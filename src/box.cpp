#include "box.h"
#include "martini.h"
#include <cmath>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <map>
#include <vector>

using namespace h5;

namespace simulation_box {

namespace npt {

// Global registry mapping DerivEngine* to BarostatState
static std::mutex g_baro_mutex;
static std::map<DerivEngine*, BarostatState> g_baro_state;
static std::atomic<NodeBoxUpdater> g_node_box_updater{nullptr};

// Attempt to read barostat settings from H5 config
static BarostatSettings read_barostat_settings(hid_t root) {
    BarostatSettings s;
    if(h5_exists(root, "/input/barostat")) {
        auto grp = open_group(root, "/input/barostat");
        s.target_p_xy = read_attribute<float>(grp.get(), ".", "target_p_xy", s.target_p_xy);
        s.target_p_z = read_attribute<float>(grp.get(), ".", "target_p_z", s.target_p_z);
        s.tau_p = read_attribute<float>(grp.get(), ".", "tau_p", s.tau_p);
        s.interval = read_attribute<int>(grp.get(), ".", "interval", s.interval);
        s.compressibility = read_attribute<float>(grp.get(), ".", "compressibility", s.compressibility);
        s.compressibility_xy = read_attribute<float>(grp.get(), ".", "compressibility_xy", s.compressibility);
        s.compressibility_z = read_attribute<float>(grp.get(), ".", "compressibility_z", s.compressibility);
        s.semi_isotropic = read_attribute<int>(grp.get(), ".", "semi_isotropic", int(s.semi_isotropic)) != 0;
        s.mc_dmax_xy = read_attribute<float>(grp.get(), ".", "mc_dmax_xy", s.mc_dmax_xy);
        s.mc_dmax_z  = read_attribute<float>(grp.get(), ".", "mc_dmax_z",  s.mc_dmax_z);
        s.mc_seed    = read_attribute<unsigned>(grp.get(), ".", "mc_seed", s.mc_seed);
    }
    return s;
}

void register_barostat_for_engine(hid_t config_root, DerivEngine& engine) {
    if(!h5_exists(config_root, "/input/barostat")) return;
    BarostatSettings s = read_barostat_settings(config_root);

    float bx = 0.f, by = 0.f, bz = 0.f;
    try {
        // Try to get box dimensions from various sources
        if(h5_exists(config_root, "/input/potential/periodic_boundary_potential")) {
            auto grp = open_group(config_root, "/input/potential/periodic_boundary_potential");
            bx = read_attribute<float>(grp.get(), ".", "x_len");
            by = read_attribute<float>(grp.get(), ".", "y_len");
            bz = read_attribute<float>(grp.get(), ".", "z_len");
        } else if(h5_exists(config_root, "/input/potential/martini_potential")) {
            auto grp = open_group(config_root, "/input/potential/martini_potential");
            bx = read_attribute<float>(grp.get(), ".", "x_len");
            by = read_attribute<float>(grp.get(), ".", "y_len");
            bz = read_attribute<float>(grp.get(), ".", "z_len");
        }
    } catch(...) {
        bx = by = bz = 0.f;
    }
    
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto& st = g_baro_state[&engine];
    st.settings = s;
    st.box_x = bx;
    st.box_y = by;
    st.box_z = bz;

    // Monte-Carlo barostat: group atoms by molecule (COM scaling) + seed the RNG
    st.mc_mol_atoms.clear();
    st.mc_init = false;
    try {
        if(h5_exists(config_root, "/input/molecule_ids")) {
            std::vector<int> mol_id;
            traverse_dset<1,int>(config_root, "/input/molecule_ids",
                [&](size_t, int m){ mol_id.push_back(m); });
            int max_id = -1;
            for(int m : mol_id) max_id = std::max(max_id, m);
            if(max_id >= 0) {
                st.mc_mol_atoms.assign(static_cast<size_t>(max_id + 1), {});
                for(size_t a = 0; a < mol_id.size(); ++a)
                    st.mc_mol_atoms[static_cast<size_t>(mol_id[a])].push_back(static_cast<int>(a));
                st.mc_rng.seed(s.mc_seed);
                st.mc_init = true;
            }
        }
    } catch(...) { st.mc_mol_atoms.clear(); st.mc_init = false; }
    if(!st.mc_init)
        throw std::string("Monte-Carlo barostat requires /input/molecule_ids for COM scaling");
}

void register_node_box_updater(NodeBoxUpdater updater) {
    g_node_box_updater.store(updater, std::memory_order_relaxed);
}

// Monte-Carlo barostat: one lateral (xy) + one normal (z) trial move per call, each accepting on the
// exact NPT weight exp(-beta[dU + P*dV]) * (V'/V)^N_mol. Molecule COMs are scaled (internal geometry
// preserved), so dU is purely intermolecular -- the modulus reflects membrane area elasticity, not bond
// stiffness. Correct fluctuations by construction, independent of any coupling constant.
static void apply_mc_barostat(BarostatState& st, DerivEngine& engine, int n_atom, float kT) {
    if(!st.mc_init || st.mc_mol_atoms.empty() || kT <= 0.f) return;
    auto& s = st.settings;
    VecArray pos = engine.pos->output;

    std::vector<unsigned char> fixed(static_cast<size_t>(std::max(0, n_atom)), 0);
    for(int a : martini_fix_rigid::get_fixed_atoms(engine))   if(a >= 0 && a < n_atom) fixed[a] = 1;
    for(int a : martini_fix_rigid::get_z_fixed_atoms(engine)) if(a >= 0 && a < n_atom) fixed[a] = 1;

    std::vector<float> saved(3 * static_cast<size_t>(std::max(0, n_atom)));
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    engine.compute(PotentialAndDerivMode);
    double U_old = engine.potential;

    for(int axis = 0; axis < 2; ++axis) {         // 0 = lateral (xy), 1 = normal (z)
        float dmax = (axis == 0) ? s.mc_dmax_xy : s.mc_dmax_z;
        if(dmax <= 0.f) continue;
        float P = (axis == 0) ? s.target_p_xy : s.target_p_z;
        float scale = std::exp(static_cast<float>(dmax * (2.0 * uni(st.mc_rng) - 1.0)));
        float sx = (axis == 0) ? scale : 1.f;
        float sy = (axis == 0) ? scale : 1.f;
        float sz = (axis == 0) ? 1.f   : scale;

        float bx0 = st.box_x, by0 = st.box_y, bz0 = st.box_z;
        double V_old = double(bx0) * by0 * bz0;

        for(int i = 0; i < n_atom; ++i) { saved[3*i] = pos(0,i); saved[3*i+1] = pos(1,i); saved[3*i+2] = pos(2,i); }
        long n_scaled = 0;
        for(auto& atoms : st.mc_mol_atoms) {
            if(atoms.empty()) continue;
            bool anyfix = false;
            for(int a : atoms) if(a < n_atom && fixed[a]) { anyfix = true; break; }
            if(anyfix) continue;
            double cx = 0, cy = 0, cz = 0;
            for(int a : atoms) { cx += pos(0,a); cy += pos(1,a); cz += pos(2,a); }
            double inv = 1.0 / atoms.size(); cx *= inv; cy *= inv; cz *= inv;
            float dx = float(cx*sx - cx), dy = float(cy*sy - cy), dz = float(cz*sz - cz);
            for(int a : atoms) { pos(0,a) += dx; pos(1,a) += dy; pos(2,a) += dz; }
            ++n_scaled;
        }
        if(n_scaled == 0) continue;

        st.box_x *= sx; st.box_y *= sy; st.box_z *= sz;
        update_node_boxes(engine, sx, sz);       // xy scale = sx (=sy), z scale = sz
        double V_new = double(st.box_x) * st.box_y * st.box_z;

        engine.compute(PotentialAndDerivMode);
        double U_new = engine.potential;

        double arg = -(1.0/kT) * ((U_new - U_old) + double(P) * (V_new - V_old))
                     + double(n_scaled) * std::log(V_new / V_old);
        st.mc_attempt[axis]++;
        if(std::log(uni(st.mc_rng)) < arg) {
            U_old = U_new; st.mc_accept[axis]++;                 // accept
        } else {
            for(int i = 0; i < n_atom; ++i) { pos(0,i) = saved[3*i]; pos(1,i) = saved[3*i+1]; pos(2,i) = saved[3*i+2]; }
            st.box_x = bx0; st.box_y = by0; st.box_z = bz0;      // restore box exactly
            update_node_boxes(engine, 1.f/sx, 1.f/sz);           // undo node-box scaling
        }
    }
    st.last_scale_xy = 1.f; st.last_scale_z = 1.f;
}

void maybe_apply_barostat(DerivEngine& engine,
                          const VecArray& mom,
                          int n_atom,
                          uint64_t round_num,
                          float dt,
                          int inner_step,
                          float temperature,
                          int verbose,
                          bool print_now) {
    (void)mom;
    (void)dt;
    (void)inner_step;
    (void)verbose;
    (void)print_now;
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(&engine);
    if(it == g_baro_state.end()) return;

    auto& st = it->second;
    auto& s = st.settings;

    if(s.interval <= 0) return;
    if(round_num == 0 || (round_num % s.interval) != 0) return;

    // Box must be initialized before scaling
    if(st.box_x == 0.f || st.box_y == 0.f || st.box_z == 0.f) return;

    apply_mc_barostat(st, engine, n_atom, temperature);
}

void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz) {
    bx = by = bz = 0.f;
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_baro_state.end()) {
        bx = it->second.box_x;
        by = it->second.box_y;
        bz = it->second.box_z;
    }
}

bool is_enabled(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    return it != g_baro_state.end();
}

void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z) {
    NodeBoxUpdater updater = g_node_box_updater.load(std::memory_order_relaxed);
    if(updater) updater(engine, scale_xy, scale_z);
}

void get_pressure(const DerivEngine& engine, float& pxy, float& pz) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_baro_state.end()) {
        pxy = it->second.last_pxy_inst;
        pz = it->second.last_pz_inst;
    } else {
        pxy = 0.0f;
        pz = 0.0f;
    }
}

float get_volume(const DerivEngine& engine) {
    float bx, by, bz;
    get_current_box(engine, bx, by, bz);
    return bx * by * bz;
}

} // namespace npt

} // namespace simulation_box
