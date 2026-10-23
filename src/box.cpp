#include "box.h"
#include "martini_internal.h"
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

// Every setting is required: the target pressures, the rounds between applications and the largest step of each move
static BarostatSettings read_barostat_settings(hid_t root) {
    BarostatSettings s;
    auto grp = open_group(root, "/input/barostat");
    s.target_p_xy = read_attribute<float>(grp.get(), ".", "target_p_xy");
    s.target_p_z  = read_attribute<float>(grp.get(), ".", "target_p_z");
    s.interval    = read_attribute<int>  (grp.get(), ".", "interval");
    s.mc_dmax_xy  = read_attribute<float>(grp.get(), ".", "mc_dmax_xy");
    s.mc_dmax_z   = read_attribute<float>(grp.get(), ".", "mc_dmax_z");
    if(s.interval < 1) throw std::string("/input/barostat@interval must be at least 1");
    if(!(s.mc_dmax_xy >= 0.f && s.mc_dmax_z >= 0.f)) throw std::string("/input/barostat@mc_dmax_* must be nonnegative");
    return s;
}

void unregister_barostat_for_engine(DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    g_baro_state.erase(&engine);
}

void register_barostat_for_engine(hid_t config_root, DerivEngine& engine, uint32_t random_seed) {
    if(!h5_exists(config_root, "/input/barostat")) return;
    BarostatSettings s = read_barostat_settings(config_root);

    if(!h5_exists(config_root, "/input/potential/martini_potential"))
        throw std::string("the barostat needs a periodic box (x_len, y_len, z_len on martini_potential)");
    auto grp = open_group(config_root, "/input/potential/martini_potential");
    float bx = read_attribute<float>(grp.get(), ".", "x_len");
    float by = read_attribute<float>(grp.get(), ".", "y_len");
    float bz = read_attribute<float>(grp.get(), ".", "z_len");
    if(!(bx > 0.f && by > 0.f && bz > 0.f))
        throw std::string("the barostat needs a positive periodic box");

    // Monte-Carlo barostat: group atoms by molecule (COM scaling)
    int n_atom = engine.pos->n_atom;
    if(!h5_exists(config_root, "/input/molecule_ids"))
        throw std::string("Monte-Carlo barostat requires /input/molecule_ids for COM scaling");
    check_size(config_root, "/input/molecule_ids", n_atom);
    std::vector<int> mol_id;
    traverse_dset<1,int>(config_root, "/input/molecule_ids", [&](size_t, int m){ mol_id.push_back(m); });
    int max_id = *std::max_element(mol_id.begin(), mol_id.end());
    if(max_id < 0 || *std::min_element(mol_id.begin(), mol_id.end()) < 0)
        throw std::string("/input/molecule_ids must be nonnegative");

    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto& st = g_baro_state[&engine];
    st.settings = s;
    st.box_x = bx;
    st.box_y = by;
    st.box_z = bz;
    st.mc_mol_atoms.assign(static_cast<size_t>(max_id + 1), {});
    for(int a = 0; a < n_atom; ++a)
        st.mc_mol_atoms[static_cast<size_t>(mol_id[a])].push_back(a);
    st.mc_seed = random_seed;
}

void register_node_box_updater(NodeBoxUpdater updater) {
    g_node_box_updater.store(updater, std::memory_order_relaxed);
}

// The hybrid position node recomputes each residue's BB proxy and backbone O from the protein carriers at
// every evaluation. Those slots are held fixed only so that the integrator leaves them alone: they move with
// their molecule and do not pin it in space.
static std::vector<unsigned char> hybrid_derived_slots(const DerivEngine& engine, int n_atom) {
    std::vector<unsigned char> derived(static_cast<size_t>(n_atom), 0);
    auto st = martini_hybrid::get_state_for_coord(*engine.pos);
    if(!st || !st->active) return derived;
    for(size_t k = 0; k < st->bb_atom_index.size() && k < st->bb_reference_runtime_atom_indices.size(); ++k) {
        int bb = st->bb_atom_index[k];
        int o  = st->bb_reference_runtime_atom_indices[k][3];
        if(bb >= 0 && bb < n_atom) derived[bb] = 1;
        if(o  >= 0 && o  < n_atom) derived[o]  = 1;
    }
    return derived;
}

// Monte-Carlo barostat: one lateral (xy) + one normal (z) trial move per call. Each move draws the log of
// the box-length scale uniformly and scales molecule centres, internal geometry preserved, so dU is purely
// intermolecular -- the modulus reflects membrane area elasticity, not bond stiffness. A move is accepted on
// the exact NPT weight exp(-beta[dU + P*dV]) * (V'/V)^(N+1): one factor V'/V for each of the N scaled
// centres and one because the proposal is uniform in ln V. A molecule holding a fixed atom (other than a
// derived hybrid slot) is not scaled; one holding a z-fixed atom is scaled laterally but not by the normal move.
// A molecule's centre is taken over its own particles. The derived slots move with their molecule but do not place
// it: the integrator never moves them, so their stored positions are stale, and counting them would make the
// centre's scaling a map whose Jacobian is not the V'/V the acceptance assumes.
static void apply_mc_barostat(BarostatState& st, DerivEngine& engine, int n_atom, float kT) {
    if(kT <= 0.f) return;
    auto& s = st.settings;
    VecArray pos = engine.pos->output;

    std::vector<unsigned char> derived = hybrid_derived_slots(engine, n_atom);
    std::vector<unsigned char> fixed(static_cast<size_t>(n_atom), 0);
    std::vector<unsigned char> z_fixed(static_cast<size_t>(n_atom), 0);
    for(int a : martini_fix_rigid::get_fixed_atoms(engine))   if(a >= 0 && a < n_atom && !derived[a]) fixed[a] = 1;
    for(int a : martini_fix_rigid::get_z_fixed_atoms(engine)) if(a >= 0 && a < n_atom) z_fixed[a] = 1;

    std::vector<float> saved(3 * static_cast<size_t>(n_atom));
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
            bool held = false;
            for(int a : atoms) if(fixed[a] || (axis == 1 && z_fixed[a])) { held = true; break; }
            if(held) continue;
            // centre of the molecule made whole: each of its particles taken at its minimum image from the first
            auto first = std::find_if(atoms.begin(), atoms.end(), [&](int a) { return !derived[a]; });
            if(first == atoms.end()) throw std::string("a molecule of /input/molecule_ids holds only derived hybrid slots");
            int a0 = *first;
            float x0 = pos(0,a0), y0 = pos(1,a0), z0 = pos(2,a0);
            double cx = 0, cy = 0, cz = 0;
            int n_centre = 0;
            for(int a : atoms) {
                if(derived[a]) continue;
                float dx = pos(0,a) - x0, dy = pos(1,a) - y0, dz = pos(2,a) - z0;
                minimum_image_scalar(dx, dy, dz, bx0, by0, bz0);
                cx += dx; cy += dy; cz += dz;
                ++n_centre;
            }
            double inv = 1.0 / n_centre;
            cx = x0 + cx*inv; cy = y0 + cy*inv; cz = z0 + cz*inv;
            // An atom stored n box lengths from its place in the whole molecule moves by (s-1)(centre + n box):
            // the scaled configuration is the same whichever periodic images the atoms are stored in.
            for(int a : atoms) {
                float nx = roundf((pos(0,a) - x0) / bx0);
                float ny = roundf((pos(1,a) - y0) / by0);
                float nz = roundf((pos(2,a) - z0) / bz0);
                pos(0,a) += float((sx - 1.0) * (cx + nx*bx0));
                pos(1,a) += float((sy - 1.0) * (cy + ny*by0));
                pos(2,a) += float((sz - 1.0) * (cz + nz*bz0));
            }
            ++n_scaled;
        }
        if(n_scaled == 0) continue;

        st.box_x = bx0 * sx; st.box_y = by0 * sy; st.box_z = bz0 * sz;
        update_node_boxes(engine, st.box_x, st.box_y, st.box_z);
        double V_new = double(st.box_x) * st.box_y * st.box_z;

        engine.compute(PotentialAndDerivMode);
        double U_new = engine.potential;

        double arg = -(1.0/kT) * ((U_new - U_old) + double(P) * (V_new - V_old))
                     + double(n_scaled + 1) * std::log(V_new / V_old);
        if(std::log(uni(st.mc_rng)) < arg) {
            U_old = U_new;                                        // accept
        } else {
            for(int i = 0; i < n_atom; ++i) { pos(0,i) = saved[3*i]; pos(1,i) = saved[3*i+1]; pos(2,i) = saved[3*i+2]; }
            st.box_x = bx0; st.box_y = by0; st.box_z = bz0;
            update_node_boxes(engine, bx0, by0, bz0);            // restore every box exactly
        }
    }
}

void maybe_apply_barostat(DerivEngine& engine,
                          int n_atom,
                          uint64_t round_num,
                          float temperature) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(&engine);
    if(it == g_baro_state.end()) return;

    auto& st = it->second;
    auto& s = st.settings;

    if(s.interval <= 0) return;
    if(round_num == 0 || (round_num % s.interval) != 0) return;

    // Seed the trial moves from the round, so that a continuation draws the moves of the continuous run.
    std::seed_seq seed{st.mc_seed, uint32_t(round_num), uint32_t(round_num >> 32)};
    st.mc_rng.seed(seed);
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

void update_node_boxes(DerivEngine& engine, float box_x, float box_y, float box_z) {
    NodeBoxUpdater updater = g_node_box_updater.load(std::memory_order_relaxed);
    if(updater) updater(engine, box_x, box_y, box_z);
}

float get_volume(const DerivEngine& engine) {
    float bx, by, bz;
    get_current_box(engine, bx, by, bz);
    return bx * by * bz;
}

} // namespace npt

} // namespace simulation_box
