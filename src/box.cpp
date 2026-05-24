#include "box.h"
#include "martini.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <map>
#include <vector>

using namespace h5;

namespace simulation_box {

void wrap_positions(VecArray pos, int n_atom, float box_x, float box_y, float box_z) {
    for(int i = 0; i < n_atom; ++i) {
        float x = pos(0, i);
        float y = pos(1, i);
        float z = pos(2, i);
        
        // Wrap into [0, box_dim)
        x = fmodf(x, box_x);
        if(x < 0) x += box_x;
        y = fmodf(y, box_y);
        if(y < 0) y += box_y;
        z = fmodf(z, box_z);
        if(z < 0) z += box_z;
        
        pos(0, i) = x;
        pos(1, i) = y;
        pos(2, i) = z;
    }
}

namespace npt {

// Global registry mapping DerivEngine* to BarostatState
static std::mutex g_baro_mutex;
static std::map<DerivEngine*, BarostatState> g_baro_state;
static std::atomic<NodeBoxUpdater> g_node_box_updater{nullptr};

static inline float volume_xyz(float x, float y, float z) { return x * y * z; }

// Attempt to read barostat settings from H5 config
static BarostatSettings read_barostat_settings(hid_t root) {
    BarostatSettings s;
    if(h5_exists(root, "/input/barostat")) {
        auto grp = open_group(root, "/input/barostat");
        s.enabled = read_attribute<int>(grp.get(), ".", "enable", int(s.enabled)) != 0;
        s.target_p_xy = read_attribute<float>(grp.get(), ".", "target_p_xy", s.target_p_xy);
        s.target_p_z = read_attribute<float>(grp.get(), ".", "target_p_z", s.target_p_z);
        s.tau_p = read_attribute<float>(grp.get(), ".", "tau_p", s.tau_p);
        s.interval = read_attribute<int>(grp.get(), ".", "interval", s.interval);
        s.compressibility = read_attribute<float>(grp.get(), ".", "compressibility", s.compressibility);
        s.compressibility_xy = read_attribute<float>(grp.get(), ".", "compressibility_xy", s.compressibility);
        s.compressibility_z = read_attribute<float>(grp.get(), ".", "compressibility_z", s.compressibility);
        s.semi_isotropic = read_attribute<int>(grp.get(), ".", "semi_isotropic", int(s.semi_isotropic)) != 0;
        int type_int = read_attribute<int>(grp.get(), ".", "type", 0);
        if(type_int == 1) {
            s.type = BarostatType::ParrinelloRahman;
        } else {
            s.type = BarostatType::Berendsen;
        }
    }
    return s;
}

// Scale positions for semi-isotropic barostat
static void apply_semi_isotropic_scaling(DerivEngine& engine, float scale_xy, float scale_z) {
    VecArray pos = engine.pos->output;
    int n_atom = engine.pos->n_elem;

    std::vector<unsigned char> fixed_mask(static_cast<size_t>(std::max(0, n_atom)), 0);
    std::vector<unsigned char> z_fixed_mask(static_cast<size_t>(std::max(0, n_atom)), 0);
    for(int atom_idx : martini_fix_rigid::get_fixed_atoms(engine)) {
        if(atom_idx >= 0 && atom_idx < n_atom) {
            fixed_mask[static_cast<size_t>(atom_idx)] = 1;
        }
    }
    for(int atom_idx : martini_fix_rigid::get_z_fixed_atoms(engine)) {
        if(atom_idx >= 0 && atom_idx < n_atom) {
            z_fixed_mask[static_cast<size_t>(atom_idx)] = 1;
        }
    }

    for(int i = 0; i < n_atom; ++i) {
        if(fixed_mask[static_cast<size_t>(i)]) continue;
        pos(0, i) *= scale_xy;
        pos(1, i) *= scale_xy;
        if(!z_fixed_mask[static_cast<size_t>(i)]) {
            pos(2, i) *= scale_z;
        }
    }
    
    update_node_boxes(engine, scale_xy, scale_z);
}

// Estimate pressure tensor from kinetic + virial contributions
static void estimate_pressure_tensor(const VecArray& mom,
                                     const VecArray& pos,
                                     const VecArray& sens,
                                     const std::vector<float>& masses,
                                     int n_atom,
                                     float V,
                                     const std::vector<unsigned char>* fixed_mask,
                                     const std::vector<unsigned char>* z_fixed_mask,
                                     float& pxx, float& pyy, float& pzz) {
    // Kinetic contribution: sum_i (p_i^2 / m_i)
    double kx = 0.0, ky = 0.0, kz = 0.0;
    // Virial contribution: -sum_i (r_i . F_i) where F = -deriv
    double wxx = 0.0, wyy = 0.0, wzz = 0.0;

    for(int i = 0; i < n_atom; ++i) {
        bool is_fixed = (fixed_mask && !fixed_mask->empty() && (*fixed_mask)[static_cast<size_t>(i)]);
        if(is_fixed) continue;
        bool is_z_fixed = (z_fixed_mask && !z_fixed_mask->empty() && (*z_fixed_mask)[static_cast<size_t>(i)]);

        float3 p = load_vec<3>(mom, i);
        float m = (masses.size() > (size_t)i) ? masses[i] : 1.0f;
        if(m <= 0.f) m = 1.0f;

        kx += double(p.x() * p.x()) / m;
        ky += double(p.y() * p.y()) / m;
        if(!is_z_fixed) kz += double(p.z() * p.z()) / m;

        // Virial: r . (-deriv) = -r . sens (sens is gradient)
        float3 r = load_vec<3>(pos, i);
        float3 F = load_vec<3>(sens, i);
        F = make_vec3(-F.x(), -F.y(), -F.z()); // Force = -gradient

        wxx += double(r.x() * F.x());
        wyy += double(r.y() * F.y());
        if(!is_z_fixed) wzz += double(r.z() * F.z());
    }

    pxx = float((kx + wxx) / V);
    pyy = float((ky + wyy) / V);
    pzz = float((kz + wzz) / V);
}

// Apply Berendsen barostat scaling
static void apply_berendsen_barostat(BarostatState& st,
                                     float pxy_inst, float pz_inst,
                                     float delta_t,
                                     float& scale_xy, float& scale_z) {
    auto& s = st.settings;
    scale_xy = 1.0f;
    scale_z = 1.0f;

    if(s.semi_isotropic) {
        float beta_xy = std::max(0.0f, s.compressibility_xy);
        float beta_z = std::max(0.0f, s.compressibility_z);
        float factor_xy = 1.0f;
        float factor_z = 1.0f;
        if(beta_xy > 0.0f && std::isfinite(pxy_inst) && std::isfinite(s.target_p_xy)) {
            factor_xy = 1.f + beta_xy * (delta_t / s.tau_p) * (pxy_inst - s.target_p_xy);
        }
        if(beta_z > 0.0f && std::isfinite(pz_inst) && std::isfinite(s.target_p_z)) {
            factor_z = 1.f + beta_z * (delta_t / s.tau_p) * (pz_inst - s.target_p_z);
        }

        // Clamp scale factors
        factor_xy = std::max(0.98f, std::min(1.02f, factor_xy));
        factor_z = std::max(0.98f, std::min(1.02f, factor_z));

        scale_xy = powf(factor_xy, 0.5f);
        scale_z = factor_z;

        scale_xy = std::max(0.995f, std::min(1.005f, scale_xy));
        scale_z = std::max(0.995f, std::min(1.005f, scale_z));

        // Monotonic shrink guard
        if(pxy_inst < s.target_p_xy) scale_xy = std::min(scale_xy, 1.0f);
        if(pz_inst < s.target_p_z) scale_z = std::min(scale_z, 1.0f);
    } else {
        float beta = std::max(0.0f, s.compressibility);
        float p_inst = (2.f * pxy_inst + pz_inst) / 3.f;
        float target_p = (2.f * s.target_p_xy + s.target_p_z) / 3.f;
        float factor = 1.0f;
        if(beta > 0.0f && std::isfinite(p_inst) && std::isfinite(target_p)) {
            factor = 1.f + beta * (delta_t / s.tau_p) * (p_inst - target_p);
        }
        factor = std::max(0.98f, std::min(1.02f, factor));
        float sc = powf(factor, 1.f / 3.f);
        sc = std::max(0.995f, std::min(1.005f, sc));
        scale_xy = sc;
        scale_z = sc;
    }

    // First application policy
    if(!st.has_applied_once && s.prefer_shrink_first) {
        scale_xy = std::min(scale_xy, 1.0f);
        scale_z = std::min(scale_z, 1.0f);
    }
}

// Apply Parrinello-Rahman barostat
static void apply_parrinello_rahman_barostat(BarostatState& st,
                                             float pxy_inst, float pz_inst,
                                             float delta_t,
                                             float& scale_xy, float& scale_z) {
    auto& s = st.settings;

    // Update box velocities using pressure difference
    // Correct Parrinello-Rahman equation: dv/dt = (1/W) * (P - P_target)
    // Note: Volume V should NOT be in the numerator
    float bx = st.box_x, by = st.box_y, bz = st.box_z;

    if(s.semi_isotropic) {
        float beta_xy = std::max(0.0f, s.compressibility_xy);
        float beta_z = std::max(0.0f, s.compressibility_z);
        float W_xy = (beta_xy > 0.0f) ? (s.tau_p * s.tau_p / beta_xy) : 0.0f;
        float W_z = (beta_z > 0.0f) ? (s.tau_p * s.tau_p / beta_z) : 0.0f;

        // Lateral (xy) and normal (z) are independent
        float dv_xy_dt = (W_xy > 0.0f && std::isfinite(pxy_inst) && std::isfinite(s.target_p_xy))
            ? ((1.0f / W_xy) * (pxy_inst - s.target_p_xy)) : 0.0f;
        float dv_z_dt = (W_z > 0.0f && std::isfinite(pz_inst) && std::isfinite(s.target_p_z))
            ? ((1.0f / W_z) * (pz_inst - s.target_p_z)) : 0.0f;

        // Update velocities
        st.box_vel_xy += dv_xy_dt * delta_t;
        st.box_vel_z += dv_z_dt * delta_t;

        // Apply damping to prevent oscillations
        float damping = 0.95f;
        st.box_vel_xy *= damping;
        st.box_vel_z *= damping;

        // Compute scaling factors from velocities
        // scale = exp(v * dt / box_dim)
        // For small v*dt, scale is approximately 1 + v*dt/box_dim.
        float avg_box_xy = 0.5f * (bx + by);
        scale_xy = 1.0f + (st.box_vel_xy * delta_t) / avg_box_xy;
        scale_z = 1.0f + (st.box_vel_z * delta_t) / bz;

        if(beta_xy == 0.0f) scale_xy = 1.0f;
        if(beta_z == 0.0f) scale_z = 1.0f;

        // Clamp to prevent instabilities
        scale_xy = std::max(0.99f, std::min(1.01f, scale_xy));
        scale_z = std::max(0.99f, std::min(1.01f, scale_z));
    } else {
        float beta = std::max(0.0f, s.compressibility);
        float W = (beta > 0.0f) ? (s.tau_p * s.tau_p / beta) : 0.0f;
        // Isotropic coupling
        float p_inst = (2.f * pxy_inst + pz_inst) / 3.f;
        float target_p = (2.f * s.target_p_xy + s.target_p_z) / 3.f;

        float dv_dt = (W > 0.0f && std::isfinite(p_inst) && std::isfinite(target_p))
            ? ((1.0f / W) * (p_inst - target_p)) : 0.0f;
        st.box_vel_xy += dv_dt * delta_t;
        st.box_vel_z = st.box_vel_xy;

        float damping = 0.95f;
        st.box_vel_xy *= damping;
        st.box_vel_z *= damping;

        float avg_box = powf(bx * by * bz, 1.0f/3.0f);
        float scale = 1.0f + (st.box_vel_xy * delta_t) / avg_box;
        scale = (beta > 0.0f) ? std::max(0.99f, std::min(1.01f, scale)) : 1.0f;
        scale_xy = scale;
        scale_z = scale;
    }
}

void register_barostat_for_engine(hid_t config_root, DerivEngine& engine) {
    BarostatSettings s = read_barostat_settings(config_root);
    if(!s.enabled) return;
    
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
    
    // Load masses if available
    try {
        if(h5_exists(config_root, "/input/mass")) {
            st.masses.clear();
            traverse_dset<1,float>(config_root, "/input/mass", 
                [&](size_t i, float m){ st.masses.push_back(m); });
        }
    } catch(...) {}
    
}

void register_node_box_updater(NodeBoxUpdater updater) {
    g_node_box_updater.store(updater, std::memory_order_relaxed);
}

void maybe_apply_barostat(DerivEngine& engine,
                          const VecArray& mom,
                          int n_atom,
                          uint64_t round_num,
                          float dt,
                          int inner_step,
                          int verbose,
                          bool print_now) {
    (void)verbose;
    (void)print_now;
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(&engine);
    if(it == g_baro_state.end()) return;
    
    auto& st = it->second;
    auto& s = st.settings;
    if(!s.enabled) return;

    if(s.interval <= 0) return;
    if(round_num == 0 || (round_num % s.interval) != 0) return;
    
    // Determine current box
    float bx = st.box_x, by = st.box_y, bz = st.box_z;
    if(bx == 0.f || by == 0.f || bz == 0.f) return;
    float V = volume_xyz(bx, by, bz);
    
    // Compute forces for virial
    engine.compute(PotentialAndDerivMode);
    VecArray pos = engine.pos->output;
    VecArray sens = engine.pos->sens;
    
    std::vector<unsigned char> fixed_mask(static_cast<size_t>(std::max(0, n_atom)), 0);
    std::vector<unsigned char> z_fixed_mask(static_cast<size_t>(std::max(0, n_atom)), 0);
    for(int atom_idx : martini_fix_rigid::get_fixed_atoms(engine)) {
        if(atom_idx >= 0 && atom_idx < n_atom) {
            fixed_mask[static_cast<size_t>(atom_idx)] = 1;
        }
    }
    for(int atom_idx : martini_fix_rigid::get_z_fixed_atoms(engine)) {
        if(atom_idx >= 0 && atom_idx < n_atom) {
            z_fixed_mask[static_cast<size_t>(atom_idx)] = 1;
        }
    }

    float pxx = 0.f, pyy = 0.f, pzz = 0.f;
    estimate_pressure_tensor(mom, pos, sens, st.masses, n_atom, V,
                             &fixed_mask, &z_fixed_mask,
                             pxx, pyy, pzz);
    
    float pxy_inst = 0.5f * (pxx + pyy);
    float pz_inst = pzz;

    if(!std::isfinite(pxy_inst) || !std::isfinite(pz_inst)) return;

    // Store for logging
    st.last_pxy_inst = pxy_inst;
    st.last_pz_inst = pz_inst;

    // Calculate time step
    float delta_t = dt * inner_step;

    // Apply appropriate barostat algorithm
    float scale_xy = 1.0f;
    float scale_z = 1.0f;

    if(s.type == BarostatType::ParrinelloRahman) {
        apply_parrinello_rahman_barostat(st, pxy_inst, pz_inst, delta_t, scale_xy, scale_z);
    } else {
        apply_berendsen_barostat(st, pxy_inst, pz_inst, delta_t, scale_xy, scale_z);
    }
    
    apply_semi_isotropic_scaling(engine, scale_xy, scale_z);
    st.box_x = bx * scale_xy;
    st.box_y = by * scale_xy;
    st.box_z = bz * scale_z;
    st.has_applied_once = true;

    // Store scale factors for logging
    st.last_scale_xy = scale_xy;
    st.last_scale_z = scale_z;
}

void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz) {
    bx = by = bz = 0.f;
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_baro_state.end() && it->second.settings.enabled) {
        bx = it->second.box_x;
        by = it->second.box_y;
        bz = it->second.box_z;
    }
}

bool is_enabled(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_baro_state.end()) {
        return it->second.settings.enabled;
    }
    return false;
}

void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z) {
    NodeBoxUpdater updater = g_node_box_updater.load(std::memory_order_relaxed);
    if(updater) updater(engine, scale_xy, scale_z);
}

void get_pressure(const DerivEngine& engine, float& pxy, float& pz) {
    std::lock_guard<std::mutex> lk(g_baro_mutex);
    auto it = g_baro_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_baro_state.end() && it->second.settings.enabled) {
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
