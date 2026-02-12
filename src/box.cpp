// box.cpp - Simulation box management for periodic boundaries and NPT/NVT ensembles
// Ported from origin/martini branch
//
// H5 configuration expected at: /input/barostat (group)
//   attrs:
//     enable            (int, 0/1)
//     target_p_xy       (float)   target lateral pressure (UP units)
//     target_p_z        (float)   target normal pressure  (UP units)
//     tau_p             (float)   barostat time constant (time units)
//     interval          (int)     apply every N integrator steps
//     compressibility_xy(float)   in-plane compressibility (1/pressure)
//     compressibility_z (float)   normal compressibility   (1/pressure)
//     compressibility   (float)   legacy isotropic compressibility fallback
//     semi_isotropic    (int)     1 = scale x,y together; z separately
//     debug             (int)     enable prints alongside existing output

#include "box.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <mutex>
#include <map>
#include <vector>

using namespace h5;

// Local helper to check if an HDF5 attribute exists
// (Mirrors the function in martini.cpp)
static inline bool attribute_exists(hid_t loc_id, const char* obj_name, const char* attr_name) {
    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

namespace simulation_box {

// ===================== PBC UTILITIES =====================

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

// ===================== NPT BAROSTAT =====================

namespace npt {

// Global registry mapping DerivEngine* to BarostatState
static std::mutex g_baro_mutex;
static std::map<DerivEngine*, BarostatState> g_baro_state;

static inline float volume_xyz(float x, float y, float z) { return x * y * z; }
static inline float pressure_rel_deviation(float value, float target) {
    const float abs_tol = 1e-6f;
    float denom = std::max(fabsf(target), abs_tol);
    return fabsf(value - target) / denom;
}

// Attempt to read barostat settings from H5 config
static BarostatSettings read_barostat_settings(hid_t root) {
    BarostatSettings s;
    try {
        if(h5_exists(root, "/input/barostat")) {
            auto grp = open_group(root, "/input/barostat");
            if(attribute_exists(grp.get(), ".", "enable"))
                s.enabled = read_attribute<int>(grp.get(), ".", "enable") != 0;
            if(attribute_exists(grp.get(), ".", "target_p_xy"))
                s.target_p_xy = read_attribute<float>(grp.get(), ".", "target_p_xy");
            if(attribute_exists(grp.get(), ".", "target_p_z"))
                s.target_p_z = read_attribute<float>(grp.get(), ".", "target_p_z");
            if(attribute_exists(grp.get(), ".", "tau_p"))
                s.tau_p = read_attribute<float>(grp.get(), ".", "tau_p");
            if(attribute_exists(grp.get(), ".", "interval"))
                s.interval = read_attribute<int>(grp.get(), ".", "interval");
            if(attribute_exists(grp.get(), ".", "compressibility"))
                s.compressibility = read_attribute<float>(grp.get(), ".", "compressibility");
            s.compressibility_xy = s.compressibility;
            s.compressibility_z = s.compressibility;
            if(attribute_exists(grp.get(), ".", "compressibility_xy"))
                s.compressibility_xy = read_attribute<float>(grp.get(), ".", "compressibility_xy");
            if(attribute_exists(grp.get(), ".", "compressibility_z"))
                s.compressibility_z = read_attribute<float>(grp.get(), ".", "compressibility_z");
            if(attribute_exists(grp.get(), ".", "semi_isotropic"))
                s.semi_isotropic = read_attribute<int>(grp.get(), ".", "semi_isotropic") != 0;
            if(attribute_exists(grp.get(), ".", "debug"))
                s.debug = read_attribute<int>(grp.get(), ".", "debug") != 0;
            // Read barostat type: 0 = Berendsen (default), 1 = Parrinello-Rahman
            if(attribute_exists(grp.get(), ".", "type")) {
                int type_int = read_attribute<int>(grp.get(), ".", "type");
                if(type_int == 1) {
                    s.type = BarostatType::ParrinelloRahman;
                } else {
                    s.type = BarostatType::Berendsen;
                }
            }
        }
    } catch(...) {
        s.enabled = false;
    }
    return s;
}

// Scale positions for semi-isotropic barostat
static void apply_semi_isotropic_scaling(DerivEngine& engine, float scale_xy, float scale_z) {
    VecArray pos = engine.pos->output;
    int n_atom = engine.pos->n_elem;
    
    for(int i = 0; i < n_atom; ++i) {
        pos(0, i) *= scale_xy;
        pos(1, i) *= scale_xy;
        pos(2, i) *= scale_z;
    }
    
    // Also update node box dimensions
    update_node_boxes(engine, scale_xy, scale_z);
}

// Estimate pressure tensor from kinetic + virial contributions
static void estimate_pressure_tensor(const VecArray& mom,
                                     const VecArray& pos,
                                     const VecArray& sens,
                                     const std::vector<float>& masses,
                                     int n_atom,
                                     float V,
                                     float& pxx, float& pyy, float& pzz) {
    // Kinetic contribution: sum_i (p_i^2 / m_i)
    double kx = 0.0, ky = 0.0, kz = 0.0;
    // Virial contribution: -sum_i (r_i . F_i) where F = -deriv
    double wxx = 0.0, wyy = 0.0, wzz = 0.0;

    for(int i = 0; i < n_atom; ++i) {
        float3 p = load_vec<3>(mom, i);
        float m = (masses.size() > (size_t)i) ? masses[i] : 1.0f;
        if(m <= 0.f) m = 1.0f;

        kx += double(p.x() * p.x()) / m;
        ky += double(p.y() * p.y()) / m;
        kz += double(p.z() * p.z()) / m;

        // Virial: r . (-deriv) = -r . sens (sens is gradient)
        float3 r = load_vec<3>(pos, i);
        float3 F = load_vec<3>(sens, i);
        F = make_vec3(-F.x(), -F.y(), -F.z()); // Force = -gradient

        wxx += double(r.x() * F.x());
        wyy += double(r.y() * F.y());
        wzz += double(r.z() * F.z());
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
        // For small v*dt, scale ≈ 1 + v*dt/box_dim
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
    
    if(s.debug) {
        const char* type_name = (s.type == BarostatType::ParrinelloRahman) ? "Parrinello-Rahman" : "Berendsen";
        printf("[NPT] Barostat registered: %s, box %.2f x %.2f x %.2f, target Pxy=%.3e Pz=%.3e, beta_xy=%.3e beta_z=%.3e, tau=%.2f, interval=%d\n",
               type_name, bx, by, bz, s.target_p_xy, s.target_p_z, s.compressibility_xy, s.compressibility_z, s.tau_p, s.interval);
    }
}

void maybe_apply_barostat(DerivEngine& engine,
                          const VecArray& mom,
                          int n_atom,
                          uint64_t round_num,
                          float dt,
                          int inner_step,
                          int verbose,
                          bool print_now) {
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
    
    float pxx = 0.f, pyy = 0.f, pzz = 0.f;
    estimate_pressure_tensor(mom, pos, sens, st.masses, n_atom, V, pxx, pyy, pzz);
    
    float pxy_inst = 0.5f * (pxx + pyy);
    float pz_inst = pzz;

    if(!std::isfinite(pxy_inst) || !std::isfinite(pz_inst)) {
        if(verbose && s.debug && print_now) {
            printf(" [NPT] skipped update due to non-finite pressure (Pxy=%.3e Pz=%.3e)\n",
                   pxy_inst, pz_inst);
        }
        return;
    }

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
    
    // Check for equilibrium
    bool at_equilibrium = false;
    if(st.has_applied_once) {
        float pxy_deviation = pressure_rel_deviation(pxy_inst, s.target_p_xy);
        float pz_deviation = pressure_rel_deviation(pz_inst, s.target_p_z);
        const float PRESSURE_TOLERANCE = 0.05f;
        
        if(pxy_deviation < PRESSURE_TOLERANCE && pz_deviation < PRESSURE_TOLERANCE) {
            st.equilibrium_count++;
        } else {
            st.equilibrium_count = 0;
        }
        
        if(st.equilibrium_count >= BarostatState::EQUILIBRIUM_THRESHOLD) {
            at_equilibrium = true;
            if(print_now && verbose) {
                printf(" [EQUILIBRIUM] Pressure converged to target\n");
            }
        }
    }
    
    st.prev_box_x = st.box_x;
    st.prev_box_y = st.box_y;
    st.prev_box_z = st.box_z;
    
    if(!at_equilibrium) {
        apply_semi_isotropic_scaling(engine, scale_xy, scale_z);
        st.box_x = bx * scale_xy;
        st.box_y = by * scale_xy;
        st.box_z = bz * scale_z;
        st.has_applied_once = true;
    }

    // Store scale factors for logging
    st.last_scale_xy = scale_xy;
    st.last_scale_z = scale_z;
    
    if(verbose && s.debug && print_now) {
        printf("[NPT] t %.3f box %.2f %.2f %.2f\n",
               double(round_num * delta_t),
               st.box_x, st.box_y, st.box_z);
        fflush(stdout);
    }
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

// Placeholder for updating node boxes - will be implemented in martini.cpp
void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z) {
    // This function should update box dimensions in MartiniPotential and other nodes
    // Implementation depends on how box dimensions are stored in these nodes
    // For now, positions are scaled directly; nodes will read from engine state
    (void)engine;
    (void)scale_xy;
    (void)scale_z;
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

// ===================== EWALD SUMMATION =====================
namespace ewald {

static std::mutex g_ewald_mutex;
static std::map<DerivEngine*, EwaldState> g_ewald_state;
struct TrigSplineTable {
    int n_grid = 0;
    std::vector<float> cos_table;
    std::vector<float> sin_table;
};
static TrigSplineTable g_trig_table;

static inline bool attribute_exists_ewald(hid_t loc_id, const char* obj_name, const char* attr_name) {
    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

static inline int wrap_idx(int i, int n) {
    int r = i % n;
    return (r < 0) ? (r + n) : r;
}

static void ensure_trig_table(int n_grid) {
    if(n_grid < 256) n_grid = 256;
    if(g_trig_table.n_grid == n_grid && !g_trig_table.cos_table.empty()) return;
    g_trig_table.n_grid = n_grid;
    g_trig_table.cos_table.resize(n_grid);
    g_trig_table.sin_table.resize(n_grid);
    const float twopi = 2.0f * float(M_PI);
    for(int i = 0; i < n_grid; ++i) {
        float a = twopi * (float(i) / float(n_grid));
        g_trig_table.cos_table[i] = cosf(a);
        g_trig_table.sin_table[i] = sinf(a);
    }
}

// Periodic cardinal cubic B-spline interpolation for sin/cos over [0, 2pi)
static inline void sincos_cardinal_bspline(float theta, float& c_out, float& s_out) {
    const int n = g_trig_table.n_grid;
    const float inv_twopi = 1.0f / (2.0f * float(M_PI));
    float u = theta * inv_twopi * float(n);
    float uf = floorf(u);
    int i = int(uf);
    float t = u - uf;

    float om = 1.0f - t;
    float t2 = t * t;
    float t3 = t2 * t;
    float om3 = om * om * om;
    // Cardinal cubic B-spline basis
    float w0 = om3 / 6.0f;
    float w1 = (3.0f * t3 - 6.0f * t2 + 4.0f) / 6.0f;
    float w2 = (-3.0f * t3 + 3.0f * t2 + 3.0f * t + 1.0f) / 6.0f;
    float w3 = t3 / 6.0f;

    int i0 = wrap_idx(i - 1, n);
    int i1 = wrap_idx(i, n);
    int i2 = wrap_idx(i + 1, n);
    int i3 = wrap_idx(i + 2, n);

    c_out = w0 * g_trig_table.cos_table[i0] + w1 * g_trig_table.cos_table[i1] +
            w2 * g_trig_table.cos_table[i2] + w3 * g_trig_table.cos_table[i3];
    s_out = w0 * g_trig_table.sin_table[i0] + w1 * g_trig_table.sin_table[i1] +
            w2 * g_trig_table.sin_table[i2] + w3 * g_trig_table.sin_table[i3];
}

// Build k-vectors and prefactors for current box dimensions
static void build_kvectors(EwaldState& st) {
    auto& s = st.settings;
    float bx = st.box_x, by = st.box_y, bz = st.box_z;
    if(bx <= 0.f || by <= 0.f || bz <= 0.f) return;

    float V = bx * by * bz;
    float alpha = s.alpha;
    float four_alpha2 = 4.0f * alpha * alpha;
    int kmax = s.kmax;

    st.kx.clear(); st.ky.clear(); st.kz.clear();
    st.k2.clear(); st.k_prefactor.clear();

    float twopi_bx = 2.0f * M_PI / bx;
    float twopi_by = 2.0f * M_PI / by;
    float twopi_bz = 2.0f * M_PI / bz;

    for(int nx = -kmax; nx <= kmax; ++nx) {
        for(int ny = -kmax; ny <= kmax; ++ny) {
            for(int nz = -kmax; nz <= kmax; ++nz) {
                if(nx == 0 && ny == 0 && nz == 0) continue;
                float gx = nx * twopi_bx;
                float gy = ny * twopi_by;
                float gz = nz * twopi_bz;
                float g2 = gx*gx + gy*gy + gz*gz;

                float prefactor = (4.0f * M_PI / V) * expf(-g2 / four_alpha2) / g2;

                st.kx.push_back(gx);
                st.ky.push_back(gy);
                st.kz.push_back(gz);
                st.k2.push_back(g2);
                st.k_prefactor.push_back(prefactor);
            }
        }
    }

    // Self-energy correction: -alpha/sqrt(pi) * sum_i q_i^2
    double q2_sum = 0.0;
    for(int i = 0; i < st.n_atom; ++i) {
        q2_sum += double(st.charges[i]) * double(st.charges[i]);
    }
    st.self_energy = -s.coulomb_k * float(alpha / sqrtf(M_PI) * q2_sum);
}

void initialize_ewald(hid_t config_root, DerivEngine& engine) {
    EwaldSettings s;
    try {
        if(h5_exists(config_root, "/input/potential/martini_potential")) {
            auto grp = open_group(config_root, "/input/potential/martini_potential");
            if(attribute_exists_ewald(grp.get(), ".", "ewald_enabled"))
                s.enabled = read_attribute<int>(grp.get(), ".", "ewald_enabled") != 0;
            if(!s.enabled) return;
            if(attribute_exists_ewald(grp.get(), ".", "ewald_alpha"))
                s.alpha = read_attribute<float>(grp.get(), ".", "ewald_alpha");
            if(attribute_exists_ewald(grp.get(), ".", "ewald_kmax"))
                s.kmax = read_attribute<int>(grp.get(), ".", "ewald_kmax");
            if(attribute_exists_ewald(grp.get(), ".", "ewald_use_cardinal_bspline"))
                s.use_cardinal_bspline = read_attribute<int>(grp.get(), ".", "ewald_use_cardinal_bspline") != 0;
            if(attribute_exists_ewald(grp.get(), ".", "ewald_bspline_grid"))
                s.bspline_grid = read_attribute<int>(grp.get(), ".", "ewald_bspline_grid");
        } else {
            return;
        }
    } catch(...) {
        return;
    }

    // Read box dimensions
    float bx = 0.f, by = 0.f, bz = 0.f;
    try {
        if(h5_exists(config_root, "/input/potential/martini_potential")) {
            auto grp = open_group(config_root, "/input/potential/martini_potential");
            bx = read_attribute<float>(grp.get(), ".", "x_len");
            by = read_attribute<float>(grp.get(), ".", "y_len");
            bz = read_attribute<float>(grp.get(), ".", "z_len");
        }
    } catch(...) {
        bx = by = bz = 0.f;
    }

    // Read per-atom charges
    int n_atom = engine.pos->n_elem;
    std::vector<float> charges(n_atom, 0.f);
    try {
        if(h5_exists(config_root, "/input/potential/martini_potential")) {
            auto grp = open_group(config_root, "/input/potential/martini_potential");
            traverse_dset<1,float>(grp.get(), "charges", [&](size_t i, float q) {
                if((int)i < n_atom) charges[i] = q;
            });
        }
    } catch(...) {}

    std::lock_guard<std::mutex> lk(g_ewald_mutex);
    auto& st = g_ewald_state[&engine];
    st.settings = s;
    st.box_x = bx; st.box_y = by; st.box_z = bz;
    st.n_atom = n_atom;
    st.charges = std::move(charges);

    build_kvectors(st);
    if(s.use_cardinal_bspline) {
        ensure_trig_table(s.bspline_grid);
    }

    printf("[EWALD] Initialized: alpha=%.4f kmax=%d box=%.2f x %.2f x %.2f n_kvec=%d self_energy=%.6f bspline=%d grid=%d\n",
           s.alpha, s.kmax, bx, by, bz, (int)st.kx.size(), st.self_energy, s.use_cardinal_bspline ? 1 : 0, s.bspline_grid);
}

void update_kvectors(DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_ewald_mutex);
    auto it = g_ewald_state.find(&engine);
    if(it == g_ewald_state.end()) return;
    auto& st = it->second;
    if(!st.settings.enabled) return;

    // Sync box from NPT barostat
    float bx, by, bz;
    npt::get_current_box(engine, bx, by, bz);
    if(bx > 0.f && by > 0.f && bz > 0.f) {
        st.box_x = bx; st.box_y = by; st.box_z = bz;
        build_kvectors(st);
    }
}

void compute_ewald_reciprocal(DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_ewald_mutex);
    auto it = g_ewald_state.find(&engine);
    if(it == g_ewald_state.end()) return;
    auto& st = it->second;
    if(!st.settings.enabled) return;
    if(st.kx.empty()) return;

    int n_atom = st.n_atom;
    int n_kvec = (int)st.kx.size();
    float coulomb_k = st.settings.coulomb_k;

    VecArray pos = engine.pos->output;
    VecArray sens = engine.pos->sens;

    // For each k-vector, compute structure factor S(k) = sum_i q_i * exp(i k.r_i)
    // Energy = coulomb_k * sum_k prefactor * |S(k)|^2  +  self_energy
    double E_recip = 0.0;

    // Temporary per-atom force accumulators
    std::vector<double> fx(n_atom, 0.0), fy(n_atom, 0.0), fz(n_atom, 0.0);
    std::vector<float> cos_cache(n_atom, 0.0f), sin_cache(n_atom, 0.0f);
    const bool use_bspline = st.settings.use_cardinal_bspline;

    for(int m = 0; m < n_kvec; ++m) {
        float gx = st.kx[m], gy = st.ky[m], gz = st.kz[m];
        float pref = st.k_prefactor[m];

        // Compute structure factor components
        double S_re = 0.0, S_im = 0.0;
        for(int i = 0; i < n_atom; ++i) {
            float qi = st.charges[i];
            if(qi == 0.f) continue;
            float kr = gx * pos(0,i) + gy * pos(1,i) + gz * pos(2,i);
            float c = 0.0f, s = 0.0f;
            if(use_bspline) {
                sincos_cardinal_bspline(kr, c, s);
            } else {
                c = cosf(kr);
                s = sinf(kr);
            }
            cos_cache[i] = c;
            sin_cache[i] = s;
            S_re += qi * c;
            S_im += qi * s;
        }

        double S2 = S_re * S_re + S_im * S_im;
        E_recip += pref * S2;

        // Force on atom i: F_i = -dE/dr_i
        // dE/dk = 2 * pref * (S_re * d(S_re)/dr_i + S_im * d(S_im)/dr_i)
        // d(S_re)/dr_i = -q_i * sin(k.r_i) * k
        // d(S_im)/dr_i =  q_i * cos(k.r_i) * k
        // => dE/dr_i = 2 * pref * q_i * (S_im * cos(k.r_i) - S_re * sin(k.r_i)) * k
        // Force = -dE/dr_i
        for(int i = 0; i < n_atom; ++i) {
            float qi = st.charges[i];
            if(qi == 0.f) continue;
            float cos_kr = cos_cache[i];
            float sin_kr = sin_cache[i];
            // gradient component = 2 * pref * qi * (S_im*cos - S_re*sin)
            double grad_scalar = 2.0 * pref * qi * (S_im * cos_kr - S_re * sin_kr);
            // This is dE/d(component), force = -gradient
            // But sens stores gradient (not force), so we ADD grad_scalar * k
            fx[i] += grad_scalar * gx;
            fy[i] += grad_scalar * gy;
            fz[i] += grad_scalar * gz;
        }
    }

    // Scale by coulomb_k and add to sens (gradient accumulator)
    // sens stores -force = gradient of potential
    for(int i = 0; i < n_atom; ++i) {
        sens(0, i) += float(coulomb_k * fx[i]);
        sens(1, i) += float(coulomb_k * fy[i]);
        sens(2, i) += float(coulomb_k * fz[i]);
    }

    st.reciprocal_energy = float(coulomb_k * E_recip) + st.self_energy;

    // Add reciprocal energy to engine potential
    engine.potential += st.reciprocal_energy;
}

bool is_enabled(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_ewald_mutex);
    auto it = g_ewald_state.find(const_cast<DerivEngine*>(&engine));
    return it != g_ewald_state.end() && it->second.settings.enabled;
}

float get_reciprocal_energy(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lk(g_ewald_mutex);
    auto it = g_ewald_state.find(const_cast<DerivEngine*>(&engine));
    if(it != g_ewald_state.end()) return it->second.reciprocal_energy;
    return 0.f;
}

} // namespace ewald

} // namespace simulation_box
