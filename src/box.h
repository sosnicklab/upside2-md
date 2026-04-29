#pragma once
#include "deriv_engine.h"
#include "h5_support.h"
#include <vector>
#include <map>
#include <mutex>

namespace simulation_box {

inline Vec<3,float> minimum_image(const Vec<3,float>& dr, float box_x, float box_y, float box_z) {
    Vec<3,float> out = dr;
    if (out.x() >  0.5f * box_x) out.x() -= box_x;
    if (out.x() < -0.5f * box_x) out.x() += box_x;
    if (out.y() >  0.5f * box_y) out.y() -= box_y;
    if (out.y() < -0.5f * box_y) out.y() += box_y;
    if (out.z() >  0.5f * box_z) out.z() -= box_z;
    if (out.z() < -0.5f * box_z) out.z() += box_z;
    return out;
}

inline void minimum_image_scalar(float& dx, float& dy, float& dz, float box_x, float box_y, float box_z) {
    if (dx >  0.5f * box_x) dx -= box_x;
    if (dx < -0.5f * box_x) dx += box_x;
    if (dy >  0.5f * box_y) dy -= box_y;
    if (dy < -0.5f * box_y) dy += box_y;
    if (dz >  0.5f * box_z) dz -= box_z;
    if (dz < -0.5f * box_z) dz += box_z;
}

void wrap_positions(VecArray pos, int n_atom, float box_x, float box_y, float box_z);

namespace npt {

enum class BarostatType {
    Berendsen,
    ParrinelloRahman
};

struct BarostatSettings {
    bool  enabled = false;
    bool  semi_isotropic = true;
    BarostatType type = BarostatType::Berendsen;
    int   interval = 0;           // steps between applications
    float target_p_xy = 1.0f;     // UP units
    float target_p_z  = 1.0f;     // UP units
    float tau_p = 5.0f;           // time constant
    float compressibility_xy = 4.5e-5f; // 1/pressure in bilayer plane
    float compressibility_z = 4.5e-5f;  // 1/pressure normal to bilayer
    float compressibility = 4.5e-5f;    // legacy fallback (kept for compatibility)
    bool  prefer_shrink_first = true; // on first application, avoid any expansion
};

struct BarostatState {
    BarostatSettings settings;
    float box_x = 0.f, box_y = 0.f, box_z = 0.f;
    std::vector<float> masses;
    bool has_applied_once = false;
    float last_pxy_inst = 0.0f;
    float last_pz_inst = 0.0f;
    float last_scale_xy = 1.0f;
    float last_scale_z = 1.0f;
    float box_vel_xy = 0.0f;  // box velocity for lateral dimensions
    float box_vel_z = 0.0f;   // box velocity for normal dimension
};

using NodeBoxUpdater = void (*)(DerivEngine& engine, float scale_xy, float scale_z);

void register_barostat_for_engine(hid_t config_root, DerivEngine& engine);

void register_node_box_updater(NodeBoxUpdater updater);

void maybe_apply_barostat(DerivEngine& engine,
                          const VecArray& mom,
                          int n_atom,
                          uint64_t round_num,
                          float dt,
                          int inner_step,
                          int verbose,
                          bool print_now);

void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz);

void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z);

bool is_enabled(const DerivEngine& engine);

void get_pressure(const DerivEngine& engine, float& pxy, float& pz);

float get_volume(const DerivEngine& engine);

} // namespace npt

namespace ewald {

struct EwaldSettings {
    bool  enabled = false;
    float alpha = 0.2f;          // Ewald screening parameter (1/Angstrom)
    int   kmax = 5;              // k-space cutoff (number of k-vectors per dimension)
    float coulomb_k = 31.775347952181f; // Coulomb constant (includes epsilon_r=15)
    bool  use_cardinal_bspline = true;  // approximate trig via periodic cubic B-spline
    int   bspline_grid = 16384;         // lookup grid size for periodic trig table
};

struct EwaldState {
    EwaldSettings settings;
    float box_x = 0.f, box_y = 0.f, box_z = 0.f;
    int n_atom = 0;
    std::vector<float> charges;
    std::vector<float> kx, ky, kz;
    std::vector<float> k2;
    std::vector<float> k_prefactor;
    float self_energy = 0.f;          // Ewald self-energy correction (constant)
    float reciprocal_energy = 0.f;    // last computed reciprocal energy
};

void initialize_ewald(hid_t config_root, DerivEngine& engine);

void update_kvectors(DerivEngine& engine);

void compute_ewald_reciprocal(DerivEngine& engine);

bool is_enabled(const DerivEngine& engine);

float get_reciprocal_energy(const DerivEngine& engine);

} // namespace ewald

} // namespace simulation_box
