#pragma once
#include "deriv_engine.h"
#include "h5_support.h"
#include <cmath>
#include <vector>
#include <random>

namespace simulation_box {

inline Vec<3,float> minimum_image(const Vec<3,float>& dr, float box_x, float box_y, float box_z) {
    Vec<3,float> out = dr;
    if(box_x > 0.f) out.x() -= box_x * roundf(out.x() / box_x);
    if(box_y > 0.f) out.y() -= box_y * roundf(out.y() / box_y);
    if(box_z > 0.f) out.z() -= box_z * roundf(out.z() / box_z);
    return out;
}

inline void minimum_image_scalar(float& dx, float& dy, float& dz, float box_x, float box_y, float box_z) {
    if(box_x > 0.f) dx -= box_x * roundf(dx / box_x);
    if(box_y > 0.f) dy -= box_y * roundf(dy / box_y);
    if(box_z > 0.f) dz -= box_z * roundf(dz / box_z);
}

namespace npt {

enum class BarostatType {
    Berendsen,
    ParrinelloRahman,
    MonteCarlo
};

struct BarostatSettings {
    bool  semi_isotropic = true;
    BarostatType type = BarostatType::Berendsen;
    int   interval = 0;           // steps between applications
    float target_p_xy = 1.0f;     // UP units
    float target_p_z  = 1.0f;     // UP units
    float tau_p = 5.0f;           // time constant
    float compressibility_xy = 4.5e-5f; // 1/pressure in bilayer plane
    float compressibility_z = 4.5e-5f;  // 1/pressure normal to bilayer
    float compressibility = 4.5e-5f;
    bool  prefer_shrink_first = true; // on first application, avoid any expansion
    // Monte-Carlo barostat: max fractional box-length step per trial move
    float    mc_dmax_xy = 0.004f;
    float    mc_dmax_z  = 0.004f;
    unsigned mc_seed    = 20240719u;
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
    // Monte-Carlo barostat: molecule grouping (COM scaling) + RNG + acceptance stats
    std::vector<std::vector<int>> mc_mol_atoms;  // atoms of each scaled molecule
    std::mt19937_64 mc_rng;
    bool mc_init = false;
    long mc_attempt[2] = {0, 0};   // [0]=lateral (xy), [1]=normal (z)
    long mc_accept[2]  = {0, 0};
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
                          float temperature,
                          int verbose,
                          bool print_now);

void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz);

void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z);

bool is_enabled(const DerivEngine& engine);

void get_pressure(const DerivEngine& engine, float& pxy, float& pz);

float get_volume(const DerivEngine& engine);

} // namespace npt

} // namespace simulation_box
