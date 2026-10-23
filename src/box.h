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

struct BarostatSettings {
    int   interval;      // rounds between applications
    float target_p_xy;   // E_up/A^3
    float target_p_z;    // E_up/A^3
    // Monte-Carlo barostat: largest log box-length step of each trial move; 0 holds that axis
    float mc_dmax_xy;
    float mc_dmax_z;
};

struct BarostatState {
    BarostatSettings settings;
    float box_x = 0.f, box_y = 0.f, box_z = 0.f;
    // Monte-Carlo barostat: molecule grouping (COM scaling) + RNG
    std::vector<std::vector<int>> mc_mol_atoms;  // atoms of each scaled molecule
    std::mt19937_64 mc_rng;
    uint32_t mc_seed = 0;
};

// Sets the periodic box of every node that holds one.
using NodeBoxUpdater = void (*)(DerivEngine& engine, float box_x, float box_y, float box_z);

void register_barostat_for_engine(hid_t config_root, DerivEngine& engine, uint32_t random_seed);
void unregister_barostat_for_engine(DerivEngine& engine);

void register_node_box_updater(NodeBoxUpdater updater);

void maybe_apply_barostat(DerivEngine& engine,
                          int n_atom,
                          uint64_t round_num,
                          float temperature);

void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz);

void update_node_boxes(DerivEngine& engine, float box_x, float box_y, float box_z);

bool is_enabled(const DerivEngine& engine);

float get_volume(const DerivEngine& engine);

} // namespace npt

} // namespace simulation_box
