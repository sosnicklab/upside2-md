#pragma once
// box.h - Simulation box management for periodic boundaries and NPT/NVT ensembles
// Ported from origin/martini branch

#include "deriv_engine.h"
#include "h5_support.h"
#include <vector>
#include <map>
#include <mutex>

namespace simulation_box {

// ===================== PBC MINIMUM IMAGE =====================
// Helper for minimum image convention in a rectangular box
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

// Minimum image for scalar components (used by dihedral/Float4 calculations)
inline void minimum_image_scalar(float& dx, float& dy, float& dz, float box_x, float box_y, float box_z) {
    if (dx >  0.5f * box_x) dx -= box_x;
    if (dx < -0.5f * box_x) dx += box_x;
    if (dy >  0.5f * box_y) dy -= box_y;
    if (dy < -0.5f * box_y) dy += box_y;
    if (dz >  0.5f * box_z) dz -= box_z;
    if (dz < -0.5f * box_z) dz += box_z;
}

// Wrap all particle positions into the primary box [0, box_dim)
void wrap_positions(VecArray pos, int n_atom, float box_x, float box_y, float box_z);

// ===================== NPT BAROSTAT =====================
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
    bool  debug = true;
    bool  prefer_shrink_first = true; // on first application, avoid any expansion
};

struct BarostatState {
    BarostatSettings settings;
    float box_x = 0.f, box_y = 0.f, box_z = 0.f;
    std::vector<float> masses;
    bool has_applied_once = false;
    // Equilibrium detection
    float prev_box_x = 0.f, prev_box_y = 0.f, prev_box_z = 0.f;
    int equilibrium_count = 0;
    static constexpr int EQUILIBRIUM_THRESHOLD = 5;
    static constexpr float EQUILIBRIUM_TOLERANCE = 0.001f;
    // Last applied values for logging
    float last_pxy_inst = 0.0f;
    float last_pz_inst = 0.0f;
    float last_scale_xy = 1.0f;
    float last_scale_z = 1.0f;
    // Parrinello-Rahman state variables
    float box_vel_xy = 0.0f;  // box velocity for lateral dimensions
    float box_vel_z = 0.0f;   // box velocity for normal dimension
};

// Initialize barostat for an engine (reads settings from H5)
void register_barostat_for_engine(hid_t config_root, DerivEngine& engine);

// Apply barostat at configured interval
void maybe_apply_barostat(DerivEngine& engine,
                          const VecArray& mom,
                          int n_atom,
                          uint64_t round_num,
                          float dt,
                          int inner_step,
                          int verbose,
                          bool print_now);

// Get current box dimensions for an engine
void get_current_box(const DerivEngine& engine, float& bx, float& by, float& bz);

// Update box dimensions in MARTINI potential nodes after NPT scaling
void update_node_boxes(DerivEngine& engine, float scale_xy, float scale_z);

// Check if barostat is enabled for this engine
bool is_enabled(const DerivEngine& engine);

// Get last calculated pressure values
void get_pressure(const DerivEngine& engine, float& pxy, float& pz);

// Get current box volume
float get_volume(const DerivEngine& engine);

} // namespace npt

// ===================== EWALD SUMMATION =====================
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
    std::vector<float> charges;       // per-atom charges
    // Pre-computed k-vectors and their |k|^2
    std::vector<float> kx, ky, kz;    // k-vector components
    std::vector<float> k2;            // |k|^2
    std::vector<float> k_prefactor;   // (4*pi/V) * exp(-k2/(4*alpha^2)) / k2
    float self_energy = 0.f;          // Ewald self-energy correction (constant)
    float reciprocal_energy = 0.f;    // last computed reciprocal energy
};

// Initialize Ewald state for an engine (reads charges + settings from H5)
void initialize_ewald(hid_t config_root, DerivEngine& engine);

// Recompute k-vector prefactors after box change (NPT)
void update_kvectors(DerivEngine& engine);

// Compute reciprocal-space energy and forces; adds to pos.sens
void compute_ewald_reciprocal(DerivEngine& engine);

// Query whether Ewald is enabled for this engine
bool is_enabled(const DerivEngine& engine);

// Get last computed reciprocal-space energy
float get_reciprocal_energy(const DerivEngine& engine);

} // namespace ewald

} // namespace simulation_box
