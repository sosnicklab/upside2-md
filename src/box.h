#pragma once

#include "deriv_engine.h"
#include "h5_support.h"
#include <vector>

namespace simulation_box {

inline Vec<3,float> minimum_image(const Vec<3,float>& dr, float box_x, float box_y, float box_z) {
    Vec<3,float> out = dr;
    if(out.x() >  0.5f * box_x) out.x() -= box_x;
    if(out.x() < -0.5f * box_x) out.x() += box_x;
    if(out.y() >  0.5f * box_y) out.y() -= box_y;
    if(out.y() < -0.5f * box_y) out.y() += box_y;
    if(out.z() >  0.5f * box_z) out.z() -= box_z;
    if(out.z() < -0.5f * box_z) out.z() += box_z;
    return out;
}

inline void minimum_image_scalar(float& dx, float& dy, float& dz, float box_x, float box_y, float box_z) {
    if(dx >  0.5f * box_x) dx -= box_x;
    if(dx < -0.5f * box_x) dx += box_x;
    if(dy >  0.5f * box_y) dy -= box_y;
    if(dy < -0.5f * box_y) dy += box_y;
    if(dz >  0.5f * box_z) dz -= box_z;
    if(dz < -0.5f * box_z) dz += box_z;
}

namespace npt {

struct BarostatSettings {
    bool enabled;
    bool semi_isotropic;
    int interval;
    float target_p_xy;
    float target_p_z;
    float tau_p;
    float compressibility_xy;
    float compressibility_z;
    bool debug;
    bool prefer_shrink_first;

    BarostatSettings():
        enabled(false),
        semi_isotropic(true),
        interval(0),
        target_p_xy(1.0f),
        target_p_z(1.0f),
        tau_p(5.0f),
        compressibility_xy(4.5e-5f),
        compressibility_z(4.5e-5f),
        debug(true),
        prefer_shrink_first(true) {}
};

struct BarostatState {
    enum { EQUILIBRIUM_THRESHOLD = 5 };

    BarostatSettings settings;
    float box_x;
    float box_y;
    float box_z;
    std::vector<float> masses;
    bool has_applied_once;
    int equilibrium_count;
    float last_pxy_inst;
    float last_pz_inst;

    BarostatState():
        box_x(0.f),
        box_y(0.f),
        box_z(0.f),
        has_applied_once(false),
        equilibrium_count(0),
        last_pxy_inst(0.f),
        last_pz_inst(0.f) {}
};

typedef void (*NodeBoxUpdater)(DerivEngine& engine, float scale_xy, float scale_z);

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

}

namespace ewald {

struct EwaldSettings {
    bool enabled;
    float alpha;
    int kmax;
    float coulomb_k;
    bool use_cardinal_bspline;
    int bspline_grid;

    EwaldSettings():
        enabled(false),
        alpha(0.2f),
        kmax(5),
        coulomb_k(31.775347952181f),
        use_cardinal_bspline(true),
        bspline_grid(16384) {}
};

struct EwaldState {
    EwaldSettings settings;
    float box_x;
    float box_y;
    float box_z;
    int n_atom;
    std::vector<float> charges;
    std::vector<float> kx;
    std::vector<float> ky;
    std::vector<float> kz;
    std::vector<float> k2;
    std::vector<float> k_prefactor;
    float self_energy;
    float reciprocal_energy;

    EwaldState():
        box_x(0.f),
        box_y(0.f),
        box_z(0.f),
        n_atom(0),
        self_energy(0.f),
        reciprocal_energy(0.f) {}
};

void initialize_ewald(hid_t config_root, DerivEngine& engine);
void update_kvectors(DerivEngine& engine);
void compute_ewald_reciprocal(DerivEngine& engine);
bool is_enabled(const DerivEngine& engine);

}

}
