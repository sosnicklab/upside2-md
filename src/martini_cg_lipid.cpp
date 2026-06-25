#include "deriv_engine.h"
#include "spline.h"
#include "box.h"
#include "martini.h"
#include "random.h"
#include "state_logger.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace h5;

namespace {

const int CG_LIPID_N_ANGULAR = 15;
const int CG_LIPID_N_RADIAL = 12;
const int CG_LIPID_N_PARAM = 54;
const float CG_LIPID_DEFAULT_KNOT_SPACING = 0.7f;
const float CG_LIPID_INV_DTHETA = 6.f;

struct QuadsplineEval {
    float value;
    float d_dr;
    float d_da1;
    float d_da2;
};

struct TargetSplineEval {
    float value;
    float d_dr;
    float d_da;
};

struct CGLPairIndex {
    size_t first;
    size_t second;
};

struct CGLBodySupport {
    bool has_geometry;
    float axis_radius;
    float perp_radius;
    float bead_cutoff;
};

static vector<int> read_int_dataset(hid_t grp, const char* name) {
    vector<hsize_t> sz = get_dset_size(1, grp, name);
    vector<int> out(sz[0]);
    traverse_dset<1, int32_t>(grp, name, [&](size_t i, int32_t v) {
        out[i] = v;
    });
    return out;
}

static vector<float> read_float_dataset(hid_t grp, const char* name) {
    vector<hsize_t> sz = get_dset_size(2, grp, name);
    vector<float> out(sz[0] * sz[1], 0.f);
    traverse_dset<2, float>(grp, name, [&](size_t i, size_t j, float v) {
        out[i * sz[1] + j] = v;
    });
    return out;
}

static vector<float> read_float_dataset_1d(hid_t grp, const char* name) {
    vector<hsize_t> sz = get_dset_size(1, grp, name);
    vector<float> out(sz[0], 0.f);
    traverse_dset<1, float>(grp, name, [&](size_t i, float v) {
        out[i] = v;
    });
    return out;
}

static vector<float> read_param_dataset(hid_t grp, int& n_type1, int& n_type2) {
    vector<hsize_t> sz = get_dset_size(3, grp, "interaction_param");
    if(sz[2] != CG_LIPID_N_PARAM)
        throw string("cg lipid interaction_param must have last dimension 54");
    n_type1 = int(sz[0]);
    n_type2 = int(sz[1]);
    vector<float> out(n_type1 * n_type2 * CG_LIPID_N_PARAM, 0.f);
    traverse_dset<3, float>(grp, "interaction_param",
            [&](size_t i, size_t j, size_t k, float v) {
        out[(i * n_type2 + j) * CG_LIPID_N_PARAM + k] = v;
    });
    return out;
}

static vector<float> read_param_dataset_any(
        hid_t grp, int& n_type1, int& n_type2, int& n_param) {
    vector<hsize_t> sz = get_dset_size(3, grp, "interaction_param");
    n_type1 = int(sz[0]);
    n_type2 = int(sz[1]);
    n_param = int(sz[2]);
    vector<float> out(n_type1 * n_type2 * n_param, 0.f);
    traverse_dset<3, float>(grp, "interaction_param",
            [&](size_t i, size_t j, size_t k, float v) {
        out[(i * n_type2 + j) * n_param + k] = v;
    });
    return out;
}

static inline const float* param_ptr(
        const vector<float>& param, int n_type2, int type1, int type2) {
    return &param[(type1 * n_type2 + type2) * CG_LIPID_N_PARAM];
}

static inline const float* param_ptr(
        const vector<float>& param, int n_type2, int n_param, int type1, int type2) {
    return &param[(type1 * n_type2 + type2) * n_param];
}

static inline float clamp_angle(float x) {
    return std::max(-1.f, std::min(1.f, x));
}

static inline float angular_spline_coord(float cos_angle) {
    float t = (clamp_angle(cos_angle) + 1.f) * CG_LIPID_INV_DTHETA + 1.f;
    return std::max(1.0001f, std::min(float(CG_LIPID_N_ANGULAR - 2) - 0.0001f, t));
}

static inline Vec<2> radial_deBoor_value_and_deriv(
        const float* bspline_coeff, const float x, int n_knot) {
    if(x < 1.f) {
        constexpr float x0 = 1.0001f;
        Vec<2> v = deBoor_value_and_deriv(bspline_coeff, x0);
        return make_vec2(v.x() + (x - x0) * v.y(), v.y());
    }
    return clamped_deBoor_value_and_deriv(bspline_coeff, x, n_knot);
}

static bool eval_quadspline(
        const float* p,
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float knot_spacing,
        float cutoff,
        float taper_width,
        QuadsplineEval& out) {

    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f) return false;

    float r = sqrtf(r2);
    if(r >= cutoff) return false;

    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};

    float a1 = -(n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2]);
    float a2 =  (n2[0] * unit[0] + n2[1] * unit[1] + n2[2] * unit[2]);

    Vec<2> ang1 = deBoor_value_and_deriv(p, angular_spline_coord(a1));
    Vec<2> ang2 = deBoor_value_and_deriv(p + CG_LIPID_N_ANGULAR, angular_spline_coord(a2));
    float radial_coord = r / knot_spacing;
    Vec<2> vr = radial_deBoor_value_and_deriv(
            p + 2 * CG_LIPID_N_ANGULAR, radial_coord, CG_LIPID_N_RADIAL);
    Vec<2> va = radial_deBoor_value_and_deriv(
            p + 2 * CG_LIPID_N_ANGULAR + CG_LIPID_N_RADIAL,
            radial_coord, CG_LIPID_N_RADIAL);

    float raw_value = vr.x() + ang1.x() * ang2.x() * va.x();
    float raw_d_dr = (vr.y() + ang1.x() * ang2.x() * va.y()) / knot_spacing;
    float raw_d_da1 = ang1.y() * CG_LIPID_INV_DTHETA * ang2.x() * va.x();
    float raw_d_da2 = ang2.y() * CG_LIPID_INV_DTHETA * ang1.x() * va.x();

    float taper = 1.f;
    float d_taper_dr = 0.f;
    taper_width = std::max(taper_width, 1e-6f);
    float taper_start = cutoff - taper_width;
    if(r > taper_start) {
        float u = (cutoff - r) / taper_width;
        u = std::max(0.f, std::min(1.f, u));
        taper = u * u * (3.f - 2.f * u);
        d_taper_dr = -6.f * u * (1.f - u) / taper_width;
    }

    out.value = taper * raw_value;
    out.d_dr = taper * raw_d_dr + raw_value * d_taper_dr;
    out.d_da1 = taper * raw_d_da1;
    out.d_da2 = taper * raw_d_da2;
    return true;
}

static inline float angular_spline_coord(float cos_angle, int n_angular) {
    float inv_dtheta = float(n_angular - 3) * 0.5f;
    float t = (clamp_angle(cos_angle) + 1.f) * inv_dtheta + 1.f;
    return std::max(1.0001f, std::min(float(n_angular - 2) - 0.0001f, t));
}

static bool eval_multimode_pair(
        const float* p,
        int n_modes,
        int n_angular,
        int n_radial,
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float knot_spacing,
        float cutoff,
        float taper_width,
        QuadsplineEval& out) {

    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f) return false;

    float r = sqrtf(r2);
    if(r >= cutoff) return false;

    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};

    float a1 = -(n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2]);
    float a2 =  (n2[0] * unit[0] + n2[1] * unit[1] + n2[2] * unit[2]);
    float a1_coord = angular_spline_coord(a1, n_angular);
    float a2_coord = angular_spline_coord(a2, n_angular);
    float inv_dtheta = float(n_angular - 3) * 0.5f;
    float radial_coord = r / knot_spacing;

    Vec<2> v0 = radial_deBoor_value_and_deriv(p, radial_coord, n_radial);
    float raw_value = v0.x();
    float raw_d_dr = v0.y() / knot_spacing;
    float raw_d_da1 = 0.f;
    float raw_d_da2 = 0.f;

    const float* mode = p + n_radial;
    for(int m = 0; m < n_modes; ++m) {
        const float* a1_ptr = mode;
        const float* a2_ptr = mode + n_angular;
        const float* v_ptr = mode + 2 * n_angular;
        Vec<2> ang1 = deBoor_value_and_deriv(a1_ptr, a1_coord);
        Vec<2> ang2 = deBoor_value_and_deriv(a2_ptr, a2_coord);
        Vec<2> vm = radial_deBoor_value_and_deriv(v_ptr, radial_coord, n_radial);
        raw_value += ang1.x() * ang2.x() * vm.x();
        raw_d_dr += ang1.x() * ang2.x() * vm.y() / knot_spacing;
        raw_d_da1 += ang1.y() * inv_dtheta * ang2.x() * vm.x();
        raw_d_da2 += ang2.y() * inv_dtheta * ang1.x() * vm.x();
        mode += 2 * n_angular + n_radial;
    }

    float taper = 1.f;
    float d_taper_dr = 0.f;
    taper_width = std::max(taper_width, 1e-6f);
    float taper_start = cutoff - taper_width;
    if(r > taper_start) {
        float u = (cutoff - r) / taper_width;
        u = std::max(0.f, std::min(1.f, u));
        taper = u * u * (3.f - 2.f * u);
        d_taper_dr = -6.f * u * (1.f - u) / taper_width;
    }

    out.value = taper * raw_value;
    out.d_dr = taper * raw_d_dr + raw_value * d_taper_dr;
    out.d_da1 = taper * raw_d_da1;
    out.d_da2 = taper * raw_d_da2;
    return true;
}

static void clamped_uniform_bspline_basis(
        float x,
        int n_knot,
        float w[4],
        float dw[4],
        int& base) {
    for(int i = 0; i < 4; ++i) {
        w[i] = 0.f;
        dw[i] = 0.f;
    }
    if(x <= 1.f) {
        base = 0;
        w[0] = 1.f / 6.f;
        w[1] = 2.f / 3.f;
        w[2] = 1.f / 6.f;
        return;
    }
    if(x >= float(n_knot - 2)) {
        base = n_knot - 3;
        w[0] = 1.f / 6.f;
        w[1] = 2.f / 3.f;
        w[2] = 1.f / 6.f;
        return;
    }

    int x_bin = int(x);
    base = x_bin - 1;
    float excess = x - float(x_bin);

    for(int ci = 0; ci < 4; ++ci) {
        float c0[4] = {0.f, 0.f, 0.f, 0.f};
        c0[ci] = 1.f;
        const float u[4] = {-2.f, -1.f, 0.f, 0.f};

        float yu1 = excess - u[0];
        float yu2 = excess - u[1];
        float yu3 = excess - u[2];

        float a11 = yu1 / 3.f;
        float a12 = yu2 / 3.f;
        float a13 = yu3 / 3.f;
        float c11 = (1.f - a11) * c0[0] + a11 * c0[1];
        float c12 = (1.f - a12) * c0[1] + a12 * c0[2];
        float c13 = (1.f - a13) * c0[2] + a13 * c0[3];

        float d11 = c0[1] - c0[0];
        float d12 = c0[2] - c0[1];
        float d13 = c0[3] - c0[2];

        float a22 = yu2 / 2.f;
        float a23 = yu3 / 2.f;
        float c22 = (1.f - a22) * c11 + a22 * c12;
        float c23 = (1.f - a23) * c12 + a23 * c13;
        float d22 = (1.f - a22) * d11 + a22 * d12;
        float d23 = (1.f - a23) * d12 + a23 * d13;

        float a33 = yu3;
        w[ci] = (1.f - a33) * c22 + a33 * c23;
        dw[ci] = (1.f - a33) * d22 + a33 * d23;
    }
}

static void radial_uniform_bspline_basis(
        float x,
        int n_knot,
        float w[4],
        float dw[4],
        int& base) {
    if(x < 1.f) {
        constexpr float x0 = 1.0001f;
        clamped_uniform_bspline_basis(x0, n_knot, w, dw, base);
        float dx = x - x0;
        for(int i = 0; i < 4; ++i)
            w[i] += dx * dw[i];
        return;
    }
    clamped_uniform_bspline_basis(x, n_knot, w, dw, base);
}

static bool eval_full_pair_tensor(
        const float* p,
        int n_angular,
        int n_radial,
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float knot_spacing,
        float cutoff,
        float taper_width,
        QuadsplineEval& out) {

    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f) return false;
    float r = sqrtf(r2);
    if(r >= cutoff) return false;

    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};
    float a1 = -(n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2]);
    float a2 =  (n2[0] * unit[0] + n2[1] * unit[1] + n2[2] * unit[2]);

    float wr[4], dwr[4], wa1[4], dwa1[4], wa2[4], dwa2[4];
    int br = 0, ba1 = 0, ba2 = 0;
    radial_uniform_bspline_basis(r / knot_spacing, n_radial, wr, dwr, br);
    clamped_uniform_bspline_basis(angular_spline_coord(a1, n_angular), n_angular, wa1, dwa1, ba1);
    clamped_uniform_bspline_basis(angular_spline_coord(a2, n_angular), n_angular, wa2, dwa2, ba2);
    float inv_dtheta = float(n_angular - 3) * 0.5f;

    float raw_value = 0.f;
    float raw_d_dr = 0.f;
    float raw_d_da1 = 0.f;
    float raw_d_da2 = 0.f;
    for(int ir = 0; ir < 4; ++ir) {
        int rr = br + ir;
        if(rr < 0 || rr >= n_radial) continue;
        for(int i1 = 0; i1 < 4; ++i1) {
            int aa1 = ba1 + i1;
            if(aa1 < 0 || aa1 >= n_angular) continue;
            for(int i2 = 0; i2 < 4; ++i2) {
                int aa2 = ba2 + i2;
                if(aa2 < 0 || aa2 >= n_angular) continue;
                float coeff = p[(rr * n_angular + aa1) * n_angular + aa2];
                raw_value += wr[ir] * wa1[i1] * wa2[i2] * coeff;
                raw_d_dr += dwr[ir] * wa1[i1] * wa2[i2] * coeff / knot_spacing;
                raw_d_da1 += wr[ir] * dwa1[i1] * wa2[i2] * coeff * inv_dtheta;
                raw_d_da2 += wr[ir] * wa1[i1] * dwa2[i2] * coeff * inv_dtheta;
            }
        }
    }

    float taper = 1.f;
    float d_taper_dr = 0.f;
    taper_width = std::max(taper_width, 1e-6f);
    float taper_start = cutoff - taper_width;
    if(r > taper_start) {
        float u = (cutoff - r) / taper_width;
        u = std::max(0.f, std::min(1.f, u));
        taper = u * u * (3.f - 2.f * u);
        d_taper_dr = -6.f * u * (1.f - u) / taper_width;
    }

    out.value = taper * raw_value;
    out.d_dr = taper * raw_d_dr + raw_value * d_taper_dr;
    out.d_da1 = taper * raw_d_da1;
    out.d_da2 = taper * raw_d_da2;
    return true;
}

static bool eval_cg_target_tensor(
        const float* p,
        int n_angular,
        int n_radial,
        const float dr[3],
        const float n_cg[3],
        float knot_spacing,
        float cutoff,
        float taper_width,
        TargetSplineEval& out) {

    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f) return false;
    float r = sqrtf(r2);
    if(r >= cutoff) return false;

    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};
    float a = n_cg[0] * unit[0] + n_cg[1] * unit[1] + n_cg[2] * unit[2];

    float wr[4], dwr[4], wa[4], dwa[4];
    int br = 0, ba = 0;
    radial_uniform_bspline_basis(r / knot_spacing, n_radial, wr, dwr, br);
    clamped_uniform_bspline_basis(angular_spline_coord(a, n_angular), n_angular, wa, dwa, ba);
    float inv_dtheta = float(n_angular - 3) * 0.5f;

    float raw_value = 0.f;
    float raw_d_dr = 0.f;
    float raw_d_da = 0.f;
    for(int ir = 0; ir < 4; ++ir) {
        int rr = br + ir;
        if(rr < 0 || rr >= n_radial) continue;
        for(int ia = 0; ia < 4; ++ia) {
            int aa = ba + ia;
            if(aa < 0 || aa >= n_angular) continue;
            float coeff = p[rr * n_angular + aa];
            raw_value += wr[ir] * wa[ia] * coeff;
            raw_d_dr += dwr[ir] * wa[ia] * coeff / knot_spacing;
            raw_d_da += wr[ir] * dwa[ia] * coeff * inv_dtheta;
        }
    }

    float taper = 1.f;
    float d_taper_dr = 0.f;
    taper_width = std::max(taper_width, 1e-6f);
    float taper_start = cutoff - taper_width;
    if(r > taper_start) {
        float u = (cutoff - r) / taper_width;
        u = std::max(0.f, std::min(1.f, u));
        taper = u * u * (3.f - 2.f * u);
        d_taper_dr = -6.f * u * (1.f - u) / taper_width;
    }

    out.value = taper * raw_value;
    out.d_dr = taper * raw_d_dr + raw_value * d_taper_dr;
    out.d_da = taper * raw_d_da;
    return true;
}

static inline void compute_cutoff_taper(
        float r,
        float cutoff,
        float taper_width,
        float& taper,
        float& d_taper_dr) {
    taper = 1.f;
    d_taper_dr = 0.f;
    taper_width = std::max(taper_width, 1e-6f);
    float taper_start = cutoff - taper_width;
    if(r > taper_start) {
        float u = (cutoff - r) / taper_width;
        u = std::max(0.f, std::min(1.f, u));
        taper = u * u * (3.f - 2.f * u);
        d_taper_dr = -6.f * u * (1.f - u) / taper_width;
    }
}

static inline float norm_dr(const float dr[3]) {
    return sqrtf(std::max(
                dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2],
                1e-12f));
}

static inline void apply_transformed_taper(
        QuadsplineEval& e,
        const float dr[3],
        float cutoff,
        float taper_width) {
    float taper = 1.f;
    float d_taper_dr = 0.f;
    compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
    float energy = e.value;
    e.value = taper * energy;
    e.d_dr = taper * e.d_dr + energy * d_taper_dr;
    e.d_da1 *= taper;
    e.d_da2 *= taper;
}

static inline void apply_transformed_taper(
        TargetSplineEval& e,
        const float dr[3],
        float cutoff,
        float taper_width) {
    float taper = 1.f;
    float d_taper_dr = 0.f;
    compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
    float energy = e.value;
    e.value = taper * energy;
    e.d_dr = taper * e.d_dr + energy * d_taper_dr;
    e.d_da *= taper;
}

static inline void apply_log1p_reduced_transform(
        QuadsplineEval& e,
        float reference_energy,
        float boltzmann_temperature,
        const float dr[3],
        float cutoff,
        float taper_width) {
    float deriv_scale = boltzmann_temperature * expf(e.value);
    e.value = reference_energy + boltzmann_temperature * expm1f(e.value);
    e.d_dr *= deriv_scale;
    e.d_da1 *= deriv_scale;
    e.d_da2 *= deriv_scale;
    apply_transformed_taper(e, dr, cutoff, taper_width);
}

static inline void apply_log1p_reduced_transform(
        TargetSplineEval& e,
        float reference_energy,
        float boltzmann_temperature,
        const float dr[3],
        float cutoff,
        float taper_width) {
    float deriv_scale = boltzmann_temperature * expf(e.value);
    e.value = reference_energy + boltzmann_temperature * expm1f(e.value);
    e.d_dr *= deriv_scale;
    e.d_da *= deriv_scale;
    apply_transformed_taper(e, dr, cutoff, taper_width);
}

static inline void apply_boltzmann_weight_transform(
        TargetSplineEval& e,
        float reference_energy,
        float boltzmann_temperature,
        float minimum_boltzmann_weight,
        const float dr[3],
        float cutoff,
        float taper_width) {
    float weight = std::max(e.value, minimum_boltzmann_weight);
    float deriv_scale = -boltzmann_temperature / weight;
    e.value = reference_energy - boltzmann_temperature * logf(weight);
    e.d_dr *= deriv_scale;
    e.d_da *= deriv_scale;
    apply_transformed_taper(e, dr, cutoff, taper_width);
}

static inline void load_vec6(VecArray data, int idx, float pos[3], float dir[3]) {
    pos[0] = data(0, idx);
    pos[1] = data(1, idx);
    pos[2] = data(2, idx);
    dir[0] = data(3, idx);
    dir[1] = data(4, idx);
    dir[2] = data(5, idx);
}

static inline void add_vec6_sens(
        VecArray sens, int idx, const float dpos[3], const float ddir[3]) {
    sens(0, idx) += dpos[0];
    sens(1, idx) += dpos[1];
    sens(2, idx) += dpos[2];
    sens(3, idx) += ddir[0];
    sens(4, idx) += ddir[1];
    sens(5, idx) += ddir[2];
}

static void accumulate_deriv(
        const float dr[3],
        const float n1[3],
        const float n2[3],
        const QuadsplineEval& e,
        float dpos1[3],
        float ddir1[3],
        float dpos2[3],
        float ddir2[3]) {

    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    float r = sqrtf(std::max(r2, 1e-12f));
    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};

    float d_dn[3] = {
        -e.d_da1 * n1[0] + e.d_da2 * n2[0],
        -e.d_da1 * n1[1] + e.d_da2 * n2[1],
        -e.d_da1 * n1[2] + e.d_da2 * n2[2],
    };
    float parallel = d_dn[0] * unit[0] + d_dn[1] * unit[1] + d_dn[2] * unit[2];
    float d_dr[3] = {
        e.d_dr * unit[0] + (d_dn[0] - parallel * unit[0]) * inv_r,
        e.d_dr * unit[1] + (d_dn[1] - parallel * unit[1]) * inv_r,
        e.d_dr * unit[2] + (d_dn[2] - parallel * unit[2]) * inv_r,
    };

    for(int k = 0; k < 3; ++k) {
        dpos1[k] -= d_dr[k];
        dpos2[k] += d_dr[k];
        ddir1[k] += -e.d_da1 * unit[k];
        ddir2[k] +=  e.d_da2 * unit[k];
    }
}

static inline void center_dr(
        VecArray first,
        int first_idx,
        VecArray second,
        int second_idx,
        float box_x,
        float box_y,
        float box_z,
        float dr[3]) {
    dr[0] = second(0, second_idx) - first(0, first_idx);
    dr[1] = second(1, second_idx) - first(1, first_idx);
    dr[2] = second(2, second_idx) - first(2, first_idx);
    if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
        simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);
}

static bool cached_positions_moved(
        VecArray pos,
        const vector<int>& index,
        const vector<float>& cached_pos,
        float max_disp2) {
    if(cached_pos.size() != index.size() * 3u)
        return true;
    for(size_t i = 0; i < index.size(); ++i) {
        int idx = index[i];
        float dx = pos(0, idx) - cached_pos[3 * i];
        float dy = pos(1, idx) - cached_pos[3 * i + 1];
        float dz = pos(2, idx) - cached_pos[3 * i + 2];
        if(dx * dx + dy * dy + dz * dz >= max_disp2)
            return true;
    }
    return false;
}

static void cache_positions(
        VecArray pos,
        const vector<int>& index,
        vector<float>& cached_pos) {
    cached_pos.resize(index.size() * 3u);
    for(size_t i = 0; i < index.size(); ++i) {
        int idx = index[i];
        cached_pos[3 * i] = pos(0, idx);
        cached_pos[3 * i + 1] = pos(1, idx);
        cached_pos[3 * i + 2] = pos(2, idx);
    }
}

static bool cached_cgl_bodies_moved(
        VecArray cg,
        const vector<int>& index,
        const vector<float>& cached_body,
        float axis_radius,
        float max_disp2) {
    if(cached_body.size() != index.size() * 6u)
        return true;
    float max_disp = sqrtf(std::max(max_disp2, 0.f));
    for(size_t i = 0; i < index.size(); ++i) {
        int idx = index[i];
        float dx = cg(0, idx) - cached_body[6 * i];
        float dy = cg(1, idx) - cached_body[6 * i + 1];
        float dz = cg(2, idx) - cached_body[6 * i + 2];
        float dn0 = cg(3, idx) - cached_body[6 * i + 3];
        float dn1 = cg(4, idx) - cached_body[6 * i + 4];
        float dn2 = cg(5, idx) - cached_body[6 * i + 5];
        float body_disp = sqrtf(dx * dx + dy * dy + dz * dz)
                        + axis_radius * sqrtf(dn0 * dn0 + dn1 * dn1 + dn2 * dn2);
        if(body_disp >= max_disp)
            return true;
    }
    return false;
}

static void cache_cgl_bodies(
        VecArray cg,
        const vector<int>& index,
        vector<float>& cached_body) {
    cached_body.resize(index.size() * 6u);
    for(size_t i = 0; i < index.size(); ++i) {
        int idx = index[i];
        for(int d = 0; d < 6; ++d)
            cached_body[6 * i + d] = cg(d, idx);
    }
}

static inline float cgl_pairlist_disp_threshold2(float cache_buffer) {
    float half_buffer = 0.5f * cache_buffer;
    return half_buffer * half_buffer;
}

static inline float dot3(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline float clamp_range(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
}

static CGLBodySupport read_cgl_body_support(hid_t grp) {
    float length_conv = read_attribute<float>(grp, ".", "length_conversion_ang_per_nm", 10.f);
    float bead_cutoff_nm = read_attribute<float>(grp, ".", "bead_nonbonded_cutoff_nm", 0.f);
    CGLBodySupport out;
    out.axis_radius = read_attribute<float>(grp, ".", "max_axis_radius_ang", 0.f);
    out.perp_radius = read_attribute<float>(grp, ".", "max_perp_radius_ang", 0.f);
    out.bead_cutoff = bead_cutoff_nm * length_conv;
    out.has_geometry = out.axis_radius > 0.f && out.perp_radius > 0.f && out.bead_cutoff > 0.f;
    return out;
}

static float cgl_point_distance2(const float dr[3], const float n[3], float axis_radius) {
    float s = clamp_range(dot3(n, dr), -axis_radius, axis_radius);
    float d0 = dr[0] - s * n[0];
    float d1 = dr[1] - s * n[1];
    float d2 = dr[2] - s * n[2];
    return d0 * d0 + d1 * d1 + d2 * d2;
}

static float cgl_cgl_distance2(
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float axis_radius) {
    float b = dot3(n1, n2);
    float d1 = dot3(n1, dr);
    float d2 = dot3(n2, dr);
    float best = std::numeric_limits<float>::max();
    auto add_candidate = [&](float s, float t) {
        float d0 = dr[0] + t * n2[0] - s * n1[0];
        float d1v = dr[1] + t * n2[1] - s * n1[1];
        float d2v = dr[2] + t * n2[2] - s * n1[2];
        best = std::min(best, d0 * d0 + d1v * d1v + d2v * d2v);
    };

    float denom = 1.f - b * b;
    if(fabsf(denom) > 1e-6f) {
        float s = (d1 - b * d2) / denom;
        float t = (b * d1 - d2) / denom;
        if(s >= -axis_radius && s <= axis_radius &&
           t >= -axis_radius && t <= axis_radius) {
            add_candidate(s, t);
        }
    }

    for(float s : {-axis_radius, axis_radius}) {
        float t = clamp_range(b * s - d2, -axis_radius, axis_radius);
        add_candidate(s, t);
    }
    for(float t : {-axis_radius, axis_radius}) {
        float s = clamp_range(d1 + b * t, -axis_radius, axis_radius);
        add_candidate(s, t);
    }
    return best;
}

static bool cgl_cgl_pairlist_candidate(
        VecArray cg,
        int idx1,
        int idx2,
        const float dr[3],
        const CGLBodySupport& support,
        float cache_buffer,
        float fallback_cutoff2) {
    if(!support.has_geometry) {
        float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        return r2 < fallback_cutoff2;
    }
    float x1[3], n1[3], x2[3], n2[3];
    load_vec6(cg, idx1, x1, n1);
    load_vec6(cg, idx2, x2, n2);
    float contact = support.bead_cutoff + cache_buffer + 2.f * support.perp_radius;
    return cgl_cgl_distance2(dr, n1, n2, support.axis_radius) < contact * contact;
}

static bool cgl_target_pairlist_candidate(
        VecArray cg,
        int idx_cg,
        const float dr[3],
        const CGLBodySupport& support,
        float cache_buffer,
        float fallback_cutoff2) {
    if(!support.has_geometry) {
        float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        return r2 < fallback_cutoff2;
    }
    float x[3], n[3];
    load_vec6(cg, idx_cg, x, n);
    float contact = support.bead_cutoff + cache_buffer + support.perp_radius;
    return cgl_point_distance2(dr, n, support.axis_radius) < contact * contact;
}

}  // namespace

struct CGLOrientationState : public CoordNode {
    float rotational_inertia;
    float thermostat_timescale;

    CGLOrientationState(hid_t grp)
        : CoordNode(int(h5::get_dset_size(2, grp, "direction")[0]), 3)
        , rotational_inertia(read_attribute<float>(grp, ".", "rotational_inertia_up", 1.f))
        , thermostat_timescale(read_attribute<float>(grp, ".", "rotational_thermostat_timescale", 5.f))
    {
        if(rotational_inertia <= 0.f)
            throw string("cgl_orientation_state rotational_inertia_up must be positive");
        if(thermostat_timescale <= 0.f)
            throw string("cgl_orientation_state rotational_thermostat_timescale must be positive");

        traverse_dset<2, float>(grp, "direction", [&](size_t i, size_t d, float v) {
            output(int(d), int(i)) = v;
        });
        for(int i = 0; i < n_elem; ++i) {
            float3 n = load_vec<3>(output, i);
            float norm = mag(n);
            if(norm <= 1e-6f)
                throw string("cgl_orientation_state direction contains a zero vector");
            store_vec(output, i, n / norm);
        }
    }

    virtual void compute_value(ComputeMode) override {}
    virtual void propagate_deriv() override {}
};

namespace {

struct CGLDynamicOrientationRuntime {
    CGLOrientationState* state = nullptr;
    VecArrayStorage angular_momentum;
    float temperature = 1.f;
    float thermostat_delta_t = 0.f;
    uint32_t random_seed = 0u;
    uint64_t thermostat_invocations = 0u;

    CGLDynamicOrientationRuntime(CGLOrientationState* state_, uint32_t random_seed_)
        : state(state_)
        , angular_momentum(3, state_ ? state_->n_elem : 1)
        , random_seed(random_seed_)
    {
        fill(angular_momentum, 0.f);
    }
};

static std::map<DerivEngine*, std::unique_ptr<CGLDynamicOrientationRuntime>>& orientation_runtime_map() {
    static std::map<DerivEngine*, std::unique_ptr<CGLDynamicOrientationRuntime>> m;
    return m;
}

static CGLDynamicOrientationRuntime* get_orientation_runtime(DerivEngine* engine) {
    auto& m = orientation_runtime_map();
    auto it = m.find(engine);
    return it == m.end() ? nullptr : it->second.get();
}

struct CGLGLERuntime {
    vector<int> atom_index;
    vector<float> memory_tau;
    vector<float> coupling;
    VecArrayStorage aux_momentum;
    float temperature = 1.f;
    float delta_t = 0.f;
    uint32_t random_seed = 0u;
    uint64_t thermostat_invocations = 0u;

    CGLGLERuntime(
            const vector<int>& atom_index_,
            uint32_t random_seed_,
            const vector<float>& memory_tau_,
            const vector<float>& coupling_)
        : atom_index(atom_index_)
        , memory_tau(memory_tau_)
        , coupling(coupling_)
        , aux_momentum(3, int(std::max<size_t>(1, atom_index_.size() * memory_tau_.size())))
        , random_seed(random_seed_)
    {
        fill(aux_momentum, 0.f);
    }

    int n_mode() const { return int(memory_tau.size()); }
    int aux_index(int mode, int lipid) const {
        return mode * int(atom_index.size()) + lipid;
    }
};

static std::map<DerivEngine*, std::unique_ptr<CGLGLERuntime>>& cgl_gle_runtime_map() {
    static std::map<DerivEngine*, std::unique_ptr<CGLGLERuntime>> m;
    return m;
}

static CGLGLERuntime* get_cgl_gle_runtime(DerivEngine* engine) {
    auto& m = cgl_gle_runtime_map();
    auto it = m.find(engine);
    return it == m.end() ? nullptr : it->second.get();
}

static float3 project_tangent(float3 v, const float3& n) {
    return v - n * dot(v, n);
}

static float3 rotate_unit_vector(const float3& n, const float3& omega, float dt) {
    float omega_mag = mag(omega);
    if(omega_mag <= 1e-8f) return n;
    float angle = omega_mag * dt;
    float3 axis = omega / omega_mag;
    float c = cosf(angle);
    float s = sinf(angle);
    float3 out = n * c + cross(axis, n) * s + axis * (dot(axis, n) * (1.f - c));
    float out_norm = mag(out);
    return out_norm > 1e-8f ? out / out_norm : n;
}

}  // namespace

struct ComposeVector6D : public CoordNode {
    CoordNode& pos;
    CoordNode* orientation_state;
    vector<int> elem_index;
    vector<int> orientation_index;
    vector<float> direction;
    bool dynamic_orientation;
    bool state_orientation;
    float box_x;
    float box_y;
    float box_z;

    ComposeVector6D(hid_t grp, CoordNode& pos_node);
    ComposeVector6D(hid_t grp, const ArgList& arguments);

    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
};

ComposeVector6D::ComposeVector6D(hid_t grp, CoordNode& pos_node)
    : CoordNode(h5::get_dset_size(1, grp, "elem_index")[0], 6)
    , pos(pos_node)
    , orientation_state(nullptr)
    , dynamic_orientation(h5_exists(grp, "orientation_index"))
    , state_orientation(false)
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
{
    check_elem_width(pos_node, 3);

    int n = n_elem;
    elem_index.resize(n);
    traverse_dset<1, int32_t>(grp, "elem_index", [&](size_t i, int32_t v) {
        elem_index[i] = v;
    });

    if(dynamic_orientation) {
        vector<hsize_t> sz = get_dset_size(1, grp, "orientation_index");
        if(int(sz[0]) != n_elem)
            throw string("compose_vector6d orientation_index size must match elem_index");
        orientation_index.resize(n);
        traverse_dset<1, int32_t>(grp, "orientation_index", [&](size_t i, int32_t v) {
            orientation_index[i] = v;
        });
    }

    direction.resize(n * 3);
    traverse_dset<2, float>(grp, "direction", [&](size_t i, size_t d, float v) {
        direction[i * 3 + d] = v;
    });

    fill(output, 0.f);
}

ComposeVector6D::ComposeVector6D(hid_t grp, const ArgList& arguments)
    : ComposeVector6D(grp, *arguments.at(0))
{
    if(arguments.size() == 2u) {
        orientation_state = arguments.at(1);
        check_elem_width(*orientation_state, 3);
        if(orientation_state->n_elem != n_elem)
            throw string("compose_vector6d orientation state size must match elem_index");
        state_orientation = true;
        dynamic_orientation = false;
    } else if(arguments.size() != 1u) {
        throw string("compose_vector6d expects pos or pos plus cgl_orientation_state");
    }
}

void ComposeVector6D::compute_value(ComputeMode) {
    VecArray posc = pos.output;
    VecArray orient = orientation_state ? orientation_state->output : VecArray();
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        output(0, i) = posc(0, ai);
        output(1, i) = posc(1, ai);
        output(2, i) = posc(2, ai);
        if(state_orientation) {
            float3 n = load_vec<3>(orient, i);
            float norm = mag(n);
            if(norm > 1e-6f) n /= norm;
            else n = make_vec3(direction[i * 3], direction[i * 3 + 1], direction[i * 3 + 2]);
            output(3, i) = n.x();
            output(4, i) = n.y();
            output(5, i) = n.z();
        } else if(dynamic_orientation) {
            int oi = orientation_index[i];
            float dr[3] = {
                posc(0, oi) - posc(0, ai),
                posc(1, oi) - posc(1, ai),
                posc(2, oi) - posc(2, ai),
            };
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
                simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);
            float norm = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
            if(norm > 1e-6f) {
                float inv_norm = 1.f / norm;
                output(3, i) = dr[0] * inv_norm;
                output(4, i) = dr[1] * inv_norm;
                output(5, i) = dr[2] * inv_norm;
            } else {
                output(3, i) = direction[i * 3];
                output(4, i) = direction[i * 3 + 1];
                output(5, i) = direction[i * 3 + 2];
            }
        } else {
            output(3, i) = direction[i * 3];
            output(4, i) = direction[i * 3 + 1];
            output(5, i) = direction[i * 3 + 2];
        }
    }
}

void ComposeVector6D::propagate_deriv() {
    VecArray pos_sens = pos.sens;
    VecArray posc = pos.output;
    VecArray orient_sens = orientation_state ? orientation_state->sens : VecArray();
    VecArray orient = orientation_state ? orientation_state->output : VecArray();
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        pos_sens(0, ai) += sens(0, i);
        pos_sens(1, ai) += sens(1, i);
        pos_sens(2, ai) += sens(2, i);

        if(state_orientation) {
            float3 n = load_vec<3>(orient, i);
            float norm2 = dot(n, n);
            if(norm2 <= 1e-12f) continue;
            float norm = sqrtf(norm2);
            n /= norm;
            float3 g = make_vec3(sens(3, i), sens(4, i), sens(5, i));
            g = (g - n * dot(g, n)) / norm;
            orient_sens(0, i) += g.x();
            orient_sens(1, i) += g.y();
            orient_sens(2, i) += g.z();
            continue;
        }

        if(!dynamic_orientation) continue;

        int oi = orientation_index[i];
        float dr[3] = {
            posc(0, oi) - posc(0, ai),
            posc(1, oi) - posc(1, ai),
            posc(2, oi) - posc(2, ai),
        };
        if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
            simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);
        float norm2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        if(norm2 <= 1e-12f) continue;
        float norm = sqrtf(norm2);
        float inv_norm = 1.f / norm;
        float n[3] = {dr[0] * inv_norm, dr[1] * inv_norm, dr[2] * inv_norm};
        float g[3] = {sens(3, i), sens(4, i), sens(5, i)};
        float g_dot_n = g[0] * n[0] + g[1] * n[1] + g[2] * n[2];
        for(int k = 0; k < 3; ++k) {
            float d_dk = (g[k] - g_dot_n * n[k]) * inv_norm;
            pos_sens(k, ai) -= d_dk;
            pos_sens(k, oi) += d_dk;
        }
    }
}

struct CGLipidPairPotential : public PotentialNode {
    CoordNode& cg_pos;
    vector<float> interaction_param;
    vector<int> index;
    vector<int> type;
    vector<int> id;
    vector<float> reference_energy_eup;
    int n_type1;
    int n_type2;
    int n_param;
    int n_modes;
    int n_radial;
    int n_angular;
    bool full_tensor;
    float box_x;
    float box_y;
    float box_z;
    float knot_spacing;
    float cutoff;
    float taper_width;
    bool log1p_reduced_transform;
    float boltzmann_temperature;
    float cache_buffer;
    CGLBodySupport body_support;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_body;
    vector<CGLPairIndex> active_pairs;

    CGLipidPairPotential(hid_t grp, CoordNode& cg_pos);
    void ensure_pairlist(VecArray cg);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == interaction_param.size()) interaction_param = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        pairlist_valid = false;
    }
};

CGLipidPairPotential::CGLipidPairPotential(hid_t grp, CoordNode& cg_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , n_type1(0)
    , n_type2(0)
    , n_param(0)
    , n_modes(read_attribute<int>(grp, ".", "n_modes", 0))
    , n_radial(read_attribute<int>(grp, ".", "n_radial", CG_LIPID_N_RADIAL))
    , n_angular(read_attribute<int>(grp, ".", "n_angular", CG_LIPID_N_ANGULAR))
    , full_tensor(false)
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", CG_LIPID_DEFAULT_KNOT_SPACING))
    , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                float(n_radial - 2) * CG_LIPID_DEFAULT_KNOT_SPACING))
    , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
    , log1p_reduced_transform(read_attribute<int>(grp, ".", "log1p_reduced_transform", 0) != 0)
    , boltzmann_temperature(read_attribute<float>(grp, ".", "boltzmann_temperature_upside", 0.f))
    , cache_buffer(std::max(0.f, read_attribute<float>(grp, ".", "pairlist_buffer_ang", 1.f)))
    , body_support(read_cgl_body_support(grp))
    , cached_box_x(0.f)
    , cached_box_y(0.f)
    , cached_box_z(0.f)
    , pairlist_valid(false)
{
    check_elem_width(cg_pos, 6);
    H5Obj pi_obj = open_group(grp, "pair_interaction");
    hid_t pi = pi_obj.get();
    interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
    int expected_n_param = n_radial + n_modes * (2 * n_angular + n_radial);
    int expected_tensor_param = n_radial * n_angular * n_angular;
    full_tensor = (n_modes == 0 && n_radial > 3 && n_angular > 3 && n_param == expected_tensor_param);
    if(!full_tensor && (n_modes <= 0 || n_radial <= 3 || n_angular <= 3 || n_param != expected_n_param))
        throw string("cg_lipid_pair requires full multimode params with matching n_modes/n_radial/n_angular attrs");
    if(H5Lexists(pi, "reference_energy_eup", H5P_DEFAULT) > 0) {
        reference_energy_eup = read_float_dataset(pi, "reference_energy_eup");
        if(int(reference_energy_eup.size()) != n_type1 * n_type2)
            throw string("cg_lipid_pair reference_energy_eup shape mismatch");
    } else {
        reference_energy_eup.assign(n_type1 * n_type2, 0.f);
    }
    if(log1p_reduced_transform && boltzmann_temperature <= 0.f)
        throw string("cg_lipid_pair log1p-reduced transform requires positive temperature");
    index = read_int_dataset(pi, "index");
    type = read_int_dataset(pi, "type");
    id = read_int_dataset(pi, "id");
    if(index.size() != type.size() || index.size() != id.size())
        throw string("cg_lipid_pair index/type/id size mismatch");
}

void CGLipidPairPotential::ensure_pairlist(VecArray cg) {
    bool rebuild = !pairlist_valid || cache_buffer <= 0.f ||
        cached_box_x != box_x || cached_box_y != box_y || cached_box_z != box_z ||
        cached_cgl_bodies_moved(cg, index, cached_body, body_support.axis_radius,
                cgl_pairlist_disp_threshold2(cache_buffer));
    if(!rebuild)
        return;

    active_pairs.clear();
    float pairlist_cutoff = cutoff + cache_buffer;
    float pairlist_cutoff2 = pairlist_cutoff * pairlist_cutoff;
    for(size_t ai = 0; ai < index.size(); ++ai) {
        int t1 = type[ai];
        if(t1 < 0 || t1 >= n_type1) continue;
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if((id[ai] >> 4) == (id[bi] >> 4)) continue;
            int t2 = type[bi];
            if(t2 < 0 || t2 >= n_type2) continue;
            float dr[3];
            center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
            if(cgl_cgl_pairlist_candidate(
                    cg, index[ai], index[bi], dr, body_support, cache_buffer, pairlist_cutoff2))
                active_pairs.push_back({ai, bi});
        }
    }
    cache_cgl_bodies(cg, index, cached_body);
    cached_box_x = box_x;
    cached_box_y = box_y;
    cached_box_z = box_z;
    pairlist_valid = true;
}

void CGLipidPairPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray cg = cg_pos.output;
    ensure_pairlist(cg);
    float total = 0.f;
    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type[ai];
        int t2 = type[bi];

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);

        QuadsplineEval e;
        bool ok = full_tensor
            ? eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e)
            : eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_modes, n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e);
        if(ok) {
            if(log1p_reduced_transform)
                apply_log1p_reduced_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, dr, cutoff, taper_width);
            total += e.value;
        }
    }
    potential = total;
}

void CGLipidPairPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    ensure_pairlist(cg);

    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type[ai];
        int t2 = type[bi];

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);

        QuadsplineEval e;
        bool ok = full_tensor
            ? eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e)
            : eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_modes, n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e);
        if(!ok)
            continue;
        if(log1p_reduced_transform)
            apply_log1p_reduced_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, dr, cutoff, taper_width);

        float dpos1[3] = {0.f, 0.f, 0.f};
        float ddir1[3] = {0.f, 0.f, 0.f};
        float dpos2[3] = {0.f, 0.f, 0.f};
        float ddir2[3] = {0.f, 0.f, 0.f};
        accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
        add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
        add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
    }
}

struct CGLipidDensityPotential : public PotentialNode {
    CoordNode& cg_pos;
    vector<int> index;
    vector<float> kernel_coeff;
    vector<float> embedding_coeff;
    int n_kernel;
    int n_embedding;
    float box_x;
    float box_y;
    float box_z;
    float kernel_spacing;
    float kernel_cutoff;
    float rho_min;
    float rho_spacing;

    CGLipidDensityPotential(hid_t grp, CoordNode& cg_pos);
    void accumulate_density(VecArray cg, vector<float>& rho) const;
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == embedding_coeff.size()) embedding_coeff = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
    }
};

CGLipidDensityPotential::CGLipidDensityPotential(hid_t grp, CoordNode& cg_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , n_kernel(read_attribute<int>(grp, ".", "kernel_n_knot", 0))
    , n_embedding(read_attribute<int>(grp, ".", "embedding_n_knot", 0))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , kernel_spacing(read_attribute<float>(grp, ".", "kernel_knot_spacing_ang", 0.f))
    , kernel_cutoff(read_attribute<float>(grp, ".", "kernel_cutoff_ang", 0.f))
    , rho_min(read_attribute<float>(grp, ".", "embedding_rho_min", 0.f))
    , rho_spacing(read_attribute<float>(grp, ".", "embedding_rho_spacing", 0.f))
{
    check_elem_width(cg_pos, 6);
    index = read_int_dataset(grp, "index");
    kernel_coeff = read_float_dataset_1d(grp, "kernel_coeff");
    embedding_coeff = read_float_dataset_1d(grp, "embedding_coeff");
    if(index.empty()) throw string("cg_lipid_density requires at least one CGL index");
    if(n_kernel <= 3) n_kernel = int(kernel_coeff.size());
    if(n_embedding <= 3) n_embedding = int(embedding_coeff.size());
    if(int(kernel_coeff.size()) != n_kernel)
        throw string("cg_lipid_density kernel_coeff size does not match kernel_n_knot");
    if(int(embedding_coeff.size()) != n_embedding)
        throw string("cg_lipid_density embedding_coeff size does not match embedding_n_knot");
    if(kernel_spacing <= 0.f || kernel_cutoff <= 0.f || rho_spacing <= 0.f)
        throw string("cg_lipid_density requires positive kernel/embedding spacing and cutoff");
}

void CGLipidDensityPotential::accumulate_density(VecArray cg, vector<float>& rho) const {
    rho.assign(index.size(), 0.f);
    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            float dr[3];
            center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
            float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if(r2 <= 1e-12f) continue;
            float r = sqrtf(r2);
            if(r >= kernel_cutoff) continue;
            Vec<2> w = radial_deBoor_value_and_deriv(
                    kernel_coeff.data(), r / kernel_spacing, n_kernel);
            rho[ai] += w.x();
            rho[bi] += w.x();
        }
    }
}

void CGLipidDensityPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray cg = cg_pos.output;
    vector<float> rho;
    accumulate_density(cg, rho);
    float total = 0.f;
    for(float r: rho) {
        float coord = 1.f + (r - rho_min) / rho_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        total += emb.x();
    }
    potential = total;
}

void CGLipidDensityPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    vector<float> rho;
    accumulate_density(cg, rho);
    vector<float> dF_drho(index.size(), 0.f);
    for(size_t i = 0; i < index.size(); ++i) {
        float coord = 1.f + (rho[i] - rho_min) / rho_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        dF_drho[i] = emb.y() / rho_spacing;
    }

    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            float dr[3];
            center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
            float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if(r2 <= 1e-12f) continue;
            float r = sqrtf(r2);
            if(r >= kernel_cutoff) continue;
            Vec<2> w = radial_deBoor_value_and_deriv(
                    kernel_coeff.data(), r / kernel_spacing, n_kernel);
            float dE_dr = (dF_drho[ai] + dF_drho[bi]) * w.y() / kernel_spacing;
            float inv_r = 1.f / r;
            float dpos1[3] = {
                -dE_dr * dr[0] * inv_r,
                -dE_dr * dr[1] * inv_r,
                -dE_dr * dr[2] * inv_r,
            };
            float dpos2[3] = {
                 dE_dr * dr[0] * inv_r,
                 dE_dr * dr[1] * inv_r,
                 dE_dr * dr[2] * inv_r,
            };
            float zero[3] = {0.f, 0.f, 0.f};
            add_vec6_sens(cg_sens, index[ai], dpos1, zero);
            add_vec6_sens(cg_sens, index[bi], dpos2, zero);
        }
    }
}

struct CGLipidSCPotential : public PotentialNode {
    CoordNode& sc_pos;
    CoordNode& cg_pos;
    vector<float> interaction_param;
    vector<int> index1;
    vector<int> type1;
    vector<int> id1;
    vector<int> index2;
    vector<int> type2;
    vector<int> id2;
    vector<float> reference_energy_eup;
    int n_type1;
    int n_type2;
    int n_param;
    int n_radial;
    int n_angular;
    float box_x;
    float box_y;
    float box_z;
    float knot_spacing;
    float cutoff;
    float taper_width;
    bool log1p_reduced_transform;
    float boltzmann_temperature;
    float cache_buffer;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_pos1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;

    CGLipidSCPotential(hid_t grp, CoordNode& sc_pos, CoordNode& cg_pos);
    void ensure_pairlist(VecArray sc, VecArray cg);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == interaction_param.size()) interaction_param = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        pairlist_valid = false;
    }
};

CGLipidSCPotential::CGLipidSCPotential(hid_t grp, CoordNode& sc_pos_, CoordNode& cg_pos_)
    : PotentialNode()
    , sc_pos(sc_pos_)
    , cg_pos(cg_pos_)
    , n_type1(0)
    , n_type2(0)
    , n_param(0)
    , n_radial(read_attribute<int>(grp, ".", "n_radial", CG_LIPID_N_RADIAL))
    , n_angular(read_attribute<int>(grp, ".", "n_angular", CG_LIPID_N_ANGULAR))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", CG_LIPID_DEFAULT_KNOT_SPACING))
    , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                (CG_LIPID_N_RADIAL - 2) * CG_LIPID_DEFAULT_KNOT_SPACING))
    , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
    , log1p_reduced_transform(read_attribute<int>(grp, ".", "log1p_reduced_transform", 0) != 0)
    , boltzmann_temperature(read_attribute<float>(grp, ".", "boltzmann_temperature_upside", 0.f))
    , cache_buffer(std::max(0.f, read_attribute<float>(grp, ".", "pairlist_buffer_ang", 1.f)))
    , cached_box_x(0.f)
    , cached_box_y(0.f)
    , cached_box_z(0.f)
    , pairlist_valid(false)
{
    check_elem_width(sc_pos, 6);
    check_elem_width(cg_pos, 6);
    H5Obj pi_obj = open_group(grp, "pair_interaction");
    hid_t pi = pi_obj.get();
    interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
    int expected_tensor_param = n_radial * n_angular * n_angular;
    if(n_radial <= 3 || n_angular <= 3 || n_param != expected_tensor_param)
        throw string("cg_lipid_sc requires full tensor params matching n_radial*n_angular*n_angular");
    if(H5Lexists(pi, "reference_energy_eup", H5P_DEFAULT) > 0) {
        reference_energy_eup = read_float_dataset(pi, "reference_energy_eup");
        if(int(reference_energy_eup.size()) != n_type1 * n_type2)
            throw string("cg_lipid_sc reference_energy_eup shape mismatch");
    } else {
        reference_energy_eup.assign(n_type1 * n_type2, 0.f);
    }
    if(log1p_reduced_transform && boltzmann_temperature <= 0.f)
        throw string("cg_lipid_sc log1p-reduced transform requires positive temperature");
    index1 = read_int_dataset(pi, "index1");
    type1 = read_int_dataset(pi, "type1");
    id1 = read_int_dataset(pi, "id1");
    index2 = read_int_dataset(pi, "index2");
    type2 = read_int_dataset(pi, "type2");
    id2 = read_int_dataset(pi, "id2");
    if(index1.size() != type1.size() || index1.size() != id1.size())
        throw string("cg_lipid_sc source1 index/type/id size mismatch");
    if(index2.size() != type2.size() || index2.size() != id2.size())
        throw string("cg_lipid_sc source2 index/type/id size mismatch");
}

void CGLipidSCPotential::ensure_pairlist(VecArray sc, VecArray cg) {
    bool rebuild = !pairlist_valid || cache_buffer <= 0.f ||
        cached_box_x != box_x || cached_box_y != box_y || cached_box_z != box_z ||
        cached_positions_moved(sc, index1, cached_pos1, cgl_pairlist_disp_threshold2(cache_buffer)) ||
        cached_positions_moved(cg, index2, cached_pos2, cgl_pairlist_disp_threshold2(cache_buffer));
    if(!rebuild)
        return;

    active_pairs.clear();
    float pairlist_cutoff = cutoff + cache_buffer;
    float pairlist_cutoff2 = pairlist_cutoff * pairlist_cutoff;
    for(size_t ai = 0; ai < index1.size(); ++ai) {
        int t1 = type1[ai];
        if(t1 < 0 || t1 >= n_type1) continue;
        for(size_t bi = 0; bi < index2.size(); ++bi) {
            if((id1[ai] >> 4) == (id2[bi] >> 4)) continue;
            int t2 = type2[bi];
            if(t2 < 0 || t2 >= n_type2) continue;
            float dr[3];
            center_dr(sc, index1[ai], cg, index2[bi], box_x, box_y, box_z, dr);
            float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if(r2 < pairlist_cutoff2)
                active_pairs.push_back({ai, bi});
        }
    }
    cache_positions(sc, index1, cached_pos1);
    cache_positions(cg, index2, cached_pos2);
    cached_box_x = box_x;
    cached_box_y = box_y;
    cached_box_z = box_z;
    pairlist_valid = true;
}

void CGLipidSCPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray sc = sc_pos.output;
    VecArray cg = cg_pos.output;
    ensure_pairlist(sc, cg);
    float total = 0.f;

    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type1[ai];
        int t2 = type2[bi];

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(sc, index1[ai], x1, n1);
        load_vec6(cg, index2[bi], x2, n2);
        center_dr(sc, index1[ai], cg, index2[bi], box_x, box_y, box_z, dr);

        QuadsplineEval e;
        bool ok = eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                n_angular, n_radial,
                dr, n1, n2, knot_spacing, cutoff,
                log1p_reduced_transform ? 0.f : taper_width, e);
        if(ok) {
            if(log1p_reduced_transform)
                apply_log1p_reduced_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, dr, cutoff, taper_width);
            total += e.value;
        }
    }
    potential = total;
}

void CGLipidSCPotential::propagate_deriv() {
    VecArray sc = sc_pos.output;
    VecArray sc_sens = sc_pos.sens;
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    ensure_pairlist(sc, cg);

    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type1[ai];
        int t2 = type2[bi];

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(sc, index1[ai], x1, n1);
        load_vec6(cg, index2[bi], x2, n2);
        center_dr(sc, index1[ai], cg, index2[bi], box_x, box_y, box_z, dr);

        QuadsplineEval e;
        bool ok = eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                n_angular, n_radial,
                dr, n1, n2, knot_spacing, cutoff,
                log1p_reduced_transform ? 0.f : taper_width, e);
        if(!ok)
            continue;
        if(log1p_reduced_transform)
            apply_log1p_reduced_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, dr, cutoff, taper_width);

        float dpos1[3] = {0.f, 0.f, 0.f};
        float ddir1[3] = {0.f, 0.f, 0.f};
        float dpos2[3] = {0.f, 0.f, 0.f};
        float ddir2[3] = {0.f, 0.f, 0.f};
        accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
        add_vec6_sens(sc_sens, index1[ai], dpos1, ddir1);
        add_vec6_sens(cg_sens, index2[bi], dpos2, ddir2);
    }
}

struct CGLipidTargetPotential : public PotentialNode {
    CoordNode& cg_pos;
    CoordNode& tgt_pos;
    vector<float> interaction_param;
    vector<int> index1;
    vector<int> type1;
    vector<int> id1;
    vector<int> index2;
    vector<int> type2;
    vector<int> id2;
    vector<float> reference_energy_eup;
    int n_type1;
    int n_type2;
    int n_param;
    int n_modes;
    int n_radial;
    int n_angular;
    float box_x;
    float box_y;
    float box_z;
    float knot_spacing;
    float cutoff;
    float taper_width;
    bool boltzmann_weight_transform;
    bool log1p_reduced_transform;
    float boltzmann_temperature;
    float minimum_boltzmann_weight;
    float cache_buffer;
    CGLBodySupport body_support;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_body1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;

    CGLipidTargetPotential(hid_t grp, CoordNode& cg_pos_, CoordNode& tgt_pos_);
    void ensure_pairlist(VecArray cg, VecArray tgt);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == interaction_param.size()) interaction_param = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        pairlist_valid = false;
    }
};

struct CGLipidSCOneBody : public CoordNode {
    CoordNode& sc_pos;
    CoordNode& cg_pos;
    vector<float> interaction_param;
    vector<int> row_type;
    vector<int> row_residue_index;
    vector<int> row_rotamer_index;
    vector<int> row_group_count;
    vector<int> row_index;
    vector<int> index2;
    vector<int> type2;
    vector<int> id2;
    vector<float> reference_energy_eup;
    int n_type1;
    int n_type2;
    int n_param;
    int n_radial;
    int n_angular;
    float box_x;
    float box_y;
    float box_z;
    float knot_spacing;
    float cutoff;
    float taper_width;
    bool log1p_reduced_transform;
    float boltzmann_temperature;
    float cache_buffer;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_pos1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;

    CGLipidSCOneBody(hid_t grp, CoordNode& sc_pos_, CoordNode& cg_pos_)
        : CoordNode(sc_pos_.n_elem, 1)
        , sc_pos(sc_pos_)
        , cg_pos(cg_pos_)
        , n_type1(0)
        , n_type2(0)
        , n_param(0)
        , n_radial(read_attribute<int>(grp, ".", "n_radial", CG_LIPID_N_RADIAL))
        , n_angular(read_attribute<int>(grp, ".", "n_angular", CG_LIPID_N_ANGULAR))
        , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
        , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
        , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
        , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", CG_LIPID_DEFAULT_KNOT_SPACING))
        , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                    (CG_LIPID_N_RADIAL - 2) * CG_LIPID_DEFAULT_KNOT_SPACING))
        , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
        , log1p_reduced_transform(read_attribute<int>(grp, ".", "log1p_reduced_transform", 0) != 0)
        , boltzmann_temperature(read_attribute<float>(grp, ".", "boltzmann_temperature_upside", 0.f))
        , cache_buffer(std::max(0.f, read_attribute<float>(grp, ".", "pairlist_buffer_ang", 1.f)))
        , cached_box_x(0.f)
        , cached_box_y(0.f)
        , cached_box_z(0.f)
        , pairlist_valid(false)
    {
        check_elem_width(sc_pos, 6);
        check_elem_width(cg_pos, 6);
        H5Obj pi_obj = open_group(grp, "pair_interaction");
        hid_t pi = pi_obj.get();
        interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
        int expected_tensor_param = n_radial * n_angular * n_angular;
        if(n_radial <= 3 || n_angular <= 3 || n_param != expected_tensor_param)
            throw string("cg_lipid_rotamer_sc requires full tensor params matching n_radial*n_angular*n_angular");
        if(H5Lexists(pi, "reference_energy_eup", H5P_DEFAULT) > 0) {
            reference_energy_eup = read_float_dataset(pi, "reference_energy_eup");
            if(int(reference_energy_eup.size()) != n_type1 * n_type2)
                throw string("cg_lipid_rotamer_sc reference_energy_eup shape mismatch");
        } else {
            reference_energy_eup.assign(n_type1 * n_type2, 0.f);
        }
        if(log1p_reduced_transform && boltzmann_temperature <= 0.f)
            throw string("cg_lipid_rotamer_sc log1p-reduced transform requires positive temperature");

        row_type = read_int_dataset(pi, "type1");
        row_residue_index = read_int_dataset(pi, "row_residue_index");
        row_rotamer_index = read_int_dataset(pi, "row_rotamer_index");
        index2 = read_int_dataset(pi, "index2");
        type2 = read_int_dataset(pi, "type2");
        id2 = read_int_dataset(pi, "id2");
        if(int(row_type.size()) != n_elem ||
            int(row_residue_index.size()) != n_elem ||
            int(row_rotamer_index.size()) != n_elem) {
            throw string("cg_lipid_rotamer_sc row datasets must match source1 n_elem");
        }
        if(index2.size() != type2.size() || index2.size() != id2.size())
            throw string("cg_lipid_rotamer_sc source2 index/type/id size mismatch");

        row_index.resize(n_elem);
        for(int i = 0; i < n_elem; ++i)
            row_index[i] = i;

        row_group_count.assign(n_elem, 1);
        unordered_map<uint64_t, int> group_counts;
        group_counts.reserve(static_cast<size_t>(n_elem) * 2u);
        for(int i = 0; i < n_elem; ++i) {
            uint64_t key = (uint64_t(uint32_t(row_residue_index[i])) << 32)
                         | uint32_t(row_rotamer_index[i]);
            group_counts[key] += 1;
        }
        for(int i = 0; i < n_elem; ++i) {
            uint64_t key = (uint64_t(uint32_t(row_residue_index[i])) << 32)
                         | uint32_t(row_rotamer_index[i]);
            auto it = group_counts.find(key);
            row_group_count[i] = (it == group_counts.end() || it->second < 1) ? 1 : it->second;
        }
    }

    void update_box_dimensions_anisotropic(float scale_xy, float scale_z) {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        pairlist_valid = false;
    }

    void ensure_pairlist(VecArray sc, VecArray cg) {
        bool rebuild = !pairlist_valid || cache_buffer <= 0.f ||
            cached_box_x != box_x || cached_box_y != box_y || cached_box_z != box_z ||
            cached_positions_moved(sc, row_index, cached_pos1, cgl_pairlist_disp_threshold2(cache_buffer)) ||
            cached_positions_moved(cg, index2, cached_pos2, cgl_pairlist_disp_threshold2(cache_buffer));
        if(!rebuild)
            return;

        active_pairs.clear();
        float pairlist_cutoff = cutoff + cache_buffer;
        float pairlist_cutoff2 = pairlist_cutoff * pairlist_cutoff;
        for(int ai = 0; ai < n_elem; ++ai) {
            int t1 = row_type[ai];
            if(t1 < 0 || t1 >= n_type1) continue;
            for(size_t bi = 0; bi < index2.size(); ++bi) {
                int t2 = type2[bi];
                if(t2 < 0 || t2 >= n_type2) continue;
                float dr[3];
                center_dr(sc, ai, cg, index2[bi], box_x, box_y, box_z, dr);
                float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if(r2 < pairlist_cutoff2)
                    active_pairs.push_back({size_t(ai), bi});
            }
        }
        cache_positions(sc, row_index, cached_pos1);
        cache_positions(cg, index2, cached_pos2);
        cached_box_x = box_x;
        cached_box_y = box_y;
        cached_box_z = box_z;
        pairlist_valid = true;
    }

    virtual void compute_value(ComputeMode mode) override {
        (void)mode;
        fill(output, 0.f);
        VecArray sc = sc_pos.output;
        VecArray cg = cg_pos.output;
        ensure_pairlist(sc, cg);

        for(const CGLPairIndex& pair : active_pairs) {
            int ai = int(pair.first);
            size_t bi = pair.second;
            int t1 = row_type[ai];
            int t2 = type2[bi];

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(sc, ai, x1, n1);
            load_vec6(cg, index2[bi], x2, n2);
            center_dr(sc, ai, cg, index2[bi], box_x, box_y, box_z, dr);

            QuadsplineEval e;
            bool ok = eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e);
            if(ok) {
                if(log1p_reduced_transform)
                    apply_log1p_reduced_transform(
                            e, reference_energy_eup[t1 * n_type2 + t2],
                            boltzmann_temperature, dr, cutoff, taper_width);
                output(0, ai) += e.value;
            }
        }
        for(int ai = 0; ai < n_elem; ++ai)
            output(0, ai) /= float(std::max(1, row_group_count[ai]));
    }

    virtual void propagate_deriv() override {
        VecArray sc = sc_pos.output;
        VecArray sc_sens = sc_pos.sens;
        VecArray cg = cg_pos.output;
        VecArray cg_sens = cg_pos.sens;
        ensure_pairlist(sc, cg);

        for(const CGLPairIndex& pair : active_pairs) {
            int ai = int(pair.first);
            size_t bi = pair.second;
            float row_scale = sens(0, ai);
            if(row_scale == 0.f) continue;
            int t1 = row_type[ai];
            int t2 = type2[bi];
            row_scale /= float(std::max(1, row_group_count[ai]));

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(sc, ai, x1, n1);
            load_vec6(cg, index2[bi], x2, n2);
            center_dr(sc, ai, cg, index2[bi], box_x, box_y, box_z, dr);

            QuadsplineEval e;
            bool ok = eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e);
            if(!ok) continue;
            if(log1p_reduced_transform)
                apply_log1p_reduced_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, dr, cutoff, taper_width);

            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
            for(int d = 0; d < 3; ++d) {
                dpos1[d] *= row_scale;
                ddir1[d] *= row_scale;
                dpos2[d] *= row_scale;
                ddir2[d] *= row_scale;
            }
            add_vec6_sens(sc_sens, ai, dpos1, ddir1);
            add_vec6_sens(cg_sens, index2[bi], dpos2, ddir2);
        }
    }
};

CGLipidTargetPotential::CGLipidTargetPotential(
        hid_t grp, CoordNode& cg_pos_, CoordNode& tgt_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , tgt_pos(tgt_pos_)
    , n_type1(0)
    , n_type2(0)
    , n_param(0)
    , n_modes(read_attribute<int>(grp, ".", "n_modes", 0))
    , n_radial(read_attribute<int>(grp, ".", "n_radial", 14))
    , n_angular(read_attribute<int>(grp, ".", "n_angular", 15))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
        , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", 1.4f))
        , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                    float(n_radial - 2) * 1.4f))
        , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
        , boltzmann_weight_transform(read_attribute<int>(grp, ".", "boltzmann_weight_transform", 0) != 0)
        , log1p_reduced_transform(read_attribute<int>(grp, ".", "log1p_reduced_transform", 0) != 0)
        , boltzmann_temperature(read_attribute<float>(grp, ".", "boltzmann_temperature_upside", 0.f))
        , minimum_boltzmann_weight(read_attribute<float>(
                    grp, ".", "minimum_boltzmann_weight", numeric_limits<float>::min()))
        , cache_buffer(std::max(0.f, read_attribute<float>(grp, ".", "pairlist_buffer_ang", 1.f)))
        , body_support(read_cgl_body_support(grp))
        , cached_box_x(0.f)
        , cached_box_y(0.f)
        , cached_box_z(0.f)
        , pairlist_valid(false)
    {
    check_elem_width(cg_pos, 6);
    // Targets may be 3D particles or 6D backbone sites; only positions are used.

    H5Obj pi_obj = open_group(grp, "pair_interaction");
    hid_t pi = pi_obj.get();
    interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
    if(n_modes != 0 || n_radial <= 3 || n_angular <= 3 || n_param != n_radial * n_angular)
        throw string("cg_lipid_target requires tensor params with n_radial*n_angular coefficients");
    if(H5Lexists(pi, "reference_energy_eup", H5P_DEFAULT) > 0) {
        reference_energy_eup = read_float_dataset(pi, "reference_energy_eup");
        if(int(reference_energy_eup.size()) != n_type1 * n_type2)
            throw string("cg_lipid_target reference_energy_eup shape mismatch");
        if(!log1p_reduced_transform && !boltzmann_weight_transform)
            boltzmann_weight_transform = true;
    } else {
        reference_energy_eup.assign(n_type1 * n_type2, 0.f);
    }
    if(log1p_reduced_transform && boltzmann_weight_transform)
        throw string("cg_lipid_target cannot enable both Boltzmann-weight and log1p-reduced transforms");
    if((boltzmann_weight_transform || log1p_reduced_transform) && boltzmann_temperature <= 0.f)
        throw string("cg_lipid_target transformed PMF requires positive temperature");
    minimum_boltzmann_weight = std::max(minimum_boltzmann_weight, numeric_limits<float>::min());
    index1 = read_int_dataset(pi, "index1");
    type1 = read_int_dataset(pi, "type1");
    id1 = read_int_dataset(pi, "id1");
    index2 = read_int_dataset(pi, "index2");
    type2 = read_int_dataset(pi, "type2");
    id2 = read_int_dataset(pi, "id2");
    if(index1.size() != type1.size() || index1.size() != id1.size())
        throw string("cg_lipid_target source1 index/type/id size mismatch");
    if(index2.size() != type2.size() || index2.size() != id2.size())
        throw string("cg_lipid_target source2 index/type/id size mismatch");
}

void CGLipidTargetPotential::ensure_pairlist(VecArray cg, VecArray tgt) {
    bool rebuild = !pairlist_valid || cache_buffer <= 0.f ||
        cached_box_x != box_x || cached_box_y != box_y || cached_box_z != box_z ||
        cached_cgl_bodies_moved(cg, index1, cached_body1, body_support.axis_radius,
                cgl_pairlist_disp_threshold2(cache_buffer)) ||
        cached_positions_moved(tgt, index2, cached_pos2, cgl_pairlist_disp_threshold2(cache_buffer));
    if(!rebuild)
        return;

    active_pairs.clear();
    float pairlist_cutoff = cutoff + cache_buffer;
    float pairlist_cutoff2 = pairlist_cutoff * pairlist_cutoff;
    for(size_t ai = 0; ai < index1.size(); ++ai) {
        int t1 = type1[ai];
        if(t1 < 0 || t1 >= n_type1) continue;
        for(size_t bi = 0; bi < index2.size(); ++bi) {
            if((id1[ai] >> 4) == (id2[bi] >> 4)) continue;
            int t2 = type2[bi];
            if(t2 < 0 || t2 >= n_type2) continue;
            float dr[3];
            center_dr(cg, index1[ai], tgt, index2[bi], box_x, box_y, box_z, dr);
            if(cgl_target_pairlist_candidate(
                    cg, index1[ai], dr, body_support, cache_buffer, pairlist_cutoff2))
                active_pairs.push_back({ai, bi});
        }
    }
    cache_cgl_bodies(cg, index1, cached_body1);
    cache_positions(tgt, index2, cached_pos2);
    cached_box_x = box_x;
    cached_box_y = box_y;
    cached_box_z = box_z;
    pairlist_valid = true;
}

void CGLipidTargetPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray cg = cg_pos.output;
    VecArray tgt = tgt_pos.output;
    ensure_pairlist(cg, tgt);
    float total = 0.f;

    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type1[ai];
        int t2 = type2[bi];

        float x1[3], n1[3], dr[3];
        load_vec6(cg, index1[ai], x1, n1);
        center_dr(cg, index1[ai], tgt, index2[bi], box_x, box_y, box_z, dr);

        TargetSplineEval e;
        if(eval_cg_target_tensor(
            param_ptr(interaction_param, n_type2, n_param, t1, t2),
            n_angular, n_radial, dr, n1, knot_spacing, cutoff,
            (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e))
        {
            if(boltzmann_weight_transform) {
                apply_boltzmann_weight_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, minimum_boltzmann_weight,
                        dr, cutoff, taper_width);
            } else if(log1p_reduced_transform) {
                apply_log1p_reduced_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, dr, cutoff, taper_width);
            }
            total += e.value;
        }
    }
    potential = total;
}

void CGLipidTargetPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    VecArray tgt = tgt_pos.output;
    VecArray tgt_sens = tgt_pos.sens;
    ensure_pairlist(cg, tgt);

    for(const CGLPairIndex& pair : active_pairs) {
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type1[ai];
        int t2 = type2[bi];

        float x1[3], n1[3], dr[3];
        load_vec6(cg, index1[ai], x1, n1);
        center_dr(cg, index1[ai], tgt, index2[bi], box_x, box_y, box_z, dr);

        TargetSplineEval e;
        if(!eval_cg_target_tensor(
            param_ptr(interaction_param, n_type2, n_param, t1, t2),
            n_angular, n_radial, dr, n1, knot_spacing, cutoff,
            (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e))
            continue;
        if(boltzmann_weight_transform) {
            apply_boltzmann_weight_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, minimum_boltzmann_weight,
                    dr, cutoff, taper_width);
        } else if(log1p_reduced_transform) {
            apply_log1p_reduced_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, dr, cutoff, taper_width);
        }

        float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        float r = sqrtf(std::max(r2, 1e-12f));
        float inv_r = 1.f / r;
        float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};
        float a = n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2];
        float dpos[3] = {
            e.d_dr * unit[0] + e.d_da * (n1[0] - a * unit[0]) * inv_r,
            e.d_dr * unit[1] + e.d_da * (n1[1] - a * unit[1]) * inv_r,
            e.d_dr * unit[2] + e.d_da * (n1[2] - a * unit[2]) * inv_r,
        };
        float ddir[3] = {
            e.d_da * unit[0],
            e.d_da * unit[1],
            e.d_da * unit[2],
        };

        cg_sens(0, index1[ai]) -= dpos[0];
        cg_sens(1, index1[ai]) -= dpos[1];
        cg_sens(2, index1[ai]) -= dpos[2];
        cg_sens(3, index1[ai]) += ddir[0];
        cg_sens(4, index1[ai]) += ddir[1];
        cg_sens(5, index1[ai]) += ddir[2];
        Vec<3> tgt_grad = make_vec3(dpos[0], dpos[1], dpos[2]);
        bool projected = martini_hybrid::project_bb_proxy_gradient_for_coord(
            tgt_pos, tgt_sens, tgt_pos.n_elem, index2[bi], tgt_grad);
        if(!projected) {
            tgt_sens(0, index2[bi]) += dpos[0];
            tgt_sens(1, index2[bi]) += dpos[1];
            tgt_sens(2, index2[bi]) += dpos[2];
        }
    }
}

namespace martini_cg_lipid {

void register_dynamic_orientation_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed) {
    if(!engine) return;
    clear_dynamic_orientation_for_engine(engine);
    int idx = engine->get_idx("cgl_orientation_state", false);
    if(idx < 0) return;
    auto* state = dynamic_cast<CGLOrientationState*>(engine->nodes[idx].computation.get());
    if(!state) throw string("cgl_orientation_state node has unexpected type");
    auto runtime = std::unique_ptr<CGLDynamicOrientationRuntime>(
        new CGLDynamicOrientationRuntime(state, random_seed));

    if(h5_exists(config_root, "/input/cgl_orientation_mom")) {
        auto shape = get_dset_size(2, config_root, "/input/cgl_orientation_mom");
        if(shape[0] != hsize_t(state->n_elem) || shape[1] != 3u)
            throw string("/input/cgl_orientation_mom shape must be n_cgl x 3");
        traverse_dset<2, float>(config_root, "/input/cgl_orientation_mom",
                [&](size_t i, size_t d, float v) {
            runtime->angular_momentum(int(d), int(i)) = v;
        });
        for(int i = 0; i < state->n_elem; ++i) {
            float3 n = normalized(load_vec<3>(state->output, i));
            float3 l = project_tangent(load_vec<3>(runtime->angular_momentum, i), n);
            store_vec(runtime->angular_momentum, i, l);
        }
    }
    orientation_runtime_map()[engine] = std::move(runtime);
}

void clear_dynamic_orientation_for_engine(DerivEngine* engine) {
    orientation_runtime_map().erase(engine);
}

bool has_dynamic_orientation(DerivEngine* engine) {
    return get_orientation_runtime(engine) != nullptr;
}

void set_dynamic_orientation_temperature(DerivEngine* engine, float temperature) {
    auto* runtime = get_orientation_runtime(engine);
    if(runtime) runtime->temperature = temperature;
}

void set_dynamic_orientation_thermostat_delta_t(DerivEngine* engine, float delta_t) {
    auto* runtime = get_orientation_runtime(engine);
    if(runtime) runtime->thermostat_delta_t = delta_t;
}

void apply_dynamic_orientation_thermostat(DerivEngine* engine) {
    auto* runtime = get_orientation_runtime(engine);
    if(!runtime || runtime->thermostat_delta_t <= 0.f) return;
    auto* state = runtime->state;
    float tau = state->thermostat_timescale;
    float mom_scale = expf(-runtime->thermostat_delta_t / tau);
    float noise_scale = sqrtf(
        state->rotational_inertia * runtime->temperature *
        std::max(0.f, 1.f - mom_scale * mom_scale));
    for(int i = 0; i < state->n_elem; ++i) {
        RandomGenerator random(
            runtime->random_seed,
            CGL_ORIENTATION_THERMOSTAT_RANDOM_STREAM,
            uint32_t(i),
            runtime->thermostat_invocations);
        float3 n = normalized(load_vec<3>(state->output, i));
        float3 l = project_tangent(load_vec<3>(runtime->angular_momentum, i), n);
        float3 kick = project_tangent(random.normal3(), n) * noise_scale;
        l = project_tangent(l * mom_scale + kick, n);
        store_vec(runtime->angular_momentum, i, l);
    }
    runtime->thermostat_invocations++;
}

void integrate_dynamic_orientation(DerivEngine* engine, float dt) {
    auto* runtime = get_orientation_runtime(engine);
    if(!runtime) return;
    auto* state = runtime->state;
    for(int i = 0; i < state->n_elem; ++i) {
        float3 n = normalized(load_vec<3>(state->output, i));
        float3 grad = load_vec<3>(state->sens, i);
        float3 torque = -cross(n, grad);
        float3 l = project_tangent(load_vec<3>(runtime->angular_momentum, i) + torque * dt, n);
        float3 omega = l / state->rotational_inertia;
        n = rotate_unit_vector(n, omega, dt);
        l = project_tangent(l, n);
        store_vec(state->output, i, n);
        store_vec(runtime->angular_momentum, i, l);
    }
}

void add_dynamic_orientation_loggers(DerivEngine* engine, H5Logger& logger, bool record_momentum) {
    auto* runtime = get_orientation_runtime(engine);
    if(!runtime) return;
    auto* state = runtime->state;
    logger.add_logger<float>("cgl_orientation", {1, state->n_elem, 3}, [state](float* buffer) {
        VecArray orient = state->output;
        for(int i = 0; i < state->n_elem; ++i)
            for(int d = 0; d < 3; ++d)
                buffer[i * 3 + d] = orient(d, i);
    });
    if(record_momentum) {
        logger.add_logger<float>("cgl_orientation_mom", {1, state->n_elem, 3}, [runtime](float* buffer) {
            VecArray mom = runtime->angular_momentum;
            for(int i = 0; i < runtime->state->n_elem; ++i)
            for(int d = 0; d < 3; ++d)
                buffer[i * 3 + d] = mom(d, i);
        });
    }
}

void register_cgl_gle_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed) {
    if(!engine) return;
    cgl_gle_runtime_map().erase(engine);
    if(!h5_exists(config_root, "/input/cgl_gle")) return;

    H5Obj gle_obj = open_group(config_root, "/input/cgl_gle");
    hid_t grp = gle_obj.get();

    vector<int> atom_index = read_int_dataset(grp, "atom_index");
    if(atom_index.empty())
        throw string("/input/cgl_gle/atom_index must contain at least one CGL atom");
    vector<float> memory_tau;
    vector<float> coupling;
    if(H5Lexists(grp, "memory_tau", H5P_DEFAULT) > 0) {
        memory_tau = read_float_dataset_1d(grp, "memory_tau");
    } else {
        memory_tau.push_back(read_attribute<float>(grp, ".", "memory_tau", 0.f));
    }
    if(H5Lexists(grp, "coupling", H5P_DEFAULT) > 0) {
        coupling = read_float_dataset_1d(grp, "coupling");
    } else {
        coupling.push_back(read_attribute<float>(grp, ".", "coupling", 0.f));
    }
    if(memory_tau.empty() || coupling.empty() || memory_tau.size() != coupling.size())
        throw string("/input/cgl_gle memory_tau and coupling must have the same positive length");
    for(size_t m = 0; m < memory_tau.size(); ++m) {
        if(memory_tau[m] <= 0.f || coupling[m] <= 0.f)
            throw string("/input/cgl_gle requires positive memory_tau and coupling values");
    }

    auto runtime = std::unique_ptr<CGLGLERuntime>(
        new CGLGLERuntime(atom_index, random_seed, memory_tau, coupling));
    if(H5Lexists(grp, "aux_momentum", H5P_DEFAULT) > 0) {
        int n_mode = runtime->n_mode();
        if(n_mode == 1) {
            auto shape = get_dset_size(2, grp, "aux_momentum");
            if(shape[0] != hsize_t(atom_index.size()) || shape[1] != 3u)
                throw string("/input/cgl_gle/aux_momentum shape must be n_cgl x 3 for single-mode GLE");
            traverse_dset<2, float>(grp, "aux_momentum", [&](size_t i, size_t d, float v) {
                runtime->aux_momentum(int(d), runtime->aux_index(0, int(i))) = v;
            });
        } else {
            auto shape = get_dset_size(3, grp, "aux_momentum");
            if(shape[0] != hsize_t(n_mode) || shape[1] != hsize_t(atom_index.size()) || shape[2] != 3u)
                throw string("/input/cgl_gle/aux_momentum shape must be n_mode x n_cgl x 3 for multi-mode GLE");
            traverse_dset<3, float>(grp, "aux_momentum", [&](size_t m, size_t i, size_t d, float v) {
                runtime->aux_momentum(int(d), runtime->aux_index(int(m), int(i))) = v;
            });
        }
    }
    cgl_gle_runtime_map()[engine] = std::move(runtime);
}

void clear_cgl_gle_for_engine(DerivEngine* engine) {
    cgl_gle_runtime_map().erase(engine);
}

bool has_cgl_gle(DerivEngine* engine) {
    return get_cgl_gle_runtime(engine) != nullptr;
}

void set_cgl_gle_temperature(DerivEngine* engine, float temperature) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(runtime) runtime->temperature = temperature;
}

void set_cgl_gle_delta_t(DerivEngine* engine, float delta_t) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(runtime) runtime->delta_t = delta_t;
}

void apply_cgl_gle_thermostat(DerivEngine* engine, VecArray mom) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(!runtime || runtime->delta_t <= 0.f) return;

    bool have_masses = martini_masses::has_masses(engine);

    for(size_t i = 0; i < runtime->atom_index.size(); ++i) {
        int atom = runtime->atom_index[i];
        if(atom < 0 || atom >= engine->pos->n_atom) continue;
        float mass = have_masses ? martini_masses::get_mass(engine, atom) : 1.f;
        float3 p = load_vec<3>(mom, atom);
        float sqrt_mass = (mass > 0.f) ? sqrtf(mass) : 1.f;
        for(int mode = 0; mode < runtime->n_mode(); ++mode) {
            float angle = runtime->coupling[mode] * runtime->delta_t;
            float c = cosf(angle);
            float s = sinf(angle);
            float ou_scale = expf(-runtime->delta_t / runtime->memory_tau[mode]);
            float noise_scale = sqrtf(
                runtime->temperature * std::max(0.f, 1.f - ou_scale * ou_scale)) * sqrt_mass;
            uint32_t random_atom = uint32_t(atom) ^ (uint32_t(mode + 1) * 0x9e3779b9u);
            RandomGenerator random(
                runtime->random_seed,
                CGL_GLE_THERMOSTAT_RANDOM_STREAM,
                random_atom,
                runtime->thermostat_invocations);

            int aux_idx = runtime->aux_index(mode, int(i));
            float3 a = load_vec<3>(runtime->aux_momentum, aux_idx);
            float3 p_rot = p * c - a * s;
            float3 a_rot = p * s + a * c;
            float3 a_next = a_rot * ou_scale + random.normal3() * noise_scale;
            p = p_rot;
            store_vec(runtime->aux_momentum, aux_idx, a_next);
        }
        store_vec(mom, atom, p);
    }
    runtime->thermostat_invocations++;
}

void add_cgl_gle_loggers(DerivEngine* engine, H5Logger& logger, bool record_aux) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(!runtime || !record_aux) return;
    logger.add_logger<float>("cgl_gle_aux_momentum", {runtime->n_mode(), int(runtime->atom_index.size()), 3},
            [runtime](float* buffer) {
        VecArray aux = runtime->aux_momentum;
        int n_cgl = int(runtime->atom_index.size());
        for(int mode = 0; mode < runtime->n_mode(); ++mode)
            for(int i = 0; i < n_cgl; ++i)
                for(int d = 0; d < 3; ++d)
                    buffer[(mode * n_cgl + i) * 3 + d] =
                        aux(d, runtime->aux_index(mode, i));
    });
}

}  // namespace martini_cg_lipid

static RegisterNodeType<CGLOrientationState, 0> _reg_cg_orientation("cgl_orientation_state");
static RegisterNodeType<ComposeVector6D, -1> _reg_cv("compose_vector6d");
static RegisterNodeType<CGLipidPairPotential, 1> _reg_cg_pair("cg_lipid_pair");
static RegisterNodeType<CGLipidDensityPotential, 1> _reg_cg_density("cg_lipid_density");
static RegisterNodeType<CGLipidSCPotential, 2> _reg_cg_sc("cg_lipid_sc");
static RegisterNodeType<CGLipidSCOneBody, 2> _reg_cg_sc_1body("cg_lipid_rotamer_sc");
static RegisterNodeType<CGLipidTargetPotential, 2> _reg_cg_target("cg_lipid_target");
