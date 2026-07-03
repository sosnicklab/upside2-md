#include "deriv_engine.h"
#include "spline.h"
#include "box.h"
#include "martini.h"
#include "random.h"
#include "state_logger.h"
#include <H5Apublic.h> // for H5Aexists

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

struct ImplicitCompactionField {
    vector<float> survival;
    vector<float> contact_field;
    vector<float> contact_slope;
    vector<float> response_slope;
    vector<float> extended;
    vector<float> pair_weight;
    vector<float> pair_distance;
    vector<float> gap_coord;
    vector<float> gap_weight_sum;
    vector<uint8_t> leaflet_side;
    vector<int> dominant_pair;
};

static inline float implicit_compaction_contact_field_from_survival(
        float survival,
        float contact_field_cap) {
    float clamped_survival = std::max(survival, 1.0e-6f);
    float contact_field = -logf(clamped_survival);
    if(contact_field_cap > 0.f)
        contact_field = std::min(contact_field, contact_field_cap);
    return contact_field;
}

static inline float implicit_compaction_probability_from_field(
        float contact_field,
        float isolated_compact_probability,
        float local_field_scale) {
    if(local_field_scale <= 0.f)
        return 1.f - expf(-contact_field);

    if(isolated_compact_probability <= 0.f)
        return 1.f - expf(-local_field_scale * contact_field);

    float p0 = std::max(1.0e-6f, std::min(1.f - 1.0e-6f, isolated_compact_probability));
    float compact_logit = logf(p0 / (1.f - p0)) + local_field_scale * contact_field;
    return 1.f / (1.f + expf(-compact_logit));
}

static inline float implicit_compaction_probability_field_slope(
        float contact_field,
        float isolated_compact_probability,
        float local_field_scale) {
    if(local_field_scale <= 0.f)
        return expf(-contact_field);

    if(isolated_compact_probability <= 0.f)
        return local_field_scale * expf(-local_field_scale * contact_field);

    float probability = implicit_compaction_probability_from_field(
            contact_field,
            isolated_compact_probability,
            local_field_scale);
    return local_field_scale * probability * (1.f - probability);
}

struct MeanFieldPairControl {
    float base_control;
    float corr_ee;
    float corr_ec;
    float corr_cc;
    float distance;
    bool valid;
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

struct TemperatureScaleTable {
    vector<float> value;
    int n_mode = 1;
};

static TemperatureScaleTable make_uniform_temperature_scale(size_t n_temp, float scale) {
    TemperatureScaleTable out;
    out.value.assign(n_temp, scale);
    out.n_mode = 1;
    return out;
}

static TemperatureScaleTable read_temperature_scale_table(
        hid_t grp, const char* name, size_t expected_n_temp, size_t expected_n_mode) {
    H5Obj dset = h5_obj(H5Dclose, H5Dopen2(grp, name, H5P_DEFAULT));
    H5Obj space = h5_obj(H5Sclose, H5Dget_space(dset.get()));
    int ndims = H5Sget_simple_extent_ndims(space.get());
    if(ndims < 0)
        throw string("hdf5 error while reading /input/cgl_gle/") + name;

    TemperatureScaleTable out;
    if(ndims == 1) {
        out.value = read_float_dataset_1d(grp, name);
        if(out.value.size() != expected_n_temp)
            throw string("/input/cgl_gle ") + name + " length must match temperature_grid length";
        out.n_mode = 1;
        return out;
    }
    if(ndims == 2) {
        vector<hsize_t> sz = get_dset_size(2, grp, name);
        if(sz[0] != expected_n_temp || sz[1] != expected_n_mode)
            throw string("/input/cgl_gle ") + name + " shape must be n_temp x n_mode";
        out.value = read_float_dataset(grp, name);
        out.n_mode = int(sz[1]);
        return out;
    }
    throw string("/input/cgl_gle ") + name + " must be a 1D or 2D float dataset";
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

static vector<float> read_named_param_dataset_any(
        hid_t grp,
        const char* name,
        int expected_type1,
        int expected_type2,
        int expected_param) {
    vector<hsize_t> sz = get_dset_size(3, grp, name);
    if(int(sz[0]) != expected_type1
            || int(sz[1]) != expected_type2
            || int(sz[2]) != expected_param)
        throw string(name) + " shape mismatch";
    vector<float> out(size_t(expected_type1) * size_t(expected_type2) * size_t(expected_param), 0.f);
    traverse_dset<3, float>(grp, name, [&](size_t i, size_t j, size_t k, float v) {
        out[(i * expected_type2 + j) * expected_param + k] = v;
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

static inline Vec<2> extrapolated_deBoor_value_and_deriv(
        const float* bspline_coeff, const float x, int n_knot) {
    constexpr float eps = 1.0e-4f;
    float left_x = 1.f + eps;
    float right_x = float(n_knot - 2) - eps;
    if(x <= 1.f) {
        Vec<2> v = deBoor_value_and_deriv(bspline_coeff, left_x);
        return make_vec2(v.x() + (x - left_x) * v.y(), v.y());
    }
    if(x >= float(n_knot - 2)) {
        Vec<2> v = deBoor_value_and_deriv(bspline_coeff, right_x);
        return make_vec2(v.x() + (x - right_x) * v.y(), v.y());
    }
    return deBoor_value_and_deriv(bspline_coeff, x);
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

static inline float clamp_probability(float p) {
    return std::max(1.0e-5f, std::min(1.0f - 1.0e-5f, p));
}

static inline float logistic(float x) {
    if(x >= 0.f) {
        float e = expf(-x);
        return 1.f / (1.f + e);
    } else {
        float e = expf(x);
        return e / (1.f + e);
    }
}

static inline float logit(float p) {
    p = clamp_probability(p);
    return logf(p / (1.f - p));
}

static inline float binary_state_relative_entropy(
        float compact_probability,
        float isolated_compact_probability,
        float boltzmann_temperature) {
    float s = clamp_probability(compact_probability);
    float p0 = clamp_probability(isolated_compact_probability);
    return boltzmann_temperature * (
            s * logf(s / p0)
            + (1.f - s) * logf((1.f - s) / (1.f - p0)));
}

static inline void compute_compaction_weights(
        float compact_i,
        float compact_j,
        float& w_ee,
        float& w_ec,
        float& w_cc) {
    w_ee = (1.f - compact_i) * (1.f - compact_j);
    w_ec = compact_i * (1.f - compact_j) + (1.f - compact_i) * compact_j;
    w_cc = compact_i * compact_j;
}

static inline float transformed_pair_scale_from_control(
        float control_value,
        float boltzmann_temperature,
        float distance,
        float cutoff,
        float taper_width) {
    float taper = 1.f;
    float d_taper_dr = 0.f;
    compute_cutoff_taper(distance, cutoff, taper_width, taper, d_taper_dr);
    return boltzmann_temperature * expf(control_value) * taper;
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

static vector<uint8_t> classify_leaflets_by_median_z(
        VecArray cg, const vector<int>& index) {
    vector<float> z_values;
    z_values.reserve(index.size());
    for(int idx: index) z_values.push_back(cg(2, idx));
    vector<float> sorted_z = z_values;
    auto mid_it = sorted_z.begin() + sorted_z.size() / 2u;
    nth_element(sorted_z.begin(), mid_it, sorted_z.end());
    float mid_z = *mid_it;

    vector<uint8_t> side(index.size(), 0u);
    for(size_t i = 0; i < index.size(); ++i)
        side[i] = uint8_t(z_values[i] > mid_z);
    return side;
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

static inline float cgl_cross_leaflet_face_weight(
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float radial_cutoff_ang,
        float face_cos_min,
        bool smooth) {
    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f) return 0.f;
    float r = sqrtf(r2);
    if(r > radial_cutoff_ang) return 0.f;
    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};
    float a1 = -(n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2]);
    float a2 =  (n2[0] * unit[0] + n2[1] * unit[1] + n2[2] * unit[2]);
    float face1 = -a1;
    float face2 = -a2;
    if(!smooth)
        return (face1 >= face_cos_min && face2 >= face_cos_min) ? 1.f : 0.f;
    float ang_denom = std::max(1.f - face_cos_min, 1e-6f);
    float w1 = clamp_range((face1 - face_cos_min) / ang_denom, 0.f, 1.f);
    float w2 = clamp_range((face2 - face_cos_min) / ang_denom, 0.f, 1.f);
    float wr = clamp_range(1.f - r / std::max(radial_cutoff_ang, 1e-6f), 0.f, 1.f);
    return wr * w1 * w2;
}

static bool eval_cgl_cross_leaflet_face_weight(
        const float dr[3],
        const float n1[3],
        const float n2[3],
        float radial_cutoff_ang,
        float face_cos_min,
        bool smooth,
        QuadsplineEval& out) {
    out.value = 0.f;
    out.d_dr = 0.f;
    out.d_da1 = 0.f;
    out.d_da2 = 0.f;
    if(!smooth)
        return false;
    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if(r2 <= 1e-12f)
        return false;
    float r = sqrtf(r2);
    if(r >= radial_cutoff_ang)
        return false;
    float inv_r = 1.f / r;
    float unit[3] = {dr[0] * inv_r, dr[1] * inv_r, dr[2] * inv_r};
    float a1 = -(n1[0] * unit[0] + n1[1] * unit[1] + n1[2] * unit[2]);
    float a2 =  (n2[0] * unit[0] + n2[1] * unit[1] + n2[2] * unit[2]);
    float face1 = -a1;
    float face2 = -a2;
    float ang_denom = std::max(1.f - face_cos_min, 1e-6f);
    float w1 = clamp_range((face1 - face_cos_min) / ang_denom, 0.f, 1.f);
    float w2 = clamp_range((face2 - face_cos_min) / ang_denom, 0.f, 1.f);
    float wr = clamp_range(1.f - r / std::max(radial_cutoff_ang, 1e-6f), 0.f, 1.f);
    if(w1 <= 0.f || w2 <= 0.f || wr <= 0.f)
        return false;
    float dw1_da1 = (w1 < 1.f) ? (-1.f / ang_denom) : 0.f;
    float dw2_da2 = (w2 < 1.f) ? (-1.f / ang_denom) : 0.f;
    float dwr_dr = (wr < 1.f) ? (-1.f / std::max(radial_cutoff_ang, 1e-6f)) : 0.f;
    out.value = wr * w1 * w2;
    out.d_dr = dwr_dr * w1 * w2;
    out.d_da1 = wr * dw1_da1 * w2;
    out.d_da2 = wr * w1 * dw2_da2;
    return true;
}

static ImplicitCompactionField build_gap_response_compaction_field(
        VecArray cg,
        const vector<int>& index,
        const vector<CGLPairIndex>& active_pairs,
        float box_x,
        float box_y,
        float box_z,
        const vector<float>& gap_response_coeff,
        int gap_response_n_knot,
        float gap_response_coord_min_ang,
        float gap_response_coord_spacing_ang,
        float gap_response_radial_cutoff_ang,
        float gap_response_face_cos_min,
        float gap_response_fallback_ang,
        bool gap_response_smooth_weight) {
    ImplicitCompactionField field;
    field.survival.assign(index.size(), 1.f);
    field.contact_field.assign(index.size(), 0.f);
    field.contact_slope.assign(index.size(), 0.f);
    field.response_slope.assign(index.size(), 0.f);
    field.extended.assign(index.size(), 0.f);
    field.pair_weight.assign(active_pairs.size(), 0.f);
    field.pair_distance.assign(active_pairs.size(), 0.f);
    field.gap_coord.assign(index.size(), gap_response_fallback_ang > 0.f ? gap_response_fallback_ang : 0.f);
    field.gap_weight_sum.assign(index.size(), 0.f);
    if(index.empty())
        return field;

    field.leaflet_side = classify_leaflets_by_median_z(cg, index);
    field.dominant_pair.assign(index.size(), -1);
    vector<float> gap_weighted_sum(index.size(), 0.f);
    for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
        const CGLPairIndex& pair = active_pairs[pidx];
        size_t ai = pair.first;
        size_t bi = pair.second;
        if(field.leaflet_side[ai] == field.leaflet_side[bi]) continue;
        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
        float w = cgl_cross_leaflet_face_weight(
                dr, n1, n2,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_smooth_weight);
        if(w <= 0.f) continue;
        float r = norm_dr(dr);
        field.pair_weight[pidx] = w;
        field.pair_distance[pidx] = r;
        gap_weighted_sum[ai] += w * r;
        gap_weighted_sum[bi] += w * r;
        field.gap_weight_sum[ai] += w;
        field.gap_weight_sum[bi] += w;
    }
    for(size_t i = 0; i < field.extended.size(); ++i) {
        float gap = gap_response_fallback_ang > 0.f ? gap_response_fallback_ang : gap_response_radial_cutoff_ang;
        if(field.gap_weight_sum[i] > 1.0e-12f)
            gap = gap_weighted_sum[i] / field.gap_weight_sum[i];
        field.gap_coord[i] = gap;
        if(gap_response_coeff.empty() || gap_response_n_knot <= 3 || gap_response_coord_spacing_ang <= 0.f)
            continue;
        float coord = 1.f + (gap - gap_response_coord_min_ang) / gap_response_coord_spacing_ang;
        Vec<2> response = clamped_deBoor_value_and_deriv(
                gap_response_coeff.data(),
                coord,
                gap_response_n_knot);
        field.extended[i] = clamp_probability(response.x());
        field.response_slope[i] = response.y() / gap_response_coord_spacing_ang;
    }
    return field;
}

static void propagate_gap_response_compaction_sens(
        VecArray cg,
        VecArray cg_sens,
        const vector<int>& index,
        const vector<CGLPairIndex>& active_pairs,
        float box_x,
        float box_y,
        float box_z,
        const ImplicitCompactionField& field,
        const vector<float>& implicit_q_sens,
        float gap_response_radial_cutoff_ang,
        float gap_response_face_cos_min,
        bool gap_response_smooth_weight) {
    if(index.empty() || field.leaflet_side.empty())
        return;
    for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
        const CGLPairIndex& pair = active_pairs[pidx];
        size_t ai = pair.first;
        size_t bi = pair.second;
        if(field.leaflet_side[ai] == field.leaflet_side[bi]) continue;
        float w = field.pair_weight[pidx];
        if(w <= 0.f) continue;
        float r = field.pair_distance[pidx];
        float lambda_r = 0.f;
        float lambda_w = 0.f;
        if(field.gap_weight_sum[ai] > 1.0e-12f) {
            float pref = implicit_q_sens[ai] * field.response_slope[ai]
                / field.gap_weight_sum[ai];
            lambda_r += pref * w;
            lambda_w += pref * (r - field.gap_coord[ai]);
        }
        if(field.gap_weight_sum[bi] > 1.0e-12f) {
            float pref = implicit_q_sens[bi] * field.response_slope[bi]
                / field.gap_weight_sum[bi];
            lambda_r += pref * w;
            lambda_w += pref * (r - field.gap_coord[bi]);
        }
        if(fabsf(lambda_r) <= 1.0e-12f
                && (!gap_response_smooth_weight || fabsf(lambda_w) <= 1.0e-12f))
            continue;

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);

        if(fabsf(lambda_r) > 1.0e-12f) {
            QuadsplineEval r_eval = {0.f, lambda_r, 0.f, 0.f};
            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, r_eval, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
        }

        if(gap_response_smooth_weight && fabsf(lambda_w) > 1.0e-12f) {
            QuadsplineEval w_eval;
            if(!eval_cgl_cross_leaflet_face_weight(
                        dr, n1, n2,
                        gap_response_radial_cutoff_ang,
                        gap_response_face_cos_min,
                        gap_response_smooth_weight,
                        w_eval))
                continue;
            w_eval.d_dr *= lambda_w;
            w_eval.d_da1 *= lambda_w;
            w_eval.d_da2 *= lambda_w;
            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, w_eval, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
        }
    }
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

struct CGLCompactionState : public CoordNode {
    float mass;
    float thermostat_timescale;
    float coord_min;
    float coord_max;

    CGLCompactionState(hid_t grp)
        : CoordNode(int(h5::get_dset_size(1, grp, "value")[0]), 1)
        , mass(read_attribute<float>(grp, ".", "mass_up", 1.f))
        , thermostat_timescale(read_attribute<float>(grp, ".", "thermostat_timescale", 5.f))
        , coord_min(read_attribute<float>(grp, ".", "self_coord_min_ang", -1e6f))
        , coord_max(read_attribute<float>(grp, ".", "self_coord_max_ang", 1e6f))
    {
        if(mass <= 0.f)
            throw string("cgl_compaction_state mass_up must be positive");
        if(thermostat_timescale <= 0.f)
            throw string("cgl_compaction_state thermostat_timescale must be positive");
        if(!(coord_max > coord_min))
            throw string("cgl_compaction_state requires self_coord_max_ang > self_coord_min_ang");

        traverse_dset<1, float>(grp, "value", [&](size_t i, float v) {
            output(0, int(i)) = std::max(coord_min, std::min(coord_max, v));
        });
    }

    virtual void compute_value(ComputeMode) override {}
    virtual void propagate_deriv() override {}
};

namespace {

struct CGLDynamicOrientationRuntime {
    CGLOrientationState* state = nullptr;
    VecArrayStorage angular_momentum;
    float temperature = 1.f;
    float hidden_state_temperature = -1.f;
    float thermostat_delta_t = 0.f;
    uint32_t random_seed = 0u;
    uint64_t thermostat_invocations = 0u;
    bool hidden_state_initialized = false;
    bool temperature_is_set = false;

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

struct CGLDynamicCompactionRuntime {
    CGLCompactionState* state = nullptr;
    VecArrayStorage momentum;
    float temperature = 1.f;
    float hidden_state_temperature = -1.f;
    float thermostat_delta_t = 0.f;
    uint32_t random_seed = 0u;
    uint64_t thermostat_invocations = 0u;
    bool hidden_state_initialized = false;
    bool temperature_is_set = false;

    CGLDynamicCompactionRuntime(CGLCompactionState* state_, uint32_t random_seed_)
        : state(state_)
        , momentum(1, state_ ? state_->n_elem : 1)
        , random_seed(random_seed_)
    {
        fill(momentum, 0.f);
    }
};

static std::map<DerivEngine*, std::unique_ptr<CGLDynamicCompactionRuntime>>& compaction_runtime_map() {
    static std::map<DerivEngine*, std::unique_ptr<CGLDynamicCompactionRuntime>> m;
    return m;
}

static CGLDynamicCompactionRuntime* get_compaction_runtime(DerivEngine* engine) {
    auto& m = compaction_runtime_map();
    auto it = m.find(engine);
    return it == m.end() ? nullptr : it->second.get();
}

struct CGLGLERuntime {
    vector<int> atom_index;
    vector<float> base_memory_tau;
    vector<float> memory_tau;
    vector<float> base_coupling;
    vector<float> coupling;
    vector<float> temperature_grid;
    TemperatureScaleTable coupling_scale;
    TemperatureScaleTable memory_tau_scale;
    VecArrayStorage aux_momentum;
    float temperature = 1.f;
    float hidden_state_temperature = -1.f;
    float delta_t = 0.f;
    uint32_t random_seed = 0u;
    uint64_t thermostat_invocations = 0u;
    vector<double> linear_step_matrix;
    vector<double> noise_cholesky_unit;
    bool hidden_state_initialized = false;
    bool temperature_is_set = false;
    bool linear_step_ready = false;

    CGLGLERuntime(
            const vector<int>& atom_index_,
            uint32_t random_seed_,
            const vector<float>& memory_tau_,
            const vector<float>& coupling_,
            const vector<float>& temperature_grid_,
            const TemperatureScaleTable& coupling_scale_,
            const TemperatureScaleTable& memory_tau_scale_)
        : atom_index(atom_index_)
        , base_memory_tau(memory_tau_)
        , memory_tau(memory_tau_)
        , base_coupling(coupling_)
        , coupling(coupling_)
        , temperature_grid(temperature_grid_)
        , coupling_scale(coupling_scale_)
        , memory_tau_scale(memory_tau_scale_)
        , aux_momentum(3, int(std::max<size_t>(1, atom_index_.size() * memory_tau_.size())))
        , random_seed(random_seed_)
    {
        fill(aux_momentum, 0.f);
    }

    int n_mode() const { return int(memory_tau.size()); }
    int linear_step_dim() const { return n_mode() + 1; }
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

static inline bool attribute_exists_local(hid_t loc_id, const char* obj_name, const char* attr_name) {
    if(!obj_name || obj_name[0] == '\0' || (obj_name[0] == '.' && obj_name[1] == '\0')) {
        htri_t exists = H5Aexists(loc_id, attr_name);
        return exists > 0;
    }

    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if(obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

static inline double& dense_square_at(vector<double>& m, int n, int row, int col) {
    return m[size_t(row) * size_t(n) + size_t(col)];
}

static inline double dense_square_at(const vector<double>& m, int n, int row, int col) {
    return m[size_t(row) * size_t(n) + size_t(col)];
}

static vector<double> dense_square_identity(int n) {
    vector<double> out(size_t(n) * size_t(n), 0.0);
    for(int i = 0; i < n; ++i)
        dense_square_at(out, n, i, i) = 1.0;
    return out;
}

static vector<double> dense_square_multiply(const vector<double>& a, const vector<double>& b, int n) {
    vector<double> out(size_t(n) * size_t(n), 0.0);
    for(int i = 0; i < n; ++i) {
        for(int k = 0; k < n; ++k) {
            double aik = dense_square_at(a, n, i, k);
            if(aik == 0.0) continue;
            for(int j = 0; j < n; ++j)
                dense_square_at(out, n, i, j) += aik * dense_square_at(b, n, k, j);
        }
    }
    return out;
}

static vector<double> dense_square_transpose(const vector<double>& a, int n) {
    vector<double> out(size_t(n) * size_t(n), 0.0);
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            dense_square_at(out, n, j, i) = dense_square_at(a, n, i, j);
    return out;
}

static double dense_square_max_abs(const vector<double>& a) {
    double out = 0.0;
    for(double x : a)
        out = std::max(out, fabs(x));
    return out;
}

static vector<double> dense_square_exponential(const vector<double>& a, int n) {
    vector<double> result = dense_square_identity(n);
    vector<double> term = dense_square_identity(n);
    for(int k = 1; k <= 128; ++k) {
        vector<double> next = dense_square_multiply(term, a, n);
        double inv_k = 1.0 / double(k);
        for(size_t i = 0; i < next.size(); ++i)
            next[i] *= inv_k;
        term.swap(next);
        for(size_t i = 0; i < result.size(); ++i)
            result[i] += term[i];
        if(dense_square_max_abs(term) < 1.0e-14)
            break;
    }
    return result;
}

static bool dense_square_cholesky_lower(vector<double>& a, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j <= i; ++j) {
            double sum = dense_square_at(a, n, i, j);
            for(int k = 0; k < j; ++k)
                sum -= dense_square_at(a, n, i, k) * dense_square_at(a, n, j, k);

            if(i == j) {
                if(sum < -1.0e-10)
                    return false;
                dense_square_at(a, n, i, j) = sqrt(std::max(sum, 1.0e-12));
            } else {
                dense_square_at(a, n, i, j) = sum / dense_square_at(a, n, j, j);
            }
        }
        for(int j = i + 1; j < n; ++j)
            dense_square_at(a, n, i, j) = 0.0;
    }
    return true;
}

static void update_cgl_gle_linear_step(CGLGLERuntime* runtime) {
    if(!runtime) return;
    // Build the exact linear-Gaussian step for the coupled physical momentum
    // plus GLE auxiliary momenta using the Van Loan block exponential.

    int n_mode = runtime->n_mode();
    int n_state = runtime->linear_step_dim();
    runtime->linear_step_matrix.assign(size_t(n_state) * size_t(n_state), 0.0);
    runtime->noise_cholesky_unit.assign(size_t(n_state) * size_t(n_state), 0.0);
    runtime->linear_step_ready = false;

    if(runtime->delta_t <= 0.f || n_mode <= 0)
        return;

    vector<double> drift(size_t(n_state) * size_t(n_state), 0.0);
    for(int mode = 0; mode < n_mode; ++mode) {
        double coupling = runtime->coupling[size_t(mode)];
        double inv_tau = 1.0 / std::max(double(runtime->memory_tau[size_t(mode)]), 1.0e-12);
        dense_square_at(drift, n_state, 0, mode + 1) = -coupling;
        dense_square_at(drift, n_state, mode + 1, 0) = coupling;
        dense_square_at(drift, n_state, mode + 1, mode + 1) = -inv_tau;
    }

    vector<double> diffusion(size_t(n_state) * size_t(n_state), 0.0);
    for(int mode = 0; mode < n_mode; ++mode) {
        double inv_tau = 1.0 / std::max(double(runtime->memory_tau[size_t(mode)]), 1.0e-12);
        dense_square_at(diffusion, n_state, mode + 1, mode + 1) = 2.0 * inv_tau;
    }

    int block_dim = 2 * n_state;
    vector<double> van_loan(size_t(block_dim) * size_t(block_dim), 0.0);
    for(int i = 0; i < n_state; ++i) {
        for(int j = 0; j < n_state; ++j) {
            double scaled_drift = double(runtime->delta_t) * dense_square_at(drift, n_state, i, j);
            dense_square_at(van_loan, block_dim, i, j) = scaled_drift;
            dense_square_at(van_loan, block_dim, i, j + n_state) =
                double(runtime->delta_t) * dense_square_at(diffusion, n_state, i, j);
            dense_square_at(van_loan, block_dim, i + n_state, j + n_state) =
                -double(runtime->delta_t) * dense_square_at(drift, n_state, j, i);
        }
    }

    vector<double> block_exp = dense_square_exponential(van_loan, block_dim);
    vector<double> propagator(size_t(n_state) * size_t(n_state), 0.0);
    vector<double> cross(size_t(n_state) * size_t(n_state), 0.0);
    for(int i = 0; i < n_state; ++i) {
        for(int j = 0; j < n_state; ++j) {
            dense_square_at(propagator, n_state, i, j) = dense_square_at(block_exp, block_dim, i, j);
            dense_square_at(cross, n_state, i, j) = dense_square_at(block_exp, block_dim, i, j + n_state);
        }
    }

    vector<double> covariance = dense_square_multiply(
        cross, dense_square_transpose(propagator, n_state), n_state);
    for(int i = 0; i < n_state; ++i)
        for(int j = i + 1; j < n_state; ++j) {
            double sym = 0.5 * (dense_square_at(covariance, n_state, i, j)
                    + dense_square_at(covariance, n_state, j, i));
            dense_square_at(covariance, n_state, i, j) = sym;
            dense_square_at(covariance, n_state, j, i) = sym;
        }

    if(!dense_square_cholesky_lower(covariance, n_state))
        throw string("cgl_gle exact linear step produced a non-positive covariance");

    runtime->linear_step_matrix.swap(propagator);
    runtime->noise_cholesky_unit.swap(covariance);
    runtime->linear_step_ready = true;
}

static float cgl_gle_transport_scale_at_temperature(
        const vector<float>& temperature_grid,
        const TemperatureScaleTable& scale,
        int mode,
        float temperature) {
    if(temperature_grid.empty() || scale.value.empty()) return 1.f;
    int scale_mode = scale.n_mode == 1 ? 0 : mode;
    if(temperature <= temperature_grid.front()) return scale.value[scale_mode];
    if(temperature >= temperature_grid.back()) return scale.value[(int(temperature_grid.size()) - 1) * scale.n_mode + scale_mode];

    auto upper = std::upper_bound(temperature_grid.begin(), temperature_grid.end(), temperature);
    int upper_idx = int(upper - temperature_grid.begin());
    int lower_idx = upper_idx - 1;
    float lower_t = temperature_grid[lower_idx];
    float upper_t = temperature_grid[upper_idx];
    float alpha = (temperature - lower_t) / (upper_t - lower_t);
    float lower_scale = scale.value[lower_idx * scale.n_mode + scale_mode];
    float upper_scale = scale.value[upper_idx * scale.n_mode + scale_mode];
    return lower_scale + alpha * (upper_scale - lower_scale);
}

static void update_cgl_gle_transport_for_temperature(CGLGLERuntime* runtime) {
    if(!runtime) return;
    for(size_t mode = 0; mode < runtime->coupling.size(); ++mode) {
        float coupling_scale = cgl_gle_transport_scale_at_temperature(
            runtime->temperature_grid, runtime->coupling_scale, int(mode), runtime->temperature);
        float memory_tau_scale = cgl_gle_transport_scale_at_temperature(
            runtime->temperature_grid, runtime->memory_tau_scale, int(mode), runtime->temperature);
        runtime->coupling[mode] = runtime->base_coupling[mode] * coupling_scale;
        runtime->memory_tau[mode] = runtime->base_memory_tau[mode] * memory_tau_scale;
    }
    update_cgl_gle_linear_step(runtime);
}

static float3 project_tangent(float3 v, const float3& n);

template <class Runtime>
static float hidden_state_reference_temperature(const Runtime* runtime) {
    if(!runtime) return -1.f;
    if(runtime->hidden_state_temperature > 0.f) return runtime->hidden_state_temperature;
    if(runtime->temperature_is_set && runtime->temperature > 0.f) return runtime->temperature;
    return -1.f;
}

static inline bool hidden_state_temperature_changed(float old_temperature, float new_temperature) {
    return old_temperature > 0.f && new_temperature > 0.f
        && fabsf(new_temperature - old_temperature) > 1e-6f;
}

static inline float hidden_state_temperature_scale(float old_temperature, float new_temperature) {
    return sqrtf(new_temperature / old_temperature);
}

static void rescale_dynamic_compaction_hidden_state(
        CGLDynamicCompactionRuntime* runtime,
        float old_temperature,
        float new_temperature) {
    if(!runtime || !runtime->state || !hidden_state_temperature_changed(old_temperature, new_temperature)) return;
    float scale = hidden_state_temperature_scale(old_temperature, new_temperature);
    for(int i = 0; i < runtime->state->n_elem; ++i)
        runtime->momentum(0, i) *= scale;
    runtime->hidden_state_temperature = new_temperature;
}

static void rescale_dynamic_orientation_hidden_state(
        CGLDynamicOrientationRuntime* runtime,
        float old_temperature,
        float new_temperature) {
    if(!runtime || !runtime->state || !hidden_state_temperature_changed(old_temperature, new_temperature)) return;
    float scale = hidden_state_temperature_scale(old_temperature, new_temperature);
    for(int i = 0; i < runtime->state->n_elem; ++i) {
        float3 n = normalized(load_vec<3>(runtime->state->output, i));
        float3 l = project_tangent(load_vec<3>(runtime->angular_momentum, i), n) * scale;
        store_vec(runtime->angular_momentum, i, l);
    }
    runtime->hidden_state_temperature = new_temperature;
}

static void rescale_cgl_gle_hidden_state(
        CGLGLERuntime* runtime,
        float old_temperature,
        float new_temperature) {
    if(!runtime || !hidden_state_temperature_changed(old_temperature, new_temperature)) return;
    float scale = hidden_state_temperature_scale(old_temperature, new_temperature);
    int n_cgl = int(runtime->atom_index.size());
    for(int mode = 0; mode < runtime->n_mode(); ++mode) {
        for(int i = 0; i < n_cgl; ++i) {
            int aux_idx = runtime->aux_index(mode, i);
            store_vec(runtime->aux_momentum, aux_idx, load_vec<3>(runtime->aux_momentum, aux_idx) * scale);
        }
    }
    runtime->hidden_state_temperature = new_temperature;
}

static float3 project_tangent(float3 v, const float3& n) {
    return v - n * dot(v, n);
}

static inline float clamp_scalar(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
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

struct CGLipidCompactionSelfPotential : public PotentialNode {
    CoordNode& compaction_state;
    vector<float> self_coeff;
    int n_knot;
    float coord_min;
    float coord_spacing;

    CGLipidCompactionSelfPotential(hid_t grp, CoordNode& compaction_state_);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
};

CGLipidCompactionSelfPotential::CGLipidCompactionSelfPotential(
        hid_t grp, CoordNode& compaction_state_)
    : PotentialNode()
    , compaction_state(compaction_state_)
    , n_knot(read_attribute<int>(grp, ".", "self_n_knot", 0))
    , coord_min(read_attribute<float>(grp, ".", "self_coord_min_ang", 0.f))
    , coord_spacing(read_attribute<float>(grp, ".", "self_coord_spacing_ang", 0.f))
{
    check_elem_width(compaction_state, 1);
    self_coeff = read_float_dataset_1d(grp, "self_coeff");
    if(n_knot <= 3) n_knot = int(self_coeff.size());
    if(int(self_coeff.size()) != n_knot)
        throw string("cg_lipid_compaction_self self_coeff size does not match self_n_knot");
    if(coord_spacing <= 0.f)
        throw string("cg_lipid_compaction_self requires positive self_coord_spacing_ang");
}

void CGLipidCompactionSelfPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray state = compaction_state.output;
    float total = 0.f;
    for(int i = 0; i < compaction_state.n_elem; ++i) {
        float coord = 1.f + (state(0, i) - coord_min) / coord_spacing;
        Vec<2> v = extrapolated_deBoor_value_and_deriv(self_coeff.data(), coord, n_knot);
        total += v.x();
    }
    potential = total;
}

void CGLipidCompactionSelfPotential::propagate_deriv() {
    VecArray state = compaction_state.output;
    VecArray sens = compaction_state.sens;
    for(int i = 0; i < compaction_state.n_elem; ++i) {
        float coord = 1.f + (state(0, i) - coord_min) / coord_spacing;
        Vec<2> v = extrapolated_deBoor_value_and_deriv(self_coeff.data(), coord, n_knot);
        sens(0, i) += v.y() / coord_spacing;
    }
}

struct CGLipidPairPotential : public PotentialNode {
    CoordNode& cg_pos;
    CoordNode* compaction_state;
    vector<float> interaction_param;
    vector<float> delta_extended_extended;
    vector<float> delta_extended_compact;
    vector<float> delta_compact_compact;
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
    bool has_compaction_correction;
    bool implicit_compaction_response;
    bool implicit_compaction_mean_field;
    bool implicit_compaction_gap_response;
    bool implicit_compaction_smooth;
    float implicit_compaction_local_field_scale;
    float implicit_compaction_contact_field_cap;
    bool implicit_compaction_max_contact;
    float compact_state_center;
    float extended_state_center;
    float compact_state_probability;
    float implicit_face_cos_min;
    float implicit_radial_cutoff_ang;
    vector<float> gap_response_coeff;
    int gap_response_n_knot;
    float gap_response_coord_min_ang;
    float gap_response_coord_spacing_ang;
    float gap_response_radial_cutoff_ang;
    float gap_response_face_cos_min;
    float gap_response_fallback_ang;
    bool gap_response_smooth_weight;
    float cache_buffer;
    CGLBodySupport body_support;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_body;
    vector<CGLPairIndex> active_pairs;

    CGLipidPairPotential(hid_t grp, CoordNode& cg_pos);
    CGLipidPairPotential(hid_t grp, const ArgList& arguments);
    void ensure_pairlist(VecArray cg);
    ImplicitCompactionField build_implicit_compaction_field(VecArray cg) const;
    vector<float> solve_mean_field_compaction(VecArray cg) const;
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
    , compaction_state(nullptr)
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
    , has_compaction_correction(false)
    , implicit_compaction_response(read_attribute<int>(grp, ".", "implicit_compaction_response", 0) != 0)
    , implicit_compaction_mean_field(read_attribute<int>(grp, ".", "implicit_compaction_mean_field", 0) != 0)
    , implicit_compaction_gap_response(read_attribute<int>(grp, ".", "implicit_compaction_gap_response", 0) != 0)
    , implicit_compaction_smooth(read_attribute<int>(grp, ".", "implicit_compaction_smooth", 0) != 0)
    , implicit_compaction_local_field_scale(
            read_attribute<float>(grp, ".", "implicit_compaction_local_field_scale", 0.f))
    , implicit_compaction_contact_field_cap(
            read_attribute<float>(grp, ".", "implicit_compaction_contact_field_cap", 0.f))
    , implicit_compaction_max_contact(
            read_attribute<int>(grp, ".", "implicit_compaction_max_contact", 0) != 0)
    , compact_state_center(read_attribute<float>(grp, ".", "compact_state_center_ang", 0.f))
    , extended_state_center(read_attribute<float>(grp, ".", "extended_state_center_ang", 0.f))
    , compact_state_probability(read_attribute<float>(grp, ".", "compact_state_probability", 0.5f))
    , implicit_face_cos_min(read_attribute<float>(grp, ".", "implicit_compaction_face_cos_min", 0.f))
    , implicit_radial_cutoff_ang(read_attribute<float>(grp, ".", "implicit_compaction_radial_cutoff_ang", 0.f))
    , gap_response_n_knot(read_attribute<int>(grp, ".", "gap_response_n_knot", 0))
    , gap_response_coord_min_ang(read_attribute<float>(grp, ".", "gap_response_coord_min_ang", 0.f))
    , gap_response_coord_spacing_ang(read_attribute<float>(grp, ".", "gap_response_coord_spacing_ang", 0.f))
    , gap_response_radial_cutoff_ang(read_attribute<float>(grp, ".", "gap_response_radial_cutoff_ang", 0.f))
    , gap_response_face_cos_min(read_attribute<float>(grp, ".", "gap_response_face_cos_min", 0.f))
    , gap_response_fallback_ang(read_attribute<float>(grp, ".", "gap_response_fallback_ang", 0.f))
    , gap_response_smooth_weight(read_attribute<int>(grp, ".", "gap_response_smooth_weight", 0) != 0)
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
    if(H5Lexists(pi, "delta_extended_extended", H5P_DEFAULT) > 0
            || H5Lexists(pi, "delta_extended_compact", H5P_DEFAULT) > 0
            || H5Lexists(pi, "delta_compact_compact", H5P_DEFAULT) > 0) {
        if(!full_tensor)
            throw string("cg_lipid_pair compaction correction requires a full tensor base table");
        if(!log1p_reduced_transform)
            throw string("cg_lipid_pair compaction correction requires log1p_reduced_transform");
        delta_extended_extended = read_float_dataset_1d(pi, "delta_extended_extended");
        delta_extended_compact = read_float_dataset_1d(pi, "delta_extended_compact");
        delta_compact_compact = read_float_dataset_1d(pi, "delta_compact_compact");
        if(int(delta_extended_extended.size()) != n_param
                || int(delta_extended_compact.size()) != n_param
                || int(delta_compact_compact.size()) != n_param)
            throw string("cg_lipid_pair compaction correction size must match interaction_param");
        if(!(fabsf(compact_state_center - extended_state_center) > 1.0e-6f))
            throw string("cg_lipid_pair requires distinct compact_state_center_ang and extended_state_center_ang");
        has_compaction_correction = true;
    }
    if(H5Lexists(pi, "gap_response_coeff", H5P_DEFAULT) > 0) {
        gap_response_coeff = read_float_dataset_1d(pi, "gap_response_coeff");
        if(gap_response_n_knot <= 3)
            gap_response_n_knot = int(gap_response_coeff.size());
        if(int(gap_response_coeff.size()) != gap_response_n_knot)
            throw string("cg_lipid_pair gap_response_coeff size does not match gap_response_n_knot");
        if(gap_response_coord_spacing_ang <= 0.f || gap_response_radial_cutoff_ang <= 0.f)
            throw string("cg_lipid_pair gap response requires positive coordinate spacing and radial cutoff");
        if(gap_response_fallback_ang <= 0.f)
            gap_response_fallback_ang = gap_response_radial_cutoff_ang;
    } else {
        implicit_compaction_gap_response = false;
    }
}

CGLipidPairPotential::CGLipidPairPotential(hid_t grp, const ArgList& arguments)
    : CGLipidPairPotential(grp, *arguments.at(0))
{
    if(arguments.size() == 2u) {
        compaction_state = arguments.at(1);
        check_elem_width(*compaction_state, 1);
        if(compaction_state->n_elem != int(index.size()))
            throw string("cg_lipid_pair compaction state size must match index size");
        if(!has_compaction_correction)
            throw string("cg_lipid_pair received cgl_compaction_state but lacks compaction correction datasets");
    } else if(arguments.size() != 1u) {
        throw string("cg_lipid_pair expects compose_vector6d or compose_vector6d plus cgl_compaction_state");
    } else if(has_compaction_correction && !implicit_compaction_response && !implicit_compaction_mean_field) {
        throw string("cg_lipid_pair compaction correction datasets require cgl_compaction_state");
    }
}

ImplicitCompactionField CGLipidPairPotential::build_implicit_compaction_field(VecArray cg) const {
    ImplicitCompactionField field;
    field.survival.assign(index.size(), 1.f);
    field.contact_field.assign(index.size(), 0.f);
    field.contact_slope.assign(index.size(), 0.f);
    field.response_slope.assign(index.size(), 0.f);
    field.extended.assign(index.size(), 0.f);
    field.pair_weight.assign(active_pairs.size(), 0.f);
    field.pair_distance.assign(active_pairs.size(), 0.f);
    field.gap_coord.assign(index.size(), gap_response_fallback_ang > 0.f ? gap_response_fallback_ang : 0.f);
    field.gap_weight_sum.assign(index.size(), 0.f);
    if(!implicit_compaction_response || !has_compaction_correction)
        return field;
    if(implicit_compaction_gap_response) {
        return build_gap_response_compaction_field(
                cg,
                index,
                active_pairs,
                box_x,
                box_y,
                box_z,
                gap_response_coeff,
                gap_response_n_knot,
                gap_response_coord_min_ang,
                gap_response_coord_spacing_ang,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_fallback_ang,
                gap_response_smooth_weight);
    }
    field.leaflet_side = classify_leaflets_by_median_z(cg, index);
    field.dominant_pair.assign(index.size(), -1);
    for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
        const CGLPairIndex& pair = active_pairs[pidx];
        size_t ai = pair.first;
        size_t bi = pair.second;
        if(field.leaflet_side[ai] == field.leaflet_side[bi]) continue;
        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
        float w = cgl_cross_leaflet_face_weight(
                dr, n1, n2,
                implicit_radial_cutoff_ang,
                implicit_face_cos_min,
                implicit_compaction_smooth);
        if(w <= 0.f) continue;
        field.pair_weight[pidx] = w;
        if(implicit_compaction_max_contact) {
            float contact = -logf(std::max(1.f - w, 1.0e-6f));
            if(contact > field.contact_field[ai]) {
                field.contact_field[ai] = contact;
                field.dominant_pair[ai] = int(pidx);
            }
            if(contact > field.contact_field[bi]) {
                field.contact_field[bi] = contact;
                field.dominant_pair[bi] = int(pidx);
            }
        } else {
            field.survival[ai] *= std::max(0.f, 1.f - w);
            field.survival[bi] *= std::max(0.f, 1.f - w);
        }
    }
    for(size_t i = 0; i < field.extended.size(); ++i) {
        float raw_contact_field = implicit_compaction_max_contact
            ? field.contact_field[i]
            : -logf(std::max(field.survival[i], 1.0e-6f));
        float contact_field = implicit_compaction_max_contact
            ? ((implicit_compaction_contact_field_cap > 0.f)
                    ? std::min(raw_contact_field, implicit_compaction_contact_field_cap)
                    : raw_contact_field)
            : implicit_compaction_contact_field_from_survival(
                    field.survival[i],
                    implicit_compaction_contact_field_cap);
        field.contact_field[i] = contact_field;
        field.extended[i] = implicit_compaction_probability_from_field(
                contact_field,
                compact_state_probability,
                implicit_compaction_local_field_scale);
        if(implicit_compaction_contact_field_cap > 0.f
                && raw_contact_field >= implicit_compaction_contact_field_cap)
            field.contact_slope[i] = 0.f;
        else
            field.contact_slope[i] = implicit_compaction_probability_field_slope(
                    contact_field,
                    compact_state_probability,
                    implicit_compaction_local_field_scale);
    }
    return field;
}

vector<float> CGLipidPairPotential::solve_mean_field_compaction(VecArray cg) const {
    vector<float> compact(index.size(), clamp_probability(compact_state_probability));
    if(!has_compaction_correction || !implicit_compaction_mean_field)
        return compact;

    vector<MeanFieldPairControl> pair_ctrl(active_pairs.size(), MeanFieldPairControl{0.f, 0.f, 0.f, 0.f, 0.f, false});
    for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
        const CGLPairIndex& pair = active_pairs[pidx];
        size_t ai = pair.first;
        size_t bi = pair.second;
        int t1 = type[ai];
        int t2 = type[bi];

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);

        QuadsplineEval base = {0.f, 0.f, 0.f, 0.f};
        bool ok = full_tensor
            ? eval_full_pair_tensor(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, base)
            : eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                    n_modes, n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, base);
        if(!ok)
            continue;

        QuadsplineEval e_ee = {0.f, 0.f, 0.f, 0.f};
        QuadsplineEval e_ec = {0.f, 0.f, 0.f, 0.f};
        QuadsplineEval e_cc = {0.f, 0.f, 0.f, 0.f};
        bool ok_ee = eval_full_pair_tensor(
                delta_extended_extended.data(), n_angular, n_radial,
                dr, n1, n2, knot_spacing, cutoff,
                log1p_reduced_transform ? 0.f : taper_width, e_ee);
        bool ok_ec = eval_full_pair_tensor(
                delta_extended_compact.data(), n_angular, n_radial,
                dr, n1, n2, knot_spacing, cutoff,
                log1p_reduced_transform ? 0.f : taper_width, e_ec);
        bool ok_cc = eval_full_pair_tensor(
                delta_compact_compact.data(), n_angular, n_radial,
                dr, n1, n2, knot_spacing, cutoff,
                log1p_reduced_transform ? 0.f : taper_width, e_cc);

        MeanFieldPairControl& ctrl = pair_ctrl[pidx];
        ctrl.base_control = base.value;
        ctrl.corr_ee = ok_ee ? e_ee.value : 0.f;
        ctrl.corr_ec = ok_ec ? e_ec.value : 0.f;
        ctrl.corr_cc = ok_cc ? e_cc.value : 0.f;
        ctrl.distance = norm_dr(dr);
        ctrl.valid = ok_ee || ok_ec || ok_cc;
    }

    float logit_p0 = logit(compact_state_probability);
    float inv_kbt = 1.f / std::max(boltzmann_temperature, 1.0e-6f);
    vector<float> field(index.size(), 0.f);
    for(int iter = 0; iter < 24; ++iter) {
        std::fill(field.begin(), field.end(), 0.f);
        for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
            const MeanFieldPairControl& ctrl = pair_ctrl[pidx];
            if(!ctrl.valid)
                continue;
            const CGLPairIndex& pair = active_pairs[pidx];
            size_t ai = pair.first;
            size_t bi = pair.second;
            float si = compact[ai];
            float sj = compact[bi];
            float w_ee = 0.f;
            float w_ec = 0.f;
            float w_cc = 0.f;
            compute_compaction_weights(si, sj, w_ee, w_ec, w_cc);
            float control_value = ctrl.base_control
                + w_ee * ctrl.corr_ee
                + w_ec * ctrl.corr_ec
                + w_cc * ctrl.corr_cc;
            float scale = log1p_reduced_transform
                ? transformed_pair_scale_from_control(
                        control_value,
                        boltzmann_temperature,
                        ctrl.distance,
                        cutoff,
                        taper_width)
                : 1.f;
            float du_dsi = (1.f - sj) * (ctrl.corr_ec - ctrl.corr_ee)
                + sj * (ctrl.corr_cc - ctrl.corr_ec);
            float du_dsj = (1.f - si) * (ctrl.corr_ec - ctrl.corr_ee)
                + si * (ctrl.corr_cc - ctrl.corr_ec);
            field[ai] += scale * du_dsi;
            field[bi] += scale * du_dsj;
        }

        float max_delta = 0.f;
        for(size_t i = 0; i < compact.size(); ++i) {
            float updated = logistic(logit_p0 - inv_kbt * field[i]);
            float blended = 0.5f * compact[i] + 0.5f * updated;
            max_delta = std::max(max_delta, fabsf(blended - compact[i]));
            compact[i] = clamp_probability(blended);
        }
        if(max_delta < 1.0e-5f)
            break;
    }
    return compact;
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
    VecArray comp = compaction_state ? compaction_state->output : VecArray();
    ensure_pairlist(cg);
    ImplicitCompactionField implicit_field = build_implicit_compaction_field(cg);
    vector<float> mean_field_compact = solve_mean_field_compaction(cg);
    float total = 0.f;
    if(has_compaction_correction && implicit_compaction_mean_field && !compaction_state) {
        for(float s : mean_field_compact)
            total += binary_state_relative_entropy(s, compact_state_probability, boltzmann_temperature);
    }
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
            if(has_compaction_correction && (compaction_state || implicit_compaction_response || implicit_compaction_mean_field)) {
                float w_ee = 0.f;
                float w_ec = 0.f;
                float w_cc = 0.f;
                if(compaction_state) {
                    float ci = comp(0, int(ai));
                    float cj = comp(0, int(bi));
                    float denom = compact_state_center - extended_state_center;
                    if(fabsf(denom) < 1.0e-6f)
                        denom = (denom >= 0.f) ? 1.0e-6f : -1.0e-6f;
                    float si = (ci - extended_state_center) / denom;
                    float sj = (cj - extended_state_center) / denom;
                    compute_compaction_weights(si, sj, w_ee, w_ec, w_cc);
                } else if(implicit_compaction_mean_field) {
                    float si = mean_field_compact[ai];
                    float sj = mean_field_compact[bi];
                    compute_compaction_weights(si, sj, w_ee, w_ec, w_cc);
                } else {
                    // Cross-leaflet face contact should raise the compact-state weight.
                    float qi = implicit_field.extended[ai];
                    float qj = implicit_field.extended[bi];
                    compute_compaction_weights(qi, qj, w_ee, w_ec, w_cc);
                }

                QuadsplineEval e_ee = {0.f, 0.f, 0.f, 0.f};
                QuadsplineEval e_ec = {0.f, 0.f, 0.f, 0.f};
                QuadsplineEval e_cc = {0.f, 0.f, 0.f, 0.f};
                bool ok_ee = eval_full_pair_tensor(
                        delta_extended_extended.data(), n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_ee);
                bool ok_ec = eval_full_pair_tensor(
                        delta_extended_compact.data(), n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_ec);
                bool ok_cc = eval_full_pair_tensor(
                        delta_compact_compact.data(), n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_cc);
                if(ok_ee) {
                    e.value += w_ee * e_ee.value;
                    e.d_dr += w_ee * e_ee.d_dr;
                    e.d_da1 += w_ee * e_ee.d_da1;
                    e.d_da2 += w_ee * e_ee.d_da2;
                }
                if(ok_ec) {
                    e.value += w_ec * e_ec.value;
                    e.d_dr += w_ec * e_ec.d_dr;
                    e.d_da1 += w_ec * e_ec.d_da1;
                    e.d_da2 += w_ec * e_ec.d_da2;
                }
                if(ok_cc) {
                    e.value += w_cc * e_cc.value;
                    e.d_dr += w_cc * e_cc.d_dr;
                    e.d_da1 += w_cc * e_cc.d_da1;
                    e.d_da2 += w_cc * e_cc.d_da2;
                }
            }
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
    VecArray comp = compaction_state ? compaction_state->output : VecArray();
    VecArray comp_sens = compaction_state ? compaction_state->sens : VecArray();
    ensure_pairlist(cg);
    ImplicitCompactionField implicit_field = build_implicit_compaction_field(cg);
    vector<float> mean_field_compact = solve_mean_field_compaction(cg);
    vector<float> implicit_q_sens(index.size(), 0.f);

    for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
        const CGLPairIndex& pair = active_pairs[pidx];
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
        float comp_ctrl_sens_ai = 0.f;
        float comp_ctrl_sens_bi = 0.f;
        float implicit_ctrl_sens_ai = 0.f;
        float implicit_ctrl_sens_bi = 0.f;
        if(has_compaction_correction && (compaction_state || implicit_compaction_response || implicit_compaction_mean_field)) {
            float w_ee = 0.f;
            float w_ec = 0.f;
            float w_cc = 0.f;
            float d_w_ee_d_si = 0.f;
            float d_w_ec_d_si = 0.f;
            float d_w_cc_d_si = 0.f;
            float d_w_ee_d_sj = 0.f;
            float d_w_ec_d_sj = 0.f;
            float d_w_cc_d_sj = 0.f;
            float dsi_dci = 0.f;
            float dsj_dcj = 0.f;
            if(compaction_state) {
                float ci = comp(0, int(ai));
                float cj = comp(0, int(bi));
                float denom = compact_state_center - extended_state_center;
                if(fabsf(denom) < 1.0e-6f)
                    denom = (denom >= 0.f) ? 1.0e-6f : -1.0e-6f;
                float si = (ci - extended_state_center) / denom;
                float sj = (cj - extended_state_center) / denom;
                dsi_dci = 1.f / denom;
                dsj_dcj = 1.f / denom;
                compute_compaction_weights(si, sj, w_ee, w_ec, w_cc);
                d_w_ee_d_si = -(1.f - sj);
                d_w_ec_d_si = 1.f - 2.f * sj;
                d_w_cc_d_si = sj;
                d_w_ee_d_sj = -(1.f - si);
                d_w_ec_d_sj = 1.f - 2.f * si;
                d_w_cc_d_sj = si;
            } else if(implicit_compaction_mean_field) {
                float si = mean_field_compact[ai];
                float sj = mean_field_compact[bi];
                compute_compaction_weights(si, sj, w_ee, w_ec, w_cc);
            } else {
                float qi = implicit_field.extended[ai];
                float qj = implicit_field.extended[bi];
                compute_compaction_weights(qi, qj, w_ee, w_ec, w_cc);
            }

            QuadsplineEval e_ee = {0.f, 0.f, 0.f, 0.f};
            QuadsplineEval e_ec = {0.f, 0.f, 0.f, 0.f};
            QuadsplineEval e_cc = {0.f, 0.f, 0.f, 0.f};
            bool ok_ee = eval_full_pair_tensor(
                    delta_extended_extended.data(), n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e_ee);
            bool ok_ec = eval_full_pair_tensor(
                    delta_extended_compact.data(), n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e_ec);
            bool ok_cc = eval_full_pair_tensor(
                    delta_compact_compact.data(), n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e_cc);
            if(ok_ee) {
                e.value += w_ee * e_ee.value;
                e.d_dr += w_ee * e_ee.d_dr;
                e.d_da1 += w_ee * e_ee.d_da1;
                e.d_da2 += w_ee * e_ee.d_da2;
            }
            if(ok_ec) {
                e.value += w_ec * e_ec.value;
                e.d_dr += w_ec * e_ec.d_dr;
                e.d_da1 += w_ec * e_ec.d_da1;
                e.d_da2 += w_ec * e_ec.d_da2;
            }
            if(ok_cc) {
                e.value += w_cc * e_cc.value;
                e.d_dr += w_cc * e_cc.d_dr;
                e.d_da1 += w_cc * e_cc.d_da1;
                e.d_da2 += w_cc * e_cc.d_da2;
            }

            float corr_ee = ok_ee ? e_ee.value : 0.f;
            float corr_ec = ok_ec ? e_ec.value : 0.f;
            float corr_cc = ok_cc ? e_cc.value : 0.f;
            if(compaction_state) {
                comp_ctrl_sens_ai = dsi_dci * (
                        d_w_ee_d_si * corr_ee
                        + d_w_ec_d_si * corr_ec
                        + d_w_cc_d_si * corr_cc);
                comp_ctrl_sens_bi = dsj_dcj * (
                        d_w_ee_d_sj * corr_ee
                        + d_w_ec_d_sj * corr_ec
                        + d_w_cc_d_sj * corr_cc);
            } else if(!implicit_compaction_mean_field) {
                float qi = implicit_field.extended[ai];
                float qj = implicit_field.extended[bi];
                implicit_ctrl_sens_ai = (
                        (qj - 1.f) * corr_ee
                        + (1.f - 2.f * qj) * corr_ec
                        + qj * corr_cc);
                implicit_ctrl_sens_bi = (
                        (qi - 1.f) * corr_ee
                        + (1.f - 2.f * qi) * corr_ec
                        + qi * corr_cc);
            }
        }
        if(log1p_reduced_transform) {
            float taper = 1.f;
            float d_taper_dr = 0.f;
            compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
            float control_value = e.value;
            float control_scale = boltzmann_temperature * expf(control_value) * taper;
            apply_log1p_reduced_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, dr, cutoff, taper_width);
            comp_ctrl_sens_ai *= control_scale;
            comp_ctrl_sens_bi *= control_scale;
            implicit_ctrl_sens_ai *= control_scale;
            implicit_ctrl_sens_bi *= control_scale;
        }

        float dpos1[3] = {0.f, 0.f, 0.f};
        float ddir1[3] = {0.f, 0.f, 0.f};
        float dpos2[3] = {0.f, 0.f, 0.f};
        float ddir2[3] = {0.f, 0.f, 0.f};
        accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
        if(compaction_state && has_compaction_correction) {
            comp_sens(0, int(ai)) += comp_ctrl_sens_ai;
            comp_sens(0, int(bi)) += comp_ctrl_sens_bi;
        } else if(has_compaction_correction && implicit_compaction_response) {
            implicit_q_sens[ai] += implicit_ctrl_sens_ai;
            implicit_q_sens[bi] += implicit_ctrl_sens_bi;
        }
        add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
        add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
    }

    if(has_compaction_correction && implicit_compaction_response && !compaction_state) {
        if(implicit_compaction_gap_response) {
            propagate_gap_response_compaction_sens(
                    cg,
                    cg_sens,
                    index,
                    active_pairs,
                    box_x,
                    box_y,
                    box_z,
                    implicit_field,
                    implicit_q_sens,
                    gap_response_radial_cutoff_ang,
                    gap_response_face_cos_min,
                    gap_response_smooth_weight);
        } else if(implicit_compaction_smooth) {
            for(size_t pidx = 0; pidx < active_pairs.size(); ++pidx) {
                const CGLPairIndex& pair = active_pairs[pidx];
                size_t ai = pair.first;
                size_t bi = pair.second;
                if(implicit_field.leaflet_side.empty()) break;
                if(implicit_field.leaflet_side[ai] == implicit_field.leaflet_side[bi]) continue;
                float w = implicit_field.pair_weight[pidx];
                if(w <= 0.f || w >= 1.f) continue;
                float denom = std::max(1.f - w, 1.0e-6f);
                float lambda = 0.f;
                if(implicit_compaction_max_contact) {
                    if(int(pidx) == implicit_field.dominant_pair[ai])
                        lambda += implicit_q_sens[ai] * (implicit_field.contact_slope[ai] / denom);
                    if(int(pidx) == implicit_field.dominant_pair[bi])
                        lambda += implicit_q_sens[bi] * (implicit_field.contact_slope[bi] / denom);
                } else {
                    lambda = implicit_q_sens[ai] * (implicit_field.contact_slope[ai] / denom)
                           + implicit_q_sens[bi] * (implicit_field.contact_slope[bi] / denom);
                }
                if(fabsf(lambda) <= 1.0e-12f) continue;

                float x1[3], n1[3], x2[3], n2[3], dr[3];
                load_vec6(cg, index[ai], x1, n1);
                load_vec6(cg, index[bi], x2, n2);
                center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
                QuadsplineEval w_eval;
                if(!eval_cgl_cross_leaflet_face_weight(
                            dr, n1, n2,
                            implicit_radial_cutoff_ang,
                            implicit_face_cos_min,
                            implicit_compaction_smooth,
                            w_eval))
                    continue;
                w_eval.d_dr *= lambda;
                w_eval.d_da1 *= lambda;
                w_eval.d_da2 *= lambda;
                float dpos1[3] = {0.f, 0.f, 0.f};
                float ddir1[3] = {0.f, 0.f, 0.f};
                float dpos2[3] = {0.f, 0.f, 0.f};
                float ddir2[3] = {0.f, 0.f, 0.f};
                accumulate_deriv(dr, n1, n2, w_eval, dpos1, ddir1, dpos2, ddir2);
                add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
                add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
            }
        }
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
    vector<uint8_t> leaflet_side = classify_leaflets_by_median_z(cg, index);
    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if(leaflet_side[ai] != leaflet_side[bi]) continue;
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
    vector<uint8_t> leaflet_side = classify_leaflets_by_median_z(cg, index);
    vector<float> dF_drho(index.size(), 0.f);
    for(size_t i = 0; i < index.size(); ++i) {
        float coord = 1.f + (rho[i] - rho_min) / rho_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        dF_drho[i] = emb.y() / rho_spacing;
    }

    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if(leaflet_side[ai] != leaflet_side[bi]) continue;
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

struct CGLipidContactEmbeddingPotential : public PotentialNode {
    CoordNode& cg_pos;
    vector<int> index;
    vector<float> embedding_coeff;
    int n_embedding;
    float box_x;
    float box_y;
    float box_z;
    float coord_min;
    float coord_spacing;
    float radial_cutoff;
    float face_cos_min;
    bool smooth_weight;

    CGLipidContactEmbeddingPotential(hid_t grp, CoordNode& cg_pos);
    void accumulate_contact_field(
            VecArray cg,
            vector<float>& survival,
            vector<float>& contact,
            vector<CGLPairIndex>* pairs = nullptr,
            vector<float>* pair_weight = nullptr) const;
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == embedding_coeff.size()) embedding_coeff = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        radial_cutoff *= scale_xy;
    }
};

CGLipidContactEmbeddingPotential::CGLipidContactEmbeddingPotential(hid_t grp, CoordNode& cg_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , n_embedding(read_attribute<int>(grp, ".", "embedding_n_knot", 0))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , coord_min(read_attribute<float>(grp, ".", "embedding_coord_min", 0.f))
    , coord_spacing(read_attribute<float>(grp, ".", "embedding_coord_spacing", 0.f))
    , radial_cutoff(read_attribute<float>(grp, ".", "contact_radial_cutoff_ang", 0.f))
    , face_cos_min(read_attribute<float>(grp, ".", "contact_face_cos_min", 0.f))
    , smooth_weight(read_attribute<int>(grp, ".", "contact_smooth_weight", 0) != 0)
{
    check_elem_width(cg_pos, 6);
    index = read_int_dataset(grp, "index");
    embedding_coeff = read_float_dataset_1d(grp, "embedding_coeff");
    if(index.empty())
        throw string("cg_lipid_contact_embedding requires at least one CGL index");
    if(n_embedding <= 3)
        n_embedding = int(embedding_coeff.size());
    if(int(embedding_coeff.size()) != n_embedding)
        throw string("cg_lipid_contact_embedding embedding_coeff size does not match embedding_n_knot");
    if(coord_spacing <= 0.f || radial_cutoff <= 0.f)
        throw string("cg_lipid_contact_embedding requires positive embedding spacing and radial cutoff");
}

void CGLipidContactEmbeddingPotential::accumulate_contact_field(
        VecArray cg,
        vector<float>& survival,
        vector<float>& contact,
        vector<CGLPairIndex>* pairs,
        vector<float>* pair_weight) const {
    survival.assign(index.size(), 1.f);
    contact.assign(index.size(), 0.f);
    vector<uint8_t> leaflet_side = classify_leaflets_by_median_z(cg, index);
    if(pairs) pairs->clear();
    if(pair_weight) pair_weight->clear();

    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if(leaflet_side[ai] == leaflet_side[bi]) continue;
            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(cg, index[ai], x1, n1);
            load_vec6(cg, index[bi], x2, n2);
            center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
            float w = cgl_cross_leaflet_face_weight(
                    dr, n1, n2,
                    radial_cutoff,
                    face_cos_min,
                    smooth_weight);
            if(w <= 0.f)
                continue;
            survival[ai] *= std::max(0.f, 1.f - w);
            survival[bi] *= std::max(0.f, 1.f - w);
            if(pairs) pairs->push_back({ai, bi});
            if(pair_weight) pair_weight->push_back(w);
        }
    }
    for(size_t i = 0; i < contact.size(); ++i)
        contact[i] = 1.f - survival[i];
}

void CGLipidContactEmbeddingPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray cg = cg_pos.output;
    vector<float> survival;
    vector<float> contact;
    accumulate_contact_field(cg, survival, contact);
    float total = 0.f;
    for(float q : contact) {
        float coord = 1.f + (q - coord_min) / coord_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        total += emb.x();
    }
    potential = total;
}

void CGLipidContactEmbeddingPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    vector<float> survival;
    vector<float> contact;
    vector<CGLPairIndex> pairs;
    vector<float> pair_weight;
    accumulate_contact_field(cg, survival, contact, &pairs, &pair_weight);
    vector<float> dF_dq(index.size(), 0.f);
    for(size_t i = 0; i < index.size(); ++i) {
        float coord = 1.f + (contact[i] - coord_min) / coord_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        dF_dq[i] = emb.y() / coord_spacing;
    }

    if(!smooth_weight)
        return;

    for(size_t pidx = 0; pidx < pairs.size(); ++pidx) {
        size_t ai = pairs[pidx].first;
        size_t bi = pairs[pidx].second;
        float w = pair_weight[pidx];
        if(w <= 0.f || w >= 1.f)
            continue;
        float denom = std::max(1.f - w, 1.0e-6f);
        float lambda = dF_dq[ai] * (survival[ai] / denom)
                     + dF_dq[bi] * (survival[bi] / denom);
        if(fabsf(lambda) <= 1.0e-12f)
            continue;

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
        QuadsplineEval w_eval;
        if(!eval_cgl_cross_leaflet_face_weight(
                    dr, n1, n2,
                    radial_cutoff,
                    face_cos_min,
                    smooth_weight,
                    w_eval))
            continue;
        w_eval.d_dr *= lambda;
        w_eval.d_da1 *= lambda;
        w_eval.d_da2 *= lambda;
        float dpos1[3] = {0.f, 0.f, 0.f};
        float ddir1[3] = {0.f, 0.f, 0.f};
        float dpos2[3] = {0.f, 0.f, 0.f};
        float ddir2[3] = {0.f, 0.f, 0.f};
        accumulate_deriv(dr, n1, n2, w_eval, dpos1, ddir1, dpos2, ddir2);
        add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
        add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
    }
}

struct CGLipidGapEmbeddingPotential : public PotentialNode {
    CoordNode& cg_pos;
    vector<int> index;
    vector<float> embedding_coeff;
    int n_embedding;
    float box_x;
    float box_y;
    float box_z;
    float coord_min;
    float coord_spacing;
    float radial_cutoff;
    float face_cos_min;
    float fallback_gap;
    bool smooth_weight;

    CGLipidGapEmbeddingPotential(hid_t grp, CoordNode& cg_pos);
    void accumulate_gap_field(
            VecArray cg,
            vector<float>& gap,
            vector<float>& gap_weight_sum,
            vector<CGLPairIndex>* pairs = nullptr,
            vector<float>* pair_weight = nullptr,
            vector<float>* pair_distance = nullptr) const;
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == embedding_coeff.size()) embedding_coeff = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
        radial_cutoff *= scale_xy;
        fallback_gap *= scale_xy;
    }
};

CGLipidGapEmbeddingPotential::CGLipidGapEmbeddingPotential(hid_t grp, CoordNode& cg_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , n_embedding(read_attribute<int>(grp, ".", "embedding_n_knot", 0))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , coord_min(read_attribute<float>(grp, ".", "embedding_coord_min", 0.f))
    , coord_spacing(read_attribute<float>(grp, ".", "embedding_coord_spacing", 0.f))
    , radial_cutoff(read_attribute<float>(grp, ".", "gap_radial_cutoff_ang", 0.f))
    , face_cos_min(read_attribute<float>(grp, ".", "gap_face_cos_min", 0.f))
    , fallback_gap(read_attribute<float>(grp, ".", "gap_fallback_ang", 0.f))
    , smooth_weight(read_attribute<int>(grp, ".", "gap_smooth_weight", 0) != 0)
{
    check_elem_width(cg_pos, 6);
    index = read_int_dataset(grp, "index");
    embedding_coeff = read_float_dataset_1d(grp, "embedding_coeff");
    if(index.empty())
        throw string("cg_lipid_gap_embedding requires at least one CGL index");
    if(n_embedding <= 3)
        n_embedding = int(embedding_coeff.size());
    if(int(embedding_coeff.size()) != n_embedding)
        throw string("cg_lipid_gap_embedding embedding_coeff size does not match embedding_n_knot");
    if(coord_spacing <= 0.f || radial_cutoff <= 0.f)
        throw string("cg_lipid_gap_embedding requires positive embedding spacing and radial cutoff");
    if(fallback_gap <= 0.f)
        fallback_gap = radial_cutoff;
}

void CGLipidGapEmbeddingPotential::accumulate_gap_field(
        VecArray cg,
        vector<float>& gap,
        vector<float>& gap_weight_sum,
        vector<CGLPairIndex>* pairs,
        vector<float>* pair_weight,
        vector<float>* pair_distance) const {
    gap.assign(index.size(), fallback_gap);
    gap_weight_sum.assign(index.size(), 0.f);
    vector<float> gap_weighted_sum(index.size(), 0.f);
    vector<uint8_t> leaflet_side = classify_leaflets_by_median_z(cg, index);
    if(pairs) pairs->clear();
    if(pair_weight) pair_weight->clear();
    if(pair_distance) pair_distance->clear();

    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if(leaflet_side[ai] == leaflet_side[bi]) continue;
            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(cg, index[ai], x1, n1);
            load_vec6(cg, index[bi], x2, n2);
            center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);
            float w = cgl_cross_leaflet_face_weight(
                    dr, n1, n2,
                    radial_cutoff,
                    face_cos_min,
                    smooth_weight);
            if(w <= 0.f)
                continue;
            float r = norm_dr(dr);
            gap_weighted_sum[ai] += w * r;
            gap_weighted_sum[bi] += w * r;
            gap_weight_sum[ai] += w;
            gap_weight_sum[bi] += w;
            if(pairs) pairs->push_back({ai, bi});
            if(pair_weight) pair_weight->push_back(w);
            if(pair_distance) pair_distance->push_back(r);
        }
    }
    for(size_t i = 0; i < gap.size(); ++i) {
        if(gap_weight_sum[i] > 1.0e-12f)
            gap[i] = gap_weighted_sum[i] / gap_weight_sum[i];
    }
}

void CGLipidGapEmbeddingPotential::compute_value(ComputeMode mode) {
    (void)mode;
    if(mode == DerivMode) {
        potential = 0.f;
        return;
    }
    VecArray cg = cg_pos.output;
    vector<float> gap;
    vector<float> gap_weight_sum;
    accumulate_gap_field(cg, gap, gap_weight_sum);
    float total = 0.f;
    for(float g : gap) {
        float coord = 1.f + (g - coord_min) / coord_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        total += emb.x();
    }
    potential = total;
}

void CGLipidGapEmbeddingPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;
    vector<float> gap;
    vector<float> gap_weight_sum;
    vector<CGLPairIndex> pairs;
    vector<float> pair_weight;
    vector<float> pair_distance;
    accumulate_gap_field(cg, gap, gap_weight_sum, &pairs, &pair_weight, &pair_distance);
    vector<float> dF_dg(index.size(), 0.f);
    for(size_t i = 0; i < index.size(); ++i) {
        float coord = 1.f + (gap[i] - coord_min) / coord_spacing;
        Vec<2> emb = clamped_deBoor_value_and_deriv(
                embedding_coeff.data(), coord, n_embedding);
        dF_dg[i] = emb.y() / coord_spacing;
    }

    for(size_t pidx = 0; pidx < pairs.size(); ++pidx) {
        size_t ai = pairs[pidx].first;
        size_t bi = pairs[pidx].second;
        float w = pair_weight[pidx];
        if(w <= 0.f)
            continue;
        float r = pair_distance[pidx];
        float lambda_r = 0.f;
        float lambda_w = 0.f;
        if(gap_weight_sum[ai] > 1.0e-12f) {
            float inv_sum = 1.f / gap_weight_sum[ai];
            lambda_r += dF_dg[ai] * (w * inv_sum);
            lambda_w += dF_dg[ai] * ((r - gap[ai]) * inv_sum);
        }
        if(gap_weight_sum[bi] > 1.0e-12f) {
            float inv_sum = 1.f / gap_weight_sum[bi];
            lambda_r += dF_dg[bi] * (w * inv_sum);
            lambda_w += dF_dg[bi] * ((r - gap[bi]) * inv_sum);
        }
        if(fabsf(lambda_r) <= 1.0e-12f && (!smooth_weight || fabsf(lambda_w) <= 1.0e-12f))
            continue;

        float x1[3], n1[3], x2[3], n2[3], dr[3];
        load_vec6(cg, index[ai], x1, n1);
        load_vec6(cg, index[bi], x2, n2);
        center_dr(cg, index[ai], cg, index[bi], box_x, box_y, box_z, dr);

        if(fabsf(lambda_r) > 1.0e-12f) {
            QuadsplineEval r_eval = {0.f, lambda_r, 0.f, 0.f};
            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, r_eval, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
        }

        if(smooth_weight && fabsf(lambda_w) > 1.0e-12f) {
            QuadsplineEval w_eval;
            if(!eval_cgl_cross_leaflet_face_weight(
                        dr, n1, n2,
                        radial_cutoff,
                        face_cos_min,
                        smooth_weight,
                        w_eval))
                continue;
            w_eval.d_dr *= lambda_w;
            w_eval.d_da1 *= lambda_w;
            w_eval.d_da2 *= lambda_w;
            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, w_eval, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
        }
    }
}

struct CGLipidSCPotential : public PotentialNode {
    CoordNode& sc_pos;
    CoordNode& cg_pos;
    vector<float> interaction_param;
    vector<float> delta_extended;
    vector<float> delta_compact;
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
    bool has_compaction_correction;
    bool implicit_compaction_gap_response;
    int gap_response_n_knot;
    float gap_response_coord_min_ang;
    float gap_response_coord_spacing_ang;
    float gap_response_radial_cutoff_ang;
    float gap_response_face_cos_min;
    float gap_response_fallback_ang;
    bool gap_response_smooth_weight;
    float cache_buffer;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_pos1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;
    vector<CGLPairIndex> active_cg_pairs;
    vector<float> gap_response_coeff;

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
    , has_compaction_correction(false)
    , implicit_compaction_gap_response(read_attribute<int>(grp, ".", "implicit_compaction_gap_response", 0) != 0)
    , gap_response_n_knot(read_attribute<int>(grp, ".", "gap_response_n_knot", 0))
    , gap_response_coord_min_ang(read_attribute<float>(grp, ".", "gap_response_coord_min_ang", 0.f))
    , gap_response_coord_spacing_ang(read_attribute<float>(grp, ".", "gap_response_coord_spacing_ang", 0.f))
    , gap_response_radial_cutoff_ang(read_attribute<float>(grp, ".", "gap_response_radial_cutoff_ang", 0.f))
    , gap_response_face_cos_min(read_attribute<float>(grp, ".", "gap_response_face_cos_min", 0.f))
    , gap_response_fallback_ang(read_attribute<float>(grp, ".", "gap_response_fallback_ang", 0.f))
    , gap_response_smooth_weight(read_attribute<int>(grp, ".", "gap_response_smooth_weight", 0) != 0)
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
    if(H5Lexists(pi, "delta_extended", H5P_DEFAULT) > 0
            || H5Lexists(pi, "delta_compact", H5P_DEFAULT) > 0) {
        if(!log1p_reduced_transform)
            throw string("cg_lipid_sc compaction correction requires log1p_reduced_transform");
        delta_extended = read_named_param_dataset_any(pi, "delta_extended", n_type1, n_type2, n_param);
        delta_compact = read_named_param_dataset_any(pi, "delta_compact", n_type1, n_type2, n_param);
        has_compaction_correction = true;
        if(H5Lexists(pi, "gap_response_coeff", H5P_DEFAULT) > 0) {
            gap_response_coeff = read_float_dataset_1d(pi, "gap_response_coeff");
            if(gap_response_n_knot <= 3)
                gap_response_n_knot = int(gap_response_coeff.size());
            if(int(gap_response_coeff.size()) != gap_response_n_knot)
                throw string("cg_lipid_sc gap_response_coeff size does not match gap_response_n_knot");
        } else {
            implicit_compaction_gap_response = false;
        }
        if(!implicit_compaction_gap_response)
            throw string("cg_lipid_sc compaction correction currently requires implicit gap response metadata");
        if(gap_response_coord_spacing_ang <= 0.f || gap_response_radial_cutoff_ang <= 0.f)
            throw string("cg_lipid_sc gap response requires positive coordinate spacing and radial cutoff");
        if(gap_response_fallback_ang <= 0.f)
            gap_response_fallback_ang = gap_response_radial_cutoff_ang;
    }
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
    active_cg_pairs.clear();
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
    if(has_compaction_correction) {
        float cg_cutoff = std::max(gap_response_radial_cutoff_ang, 0.f) + cache_buffer;
        float cg_cutoff2 = cg_cutoff * cg_cutoff;
        for(size_t ai = 0; ai < index2.size(); ++ai) {
            for(size_t bi = ai + 1; bi < index2.size(); ++bi) {
                float dr[3];
                center_dr(cg, index2[ai], cg, index2[bi], box_x, box_y, box_z, dr);
                float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if(r2 < cg_cutoff2)
                    active_cg_pairs.push_back({ai, bi});
            }
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
    ImplicitCompactionField implicit_field;
    if(has_compaction_correction) {
        implicit_field = build_gap_response_compaction_field(
                cg,
                index2,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                gap_response_coeff,
                gap_response_n_knot,
                gap_response_coord_min_ang,
                gap_response_coord_spacing_ang,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_fallback_ang,
                gap_response_smooth_weight);
    }

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
            if(has_compaction_correction) {
                QuadsplineEval e_ext = {0.f, 0.f, 0.f, 0.f};
                QuadsplineEval e_comp = {0.f, 0.f, 0.f, 0.f};
                bool ok_ext = eval_full_pair_tensor(
                        param_ptr(delta_extended, n_type2, n_param, t1, t2),
                        n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_ext);
                bool ok_comp = eval_full_pair_tensor(
                        param_ptr(delta_compact, n_type2, n_param, t1, t2),
                        n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_comp);
                float q = implicit_field.extended[bi];
                if(ok_ext) {
                    float w = 1.f - q;
                    e.value += w * e_ext.value;
                    e.d_dr += w * e_ext.d_dr;
                    e.d_da1 += w * e_ext.d_da1;
                    e.d_da2 += w * e_ext.d_da2;
                }
                if(ok_comp) {
                    e.value += q * e_comp.value;
                    e.d_dr += q * e_comp.d_dr;
                    e.d_da1 += q * e_comp.d_da1;
                    e.d_da2 += q * e_comp.d_da2;
                }
            }
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
    ImplicitCompactionField implicit_field;
    vector<float> implicit_q_sens(index2.size(), 0.f);
    if(has_compaction_correction) {
        implicit_field = build_gap_response_compaction_field(
                cg,
                index2,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                gap_response_coeff,
                gap_response_n_knot,
                gap_response_coord_min_ang,
                gap_response_coord_spacing_ang,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_fallback_ang,
                gap_response_smooth_weight);
    }

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
        float q_ctrl_sens = 0.f;
        if(has_compaction_correction) {
            QuadsplineEval e_ext = {0.f, 0.f, 0.f, 0.f};
            QuadsplineEval e_comp = {0.f, 0.f, 0.f, 0.f};
            bool ok_ext = eval_full_pair_tensor(
                    param_ptr(delta_extended, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e_ext);
            bool ok_comp = eval_full_pair_tensor(
                    param_ptr(delta_compact, n_type2, n_param, t1, t2),
                    n_angular, n_radial,
                    dr, n1, n2, knot_spacing, cutoff,
                    log1p_reduced_transform ? 0.f : taper_width, e_comp);
            float q = implicit_field.extended[bi];
            float corr_ext = ok_ext ? e_ext.value : 0.f;
            float corr_comp = ok_comp ? e_comp.value : 0.f;
            if(ok_ext) {
                float w = 1.f - q;
                e.value += w * e_ext.value;
                e.d_dr += w * e_ext.d_dr;
                e.d_da1 += w * e_ext.d_da1;
                e.d_da2 += w * e_ext.d_da2;
            }
            if(ok_comp) {
                e.value += q * e_comp.value;
                e.d_dr += q * e_comp.d_dr;
                e.d_da1 += q * e_comp.d_da1;
                e.d_da2 += q * e_comp.d_da2;
            }
            q_ctrl_sens = corr_comp - corr_ext;
        }
        if(log1p_reduced_transform)
        {
            float taper = 1.f;
            float d_taper_dr = 0.f;
            compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
            float control_scale = boltzmann_temperature * expf(e.value) * taper;
            q_ctrl_sens *= control_scale;
            apply_log1p_reduced_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, dr, cutoff, taper_width);
        }

        float dpos1[3] = {0.f, 0.f, 0.f};
        float ddir1[3] = {0.f, 0.f, 0.f};
        float dpos2[3] = {0.f, 0.f, 0.f};
        float ddir2[3] = {0.f, 0.f, 0.f};
        accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
        add_vec6_sens(sc_sens, index1[ai], dpos1, ddir1);
        add_vec6_sens(cg_sens, index2[bi], dpos2, ddir2);
        if(has_compaction_correction)
            implicit_q_sens[bi] += q_ctrl_sens;
    }
    if(has_compaction_correction)
        propagate_gap_response_compaction_sens(
                cg,
                cg_sens,
                index2,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                implicit_field,
                implicit_q_sens,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_smooth_weight);
}

struct CGLipidTargetPotential : public PotentialNode {
    CoordNode& cg_pos;
    CoordNode& tgt_pos;
    CoordNode* compaction_state;
    vector<float> interaction_param;
    vector<float> delta_extended;
    vector<float> delta_compact;
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
    bool has_compaction_correction;
    float compact_state_center;
    float extended_state_center;
    bool implicit_compaction_gap_response;
    int gap_response_n_knot;
    float gap_response_coord_min_ang;
    float gap_response_coord_spacing_ang;
    float gap_response_radial_cutoff_ang;
    float gap_response_face_cos_min;
    float gap_response_fallback_ang;
    bool gap_response_smooth_weight;
    float cache_buffer;
    CGLBodySupport body_support;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_body1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;
    vector<CGLPairIndex> active_cg_pairs;
    vector<float> gap_response_coeff;

    CGLipidTargetPotential(hid_t grp, CoordNode& cg_pos_, CoordNode& tgt_pos_);
    CGLipidTargetPotential(hid_t grp, const ArgList& arguments);
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
    CoordNode* compaction_state;
    vector<float> interaction_param;
    vector<float> delta_extended;
    vector<float> delta_compact;
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
    bool has_compaction_correction;
    bool implicit_compaction_gap_response;
    int gap_response_n_knot;
    float gap_response_coord_min_ang;
    float gap_response_coord_spacing_ang;
    float gap_response_radial_cutoff_ang;
    float gap_response_face_cos_min;
    float gap_response_fallback_ang;
    bool gap_response_smooth_weight;
    float cache_buffer;
    float cached_box_x;
    float cached_box_y;
    float cached_box_z;
    bool pairlist_valid;
    vector<float> cached_pos1;
    vector<float> cached_pos2;
    vector<CGLPairIndex> active_pairs;
    vector<CGLPairIndex> active_cg_pairs;
    vector<float> gap_response_coeff;

    CGLipidSCOneBody(hid_t grp, CoordNode& sc_pos_, CoordNode& cg_pos_)
        : CoordNode(sc_pos_.n_elem, 1)
        , sc_pos(sc_pos_)
        , cg_pos(cg_pos_)
        , compaction_state(nullptr)
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
        , has_compaction_correction(false)
        , implicit_compaction_gap_response(read_attribute<int>(grp, ".", "implicit_compaction_gap_response", 0) != 0)
        , gap_response_n_knot(read_attribute<int>(grp, ".", "gap_response_n_knot", 0))
        , gap_response_coord_min_ang(read_attribute<float>(grp, ".", "gap_response_coord_min_ang", 0.f))
        , gap_response_coord_spacing_ang(read_attribute<float>(grp, ".", "gap_response_coord_spacing_ang", 0.f))
        , gap_response_radial_cutoff_ang(read_attribute<float>(grp, ".", "gap_response_radial_cutoff_ang", 0.f))
        , gap_response_face_cos_min(read_attribute<float>(grp, ".", "gap_response_face_cos_min", 0.f))
        , gap_response_fallback_ang(read_attribute<float>(grp, ".", "gap_response_fallback_ang", 0.f))
        , gap_response_smooth_weight(read_attribute<int>(grp, ".", "gap_response_smooth_weight", 0) != 0)
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
        if(H5Lexists(pi, "delta_extended", H5P_DEFAULT) > 0
                || H5Lexists(pi, "delta_compact", H5P_DEFAULT) > 0) {
            if(!log1p_reduced_transform)
                throw string("cg_lipid_rotamer_sc compaction correction requires log1p_reduced_transform");
            delta_extended = read_named_param_dataset_any(pi, "delta_extended", n_type1, n_type2, n_param);
            delta_compact = read_named_param_dataset_any(pi, "delta_compact", n_type1, n_type2, n_param);
            has_compaction_correction = true;
            if(H5Lexists(pi, "gap_response_coeff", H5P_DEFAULT) > 0) {
                gap_response_coeff = read_float_dataset_1d(pi, "gap_response_coeff");
                if(gap_response_n_knot <= 3)
                    gap_response_n_knot = int(gap_response_coeff.size());
                if(int(gap_response_coeff.size()) != gap_response_n_knot)
                    throw string("cg_lipid_rotamer_sc gap_response_coeff size does not match gap_response_n_knot");
            } else {
                implicit_compaction_gap_response = false;
            }
            if(implicit_compaction_gap_response
                    && (gap_response_coord_spacing_ang <= 0.f || gap_response_radial_cutoff_ang <= 0.f))
                throw string("cg_lipid_rotamer_sc gap response requires positive coordinate spacing and radial cutoff");
            if(implicit_compaction_gap_response && gap_response_fallback_ang <= 0.f)
                gap_response_fallback_ang = gap_response_radial_cutoff_ang;
        }

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

    CGLipidSCOneBody(hid_t grp, const ArgList& arguments)
        : CGLipidSCOneBody(grp, *arguments.at(0), *arguments.at(1))
    {
        if(arguments.size() == 3u) {
            compaction_state = arguments.at(2);
            check_elem_width(*compaction_state, 1);
            if(compaction_state->n_elem != int(index2.size()))
                throw string("cg_lipid_rotamer_sc compaction state size must match CGL count");
            if(!has_compaction_correction)
                throw string("cg_lipid_rotamer_sc received cgl_compaction_state but lacks compaction correction datasets");
        } else if(arguments.size() != 2u) {
            throw string(
                    "cg_lipid_rotamer_sc expects placement_fixed_point_vector_only plus compose_vector6d, "
                    "optionally with cgl_compaction_state");
        } else if(has_compaction_correction && !implicit_compaction_gap_response) {
            throw string("cg_lipid_rotamer_sc compaction correction currently requires implicit gap response metadata");
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
        active_cg_pairs.clear();
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
        if(has_compaction_correction && !compaction_state) {
            float cg_cutoff = std::max(gap_response_radial_cutoff_ang, 0.f) + cache_buffer;
            float cg_cutoff2 = cg_cutoff * cg_cutoff;
            for(size_t ai = 0; ai < index2.size(); ++ai) {
                for(size_t bi = ai + 1; bi < index2.size(); ++bi) {
                    float dr[3];
                    center_dr(cg, index2[ai], cg, index2[bi], box_x, box_y, box_z, dr);
                    float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                    if(r2 < cg_cutoff2)
                        active_cg_pairs.push_back({ai, bi});
                }
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
        VecArray comp = compaction_state ? compaction_state->output : VecArray();
        ensure_pairlist(sc, cg);
        ImplicitCompactionField implicit_field;
        if(has_compaction_correction && !compaction_state) {
            implicit_field = build_gap_response_compaction_field(
                    cg,
                    index2,
                    active_cg_pairs,
                    box_x,
                    box_y,
                    box_z,
                    gap_response_coeff,
                    gap_response_n_knot,
                    gap_response_coord_min_ang,
                    gap_response_coord_spacing_ang,
                    gap_response_radial_cutoff_ang,
                    gap_response_face_cos_min,
                    gap_response_fallback_ang,
                    gap_response_smooth_weight);
        }

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
                if(has_compaction_correction) {
                    QuadsplineEval e_ext = {0.f, 0.f, 0.f, 0.f};
                    QuadsplineEval e_comp = {0.f, 0.f, 0.f, 0.f};
                    bool ok_ext = eval_full_pair_tensor(
                            param_ptr(delta_extended, n_type2, n_param, t1, t2),
                            n_angular, n_radial,
                            dr, n1, n2, knot_spacing, cutoff,
                            log1p_reduced_transform ? 0.f : taper_width, e_ext);
                    bool ok_comp = eval_full_pair_tensor(
                            param_ptr(delta_compact, n_type2, n_param, t1, t2),
                            n_angular, n_radial,
                            dr, n1, n2, knot_spacing, cutoff,
                            log1p_reduced_transform ? 0.f : taper_width, e_comp);
                    float q = compaction_state
                        ? comp(0, index2[bi])
                        : implicit_field.extended[bi];
                    if(ok_ext) {
                        float w = 1.f - q;
                        e.value += w * e_ext.value;
                        e.d_dr += w * e_ext.d_dr;
                        e.d_da1 += w * e_ext.d_da1;
                        e.d_da2 += w * e_ext.d_da2;
                    }
                    if(ok_comp) {
                        e.value += q * e_comp.value;
                        e.d_dr += q * e_comp.d_dr;
                        e.d_da1 += q * e_comp.d_da1;
                        e.d_da2 += q * e_comp.d_da2;
                    }
                }
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
        VecArray comp = compaction_state ? compaction_state->output : VecArray();
        VecArray comp_sens = compaction_state ? compaction_state->sens : VecArray();
        ensure_pairlist(sc, cg);
        ImplicitCompactionField implicit_field;
        vector<float> implicit_q_sens(index2.size(), 0.f);
        if(has_compaction_correction && !compaction_state) {
            implicit_field = build_gap_response_compaction_field(
                    cg,
                    index2,
                    active_cg_pairs,
                    box_x,
                    box_y,
                    box_z,
                    gap_response_coeff,
                    gap_response_n_knot,
                    gap_response_coord_min_ang,
                    gap_response_coord_spacing_ang,
                    gap_response_radial_cutoff_ang,
                    gap_response_face_cos_min,
                    gap_response_fallback_ang,
                    gap_response_smooth_weight);
        }

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
            float q_ctrl_sens = 0.f;
            if(has_compaction_correction) {
                QuadsplineEval e_ext = {0.f, 0.f, 0.f, 0.f};
                QuadsplineEval e_comp = {0.f, 0.f, 0.f, 0.f};
                bool ok_ext = eval_full_pair_tensor(
                        param_ptr(delta_extended, n_type2, n_param, t1, t2),
                        n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_ext);
                bool ok_comp = eval_full_pair_tensor(
                        param_ptr(delta_compact, n_type2, n_param, t1, t2),
                        n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff,
                        log1p_reduced_transform ? 0.f : taper_width, e_comp);
                float q = compaction_state
                    ? comp(0, index2[bi])
                    : implicit_field.extended[bi];
                float corr_ext = ok_ext ? e_ext.value : 0.f;
                float corr_comp = ok_comp ? e_comp.value : 0.f;
                if(ok_ext) {
                    float w = 1.f - q;
                    e.value += w * e_ext.value;
                    e.d_dr += w * e_ext.d_dr;
                    e.d_da1 += w * e_ext.d_da1;
                    e.d_da2 += w * e_ext.d_da2;
                }
                if(ok_comp) {
                    e.value += q * e_comp.value;
                    e.d_dr += q * e_comp.d_dr;
                    e.d_da1 += q * e_comp.d_da1;
                    e.d_da2 += q * e_comp.d_da2;
                }
                q_ctrl_sens = corr_comp - corr_ext;
            }
            if(log1p_reduced_transform)
            {
                float taper = 1.f;
                float d_taper_dr = 0.f;
                compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
                float control_scale = boltzmann_temperature * expf(e.value) * taper;
                q_ctrl_sens *= control_scale;
                apply_log1p_reduced_transform(
                        e, reference_energy_eup[t1 * n_type2 + t2],
                        boltzmann_temperature, dr, cutoff, taper_width);
            }

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
            if(has_compaction_correction && compaction_state)
                comp_sens(0, index2[bi]) += row_scale * q_ctrl_sens;
            else if(has_compaction_correction)
                implicit_q_sens[bi] += row_scale * q_ctrl_sens;
        }
        if(has_compaction_correction && !compaction_state)
            propagate_gap_response_compaction_sens(
                    cg,
                    cg_sens,
                    index2,
                    active_cg_pairs,
                    box_x,
                    box_y,
                    box_z,
                    implicit_field,
                    implicit_q_sens,
                    gap_response_radial_cutoff_ang,
                    gap_response_face_cos_min,
                    gap_response_smooth_weight);
    }
};

CGLipidTargetPotential::CGLipidTargetPotential(
        hid_t grp, CoordNode& cg_pos_, CoordNode& tgt_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , tgt_pos(tgt_pos_)
    , compaction_state(nullptr)
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
        , has_compaction_correction(false)
        , compact_state_center(read_attribute<float>(grp, ".", "compact_state_center_ang", 0.f))
        , extended_state_center(read_attribute<float>(grp, ".", "extended_state_center_ang", 0.f))
        , implicit_compaction_gap_response(read_attribute<int>(grp, ".", "implicit_compaction_gap_response", 0) != 0)
        , gap_response_n_knot(read_attribute<int>(grp, ".", "gap_response_n_knot", 0))
        , gap_response_coord_min_ang(read_attribute<float>(grp, ".", "gap_response_coord_min_ang", 0.f))
        , gap_response_coord_spacing_ang(read_attribute<float>(grp, ".", "gap_response_coord_spacing_ang", 0.f))
        , gap_response_radial_cutoff_ang(read_attribute<float>(grp, ".", "gap_response_radial_cutoff_ang", 0.f))
        , gap_response_face_cos_min(read_attribute<float>(grp, ".", "gap_response_face_cos_min", 0.f))
        , gap_response_fallback_ang(read_attribute<float>(grp, ".", "gap_response_fallback_ang", 0.f))
        , gap_response_smooth_weight(read_attribute<int>(grp, ".", "gap_response_smooth_weight", 0) != 0)
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
    if(H5Lexists(pi, "delta_extended", H5P_DEFAULT) > 0
            || H5Lexists(pi, "delta_compact", H5P_DEFAULT) > 0) {
        if(!log1p_reduced_transform)
            throw string("cg_lipid_target compaction correction requires log1p_reduced_transform");
        delta_extended = read_named_param_dataset_any(pi, "delta_extended", n_type1, n_type2, n_param);
        delta_compact = read_named_param_dataset_any(pi, "delta_compact", n_type1, n_type2, n_param);
        has_compaction_correction = true;
        if(H5Lexists(pi, "gap_response_coeff", H5P_DEFAULT) > 0) {
            gap_response_coeff = read_float_dataset_1d(pi, "gap_response_coeff");
            if(gap_response_n_knot <= 3)
                gap_response_n_knot = int(gap_response_coeff.size());
            if(int(gap_response_coeff.size()) != gap_response_n_knot)
                throw string("cg_lipid_target gap_response_coeff size does not match gap_response_n_knot");
        } else {
            implicit_compaction_gap_response = false;
        }
        if(gap_response_coord_spacing_ang <= 0.f || gap_response_radial_cutoff_ang <= 0.f)
            implicit_compaction_gap_response = false;
        if(gap_response_fallback_ang <= 0.f)
            gap_response_fallback_ang = gap_response_radial_cutoff_ang;
    }
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

CGLipidTargetPotential::CGLipidTargetPotential(hid_t grp, const ArgList& arguments)
    : CGLipidTargetPotential(grp, *arguments.at(0), *arguments.at(1))
{
    if(arguments.size() == 3u) {
        compaction_state = arguments.at(2);
        check_elem_width(*compaction_state, 1);
        if(compaction_state->n_elem != int(index1.size()))
            throw string("cg_lipid_target compaction state size must match CGL count");
        if(!has_compaction_correction)
            throw string("cg_lipid_target received cgl_compaction_state but lacks compaction correction datasets");
    } else if(arguments.size() != 2u) {
        throw string("cg_lipid_target expects compose_vector6d plus pos, optionally with cgl_compaction_state");
    } else if(has_compaction_correction && !implicit_compaction_gap_response) {
        throw string("cg_lipid_target compaction correction currently requires implicit gap response metadata");
    }
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
    active_cg_pairs.clear();
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
    if(has_compaction_correction) {
        float cg_cutoff = std::max(gap_response_radial_cutoff_ang, 0.f) + cache_buffer;
        float cg_cutoff2 = cg_cutoff * cg_cutoff;
        for(size_t ai = 0; ai < index1.size(); ++ai) {
            for(size_t bi = ai + 1; bi < index1.size(); ++bi) {
                float dr[3];
                center_dr(cg, index1[ai], cg, index1[bi], box_x, box_y, box_z, dr);
                float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if(r2 < cg_cutoff2)
                    active_cg_pairs.push_back({ai, bi});
            }
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
    VecArray comp = compaction_state ? compaction_state->output : VecArray();
    VecArray tgt = tgt_pos.output;
    ensure_pairlist(cg, tgt);
    float total = 0.f;
    ImplicitCompactionField implicit_field;
    if(has_compaction_correction && !compaction_state) {
        implicit_field = build_gap_response_compaction_field(
                cg,
                index1,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                gap_response_coeff,
                gap_response_n_knot,
                gap_response_coord_min_ang,
                gap_response_coord_spacing_ang,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_fallback_ang,
                gap_response_smooth_weight);
    }

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
            if(has_compaction_correction) {
                TargetSplineEval e_ext = {0.f, 0.f, 0.f};
                TargetSplineEval e_comp = {0.f, 0.f, 0.f};
                bool ok_ext = eval_cg_target_tensor(
                        param_ptr(delta_extended, n_type2, n_param, t1, t2),
                        n_angular, n_radial, dr, n1, knot_spacing, cutoff,
                        (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e_ext);
                bool ok_comp = eval_cg_target_tensor(
                        param_ptr(delta_compact, n_type2, n_param, t1, t2),
                        n_angular, n_radial, dr, n1, knot_spacing, cutoff,
                        (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e_comp);
                float q = implicit_field.extended[ai];
                if(compaction_state) {
                    // The dynamic compaction state is the normalized hidden
                    // compactness coordinate on [0, 1], not the physical tail
                    // extension used to pick representative training ensembles.
                    q = comp(0, int(ai));
                }
                if(ok_ext) {
                    float w = 1.f - q;
                    e.value += w * e_ext.value;
                    e.d_dr += w * e_ext.d_dr;
                    e.d_da += w * e_ext.d_da;
                }
                if(ok_comp) {
                    e.value += q * e_comp.value;
                    e.d_dr += q * e_comp.d_dr;
                    e.d_da += q * e_comp.d_da;
                }
            }
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
    VecArray comp = compaction_state ? compaction_state->output : VecArray();
    VecArray comp_sens = compaction_state ? compaction_state->sens : VecArray();
    VecArray tgt = tgt_pos.output;
    VecArray tgt_sens = tgt_pos.sens;
    ensure_pairlist(cg, tgt);
    ImplicitCompactionField implicit_field;
    vector<float> implicit_q_sens(index1.size(), 0.f);
    if(has_compaction_correction && !compaction_state) {
        implicit_field = build_gap_response_compaction_field(
                cg,
                index1,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                gap_response_coeff,
                gap_response_n_knot,
                gap_response_coord_min_ang,
                gap_response_coord_spacing_ang,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_fallback_ang,
                gap_response_smooth_weight);
    }

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
        float q_ctrl_sens = 0.f;
        if(has_compaction_correction) {
            TargetSplineEval e_ext = {0.f, 0.f, 0.f};
            TargetSplineEval e_comp = {0.f, 0.f, 0.f};
            bool ok_ext = eval_cg_target_tensor(
                    param_ptr(delta_extended, n_type2, n_param, t1, t2),
                    n_angular, n_radial, dr, n1, knot_spacing, cutoff,
                    (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e_ext);
            bool ok_comp = eval_cg_target_tensor(
                    param_ptr(delta_compact, n_type2, n_param, t1, t2),
                    n_angular, n_radial, dr, n1, knot_spacing, cutoff,
                    (boltzmann_weight_transform || log1p_reduced_transform) ? 0.f : taper_width, e_comp);
            float q = implicit_field.extended[ai];
            float q_to_compaction_sens = 0.f;
            if(compaction_state) {
                q = comp(0, int(ai));
                q_to_compaction_sens = 1.f;
            }
            float corr_ext = ok_ext ? e_ext.value : 0.f;
            float corr_comp = ok_comp ? e_comp.value : 0.f;
            if(ok_ext) {
                float w = 1.f - q;
                e.value += w * e_ext.value;
                e.d_dr += w * e_ext.d_dr;
                e.d_da += w * e_ext.d_da;
            }
            if(ok_comp) {
                e.value += q * e_comp.value;
                e.d_dr += q * e_comp.d_dr;
                e.d_da += q * e_comp.d_da;
            }
            q_ctrl_sens = corr_comp - corr_ext;
            if(compaction_state)
                q_ctrl_sens *= q_to_compaction_sens;
        }
        if(boltzmann_weight_transform) {
            apply_boltzmann_weight_transform(
                    e, reference_energy_eup[t1 * n_type2 + t2],
                    boltzmann_temperature, minimum_boltzmann_weight,
                    dr, cutoff, taper_width);
        } else if(log1p_reduced_transform) {
            float taper = 1.f;
            float d_taper_dr = 0.f;
            compute_cutoff_taper(norm_dr(dr), cutoff, taper_width, taper, d_taper_dr);
            float control_scale = boltzmann_temperature * expf(e.value) * taper;
            q_ctrl_sens *= control_scale;
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
        if(has_compaction_correction && compaction_state)
            comp_sens(0, int(ai)) += q_ctrl_sens;
        else if(has_compaction_correction)
            implicit_q_sens[ai] += q_ctrl_sens;
    }
    if(has_compaction_correction && !compaction_state)
        propagate_gap_response_compaction_sens(
                cg,
                cg_sens,
                index1,
                active_cg_pairs,
                box_x,
                box_y,
                box_z,
                implicit_field,
                implicit_q_sens,
                gap_response_radial_cutoff_ang,
                gap_response_face_cos_min,
                gap_response_smooth_weight);
}

namespace martini_cg_lipid {

void register_dynamic_compaction_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed) {
    if(!engine) return;
    clear_dynamic_compaction_for_engine(engine);
    int idx = engine->get_idx("cgl_compaction_state", false);
    if(idx < 0) return;
    auto* state = dynamic_cast<CGLCompactionState*>(engine->nodes[idx].computation.get());
    if(!state) throw string("cgl_compaction_state node has unexpected type");
    auto runtime = std::unique_ptr<CGLDynamicCompactionRuntime>(
        new CGLDynamicCompactionRuntime(state, random_seed));

    if(h5_exists(config_root, "/input/cgl_compaction_mom")) {
        auto shape = get_dset_size(1, config_root, "/input/cgl_compaction_mom");
        if(shape[0] != hsize_t(state->n_elem))
            throw string("/input/cgl_compaction_mom shape must be n_cgl");
        traverse_dset<1, float>(config_root, "/input/cgl_compaction_mom",
                [&](size_t i, float v) {
            runtime->momentum(0, int(i)) = v;
        });
        runtime->hidden_state_initialized = true;
        runtime->hidden_state_temperature = read_attribute<float>(
            config_root, "/input/cgl_compaction_mom", "restart_temperature", -1.f);
    }
    compaction_runtime_map()[engine] = std::move(runtime);
}

void clear_dynamic_compaction_for_engine(DerivEngine* engine) {
    compaction_runtime_map().erase(engine);
}

bool has_dynamic_compaction(DerivEngine* engine) {
    return get_compaction_runtime(engine) != nullptr;
}

void set_dynamic_compaction_temperature(DerivEngine* engine, float temperature) {
    auto* runtime = get_compaction_runtime(engine);
    if(!runtime) return;
    if(temperature > 0.f) {
        if(runtime->hidden_state_initialized) {
            float old_temperature = hidden_state_reference_temperature(runtime);
            if(hidden_state_temperature_changed(old_temperature, temperature))
                rescale_dynamic_compaction_hidden_state(runtime, old_temperature, temperature);
            runtime->hidden_state_temperature = temperature;
        }
    }
    runtime->temperature = temperature;
    runtime->temperature_is_set = true;
}

void set_dynamic_compaction_thermostat_delta_t(DerivEngine* engine, float delta_t) {
    auto* runtime = get_compaction_runtime(engine);
    if(runtime) runtime->thermostat_delta_t = delta_t;
}

void apply_dynamic_compaction_thermostat(DerivEngine* engine) {
    auto* runtime = get_compaction_runtime(engine);
    if(!runtime || runtime->thermostat_delta_t <= 0.f) return;
    auto* state = runtime->state;
    float mom_scale = expf(-runtime->thermostat_delta_t / state->thermostat_timescale);
    float noise_scale = sqrtf(
        state->mass * runtime->temperature *
        std::max(0.f, 1.f - mom_scale * mom_scale));
    for(int i = 0; i < state->n_elem; ++i) {
        RandomGenerator random(
            runtime->random_seed,
            CGL_COMPACTION_THERMOSTAT_RANDOM_STREAM,
            uint32_t(i),
            runtime->thermostat_invocations);
        runtime->momentum(0, i) = runtime->momentum(0, i) * mom_scale
            + noise_scale * random.normal().x();
    }
    if(!runtime->hidden_state_initialized && runtime->temperature > 0.f) {
        runtime->hidden_state_initialized = true;
        runtime->hidden_state_temperature = runtime->temperature;
    }
    runtime->thermostat_invocations++;
}

void integrate_dynamic_compaction(DerivEngine* engine, float dt) {
    auto* runtime = get_compaction_runtime(engine);
    if(!runtime) return;
    auto* state = runtime->state;
    for(int i = 0; i < state->n_elem; ++i) {
        float p = runtime->momentum(0, i) - state->sens(0, i) * dt;
        float c = state->output(0, i) + (p / state->mass) * dt;
        if(c < state->coord_min) {
            c = state->coord_min;
            if(p < 0.f) p = 0.f;
        } else if(c > state->coord_max) {
            c = state->coord_max;
            if(p > 0.f) p = 0.f;
        }
        state->output(0, i) = c;
        runtime->momentum(0, i) = p;
    }
}

void add_dynamic_compaction_loggers(DerivEngine* engine, H5Logger& logger, bool record_momentum) {
    auto* runtime = get_compaction_runtime(engine);
    if(!runtime) return;
    auto* state = runtime->state;
    logger.add_logger<float>("cgl_compaction", {1, state->n_elem}, [state](float* buffer) {
        VecArray coord = state->output;
        for(int i = 0; i < state->n_elem; ++i)
            buffer[i] = coord(0, i);
    });
    if(record_momentum) {
        logger.add_logger<float>("cgl_compaction_mom", {1, state->n_elem}, [runtime](float* buffer) {
            VecArray mom = runtime->momentum;
            for(int i = 0; i < runtime->state->n_elem; ++i)
                buffer[i] = mom(0, i);
        });
    }
}

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
        runtime->hidden_state_initialized = true;
        runtime->hidden_state_temperature = read_attribute<float>(
            config_root, "/input/cgl_orientation_mom", "restart_temperature", -1.f);
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
    if(!runtime) return;
    if(temperature > 0.f) {
        if(runtime->hidden_state_initialized) {
            float old_temperature = hidden_state_reference_temperature(runtime);
            if(hidden_state_temperature_changed(old_temperature, temperature))
                rescale_dynamic_orientation_hidden_state(runtime, old_temperature, temperature);
            runtime->hidden_state_temperature = temperature;
        }
    }
    runtime->temperature = temperature;
    runtime->temperature_is_set = true;
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
    if(!runtime->hidden_state_initialized && runtime->temperature > 0.f) {
        runtime->hidden_state_initialized = true;
        runtime->hidden_state_temperature = runtime->temperature;
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
    vector<float> temperature_grid;
    TemperatureScaleTable coupling_scale;
    TemperatureScaleTable memory_tau_scale;
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
    bool have_temperature_grid = H5Lexists(grp, "temperature_grid", H5P_DEFAULT) > 0;
    bool have_coupling_scale = H5Lexists(grp, "coupling_scale", H5P_DEFAULT) > 0;
    bool have_memory_tau_scale = H5Lexists(grp, "memory_tau_scale", H5P_DEFAULT) > 0;
    if(!have_temperature_grid && (have_coupling_scale || have_memory_tau_scale))
        throw string("/input/cgl_gle temperature_grid is required when coupling_scale or memory_tau_scale is present");
    if(have_temperature_grid) {
        temperature_grid = read_float_dataset_1d(grp, "temperature_grid");
        coupling_scale = have_coupling_scale
            ? read_temperature_scale_table(grp, "coupling_scale", temperature_grid.size(), coupling.size())
            : make_uniform_temperature_scale(temperature_grid.size(), 1.f);
        memory_tau_scale = have_memory_tau_scale
            ? read_temperature_scale_table(grp, "memory_tau_scale", temperature_grid.size(), memory_tau.size())
            : make_uniform_temperature_scale(temperature_grid.size(), 1.f);
        if(temperature_grid.size() < 2u)
            throw string("/input/cgl_gle temperature_grid must have length at least 2");
        for(size_t i = 0; i < temperature_grid.size(); ++i) {
            if(temperature_grid[i] <= 0.f)
                throw string("/input/cgl_gle temperature_grid, coupling_scale, and memory_tau_scale must contain positive values");
            if(i && temperature_grid[i] <= temperature_grid[i-1])
                throw string("/input/cgl_gle temperature_grid must be strictly increasing");
        }
        for(float scale : coupling_scale.value)
            if(scale <= 0.f)
                throw string("/input/cgl_gle temperature_grid, coupling_scale, and memory_tau_scale must contain positive values");
        for(float scale : memory_tau_scale.value)
            if(scale <= 0.f)
                throw string("/input/cgl_gle temperature_grid, coupling_scale, and memory_tau_scale must contain positive values");
    }

    auto runtime = std::unique_ptr<CGLGLERuntime>(
        new CGLGLERuntime(
            atom_index, random_seed, memory_tau, coupling, temperature_grid, coupling_scale, memory_tau_scale));
    bool aux_momentum_nonzero = false;
    if(H5Lexists(grp, "aux_momentum", H5P_DEFAULT) > 0) {
        int n_mode = runtime->n_mode();
        if(n_mode == 1) {
            auto shape = get_dset_size(2, grp, "aux_momentum");
            if(shape[0] != hsize_t(atom_index.size()) || shape[1] != 3u)
                throw string("/input/cgl_gle/aux_momentum shape must be n_cgl x 3 for single-mode GLE");
            traverse_dset<2, float>(grp, "aux_momentum", [&](size_t i, size_t d, float v) {
                runtime->aux_momentum(int(d), runtime->aux_index(0, int(i))) = v;
                aux_momentum_nonzero = aux_momentum_nonzero || fabsf(v) > 0.f;
            });
        } else {
            auto shape = get_dset_size(3, grp, "aux_momentum");
            if(shape[0] != hsize_t(n_mode) || shape[1] != hsize_t(atom_index.size()) || shape[2] != 3u)
                throw string("/input/cgl_gle/aux_momentum shape must be n_mode x n_cgl x 3 for multi-mode GLE");
            traverse_dset<3, float>(grp, "aux_momentum", [&](size_t m, size_t i, size_t d, float v) {
                runtime->aux_momentum(int(d), runtime->aux_index(int(m), int(i))) = v;
                aux_momentum_nonzero = aux_momentum_nonzero || fabsf(v) > 0.f;
            });
        }
        bool have_restart_attr = attribute_exists_local(grp, ".", "aux_momentum_restart_source");
        if(have_restart_attr || aux_momentum_nonzero) {
            runtime->hidden_state_initialized = true;
            runtime->hidden_state_temperature = read_attribute<float>(
                grp, ".", "aux_momentum_restart_temperature", -1.f);
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
    if(!runtime) return;
    if(temperature > 0.f) {
        if(runtime->hidden_state_initialized) {
            float old_temperature = hidden_state_reference_temperature(runtime);
            if(hidden_state_temperature_changed(old_temperature, temperature))
                rescale_cgl_gle_hidden_state(runtime, old_temperature, temperature);
            runtime->hidden_state_temperature = temperature;
        }
    }
    runtime->temperature = temperature;
    runtime->temperature_is_set = true;
    update_cgl_gle_transport_for_temperature(runtime);
}

void set_cgl_gle_delta_t(DerivEngine* engine, float delta_t) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(!runtime) return;
    runtime->delta_t = delta_t;
    update_cgl_gle_linear_step(runtime);
}

void apply_cgl_gle_thermostat(DerivEngine* engine, VecArray mom) {
    auto* runtime = get_cgl_gle_runtime(engine);
    if(!runtime || runtime->delta_t <= 0.f || !runtime->linear_step_ready) return;

    bool have_masses = martini_masses::has_masses(engine);
    int n_state = runtime->linear_step_dim();
    vector<double> state(size_t(n_state), 0.0);
    vector<double> noise(size_t(n_state), 0.0);
    vector<double> updated(size_t(n_state), 0.0);

    for(size_t i = 0; i < runtime->atom_index.size(); ++i) {
        int atom = runtime->atom_index[i];
        if(atom < 0 || atom >= engine->pos->n_atom) continue;
        float mass = have_masses ? martini_masses::get_mass(engine, atom) : 1.f;
        float3 p = load_vec<3>(mom, atom);
        double noise_scale = sqrt(std::max(0.f, runtime->temperature) * std::max(mass, 0.f));
        uint32_t random_atom = uint32_t(atom);
        RandomGenerator random(
            runtime->random_seed,
            CGL_GLE_THERMOSTAT_RANDOM_STREAM,
            random_atom,
            runtime->thermostat_invocations);

        float3 next_p = p;
        for(int d = 0; d < 3; ++d) {
            state[0] = p[d];
            for(int mode = 0; mode < runtime->n_mode(); ++mode)
                state[size_t(mode + 1)] = runtime->aux_momentum(d, runtime->aux_index(mode, int(i)));

            int noise_pos = 0;
            while(noise_pos < n_state) {
                float4 z = random.normal();
                noise[size_t(noise_pos++)] = z.x();
                if(noise_pos < n_state) noise[size_t(noise_pos++)] = z.y();
                if(noise_pos < n_state) noise[size_t(noise_pos++)] = z.z();
                if(noise_pos < n_state) noise[size_t(noise_pos++)] = z.w();
            }

            for(int row = 0; row < n_state; ++row) {
                double value = 0.0;
                for(int col = 0; col < n_state; ++col)
                    value += dense_square_at(runtime->linear_step_matrix, n_state, row, col) * state[size_t(col)];
                double stochastic = 0.0;
                for(int col = 0; col <= row; ++col)
                    stochastic += dense_square_at(runtime->noise_cholesky_unit, n_state, row, col) * noise[size_t(col)];
                updated[size_t(row)] = value + noise_scale * stochastic;
            }

            next_p[d] = float(updated[0]);
            for(int mode = 0; mode < runtime->n_mode(); ++mode)
                runtime->aux_momentum(d, runtime->aux_index(mode, int(i))) = float(updated[size_t(mode + 1)]);
        }
        p = next_p;
        store_vec(mom, atom, p);
    }
    if(!runtime->hidden_state_initialized && runtime->temperature > 0.f) {
        runtime->hidden_state_initialized = true;
        runtime->hidden_state_temperature = runtime->temperature;
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

static RegisterNodeType<CGLCompactionState, 0> _reg_cg_compaction("cgl_compaction_state");
static RegisterNodeType<CGLOrientationState, 0> _reg_cg_orientation("cgl_orientation_state");
static RegisterNodeType<ComposeVector6D, -1> _reg_cv("compose_vector6d");
static RegisterNodeType<CGLipidCompactionSelfPotential, 1> _reg_cg_compaction_self("cg_lipid_compaction_self");
static RegisterNodeType<CGLipidPairPotential, -1> _reg_cg_pair("cg_lipid_pair");
static RegisterNodeType<CGLipidDensityPotential, 1> _reg_cg_density("cg_lipid_density");
static RegisterNodeType<CGLipidContactEmbeddingPotential, 1> _reg_cg_contact_embedding("cg_lipid_contact_embedding");
static RegisterNodeType<CGLipidGapEmbeddingPotential, 1> _reg_cg_gap_embedding("cg_lipid_gap_embedding");
static RegisterNodeType<CGLipidSCPotential, 2> _reg_cg_sc("cg_lipid_sc");
static RegisterNodeType<CGLipidSCOneBody, -1> _reg_cg_sc_1body("cg_lipid_rotamer_sc");
static RegisterNodeType<CGLipidTargetPotential, -1> _reg_cg_target("cg_lipid_target");
