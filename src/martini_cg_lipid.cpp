#include "deriv_engine.h"
#include "spline.h"
#include "box.h"
#include "state_logger.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
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

static vector<int> read_int_dataset(hid_t grp, const char* name) {
    vector<hsize_t> sz = get_dset_size(1, grp, name);
    vector<int> out(sz[0]);
    traverse_dset<1, int32_t>(grp, name, [&](size_t i, int32_t v) {
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
    Vec<2> vr = clamped_deBoor_value_and_deriv(
            p + 2 * CG_LIPID_N_ANGULAR, radial_coord, CG_LIPID_N_RADIAL);
    Vec<2> va = clamped_deBoor_value_and_deriv(
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

    Vec<2> v0 = clamped_deBoor_value_and_deriv(p, radial_coord, n_radial);
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
        Vec<2> vm = clamped_deBoor_value_and_deriv(v_ptr, radial_coord, n_radial);
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

}  // namespace

// ===================================================================
// ComposeVector6D — Combines 3D positions with static 3D directions
// ===================================================================
struct ComposeVector6D : public CoordNode {
    CoordNode& pos;
    vector<int> elem_index;
    vector<int> orientation_index;
    vector<float> direction;
    bool dynamic_orientation;
    float box_x;
    float box_y;
    float box_z;

    ComposeVector6D(hid_t grp, CoordNode& pos_node);

    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
};

ComposeVector6D::ComposeVector6D(hid_t grp, CoordNode& pos_node)
    : CoordNode(h5::get_dset_size(1, grp, "elem_index")[0], 6)
    , pos(pos_node)
    , dynamic_orientation(h5_exists(grp, "orientation_index"))
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

void ComposeVector6D::compute_value(ComputeMode) {
    VecArray posc = pos.output;
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        output(0, i) = posc(0, ai);
        output(1, i) = posc(1, ai);
        output(2, i) = posc(2, ai);
        if(dynamic_orientation) {
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
    for (int i = 0; i < n_elem; i++) {
        int ai = elem_index[i];
        pos_sens(0, ai) += sens(0, i);
        pos_sens(1, ai) += sens(1, i);
        pos_sens(2, ai) += sens(2, i);

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

// ===================================================================
// CGLipidPairPotential — directional residual CG lipid ↔ CG lipid term
// ===================================================================
struct CGLipidPairPotential : public PotentialNode {
    CoordNode& cg_pos;
    vector<float> interaction_param;
    vector<int> index;
    vector<int> type;
    vector<int> id;
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

    CGLipidPairPotential(hid_t grp, CoordNode& cg_pos);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == interaction_param.size()) interaction_param = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
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
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", CG_LIPID_DEFAULT_KNOT_SPACING))
    , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                float(n_radial - 2) * CG_LIPID_DEFAULT_KNOT_SPACING))
    , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
{
    check_elem_width(cg_pos, 6);
    H5Obj pi_obj = open_group(grp, "pair_interaction");
    hid_t pi = pi_obj.get();
    interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
    int expected_n_param = n_radial + n_modes * (2 * n_angular + n_radial);
    if(n_modes <= 0 || n_radial <= 3 || n_angular <= 3 || n_param != expected_n_param)
        throw string("cg_lipid_pair requires full multimode params with matching n_modes/n_radial/n_angular attrs");
    index = read_int_dataset(pi, "index");
    type = read_int_dataset(pi, "type");
    id = read_int_dataset(pi, "id");
    if(index.size() != type.size() || index.size() != id.size())
        throw string("cg_lipid_pair index/type/id size mismatch");
}

void CGLipidPairPotential::compute_value(ComputeMode mode) {
    (void)mode;
    VecArray cg = cg_pos.output;
    float total = 0.f;
    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if((id[ai] >> 4) == (id[bi] >> 4)) continue;
            int t1 = type[ai];
            int t2 = type[bi];
            if(t1 < 0 || t1 >= n_type1 || t2 < 0 || t2 >= n_type2) continue;

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(cg, index[ai], x1, n1);
            load_vec6(cg, index[bi], x2, n2);
            dr[0] = x2[0] - x1[0];
            dr[1] = x2[1] - x1[1];
            dr[2] = x2[2] - x1[2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
                simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);

            QuadsplineEval e;
            if(eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                        n_modes, n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e))
                total += e.value;
        }
    }
    potential = total;
}

void CGLipidPairPotential::propagate_deriv() {
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;

    for(size_t ai = 0; ai < index.size(); ++ai) {
        for(size_t bi = ai + 1; bi < index.size(); ++bi) {
            if((id[ai] >> 4) == (id[bi] >> 4)) continue;
            int t1 = type[ai];
            int t2 = type[bi];
            if(t1 < 0 || t1 >= n_type1 || t2 < 0 || t2 >= n_type2) continue;

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(cg, index[ai], x1, n1);
            load_vec6(cg, index[bi], x2, n2);
            dr[0] = x2[0] - x1[0];
            dr[1] = x2[1] - x1[1];
            dr[2] = x2[2] - x1[2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
                simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);

            QuadsplineEval e;
            if(!eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                        n_modes, n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e))
                continue;

            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(cg_sens, index[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index[bi], dpos2, ddir2);
        }
    }
}

// ===================================================================
// CGLipidSCPotential — directional sidechain ↔ CG lipid term
// ===================================================================
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

    CGLipidSCPotential(hid_t grp, CoordNode& sc_pos, CoordNode& cg_pos);
    virtual void compute_value(ComputeMode mode) override;
    virtual void propagate_deriv() override;
    virtual void set_param(const vector<float>& new_param) override {
        if(new_param.size() == interaction_param.size()) interaction_param = new_param;
    }
    virtual void update_box_dimensions_anisotropic(float scale_xy, float scale_z) override {
        box_x *= scale_xy;
        box_y *= scale_xy;
        box_z *= scale_z;
    }
};

CGLipidSCPotential::CGLipidSCPotential(hid_t grp, CoordNode& sc_pos_, CoordNode& cg_pos_)
    : PotentialNode()
    , sc_pos(sc_pos_)
    , cg_pos(cg_pos_)
    , n_type1(0)
    , n_type2(0)
    , n_param(0)
    , n_modes(read_attribute<int>(grp, ".", "n_modes", 0))
    , n_radial(read_attribute<int>(grp, ".", "n_radial", CG_LIPID_N_RADIAL))
    , n_angular(read_attribute<int>(grp, ".", "n_angular", CG_LIPID_N_ANGULAR))
    , box_x(read_attribute<float>(grp, ".", "x_len", 0.f))
    , box_y(read_attribute<float>(grp, ".", "y_len", 0.f))
    , box_z(read_attribute<float>(grp, ".", "z_len", 0.f))
    , knot_spacing(read_attribute<float>(grp, ".", "knot_spacing_ang", CG_LIPID_DEFAULT_KNOT_SPACING))
    , cutoff(read_attribute<float>(grp, ".", "cutoff_ang",
                (CG_LIPID_N_RADIAL - 2) * CG_LIPID_DEFAULT_KNOT_SPACING))
    , taper_width(read_attribute<float>(grp, ".", "taper_width_ang", knot_spacing))
{
    check_elem_width(sc_pos, 6);
    check_elem_width(cg_pos, 6);
    H5Obj pi_obj = open_group(grp, "pair_interaction");
    hid_t pi = pi_obj.get();
    interaction_param = read_param_dataset_any(pi, n_type1, n_type2, n_param);
    int expected_n_param = n_radial + n_modes * (2 * n_angular + n_radial);
    if(n_modes > 0) {
        if(n_radial <= 3 || n_angular <= 3 || n_param != expected_n_param)
            throw string("cg_lipid_sc full multimode params require matching n_modes/n_radial/n_angular attrs");
    } else if(n_param != CG_LIPID_N_PARAM) {
        throw string("cg_lipid_sc legacy params must have last dimension 54");
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

void CGLipidSCPotential::compute_value(ComputeMode mode) {
    (void)mode;
    VecArray sc = sc_pos.output;
    VecArray cg = cg_pos.output;
    float total = 0.f;

    for(size_t ai = 0; ai < index1.size(); ++ai) {
        for(size_t bi = 0; bi < index2.size(); ++bi) {
            if((id1[ai] >> 4) == (id2[bi] >> 4)) continue;
            int t1 = type1[ai];
            int t2 = type2[bi];
            if(t1 < 0 || t1 >= n_type1 || t2 < 0 || t2 >= n_type2) continue;

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(sc, index1[ai], x1, n1);
            load_vec6(cg, index2[bi], x2, n2);
            dr[0] = x2[0] - x1[0];
            dr[1] = x2[1] - x1[1];
            dr[2] = x2[2] - x1[2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
                simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);

            QuadsplineEval e;
            bool ok = false;
            if(n_modes > 0) {
                ok = eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                        n_modes, n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e);
            } else {
                ok = eval_quadspline(param_ptr(interaction_param, n_type2, t1, t2),
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e);
            }
            if(ok) total += e.value;
        }
    }
    potential = total;
}

void CGLipidSCPotential::propagate_deriv() {
    VecArray sc = sc_pos.output;
    VecArray sc_sens = sc_pos.sens;
    VecArray cg = cg_pos.output;
    VecArray cg_sens = cg_pos.sens;

    for(size_t ai = 0; ai < index1.size(); ++ai) {
        for(size_t bi = 0; bi < index2.size(); ++bi) {
            if((id1[ai] >> 4) == (id2[bi] >> 4)) continue;
            int t1 = type1[ai];
            int t2 = type2[bi];
            if(t1 < 0 || t1 >= n_type1 || t2 < 0 || t2 >= n_type2) continue;

            float x1[3], n1[3], x2[3], n2[3], dr[3];
            load_vec6(sc, index1[ai], x1, n1);
            load_vec6(cg, index2[bi], x2, n2);
            dr[0] = x2[0] - x1[0];
            dr[1] = x2[1] - x1[1];
            dr[2] = x2[2] - x1[2];
            if(box_x > 0.f && box_y > 0.f && box_z > 0.f)
                simulation_box::minimum_image_scalar(dr[0], dr[1], dr[2], box_x, box_y, box_z);

            QuadsplineEval e;
            bool ok = false;
            if(n_modes > 0) {
                ok = eval_multimode_pair(param_ptr(interaction_param, n_type2, n_param, t1, t2),
                        n_modes, n_angular, n_radial,
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e);
            } else {
                ok = eval_quadspline(param_ptr(interaction_param, n_type2, t1, t2),
                        dr, n1, n2, knot_spacing, cutoff, taper_width, e);
            }
            if(!ok)
                continue;

            float dpos1[3] = {0.f, 0.f, 0.f};
            float ddir1[3] = {0.f, 0.f, 0.f};
            float dpos2[3] = {0.f, 0.f, 0.f};
            float ddir2[3] = {0.f, 0.f, 0.f};
            accumulate_deriv(dr, n1, n2, e, dpos1, ddir1, dpos2, ddir2);
            add_vec6_sens(sc_sens, index1[ai], dpos1, ddir1);
            add_vec6_sens(cg_sens, index2[bi], dpos2, ddir2);
        }
    }
}

// ===================================================================
// OrientationSpring — penalizes angular deviation from reference direction
// ===================================================================
struct OrientationSpring : public PotentialNode {
    CoordNode& cg_pos;
    int n_elem;
    vector<float> ref_dir;
    float k_orient;

    OrientationSpring(hid_t grp, CoordNode& cg_pos_);
    virtual void compute_value(ComputeMode mode) override;
};

OrientationSpring::OrientationSpring(hid_t grp, CoordNode& cg_pos_)
    : PotentialNode()
    , cg_pos(cg_pos_)
    , n_elem(get_dset_size(2, grp, "ref_dir")[0])
    , ref_dir(n_elem * 3)
    , k_orient(read_attribute<float>(grp, ".", "k_orient", 0.f))
{
    check_elem_width(cg_pos, 6);
    traverse_dset<2, float>(grp, "ref_dir", [&](size_t i, size_t d, float v) {
        ref_dir[i * 3 + d] = v;
    });
}

void OrientationSpring::compute_value(ComputeMode mode) {
    float* pot = (mode == PotentialAndDerivMode) ? &potential : nullptr;
    if(pot) *pot = 0.f;

    for(int i = 0; i < n_elem; ++i) {
        float n[3] = {cg_pos.output(3,i), cg_pos.output(4,i), cg_pos.output(5,i)};
        float nr[3] = {ref_dir[i*3], ref_dir[i*3+1], ref_dir[i*3+2]};
        float dot = n[0]*nr[0] + n[1]*nr[1] + n[2]*nr[2];
        if(pot) *pot += k_orient * (1.f - dot);
        // dE/dn = -k_orient * n_ref  (energy does not depend on position)
        cg_pos.sens(3,i) += -k_orient * nr[0];
        cg_pos.sens(4,i) += -k_orient * nr[1];
        cg_pos.sens(5,i) += -k_orient * nr[2];
    }
}

// ===================================================================
// Node registration
// ===================================================================
static RegisterNodeType<ComposeVector6D, 1> _reg_cv("compose_vector6d");
static RegisterNodeType<CGLipidPairPotential, 1> _reg_cg_pair("cg_lipid_pair");
static RegisterNodeType<CGLipidSCPotential, 2> _reg_cg_sc("cg_lipid_sc");
static RegisterNodeType<OrientationSpring, 1> _reg_orient("cg_lipid_orientation_spring");
