#include "martini.h"
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <mutex>
#include "spline.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists
#include <cmath> // For pow, cosf, sinf, acosf
#include <cctype>
#include <cstdint>
#include <set> // For std::set
#include <array>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <limits>
#include "box.h" // For PBC minimum_image function

using namespace h5;
using namespace std;

namespace martini_fix_rigid {

// Global registry for fix rigid constraints per engine
static std::mutex g_fix_rigid_mutex;
static std::map<DerivEngine*, std::vector<int>> g_user_fixed_atoms;
static std::map<DerivEngine*, std::vector<int>> g_dynamic_fixed_atoms;
static std::map<DerivEngine*, std::vector<int>> g_dynamic_z_fixed_atoms;
static std::map<DerivEngine*, std::vector<int>> g_fixed_atoms;
static std::map<DerivEngine*, std::vector<int>> g_z_fixed_atoms;
static std::map<DerivEngine*, std::vector<float>> g_atom_masses;

struct RigidBodyGroup {
    std::vector<int> atom_indices;
    std::vector<float> masses;
    std::vector<std::array<float,3>> reference_rel;
    bool initialized = false;
};
static std::map<DerivEngine*, std::vector<std::vector<int>>> g_dynamic_rigid_groups;
static std::map<DerivEngine*, std::vector<RigidBodyGroup>> g_rigid_groups;
static void rebuild_rigid_groups(DerivEngine& engine);

static void normalize_atom_list(std::vector<int>& atoms) {
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
}

static void rebuild_fixed_atoms(DerivEngine& engine) {
    auto& merged = g_fixed_atoms[&engine];
    merged.clear();
    auto uit = g_user_fixed_atoms.find(&engine);
    if(uit != g_user_fixed_atoms.end()) {
        merged.insert(merged.end(), uit->second.begin(), uit->second.end());
    }
    auto dit = g_dynamic_fixed_atoms.find(&engine);
    if(dit != g_dynamic_fixed_atoms.end()) {
        merged.insert(merged.end(), dit->second.begin(), dit->second.end());
    }
    normalize_atom_list(merged);

    auto& merged_z = g_z_fixed_atoms[&engine];
    merged_z.clear();
    auto zdit = g_dynamic_z_fixed_atoms.find(&engine);
    if(zdit != g_dynamic_z_fixed_atoms.end()) {
        merged_z.insert(merged_z.end(), zdit->second.begin(), zdit->second.end());
    }
    normalize_atom_list(merged_z);
    if(!merged.empty() && !merged_z.empty()) {
        std::vector<int> filtered;
        filtered.reserve(merged_z.size());
        for(int atom_idx : merged_z) {
            if(!std::binary_search(merged.begin(), merged.end(), atom_idx)) {
                filtered.push_back(atom_idx);
            }
        }
        merged_z.swap(filtered);
    }
}

static void merge_fixed_atoms(DerivEngine& engine, const std::vector<int>& extra_atoms) {
    if(extra_atoms.empty()) return;
    auto& fixed_atoms = g_user_fixed_atoms[&engine];
    fixed_atoms.insert(fixed_atoms.end(), extra_atoms.begin(), extra_atoms.end());
    normalize_atom_list(fixed_atoms);
    rebuild_fixed_atoms(engine);
}

static bool invert_3x3(const float a[3][3], float inv_out[3][3]) {
    float det =
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if(std::fabs(det) < 1.0e-12f) return false;
    float inv_det = 1.0f / det;
    inv_out[0][0] =  (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det;
    inv_out[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) * inv_det;
    inv_out[0][2] =  (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det;
    inv_out[1][0] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) * inv_det;
    inv_out[1][1] =  (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det;
    inv_out[1][2] = -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) * inv_det;
    inv_out[2][0] =  (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det;
    inv_out[2][1] = -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) * inv_det;
    inv_out[2][2] =  (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det;
    return true;
}

static std::array<float,3> cross3(
        const std::array<float,3>& a,
        const std::array<float,3>& b) {
    return std::array<float,3>{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

static std::array<float,3> mat3_mul_vec(
        const float m[3][3],
        const std::array<float,3>& v) {
    return std::array<float,3>{
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    };
}

static void quat_to_rotmat(const float q[4], float R[3][3]) {
    float w = q[0], x = q[1], y = q[2], z = q[3];
    float ww = w*w, xx = x*x, yy = y*y, zz = z*z;
    float wx = w*x, wy = w*y, wz = w*z;
    float xy = x*y, xz = x*z, yz = y*z;
    R[0][0] = ww + xx - yy - zz;
    R[0][1] = 2.f * (xy - wz);
    R[0][2] = 2.f * (xz + wy);
    R[1][0] = 2.f * (xy + wz);
    R[1][1] = ww - xx + yy - zz;
    R[1][2] = 2.f * (yz - wx);
    R[2][0] = 2.f * (xz - wy);
    R[2][1] = 2.f * (yz + wx);
    R[2][2] = ww - xx - yy + zz;
}

static std::array<float,3> apply_rot3(const float R[3][3], const std::array<float,3>& v) {
    return std::array<float,3>{
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
    };
}

static bool best_fit_rotation_from_reference(
        const std::vector<std::array<float,3>>& reference_rel,
        const std::vector<std::array<float,3>>& current_rel,
        const std::vector<float>& masses,
        float R[3][3]) {
    if(reference_rel.size() != current_rel.size() || reference_rel.size() != masses.size()) {
        return false;
    }
    if(reference_rel.size() < 2) {
        R[0][0] = 1.f; R[0][1] = 0.f; R[0][2] = 0.f;
        R[1][0] = 0.f; R[1][1] = 1.f; R[1][2] = 0.f;
        R[2][0] = 0.f; R[2][1] = 0.f; R[2][2] = 1.f;
        return true;
    }

    float Sxx = 0.f, Sxy = 0.f, Sxz = 0.f;
    float Syx = 0.f, Syy = 0.f, Syz = 0.f;
    float Szx = 0.f, Szy = 0.f, Szz = 0.f;
    float total_weight = 0.f;
    for(size_t i = 0; i < reference_rel.size(); ++i) {
        float w = masses[i];
        if(!(w > 0.f)) w = 1.f;
        const auto& r = reference_rel[i];
        const auto& c = current_rel[i];
        total_weight += w;
        Sxx += w * c[0] * r[0];
        Sxy += w * c[0] * r[1];
        Sxz += w * c[0] * r[2];
        Syx += w * c[1] * r[0];
        Syy += w * c[1] * r[1];
        Syz += w * c[1] * r[2];
        Szx += w * c[2] * r[0];
        Szy += w * c[2] * r[1];
        Szz += w * c[2] * r[2];
    }
    if(!(total_weight > 0.f)) return false;

    float N[4][4];
    N[0][0] = Sxx + Syy + Szz;
    N[0][1] = Syz - Szy;
    N[0][2] = Szx - Sxz;
    N[0][3] = Sxy - Syx;
    N[1][0] = Syz - Szy;
    N[1][1] = Sxx - Syy - Szz;
    N[1][2] = Sxy + Syx;
    N[1][3] = Szx + Sxz;
    N[2][0] = Szx - Sxz;
    N[2][1] = Sxy + Syx;
    N[2][2] = -Sxx + Syy - Szz;
    N[2][3] = Syz + Szy;
    N[3][0] = Sxy - Syx;
    N[3][1] = Szx + Sxz;
    N[3][2] = Syz + Szy;
    N[3][3] = -Sxx - Syy + Szz;

    float q[4] = {1.f, 0.f, 0.f, 0.f};
    for(int iter = 0; iter < 16; ++iter) {
        float qn[4] = {
            N[0][0]*q[0] + N[0][1]*q[1] + N[0][2]*q[2] + N[0][3]*q[3],
            N[1][0]*q[0] + N[1][1]*q[1] + N[1][2]*q[2] + N[1][3]*q[3],
            N[2][0]*q[0] + N[2][1]*q[1] + N[2][2]*q[2] + N[2][3]*q[3],
            N[3][0]*q[0] + N[3][1]*q[1] + N[3][2]*q[2] + N[3][3]*q[3],
        };
        float norm = sqrtf(qn[0]*qn[0] + qn[1]*qn[1] + qn[2]*qn[2] + qn[3]*qn[3]);
        if(norm <= 1.0e-12f) return false;
        q[0] = qn[0] / norm;
        q[1] = qn[1] / norm;
        q[2] = qn[2] / norm;
        q[3] = qn[3] / norm;
    }
    quat_to_rotmat(q, R);
    return true;
}

static void initialize_rigid_group_if_needed(DerivEngine& engine, RigidBodyGroup& group, VecArray pos) {
    if(group.initialized) return;
    float total_mass = 0.f;
    std::array<float,3> com{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        total_mass += m;
        com[0] += m * pos(0, atom);
        com[1] += m * pos(1, atom);
        com[2] += m * pos(2, atom);
    }
    if(!(total_mass > 0.f)) return;
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;
    group.reference_rel.assign(group.atom_indices.size(), std::array<float,3>{{0.f,0.f,0.f}});
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        group.reference_rel[i][0] = pos(0, atom) - com[0];
        group.reference_rel[i][1] = pos(1, atom) - com[1];
        group.reference_rel[i][2] = pos(2, atom) - com[2];
    }
    group.initialized = true;
}

static void apply_rigid_group_constraints(
        DerivEngine& engine,
        RigidBodyGroup& group,
        VecArray pos,
        VecArray deriv,
        VecArray mom) {
    if(group.atom_indices.size() < 2) return;
    initialize_rigid_group_if_needed(engine, group, pos);
    if(!group.initialized || group.reference_rel.size() != group.atom_indices.size()) return;

    float total_mass = 0.f;
    std::array<float,3> com{{0.f, 0.f, 0.f}};
    std::array<float,3> total_force{{0.f, 0.f, 0.f}};
    std::array<float,3> total_momentum{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        total_mass += m;
        com[0] += m * pos(0, atom);
        com[1] += m * pos(1, atom);
        com[2] += m * pos(2, atom);
        total_force[0] += -deriv(0, atom);
        total_force[1] += -deriv(1, atom);
        total_force[2] += -deriv(2, atom);
        total_momentum[0] += mom(0, atom);
        total_momentum[1] += mom(1, atom);
        total_momentum[2] += mom(2, atom);
    }
    if(!(total_mass > 0.f)) return;
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;

    std::vector<std::array<float,3>> current_rel(group.atom_indices.size(), std::array<float,3>{{0.f,0.f,0.f}});
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        current_rel[i][0] = pos(0, atom) - com[0];
        current_rel[i][1] = pos(1, atom) - com[1];
        current_rel[i][2] = pos(2, atom) - com[2];
    }
    float Rfit[3][3];
    bool have_rotation = best_fit_rotation_from_reference(group.reference_rel, current_rel, group.masses, Rfit);
    if(have_rotation) {
        for(size_t i = 0; i < group.atom_indices.size(); ++i) {
            int atom = group.atom_indices[i];
            if(atom < 0 || atom >= engine.pos->n_atom) continue;
            auto rfit = apply_rot3(Rfit, group.reference_rel[i]);
            pos(0, atom) = com[0] + rfit[0];
            pos(1, atom) = com[1] + rfit[1];
            pos(2, atom) = com[2] + rfit[2];
            current_rel[i] = rfit;
        }
    }

    float inertia[3][3] = {{0.f,0.f,0.f},{0.f,0.f,0.f},{0.f,0.f,0.f}};
    std::array<float,3> torque{{0.f, 0.f, 0.f}};
    std::array<float,3> ang_momentum{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        std::array<float,3> r = current_rel[i];
        std::array<float,3> f{
            -deriv(0, atom),
            -deriv(1, atom),
            -deriv(2, atom),
        };
        std::array<float,3> p{
            mom(0, atom),
            mom(1, atom),
            mom(2, atom),
        };
        auto rx_f = cross3(r, f);
        auto rx_p = cross3(r, p);
        torque[0] += rx_f[0];
        torque[1] += rx_f[1];
        torque[2] += rx_f[2];
        ang_momentum[0] += rx_p[0];
        ang_momentum[1] += rx_p[1];
        ang_momentum[2] += rx_p[2];

        float r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        inertia[0][0] += m * (r2 - r[0]*r[0]);
        inertia[1][1] += m * (r2 - r[1]*r[1]);
        inertia[2][2] += m * (r2 - r[2]*r[2]);
        inertia[0][1] -= m * r[0] * r[1];
        inertia[1][0] -= m * r[0] * r[1];
        inertia[0][2] -= m * r[0] * r[2];
        inertia[2][0] -= m * r[0] * r[2];
        inertia[1][2] -= m * r[1] * r[2];
        inertia[2][1] -= m * r[1] * r[2];
    }

    float inertia_inv[3][3];
    bool invertible = invert_3x3(inertia, inertia_inv);
    std::array<float,3> alpha{{0.f, 0.f, 0.f}};
    std::array<float,3> omega{{0.f, 0.f, 0.f}};
    if(invertible) {
        alpha = mat3_mul_vec(inertia_inv, torque);
        omega = mat3_mul_vec(inertia_inv, ang_momentum);
    }
    std::array<float,3> vcm{
        total_momentum[0] / total_mass,
        total_momentum[1] / total_mass,
        total_momentum[2] / total_mass,
    };
    std::array<float,3> acom{
        total_force[0] / total_mass,
        total_force[1] / total_mass,
        total_force[2] / total_mass,
    };

    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        std::array<float,3> r = current_rel[i];
        auto alpha_x_r = cross3(alpha, r);
        auto omega_x_r = cross3(omega, r);
        std::array<float,3> rigid_force{
            m * (acom[0] + alpha_x_r[0]),
            m * (acom[1] + alpha_x_r[1]),
            m * (acom[2] + alpha_x_r[2]),
        };
        std::array<float,3> rigid_mom{
            m * (vcm[0] + omega_x_r[0]),
            m * (vcm[1] + omega_x_r[1]),
            m * (vcm[2] + omega_x_r[2]),
        };
        deriv(0, atom) = -rigid_force[0];
        deriv(1, atom) = -rigid_force[1];
        deriv(2, atom) = -rigid_force[2];
        mom(0, atom) = rigid_mom[0];
        mom(1, atom) = rigid_mom[1];
        mom(2, atom) = rigid_mom[2];
    }
}

static void apply_rigid_group_minimization(
        DerivEngine& engine,
        RigidBodyGroup& group,
        VecArray pos,
        VecArray deriv) {
    if(group.atom_indices.size() < 2) return;
    initialize_rigid_group_if_needed(engine, group, pos);
    if(!group.initialized) return;

    float total_mass = 0.f;
    std::array<float,3> com{{0.f, 0.f, 0.f}};
    std::array<float,3> total_force{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        total_mass += m;
        com[0] += m * pos(0, atom);
        com[1] += m * pos(1, atom);
        com[2] += m * pos(2, atom);
        total_force[0] += -deriv(0, atom);
        total_force[1] += -deriv(1, atom);
        total_force[2] += -deriv(2, atom);
    }
    if(!(total_mass > 0.f)) return;
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;

    std::vector<std::array<float,3>> current_rel(group.atom_indices.size(), std::array<float,3>{{0.f,0.f,0.f}});
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        current_rel[i][0] = pos(0, atom) - com[0];
        current_rel[i][1] = pos(1, atom) - com[1];
        current_rel[i][2] = pos(2, atom) - com[2];
    }
    float Rfit[3][3];
    bool have_rotation = best_fit_rotation_from_reference(group.reference_rel, current_rel, group.masses, Rfit);
    if(have_rotation) {
        for(size_t i = 0; i < group.atom_indices.size(); ++i) {
            int atom = group.atom_indices[i];
            if(atom < 0 || atom >= engine.pos->n_atom) continue;
            auto rfit = apply_rot3(Rfit, group.reference_rel[i]);
            pos(0, atom) = com[0] + rfit[0];
            pos(1, atom) = com[1] + rfit[1];
            pos(2, atom) = com[2] + rfit[2];
            current_rel[i] = rfit;
        }
    }

    float inertia[3][3] = {{0.f,0.f,0.f},{0.f,0.f,0.f},{0.f,0.f,0.f}};
    std::array<float,3> torque{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        std::array<float,3> r = current_rel[i];
        std::array<float,3> f{
            -deriv(0, atom),
            -deriv(1, atom),
            -deriv(2, atom),
        };
        auto rx_f = cross3(r, f);
        torque[0] += rx_f[0];
        torque[1] += rx_f[1];
        torque[2] += rx_f[2];
        float r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        inertia[0][0] += m * (r2 - r[0]*r[0]);
        inertia[1][1] += m * (r2 - r[1]*r[1]);
        inertia[2][2] += m * (r2 - r[2]*r[2]);
        inertia[0][1] -= m * r[0] * r[1];
        inertia[1][0] -= m * r[0] * r[1];
        inertia[0][2] -= m * r[0] * r[2];
        inertia[2][0] -= m * r[0] * r[2];
        inertia[1][2] -= m * r[1] * r[2];
        inertia[2][1] -= m * r[1] * r[2];
    }
    float inertia_inv[3][3];
    std::array<float,3> alpha{{0.f, 0.f, 0.f}};
    if(invert_3x3(inertia, inertia_inv)) {
        alpha = mat3_mul_vec(inertia_inv, torque);
    }
    std::array<float,3> acom{
        total_force[0] / total_mass,
        total_force[1] / total_mass,
        total_force[2] / total_mass,
    };
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        std::array<float,3> r = current_rel[i];
        auto alpha_x_r = cross3(alpha, r);
        deriv(0, atom) = -m * (acom[0] + alpha_x_r[0]);
        deriv(1, atom) = -m * (acom[1] + alpha_x_r[1]);
        deriv(2, atom) = -m * (acom[2] + alpha_x_r[2]);
    }
}

static void project_rigid_group_geometry(
        DerivEngine& engine,
        RigidBodyGroup& group,
        VecArray pos) {
    if(group.atom_indices.size() < 2) return;
    initialize_rigid_group_if_needed(engine, group, pos);
    if(!group.initialized) return;

    float total_mass = 0.f;
    std::array<float,3> com{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        float m = group.masses[i];
        if(!(m > 0.f)) m = 1.f;
        total_mass += m;
        com[0] += m * pos(0, atom);
        com[1] += m * pos(1, atom);
        com[2] += m * pos(2, atom);
    }
    if(!(total_mass > 0.f)) return;
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;

    std::vector<std::array<float,3>> current_rel(group.atom_indices.size(), std::array<float,3>{{0.f,0.f,0.f}});
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        current_rel[i][0] = pos(0, atom) - com[0];
        current_rel[i][1] = pos(1, atom) - com[1];
        current_rel[i][2] = pos(2, atom) - com[2];
    }

    float Rfit[3][3];
    bool have_rotation = best_fit_rotation_from_reference(group.reference_rel, current_rel, group.masses, Rfit);
    if(!have_rotation) return;
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        if(atom < 0 || atom >= engine.pos->n_atom) continue;
        auto rfit = apply_rot3(Rfit, group.reference_rel[i]);
        pos(0, atom) = com[0] + rfit[0];
        pos(1, atom) = com[1] + rfit[1];
        pos(2, atom) = com[2] + rfit[2];
    }
}

// Read fix rigid settings from H5 configuration
std::vector<int> read_fix_rigid_settings(hid_t root) {
    std::vector<int> fixed_atoms;
    try {
        if(h5_exists(root, "/input/fix_rigid")) {
            auto grp = open_group(root, "/input/fix_rigid");
            int enable = read_attribute<int>(grp.get(), ".", "enable", 0);
            if(enable) {
                // Read atom indices to fix
                if(h5_exists(grp.get(), "atom_indices")) {
                    traverse_dset<1,int>(grp.get(), "atom_indices", [&](size_t i, int atom_idx) {
                        fixed_atoms.push_back(atom_idx);
                    });
                }
            }
        }
    } catch(...) { 
        // Return empty vector if no fix rigid settings found
    }
    return fixed_atoms;
}

std::vector<int> read_martini_backbone_hold(hid_t root, const std::string& atom_role_name) {
    std::vector<int> fixed_atoms;
    try {
        if(atom_role_name.empty()) return fixed_atoms;
        if(!h5_exists(root, "/input/atom_roles")) {
            return fixed_atoms;
        }
        traverse_string_dset<1>(root, "/input/atom_roles", [&](size_t i, const std::string& name) {
            if(name == atom_role_name) {
                fixed_atoms.push_back(static_cast<int>(i));
            }
        });
    } catch(...) {
        return std::vector<int>();
    }
    return fixed_atoms;
}

// Register fix rigid constraints for an engine
void register_fix_rigid_for_engine(hid_t config_root, DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto& masses = g_atom_masses[&engine];
    masses.assign(std::max(0, engine.pos->n_atom), 1.f);
    if(h5_exists(config_root, "/input/mass")) {
        traverse_dset<1,float>(config_root, "/input/mass", [&](size_t i, float m) {
            if(i < masses.size() && m > 0.f) masses[i] = m;
        });
    }
    auto fixed_atoms = read_fix_rigid_settings(config_root);
    merge_fixed_atoms(engine, fixed_atoms);
    rebuild_rigid_groups(engine);
}

void register_fix_rigid_backbone_for_engine(hid_t config_root, DerivEngine& engine, const std::string& atom_name) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    if(h5_exists(config_root, "/input/potential/martini_potential")) {
        if(!h5_exists(config_root, "/input/atom_roles")) {
            throw string("MARTINI backbone hold requires /input/atom_roles");
        }
    }
    auto fixed_atoms = read_martini_backbone_hold(config_root, atom_name);
    merge_fixed_atoms(engine, fixed_atoms);
}

void set_dynamic_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto& dyn = g_dynamic_fixed_atoms[&engine];
    dyn = atom_indices;
    normalize_atom_list(dyn);
    rebuild_fixed_atoms(engine);
}

void clear_dynamic_fixed_atoms(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    g_dynamic_fixed_atoms.erase(&engine);
    rebuild_fixed_atoms(engine);
}

void set_dynamic_z_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto& dyn = g_dynamic_z_fixed_atoms[&engine];
    dyn = atom_indices;
    normalize_atom_list(dyn);
    rebuild_fixed_atoms(engine);
}

void clear_dynamic_z_fixed_atoms(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    g_dynamic_z_fixed_atoms.erase(&engine);
    rebuild_fixed_atoms(engine);
}

static void rebuild_rigid_groups(DerivEngine& engine) {
    auto it = g_dynamic_rigid_groups.find(&engine);
    if(it == g_dynamic_rigid_groups.end()) {
        g_rigid_groups.erase(&engine);
        return;
    }
    const auto mit = g_atom_masses.find(&engine);
    const std::vector<float> empty_masses;
    const auto& masses = (mit != g_atom_masses.end()) ? mit->second : empty_masses;
    std::vector<RigidBodyGroup> groups;
    for(const auto& raw_group : it->second) {
        std::vector<int> atoms = raw_group;
        normalize_atom_list(atoms);
        if(atoms.size() < 2) continue;
        RigidBodyGroup grp;
        grp.atom_indices = atoms;
        grp.masses.resize(atoms.size(), 1.f);
        for(size_t i = 0; i < atoms.size(); ++i) {
            int atom = atoms[i];
            if(atom >= 0 && atom < static_cast<int>(masses.size()) && masses[atom] > 0.f) {
                grp.masses[i] = masses[atom];
            }
        }
        groups.push_back(std::move(grp));
    }
    if(groups.empty()) g_rigid_groups.erase(&engine);
    else g_rigid_groups[&engine] = std::move(groups);
}

void set_dynamic_rigid_groups(DerivEngine& engine, const std::vector<std::vector<int>>& groups) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    if(groups.empty()) {
        g_dynamic_rigid_groups.erase(&engine);
        g_rigid_groups.erase(&engine);
        return;
    }
    g_dynamic_rigid_groups[&engine] = groups;
    rebuild_rigid_groups(engine);
}

void clear_dynamic_rigid_groups(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    g_dynamic_rigid_groups.erase(&engine);
    g_rigid_groups.erase(&engine);
}

// Apply fix rigid constraints during minimization
void apply_fix_rigid_minimization(DerivEngine& engine, VecArray pos, VecArray deriv) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_fixed_atoms.find(&engine);
    if(it != g_fixed_atoms.end()) {
        const auto& fixed_atoms = it->second;
        for(int atom_idx : fixed_atoms) {
            if(atom_idx >= 0 && atom_idx < engine.pos->n_atom) {
                // Zero out forces on fixed atoms
                deriv(0, atom_idx) = 0.0f;
                deriv(1, atom_idx) = 0.0f;
                deriv(2, atom_idx) = 0.0f;
            }
        }
    }
    auto zit = g_z_fixed_atoms.find(&engine);
    if(zit != g_z_fixed_atoms.end()) {
        const auto& z_fixed_atoms = zit->second;
        for(int atom_idx : z_fixed_atoms) {
            if(atom_idx >= 0 && atom_idx < engine.pos->n_atom) {
                deriv(2, atom_idx) = 0.0f;
            }
        }
    }
    auto rit = g_rigid_groups.find(&engine);
    if(rit != g_rigid_groups.end()) {
        for(auto& group : rit->second) {
            apply_rigid_group_minimization(engine, group, pos, deriv);
        }
    }
}

void apply_fix_rigid_projection(DerivEngine& engine, VecArray pos) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto rit = g_rigid_groups.find(&engine);
    if(rit == g_rigid_groups.end()) return;
    for(auto& group : rit->second) {
        project_rigid_group_geometry(engine, group, pos);
    }
}

// Apply fix rigid constraints during MD (zero forces and velocities)
void apply_fix_rigid_md(DerivEngine& engine, VecArray pos, VecArray deriv, VecArray mom) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_fixed_atoms.find(&engine);
    if(it != g_fixed_atoms.end()) {
        const auto& fixed_atoms = it->second;
        for(int atom_idx : fixed_atoms) {
            if(atom_idx >= 0 && atom_idx < engine.pos->n_atom) {
                // Zero out forces on fixed atoms
                deriv(0, atom_idx) = 0.0f;
                deriv(1, atom_idx) = 0.0f;
                deriv(2, atom_idx) = 0.0f;
                
                // Zero out velocities on fixed atoms
                if(mom.row_width > 0) {
                    mom(0, atom_idx) = 0.0f;
                    mom(1, atom_idx) = 0.0f;
                    mom(2, atom_idx) = 0.0f;
                }
            }
        }
    }
    auto zit = g_z_fixed_atoms.find(&engine);
    if(zit != g_z_fixed_atoms.end()) {
        const auto& z_fixed_atoms = zit->second;
        for(int atom_idx : z_fixed_atoms) {
            if(atom_idx >= 0 && atom_idx < engine.pos->n_atom) {
                deriv(2, atom_idx) = 0.0f;
                if(mom.row_width > 2) {
                    mom(2, atom_idx) = 0.0f;
                }
            }
        }
    }
    auto rit = g_rigid_groups.find(&engine);
    if(rit != g_rigid_groups.end()) {
        for(auto& group : rit->second) {
            apply_rigid_group_constraints(engine, group, pos, deriv, mom);
        }
    }
}

// Check if an atom is fixed
bool is_atom_fixed(const DerivEngine& engine, int atom_idx) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_fixed_atoms.find(const_cast<DerivEngine*>(&engine));
    if(it != g_fixed_atoms.end()) {
        const auto& fixed_atoms = it->second;
        return std::find(fixed_atoms.begin(), fixed_atoms.end(), atom_idx) != fixed_atoms.end();
    }
    return false;
}

// Get list of fixed atoms for an engine
std::vector<int> get_fixed_atoms(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_fixed_atoms.find(const_cast<DerivEngine*>(&engine));
    if(it != g_fixed_atoms.end()) {
        return it->second;
    }
    return std::vector<int>();
}

std::vector<int> get_z_fixed_atoms(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_z_fixed_atoms.find(const_cast<DerivEngine*>(&engine));
    if(it != g_z_fixed_atoms.end()) {
        return it->second;
    }
    return std::vector<int>();
}

} // namespace martini_fix_rigid
