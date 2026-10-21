#include "martini.h"
#include "deriv_engine.h"
#include <mutex>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

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

struct RigidBodyGroup {
    std::vector<int> atom_indices;
    std::vector<float> masses;
    std::vector<std::array<float,3>> reference_rel;
    bool initialized = false;
};
static std::map<DerivEngine*, std::vector<std::vector<int>>> g_dynamic_rigid_groups;
static std::map<DerivEngine*, std::vector<RigidBodyGroup>> g_rigid_groups;

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
    for(auto* atoms: {&merged, &merged_z})
        if(!atoms->empty() && (atoms->front() < 0 || atoms->back() >= engine.pos->n_atom))
            throw std::string("fixed atom index out of range");
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

static void invert_3x3(const float a[3][3], float inv_out[3][3]) {
    float det =
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if(std::fabs(det) < 1.0e-12f) throw string("rigid group has a singular inertia tensor (collinear atoms)");
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

// Rotation R taking the reference onto the current positions, R*reference_rel[i] ~ current_rel[i]: the
// quaternion is the eigenvector of the largest eigenvalue of Horn's matrix built from
// S(a,b) = sum_i m_i reference_a current_b, as in AffineAlignment (eig.cpp).
static void best_fit_rotation_from_reference(
        const std::vector<std::array<float,3>>& reference_rel,
        const std::vector<std::array<float,3>>& current_rel,
        const std::vector<float>& masses,
        float R[3][3]) {
    Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
    for(size_t i = 0; i < reference_rel.size(); ++i)
        for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
                S(a,b) += double(masses[i]) * reference_rel[i][a] * current_rel[i][b];

    Eigen::Matrix4d N;
    N << S(0,0)+S(1,1)+S(2,2), S(1,2)-S(2,1),         S(2,0)-S(0,2),         S(0,1)-S(1,0),
         S(1,2)-S(2,1),        S(0,0)-S(1,1)-S(2,2),  S(0,1)+S(1,0),         S(0,2)+S(2,0),
         S(2,0)-S(0,2),        S(0,1)+S(1,0),        -S(0,0)+S(1,1)-S(2,2),  S(1,2)+S(2,1),
         S(0,1)-S(1,0),        S(0,2)+S(2,0),         S(1,2)+S(2,1),        -S(0,0)-S(1,1)+S(2,2);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eig(N);
    Eigen::Vector4d q = eig.eigenvectors().col(3);   // eigenvalues are in increasing order
    float qf[4] = {float(q[0]), float(q[1]), float(q[2]), float(q[3])};
    quat_to_rotmat(qf, R);
}

static std::array<float,3> center_of_mass(const RigidBodyGroup& group, VecArray pos) {
    float total_mass = 0.f;
    std::array<float,3> com{{0.f, 0.f, 0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        float m = group.masses[i];
        total_mass += m;
        com[0] += m * pos(0, atom);
        com[1] += m * pos(1, atom);
        com[2] += m * pos(2, atom);
    }
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;
    return com;
}

// Moves the group onto the rigid placement of its reference that best fits the current positions, and returns
// the atoms' positions relative to the centre of mass.  The reference is the group's geometry when it was built.
static std::vector<std::array<float,3>> project_rigid_group_geometry(RigidBodyGroup& group, VecArray pos) {
    std::vector<std::array<float,3>> rel(group.atom_indices.size());
    auto com = center_of_mass(group, pos);
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        for(int d = 0; d < 3; ++d) rel[i][d] = pos(d, atom) - com[d];
    }
    if(!group.initialized) {
        group.reference_rel = rel;
        group.initialized = true;
        return rel;
    }

    float Rfit[3][3];
    best_fit_rotation_from_reference(group.reference_rel, rel, group.masses, Rfit);
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        rel[i] = apply_rot3(Rfit, group.reference_rel[i]);
        for(int d = 0; d < 3; ++d) pos(d, atom) = com[d] + rel[i][d];
    }
    return rel;
}

static void add_inertia(float inertia[3][3], float m, const std::array<float,3>& r) {
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

// Places the group rigidly and keeps only the rigid-body part of its momenta: the centre-of-mass velocity and
// the angular velocity I^-1 L about the centre of mass.
static void apply_rigid_group_constraints(RigidBodyGroup& group, VecArray pos, VecArray mom) {
    auto rel = project_rigid_group_geometry(group, pos);

    float total_mass = 0.f;
    std::array<float,3> total_momentum{{0.f, 0.f, 0.f}};
    std::array<float,3> ang_momentum{{0.f, 0.f, 0.f}};
    float inertia[3][3] = {{0.f,0.f,0.f},{0.f,0.f,0.f},{0.f,0.f,0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        float m = group.masses[i];
        std::array<float,3> p{mom(0, atom), mom(1, atom), mom(2, atom)};
        total_mass += m;
        total_momentum[0] += p[0];
        total_momentum[1] += p[1];
        total_momentum[2] += p[2];
        auto rx_p = cross3(rel[i], p);
        ang_momentum[0] += rx_p[0];
        ang_momentum[1] += rx_p[1];
        ang_momentum[2] += rx_p[2];
        add_inertia(inertia, m, rel[i]);
    }

    float inertia_inv[3][3];
    invert_3x3(inertia, inertia_inv);
    auto omega = mat3_mul_vec(inertia_inv, ang_momentum);
    std::array<float,3> vcm{
        total_momentum[0] / total_mass,
        total_momentum[1] / total_mass,
        total_momentum[2] / total_mass,
    };

    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        float m = group.masses[i];
        auto omega_x_r = cross3(omega, rel[i]);
        mom(0, atom) = m * (vcm[0] + omega_x_r[0]);
        mom(1, atom) = m * (vcm[1] + omega_x_r[1]);
        mom(2, atom) = m * (vcm[2] + omega_x_r[2]);
    }
}

// Places the group rigidly and replaces its gradient by the rigid-body motion the net force and torque drive.
static void apply_rigid_group_minimization(RigidBodyGroup& group, VecArray pos, VecArray deriv) {
    auto rel = project_rigid_group_geometry(group, pos);

    float total_mass = 0.f;
    std::array<float,3> total_force{{0.f, 0.f, 0.f}};
    std::array<float,3> torque{{0.f, 0.f, 0.f}};
    float inertia[3][3] = {{0.f,0.f,0.f},{0.f,0.f,0.f},{0.f,0.f,0.f}};
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        float m = group.masses[i];
        std::array<float,3> f{-deriv(0, atom), -deriv(1, atom), -deriv(2, atom)};
        total_mass += m;
        total_force[0] += f[0];
        total_force[1] += f[1];
        total_force[2] += f[2];
        auto rx_f = cross3(rel[i], f);
        torque[0] += rx_f[0];
        torque[1] += rx_f[1];
        torque[2] += rx_f[2];
        add_inertia(inertia, m, rel[i]);
    }

    float inertia_inv[3][3];
    invert_3x3(inertia, inertia_inv);
    auto alpha = mat3_mul_vec(inertia_inv, torque);
    std::array<float,3> acom{
        total_force[0] / total_mass,
        total_force[1] / total_mass,
        total_force[2] / total_mass,
    };
    for(size_t i = 0; i < group.atom_indices.size(); ++i) {
        int atom = group.atom_indices[i];
        float m = group.masses[i];
        auto alpha_x_r = cross3(alpha, rel[i]);
        deriv(0, atom) = -m * (acom[0] + alpha_x_r[0]);
        deriv(1, atom) = -m * (acom[1] + alpha_x_r[1]);
        deriv(2, atom) = -m * (acom[2] + alpha_x_r[2]);
    }
}

void unregister_fix_rigid_for_engine(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    for(auto* m: {&g_user_fixed_atoms, &g_dynamic_fixed_atoms, &g_dynamic_z_fixed_atoms, &g_fixed_atoms, &g_z_fixed_atoms})
        m->erase(&engine);
    g_dynamic_rigid_groups.erase(&engine);
    g_rigid_groups.erase(&engine);
}

void register_fix_rigid_for_engine(hid_t config_root, DerivEngine& engine) {
    unregister_fix_rigid_for_engine(engine);
    if(!h5_exists(config_root, "/input/fix_rigid")) return;

    std::vector<int> fixed_atoms;
    int n_atom = engine.pos->n_atom;
    traverse_dset<1,int>(config_root, "/input/fix_rigid/atom_indices", [&](size_t i, int atom_idx) {
        if(atom_idx < 0 || atom_idx >= n_atom) throw string("/input/fix_rigid/atom_indices out of range");
        fixed_atoms.push_back(atom_idx);});
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    merge_fixed_atoms(engine, fixed_atoms);
}

void set_dynamic_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto& dyn = g_dynamic_fixed_atoms[&engine];
    dyn = atom_indices;
    normalize_atom_list(dyn);
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
    std::vector<RigidBodyGroup> groups;
    for(const auto& raw_group : it->second) {
        std::vector<int> atoms = raw_group;
        normalize_atom_list(atoms);
        if(atoms.size() < 2) continue;
        RigidBodyGroup grp;
        grp.atom_indices = atoms;
        for(int atom : atoms) {
            if(atom < 0 || atom >= engine.pos->n_atom) throw string("rigid group atom index out of range");
            grp.masses.push_back(martini_masses::get_mass(&engine, atom));
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
            apply_rigid_group_minimization(group, pos, deriv);
        }
    }
}

void apply_fix_rigid_projection(DerivEngine& engine, VecArray pos) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto rit = g_rigid_groups.find(&engine);
    if(rit == g_rigid_groups.end()) return;
    for(auto& group : rit->second) {
        project_rigid_group_geometry(group, pos);
    }
}

// Apply fix rigid constraints during MD (zero forces and velocities)
// Holds fixed atoms at rest, z-fixed atoms at rest along z, and each rigid group to rigid-body motion.
void apply_fix_rigid_md(DerivEngine& engine, VecArray pos, VecArray mom) {
    std::lock_guard<std::mutex> lock(g_fix_rigid_mutex);
    auto it = g_fixed_atoms.find(&engine);
    if(it != g_fixed_atoms.end())
        for(int atom_idx : it->second)
            for(int d = 0; d < 3; ++d) mom(d, atom_idx) = 0.f;
    auto zit = g_z_fixed_atoms.find(&engine);
    if(zit != g_z_fixed_atoms.end())
        for(int atom_idx : zit->second) mom(2, atom_idx) = 0.f;
    auto rit = g_rigid_groups.find(&engine);
    if(rit != g_rigid_groups.end())
        for(auto& group : rit->second) apply_rigid_group_constraints(group, pos, mom);
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
