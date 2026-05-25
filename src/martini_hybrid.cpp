#include "martini_internal.h"
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
#include <limits>
#include "box.h" // For PBC minimum_image function

using namespace h5;
using namespace std;

namespace martini_hybrid {
static std::mutex g_hybrid_mutex;

static std::map<DerivEngine*, std::shared_ptr<HybridRuntimeState>> g_hybrid_state;
static std::map<const CoordNode*, std::shared_ptr<HybridRuntimeState>> g_hybrid_state_by_coord;

static std::vector<int> active_protein_proxy_fixed_atoms(const HybridRuntimeState& st) {
    std::vector<int> atoms;
    if(!st.enabled || !st.active) return atoms;
    atoms.reserve(st.protein_membership.size());
    for(size_t atom_idx = 0; atom_idx < st.protein_membership.size(); ++atom_idx) {
        if(st.protein_membership[atom_idx] < 0) continue;
        if(atom_idx >= st.atom_role_class.size()) continue;
        auto role = st.atom_role_class[atom_idx];
        bool is_backbone_carrier =
            atom_idx < st.atom_backbone_carrier_mask.size() &&
            st.atom_backbone_carrier_mask[atom_idx] != 0u;
        if((role == ROLE_BB || role == ROLE_SC) && !is_backbone_carrier) {
            atoms.push_back(static_cast<int>(atom_idx));
        }
    }
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    return atoms;
}

static inline int direct_ca_atom_for_bb_proxy(const HybridRuntimeState& st, int bb_proxy_atom) {
    if(bb_proxy_atom < 0 || bb_proxy_atom >= static_cast<int>(st.bb_proxy_to_ca_atom.size())) return -1;
    return st.bb_proxy_to_ca_atom[static_cast<size_t>(bb_proxy_atom)];
}

int bb_map_index_for_proxy(const HybridRuntimeState& st, int bb_proxy_atom) {
    if(bb_proxy_atom < 0 || bb_proxy_atom >= static_cast<int>(st.bb_proxy_to_map_index.size())) return -1;
    return st.bb_proxy_to_map_index[static_cast<size_t>(bb_proxy_atom)];
}

void project_bb_proxy_gradient_if_active(
        const HybridRuntimeState& st,
        VecArray pos_sens,
        int n_atom,
        int bb_proxy_atom,
        const Vec<3>& grad) {
    int map_idx = bb_map_index_for_proxy(st, bb_proxy_atom);
    if(map_idx < 0 || map_idx >= static_cast<int>(st.n_bb)) return;
    bool site_is_carrier =
        bb_proxy_atom >= 0 &&
        bb_proxy_atom < static_cast<int>(st.atom_backbone_carrier_mask.size()) &&
        st.atom_backbone_carrier_mask[static_cast<size_t>(bb_proxy_atom)] != 0u;

    for(int d = 0; d < 4; ++d) {
        if(st.atom_mask[static_cast<size_t>(map_idx)][d] == 0) continue;
        int atom_idx = st.atom_indices[static_cast<size_t>(map_idx)][d];
        float w = st.weights[static_cast<size_t>(map_idx)][d];
        if(atom_idx < 0 || atom_idx >= n_atom || w == 0.f) continue;
        if(!site_is_carrier && atom_idx == bb_proxy_atom) continue;
        update_vec<3>(pos_sens, atom_idx, w * grad);
    }
}

static inline bool mapped_bb_site_position_if_active(
        const HybridRuntimeState& st,
        VecArray pos,
        int n_atom,
        int bb_site_atom,
        Vec<3>& site_pos) {
    int map_idx = bb_map_index_for_proxy(st, bb_site_atom);
    if(map_idx < 0 || map_idx >= static_cast<int>(st.n_bb)) return false;
    Vec<3> com = make_zero<3>();
    float wsum = 0.f;
    for(int d = 0; d < 4; ++d) {
        if(st.atom_mask[static_cast<size_t>(map_idx)][d] == 0) continue;
        int atom_idx = st.atom_indices[static_cast<size_t>(map_idx)][d];
        float w = st.weights[static_cast<size_t>(map_idx)][d];
        if(atom_idx < 0 || atom_idx >= n_atom || w == 0.f) continue;
        com += w * load_vec<3>(pos, atom_idx);
        wsum += w;
    }
    if(!(wsum > 0.f)) return false;
    if(fabsf(wsum - 1.f) > 1e-6f) com *= (1.f / wsum);
    site_pos = com;
    return true;
}

static void refresh_backbone_o_positions_if_active(const HybridRuntimeState& st, VecArray pos, int n_atom);

void refresh_bb_positions_if_active(
        const HybridRuntimeState& st,
        VecArray pos,
        int n_atom) {
    if(!st.enabled || !st.active) return;
    refresh_backbone_o_positions_if_active(st, pos, n_atom);
    for(size_t k = 0; k < st.n_bb; ++k) {
        int bb = st.bb_atom_index[k];
        if(bb < 0 || bb >= n_atom) continue;
        Vec<3> com = make_zero<3>();
        float wsum = 0.f;
        for(int d = 0; d < 4; ++d) {
            if(st.atom_mask[k][d] == 0) continue;
            int ai = st.atom_indices[k][d];
            float w = st.weights[k][d];
            if(ai < 0 || ai >= n_atom || w == 0.f) continue;
            if(ai == bb) continue;
            com += w * load_vec<3>(pos, ai);
            wsum += w;
        }
        if(wsum > 0.f) {
            if(fabsf(wsum - 1.f) > 1e-6f) com *= (1.f / wsum);
            store_vec<3>(pos, bb, com);
        }
    }
}

static inline std::array<float,3> apply_rot(const float R[3][3], const std::array<float,3>& v) {
    return std::array<float,3>{
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
    };
}

static inline std::array<float,3> vec_sub(const std::array<float,3>& a, const std::array<float,3>& b) {
    return std::array<float,3>{a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

static inline std::array<float,3> vec_add(const std::array<float,3>& a, const std::array<float,3>& b) {
    return std::array<float,3>{a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

static inline std::array<float,3> vec_scale(const std::array<float,3>& a, float s) {
    return std::array<float,3>{a[0]*s, a[1]*s, a[2]*s};
}

static inline float vec_dot(const std::array<float,3>& a, const std::array<float,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline std::array<float,3> vec_cross(const std::array<float,3>& a, const std::array<float,3>& b) {
    return std::array<float,3>{
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

static inline float vec_norm(const std::array<float,3>& a) {
    return sqrtf(vec_dot(a, a));
}

static inline std::array<float,3> vec_normalize(const std::array<float,3>& a) {
    float n = vec_norm(a);
    if(n <= 1e-8f) return std::array<float,3>{0.f,0.f,0.f};
    return vec_scale(a, 1.f/n);
}

static inline bool build_frame_from_three(
        const std::array<float,3>& p0,
        const std::array<float,3>& p1,
        const std::array<float,3>& p2,
        float F[3][3]) {
    auto e1 = vec_normalize(vec_sub(p1, p0));
    auto v2 = vec_sub(p2, p0);
    auto e3 = vec_cross(e1, v2);
    float n3 = vec_norm(e3);
    if(n3 <= 1e-8f) return false;
    e3 = vec_scale(e3, 1.f/n3);
    auto e2 = vec_cross(e3, e1);

    // Columns are basis vectors.
    F[0][0] = e1[0]; F[1][0] = e1[1]; F[2][0] = e1[2];
    F[0][1] = e2[0]; F[1][1] = e2[1]; F[2][1] = e2[2];
    F[0][2] = e3[0]; F[1][2] = e3[1]; F[2][2] = e3[2];
    return true;
}

static inline void mat_mul(const float A[3][3], const float B[3][3], float C[3][3]) {
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
        C[i][j] = 0.f;
        for(int k=0;k<3;++k) C[i][j] += A[i][k]*B[k][j];
    }
}

static inline void mat_transpose(const float A[3][3], float AT[3][3]) {
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) AT[i][j] = A[j][i];
}

static void refresh_backbone_o_positions_if_active(const HybridRuntimeState& st, VecArray pos, int n_atom) {
    if(!st.enabled || !st.active) return;
    if(st.bb_reference_atom_coords.size() != st.n_bb) return;
    if(st.bb_reference_runtime_atom_indices.size() != st.n_bb) return;

    for(size_t k = 0; k < st.n_bb; ++k) {
        const auto& atom_idx = st.bb_reference_runtime_atom_indices[k];
        const int n_idx = atom_idx[0];
        const int ca_idx = atom_idx[1];
        const int c_idx = atom_idx[2];
        const int o_idx = atom_idx[3];
        if(n_idx < 0 || n_idx >= n_atom ||
           ca_idx < 0 || ca_idx >= n_atom ||
           c_idx < 0 || c_idx >= n_atom ||
           o_idx < 0 || o_idx >= n_atom) {
            continue;
        }

        const auto& ref = st.bb_reference_atom_coords[k];
        float F_ref[3][3], F_cur[3][3], F_ref_T[3][3], R[3][3];
        if(!build_frame_from_three(ref[0], ref[1], ref[2], F_ref)) continue;

        auto cur_n = load_vec<3>(pos, n_idx);
        auto cur_ca = load_vec<3>(pos, ca_idx);
        auto cur_c = load_vec<3>(pos, c_idx);
        auto cur_n_arr = std::array<float,3>{cur_n[0], cur_n[1], cur_n[2]};
        auto cur_ca_arr = std::array<float,3>{cur_ca[0], cur_ca[1], cur_ca[2]};
        auto cur_c_arr = std::array<float,3>{cur_c[0], cur_c[1], cur_c[2]};
        if(!build_frame_from_three(cur_n_arr, cur_ca_arr, cur_c_arr, F_cur)) continue;

        mat_transpose(F_ref, F_ref_T);
        mat_mul(F_cur, F_ref_T, R);

        auto ref_local_o = vec_sub(ref[3], ref[1]);
        auto mapped_o = vec_add(apply_rot(R, ref_local_o), cur_ca_arr);
        store_vec<3>(pos, o_idx, make_vec3(mapped_o[0], mapped_o[1], mapped_o[2]));
    }
}
static std::string trim_h5_string(const std::string& in);

static std::vector<std::string> split_csv_tokens(const std::string& s) {
    auto trim_token = [](const std::string& in) {
        size_t begin = 0;
        while(begin < in.size() && (in[begin] == '\0' || std::isspace(static_cast<unsigned char>(in[begin])))) {
            ++begin;
        }
        size_t end = in.size();
        while(end > begin && (in[end - 1] == '\0' || std::isspace(static_cast<unsigned char>(in[end - 1])))) {
            --end;
        }
        return in.substr(begin, end - begin);
    };
    std::vector<std::string> out;
    std::string cur;
    for(char c : s) {
        if(c == ',') {
            auto t = trim_token(cur);
            if(!t.empty()) out.push_back(t);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    auto t = trim_token(cur);
    if(!t.empty()) out.push_back(t);
    return out;
}

static std::string normalize_role_token(const std::string& raw) {
    std::string out = trim_h5_string(raw);
    for(char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

static std::string normalize_mode_token(const std::string& raw) {
    std::string out = trim_h5_string(raw);
    for(char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

static inline unsigned char classify_atom_role_name(const std::string& raw_name) {
    size_t begin = 0;
    while(begin < raw_name.size() &&
          (raw_name[begin] == '\0' || std::isspace(static_cast<unsigned char>(raw_name[begin])))) {
        ++begin;
    }
    size_t end = raw_name.size();
    while(end > begin &&
          (raw_name[end - 1] == '\0' || std::isspace(static_cast<unsigned char>(raw_name[end - 1])))) {
        --end;
    }
    std::string name = raw_name.substr(begin, end - begin);
    for(char& c : name) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if(name == "BB") return ROLE_BB;
    if(name.size() >= 2 && name[0] == 'S' && name[1] == 'C') return ROLE_SC;
    return ROLE_OTHER;
}

static inline bool is_backbone_carrier_role_name(const std::string& raw_name) {
    std::string name = normalize_role_token(raw_name);
    return name == "N" || name == "CA" || name == "C" || name == "O";
}

static inline bool is_env_po4_atom(const HybridRuntimeState& st, int atom) {
    return atom >= 0 &&
           atom < static_cast<int>(st.sc_env_po4_env_mask.size()) &&
           st.sc_env_po4_env_mask[atom] != 0;
}

static inline bool active_sc_env_backbone_hold_enabled(const HybridRuntimeState& st) {
    return st.enabled &&
           st.active &&
           st.sc_env_backbone_hold_steps > 0 &&
           st.sc_env_transition_step < static_cast<uint64_t>(st.sc_env_backbone_hold_steps);
}

float compute_sc_backbone_feedback_mix(const HybridRuntimeState& st) {
    if(!active_sc_env_backbone_hold_enabled(st)) return 1.f;
    if(st.sc_env_backbone_hold_steps <= 1) return 1.f;

    uint64_t final_hold_step = static_cast<uint64_t>(st.sc_env_backbone_hold_steps - 1);
    if(st.sc_env_transition_step >= final_hold_step) return 1.f;
    return float(st.sc_env_transition_step) / float(final_hold_step);
}

static inline bool active_sc_env_force_ramp_enabled(const HybridRuntimeState& st) {
    return st.enabled &&
           st.active &&
           st.sc_env_relax_steps > 0 &&
           st.sc_env_transition_step < static_cast<uint64_t>(st.sc_env_relax_steps);
}

static inline bool active_sc_env_po4_z_hold_enabled(const HybridRuntimeState& st) {
    return st.enabled &&
           st.active &&
           st.sc_env_po4_z_clamp_enabled &&
           st.sc_env_po4_z_hold_steps > 0 &&
           st.sc_env_transition_step < static_cast<uint64_t>(st.sc_env_po4_z_hold_steps);
}

static inline bool active_cross_interface_pair(
        const HybridRuntimeState& st,
        bool i_is_protein,
        bool j_is_protein) {
    return st.enabled && st.active && (i_is_protein != j_is_protein);
}

float active_interface_interaction_scale(
        const HybridRuntimeState& st,
        bool i_is_protein,
        bool j_is_protein) {
    return active_cross_interface_pair(st, i_is_protein, j_is_protein)
               ? st.protein_env_interface_scale
               : 1.f;
}

bool deterministic_startup_pair_cap_enabled(
        const HybridRuntimeState& st,
        bool i_is_protein,
        bool j_is_protein,
        unsigned char i_role,
        unsigned char j_role) {
    if(!active_sc_env_force_ramp_enabled(st)) return false;
    if(i_is_protein == j_is_protein) return false;
    unsigned char protein_role = i_is_protein ? i_role : j_role;
    return protein_role == ROLE_BB || protein_role == ROLE_SC;
}

void initialize_sc_env_po4_z_reference(
        HybridRuntimeState& st,
        VecArray pos,
        int n_atom) {
    if(!active_sc_env_po4_z_hold_enabled(st)) return;
    if(st.sc_env_po4_z_reference_initialized &&
       static_cast<int>(st.sc_env_po4_z_reference.size()) == n_atom) {
        return;
    }

    st.sc_env_po4_z_reference.assign(static_cast<size_t>(n_atom), 0.f);
    for(int i = 0; i < n_atom; ++i) {
        if(!is_env_po4_atom(st, i)) continue;
        st.sc_env_po4_z_reference[i] = pos(2, i);
    }
    st.sc_env_po4_z_reference_initialized = true;
}

static inline bool same_protein_membership_pair(const HybridRuntimeState& st, int i, int j) {
    if(i < 0 || j < 0) return false;
    if(i >= (int)st.protein_membership.size() || j >= (int)st.protein_membership.size()) return false;
    int pi = st.protein_membership[i];
    int pj = st.protein_membership[j];
    return (pi >= 0 && pj >= 0 && pi == pj);
}

unsigned char atom_role_class_at(const HybridRuntimeState& st, int i) {
    if(i < 0 || i >= (int)st.atom_role_class.size()) return ROLE_OTHER;
    return st.atom_role_class[i];
}

bool atom_is_backbone_carrier_at(const HybridRuntimeState& st, int i) {
    return i >= 0 &&
           i < static_cast<int>(st.atom_backbone_carrier_mask.size()) &&
           st.atom_backbone_carrier_mask[static_cast<size_t>(i)] != 0;
}

static inline bool allow_protein_pair_by_rule(const HybridRuntimeState& st, int i, int j) {
    (void)st;
    (void)i;
    (void)j;
    // In the direct-Upside active stage, protein internal MARTINI proxy terms
    // are bookkeeping only. Keeping any proxy-proxy protein interaction live
    // leaks legacy bonded/nonbonded energy into the production bucket without
    // feeding that force back through the Upside carrier path.
    return false;
}

bool allow_intra_protein_pair_if_active(const HybridRuntimeState& st, int i, int j) {
    if(!st.enabled || !st.active) return true;
    if(!same_protein_membership_pair(st, i, j)) return true;
    return allow_protein_pair_by_rule(st, i, j);
}

bool allow_multibody_term_if_active(const HybridRuntimeState& st, const index_t* atoms, int n_atom_dep) {
    if(!st.enabled || !st.active) return true;
    for(int a = 0; a < n_atom_dep; ++a) {
        for(int b = a + 1; b < n_atom_dep; ++b) {
            int ia = atoms[a];
            int ib = atoms[b];
            if(!same_protein_membership_pair(st, ia, ib)) continue;
            if(!allow_protein_pair_by_rule(st, ia, ib)) return false;
        }
    }
    return true;
}

static inline bool attribute_exists_hybrid(hid_t loc_id, const char* obj_name, const char* attr_name) {
    // H5Oopen(loc, ".") may alias the same underlying ID in some HDF5 builds.
    // Closing that alias can invalidate the caller's group handle and cause
    // downstream "error -1" / close-reference failures.
    if(!obj_name || obj_name[0] == '\0' || (obj_name[0] == '.' && obj_name[1] == '\0')) {
        htri_t exists = H5Aexists(loc_id, attr_name);
        return exists > 0;
    }

    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id < 0) return false;
    htri_t exists = H5Aexists(obj_id, attr_name);
    H5Oclose(obj_id);
    return exists > 0;
}

static std::string trim_h5_string(const std::string& in) {
    size_t begin = 0;
    while(begin < in.size() && (in[begin] == '\0' || std::isspace(static_cast<unsigned char>(in[begin])))) {
        ++begin;
    }
    size_t end = in.size();
    while(end > begin && (in[end - 1] == '\0' || std::isspace(static_cast<unsigned char>(in[end - 1])))) {
        --end;
    }
    return in.substr(begin, end - begin);
}

std::string read_string_attribute_or_default(hid_t group, const char* attr_name, const std::string& fallback) {
    if(!attribute_exists_hybrid(group, ".", attr_name)) return fallback;

    auto attr = h5_obj(H5Aclose, H5Aopen(group, attr_name, H5P_DEFAULT));
    auto dtype = h5_obj(H5Tclose, H5Aget_type(attr.get()));
    if(H5Tis_variable_str(dtype.get())) {
        char* raw = nullptr;
        if(H5Aread(attr.get(), dtype.get(), &raw) < 0) return fallback;
        std::string out = raw ? trim_h5_string(std::string(raw)) : fallback;
        if(raw) H5free_memory(raw);
        return out.empty() ? fallback : out;
    }
    size_t nchar = H5Tget_size(dtype.get());
    if(nchar == 0) return fallback;

    std::vector<char> data(nchar + 1, '\0');
    h5_noerr(H5Aread(attr.get(), dtype.get(), data.data()));
    auto out = trim_h5_string(std::string(data.data(), nchar));
    return out.empty() ? fallback : out;
}

static HybridRuntimeState read_hybrid_settings(hid_t root, int n_atom) {
    HybridRuntimeState out;

    if(!h5_exists(root, "/input/hybrid_control")) {
        return out;
    }
    out.has_config = true;

    auto ctrl = open_group(root, "/input/hybrid_control");
    out.enabled = (read_attribute<int>(ctrl.get(), ".", "enable", 0) != 0);
    out.activation_stage = read_string_attribute_or_default(ctrl.get(), "activation_stage", "production");
    out.preprod_mode = read_string_attribute_or_default(ctrl.get(), "preprod_protein_mode", "rigid");
    out.exclude_intra_protein_martini =
        (read_attribute<int>(ctrl.get(), ".", "exclude_intra_protein_martini", 1) != 0);
    out.preprod_rigid = (out.preprod_mode == "rigid" || out.preprod_mode == "rigid_body");
    out.protein_env_interface_scale =
        read_attribute<float>(ctrl.get(), ".", "protein_env_interface_scale", out.protein_env_interface_scale);
    out.sc_env_lj_force_cap = read_attribute<float>(ctrl.get(), ".", "sc_env_lj_force_cap", out.sc_env_lj_force_cap);
    out.sc_env_coul_force_cap = read_attribute<float>(ctrl.get(), ".", "sc_env_coul_force_cap", out.sc_env_coul_force_cap);
    out.sc_env_relax_steps = read_attribute<int>(ctrl.get(), ".", "sc_env_relax_steps", out.sc_env_relax_steps);
    out.sc_env_backbone_hold_steps =
        read_attribute<int>(ctrl.get(), ".", "sc_env_backbone_hold_steps", out.sc_env_backbone_hold_steps);
    out.sc_env_po4_z_hold_steps =
        read_attribute<int>(ctrl.get(), ".", "sc_env_po4_z_hold_steps", out.sc_env_po4_z_hold_steps);
    out.sc_env_po4_z_clamp_enabled =
        (read_attribute<int>(ctrl.get(), ".", "sc_env_po4_z_clamp_enabled", out.sc_env_po4_z_clamp_enabled ? 1 : 0) != 0);
    int transition_start = read_attribute<int>(ctrl.get(), ".", "sc_env_transition_step_start", 0);
    out.sc_env_transition_step_start = static_cast<uint64_t>(std::max(0, transition_start));
    out.sc_env_transition_step = out.sc_env_transition_step_start;
    if(out.sc_env_lj_force_cap < 0.f) out.sc_env_lj_force_cap = 0.f;
    if(out.sc_env_coul_force_cap < 0.f) out.sc_env_coul_force_cap = 0.f;
    if(out.sc_env_relax_steps < 1) out.sc_env_relax_steps = 1;
    if(out.sc_env_backbone_hold_steps < 0) out.sc_env_backbone_hold_steps = 0;
    if(out.sc_env_po4_z_hold_steps < 0) out.sc_env_po4_z_hold_steps = 0;
    if(!(out.protein_env_interface_scale > 0.f) || !std::isfinite(out.protein_env_interface_scale)) {
        throw string("hybrid_control/protein_env_interface_scale must be finite and > 0");
    }

    if(!out.enabled) {
        return out;
    }

    if(!h5_exists(root, "/input/hybrid_bb_map")) {
        throw string("Hybrid mode enabled but /input/hybrid_bb_map is missing");
    }
    {
        auto bb = open_group(root, "/input/hybrid_bb_map");
        auto atom_idx_shape = get_dset_size(2, bb.get(), "atom_indices");
        auto atom_mask_shape = get_dset_size(2, bb.get(), "atom_mask");
        auto weights_shape = get_dset_size(2, bb.get(), "weights");
        auto res_shape = get_dset_size(1, bb.get(), "bb_residue_index");
        auto pid_shape = get_dset_size(1, bb.get(), "protein_id");
        if(atom_idx_shape[1] != 4 || atom_mask_shape[1] != 4 || weights_shape[1] != 4) {
            throw string("Hybrid BB map must use (n_bb,4) layout for atom_indices/atom_mask/weights");
        }
        if(atom_idx_shape[0] != atom_mask_shape[0] ||
           atom_idx_shape[0] != weights_shape[0] ||
           atom_idx_shape[0] != res_shape[0] ||
           atom_idx_shape[0] != pid_shape[0]) {
            throw string("Hybrid BB map datasets have inconsistent n_bb sizes");
        }
        out.n_bb = atom_idx_shape[0];
        out.bb_residue_index.assign(out.n_bb, -1);
        out.bb_atom_index.assign(out.n_bb, -1);
        out.bb_ca_atom_index.assign(out.n_bb, -1);
        out.atom_indices.assign(out.n_bb, std::array<int,4>{{-1,-1,-1,-1}});
        out.atom_mask.assign(out.n_bb, std::array<int,4>{{0,0,0,0}});
        out.weights.assign(out.n_bb, std::array<float,4>{{0.f,0.f,0.f,0.f}});
        out.bb_reference_runtime_atom_indices.assign(out.n_bb, std::array<int,4>{{-1,-1,-1,-1}});
        out.bb_reference_atom_coords.assign(
            out.n_bb,
            std::array<std::array<float,3>,4>{{
                std::array<float,3>{{0.f,0.f,0.f}},
                std::array<float,3>{{0.f,0.f,0.f}},
                std::array<float,3>{{0.f,0.f,0.f}},
                std::array<float,3>{{0.f,0.f,0.f}},
            }});

        traverse_dset<1,int>(bb.get(), "bb_residue_index", [&](size_t i, int v) {
            out.bb_residue_index[i] = v;
        });
        traverse_dset<2,int>(bb.get(), "atom_indices", [&](size_t i, size_t j, int v) {
            out.atom_indices[i][j] = v;
        });
        traverse_dset<2,int>(bb.get(), "atom_mask", [&](size_t i, size_t j, int v) {
            out.atom_mask[i][j] = v;
        });
        traverse_dset<2,float>(bb.get(), "weights", [&](size_t i, size_t j, float v) {
            out.weights[i][j] = v;
        });
        out.bb_reference_runtime_atom_indices = out.atom_indices;
        if(h5_exists(bb.get(), "reference_atom_coords")) {
            auto ref_shape = get_dset_size(3, bb.get(), "reference_atom_coords");
            if(ref_shape[0] != out.n_bb || ref_shape[1] != 4 || ref_shape[2] != 3) {
                throw string("Hybrid BB reference_atom_coords must have shape (n_bb,4,3)");
            }
            traverse_dset<3,float>(bb.get(), "reference_atom_coords", [&](size_t i, size_t j, size_t d, float v) {
                out.bb_reference_atom_coords[i][j][d] = v;
            });
        }
        if(h5_exists(bb.get(), "bb_atom_index")) {
            check_size(bb.get(), "bb_atom_index", out.n_bb);
            traverse_dset<1,int>(bb.get(), "bb_atom_index", [&](size_t i, int v) {
                out.bb_atom_index[i] = v;
            });
        }
    }

    if(h5_exists(root, "/input/hybrid_env_topology")) {
        auto env = open_group(root, "/input/hybrid_env_topology");
        auto env_idx_shape = get_dset_size(1, env.get(), "env_atom_indices");
        auto member_shape = get_dset_size(1, env.get(), "protein_membership");
        if(static_cast<int>(member_shape[0]) != n_atom) {
            throw string("Hybrid env topology protein_membership length must match n_atom");
        }
        out.n_env = env_idx_shape[0];
        out.protein_membership.assign(member_shape[0], -1);
        traverse_dset<1,int>(env.get(), "protein_membership", [&](size_t i, int v) {
            out.protein_membership[i] = v;
        });
    } else {
        throw string("Hybrid mode enabled but /input/hybrid_env_topology is missing");
    }
    out.atom_role_class.assign(static_cast<size_t>(n_atom), ROLE_OTHER);
    out.atom_backbone_carrier_mask.assign(static_cast<size_t>(n_atom), 0u);
    if(h5_exists(root, "/input/atom_roles")) {
        traverse_string_dset<1>(root, "/input/atom_roles", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) {
                out.atom_role_class[i] = classify_atom_role_name(v);
                out.atom_backbone_carrier_mask[i] = is_backbone_carrier_role_name(v) ? 1u : 0u;
            }
        });
    } else if(h5_exists(root, "/input/atom_names")) {
        traverse_string_dset<1>(root, "/input/atom_names", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) {
                out.atom_role_class[i] = classify_atom_role_name(v);
                out.atom_backbone_carrier_mask[i] = is_backbone_carrier_role_name(v) ? 1u : 0u;
            }
        });
    }

    std::vector<std::string> atom_roles(n_atom, "");
    bool has_roles = false;
    if(h5_exists(root, "/input/atom_roles")) {
        has_roles = true;
        traverse_string_dset<1>(root, "/input/atom_roles", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) atom_roles[i] = normalize_role_token(v);
        });
    } else if(h5_exists(root, "/input/atom_names")) {
        has_roles = true;
        traverse_string_dset<1>(root, "/input/atom_names", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) atom_roles[i] = normalize_role_token(v);
        });
    }

    if(out.sc_env_po4_z_clamp_enabled && !out.protein_membership.empty()) {
        out.sc_env_po4_env_mask.assign(static_cast<size_t>(n_atom), 0u);
        if(has_roles) {
            for(int i = 0; i < n_atom; ++i) {
                bool is_protein = (i < (int)out.protein_membership.size() && out.protein_membership[i] >= 0);
                if(is_protein) continue;
                if(i >= (int)atom_roles.size()) continue;
                if(atom_roles[i] == "PO4") {
                    out.sc_env_po4_env_mask[i] = 1u;
                }
            }
        }
    }

    if(out.preprod_rigid && !out.protein_membership.empty()) {
        std::set<std::string> lipid_head_roles{"PO4"};
        std::string role_csv = read_string_attribute_or_default(
            ctrl.get(), "preprod_lipid_headgroup_roles", "PO4");
        auto role_tokens = split_csv_tokens(role_csv);
        if(!role_tokens.empty()) {
            lipid_head_roles.clear();
            for(const auto& tok : role_tokens) lipid_head_roles.insert(normalize_role_token(tok));
        }

        out.preprod_fixed_atom_indices.reserve(static_cast<size_t>(n_atom));
        out.preprod_z_fixed_atom_indices.reserve(static_cast<size_t>(n_atom));
        for(int i = 0; i < n_atom; ++i) {
            bool is_protein = (i < (int)out.protein_membership.size() && out.protein_membership[i] >= 0);
            bool is_head = false;
            if(has_roles && i < (int)atom_roles.size()) {
                is_head = (lipid_head_roles.count(atom_roles[i]) > 0);
            }
            if(is_protein) {
                out.preprod_fixed_atom_indices.push_back(i);
            } else if(is_head) {
                out.preprod_z_fixed_atom_indices.push_back(i);
            }
        }
    }

    if(!out.sc_env_po4_env_mask.empty()) {
        out.sc_env_po4_z_hold_atom_indices.reserve(out.sc_env_po4_env_mask.size());
        for(int i = 0; i < n_atom; ++i) {
            if(is_env_po4_atom(out, i)) {
                out.sc_env_po4_z_hold_atom_indices.push_back(i);
            }
        }
    }

    // Fallback BB atom inference if bb_atom_index is absent.
    bool need_infer_bb = false;
    for(int idx : out.bb_atom_index) if(idx < 0) { need_infer_bb = true; break; }
    if(need_infer_bb && h5_exists(root, "/input/residue_ids")) {
        auto residue_shape = get_dset_size(1, root, "/input/residue_ids");
        if(static_cast<int>(residue_shape[0]) == n_atom) {
            std::vector<int> residue_ids(n_atom, -1);
            traverse_dset<1,int>(root, "/input/residue_ids", [&](size_t i, int v) { residue_ids[i] = v; });

            std::vector<std::string> atom_names(n_atom, "");
            bool has_names = false;
            if(h5_exists(root, "/input/atom_names")) {
                has_names = true;
                traverse_string_dset<1>(root, "/input/atom_names", [&](size_t i, const std::string& v) {
                    atom_names[i] = trim_h5_string(v);
                });
            } else if(h5_exists(root, "/input/atom_roles")) {
                has_names = true;
                traverse_string_dset<1>(root, "/input/atom_roles", [&](size_t i, const std::string& v) {
                    atom_names[i] = trim_h5_string(v);
                });
            }

            if(has_names) {
                for(size_t k = 0; k < out.n_bb; ++k) {
                    if(out.bb_atom_index[k] >= 0) continue;
                    int resid = out.bb_residue_index[k];
                    for(int i = 0; i < n_atom; ++i) {
                        if(residue_ids[i] == resid && atom_names[i] == "BB") {
                            out.bb_atom_index[k] = i;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Final hybrid mapping sanity checks in runtime index space.
    if(!out.protein_membership.empty()) {
        for(size_t k = 0; k < out.n_bb; ++k) {
            int bb = out.bb_atom_index[k];
            int ca = (out.atom_mask[k][1] != 0) ? out.atom_indices[k][1] : -1;
            out.bb_ca_atom_index[k] = ca;
            if(bb >= 0) {
                if(bb >= n_atom) {
                    throw string("Hybrid BB proxy index out of bounds");
                }
                if(out.protein_membership[bb] < 0) {
                    throw string("Hybrid BB proxy index must be protein atom");
                }
            }
            if(ca < 0 || ca >= n_atom) {
                throw string("Hybrid BB mapping requires a valid CA carrier index");
            }
            if(out.protein_membership[ca] < 0) {
                throw string("Hybrid BB CA carrier index must be protein atom");
            }
            for(int d = 0; d < 4; ++d) {
                if(out.atom_mask[k][d] == 0) continue;
                int ai = out.atom_indices[k][d];
                if(ai < 0 || ai >= n_atom) {
                    throw string("Hybrid BB target index out of bounds");
                }
                if(out.protein_membership[ai] < 0) {
                    throw string("Hybrid BB target index must be protein atom");
                }
            }
        }

        out.bb_proxy_to_ca_atom.assign(static_cast<size_t>(n_atom), -1);
        out.bb_proxy_to_map_index.assign(static_cast<size_t>(n_atom), -1);
        for(size_t k = 0; k < out.n_bb; ++k) {
            int bb = out.bb_atom_index[k];
            int ca = out.bb_ca_atom_index[k];
            if(bb >= 0 && bb < n_atom && ca >= 0 && ca < n_atom) {
                out.bb_proxy_to_ca_atom[static_cast<size_t>(bb)] = ca;
            }
            if(bb >= 0 && bb < n_atom) {
                out.bb_proxy_to_map_index[static_cast<size_t>(bb)] = static_cast<int>(k);
            }
            if(ca >= 0 && ca < n_atom) {
                if(out.bb_proxy_to_ca_atom[static_cast<size_t>(ca)] < 0) {
                    out.bb_proxy_to_ca_atom[static_cast<size_t>(ca)] = ca;
                }
                if(out.bb_proxy_to_map_index[static_cast<size_t>(ca)] < 0) {
                    out.bb_proxy_to_map_index[static_cast<size_t>(ca)] = static_cast<int>(k);
                }
            }
            for(int d = 0; d < 4; ++d) {
                if(out.atom_mask[k][d] == 0) continue;
                int ai = out.atom_indices[k][d];
                if(ai < 0 || ai >= n_atom) continue;
                if(out.bb_proxy_to_map_index[static_cast<size_t>(ai)] < 0) {
                    out.bb_proxy_to_map_index[static_cast<size_t>(ai)] = static_cast<int>(k);
                }
            }
        }

    }

    return out;
}

void update_stage_for_engine(DerivEngine* engine, const std::string& stage) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(engine);
    if(it == g_hybrid_state.end()) return;
    auto st = it->second;
    if(!st) return;
    bool was_active = st->active;
    if(!st->enabled) {
        st->active = false;
        st->has_prev_bb = false;
        st->prev_bb_pos.clear();
        st->sc_env_po4_z_reference.clear();
        st->sc_env_po4_z_reference_initialized = false;
        st->bb_env_interface_potential = 0.f;
        st->sc_env_transition_step = 0;
        martini_fix_rigid::clear_dynamic_fixed_atoms(*engine);
        martini_fix_rigid::clear_dynamic_z_fixed_atoms(*engine);
        martini_fix_rigid::clear_dynamic_rigid_groups(*engine);
        return;
    }
    st->active = (stage == st->activation_stage);
    if(st->active != was_active) {
        st->has_prev_bb = false;
        st->prev_bb_pos.clear();
        st->sc_env_po4_z_reference.clear();
        st->sc_env_po4_z_reference_initialized = false;
        st->bb_env_interface_potential = 0.f;
        st->sc_env_transition_step = st->sc_env_transition_step_start;
    }
    bool enforce_preprod_rigid = st->preprod_rigid && (stage != "production");
    if(enforce_preprod_rigid) {
        martini_fix_rigid::clear_dynamic_fixed_atoms(*engine);
        martini_fix_rigid::set_dynamic_rigid_groups(*engine, {st->preprod_fixed_atom_indices});
        martini_fix_rigid::set_dynamic_z_fixed_atoms(*engine, st->preprod_z_fixed_atom_indices);
    } else {
        martini_fix_rigid::clear_dynamic_rigid_groups(*engine);
        martini_fix_rigid::set_dynamic_fixed_atoms(*engine, active_protein_proxy_fixed_atoms(*st));
        if(active_sc_env_po4_z_hold_enabled(*st)) {
            martini_fix_rigid::set_dynamic_z_fixed_atoms(*engine, st->sc_env_po4_z_hold_atom_indices);
        } else {
            martini_fix_rigid::clear_dynamic_z_fixed_atoms(*engine);
        }
    }
}

void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine) {
    auto parsed = read_hybrid_settings(config_root, engine.pos->n_elem);
    auto current_stage = martini_stage_params::get_current_stage(&engine);
    parsed.active = parsed.enabled && (current_stage == parsed.activation_stage);

    auto st = std::make_shared<HybridRuntimeState>(std::move(parsed));

    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    g_hybrid_state[&engine] = st;
    g_hybrid_state_by_coord[static_cast<const CoordNode*>(engine.pos)] = st;
    bool enforce_preprod_rigid = st->enabled && st->preprod_rigid && (current_stage != "production");
    if(enforce_preprod_rigid) {
        martini_fix_rigid::clear_dynamic_fixed_atoms(engine);
        martini_fix_rigid::set_dynamic_rigid_groups(engine, {st->preprod_fixed_atom_indices});
        martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st->preprod_z_fixed_atom_indices);
    } else {
        martini_fix_rigid::clear_dynamic_rigid_groups(engine);
        martini_fix_rigid::set_dynamic_fixed_atoms(engine, active_protein_proxy_fixed_atoms(*st));
        if(active_sc_env_po4_z_hold_enabled(*st)) {
            martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st->sc_env_po4_z_hold_atom_indices);
        } else {
            martini_fix_rigid::clear_dynamic_z_fixed_atoms(engine);
        }
    }
}

void refresh_transition_holds_for_engine(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(&engine);
    if(it == g_hybrid_state.end() || !it->second) return;
    auto st = it->second;
    std::string current_stage = martini_stage_params::get_current_stage(&engine);
    bool enforce_preprod_rigid = st->preprod_rigid && (current_stage != "production");
    if(enforce_preprod_rigid) {
        martini_fix_rigid::clear_dynamic_fixed_atoms(engine);
        martini_fix_rigid::set_dynamic_rigid_groups(engine, {st->preprod_fixed_atom_indices});
        martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st->preprod_z_fixed_atom_indices);
        return;
    }
    martini_fix_rigid::clear_dynamic_rigid_groups(engine);
    martini_fix_rigid::set_dynamic_fixed_atoms(engine, active_protein_proxy_fixed_atoms(*st));
    if(st->enabled && st->active &&
       st->sc_env_transition_step < std::numeric_limits<uint64_t>::max()) {
        st->sc_env_transition_step += 1;
    }
    if(active_sc_env_po4_z_hold_enabled(*st)) {
        martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st->sc_env_po4_z_hold_atom_indices);
    } else {
        martini_fix_rigid::clear_dynamic_z_fixed_atoms(engine);
    }
}

bool is_hybrid_enabled(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    return it != g_hybrid_state.end() && it->second && it->second->enabled;
}

bool is_hybrid_active(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    return it != g_hybrid_state.end() && it->second && it->second->active;
}

double get_last_bb_env_interface_potential(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    if(it == g_hybrid_state.end() || !it->second) return 0.0;
    return static_cast<double>(it->second->bb_env_interface_potential);
}

float compute_sc_force_uncap_mix(const HybridRuntimeState& st) {
    if(st.sc_env_relax_steps <= 1) return 1.f;

    uint64_t final_ramp_step = static_cast<uint64_t>(st.sc_env_relax_steps - 1);
    if(st.sc_env_transition_step >= final_ramp_step) return 1.f;

    return float(st.sc_env_transition_step) / float(final_ramp_step);
}

void cap_force_vector(Vec<3>& force, float cap_mag) {
    if(!(cap_mag > 0.f)) return;
    float force2 = mag2(force);
    float cap2 = cap_mag * cap_mag;
    if(force2 > cap2 && force2 > 0.f) {
        force *= cap_mag / sqrtf(force2);
    }
}

bool preproduction_requires_rigid(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    if(it == g_hybrid_state.end() || !it->second) return false;
    return it->second->preprod_rigid;
}

std::shared_ptr<const HybridRuntimeState> get_state_for_coord(const CoordNode& coord) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state_by_coord.find(&coord);
    if(it == g_hybrid_state_by_coord.end()) return nullptr;
    return it->second;
}

bool project_bb_proxy_gradient_for_coord(
        const CoordNode& coord,
        VecArray pos_sens,
        int n_atom,
        int atom_idx,
        const Vec<3>& grad) {
    auto st = get_state_for_coord(coord);
    if(!st || !st->enabled || !st->active) return false;
    if(atom_role_class_at(*st, atom_idx) != ROLE_BB) return false;
    if(atom_is_backbone_carrier_at(*st, atom_idx)) return false;
    if(bb_map_index_for_proxy(*st, atom_idx) < 0) return false;

    Vec<3> routed_grad = compute_sc_backbone_feedback_mix(*st) * grad;
    project_bb_proxy_gradient_if_active(*st, pos_sens, n_atom, atom_idx, routed_grad);
    return true;
}

bool skip_pair_if_intra_protein(const HybridRuntimeState& st, int i, int j) {
    if(!st.enabled || !st.active || !st.exclude_intra_protein_martini) return false;
    return !allow_intra_protein_pair_if_active(st, i, j);
}

void clear_hybrid_for_engine(DerivEngine* engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    if(engine) {
        martini_fix_rigid::clear_dynamic_fixed_atoms(*engine);
        martini_fix_rigid::clear_dynamic_z_fixed_atoms(*engine);
        martini_fix_rigid::clear_dynamic_rigid_groups(*engine);
    }
    if(engine && engine->pos) {
        g_hybrid_state_by_coord.erase(static_cast<const CoordNode*>(engine->pos));
    }
    g_hybrid_state.erase(engine);
}
} // namespace martini_hybrid
