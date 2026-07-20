#include "martini_internal.h"
#include "deriv_engine.h"
#include <mutex>
#include <H5Apublic.h>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <set>
#include <array>
#include <vector>
#include <algorithm>
#include <limits>

using namespace h5;
using namespace std;

namespace martini_hybrid {
static std::mutex g_hybrid_mutex;

static std::map<DerivEngine*, std::shared_ptr<HybridRuntimeState>> g_hybrid_state;
static std::map<const CoordNode*, std::shared_ptr<HybridRuntimeState>> g_hybrid_state_by_coord;

static std::vector<int> active_protein_proxy_fixed_atoms(const HybridRuntimeState& st) {
    std::vector<int> atoms;
    if(!st.active) return atoms;
    atoms.reserve(st.protein_membership.size());
    for(size_t atom_idx = 0; atom_idx < st.protein_membership.size(); ++atom_idx) {
        if(st.protein_membership[atom_idx] < 0) continue;
        if(atom_idx >= st.atom_role_class.size()) continue;
        auto role = st.atom_role_class[atom_idx];
        bool is_bb_source =
            atom_idx < st.atom_bb_source_mask.size() &&
            st.atom_bb_source_mask[atom_idx] != 0u;
        if((role == ROLE_BB || role == ROLE_SC) && !is_bb_source) {
            atoms.push_back(static_cast<int>(atom_idx));
        }
    }
    for(const auto& atom_idx : st.bb_reference_runtime_atom_indices) {
        if(atom_idx[3] >= 0) atoms.push_back(atom_idx[3]);
    }
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    return atoms;
}

int bb_map_index_for_proxy(const HybridRuntimeState& st, int bb_proxy_atom) {
    if(bb_proxy_atom < 0 || bb_proxy_atom >= static_cast<int>(st.bb_proxy_to_map_index.size())) return -1;
    return st.bb_proxy_to_map_index[static_cast<size_t>(bb_proxy_atom)];
}

static inline std::array<float,3> vec_sub(const std::array<float,3>& a, const std::array<float,3>& b) {
    return std::array<float,3>{a[0]-b[0], a[1]-b[1], a[2]-b[2]};
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

struct Dual9 {
    double value = 0.0;
    std::array<double,9> deriv{};
};

using DualVec3 = std::array<Dual9,3>;

static Dual9 dual_constant(double value) {
    Dual9 out;
    out.value = value;
    return out;
}

static Dual9 dual_variable(double value, int index) {
    Dual9 out = dual_constant(value);
    out.deriv[static_cast<size_t>(index)] = 1.0;
    return out;
}

static Dual9 dual_add(const Dual9& a, const Dual9& b) {
    Dual9 out;
    out.value = a.value + b.value;
    for(int i = 0; i < 9; ++i) out.deriv[i] = a.deriv[i] + b.deriv[i];
    return out;
}

static Dual9 dual_sub(const Dual9& a, const Dual9& b) {
    Dual9 out;
    out.value = a.value - b.value;
    for(int i = 0; i < 9; ++i) out.deriv[i] = a.deriv[i] - b.deriv[i];
    return out;
}

static Dual9 dual_mul(const Dual9& a, const Dual9& b) {
    Dual9 out;
    out.value = a.value * b.value;
    for(int i = 0; i < 9; ++i) {
        out.deriv[i] = a.deriv[i]*b.value + a.value*b.deriv[i];
    }
    return out;
}

static Dual9 dual_scale(const Dual9& a, double scale) {
    Dual9 out;
    out.value = a.value * scale;
    for(int i = 0; i < 9; ++i) out.deriv[i] = a.deriv[i] * scale;
    return out;
}

static Dual9 dual_div(const Dual9& a, const Dual9& b) {
    Dual9 out;
    out.value = a.value / b.value;
    double inv_b2 = 1.0 / (b.value*b.value);
    for(int i = 0; i < 9; ++i) {
        out.deriv[i] = (a.deriv[i]*b.value - a.value*b.deriv[i]) * inv_b2;
    }
    return out;
}

static Dual9 dual_sqrt(const Dual9& a) {
    Dual9 out;
    out.value = std::sqrt(a.value);
    double scale = 0.5 / out.value;
    for(int i = 0; i < 9; ++i) out.deriv[i] = a.deriv[i] * scale;
    return out;
}

static DualVec3 dual_vec_sub(const DualVec3& a, const DualVec3& b) {
    return DualVec3{{dual_sub(a[0], b[0]), dual_sub(a[1], b[1]), dual_sub(a[2], b[2])}};
}

static Dual9 dual_vec_dot(const DualVec3& a, const DualVec3& b) {
    Dual9 out = dual_constant(0.0);
    for(int d = 0; d < 3; ++d) out = dual_add(out, dual_mul(a[d], b[d]));
    return out;
}

static DualVec3 dual_vec_cross(const DualVec3& a, const DualVec3& b) {
    return DualVec3{{
        dual_sub(dual_mul(a[1], b[2]), dual_mul(a[2], b[1])),
        dual_sub(dual_mul(a[2], b[0]), dual_mul(a[0], b[2])),
        dual_sub(dual_mul(a[0], b[1]), dual_mul(a[1], b[0]))
    }};
}

static bool dual_vec_normalize(const DualVec3& value, DualVec3& normalized) {
    Dual9 norm2 = dual_vec_dot(value, value);
    if(norm2.value <= 1e-16) return false;
    Dual9 norm = dual_sqrt(norm2);
    for(int d = 0; d < 3; ++d) normalized[d] = dual_div(value[d], norm);
    return true;
}

static bool map_backbone_sites(
        const HybridRuntimeState& st,
        size_t k,
        const std::array<std::array<float,3>,3>& carrier,
        std::array<float,3>& mapped_o,
        std::array<float,3>& mapped_bb,
        std::array<float,27>& jacobian) {
    const auto& ref = st.bb_reference_atom_coords[k];
    float F_ref[3][3];
    if(!build_frame_from_three(ref[0], ref[1], ref[2], F_ref)) return false;

    auto ref_local_o = vec_sub(ref[3], ref[1]);
    std::array<double,3> local_o{0.0, 0.0, 0.0};
    for(int basis = 0; basis < 3; ++basis) {
        for(int d = 0; d < 3; ++d) {
            local_o[basis] += F_ref[d][basis] * ref_local_o[d];
        }
    }

    std::array<DualVec3,3> carrier_dual;
    for(int atom = 0; atom < 3; ++atom) {
        for(int d = 0; d < 3; ++d) {
            carrier_dual[atom][d] = dual_variable(carrier[atom][d], 3*atom + d);
        }
    }
    DualVec3 e1;
    if(!dual_vec_normalize(dual_vec_sub(carrier_dual[1], carrier_dual[0]), e1)) return false;
    DualVec3 v2 = dual_vec_sub(carrier_dual[2], carrier_dual[0]);
    DualVec3 e3;
    if(!dual_vec_normalize(dual_vec_cross(e1, v2), e3)) return false;
    DualVec3 e2 = dual_vec_cross(e3, e1);

    DualVec3 mapped_o_dual;
    for(int d = 0; d < 3; ++d) {
        mapped_o_dual[d] = carrier_dual[1][d];
        mapped_o_dual[d] = dual_add(mapped_o_dual[d], dual_scale(e1[d], local_o[0]));
        mapped_o_dual[d] = dual_add(mapped_o_dual[d], dual_scale(e2[d], local_o[1]));
        mapped_o_dual[d] = dual_add(mapped_o_dual[d], dual_scale(e3[d], local_o[2]));
    }

    std::array<DualVec3,4> sites{{
        carrier_dual[0], carrier_dual[1], carrier_dual[2], mapped_o_dual
    }};
    DualVec3 mapped_bb_dual{{
        dual_constant(0.0), dual_constant(0.0), dual_constant(0.0)
    }};
    float wsum = 0.f;
    for(int d = 0; d < 4; ++d) {
        if(st.atom_mask[k][d] == 0) continue;
        float w = st.weights[k][d];
        if(w == 0.f) continue;
        for(int coord = 0; coord < 3; ++coord) {
            mapped_bb_dual[coord] =
                dual_add(mapped_bb_dual[coord], dual_scale(sites[d][coord], w));
        }
        wsum += w;
    }
    if(wsum <= 0.f) return false;
    if(fabsf(wsum - 1.f) > 1e-6f) {
        for(int d = 0; d < 3; ++d) mapped_bb_dual[d] = dual_scale(mapped_bb_dual[d], 1.f/wsum);
    }
    for(int output = 0; output < 3; ++output) {
        mapped_o[output] = static_cast<float>(mapped_o_dual[output].value);
        mapped_bb[output] = static_cast<float>(mapped_bb_dual[output].value);
        for(int input = 0; input < 9; ++input) {
            jacobian[9*output + input] =
                static_cast<float>(mapped_bb_dual[output].deriv[input]);
        }
    }
    return true;
}

struct HybridPositionNode : public CoordNode {
    CoordNode& source;
    std::vector<std::array<float,27>> jacobian;

    HybridPositionNode(hid_t, CoordNode& source_):
        CoordNode(source_.n_elem, 3),
        source(source_) {
        check_elem_width_lower_bound(source, 3);
    }

    virtual void compute_value(ComputeMode) override {
        VecArray source_pos = source.output;
        for(int atom = 0; atom < n_elem; ++atom) {
            store_vec<3>(output, atom, load_vec<3>(source_pos, atom));
        }

        auto st = get_state_for_coord(source);
        if(!st || !st->active) return;
        if(st->bb_reference_atom_coords.size() != st->n_bb ||
           st->bb_reference_runtime_atom_indices.size() != st->n_bb) {
            throw string("Hybrid position node requires complete BB reference data");
        }
        jacobian.resize(st->n_bb);
        for(size_t k = 0; k < st->n_bb; ++k) {
            const auto& atom_idx = st->bb_reference_runtime_atom_indices[k];
            int o_idx = atom_idx[3];
            int bb_idx = st->bb_atom_index[k];
            bool valid = o_idx >= 0 && o_idx < n_elem && bb_idx >= 0 && bb_idx < n_elem;
            for(int d = 0; d < 3; ++d) valid = valid && atom_idx[d] >= 0 && atom_idx[d] < n_elem;
            if(!valid) throw string("Hybrid position node has an invalid BB mapping index");

            std::array<std::array<float,3>,3> carrier;
            for(int atom = 0; atom < 3; ++atom) {
                auto value = load_vec<3>(source_pos, atom_idx[atom]);
                carrier[atom] = std::array<float,3>{{value[0], value[1], value[2]}};
            }
            std::array<float,3> mapped_o, mapped_bb;
            if(!map_backbone_sites(*st, k, carrier, mapped_o, mapped_bb, jacobian[k])) {
                throw string("Hybrid position node encountered a singular BB frame");
            }
            store_vec<3>(output, o_idx, make_vec3(mapped_o[0], mapped_o[1], mapped_o[2]));
            store_vec<3>(output, bb_idx, make_vec3(mapped_bb[0], mapped_bb[1], mapped_bb[2]));
        }
    }

    virtual void propagate_deriv() override {
        auto st = get_state_for_coord(source);
        if(!st || !st->active) {
            for(int atom = 0; atom < n_elem; ++atom) {
                update_vec<3>(source.sens, atom, load_vec<3>(sens, atom));
            }
            return;
        }
        if(jacobian.size() != st->n_bb) {
            throw string("Hybrid position node Jacobian is not initialized");
        }

        std::vector<unsigned char> derived(static_cast<size_t>(n_elem), 0u);
        for(size_t k = 0; k < st->n_bb; ++k) {
            int o_idx = st->bb_reference_runtime_atom_indices[k][3];
            int bb_idx = st->bb_atom_index[k];
            if(o_idx >= 0 && o_idx < n_elem) derived[static_cast<size_t>(o_idx)] = 1u;
            if(bb_idx >= 0 && bb_idx < n_elem) derived[static_cast<size_t>(bb_idx)] = 1u;
        }
        for(int atom = 0; atom < n_elem; ++atom) {
            if(!derived[static_cast<size_t>(atom)]) {
                update_vec<3>(source.sens, atom, load_vec<3>(sens, atom));
            }
        }
        for(size_t k = 0; k < st->n_bb; ++k) {
            int bb_idx = st->bb_atom_index[k];
            Vec<3> grad = load_vec<3>(sens, bb_idx);
            const auto& jac = jacobian[k];
            for(int carrier = 0; carrier < 3; ++carrier) {
                int atom_idx = st->atom_indices[k][carrier];
                Vec<3> projected = make_zero<3>();
                for(int input_dim = 0; input_dim < 3; ++input_dim) {
                    int input = 3*carrier + input_dim;
                    for(int output_dim = 0; output_dim < 3; ++output_dim) {
                        projected[input_dim] += jac[9*output_dim + input] * grad[output_dim];
                    }
                }
                update_vec<3>(source.sens, atom_idx, projected);
            }
        }
    }
};

static RegisterNodeType<HybridPositionNode, 1> hybrid_position_node("martini_hybrid_position");

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

static inline bool is_bb_source_role_name(const std::string& raw_name) {
    std::string name = normalize_role_token(raw_name);
    return name == "N" || name == "CA" || name == "C" || name == "O";
}

static inline bool is_env_po4_atom(const HybridRuntimeState& st, int atom) {
    return atom >= 0 &&
           atom < static_cast<int>(st.sc_env_po4_env_mask.size()) &&
           st.sc_env_po4_env_mask[atom] != 0;
}

static inline bool active_sc_env_backbone_hold_enabled(const HybridRuntimeState& st) {
    return st.active &&
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

static inline bool active_sc_env_po4_z_hold_enabled(const HybridRuntimeState& st) {
    return st.active &&
           st.sc_env_po4_z_hold_steps > 0 &&
           st.sc_env_transition_step < static_cast<uint64_t>(st.sc_env_po4_z_hold_steps);
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

bool atom_is_bb_source_at(const HybridRuntimeState& st, int i) {
    return i >= 0 &&
           i < static_cast<int>(st.atom_bb_source_mask.size()) &&
           st.atom_bb_source_mask[static_cast<size_t>(i)] != 0;
}

static inline bool allow_protein_pair_by_rule(const HybridRuntimeState& st, int i, int j) {
    (void)st;
    (void)i;
    (void)j;
    // In the direct-Upside active stage, protein internal MARTINI proxy terms
    // are bookkeeping only. Keeping any proxy-proxy protein interaction live
    // adds proxy bonded/nonbonded energy to the production bucket without
    // feeding that force back through the Upside carrier path.
    return false;
}

bool allow_intra_protein_pair_if_active(const HybridRuntimeState& st, int i, int j) {
    if(!st.active) return true;
    if(!same_protein_membership_pair(st, i, j)) return true;
    return allow_protein_pair_by_rule(st, i, j);
}

bool allow_multibody_term_if_active(const HybridRuntimeState& st, const index_t* atoms, int n_atom_dep) {
    if(!st.active) return true;
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
    out.activation_stage = read_string_attribute_or_default(ctrl.get(), "activation_stage", "production");
    out.preprod_mode = read_string_attribute_or_default(ctrl.get(), "preprod_protein_mode", "rigid");
    out.exclude_intra_protein_martini =
        (read_attribute<int>(ctrl.get(), ".", "exclude_intra_protein_martini", 1) != 0);
    out.preprod_rigid = (out.preprod_mode == "rigid" || out.preprod_mode == "rigid_body");
    out.sc_env_backbone_hold_steps =
        read_attribute<int>(ctrl.get(), ".", "sc_env_backbone_hold_steps", out.sc_env_backbone_hold_steps);
    out.sc_env_po4_z_hold_steps =
        read_attribute<int>(ctrl.get(), ".", "sc_env_po4_z_hold_steps", out.sc_env_po4_z_hold_steps);
    int transition_start = read_attribute<int>(ctrl.get(), ".", "sc_env_transition_step_start", 0);
    out.sc_env_transition_step_start = static_cast<uint64_t>(std::max(0, transition_start));
    out.sc_env_transition_step = out.sc_env_transition_step_start;
    if(out.sc_env_backbone_hold_steps < 0) out.sc_env_backbone_hold_steps = 0;
    if(out.sc_env_po4_z_hold_steps < 0) out.sc_env_po4_z_hold_steps = 0;

    if(!h5_exists(root, "/input/hybrid_bb_map")) {
        throw string("Hybrid control requires /input/hybrid_bb_map");
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
        throw string("Hybrid control requires /input/hybrid_env_topology");
    }
    out.atom_role_class.assign(static_cast<size_t>(n_atom), ROLE_OTHER);
    out.atom_bb_source_mask.assign(static_cast<size_t>(n_atom), 0u);
    if(h5_exists(root, "/input/atom_roles")) {
        traverse_string_dset<1>(root, "/input/atom_roles", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) {
                out.atom_role_class[i] = classify_atom_role_name(v);
                out.atom_bb_source_mask[i] = is_bb_source_role_name(v) ? 1u : 0u;
            }
        });
    } else if(h5_exists(root, "/input/atom_names")) {
        traverse_string_dset<1>(root, "/input/atom_names", [&](size_t i, const std::string& v) {
            if(static_cast<int>(i) < n_atom) {
                out.atom_role_class[i] = classify_atom_role_name(v);
                out.atom_bb_source_mask[i] = is_bb_source_role_name(v) ? 1u : 0u;
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

    if(!out.protein_membership.empty()) {
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
            bool is_carrier = false;
            if(is_protein && has_roles && i < static_cast<int>(atom_roles.size())) {
                const std::string& role = atom_roles[i];
                is_carrier = role == "N" || role == "CA" || role == "C";
            }
            if(is_carrier) {
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

static inline bool enforce_preprod_rigid_stage(
        const HybridRuntimeState& st,
        const std::string& stage) {
    return st.preprod_rigid && (stage != "production");
}

static void apply_stage_fixing(
        DerivEngine& engine,
        HybridRuntimeState& st,
        const std::string& stage) {
    if(enforce_preprod_rigid_stage(st, stage)) {
        martini_fix_rigid::set_dynamic_fixed_atoms(engine, active_protein_proxy_fixed_atoms(st));
        martini_fix_rigid::set_dynamic_rigid_groups(engine, {st.preprod_fixed_atom_indices});
        martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st.preprod_z_fixed_atom_indices);
        return;
    }

    martini_fix_rigid::clear_dynamic_rigid_groups(engine);
    martini_fix_rigid::set_dynamic_fixed_atoms(engine, active_protein_proxy_fixed_atoms(st));
    if(active_sc_env_po4_z_hold_enabled(st)) {
        martini_fix_rigid::set_dynamic_z_fixed_atoms(engine, st.sc_env_po4_z_hold_atom_indices);
    } else {
        martini_fix_rigid::clear_dynamic_z_fixed_atoms(engine);
    }
}

void update_stage_for_engine(DerivEngine* engine, const std::string& stage) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(engine);
    if(it == g_hybrid_state.end()) return;
    auto st = it->second;
    if(!st) return;
    bool was_active = st->active;
    if(!st->has_config) {
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
    apply_stage_fixing(*engine, *st, stage);
}

void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine) {
    auto parsed = read_hybrid_settings(config_root, engine.pos->n_elem);
    auto current_stage = martini_stage_params::get_current_stage(&engine);
    parsed.active = parsed.has_config && (current_stage == parsed.activation_stage);

    auto st = std::make_shared<HybridRuntimeState>(std::move(parsed));

    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    g_hybrid_state[&engine] = st;
    g_hybrid_state_by_coord[static_cast<const CoordNode*>(engine.pos)] = st;
    for(auto& node : engine.nodes) {
        if(node.name != "martini_hybrid_position") continue;
        auto* coord = dynamic_cast<CoordNode*>(node.computation.get());
        if(coord) g_hybrid_state_by_coord[coord] = st;
    }
    apply_stage_fixing(engine, *st, current_stage);
}

void refresh_transition_holds_for_engine(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(&engine);
    if(it == g_hybrid_state.end() || !it->second) return;
    auto st = it->second;
    std::string current_stage = martini_stage_params::get_current_stage(&engine);
    if(!enforce_preprod_rigid_stage(*st, current_stage) &&
       st->active &&
       st->sc_env_transition_step < std::numeric_limits<uint64_t>::max()) {
        st->sc_env_transition_step += 1;
    }
    apply_stage_fixing(engine, *st, current_stage);
}

bool is_hybrid_enabled(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    return it != g_hybrid_state.end() && it->second && it->second->has_config;
}

double get_last_bb_env_interface_potential(const DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(const_cast<DerivEngine*>(&engine));
    if(it == g_hybrid_state.end() || !it->second) return 0.0;
    return static_cast<double>(it->second->bb_env_interface_potential);
}

void copy_positions_for_output(DerivEngine& engine, float* buffer) {
    for(auto& node : engine.nodes) {
        if(node.name != "martini_hybrid_position") continue;
        auto* coord = dynamic_cast<CoordNode*>(node.computation.get());
        if(!coord) break;
        copy_vec_array_to_buffer(coord->output, coord->n_elem, 3, buffer);
        return;
    }
    copy_vec_array_to_buffer(engine.pos->output, engine.pos->n_elem, 3, buffer);
}

std::shared_ptr<const HybridRuntimeState> get_state_for_coord(const CoordNode& coord) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state_by_coord.find(&coord);
    if(it == g_hybrid_state_by_coord.end()) return nullptr;
    return it->second;
}

bool skip_pair_if_intra_protein(const HybridRuntimeState& st, int i, int j) {
    if(!st.active || !st.exclude_intra_protein_martini) return false;
    return !allow_intra_protein_pair_if_active(st, i, j);
}

} // namespace martini_hybrid
