#include "martini_internal.h"
#include "deriv_engine.h"
#include <mutex>
#include <H5Apublic.h>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <array>
#include <vector>
#include <unordered_map>
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

struct HybridPositionNode : public CoordNode {
    CoordNode& source;
    CoordNode& hbond;   // infer_H_O: supplies the backbone O this proxy is partly built from

    HybridPositionNode(hid_t, CoordNode& source_, CoordNode& hbond_):
        CoordNode(source_.n_elem, 3),
        source(source_),
        hbond(hbond_) {
        check_elem_width_lower_bound(source, 3);
        check_elem_width_lower_bound(hbond, 3);
    }

    virtual void compute_value(ComputeMode) override {
        VecArray source_pos = source.output;
        for(int atom = 0; atom < n_elem; ++atom) {
            store_vec<3>(output, atom, load_vec<3>(source_pos, atom));
        }

        auto st = get_state_for_coord(source);
        if(!st || !st->active) return;
        if(st->bb_reference_runtime_atom_indices.size() != st->n_bb) {
            throw string("Hybrid position node requires a complete BB carrier map");
        }
        VecArray o_pos = hbond.output;
        for(size_t k = 0; k < st->n_bb; ++k) {
            const auto& atom_idx = st->bb_reference_runtime_atom_indices[k];
            int o_idx  = atom_idx[3];
            int bb_idx = st->bb_atom_index[k];
            int o_elem = st->bb_o_hbond_elem[k];
            bool valid = o_idx >= 0 && o_idx < n_elem && bb_idx >= 0 && bb_idx < n_elem &&
                         o_elem < hbond.n_elem;
            for(int d = 0; d < 3; ++d) valid = valid && atom_idx[d] >= 0 && atom_idx[d] < n_elem;
            if(!valid) throw string("Hybrid position node has an invalid BB mapping index");

            // O is Upside's own derived site; take it rather than rebuilding it from a local frame.
            // The C-terminal residue has no acceptor -- infer_H_O builds a carbonyl O from the NEXT
            // residue's N, which does not exist there -- so that one BB is the N/CA/C mass centre with the
            // weights renormalised, and its O slot is left as it came in (no MARTINI term reads it).
            const bool have_o = o_elem >= 0;
            auto o = have_o ? load_vec<3>(o_pos, o_elem) : make_vec3(0.f, 0.f, 0.f);

            // BB is the mass-weighted centre of N/CA/C/O. Because O is a node in its own right, this is a
            // plain linear combination and its derivative needs no placement Jacobian: the O share is
            // handed to infer_H_O, which already carries the frame's rotational dependence exactly.
            auto bb = make_vec3(0.f, 0.f, 0.f);
            float wsum = 0.f;
            for(int d = 0; d < 4; ++d) {
                if(st->atom_mask[k][d] == 0) continue;
                if(d == 3 && !have_o) continue;
                float w = st->weights[k][d];
                if(w == 0.f) continue;
                bb += w * ((d == 3) ? o : load_vec<3>(source_pos, atom_idx[d]));
                wsum += w;
            }
            if(!(wsum > 0.f)) throw string("Hybrid BB weights sum to zero");
            if(fabsf(wsum - 1.f) > 1e-6f) bb *= 1.f/wsum;

            if(have_o) store_vec<3>(output, o_idx, o);
            store_vec<3>(output, bb_idx, bb);
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

        // Both derived sites are linear in nodes the engine already differentiates, so their sensitivity
        // splits by the same constant weights that build them: the N/CA/C shares go straight to pos, and
        // the O share goes to infer_H_O, whose own propagate_deriv carries it back to CA, C and the next
        // residue's N. That is where the frame's rotational dependence lives, so nothing is dropped and no
        // placement Jacobian is needed here.
        VecArray o_sens = hbond.sens;
        for(size_t k = 0; k < st->n_bb; ++k) {
            const auto& atom_idx = st->bb_reference_runtime_atom_indices[k];
            int o_idx  = atom_idx[3];
            int bb_idx = st->bb_atom_index[k];
            int o_elem = st->bb_o_hbond_elem[k];
            if(o_elem >= hbond.n_elem) {
                throw string("Hybrid position node has an out-of-range infer_H_O element for a BB site");
            }
            const bool have_o = o_elem >= 0;

            float wsum = 0.f;
            for(int d = 0; d < 4; ++d) {
                if(st->atom_mask[k][d] == 0) continue;
                if(d == 3 && !have_o) continue;
                wsum += st->weights[k][d];
            }
            if(!(wsum > 0.f)) throw string("Hybrid BB weights sum to zero");
            const float inv_wsum = 1.f/wsum;

            auto bb_s = load_vec<3>(sens, bb_idx);

            // O carries its own sensitivity plus BB's O-weighted share
            auto o_total = have_o ? load_vec<3>(sens, o_idx) : make_vec3(0.f, 0.f, 0.f);
            for(int d = 0; d < 4; ++d) {
                if(st->atom_mask[k][d] == 0) continue;
                if(d == 3 && !have_o) continue;
                float w = st->weights[k][d]*inv_wsum;
                if(w == 0.f) continue;
                if(d == 3) o_total += w*bb_s;
                else       update_vec<3>(source.sens, atom_idx[d], w*bb_s);
            }
            if(have_o) update_vec<3>(o_sens, o_elem, o_total);
            else       update_vec<3>(source.sens, o_idx, load_vec<3>(sens, o_idx));
        }
    }
};

static RegisterNodeType<HybridPositionNode, 2> hybrid_position_node("martini_hybrid_position");

static std::string trim_h5_string(const std::string& in);

static std::string normalize_role_token(const std::string& raw) {
    std::string out = trim_h5_string(raw);
    for(char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
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

// With the interface on, the protein is Upside's alone: no MARTINI term acts within it.
bool allow_intra_protein_pair_if_active(const HybridRuntimeState& st, int i, int j) {
    return !st.active || !same_protein_membership_pair(st, i, j);
}

bool allow_multibody_term_if_active(const HybridRuntimeState& st, const index_t* atoms, int n_atom_dep) {
    if(!st.active) return true;
    for(int a = 0; a < n_atom_dep; ++a)
        for(int b = a + 1; b < n_atom_dep; ++b)
            if(same_protein_membership_pair(st, atoms[a], atoms[b])) return false;
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
    out.preprod_rigid = (out.preprod_mode == "rigid" || out.preprod_mode == "rigid_body");
    out.sc_env_po4_z_hold_steps =
        read_attribute<int>(ctrl.get(), ".", "sc_env_po4_z_hold_steps", out.sc_env_po4_z_hold_steps);
    int transition_start = read_attribute<int>(ctrl.get(), ".", "sc_env_transition_step_start", 0);
    out.sc_env_transition_step_start = static_cast<uint64_t>(std::max(0, transition_start));
    out.sc_env_transition_step = out.sc_env_transition_step_start;
    if(out.sc_env_po4_z_hold_steps < 0) out.sc_env_po4_z_hold_steps = 0;

    if(!h5_exists(root, "/input/hybrid_bb_map")) {
        throw string("Hybrid control requires /input/hybrid_bb_map");
    }
    {
        auto bb = open_group(root, "/input/hybrid_bb_map");
        auto atom_idx_shape = get_dset_size(2, bb.get(), "atom_indices");
        auto atom_mask_shape = get_dset_size(2, bb.get(), "atom_mask");
        auto weights_shape = get_dset_size(2, bb.get(), "weights");
        if(atom_idx_shape[1] != 4 || atom_mask_shape[1] != 4 || weights_shape[1] != 4) {
            throw string("Hybrid BB map must use (n_bb,4) layout for atom_indices/atom_mask/weights");
        }
        if(atom_idx_shape[0] != atom_mask_shape[0] ||
           atom_idx_shape[0] != weights_shape[0]) {
            throw string("Hybrid BB map datasets have inconsistent n_bb sizes");
        }
        out.n_bb = atom_idx_shape[0];
        out.bb_atom_index.assign(out.n_bb, -1);
        out.atom_indices.assign(out.n_bb, std::array<int,4>{{-1,-1,-1,-1}});
        out.atom_mask.assign(out.n_bb, std::array<int,4>{{0,0,0,0}});
        out.weights.assign(out.n_bb, std::array<float,4>{{0.f,0.f,0.f,0.f}});
        out.bb_reference_runtime_atom_indices.assign(out.n_bb, std::array<int,4>{{-1,-1,-1,-1}});

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
        check_size(bb.get(), "bb_atom_index", out.n_bb);
        traverse_dset<1,int>(bb.get(), "bb_atom_index", [&](size_t i, int v) {
            out.bb_atom_index[i] = v;
        });
    }

    if(h5_exists(root, "/input/hybrid_env_topology")) {
        auto env = open_group(root, "/input/hybrid_env_topology");
        auto member_shape = get_dset_size(1, env.get(), "protein_membership");
        if(static_cast<int>(member_shape[0]) != n_atom) {
            throw string("Hybrid env topology protein_membership length must match n_atom");
        }
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

    // The roles name the protein carriers (N, CA, C) and the lipid phosphates (PO4).
    if(!out.protein_membership.empty()) {
        const char* role_dset = h5_exists(root, "/input/atom_roles") ? "/input/atom_roles" : "/input/atom_names";
        check_size(root, role_dset, n_atom);
        std::vector<std::string> atom_roles;
        traverse_string_dset<1>(root, role_dset, [&](size_t i, const std::string& v) {
            atom_roles.push_back(normalize_role_token(v));});

        out.sc_env_po4_env_mask.assign(static_cast<size_t>(n_atom), 0u);
        for(int i = 0; i < n_atom; ++i) {
            const std::string& role = atom_roles[i];
            if(out.protein_membership[i] < 0) {
                if(role == "PO4") out.sc_env_po4_env_mask[i] = 1u;
            } else if(out.preprod_rigid && (role == "N" || role == "CA" || role == "C")) {
                out.preprod_fixed_atom_indices.push_back(i);
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

    // Final hybrid mapping sanity checks in runtime index space.
    if(!out.protein_membership.empty()) {
        for(size_t k = 0; k < out.n_bb; ++k) {
            int bb = out.bb_atom_index[k];
            int ca = (out.atom_mask[k][1] != 0) ? out.atom_indices[k][1] : -1;
            if(bb < 0 || bb >= n_atom) {
                throw string("Hybrid BB proxy index out of bounds");
            }
            if(out.protein_membership[bb] < 0) {
                throw string("Hybrid BB proxy index must be protein atom");
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

        out.bb_proxy_to_map_index.assign(static_cast<size_t>(n_atom), -1);
        for(size_t k = 0; k < out.n_bb; ++k) {
            int bb = out.bb_atom_index[k];
            int ca = (out.atom_mask[k][1] != 0) ? out.atom_indices[k][1] : -1;
            if(bb >= 0 && bb < n_atom) {
                out.bb_proxy_to_map_index[static_cast<size_t>(bb)] = static_cast<int>(k);
            }
            if(ca >= 0 && ca < n_atom && out.bb_proxy_to_map_index[static_cast<size_t>(ca)] < 0) {
                out.bb_proxy_to_map_index[static_cast<size_t>(ca)] = static_cast<int>(k);
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

// A hybrid configuration is built for preparation (activation_stage "minimization") or for production
// (activation_stage "production") and runs only in the stages built for it, every one with the SC-env and
// BB-env interface on:
//   minimization        preparation (stages 6.x): protein N/CA/C move as one rigid body
//   production_handoff  production minimization and burn-in (stage 7.0): protein rigid
//   production          protein flexible; PO4 z held for the first sc_env_po4_z_hold_steps rounds
// The lipids are free in the first two: the preparation restrains the phosphate heights with the harmonic
// restraint_position node the stage file carries, as CHARMM-GUI does, never with a hard freeze.
static void check_hybrid_stage(const HybridRuntimeState& st, const std::string& stage) {
    bool preparation = st.activation_stage == "minimization" && stage == "minimization";
    bool production  = st.activation_stage == "production" &&
                       (stage == "production_handoff" || stage == "production");
    if(!preparation && !production)
        throw string("hybrid stage '") + stage + "' is not a stage of activation_stage '" + st.activation_stage +
              "' (minimization runs in minimization; production in production_handoff and production)";
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
        martini_fix_rigid::clear_dynamic_z_fixed_atoms(engine);
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

// Locate each residue's backbone O inside the infer_H_O node. That node emits donors first, then
// acceptors, and an acceptor is identified by the C atom it hangs off, which is exactly atom_indices[k][2]
// of the BB map. Resolving the element here keeps the position node free of h5 access.
static void resolve_bb_o_hbond_elements(hid_t root, HybridRuntimeState& st) {
    st.bb_o_hbond_elem.assign(st.n_bb, -1);
    if(!st.has_config || !h5_exists(root, "/input/potential/infer_H_O")) return;

    auto hb  = open_group(root, "/input/potential/infer_H_O");
    auto don = open_group(hb.get(), "donors");
    auto acc = open_group(hb.get(), "acceptors");
    int n_donor    = int(get_dset_size(2, don.get(), "id")[0]);
    int n_acceptor = int(get_dset_size(2, acc.get(), "id")[0]);

    // acceptor curr atom (column 1) is the carbonyl C
    std::unordered_map<int,int> c_atom_to_elem;
    traverse_dset<2,int>(acc.get(), "id", [&](size_t i, size_t j, int x) {
        if(j == 1) c_atom_to_elem[x] = n_donor + int(i);
    });
    if(int(c_atom_to_elem.size()) != n_acceptor) {
        throw string("infer_H_O acceptors do not have distinct carbonyl C atoms");
    }
    for(size_t k = 0; k < st.n_bb; ++k) {
        auto it = c_atom_to_elem.find(st.atom_indices[k][2]);
        if(it != c_atom_to_elem.end()) st.bb_o_hbond_elem[k] = it->second;
    }
}

void unregister_hybrid_for_engine(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_hybrid_mutex);
    auto it = g_hybrid_state.find(&engine);
    if(it == g_hybrid_state.end()) return;
    for(auto c = g_hybrid_state_by_coord.begin(); c != g_hybrid_state_by_coord.end(); )
        c = (c->second == it->second) ? g_hybrid_state_by_coord.erase(c) : next(c);
    g_hybrid_state.erase(it);
}

void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine) {
    unregister_hybrid_for_engine(engine);
    auto parsed = read_hybrid_settings(config_root, engine.pos->n_elem);
    resolve_bb_o_hbond_elements(config_root, parsed);
    auto current_stage = martini_stage_params::get_current_stage(&engine);
    if(parsed.has_config) check_hybrid_stage(parsed, current_stage);
    parsed.active = parsed.has_config;

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
    return !allow_intra_protein_pair_if_active(st, i, j);
}

} // namespace martini_hybrid
