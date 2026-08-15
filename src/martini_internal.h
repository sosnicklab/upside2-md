#pragma once

#include "martini.h"
#include <array>
#include <cstdint>
#include <memory>

namespace martini_hybrid {

enum AtomRoleClass : unsigned char {
    ROLE_OTHER = 0,
    ROLE_BB = 1,
    ROLE_SC = 2
};

struct HybridRuntimeState {
    bool has_config = false;
    bool active = false;
    bool exclude_intra_protein_martini = true;
    bool preprod_rigid = true;
    std::string activation_stage = "production";
    std::string preprod_mode = "rigid_body";
    size_t n_bb = 0;
    size_t n_env = 0;
    std::vector<int> bb_residue_index;
    std::vector<int> bb_atom_index;
    std::vector<int> bb_proxy_to_map_index;
    std::vector<std::array<int,4>> atom_indices;
    std::vector<std::array<int,4>> atom_mask;
    std::vector<std::array<float,4>> weights;
    std::vector<std::array<int,4>> bb_reference_runtime_atom_indices;
    // Element of the infer_H_O node holding each residue's backbone O. The BB proxy is the mass-weighted
    // centre of N/CA/C/O, and O is Upside's own derived site, so taking it from that node makes the proxy
    // an exact linear combination of nodes the engine already differentiates.
    std::vector<int> bb_o_hbond_elem;
    std::vector<int> protein_membership;
    std::vector<unsigned char> atom_role_class;
    std::vector<unsigned char> atom_bb_source_mask;
    int sc_env_po4_z_hold_steps = 150;
    uint64_t sc_env_transition_step = 0;
    uint64_t sc_env_transition_step_start = 0;
    std::vector<unsigned char> sc_env_po4_env_mask;
    std::vector<float> sc_env_po4_z_reference;
    bool sc_env_po4_z_reference_initialized = false;
    std::vector<int> sc_env_po4_z_hold_atom_indices;
    float bb_env_interface_potential = 0.f;
    std::vector<int> preprod_fixed_atom_indices;
    std::vector<int> preprod_z_fixed_atom_indices;
};

int bb_map_index_for_proxy(const HybridRuntimeState& st, int bb_proxy_atom);
void initialize_sc_env_po4_z_reference(HybridRuntimeState& st, VecArray pos, int n_atom);
unsigned char atom_role_class_at(const HybridRuntimeState& st, int i);
bool atom_is_bb_source_at(const HybridRuntimeState& st, int i);
bool allow_intra_protein_pair_if_active(const HybridRuntimeState& st, int i, int j);
bool allow_multibody_term_if_active(const HybridRuntimeState& st, const index_t* atoms, int n_atom_dep);
std::shared_ptr<const HybridRuntimeState> get_state_for_coord(const CoordNode& coord);
bool skip_pair_if_intra_protein(const HybridRuntimeState& st, int i, int j);
std::string read_string_attribute_or_default(
        hid_t group,
        const char* attr_name,
        const std::string& fallback);
void update_stage_for_engine(DerivEngine* engine, const std::string& stage);

} // namespace martini_hybrid
