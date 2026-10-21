#pragma once

#include "deriv_engine.h"
#include "h5_support.h"
#include "martini_brownian.h"
#include <string>
#include <vector>

struct H5Logger;

namespace martini_fix_rigid {
void register_fix_rigid_for_engine(hid_t config_root, DerivEngine& engine);
void unregister_fix_rigid_for_engine(DerivEngine& engine);
void apply_fix_rigid_minimization(DerivEngine& engine, VecArray pos, VecArray deriv);
void apply_fix_rigid_projection(DerivEngine& engine, VecArray pos);
void apply_fix_rigid_md(DerivEngine& engine, VecArray pos, VecArray mom);
void set_dynamic_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices);
void set_dynamic_z_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices);
void clear_dynamic_z_fixed_atoms(DerivEngine& engine);
void set_dynamic_rigid_groups(DerivEngine& engine, const std::vector<std::vector<int>>& groups);
void clear_dynamic_rigid_groups(DerivEngine& engine);
std::vector<int> get_fixed_atoms(const DerivEngine& engine);
std::vector<int> get_z_fixed_atoms(const DerivEngine& engine);
}

namespace martini_stage_params {
void register_stage_params_for_engine(DerivEngine* engine, hid_t root);
void unregister_stage_params_for_engine(DerivEngine* engine);
std::string get_current_stage(DerivEngine* engine);
}

namespace martini_hybrid {
void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine);
void unregister_hybrid_for_engine(DerivEngine& engine);
bool is_hybrid_enabled(const DerivEngine& engine);
double get_last_bb_env_interface_potential(const DerivEngine& engine);
void refresh_transition_holds_for_engine(DerivEngine& engine);
void copy_positions_for_output(DerivEngine& engine, float* buffer);
}

namespace martini_masses {
void load_masses_for_engine(DerivEngine* engine, hid_t config_root);
void unregister_masses_for_engine(DerivEngine* engine);
float get_mass(DerivEngine* engine, int atom_index);
bool has_masses(DerivEngine* engine);
}

// The MARTINI modules keep per-engine state in maps keyed by the engine's address; an engine's entries must be
// dropped when it is destroyed, or a later engine at the same address inherits them.
void martini_unregister_engine(DerivEngine& engine);

void martini_run_minimization(DerivEngine& engine,
                              int max_iterations,
                              double energy_tolerance,
                              double force_tolerance,
                              double initial_step_size,
                              int verbose);
