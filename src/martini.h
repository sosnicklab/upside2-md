#pragma once

#include "deriv_engine.h"
#include "h5_support.h"
#include <string>
#include <vector>

struct H5Logger;

namespace martini_fix_rigid {
void register_fix_rigid_for_engine(hid_t config_root, DerivEngine& engine);
void apply_fix_rigid_minimization(DerivEngine& engine, VecArray pos, VecArray deriv);
void apply_fix_rigid_projection(DerivEngine& engine, VecArray pos);
void apply_fix_rigid_md(DerivEngine& engine, VecArray pos, VecArray deriv, VecArray mom);
void set_dynamic_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices);
void clear_dynamic_fixed_atoms(DerivEngine& engine);
void set_dynamic_z_fixed_atoms(DerivEngine& engine, const std::vector<int>& atom_indices);
void clear_dynamic_z_fixed_atoms(DerivEngine& engine);
void set_dynamic_rigid_groups(DerivEngine& engine, const std::vector<std::vector<int>>& groups);
void clear_dynamic_rigid_groups(DerivEngine& engine);
std::vector<int> get_fixed_atoms(const DerivEngine& engine);
std::vector<int> get_z_fixed_atoms(const DerivEngine& engine);
}

namespace martini_stage_params {
void register_stage_params_for_engine(DerivEngine* engine, hid_t root);
void switch_simulation_stage(DerivEngine* engine, const std::string& new_stage);
std::string get_current_stage(DerivEngine* engine);
void apply_stage_bond_params(DerivEngine& engine);
void apply_stage_angle_params(DerivEngine& engine);
void clear_stage_params_for_engine(DerivEngine* engine);
}

namespace martini_hybrid {
void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine);
bool is_hybrid_enabled(const DerivEngine& engine);
bool is_hybrid_active(const DerivEngine& engine);
double get_last_bb_env_interface_potential(const DerivEngine& engine);
bool preproduction_requires_rigid(const DerivEngine& engine);
void refresh_transition_holds_for_engine(DerivEngine& engine);
void clear_hybrid_for_engine(DerivEngine* engine);
bool project_bb_proxy_gradient_for_coord(
        const CoordNode& coord,
        VecArray pos_sens,
        int n_atom,
        int atom_idx,
        const Vec<3>& grad);
}

namespace martini_masses {
void load_masses_for_engine(DerivEngine* engine, hid_t config_root);
float get_mass(DerivEngine* engine, int atom_index);
void clear_masses_for_engine(DerivEngine* engine);
bool has_masses(DerivEngine* engine);
}

namespace martini_cg_lipid {
void register_dynamic_compaction_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);
void clear_dynamic_compaction_for_engine(DerivEngine* engine);
bool has_dynamic_compaction(DerivEngine* engine);
void set_dynamic_compaction_temperature(DerivEngine* engine, float temperature);
void set_dynamic_compaction_thermostat_delta_t(DerivEngine* engine, float delta_t);
void apply_dynamic_compaction_thermostat(DerivEngine* engine);
void integrate_dynamic_compaction(DerivEngine* engine, float dt);
void add_dynamic_compaction_loggers(DerivEngine* engine, H5Logger& logger, bool record_momentum);
void register_dynamic_orientation_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);
void clear_dynamic_orientation_for_engine(DerivEngine* engine);
bool has_dynamic_orientation(DerivEngine* engine);
void set_dynamic_orientation_temperature(DerivEngine* engine, float temperature);
void set_dynamic_orientation_thermostat_delta_t(DerivEngine* engine, float delta_t);
void apply_dynamic_orientation_thermostat(DerivEngine* engine);
void integrate_dynamic_orientation(DerivEngine* engine, float dt);
void add_dynamic_orientation_loggers(DerivEngine* engine, H5Logger& logger, bool record_momentum);
void register_cgl_gle_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);
void clear_cgl_gle_for_engine(DerivEngine* engine);
bool has_cgl_gle(DerivEngine* engine);
void set_cgl_gle_temperature(DerivEngine* engine, float temperature);
void set_cgl_gle_delta_t(DerivEngine* engine, float delta_t);
void apply_cgl_gle_thermostat(DerivEngine* engine, VecArray mom);
void add_cgl_gle_loggers(DerivEngine* engine, H5Logger& logger, bool record_aux);
}

void martini_run_minimization(DerivEngine& engine,
                              int max_iterations,
                              double energy_tolerance,
                              double force_tolerance,
                              double initial_step_size,
                              int verbose);
