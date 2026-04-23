#pragma once

#include "deriv_engine.h"
#include <string>

extern "C" {
    int upside_main(int argc, const char* const * argv, int verbose);
}

namespace martini_fix_rigid {
    void register_fix_rigid_for_engine(hid_t config_root, DerivEngine& engine);
    void apply_fix_rigid_minimization(DerivEngine& engine, VecArray pos, VecArray deriv);
    void apply_fix_rigid_md(DerivEngine& engine, VecArray pos, VecArray deriv, VecArray mom);
    std::vector<int> get_rigid_body_atoms(const DerivEngine& engine);
}

namespace martini_stage_params {
    void register_stage_params_for_engine(DerivEngine* engine, hid_t root);
    void switch_simulation_stage(DerivEngine* engine, const std::string& new_stage);
    std::string get_current_stage(DerivEngine* engine);
}

namespace martini_hybrid {
    void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine);
    bool is_hybrid_enabled(const DerivEngine& engine);
    bool is_hybrid_active(const DerivEngine& engine);
    bool preproduction_requires_rigid(const DerivEngine& engine);
    void refresh_transition_holds_for_engine(DerivEngine& engine);
    void clear_hybrid_for_engine(DerivEngine* engine);
}

namespace martini_masses {
    void load_masses_for_engine(DerivEngine* engine, hid_t config_root);
    float get_mass(DerivEngine* engine, int atom_index);
    bool has_masses(DerivEngine* engine);
}

void martini_run_minimization(DerivEngine& engine,
        int max_iter,
        double energy_tol,
        double force_tol,
        double initial_step,
        int verbose);
