#include <string>

extern "C" {
    int upside_main(int argc, const char* const * argv, int verbose);
}

// Forward declarations
struct DerivEngine;
typedef int64_t hid_t;

// Stage-specific parameter declarations
namespace martini_stage_params {
    void register_stage_params_for_engine(DerivEngine* engine, hid_t root);
    void switch_simulation_stage(DerivEngine* engine, const std::string& new_stage);
    std::string get_current_stage(DerivEngine* engine);
}

namespace martini_hybrid {
    void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine);
    bool is_hybrid_enabled(const DerivEngine& engine);
    bool is_hybrid_active(const DerivEngine& engine);
    bool is_sc_env_energy_dump_enabled(const DerivEngine& engine);
    bool sample_sc_env_energy_for_logging(DerivEngine& engine, float& total, float& lj, float& coul);
    bool preproduction_requires_rigid(const DerivEngine& engine);
    void refresh_transition_holds_for_engine(DerivEngine& engine);
    void clear_hybrid_for_engine(DerivEngine* engine);
}

// Mass storage declarations for MARTINI integrators
namespace martini_masses {
    void load_masses_for_engine(DerivEngine* engine, hid_t config_root);
    float get_mass(DerivEngine* engine, int atom_index);
    bool has_masses(DerivEngine* engine);
}
