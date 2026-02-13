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
    void apply_stage_bond_params(DerivEngine& engine);
    void apply_stage_angle_params(DerivEngine& engine);
    void clear_stage_params_for_engine(DerivEngine* engine);
}

namespace martini_hybrid {
    void register_hybrid_for_engine(hid_t config_root, DerivEngine& engine);
    bool is_hybrid_enabled(const DerivEngine& engine);
    bool is_hybrid_active(const DerivEngine& engine);
    bool preproduction_requires_rigid(const DerivEngine& engine);
    void clear_hybrid_for_engine(DerivEngine* engine);
}

// Mass storage declarations for MARTINI integrators
namespace martini_masses {
    void load_masses_for_engine(DerivEngine* engine, hid_t config_root);
    float get_mass(DerivEngine* engine, int atom_index);
    void clear_masses_for_engine(DerivEngine* engine);
    bool has_masses(DerivEngine* engine);
}
