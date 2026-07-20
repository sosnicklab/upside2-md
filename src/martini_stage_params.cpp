#include "martini_internal.h"
#include "deriv_engine.h"
#include <mutex>
#include <vector>

using namespace h5;
using namespace std;

namespace martini_stage_params {

// Global registry for stage-specific parameters per engine
static std::mutex g_stage_mutex;
static std::map<DerivEngine*, std::string> g_current_stage;
static std::map<DerivEngine*, std::map<std::string, std::vector<float>>> g_stage_bond_params;
static std::map<DerivEngine*, std::map<std::string, std::vector<float>>> g_stage_angle_params;

// Read stage-specific parameter settings from H5 configuration
struct StageParamData {
    std::string stage;
    std::map<std::string, std::vector<float>> bond_params;
    std::map<std::string, std::vector<float>> angle_params;
    bool has_config;
};

StageParamData read_stage_param_settings(hid_t root) {
    StageParamData data;
    data.has_config = false;
    
    try {
        if(h5_exists(root, "/input/stage_parameters")) {
            auto grp = open_group(root, "/input/stage_parameters");
            data.has_config = true;

            // Set default stage and allow override from H5.
            data.stage = "production";
            data.stage = martini_hybrid::read_string_attribute_or_default(
                grp.get(), "current_stage", data.stage);

            // Read bond parameters for different stages
            if(h5_exists(grp.get(), "minimization_bonds")) {
                auto min_grp = open_group(grp.get(), "minimization_bonds");
                // Read bond force constants
                if(h5_exists(min_grp.get(), "force_constants")) {
                    traverse_dset<1,float>(min_grp.get(), "force_constants", [&](size_t i, float fc) {
                        data.bond_params["minimization"].push_back(fc);
                    });
                }
            }

            if(h5_exists(grp.get(), "production_bonds")) {
                auto prod_grp = open_group(grp.get(), "production_bonds");
                if(h5_exists(prod_grp.get(), "force_constants")) {
                    traverse_dset<1,float>(prod_grp.get(), "force_constants", [&](size_t i, float fc) {
                        data.bond_params["production"].push_back(fc);
                    });
                }
            }

            // Read angle parameters for different stages
            if(h5_exists(grp.get(), "minimization_angles")) {
                auto min_grp = open_group(grp.get(), "minimization_angles");
                if(h5_exists(min_grp.get(), "force_constants")) {
                    traverse_dset<1,float>(min_grp.get(), "force_constants", [&](size_t i, float fc) {
                        data.angle_params["minimization"].push_back(fc);
                    });
                }
            }

            if(h5_exists(grp.get(), "production_angles")) {
                auto prod_grp = open_group(grp.get(), "production_angles");
                if(h5_exists(prod_grp.get(), "force_constants")) {
                    traverse_dset<1,float>(prod_grp.get(), "force_constants", [&](size_t i, float fc) {
                        data.angle_params["production"].push_back(fc);
                    });
                }
            }
        }
    } catch(...) { 
        // Return empty data if no stage parameter settings found
    }
    return data;
}

// Register stage-specific parameters for an engine
void register_stage_params_for_engine(DerivEngine* engine, hid_t root) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    
    auto data = read_stage_param_settings(root);
    if(data.has_config) {
        g_current_stage[engine] = data.stage;
        g_stage_bond_params[engine] = data.bond_params;
        g_stage_angle_params[engine] = data.angle_params;
    }
}

// Switch to a different simulation stage
void switch_simulation_stage(DerivEngine* engine, const std::string& new_stage) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    
    auto it = g_current_stage.find(engine);
    if(it != g_current_stage.end()) {
        it->second = new_stage;
        martini_hybrid::update_stage_for_engine(engine, new_stage);
    }
}

// Get current simulation stage
std::string get_current_stage(DerivEngine* engine) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    
    auto it = g_current_stage.find(engine);
    if(it != g_current_stage.end()) {
        return it->second;
    }
    return "production"; // Default stage
}

// Apply stage-specific bond parameters
void apply_stage_bond_params(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    
    auto stage_it = g_current_stage.find(&engine);
    if(stage_it == g_current_stage.end()) return;
    
    std::string current_stage = stage_it->second;
    auto bond_it = g_stage_bond_params.find(&engine);
    if(bond_it == g_stage_bond_params.end()) return;
    
    auto& stage_bonds = bond_it->second;
    auto param_it = stage_bonds.find(current_stage);
    if(param_it == stage_bonds.end()) return;
    
}

// Apply stage-specific angle parameters  
void apply_stage_angle_params(DerivEngine& engine) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    
    auto stage_it = g_current_stage.find(&engine);
    if(stage_it == g_current_stage.end()) return;
    
    std::string current_stage = stage_it->second;
    auto angle_it = g_stage_angle_params.find(&engine);
    if(angle_it == g_stage_angle_params.end()) return;
    
    auto& stage_angles = angle_it->second;
    auto param_it = stage_angles.find(current_stage);
    if(param_it == stage_angles.end()) return;
    
}

} // namespace martini_stage_params
