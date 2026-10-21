#include "martini_internal.h"
#include "deriv_engine.h"
#include <mutex>

using namespace h5;
using namespace std;

namespace martini_stage_params {

// Global registry of the current simulation stage per engine
static std::mutex g_stage_mutex;
static std::map<DerivEngine*, std::string> g_current_stage;

void unregister_stage_params_for_engine(DerivEngine* engine) {
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    g_current_stage.erase(engine);
}

// Read /input/stage_parameters@current_stage ("production" if the attribute is absent)
void register_stage_params_for_engine(DerivEngine* engine, hid_t root) {
    unregister_stage_params_for_engine(engine);
    if(!h5_exists(root, "/input/stage_parameters")) return;
    auto grp = open_group(root, "/input/stage_parameters");
    std::string stage = martini_hybrid::read_string_attribute_or_default(grp.get(), "current_stage", "production");
    std::lock_guard<std::mutex> lock(g_stage_mutex);
    g_current_stage[engine] = stage;
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

} // namespace martini_stage_params
