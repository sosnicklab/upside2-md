#include "martini.h"
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include <mutex>
#include "spline.h"
#include <iostream>
#include <H5Apublic.h> // for H5Aexists
#include <cmath> // For pow, cosf, sinf, acosf
#include <cctype>
#include <cstdint>
#include <set> // For std::set
#include <array>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <limits>
#include "box.h" // For PBC minimum_image function

using namespace h5;
using namespace std;

namespace martini_masses {
    static std::mutex g_mass_mutex;
    static std::map<DerivEngine*, std::vector<float>> g_masses;

    void load_masses_for_engine(DerivEngine* engine, hid_t config_root) {
        std::lock_guard<std::mutex> lk(g_mass_mutex);
        auto& masses = g_masses[engine];
        masses.clear();
        
        try {
            if(h5_exists(config_root, "/input/mass")) {
                traverse_dset<1,float>(config_root, "/input/mass", [&](size_t i, float m){ masses.push_back(m); });
            }
        } catch(...) {
            // Silently handle errors
        }
    }
    
    // Get mass for a specific atom
    float get_mass(DerivEngine* engine, int atom_index) {
        std::lock_guard<std::mutex> lk(g_mass_mutex);
        auto it = g_masses.find(engine);
        if(it == g_masses.end()) return 1.0f;
        if(atom_index < 0 || atom_index >= (int)it->second.size()) return 1.0f;
        return it->second[atom_index];
    }
    
    // Clean up masses for an engine
    void clear_masses_for_engine(DerivEngine* engine) {
        std::lock_guard<std::mutex> lk(g_mass_mutex);
        g_masses.erase(engine);
    }
    
    bool has_masses(DerivEngine* engine) {
        std::lock_guard<std::mutex> lk(g_mass_mutex);
        auto it = g_masses.find(engine);
        return (it != g_masses.end() && !it->second.empty());
    }
}
