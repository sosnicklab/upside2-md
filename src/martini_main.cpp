#include "main.h"
#include "monte_carlo_sampler.h"
#include "h5_support.h"
#include <tclap/CmdLine.h>
#include "deriv_engine.h"
#include "timing.h"
#include "thermostat.h"
#include <chrono>
#include <algorithm>
#include <set>
#include "random.h"
#include <random>
#include "state_logger.h"
#include <csignal>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>


#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
using namespace h5;

struct MD {
    vector<int> replica_indices;  // Optional: track identity of replicas

    MD(vector<System>& systems) {
        int n_system = systems.size();
        for (int i = 0; i < n_system; ++i) {
            replica_indices.push_back(i);
            auto logger = systems[i].logger;
            if (!logger) continue;
            if (static_cast<int>(logger->level) < static_cast<int>(LOG_BASIC)) continue;

            auto* replica_ptr = &replica_indices[i];
            logger->add_logger<int>("replica_index", {1}, 
                [replica_ptr](int* buffer){ buffer[0] = *replica_ptr; });
        }
    }

    void run_md_step(vector<System>& systems, float timestep) {
        for (auto& sys : systems) {
            sys.engine.compute(PotentialAndDerivMode);
            sys.engine.integrate(timestep);
        }
    }
};