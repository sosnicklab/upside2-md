#include "deriv_engine.h"
#include "metadynamics.h"
#include "timing.h"
#include "state_logger.h"
#include <vector>
#include <cmath>

using namespace h5;
using namespace std;

struct Metadynamics : public PotentialNode
{
    CoordNode& cv;
    int   cv_dim;          // dimension of cv.output to bias
    int   cv_elem;         // element of cv.output to bias
    float sigma;           // Gaussian width (CV units)
    float height0;         // hill height (E_up)
    float bias_factor;     // well-tempered gamma (>1); <=1 means standard metadynamics
    int   deposit_period;  // rounds between deposits
    float grid_min, grid_max;  // CV range for the logged bias profile
    int   grid_bins;

    vector<float> hill_center;
    vector<float> hill_height;
    uint64_t last_deposit;
    float last_cv, last_bias;

    Metadynamics(hid_t grp, CoordNode& cv_):
        PotentialNode(),
        cv(cv_),
        cv_dim(         read_attribute<int>  (grp, ".", "cv_dim",  0)),
        cv_elem(        read_attribute<int>  (grp, ".", "cv_elem", 0)),
        sigma(          read_attribute<float>(grp, ".", "sigma")),
        height0(        read_attribute<float>(grp, ".", "height")),
        bias_factor(    read_attribute<float>(grp, ".", "bias_factor", -1.f)),
        deposit_period( read_attribute<int>  (grp, ".", "deposit_period")),
        grid_min(       read_attribute<float>(grp, ".", "grid_min")),
        grid_max(       read_attribute<float>(grp, ".", "grid_max")),
        grid_bins(      read_attribute<int>  (grp, ".", "grid_bins")),
        last_deposit(0), last_cv(0.f), last_bias(0.f)
    {
        if(cv_dim  >= cv.elem_width) throw string("metadynamics cv_dim out of range");
        if(cv_elem >= cv.n_elem)     throw string("metadynamics cv_elem out of range");

        // Restart: reload previously deposited hills if present.
        if(h5_exists(grp, "hill_center") && h5_exists(grp, "hill_height")) {
            traverse_dset<1,float>(grp, "hill_center", [&](size_t, float x){ hill_center.push_back(x); });
            traverse_dset<1,float>(grp, "hill_height", [&](size_t, float x){ hill_height.push_back(x); });
        }

        if(logging(LOG_BASIC)) {
            default_logger->add_logger<float>("metadynamics_cv", {2}, [&](float* buffer) {
                buffer[0] = last_cv;
                buffer[1] = last_bias;
            });
            // The last frame of this grid IS the converged bias -> FES = -(gamma/(gamma-1)) * bias.
            default_logger->add_logger<float>("metadynamics_bias", {grid_bins}, [&](float* buffer) {
                for(int j=0; j<grid_bins; ++j) {
                    float s = grid_min + (grid_max-grid_min)*(j+0.5f)/grid_bins;
                    double V, dV; eval_bias(s, V, dV);
                    buffer[j] = float(V);
                }
            });
        }
    }

    void eval_bias(float s, double& V, double& dV) const {
        V = 0.0; dV = 0.0;
        if(hill_center.empty()) return;
        const double inv_s2 = 1.0/(double(sigma)*double(sigma));
        for(size_t i=0; i<hill_center.size(); ++i) {
            double d = double(s) - hill_center[i];
            double e = hill_height[i] * exp(-0.5*d*d*inv_s2);
            V  += e;
            dV += e * (-d*inv_s2);
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("metadynamics"));
        float s = cv.output(cv_dim, cv_elem);
        double V, dV;
        eval_bias(s, V, dV);
        if(mode == PotentialAndDerivMode) potential = float(V);
        cv.sens(cv_dim, cv_elem) += float(dV);
        last_cv = s; last_bias = float(V);
    }

    void deposit(uint64_t round_num, float temperature) {
        if(deposit_period <= 0) return;
        if(round_num == 0 || (round_num % (uint64_t)deposit_period) != 0) return;
        if(round_num == last_deposit) return;
        last_deposit = round_num;
        float s = cv.output(cv_dim, cv_elem);
        double V, dV;
        eval_bias(s, V, dV);
        float w = height0;
        if(bias_factor > 1.f && temperature > 0.f)          // well-tempered scaling
            w = height0 * expf(-float(V)/(temperature*(bias_factor-1.f)));
        hill_center.push_back(s);
        hill_height.push_back(w);
    }
};
static RegisterNodeType<Metadynamics,1> metadynamics_node("metadynamics");

namespace metadynamics {
    void maybe_deposit(DerivEngine& engine, uint64_t round_num, float temperature) {
        for(auto& node : engine.nodes) {
            auto* m = dynamic_cast<Metadynamics*>(node.computation.get());
            if(m) m->deposit(round_num, temperature);
        }
    }
}
