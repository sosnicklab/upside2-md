#include "deriv_engine.h"
#include "main.h"
#include "thermostat.h"
#include "timing.h"
#include "random.h"
#include <string>

using namespace std;

void OrnsteinUhlenbeckThermostat::apply(VecArray mom, int n_atom, DerivEngine* engine) {
    Timer timer(string("thermostat"));

    const bool use_mass_aware_noise = engine && martini_masses::has_masses(engine);
    for(int na=0; na<n_atom; ++na) {
        RandomGenerator random(random_seed, THERMOSTAT_RANDOM_STREAM, na, n_invocations);
        auto p = load_vec<3>(mom, na);
        float atom_noise_scale = noise_scale;
        if(use_mass_aware_noise) {
            float mass = martini_masses::get_mass(engine, na);
            atom_noise_scale *= (mass > 0.f) ? sqrtf(mass) : 1.f;
        }
        store_vec(mom, na, mom_scale*p + atom_noise_scale*random.normal3());
    }
    n_invocations++;
}
