#include "deriv_engine.h"
#include "martini.h"
#include "thermostat.h"
#include "timing.h"
#include "random.h"
#include <string>

using namespace std;

void OrnsteinUhlenbeckThermostat::apply(
        VecArray mom,
        int n_atom,
        DerivEngine* engine,
        bool include_brownian,
        bool full_refresh) {
    Timer timer(string("thermostat"));

    const bool use_mass_aware_noise = engine && martini_masses::has_masses(engine);
    const bool use_atom_timescale = atom_timescale.size() == static_cast<size_t>(n_atom);
    const auto& stochastic_mask = engine ? martini_brownian::stochastic_mask(engine)
                                         : std::vector<unsigned char>();
    const bool use_stochastic_mask = stochastic_mask.size() == static_cast<size_t>(n_atom);
    std::vector<unsigned char> fixed_mask(static_cast<size_t>(n_atom), 0u);
    if(engine) {
        for(int atom : martini_fix_rigid::get_fixed_atoms(*engine)) {
            if(atom >= 0 && atom < n_atom) fixed_mask[static_cast<size_t>(atom)] = 1u;
        }
    }
    for(int na=0; na<n_atom; ++na) {
        if(fixed_mask[static_cast<size_t>(na)]) continue;
        if(!include_brownian && use_stochastic_mask && stochastic_mask[static_cast<size_t>(na)]) continue;
        RandomGenerator random(random_seed, THERMOSTAT_RANDOM_STREAM, na, n_invocations);
        auto p = load_vec<3>(mom, na);
        float atom_mom_scale = full_refresh ? 0.f : mom_scale;
        float atom_noise_scale = full_refresh ? sqrtf(temp) : noise_scale;
        if(use_atom_timescale && !full_refresh) {
            float tau = atom_timescale[static_cast<size_t>(na)];
            atom_mom_scale = expf(-delta_t/tau);
            atom_noise_scale = sqrtf(temp * (1.f-atom_mom_scale*atom_mom_scale));
        }
        if(use_mass_aware_noise) {
            float mass = martini_masses::get_mass(engine, na);
            atom_noise_scale *= (mass > 0.f) ? sqrtf(mass) : 1.f;
        }
        store_vec(mom, na, atom_mom_scale*p + atom_noise_scale*random.normal3());
    }
    n_invocations++;
}
