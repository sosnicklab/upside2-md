#pragma once

#include "deriv_engine.h"
#include <cstdint>
#include <vector>

// Single-step inertial g-JF Langevin integrator for a MARTINI system. Every dynamic particle is
// listed in /input/brownian; zero-friction protein carriers reduce exactly to velocity Verlet and
// positive-friction environment particles receive the g-JF bath.
namespace martini_brownian {

// Read /input/brownian (atom_index and absolute friction datasets plus the numerical timestep).
// No-op if the group is absent.
void register_brownian_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);

// True if this engine has g-JF atoms.
bool has_brownian(DerivEngine* engine);

// Per-atom mask for particles that receive stochastic noise in this integrator. A zero-friction
// g-JF particle is thermostatted by the global OU bath instead.
const std::vector<unsigned char>& stochastic_mask(DerivEngine* engine);

// Set the target temperature (kT in E_up) used for the Brownian noise.
void set_brownian_temperature(DerivEngine* engine, float temperature);

// Advance every dynamic atom by one g-JF step of size dt using pos->sens (= dU/dx).
void apply_langevin_step(DerivEngine* engine, VecArray mom, float dt);

}
