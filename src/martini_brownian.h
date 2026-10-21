#pragma once

#include "deriv_engine.h"
#include <cstdint>
#include <vector>

// Inertial g-JF Langevin integrator for a MARTINI system. Every dynamic particle is listed in
// /input/brownian; zero-friction protein carriers reduce exactly to velocity Verlet and
// positive-friction environment particles receive the g-JF bath.
//
// Two attributes on /input/brownian control the discretization and the mobility:
//   inner_steps             N substeps of dt/N per outer step (default 1). The outer dt is what the
//                           friction and the ps-per-step clock are calibrated against, so it stays
//                           fixed; N only makes the integration of that same Langevin equation more
//                           accurate. Needed because mass-1 protein sites at dt = 0.009 pick up a
//                           kinetic-temperature bias against the steep MARTINI core (findings 3.10).
//   reference_temperature_up  temperature the friction was calibrated at, as friction = kT_ref/D_target.
//                           Friction is scaled by T/T_ref at runtime so the target diffusion is
//                           delivered at whatever temperature the replica actually runs. Absent or
//                           zero disables the scaling.
namespace martini_brownian {

// Read /input/brownian (atom_index and absolute friction datasets, the numerical timestep, and the
// optional inner_steps / reference_temperature_up attributes). No-op if the group is absent.
void register_brownian_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);

// True if this engine has g-JF atoms.
bool has_brownian(DerivEngine* engine);

// Per-atom mask for particles that receive stochastic noise in this integrator. A zero-friction
// g-JF particle is thermostatted by the global OU bath instead.
const std::vector<unsigned char>& stochastic_mask(DerivEngine* engine);

// Set the target temperature (kT in E_up) used for the Brownian noise.
void set_brownian_temperature(DerivEngine* engine, float temperature);

// Advance every dynamic atom by one outer g-JF step of size dt using pos->sens (= dU/dx), as
// inner_steps substeps of dt/inner_steps. dt must equal /input/brownian numerical_time_step.
void apply_langevin_step(DerivEngine* engine, VecArray mom, float dt);

}
