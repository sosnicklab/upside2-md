#pragma once

#include "deriv_engine.h"
#include <cstdint>
#include <vector>

// Langevin integrator for MARTINI beads: overdamped Brownian with RESPA sub-stepping (default)
// or single-step inertial g-JF (Gronbech-Jensen & Farago). Beads listed in /input/brownian are
// advanced here; all other atoms keep their Verlet dynamics (the Verlet loop must skip them).
namespace martini_brownian {

// Read /input/brownian (atom_index dataset; gamma_ref, t_ref, activation_energy, n_substep
// attributes). No-op if the group is absent.
void register_brownian_for_engine(DerivEngine* engine, hid_t config_root, uint32_t random_seed);

// True if this engine has Brownian atoms.
bool has_brownian(DerivEngine* engine);

// Per-atom mask (length n_atom): 1 where the atom is integrated by Brownian dynamics.
const std::vector<unsigned char>& brownian_mask(DerivEngine* engine);

// Set the target temperature (kT in E_up) used for the Brownian noise.
void set_brownian_temperature(DerivEngine* engine, float temperature);

// Advance the Brownian atoms by one overdamped step of size dt using pos->sens (= dU/dx).
// Call exactly once per integration cycle; the Verlet loop must skip these atoms.
void apply_brownian_step(DerivEngine* engine, VecArray mom, float dt);

}
