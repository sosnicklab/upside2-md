# Ewald Summation Implementation Plan

## Overview
Implement standard Ewald summation for periodic long-range electrostatics across 4 files:
`src/box.h`, `src/box.cpp`, `src/martini.cpp`, `src/main.cpp`, and `example/16.MARTINI/run_sim_bilayer.sh`.

## Step 1: `src/box.h` — Add EwaldState struct and function declarations

Add to `simulation_box::npt` namespace:
- `EwaldSettings` struct: `enabled`, `alpha`, `kmax`, `coulomb_k` (31.775347952181)
- `EwaldState` struct: `settings`, `charges` (per-atom), `kvectors`, `self_energy`
- Function declarations:
  - `void initialize_ewald(hid_t config_root, DerivEngine& engine);`
  - `void apply_ewald_reciprocal(DerivEngine& engine);`
  - `float get_ewald_self_energy(const DerivEngine& engine);`
  - `float get_ewald_reciprocal_energy(const DerivEngine& engine);`

## Step 2: `src/box.cpp` — Implement reciprocal space logic

Add global registry `g_ewald_state` (same pattern as `g_baro_state`).

### `initialize_ewald()`
- Read `ewald_enabled` and `ewald_alpha` from `/input/potential/martini_potential`
- Read per-atom charges from `/input/potential/martini_potential/charges`
- Read box dimensions (same source as barostat)
- Pre-compute self-energy: `V_self = -kC * (alpha/sqrt(pi)) * sum(q_i^2)`
- Pre-generate k-vectors up to kmax (integer triplets where |k| <= kmax)

### `apply_ewald_reciprocal()`
- Get current box dimensions from barostat state (handles NPT box changes)
- For each k-vector (nx, ny, nz) != (0,0,0):
  - Compute k = 2*pi*(nx/Lx, ny/Ly, nz/Lz)
  - Compute k^2 = kx^2 + ky^2 + kz^2
  - Compute structure factor: rho(k) = sum_i q_i * exp(i*k.r_i)
  - Accumulate energy: (kC / (2*pi*V)) * (4*pi^2/k^2) * exp(-k^2/(4*alpha^2)) * |rho(k)|^2
  - Accumulate forces on each atom i: gradient of the above w.r.t. r_i
  - Accumulate virial contribution for pressure tensor
- Store reciprocal energy for logging
- Add forces to `engine.pos->sens` (as negative gradient, matching convention)

## Step 3: `src/martini.cpp` — Screen real-space Coulomb

In `MartiniPotential`:
- Add members: `bool ewald_enabled; float ewald_alpha;`
- In constructor: read `ewald_enabled` and `ewald_alpha` from HDF5 attributes
- In Coulomb spline generation: when `ewald_enabled`, multiply potential by `erfc(alpha*r)`:
  `V_real(r) = kC * qq * erfc(alpha*r) / r`
- In analytical fallback: same screening
- Add per-atom charge array member for Ewald to read

## Step 4: `src/main.cpp` — Hook Ewald into simulation loop

- Add forward declarations for `initialize_ewald` and `apply_ewald_reciprocal`
- Call `initialize_ewald` after `register_barostat_for_engine` during engine setup (~line 779)
- Call `apply_ewald_reciprocal` after `maybe_apply_barostat` in the simulation loop (~line 1128)
- Add loggers for `ewald_reciprocal_energy` and `ewald_self_energy`
- Update potential energy print to include Ewald terms

## Step 5: `example/16.MARTINI/run_sim_bilayer.sh` — Enable via environment

- Add `UPSIDE_EWALD_ENABLE` and `UPSIDE_EWALD_ALPHA` environment variables
- Document in the header comment block

## Step 6: `example/16.MARTINI/prepare_martini.py` — Write HDF5 config

- Add `ewald_enabled` and `ewald_alpha` attributes to martini_potential group
- Read from `UPSIDE_EWALD_ENABLE` and `UPSIDE_EWALD_ALPHA` environment variables

## Key Design Decisions
- Direct k-space summation O(N*Nk) — no FFT dependency
- Reuse barostat's box dimensions (no duplication)
- Thread-safe global registry (same pattern as barostat)
- erfc screening applied at spline generation time (not per-step)
- Self-energy is constant, computed once at initialization
- Reciprocal sum re-evaluated every step (box may change under NPT)
