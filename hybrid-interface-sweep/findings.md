# Findings

## 2026-04-13 (User Correction)
- The current `hybrid-interface-sweep` is scientifically mis-scoped for the stated goal.
  - Sweeping `PROTEIN_ENV_INTERFACE_SCALE` through the full hybrid protein+bilayer workflow answers a protein/bilayer interaction-strength question.
  - It does not answer the requested bilayer-only softening-calibration question.
- The corrected workflow target is:
  - run bilayer-only dry-MARTINI tests,
  - vary the softening applied to Lennard-Jones and Coulomb interactions,
  - measure bilayer lateral diffusion and a viscosity proxy,
  - then use that calibration to choose softening intended for later protein/bilayer interface use.

## 2026-04-13 (Implementation Facts)
- `py/martini_prepare_system.py` already supports `--mode bilayer`, so the corrected workflow can prepare bilayer-only stage inputs without inventing a new packer.
- Production softening is represented directly on `input/potential/martini_potential` by these attrs:
  - `lj_soften`
  - `lj_soften_alpha`
  - `coulomb_soften`
  - `slater_alpha`
- `src/martini.cpp` reads those attrs at runtime and uses them to generate the softened LJ and Coulomb spline tables for the active MARTINI Hamiltonian.
- The existing repo context does not provide a ready-made bilayer viscosity observable, so the workflow should report diffusion plus a clearly labeled reciprocal-diffusion viscosity proxy instead of pretending to compute a direct viscosity.

## Lessons
- When the user says a sweep should calibrate a future interface control through bilayer-only tests, do not route the workflow through the actual hybrid protein+bilayer simulation path.
- When the user asks for softening calibration, distinguish interaction-strength scaling from Hamiltonian softening parameters before writing the sweep surface.
- If a workflow’s scientific observable changes, rewrite the tracker files first so implementation and verification are aligned with the corrected question.
