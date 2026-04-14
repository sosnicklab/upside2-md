# Findings

## 2026-04-13 (Current Workflow Target)
- `hybrid-interface-sweep/` must be a bilayer-only calibration workflow.
- The factor being calibrated is a simple scalar strength multiplier `f`, not a softened-potential shape parameter.
- The requested change belongs only in `hybrid-interface-sweep/`.
  - `example/16.MARTINI/` is not part of this correction.

## 2026-04-13 (Control Semantics)
- The scientifically correct sweep surface for this workflow is one shared nonbonded interaction factor.
- For the staged MARTINI coefficient table this is implemented as:
  - `epsilon *= f` for LJ,
  - `q -> sqrt(f) * q` for Coulomb.
- This is equivalent to multiplying evaluated LJ and Coulomb pair energies and forces by `f` at runtime.
- Consequence:
  - `lj_soften_alpha` and `slater_alpha` are the wrong sweep variables for this study,
  - the earlier downloaded `analysis/` tree is not valid for selecting the requested factor.

## 2026-04-13 (Stage-7 Rewrite Surface)
- The bilayer prep path already writes all required nonbonded information into `input/potential/martini_potential`.
- Runtime pair evaluation uses the `coefficients` table entries:
  - column `0`: LJ `epsilon`,
  - column `1`: LJ `sigma`,
  - columns `2` and `3`: per-pair Coulomb charges.
- Rewriting only the staged `7.0` coefficient table is sufficient to change the production nonbonded forces seen by the bilayer-only diffusion measurement.
- The corrected workflow should also force:
  - `lj_soften = 0`
  - `coulomb_soften = 0`
  on the production-stage file so the scalar factor is not mixed with soft-core/Slater shaping.

## 2026-04-13 (Analysis Surface)
- The existing bilayer-only analysis path remains the right observable family for this calibration:
  - `PO4` lateral diffusion relative to bilayer COM,
  - physical-unit diffusion conversion via the configured `ps/step` assumption,
  - reciprocal-diffusion viscosity proxy,
  - optional target-diffusion ranking.
- Only the grouping axis and recommendation surface need to change:
  - from `(lj_alpha, slater_alpha)`
  - to `interaction_scale`.

## 2026-04-13 (Verification Facts)
- The corrected scalar-factor workflow was smoke-tested at `interaction_scale = 0.85` under `/tmp/hybrid_interface_scale_smoke_v2`.
- Direct comparison of `bilayer.stage_7.0.prepared.up` and `bilayer.stage_7.0.up` confirmed:
  - `epsilon` ratios are exactly `0.85` within float32 precision,
  - both charge columns scale by `0.921954`, matching `sqrt(0.85)`,
  - `lj_soften = 0`,
  - `coulomb_soften = 0`,
  - provenance attrs record `interaction_scale` and `coulomb_charge_scale`.
- The corrected analysis path assembles cleanly on the smoke tree with `interaction_scale` as the only condition axis.

## Lessons
- When the user says a softening factor is “a simple number,” implement a scalar interaction-strength factor, not a softened Hamiltonian shape.
- When the user explicitly says not to modify a reference workflow, keep the fix inside the named target workflow even if the reference path influenced earlier designs.
- When a sweep control changes scientifically, bump manifest/result schemas so stale run directories fail loudly instead of silently reusing old task definitions.
