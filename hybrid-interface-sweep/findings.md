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

## 2026-04-14 (Downloaded Scalar-Factor Analysis Review)
- The downloaded scalar-factor analysis tree is numerically complete:
  - `15 / 15` analysis tasks succeeded,
  - `0` failed tasks were assembled,
  - all `5` tested scales have full `3 / 3` replicate coverage.
- Internal consistency is strong across the full sweep:
  - every task uses `160` post-burn-in frames,
  - every task uses `102` `PO4` beads,
  - all diffusion values are finite,
  - fit quality stays tight at `R^2 = 0.9882 -> 0.9943`.
- The sweep response is monotonic over the tested range:
  - `interaction_scale = 1.0 -> 0.6`
  - `PO4` diffusion rises from `0.824 -> 1.128 um^2/s`.
- The saved recommendation `interaction_scale = 0.6` is justified by the current selector:
  - it has the highest mean diffusion in-grid,
  - its replicate spread remains modest (`CV = 0.035`),
  - and it also has the best minimum fit quality in-grid (`min R^2 = 0.9915`).
- Relative to the existing `40 ps/step` temperature-specific proxy target:
  - workflow temperature `T = 0.8647` corresponds to about `303.2 K`,
  - linear interpolation between the `0.8` and `0.9` reference rows in `hybrid_timescale.md` implies a target near `2.892 um^2/s`,
  - no tested scale reaches that target.
- Consequence:
  - if the goal is simply “best within this tested scalar grid,” choose `interaction_scale = 0.6`,
  - if the goal is to match the provisional `40 ps/step` physical target, the current grid is still unresolved and should be extended below `0.6`.

## 2026-04-14 (Default Grid Extension)
- The workflow default scale list should no longer stop at `0.6`.
- The next default run should probe substantially lower interaction strength without requiring a manual env override.
- The widened default grid is:
  - `1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1`
- This is an exploration extension, not a claim that the target lies inside that range.

## 2026-04-14 (Downloaded Extended Scalar-Factor Analysis Review)
- The downloaded extended scalar-factor analysis tree is numerically complete:
  - `36 / 36` analysis tasks succeeded,
  - `0` failed tasks were assembled,
  - all `12` tested scales have full `3 / 3` replicate coverage.
- Internal consistency remains strong across the extended sweep:
  - every task uses `160` post-burn-in frames,
  - every task uses `102` `PO4` beads,
  - all diffusion values are finite,
  - fit quality stays tight at `R^2 = 0.9882 -> 0.9967`.
- The response remains broadly monotone toward smaller `interaction_scale`:
  - `1.0 -> 0.1` raises mean diffusion from `0.824 -> 1.598 um^2/s`,
  - there is only a mild local dip at `0.20` versus `0.25`, not a broader instability signature.
- The saved recommendation `interaction_scale = 0.1` is justified by the current selector:
  - it has the highest mean diffusion in-grid,
  - one of the lowest replicate spreads (`CV = 0.020`),
  - and the best minimum fit quality in-grid (`min R^2 = 0.9950`).
- Relative to the same provisional `40 ps/step` target proxy:
  - the workflow target remains about `2.892 um^2/s` at `T = 0.8647`,
  - the best tested point `0.1` reaches `1.598 um^2/s`,
  - so the extended grid still does not achieve the target and reaches only about `55%` of it.
- Consequence:
  - if the goal is “best within the tested scalar grid,” choose `interaction_scale = 0.1`,
  - if the goal is to match the provisional physical target, the study remains unresolved even after the low-end extension.

## 2026-04-14 (Low-End Default Grid Refinement)
- The workflow default scale list should no longer stop at `0.1`.
- Since `0.1` was the best tested point in the extended sweep, the next default run should sample the boundary region below `0.1` more densely.
- The refined default low-end addition is:
  - `0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01`
- This keeps the sweep inside the same scalar-factor semantics and only increases resolution near the current best boundary.

## 2026-04-14 (Downloaded Refined Low-End Scalar-Factor Analysis Review)
- The downloaded refined low-end analysis tree is numerically complete:
  - `63 / 63` analysis tasks succeeded,
  - `0` failed tasks were assembled,
  - all `21` tested scales have full `3 / 3` replicate coverage.
- Internal consistency remains strong across the refined sweep:
  - every task uses `160` post-burn-in frames,
  - every task uses `102` `PO4` beads,
  - all diffusion values are finite,
  - fit quality stays tight at `R^2 = 0.9882 -> 0.9979`.
- The response stays mostly monotone toward smaller `interaction_scale`, with only mild local dips at:
  - `0.25 -> 0.20`,
  - `0.10 -> 0.09`,
  - `0.05 -> 0.04`,
  - `0.02 -> 0.01`.
- The saved recommendation `interaction_scale = 0.02` is justified by the current selector:
  - it has the highest mean diffusion in-grid at `1.8275 um^2/s`,
  - low replicate spread (`CV = 0.015`),
  - and strong fit quality (`min R^2 = 0.9966`).
- The neighboring `0.01` point is statistically indistinguishable from `0.02` within replicate uncertainty, but it does not improve the mean.
  - Because the workflow tie-break rule prefers the larger scale when performance is effectively tied, `0.02` is the conservative best-tested choice for later interface use.
- Relative to the same provisional `40 ps/step` target proxy:
  - the workflow target remains about `2.892 um^2/s` at `T = 0.8647`,
  - the best tested point `0.02` reaches about `63%` of that target,
  - so the refined grid still does not demonstrate a true target match.
- Consequence:
  - if the goal is “best stable scalar factor supported by this bilayer sweep,” choose `interaction_scale = 0.02`,
  - for the actual hybrid run surface, map that to `PROTEIN_ENV_INTERFACE_SCALE = 0.02`,
  - if the goal is a strict physical-timescale match, the calibration remains unresolved.

## 2026-04-14 (Trend-Line Fit On Refined Sweep)
- A simple weighted linear fit across the `21` condition means is a reasonable global trend summary:
  - `D_um2_per_s = -0.9807 * interaction_scale + 1.7549`
  - `R^2 = 0.9707`
- This fit is useful as a sweep-level trend, but it is not a perfect local model of the low-end region:
  - at the best tested point `interaction_scale = 0.02`,
  - the fit predicts `1.7353 um^2/s`,
  - while the observed mean is `1.8275 um^2/s`.
- Extrapolation from the fitted line does not rescue the provisional target:
  - predicted `D(0) = 1.7549 um^2/s`,
  - fitted `interaction_scale` needed for `D = 2.892 um^2/s` is `-1.159`,
  - which lies outside the physically allowed `[0, 1]` range.
- Consequence:
  - the fitted line supports the same scientific conclusion as the raw sweep,
  - the bilayer calibration trend is real,
  - but the current physical-target mapping is still unresolved rather than just hidden between sampled points.

## Lessons
- When the user says a softening factor is “a simple number,” implement a scalar interaction-strength factor, not a softened Hamiltonian shape.
- When the user explicitly says not to modify a reference workflow, keep the fix inside the named target workflow even if the reference path influenced earlier designs.
- When a sweep control changes scientifically, bump manifest/result schemas so stale run directories fail loudly instead of silently reusing old task definitions.
