# Hybrid Bilayer Timescale Settings

This file reflects the refreshed boundary-focused dryMARTINI bilayer scan. It supersedes the earlier broad joint-scan summary.

## Hybrid Interface Scale Calibration

This file also records the separate protein-interface calibration that was used to choose the hybrid `1rkl` production scale.

- Goal:
  - choose the Upside/dryMARTINI interface factor that reproduces the protein fluctuation amplitude of the implicit membrane reference
- Reference method:
  - `example/08.MembraneSimulation`
  - `1rkl`
  - fixed-curvature implicit membrane
  - `membrane_thickness = 24.8`
  - `curvature_radius = 120.0`
  - `curvature_sign = 1`
  - `temperature = 0.80`
- Hybrid method:
  - `example/16.MARTINI/run_sim_1rkl.sh`
  - sweep / production control:
    - `PROTEIN_ENV_INTERFACE_SCALE`

### Calibration Method

- Compare protein-only fluctuations, not the full hybrid protein+bilayer dump.
- Use only the membrane-embedded protein region derived from the reference membrane geometry.
- Embedded reference region used in the latest rerun:
  - residues `7-24`
- Primary metric:
  - embedded-region RMSD amplitude difference to the reference
- Primary presentation artifact:
  - direct reference-vs-hybrid RMSF overlays for each stable sampled scale
- Secondary artifact only:
  - fitted trend line versus `interface_scale`

### What We Ran

- Replaced the older broad exploratory default with a focused high-end sweep in `hybrid-interface-sweep/`.
- Lengthened the default trajectories to reduce RMSF noise:
  - reference `REFERENCE_DURATION = 200001`
  - hybrid `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1000`
  - hybrid `PROD_70_NSTEPS = 50000`
- Raised the default Slurm run walltime to:
  - `36:00:00`
- Switched sweep initialization to fresh random seeds by default.
- Filtered destroyed trajectories from the analysis:
  - explicit non-finite backbone-frame failures are excluded
  - structurally unstable protein trajectories are excluded by reference-relative RMSF / geometry thresholds

### Latest Rerun Result

- Downloaded rerun coverage:
  - `20 / 20` run tasks succeeded
  - `19 / 20` analysis tasks succeeded
  - `1` analysis task failed:
    - `scale1p1_r01`
    - `Too few finite backbone frames remain after filtering for RMSF analysis: 0`
- Stable sampled scales that survived the filter:
  - `0.825`
  - `0.925`
  - `1.0`
  - `1.075`
  - `1.15`
  - `1.20`
- Best sampled stable scale:
  - `interface_scale = 1.15`
  - embedded-region RMSD delta `= 0.05765 A`
- Fitted trend line on the same rerun:
  - recommended `0.915`
  - `R^2 = 0.268`
- Decision:
  - the fit is too weak to use as the main argument
  - the direct overlay comparison was used instead
  - chosen production interface scale:
    - `PROTEIN_ENV_INTERFACE_SCALE = 1.15`

### Current Project Default

- `hybrid-interface-sweep/workflow.py` now defaults to:
  - `interface_scale = 1.15`
- `example/16.MARTINI/run_sim_1rkl.sh` now defaults to:
  - `PROTEIN_ENV_INTERFACE_SCALE = 1.15`

## Unit Contract

- Upside temperature conversion from [AGENTS.md](/Users/yinhan/Documents/upside2-md-martini/AGENTS.md#L31): `1.0 T_up = 350.588235 K`
- Upside timestep calibration used for the downloaded analysis: `40 ps` per integration step
- Workflow integration step used in the bilayer scan: `0.01`
- Therefore:
  - `1` workflow time unit = `4 ns`
  - `nm^2/ns = diffusion_mean_nm2_per_time / 4`
  - `um^2/s = diffusion_mean_nm2_per_time * 250`
- Caveat:
  - workflow `16` now matches the standard example thermostat default again,
  - but its hybrid stages still use smaller MD timesteps than the examples behind the borrowed `40 ps/step` calibration,
  - so the absolute mapping should still be treated as provisional until it is validated directly for the hybrid stage-7 setup.

## Boundary Scan

- Downloaded `stage7-analysis` counts:
  - `630` stage-7 files
  - `630` successful task rows
  - `210` condition rows
  - `4` tasks below `R^2 < 0.95`
- Scan grid:
  - `temperature = 0.7, 0.8, 0.9, 1.0, 1.1`
  - `tau = 5, 8, 10, 12, 14, 16`
  - `mass = 0.02, 0.015, 0.01, 0.0075, 0.005, 0.003, 0.002`
- Raw-grid caveats:
  - `16` task rows have `NaN` diffusion outputs
  - `9` condition rows have non-finite diffusion
  - `70` conditions exceed `CV > 0.25`

## Robust Subset

The saved recommendations use a robust subset only:

- completed replicates only
- `R^2 >= 0.95`
- `CV <= 0.25`
- continuity-stable mass progression, with the first `> 4x` adjacent-mass jump treated as the start of an unstable lighter-mass branch

That leaves `125 / 210` usable conditions.

Within that robust subset:

- temperature effect: `3.02x`
- damping effect: `1.33x`
- mass effect: `1.30x`
- best robust point: `T = 1.1`, `tau = 5`, `mass = 0.005`, `D = 0.394 um^2/s`
- no robust condition reaches `1.0 um^2/s`

## What The Boundary Scan Means

- This scan improved the stability picture:
  - the robust subset grew from `49` to `125` conditions,
  - the raw-grid NaN count dropped sharply compared with the earlier joint scan.
- But it did not improve the calibration outcome:
  - the best robust diffusion stayed below `0.4 um^2/s`,
  - the best robust settings are still at `tau = 5`,
  - higher damping timescales `10-16` do not beat the `tau = 5` branch in the trustworthy subset.

So the boundary scan says:

- `tau = 5` remains the best damping anchor
- `tau = 8` is still usable but slower on average
- `tau = 10-16` mostly traces the approach to the instability frontier rather than a better production regime
- changing mass inside `0.02 -> 0.002` still gives only modest gain inside the trustworthy region

## Instability Frontier

The full grid does reach much faster diffusion, but only on excluded branches. The onset is now clearer:

- mostly `tau = 12-16`
- most obvious at `mass = 0.003` and `0.002`
- characterized by abrupt mass-to-mass jumps, high CV, or non-finite outputs

Examples from the excluded frontier include:

- `T = 0.7`, `tau = 16`, `mass = 0.002`, `D = 121.8 um^2/s`, `CV = 0.45`
- `T = 0.8`, `tau = 16`, `mass = 0.002`, `D = 114.6 um^2/s`, continuity-failed
- `T = 1.1`, `tau = 16`, `mass = 0.003`, `D = 53.4 um^2/s`, `CV = 0.80`

These are instability-frontier points, not calibration candidates.

## Temperature-Specific Robust Summary

Using the physical DOPC target curve as the timescale proxy, the nearest robust in-grid points are still all too slow:

| Upside temperature | Kelvin | target D | nearest robust `(mass, tau, D)` | Result |
| --- | ---: | ---: | --- | --- |
| `0.7` | `245.4 K` | `0.584 um^2/s` | `(0.002, 5, 0.154)` | unresolved |
| `0.8` | `280.5 K` | `1.621 um^2/s` | `(0.002, 5, 0.254)` | unresolved |
| `0.9` | `315.5 K` | `3.585 um^2/s` | `(0.0075, 5, 0.297)` | unresolved |
| `1.0` | `350.6 K` | `6.765 um^2/s` | `(0.002, 5, 0.373)` | unresolved |
| `1.1` | `385.6 K` | `11.374 um^2/s` | `(0.005, 5, 0.394)` | unresolved |

`Unresolved` here means the robust subset never reaches the target at that temperature.

## Practical Recommendation

- If you need a defensible stable bilayer setting now, stay at `tau = 5`.
- The best current robust mass window is roughly `0.002 -> 0.0075`, depending on temperature.
- Do not promote `tau = 10-16` into production as a speedup mechanism. This scan says those values mostly move you toward the noisy frontier, not toward a better robust match.
- The current boundary scan therefore strengthens the earlier conclusion:
  - matching dryMARTINI bilayer timescale to Upside is not being solved by damping retuning,
  - mass retuning in the currently tested range is still too weak,
  - a trustworthy timescale match is still not demonstrated.

## Source Artifacts

- General report: [report.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/report.md)
- Recommendation summary: [recommendation_summary.json](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/recommendation_summary.json)
- Robust condition table: [stable_condition_summary_physical_units.csv](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/stable_condition_summary_physical_units.csv)
- Temperature-specific correction note: [temperature_specific_correction.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.md)
- Temperature-specific correction table: [temperature_specific_corrections.csv](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_corrections.csv)
