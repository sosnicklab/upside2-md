# Hybrid Bilayer Timescale Settings

This file reflects the focused post-fix mass-only bilayer rerun at fixed `tau = 5`. The old strong mass dependence from the pre-fix scan was an artifact of the unit-mass OU noise bug and should not be used anymore.

## Unit Contract

- Upside temperature conversion from [AGENTS.md](/Users/yinhan/Documents/upside2-md-martini/AGENTS.md#L31): `1.0 T_up = 350.588235 K`
- Upside timestep calibration used for the downloaded analysis: `40 ps` per integration step
- Caveat:
  - the repo examples that motivate that `40 ps` calibration use the standard Upside defaults, including the default thermostat timescale and default `dt = 0.009`,
  - workflow `16` now matches the standard example thermostat default again, but its hybrid stages still run at `EQ_TIME_STEP=0.010` and `PROD_TIME_STEP=0.002`,
  - so the borrowed `40 ps/step` calibration should still be treated as provisional for hybrid stage `7.0` until it is validated for that smaller hybrid timestep.
- Workflow integration step used in the bilayer scan: `0.01`
- Therefore:
  - `1` workflow time unit = `4 ns`
  - `nm^2/ns = diffusion_mean_nm2_per_time / 4`
  - `um^2/s = diffusion_mean_nm2_per_time * 250`

## Focused Tau-5 Scan Summary

- Refreshed `stage7-analysis` validation passed:
  - `90` stage-7 files
  - `90` successful task rows
  - `30` condition rows
  - no tasks below `R^2 < 0.95`
- This rerun kept damping fixed at `tau = 5` and only pushed mass lower:
  - `mass = 0.1, 0.05, 0.02, 0.01, 0.005, 0.002`
  - `temperature = 0.7, 0.8, 0.9, 1.0, 1.1`
- In this focused grid:
  - temperature effect: about `3.20x`
  - mass effect: about `1.61x`
  - damping effect: fixed at `1.00x` because `tau` was held constant
- Interpretation:
  - lowering mass does accelerate the bilayer, but temperature still changes diffusion more strongly than mass across this tested range,
  - even at `mass = 0.002`, the bilayer remains far below the physical DOPC target at the same temperature.

## Best Robust Settings In The Current Grid

These are the best low-variance (`CV <= 0.25`) points currently available at fixed `tau = 5`. They are the most defensible operational choices from the rerun, but they are still much too slow for a physical DOPC match.

| Upside temperature | Kelvin | Fixed damping | Best robust mass | Observed diffusion | Notes |
| --- | ---: | ---: | ---: | ---: | --- |
| `0.7` | `245.4 K` | `5` | `0.002` | `0.120 um^2/s` | lowest mass is usable but only modestly faster |
| `0.8` | `280.5 K` | `5` | `0.002` | `0.225 um^2/s` | lowest mass is usable |
| `0.9` | `315.5 K` | `5` | `0.005` | `0.284 um^2/s` | `0.002` is noisier and should not anchor calibration |
| `1.0` | `350.6 K` | `5` | `0.010` | `0.311 um^2/s` | best robust point at this temperature |
| `1.1` | `385.6 K` | `5` | `0.005` | `0.365 um^2/s` | best robust point in the whole grid |

## Closest Stable Points To The Physical Target

Using the temperature-specific calibration report, the closest in-grid points to the physical DOPC target are still far too slow:

| Upside temperature | Kelvin | Fixed damping | Nearest mass | Observed diffusion | Sim / target |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0.7` | `245.4 K` | `5` | `0.002` | `0.120 um^2/s` | `0.206` |
| `0.8` | `280.5 K` | `5` | `0.002` | `0.225 um^2/s` | `0.139` |
| `0.9` | `315.5 K` | `5` | `0.002` | `0.410 um^2/s` | `0.114` |
| `1.0` | `350.6 K` | `5` | `0.010` | `0.311 um^2/s` | `0.046` |
| `1.1` | `385.6 K` | `5` | `0.002` | `0.432 um^2/s` | `0.038` |

## What This Means

- The focused mass-only rerun improved bilayer fluidity, but it still did not resolve a reliable timescale match.
- No reportable mass correction is resolved by the current grid; the calibration file still marks every temperature as `unresolved (< 0.01)`.
- The practical recommendation is now:
  - if you need the best currently tested stable settings at fixed Upside damping, use the `Best Robust Settings` table above,
  - if you need a quantitative bilayer timescale match while keeping `tau = 5`, you will need either even lighter effective bilayer masses or a different acceleration model.

## Source Artifacts

- General report: [report.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/report.md)
- Fixed-temperature calibration note: [temperature_specific_correction.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.md)
- Fixed-temperature calibration table: [temperature_specific_corrections.csv](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_corrections.csv)
- Fixed-temperature calibration plot: [temperature_specific_correction.png](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.png)
