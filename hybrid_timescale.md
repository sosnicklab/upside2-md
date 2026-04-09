# Hybrid Bilayer Timescale Settings

This file reflects the post-thermostat-fix bilayer rerun. The old strong mass dependence was an artifact of the unit-mass OU noise bug and should not be used anymore.

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

## Post-Fix Scan Summary

- Refreshed `stage7-analysis` validation passed:
  - `525` stage-7 files
  - `525` successful task rows
  - `175` condition rows
  - no tasks below `R^2 < 0.95`
- Within the stable damping window, the corrected scan shows:
  - temperature effect: about `2.87x`
  - mass effect: about `1.36x`
  - damping effect: about `1.05x`
- Interpretation:
  - after fixing the thermostat, temperature is a stronger control than mass inside the tested `mass = 0.1 -> 1.0` range,
  - the current grid does not provide enough leverage to fit a trustworthy bilayer mass correction.

## Closest Stable Settings In The Current Grid

These are the closest stable scanned points to the physical DOPC target at fixed Upside temperature. They are the best currently available settings from the rerun, but they are still much too slow.

| Upside temperature | Kelvin | Stable damping | Best tested mass | Observed diffusion | Sim / target |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0.7` | `245.4 K` | `12` | `0.1` | `0.089 um^2/s` | `0.153` |
| `0.8` | `280.5 K` | `12` | `0.1` | `0.149 um^2/s` | `0.092` |
| `0.9` | `315.5 K` | `12` | `0.1` | `0.201 um^2/s` | `0.056` |
| `1.0` | `350.6 K` | `16` | `0.1` | `0.253 um^2/s` | `0.037` |
| `1.1` | `385.6 K` | `16` | `0.1` | `0.320 um^2/s` | `0.028` |

## What This Means

- No reliable mass correction is resolved by the current post-fix scan.
- Any literal extrapolation from this grid wants a mass far below the tested range, which is not trustworthy enough to save as a calibrated setting.
- The practical post-fix recommendation is:
  - if you need the best currently tested stable settings, use the table above,
  - if you need a quantitative bilayer timescale match, run a new scan with a broader acceleration range instead of fitting a precise correction from the current grid.

## Source Artifacts

- General report: [report.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/report.md)
- Fixed-temperature calibration note: [temperature_specific_correction.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.md)
- Fixed-temperature calibration table: [temperature_specific_corrections.csv](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_corrections.csv)
- Fixed-temperature calibration plot: [temperature_specific_correction.png](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.png)
