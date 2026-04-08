# Hybrid Bilayer Timescale Settings

Warning: these settings were derived before fixing a thermostat bug where momentum noise ignored per-atom mass. They are provisional only and should be regenerated after rerunning the bilayer diffusion scan with the corrected thermostat.

These settings keep the scanned Upside temperature fixed and apply a bilayer-only mass / damping correction to match a physical DOPC diffusion target as closely as possible.

## Unit Contract

- Upside temperature conversion from [AGENTS.md](/Users/yinhan/Documents/upside2-md-martini/AGENTS.md#L31): `1.0 T_up = 350.588235 K`
- Physical DOPC target curve from:
  - Bag N, Yap DHX, Wohland T. *Temperature dependence of diffusion in model and live cell membranes characterized by imaging fluorescence correlation spectroscopy.* Biochim Biophys Acta Biomembr. 2014;1838(3):802-813.
  - DOI: `10.1016/j.bbamem.2013.10.009`

## Recommended Bilayer Corrections

| Upside temperature | Kelvin | Recommended damping | Recommended mass scale | Closest stable scanned point |
| --- | ---: | ---: | ---: | --- |
| `0.7` | `245.4 K` | `16` | `0.0676` | `mass=0.1, tau=4, D=0.311 um^2/s` |
| `0.8` | `280.5 K` | `4` | `0.0419` | `mass=0.1, tau=12, D=0.374 um^2/s` |
| `0.9` | `315.5 K` | `6` | `0.0324` | `mass=0.1, tau=16, D=0.454 um^2/s` |
| `1.0` | `350.6 K` | `8` | `0.0303` | `mass=0.1, tau=8, D=0.507 um^2/s` |
| `1.1` | `385.6 K` | `6` | `0.0296` | `mass=0.1, tau=4, D=0.758 um^2/s` |

## Interpretation

- Damping is a secondary knob in this scan. Mass scaling is the main correction.
- Every temperature requires `mass < 0.1`, so the current scan minimum mass is still too heavy.
- A practical next scan would keep temperature fixed, keep damping in a stable band around `tau = 6-8`, and extend the bilayer mass scale below `0.1`.
- The source experiment directly constrains `298-313 K`. The rows for `T_up = 0.7`, `0.8`, `1.0`, and `1.1` are therefore extrapolations under the Upside temperature mapping.

## Source Artifacts

- Detailed calibration table: [temperature_specific_corrections.csv](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_corrections.csv)
- Detailed calibration note: [temperature_specific_correction.md](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.md)
- Calibration plot: [temperature_specific_correction.png](/Users/yinhan/Documents/upside2-md-martini/bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.png)
