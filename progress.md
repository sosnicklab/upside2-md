# Progress Log

## Current Task: Persisted CGL Output Geometry Issue

- Actions taken:
  - Inspected the named 1AFO/1RKL H5/VTF outputs and confirmed salts are
    present in the current prep summaries.
  - Found that VTF rod mode was moving the original CGL atom to a PO4 display
    endpoint instead of preserving the physical CGL center.
  - Found that 1RKL CGL center outliers are introduced during stage-6
    minimization and inherited by stage 7.
  - Fixed rod-mode extraction so the physical CGL center atom is preserved in
    VTF output.
  - Reverted the mass-default and CGL z-restraint changes after user correction:
    they suppress CGL z motion by destroying the calibrated timescale and are
    invalid fixes.
  - Identified the remaining 1RKL minimizer ejection as seed/basin-sensitive:
    seed `2026` fails, while `11`, `22`, and `2027` pass under calibrated CGL
    dynamics.
  - Added bounded coarse-CGL stage-6 geometry rejection/retry that regenerates
    stage-0 packing and reruns stage 6 with incremented seeds.
  - Confirmed no leftover unneeded diagnostic scripts remain in
    `example/16.MARTINI`.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `findings.md`
  - `py/martini_extract_vtf.py`
  - `py/martini_prepare_system.py`
- Verification:
  - Python syntax checks passed for `py/martini_prepare_system.py`,
    `py/martini_prepare_system_lib.py`, and `py/martini_extract_vtf.py`.
  - Shell syntax checks passed for `run_sim_hybrid.sh`, `run_sim_1rkl.sh`,
    `run_sim_1afo.sh`, `run_sim_1rkl_full.sh`, and `run_sim_1afo_full.sh`.
  - The shortened geometry run that relied on mass/restraint changes is invalid
    evidence and should not be used.
  - 1RKL retry test with bad seeds `2026/2026` rejected the first stage-6
    geometry and passed after retrying with `2027/2027`; final CGL mass stayed
    at calibrated `1.008` Upside units and wrapped CGL z outliers were `0`.
  - 1AFO retry-enabled test with `2026/2026` passed without retry; final CGL
    mass stayed at `1.008` and wrapped CGL z outliers were `0`.
