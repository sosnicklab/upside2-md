# Task Plan

## Project Goal
- Remove the hybrid RMSD-align feature from the MARTINI runtime and workflow, then verify the hybrid example still runs without it.

## Architecture & Key Decisions
- Delete the runtime bookkeeping path entirely instead of keeping a dormant flag.
- Remove the corresponding workflow attrs and shell knob so the example no longer carries dead RMSD-align configuration.
- Do not introduce any replacement physics; this change only removes the unused whole-protein RMSD-align bookkeeping path.
- Verify with a reduced end-to-end hybrid workflow run after rebuilding `obj/upside`.
- Align workflow `16` thermostat defaults with the standard Upside examples when comparing against repo-wide time calibrations.
- Keep workflow `16`'s smaller hybrid timesteps for stability, but document that the borrowed `40 ps/step` calibration from the standard examples is not automatically transferable to stage `7.0`.

## Execution Phases
- [x] Phase 1: Recreate task tracker files for the current MARTINI example directory and audit the current RMSD-align path.
- [x] Phase 2: Remove the RMSD-align runtime state, parser fields, helper function, and integrator call sites from C++.
- [x] Phase 3: Remove stale workflow-level RMSD-align knobs and metadata writes.
- [x] Phase 4: Rebuild and run a reduced hybrid workflow verification.
- [x] Phase 5: Audit thermostat defaults across the examples and align workflow `16` with the standard example thermostat setting.

## Known Errors / Blockers
- No blocker identified.

## Review
- Removed the feature from C++ runtime code:
  - deleted the `align_active_protein_coordinates(...)` path and its five integrator call sites,
  - removed the RMSD-align parser fields, runtime state, previous-frame storage, and banner output from `src/martini.cpp`.
- Removed stale workflow plumbing:
  - deleted the `INTEGRATION_RMSD_ALIGN_ENABLE` shell knob,
  - removed `integration_rmsd_align_enable`, `coupling_align_debug`, and `coupling_align_interval` metadata writes,
  - simplified production-stage hybrid-control updates to only set active knobs that still exist.
- Verification:
  - `rg` over `src`, `py/martini_prepare_system_lib.py`, and `example/16.MARTINI/run_sim_1rkl.sh` found no remaining code references to the RMSD-align feature,
  - `python3 -m py_compile py/martini_prepare_system_lib.py`,
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`,
  - `cmake --build obj -j4`,
  - reduced end-to-end workflow run in `outputs/martini_test_1rkl_noalign_removed` completed through stage `7.0`.
- Thermostat audit:
  - standard examples such as `example/01.GettingStarted/0.run.py` and `example/08.MembraneSimulation/0.normal.run.py` do not pass `--thermostat-timescale`, so they inherit the engine default `tau = 5.0`,
  - workflow `16` had been the outlier with `THERMOSTAT_TIMESCALE=4.0`,
  - workflow `16` should therefore default to `THERMOSTAT_TIMESCALE=5.0` if it is meant to match the standard example damping.
- Scope limit:
  - this confirms the feature can be removed cleanly without breaking the example workflow,
  - it does not establish long-timescale physics equivalence.
  - workflow `16` still uses `EQ_TIME_STEP=0.010` and `PROD_TIME_STEP=0.002`, so the repo-wide `40 ps/step` calibration borrowed from the default examples still needs separate validation before using it literally for hybrid stage `7.0`.
