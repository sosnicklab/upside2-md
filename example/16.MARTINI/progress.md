# Progress Log

## 2026-04-06 (Audit: Does Upside RMSD-Align Protein Backbone Each Step?)
- Recreated the required tracker files in `example/16.MARTINI` because they had been removed during prior directory cleanup.
- Audited runtime references related to RMSD alignment across:
  - `src/martini.cpp`
  - `src/deriv_engine.cpp`
  - `example/16.MARTINI/run_sim_1rkl.sh`
  - supporting config code in `py/martini_prepare_system_lib.py`
- Confirmed workflow wiring:
  - `INTEGRATION_RMSD_ALIGN_ENABLE` defaults to `1` in `run_sim_1rkl.sh`,
  - production stage setup writes `integration_rmsd_align_enable` into `/input/hybrid_control`.
- Confirmed runtime behavior:
  - `martini_hybrid::align_active_protein_coordinates(...)` is called repeatedly from the integrator before derivative evaluations,
  - but the function only computes/stores RMSD-alignment bookkeeping and optional debug output,
  - it does not modify integrated coordinates or momenta.
- Verification:
  - `rg -n "rmsd|align|coupling_align|integration_rmsd_align|affine_alignment" src py example/16.MARTINI/run_sim_1rkl.sh`
  - inspected `src/martini.cpp` around `integration_rmsd_align_enable` parsing and `align_active_protein_coordinates(...)`
  - inspected `src/deriv_engine.cpp` call sites for `martini_hybrid::align_active_protein_coordinates(...)`
  - inspected `example/16.MARTINI/run_sim_1rkl.sh` production-control wiring

## 2026-04-06 (Audit: Would Removing Whole-Protein Rigid-Body Momentum Preserve Results?)
- Extended the runtime audit to inspect likely pure-Upside rolling sources and the effect of possible fixes.
- Inspected:
  - `src/thermostat.cpp` and `src/thermostat.h` for OU momentum injection,
  - `src/main.cpp` for thermostat application and recentering hooks,
  - `src/deriv_engine.cpp` for position recentering behavior.
- Findings:
  - OU thermostat applies atomwise stochastic momentum updates with no subsequent rigid-body angular-momentum removal.
  - Upside does have position recentering, but only translational recentering of coordinates, not momentum or angular-velocity removal.
  - Whole-protein angular-momentum removal would change rotational dynamics, so it is not a no-change substitute for the current bookkeeping-only RMSD alignment hook.
- Conclusion recorded for decision-making:
  - if the goal is no change to Upside results, do not default to angular-momentum removal;
  - diagnose net force / torque / linear momentum / angular momentum first to separate force bugs from expected rigid-body thermostat motion.

## 2026-04-06 (Rollback Test: Hybrid RMSD-Align Disabled)
- First changed `example/16.MARTINI/run_sim_1rkl.sh` so `INTEGRATION_RMSD_ALIGN_ENABLE` defaults to `0`.
- Ran an initial reduced workflow smoke test and found stage `6.0` still reported `integration_rmsd_align=1`.
- Traced the issue to generated hybrid-control metadata:
  - `py/martini_prepare_system_lib.py::write_hybrid_mapping_h5(...)` created `/input/hybrid_control` without `integration_rmsd_align_enable`,
  - `src/martini.cpp` therefore fell back to runtime default `1`.
- Fixed the workflow-level source of truth by adding `ctrl.attrs["integration_rmsd_align_enable"] = np.int8(0)` in `py/martini_prepare_system_lib.py`.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system_lib.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - reduced end-to-end run:
    - `RUN_DIR=outputs/martini_test_1rkl_noalign_smoke2`
    - `MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1`
    - `EQ_62_NSTEPS=0 EQ_63_NSTEPS=0 EQ_64_NSTEPS=0 EQ_65_NSTEPS=0 EQ_66_NSTEPS=0`
    - `PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1`
- Observed result:
  - stages `6.0` through `7.0` completed successfully,
  - runtime banners for stages `6.0`, `6.1`, `6.2`, `6.3`, `6.4`, `6.5`, `6.6`, and `7.0` all reported `integration_rmsd_align=0`,
  - workflow finished and wrote the expected checkpoint and VTF outputs under `outputs/martini_test_1rkl_noalign_smoke2`.

## 2026-04-06 (Implementation: Remove Hybrid RMSD-Align from Runtime and Workflow)
- Deleted the C++ runtime feature:
  - removed `martini_hybrid::align_active_protein_coordinates(...)` from `src/martini.cpp`,
  - removed the associated Kabsch/Horn quaternion helpers that were only used by that path,
  - removed five integrator call sites plus the forward declaration from `src/deriv_engine.cpp`,
  - removed the RMSD-align fields from `HybridRuntimeState`, hybrid-control parsing, stage-reset logic, and runtime banner text in `src/martini.cpp`.
- Deleted stale workflow plumbing:
  - removed `INTEGRATION_RMSD_ALIGN_ENABLE` from `example/16.MARTINI/run_sim_1rkl.sh`,
  - simplified `set_hybrid_production_controls()` to update only `production_nonprotein_hard_sphere`,
  - removed `integration_rmsd_align_enable`, `coupling_align_debug`, and `coupling_align_interval` attrs from `py/martini_prepare_system_lib.py`.
- Verification:
  - `rg -n "integration_rmsd_align|align_active_protein_coordinates|prev_bb_pos_rmsd|has_prev_bb_rmsd|coupling_align_debug|coupling_align_interval|integration_align_step" src example/16.MARTINI/run_sim_1rkl.sh py/martini_prepare_system_lib.py`
  - `python3 -m py_compile py/martini_prepare_system_lib.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `cmake --build obj -j4`
  - reduced end-to-end workflow run:
    - `RUN_DIR=outputs/martini_test_1rkl_noalign_removed`
    - `MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1`
    - `EQ_62_NSTEPS=0 EQ_63_NSTEPS=0 EQ_64_NSTEPS=0 EQ_65_NSTEPS=0 EQ_66_NSTEPS=0`
    - `PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1`
- Observed result:
  - rebuild succeeded,
  - runtime banners no longer print any RMSD-align field,
  - stages `6.0` through `7.0` completed successfully under `outputs/martini_test_1rkl_noalign_removed`.

## 2026-04-09 (Thermostat Consistency Audit)
- Audited thermostat usage across the example workflows:
  - `example/01.GettingStarted/0.run.py`,
  - `example/08.MembraneSimulation/0.normal.run.py`,
  - `example/12.MultistepIntegrator/0.run.py`,
  - `example/16.MARTINI/run_sim_1rkl.sh`,
  - CLI defaults in `src/main.cpp`.
- Confirmed that the standard examples inherit the engine default `--thermostat-timescale=5.0`, while workflow `16` was the outlier with `THERMOSTAT_TIMESCALE=4.0`.
- Updated `example/16.MARTINI/run_sim_1rkl.sh` so workflow `16` now defaults to `THERMOSTAT_TIMESCALE=5.0`.
- Updated `example/16.MARTINI/readme.md` to document:
  - the thermostat-default alignment with the other examples,
  - the remaining caveat that workflow `16` still uses smaller hybrid timesteps, so the borrowed `40 ps/step` calibration still needs separate validation for stage `7.0`.

## 2026-04-14 (Packing-First Calibration Tooling)
- Opened a new local task to address membrane packing before further timescale tuning and wrote `example/16.MARTINI/plan.md`.
- Audited the current hybrid workflow and confirmed:
  - equilibration stages `6.0 -> 6.6` already run with NPT enabled,
  - production stage `7.0` has NPT plumbing but defaults to `PROD_70_NPT_ENABLE=0`,
  - stage outputs contain the box and lipid metadata needed for direct packing analysis.
- Added `example/16.MARTINI/analyze_packing.py`:
  - reads completed `stage_7.0.up` files,
  - computes XY area,
  - computes leaflet-resolved and mean APL from `PO4`-identified lipid molecules,
  - computes `PO4` leaflet thickness,
  - computes tail-bond orientational order relative to the bilayer normal,
  - writes optional JSON summary and per-frame CSV outputs.
- Updated `example/16.MARTINI/run.py` to expose:
  - `PROD_70_NPT_ENABLE`
  - `PROD_70_BAROSTAT_TYPE`
  in the top-level wrapper settings.
- Updated `example/16.MARTINI/readme.md` with:
  - production-NPT guidance for packing tests,
  - usage for `analyze_packing.py`,
  - a packing-first calibration loop.
- Verification:
  - `.venv/bin/python -m py_compile example/16.MARTINI/analyze_packing.py example/16.MARTINI/run.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `.venv/bin/python example/16.MARTINI/analyze_packing.py example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up --json-out /tmp/packing_summary.json --csv-out /tmp/packing_timeseries.csv`
- Observed baseline on the current `martini_test_1rkl_hybrid` output:
  - `APL = 89.01 A^2`
  - `PO4 thickness = 32.57 A`
  - `tail order = 0.238`
  - `n_lipid_molecules = 280`
  - no dynamic `output/box` frames were present, so area stayed fixed under the expected NVT fallback path.
