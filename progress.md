# Progress Log

## 2026-06-16 1RKL MARTINI Secondary-Structure Regression

- Actions taken:
  - Compared uncommitted workflow/runtime changes against `HEAD`.
  - Read the generated stage-7 logs for both reported outputs.
  - Quantified protein H-bond/Rg/RMSD metrics from the current CGL and full-lipid stage-7 HDF5 files.
  - Ran controlled copied-checkpoint continuations at `T=0.8647` and `T=1.2` with the same seed and current integrator settings.
  - Restored the wrapper default temperature to the committed `0.8647`.
  - Changed hybrid verbose progress logging to print actual MD step counts.
- Files modified:
  - `example/16.MARTINI/run_sim_hybrid.sh`
  - `src/main.cpp`
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - Bad generated logs both used `--temperature 1.2`.
  - Bad CGL output: hbond sum `23.92 -> 14.28`, CA Rg `12.73 -> 14.02 A`.
  - Bad full-lipid output: hbond sum `26.38 -> 15.45`, CA Rg `12.69 -> 13.73 A`.
  - Controlled `T=0.8647` continuation: hbond final `33.48`, last-20 mean `32.15`, CA RMSD max `1.03 A`.
  - Controlled `T=1.2` continuation: hbond final `22.30`, last-20 mean `21.55`, CA RMSD max `1.34 A`.
  - `make -C obj -j4` passed after the logging change; remaining warnings are pre-existing random/MonteCarlo warnings.
  - Copied 120-step smoke run printed `step 0 / 120`, `step 50 / 120`, and `step 100 / 120`.
- Residual notes:
  - The existing reported VTF/HDF5 outputs were generated before the fix and should be regenerated if those exact artifact paths need to be clean.

## 2026-06-16 CGL-Only Lateral Diffusion Timescale

- Actions taken:
  - Replaced hidden `3x` MARTINI/Upside time accounting with explicit
    `--integrator mv --inner-step 1` MD execution in the main hybrid runner and
    CGL calibration/validation helpers.
  - Updated `analyze_lipid_diffusion.py` to label the observable as lateral `xy` diffusion and report real-DOPC-reference effective timescale fields.
  - Changed `sweep_cgl_dissipation.py` to default to the CGL-only bilayer input instead of the 1RKL protein-bilayer checkpoint.
  - Changed `stable_cgl_mass_boundary.py` so diffusion/timescale fields are measured from the CGL-only bilayer run; protein-bilayer runs remain stability checks only.
  - Ran a native-mass CGL-only bilayer measurement at `tau=5`, `time_step=0.002`, `nsteps=20000`.
  - Ran explicit 1:1 inner-step validation on the CGL-only bilayer and the active 1RKL protein-bilayer checkpoint.
  - Cleaned up the 1:1 stepping implementation by removing public `INNER_STEP`/`--inner-step` workflow knobs and using fixed MARTINI MD constants instead.
- Files modified:
  - `example/16.MARTINI/analyze_lipid_diffusion.py`
  - `example/16.MARTINI/sweep_cgl_dissipation.py`
  - `example/16.MARTINI/stable_cgl_mass_boundary.py`
  - `example/16.MARTINI/validate_cgl_mass_stability.py`
  - `example/16.MARTINI/run_sim_hybrid.sh`
  - `example/16.MARTINI/run.py`
  - `example/16.MARTINI/test_cg_bilayer/run_test.py`
  - `py/martini_prepare_system.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py example/16.MARTINI/run.py example/16.MARTINI/sweep_cgl_dissipation.py example/16.MARTINI/validate_cgl_mass_stability.py example/16.MARTINI/stable_cgl_mass_boundary.py example/16.MARTINI/test_cg_bilayer/run_test.py example/16.MARTINI/analyze_lipid_diffusion.py` passed.
  - Native CGL-only bilayer run completed: `outputs/cgl_bilayer_timescale_native_tau5/test.input.native_mass.tau5p0.up`.
  - Lateral diffusion fit: `D_xy=0.017319996 A^2/native_time`, `0.008659998 um^2/s` under the Upside `40 ps/step` mapping, fit `R^2=0.97935`.
  - Effective dry-MARTINI step from real-DOPC references: `0.041238 ps` for `8.4 um^2/s`, `0.028867 ps` for `12.0 um^2/s`, and `0.024743 ps` for `14.0 um^2/s`.
  - Explicit 1:1 CGL-only bilayer run completed under `outputs/cgl_bilayer_inner1_native_tau5/`; geometry passed with no leaflet crossings, output `time_final=40.0`, and measured `D_xy=0.0148313 um^2/s` under the Upside `40 ps/step` mapping.
  - Explicit 1:1 1RKL protein-bilayer run completed under `outputs/inner_step1_native_stability/`; CGL geometry passed with no leaflet crossings, protein passed with CA Rg `12.74 -> 13.02 A`, final CA RMSD `1.65 A`, and output `time_final=40.0`.
  - Help output for `martini_prepare_system.py run-hybrid-workflow`, `sweep_cgl_dissipation.py`, and `validate_cgl_mass_stability.py` no longer includes an inner-step option.
- Residual notes:
  - Do not use 1RKL/protein-bilayer diffusion to infer the dry-MARTINI lipid timescale. Use CGL-only bilayer lateral diffusion; use protein-bilayer trajectories only for stability checks.

## 2026-06-15 Stable-Gated CGL Mass Boundary Sweep

- Actions taken:
  - Added `example/16.MARTINI/stable_cgl_mass_boundary.py` to run bilayer-only validation, protein-bilayer validation, and diffusion fitting as one gated workflow.
  - Updated the stability validator so dense protein-bilayer initial contacts are judged relative to the starting geometry instead of a standalone bilayer-only absolute threshold.
  - Ran low and moderate mass grids at uniform `--time-step 0.002`.
- Files modified:
  - `example/16.MARTINI/stable_cgl_mass_boundary.py`
  - `example/16.MARTINI/validate_cgl_mass_stability.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `mass=0.02,0.03,0.05,0.1,0.2,0.5 m_up` failed standalone bilayer geometry for `tau=5,10,20,50`.
  - `mass=1.0,tau=5` passed standalone bilayer geometry and protein stability, but failed protein-bilayer CGL geometry; `mass=1.0,tau=10,20,50` failed standalone bilayer geometry.
  - No row passed the full stable gate, so diffusion fitting was correctly skipped for every tested point.
  - Stable-boundary artifacts were written under `example/16.MARTINI/outputs/cgl_stable_mass_boundary_40ps_v2/` and `example/16.MARTINI/outputs/cgl_stable_mass_boundary_40ps_extended/`.
- Residual notes:
  - Do not promote any low-mass CGL setting from this sweep. The target crossing seen at `mass=0.015 m_up` came from an unstable regime, and stable-gated candidates up to `1.0 m_up` did not survive validation.

## 2026-06-15 CGL Mass=0.015 Stability and Stokes-Einstein Check

- Actions taken:
  - Added `example/16.MARTINI/validate_cgl_mass_stability.py` for copied-input mass/tau validation.
  - Ran 20k-step bilayer-only and 1RKL protein-bilayer validations with `CGL/CGLD mass=0.015 m_up`, `tau=5.0635`, `gamma=0.19749`, and uniform `--time-step 0.002`.
  - Ran the native-mass bilayer-only test for the same 20k-step duration as a control.
  - Generated a short-run Stokes-Einstein plot/fits for the `mass=0.015` sweep points.
- Files modified:
  - `example/16.MARTINI/validate_cgl_mass_stability.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python3 -m py_compile example/16.MARTINI/validate_cgl_mass_stability.py example/16.MARTINI/sweep_cgl_dissipation.py example/16.MARTINI/analyze_lipid_diffusion.py` passed.
  - Native-mass bilayer-only baseline passed 20k steps.
  - `mass=0.015` bilayer-only validation failed: leaflet crossings `14/35`, `bad_parallel=26`, `bad_flip=45`, same-leaflet nearest-neighbor minimum `0.040 A`.
  - `mass=0.015` 1RKL protein-bilayer validation failed CGL geometry: leaflet crossings `30/4`, `bad_parallel=8`, `bad_flip=2`, same-leaflet nearest-neighbor minimum `0.120 A`.
  - Protein in the failed 1RKL validation stayed finite but drifted: CA Rg `12.71 -> 14.33 A`, CA RMSD final `5.52 A`, hbond last-20 mean `24.19`.
  - Stokes-Einstein check does not support linear mobility scaling: `D` versus `1/gamma` fit `R^2=0.507` with negative slope.
- Residual notes:
  - Do not promote `mass=0.015 m_up` as a production setting. It matched the target diffusion only in short trajectories and is unstable at the uniform production timestep.

## 2026-06-15 Smaller-Mass CGL Diffusion Sweep

- Actions taken:
  - Extended the sweep driver to scan multiple CGL/CGLD mass overrides in one run.
  - Kept the hybrid simulation timestep uniform at `--time-step 0.002` for every atom and every sweep point.
  - Ran two smaller-mass grids from the active 1RKL stage-7 checkpoint.
- Files modified:
  - `example/16.MARTINI/sweep_cgl_dissipation.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python3 -m py_compile example/16.MARTINI/sweep_cgl_dissipation.py example/16.MARTINI/analyze_lipid_diffusion.py` passed.
  - `mass=0.5,0.25,0.1,0.05,0.02 m_up` sweep completed 25/25 runs with no failures but stayed below the `0.5 um^2/s` target.
  - `mass=0.015,0.01,0.005 m_up` sweep completed 12/12 runs with no failures and crossed the target.
  - Best direct target match: `mass=0.015 m_up`, `tau=5`, `gamma=0.2`, `D=0.500349 um^2/s`.
  - Interpolated short-run recommendation for `mass=0.015 m_up`: `tau=5.0635`, `gamma=0.19749`.
- Residual notes:
  - The recommendation is based on short 5000-step calibration trajectories. Run a longer validation before treating it as a production default.

## 2026-06-15 CGL Mass=1 Diffusion Diagnostic

- Actions taken:
  - Added a copied-input `--cgl-mass-up` override to the CGL dissipation sweep driver.
  - Ran a 5000-step `CGL/CGLD mass=1.0 m_up`, `tau=20` stability check.
  - Ran the full `mass=1.0 m_up` dissipation sweep for `tau=5,10,20,50,100,200,500`.
- Files modified:
  - `example/16.MARTINI/sweep_cgl_dissipation.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python3 -m py_compile example/16.MARTINI/analyze_lipid_diffusion.py example/16.MARTINI/sweep_cgl_dissipation.py` passed.
  - Stability check completed with finite positions; `CGL` and `CGLD` masses were both patched to `1.0 m_up`.
  - Full mass=1 sweep completed all seven tau points without numerical failure.
  - Best measured mass=1 point: `tau=10`, `gamma=0.1`, `D=0.101894 um^2/s`; target remains `0.5 um^2/s`.
- Residual notes:
  - Mass=1 improves CGL diffusion by about `6.5x` over native mass but remains about `4.9x` too slow for the current target.

## 2026-06-15 Cross-System CGL Cutoff and Speed Check

- Actions taken:
  - Compared active 1RKL and 1AFO coarse HDF5 cutoff metadata for CGL-CGL, CGL-target, and SC-CGL runtime nodes.
  - Checked the 1RKL and 1AFO wrapper scripts for workflow-specific CGL cutoff overrides.
  - Compared active stage-7 coarse/full timings for both systems.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - 1RKL and 1AFO both use CGL-CGL `cutoff_ang=41.3`, CGL-target `cutoff_ang=26.6`, SC-CGL `cutoff_ang=26.6`, bead nonbonded cutoff `1.2 nm`, DOPC axis radius `14.55848 A`, and DOPC perpendicular radius `4.114344 A`.
  - 1RKL timing: CGL `1707.99 us/systems/step`; full lipid `3328.29 us/systems/step`.
  - 1AFO timing: CGL `1923.47 us/systems/step`; full lipid `2522.99 us/systems/step`.
  - Final CGL orientation checks have no bad-parallel or flipped rods for either active system.
- Residual notes:
  - No code change was required because the active outputs already use the same cutoff metadata and satisfy the CGL-faster-than-full requirement for both systems.

## 2026-06-15 Active 1RKL CGL Orientation and Timing Check

- Actions taken:
  - Audited the reported active `1rkl.stage_7.0.vtf` issue from the HDF5 CGL/CGLD vectors rather than the display file alone.
  - Confirmed the old active stage-7 trajectory had real CGL orientation outliers, already present at the stage-7 production handoff.
  - Regenerated a default 1RKL coarse hybrid trajectory with the current directional CGL pairlist and strict stage handoff.
  - Replaced `example/16.MARTINI/outputs/martini_1rkl_hybrid` with the validated regenerated output.
  - Preserved the previous active output as `example/16.MARTINI/outputs/martini_1rkl_hybrid.pre_directional_pairlist_backup_20260615_140424`.
  - Compared current CGL and full-resolution timings for 1RKL hybrid and same-size DOPC-only bilayer tests.
- Files and outputs modified:
  - `example/16.MARTINI/outputs/martini_1rkl_hybrid`
  - `example/16.MARTINI/outputs/martini_1rkl_hybrid.pre_directional_pairlist_backup_20260615_140424`
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - Active replacement `checkpoints/1rkl.stage_7.0.up` final CGL metrics: aligned-z min/p05/mean `0.770/0.849/0.949`, `bad_parallel=0`, `bad_flip=0`, leaflet crossings `0/0`.
  - Active replacement `logs/stage_7.0.log`: coarse CGL hybrid `1707.99 us/systems/step`.
  - Existing active full-resolution comparator `logs/stage_7.0.log`: `3328.29 us/systems/step`.
  - Same-size DOPC-only bilayer timings measured from command output: CGL `350.29 us/systems/step`; full-resolution lipid `797.32 us/systems/step`.
- Residual notes:
  - No additional performance code change was made in this phase because the current measured CGL paths are already faster than the comparable full-resolution paths.

## 2026-06-15 CGL Directional Pairlist Support

- Actions taken:
  - Added DOPC-derived `max_axis_radius_ang` metadata.
  - Propagated CGL body geometry and dry-MARTINI bead cutoff attrs into runtime CGL-CGL and CGL-target nodes.
  - Replaced broad spherical CGL-CGL and CGL-target candidate filters with conservative capsule-envelope filters using the dry-MARTINI 1.2 nm bead cutoff plus pairlist buffer.
  - Kept SC-CGL pairlists on their existing spherical support until side-chain body extents are explicit enough for a conservative directional filter.
  - Rebuilt installed MARTINI parameter files, including `parameters/dryMARTINI/dopc.h5`.
- Files modified:
  - `src/martini_cg_lipid.cpp`
  - `py/martini_cg_lipid_params.py`
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system_lib.py`
  - `parameters/dryMARTINI/dopc.h5`
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - `python3 -m py_compile py/martini_cg_lipid_params.py py/martini_build_tables.py py/martini_prepare_system_lib.py` passed.
  - `make -C obj -j4` passed.
  - `example/16.MARTINI/build_martini_h5_m1.sh` completed successfully.
  - Rebuilt `dopc.h5` records `max_axis_radius_ang=14.55848`, `max_perp_radius_ang=4.114344`, and `bead_nonbonded_cutoff_nm=1.2` on CGL table groups.
  - `example/16.MARTINI/test_cg_bilayer/run_test.sh` passed with `validation=PASS`.
  - Fresh short `1rkl` coarse hybrid smoke under `outputs/phase_cgl_directional_1rkl_smoke` completed through stage 7.0.
  - Candidate audit on the fresh 1RKL stage: CGL-CGL candidates reduced from `14678` spherical pairs to `4579` capsule pairs; CGL-target candidates reduced from `1238` to `679`.
- Residual notes:
  - Old hybrid `.up` files without `compose_vector6d/max_axis_radius_ang` are intentionally rejected during reinjection with rebuilt CGL tables. Regenerate stages instead of mixing old geometry metadata with new tables.

## 2026-06-15 CGL Diffusion Calibration

- Actions taken:
  - Added an optional per-atom OU thermostat timescale dataset read by the engine.
  - Added a CGL-only thermostat-timescale workflow knob for coarse MARTINI hybrid stage files.
  - Added a local lateral lipid diffusion analyzer for coarse CGL centers and full DOPC COMs.
  - Ran a one-step 1RKL coarse smoke with `CG_LIPID_THERMOSTAT_TIMESCALE=8.0`.
  - Corrected the target after user feedback: full DOPC is not the target; the calibration target is the CGL diffusion implied by the Upside `40 ps/step` mapping.
  - Added and ran a CGL dissipation sweep/plot driver.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `src/main.cpp`
  - `src/thermostat.cpp`
  - `src/thermostat.h`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/run_sim_hybrid.sh`
  - `example/16.MARTINI/analyze_lipid_diffusion.py`
  - `example/16.MARTINI/sweep_cgl_dissipation.py`
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py example/16.MARTINI/analyze_lipid_diffusion.py` passed.
  - `make -C obj -j4` passed.
  - Short 1RKL coarse smoke completed through stage 7.0.
  - Smoke stage-7 HDF5 has `248` atoms at global timescale `5.0` and `550` CGL/CGLD atoms at timescale `8.0`.
  - `python3 -m py_compile example/16.MARTINI/analyze_lipid_diffusion.py example/16.MARTINI/sweep_cgl_dissipation.py` passed.
  - Sweep completed for `tau=5,10,20,50,100,200,500` from the active 1RKL coarse stage-7 file.
  - Plot/table written under `example/16.MARTINI/outputs/cgl_dissipation_sweep_40ps/`.
  - Best measured point: `tau=20`, `gamma=0.05`, `D=0.015631 um^2/s`; target is `0.5 um^2/s`.
- Residual notes:
  - No production dissipation constant is recommended from this sweep. The measured range is about `32x` too slow at best, so dissipation-only tuning did not solve the 40 ps/step CGL timescale calibration.
