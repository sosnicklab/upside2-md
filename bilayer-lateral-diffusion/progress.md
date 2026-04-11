# Progress Log

## 2026-04-06 (Initialization)
- Created a dedicated root workflow folder for bilayer lateral diffusion screening.
- Recorded the implementation goal, design choices, and execution phases in the local tracker files before code changes.

## 2026-04-06 (Implementation)
- Added `bilayer-lateral-diffusion/workflow.py` with:
  - run initialization and manifest generation,
  - per-task stage preparation and execution,
  - production-stage MSD analysis,
  - result aggregation,
  - Slurm array and collector staging.
- Added `bilayer-lateral-diffusion/run_local.sh`, `bilayer-lateral-diffusion/submit_remote_round.sh`, and `bilayer-lateral-diffusion/README.md`.
- Marked the shell wrappers executable so the documented `./...sh` entrypoints work directly.

## 2026-04-06 (Verification And Fixes)
- Ran `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`.
- Ran `bash -n bilayer-lateral-diffusion/run_local.sh`.
- Ran `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`.
- Initialized a reduced smoke run in `/tmp/bilayer_diff_smoke`.
- Verified `submit-slurm --no-submit` writes Slurm manifests and scripts for the smoke run.
- Ran one smoke array task end-to-end and fixed two analysis issues discovered during verification:
  - lipid selection cannot rely on `particle_class == LIPID` for the current bilayer prep path,
  - leaflet assignment should use an early-stage reference geometry instead of the production-stage first frame.
- Re-ran the smoke task after both fixes and verified it writes a successful task JSON and aggregated CSV/JSON outputs.

## 2026-04-06 (Wrapper Portability Fix)
- Investigated a `submit_remote_round.sh` failure reported from the Slurm environment during `init-run`.
- Identified two portability issues:
  - the shell wrappers and generated Slurm scripts could run `workflow.py` under the module-loaded `python3` instead of the project venv,
  - workflow child processes inherited ambient `UPSIDE_HOME` and related repo-root variables instead of forcing the current repository root.
- Updated `run_local.sh`, `submit_remote_round.sh`, and generated Slurm script content in `workflow.py` so `PYTHON_BIN` is resolved after activation and prefers `${VIRTUAL_ENV}/bin/python3`.
- Updated `workflow.py` runtime environment setup so `UPSIDE_HOME`, MARTINI/NPT variables, `PATH`, and `PYTHONPATH` are set deterministically for child commands.
- Updated child-command failure reporting to include the log path and the tail of the failing log.
- Re-verified:
  - `bash -n bilayer-lateral-diffusion/run_local.sh`
  - `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`
  - `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`
  - `BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_smoke BILAYER_DIFF_NO_SUBMIT=1 ./bilayer-lateral-diffusion/submit_remote_round.sh`
  - `BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_smoke BILAYER_DIFF_MAX_TASKS=0 ./bilayer-lateral-diffusion/run_local.sh`
  - `source .venv/bin/activate && BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_smoke BILAYER_DIFF_NO_SUBMIT=1 ./bilayer-lateral-diffusion/submit_remote_round.sh`
  - `source .venv/bin/activate && BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_smoke BILAYER_DIFF_MAX_TASKS=0 ./bilayer-lateral-diffusion/run_local.sh`

## 2026-04-06 (Cluster Python Compatibility Fix)
- Investigated a new cluster traceback from `martini_prepare_system_lib.py` during `init-run`.
- Root cause: `py/martini_prepare_system_lib.py` used PEP 604 annotations like `Path | None` without `from __future__ import annotations`, which breaks import-time annotation evaluation under Python 3.9.
- Added `from __future__ import annotations` at the top of `py/martini_prepare_system_lib.py`.
- Re-verified:
  - `python3 -m py_compile py/martini_prepare_system_lib.py py/martini_prepare_system.py bilayer-lateral-diffusion/workflow.py`
  - normal import with `PYTHONPATH=py python3 -c 'import martini_prepare_system_lib'`

## 2026-04-06 (Checkout Root Preservation Fix)
- Investigated the cluster command path still jumping from the invoked checkout under `/home/...` to `/project/...`.
- Root cause: `bilayer-lateral-diffusion/workflow.py` used `Path(__file__).resolve()`, which canonicalizes symlinks / mounted paths and can redirect the workflow into a different checkout tree.
- Replaced that with `Path(os.path.abspath(__file__))` so the workflow stays rooted in the checkout path that was actually launched.
- Re-verified with `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`.

## 2026-04-06 (Cluster HDF5 Runtime Fix)
- Investigated a Slurm array-task runtime failure from `obj/upside`:
  - `error while loading shared libraries: libhdf5.so.310`
- User provided the cluster-specific correction that `hdf5/1.14.3` must be loaded with the other modules.
- Updated:
  - `bilayer-lateral-diffusion/submit_remote_round.sh`
  - `bilayer-lateral-diffusion/run_local.sh`
  - generated Slurm train/collector scripts in `bilayer-lateral-diffusion/workflow.py`
- Added `BILAYER_DIFF_HDF5_MODULE` with default `hdf5/1.14.3`.
- Re-verified:
  - `bash -n bilayer-lateral-diffusion/run_local.sh`
  - `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`
  - `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`
  - regenerated Slurm scripts contain `module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true`

## 2026-04-08 (Stage-7 Post-Analysis Workflow)
- Added new post hoc commands to `bilayer-lateral-diffusion/workflow.py`:
  - `init-stage7-analysis`
  - `run-stage7-analysis-local`
  - `run-stage7-analysis-task`
  - `assemble-stage7-analysis`
  - `submit-stage7-analysis-slurm`
  - `pack-stage7-results`
- Added wrapper scripts:
  - `bilayer-lateral-diffusion/submit_stage7_analysis.sh`
  - `bilayer-lateral-diffusion/pack_stage7_results.sh`
- Updated `bilayer-lateral-diffusion/README.md` with stage-7 post-analysis and packing usage.
- Verified:
  - `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`
  - `bash -n bilayer-lateral-diffusion/submit_stage7_analysis.sh`
  - `bash -n bilayer-lateral-diffusion/pack_stage7_results.sh`
  - `python3 bilayer-lateral-diffusion/workflow.py init-stage7-analysis --base-dir bilayer-lateral-diffusion/runs/default`
  - `python3 bilayer-lateral-diffusion/workflow.py run-stage7-analysis-local --base-dir bilayer-lateral-diffusion/runs/default --max-tasks 1 --overwrite`
  - `python3 bilayer-lateral-diffusion/workflow.py submit-stage7-analysis-slurm --base-dir bilayer-lateral-diffusion/runs/default --no-submit`
  - `python3 bilayer-lateral-diffusion/workflow.py pack-stage7-results --base-dir bilayer-lateral-diffusion/runs/default --metadata-only --output /tmp/bilayer-stage7-metadata.tar.gz`
  - wrapper entrypoints:
    - `BILAYER_DIFF_BASE_DIR=$PWD/bilayer-lateral-diffusion/runs/default BILAYER_DIFF_NO_SUBMIT=1 ./bilayer-lateral-diffusion/submit_stage7_analysis.sh`
    - `BILAYER_DIFF_BASE_DIR=$PWD/bilayer-lateral-diffusion/runs/default BILAYER_DIFF_PACK_METADATA_ONLY=1 BILAYER_DIFF_PACK_OUTPUT=/tmp/bilayer-stage7-wrapper.tar.gz ./bilayer-lateral-diffusion/pack_stage7_results.sh`
- Explicitly verified the stage-7-only fallback path by calling `_analyze_production(stage7_file, reference_stage_file=None)` on one completed `stage_7.0.up`; it used the stage-7 file as its own leaflet reference and completed successfully.
- The current local completed run used for verification contains `65` surviving `stage_7.0.up` files; the stage-7 analysis manifest records that discovered set.

## 2026-04-08 (Downloaded Stage-7 Report And Calibration Check)
- Validated the downloaded `bilayer-lateral-diffusion/stage7-analysis/` tree for internal consistency:
  - manifest stage-7 count = `525`,
  - summary total / successful tasks = `525 / 525`,
  - `task_results.csv` rows = `525`,
  - `condition_summary.csv` rows = `175`,
  - task result JSON count = `525`.
- Measured fit-quality health on the downloaded analysis:
  - median `fit_r2_mean = 0.9904`,
  - `2` tasks below `R^2 < 0.9`,
  - `10` tasks below `R^2 < 0.95`.
- Added `bilayer-lateral-diffusion/report_stage7_analysis.py` to:
  - validate a downloaded `stage7-analysis/` tree,
  - convert diffusion into `nm^2/ns` and `um^2/s` using a user-specified Upside timestep calibration,
  - generate report plots and recommendation tables under `stage7-analysis/report/`.
- Updated `bilayer-lateral-diffusion/README.md` with offline report usage.
- Verified:
  - `python3 -m py_compile bilayer-lateral-diffusion/report_stage7_analysis.py`
  - `python3 bilayer-lateral-diffusion/report_stage7_analysis.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- Generated report outputs under `bilayer-lateral-diffusion/stage7-analysis/report/`, including:
  - `diffusion_heatmaps_um2_s.png`
  - `main_effects_um2_s.png`
  - `interaction_trends_um2_s.png`
  - `fit_quality.png`
  - `report.md`
  - `validation_summary.json`
  - `recommendation_summary.json`

## 2026-04-08 (Temperature-Specific Fixed-Temperature Calibration)
- User corrected the calibration goal: keep the scanned Upside temperature unchanged and use the repository’s real-temperature conversion from `AGENTS.md` instead of changing temperature to compensate for bilayer slowness.
- Added `bilayer-lateral-diffusion/temperature_specific_calibration.py`.
- The new script:
  - reads `stage7-analysis/report/condition_summary_physical_units.csv`,
  - converts `T_up -> K` with `1.0 T_up = 350.588235 K`,
  - fits an Arrhenius target curve from the DOPC diffusion values reported by Bag et al. 2014,
  - derives per-temperature mass corrections at stable damping values.
- Verified:
  - `python3 -m py_compile bilayer-lateral-diffusion/temperature_specific_calibration.py`
  - `python3 bilayer-lateral-diffusion/temperature_specific_calibration.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- Generated new outputs:
  - `bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_corrections.csv`
  - `bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.md`
  - `bilayer-lateral-diffusion/stage7-analysis/report/temperature_specific_correction.png`
- Result:
  - damping remains a weak correction knob,
  - the current scan minimum mass `0.1` is still too heavy at every temperature,
  - target-matching mass scales are approximately `0.068`, `0.042`, `0.032`, `0.030`, and `0.030` for `T_up = 0.7, 0.8, 0.9, 1.0, 1.1`.

## 2026-04-08 (Saved Hybrid Settings Summary)
- Added root-level `hybrid_timescale.md`.
- Saved the fixed-temperature bilayer correction table there so the recommended `(temperature, damping, mass scale)` settings are easy to reuse without reopening the generated report directory.

## 2026-04-08 (Thermostat Einstein-Relation Audit)
- Investigated the user-reported concern that the damping control should satisfy the Einstein relation and that the observed dependence was not behaving linearly.
- Inspected:
  - `src/thermostat.h`
  - `src/thermostat.cpp`
  - `src/main.cpp`
  - `src/deriv_engine.cpp`
  - `src/martini.cpp`
- Found a real bug:
  - the integrator already uses `/input/mass` in the position update,
  - but the OU thermostat still injected the same momentum-noise amplitude for all atoms, independent of mass.
- Patched the thermostat to accept the engine pointer and scale noise by `sqrt(mass)` when MARTINI masses are loaded.
- Patched kinetic-energy logging so it uses `p^2 / m` when MARTINI masses are active.
- Rebuilt successfully with `cmake --build obj --target upside`.
- Marked `hybrid_timescale.md` as provisional because the old bilayer diffusion scan was run with the buggy thermostat; the scan should be repeated before trusting any previous mass/damping recommendations.

## 2026-04-09 (Post-Fix Reanalysis)
- Revalidated the refreshed downloaded `bilayer-lateral-diffusion/stage7-analysis/` tree after the thermostat fix:
  - `525` stage-7 files,
  - `525` successful tasks,
  - `175` condition rows,
  - no tasks below `R^2 < 0.95`.
- Re-ran:
  - `python3 bilayer-lateral-diffusion/report_stage7_analysis.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
  - `python3 bilayer-lateral-diffusion/temperature_specific_calibration.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- Found that the refreshed post-fix scan differs qualitatively from the provisional pre-fix one:
  - mass effect is now weak across the tested range,
  - temperature is the larger lever inside the current grid,
  - `tau=64` remains a noisy branch and should not anchor calibration.
- Updated `report_stage7_analysis.py` so its recommendation section is data-driven for the post-fix regime instead of repeating the old pre-fix interpretation.
- Hardened `temperature_specific_calibration.py` so it reports `unresolved` when the current grid is too weakly mass-dependent to justify a literal extrapolated mass correction.
- Rewrote root-level `hybrid_timescale.md` with the post-fix conclusion:
  - no reliable mass correction is resolved by the current grid,
  - the best currently tested stable settings are `mass = 0.1` with `tau = 12` for `T_up = 0.7-0.9` and `tau = 16` for `T_up = 1.0-1.1`.

## 2026-04-09 (Focused Tau-5 Mass-Only Rerun Surface)
- Audited the standard example workflows and confirmed the repo-wide example convention is to inherit the engine thermostat default `tau = 5.0`, while the prior bilayer scan and workflow `16` had been using explicit non-default damping values.
- Concluded that the next bilayer rerun should:
  - keep the Upside-consistent damping fixed at `tau = 5.0`,
  - keep the scanned temperatures unchanged,
  - push only the mass grid lower because the refreshed post-fix scan always hit its fastest point at `mass = 0.1`.
- Updated wrapper scripts:
  - `bilayer-lateral-diffusion/submit_remote_round.sh`
  - `bilayer-lateral-diffusion/run_local.sh`
  so they now pass optional init-run settings from environment variables when creating a manifest.
- Added `bilayer-lateral-diffusion/submit_tau5_mass_only_round.sh` as a convenience preset for the next rerun.
- Updated `bilayer-lateral-diffusion/README.md` with:
  - the recommended next scan,
  - the new wrapper environment variables,
  - one-command usage for the focused `tau=5` mass-only rerun.
- Verification:
  - `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`
  - `bash -n bilayer-lateral-diffusion/run_local.sh`
  - `bash -n bilayer-lateral-diffusion/submit_tau5_mass_only_round.sh`
  - `python3 -m py_compile bilayer-lateral-diffusion/workflow.py bilayer-lateral-diffusion/report_stage7_analysis.py bilayer-lateral-diffusion/temperature_specific_calibration.py`
  - staged a no-submit run with:
    - `BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_tau5_mass_only`
    - `BILAYER_DIFF_NO_SUBMIT=1`
    - `BILAYER_DIFF_FORCE_INIT=1`
    - `./bilayer-lateral-diffusion/submit_tau5_mass_only_round.sh`
  - inspected the generated manifest and confirmed:
    - `damping_values = [5.0]`
    - `mass_scales = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002]`
  - `temperature_values = [0.7, 0.8, 0.9, 1.0, 1.1]`
  - `replicates = 3`
  - `task_count = 90`

## 2026-04-10 (Stage-7 Wrapper Alignment For Tau-5 Preset)
- Updated `bilayer-lateral-diffusion/submit_stage7_analysis.sh` so its default base directory can follow `BILAYER_DIFF_RUN_PRESET`, instead of always assuming `runs/default`.
- Added `bilayer-lateral-diffusion/submit_tau5_mass_only_stage7_analysis.sh` as the matching stage-7 analysis entrypoint for the focused `tau=5` mass-only rerun.
- Updated `bilayer-lateral-diffusion/README.md` with:
  - the new `BILAYER_DIFF_RUN_PRESET` environment variable for stage-7 analysis,
  - the one-command `submit_tau5_mass_only_stage7_analysis.sh` usage.
- Verification:
  - `bash -n bilayer-lateral-diffusion/submit_stage7_analysis.sh`
  - `bash -n bilayer-lateral-diffusion/submit_tau5_mass_only_stage7_analysis.sh`
  - created a fake surviving checkpoint under `/tmp/bilayer_diff_tau5_mass_only/tasks/tau5_m0p1_t0p7_r01/checkpoints/bilayer_diffusion.stage_7.0.up`
  - `source .venv/bin/activate && source source.sh && BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_tau5_mass_only BILAYER_DIFF_NO_SUBMIT=1 ./bilayer-lateral-diffusion/submit_tau5_mass_only_stage7_analysis.sh`
  - confirmed it wrote:
    - `stage7-analysis/analysis_manifest.json`
    - `stage7-analysis/slurm/round_manifest.json`
    - `stage7-analysis/slurm/analyze_stage7_array.sbatch`
    - `stage7-analysis/slurm/collect_stage7_analysis.sbatch`

## 2026-04-10 (Focused Tau-5 Mass-Only Reanalysis)
- Revalidated the refreshed downloaded `bilayer-lateral-diffusion/stage7-analysis/` tree for the focused mass-only rerun:
  - `90` stage-7 files,
  - `90` task rows,
  - `30` condition rows,
  - no tasks below `R^2 < 0.95`.
- Re-ran:
  - `python3 bilayer-lateral-diffusion/report_stage7_analysis.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
  - `python3 bilayer-lateral-diffusion/temperature_specific_calibration.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- Found and fixed two analysis issues exposed by the focused rerun:
  - `temperature_specific_calibration.py` assumed the older stable damping list `4,6,8,12,16`; it now falls back to the damping values actually present in the current report, which lets the fixed-`tau=5` scan analyze cleanly.
  - `report_stage7_analysis.py` contained stale broad-scan wording (`tau = 4 to 64`, `mass = 0.1 -> 1.0`) and now writes scan-aware markdown for fixed-damping reruns.
- Updated `hybrid_timescale.md` to reflect the focused `tau=5` mass-only rerun:
  - temperature effect `~3.20x`,
  - mass effect `~1.61x`,
  - best robust point `T = 1.1`, `tau = 5`, `mass = 0.005`, `D = 0.365 um^2/s`,
  - lowest tested mass `0.002` is already noisy at `T = 0.9` and `1.1`,
  - all temperatures remain well below the physical DOPC target.
- Verification:
  - `python3 -m py_compile bilayer-lateral-diffusion/report_stage7_analysis.py`
  - `python3 -m py_compile bilayer-lateral-diffusion/temperature_specific_calibration.py`
  - report regeneration completed with validation status `valid`
  - temperature-specific calibration regenerated `temperature_specific_corrections.csv`, `temperature_specific_correction.md`, and `temperature_specific_correction.png`

## 2026-04-10 (Joint Mass+Damping Follow-Up Preset)

## 2026-04-10 (Refreshed Joint-Scan Reanalysis)
- Revalidated the refreshed downloaded `bilayer-lateral-diffusion/stage7-analysis/` tree for the joint mass+damping follow-up scan:
  - `375` stage-7 files,
  - `375` successful task rows,
  - `125` condition rows,
  - `1` task below `R^2 < 0.95`,
  - `62` task rows with `NaN` diffusion outputs in the raw assembled table.
- Re-ran:
  - `python3 bilayer-lateral-diffusion/report_stage7_analysis.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
  - `python3 bilayer-lateral-diffusion/temperature_specific_calibration.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- Found that the first robust-subset hardening pass was still too permissive:
  - a low-mass `tau = 24` branch with internally smooth but obviously nonphysical diffusion still survived the old adjacent-mass filter,
  - the temperature-specific calibration was still selecting high-CV conditions because it only enforced `R^2`.
- Updated `bilayer-lateral-diffusion/report_stage7_analysis.py`:
  - changed the continuity filter to treat the first large adjacent-mass jump as the start of an unstable lighter-mass branch,
  - tightened the default jump threshold from `10x` to `4x`,
  - added `stable_condition_summary_physical_units.csv`,
  - added raw-grid versus robust-subset max diffusion to `recommendation_summary.json`,
  - rewrote the markdown recommendation text for the joint scan.
- Updated `bilayer-lateral-diffusion/temperature_specific_calibration.py`:
  - default stable damping window is now `5,8,12,16`,
  - added a robust-subset filter using completed replicates, `R^2 >= 0.95`, `CV <= 0.25`, and `continuity_stable` when present,
  - rewrote the markdown so it explicitly states that target-like diffusion appears only on excluded branches.
- Re-verified:
  - `python3 -m py_compile bilayer-lateral-diffusion/report_stage7_analysis.py`
  - `python3 -m py_compile bilayer-lateral-diffusion/temperature_specific_calibration.py`
  - regenerated `report.md`, `recommendation_summary.json`, `stable_condition_summary_physical_units.csv`, `temperature_specific_correction.md`, `temperature_specific_corrections.csv`, and the updated plots.
- Final refreshed conclusion from the joint scan:
  - only `49 / 125` conditions remain in the robust subset,
  - robust-subset temperature effect is about `2.95x`,
  - robust-subset damping effect is about `1.86x`,
  - robust-subset mass effect is about `1.32x`,
  - best robust point is `T = 1.1`, `tau = 5`, `mass = 0.01`, `D = 0.435 um^2/s`,
  - no robust condition reaches `1.0 um^2/s`,
  - target-like diffusion appears only on excluded `tau = 16-24` low-mass branches.

## 2026-04-10 (Boundary-Focused Next Scan Preset)
- Added `bilayer-lateral-diffusion/submit_mass_damping_boundary_round.sh` as the next focused rerun entrypoint.
- Added `bilayer-lateral-diffusion/submit_mass_damping_boundary_stage7_analysis.sh` as the matching post-run stage-7 analysis entrypoint.
- Chosen default grid:
  - `tau = 5, 8, 10, 12, 14, 16`
  - `mass = 0.02, 0.015, 0.01, 0.0075, 0.005, 0.003, 0.002`
  - `temperature = 0.70, 0.80, 0.90, 1.00, 1.10`
  - `replicates = 3`
- Updated `bilayer-lateral-diffusion/README.md` with:
  - the new recommended boundary scan section,
  - launch instructions for the new preset,
  - launch instructions for the new stage-7 analysis wrapper.
- Verification:
  - `bash -n bilayer-lateral-diffusion/submit_mass_damping_boundary_round.sh`
  - `bash -n bilayer-lateral-diffusion/submit_mass_damping_boundary_stage7_analysis.sh`
  - `source .venv/bin/activate && source source.sh && BILAYER_DIFF_BASE_DIR=/tmp/bilayer_diff_mass_damping_boundary BILAYER_DIFF_NO_SUBMIT=1 BILAYER_DIFF_FORCE_INIT=1 ./bilayer-lateral-diffusion/submit_mass_damping_boundary_round.sh`
  - confirmed `Task count: 630`
  - staged `submit_mass_damping_boundary_stage7_analysis.sh` successfully after placing one fake `stage_7.0.up` file under the temp preset tree, which confirmed the preset pathing and Slurm-script generation.

## 2026-04-11 (Boundary-Scan Reanalysis)
- Revalidated the refreshed downloaded `bilayer-lateral-diffusion/stage7-analysis/` tree for the boundary-focused rerun:
  - `630` stage-7 files,
  - `630` successful task rows,
  - `210` condition rows,
  - `16` task rows with `NaN` diffusion outputs,
  - `9` condition rows with non-finite diffusion,
  - `4` tasks below `R^2 < 0.95`, with minimum `R^2 = 0.799`.
- Re-ran:
  - `python3 bilayer-lateral-diffusion/report_stage7_analysis.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
  - `python3 bilayer-lateral-diffusion/temperature_specific_calibration.py --analysis-dir bilayer-lateral-diffusion/stage7-analysis`
- The report remained `invalid` at the full-grid level because of the raw NaNs and weak-fit outliers, but the existing robust-subset filter still behaved correctly on the new grid; no additional analysis-script patch was needed.
- Regenerated outputs:
  - `report.md`
  - `recommendation_summary.json`
  - `stable_condition_summary_physical_units.csv`
  - `temperature_specific_correction.md`
  - `temperature_specific_corrections.csv`
  - updated plots under `stage7-analysis/report/`
- Main refreshed conclusion from the boundary scan:
  - robust subset now contains `125 / 210` conditions,
  - robust-subset temperature effect is about `3.02x`,
  - robust-subset damping effect is about `1.33x`,
  - robust-subset mass effect is about `1.30x`,
  - best robust point is `T = 1.1`, `tau = 5`, `mass = 0.005`, `D = 0.394 um^2/s`,
  - no robust condition reaches `1.0 um^2/s`,
  - `tau = 10-16` does not beat the `tau = 5` branch in the trustworthy subset and mostly maps the approach to the instability frontier.
- Updated `hybrid_timescale.md` so the root summary now reflects the boundary scan instead of the earlier broad joint scan.
- Chose the next joint follow-up grid based on the focused `tau=5` rerun plus the earlier broad damping scan:
  - keep temperature fixed to the same Upside temperature,
  - keep `tau = 5` as the aligned baseline,
  - probe a moderate weaker-damping window `5, 8, 12, 16, 24`,
  - focus mass on `0.02, 0.01, 0.005, 0.002, 0.001`.
- Added new preset wrappers:
  - `bilayer-lateral-diffusion/submit_mass_damping_followup_round.sh`
  - `bilayer-lateral-diffusion/submit_mass_damping_followup_stage7_analysis.sh`
- Updated `bilayer-lateral-diffusion/README.md` with:
  - the rationale for the joint scan,
  - the `375`-task count,
  - the new run and stage-7 analysis entrypoints.

## 2026-04-11 (Integrator/Time-Step Diagnostic Preset)
- Updated `bilayer-lateral-diffusion/workflow.py` so `init-run` accepts `--integrator` and `--max-force`, stores them in the manifest, and uses them for MD stages while keeping minimization on `integrator = v`.
- Updated `bilayer-lateral-diffusion/submit_remote_round.sh` and `bilayer-lateral-diffusion/run_local.sh` to forward `BILAYER_DIFF_INTEGRATOR` and `BILAYER_DIFF_MAX_FORCE`.
- Added `bilayer-lateral-diffusion/submit_integrator_timestep_boundary_round.sh`.
- Added `bilayer-lateral-diffusion/submit_integrator_timestep_boundary_stage7_analysis.sh`.
- Updated `bilayer-lateral-diffusion/README.md` with the diagnostic matrix:
  - `integrator = v, nvtc`
  - `dt = 0.01, 0.0075, 0.005`
  - `damping = 5, 12, 14, 16`
  - `mass_scale = 0.005, 0.003, 0.002`
  - `temperature = 0.70, 0.90, 1.10`
  - `replicates = 3`
  - `648` tasks across six subruns
- Verification:
  - `python3 -m py_compile bilayer-lateral-diffusion/workflow.py`
  - `bash -n bilayer-lateral-diffusion/submit_integrator_timestep_boundary_round.sh`
  - `bash -n bilayer-lateral-diffusion/submit_integrator_timestep_boundary_stage7_analysis.sh`
  - `source .venv/bin/activate && source source.sh && BILAYER_DIFF_PROBE_BASE_ROOT=/tmp/bilayer_diff_integrator_timestep_boundary BILAYER_DIFF_PROBE_INTEGRATORS=nvtc BILAYER_DIFF_PROBE_TIME_STEPS=0.0075 BILAYER_DIFF_NO_SUBMIT=1 BILAYER_DIFF_FORCE_INIT=1 ./bilayer-lateral-diffusion/submit_integrator_timestep_boundary_round.sh`
  - confirmed the staged subrun manifest reports `integrator = nvtc`, `dt = 0.0075`, `max_force = 0.0`, and `108` tasks
  - placed one fake `tasks/.../checkpoints/bilayer_diffusion.stage_7.0.up` file under the temp subrun tree
  - `source .venv/bin/activate && source source.sh && BILAYER_DIFF_PROBE_BASE_ROOT=/tmp/bilayer_diff_integrator_timestep_boundary BILAYER_DIFF_PROBE_INTEGRATORS=nvtc BILAYER_DIFF_PROBE_TIME_STEPS=0.0075 BILAYER_DIFF_NO_SUBMIT=1 ./bilayer-lateral-diffusion/submit_integrator_timestep_boundary_stage7_analysis.sh`
  - confirmed the stage-7 wrapper wrote `analysis_manifest.json`, `round_manifest.json`, `analyze_stage7_array.sbatch`, and `collect_stage7_analysis.sbatch`
