# bilayer-lateral-diffusion

This folder contains a bilayer-only dry-MARTINI workflow that measures DOPC lateral diffusion across damping, mass, and temperature conditions.

## Scope

- The workflow prepares one shared DOPC bilayer system and reuses it across all scan tasks.
- Each task runs a fixed stage ladder with minimization, NPT equilibration, and one production measurement stage.
- Bilayer fluidity is measured from leaflet-resolved DOPC `PO4` lateral MSD on the production `.up` output.
- The workflow supports local execution and Slurm array staging/submission.
- The current bundled DOPC system is used as a relative screening benchmark, not an absolute diffusion calibration benchmark.
- After simulation completion, a separate stage-7-only post-analysis path can re-measure diffusion from surviving `stage_7.0.up` checkpoints even if stage-6 files were deleted.

## Default Scan

- Damping values: `4, 6, 8, 12, 16, 32, 64`
- Mass scales: `1.0, 0.75, 0.5, 0.25, 0.1`
- Temperatures: `0.70, 0.80, 0.90, 1.00, 1.10`
- Replicates: `3`
- Total tasks: `525`

## Recommended Next Scan

- Keep the thermostat consistent with the standard Upside examples: `damping = 5.0`
- Keep the scanned Upside temperatures unchanged: `0.70, 0.80, 0.90, 1.00, 1.10`
- Push only the bilayer mass lower than the previous scan:
  - `mass_scale = 0.1, 0.05, 0.02, 0.01, 0.005, 0.002`
- Replicates: `3`
- Total tasks: `90`

Launch that focused rerun with:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_tau5_mass_only_round.sh
```

The default output directory for that preset is `bilayer-lateral-diffusion/runs/tau5_mass_only`.

## Recommended Joint Mass+Damping Follow-Up

- Keep dryMARTINI at the same temperature as Upside for each condition.
- Probe a moderate weaker-damping window around the aligned Upside baseline:
  - `damping = 5, 8, 12, 16, 24`
  - `mass_scale = 0.02, 0.01, 0.005, 0.002, 0.001`
  - `temperature = 0.70, 0.80, 0.90, 1.00, 1.10`
  - `replicates = 3`
- This is `375` tasks total.
- Why this grid:
  - `tau = 5` is the current Upside-consistent baseline,
  - `8 -> 24` tests weaker damping without jumping straight to the previously noisy extreme-`tau` regime,
  - `0.02 -> 0.001` focuses on the low-mass region that already sped the bilayer up, while pushing one decade lower than the last robust hot-temperature points.

Launch it with:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_mass_damping_followup_round.sh
```

The default output directory for that preset is `bilayer-lateral-diffusion/runs/mass_damping_followup`.

## Recommended Boundary Scan

- This is the next scan to run after the refreshed joint analysis.
- Keep dryMARTINI at the same temperature as Upside for each condition.
- Focus on the transition between the robust `tau = 5-8` region and the unstable `tau = 16-24` region:
  - `damping = 5, 8, 10, 12, 14, 16`
  - `mass_scale = 0.02, 0.015, 0.01, 0.0075, 0.005, 0.003, 0.002`
  - `temperature = 0.70, 0.80, 0.90, 1.00, 1.10`
  - `replicates = 3`
- This is `630` tasks total.
- Why this grid:
  - keep `tau = 5` and `8` as the robust anchors already observed in the data,
  - add `10, 12, 14, 16` to locate the onset of the instability frontier,
  - drop `tau = 24` because it is already clearly pathological,
  - drop `mass = 0.001` because it repeatedly produced NaNs and exploded diffusion,
  - refine the useful mass window between `0.015` and `0.002` where the speed / stability tradeoff is actually happening.

Launch it with:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_mass_damping_boundary_round.sh
```

The default output directory for that preset is `bilayer-lateral-diffusion/runs/mass_damping_boundary`.

## Recommended Integrator/Time-Step Boundary Diagnostic

- This is the next diagnostic workflow if you want to test whether the dryMARTINI instability frontier is mainly a time-integration problem.
- It runs separate subdirectories for each `(integrator, dt)` pair instead of mixing those settings inside one manifest.
- Default grid:
  - `integrator = v, nvtc`
  - `dt = 0.01, 0.0075, 0.005`
  - `damping = 5, 12, 14, 16`
  - `mass_scale = 0.005, 0.003, 0.002`
  - `temperature = 0.70, 0.90, 1.10`
  - `replicates = 3`
- This is `108` tasks per `(integrator, dt)` subrun and `648` tasks total across the full matrix.
- Why this workflow:
  - `tau = 5` is the current robust anchor,
  - `12, 14, 16` sit on the onset region of the instability frontier,
  - `0.005 -> 0.002` is the low-mass window where the blow-up behavior is actually happening,
  - comparing `v` against `nvtc` while shrinking `dt` tells you whether the frontier is mainly a timestep-resolution issue or something deeper.
- Default directory layout:
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/v_dt0p01`
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/v_dt0p0075`
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/v_dt0p005`
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/nvtc_dt0p01`
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/nvtc_dt0p0075`
  - `bilayer-lateral-diffusion/runs/integrator_timestep_boundary/nvtc_dt0p005`

Launch the full diagnostic matrix with:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_integrator_timestep_boundary_round.sh
```

Optional environment variables for this diagnostic:

- `BILAYER_DIFF_PROBE_BASE_ROOT`
- `BILAYER_DIFF_PROBE_INTEGRATORS`
- `BILAYER_DIFF_PROBE_TIME_STEPS`
- `BILAYER_DIFF_PROBE_NVT_MAX_FORCE`: default `0.0`
- `BILAYER_DIFF_DAMPING_VALUES`
- `BILAYER_DIFF_MASS_SCALES`
- `BILAYER_DIFF_TEMPERATURE_VALUES`
- `BILAYER_DIFF_REPLICATES`

By default `nvtc` is tested with `max-force = 0.0`, so this compares integrator behavior without turning on force clipping. If you specifically want to test force clipping, override `BILAYER_DIFF_PROBE_NVT_MAX_FORCE`.

## Measurement

- Main metric: DOPC `PO4` 2D lateral diffusion coefficient
- Analysis details:
  - identify lipid molecules from `PO4`-containing molecules in the generated `.up` file,
  - assign leaflets from the early-stage bilayer geometry with a wrapped-`z` split and keep that assignment fixed during analysis,
  - unwrap XY displacements with periodic minimum image,
  - remove whole-bilayer lateral COM drift,
  - discard the first `20%` of production frames,
  - fit `MSD = 4Dt + b` over the `10% -> 40%` lag window of the retained trajectory.

## Stage-7 Post-Analysis

- The post-analysis flow discovers surviving `tasks/*/checkpoints/bilayer_diffusion.stage_7.0.up` files.
- If a matching `stage_6.0.up` file still exists it is used as the leaflet reference; otherwise the `stage_7.0.up` file is used as its own reference.
- Post-analysis outputs are written under `stage7-analysis/` so they do not overwrite the original run-time results.

Submit the stage-7 analysis on Slurm:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_stage7_analysis.sh
```

Optional environment variables:

- `BILAYER_DIFF_BASE_DIR`
- `BILAYER_DIFF_RUN_PRESET`: default `default`; for preset runs this resolves to `bilayer-lateral-diffusion/runs/<preset>`
- `BILAYER_DIFF_NO_SUBMIT=1`
- `BILAYER_DIFF_OVERWRITE=1`
- `BILAYER_DIFF_PYTHON`
- `BILAYER_DIFF_HDF5_MODULE`
- `BILAYER_DIFF_ANALYSIS_WALLTIME`
- `BILAYER_DIFF_ANALYSIS_COLLECT_WALLTIME`
- `BILAYER_DIFF_ANALYSIS_CPUS_PER_TASK`

For the focused `tau=5` mass-only rerun, use:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_tau5_mass_only_stage7_analysis.sh
```

For the joint mass+damping follow-up rerun, use:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_mass_damping_followup_stage7_analysis.sh
```

For the boundary-focused rerun, use:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_mass_damping_boundary_stage7_analysis.sh
```

For the integrator/time-step diagnostic matrix, use:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_integrator_timestep_boundary_stage7_analysis.sh
```

This wrapper walks the `integrator_timestep_boundary` subruns and stages analysis only for those that already contain `stage_7.0.up` outputs.

Pack the finished results for download:

```bash
source .venv/bin/activate
./bilayer-lateral-diffusion/pack_stage7_results.sh
```

Optional environment variables:

- `BILAYER_DIFF_BASE_DIR`
- `BILAYER_DIFF_PACK_OUTPUT`
- `BILAYER_DIFF_PACK_METADATA_ONLY=1`
- `BILAYER_DIFF_PYTHON`

## Offline Report For Downloaded Analysis

If you have already downloaded a completed `stage7-analysis/` tree, generate a local report with plots and a calibration summary:

```bash
source .venv/bin/activate
source source.sh
python3 bilayer-lateral-diffusion/report_stage7_analysis.py \
  --analysis-dir bilayer-lateral-diffusion/stage7-analysis
```

This writes outputs under `bilayer-lateral-diffusion/stage7-analysis/report/`, including:

- validation summary JSON,
- converted CSV tables in physical units,
- heatmaps and main-effect plots,
- a markdown recommendation summary tied to the chosen Upside timestep calibration.

Generate a temperature-specific bilayer correction that keeps the scanned Upside temperature fixed and uses the repo unit contract plus the DOPC experimental target curve:

```bash
source .venv/bin/activate
source source.sh
python3 bilayer-lateral-diffusion/temperature_specific_calibration.py \
  --analysis-dir bilayer-lateral-diffusion/stage7-analysis
```

This adds:

- `stage7-analysis/report/temperature_specific_corrections.csv`
- `stage7-analysis/report/temperature_specific_correction.md`
- `stage7-analysis/report/temperature_specific_correction.png`

## Local Usage

From the project root:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/run_local.sh
```

Optional environment variables:

- `BILAYER_DIFF_BASE_DIR`: run directory, default `bilayer-lateral-diffusion/runs/default`
- `BILAYER_DIFF_FORCE_INIT=1`: regenerate the manifest and shared system even if they already exist
- `BILAYER_DIFF_BILAYER_PDB`
- `BILAYER_DIFF_DAMPING_VALUES`
- `BILAYER_DIFF_MASS_SCALES`
- `BILAYER_DIFF_TEMPERATURE_VALUES`
- `BILAYER_DIFF_REPLICATES`
- `BILAYER_DIFF_SEED`
- `BILAYER_DIFF_XY_SCALE`
- `BILAYER_DIFF_BOX_PADDING_Z`
- `BILAYER_DIFF_SALT_MOLAR`
- `BILAYER_DIFF_ION_CUTOFF`
- `BILAYER_DIFF_MIN_60_MAX_ITER`
- `BILAYER_DIFF_MIN_61_MAX_ITER`
- `BILAYER_DIFF_EQUILIBRATION_STEPS`
- `BILAYER_DIFF_PRODUCTION_STEPS`
- `BILAYER_DIFF_TIME_STEP`
- `BILAYER_DIFF_INTEGRATOR`
- `BILAYER_DIFF_MAX_FORCE`
- `BILAYER_DIFF_EQUILIBRATION_FRAME_STEPS`
- `BILAYER_DIFF_PRODUCTION_FRAME_STEPS`
- `BILAYER_DIFF_MAX_TASKS`: limit local execution to the first `N` tasks
- `BILAYER_DIFF_START_TASK`: start local execution at task index `N`
- `BILAYER_DIFF_OVERWRITE=1`: rerun completed tasks
- `BILAYER_DIFF_NO_ASSEMBLE=1`: skip aggregation after local execution
- `BILAYER_DIFF_PYTHON`: Python interpreter override
- `BILAYER_DIFF_HDF5_MODULE`: HDF5 module name, default `hdf5/1.14.3`

## Slurm Usage

From the project root on a Slurm system:

```bash
source .venv/bin/activate
source source.sh
./bilayer-lateral-diffusion/submit_remote_round.sh
```

Optional environment variables:

- `BILAYER_DIFF_BASE_DIR`
- `BILAYER_DIFF_FORCE_INIT=1`
- `BILAYER_DIFF_BILAYER_PDB`
- `BILAYER_DIFF_DAMPING_VALUES`
- `BILAYER_DIFF_MASS_SCALES`
- `BILAYER_DIFF_TEMPERATURE_VALUES`
- `BILAYER_DIFF_REPLICATES`
- `BILAYER_DIFF_SEED`
- `BILAYER_DIFF_XY_SCALE`
- `BILAYER_DIFF_BOX_PADDING_Z`
- `BILAYER_DIFF_SALT_MOLAR`
- `BILAYER_DIFF_ION_CUTOFF`
- `BILAYER_DIFF_MIN_60_MAX_ITER`
- `BILAYER_DIFF_MIN_61_MAX_ITER`
- `BILAYER_DIFF_EQUILIBRATION_STEPS`
- `BILAYER_DIFF_PRODUCTION_STEPS`
- `BILAYER_DIFF_TIME_STEP`
- `BILAYER_DIFF_INTEGRATOR`
- `BILAYER_DIFF_MAX_FORCE`
- `BILAYER_DIFF_EQUILIBRATION_FRAME_STEPS`
- `BILAYER_DIFF_PRODUCTION_FRAME_STEPS`
- `BILAYER_DIFF_NO_SUBMIT=1`
- `BILAYER_DIFF_PYTHON`
- `BILAYER_DIFF_HDF5_MODULE`
- `BILAYER_DIFF_SBATCH_PARTITION`
- `BILAYER_DIFF_SBATCH_ACCOUNT`
- `BILAYER_DIFF_SBATCH_QOS`
- `BILAYER_DIFF_SBATCH_CONSTRAINT`
- `BILAYER_DIFF_SBATCH_MEM`
- `BILAYER_DIFF_TRAIN_WALLTIME`
- `BILAYER_DIFF_COLLECT_WALLTIME`

To stage scripts without submitting:

```bash
python3 bilayer-lateral-diffusion/workflow.py submit-slurm --base-dir bilayer-lateral-diffusion/runs/default --no-submit
```

To rerun the focused `tau=5` mass-only scan under a different output directory or with a tweaked mass grid:

```bash
BILAYER_DIFF_BASE_DIR=bilayer-lateral-diffusion/runs/tau5_mass_only_retry \
BILAYER_DIFF_MASS_SCALES=0.1,0.05,0.02,0.01,0.005,0.002 \
./bilayer-lateral-diffusion/submit_tau5_mass_only_round.sh
```

## Main Commands

```bash
python3 bilayer-lateral-diffusion/workflow.py init-run --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py run-local --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py run-array-task --round-manifest bilayer-lateral-diffusion/runs/default/slurm/round_manifest.json --task-id 0
python3 bilayer-lateral-diffusion/workflow.py assemble-results --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py submit-slurm --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py init-stage7-analysis --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py run-stage7-analysis-local --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py run-stage7-analysis-task --round-manifest bilayer-lateral-diffusion/runs/default/stage7-analysis/slurm/round_manifest.json --task-id 0
python3 bilayer-lateral-diffusion/workflow.py assemble-stage7-analysis --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py submit-stage7-analysis-slurm --base-dir bilayer-lateral-diffusion/runs/default
python3 bilayer-lateral-diffusion/workflow.py pack-stage7-results --base-dir bilayer-lateral-diffusion/runs/default
```

## Outputs

Each run directory contains:

- `diffusion_manifest.json`: task manifest and run settings
- `system/`: shared prepared bilayer assets
- `tasks/<task-code>/`: per-task stage files and logs
- `results/tasks/*.json`: one JSON result per task
- `assembled/task_results.csv`: completed task table
- `assembled/condition_summary.csv`: replicate-aggregated condition summary
- `slurm/`: generated Slurm array and collector scripts
- `stage7-analysis/`: post hoc stage-7 analysis manifest, task JSON files, assembled CSV/JSON outputs, and Slurm scripts
- `download/*.tar.gz`: packed download archives created by `pack-stage7-results`
