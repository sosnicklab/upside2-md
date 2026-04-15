# hybrid-interface-sweep

This folder now runs an interface-scale RMSF calibration for `1rkl`.

It compares two methods:

- Reference: `example/08.MembraneSimulation/3.fixed_curvature.run.py` physics for `1rkl`
- Swept method: `example/16.MARTINI/run_sim_1rkl.sh`

The sweep variable is `PROTEIN_ENV_INTERFACE_SCALE`. The analysis computes residue-wise backbone RMSF for both methods, fits a trend line over RMSF error versus `interface_scale`, and writes plot-ready CSV files plus a recommended scale.

## Scope

- Each `reference` task runs one independent example-08-style fixed-curvature replicate for `1rkl`.
- Each `hybrid` task runs one example-16 hybrid replicate for one `interface_scale`.
- The workflow supports local execution and Slurm array staging/submission.
- The workflow also supports a post-run RMSF analysis phase over completed outputs.

## Default Sweep

- `interface_scale = 1.00, 0.85, 0.70, 0.55, 0.40, 0.25`
- `reference_replicates = 4`
- `hybrid_replicates = 3`
- `pdb_id = 1rkl`
- RMSF burn-in fraction: `0.20`
- trend-line sample count: `201`

The default output directory is `hybrid-interface-sweep/runs/default`.

Old softening or diffusion manifests are intentionally incompatible with this workflow. Reinitialize with `HYBRID_SWEEP_FORCE_INIT=1` or use a fresh `HYBRID_SWEEP_BASE_DIR`.

## Run Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_HYBRID_REPLICATES`
- `HYBRID_SWEEP_REFERENCE_REPLICATES`
- `HYBRID_SWEEP_REPLICATES`
  - fallback that applies the same count to both task families
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_BURN_IN_FRACTION`
- `HYBRID_SWEEP_TRENDLINE_SAMPLES`
- `HYBRID_SWEEP_FORCE_INIT=1`
- `HYBRID_SWEEP_MAX_TASKS`
- `HYBRID_SWEEP_START_TASK`
- `HYBRID_SWEEP_OVERWRITE=1`

## Submit To Slurm

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/submit_remote_round.sh
```

This stages and submits a Slurm array with one array element per task from `sweep_manifest.json`, followed by a dependent collector job.

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_HYBRID_REPLICATES`
- `HYBRID_SWEEP_REFERENCE_REPLICATES`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_BURN_IN_FRACTION`
- `HYBRID_SWEEP_TRENDLINE_SAMPLES`
- `HYBRID_SWEEP_FORCE_INIT=1`
- `HYBRID_SWEEP_NO_SUBMIT=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

Optional Slurm resource environment variables:

- `HYBRID_SWEEP_TRAIN_WALLTIME`
- `HYBRID_SWEEP_COLLECT_WALLTIME`
- `HYBRID_SWEEP_CPUS_PER_TASK`
- `HYBRID_SWEEP_SBATCH_PARTITION`
- `HYBRID_SWEEP_SBATCH_ACCOUNT`
- `HYBRID_SWEEP_SBATCH_QOS`
- `HYBRID_SWEEP_SBATCH_CONSTRAINT`
- `HYBRID_SWEEP_SBATCH_MEM`

## Runtime Overrides

### Reference Method

The example-08 reference side is configured from these environment variables at `init-run` time:

- `REFERENCE_TEMPERATURE`
- `REFERENCE_DURATION`
- `REFERENCE_FRAME_INTERVAL`
- `REFERENCE_MEMBRANE_THICKNESS`
- `REFERENCE_USE_CURVATURE`
- `REFERENCE_CURVATURE_RADIUS`
- `REFERENCE_CURVATURE_SIGN`

Defaults match `example/08.MembraneSimulation/3.fixed_curvature.run.py` for `1rkl`.

### Hybrid Method

The hybrid side replays a whitelist of `example/16.MARTINI/run_sim_1rkl.sh` overrides captured at `init-run` time.

Examples:

- `TEMPERATURE`
- `THERMOSTAT_TIMESCALE`
- `MIN_60_MAX_ITER`
- `EQ_62_NSTEPS` through `EQ_66_NSTEPS`
- `PROD_70_NSTEPS`
- `EQ_FRAME_STEPS`
- `PROD_FRAME_STEPS`
- `PROD_70_NPT_ENABLE`
- `PROD_70_BAROSTAT_TYPE`
- `PROTEIN_LIPID_CUTOFF`
- `SALT_MOLAR`
- `SC_ENV_LJ_FORCE_CAP`
- `SC_ENV_COUL_FORCE_CAP`

To capture additional hybrid environment keys, set:

- `HYBRID_SWEEP_EXTRA_ENV_KEYS=KEY1,KEY2,...`

To capture additional reference-side environment keys, set:

- `HYBRID_SWEEP_EXTRA_REFERENCE_KEYS=KEY1,KEY2,...`

## Analysis Outputs

The post-run analysis compares residue-wise backbone RMSF from the hybrid runs against the aggregated example-08 reference profile.

The analysis uses:

- rigid 3D alignment over the mapped backbone carrier atoms
- residue RMSF from the mapped CA-like backbone atom
- reference profile mean/std across reference replicates
- condition-level RMSF error versus `interface_scale`
- a fitted trend line over `condition_profile_rmse_angstrom`

The analysis writes under `BASE_DIR/analysis/`:

- `analysis_manifest.json`
- `results/tasks/*.json`
- `assembled/task_results.csv`
- `assembled/residue_rmsf_profiles.csv`
- `assembled/reference_profile.csv`
- `assembled/condition_profiles.csv`
- `assembled/condition_summary.csv`
- `assembled/trendline_points.csv`
- `assembled/failed_tasks.csv`
- `assembled/recommendation_summary.json`
- `assembled/summary.json`

The main plot-ready files are:

- observed condition points: `assembled/condition_summary.csv`
- fitted trend line: `assembled/trendline_points.csv`

The recommended scale is saved in:

- `assembled/recommendation_summary.json`

## Run Analysis Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_analysis_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_FORCE_ANALYSIS_INIT=1`
- `HYBRID_SWEEP_ANALYSIS_MAX_TASKS`
- `HYBRID_SWEEP_ANALYSIS_START_TASK`
- `HYBRID_SWEEP_ANALYSIS_NO_ASSEMBLE=1`
- `HYBRID_SWEEP_OVERWRITE=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

## Submit Analysis To Slurm

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/submit_analysis.sh
```

This follows the same pattern as the run phase:

- `init-analysis` discovers successful reference and hybrid outputs,
- `submit-analysis-slurm` stages an analysis Slurm array,
- a dependent collector job assembles the RMSF tables and trend-line outputs.

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_FORCE_ANALYSIS_INIT=1`
- `HYBRID_SWEEP_NO_SUBMIT=1`
- `HYBRID_SWEEP_OVERWRITE=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

Optional analysis Slurm resource environment variables:

- `HYBRID_SWEEP_ANALYSIS_WALLTIME`
- `HYBRID_SWEEP_ANALYSIS_COLLECT_WALLTIME`
- `HYBRID_SWEEP_ANALYSIS_CPUS_PER_TASK`
- `HYBRID_SWEEP_SBATCH_PARTITION`
- `HYBRID_SWEEP_SBATCH_ACCOUNT`
- `HYBRID_SWEEP_SBATCH_QOS`
- `HYBRID_SWEEP_SBATCH_CONSTRAINT`
- `HYBRID_SWEEP_SBATCH_MEM`
