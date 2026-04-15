# hybrid-interface-sweep

This folder stages and submits bilayer-only dry-MARTINI softening sweeps for later hybrid interface calibration.

## Scope

- Each task runs one full bilayer-only `6.0 -> 7.0` membrane workflow instance for one `(lj_alpha, slater_alpha, replicate)` condition.
- Production-stage `7.0` softening is applied only inside this sweep workflow by rewriting the staged `martini_potential` attrs:
  - `lj_soften`
  - `lj_soften_alpha`
  - `coulomb_soften`
  - `slater_alpha`
- The canonical `example/16.MARTINI/` workflow is not modified by this sweep.
- The workflow measures bilayer `PO4` lateral diffusion and reports a reciprocal-diffusion viscosity proxy.
- The workflow supports local execution and Slurm array staging/submission.
- The workflow also supports a post-run analysis phase over completed `stage_7.0.up` files.

## Default Sweep

- `lj_alpha = 0.00, 0.025, 0.05, 0.10, 0.20`
- `slater_alpha = 0.00, 0.25, 0.50, 1.00, 2.00`
- `replicates = 3`
- `pdb_id = bilayer`
- physical conversion assumption: `40 ps` per integrator step

The default output directory is `hybrid-interface-sweep/runs/default`.

If a base directory still contains a scalar `interaction_scale` manifest, reinitialize it with `HYBRID_SWEEP_FORCE_INIT=1` or use a fresh `HYBRID_SWEEP_BASE_DIR`.

## Run Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_LJ_ALPHAS`
- `HYBRID_SWEEP_SLATER_ALPHAS`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_INTEGRATION_PS_PER_STEP`
- `HYBRID_SWEEP_TARGET_DIFFUSION_UM2_S`
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

This stages and submits a Slurm array with one array element per sweep task from `sweep_manifest.json`, followed by a dependent collector job.

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_LJ_ALPHAS`
- `HYBRID_SWEEP_SLATER_ALPHAS`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_INTEGRATION_PS_PER_STEP`
- `HYBRID_SWEEP_TARGET_DIFFUSION_UM2_S`
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

## Runtime Settings

At `init-run` time the workflow resolves and stores the bilayer runtime settings so every task uses the same preparation and simulation surface.

Examples:

- `TEMPERATURE`
- `THERMOSTAT_TIMESCALE`
- `MIN_60_MAX_ITER`
- `EQ_62_NSTEPS` through `EQ_66_NSTEPS`
- `PROD_70_NSTEPS`
- `MIN_TIME_STEP`
- `EQ_TIME_STEP`
- `PROD_TIME_STEP`
- `EQ_FRAME_STEPS`
- `PROD_FRAME_STEPS`
- `PROD_70_NPT_ENABLE`
- `PROD_70_BAROSTAT_TYPE`
- `BILAYER_PDB`
- `SALT_MOLAR`
- `ION_CUTOFF`
- `BOX_PADDING_XY`
- `BOX_PADDING_Z`
- `UPSIDE_EWALD_ENABLE`
- `UPSIDE_EWALD_ALPHA`
- `UPSIDE_EWALD_KMAX`

To capture additional environment keys, set:

- `HYBRID_SWEEP_EXTRA_ENV_KEYS=KEY1,KEY2,...`

## Analysis Outputs

The post-run analysis treats bilayer `PO4` lateral diffusion as the main fluidity signal and groups results by `(lj_alpha, slater_alpha)`.

- Main signal: `PO4` lateral diffusion relative to bilayer COM.
- Physical-unit conversion:
  - uses the task production `dt`,
  - assumes `40 ps` per integrator step unless overridden at `init-run`.
- Viscosity proxy:
  - `1 / diffusion_um2_per_s`
  - reported explicitly as a proxy, not a direct viscosity measurement.
- Optional recommendation:
  - if `HYBRID_SWEEP_TARGET_DIFFUSION_UM2_S` is provided at `init-run`, the assembled analysis ranks conditions by absolute error to that target.

The analysis writes under `BASE_DIR/analysis/`:

- `analysis_manifest.json`
- `results/tasks/*.json`
- `assembled/task_results.csv`
- `assembled/condition_summary.csv`
- `assembled/failed_tasks.csv`
- `assembled/recommendation_summary.json`
- `assembled/summary.json`

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

This follows the same pattern as the sweep run path:

- `init-analysis` discovers completed `stage_7.0.up` files from successful sweep results,
- `submit-analysis-slurm` stages an analysis Slurm array with one array element per discovered checkpoint,
- a dependent collector job assembles the per-task analysis JSON files into CSV and recommendation outputs.

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
