# hybrid-interface-sweep

This folder stages and submits sweeps over `PROTEIN_ENV_INTERFACE_SCALE` using the real hybrid `example/16.MARTINI/run_sim_1rkl.sh` workflow.

## Scope

- Each task runs one full hybrid workflow instance for one `(interface_scale, replicate)` pair.
- Every task gets its own task-local `RUN_DIR`, logs, checkpoints, and result JSON.
- The workflow supports local execution and Slurm array staging/submission.
- The workflow records execution status only in v1; it does not add a new scientific analysis pipeline.

## Default Sweep

- `interface_scale = 1.00, 0.85, 0.70, 0.55, 0.40, 0.25`
- `replicates = 3`
- `pdb_id = 1rkl`

The default output directory is `hybrid-interface-sweep/runs/default`.

## Run Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
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

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
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

## Captured Hybrid Overrides

At `init-run` time, the workflow captures a whitelist of relevant `run_sim_1rkl.sh` environment overrides into the manifest so they are replayed consistently on local runs and Slurm array tasks.

Examples:

- `TEMPERATURE`
- `THERMOSTAT_TIMESCALE`
- `MIN_60_MAX_ITER`
- `EQ_62_NSTEPS` through `EQ_66_NSTEPS`
- `PROD_70_NSTEPS`
- `EQ_FRAME_STEPS`
- `PROD_FRAME_STEPS`
- `PROD_70_NPT_ENABLE`
- `PROD_70_BACKBONE_FIX_RIGID_ENABLE`
- `PRODUCTION_NONPROTEIN_HARD_SPHERE`

To capture additional environment keys, set:

- `HYBRID_SWEEP_EXTRA_ENV_KEYS=KEY1,KEY2,...`
