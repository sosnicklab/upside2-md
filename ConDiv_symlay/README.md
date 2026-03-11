# ConDiv_symlay Membrane Workflow

`ConDiv_symlay/` is a symmetric multilayer clone of `ConDiv/`.

It keeps the same checkpoint/restart workflow surface:
- `ConDiv_mem.py`
- `run_init.sh`
- `run_local.sh`
- `run_remote.sh`
- `run_validate_rounds.sh`
- `check_membrane_gradient.py`

and adds the symmetric-layer tooling:
- `layer_template.json`
- `build_layer_manifest.py`
- `validate_symlay_constraints.py`
- `symlay_utils.py`

The training target still uses the standard `membrane.h5` schema. The difference is that initialization and every optimizer update are projected into a symmetric topology-slot subspace derived from the DOPC bilayer structure.

This workflow now also uses one fixed DOPC membrane thickness for every training target. It does not read per-target `.thickness` files during ConDiv training.

The design rationale for the training potential and constraints is documented in [TRAINING_POTENTIAL_DESIGN.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/TRAINING_POTENTIAL_DESIGN.md).

## Environment

Run from the project root:

```bash
cd /Users/yinhan/Documents/upside2-md
source .venv/bin/activate
source source.sh
```

The runner scripts do this automatically.

## Topology-Slot Model

The hard constrained full-sequence type template is:

```text
Q0-Qa-Na-C1-C3-C1-C1-C1-C1-C3-C1-Na-Qa-Q0
```

The slot manifest is built from:
- `example/16.MARTINI/pdb/bilayer.MARTINI.pdb`
- `example/16.MARTINI/ff_dry/dry_martini_v2.1_lipids.itp`

and written into the active run directory as:
- `layer_manifest.json`
- `layer_manifest.csv`
- `layer_manifest.png`

## Initialization

Default initialization starts from the current published force field:

```bash
cd /Users/yinhan/Documents/upside2-md/ConDiv_symlay
./run_init.sh
```

This does four things:
1. builds the DOPC topology-slot manifest
2. seeds a symmetric constrained `membrane.h5` under `<base_dir>/seed_forcefield/`
3. creates `initial_checkpoint.pkl`
4. validates the seeded checkpoint with `validate_symlay_constraints.py`

By default this creates the run under `ConDiv_symlay/test_dimer3`. If that directory already contains an initialized training run, `./run_init.sh` now refuses to overwrite it and tells you to resume with `./submit_remote_round.sh` instead.

After a successful init, `run_init.sh` records the resolved run directory in `ConDiv_symlay/.condiv_current_run_dir`. The zero-argument `./submit_remote_round.sh` and `sbatch run_remote.sh` paths both use that recorded run directory by default, so stale inherited `BASE_DIR` values do not silently redirect training to a different checkout. If that record is missing, the scripts fall back to discovering initialized run directories under the current `ConDiv_symlay` checkout before they use the built-in `ConDiv_symlay/test_dimer3` default.

Useful overrides:

```bash
export PROFILE=dimer3
export INIT_FORCEFIELD_DIR=/Users/yinhan/Documents/upside2-md/parameters/ff_2.1
export WORKER_LAUNCH=auto
export CONDIV_MINIBATCH_SIZE=15
export CONDIV_DOPC_THICKNESS=30.2
export CONDIV_SYMLAY_DENSE_GRID_SIZE=801
export CONDIV_SYMLAY_SUPPORT_MARGIN=0.5
```

`CONDIV_DOPC_THICKNESS` is the global membrane thickness in Angstrom used for every target in the run. The default is `30.2`, matching the median of the current `upside_input2/*.thickness` values while removing target-specific variation.

`CONDIV_MINIBATCH_SIZE` is fixed when `initial_checkpoint.pkl` is created. If you want larger Slurm fanout, set it before `./run_init.sh` and reinitialize the run directory. For the current default `pdb_list2`, there are 45 proteins total, so the practical upper bound is 45.

## Local Restart

```bash
./run_local.sh 20
```

## Separate-Job Slurm Restart

```bash
./submit_remote_round.sh
```

This is the recommended remote workflow when you want each protein simulation to receive its own Slurm job ID. The submitter stages the next minibatch round under:

- `<base_dir>/epoch_<epoch>_minibatch_<minibatch>/slurm/worker_specs/*.json`
- `<base_dir>/epoch_<epoch>_minibatch_<minibatch>/slurm/simulate_array.sbatch`
- `<base_dir>/epoch_<epoch>_minibatch_<minibatch>/slurm/round_manifest.json`

Then it:

- submits one Slurm job array for the whole minibatch round
- maps one protein simulation to one array task
- submits `run_remote_update.sh` with an `afterany` dependency on that array job
- lets the update job compute the FF update, write `checkpoint.pkl`, append `training_progress.jsonl`, update `training_status.json`, and submit the next round when convergence has not been reached

The simulation jobs use:

- `1` Slurm array task per protein
- `1` Slurm task inside each array element
- `CONDIV_OMP_THREADS` CPUs per task
- `CONDIV_N_REPLICA` local Upside replicas inside that job

The generated array script carries its own `#SBATCH --array=0-(n-1)` line, so there is only one simulation `sbatch` file per round.

By default `submit_remote_round.sh` uses `RUN_STEPS=0`, which means no submitter-imposed chain limit. If you pass one positional integer, it is treated as the remaining number of rounds to auto-submit in that chain:

```bash
./submit_remote_round.sh 10
```

Useful Slurm overrides for the separate-job workflow:

```bash
export CONDIV_SIM_WALLTIME=24:00:00
export CONDIV_UPDATE_WALLTIME=02:00:00
export CONDIV_SIM_MEM=16G
export CONDIV_UPDATE_MEM=8G
export CONDIV_SBATCH_PARTITION=caslake
export CONDIV_SBATCH_ACCOUNT=<account>
```

`run_remote_update.sh` is the single Slurm-side update driver. It is normally submitted automatically by `submit_remote_round.sh`; manual `sbatch` is only needed for debugging or recovery, and requires `CONDIV_ROUND_MANIFEST` plus `CONDIV_PROJECT_ROOT` in the environment.

## Legacy Single-Allocation Slurm Restart

```bash
sbatch run_remote.sh
```

`run_remote.sh` keeps the older model where one Slurm allocation owns the entire minibatch and fan-out happens inside that allocation via `srun`. That mode is still available, but it does not give each protein simulation a separate Slurm job ID.

`run_remote.sh` uses the same restart surface as `run_local.sh`, requests the full node memory by default with `#SBATCH --mem=0`, resolves the real workflow directory from `SLURM_SUBMIT_DIR` / `CONDIV_PROJECT_ROOT` so Slurm spool-copy execution still finds the repo checkout, prefers `ConDiv_symlay/venv/bin/activate` when present (falling back to the project-root `.venv` in this checkout), sources `source.sh` from the project root, loads `python/3.11.9`, `cmake`, and `openmpi` through the module system, and caps BLAS/OpenMP threads to keep worker launches predictable.

With no argument, `run_remote.sh` defaults to `RUN_STEPS=20`. You can still pass a single positional step count if you want to override that default.

Inside Slurm, `WORKER_LAUNCH=auto` now resolves to `srun`. The launch model matches the reference `/Users/yinhan/Documents/Train` workflow:
- one `srun --ntasks=1 --cpus-per-task=<CONDIV_OMP_THREADS>` step per protein worker
- each worker launches one local `upside` replica bundle
- the optimizer updates once after the whole minibatch finishes

The wrapper defaults:
- `CONDIV_N_REPLICA=8`
- `CONDIV_OMP_THREADS=8`
- `CONDIV_MAX_PARALLEL_WORKERS=6`

The Slurm wrapper is now hard-coded to a fixed 48-slot / 48-CPU layout:
- `#SBATCH --ntasks-per-node=48`
- `#SBATCH --cpus-per-task=1`
- interpreted as `48` CPU slots total
- `6` protein workers in parallel
- `8` CPUs per worker
- `8` replicas per worker
- worker steps use `srun --cpu-bind=cores`

It now fails fast if the live Slurm allocation does not match that `48 slots -> 6 x 8 = 48 CPUs` layout. This is intentional; the wrapper no longer tries to infer or adapt the CPU layout from the allocation.

After each Slurm job, the wrapper:
- appends one record to `<base_dir>/training_progress.jsonl`
- writes the latest convergence summary to `<base_dir>/training_status.json`
- resubmits itself with `afterok` until the convergence window passes

Default convergence settings:

```bash
export CONDIV_CONVERGENCE_GRAD_NORM=1.0
export CONDIV_CONVERGENCE_UPDATE_NORM=0.05
export CONDIV_CONVERGENCE_PATIENCE=3
export CONDIV_AUTO_RESUBMIT=1
export CONDIV_AUTO_MAX_SUBMISSIONS=0
```

`CONDIV_AUTO_MAX_SUBMISSIONS=0` means no wrapper-imposed submission limit. Convergence is checked on the trailing `CONDIV_CONVERGENCE_PATIENCE` progress records, requiring every record in that window to satisfy both the total gradient-norm and total update-norm thresholds.

## Validation

Run gradient + constraint validation rounds:

```bash
./run_validate_rounds.sh
```

Useful overrides:

```bash
export ROUNDS=3
export STEPS_PER_ROUND=2
export FD_SAMPLES=24
export FD_EPS=1e-3
```

Per-round reports:
- `<base_dir>/gradient_round_<round>.json`
- `<base_dir>/constraint_round_<round>.json`

## Notes

- The workflow copies `ConDiv_mem.py` and `symlay_utils.py` into the run directory so checkpoints remain restartable.
- `ConDiv_symlay/` intentionally does not carry the nested `ConDiv/venv/`; it uses the project-root `.venv`.
- Output membrane files remain standard `membrane.h5` files, so existing readers continue to work.
