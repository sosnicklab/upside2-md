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

By default this creates the run under `ConDiv_symlay/test_dimer3`. If that directory already contains an initialized training run, `./run_init.sh` now refuses to overwrite it and tells you to resume with `sbatch run_remote.sh` instead.

After a successful init, `run_init.sh` records the resolved run directory in `ConDiv_symlay/.condiv_current_run_dir`. The zero-argument `sbatch run_remote.sh` path uses that recorded run directory by default, so stale inherited `BASE_DIR` values do not silently redirect training to a different checkout. If that record is missing, `run_local.sh` and `run_remote.sh` fall back to discovering initialized run directories under the current `ConDiv_symlay` checkout before they use the built-in `ConDiv_symlay/test_dimer3` default.

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

## Slurm Restart

```bash
sbatch run_remote.sh
```

`run_remote.sh` uses the same restart surface as `run_local.sh`, requests the full node memory by default with `#SBATCH --mem=0`, resolves the real workflow directory from `SLURM_SUBMIT_DIR` / `CONDIV_PROJECT_ROOT` so Slurm spool-copy execution still finds the repo checkout, prefers `ConDiv_symlay/venv/bin/activate` when present (falling back to the project-root `.venv` in this checkout), sources `source.sh` from the project root, loads `python/3.11.9`, `cmake`, and `openmpi` through the module system, and caps BLAS/OpenMP threads to keep worker launches predictable.

With no argument, `run_remote.sh` defaults to `RUN_STEPS=20`. You can still pass a single positional step count if you want to override that default.

Inside Slurm, `WORKER_LAUNCH=auto` now resolves to `srun`. The launch model matches the reference `/Users/yinhan/Documents/Train` workflow:
- one `srun --ntasks=1 --cpus-per-task=<CONDIV_OMP_THREADS>` step per protein worker
- each worker launches one local `upside` replica bundle
- the optimizer updates once after the whole minibatch finishes

The wrapper defaults:
- `CONDIV_N_REPLICA=8`
- `CONDIV_OMP_THREADS=CONDIV_N_REPLICA`
- `CONDIV_MAX_PARALLEL_WORKERS=0`

`CONDIV_MAX_PARALLEL_WORKERS=0` means the Python layer does not cap the minibatch fanout. It launches the whole minibatch and lets Slurm determine how many protein workers can run concurrently inside the allocation.

With the default `#SBATCH --ntasks-per-node=48`, `CONDIV_N_REPLICA=8`, and `CONDIV_OMP_THREADS=8`, the allocation can sustain about `6` protein workers at once. If you want more proteins running concurrently, reduce `CONDIV_OMP_THREADS` / `CONDIV_N_REPLICA` or request more CPUs.

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
