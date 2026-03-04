# ConDiv Membrane Workflow (Modernized)

This folder contains a Python 3 membrane ConDiv workflow ported from `/Users/yinhan/Documents/Train`, adapted for both:
- local M1 execution
- Slurm execution

It also includes strict gradient validation (finite-difference + sanity checks) for multi-round training.

## Contents

- `ConDiv_mem.py`: modernized membrane training script (`initialize`, `restart`, `worker` modes)
- `check_membrane_gradient.py`: strict membrane gradient validator
- `run_init.sh`: initialize a run directory
- `run_local.sh`: resume training locally
- `run_remote.sh`: Slurm resume script
- `run_validate_rounds.sh`: run several rounds + gradient validation each round
- `setup_venv.sh`: optional standalone venv setup for this folder
- copied membrane inputs from `/Users/yinhan/Documents/Train`:
  - `param0/ param1/ param2/ param3/`
  - `upside_input/ upside_input2/`
  - `pdb_list pdb_list2`
  - helper scripts (`get_dmemb.py`, `get_potential_3b.py`, `get_potential_3d.py`, `show_memb*.py`)

## Environment

Preferred environment is project root:

```bash
cd /Users/yinhan/Documents/upside2-md
source .venv/bin/activate
source source.sh
```

All runner scripts do this automatically.

## Profiles

- `dimer3` (default): `upside_input2 + pdb_list2`
- `legacy60`: `upside_input + pdb_list`

Set with:

```bash
export PROFILE=dimer3
# or
export PROFILE=legacy60
```

## Initialization

Default initialization uses the current force field:

`/Users/yinhan/Documents/upside2-md/parameters/ff_2.1`

Run:

```bash
cd /Users/yinhan/Documents/upside2-md/ConDiv
./run_init.sh
```

Useful overrides:

```bash
export INIT_FORCEFIELD_DIR=/Users/yinhan/Documents/upside2-md/parameters/ff_2.1
export BASE_DIR=/Users/yinhan/Documents/upside2-md/ConDiv/test_dimer3
export WORKER_LAUNCH=auto   # auto|local|srun
```

## Local Restart

```bash
./run_local.sh 20
```

`20` is restart iterations. If omitted, `RUN_STEPS` or default `20` is used.

## Slurm Restart

```bash
sbatch run_remote.sh 20
```

## Multi-round Training + Gradient Validation

Run several rounds and validate gradients after each round:

```bash
./run_validate_rounds.sh
```

Useful overrides:

```bash
export ROUNDS=3
export STEPS_PER_ROUND=2
export FD_SAMPLES=24
export FD_EPS=1e-3
export FD_REL_MEDIAN_THRESHOLD=5e-2
export FD_REL_MAX_THRESHOLD=2.5e-1
```

Per-round reports are written to:

- `${BASE_DIR}/gradient_round_<round>.json`

The script exits non-zero if any round fails gradient validation.

## Notes

- `WORKER_LAUNCH=auto` defaults to local worker subprocesses for cross-platform parity.
- `WORKER_LAUNCH=srun` enables per-worker `srun` launch in cluster environments.
- Full production-scale convergence is not part of smoke validation; use reduced rounds/steps for portability checks first.
