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

Useful overrides:

```bash
export PROFILE=dimer3
export BASE_DIR=/Users/yinhan/Documents/upside2-md/ConDiv_symlay/test_dimer3
export INIT_FORCEFIELD_DIR=/Users/yinhan/Documents/upside2-md/parameters/ff_2.1
export WORKER_LAUNCH=auto
export CONDIV_SYMLAY_DENSE_GRID_SIZE=801
export CONDIV_SYMLAY_SUPPORT_MARGIN=0.5
```

## Local Restart

```bash
./run_local.sh 20
```

## Slurm Restart

```bash
sbatch run_remote.sh 20
```

`run_remote.sh` uses the same restart surface as `run_local.sh`, requests `16G` memory by default, and caps BLAS/OpenMP threads to keep worker launches predictable.

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
