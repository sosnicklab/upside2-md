This example runs hybrid Upside/dry-MARTINI membrane systems with
full-resolution 14-bead DOPC lipids.

Setup:

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the consolidated dry-MARTINI force-field file when the force-field inputs
change (writes `parameters/ff_2.1/martini.h5`):

```bash
bash build_martini_h5_m1.sh
```

Run a full-resolution 1RKL system:

```bash
bash run_sim_1rkl_full.sh
```

The shell wrappers call `run_sim_hybrid.sh`. Common overrides are environment
variables:

```bash
RUN_DIR=outputs/test_1rkl \
RUNTIME_PDB_ID=test_1rkl \
EQ_60_NSTEPS=500 \
PROD_70_NSTEPS=10000 \
bash run_sim_1rkl_full.sh
```

The Python wrapper exposes the same workflow through editable settings near the
top of `run.py`:

```bash
python run.py
```

Outputs are written under `outputs/`. Each run contains stage checkpoints, logs,
and VTF files.
