This example runs hybrid Upside/dry-MARTINI membrane systems.

Setup:

```bash
source ../../.venv/bin/activate
source ../../source.sh
```

Build the dry-MARTINI HDF5 tables when the force-field inputs change:

```bash
bash build_martini_h5_m1.sh
```

Run a coarse-lipid 1RKL system:

```bash
bash run_sim_1rkl.sh
```

Run the same system with full-resolution DOPC:

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
bash run_sim_1rkl.sh
```

The Python wrapper exposes the same workflow through editable settings near the
top of `run.py`:

```bash
python run.py
```

Outputs are written under `outputs/`. Each run contains stage checkpoints, logs,
and VTF files. In coarse-lipid mode, VTF output represents each CGL particle as a
head endpoint and a tail endpoint connected by one display bond.
