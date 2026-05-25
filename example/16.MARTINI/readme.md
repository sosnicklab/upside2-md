This example shows:
1. how to run the hybrid dry-MARTINI + Upside membrane workflow
2. how to launch the workflow through the python wrapper `run.py`
3. how to generate VTF output for visualization from the stage workflow

Similar to the other examples, you can use the python script `run.py` in this folder:

source ../../.venv/bin/activate
source ../../source.sh
python run.py

Edit the settings near the top of `run.py` to change the PDB id, run
directory, salt concentration, lipid representation, or stage lengths.  The
MARTINI workflow helper scripts live under `../../py/`.

The workflow now defaults to the same thermostat timescale as the standard Upside examples (`tau = 5.0`) unless you override `THERMOSTAT_TIMESCALE`.
Its hybrid stages still use smaller explicit timesteps than the standard examples for stability, so any repo-wide physical time calibration borrowed from those examples should be revalidated before applying it literally to stage `7.0`.
Fresh workflow runs now execute a rigid-protein `6.0` NPT box-relaxation stage before handing coordinates and box dimensions into `7.0`.  Adjust its length with `EQ_60_NSTEPS` when you need a shorter or longer packing pass.
Production stage `7.0` defaults to NVT (`PROD_70_NPT_ENABLE=0`).  That is
useful for stability comparisons, but membrane packing can drift away from the
dry-MARTINI reference when production NPT is disabled.

You can also run the shell workflow directly:

bash run_sim_1rkl.sh

The default results are stored in `outputs/martini_1rkl_hybrid`.  The workflow
writes stage checkpoints, logs, and VTF files there.

To test packing with production NPT enabled, override:

`PROD_70_NPT_ENABLE=1`

The workflow uses the semi-isotropic barostat settings from
`run_sim_hybrid.sh`, so this keeps `x/y` coupled while using the configured
normal-direction compressibility.

Generated VTF files render each CG lipid as two visible particles connected by one bond: a hydrophilic endpoint and a hydrophobic endpoint.  Their separation follows the source DOPC head-to-tail reference span rather than the hidden runtime orientation-site bond length.
