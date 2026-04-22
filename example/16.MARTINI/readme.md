This example shows:
1. how to run the hybrid dry-MARTINI + Upside membrane workflow for `1rkl`
2. how to launch the workflow through the python wrapper `run.py`
3. how to generate VTF output for visualization from the stage workflow

Similar to the other examples, you can use the python script `run.py` in this folder:

source ../../source.sh
python run.py

The script `run.py` is a thin wrapper around `run_sim_1rkl.sh`.  Edit the settings near the top of `run.py` if you want to change the PDB id, run directory, salt concentration, or stage lengths.
The MARTINI workflow helper scripts used by `run_sim_1rkl.sh` live under `../../py/`.
For the bilayer input, use a lipid bilayer structure generated from CHARMM-GUI; note that while `DOPC.pdb` is used as an example, any valid bilayer structure from CHARMM-GUI is compatible.

The workflow now defaults to the same thermostat timescale as the standard Upside examples (`tau = 5.0`) unless you override `THERMOSTAT_TIMESCALE`.
Its hybrid stages still use smaller explicit timesteps than the standard examples for stability, so any repo-wide physical time calibration borrowed from those examples should be revalidated before applying it literally to stage `7.0`.
Production stage `7.0` still defaults to NVT (`PROD_70_NPT_ENABLE=0`).  That is useful for stability comparisons, but it also means membrane packing can drift away from the dry-MARTINI reference even when the interface coupling is otherwise reasonable.

You can also run the shell workflow directly:

bash run_sim_1rkl.sh

The default results are stored in `outputs/martini_test_1rkl_hybrid`.  The workflow writes stage checkpoints, logs, and VTF files there.

To test packing with production NPT enabled, override:

`PROD_70_NPT_ENABLE=1`

The workflow already uses the semi-isotropic barostat settings defined in `run_sim_1rkl.sh`, so this keeps `x/y` coupled while leaving `z` compressibility at the configured normal setting.

## Packing Analysis

Use `analyze_packing.py` on a completed `stage_7.0.up` file to measure:

- XY box area
- area per lipid (APL) from `PO4`-identified lipid molecules
- `PO4` leaflet thickness
- tail-bond orientational order relative to the bilayer normal

Example:

`.venv/bin/python analyze_packing.py outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up --json-out outputs/martini_test_1rkl_hybrid/packing_summary.json --csv-out outputs/martini_test_1rkl_hybrid/packing_timeseries.csv`

Recommended calibration loop:

1. Run the current hybrid workflow with the intended production settings.
2. Analyze `stage_7.0.up` with `analyze_packing.py`.
3. Compare APL / thickness / order against the dry-MARTINI reference.
4. Only after packing is matched, revisit diffusion and interface-scale tuning.
