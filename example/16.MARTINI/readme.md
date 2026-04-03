This example shows:
1. how to run the hybrid dry-MARTINI + Upside membrane workflow for `1rkl`
2. how to launch the workflow through the python wrapper `run.py`
3. how to generate VTF output for visualization from the stage workflow

Similar to the other examples, you can use the python script `run.py` in this folder:

source ../../source.sh
python run.py

The script `run.py` is a thin wrapper around `run_sim_1rkl.sh`.  Edit the settings near the top of `run.py` if you want to change the PDB id, run directory, salt concentration, or stage lengths.

You can also run the shell workflow directly:

bash run_sim_1rkl.sh

The default results are stored in `outputs/martini_test_1rkl_hybrid`.  The workflow writes stage checkpoints, logs, and VTF files there.
