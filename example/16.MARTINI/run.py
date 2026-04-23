import os
import subprocess as sp
from pathlib import Path


#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id = "1rkl"
runtime_pdb_id = f"{pdb_id}_aabb"
run_dir = "./outputs/martini_test_1rkl_aabb"

salt_molar = 0.15
protein_lipid_cutoff = 4.5
xy_scale = 1.0
box_padding_xy = 0.0
box_padding_z = 20.0
temperature = 0.8647

min_60_max_iter = 500
min_61_max_iter = 500
eq_62_nsteps = 500
eq_63_nsteps = 500
eq_64_nsteps = 500
eq_65_nsteps = 500
eq_66_nsteps = 500
prod_70_nsteps = 10000

prod_frame_steps = 50
prod_70_npt_enable = 0


#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

script_dir = Path(__file__).resolve().parent
workflow = script_dir / "run_sim_1rkl.sh"

if not workflow.exists():
    raise FileNotFoundError(f"Workflow script not found: {workflow}")

env_updates = {
    "RUNTIME_PDB_ID": runtime_pdb_id,
    "RUN_DIR": run_dir,
    "SALT_MOLAR": salt_molar,
    "PROTEIN_LIPID_CUTOFF": protein_lipid_cutoff,
    "XY_SCALE": xy_scale,
    "BOX_PADDING_XY": box_padding_xy,
    "BOX_PADDING_Z": box_padding_z,
    "TEMPERATURE": temperature,
    "MIN_60_MAX_ITER": min_60_max_iter,
    "MIN_61_MAX_ITER": min_61_max_iter,
    "EQ_62_NSTEPS": eq_62_nsteps,
    "EQ_63_NSTEPS": eq_63_nsteps,
    "EQ_64_NSTEPS": eq_64_nsteps,
    "EQ_65_NSTEPS": eq_65_nsteps,
    "EQ_66_NSTEPS": eq_66_nsteps,
    "PROD_70_NSTEPS": prod_70_nsteps,
    "PROD_FRAME_STEPS": prod_frame_steps,
    "PROD_70_NPT_ENABLE": prod_70_npt_enable,
}

run_env = os.environ.copy()
run_env.update({key: str(value) for key, value in env_updates.items()})

cmd = ["bash", str(workflow), f"PDB_ID={pdb_id}"]

print("Running hybrid dry-MARTINI + Upside example...")
print("Command:")
print(" ".join(cmd))
print("Environment overrides:")
for key in sorted(env_updates):
    print(f"  {key}={env_updates[key]}")

sp.check_call(cmd, cwd=str(script_dir), env=run_env)
