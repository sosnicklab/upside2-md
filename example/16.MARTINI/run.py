import os, sys
import subprocess as sp
from pathlib import Path


upside_path = Path(os.environ.get("UPSIDE_HOME", Path(__file__).resolve().parents[2])).resolve()
upside_utils_dir = upside_path / "py"
sys.path.insert(0, str(upside_utils_dir))

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id             = "1rkl"
pdb_dir            = "./pdb"
opm_pdb            = "{}/{}.pdb".format(pdb_dir, pdb_id)
sim_id             = "martini_{}_hybrid".format(pdb_id)
work_dir           = "./"

lipid_model        = "single_particle"  # "single_particle" for CGL/CGLD, "cg_lipid" for 14-bead DOPC
lipid_composition  = dict(DOPC=1.0)

salt_molar         = 0.15
protein_lipid_cutoff = 4.5
protein_lipid_min_gap = 4.5
ion_cutoff         = 10.0
xy_scale           = 1.0
box_padding_xy     = 0.0
box_padding_z      = 20.0
protein_placement_mode = "embed"
protein_orientation_mode = "input"
protein_surface_gap = 6.0

temperature        = 0.8647
eq_time_step       = 0.010
prod_time_step     = 0.002
eq_frame_steps     = 1000
prod_frame_steps   = 50
prod_70_npt_enable = 0

min_60_max_iter    = 500
min_61_max_iter    = 0
min_70_max_iter    = 500
eq_60_nsteps       = 500
eq_62_nsteps       = 500
eq_63_nsteps       = 500
eq_64_nsteps       = 500
eq_65_nsteps       = 500
eq_66_nsteps       = 500
prod_70_burnin_nsteps = 40000
prod_70_nsteps     = 10000

prep_seed          = 2026
randomseed         = 1

submit_slurm       = False
account            = ""
partition          = ""
run_time           = "36:00:00"
ntasks_per_node    = 1


#----------------------------------------------------------------------
# Set the path and filename
#----------------------------------------------------------------------

script_dir = Path(__file__).resolve().parent
work_dir_path = Path(work_dir).expanduser()
if not work_dir_path.is_absolute():
    work_dir_path = script_dir / work_dir_path
run_dir = work_dir_path / "outputs" / sim_id
log_dir = run_dir / "logs"
slurm_dir = run_dir / "slurm"

for direc in (run_dir, log_dir, slurm_dir):
    if not direc.exists():
        direc.mkdir(parents=True)

lipid_resolution_by_model = {
    "single_particle": "coarse",
    "cg_lipid": "full",
}
if lipid_model not in lipid_resolution_by_model:
    raise ValueError("lipid_model must be 'single_particle' or 'cg_lipid'")
lipid_resolution = lipid_resolution_by_model[lipid_model]

if set(lipid_composition) != {"DOPC"} or abs(lipid_composition["DOPC"] - 1.0) > 1e-6:
    raise ValueError("This example currently supports only 100% DOPC lipid composition")

runtime_pdb_id = "{}_{}_hybrid".format(pdb_id, lipid_model)
bilayer_pdb = upside_path / "parameters" / "dryMARTINI" / "DOPC.pdb"
prep_script = upside_path / "py" / "martini_prepare_system.py"
param_script = upside_path / "py" / "martini_gen_params.py"
extract_vtf_script = upside_path / "py" / "martini_extract_vtf.py"

opm_pdb_path = (script_dir / opm_pdb).resolve()
if not opm_pdb_path.exists():
    raise FileNotFoundError("OPM-style input PDB not found: {}".format(opm_pdb_path))


#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

python = upside_path / ".venv" / "bin" / "python3"
if not python.exists():
    python = Path(sys.executable)

env = os.environ.copy()
env["UPSIDE_HOME"] = str(upside_path)
env["PATH"] = "{}{}{}".format(upside_path / "obj", os.pathsep, env.get("PATH", ""))
env["PYTHONPATH"] = "{}{}{}".format(upside_utils_dir, os.pathsep, env.get("PYTHONPATH", ""))

required_params = [
    upside_path / "parameters" / "dryMARTINI" / "particle.h5",
    upside_path / "parameters" / "dryMARTINI" / "sidechain.h5",
    upside_path / "parameters" / "dryMARTINI" / "interlipid.h5",
]
if lipid_resolution == "coarse":
    required_params.append(upside_path / "parameters" / "dryMARTINI" / "dopc.h5")

workflow_cmd = [
    str(python), str(prep_script), "run-hybrid-workflow",
    "--pdb-id", pdb_id,
    "--runtime-pdb-id", runtime_pdb_id,
    "--upside-home", str(upside_path),
    "--run-dir", str(run_dir),
    "--protein-aa-pdb", str(opm_pdb_path),
    "--bilayer-pdb", str(bilayer_pdb),
    "--extract-vtf-script", str(extract_vtf_script),
    "--salt-molar", str(salt_molar),
    "--protein-lipid-cutoff", str(protein_lipid_cutoff),
    "--protein-lipid-min-gap", str(protein_lipid_min_gap),
    "--ion-cutoff", str(ion_cutoff),
    "--xy-scale", str(xy_scale),
    "--box-padding-xy", str(box_padding_xy),
    "--box-padding-z", str(box_padding_z),
    "--protein-placement-mode", protein_placement_mode,
    "--protein-orientation-mode", protein_orientation_mode,
    "--protein-surface-gap", str(protein_surface_gap),
    "--temperature", str(temperature),
    "--min-60-max-iter", str(min_60_max_iter),
    "--min-61-max-iter", str(min_61_max_iter),
    "--min-70-max-iter", str(min_70_max_iter),
    "--eq-60-nsteps", str(eq_60_nsteps),
    "--eq-62-nsteps", str(eq_62_nsteps),
    "--eq-63-nsteps", str(eq_63_nsteps),
    "--eq-64-nsteps", str(eq_64_nsteps),
    "--eq-65-nsteps", str(eq_65_nsteps),
    "--eq-66-nsteps", str(eq_66_nsteps),
    "--prod-70-burnin-nsteps", str(prod_70_burnin_nsteps),
    "--prod-70-nsteps", str(prod_70_nsteps),
    "--eq-time-step", str(eq_time_step),
    "--prod-time-step", str(prod_time_step),
    "--eq-frame-steps", str(eq_frame_steps),
    "--prod-frame-steps", str(prod_frame_steps),
    "--prod-70-npt-enable", str(prod_70_npt_enable),
    "--prep-seed", str(prep_seed),
    "--seed", str(randomseed),
    "--lipid-resolution", lipid_resolution,
]

print("MARTINI hybrid example")
print("  OPM PDB: {}".format(opm_pdb_path))
print("  lipid composition: {}".format(lipid_composition))
print("  lipid model: {} ({})".format(lipid_model, lipid_resolution))
print("  run directory: {}".format(run_dir))

if submit_slurm:
    slurm_script = slurm_dir / "{}.sbatch".format(sim_id)
    slurm_log = log_dir / "{}.slurm.out".format(sim_id)
    with slurm_script.open("w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name={}\n".format(sim_id))
        f.write("#SBATCH --output={}\n".format(slurm_log.resolve()))
        f.write("#SBATCH --time={}\n".format(run_time))
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks-per-node={}\n".format(ntasks_per_node))
        if account:
            f.write("#SBATCH --account={}\n".format(account))
        if partition:
            f.write("#SBATCH --partition={}\n".format(partition))
        f.write("set -euo pipefail\n\n")
        f.write("PROJECT_ROOT={}\n".format(upside_path))
        f.write("if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi\n")
        f.write("if command -v module >/dev/null 2>&1; then\n")
        f.write("  module load python/3.11.9 || true\n")
        f.write("  module load cmake || true\n")
        f.write("  module load openmpi || true\n")
        f.write("  module load hdf5/1.14.3 || true\n")
        f.write("fi\n")
        f.write("if [ -f \"$PROJECT_ROOT/.venv/bin/activate\" ]; then source \"$PROJECT_ROOT/.venv/bin/activate\"; fi\n")
        f.write("export UPSIDE_HOME=\"$PROJECT_ROOT\"\n")
        f.write("export PATH=\"$PROJECT_ROOT/obj:$PATH\"\n")
        f.write("export PYTHONPATH=\"$PROJECT_ROOT/py${PYTHONPATH:+:$PYTHONPATH}\"\n")
        f.write("export UPSIDE_SKIP_SOURCE_SH=1\n")
        f.write("cd {}\n\n".format(script_dir))
        f.write("missing=0\n")
        for param in required_params:
            f.write("[ -f {} ] || missing=1\n".format(param))
        f.write("if [ \"$missing\" = \"1\" ]; then\n")
        f.write("  python3 {} --upside-home \"$PROJECT_ROOT\"\n".format(param_script))
        f.write("fi\n")
        f.write("{}\n".format(" ".join(workflow_cmd)))

    print("Submitting {}".format(slurm_script))
    sp.check_call(["sbatch", str(slurm_script)])
else:
    if any(not path.exists() for path in required_params):
        print("Generating missing MARTINI parameter files")
        sp.check_call([str(python), str(param_script), "--upside-home", str(upside_path)], env=env)

    print("Running locally")
    print(" ".join(workflow_cmd))
    sp.check_call(workflow_cmd, cwd=str(script_dir), env=env)
