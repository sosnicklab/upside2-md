import os
import sys
import numpy as np

# --- CONFIGURATION ---
base_dir = "/home/yinhanw/project/yinhan/upside2-md/example/16.MARTINI/water_diffusion"
input_file = "/project/trsosnic/yinhan/upside2-md/example/16.MARTINI/inputs/water.up"
upside_exec = "/home/yinhanw/project/yinhan/upside2-md/obj/upside"
source_script = "/home/yinhanw/project/yinhan/upside2-md/source.sh"
venv_activate = "/home/yinhanw/project/yinhan/upside2-md/.venv/bin/activate"

# --- PARAMETERS ---
# (Restore your full parameter ranges here)
temperatures = np.arange(0.500, 1.000, 0.020)
taus = [0.000, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000, 2.000, 5.000, 10.000]

# --- SCRIPT GENERATION ---

if not os.path.exists(base_dir):
    os.makedirs(base_dir)

run_scripts = []

for T in temperatures:
    for tau in taus:
        # Folder setup
        folder_name = f"T{T:.3f}_tau{tau:.3f}".replace('.', 'p')
        run_dir = os.path.join(base_dir, folder_name)
        
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
            
        # 1. Generate the individual simulation script
        # We inject {input_file} into the Python block to fix the "NoSuchNodeError"
        script_content = f"""#!/bin/bash
# Simulation run script for T={T:.3f}, tau={tau:.3f}

set -e

# Source environment
source {source_script}
source {venv_activate}

cd "$(dirname "$0")"

# Run Simulation
{upside_exec} --duration 10000 --frame-interval 50 --temperature {T:.3f} --thermostat-timescale {tau:.3f} --output water.run.up {input_file}

# Calculate Diffusion
python - << 'EOF'
import os
import sys
sys.path.append(os.path.abspath('../../../py'))
import mdtraj as md
import mdtraj_upside as mu
import numpy as np

# Load trajectory
traj_file = "water.run.up"
top_file = "{input_file}"  # <--- FIX: Point to input file for topology

# Load Upside trajectory
traj = mu.load_upside_traj(traj_file, top_file)

# Select water
water_sel = traj.top.select("resname W")
water_traj = traj.atom_slice(water_sel)

# Filter bulk (Non-periodic logic)
bulk_traj = water_traj
if len(water_traj) > 0 and water_traj.n_atoms > 1:
    neighbor_cutoff = 0.5 
    bulk_neighbor_threshold = 3 
    
    # Neighbor calculation
    neighbors = md.compute_neighbors(water_traj, cutoff=neighbor_cutoff,
                                     query_indices=range(water_traj.n_atoms),
                                     haystack_indices=None, frame=0)
    neighbor_counts = [len(nbrs) for nbrs in neighbors]
    bulk_indices = [i for i, count in enumerate(neighbor_counts) if count >= bulk_neighbor_threshold]

    if len(bulk_indices) > 0:
        bulk_traj = water_traj.atom_slice(bulk_indices)
        print(f"Selected {{len(bulk_indices)}} bulk water atoms.")
    else:
        print("Warning: No bulk water atoms found, using all atoms")

# MSD and Diffusion
msd = md.compute_msd(bulk_traj, select="all", window=100)
times = traj.time

if len(times) >= 2:
    slope, intercept = np.polyfit(times[1:], msd[1:], 1)
    D_cm2s = (slope / 6.0) * 1e-5  # nm^2/ps to cm^2/s

    with open("diffusion_rate.txt", "w") as f:
        f.write(f"Temperature: {T:.3f} reduced units\\n")
        f.write(f"ThermostatTimescale: {tau:.3f} reduced units\\n")
        f.write(f"DiffusionRate: {{D_cm2s:.9f}} cm^2/s\\n")
        f.write(f"MSDSlope: {{slope:.6f}} nm^2/ps\\n")
    print(f"Completed T={T:.3f}, tau={tau:.3f}: D={{D_cm2s:.9f}} cm^2/s")
EOF
"""
        
        script_path = os.path.join(run_dir, "run_simulation.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        run_scripts.append(script_path)

print(f"Generated {len(run_scripts)} simulation scripts.")

# --- 2. GENERATE ROBUST SLURM SUBMISSION ---

# Write the list of all tasks to a text file
tasks_file = os.path.join(base_dir, "tasks.txt")
with open(tasks_file, "w") as f:
    for script in run_scripts:
        f.write(f"{script}\n")

# Generate the SLURM script that reads from tasks.txt
# Array is 1-based to match 'sed' line numbers
slurm_content = f"""#!/bin/bash
#SBATCH --job-name=water_diff_scan
#SBATCH --output=water_scan_%A_%a.out
#SBATCH --error=water_scan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --array=1-{len(run_scripts)}

# Set the file containing the list of scripts
TASK_LIST="{tasks_file}"

# Get the script path for this task ID (read specific line from file)
SCRIPT_PATH=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" "$TASK_LIST")

# Check if script exists and execute
if [ -f "$SCRIPT_PATH" ]; then
    echo "Running task $SLURM_ARRAY_TASK_ID: $SCRIPT_PATH"
    "$SCRIPT_PATH"
else
    echo "Error: Script not found at $SCRIPT_PATH"
    exit 1
fi
"""

slurm_script_path = os.path.join(base_dir, "run_scan.slurm")
with open(slurm_script_path, "w") as f:
    f.write(slurm_content)

print(f"Generated master SLURM script: {slurm_script_path}")
print(f"Task list saved to: {tasks_file}")