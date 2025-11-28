import os
import sys
import numpy as np

# --- CONFIGURATION ---
# PLEASE VERIFY THESE PATHS MATCH YOUR CLUSTER EXACTLY
base_dir = "/home/yinhanw/project/yinhan/upside2-md/example/16.MARTINI/water_diffusion"
input_file = "/project/trsosnic/yinhan/upside2-md/example/16.MARTINI/inputs/water.up" 
upside_exec = "/home/yinhanw/project/yinhan/upside2-md/obj/upside"
source_script = "/home/yinhanw/project/yinhan/upside2-md/source.sh"
venv_activate = "/home/yinhanw/project/yinhan/upside2-md/.venv/bin/activate"

# --- PARAMETERS ---
temperatures = np.arange(0.500, 1.000, 0.020)
taus = [0.000, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000, 2.000, 5.000, 10.000]

# --- SCRIPT GENERATION ---

if not os.path.exists(base_dir):
    os.makedirs(base_dir)

run_scripts = []

for T in temperatures:
    for tau in taus:
        folder_name = f"T{T:.3f}_tau{tau:.3f}".replace('.', 'p')
        run_dir = os.path.join(base_dir, folder_name)
        
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
            
        # Embedded Python script with 10000 steps + 10% Skip Logic
        script_content = f"""#!/bin/bash
# Simulation run script for T={T:.3f}, tau={tau:.3f}

set -e

# Source environment
source {source_script}
source {venv_activate}

cd "$(dirname "$0")"

# Run Simulation (Duration restored to 10000)
{upside_exec} --duration 10000 --frame-interval 50 --temperature {T:.3f} --thermostat-timescale {tau:.3f} --output water.run.up {input_file}

# Calculate Diffusion
python - << 'EOF'
import os
import sys
import tables
import mdtraj as md
import numpy as np
import warnings

# Suppress harmless MDTraj warnings
warnings.filterwarnings("ignore", category=UserWarning)

traj_file = "water.run.up"

# --- HELPER: CUSTOM MSD ---
def calculate_msd(traj):
    xyz = traj.xyz
    n_frames = xyz.shape[0]
    msd = np.zeros(n_frames)
    # Loop over time lags (tau)
    for tau in range(1, n_frames):
        diff = xyz[tau:] - xyz[:-tau]
        sq_dist = np.sum(diff**2, axis=-1)
        msd[tau] = np.mean(sq_dist)
    return msd

# --- MAIN LOADING ---
try:
    with tables.open_file(traj_file, 'r') as f:
        # 1. Load Coords (try pos, then xyz)
        if hasattr(f.root.output, 'pos'):
            xyz_raw = f.root.output.pos.read()
        elif hasattr(f.root.output, 'xyz'):
            xyz_raw = f.root.output.xyz.read()
        else:
            raise ValueError("No coordinates found (pos/xyz).")

        # 2. Fix Shape (Frames, Replicas, Atoms, 3) -> (Frames, Atoms, 3)
        if xyz_raw.ndim == 4:
            xyz = xyz_raw[:, 0, :, :]
        elif xyz_raw.ndim == 3:
            xyz = xyz_raw
        else:
            raise ValueError(f"Unexpected shape: {{xyz_raw.shape}}")

        # 3. Load Time
        if hasattr(f.root.output, 'time'):
            time = f.root.output.time.read()
        else:
            time = np.arange(xyz.shape[0])

    # 4. Manual Topology (All atoms = Water)
    n_atoms = xyz.shape[1]
    top = md.Topology()
    chain = top.add_chain()
    for _ in range(n_atoms):
        res = top.add_residue("W", chain)
        top.add_atom("W", md.element.oxygen, res)

    traj = md.Trajectory(xyz, top)
    traj.time = time

except Exception as e:
    print(f"CRITICAL ERROR loading: {{e}}")
    sys.exit(1)

# --- ANALYSIS: BULK FILTERING ---
water_traj = traj
bulk_traj = None

if len(water_traj) > 0 and water_traj.n_atoms > 1:
    # Strategy A: Neighbor Density
    neighbor_cutoff = 0.8 # nm
    bulk_threshold = 3
    
    # Slice [0] to get the first frame for neighbor calc
    neighbors_list = md.compute_neighbors(water_traj[0], cutoff=neighbor_cutoff, 
                                          query_indices=range(water_traj.n_atoms))
    neighbors = neighbors_list[0]
    neighbor_counts = [len(nbrs) for nbrs in neighbors]
    bulk_indices = [i for i, count in enumerate(neighbor_counts) if count >= bulk_threshold]

    if len(bulk_indices) > 10:
        bulk_traj = water_traj.atom_slice(bulk_indices)
        print(f"Selected {{len(bulk_indices)}} bulk atoms (Neighbor method).")
    else:
        # Strategy B: Center of Mass Core (Fallback)
        xyz0 = water_traj.xyz[0]
        com = np.mean(xyz0, axis=0)
        dists = np.linalg.norm(xyz0 - com, axis=1)
        
        # Keep inner 80%
        sorted_indices = np.argsort(dists)
        n_keep = int(0.8 * len(dists))
        bulk_indices = sorted_indices[:n_keep]
        
        bulk_traj = water_traj.atom_slice(bulk_indices)
        print(f"Selected {{len(bulk_indices)}} core atoms (Center-of-Mass method).")

if bulk_traj is None:
    print("Warning: Bulk selection failed. Using all atoms.")
    bulk_traj = water_traj

# --- MSD & DIFFUSION ---
msd = calculate_msd(bulk_traj)

# Determine start frame (Skip first 10%)
n_frames = len(traj.time)
start_frame = int(n_frames * 0.10)

if n_frames > (start_frame + 2):
    # Linear fit on the latter 90%
    fit_time = traj.time[start_frame:]
    fit_msd = msd[start_frame:]
    
    slope, intercept = np.polyfit(fit_time, fit_msd, 1)
    
    # D = Slope / 6
    # Conversion: nm^2/ps -> cm^2/s using 1e-5 scaling
    D_cm2s = (slope / 6.0) * 1e-5 

    with open("diffusion_rate.txt", "w") as f:
        f.write(f"Temperature: {T:.3f} reduced units\\n")
        f.write(f"ThermostatTimescale: {tau:.3f} reduced units\\n")
        f.write(f"DiffusionRate: {{D_cm2s:.9f}} cm^2/s\\n")
        f.write(f"MSDSlope: {{slope:.6f}} nm^2/ps\\n")
        f.write(f"AnalysisStartFrame: {{start_frame}}\\n")
        f.write(f"BulkMethod: {{'Neighbors' if len(bulk_indices) > 10 and 'neighbor_counts' in locals() and len([x for x in neighbor_counts if x>=3]) > 10 else 'CenterOfMass'}}\\n")
    print(f"Completed T={T:.3f}, tau={tau:.3f}: D={{D_cm2s:.9f}} cm^2/s (Skipped first {{start_frame}} frames)")
else:
    print("Not enough frames.")
EOF
"""
        
        # Write script file
        script_path = os.path.join(run_dir, "run_simulation.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        run_scripts.append(script_path)

print(f"Generated {len(run_scripts)} simulation scripts.")

# --- SLURM SUBMISSION GENERATOR ---
tasks_file = os.path.join(base_dir, "tasks.txt")
with open(tasks_file, "w") as f:
    for script in run_scripts:
        f.write(f"{script}\n")

slurm_content = f"""#!/bin/bash
#SBATCH --job-name=water_diff_scan
#SBATCH --output=water_scan_%A_%a.out
#SBATCH --error=water_scan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --array=1-{len(run_scripts)}

TASK_LIST="{tasks_file}"
SCRIPT_PATH=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" "$TASK_LIST")

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