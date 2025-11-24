import os
import subprocess
import numpy as np

# Configuration
base_dir = "/home/yinhanw/project/yinhan/upside2-md/example/16.MARTINI/water_diffusion"
input_file = "/project/trsosnic/yinhan/upside2-md/example/16.MARTINI/inputs/water.up"
upside_exec = "/home/yinhanw/project/yinhan/upside2-md/obj/upside"
source_script = "/home/yinhanw/project/yinhan/upside2-md/source.sh"
venv_activate = "/home/yinhanw/project/yinhan/upside2-md/.venv/bin/activate"

# Simulation parameters
temperatures = np.arange(0.500, 1.000, 0.020)
taus = [0.000, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000, 2.000, 5.000, 10.000]

# Create base directory if it doesn't exist
if not os.path.exists(base_dir):
    os.makedirs(base_dir)

for T in temperatures:
    for tau in taus:
        # Format folder name
        folder_name = f"T{T:.3f}_tau{tau:.3f}".replace('.', 'p')
        run_dir = os.path.join(base_dir, folder_name)
        
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
            
        # Define the script content
        script_content = f"""#!/bin/bash
# Simulation run script for T={T:.3f}, tau={tau:.3f}

set -e

# Source the environment setup and activate virtual environment
source {source_script}
source {venv_activate}

# Change to the run directory
cd "$(dirname "$0")"

# Run the simulation
{upside_exec} --duration 10000 --frame-interval 50 --temperature {T:.3f} --thermostat-timescale {tau:.3f} --output water.run.up {input_file}

# Calculate diffusion rate
python - << 'EOF'
import os
import sys
sys.path.append(os.path.abspath('../../../py'))
import mdtraj as md
import mdtraj_upside as mu
import numpy as np

# Load trajectory
traj_file = "water.run.up"
# FIX: Use the input file for topology to avoid NoSuchNodeError
top_file = "{input_file}" 

# Load Upside trajectory (top_file provides the 'input' group/sequence)
traj = mu.load_upside_traj(traj_file, top_file)

# Select water atoms
water_sel = traj.top.select("resname W")
water_traj = traj.atom_slice(water_sel)

# Filter out edge particles (system is non-periodic, so edge particles have fewer neighbors)
bulk_traj = water_traj
if len(water_traj) > 0 and water_traj.n_atoms > 1:
    # Parameters for bulk water selection
    neighbor_cutoff = 0.5  # nm (5 Å) - appropriate for MARTINI water
    bulk_neighbor_threshold = 3  # At least 3 neighbors to be considered bulk

    # Compute neighbors for all atoms in the first frame
    neighbors = md.compute_neighbors(water_traj, cutoff=neighbor_cutoff,
                                     query_indices=range(water_traj.n_atoms),
                                     haystack_indices=None, frame=0)

    # Count neighbors for each atom
    neighbor_counts = [len(nbrs) for nbrs in neighbors]

    # Select bulk atoms (atoms with sufficient neighbors)
    bulk_indices = [i for i, count in enumerate(neighbor_counts) if count >= bulk_neighbor_threshold]

    if len(bulk_indices) > 0:
        bulk_traj = water_traj.atom_slice(bulk_indices)
        print(f"Selected {{len(bulk_indices)}} bulk water atoms (excluded {{water_traj.n_atoms - len(bulk_indices)}} edge atoms)")
    else:
        print("Warning: No bulk water atoms found, using all atoms")

# Calculate MSD
msd = md.compute_msd(bulk_traj, select="all", window=100)
times = traj.time

# Calculate diffusion rate using Einstein relation (D = slope/6)
if len(times) >= 2:
    slope, intercept = np.polyfit(times[1:], msd[1:], 1)
    D_cm2s = (slope / 6.0) * 1e-5  # Convert nm²/ps to cm²/s

    # Write result
    with open("diffusion_rate.txt", "w") as f:
        f.write(f"Temperature: {T:.3f} reduced units\\n")
        f.write(f"ThermostatTimescale: {tau:.3f} reduced units\\n")
        f.write(f"DiffusionRate: {{D_cm2s:.9f}} cm²/s\\n")
        f.write(f"MSDSlope: {{slope:.6f}} nm²/ps\\n")
        f.write(f"MSDIntercept: {{intercept:.6f}} nm²\\n")
        f.write(f"BulkWaterAtoms: {{bulk_traj.n_atoms}}\\n")
        f.write(f"TotalWaterAtoms: {{water_traj.n_atoms}}\\n")
    print(f"Completed simulation T={T:.3f}, tau={tau:.3f}: D={{D_cm2s:.9f}} cm²/s")
EOF
"""
        
        # Write the script
        script_path = os.path.join(run_dir, "run_simulation.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_path, 0o755)
        
        print(f"Generated script for T={T:.3f}, tau={tau:.3f} in {run_dir}")
        
        # Optional: Submit to slurm or run locally
        # subprocess.run(["sbatch", script_path])