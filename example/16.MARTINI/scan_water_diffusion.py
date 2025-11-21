#!/usr/bin/env python3
"""
Parameter Scanning Script for Water Diffusion Rate
Scans temperature (0.6-1.2, 0.05 step) and thermostat timescale (0.0-0.2, 0.02 step)
"""

import os
import sys
import itertools
import numpy as np
import subprocess
import shutil
from datetime import datetime

# Default configuration
UPSIDE_HOME = os.environ.get('UPSIDE_HOME', '')
if not UPSIDE_HOME:
    print("Error: UPSIDE_HOME environment variable not set!")
    sys.exit(1)

# Parameter ranges
TEMPERATURES = np.arange(0.6, 1.21, 0.05)  # 0.6 to 1.2 with 0.05 step

# Thermostat timescale: 0 to 0.2 with 0.1 gap, include 0.135
# Base range with 0.1 gap: 0.0, 0.1, 0.2 (3 points)
THERMOSTAT_TIMESCALES = np.array([0.0, 0.1, 0.2])
# Check if 0.135 is already included (handle floating-point precision)
target = 0.135
tolerance = 1e-9
is_present = any(np.abs(tau - target) < tolerance for tau in THERMOSTAT_TIMESCALES)
if not is_present:
    THERMOSTAT_TIMESCALES = np.sort(np.append(THERMOSTAT_TIMESCALES, target))
DURATION = 10000  # Simulation duration (steps)
FRAME_INTERVAL = 50  # Frame saving interval
PDB_FILE = "water"  # PDB ID (for prepare_martini.py, it will look for pdb/water.MARTINI.pdb)
RUN_PREFIX = "water_diffusion_scan"  # Prefix for output directories

def run_command(cmd, cwd=None):
    """Run a command and return output"""
    print(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Command failed with error:")
        print(result.stderr)
        sys.exit(1)
    return result.stdout

def prepare_simulation_input():
    """Prepare the initial simulation input file if it doesn't exist"""
    input_file = "inputs/water.up"

    if not os.path.exists(input_file):
        # Create inputs directory if it doesn't exist
        os.makedirs("inputs", exist_ok=True)

        # prepare_martini.py requires only PDB ID as argument (it looks for pdb/{PDB_ID}.MARTINI.pdb internally)
        print(f"Preparing initial input file: {input_file}")
        cmd = f"python3 prepare_martini.py {PDB_FILE}"
        run_command(cmd)

        # Move the generated input file to the expected location
        generated_file = f"outputs/martini_test/{PDB_FILE}.MARTINI.pdb.up"  # prepare_martini.py generates this
        if os.path.exists(generated_file):
            os.rename(generated_file, input_file)
        else:
            # Check if the output file name is different
            print(f"Warning: Expected generated file not found: {generated_file}")
            print("Checking for possible output files...")
            generated_files = [f for f in os.listdir("outputs/martini_test") if f.endswith(".up")]
            if generated_files:
                os.rename(f"outputs/martini_test/{generated_files[0]}", input_file)
                print(f"Found and renamed: {generated_files[0]}")
            else:
                raise FileNotFoundError(f"Could not find the generated .up file")

    return input_file

def generate_simulation_script(temperature, tau):
    """Generate a bash script for a single simulation run"""
    # Create a unique directory name
    temp_str = f"{temperature:.3f}".replace(".", "p")
    tau_str = f"{tau:.3f}".replace(".", "p")
    run_dir = f"{RUN_PREFIX}/T{temp_str}_tau{tau_str}"

    # Create run directory
    os.makedirs(run_dir, exist_ok=True)

    # Simulation command
    upside_binary = os.path.join(UPSIDE_HOME, "obj", "upside")
    input_file = os.path.abspath("inputs/water.up")
    # Output file will be written in the current directory (run_dir)
    output_file = "water.run.up"

    # Build command with thermostat timescale
    cmd = [
        f"{upside_binary}",
        f"--duration {DURATION}",
        f"--frame-interval {FRAME_INTERVAL}",
        f"--temperature {temperature:.3f}",
        f"--thermostat-timescale {tau:.3f}",
        f"--output {output_file}",
        f"{input_file}"
    ]

    # Write the run script
    script_content = f"""#!/bin/bash
# Simulation run script for T={temperature:.3f}, tau={tau:.3f}

set -e

# Source the environment setup and activate virtual environment
source {UPSIDE_HOME}/source.sh
source {UPSIDE_HOME}/.venv/bin/activate

# Change to the run directory
cd "$(dirname "$0")"

# Run the simulation
{' '.join(cmd)}

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
traj = mu.load_upside_traj(traj_file, traj_file)

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
        f.write(f"Temperature: {temperature:.3f} reduced units\\n")
        f.write(f"ThermostatTimescale: {tau:.3f} reduced units\\n")
        f.write(f"DiffusionRate: {{D_cm2s:.9f}} cm²/s\\n")
        f.write(f"MSDSlope: {{slope:.6f}} nm²/ps\\n")
        f.write(f"MSDIntercept: {{intercept:.6f}} nm²\\n")
        f.write(f"BulkWaterAtoms: {{bulk_traj.n_atoms}}\\n")
        f.write(f"TotalWaterAtoms: {{water_traj.n_atoms}}\\n")
    print(f"Completed simulation T={temperature:.3f}, tau={tau:.3f}: D={{D_cm2s:.9f}} cm²/s")
EOF

"""

    script_name = os.path.join(run_dir, "run_simulation.sh")
    with open(script_name, "w") as f:
        f.write(script_content)

    # Make the script executable
    os.chmod(script_name, 0o755)

    return run_dir, script_name

def generate_slurm_script(all_scripts):
    """Generate a SLURM job array script to run all simulations in parallel"""
    # Create the scripts array line with absolute paths
    import os
    all_scripts_absolute = [os.path.abspath(script) for script in all_scripts]
    scripts_array_line = 'scripts=(' + ' '.join(all_scripts_absolute) + ')'

    # Now build the full SLURM script content
    slurm_content = f"""#!/bin/bash
# SLURM job array for water diffusion parameter scan

#SBATCH --job-name=water_diffusion_scan
#SBATCH --account=pi-trsosnic
#SBATCH --partition=caslake
#SBATCH --output=slurm_output/%A_%a.out
#SBATCH --error=slurm_output/%A_%a.err
#SBATCH --array=0-{len(all_scripts)-1}

# Change to the working directory
cd "$(dirname "$0")"

# Execute the corresponding simulation script
{scripts_array_line}
script=${{scripts[$SLURM_ARRAY_TASK_ID]}}
bash "$script"
"""

    # Create slurm_output directory
    os.makedirs("slurm_output", exist_ok=True)

    slurm_file = "run_scan.slurm"
    with open(slurm_file, "w") as f:
        f.write(slurm_content)

    os.chmod(slurm_file, 0o755)

    return slurm_file

def main():
    print("=== WATER DIFFUSION PARAMETER SCAN ===")
    print(f"Temperature range: {TEMPERATURES[0]:.3f} to {TEMPERATURES[-1]:.3f} (step 0.05)")
    # Show the range with special note about 0.135
    print(f"Thermostat timescale range: {THERMOSTAT_TIMESCALES[0]:.3f} to {THERMOSTAT_TIMESCALES[-1]:.3f} (main gap 0.1, including 0.135)")
    print(f"Total combinations: {len(TEMPERATURES) * len(THERMOSTAT_TIMESCALES)}")
    print()

    # Prepare initial input
    input_file = prepare_simulation_input()
    print(f"Prepared input file: {input_file}")
    print()

    # Generate all simulation scripts
    print("Generating simulation scripts...")
    all_scripts = []
    for temperature in TEMPERATURES:
        for tau in THERMOSTAT_TIMESCALES:
            run_dir, script = generate_simulation_script(temperature, tau)
            all_scripts.append(script)

    print(f"Generated {len(all_scripts)} simulation scripts")
    print()

    # Generate SLURM job array script
    print("Generating SLURM job array script...")
    slurm_script = generate_slurm_script(all_scripts)
    print(f"Generated SLURM script: {slurm_script}")
    print()

    print("=== INSTRUCTIONS ===")
    print("1. To run all simulations in parallel on SLURM:")
    print(f"   sbatch {slurm_script}")
    print()
    print("2. To run a single simulation (for testing):")
    print(f"   bash {all_scripts[0]}")
    print()
    print("3. After all simulations complete, run the plotting script:")
    print("   python plot_diffusion_results.py")
    print()
    print("Done!")

if __name__ == "__main__":
    main()
