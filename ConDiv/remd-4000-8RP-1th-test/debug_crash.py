import os
import sys
import subprocess as sp
import numpy as np

# --- CONFIGURATION ---
# Use one of the proteins that failed
pdb_code = "1s29" 
# Path to your upside_config script
config_script = os.path.join(os.getcwd(), "../../py/upside_config.py")
# Path to upside binary
upside_bin = os.path.join(os.getcwd(), "../../obj/upside")

# Define paths (Adjust if your folder structure is different)
input_fasta = f"upside_input/{pdb_code}.fasta"
input_init  = f"upside_input/{pdb_code}.initial.pkl"
output_h5   = f"debug_{pdb_code}.h5"
output_log  = f"debug_{pdb_code}.log"

# Mock parameters from your ff_2.1 folder
param_dir = "../../parameters/ff_2.1"
common_dir = "../../parameters/common"

print(f"--- STARTING DEBUG RUN FOR {pdb_code} ---")

# 1. Generate Configuration (mimicking ConDiv.py)
print("Generating configuration...")
cmd = [
    sys.executable, config_script,
    f"--fasta={input_fasta}",
    f"--output={output_h5}",
    f"--initial-structure={input_init}",
    
    # Force Field Params
    f"--hbond-energy={param_dir}/hbond.h5",
    f"--rama-sheet-mixing-energy={param_dir}/sheet",
    f"--environment-potential={param_dir}/environment.h5",
    f"--rotamer-placement={param_dir}/sidechain.h5",
    f"--rotamer-interaction={param_dir}/sidechain.h5",
    "--dynamic-rotamer-1body",
    
    # Common Params
    f"--rama-library={common_dir}/rama.dat",
    f"--reference-state-rama={common_dir}/rama_reference.pkl",
    
    # THE RESTRAINTS (The suspected cause)
    "--restraint-group=0-50", # Restrain first 50 residues
    "--restraint-spring=0.111"
]

try:
    sp.check_output(cmd, stderr=sp.STDOUT)
    print("Configuration generated successfully.")
except sp.CalledProcessError as e:
    print("FAILED to generate config:")
    print(e.output.decode())
    sys.exit(1)

# 2. Run Simulation
print("Running Upside Engine...")
run_cmd = [
    upside_bin,
    "--duration", "100",      # Short run to test crash
    "--frame-interval", "10",
    "--temperature", "0.8",
    output_h5
]

with open(output_log, "w") as log:
    # We pipe stdout to file but verify return code
    p = sp.Popen(run_cmd, stdout=log, stderr=sp.STDOUT)
    ret = p.wait()

print(f"Engine finished with exit code: {ret}")

# 3. Check Log
print("\n--- LAST 20 LINES OF LOG ---")
os.system(f"tail -n 20 {output_log}")

if ret == -11:
    print("\nRESULT: Segmentation Fault (SIGSEGV). This confirms memory corruption.")
    print("Likely cause: Float64/Float32 mismatch in H5 file.")
elif ret == 0:
    print("\nRESULT: Success! The debug run worked.")
    print("If this worked but ConDiv fails, the issue is ConDiv's parallelization.")
else:
    print(f"\nRESULT: Crashed with code {ret}")
