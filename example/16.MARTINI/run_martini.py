import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb
import h5py as h5

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

# Get the Python interpreter path
python_path = sys.executable

# Import UPSIDE utilities
import run_upside as ru
from advanced_config import write_external_pairs_potential

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

sim_id         = 'martini_test'
T              = 0.8  # Original temperature
duration       = 1000
frame_interval = 50
base_dir       = './'

#----------------------------------------------------------------------
## Initialization
#----------------------------------------------------------------------

input_dir  = "{}/inputs".format(base_dir)
output_dir = "{}/outputs".format(base_dir)
run_dir    = "{}/{}".format(output_dir, sim_id) 

make_dirs = [input_dir, output_dir, run_dir]
for direc in make_dirs:
    if not os.path.exists(direc):
        os.makedirs(direc)

#----------------------------------------------------------------------
## Convert MARTINI PDB to UPSIDE format
#----------------------------------------------------------------------

print("Converting MARTINI PDB to UPSIDE format...")
cmd = f"python martini_to_upside.py input.pdb {input_dir}/martini"
try:
    sp.check_call(cmd, shell=True)
except sp.CalledProcessError as e:
    print(f"Error converting MARTINI PDB: {e}")
    sys.exit(1)

#----------------------------------------------------------------------
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + 'ff_2.1/'  # Using standard force field

# options
print("Configuring...")
fasta = "{}/martini.fasta".format(input_dir)  # You'll need to create this

# Create a minimal fasta file with the correct number of residues
with open(fasta, 'w') as f:
    f.write("> MARTINI\n")
    # Add appropriate number of 'A' residues based on your system
    f.write("A" * 8)  # 8 water beads

# First create the basic file structure using upside_config.py
temp_file = "{}/martini_temp.up".format(input_dir)

config_args = [
    python_path,
    os.path.join(upside_utils_dir, 'upside_config.py'),
    '--fasta', fasta,
    '--output', temp_file,
    '--rama-library', param_dir_common + "rama.dat",
    '--rama-sheet-mixing-energy', param_dir_ff + "sheet",
    '--reference-state-rama', param_dir_common + "rama_reference.pkl",
    '--hbond-energy', param_dir_ff + "hbond.h5",  # Use the standard hbond parameter file
    '--rotamer-placement', param_dir_ff + "sidechain.h5",
    '--dynamic-rotamer-1body',
    '--rotamer-interaction', param_dir_ff + "sidechain.h5",
    '--environment-potential', param_dir_ff + "environment.h5",
    '--bb-environment-potential', param_dir_ff + "bb_env.dat",
    '--no-backbone',  # Disable backbone interactions since we only have water beads
]

try:
    sp.check_call(config_args)
except sp.CalledProcessError as e:
    print(f"Error running upside_config: {e}")
    sys.exit(1)

# Read initial positions from PDB
initial_positions = []
with open("input.pdb", 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
initial_positions = np.array(initial_positions)
print(f"Initial positions shape: {initial_positions.shape}")
print(f"Initial positions:\n{initial_positions}")

# Open the temporary file for modification
with tb.open_file(temp_file, 'r+') as t:
    # Create martini_potential group
    martini_group = t.create_group(t.root.input.potential, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'pair'
    martini_group._v_attrs.epsilon = 1.715293655572  # kJ/mol (unit converted for UPSIDE)
    martini_group._v_attrs.sigma = 4.7    # Å (converted from 0.47 nm)
    martini_group._v_attrs.lj_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.coul_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.dielectric = 15.0  # MARTINI water dielectric
    martini_group._v_attrs.n_types = 1  # Only one type of particle (water)
    martini_group._v_attrs.n_params = 4  # [epsilon, sigma, charge1, charge2]
    martini_group._v_attrs.cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.cache_buffer = 1.0  # Default cache buffer
    
    # Get the expected number of atoms from the FASTA sequence
    with open(fasta, 'r') as f:
        fasta_content = f.read().split('\n')[1]  # Skip header line
    n_residues = len(fasta_content)
    n_atoms = n_residues * 3  # Each residue has 3 atoms (N, CA, C)
    
    # Add extra particles for water beads
    n_water_beads = 8  # Number of water beads
    n_total_atoms = n_atoms + n_water_beads
    
    # Create a mapping between MARTINI atoms and UPSIDE atoms
    # For now, we'll map each MARTINI atom to the nearest UPSIDE atom
    martini_pos = initial_positions
    upside_pos = np.zeros((n_total_atoms, 3))
    
    # Create a simple linear chain for UPSIDE atoms
    for i in range(n_atoms):
        upside_pos[i] = [i * 3.8, 0, 0]  # 3.8Å spacing between atoms
    
    # Add water bead positions
    for i in range(n_water_beads):
        upside_pos[n_atoms + i] = martini_pos[i]
    
    print(f"Final positions shape: {upside_pos.shape}")
    print(f"Final positions:\n{upside_pos}")
    
    # Update the position array in the input group
    if hasattr(t.root.input, 'pos'):
        t.remove_node(t.root.input, 'pos')
    
    # Create position array with correct shape and attributes
    pos = np.zeros((n_total_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = upside_pos
    
    # Create the position array with proper attributes
    pos_array = t.create_array(t.root.input, 'pos', obj=pos)
    pos_array._v_attrs.arguments = np.array([b'pos'])
    pos_array._v_attrs.shape = pos.shape
    
    # Create arrays in the martini_potential group
    # All atoms are water beads (type 0)
    atoms = np.zeros(n_total_atoms, dtype=int)  # Changed to include all atoms
    t.create_array(martini_group, 'atom_indices', obj=atoms)
    t.create_array(martini_group, 'charges', obj=np.zeros(n_total_atoms))  # Water is neutral
    
    # Add MARTINI force field parameters
    # Create pairs array for MARTINI interactions
    # Create pairs between all atoms
    n_pairs = 0
    pairs_list = []
    for i in range(n_total_atoms):
        for j in range(i+1, n_total_atoms):
            pairs_list.append([i, j])
            n_pairs += 1
    
    pairs_array = np.array(pairs_list, dtype=int)
    print(f"Number of pairs: {n_pairs}")
    print(f"First few pairs:\n{pairs_array[:5]}")
    
    # Create coefficients array for MARTINI interactions
    coeff_array = np.zeros((n_pairs, 4))  # [epsilon, sigma, charge1, charge2]
    for i in range(n_pairs):
        # Use MARTINI water parameters for all interactions
        coeff_array[i] = [1.715293655572, 4.7, 0.0, 0.0]  # MARTINI water parameters (unit converted epsilon)
    
    print(f"Coefficients shape: {coeff_array.shape}")
    print(f"First few coefficients:\n{coeff_array[:5]}")
    
    # Add these arrays to the martini_potential group
    t.create_array(martini_group, 'pairs', obj=pairs_array)
    t.create_array(martini_group, 'coefficients', obj=coeff_array)
    
    # Remove all other potential groups except martini_potential
    groups_to_remove = [group for group in t.root.input.potential._v_groups 
                       if group not in ['martini_potential']]
    for group in groups_to_remove:
        t.remove_node(t.root.input.potential, group, recursive=True)

# Create wall-const-xyz.dat file
wall_file = "{}/wall-const-xyz.dat".format(input_dir)
with open(wall_file, 'w') as f:
    f.write("residue radius spring_const wall_type x0 y0 z0\n")
    f.write("0 50.0 100.0 1 0 0 0\n")  # Wall at origin with radius 50Å and force constant 100 kJ/mol/Å^2

# Add wall potential using advanced_config
kwargs = dict(
    fixed_wall = wall_file
)
config_stdout = ru.advanced_config("{}/martini.up".format(input_dir), **kwargs)
print("Advanced Config commandline options:")
print(config_stdout)

# Now copy the temporary file to the final location
import shutil
shutil.copy2(temp_file, "{}/martini.up".format(input_dir))

# Inspect the generated file
print("\nInspecting generated .up file:")
with tb.open_file("{}/martini.up".format(input_dir), 'r') as t:
    def print_structure(name, obj):
        if isinstance(obj, tb.Group):
            print(f"\nGroup: {name}")
            for key, value in obj._v_attrs.items():
                print(f"  Attribute: {key} = {value}")
        elif isinstance(obj, tb.Array):
            print(f"\nDataset: {name}")
            print(f"  Shape: {obj.shape}")
            print(f"  Type: {obj.dtype}")
            for key, value in obj._v_attrs.items():
                print(f"  Attribute: {key} = {value}")
    
    t.walk_nodes('/', print_structure)

print("\nDone generating .up file. Please inspect the output above to verify the structure.")

#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

randomseed = 1  # Fixed seed for reproducibility

upside_opts = (
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--seed {} "
)
upside_opts = upside_opts.format(duration, frame_interval, T, randomseed)

h5_file  = "{}/martini.run.up".format(run_dir)
log_file = "{}/martini.run.log".format(run_dir)

# Copy the input file to the run directory using PyTables
try:
    # Open source file
    with tb.open_file("{}/martini.up".format(input_dir), 'r') as src:
        # Create destination file
        with tb.open_file(h5_file, 'w') as dst:
            # Copy all groups and arrays
            src.copy_children(src.root, dst.root, recursive=True)
except Exception as e:
    print(f"Error copying HDF5 file: {e}")
    sys.exit(1)

print("Running...")

# Check if we're on a SLURM system
is_slurm = 'SLURM_JOB_ID' in os.environ

if is_slurm:
    # SLURM system - use srun
    cmd = "srun {}/obj/upside {} {} | tee {}".format(upside_path, upside_opts, h5_file, log_file)
else:
    # Non-SLURM system (like Mac) - run directly
    cmd = "{}/obj/upside {} {} | tee {}".format(upside_path, upside_opts, h5_file, log_file)

# Execute the command
try:
    sp.check_call(cmd, shell=True)
except sp.CalledProcessError as e:
    print(f"Error running UPSIDE: {e}")
    sys.exit(1)
except FileNotFoundError as e:
    print(f"Error: Could not find UPSIDE executable. Make sure UPSIDE_HOME is set correctly.")
    print(f"Current UPSIDE_HOME: {upside_path}")
    sys.exit(1)
