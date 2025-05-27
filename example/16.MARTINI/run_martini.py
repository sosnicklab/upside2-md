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
duration       = 10000  # Increased from 1 to 10000 steps
frame_interval = 50   # Increased from 1 to 100 steps
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

# Open the temporary file for modification
with tb.open_file(temp_file, 'r+') as t:
    # Create martini_potential group
    martini_group = t.create_group(t.root.input.potential, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'pair'
    martini_group._v_attrs.epsilon = 1.715293655572  # Original MARTINI epsilon (kJ/mol)
    martini_group._v_attrs.sigma = 4.7    # Original MARTINI sigma (Å)
    martini_group._v_attrs.lj_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.coul_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.dielectric = 15.0  # MARTINI water dielectric
    martini_group._v_attrs.n_types = 1  # Only one type of particle (water)
    martini_group._v_attrs.n_params = 4  # [epsilon, sigma, charge1, charge2]
    martini_group._v_attrs.cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.cache_buffer = 1.0  # Default cache buffer
    martini_group._v_attrs.initialized = True  # Mark as initialized
    
    # Create arrays in the martini_potential group
    # All atoms are water beads (type 0)
    n_atoms = len(initial_positions)
    atoms = np.arange(n_atoms, dtype=int)  # Use actual indices 0 to n_atoms-1
    t.create_array(martini_group, 'atom_indices', obj=atoms)
    t.create_array(martini_group, 'charges', obj=np.zeros(n_atoms))  # Water is neutral
    
    # Add MARTINI force field parameters
    # Create pairs array for MARTINI interactions
    # Create pairs between all atoms
    n_pairs = 0
    pairs_list = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])
            n_pairs += 1
    
    pairs_array = np.array(pairs_list, dtype=int)
    
    # Create coefficients array for MARTINI interactions
    coeff_array = np.zeros((n_pairs, 4))  # [epsilon, sigma, charge1, charge2]
    for i in range(n_pairs):
        # Use MARTINI water parameters for all interactions
        coeff_array[i] = [1.715293655572, 4.7, 0.0, 0.0]  # MARTINI water parameters
    
    # Add these arrays to the martini_potential group
    pairs_array = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_array._v_attrs.initialized = True
    coeff_array = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_array._v_attrs.initialized = True
    
    # Remove all other potential groups except martini_potential
    groups_to_remove = [group for group in t.root.input.potential._v_groups 
                       if group not in ['martini_potential']]
    for group in groups_to_remove:
        t.remove_node(t.root.input.potential, group, recursive=True)
    
    # Update the position array in the input group
    if hasattr(t.root.input, 'pos'):
        t.remove_node(t.root.input, 'pos')
    
    # Create position array with correct shape and attributes
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')  # Shape: (n_atoms, 3, 1)
    pos[:,:,0] = initial_positions  # Set initial positions in first frame
    
    # Create the position array with proper attributes
    pos_array = t.create_array(t.root.input, 'pos', obj=pos)
    pos_array._v_attrs.arguments = np.array([b'pos'])
    pos_array._v_attrs.shape = pos.shape
    pos_array._v_attrs.n_atoms = n_atoms
    pos_array._v_attrs.n_frames = 1
    pos_array._v_attrs.dim = 3
    pos_array._v_attrs.initialized = True  # Mark as initialized
    
    # Add initial velocities
    vel = np.zeros((n_atoms, 3), dtype='f4')
    # Add small random velocities to help system equilibrate
    vel = np.random.normal(0, 0.001, (n_atoms, 3))  # Very small initial velocities
    vel_array = t.create_array(t.root.input, 'vel', obj=vel)
    vel_array._v_attrs.arguments = np.array([b'vel'])
    vel_array._v_attrs.shape = vel.shape
    vel_array._v_attrs.n_atoms = n_atoms
    vel_array._v_attrs.dim = 3
    vel_array._v_attrs.initialized = True  # Mark as initialized

    # Add mass array
    mass = np.ones(n_atoms, dtype='f4') * 1.0  # Set mass to 1.0 for all particles
    mass_array = t.create_array(t.root.input, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True

    # Add type array
    type_array = t.create_array(t.root.input, 'type', obj=atoms)
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atoms.shape
    type_array._v_attrs.n_atoms = n_atoms
    type_array._v_attrs.initialized = True

# Create wall-const-xyz.dat file
wall_file = "{}/wall-const-xyz.dat".format(input_dir)
with open(wall_file, 'w') as f:
    f.write("residue radius spring_const wall_type x0 y0 z0\n")
    f.write("0 50.0 10.0 1 0 0 0\n")  # Restored original wall potential

# Add wall potential using advanced_config
kwargs = dict(
    fixed_wall = wall_file
)
config_stdout = ru.advanced_config("{}/martini.up".format(input_dir), **kwargs)

# Now copy the temporary file to the final location
import shutil
shutil.copy2(temp_file, "{}/martini.up".format(input_dir))

# Comment out or remove debug/inspection print statements
# print(f"Initial positions shape: {initial_positions.shape}")
# print(f"Initial positions:\n{initial_positions}")
# print(f"Number of pairs: {n_pairs}")
# print(f"First few pairs:\n{pairs_array[:5]}")
# print(f"Coefficients shape: {coeff_array.shape}")
# print(f"First few coefficients:\n{coeff_array[:5]}")
# print("Advanced Config commandline options:")
# print(config_stdout)
# print("\nInspecting generated .up file:")
# with tb.open_file("{}/martini.up".format(input_dir), 'r') as t:
#     def print_structure(name, obj):
#         if isinstance(obj, tb.Group):
#             print(f"\nGroup: {name}")
#             for key, value in obj._v_attrs.items():
#                 print(f"  Attribute: {key} = {value}")
#         elif isinstance(obj, tb.Array):
#             print(f"\nDataset: {name}")
#             print(f"  Shape: {obj.shape}")
#             print(f"  Type: {obj.dtype}")
#             for key, value in obj._v_attrs.items():
#                 print(f"  Attribute: {key} = {value}")
#     t.walk_nodes('/', print_structure)
# print("\nDone generating .up file. Please inspect the output above to verify the structure.")
# print("\nDetailed inspection of input .up file:")
# try:
#     with tb.open_file(h5_file, 'r') as t:
#         def print_structure(name, obj):
#             if isinstance(obj, tb.Group):
#                 print(f"\nGroup: {name}")
#                 print("  Attributes:")
#                 for key, value in obj._v_attrs.items():
#                     print(f"    {key} = {value}")
#                 print("  Children:")
#                 for child in obj._v_children:
#                     print(f"    {child}")
#             elif isinstance(obj, tb.Array):
#                 print(f"\nDataset: {name}")
#                 print(f"  Shape: {obj.shape}")
#                 print(f"  Type: {obj.dtype}")
#                 print("  Attributes:")
#                 for key, value in obj._v_attrs.items():
#                     print(f"    {key} = {value}")
#                 print("  Data (first few elements):")
#                 print(f"    {obj[:5]}")
#         t.walk_nodes('/', print_structure)
# except Exception as e:
#     print(f"Error inspecting HDF5 file: {e}")
#     sys.exit(1)
# print(f"Command: {upside_path}/obj/upside {upside_opts}")
# print("\nAnalyzing trajectory...")

#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

randomseed = 1  # Fixed seed for reproducibility

# Define file paths first
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

# Put the file path first, then the options
upside_opts = (
    "{} "  # File path first
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--seed {}"
)
upside_opts = upside_opts.format(h5_file, duration, frame_interval, T, randomseed)

print("\nRunning...")
print(f"Command: {upside_path}/obj/upside {upside_opts}")

# Check if we're on a SLURM system
is_slurm = 'SLURM_JOB_ID' in os.environ

if is_slurm:
    # SLURM system - use srun
    cmd = "srun {}/obj/upside {} | tee {}".format(upside_path, upside_opts, log_file)
else:
    # Non-SLURM system (like Mac) - run directly
    cmd = "{}/obj/upside {} | tee {}".format(upside_path, upside_opts, log_file)

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

# Analyze the trajectory
print("\nAnalyzing trajectory...")
try:
    with tb.open_file(h5_file, 'r') as t:
        # Get the number of frames
        if not hasattr(t.root, 'output') or not hasattr(t.root.output, 'pos'):
            print("No output data found in trajectory file")
            sys.exit(1)
        
        n_frames = t.root.output.pos.shape[0]
        print(f"Total number of frames: {n_frames}")
        
        # Get the last 5 frames (or all frames if less than 5)
        last_frames = min(5, n_frames)
        print(f"\nAnalyzing last {last_frames} frames:")
        
        for frame in range(n_frames - last_frames, n_frames):
            pos = t.root.output.pos[frame, 0, :, :]  # shape: (n_atoms, 3)
            print(f"\nFrame {frame}:")
            print(f"Shape: {pos.shape}")
            print("First few positions:")
            print(pos[:5])  # Print first 5 positions
            print("Last few positions:")
            print(pos[-5:])  # Print last 5 positions
            
            # Calculate distances between consecutive particles
            if len(pos) > 1:
                distances = np.linalg.norm(pos[1:] - pos[:-1], axis=1)
                print(f"Average distance between consecutive particles: {np.mean(distances):.2f} Å")
                print(f"Min distance: {np.min(distances):.2f} Å")
                print(f"Max distance: {np.max(distances):.2f} Å")
            
            # Print potential energy if available
            if hasattr(t.root.output, 'potential'):
                try:
                    pot = t.root.output.potential[frame, 0]
                    print(f"Potential energy: {pot}")
                except IndexError:
                    print("Potential energy not available for this frame")
            
            # Print kinetic energy if available
            if hasattr(t.root.output, 'kinetic'):
                try:
                    kin = t.root.output.kinetic[frame, 0]
                    print(f"Kinetic energy: {kin}")
                except IndexError:
                    print("Kinetic energy not available for this frame")
            
            # Print total energy if available
            if hasattr(t.root.output, 'total'):
                try:
                    tot = t.root.output.total[frame, 0]
                    print(f"Total energy: {tot}")
                except IndexError:
                    print("Total energy not available for this frame")
except Exception as e:
    print(f"Error analyzing trajectory: {e}")
    sys.exit(1)
