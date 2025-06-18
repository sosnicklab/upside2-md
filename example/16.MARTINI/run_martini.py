import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

sim_id         = 'martini_test'
base_dir       = './'

# MARTINI parameters
martini_epsilon = 1.715293655572  # Test with overflow protection
martini_sigma   = 4.7  # MARTINI water sigma (Angstroms)

print(f"Testing epsilon = {martini_epsilon:.3f} UPSIDE units = {martini_epsilon * 2.332:.3f} kJ/mol with overflow protection")

# Simulation parameters
T              = 0.8  # Temperature (UPSIDE units)
duration       = 1000  # Total simulation steps  
frame_interval = 50   # Output every N steps
dt             = 0.001  # Time step

# Wall box size (Angstroms) - should contain centered particles  
wall_box_size = 25

#----------------------------------------------------------------------
## Setup directories
#----------------------------------------------------------------------

input_dir = "{}/inputs".format(base_dir)
output_dir = "{}/outputs".format(base_dir)
run_dir = "{}/{}".format(output_dir, sim_id)

make_dirs = [input_dir, output_dir, run_dir]
for direc in make_dirs:
    if not os.path.exists(direc):
        os.makedirs(direc)

#----------------------------------------------------------------------
## Create input structure
#----------------------------------------------------------------------

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
n_atoms = len(initial_positions)

print(f"Loaded {n_atoms} atoms from PDB")

# Check initial particle range
pos_min = np.min(initial_positions, axis=0)
pos_max = np.max(initial_positions, axis=0)
pos_range = pos_max - pos_min
print(f"Initial position range: X=[{pos_min[0]:.1f}, {pos_max[0]:.1f}], Y=[{pos_min[1]:.1f}, {pos_max[1]:.1f}], Z=[{pos_min[2]:.1f}, {pos_max[2]:.1f}]")
print(f"Position range: X={pos_range[0]:.1f}, Y={pos_range[1]:.1f}, Z={pos_range[2]:.1f} Angstroms")

# Center the particles around origin
center = (pos_max + pos_min) / 2
initial_positions -= center
print(f"Centered particles around origin. New center: {np.mean(initial_positions, axis=0)}")

# Check if particles fit within wall boundaries
new_min = np.min(initial_positions, axis=0)
new_max = np.max(initial_positions, axis=0)
max_coord = np.max(np.abs([new_min, new_max]))
print(f"After centering, max coordinate magnitude: {max_coord:.1f} Angstroms")

if max_coord > wall_box_size:
    print(f"WARNING: Particles extend beyond wall boundaries (±{wall_box_size} Angstroms)")
    print(f"Consider increasing wall_box_size to at least {max_coord + 2:.1f} Angstroms")
else:
    print(f"Particles fit within wall boundaries (±{wall_box_size} Angstroms)")

# Create HDF5 input file
input_file = "{}/test.up".format(input_dir)

with tb.open_file(input_file, 'w') as t:
    input_grp = t.create_group(t.root, 'input')
    
    # Position array
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = initial_positions
    
    pos_array = t.create_array(input_grp, 'pos', obj=pos)
    pos_array._v_attrs.arguments = np.array([b'pos'])
    pos_array._v_attrs.shape = pos.shape
    pos_array._v_attrs.n_atoms = n_atoms
    pos_array._v_attrs.n_frames = 1
    pos_array._v_attrs.dim = 3
    pos_array._v_attrs.initialized = True
    
    # Velocity array
    velocity = np.zeros((n_atoms, 3), dtype='f4')
    vel_array = t.create_array(input_grp, 'vel', obj=velocity)
    vel_array._v_attrs.arguments = np.array([b'vel'])
    vel_array._v_attrs.shape = velocity.shape
    vel_array._v_attrs.n_atoms = n_atoms
    vel_array._v_attrs.dim = 3
    vel_array._v_attrs.initialized = True

    # Mass array
    mass = np.ones(n_atoms, dtype='f4') * 1.0
    mass_array = t.create_array(input_grp, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True

    # Type array
    atoms = np.arange(n_atoms, dtype=int)
    type_array = t.create_array(input_grp, 'type', obj=atoms)
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atoms.shape
    type_array._v_attrs.n_atoms = n_atoms
    type_array._v_attrs.initialized = True

    # Potential group
    potential_grp = t.create_group(input_grp, 'potential')
    
    # Create MARTINI potential group
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = martini_epsilon
    martini_group._v_attrs.sigma = martini_sigma
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = 15.0
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    
    # Create wall potential group
    wall_group = t.create_group(potential_grp, 'periodic_boundary_potential')
    wall_group._v_attrs.arguments = np.array([b'pos'])
    wall_group._v_attrs.wall_xlo = -wall_box_size
    wall_group._v_attrs.wall_xhi = wall_box_size
    wall_group._v_attrs.wall_ylo = -wall_box_size
    wall_group._v_attrs.wall_yhi = wall_box_size
    wall_group._v_attrs.wall_zlo = -wall_box_size
    wall_group._v_attrs.wall_zhi = wall_box_size
    wall_group._v_attrs.initialized = True
    
    # Add atom indices and charges for MARTINI
    t.create_array(martini_group, 'atom_indices', obj=atoms)
    t.create_array(martini_group, 'charges', obj=np.zeros(n_atoms))
    
    # Create pairs and coefficients
    pairs_list = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])

    pairs_array = np.array(pairs_list, dtype=int)
    n_pairs = len(pairs_list)

    coeff_array = np.zeros((n_pairs, 4))
    for i in range(n_pairs):
        coeff_array[i] = [martini_epsilon, martini_sigma, 0.0, 0.0]
    
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True

print(f"Created input with {n_pairs} LJ pairs")

#----------------------------------------------------------------------
## Run MD Simulation 
#----------------------------------------------------------------------

h5_file = "{}/test.run.up".format(run_dir)
log_file = "{}/test.run.log".format(run_dir)

# Copy input to run directory
with tb.open_file(input_file, 'r') as src:
    with tb.open_file(h5_file, 'w') as dst:
        src.copy_children(src.root, dst.root, recursive=True)

print(f"\nRunning MD simulation:")
print(f"  Temperature: {T:.2f} UPSIDE units = {T * 350.59:.1f} K")
print(f"  Duration: {duration} steps")
print(f"  Time step: {dt}")
print(f"  Frame interval: {frame_interval}")

# Run simulation
upside_opts = (
    "{} "
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--time-step {} "
    "--seed {}"
)
upside_opts = upside_opts.format(h5_file, duration, frame_interval, T, dt, 12345)

cmd = "{}/obj/upside {}".format(upside_path, upside_opts)

print(f"Command: {cmd}")
result = sp.run(cmd, shell=True)

if result.returncode != 0:
    print("Simulation failed!")
    print("Return code:", result.returncode)
else:
    print("Simulation completed successfully!")

#----------------------------------------------------------------------
## Convert trajectory to VTF format for visualization
#----------------------------------------------------------------------

print("\nConverting trajectory to VTF format...")

vtf_script = "extract_martini_vtf.py"
input_traj = h5_file
output_vtf = "{}/martini_trajectory.vtf".format(run_dir)

# Run VTF conversion
convert_cmd = "python {} {} {}".format(vtf_script, input_traj, output_vtf)
print(f"VTF conversion command: {convert_cmd}")

vtf_result = sp.run(convert_cmd, shell=True)

if vtf_result.returncode != 0:
    print("VTF conversion failed!")
    print("Return code:", vtf_result.returncode)
else:
    print(f"VTF trajectory saved to: {output_vtf}")

print("Done!")
