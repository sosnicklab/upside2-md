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
#martini_sigma   = 0.47  # MARTINI water sigma (Angstroms)

dielectric_constant = 15.0  # Default MARTINI dielectric, user can change

# Softened LJ parameters for minimization
soft_epsilon = martini_epsilon * 0.1  # 10x softer
soft_sigma = martini_sigma * 1.1      # 10% larger

# Minimization parameters
minimization_steps = 1000
minimization_T = 0.01  # Very low temperature for minimization
minimization_dt = 0.01
minimization_frame_interval = 100

print(f"Testing epsilon = {martini_epsilon:.3f} UPSIDE units = {martini_epsilon * 2.332:.3f} kJ/mol with overflow protection")

# Simulation parameters
T              = 0.8  # Temperature (UPSIDE units)
duration       = 1000  # Total simulation steps  
frame_interval = 50   # Output every N steps
dt             = 0.01  # Time step

# Wall box size (Angstroms) - should contain centered particles  
wall_box_size = 23


martini_table = {'Qda': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qd': {'Qda': 5.6, 'Qd': 5.0, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 4.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qa': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Q0': {'Qda': 4.5, 'Qd': 4.5, 'Qa': 4.5, 'Q0': 3.5, 'P5': 5.0, 'P4': 5.6, 'P3': 5.0, 'P2': 4.5, 'P1': 4.0, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'P5': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.6, 'P1': 5.6, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P4': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.6, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P3': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P2': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P1': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'Nda': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Nd': {'Qda': 5.0, 'Qd': 4.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.0, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Na': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 4.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.0, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'N0': {'Qda': 3.5, 'Qd': 3.5, 'Qa': 3.5, 'Q0': 3.5, 'P5': 3.5, 'P4': 3.5, 'P3': 3.5, 'P2': 4.0, 'P1': 4.0, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'C5': {'Qda': 3.1, 'Qd': 3.1, 'Qa': 3.1, 'Q0': 3.1, 'P5': 3.1, 'P4': 3.1, 'P3': 3.5, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C4': {'Qda': 2.7, 'Qd': 2.7, 'Qa': 2.7, 'Q0': 2.7, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.1, 'Nd': 3.1, 'Na': 3.1, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C3': {'Qda': 2.3, 'Qd': 2.3, 'Qa': 2.3, 'Q0': 2.3, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.1, 'P1': 3.5, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C2': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.3, 'P4': 2.3, 'P3': 2.7, 'P2': 2.7, 'P1': 3.1, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.1, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C1': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.0, 'P4': 2.0, 'P3': 2.3, 'P2': 2.3, 'P1': 2.7, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 2.7, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}}

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

# Read initial positions and atom types from PDB
initial_positions = []
atom_types = []
charges = []
# Mapping: W -> P4, NA -> Qd, CL -> Qa
pdb_to_martini = {'W': 'P4', 'NA': 'Qd', 'CL': 'Qa'}
martini_charges = {'P4': 0.0, 'Qd': 1.0, 'Qa': -1.0}
with open("input.pdb", 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            pdb_name = line[12:16].strip().upper()
            # Assign MARTINI type
            if pdb_name in pdb_to_martini:
                mtype = pdb_to_martini[pdb_name]
            else:
                mtype = 'P4'  # Default to P4 (water)
            atom_types.append(mtype)
            charges.append(martini_charges.get(mtype, 0.0))
initial_positions = np.array(initial_positions)
n_atoms = len(initial_positions)
atom_types = np.array(atom_types)
charges = np.array(charges)

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

# Create HDF5 input file for minimization
min_input_file = f"{input_dir}/minimize.up"

with tb.open_file(min_input_file, 'w') as t:
    input_grp = t.create_group(t.root, 'input')
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = initial_positions
    pos_array = t.create_array(input_grp, 'pos', obj=pos)
    pos_array._v_attrs.arguments = np.array([b'pos'])
    pos_array._v_attrs.shape = pos.shape
    pos_array._v_attrs.n_atoms = n_atoms
    pos_array._v_attrs.n_frames = 1
    pos_array._v_attrs.dim = 3
    pos_array._v_attrs.initialized = True
    velocity = np.zeros((n_atoms, 3), dtype='f4')
    vel_array = t.create_array(input_grp, 'vel', obj=velocity)
    vel_array._v_attrs.arguments = np.array([b'vel'])
    vel_array._v_attrs.shape = velocity.shape
    vel_array._v_attrs.n_atoms = n_atoms
    vel_array._v_attrs.dim = 3
    vel_array._v_attrs.initialized = True
    mass = np.ones(n_atoms, dtype='f4') * 1.0
    mass_array = t.create_array(input_grp, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True
    # Store MARTINI type as string array
    type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atom_types.shape
    type_array._v_attrs.n_atoms = n_atoms
    type_array._v_attrs.initialized = True
    # Store charges
    charge_array = t.create_array(input_grp, 'charges', obj=charges)
    charge_array._v_attrs.arguments = np.array([b'charges'])
    charge_array._v_attrs.shape = charges.shape
    charge_array._v_attrs.n_atoms = n_atoms
    charge_array._v_attrs.initialized = True
    potential_grp = t.create_group(input_grp, 'potential')
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = soft_epsilon
    martini_group._v_attrs.sigma = soft_sigma
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 1  # Enable force capping for minimization
    wall_group = t.create_group(potential_grp, 'periodic_boundary_potential')
    wall_group._v_attrs.arguments = np.array([b'pos'])
    wall_group._v_attrs.wall_xlo = -wall_box_size
    wall_group._v_attrs.wall_xhi = wall_box_size
    wall_group._v_attrs.wall_ylo = -wall_box_size
    wall_group._v_attrs.wall_yhi = wall_box_size
    wall_group._v_attrs.wall_zlo = -wall_box_size
    wall_group._v_attrs.wall_zhi = wall_box_size
    wall_group._v_attrs.initialized = True
    t.create_array(martini_group, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_group, 'charges', obj=charges)
    # Create pairs and coefficients (LJ and Coulomb)
    pairs_list = []
    coeff_array = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon from table and convert to simulation units
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            # Use softened params for minimization
            if 'minimize' in min_input_file:
                epsilon = epsilon_sim * 0.1
                sigma_val = 4.7 * 1.1
            else:
                epsilon = epsilon_sim
                sigma_val = 4.7
            q1 = charges[i]
            q2 = charges[j]
            # [epsilon, sigma, q1, q2] for both LJ and Coulomb
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True

print(f"Created minimization input with {len(pairs_list)} LJ pairs")

# Run minimization
min_run_dir = f"{output_dir}/minimize"
if not os.path.exists(min_run_dir):
    os.makedirs(min_run_dir)
min_h5_file = f"{min_run_dir}/minimize.run.up"
with tb.open_file(min_input_file, 'r') as src:
    with tb.open_file(min_h5_file, 'w') as dst:
        src.copy_children(src.root, dst.root, recursive=True)
min_cmd = f"{upside_path}/obj/upside {min_h5_file} --duration {minimization_steps} --frame-interval {minimization_frame_interval} --temperature {minimization_T} --time-step {minimization_dt} --seed 12345"
print(f"\nRunning minimization:\n  Command: {min_cmd}")
min_result = sp.run(min_cmd, shell=True)
if min_result.returncode != 0:
    print("Minimization failed!")
    sys.exit(1)
else:
    print("Minimization completed successfully!")

# Extract minimized positions from minimization output
with tb.open_file(min_h5_file, 'r') as t:
    min_pos = t.root.input.pos[:,:,-1]  # Last frame

# Create HDF5 input file for production MD (overwrite initial_positions with min_pos)
input_file = "{}/test.up".format(input_dir)
with tb.open_file(input_file, 'w') as t:
    input_grp = t.create_group(t.root, 'input')
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = min_pos  # Use minimized positions
    
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
    type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atom_types.shape
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
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 0  # Disable force capping for normal MD
    
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
    t.create_array(martini_group, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_group, 'charges', obj=charges)
    
    # Create pairs and coefficients (LJ and Coulomb)
    pairs_list = []
    coeff_array = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon from table and convert to simulation units
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            # Use softened params for minimization
            if 'minimize' in min_input_file:
                epsilon = epsilon_sim * 0.1
                sigma_val = 4.7 * 1.1
            else:
                epsilon = epsilon_sim
                sigma_val = 4.7
            q1 = charges[i]
            q2 = charges[j]
            # [epsilon, sigma, q1, q2] for both LJ and Coulomb
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True

print(f"Created input with {len(pairs_list)} LJ pairs")

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
