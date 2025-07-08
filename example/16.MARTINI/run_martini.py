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

# Box dimensions (Angstroms) - now configurable
x_len = 50.0  # Box length in X direction (Angstroms)
y_len = 50.0  # Box length in Y direction (Angstroms)
z_len = 50.0  # Box length in Z direction (Angstroms)

# Input PDB file
input_pdb_file = 'input.pdb'

# MARTINI parameters - will be calculated from table below
martini_sigma   = 4.7  # MARTINI water sigma (Angstroms)
dielectric_constant = 15.0  # Default MARTINI dielectric, user can change

# Minimization parameters
minimization_steps = 1000
minimization_T = 0.01  # Very low temperature for minimization
minimization_dt = 0.01
minimization_frame_interval = 100

# Conjugate gradient minimizer parameters
minimizer_max_iterations = 1000
minimizer_energy_tolerance = 1e-6
minimizer_force_tolerance = 1e-6
minimizer_step_size = 0.01  # Lowered from 0.1 for stability
minimizer_verbose = True

# Simulation parameters
T              = 0.8  # Temperature (UPSIDE units)
duration       = 1000  # Total simulation steps
frame_interval = 50   # Output every N steps
dt             = 0.005  # Reduced time step for stability (was 0.01)
thermostat_timescale = 0.0  # Thermostat timescale (10x from default 5.0 for increased Langevin factor)

martini_table = {'Qda': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qd': {'Qda': 5.6, 'Qd': 5.0, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 4.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qa': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Q0': {'Qda': 4.5, 'Qd': 4.5, 'Qa': 4.5, 'Q0': 3.5, 'P5': 5.0, 'P4': 5.6, 'P3': 5.0, 'P2': 4.5, 'P1': 4.0, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'P5': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.6, 'P1': 5.6, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P4': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.6, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P3': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P2': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P1': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'Nda': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Nd': {'Qda': 5.0, 'Qd': 4.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.0, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Na': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 4.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.0, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'N0': {'Qda': 3.5, 'Qd': 3.5, 'Qa': 3.5, 'Q0': 3.5, 'P5': 3.5, 'P4': 3.5, 'P3': 3.5, 'P2': 4.0, 'P1': 4.0, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'C5': {'Qda': 3.1, 'Qd': 3.1, 'Qa': 3.1, 'Q0': 3.1, 'P5': 3.1, 'P4': 3.1, 'P3': 3.5, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C4': {'Qda': 2.7, 'Qd': 2.7, 'Qa': 2.7, 'Q0': 2.7, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.1, 'Nd': 3.1, 'Na': 3.1, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C3': {'Qda': 2.3, 'Qd': 2.3, 'Qa': 2.3, 'Q0': 2.3, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.1, 'P1': 3.5, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C2': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.3, 'P4': 2.3, 'P3': 2.7, 'P2': 2.7, 'P1': 3.1, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.1, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C1': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.0, 'P4': 2.0, 'P3': 2.3, 'P2': 2.3, 'P1': 2.7, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 2.7, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}}

# Calculate MARTINI epsilon from table for P4 (water) interactions
martini_epsilon_table = martini_table['P4']['P4']  # P4-P4 interaction strength (kJ/mol)
martini_epsilon = martini_epsilon_table / 2.914952774272  # Convert to UPSIDE units

# Softened LJ parameters for minimization
soft_epsilon = martini_epsilon * 0.1  # 10x softer
soft_sigma = martini_sigma * 1.1      # 10% larger

print(f"Box dimensions: X={x_len:.1f}, Y={y_len:.1f}, Z={z_len:.1f} Angstroms")
print(f"P4-P4 epsilon from MARTINI table: {martini_epsilon_table:.1f} kJ/mol")
print(f"Converted to UPSIDE units: {martini_epsilon:.6f}")
print(f"Equivalent in kJ/mol: {martini_epsilon * 2.332:.3f} kJ/mol")
print(f"Thermostat timescale: {thermostat_timescale:.1f} (10x from default 5.0)")

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
# All particles are P4 (water) with zero charge
pdb_to_martini = {'W': 'P4', 'NA': 'P4', 'CL': 'P4'}  # All mapped to P4
martini_charges = {'P4': 0.0}  # All particles have zero charge

if not os.path.exists(input_pdb_file):
    print(f"Error: Input PDB file '{input_pdb_file}' not found!")
    sys.exit(1)

with open(input_pdb_file, 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            pdb_name = line[12:16].strip().upper()
            # All particles are P4 (water)
            mtype = 'P4'
            atom_types.append(mtype)
            charges.append(0.0)  # All particles have zero charge

if not initial_positions:
    print("Error: No atoms found in PDB file!")
    sys.exit(1)

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

# Center the particles properly for periodic boundary conditions
# Move particles to range [0, box_dimension] instead of centering around origin
center = (pos_max + pos_min) / 2
initial_positions -= center  # First center around origin

# Now shift to box center so coordinates are in [0, box_dimension] range
box_center = np.array([x_len/2, y_len/2, z_len/2])
initial_positions += box_center  # Shift to box center

print(f"Positioned particles for periodic boundaries. New center: {np.mean(initial_positions, axis=0)}")

# Check if particles fit within box boundaries [0, box_dimension]
new_min = np.min(initial_positions, axis=0)
new_max = np.max(initial_positions, axis=0)
box_dims = np.array([x_len, y_len, z_len])
print(f"After repositioning: min=({new_min[0]:.1f}, {new_min[1]:.1f}, {new_min[2]:.1f}), max=({new_max[0]:.1f}, {new_max[1]:.1f}, {new_max[2]:.1f})")

if np.any(new_min < 0) or np.any(new_max > box_dims):
    print(f"WARNING: Particles extend beyond box boundaries [0, {box_dims}]")
    print(f"Consider increasing box dimensions")
else:
    print(f"Particles properly positioned within box boundaries [0, {box_dims}]")

# --- MINIMIZATION: Create minimization input with softened parameters ---
print("\nCreating minimization input with softened LJ parameters...")
print(f"Soft epsilon: {soft_epsilon:.6f} (10x softer)")
print(f"Soft sigma: {soft_sigma:.2f} Angstroms (10% larger)")

min_input_file = f"{input_dir}/minimize.up"
with tb.open_file(min_input_file, 'w') as t:
    input_grp = t.create_group(t.root, 'input')
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = initial_positions  # Use centered positions from PDB
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
    mass = np.ones(n_atoms, dtype='f4') * 6.0  # Changed from 72.0 to 6.0
    mass_array = t.create_array(input_grp, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True
    type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atom_types.shape
    type_array._v_attrs.n_atoms = n_atoms
    type_array._v_attrs.initialized = True
    charge_array = t.create_array(input_grp, 'charges', obj=charges)
    charge_array._v_attrs.arguments = np.array([b'charges'])
    charge_array._v_attrs.shape = charges.shape
    charge_array._v_attrs.n_atoms = n_atoms
    charge_array._v_attrs.initialized = True
    potential_grp = t.create_group(input_grp, 'potential')
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = soft_epsilon  # Use softened epsilon for minimization
    martini_group._v_attrs.sigma = soft_sigma      # Use larger sigma for minimization
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 1  # Enable force capping to prevent blowup
    martini_group._v_attrs.debug_mode = 0  # Disable debug mode
    martini_group._v_attrs.x_len = x_len
    martini_group._v_attrs.y_len = y_len
    martini_group._v_attrs.z_len = z_len
    wall_group = t.create_group(potential_grp, 'periodic_boundary_potential')
    wall_group._v_attrs.arguments = np.array([b'pos'])
    wall_group._v_attrs.x_len = x_len
    wall_group._v_attrs.y_len = y_len
    wall_group._v_attrs.z_len = z_len
    wall_group._v_attrs.initialized = True
    t.create_array(martini_group, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_group, 'charges', obj=charges)
    pairs_list = []
    coeff_array = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            epsilon = epsilon_sim * 0.1  # Apply softening factor for minimization
            sigma_val = soft_sigma  # Use softened sigma
            q1 = charges[i]
            q2 = charges[j]
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True

print(f"Created minimization input with {len(pairs_list)} softened LJ pairs")

# --- RUN MINIMIZATION ---
print(f"\nRunning minimization:")
print(f"  Steps: {minimization_steps}")
print(f"  Temperature: {minimization_T}")
print(f"  Time step: {minimization_dt}")
print(f"  Thermostat timescale: {thermostat_timescale:.1f}")

min_h5_file = "{}/minimize.run.up".format(run_dir)

# Copy minimization input to run directory
with tb.open_file(min_input_file, 'r') as src:
    with tb.open_file(min_h5_file, 'w') as dst:
        src.copy_children(src.root, dst.root, recursive=True)

# Run minimization simulation
min_opts = (
    "{} "
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--time-step {} "
    "--thermostat-timescale {} "
    "--seed {}"
)
min_opts = min_opts.format(min_h5_file, minimization_steps, minimization_frame_interval, minimization_T, minimization_dt, thermostat_timescale, 12345)

cmd = "{}/obj/upside {}".format(upside_path, min_opts)

print(f"Minimization command: {cmd}")
result = sp.run(cmd, shell=True)

if result.returncode != 0:
    print("Minimization failed!")
    print("Return code:", result.returncode)
    sys.exit(1)
else:
    print("Minimization completed successfully!")

# --- EXTRACT MINIMIZED POSITIONS ---
print("Extracting minimized positions...")
with tb.open_file(min_h5_file, 'r') as t:
    final_pos = t.root.output.pos[-1]  # Get final frame positions
    minimized_positions = final_pos.copy()
    
print(f"Extracted minimized positions with shape: {minimized_positions.shape}")

# Check position drift during minimization
pos_drift = np.linalg.norm(minimized_positions - initial_positions, axis=1)
max_drift = np.max(pos_drift)
mean_drift = np.mean(pos_drift)
print(f"Position drift during minimization: max={max_drift:.3f}, mean={mean_drift:.3f} Angstroms")

# --- CREATE MD INPUT WITH MINIMIZED POSITIONS ---
print("Creating MD input with minimized positions and full MARTINI parameters...")
input_file = f"{input_dir}/test.up"
with tb.open_file(input_file, 'w') as t:
    input_grp = t.create_group(t.root, 'input')
    pos = np.zeros((n_atoms, 3, 1), dtype='f4')
    pos[:,:,0] = minimized_positions  # Use minimized positions for MD
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
    mass = np.ones(n_atoms, dtype='f4') * 6.0  # Changed from 72.0 to 6.0
    mass_array = t.create_array(input_grp, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True
    type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
    type_array._v_attrs.arguments = np.array([b'type'])
    type_array._v_attrs.shape = atom_types.shape
    type_array._v_attrs.n_atoms = n_atoms
    type_array._v_attrs.initialized = True
    charge_array = t.create_array(input_grp, 'charges', obj=charges)
    charge_array._v_attrs.arguments = np.array([b'charges'])
    charge_array._v_attrs.shape = charges.shape
    charge_array._v_attrs.n_atoms = n_atoms
    charge_array._v_attrs.initialized = True
    potential_grp = t.create_group(input_grp, 'potential')
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = martini_epsilon  # Use full MARTINI epsilon for MD
    martini_group._v_attrs.sigma = martini_sigma      # Use standard sigma for MD
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 1  # Enable force capping to prevent blowup
    martini_group._v_attrs.debug_mode = 0  # Disable debug mode
    martini_group._v_attrs.x_len = x_len
    martini_group._v_attrs.y_len = y_len
    martini_group._v_attrs.z_len = z_len
    wall_group = t.create_group(potential_grp, 'periodic_boundary_potential')
    wall_group._v_attrs.arguments = np.array([b'pos'])
    wall_group._v_attrs.x_len = x_len
    wall_group._v_attrs.y_len = y_len
    wall_group._v_attrs.z_len = z_len
    wall_group._v_attrs.initialized = True
    t.create_array(martini_group, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_group, 'charges', obj=charges)
    pairs_list = []
    coeff_array = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            epsilon = epsilon_sim  # Use full MARTINI epsilon for MD
            sigma_val = martini_sigma  # Use standard sigma for MD
            q1 = charges[i]
            q2 = charges[j]
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True
print(f"Created MD input with {len(pairs_list)} LJ pairs using minimized positions")

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
print(f"  Thermostat timescale: {thermostat_timescale:.1f}")

# Run simulation
upside_opts = (
    "{} "
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--time-step {} "
    "--thermostat-timescale {} "
    "--seed {}"
)
upside_opts = upside_opts.format(h5_file, duration, frame_interval, T, dt, thermostat_timescale, 12345)

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