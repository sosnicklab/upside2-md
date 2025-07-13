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

# Input PDB file
input_pdb_file = 'input.pdb'

# Box dimensions (Angstroms) - will be read from input.pdb CRYST1 record
x_len = None  # Will be set from CRYST1 record
y_len = None  # Will be set from CRYST1 record  
z_len = None  # Will be set from CRYST1 record

# MARTINI parameters - will be calculated from table below
martini_sigma   = 4.7  # MARTINI water sigma (Angstroms)
dielectric_constant = 15.0  # MARTINI dielectric constant (standard value)

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

# Simulation parameters - optimized for stability
T              = 0.86  # Temperature (UPSIDE units) - proper MARTINI temperature ~300K  
duration       = 100  # Total simulation steps (production run)
frame_interval = 20   # Output every N steps (more frames for better trajectory)
dt             = 0.001  # Time step - reduced for stability with lipids
thermostat_timescale = 5.0  # Thermostat timescale (default Langevin damping)

martini_table = {'Qda': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qd': {'Qda': 5.6, 'Qd': 5.0, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 4.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qa': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Q0': {'Qda': 4.5, 'Qd': 4.5, 'Qa': 4.5, 'Q0': 3.5, 'P5': 5.0, 'P4': 5.6, 'P3': 5.0, 'P2': 4.5, 'P1': 4.0, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'P5': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.6, 'P1': 5.6, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P4': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.6, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P3': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P2': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P1': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'Nda': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Nd': {'Qda': 5.0, 'Qd': 4.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.0, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Na': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 4.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.0, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'N0': {'Qda': 3.5, 'Qd': 3.5, 'Qa': 3.5, 'Q0': 3.5, 'P5': 3.5, 'P4': 3.5, 'P3': 3.5, 'P2': 4.0, 'P1': 4.0, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'C5': {'Qda': 3.1, 'Qd': 3.1, 'Qa': 3.1, 'Q0': 3.1, 'P5': 3.1, 'P4': 3.1, 'P3': 3.5, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C4': {'Qda': 2.7, 'Qd': 2.7, 'Qa': 2.7, 'Q0': 2.7, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.1, 'Nd': 3.1, 'Na': 3.1, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C3': {'Qda': 2.3, 'Qd': 2.3, 'Qa': 2.3, 'Q0': 2.3, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.1, 'P1': 3.5, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C2': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.3, 'P4': 2.3, 'P3': 2.7, 'P2': 2.7, 'P1': 3.1, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.1, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C1': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.0, 'P4': 2.0, 'P3': 2.3, 'P2': 2.3, 'P1': 2.7, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 2.7, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}}

# Calculate MARTINI epsilon from table for P4 (water) interactions
martini_epsilon_table = martini_table['P4']['P4']  # P4-P4 interaction strength (kJ/mol)
martini_epsilon = martini_epsilon_table / 2.914952774272  # Convert to UPSIDE units

# Unit conversion: set 4πϵ₀ = 1
# Each unit charge (1 or -1) corresponds to 21.831807297541 after conversion
charge_conversion_factor = 21.831807297541
coulomb_constant_upside = charge_conversion_factor**2  # This gives us the effective Coulomb constant

# Soft MARTINI LJ parameters for initial minimization (even after pre-separation)
soft_epsilon = martini_epsilon * 0.1  # 10x softer for initial minimization
soft_sigma = martini_sigma * 1.2      # 20% larger for initial minimization

# DOPC lipid topology from MARTINI v2.0
# DOPC bead types: NC3(Q0), PO4(Qa), GL1(Na), GL2(Na), C1A(C1), D2A(C3), C3A(C1), C4A(C1), C1B(C1), D2B(C3), C3B(C1), C4B(C1)
dopc_bead_types = ['Q0', 'Qa', 'Na', 'Na', 'C1', 'C3', 'C1', 'C1', 'C1', 'C3', 'C1', 'C1']
dopc_charges = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# DOPC bonds (1-indexed, convert to 0-indexed)
dopc_bonds = [
    (0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (5, 6), (6, 7), (3, 8), (8, 9), (9, 10), (10, 11)
]

# Bond parameters (length in nm, force constant in kJ/mol/nm^2)
bond_lengths_nm = [0.47, 0.47, 0.37, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47]  # nm
bond_force_constants_martini = [1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250]  # kJ/mol/nm^2

# Convert MARTINI units to UPSIDE units:
# Bond lengths: nm -> Angstroms (multiply by 10)
bond_lengths = [L * 10.0 for L in bond_lengths_nm]  # Convert nm to Angstroms

# Bond force constants: kJ/(mol·nm^2) -> E_up/Å^2
# Need to convert both energy (kJ/mol -> E_up) and length (nm^2 -> Å^2)
# Since 1 nm = 10 Å, then 1 nm^2 = 100 Å^2
# Conversion: kJ/(mol·nm^2) = kJ/(mol·100Å^2) = (kJ/mol)/100Å^2
# Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)
bond_force_constants = [k / (2.914952774272 * 100.0) for k in bond_force_constants_martini]  # E_up/Å^2

# MARTINI angle parameters for DOPC (from topology file)
# Format: [atom1, atom2, atom3] where atom2 is the central atom
angle_atoms_martini = [
    [1, 2, 3],   # PO4-GL1-GL2 (120.0°)
    [1, 2, 4],   # PO4-GL1-C1A (180.0°)
    [2, 4, 5],   # GL1-C1A-D2A (180.0°)
    [4, 5, 6],   # C1A-D2A-C3A (120.0°)
    [5, 6, 7],   # D2A-C3A-C4A (180.0°)
    [3, 8, 9],   # GL2-C1B-D2B (180.0°)
    [8, 9, 10],  # C1B-D2B-C3B (120.0°)
    [9, 10, 11]  # D2B-C3B-C4B (180.0°)
]
angle_values_martini = [120.0, 180.0, 180.0, 120.0, 180.0, 180.0, 120.0, 180.0]  # degrees
angle_force_constants_martini = [25.0, 25.0, 25.0, 45.0, 25.0, 25.0, 45.0, 25.0]  # kJ/mol/rad²

# Convert angle parameters to UPSIDE units
# Note: AngleSpring uses dot product (cosine) as equilibrium value and implements E = 1/2*K*(cos(theta)-cos(theta_0))^2
import math
angle_equil_dp = [math.cos(math.radians(angle)) for angle in angle_values_martini]  # cosine of angle

# Angle force constants: kJ/(mol·rad^2) -> E_up/rad^2
# Only need energy conversion (angles are dimensionless)
# Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)
angle_force_constants = [k / 2.914952774272 for k in angle_force_constants_martini]  # E_up/rad²

# Print unit conversion summary for verification
print(f"\n=== MARTINI Unit Conversions ===")
print(f"Bond lengths (nm -> Å): {bond_lengths_nm[0]:.2f} nm -> {bond_lengths[0]:.1f} Å")
print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): {bond_force_constants_martini[0]} -> {bond_force_constants[0]:.6f}")
print(f"Angle equilibrium (degrees -> cosine): {angle_values_martini[0]}° -> {angle_equil_dp[0]:.6f}")
print(f"Angle force constants (kJ/mol/rad² -> E_up/rad²): {angle_force_constants_martini[0]} -> {angle_force_constants[0]:.6f}")
print(f"Coulomb constant (unit conversion): {coulomb_constant_upside:.1f} (charge_factor² = {charge_conversion_factor:.1f}²)")
print(f"Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)")
print(f"Length conversion: 1 nm = 10 Å, so 1 nm² = 100 Å²")

# Two-stage minimization option to preserve bilayer structure
skip_minimization = False  # Enable minimization for better structure

# Enable all interactions now that PBC wrapping is fixed
skip_bonds_and_angles = False  # Re-enable bonded interactions

# Stage 1: Hold lipids, minimize water/ions
stage1_steps = 100
stage1_T = 0.1   # Much higher temperature for better equilibration
stage1_dt = 0.001  # Larger time step for faster minimization
stage1_frame_interval = 10  # Reasonable frame interval

# Stage 2: Minimize everything gently
stage2_steps = 200
stage2_T = 0.05   # Gentle but not too low
stage2_dt = 0.001  # Larger time step for faster minimization
stage2_frame_interval = 20  # Reasonable frame interval

# These prints will be moved to after box dimensions are read from PDB

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
residue_ids = []
atom_names = []

# Enhanced mapping for different molecule types
pdb_to_martini = {
    # DOPC lipid beads
    'NC3': 'Q0', 'PO4': 'Qa', 'GL1': 'Na', 'GL2': 'Na', 
    'C1A': 'C1', 'D2A': 'C3', 'C3A': 'C1', 'C4A': 'C1',
    'C1B': 'C1', 'D2B': 'C3', 'C3B': 'C1', 'C4B': 'C1',
    # Water and ions
    'W': 'P4', 'NA': 'Qd', 'CL': 'Qa'
}

martini_charges = {
    'Q0': 1.0, 'Qa': -1.0, 'Na': 0.0, 'C1': 0.0, 'C3': 0.0,
    'P4': 0.0, 'Qd': 1.0
}

if not os.path.exists(input_pdb_file):
    print(f"Error: Input PDB file '{input_pdb_file}' not found!")
    sys.exit(1)

# Read box dimensions from CRYST1 record in PDB file
print(f"Reading box dimensions from {input_pdb_file}...")
with open(input_pdb_file, 'r') as f:
    for line in f:
        if line.startswith('CRYST1'):
            # CRYST1 format: CRYST1    a    b    c  alpha  beta  gamma spacegroup
            # Extract a, b, c (box dimensions in Angstroms)
            fields = line.split()
            if len(fields) >= 4:
                x_len = float(fields[1])  # a dimension
                y_len = float(fields[2])  # b dimension  
                z_len = float(fields[3])  # c dimension
                print(f"Found CRYST1 record: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
                break
    else:
        # No CRYST1 record found - use default values
        x_len = 50.0
        y_len = 50.0
        z_len = 50.0
        print(f"WARNING: No CRYST1 record found in PDB file!")
        print(f"Using default box dimensions: X={x_len:.1f}, Y={y_len:.1f}, Z={z_len:.1f} Angstroms")

# Verify box dimensions are set
if x_len is None or y_len is None or z_len is None:
    print("ERROR: Failed to set box dimensions!")
    sys.exit(1)

# Print system parameters now that box dimensions are known
print(f"\n=== System Parameters ===")
print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
print(f"Box volume: {x_len * y_len * z_len:.1f} Å³")
print(f"P4-P4 epsilon from MARTINI table: {martini_epsilon_table:.1f} kJ/mol")
print(f"Converted to UPSIDE units: {martini_epsilon:.6f}")
print(f"Equivalent in kJ/mol: {martini_epsilon * 2.332:.3f} kJ/mol")
print(f"Temperature: {T:.2f} UPSIDE units = {T * 350.59:.1f} K")
print(f"Thermostat timescale: {thermostat_timescale:.1f} (Langevin damping)")

# Check if box dimensions are appropriate for bilayer systems
if x_len > 60 and y_len > 60:
    print("Large box detected - suitable for extended bilayer systems")
elif x_len < 30 or y_len < 30:
    print("Small box detected - suitable for small systems or micelles")
else:
    print("Medium box detected - suitable for standard bilayer patches")

with open(input_pdb_file, 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            
            pdb_name = line[12:16].strip().upper()
            residue_name = line[17:21].strip()
            residue_id = int(line[22:26])
            
            atom_names.append(pdb_name)
            residue_ids.append(residue_id)
            
            # Assign MARTINI type
            if pdb_name in pdb_to_martini:
                mtype = pdb_to_martini[pdb_name]
            else:
                mtype = 'P4'  # Default to P4 (water)
            atom_types.append(mtype)
            # Apply charge conversion factor to convert to UPSIDE units
            raw_charge = martini_charges.get(mtype, 0.0)
            converted_charge = raw_charge * charge_conversion_factor
            charges.append(converted_charge)

if not initial_positions:
    print("Error: No atoms found in PDB file!")
    sys.exit(1)

initial_positions = np.array(initial_positions)
n_atoms = len(initial_positions)
atom_types = np.array(atom_types)
charges = np.array(charges)
residue_ids = np.array(residue_ids)
atom_names = np.array(atom_names)

print(f"Loaded {n_atoms} atoms from PDB")
print(f"Found {len(set(residue_ids))} residues")

# Count different molecule types based on residue analysis
dopc_residues = []
water_residues = []
ion_residues = []

# Classify residues properly
for rid in set(residue_ids):
    # Get all atoms for this residue
    residue_atoms = np.where(residue_ids == rid)[0]
    residue_atom_names = [atom_names[i] for i in residue_atoms]
    
    # Check if this residue contains DOPC atoms
    if any(name in ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'D2A', 'C3A', 'C4A', 'C1B', 'D2B', 'C3B', 'C4B'] for name in residue_atom_names):
        dopc_residues.append(rid)
    elif all(name == 'W' for name in residue_atom_names):
        water_residues.append(rid)
    elif all(name in ['NA', 'CL'] for name in residue_atom_names):
        ion_residues.append(rid)

print(f"DOPC lipids: {len(dopc_residues)}")
print(f"Water molecules: {len(water_residues)}")
print(f"Ions: {len(ion_residues)}")

# Check initial particle range
pos_min = np.min(initial_positions, axis=0)
pos_max = np.max(initial_positions, axis=0)
pos_range = pos_max - pos_min
print(f"Initial position range: X=[{pos_min[0]:.1f}, {pos_max[0]:.1f}], Y=[{pos_min[1]:.1f}, {pos_max[1]:.1f}], Z=[{pos_min[2]:.1f}, {pos_max[2]:.1f}]")
print(f"Position range: X={pos_range[0]:.1f}, Y={pos_range[1]:.1f}, Z={pos_range[2]:.1f} Angstroms")

# Apply periodic boundary conditions to initial positions
# Wrap particles that are outside the box back into the box
print(f"Box dimensions: X=[0, {x_len:.1f}], Y=[0, {y_len:.1f}], Z=[0, {z_len:.1f}] Angstroms")

# Check initial particle range before wrapping
pos_min = np.min(initial_positions, axis=0)
pos_max = np.max(initial_positions, axis=0)
print(f"Initial particle range BEFORE wrapping:")
print(f"  X=[{pos_min[0]:.1f}, {pos_max[0]:.1f}] (range: {pos_max[0]-pos_min[0]:.1f} Å)")
print(f"  Y=[{pos_min[1]:.1f}, {pos_max[1]:.1f}] (range: {pos_max[1]-pos_min[1]:.1f} Å)")
print(f"  Z=[{pos_min[2]:.1f}, {pos_max[2]:.1f}] (range: {pos_max[2]-pos_min[2]:.1f} Å)")

# Apply periodic boundary conditions to wrap particles into the box
# For a box from [0, L], wrap coordinates using: x = x - L * floor(x / L)
wrapped_positions = initial_positions.copy()

# X dimension
wrapped_positions[:, 0] = wrapped_positions[:, 0] - x_len * np.floor(wrapped_positions[:, 0] / x_len)
# Y dimension  
wrapped_positions[:, 1] = wrapped_positions[:, 1] - y_len * np.floor(wrapped_positions[:, 1] / y_len)

# Z dimension - preserve bilayer structure by simple shift instead of wrapping
# The PDB has correct bilayer structure, just shift to positive coordinates
z_min = np.min(initial_positions[:, 2])
if z_min < 0:
    z_shift = -z_min + 1.0  # Shift by 1 Å to ensure positive coordinates
    wrapped_positions[:, 2] = initial_positions[:, 2] + z_shift
    print(f"Shifted Z coordinates by {z_shift:.1f} Å to preserve bilayer structure")
else:
    wrapped_positions[:, 2] = initial_positions[:, 2]

# Update initial_positions with wrapped coordinates
initial_positions = wrapped_positions

# Check particle range after wrapping
pos_min_wrapped = np.min(initial_positions, axis=0)
pos_max_wrapped = np.max(initial_positions, axis=0)
print(f"Particle range AFTER wrapping:")
print(f"  X=[{pos_min_wrapped[0]:.1f}, {pos_max_wrapped[0]:.1f}] (range: {pos_max_wrapped[0]-pos_min_wrapped[0]:.1f} Å)")
print(f"  Y=[{pos_min_wrapped[1]:.1f}, {pos_max_wrapped[1]:.1f}] (range: {pos_max_wrapped[1]-pos_min_wrapped[1]:.1f} Å)")
print(f"  Z=[{pos_min_wrapped[2]:.1f}, {pos_max_wrapped[2]:.1f}] (range: {pos_max_wrapped[2]-pos_min_wrapped[2]:.1f} Å)")

# Check for any particles still outside the box (should not happen after wrapping)
outside_x = np.any((initial_positions[:, 0] < 0) | (initial_positions[:, 0] >= x_len))
outside_y = np.any((initial_positions[:, 1] < 0) | (initial_positions[:, 1] >= y_len))
outside_z = np.any((initial_positions[:, 2] < 0) | (initial_positions[:, 2] >= z_len))

if outside_x or outside_y or outside_z:
    print(f"ERROR: Some particles still outside box after wrapping!")
else:
    print(f"All particles now within box boundaries")

# Check for severe overlaps and pre-separate particles
print("\n=== OVERLAP ANALYSIS AND PRE-SEPARATION ===")
min_distances = []
severe_overlaps = []
overlap_threshold = 2.0  # Severe overlap if distance < 2.0 Å

for i in range(n_atoms):
    for j in range(i+1, n_atoms):
        dx = initial_positions[i, 0] - initial_positions[j, 0]
        dy = initial_positions[i, 1] - initial_positions[j, 1] 
        dz = initial_positions[i, 2] - initial_positions[j, 2]
        
        # Apply minimum image convention for distance calculation
        if dx > x_len/2: dx -= x_len
        elif dx < -x_len/2: dx += x_len
        if dy > y_len/2: dy -= y_len
        elif dy < -y_len/2: dy += y_len
        if dz > z_len/2: dz -= z_len
        elif dz < -z_len/2: dz += z_len
        
        distance = np.sqrt(dx*dx + dy*dy + dz*dz)
        min_distances.append(distance)
        
        if distance < overlap_threshold:
            severe_overlaps.append((i, j, distance))

print(f"Minimum inter-particle distance: {min(min_distances):.3f} Å")
print(f"Number of severe overlaps (< {overlap_threshold} Å): {len(severe_overlaps)}")

if severe_overlaps:
    print(f"WARNING: {len(severe_overlaps)} severe overlaps detected - will perform pre-separation")
    print("  First 10 overlaps:")
    for i, j, dist in severe_overlaps[:10]:
        print(f"    Particles {i} ({atom_types[i]}) and {j} ({atom_types[j]}): {dist:.3f} Å")
else:
    print("No severe overlaps detected - system should be stable")

# Pre-separate overlapping particles if any severe overlaps were found
# DISABLED: Overlap prevention is turned off
if severe_overlaps:
    print("\n=== OVERLAP PREVENTION DISABLED ===")
    print(f"WARNING: {len(severe_overlaps)} severe overlaps detected but pre-separation is disabled")
    print("  First 10 overlaps:")
    for i, j, dist in severe_overlaps[:10]:
        print(f"    Particles {i} ({atom_types[i]}) and {j} ({atom_types[j]}): {dist:.3f} Å")
    print("  System will proceed with overlaps - this may cause instability!")
    print("  Consider enabling pre-separation if simulation fails")
else:
    print("No severe overlaps detected - system should be stable")

# Bilayer structure analysis if lipids are present
if dopc_residues:
    print("\n=== Bilayer Structure Analysis ===")
    nc3_atoms = [i for i, name in enumerate(atom_names) if name == 'NC3']
    po4_atoms = [i for i, name in enumerate(atom_names) if name == 'PO4']
    
    if nc3_atoms and po4_atoms:
        nc3_z = [initial_positions[i][2] for i in nc3_atoms]
        po4_z = [initial_positions[i][2] for i in po4_atoms]
        
        print(f"NC3 headgroups: {len(nc3_atoms)} atoms, Z range: [{min(nc3_z):.1f}, {max(nc3_z):.1f}] Å")
        print(f"PO4 headgroups: {len(po4_atoms)} atoms, Z range: [{min(po4_z):.1f}, {max(po4_z):.1f}] Å")
        
        # Check if bilayer is properly formed
        bilayer_center = np.mean(initial_positions[:, 2])
        nc3_center = np.mean(nc3_z)
        po4_center = np.mean(po4_z)
        
        print(f"System center: {bilayer_center:.1f} Å")
        print(f"NC3 center: {nc3_center:.1f} Å")
        print(f"PO4 center: {po4_center:.1f} Å")
        
        # Estimate bilayer thickness
        bilayer_thickness = abs(nc3_center - po4_center) + 10  # Approximate
        print(f"Estimated bilayer thickness: {bilayer_thickness:.1f} Å")
        
        if bilayer_thickness > 30 and bilayer_thickness < 60:
            print("Bilayer structure appears reasonable")
        else:
            print("WARNING: Bilayer structure may need attention")
            print("   Consider setting skip_minimization = True")
    else:
        print("WARNING: Could not identify bilayer headgroups")

# Create bonds list for DOPC lipids
bonds_list = []
bond_lengths_list = []
bond_force_constants_list = []

for lipid_id in dopc_residues:
    # Find atoms belonging to this lipid
    lipid_atoms = np.where(residue_ids == lipid_id)[0]
    # Filter to only include DOPC atoms in the right order
    dopc_atoms = [i for i in lipid_atoms if atom_names[i] in ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'D2A', 'C3A', 'C4A', 'C1B', 'D2B', 'C3B', 'C4B']]
    if len(dopc_atoms) == 12:  # DOPC has 12 beads
        # Add bonds for this lipid
        for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
            atom1 = dopc_atoms[bond_idx1]
            atom2 = dopc_atoms[bond_idx2]
            bonds_list.append([atom1, atom2])
            bond_lengths_list.append(bond_lengths[i])
            bond_force_constants_list.append(bond_force_constants[i])

print(f"Created {len(bonds_list)} bonds for {len(dopc_residues)} DOPC lipids")

# Create angles list for DOPC lipids
angles_list = []
angle_equil_dp_list = []
angle_force_constants_list = []

for lipid_id in dopc_residues:
    # Find atoms belonging to this lipid
    lipid_atoms = np.where(residue_ids == lipid_id)[0]
    # Filter to only include DOPC atoms in the right order
    dopc_atoms = [i for i in lipid_atoms if atom_names[i] in ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'D2A', 'C3A', 'C4A', 'C1B', 'D2B', 'C3B', 'C4B']]
    if len(dopc_atoms) == 12:  # DOPC has 12 beads
        # Add angles for this lipid
        for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini):
            atom1 = dopc_atoms[angle_idx1]
            atom2 = dopc_atoms[angle_idx2]  # Central atom
            atom3 = dopc_atoms[angle_idx3]
            angles_list.append([atom1, atom2, atom3])
            angle_equil_dp_list.append(angle_equil_dp[i])
            angle_force_constants_list.append(angle_force_constants[i])

print(f"Created {len(angles_list)} angles for {len(dopc_residues)} DOPC lipids")

# --- MINIMIZATION: Create minimization input with softened parameters ---
print("\nCreating minimization input with pre-separated positions...")
print(f"Soft epsilon: {soft_epsilon:.6f} (10x softer than {martini_epsilon:.6f})")
print(f"Soft sigma: {soft_sigma:.2f} Angstroms (20% larger than {martini_sigma:.2f})")
print(f"Coulomb interactions: ENABLED")
if skip_bonds_and_angles:
    print("DEBUG: Using only non-bonded interactions (bonds/angles skipped)")
else:
    print("Using full MARTINI force field (LJ + Coulomb + bonds + angles)")

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
    mass = np.ones(n_atoms, dtype='f4') * 1.0  # Standard MARTINI mass
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
    
    # Store residue IDs for VTF output and topology handling
    residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
    residue_array._v_attrs.arguments = np.array([b'residue_ids'])
    residue_array._v_attrs.shape = residue_ids.shape
    residue_array._v_attrs.n_atoms = n_atoms
    residue_array._v_attrs.initialized = True
    potential_grp = t.create_group(input_grp, 'potential')
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = soft_epsilon  # Use soft epsilon for minimization
    martini_group._v_attrs.sigma = soft_sigma      # Use soft sigma for minimization
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.coulomb_constant = 1.0  # 4πε₀ = 1 in UPSIDE units
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 0  # Disable force capping for gentler minimization
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
    
    # Create pairs and coefficients for non-bonded interactions
    # Exclude 1-2, 1-3, 1-4 interactions as per MARTINI nrexcl=1
    pairs_list = []
    coeff_array = []
    
    # Create sets of bonded pairs for exclusions
    bonded_pairs_12 = set()  # Directly bonded (1-2)
    bonded_pairs_13 = set()  # 1-3 interactions  
    bonded_pairs_14 = set()  # 1-4 interactions
    
    # Add 1-2 exclusions from bond list
    if not skip_bonds_and_angles:
        for bond in bonds_list:
            bonded_pairs_12.add((min(bond[0], bond[1]), max(bond[0], bond[1])))
    
    # Add 1-3 exclusions from angle list
    if not skip_bonds_and_angles:
        for angle in angles_list:
            # Add 1-3 pairs from angles
            bonded_pairs_13.add((min(angle[0], angle[2]), max(angle[0], angle[2])))
    
    # For particles without explicit bonds, exclude very close particles as 1-2
    close_threshold = 2.5  # Å - particles closer than this are considered "bonded"
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Skip if already in bonded sets
            if (i, j) in bonded_pairs_12 or (i, j) in bonded_pairs_13 or (i, j) in bonded_pairs_14:
                continue
            
            # Check distance for implicit 1-2 exclusions
            dx = initial_positions[i, 0] - initial_positions[j, 0]
            dy = initial_positions[i, 1] - initial_positions[j, 1] 
            dz = initial_positions[i, 2] - initial_positions[j, 2]
            
            # Apply minimum image convention
            if dx > x_len/2: dx -= x_len
            elif dx < -x_len/2: dx += x_len
            if dy > y_len/2: dy -= y_len
            elif dy < -y_len/2: dy += y_len
            if dz > z_len/2: dz -= z_len
            elif dz < -z_len/2: dz += z_len
            
            distance = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            # Exclude very close particles as 1-2 interactions
            if distance < close_threshold:
                bonded_pairs_12.add((i, j))
                continue
            
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon from table and convert to simulation units
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            # Use soft parameters for initial minimization
            epsilon = epsilon_sim * 0.1  # Use soft epsilon
            sigma_val = soft_sigma  # Use soft sigma
            q1 = charges[i]  # Already converted to UPSIDE units
            q2 = charges[j]  # Already converted to UPSIDE units
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True
    
    # Add bond potential for DOPC lipids using MARTINI dist_spring
    if bonds_list and not skip_bonds_and_angles:
        bond_group = t.create_group(potential_grp, 'dist_spring')
        bond_group._v_attrs.arguments = np.array([b'pos'])
        bond_group._v_attrs.initialized = True
        bond_group._v_attrs.x_len = x_len
        bond_group._v_attrs.y_len = y_len
        bond_group._v_attrs.z_len = z_len
        
        bonds_array = np.array(bonds_list, dtype=int)
        bond_lengths_array = np.array(bond_lengths_list, dtype='f4')
        bond_force_constants_array = np.array(bond_force_constants_list, dtype='f4')
        
        # Create bonded_atoms array (1 for bonded, 0 for non-bonded)
        bonded_atoms = np.ones(len(bonds_list), dtype=int)
        
        t.create_array(bond_group, 'id', obj=bonds_array)
        t.create_array(bond_group, 'equil_dist', obj=bond_lengths_array)
        t.create_array(bond_group, 'spring_const', obj=bond_force_constants_array)
        t.create_array(bond_group, 'bonded_atoms', obj=bonded_atoms)
    
    # Add angle potential for DOPC lipids using MARTINI angle_spring
    if angles_list and not skip_bonds_and_angles:
        angle_group = t.create_group(potential_grp, 'angle_spring')
        angle_group._v_attrs.arguments = np.array([b'pos'])
        angle_group._v_attrs.initialized = True
        angle_group._v_attrs.x_len = x_len
        angle_group._v_attrs.y_len = y_len
        angle_group._v_attrs.z_len = z_len
        
        angles_array = np.array(angles_list, dtype=int)
        angle_equil_dp_array = np.array(angle_equil_dp_list, dtype='f4')
        angle_force_constants_array = np.array(angle_force_constants_list, dtype='f4')
        
        t.create_array(angle_group, 'id', obj=angles_array)
        t.create_array(angle_group, 'equil_dist', obj=angle_equil_dp_array)
        t.create_array(angle_group, 'spring_const', obj=angle_force_constants_array)

# Count how many bonded interactions were actually added
bonds_added = len(bonds_list) if bonds_list and not skip_bonds_and_angles else 0
angles_added = len(angles_list) if angles_list and not skip_bonds_and_angles else 0

if skip_bonds_and_angles:
    print(f"DEBUG: Skipped all bonded interactions for debugging")

print(f"Created minimization input with {len(pairs_list)} non-bonded pairs and {bonds_added} bonds and {angles_added} angles")

# --- RUN MINIMIZATION ---
if skip_minimization:
    print("\nSkipping minimization as requested")
    minimized_positions = initial_positions  # Use initial positions directly
    min_h5_file = None  # No minimization file to convert
else:
    min_run_dir = f"{output_dir}/minimize"
    if not os.path.exists(min_run_dir):
        os.makedirs(min_run_dir)
    
    if dopc_residues:
        # Two-stage minimization for bilayer systems
        print("\n=== Two-Stage Minimization for Bilayer System ===")
        
        # Stage 1: Hold lipids fixed, minimize water and ions only
        print(f"\nStage 1: Minimizing water/ions (lipids held fixed)")
        print(f"  Steps: {stage1_steps}")
        print(f"  Temperature: {stage1_T}")
        print(f"  Time step: {stage1_dt}")
        
        # Create stage 1 input with position restraints for lipids
        stage1_input_file = f"{input_dir}/minimize_stage1.up"
        with tb.open_file(min_input_file, 'r') as src:
            with tb.open_file(stage1_input_file, 'w') as dst:
                src.copy_children(src.root, dst.root, recursive=True)
                
                # Add position restraints for lipid atoms
                potential_grp = dst.root.input.potential
                
                # Create position restraint potential for lipids
                restraint_group = dst.create_group(potential_grp, 'position_restraint')
                restraint_group._v_attrs.arguments = np.array([b'pos'])
                restraint_group._v_attrs.initialized = True
                restraint_group._v_attrs.x_len = x_len
                restraint_group._v_attrs.y_len = y_len
                restraint_group._v_attrs.z_len = z_len
                
                # Create arrays for position restraints
                lipid_atoms = []
                restraint_positions = []
                restraint_force_constants = []
                
                for i in range(len(residue_ids)):
                    if residue_ids[i] in dopc_residues:
                        lipid_atoms.append(i)
                        restraint_positions.append(initial_positions[i])
                        restraint_force_constants.append(1000.0)  # Strong restraint
                
                if lipid_atoms:
                    lipid_atoms_array = np.array(lipid_atoms, dtype=int)
                    restraint_positions_array = np.array(restraint_positions, dtype='f4')
                    restraint_force_constants_array = np.array(restraint_force_constants, dtype='f4')
                    
                    dst.create_array(restraint_group, 'atom_indices', obj=lipid_atoms_array)
                    dst.create_array(restraint_group, 'equilibrium_positions', obj=restraint_positions_array)
                    dst.create_array(restraint_group, 'force_constants', obj=restraint_force_constants_array)
                
                # Print debug info about freezing
                frozen_count = len(lipid_atoms)
                print(f"Stage 1: Position-restrained {frozen_count} lipid atoms with force constant 1000.0")
                print(f"Stage 1: {len(residue_ids) - frozen_count} water/ion atoms free to move")
        
        stage1_h5_file = f"{min_run_dir}/stage1.run.up"
        with tb.open_file(stage1_input_file, 'r') as src:
            with tb.open_file(stage1_h5_file, 'w') as dst:
                src.copy_children(src.root, dst.root, recursive=True)
        
        # Run stage 1
        stage1_opts = f"{stage1_h5_file} --duration {stage1_steps} --frame-interval {stage1_frame_interval} --temperature {stage1_T} --time-step {stage1_dt} --thermostat-timescale {thermostat_timescale} --seed 12345"
        cmd = f"{upside_path}/obj/upside {stage1_opts}"
        print(f"Stage 1 command: {cmd}")
        result = sp.run(cmd, shell=True)
        
        if result.returncode != 0:
            print("Stage 1 minimization failed!")
            sys.exit(1)
        else:
            print("Stage 1 completed successfully!")
        
        # Extract stage 1 positions
        with tb.open_file(stage1_h5_file, 'r') as t:
            stage1_positions = t.root.output.pos[-1].copy()
            # Reshape to (n_atoms, 3) if needed
            if len(stage1_positions.shape) == 3:
                stage1_positions = stage1_positions[0].reshape(len(residue_ids), 3)
        
        # Verify bilayer structure is preserved in stage 1
        print(f"\n=== Stage 1 Bilayer Structure Verification ===")
        lipid_drift = []
        water_drift = []
        
        for i in range(len(residue_ids)):
            initial_pos = initial_positions[i]
            stage1_pos = stage1_positions[i]
            drift = np.linalg.norm(stage1_pos - initial_pos)
            
            if residue_ids[i] in dopc_residues:
                lipid_drift.append(drift)
            else:
                water_drift.append(drift)
        
        if lipid_drift:
            max_lipid_drift = max(lipid_drift)
            mean_lipid_drift = np.mean(lipid_drift)
            print(f"Lipid atoms: max drift = {max_lipid_drift:.6f} Å, mean drift = {mean_lipid_drift:.6f} Å")
            if max_lipid_drift > 0.001:
                print(f"WARNING: Lipid atoms moved during stage 1! This indicates freezing failed.")
            else:
                print(f"✓ Lipid atoms properly frozen (drift < 0.001 Å)")
        
        if water_drift:
            max_water_drift = max(water_drift)
            mean_water_drift = np.mean(water_drift)
            print(f"Water/ion atoms: max drift = {max_water_drift:.3f} Å, mean drift = {mean_water_drift:.3f} Å")
            print(f"✓ Water/ions relaxed as expected")
        
        # Stage 2: Gentle minimization of everything
        print(f"\nStage 2: Gentle minimization of entire system")
        print(f"  Steps: {stage2_steps}")
        print(f"  Temperature: {stage2_T}")
        print(f"  Time step: {stage2_dt}")
        
        # Create stage 2 input with stage 1 positions and normal masses
        stage2_input_file = f"{input_dir}/minimize_stage2.up"
        with tb.open_file(min_input_file, 'r') as src:
            with tb.open_file(stage2_input_file, 'w') as dst:
                src.copy_children(src.root, dst.root, recursive=True)
                # Update positions with stage 1 results
                pos_data = dst.root.input.pos[:]
                # Ensure stage1_positions is in the right format (n_atoms, 3)
                if len(stage1_positions.shape) == 3:
                    stage1_pos_reshaped = stage1_positions[0].reshape(len(residue_ids), 3)
                else:
                    stage1_pos_reshaped = stage1_positions
                pos_data[:, :, 0] = stage1_pos_reshaped
                dst.root.input.pos[:] = pos_data
        
        stage2_h5_file = f"{min_run_dir}/stage2.run.up"
        with tb.open_file(stage2_input_file, 'r') as src:
            with tb.open_file(stage2_h5_file, 'w') as dst:
                src.copy_children(src.root, dst.root, recursive=True)
        
        # Run stage 2
        stage2_opts = f"{stage2_h5_file} --duration {stage2_steps} --frame-interval {stage2_frame_interval} --temperature {stage2_T} --time-step {stage2_dt} --thermostat-timescale {thermostat_timescale} --seed 12346"
        cmd = f"{upside_path}/obj/upside {stage2_opts}"
        print(f"Stage 2 command: {cmd}")
        result = sp.run(cmd, shell=True)
        
        if result.returncode != 0:
            print("Stage 2 minimization failed!")
            sys.exit(1)
        else:
            print("Stage 2 completed successfully!")
        
        # Extract final minimized positions
        with tb.open_file(stage2_h5_file, 'r') as t:
            minimized_positions = t.root.output.pos[-1].copy()
            # Reshape to (n_atoms, 3) if needed
            if len(minimized_positions.shape) == 3:
                minimized_positions = minimized_positions[0].reshape(len(residue_ids), 3)
        
        print("Two-stage minimization completed!")
        
        # For VTF conversion, we'll use the final stage file
        min_h5_file = stage2_h5_file
        
    else:
        # Single-stage minimization for simple systems
        print(f"\nRunning single-stage minimization:")
        print(f"  Steps: {minimization_steps}")
        print(f"  Temperature: {minimization_T}")
        print(f"  Time step: {minimization_dt}")
        
        min_h5_file = f"{min_run_dir}/minimize.run.up"
        
        # Copy minimization input to run directory
        with tb.open_file(min_input_file, 'r') as src:
            with tb.open_file(min_h5_file, 'w') as dst:
                src.copy_children(src.root, dst.root, recursive=True)
        
        # Run minimization simulation
        min_opts = f"{min_h5_file} --duration {minimization_steps} --frame-interval {minimization_frame_interval} --temperature {minimization_T} --time-step {minimization_dt} --thermostat-timescale {thermostat_timescale} --seed 12345"
        cmd = f"{upside_path}/obj/upside {min_opts}"
        
        print(f"Minimization command: {cmd}")
        result = sp.run(cmd, shell=True)
        
        if result.returncode != 0:
            print("Minimization failed!")
            sys.exit(1)
        else:
            print("Minimization completed successfully!")
        
        # Extract minimized positions
        with tb.open_file(min_h5_file, 'r') as t:
            minimized_positions = t.root.output.pos[-1].copy()
            # Reshape to (n_atoms, 3) if needed
            if len(minimized_positions.shape) == 3:
                minimized_positions = minimized_positions[0].reshape(len(residue_ids), 3)
    
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
    mass = np.ones(n_atoms, dtype='f4') * 1.0  # Standard MARTINI mass
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
    
    # Store residue IDs for VTF output and topology handling
    residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
    residue_array._v_attrs.arguments = np.array([b'residue_ids'])
    residue_array._v_attrs.shape = residue_ids.shape
    residue_array._v_attrs.n_atoms = n_atoms
    residue_array._v_attrs.initialized = True
    
    potential_grp = t.create_group(input_grp, 'potential')
    martini_group = t.create_group(potential_grp, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = martini_epsilon  # Use full MARTINI epsilon for MD
    martini_group._v_attrs.sigma = martini_sigma      # Use standard sigma for MD
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.coulomb_constant = 1.0  # 4πε₀ = 1 in UPSIDE units
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 0  # Disable force capping for MD
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
    
    # Create pairs and coefficients for non-bonded interactions
    # Exclude 1-2, 1-3, 1-4 interactions as per MARTINI nrexcl=1
    pairs_list = []
    coeff_array = []
    
    # Create sets of bonded pairs for exclusions
    bonded_pairs_12 = set()  # Directly bonded (1-2)
    bonded_pairs_13 = set()  # 1-3 interactions  
    bonded_pairs_14 = set()  # 1-4 interactions
    
    # Add 1-2 exclusions from bond list
    if not skip_bonds_and_angles:
        for bond in bonds_list:
            bonded_pairs_12.add((min(bond[0], bond[1]), max(bond[0], bond[1])))
    
    # Add 1-3 exclusions from angle list
    if not skip_bonds_and_angles:
        for angle in angles_list:
            # Add 1-3 pairs from angles
            bonded_pairs_13.add((min(angle[0], angle[2]), max(angle[0], angle[2])))
    
    # For particles without explicit bonds, exclude very close particles as 1-2
    close_threshold = 2.5  # Å - particles closer than this are considered "bonded"
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Skip if already in bonded sets
            if (i, j) in bonded_pairs_12 or (i, j) in bonded_pairs_13 or (i, j) in bonded_pairs_14:
                continue
            
            # Check distance for implicit 1-2 exclusions
            dx = initial_positions[i, 0] - initial_positions[j, 0]
            dy = initial_positions[i, 1] - initial_positions[j, 1] 
            dz = initial_positions[i, 2] - initial_positions[j, 2]
            
            # Apply minimum image convention
            if dx > x_len/2: dx -= x_len
            elif dx < -x_len/2: dx += x_len
            if dy > y_len/2: dy -= y_len
            elif dy < -y_len/2: dy += y_len
            if dz > z_len/2: dz -= z_len
            elif dz < -z_len/2: dz += z_len
            
            distance = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            # Exclude very close particles as 1-2 interactions
            if distance < close_threshold:
                bonded_pairs_12.add((i, j))
                continue
            
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon from table and convert to simulation units
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / 2.914952774272
            # Use full MARTINI parameters for production MD
            epsilon = epsilon_sim  # Use full MARTINI epsilon
            sigma_val = martini_sigma  # Use full MARTINI sigma
            q1 = charges[i]  # Already converted to UPSIDE units
            q2 = charges[j]  # Already converted to UPSIDE units
            coeff_array.append([epsilon, sigma_val, q1, q2])
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True
    
    # Add bond potential for DOPC lipids using MARTINI dist_spring
    if bonds_list and not skip_bonds_and_angles:
        bond_group = t.create_group(potential_grp, 'dist_spring')
        bond_group._v_attrs.arguments = np.array([b'pos'])
        bond_group._v_attrs.initialized = True
        bond_group._v_attrs.x_len = x_len
        bond_group._v_attrs.y_len = y_len
        bond_group._v_attrs.z_len = z_len
        
        bonds_array = np.array(bonds_list, dtype=int)
        bond_lengths_array = np.array(bond_lengths_list, dtype='f4')
        bond_force_constants_array = np.array(bond_force_constants_list, dtype='f4')
        
        # Create bonded_atoms array (1 for bonded, 0 for non-bonded)
        bonded_atoms = np.ones(len(bonds_list), dtype=int)
        
        t.create_array(bond_group, 'id', obj=bonds_array)
        t.create_array(bond_group, 'equil_dist', obj=bond_lengths_array)
        t.create_array(bond_group, 'spring_const', obj=bond_force_constants_array)
        t.create_array(bond_group, 'bonded_atoms', obj=bonded_atoms)
    
    # Add angle potential for DOPC lipids using MARTINI angle_spring
    if angles_list and not skip_bonds_and_angles:
        angle_group = t.create_group(potential_grp, 'angle_spring')
        angle_group._v_attrs.arguments = np.array([b'pos'])
        angle_group._v_attrs.initialized = True
        angle_group._v_attrs.x_len = x_len
        angle_group._v_attrs.y_len = y_len
        angle_group._v_attrs.z_len = z_len
        
        angles_array = np.array(angles_list, dtype=int)
        angle_equil_dp_array = np.array(angle_equil_dp_list, dtype='f4')
        angle_force_constants_array = np.array(angle_force_constants_list, dtype='f4')
        
        t.create_array(angle_group, 'id', obj=angles_array)
        t.create_array(angle_group, 'equil_dist', obj=angle_equil_dp_array)
        t.create_array(angle_group, 'spring_const', obj=angle_force_constants_array)

# Count how many bonded interactions were actually added for MD
bonds_added_md = len(bonds_list) if bonds_list and not skip_bonds_and_angles else 0
angles_added_md = len(angles_list) if angles_list and not skip_bonds_and_angles else 0

print(f"Created MD input with {len(pairs_list)} non-bonded pairs and {bonds_added_md} bonds and {angles_added_md} angles using minimized positions")

# Debug: Print first 10 atom types and Z positions in MD input
with tb.open_file(input_file, 'r') as t:
    types = t.root.input.type[:].astype(str)
    pos = t.root.input.pos[:,:,0]
    print("\n[DEBUG] First 10 atoms in MD input (.up):")
    for i in range(10):
        print(f"  Atom {i}: type={types[i]}, Z={pos[i,2]:.2f} Å")

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
if dopc_residues:
    print(f"  Lipid system: {len(dopc_residues)} DOPC lipids with bonds and angles")
else:
    print(f"  Simple water system: {len(water_residues)} water molecules")

# Run production MD simulation
print(f"Running production MD simulation...")
# Use the original duration setting for production

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

#----------------------------------------------------------------------
## Convert trajectories to VTF format for visualization
#----------------------------------------------------------------------

print("\nConverting trajectories to VTF format...")

vtf_script = "extract_martini_vtf.py"

# Convert minimization trajectory if it exists
if min_h5_file is not None:
    print(f"\nConverting minimization trajectory...")
    min_output_vtf = "{}/minimization_trajectory.vtf".format(run_dir)
    min_convert_cmd = "python {} {} {}".format(vtf_script, min_h5_file, min_output_vtf)
    print(f"Minimization VTF conversion command: {min_convert_cmd}")
    
    min_vtf_result = sp.run(min_convert_cmd, shell=True)
    
    if min_vtf_result.returncode != 0:
        print("Minimization VTF conversion failed!")
        print("Return code:", min_vtf_result.returncode)
    else:
        print(f"Minimization VTF trajectory saved to: {min_output_vtf}")

# Convert production MD trajectory
print(f"\nConverting production MD trajectory...")
input_traj = h5_file
output_vtf = "{}/martini_trajectory.vtf".format(run_dir)

# Run VTF conversion
convert_cmd = "python {} {} {}".format(vtf_script, input_traj, output_vtf)
print(f"Production VTF conversion command: {convert_cmd}")

vtf_result = sp.run(convert_cmd, shell=True)

if vtf_result.returncode != 0:
    print("Production VTF conversion failed!")
    print("Return code:", vtf_result.returncode)
else:
    print(f"Production VTF trajectory saved to: {output_vtf}")

print("Done!") 