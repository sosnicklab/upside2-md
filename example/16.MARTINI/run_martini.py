import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb
from collections import Counter

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
dt             = 0.001  # Time step - further reduced for better stability
thermostat_timescale = 5.0  # Thermostat timescale (default Langevin damping)

# Enable thermostat for stable equilibrium
thermostat_interval = 1  # Apply thermostat every step for stability

martini_table = {'Qda': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qd': {'Qda': 5.6, 'Qd': 5.0, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 4.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qa': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Q0': {'Qda': 4.5, 'Qd': 4.5, 'Qa': 4.5, 'Q0': 3.5, 'P5': 5.0, 'P4': 5.6, 'P3': 5.0, 'P2': 4.5, 'P1': 4.0, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'P5': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.6, 'P1': 5.6, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P4': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.6, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P3': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P2': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P1': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'Nda': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Nd': {'Qda': 5.0, 'Qd': 4.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.0, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Na': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 4.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.0, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'N0': {'Qda': 3.5, 'Qd': 3.5, 'Qa': 3.5, 'Q0': 3.5, 'P5': 3.5, 'P4': 3.5, 'P3': 3.5, 'P2': 4.0, 'P1': 4.0, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'C5': {'Qda': 3.1, 'Qd': 3.1, 'Qa': 3.1, 'Q0': 3.1, 'P5': 3.1, 'P4': 3.1, 'P3': 3.5, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C4': {'Qda': 2.7, 'Qd': 2.7, 'Qa': 2.7, 'Q0': 2.7, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.1, 'Nd': 3.1, 'Na': 3.1, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C3': {'Qda': 2.3, 'Qd': 2.3, 'Qa': 2.3, 'Q0': 2.3, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.1, 'P1': 3.5, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C2': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.3, 'P4': 2.3, 'P3': 2.7, 'P2': 2.7, 'P1': 3.1, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.1, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C1': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.0, 'P4': 2.0, 'P3': 2.3, 'P2': 2.3, 'P1': 2.7, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 2.7, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}}

# Calculate MARTINI epsilon from table for P4 (water) interactions

energy_conversion_factor = 2.914952774272
martini_epsilon_table = martini_table['P4']['P4']  # P4-P4 interaction strength (kJ/mol)
martini_epsilon = martini_epsilon_table / energy_conversion_factor  # Convert to UPSIDE units

# Unit conversion: set 4πϵ₀ = 1
# Each unit charge (1 or -1) corresponds to 21.831807297541 after conversion
charge_conversion_factor = 21.831807297541
coulomb_constant_upside = charge_conversion_factor**2  # This gives us the effective Coulomb constant

# Use full MARTINI parameters (no softening needed since we skip minimization)

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
bond_force_constants = [k / (energy_conversion_factor * 100.0) for k in bond_force_constants_martini]  # E_up/Å^2

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
angle_force_constants_martini = [25.0, 25.0, 25.0, 45.0, 25.0, 25.0, 45.0, 25.0]  # kJ/mol

# Convert angle parameters to UPSIDE units
# Note: Updated AngleSpring now uses degrees directly, not cosines
import math
angle_equil_deg = angle_values_martini  # Keep angles in degrees

# Angle force constants: kJ/(mol·deg^2) -> E_up/deg^2
# Apply energy conversion factor (angles stay in degrees)
# Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)
angle_force_constants = [k / energy_conversion_factor for k in angle_force_constants_martini]  # E_up/deg²

# Print unit conversion summary for verification
print(f"\n=== MARTINI Unit Conversions ===")
print(f"Bond lengths (nm -> Å): {bond_lengths_nm[0]:.2f} nm -> {bond_lengths[0]:.1f} Å")
print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): {bond_force_constants_martini[0]} -> {bond_force_constants[0]:.6f}")
print(f"Angle equilibrium (degrees): {angle_values_martini[0]}°")
print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): {angle_force_constants_martini[0]} -> {angle_force_constants[0]:.6f}")
print(f"Coulomb constant: 476.6 (standard value for unit conversion)")
print(f"Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)")
print(f"Length conversion: 1 nm = 10 Å, so 1 nm² = 100 Å²")

# Skip minimization entirely
skip_minimization = True  # Disable minimization - use original structure

# Enable all interactions now that PBC wrapping is fixed
skip_bonds_and_angles = False  # Re-enable bonded interactions

# Stage 1: Hold lipids, minimize water/ions
stage1_steps = 100
stage1_T = 0.1   # Much higher temperature for better equilibration
stage1_dt = 0.001  # Smaller time step for better stability
stage1_frame_interval = 10  # Reasonable frame interval

# Stage 2: Minimize everything gently
stage2_steps = 200
stage2_T = 0.05   # Gentle but not too low
stage2_dt = 0.001  # Smaller time step for better stability
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

# === Read ATOM records from PDB and populate arrays ===
with open(input_pdb_file, 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            atom_name = line[12:16].strip().upper()
            atom_names.append(atom_name)
            residue_id = int(line[22:26])
            residue_ids.append(residue_id)
            # Map to MARTINI type
            martini_type = pdb_to_martini.get(atom_name, 'P4')  # Default to water if not found
            atom_types.append(martini_type)
            charge = martini_charges.get(martini_type, 0.0)
            charges.append(charge)

# Convert lists to numpy arrays
initial_positions = np.array(initial_positions, dtype=float)
atom_types = np.array(atom_types)
# Keep charges in standard units (1.0, -1.0) - coulomb_constant will handle the conversion
charges = np.array(charges, dtype=float)
residue_ids = np.array(residue_ids, dtype=int)
atom_names = np.array(atom_names)
n_atoms = len(initial_positions)

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

# --- Molecule extraction based on residue ID and residue name ---
# Parse PDB and group atoms into molecules by residue ID and residue name
molecules = []
current_resid = None
current_resname = None
current_atoms = []
current_indices = []
atom_idx = 0

with open(input_pdb_file) as f:
    for line in f:
        if not line.startswith('ATOM'):
            continue
        
        resname = line[17:20].strip().upper()
        resid = int(line[22:26])
        atom_name = line[12:16].strip().upper()
        
        # If we encounter a new residue, save the previous molecule and start a new one
        if current_resid is not None and (resid != current_resid or resname != current_resname):
            # Determine molecule type based on residue name
            if current_resname == 'DOPC':
                moltype = 'DOPC'
            elif current_resname == 'W':
                moltype = 'WATER'
            elif current_resname == 'NA':
                moltype = 'NA'
            elif current_resname == 'CL':
                moltype = 'CL'
            else:
                moltype = 'UNKNOWN'
            
            molecules.append((moltype, current_atoms.copy(), current_indices.copy()))
            current_atoms = []
            current_indices = []
        
        # Add current atom to the molecule
        current_atoms.append(atom_name)
        current_indices.append(atom_idx)
        current_resid = resid
        current_resname = resname
        atom_idx += 1

# Don't forget the last molecule
if current_atoms:
    if current_resname == 'DOPC':
        moltype = 'DOPC'
    elif current_resname == 'W':
        moltype = 'WATER'
    elif current_resname == 'NA':
        moltype = 'NA'
    elif current_resname == 'CL':
        moltype = 'CL'
    else:
        moltype = 'UNKNOWN'
    
    molecules.append((moltype, current_atoms.copy(), current_indices.copy()))

# Now, molecules is a list of (moltype, [atom_names], [pdb_line_indices])
# DOPC extraction: find DOPC molecules
dopc_molecules = [m for m in molecules if m[0] == 'DOPC']

# Count different molecule types based on new molecule parsing
dopc_count = len([m for m in molecules if m[0] == 'DOPC'])
water_count = len([m for m in molecules if m[0] == 'WATER'])
na_count = len([m for m in molecules if m[0] == 'NA'])
cl_count = len([m for m in molecules if m[0] == 'CL'])

print(f"DOPC lipids: {dopc_count}")
print(f"Water molecules: {water_count}")
print(f"Sodium ions: {na_count}")
print(f"Chloride ions: {cl_count}")

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

# Skip overlap checking and pre-separation
print("\n=== SKIPPING OVERLAP CHECKING AND PRE-SEPARATION ===")
print("Using original PDB structure without modification")

# Bilayer structure analysis if lipids are present
if dopc_count > 0:
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

# DOPC bead names (for topology)
DOPC_NAMES = ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'D2A', 'C3A', 'C4A', 'C1B', 'D2B', 'C3B', 'C4B']

# Create bonds and angles for DOPC lipids
bonds_list = []
bond_lengths_list = []
bond_force_constants_list = []
angles_list = []
angle_equil_deg_list = []
angle_force_constants_list = []

# When building bonds/angles, iterate over dopc_molecules
for mol in dopc_molecules:
    _, atom_names_mol, atom_indices = mol
    # Map DOPC_NAMES to atom_indices for this molecule
    name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}
    # Create bonds for this lipid
    for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
        atom1_name = DOPC_NAMES[bond_idx1]
        atom2_name = DOPC_NAMES[bond_idx2]
        atom1 = name_to_idx[atom1_name]
        atom2 = name_to_idx[atom2_name]
        bonds_list.append([atom1, atom2])
        bond_lengths_list.append(bond_lengths[i])
        bond_force_constants_list.append(bond_force_constants[i])
    # Create angles for this lipid
    for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini):
        atom1_name = DOPC_NAMES[angle_idx1]
        atom2_name = DOPC_NAMES[angle_idx2]
        atom3_name = DOPC_NAMES[angle_idx3]
        atom1 = name_to_idx[atom1_name]
        atom2 = name_to_idx[atom2_name]
        atom3 = name_to_idx[atom3_name]
        angles_list.append([atom1, atom2, atom3])
        angle_equil_deg_list.append(angle_equil_deg[i])
        angle_force_constants_list.append(angle_force_constants[i])

print(f"Created {len(bonds_list)} bonds for {dopc_count} DOPC lipids")
print(f"Created {len(angles_list)} angles for {dopc_count} DOPC lipids")

# Show bond and angle structure for first DOPC molecule as example
if dopc_molecules:
    print(f"\n=== DOPC Bond/Angle Structure (example from first molecule) ===")
    first_dopc = dopc_molecules[0]
    _, atom_names_mol, atom_indices = first_dopc
    name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}
    
    print("Bonds:")
    for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds[:5]):  # Show first 5 bonds
        atom1_name = DOPC_NAMES[bond_idx1]
        atom2_name = DOPC_NAMES[bond_idx2]
        atom1 = name_to_idx[atom1_name]
        atom2 = name_to_idx[atom2_name]
        print(f"  {i+1}: {atom1_name}({atom1}) - {atom2_name}({atom2})")
    
    print("Angles:")
    for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini[:5]):  # Show first 5 angles
        atom1_name = DOPC_NAMES[angle_idx1]
        atom2_name = DOPC_NAMES[angle_idx2]
        atom3_name = DOPC_NAMES[angle_idx3]
        atom1 = name_to_idx[atom1_name]
        atom2 = name_to_idx[atom2_name]
        atom3 = name_to_idx[atom3_name]
        print(f"  {i+1}: {atom1_name}({atom1}) - {atom2_name}({atom2}) - {atom3_name}({atom3})")

# Skip minimization entirely - use original structure
print("\n=== SKIPPING MINIMIZATION ===")
print("Using original PDB structure for MD simulation")

# Skip minimization entirely
print("\nSkipping minimization as requested")
minimized_positions = initial_positions  # Use initial positions directly
min_h5_file = None  # No minimization file to convert

# --- Center and wrap positions just before writing the .up file ---
center = np.mean(initial_positions, axis=0)
centered_positions = initial_positions - center
half_box = np.array([x_len/2, y_len/2, z_len/2])
centered_positions = (centered_positions + half_box) % (2*half_box) - half_box
minimized_positions = centered_positions  # Use these for writing to .up file
print("First 10 positions to be written to .up file (centered and wrapped):")
for i in range(min(10, n_atoms)):
    print(f"  Atom {i}: {minimized_positions[i]}")

# --- Assign unique residue IDs for each molecule (for .up file) ---
residue_ids = np.zeros(len(atom_names), dtype=int)
for unique_id, mol in enumerate(molecules):
    _, _, atom_indices = mol
    for idx in atom_indices:
        residue_ids[idx] = unique_id
# Now residue_ids is a per-atom array, with a unique value for each molecule

# Print molecule summary for verification
print(f"\n=== Molecule Summary ===")
mol_counts = {}
for moltype, atom_names_mol, _ in molecules:
    if moltype not in mol_counts:
        mol_counts[moltype] = 0
    mol_counts[moltype] += 1

for moltype, count in mol_counts.items():
    print(f"{moltype}: {count} molecules")

# Compact molecule structure debug output
print(f"\n=== Molecule Structures ===")
for i, (moltype, atom_names_mol, atom_indices) in enumerate(molecules):
    print(f"Molecule {i+1}: {moltype} - Atoms: {atom_names_mol}")

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
    
    # Add mom dataset with zero momenta (explicitly set to zero)
    momentum = np.zeros((n_atoms, 3, 1), dtype='f4')  # Zero momenta for all atoms, shape (n_atoms, 3, 1)
    mom_array = t.create_array(input_grp, 'mom', obj=momentum)
    mom_array._v_attrs.arguments = np.array([b'mom'])
    mom_array._v_attrs.shape = momentum.shape
    mom_array._v_attrs.n_atoms = n_atoms
    mom_array._v_attrs.dim = 3
    mom_array._v_attrs.initialized = True
    mass = np.ones(n_atoms, dtype='f4') * 72.0  # Set mass to 72
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
    martini_group._v_attrs.coulomb_constant = 476.6  # Standard Coulomb constant for unit conversion
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
    # Exclude only 1-2 interactions as per MARTINI nrexcl=1
    pairs_list = []
    coeff_array = []

    # Create set of bonded pairs for 1-2 exclusions
    bonded_pairs_12 = set()  # Directly bonded (1-2)

    # Add 1-2 exclusions from bond list
    if not skip_bonds_and_angles:
        for bond in bonds_list:
            bonded_pairs_12.add((min(bond[0], bond[1]), max(bond[0], bond[1])))

    # Generate all unique pairs (i < j), excluding only 1-2 bonded pairs
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Skip if already in bonded set (1-2 exclusion only)
            if (i, j) in bonded_pairs_12:
                continue
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon from table and convert to simulation units
            epsilon_table = martini_table[type_i][type_j]
            epsilon_sim = epsilon_table / energy_conversion_factor
            # Use full MARTINI parameters for production MD
            epsilon = epsilon_sim  # Use full MARTINI epsilon
            sigma_val = martini_sigma  # Use full MARTINI sigma
            q1 = charges[i]  # Standard units (1.0, -1.0)
            q2 = charges[j]  # Standard units (1.0, -1.0)
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
        angle_equil_deg_array = np.array(angle_equil_deg_list, dtype='f4')
        angle_force_constants_array = np.array(angle_force_constants_list, dtype='f4')
        
        t.create_array(angle_group, 'id', obj=angles_array)
        t.create_array(angle_group, 'equil_angle_deg', obj=angle_equil_deg_array)
        t.create_array(angle_group, 'spring_const', obj=angle_force_constants_array)

# Count how many bonded interactions were actually added for MD
bonds_added_md = len(bonds_list) if bonds_list and not skip_bonds_and_angles else 0
angles_added_md = len(angles_list) if angles_list and not skip_bonds_and_angles else 0

print(f"Created MD input with {len(pairs_list)} non-bonded pairs and {bonds_added_md} bonds and {angles_added_md} angles using minimized positions")


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
if dopc_count > 0:
    print(f"  Lipid system: {dopc_count} DOPC lipids with bonds and angles")
else:
    print(f"  Simple water system: {water_count} water molecules")

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
    "--thermostat-interval {} "
    "--seed {} "
    "--integrator v "  # Changed from 'verlet' to 'v' for standard Verlet
    "--disable-initial-thermalization "  # Ensure zero initial velocities
    "--disable-thermostat "  # Disable thermostat entirely to keep zero velocities
    "--restart-using-momentum"  # Force reading input/mom for initial momenta
)
upside_opts = upside_opts.format(h5_file, duration, frame_interval, T, dt, thermostat_timescale, thermostat_interval, 12345)

cmd = "{}/obj/upside {}".format(upside_path, upside_opts)

print(f"Command: {cmd}")
result = sp.run(cmd, shell=True)

#----------------------------------------------------------------------
## Convert trajectories to VTF format for visualization
#----------------------------------------------------------------------

print("\nConverting trajectories to VTF format...")

vtf_script = "extract_martini_vtf.py"

# Skip minimization VTF conversion since we're not doing minimization

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