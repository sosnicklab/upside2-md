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
# martini_sigma will be set from MARTINI 3.00 table
dielectric_constant = 15.0  # MARTINI dielectric constant (standard value)

# Optimization settings - All MARTINI interactions use spline interpolation for maximum performance

# Minimization parameters
minimization_steps = 5000
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
T              = 0.8  # Temperature (UPSIDE units) - increased to see temperature effects  
duration       = 1000  # Total simulation steps (production run)
frame_interval = 20   # Output every N steps (more frames for better trajectory)
dt             = 0.01  # Time step - further reduced for better stability
thermostat_timescale = 0.135  # Thermostat timescale (default Langevin damping)
#thermostat_timescale = 0.001  # Thermostat timescale (default Langevin damping)

# NPT ensemble parameters for lipid bilayer simulation
pressure = 0.0  # Target pressure (UPSIDE units) - 0 for atmospheric pressure
barostat_timescale = 0.001  # Very conservative barostat timescale for pressure coupling
barostat_interval = 1.0  # Less frequent barostat application for stability

# Option to disable barostat for testing (set to True to disable)
disable_barostat = False  # Set to True to run NVT instead of NPT

# Enable thermostat for stable equilibrium
thermostat_interval = 1  # Apply thermostat every step for stability

# Import topology reader for MARTINI 3.00
from read_martini3_topology import read_martini3_atoms, read_martini3_bonds, read_martini3_angles, read_martini3_nonbond_params, read_martini3_masses

# Read MARTINI 3.00 nonbonded parameters
main_itp_file = "ff3.00/martini_v3.0.0.itp"
martini3_table = read_martini3_nonbond_params(main_itp_file)

# Calculate MARTINI 3.00 epsilon and sigma from table for W (water) interactions

energy_conversion_factor = 2.914952774272
# Get W-W interaction from MARTINI 3.00 table (sigma, epsilon)
if ('W', 'W') in martini3_table:
    martini_sigma_table, martini_epsilon_table = martini3_table[('W', 'W')]
    martini_epsilon = martini_epsilon_table / energy_conversion_factor  # Convert to UPSIDE units
    martini_sigma = martini_sigma_table * 10.0  # Convert nm to Angstroms
    print(f"Found W-W interaction in MARTINI 3.00 table: sigma={martini_sigma_table:.3f} nm ({martini_sigma:.1f} Å), epsilon={martini_epsilon_table:.3f} kJ/mol")
else:
    # Fallback to default values if W-W not found
    martini_epsilon_table = 4.65  # Default W-W epsilon from MARTINI 3.00
    martini_sigma_table = 0.47  # Default W-W sigma from MARTINI 3.00 (nm)
    martini_epsilon = martini_epsilon_table / energy_conversion_factor
    martini_sigma = martini_sigma_table * 10.0  # Convert nm to Angstroms
    print(f"W-W interaction not found in table, using default sigma={martini_sigma_table:.3f} nm ({martini_sigma:.1f} Å), epsilon={martini_epsilon_table:.3f} kJ/mol")

# Unit conversion: set 4πϵ₀ = 1
# Each unit charge (1 or -1) corresponds to 21.831807297541 after conversion
charge_conversion_factor = 21.831807297541
coulomb_constant_upside = charge_conversion_factor**2  # This gives us the effective Coulomb constant

# Use full MARTINI parameters (no softening needed since we skip minimization)

# Read DOPC topology from MARTINI 3.00 phospholipids file
itp_file = "ff3.00/martini_v3.0.0_phospholipids_v1.itp"
dopc_bead_types, dopc_charges = read_martini3_atoms(itp_file, "DOPC")

print(f"Read DOPC topology from {itp_file}:")
print(f"  Bead types: {dopc_bead_types}")
print(f"  Charges: {dopc_charges}")

# Read ion topologies from MARTINI 3.00 ions file
ions_itp_file = "ff3.00/martini_v3.0.0_ions_v1.itp"
na_bead_types, na_charges = read_martini3_atoms(ions_itp_file, "NA")
cl_bead_types, cl_charges = read_martini3_atoms(ions_itp_file, "CL")

print(f"Read ion topologies from {ions_itp_file}:")
print(f"  NA bead types: {na_bead_types}, charges: {na_charges}")
print(f"  CL bead types: {cl_bead_types}, charges: {cl_charges}")

# Read water topology from MARTINI 3.00 solvents file
solvents_itp_file = "ff3.00/martini_v3.0.0_solvents_v1.itp"
water_bead_types, water_charges = read_martini3_atoms(solvents_itp_file, "W")

print(f"Read water topology from {solvents_itp_file}:")
print(f"  W bead types: {water_bead_types}, charges: {water_charges}")

# Read bead masses from MARTINI 3.00 main file
martini_masses = read_martini3_masses(main_itp_file)
print(f"Read bead masses from {main_itp_file}:")
for bead_type in ['W', 'TQ5', 'Q1', 'Q5', 'SN4a', 'N4a', 'C1', 'C4h']:
    if bead_type in martini_masses:
        print(f"  {bead_type}: {martini_masses[bead_type]}")

# Read DOPC bond topology from MARTINI 3.00 phospholipids file
dopc_bonds_1indexed, bond_lengths_nm, bond_force_constants_martini = read_martini3_bonds(itp_file, "DOPC")

# Convert 1-indexed to 0-indexed for internal use
dopc_bonds = [(b[0]-1, b[1]-1) for b in dopc_bonds_1indexed]

print(f"Read {len(dopc_bonds)} bonds from MARTINI parameter file")
print(f"Bonds (0-indexed): {dopc_bonds}")

# Bond parameters are now read from the parameter file above

# Convert MARTINI units to UPSIDE units:
# Bond lengths: nm -> Angstroms (multiply by 10)
bond_lengths = [L * 10.0 for L in bond_lengths_nm]  # Convert nm to Angstroms

# Bond force constants: kJ/(mol·nm^2) -> E_up/Å^2
# Need to convert both energy (kJ/mol -> E_up) and length (nm^2 -> Å^2)
# Since 1 nm = 10 Å, then 1 nm^2 = 100 Å^2
# Conversion: kJ/(mol·nm^2) = kJ/(mol·100Å^2) = (kJ/mol)/100Å^2
# Energy conversion factor: 2.914952774272 (kJ/mol -> E_up)
bond_force_constants = [k / (energy_conversion_factor * 100.0) for k in bond_force_constants_martini]  # E_up/Å^2

# Read DOPC angle topology from MARTINI 3.00 phospholipids file
dopc_angles_1indexed, angle_values_martini, angle_force_constants_martini = read_martini3_angles(itp_file, "DOPC")

# Convert 1-indexed to 0-indexed for internal use
angle_atoms_martini = [(a[0]-1, a[1]-1, a[2]-1) for a in dopc_angles_1indexed]

print(f"Read {len(angle_atoms_martini)} angles from MARTINI parameter file")
print(f"Angles (0-indexed): {angle_atoms_martini}")

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
if bond_lengths_nm:
    print(f"Bond lengths (nm -> Å): {bond_lengths_nm[0]:.2f} nm -> {bond_lengths[0]:.1f} Å")
    print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): {bond_force_constants_martini[0]} -> {bond_force_constants[0]:.6f}")
else:
    print("No bonds found - bond parameters not available")
if angle_values_martini:
    print(f"Angle equilibrium (degrees): {angle_values_martini[0]}°")
    print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): {angle_force_constants_martini[0]} -> {angle_force_constants[0]:.6f}")
else:
    print("No angles found - angle parameters not available")
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

# Enhanced mapping for different molecule types (MARTINI 3.00)
pdb_to_martini = {
    # DOPC lipid beads (MARTINI 3.00 bead types)
    'NC3': 'Q1', 'PO4': 'Q5', 'GL1': 'SN4a', 'GL2': 'N4a', 
    'C1A': 'C1', 'D2A': 'C4h', 'C3A': 'C1', 'C4A': 'C1',
    'C1B': 'C1', 'D2B': 'C4h', 'C3B': 'C1', 'C4B': 'C1',
    # Water and ions (MARTINI 3.00) - updated to correct bead types
    'W': 'W', 'NA': 'TQ5', 'CL': 'TQ5'
}

martini_charges = {
    'Q1': 1.0, 'Q5': -1.0, 'SN4a': 0.0, 'N4a': 0.0, 'C1': 0.0, 'C4h': 0.0,
    'W': 0.0, 'TQ5': 0.0  # Charges will be set from topology files
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
            residue_name = line[17:20].strip().upper()
            
            # Map to MARTINI type based on residue name and atom name
            if residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file
                if atom_name in dopc_bead_types:
                    atom_idx = dopc_bead_types.index(atom_name)
                    martini_type = dopc_bead_types[atom_idx]
                    charge = dopc_charges[atom_idx]
                else:
                    # Fallback to mapping
                    martini_type = pdb_to_martini.get(atom_name, 'W')
                    charge = 0.0
            elif residue_name == 'W':
                # Water - use topology from parameter file
                martini_type = water_bead_types[0] if water_bead_types else 'W'
                charge = water_charges[0] if water_charges else 0.0
            elif residue_name == 'NA':
                # Sodium ion - use topology from parameter file
                martini_type = na_bead_types[0] if na_bead_types else 'TQ5'
                charge = na_charges[0] if na_charges else 1.0
            elif residue_name == 'CL':
                # Chloride ion - use topology from parameter file
                martini_type = cl_bead_types[0] if cl_bead_types else 'TQ5'
                charge = cl_charges[0] if cl_charges else -1.0
            else:
                # Fallback to mapping
                martini_type = pdb_to_martini.get(atom_name, 'W')
                charge = martini_charges.get(martini_type, 0.0)
            
            atom_types.append(martini_type)
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
print(f"W-W sigma from MARTINI 3.00 table: {martini_sigma_table:.3f} nm ({martini_sigma:.1f} Å)")
print(f"W-W epsilon from MARTINI 3.00 table: {martini_epsilon_table:.1f} kJ/mol")
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
            if current_resname == 'DOPC' or current_resname == 'DOP':  # Handle truncated residue names
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
    if current_resname == 'DOPC' or current_resname == 'DOP':  # Handle truncated residue names
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

# Create bonds and angles for DOPC lipids using dynamic topology from parameter files
bonds_list = []
bond_lengths_list = []
bond_force_constants_list = []
angles_list = []
angle_equil_deg_list = []
angle_force_constants_list = []

# When building bonds/angles, iterate over dopc_molecules
for mol in dopc_molecules:
    _, atom_names_mol, atom_indices = mol
    # Map atom names to indices for this molecule
    name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}
    
    # Create bonds for this lipid using the actual topology from parameter files
    for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
        # Use 0-indexed bond indices to get atom names from the topology
        if bond_idx1 < len(dopc_bead_types) and bond_idx2 < len(dopc_bead_types):
            atom1_name = atom_names_mol[bond_idx1]  # Use actual atom names from PDB
            atom2_name = atom_names_mol[bond_idx2]  # Use actual atom names from PDB
            atom1 = name_to_idx[atom1_name]
            atom2 = name_to_idx[atom2_name]
            bonds_list.append([atom1, atom2])
            bond_lengths_list.append(bond_lengths[i])
            bond_force_constants_list.append(bond_force_constants[i])
    
    # Create angles for this lipid using the actual topology from parameter files
    for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini):
        # Use 0-indexed angle indices to get atom names from the topology
        if (angle_idx1 < len(dopc_bead_types) and angle_idx2 < len(dopc_bead_types) and 
            angle_idx3 < len(dopc_bead_types)):
            atom1_name = atom_names_mol[angle_idx1]  # Use actual atom names from PDB
            atom2_name = atom_names_mol[angle_idx2]  # Use actual atom names from PDB
            atom3_name = atom_names_mol[angle_idx3]  # Use actual atom names from PDB
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
        if bond_idx1 < len(atom_names_mol) and bond_idx2 < len(atom_names_mol):
            atom1_name = atom_names_mol[bond_idx1]
            atom2_name = atom_names_mol[bond_idx2]
            atom1 = name_to_idx[atom1_name]
            atom2 = name_to_idx[atom2_name]
            print(f"  {i+1}: {atom1_name}({atom1}) - {atom2_name}({atom2})")
    
    print("Angles:")
    for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini[:5]):  # Show first 5 angles
        if (angle_idx1 < len(atom_names_mol) and angle_idx2 < len(atom_names_mol) and 
            angle_idx3 < len(atom_names_mol)):
            atom1_name = atom_names_mol[angle_idx1]
            atom2_name = atom_names_mol[angle_idx2]
            atom3_name = atom_names_mol[angle_idx3]
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
    # Assign masses based on bead types from MARTINI 3.00 topology
    # UPSIDE expects masses relative to a reference mass (72.0 for MARTINI)
    mass = np.zeros(n_atoms, dtype='f4')
    reference_mass = 72.0  # Reference mass for MARTINI
    for i, bead_type in enumerate(atom_types):
        if bead_type in martini_masses:
            # Normalize mass relative to reference mass
            mass[i] = martini_masses[bead_type] / reference_mass
        else:
            # Fallback to default mass if bead type not found
            mass[i] = 1.0  # Normalized to reference mass
            print(f"WARNING: Bead type '{bead_type}' not found in mass table, using default normalized mass 1.0")
    mass_array = t.create_array(input_grp, 'mass', obj=mass)
    mass_array._v_attrs.arguments = np.array([b'mass'])
    mass_array._v_attrs.shape = mass.shape
    mass_array._v_attrs.n_atoms = n_atoms
    mass_array._v_attrs.initialized = True
    
    # Print mass distribution for verification
    unique_masses, mass_counts = np.unique(mass, return_counts=True)
    print(f"Mass distribution:")
    for m, count in zip(unique_masses, mass_counts):
        print(f"  Mass {m}: {count} atoms")
    
    # Print mass for each bead type
    print(f"\nMass per bead type:")
    bead_type_masses = {}
    for i, bead_type in enumerate(atom_types):
        if bead_type not in bead_type_masses:
            bead_type_masses[bead_type] = mass[i]
    
    for bead_type in sorted(bead_type_masses.keys()):
        normalized_mass = bead_type_masses[bead_type]
        actual_mass = normalized_mass * reference_mass
        print(f"  {bead_type}: {normalized_mass:.3f} (normalized), {actual_mass:.1f} (actual mass units)")
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
    
    # Create MARTINI potential with spline interpolation
    martini_group = t.create_group(potential_grp, 'martini_potential')
    print("Using MARTINI potential with spline interpolation for LJ and Coulomb calculations")
    
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'lj_coulomb'
    martini_group._v_attrs.epsilon = martini_epsilon  # Use full MARTINI epsilon for MD
    martini_group._v_attrs.sigma = martini_sigma      # Use standard sigma for MD
    martini_group._v_attrs.lj_cutoff = 12.0
    martini_group._v_attrs.coul_cutoff = 12.0
    martini_group._v_attrs.dielectric = dielectric_constant
    martini_group._v_attrs.coulomb_constant = 476.627809876965  # Exact Coulomb constant as specified
    martini_group._v_attrs.n_types = 1
    martini_group._v_attrs.n_params = 4
    martini_group._v_attrs.cutoff = 12.0
    martini_group._v_attrs.cache_buffer = 1.0
    martini_group._v_attrs.initialized = True
    martini_group._v_attrs.force_cap = 0  # Disable force capping for MD
    martini_group._v_attrs.mass_scale = 1.0 / 72.0  # Mass scaling for MARTINI (reference mass = 72)
    martini_group._v_attrs.x_len = x_len
    martini_group._v_attrs.y_len = y_len
    martini_group._v_attrs.z_len = z_len
    martini_group._v_attrs.debug_mode = 1  # Enable spline debug output
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
            sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
            bonded_pairs_12.add(sorted_bond)
    else:
        print("WARNING: Bonds and angles are disabled - no 1-2 exclusions will be applied!")

    # Generate all unique pairs (i < j), excluding only 1-2 bonded pairs
    excluded_count = 0
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Skip if already in bonded set (1-2 exclusion only)
            if (i, j) in bonded_pairs_12:
                excluded_count += 1
                continue
            pairs_list.append([i, j])
            type_i = atom_types[i]
            type_j = atom_types[j]
            # Get epsilon and sigma from MARTINI 3.00 table and convert to simulation units
            if (type_i, type_j) in martini3_table:
                sigma_table, epsilon_table = martini3_table[(type_i, type_j)]
                # Convert sigma from nm to Angstroms
                sigma_val = sigma_table * 10.0
                epsilon_sim = epsilon_table / energy_conversion_factor
            else:
                # Fallback to default values if interaction not found
                print(f"WARNING: Interaction {type_i}-{type_j} not found in MARTINI 3.00 table, using defaults")
                sigma_val = martini_sigma  # Use default sigma
                epsilon_sim = martini_epsilon  # Use default epsilon
            
            # Use full MARTINI 3.00 parameters for production MD
            epsilon = epsilon_sim
            q1 = charges[i]  # Standard units (1.0, -1.0)
            q2 = charges[j]  # Standard units (1.0, -1.0)
            coeff_array.append([epsilon, sigma_val, q1, q2])
    
    pairs_array = np.array(pairs_list, dtype=int)
    coeff_array = np.array(coeff_array)
    print(f"Excluded {excluded_count} 1-2 bonded pairs from non-bonded interactions")
    pairs_data = t.create_array(martini_group, 'pairs', obj=pairs_array)
    pairs_data._v_attrs.initialized = True
    coeff_data = t.create_array(martini_group, 'coefficients', obj=coeff_array)
    coeff_data._v_attrs.initialized = True
    
    # Add bond potential for DOPC lipids using MARTINI dist_spring
    if bonds_list and not skip_bonds_and_angles:
        bond_group = t.create_group(potential_grp, 'dist_spring')
        print("Using bond potential with spline interpolation")
            
        bond_group._v_attrs.arguments = np.array([b'pos'])
        bond_group._v_attrs.initialized = True
        bond_group._v_attrs.mass_scale = 1.0 / 72.0  # Mass scaling for MARTINI (reference mass = 72)
        bond_group._v_attrs.x_len = x_len
        bond_group._v_attrs.y_len = y_len
        bond_group._v_attrs.z_len = z_len
        bond_group._v_attrs.debug_mode = 1  # Enable spline debug output
        
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
        print("Using angle potential with spline interpolation")
            
        angle_group._v_attrs.arguments = np.array([b'pos'])
        angle_group._v_attrs.initialized = True
        angle_group._v_attrs.mass_scale = 1.0 / 72.0  # Mass scaling for MARTINI (reference mass = 72)
        angle_group._v_attrs.x_len = x_len
        angle_group._v_attrs.y_len = y_len
        angle_group._v_attrs.z_len = z_len
        angle_group._v_attrs.debug_mode = 1  # Enable spline debug output
        
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

print(f"\nRunning NPT MD simulation:")
print(f"  Temperature: {T:.2f} UPSIDE units = {T * 350.59:.1f} K")
print(f"  Pressure: {pressure:.2f} UPSIDE units")
print(f"  Duration: {duration} steps")
print(f"  Time step: {dt}")
print(f"  Frame interval: {frame_interval}")
print(f"  Thermostat timescale: {thermostat_timescale:.1f}")
print(f"  Barostat timescale: {barostat_timescale:.1f}")
print(f"  Barostat interval: {barostat_interval:.1f}")
if dopc_count > 0:
    print(f"  Lipid system: {dopc_count} DOPC lipids with bonds and angles")
else:
    print(f"  Simple water system: {water_count} water molecules")

# Run production MD simulation
print(f"Running production NPT MD simulation...")
# Use the original duration setting for production

# MPI support
use_mpi = False  # Set to False for serial execution (development)
mpi_ranks = 2   # Number of MPI ranks to use (reduced for testing)

# Run simulation
# Choose integrator based on barostat setting
if disable_barostat:
    integrator = "vv"  # Use Velocity Verlet for NVT
    print("Using NVT ensemble (Velocity Verlet integrator)")
else:
    integrator = "npt"  # Use NPT ensemble
    print("Using NPT ensemble (NPT integrator)")

# Choose execution method
if use_mpi:
    print(f"Using MPI with {mpi_ranks} ranks")
    mpi_cmd = f"mpirun -np {mpi_ranks}"
    
    # Create multiple copies of the system for MPI parallelization
    import shutil
    import os
    
    # Create multiple system files
    system_files = []
    for i in range(mpi_ranks):
        copy_file = f".//outputs/martini_test/test.run.{i}.up"
        shutil.copy2(h5_file, copy_file)
        system_files.append(copy_file)
    
    h5_file = " ".join(system_files)
    print(f"Created {mpi_ranks} system copies for MPI parallelization")
else:
    print("Using single-process execution")
    mpi_cmd = ""

upside_opts = (
    "{} "
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--time-step {} "
    "--thermostat-timescale {} "
    "--thermostat-interval {} "
    "--pressure {} "
    "--barostat-timescale {} "
    "--barostat-interval {} "
    "--seed {} "
    "--integrator {} "  # Use selected integrator
    "--disable-initial-thermalization "  # Ensure zero initial velocities
    "--restart-using-momentum "  # Force reading input/mom for initial momenta
)
upside_opts = upside_opts.format(h5_file, duration, frame_interval, T, dt, thermostat_timescale, thermostat_interval, pressure, barostat_timescale, barostat_interval, 12345, integrator)

if use_mpi:
    cmd = "{} {}/obj/upside {}".format(mpi_cmd, upside_path, upside_opts)
else:
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
if use_mpi and mpi_ranks > 1:
    # For MPI, convert only the first system's trajectory
    input_traj = system_files[0] if 'system_files' in locals() else h5_file.split()[0]
else:
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