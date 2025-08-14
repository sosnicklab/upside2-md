#!/usr/bin/env python3

import os

def read_martini3_nonbond_params(itp_file):
    """
    Read MARTINI 3.00 nonbonded parameters from the main .itp file
    Returns a dictionary mapping (type1, type2) tuples to (sigma, epsilon) values
    """
    martini_table = {}
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI 3.00 parameter file '{itp_file}' not found!")
        return martini_table
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the nonbond_params section
    in_nonbond_params = False
    param_count = 0
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Check for nonbond_params section start
        if line == '[ nonbond_params ]' or line == '[nonbond_params]':
            in_nonbond_params = True
            continue
        elif line.startswith('[') and line.endswith(']') and line != '[ nonbond_params ]' and line != '[nonbond_params]':
            if in_nonbond_params:
                break
            in_nonbond_params = False
            continue
        
        # Parse nonbond_params lines
        if in_nonbond_params and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                type1 = parts[0]
                type2 = parts[1]
                func = int(parts[2])
                sigma = float(parts[3])  # nm
                epsilon = float(parts[4])  # kJ/mol
                
                # Store both orientations for easy lookup
                martini_table[(type1, type2)] = (sigma, epsilon)
                martini_table[(type2, type1)] = (sigma, epsilon)
                param_count += 1
                

    
    print(f"Read {len(martini_table)//2} unique nonbonded parameter pairs from MARTINI 3.00")
    return martini_table

def read_martini3_atoms(itp_file, molecule_name):
    """
    Read MARTINI 3.00 atom types and charges for a specific molecule
    Returns (bead_types, charges)
    """
    bead_types = []
    charges = []
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI 3.00 topology file '{itp_file}' not found!")
        return bead_types, charges
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the molecule section
    in_molecule = False
    in_atoms = False
    found_molecule = False
    
    for line in lines:
        line = line.strip()
        
        # Check for molecule start
        if line == f'[ moleculetype ]' or line == f'[ moleculetype ]':
            if found_molecule:
                break  # Stop after reading the first occurrence of our molecule
            in_molecule = False
            in_atoms = False
        
        # Check if this is our molecule (handle different formats)
        if line == molecule_name or line.startswith(molecule_name + " ") or line.strip().startswith(molecule_name):
            in_molecule = True
            found_molecule = True
            continue
        
        # Stop if we encounter a different molecule name (after finding our molecule)
        if found_molecule and in_molecule and line and not line.startswith(';') and not line.startswith('[') and not line.startswith('@'):
            # Check if this line starts with a different molecule name
            parts = line.split()
            if len(parts) >= 1 and parts[0] != molecule_name and not parts[0].isdigit():
                break  # Stop when we encounter a different molecule
        
        # Check for atoms section within our molecule
        if in_molecule and (line == '[ atoms ]' or line == '[atoms]'):
            in_atoms = True
            continue
        elif in_molecule and line.startswith('[') and line.endswith(']') and line != '[ atoms ]' and line != '[atoms]':
            in_atoms = False
            continue
        
        # Parse atom lines
        if in_molecule and in_atoms and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 7:
                bead_type = parts[1]
                charge = float(parts[6])
                bead_types.append(bead_type)
                charges.append(charge)
                # Read all atoms for the molecule - no hardcoded limits
    
    print(f"Read {len(bead_types)} atoms for {molecule_name} from MARTINI 3.00")
    return bead_types, charges

def read_martini3_bonds(itp_file, molecule_name):
    """
    Read MARTINI 3.00 bonds for a specific molecule
    Returns (bonds_1indexed, bond_lengths_nm, bond_force_constants)
    """
    bonds_1indexed = []
    bond_lengths_nm = []
    bond_force_constants = []
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI 3.00 topology file '{itp_file}' not found!")
        return bonds_1indexed, bond_lengths_nm, bond_force_constants
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the molecule section
    in_molecule = False
    in_bonds = False
    found_molecule = False
    
    for line in lines:
        line = line.strip()
        
        # Check for molecule start
        if line == f'[ moleculetype ]' or line == f'[ moleculetype ]':
            if found_molecule:
                break  # Stop after reading the first occurrence of our molecule
            in_molecule = False
            in_bonds = False
        
        # Check if this is our molecule (handle different formats)
        if line == molecule_name or line.startswith(molecule_name + " ") or line.strip().startswith(molecule_name):
            in_molecule = True
            found_molecule = True
            continue
        
        # Stop if we encounter a different molecule name (after finding our molecule)
        if found_molecule and in_molecule and line and not line.startswith(';') and not line.startswith('[') and not line.startswith('@'):
            # Check if this line starts with a different molecule name
            parts = line.split()
            if len(parts) >= 1 and parts[0] != molecule_name and not parts[0].isdigit():
                break  # Stop when we encounter a different molecule
        
        # Check for bonds section within our molecule
        if in_molecule and (line == '[ bonds ]' or line == '[bonds]'):
            in_bonds = True
            continue
        elif in_molecule and line.startswith('[') and line.endswith(']') and line != '[ bonds ]' and line != '[bonds]':
            in_bonds = False
            continue
        
        # Parse bond lines
        if in_molecule and in_bonds and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                func = int(parts[2])
                length = float(parts[3])
                force_const = float(parts[4])
                
                bonds_1indexed.append([atom1, atom2])
                bond_lengths_nm.append(length)
                bond_force_constants.append(force_const)
                # Read all bonds for the molecule - no hardcoded limits
    
    print(f"Read {len(bonds_1indexed)} bonds for {molecule_name} from MARTINI 3.00")
    return bonds_1indexed, bond_lengths_nm, bond_force_constants

def read_martini3_angles(itp_file, molecule_name):
    """
    Read MARTINI 3.00 angles for a specific molecule
    Returns (angles_1indexed, angle_values_deg, angle_force_constants)
    """
    angles_1indexed = []
    angle_values_deg = []
    angle_force_constants = []
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI 3.00 topology file '{itp_file}' not found!")
        return angles_1indexed, angle_values_deg, angle_force_constants
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the molecule section
    in_molecule = False
    in_angles = False
    found_molecule = False
    
    for line in lines:
        line = line.strip()
        
        # Check for molecule start
        if line == f'[ moleculetype ]' or line == f'[ moleculetype ]':
            if found_molecule:
                break  # Stop after reading the first occurrence of our molecule
            in_molecule = False
            in_angles = False
        
        # Check if this is our molecule (handle different formats)
        if line == molecule_name or line.startswith(molecule_name + " ") or line.strip().startswith(molecule_name):
            in_molecule = True
            found_molecule = True
            continue
        
        # Stop if we encounter a different molecule name (after finding our molecule)
        if found_molecule and in_molecule and line and not line.startswith(';') and not line.startswith('[') and not line.startswith('@'):
            # Check if this line starts with a different molecule name
            parts = line.split()
            if len(parts) >= 1 and parts[0] != molecule_name and not parts[0].isdigit():
                break  # Stop when we encounter a different molecule
        
        # Check for angles section within our molecule
        if in_molecule and (line == '[ angles ]' or line == '[angles]'):
            in_angles = True
            continue
        elif in_molecule and line.startswith('[') and line.endswith(']') and line != '[ angles ]' and line != '[angles]':
            in_angles = False
            continue
        
        # Parse angle lines
        if in_molecule and in_angles and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 6:
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                atom3 = int(parts[2])
                func = int(parts[3])
                angle = float(parts[4])
                force_const = float(parts[5])
                
                angles_1indexed.append([atom1, atom2, atom3])
                angle_values_deg.append(angle)
                angle_force_constants.append(force_const)
                # Read all angles for the molecule - no hardcoded limits
    
    print(f"Read {len(angles_1indexed)} angles for {molecule_name} from MARTINI 3.00")
    return angles_1indexed, angle_values_deg, angle_force_constants
