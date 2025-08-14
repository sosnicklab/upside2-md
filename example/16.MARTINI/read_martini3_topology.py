import os
import re

def read_martini3_nonbond_params(itp_file):
    """
    Read MARTINI 3.00 nonbonded parameters from the main itp file
    Returns a dictionary mapping (type1, type2) -> (sigma, epsilon)
    """
    martini3_table = {}
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI 3.00 parameter file '{itp_file}' not found!")
        return martini3_table
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the nonbond_params section
    in_nonbond_section = False
    for line in lines:
        line = line.strip()
        if line == '[ nonbond_params ]':
            in_nonbond_section = True
            continue
        elif line.startswith('[') and line.endswith(']'):
            in_nonbond_section = False
            continue
        
        if in_nonbond_section and line and not line.startswith(';'):
            # Parse line: type1 type2 func sigma epsilon
            parts = line.split()
            if len(parts) >= 5:
                type1 = parts[0]
                type2 = parts[1]
                sigma = float(parts[3])
                epsilon = float(parts[4])
                
                # Store both directions for symmetric interactions
                martini3_table[(type1, type2)] = (sigma, epsilon)
                martini3_table[(type2, type1)] = (sigma, epsilon)
    
    print(f"Read {len(martini3_table)//2} unique nonbonded parameter pairs from MARTINI 3.00")
    return martini3_table

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
                # For DOPC, only read the first 12 atoms (standard DOPC structure)
                if molecule_name == "DOPC" and len(bead_types) >= 12:
                    break
    
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
                # For DOPC, only read bonds for the first 12 atoms (standard DOPC structure)
                if molecule_name == "DOPC" and (atom1 > 12 or atom2 > 12):
                    break
                # For DOPC, stop after reading 11 bonds (standard DOPC structure)
                if molecule_name == "DOPC" and len(bonds_1indexed) >= 11:
                    break
    
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
                # For DOPC, only read angles for the first 12 atoms (standard DOPC structure)
                if molecule_name == "DOPC" and (atom1 > 12 or atom2 > 12 or atom3 > 12):
                    break
                # For DOPC, stop after reading 8 angles (standard DOPC structure)
                if molecule_name == "DOPC" and len(angles_1indexed) >= 8:
                    break
    
    print(f"Read {len(angles_1indexed)} angles for {molecule_name} from MARTINI 3.00")
    return angles_1indexed, angle_values_deg, angle_force_constants
