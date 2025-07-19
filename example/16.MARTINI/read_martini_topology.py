#!/usr/bin/env python3

def read_martini_bonds(itp_file, molecule_name):
    """
    Read bond topology from MARTINI .itp file
    
    Args:
        itp_file (str): Path to the .itp file
        molecule_name (str): Name of the molecule (e.g., 'DOPC')
    
    Returns:
        list: List of bond tuples [(atom1_idx, atom2_idx), ...] (1-indexed)
        list: List of bond lengths in nm
        list: List of bond force constants in kJ/mol/nm^2
    """
    bonds = []
    bond_lengths = []
    bond_force_constants = []
    
    in_molecule = False
    in_bonds_section = False
    
    with open(itp_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith(';'):
                continue
            
            # Check for molecule section
            if line.startswith('[moleculetype]'):
                in_molecule = False
                in_bonds_section = False
                continue
            elif line.startswith('[bonds]'):
                if in_molecule:
                    in_bonds_section = True
                continue
            elif line.startswith('['):
                in_bonds_section = False
                continue
            
            # Check if we're in the right molecule
            if not in_molecule and not in_bonds_section:
                parts = line.split()
                if parts and parts[0] == molecule_name:
                    in_molecule = True
                    print(f"Found molecule {molecule_name} at line {line_num}")
                continue
            
            # Read bond data
            if in_bonds_section and in_molecule:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        length = float(parts[3])
                        force_const = float(parts[4])
                        
                        bonds.append((atom1, atom2))
                        bond_lengths.append(length)
                        bond_force_constants.append(force_const)
                    except (ValueError, IndexError):
                        continue
    
    return bonds, bond_lengths, bond_force_constants

def read_martini_angles(itp_file, molecule_name):
    """
    Read angle topology from MARTINI .itp file
    
    Args:
        itp_file (str): Path to the .itp file
        molecule_name (str): Name of the molecule (e.g., 'DOPC')
    
    Returns:
        list: List of angle tuples [(atom1_idx, atom2_idx, atom3_idx), ...] (1-indexed)
        list: List of angle values in degrees
        list: List of angle force constants in kJ/mol/deg^2
    """
    angles = []
    angle_values = []
    angle_force_constants = []
    
    in_molecule = False
    in_angles_section = False
    
    with open(itp_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith(';'):
                continue
            
            # Check for molecule section
            if line.startswith('[moleculetype]'):
                in_molecule = False
                in_angles_section = False
                continue
            elif line.startswith('[angles]'):
                if in_molecule:
                    in_angles_section = True
                continue
            elif line.startswith('['):
                in_angles_section = False
                continue
            
            # Check if we're in the right molecule
            if not in_molecule and not in_angles_section:
                parts = line.split()
                if parts and parts[0] == molecule_name:
                    in_molecule = True
                    print(f"Found molecule {molecule_name} at line {line_num}")
                continue
            
            # Read angle data
            if in_angles_section and in_molecule:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        atom3 = int(parts[2])
                        angle_val = float(parts[4])
                        force_const = float(parts[5])
                        
                        angles.append((atom1, atom2, atom3))
                        angle_values.append(angle_val)
                        angle_force_constants.append(force_const)
                    except (ValueError, IndexError):
                        continue
    
    return angles, angle_values, angle_force_constants

if __name__ == "__main__":
    # Test the function
    itp_file = "martini_v2.0_lipids_all_201506.itp"
    molecule_name = "DOPC"
    
    print("Testing bond reading...")
    bonds, bond_lengths, bond_force_constants = read_martini_bonds(itp_file, molecule_name)
    
    print("Testing angle reading...")
    angles, angle_values, angle_force_constants = read_martini_angles(itp_file, molecule_name)
    
    print(f"\nBonds for {molecule_name}:")
    for i, (bond, length, force_const) in enumerate(zip(bonds, bond_lengths, bond_force_constants)):
        print(f"  {i+1}: {bond} - length: {length} nm, force_const: {force_const} kJ/mol/nm^2")
    
    print(f"\nAngles for {molecule_name}:")
    for i, (angle, angle_val, force_const) in enumerate(zip(angles, angle_values, angle_force_constants)):
        print(f"  {i+1}: {angle} - angle: {angle_val}°, force_const: {force_const} kJ/mol/deg^2") 