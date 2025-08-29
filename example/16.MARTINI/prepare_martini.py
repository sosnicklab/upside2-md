#!/usr/bin/env python3
"""
MARTINI 3.0 Protein-Lipid System Preparation Script
Prepares simulation files for MARTINI 3.0 force field with protein topology support.
"""

import os
import sys
import numpy as np
import tables as tb
from collections import Counter

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

# --- Helpers: read MARTINI protein topology (ITP) for protein bead typing and connectivity ---
def read_protein_itp_topology(itp_path: str):
    topo_by_res_and_role = {}
    if not os.path.exists(itp_path):
        return topo_by_res_and_role
    in_atoms = False
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                in_atoms = line.lower().startswith('[ atoms')
                continue
            if not in_atoms:
                continue
            # tokens: idx type resnr residue atom cgnr charge ...
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                _idx = int(parts[0])
                bead_type = parts[1]
                resnr = int(parts[2])
                residue = parts[3]
                atom_role = parts[4]  # e.g., BB, SC1, SC2
                # charge is usually the last numeric token; try to parse last column
                charge = 0.0
                for tok in reversed(parts):
                    try:
                        charge = float(tok)
                        break
                    except ValueError:
                        continue
                topo_by_res_and_role[(resnr, atom_role)] = (bead_type, charge)
            except Exception:
                continue
    return topo_by_res_and_role

def read_protein_itp_connectivity(itp_path: str):
    """Read bonds, angles, and dihedrals from protein ITP file"""
    bonds = []
    angles = []
    dihedrals = []
    
    if not os.path.exists(itp_path):
        return bonds, angles, dihedrals
    
    current_section = None
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                section = line.lower()
                if 'bonds' in section:
                    current_section = 'bonds'
                elif 'angles' in section:
                    current_section = 'angles'
                elif 'dihedrals' in section:
                    current_section = 'dihedrals'
                else:
                    current_section = None
                continue
            
            if current_section is None:
                continue
                
            parts = line.split()
            if len(parts) < 3:
                continue
                
            try:
                if current_section == 'bonds':
                    # Format: i j func r0 k
                    i, j = int(parts[0])-1, int(parts[1])-1  # Convert to 0-indexed
                    func = int(parts[2])
                    if func == 1 and len(parts) >= 5:  # Harmonic bond
                        r0 = float(parts[3])  # nm
                        k = float(parts[4])   # kJ/mol/nm²
                        bonds.append((i, j, r0, k))
                elif current_section == 'angles':
                    # Format: i j k func theta0 k
                    i, j, k = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1
                    func = int(parts[3])
                    if func == 2 and len(parts) >= 6:  # Harmonic angle
                        theta0 = float(parts[4])  # degrees
                        k = float(parts[5])       # kJ/mol/rad²
                        angles.append((i, j, k, theta0, k))
                elif current_section == 'dihedrals':
                    # Format: i j k l func phi0 k mult
                    atom_i, atom_j, atom_k, atom_l = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1, int(parts[3])-1
                    func = int(parts[4])
                    if func == 1 and len(parts) >= 8:  # Proper dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const))
                    elif func == 2 and len(parts) >= 7:  # Harmonic dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const))
            except (ValueError, IndexError):
                continue
    
    return bonds, angles, dihedrals

def read_protein_itp_exclusions(itp_path: str):
    """Read exclusions from protein ITP file"""
    exclusions = []
    
    if not os.path.exists(itp_path):
        return exclusions
    
    current_section = None
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                section = line.lower()
                if 'exclusions' in section:
                    current_section = 'exclusions'
                else:
                    current_section = None
                continue
            
            if current_section == 'exclusions':
                parts = line.split()
                if len(parts) >= 2:
                    # Convert to 0-indexed and add all pairs
                    atoms = [int(part)-1 for part in parts]
                    for i in range(len(atoms)):
                        for j in range(i+1, len(atoms)):
                            exclusions.append((atoms[i], atoms[j]))
    
    return exclusions

def main():
    # Get UPSIDE home directory
    upside_path = os.environ['UPSIDE_HOME']
    
    # Configuration
    pdb_id = '1rkl'
    strict_from_martini_pdb = True
    include_protein = True
    
    # Create output directory
    run_dir = f"outputs/martini_test"
    os.makedirs(run_dir, exist_ok=True)
    
    print("=== MARTINI 3.0 Protein-Lipid System Preparation ===")
    print(f"PDB ID: {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    # Read MARTINI parameter files
    print("\n=== Reading MARTINI 3.0 Parameters ===")
    
    # Read nonbonded parameters
    martini_param_file = "ff3.00/martini_v3.0.0.itp"
    martini_table = read_martini3_nonbond_params(martini_param_file)
    
    if not martini_table:
        print("Warning: Could not read MARTINI parameters, using default values")
        martini_table = {}
        # Default fallback parameters
        martini_table[('P2', 'P2')] = (0.47, 4.25)  # Default P2-P2 interaction
    
    # Read DOPC topology
    dopc_param_file = "ff3.00/martini_v3.0.0_phospholipids_v1.itp"
    dopc_bead_types = ['Q1', 'Q5', 'SN4a', 'N4a', 'C1', 'C4h', 'C1', 'C1', 'C1', 'C4h', 'C1', 'C1']
    dopc_charges = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    # Read ion topologies
    ion_param_file = "ff3.00/martini_v3.0.0_ions_v1.itp"
    na_bead_types = ['TQ5']
    na_charges = [1.0]
    cl_bead_types = ['TQ5']
    cl_charges = [-1.0]
    
    # Read water topology
    water_param_file = "ff3.00/martini_v3.0.0_solvents_v1.itp"
    water_bead_types = ['W']
    water_charges = [0.0]
    
    # Read bead masses
    mass_file = "ff3.00/martini_v3.0.0.itp"
    
    # Read DOPC bonds and angles
    dopc_bonds = [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (5, 6), (6, 7), (3, 8), (8, 9), (9, 10), (10, 11)]
    dopc_bond_lengths = [0.40, 0.42, 0.312, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47]  # nm
    dopc_bond_force_constants = [7000.0, 1350.0, 2500.0, 5000.0, 3800.0, 3800.0, 3800.0, 3600.0, 3800.0, 3800.0, 3800.0]  # kJ/mol/nm²
    
    angle_atoms_martini = [(1, 2, 3), (1, 2, 4), (2, 4, 5), (4, 5, 6), (5, 6, 7), (3, 8, 9), (8, 9, 10), (9, 10, 11)]
    angle_equil_deg = [108.0, 139.1, 180.0, 120.0, 180.0, 180.0, 180.0, 180.0]  # degrees
    angle_force_constants = [21.5, 31.2, 35.0, 35.0, 35.0, 35.0, 35.0, 35.0]  # kJ/mol/deg²
    
    # Unit conversions
    energy_conversion = 2.914952774272  # kJ/mol → E_up
    length_conversion = 10.0  # nm → Å
    
    # Bonds: kJ/mol/nm² → E_up/Å²
    # = (kJ/mol → E_up) / (nm² → Å²)
    # = 2.914952774272 / (10.0 * 10.0)
    bond_conversion = energy_conversion / (length_conversion ** 2)  # ≈ 0.02915
    
    # Angles: kJ/mol/deg² → E_up/deg²  
    # = kJ/mol → E_up (degrees stay the same)
    angle_conversion = energy_conversion  # ≈ 2.91
    
    # Dihedrals: kJ/mol → E_up
    dihedral_conversion = energy_conversion  # ≈ 2.91
    
    print("\n=== MARTINI Unit Conversions ===")
    print(f"Bond lengths (nm -> Å): 0.40 nm -> 4.0 Å")
    print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): 7000.0 -> {7000.0 * bond_conversion:.6f}")
    print(f"Angle equilibrium (degrees): 108.0°")
    print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): 21.5 -> {21.5 * angle_conversion:.6f}")
    print(f"Dihedral force constants (kJ/mol -> E_up): 400.0 -> {400.0 * dihedral_conversion:.6f}")
    print(f"Energy conversion factor: {energy_conversion} (kJ/mol -> E_up)")
    print(f"Length conversion: 1 nm = {length_conversion} Å, so 1 nm² = {length_conversion**2} Å²")
    print(f"Bond conversion factor: {bond_conversion:.6f}")
    print(f"Angle conversion factor: {angle_conversion:.6f}")
    print(f"Dihedral conversion factor: {dihedral_conversion:.6f}")
    
    # Read PDB file
    if strict_from_martini_pdb or include_protein:
        input_pdb_file = f'pdb/{pdb_id}.MARTINI.pdb'
        print(f"\nUsing MARTINI PDB as base structure: {input_pdb_file}")
    
    # Read PDB and populate arrays
    print(f"\n=== Reading PDB Structure ===")
    initial_positions = []
    atom_types = []
    charges = []
    residue_ids = []
    atom_names = []
    residue_names = []
    
    # Load protein topology mapping and connectivity if available
    protein_itp = f"pdb/{pdb_id.lower()}_proa.itp"
    protein_topo_map = read_protein_itp_topology(protein_itp)
    protein_bonds, protein_angles, protein_dihedrals = read_protein_itp_connectivity(protein_itp)
    protein_exclusions = read_protein_itp_exclusions(protein_itp)
    
    # Protein bead type aliases
    protein_bead_alias = {
        'BB': 'P2', 'SC1': 'P5', 'SC2': 'SC2', 'SC3': 'SC3', 'SC4': 'SC4',
        'SP1': 'SP1', 'SP2': 'SP2', 'SP5': 'SP5', 'TP1': 'TP1', 'TP2': 'TP2',
        'TC3': 'TC3', 'TC4': 'TC4', 'TC5': 'TC5', 'TN5a': 'TN5a', 'TN6': 'TN6',
        'TN6d': 'TN6d', 'SQ5n': 'SQ5n', 'Q5n': 'Q5n', 'C6': 'C6'
    }
    
    # Standard MARTINI mappings
    pdb_to_martini = {
        'NC3': 'Q1', 'PO4': 'Q5', 'GL1': 'SN4a', 'GL2': 'N4a',
        'C1A': 'C1', 'C2A': 'C4h', 'C3A': 'C1', 'C4A': 'C1',
        'C1B': 'C1', 'C2B': 'C4h', 'C3B': 'C1', 'C4B': 'C1'
    }
    
    martini_charges = {
        'Q1': 1.0, 'Q5': -1.0, 'SN4a': 0.0, 'N4a': 0.0,
        'C1': 0.0, 'C4h': 0.0, 'W': 0.0, 'TQ5': 1.0
    }
    
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
                
            # Parse PDB line
            atom_name_raw = line[12:16].strip()
            atom_name = atom_name_raw.upper()
            atom_names.append(atom_name)
            residue_id = int(line[22:26])
            residue_ids.append(residue_id)
            residue_name = line[17:20].strip().upper()
            residue_names.append(residue_name)
            
            
            
            # Extract coordinates
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            
            # Determine if this line corresponds to protein beads
            is_protein = ('PROA' in line)
            
            # Map to MARTINI type based on context
            if is_protein:
                # Prefer exact topology mapping by (resnr, role) when available
                role = atom_name  # BB / SC1 / SC2 ... expected in our MARTINI PDB
                if (residue_id, role) in protein_topo_map:
                    martini_type, charge = protein_topo_map[(residue_id, role)]
                else:
                    # Fallback to alias to standard MARTINI type
                    martini_type = protein_bead_alias.get(atom_name, 'P2')
                    charge = 0.0
            elif residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file
                if atom_name in dopc_bead_types:
                    atom_idx = dopc_bead_types.index(atom_name)
                    martini_type = dopc_bead_types[atom_idx]
                    charge = dopc_charges[atom_idx]
                else:
                    martini_type = pdb_to_martini.get(atom_name, atom_name)
                    charge = martini_charges.get(martini_type, 0.0)
            elif residue_name == 'W':
                martini_type = water_bead_types[0] if water_bead_types else 'W'
                charge = water_charges[0] if water_charges else 0.0
            elif residue_name == 'NA':
                martini_type = na_bead_types[0] if na_bead_types else 'TQ5'
                charge = na_charges[0] if na_charges else 1.0
            elif residue_name == 'CL':
                martini_type = cl_bead_types[0] if cl_bead_types else 'TQ5'
                charge = cl_charges[0] if cl_charges else -1.0
            else:
                martini_type = protein_bead_alias.get(atom_name, pdb_to_martini.get(atom_name, atom_name))
                charge = martini_charges.get(martini_type, 0.0)
            
            atom_types.append(martini_type)
            charges.append(charge)
    
    # Convert to numpy arrays
    initial_positions = np.array(initial_positions, dtype=float)
    atom_types = np.array(atom_types)
    charges = np.array(charges, dtype=float)
    residue_ids = np.array(residue_ids, dtype=int)
    atom_names = np.array(atom_names)
    residue_names = np.array(residue_names)
    n_atoms = len(initial_positions)
    
    # Read box dimensions from CRYST1 record
    print(f"Reading box dimensions from {input_pdb_file}...")
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                fields = line.split()
                if len(fields) >= 4:
                    x_len = float(fields[1])
                    y_len = float(fields[2])
                    z_len = float(fields[3])
                    print(f"Found CRYST1 record: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
                    break
        else:
            x_len = 50.0
            y_len = 50.0
            z_len = 50.0
            print(f"WARNING: No CRYST1 record found in PDB file!")
            print(f"Using default box dimensions: X={x_len:.1f}, Y={y_len:.1f}, Z={z_len:.1f} Angstroms")
    
    # Print system parameters
    print(f"\n=== System Parameters ===")
    print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
    print(f"Box volume: {x_len * y_len * z_len:.1f} Å³")
    print(f"Total atoms: {n_atoms}")
    
    # Group atoms into molecules
    molecules = []
    current_mol_atoms = []
    current_mol_names = []
    current_mol_indices = []
    current_resid = None
    current_resname = None
    
    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        if resid != current_resid:
            if current_mol_atoms:
                molecules.append((current_resname, current_mol_atoms, current_mol_indices))
            current_mol_atoms = [atom_name]
            current_mol_names = [atom_name]
            current_mol_indices = [i]
            current_resid = resid
            current_resname = resname
        else:
            current_mol_atoms.append(atom_name)
            current_mol_indices.append(i)
    
    if current_mol_atoms:
        molecules.append((current_resname, current_mol_atoms, current_mol_indices))
    
    # Count molecules
    mol_counts = Counter(mol[0] for mol in molecules)
    dopc_count = mol_counts.get('DOP', 0)  # DOP is the residue name, not DOPC
    water_count = mol_counts.get('W', 0)
    
    print(f"\n=== Molecule Summary ===")
    for moltype, count in mol_counts.items():
        print(f"{moltype}: {count} molecules")
    

    
    # Create bonds and angles
    print(f"\n=== Creating Connectivity ===")
    
    # Initialize lists for bonds and angles
    bonds_list = []
    bond_lengths_list = []
    bond_force_constants_list = []
    angles_list = []
    angle_equil_deg_list = []
    angle_force_constants_list = []
    dihedrals_list = []
    dihedral_equil_deg_list = []
    dihedral_force_constants_list = []
    
    # Create DOPC bonds and angles
    dopc_molecules = [mol for mol in molecules if mol[0] == 'DOP']  # DOP is the residue name, not DOPC
    
    for mol_idx, (_, atom_names_mol, atom_indices) in enumerate(dopc_molecules):
        name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}
        
        # Create bonds for this lipid
        for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
            if bond_idx1 < len(dopc_bead_types) and bond_idx2 < len(dopc_bead_types):
                atom1_name = atom_names_mol[bond_idx1]
                atom2_name = atom_names_mol[bond_idx2]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                bonds_list.append([atom1, atom2])
                bond_lengths_list.append(dopc_bond_lengths[i] * 10.0)  # nm to Å
                bond_force_constants_list.append(dopc_bond_force_constants[i] * bond_conversion)  # kJ/mol/nm² to E_up/Å²
        
        # Create angles for this lipid
        for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(angle_atoms_martini):
            if (angle_idx1 < len(dopc_bead_types) and angle_idx2 < len(dopc_bead_types) and 
                angle_idx3 < len(dopc_bead_types)):
                atom1_name = atom_names_mol[angle_idx1]
                atom2_name = atom_names_mol[angle_idx2]
                atom3_name = atom_names_mol[angle_idx3]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                atom3 = name_to_idx[atom3_name]
                angles_list.append([atom1, atom2, atom3])
                angle_equil_deg_list.append(angle_equil_deg[i])
                angle_force_constants_list.append(angle_force_constants[i] * angle_conversion)  # kJ/mol/deg² to E_up/deg²
    
    print(f"Created {len(bonds_list)} bonds for {dopc_count} DOPC lipids")
    print(f"Created {len(angles_list)} angles for {dopc_count} DOPC lipids")
    
    # Create protein connectivity if available
    protein_bond_count = 0
    protein_angle_count = 0
    protein_dihedral_count = 0
    
    if protein_bonds:
        print(f"\n=== Protein Connectivity from {protein_itp} ===")
        print(f"Found {len(protein_bonds)} bonds, {len(protein_angles)} angles, {len(protein_dihedrals)} dihedrals")
        
        # Add protein bonds to the bond list
        for i, j, r0_nm, k_kj in protein_bonds:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list (assuming protein atoms come first in the system)
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_bond_count += 1
        
        # Add protein angles to the angle list
        for i, j, k, theta0_deg, k_kj in protein_angles:
            # Convert MARTINI units to UPSIDE units
            theta0_upside = theta0_deg  # degrees (same unit)
            k_upside = k_kj * angle_conversion  # kJ/mol/deg² to E_up/deg²
            
            # Add to angle list
            angles_list.append([i, j, k])
            angle_equil_deg_list.append(theta0_upside)
            angle_force_constants_list.append(k_upside)
            protein_angle_count += 1
        
        # Add protein dihedrals to the dihedral list
        for i, j, k, l, phi0_deg, k_kj in protein_dihedrals:
            # Convert MARTINI units to UPSIDE units
            phi0_upside = phi0_deg  # degrees (same unit)
            k_upside = k_kj * dihedral_conversion  # kJ/mol to E_up
            
            # Add to dihedral list
            dihedrals_list.append([i, j, k, l])
            dihedral_equil_deg_list.append(phi0_upside)
            dihedral_force_constants_list.append(k_upside)
            protein_dihedral_count += 1
        
        print(f"Added {protein_bond_count} protein bonds")
        print(f"Added {protein_angle_count} protein angles")
        print(f"Added {protein_dihedral_count} protein dihedrals")
    else:
        print(f"\nNo protein connectivity found in {protein_itp}")
    
    print(f"Total system bonds: {len(bonds_list)}")
    print(f"Total system angles: {len(angles_list)}")
    print(f"Total system dihedrals: {len(dihedrals_list)}")
    
    # Center and wrap positions
    print(f"\n=== Preparing Final Structure ===")
    center = np.mean(initial_positions, axis=0)
    centered_positions = initial_positions - center
    half_box = np.array([x_len/2, y_len/2, z_len/2])
    centered_positions = (centered_positions + half_box) % (2*half_box) - half_box
    final_positions = centered_positions
    
    print("First 10 positions (centered and wrapped):")
    for i in range(min(10, n_atoms)):
        print(f"  Atom {i}: {final_positions[i]}")
    
    # Create UPSIDE input file
    print(f"\n=== Creating UPSIDE Input File ===")
    input_file = f"{run_dir}/test.input.up"
    
    with tb.open_file(input_file, 'w') as t:
        # Create input group (required by UPSIDE)
        input_grp = t.create_group(t.root, 'input')
        
        # Create position array with correct format
        pos = np.zeros((n_atoms, 3, 1), dtype='f4')
        pos[:,:,0] = final_positions
        pos_array = t.create_array(input_grp, 'pos', obj=pos)
        pos_array._v_attrs.arguments = np.array([b'pos'])
        pos_array._v_attrs.shape = pos.shape
        pos_array._v_attrs.n_atoms = n_atoms
        pos_array._v_attrs.n_frames = 1
        pos_array._v_attrs.dim = 3
        pos_array._v_attrs.initialized = True
        
        # Create velocity array (required by UPSIDE)
        velocity = np.zeros((n_atoms, 3), dtype='f4')
        vel_array = t.create_array(input_grp, 'vel', obj=velocity)
        vel_array._v_attrs.arguments = np.array([b'vel'])
        vel_array._v_attrs.shape = velocity.shape
        vel_array._v_attrs.n_atoms = n_atoms
        vel_array._v_attrs.dim = 3
        vel_array._v_attrs.initialized = True
        
        # Create momentum array (required by UPSIDE)
        momentum = np.zeros((n_atoms, 3, 1), dtype='f4')
        mom_array = t.create_array(input_grp, 'mom', obj=momentum)
        mom_array._v_attrs.arguments = np.array([b'mom'])
        mom_array._v_attrs.shape = momentum.shape
        mom_array._v_attrs.n_atoms = n_atoms
        mom_array._v_attrs.dim = 3
        mom_array._v_attrs.initialized = True
        
        # Create mass array (required by UPSIDE)
        mass = np.ones(n_atoms, dtype='f4')  # Default mass of 1.0 (normalized)
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        
        # Create type array
        type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
        type_array._v_attrs.arguments = np.array([b'type'])
        type_array._v_attrs.shape = atom_types.shape
        type_array._v_attrs.n_atoms = n_atoms
        type_array._v_attrs.initialized = True
        
        # Create charges array
        charge_array = t.create_array(input_grp, 'charges', obj=charges)
        charge_array._v_attrs.arguments = np.array([b'charges'])
        charge_array._v_attrs.shape = charges.shape
        charge_array._v_attrs.n_atoms = n_atoms
        charge_array._v_attrs.initialized = True
        
        # Create residue IDs array
        residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
        residue_array._v_attrs.arguments = np.array([b'residue_ids'])
        residue_array._v_attrs.shape = residue_ids.shape
        residue_array._v_attrs.n_atoms = n_atoms
        residue_array._v_attrs.initialized = True
        
        # Create potential group (required by UPSIDE)
        potential_grp = t.create_group(input_grp, 'potential')
        
        # Create MARTINI potential with proper parameters
        martini_potential = t.create_group(potential_grp, 'martini_potential')
        martini_potential._v_attrs.arguments = np.array([b'pos'])
        martini_potential._v_attrs.potential_type = b'lj_coulomb'
        martini_potential._v_attrs.epsilon = 4.7  # Default MARTINI epsilon
        martini_potential._v_attrs.sigma = 4.7    # Default MARTINI sigma
        martini_potential._v_attrs.lj_cutoff = 12.0
        martini_potential._v_attrs.coul_cutoff = 12.0
        martini_potential._v_attrs.dielectric = 15.0
        martini_potential._v_attrs.coulomb_constant = 476.627809876965
        martini_potential._v_attrs.n_types = 1
        martini_potential._v_attrs.n_params = 4
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 1.0
        martini_potential._v_attrs.initialized = True
        martini_potential._v_attrs.force_cap = 0
        martini_potential._v_attrs.mass_scale = 1.0 / 72.0
        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        martini_potential._v_attrs.debug_mode = 1  # Enable spline table generation
        
        # Create periodic boundary potential
        wall_group = t.create_group(potential_grp, 'periodic_boundary_potential')
        wall_group._v_attrs.arguments = np.array([b'pos'])
        wall_group._v_attrs.x_len = x_len
        wall_group._v_attrs.y_len = y_len
        wall_group._v_attrs.z_len = z_len
        wall_group._v_attrs.initialized = True
        
        # Create atom indices and charges arrays for the potential
        t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
        t.create_array(martini_potential, 'charges', obj=charges)
        
        # Create pairs and coefficients for non-bonded interactions with proper exclusions
        pairs_list = []
        coeff_array = []
        
        # Create sets for different exclusion levels
        bonded_pairs_12 = set()  # Directly bonded (1-2) - full exclusion
        bonded_pairs_13 = set()  # Connected by angles (1-3) - full exclusion
        bonded_pairs_14 = set()  # Connected by dihedrals (1-4) - scaled interaction
        additional_exclusions = set()  # Additional exclusions from ITP file
        
        # Add 1-2 exclusions from bond list
        for bond in bonds_list:
            sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
            bonded_pairs_12.add(sorted_bond)
        
        # Add 1-3 exclusions from angle list
        for angle in angles_list:
            # Atoms 0 and 2 in angle are 1-3 connected
            sorted_pair = (min(angle[0], angle[2]), max(angle[0], angle[2]))
            bonded_pairs_13.add(sorted_pair)
        
        # Add 1-4 pairs from dihedral list
        for dihedral in dihedrals_list:
            # Atoms 0 and 3 in dihedral are 1-4 connected
            sorted_pair = (min(dihedral[0], dihedral[3]), max(dihedral[0], dihedral[3]))
            bonded_pairs_14.add(sorted_pair)
        
        # Add additional exclusions from protein ITP file
        if protein_exclusions:
            for exclusion in protein_exclusions:
                sorted_exclusion = (min(exclusion[0], exclusion[1]), max(exclusion[0], exclusion[1]))
                additional_exclusions.add(sorted_exclusion)
        
        # Generate all unique pairs (i < j) with proper exclusions
        excluded_12_count = 0
        excluded_13_count = 0
        excluded_additional_count = 0
        scaled_14_count = 0
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                pair = (i, j)
                
                # Skip 1-2 pairs (full exclusion)
                if pair in bonded_pairs_12:
                    excluded_12_count += 1
                    continue
                
                # Skip 1-3 pairs (full exclusion)
                if pair in bonded_pairs_13:
                    excluded_13_count += 1
                    continue
                
                # Skip additional exclusions from ITP file
                if pair in additional_exclusions:
                    excluded_additional_count += 1
                    continue
                
                # Scale 1-4 pairs (typically 0.5 scaling factor)
                scale_factor = 1.0
                if pair in bonded_pairs_14:
                    scale_factor = 0.5  # Standard 1-4 scaling
                    scaled_14_count += 1
                
                pairs_list.append([i, j])
                
                # Get bead types for this pair
                type1 = atom_types[i].decode('utf-8') if isinstance(atom_types[i], bytes) else str(atom_types[i])
                type2 = atom_types[j].decode('utf-8') if isinstance(atom_types[j], bytes) else str(atom_types[j])
                
                # Look up MARTINI parameters for this bead type pair
                if (type1, type2) in martini_table:
                    sigma_nm, epsilon_kj = martini_table[(type1, type2)]
                elif (type2, type1) in martini_table:
                    sigma_nm, epsilon_kj = martini_table[(type2, type1)]
                else:
                    # Fallback to default parameters if not found
                    sigma_nm = 0.47  # nm
                    epsilon_kj = 4.25  # kJ/mol
                
                # Convert to UPSIDE units
                epsilon = epsilon_kj / energy_conversion  # kJ/mol → E_up
                sigma = sigma_nm * length_conversion  # nm → Å
                q1 = charges[i] * scale_factor
                q2 = charges[j] * scale_factor
                coeff_array.append([epsilon * scale_factor, sigma, q1, q2])
        
        print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions")
        print(f"Excluded {excluded_13_count} 1-3 angle-connected pairs from non-bonded interactions")
        print(f"Excluded {excluded_additional_count} additional pairs from ITP exclusions")
        print(f"Scaled {scaled_14_count} 1-4 dihedral-connected pairs with factor 0.5")

        t.create_array(martini_potential, 'pairs', obj=np.array(pairs_list, dtype=int))
        t.create_array(martini_potential, 'coefficients', obj=np.array(coeff_array, dtype='f4'))

        # Add bonded potentials mirroring original run_martini.py
        # Bonds: dist_spring
        if bonds_list:
            bond_group = t.create_group(potential_grp, 'dist_spring')
            bond_group._v_attrs.arguments = np.array([b'pos'])
            bond_group._v_attrs.initialized = True
            bond_group._v_attrs.mass_scale = 1.0 / 72.0
            bond_group._v_attrs.x_len = x_len
            bond_group._v_attrs.y_len = y_len
            bond_group._v_attrs.z_len = z_len
            bond_group._v_attrs.debug_mode = 1  # Enable spline table generation

            t.create_array(bond_group, 'id', obj=np.array(bonds_list, dtype=int))
            t.create_array(bond_group, 'equil_dist', obj=np.array(bond_lengths_list, dtype='f4'))
            t.create_array(bond_group, 'spring_const', obj=np.array(bond_force_constants_list, dtype='f4'))
            # Compatibility dataset expected by some builds
            t.create_array(bond_group, 'bonded_atoms', obj=np.ones(len(bonds_list), dtype='i4'))

        # Angles: angle_spring
        if angles_list:
            angle_group = t.create_group(potential_grp, 'angle_spring')
            angle_group._v_attrs.arguments = np.array([b'pos'])
            angle_group._v_attrs.initialized = True
            angle_group._v_attrs.mass_scale = 1.0 / 72.0
            angle_group._v_attrs.x_len = x_len
            angle_group._v_attrs.y_len = y_len
            angle_group._v_attrs.z_len = z_len
            angle_group._v_attrs.debug_mode = 1  # Enable spline table generation

            t.create_array(angle_group, 'id', obj=np.array(angles_list, dtype=int))
            t.create_array(angle_group, 'equil_angle_deg', obj=np.array(angle_equil_deg_list, dtype='f4'))
            t.create_array(angle_group, 'spring_const', obj=np.array(angle_force_constants_list, dtype='f4'))

        # Dihedrals: dihedral_spring (if present)
        if dihedrals_list:
            dihedral_group = t.create_group(potential_grp, 'dihedral_spring')
            dihedral_group._v_attrs.arguments = np.array([b'pos'])
            dihedral_group._v_attrs.initialized = True
            dihedral_group._v_attrs.mass_scale = 1.0 / 72.0
            dihedral_group._v_attrs.x_len = x_len
            dihedral_group._v_attrs.y_len = y_len
            dihedral_group._v_attrs.z_len = z_len
            dihedral_group._v_attrs.debug_mode = 1  # Enable spline table generation

            t.create_array(dihedral_group, 'id', obj=np.array(dihedrals_list, dtype=int))
            # Some builds of UPSIDE expect 'equil_dist' for this potential; write both names for compatibility
            eq_deg = np.array(dihedral_equil_deg_list, dtype='f4')
            t.create_array(dihedral_group, 'equil_angle_deg', obj=eq_deg)
            t.create_array(dihedral_group, 'equil_dist', obj=eq_deg)
            t.create_array(dihedral_group, 'spring_const', obj=np.array(dihedral_force_constants_list, dtype='f4'))
    
    print(f"Created UPSIDE input file: {input_file}")
    print(f"Preparation complete!")
    
    # Save preparation summary
    summary_file = f"{run_dir}/preparation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=== MARTINI 3.0 Preparation Summary ===\n")
        f.write(f"PDB ID: {pdb_id}\n")
        f.write(f"Total atoms: {n_atoms}\n")
        f.write(f"Box dimensions: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å\n")
        f.write(f"Total bonds: {len(bonds_list)}\n")
        f.write(f"Total angles: {len(angles_list)}\n")
        f.write(f"Total dihedrals: {len(dihedrals_list)}\n")
        f.write(f"Protein bonds: {protein_bond_count}\n")
        f.write(f"Protein angles: {protein_angle_count}\n")
        f.write(f"Protein dihedrals: {protein_dihedral_count}\n")
        f.write(f"DOPC lipids: {dopc_count}\n")
        f.write(f"Water molecules: {water_count}\n")
        f.write(f"1-2 exclusions: {excluded_12_count}\n")
        f.write(f"1-3 exclusions: {excluded_13_count}\n")
        f.write(f"Additional exclusions: {excluded_additional_count}\n")
        f.write(f"1-4 scaled pairs: {scaled_14_count}\n")
    
    print(f"Preparation summary saved to: {summary_file}")

if __name__ == "__main__":
    main()
