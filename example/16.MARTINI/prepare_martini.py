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

def read_protein_itp_connectivity(itp_path: str, simulation_stage='minimization'):
    """
    Read bonds, angles, dihedrals, constraints, and position restraints from protein ITP file
    simulation_stage: 'minimization' or 'production'
    """
    bonds = []
    angles = []
    dihedrals = []
    constraints = []
    position_restraints = []
    
    if not os.path.exists(itp_path):
        return bonds, angles, dihedrals, constraints, position_restraints
    
    current_section = None
    in_normang_section = False
    in_nonnormang_section = False
    in_flexible_section = False
    in_posres_section = False
    
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            
            # Handle preprocessor directives
            if line.startswith('#ifdef NORMANG'):
                in_normang_section = True
                in_nonnormang_section = False
                continue
            elif line.startswith('#ifndef NORMANG'):
                in_normang_section = False
                in_nonnormang_section = True
                continue
            elif line.startswith('#ifdef FLEXIBLE'):
                in_flexible_section = True
                continue
            elif line.startswith('#ifndef FLEXIBLE'):
                in_flexible_section = False
                continue
            elif line.startswith('#ifdef POSRES'):
                in_posres_section = True
                continue
            elif line.startswith('#ifndef POSRES'):
                # Don't reset in_posres_section here - let #endif handle it
                continue
            elif line.startswith('#ifndef POSRES_FC'):
                # This is a different directive - don't reset in_posres_section
                continue
            elif line.startswith('#endif'):
                # Only reset flags if we're not in a position restraints section
                if current_section != 'position_restraints':
                    in_normang_section = False
                    in_nonnormang_section = False
                    in_flexible_section = False
                    in_posres_section = False
                continue
            
            if line.startswith('['):
                section = line.lower()
                if 'bonds' in section:
                    current_section = 'bonds'
                elif 'angles' in section:
                    current_section = 'angles'
                elif 'dihedrals' in section:
                    current_section = 'dihedrals'
                elif 'constraints' in section:
                    current_section = 'constraints'
                elif 'position_restraints' in section:
                    current_section = 'position_restraints'
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
                        
                        # For minimization: use FLEXIBLE bonds (large spring constants)
                        # For production: use regular bonds
                        if simulation_stage == 'minimization' and in_flexible_section:
                            k = 1000000.0  # Large spring constant for minimization
                        
                        bonds.append((i, j, r0, k))
                elif current_section == 'angles':
                    # Format: i j k func theta0 k
                    # For minimization: use NORMANG section
                    # For production: use regular angles section
                    if simulation_stage == 'minimization' and in_normang_section and not in_nonnormang_section:
                        i, j, k = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1
                        func = int(parts[3])
                        if func == 2 and len(parts) >= 6:  # Harmonic angle
                            theta0 = float(parts[4])  # degrees
                            k = float(parts[5])       # kJ/mol/rad²
                            angles.append((i, j, k, theta0, k))
                    elif simulation_stage == 'production' and not in_normang_section and not in_nonnormang_section:
                        # Regular angles section for production
                        i, j, k = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1
                        func = int(parts[3])
                        if func == 2 and len(parts) >= 6:  # Harmonic angle
                            theta0 = float(parts[4])  # degrees
                            k = float(parts[5])       # kJ/mol/rad²
                            angles.append((i, j, k, theta0, k))
                elif current_section == 'constraints':
                    # Format: i j func r0
                    # These should be treated as bonds with large spring constants
                    i, j = int(parts[0])-1, int(parts[1])-1  # Convert to 0-indexed
                    func = int(parts[2])
                    if func == 1 and len(parts) >= 4:  # Constraint
                        r0 = float(parts[3])  # nm
                        k = 1000000.0  # Large spring constant as requested
                        constraints.append((i, j, r0, k))
                elif current_section == 'position_restraints':
                    # Format: i func fx fy fz
                    # Only process during minimization stage
                    if simulation_stage == 'minimization' and in_posres_section:
                        i = int(parts[0])-1  # Convert to 0-indexed
                        func = int(parts[1])
                        if func == 1 and len(parts) >= 5:  # Position restraint
                            fx = float(parts[2]) if parts[2] != 'POSRES_FC' else 1000.0
                            fy = float(parts[3]) if parts[3] != 'POSRES_FC' else 1000.0
                            fz = float(parts[4]) if parts[4] != 'POSRES_FC' else 1000.0
                            position_restraints.append((i, fx, fy, fz))
                elif current_section == 'dihedrals':
                    # Format: i j k l func phi0 k mult
                    atom_i, atom_j, atom_k, atom_l = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1, int(parts[3])-1
                    func = int(parts[4])
                    if func == 1 and len(parts) >= 8:  # Periodic dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const, func))
                    elif func == 2 and len(parts) >= 7:  # Harmonic dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const, func))
            except (ValueError, IndexError):
                continue
    
    return bonds, angles, dihedrals, constraints, position_restraints

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

def parse_itp_file(itp_file, target_molecule=None):
    """
    Universal ITP parser to read MARTINI topology files.
    Returns a dictionary with parsed sections: atoms, bonds, angles, dihedrals, etc.
    If target_molecule is specified, only returns data for that specific molecule type.
    """
    topology = {
        'atoms': [],
        'bonds': [], 
        'angles': [],
        'dihedrals': [],
        'exclusions': [],
        'moleculetype': None,
        'molecules': {}  # Store multiple molecule types
    }
    
    if not os.path.exists(itp_file):
        print(f"Warning: ITP file {itp_file} not found")
        return topology
    
    current_section = None
    current_molecule = None
    current_mol_data = None
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith(';'):
                continue
            
            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                section_name = line[1:-1].strip().lower()
                current_section = section_name
                continue
            
            # Parse based on current section
            if current_section == 'moleculetype':
                parts = line.split()
                if len(parts) >= 1:
                    current_molecule = parts[0]
                    if topology['moleculetype'] is None:
                        topology['moleculetype'] = current_molecule
                    
                    # Initialize data structure for this molecule
                    current_mol_data = {
                        'atoms': [],
                        'bonds': [],
                        'angles': [],
                        'dihedrals': [],
                        'exclusions': []
                    }
                    topology['molecules'][current_molecule] = current_mol_data
            
            elif current_section == 'atoms' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        atom_data = {
                            'id': int(parts[0]),
                            'type': parts[1],
                            'resnr': int(parts[2]),
                            'residue': parts[3],
                            'atom': parts[4],
                            'cgnr': int(parts[5]),
                            'charge': 0.0,
                            'mass': 0.0
                        }
                        # Parse charge and mass more carefully
                        # Charge is typically in column 6 (0-indexed), mass in column 7
                        if len(parts) >= 7:
                            try:
                                atom_data['charge'] = float(parts[6])
                            except ValueError:
                                pass
                        if len(parts) >= 8:
                            try:
                                atom_data['mass'] = float(parts[7])
                            except ValueError:
                                pass
                        current_mol_data['atoms'].append(atom_data)
                        topology['atoms'].append(atom_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'bonds' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        bond_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'func': int(parts[2]),
                            'r0': 0.0,
                            'k': 0.0
                        }
                        if len(parts) >= 5:
                            bond_data['r0'] = float(parts[3])  # nm
                            bond_data['k'] = float(parts[4])   # kJ/mol/nm²
                        current_mol_data['bonds'].append(bond_data)
                        topology['bonds'].append(bond_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'angles' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 6:  # Need at least i j k func theta0 force_k
                    try:
                        angle_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'func': int(parts[3]),
                            'theta0': float(parts[4]),  # degrees
                            'force_k': float(parts[5])  # kJ/mol/rad² (renamed to avoid conflict)
                        }
                        current_mol_data['angles'].append(angle_data)
                        topology['angles'].append(angle_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'dihedrals' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        dihedral_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'l': int(parts[3]) - 1,
                            'func': int(parts[4]),
                            'phi0': 0.0,
                            'k': 0.0,
                            'mult': 1
                        }
                        if len(parts) >= 7:
                            dihedral_data['phi0'] = float(parts[5])  # degrees
                            dihedral_data['k'] = float(parts[6])     # kJ/mol
                        if len(parts) >= 8:
                            dihedral_data['mult'] = int(parts[7])
                        current_mol_data['dihedrals'].append(dihedral_data)
                        topology['dihedrals'].append(dihedral_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'exclusions' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # Convert to 0-indexed and add all pairs
                        atoms = [int(part)-1 for part in parts]
                        for i in range(len(atoms)):
                            for j in range(i+1, len(atoms)):
                                exclusion = (atoms[i], atoms[j])
                                current_mol_data['exclusions'].append(exclusion)
                                topology['exclusions'].append(exclusion)
                    except ValueError:
                        continue
    
    # If target_molecule is specified, return only that molecule's data
    if target_molecule and target_molecule in topology['molecules']:
        mol_data = topology['molecules'][target_molecule]
        return {
            'atoms': mol_data['atoms'],
            'bonds': mol_data['bonds'],
            'angles': mol_data['angles'],
            'dihedrals': mol_data['dihedrals'],
            'exclusions': mol_data['exclusions'],
            'moleculetype': target_molecule,
            'molecules': {target_molecule: mol_data}
        }
    
    return topology

def read_martini_masses(ff_file):
    """Read atom type masses from MARTINI force field file"""
    masses = {}
    if not os.path.exists(ff_file):
        raise ValueError(f"FATAL ERROR: Force field file '{ff_file}' not found.\n"
                        f"  This file is required for atom type masses.\n"
                        f"  Please ensure the MARTINI force field file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    in_atomtypes = False
    with open(ff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[ atomtypes ]'):
                in_atomtypes = True
                continue
            elif line.startswith('['):
                in_atomtypes = False
                continue
            
            if in_atomtypes and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    atom_type = parts[0]
                    try:
                        mass = float(parts[1])
                        masses[atom_type] = mass
                    except ValueError:
                        continue
    
    return masses

def main(stage='minimization', run_dir=None):
    """
    Main preparation function with stage-specific parameterization
    
    Args:
        stage: Simulation stage ('minimization', 'npt_equil', 'npt_equil_reduced', 'npt_prod')
        run_dir: Optional run directory (default: outputs/martini_test)
    """

    # Get UPSIDE home directory
    upside_path = os.environ['UPSIDE_HOME']

    # Check command line arguments
    if len(sys.argv) > 1:
        pdb_id = sys.argv[1]
    else:
        raise ValueError("FATAL ERROR: No PDB ID provided as command line argument.\n"
                        f"  Usage: python {sys.argv[0]} <pdb_id>\n"
                        f"  Example: python {sys.argv[0]} 1rkl\n"
                        f"  Aborting to prevent incorrect simulation results.")

    # Get run directory from function parameter or use default
    if run_dir is None:
        if len(sys.argv) > 2 and sys.argv[2] != '--stage':
            run_dir = sys.argv[2]
        else:
            run_dir = "outputs/martini_test"
    os.makedirs(run_dir, exist_ok=True)

    # Get stage from environment variable or use default
    stage = os.environ.get('UPSIDE_SIMULATION_STAGE', stage)
    print(f"Preparing for stage: {stage}")

    # Stage-specific parameterization
    stage_params = {
        'minimization': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
            'barostat_type': 0  # Berendsen
        },
        'npt_equil': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
            'barostat_type': 0  # Berendsen
        },
        'npt_equil_reduced': {
            'lj_soften': 1,
            'lj_alpha': 0.05,
            'coulomb_soften': 1,
            'slater_alpha': 0.5,
            'barostat_type': 0  # Berendsen
        },
        'npt_prod': {
            'lj_soften': 0,
            'lj_alpha': 0.0,
            'coulomb_soften': 0,
            'slater_alpha': 0.0,
            'barostat_type': 1  # Parrinello-Rahman
        }
    }

    # Get parameters for current stage
    params = stage_params.get(stage, stage_params['npt_prod'])


    # Configuration
    strict_from_martini_pdb = True
    include_protein = True
    
    print("=== MARTINI 3.0 Protein-Lipid System Preparation ===")
    print(f"PDB ID: {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    # Read MARTINI parameter files
    print("\n=== Reading MARTINI 3.0 Parameters ===")
    
    # Read nonbonded parameters
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    martini_param_file = os.path.join(SCRIPT_DIR, "ff3.00/martini_v3.0.0.itp")
    martini_table = read_martini3_nonbond_params(martini_param_file)
    
    if not martini_table:
        raise ValueError(f"FATAL ERROR: Could not read MARTINI parameters from '{martini_param_file}'\n"
                        f"  This file is required for proper force field parameterization.\n"
                        f"  Please ensure the MARTINI 3.0 parameter file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    # Determine system type and load appropriate topology
    # Check if this is a protein system by looking for protein ITP file
    protein_itp = os.path.join(SCRIPT_DIR, f"pdb/{pdb_id}_proa.itp")
    has_protein = os.path.exists(protein_itp)
    
    # Check if this is a mixed protein-lipid system by looking for both protein and lipid residues
    # This will be determined during PDB parsing, but we can prepare for both cases
    
    if has_protein:
        print("=== Mixed Protein-Lipid System Detected ===")
        print(f"Using protein topology from: {protein_itp}")
        print("System contains both protein and lipid components")
        print("Using MARTINI 3.0 force field parameters for both protein and lipid components")
    else:
        print("=== Lipid System Detected ===")
    
    # For both protein and lipid systems, we need DOPC parameters
    # Parse DOPC topology from ITP file
    dopc_param_file = os.path.join(SCRIPT_DIR, "ff3.00/martini_v3.0.0_phospholipids_v1.itp")
    # First check what molecules are available
    full_topology = parse_itp_file(dopc_param_file)
    print(f"Available molecules in {dopc_param_file}: {list(full_topology['molecules'].keys())}")
    
    # Try to find DOPC or similar molecule
    dopc_molecule = None
    for mol_name in full_topology['molecules'].keys():
        if 'DOPC' in mol_name.upper() or 'DOP' in mol_name.upper():
            dopc_molecule = mol_name
            break
    
    if dopc_molecule:
        dopc_topology = parse_itp_file(dopc_param_file, dopc_molecule)
        print(f"Using molecule type: {dopc_molecule}")
    else:
        available_molecules = list(full_topology['molecules'].keys())
        raise ValueError(f"FATAL ERROR: DOPC molecule not found in '{dopc_param_file}'.\n"
                        f"  Available molecules: {available_molecules}\n"
                        f"  Please ensure DOPC is defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    dopc_bead_types = [atom['type'] for atom in dopc_topology['atoms']]
    dopc_charges = [atom['charge'] for atom in dopc_topology['atoms']]
    # Create mapping from atom names to bead types and charges
    dopc_atom_to_type = {atom['atom']: atom['type'] for atom in dopc_topology['atoms']}
    dopc_atom_to_charge = {atom['atom']: atom['charge'] for atom in dopc_topology['atoms']}
    print(f"Read DOPC topology: {len(dopc_bead_types)} bead types from {dopc_param_file}")
    print(f"DOPC atom name mapping: {dopc_atom_to_type}")
    
    # Parse ion topologies from ITP file
    ion_param_file = os.path.join(SCRIPT_DIR, "ff3.00/martini_v3.0.0_ions_v1.itp")
    ion_topology = parse_itp_file(ion_param_file)
    # Extract NA and CL atoms specifically
    # In MARTINI ion ITP files, residue name is "ION" and atom name is "NA" or "CL"
    na_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'NA']
    cl_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'CL']
    
    # For standard ions, use the first occurrence (standard chloride is TQ5, not SQ5n which is for acetate)
    na_bead_types = [na_atoms[0]['type']] if na_atoms else []
    na_charges = [na_atoms[0]['charge']] if na_atoms else []
    cl_bead_types = [cl_atoms[0]['type']] if cl_atoms else []  # Use first CL (TQ5), not acetate CL (SQ5n)
    cl_charges = [cl_atoms[0]['charge']] if cl_atoms else []
    print(f"Read ion topology: NA={len(na_bead_types)} types, CL={len(cl_bead_types)} types from {ion_param_file}")
    
    # Parse water topology from ITP file
    water_param_file = os.path.join(SCRIPT_DIR, "ff3.00/martini_v3.0.0_solvents_v1.itp")
    water_topology = parse_itp_file(water_param_file)
    water_atoms = [atom for atom in water_topology['atoms'] if atom['residue'].upper() == 'W']
    water_bead_types = [atom['type'] for atom in water_atoms]
    water_charges = [atom['charge'] for atom in water_atoms]
    print(f"Read water topology: {len(water_bead_types)} bead types from {water_param_file}")
    
    # Read bead masses from force field file
    mass_file = "ff3.00/martini_v3.0.0.itp"
    martini_masses = read_martini_masses(mass_file)
    print(f"Read {len(martini_masses)} atom type masses from force field file")
    
    # Read DOPC bonds and angles from parsed topology (for both lipid and mixed systems)
    dopc_bonds = [(bond['i'], bond['j']) for bond in dopc_topology['bonds']]
    dopc_bond_lengths = [bond['r0'] for bond in dopc_topology['bonds']]  # nm
    dopc_bond_force_constants = [bond['k'] for bond in dopc_topology['bonds']]  # kJ/mol/nm²
    
    dopc_angles = [(angle['i'], angle['j'], angle['k']) for angle in dopc_topology['angles']]
    dopc_angle_equil_deg = [angle['theta0'] for angle in dopc_topology['angles']]  # degrees
    dopc_angle_force_constants = [angle['force_k'] for angle in dopc_topology['angles']]  # kJ/mol/rad²
    
    print(f"Read DOPC connectivity: {len(dopc_bonds)} bonds, {len(dopc_angles)} angles")
    
    # Validate that required topology data was found
    if not dopc_bonds:
        raise ValueError(f"FATAL ERROR: No DOPC bonds found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC bonds are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    if not dopc_angles:
        raise ValueError(f"FATAL ERROR: No DOPC angles found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC angles are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    # Unit conversions
    energy_conversion = 2.914952774272  # kJ/mol → E_up (kT at T=350.59K)
    length_conversion = 10.0  # nm → Å

    # Pressure conversion (for NPT simulations):
    # 1 atm = 101325 Pa = 101.325 kJ/m³
    # 1 atm = 1.01325e-28 kJ/Å³
    # 1 atm = (1.01325e-28 kJ/Å³) × (6.02214e23 mol⁻¹) / (2.914952774272 kJ/mol)
    # 1 atm = 0.000020933 E_up/Å³
    # 1 bar = 0.986923 atm = 0.000020659 E_up/Å³
    # For LAMMPS real units compatibility, use 1 bar = 0.000020659 E_up/Å³
    pressure_conversion_bar_to_eup = 0.000020659  # bar → E_up/Å³

    # Bonds: kJ/mol/nm² → E_up/Å²
    # = (kJ/mol → E_up) / (nm² → Å²)
    # = (1/energy_conversion) / (length_conversion²)
    bond_conversion = 1.0 / (energy_conversion * length_conversion ** 2)  # ≈ 0.003406

    # Angles: kJ/mol/deg² → E_up/deg²
    # = (kJ/mol → E_up) (degrees stay the same)
    angle_conversion = 1.0 / energy_conversion  # ≈ 0.343

    # Dihedrals: kJ/mol → E_up
    dihedral_conversion = 1.0 / energy_conversion  # ≈ 0.343

    print("\n=== MARTINI Unit Conversions ===")
    print(f"Bond lengths (nm -> Å): 0.40 nm -> 4.0 Å")
    print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): 7000.0 -> {7000.0 * bond_conversion:.3f}")
    print(f"Angle equilibrium (degrees): 108.0°")
    print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): 21.5 -> {21.5 * angle_conversion:.6f}")
    print(f"Dihedral force constants (kJ/mol -> E_up): 400.0 -> {400.0 * dihedral_conversion:.6f}")
    print(f"Pressure (bar -> E_up/Å³): 1 bar -> {pressure_conversion_bar_to_eup:.9f}")
    print(f"Energy conversion factor: {energy_conversion} (kJ/mol -> E_up)")
    print(f"Length conversion: 1 nm = {length_conversion} Å, so 1 nm² = {length_conversion**2} Å²")
    print(f"Bond conversion factor: {bond_conversion:.6f} (divide by energy_conv × length_conv²)")
    print(f"Angle conversion factor: {angle_conversion:.6f} (divide by energy_conv)")
    print(f"Dihedral conversion factor: {dihedral_conversion:.6f} (divide by energy_conv)")
    
    # Read PDB file
    if strict_from_martini_pdb or include_protein:
        input_pdb_file = os.path.join(SCRIPT_DIR, f'pdb/{pdb_id}.MARTINI.pdb')
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
    protein_itp = f"pdb/{pdb_id}_proa.itp"
    protein_topo_map = read_protein_itp_topology(protein_itp)
    
    # Parse protein connectivity for minimization stage (uses FLEXIBLE, NORMANG, POSRES sections)
    protein_bonds, protein_angles, protein_dihedrals, protein_constraints, protein_position_restraints = read_protein_itp_connectivity(protein_itp, 'minimization')
    protein_exclusions = read_protein_itp_exclusions(protein_itp)
    
    print(f"\n=== Protein Connectivity for Minimization Stage ===")
    print(f"Bonds: {len(protein_bonds)} (including FLEXIBLE bonds with large spring constants)")
    print(f"Angles: {len(protein_angles)} (from NORMANG section)")
    print(f"Dihedrals: {len(protein_dihedrals)}")
    print(f"Constraints: {len(protein_constraints)} (as bonds with large spring constants)")
    print(f"Position restraints: {len(protein_position_restraints)} (POSRES atoms - will be ignored)")
    
    
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
                    # Raise error for unknown protein atom
                    raise ValueError(f"FATAL ERROR: Unknown protein atom '{atom_name}' in residue '{residue_name}'.\n"
                                   f"  This indicates incomplete protein topology mapping in the ITP file.\n"
                                   f"  Please ensure the protein topology file '{protein_itp}' contains proper mapping for this atom.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
            elif residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file (for both lipid and mixed systems)
                if atom_name in dopc_atom_to_type:
                    martini_type = dopc_atom_to_type[atom_name]
                    charge = dopc_atom_to_charge[atom_name]
                else:
                    available_atom_names = sorted(dopc_atom_to_type.keys())
                    raise ValueError(f"FATAL ERROR: Unknown DOPC atom '{atom_name}' in residue '{residue_name}'.\n"
                                   f"  Available DOPC atom names: {available_atom_names}\n"
                                   f"  This indicates incomplete DOPC topology mapping.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
            elif residue_name == 'W':
                if not water_bead_types:
                    raise ValueError(f"FATAL ERROR: No water bead types found in topology.\n"
                                   f"  This indicates incomplete water parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = water_bead_types[0]
                charge = water_charges[0]
            elif residue_name == 'NA':
                if not na_bead_types:
                    raise ValueError(f"FATAL ERROR: No sodium bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = na_bead_types[0]
                charge = na_charges[0]
            elif residue_name == 'CL':
                if not cl_bead_types:
                    raise ValueError(f"FATAL ERROR: No chloride bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = cl_bead_types[0]
                charge = cl_charges[0]
            else:
                raise ValueError(f"FATAL ERROR: Unknown residue type '{residue_name}' for atom '{atom_name}'.\n"
                               f"  Supported residue types: PROTEIN, DOPC, W, NA, CL\n"
                               f"  This indicates incomplete system definition.\n"
                               f"  Aborting to prevent incorrect simulation results.")
            
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
            raise ValueError(f"FATAL ERROR: No CRYST1 record found in PDB file '{input_pdb_file}'.\n"
                           f"  Box dimensions are required for proper simulation setup.\n"
                           f"  Please ensure the PDB file contains a CRYST1 record with box dimensions.\n"
                           f"  Aborting to prevent incorrect simulation results.")
    
    # Print system parameters
    print(f"\n=== System Parameters ===")
    print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
    print(f"Box volume: {x_len * y_len * z_len:.1f} Å³")
    print(f"Total atoms: {n_atoms}")
    
    # Group atoms into molecules by residue ID and residue name
    molecules = []
    current_mol_atoms = []
    current_mol_names = []
    current_mol_indices = []
    current_resid = None
    current_resname = None
    
    # Define protein residue names
    protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
    
    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        # Determine molecule type
        if resname in protein_residues:
            mol_type = 'PROTEIN'
        else:
            # Normalize DOP (resname in PDB) to DOPC for reporting and selection
            mol_type = 'DOPC' if resname == 'DOP' else resname
            
        # Start new molecule if residue ID OR residue name changes
        if resid != current_resid or resname != current_resname:
            if current_mol_atoms:
                molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))
            current_mol_atoms = [atom_name]
            current_mol_names = [atom_name]
            current_mol_indices = [i]
            current_resid = resid
            current_resname = resname
            current_mol_type = mol_type
        else:
            current_mol_atoms.append(atom_name)
            current_mol_indices.append(i)
    
    if current_mol_atoms:
        molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))
    
    # Validate that each molecule contains only atoms of the same residue type
    print(f"\n=== Validating Molecule Integrity ===")
    for i, (mol_type, atoms, indices) in enumerate(molecules):
        residue_names_in_mol = [residue_names[idx] for idx in indices]
        unique_resnames = set(residue_names_in_mol)
        if len(unique_resnames) > 1:
            raise ValueError(f"FATAL ERROR: Molecule {i} contains mixed residue types: {unique_resnames}\n"
                           f"  This indicates incorrect molecule grouping.\n"
                           f"  Molecule atoms: {atoms}\n"
                           f"  Residue names: {residue_names_in_mol}\n"
                           f"  Aborting to prevent incorrect simulation results.")
        print(f"Molecule {i} ({mol_type}): {len(atoms)} atoms, residue type: {list(unique_resnames)[0]} ✓")
    
    # Count molecules by type, but group all protein residues together
    mol_counts = Counter()
    protein_residue_count = 0
    
    for mol_type, _, _ in molecules:
        if mol_type == 'PROTEIN':
            protein_residue_count += 1
        else:
            mol_counts[mol_type] += 1
    
    if protein_residue_count > 0:
        mol_counts['PROTEIN'] = f"1 chain ({protein_residue_count} residues)"
    
    dopc_count = mol_counts.get('DOPC', 0)
    water_count = mol_counts.get('W', 0)
    
    print(f"\n=== Molecule Summary ===")
    for moltype, count in mol_counts.items():
        print(f"{moltype}: {count} molecules")
    
    # Debug: Show first few molecules to verify proper separation
    print(f"\n=== Molecule Details (first 10) ===")
    for i, (mol_type, atoms, indices) in enumerate(molecules[:10]):
        print(f"Molecule {i}: {mol_type} with {len(atoms)} atoms (indices {indices[0]}-{indices[-1]})")
        print(f"  Atoms: {atoms}")
        print(f"  Residue IDs: {[residue_ids[idx] for idx in indices]}")
        print(f"  Residue Names: {[residue_names[idx] for idx in indices]}")
        print()
    
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
    dihedral_type_list = []
    
    # Create DOPC bonds and angles (for both lipid and mixed systems)
    dopc_molecules = [mol for mol in molecules if mol[0] == 'DOPC']  # unified label
    
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
        for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(dopc_angles):
            if (angle_idx1 < len(dopc_bead_types) and angle_idx2 < len(dopc_bead_types) and 
                angle_idx3 < len(dopc_bead_types)):
                atom1_name = atom_names_mol[angle_idx1]
                atom2_name = atom_names_mol[angle_idx2]
                atom3_name = atom_names_mol[angle_idx3]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                atom3 = name_to_idx[atom3_name]
                angles_list.append([atom1, atom2, atom3])
                angle_equil_deg_list.append(dopc_angle_equil_deg[i])
                angle_force_constants_list.append(dopc_angle_force_constants[i] * angle_conversion)  # kJ/mol/deg² to E_up/deg²
    
    print(f"Created {len(bonds_list)} bonds for {dopc_count} DOPC lipids")
    print(f"Created {len(angles_list)} angles for {dopc_count} DOPC lipids")
    
    # Create protein connectivity if available
    protein_bond_count = 0
    protein_angle_count = 0
    protein_dihedral_count = 0
    protein_constraint_count = 0
    
    if protein_bonds or protein_constraints:
        print(f"\n=== Protein Connectivity from {protein_itp} ===")
        print(f"Found {len(protein_bonds)} bonds, {len(protein_angles)} angles, {len(protein_dihedrals)} dihedrals")
        print(f"Found {len(protein_constraints)} constraints, {len(protein_position_restraints)} position restraints")
        
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
        
        # Add protein constraints as bonds with large spring constants
        for i, j, r0_nm, k_kj in protein_constraints:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_constraint_count += 1
        
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
        for i, j, k, l, phi0_deg, k_kj, func_type in protein_dihedrals:
            # Convert MARTINI units to UPSIDE units
            phi0_upside = phi0_deg  # degrees (same unit)
            k_upside = k_kj * dihedral_conversion  # kJ/mol to E_up
            
            # Add to dihedral list
            dihedrals_list.append([i, j, k, l])
            dihedral_equil_deg_list.append(phi0_upside)
            dihedral_force_constants_list.append(k_upside)
            dihedral_type_list.append(func_type)
            protein_dihedral_count += 1
        
        print(f"Added {protein_bond_count} protein bonds")
        print(f"Added {protein_constraint_count} protein constraints (as bonds with large spring constants)")
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
        mass = np.zeros(n_atoms, dtype='f4')
        for i, atom_type in enumerate(atom_types):
            # Get mass from force field file, raise error if not found
            if atom_type not in martini_masses:
                raise ValueError(f"FATAL ERROR: Mass not found for atom type '{atom_type}' (atom index {i}).\n"
                                f"  Available atom types with masses: {sorted(martini_masses.keys())}\n"
                                f"  This indicates incomplete force field parameters.\n"
                                f"  Aborting to prevent incorrect simulation results.")
            # Divide by 12.0 for reduced mass units (1 unit = 12 g/mol)
            mass[i] = martini_masses[atom_type] / 12.0
        
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        
        # Protein should NOT be held rigid during minimization
        # Allow the protein to relax and minimize its energy
        print("Protein atoms are free to move during minimization (no rigid constraints)")
        
        # Create stage-specific parameters group (always create this)
        stage_grp = t.create_group(input_grp, 'stage_parameters')
        stage_grp._v_attrs.enable = 1
        stage_grp._v_attrs.current_stage = b'minimization'
        
        # Store minimization stage bond parameters (large spring constants)
        min_bonds_grp = t.create_group(stage_grp, 'minimization_bonds')
        min_bond_fc = np.array([1000000.0] * len(protein_bonds), dtype='f4')
        t.create_array(min_bonds_grp, 'force_constants', obj=min_bond_fc)
        
        # Store production stage bond parameters (regular spring constants)
        prod_bonds_grp = t.create_group(stage_grp, 'production_bonds')
        prod_bond_fc = np.array([bond[3] for bond in protein_bonds], dtype='f4')  # Regular k values
        t.create_array(prod_bonds_grp, 'force_constants', obj=prod_bond_fc)
        
        # Store minimization stage angle parameters (NORMANG angles)
        min_angles_grp = t.create_group(stage_grp, 'minimization_angles')
        min_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(min_angles_grp, 'force_constants', obj=min_angle_fc)
        
        # Store production stage angle parameters (regular angles)
        prod_angles_grp = t.create_group(stage_grp, 'production_angles')
        # For production, we'd use different angle parameters (not NORMANG)
        # For now, use the same as minimization but this could be different
        prod_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(prod_angles_grp, 'force_constants', obj=prod_angle_fc)
        
        print(f"Stage-specific parameters: minimization bonds={len(min_bond_fc)}, production bonds={len(prod_bond_fc)}")
        print(f"Stage-specific parameters: minimization angles={len(min_angle_fc)}, production angles={len(prod_angle_fc)}")
        
        # ===================== NPT BAROSTAT CONFIGURATION =====================
        # Create barostat configuration group for NPT simulations
        # Settings are read from environment variables (set by run_sim_bilayer.sh)
        barostat_enable = int(os.environ.get('UPSIDE_NPT_ENABLE', '0'))
        if barostat_enable:
            print(f"\n=== Creating NPT Barostat Configuration ===")
            barostat_grp = t.create_group(input_grp, 'barostat')
            barostat_grp._v_attrs.enable = barostat_enable
            # Default: 1 bar = 0.000020659 E_up/Angstrom^3 (from 1 atm = 0.000020933215)
            barostat_grp._v_attrs.target_p_xy = float(os.environ.get('UPSIDE_NPT_TARGET_PXY', '0.000020659'))
            barostat_grp._v_attrs.target_p_z = float(os.environ.get('UPSIDE_NPT_TARGET_PZ', '0.000020659'))
            barostat_grp._v_attrs.tau_p = float(os.environ.get('UPSIDE_NPT_TAU', '1.0'))
            barostat_grp._v_attrs.compressibility = float(os.environ.get('UPSIDE_NPT_COMPRESSIBILITY', '14.521180763676'))
            barostat_grp._v_attrs.interval = int(os.environ.get('UPSIDE_NPT_INTERVAL', '10'))
            barostat_grp._v_attrs.semi_isotropic = int(os.environ.get('UPSIDE_NPT_SEMI', '1'))
            barostat_grp._v_attrs.debug = int(os.environ.get('UPSIDE_NPT_DEBUG', '1'))
            # Barostat type: 0 = Berendsen (default), 1 = Parrinello-Rahman
            barostat_grp._v_attrs.type = params['barostat_type']
            print(f"  Enabled: {barostat_enable}")
            print(f"  Type: {'Parrinello-Rahman' if barostat_grp._v_attrs.type == 1 else 'Berendsen'}")
            print(f"  Target Pxy: {barostat_grp._v_attrs.target_p_xy} E_up/Angstrom^3 (~1 bar)")
            print(f"  Target Pz: {barostat_grp._v_attrs.target_p_z} E_up/Angstrom^3 (~1 bar)")
            print(f"  Tau_p: {barostat_grp._v_attrs.tau_p}")
            print(f"  Interval: {barostat_grp._v_attrs.interval} steps")
        else:
            print(f"\n=== NPT Barostat Disabled (NVT mode) ===")
        
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
        # Note: epsilon and sigma attributes are required by C++ interface but not used in computation
        # They will be calculated from the coefficients array after it's created
        martini_potential._v_attrs.lj_cutoff = 12.0
        martini_potential._v_attrs.coul_cutoff = 12.0
        # dielectric constant is now included in the Coulomb k constant (31.775347952181)
        # coulomb_constant is now hardcoded as 31.775347952181 in the C++ code
        martini_potential._v_attrs.n_types = 1
        martini_potential._v_attrs.n_params = 4
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 1.0
        martini_potential._v_attrs.initialized = True
        force_cap = float(os.environ.get('UPSIDE_FORCE_CAP', '0'))
        martini_potential._v_attrs.force_cap = force_cap
        
        # PME configuration for long-range Coulomb interactions
        use_pme = int(os.environ.get('UPSIDE_USE_PME', '1'))
        pme_alpha = float(os.environ.get('UPSIDE_PME_ALPHA', '0.01'))
        pme_rcut = float(os.environ.get('UPSIDE_PME_RCUT', '10.0'))
        pme_nx = int(os.environ.get('UPSIDE_PME_NX', '32'))
        pme_ny = int(os.environ.get('UPSIDE_PME_NY', '32'))
        pme_nz = int(os.environ.get('UPSIDE_PME_NZ', '32'))
        pme_order = int(os.environ.get('UPSIDE_PME_ORDER', '4'))
        
        martini_potential._v_attrs.use_pme = use_pme
        martini_potential._v_attrs.pme_alpha = pme_alpha
        martini_potential._v_attrs.pme_rcut = pme_rcut
        martini_potential._v_attrs.pme_nx = pme_nx
        martini_potential._v_attrs.pme_ny = pme_ny
        martini_potential._v_attrs.pme_nz = pme_nz
        martini_potential._v_attrs.pme_order = pme_order

        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        # Optional softening/overwrite controls via environment variables
        # UPSIDE_SOFTEN_COULOMB: 1 to enable Slater softening for Coulomb
        # UPSIDE_SLATER_ALPHA: float value for Slater alpha (1/Angstrom)
        # UPSIDE_SOFTEN_LJ: 1 to enable soft-core LJ
        # UPSIDE_LJ_ALPHA: float value for LJ softening alpha (dimensionless)
        # UPSIDE_OVERWRITE_SPLINES: 1 to truncate spline debug files before writing
        soften_coul = int(os.environ.get('UPSIDE_SOFTEN_COULOMB', '0'))
        slater_alpha = float(os.environ.get('UPSIDE_SLATER_ALPHA', '1.0'))
        soften_lj = int(os.environ.get('UPSIDE_SOFTEN_LJ', '0'))
        lj_alpha = float(os.environ.get('UPSIDE_LJ_ALPHA', '1.0'))
        overwrite_splines = int(os.environ.get('UPSIDE_OVERWRITE_SPLINES', '0'))
        
        # PME configuration via environment variables
        # UPSIDE_USE_PME: 1 to enable Particle Mesh Ewald for long-range Coulomb
        # UPSIDE_PME_ALPHA: PME screening parameter (default: 0.2)
        # UPSIDE_PME_RCUT: Real space cutoff in Angstroms (default: 10.0)
        # UPSIDE_PME_NX/NY/NZ: Grid dimensions (default: 32, should be powers of 2)
        # UPSIDE_PME_ORDER: B-spline interpolation order (default: 4)
        use_pme = int(os.environ.get('UPSIDE_USE_PME', '0'))
        pme_alpha = float(os.environ.get('UPSIDE_PME_ALPHA', '0.01'))
        pme_rcut = float(os.environ.get('UPSIDE_PME_RCUT', '10.0'))
        pme_nx = int(os.environ.get('UPSIDE_PME_NX', '32'))
        pme_ny = int(os.environ.get('UPSIDE_PME_NY', '32'))
        pme_nz = int(os.environ.get('UPSIDE_PME_NZ', '32'))
        pme_order = int(os.environ.get('UPSIDE_PME_ORDER', '4'))

        # Set stage-specific softening parameters
        martini_potential._v_attrs.coulomb_soften = params['coulomb_soften']
        if params['coulomb_soften']:
            martini_potential._v_attrs.slater_alpha = params['slater_alpha']
        martini_potential._v_attrs.lj_soften = params['lj_soften']
        if params['lj_soften']:
            martini_potential._v_attrs.lj_soften_alpha = params['lj_alpha']
        martini_potential._v_attrs.overwrite_spline_tables = overwrite_splines
        
        # PME configuration
        martini_potential._v_attrs.use_pme = use_pme
        if use_pme:
            martini_potential._v_attrs.pme_alpha = pme_alpha
            martini_potential._v_attrs.pme_rcut = pme_rcut
            print(f"PME enabled: alpha={pme_alpha}, rcut={pme_rcut}, grid={pme_nx}x{pme_ny}x{pme_nz}, order={pme_order}")
        else:
            print("PME disabled: using standard Coulomb cutoff")

        martini_potential._v_attrs.debug_mode = 1  # Enable spline table generation
        martini_potential._v_attrs.force_debug_mode = 1  # Enable force debugging for charged particles
        
        # Periodic boundary potential removed - using NVT ensemble without boundaries
        
        # PME node creation removed - using Coulomb spline tables instead
        
        # NPT barostat configuration removed - using NVT ensemble without boundaries

        # Create atom indices and charges arrays for the potential
        t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
        t.create_array(martini_potential, 'charges', obj=charges)
        
        # Create pairs and coefficients for non-bonded interactions with proper exclusions
        pairs_list = []
        coeff_array = []
        
        # Create sets for exclusions (MARTINI uses nrexcl=1, so only 1-2 exclusions)
        bonded_pairs_12 = set()  # Directly bonded (1-2) - full exclusion
        additional_exclusions = set()  # Additional exclusions from ITP file
        
        # Add 1-2 exclusions from bond list
        for bond in bonds_list:
            sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
            bonded_pairs_12.add(sorted_bond)
        
        # Add additional exclusions from protein ITP file
        if protein_exclusions:
            for exclusion in protein_exclusions:
                sorted_exclusion = (min(exclusion[0], exclusion[1]), max(exclusion[0], exclusion[1]))
                additional_exclusions.add(sorted_exclusion)
        
        # Generate all unique pairs (i < j) with proper exclusions (nrexcl=1)
        excluded_12_count = 0
        excluded_additional_count = 0
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                pair = (i, j)
                
                # Skip 1-2 pairs (full exclusion) - nrexcl=1
                if pair in bonded_pairs_12:
                    excluded_12_count += 1
                    continue
                
                # Skip additional exclusions from ITP file
                if pair in additional_exclusions:
                    excluded_additional_count += 1
                    continue
                
                # No scaling needed for MARTINI (nrexcl=1)
                scale_factor = 1.0
                
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
                    # Raise error for missing interaction parameters
                    available_types = sorted(set([t[0] for t in martini_table.keys()] + [t[1] for t in martini_table.keys()]))
                    raise ValueError(f"FATAL ERROR: Missing interaction parameters for bead type pair ({type1}, {type2})\n"
                                   f"  Atom indices: {i} ({type1}) - {j} ({type2})\n"
                                   f"  This indicates incomplete MARTINI force field parameters.\n"
                                   f"  Available bead types in parameter table: {available_types}\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                
                # Convert to UPSIDE units
                epsilon = epsilon_kj / energy_conversion  # kJ/mol → E_up
                sigma = sigma_nm * length_conversion  # nm → Å
                q1 = charges[i] * scale_factor
                q2 = charges[j] * scale_factor
                coeff_array.append([epsilon * scale_factor, sigma, q1, q2])
        
        print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions (nrexcl=1)")
        print(f"Excluded {excluded_additional_count} additional pairs from ITP exclusions")

        t.create_array(martini_potential, 'pairs', obj=np.array(pairs_list, dtype=int))
        t.create_array(martini_potential, 'coefficients', obj=np.array(coeff_array, dtype='f4'))
        
        # Calculate representative epsilon and sigma from the coefficients array
        # These are required by the C++ interface but not used in computation
        if coeff_array:
            # Use the most common epsilon and sigma values from the coefficients
            epsilon_values = [coeff[0] for coeff in coeff_array]  # epsilon values
            sigma_values = [coeff[1] for coeff in coeff_array]    # sigma values
            
            # Use the median values as representative
            epsilon_values.sort()
            sigma_values.sort()
            median_epsilon = epsilon_values[len(epsilon_values)//2]
            median_sigma = sigma_values[len(sigma_values)//2]
            
            martini_potential._v_attrs.epsilon = median_epsilon
            martini_potential._v_attrs.sigma = median_sigma
        else:
            raise ValueError("FATAL ERROR: No interaction coefficients found.\n"
                           f"  This indicates no non-bonded interactions were defined.\n"
                           f"  Aborting to prevent incorrect simulation results.")

        # Add bonded potentials mirroring original run_martini.py
        # Bonds: dist_spring
        if bonds_list:
            bond_group = t.create_group(potential_grp, 'dist_spring')
            bond_group._v_attrs.arguments = np.array([b'pos'])
            bond_group._v_attrs.initialized = True

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
            # Store dihedral type information (1=periodic, 2=harmonic)
            t.create_array(dihedral_group, 'dihedral_type', obj=np.array(dihedral_type_list, dtype=int))
    
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
        f.write(f"Protein constraints: {protein_constraint_count}\n")
        f.write(f"Protein angles: {protein_angle_count}\n")
        f.write(f"Protein dihedrals: {protein_dihedral_count}\n")
        f.write(f"DOPC lipids: {dopc_count}\n")
        f.write(f"Water molecules: {water_count}\n")
        f.write(f"1-2 exclusions (nrexcl=1): {excluded_12_count}\n")
        f.write(f"Additional exclusions: {excluded_additional_count}\n")
    
    print(f"Preparation summary saved to: {summary_file}")

def create_production_input(input_file, pdb_id):
    """Create production input file with production-stage parameters"""
    print(f"\n=== Creating Production Input File ===")
    
    # Parse protein connectivity for production stage (uses regular sections)
    protein_itp = f"pdb/{pdb_id}_proa.itp"
    protein_bonds_prod, protein_angles_prod, protein_dihedrals_prod, protein_constraints_prod, protein_position_restraints_prod = read_protein_itp_connectivity(protein_itp, 'production')
    
    print(f"Production stage parameters:")
    print(f"  Bonds: {len(protein_bonds_prod)} (regular bonds)")
    print(f"  Angles: {len(protein_angles_prod)} (regular angles)")
    print(f"  Dihedrals: {len(protein_dihedrals_prod)}")
    print(f"  Constraints: {len(protein_constraints_prod)}")
    print(f"  Position restraints: {len(protein_position_restraints_prod)} (none for production)")
    
    # Create production input file
    production_file = f"outputs/martini_test/production.input.up"
    
    # Copy the minimization file and modify for production
    import shutil
    shutil.copy2(input_file, production_file)
    
    # Remove position restraints from production file
    import h5py
    with h5py.File(production_file, 'r+') as f:
        if 'input' in f and 'position_restraints' in f['input']:
            del f['input']['position_restraints']
            print("Removed position restraints from production file")
    
    print(f"Production input file created: {production_file}")
    print("Note: Production file uses regular parameters without position restraints")

def main_always_fixed(pdb_id):
    """Create input file with protein always fixed rigid throughout simulation"""
    print(f"\n=== Creating Input with Always-Fixed Protein ===")
    print(f"PDB ID: {pdb_id}")
    print("Output directory: outputs/martini_test")
    
    # Use the same setup as main() but modify the fix rigid configuration
    # First, run the normal preparation
    main()
    
    # Then modify the H5 file to ensure fix rigid is always enabled
    import h5py
    input_file = f"outputs/martini_test/test.input.up"
    
    print(f"\n=== Modifying H5 File for Always-Fixed Protein ===")
    with h5py.File(input_file, 'r+') as f:
        if 'input' in f and 'fix_rigid' in f['input']:
            fix_rigid_grp = f['input']['fix_rigid']
            # Ensure fix rigid is always enabled
            fix_rigid_grp.attrs['enable'] = 1
            print("Fix rigid enabled for entire simulation")
            
            # Remove stage-specific parameters since we want fix rigid throughout
            if 'stage_parameters' in f['input']:
                del f['input']['stage_parameters']
                print("Removed stage-specific parameters (using fix rigid throughout)")
        else:
            print("WARNING: No fix rigid group found in H5 file")
    
    print(f"Modified input file: {input_file}")
    print("Note: Protein will be fixed rigid throughout entire simulation")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2 and sys.argv[2] == '--production':
        # Create production input file
        create_production_input("outputs/martini_test/test.input.up", sys.argv[1])
    elif len(sys.argv) > 2 and sys.argv[2] == '--always-fix-protein':
        # Create input file with protein always fixed rigid
        print("=== Creating Input with Always-Fixed Protein ===")
        pdb_id = sys.argv[1]
        main_always_fixed(pdb_id)
    elif len(sys.argv) > 2 and sys.argv[2] == '--stage':
        # Create stage-specific input file
        if len(sys.argv) < 4:
            print("Usage: python prepare_martini.py <pdb_id> --stage <stage_name> [run_dir]")
            print("  Stages: minimization, npt_equil, npt_equil_reduced, npt_prod")
            sys.exit(1)
        stage = sys.argv[3]
        run_dir = sys.argv[4] if len(sys.argv) > 4 else "outputs/martini_test"
        main(stage=stage, run_dir=run_dir)
    else:
        # Create minimization input file (default)
        main()
