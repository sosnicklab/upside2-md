#!/usr/bin/env python3

import sys
import os
import numpy as np
import tables as tb
from collections import Counter, defaultdict

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

def analyze_molecules(input_file):
    """Analyze molecules in the input file"""
    print("=" * 80)
    print("1. MOLECULE ANALYSIS")
    print("=" * 80)
    
    with tb.open_file(input_file, 'r') as t:
        # Read atom types and residue IDs
        atom_types = t.root.input.type[:]
        residue_ids = t.root.input.residue_ids[:]
        charges = t.root.input.charges[:]
        pos = t.root.input.pos[:, :, 0]  # First frame positions
        
        # Convert bytes to strings
        atom_types = [at.decode('utf-8').strip() for at in atom_types]
        
        # Group atoms by residue ID
        molecules = defaultdict(list)
        for i, resid in enumerate(residue_ids):
            molecules[resid].append({
                'index': i,
                'type': atom_types[i],
                'charge': charges[i],
                'pos': pos[i]
            })
        
        # Analyze molecule types
        molecule_types = defaultdict(list)
        for resid, atoms in molecules.items():
            # Determine molecule type based on atom types
            atom_type_list = [atom['type'] for atom in atoms]
            
            # Count different atom types
            type_counts = Counter(atom_type_list)
            
            # Determine molecule type
            if 'Q1' in type_counts and 'Q5' in type_counts and 'SN4a' in type_counts:
                mol_type = 'DOPC'
            elif 'P4' in type_counts and len(atoms) == 1:
                mol_type = 'WATER'
            elif 'Qd' in type_counts and len(atoms) == 1:
                mol_type = 'SODIUM'
            elif 'Qa' in type_counts and len(atoms) == 1:
                mol_type = 'CHLORIDE'
            else:
                mol_type = 'UNKNOWN'
            
            molecule_types[mol_type].append({
                'resid': resid,
                'atoms': atoms,
                'type_counts': type_counts
            })
        
        # Print molecule analysis
        print(f"Total number of atoms: {len(atom_types)}")
        print(f"Total number of molecules: {len(molecules)}")
        print(f"Number of molecule types: {len(molecule_types)}")
        print()
        
        for mol_type, mol_list in molecule_types.items():
            print(f"Molecule type: {mol_type}")
            print(f"  Count: {len(mol_list)}")
            print(f"  Atom indices:")
            
            for i, mol in enumerate(mol_list):
                atom_indices = [atom['index'] for atom in mol['atoms']]
                print(f"    {mol_type} {i+1} (resid {mol['resid']}): {min(atom_indices)}-{max(atom_indices)} ({len(atom_indices)} atoms)")
                
                # Show atom details for first few molecules
                if i < 2:
                    print(f"      Atoms: {[f'{atom['index']}({atom['type']})' for atom in mol['atoms']]}")
            
            print()

def analyze_nonbonded_interactions(input_file):
    """Analyze nonbonded interactions"""
    print("=" * 80)
    print("2. NONBONDED INTERACTION TABLES")
    print("=" * 80)
    
    with tb.open_file(input_file, 'r') as t:
        # Check if martini_potential exists
        if 'martini_potential' in t.root.input.potential:
            martini_group = t.root.input.potential.martini_potential
            
            # Read pairs and coefficients
            pairs = martini_group.pairs[:]
            coefficients = martini_group.coefficients[:]
            
            print(f"Total nonbonded pairs: {len(pairs)}")
            print(f"Total coefficient sets: {len(coefficients)}")
            print()
            
            # Create epsilon and sigma tables
            epsilon_table = {}
            sigma_table = {}
            
            for i, (pair, coeff) in enumerate(zip(pairs, coefficients)):
                atom1_idx, atom2_idx = pair
                epsilon, sigma, q1, q2 = coeff
                
                # Get atom types
                atom_types = t.root.input.type[:]
                atom_types = [at.decode('utf-8').strip() for at in atom_types]
                type1, type2 = atom_types[atom1_idx], atom_types[atom2_idx]
                
                # Store in tables (use sorted key for consistency)
                key = tuple(sorted([type1, type2]))
                epsilon_table[key] = epsilon
                sigma_table[key] = sigma
            
            # Print epsilon table
            print("EPSILON TABLE (E_up units):")
            print("-" * 50)
            unique_keys = sorted(epsilon_table.keys())
            for key in unique_keys[:10]:  # Show first 10
                print(f"  {key[0]}-{key[1]}: {epsilon_table[key]:.6f}")
            if len(unique_keys) > 10:
                print(f"  ... and {len(unique_keys) - 10} more interactions")
            print()
            
            # Print sigma table
            print("SIGMA TABLE (Angstroms):")
            print("-" * 50)
            for key in unique_keys[:10]:  # Show first 10
                print(f"  {key[0]}-{key[1]}: {sigma_table[key]:.3f}")
            if len(unique_keys) > 10:
                print(f"  ... and {len(unique_keys) - 10} more interactions")
            print()
            
            # Show some example pairs
            print("EXAMPLE NONBONDED PAIRS:")
            print("-" * 50)
            for i in range(min(5, len(pairs))):
                pair = pairs[i]
                coeff = coefficients[i]
                atom_types = t.root.input.type[:]
                atom_types = [at.decode('utf-8').strip() for at in atom_types]
                type1, type2 = atom_types[pair[0]], atom_types[pair[1]]
                print(f"  Pair {i+1}: atoms {pair[0]}({type1})-{pair[1]}({type2})")
                print(f"    epsilon: {coeff[0]:.6f}, sigma: {coeff[1]:.3f}, charges: {coeff[2]:.1f}, {coeff[3]:.1f}")
            print()
        else:
            print("No martini_potential found in input file")

def analyze_bonded_interactions(input_file):
    """Analyze bonded interactions"""
    print("=" * 80)
    print("3. BONDED INTERACTION TABLES")
    print("=" * 80)
    
    with tb.open_file(input_file, 'r') as t:
        potential_group = t.root.input.potential
        
        # Check for bonds (dist_spring)
        if 'dist_spring' in potential_group:
            bond_group = potential_group.dist_spring
            bonds = bond_group.id[:]
            equil_dist = bond_group.equil_dist[:]
            spring_const = bond_group.spring_const[:]
            
            print("BONDS:")
            print("-" * 30)
            print(f"Total bonds: {len(bonds)}")
            print()
            
            # Get atom types for bond analysis
            atom_types = t.root.input.type[:]
            atom_types = [at.decode('utf-8').strip() for at in atom_types]
            
            for i in range(min(10, len(bonds))):
                bond = bonds[i]
                atom1_idx, atom2_idx = bond
                type1, type2 = atom_types[atom1_idx], atom_types[atom2_idx]
                print(f"  Bond {i+1}: atoms {atom1_idx}({type1})-{atom2_idx}({type2})")
                print(f"    equilibrium distance: {equil_dist[i]:.3f} Å")
                print(f"    spring constant: {spring_const[i]:.6f} E_up/Å²")
            if len(bonds) > 10:
                print(f"  ... and {len(bonds) - 10} more bonds")
            print()
        else:
            print("No bonds (dist_spring) found")
        
        # Check for angles (angle_spring)
        if 'angle_spring' in potential_group:
            angle_group = potential_group.angle_spring
            angles = angle_group.id[:]
            equil_angle = angle_group.equil_angle_deg[:]
            spring_const = angle_group.spring_const[:]
            
            print("ANGLES:")
            print("-" * 30)
            print(f"Total angles: {len(angles)}")
            print()
            
            # Get atom types for angle analysis
            atom_types = t.root.input.type[:]
            atom_types = [at.decode('utf-8').strip() for at in atom_types]
            
            for i in range(min(10, len(angles))):
                angle = angles[i]
                atom1_idx, atom2_idx, atom3_idx = angle
                type1, type2, type3 = atom_types[atom1_idx], atom_types[atom2_idx], atom_types[atom3_idx]
                print(f"  Angle {i+1}: atoms {atom1_idx}({type1})-{atom2_idx}({type2})-{atom3_idx}({type3})")
                print(f"    equilibrium angle: {equil_angle[i]:.1f}°")
                print(f"    spring constant: {spring_const[i]:.6f} E_up/deg²")
            if len(angles) > 10:
                print(f"  ... and {len(angles) - 10} more angles")
            print()
        else:
            print("No angles (angle_spring) found")
        
        # Check for dihedrals (torsion_spring)
        if 'torsion_spring' in potential_group:
            torsion_group = potential_group.torsion_spring
            torsions = torsion_group.id[:]
            equil_torsion = torsion_group.equil_torsion_deg[:]
            spring_const = torsion_group.spring_const[:]
            
            print("DIHEDRALS:")
            print("-" * 30)
            print(f"Total dihedrals: {len(torsions)}")
            print()
            
            # Get atom types for torsion analysis
            atom_types = t.root.input.type[:]
            atom_types = [at.decode('utf-8').strip() for at in atom_types]
            
            for i in range(min(10, len(torsions))):
                torsion = torsions[i]
                atom1_idx, atom2_idx, atom3_idx, atom4_idx = torsion
                type1, type2, type3, type4 = atom_types[atom1_idx], atom_types[atom2_idx], atom_types[atom3_idx], atom_types[atom4_idx]
                print(f"  Dihedral {i+1}: atoms {atom1_idx}({type1})-{atom2_idx}({type2})-{atom3_idx}({type3})-{atom4_idx}({type4})")
                print(f"    equilibrium torsion: {equil_torsion[i]:.1f}°")
                print(f"    spring constant: {spring_const[i]:.6f} E_up/deg²")
            if len(torsions) > 10:
                print(f"  ... and {len(torsions) - 10} more dihedrals")
            print()
        else:
            print("No dihedrals (torsion_spring) found")

def analyze_spline_tables(input_file):
    """Analyze spline tables for interactions"""
    print("=" * 80)
    print("4. SPLINE TABLES")
    print("=" * 80)
    
    print("Note: Spline tables are generated internally by UPSIDE during simulation.")
    print("The actual spline data is not stored in the input file but computed on-the-fly.")
    print("However, we can show the parameters that would be used to generate splines:")
    print()
    
    with tb.open_file(input_file, 'r') as t:
        potential_group = t.root.input.potential
        
        # Nonbonded spline parameters
        if 'martini_potential' in potential_group:
            martini_group = potential_group.martini_potential
            print("NONBONDED SPLINE PARAMETERS:")
            print("-" * 40)
            print(f"  LJ cutoff: {martini_group._v_attrs.lj_cutoff:.1f} Å")
            print(f"  Coulomb cutoff: {martini_group._v_attrs.coul_cutoff:.1f} Å")
            print(f"  Cache buffer: {martini_group._v_attrs.cache_buffer:.1f} Å")
            print(f"  Debug mode: {martini_group._v_attrs.debug_mode}")
            print()
            
            # Show some example coefficients that would be used for splines
            coefficients = martini_group.coefficients[:]
            print("EXAMPLE NONBONDED COEFFICIENTS (for spline generation):")
            print("-" * 60)
            for i in range(min(5, len(coefficients))):
                coeff = coefficients[i]
                print(f"  Interaction {i+1}: epsilon={coeff[0]:.6f}, sigma={coeff[1]:.3f}, q1={coeff[2]:.1f}, q2={coeff[3]:.1f}")
            print()
        
        # Bonded spline parameters
        if 'dist_spring' in potential_group:
            bond_group = potential_group.dist_spring
            print("BOND SPLINE PARAMETERS:")
            print("-" * 40)
            print(f"  Debug mode: {bond_group._v_attrs.debug_mode}")
            print()
            
            # Show some example bond parameters
            equil_dist = bond_group.equil_dist[:]
            spring_const = bond_group.spring_const[:]
            print("EXAMPLE BOND PARAMETERS (for spline generation):")
            print("-" * 50)
            for i in range(min(5, len(equil_dist))):
                print(f"  Bond {i+1}: equil_dist={equil_dist[i]:.3f} Å, spring_const={spring_const[i]:.6f} E_up/Å²")
            print()
        
        if 'angle_spring' in potential_group:
            angle_group = potential_group.angle_spring
            print("ANGLE SPLINE PARAMETERS:")
            print("-" * 40)
            print(f"  Debug mode: {angle_group._v_attrs.debug_mode}")
            print()
            
            # Show some example angle parameters
            equil_angle = angle_group.equil_angle_deg[:]
            spring_const = angle_group.spring_const[:]
            print("EXAMPLE ANGLE PARAMETERS (for spline generation):")
            print("-" * 50)
            for i in range(min(5, len(equil_angle))):
                print(f"  Angle {i+1}: equil_angle={equil_angle[i]:.1f}°, spring_const={spring_const[i]:.6f} E_up/deg²")
            print()

def main():
    """Main function to generate debug output"""
    input_file = "inputs/test.up"
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found!")
        print("Please run the main simulation first to generate the input file.")
        return
    
    print("MARTINI 3.00 SIMULATION DEBUG OUTPUT")
    print("=" * 80)
    print()
    
    # Run all analyses
    analyze_molecules(input_file)
    analyze_nonbonded_interactions(input_file)
    analyze_bonded_interactions(input_file)
    analyze_spline_tables(input_file)
    
    print("=" * 80)
    print("DEBUG OUTPUT COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
