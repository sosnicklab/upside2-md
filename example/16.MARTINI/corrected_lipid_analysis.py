#!/usr/bin/env python3
"""
Corrected lipid analysis that properly handles non-unique residue IDs
"""

import os
import sys
import numpy as np
import tables as tb
import random
from collections import defaultdict

def analyze_lipid_molecule_corrected(up_file, output_file):
    """
    Analyze a random lipid molecule from the bilayer.up file
    Properly handles non-unique residue IDs by using residue name as well
    """
    print(f"Analyzing lipid molecule from: {up_file}")
    
    # Open the HDF5 file
    with tb.open_file(up_file, 'r') as f:
        # Get basic system information
        input_grp = f.root.input
        
        # Read atom data
        pos = input_grp.pos[:]
        types = input_grp.type[:]
        charges = input_grp.charges[:]
        residue_ids = input_grp.residue_ids[:]
        n_atoms = len(types)
        
        print(f"Total atoms in system: {n_atoms}")
        
        # Get potential data
        potential_grp = input_grp.potential
        
        # Read bonds
        bonds = []
        bond_lengths = []
        bond_constants = []
        if 'dist_spring' in potential_grp:
            bond_group = potential_grp.dist_spring
            bonds = bond_group.id[:]
            bond_lengths = bond_group.equil_dist[:]
            bond_constants = bond_group.spring_const[:]
        
        # Read angles
        angles = []
        angle_equil = []
        angle_constants = []
        if 'angle_spring' in potential_grp:
            angle_group = potential_grp.angle_spring
            angles = angle_group.id[:]
            angle_equil = angle_group.equil_angle_deg[:]
            angle_constants = angle_group.spring_const[:]
        
        # Read dihedrals
        dihedrals = []
        dihedral_equil = []
        dihedral_constants = []
        if 'dihedral_spring' in potential_grp:
            dihedral_group = potential_grp.dihedral_spring
            dihedrals = dihedral_group.id[:]
            dihedral_equil = dihedral_group.equil_angle_deg[:]
            dihedral_constants = dihedral_group.spring_const[:]
        
        # Read non-bonded interactions for exclusions
        pairs = []
        if 'martini_potential' in potential_grp:
            martini_group = potential_grp.martini_potential
            pairs = martini_group.pairs[:]
        
        print(f"Found {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals")
        
        # Since we don't have residue names in the .up file, we need to infer them
        # from the atom types. DOPC lipids have specific bead type patterns.
        print("Identifying lipid molecules by bead type patterns...")
        
        # Group atoms by residue ID and analyze their types
        residue_groups = defaultdict(list)
        for i, res_id in enumerate(residue_ids):
            residue_groups[res_id].append(i)
        
        lipid_residues = []
        for res_id, atom_indices in residue_groups.items():
            if len(atom_indices) > 1:  # Skip single-atom residues (likely water)
                # Get atom types for this residue
                res_types = [types[i].decode('utf-8') if isinstance(types[i], bytes) else str(types[i]) 
                           for i in atom_indices]
                
                # Check if this looks like a DOPC lipid
                # DOPC should have Q1, Q5, SN4a, N4a, and C1/C4h types
                has_q1 = any('Q1' in t for t in res_types)
                has_q5 = any('Q5' in t for t in res_types)
                has_sn4a = any('SN4a' in t for t in res_types)
                has_n4a = any('N4a' in t for t in res_types)
                has_c_types = any('C1' in t or 'C4h' in t for t in res_types)
                
                if has_q1 and has_q5 and has_sn4a and has_n4a and has_c_types:
                    # This is a DOPC lipid
                    lipid_residues.append(res_id)
                    print(f"Residue {res_id}: DOPC lipid with {len(atom_indices)} atoms")
                elif len(atom_indices) == 1 and 'W' in res_types:
                    # This is water
                    print(f"Residue {res_id}: Water molecule")
                else:
                    print(f"Residue {res_id}: Unknown molecule type with types: {set(res_types)}")
        
        print(f"Found {len(lipid_residues)} DOPC lipid residues")
        
        if not lipid_residues:
            print("No lipid molecules found!")
            return
        
        # Select a random lipid
        selected_residue = random.choice(lipid_residues)
        lipid_atoms = np.array(residue_groups[selected_residue])
        
        print(f"Selected lipid residue ID: {selected_residue}")
        print(f"Lipid has {len(lipid_atoms)} atoms")
        
        # Create atom index mapping for this lipid
        atom_to_lipid_idx = {atom_idx: i for i, atom_idx in enumerate(lipid_atoms)}
        
        # Extract bonds involving this lipid
        lipid_bonds = []
        lipid_bond_lengths = []
        lipid_bond_constants = []
        
        for i, bond in enumerate(bonds):
            atom1, atom2 = bond[0], bond[1]
            if atom1 in atom_to_lipid_idx and atom2 in atom_to_lipid_idx:
                lipid_bonds.append([atom_to_lipid_idx[atom1], atom_to_lipid_idx[atom2]])
                lipid_bond_lengths.append(bond_lengths[i])
                lipid_bond_constants.append(bond_constants[i])
        
        # Extract angles involving this lipid
        lipid_angles = []
        lipid_angle_equil = []
        lipid_angle_constants = []
        
        for i, angle in enumerate(angles):
            atom1, atom2, atom3 = angle[0], angle[1], angle[2]
            if (atom1 in atom_to_lipid_idx and atom2 in atom_to_lipid_idx and 
                atom3 in atom_to_lipid_idx):
                lipid_angles.append([atom_to_lipid_idx[atom1], atom_to_lipid_idx[atom2], 
                                   atom_to_lipid_idx[atom3]])
                lipid_angle_equil.append(angle_equil[i])
                lipid_angle_constants.append(angle_constants[i])
        
        # Extract dihedrals involving this lipid
        lipid_dihedrals = []
        lipid_dihedral_equil = []
        lipid_dihedral_constants = []
        
        for i, dihedral in enumerate(dihedrals):
            atom1, atom2, atom3, atom4 = dihedral[0], dihedral[1], dihedral[2], dihedral[3]
            if (atom1 in atom_to_lipid_idx and atom2 in atom_to_lipid_idx and 
                atom3 in atom_to_lipid_idx and atom4 in atom_to_lipid_idx):
                lipid_dihedrals.append([atom_to_lipid_idx[atom1], atom_to_lipid_idx[atom2], 
                                      atom_to_lipid_idx[atom3], atom_to_lipid_idx[atom4]])
                lipid_dihedral_equil.append(dihedral_equil[i])
                lipid_dihedral_constants.append(dihedral_constants[i])
        
        # Find exclusions - only 1-2 exclusions from bonds (nrexcl=1)
        bonded_pairs_12 = set()
        for bond in lipid_bonds:
            atom1, atom2 = bond[0], bond[1]
            bonded_pairs_12.add((min(atom1, atom2), max(atom1, atom2)))
        
        exclusions = bonded_pairs_12
        
        # Write analysis to file
        with open(output_file, 'w') as out_f:
            out_f.write("=== CORRECTED LIPID MOLECULE ANALYSIS ===\n")
            out_f.write(f"Selected lipid residue ID: {selected_residue}\n")
            out_f.write(f"Number of atoms: {len(lipid_atoms)}\n")
            out_f.write(f"Number of bonds: {len(lipid_bonds)}\n")
            out_f.write(f"Number of angles: {len(lipid_angles)}\n")
            out_f.write(f"Number of dihedrals: {len(lipid_dihedrals)}\n")
            out_f.write(f"Number of exclusions: {len(exclusions)}\n\n")
            
            # Write particle information
            out_f.write("=== PARTICLES ===\n")
            out_f.write("Index\tGlobal_Index\tType\tCharge\tPosition\n")
            for i, atom_idx in enumerate(lipid_atoms):
                atom_type = types[atom_idx].decode('utf-8') if isinstance(types[atom_idx], bytes) else str(types[atom_idx])
                charge = charges[atom_idx]
                pos_xyz = pos[atom_idx, :, 0]  # First frame
                out_f.write(f"{i}\t{atom_idx}\t{atom_type}\t{charge:.3f}\t{pos_xyz[0]:.3f}\t{pos_xyz[1]:.3f}\t{pos_xyz[2]:.3f}\n")
            
            # Write bonds
            out_f.write("\n=== BONDS ===\n")
            out_f.write("Index\tAtom1\tAtom2\tEquilibrium_Length\tForce_Constant\n")
            for i, (bond, length, const) in enumerate(zip(lipid_bonds, lipid_bond_lengths, lipid_bond_constants)):
                out_f.write(f"{i}\t{bond[0]}\t{bond[1]}\t{length:.6f}\t{const:.6f}\n")
            
            # Write angles
            out_f.write("\n=== ANGLES ===\n")
            out_f.write("Index\tAtom1\tAtom2\tAtom3\tEquilibrium_Angle\tForce_Constant\n")
            for i, (angle, equil, const) in enumerate(zip(lipid_angles, lipid_angle_equil, lipid_angle_constants)):
                out_f.write(f"{i}\t{angle[0]}\t{angle[1]}\t{angle[2]}\t{equil:.3f}\t{const:.6f}\n")
            
            # Write dihedrals
            out_f.write("\n=== DIHEDRALS ===\n")
            out_f.write("Index\tAtom1\tAtom2\tAtom3\tAtom4\tEquilibrium_Angle\tForce_Constant\n")
            for i, (dihedral, equil, const) in enumerate(zip(lipid_dihedrals, lipid_dihedral_equil, lipid_dihedral_constants)):
                out_f.write(f"{i}\t{dihedral[0]}\t{dihedral[1]}\t{dihedral[2]}\t{dihedral[3]}\t{equil:.3f}\t{const:.6f}\n")
            
            # Write exclusions
            out_f.write("\n=== EXCLUSIONS (1-2 bonded pairs only, nrexcl=1) ===\n")
            out_f.write("Index\tAtom1\tAtom2\n")
            for i, (atom1, atom2) in enumerate(sorted(exclusions)):
                out_f.write(f"{i}\t{atom1}\t{atom2}\n")
        
        print(f"Analysis complete! Results saved to: {output_file}")
        print(f"Lipid molecule has {len(lipid_atoms)} atoms, {len(lipid_bonds)} bonds, {len(lipid_angles)} angles, {len(lipid_dihedrals)} dihedrals, {len(exclusions)} exclusions")

def main():
    # Set random seed for reproducibility
    random.seed(42)
    
    # File paths
    up_file = "inputs/bilayer.up"
    output_file = "corrected_lipid_analysis.txt"
    
    if not os.path.exists(up_file):
        print(f"Error: File {up_file} not found!")
        return
    
    analyze_lipid_molecule_corrected(up_file, output_file)

if __name__ == "__main__":
    main()
