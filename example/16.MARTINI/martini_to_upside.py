#!/usr/bin/env python

import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb
import argparse

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

# Import UPSIDE utilities
import run_upside as ru

def read_martini_pdb(pdb_file):
    """Read MARTINI PDB file and return positions of water beads."""
    positions = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                positions.append([x, y, z])
    return np.array(positions)

def create_martini_fasta(n_beads, output_file):
    """Create a minimal FASTA file with the correct number of residues."""
    with open(output_file, 'w') as f:
        f.write("> MARTINI\n")
        f.write("A" * n_beads)  # Use 'A' for all residues

def main():
    parser = argparse.ArgumentParser(description='Convert MARTINI PDB to UPSIDE format')
    parser.add_argument('pdb', help='Input MARTINI PDB file')
    parser.add_argument('output', help='Output base name (without extension)')
    args = parser.parse_args()

    # Create output directories
    input_dir = os.path.dirname(args.output) or '.'
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)

    # Read MARTINI positions
    positions = read_martini_pdb(args.pdb)
    n_beads = len(positions)

    # Create FASTA file
    fasta_file = f"{args.output}.fasta"
    create_martini_fasta(n_beads, fasta_file)

    # Create initial structure file
    # UPSIDE expects (n_atoms, 3) shape for initial structure
    # Each residue has 3 atoms (N, CA, C)
    n_atoms = n_beads * 3
    initial_pos = np.zeros((n_atoms, 3))
    
    # For each water bead, create 3 atoms in a line
    for i in range(n_beads):
        # Get the water bead position
        pos = positions[i]
        # Create 3 atoms in a line with 3.8Å spacing
        initial_pos[i*3] = pos  # N atom
        initial_pos[i*3+1] = pos + np.array([3.8, 0, 0])  # CA atom
        initial_pos[i*3+2] = pos + np.array([7.6, 0, 0])  # C atom

    # Save initial structure
    initial_file = f"{args.output}.initial.npy"
    np.save(initial_file, initial_pos)

    # Create UPSIDE configuration
    config_file = f"{args.output}.up"
    
    # Parameters
    param_dir_base = os.path.expanduser(upside_path+"/parameters/")
    param_dir_common = param_dir_base + "common/"
    param_dir_ff = param_dir_base + 'ff_2.1/'

    # Configure UPSIDE
    kwargs = dict(
        rama_library=param_dir_common + "rama.dat",
        rama_sheet_mix_energy=param_dir_ff + "sheet",
        reference_state_rama=param_dir_common + "rama_reference.pkl",
        hbond_energy=param_dir_ff + "hbond.h5",
        rotamer_placement=param_dir_ff + "sidechain.h5",
        dynamic_rotamer_1body=True,
        rotamer_interaction=param_dir_ff + "sidechain.h5",
        environment_potential=param_dir_ff + "environment.h5",
        bb_environment_potential=param_dir_ff + "bb_env.dat",
        initial_structure=initial_file,
    )

    # Generate base configuration
    config_stdout = ru.upside_config(fasta_file, config_file, **kwargs)
    print("Config commandline options:")
    print(config_stdout)

    # Add MARTINI potential
    with tb.open_file(config_file, 'r+') as t:
        # Create martini_potential group
        martini_group = t.create_group(t.root.input.potential, 'martini_potential')
        martini_group._v_attrs.arguments = np.array([b'pos'])
        martini_group._v_attrs.epsilon = 5.6  # kJ/mol
        martini_group._v_attrs.sigma = 0.47   # nm
        martini_group._v_attrs.lj_cutoff = 1.2  # nm
        martini_group._v_attrs.coul_cutoff = 1.2  # nm
        martini_group._v_attrs.dielectric = 15.0  # MARTINI water dielectric

        # Create arrays in the martini_potential group
        atoms = np.zeros(n_atoms, dtype=int)  # All atoms are water beads (type 0)
        t.create_array(martini_group, 'atom_indices', obj=atoms)
        t.create_array(martini_group, 'charges', obj=np.zeros(n_atoms))  # Water is neutral

        # Add MARTINI force field parameters
        # Create pairs array for MARTINI interactions
        # Only create pairs between water beads (every 3rd atom)
        n_pairs = 0
        pairs_list = []
        for i in range(0, n_atoms, 3):
            for j in range(i+3, n_atoms, 3):
                pairs_list.append([i, j])
                n_pairs += 1

        pairs_array = np.array(pairs_list, dtype=int)

        # Create coefficients array for MARTINI interactions
        coeff_array = np.zeros((n_pairs, 4))  # [epsilon, sigma, charge1, charge2]
        for i in range(n_pairs):
            coeff_array[i] = [5.6, 0.47, 0.0, 0.0]  # MARTINI water parameters

        # Add these arrays to the martini_potential group
        t.create_array(martini_group, 'pairs', obj=pairs_array)
        t.create_array(martini_group, 'coefficients', obj=coeff_array)

    print(f"\nDone! Created {config_file}")

if __name__ == '__main__':
    main()
