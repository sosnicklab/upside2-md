#!/usr/bin/env python

import numpy as np
import tables as tb
import sys
import os
from upside_nodes import *
import argparse
import prody  # for PDB reading
import re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def bstring(string):
    return bytes(string, encoding="ascii")

# Default filter for HDF5 compression
default_filter = tb.Filters(complib='zlib', complevel=5, fletcher32=True)

def create_array(grp, nm, obj=None):
    return tb.create_earray(grp, nm, obj=obj, filters=default_filter)

def read_pdb_coordinates(pdb_file):
    """
    Read PDB file and extract coarse-grained particle coordinates
    Args:
        pdb_file: path to PDB file
    Returns:
        - coordinates: dictionary mapping particle names to coordinates
        - sequence: protein sequence
    """
    structure = prody.parsePDB(pdb_file)
    coordinates = {}
    sequence = []
    
    for chain in structure.iterChains():
        for residue in chain.iterResidues():
            resname = residue.getResname()
            sequence.append(resname)
            
            for atom in residue.iterAtoms():
                atom_name = atom.getName()
                coords = atom.getCoords()
                key = f"{resname}_{atom_name}"
                coordinates[key] = coords
    
    return coordinates, sequence

def read_itp_files(itp_files):
    """
    Read MARTINI ITP files for topology and interactions
    Args:
        itp_files: list of paths to ITP files
    Returns:
        - topology: dictionary containing topology information
        - interactions: dictionary containing interaction parameters
    """
    topology = {}
    interactions = {}
    
    for itp_file in itp_files:
        with open(itp_file, 'r') as f:
            content = f.read()
            
            # Parse [moleculetype] section
            moltype_match = re.search(r'\[ moleculetype \]\s*(.*?)(?=\[|$)', content, re.DOTALL)
            if moltype_match:
                topology['moleculetype'] = parse_moleculetype(moltype_match.group(1))
            
            # Parse [atoms] section
            atoms_match = re.search(r'\[ atoms \]\s*(.*?)(?=\[|$)', content, re.DOTALL)
            if atoms_match:
                topology['atoms'] = parse_atoms(atoms_match.group(1))
            
            # Parse [bonds] section
            bonds_match = re.search(r'\[ bonds \]\s*(.*?)(?=\[|$)', content, re.DOTALL)
            if bonds_match:
                interactions['bonds'] = parse_bonds(bonds_match.group(1))
            
            # Parse [angles] section
            angles_match = re.search(r'\[ angles \]\s*(.*?)(?=\[|$)', content, re.DOTALL)
            if angles_match:
                interactions['angles'] = parse_angles(angles_match.group(1))
            
            # Parse [dihedrals] section
            dihedrals_match = re.search(r'\[ dihedrals \]\s*(.*?)(?=\[|$)', content, re.DOTALL)
            if dihedrals_match:
                interactions['dihedrals'] = parse_dihedrals(dihedrals_match.group(1))
    
    return topology, interactions

def parse_moleculetype(content):
    """Parse [moleculetype] section from ITP file"""
    lines = content.strip().split('\n')
    return {
        'name': lines[0].split()[0],
        'nrexcl': int(lines[0].split()[1])
    }

def parse_atoms(content):
    """Parse [atoms] section from ITP file"""
    atoms = []
    for line in content.strip().split('\n'):
        if line and not line.startswith(';'):
            fields = line.split()
            atoms.append({
                'nr': int(fields[0]),
                'type': fields[1],
                'resnr': int(fields[2]),
                'residue': fields[3],
                'atom': fields[4],
                'cgnr': int(fields[5]),
                'charge': float(fields[6]),
                'mass': float(fields[7])
            })
    return atoms

def parse_bonds(content):
    """Parse [bonds] section from ITP file"""
    bonds = []
    for line in content.strip().split('\n'):
        if line and not line.startswith(';'):
            fields = line.split()
            bonds.append({
                'i': int(fields[0]),
                'j': int(fields[1]),
                'func': int(fields[2]),
                'parameters': [float(x) for x in fields[3:]]
            })
    return bonds

def parse_angles(content):
    """Parse [angles] section from ITP file"""
    angles = []
    for line in content.strip().split('\n'):
        if line and not line.startswith(';'):
            fields = line.split()
            angles.append({
                'i': int(fields[0]),
                'j': int(fields[1]),
                'k': int(fields[2]),
                'func': int(fields[3]),
                'parameters': [float(x) for x in fields[4:]]
            })
    return angles

def parse_dihedrals(content):
    """Parse [dihedrals] section from ITP file"""
    dihedrals = []
    for line in content.strip().split('\n'):
        if line and not line.startswith(';'):
            fields = line.split()
            dihedrals.append({
                'i': int(fields[0]),
                'j': int(fields[1]),
                'k': int(fields[2]),
                'l': int(fields[3]),
                'func': int(fields[4]),
                'parameters': [float(x) for x in fields[5:]]
            })
    return dihedrals

def convert_to_upside_format(pdb_data, itp_data):
    """
    Convert MARTINI data to Upside format
    Args:
        pdb_data: tuple of (coordinates, sequence) from PDB
        itp_data: tuple of (topology, interactions) from ITP files
    Returns:
        - fasta: protein sequence in FASTA format
        - initial_structure: backbone coordinates
        - interactions: converted interaction parameters
    """
    coordinates, sequence = pdb_data
    topology, interactions = itp_data
    
    # Convert sequence to FASTA format
    fasta = ''.join(three_letter_to_one_letter.get(res, 'X') for res in sequence)
    
    # Convert coordinates to Upside format
    initial_structure = convert_coordinates(coordinates, topology)
    
    # Convert interactions to Upside format
    upside_interactions = convert_interactions(interactions, topology)
    
    return {
        'fasta': fasta,
        'initial_structure': initial_structure,
        'interactions': upside_interactions
    }

def write_upside_config(output_file, fasta, initial_structure, interactions):
    """
    Write Upside configuration file
    Args:
        output_file: path to output .up file
        fasta: protein sequence
        initial_structure: backbone coordinates
        interactions: converted interaction parameters
    """
    t = tb.open_file(output_file, 'w')
    input_group = t.create_group(t.root, 'input')
    
    # Write basic structure
    write_basic_structure(input_group, fasta, initial_structure)
    
    # Write interactions
    write_interactions(input_group, interactions)
    
    t.close()

def write_basic_structure(group, fasta, initial_structure):
    """
    Write basic structure information to HDF5 group
    """
    # Write sequence
    create_array(group, 'sequence', obj=np.array([ord(c) for c in fasta]))
    
    # Write initial structure
    create_array(group, 'pos', obj=initial_structure)

def write_potentials(group, fasta):
    """
    Write various potential terms to HDF5 group
    """
    # Create potential group
    potential_group = group.create_group(group, 'potential')
    
    # Write backbone potentials
    write_backbone_potentials(potential_group, fasta)
    
    # Write sidechain potentials
    write_sidechain_potentials(potential_group, fasta)
    
    # Write environment potentials
    write_environment_potentials(potential_group, fasta)

def write_membrane_potential(group, fasta, membrane_params):
    """
    Write membrane potential to HDF5 group
    """
    # TODO: Implement membrane potential writing
    pass

def main():
    parser = argparse.ArgumentParser(description='Convert MARTINI input to Upside format')
    parser.add_argument('pdb_file', help='Input PDB file with coarse-grained coordinates')
    parser.add_argument('itp_files', nargs='+', help='Input ITP file(s) for topology and interactions')
    parser.add_argument('output', help='Output .up file')
    args = parser.parse_args()
    
    # Read PDB coordinates
    pdb_data = read_pdb_coordinates(args.pdb_file)
    
    # Read ITP files
    itp_data = read_itp_files(args.itp_files)
    
    # Convert to Upside format
    upside_data = convert_to_upside_format(pdb_data, itp_data)
    
    # Write Upside configuration
    write_upside_config(args.output, 
                       upside_data['fasta'],
                       upside_data['initial_structure'],
                       upside_data['interactions'])

if __name__ == '__main__':
    main()
