#!.venv/bin/python
import numpy as np
import tables as tb
import argparse
import sys
import re

def read_pdb_line(line):
    """Parse a single PDB line and return atom information."""
    if not line.startswith('ATOM') and not line.startswith('HETATM'):
        return None
        
    try:
        # Extract atom information from PDB line
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21]
        res_num = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        
        # Read occupancy (charge) with default 0.0
        try:
            occupancy = float(line[54:60])  # Used for charge in MARTINI
        except ValueError:
            occupancy = 0.0  # Default charge if not specified
            
        try:
            temp_factor = float(line[60:66])  # Used for LJ parameters in MARTINI
        except ValueError:
            temp_factor = 0.0  # Default temp factor if not specified
        
        return {
            'atom_name': atom_name,
            'res_name': res_name,
            'chain_id': chain_id,
            'res_num': res_num,
            'coords': np.array([x, y, z]),
            'occupancy': occupancy,  # Charge
            'temp_factor': temp_factor  # LJ parameter
        }
    except (ValueError, IndexError):
        return None

def read_martini_pdb(pdb_file):
    """Read MARTINI PDB file and extract coordinates and atom information."""
    coords = []
    atom_indices = []
    charges = []  # Will be filled with 0.0 for all atoms
    atom_types = []
    res_names = []
    chain_ids = []
    res_nums = []
    
    with open(pdb_file, 'r') as f:
        for i, line in enumerate(f):
            atom_info = read_pdb_line(line)
            if atom_info:
                coords.append(atom_info['coords'])
                atom_indices.append(i)
                charges.append(0.0)  # Set all charges to 0.0
                atom_types.append(atom_info['atom_name'])
                res_names.append(atom_info['res_name'])
                chain_ids.append(atom_info['chain_id'])
                res_nums.append(atom_info['res_num'])
    
    return (np.array(coords), np.array(atom_indices), np.array(charges), 
            atom_types, res_names, chain_ids, res_nums)

def create_upside_h5(coords, atom_indices, charges, atom_types, res_names, 
                     chain_ids, res_nums, output_file, epsilon=1.0, sigma=1.0, 
                     lj_cutoff=10.0, coul_cutoff=12.0, dielectric=80.0):
    """Create UPSIDE h5 file with pairs potential parameters."""
    
    with tb.open_file(output_file, 'w') as f:
        # Create root group
        root = f.root
        
        # Create all possible pairs of atoms
        n_atoms = len(atom_indices)
        pairs = []
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                pairs.append([i, j])
        pairs = np.array(pairs)
        
        # Create coefficients array (using LJ parameters)
        coefficients = np.ones(len(pairs)) * epsilon
        
        # Create atoms array (0 for CB, 1 for CA)
        atoms = np.zeros((len(pairs), 2), dtype=int)
        for i, (a1, a2) in enumerate(pairs):
            if atom_types[a1] == 'CA':
                atoms[i, 0] = 1
            if atom_types[a2] == 'CA':
                atoms[i, 1] = 1
        
        # Create centers array (using initial coordinates)
        centers = np.zeros((len(pairs), 3))
        for i, (a1, a2) in enumerate(pairs):
            centers[i] = (coords[a1] + coords[a2]) / 2
        
        # Store all arrays in root group
        f.create_array(root, 'pairs', pairs)
        f.create_array(root, 'coefficients', coefficients)
        f.create_array(root, 'atoms', atoms)
        f.create_array(root, 'centers', centers)
        
        # Store metadata
        root._v_attrs.num_atoms = len(atom_indices)
        root._v_attrs.num_chains = len(set(chain_ids))
        root._v_attrs.epsilon = epsilon
        root._v_attrs.sigma = sigma
        root._v_attrs.lj_cutoff = lj_cutoff
        root._v_attrs.coul_cutoff = coul_cutoff
        root._v_attrs.dielectric = dielectric

def main():
    parser = argparse.ArgumentParser(description='Convert MARTINI PDB to UPSIDE h5 input file')
    parser.add_argument('pdb_file', help='Input MARTINI PDB file')
    parser.add_argument('output_file', help='Output UPSIDE h5 file')
    parser.add_argument('--epsilon', type=float, default=1.0, help='LJ epsilon parameter')
    parser.add_argument('--sigma', type=float, default=1.0, help='LJ sigma parameter')
    parser.add_argument('--lj-cutoff', type=float, default=10.0, help='LJ cutoff distance')
    parser.add_argument('--coul-cutoff', type=float, default=12.0, help='Coulombic cutoff distance')
    parser.add_argument('--dielectric', type=float, default=80.0, help='Dielectric constant')
    
    args = parser.parse_args()
    
    try:
        # Read MARTINI PDB
        print(f"Reading MARTINI PDB file: {args.pdb_file}")
        coords, atom_indices, charges, atom_types, res_names, chain_ids, res_nums = read_martini_pdb(args.pdb_file)
        
        # Create UPSIDE h5 file
        print(f"Creating UPSIDE h5 file: {args.output_file}")
        create_upside_h5(
            coords, atom_indices, charges, atom_types, res_names, 
            chain_ids, res_nums, args.output_file,
            epsilon=args.epsilon,
            sigma=args.sigma,
            lj_cutoff=args.lj_cutoff,
            coul_cutoff=args.coul_cutoff,
            dielectric=args.dielectric
        )
        
        print("Conversion completed successfully!")
        print(f"Number of atoms processed: {len(atom_indices)}")
        print(f"Number of chains: {len(set(chain_ids))}")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
