import numpy as np
import argparse
import importlib.util

# Simple PDB parser
def parse_pdb(filename):
    atoms = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                atom = {
                    'record': line[:6].strip(),
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'resName': line[17:20].strip(),
                    'chainID': line[21],
                    'resSeq': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'line': line
                }
                atoms.append(atom)
    return atoms

# Load mapping file provided by user
def load_mapping(mapping_file):
    spec = importlib.util.spec_from_file_location("MAP", mapping_file)
    MAP = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(MAP)
    return MAP.mapping

# Write PDB file
def write_pdb(atoms, filename):
    with open(filename, 'w') as file:
        for atom in atoms:
            file.write(atom['line'][:30] + f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}" + atom['line'][54:])

# RMSD alignment
def rmsd_align(ref_coords, target_coords):
    ref_centroid = np.mean(ref_coords, axis=0)
    target_centroid = np.mean(target_coords, axis=0)

    ref_centered = ref_coords - ref_centroid
    target_centered = target_coords - target_centroid

    H = np.dot(target_centered.T, ref_centered)
    U, _, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)

    aligned_coords = np.dot(target_centered, R) + ref_centroid
    return aligned_coords

# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='RMSD align MARTINI and all-atom structures or load backbone from npy.')
    parser.add_argument('-aa', '--all_atom', help='All-atom pdb file.')
    parser.add_argument('-cg', '--martini', required=True, help='MARTINI pdb file.')
    parser.add_argument('-map', '--mapping_file', required=True, help='Python mapping file (MAP.py).')
    parser.add_argument('-o', '--output', required=True, help='Output pdb file.')
    parser.add_argument('--to_aa', action='store_true', help='If set, update AA coordinates from MARTINI; otherwise, update MARTINI coordinates from AA.')
    parser.add_argument('--npy', help='Optional numpy file with backbone coordinates.')
    return parser.parse_args()

# Main alignment function
def align_and_update(aa_file, cg_file, mapping, output_file, to_aa, npy_file):
    cg_atoms = parse_pdb(cg_file)

    cg_names = list(mapping.keys())
    cg_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in cg_atoms if atom['name'] in cg_names])

    if npy_file:
        aa_coords = np.load(npy_file)
        updated_coords = rmsd_align(aa_coords, cg_coords)
        j = 0
        for atom in cg_atoms:
            if atom['name'] in cg_names:
                atom['x'], atom['y'], atom['z'] = updated_coords[j]
                j += 1
        write_pdb(cg_atoms, output_file)
    else:
        aa_atoms = parse_pdb(aa_file)
        aa_indices = [idx for indices in mapping.values() for idx in indices]
        aa_coords = np.array([[aa_atoms[i]['x'], aa_atoms[i]['y'], aa_atoms[i]['z']] for i in aa_indices])

        if to_aa:
            updated_coords = rmsd_align(cg_coords, aa_coords)
            for i, idx in enumerate(aa_indices):
                aa_atoms[idx]['x'], aa_atoms[idx]['y'], aa_atoms[idx]['z'] = updated_coords[i]
            write_pdb(aa_atoms, output_file)
        else:
            updated_coords = rmsd_align(aa_coords, cg_coords)
            j = 0
            for atom in cg_atoms:
                if atom['name'] in cg_names:
                    atom['x'], atom['y'], atom['z'] = updated_coords[j]
                    j += 1
            write_pdb(cg_atoms, output_file)

# Run the script
if __name__ == '__main__':
    args = parse_args()
    mapping = load_mapping(args.mapping_file)
    align_and_update(args.all_atom, args.martini, mapping, args.output, args.to_aa, args.npy)
