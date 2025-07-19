import sys
import h5py
import numpy as np
from collections import defaultdict

# Charges are now in standard units (1.0, -1.0) - no conversion needed
CHARGE_CONVERSION_FACTOR = 1.0

# DOPC bead names (for topology)
DOPC_NAMES = ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'D2A', 'C3A', 'C4A', 'C1B', 'D2B', 'C3B', 'C4B']

# Helper to infer molecule type from atom types
def infer_molecule_type(atom_types):
    if sorted(atom_types) == sorted(DOPC_NAMES):
        return 'DOPC'
    elif atom_types == ['W']:
        return 'WATER'
    elif atom_types == ['NA']:
        return 'NA'
    elif atom_types == ['CL']:
        return 'CL'
    else:
        return 'UNKNOWN'

def print_molecules_and_particles(h5):
    print("=== Molecules and Particles ===")
    residue_ids = h5['input/residue_ids'][:]
    atom_types = h5['input/type'][:].astype(str)
    molecules = defaultdict(list)
    for idx, resid in enumerate(residue_ids):
        molecules[resid].append(idx)
    for resid, particles in molecules.items():
        types = [atom_types[i] for i in particles]
        print(f"Molecule (residue_id={resid}): Particles={particles}, Types={types}")
    print(f"Total molecules: {len(molecules)}\n")


def print_pair_interactions(h5):
    print("=== Pair (Non-bonded) Interactions ===")
    pairs = h5['input/potential/martini_potential/pairs'][:]
    coeffs = h5['input/potential/martini_potential/coefficients'][:]
    for i, (a, b) in enumerate(pairs):
        epsilon, sigma, q1, q2 = coeffs[i]
        # Charges are already in standard units (1.0, -1.0)
        q1_real = q1
        q2_real = q2
        print(f"Pair: ({a}, {b}) | epsilon={epsilon:.3f}, sigma={sigma:.3f}, q1={q1_real:.2f}e, q2={q2_real:.2f}e")
    print(f"Total pair interactions: {len(pairs)}\n")


def print_bonded_interactions(h5):
    print("=== Bonded Interactions ===")
    # Bonds
    if 'input/potential/dist_spring' in h5:
        print("-- Bonds --")
        bonds = h5['input/potential/dist_spring/id'][:]
        equil_dist = h5['input/potential/dist_spring/equil_dist'][:]
        spring_const = h5['input/potential/dist_spring/spring_const'][:]
        for i, (a, b) in enumerate(bonds):
            print(f"Bond: ({a}, {b}) | length={equil_dist[i]:.2f}, k={spring_const[i]:.2f}")
        print(f"Total bonds: {len(bonds)}\n")
    else:
        print("No bonds found.\n")
    # Angles
    if 'input/potential/angle_spring' in h5:
        print("-- Angles --")
        angles = h5['input/potential/angle_spring/id'][:]
        equil_angle = h5['input/potential/angle_spring/equil_angle_deg'][:]
        spring_const = h5['input/potential/angle_spring/spring_const'][:]
        for i, (a, b, c) in enumerate(angles):
            print(f"Angle: ({a}, {b}, {c}) | angle={equil_angle[i]:.1f}°, k={spring_const[i]:.2f}")
        print(f"Total angles: {len(angles)}\n")
    else:
        print("No angles found.\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python visualize_up_setting.py <input.up>")
        sys.exit(1)
    up_file = sys.argv[1]
    with h5py.File(up_file, 'r') as h5:
        print_molecules_and_particles(h5)
        print_pair_interactions(h5)
        print_bonded_interactions(h5)

if __name__ == "__main__":
    main() 