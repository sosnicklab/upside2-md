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

def detect_potential_type(h5):
    """Detect which MARTINI potential is being used"""
    if 'input/potential/martini_potential' in h5:
        return 'martini_potential'
    else:
        return None

def detect_bonded_potential_types(h5):
    """Detect which bonded potentials are being used"""
    bond_type = None
    angle_type = None
    
    if 'input/potential/dist_spring' in h5:
        bond_type = 'dist_spring'
        
    if 'input/potential/angle_spring' in h5:
        angle_type = 'angle_spring'
        
    return bond_type, angle_type

def print_potential_info(h5, potential_type):
    """Print information about the MARTINI potential being used"""
    print("=== MARTINI Potential Configuration ===")
    
    if potential_type == 'martini_potential':
        print("🔧 Using MARTINI potential with spline interpolation")
        print("   - LJ interactions: Precomputed spline tables (1000 knots)")
        print("   - Coulomb interactions: Precomputed spline tables (1000 knots)")
        print("   - Expected performance improvement for LJ and Coulomb calculations")
        
        # Check if bonded potentials are also used
        bond_type, angle_type = detect_bonded_potential_types(h5)
        if bond_type == 'dist_spring':
            print("   - Bond interactions: Precomputed spline tables (1000 knots)")
        if angle_type == 'angle_spring':
            print("   - Angle interactions: Precomputed spline tables (1000 knots)")
    else:
        print("❌ No MARTINI potential found!")
        return
    
    # Get potential parameters
    pot_group = h5[f'input/potential/{potential_type}']
    
    # Display key parameters
    epsilon = pot_group.attrs.get('epsilon', 'N/A')
    sigma = pot_group.attrs.get('sigma', 'N/A')
    lj_cutoff = pot_group.attrs.get('lj_cutoff', 'N/A')
    coul_cutoff = pot_group.attrs.get('coul_cutoff', 'N/A')
    dielectric = pot_group.attrs.get('dielectric', 'N/A')
    coulomb_constant = pot_group.attrs.get('coulomb_constant', 'N/A')
    
    print(f"   - Epsilon: {epsilon}")
    print(f"   - Sigma: {sigma} Å")
    print(f"   - LJ cutoff: {lj_cutoff} Å")
    print(f"   - Coulomb cutoff: {coul_cutoff} Å")
    print(f"   - Dielectric constant: {dielectric}")
    print(f"   - Coulomb constant: {coulomb_constant}")
    
    # Display box dimensions if available
    if 'x_len' in pot_group.attrs:
        x_len = pot_group.attrs['x_len']
        y_len = pot_group.attrs['y_len']
        z_len = pot_group.attrs['z_len']
        print(f"   - Box dimensions: {x_len:.1f} × {y_len:.1f} × {z_len:.1f} Å")
    
    print()

def print_pair_interactions(h5, potential_type):
    print("=== Pair (Non-bonded) Interactions ===")
    
    if potential_type is None:
        print("❌ No MARTINI potential found!")
        return
    
    pairs = h5[f'input/potential/{potential_type}/pairs'][:]
    coeffs = h5[f'input/potential/{potential_type}/coefficients'][:]
    
    # Count different interaction types
    lj_interactions = 0
    coulomb_interactions = 0
    both_interactions = 0
    
    for i, (a, b) in enumerate(pairs):
        epsilon, sigma, q1, q2 = coeffs[i]
        # Charges are already in standard units (1.0, -1.0)
        q1_real = q1
        q2_real = q2
        
        # Count interaction types
        has_lj = epsilon != 0.0 and sigma != 0.0
        has_coulomb = q1 != 0.0 and q2 != 0.0
        
        if has_lj and has_coulomb:
            both_interactions += 1
        elif has_lj:
            lj_interactions += 1
        elif has_coulomb:
            coulomb_interactions += 1
        
        # Only print first few interactions to avoid overwhelming output
        if i < 10:
            interaction_types = []
            if has_lj:
                interaction_types.append("LJ")
            if has_coulomb:
                interaction_types.append("Coulomb")
            types_str = "+".join(interaction_types) if interaction_types else "None"
            
            print(f"Pair {i}: ({a}, {b}) | Types: {types_str} | "
                  f"ε={epsilon:.3f}, σ={sigma:.3f}, q1={q1_real:.2f}e, q2={q2_real:.2f}e")
    
    if len(pairs) > 10:
        print(f"... and {len(pairs) - 10} more pairs")
    
    print(f"\nInteraction Summary:")
    print(f"  - Total pairs: {len(pairs)}")
    print(f"  - LJ only: {lj_interactions}")
    print(f"  - Coulomb only: {coulomb_interactions}")
    print(f"  - LJ + Coulomb: {both_interactions}")
    print()

def print_bonded_interactions(h5):
    print("=== Bonded Interactions ===")
    
    # Detect bonded potential types
    bond_type, angle_type = detect_bonded_potential_types(h5)
    
    # Bonds
    if bond_type:
        print("-- Bonds --")
        print("🔧 Using bond potential with spline interpolation")
            
        bonds = h5[f'input/potential/{bond_type}/id'][:]
        equil_dist = h5[f'input/potential/{bond_type}/equil_dist'][:]
        spring_const = h5[f'input/potential/{bond_type}/spring_const'][:]
        
        # Show first few bonds
        for i, (a, b) in enumerate(bonds[:10]):
            print(f"Bond {i}: ({a}, {b}) | length={equil_dist[i]:.2f} Å, k={spring_const[i]:.2f}")
        
        if len(bonds) > 10:
            print(f"... and {len(bonds) - 10} more bonds")
        print(f"Total bonds: {len(bonds)}\n")
    else:
        print("No bonds found.\n")
    
    # Angles
    if angle_type:
        print("-- Angles --")
        print("🔧 Using angle potential with spline interpolation")
            
        angles = h5[f'input/potential/{angle_type}/id'][:]
        equil_angle = h5[f'input/potential/{angle_type}/equil_angle_deg'][:]
        spring_const = h5[f'input/potential/{angle_type}/spring_const'][:]
        
        # Show first few angles
        for i, (a, b, c) in enumerate(angles[:10]):
            print(f"Angle {i}: ({a}, {b}, {c}) | angle={equil_angle[i]:.1f}°, k={spring_const[i]:.2f}")
        
        if len(angles) > 10:
            print(f"... and {len(angles) - 10} more angles")
        print(f"Total angles: {len(angles)}\n")
    else:
        print("No angles found.\n")

def print_system_summary(h5):
    """Print a summary of the system"""
    print("=== System Summary ===")
    
    # Count atoms by type
    atom_types = h5['input/type'][:].astype(str)
    type_counts = defaultdict(int)
    for atom_type in atom_types:
        type_counts[atom_type] += 1
    
    print("Atom type counts:")
    for atom_type, count in sorted(type_counts.items()):
        print(f"  {atom_type}: {count}")
    
    total_atoms = len(atom_types)
    print(f"\nTotal atoms: {total_atoms}")
    
    # Count molecules
    residue_ids = h5['input/residue_ids'][:]
    unique_residues = len(set(residue_ids))
    print(f"Total molecules: {unique_residues}")
    
    # Check for periodic boundary conditions
    if 'input/potential/periodic_boundary_potential' in h5:
        print("Periodic boundary conditions: Enabled")
    else:
        print("Periodic boundary conditions: Disabled")
    
    print()

def main():
    if len(sys.argv) != 2:
        print("Usage: python visualize_up_setting.py <input.up>")
        print("\nThis script visualizes the configuration of a UPSIDE input file,")
        print("including MARTINI potential settings, interactions, and system summary.")
        sys.exit(1)
    
    up_file = sys.argv[1]
    
    try:
        with h5py.File(up_file, 'r') as h5:
            print(f"📁 Analyzing UPSIDE input file: {up_file}")
            print("=" * 60)
            
            # Detect potential type
            potential_type = detect_potential_type(h5)
            
            # Print system information
            print_system_summary(h5)
            print_molecules_and_particles(h5)
            
            # Print potential information
            print_potential_info(h5, potential_type)
            
            # Print interactions
            print_pair_interactions(h5, potential_type)
            print_bonded_interactions(h5)
            
            print("✅ Analysis complete!")
            
    except FileNotFoundError:
        print(f"❌ Error: File '{up_file}' not found!")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 