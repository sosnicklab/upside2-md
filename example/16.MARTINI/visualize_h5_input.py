#!/usr/bin/env python3
"""
Comprehensive HDF5 Input File Visualizer for MARTINI Simulations
Shows all particles, types, interactions, and spline tables
"""

import h5py
import numpy as np
import sys
import os

def format_coord(coord):
    """Format coordinate with consistent spacing"""
    return f"({float(coord[0]):8.3f}, {float(coord[1]):8.3f}, {float(coord[2]):8.3f})"

def format_charge(charge):
    """Format charge with consistent spacing"""
    return f"{float(charge):6.2f}"

def get_particle_type(particle_id, atom_types):
    """Get particle type from atom types"""
    if particle_id < len(atom_types):
        atom_type = atom_types[particle_id].decode('utf-8') if isinstance(atom_types[particle_id], bytes) else atom_types[particle_id]
        return atom_type
    return "Unknown"

def analyze_h5_file(h5_file_path, output_file=None, spline_dir="."):
    """Analyze HDF5 file and generate comprehensive report"""
    
    if output_file is None:
        output_file = "h5_visualization_report.txt"
    
    with open(output_file, 'w') as f_out:
        f_out.write("=" * 80 + "\n")
        f_out.write("COMPREHENSIVE HDF5 INPUT FILE VISUALIZATION REPORT\n")
        f_out.write("=" * 80 + "\n\n")
        
        with h5py.File(h5_file_path, 'r') as h5:
            
            # ====================================================================
            # 1. SYSTEM OVERVIEW
            # ====================================================================
            f_out.write("1. SYSTEM OVERVIEW\n")
            f_out.write("-" * 40 + "\n")
            
            # Box dimensions
            if 'input/box' in h5:
                box = h5['input/box'][:]
                f_out.write(f"Box dimensions: {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} Å\n")
                f_out.write(f"Box volume: {box[0] * box[1] * box[2]:.3f} Å³\n\n")
            
            # Particle count
            if 'input/pos' in h5:
                n_particles = len(h5['input/pos'])
                f_out.write(f"Total particles: {n_particles}\n\n")
            
            # ====================================================================
            # 2. ALL PARTICLES DETAILED LISTING
            # ====================================================================
            f_out.write("2. ALL PARTICLES DETAILED LISTING\n")
            f_out.write("-" * 40 + "\n")
            
            if 'input/pos' in h5:
                positions = h5['input/pos'][:]
                # Get charges - try both possible locations
                charges = None
                if 'input/charges' in h5:
                    charges = h5['input/charges'][:]
                elif 'input/potential/martini_potential/charges' in h5:
                    charges = h5['input/potential/martini_potential/charges'][:]
                else:
                    charges = np.zeros(len(positions))
                
                # Get atom types if available
                atom_types = h5['input/type'][:] if 'input/type' in h5 else [b'UNK'] * len(positions)
                
                f_out.write(f"{'ID':>4} {'Type':>10} {'Position':>30} {'Charge':>8}\n")
                f_out.write("-" * 60 + "\n")
                
                for i in range(len(positions)):
                    particle_type = get_particle_type(i, atom_types)
                    pos_str = format_coord(positions[i])
                    charge_str = format_charge(charges[i])
                    f_out.write(f"{i:4d} {particle_type:>10} {pos_str:>30} {charge_str:>8}\n")
                
                f_out.write("\n")
            
            # ====================================================================
            # 3. BONDED INTERACTIONS
            # ====================================================================
            f_out.write("3. BONDED INTERACTIONS\n")
            f_out.write("-" * 40 + "\n")
            
            # Bonds
            if 'input/potential/dist_spring' in h5:
                bond_group = h5['input/potential/dist_spring']
                if 'id' in bond_group:
                    bond_ids = bond_group['id'][:]
                    bond_equil = bond_group['equil_dist'][:] if 'equil_dist' in bond_group else np.zeros(len(bond_ids))
                    bond_force = bond_group['spring_const'][:] if 'spring_const' in bond_group else np.zeros(len(bond_ids))
                    
                    f_out.write(f"BONDS ({len(bond_ids)} total):\n")
                    f_out.write("-" * 20 + "\n")
                    f_out.write(f"{'ID':>4} {'Atoms':>10} {'Length':>10} {'Force':>12}\n")
                    f_out.write("-" * 40 + "\n")
                    
                    for i in range(min(20, len(bond_ids))):  # Show first 20
                        atoms = f"{bond_ids[i][0]:3d}-{bond_ids[i][1]:3d}"
                        f_out.write(f"{i:4d} {atoms:>10} {bond_equil[i]:8.3f}Å {bond_force[i]:10.3f}\n")
                    
                    if len(bond_ids) > 20:
                        f_out.write(f"... and {len(bond_ids) - 20} more bonds\n")
                    f_out.write("\n")
            
            # Angles
            if 'input/potential/angle_spring' in h5:
                angle_group = h5['input/potential/angle_spring']
                if 'id' in angle_group:
                    angle_ids = angle_group['id'][:]
                    angle_equil = angle_group['equil_angle_deg'][:] if 'equil_angle_deg' in angle_group else np.zeros(len(angle_ids))
                    angle_force = angle_group['spring_const'][:] if 'spring_const' in angle_group else np.zeros(len(angle_ids))
                    
                    f_out.write(f"ANGLES ({len(angle_ids)} total):\n")
                    f_out.write("-" * 20 + "\n")
                    f_out.write(f"{'ID':>4} {'Atoms':>15} {'Angle':>10} {'Force':>12}\n")
                    f_out.write("-" * 45 + "\n")
                    
                    for i in range(min(20, len(angle_ids))):  # Show first 20
                        atoms = f"{angle_ids[i][0]:3d}-{angle_ids[i][1]:3d}-{angle_ids[i][2]:3d}"
                        f_out.write(f"{i:4d} {atoms:>15} {angle_equil[i]:8.1f}° {angle_force[i]:10.3f}\n")
                    
                    if len(angle_ids) > 20:
                        f_out.write(f"... and {len(angle_ids) - 20} more angles\n")
                    f_out.write("\n")
            
            # Dihedrals
            if 'input/potential/dihedral_spring' in h5:
                dihedral_group = h5['input/potential/dihedral_spring']
                if 'id' in dihedral_group:
                    dihedral_ids = dihedral_group['id'][:]
                    dihedral_equil = dihedral_group['equil_angle_deg'][:] if 'equil_angle_deg' in dihedral_group else np.zeros(len(dihedral_ids))
                    dihedral_force = dihedral_group['spring_const'][:] if 'spring_const' in dihedral_group else np.zeros(len(dihedral_ids))
                    
                    f_out.write(f"DIHEDRALS ({len(dihedral_ids)} total):\n")
                    f_out.write("-" * 20 + "\n")
                    f_out.write(f"{'ID':>4} {'Atoms':>20} {'Angle':>10} {'Force':>12}\n")
                    f_out.write("-" * 50 + "\n")
                    
                    for i in range(min(20, len(dihedral_ids))):  # Show first 20
                        atoms = f"{dihedral_ids[i][0]:3d}-{dihedral_ids[i][1]:3d}-{dihedral_ids[i][2]:3d}-{dihedral_ids[i][3]:3d}"
                        f_out.write(f"{i:4d} {atoms:>20} {dihedral_equil[i]:8.1f}° {dihedral_force[i]:10.3f}\n")
                    
                    if len(dihedral_ids) > 20:
                        f_out.write(f"... and {len(dihedral_ids) - 20} more dihedrals\n")
                    f_out.write("\n")
            
            # ====================================================================
            # 4. NON-BONDED INTERACTIONS
            # ====================================================================
            f_out.write("4. NON-BONDED INTERACTIONS\n")
            f_out.write("-" * 40 + "\n")
            
            # Check for non-bonded interactions in both possible locations
            pairs = None
            coeffs = None
            
            if 'input/potential/pair' in h5:
                pair_group = h5['input/potential/pair']
                if 'pairs' in pair_group and 'coefficients' in pair_group:
                    pairs = pair_group['pairs'][:]
                    coeffs = pair_group['coefficients'][:]
            elif 'input/potential/martini_potential' in h5:
                martini_group = h5['input/potential/martini_potential']
                if 'pairs' in martini_group and 'coefficients' in martini_group:
                    pairs = martini_group['pairs'][:]
                    coeffs = martini_group['coefficients'][:]
            
            if pairs is not None and coeffs is not None:
                f_out.write(f"NON-BONDED PAIRS ({len(pairs)} total):\n")
                f_out.write("-" * 25 + "\n")
                f_out.write(f"{'ID':>4} {'Atoms':>10} {'Epsilon':>10} {'Sigma':>8} {'Q1':>8} {'Q2':>8}\n")
                f_out.write("-" * 55 + "\n")
                
                for i in range(min(20, len(pairs))):  # Show first 20
                    atom_pair = f"{pairs[i][0]:3d}-{pairs[i][1]:3d}"
                    epsilon, sigma, q1, q2 = coeffs[i]
                    f_out.write(f"{i:4d} {atom_pair:>10} {epsilon:8.3f} {sigma:8.3f} {q1:8.3f} {q2:8.3f}\n")
                
                if len(pairs) > 20:
                    f_out.write(f"... and {len(pairs) - 20} more pairs\n")
                f_out.write("\n")
            
            # ====================================================================
            # 5. SPLINE TABLES AND INTERACTION PARAMETERS
            # ====================================================================
            f_out.write("5. SPLINE TABLES AND INTERACTION PARAMETERS\n")
            f_out.write("-" * 40 + "\n")
            
            f_out.write("Note: Spline tables are generated internally by UPSIDE during simulation.\n")
            f_out.write("The actual spline data is not stored in the input file but computed on-the-fly.\n")
            f_out.write("Below are the parameters that would be used to generate splines:\n\n")
            
            # Look for spline tables in the potential group (in case they exist)
            spline_found = False
            if 'input/potential' in h5:
                potential_group = h5['input/potential']
                
                for key in potential_group.keys():
                    if 'spline' in key.lower() or 'table' in key.lower():
                        spline_found = True
                        f_out.write(f"\nSPLINE TABLE: {key}\n")
                        f_out.write("-" * 20 + "\n")
                        
                        try:
                            table_data = potential_group[key][:]
                            if len(table_data.shape) == 2:
                                # Show first few and last few rows
                                f_out.write(f"Table shape: {table_data.shape}\n")
                                f_out.write("First 5 rows:\n")
                                for i in range(min(5, len(table_data))):
                                    f_out.write(f"  Row {i}: {table_data[i]}\n")
                                
                                if len(table_data) > 10:
                                    f_out.write("...\n")
                                    f_out.write("Last 5 rows:\n")
                                    for i in range(max(5, len(table_data) - 5), len(table_data)):
                                        f_out.write(f"  Row {i}: {table_data[i]}\n")
                                elif len(table_data) > 5:
                                    f_out.write("Remaining rows:\n")
                                    for i in range(5, len(table_data)):
                                        f_out.write(f"  Row {i}: {table_data[i]}\n")
                            else:
                                f_out.write(f"1D array with {len(table_data)} elements\n")
                                f_out.write(f"First 10: {table_data[:10]}\n")
                                if len(table_data) > 10:
                                    f_out.write(f"Last 10: {table_data[-10:]}\n")
                        except Exception as e:
                            f_out.write(f"Error reading table: {e}\n")
            
            if not spline_found:
                f_out.write("No pre-computed spline tables found in input file.\n\n")
            
            # ====================================================================
            # 6. EXTERNAL SPLINE FILES (if they exist)
            # ====================================================================
            f_out.write("6. EXTERNAL SPLINE FILES\n")
            f_out.write("-" * 40 + "\n")
            
            # Check for spline files that might have been generated by UPSIDE
            spline_files = []
            for filename in ['all_splines.txt', 'bond_splines.txt', 'angle_splines.txt', 'dihedral_splines.txt']:
                spline_path = os.path.join(spline_dir, filename)
                if os.path.exists(spline_path):
                    spline_files.append(spline_path)
            
            if spline_files:
                f_out.write(f"Found {len(spline_files)} spline files generated by UPSIDE simulation:\n")
                for filename in spline_files:
                    f_out.write(f"  - {filename}\n")
                f_out.write("\n")
                
                # Read and display content from each spline file
                for spline_path in spline_files:
                    filename = os.path.basename(spline_path)
                    f_out.write(f"SPLINE FILE: {filename}\n")
                    f_out.write("-" * 30 + "\n")
                    
                    try:
                        with open(spline_path, 'r') as spline_file:
                            lines = spline_file.readlines()
                            
                        f_out.write(f"File size: {len(lines)} lines\n")
                        
                        # Show first 20 lines
                        f_out.write("First 20 lines:\n")
                        for i, line in enumerate(lines[:20]):
                            f_out.write(f"  {i+1:3d}: {line.rstrip()}\n")
                        
                        if len(lines) > 20:
                            f_out.write(f"  ... ({len(lines) - 20} more lines)\n")
                            
                            # Show last 10 lines if file is long
                            if len(lines) > 30:
                                f_out.write("Last 10 lines:\n")
                                for i, line in enumerate(lines[-10:], len(lines) - 9):
                                    f_out.write(f"  {i:3d}: {line.rstrip()}\n")
                        
                        f_out.write("\n")
                        
                    except Exception as e:
                        f_out.write(f"Error reading spline file {filename}: {e}\n\n")
            else:
                f_out.write("No external spline files found.\n")
                f_out.write("Spline files are generated by UPSIDE during simulation runs.\n")
                f_out.write("Run a simulation to generate spline tables.\n\n")
            
            # Show interaction parameters that would be used for spline generation
            f_out.write("INTERACTION PARAMETERS (for spline generation):\n")
            f_out.write("-" * 40 + "\n")
            
            # Non-bonded parameters
            if 'input/potential/martini_potential' in h5:
                martini_group = h5['input/potential/martini_potential']
                if 'coefficients' in martini_group:
                    coeffs = martini_group['coefficients'][:]
                    f_out.write("MARTINI NON-BONDED PARAMETERS:\n")
                    f_out.write("Format: [epsilon, sigma, q1, q2] for each pair\n")
                    f_out.write("First 10 coefficient sets:\n")
                    for i in range(min(10, len(coeffs))):
                        epsilon, sigma, q1, q2 = coeffs[i]
                        f_out.write(f"  Pair {i}: ε={epsilon:.6f}, σ={sigma:.6f}, q1={q1:.3f}, q2={q2:.3f}\n")
                    if len(coeffs) > 10:
                        f_out.write(f"  ... and {len(coeffs) - 10} more coefficient sets\n")
                    f_out.write("\n")
            
            # Bond parameters
            if 'input/potential/dist_spring' in h5:
                bond_group = h5['input/potential/dist_spring']
                if 'spring_const' in bond_group and 'equil_dist' in bond_group:
                    spring_consts = bond_group['spring_const'][:]
                    equil_dists = bond_group['equil_dist'][:]
                    f_out.write("BOND PARAMETERS:\n")
                    f_out.write("Format: [equilibrium_distance, force_constant] for each bond\n")
                    f_out.write("First 10 bond parameter sets:\n")
                    for i in range(min(10, len(spring_consts))):
                        f_out.write(f"  Bond {i}: r0={equil_dists[i]:.6f}Å, k={spring_consts[i]:.6f} E_up/Å²\n")
                    if len(spring_consts) > 10:
                        f_out.write(f"  ... and {len(spring_consts) - 10} more bond parameter sets\n")
                    f_out.write("\n")
            
            # Angle parameters
            if 'input/potential/angle_spring' in h5:
                angle_group = h5['input/potential/angle_spring']
                if 'spring_const' in angle_group and 'equil_angle' in angle_group:
                    spring_consts = angle_group['spring_const'][:]
                    equil_angles = angle_group['equil_angle'][:]
                    f_out.write("ANGLE PARAMETERS:\n")
                    f_out.write("Format: [equilibrium_angle, force_constant] for each angle\n")
                    f_out.write("First 10 angle parameter sets:\n")
                    for i in range(min(10, len(spring_consts))):
                        f_out.write(f"  Angle {i}: θ0={equil_angles[i]:.6f}°, k={spring_consts[i]:.6f} E_up/deg²\n")
                    if len(spring_consts) > 10:
                        f_out.write(f"  ... and {len(spring_consts) - 10} more angle parameter sets\n")
                    f_out.write("\n")
            
            # Dihedral parameters
            if 'input/potential/dihedral_spring' in h5:
                dihedral_group = h5['input/potential/dihedral_spring']
                if 'spring_const' in dihedral_group and 'equil_angle_deg' in dihedral_group:
                    spring_consts = dihedral_group['spring_const'][:]
                    equil_angles = dihedral_group['equil_angle_deg'][:]
                    f_out.write("DIHEDRAL PARAMETERS:\n")
                    f_out.write("Format: [equilibrium_angle, force_constant] for each dihedral\n")
                    f_out.write("First 10 dihedral parameter sets:\n")
                    for i in range(min(10, len(spring_consts))):
                        f_out.write(f"  Dihedral {i}: φ0={equil_angles[i]:.6f}°, k={spring_consts[i]:.6f} E_up\n")
                    if len(spring_consts) > 10:
                        f_out.write(f"  ... and {len(spring_consts) - 10} more dihedral parameter sets\n")
                    f_out.write("\n")
            
            # ====================================================================
            # 7. INTERACTION SUMMARY
            # ====================================================================
            f_out.write("\n7. INTERACTION SUMMARY\n")
            f_out.write("-" * 40 + "\n")
            
            # Count different types of interactions
            bond_count = len(h5['input/potential/dist_spring']['id']) if 'input/potential/dist_spring' in h5 and 'id' in h5['input/potential/dist_spring'] else 0
            angle_count = len(h5['input/potential/angle_spring']['id']) if 'input/potential/angle_spring' in h5 and 'id' in h5['input/potential/angle_spring'] else 0
            dihedral_count = len(h5['input/potential/dihedral_spring']['id']) if 'input/potential/dihedral_spring' in h5 and 'id' in h5['input/potential/dihedral_spring'] else 0
            pair_count = len(h5['input/potential/pair']['pairs']) if 'input/potential/pair' in h5 and 'pairs' in h5['input/potential/pair'] else 0
            
            f_out.write(f"Total bonds: {bond_count}\n")
            f_out.write(f"Total angles: {angle_count}\n")
            f_out.write(f"Total dihedrals: {dihedral_count}\n")
            f_out.write(f"Total non-bonded pairs: {pair_count}\n")
            
            # Charge summary
            charges = None
            if 'input/charges' in h5:
                charges = h5['input/charges'][:]
            elif 'input/potential/martini_potential/charges' in h5:
                charges = h5['input/potential/martini_potential/charges'][:]
            
            if charges is not None:
                total_charge = np.sum(charges)
                charged_particles = np.sum(charges != 0)
                f_out.write(f"\nCharge summary:\n")
                f_out.write(f"Total charge: {total_charge:.6f}\n")
                f_out.write(f"Charged particles: {charged_particles}/{len(charges)}\n")
                
                # Show unique charges
                unique_charges = np.unique(charges)
                f_out.write(f"Unique charges: {unique_charges}\n")
            
            f_out.write("\n" + "=" * 80 + "\n")
            f_out.write("VISUALIZATION COMPLETE\n")
            f_out.write("=" * 80 + "\n")
    
    print(f"Comprehensive HDF5 visualization report saved to: {output_file}")
    return output_file

def main():
    """Main function"""
    if len(sys.argv) < 2:
        print("Usage: python visualize_h5_input.py <h5_file> [output_file] [spline_dir]")
        print("Example: python visualize_h5_input.py outputs/martini_test/test.input.up")
        print("Example: python visualize_h5_input.py outputs/martini_test/test.input.up report.txt bak/")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    spline_dir = sys.argv[3] if len(sys.argv) > 3 else "."
    
    if not os.path.exists(h5_file):
        print(f"Error: HDF5 file '{h5_file}' not found!")
        sys.exit(1)
    
    try:
        output_file = analyze_h5_file(h5_file, output_file, spline_dir)
        print(f"✅ Visualization complete! Report saved to: {output_file}")
    except Exception as e:
        print(f"❌ Error during visualization: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
