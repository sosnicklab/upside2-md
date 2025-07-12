#!/usr/bin/env python3

import sys, os
import numpy as np
import h5py
from datetime import datetime

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

def write_vtf_frame(f, pos, frame_num):
    """Write a single frame in VTF format"""
    f.write("\ntimestep ordered\n")
    for i in range(len(pos)):
        x, y, z = pos[i]
        # Use coordinates as-is (already in Angstroms)
        f.write(f"{x:.3f} {y:.3f} {z:.3f}\n")

def write_pdb_frame(f, pos, frame_num, atom_types, residue_ids, x_len, y_len, z_len):
    """Write a single frame in PDB format with proper atom names and residue info"""
    f.write(f"MODEL     {frame_num+1:4d}\n")
    # Add CRYST1 record for each frame with current box dimensions
    f.write(f"CRYST1{x_len:9.3f}{y_len:9.3f}{z_len:9.3f}  90.00  90.00  90.00 P 1           1\n")
    
    for i in range(len(pos)):
        x, y, z = pos[i]
        mtype = atom_types[i]
        resid = residue_ids[i] if i < len(residue_ids) else 1
        
        # Map MARTINI types to PDB atom names
        if mtype == 'Q0':  # NC3 (choline)
            atom_name = 'NC3'
            resname = 'DOPC'
        elif mtype == 'Qa':  # PO4 (phosphate)
            atom_name = 'PO4'
            resname = 'DOPC'
        elif mtype == 'Na':  # GL1/GL2 (glycerol)
            atom_name = 'GL1' if i % 12 < 6 else 'GL2'
            resname = 'DOPC'
        elif mtype == 'C1':  # Carbon beads
            if 'C1A' in str(i) or i % 12 == 4:
                atom_name = 'C1A'
            elif 'C1B' in str(i) or i % 12 == 8:
                atom_name = 'C1B'
            elif 'C3A' in str(i) or i % 12 == 6:
                atom_name = 'C3A'
            elif 'C4A' in str(i) or i % 12 == 7:
                atom_name = 'C4A'
            elif 'C3B' in str(i) or i % 12 == 10:
                atom_name = 'C3B'
            elif 'C4B' in str(i) or i % 12 == 11:
                atom_name = 'C4B'
            else:
                atom_name = 'C1'
            resname = 'DOPC'
        elif mtype == 'C3':  # D2A/D2B (double bond)
            if 'D2A' in str(i) or i % 12 == 5:
                atom_name = 'D2A'
            elif 'D2B' in str(i) or i % 12 == 9:
                atom_name = 'D2B'
            else:
                atom_name = 'D2'
            resname = 'DOPC'
        elif mtype == 'P4':  # Water
            atom_name = 'W'
            resname = 'WAT'
        elif mtype == 'Qd':  # Sodium
            atom_name = 'NA'
            resname = 'NA'
        elif mtype == 'Qa' and resid > 200:  # Chloride
            atom_name = 'CL'
            resname = 'CL'
        else:
            atom_name = mtype
            resname = 'UNK'
        
        f.write(f"ATOM  {i+1:5d} {atom_name:>4} {resname:3} A{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n")
    f.write("ENDMDL\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_martini_vtf.py <input.up> <output.vtf/pdb>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    output_format = output_file.split('.')[-1].lower()

    if output_format not in ['vtf', 'pdb']:
        print("Error: Output file must have .vtf or .pdb extension")
        sys.exit(1)

    # Open the HDF5 file
    with h5py.File(input_file, 'r') as t:
        # Get the number of frames
        n_frame = len(t['output/pos'])
        print(f"Number of frames: {n_frame}")
        
        # Get position data shape
        pos = t['output/pos'][0]
        n_particles = pos.shape[1]
        print(f"Number of particles: {n_particles}")
        print(f"Position data shape: {pos.shape}")
        
        # Print debug info for first frame
        print("\nFirst frame statistics:")
        print(f"Min position: {np.min(pos)}")
        print(f"Max position: {np.max(pos)}")
        print(f"Mean position: {np.mean(pos)}")
        print(f"Std position: {np.std(pos)}")
        
        # Read box dimensions from HDF5 file
        # Check different potential groups for box dimensions
        x_len = y_len = z_len = None
        
        # Try to read from martini_potential group
        if 'input/potential/martini_potential' in t:
            potential_group = t['input/potential/martini_potential']
            if 'x_len' in potential_group.attrs:
                x_len = potential_group.attrs['x_len']
                y_len = potential_group.attrs['y_len'] 
                z_len = potential_group.attrs['z_len']
                print(f"Found box dimensions from martini_potential: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
        
        # Try to read from periodic_boundary_potential group as fallback
        if x_len is None and 'input/potential/periodic_boundary_potential' in t:
            pbc_group = t['input/potential/periodic_boundary_potential']
            if 'x_len' in pbc_group.attrs:
                x_len = pbc_group.attrs['x_len']
                y_len = pbc_group.attrs['y_len']
                z_len = pbc_group.attrs['z_len']
                print(f"Found box dimensions from periodic_boundary_potential: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
        
        # Use default values if not found
        if x_len is None:
            x_len = y_len = z_len = 50.0
            print(f"WARNING: No box dimensions found in HDF5 file! Using defaults: {x_len:.1f} x {y_len:.1f} x {z_len:.1f} Å")
        
        # Read atom types and residue IDs from HDF5 input
        atom_types = t['input/type'][:].astype(str)
        residue_ids = t['input/residue_ids'][:] if 'residue_ids' in t['input'] else np.ones(n_particles, dtype=int)
        
        print(f"Atom types found: {set(atom_types)}")
        print(f"Number of unique residues: {len(set(residue_ids))}")
        print(f"Residue ID range: {min(residue_ids)} to {max(residue_ids)}")
        
        # Count DOPC atoms by residue
        dopc_residues = set()
        for i, mtype in enumerate(atom_types):
            if mtype in ['Q0', 'Qa', 'Na', 'C1', 'C3']:  # DOPC bead types
                dopc_residues.add(residue_ids[i])
        
        print(f"DOPC residues found: {sorted(dopc_residues)}")
        print(f"Number of DOPC residues: {len(dopc_residues)}")
        
        # MARTINI to VTF mapping with colors and radii
        martini_to_vtf = {
            # DOPC lipid beads
            'Q0': 'NC3',   # Choline (blue)
            'Qa': 'PO4',   # Phosphate (red)
            'Na': 'GL',    # Glycerol (green)
            'C1': 'C1',    # Carbon (gray)
            'C3': 'C3',    # Double bond carbon (yellow)
            # Water and ions
            'P4': 'W',     # Water (cyan)
            'Qd': 'NA',    # Sodium (orange)
        }
        
        # Color mapping for VTF
        color_map = {
            'NC3': 'blue',
            'PO4': 'red', 
            'GL': 'green',
            'C1': 'gray',
            'C3': 'yellow',
            'W': 'cyan',
            'NA': 'orange',
            'CL': 'purple'
        }
        
        # Radius mapping for VTF
        radius_map = {
            'NC3': 2.5,   # Choline head
            'PO4': 2.5,   # Phosphate head
            'GL': 2.0,    # Glycerol
            'C1': 1.8,    # Carbon
            'C3': 1.8,    # Double bond carbon
            'W': 1.5,     # Water
            'NA': 2.0,    # Sodium
            'CL': 2.0     # Chloride
        }
        
        # Create output file
        with open(output_file, 'w') as f:
            if output_format == 'vtf':
                # Write VTF structure section with atom types, colors, and radii
                f.write("# VTF file generated from MARTINI bilayer simulation\n")
                f.write("# Structure section\n")
                
                for i in range(n_particles):
                    mtype = atom_types[i]
                    vtf_name = martini_to_vtf.get(mtype, mtype)
                    color = color_map.get(vtf_name, 'white')
                    radius = radius_map.get(vtf_name, 2.0)
                    resid = residue_ids[i] if i < len(residue_ids) else 1
                    
                    # Determine residue name based on atom type FIRST, then residue ID
                    # This prevents water molecules (P4) from being misidentified as lipids due to residue ID overlap
                    if mtype == 'P4':  # Water - check atom type FIRST
                        resname = 'WAT'
                    elif resid in dopc_residues:  # DOPC lipids
                        resname = 'DOPC'
                    else:  # Ions
                        resname = 'ION'
                    
                    f.write(f"atom {i} name {vtf_name} type {vtf_name} radius {radius:.1f} resid {resid} resname {resname} segid s0\n")
                
                # Add periodic boundary conditions to VTF
                # VTF format: pbc <a> <b> <c> [alpha] [beta] [gamma]
                # For orthorhombic box: alpha=beta=gamma=90.0
                f.write(f"\npbc {x_len:.3f} {y_len:.3f} {z_len:.3f} 90.0 90.0 90.0\n")
                f.write("# Periodic boundary conditions enabled\n")
                
                # Add bonds for DOPC lipids
                f.write("\n# Bonds section\n")
                bond_count = 0
                for resid in dopc_residues:
                    # Find atoms belonging to this lipid
                    lipid_atoms = np.where(residue_ids == resid)[0]
                    # Sort atoms by their index to ensure proper order
                    lipid_atoms = sorted(lipid_atoms)
                    
                    if len(lipid_atoms) >= 12:  # DOPC has 12 beads
                        # Take first 12 atoms if there are more
                        lipid_atoms = lipid_atoms[:12]
                        
                        # DOPC bonds: NC3-PO4, PO4-GL1, GL1-GL2, GL1-C1A, C1A-D2A, D2A-C3A, C3A-C4A, GL2-C1B, C1B-D2B, D2B-C3B, C3B-C4B
                        bonds = [
                            (lipid_atoms[0], lipid_atoms[1]),   # NC3-PO4
                            (lipid_atoms[1], lipid_atoms[2]),   # PO4-GL1
                            (lipid_atoms[2], lipid_atoms[3]),   # GL1-GL2
                            (lipid_atoms[2], lipid_atoms[4]),   # GL1-C1A
                            (lipid_atoms[4], lipid_atoms[5]),   # C1A-D2A
                            (lipid_atoms[5], lipid_atoms[6]),   # D2A-C3A
                            (lipid_atoms[6], lipid_atoms[7]),   # C3A-C4A
                            (lipid_atoms[3], lipid_atoms[8]),   # GL2-C1B
                            (lipid_atoms[8], lipid_atoms[9]),   # C1B-D2B
                            (lipid_atoms[9], lipid_atoms[10]),  # D2B-C3B
                            (lipid_atoms[10], lipid_atoms[11])  # C3B-C4B
                        ]
                        for atom1, atom2 in bonds:
                            f.write(f"bond {atom1}:{atom2}\n")
                            bond_count += 1
                    else:
                        print(f"Warning: Lipid residue {resid} has {len(lipid_atoms)} atoms, expected 12")
                
                print(f"Added {bond_count} bonds for lipid molecules")
                
            else:  # PDB format
                # Write PDB header
                f.write("TITLE     MARTINI LIPID BILAYER SIMULATION\n")
                f.write("REMARK    GENERATED BY UPSIDE\n")
                f.write(f"REMARK    DATE: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("REMARK    CONTAINS DOPC LIPIDS, WATER, AND IONS\n")
                f.write(f"REMARK    BOX DIMENSIONS: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Angstroms\n")
            
            # Write frames
            for frame in range(n_frame):
                pos = t['output/pos'][frame]
                
                # Reshape position data to (n_particles, 3)
                frame_pos = pos[0].reshape(n_particles, 3)
                
                # Check for NaN values and replace with last valid frame
                if np.isnan(frame_pos).any():
                    if frame > 0:
                        valid_pos = t['output/pos'][frame-1][0].reshape(n_particles, 3)
                        frame_pos = np.where(np.isnan(frame_pos), valid_pos, frame_pos)
                    else:
                        print(f"Warning: NaN values in first frame, replacing with zeros")
                        frame_pos = np.where(np.isnan(frame_pos), 0.0, frame_pos)
                
                # Write frame in appropriate format
                if output_format == 'vtf':
                    write_vtf_frame(f, frame_pos, frame)
                else:  # PDB format
                    write_pdb_frame(f, frame_pos, frame, atom_types, residue_ids, x_len, y_len, z_len)
                
                # Print progress every 100 frames
                if frame % 100 == 0:
                    print(f"Processed frame {frame}/{n_frame}")

if __name__ == "__main__":
    main() 