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

def write_pdb_frame(f, pos, frame_num, atom_types, residue_ids, x_len, y_len, z_len, pdb_atom_names=None, pdb_residue_names=None, pdb_residue_ids=None):
    """Write a single frame in PDB format with proper atom names and residue info"""
    f.write(f"MODEL     {frame_num+1:4d}\n")
    # Add CRYST1 record for each frame with current box dimensions
    f.write(f"CRYST1{x_len:9.3f}{y_len:9.3f}{z_len:9.3f}  90.00  90.00  90.00 P 1           1\n")
    
    for i in range(len(pos)):
        x, y, z = pos[i]
        mtype = atom_types[i]
        resid = residue_ids[i] if i < len(residue_ids) else 1
        
        # Use PDB information if available, otherwise fall back to MARTINI mapping
        if pdb_atom_names is not None and i < len(pdb_atom_names):
            # Use original PDB atom names and residue information
            atom_name = pdb_atom_names[i]
            resname = pdb_residue_names[i]
            resid = pdb_residue_ids[i]
        else:
            # Fallback to MARTINI type mapping
            if mtype == 'Q0':  # NC3 (choline)
                atom_name = 'NC3'
                resname = 'DOPC'
            elif mtype == 'Qa':  # PO4 (phosphate) or CL (chloride)
                if resid > 200:  # Chloride (high residue IDs)
                    atom_name = 'CL'
                    resname = 'CL'
                else:  # Phosphate in DOPC
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

    # Read original PDB file to get proper atom names and residue information
    pdb_file = 'input.pdb'
    if not os.path.exists(pdb_file):
        print(f"Warning: Original PDB file '{pdb_file}' not found. Using MARTINI type mapping.")
        pdb_atom_names = None
        pdb_residue_names = None
        pdb_residue_ids = None
    else:
        print(f"Reading original PDB file: {pdb_file}")
        pdb_atom_names = []
        pdb_residue_names = []
        pdb_residue_ids = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    atom_name = line[12:16].strip()
                    residue_name = line[17:21].strip()
                    residue_id = int(line[22:26])
                    pdb_atom_names.append(atom_name)
                    pdb_residue_names.append(residue_name)
                    pdb_residue_ids.append(residue_id)
        
        print(f"Read {len(pdb_atom_names)} atoms from PDB file")
        print(f"Residue types found: {set(pdb_residue_names)}")
        print(f"Residue ID range: {min(pdb_residue_ids)} to {max(pdb_residue_ids)}")

    # Open the HDF5 file
    with h5py.File(input_file, 'r') as t:
        # Get the number of frames
        n_frame = len(t['output/pos'])
        print(f"Number of frames: {n_frame}")
        
        # Get position data shape from input (first frame should use input positions)
        input_pos = t['input/pos'][:]
        if len(input_pos.shape) == 3:  # (n_atoms, 3, n_frames)
            input_pos = input_pos[:, :, 0]  # Take first frame
        n_particles = input_pos.shape[0]
        print(f"Number of particles: {n_particles}")
        print(f"Input position data shape: {input_pos.shape}")
        
        # Print debug info for input frame
        print("\nInput frame statistics:")
        print(f"Min position: {np.min(input_pos)}")
        print(f"Max position: {np.max(input_pos)}")
        print(f"Mean position: {np.mean(input_pos)}")
        print(f"Std position: {np.std(input_pos)}")
        
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
        
        # Debug: Print first 10 atom types and Z positions as read from .up
        print("\n[DEBUG] First 10 atoms in VTF input (.up):")
        for i in range(10):
            print(f"  Atom {i}: type={atom_types[i]}, Z={input_pos[i,2]:.2f} Å")
        
        print(f"Atom types found: {set(atom_types)}")
        print(f"Number of unique residues: {len(set(residue_ids))}")
        print(f"Residue ID range: {min(residue_ids)} to {max(residue_ids)}")
        
        # Count DOPC atoms by residue
        dopc_residues = set()
        if pdb_residue_names is not None:
            # Use PDB residue names to identify DOPC
            for i, resname in enumerate(pdb_residue_names):
                if resname == 'DOPC':
                    dopc_residues.add(pdb_residue_ids[i])
        else:
            # Fallback to MARTINI type-based identification
            for i, mtype in enumerate(atom_types):
                if mtype in ['Q0', 'Qa', 'Na', 'C1', 'C3']:  # DOPC bead types
                    dopc_residues.add(residue_ids[i])
        
        print(f"DOPC residues found: {sorted(dopc_residues)}")
        print(f"Number of DOPC residues: {len(dopc_residues)}")
        
        # MARTINI to VTF mapping with colors and radii
        martini_to_vtf = {
            # DOPC lipid beads
            'Q0': 'NC3',   # Choline (blue)
            'Qa': 'PO4',   # Phosphate (red) - will be overridden for chloride
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
                    resid = residue_ids[i] if i < len(residue_ids) else 1
                    
                    # Use MARTINI type mapping from HDF5 file (correct approach)
                    if mtype == 'P4':  # Water
                        vtf_name = 'W'
                        resname = 'WAT'
                    elif mtype == 'Qd':  # Sodium
                        vtf_name = 'NA'
                        resname = 'NA'
                    elif mtype == 'Qa' and resid > 200:  # Chloride (high residue IDs)
                        vtf_name = 'CL'
                        resname = 'CL'
                    elif mtype == 'Qa' and resid in dopc_residues:  # Phosphate in DOPC
                        vtf_name = 'PO4'
                        resname = 'DOPC'
                    elif resid in dopc_residues:  # Other DOPC atoms
                        vtf_name = martini_to_vtf.get(mtype, mtype)
                        resname = 'DOPC'
                    else:  # Unknown ions or other molecules
                        vtf_name = martini_to_vtf.get(mtype, mtype)
                        resname = 'ION'
                    
                    # Map atom names to colors and radii
                    color = color_map.get(vtf_name, 'white')
                    radius = radius_map.get(vtf_name, 2.0)
                    
                    f.write(f"atom {i} name {vtf_name} type {vtf_name} radius {radius:.1f} resid {resid} resname {resname} segid s0\n")
                
                # Add periodic boundary conditions to VTF
                # VTF format: pbc <a> <b> <c> [alpha] [beta] [gamma]
                # For orthorhombic box: alpha=beta=gamma=90.0
                f.write(f"\npbc {x_len:.3f} {y_len:.3f} {z_len:.3f} 90.0 90.0 90.0\n")
                f.write("# Periodic boundary conditions enabled\n")
                
                # Add bonds for DOPC lipids
                f.write("\n# Bonds section\n")
                bond_count = 0
                
                # Use PDB residue IDs if available, otherwise use HDF5 residue IDs
                bond_residue_ids = pdb_residue_ids if pdb_residue_ids is not None else residue_ids
                
                # Ensure bond_residue_ids is a numpy array
                if not isinstance(bond_residue_ids, np.ndarray):
                    bond_residue_ids = np.array(bond_residue_ids)
                
                for resid in dopc_residues:
                    # Find atoms belonging to this lipid
                    lipid_atoms = np.where(bond_residue_ids == resid)[0]
                    
                    if len(lipid_atoms) >= 12:  # DOPC has 12 beads
                        # Sort atoms by their index to ensure proper order
                        lipid_atoms = sorted(lipid_atoms)
                        
                        # Take first 12 atoms if there are more
                        lipid_atoms = lipid_atoms[:12]
                        
                        # Verify we have the right atom types for DOPC
                        lipid_types = [atom_types[i] for i in lipid_atoms]
                        expected_types = ['Q0', 'Qa', 'Na', 'Na', 'C1', 'C3', 'C1', 'C1', 'C1', 'C3', 'C1', 'C1']
                        
                        # Check if types match (allowing for some flexibility)
                        type_matches = sum(1 for i, (actual, expected) in enumerate(zip(lipid_types, expected_types)) if actual == expected)
                        
                        if type_matches >= 8:  # At least 8 out of 12 types should match
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
                            print(f"Warning: Lipid residue {resid} has unexpected atom types: {lipid_types}")
                            print(f"  Expected: {expected_types}")
                            print(f"  Only {type_matches}/12 types match")
                    else:
                        print(f"Warning: Lipid residue {resid} has {len(lipid_atoms)} atoms, expected 12")
                
                print(f"Added {bond_count} bonds for lipid molecules")
                
                # Write frames starting with initial PDB structure
                for frame in range(n_frame + 1):  # +1 to include initial structure
                    if frame == 0:
                        # Frame 0: Initial PDB structure (before any simulation)
                        print(f"Writing initial PDB structure as frame 0...")
                        
                        # Read PDB positions and apply same processing as run_martini.py
                        pdb_positions = []
                        with open('input.pdb', 'r') as pdb_f:
                            for line in pdb_f:
                                if line.startswith('ATOM'):
                                    x = float(line[30:38])
                                    y = float(line[38:46])
                                    z = float(line[46:54])
                                    pdb_positions.append([x, y, z])
                        
                        pdb_positions = np.array(pdb_positions)
                        
                        # Shift initial PDB positions to centered box convention
                        center_shift = np.array([x_len/2, y_len/2, z_len/2])
                        frame_pos = pdb_positions - center_shift
                        # Debug print for frame 0
                        print(f"[DEBUG] Frame 0 (PDB after centering): mean={frame_pos.mean(axis=0)}, min={frame_pos.min(axis=0)}, max={frame_pos.max(axis=0)}")
                    else:
                        # Use output positions for simulation frames (frame-1 because we added initial frame)
                        pos = t['output/pos'][frame-1]
                        
                        # Reshape position data to (n_particles, 3)
                        frame_pos = pos[0].reshape(n_particles, 3)
                        
                        # Apply periodic boundary conditions to wrap all atoms into the box
                        frame_pos[:, 0] = frame_pos[:, 0] - x_len * np.floor(frame_pos[:, 0] / x_len)
                        frame_pos[:, 1] = frame_pos[:, 1] - y_len * np.floor(frame_pos[:, 1] / y_len)
                        frame_pos[:, 2] = frame_pos[:, 2] - z_len * np.floor(frame_pos[:, 2] / z_len)
                        # Shift to center at origin
                        frame_pos -= np.array([x_len/2, y_len/2, z_len/2])
                        
                        # Check for NaN values and replace with last valid frame
                        if np.isnan(frame_pos).any():
                            if frame > 1:  # frame > 1 because we added initial frame
                                valid_pos = t['output/pos'][frame-2][0].reshape(n_particles, 3)
                                frame_pos = np.where(np.isnan(frame_pos), valid_pos, frame_pos)
                            else:
                                print(f"Warning: NaN values in first simulation frame, replacing with zeros")
                                frame_pos = np.where(np.isnan(frame_pos), 0.0, frame_pos)
                    
                    # Write frame in VTF format
                    write_vtf_frame(f, frame_pos, frame)
                    
                    # Print progress every 100 frames
                    if frame % 100 == 0:
                        print(f"Processed frame {frame}/{n_frame}")
                
            else:  # PDB format
                # Write PDB header
                f.write("TITLE     MARTINI LIPID BILAYER SIMULATION\n")
                f.write("REMARK    GENERATED BY UPSIDE\n")
                f.write(f"REMARK    DATE: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("REMARK    CONTAINS DOPC LIPIDS, WATER, AND IONS\n")
                f.write(f"REMARK    BOX DIMENSIONS: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Angstroms\n")
                
                # Write frames for PDB format
                for frame in range(n_frame + 1):  # +1 to include initial structure
                    if frame == 0:
                        # Frame 0: Initial PDB structure (before any simulation)
                        print(f"Writing initial PDB structure as frame 0...")
                        
                        # Read PDB positions and apply same processing as run_martini.py
                        pdb_positions = []
                        with open('input.pdb', 'r') as pdb_f:
                            for line in pdb_f:
                                if line.startswith('ATOM'):
                                    x = float(line[30:38])
                                    y = float(line[38:46])
                                    z = float(line[46:54])
                                    pdb_positions.append([x, y, z])
                        
                        pdb_positions = np.array(pdb_positions)
                        
                        # Apply PBC wrapping (same as in run_martini.py)
                        wrapped_positions = pdb_positions.copy()
                        
                        # X dimension wrapping
                        wrapped_positions[:, 0] = wrapped_positions[:, 0] - x_len * np.floor(wrapped_positions[:, 0] / x_len)
                        # Y dimension wrapping  
                        wrapped_positions[:, 1] = wrapped_positions[:, 1] - y_len * np.floor(wrapped_positions[:, 1] / y_len)
                        
                        # Z dimension - simple shift (no wrapping)
                        z_min = np.min(pdb_positions[:, 2])
                        if z_min < 0:
                            z_shift = -z_min + 1.0
                            wrapped_positions[:, 2] = pdb_positions[:, 2] + z_shift
                        else:
                            wrapped_positions[:, 2] = pdb_positions[:, 2]
                        
                        frame_pos = wrapped_positions
                    else:
                        # Use output positions for simulation frames (frame-1 because we added initial frame)
                        pos = t['output/pos'][frame-1]
                        
                        # Reshape position data to (n_particles, 3)
                        frame_pos = pos[0].reshape(n_particles, 3)
                        
                        # Apply PBC wrapping (same as in run_martini.py)
                        frame_pos[:, 0] = frame_pos[:, 0] - x_len * np.floor(frame_pos[:, 0] / x_len)
                        frame_pos[:, 1] = frame_pos[:, 1] - y_len * np.floor(frame_pos[:, 1] / y_len)
                        frame_pos[:, 2] = frame_pos[:, 2] - z_len * np.floor(frame_pos[:, 2] / z_len)
                        # Shift to center at origin
                        frame_pos -= np.array([x_len/2, y_len/2, z_len/2])
                        
                        # Check for NaN values and replace with last valid frame
                        if np.isnan(frame_pos).any():
                            if frame > 1:  # frame > 1 because we added initial frame
                                valid_pos = t['output/pos'][frame-2][0].reshape(n_particles, 3)
                                frame_pos = np.where(np.isnan(frame_pos), valid_pos, frame_pos)
                            else:
                                print(f"Warning: NaN values in first simulation frame, replacing with zeros")
                                frame_pos = np.where(np.isnan(frame_pos), 0.0, frame_pos)
                    
                    # Write frame in PDB format
                    write_pdb_frame(f, frame_pos, frame, atom_types, residue_ids, x_len, y_len, z_len, pdb_atom_names, pdb_residue_names, pdb_residue_ids)
                    
                    # Print progress every 100 frames
                    if frame % 100 == 0:
                        print(f"Processed frame {frame}/{n_frame}")

if __name__ == "__main__":
    main() 