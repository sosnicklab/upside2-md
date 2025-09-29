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
            # If no PDB file available, use a simple fallback
            atom_name = mtype
            resname = 'UNK'
        
        f.write(f"ATOM  {i+1:5d} {atom_name:>4} {resname:3} A{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n")
    f.write("ENDMDL\n")

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("Usage: python extract_martini_vtf.py <input.up> <output.vtf/pdb> [input_file_for_structure] [pdb_id]")
        print("  If input_file_for_structure is not provided, will try to use input.up for structure data")
        print("  If pdb_id is not provided, will use default '1rkl'")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    input_file_for_structure = sys.argv[3] if len(sys.argv) >= 4 else input_file

    # Auto-detect PDB file based on input file name
    if len(sys.argv) >= 5:
        pdb_id = sys.argv[4]
    else:
        # Try to infer PDB ID from input file name
        input_basename = os.path.basename(input_file)
        if 'bilayer' in input_basename.lower():
            pdb_id = 'bilayer'
        else:
            pdb_id = '1rkl'  # fallback default
    output_format = output_file.split('.')[-1].lower()

    if output_format not in ['vtf', 'pdb']:
        print("Error: Output file must have .vtf or .pdb extension")
        sys.exit(1)

    # Read original PDB file to get proper atom names and residue information
    pdb_file = f'pdb/{pdb_id}.MARTINI.pdb'
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

    # Open the HDF5 files
    with h5py.File(input_file, 'r') as t, h5py.File(input_file_for_structure, 'r') as t_struct:
        # Count frames exclusively from PRODUCTION (regular) stage
        n_frame_regular = len(t['output/pos']) if 'output/pos' in t else 0
        n_frame_total = n_frame_regular

        print(f"Number of production (regular) frames: {n_frame_regular}")
        
        # Get position data shape from input (first frame should use input positions)
        input_pos = t_struct['input/pos'][:]
        print(f"Raw input position data shape: {input_pos.shape}")
        
        # Handle different position data formats
        if len(input_pos.shape) == 3:  # (n_atoms, 3, n_frames)
            input_pos = input_pos[:, :, 0]  # Take first frame
        elif len(input_pos.shape) == 2:  # (n_atoms, 3)
            input_pos = input_pos
        else:
            raise ValueError(f"Unexpected input position shape: {input_pos.shape}")
            
        n_particles = input_pos.shape[0]
        print(f"Number of particles: {n_particles}")
        print(f"Input position data shape: {input_pos.shape}")
        
        # Print debug info for input frame
        print("\nInput frame statistics:")
        print(f"Min position: {np.min(input_pos)}")
        print(f"Max position: {np.max(input_pos)}")
        print(f"Mean position: {np.mean(input_pos)}")
        print(f"Std position: {np.std(input_pos)}")
        
        # Read equilibrated box dimensions if present; otherwise fall back to input attrs
        x_len = y_len = z_len = None
        if 'equilibrated_box' in t:
            boxgrp = t['equilibrated_box']
            if all(k in boxgrp.attrs for k in ('x_len','y_len','z_len')):
                x_len = float(boxgrp.attrs['x_len'])
                y_len = float(boxgrp.attrs['y_len'])
                z_len = float(boxgrp.attrs['z_len'])
                print(f"Using equilibrated box: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")

        if x_len is None and 'input/potential/martini_potential' in t_struct:
            potential_group = t_struct['input/potential/martini_potential']
            if all(k in potential_group.attrs for k in ('x_len','y_len','z_len')):
                x_len = float(potential_group.attrs['x_len'])
                y_len = float(potential_group.attrs['y_len'])
                z_len = float(potential_group.attrs['z_len'])
                print(f"Found box dimensions from martini_potential: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")

        if x_len is None and 'input/potential/periodic_boundary_potential' in t_struct:
            pbc_group = t_struct['input/potential/periodic_boundary_potential']
            if all(k in pbc_group.attrs for k in ('x_len','y_len','z_len')):
                x_len = float(pbc_group.attrs['x_len'])
                y_len = float(pbc_group.attrs['y_len'])
                z_len = float(pbc_group.attrs['z_len'])
                print(f"Found box dimensions from periodic_boundary_potential: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")

        if x_len is None:
            x_len = y_len = z_len = 50.0
            print(f"WARNING: No box dimensions found in HDF5 file! Using defaults: {x_len:.1f} x {y_len:.1f} x {z_len:.1f} Å")
        
        # Read atom types and residue IDs from HDF5 input
        atom_types = t_struct['input/type'][:].astype(str)
        residue_ids = t_struct['input/residue_ids'][:] if 'residue_ids' in t_struct['input'] else np.ones(n_particles, dtype=int)
        
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
        
        # Use PDB atom names directly for VTF mapping
        # This ensures we use the actual topology from the input PDB file
        
        # Color mapping for VTF - use actual PDB atom names
        color_map = {}
        radius_map = {}
        
        # If we have PDB data, build color/radius maps from actual atom names
        if pdb_atom_names is not None:
            unique_atom_names = set(pdb_atom_names)
            for atom_name in unique_atom_names:
                if atom_name == 'NC3':
                    color_map[atom_name] = 'blue'
                    radius_map[atom_name] = 2.5
                elif atom_name == 'PO4':
                    color_map[atom_name] = 'red'
                    radius_map[atom_name] = 2.5
                elif atom_name in ['GL1', 'GL2']:
                    color_map[atom_name] = 'green'
                    radius_map[atom_name] = 2.0
                elif atom_name in ['C1A', 'C1B', 'C3A', 'C3B', 'C4A', 'C4B']:
                    color_map[atom_name] = 'gray'
                    radius_map[atom_name] = 1.8
                elif atom_name in ['D2A', 'D2B']:
                    color_map[atom_name] = 'yellow'
                    radius_map[atom_name] = 1.8
                elif atom_name == 'W':
                    color_map[atom_name] = 'cyan'
                    radius_map[atom_name] = 1.5
                elif atom_name == 'NA':
                    color_map[atom_name] = 'orange'
                    radius_map[atom_name] = 2.0
                elif atom_name == 'CL':
                    color_map[atom_name] = 'purple'
                    radius_map[atom_name] = 2.0
                else:
                    # Default for unknown atom types
                    color_map[atom_name] = 'white'
                    radius_map[atom_name] = 1.5
        else:
            # Fallback for when no PDB file is available
            color_map = {
                'Q0': 'blue', 'Qa': 'red', 'Na': 'green', 'C1': 'gray', 'C3': 'yellow',
                'P4': 'cyan', 'Qd': 'orange'
            }
            radius_map = {
                'Q0': 2.5, 'Qa': 2.5, 'Na': 2.0, 'C1': 1.8, 'C3': 1.8,
                'P4': 1.5, 'Qd': 2.0
            }
        
        # Create output file
        with open(output_file, 'w') as f:
            if output_format == 'vtf':
                # Write VTF structure section with atom types, colors, and radii
                f.write("# VTF file generated from MARTINI bilayer simulation\n")
                f.write("# Structure section\n")
                
                # Write atom types using PDB atom names if available, otherwise use MARTINI types
                for i in range(n_particles):
                    if pdb_atom_names is not None and i < len(pdb_atom_names):
                        atom_name = pdb_atom_names[i]
                    else:
                        atom_name = atom_types[i]
                    f.write(f"atom {i} name {atom_name}\n")
                
                # Write bonds if available
                if 'input/potential/dist_spring' in t:
                    bond_group = t['input/potential/dist_spring']
                    if 'id' in bond_group:
                        bonds = bond_group['id'][:]
                        for bond in bonds:
                            f.write(f"bond {bond[0]}:{bond[1]}\n")
                
                # Write box information
                f.write(f"pbc {x_len} {y_len} {z_len}\n")
                
                half_box = np.array([x_len/2, y_len/2, z_len/2])

                # Function to get frame positions from production (regular) stage only
                def get_frame_pos(frame_idx):
                    pos = t['output/pos'][frame_idx]

                    # Handle different output position data formats
                    if len(pos.shape) == 3:  # (n_frames, n_atoms, 3) or (n_atoms, 3, n_frames)
                        if pos.shape[0] == 1:  # (1, n_atoms, 3) - single frame
                            frame_pos = pos[0]  # Take first frame
                        elif pos.shape[2] == 3:  # (n_atoms, 3, n_frames)
                            frame_pos = pos[:, :, 0]  # Take first frame
                        else:
                            raise ValueError(f"Unexpected 3D position shape: {pos.shape}")
                    elif len(pos.shape) == 2:  # (n_atoms, 3)
                        frame_pos = pos
                    elif len(pos.shape) == 1:  # Flattened array
                        frame_pos = pos.reshape(n_particles, 3)
                    else:
                        raise ValueError(f"Unexpected output position shape: {pos.shape}")

                    return frame_pos

                # Write only production frames
                for frame in range(n_frame_total):
                    frame_pos = get_frame_pos(frame)

                    # Apply PBC wrapping to centered box convention [-L/2, L/2]
                    frame_pos = (frame_pos + half_box) % (2*half_box) - half_box

                    if np.isnan(frame_pos).any():
                        if frame > 0:
                            # Get previous valid frame
                            prev_frame_pos = get_frame_pos(frame - 1)
                            prev_frame_pos = (prev_frame_pos + half_box) % (2*half_box) - half_box
                            frame_pos = np.where(np.isnan(frame_pos), prev_frame_pos, frame_pos)
                        else:
                            print(f"Warning: NaN values in frame {frame}, replacing with zeros")
                            frame_pos = np.where(np.isnan(frame_pos), 0.0, frame_pos)

                    write_vtf_frame(f, frame_pos, frame)
                    if frame % 100 == 0:
                        print(f"Processed frame {frame}/{n_frame_total - 1}")

            else:  # PDB format
                # Write PDB header
                f.write("TITLE     MARTINI LIPID BILAYER SIMULATION\n")
                f.write("REMARK    GENERATED BY UPSIDE\n")
                f.write(f"REMARK    DATE: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("REMARK    CONTAINS DOPC LIPIDS, WATER, AND IONS\n")
                f.write(f"REMARK    BOX DIMENSIONS: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Angstroms\n")
                
                # Function to get frame positions from production (regular) stage only (for PDB)
                def get_frame_pos_pdb(frame_idx):
                    pos = t['output/pos'][frame_idx]

                    # Handle different output position data formats
                    if len(pos.shape) == 3:  # (n_frames, n_atoms, 3) or (n_atoms, 3, n_frames)
                        if pos.shape[0] == 1:  # (1, n_atoms, 3) - single frame
                            frame_pos = pos[0]  # Take first frame
                        elif pos.shape[2] == 3:  # (n_atoms, 3, n_frames)
                            frame_pos = pos[:, :, 0]  # Take first frame
                        else:
                            raise ValueError(f"Unexpected 3D position shape: {pos.shape}")
                    elif len(pos.shape) == 2:  # (n_atoms, 3)
                        frame_pos = pos
                    elif len(pos.shape) == 1:  # Flattened array
                        frame_pos = pos.reshape(n_particles, 3)
                    else:
                        raise ValueError(f"Unexpected output position shape: {pos.shape}")

                    return frame_pos

                # Write only production frames
                for frame in range(n_frame_total):
                    frame_pos = get_frame_pos_pdb(frame)

                    # Apply PBC wrapping to centered box convention [-L/2, L/2]
                    half_box = np.array([x_len/2, y_len/2, z_len/2])
                    frame_pos = (frame_pos + half_box) % (2*half_box) - half_box

                    # Check for NaN values and replace with last valid frame
                    if np.isnan(frame_pos).any():
                        if frame > 0:
                            # Get previous valid frame from the same stage
                            prev_frame_pos = get_frame_pos_pdb(frame - 1)
                            prev_frame_pos = (prev_frame_pos + half_box) % (2*half_box) - half_box
                            frame_pos = np.where(np.isnan(frame_pos), prev_frame_pos, frame_pos)
                        else:
                            print(f"Warning: NaN values in first simulation frame, replacing with zeros")
                            frame_pos = np.where(np.isnan(frame_pos), 0.0, frame_pos)

                    # Write frame in PDB format
                    write_pdb_frame(f, frame_pos, frame, atom_types, residue_ids, x_len, y_len, z_len, pdb_atom_names, pdb_residue_names, pdb_residue_ids)

                    # Print progress every 100 frames
                    if frame % 100 == 0:
                        print(f"Processed frame {frame}/{n_frame_total}")

if __name__ == "__main__":
    main() 