#!/usr/bin/env python3

import sys, os
import numpy as np
import h5py
import re
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

def centralize_system(frame_pos, atom_types, residue_ids, x_len, y_len, z_len):
    """
    Centralize the system by keeping the lipid bilayer centered in Z direction.
    Handle water molecules that flow over box boundaries by redistributing them evenly.
    """
    # Find lipid atoms (DOPC residues)
    lipid_mask = np.array([res_name == 'DOPC' for res_name in residue_ids])
    water_mask = np.array([res_name == 'W' for res_name in residue_ids])
    
    
    if np.any(lipid_mask):
        # Calculate the center of mass of lipids in Z direction
        lipid_z_positions = frame_pos[lipid_mask, 2]
        lipid_z_center = np.mean(lipid_z_positions)
        # Calculate the shift needed to center lipids at Z=0
        z_shift = -lipid_z_center
        
        # Apply shift to all atoms
        frame_pos[:, 2] += z_shift
        
        # Improve water distribution around the centered bilayer
        half_z = z_len / 2.0
        
        # Get all water molecule indices
        water_indices = np.where(water_mask)[0]
        n_water = len(water_indices)
        
        if n_water > 0:
            # Calculate the lipid bilayer thickness (approximate)
            lipid_z_positions = frame_pos[lipid_mask, 2]
            bilayer_thickness = np.max(lipid_z_positions) - np.min(lipid_z_positions)
            bilayer_margin = bilayer_thickness * 0.5  # Add some margin around bilayer
            
            # Define water regions: above and below the bilayer
            water_above_center = bilayer_margin
            water_below_center = -bilayer_margin
            
            # Split water molecules evenly between above and below the bilayer
            n_water_above = n_water // 2
            n_water_below = n_water - n_water_above
            
            # Redistribute water molecules
            if n_water_above > 0:
                # Place water above the bilayer
                above_indices = water_indices[:n_water_above]
                frame_pos[above_indices, 2] = np.random.uniform(water_above_center, half_z, n_water_above)
            
            if n_water_below > 0:
                # Place water below the bilayer
                below_indices = water_indices[n_water_above:]
                frame_pos[below_indices, 2] = np.random.uniform(-half_z, water_below_center, n_water_below)
        
        # Apply periodic boundary conditions to keep all atoms within the box
        frame_pos[:, 0] = ((frame_pos[:, 0] + x_len/2) % x_len) - x_len/2
        frame_pos[:, 1] = ((frame_pos[:, 1] + y_len/2) % y_len) - y_len/2
        frame_pos[:, 2] = ((frame_pos[:, 2] + z_len/2) % z_len) - z_len/2
        
    
    return frame_pos

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
    
    print(f"Extracting trajectory from: {input_file}")
    print(f"Output file: {output_file}")
    print(f"Structure source: {input_file_for_structure}")
    print(f"PDB ID: {pdb_id}")
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
        
        # First, try to read from simulation output (equilibrated box dimensions)
        if 'output/box' in t:
            box_data = t['output/box'][:]
            if len(box_data) > 0:
                # Get the last (most recent) box dimensions from the simulation
                last_box = box_data[-1]
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:  # (time, 3) format
                    x_len = float(last_box[0, 0])  # x dimension
                    y_len = float(last_box[0, 1])  # y dimension  
                    z_len = float(last_box[0, 2])  # z dimension
                    # Check if box dimensions are valid (non-zero)
                    if x_len > 0 and y_len > 0 and z_len > 0:
                        print(f"Using equilibrated box from simulation output: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
                    else:
                        print(f"Invalid box dimensions from simulation output: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å (likely NVT simulation)")
                        x_len = y_len = z_len = None
                elif len(last_box.shape) == 1 and len(last_box) == 3:  # (3,) format
                    x_len = float(last_box[0])
                    y_len = float(last_box[1])
                    z_len = float(last_box[2])
                    # Check if box dimensions are valid (non-zero)
                    if x_len > 0 and y_len > 0 and z_len > 0:
                        print(f"Using equilibrated box from simulation output: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
                    else:
                        print(f"Invalid box dimensions from simulation output: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å (likely NVT simulation)")
                        x_len = y_len = z_len = None
        
        # If no box data in HDF5, try to extract from log file
        if x_len is None:
            # Look for log file in the same directory as the input file
            log_file = input_file.replace('.up', '.log')
            if not os.path.exists(log_file):
                # Try alternative log file names
                log_file = input_file.replace('inputs/', 'outputs/martini_test/').replace('.up', '.log')
                if not os.path.exists(log_file):
                    log_file = input_file.replace('inputs/', 'outputs/martini_test/').replace('.up', '.optimized.log')
                    if not os.path.exists(log_file):
                        # Try the test.run.optimized.log pattern
                        log_file = input_file.replace('inputs/', 'outputs/martini_test/').replace('bilayer.up', 'test.run.optimized.log')
            
            if os.path.exists(log_file):
                print(f"Reading box dimensions from log file: {log_file}")
                try:
                    with open(log_file, 'r') as f:
                        lines = f.readlines()
                    
                    # Look for the last NPT box line in the log
                    last_box_line = None
                    for line in reversed(lines):
                        if '[NPT]' in line and 'box' in line:
                            last_box_line = line
                            break
                    
                    if last_box_line:
                        # Parse box dimensions from log line like:
                        # [NPT] t 1980.000 scale_xy 1.0000 scale_z 0.9965 | Pxy -4.147e-04 tgt 2.066e-05, Pz 1.198e-02 tgt 2.066e-05 | box 28.36 28.36 70.56
                        box_match = re.search(r'box\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', last_box_line)
                        if box_match:
                            x_len = float(box_match.group(1))
                            y_len = float(box_match.group(2))
                            z_len = float(box_match.group(3))
                            print(f"Using equilibrated box from log file: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
                except Exception as e:
                    print(f"Warning: Could not read box dimensions from log file: {e}")
        
        # Fallback to equilibrated_box group if present
        if x_len is None and 'equilibrated_box' in t:
            boxgrp = t['equilibrated_box']
            if all(k in boxgrp.attrs for k in ('x_len','y_len','z_len')):
                x_len = float(boxgrp.attrs['x_len'])
                y_len = float(boxgrp.attrs['y_len'])
                z_len = float(boxgrp.attrs['z_len'])
                print(f"Using equilibrated box: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")

        # Fallback to input potential attributes
        if x_len is None and 'input/potential/martini_potential' in t_struct:
            potential_group = t_struct['input/potential/martini_potential']
            if all(k in potential_group.attrs for k in ('x_len','y_len','z_len')):
                x_len = float(potential_group.attrs['x_len'])
                y_len = float(potential_group.attrs['y_len'])
                z_len = float(potential_group.attrs['z_len'])
                print(f"Found box dimensions from martini_potential: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")

        # Periodic boundary potential removed - using NVT ensemble without boundaries

        # Final fallback: try to read from PDB file
        if x_len is None and pdb_file is not None:
            try:
                with open(pdb_file, 'r') as f:
                    for line in f:
                        if line.startswith('CRYST1'):
                            parts = line.split()
                            if len(parts) >= 4:
                                x_len = float(parts[1])
                                y_len = float(parts[2])
                                z_len = float(parts[3])
                                print(f"Using original box dimensions from PDB file: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å")
                                break
            except Exception as e:
                print(f"Warning: Could not read box dimensions from PDB file: {e}")
        
        if x_len is None:
            x_len = y_len = z_len = 50.0
            print(f"WARNING: No box dimensions found! Using defaults: {x_len:.1f} x {y_len:.1f} x {z_len:.1f} Å")
        
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

                    # Centralize the system: keep lipid bilayer centered in Z
                    frame_pos = centralize_system(frame_pos, atom_types, pdb_residue_names, x_len, y_len, z_len)

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

                    # Centralize the system: keep lipid bilayer centered in Z
                    frame_pos = centralize_system(frame_pos, atom_types, pdb_residue_names, x_len, y_len, z_len)

                    # Write frame in PDB format
                    write_pdb_frame(f, frame_pos, frame, atom_types, residue_ids, x_len, y_len, z_len, pdb_atom_names, pdb_residue_names, pdb_residue_ids)

                    # Print progress every 100 frames
                    if frame % 100 == 0:
                        print(f"Processed frame {frame}/{n_frame_total}")

if __name__ == "__main__":
    main() 