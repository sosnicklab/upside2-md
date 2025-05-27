#!/usr/bin/env python3

import sys, os
import numpy as np
import h5py

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_martini_vtf.py <input.up> <output.vtf>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

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
        
        # Get MARTINI data
        with h5py.File("inputs/martini.h5", 'r') as m:
            atoms = m['atoms'][:]  # Shape: (28, 2)
            centers = m['centers'][:]  # Shape: (28, 3)
            pairs = m['pairs'][:]  # Shape: (28, 2)
            print(f"MARTINI atoms shape: {atoms.shape}")
            print(f"MARTINI centers shape: {centers.shape}")
            print(f"MARTINI pairs shape: {pairs.shape}")
        
        # Create VTF file
        with open(output_file, 'w') as f:
            # Write structure
            f.write("structure {\n")
            
            # Write all particles as water (W)
            for i in range(n_particles):
                f.write(f"atom {i} name W resname WAT resid {i+1} reschain W\n")
            
            # Write bonds from pairs
            for i, j in pairs:
                f.write(f"bond {i}:{j}\n")
            
            f.write("}\n\n")
            
            # Find the last valid frame (no NaN values)
            last_valid_frame = None
            for frame in range(n_frame-1, -1, -1):
                pos = t['output/pos'][frame]
                if not np.isnan(pos).any():
                    last_valid_frame = frame
                    break
            
            if last_valid_frame is None:
                print("Error: No valid frames found in trajectory")
                sys.exit(1)
            
            print(f"Using frame {last_valid_frame} as reference for NaN positions")
            
            # Write frames
            for frame in range(n_frame):
                pos = t['output/pos'][frame]
                if frame == 0:
                    print(f"First frame sample: {pos[0,:5]}")
                    print(f"Position data type: {pos.dtype}")
                    print(f"Position data min: {np.min(pos)}")
                    print(f"Position data max: {np.max(pos)}")
                
                f.write(f"timestep {frame}\n")
                
                # Write all particle positions
                for i in range(n_particles):
                    # Get position and convert to float64 to avoid any precision issues
                    x, y, z = pos[0, i].astype(np.float64)
                    # Check for NaN values
                    if np.isnan(x) or np.isnan(y) or np.isnan(z):
                        # Use the last valid frame's position
                        valid_pos = t['output/pos'][last_valid_frame]
                        x, y, z = valid_pos[0, i].astype(np.float64)
                    # Format with fixed width for better compatibility
                    f.write(f"{x:12.6f} {y:12.6f} {z:12.6f}\n")
                
                # Add a newline after each timestep
                f.write("\n")

if __name__ == "__main__":
    main() 