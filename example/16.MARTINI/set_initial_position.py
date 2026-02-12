
#!/usr/bin/env python
import sys
import h5py
import numpy as np

def set_initial_position(input_file, output_file):
    # Open input file and get last frame's position and box dimensions
    with h5py.File(input_file, 'r') as f:
        # Try to read from output/pos (if simulation has already produced frames)
        if '/output/pos' in f and f['/output/pos'].shape[0] > 0:
            last_pos = f['/output/pos'][-1, 0, :, :]  # Last frame, first replica
            last_pos = last_pos[:, :, np.newaxis]  # Add replica dimension (1)
        else:
            # Fallback when output exists but is empty (e.g. minimize-only runs)
            last_pos = f['/input/pos'][:, :, 0]
            last_pos = last_pos[:, :, np.newaxis]

        # Get last frame's box dimensions if available
        last_box = None
        if '/output/box' in f:
            box_data = f['/output/box'][:]
            if box_data.size > 0:
                last_box = box_data[-1]
                # Handle different box data formats
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:
                    last_box = last_box[0]
        # If no output/box, try to get box from input potential attributes
        if last_box is None:
            if '/input/potential/martini_potential' in f:
                pot_grp = f['/input/potential/martini_potential']
                if all(k in pot_grp.attrs for k in ('x_len', 'y_len', 'z_len')):
                    last_box = np.array([
                        pot_grp.attrs['x_len'],
                        pot_grp.attrs['y_len'],
                        pot_grp.attrs['z_len']
                    ])

    # Open output file and write last frame's position to input/pos
    with h5py.File(output_file, 'r+') as f:
        if '/input/pos' in f:
            del f['/input/pos']
        f.create_dataset('/input/pos', data=last_pos)

        # Update box dimensions in martini_potential if available
        if last_box is not None and '/input/potential/martini_potential' in f:
            pot_grp = f['/input/potential/martini_potential']
            pot_grp.attrs['x_len'] = float(last_box[0])
            pot_grp.attrs['y_len'] = float(last_box[1])
            pot_grp.attrs['z_len'] = float(last_box[2])
            print(f"Updated box dimensions: x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 set_initial_position.py <input_file> <output_file>")
        sys.exit(1)
    set_initial_position(sys.argv[1], sys.argv[2])
