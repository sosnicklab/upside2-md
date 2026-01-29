
#!/usr/bin/env python
import sys
import h5py
import numpy as np

def set_initial_position(input_file, output_file):
    # Open input file and get last frame's position
    with h5py.File(input_file, 'r') as f:
        last_pos = f['/output/pos'][-1, 0, :, :]  # Last frame, first replica, all atoms, all coordinates
        last_pos = last_pos[:, :, np.newaxis]  # Add replica dimension (1)

    # Open output file and write last frame's position to input/pos
    with h5py.File(output_file, 'r+') as f:
        if '/input/pos' in f:
            del f['/input/pos']
        f.create_dataset('/input/pos', data=last_pos)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 set_initial_position.py <input_file> <output_file>")
        sys.exit(1)
    set_initial_position(sys.argv[1], sys.argv[2])
