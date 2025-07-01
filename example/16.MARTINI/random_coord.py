import re

min_coords = [float('inf')] * 3
max_coords = [float('-inf')] * 3
atoms = []

# Read and extract atom lines
with open("ions.pdb", "r") as f:
    for line in f:
        if line.startswith("ATOM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((line, x, y, z))
            for i, val in enumerate([x, y, z]):
                min_coords[i] = min(min_coords[i], val)
                max_coords[i] = max(max_coords[i], val)

# Compute center and scale
center = [(minc + maxc) / 2 for minc, maxc in zip(min_coords, max_coords)]
range_max = max(maxc - minc for minc, maxc in zip(min_coords, max_coords))
scale = 46.0 / range_max  # to fit in [-23, 23]

# Apply transform and write output
with open("ions_out.pdb", "w") as f:
    for line, x, y, z in atoms:
        x_new = (x - center[0]) * scale
        y_new = (y - center[1]) * scale
        z_new = (z - center[2]) * scale
        newline = line[:30] + f"{x_new:8.3f}{y_new:8.3f}{z_new:8.3f}" + line[54:]
        f.write(newline)