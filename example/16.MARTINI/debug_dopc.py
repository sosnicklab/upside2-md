#!/usr/bin/env python3

from read_martini3_topology import read_martini3_atoms, read_martini3_bonds, read_martini3_angles

# Test DOPC parsing
itp_file = "ff3.00/martini_v3.0.0_phospholipids_v1.itp"

print("Testing DOPC parsing with fixed functions...")

# Test atoms
bead_types, charges = read_martini3_atoms(itp_file, "DOPC")
print(f"Atoms: {len(bead_types)} found")
print(f"Bead types: {bead_types}")
print(f"Charges: {charges}")

# Test bonds
bonds, lengths, forces = read_martini3_bonds(itp_file, "DOPC")
print(f"Bonds: {len(bonds)} found")
if bonds:
    print(f"First few bonds: {bonds[:5]}")
    print(f"Bond lengths (nm): {lengths[:5]}")

# Test angles
angles, values, forces = read_martini3_angles(itp_file, "DOPC")
print(f"Angles: {len(angles)} found")
if angles:
    print(f"First few angles: {angles[:5]}")
    print(f"Angle values (deg): {values[:5]}")
