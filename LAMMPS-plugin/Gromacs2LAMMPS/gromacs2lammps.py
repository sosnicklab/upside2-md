import sys
import re

def extract_section(lines, section):
    """Extract [ section ] from ITP file."""
    in_section = False
    section_data = []
    for line in lines:
        line = line.strip()
        if line.startswith(f"[ {section} ]"):
            in_section = True
            continue
        if in_section:
            if line.startswith("[") or line == "":
                break
            if line.startswith(";"):
                continue
            section_data.append(line.split())
    return section_data

def parse_pdb(pdb_file):
    atoms = []
    box = [0.0, 0.0, 0.0]
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_id = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "id": atom_id,
                    "name": atom_name,
                    "resname": res_name,
                    "resid": res_id,
                    "x": x, "y": y, "z": z
                })
            elif line.startswith("CRYST1"):
                box = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
    return atoms, box

def write_lammps_data(filename, pdb_atoms, itp_atoms, bonds, angles, dihedrals, box):
    type_map = {}
    type_counter = 1
    mass_map = {}
    charge_map = {}
    atom_lines = []

    for i, pdb_atom in enumerate(pdb_atoms):
        itp_atom = itp_atoms[i]
        atype = itp_atom["type"]
        if atype not in type_map:
            type_map[atype] = type_counter
            mass_map[type_counter] = itp_atom["mass"]
            charge_map[type_counter] = itp_atom["charge"]
            type_counter += 1
        atom_lines.append((pdb_atom["id"], pdb_atom["resid"], type_map[atype],
                           itp_atom["charge"], pdb_atom["x"], pdb_atom["y"], pdb_atom["z"]))

    with open(filename, 'w') as f:
        f.write("LAMMPS data file generated from PDB and ITP\n\n")
        f.write(f"{len(pdb_atoms)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write(f"{len(dihedrals)} dihedrals\n\n")
        f.write(f"{len(type_map)} atom types\n")
        f.write("1 bond types\n1 angle types\n1 dihedral types\n\n")  # simple case

        f.write(f"0.0 {box[0]:.3f} xlo xhi\n")
        f.write(f"0.0 {box[1]:.3f} ylo yhi\n")
        f.write(f"0.0 {box[2]:.3f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for tid, mass in mass_map.items():
            f.write(f"{tid} {mass:.4f}\n")

        f.write("\nAtoms\n\n")
        for a in atom_lines:
            f.write(f"{a[0]} {a[1]} {a[2]} {a[3]:.4f} {a[4]:.4f} {a[5]:.4f} {a[6]:.4f}\n")

        f.write("\nBonds\n\n")
        for i, (a1, a2) in enumerate(bonds, 1):
            f.write(f"{i} 1 {a1} {a2}\n")

        f.write("\nAngles\n\n")
        for i, (a1, a2, a3) in enumerate(angles, 1):
            f.write(f"{i} 1 {a1} {a2} {a3}\n")

        f.write("\nDihedrals\n\n")
        for i, (a1, a2, a3, a4) in enumerate(dihedrals, 1):
            f.write(f"{i} 1 {a1} {a2} {a3} {a4}\n")

    print(f"✅ LAMMPS data written to: {filename}")

def parse_itp_atoms(lines):
    raw_section = extract_section(lines, "atoms")
    atoms = []
    for tokens in raw_section:
        # Remove inline comments (strip anything after ';')
        if ';' in tokens:
            tokens = tokens[:tokens.index(';')]
        if len(tokens) < 7:
            continue
        try:
            atoms.append({
                "id": int(tokens[0]),
                "type": tokens[1],
                "resid": int(tokens[2]),
                "resname": tokens[3],
                "atomname": tokens[4],
                "charge": float(tokens[6]),
                "mass": float(tokens[7]) if len(tokens) >= 8 else 0.0
            })
        except ValueError:
            print(f"⚠ Skipping atom line due to conversion error: {tokens}")
            continue
    return atoms

def parse_itp_pairs(lines, section):
    data = extract_section(lines, section)
    if section == "bonds":
        return [(int(a[0]), int(a[1])) for a in data]
    elif section == "angles":
        return [(int(a[0]), int(a[1]), int(a[2])) for a in data]
    elif section == "dihedrals":
        return [(int(a[0]), int(a[1]), int(a[2]), int(a[3])) for a in data]
    return []

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python pdb_itp_to_lammps_full.py <input.pdb> <output.data> <top.itp> [extra.itp ...]")
        sys.exit(1)

    pdb_file = sys.argv[1]
    out_file = sys.argv[2]
    itp_files = sys.argv[3:]

    # Use the first ITP as the main one for atoms, rest are support
    with open(itp_files[0], 'r') as f:
        itp_lines = f.readlines()
        atoms = parse_itp_atoms(itp_lines)
        bonds = parse_itp_pairs(itp_lines, "bonds")
        angles = parse_itp_pairs(itp_lines, "angles")
        dihedrals = parse_itp_pairs(itp_lines, "dihedrals")

    pdb_atoms, box = parse_pdb(pdb_file)

    if len(pdb_atoms) != len(atoms):
        print(f"⚠ Warning: Atom count mismatch (PDB={len(pdb_atoms)} vs ITP={len(atoms)})")

    if len(pdb_atoms) > len(itp_atoms):
        print(f"⚠ Truncating PDB atoms from {len(pdb_atoms)} to match ITP atoms ({len(itp_atoms)})")
        pdb_atoms = pdb_atoms[:len(itp_atoms)]
    elif len(pdb_atoms) < len(itp_atoms):
        raise ValueError("ITP defines more atoms than available in the PDB file.")

    write_lammps_data(out_file, pdb_atoms, atoms, bonds, angles, dihedrals, box)