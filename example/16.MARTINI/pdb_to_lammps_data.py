import io

def convert_pdb_to_lammps_data(pdb_content, box_dimensions_angstroms):
    """
    Converts PDB content (as a string) to a LAMMPS data file format.

    Args:
        pdb_content (str): The content of the PDB file as a string.
        box_dimensions_angstroms (dict): A dictionary specifying the box dimensions.
                                         Example: {'x': [-23, 23], 'y': [-23, 23], 'z': [-23, 23]}

    Returns:
        str: A string containing the LAMMPS data file content.
    """
    atoms_data = []
    # In this specific case, all atoms are 'W' and map to LAMMPS atom type 1.
    atom_type_map = {'W': 1}
    num_atom_types = 1 # Assuming only 'W' type is present in this PDB snippet

    # Parse the PDB content
    for line in pdb_content.splitlines():
        if line.startswith("ATOM"):
            try:
                # Atom ID from PDB (column 7-11, 1-indexed in PDB spec)
                atom_id = int(line[6:11].strip())
                # Atom name (column 13-16)
                atom_name = line[12:16].strip()
                # Coordinates (X, Y, Z from column 31-38, 39-46, 47-54)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                lammps_atom_type = atom_type_map.get(atom_name)
                if lammps_atom_type is None:
                    print(f"Warning: Atom type '{atom_name}' not found in map. Defaulting to 1.")
                    lammps_atom_type = 1 # Default to 1 if not in map

                atoms_data.append({
                    'id': atom_id,
                    'type': lammps_atom_type,
                    'x': x,
                    'y': y,
                    'z': z
                })
            except ValueError as e:
                print(f"Error parsing line: {line.strip()} - {e}")
                continue

    num_atoms = len(atoms_data)

    # Prepare LAMMPS data file content
    lammps_data_output = io.StringIO()

    lammps_data_output.write("LAMMPS data file generated from PDB\n\n")
    lammps_data_output.write(f"{num_atoms} atoms\n")
    lammps_data_output.write(f"{num_atom_types} atom types\n\n")

    # Box dimensions
    xlo = box_dimensions_angstroms['x'][0]
    xhi = box_dimensions_angstroms['x'][1]
    ylo = box_dimensions_angstroms['y'][0]
    yhi = box_dimensions_angstroms['y'][1]
    zlo = box_dimensions_angstroms['z'][0]
    zhi = box_dimensions_angstroms['z'][1]

    lammps_data_output.write(f"{xlo:.3f} {xhi:.3f} xlo xhi\n")
    lammps_data_output.write(f"{ylo:.3f} {yhi:.3f} ylo yhi\n")
    lammps_data_output.write(f"{zlo:.3f} {zhi:.3f} zlo zhi\n\n")

    # Masses (assuming mass 72.0 for atom type 1 as per your LAMMPS script)
    lammps_data_output.write("Masses\n\n")
    lammps_data_output.write("1 72.0\n\n") # Assuming only one atom type 'W' with mass 72.0

    # Atoms section
    lammps_data_output.write("Atoms\n\n")
    for atom in atoms_data:
        lammps_data_output.write(f"{atom['id']} {atom['type']} {atom['x']:.3f} {atom['y']:.3f} {atom['z']:.3f}\n")

    return lammps_data_output.getvalue()

# Example usage (you would replace this with your actual PDB content and desired box dimensions)
pdb_text = open("input.pdb").read()
box_dims = {'x': [-23, 23], 'y': [-23, 23], 'z': [-23, 23]}
lammps_data = convert_pdb_to_lammps_data(pdb_text, box_dims)
print(lammps_data)
