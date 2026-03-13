import collections
import glob
import math
import os

import numpy as np


def _is_float(token):
    try:
        float(token)
    except ValueError:
        return False
    return True


def _strip_comment(line):
    return line.split("!", 1)[0].rstrip()


def _read_records(paths):
    records = []
    for path in paths:
        with open(path) as handle:
            carry = ""
            for raw_line in handle:
                line = _strip_comment(raw_line.rstrip("\n"))
                if not line:
                    continue
                if line.lstrip().startswith("*"):
                    continue
                line = line.strip()
                if carry:
                    line = carry + " " + line
                    carry = ""
                if line.endswith("-"):
                    carry = line[:-1].rstrip()
                    continue
                records.append(line)
            if carry:
                records.append(carry)
    return records


class ResidueTemplate:
    def __init__(self, name, total_charge, atoms, bonds, impropers):
        self.name = name
        self.total_charge = float(total_charge)
        self.atoms = atoms
        self.bonds = bonds
        self.impropers = impropers


class CharmmForceField:
    def __init__(self):
        self.residues = {}
        self.nonbonded = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = collections.defaultdict(list)
        self.impropers = []

    def bond_param(self, type1, type2):
        key = frozenset((type1, type2))
        if key not in self.bonds:
            raise KeyError("Missing bond parameter for %s-%s" % (type1, type2))
        return self.bonds[key]

    def angle_param(self, type1, type2, type3):
        key = (type1, type2, type3)
        rev_key = (type3, type2, type1)
        if key in self.angles:
            return self.angles[key]
        if rev_key in self.angles:
            return self.angles[rev_key]
        raise KeyError("Missing angle parameter for %s-%s-%s" % (type1, type2, type3))

    def dihedral_terms(self, type1, type2, type3, type4):
        direct = (type1, type2, type3, type4)
        reverse = (type4, type3, type2, type1)
        terms = []
        for pattern, pattern_terms in self.dihedrals.items():
            if _match_pattern(direct, pattern) or _match_pattern(reverse, pattern):
                terms.extend(pattern_terms)
        if not terms:
            raise KeyError("Missing dihedral parameter for %s-%s-%s-%s" % direct)
        return terms

    def improper_param(self, type1, type2, type3, type4):
        direct = (type1, type2, type3, type4)
        reverse = (type4, type3, type2, type1)
        for pattern, param in self.impropers:
            if _match_pattern(direct, pattern) or _match_pattern(reverse, pattern):
                return param
        raise KeyError("Missing improper parameter for %s-%s-%s-%s" % direct)


def _match_pattern(types, pattern):
    return all(p == "X" or p == t for p, t in zip(pattern, types))


def _resolve_residue_name(forcefield, residue_name):
    if residue_name in forcefield.residues:
        return residue_name
    if residue_name == "HIS":
        for candidate in ("HSD", "HSE", "HSP"):
            if candidate in forcefield.residues:
                return candidate
    raise KeyError("Missing residue template for %s" % residue_name)


def _atom_alias_candidates(residue_name, atom_name):
    aliases = {
        ("ILE", "CD1"): ("CD1", "CD"),
    }
    return aliases.get((residue_name, atom_name), (atom_name,))


def _finalize_residue(name, total_charge, atom_defs, bonds_raw, impropers_raw):
    atoms = collections.OrderedDict(atom_defs)
    atom_names = set(atoms)
    bonds = []
    for atom1, atom2 in bonds_raw:
        if atom1 in atom_names and atom2 in atom_names:
            bonds.append((atom1, atom2))
    impropers = []
    for atom1, atom2, atom3, atom4 in impropers_raw:
        if atom1 in atom_names and atom2 in atom_names and atom3 in atom_names and atom4 in atom_names:
            impropers.append((atom1, atom2, atom3, atom4))
    return ResidueTemplate(name, total_charge, atoms, bonds, impropers)


def _parse_topology(paths, forcefield):
    current = None
    total_charge = 0.0
    atom_defs = []
    bonds_raw = []
    impropers_raw = []

    def finish_current():
        nonlocal current, total_charge, atom_defs, bonds_raw, impropers_raw
        if current is not None:
            forcefield.residues[current] = _finalize_residue(
                current,
                total_charge,
                atom_defs,
                bonds_raw,
                impropers_raw,
            )
        current = None
        total_charge = 0.0
        atom_defs = []
        bonds_raw = []
        impropers_raw = []

    for line in _read_records(paths):
        tokens = line.split()
        key = tokens[0].upper()

        if key == "RESI":
            finish_current()
            current = tokens[1].upper()
            total_charge = float(tokens[2])
            continue
        if key in ("PRES", "END", "RETURN"):
            finish_current()
            continue
        if current is None:
            continue
        if key in (
            "GROUP",
            "MASS",
            "DECL",
            "DEFA",
            "AUTO",
            "PATCH",
            "PATCHING",
            "READ",
            "IC",
            "DELETE",
            "DONOR",
            "ACCEPTOR",
            "BILD",
            "CMAP",
        ):
            continue
        if key == "ATOM":
            atom_defs.append((tokens[1], {"type": tokens[2], "charge": float(tokens[3])}))
        elif key in ("BOND", "DOUBLE"):
            fields = tokens[1:]
            for idx in range(0, len(fields), 2):
                atom1 = fields[idx]
                atom2 = fields[idx + 1]
                if atom1.startswith(("+", "-")) or atom2.startswith(("+", "-")):
                    continue
                bonds_raw.append((atom1, atom2))
        elif key == "IMPR":
            fields = tokens[1:]
            for idx in range(0, len(fields), 4):
                atoms = tuple(fields[idx: idx + 4])
                if any(atom.startswith(("+", "-")) for atom in atoms):
                    continue
                impropers_raw.append(atoms)

    finish_current()


def _parse_parameters(paths, forcefield):
    section = None
    for line in _read_records(paths):
        tokens = line.split()
        key = tokens[0].upper()

        if key in ("BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "NONBONDED"):
            section = key
            continue
        if key in ("CMAP", "HBOND", "NBFIX", "END", "RETURN"):
            section = None
            continue
        if key in ("READ", "MASS"):
            continue
        if section is None:
            continue

        if section == "BONDS":
            if len(tokens) < 4:
                continue
            forcefield.bonds[frozenset((tokens[0], tokens[1]))] = (float(tokens[2]), float(tokens[3]))
        elif section == "ANGLES":
            if len(tokens) < 5:
                continue
            forcefield.angles[(tokens[0], tokens[1], tokens[2])] = (float(tokens[3]), float(tokens[4]))
        elif section == "DIHEDRALS":
            if len(tokens) < 7:
                continue
            pattern = (tokens[0], tokens[1], tokens[2], tokens[3])
            forcefield.dihedrals[pattern].append((float(tokens[4]), float(tokens[5]), math.radians(float(tokens[6]))))
        elif section == "IMPROPER":
            if len(tokens) < 6:
                continue
            pattern = (tokens[0], tokens[1], tokens[2], tokens[3])
            forcefield.impropers.append((pattern, (float(tokens[4]), math.radians(float(tokens[5])))))
        elif section == "NONBONDED":
            if key in ("NBXMOD", "CUTNB", "CTONNB", "CTOFNB", "EPS", "E14FAC", "WMIN", "ATOM", "CDIEL", "SHIFT", "VATOM", "VDISTANCE", "VSWITCH"):
                continue
            if len(tokens) < 4 or not _is_float(tokens[1]) or not _is_float(tokens[2]) or not _is_float(tokens[3]):
                continue
            forcefield.nonbonded[tokens[0]] = {
                "epsilon": abs(float(tokens[2])),
                "rmin_half": float(tokens[3]),
            }


def _unique_paths(paths):
    ordered = []
    seen = set()
    for path in paths:
        path = os.path.abspath(path)
        if path not in seen:
            seen.add(path)
            ordered.append(path)
    return ordered


def _split_env_paths(value):
    if not value:
        return []
    return [field for field in value.split(os.pathsep) if field]


def _search_existing_path(spec, search_roots):
    if os.path.isabs(spec):
        if os.path.exists(spec):
            return os.path.abspath(spec)
        raise IOError("CHARMM input path does not exist: %s" % spec)

    direct_hits = []
    for root in search_roots:
        for candidate in (
            os.path.join(root, spec),
            os.path.join(root, "toppar", spec),
        ):
            if os.path.exists(candidate):
                direct_hits.append(os.path.abspath(candidate))
    direct_hits = _unique_paths(direct_hits)
    if direct_hits:
        return direct_hits[0]

    recursive_hits = []
    if os.sep not in spec:
        for root in search_roots:
            recursive_hits.extend(
                glob.glob(os.path.join(root, "**", spec), recursive=True)
            )
    recursive_hits = _unique_paths(recursive_hits)
    if len(recursive_hits) == 1:
        return recursive_hits[0]
    if len(recursive_hits) > 1:
        raise IOError("Ambiguous CHARMM input path '%s': %s" % (spec, ", ".join(recursive_hits)))
    raise IOError("Unable to resolve CHARMM input path: %s" % spec)


def _default_topology_specs(toppar_dir, bundled_dir):
    if not toppar_dir:
        return ["protein_subset.rtf", "benzene.rtf"]
    specs = ["top_all36_prot.rtf", "benzene.rtf"]
    cgenff_top = None
    try:
        cgenff_top = _search_existing_path("top_all36_cgenff.rtf", [toppar_dir])
    except IOError:
        cgenff_top = None
    if cgenff_top is not None:
        specs.insert(1, "top_all36_cgenff.rtf")
    return specs


def _default_parameter_specs(toppar_dir, bundled_dir):
    if not toppar_dir:
        return ["protein_subset.prm", "benzene.prm"]
    specs = []
    for candidate in ("par_all36m_prot.prm", "par_all36_prot.prm"):
        try:
            _search_existing_path(candidate, [toppar_dir])
            specs.append(candidate)
            break
        except IOError:
            continue
    if not specs:
        raise IOError(
            "Unable to locate a protein parameter file under %s; set UPSIDE_CHARMM_PARAMETER_FILES"
            % toppar_dir
        )
    try:
        _search_existing_path("par_all36_cgenff.prm", [toppar_dir])
        specs.append("par_all36_cgenff.prm")
    except IOError:
        pass
    specs.append("benzene.prm")
    return specs


def resolve_forcefield_inputs(example_root):
    bundled_dir = os.path.join(example_root, "charmm")
    toppar_dir = os.environ.get("UPSIDE_CHARMM_TOPPAR_DIR")
    if not toppar_dir:
        default_toppar = os.path.join(example_root, "toppar")
        if os.path.isdir(default_toppar):
            toppar_dir = default_toppar
    topology_specs = _split_env_paths(os.environ.get("UPSIDE_CHARMM_TOPOLOGY_FILES"))
    parameter_specs = _split_env_paths(os.environ.get("UPSIDE_CHARMM_PARAMETER_FILES"))
    stream_specs = _split_env_paths(os.environ.get("UPSIDE_CHARMM_STREAM_FILES"))

    search_roots = [bundled_dir]
    if toppar_dir:
        search_roots.insert(0, os.path.abspath(toppar_dir))

    if not topology_specs:
        topology_specs = _default_topology_specs(toppar_dir, bundled_dir)
    if not parameter_specs:
        parameter_specs = _default_parameter_specs(toppar_dir, bundled_dir)

    return {
        "topology_paths": [_search_existing_path(spec, search_roots) for spec in topology_specs],
        "parameter_paths": [_search_existing_path(spec, search_roots) for spec in parameter_specs],
        "stream_paths": [_search_existing_path(spec, search_roots) for spec in stream_specs],
    }


def load_forcefield(topology_paths, parameter_paths, stream_paths=()):
    forcefield = CharmmForceField()
    topology_paths = _unique_paths(list(topology_paths) + list(stream_paths))
    parameter_paths = _unique_paths(list(parameter_paths) + list(stream_paths))
    _parse_topology(topology_paths, forcefield)
    _parse_parameters(parameter_paths, forcefield)
    return forcefield


def atom_params_for_residue(forcefield, residue_name, atom_names):
    residue_name = _resolve_residue_name(forcefield, residue_name)
    residue = forcefield.residues[residue_name]
    params = []
    for atom_name in atom_names:
        atom = None
        resolved_name = None
        for candidate in _atom_alias_candidates(residue_name, atom_name):
            atom = residue.atoms.get(candidate)
            if atom is not None:
                resolved_name = candidate
                break
        if atom is None:
            raise KeyError("Missing atom %s for residue %s" % (atom_name, residue_name))
        nb = forcefield.nonbonded[atom["type"]]
        params.append({
            "name": atom_name,
            "element": atom_name[0],
            "charge": atom["charge"],
            "rmin_half": nb["rmin_half"],
            "epsilon": nb["epsilon"],
            "type": atom["type"],
            "source_name": resolved_name,
        })
    return params


def first_atom_params_for_residue(forcefield, residue_name, atom_name_candidates):
    residue_name = _resolve_residue_name(forcefield, residue_name)
    residue = forcefield.residues[residue_name]
    for atom_name in atom_name_candidates:
        for candidate in _atom_alias_candidates(residue_name, atom_name):
            atom = residue.atoms.get(candidate)
            if atom is None:
                continue
            nb = forcefield.nonbonded[atom["type"]]
            return {
                "name": atom_name,
                "element": atom_name[0],
                "charge": atom["charge"],
                "rmin_half": nb["rmin_half"],
                "epsilon": nb["epsilon"],
                "type": atom["type"],
                "source_name": candidate,
            }
    raise KeyError(
        "Missing atom %s for residue %s" % ("/".join(atom_name_candidates), residue_name)
    )


def residue_sidechain_bonds(forcefield, residue_name, sidechain_names):
    residue_name = _resolve_residue_name(forcefield, residue_name)
    sidechain_names = list(sidechain_names)
    index = {}
    for idx, name in enumerate(sidechain_names):
        for candidate in _atom_alias_candidates(residue_name, name):
            index[candidate] = idx
    bonds = []
    for atom1, atom2 in forcefield.residues[residue_name].bonds:
        if atom1 in index and atom2 in index:
            bonds.append((index[atom1], index[atom2]))
    return np.asarray(bonds, dtype="i4").reshape((-1, 2))


def protein_backbone_params(forcefield, sequence, pocket_residue):
    protein_index = []
    protein_rmin_half = []
    protein_epsilon = []
    protein_charge = []

    for residue in pocket_residue:
        residue = int(residue)
        residue_name = "PRO" if sequence[residue] == "CPR" else sequence[residue]
        residue_params = {p["name"]: p for p in atom_params_for_residue(forcefield, residue_name, ("N", "CA", "C"))}
        for offset, atom_name in enumerate(("N", "CA", "C")):
            atom = residue_params[atom_name]
            protein_index.append(3 * residue + offset)
            protein_rmin_half.append(atom["rmin_half"])
            protein_epsilon.append(atom["epsilon"])
            protein_charge.append(atom["charge"])

    return (
        np.asarray(protein_index, dtype="i4"),
        np.asarray(protein_rmin_half, dtype="f4"),
        np.asarray(protein_epsilon, dtype="f4"),
        np.asarray(protein_charge, dtype="f4"),
    )


def protein_virtual_backbone_params(forcefield, infer_group, sequence, pocket_residue):
    pocket_residue = set(int(x) for x in pocket_residue)
    protein_index = []
    protein_rmin_half = []
    protein_epsilon = []
    protein_charge = []

    donor_residue = infer_group.donors.residue[:]
    acceptor_residue = infer_group.acceptors.residue[:]
    for idx, residue in enumerate(donor_residue):
        residue = int(residue)
        if residue not in pocket_residue:
            continue
        residue_name = "PRO" if sequence[residue] == "CPR" else sequence[residue]
        atom = first_atom_params_for_residue(forcefield, residue_name, ("H", "HN", "HN1", "HNT", "HT1"))
        protein_index.append(int(idx))
        protein_rmin_half.append(atom["rmin_half"])
        protein_epsilon.append(atom["epsilon"])
        protein_charge.append(atom["charge"])

    acceptor_offset = len(donor_residue)
    for idx, residue in enumerate(acceptor_residue):
        residue = int(residue)
        if residue not in pocket_residue:
            continue
        residue_name = "PRO" if sequence[residue] == "CPR" else sequence[residue]
        atom = first_atom_params_for_residue(forcefield, residue_name, ("O", "OT1"))
        protein_index.append(acceptor_offset + int(idx))
        protein_rmin_half.append(atom["rmin_half"])
        protein_epsilon.append(atom["epsilon"])
        protein_charge.append(atom["charge"])

    return (
        np.asarray(protein_index, dtype="i4"),
        np.asarray(protein_rmin_half, dtype="f4"),
        np.asarray(protein_epsilon, dtype="f4"),
        np.asarray(protein_charge, dtype="f4"),
    )


def _build_adjacency(atom_names, bond_pairs):
    adjacency = {name: set() for name in atom_names}
    for atom1, atom2 in bond_pairs:
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)
    return adjacency


def _enumerate_angles(atom_names, adjacency):
    seen = set()
    angles = []
    for center in atom_names:
        neighbors = sorted(adjacency[center])
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                angle = (neighbors[i], center, neighbors[j])
                reverse = (neighbors[j], center, neighbors[i])
                canonical = angle if angle < reverse else reverse
                if canonical not in seen:
                    seen.add(canonical)
                    angles.append(canonical)
    return angles


def _enumerate_dihedrals(atom_names, adjacency):
    seen = set()
    dihedrals = []
    for atom2 in atom_names:
        for atom3 in adjacency[atom2]:
            if atom2 >= atom3:
                continue
            left_neighbors = sorted(neighbor for neighbor in adjacency[atom2] if neighbor != atom3)
            right_neighbors = sorted(neighbor for neighbor in adjacency[atom3] if neighbor != atom2)
            for atom1 in left_neighbors:
                for atom4 in right_neighbors:
                    if len({atom1, atom2, atom3, atom4}) < 4:
                        continue
                    dihedral = (atom1, atom2, atom3, atom4)
                    reverse = dihedral[::-1]
                    canonical = dihedral if dihedral < reverse else reverse
                    if canonical not in seen:
                        seen.add(canonical)
                        dihedrals.append(dihedral)
    return dihedrals


def _nonbonded_pairs(atom_names, bonds, angles):
    excluded = set()
    for atom1, atom2 in bonds:
        excluded.add(tuple(sorted((atom1, atom2))))
    for atom1, _center, atom3 in angles:
        excluded.add(tuple(sorted((atom1, atom3))))

    pairs = []
    for i, atom1 in enumerate(atom_names):
        for atom2 in atom_names[i + 1:]:
            pair = tuple(sorted((atom1, atom2)))
            if pair not in excluded:
                pairs.append((atom1, atom2))
    return pairs


def build_ligand_asset(forcefield, residue_name, ligand_atoms):
    residue = forcefield.residues[residue_name]
    atom_names = [atom["atom_name"] for atom in ligand_atoms]
    elements = [atom["element"] for atom in ligand_atoms]
    coords = np.asarray([atom["coord"] for atom in ligand_atoms], dtype="f4")
    index = {name: idx for idx, name in enumerate(atom_names)}

    residue_params = atom_params_for_residue(forcefield, residue_name, atom_names)
    bonds = [(atom1, atom2) for atom1, atom2 in residue.bonds if atom1 in index and atom2 in index]
    adjacency = _build_adjacency(atom_names, bonds)
    angles = _enumerate_angles(atom_names, adjacency)
    dihedrals = _enumerate_dihedrals(atom_names, adjacency)
    nonbonded_pairs = _nonbonded_pairs(atom_names, bonds, angles)

    bond_pairs = []
    bond_length = []
    bond_spring = []
    for atom1, atom2 in bonds:
        idx1 = index[atom1]
        idx2 = index[atom2]
        atom_type1 = residue.atoms[atom1]["type"]
        atom_type2 = residue.atoms[atom2]["type"]
        spring, length = forcefield.bond_param(atom_type1, atom_type2)
        bond_pairs.append((idx1, idx2))
        bond_length.append(length)
        bond_spring.append(spring)

    angle_triples = []
    angle_theta0 = []
    angle_spring = []
    for atom1, atom2, atom3 in angles:
        idx1 = index[atom1]
        idx2 = index[atom2]
        idx3 = index[atom3]
        atom_type1 = residue.atoms[atom1]["type"]
        atom_type2 = residue.atoms[atom2]["type"]
        atom_type3 = residue.atoms[atom3]["type"]
        spring, theta0 = forcefield.angle_param(atom_type1, atom_type2, atom_type3)
        angle_triples.append((idx1, idx2, idx3))
        angle_theta0.append(math.radians(theta0))
        angle_spring.append(spring)

    torsion_quads = []
    torsion_periodicity = []
    torsion_phase = []
    torsion_amplitude = []
    for atom1, atom2, atom3, atom4 in dihedrals:
        idx1 = index[atom1]
        idx2 = index[atom2]
        idx3 = index[atom3]
        idx4 = index[atom4]
        atom_type1 = residue.atoms[atom1]["type"]
        atom_type2 = residue.atoms[atom2]["type"]
        atom_type3 = residue.atoms[atom3]["type"]
        atom_type4 = residue.atoms[atom4]["type"]
        for amplitude, periodicity, phase in forcefield.dihedral_terms(atom_type1, atom_type2, atom_type3, atom_type4):
            torsion_quads.append((idx1, idx2, idx3, idx4))
            torsion_periodicity.append(periodicity)
            torsion_phase.append(phase)
            torsion_amplitude.append(amplitude)

    improper_quads = []
    improper_theta0 = []
    improper_spring = []
    for atom1, atom2, atom3, atom4 in residue.impropers:
        if atom1 not in index or atom2 not in index or atom3 not in index or atom4 not in index:
            continue
        atom_type1 = residue.atoms[atom1]["type"]
        atom_type2 = residue.atoms[atom2]["type"]
        atom_type3 = residue.atoms[atom3]["type"]
        atom_type4 = residue.atoms[atom4]["type"]
        spring, theta0 = forcefield.improper_param(atom_type1, atom_type2, atom_type3, atom_type4)
        improper_quads.append((index[atom1], index[atom2], index[atom3], index[atom4]))
        improper_theta0.append(theta0)
        improper_spring.append(spring)

    pair_index = np.asarray([(index[a], index[b]) for a, b in nonbonded_pairs], dtype="i4").reshape((-1, 2))

    return {
        "resname": residue_name,
        "atom_names": atom_names,
        "elements": elements,
        "coords": coords,
        "charmm_rmin_half": np.asarray([p["rmin_half"] for p in residue_params], dtype="f4"),
        "charmm_epsilon": np.asarray([p["epsilon"] for p in residue_params], dtype="f4"),
        "charmm_charge": np.asarray([p["charge"] for p in residue_params], dtype="f4"),
        "bond_pairs": np.asarray(bond_pairs, dtype="i4").reshape((-1, 2)),
        "bond_length": np.asarray(bond_length, dtype="f4"),
        "bond_spring": np.asarray(bond_spring, dtype="f4"),
        "angle_triples": np.asarray(angle_triples, dtype="i4").reshape((-1, 3)),
        "angle_theta0": np.asarray(angle_theta0, dtype="f4"),
        "angle_spring": np.asarray(angle_spring, dtype="f4"),
        "torsion_quads": np.asarray(torsion_quads, dtype="i4").reshape((-1, 4)),
        "torsion_periodicity": np.asarray(torsion_periodicity, dtype="f4"),
        "torsion_phase": np.asarray(torsion_phase, dtype="f4"),
        "torsion_amplitude": np.asarray(torsion_amplitude, dtype="f4"),
        "improper_quads": np.asarray(improper_quads, dtype="i4").reshape((-1, 4)),
        "improper_theta0": np.asarray(improper_theta0, dtype="f4"),
        "improper_spring": np.asarray(improper_spring, dtype="f4"),
        "nonbonded_pairs": pair_index,
    }
