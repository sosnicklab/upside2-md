#!/usr/bin/env python3

import math
import os
import subprocess as sp
import sys

import numpy as np
import tables as tb

import charmm_import as ci
from pocket_library import build_sidechain_geometry


POCKET_CUTOFF_ANGSTROM = 8.0
SIM_ID = "smoke"
FF_NAME = "ff_2.1"
LIGAND_POCKET_RESTRAINT = 10.0
N_ANCHOR_PER_LIGAND_ATOM = 3
N_BIT_ROTAMER = 4

def repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def upside_home():
    return os.environ.get("UPSIDE_HOME", repo_root())


def py_dir():
    return os.path.join(upside_home(), "py")


def example_root():
    return os.path.abspath(os.path.dirname(__file__))


def input_dir():
    return os.path.join(example_root(), "inputs")


def output_dir():
    return os.path.join(example_root(), "outputs", SIM_ID)


def results_dir():
    return os.path.join(example_root(), "results")


def ensure_loader_compat():
    obj_dir = os.path.join(upside_home(), "obj")
    dylib = os.path.join(obj_dir, "libupside.dylib")
    so = os.path.join(obj_dir, "libupside.so")
    if os.path.exists(dylib) and not os.path.exists(so):
        os.symlink(dylib, so)


def load_upside_engine():
    ensure_loader_compat()
    py_path = py_dir()
    if py_path not in sys.path:
        sys.path.insert(0, py_path)
    import upside_engine as ue

    return ue


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def bstrings(items):
    return np.array([item.encode("ascii") for item in items])


def parse_pdb_atom(line):
    element = line[76:78].strip()
    return {
        "record": line[:6].strip(),
        "atom_name": line[12:16].strip(),
        "altloc": line[16],
        "resname": line[17:20].strip(),
        "chain": line[21].strip(),
        "resseq": line[22:26].strip(),
        "icode": line[26].strip(),
        "coord": np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype="f4",
        ),
        "element": element if element else line[12:16].strip()[0],
        "line": line.rstrip("\n"),
    }


def load_backbone_coords(pdb_path):
    coords = []
    with open(pdb_path) as handle:
        for raw_line in handle:
            if not raw_line.startswith("ATOM"):
                continue
            atom = parse_pdb_atom(raw_line)
            if atom["atom_name"] in ("N", "CA", "C"):
                coords.append(atom["coord"])
    return np.asarray(coords, dtype="f4")


def kabsch_row(mobile, target):
    center_mobile = mobile.mean(axis=0)
    center_target = target.mean(axis=0)
    mobile_centered = mobile - center_mobile
    target_centered = target - center_target
    cov = mobile_centered.T.dot(target_centered)
    u, _s, vt = np.linalg.svd(cov)
    rot = u.dot(vt)
    if np.linalg.det(rot) < 0.0:
        u[:, -1] *= -1.0
        rot = u.dot(vt)
    return rot.astype("f4"), center_mobile.astype("f4"), center_target.astype("f4")


def apply_rigid_transform(coords, rot, center_mobile, center_target):
    return (coords - center_mobile[None, :]).dot(rot) + center_target[None, :]


def transform_ligand_to_initial_frame(receptor_pdb_path, initial_structure_path, ligand_coords):
    receptor_backbone = load_backbone_coords(receptor_pdb_path)
    initial_backbone = np.load(initial_structure_path).astype("f4")
    rot, center_mobile, center_target = kabsch_row(receptor_backbone, initial_backbone)
    fitted_backbone = apply_rigid_transform(receptor_backbone, rot, center_mobile, center_target)
    fit_rms = np.sqrt(np.mean(np.sum((fitted_backbone - initial_backbone) ** 2, axis=1)))
    if fit_rms > 1e-3:
        raise RuntimeError("Protein frame alignment RMS is too large: %.6f" % fit_rms)
    return apply_rigid_transform(ligand_coords, rot, center_mobile, center_target).astype("f4")


def cleaned_protein_line(line):
    if line[16] == "A":
        line = line[:16] + " " + line[17:]
    return line


def load_complex(raw_pdb_path):
    protein_lines = []
    protein_atoms = []
    ligand_atoms = []
    residue_order = []
    residue_index = {}

    with open(raw_pdb_path) as handle:
        for raw_line in handle:
            if not (raw_line.startswith("ATOM") or raw_line.startswith("HETATM")):
                continue
            atom = parse_pdb_atom(raw_line)
            if atom["chain"] != "A":
                continue
            if atom["record"] == "ATOM":
                if atom["altloc"] not in (" ", "A"):
                    continue
                protein_lines.append(cleaned_protein_line(atom["line"]))
                protein_atoms.append(atom)
                key = (atom["chain"], atom["resseq"], atom["icode"])
                if key not in residue_index:
                    residue_index[key] = len(residue_order)
                    residue_order.append(key)
            elif atom["resname"] == "BNZ":
                ligand_atoms.append(atom)

    if not protein_lines:
        raise RuntimeError("No chain A protein atoms found in sample PDB")
    if len(ligand_atoms) != 6:
        raise RuntimeError("Expected 6 benzene atoms, found %d" % len(ligand_atoms))

    return protein_lines, protein_atoms, ligand_atoms, residue_index


def pocket_residues(protein_atoms, ligand_atoms, residue_index, cutoff):
    cutoff2 = cutoff * cutoff
    picked = set()
    for atom in protein_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        resid = residue_index[key]
        for lig in ligand_atoms:
            diff = atom["coord"] - lig["coord"]
            if np.dot(diff, diff) <= cutoff2:
                picked.add(resid)
                break
    return np.array(sorted(picked), dtype="i4")


def load_ligand_atoms(pdb_path):
    atoms = []
    with open(pdb_path) as handle:
        for raw_line in handle:
            if not (raw_line.startswith("ATOM") or raw_line.startswith("HETATM")):
                continue
            atoms.append(parse_pdb_atom(raw_line))
    if not atoms:
        raise RuntimeError("No ligand atoms found in %s" % pdb_path)
    return atoms


def write_sample_inputs(raw_pdb_path, receptor_path, ligand_path, pocket_path):
    protein_lines, protein_atoms, ligand_atoms, residue_index = load_complex(raw_pdb_path)
    with open(receptor_path, "w") as handle:
        for line in protein_lines:
            handle.write(line + "\n")
        handle.write("END\n")

    with open(ligand_path, "w") as handle:
        for atom in ligand_atoms:
            handle.write(atom["line"] + "\n")
        handle.write("END\n")

    picked = pocket_residues(protein_atoms, ligand_atoms, residue_index, POCKET_CUTOFF_ANGSTROM)
    np.savetxt(pocket_path, picked, fmt="%d")
    return picked


def run_python(args):
    print("Running:", " ".join(args))
    sp.check_call(args)


def load_pocket_residue_list(pocket_residue_path):
    pocket_residue = np.loadtxt(pocket_residue_path, dtype="i4")
    return np.atleast_1d(pocket_residue).astype("i4")


def ensure_infer_h_o(h5file, potential, sequence):
    if "infer_H_O" in potential:
        return

    donor_residue = np.array([i for i, aa in enumerate(sequence) if i > 0 and aa != "PRO"], dtype="i4")
    acceptor_residue = np.array([i for i in range(len(sequence) - 1)], dtype="i4")

    grp = create_group(h5file, potential, "infer_H_O", ["pos"])
    donors = h5file.create_group(grp, "donors")
    acceptors = h5file.create_group(grp, "acceptors")
    h5file.create_array(donors, "residue", donor_residue)
    h5file.create_array(acceptors, "residue", acceptor_residue)
    h5file.create_array(donors, "bond_length", np.full(len(donor_residue), 0.88, dtype="f4"))
    h5file.create_array(acceptors, "bond_length", np.full(len(acceptor_residue), 1.24, dtype="f4"))
    if len(donor_residue):
        donor_id = (np.array([-1, 0, 1], dtype="i4")[None, :] + 3 * donor_residue[:, None]).astype("i4")
    else:
        donor_id = np.zeros((0, 3), dtype="i4")
    if len(acceptor_residue):
        acceptor_id = (np.array([1, 2, 3], dtype="i4")[None, :] + 3 * acceptor_residue[:, None]).astype("i4")
    else:
        acceptor_id = np.zeros((0, 3), dtype="i4")
    h5file.create_array(donors, "id", donor_id)
    h5file.create_array(acceptors, "id", acceptor_id)


def canonical_restype(restype):
    return "PRO" if restype == "CPR" else restype


def load_rotamer_chi(sequence, affine_residue, layer_index, placement_library_path):
    with tb.open_file(placement_library_path) as data:
        restype_order = [x.decode("ascii") for x in data.root.restype_order[:]]
        restype_num = dict((aa, i) for i, aa in enumerate(restype_order))
        start_stop = data.root.rotamer_start_stop_bead[:]
        chi_lookup = {}
        for restype, chi1, chi2, state in data.root.restype_and_chi_and_state[:]:
            chi_lookup[(int(restype), int(state))] = (float(chi1), float(chi2))

    mapped_sequence = [("PRO" if aa == "CPR" else aa) for aa in sequence]
    restype_index = np.array([restype_num[mapped_sequence[int(residue)]] for residue in affine_residue], dtype="i4")
    chi = np.zeros((len(affine_residue), 4), dtype="f4")
    for idx, (restype, layer) in enumerate(zip(restype_index, layer_index)):
        start, stop, n_bead = [int(x) for x in start_stop[int(restype)]]
        if not (start <= int(layer) < stop):
            raise RuntimeError("Rotamer layer %d does not match residue type %d" % (int(layer), int(restype)))
        state = (int(layer) - start) // n_bead
        chi1, chi2 = chi_lookup.get((int(restype), int(state)), (0.0, 0.0))
        chi[idx, 0] = chi1
        chi[idx, 1] = chi2
    return restype_index, chi


def build_rotamer_atom_library(sequence, affine_residue, chi, forcefield):
    atom_start = np.zeros(len(affine_residue), dtype="i4")
    atom_count = np.zeros(len(affine_residue), dtype="i4")
    bond_start = np.zeros(len(affine_residue), dtype="i4")
    bond_count = np.zeros(len(affine_residue), dtype="i4")

    atom_name = []
    atom_element = []
    atom_local_pos = []
    atom_rmin_half = []
    atom_epsilon = []
    atom_charge = []
    bond_index = []

    for idx, (residue, chi_vec) in enumerate(zip(affine_residue, chi)):
        residue_name = canonical_restype(sequence[int(residue)])
        atom_start[idx] = len(atom_name)
        bond_start[idx] = len(bond_index)
        atom_block = build_sidechain_geometry(residue_name, chi_vec)
        imported_params = ci.atom_params_for_residue(forcefield, residue_name, atom_block["atom_name"])
        imported_bonds = ci.residue_sidechain_bonds(forcefield, residue_name, atom_block["atom_name"])

        atom_name.extend(atom_block["atom_name"])
        atom_element.extend(atom_block["atom_element"])
        atom_local_pos.extend(atom_block["atom_local_pos"].tolist())
        atom_rmin_half.extend([param["rmin_half"] for param in imported_params])
        atom_epsilon.extend([param["epsilon"] for param in imported_params])
        atom_charge.extend([param["charge"] for param in imported_params])
        bond_index.extend(imported_bonds.tolist())

        atom_count[idx] = len(atom_block["atom_name"])
        bond_count[idx] = imported_bonds.shape[0]

    return {
        "atom_start": atom_start,
        "atom_count": atom_count,
        "bond_start": bond_start,
        "bond_count": bond_count,
        "atom_name": bstrings(atom_name),
        "atom_element": bstrings(atom_element),
        "atom_local_pos": np.asarray(atom_local_pos, dtype="f4").reshape((-1, 3)),
        "atom_rmin_half": np.asarray(atom_rmin_half, dtype="f4"),
        "atom_epsilon": np.asarray(atom_epsilon, dtype="f4"),
        "atom_charge": np.asarray(atom_charge, dtype="f4"),
        "bond_index": np.asarray(bond_index, dtype="i4").reshape((-1, 2)),
    }


def create_group(h5file, parent, name, arguments):
    group = h5file.create_group(parent, name)
    group._v_attrs.arguments = bstrings(arguments)
    return group


def select_sidechain_anchor_pairs(config_path, ligand_coords, pocket_residues):
    ue = load_upside_engine()
    with tb.open_file(config_path) as t:
        pos = t.root.input.pos[:, :, 0].astype("f4")
        sc_node = t.root.input.potential.placement_fixed_point_vector_only
        affine_residue = sc_node.affine_residue[:]
        id_seq = sc_node.id_seq[:]

    engine = ue.Upside(config_path)
    engine.energy(pos)
    placement_pos = engine.get_output("placement_fixed_point_vector_only")[:, :3]
    placement_prob = engine.get_value_by_name((len(affine_residue),), "rotamer", "placement_marginal")
    selector = (1 << N_BIT_ROTAMER) - 1
    rotamer_id = id_seq & selector

    candidate_index = []
    for residue in pocket_residues:
        residue_mask = affine_residue == residue
        if not np.any(residue_mask):
            continue
        best_rot = None
        best_prob = -1.0
        for rot in sorted(set(int(x) for x in rotamer_id[residue_mask])):
            rot_mask = residue_mask & (rotamer_id == rot)
            rot_prob = float(np.mean(placement_prob[rot_mask]))
            if rot_prob > best_prob:
                best_prob = rot_prob
                best_rot = rot
        if best_rot is None:
            continue
        candidate_index.extend(np.flatnonzero(residue_mask & (rotamer_id == best_rot)).tolist())

    if not candidate_index:
        raise RuntimeError("Failed to identify pocket sidechain anchor placements")

    candidate_index = np.asarray(candidate_index, dtype="i4")
    candidate_xyz = placement_pos[candidate_index]
    pair_ids = []
    pair_dists = []
    for ligand_atom, ligand_xyz in enumerate(ligand_coords):
        deltas = candidate_xyz - ligand_xyz[None, :]
        dist2 = np.sum(deltas * deltas, axis=1)
        ranked = np.argsort(dist2)[:min(N_ANCHOR_PER_LIGAND_ATOM, len(candidate_index))]
        for best in ranked:
            pair_ids.append([ligand_atom, int(candidate_index[int(best)])])
            pair_dists.append(np.sqrt(dist2[int(best)]))
    pair_ids = np.asarray(pair_ids, dtype="i4")
    pair_dists = np.asarray(pair_dists, dtype="f4")
    return pair_ids, pair_dists


def add_ligand_pocket_restraint(h5file, potential, pair_ids, pair_dists):
    anchor_dist = create_group(
        h5file,
        potential,
        "Distance3D_ligand_anchor",
        ["slice_ligand", "placement_fixed_point_vector_only"],
    )
    h5file.create_array(anchor_dist, "id", pair_ids)

    anchor_spring = create_group(h5file, potential, "Spring_ligand_anchor", ["Distance3D_ligand_anchor"])
    anchor_spring._v_attrs.dim1 = 0
    anchor_spring._v_attrs.pbc = 0
    anchor_spring._v_attrs.box_len = 0.0
    h5file.create_array(anchor_spring, "id", np.arange(pair_ids.shape[0], dtype="i4"))
    h5file.create_array(anchor_spring, "equil_dist", pair_dists)
    h5file.create_array(
        anchor_spring,
        "spring_const",
        np.full(pair_ids.shape[0], LIGAND_POCKET_RESTRAINT, dtype="f4"),
    )


def patch_config(config_path, ligand_asset, pocket_residue_path, forcefield):
    with tb.open_file(config_path, "a") as t:
        protein_pos = t.root.input.pos[:, :, 0]
        n_protein_atom = protein_pos.shape[0]
        n_ligand_atom = ligand_asset["coords"].shape[0]
        pocket_residue = load_pocket_residue_list(pocket_residue_path)
        new_pos = np.concatenate((protein_pos, ligand_asset["coords"]), axis=0)[:, :, None]

        t.remove_node("/input/pos")
        t.create_array(t.root.input, "pos", new_pos)

        if "ligand" in t.root.input:
            t.remove_node("/input/ligand", recursive=True)
        ligand_grp = t.create_group(t.root.input, "ligand")
        ligand_grp._v_attrs.resname = ligand_asset["resname"]
        ligand_grp._v_attrs.atom_start = n_protein_atom
        ligand_grp._v_attrs.pocket_cutoff_angstrom = POCKET_CUTOFF_ANGSTROM
        t.create_array(ligand_grp, "initial_pos", ligand_asset["coords"])
        t.create_array(ligand_grp, "atom_name", bstrings(ligand_asset["atom_names"]))
        t.create_array(ligand_grp, "element", bstrings(ligand_asset["elements"]))
        t.create_array(ligand_grp, "bond_index", ligand_asset["bond_pairs"])
        t.create_array(ligand_grp, "pocket_residue", pocket_residue)

        potential = t.root.input.potential

        slice_grp = create_group(t, potential, "slice_ligand", ["pos"])
        t.create_array(slice_grp, "id", np.arange(n_protein_atom, n_protein_atom + n_ligand_atom, dtype="i4"))

        bond_grp = create_group(t, potential, "ligand_harmonic_bond", ["slice_ligand"])
        t.create_array(bond_grp, "atom1", ligand_asset["bond_pairs"][:, 0])
        t.create_array(bond_grp, "atom2", ligand_asset["bond_pairs"][:, 1])
        t.create_array(bond_grp, "length", ligand_asset["bond_length"])
        t.create_array(bond_grp, "spring", ligand_asset["bond_spring"])

        angle_grp = create_group(t, potential, "ligand_angle", ["slice_ligand"])
        t.create_array(angle_grp, "atom1", ligand_asset["angle_triples"][:, 0])
        t.create_array(angle_grp, "atom2", ligand_asset["angle_triples"][:, 1])
        t.create_array(angle_grp, "atom3", ligand_asset["angle_triples"][:, 2])
        t.create_array(angle_grp, "theta0", ligand_asset["angle_theta0"])
        t.create_array(angle_grp, "spring", ligand_asset["angle_spring"])

        torsion_grp = create_group(t, potential, "ligand_torsion", ["slice_ligand"])
        t.create_array(torsion_grp, "atom1", ligand_asset["torsion_quads"][:, 0])
        t.create_array(torsion_grp, "atom2", ligand_asset["torsion_quads"][:, 1])
        t.create_array(torsion_grp, "atom3", ligand_asset["torsion_quads"][:, 2])
        t.create_array(torsion_grp, "atom4", ligand_asset["torsion_quads"][:, 3])
        t.create_array(torsion_grp, "periodicity", ligand_asset["torsion_periodicity"])
        t.create_array(torsion_grp, "phase", ligand_asset["torsion_phase"])
        t.create_array(torsion_grp, "amplitude", ligand_asset["torsion_amplitude"])

        if len(ligand_asset["improper_quads"]):
            improper_grp = create_group(t, potential, "ligand_improper", ["slice_ligand"])
            t.create_array(improper_grp, "atom1", ligand_asset["improper_quads"][:, 0])
            t.create_array(improper_grp, "atom2", ligand_asset["improper_quads"][:, 1])
            t.create_array(improper_grp, "atom3", ligand_asset["improper_quads"][:, 2])
            t.create_array(improper_grp, "atom4", ligand_asset["improper_quads"][:, 3])
            t.create_array(improper_grp, "theta0", ligand_asset["improper_theta0"])
            t.create_array(improper_grp, "spring", ligand_asset["improper_spring"])

        pair_index = ligand_asset["nonbonded_pairs"]
        pair_rmin = []
        pair_epsilon = []
        pair_charge = []
        for left, right in pair_index:
            pair_rmin.append(ligand_asset["charmm_rmin_half"][left] + ligand_asset["charmm_rmin_half"][right])
            pair_epsilon.append(math.sqrt(ligand_asset["charmm_epsilon"][left] * ligand_asset["charmm_epsilon"][right]))
            pair_charge.append(ligand_asset["charmm_charge"][left] * ligand_asset["charmm_charge"][right])

        nb_grp = create_group(t, potential, "charmm_ligand_nonbonded", ["slice_ligand"])
        nb_grp._v_attrs.cutoff = 12.0
        t.create_array(nb_grp, "atom1", pair_index[:, 0])
        t.create_array(nb_grp, "atom2", pair_index[:, 1])
        t.create_array(nb_grp, "rmin", np.asarray(pair_rmin, dtype="f4"))
        t.create_array(nb_grp, "epsilon", np.asarray(pair_epsilon, dtype="f4"))
        t.create_array(nb_grp, "charge", np.asarray(pair_charge, dtype="f4"))

        sequence = [x.decode("ascii") for x in t.root.input.sequence[:]]
        ensure_infer_h_o(t, potential, sequence)

        protein_index, protein_rmin_half, protein_epsilon, protein_charge = ci.protein_backbone_params(
            forcefield, sequence, pocket_residue
        )
        bb_grp = create_group(t, potential, "ligand_backbone_nonbonded", ["slice_ligand", "pos"])
        bb_grp._v_attrs.cutoff = 12.0
        t.create_array(bb_grp, "protein_index", protein_index)
        t.create_array(bb_grp, "ligand_rmin_half", ligand_asset["charmm_rmin_half"])
        t.create_array(bb_grp, "ligand_epsilon", ligand_asset["charmm_epsilon"])
        t.create_array(bb_grp, "ligand_charge", ligand_asset["charmm_charge"])
        t.create_array(bb_grp, "protein_rmin_half", protein_rmin_half)
        t.create_array(bb_grp, "protein_epsilon", protein_epsilon)
        t.create_array(bb_grp, "protein_charge", protein_charge)

        infer_group = t.root.input.potential.infer_H_O
        virtual_index, virtual_rmin_half, virtual_epsilon, virtual_charge = ci.protein_virtual_backbone_params(
            forcefield, infer_group, sequence, pocket_residue
        )
        virtual_grp = create_group(t, potential, "ligand_virtual_backbone_nonbonded", ["slice_ligand", "infer_H_O"])
        virtual_grp._v_attrs.cutoff = 12.0
        t.create_array(virtual_grp, "protein_index", virtual_index)
        t.create_array(virtual_grp, "ligand_rmin_half", ligand_asset["charmm_rmin_half"])
        t.create_array(virtual_grp, "ligand_epsilon", ligand_asset["charmm_epsilon"])
        t.create_array(virtual_grp, "ligand_charge", ligand_asset["charmm_charge"])
        t.create_array(virtual_grp, "protein_rmin_half", virtual_rmin_half)
        t.create_array(virtual_grp, "protein_epsilon", virtual_epsilon)
        t.create_array(virtual_grp, "protein_charge", virtual_charge)

        sc_node = t.root.input.potential.placement_fixed_point_vector_only
        restype_index, chi = load_rotamer_chi(
            sequence,
            sc_node.affine_residue[:],
            sc_node.layer_index[:],
            os.path.join(upside_home(), "parameters", FF_NAME, "sidechain.h5"),
        )
        rotamer_grp = create_group(
            t,
            potential,
            "ligand_rotamer_charmm_1body",
            ["affine_alignment", "slice_ligand"],
        )
        rotamer_grp._v_attrs.cutoff = 12.0
        t.create_array(rotamer_grp, "affine_residue", sc_node.affine_residue[:])
        t.create_array(rotamer_grp, "restype_index", restype_index)
        t.create_array(rotamer_grp, "chi", chi)
        t.create_array(rotamer_grp, "ligand_rmin_half", ligand_asset["charmm_rmin_half"])
        t.create_array(rotamer_grp, "ligand_epsilon", ligand_asset["charmm_epsilon"])
        t.create_array(rotamer_grp, "ligand_charge", ligand_asset["charmm_charge"])
        atom_library = build_rotamer_atom_library(sequence, sc_node.affine_residue[:], chi, forcefield)
        t.create_array(rotamer_grp, "atom_start", atom_library["atom_start"])
        t.create_array(rotamer_grp, "atom_count", atom_library["atom_count"])
        t.create_array(rotamer_grp, "bond_start", atom_library["bond_start"])
        t.create_array(rotamer_grp, "bond_count", atom_library["bond_count"])
        t.create_array(rotamer_grp, "atom_name", atom_library["atom_name"])
        t.create_array(rotamer_grp, "atom_element", atom_library["atom_element"])
        t.create_array(rotamer_grp, "atom_local_pos", atom_library["atom_local_pos"])
        t.create_array(rotamer_grp, "atom_rmin_half", atom_library["atom_rmin_half"])
        t.create_array(rotamer_grp, "atom_epsilon", atom_library["atom_epsilon"])
        t.create_array(rotamer_grp, "atom_charge", atom_library["atom_charge"])
        t.create_array(rotamer_grp, "bond_index", atom_library["bond_index"])

        rotamer_args = [arg.decode("ascii") if isinstance(arg, bytes) else str(arg) for arg in t.root.input.potential.rotamer._v_attrs.arguments]
        if "ligand_rotamer_charmm_1body" not in rotamer_args:
            rotamer_args.append("ligand_rotamer_charmm_1body")
            t.root.input.potential.rotamer._v_attrs.arguments = bstrings(rotamer_args)


def main():
    ensure_dir(input_dir())
    ensure_dir(output_dir())
    ensure_dir(results_dir())

    raw_pdb = os.path.join(example_root(), "pdb", "4W52.pdb")
    receptor_pdb = os.path.join(input_dir(), "4w52_receptor.pdb")
    ligand_pdb = os.path.join(input_dir(), "4w52_benzene.pdb")
    pocket_path = os.path.join(input_dir(), "pocket_residues.txt")
    picked = write_sample_inputs(raw_pdb, receptor_pdb, ligand_pdb, pocket_path)

    print("Pocket residues:", " ".join(str(x) for x in picked))

    base = os.path.join(input_dir(), "4w52_receptor")
    run_python([
        sys.executable,
        os.path.join(py_dir(), "PDB_to_initial_structure.py"),
        receptor_pdb,
        base,
        "--record-chain-breaks",
    ])

    config_path = os.path.join(input_dir(), "4w52_ligand.up")
    param_base = os.path.join(upside_home(), "parameters")
    ff_dir = os.path.join(param_base, FF_NAME)
    common_dir = os.path.join(param_base, "common")

    run_python([
        sys.executable,
        os.path.join(py_dir(), "upside_config.py"),
        "--output=%s" % config_path,
        "--fasta=%s.fasta" % base,
        "--initial-structure=%s.initial.npy" % base,
        "--rama-library=%s" % os.path.join(common_dir, "rama.dat"),
        "--rama-sheet-mixing-energy=%s" % os.path.join(ff_dir, "sheet"),
        "--reference-state-rama=%s" % os.path.join(common_dir, "rama_reference.pkl"),
        "--hbond-energy=%s" % os.path.join(ff_dir, "hbond.h5"),
        "--dynamic-rotamer-1body",
        "--rotamer-placement=%s" % os.path.join(ff_dir, "sidechain.h5"),
        "--rotamer-interaction=%s" % os.path.join(ff_dir, "sidechain.h5"),
        "--environment-potential=%s" % os.path.join(ff_dir, "environment.h5"),
        "--bb-environment-potential=%s" % os.path.join(ff_dir, "bb_env.dat"),
    ])

    ff_inputs = ci.resolve_forcefield_inputs(example_root())
    print("Using CHARMM topology files:")
    for path in ff_inputs["topology_paths"]:
        print("  ", path)
    print("Using CHARMM parameter files:")
    for path in ff_inputs["parameter_paths"]:
        print("  ", path)
    if ff_inputs["stream_paths"]:
        print("Using CHARMM stream files:")
        for path in ff_inputs["stream_paths"]:
            print("  ", path)

    forcefield = ci.load_forcefield(
        ff_inputs["topology_paths"],
        ff_inputs["parameter_paths"],
        ff_inputs["stream_paths"],
    )
    ligand_asset = ci.build_ligand_asset(forcefield, "BNZ", load_ligand_atoms(ligand_pdb))
    ligand_asset["coords"] = transform_ligand_to_initial_frame(
        receptor_pdb,
        "%s.initial.npy" % base,
        ligand_asset["coords"],
    )
    patch_config(config_path, ligand_asset, pocket_path, forcefield)
    pair_ids, pair_dists = select_sidechain_anchor_pairs(config_path, ligand_asset["coords"], picked)
    with tb.open_file(config_path, "a") as t:
        add_ligand_pocket_restraint(t, t.root.input.potential, pair_ids, pair_dists)
    print("Prepared ligand-aware config:", config_path)


if __name__ == "__main__":
    main()
