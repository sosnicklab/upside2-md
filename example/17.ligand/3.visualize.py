#!/usr/bin/env python3

import os
import sys

import numpy as np
import tables as tb


N_BIT_ROTAMER = 4
SIM_ID = "smoke"


def repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def upside_home():
    return os.environ.get("UPSIDE_HOME", repo_root())


def ensure_loader_compat():
    obj_dir = os.path.join(upside_home(), "obj")
    dylib = os.path.join(obj_dir, "libupside.dylib")
    so = os.path.join(obj_dir, "libupside.so")
    if os.path.exists(dylib) and not os.path.exists(so):
        os.symlink(dylib, so)


ensure_loader_compat()


sys.path.insert(0, os.path.join(upside_home(), "py"))
import upside_engine as ue  # noqa: E402


def example_root():
    return os.path.abspath(os.path.dirname(__file__))


def run_path():
    return os.path.join(example_root(), "outputs", SIM_ID, "4w52_ligand.run.up")


def results_dir():
    return os.path.join(example_root(), "results")


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def latest_output_positions(handle):
    chunks = []
    idx = 0
    while "output_previous_%d" % idx in handle.root:
        group = handle.get_node("/output_previous_%d" % idx)
        chunks.append(group.pos[:, 0] if idx == 0 else group.pos[1:, 0])
        idx += 1
    if "output" in handle.root:
        group = handle.get_node("/output")
        chunks.append(group.pos[:, 0] if not chunks else group.pos[1:, 0])
    if not chunks:
        raise RuntimeError("No trajectory positions found")
    return np.concatenate(chunks, axis=0)


def backbone_bonds(n_res):
    bonds = []
    for resid in range(n_res):
        base = resid * 3
        bonds.append((base + 0, base + 1))
        bonds.append((base + 1, base + 2))
        if resid + 1 < n_res:
            bonds.append((base + 2, base + 3))
    return bonds


def quat_to_rot(quat):
    a, b, c, d = quat
    return np.array([
        [a * a + b * b - c * c - d * d, 2.0 * b * c - 2.0 * a * d, 2.0 * b * d + 2.0 * a * c],
        [2.0 * b * c + 2.0 * a * d, a * a - b * b + c * c - d * d, 2.0 * c * d - 2.0 * a * b],
        [2.0 * b * d - 2.0 * a * c, 2.0 * c * d + 2.0 * a * b, a * a - b * b - c * c + d * d],
    ], dtype="f4")


def apply_affine(affine, local_pos):
    rot = quat_to_rot(affine[3:])
    return local_pos.dot(rot.T) + affine[:3][None, :]


def write_vtf(path, atom_records, bonds, frames):
    with open(path, "w") as handle:
        for atom_id, record in enumerate(atom_records):
            handle.write(
                "atom %d name %s resid %d resname %s segid %s chain %s\n" % (
                    atom_id,
                    record["name"],
                    record["resid"],
                    record["resname"],
                    record["segid"],
                    record["chain"],
                )
            )
        for atom_a, atom_b in bonds:
            handle.write("bond %d:%d\n" % (atom_a, atom_b))
        for frame in frames:
            handle.write("\ntimestep ordered\n")
            for xyz in frame:
                handle.write("%.3f %.3f %.3f\n" % (xyz[0], xyz[1], xyz[2]))


def pdb_atom_line(record_type, serial, name, resname, chain, resid, xyz, element):
    return (
        "%-6s%5d %-4s %-3s %1s%4d    "
        "%8.3f%8.3f%8.3f  1.00  0.00           %2s"
    ) % (record_type, serial, name[:4], resname[:3], chain[:1], resid, xyz[0], xyz[1], xyz[2], element[:2])


def write_pdb(path, atom_records, frames):
    with open(path, "w") as handle:
        for model_id, frame in enumerate(frames, start=1):
            handle.write("MODEL     %4d\n" % model_id)
            for serial, (record, xyz) in enumerate(zip(atom_records, frame), start=1):
                handle.write(
                    pdb_atom_line(
                        record["record"],
                        serial,
                        record["name"],
                        record["resname"],
                        record["chain"],
                        record["resid"] + 1,
                        xyz,
                        record["element"],
                    ) + "\n"
                )
            handle.write("ENDMDL\n")
        handle.write("END\n")


def main():
    ensure_dir(results_dir())

    with tb.open_file(run_path()) as t:
        seq = [x.decode("ascii") for x in t.root.input.sequence[:]]
        ligand_group = t.root.input.ligand
        ligand_start = int(ligand_group._v_attrs.atom_start)
        ligand_names = [x.decode("ascii") for x in ligand_group.atom_name[:]]
        ligand_elements = [x.decode("ascii") for x in ligand_group.element[:]]
        ligand_bonds = [tuple(x) for x in ligand_group.bond_index[:]]
        pocket_residues = set(int(x) for x in ligand_group.pocket_residue[:])

        sc_node = t.root.input.potential.placement_fixed_point_vector_only
        affine_residue = sc_node.affine_residue[:]
        id_seq = sc_node.id_seq[:]
        rotamer_node = t.root.input.potential.ligand_rotamer_charmm_1body
        atom_start = rotamer_node.atom_start[:]
        atom_count = rotamer_node.atom_count[:]
        bond_start = rotamer_node.bond_start[:]
        bond_count = rotamer_node.bond_count[:]
        atom_name = [x.decode("ascii") for x in rotamer_node.atom_name[:]]
        atom_element = [x.decode("ascii") for x in rotamer_node.atom_element[:]]
        atom_local_pos = rotamer_node.atom_local_pos[:]
        bond_index = rotamer_node.bond_index[:]

        frames = latest_output_positions(t)

    engine = ue.Upside(run_path())

    n_res = len(seq)
    selector = (1 << N_BIT_ROTAMER) - 1
    pocket_slots = []
    for residue in sorted(pocket_residues):
        mask = affine_residue == residue
        residue_ids = id_seq[mask]
        rot_ids = residue_ids & selector
        unique_rot = sorted(set(int(x) for x in rot_ids))
        if not unique_rot:
            continue
        rep_index = int(np.flatnonzero(mask)[0])
        rep_start = int(atom_start[rep_index])
        rep_count = int(atom_count[rep_index])
        rep_bond_start = int(bond_start[rep_index])
        rep_bond_count = int(bond_count[rep_index])
        slot = {
            "residue": residue,
            "mask": mask,
            "rot_ids": rot_ids,
            "affine_residue": int(affine_residue[rep_index]),
            "atom_names": atom_name[rep_start: rep_start + rep_count],
            "atom_elements": atom_element[rep_start: rep_start + rep_count],
            "atom_local_pos": atom_local_pos[rep_start: rep_start + rep_count],
            "bond_index": bond_index[rep_bond_start: rep_bond_start + rep_bond_count],
        }
        residue_index = np.flatnonzero(mask)
        for idx in residue_index[1:]:
            idx = int(idx)
            start = int(atom_start[idx])
            count = int(atom_count[idx])
            if count != rep_count:
                raise RuntimeError("Atomistic sidechain library is not residue-consistent for residue %d" % residue)
            if atom_name[start: start + count] != slot["atom_names"]:
                raise RuntimeError("Atom naming mismatch across placements for residue %d" % residue)
        pocket_slots.append(slot)

    atom_records = []
    for resid, resname in enumerate(seq):
        atom_records.append({"record": "ATOM", "name": "N", "resname": resname, "chain": "A", "resid": resid, "segid": "PROT", "element": "N"})
        atom_records.append({"record": "ATOM", "name": "CA", "resname": resname, "chain": "A", "resid": resid, "segid": "PROT", "element": "C"})
        atom_records.append({"record": "ATOM", "name": "C", "resname": resname, "chain": "A", "resid": resid, "segid": "PROT", "element": "C"})

    bonds = backbone_bonds(n_res)
    for slot in pocket_slots:
        slot["output_atom_start"] = len(atom_records)
        for atom_name_i, element_i in zip(slot["atom_names"], slot["atom_elements"]):
            atom_records.append({
                "record": "ATOM",
                "name": atom_name_i,
                "resname": seq[slot["residue"]],
                "chain": "A",
                "resid": slot["residue"],
                "segid": "PCKT",
                "element": element_i,
            })
        if slot["atom_names"]:
            bonds.append((slot["residue"] * 3 + 1, slot["output_atom_start"]))
        for left, right in slot["bond_index"]:
            bonds.append((slot["output_atom_start"] + int(left), slot["output_atom_start"] + int(right)))

    ligand_offset = len(atom_records)
    for atom_name_i, element in zip(ligand_names, ligand_elements):
        atom_records.append({"record": "HETATM", "name": atom_name_i, "resname": "BNZ", "chain": "A", "resid": n_res, "segid": "LIG", "element": element})

    for left, right in ligand_bonds:
        bonds.append((ligand_offset + left, ligand_offset + right))

    out_frames = []
    for frame in frames:
        engine.energy(frame.astype("f4"))
        placement_prob = engine.get_value_by_name((len(affine_residue),), "rotamer", "placement_marginal")
        affine_output = engine.get_output("affine_alignment")

        coords = []
        coords.extend(frame[: 3 * n_res])

        for slot in pocket_slots:
            scores = []
            for rot_id in sorted(set(int(x) for x in slot["rot_ids"])):
                rot_mask = slot["mask"] & ((id_seq & selector) == rot_id)
                rot_prob = float(np.mean(placement_prob[rot_mask]))
                scores.append((rot_prob, rot_mask))
            _best_score, best_mask = max(scores, key=lambda item: item[0])
            rep_index = int(np.flatnonzero(best_mask)[0])
            aff = affine_output[int(affine_residue[rep_index])]
            if slot["atom_local_pos"].size:
                coords.extend(apply_affine(aff, slot["atom_local_pos"]))

        coords.extend(frame[ligand_start:])
        out_frames.append(np.asarray(coords, dtype="f4"))

    vtf_path = os.path.join(results_dir(), "4w52_ligand.vtf")
    pdb_path = os.path.join(results_dir(), "4w52_ligand.pdb")
    write_vtf(vtf_path, atom_records, bonds, out_frames)
    write_pdb(pdb_path, atom_records, out_frames)
    print("Wrote", vtf_path)
    print("Wrote", pdb_path)


if __name__ == "__main__":
    main()
