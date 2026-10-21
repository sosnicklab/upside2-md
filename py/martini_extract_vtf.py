#!/usr/bin/env python

import argparse
import os
import re
from pathlib import Path

import h5py
import numpy as np

PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"
OUTPUT_PREVIOUS_RE = re.compile(r"^output_previous_(\d+)$")

ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "P": 15,
    "S": 16,
    "CL": 17,
    "NA": 11,
}


def infer_atomic_number(atom_name, atom_type=None):
    for raw in (atom_name, atom_type):
        token = str(raw).strip().upper()
        if not token:
            continue
        if token in ATOMIC_NUMBERS:
            return ATOMIC_NUMBERS[token]
        if len(token) >= 2 and token[:2] in ATOMIC_NUMBERS:
            return ATOMIC_NUMBERS[token[:2]]
        if token[0] in ATOMIC_NUMBERS:
            return ATOMIC_NUMBERS[token[0]]
    return 6


def decode_str_array(dset):
    arr = dset[:]
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8").strip())
        else:
            out.append(str(x).strip())
    return np.array(out, dtype=object)


def find_martini_metadata_pdb(pdb_id, *up_files):
    """Locate <pdb_id>.MARTINI.pdb, which supplies lipid and ion residue labels for the VTF.

    The workflow writes it to <run_dir>/hybrid_prep/, not to the example's pdb/ directory, so looking only
    in pdb/ silently drops every environment resname to UNK and makes the trajectory unselectable by lipid
    in a viewer. Search the run directories implied by the .up files actually passed in, then pdb/. Every
    candidate is derived from an explicit argument, so no system's metadata can be attached to another's.
    """
    if not pdb_id:
        return ""
    name = f"{pdb_id}.MARTINI.pdb"
    candidates = []
    for up in up_files:
        if not up:
            continue
        up_dir = Path(up).expanduser().resolve().parent
        candidates.append(up_dir.parent / "hybrid_prep" / name)
        candidates.append(up_dir / "hybrid_prep" / name)
        candidates.append(up_dir / name)
    candidates.append(WORKFLOW_DIR / "pdb" / name)
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return ""


def normalize_frame(pos, n_particles):
    arr = np.asarray(pos)
    while arr.ndim > 2:
        arr = arr[0]
    if arr.ndim == 1:
        arr = arr.reshape(n_particles, 3)
    elif arr.ndim == 2:
        if arr.shape == (3, n_particles):
            arr = arr.T
        elif arr.shape != (n_particles, 3):
            raise ValueError(f"Unexpected 2D position shape: {arr.shape}")
    else:
        raise ValueError(f"Unexpected position shape: {arr.shape}")
    return np.asarray(arr, dtype=np.float32)


def normalize_input_positions(input_pos):
    arr = np.asarray(input_pos)
    if arr.ndim == 2:
        if arr.shape[1] == 3:
            return np.asarray(arr, dtype=np.float32)
        if arr.shape[0] == 3:
            return np.asarray(arr.T, dtype=np.float32)
    if arr.ndim == 3:
        # Common H5 layout from UPSIDE input: (n_atom, 3, 1)
        if arr.shape[1] == 3:
            return np.asarray(arr[:, :, 0], dtype=np.float32)
        # Alternate layout: (n_atom, 1, 3)
        if arr.shape[2] == 3:
            return np.asarray(arr[:, 0, :], dtype=np.float32)
        # Alternate layout: (3, n_atom, n_frame)
        if arr.shape[0] == 3:
            return np.asarray(arr[:, :, 0].T, dtype=np.float32)
    raise ValueError(f"Unexpected input/pos shape: {arr.shape}")


def read_martini_pdb_metadata(pdb_file):
    atom_names = []
    residue_names = []
    residue_ids = []
    chain_ids = []
    if not os.path.exists(pdb_file):
        return None, None, None, None

    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_names.append(line[12:16].strip())
            residue_names.append(line[17:21].strip())
            residue_ids.append(int(line[22:26]))
            chain_ids.append((line[21:22].strip() or "X"))

    return (
        np.array(atom_names, dtype=object),
        np.array(residue_names, dtype=object),
        np.array(residue_ids, dtype=int),
        np.array(chain_ids, dtype=object),
    )


def residue_names_from_class(atom_names, particle_class):
    """Residue names for a structure file that records none, from the particle class."""
    out = []
    for aname, pclass in zip(atom_names, particle_class):
        cls = str(pclass).upper()
        an = str(aname).upper()
        if cls == "PROTEIN":
            out.append("PRO")
        elif cls == "LIPID":
            out.append("LIP")
        elif cls == "ION":
            out.append(an if an in ("NA", "CL") else "ION")
        else:
            out.append("UNK")
    return np.array(out, dtype=object)


def chain_ids_from_class(particle_class):
    chain_by_class = {"PROTEIN": "A", "LIPID": "L", "ION": "I"}
    return np.array([chain_by_class.get(str(pclass).upper(), "X") for pclass in particle_class], dtype=object)


def group_boxes(traj_h5, struct_h5, output_group, n_frame):
    """The box of every frame: logged per frame under a barostat, otherwise the fixed box of the structure."""
    box_path = f"{output_group}/box"
    if output_group != "input" and box_path in traj_h5:
        return np.asarray(traj_h5[box_path][:], dtype=float).reshape(n_frame, 3)
    g = struct_h5["input/potential/martini_potential"] if "input/potential/martini_potential" in struct_h5 else None
    if g is None or not all(k in g.attrs for k in ("x_len", "y_len", "z_len")):
        raise ValueError(f"no box for {output_group}: neither {box_path} nor "
                         "input/potential/martini_potential x_len/y_len/z_len exists")
    return np.tile([float(g.attrs["x_len"]), float(g.attrs["y_len"]), float(g.attrs["z_len"])], (n_frame, 1))


def build_molecule_topology(n_output, bonds):
    """Bond-walk order and connected-component labels for the output particles.

    The walk is (child, parent) pairs with every parent placed before its children, so applying a
    minimum-image displacement in that order rebuilds a whole molecule from a single anchor atom.
    The labels group the particles into molecules; an unbonded particle is a molecule of one."""
    adjacency = [[] for _ in range(n_output)]
    for i, j in bonds:
        if 0 <= i < n_output and 0 <= j < n_output:
            adjacency[i].append(j)
            adjacency[j].append(i)

    walk = []
    molecule_id = np.full(n_output, -1, dtype=np.int64)
    n_molecules = 0
    for root in range(n_output):
        if molecule_id[root] >= 0:
            continue
        molecule_id[root] = n_molecules
        stack = [root]
        while stack:
            parent = stack.pop()
            for child in adjacency[parent]:
                if molecule_id[child] < 0:
                    molecule_id[child] = n_molecules
                    walk.append((child, parent))
                    stack.append(child)
        n_molecules += 1
    return walk, molecule_id, n_molecules


def unwrap_molecules(frame_pos, walk, x_len, y_len, z_len):
    """Make every bonded molecule whole again.

    The integrator stores each particle wrapped into the cell, so a molecule straddling a periodic
    face comes out torn and its bonds span most of the box when rendered. Walking the bond topology
    with a minimum-image displacement reassembles each molecule around its anchor atom."""
    box = np.array([x_len, y_len, z_len], dtype=float)
    out = np.array(frame_pos, copy=True)
    for child, parent in walk:
        out[child] -= box * np.round((out[child] - out[parent]) / box)
    return out


def centralize_system(frame_pos, residue_names, molecule_id, n_molecules, x_len, y_len, z_len):
    """Put the protein at the origin and bring every molecule back into the central cell.

    The molecules must already be whole (`unwrap_molecules`). Wrapping is then applied per molecule
    through its centroid, never per particle. Wrapping particles individually tears every molecule
    that crosses a face, and rebuilding them afterwards with the bond walk re-seats each one on
    whichever anchor atom happened to land inside: for glpG that anchor is the floppy N-terminal
    amide, which threw the entire protein a full box length out of the bilayer in 159 of 1822
    frames and left a protein-shaped hole behind."""
    out = np.array(frame_pos, copy=True)
    box = np.array([x_len, y_len, z_len], dtype=np.float64)

    protein_mask = None
    if residue_names is not None:
        residue_names = np.asarray(residue_names, dtype=object)
        protein_residues = {
            "PRO", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL",
            "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX",
        }
        protein_mask = np.array([str(name).upper() in protein_residues for name in residue_names], dtype=bool)

    if protein_mask is not None and np.any(protein_mask):
        out -= out[protein_mask].mean(axis=0)[None, :]

    counts = np.bincount(molecule_id, minlength=n_molecules).astype(np.float64)
    centroid = np.stack(
        [np.bincount(molecule_id, weights=out[:, ax], minlength=n_molecules) / counts for ax in range(3)],
        axis=1,
    )
    return out - box * np.round(centroid[molecule_id] / box)


def write_vtf_frame(fh, pos, box=None):
    fh.write("\ntimestep ordered\n")
    if box is not None:
        fh.write(f"pbc {box[0]} {box[1]} {box[2]}\n")
    for x, y, z in pos:
        fh.write(f"{x:.3f} {y:.3f} {z:.3f}\n")


def list_output_groups(traj_h5):
    groups = []
    for name in traj_h5.keys():
        match = OUTPUT_PREVIOUS_RE.fullmatch(str(name))
        if match and isinstance(traj_h5[name], h5py.Group):
            groups.append((int(match.group(1)), str(name)))
    groups.sort()
    ordered = [name for _, name in groups]
    if "output" in traj_h5 and isinstance(traj_h5["output"], h5py.Group):
        ordered.append("output")
    return ordered


def build_segment_output_path(output_file, segment_index):
    path = Path(output_file)
    return str(path.with_name(f"{path.stem}.segment_{segment_index}{path.suffix}"))


def collect_dist_spring_bonds(struct_h5):
    path = "input/potential/dist_spring/id"
    if path not in struct_h5:
        return np.zeros((0, 2), dtype=int)
    bonds = np.asarray(struct_h5[path][:], dtype=int)
    if bonds.ndim != 2 or bonds.shape[1] < 2:
        return np.zeros((0, 2), dtype=int)
    return bonds[:, :2]


def remap_bonds(bonds, index_map):
    out = []
    for i, j in bonds:
        ni = index_map.get(int(i))
        nj = index_map.get(int(j))
        if ni is None or nj is None:
            continue
        out.append((ni, nj))
    return out


def build_backbone_projection_map(struct_h5, input_pos):
    """The N/CA/C/O carriers of every residue, from /input/hybrid_bb_map (row i is Upside residue i).

    A file without a protein has no map and gets None; a map that does not resolve all four carriers of every
    residue is a broken build and is refused."""
    if "input/hybrid_bb_map" not in struct_h5:
        return None

    bb = struct_h5["input/hybrid_bb_map"]
    carriers = np.asarray(bb["atom_indices"][:], dtype=int)
    if carriers.ndim != 2 or carriers.shape[1] != 4:
        raise ValueError("hybrid_bb_map/atom_indices must have shape (n_residue, 4)")
    if np.any(carriers < 0) or np.any(carriers >= input_pos.shape[0]):
        raise ValueError("hybrid_bb_map/atom_indices contains an index outside input/pos")
    ref_atom_names = decode_str_array(bb["reference_atom_names"])
    atom_roles = decode_str_array(struct_h5["input/atom_roles"])
    if not np.all(atom_roles[carriers] == np.asarray(ref_atom_names, dtype=object).reshape(1, 4)):
        raise ValueError(
            "hybrid_bb_map/atom_indices do not point at %s for every residue" % list(ref_atom_names))
    sequence = decode_str_array(struct_h5["input/sequence"])
    if sequence.shape[0] != carriers.shape[0]:
        raise ValueError("hybrid_bb_map has %d residues but input/sequence has %d" % (
            carriers.shape[0], sequence.shape[0]))

    # The engine writes a residue's O slot only where infer_H_O has an acceptor on that residue's
    # carbonyl C (martini_hybrid.cpp, resolve_bb_o_hbond_elements). The C-terminal residue has none, so
    # its O slot keeps the input coordinate for the whole run and is not part of the trajectory.
    if "input/potential/infer_H_O/acceptors/id" not in struct_h5:
        raise ValueError("hybrid_bb_map requires /input/potential/infer_H_O/acceptors, which supplies the O slots")
    carbonyl_c = np.asarray(struct_h5["input/potential/infer_H_O/acceptors/id"][:, 1], dtype=int)
    carrier_present = np.ones(carriers.shape, dtype=bool)
    carrier_present[:, 3] = np.isin(carriers[:, 2], carbonyl_c)

    return {
        "bb_residue_index": np.asarray(bb["bb_residue_index"][:], dtype=int),
        "bb_residue_names": sequence,
        "bb_chain_ids": decode_str_array(bb["bb_chain_id"]),
        "ref_atom_names": np.array(ref_atom_names, dtype=object),
        "runtime_carrier_index": carriers,
        "carrier_present": carrier_present,
    }


def backbone_output_atoms(bb_map):
    """(atom name, residue name, residue id, chain id) of every backbone carrier the trajectory updates,
    in the order reconstruct_backbone_aa writes their coordinates."""
    for row, (resid, rname, chain_id) in enumerate(
        zip(bb_map["bb_residue_index"], bb_map["bb_residue_names"], bb_map["bb_chain_ids"])
    ):
        for k, aname in enumerate(bb_map["ref_atom_names"]):
            if bb_map["carrier_present"][row, k]:
                yield str(aname), str(rname), int(resid), str(chain_id).strip() or "A"


def build_mode1_mapping(
    struct_h5,
    input_pos,
    atom_names,
    residue_names,
    residue_ids,
    particle_class,
    pdb_metadata=None,
):
    # The backbone carriers are written once, after the MARTINI particles, by backbone_output_atoms.
    bb_map = build_backbone_projection_map(struct_h5, input_pos)
    keep = np.ones(input_pos.shape[0], dtype=bool)
    if bb_map is not None:
        keep[bb_map["runtime_carrier_index"].reshape(-1)] = False
    martini_indices = np.where(keep)[0]

    if (
        pdb_metadata is not None
        and pdb_metadata[0] is not None
        and pdb_metadata[0].shape[0] == martini_indices.shape[0]
    ):
        mart_atom_names, mart_res_names, mart_res_ids, mart_chain_ids = pdb_metadata
    else:
        mart_atom_names = np.asarray(atom_names, dtype=object)[martini_indices]
        mart_res_names = np.asarray(residue_names, dtype=object)[martini_indices]
        mart_res_ids = np.asarray(residue_ids, dtype=int)[martini_indices]
        mart_chain_ids = chain_ids_from_class(np.asarray(particle_class, dtype=object)[martini_indices])

    out_atom_names = list(mart_atom_names)
    out_atom_types = [str(name).strip() or "X" for name in mart_atom_names]
    out_res_names = list(mart_res_names)
    out_res_ids = [int(x) for x in np.asarray(mart_res_ids, dtype=int)]
    out_chain_ids = [str(x).strip() or "X" for x in np.asarray(mart_chain_ids, dtype=object)]
    out_atomic_numbers = [
        infer_atomic_number(aname, atype)
        for aname, atype in zip(out_atom_names, out_atom_types)
    ]
    include_aa_backbone = bb_map is not None
    if include_aa_backbone:
        for aname, rname, resid, chain_id in backbone_output_atoms(bb_map):
            out_atom_names.append(aname)
            out_atom_types.append(aname)
            out_res_names.append(rname)
            out_res_ids.append(resid)
            out_chain_ids.append(chain_id)
            out_atomic_numbers.append(infer_atomic_number(aname))

    return {
        "mode": 1,
        "martini_indices": martini_indices,
        "orig_to_output": {int(orig): int(out) for out, orig in enumerate(martini_indices)},
        "include_aa_backbone": include_aa_backbone,
        "bb_map": bb_map,
        "output_atom_names": np.array(out_atom_names, dtype=object),
        "output_atom_types": np.array(out_atom_types, dtype=object),
        "output_residue_names": np.array(out_res_names, dtype=object),
        "output_residue_ids": np.array(out_res_ids, dtype=int),
        "output_chain_ids": np.array(out_chain_ids, dtype=object),
        "output_atomic_numbers": np.array(out_atomic_numbers, dtype=int),
    }


def build_mode2_mapping(
    struct_h5,
    input_pos,
    atom_names,
    residue_names,
    residue_ids,
    particle_class,
    pdb_metadata=None,
):
    if "input/hybrid_env_topology/protein_membership" not in struct_h5:
        raise ValueError("Mode 2 requires /input/hybrid_env_topology/protein_membership")

    protein_membership = np.asarray(struct_h5["input/hybrid_env_topology/protein_membership"][:], dtype=int)
    n_particles = input_pos.shape[0]
    if protein_membership.shape[0] != n_particles:
        raise ValueError("protein_membership length mismatch")

    non_protein_idx = np.where(protein_membership < 0)[0]
    bb_map = build_backbone_projection_map(struct_h5, input_pos)
    if bb_map is None:
        raise ValueError("Mode 2 requires /input/hybrid_bb_map with valid backbone entries")

    if (
        pdb_metadata is not None
        and pdb_metadata[0] is not None
        and pdb_metadata[0].shape[0] == non_protein_idx.shape[0]
    ):
        env_atom_names, env_res_names, env_res_ids, env_chain_ids = pdb_metadata
    else:
        env_atom_names = np.asarray(atom_names, dtype=object)[non_protein_idx]
        env_res_names = np.asarray(residue_names, dtype=object)[non_protein_idx]
        env_res_ids = np.asarray(residue_ids, dtype=int)[non_protein_idx]
        env_chain_ids = chain_ids_from_class(np.asarray(particle_class, dtype=object)[non_protein_idx])

    out_atom_names = list(env_atom_names)
    out_atom_types = [str(name).strip() or "X" for name in env_atom_names]
    out_res_names = list(env_res_names)
    out_res_ids = [int(x) for x in env_res_ids]
    out_chain_ids = [str(x).strip() or "X" for x in np.asarray(env_chain_ids, dtype=object)]
    out_atomic_numbers = [
        infer_atomic_number(aname, atype)
        for aname, atype in zip(out_atom_names, out_atom_types)
    ]

    for aname, rname, resid, chain_id in backbone_output_atoms(bb_map):
        out_atom_names.append(aname)
        out_atom_types.append(aname)
        out_res_names.append(rname)
        out_res_ids.append(resid)
        out_chain_ids.append(chain_id)
        out_atomic_numbers.append(infer_atomic_number(aname))

    return {
        "mode": 2,
        "non_protein_idx": non_protein_idx,
        "orig_to_output": {int(orig): int(out) for out, orig in enumerate(non_protein_idx)},
        "bb_map": bb_map,
        "output_atom_names": np.array(out_atom_names, dtype=object),
        "output_atom_types": np.array(out_atom_types, dtype=object),
        "output_residue_names": np.array(out_res_names, dtype=object),
        "output_residue_ids": np.array(out_res_ids, dtype=int),
        "output_chain_ids": np.array(out_chain_ids, dtype=object),
        "output_atomic_numbers": np.array(out_atomic_numbers, dtype=int),
    }


def backbone_bonds(start_idx, bb_map):
    """N-CA, CA-C and C-O within each residue and C-N to the next residue of the same chain, indexed
    over the carriers backbone_output_atoms writes (a residue without an O slot has no C-O bond)."""
    present = bb_map["carrier_present"]
    out_index = np.full(present.shape, -1, dtype=int)
    out_index[present] = start_idx + np.arange(int(present.sum()))
    chain_ids = [str(c).strip() or " " for c in bb_map["bb_chain_ids"]]
    bonds = []
    for r in range(present.shape[0]):
        n, ca, c, o = out_index[r]
        bonds.append((n, ca))
        bonds.append((ca, c))
        if o >= 0:
            bonds.append((c, o))
        if r + 1 < present.shape[0] and chain_ids[r] == chain_ids[r + 1]:
            bonds.append((c, out_index[r + 1, 0]))
    return bonds


def reconstruct_backbone_aa(frame_pos, bb_map):
    """Backbone N/CA/C/O for the VTF, read straight out of the frame.

    These are real entries in the trajectory, so there is nothing to reconstruct. The previous fallback
    rigidly translated a stored reference geometry by the BB proxy's displacement, which could not
    represent any backbone rotation and only existed because the proxy predated the carriers being
    addressable. A build whose BB map does not resolve all four carriers is a broken build, not a case to
    approximate around. Only the carriers the trajectory updates are written (carrier_present).
    """
    carriers = bb_map["runtime_carrier_index"][bb_map["carrier_present"]]
    return np.asarray(frame_pos[carriers], dtype=np.float32).reshape(-1, 3)


def assemble_mode1_frame(frame_pos, mapping):
    martini_pos = frame_pos[mapping["martini_indices"]]
    if not mapping.get("include_aa_backbone", False):
        return martini_pos
    aa_pos = reconstruct_backbone_aa(frame_pos, mapping["bb_map"])
    return np.vstack((martini_pos, aa_pos))


def assemble_mode2_frame(frame_pos, mapping):
    env_pos = frame_pos[mapping["non_protein_idx"]]
    aa_pos = reconstruct_backbone_aa(frame_pos, mapping["bb_map"])
    return np.vstack((env_pos, aa_pos))


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract MARTINI trajectory into VTF.\n"
            "mode 1: all MARTINI particles + protein all-atom backbone (N,CA,C,O)\n"
            "mode 2: non-protein MARTINI particles + protein all-atom backbone (N,CA,C,O)"
        )
    )
    parser.add_argument("input_up", help="trajectory .up file")
    parser.add_argument("output_file", help="output .vtf file")
    parser.add_argument("structure_up", nargs="?", default=None, help="structure source .up (default: input_up)")
    parser.add_argument("pdb_id", nargs="?", default=None, help="PDB id for pdb/<id>.MARTINI.pdb metadata")
    parser.add_argument(
        "--mode",
        type=int,
        choices=(1, 2),
        default=None,
        help=(
            "output mode (default: auto). "
            "Auto selects mode 2 when the structure file contains "
            "input/hybrid_env_topology/protein_membership (hybrid systems), "
            "and mode 1 otherwise."
        ),
    )
    parser.add_argument(
        "--split-segments",
        action="store_true",
        help="Write one output file per output_previous_* / output trajectory segment instead of joining them "
             "into one trajectory.",
    )
    return parser.parse_args()


def extract_trajectory(traj_h5, struct_h5, input_pos, mapping, out_bonds, input_file, structure_file,
                       output_file, pdb_id, mode, output_groups):
    """Write the frames of `output_groups`, joined in order, to one VTF.

    The first frame of every group after the first repeats the frame that ended the group before it (a
    restart logs its starting state) and is dropped, as mdtraj_upside.load_upside_traj does. Each frame is
    made whole and wrapped in its own box; a barostat run writes that box into every timestep."""
    if output_file.split(".")[-1].lower() != "vtf":
        raise ValueError("Output file must be .vtf")
    n_particles = input_pos.shape[0]

    print(f"Extracting trajectory from: {input_file}")
    print(f"Structure source: {structure_file}")
    print(f"Output: {output_file}")
    print(f"Mode: {mode}")
    print(f"PDB ID: {pdb_id if pdb_id else '(none - no PDB metadata)'}")
    print(f"Trajectory groups: {', '.join(output_groups)}")

    frames = []   # (group, frame index, box)
    per_frame_box = False
    for g_no, output_group in enumerate(output_groups):
        if output_group == "input":
            n_frame = 1
        else:
            if "pos" not in traj_h5[output_group]:
                raise ValueError(f"Trajectory group has no positions: {output_group}")
            n_frame = int(traj_h5[output_group]["pos"].shape[0])
            per_frame_box |= f"{output_group}/box" in traj_h5
        boxes = group_boxes(traj_h5, struct_h5, output_group, n_frame)
        frames.extend((output_group, k, boxes[k]) for k in range(1 if g_no else 0, n_frame))

    n_output_particles = mapping['output_atom_names'].shape[0]
    print(f"Particles (input): {n_particles}")
    print(f"Particles (output): {n_output_particles}")
    print(f"Frames: {len(frames)}")
    print(f"Box: {frames[0][2][0]:.3f} {frames[0][2][1]:.3f} {frames[0][2][2]:.3f}")

    bond_walk, molecule_id, n_molecules = build_molecule_topology(n_output_particles, out_bonds)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# VTF extracted from UPSIDE MARTINI trajectory\n")
        f.write(f"# mode {mode}\n")
        f.write(f"# groups {' '.join(output_groups)}\n")
        for i, (aname, atype, rname, resid, chain_id, atomic_number) in enumerate(
            zip(
                mapping["output_atom_names"],
                mapping["output_atom_types"],
                mapping["output_residue_names"],
                mapping["output_residue_ids"],
                mapping["output_chain_ids"],
                mapping["output_atomic_numbers"],
            )
        ):
            chain = (str(chain_id).strip() or "X")[0]
            segid = f"s{chain}"
            f.write(
                f"atom {i} name {aname} type {atype} resid {int(resid)} "
                f"resname {str(rname)} segid {segid} chain {chain} "
                f"atomicnumber {int(atomic_number)}\n"
            )
        for i, j in out_bonds:
            f.write(f"bond {i}:{j}\n")
        first_box = frames[0][2]
        f.write(f"pbc {first_box[0]} {first_box[1]} {first_box[2]}\n")

        for frame_idx, (output_group, k, frame_box) in enumerate(frames):
            if output_group == "input":
                frame = input_pos
            else:
                frame = normalize_frame(traj_h5[output_group]["pos"][k], n_particles)

            if mode == 1:
                out_frame = assemble_mode1_frame(frame, mapping)
            else:
                out_frame = assemble_mode2_frame(frame, mapping)

            out_frame = unwrap_molecules(out_frame, bond_walk, *frame_box)
            out_frame = centralize_system(
                out_frame,
                mapping["output_residue_names"][:n_output_particles],
                molecule_id,
                n_molecules,
                *frame_box,
            )
            if np.isnan(out_frame).any():
                raise ValueError(f"NaN coordinates found in frame {frame_idx} ({output_group} frame {k})")
            write_vtf_frame(f, out_frame, frame_box if per_frame_box else None)
            if frame_idx % 100 == 0:
                print(f"Processed frame {frame_idx}/{len(frames) - 1}")


def main():
    args = parse_args()

    input_file = args.input_up
    output_file = args.output_file
    structure_file = args.structure_up or input_file
    pdb_id = args.pdb_id

    pdb_file = find_martini_metadata_pdb(pdb_id, input_file, structure_file)

    with h5py.File(input_file, "r") as t, h5py.File(structure_file, "r") as s:
        input_pos = normalize_input_positions(s["input/pos"][:])
        n_particles = input_pos.shape[0]

        atom_names = decode_str_array(s["input/atom_names"])
        residue_ids = np.asarray(s["input/residue_ids"][:], dtype=int)
        particle_class = decode_str_array(s["input/particle_class"])
        for name, values in (("atom_names", atom_names), ("residue_ids", residue_ids),
                             ("particle_class", particle_class)):
            if values.shape[0] != n_particles:
                raise ValueError(f"input/{name} length mismatch")
        if "input/residue_names" in s:
            residue_names = decode_str_array(s["input/residue_names"])
        else:
            residue_names = residue_names_from_class(atom_names, particle_class)

        pdb_atom_names, pdb_res_names, pdb_res_ids, pdb_chain_ids = read_martini_pdb_metadata(pdb_file)
        pdb_metadata = (pdb_atom_names, pdb_res_names, pdb_res_ids, pdb_chain_ids)

        is_hybrid = "input/hybrid_env_topology/protein_membership" in s
        if args.mode is None:
            mode = 2 if is_hybrid else 1
            print(f"Mode: auto-detected {'2 (hybrid)' if is_hybrid else '1 (pure MARTINI)'}")
        else:
            mode = args.mode

        build_mapping = build_mode1_mapping if mode == 1 else build_mode2_mapping
        mapping = build_mapping(s, input_pos, atom_names, residue_names, residue_ids, particle_class,
                                pdb_metadata=pdb_metadata)

        dist_bonds = collect_dist_spring_bonds(s)
        kept = mapping["martini_indices"] if mode == 1 else mapping["non_protein_idx"]
        idx_map = {int(old): int(new) for new, old in enumerate(kept)}
        out_bonds = remap_bonds(dist_bonds, idx_map)
        if mapping["bb_map"] is not None:
            out_bonds.extend(backbone_bonds(len(kept), mapping["bb_map"]))

        output_groups = list_output_groups(t) or ["input"]
        if args.split_segments and len(output_groups) > 1:
            print(f"Found {len(output_groups)} trajectory segments; writing one file per segment.")
            for segment_index, output_group in enumerate(output_groups):
                segment_output_file = build_segment_output_path(output_file, segment_index)
                extract_trajectory(t, s, input_pos, mapping, out_bonds, input_file, structure_file,
                                   segment_output_file, pdb_id, mode, [output_group])
                print(f"Wrote segment {segment_index}: {segment_output_file} ({output_group})")
        else:
            extract_trajectory(t, s, input_pos, mapping, out_bonds, input_file, structure_file,
                               output_file, pdb_id, mode, output_groups)


if __name__ == "__main__":
    main()
