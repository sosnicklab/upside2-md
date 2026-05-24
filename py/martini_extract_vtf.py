#!/usr/bin/env python3

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


def infer_pdb_id(input_file, user_pdb_id):
    if user_pdb_id:
        return user_pdb_id
    base = os.path.basename(input_file).lower()
    if "bilayer" in base:
        return "bilayer"
    return "1rkl"


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


def infer_residue_names_from_class(atom_names, particle_class):
    if particle_class is None:
        return np.array(["UNK"] * len(atom_names), dtype=object)

    lipid_atom_names = {
        "NC3", "PO4", "GL1", "GL2", "C1A", "C2A", "D3A", "C4A", "C5A",
        "C1B", "C2B", "D3B", "C4B", "C5B",
    }
    out = []
    for aname, pclass in zip(atom_names, particle_class):
        cls = str(pclass).upper()
        an = str(aname).upper()
        if cls == "PROTEIN":
            out.append("PRO")
        elif cls == "ION":
            if an == "NA":
                out.append("NA")
            elif an == "CL":
                out.append("CL")
            else:
                out.append("ION")
        elif cls == "OTHER" and an in lipid_atom_names:
            out.append("DOPC")
        else:
            out.append("UNK")
    return np.array(out, dtype=object)


def infer_chain_ids_from_class(atom_names, particle_class):
    if particle_class is None:
        return np.array(["X"] * len(atom_names), dtype=object)

    out = []
    for aname, pclass in zip(atom_names, particle_class):
        cls = str(pclass).upper()
        an = str(aname).upper()
        if cls == "PROTEIN" or cls == "PROTEINAA":
            out.append("A")
        elif cls == "ION":
            out.append("I")
        elif cls == "OTHER" and an in {"NC3", "PO4", "GL1", "GL2", "C1A", "C2A", "D3A", "C4A", "C5A", "C1B", "C2B", "D3B", "C4B", "C5B"}:
            out.append("L")
        else:
            out.append("X")
    return np.array(out, dtype=object)


def infer_box_lengths(traj_h5, struct_h5, pdb_file, input_file, output_group="output"):
    x_len = y_len = z_len = None

    box_path = f"{output_group}/box" if output_group else None
    if box_path and box_path in traj_h5:
        box_data = np.asarray(traj_h5[box_path][:])
        if box_data.size >= 3:
            last = box_data[-1]
            while np.asarray(last).ndim > 1:
                last = last[0]
            last = np.asarray(last).reshape(-1)
            if last.size >= 3:
                bx, by, bz = float(last[0]), float(last[1]), float(last[2])
                if bx > 0 and by > 0 and bz > 0:
                    x_len, y_len, z_len = bx, by, bz

    if x_len is None:
        log_file = input_file.replace(".up", ".log")
        if os.path.exists(log_file):
            with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
            for line in reversed(lines):
                if "[NPT]" in line and "box" in line:
                    m = re.search(r"box\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
                    if m:
                        x_len, y_len, z_len = float(m.group(1)), float(m.group(2)), float(m.group(3))
                    break

    if x_len is None and "input/potential/martini_potential" in struct_h5:
        g = struct_h5["input/potential/martini_potential"]
        if all(k in g.attrs for k in ("x_len", "y_len", "z_len")):
            x_len = float(g.attrs["x_len"])
            y_len = float(g.attrs["y_len"])
            z_len = float(g.attrs["z_len"])

    if x_len is None and os.path.exists(pdb_file):
        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("CRYST1"):
                    parts = line.split()
                    if len(parts) >= 4:
                        x_len = float(parts[1])
                        y_len = float(parts[2])
                        z_len = float(parts[3])
                    break

    if x_len is None:
        x_len = y_len = z_len = 50.0

    return x_len, y_len, z_len


def centralize_system(frame_pos, residue_names, x_len, y_len, z_len):
    out = np.array(frame_pos, copy=True)
    box = np.array([x_len, y_len, z_len], dtype=np.float32)
    half_box = 0.5 * box

    protein_mask = None
    if residue_names is not None:
        residue_names = np.asarray(residue_names, dtype=object)
        protein_residues = {
            "PRO", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL",
            "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX",
        }
        protein_mask = np.array([str(name).upper() in protein_residues for name in residue_names], dtype=bool)

    # Mode 2 output includes backmapped protein backbone residues by their
    # sequence-derived residue names. Centering by periodic protein COM avoids
    # boundary-split rendering artifacts.
    if protein_mask is not None and np.any(protein_mask):
        prot = np.mod(out[protein_mask], box[None, :])
        protein_center = np.zeros(3, dtype=np.float64)
        for ax in range(3):
            angles = 2.0 * np.pi * prot[:, ax] / float(box[ax])
            s = np.mean(np.sin(angles))
            c = np.mean(np.cos(angles))
            if abs(s) < 1e-12 and abs(c) < 1e-12:
                protein_center[ax] = float(np.mean(prot[:, ax]))
            else:
                a = np.arctan2(s, c)
                if a < 0.0:
                    a += 2.0 * np.pi
                protein_center[ax] = float(box[ax]) * a / (2.0 * np.pi)
        out -= protein_center[None, :]
    else:
        if residue_names is not None:
            lipid_mask = np.array([name == "DOPC" for name in residue_names], dtype=bool)
            if np.any(lipid_mask):
                out[:, 2] -= np.mean(out[lipid_mask, 2])

    out = (out + half_box) % box - half_box
    return out


def minimum_image_delta(delta, box_lengths):
    out = np.asarray(delta, dtype=np.float32)
    box = np.asarray(box_lengths, dtype=np.float32).reshape(1, 3)
    safe_box = np.where(box > 0.0, box, 1.0)
    return out - safe_box * np.round(out / safe_box)


def write_vtf_frame(fh, pos):
    fh.write("\ntimestep ordered\n")
    for x, y, z in pos:
        fh.write(f"{x:.3f} {y:.3f} {z:.3f}\n")


def build_cg_lipid_vector_info(struct_h5):
    """Read compose_vector6d data and build CG lipid display-vector info.

    Returns None if no CG lipids are present, else a dict with:
      - elem_indices: (n_cg,) int, original HDF5 indices of CG lipids
      - orientation_indices: (n_cg,) int or None, original HDF5 CGLD indices
      - direction: (n_cg, 3) float, static fallback unit direction vectors
      - head_offsets: (n_cg,) float, hydrophilic outer offsets in Angstrom
      - tail_offsets: (n_cg,) float, hydrophobic outer offsets in Angstrom
    """
    cv6_path = "input/potential/compose_vector6d"
    if cv6_path not in struct_h5:
        return None

    cv6 = struct_h5[cv6_path]
    if "elem_index" not in cv6 or "direction" not in cv6:
        return None

    elem_index = np.asarray(cv6["elem_index"][:], dtype=int)  # maps CG to pos atom index
    direction = np.asarray(cv6["direction"][:], dtype=np.float32)
    orientation_index = (
        np.asarray(cv6["orientation_index"][:], dtype=int)
        if "orientation_index" in cv6
        else None
    )
    fallback_span = float(cv6.attrs.get("orientation_length_ang", 5.0))
    metadata_head_offsets = (
        np.asarray(cv6["display_head_offset_ang"][:], dtype=np.float32)
        if "display_head_offset_ang" in cv6
        else np.full(elem_index.shape[0], -fallback_span, dtype=np.float32)
    )
    metadata_tail_offsets = (
        np.asarray(cv6["display_tail_offset_ang"][:], dtype=np.float32)
        if "display_tail_offset_ang" in cv6
        else np.full(elem_index.shape[0], fallback_span, dtype=np.float32)
    )

    if elem_index.size == 0:
        return None
    if metadata_head_offsets.shape[0] != elem_index.shape[0] or metadata_tail_offsets.shape[0] != elem_index.shape[0]:
        raise ValueError("CG lipid VTF endpoint-offset metadata length mismatch")
    display_span = float(np.mean(metadata_tail_offsets - metadata_head_offsets))
    if not np.isfinite(display_span) or display_span <= 0.0:
        display_span = float(cv6.attrs.get("head_tail_span_ang", 0.0))
    if not np.isfinite(display_span) or display_span <= 0.0:
        display_span = fallback_span

    return {
        "elem_indices": elem_index,
        "orientation_indices": orientation_index,
        "direction": direction,
        "head_offsets": metadata_head_offsets,
        "tail_offsets": metadata_tail_offsets,
        "display_span": display_span,
    }


def extend_with_lipid_vector_atoms(mapping, vector_info, out_bonds):
    """Convert visible CGL centers into two colored half-rods."""
    elem_idx = vector_info["elem_indices"]
    orig_to_output = mapping.get("orig_to_output", {})
    output_cg_idx = []
    source_ordinals = []
    for ordinal, original_idx in enumerate(elem_idx):
        output_idx = orig_to_output.get(int(original_idx))
        if output_idx is None:
            continue
        output_cg_idx.append(int(output_idx))
        source_ordinals.append(int(ordinal))
    output_cg_idx = np.asarray(output_cg_idx, dtype=int)
    source_ordinals = np.asarray(source_ordinals, dtype=int)
    n_cg = int(output_cg_idx.size)
    if n_cg == 0:
        vector_info["_tail_start"] = None
        vector_info["_n_vectors"] = 0
        return

    n_orig = mapping["output_atom_names"].shape[0]

    hydrophilic_names = np.array(["PO4"] * n_cg, dtype=object)
    hydrophilic_types = np.array(["Qa"] * n_cg, dtype=object)
    hydrophilic_resnames = np.array(["LIPH"] * n_cg, dtype=object)
    hydrophilic_atomic_numbers = np.full(n_cg, 15, dtype=int)
    hydrophobic_names = np.array(["C1A"] * n_cg, dtype=object)
    hydrophobic_types = np.array(["C1"] * n_cg, dtype=object)
    hydrophobic_resnames = np.array(["LIPT"] * n_cg, dtype=object)
    hydrophobic_atomic_numbers = np.full(n_cg, 6, dtype=int)
    lipid_resids = np.array(mapping["output_residue_ids"], dtype=int)[output_cg_idx]
    lipid_chain_ids = np.array(mapping["output_chain_ids"], dtype=object)[output_cg_idx]

    mapping["output_atom_names"][output_cg_idx] = hydrophilic_names
    mapping["output_atom_types"][output_cg_idx] = hydrophilic_types
    mapping["output_residue_names"][output_cg_idx] = hydrophilic_resnames
    mapping["output_atomic_numbers"][output_cg_idx] = hydrophilic_atomic_numbers

    mapping["output_atom_names"] = np.concatenate([
        mapping["output_atom_names"], hydrophilic_names, hydrophobic_names, hydrophobic_names,
    ])
    mapping["output_atom_types"] = np.concatenate([
        mapping["output_atom_types"], hydrophilic_types, hydrophobic_types, hydrophobic_types,
    ])
    mapping["output_residue_names"] = np.concatenate([
        mapping["output_residue_names"], hydrophilic_resnames, hydrophobic_resnames, hydrophobic_resnames,
    ])
    mapping["output_residue_ids"] = np.concatenate([
        mapping["output_residue_ids"], lipid_resids, lipid_resids, lipid_resids,
    ])
    mapping["output_chain_ids"] = np.concatenate([
        mapping["output_chain_ids"], lipid_chain_ids, lipid_chain_ids, lipid_chain_ids,
    ])
    mapping["output_atomic_numbers"] = np.concatenate([
        mapping["output_atomic_numbers"],
        hydrophilic_atomic_numbers,
        hydrophobic_atomic_numbers,
        hydrophobic_atomic_numbers,
    ])

    hydrophilic_center_start = n_orig
    hydrophobic_center_start = n_orig + n_cg
    hydrophobic_tail_start = n_orig + 2 * n_cg

    for gi in range(n_cg):
        head_atom = int(output_cg_idx[gi])
        hydrophilic_center_atom = hydrophilic_center_start + gi
        hydrophobic_center_atom = hydrophobic_center_start + gi
        tail_atom = hydrophobic_tail_start + gi
        head_name = str(mapping["output_atom_names"][head_atom]).strip().upper()
        if head_name != "PO4":
            raise ValueError(
                f"CG lipid hydrophilic endpoint remapped to unexpected atom {head_atom} ({head_name})"
            )
        out_bonds.append((head_atom, hydrophilic_center_atom))
        out_bonds.append((hydrophobic_center_atom, tail_atom))

    vector_info["_hydrophilic_center_start"] = hydrophilic_center_start
    vector_info["_hydrophobic_center_start"] = hydrophobic_center_start
    vector_info["_tail_start"] = hydrophobic_tail_start
    vector_info["_n_vectors"] = n_cg
    vector_info["_output_head_indices"] = output_cg_idx
    vector_info["_source_ordinals"] = source_ordinals

    print(f"CG lipid vectors: {n_cg} hydrophilic/hydrophobic half-rod pairs emitted")


def extend_frame_with_lipid_vectors(frame, vector_info, source_frame=None, box_lengths=None):
    """Replace CGL center positions with side-colored display rod atoms."""
    if vector_info is None or "_tail_start" not in vector_info:
        return frame
    if vector_info.get("_n_vectors", 0) == 0:
        return frame

    output_head_idx = vector_info["_output_head_indices"]
    source_ordinals = vector_info["_source_ordinals"]
    elem_idx = vector_info["elem_indices"][source_ordinals]
    orientation_idx = vector_info.get("orientation_indices")
    direction = vector_info["direction"][source_ordinals]
    head_offsets = vector_info["head_offsets"][source_ordinals]
    tail_offsets = vector_info["tail_offsets"][source_ordinals]

    if source_frame is not None and orientation_idx is not None:
        orient_idx = orientation_idx[source_ordinals]
        if np.max(orient_idx) < source_frame.shape[0] and np.max(elem_idx) < source_frame.shape[0]:
            direction = source_frame[orient_idx] - source_frame[elem_idx]
            if box_lengths is not None:
                direction = minimum_image_delta(direction, box_lengths)
        else:
            direction = np.asarray(direction, dtype=np.float32)
        norm = np.linalg.norm(direction, axis=1)
        mask = norm > 1e-12
        if np.any(mask):
            direction = direction.copy()
            direction[mask] /= norm[mask, None]

    out = np.asarray(frame, dtype=np.float32).copy()
    centers = out[output_head_idx].copy()
    out[output_head_idx] = centers + direction * head_offsets[:, None]
    hydrophilic_center_pos = centers.copy()
    hydrophobic_center_pos = centers.copy()
    tail_pos = centers + direction * tail_offsets[:, None]
    return np.concatenate([out, hydrophilic_center_pos, hydrophobic_center_pos, tail_pos], axis=0)


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
    if "input/hybrid_bb_map" not in struct_h5:
        return None

    bb = struct_h5["input/hybrid_bb_map"]
    if "bb_atom_index" not in bb or "bb_residue_index" not in bb or "reference_atom_coords" not in bb:
        return None

    bb_atom_index_raw = np.asarray(bb["bb_atom_index"][:], dtype=int)
    bb_residue_index = np.asarray(bb["bb_residue_index"][:], dtype=int)
    ref_coords = np.asarray(bb["reference_atom_coords"][:], dtype=np.float32)
    if ref_coords.ndim != 3 or ref_coords.shape[1:] != (4, 3):
        return None
    if bb_atom_index_raw.ndim != 1:
        return None
    n_bb = bb_atom_index_raw.shape[0]
    if bb_residue_index.shape[0] != n_bb or ref_coords.shape[0] != n_bb:
        return None

    if "reference_atom_names" in bb:
        ref_atom_names = decode_str_array(bb["reference_atom_names"])
    else:
        ref_atom_names = np.array(["N", "CA", "C", "O"], dtype=object)
    if "bb_chain_id" in bb:
        bb_chain_ids = decode_str_array(bb["bb_chain_id"])
    else:
        bb_chain_ids = np.array(["A"] * n_bb, dtype=object)
    bb_chain_ids = np.asarray(bb_chain_ids, dtype=object)
    if bb_chain_ids.shape[0] != n_bb:
        bb_chain_ids = np.array(["A"] * n_bb, dtype=object)

    n_particles = input_pos.shape[0]
    valid_proxy = (bb_atom_index_raw >= 0) & (bb_atom_index_raw < n_particles)

    atom_indices_full = None
    valid_carrier = np.zeros(n_bb, dtype=bool)
    if "atom_indices" in bb:
        atom_indices_candidate = np.asarray(bb["atom_indices"][:], dtype=int)
        if atom_indices_candidate.ndim == 2 and atom_indices_candidate.shape == (n_bb, 4):
            atom_indices_full = atom_indices_candidate
            valid_carrier = np.all((atom_indices_full >= 0) & (atom_indices_full < n_particles), axis=1)

    valid = valid_proxy | valid_carrier
    if not np.any(valid):
        return None

    if atom_indices_full is not None:
        fallback_ca = atom_indices_full[:, 1]
    else:
        fallback_ca = np.full(n_bb, -1, dtype=int)
    bb_atom_index = np.where(valid_proxy, bb_atom_index_raw, fallback_ca)
    bb_atom_index = bb_atom_index[valid]
    bb_residue_index = bb_residue_index[valid]
    ref_coords = ref_coords[valid]
    bb_chain_ids = bb_chain_ids[valid]

    bb_residue_names = np.array(["UNK"] * bb_residue_index.shape[0], dtype=object)
    if "input/sequence" in struct_h5:
        sequence = decode_str_array(struct_h5["input/sequence"])
        unique_residue_ids = []
        for resid in bb_residue_index.tolist():
            resid = int(resid)
            if resid not in unique_residue_ids:
                unique_residue_ids.append(resid)
        if len(sequence) == len(unique_residue_ids):
            residue_name_by_id = {
                int(resid): str(resname)
                for resid, resname in zip(unique_residue_ids, sequence.tolist())
            }
            bb_residue_names = np.array(
                [residue_name_by_id.get(int(resid), "UNK") for resid in bb_residue_index.tolist()],
                dtype=object,
            )

    if "weights" in bb:
        weights = np.asarray(bb["weights"][:], dtype=np.float32)
        if weights.ndim != 2 or weights.shape != (n_bb, 4):
            weights = np.full((bb_atom_index.shape[0], 4), 0.25, dtype=np.float32)
        else:
            weights = weights[valid]
    else:
        weights = np.full((bb_atom_index.shape[0], 4), 0.25, dtype=np.float32)

    wsum = np.sum(weights, axis=1, keepdims=True)
    safe_wsum = np.where(wsum > 1.0e-8, wsum, 1.0)
    weights = weights / safe_wsum
    ref_anchor = np.sum(ref_coords * weights[:, :, None], axis=1)

    runtime_carrier_index = None
    use_runtime_carriers = False
    if atom_indices_full is not None:
        atom_indices = atom_indices_full[valid]
        if np.all((atom_indices >= 0) & (atom_indices < n_particles)):
            runtime_carrier_index = atom_indices
            if "input/atom_roles" in struct_h5:
                atom_roles = decode_str_array(struct_h5["input/atom_roles"])
                runtime_roles = atom_roles[runtime_carrier_index]
                ref_role_row = np.asarray(ref_atom_names, dtype=object).reshape(1, 4)
                use_runtime_carriers = bool(np.all(runtime_roles == ref_role_row))
            else:
                use_runtime_carriers = True

    return {
        "bb_atom_index": bb_atom_index,
        "bb_residue_index": bb_residue_index,
        "bb_residue_names": bb_residue_names,
        "bb_chain_ids": bb_chain_ids,
        "ref_coords": ref_coords,
        "ref_atom_names": np.array(ref_atom_names, dtype=object),
        "weights": weights,
        "ref_anchor": ref_anchor,
        "runtime_carrier_index": runtime_carrier_index,
        "use_runtime_carriers": use_runtime_carriers,
    }


def build_mode1_mapping(
    struct_h5,
    input_pos,
    atom_names,
    residue_names,
    residue_ids,
    particle_class,
    pdb_metadata=None,
):
    n_particles = input_pos.shape[0]
    martini_indices = np.arange(n_particles, dtype=int)
    if particle_class is not None:
        keep_mask = np.array([pc != "PROTEINAA" for pc in particle_class], dtype=bool)
        if np.any(~keep_mask):
            martini_indices = np.where(keep_mask)[0]
    visible_mask = np.array(
        [str(atom_names[idx]).strip().upper() != "CGLD" for idx in martini_indices],
        dtype=bool,
    )
    if np.any(~visible_mask):
        martini_indices = martini_indices[visible_mask]

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
        mart_chain_ids = infer_chain_ids_from_class(
            np.asarray(atom_names, dtype=object)[martini_indices],
            None if particle_class is None else np.asarray(particle_class, dtype=object)[martini_indices],
        )

    bb_map = build_backbone_projection_map(struct_h5, input_pos)

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
        for resid, rname, chain_id in zip(
            bb_map["bb_residue_index"], bb_map["bb_residue_names"], bb_map["bb_chain_ids"]
        ):
            for aname in bb_map["ref_atom_names"]:
                out_atom_names.append(str(aname))
                out_atom_types.append(str(aname))
                out_res_names.append(str(rname))
                out_res_ids.append(int(resid))
                out_chain_ids.append(str(chain_id).strip() or "A")
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
    visible_mask = np.array(
        [str(atom_names[idx]).strip().upper() != "CGLD" for idx in non_protein_idx],
        dtype=bool,
    )
    if np.any(~visible_mask):
        non_protein_idx = non_protein_idx[visible_mask]
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
        env_chain_ids = infer_chain_ids_from_class(
            np.asarray(atom_names, dtype=object)[non_protein_idx],
            None if particle_class is None else np.asarray(particle_class, dtype=object)[non_protein_idx],
        )

    out_atom_names = list(env_atom_names)
    out_atom_types = [str(name).strip() or "X" for name in env_atom_names]
    out_res_names = list(env_res_names)
    out_res_ids = [int(x) for x in env_res_ids]
    out_chain_ids = [str(x).strip() or "X" for x in np.asarray(env_chain_ids, dtype=object)]
    out_atomic_numbers = [
        infer_atomic_number(aname, atype)
        for aname, atype in zip(out_atom_names, out_atom_types)
    ]

    for resid, rname, chain_id in zip(
        bb_map["bb_residue_index"], bb_map["bb_residue_names"], bb_map["bb_chain_ids"]
    ):
        for aname in bb_map["ref_atom_names"]:
            out_atom_names.append(str(aname))
            out_atom_types.append(str(aname))
            out_res_names.append(str(rname))
            out_res_ids.append(int(resid))
            out_chain_ids.append(str(chain_id).strip() or "A")
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


def mode2_backbone_bonds(start_idx, n_bb, bb_chain_ids=None):
    bonds = []
    for r in range(n_bb):
        b = start_idx + 4 * r
        bonds.append((b + 0, b + 1))
        bonds.append((b + 1, b + 2))
        bonds.append((b + 2, b + 3))
    for r in range(n_bb - 1):
        if bb_chain_ids is not None:
            chain_left = str(bb_chain_ids[r]).strip() or " "
            chain_right = str(bb_chain_ids[r + 1]).strip() or " "
            if chain_left != chain_right:
                continue
        b0 = start_idx + 4 * r
        b1 = start_idx + 4 * (r + 1)
        bonds.append((b0 + 2, b1 + 0))
    return bonds


def reconstruct_backbone_aa(frame_pos, bb_map, box_lengths=None):
    runtime_carrier_index = bb_map.get("runtime_carrier_index")
    if bb_map.get("use_runtime_carriers", False) and runtime_carrier_index is not None:
        aa_pos = frame_pos[runtime_carrier_index.reshape(-1)]
        return np.asarray(aa_pos, dtype=np.float32).reshape(-1, 3)

    bb_cur = frame_pos[bb_map["bb_atom_index"]]
    delta = bb_cur - bb_map["ref_anchor"]
    if box_lengths is not None:
        box = np.asarray(box_lengths, dtype=np.float32).reshape(1, 3)
        safe_box = np.where(box > 0.0, box, 1.0)
        delta = delta - safe_box * np.round(delta / safe_box)
    aa_pos = bb_map["ref_coords"] + delta[:, None, :]
    aa_pos = aa_pos.reshape(-1, 3)
    return aa_pos


def assemble_mode1_frame(frame_pos, mapping, box_lengths=None):
    martini_pos = frame_pos[mapping["martini_indices"]]
    if not mapping.get("include_aa_backbone", False):
        return martini_pos
    aa_pos = reconstruct_backbone_aa(frame_pos, mapping["bb_map"], box_lengths=box_lengths)
    return np.vstack((martini_pos, aa_pos))


def assemble_mode2_frame(frame_pos, mapping, box_lengths=None):
    env_pos = frame_pos[mapping["non_protein_idx"]]
    aa_pos = reconstruct_backbone_aa(frame_pos, mapping["bb_map"], box_lengths=box_lengths)
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
    parser.add_argument("--mode", type=int, choices=(1, 2), default=1, help="output mode (default: 1)")
    parser.add_argument(
        "--split-segments",
        action="store_true",
        help="Write one output file per output_previous_* / output trajectory segment instead of combining them.",
    )
    return parser.parse_args()


def extract_trajectory(
    traj_h5,
    struct_h5,
    input_pos,
    n_particles,
    mapping,
    out_bonds,
    input_file,
    structure_file,
    output_file,
    pdb_file,
    pdb_id,
    mode,
    output_group,
):
    fmt = output_file.split(".")[-1].lower()
    if fmt != "vtf":
        raise ValueError("Output file must be .vtf")

    print(f"Extracting trajectory from: {input_file}")
    print(f"Structure source: {structure_file}")
    print(f"Output: {output_file}")
    print(f"Mode: {mode}")
    print(f"PDB ID: {pdb_id}")
    print(f"Trajectory group: {output_group}")

    if output_group == "input":
        pos_data = None
        n_frame_total = 1
    else:
        if output_group not in traj_h5:
            raise ValueError(f"Trajectory group not found: {output_group}")
        group = traj_h5[output_group]
        if "pos" not in group:
            raise ValueError(f"Trajectory group has no positions: {output_group}")
        pos_data = group["pos"]
        n_frame_total = int(pos_data.shape[0])

    x_len, y_len, z_len = infer_box_lengths(
        traj_h5,
        struct_h5,
        pdb_file,
        input_file,
        output_group=output_group,
    )

    print(f"Particles (input): {n_particles}")
    n_output_particles = mapping['output_atom_names'].shape[0]
    print(f"Particles (output): {n_output_particles}")

    vector_info = build_cg_lipid_vector_info(struct_h5)
    if vector_info is not None:
        extend_with_lipid_vector_atoms(mapping, vector_info, out_bonds)
        print(f"Particles (output+lipid-vectors): {mapping['output_atom_names'].shape[0]}")
        span = vector_info["tail_offsets"] - vector_info["head_offsets"]
        print(
            "CG lipid display span: "
            f"mean={float(np.mean(span)):.3f} A "
            f"min={float(np.min(span)):.3f} A "
            f"max={float(np.max(span)):.3f} A"
        )

    print(f"Frames: {n_frame_total}")
    print(f"Box: {x_len:.3f} {y_len:.3f} {z_len:.3f}")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# VTF extracted from UPSIDE MARTINI trajectory\n")
        f.write(f"# mode {mode}\n")
        f.write(f"# group {output_group}\n")
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
        f.write(f"pbc {x_len} {y_len} {z_len}\n")

        for frame_idx in range(n_frame_total):
            if pos_data is None:
                frame = input_pos
            else:
                frame = normalize_frame(pos_data[frame_idx], n_particles)

            if mode == 1:
                out_frame = assemble_mode1_frame(frame, mapping, box_lengths=(x_len, y_len, z_len))
            else:
                out_frame = assemble_mode2_frame(frame, mapping, box_lengths=(x_len, y_len, z_len))

            out_frame = centralize_system(
                out_frame,
                mapping["output_residue_names"][:n_output_particles],
                x_len,
                y_len,
                z_len,
            )
            if vector_info is not None:
                out_frame = extend_frame_with_lipid_vectors(
                    out_frame,
                    vector_info,
                    source_frame=frame,
                    box_lengths=(x_len, y_len, z_len),
                )
            if np.isnan(out_frame).any():
                raise ValueError(f"NaN coordinates found in frame {frame_idx} for group {output_group}")
            write_vtf_frame(f, out_frame)
            if frame_idx % 100 == 0:
                print(f"Processed frame {frame_idx}/{n_frame_total - 1}")


def main():
    args = parse_args()

    input_file = args.input_up
    output_file = args.output_file
    structure_file = args.structure_up or input_file
    pdb_id = infer_pdb_id(input_file, args.pdb_id)
    mode = args.mode

    pdb_file = str(WORKFLOW_DIR / "pdb" / f"{pdb_id}.MARTINI.pdb")

    with h5py.File(input_file, "r") as t, h5py.File(structure_file, "r") as s:
        input_pos = normalize_input_positions(s["input/pos"][:])
        n_particles = input_pos.shape[0]

        if "input/atom_names" in s:
            atom_names = decode_str_array(s["input/atom_names"])
        else:
            atom_names = decode_str_array(s["input/type"])
        if atom_names.shape[0] != n_particles:
            raise ValueError("input/atom_names length mismatch")

        if "input/residue_ids" in s:
            residue_ids = np.asarray(s["input/residue_ids"][:], dtype=int)
            if residue_ids.shape[0] != n_particles:
                residue_ids = np.arange(1, n_particles + 1, dtype=int)
        else:
            residue_ids = np.arange(1, n_particles + 1, dtype=int)

        if "input/particle_class" in s:
            particle_class = decode_str_array(s["input/particle_class"])
            if particle_class.shape[0] != n_particles:
                particle_class = None
        else:
            particle_class = None

        inferred_res_names = infer_residue_names_from_class(atom_names, particle_class)

        pdb_atom_names, pdb_res_names, pdb_res_ids, pdb_chain_ids = read_martini_pdb_metadata(pdb_file)
        pdb_metadata = (pdb_atom_names, pdb_res_names, pdb_res_ids, pdb_chain_ids)

        if mode == 1:
            mapping = build_mode1_mapping(
                s,
                input_pos,
                atom_names,
                inferred_res_names,
                residue_ids,
                particle_class,
                pdb_metadata=pdb_metadata,
            )
        else:
            mapping = build_mode2_mapping(
                s,
                input_pos,
                atom_names,
                inferred_res_names,
                residue_ids,
                particle_class,
                pdb_metadata=pdb_metadata,
            )

        dist_bonds = collect_dist_spring_bonds(s)
        if mode == 1:
            mart_idx = mapping["martini_indices"]
            idx_map = {int(old): int(new) for new, old in enumerate(mart_idx)}
            out_bonds = remap_bonds(dist_bonds, idx_map)
            if mapping.get("include_aa_backbone", False):
                bb_start = len(mart_idx)
                out_bonds.extend(
                    mode2_backbone_bonds(
                        bb_start,
                        mapping["bb_map"]["bb_atom_index"].shape[0],
                        mapping["bb_map"].get("bb_chain_ids"),
                    )
                )
        else:
            non_idx = mapping["non_protein_idx"]
            idx_map = {int(old): int(new) for new, old in enumerate(non_idx)}
            out_bonds = remap_bonds(dist_bonds, idx_map)
            bb_start = len(non_idx)
            out_bonds.extend(
                mode2_backbone_bonds(
                    bb_start,
                    mapping["bb_map"]["bb_atom_index"].shape[0],
                    mapping["bb_map"].get("bb_chain_ids"),
                )
            )

        if args.split_segments:
            output_groups = list_output_groups(t)
            if len(output_groups) <= 1:
                target_group = output_groups[0] if output_groups else "input"
                extract_trajectory(
                    t,
                    s,
                    input_pos,
                    n_particles,
                    mapping,
                    out_bonds,
                    input_file,
                    structure_file,
                    output_file,
                    pdb_file,
                    pdb_id,
                    mode,
                    output_group=target_group,
                )
            else:
                print(f"Found {len(output_groups)} trajectory segments; writing one file per segment.")
                for segment_index, output_group in enumerate(output_groups):
                    segment_output_file = build_segment_output_path(output_file, segment_index)
                    extract_trajectory(
                        t,
                        s,
                        input_pos,
                        n_particles,
                        mapping,
                        out_bonds,
                        input_file,
                        structure_file,
                        segment_output_file,
                        pdb_file,
                        pdb_id,
                        mode,
                        output_group=output_group,
                    )
                    print(f"Wrote segment {segment_index}: {segment_output_file} ({output_group})")
        else:
            target_group = "output" if "output" in t else "input"
            extract_trajectory(
                t,
                s,
                input_pos,
                n_particles,
                mapping,
                out_bonds,
                input_file,
                structure_file,
                output_file,
                pdb_file,
                pdb_id,
                mode,
                output_group=target_group,
            )


if __name__ == "__main__":
    main()
