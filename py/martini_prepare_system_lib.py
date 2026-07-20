#!/usr/bin/env python

import _pickle as cPickle
import importlib.util
import os
import sys
import types
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path

import math

import h5py
import numpy as np
import tables as tb

from martini_itp_reader import parse_dry_forcefield, parse_itp_atomtype_masses, parse_itp_file

PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

NA_AVOGADRO = 6.02214076e23
BB_COMPONENT_NAMES = ("N", "CA", "C", "O")
BB_COMPONENT_MASSES = (14.0, 12.0, 12.0, 16.0)
MARTINIZE_BB_SS_ORDER = ("F", "E", "H", "1", "2", "3", "T", "S", "C")
MARTINIZE_BB_DEFAULT = ("N0", "Nda", "N0", "Nd", "Na", "Nda", "Nda", "P5", "P5")
MARTINIZE_BB_RESIDUE_OVERRIDE = {
    "ALA": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "P4", "P4"),
    "PRO": ("C5", "N0", "C5", "N0", "Na", "N0", "N0", "Na", "Na"),
    "HYP": ("C5", "N0", "C5", "N0", "N0", "N0", "N0", "Na", "Na"),
}
BB_TYPE_CHARGE = {
    "Qd": 1.0,
    "Qa": -1.0,
    "SQd": 1.0,
    "SQa": -1.0,
    "RQd": 1.0,
    "AQa": -1.0,
}
TWOPI = 2.0 * np.pi


CANONICAL_AFFINE_REF = np.array(
    [
        [-1.19280531, -0.83127186, 0.0],
        [0.0, 0.0, 0.0],
        [1.25222632, -0.87268266, 0.0],
    ],
    dtype=np.float32,
)
CANONICAL_AFFINE_REF -= CANONICAL_AFFINE_REF.mean(axis=0, keepdims=True)
STAGE_PARAMS = {
    "minimization": {
        "lj_soften": 1,
        "lj_alpha": 0.2,
        "coulomb_soften": 1,
        "slater_alpha": 2.0,
        "barostat_type": 0,
    },
    "npt_equil": {
        "lj_soften": 1,
        "lj_alpha": 0.2,
        "coulomb_soften": 1,
        "slater_alpha": 2.0,
        "barostat_type": 0,
    },
    "npt_equil_reduced": {
        "lj_soften": 1,
        "lj_alpha": 0.05,
        "coulomb_soften": 1,
        "slater_alpha": 0.5,
        "barostat_type": 0,
    },
    "npt_prod": {
        "lj_soften": 0,
        "lj_alpha": 0.0,
        "coulomb_soften": 0,
        "slater_alpha": 0.0,
        "barostat_type": 1,
    },
}


CB_PLACEMENT = np.array([[0.0, 0.94375626, 1.2068012]], dtype=np.float32)
CB_VECTOR = CB_PLACEMENT / np.linalg.norm(CB_PLACEMENT, axis=1, keepdims=True)
N_BIT_ROTAMER = 4
LEGACY_STAGE7_NODES = [
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
    "martini_sc_table_potential",
    "martini_sc_table_1body",
]
BACKBONE_NODES = [
    "Distance3D",
    "Angle",
    "Dihedral_omega",
    "Dihedral_phi",
    "Dihedral_psi",
    "Spring_bond",
    "Spring_angle",
    "Spring_omega",
    "rama_coord",
    "rama_map_pot",
    "rama_map_pot_ref",
    "infer_H_O",
    "protein_hbond",
    "hbond_energy",
    "backbone_pairs",
]


def parse_pdb(path: Path):
    atoms = []
    cryst1 = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            record = line[:6].strip()
            if record == "CRYST1":
                cryst1 = (
                    float(line[6:15]),
                    float(line[15:24]),
                    float(line[24:33]),
                )
                continue
            if record not in {"ATOM", "HETATM"}:
                continue

            atom = {
                "record": record,
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:21].strip(),
                "chain": line[21].strip(),
                "resseq": int(line[22:26]),
                "icode": line[26].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "occ": float(line[54:60]) if line[54:60].strip() else 1.0,
                "bfac": float(line[60:66]) if line[60:66].strip() else 0.0,
                "segid": line[72:76].strip() if len(line) >= 76 else "",
                "element": line[76:78].strip() if len(line) >= 78 else "",
                "charge": line[78:80].strip() if len(line) >= 80 else "",
            }
            atoms.append(atom)
    return atoms, cryst1


def write_pdb(path: Path, atoms, box_lengths):
    with path.open("w", encoding="utf-8") as f:
        if box_lengths is not None:
            f.write(
                f"CRYST1{box_lengths[0]:9.3f}{box_lengths[1]:9.3f}{box_lengths[2]:9.3f}"
                "  90.00  90.00  90.00 P 1           1\n"
            )
        for idx, atom in enumerate(atoms, start=1):
            chain = atom["chain"] if atom["chain"] else " "
            segid = atom["segid"][:4].ljust(4) if atom["segid"] else "    "
            f.write(
                f"{atom['record']:<6}{idx:5d} {atom['name'][:4]:>4} {atom['resname'][:4]:>4}{chain:1}"
                f"{atom['resseq']:4d}{atom['icode'][:1]:1}   "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"{atom['occ']:6.2f}{atom['bfac']:6.2f}      {segid}"
                f"{atom['element'][:2]:>2}{atom['charge'][:2]:>2}\n"
            )
        f.write("END\n")


def _minimum_image_delta(delta, box_lengths):
    delta = np.asarray(delta, dtype=np.float64).copy()
    if box_lengths is None:
        return delta
    box = np.asarray(box_lengths, dtype=np.float64)
    valid = box > 0.0
    delta[..., valid] -= box[valid] * np.round(delta[..., valid] / box[valid])
    return delta


def coords(atoms):
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def set_coords(atoms, xyz):
    for atom, c in zip(atoms, xyz):
        atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])


def center_of_mass(xyz):
    return np.mean(xyz, axis=0)


def lipid_resname(resname: str) -> bool:
    return resname.upper() in {"DOP", "DOPC"}


def canonical_lipid_resname(resname: str) -> str:
    if lipid_resname(resname):
        return "DOPC"
    return resname


def infer_protein_charge_from_residues(protein_atoms):
    charged_res = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
    seen = set()
    total = 0
    for atom in protein_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen:
            continue
        seen.add(key)
        total += charged_res.get(atom["resname"].upper(), 0)
    return total


def residue_key(atom):
    return (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].upper())


def extract_protein_backbone_atoms_from_aa(aa_atoms):
    protein_res = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "CYX",
    }
    backbone_roles = {"N", "CA", "C", "O"}
    by_residue = defaultdict(dict)
    for atom in aa_atoms:
        resname = atom["resname"].upper()
        if resname not in protein_res:
            continue
        role = atom["name"].strip().upper()
        if role not in backbone_roles:
            continue
        key = residue_key(atom)
        if role not in by_residue[key]:
            by_residue[key][role] = deepcopy(atom)

    ordered_keys = sorted(by_residue.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    out = []
    residue_index = []
    for rid, key in enumerate(ordered_keys, start=1):
        role_map = by_residue[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        for role in BB_COMPONENT_NAMES:
            atom = role_map[role]
            atom["name"] = role
            atom["segid"] = "PROA"
            out.append(atom)
            residue_index.append(rid)
    if not out:
        raise ValueError("No complete protein backbone residues (N/CA/C/O) found in AA PDB.")
    return out, np.asarray(residue_index, dtype=np.int32)


def _complete_backbone_residue_groups(backbone_atoms):
    residue_groups = defaultdict(dict)
    residue_order = []
    for atom in backbone_atoms:
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        if key not in residue_groups:
            residue_order.append(key)
        residue_groups[key][role] = atom

    valid_order = []
    valid_groups = {}
    for key in residue_order:
        group = residue_groups[key]
        if not all(role in group for role in BB_COMPONENT_NAMES):
            continue
        valid_order.append(key)
        valid_groups[key] = group
    if not valid_order:
        raise ValueError("Cannot map BB types: no complete N/CA/C/O residue rows.")
    return valid_order, valid_groups


def _atom_xyz(atom):
    return np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=np.float64)


def _dihedral_degrees(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    norm_b1 = np.linalg.norm(b1)
    if norm_b1 < 1.0e-8:
        return None
    b1 = b1 / norm_b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    norm_v = np.linalg.norm(v)
    norm_w = np.linalg.norm(w)
    if norm_v < 1.0e-8 or norm_w < 1.0e-8:
        return None
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def _same_chain_adjacent(prev_key, curr_key):
    return curr_key[0] == prev_key[0] and curr_key[1] == prev_key[1] + 1


def _raw_secondary_from_geometry(residue_order, residue_groups):
    raw = {key: "C" for key in residue_order}
    ca_pos = {key: _atom_xyz(residue_groups[key]["CA"]) for key in residue_order}
    raw_list = ["C"] * len(residue_order)

    for i, key in enumerate(residue_order):
        helix_like = False
        if 0 < i < len(residue_order) - 1:
            prev_key = residue_order[i - 1]
            next_key = residue_order[i + 1]
            if _same_chain_adjacent(prev_key, key) and _same_chain_adjacent(key, next_key):
                prev_group = residue_groups[prev_key]
                curr_group = residue_groups[key]
                next_group = residue_groups[next_key]
                phi = _dihedral_degrees(
                    _atom_xyz(prev_group["C"]),
                    _atom_xyz(curr_group["N"]),
                    _atom_xyz(curr_group["CA"]),
                    _atom_xyz(curr_group["C"]),
                )
                psi = _dihedral_degrees(
                    _atom_xyz(curr_group["N"]),
                    _atom_xyz(curr_group["CA"]),
                    _atom_xyz(curr_group["C"]),
                    _atom_xyz(next_group["N"]),
                )
                if phi is not None and psi is not None and -100.0 <= phi <= -30.0 and -85.0 <= psi <= 15.0:
                    helix_like = True

        for offset in (-4, 4):
            j = i + offset
            if helix_like or j < 0 or j >= len(residue_order):
                continue
            other_key = residue_order[j]
            if other_key[0] != key[0]:
                continue
            dist = float(np.linalg.norm(ca_pos[key] - ca_pos[other_key]))
            if 4.8 <= dist <= 6.8:
                helix_like = True

        if helix_like:
            raw_list[i] = "H"

    for i in range(1, len(raw_list) - 1):
        if (
            raw_list[i] == "C"
            and raw_list[i - 1] == "H"
            and raw_list[i + 1] == "H"
            and residue_order[i - 1][0] == residue_order[i][0] == residue_order[i + 1][0]
        ):
            raw_list[i] = "H"

    i = 0
    while i < len(raw_list):
        if raw_list[i] != "H":
            i += 1
            continue
        j = i + 1
        while j < len(raw_list) and raw_list[j] == "H" and residue_order[j][0] == residue_order[i][0]:
            j += 1
        if j - i < 4:
            for k in range(i, j):
                raw_list[k] = "C"
        i = j

    for key, ss in zip(residue_order, raw_list):
        raw[key] = ss
    return raw


def _apply_martinize_helix_labels(residue_order, raw_secondary):
    labels = {key: raw_secondary.get(key, "C") for key in residue_order}
    i = 0
    while i < len(residue_order):
        key = residue_order[i]
        if raw_secondary.get(key, "C") != "H":
            i += 1
            continue
        j = i + 1
        while (
            j < len(residue_order)
            and residue_order[j][0] == residue_order[i][0]
            and raw_secondary.get(residue_order[j], "C") == "H"
        ):
            j += 1
        run_len = j - i
        if run_len < 4:
            for k in range(i, j):
                labels[residue_order[k]] = "C"
        else:
            for offset, k in enumerate(range(i, j)):
                if offset < min(4, run_len):
                    labels[residue_order[k]] = "1"
                elif run_len - offset <= 4:
                    labels[residue_order[k]] = "2"
                else:
                    labels[residue_order[k]] = "H"
        i = j
    return labels


def _martinize_bb_type_for_residue(resname, secondary_structure):
    table = MARTINIZE_BB_RESIDUE_OVERRIDE.get(resname.upper(), MARTINIZE_BB_DEFAULT)
    try:
        return table[MARTINIZE_BB_SS_ORDER.index(str(secondary_structure).upper())]
    except ValueError:
        return table[MARTINIZE_BB_SS_ORDER.index("C")]


def map_backbone_types_from_structure(backbone_atoms):
    residue_order, residue_groups = _complete_backbone_residue_groups(backbone_atoms)
    raw_secondary = _raw_secondary_from_geometry(residue_order, residue_groups)
    secondary_by_residue = _apply_martinize_helix_labels(residue_order, raw_secondary)
    type_by_residue = {}
    source_by_residue = {}
    for key in residue_order:
        ss = secondary_by_residue[key]
        type_by_residue[key] = _martinize_bb_type_for_residue(key[3], ss)
        source_by_residue[key] = "structure_geometry" if ss in {"H", "1", "2", "3"} else "coil_geometry_fallback"

    if residue_order:
        terminal_overrides = [(residue_order[0], "Qd"), (residue_order[-1], "Qa")]
        for i in range(1, len(residue_order)):
            prev = residue_order[i - 1]
            curr = residue_order[i]
            chain_break = curr[0] != prev[0]
            seq_break = curr[1] != (prev[1] + 1)
            if chain_break or seq_break:
                terminal_overrides.append((curr, "Qd"))
                terminal_overrides.append((prev, "Qa"))
        for key, bb_type in terminal_overrides:
            type_by_residue[key] = bb_type
            source_by_residue[key] = "terminal_charge_override"

    return type_by_residue, secondary_by_residue, source_by_residue


def map_backbone_types_from_martinize_fallback(backbone_atoms):
    type_by_residue, _, _ = map_backbone_types_from_structure(backbone_atoms)
    return type_by_residue


def compute_lipid_residue_indices(bilayer_atoms):
    lipid_residues = defaultdict(list)
    keep_nonlipid = []
    for idx, atom in enumerate(bilayer_atoms):
        if lipid_resname(atom["resname"]):
            key = (atom["chain"], atom["resseq"], atom["icode"])
            lipid_residues[key].append(idx)
        else:
            keep_nonlipid.append(idx)
    return lipid_residues, keep_nonlipid


def residue_group_atoms(atoms):
    groups = defaultdict(list)
    for atom in atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].upper())
        groups[key].append(atom)
    return groups


def tile_and_crop_bilayer_lipids(
    bilayer_atoms,
    bilayer_box,
    target_xy_min,
    target_xy_max,
):
    # Keep only lipid residues for tiling; ions are re-generated after packing.
    lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not lipid_atoms:
        raise ValueError("No lipid residues found in bilayer input for tiling/cropping.")

    lip_xyz = coords(lipid_atoms)
    lip_center = center_of_mass(lip_xyz)

    if bilayer_box is not None:
        tile_x, tile_y = float(bilayer_box[0]), float(bilayer_box[1])
    else:
        lip_min = lip_xyz.min(axis=0)
        lip_max = lip_xyz.max(axis=0)
        span = lip_max - lip_min
        tile_x, tile_y = float(span[0]), float(span[1])
    if tile_x <= 0.0 or tile_y <= 0.0:
        raise ValueError("Invalid bilayer tile dimensions for XY tiling.")

    target_xy_min = np.array(target_xy_min, dtype=float)
    target_xy_max = np.array(target_xy_max, dtype=float)
    target_center = 0.5 * (target_xy_min + target_xy_max)
    target_span = target_xy_max - target_xy_min

    # Symmetric tiling around the bilayer center.
    # +2 guard ensures crop coverage despite residue COM discretization.
    nx = int(np.ceil(target_span[0] / tile_x)) + 2
    ny = int(np.ceil(target_span[1] / tile_y)) + 2

    base_groups = residue_group_atoms(lipid_atoms)
    tiled_groups = []
    for ix in range(-nx, nx + 1):
        for iy in range(-ny, ny + 1):
            shift = np.array([ix * tile_x, iy * tile_y, 0.0], dtype=float)
            for group in base_groups.values():
                gcopy = []
                for atom in group:
                    a = deepcopy(atom)
                    a["resname"] = canonical_lipid_resname(a["resname"])
                    a["x"] = float(a["x"] + shift[0])
                    a["y"] = float(a["y"] + shift[1])
                    gcopy.append(a)
                tiled_groups.append(gcopy)

    # Recenter tiled lattice so its center matches target center.
    tiled_xyz = np.array(
        [[a["x"], a["y"], a["z"]] for group in tiled_groups for a in group],
        dtype=float,
    )
    tiled_center_xy = center_of_mass(tiled_xyz)[:2]
    recenter_shift_xy = target_center - tiled_center_xy
    for group in tiled_groups:
        for atom in group:
            atom["x"] = float(atom["x"] + recenter_shift_xy[0])
            atom["y"] = float(atom["y"] + recenter_shift_xy[1])

    # Crop by residue center to preserve complete molecules.
    cropped_groups = []
    for group in tiled_groups:
        gxyz = np.array([[a["x"], a["y"]] for a in group], dtype=float)
        gcom = gxyz.mean(axis=0)
        if (
            target_xy_min[0] <= gcom[0] <= target_xy_max[0]
            and target_xy_min[1] <= gcom[1] <= target_xy_max[1]
        ):
            cropped_groups.append(group)

    if not cropped_groups:
        raise RuntimeError("Bilayer tiling/cropping produced no lipid residues.")

    # Renumber residues for unique molecule ids downstream.
    out_atoms = []
    next_resseq = 1
    for group in cropped_groups:
        chain = group[0]["chain"]
        icode = group[0]["icode"]
        for atom in group:
            atom["chain"] = chain
            atom["icode"] = icode
            atom["resseq"] = next_resseq
            atom["resname"] = canonical_lipid_resname(atom["resname"])
            atom["segid"] = "MEMB"
            out_atoms.append(atom)
        next_resseq += 1
    return out_atoms


def remove_overlapping_lipids(
    bilayer_atoms,
    protein_atoms,
    lipid_residues,
    keep_nonlipid,
    cutoff,
):
    protein_xyz = coords(protein_atoms)
    cutoff2 = cutoff * cutoff
    keep_atom_idx = set(keep_nonlipid)
    removed_residues = 0
    for res_key, atom_idx_list in lipid_residues.items():
        lipid_xyz = np.array(
            [[bilayer_atoms[i]["x"], bilayer_atoms[i]["y"], bilayer_atoms[i]["z"]] for i in atom_idx_list],
            dtype=float,
        )
        delta = lipid_xyz[:, None, :] - protein_xyz[None, :, :]
        dist2 = np.sum(delta * delta, axis=2)
        if np.min(dist2) < cutoff2:
            removed_residues += 1
            continue
        for i in atom_idx_list:
            keep_atom_idx.add(i)
    kept = [bilayer_atoms[i] for i in sorted(keep_atom_idx)]
    return kept, removed_residues


def set_box_from_lipid_xy(
    all_atoms,
    lipid_atoms,
    pad_z,
    force_square_xy=True,
    min_box_z=None,
    center_lipid_in_z=True,
):
    all_xyz = coords(all_atoms)
    lip_xyz = coords(lipid_atoms)
    lip_min = lip_xyz.min(axis=0)
    lip_max = lip_xyz.max(axis=0)

    # Force XY box edges to coincide with bilayer edges.
    span_x = float(lip_max[0] - lip_min[0])
    span_y = float(lip_max[1] - lip_min[1])
    box_x = span_x
    box_y = span_y
    if force_square_xy:
        square_side = max(span_x, span_y)
        box_x = square_side
        box_y = square_side
    if span_x <= 0.0 or span_y <= 0.0:
        raise ValueError("Invalid lipid XY span while defining box edges.")

    min_z = float(all_xyz[:, 2].min())
    max_z = float(all_xyz[:, 2].max())
    lip_mid_z = 0.5 * float(lip_min[2] + lip_max[2])
    span_z = max_z - min_z
    box_z = float(span_z + 2.0 * pad_z)
    if min_box_z is not None:
        box_z = float(max(box_z, float(min_box_z)))
    if center_lipid_in_z:
        # Keep lipid bilayer centered around box_z/2 while ensuring all atoms remain inside [0, box_z].
        needed_z = 2.0 * max(lip_mid_z - min_z, max_z - lip_mid_z)
        box_z = float(max(box_z, needed_z))

    shift = np.array([0.0, 0.0, 0.0], dtype=float)
    shift[0] = -lip_min[0]
    shift[1] = -lip_min[1]
    if center_lipid_in_z:
        shift[2] = 0.5 * box_z - lip_mid_z
    else:
        shift[2] = pad_z - min_z

    shifted = all_xyz + shift
    set_coords(all_atoms, shifted)
    return np.array([box_x, box_y, box_z], dtype=float)


def estimate_salt_pairs(box_lengths, salt_molar, effective_volume_fraction=1.0):
    volume_a3 = float(box_lengths[0] * box_lengths[1] * box_lengths[2])
    volume_a3 *= float(max(0.0, min(1.0, effective_volume_fraction)))
    volume_l = volume_a3 * 1e-27
    pairs = int(round(salt_molar * NA_AVOGADRO * volume_l))
    return max(0, pairs)


def infer_effective_ion_volume_fraction_from_template(
    bilayer_atoms,
    bilayer_box,
    salt_molar,
):
    # Calibrate effective ion-accessible volume from the template system ions.
    # This avoids overestimating ion count by using full geometric box volume.
    if bilayer_box is None or salt_molar <= 0.0:
        return 1.0

    n_na = sum(1 for a in bilayer_atoms if a["resname"].upper() == "NA")
    n_cl = sum(1 for a in bilayer_atoms if a["resname"].upper() == "CL")
    base_pairs = min(n_na, n_cl)
    if base_pairs <= 0:
        return 1.0

    base_volume_l = float(bilayer_box[0] * bilayer_box[1] * bilayer_box[2]) * 1e-27
    expected_pairs = salt_molar * NA_AVOGADRO * base_volume_l
    if expected_pairs <= 0.0:
        return 1.0

    frac = float(base_pairs) / float(expected_pairs)
    return float(max(0.0, min(1.0, frac)))


def place_ions(atoms, box_lengths, n_na, n_cl, cutoff, rng):
    existing = coords(atoms)
    placed = []
    cutoff2 = cutoff * cutoff
    box = np.asarray(box_lengths, dtype=float)
    types = [("NA", "NA", 1.0, "IONS")] * n_na + [("CL", "CL", -1.0, "IONS")] * n_cl
    for name, resname, _charge, segid in types:
        accepted = False
        for _ in range(20000):
            trial = rng.uniform([0, 0, 0], box)
            if existing.size:
                d2 = np.sum(_minimum_image_delta(existing - trial, box) ** 2, axis=1)
                if np.min(d2) < cutoff2:
                    continue
            if placed:
                placed_xyz = np.array([[a["x"], a["y"], a["z"]] for a in placed], dtype=float)
                d2_placed = np.sum(_minimum_image_delta(placed_xyz - trial, box) ** 2, axis=1)
                if np.min(d2_placed) < cutoff2:
                    continue
            placed.append(
                {
                    "record": "HETATM",
                    "serial": 0,
                    "name": name,
                    "resname": resname,
                    "chain": "I",
                    "resseq": len(placed) + 1,
                    "icode": "",
                    "x": float(trial[0]),
                    "y": float(trial[1]),
                    "z": float(trial[2]),
                    "occ": 1.0,
                    "bfac": 0.0,
                    "segid": segid,
                    "element": name[:2],
                    "charge": "",
                }
            )
            accepted = True
            break
        if not accepted:
            raise RuntimeError("Failed to place ions without overlaps; relax cutoff or enlarge box.")
    return placed


def build_backbone_with_virtual_bb(
    backbone_atoms,
    bb_type_by_residue,
    bb_secondary_by_residue=None,
    bb_type_source_by_residue=None,
):
    residue_groups = defaultdict(dict)
    for atom_idx, atom in enumerate(backbone_atoms):
        key = residue_key(atom)
        role = atom["name"].strip().upper()
        residue_groups[key][role] = atom_idx

    ordered_keys = sorted(residue_groups.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    protein_atoms = [deepcopy(atom) for atom in backbone_atoms]
    bb_entries = []
    for seq_idx, key in enumerate(ordered_keys, start=1):
        role_map = residue_groups[key]
        if not all(role in role_map for role in BB_COMPONENT_NAMES):
            continue
        idxs = [int(role_map[role]) for role in BB_COMPONENT_NAMES]
        masses = np.asarray(BB_COMPONENT_MASSES, dtype=np.float32)
        weights = (masses / masses.sum()).tolist()
        coords_row = []
        com = np.zeros(3, dtype=np.float64)
        wsum = 0.0
        for role, weight in zip(BB_COMPONENT_NAMES, weights):
            atom = backbone_atoms[role_map[role]]
            xyz = np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=np.float64)
            coords_row.append([float(xyz[0]), float(xyz[1]), float(xyz[2])])
            com += float(weight) * xyz
            wsum += float(weight)
        if wsum <= 0.0:
            raise ValueError("Backbone BB COM weights sum to zero.")
        if abs(wsum - 1.0) > 1.0e-8:
            com /= wsum

        bb_type = str(bb_type_by_residue.get(key, "P5")).strip()
        bb_secondary = str((bb_secondary_by_residue or {}).get(key, "C")).strip() or "C"
        bb_type_source = str((bb_type_source_by_residue or {}).get(key, "structure_default")).strip()
        proxy = deepcopy(backbone_atoms[role_map["CA"]])
        proxy["name"] = "BB"
        proxy["x"] = float(com[0])
        proxy["y"] = float(com[1])
        proxy["z"] = float(com[2])
        proxy["segid"] = "PROA"
        proxy["element"] = "C"
        bb_atom_index = len(protein_atoms)
        protein_atoms.append(proxy)

        bb_entries.append(
            {
                "bb_residue_index": seq_idx,
                "bb_resseq": int(key[1]),
                "bb_chain": str(key[0]),
                "bb_icode": str(key[2]),
                "bb_atom_index": bb_atom_index,
                "bb_type": bb_type,
                "bb_secondary_structure": bb_secondary,
                "bb_type_source": bb_type_source,
                "atom_indices": idxs,
                "atom_mask": [1, 1, 1, 1],
                "weights": weights,
                "reference_atom_indices": idxs,
                "reference_atom_coords": coords_row,
                "bb_comment": (
                    f"Backbone virtual BB map resid={seq_idx} resseq={key[1]} chain={key[0]} "
                    f"type={bb_type} ss={bb_secondary} source={bb_type_source} "
                    f"bb_atom_index={bb_atom_index} atom_indices={idxs}"
                ),
            }
        )
    if not bb_entries:
        raise ValueError("No complete backbone mapping entries could be generated.")
    return protein_atoms, bb_entries


def derive_chain_break_metadata(bb_entries):
    if not bb_entries:
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

    chain_first_residue = []
    chain_counts = [1]
    prev_chain = str(bb_entries[0].get("bb_chain", "")).strip() or " "
    prev_resseq = int(bb_entries[0].get("bb_resseq", 0))
    for residue_index, entry in enumerate(bb_entries[1:], start=1):
        chain = str(entry.get("bb_chain", "")).strip() or " "
        resseq = int(entry.get("bb_resseq", residue_index))
        chain_break = chain != prev_chain
        seq_break = (not chain_break) and resseq != prev_resseq + 1
        if chain_break or seq_break:
            chain_first_residue.append(residue_index)
            if chain_break:
                chain_counts.append(1)
            else:
                chain_counts[-1] += 1
        prev_chain = chain
        prev_resseq = resseq
    return np.asarray(chain_first_residue, dtype=np.int32), np.asarray(chain_counts, dtype=np.int32)


def write_hybrid_mapping_h5(
    path: Path,
    bb_entries,
    total_martini_atoms,
    env_atom_indices,
    n_protein_atoms,
):
    with h5py.File(path, "w") as h5:
        inp = h5.create_group("input")

        ctrl = inp.create_group("hybrid_control")
        ctrl.attrs["activation_stage"] = b"minimization"
        ctrl.attrs["preprod_protein_mode"] = b"rigid_body"
        ctrl.attrs["preprod_lipid_headgroup_roles"] = b"PO4"
        ctrl.attrs["exclude_intra_protein_martini"] = np.int8(1)
        ctrl.attrs["virtual_backbone_com_mode"] = np.int8(1)

        bb_grp = inp.create_group("hybrid_bb_map")
        bb_grp.attrs["atom_index_space"] = b"protein_aa_pdb_0based"
        bb_grp.attrs["runtime_index_space"] = b"stage_runtime_after_injection"
        bb_grp.create_dataset(
            "bb_residue_index",
            data=np.array([b.get("bb_residue_index", b["bb_resseq"]) for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_resseq",
            data=np.array([b["bb_resseq"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_chain_id",
            data=np.array(
                [str(b.get("bb_chain", "")).strip() or " " for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_icode",
            data=np.array(
                [str(b.get("bb_icode", "")).strip() or " " for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_atom_index",
            data=np.array([b["bb_atom_index"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_type",
            data=np.array(
                [str(b.get("bb_type", "P5")).strip() for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_secondary_structure",
            data=np.array(
                [str(b.get("bb_secondary_structure", "C")).strip() or "C" for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "bb_type_source",
            data=np.array(
                [str(b.get("bb_type_source", "structure_default")).strip() for b in bb_entries],
                dtype=h5py.string_dtype(encoding="utf-8"),
            ),
        )
        bb_grp.create_dataset(
            "atom_indices",
            data=np.array([b["atom_indices"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "atom_mask",
            data=np.array([b["atom_mask"] for b in bb_entries], dtype=np.int8),
        )
        bb_grp.create_dataset(
            "weights",
            data=np.array([b["weights"] for b in bb_entries], dtype=np.float32),
        )
        bb_grp.create_dataset(
            "protein_id",
            data=np.zeros(len(bb_entries), dtype=np.int32),
        )
        bb_grp.attrs["virtual_backbone_com_mode"] = np.int8(1)
        # Optional reference backbone metadata for auditability in .up/.h5 artifacts.
        bb_grp.create_dataset(
            "reference_atom_names",
            data=np.array(["N", "CA", "C", "O"], dtype=h5py.string_dtype(encoding="utf-8")),
        )
        bb_grp.create_dataset(
            "reference_atom_indices",
            data=np.array([b["reference_atom_indices"] for b in bb_entries], dtype=np.int32).reshape(
                len(bb_entries), 4
            ),
        )
        bb_grp.create_dataset(
            "reference_atom_coords",
            data=np.array([b["reference_atom_coords"] for b in bb_entries], dtype=np.float32).reshape(
                len(bb_entries), 4, 3
            ),
        )
        bb_grp.create_dataset(
            "bb_comment",
            data=np.array([b["bb_comment"] for b in bb_entries], dtype=h5py.string_dtype(encoding="utf-8")),
        )

        env_grp = inp.create_group("hybrid_env_topology")
        env_grp.create_dataset(
            "env_atom_indices",
            data=np.array(env_atom_indices, dtype=np.int32),
        )
        membership = np.full(total_martini_atoms, -1, dtype=np.int32)
        membership[:n_protein_atoms] = 0
        env_grp.create_dataset(
            "protein_membership",
            data=membership,
        )

        chain_first_residue, chain_counts = derive_chain_break_metadata(bb_entries)
        if chain_first_residue.size:
            break_grp = inp.create_group("chain_break")
            break_grp.create_dataset("chain_first_residue", data=chain_first_residue, dtype=np.int32)
            break_grp.create_dataset("chain_counts", data=chain_counts, dtype=np.int32)


def runtime_input_pdb_path(script_dir, pdb_id):
    """Resolve runtime protein+bilayer PDB path with optional override."""
    override = os.environ.get("UPSIDE_RUNTIME_PDB_FILE", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(script_dir, f"pdb/{pdb_id}.pdb")


def read_martini3_nonbond_params(itp_file):
    _, pair_params = parse_dry_forcefield(itp_file)
    martini_table = {
        key: (value["sigma_nm"], value["epsilon_kj_mol"])
        for key, value in pair_params.items()
    }
    print(f"Read {len(martini_table)//2} unique nonbonded parameter pairs")
    return martini_table


def read_martini_masses(ff_file):
    ff_file_path = Path(ff_file).expanduser()
    if not ff_file_path.is_absolute():
        ff_file_path = (WORKFLOW_DIR / ff_file_path).resolve()
    else:
        ff_file_path = ff_file_path.resolve()
    
    if not ff_file_path.exists():
        raise ValueError(
            f"Force-field file '{ff_file}' not found; expected {ff_file_path}"
        )
    
    return parse_itp_atomtype_masses(ff_file_path)


def build_hybrid_mass_array_up(atom_types, particle_class, martini_masses):
    mass = np.empty(len(atom_types), dtype=np.float32)
    for index, (atom_type, class_value) in enumerate(
        zip(atom_types, particle_class)
    ):
        class_name = (
            class_value.decode("ascii")
            if isinstance(class_value, (bytes, np.bytes_))
            else str(class_value)
        ).strip().upper()
        if class_name == "PROTEIN":
            mass[index] = 1.0
            continue
        if atom_type not in martini_masses:
            raise ValueError(
                f"Mass not found for atom type '{atom_type}' (atom index {index}); "
                f"available atom types: {sorted(martini_masses.keys())}"
            )
        mass[index] = martini_masses[atom_type] / 12.0
    return mass


def _load_martini_topology(ff_path, ff_files, stage_lipidhead_fc):
    """Read the dry MARTINI force field, DOPC/ion/water topologies and masses."""
    def pick_ff_file(name, required=True):
        path = os.path.join(ff_path, name)
        if os.path.exists(path):
            return path
        if required:
            raise ValueError(
                f"Required dry MARTINI file '{name}' not found in '{ff_path}'. Available: {ff_files}"
            )
        return None

    martini_param_file = pick_ff_file("dry_martini_v2.1.itp")
    martini_table = read_martini3_nonbond_params(martini_param_file)

    if not martini_table:
        raise ValueError(f"Could not read MARTINI parameters from '{martini_param_file}'")

    print("Protein runtime representation: N/CA/C carriers; regenerated O/BB virtual sites")

    dopc_param_file = pick_ff_file("dry_martini_v2.1_lipids.itp")
    lipid_preproc_defs = {}
    if stage_lipidhead_fc > 0.0:
        lipid_preproc_defs['BILAYER_LIPIDHEAD_FC'] = stage_lipidhead_fc
    full_topology = parse_itp_file(dopc_param_file, preprocessor_defines=lipid_preproc_defs)

    dopc_molecule = None
    for mol_name in full_topology['molecules'].keys():
        if 'DOPC' in mol_name.upper() or 'DOP' in mol_name.upper():
            dopc_molecule = mol_name
            break

    if dopc_molecule:
        dopc_topology = parse_itp_file(
            dopc_param_file, dopc_molecule, preprocessor_defines=lipid_preproc_defs
        )
    else:
        available_molecules = list(full_topology['molecules'].keys())
        raise ValueError(
            f"DOPC molecule not found in '{dopc_param_file}'. Available molecules: {available_molecules}"
        )

    dopc_bead_types = [atom['type'] for atom in dopc_topology['atoms']]
    dopc_charges = [atom['charge'] for atom in dopc_topology['atoms']]
    dopc_atom_to_type = {atom['atom']: atom['type'] for atom in dopc_topology['atoms']}
    dopc_atom_to_charge = {atom['atom']: atom['charge'] for atom in dopc_topology['atoms']}
    print(f"Read DOPC topology: {len(dopc_bead_types)} bead types from {dopc_param_file}")

    ion_param_file = pick_ff_file("dry_martini_v2.1_ions.itp", required=False)
    if ion_param_file:
        ion_topology = parse_itp_file(ion_param_file)
        na_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'NA']
        cl_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'CL']

        na_bead_types = [na_atoms[0]['type']] if na_atoms else []
        na_charges = [na_atoms[0]['charge']] if na_atoms else []
        cl_bead_types = [cl_atoms[0]['type']] if cl_atoms else []
        cl_charges = [cl_atoms[0]['charge']] if cl_atoms else []
        print(f"Ion topology loaded: NA={len(na_bead_types)} type(s), CL={len(cl_bead_types)} type(s)")
    else:
        na_bead_types, na_charges = [], []
        cl_bead_types, cl_charges = [], []
        print("Ion topology file not found in selected FF")

    water_param_file = pick_ff_file("dry_martini_v2.1_solvents.itp", required=False)
    if water_param_file:
        water_topology = parse_itp_file(water_param_file)
        water_atoms = [atom for atom in water_topology['atoms'] if atom['residue'].upper() == 'W']
        water_bead_types = [atom['type'] for atom in water_atoms]
        water_charges = [atom['charge'] for atom in water_atoms]
        print(f"Read water topology: {len(water_bead_types)} bead types from {water_param_file}")
    else:
        water_bead_types, water_charges = [], []
        print("Water topology file not found in selected FF")

    martini_masses = read_martini_masses(martini_param_file)
    print(f"Read {len(martini_masses)} atom type masses from force field file")

    dopc_bonds = [(bond['i'], bond['j']) for bond in dopc_topology['bonds']]
    dopc_bond_lengths = [bond['r0'] for bond in dopc_topology['bonds']]
    dopc_bond_force_constants = [bond['k'] for bond in dopc_topology['bonds']]

    dopc_angles = [(angle['i'], angle['j'], angle['k']) for angle in dopc_topology['angles']]
    dopc_angle_equil_deg = [angle['theta0'] for angle in dopc_topology['angles']]
    dopc_angle_force_constants = [angle['force_k'] for angle in dopc_topology['angles']]
    dopc_position_restraints = dopc_topology.get('position_restraints', [])

    print(f"Read DOPC connectivity: {len(dopc_bonds)} bonds, {len(dopc_angles)} angles")

    if not dopc_bonds:
        raise ValueError(f"No DOPC bonds found in topology from '{dopc_param_file}'")

    if not dopc_angles:
        raise ValueError(f"No DOPC angles found in topology from '{dopc_param_file}'")

    return {
        'martini_table': martini_table,
        'martini_masses': martini_masses,
        'dopc_bead_types': dopc_bead_types,
        'dopc_charges': dopc_charges,
        'dopc_atom_to_type': dopc_atom_to_type,
        'dopc_atom_to_charge': dopc_atom_to_charge,
        'dopc_bonds': dopc_bonds,
        'dopc_bond_lengths': dopc_bond_lengths,
        'dopc_bond_force_constants': dopc_bond_force_constants,
        'dopc_angles': dopc_angles,
        'dopc_angle_equil_deg': dopc_angle_equil_deg,
        'dopc_angle_force_constants': dopc_angle_force_constants,
        'dopc_position_restraints': dopc_position_restraints,
        'na_bead_types': na_bead_types,
        'na_charges': na_charges,
        'cl_bead_types': cl_bead_types,
        'cl_charges': cl_charges,
        'water_bead_types': water_bead_types,
        'water_charges': water_charges,
    }


def _parse_runtime_pdb_atoms(input_pdb_file, topo):
    """Parse the runtime PDB, map each atom to its MARTINI bead type, read the box."""
    dopc_atom_to_type = topo['dopc_atom_to_type']
    dopc_atom_to_charge = topo['dopc_atom_to_charge']
    water_bead_types = topo['water_bead_types']
    water_charges = topo['water_charges']
    na_bead_types = topo['na_bead_types']
    na_charges = topo['na_charges']
    cl_bead_types = topo['cl_bead_types']
    cl_charges = topo['cl_charges']
    dopc_bead_types = topo['dopc_bead_types']

    print("Reading PDB structure")
    initial_positions = []
    atom_types = []
    charges = []
    residue_ids = []
    atom_names = []
    residue_names = []
    chain_ids = []
    seg_ids = []

    protein_atoms_for_mapping, _ = parse_pdb(Path(input_pdb_file))
    protein_like_atoms = [
        atom for atom in protein_atoms_for_mapping
        if atom["resname"].upper() in {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
        }
    ]
    if protein_like_atoms:
        protein_backbone_atoms_for_mapping, _ = extract_protein_backbone_atoms_from_aa(
            protein_atoms_for_mapping
        )
        bb_type_by_residue = map_backbone_types_from_martinize_fallback(
            protein_backbone_atoms_for_mapping
        )
    else:
        bb_type_by_residue = {}

    protein_residue_names = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
        'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX'
    }
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            atom_name_raw = line[12:16].strip()
            atom_name = atom_name_raw.upper()
            atom_names.append(atom_name)
            residue_id = int(line[22:26])
            residue_ids.append(residue_id)
            # Read 4-character residue field to support DOPC while preserving
            # standard 3-letter protein residues (e.g., " ASN" -> "ASN").
            residue_name = line[17:21].strip().upper()
            residue_names.append(residue_name)
            chain_id = line[21:22].strip()
            chain_ids.append(chain_id)
            seg_id = line[72:76].strip()
            seg_ids.append(seg_id)

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])

            is_protein = (residue_name in protein_residue_names)

            # Map to MARTINI type based on context
            if is_protein:
                if atom_name not in {"N", "CA", "C", "O", "BB"}:
                    raise ValueError(
                        f"Protein atom '{atom_name}' in residue '{residue_name}' is not one of N/CA/C/O/BB"
                    )
                icode = line[26:27].strip()
                res_key = (chain_id, residue_id, icode, residue_name)
                martini_type = bb_type_by_residue.get(res_key, "P5")
                charge = float(BB_TYPE_CHARGE.get(martini_type, 0.0))
            elif residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file (for both lipid and mixed systems)
                if atom_name in dopc_atom_to_type:
                    martini_type = dopc_atom_to_type[atom_name]
                    charge = dopc_atom_to_charge[atom_name]
                else:
                    available_atom_names = sorted(dopc_atom_to_type.keys())
                    raise ValueError(
                        f"Unknown DOPC atom '{atom_name}' in residue '{residue_name}'. "
                        f"Available DOPC atom names: {available_atom_names}"
                    )
            elif residue_name == 'W':
                if not water_bead_types:
                    raise ValueError("No water bead types found in topology")
                martini_type = water_bead_types[0]
                charge = water_charges[0]
            elif residue_name == 'NA':
                if not na_bead_types:
                    raise ValueError("No sodium bead types found in topology")
                martini_type = na_bead_types[0]
                charge = na_charges[0]
            elif residue_name == 'CL':
                if not cl_bead_types:
                    raise ValueError("No chloride bead types found in topology")
                martini_type = cl_bead_types[0]
                charge = cl_charges[0]
            else:
                raise ValueError(
                    f"Unknown residue type '{residue_name}' for atom '{atom_name}'. "
                    "Supported residue types: PROTEIN, DOPC, W, NA, CL"
                )

            atom_types.append(martini_type)
            charges.append(charge)

    initial_positions = np.array(initial_positions, dtype=float)
    atom_types = np.array(atom_types)
    charges = np.array(charges, dtype=float)
    residue_ids = np.array(residue_ids, dtype=int)
    atom_names = np.array(atom_names)
    residue_names = np.array(residue_names)
    chain_ids = np.array(chain_ids)
    seg_ids = np.array(seg_ids)
    n_atoms = len(initial_positions)

    n_lipid_atoms = int(np.sum(np.isin(residue_names, ["DOPC", "DOP"])))
    n_lipid_mols = n_lipid_atoms // len(dopc_bead_types) if dopc_bead_types else 0
    print(f"Full-resolution mode: keeping all {n_lipid_mols} DOPC lipids "
          f"({n_lipid_atoms} beads) as individual particles")

    print(f"Reading box dimensions from {input_pdb_file}...")
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                fields = line.split()
                if len(fields) >= 4:
                    x_len = float(fields[1])
                    y_len = float(fields[2])
                    z_len = float(fields[3])
                    print(f"Found CRYST1 record: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
                    break
        else:
            raise ValueError(f"No CRYST1 record found in PDB file '{input_pdb_file}'")

    print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
    print(f"Box volume: {x_len * y_len * z_len:.1f} A^3")
    print(f"Total atoms: {n_atoms}")

    return (initial_positions, atom_types, charges, residue_ids, atom_names,
            residue_names, chain_ids, seg_ids, n_atoms, x_len, y_len, z_len)


def _group_atoms_into_molecules(residue_ids, residue_names, atom_names, chain_ids, seg_ids, n_atoms):
    """Group atoms into molecules by protein chain or by residue for everything else."""
    molecules = []
    current_mol_atoms = []
    current_mol_indices = []
    current_resid = None
    current_resname = None
    current_chain = None

    protein_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
    }

    def pick_chain_id(idx):
        # Prefer PDB chain ID to preserve distinct protein chains in multi-chain
        # systems where segid may be shared (e.g., "PROA" for all protein atoms).
        if chain_ids[idx]:
            return chain_ids[idx]
        return seg_ids[idx]

    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        if resname in protein_residues:
            mol_type = 'PROTEIN'
        else:
            # Normalize DOP (resname in PDB) to DOPC for reporting and selection
            mol_type = 'DOPC' if resname == 'DOP' else resname

        # Start new molecule if protein chain changes or residue ID/name changes for non-proteins
        chain_id = pick_chain_id(i) if mol_type == 'PROTEIN' else ""
        start_new = False
        if mol_type == 'PROTEIN':
            if chain_id != current_chain:
                start_new = True
        else:
            if resid != current_resid or resname != current_resname:
                start_new = True

        if start_new:
            if current_mol_atoms:
                molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))
            current_mol_atoms = [atom_name]
            current_mol_indices = [i]
            current_resid = resid
            current_resname = resname
            current_mol_type = mol_type
            current_chain = chain_id
        else:
            current_mol_atoms.append(atom_name)
            current_mol_indices.append(i)

    if current_mol_atoms:
        molecules.append((current_mol_type, current_mol_atoms, current_mol_indices))

    # Build per-atom molecule indices
    molecule_ids = np.zeros(n_atoms, dtype=np.int32)
    for mol_idx, (_, _, atom_indices) in enumerate(molecules):
        for atom_idx in atom_indices:
            molecule_ids[atom_idx] = mol_idx
    # Validate that each non-protein molecule contains only atoms of the same residue type
    for i, (mol_type, atoms, indices) in enumerate(molecules):
        residue_names_in_mol = [residue_names[idx] for idx in indices]
        unique_resnames = set(residue_names_in_mol)
        if mol_type != 'PROTEIN':
            if len(unique_resnames) > 1:
                raise ValueError(
                    f"Molecule {i} contains mixed residue types: {unique_resnames}; "
                    f"atoms={atoms}, residue_names={residue_names_in_mol}"
                )

    # Count molecules by type, but group all protein residues together
    mol_counts = Counter()
    protein_sequence = []
    protein_sequence_seen = set()
    for idx, (resid, resname) in enumerate(zip(residue_ids, residue_names)):
        if resname not in protein_residues:
            continue
        key = (pick_chain_id(idx), int(resid))
        if key in protein_sequence_seen:
            continue
        protein_sequence_seen.add(key)
        protein_sequence.append(normalize_resname(str(resname)))

    protein_residue_keys = set()
    for resid, resname, idx in zip(residue_ids, residue_names, range(n_atoms)):
        if resname in protein_residues:
            protein_residue_keys.add((pick_chain_id(idx), resid))

    for mol_type, _, _ in molecules:
        mol_counts[mol_type] += 1

    protein_residue_count = len(protein_residue_keys)
    if protein_residue_count > 0:
        protein_chains = mol_counts.get('PROTEIN', 0)
        mol_counts['PROTEIN'] = f"{protein_chains} chain(s) ({protein_residue_count} residues)"

    dopc_count = mol_counts.get('DOPC', 0)
    water_count = mol_counts.get('W', 0)

    print("\nMolecule summary")
    for moltype, count in mol_counts.items():
        print(f"{moltype}: {count} molecules")

    return molecules, molecule_ids, protein_sequence, dopc_count, water_count


def _build_bonded_terms(molecules, topo, initial_positions,
                        bond_conversion, angle_conversion, dihedral_conversion,
                        protein_bonds, protein_angles, protein_dihedrals,
                        protein_constraints, protein_position_restraints):
    """Build DOPC and protein bonds/angles/dihedrals plus lipid position restraints."""
    dopc_bead_types = topo['dopc_bead_types']
    dopc_bonds = topo['dopc_bonds']
    dopc_bond_lengths = topo['dopc_bond_lengths']
    dopc_bond_force_constants = topo['dopc_bond_force_constants']
    dopc_angles = topo['dopc_angles']
    dopc_angle_equil_deg = topo['dopc_angle_equil_deg']
    dopc_angle_force_constants = topo['dopc_angle_force_constants']
    dopc_position_restraints = topo['dopc_position_restraints']

    print("\nCreating connectivity")

    bonds_list = []
    bond_lengths_list = []
    bond_force_constants_list = []
    angles_list = []
    angle_equil_deg_list = []
    angle_force_constants_list = []
    dihedrals_list = []
    dihedral_equil_deg_list = []
    dihedral_force_constants_list = []
    dihedral_type_list = []
    lipid_restraint_indices = []
    lipid_restraint_ref_pos = []
    lipid_restraint_spring_xyz = []

    dopc_molecules = [mol for mol in molecules if mol[0] == 'DOPC']
    for mol_idx, (_, atom_names_mol, atom_indices) in enumerate(dopc_molecules):
        name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}

        for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
            if bond_idx1 < len(dopc_bead_types) and bond_idx2 < len(dopc_bead_types):
                atom1_name = atom_names_mol[bond_idx1]
                atom2_name = atom_names_mol[bond_idx2]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                bonds_list.append([atom1, atom2])
                bond_lengths_list.append(dopc_bond_lengths[i] * 10.0)
                bond_force_constants_list.append(
                    dopc_bond_force_constants[i] * bond_conversion
                )

        for i, (angle_idx1, angle_idx2, angle_idx3) in enumerate(dopc_angles):
            if (angle_idx1 < len(dopc_bead_types) and angle_idx2 < len(dopc_bead_types) and
                angle_idx3 < len(dopc_bead_types)):
                atom1_name = atom_names_mol[angle_idx1]
                atom2_name = atom_names_mol[angle_idx2]
                atom3_name = atom_names_mol[angle_idx3]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                atom3 = name_to_idx[atom3_name]
                angles_list.append([atom1, atom2, atom3])
                angle_equil_deg_list.append(dopc_angle_equil_deg[i])
                angle_force_constants_list.append(
                    dopc_angle_force_constants[i] * angle_conversion
                )

        for restraint in dopc_position_restraints:
            local_idx = restraint['i']
            if 0 <= local_idx < len(atom_indices):
                atom_idx = atom_indices[local_idx]
                lipid_restraint_indices.append(atom_idx)
                lipid_restraint_ref_pos.append(initial_positions[atom_idx].tolist())
                lipid_restraint_spring_xyz.append([
                    restraint['fx'] * bond_conversion,
                    restraint['fy'] * bond_conversion,
                    restraint['fz'] * bond_conversion,
                ])

    n_lipid_bonds = len(bonds_list)
    n_lipid_angles = len(angles_list)
    print(f"Full-resolution mode: created {n_lipid_bonds} DOPC bonds, "
          f"{n_lipid_angles} angles, {len(lipid_restraint_indices)} restraints "
          f"across {len(dopc_molecules)} lipids")

    protein_bond_count = 0
    protein_angle_count = 0
    protein_dihedral_count = 0
    protein_constraint_count = 0

    if protein_bonds or protein_constraints:
        print("\nProtein connectivity")
        print(f"Found {len(protein_bonds)} bonds, {len(protein_angles)} angles, {len(protein_dihedrals)} dihedrals")
        print(f"Found {len(protein_constraints)} constraints, {len(protein_position_restraints)} position restraints")

        for i, j, r0_nm, k_kj in protein_bonds:
            r0_angstrom = r0_nm * 10.0
            k_upside = k_kj * bond_conversion

            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_bond_count += 1

        for i, j, r0_nm, k_kj in protein_constraints:
            r0_angstrom = r0_nm * 10.0
            k_upside = k_kj * bond_conversion

            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_constraint_count += 1

        for i, j, k, theta0_deg, k_kj in protein_angles:
            theta0_upside = theta0_deg
            k_upside = k_kj * angle_conversion

            angles_list.append([i, j, k])
            angle_equil_deg_list.append(theta0_upside)
            angle_force_constants_list.append(k_upside)
            protein_angle_count += 1

        for i, j, k, l, phi0_deg, k_kj, func_type in protein_dihedrals:
            phi0_upside = phi0_deg
            k_upside = k_kj * dihedral_conversion  # kJ/mol to E_up

            dihedrals_list.append([i, j, k, l])
            dihedral_equil_deg_list.append(phi0_upside)
            dihedral_force_constants_list.append(k_upside)
            dihedral_type_list.append(func_type)
            protein_dihedral_count += 1

        print(f"Added {protein_bond_count} protein bonds")
        print(f"Added {protein_constraint_count} protein constraints (as bonds with large spring constants)")
        print(f"Added {protein_angle_count} protein angles")
        print(f"Added {protein_dihedral_count} protein dihedrals")
    else:
        print("\nNo explicit protein MARTINI connectivity is used for AA backbone runtime mode")

    print(f"Total system bonds: {len(bonds_list)}")
    print(f"Total system angles: {len(angles_list)}")
    print(f"Total system dihedrals: {len(dihedrals_list)}")

    return (bonds_list, bond_lengths_list, bond_force_constants_list,
            angles_list, angle_equil_deg_list, angle_force_constants_list,
            dihedrals_list, dihedral_equil_deg_list, dihedral_force_constants_list,
            dihedral_type_list, lipid_restraint_indices, lipid_restraint_ref_pos,
            lipid_restraint_spring_xyz, protein_bond_count, protein_constraint_count,
            protein_angle_count, protein_dihedral_count)


def _write_nonbonded_pairs(t, martini_potential, n_atoms, atom_types, charges,
                           bonds_list, protein_exclusions, martini_table,
                           energy_conversion, length_conversion):
    """Write the optimized non-bonded pair table (pairs + coefficient indices)."""
    # Create atom indices and charges arrays for the potential
    t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
    t.create_array(martini_potential, 'charges', obj=charges)

    # Create pairs and optimized coefficient indices for non-bonded interactions.
    # The full coefficient table is O(N^2) and can OOM for large membrane boxes;
    # store unique coefficient rows plus one int index per pair instead.
    # Create sets for exclusions (MARTINI uses nrexcl=1, so only 1-2 exclusions)
    bonded_pairs_12 = set()  # Directly bonded (1-2) - full exclusion
    additional_exclusions = set()  # Additional exclusions from ITP file

    # Add 1-2 exclusions from bond list
    for bond in bonds_list:
        sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
        bonded_pairs_12.add(sorted_bond)

    # Add additional exclusions from protein ITP file
    if protein_exclusions:
        for exclusion in protein_exclusions:
            sorted_exclusion = (min(exclusion[0], exclusion[1]), max(exclusion[0], exclusion[1]))
            additional_exclusions.add(sorted_exclusion)

    # Generate all unique pairs (i < j) with proper exclusions (nrexcl=1)
    excluded_12_count = 0
    excluded_additional_count = 0
    total_candidate_pairs = n_atoms * (n_atoms - 1) // 2
    total_pairs_written = 0
    unique_coeffs = []
    coeff_to_index = {}
    atom_type_strings = [
        atom_type.decode('utf-8') if isinstance(atom_type, bytes) else str(atom_type)
        for atom_type in atom_types
    ]
    atom_classes = [(atom_type_strings[i], float(charges[i])) for i in range(n_atoms)]
    class_to_id = {}
    atom_class_ids = np.empty(n_atoms, dtype=np.int32)
    class_values = []
    for i, cls in enumerate(atom_classes):
        class_id = class_to_id.get(cls)
        if class_id is None:
            class_id = len(class_values)
            class_to_id[cls] = class_id
            class_values.append(cls)
        atom_class_ids[i] = class_id

    def coefficient_index_for_classes(class_i, class_j, atom_i):
        type1, q1_raw = class_values[class_i]
        type2, q2_raw = class_values[class_j]
        if (type1, type2) in martini_table:
            sigma_nm, epsilon_kj = martini_table[(type1, type2)]
        elif (type2, type1) in martini_table:
            sigma_nm, epsilon_kj = martini_table[(type2, type1)]
        else:
            available_types = sorted(set([t[0] for t in martini_table.keys()] + [t[1] for t in martini_table.keys()]))
            raise ValueError(
                f"Missing interaction parameters for bead type pair ({type1}, {type2}) "
                f"at atom index {atom_i}; available bead types: {available_types}"
            )

        epsilon = np.float32(epsilon_kj / energy_conversion)
        sigma = np.float32(sigma_nm * length_conversion)
        q1 = np.float32(q1_raw)
        q2 = np.float32(q2_raw)
        coeff_key = (float(epsilon), float(sigma), float(q1), float(q2))
        coeff_index = coeff_to_index.get(coeff_key)
        if coeff_index is None:
            coeff_index = len(unique_coeffs)
            coeff_to_index[coeff_key] = coeff_index
            unique_coeffs.append((float(epsilon), float(sigma), float(q1), float(q2)))
        return coeff_index

    pairs_array = t.create_earray(
        martini_potential,
        'pairs',
        atom=tb.Int32Atom(),
        shape=(0, 2),
        expectedrows=total_candidate_pairs,
    )
    coeff_index_array = t.create_earray(
        martini_potential,
        'coefficient_indices',
        atom=tb.Int64Atom(),
        shape=(0,),
        expectedrows=total_candidate_pairs,
    )

    bonded_by_i = defaultdict(set)
    for i, j in bonded_pairs_12:
        bonded_by_i[i].add(j)
    additional_by_i = defaultdict(set)
    for i, j in additional_exclusions:
        additional_by_i[i].add(j)

    for i in range(n_atoms):
        if i + 1 >= n_atoms:
            continue
        js = np.arange(i + 1, n_atoms, dtype=np.int32)
        bonded_js = bonded_by_i.get(i)
        additional_js = additional_by_i.get(i)
        if bonded_js:
            bonded_mask = np.isin(js, np.fromiter(bonded_js, dtype=np.int32), assume_unique=True)
            excluded_12_count += int(np.count_nonzero(bonded_mask))
            js = js[~bonded_mask]
        if additional_js and js.size:
            additional_candidates = np.fromiter(
                (j for j in additional_js if not bonded_js or j not in bonded_js),
                dtype=np.int32,
            )
            if additional_candidates.size:
                additional_mask = np.isin(js, additional_candidates, assume_unique=True)
                excluded_additional_count += int(np.count_nonzero(additional_mask))
                js = js[~additional_mask]
        if js.size == 0:
            continue

        pairs_chunk = np.empty((js.size, 2), dtype=np.int32)
        pairs_chunk[:, 0] = i
        pairs_chunk[:, 1] = js

        coeff_indices = np.empty(js.size, dtype=np.int64)
        class_i = int(atom_class_ids[i])
        js_class_ids = atom_class_ids[js]
        for class_j in np.unique(js_class_ids):
            class_mask = js_class_ids == class_j
            count = int(np.count_nonzero(class_mask))
            coeff_index = coefficient_index_for_classes(class_i, int(class_j), i)
            coeff_indices[class_mask] = coeff_index

        pairs_array.append(pairs_chunk)
        coeff_index_array.append(coeff_indices)
        total_pairs_written += int(js.size)

    print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions (nrexcl=1)")
    print(f"Excluded {excluded_additional_count} additional pairs from ITP exclusions")
    print(f"Wrote {total_pairs_written} non-bonded pairs with {len(unique_coeffs)} unique coefficient rows")
    martini_potential._v_attrs.optimized_format = 1
    coeff_array = np.array(unique_coeffs, dtype='f4')
    if coeff_array.size == 0:
        coeff_array = np.zeros((0, 4), dtype='f4')
    t.create_array(martini_potential, 'coefficients', obj=coeff_array)

    return excluded_12_count, excluded_additional_count


def convert_stage(pdb_id=None, stage='minimization', run_dir=None):
    """Build the HDF5 input for one MARTINI workflow stage."""
    if 'UPSIDE_HOME' not in os.environ:
        raise ValueError("UPSIDE_HOME is required")

    if pdb_id is None:
        if len(sys.argv) > 1:
            pdb_id = sys.argv[1]
        else:
            raise ValueError(f"No PDB ID provided; usage: python {sys.argv[0]} <pdb_id>")

    if run_dir is None:
        if len(sys.argv) > 2 and sys.argv[2] != '--stage':
            run_dir = sys.argv[2]
        else:
            run_dir = "outputs/martini_test"
    os.makedirs(run_dir, exist_ok=True)

    stage = os.environ.get('UPSIDE_SIMULATION_STAGE', stage)

    if stage not in STAGE_PARAMS:
        valid = ", ".join(sorted(STAGE_PARAMS))
        raise ValueError(f"Unknown MARTINI preparation stage {stage!r}; expected one of: {valid}")
    params = STAGE_PARAMS[stage]
    stage_lipidhead_fc = float(os.environ.get('UPSIDE_BILAYER_LIPIDHEAD_FC', '0'))

    print(f"Preparing MARTINI stage {stage} for {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    workflow_dir = str(WORKFLOW_DIR)
    print("Reading dry MARTINI parameters")
    ff_dir = Path(
        os.environ.get('UPSIDE_MARTINI_FF_DIR', str(REPO_ROOT / "parameters" / "dryMARTINI"))
    ).expanduser()
    if not ff_dir.is_absolute():
        ff_dir = (REPO_ROOT / ff_dir).resolve()
    else:
        ff_dir = ff_dir.resolve()
    ff_path = str(ff_dir)

    if not os.path.isdir(ff_path):
        raise ValueError(
            f"Force-field directory '{ff_path}' not found; set UPSIDE_MARTINI_FF_DIR"
        )

    ff_files = sorted(os.listdir(ff_path))
    print(f"Using force-field directory: {ff_path}")

    topo = _load_martini_topology(ff_path, ff_files, stage_lipidhead_fc)
    martini_masses = topo['martini_masses']
    
    energy_conversion_raw = os.environ.get('UPSIDE_MARTINI_ENERGY_CONVERSION', '').strip()
    length_conversion_raw = os.environ.get('UPSIDE_MARTINI_LENGTH_CONVERSION', '').strip()
    if not energy_conversion_raw:
        raise ValueError("Missing required environment variable UPSIDE_MARTINI_ENERGY_CONVERSION")
    if not length_conversion_raw:
        raise ValueError("Missing required environment variable UPSIDE_MARTINI_LENGTH_CONVERSION")
    energy_conversion = float(energy_conversion_raw)
    length_conversion = float(length_conversion_raw)
    coulomb_constant_native = float(os.environ.get('UPSIDE_MARTINI_COULOMB_CONSTANT_NATIVE', str(138.935458 / 15.0)))
    if energy_conversion <= 0.0:
        raise ValueError("UPSIDE_MARTINI_ENERGY_CONVERSION must be positive")
    if length_conversion <= 0.0:
        raise ValueError("UPSIDE_MARTINI_LENGTH_CONVERSION must be positive")

    pressure_conversion_bar_to_eup = 0.000020659
    bond_conversion = 1.0 / (energy_conversion * length_conversion ** 2)
    angle_conversion = 1.0 / energy_conversion
    dihedral_conversion = 1.0 / energy_conversion

    print(
        "Unit conversions: "
        f"energy={energy_conversion:g} kJ/mol per E_up, "
        f"length={length_conversion:g} A/nm, "
        f"pressure={pressure_conversion_bar_to_eup:.9f} E_up/A^3 per bar"
    )
    
    input_pdb_file = runtime_input_pdb_path(workflow_dir, pdb_id)
    print(f"Using MARTINI PDB as base structure: {input_pdb_file}")
    
    protein_bonds, protein_angles, protein_dihedrals, protein_constraints = [], [], [], []
    protein_exclusions = []
    protein_position_restraints = []

    (initial_positions, atom_types, charges, residue_ids, atom_names,
     residue_names, chain_ids, seg_ids, n_atoms, x_len, y_len, z_len) = \
        _parse_runtime_pdb_atoms(input_pdb_file, topo)
    
    molecules, molecule_ids, protein_sequence, dopc_count, water_count = \
        _group_atoms_into_molecules(
            residue_ids, residue_names, atom_names, chain_ids, seg_ids, n_atoms
        )
    
    (bonds_list, bond_lengths_list, bond_force_constants_list,
     angles_list, angle_equil_deg_list, angle_force_constants_list,
     dihedrals_list, dihedral_equil_deg_list, dihedral_force_constants_list,
     dihedral_type_list, lipid_restraint_indices, lipid_restraint_ref_pos,
     lipid_restraint_spring_xyz, protein_bond_count, protein_constraint_count,
     protein_angle_count, protein_dihedral_count) = _build_bonded_terms(
        molecules, topo, initial_positions,
        bond_conversion, angle_conversion, dihedral_conversion,
        protein_bonds, protein_angles, protein_dihedrals,
        protein_constraints, protein_position_restraints,
    )
    
    print("\nPreparing final structure")
    center_mask = np.ones(n_atoms, dtype=bool)
    center = np.mean(initial_positions[center_mask], axis=0)
    centered_positions = initial_positions - center
    half_box = np.array([x_len/2, y_len/2, z_len/2])
    centered_positions = (centered_positions + half_box) % (2*half_box) - half_box
    final_positions = centered_positions

    protein_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
    }
    particle_class = np.empty(n_atoms, dtype='S10')
    for i, resname in enumerate(residue_names):
        if resname in protein_residues:
            particle_class[i] = b"PROTEIN"
        elif resname == 'DOP' or resname == 'DOPC':
            particle_class[i] = b"LIPID"
        elif resname == 'W':
            particle_class[i] = b"WATER"
        elif resname in ('NA', 'CL'):
            particle_class[i] = b"ION"
        else:
            particle_class[i] = b"OTHER"

    print("\nCreating UPSIDE input file")
    input_file = f"{run_dir}/test.input.up"
    
    with tb.open_file(input_file, 'w') as t:
        # Create input group (required by UPSIDE)
        input_grp = t.create_group(t.root, 'input')
        
        # Create position array with correct format
        pos = np.zeros((n_atoms, 3, 1), dtype='f4')
        pos[:,:,0] = final_positions
        pos_array = t.create_array(input_grp, 'pos', obj=pos)
        pos_array._v_attrs.arguments = np.array([b'pos'])
        pos_array._v_attrs.shape = pos.shape
        pos_array._v_attrs.n_atoms = n_atoms
        pos_array._v_attrs.n_frames = 1
        pos_array._v_attrs.dim = 3
        pos_array._v_attrs.initialized = True
        
        # Create velocity array (required by UPSIDE)
        velocity = np.zeros((n_atoms, 3), dtype='f4')
        vel_array = t.create_array(input_grp, 'vel', obj=velocity)
        vel_array._v_attrs.arguments = np.array([b'vel'])
        vel_array._v_attrs.shape = velocity.shape
        vel_array._v_attrs.n_atoms = n_atoms
        vel_array._v_attrs.dim = 3
        vel_array._v_attrs.initialized = True
        
        # Create momentum array (required by UPSIDE)
        momentum = np.zeros((n_atoms, 3, 1), dtype='f4')
        mom_array = t.create_array(input_grp, 'mom', obj=momentum)
        mom_array._v_attrs.arguments = np.array([b'mom'])
        mom_array._v_attrs.shape = momentum.shape
        mom_array._v_attrs.n_atoms = n_atoms
        mom_array._v_attrs.dim = 3
        mom_array._v_attrs.initialized = True
        
        # Create mass array (required by UPSIDE)
        mass = build_hybrid_mass_array_up(
            atom_types,
            particle_class,
            martini_masses,
        )
        
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        mass_array._v_attrs.protein_mass_up = np.float32(1.0)
        mass_array._v_attrs.protein_mass_source = b"upside_core_unit_mass"
        mass_array._v_attrs.environment_mass_source = b"native_dry_martini_mass_div_12"

        print("Hybrid stage files use N/CA/C protein carriers with regenerated O/BB virtual sites")
        
        # Create stage-specific parameters group (always create this)
        stage_grp = t.create_group(input_grp, 'stage_parameters')
        stage_grp._v_attrs.current_stage = b'minimization'
        
        # Store minimization stage bond parameters (large spring constants)
        min_bonds_grp = t.create_group(stage_grp, 'minimization_bonds')
        min_bond_fc = np.array([1000000.0] * len(protein_bonds), dtype='f4')
        t.create_array(min_bonds_grp, 'force_constants', obj=min_bond_fc)
        
        # Store production stage bond parameters (regular spring constants)
        prod_bonds_grp = t.create_group(stage_grp, 'production_bonds')
        prod_bond_fc = np.array([bond[3] for bond in protein_bonds], dtype='f4')  # Regular k values
        t.create_array(prod_bonds_grp, 'force_constants', obj=prod_bond_fc)
        
        # Store minimization stage angle parameters (NORMANG angles)
        min_angles_grp = t.create_group(stage_grp, 'minimization_angles')
        min_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(min_angles_grp, 'force_constants', obj=min_angle_fc)
        
        prod_angles_grp = t.create_group(stage_grp, 'production_angles')
        prod_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(prod_angles_grp, 'force_constants', obj=prod_angle_fc)
        
        print(f"Stage-specific parameters: minimization bonds={len(min_bond_fc)}, production bonds={len(prod_bond_fc)}")
        print(f"Stage-specific parameters: minimization angles={len(min_angle_fc)}, production angles={len(prod_angle_fc)}")
        
        barostat_enable = int(os.environ.get('UPSIDE_NPT_ENABLE', '0'))
        if barostat_enable:
            print("\nCreating NPT barostat configuration")
            barostat_grp = t.create_group(input_grp, 'barostat')
            barostat_grp._v_attrs.target_p_xy = float(os.environ.get('UPSIDE_NPT_TARGET_PXY', '0.000020659'))
            barostat_grp._v_attrs.target_p_z = float(os.environ.get('UPSIDE_NPT_TARGET_PZ', '0.000020659'))
            barostat_grp._v_attrs.tau_p = float(os.environ.get('UPSIDE_NPT_TAU', '1.0'))
            base_compressibility = float(os.environ.get('UPSIDE_NPT_COMPRESSIBILITY', '14.521180763676'))
            barostat_grp._v_attrs.compressibility = base_compressibility
            barostat_grp._v_attrs.compressibility_xy = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_XY', str(base_compressibility))
            )
            barostat_grp._v_attrs.compressibility_z = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_Z', str(base_compressibility))
            )
            barostat_grp._v_attrs.interval = int(os.environ.get('UPSIDE_NPT_INTERVAL', '10'))
            barostat_grp._v_attrs.semi_isotropic = int(os.environ.get('UPSIDE_NPT_SEMI', '1'))
            barostat_grp._v_attrs.type = params['barostat_type']
            print(f"  Enabled: {barostat_enable}")
            print(f"  Type: {'Parrinello-Rahman' if barostat_grp._v_attrs.type == 1 else 'Berendsen'}")
            print(f"  Target Pxy: {barostat_grp._v_attrs.target_p_xy} E_up/Angstrom^3 (~1 bar)")
            print(f"  Target Pz: {barostat_grp._v_attrs.target_p_z} E_up/Angstrom^3 (~1 bar)")
            print(f"  Compressibility XY: {barostat_grp._v_attrs.compressibility_xy} Angstrom^3/E_up")
            print(f"  Compressibility Z: {barostat_grp._v_attrs.compressibility_z} Angstrom^3/E_up")
            print(f"  Tau_p: {barostat_grp._v_attrs.tau_p}")
            print(f"  Interval: {barostat_grp._v_attrs.interval} steps")
        else:
            print("\nNPT barostat disabled (NVT mode)")
        
        # Create type array
        type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
        type_array._v_attrs.arguments = np.array([b'type'])
        type_array._v_attrs.shape = atom_types.shape
        type_array._v_attrs.n_atoms = n_atoms
        type_array._v_attrs.initialized = True
        type_array._v_attrs.description = b"Interaction matrix type names (e.g., C1, C2, Qa, Qd)"

        particle_class_array = t.create_array(input_grp, 'particle_class', obj=particle_class)
        particle_class_array._v_attrs.arguments = np.array([b'particle_class'])
        particle_class_array._v_attrs.shape = particle_class.shape
        particle_class_array._v_attrs.n_atoms = n_atoms
        particle_class_array._v_attrs.initialized = True
        particle_class_array._v_attrs.description = b"Per-atom class: PROTEIN, LIPID, WATER, ION, OTHER"

        # Every physical degree of freedom uses one g-JF step. Regenerated O and BB sites are excluded:
        # their positions and BB-env derivatives are projected through N/CA/C. Production uses the requested
        # bare-particle mobility clock. Protein drag is the sum of the bead frictions inside the existing
        # protein--environment spline range.
        lipid_beads = np.where(particle_class == b"LIPID")[0].astype('i4')
        if lipid_beads.size:
            protein_carrier = (
                (particle_class == b"PROTEIN")
                & np.isin(atom_names, np.array(["N", "CA", "C"], dtype=atom_names.dtype))
            )
            dynamic_atom = (particle_class != b"PROTEIN") | protein_carrier
            brownian_atom_index = np.where(dynamic_atom)[0].astype('i4')
            numerical_time_step = float(os.environ.get("UPSIDE_MARTINI_TIME_STEP_UP", "0.009"))
            protein_time_ps = float(os.environ.get("UPSIDE_PROTEIN_TIME_PS_PER_STEP", "40.0"))
            martini_time_factor = float(os.environ.get("UPSIDE_MARTINI_TIME_FACTOR", "4.0"))
            target_diffusion = float(os.environ.get("UPSIDE_DOPC_TARGET_DIFFUSION_UM2_S", "11.5"))
            reference_temperature = float(os.environ.get("UPSIDE_DOPC_REFERENCE_TEMPERATURE_UP", "0.8647"))
            raw_relaxation_ps = float(os.environ.get("UPSIDE_DRY_MARTINI_RELAXATION_PS", "4.0"))
            dynamics_phase = os.environ.get("UPSIDE_MARTINI_DYNAMICS_PHASE", "production")
            if numerical_time_step <= 0.0 or protein_time_ps <= 0.0 or martini_time_factor <= 0.0:
                raise ValueError("MARTINI clock parameters must be positive")
            if target_diffusion <= 0.0 or reference_temperature <= 0.0 or raw_relaxation_ps <= 0.0:
                raise ValueError("DOPC transport target and reference temperature must be positive")
            if dynamics_phase not in {"early_equilibration", "overlap_settling", "production"}:
                raise ValueError(f"Unknown MARTINI dynamics phase '{dynamics_phase}'")
            raw_martini_time_ps = protein_time_ps / martini_time_factor
            raw_target_diffusion = target_diffusion * martini_time_factor
            time_unit_ps = 1.0e12 * math.sqrt(
                (0.012 * 1.0e-20) / (energy_conversion * 1000.0))
            if dynamics_phase == "production":
                raw_target_diffusion_A2_ps = raw_target_diffusion * 1.0e-4
                particle_diffusion_up = (
                    raw_target_diffusion_A2_ps * raw_martini_time_ps / numerical_time_step
                )
                particle_friction = reference_temperature / particle_diffusion_up
                friction = np.full(brownian_atom_index.size, particle_friction, dtype=np.float64)
            else:
                if dynamics_phase == "overlap_settling":
                    relaxation_up = numerical_time_step * raw_relaxation_ps / raw_martini_time_ps
                else:
                    relaxation_up = raw_relaxation_ps / time_unit_ps
                friction = mass[brownian_atom_index] / relaxation_up

            interface_cutoff = 12.0
            contact_count = np.zeros(n_atoms, dtype=np.int32)
            if dynamics_phase == "production" and np.any(protein_carrier):
                box = np.array([x_len, y_len, z_len], dtype=np.float64)
                carrier_index = np.where(protein_carrier)[0]
                for start in range(0, lipid_beads.size, 512):
                    delta = (
                        final_positions[carrier_index, None, :]
                        - final_positions[lipid_beads[start:start + 512]][None, :, :]
                    )
                    delta -= box * np.rint(delta / box)
                    contact_count[carrier_index] += np.sum(
                        np.sum(delta * delta, axis=2) < interface_cutoff * interface_cutoff,
                        axis=1,
                    )
            protein_friction = protein_carrier[brownian_atom_index]
            if dynamics_phase == "production":
                friction[protein_friction] = (
                    particle_friction * contact_count[brownian_atom_index[protein_friction]]
                )
            else:
                friction[protein_friction] = 0.0
            brownian_grp = t.create_group(input_grp, 'brownian')
            t.create_array(brownian_grp, 'atom_index', obj=brownian_atom_index)
            t.create_array(brownian_grp, 'friction', obj=friction.astype(np.float32))
            brownian_grp._v_attrs.numerical_time_step = np.float64(numerical_time_step)
            brownian_grp._v_attrs.protein_time_ps_per_step = np.float64(protein_time_ps)
            brownian_grp._v_attrs.martini_time_factor = np.float64(martini_time_factor)
            brownian_grp._v_attrs.raw_martini_time_ps_per_step = np.float64(raw_martini_time_ps)
            brownian_grp._v_attrs.target_dopc_diffusion_um2_s = np.float64(target_diffusion)
            brownian_grp._v_attrs.raw_target_diffusion_um2_s = np.float64(raw_target_diffusion)
            brownian_grp._v_attrs.reference_temperature_up = np.float64(reference_temperature)
            brownian_grp._v_attrs.dynamics_phase = dynamics_phase.encode()
            if dynamics_phase == "production":
                contact_atom_index = np.where(contact_count > 0)[0].astype('i4')
                t.create_array(brownian_grp, 'protein_contact_atom_index', obj=contact_atom_index)
                t.create_array(
                    brownian_grp,
                    'protein_lipid_contact_count',
                    obj=contact_count[contact_atom_index],
                )
                brownian_grp._v_attrs.bare_particle_diffusion_up = np.float64(particle_diffusion_up)
                brownian_grp._v_attrs.bare_particle_friction_up = np.float64(particle_friction)
                brownian_grp._v_attrs.interface_friction_cutoff_A = np.float64(interface_cutoff)
                brownian_grp._v_attrs.transport_observable = b"bare_martini_particle_lateral_diffusion"
            else:
                brownian_grp._v_attrs.raw_langevin_relaxation_ps = np.float64(raw_relaxation_ps)
                brownian_grp._v_attrs.relaxation_time_up = np.float64(relaxation_up)
                brownian_grp._v_attrs.time_unit_ps = np.float64(time_unit_ps)
                if dynamics_phase == "overlap_settling":
                    brownian_grp._v_attrs.transport_observable = b"overlap_settling_damping"
                else:
                    brownian_grp._v_attrs.transport_observable = b"native_dry_martini_relaxation"

        # Create charges array
        charge_array = t.create_array(input_grp, 'charges', obj=charges)
        charge_array._v_attrs.arguments = np.array([b'charges'])
        charge_array._v_attrs.shape = charges.shape
        charge_array._v_attrs.n_atoms = n_atoms
        charge_array._v_attrs.initialized = True
        
        # Create residue IDs array
        residue_array = t.create_array(input_grp, 'residue_ids', obj=residue_ids)
        residue_array._v_attrs.arguments = np.array([b'residue_ids'])
        residue_array._v_attrs.shape = residue_ids.shape
        residue_array._v_attrs.n_atoms = n_atoms
        residue_array._v_attrs.initialized = True
        residue_array._v_attrs.description = b"PDB residue indices per atom"

        # Create molecule IDs array (per-atom molecule index)
        molecule_array = t.create_array(input_grp, 'molecule_ids', obj=molecule_ids)
        molecule_array._v_attrs.arguments = np.array([b'molecule_ids'])
        molecule_array._v_attrs.shape = molecule_ids.shape
        molecule_array._v_attrs.n_atoms = n_atoms
        molecule_array._v_attrs.initialized = True
        molecule_array._v_attrs.description = b"Per-atom molecule index (0-based, contiguous)"

        # Store atom names for MARTINI-specific selections (e.g., BB backbone)
        atom_name_array = t.create_array(input_grp, 'atom_names', obj=atom_names.astype('S4'))
        atom_name_array._v_attrs.arguments = np.array([b'atom_names'])
        atom_name_array._v_attrs.shape = atom_names.shape
        atom_name_array._v_attrs.n_atoms = n_atoms
        atom_name_array._v_attrs.initialized = True

        # Store atom role names from PDB (BB, SC1, SC2, W, etc.)
        atom_roles = atom_names.astype('S4')
        atom_role_array = t.create_array(input_grp, 'atom_roles', obj=atom_roles)
        atom_role_array._v_attrs.arguments = np.array([b'atom_roles'])
        atom_role_array._v_attrs.shape = atom_roles.shape
        atom_role_array._v_attrs.n_atoms = n_atoms
        atom_role_array._v_attrs.initialized = True
        atom_role_array._v_attrs.description = b"PDB atom role names (BB, SC1, SC2, W, etc.)"

        if protein_sequence:
            t.create_array(input_grp, 'sequence', obj=np.asarray([x.encode('ascii') for x in protein_sequence], dtype='S3'))
        
        # Create potential group (required by UPSIDE)
        potential_grp = t.create_group(input_grp, 'potential')

        hybrid_position = t.create_group(potential_grp, 'martini_hybrid_position')
        hybrid_position._v_attrs.arguments = np.array([b'pos'])
        
        # Create MARTINI potential with proper parameters
        martini_potential = t.create_group(potential_grp, 'martini_potential')
        martini_potential._v_attrs.arguments = np.array([b'martini_hybrid_position'])
        martini_potential._v_attrs.potential_type = b'lj_coulomb'
        # coulomb_constant feeds the Python soft-core builder for the equilibration
        # stages; the runtime spline table already has Coulomb baked in.
        martini_potential._v_attrs.coulomb_constant = coulomb_constant_native * (length_conversion / energy_conversion)
        martini_potential._v_attrs.n_types = 1
        martini_potential._v_attrs.n_params = 4
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 2.0   # Verlet-list skin; rebuild is O(N) cell list, so keep the active list lean
        martini_potential._v_attrs.initialized = True

        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        martini_potential._v_attrs.coulomb_soften = params['coulomb_soften']
        if params['coulomb_soften']:
            martini_potential._v_attrs.slater_alpha = params['slater_alpha']
        martini_potential._v_attrs.lj_soften = params['lj_soften']
        if params['lj_soften']:
            martini_potential._v_attrs.lj_soften_alpha = params['lj_alpha']
        
        # Periodic boundary potential removed - using NVT ensemble without boundaries


        excluded_12_count, excluded_additional_count = _write_nonbonded_pairs(
            t, martini_potential, n_atoms, atom_types, charges,
            bonds_list, protein_exclusions, topo['martini_table'],
            energy_conversion, length_conversion,
        )

        # Add bonded potentials mirroring original run_martini.py
        # Bonds: dist_spring
        if bonds_list:
            bond_group = t.create_group(potential_grp, 'dist_spring')
            bond_group._v_attrs.arguments = np.array([b'pos'])
            bond_group._v_attrs.initialized = True

            bond_group._v_attrs.x_len = x_len
            bond_group._v_attrs.y_len = y_len
            bond_group._v_attrs.z_len = z_len
            t.create_array(bond_group, 'id', obj=np.array(bonds_list, dtype=int))
            t.create_array(bond_group, 'equil_dist', obj=np.array(bond_lengths_list, dtype='f4'))
            t.create_array(bond_group, 'spring_const', obj=np.array(bond_force_constants_list, dtype='f4'))
            # Compatibility dataset expected by some builds
            t.create_array(bond_group, 'bonded_atoms', obj=np.ones(len(bonds_list), dtype='i4'))

        # Angles: angle_spring
        if angles_list:
            angle_group = t.create_group(potential_grp, 'angle_spring')
            angle_group._v_attrs.arguments = np.array([b'pos'])
            angle_group._v_attrs.initialized = True

            angle_group._v_attrs.x_len = x_len
            angle_group._v_attrs.y_len = y_len
            angle_group._v_attrs.z_len = z_len
            t.create_array(angle_group, 'id', obj=np.array(angles_list, dtype=int))
            t.create_array(angle_group, 'equil_angle_deg', obj=np.array(angle_equil_deg_list, dtype='f4'))
            t.create_array(angle_group, 'spring_const', obj=np.array(angle_force_constants_list, dtype='f4'))

        # Dihedrals: dihedral_spring (if present)
        if dihedrals_list:
            dihedral_group = t.create_group(potential_grp, 'dihedral_spring')
            dihedral_group._v_attrs.arguments = np.array([b'pos'])
            dihedral_group._v_attrs.initialized = True

            dihedral_group._v_attrs.x_len = x_len
            dihedral_group._v_attrs.y_len = y_len
            dihedral_group._v_attrs.z_len = z_len
            t.create_array(dihedral_group, 'id', obj=np.array(dihedrals_list, dtype=int))
            eq_deg = np.array(dihedral_equil_deg_list, dtype='f4')
            t.create_array(dihedral_group, 'equil_angle_deg', obj=eq_deg)
            t.create_array(dihedral_group, 'spring_const', obj=np.array(dihedral_force_constants_list, dtype='f4'))
            t.create_array(dihedral_group, 'dihedral_type', obj=np.array(dihedral_type_list, dtype=int))

        # Position restraints (dry MARTINI lipid-head ramp before production)
        if lipid_restraint_indices:
            restraint_group = t.create_group(potential_grp, 'restraint_position')
            restraint_group._v_attrs.arguments = np.array([b'pos'])
            restraint_group._v_attrs.initialized = True
            t.create_array(
                restraint_group, 'restraint_indices', obj=np.array(lipid_restraint_indices, dtype='i4')
            )
            t.create_array(
                restraint_group, 'ref_pos', obj=np.array(lipid_restraint_ref_pos, dtype='f4')
            )
            spring_xyz = np.array(lipid_restraint_spring_xyz, dtype='f4')
            t.create_array(restraint_group, 'spring_const_xyz', obj=spring_xyz)

    print(f"Created UPSIDE input file: {input_file}")
    print("Preparation complete")
    
    summary_file = f"{run_dir}/preparation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("MARTINI preparation summary\n")
        f.write(f"PDB ID: {pdb_id}\n")
        f.write(f"Total atoms: {n_atoms}\n")
        f.write(f"Box dimensions: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} A\n")
        f.write(f"Total bonds: {len(bonds_list)}\n")
        f.write(f"Total angles: {len(angles_list)}\n")
        f.write(f"Total dihedrals: {len(dihedrals_list)}\n")
        f.write(f"Protein bonds: {protein_bond_count}\n")
        f.write(f"Protein constraints: {protein_constraint_count}\n")
        f.write(f"Protein angles: {protein_angle_count}\n")
        f.write(f"Protein dihedrals: {protein_dihedral_count}\n")
        f.write(f"DOPC lipids: {dopc_count}\n")
        f.write(f"Water molecules: {water_count}\n")
        f.write(f"1-2 exclusions (nrexcl=1): {excluded_12_count}\n")
        f.write(f"Additional exclusions: {excluded_additional_count}\n")

    print(f"Preparation summary saved to: {summary_file}")



def require_group(root, path):
    if path not in root:
        raise ValueError(f"Missing group: {path}")
    return root[path]


def require_dataset(group, name):
    if name not in group:
        raise ValueError(f"Missing dataset: {group.name}/{name}")
    return group[name]


def decode_attr_string(value, default=""):
    if value is None:
        return default
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def validate_hybrid_mapping(mapping_h5: Path, n_atom: int | None = None):
    path = Path(mapping_h5).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)

    with h5py.File(path, "r") as h5:
        inp = require_group(h5, "/input")
        ctrl = require_group(inp, "hybrid_control")
        bb = require_group(inp, "hybrid_bb_map")
        env = require_group(inp, "hybrid_env_topology")

        for attr in [
            "activation_stage",
            "preprod_protein_mode",
            "exclude_intra_protein_martini",
        ]:
            if attr not in ctrl.attrs:
                raise ValueError(f"Missing hybrid_control attr: {attr}")

        bb_atom_idx = require_dataset(bb, "bb_atom_index")[:]
        bb_atom_map = require_dataset(bb, "atom_indices")[:]
        bb_mask = require_dataset(bb, "atom_mask")[:]
        bb_w = require_dataset(bb, "weights")[:]

        if bb_atom_map.ndim != 2 or bb_atom_map.shape[1] != 4:
            raise ValueError("hybrid_bb_map/atom_indices must have shape (n_bb,4)")
        if bb_mask.shape != bb_atom_map.shape or bb_w.shape != bb_atom_map.shape:
            raise ValueError("hybrid_bb_map mask/weights shapes must match atom_indices")

        n_bb = bb_atom_map.shape[0]
        if bb_atom_idx.shape != (n_bb,):
            raise ValueError("hybrid_bb_map/bb_atom_index must have shape (n_bb,)")
        bb_index_space = decode_attr_string(bb.attrs.get("atom_index_space", "runtime_n_atom"))
        bb_reference_space = bb_index_space == "protein_aa_pdb_0based"

        if "reference_atom_names" in bb and bb["reference_atom_names"][:].shape != (4,):
            raise ValueError("hybrid_bb_map/reference_atom_names must have shape (4,)")
        if "reference_atom_indices" in bb and bb["reference_atom_indices"][:].shape != (n_bb, 4):
            raise ValueError("hybrid_bb_map/reference_atom_indices must have shape (n_bb,4)")
        if "reference_atom_coords" in bb and bb["reference_atom_coords"][:].shape != (n_bb, 4, 3):
            raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")
        if "bb_comment" in bb and bb["bb_comment"][:].shape != (n_bb,):
            raise ValueError("hybrid_bb_map/bb_comment must have shape (n_bb,)")

        for i in range(n_bb):
            mask = bb_mask[i].astype(bool)
            if mask.any():
                wsum = float(bb_w[i][mask].sum())
                if abs(wsum - 1.0) > 1e-4:
                    raise ValueError(f"BB weights do not sum to 1 for row {i}: {wsum}")

        membership = require_dataset(env, "protein_membership")[:]
        if membership.ndim != 1:
            raise ValueError("hybrid_env_topology/protein_membership must be 1D")
        if n_atom is not None and membership.shape[0] != n_atom:
            raise ValueError(
                f"protein_membership length ({membership.shape[0]}) != expected n_atom ({n_atom})"
            )
        n_atom_runtime = int(membership.shape[0])

        if "chain_break" in inp:
            chain_break = require_group(inp, "chain_break")
            chain_first_residue = require_dataset(chain_break, "chain_first_residue")[:]
            chain_counts = require_dataset(chain_break, "chain_counts")[:]
            if chain_first_residue.ndim != 1:
                raise ValueError("chain_break/chain_first_residue must be 1D")
            if chain_counts.ndim != 1:
                raise ValueError("chain_break/chain_counts must be 1D")
            if chain_counts.shape[0] != chain_first_residue.shape[0] + 1:
                raise ValueError("chain_break/chain_counts length must equal number of chains")
            if chain_first_residue.size and np.any(chain_first_residue[:-1] >= chain_first_residue[1:]):
                raise ValueError("chain_break/chain_first_residue must be strictly increasing")
            if chain_first_residue.size and (
                np.any(chain_first_residue <= 0) or np.any(chain_first_residue >= n_bb)
            ):
                raise ValueError("chain_break/chain_first_residue entries must be within (0, n_bb)")
            if np.any(chain_counts <= 0):
                raise ValueError("chain_break/chain_counts entries must be positive")

        for i in range(n_bb):
            bb_i = int(bb_atom_idx[i])
            if bb_i >= 0:
                if bb_i >= n_atom_runtime:
                    raise ValueError(f"BB proxy index out of bounds at row {i}: {bb_i}")
                if membership[bb_i] < 0:
                    raise ValueError(f"BB proxy index is not protein atom at row {i}: {bb_i}")
            elif bb_i != -1:
                raise ValueError(f"BB proxy index must be >=0 or -1 sentinel at row {i}: {bb_i}")

            if int(bb_mask[i, 1]) == 0:
                raise ValueError(f"BB map row {i} must include a CA carrier (atom_mask[:,1]==1)")
            for j in range(4):
                if int(bb_mask[i, j]) == 0:
                    continue
                ai = int(bb_atom_map[i, j])
                if bb_reference_space:
                    if ai < 0:
                        raise ValueError(
                            f"BB reference target index invalid at row {i}, col {j}: {ai}"
                        )
                else:
                    if ai < 0 or ai >= n_atom_runtime:
                        raise ValueError(f"BB target index out of bounds at row {i}, col {j}: {ai}")
                    if membership[ai] < 0:
                        raise ValueError(f"BB target index is not protein atom at row {i}, col {j}: {ai}")

        n_protein = int(np.sum(membership >= 0))
        n_env = int(np.sum(membership < 0))
        print(f"OK: {path}")
        print(f"  n_atom={membership.shape[0]} n_protein={n_protein} n_env={n_env}")
        print(f"  n_bb={n_bb} bb_index_space={bb_index_space}")


def decode_stage_label(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip().lower()
    return str(value).strip().lower()


def env_bool(name, default=False):
    raw = os.getenv(name)
    if raw is None:
        return default
    value = str(raw).strip().lower()
    return value not in ("0", "false", "no", "off", "")


def env_float(name, default):
    raw = os.getenv(name)
    if raw is None:
        return float(default)
    value = str(raw).strip()
    if not value:
        return float(default)
    return float(value)


def recenter_protein_for_production(h5f, pos, box_lengths):
    if box_lengths is None:
        return pos
    if "/input/stage_parameters" not in h5f:
        return pos
    stage = decode_stage_label(h5f["/input/stage_parameters"].attrs.get("current_stage", b""))
    if stage != "production":
        return pos
    if "/input/hybrid_bb_map/bb_atom_index" not in h5f:
        return pos

    box = np.asarray(box_lengths, dtype=np.float64).reshape(3)
    if not np.all(np.isfinite(box)) or np.any(box <= 0.0):
        return pos

    bb_idx = h5f["/input/hybrid_bb_map/bb_atom_index"][:].astype(np.int32)
    bb_idx = bb_idx[(bb_idx >= 0) & (bb_idx < pos.shape[0])]
    if bb_idx.size == 0:
        return pos

    bb_xyz = pos[bb_idx, :, 0].astype(np.float64)
    bb_center = np.zeros(3, dtype=np.float64)
    for ax in range(3):
        coords_ax = np.mod(bb_xyz[:, ax], box[ax])
        angles = TWOPI * coords_ax / box[ax]
        s = np.mean(np.sin(angles))
        c = np.mean(np.cos(angles))
        if abs(s) < 1e-12 and abs(c) < 1e-12:
            bb_center[ax] = float(np.mean(coords_ax))
        else:
            ang = np.arctan2(s, c)
            if ang < 0.0:
                ang += TWOPI
            bb_center[ax] = box[ax] * ang / TWOPI

    target_center = 0.5 * box
    delta = target_center - bb_center
    pos[:, :, 0] = pos[:, :, 0] + delta.astype(pos.dtype, copy=False)[None, :]
    print(
        "Recentered production protein BB to box center: "
        f"before=({bb_center[0]:.3f}, {bb_center[1]:.3f}, {bb_center[2]:.3f}) "
        f"target=({target_center[0]:.3f}, {target_center[1]:.3f}, {target_center[2]:.3f})"
    )
    return pos


def set_initial_position(input_file, output_file):
    strict_copy = env_bool("UPSIDE_SET_INITIAL_STRICT_COPY", False)
    apply_recenter_production = env_bool(
        "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION", default=(not strict_copy)
    )
    preserve_hybrid_transition = env_bool("UPSIDE_SET_INITIAL_PRESERVE_HYBRID_TRANSITION", False)
    time_step = float(os.environ.get("UPSIDE_SET_INITIAL_TIME_STEP", "0") or "0")

    with h5py.File(input_file, "r") as f:
        if "/output/pos" in f and f["/output/pos"].shape[0] > 0:
            last_pos = f["/output/pos"][-1, 0, :, :]
            last_pos = last_pos[:, :, np.newaxis]
        else:
            last_pos = f["/input/pos"][:, :, 0]
            last_pos = last_pos[:, :, np.newaxis]

        last_box = None
        if "/output/box" in f:
            box_data = f["/output/box"][:]
            if box_data.size > 0:
                last_box = box_data[-1]
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:
                    last_box = last_box[0]
        if last_box is None and "/input/potential/martini_potential" in f:
            pot_grp = f["/input/potential/martini_potential"]
            if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
                last_box = np.array(
                    [pot_grp.attrs["x_len"], pot_grp.attrs["y_len"], pot_grp.attrs["z_len"]]
                )
        last_mom = None
        if "/output/mom" in f and f["/output/mom"].shape[0] > 0:
            last_mom = f["/output/mom"][-1, 0, :, :]
            last_mom = last_mom[:, :, np.newaxis]

        source_stage = ""
        if "/input/stage_parameters" in f:
            source_stage = f["/input/stage_parameters"].attrs.get("current_stage", b"")
            if isinstance(source_stage, (bytes, np.bytes_)):
                source_stage = source_stage.decode("utf-8", errors="ignore")
            else:
                source_stage = str(source_stage)
            source_stage = source_stage.strip()
        source_transition_start = 0
        if "/input/hybrid_control" in f:
            source_transition_start = int(
                f["/input/hybrid_control"].attrs.get("sc_env_transition_step_start", 0)
            )
        source_elapsed_steps = 0
        if preserve_hybrid_transition and source_stage == "production" and time_step > 0.0:
            if "/output/time" in f and f["/output/time"].shape[0] > 0:
                last_time = float(f["/output/time"][-1])
                source_elapsed_steps = max(0, int(round(last_time / time_step)))
        exact_production_restart = preserve_hybrid_transition and source_stage == "production" and last_mom is not None

    with h5py.File(output_file, "r+") as f:
        target_stage = ""
        if "/input/stage_parameters" in f:
            target_stage = f["/input/stage_parameters"].attrs.get("current_stage", b"")
            if isinstance(target_stage, (bytes, np.bytes_)):
                target_stage = target_stage.decode("utf-8", errors="ignore")
            else:
                target_stage = str(target_stage)
            target_stage = target_stage.strip()

        target_pos = f["/input/pos"][:]
        target_n = target_pos.shape[0]
        source_n = last_pos.shape[0]
        if source_n != target_n:
            merged = target_pos.copy()
            merged[: min(source_n, target_n), :, :] = last_pos[: min(source_n, target_n), :, :]
            last_pos = merged

        if strict_copy and not apply_recenter_production:
            print("Strict handoff mode: preserving exact coordinates from previous stage output.")
        if exact_production_restart and apply_recenter_production:
            raise ValueError(
                "Production restart must preserve saved coordinates exactly when saved momentum is reused."
            )

        release_from_preprod = (source_stage != "production" and target_stage == "production")
        if apply_recenter_production:
            last_pos = recenter_protein_for_production(f, last_pos, last_box)

        if "/input/pos" in f:
            del f["/input/pos"]
        f.create_dataset("/input/pos", data=last_pos)

        if release_from_preprod:
            if "/input/mom" in f:
                del f["/input/mom"]
            mom_dtype = last_mom.dtype if last_mom is not None else target_pos.dtype
            mom_ds = f.create_dataset(
                "/input/mom",
                data=np.zeros((target_n, 3, 1), dtype=mom_dtype),
            )
            mom_ds.attrs["restart_valid"] = np.int8(0)
        elif last_mom is not None:
            if last_mom.shape[0] != target_n:
                merged_mom = f["/input/mom"][:] if "/input/mom" in f else np.zeros((target_n, 3, 1), dtype=last_mom.dtype)
                merged_mom[: min(last_mom.shape[0], target_n), :, :] = last_mom[: min(last_mom.shape[0], target_n), :, :]
                last_mom = merged_mom
            if "/input/mom" in f:
                del f["/input/mom"]
            mom_ds = f.create_dataset("/input/mom", data=last_mom)
            mom_ds.attrs["restart_valid"] = np.int8(1)
        elif "/input/mom" in f:
            f["/input/mom"].attrs["restart_valid"] = np.int8(0)

        if preserve_hybrid_transition and source_stage == "production" and "/input/hybrid_control" in f:
            transition_step = source_transition_start + source_elapsed_steps
            f["/input/hybrid_control"].attrs["sc_env_transition_step_start"] = np.int32(transition_step)

        if last_box is not None and "/input/potential" in f:
            n_updated = 0
            for pot_grp in f["/input/potential"].values():
                if not isinstance(pot_grp, h5py.Group):
                    continue
                if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
                    pot_grp.attrs["x_len"] = float(last_box[0])
                    pot_grp.attrs["y_len"] = float(last_box[1])
                    pot_grp.attrs["z_len"] = float(last_box[2])
                    n_updated += 1
            print(
                f"Updated box dimensions for {n_updated} potential nodes: "
                f"x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}"
            )


def normalize_resname(name: str) -> str:
    name = name.upper()
    aliases = {
        "HSD": "HIS",
        "HSE": "HIS",
        "HSP": "HIS",
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "CYX": "CYS",
    }
    return aliases.get(name, name)


def decode_string_array(dataset):
    return [x.decode("ascii") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in dataset[:]]


def unique_preserving_order(values):
    out = []
    seen = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def parse_itp_residue_names(path: Path):
    resnames = []
    seen = set()
    in_atoms = False
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low in {"[ atoms ]", "[atoms]"}:
                in_atoms = True
                continue
            if line.startswith("[") and line.endswith("]"):
                in_atoms = False
                continue
            if not in_atoms:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                resnr = int(parts[2])
            except ValueError:
                continue
            if parts[4].strip().upper() != "BB" or resnr in seen:
                continue
            seen.add(resnr)
            resnames.append(parts[3].strip().upper())
    return resnames


def resolve_sequence(inp, residue_count: int, protein_itp: Path | None):
    mismatch_notes = []
    if "sequence" in inp:
        sequence = [normalize_resname(x) for x in decode_string_array(inp["sequence"])]
        if len(sequence) == residue_count:
            return sequence
        mismatch_notes.append(f"/input/sequence has {len(sequence)} residues")

    if protein_itp is not None:
        if not protein_itp.exists():
            raise ValueError(f"Protein ITP not found for sequence fallback: {protein_itp}")
        sequence = [normalize_resname(x) for x in parse_itp_residue_names(protein_itp)]
        if len(sequence) == residue_count:
            return sequence
        mismatch_notes.append(f"{protein_itp} has {len(sequence)} residues")

    detail = ", ".join(mismatch_notes) if mismatch_notes else "no sequence source was available"
    raise ValueError(
        f"Could not resolve a residue sequence matching hybrid_bb_map residue count {residue_count}: {detail}"
    )


def build_affine_atoms(inp):
    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map in stage file")
    bb_grp = inp["hybrid_bb_map"]
    if "atom_indices" not in bb_grp:
        raise ValueError("hybrid_bb_map/atom_indices is required for stage-7 CB placement")

    bb_residue_raw = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_atom_idx = bb_grp["atom_indices"][:].astype(np.int32)
    residue_ids = unique_preserving_order(int(x) for x in bb_residue_raw.tolist())

    residue_to_ncac = {}
    n_atom = int(inp["pos"].shape[0])
    for resid, atom_row in zip(bb_residue_raw.tolist(), bb_atom_idx.tolist()):
        rid = int(resid)
        n_idx, ca_idx, c_idx = [int(atom_row[0]), int(atom_row[1]), int(atom_row[2])]
        if n_idx < 0 or ca_idx < 0 or c_idx < 0:
            continue
        for idx in (n_idx, ca_idx, c_idx):
            if idx < 0 or idx >= n_atom:
                raise ValueError(
                    f"Backbone carrier index out of bounds for residue {rid}: idx={idx}, n_atom={n_atom}"
                )
        residue_to_ncac.setdefault(rid, (n_idx, ca_idx, c_idx))

    missing = [rid for rid in residue_ids if rid not in residue_to_ncac]
    if missing:
        raise ValueError(f"Missing N/CA/C runtime mapping for residues: {missing[:8]}")

    affine_atoms = np.zeros((len(residue_ids), 3), dtype=np.int32)
    for seq_idx, resid in enumerate(residue_ids):
        affine_atoms[seq_idx, :] = residue_to_ncac[resid]
    return residue_ids, affine_atoms


def load_sidechain_rotamer_payload(sidechain_lib: Path, sequence, restype_to_index):
    if not sidechain_lib.exists():
        raise ValueError(f"Sidechain library not found: {sidechain_lib}")

    with h5py.File(sidechain_lib, "r") as sclib:
        sc_restype_order = decode_string_array(sclib["restype_order"])
        sc_restype_num = {name: i for i, name in enumerate(sc_restype_order)}
        start_stop_bead = sclib["rotamer_start_stop_bead"][:].astype(np.int32)
        rot_center_fixed = sclib["rotamer_center_fixed"][:, :6].astype(np.float32)
        if "rotamer_prob_fixed" in sclib:
            rot_energy_fixed = sclib["rotamer_prob_fixed"][:].astype(np.float32).reshape(-1)
        elif "rotamer_prob" in sclib:
            rot_prob = sclib["rotamer_prob"][:].astype(np.float64)
            if rot_prob.ndim != 3:
                raise ValueError(f"Unsupported rotamer_prob shape in {sidechain_lib}: {rot_prob.shape}")
            rot_prob_mean = np.clip(rot_prob.mean(axis=(0, 1)), 1.0e-12, None)
            rot_energy_fixed = (-np.log(rot_prob_mean)).astype(np.float32)
        else:
            raise ValueError(
                f"Missing rotamer probability tables in {sidechain_lib}: need rotamer_prob_fixed or rotamer_prob"
            )
        bead_order = decode_string_array(sclib["bead_order"])
        bead_num = {name: i for i, name in enumerate(bead_order)}
        pair_interaction = sclib["pair_interaction"][:].astype(np.float32)

    count_by_n_rot = {}
    affine_residue = []
    layer_index = []
    beadtype_seq = []
    bead_type_index = []
    id_seq = []
    row_rotamer_index = []
    row_residue_table_index = []
    skipped = []

    for seq_idx, raw_resname in enumerate(sequence):
        resname = normalize_resname(raw_resname)
        restype_idx = sc_restype_num.get(resname)
        residue_table_idx = restype_to_index.get(resname)
        if restype_idx is None or residue_table_idx is None:
            skipped.append((seq_idx, resname))
            continue

        start, stop, n_bead = [int(x) for x in start_stop_bead[restype_idx]]
        if n_bead <= 0 or stop <= start:
            skipped.append((seq_idx, resname))
            continue
        n_rot = (stop - start) // n_bead
        if n_rot <= 0:
            skipped.append((seq_idx, resname))
            continue

        if n_rot not in count_by_n_rot:
            count_by_n_rot[n_rot] = 0
        base_id = (count_by_n_rot[n_rot] << N_BIT_ROTAMER) + n_rot
        count_by_n_rot[n_rot] += 1

        for rel in range(stop - start):
            lid = start + rel
            rot_idx = rel // n_bead
            bead_name = f"{resname}_{rel % n_bead}"
            if bead_name not in bead_num:
                raise ValueError(f"Missing bead '{bead_name}' in sidechain library bead_order")
            affine_residue.append(seq_idx)
            layer_index.append(lid)
            beadtype_seq.append(bead_name)
            bead_type_index.append(bead_num[bead_name])
            id_seq.append(rot_idx + (base_id << N_BIT_ROTAMER))
            row_rotamer_index.append(rot_idx)
            row_residue_table_index.append(residue_table_idx)

    if not layer_index:
        raise ValueError("No production rotamer rows could be generated from the sidechain library")

    layer_index_arr = np.asarray(layer_index, dtype=np.int32)
    return {
        "pair_interaction": pair_interaction,
        "affine_residue": np.asarray(affine_residue, dtype=np.int32),
        "layer_index": layer_index_arr,
        "beadtype_seq": np.asarray([np.bytes_(x) for x in beadtype_seq], dtype="S16"),
        "bead_type_index": np.asarray(bead_type_index, dtype=np.int32),
        "id_seq": np.asarray(id_seq, dtype=np.int32),
        "placement_data": rot_center_fixed[layer_index_arr].astype(np.float32),
        "placement_scalar_data": rot_energy_fixed[layer_index_arr][:, None].astype(np.float32),
        "row_rotamer_index": np.asarray(row_rotamer_index, dtype=np.int32),
        "row_residue_table_index": np.asarray(row_residue_table_index, dtype=np.int32),
        "skipped": skipped,
    }


def build_env_rows(inp, target_to_index):
    if "hybrid_env_topology" not in inp:
        raise ValueError("Missing /input/hybrid_env_topology in stage file")
    env_grp = inp["hybrid_env_topology"]
    if "protein_membership" not in env_grp:
        raise ValueError("Missing hybrid_env_topology/protein_membership in stage file")

    protein_membership = env_grp["protein_membership"][:].astype(np.int32)
    atom_types = decode_string_array(inp["type"])

    env_atom_index = []
    env_target_index = []
    for atom_idx, (membership, atom_type) in enumerate(zip(protein_membership.tolist(), atom_types)):
        if membership >= 0:
            continue
        target_idx = target_to_index.get(atom_type)
        if target_idx is None:
            continue
        env_atom_index.append(atom_idx)
        env_target_index.append(target_idx)

    return np.asarray(env_atom_index, dtype=np.int32), np.asarray(env_target_index, dtype=np.int32)


def recreate_group(parent, name):
    if name in parent:
        del parent[name]
    return parent.create_group(name)


def load_upside_config(upside_home: Path):
    upside_config_py = upside_home / "py" / "upside_config.py"
    if not upside_config_py.exists():
        raise ValueError(f"upside_config.py not found: {upside_config_py}")
    spec = importlib.util.spec_from_file_location("upside_config_runtime", str(upside_config_py))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def remap_atom_index_array(arr, atom_map):
    out = arr.copy()
    if out.size == 0:
        return out
    finite_mask = np.isfinite(out)
    if not np.any(finite_mask):
        return out
    nonneg_mask = finite_mask & (out >= 0)
    if not np.any(nonneg_mask):
        return out
    idx = out[nonneg_mask].astype(np.int64)
    if np.any(idx >= atom_map.shape[0]):
        raise ValueError(
            f"Backbone atom remap index exceeds reference map size: max={int(idx.max())}, map_size={atom_map.shape[0]}"
        )
    mapped = atom_map[idx]
    if np.any(mapped < 0):
        bad = np.where(mapped < 0)[0][0]
        raise ValueError(
            f"Backbone atom remap produced negative target for reference index {int(idx[bad])}"
        )
    out[nonneg_mask] = mapped.astype(out.dtype, copy=False)
    return out


def inject_backbone_nodes(
    up_file: Path,
    sequence,
    affine_atoms,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
    upside_home: Path,
):
    uc = load_upside_config(upside_home)
    fasta_seq = np.asarray(sequence)
    ref_n_atom = 3 * len(sequence)
    spring_args = types.SimpleNamespace(bond_stiffness=48.0, angle_stiffness=175.0, omega_stiffness=30.0)
    chain_first_residue = np.array([], dtype=np.int32)

    with h5py.File(up_file, "r") as up:
        inp = up["input"]
        if "chain_break" in inp and "chain_first_residue" in inp["chain_break"]:
            chain_first_residue = inp["chain_break"]["chain_first_residue"][:].astype(np.int32)
        elif "hybrid_bb_map" in inp and "bb_chain_id" in inp["hybrid_bb_map"]:
            bb_chain_ids = decode_string_array(inp["hybrid_bb_map"]["bb_chain_id"])
            chain_first_residue, _ = derive_chain_break_metadata(
                [{"bb_chain": chain_id} for chain_id in bb_chain_ids]
            )
    n_chains = int(chain_first_residue.size) + 1
    chain_starts = np.append(np.array([0], dtype=np.int32), chain_first_residue) * 3
    # Mirror upside_config chain-break handling: exclude residues adjacent to
    # each chain boundary from inferred H-bond donor/acceptor construction.
    # For each chain first residue i (except chain 0), exclude residues i-1 and i.
    hbond_excluded_residues = np.array(
        sorted({int(i + j) for i in chain_first_residue.tolist() for j in (-1, 0)}),
        dtype=np.int32,
    )

    with tb.open_file(str(up_file), mode="a") as tf:
        uc.t = tf
        uc.potential = tf.root.input.potential
        uc.n_atom = ref_n_atom
        uc.n_chains = n_chains
        uc.chain_first_residue = chain_first_residue.astype(np.int32, copy=False)
        uc.chain_starts = chain_starts.astype(np.int32, copy=False)
        uc.rl_chains = None
        uc.use_intensive_memory = False

        uc.write_dist_spring(spring_args)
        uc.write_angle_spring(spring_args)
        uc.write_omega_spring1(spring_args, fasta_seq)
        uc.write_rama_map_pot(
            fasta_seq,
            str(rama_library),
            str(rama_sheet_mixing),
            secstr_bias="",
            mode="mixture",
            param_deriv=False,
        )

        with open(reference_state_rama, "rb") as fh:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"dtype\(\): align should be passed as Python or NumPy boolean.*",
                    category=np.exceptions.VisibleDeprecationWarning,
                )
                ref_state_raw = cPickle.load(fh, encoding="latin1")
        ref_state_cor = np.log(np.asarray(ref_state_raw, dtype=np.float64))
        ref_state_cor -= ref_state_cor.mean()
        grp = tf.create_group(tf.root.input.potential, "rama_map_pot_ref")
        grp._v_attrs.arguments = np.array([b"rama_coord"])
        grp._v_attrs.log_pot = 0
        uc.create_array(grp, "residue_id", obj=np.arange(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_map_id", obj=np.zeros(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_pot", obj=ref_state_cor[None])

        uc.write_infer_H_O(fasta_seq, hbond_excluded_residues)
        uc.write_count_hbond(fasta_seq, False)
        uc.write_short_hbond(fasta_seq, str(hbond_energy))
        uc.write_rama_coord2()
        uc.write_backbone_pair(fasta_seq)

    atom_map = np.full((ref_n_atom,), -1, dtype=np.int64)
    for seq_idx, (n_idx, ca_idx, c_idx) in enumerate(affine_atoms.tolist()):
        atom_map[3 * seq_idx + 0] = n_idx
        atom_map[3 * seq_idx + 1] = ca_idx
        atom_map[3 * seq_idx + 2] = c_idx

    remap_datasets = [
        ("Distance3D", "id"),
        ("Angle", "id"),
        ("Dihedral_omega", "id"),
        ("Dihedral_phi", "id"),
        ("Dihedral_psi", "id"),
        ("rama_coord", "id"),
        ("infer_H_O/donors", "id"),
        ("infer_H_O/acceptors", "id"),
    ]
    with h5py.File(up_file, "r+") as up:
        pot = up["/input/potential"]
        for node_name in BACKBONE_NODES:
            if node_name not in pot:
                raise ValueError(f"Missing generated backbone node in stage file {up_file}: {node_name}")
        for grp_name, dset_name in remap_datasets:
            ds = pot[grp_name][dset_name]
            ds[...] = remap_atom_index_array(ds[:], atom_map).astype(ds.dtype, copy=False)


def inject_particles_table(up_file: Path, martini_h5: Path):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()

    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")

    with h5py.File(martini_h5, "r") as src:
        pg = src["particles"]
        eps_arr = pg["unique_eps_eup"][:].astype(np.float64)
        sig_arr = pg["unique_sig_ang"][:].astype(np.float64)
        qq_arr = pg["unique_charge_product"][:].astype(np.float64)
        combined_grids = pg["combined_energy_grids"][:].astype(np.float64)
        n_triples = len(eps_arr)

    has_zero_triple = np.any(
        (np.abs(eps_arr) < 1e-12)
        & (np.abs(sig_arr) < 1e-12)
        & (np.abs(qq_arr) < 1e-12)
    )
    if not has_zero_triple:
        eps_arr = np.concatenate([eps_arr, np.array([0.0], dtype=np.float64)])
        sig_arr = np.concatenate([sig_arr, np.array([0.0], dtype=np.float64)])
        qq_arr = np.concatenate([qq_arr, np.array([0.0], dtype=np.float64)])
        combined_grids = np.vstack([
            combined_grids,
            np.zeros((1, combined_grids.shape[1]), dtype=np.float64),
        ])
        n_triples = len(eps_arr)

    GRID_N = 1000
    R_MIN = 0.0
    R_MAX = 12.0

    with h5py.File(up_file, "r+") as up:
        g = up["/input/potential/martini_potential"]

        lj_soften = bool(int(g.attrs.get("lj_soften", 0)))
        lj_soften_alpha = float(g.attrs.get("lj_soften_alpha", 0.0))
        coulomb_soften = bool(int(g.attrs.get("coulomb_soften", 0)))
        slater_alpha = float(g.attrs.get("slater_alpha", 0.0))
        has_softening = (lj_soften and lj_soften_alpha > 0.0) or coulomb_soften

        if has_softening:
            coulomb_k = float(g.attrs["coulomb_constant"])
            lj_alpha = float(lj_soften_alpha) if lj_soften else 0.0
            slater = float(slater_alpha) if coulomb_soften else 0.0

            for idx in range(n_triples):
                eps = float(eps_arr[idx])
                sig = float(sig_arr[idx])
                qq = float(qq_arr[idx])
                for i in range(GRID_N):
                    r = R_MIN + i * (R_MAX - R_MIN) / (GRID_N - 1)
                    if r == 0.0:
                        r = 1.0e-6

                    # LJ component with optional soft-core
                    if eps != 0.0 and sig != 0.0:
                        if lj_soften and lj_alpha > 0.0:
                            x = r / sig
                            t = (x * x) ** 3 + lj_alpha
                            inv_t = 1.0 / t
                            lj = 4.0 * eps * (inv_t * inv_t - inv_t)
                        else:
                            x = sig / r
                            x6 = (x * x) ** 3
                            lj = 4.0 * eps * (x6 * x6 - x6)
                    else:
                        lj = 0.0

                    # Coulomb component with optional softening
                    if abs(qq) > 1e-10:
                        coul = coulomb_k * qq / r
                        if coulomb_soften:
                            alpha_r = slater * r
                            coul *= 1.0 - (1.0 + alpha_r * 0.5) * math.exp(-alpha_r)
                    else:
                        coul = 0.0

                    combined_grids[idx, i] = lj + coul

        if "coefficients" in g:
            coeffs = g["coefficients"][:].astype(np.float64)
            if coeffs.size:
                if coeffs.ndim != 2 or coeffs.shape[1] < 4:
                    raise ValueError(
                        f"{up_file}: martini_potential/coefficients must have shape (n,4)"
                    )
                local_eps = np.empty(coeffs.shape[0], dtype=np.float64)
                local_sig = np.empty(coeffs.shape[0], dtype=np.float64)
                local_qq = np.empty(coeffs.shape[0], dtype=np.float64)
                local_grids = np.empty((coeffs.shape[0], combined_grids.shape[1]), dtype=np.float64)
                for coeff_i, (eps, sig, q1, q2) in enumerate(coeffs[:, :4]):
                    qq = float(q1) * float(q2)
                    match = np.where(
                        np.isclose(eps_arr, eps, rtol=1e-5, atol=1e-7)
                        & np.isclose(sig_arr, sig, rtol=1e-5, atol=1e-6)
                        & np.isclose(qq_arr, qq, rtol=1e-5, atol=1e-7)
                    )[0]
                    if match.size == 0:
                        raise ValueError(
                            f"{up_file}: no particle spline row for local coefficient "
                            f"{coeff_i}: eps={eps:g}, sig={sig:g}, q1*q2={qq:g}"
                        )
                    src_i = int(match[0])
                    local_eps[coeff_i] = eps_arr[src_i]
                    local_sig[coeff_i] = sig_arr[src_i]
                    local_qq[coeff_i] = qq_arr[src_i]
                    local_grids[coeff_i, :] = combined_grids[src_i, :]
                eps_arr = local_eps
                sig_arr = local_sig
                qq_arr = local_qq
                combined_grids = local_grids
                n_triples = len(eps_arr)
                g.attrs["coefficient_indices_reference"] = "coefficients"

        for name in (
            "unique_eps_eup",
            "unique_sig_ang",
            "unique_charge_product",
            "combined_energy_grids",
        ):
            if name in g:
                del g[name]
        g.create_dataset("unique_eps_eup", data=eps_arr, dtype=np.float64)
        g.create_dataset("unique_sig_ang", data=sig_arr, dtype=np.float64)
        g.create_dataset("unique_charge_product", data=qq_arr, dtype=np.float64)
        g.create_dataset("combined_energy_grids", data=combined_grids, dtype=np.float64)
        g.attrs["n_combined_triples"] = np.int32(n_triples)
        g.attrs["r_min_ang"] = np.float32(R_MIN)
        g.attrs["r_max_ang"] = np.float32(R_MAX)
        g.attrs["n_points"] = np.int32(GRID_N)


def inject_stage7_sc_table_nodes(
    up_file: Path,
    martini_h5: Path,
    upside_home: Path,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
    protein_itp: Path | None = None,
):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()
    upside_home = Path(upside_home).expanduser().resolve()
    rama_library = Path(rama_library).expanduser().resolve()
    rama_sheet_mixing = Path(rama_sheet_mixing).expanduser().resolve()
    hbond_energy = Path(hbond_energy).expanduser().resolve()
    reference_state_rama = Path(reference_state_rama).expanduser().resolve()
    protein_itp = Path(protein_itp).expanduser().resolve() if protein_itp else None
    sidechain_lib = (upside_home / "parameters" / "ff_2.1" / "sidechain.h5").resolve()

    if not up_file.exists():
        raise SystemExit(f"ERROR: stage file not found: {up_file}")
    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")
    for path in [rama_library, rama_sheet_mixing, hbond_energy, reference_state_rama, sidechain_lib]:
        if not path.exists():
            raise SystemExit(f"ERROR: required Upside input not found: {path}")

    with h5py.File(martini_h5, "r") as sc_lib:
        sc_grp = sc_lib["sc_table"] if "sc_table" in sc_lib else sc_lib
        required_sc_datasets = [
            "restype_order",
            "target_order",
            "grid_ang",
            "cos_theta_grid",
            "rotamer_count",
            "rotamer_probability_fixed",
            "rotamer_radial_energy_eup",
            "rotamer_angular_energy_eup",
            "rotamer_angular_profile",
            "rotamer_full_energy_eup",
        ]
        missing_sc_datasets = [name for name in required_sc_datasets if name not in sc_grp]
        if missing_sc_datasets:
            missing_text = ", ".join(missing_sc_datasets)
            raise SystemExit(
                f"ERROR: {martini_h5} is missing required rotamer-resolved SC datasets: {missing_text}. "
                "Regenerate martini.h5 via build_martini_h5()."
            )
        restype_order = decode_string_array(sc_grp["restype_order"])
        target_order = decode_string_array(sc_grp["target_order"])
        grid_ang = sc_grp["grid_ang"][:].astype(np.float32)
        cos_theta_grid = sc_grp["cos_theta_grid"][:].astype(np.float32)
        rotamer_count = sc_grp["rotamer_count"][:].astype(np.int32)
        rotamer_probability_fixed = sc_grp["rotamer_probability_fixed"][:].astype(np.float32)
        rotamer_radial_energy_eup = sc_grp["rotamer_radial_energy_eup"][:].astype(np.float32)
        rotamer_angular_energy_eup = sc_grp["rotamer_angular_energy_eup"][:].astype(np.float32)
        rotamer_angular_profile = sc_grp["rotamer_angular_profile"][:].astype(np.float32)
        rotamer_full_energy_eup = sc_grp["rotamer_full_energy_eup"][:].astype(np.float32)

    restype_to_index = {name: i for i, name in enumerate(restype_order)}
    target_to_index = {name: i for i, name in enumerate(target_order)}

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]

        residue_ids, affine_atoms = build_affine_atoms(inp)
        sequence = resolve_sequence(inp, len(residue_ids), protein_itp)
        env_atom_index, env_target_index = build_env_rows(inp, target_to_index)
        rotamer_payload = load_sidechain_rotamer_payload(sidechain_lib, sequence, restype_to_index)

        if "hybrid_sc_map" in inp:
            del inp["hybrid_sc_map"]
        if "sequence" in inp:
            del inp["sequence"]
        inp.create_dataset("sequence", data=np.asarray([np.bytes_(x) for x in sequence], dtype="S3"))

        for node_name in [*LEGACY_STAGE7_NODES, *BACKBONE_NODES]:
            if node_name in pot:
                del pot[node_name]

        g_aff = recreate_group(pot, "affine_alignment")
        g_aff.attrs["arguments"] = np.asarray([np.bytes_("pos")])
        g_aff.create_dataset("atoms", data=affine_atoms, dtype=np.int32)
        g_aff.create_dataset(
            "ref_geom",
            data=np.repeat(CANONICAL_AFFINE_REF[None, :, :], len(sequence), axis=0),
            dtype=np.float32,
        )

        g_sc_place = recreate_group(pot, "placement_fixed_point_vector_only")
        g_sc_place.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_sc_place.create_dataset("rama_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc_place.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc_place.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_sc_place.create_dataset("placement_data", data=rotamer_payload["placement_data"], dtype=np.float32)
        g_sc_place.create_dataset("beadtype_seq", data=rotamer_payload["beadtype_seq"])
        g_sc_place.create_dataset("id_seq", data=rotamer_payload["id_seq"], dtype=np.int32)
        g_sc_place.create_dataset("fix_rotamer", data=np.zeros((0, 2), dtype=np.int32), dtype=np.int32)

        g_pl = recreate_group(pot, "placement_fixed_scalar")
        g_pl.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_pl.create_dataset("rama_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("affine_residue", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_pl.create_dataset("layer_index", data=rotamer_payload["layer_index"], dtype=np.int32)
        g_pl.create_dataset("placement_data", data=rotamer_payload["placement_scalar_data"], dtype=np.float32)

        g_rot = recreate_group(pot, "rotamer")
        g_rot.attrs["arguments"] = np.asarray(
            [
                np.bytes_("placement_fixed_point_vector_only"),
                np.bytes_("placement_fixed_scalar"),
                np.bytes_("martini_sc_table_1body"),
            ]
        )
        g_rot.attrs["integrator_level"] = np.int32(1)
        g_rot.attrs["max_iter"] = np.int32(1000)
        g_rot.attrs["tol"] = np.float32(1.0e-3)
        g_rot.attrs["damping"] = np.float32(0.4)
        g_rot.attrs["iteration_chunk_size"] = np.int32(2)

        g_pair = g_rot.create_group("pair_interaction")
        g_pair.create_dataset("interaction_param", data=rotamer_payload["pair_interaction"], dtype=np.float32)
        g_pair.create_dataset(
            "index", data=np.arange(len(rotamer_payload["id_seq"]), dtype=np.int32), dtype=np.int32
        )
        g_pair.create_dataset("type", data=rotamer_payload["bead_type_index"], dtype=np.int32)
        g_pair.create_dataset("id", data=rotamer_payload["id_seq"], dtype=np.int32)

        g_cb = recreate_group(pot, "placement_fixed_point_vector_only_CB")
        g_cb.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_cb.create_dataset("affine_residue", data=np.arange(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset("layer_index", data=np.zeros(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset(
            "placement_data", data=np.concatenate([CB_PLACEMENT, CB_VECTOR], axis=1), dtype=np.float32
        )

    inject_backbone_nodes(
        up_file=up_file,
        sequence=sequence,
        affine_atoms=affine_atoms,
        rama_library=rama_library,
        rama_sheet_mixing=rama_sheet_mixing,
        hbond_energy=hbond_energy,
        reference_state_rama=reference_state_rama,
        upside_home=upside_home,
    )

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]
        g_sc = recreate_group(pot, "martini_sc_table_1body")
        g_sc.attrs["arguments"] = np.asarray([np.bytes_("pos"), np.bytes_("placement_fixed_point_vector_only_CB")])
        g_sc.attrs["x_len"] = np.float32(martini_potential.attrs["x_len"])
        g_sc.attrs["y_len"] = np.float32(martini_potential.attrs["y_len"])
        g_sc.attrs["z_len"] = np.float32(martini_potential.attrs["z_len"])
        g_sc.create_dataset("row_residue_index", data=rotamer_payload["affine_residue"], dtype=np.int32)
        g_sc.create_dataset("row_rotamer_index", data=rotamer_payload["row_rotamer_index"], dtype=np.int32)
        g_sc.create_dataset(
            "row_residue_table_index", data=rotamer_payload["row_residue_table_index"], dtype=np.int32
        )
        g_sc.create_dataset("env_atom_index", data=env_atom_index, dtype=np.int32)
        g_sc.create_dataset("env_target_index", data=env_target_index, dtype=np.int32)
        g_sc.create_dataset("grid_ang", data=grid_ang, dtype=np.float32)
        g_sc.create_dataset("cos_theta_grid", data=cos_theta_grid, dtype=np.float32)
        g_sc.create_dataset("rotamer_count", data=rotamer_count.astype(np.int32), dtype=np.int32)
        g_sc.create_dataset("rotamer_probability_fixed", data=rotamer_probability_fixed, dtype=np.float32)
        g_sc.create_dataset("rotamer_radial_energy_eup", data=rotamer_radial_energy_eup, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_energy_eup", data=rotamer_angular_energy_eup, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_profile", data=rotamer_angular_profile, dtype=np.float32)
        g_sc.create_dataset("rotamer_full_energy_eup", data=rotamer_full_energy_eup, dtype=np.float32)
        g_sc.create_dataset("restype_order", data=np.asarray([np.bytes_(x) for x in restype_order], dtype="S4"))
        g_sc.create_dataset("target_order", data=np.asarray([np.bytes_(x) for x in target_order], dtype="S8"))

    print(
        f"Injected rotamer-weighted martini_sc_table_1body into {up_file}: "
        f"n_rows={len(rotamer_payload['id_seq'])} skipped={len(rotamer_payload['skipped'])} "
        f"n_env={len(env_atom_index)} n_restypes={len(restype_order)} n_targets={len(target_order)}"
    )


if __name__ == "__main__":
    raise SystemExit("Use martini_prepare_system.py as the workflow preparation entrypoint.")
