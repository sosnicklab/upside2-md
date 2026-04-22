#!/usr/bin/env python3
from __future__ import annotations

import _pickle as cPickle
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
import tables as tb

PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

NA_AVOGADRO = 6.02214076e23
BB_COMPONENT_NAMES = ("N", "CA", "C", "O")
PROTEIN_AA_RESNAMES = {
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
AA_BACKBONE_PLACEHOLDER_TYPE = "AABB"
AA_BACKBONE_ATOM_MASS = {
    "N": 14.0,
    "CA": 12.0,
    "C": 12.0,
    "O": 16.0,
}
NONINTERACTING_ATOM_TYPES = {AA_BACKBONE_PLACEHOLDER_TYPE}
DEFAULT_SC_TABLE_JSON = Path(
    "/Users/yinhan/Documents/upside2-md-martini/SC-training/runs/default/results/assembled/sc_table.json"
)
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
CB_PLACEMENT = np.array([[0.0, 0.94375626, 1.2068012]], dtype=np.float32)
CB_VECTOR = CB_PLACEMENT / np.linalg.norm(CB_PLACEMENT, axis=1, keepdims=True)
N_BIT_ROTAMER = 4
LEGACY_STAGE7_NODES = [
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
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
EXPLICIT_BACKBONE_OXYGEN_BOND_STIFFNESS = 48.0
EXPLICIT_BACKBONE_OXYGEN_ANGLE_STIFFNESS = 175.0

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


def extract_protein_aa_atoms(aa_atoms):
    protein_atoms = [deepcopy(atom) for atom in aa_atoms if atom["resname"].strip().upper() in PROTEIN_AA_RESNAMES]
    if not protein_atoms:
        raise ValueError("No protein AA atoms found in protein AA PDB.")
    return protein_atoms


def extract_protein_aa_backbone_atoms(aa_atoms):
    protein_atoms = extract_protein_aa_atoms(aa_atoms)
    backbone = []
    for atom in protein_atoms:
        atom_name = atom["name"].strip().upper()
        if atom_name in BB_COMPONENT_NAMES:
            backbone.append(atom)
    if not backbone:
        raise ValueError("No AA backbone atoms (N/CA/C/O) found in protein AA PDB.")
    return backbone


def infer_protein_charge_from_aa(protein_aa_atoms):
    charged_res = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
    seen = set()
    total = 0
    for atom in protein_aa_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen:
            continue
        seen.add(key)
        total += charged_res.get(atom["resname"].upper(), 0)
    return total


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


def build_sidechain_centroid_proxy_atoms(protein_aa_atoms):
    proxies = []
    for atoms in residue_group_atoms(protein_aa_atoms).values():
        sidechain_atoms = [atom for atom in atoms if atom["name"].strip().upper() not in BB_COMPONENT_NAMES]
        if not sidechain_atoms:
            continue
        centroid = center_of_mass(coords(sidechain_atoms))
        proxy = deepcopy(sidechain_atoms[0])
        proxy["name"] = "SC"
        proxy["x"], proxy["y"], proxy["z"] = (float(centroid[0]), float(centroid[1]), float(centroid[2]))
        proxies.append(proxy)
    return proxies


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


def resize_and_shift(all_atoms, pad_xy, pad_z):
    xyz = coords(all_atoms)
    min_c = xyz.min(axis=0)
    max_c = xyz.max(axis=0)
    spans = max_c - min_c
    box = np.array(
        [spans[0] + 2.0 * pad_xy, spans[1] + 2.0 * pad_xy, spans[2] + 2.0 * pad_z],
        dtype=float,
    )
    shift = np.array([pad_xy, pad_xy, pad_z], dtype=float) - min_c
    shifted = xyz + shift
    set_coords(all_atoms, shifted)
    return box


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
    types = [("NA", "NA", 1.0, "IONS")] * n_na + [("CL", "CL", -1.0, "IONS")] * n_cl
    for name, resname, _charge, segid in types:
        accepted = False
        for _ in range(20000):
            trial = rng.uniform([0, 0, 0], box_lengths)
            if existing.size:
                d2 = np.sum((existing - trial) ** 2, axis=1)
                if np.min(d2) < cutoff2:
                    continue
            if placed:
                placed_xyz = np.array([[a["x"], a["y"], a["z"]] for a in placed], dtype=float)
                d2_placed = np.sum((placed_xyz - trial) ** 2, axis=1)
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


def extract_backbone_sequence(backbone_atoms):
    sequence = []
    seen = set()
    for atom in backbone_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in seen:
            continue
        seen.add(key)
        sequence.append(atom["resname"].strip().upper())
    return sequence


def collect_aa_backbone_map(backbone_atoms):
    by_residue = defaultdict(dict)
    residue_order = []
    for idx, atom in enumerate(backbone_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"].strip().upper())
        if key not in by_residue:
            residue_order.append(key)
        by_residue[key][atom["name"].strip().upper()] = idx

    bb_entries = []
    for chain, resseq, icode, resname in residue_order:
        atom_map = by_residue[(chain, resseq, icode, resname)]
        missing = [name for name in BB_COMPONENT_NAMES if name not in atom_map]
        if missing:
            raise ValueError(
                f"AA backbone residue {resname} {chain}{resseq}{icode} is missing atoms: {missing}"
            )

        idxs = [int(atom_map[name]) for name in BB_COMPONENT_NAMES]
        weights_raw = [float(AA_BACKBONE_ATOM_MASS[name]) for name in BB_COMPONENT_NAMES]
        weight_sum = float(sum(weights_raw))
        weights = [weight / weight_sum for weight in weights_raw]
        coords_row = []
        for atom_idx in idxs:
            atom = backbone_atoms[atom_idx]
            coords_row.append([float(atom["x"]), float(atom["y"]), float(atom["z"])])

        bb_entries.append(
            {
                "bb_resseq": int(resseq),
                "bb_chain": chain,
                "bb_icode": icode,
                "bb_atom_index": int(atom_map["CA"]),
                "atom_indices": idxs,
                "atom_mask": [1, 1, 1, 1],
                "weights": weights,
                "reference_atom_indices": idxs,
                "reference_atom_coords": coords_row,
                "bb_comment": (
                    f"AA backbone residue {resname} {chain}{resseq}{icode}; "
                    f"runtime N/CA/C/O idx={idxs}; index_space=stage_runtime; w={weights}"
                ),
            }
        )
    return bb_entries


def write_backbone_metadata_h5(
    path: Path,
    bb_entries,
    total_atoms,
    env_atom_indices,
    n_protein_atoms,
    sequence,
):
    with h5py.File(path, "w") as h5:
        inp = h5.create_group("input")
        inp.create_dataset(
            "sequence",
            data=np.asarray([np.bytes_(resname) for resname in sequence], dtype="S3"),
        )

        bb_grp = inp.create_group("hybrid_bb_map")
        bb_grp.attrs["atom_index_space"] = b"stage_runtime"
        bb_grp.attrs["runtime_index_space"] = b"stage_runtime"
        bb_grp.attrs["reference_index_space"] = b"stage_runtime"
        bb_grp.attrs["reference_index_offset"] = np.int32(0)
        bb_grp.attrs["reference_index_count"] = np.int32(n_protein_atoms)
        bb_grp.create_dataset(
            "bb_residue_index",
            data=np.array([b["bb_resseq"] for b in bb_entries], dtype=np.int32),
        )
        bb_grp.create_dataset(
            "bb_atom_index",
            data=np.array([b["bb_atom_index"] for b in bb_entries], dtype=np.int32),
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
        env_grp.create_dataset("env_atom_indices", data=np.array(env_atom_indices, dtype=np.int32))
        membership = np.full(total_atoms, -1, dtype=np.int32)
        membership[:n_protein_atoms] = 0
        env_grp.create_dataset("protein_membership", data=membership)


def write_summary(path: Path, payload: Dict):
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")



# -----------------------------------------------------------------------------
# Stage Conversion Helpers
# -----------------------------------------------------------------------------
def runtime_input_pdb_path(script_dir, pdb_id):
    """Resolve MARTINI PDB path with optional runtime override."""
    override = os.environ.get("UPSIDE_RUNTIME_PDB_FILE", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(script_dir, f"pdb/{pdb_id}.MARTINI.pdb")


def read_martini3_nonbond_params(itp_file):
    """
    Read MARTINI nonbonded parameters from the main .itp file
    Returns a dictionary mapping (type1, type2) tuples to (sigma, epsilon) values
    """
    martini_table = {}
    
    if not os.path.exists(itp_file):
        print(f"Error: MARTINI parameter file '{itp_file}' not found!")
        return martini_table
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    # Parse #define macros for dry MARTINI style nonbond_params entries
    macro_params = {}
    for raw in lines:
        line = raw.split(';', 1)[0].strip()
        if not line.startswith('#define'):
            continue
        parts = line.split()
        # Format: #define m_VI 0.47 2.700
        if len(parts) >= 4:
            macro_name = parts[1]
            try:
                sigma = float(parts[2])
                epsilon = float(parts[3])
                macro_params[macro_name] = (sigma, epsilon)
            except ValueError:
                continue

    # Find the nonbond_params section
    in_nonbond_params = False
    param_count = 0
    
    for i, line in enumerate(lines):
        line = line.split(';', 1)[0].strip()
        
        # Check for nonbond_params section start
        if line == '[ nonbond_params ]' or line == '[nonbond_params]':
            in_nonbond_params = True
            continue
        elif line.startswith('[') and line.endswith(']') and line != '[ nonbond_params ]' and line != '[nonbond_params]':
            if in_nonbond_params:
                break
            in_nonbond_params = False
            continue
        
        # Parse nonbond_params lines
        if in_nonbond_params and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                type1 = parts[0]
                type2 = parts[1]
                try:
                    _func = int(parts[2])
                except ValueError:
                    continue
                try:
                    sigma = float(parts[3])  # nm
                    epsilon = float(parts[4])  # kJ/mol
                except ValueError:
                    # Dry MARTINI often uses a macro token in column 4
                    macro = parts[3]
                    if macro not in macro_params:
                        continue
                    sigma, epsilon = macro_params[macro]

                # Store both orientations for easy lookup
                martini_table[(type1, type2)] = (sigma, epsilon)
                martini_table[(type2, type1)] = (sigma, epsilon)
                param_count += 1
            elif len(parts) >= 4:
                # Alternative dry MARTINI format: type1 type2 func macro
                type1 = parts[0]
                type2 = parts[1]
                try:
                    _func = int(parts[2])
                except ValueError:
                    continue
                macro = parts[3]
                if macro not in macro_params:
                    continue
                sigma, epsilon = macro_params[macro]
                martini_table[(type1, type2)] = (sigma, epsilon)
                martini_table[(type2, type1)] = (sigma, epsilon)
                param_count += 1
    
    print(f"Read {len(martini_table)//2} unique nonbonded parameter pairs")
    return martini_table

def parse_itp_file(itp_file, target_molecule=None, preprocessor_defines=None):
    """
    Universal ITP parser to read MARTINI topology files.
    Returns a dictionary with parsed sections: atoms, bonds, angles, dihedrals, etc.
    If target_molecule is specified, only returns data for that specific molecule type.
    """
    topology = {
        'atoms': [],
        'bonds': [], 
        'angles': [],
        'dihedrals': [],
        'position_restraints': [],
        'exclusions': [],
        'moleculetype': None,
        'molecules': {}  # Store multiple molecule types
    }
    
    if not os.path.exists(itp_file):
        print(f"Warning: ITP file {itp_file} not found")
        return topology
    
    current_section = None
    current_molecule = None
    current_mol_data = None
    macro_defs = {}
    if preprocessor_defines:
        for macro_name, macro_value in preprocessor_defines.items():
            if isinstance(macro_value, (list, tuple)):
                macro_defs[macro_name] = [float(x) for x in macro_value]
            else:
                macro_defs[macro_name] = [float(macro_value)]
    pp_stack = []
    current_active = True

    def parse_macro_value(tok, macro_index=None):
        # Try direct float
        try:
            return float(tok)
        except ValueError:
            pass
        # Try macro lookup
        values = macro_defs.get(tok)
        if values is None:
            raise ValueError(f"Unknown macro '{tok}'")
        if macro_index is None:
            if not values:
                raise ValueError(f"Macro '{tok}' has no values")
            return values[0]
        if macro_index < len(values):
            return values[macro_index]
        if values:
            return values[-1]
        raise ValueError(f"Macro '{tok}' has no values")
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.split(';', 1)[0].strip()
            
            # Skip empty lines and comments
            if not line:
                continue

            # Handle a minimal preprocessor subset so dry MARTINI #ifndef/#else blocks
            # resolve to a single active topology branch.
            if line.startswith('#ifdef') or line.startswith('#ifndef'):
                parts = line.split()
                macro_name = parts[1] if len(parts) >= 2 else ""
                is_defined = macro_name in macro_defs
                cond = is_defined if line.startswith('#ifdef') else (not is_defined)
                pp_stack.append((current_active, cond))
                current_active = current_active and cond
                continue
            if line.startswith('#else'):
                if pp_stack:
                    parent_active, cond = pp_stack[-1]
                    cond = not cond
                    pp_stack[-1] = (parent_active, cond)
                    current_active = parent_active and cond
                continue
            if line.startswith('#endif'):
                if pp_stack:
                    parent_active, _cond = pp_stack.pop()
                    current_active = parent_active
                continue
            if not current_active:
                continue

            # Parse and store macro definitions used in dry MARTINI lipids
            if line.startswith('#define'):
                parts = line.split()
                if len(parts) >= 3:
                    macro_name = parts[1]
                    macro_vals = []
                    for tok in parts[2:]:
                        try:
                            macro_vals.append(float(tok))
                        except ValueError:
                            break
                    if macro_vals:
                        macro_defs[macro_name] = macro_vals
                continue

            # Ignore preprocessor directives like #ifdef/#ifndef/#endif/#include
            if line.startswith('#'):
                continue
            
            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                section_name = line[1:-1].strip().lower()
                current_section = section_name
                continue
            
            # Parse based on current section
            if current_section == 'moleculetype':
                parts = line.split()
                if len(parts) >= 1:
                    current_molecule = parts[0]
                    if topology['moleculetype'] is None:
                        topology['moleculetype'] = current_molecule
                    
                    # Initialize data structure for this molecule
                    current_mol_data = {
                        'atoms': [],
                        'bonds': [],
                        'angles': [],
                        'dihedrals': [],
                        'position_restraints': [],
                        'exclusions': []
                    }
                    topology['molecules'][current_molecule] = current_mol_data
            
            elif current_section == 'atoms' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        atom_data = {
                            'id': int(parts[0]),
                            'type': parts[1],
                            'resnr': int(parts[2]),
                            'residue': parts[3],
                            'atom': parts[4],
                            'cgnr': int(parts[5]),
                            'charge': 0.0,
                            'mass': 0.0
                        }
                        # Parse charge and mass more carefully
                        # Charge is typically in column 6 (0-indexed), mass in column 7
                        if len(parts) >= 7:
                            try:
                                atom_data['charge'] = float(parts[6])
                            except ValueError:
                                pass
                        if len(parts) >= 8:
                            try:
                                atom_data['mass'] = float(parts[7])
                            except ValueError:
                                pass
                        current_mol_data['atoms'].append(atom_data)
                        topology['atoms'].append(atom_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'bonds' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        bond_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'func': int(parts[2]),
                            'r0': 0.0,
                            'k': 0.0
                        }
                        if len(parts) >= 5:
                            bond_data['r0'] = parse_macro_value(parts[3], macro_index=0)  # nm
                            bond_data['k'] = parse_macro_value(parts[4], macro_index=1)   # kJ/mol/nm²
                        elif len(parts) >= 4:
                            # Dry MARTINI alias form: i j func macro
                            bond_data['r0'] = parse_macro_value(parts[3], macro_index=0)
                            bond_data['k'] = parse_macro_value(parts[3], macro_index=1)
                        current_mol_data['bonds'].append(bond_data)
                        topology['bonds'].append(bond_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'angles' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        angle_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'func': int(parts[3]),
                            'theta0': 0.0,  # degrees
                            'force_k': 0.0  # kJ/mol/rad²
                        }
                        if len(parts) >= 6:
                            angle_data['theta0'] = parse_macro_value(parts[4], macro_index=0)
                            angle_data['force_k'] = parse_macro_value(parts[5], macro_index=1)
                        else:
                            # Dry MARTINI alias form: i j k func macro
                            angle_data['theta0'] = parse_macro_value(parts[4], macro_index=0)
                            angle_data['force_k'] = parse_macro_value(parts[4], macro_index=1)
                        current_mol_data['angles'].append(angle_data)
                        topology['angles'].append(angle_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'dihedrals' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        dihedral_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'j': int(parts[1]) - 1,
                            'k': int(parts[2]) - 1,
                            'l': int(parts[3]) - 1,
                            'func': int(parts[4]),
                            'phi0': 0.0,
                            'k': 0.0,
                            'mult': 1
                        }
                        if len(parts) >= 7:
                            dihedral_data['phi0'] = float(parts[5])  # degrees
                            dihedral_data['k'] = float(parts[6])     # kJ/mol
                        if len(parts) >= 8:
                            dihedral_data['mult'] = int(parts[7])
                        current_mol_data['dihedrals'].append(dihedral_data)
                        topology['dihedrals'].append(dihedral_data)
                    except (ValueError, IndexError):
                        continue

            elif current_section == 'position_restraints' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        restraint_data = {
                            'i': int(parts[0]) - 1,  # Convert to 0-indexed
                            'func': int(parts[1]),
                            'fx': parse_macro_value(parts[2]),
                            'fy': parse_macro_value(parts[3]),
                            'fz': parse_macro_value(parts[4]),
                        }
                        current_mol_data['position_restraints'].append(restraint_data)
                        topology['position_restraints'].append(restraint_data)
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == 'exclusions' and current_mol_data is not None:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # Convert to 0-indexed and add all pairs
                        atoms = [int(part)-1 for part in parts]
                        for i in range(len(atoms)):
                            for j in range(i+1, len(atoms)):
                                exclusion = (atoms[i], atoms[j])
                                current_mol_data['exclusions'].append(exclusion)
                                topology['exclusions'].append(exclusion)
                    except ValueError:
                        continue
    
    # If target_molecule is specified, return only that molecule's data
    if target_molecule and target_molecule in topology['molecules']:
        mol_data = topology['molecules'][target_molecule]
        return {
            'atoms': mol_data['atoms'],
            'bonds': mol_data['bonds'],
            'angles': mol_data['angles'],
            'dihedrals': mol_data['dihedrals'],
            'position_restraints': mol_data['position_restraints'],
            'exclusions': mol_data['exclusions'],
            'moleculetype': target_molecule,
            'molecules': {target_molecule: mol_data}
        }
    
    return topology

def read_martini_masses(ff_file):
    """Read atom type masses from MARTINI force field file"""
    masses = {}
    ff_file_path = Path(ff_file).expanduser()
    if not ff_file_path.is_absolute():
        ff_file_path = (WORKFLOW_DIR / ff_file_path).resolve()
    else:
        ff_file_path = ff_file_path.resolve()
    
    if not ff_file_path.exists():
        raise ValueError(f"FATAL ERROR: Force field file '{ff_file}' not found.\n"
                        f"  Full path: {ff_file_path}\n"
                        f"  This file is required for atom type masses.\n"
                        f"  Please ensure the MARTINI force field file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    in_atomtypes = False
    with ff_file_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[ atomtypes ]'):
                in_atomtypes = True
                continue
            elif line.startswith('['):
                in_atomtypes = False
                continue
            
            if in_atomtypes and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    atom_type = parts[0]
                    try:
                        mass = float(parts[1])
                        masses[atom_type] = mass
                    except ValueError:
                        continue
    
    return masses

def convert_stage(pdb_id=None, stage='minimization', run_dir=None):
    """
    Main preparation function with stage-specific parameterization
    
    Args:
        stage: Simulation stage ('minimization', 'npt_equil', 'npt_equil_reduced', 'npt_prod')
        run_dir: Optional run directory (default: outputs/martini_test)
    """

    # Get UPSIDE home directory
    upside_path = os.environ['UPSIDE_HOME']

    # Resolve PDB ID from explicit argument first, then legacy CLI fallback.
    if pdb_id is None:
        if len(sys.argv) > 1:
            pdb_id = sys.argv[1]
        else:
            raise ValueError("FATAL ERROR: No PDB ID provided.\n"
                            f"  Usage: python {sys.argv[0]} <pdb_id>\n"
                            f"  Example: python {sys.argv[0]} 1rkl\n"
                            f"  Aborting to prevent incorrect simulation results.")

    # Get run directory from function parameter or use default
    if run_dir is None:
        if len(sys.argv) > 2 and sys.argv[2] != '--stage':
            run_dir = sys.argv[2]
        else:
            run_dir = "outputs/martini_test"
    os.makedirs(run_dir, exist_ok=True)

    # Get stage from environment variable or use default
    stage = os.environ.get('UPSIDE_SIMULATION_STAGE', stage)
    print(f"Preparing for stage: {stage}")

    # Stage-specific parameterization
    stage_params = {
        'minimization': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
        },
        'npt_equil': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
        },
        'npt_equil_reduced': {
            'lj_soften': 1,
            'lj_alpha': 0.05,
            'coulomb_soften': 1,
            'slater_alpha': 0.5,
        },
        'npt_prod': {
            'lj_soften': 0,
            'lj_alpha': 0.0,
            'coulomb_soften': 0,
            'slater_alpha': 0.0,
        }
    }

    # Get parameters for current stage
    params = stage_params.get(stage, stage_params['npt_prod'])
    stage_lipidhead_fc = float(os.environ.get('UPSIDE_BILAYER_LIPIDHEAD_FC', '0'))


    # Configuration
    strict_from_martini_pdb = True
    include_protein = True
    
    print("=== AA-Backbone Protein-Lipid System Preparation ===")
    print(f"PDB ID: {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    workflow_dir = str(WORKFLOW_DIR)
    # Read dry MARTINI parameter files
    print("\n=== Reading Dry MARTINI Parameters ===")
    ff_dir = Path(
        os.environ.get('UPSIDE_MARTINI_FF_DIR', str(REPO_ROOT / "parameters" / "dryMARTINI"))
    ).expanduser()
    if not ff_dir.is_absolute():
        ff_dir = (REPO_ROOT / ff_dir).resolve()
    else:
        ff_dir = ff_dir.resolve()
    ff_path = str(ff_dir)

    if not os.path.isdir(ff_path):
        raise ValueError(f"FATAL ERROR: Force-field directory '{ff_path}' not found.\n"
                        "  Set UPSIDE_MARTINI_FF_DIR to a valid dry-MARTINI force-field directory.")

    ff_files = sorted(os.listdir(ff_path))
    print(f"Using force-field directory: {ff_path}")

    def pick_ff_file(name, required=True):
        path = os.path.join(ff_path, name)
        if os.path.exists(path):
            return path
        if required:
            raise ValueError(f"FATAL ERROR: Required dry MARTINI file '{name}' not found in '{ff_path}'.\n"
                            f"  Available: {ff_files}")
        return None

    # Read nonbonded parameters
    martini_param_file = pick_ff_file("dry_martini_v2.1.itp")
    martini_table = read_martini3_nonbond_params(martini_param_file)
    
    if not martini_table:
        raise ValueError(f"FATAL ERROR: Could not read MARTINI parameters from '{martini_param_file}'\n"
                        f"  This file is required for proper force field parameterization.\n"
                        f"  Please ensure the dry MARTINI parameter file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    print("=== AA-Backbone Protein-Lipid System Detected ===")
    print("System contains DOPC environment plus direct AA backbone atoms.")
    
    # For both protein and lipid systems, we need DOPC parameters
    # Parse DOPC topology from ITP file
    dopc_param_file = pick_ff_file("dry_martini_v2.1_lipids.itp")
    lipid_preproc_defs = {}
    if stage_lipidhead_fc > 0.0:
        lipid_preproc_defs['BILAYER_LIPIDHEAD_FC'] = stage_lipidhead_fc
    full_topology = parse_itp_file(dopc_param_file, preprocessor_defines=lipid_preproc_defs)
    
    # Try to find DOPC or similar molecule
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
        raise ValueError(f"FATAL ERROR: DOPC molecule not found in '{dopc_param_file}'.\n"
                        f"  Available molecules: {available_molecules}\n"
                        f"  Please ensure DOPC is defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    dopc_bead_types = [atom['type'] for atom in dopc_topology['atoms']]
    dopc_charges = [atom['charge'] for atom in dopc_topology['atoms']]
    # Create mapping from atom names to bead types and charges
    dopc_atom_to_type = {atom['atom']: atom['type'] for atom in dopc_topology['atoms']}
    dopc_atom_to_charge = {atom['atom']: atom['charge'] for atom in dopc_topology['atoms']}
    print(f"Read DOPC topology: {len(dopc_bead_types)} bead types from {dopc_param_file}")
    
    # Parse ion topologies from ITP file
    ion_param_file = pick_ff_file("dry_martini_v2.1_ions.itp", required=False)
    if ion_param_file:
        ion_topology = parse_itp_file(ion_param_file)
        # Extract NA and CL atoms specifically
        # In MARTINI ion ITP files, residue name is "ION" and atom name is "NA" or "CL"
        na_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'NA']
        cl_atoms = [atom for atom in ion_topology['atoms'] if atom['atom'].upper() == 'CL']
        
        # For standard ions, use the first occurrence (standard chloride is TQ5, not SQ5n which is for acetate)
        na_bead_types = [na_atoms[0]['type']] if na_atoms else []
        na_charges = [na_atoms[0]['charge']] if na_atoms else []
        cl_bead_types = [cl_atoms[0]['type']] if cl_atoms else []  # Use first CL (TQ5), not acetate CL (SQ5n)
        cl_charges = [cl_atoms[0]['charge']] if cl_atoms else []
        print(f"Ion topology loaded: NA={len(na_bead_types)} type(s), CL={len(cl_bead_types)} type(s)")
    else:
        na_bead_types, na_charges = [], []
        cl_bead_types, cl_charges = [], []
        print("Ion topology file not found in selected FF")
    
    # Parse water topology from ITP file
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
    
    # Read bead masses from force field file
    mass_file = martini_param_file
    martini_masses = read_martini_masses(mass_file)
    print(f"Read {len(martini_masses)} atom type masses from force field file")
    
    # Read DOPC bonds and angles from parsed topology (for both lipid and mixed systems)
    dopc_bonds = [(bond['i'], bond['j']) for bond in dopc_topology['bonds']]
    dopc_bond_lengths = [bond['r0'] for bond in dopc_topology['bonds']]  # nm
    dopc_bond_force_constants = [bond['k'] for bond in dopc_topology['bonds']]  # kJ/mol/nm²
    
    dopc_angles = [(angle['i'], angle['j'], angle['k']) for angle in dopc_topology['angles']]
    dopc_angle_equil_deg = [angle['theta0'] for angle in dopc_topology['angles']]  # degrees
    dopc_angle_force_constants = [angle['force_k'] for angle in dopc_topology['angles']]  # kJ/mol/rad²
    dopc_position_restraints = dopc_topology.get('position_restraints', [])
    
    print(f"Read DOPC connectivity: {len(dopc_bonds)} bonds, {len(dopc_angles)} angles")
    
    # Validate that required topology data was found
    if not dopc_bonds:
        raise ValueError(f"FATAL ERROR: No DOPC bonds found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC bonds are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    if not dopc_angles:
        raise ValueError(f"FATAL ERROR: No DOPC angles found in topology from '{dopc_param_file}'.\n"
                        f"  This indicates incomplete molecule definition.\n"
                        f"  Please ensure DOPC angles are properly defined in the phospholipid parameter file.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    # Unit conversions for mapping native MARTINI units into the active
    # simulation unit system. These must be provided explicitly rather than
    # baked into the generator.
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

    # Pressure conversion (for NPT simulations):
    # 1 atm = 101325 Pa = 101.325 kJ/m³
    # 1 atm = 1.01325e-28 kJ/Å³
    # 1 atm = (1.01325e-28 kJ/Å³) × (6.02214e23 mol⁻¹) / (2.914952774272 kJ/mol)
    # 1 atm = 0.000020933 E_up/Å³
    # 1 bar = 0.986923 atm = 0.000020659 E_up/Å³
    # For LAMMPS real units compatibility, use 1 bar = 0.000020659 E_up/Å³
    pressure_conversion_bar_to_eup = 0.000020659  # bar → E_up/Å³

    # Bonds: kJ/mol/nm² → E_up/Å²
    # = (kJ/mol → E_up) / (nm² → Å²)
    # = (1/energy_conversion) / (length_conversion²)
    bond_conversion = 1.0 / (energy_conversion * length_conversion ** 2)  # ≈ 0.003406

    # Angles: kJ/mol/deg² → E_up/deg²
    # = (kJ/mol → E_up) (degrees stay the same)
    angle_conversion = 1.0 / energy_conversion  # ≈ 0.343

    # Dihedrals: kJ/mol → E_up
    dihedral_conversion = 1.0 / energy_conversion  # ≈ 0.343

    print("\n=== MARTINI Unit Conversions ===")
    print(f"Bond lengths (nm -> Å): 0.40 nm -> 4.0 Å")
    print(f"Bond force constants (kJ/mol/nm² -> E_up/Å²): 7000.0 -> {7000.0 * bond_conversion:.3f}")
    print(f"Angle equilibrium (degrees): 108.0°")
    print(f"Angle force constants (kJ/mol/deg² -> E_up/deg²): 21.5 -> {21.5 * angle_conversion:.6f}")
    print(f"Dihedral force constants (kJ/mol -> E_up): 400.0 -> {400.0 * dihedral_conversion:.6f}")
    print(f"Pressure (bar -> E_up/Å³): 1 bar -> {pressure_conversion_bar_to_eup:.9f}")
    print(f"Energy conversion factor: {energy_conversion} (kJ/mol -> E_up)")
    print(f"Length conversion: 1 nm = {length_conversion} Å, so 1 nm² = {length_conversion**2} Å²")
    print(f"Bond conversion factor: {bond_conversion:.6f} (divide by energy_conv × length_conv²)")
    print(f"Angle conversion factor: {angle_conversion:.6f} (divide by energy_conv)")
    print(f"Dihedral conversion factor: {dihedral_conversion:.6f} (divide by energy_conv)")
    
    # Read PDB file
    if strict_from_martini_pdb or include_protein:
        input_pdb_file = runtime_input_pdb_path(workflow_dir, pdb_id)
        print(f"\nUsing MARTINI PDB as base structure: {input_pdb_file}")
    
    # Read PDB and populate arrays
    print(f"\n=== Reading PDB Structure ===")
    initial_positions = []
    atom_types = []
    charges = []
    residue_ids = []
    atom_names = []
    residue_names = []
    chain_ids = []
    seg_ids = []
    
    protein_bonds = []
    protein_angles = []
    protein_dihedrals = []
    protein_constraints = []
    protein_position_restraints = []
    protein_exclusions = []

    print(f"\n=== Protein Connectivity ===")
    print("Protein MARTINI topology is retired for this workflow.")
    print("Protein bonded geometry will be supplied by injected Upside backbone nodes.")
    
    
    protein_residue_names = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
        'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX'
    }
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
                
            # Parse PDB line
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
            
            
            
            # Extract coordinates
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            initial_positions.append([x, y, z])
            
            is_protein = residue_name in protein_residue_names
            
            # Map to MARTINI type based on context
            if is_protein:
                if atom_name not in AA_BACKBONE_ATOM_MASS:
                    raise ValueError(
                        f"FATAL ERROR: Unexpected protein atom '{atom_name}' in residue '{residue_name}'.\n"
                        "  The AA-backbone workflow expects runtime protein atoms to contain only N/CA/C/O.\n"
                        "  Regenerate the runtime structure with martini_prepare_system.py before stage conversion.\n"
                        "  Aborting to prevent incorrect simulation results."
                    )
                martini_type = AA_BACKBONE_PLACEHOLDER_TYPE
                charge = 0.0
            elif residue_name == 'DOPC' or residue_name == 'DOP':
                # For DOPC, use the topology from parameter file (for both lipid and mixed systems)
                if atom_name in dopc_atom_to_type:
                    martini_type = dopc_atom_to_type[atom_name]
                    charge = dopc_atom_to_charge[atom_name]
                else:
                    available_atom_names = sorted(dopc_atom_to_type.keys())
                    raise ValueError(f"FATAL ERROR: Unknown DOPC atom '{atom_name}' in residue '{residue_name}'.\n"
                                   f"  Available DOPC atom names: {available_atom_names}\n"
                                   f"  This indicates incomplete DOPC topology mapping.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
            elif residue_name == 'W':
                if not water_bead_types:
                    raise ValueError(f"FATAL ERROR: No water bead types found in topology.\n"
                                   f"  This indicates incomplete water parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = water_bead_types[0]
                charge = water_charges[0]
            elif residue_name == 'NA':
                if not na_bead_types:
                    raise ValueError(f"FATAL ERROR: No sodium bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = na_bead_types[0]
                charge = na_charges[0]
            elif residue_name == 'CL':
                if not cl_bead_types:
                    raise ValueError(f"FATAL ERROR: No chloride bead types found in topology.\n"
                                   f"  This indicates incomplete ion parameter file.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                martini_type = cl_bead_types[0]
                charge = cl_charges[0]
            else:
                raise ValueError(f"FATAL ERROR: Unknown residue type '{residue_name}' for atom '{atom_name}'.\n"
                               f"  Supported residue types: PROTEIN, DOPC, W, NA, CL\n"
                               f"  This indicates incomplete system definition.\n"
                               f"  Aborting to prevent incorrect simulation results.")
            
            atom_types.append(martini_type)
            charges.append(charge)
    
    # Convert to numpy arrays
    initial_positions = np.array(initial_positions, dtype=float)
    atom_types = np.array(atom_types)
    charges = np.array(charges, dtype=float)
    residue_ids = np.array(residue_ids, dtype=int)
    atom_names = np.array(atom_names)
    residue_names = np.array(residue_names)
    chain_ids = np.array(chain_ids)
    seg_ids = np.array(seg_ids)
    n_atoms = len(initial_positions)
    
    # Read box dimensions from CRYST1 record
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
            raise ValueError(f"FATAL ERROR: No CRYST1 record found in PDB file '{input_pdb_file}'.\n"
                           f"  Box dimensions are required for proper simulation setup.\n"
                           f"  Please ensure the PDB file contains a CRYST1 record with box dimensions.\n"
                           f"  Aborting to prevent incorrect simulation results.")
    
    # Print system parameters
    print(f"Box dimensions: X={x_len:.3f}, Y={y_len:.3f}, Z={z_len:.3f} Angstroms")
    print(f"Box volume: {x_len * y_len * z_len:.1f} Å³")
    print(f"Total atoms: {n_atoms}")
    
    # Group atoms into molecules by chain (for proteins) or residue ID (for other molecules)
    molecules = []
    current_mol_atoms = []
    current_mol_names = []
    current_mol_indices = []
    current_resid = None
    current_resname = None
    current_chain = None
    
    # Define protein residue names
    protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
    
    def pick_chain_id(idx):
        if seg_ids[idx]:
            return seg_ids[idx]
        return chain_ids[idx]

    for i, (resid, resname, atom_name) in enumerate(zip(residue_ids, residue_names, atom_names)):
        # Determine molecule type
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
            current_mol_names = [atom_name]
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
                raise ValueError(f"FATAL ERROR: Molecule {i} contains mixed residue types: {unique_resnames}\n"
                               f"  This indicates incorrect molecule grouping.\n"
                               f"  Molecule atoms: {atoms}\n"
                               f"  Residue names: {residue_names_in_mol}\n"
                               f"  Aborting to prevent incorrect simulation results.")
        else:
            pass
    
    # Count molecules by type, but group all protein residues together
    mol_counts = Counter()
    protein_residue_count = 0
    
    protein_residue_keys = set()
    for resid, resname, idx in zip(residue_ids, residue_names, range(n_atoms)):
        if resname in protein_residues:
            protein_residue_keys.add((pick_chain_id(idx), resid))

    for mol_type, _, _ in molecules:
        if mol_type == 'PROTEIN':
            mol_counts[mol_type] += 1
        else:
            mol_counts[mol_type] += 1
    
    protein_residue_count = len(protein_residue_keys)
    if protein_residue_count > 0:
        protein_chains = mol_counts.get('PROTEIN', 0)
        mol_counts['PROTEIN'] = f"{protein_chains} chain(s) ({protein_residue_count} residues)"
    
    dopc_count = mol_counts.get('DOPC', 0)
    water_count = mol_counts.get('W', 0)
    
    print(f"\n=== Molecule Summary ===")
    for moltype, count in mol_counts.items():
        print(f"{moltype}: {count} molecules")
    
    # Create bonds and angles
    print(f"\n=== Creating Connectivity ===")
    
    # Initialize lists for bonds and angles
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
    
    # Create DOPC bonds and angles (for both lipid and mixed systems)
    dopc_molecules = [mol for mol in molecules if mol[0] == 'DOPC']  # unified label
    
    for mol_idx, (_, atom_names_mol, atom_indices) in enumerate(dopc_molecules):
        name_to_idx = {name: idx for name, idx in zip(atom_names_mol, atom_indices)}
        
        # Create bonds for this lipid
        for i, (bond_idx1, bond_idx2) in enumerate(dopc_bonds):
            if bond_idx1 < len(dopc_bead_types) and bond_idx2 < len(dopc_bead_types):
                atom1_name = atom_names_mol[bond_idx1]
                atom2_name = atom_names_mol[bond_idx2]
                atom1 = name_to_idx[atom1_name]
                atom2 = name_to_idx[atom2_name]
                bonds_list.append([atom1, atom2])
                bond_lengths_list.append(dopc_bond_lengths[i] * 10.0)  # nm to Å
                bond_force_constants_list.append(dopc_bond_force_constants[i] * bond_conversion)  # kJ/mol/nm² to E_up/Å²
        
        # Create angles for this lipid
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
                angle_force_constants_list.append(dopc_angle_force_constants[i] * angle_conversion)  # kJ/mol/deg² to E_up/deg²

        # Apply stage-specific lipid head-group restraints from dry MARTINI topology
        # (e.g., BILAYER_LIPIDHEAD_FC=200/100/50/20/10 for stages 6.2-6.6).
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
    
    print(f"Created {len(bonds_list)} bonds for {dopc_count} DOPC lipids")
    print(f"Created {len(angles_list)} angles for {dopc_count} DOPC lipids")
    print(f"Created {len(lipid_restraint_indices)} lipid position restraints (BILAYER_LIPIDHEAD_FC={stage_lipidhead_fc:g})")
    
    # Create protein connectivity if available
    protein_bond_count = 0
    protein_angle_count = 0
    protein_dihedral_count = 0
    protein_constraint_count = 0
    
    if protein_bonds or protein_constraints:
        print(f"\n=== Protein Connectivity from Upside Backbone Nodes ===")
        print("No MARTINI protein bonds, angles, or dihedrals are created in stage conversion.")
        
        # Add protein bonds to the bond list
        for i, j, r0_nm, k_kj in protein_bonds:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list (assuming protein atoms come first in the system)
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_bond_count += 1
        
        # Add protein constraints as bonds with large spring constants
        for i, j, r0_nm, k_kj in protein_constraints:
            # Convert MARTINI units to UPSIDE units
            r0_angstrom = r0_nm * 10.0  # nm to Å
            k_upside = k_kj * bond_conversion  # kJ/mol/nm² to E_up/Å²
            
            # Add to bond list
            bonds_list.append([i, j])
            bond_lengths_list.append(r0_angstrom)
            bond_force_constants_list.append(k_upside)
            protein_constraint_count += 1
        
        # Add protein angles to the angle list
        for i, j, k, theta0_deg, k_kj in protein_angles:
            # Convert MARTINI units to UPSIDE units
            theta0_upside = theta0_deg  # degrees (same unit)
            k_upside = k_kj * angle_conversion  # kJ/mol/deg² to E_up/deg²
            
            # Add to angle list
            angles_list.append([i, j, k])
            angle_equil_deg_list.append(theta0_upside)
            angle_force_constants_list.append(k_upside)
            protein_angle_count += 1
        
        # Add protein dihedrals to the dihedral list
        for i, j, k, l, phi0_deg, k_kj, func_type in protein_dihedrals:
            # Convert MARTINI units to UPSIDE units
            phi0_upside = phi0_deg  # degrees (same unit)
            k_upside = k_kj * dihedral_conversion  # kJ/mol to E_up
            
            # Add to dihedral list
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
        print("\nNo MARTINI protein connectivity is expected in the AA-backbone workflow")
    
    print(f"Total system bonds: {len(bonds_list)}")
    print(f"Total system angles: {len(angles_list)}")
    print(f"Total system dihedrals: {len(dihedrals_list)}")
    
    # Center and wrap positions
    print(f"\n=== Preparing Final Structure ===")
    center = np.mean(initial_positions, axis=0)
    centered_positions = initial_positions - center
    half_box = np.array([x_len/2, y_len/2, z_len/2])
    centered_positions = (centered_positions + half_box) % (2*half_box) - half_box
    final_positions = centered_positions
    
    # Create UPSIDE input file
    print(f"\n=== Creating UPSIDE Input File ===")
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
        mass = np.zeros(n_atoms, dtype='f4')
        for i, atom_type in enumerate(atom_types):
            if atom_type in NONINTERACTING_ATOM_TYPES:
                atom_name = atom_names[i].decode('utf-8') if isinstance(atom_names[i], bytes) else str(atom_names[i])
                if atom_name not in AA_BACKBONE_ATOM_MASS:
                    raise ValueError(
                        f"FATAL ERROR: No AA backbone mass is defined for atom '{atom_name}' (atom index {i})."
                    )
                mass[i] = AA_BACKBONE_ATOM_MASS[atom_name] / 12.0
                continue
            if atom_type not in martini_masses:
                raise ValueError(f"FATAL ERROR: Mass not found for atom type '{atom_type}' (atom index {i}).\n"
                                f"  Available atom types with masses: {sorted(martini_masses.keys())}\n"
                                f"  This indicates incomplete force field parameters.\n"
                                f"  Aborting to prevent incorrect simulation results.")
            mass[i] = martini_masses[atom_type] / 12.0
        
        mass_array = t.create_array(input_grp, 'mass', obj=mass)
        mass_array._v_attrs.arguments = np.array([b'mass'])
        mass_array._v_attrs.shape = mass.shape
        mass_array._v_attrs.n_atoms = n_atoms
        mass_array._v_attrs.initialized = True
        
        # Protein should NOT be held rigid during minimization
        # Allow the protein to relax and minimize its energy
        print("Protein atoms are free to move during minimization (no rigid constraints)")
        
        # Create stage-specific parameters group (always create this)
        stage_grp = t.create_group(input_grp, 'stage_parameters')
        stage_grp._v_attrs.enable = 1
        stage_grp._v_attrs.current_stage = b'minimization'
        print(f"Stage parameters initialized: current_stage={stage_grp._v_attrs.current_stage.decode()}")
        
        # ===================== NPT BAROSTAT CONFIGURATION =====================
        # Create barostat configuration group for NPT simulations
        # Settings are read from environment variables (set by run_sim_bilayer.sh)
        barostat_enable = int(os.environ.get('UPSIDE_NPT_ENABLE', '0'))
        if barostat_enable:
            print(f"\n=== Creating NPT Barostat Configuration ===")
            barostat_grp = t.create_group(input_grp, 'barostat')
            barostat_grp._v_attrs.enable = barostat_enable
            # Default: 1 bar = 0.000020659 E_up/Angstrom^3 (from 1 atm = 0.000020933215)
            barostat_grp._v_attrs.target_p_xy = float(os.environ.get('UPSIDE_NPT_TARGET_PXY', '0.000020659'))
            barostat_grp._v_attrs.target_p_z = float(os.environ.get('UPSIDE_NPT_TARGET_PZ', '0.000020659'))
            barostat_grp._v_attrs.tau_p = float(os.environ.get('UPSIDE_NPT_TAU', '1.0'))
            legacy_compressibility = float(os.environ.get('UPSIDE_NPT_COMPRESSIBILITY', '14.521180763676'))
            # Axis-specific compressibility for semi-isotropic membrane coupling
            barostat_grp._v_attrs.compressibility_xy = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_XY', str(legacy_compressibility))
            )
            barostat_grp._v_attrs.compressibility_z = float(
                os.environ.get('UPSIDE_NPT_COMPRESSIBILITY_Z', str(legacy_compressibility))
            )
            barostat_grp._v_attrs.interval = int(os.environ.get('UPSIDE_NPT_INTERVAL', '10'))
            barostat_grp._v_attrs.semi_isotropic = int(os.environ.get('UPSIDE_NPT_SEMI', '1'))
            barostat_grp._v_attrs.debug = int(os.environ.get('UPSIDE_NPT_DEBUG', '1'))
            print(f"  Enabled: {barostat_enable}")
            print("  Type: Berendsen")
            print(f"  Target Pxy: {barostat_grp._v_attrs.target_p_xy} E_up/Angstrom^3 (~1 bar)")
            print(f"  Target Pz: {barostat_grp._v_attrs.target_p_z} E_up/Angstrom^3 (~1 bar)")
            print(f"  Compressibility XY: {barostat_grp._v_attrs.compressibility_xy} Angstrom^3/E_up")
            print(f"  Compressibility Z: {barostat_grp._v_attrs.compressibility_z} Angstrom^3/E_up")
            print(f"  Tau_p: {barostat_grp._v_attrs.tau_p}")
            print(f"  Interval: {barostat_grp._v_attrs.interval} steps")
        else:
            print(f"\n=== NPT Barostat Disabled (NVT mode) ===")
        
        # Create type array
        type_array = t.create_array(input_grp, 'type', obj=atom_types.astype('S4'))
        type_array._v_attrs.arguments = np.array([b'type'])
        type_array._v_attrs.shape = atom_types.shape
        type_array._v_attrs.n_atoms = n_atoms
        type_array._v_attrs.initialized = True
        type_array._v_attrs.description = b"Interaction matrix type names (e.g., C1, C2, Qa, Qd)"

        # Store particle class per atom (protein/lipid/water/ion/other)
        particle_class = np.empty(n_atoms, dtype='S10')
        for i, resname in enumerate(residue_names):
            if resname in protein_residues:
                particle_class[i] = b"PROTEINAA"
            elif resname == 'DOP':
                particle_class[i] = b"LIPID"
            elif resname == 'W':
                particle_class[i] = b"WATER"
            elif resname in ('NA', 'CL'):
                particle_class[i] = b"ION"
            else:
                particle_class[i] = b"OTHER"
        particle_class_array = t.create_array(input_grp, 'particle_class', obj=particle_class)
        particle_class_array._v_attrs.arguments = np.array([b'particle_class'])
        particle_class_array._v_attrs.shape = particle_class.shape
        particle_class_array._v_attrs.n_atoms = n_atoms
        particle_class_array._v_attrs.initialized = True
        particle_class_array._v_attrs.description = b"Per-atom class: PROTEIN, LIPID, WATER, ION, OTHER"
        
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
        
        # Create potential group (required by UPSIDE)
        potential_grp = t.create_group(input_grp, 'potential')
        
        # Create MARTINI potential with proper parameters
        martini_potential = t.create_group(potential_grp, 'martini_potential')
        martini_potential._v_attrs.arguments = np.array([b'pos'])
        martini_potential._v_attrs.potential_type = b'lj_coulomb'
        # Note: epsilon and sigma attributes are required by C++ interface but not used in computation
        # They will be calculated from the coefficients array after it's created
        martini_potential._v_attrs.lj_cutoff = 12.0
        martini_potential._v_attrs.coul_cutoff = 12.0
        martini_potential._v_attrs.energy_conversion_kj_per_eup = energy_conversion
        martini_potential._v_attrs.length_conversion_angstrom_per_nm = length_conversion
        martini_potential._v_attrs.coulomb_constant_native_kj_mol_nm_e2 = coulomb_constant_native
        martini_potential._v_attrs.n_types = 1
        martini_potential._v_attrs.n_params = 4
        martini_potential._v_attrs.cutoff = 12.0
        martini_potential._v_attrs.cache_buffer = 1.0
        martini_potential._v_attrs.initialized = True
        force_cap = float(os.environ.get('UPSIDE_FORCE_CAP', '0'))
        martini_potential._v_attrs.force_cap = force_cap
        
        # PME configuration for long-range Coulomb interactions
        use_pme = int(os.environ.get('UPSIDE_USE_PME', '1'))
        pme_alpha = float(os.environ.get('UPSIDE_PME_ALPHA', '0.01'))
        pme_rcut = float(os.environ.get('UPSIDE_PME_RCUT', '10.0'))
        pme_nx = int(os.environ.get('UPSIDE_PME_NX', '32'))
        pme_ny = int(os.environ.get('UPSIDE_PME_NY', '32'))
        pme_nz = int(os.environ.get('UPSIDE_PME_NZ', '32'))
        pme_order = int(os.environ.get('UPSIDE_PME_ORDER', '4'))
        
        martini_potential._v_attrs.use_pme = use_pme
        martini_potential._v_attrs.pme_alpha = pme_alpha
        martini_potential._v_attrs.pme_rcut = pme_rcut
        martini_potential._v_attrs.pme_nx = pme_nx
        martini_potential._v_attrs.pme_ny = pme_ny
        martini_potential._v_attrs.pme_nz = pme_nz
        martini_potential._v_attrs.pme_order = pme_order

        martini_potential._v_attrs.x_len = x_len
        martini_potential._v_attrs.y_len = y_len
        martini_potential._v_attrs.z_len = z_len
        # Optional softening controls via environment variables
        # UPSIDE_SOFTEN_COULOMB: 1 to enable Slater softening for Coulomb
        # UPSIDE_SLATER_ALPHA: float value for Slater alpha (1/Angstrom)
        # UPSIDE_SOFTEN_LJ: 1 to enable soft-core LJ
        # UPSIDE_LJ_ALPHA: float value for LJ softening alpha (dimensionless)
        soften_coul = int(os.environ.get('UPSIDE_SOFTEN_COULOMB', '0'))
        slater_alpha = float(os.environ.get('UPSIDE_SLATER_ALPHA', '1.0'))
        soften_lj = int(os.environ.get('UPSIDE_SOFTEN_LJ', '0'))
        lj_alpha = float(os.environ.get('UPSIDE_LJ_ALPHA', '1.0'))
        
        # PME configuration via environment variables
        # UPSIDE_USE_PME: 1 to enable Particle Mesh Ewald for long-range Coulomb
        # UPSIDE_PME_ALPHA: PME screening parameter (default: 0.2)
        # UPSIDE_PME_RCUT: Real space cutoff in Angstroms (default: 10.0)
        # UPSIDE_PME_NX/NY/NZ: Grid dimensions (default: 32, should be powers of 2)
        # UPSIDE_PME_ORDER: B-spline interpolation order (default: 4)
        use_pme = int(os.environ.get('UPSIDE_USE_PME', '0'))
        pme_alpha = float(os.environ.get('UPSIDE_PME_ALPHA', '0.01'))
        pme_rcut = float(os.environ.get('UPSIDE_PME_RCUT', '10.0'))
        pme_nx = int(os.environ.get('UPSIDE_PME_NX', '32'))
        pme_ny = int(os.environ.get('UPSIDE_PME_NY', '32'))
        pme_nz = int(os.environ.get('UPSIDE_PME_NZ', '32'))
        pme_order = int(os.environ.get('UPSIDE_PME_ORDER', '4'))

        # Set stage-specific softening parameters
        martini_potential._v_attrs.coulomb_soften = params['coulomb_soften']
        if params['coulomb_soften']:
            martini_potential._v_attrs.slater_alpha = params['slater_alpha']
        martini_potential._v_attrs.lj_soften = params['lj_soften']
        if params['lj_soften']:
            martini_potential._v_attrs.lj_soften_alpha = params['lj_alpha']
        
        # PME configuration
        martini_potential._v_attrs.use_pme = use_pme
        if use_pme:
            martini_potential._v_attrs.pme_alpha = pme_alpha
            martini_potential._v_attrs.pme_rcut = pme_rcut
            print(f"PME enabled: alpha={pme_alpha}, rcut={pme_rcut}, grid={pme_nx}x{pme_ny}x{pme_nz}, order={pme_order}")
        else:
            print("PME disabled: using standard Coulomb cutoff")

        # Ewald summation configuration via environment variables
        # UPSIDE_EWALD_ENABLE: 1 to enable Ewald summation for long-range Coulomb
        # UPSIDE_EWALD_ALPHA: Ewald screening parameter in 1/Angstrom (default: 0.2)
        # UPSIDE_EWALD_KMAX: k-space cutoff (default: 5)
        # UPSIDE_EWALD_USE_CARDINAL_BSPLINE: 1 to approximate trig via periodic cardinal cubic B-spline
        # UPSIDE_EWALD_BSPLINE_GRID: lookup grid size for the periodic trig table
        ewald_enabled = int(os.environ.get('UPSIDE_EWALD_ENABLE', '0'))
        ewald_alpha = float(os.environ.get('UPSIDE_EWALD_ALPHA', '0.2'))
        ewald_kmax = int(os.environ.get('UPSIDE_EWALD_KMAX', '5'))
        ewald_use_cardinal_bspline = int(os.environ.get('UPSIDE_EWALD_USE_CARDINAL_BSPLINE', '1'))
        ewald_bspline_grid = int(os.environ.get('UPSIDE_EWALD_BSPLINE_GRID', '16384'))

        martini_potential._v_attrs.ewald_enabled = ewald_enabled
        if ewald_enabled:
            martini_potential._v_attrs.ewald_alpha = ewald_alpha
            martini_potential._v_attrs.ewald_kmax = ewald_kmax
            martini_potential._v_attrs.ewald_use_cardinal_bspline = ewald_use_cardinal_bspline
            martini_potential._v_attrs.ewald_bspline_grid = ewald_bspline_grid
            print(f"Ewald summation enabled: alpha={ewald_alpha} A^-1, kmax={ewald_kmax}, cardinal_bspline={ewald_use_cardinal_bspline}, grid={ewald_bspline_grid}")
        else:
            print("Ewald summation disabled")

        # Periodic boundary potential removed - using NVT ensemble without boundaries
        
        # PME node creation removed - using Coulomb spline tables instead
        
        # NPT barostat configuration removed - using NVT ensemble without boundaries

        # Create atom indices and charges arrays for the potential
        t.create_array(martini_potential, 'atom_indices', obj=np.arange(n_atoms))
        t.create_array(martini_potential, 'charges', obj=charges)
        
        # Create pairs and coefficients for non-bonded interactions with proper exclusions
        pairs_list = []
        coeff_array = []
        
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
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                pair = (i, j)
                
                # Skip 1-2 pairs (full exclusion) - nrexcl=1
                if pair in bonded_pairs_12:
                    excluded_12_count += 1
                    continue
                
                # Skip additional exclusions from ITP file
                if pair in additional_exclusions:
                    excluded_additional_count += 1
                    continue
                
                # No scaling needed for MARTINI (nrexcl=1)
                scale_factor = 1.0
                
                pairs_list.append([i, j])
                
                # Get bead types for this pair
                type1 = atom_types[i].decode('utf-8') if isinstance(atom_types[i], bytes) else str(atom_types[i])
                type2 = atom_types[j].decode('utf-8') if isinstance(atom_types[j], bytes) else str(atom_types[j])
                
                if type1 in NONINTERACTING_ATOM_TYPES or type2 in NONINTERACTING_ATOM_TYPES:
                    coeff_array.append([0.0, 0.0, 0.0, 0.0])
                    continue

                # Look up MARTINI parameters for this bead type pair
                if (type1, type2) in martini_table:
                    sigma_nm, epsilon_kj = martini_table[(type1, type2)]
                elif (type2, type1) in martini_table:
                    sigma_nm, epsilon_kj = martini_table[(type2, type1)]
                else:
                    # Raise error for missing interaction parameters
                    available_types = sorted(set([t[0] for t in martini_table.keys()] + [t[1] for t in martini_table.keys()]))
                    raise ValueError(f"FATAL ERROR: Missing interaction parameters for bead type pair ({type1}, {type2})\n"
                                   f"  Atom indices: {i} ({type1}) - {j} ({type2})\n"
                                   f"  This indicates incomplete MARTINI force field parameters.\n"
                                   f"  Available bead types in parameter table: {available_types}\n"
                                   f"  Aborting to prevent incorrect simulation results.")
                
                # Convert to UPSIDE units
                epsilon = epsilon_kj / energy_conversion  # kJ/mol → E_up
                sigma = sigma_nm * length_conversion  # nm → Å
                q1 = charges[i] * scale_factor
                q2 = charges[j] * scale_factor
                coeff_array.append([epsilon * scale_factor, sigma, q1, q2])
        
        print(f"Excluded {excluded_12_count} 1-2 bonded pairs from non-bonded interactions (nrexcl=1)")
        print(f"Excluded {excluded_additional_count} additional pairs from ITP exclusions")

        t.create_array(martini_potential, 'pairs', obj=np.array(pairs_list, dtype=int))
        t.create_array(martini_potential, 'coefficients', obj=np.array(coeff_array, dtype='f4'))
        
        # Calculate representative epsilon and sigma from the coefficients array
        # These are required by the C++ interface but not used in computation
        if coeff_array:
            # Use the most common epsilon and sigma values from the coefficients
            epsilon_values = [coeff[0] for coeff in coeff_array]  # epsilon values
            sigma_values = [coeff[1] for coeff in coeff_array]    # sigma values
            
            # Use the median values as representative
            epsilon_values.sort()
            sigma_values.sort()
            median_epsilon = epsilon_values[len(epsilon_values)//2]
            median_sigma = sigma_values[len(sigma_values)//2]
            
            martini_potential._v_attrs.epsilon = median_epsilon
            martini_potential._v_attrs.sigma = median_sigma
        else:
            raise ValueError("FATAL ERROR: No interaction coefficients found.\n"
                           f"  This indicates no non-bonded interactions were defined.\n"
                           f"  Aborting to prevent incorrect simulation results.")

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
            # Some builds of UPSIDE expect 'equil_dist' for this potential; write both names for compatibility
            eq_deg = np.array(dihedral_equil_deg_list, dtype='f4')
            t.create_array(dihedral_group, 'equil_angle_deg', obj=eq_deg)
            t.create_array(dihedral_group, 'equil_dist', obj=eq_deg)
            t.create_array(dihedral_group, 'spring_const', obj=np.array(dihedral_force_constants_list, dtype='f4'))
            # Store dihedral type information (1=periodic, 2=harmonic)
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
            # Backward-compatible scalar spring constant for older readers.
            t.create_array(restraint_group, 'spring_const', obj=np.max(spring_xyz, axis=1).astype('f4'))
    
    print(f"Created UPSIDE input file: {input_file}")
    print(f"Preparation complete!")
    
    # Save preparation summary
    summary_file = f"{run_dir}/preparation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=== MARTINI 3.0 Preparation Summary ===\n")
        f.write(f"PDB ID: {pdb_id}\n")
        f.write(f"Total atoms: {n_atoms}\n")
        f.write(f"Box dimensions: {x_len:.3f} x {y_len:.3f} x {z_len:.3f} Å\n")
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

def require_h5import():
    exe = shutil.which("h5import")
    if not exe:
        raise SystemExit("ERROR: h5import was not found in PATH")
    return exe


def read_sc_table(path: Path):
    if not path.exists():
        raise SystemExit(f"ERROR: SC table JSON not found: {path}")
    return json.loads(path.read_text())


def normalize_string(value):
    return "" if value is None else str(value)


def build_factorized_sc_table(sc_table):
    tables_by_residue = sc_table.get("tables_by_residue")
    if not isinstance(tables_by_residue, dict) or not tables_by_residue:
        raise SystemExit("ERROR: sc_table.json is missing non-empty tables_by_residue")

    residues = list(tables_by_residue)
    first_residue = residues[0]
    targets = list(tables_by_residue[first_residue])
    if not targets:
        raise SystemExit("ERROR: first residue has no target tables")

    reference_grid = np.asarray(
        tables_by_residue[first_residue][targets[0]]["grid_nm"], dtype=np.float32
    )
    if reference_grid.ndim != 1 or reference_grid.size < 2:
        raise SystemExit("ERROR: invalid radial grid in sc_table.json")
    reference_cos_grid = np.asarray(
        tables_by_residue[first_residue][targets[0]]["cos_theta_grid"], dtype=np.float32
    )
    if reference_cos_grid.ndim != 1 or reference_cos_grid.size < 2:
        raise SystemExit("ERROR: invalid cos(theta) grid in sc_table.json")

    max_rotamer = 0
    for residue in residues:
        for target in targets:
            entry = tables_by_residue[residue][target]
            max_rotamer = max(max_rotamer, int(entry.get("rotamer_count", 0)))
    if max_rotamer < 1:
        raise SystemExit("ERROR: no rotamer-resolved SC tables found in sc_table.json")

    radial_energy = np.zeros((len(residues), len(targets), reference_grid.size), dtype=np.float32)
    angular_energy = np.zeros((len(residues), len(targets), reference_grid.size), dtype=np.float32)
    angular_profile = np.zeros(
        (len(residues), len(targets), reference_cos_grid.size), dtype=np.float32
    )
    rotamer_count = np.zeros((len(residues),), dtype=np.float32)
    rotamer_probability_fixed = np.zeros((len(residues), max_rotamer), dtype=np.float32)
    rotamer_radial_energy = np.zeros(
        (len(residues), max_rotamer, len(targets), reference_grid.size), dtype=np.float32
    )
    rotamer_angular_energy = np.zeros(
        (len(residues), max_rotamer, len(targets), reference_grid.size), dtype=np.float32
    )
    rotamer_angular_profile = np.zeros(
        (len(residues), max_rotamer, len(targets), reference_cos_grid.size), dtype=np.float32
    )

    for ri, residue in enumerate(residues):
        residue_tables = tables_by_residue[residue]
        missing_targets = sorted(set(targets) - set(residue_tables))
        extra_targets = sorted(set(residue_tables) - set(targets))
        if missing_targets or extra_targets:
            raise SystemExit(
                f"ERROR: target mismatch for residue {residue}: "
                f"missing={missing_targets} extra={extra_targets}"
            )

        for ti, target in enumerate(targets):
            entry = residue_tables[target]
            entry_rot_count = int(entry.get("rotamer_count", 0))
            if entry_rot_count < 1 or entry_rot_count > max_rotamer:
                raise SystemExit(
                    f"ERROR: invalid rotamer_count for residue {residue} target {target}: "
                    f"{entry_rot_count}"
                )
            if ti == 0:
                rotamer_count[ri] = float(entry_rot_count)
                prob_row = np.asarray(entry.get("rotamer_probability_fixed", []), dtype=np.float32)
                if prob_row.shape != (entry_rot_count,):
                    raise SystemExit(
                        f"ERROR: rotamer_probability_fixed shape mismatch for residue {residue}: "
                        f"{prob_row.shape} vs {(entry_rot_count,)}"
                    )
                rotamer_probability_fixed[ri, :entry_rot_count] = prob_row
            elif int(rotamer_count[ri]) != entry_rot_count:
                raise SystemExit(
                    f"ERROR: inconsistent rotamer_count for residue {residue}: "
                    f"{int(rotamer_count[ri])} vs {entry_rot_count}"
                )

            grid = np.asarray(entry["grid_nm"], dtype=np.float32)
            if grid.shape != reference_grid.shape or not np.allclose(grid, reference_grid):
                raise SystemExit(
                    f"ERROR: non-uniform radial grid for residue {residue} target {target}"
                )
            cos_grid = np.asarray(entry["cos_theta_grid"], dtype=np.float32)
            if cos_grid.shape != reference_cos_grid.shape or not np.allclose(cos_grid, reference_cos_grid):
                raise SystemExit(
                    f"ERROR: non-uniform cos(theta) grid for residue {residue} target {target}"
                )

            radial_row = np.asarray(entry["radial_energy_kj_mol"], dtype=np.float32)
            angular_row = np.asarray(entry["angular_energy_kj_mol"], dtype=np.float32)
            profile_row = np.asarray(entry["angular_profile"], dtype=np.float32)
            rot_radial = np.asarray(entry["rotamer_radial_energy_kj_mol"], dtype=np.float32)
            rot_angular = np.asarray(entry["rotamer_angular_energy_kj_mol"], dtype=np.float32)
            rot_profile = np.asarray(entry["rotamer_angular_profile"], dtype=np.float32)

            if radial_row.shape != reference_grid.shape:
                raise SystemExit(
                    f"ERROR: radial_energy_kj_mol shape mismatch for residue {residue} "
                    f"target {target}: {radial_row.shape} vs {reference_grid.shape}"
                )
            if angular_row.shape != reference_grid.shape:
                raise SystemExit(
                    f"ERROR: angular_energy_kj_mol shape mismatch for residue {residue} "
                    f"target {target}: {angular_row.shape} vs {reference_grid.shape}"
                )
            if profile_row.shape != reference_cos_grid.shape:
                raise SystemExit(
                    f"ERROR: angular_profile shape mismatch for residue {residue} target {target}: "
                    f"{profile_row.shape} vs {reference_cos_grid.shape}"
                )
            if rot_radial.shape != (entry_rot_count, reference_grid.size):
                raise SystemExit(
                    f"ERROR: rotamer_radial_energy_kj_mol shape mismatch for residue {residue} "
                    f"target {target}: {rot_radial.shape} vs {(entry_rot_count, reference_grid.size)}"
                )
            if rot_angular.shape != (entry_rot_count, reference_grid.size):
                raise SystemExit(
                    f"ERROR: rotamer_angular_energy_kj_mol shape mismatch for residue {residue} "
                    f"target {target}: {rot_angular.shape} vs {(entry_rot_count, reference_grid.size)}"
                )
            if rot_profile.shape != (entry_rot_count, reference_cos_grid.size):
                raise SystemExit(
                    f"ERROR: rotamer_angular_profile shape mismatch for residue {residue} "
                    f"target {target}: {rot_profile.shape} vs {(entry_rot_count, reference_cos_grid.size)}"
                )

            radial_energy[ri, ti, :] = radial_row
            angular_energy[ri, ti, :] = angular_row
            angular_profile[ri, ti, :] = profile_row
            rotamer_radial_energy[ri, :entry_rot_count, ti, :] = rot_radial
            rotamer_angular_energy[ri, :entry_rot_count, ti, :] = rot_angular
            rotamer_angular_profile[ri, :entry_rot_count, ti, :] = rot_profile

    return (
        residues,
        targets,
        reference_grid,
        reference_cos_grid,
        radial_energy,
        angular_energy,
        angular_profile,
        rotamer_count,
        rotamer_probability_fixed,
        rotamer_radial_energy,
        rotamer_angular_energy,
        rotamer_angular_profile,
    )


def write_text_h5_dataset(h5import_exe: str, tmpdir: Path, output_h5: Path, dataset_path: str, values):
    txt_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.txt"
    cfg_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.cfg"
    txt_path.write_text("".join(f"{normalize_string(v)}\n" for v in values))
    cfg_path.write_text(f"PATH {dataset_path}\nINPUT-CLASS STR\n")
    subprocess.run(
        [h5import_exe, str(txt_path), "-c", str(cfg_path), "-o", str(output_h5)],
        check=True,
    )


def write_float_h5_dataset(
    h5import_exe: str,
    tmpdir: Path,
    output_h5: Path,
    dataset_path: str,
    array: np.ndarray,
):
    arr = np.asarray(array, dtype=np.float32)
    bin_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.bin"
    cfg_path = tmpdir / f"{dataset_path.strip('/').replace('/', '__')}.cfg"
    arr.tofile(bin_path)
    dims = " ".join(str(x) for x in arr.shape)
    cfg_path.write_text(
        "\n".join(
            [
                f"PATH {dataset_path}",
                "INPUT-CLASS FP",
                "INPUT-SIZE 32",
                "INPUT-BYTE-ORDER LE",
                f"RANK {arr.ndim}",
                f"DIMENSION-SIZES {dims}",
                "OUTPUT-CLASS FP",
                "OUTPUT-SIZE 32",
                "OUTPUT-ARCHITECTURE NATIVE",
                "",
            ]
        )
    )
    subprocess.run(
        [h5import_exe, str(bin_path), "-c", str(cfg_path), "-o", str(output_h5)],
        check=True,
    )


def build_sc_martini_h5(sc_table_json: Path, output_h5: Path):
    h5import_exe = require_h5import()
    sc_table_path = Path(sc_table_json).expanduser().resolve()
    output_h5 = Path(output_h5).expanduser().resolve()
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    if output_h5.exists():
        output_h5.unlink()

    sc_table = read_sc_table(sc_table_path)
    (
        residues,
        targets,
        grid_nm,
        cos_theta_grid,
        radial_energy_kj_mol,
        angular_energy_kj_mol,
        angular_profile,
        rotamer_count,
        rotamer_probability_fixed,
        rotamer_radial_energy_kj_mol,
        rotamer_angular_energy_kj_mol,
        rotamer_angular_profile,
    ) = build_factorized_sc_table(sc_table)

    with tempfile.TemporaryDirectory(prefix="build_sc_martini_h5.") as tmp:
        tmpdir = Path(tmp)
        write_text_h5_dataset(h5import_exe, tmpdir, output_h5, "/schema", [sc_table.get("schema", "")])
        write_text_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/forcefield_name", [sc_table.get("forcefield_name", "")]
        )
        write_text_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/created_at_utc", [sc_table.get("created_at_utc", "")]
        )
        write_text_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/manifest_path", [sc_table.get("manifest_path", "")]
        )
        write_text_h5_dataset(h5import_exe, tmpdir, output_h5, "/restype_order", residues)
        write_text_h5_dataset(h5import_exe, tmpdir, output_h5, "/target_order", targets)
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/grid_nm", grid_nm)
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/cos_theta_grid", cos_theta_grid)
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/rotamer_count", rotamer_count)
        write_float_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/rotamer_probability_fixed", rotamer_probability_fixed
        )
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/radial_energy_kj_mol", radial_energy_kj_mol)
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/angular_energy_kj_mol", angular_energy_kj_mol)
        write_float_h5_dataset(h5import_exe, tmpdir, output_h5, "/angular_profile", angular_profile)
        write_float_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/rotamer_radial_energy_kj_mol", rotamer_radial_energy_kj_mol
        )
        write_float_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/rotamer_angular_energy_kj_mol", rotamer_angular_energy_kj_mol
        )
        write_float_h5_dataset(
            h5import_exe, tmpdir, output_h5, "/rotamer_angular_profile", rotamer_angular_profile
        )

    print(
        f"Built {output_h5} from {sc_table_path} with {len(residues)} residues, "
        f"{len(targets)} targets, {cos_theta_grid.size} angular points, {grid_nm.size} radial points"
    )


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
        bb = require_group(inp, "hybrid_bb_map")
        env = require_group(inp, "hybrid_env_topology")

        if "hybrid_control" in inp:
            ctrl = require_group(inp, "hybrid_control")
            for attr in [
                "enable",
                "activation_stage",
                "preprod_protein_mode",
                "exclude_intra_protein_martini",
                "schema_version",
            ]:
                if attr not in ctrl.attrs:
                    raise ValueError(f"Missing hybrid_control attr: {attr}")

            if "protein_env_interface_scale" in ctrl.attrs:
                interface_scale = float(ctrl.attrs["protein_env_interface_scale"])
                if not np.isfinite(interface_scale) or interface_scale <= 0.0:
                    raise ValueError(
                        "hybrid_control/protein_env_interface_scale must be finite and > 0"
                    )

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

        for i in range(n_bb):
            bb_i = int(bb_atom_idx[i])
            if bb_i < 0 or bb_i >= n_atom_runtime:
                raise ValueError(f"BB proxy index out of bounds at row {i}: {bb_i}")
            if membership[bb_i] < 0:
                raise ValueError(f"BB proxy index is not protein atom at row {i}: {bb_i}")
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


def refresh_hybrid_reference_carriers(h5f, pos):
    if "/input/hybrid_bb_map" not in h5f:
        return pos

    grp = h5f["/input/hybrid_bb_map"]
    required = ("bb_atom_index", "atom_indices", "weights", "reference_atom_coords")
    if not all((f"/input/hybrid_bb_map/{k}" in h5f) for k in required):
        return pos

    bb_idx = grp["bb_atom_index"][:].astype(np.int32)
    comp_idx = grp["atom_indices"][:].astype(np.int32)
    comp_w = grp["weights"][:].astype(np.float64)
    ref_xyz = grp["reference_atom_coords"][:].astype(np.float64)

    if comp_idx.ndim != 2 or comp_idx.shape[1] != 4:
        return pos
    if ref_xyz.shape != (comp_idx.shape[0], 4, 3):
        return pos

    n_atom = pos.shape[0]
    for k in range(comp_idx.shape[0]):
        bb = int(bb_idx[k]) if k < bb_idx.shape[0] else -1
        if bb < 0 or bb >= n_atom:
            continue

        valid = (comp_idx[k] >= 0) & (comp_idx[k] < n_atom) & (comp_w[k] > 0.0)
        if not np.any(valid):
            continue

        w = comp_w[k][valid]
        wsum = float(np.sum(w))
        if wsum <= 0.0:
            continue
        w = w / wsum

        ref_pts = ref_xyz[k][valid]
        ref_com = np.sum(ref_pts * w[:, None], axis=0)
        bb_pos = pos[bb, :, 0].astype(np.float64)
        shift = bb_pos - ref_com

        target_idx = comp_idx[k][valid]
        aligned = ref_pts + shift[None, :]
        for ai, pxyz in zip(target_idx, aligned):
            pos[int(ai), :, 0] = pxyz.astype(pos.dtype)

    return pos


def set_initial_position(input_file, output_file):
    strict_copy = env_bool("UPSIDE_SET_INITIAL_STRICT_COPY", False)
    apply_refresh_hybrid = env_bool(
        "UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS", default=(not strict_copy)
    )
    apply_recenter_production = env_bool(
        "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION", default=(not strict_copy)
    )

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

    with h5py.File(output_file, "r+") as f:
        target_pos = f["/input/pos"][:]
        target_n = target_pos.shape[0]
        source_n = last_pos.shape[0]
        if source_n != target_n:
            merged = target_pos.copy()
            merged[: min(source_n, target_n), :, :] = last_pos[: min(source_n, target_n), :, :]
            last_pos = merged

        if strict_copy and not apply_refresh_hybrid and not apply_recenter_production:
            print("Strict handoff mode: preserving exact coordinates from previous stage output.")

        if apply_refresh_hybrid:
            last_pos = refresh_hybrid_reference_carriers(f, last_pos)
        if apply_recenter_production:
            last_pos = recenter_protein_for_production(f, last_pos, last_box)

        if "/input/pos" in f:
            del f["/input/pos"]
        f.create_dataset("/input/pos", data=last_pos)

        if last_box is not None and "/input/potential/martini_potential" in f:
            pot_grp = f["/input/potential/martini_potential"]
            pot_grp.attrs["x_len"] = float(last_box[0])
            pot_grp.attrs["y_len"] = float(last_box[1])
            pot_grp.attrs["z_len"] = float(last_box[2])
            print(
                f"Updated box dimensions: x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}"
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


def resolve_sequence(inp, residue_count: int):
    if "sequence" in inp:
        sequence = [normalize_resname(x) for x in decode_string_array(inp["sequence"])]
        if len(sequence) == residue_count:
            return sequence
    raise ValueError(
        "Missing or inconsistent /input/sequence for AA-backbone sidechain injection: "
        f"expected {residue_count} residues"
    )


def build_affine_atoms(inp):
    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map in stage file")
    bb_grp = inp["hybrid_bb_map"]
    if "reference_atom_indices" not in bb_grp:
        raise ValueError("hybrid_bb_map/reference_atom_indices is required for stage-7 CB placement")
    ref_offset = int(bb_grp.attrs.get("reference_index_offset", -1))
    if ref_offset < 0:
        raise ValueError("hybrid_bb_map/reference_index_offset is missing or invalid")

    bb_residue_raw = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_ref_idx = bb_grp["reference_atom_indices"][:].astype(np.int32)
    residue_ids = unique_preserving_order(int(x) for x in bb_residue_raw.tolist())

    residue_to_ncac = {}
    n_atom = int(inp["pos"].shape[0])
    for resid, ref_row in zip(bb_residue_raw.tolist(), bb_ref_idx.tolist()):
        rid = int(resid)
        n_idx, ca_idx, c_idx = [int(ref_row[0]), int(ref_row[1]), int(ref_row[2])]
        if n_idx < 0 or ca_idx < 0 or c_idx < 0:
            continue
        n_rt = ref_offset + n_idx
        ca_rt = ref_offset + ca_idx
        c_rt = ref_offset + c_idx
        for idx in (n_rt, ca_rt, c_rt):
            if idx < 0 or idx >= n_atom:
                raise ValueError(
                    f"Backbone carrier index out of bounds for residue {rid}: idx={idx}, n_atom={n_atom}"
                )
        residue_to_ncac.setdefault(rid, (n_rt, ca_rt, c_rt))

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

    if not env_atom_index:
        raise ValueError("No non-protein dry-MARTINI particles matched martini.h5 target_order")

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


def _decode_bytes_array(arr):
    return np.asarray(
        [x.decode("utf-8", errors="ignore") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]
    )


def _replace_h5_dataset(group, name: str, data, dtype):
    if name in group:
        del group[name]
    group.create_dataset(name, data=np.asarray(data, dtype=dtype))


def _cosine_from_triplets(pos, atom0, atom1, vertex):
    vec0 = pos[np.asarray(atom0, dtype=np.int64)] - pos[np.asarray(vertex, dtype=np.int64)]
    vec1 = pos[np.asarray(atom1, dtype=np.int64)] - pos[np.asarray(vertex, dtype=np.int64)]
    norm0 = np.linalg.norm(vec0, axis=1)
    norm1 = np.linalg.norm(vec1, axis=1)
    if np.any(norm0 <= 0.0) or np.any(norm1 <= 0.0):
        raise ValueError("Cannot build explicit backbone oxygen angle with zero-length vectors")
    cosine = np.einsum("ij,ij->i", vec0, vec1) / (norm0 * norm1)
    return np.clip(cosine, -1.0, 1.0)


def append_explicit_backbone_oxygen_geometry(
    up_file: Path,
    bond_stiffness: float = EXPLICIT_BACKBONE_OXYGEN_BOND_STIFFNESS,
    angle_stiffness: float = EXPLICIT_BACKBONE_OXYGEN_ANGLE_STIFFNESS,
):
    up_file = Path(up_file).expanduser().resolve()
    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        bb = inp["hybrid_bb_map"]

        reference_names = _decode_bytes_array(bb["reference_atom_names"][:])
        column_by_name = {name.upper(): idx for idx, name in enumerate(reference_names.tolist())}
        missing = [name for name in BB_COMPONENT_NAMES if name not in column_by_name]
        if missing:
            raise ValueError(f"Missing explicit backbone oxygen mapping columns: {missing}")

        atom_indices = bb["atom_indices"][:].astype(np.int64)
        atom_mask = bb["atom_mask"][:].astype(bool)
        required_mask = np.logical_and.reduce(
            [atom_mask[:, column_by_name[name]] for name in BB_COMPONENT_NAMES]
        )
        if not np.any(required_mask):
            raise ValueError("No residues with complete N/CA/C/O mapping found for explicit oxygen geometry")

        n_idx = atom_indices[required_mask, column_by_name["N"]]
        ca_idx = atom_indices[required_mask, column_by_name["CA"]]
        c_idx = atom_indices[required_mask, column_by_name["C"]]
        o_idx = atom_indices[required_mask, column_by_name["O"]]
        pos = inp["pos"][:, :, 0]

        bond_rows = np.column_stack([c_idx, o_idx]).astype(np.int64, copy=False)
        bond_equil = np.linalg.norm(pos[c_idx] - pos[o_idx], axis=1)

        angle_rows = [np.column_stack([ca_idx, o_idx, c_idx]).astype(np.int64, copy=False)]
        angle_equil = [_cosine_from_triplets(pos, ca_idx, o_idx, c_idx)]
        if len(n_idx) > 1:
            angle_rows.append(
                np.column_stack([o_idx[:-1], n_idx[1:], c_idx[:-1]]).astype(np.int64, copy=False)
            )
            angle_equil.append(_cosine_from_triplets(pos, o_idx[:-1], n_idx[1:], c_idx[:-1]))

        dist_group = pot["Distance3D"]
        spring_bond_group = pot["Spring_bond"]
        dist_ids = dist_group["id"][:]
        spring_bond_ids = spring_bond_group["id"][:]
        spring_bond_equil = spring_bond_group["equil_dist"][:]
        spring_bond_const = spring_bond_group["spring_const"][:]

        existing_bonds = {tuple(sorted(map(int, row.tolist()))): i for i, row in enumerate(dist_ids)}
        new_bond_rows = []
        new_bond_equil = []
        for row, equil in zip(bond_rows, bond_equil):
            key = tuple(sorted(map(int, row.tolist())))
            if key in existing_bonds:
                continue
            new_bond_rows.append(row.tolist())
            new_bond_equil.append(float(equil))
        if new_bond_rows:
            new_bond_rows = np.asarray(new_bond_rows, dtype=dist_ids.dtype)
            new_bond_equil = np.asarray(new_bond_equil, dtype=spring_bond_equil.dtype)
            bond_row_start = dist_ids.shape[0]
            dist_ids = np.concatenate([dist_ids, new_bond_rows], axis=0)
            spring_bond_ids = np.concatenate(
                [
                    spring_bond_ids,
                    np.arange(bond_row_start, bond_row_start + len(new_bond_rows), dtype=spring_bond_ids.dtype),
                ]
            )
            spring_bond_equil = np.concatenate([spring_bond_equil, new_bond_equil])
            spring_bond_const = np.concatenate(
                [
                    spring_bond_const,
                    np.full(len(new_bond_rows), bond_stiffness, dtype=spring_bond_const.dtype),
                ]
            )
            _replace_h5_dataset(dist_group, "id", dist_ids, dist_ids.dtype)
            _replace_h5_dataset(spring_bond_group, "id", spring_bond_ids, spring_bond_ids.dtype)
            _replace_h5_dataset(
                spring_bond_group, "equil_dist", spring_bond_equil, spring_bond_equil.dtype
            )
            _replace_h5_dataset(
                spring_bond_group, "spring_const", spring_bond_const, spring_bond_const.dtype
            )

        angle_group = pot["Angle"]
        spring_angle_group = pot["Spring_angle"]
        angle_ids = angle_group["id"][:]
        spring_angle_ids = spring_angle_group["id"][:]
        spring_angle_equil = spring_angle_group["equil_dist"][:]
        spring_angle_const = spring_angle_group["spring_const"][:]

        existing_angles = {
            (min(int(row[0]), int(row[1])), max(int(row[0]), int(row[1])), int(row[2])): i
            for i, row in enumerate(angle_ids)
        }
        pending_angle_rows = []
        pending_angle_equil = []
        for rows, equil_vals in zip(angle_rows, angle_equil):
            for row, equil in zip(rows, equil_vals):
                key = (min(int(row[0]), int(row[1])), max(int(row[0]), int(row[1])), int(row[2]))
                if key in existing_angles:
                    continue
                pending_angle_rows.append(row.tolist())
                pending_angle_equil.append(float(equil))
        if pending_angle_rows:
            pending_angle_rows = np.asarray(pending_angle_rows, dtype=angle_ids.dtype)
            pending_angle_equil = np.asarray(pending_angle_equil, dtype=spring_angle_equil.dtype)
            angle_row_start = angle_ids.shape[0]
            angle_ids = np.concatenate([angle_ids, pending_angle_rows], axis=0)
            spring_angle_ids = np.concatenate(
                [
                    spring_angle_ids,
                    np.arange(
                        angle_row_start,
                        angle_row_start + len(pending_angle_rows),
                        dtype=spring_angle_ids.dtype,
                    ),
                ]
            )
            spring_angle_equil = np.concatenate([spring_angle_equil, pending_angle_equil])
            spring_angle_const = np.concatenate(
                [
                    spring_angle_const,
                    np.full(len(pending_angle_rows), angle_stiffness, dtype=spring_angle_const.dtype),
                ]
            )
            _replace_h5_dataset(angle_group, "id", angle_ids, angle_ids.dtype)
            _replace_h5_dataset(spring_angle_group, "id", spring_angle_ids, spring_angle_ids.dtype)
            _replace_h5_dataset(
                spring_angle_group, "equil_dist", spring_angle_equil, spring_angle_equil.dtype
            )
            _replace_h5_dataset(
                spring_angle_group, "spring_const", spring_angle_const, spring_angle_const.dtype
            )


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

    with tb.open_file(str(up_file), mode="a") as tf:
        uc.t = tf
        uc.potential = tf.root.input.potential
        uc.n_atom = ref_n_atom
        uc.n_chains = 1
        uc.chain_starts = np.array([0], dtype=np.int32)
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

        uc.write_infer_H_O(fasta_seq, np.array([], dtype=np.int32))
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
    append_explicit_backbone_oxygen_geometry(up_file)


def inject_stage7_sc_table_nodes(
    up_file: Path,
    martini_h5: Path,
    upside_home: Path,
    rama_library: Path,
    rama_sheet_mixing: Path,
    hbond_energy: Path,
    reference_state_rama: Path,
):
    up_file = Path(up_file).expanduser().resolve()
    martini_h5 = Path(martini_h5).expanduser().resolve()
    upside_home = Path(upside_home).expanduser().resolve()
    rama_library = Path(rama_library).expanduser().resolve()
    rama_sheet_mixing = Path(rama_sheet_mixing).expanduser().resolve()
    hbond_energy = Path(hbond_energy).expanduser().resolve()
    reference_state_rama = Path(reference_state_rama).expanduser().resolve()
    sidechain_lib = (upside_home / "parameters" / "ff_2.1" / "sidechain.h5").resolve()

    if not up_file.exists():
        raise SystemExit(f"ERROR: stage file not found: {up_file}")
    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")
    for path in [rama_library, rama_sheet_mixing, hbond_energy, reference_state_rama, sidechain_lib]:
        if not path.exists():
            raise SystemExit(f"ERROR: required Upside input not found: {path}")

    with h5py.File(martini_h5, "r") as sc_lib:
        required_sc_datasets = [
            "restype_order",
            "target_order",
            "grid_nm",
            "cos_theta_grid",
            "rotamer_count",
            "rotamer_probability_fixed",
            "rotamer_radial_energy_kj_mol",
            "rotamer_angular_energy_kj_mol",
            "rotamer_angular_profile",
        ]
        missing_sc_datasets = [name for name in required_sc_datasets if name not in sc_lib]
        if missing_sc_datasets:
            missing_text = ", ".join(missing_sc_datasets)
            raise SystemExit(
                f"ERROR: {martini_h5} is missing required rotamer-resolved SC datasets: {missing_text}. "
                "Rebuild martini.h5 with the integrated martini_prepare_system.py build-sc-martini-h5 command."
            )
        restype_order = decode_string_array(sc_lib["restype_order"])
        target_order = decode_string_array(sc_lib["target_order"])
        grid_nm = sc_lib["grid_nm"][:].astype(np.float32)
        cos_theta_grid = sc_lib["cos_theta_grid"][:].astype(np.float32)
        rotamer_count = sc_lib["rotamer_count"][:].astype(np.int32)
        rotamer_probability_fixed = sc_lib["rotamer_probability_fixed"][:].astype(np.float32)
        rotamer_radial_energy_kj_mol = sc_lib["rotamer_radial_energy_kj_mol"][:].astype(np.float32)
        rotamer_angular_energy_kj_mol = sc_lib["rotamer_angular_energy_kj_mol"][:].astype(np.float32)
        rotamer_angular_profile = sc_lib["rotamer_angular_profile"][:].astype(np.float32)

    restype_to_index = {name: i for i, name in enumerate(restype_order)}
    target_to_index = {name: i for i, name in enumerate(target_order)}

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]

        residue_ids, affine_atoms = build_affine_atoms(inp)
        sequence = resolve_sequence(inp, len(residue_ids))
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
        g_sc.attrs["energy_conversion_kj_per_eup"] = np.float32(martini_potential.attrs["energy_conversion_kj_per_eup"])
        g_sc.attrs["length_conversion_angstrom_per_nm"] = np.float32(
            martini_potential.attrs["length_conversion_angstrom_per_nm"]
        )
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
        g_sc.create_dataset("grid_nm", data=grid_nm, dtype=np.float32)
        g_sc.create_dataset("cos_theta_grid", data=cos_theta_grid, dtype=np.float32)
        g_sc.create_dataset("rotamer_count", data=rotamer_count.astype(np.int32), dtype=np.int32)
        g_sc.create_dataset("rotamer_probability_fixed", data=rotamer_probability_fixed, dtype=np.float32)
        g_sc.create_dataset("rotamer_radial_energy_kj_mol", data=rotamer_radial_energy_kj_mol, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_energy_kj_mol", data=rotamer_angular_energy_kj_mol, dtype=np.float32)
        g_sc.create_dataset("rotamer_angular_profile", data=rotamer_angular_profile, dtype=np.float32)
        g_sc.create_dataset("restype_order", data=np.asarray([np.bytes_(x) for x in restype_order], dtype="S4"))
        g_sc.create_dataset("target_order", data=np.asarray([np.bytes_(x) for x in target_order], dtype="S8"))

    print(
        f"Injected rotamer-weighted martini_sc_table_1body into {up_file}: "
        f"n_rows={len(rotamer_payload['id_seq'])} skipped={len(rotamer_payload['skipped'])} "
        f"n_env={len(env_atom_index)} n_restypes={len(restype_order)} n_targets={len(target_order)}"
    )


if __name__ == "__main__":
    raise SystemExit("Use martini_prepare_system.py as the workflow preparation entrypoint.")
