#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
import tables as tb

NA_AVOGADRO = 6.02214076e23
BB_COMPONENT_NAMES = ("N", "CA", "C", "O")
BB_COMPONENT_MASSES = (14.0, 12.0, 12.0, 16.0)


@dataclass
class Config:
    protein_pdb: Path
    bilayer_pdb: Path
    output_dir: Path
    protein_cg_pdb: Optional[Path]
    protein_itp: Optional[Path]
    martinize_cmd: Optional[str]
    salt_molar: float
    protein_lipid_cutoff: float
    ion_cutoff: float
    box_padding_xy: float
    box_padding_z: float
    seed: int
    protein_net_charge: Optional[int]


def parse_args() -> Config:
    parser = argparse.ArgumentParser(
        description=(
            "Pack OPM protein into MARTINI DOPC bilayer and export hybrid "
            "mapping artifacts for Upside + dry MARTINI integration."
        )
    )
    parser.add_argument("--protein-pdb", default="pdb/1rkl.pdb")
    parser.add_argument("--bilayer-pdb", default="pdb/bilayer.MARTINI.pdb")
    parser.add_argument("--output-dir", default="outputs/hybrid_1rkl")
    parser.add_argument(
        "--protein-cg-pdb",
        default=None,
        help="Optional existing MARTINI protein PDB. If not provided, --martinize-cmd is required.",
    )
    parser.add_argument(
        "--protein-itp",
        default=None,
        help=(
            "Optional MARTINI protein ITP from martinize.py. If provided, sidechain-to-backbone "
            "force-transfer mapping is derived from bonded topology and force constants."
        ),
    )
    parser.add_argument(
        "--martinize-cmd",
        default=None,
        help=(
            "Optional martinize command template used when --protein-cg-pdb is not given. "
            "Use placeholders {input} and {output}."
        ),
    )
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--protein-lipid-cutoff", type=float, default=3.0)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--box-padding-xy", type=float, default=0.0)
    parser.add_argument("--box-padding-z", type=float, default=20.0)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--protein-net-charge",
        type=int,
        default=None,
        help="Override inferred protein charge for ion placement.",
    )

    args = parser.parse_args()
    return Config(
        protein_pdb=Path(args.protein_pdb).expanduser().resolve(),
        bilayer_pdb=Path(args.bilayer_pdb).expanduser().resolve(),
        output_dir=Path(args.output_dir).expanduser().resolve(),
        protein_cg_pdb=Path(args.protein_cg_pdb).expanduser().resolve()
        if args.protein_cg_pdb
        else None,
        protein_itp=Path(args.protein_itp).expanduser().resolve()
        if args.protein_itp
        else None,
        martinize_cmd=args.martinize_cmd,
        salt_molar=args.salt_molar,
        protein_lipid_cutoff=args.protein_lipid_cutoff,
        ion_cutoff=args.ion_cutoff,
        box_padding_xy=args.box_padding_xy,
        box_padding_z=args.box_padding_z,
        seed=args.seed,
        protein_net_charge=args.protein_net_charge,
    )


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


def infer_protein_charge_from_cg(protein_atoms):
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


def run_martinize(config: Config, out_path: Path):
    if not config.martinize_cmd:
        raise ValueError(
            "protein_cg_pdb not provided and martinize_cmd missing. "
            "Provide --protein-cg-pdb or --martinize-cmd."
        )
    cmd_str = config.martinize_cmd.format(
        input=str(config.protein_pdb),
        output=str(out_path),
    )
    cmd = shlex.split(cmd_str)
    subprocess.run(cmd, check=True)
    if not out_path.exists():
        raise FileNotFoundError(f"martinize output not found: {out_path}")


def choose_protein_cg(config: Config) -> Path:
    if config.protein_cg_pdb is not None:
        return config.protein_cg_pdb

    out_path = config.output_dir / "protein.martini.pdb"
    run_martinize(config, out_path)
    return out_path


def extract_protein_cg_atoms(cg_atoms):
    aa_res = {
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
    }
    protein_like = []
    for atom in cg_atoms:
        seg = atom["segid"].upper()
        resname = atom["resname"].upper()
        if seg.startswith("PRO") or resname in aa_res:
            protein_like.append(atom)
    if not protein_like:
        raise ValueError("No protein-like MARTINI atoms found in CG PDB.")
    return protein_like


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


def collect_backbone_com_pairs(protein_aa_atoms, protein_cg_atoms):
    aa_by_res = defaultdict(dict)
    for idx, atom in enumerate(protein_aa_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"])
        aa_by_res[key][atom["name"].strip().upper()] = idx

    aa_points = []
    cg_points = []
    for atom in protein_cg_atoms:
        if atom["name"].strip().upper() != "BB":
            continue

        found_key = None
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in aa_by_res:
            found_key = key
        else:
            cands = [k for k in aa_by_res.keys() if k[1] == atom["resseq"]]
            if len(cands) == 1:
                found_key = cands[0]
        if found_key is None:
            continue

        idxs = [aa_by_res[found_key].get(nm, -1) for nm in BB_COMPONENT_NAMES]
        valid = [int(i) for i in idxs if i is not None and int(i) >= 0]
        if not valid:
            continue

        aa_xyz = np.array(
            [[protein_aa_atoms[i]["x"], protein_aa_atoms[i]["y"], protein_aa_atoms[i]["z"]] for i in valid],
            dtype=float,
        )
        aa_com = np.mean(aa_xyz, axis=0)
        cg_bb = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
        aa_points.append(aa_com)
        cg_points.append(cg_bb)

    if not aa_points:
        return np.zeros((0, 3), dtype=float), np.zeros((0, 3), dtype=float)

    return np.array(aa_points, dtype=float), np.array(cg_points, dtype=float)


def compute_backbone_com_alignment(protein_aa_atoms, protein_cg_atoms):
    p, q = collect_backbone_com_pairs(protein_aa_atoms, protein_cg_atoms)
    if p.shape[0] == 0:
        return np.eye(3, dtype=float), np.zeros(3, dtype=float), 0.0, 0

    p_center = np.mean(p, axis=0)
    q_center = np.mean(q, axis=0)

    # Use Kabsch when at least 3 matched COM points exist; otherwise
    # keep identity rotation and use translation-only alignment.
    if p.shape[0] >= 3:
        p0 = p - p_center
        q0 = q - q_center
        cov = p0.T @ q0
        u, _, vt = np.linalg.svd(cov)
        rot = u @ vt
        if np.linalg.det(rot) < 0.0:
            vt[-1, :] *= -1.0
            rot = u @ vt
    else:
        rot = np.eye(3, dtype=float)

    trans = q_center - p_center @ rot
    aligned = p @ rot + trans
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - q) ** 2, axis=1))))
    return rot, trans, rmsd, int(p.shape[0])


def summarize_backbone_frame_alignment(protein_aa_atoms, protein_cg_atoms):
    p, q = collect_backbone_com_pairs(protein_aa_atoms, protein_cg_atoms)
    matched = int(p.shape[0])
    summary = {
        "matched_residues": matched,
        "raw_rmsd_angstrom": 0.0,
        "translation_only_rmsd_angstrom": 0.0,
        "rigid_rmsd_angstrom": 0.0,
        "centroid_shift_angstrom": 0.0,
        "rotation_angle_deg": 0.0,
        "radius_gyration_ratio": 1.0,
        "frame_offset_detected": False,
    }
    if matched == 0:
        return summary

    p_center = np.mean(p, axis=0)
    q_center = np.mean(q, axis=0)
    shift = q_center - p_center
    p_shift = p + shift

    raw_rmsd = float(np.sqrt(np.mean(np.sum((p - q) ** 2, axis=1))))
    trans_rmsd = float(np.sqrt(np.mean(np.sum((p_shift - q) ** 2, axis=1))))
    rot, _trans, rigid_rmsd, _n = compute_backbone_com_alignment(protein_aa_atoms, protein_cg_atoms)

    trace_val = float(np.trace(rot))
    cos_theta = max(-1.0, min(1.0, 0.5 * (trace_val - 1.0)))
    rotation_angle = float(np.degrees(np.arccos(cos_theta)))
    rg_aa = float(np.sqrt(np.mean(np.sum((p - p_center) ** 2, axis=1))))
    rg_cg = float(np.sqrt(np.mean(np.sum((q - q_center) ** 2, axis=1))))
    rg_ratio = float(rg_aa / rg_cg) if rg_cg > 1.0e-12 else 1.0
    centroid_shift = float(np.linalg.norm(shift))

    summary.update(
        {
            "raw_rmsd_angstrom": raw_rmsd,
            "translation_only_rmsd_angstrom": trans_rmsd,
            "rigid_rmsd_angstrom": float(rigid_rmsd),
            "centroid_shift_angstrom": centroid_shift,
            "rotation_angle_deg": rotation_angle,
            "radius_gyration_ratio": rg_ratio,
            "frame_offset_detected": bool(
                centroid_shift > 1.0
                or rotation_angle > 5.0
                or (trans_rmsd - float(rigid_rmsd)) > 1.0
            ),
        }
    )
    return summary


def validate_backbone_reference_frame(
    protein_aa_atoms,
    protein_cg_atoms,
    min_matched_residues: int = 8,
    max_rigid_rmsd: float = 1.5,
    context: str = "",
):
    diag = summarize_backbone_frame_alignment(protein_aa_atoms, protein_cg_atoms)
    matched = int(diag["matched_residues"])
    label = context if context else "AA->MARTINI BB"

    if matched < int(min_matched_residues):
        raise ValueError(
            "Backbone frame preflight failed: insufficient matched residues between "
            f"AA backbone and MARTINI BB for {label}. "
            f"Matched={matched}, required>={int(min_matched_residues)}. "
            "Check chain/residue numbering and ensure AA/CG structures correspond to the same protein."
        )

    rigid_rmsd = float(diag["rigid_rmsd_angstrom"])
    if rigid_rmsd > float(max_rigid_rmsd):
        raise ValueError(
            "Backbone frame preflight failed: AA backbone and MARTINI BB are inconsistent "
            f"after optimal rigid alignment for {label}. "
            f"rigid_rmsd={rigid_rmsd:.4f} Angstrom, limit={float(max_rigid_rmsd):.4f} Angstrom, "
            f"translation_only_rmsd={float(diag['translation_only_rmsd_angstrom']):.4f} Angstrom, "
            f"centroid_shift={float(diag['centroid_shift_angstrom']):.4f} Angstrom. "
            "Check that AA and CG inputs are generated from the same structure and residue indexing."
        )

    print(
        "Backbone frame preflight (AA N/CA/C/O vs MARTINI BB): "
        f"context={label}; matched={matched}; "
        f"rigid_rmsd={rigid_rmsd:.4f} Angstrom; "
        f"translation_only_rmsd={float(diag['translation_only_rmsd_angstrom']):.4f} Angstrom; "
        f"centroid_shift={float(diag['centroid_shift_angstrom']):.4f} Angstrom; "
        f"rotation={float(diag['rotation_angle_deg']):.2f} deg"
    )
    if bool(diag["frame_offset_detected"]):
        print(
            "Detected AA/BB frame offset (shift and/or rotation); "
            "preparation will keep using rigid AA->BB alignment for reference backbone coordinates."
        )
    return diag


def collect_bb_map(protein_aa_atoms, protein_cg_atoms):
    aa_by_res = defaultdict(dict)
    for idx, atom in enumerate(protein_aa_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"])
        aa_by_res[key][atom["name"].strip().upper()] = idx

    align_rot, align_trans, align_rmsd, align_n = compute_backbone_com_alignment(
        protein_aa_atoms, protein_cg_atoms
    )
    if align_n > 0:
        print(
            "Backbone COM RMSD alignment (AA->MARTINI BB): "
            f"n={align_n}, rmsd={align_rmsd:.4f} Å"
        )

    bb_entries = []
    for cg_idx, atom in enumerate(protein_cg_atoms):
        if atom["name"].upper() != "BB":
            continue
        found_key = None
        # Primary matching by chain+resseq+icode.
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key in aa_by_res:
            found_key = key
        else:
            # Fallback: unique residue number match.
            cands = [k for k in aa_by_res.keys() if k[1] == atom["resseq"]]
            if len(cands) == 1:
                found_key = cands[0]
        if found_key is None:
            continue

        # Reference all-atom backbone atom mapping (for audit/comment metadata).
        aa_idxs = [int(aa_by_res[found_key].get(nm, -1)) for nm in BB_COMPONENT_NAMES]
        aa_coords = []
        for ai in aa_idxs:
            if ai < 0:
                aa_coords.append([0.0, 0.0, 0.0])
            else:
                aa_atom = protein_aa_atoms[ai]
                aa_vec = np.array([float(aa_atom["x"]), float(aa_atom["y"]), float(aa_atom["z"])], dtype=float)
                aa_aligned = aa_vec @ align_rot + align_trans
                aa_coords.append([float(aa_aligned[0]), float(aa_aligned[1]), float(aa_aligned[2])])

        # Active mapping is kept in protein-AA PDB index space.
        # Runtime conversion to stage-local indices is done during stage-file injection.
        idxs = [-1, -1, -1, -1]
        mask = [0, 0, 0, 0]
        raw_weights = [0.0, 0.0, 0.0, 0.0]
        for d, (ai, m) in enumerate(zip(aa_idxs, BB_COMPONENT_MASSES)):
            if ai < 0:
                continue
            idxs[d] = int(ai)
            mask[d] = 1
            raw_weights[d] = float(m)

        weights = [0.0, 0.0, 0.0, 0.0]
        wsum = float(sum(raw_weights))
        if wsum > 0.0:
            weights = [w / wsum for w in raw_weights]
        bb_comment = (
            f"BB residue {atom['resseq']} chain '{atom['chain']}' "
            f"ref N/CA/C/O idx={aa_idxs}; index_space=protein_aa_pdb_0based; "
            f"align_rmsd={align_rmsd:.4f}; w={weights}"
        )
        bb_entries.append(
            {
                "bb_resseq": atom["resseq"],
                "bb_chain": atom["chain"],
                "bb_icode": atom["icode"],
                "bb_atom_index": cg_idx,
                "atom_indices": idxs,
                "atom_mask": mask,
                "weights": weights,
                "reference_atom_indices": aa_idxs,
                "reference_atom_coords": aa_coords,
                "bb_comment": bb_comment,
            }
        )
    return bb_entries


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
        ctrl.attrs["enable"] = np.int8(1)
        ctrl.attrs["activation_stage"] = b"production"
        ctrl.attrs["preprod_protein_mode"] = b"rigid"
        ctrl.attrs["preprod_lipid_headgroup_roles"] = b"PO4"
        ctrl.attrs["exclude_intra_protein_martini"] = np.int8(1)
        ctrl.attrs["production_nonprotein_hard_sphere"] = np.int8(0)
        ctrl.attrs["coupling_align_debug"] = np.int8(0)
        ctrl.attrs["coupling_align_interval"] = np.int32(100)
        ctrl.attrs["schema_version"] = np.int32(1)

        bb_grp = inp.create_group("hybrid_bb_map")
        bb_grp.attrs["atom_index_space"] = b"protein_aa_pdb_0based"
        bb_grp.attrs["runtime_index_space"] = b"stage_runtime_after_injection"
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


def write_summary(path: Path, payload: Dict):
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def main():
    config = parse_args()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)

    protein_aa_atoms, _ = parse_pdb(config.protein_pdb)
    bilayer_atoms, bilayer_box = parse_pdb(config.bilayer_pdb)
    protein_cg_path = choose_protein_cg(config)
    protein_cg_atoms_raw, _ = parse_pdb(protein_cg_path)
    protein_cg_atoms = extract_protein_cg_atoms(protein_cg_atoms_raw)

    bilayer_lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    bilayer_xyz = coords(bilayer_lipid_atoms if bilayer_lipid_atoms else bilayer_atoms)
    protein_xyz = coords(protein_cg_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    protein_center = center_of_mass(protein_xyz)
    shift = bilayer_center - protein_center
    translated = protein_xyz + shift
    set_coords(protein_cg_atoms, translated)
    if protein_aa_atoms:
        protein_aa_xyz = coords(protein_aa_atoms)
        set_coords(protein_aa_atoms, protein_aa_xyz + shift)

    protein_xyz = coords(protein_cg_atoms)
    pmin = protein_xyz.min(axis=0)
    pmax = protein_xyz.max(axis=0)
    pspan = pmax - pmin
    pcenter_xy = 0.5 * (pmin[:2] + pmax[:2])
    min_required_xy = 3.0 * pspan[:2] + 2.0 * config.box_padding_xy
    square_side = float(np.max(min_required_xy))
    target_xy_min = np.array(
        [pcenter_xy[0] - 0.5 * square_side, pcenter_xy[1] - 0.5 * square_side],
        dtype=float,
    )
    target_xy_max = np.array(
        [pcenter_xy[0] + 0.5 * square_side, pcenter_xy[1] + 0.5 * square_side],
        dtype=float,
    )

    bilayer_lipids = tile_and_crop_bilayer_lipids(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        target_xy_min=target_xy_min,
        target_xy_max=target_xy_max,
    )

    lipid_residues, keep_nonlipid = compute_lipid_residue_indices(bilayer_lipids)
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_lipids,
        protein_atoms=protein_cg_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=config.protein_lipid_cutoff,
    )

    packed_atoms = protein_cg_atoms + bilayer_kept
    min_box_z_target = float(3.0 * pspan[2])

    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        pad_z=config.box_padding_z,
        min_box_z=min_box_z_target,
    )
    bilayer_xyz_post = coords(bilayer_kept)
    lipid_mid_z = float(0.5 * (bilayer_xyz_post[:, 2].min() + bilayer_xyz_post[:, 2].max()))

    effective_vol_frac = infer_effective_ion_volume_fraction_from_template(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        salt_molar=config.salt_molar,
    )

    protein_charge = (
        config.protein_net_charge
        if config.protein_net_charge is not None
        else infer_protein_charge_from_cg(protein_cg_atoms)
    )
    salt_pairs = estimate_salt_pairs(
        box_lengths,
        config.salt_molar,
        effective_volume_fraction=effective_vol_frac,
    )
    n_na = salt_pairs + max(0, -protein_charge)
    n_cl = salt_pairs + max(0, protein_charge)
    ion_atoms = place_ions(
        atoms=packed_atoms,
        box_lengths=box_lengths,
        n_na=n_na,
        n_cl=n_cl,
        cutoff=config.ion_cutoff,
        rng=rng,
    )

    all_atoms = packed_atoms + ion_atoms

    packed_pdb = config.output_dir / "hybrid_packed.MARTINI.pdb"
    write_pdb(packed_pdb, all_atoms, box_lengths)

    bb_entries = collect_bb_map(protein_aa_atoms, protein_cg_atoms)
    mapping_h5 = config.output_dir / "hybrid_mapping.h5"
    env_atom_indices = list(range(len(protein_cg_atoms), len(all_atoms)))
    write_hybrid_mapping_h5(
        mapping_h5,
        bb_entries=bb_entries,
        total_martini_atoms=len(all_atoms),
        env_atom_indices=env_atom_indices,
        n_protein_atoms=len(protein_cg_atoms),
    )

    mapping_json = config.output_dir / "hybrid_bb_map.json"
    write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})

    summary = {
        "protein_aa_pdb": str(config.protein_pdb),
        "protein_cg_pdb": str(protein_cg_path),
        "bilayer_pdb": str(config.bilayer_pdb),
        "output_pdb": str(packed_pdb),
        "mapping_h5": str(mapping_h5),
        "box_angstrom": [float(v) for v in box_lengths],
        "lipid_mid_z": float(lipid_mid_z),
        "box_half_z": float(0.5 * box_lengths[2]),
        "protein_z_span": float(pspan[2]),
        "min_box_z_target": float(min_box_z_target),
        "protein_charge_used": int(protein_charge),
        "salt_molar": float(config.salt_molar),
        "salt_pairs_target": int(salt_pairs),
        "ion_effective_volume_fraction": float(effective_vol_frac),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "protein_atoms_cg": int(len(protein_cg_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "lipid_residues_removed": int(removed_lipids),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "bb_map_entries": int(len(bb_entries)),
    }
    write_summary(config.output_dir / "hybrid_prep_summary.json", summary)

    print(f"Packed system written to: {packed_pdb}")
    print(f"Hybrid mapping HDF5 written to: {mapping_h5}")
    print(f"Summary written to: {config.output_dir / 'hybrid_prep_summary.json'}")



# -----------------------------------------------------------------------------
# Stage Conversion Helpers
# -----------------------------------------------------------------------------
def runtime_input_pdb_path(script_dir, pdb_id):
    """Resolve MARTINI PDB path with optional runtime override."""
    override = os.environ.get("UPSIDE_RUNTIME_PDB_FILE", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(script_dir, f"pdb/{pdb_id}.MARTINI.pdb")


def runtime_protein_itp_path(script_dir, pdb_id):
    """Resolve protein ITP path with optional runtime override."""
    override = os.environ.get("UPSIDE_RUNTIME_ITP_FILE", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(script_dir, f"pdb/{pdb_id}_proa.itp")


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

# --- Helpers: read MARTINI protein topology (ITP) for protein bead typing and connectivity ---
def read_protein_itp_topology(itp_path: str):
    topo_by_res_and_role = {}
    if not os.path.exists(itp_path):
        return topo_by_res_and_role
    in_atoms = False
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                in_atoms = line.lower().startswith('[ atoms')
                continue
            if not in_atoms:
                continue
            # tokens: idx type resnr residue atom cgnr charge ...
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                _idx = int(parts[0])
                bead_type = parts[1]
                resnr = int(parts[2])
                residue = parts[3]
                atom_role = parts[4]  # e.g., BB, SC1, SC2
                # charge is usually the last numeric token; try to parse last column
                charge = 0.0
                for tok in reversed(parts):
                    try:
                        charge = float(tok)
                        break
                    except ValueError:
                        continue
                topo_by_res_and_role[(resnr, atom_role)] = (bead_type, charge)
            except Exception:
                continue
    return topo_by_res_and_role

def read_protein_itp_connectivity(itp_path: str, simulation_stage='minimization'):
    """
    Read bonds, angles, dihedrals, constraints, and position restraints from protein ITP file
    simulation_stage: 'minimization' or 'production'
    """
    bonds = []
    angles = []
    dihedrals = []
    constraints = []
    position_restraints = []
    
    if not os.path.exists(itp_path):
        return bonds, angles, dihedrals, constraints, position_restraints
    
    current_section = None
    in_normang_section = False
    in_nonnormang_section = False
    in_flexible_section = False
    in_posres_section = False
    
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            
            # Handle preprocessor directives
            if line.startswith('#ifdef NORMANG'):
                in_normang_section = True
                in_nonnormang_section = False
                continue
            elif line.startswith('#ifndef NORMANG'):
                in_normang_section = False
                in_nonnormang_section = True
                continue
            elif line.startswith('#ifdef FLEXIBLE'):
                in_flexible_section = True
                continue
            elif line.startswith('#ifndef FLEXIBLE'):
                in_flexible_section = False
                continue
            elif line.startswith('#ifdef POSRES'):
                in_posres_section = True
                continue
            elif line.startswith('#ifndef POSRES'):
                # Don't reset in_posres_section here - let #endif handle it
                continue
            elif line.startswith('#ifndef POSRES_FC'):
                # This is a different directive - don't reset in_posres_section
                continue
            elif line.startswith('#endif'):
                # Only reset flags if we're not in a position restraints section
                if current_section != 'position_restraints':
                    in_normang_section = False
                    in_nonnormang_section = False
                    in_flexible_section = False
                    in_posres_section = False
                continue
            
            if line.startswith('['):
                section = line.lower()
                if 'bonds' in section:
                    current_section = 'bonds'
                elif 'angles' in section:
                    current_section = 'angles'
                elif 'dihedrals' in section:
                    current_section = 'dihedrals'
                elif 'constraints' in section:
                    current_section = 'constraints'
                elif 'position_restraints' in section:
                    current_section = 'position_restraints'
                else:
                    current_section = None
                continue
            
            if current_section is None:
                continue
                
            parts = line.split()
            if len(parts) < 3:
                continue
                
            try:
                if current_section == 'bonds':
                    # Format: i j func r0 k
                    i, j = int(parts[0])-1, int(parts[1])-1  # Convert to 0-indexed
                    func = int(parts[2])
                    if func == 1 and len(parts) >= 5:  # Harmonic bond
                        r0 = float(parts[3])  # nm
                        k = float(parts[4])   # kJ/mol/nm²
                        
                        # For minimization: use FLEXIBLE bonds (large spring constants)
                        # For production: use regular bonds
                        if simulation_stage == 'minimization' and in_flexible_section:
                            k = 1000000.0  # Large spring constant for minimization
                        
                        bonds.append((i, j, r0, k))
                elif current_section == 'angles':
                    # Format: i j k func theta0 k
                    # For minimization: use NORMANG section
                    # For production: use regular angles section
                    if simulation_stage == 'minimization' and in_normang_section and not in_nonnormang_section:
                        i, j, k = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1
                        func = int(parts[3])
                        if func == 2 and len(parts) >= 6:  # Harmonic angle
                            theta0 = float(parts[4])  # degrees
                            k = float(parts[5])       # kJ/mol/rad²
                            angles.append((i, j, k, theta0, k))
                    elif simulation_stage == 'production' and not in_normang_section and not in_nonnormang_section:
                        # Regular angles section for production
                        i, j, k = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1
                        func = int(parts[3])
                        if func == 2 and len(parts) >= 6:  # Harmonic angle
                            theta0 = float(parts[4])  # degrees
                            k = float(parts[5])       # kJ/mol/rad²
                            angles.append((i, j, k, theta0, k))
                elif current_section == 'constraints':
                    # Format: i j func r0
                    # These should be treated as bonds with large spring constants
                    i, j = int(parts[0])-1, int(parts[1])-1  # Convert to 0-indexed
                    func = int(parts[2])
                    if func == 1 and len(parts) >= 4:  # Constraint
                        r0 = float(parts[3])  # nm
                        k = 1000000.0  # Large spring constant as requested
                        constraints.append((i, j, r0, k))
                elif current_section == 'position_restraints':
                    # Format: i func fx fy fz
                    # Only process during minimization stage
                    if simulation_stage == 'minimization' and in_posres_section:
                        i = int(parts[0])-1  # Convert to 0-indexed
                        func = int(parts[1])
                        if func == 1 and len(parts) >= 5:  # Position restraint
                            fx = float(parts[2]) if parts[2] != 'POSRES_FC' else 1000.0
                            fy = float(parts[3]) if parts[3] != 'POSRES_FC' else 1000.0
                            fz = float(parts[4]) if parts[4] != 'POSRES_FC' else 1000.0
                            position_restraints.append((i, fx, fy, fz))
                elif current_section == 'dihedrals':
                    # Format: i j k l func phi0 k mult
                    atom_i, atom_j, atom_k, atom_l = int(parts[0])-1, int(parts[1])-1, int(parts[2])-1, int(parts[3])-1
                    func = int(parts[4])
                    if func == 1 and len(parts) >= 8:  # Periodic dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const, func))
                    elif func == 2 and len(parts) >= 7:  # Harmonic dihedral
                        phi0 = float(parts[5])  # degrees
                        force_const = float(parts[6])     # kJ/mol/rad²
                        dihedrals.append((atom_i, atom_j, atom_k, atom_l, phi0, force_const, func))
            except (ValueError, IndexError):
                continue
    
    return bonds, angles, dihedrals, constraints, position_restraints

def read_protein_itp_exclusions(itp_path: str):
    """Read exclusions from protein ITP file"""
    exclusions = []
    
    if not os.path.exists(itp_path):
        return exclusions
    
    current_section = None
    with open(itp_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                section = line.lower()
                if 'exclusions' in section:
                    current_section = 'exclusions'
                else:
                    current_section = None
                continue
            
            if current_section == 'exclusions':
                parts = line.split()
                if len(parts) >= 2:
                    # Convert to 0-indexed and add all pairs
                    atoms = [int(part)-1 for part in parts]
                    for i in range(len(atoms)):
                        for j in range(i+1, len(atoms)):
                            exclusions.append((atoms[i], atoms[j]))
    
    return exclusions

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
    # Get script directory to ensure correct path resolution
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ff_file_path = os.path.join(SCRIPT_DIR, ff_file)
    
    if not os.path.exists(ff_file_path):
        raise ValueError(f"FATAL ERROR: Force field file '{ff_file}' not found.\n"
                        f"  Full path: {ff_file_path}\n"
                        f"  This file is required for atom type masses.\n"
                        f"  Please ensure the MARTINI force field file exists and is readable.\n"
                        f"  Aborting to prevent incorrect simulation results.")
    
    in_atomtypes = False
    with open(ff_file_path, 'r') as f:
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
            'barostat_type': 0  # Berendsen
        },
        'npt_equil': {
            'lj_soften': 1,
            'lj_alpha': 0.2,
            'coulomb_soften': 1,
            'slater_alpha': 2.0,
            'barostat_type': 0  # Berendsen
        },
        'npt_equil_reduced': {
            'lj_soften': 1,
            'lj_alpha': 0.05,
            'coulomb_soften': 1,
            'slater_alpha': 0.5,
            'barostat_type': 0  # Berendsen
        },
        'npt_prod': {
            'lj_soften': 0,
            'lj_alpha': 0.0,
            'coulomb_soften': 0,
            'slater_alpha': 0.0,
            'barostat_type': 1  # Parrinello-Rahman
        }
    }

    # Get parameters for current stage
    params = stage_params.get(stage, stage_params['npt_prod'])
    stage_lipidhead_fc = float(os.environ.get('UPSIDE_BILAYER_LIPIDHEAD_FC', '0'))


    # Configuration
    strict_from_martini_pdb = True
    include_protein = True
    
    print("=== Dry MARTINI Protein-Lipid System Preparation ===")
    print(f"PDB ID: {pdb_id}")
    print(f"Output directory: {run_dir}")
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    # Read dry MARTINI parameter files
    print("\n=== Reading Dry MARTINI Parameters ===")
    ff_dir = os.environ.get('UPSIDE_MARTINI_FF_DIR', 'ff_dry')
    ff_path = os.path.join(SCRIPT_DIR, ff_dir)

    if not os.path.isdir(ff_path):
        raise ValueError(f"FATAL ERROR: Force-field directory '{ff_path}' not found.\n"
                        f"  Set UPSIDE_MARTINI_FF_DIR to a valid directory under {SCRIPT_DIR}.")

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
    
    # Determine system type and load appropriate topology
    # Check if this is a protein system by looking for protein ITP file
    protein_itp = runtime_protein_itp_path(SCRIPT_DIR, pdb_id)
    has_protein = os.path.exists(protein_itp)
    
    # Check if this is a mixed protein-lipid system by looking for both protein and lipid residues
    # This will be determined during PDB parsing, but we can prepare for both cases
    
    if has_protein:
        print("=== Mixed Protein-Lipid System Detected ===")
        print(f"Using protein topology from: {protein_itp}")
        print("System contains both protein and lipid components")
        print("Using dry MARTINI force field parameters for both protein and lipid components")
    else:
        print("=== Lipid System Detected ===")
    
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
        input_pdb_file = runtime_input_pdb_path(SCRIPT_DIR, pdb_id)
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
    
    # Load protein topology mapping and connectivity if available
    protein_itp = runtime_protein_itp_path(SCRIPT_DIR, pdb_id)
    protein_topo_map = read_protein_itp_topology(protein_itp)
    
    # Parse protein connectivity for minimization stage (uses FLEXIBLE, NORMANG, POSRES sections)
    protein_bonds, protein_angles, protein_dihedrals, protein_constraints, protein_position_restraints = read_protein_itp_connectivity(protein_itp, 'minimization')
    protein_exclusions = read_protein_itp_exclusions(protein_itp)
    
    print(f"\n=== Protein Connectivity for Minimization Stage ===")
    print(f"Bonds: {len(protein_bonds)} (including FLEXIBLE bonds with large spring constants)")
    print(f"Angles: {len(protein_angles)} (from NORMANG section)")
    print(f"Dihedrals: {len(protein_dihedrals)}")
    print(f"Constraints: {len(protein_constraints)} (as bonds with large spring constants)")
    print(f"Position restraints: {len(protein_position_restraints)} (POSRES atoms - will be ignored)")
    
    
    protein_residue_names = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
        'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX'
    }
    # martinize outputs can renumber residues (e.g., start from 1) while the packed
    # PDB can preserve original residue ids. Track sequential protein residues from
    # the PDB stream so we can map roles against ITP resnr robustly.
    protein_resseq_to_seqidx = {}
    protein_seq_next = 1

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
            
            # Determine if this line corresponds to protein beads.
            # Accept both legacy PROA-tagged MARTINI PDBs and martinize.py outputs
            # that identify protein only by residue names/atom roles.
            is_protein = (
                ('PROA' in line) or
                (residue_name in protein_residue_names) or
                ((residue_id, atom_name) in protein_topo_map)
            )
            
            # Map to MARTINI type based on context
            if is_protein:
                # Prefer exact topology mapping by (resnr, role) when available
                role = atom_name  # BB / SC1 / SC2 ... expected in our MARTINI PDB
                res_key = (chain_id, residue_id, residue_name)
                if res_key not in protein_resseq_to_seqidx:
                    protein_resseq_to_seqidx[res_key] = protein_seq_next
                    protein_seq_next += 1
                seq_resnr = protein_resseq_to_seqidx[res_key]

                if (residue_id, role) in protein_topo_map:
                    martini_type, charge = protein_topo_map[(residue_id, role)]
                elif (seq_resnr, role) in protein_topo_map:
                    martini_type, charge = protein_topo_map[(seq_resnr, role)]
                else:
                    # Raise error for unknown protein atom
                    raise ValueError(f"FATAL ERROR: Unknown protein atom '{atom_name}' in residue '{residue_name}'.\n"
                                   f"  This indicates incomplete protein topology mapping in the ITP file.\n"
                                   f"  Please ensure the protein topology file '{protein_itp}' contains proper mapping for this atom.\n"
                                   f"  Aborting to prevent incorrect simulation results.")
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
        print(f"\n=== Protein Connectivity from {protein_itp} ===")
        print(f"Found {len(protein_bonds)} bonds, {len(protein_angles)} angles, {len(protein_dihedrals)} dihedrals")
        print(f"Found {len(protein_constraints)} constraints, {len(protein_position_restraints)} position restraints")
        
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
        print(f"\nNo protein connectivity found in {protein_itp}")
    
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
            # Get mass from force field file, raise error if not found
            if atom_type not in martini_masses:
                raise ValueError(f"FATAL ERROR: Mass not found for atom type '{atom_type}' (atom index {i}).\n"
                                f"  Available atom types with masses: {sorted(martini_masses.keys())}\n"
                                f"  This indicates incomplete force field parameters.\n"
                                f"  Aborting to prevent incorrect simulation results.")
            # Divide by 12.0 for reduced mass units (1 unit = 12 g/mol)
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
        
        # Store production stage angle parameters (regular angles)
        prod_angles_grp = t.create_group(stage_grp, 'production_angles')
        # For production, we'd use different angle parameters (not NORMANG)
        # For now, use the same as minimization but this could be different
        prod_angle_fc = np.array([angle[4] for angle in protein_angles], dtype='f4')
        t.create_array(prod_angles_grp, 'force_constants', obj=prod_angle_fc)
        
        print(f"Stage-specific parameters: minimization bonds={len(min_bond_fc)}, production bonds={len(prod_bond_fc)}")
        print(f"Stage-specific parameters: minimization angles={len(min_angle_fc)}, production angles={len(prod_angle_fc)}")
        
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
            # Legacy isotropic compressibility (kept for compatibility)
            legacy_compressibility = float(os.environ.get('UPSIDE_NPT_COMPRESSIBILITY', '14.521180763676'))
            barostat_grp._v_attrs.compressibility = legacy_compressibility
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
            # Barostat type: 0 = Berendsen (default), 1 = Parrinello-Rahman
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
                particle_class[i] = b"PROTEIN"
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
        # Optional softening/overwrite controls via environment variables
        # UPSIDE_SOFTEN_COULOMB: 1 to enable Slater softening for Coulomb
        # UPSIDE_SLATER_ALPHA: float value for Slater alpha (1/Angstrom)
        # UPSIDE_SOFTEN_LJ: 1 to enable soft-core LJ
        # UPSIDE_LJ_ALPHA: float value for LJ softening alpha (dimensionless)
        # UPSIDE_OVERWRITE_SPLINES: 1 to truncate spline debug files before writing
        soften_coul = int(os.environ.get('UPSIDE_SOFTEN_COULOMB', '0'))
        slater_alpha = float(os.environ.get('UPSIDE_SLATER_ALPHA', '1.0'))
        soften_lj = int(os.environ.get('UPSIDE_SOFTEN_LJ', '0'))
        lj_alpha = float(os.environ.get('UPSIDE_LJ_ALPHA', '1.0'))
        overwrite_splines = int(os.environ.get('UPSIDE_OVERWRITE_SPLINES', '0'))
        
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
        martini_potential._v_attrs.overwrite_spline_tables = overwrite_splines
        
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

        martini_potential._v_attrs.debug_mode = 1  # Enable spline table generation
        martini_potential._v_attrs.force_debug_mode = 1  # Enable force debugging for charged particles
        
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
            bond_group._v_attrs.debug_mode = 1  # Enable spline table generation

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
            angle_group._v_attrs.debug_mode = 1  # Enable spline table generation

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
            dihedral_group._v_attrs.debug_mode = 1  # Enable spline table generation

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

def create_production_input(input_file, pdb_id):
    """Create production input file with production-stage parameters"""
    print(f"\n=== Creating Production Input File ===")
    
    # Parse protein connectivity for production stage (uses regular sections)
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    protein_itp = runtime_protein_itp_path(SCRIPT_DIR, pdb_id)
    protein_bonds_prod, protein_angles_prod, protein_dihedrals_prod, protein_constraints_prod, protein_position_restraints_prod = read_protein_itp_connectivity(protein_itp, 'production')
    
    print(f"Production stage parameters:")
    print(f"  Bonds: {len(protein_bonds_prod)} (regular bonds)")
    print(f"  Angles: {len(protein_angles_prod)} (regular angles)")
    print(f"  Dihedrals: {len(protein_dihedrals_prod)}")
    print(f"  Constraints: {len(protein_constraints_prod)}")
    print(f"  Position restraints: {len(protein_position_restraints_prod)} (none for production)")
    
    # Create production input file
    production_file = f"outputs/martini_test/production.input.up"
    
    # Copy the minimization file and modify for production
    import shutil
    shutil.copy2(input_file, production_file)
    
    # Remove position restraints from production file
    import h5py
    with h5py.File(production_file, 'r+') as f:
        if 'input' in f and 'position_restraints' in f['input']:
            del f['input']['position_restraints']
            print("Removed position restraints from production file")
    
    print(f"Production input file created: {production_file}")
    print("Note: Production file uses regular parameters without position restraints")

def main_always_fixed(pdb_id):
    """Create input file with protein always fixed rigid throughout simulation"""
    print(f"\n=== Creating Input with Always-Fixed Protein ===")
    print(f"PDB ID: {pdb_id}")
    print("Output directory: outputs/martini_test")
    
    # Use the same setup as conversion path but modify the fix rigid configuration
    # First, run the normal preparation
    convert_stage(pdb_id=pdb_id)
    
    # Then modify the H5 file to ensure fix rigid is always enabled
    import h5py
    input_file = f"outputs/martini_test/test.input.up"
    
    print(f"\n=== Modifying H5 File for Always-Fixed Protein ===")
    with h5py.File(input_file, 'r+') as f:
        if 'input' in f and 'fix_rigid' in f['input']:
            fix_rigid_grp = f['input']['fix_rigid']
            # Ensure fix rigid is always enabled
            fix_rigid_grp.attrs['enable'] = 1
            print("Fix rigid enabled for entire simulation")
            
            # Remove stage-specific parameters since we want fix rigid throughout
            if 'stage_parameters' in f['input']:
                del f['input']['stage_parameters']
                print("Removed stage-specific parameters (using fix rigid throughout)")
        else:
            print("WARNING: No fix rigid group found in H5 file")
    
    print(f"Modified input file: {input_file}")
    print("Note: Protein will be fixed rigid throughout entire simulation")


if __name__ == "__main__":
    raise SystemExit("Use prepare_system.py as the workflow preparation entrypoint.")
