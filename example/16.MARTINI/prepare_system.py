#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import json
import os
import shutil
import sys
import _pickle as cPickle
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import tables as tb

from lib import (
    DRY_NONBOND_SCHEMA,
    center_of_mass,
    collect_bb_map,
    collect_sc_map,
    convert_stage,
    compute_lipid_residue_indices,
    coords,
    estimate_salt_pairs,
    extract_protein_cg_atoms,
    infer_effective_ion_volume_fraction_from_template,
    infer_protein_charge_from_cg,
    lipid_resname,
    load_dry_nonbond_artifact,
    parse_pdb,
    place_ions,
    read_martini3_nonbond_params,
    remove_overlapping_lipids,
    set_box_from_lipid_xy,
    tile_and_crop_bilayer_lipids,
    validate_backbone_reference_frame,
    wrap_xyz_array_to_box,
    write_dry_nonbond_artifact,
    write_hybrid_mapping_h5,
    write_pdb,
)


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Unified preparation script for bilayer-only, protein-only, or mixed "
            "protein+bilayer systems. Can optionally convert prepared structure "
            "to UPSIDE input."
        )
    )
    parser.add_argument("--mode", choices=["bilayer", "protein", "both"], required=True)
    parser.add_argument("--pdb-id", required=True, help="Runtime PDB id for stage conversion")
    parser.add_argument("--runtime-pdb-output", default=None)
    parser.add_argument("--runtime-itp-output", default=None)
    parser.add_argument("--prepare-structure", type=int, default=1, choices=[0, 1])
    parser.add_argument("--stage", default=None, help="stage name for UPSIDE input conversion")
    parser.add_argument("--run-dir", default="outputs/martini_test")

    parser.add_argument("--bilayer-pdb", default="pdb/bilayer.MARTINI.pdb")
    parser.add_argument("--protein-cg-pdb", default=None)
    parser.add_argument("--protein-aa-pdb", default=None)
    parser.add_argument("--protein-itp", default=None)
    parser.add_argument("--hybrid-mapping-output", default=None)
    parser.add_argument("--hybrid-bb-map-json-output", default=None)

    parser.add_argument("--xy-scale", type=float, default=1.0)
    parser.add_argument("--box-padding-xy", type=float, default=0.0)
    parser.add_argument("--box-padding-z", type=float, default=20.0)
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=2026)

    parser.add_argument("--protein-lipid-cutoff", type=float, default=3.0)
    parser.add_argument("--protein-net-charge", type=int, default=None)
    parser.add_argument(
        "--bb-aa-min-matched-residues",
        type=int,
        default=8,
        help=(
            "Minimum matched residues required between AA backbone (N/CA/C/O) "
            "and MARTINI BB for hybrid-mapping frame preflight."
        ),
    )
    parser.add_argument(
        "--bb-aa-max-rigid-rmsd",
        type=float,
        default=1.5,
        help=(
            "Maximum allowed rigid-fit RMSD (Angstrom) between AA backbone COM "
            "and MARTINI BB in hybrid-mapping frame preflight."
        ),
    )
    parser.add_argument("--summary-json", default=None)
    return parser.parse_args(argv)


def runtime_paths(args):
    runtime_pdb = (
        Path(args.runtime_pdb_output).expanduser().resolve()
        if args.runtime_pdb_output
        else (SCRIPT_DIR / "pdb" / f"{args.pdb_id}.MARTINI.pdb")
    )
    runtime_itp = (
        Path(args.runtime_itp_output).expanduser().resolve()
        if args.runtime_itp_output
        else (SCRIPT_DIR / "pdb" / f"{args.pdb_id}_proa.itp")
    )
    return runtime_pdb, runtime_itp


def runtime_pdb_has_martini_protein(runtime_pdb: Path) -> bool:
    protein_residue_names = {
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
    with runtime_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            residue_name = line[17:21].strip().upper()
            seg_id = line[72:76].strip().upper() if len(line) >= 76 else ""
            if residue_name in protein_residue_names or seg_id.startswith("PRO"):
                return True
    return False


def copy_if_different(src: Path, dst: Path):
    src_resolved = src.expanduser().resolve()
    dst_resolved = dst.expanduser().resolve()
    if src_resolved == dst_resolved:
        return
    shutil.copy2(src_resolved, dst_resolved)


def infer_protein_charge_from_atoms(protein_atoms):
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


def extract_protein_aa_atoms(aa_atoms):
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
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "CYX",
    }
    protein_atoms = []
    for atom in aa_atoms:
        if atom["record"] != "ATOM":
            continue
        if atom["resname"].upper() not in aa_res:
            continue
        protein_atoms.append(atom)
    if not protein_atoms:
        raise ValueError("No protein amino-acid atoms found in AA PDB.")
    return protein_atoms


def collect_aa_backbone_entries(protein_aa_atoms):
    comp_names = ("N", "CA", "C", "O")
    comp_mass = {"N": 14.0, "CA": 12.0, "C": 12.0, "O": 16.0}
    aa_by_res = defaultdict(dict)
    residue_order = []
    seen_res = set()

    for atom in protein_aa_atoms:
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key not in seen_res:
            residue_order.append(key)
            seen_res.add(key)
        aa_by_res[key][atom["name"].strip().upper()] = atom

    entries = []
    next_ref_idx = 0
    for chain, resseq, icode in residue_order:
        row_idx = []
        row_mask = []
        row_coords = []
        row_mass = []
        res_atoms = aa_by_res[(chain, resseq, icode)]
        for name in comp_names:
            atom = res_atoms.get(name)
            if atom is None:
                row_idx.append(-1)
                row_mask.append(0)
                row_coords.append([0.0, 0.0, 0.0])
                row_mass.append(0.0)
                continue
            row_idx.append(next_ref_idx)
            row_mask.append(1)
            row_coords.append([float(atom["x"]), float(atom["y"]), float(atom["z"])])
            row_mass.append(comp_mass[name])
            next_ref_idx += 1

        mass_sum = sum(row_mass)
        if mass_sum > 0.0:
            weights = [float(m / mass_sum) for m in row_mass]
        else:
            weights = [0.0, 0.0, 0.0, 0.0]

        entries.append(
            {
                "bb_resseq": int(resseq),
                "bb_atom_index": -1,
                "atom_indices": list(row_idx),
                "atom_mask": list(row_mask),
                "weights": weights,
                "reference_atom_indices": list(row_idx),
                "reference_atom_coords": row_coords,
                "bb_comment": "AA-only backbone carrier row (no MARTINI protein proxy)",
            }
        )

    return entries


def parse_itp_atom_types(itp_path: Path):
    atom_types = set()
    in_atoms = False
    with itp_path.open("r", encoding="utf-8", errors="ignore") as fh:
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
            if in_atoms:
                parts = line.split()
                if len(parts) >= 2:
                    atom_types.add(parts[1])
    return sorted(atom_types)


def parse_ff_mass_atom_types(ff_file: Path):
    mass_types = set()
    in_atomtypes = False
    with ff_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low in {"[ atomtypes ]", "[atomtypes]"}:
                in_atomtypes = True
                continue
            if line.startswith("[") and line.endswith("]"):
                in_atomtypes = False
                continue
            if not in_atomtypes:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                float(parts[1])
            except ValueError:
                continue
            mass_types.add(parts[0])
    return sorted(mass_types)


def assert_protein_itp_mass_compatibility(runtime_itp: Path):
    ff_dir = os.environ.get("UPSIDE_MARTINI_FF_DIR", "ff_dry")
    ff_file = SCRIPT_DIR / ff_dir / "dry_martini_v2.1.itp"
    if not ff_file.exists():
        raise FileNotFoundError(
            f"Dry MARTINI mass table not found: {ff_file}. "
            "Set UPSIDE_MARTINI_FF_DIR to a valid force-field directory."
        )

    itp_types = parse_itp_atom_types(runtime_itp)
    ff_types = set(parse_ff_mass_atom_types(ff_file))
    missing = sorted(t for t in itp_types if t not in ff_types)
    if missing:
        raise ValueError(
            "Protein ITP contains bead types missing in dry-MARTINI mass table.\n"
            f"  ITP: {runtime_itp}\n"
            f"  FF:  {ff_file}\n"
            f"  Missing types: {missing}\n"
            "Use a dry-MARTINI-compatible protein ITP (e.g., from the hybrid path or "
            "martinize settings that match the selected force field)."
        )


def prepare_bilayer_structure(args, runtime_pdb):
    bilayer_pdb = Path(args.bilayer_pdb).expanduser().resolve()
    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not lipid_atoms:
        raise ValueError("No lipid residues detected in bilayer template.")

    if args.xy_scale < 1.0:
        raise ValueError("--xy-scale must be >= 1.0")

    lip_xyz = coords(lipid_atoms)
    lip_min = lip_xyz.min(axis=0)
    lip_max = lip_xyz.max(axis=0)
    lip_center = center_of_mass(lip_xyz)

    span_x = float(lip_max[0] - lip_min[0])
    span_y = float(lip_max[1] - lip_min[1])
    base_side = max(span_x, span_y)
    target_side = base_side * float(args.xy_scale)
    target_xy_min = np.array(
        [lip_center[0] - 0.5 * target_side, lip_center[1] - 0.5 * target_side],
        dtype=float,
    )
    target_xy_max = np.array(
        [lip_center[0] + 0.5 * target_side, lip_center[1] + 0.5 * target_side],
        dtype=float,
    )

    bilayer_lipids = tile_and_crop_bilayer_lipids(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        target_xy_min=target_xy_min,
        target_xy_max=target_xy_max,
    )

    box_lengths = set_box_from_lipid_xy(
        all_atoms=bilayer_lipids,
        lipid_atoms=bilayer_lipids,
        pad_z=float(args.box_padding_z),
        force_square_xy=True,
        min_box_z=None,
        center_lipid_in_z=True,
    )

    effective_vol_frac = infer_effective_ion_volume_fraction_from_template(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        salt_molar=float(args.salt_molar),
    )
    salt_pairs = estimate_salt_pairs(
        box_lengths=box_lengths,
        salt_molar=float(args.salt_molar),
        effective_volume_fraction=effective_vol_frac,
    )

    rng = np.random.default_rng(int(args.seed))
    ion_atoms = place_ions(
        atoms=bilayer_lipids,
        box_lengths=box_lengths,
        n_na=salt_pairs,
        n_cl=salt_pairs,
        cutoff=float(args.ion_cutoff),
        rng=rng,
    )
    all_atoms = bilayer_lipids + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)

    return {
        "mode": "bilayer",
        "input_bilayer_pdb": str(bilayer_pdb),
        "runtime_pdb": str(runtime_pdb),
        "xy_scale": float(args.xy_scale),
        "base_xy_side_angstrom": float(base_side),
        "target_xy_side_angstrom": float(target_side),
        "box_angstrom": [float(v) for v in box_lengths],
        "salt_molar": float(args.salt_molar),
        "ion_effective_volume_fraction": float(effective_vol_frac),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(salt_pairs),
        "cl_added": int(salt_pairs),
        "lipid_atoms": int(len(bilayer_lipids)),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
    }


def prepare_protein_structure(args, runtime_pdb, runtime_itp):
    if not args.protein_cg_pdb:
        raise ValueError("--protein-cg-pdb is required for mode=protein")
    protein_cg_pdb = Path(args.protein_cg_pdb).expanduser().resolve()
    if not protein_cg_pdb.exists():
        raise FileNotFoundError(f"Protein CG PDB not found: {protein_cg_pdb}")

    protein_atoms_raw, protein_box = parse_pdb(protein_cg_pdb)
    protein_atoms = extract_protein_cg_atoms(protein_atoms_raw)
    write_pdb(runtime_pdb, protein_atoms, protein_box)
    if args.protein_itp:
        protein_itp = Path(args.protein_itp).expanduser().resolve()
        if not protein_itp.exists():
            raise FileNotFoundError(f"Protein ITP not found: {protein_itp}")
        copy_if_different(protein_itp, runtime_itp)

    return {
        "mode": "protein",
        "input_protein_cg_pdb": str(protein_cg_pdb),
        "runtime_pdb": str(runtime_pdb),
        "runtime_itp": str(runtime_itp) if args.protein_itp else None,
        "total_atoms": int(len(protein_atoms)),
        "box_angstrom": [float(v) for v in protein_box] if protein_box else None,
    }


def prepare_mixed_structure(args, runtime_pdb, runtime_itp):
    if not args.bilayer_pdb:
        raise ValueError("--bilayer-pdb is required for mode=both")
    if not args.protein_aa_pdb:
        raise ValueError("--protein-aa-pdb is required for mode=both")

    bilayer_pdb = Path(args.bilayer_pdb).expanduser().resolve()
    if not bilayer_pdb.exists():
        raise FileNotFoundError(f"Bilayer PDB not found: {bilayer_pdb}")
    protein_aa_pdb = Path(args.protein_aa_pdb).expanduser().resolve()
    if not protein_aa_pdb.exists():
        raise FileNotFoundError(f"Protein AA PDB not found: {protein_aa_pdb}")

    protein_aa_atoms_raw, _ = parse_pdb(protein_aa_pdb)
    protein_aa_atoms = extract_protein_aa_atoms(protein_aa_atoms_raw)
    if not protein_aa_atoms:
        raise ValueError(f"No atoms found in protein AA PDB: {protein_aa_pdb}")

    protein_cg_pdb = None
    protein_atoms = None
    runtime_protein_atoms = []
    uses_martini_protein = bool(args.protein_cg_pdb)
    if uses_martini_protein:
        protein_cg_pdb = Path(args.protein_cg_pdb).expanduser().resolve()
        if not protein_cg_pdb.exists():
            raise FileNotFoundError(f"Protein CG PDB not found: {protein_cg_pdb}")
        protein_atoms_raw, _ = parse_pdb(protein_cg_pdb)
        protein_atoms = extract_protein_cg_atoms(protein_atoms_raw)
        runtime_protein_atoms = protein_atoms
        if not protein_atoms:
            raise ValueError("No atoms found in protein CG PDB.")
    else:
        protein_atoms = protein_aa_atoms

    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    bilayer_lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not bilayer_lipid_atoms:
        raise ValueError("No lipid residues found in bilayer template.")

    protein_xyz = coords(protein_atoms)
    bilayer_xyz = coords(bilayer_lipid_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    protein_center = center_of_mass(protein_xyz)
    bilayer_alignment_shift = protein_center - bilayer_center
    bilayer_all_xyz = coords(bilayer_atoms)
    bilayer_all_xyz = bilayer_all_xyz + bilayer_alignment_shift
    for atom, c in zip(bilayer_atoms, bilayer_all_xyz):
        atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])

    pmin = protein_xyz.min(axis=0)
    pmax = protein_xyz.max(axis=0)
    pspan = pmax - pmin
    pcenter_xy = 0.5 * (pmin[:2] + pmax[:2])

    bmin = bilayer_xyz.min(axis=0)
    bmax = bilayer_xyz.max(axis=0)
    base_side = max(float(bmax[0] - bmin[0]), float(bmax[1] - bmin[1]))
    min_required_xy = 3.0 * pspan[:2] + 2.0 * float(args.box_padding_xy)
    target_side = max(base_side * float(max(1.0, args.xy_scale)), float(np.max(min_required_xy)))

    target_xy_min = np.array(
        [pcenter_xy[0] - 0.5 * target_side, pcenter_xy[1] - 0.5 * target_side],
        dtype=float,
    )
    target_xy_max = np.array(
        [pcenter_xy[0] + 0.5 * target_side, pcenter_xy[1] + 0.5 * target_side],
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
        protein_atoms=protein_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=float(args.protein_lipid_cutoff),
    )

    packed_atoms = runtime_protein_atoms + bilayer_kept
    packed_xyz_before_box = coords(packed_atoms)
    min_box_z_target = float(3.0 * pspan[2])
    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        pad_z=float(args.box_padding_z),
        force_square_xy=True,
        min_box_z=min_box_z_target,
        center_lipid_in_z=True,
    )
    packed_xyz_after_box = coords(packed_atoms)
    if packed_xyz_before_box.shape == packed_xyz_after_box.shape and packed_xyz_after_box.size:
        box_shift = packed_xyz_after_box[0] - packed_xyz_before_box[0]
        protein_aa_xyz = coords(protein_aa_atoms)
        protein_aa_xyz = protein_aa_xyz + box_shift
        for atom, c in zip(protein_aa_atoms, protein_aa_xyz):
            atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])

    protein_aa_xyz = coords(protein_aa_atoms)
    lipid_kept_xyz = coords(bilayer_kept)
    min_protein_lipid_distance = float("nan")
    if protein_aa_xyz.size and lipid_kept_xyz.size:
        min_d2 = float("inf")
        chunk = 2000
        for i in range(0, lipid_kept_xyz.shape[0], chunk):
            block = lipid_kept_xyz[i:i + chunk]
            d = block[:, None, :] - protein_aa_xyz[None, :, :]
            d2 = np.einsum("ijk,ijk->ij", d, d)
            block_min = float(d2.min())
            if block_min < min_d2:
                min_d2 = block_min
        min_protein_lipid_distance = float(min_d2 ** 0.5)

    effective_vol_frac = infer_effective_ion_volume_fraction_from_template(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        salt_molar=float(args.salt_molar),
    )
    protein_charge = (
        int(args.protein_net_charge)
        if args.protein_net_charge is not None
        else int(infer_protein_charge_from_atoms(protein_aa_atoms))
    )
    salt_pairs = estimate_salt_pairs(
        box_lengths=box_lengths,
        salt_molar=float(args.salt_molar),
        effective_volume_fraction=effective_vol_frac,
    )
    n_na = int(salt_pairs + max(0, -protein_charge))
    n_cl = int(salt_pairs + max(0, protein_charge))

    rng = np.random.default_rng(int(args.seed))
    ion_atoms = place_ions(
        atoms=packed_atoms,
        box_lengths=box_lengths,
        n_na=n_na,
        n_cl=n_cl,
        cutoff=float(args.ion_cutoff),
        rng=rng,
    )

    all_atoms = packed_atoms + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)

    if uses_martini_protein and args.protein_itp:
        protein_itp = Path(args.protein_itp).expanduser().resolve()
        if not protein_itp.exists():
            raise FileNotFoundError(f"Protein ITP not found: {protein_itp}")
        copy_if_different(protein_itp, runtime_itp)

    mapping_summary = {}
    if args.hybrid_mapping_output:
        protein_itp_for_map = None
        if uses_martini_protein and args.protein_itp:
            protein_itp_for_map = Path(args.protein_itp).expanduser().resolve()
        elif uses_martini_protein and runtime_itp.exists():
            protein_itp_for_map = runtime_itp

        frame_diag = None
        if uses_martini_protein:
            frame_diag = validate_backbone_reference_frame(
                protein_aa_atoms,
                protein_atoms,
                min_matched_residues=max(1, int(args.bb_aa_min_matched_residues)),
                max_rigid_rmsd=float(args.bb_aa_max_rigid_rmsd),
                context=f"{protein_aa_pdb.name} -> {protein_cg_pdb.name}",
            )
            bb_entries = collect_bb_map(protein_aa_atoms, protein_atoms)
            sc_entries = collect_sc_map(
                protein_aa_atoms,
                protein_atoms,
                protein_itp_path=protein_itp_for_map,
            )
            env_atom_indices = list(range(len(runtime_protein_atoms), len(all_atoms)))
            n_protein_atoms = len(runtime_protein_atoms)
        else:
            bb_entries = collect_aa_backbone_entries(protein_aa_atoms)
            sc_entries = []
            env_atom_indices = list(range(len(all_atoms)))
            n_protein_atoms = 0

        mapping_h5 = Path(args.hybrid_mapping_output).expanduser().resolve()
        mapping_h5.parent.mkdir(parents=True, exist_ok=True)
        write_hybrid_mapping_h5(
            mapping_h5,
            bb_entries=bb_entries,
            sc_entries=sc_entries,
            total_martini_atoms=len(all_atoms),
            env_atom_indices=env_atom_indices,
            n_protein_atoms=n_protein_atoms,
        )

        mapping_summary["protein_aa_pdb"] = str(protein_aa_pdb)
        mapping_summary["mapping_h5"] = str(mapping_h5)
        if frame_diag is not None:
            mapping_summary["backbone_frame_check"] = frame_diag
        mapping_summary["bb_map_entries"] = int(len(bb_entries))
        mapping_summary["sc_map_entries"] = int(len(sc_entries))

        if args.hybrid_bb_map_json_output:
            mapping_json = Path(args.hybrid_bb_map_json_output).expanduser().resolve()
            write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})
            mapping_summary["mapping_json"] = str(mapping_json)

    summary = {
        "mode": "both",
        "protein_representation": "martini_cg_plus_upside_mapping" if uses_martini_protein else "upside_backbone_only",
        "input_protein_cg_pdb": str(protein_cg_pdb) if protein_cg_pdb is not None else None,
        "input_protein_aa_pdb": str(protein_aa_pdb),
        "input_bilayer_pdb": str(bilayer_pdb),
        "runtime_pdb": str(runtime_pdb),
        "runtime_itp": str(runtime_itp) if (uses_martini_protein and args.protein_itp) else None,
        "bilayer_alignment_shift_angstrom": [float(v) for v in bilayer_alignment_shift],
        "input_protein_center_angstrom": [float(v) for v in protein_center],
        "input_bilayer_center_angstrom": [float(v) for v in bilayer_center],
        "xy_scale": float(args.xy_scale),
        "base_xy_side_angstrom": float(base_side),
        "target_xy_side_angstrom": float(target_side),
        "box_angstrom": [float(v) for v in box_lengths],
        "protein_charge_used": int(protein_charge),
        "min_protein_lipid_distance": float(min_protein_lipid_distance),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "protein_atoms_runtime": int(len(runtime_protein_atoms)),
        "protein_atoms_reference": int(len(protein_aa_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "lipid_residues_removed": int(removed_lipids),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
    }
    summary.update(mapping_summary)
    return summary


def run_stage_conversion(args, runtime_pdb: Path, runtime_itp: Path):
    prev_pdb = os.environ.get("UPSIDE_RUNTIME_PDB_FILE")
    prev_itp = os.environ.get("UPSIDE_RUNTIME_ITP_FILE")
    uses_runtime_martini_protein = runtime_pdb_has_martini_protein(runtime_pdb)

    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(runtime_pdb)
    if uses_runtime_martini_protein:
        os.environ["UPSIDE_RUNTIME_ITP_FILE"] = str(runtime_itp)
    else:
        os.environ.pop("UPSIDE_RUNTIME_ITP_FILE", None)

    try:
        convert_stage(
            pdb_id=args.pdb_id,
            stage=args.stage,
            run_dir=args.run_dir,
        )
    finally:
        if prev_pdb is None:
            os.environ.pop("UPSIDE_RUNTIME_PDB_FILE", None)
        else:
            os.environ["UPSIDE_RUNTIME_PDB_FILE"] = prev_pdb

        if prev_itp is None:
            os.environ.pop("UPSIDE_RUNTIME_ITP_FILE", None)
        else:
            os.environ["UPSIDE_RUNTIME_ITP_FILE"] = prev_itp


def write_summary(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


BACKBONE_NODES = [
    "affine_alignment",
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

SIDECHAIN_NODES = [
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
    "placement_point_vector_only",
]

REQUIRED_BACKBONE_NODES = tuple(BACKBONE_NODES)
FORBIDDEN_SIDECHAIN_NODES = tuple(SIDECHAIN_NODES)
FORBIDDEN_MEMBRANE_RUNTIME_NODES = (
    "membrane_potential",
    "cb_membrane_potential",
    "hb_membrane_potential",
    "cb_surf_membrane_potential",
    "hb_surf_membrane_potential",
    "membranelateral_potential",
    "surface",
)
BACKBONE_CLASSES = ("N", "CA", "C", "O")
DOPC_COVERED_TYPES = ("Q0", "Qa", "Na", "C1", "C3")
SIDECHAIN_DATASET_CANDIDATES = ("sc_energy", "sidechain_energy", "sidechain")
SHEET_DATASET_CANDIDATES = ("sheet_energy", "sheet")
ROLE_PROXIES = {"N": "Qd", "CA": "C1", "C": "C2", "O": "Qa"}
ROLE_TARGET_FIELDS = {
    "N": "mean_hb_donor_attraction",
    "CA": "mean_cb_value",
    "C": "mean_cb_value",
    "O": "mean_hb_acceptor_attraction",
}
ROLE_TARGET_DESCRIPTIONS = {
    "N": "donor_hbond_mean",
    "CA": "backbone_burial_mean",
    "C": "backbone_burial_mean",
    "O": "acceptor_hbond_mean",
}
ROLE_DESCRIPTOR_NAME = "role_proxy_pair_epsilon_v1"
ROLE_RMIN_NAME = "role_proxy_pair_sigma_v1"
CROSS_CALIBRATION_MODE = "role_specific_proxy_pair_rows_v1"
RBM_CROSS_SCHEMA = "martini_backbone_cross_v1"
SPLINE_CROSS_SCHEMA = "martini_backbone_spline_cross_v1"
RBM_CROSS_NODE_NAME = "martini_rbm_cross_potential"
SPLINE_CROSS_NODE_NAME = "martini_spline_cross_potential"
REQUIRED_DEPTH_TABLE_COLUMNS = (
    "bead_name",
    "bead_type",
    "mean_depth",
    "std_depth",
    "count",
    "cb_backbone_burial0_mean",
    "cb_backbone_burial1_mean",
    "cb_backbone_mean",
    "hb_donor_unbound",
    "hb_donor_bound",
    "hb_acceptor_unbound",
    "hb_acceptor_bound",
    "hb_mean",
)
REQUIRED_CROSS_TABLE_COLUMNS = (
    "protein_class",
    "env_type",
    "proxy_type",
    "r_angstrom",
    "target_potential",
    "fitted_potential",
    "r_min_angstrom",
    "well_depth",
    "base_depth",
    "role_scale",
    "sigma_nm",
    "epsilon_pair",
    "descriptor_value",
    "fit_rmse",
)


def require_existing_file(path):
    file_path = Path(path).expanduser().resolve()
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    return file_path


def decode_bytes_array(arr):
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


def env_bool(name, default=False):
    raw = os.getenv(name)
    if raw is None:
        return default
    value = str(raw).strip().lower()
    return value not in ("0", "false", "no", "off", "")


def refresh_backbone_reference_carriers(h5f, pos):
    if "/input/hybrid_bb_map" not in h5f:
        return pos

    grp = h5f["/input/hybrid_bb_map"]
    required = ("bb_atom_index", "atom_indices", "weights", "reference_atom_coords")
    if not all((f"/input/hybrid_bb_map/{key}" in h5f) for key in required):
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
    for row_idx in range(comp_idx.shape[0]):
        bb = int(bb_idx[row_idx]) if row_idx < bb_idx.shape[0] else -1
        if bb < 0 or bb >= n_atom:
            continue

        valid = (comp_idx[row_idx] >= 0) & (comp_idx[row_idx] < n_atom) & (comp_w[row_idx] > 0.0)
        if not np.any(valid):
            continue

        weights = comp_w[row_idx][valid]
        weight_sum = float(np.sum(weights))
        if weight_sum <= 0.0:
            continue
        weights = weights / weight_sum

        ref_pts = ref_xyz[row_idx][valid]
        ref_com = np.sum(ref_pts * weights[:, None], axis=0)
        bb_pos = pos[bb, :, 0].astype(np.float64)
        shift = bb_pos - ref_com

        target_idx = comp_idx[row_idx][valid]
        aligned = ref_pts + shift[None, :]
        for atom_idx, atom_xyz in zip(target_idx, aligned):
            pos[int(atom_idx), :, 0] = atom_xyz.astype(pos.dtype)

    return pos


def set_initial_position(input_file, output_file):
    strict_copy = env_bool("UPSIDE_SET_INITIAL_STRICT_COPY", False)
    apply_refresh_backbone = env_bool(
        "UPSIDE_SET_INITIAL_REFRESH_BACKBONE_CARRIERS",
        default=(not strict_copy),
    )

    with h5py.File(input_file, "r") as src:
        if "/output/pos" in src and src["/output/pos"].shape[0] > 0:
            last_pos = src["/output/pos"][-1, 0, :, :]
            last_pos = last_pos[:, :, np.newaxis]
        else:
            last_pos = src["/input/pos"][:, :, 0]
            last_pos = last_pos[:, :, np.newaxis]

        last_box = None
        if "/output/box" in src:
            box_data = src["/output/box"][:]
            if box_data.size > 0:
                last_box = box_data[-1]
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:
                    last_box = last_box[0]
        if last_box is None and "/input/potential/martini_potential" in src:
            pot_grp = src["/input/potential/martini_potential"]
            if all(key in pot_grp.attrs for key in ("x_len", "y_len", "z_len")):
                last_box = np.array(
                    [
                        pot_grp.attrs["x_len"],
                        pot_grp.attrs["y_len"],
                        pot_grp.attrs["z_len"],
                    ]
                )

    with h5py.File(output_file, "r+") as dst:
        target_pos = dst["/input/pos"][:]
        target_n = target_pos.shape[0]
        source_n = last_pos.shape[0]
        if source_n != target_n:
            merged = target_pos.copy()
            merged[: min(source_n, target_n), :, :] = last_pos[: min(source_n, target_n), :, :]
            last_pos = merged

        if strict_copy and not apply_refresh_backbone:
            print("Strict handoff mode: preserving exact coordinates from previous stage output.")

        if apply_refresh_backbone:
            last_pos = refresh_backbone_reference_carriers(dst, last_pos)

        if last_box is None and "/input/potential/martini_potential" in dst:
            pot_grp = dst["/input/potential/martini_potential"]
            if all(key in pot_grp.attrs for key in ("x_len", "y_len", "z_len")):
                last_box = np.array(
                    [
                        pot_grp.attrs["x_len"],
                        pot_grp.attrs["y_len"],
                        pot_grp.attrs["z_len"],
                    ]
                )

        if last_box is not None:
            last_pos = wrap_xyz_array_to_box(last_pos, last_box)

        if "/input/pos" in dst:
            del dst["/input/pos"]
        dst.create_dataset("/input/pos", data=last_pos)

        if last_box is not None and "/input/potential/martini_potential" in dst:
            pot_grp = dst["/input/potential/martini_potential"]
            pot_grp.attrs["x_len"] = float(last_box[0])
            pot_grp.attrs["y_len"] = float(last_box[1])
            pot_grp.attrs["z_len"] = float(last_box[2])
            print(
                f"Updated box dimensions: x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}"
            )


def normalize_resname(name):
    aliases = {
        "HSD": "HIS",
        "HSE": "HIS",
        "HSP": "HIS",
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "CYX": "CYS",
    }
    normalized = name.strip().upper()
    return aliases.get(normalized, normalized)


def parse_backbone_itp_residue_names(path):
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
            role = parts[4].strip().upper()
            if role != "BB" or resnr in seen:
                continue
            seen.add(resnr)
            resnames.append(normalize_resname(parts[3]))
    return resnames


def parse_backbone_pdb_residue_names(path):
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
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "CYX",
    }
    resnames = []
    seen = set()
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if not raw.startswith("ATOM"):
                continue
            resname = normalize_resname(raw[17:21])
            if resname not in aa_res:
                continue
            key = (raw[21].strip(), int(raw[22:26]), raw[26].strip())
            if key in seen:
                continue
            seen.add(key)
            resnames.append(resname)
    return resnames


def parse_backbone_source_residue_names(path):
    suffix = Path(path).suffix.lower()
    if suffix == ".pdb":
        return parse_backbone_pdb_residue_names(path)
    return parse_backbone_itp_residue_names(path)


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
            f"Backbone remap index exceeds map size: max={int(idx.max())}, size={atom_map.shape[0]}"
        )
    mapped = atom_map[idx]
    if np.any(mapped < 0):
        bad = np.where(mapped < 0)[0][0]
        raise ValueError(f"Backbone remap produced negative target for ref index {int(idx[bad])}")
    out[nonneg_mask] = mapped.astype(out.dtype, copy=False)
    return out


def load_upside_config(upside_home):
    upside_config_py = Path(upside_home) / "py" / "upside_config.py"
    require_existing_file(upside_config_py)
    spec = importlib.util.spec_from_file_location("upside_config_runtime", str(upside_config_py))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_runtime_ncac_mapping(up_file):
    with h5py.File(up_file, "r") as h5:
        inp = h5["/input"]
        if "hybrid_bb_map" not in inp:
            raise ValueError("Missing /input/hybrid_bb_map; cannot map runtime N/CA/C carriers")
        bb_grp = inp["hybrid_bb_map"]
        if "bb_residue_index" not in bb_grp:
            raise ValueError("Missing /input/hybrid_bb_map/bb_residue_index")
        if "reference_atom_indices" not in bb_grp:
            raise ValueError("Missing /input/hybrid_bb_map/reference_atom_indices")
        ref_offset = int(bb_grp.attrs.get("reference_index_offset", -1))
        if ref_offset < 0:
            raise ValueError("Missing/invalid hybrid_bb_map reference_index_offset")

        bb_residue_raw = bb_grp["bb_residue_index"][:].astype(np.int32)
        bb_ref_idx = bb_grp["reference_atom_indices"][:].astype(np.int32)
        if bb_ref_idx.ndim != 2 or bb_ref_idx.shape[1] < 3:
            raise ValueError("hybrid_bb_map/reference_atom_indices must have shape (n_bb,4)")

        n_atom = int(inp["pos"].shape[0])
        residue_ids = []
        residue_to_ncac = {}
        seen = set()
        for resid, ref_row in zip(bb_residue_raw.tolist(), bb_ref_idx.tolist()):
            residue_id = int(resid)
            if residue_id not in seen:
                residue_ids.append(residue_id)
                seen.add(residue_id)
            n_idx, ca_idx, c_idx = [int(ref_row[0]), int(ref_row[1]), int(ref_row[2])]
            if n_idx < 0 or ca_idx < 0 or c_idx < 0:
                continue
            runtime_indices = (ref_offset + n_idx, ref_offset + ca_idx, ref_offset + c_idx)
            for atom_idx in runtime_indices:
                if atom_idx < 0 or atom_idx >= n_atom:
                    raise ValueError(
                        f"Backbone carrier index out of bounds for residue {residue_id}: "
                        f"idx={atom_idx}, n_atom={n_atom}"
                    )
            residue_to_ncac.setdefault(residue_id, runtime_indices)

    return residue_ids, residue_to_ncac


def write_backbone_nodes(
    upside_config,
    up_file,
    residue_ids,
    residue_to_ncac,
    residue_names,
    rama_library,
    rama_sheet_mixing,
    hbond_energy,
    reference_state_rama,
):
    if len(residue_names) != len(residue_ids):
        raise ValueError(
            f"Residue count mismatch: source has {len(residue_names)} residues, "
            f"bb map has {len(residue_ids)}"
        )

    ref_n_atom = 3 * len(residue_ids)
    atom_map = np.full((ref_n_atom,), -1, dtype=np.int64)
    for i, residue_id in enumerate(residue_ids):
        if residue_id not in residue_to_ncac:
            raise ValueError(f"Missing N/CA/C runtime mapping for residue {residue_id}")
        n_idx, ca_idx, c_idx = residue_to_ncac[residue_id]
        atom_map[3 * i + 0] = n_idx
        atom_map[3 * i + 1] = ca_idx
        atom_map[3 * i + 2] = c_idx

    with h5py.File(up_file, "r+") as up:
        inp = up["/input"]
        pot = inp["potential"]
        for node_name in SIDECHAIN_NODES + BACKBONE_NODES:
            if node_name in pot:
                del pot[node_name]
        if "sequence" in inp:
            del inp["sequence"]
        seq_data = np.asarray([np.bytes_(name) for name in residue_names], dtype="S3")
        inp.create_dataset("sequence", data=seq_data)

    fasta_seq = np.array(residue_names)
    spring_args = type(
        "SpringArgs",
        (),
        {"bond_stiffness": 48.0, "angle_stiffness": 175.0, "omega_stiffness": 30.0},
    )()

    with tb.open_file(up_file, mode="a") as tf:
        upside_config.t = tf
        upside_config.potential = tf.root.input.potential
        upside_config.n_atom = ref_n_atom
        upside_config.n_chains = 1
        upside_config.chain_starts = np.array([0], dtype=np.int32)
        upside_config.use_intensive_memory = False

        upside_config.write_dist_spring(spring_args)
        upside_config.write_angle_spring(spring_args)
        upside_config.write_omega_spring1(spring_args, fasta_seq)
        upside_config.write_rama_map_pot(
            fasta_seq,
            rama_library,
            rama_sheet_mixing,
            secstr_bias="",
            mode="mixture",
            param_deriv=False,
        )
        upside_config.write_affine_alignment(len(fasta_seq))

        ref_state_cor = np.log(cPickle.load(open(reference_state_rama, "rb"), encoding="latin1"))
        ref_state_cor -= ref_state_cor.mean()
        grp = tf.create_group(tf.root.input.potential, "rama_map_pot_ref")
        grp._v_attrs.arguments = np.array([b"rama_coord"])
        grp._v_attrs.log_pot = 0
        upside_config.create_array(grp, "residue_id", obj=np.arange(len(fasta_seq), dtype=np.int32))
        upside_config.create_array(grp, "rama_map_id", obj=np.zeros(len(fasta_seq), dtype=np.int32))
        upside_config.create_array(grp, "rama_pot", obj=ref_state_cor[None])

        upside_config.write_infer_H_O(fasta_seq, np.array([], dtype=np.int32))
        upside_config.write_count_hbond(fasta_seq, False)
        upside_config.write_short_hbond(fasta_seq, hbond_energy)
        upside_config.write_rama_coord2()
        upside_config.write_backbone_pair(fasta_seq)

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
                raise ValueError(f"Missing generated backbone node: {node_name}")
        for node_name in SIDECHAIN_NODES:
            if node_name in pot:
                raise ValueError(f"Unexpected sidechain node still present: {node_name}")
        for group_name, dataset_name in remap_datasets:
            ds = pot[group_name][dataset_name]
            ds[...] = remap_atom_index_array(ds[:], atom_map).astype(ds.dtype, copy=False)

    return ref_n_atom


def main_prepare(argv=None):
    args = parse_args(argv)
    runtime_pdb, runtime_itp = runtime_paths(args)
    runtime_pdb.parent.mkdir(parents=True, exist_ok=True)
    runtime_itp.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "mode": args.mode,
        "pdb_id": args.pdb_id,
        "prepare_structure": bool(args.prepare_structure),
        "stage": args.stage,
        "run_dir": args.run_dir,
    }

    if args.prepare_structure:
        if args.mode == "bilayer":
            summary.update(prepare_bilayer_structure(args, runtime_pdb))
        elif args.mode == "protein":
            summary.update(prepare_protein_structure(args, runtime_pdb, runtime_itp))
        else:
            summary.update(prepare_mixed_structure(args, runtime_pdb, runtime_itp))
    else:
        if not runtime_pdb.exists():
            raise FileNotFoundError(
                f"Runtime PDB not found for stage conversion: {runtime_pdb}. "
                "Run with --prepare-structure 1 first."
            )

    if args.stage:
        if runtime_pdb_has_martini_protein(runtime_pdb):
            if not runtime_itp.exists():
                raise FileNotFoundError(
                    f"Protein ITP required for MARTINI-protein stage conversion: {runtime_itp}. "
                    "Provide --protein-itp when the runtime PDB contains MARTINI protein atoms."
                )
            assert_protein_itp_mass_compatibility(runtime_itp)
        run_stage_conversion(args, runtime_pdb, runtime_itp)
        summary["upside_input"] = str(Path(args.run_dir).expanduser().resolve() / "test.input.up")

    if args.summary_json:
        summary_path = Path(args.summary_json).expanduser().resolve()
    else:
        summary_path = runtime_pdb.with_suffix(".prep_summary.json")
    write_summary(summary_path, summary)
    print(f"Preparation summary written to: {summary_path}")


def validate_stage_potential_nodes(path: Path):
    with h5py.File(path, "r") as h5:
        if "/input" not in h5:
            raise ValueError("Missing /input group")
        inp = h5["/input"]
        if "potential" not in inp:
            raise ValueError("Missing /input/potential group")
        pot = inp["potential"]
        forbidden_nodes = [name for name in FORBIDDEN_MEMBRANE_RUNTIME_NODES if name in pot]
        if forbidden_nodes:
            raise ValueError(
                "Forbidden membrane/surface runtime nodes present under /input/potential: "
                f"{forbidden_nodes}"
            )
        return {"n_potential_nodes": len(pot.keys())}


def validate_up_file(path: Path):
    with h5py.File(path, "r") as h5:
        if "/input" not in h5:
            raise ValueError("Missing /input group")
        inp = h5["/input"]
        if "potential" not in inp:
            raise ValueError("Missing /input/potential group")
        if "pos" not in inp:
            raise ValueError("Missing /input/pos dataset")
        if "sequence" not in inp:
            raise ValueError("Missing /input/sequence dataset")

        n_atom = int(inp["pos"].shape[0])
        pot = inp["potential"]
        forbidden_nodes = [name for name in FORBIDDEN_MEMBRANE_RUNTIME_NODES if name in pot]
        if forbidden_nodes:
            raise ValueError(
                "Forbidden membrane/surface runtime nodes present under /input/potential: "
                f"{forbidden_nodes}"
            )

        for node in REQUIRED_BACKBONE_NODES:
            if node not in pot:
                raise ValueError(f"Missing required backbone node: /input/potential/{node}")
        for node in FORBIDDEN_SIDECHAIN_NODES:
            if node in pot:
                raise ValueError(f"Forbidden sidechain node present: /input/potential/{node}")

        remap_targets = [
            ("Distance3D", "id"),
            ("Angle", "id"),
            ("Dihedral_omega", "id"),
            ("Dihedral_phi", "id"),
            ("Dihedral_psi", "id"),
            ("rama_coord", "id"),
            ("infer_H_O/donors", "id"),
            ("infer_H_O/acceptors", "id"),
        ]
        for group_name, dataset_name in remap_targets:
            arr = np.asarray(pot[group_name][dataset_name][:])
            finite = np.isfinite(arr)
            active = finite & (arr >= 0)
            if not np.any(active):
                continue
            idx = arr[active].astype(np.int64)
            if idx.min() < 0 or idx.max() >= n_atom:
                raise ValueError(
                    f"Index range violation in /input/potential/{group_name}/{dataset_name}: "
                    f"min={int(idx.min())} max={int(idx.max())} n_atom={n_atom}"
                )

        cross_nodes = [name for name in (RBM_CROSS_NODE_NAME, SPLINE_CROSS_NODE_NAME) if name in pot]
        if len(cross_nodes) > 1:
            raise ValueError(
                "Only one backbone cross-potential node is allowed under /input/potential; "
                f"found {cross_nodes}"
            )

        if RBM_CROSS_NODE_NAME in pot:
            cross = pot[RBM_CROSS_NODE_NAME]
            required_cross = (
                "protein_atom_indices",
                "protein_class_index",
                "env_atom_indices",
                "env_class_index",
                "radial_centers",
                "radial_widths",
                "weights",
            )
            for name in required_cross:
                if name not in cross:
                    raise ValueError(
                        f"Missing cross-potential dataset: /input/potential/{RBM_CROSS_NODE_NAME}/{name}"
                    )

            protein_atom_indices = cross["protein_atom_indices"][:].astype(np.int64)
            protein_class_index = cross["protein_class_index"][:].astype(np.int64)
            env_atom_indices = cross["env_atom_indices"][:].astype(np.int64)
            env_class_index = cross["env_class_index"][:].astype(np.int64)
            radial_centers = cross["radial_centers"][:]
            radial_widths = cross["radial_widths"][:]
            weights = cross["weights"][:]

            if protein_atom_indices.ndim != 1 or protein_class_index.ndim != 1:
                raise ValueError("Cross-potential protein index/class datasets must be rank-1")
            if env_atom_indices.ndim != 1 or env_class_index.ndim != 1:
                raise ValueError("Cross-potential environment index/class datasets must be rank-1")
            if radial_centers.ndim != 1 or radial_widths.ndim != 1:
                raise ValueError("Cross-potential radial centers/widths must be rank-1")
            if weights.ndim != 3:
                raise ValueError("Cross-potential weights must be rank-3")
            if protein_atom_indices.shape[0] != protein_class_index.shape[0]:
                raise ValueError("Cross-potential protein index/class length mismatch")
            if env_atom_indices.shape[0] != env_class_index.shape[0]:
                raise ValueError("Cross-potential environment index/class length mismatch")
            if radial_centers.shape[0] != radial_widths.shape[0]:
                raise ValueError("Cross-potential radial centers/widths length mismatch")
            if np.any(protein_atom_indices < 0) or np.any(protein_atom_indices >= n_atom):
                raise ValueError("Cross-potential protein_atom_indices out of range")
            if np.any(env_atom_indices < 0) or np.any(env_atom_indices >= n_atom):
                raise ValueError("Cross-potential env_atom_indices out of range")
            if np.any(protein_class_index < 0) or np.any(env_class_index < 0):
                raise ValueError("Cross-potential class indices must be non-negative")
            if (
                np.any(~np.isfinite(radial_centers))
                or np.any(~np.isfinite(radial_widths))
                or np.any(radial_widths <= 0)
            ):
                raise ValueError("Cross-potential radial centers/widths must be finite and widths positive")
            expected_shape = (
                int(protein_class_index.max()) + 1 if protein_class_index.size else 0,
                int(env_class_index.max()) + 1 if env_class_index.size else 0,
                int(radial_centers.shape[0]),
            )
            if weights.shape != expected_shape:
                raise ValueError(
                    f"Cross-potential weights shape mismatch: {weights.shape} vs expected {expected_shape}"
                )

        if SPLINE_CROSS_NODE_NAME in pot:
            cross = pot[SPLINE_CROSS_NODE_NAME]
            required_cross = (
                "protein_atom_indices",
                "protein_class_index",
                "env_atom_indices",
                "env_class_index",
                "radial_grid",
                "table_values",
            )
            for name in required_cross:
                if name not in cross:
                    raise ValueError(
                        f"Missing cross-potential dataset: /input/potential/{SPLINE_CROSS_NODE_NAME}/{name}"
                    )
            for attr_name in ("cutoff", "radial_min", "radial_step"):
                if attr_name not in cross.attrs:
                    raise ValueError(
                        f"Missing cross-potential attribute: /input/potential/{SPLINE_CROSS_NODE_NAME} attr {attr_name}"
                    )

            protein_atom_indices = cross["protein_atom_indices"][:].astype(np.int64)
            protein_class_index = cross["protein_class_index"][:].astype(np.int64)
            env_atom_indices = cross["env_atom_indices"][:].astype(np.int64)
            env_class_index = cross["env_class_index"][:].astype(np.int64)
            radial_grid = cross["radial_grid"][:]
            table_values = cross["table_values"][:]
            radial_min = float(cross.attrs["radial_min"])
            radial_step = float(cross.attrs["radial_step"])
            cutoff = float(cross.attrs["cutoff"])

            if protein_atom_indices.ndim != 1 or protein_class_index.ndim != 1:
                raise ValueError("Spline cross-potential protein index/class datasets must be rank-1")
            if env_atom_indices.ndim != 1 or env_class_index.ndim != 1:
                raise ValueError("Spline cross-potential environment index/class datasets must be rank-1")
            if radial_grid.ndim != 1:
                raise ValueError("Spline cross-potential radial_grid must be rank-1")
            if table_values.ndim != 3:
                raise ValueError("Spline cross-potential table_values must be rank-3")
            if protein_atom_indices.shape[0] != protein_class_index.shape[0]:
                raise ValueError("Spline cross-potential protein index/class length mismatch")
            if env_atom_indices.shape[0] != env_class_index.shape[0]:
                raise ValueError("Spline cross-potential environment index/class length mismatch")
            if np.any(protein_atom_indices < 0) or np.any(protein_atom_indices >= n_atom):
                raise ValueError("Spline cross-potential protein_atom_indices out of range")
            if np.any(env_atom_indices < 0) or np.any(env_atom_indices >= n_atom):
                raise ValueError("Spline cross-potential env_atom_indices out of range")
            if np.any(protein_class_index < 0) or np.any(env_class_index < 0):
                raise ValueError("Spline cross-potential class indices must be non-negative")
            if not np.all(np.isfinite(radial_grid)) or not np.all(np.isfinite(table_values)):
                raise ValueError("Spline cross-potential contains non-finite radial grid or table values")
            if radial_grid.size == 0:
                raise ValueError("Spline cross-potential radial_grid is empty")
            if not np.isfinite(radial_min) or not np.isfinite(radial_step) or radial_step <= 0.0:
                raise ValueError("Spline cross-potential radial_min/radial_step must be finite and radial_step positive")
            if not np.isfinite(cutoff) or cutoff <= 0.0:
                raise ValueError("Spline cross-potential cutoff must be finite and positive")
            expected_grid = radial_min + radial_step * np.arange(radial_grid.shape[0], dtype=np.float64)
            if not np.allclose(radial_grid, expected_grid, atol=1e-5):
                raise ValueError("Spline cross-potential radial_grid must match radial_min + radial_step * i")
            expected_shape = (
                int(protein_class_index.max()) + 1 if protein_class_index.size else 0,
                int(env_class_index.max()) + 1 if env_class_index.size else 0,
                int(radial_grid.shape[0]),
            )
            if table_values.shape != expected_shape:
                raise ValueError(
                    f"Spline cross-potential table_values shape mismatch: {table_values.shape} vs expected {expected_shape}"
                )

        seq = inp["sequence"][:]
        if seq.ndim != 1:
            raise ValueError("/input/sequence must be rank-1")
        if seq.shape[0] == 0:
            raise ValueError("/input/sequence is empty")

    return {"n_atom": n_atom, "n_res": int(seq.shape[0])}


def validate_table_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header row")
        missing = [col for col in REQUIRED_DEPTH_TABLE_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"CSV missing required columns: {missing}")
        rows = list(reader)
        if not rows:
            raise ValueError("CSV has no data rows")
        bead_names = set()
        for i, row in enumerate(rows):
            bead = row["bead_name"].strip()
            if not bead:
                raise ValueError(f"CSV row {i} has empty bead_name")
            if bead in bead_names:
                raise ValueError(f"CSV duplicate bead_name: {bead}")
            bead_names.add(bead)
            for key in REQUIRED_DEPTH_TABLE_COLUMNS:
                if key in ("bead_name", "bead_type"):
                    continue
                value = float(row[key])
                if not np.isfinite(value):
                    raise ValueError(f"CSV row {i} field {key} is not finite")
    return {"n_rows": len(rows)}


def validate_cross_table_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("Cross-table CSV has no header row")
        missing = [col for col in REQUIRED_CROSS_TABLE_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"Cross-table CSV missing required columns: {missing}")
        rows = list(reader)
        if not rows:
            raise ValueError("Cross-table CSV has no data rows")
        seen = set()
        protein_classes = set()
        env_types = set()
        for i, row in enumerate(rows):
            protein_class = row["protein_class"].strip()
            env_type = row["env_type"].strip()
            proxy_type = row["proxy_type"].strip()
            if not protein_class or not env_type or not proxy_type:
                raise ValueError(f"Cross-table CSV row {i} has empty class/type/proxy")
            protein_classes.add(protein_class)
            env_types.add(env_type)
            key = (protein_class, env_type, row["r_angstrom"])
            if key in seen:
                raise ValueError(f"Cross-table CSV duplicate row key: {key}")
            seen.add(key)
            for field in REQUIRED_CROSS_TABLE_COLUMNS:
                if field in ("protein_class", "env_type", "proxy_type"):
                    continue
                value = float(row[field])
                if not np.isfinite(value):
                    raise ValueError(f"Cross-table CSV row {i} field {field} is not finite")
    return {
        "n_rows": len(rows),
        "n_protein_class": len(protein_classes),
        "n_env_type": len(env_types),
    }


def validate_table_meta(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if "depth_definition" not in data:
        raise ValueError("Metadata missing depth_definition")
    depth_definition = str(data["depth_definition"])
    if "no PO4 anchor" not in depth_definition:
        raise ValueError("Metadata depth_definition must explicitly indicate no PO4 anchor")
    if "term_presence" not in data:
        raise ValueError("Metadata missing term_presence")
    term_presence = data["term_presence"]
    for key in ("backbone_cb_energy", "hbond_hb_energy"):
        if key not in term_presence or not bool(term_presence[key]):
            raise ValueError(f"Metadata term_presence indicates missing required term: {key}")
    if data.get("aggregation_mode") == "empirical_depth_samples":
        if "type_sample_stats" not in data or not isinstance(data["type_sample_stats"], dict):
            raise ValueError("Refined depth metadata missing type_sample_stats")
        if "symmetry_report" not in data or not isinstance(data["symmetry_report"], dict):
            raise ValueError("Refined depth metadata missing symmetry_report")
    return {"term_presence": term_presence}


def validate_cross_table_meta(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    for key in (
        "sources",
        "protein_classes",
        "environment_types",
        "radial_grid",
        "calibration",
        "role_proxies",
        "pair_summary",
        "fit_rmse_summary",
        "overlap_report",
    ):
        if key not in data:
            raise ValueError(f"Cross-table metadata missing {key}")
    protein_classes = [str(x).strip() for x in data["protein_classes"]]
    if protein_classes != list(BACKBONE_CLASSES):
        raise ValueError(f"Cross-table metadata protein_classes mismatch: {protein_classes}")
    env_types = [str(x).strip() for x in data["environment_types"]]
    if len(env_types) != 38:
        raise ValueError(f"Cross-table metadata expected 38 environment types, got {len(env_types)}")
    calibration = data["calibration"]
    if "covered_depth_types" not in calibration:
        raise ValueError("Cross-table metadata calibration missing covered_depth_types")
    if calibration.get("calibration_mode", "") != CROSS_CALIBRATION_MODE:
        raise ValueError(
            f"Cross-table metadata calibration_mode must be {CROSS_CALIBRATION_MODE}, "
            f"got {calibration.get('calibration_mode', '')!r}"
        )
    if calibration.get("descriptor_name", "") != ROLE_DESCRIPTOR_NAME:
        raise ValueError(
            f"Cross-table metadata descriptor_name must be {ROLE_DESCRIPTOR_NAME}, "
            f"got {calibration.get('descriptor_name', '')!r}"
        )
    for key in DOPC_COVERED_TYPES:
        if key not in calibration["covered_depth_types"]:
            raise ValueError(f"Cross-table metadata missing covered type: {key}")
    fit_summary = data["fit_rmse_summary"]
    for key in ("min", "median", "max"):
        if key not in fit_summary or not np.isfinite(float(fit_summary[key])):
            raise ValueError(f"Cross-table metadata fit_rmse_summary missing/invalid {key}")
    if "symmetry_report" in data and not isinstance(data["symmetry_report"], dict):
        raise ValueError("Cross-table metadata symmetry_report must be a dict when present")
    return {
        "n_env_type": len(env_types),
        "fit_rmse_summary": fit_summary,
        "calibration_mode": calibration["calibration_mode"],
        "descriptor_name": calibration["descriptor_name"],
    }


def validate_dry_nonbond_artifact(path: Path, expected_env_classes=None):
    artifact = load_dry_nonbond_artifact(path)
    type_names = artifact["type_names"]
    sigma_angstrom = artifact["sigma_angstrom"]
    epsilon_eup = artifact["epsilon_eup"]
    if len(type_names) != 38:
        raise ValueError(f"Dry nonbond artifact expected 38 type names, got {len(type_names)}")
    if expected_env_classes is not None and list(expected_env_classes) != list(type_names):
        raise ValueError(
            "Dry nonbond artifact type_names must match /classes/environment exactly; "
            f"got {type_names}"
        )
    if np.any(sigma_angstrom <= 0.0):
        raise ValueError("Dry nonbond artifact sigma_angstrom must be positive")
    if not np.all(np.isfinite(sigma_angstrom)) or not np.all(np.isfinite(epsilon_eup)):
        raise ValueError("Dry nonbond artifact contains non-finite values")
    return {
        "schema": artifact["schema"],
        "n_type": len(type_names),
    }


def validate_cross_artifact(path: Path):
    with h5py.File(path, "r") as h5:
        required = ("/classes/protein", "/classes/environment", "/meta")
        for key in required:
            if key not in h5:
                raise ValueError(f"Cross artifact missing {key}")
        protein_classes = decode_bytes_array(h5["/classes/protein"][:]).tolist()
        env_classes = decode_bytes_array(h5["/classes/environment"][:]).tolist()
        schema = str(h5["/meta"].attrs.get("schema", ""))
        calibration_mode = str(h5["/meta"].attrs.get("calibration_mode", ""))
        descriptor_name = str(h5["/meta"].attrs.get("descriptor_name", ""))
        if protein_classes != list(BACKBONE_CLASSES):
            raise ValueError(f"Cross artifact protein classes mismatch: {protein_classes}")
        if len(env_classes) != 38:
            raise ValueError(f"Cross artifact expected 38 environment classes, got {len(env_classes)}")
        dry_nonbond_summary = validate_dry_nonbond_artifact(path, expected_env_classes=env_classes)
        if schema == RBM_CROSS_SCHEMA:
            required = (
                "/rbm/cross/radial_centers",
                "/rbm/cross/radial_widths",
                "/rbm/cross/weights",
            )
            for key in required:
                if key not in h5:
                    raise ValueError(f"Cross artifact missing {key}")
            centers = h5["/rbm/cross/radial_centers"][:]
            widths = h5["/rbm/cross/radial_widths"][:]
            weights = h5["/rbm/cross/weights"][:]
            if centers.ndim != 1 or widths.ndim != 1:
                raise ValueError("Cross artifact radial centers/widths must be rank-1")
            if weights.shape != (len(BACKBONE_CLASSES), len(env_classes), centers.shape[0]):
                raise ValueError(f"Cross artifact weights shape mismatch: {weights.shape}")
            if not np.all(np.isfinite(centers)) or not np.all(np.isfinite(widths)) or not np.all(np.isfinite(weights)):
                raise ValueError("Cross artifact contains non-finite values")
            if np.any(widths <= 0.0):
                raise ValueError("Cross artifact radial_widths must be positive")
            return {
                "schema": schema,
                "n_env_type": len(env_classes),
                "n_radial": int(centers.shape[0]),
                "dry_schema": dry_nonbond_summary["schema"],
                "n_dry_type": dry_nonbond_summary["n_type"],
            }
        if schema == SPLINE_CROSS_SCHEMA:
            required = ("/spline/cross/radial_grid", "/spline/cross/table_values")
            for key in required:
                if key not in h5:
                    raise ValueError(f"Cross artifact missing {key}")
            if calibration_mode != CROSS_CALIBRATION_MODE:
                raise ValueError(
                    f"Cross artifact calibration_mode must be {CROSS_CALIBRATION_MODE}, "
                    f"got {calibration_mode!r}"
                )
            if descriptor_name != ROLE_DESCRIPTOR_NAME:
                raise ValueError(
                    f"Cross artifact descriptor_name must be {ROLE_DESCRIPTOR_NAME}, "
                    f"got {descriptor_name!r}"
                )
            cross = h5["/spline/cross"]
            radial_grid = cross["radial_grid"][:]
            table_values = cross["table_values"][:]
            for attr_name in ("cutoff", "radial_min", "radial_step"):
                if attr_name not in cross.attrs:
                    raise ValueError(f"Cross artifact spline/cross missing attr {attr_name}")
            cutoff = float(cross.attrs["cutoff"])
            radial_min = float(cross.attrs["radial_min"])
            radial_step = float(cross.attrs["radial_step"])
            if radial_grid.ndim != 1:
                raise ValueError("Cross artifact spline radial_grid must be rank-1")
            if table_values.shape != (len(BACKBONE_CLASSES), len(env_classes), radial_grid.shape[0]):
                raise ValueError(f"Cross artifact table_values shape mismatch: {table_values.shape}")
            if not np.all(np.isfinite(radial_grid)) or not np.all(np.isfinite(table_values)):
                raise ValueError("Cross artifact spline table contains non-finite values")
            if not np.isfinite(cutoff) or cutoff <= 0.0:
                raise ValueError("Cross artifact spline cutoff must be finite and positive")
            if not np.isfinite(radial_step) or radial_step <= 0.0:
                raise ValueError("Cross artifact spline radial_step must be finite and positive")
            expected_grid = radial_min + radial_step * np.arange(radial_grid.shape[0], dtype=np.float64)
            if not np.allclose(radial_grid, expected_grid, atol=1e-5):
                raise ValueError("Cross artifact spline radial_grid must match radial_min + radial_step * i")
            return {
                "schema": schema,
                "n_env_type": len(env_classes),
                "n_radial": int(radial_grid.shape[0]),
                "calibration_mode": calibration_mode,
                "descriptor_name": descriptor_name,
                "dry_schema": dry_nonbond_summary["schema"],
                "n_dry_type": dry_nonbond_summary["n_type"],
            }
        raise ValueError(f"Unexpected cross artifact schema: {schema}")


def parse_lipid_itp_atom_types(path: Path, allowed_molecules):
    allowed = {x.upper() for x in allowed_molecules}
    atom_to_type = {}
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    current_section = ""
    current_molecule = None
    i = 0
    while i < len(lines):
        line = lines[i].split(";", 1)[0].strip()
        i += 1
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current_section = line.strip("[]").strip().lower()
            if current_section == "moleculetype":
                current_molecule = None
            continue
        if current_section == "moleculetype":
            if current_molecule is None:
                current_molecule = line.split()[0].upper()
            continue
        if current_section != "atoms" or current_molecule not in allowed:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        atom_to_type.setdefault(parts[4].strip().upper(), parts[1].strip())
    if not atom_to_type:
        raise ValueError(f"No atoms parsed from {path} for molecules {sorted(allowed)}")
    return atom_to_type


def parse_bilayer_depths(path: Path, lipid_resnames):
    allowed = {x.upper() for x in lipid_resnames}
    bead_z = defaultdict(list)
    all_lipid_z = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if resname not in allowed:
                continue
            atom_name = line[12:16].strip().upper()
            try:
                z = float(line[46:54])
            except ValueError:
                continue
            bead_z[atom_name].append(z)
            all_lipid_z.append(z)
    if not all_lipid_z:
        raise ValueError(f"No lipid atoms found in {path} for resnames {sorted(allowed)}")
    z_center = float(np.mean(np.asarray(all_lipid_z, dtype=np.float64)))
    stats = {}
    signed_samples = {}
    for bead_name, values in bead_z.items():
        signed = np.asarray(values, dtype=np.float64) - z_center
        arr = np.abs(signed)
        stats[bead_name] = {
            "count": int(arr.size),
            "mean_depth": float(np.mean(arr)),
            "std_depth": float(np.std(arr)),
        }
        signed_samples[bead_name] = signed
    return z_center, stats, signed_samples


def load_membrane_channels(path: Path):
    with h5py.File(path, "r") as h5:
        if "cb_energy" not in h5 or "hb_energy" not in h5:
            raise ValueError(f"Missing cb_energy/hb_energy in {path}")
        cb_energy = h5["cb_energy"][:].astype(np.float64)
        hb_energy = h5["hb_energy"][:].astype(np.float64)
        cb_z = np.linspace(float(h5.attrs["cb_z_min"]), float(h5.attrs["cb_z_max"]), cb_energy.shape[2])
        hb_z = np.linspace(float(h5.attrs["hb_z_min"]), float(h5.attrs["hb_z_max"]), hb_energy.shape[2])
        cb_names = [x.decode("ascii") if isinstance(x, bytes) else str(x) for x in h5["names"][:]]
        term_presence = {
            "backbone_cb_energy": True,
            "hbond_hb_energy": True,
            "sidechain_term": any(name in h5 for name in SIDECHAIN_DATASET_CANDIDATES),
            "sheet_term": any(name in h5 for name in SHEET_DATASET_CANDIDATES),
        }
    return cb_z, hb_z, cb_energy, hb_energy, cb_names, term_presence


def interp_raw(z_grid, values, z_value):
    return float(np.interp(float(z_value), z_grid, values, left=float(values[0]), right=float(values[-1])))


def interp_symmetric(z_grid, values, z_value):
    depth = abs(float(z_value))
    yp = interp_raw(z_grid, values, depth)
    yn = interp_raw(z_grid, values, -depth)
    return 0.5 * (yp + yn)


def evaluate_at_signed_depth(z_value, cb_z, hb_z, cb_energy, hb_energy, symmetrize=True):
    interp = interp_symmetric if symmetrize else interp_raw
    cb_lvl0 = np.array([interp(cb_z, cb_energy[i, 0, :], z_value) for i in range(cb_energy.shape[0])])
    cb_lvl1 = np.array([interp(cb_z, cb_energy[i, 1, :], z_value) for i in range(cb_energy.shape[0])])
    hb_values = np.array(
        [
            interp(hb_z, hb_energy[0, 0, :], z_value),
            interp(hb_z, hb_energy[0, 1, :], z_value),
            interp(hb_z, hb_energy[1, 0, :], z_value),
            interp(hb_z, hb_energy[1, 1, :], z_value),
        ],
        dtype=np.float64,
    )
    return {
        "cb_backbone_burial0_mean": float(np.mean(cb_lvl0)),
        "cb_backbone_burial1_mean": float(np.mean(cb_lvl1)),
        "cb_backbone_mean": float(np.mean(np.concatenate([cb_lvl0, cb_lvl1]))),
        "hb_donor_unbound": float(hb_values[0]),
        "hb_donor_bound": float(hb_values[1]),
        "hb_acceptor_unbound": float(hb_values[2]),
        "hb_acceptor_bound": float(hb_values[3]),
        "hb_mean": float(np.mean(hb_values)),
    }


def evaluate_at_depth(depth, cb_z, hb_z, cb_energy, hb_energy):
    return evaluate_at_signed_depth(abs(float(depth)), cb_z, hb_z, cb_energy, hb_energy, symmetrize=True)


def evaluate_at_signed_samples(z_values, cb_z, hb_z, cb_energy, hb_energy, symmetrize=True):
    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.size == 0:
        raise ValueError("Cannot evaluate membrane channels on an empty sample set")
    acc = defaultdict(list)
    for z_value in z_values.tolist():
        sample_values = evaluate_at_signed_depth(z_value, cb_z, hb_z, cb_energy, hb_energy, symmetrize=symmetrize)
        for key, value in sample_values.items():
            acc[key].append(float(value))
    return {key: float(np.mean(np.asarray(values, dtype=np.float64))) for key, values in acc.items()}


def pool_depth_samples_by_type(atom_to_type, bead_signed_samples):
    pooled = defaultdict(list)
    for bead_name, signed in bead_signed_samples.items():
        bead_type = atom_to_type.get(bead_name, "UNKNOWN")
        pooled[bead_type].extend(np.asarray(signed, dtype=np.float64).tolist())
    return {bead_type: np.asarray(values, dtype=np.float64) for bead_type, values in pooled.items()}


def summarize_depth_samples(sample_map, cb_z, hb_z, cb_energy, hb_energy):
    summary = {}
    for key, signed in sample_map.items():
        signed = np.asarray(signed, dtype=np.float64)
        if signed.size == 0:
            continue
        summary[key] = {
            "n_samples": int(signed.size),
            "mean_depth": float(np.mean(np.abs(signed))),
            "std_depth": float(np.std(np.abs(signed))),
            "mean_signed_depth": float(np.mean(signed)),
            **evaluate_at_signed_samples(signed, cb_z, hb_z, cb_energy, hb_energy, symmetrize=True),
        }
    return summary


def audit_curve_symmetry(z_grid, values, label):
    overlap_depth = min(abs(float(z_grid[0])), abs(float(z_grid[-1])))
    if overlap_depth <= 0.0:
        return {
            "label": label,
            "max_abs_diff": 0.0,
            "worst_depth": 0.0,
            "value_pos": float(values[0]),
            "value_neg": float(values[0]),
        }
    depths = np.linspace(0.0, overlap_depth, 151, dtype=np.float64)
    diffs = []
    values_pos = []
    values_neg = []
    for depth in depths.tolist():
        pos = interp_raw(z_grid, values, depth)
        neg = interp_raw(z_grid, values, -depth)
        values_pos.append(pos)
        values_neg.append(neg)
        diffs.append(abs(pos - neg))
    diffs = np.asarray(diffs, dtype=np.float64)
    best_idx = int(np.argmax(diffs))
    return {
        "label": label,
        "max_abs_diff": float(diffs[best_idx]),
        "worst_depth": float(depths[best_idx]),
        "value_pos": float(values_pos[best_idx]),
        "value_neg": float(values_neg[best_idx]),
    }


def audit_membrane_symmetry(cb_z, hb_z, cb_energy, hb_energy, cb_names):
    cb_rows = []
    hb_rows = []
    for residue_idx, residue_name in enumerate(cb_names):
        for burial_idx in range(cb_energy.shape[1]):
            cb_rows.append(
                audit_curve_symmetry(
                    cb_z,
                    cb_energy[residue_idx, burial_idx, :],
                    f"cb:{residue_name}:burial{burial_idx}",
                )
            )
    hb_labels = (
        ("donor_unbound", hb_energy[0, 0, :]),
        ("donor_bound", hb_energy[0, 1, :]),
        ("acceptor_unbound", hb_energy[1, 0, :]),
        ("acceptor_bound", hb_energy[1, 1, :]),
    )
    for label, values in hb_labels:
        hb_rows.append(audit_curve_symmetry(hb_z, values, f"hb:{label}"))

    all_rows = cb_rows + hb_rows
    global_worst = max(all_rows, key=lambda row: row["max_abs_diff"]) if all_rows else {
        "label": "none",
        "max_abs_diff": 0.0,
        "worst_depth": 0.0,
        "value_pos": 0.0,
        "value_neg": 0.0,
    }

    def summarize(rows):
        diffs = np.asarray([row["max_abs_diff"] for row in rows], dtype=np.float64)
        return {
            "n_curves": int(len(rows)),
            "mean_max_abs_diff": float(np.mean(diffs)) if diffs.size else 0.0,
            "median_max_abs_diff": float(np.median(diffs)) if diffs.size else 0.0,
            "max_abs_diff": float(np.max(diffs)) if diffs.size else 0.0,
        }

    return {
        "policy": "symmetrize_for_build_only",
        "cb_summary": summarize(cb_rows),
        "hb_summary": summarize(hb_rows),
        "global_worst": global_worst,
    }


def repo_relative(path: Path, repo_root: Path):
    try:
        return str(path.resolve().relative_to(repo_root))
    except ValueError:
        return str(path.resolve())


def parse_dry_atomtypes(path: Path):
    atomtypes = []
    in_atomtypes = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("[") and line.endswith("]"):
            in_atomtypes = line.strip("[]").strip().lower() == "atomtypes"
            continue
        if not in_atomtypes or not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 6 and parts[3] == "A":
            atomtypes.append(parts[0])
    if not atomtypes:
        raise ValueError(f"No atom types parsed from {path}")
    return atomtypes


def aggregate_depth_table(path: Path):
    per_type = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            bead_type = row["bead_type"].strip()
            entry = per_type.setdefault(
                bead_type,
                {
                    "count": 0,
                    "depth_sum": 0.0,
                    "hb_sum": 0.0,
                    "hb_donor_sum": 0.0,
                    "hb_acceptor_sum": 0.0,
                    "cb_sum": 0.0,
                },
            )
            entry["count"] += 1
            entry["depth_sum"] += float(row["mean_depth"])
            entry["hb_sum"] += -float(row["hb_mean"])
            entry["hb_donor_sum"] += -0.5 * (float(row["hb_donor_unbound"]) + float(row["hb_donor_bound"]))
            entry["hb_acceptor_sum"] += -0.5 * (float(row["hb_acceptor_unbound"]) + float(row["hb_acceptor_bound"]))
            entry["cb_sum"] += float(row["cb_backbone_mean"])
    if not per_type:
        raise ValueError(f"No rows parsed from {path}")
    out = {}
    for bead_type, entry in per_type.items():
        count = float(entry["count"])
        out[bead_type] = {
            "n_rows": int(entry["count"]),
            "mean_depth": entry["depth_sum"] / count,
            "mean_hb_attraction": entry["hb_sum"] / count,
            "mean_hb_donor_attraction": entry["hb_donor_sum"] / count,
            "mean_hb_acceptor_attraction": entry["hb_acceptor_sum"] / count,
            "mean_cb_value": entry["cb_sum"] / count,
        }
    return out


def pair_param(param_table, type1, type2):
    if (type1, type2) in param_table:
        return param_table[(type1, type2)]
    if (type2, type1) in param_table:
        return param_table[(type2, type1)]
    raise KeyError(f"Missing nonbonded parameter pair ({type1}, {type2})")


def load_depth_targets(depth_table_csv: Path, depth_table_meta: Path | None):
    depth_by_type = aggregate_depth_table(depth_table_csv)
    if depth_table_meta is None:
        return depth_by_type, {"source": "csv_only", "legacy_depth_by_type": depth_by_type}

    meta = json.loads(depth_table_meta.read_text(encoding="utf-8"))
    type_sample_stats = meta.get("type_sample_stats", {})
    if not isinstance(type_sample_stats, dict) or not type_sample_stats:
        return depth_by_type, {"source": "csv_only", "legacy_depth_by_type": depth_by_type}

    refined = {}
    comparison = {}
    for bead_type, legacy in depth_by_type.items():
        sample_stats = type_sample_stats.get(bead_type)
        if not isinstance(sample_stats, dict):
            refined[bead_type] = legacy
            continue
        refined_entry = {
            "n_rows": int(legacy.get("n_rows", 0)),
            "n_samples": int(sample_stats.get("n_samples", 0)),
            "mean_depth": float(sample_stats.get("mean_depth", legacy["mean_depth"])),
            "mean_hb_attraction": float(-sample_stats.get("hb_mean", -legacy["mean_hb_attraction"])),
            "mean_hb_donor_attraction": float(
                -0.5
                * (
                    sample_stats.get(
                        "hb_donor_unbound",
                        -legacy.get("mean_hb_donor_attraction", legacy["mean_hb_attraction"]),
                    )
                    + sample_stats.get(
                        "hb_donor_bound",
                        -legacy.get("mean_hb_donor_attraction", legacy["mean_hb_attraction"]),
                    )
                )
            ),
            "mean_hb_acceptor_attraction": float(
                -0.5
                * (
                    sample_stats.get(
                        "hb_acceptor_unbound",
                        -legacy.get("mean_hb_acceptor_attraction", legacy["mean_hb_attraction"]),
                    )
                    + sample_stats.get(
                        "hb_acceptor_bound",
                        -legacy.get("mean_hb_acceptor_attraction", legacy["mean_hb_attraction"]),
                    )
                )
            ),
            "mean_cb_value": float(sample_stats.get("cb_backbone_mean", legacy["mean_cb_value"])),
        }
        refined[bead_type] = refined_entry
        comparison[bead_type] = {
            "legacy_mean_depth": float(legacy["mean_depth"]),
            "refined_mean_depth": float(refined_entry["mean_depth"]),
            "legacy_mean_hb_attraction": float(legacy["mean_hb_attraction"]),
            "refined_mean_hb_attraction": float(refined_entry["mean_hb_attraction"]),
            "legacy_mean_hb_donor_attraction": float(legacy.get("mean_hb_donor_attraction", legacy["mean_hb_attraction"])),
            "refined_mean_hb_donor_attraction": float(refined_entry["mean_hb_donor_attraction"]),
            "legacy_mean_hb_acceptor_attraction": float(legacy.get("mean_hb_acceptor_attraction", legacy["mean_hb_attraction"])),
            "refined_mean_hb_acceptor_attraction": float(refined_entry["mean_hb_acceptor_attraction"]),
            "legacy_mean_cb_value": float(legacy["mean_cb_value"]),
            "refined_mean_cb_value": float(refined_entry["mean_cb_value"]),
        }
    return refined, {
        "source": "empirical_samples",
        "legacy_depth_by_type": depth_by_type,
        "legacy_vs_refined": comparison,
        "symmetry_report": meta.get("symmetry_report"),
    }


def compute_role_pair_descriptor(proxy_type, env_type, param_table):
    _sigma, epsilon = pair_param(param_table, proxy_type, env_type)
    return float(epsilon)


def merge_anchor_points_with_labels(x, y, labels):
    order = np.argsort(x)
    x_sorted = np.asarray(x, dtype=np.float64)[order]
    y_sorted = np.asarray(y, dtype=np.float64)[order]
    labels_sorted = [labels[i] for i in order.tolist()]
    merged_x = []
    merged_y = []
    merged_labels = []
    i = 0
    while i < x_sorted.shape[0]:
        j = i + 1
        values = [float(y_sorted[i])]
        group_labels = [labels_sorted[i]]
        while j < x_sorted.shape[0] and abs(x_sorted[j] - x_sorted[i]) < 1e-8:
            values.append(float(y_sorted[j]))
            group_labels.append(labels_sorted[j])
            j += 1
        merged_x.append(float(x_sorted[i]))
        merged_y.append(float(np.mean(np.asarray(values, dtype=np.float64))))
        merged_labels.append("+".join(group_labels))
        i = j
    if len(merged_x) < 2:
        raise ValueError("Need at least two unique anchor descriptors")
    return np.asarray(merged_x, dtype=np.float64), np.asarray(merged_y, dtype=np.float64), merged_labels


def build_anchor_model(depth_by_type, param_table, proxy_type, target_field):
    required = ("Q0", "Qa", "Na", "C1", "C3")
    missing = [key for key in required if key not in depth_by_type]
    if missing:
        raise ValueError(f"Depth table is missing required DOPC dry types: {missing}")
    raw_anchor_points = []
    for bead_type in DOPC_COVERED_TYPES:
        sigma_nm, epsilon_pair = pair_param(param_table, proxy_type, bead_type)
        raw_anchor_points.append(
            {
                "env_type": bead_type,
                "sigma_nm": float(sigma_nm),
                "descriptor": float(epsilon_pair),
                "target_value": float(depth_by_type[bead_type][target_field]),
            }
        )
    x = np.asarray([item["descriptor"] for item in raw_anchor_points], dtype=np.float64)
    y = np.asarray([item["target_value"] for item in raw_anchor_points], dtype=np.float64)
    labels = [item["env_type"] for item in raw_anchor_points]
    merged_x, merged_y, merged_labels = merge_anchor_points_with_labels(x, y, labels)
    return merged_x, merged_y, merged_labels, raw_anchor_points


def build_role_anchor_models(depth_by_type, param_table):
    anchor_models = {}
    for role in BACKBONE_CLASSES:
        target_field = ROLE_TARGET_FIELDS[role]
        proxy_type = ROLE_PROXIES[role]
        anchor_x, anchor_y, anchor_labels, raw_anchor_points = build_anchor_model(
            depth_by_type,
            param_table,
            proxy_type,
            target_field,
        )
        anchor_models[role] = {
            "proxy_type": proxy_type,
            "target_field": target_field,
            "target_description": ROLE_TARGET_DESCRIPTIONS[role],
            "descriptor_name": ROLE_DESCRIPTOR_NAME,
            "rmin_name": ROLE_RMIN_NAME,
            "x": anchor_x,
            "y": anchor_y,
            "labels": anchor_labels,
            "raw_anchor_points": raw_anchor_points,
        }
    return anchor_models


def pchip_slopes(x, y):
    n = x.shape[0]
    h = np.diff(x)
    delta = np.diff(y) / h
    m = np.zeros(n, dtype=np.float64)
    if n == 2:
        m[:] = delta[0]
        return m
    for k in range(1, n - 1):
        if delta[k - 1] == 0.0 or delta[k] == 0.0 or np.sign(delta[k - 1]) != np.sign(delta[k]):
            m[k] = 0.0
        else:
            w1 = 2.0 * h[k] + h[k - 1]
            w2 = h[k] + 2.0 * h[k - 1]
            m[k] = (w1 + w2) / ((w1 / delta[k - 1]) + (w2 / delta[k]))
    m0 = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    if np.sign(m0) != np.sign(delta[0]):
        m0 = 0.0
    elif np.sign(delta[0]) != np.sign(delta[1]) and abs(m0) > 3.0 * abs(delta[0]):
        m0 = 3.0 * delta[0]
    m[0] = m0
    mn = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    if np.sign(mn) != np.sign(delta[-1]):
        mn = 0.0
    elif np.sign(delta[-1]) != np.sign(delta[-2]) and abs(mn) > 3.0 * abs(delta[-1]):
        mn = 3.0 * delta[-1]
    m[-1] = mn
    return m


def pchip_interpolate(x, y, xq):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xq = np.asarray(xq, dtype=np.float64)
    m = pchip_slopes(x, y)
    idx = np.searchsorted(x, xq, side="right") - 1
    idx = np.clip(idx, 0, x.shape[0] - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    m0 = m[idx]
    m1 = m[idx + 1]
    h = x1 - x0
    t = (xq - x0) / h
    h00 = (2.0 * t**3) - (3.0 * t**2) + 1.0
    h10 = (t**3) - (2.0 * t**2) + t
    h01 = (-2.0 * t**3) + (3.0 * t**2)
    h11 = (t**3) - (t**2)
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def interpolate_depth(epsilon_value, anchor_x, anchor_y):
    epsilon_value = float(epsilon_value)
    if epsilon_value <= anchor_x[0]:
        return float(anchor_y[0])
    if epsilon_value <= anchor_x[-1]:
        return float(pchip_interpolate(anchor_x, anchor_y, np.array([epsilon_value]))[0])
    slope = (anchor_y[-1] - anchor_y[-2]) / (anchor_x[-1] - anchor_x[-2])
    value = anchor_y[-1] + slope * (epsilon_value - anchor_x[-1])
    return float(min(value, 1.25 * anchor_y[-1]))


def merge_control_points(x, y):
    order = np.argsort(x)
    x_sorted = np.asarray(x, dtype=np.float64)[order]
    y_sorted = np.asarray(y, dtype=np.float64)[order]
    merged_x = []
    merged_y = []
    i = 0
    while i < x_sorted.shape[0]:
        j = i + 1
        values = [y_sorted[i]]
        while j < x_sorted.shape[0] and abs(x_sorted[j] - x_sorted[i]) < 1e-8:
            values.append(y_sorted[j])
            j += 1
        merged_x.append(float(x_sorted[i]))
        merged_y.append(float(np.mean(np.asarray(values, dtype=np.float64))))
        i = j
    if len(merged_x) < 2:
        raise ValueError("Need at least two unique control points")
    return np.asarray(merged_x, dtype=np.float64), np.asarray(merged_y, dtype=np.float64)


def radial_control_points(r_min, depth, radial_min, radial_max):
    x = np.asarray(
        [
            max(radial_min, 0.72 * r_min),
            max(radial_min + 0.2, 0.90 * r_min),
            r_min,
            min(radial_max, 1.35 * r_min),
            radial_max,
        ],
        dtype=np.float64,
    )
    y = np.asarray([3.0 * depth, 1.2 * depth, -depth, -0.30 * depth, 0.0], dtype=np.float64)
    return merge_control_points(x, y)


def gaussian_basis(r, centers, widths):
    delta = (r[:, None] - centers[None, :]) / widths[None, :]
    return np.exp(-0.5 * delta * delta)


def fit_gaussian_weights(target_values, centers, widths, radial_grid, ridge):
    basis = gaussian_basis(radial_grid, centers, widths)
    gram = basis.T @ basis
    rhs = basis.T @ target_values
    reg = ridge * np.eye(centers.shape[0], dtype=np.float64)
    weights = np.linalg.solve(gram + reg, rhs)
    fitted = basis @ weights
    rmse = float(np.sqrt(np.mean((fitted - target_values) ** 2)))
    return weights, fitted, rmse


def decode_ring_scale(atom_type):
    return 0.75 if atom_type.startswith("S") else 1.0


def build_cross_type_records(atomtypes, param_table, anchor_models_by_role, radial_grid):
    rows = []
    table_tensor = np.zeros((len(BACKBONE_CLASSES), len(atomtypes), len(radial_grid)), dtype=np.float32)
    pair_summaries = {}
    for env_index, env_type in enumerate(atomtypes):
        for role_index, role in enumerate(BACKBONE_CLASSES):
            anchor_model = anchor_models_by_role[role]
            proxy_type = anchor_model["proxy_type"]
            sigma_base, epsilon_base = pair_param(param_table, proxy_type, env_type)
            descriptor = compute_role_pair_descriptor(proxy_type, env_type, param_table)
            r_min = float((2.0 ** (1.0 / 6.0)) * sigma_base * 10.0)
            base_depth = interpolate_depth(descriptor, anchor_model["x"], anchor_model["y"])
            ctrl_x, ctrl_y = radial_control_points(r_min, base_depth, float(radial_grid[0]), float(radial_grid[-1]))
            target = pchip_interpolate(ctrl_x, ctrl_y, radial_grid)
            well_depth = float(base_depth)
            table_tensor[role_index, env_index, :] = target.astype(np.float32)
            pair_summaries[f"{role}:{env_type}"] = {
                "proxy_type": proxy_type,
                "r_min_angstrom": r_min,
                "well_depth": well_depth,
                "base_depth": float(base_depth),
                "role_scale": 1.0,
                "target_field": anchor_model["target_field"],
                "target_description": anchor_model["target_description"],
                "descriptor_name": anchor_model["descriptor_name"],
                "descriptor_value": float(descriptor),
                "sigma_nm": float(sigma_base),
                "epsilon_pair": float(epsilon_base),
                "fit_rmse": 0.0,
                "control_r": [float(x) for x in ctrl_x],
                "control_v": [float(y) for y in ctrl_y],
            }
            for r_value, target_value in zip(radial_grid, target):
                rows.append(
                    {
                        "protein_class": role,
                        "env_type": env_type,
                        "proxy_type": proxy_type,
                        "r_angstrom": float(r_value),
                        "target_potential": float(target_value),
                        "fitted_potential": float(target_value),
                        "r_min_angstrom": r_min,
                        "well_depth": well_depth,
                        "base_depth": float(base_depth),
                        "role_scale": 1.0,
                        "sigma_nm": float(sigma_base),
                        "epsilon_pair": float(epsilon_base),
                        "descriptor_value": float(descriptor),
                        "fit_rmse": 0.0,
                    }
                )
    return rows, table_tensor, pair_summaries


def build_overlap_report(depth_by_type, pair_summaries):
    overlap_rows = []
    predicted = []
    observed = []
    role_correlations = {}
    role_predicted = {role: [] for role in BACKBONE_CLASSES}
    role_observed = {role: [] for role in BACKBONE_CLASSES}
    for env_type in DOPC_COVERED_TYPES:
        role_depths = [pair_summaries[f"{role}:{env_type}"]["well_depth"] for role in BACKBONE_CLASSES]
        predicted_value = float(np.mean(np.asarray(role_depths, dtype=np.float64)))
        observed_value = float(depth_by_type[env_type]["mean_hb_attraction"])
        row = {
            "env_type": env_type,
            "observed_hb_attraction": observed_value,
            "predicted_role_mean_well_depth": predicted_value,
            "mean_depth": float(depth_by_type[env_type]["mean_depth"]),
            "n_rows": int(depth_by_type[env_type]["n_rows"]),
            "n_samples": int(depth_by_type[env_type].get("n_samples", depth_by_type[env_type]["n_rows"])),
        }
        for role in BACKBONE_CLASSES:
            target_field = ROLE_TARGET_FIELDS[role]
            role_pred = float(pair_summaries[f"{role}:{env_type}"]["well_depth"])
            role_obs = float(depth_by_type[env_type][target_field])
            row[f"observed_{role}_target"] = role_obs
            row[f"predicted_{role}_well_depth"] = role_pred
            role_predicted[role].append(role_pred)
            role_observed[role].append(role_obs)
        overlap_rows.append(row)
        predicted.append(predicted_value)
        observed.append(observed_value)
    predicted_arr = np.asarray(predicted, dtype=np.float64)
    observed_arr = np.asarray(observed, dtype=np.float64)
    corr = float(np.corrcoef(predicted_arr, observed_arr)[0, 1]) if predicted_arr.size >= 2 else 1.0
    for role in BACKBONE_CLASSES:
        pred_arr = np.asarray(role_predicted[role], dtype=np.float64)
        obs_arr = np.asarray(role_observed[role], dtype=np.float64)
        role_correlations[role] = float(np.corrcoef(pred_arr, obs_arr)[0, 1]) if pred_arr.size >= 2 else 1.0
    return {"covered_types": overlap_rows, "pearson_r": corr, "role_pearson_r": role_correlations}


def write_csv_rows(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_cross_artifact(path: Path, protein_classes, env_types, radial_grid, table_values, metadata):
    path.parent.mkdir(parents=True, exist_ok=True)
    str_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["schema"] = SPLINE_CROSS_SCHEMA
        meta.attrs["calibration_mode"] = metadata["calibration"]["calibration_mode"]
        meta.attrs["descriptor_name"] = metadata["calibration"]["descriptor_name"]
        meta.attrs["radial_min"] = np.float32(metadata["radial_grid"]["min"])
        meta.attrs["radial_max"] = np.float32(metadata["radial_grid"]["max"])
        meta.attrs["sample_step"] = np.float32(metadata["radial_grid"]["step"])
        meta.attrs["n_env_type"] = np.int32(len(env_types))
        meta.attrs["n_protein_class"] = np.int32(len(protein_classes))
        for role in protein_classes:
            meta.attrs[f"proxy_{role}"] = metadata["role_proxies"][role]

        sources = h5.create_group("sources")
        sources.attrs["ff_itp"] = metadata["sources"]["ff_itp"]
        sources.attrs["depth_table_csv"] = metadata["sources"]["depth_table_csv"]

        classes = h5.create_group("classes")
        classes.create_dataset("protein", data=np.asarray(protein_classes, dtype=str_dtype))
        classes.create_dataset("environment", data=np.asarray(env_types, dtype=str_dtype))

        cross = h5.create_group("spline").create_group("cross")
        cross.attrs["cutoff"] = np.float32(metadata["radial_grid"]["max"])
        cross.attrs["radial_min"] = np.float32(metadata["radial_grid"]["min"])
        cross.attrs["radial_step"] = np.float32(metadata["radial_grid"]["step"])
        cross.create_dataset("radial_grid", data=np.asarray(radial_grid, dtype=np.float32))
        cross.create_dataset("table_values", data=np.asarray(table_values, dtype=np.float32))


def write_dry_nonbond_forcefield_artifact(path: Path, ff_itp: Path):
    param_table = read_martini3_nonbond_params(str(ff_itp))
    if not param_table:
        raise ValueError(f"Failed to read dry MARTINI nonbond parameters from {ff_itp}")
    type_order = None
    with h5py.File(path, "r") as h5:
        if "/classes/environment" in h5:
            type_order = decode_bytes_array(h5["/classes/environment"][:]).tolist()
    write_dry_nonbond_artifact(path, param_table, ff_itp, type_order=type_order)


def load_manifest_paths(manifest_path: Path):
    with manifest_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        if isinstance(data.get("systems"), list):
            items = data["systems"]
        elif isinstance(data.get("up_files"), list):
            items = data["up_files"]
        else:
            raise ValueError(f"Unsupported manifest schema: {manifest_path}")
    elif isinstance(data, list):
        items = data
    else:
        raise ValueError(f"Unsupported manifest top-level type: {type(data)}")
    paths = []
    for item in items:
        raw_path = item if isinstance(item, str) else item.get("up_file") or item.get("upside_input")
        if not raw_path:
            continue
        candidate = Path(raw_path).expanduser()
        candidate = (manifest_path.parent / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
        paths.append(candidate)
    return paths


def collect_cross_up_paths(manifest_path, up_files):
    paths = []
    if manifest_path is not None:
        manifest = require_existing_file(manifest_path)
        paths.extend(load_manifest_paths(manifest))
    for raw in up_files:
        paths.append(Path(raw).expanduser().resolve())
    unique = []
    seen = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        unique.append(path)
    if not unique:
        raise ValueError("Provide at least one .up file via --up-files and/or --manifest")
    missing = [path for path in unique if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing .up files:\n" + "\n".join(str(path) for path in missing))
    return unique


def load_cross_artifact(path: Path):
    with h5py.File(path, "r") as h5:
        schema = str(h5["/meta"].attrs.get("schema", ""))
        protein_classes = decode_bytes_array(h5["/classes/protein"][:]).tolist()
        env_classes = decode_bytes_array(h5["/classes/environment"][:]).tolist()
        if schema == RBM_CROSS_SCHEMA:
            centers = h5["/rbm/cross/radial_centers"][:].astype(np.float32)
            widths = h5["/rbm/cross/radial_widths"][:].astype(np.float32)
            weights = h5["/rbm/cross/weights"][:].astype(np.float32)
            cutoff = float(h5["/rbm/cross"].attrs["cutoff"])
            return {
                "schema": schema,
                "protein_classes": protein_classes,
                "environment_classes": env_classes,
                "radial_centers": centers,
                "radial_widths": widths,
                "weights": weights,
                "cutoff": cutoff,
            }
        if schema == SPLINE_CROSS_SCHEMA:
            cross = h5["/spline/cross"]
            return {
                "schema": schema,
                "calibration_mode": str(h5["/meta"].attrs.get("calibration_mode", "")),
                "descriptor_name": str(h5["/meta"].attrs.get("descriptor_name", "")),
                "protein_classes": protein_classes,
                "environment_classes": env_classes,
                "radial_grid": cross["radial_grid"][:].astype(np.float32),
                "table_values": cross["table_values"][:].astype(np.float32),
                "cutoff": float(cross.attrs["cutoff"]),
                "radial_min": float(cross.attrs["radial_min"]),
                "radial_step": float(cross.attrs["radial_step"]),
            }
    raise ValueError(f"Unsupported cross artifact schema in {path}: {schema}")


def collect_cross_indices(up_path: Path, artifact):
    protein_class_to_idx = {name: i for i, name in enumerate(artifact["protein_classes"])}
    env_class_to_idx = {name: i for i, name in enumerate(artifact["environment_classes"])}
    with h5py.File(up_path, "r") as h5:
        roles = decode_bytes_array(h5["/input/atom_roles"][:])
        types = decode_bytes_array(h5["/input/type"][:])
        membership = h5["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    if roles.shape[0] != types.shape[0] or roles.shape[0] != membership.shape[0]:
        raise ValueError(f"{up_path}: role/type/membership length mismatch")
    protein_mask = membership >= 0
    env_mask = ~protein_mask
    protein_indices = []
    protein_class_indices = []
    for idx, role in enumerate(roles.tolist()):
        role_name = role.strip()
        if protein_mask[idx] and role_name in protein_class_to_idx:
            protein_indices.append(idx)
            protein_class_indices.append(protein_class_to_idx[role_name])
    if not protein_indices:
        raise ValueError(f"{up_path}: no explicit N/CA/C/O backbone carriers found for injection")
    env_indices = []
    env_class_indices = []
    missing_types = set()
    for idx, atom_type in enumerate(types.tolist()):
        type_name = atom_type.strip()
        if not env_mask[idx]:
            continue
        if type_name not in env_class_to_idx:
            missing_types.add(type_name)
            continue
        env_indices.append(idx)
        env_class_indices.append(env_class_to_idx[type_name])
    if missing_types:
        raise ValueError(f"{up_path}: artifact is missing environment classes for types {sorted(missing_types)}")
    if not env_indices:
        raise ValueError(f"{up_path}: no environment atoms found for injection")
    return (
        np.asarray(protein_indices, dtype=np.int32),
        np.asarray(protein_class_indices, dtype=np.int32),
        np.asarray(env_indices, dtype=np.int32),
        np.asarray(env_class_indices, dtype=np.int32),
    )


def inject_cross_node(up_path, artifact, protein_indices, protein_class_indices, env_indices, env_class_indices, overwrite_existing):
    with h5py.File(up_path, "r+") as h5:
        pot_grp = h5.require_group("input").require_group("potential")
        if artifact["schema"] == RBM_CROSS_SCHEMA:
            node_name = RBM_CROSS_NODE_NAME
            other_node_name = SPLINE_CROSS_NODE_NAME
        elif artifact["schema"] == SPLINE_CROSS_SCHEMA:
            node_name = SPLINE_CROSS_NODE_NAME
            other_node_name = RBM_CROSS_NODE_NAME
        else:
            raise ValueError(f"Unsupported cross artifact schema for injection: {artifact['schema']}")
        if node_name in pot_grp:
            if not overwrite_existing:
                raise ValueError(f"{up_path}: {node_name} already exists and overwrite is disabled")
            del pot_grp[node_name]
        if other_node_name in pot_grp:
            if not overwrite_existing:
                raise ValueError(
                    f"{up_path}: {other_node_name} already exists and overwrite is disabled; "
                    f"remove the legacy cross node before injecting {node_name}"
                )
            del pot_grp[other_node_name]
        used_protein_classes = sorted(set(int(x) for x in protein_class_indices.tolist()))
        used_env_classes = sorted(set(int(x) for x in env_class_indices.tolist()))
        protein_reindex = {old: new for new, old in enumerate(used_protein_classes)}
        env_reindex = {old: new for new, old in enumerate(used_env_classes)}
        local_protein_class_indices = np.asarray(
            [protein_reindex[int(x)] for x in protein_class_indices.tolist()],
            dtype=np.int32,
        )
        local_env_class_indices = np.asarray(
            [env_reindex[int(x)] for x in env_class_indices.tolist()],
            dtype=np.int32,
        )
        grp = pot_grp.create_group(node_name)
        grp.attrs["arguments"] = np.asarray([b"pos"])
        grp.attrs["initialized"] = np.int8(1)
        grp.attrs["cutoff"] = np.float32(artifact["cutoff"])
        grp.create_dataset("protein_atom_indices", data=protein_indices)
        grp.create_dataset("protein_class_index", data=local_protein_class_indices)
        grp.create_dataset("env_atom_indices", data=env_indices)
        grp.create_dataset("env_class_index", data=local_env_class_indices)
        if artifact["schema"] == RBM_CROSS_SCHEMA:
            local_weights = artifact["weights"][used_protein_classes][:, used_env_classes, :].astype(np.float32, copy=False)
            grp.create_dataset("radial_centers", data=artifact["radial_centers"])
            grp.create_dataset("radial_widths", data=artifact["radial_widths"])
            grp.create_dataset("weights", data=local_weights)
        else:
            local_tables = artifact["table_values"][used_protein_classes][:, used_env_classes, :].astype(np.float32, copy=False)
            grp.attrs["radial_min"] = np.float32(artifact["radial_min"])
            grp.attrs["radial_step"] = np.float32(artifact["radial_step"])
            grp.create_dataset("radial_grid", data=artifact["radial_grid"])
            grp.create_dataset("table_values", data=local_tables)
        grp.create_dataset(
            "protein_classes",
            data=np.asarray([artifact["protein_classes"][idx] for idx in used_protein_classes], dtype=h5py.string_dtype("utf-8")),
        )
        grp.create_dataset(
            "environment_classes",
            data=np.asarray([artifact["environment_classes"][idx] for idx in used_env_classes], dtype=h5py.string_dtype("utf-8")),
        )


def run_handoff_command(argv):
    parser = argparse.ArgumentParser(description="Copy the last frame from one stage into the next stage input.")
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args(argv)
    set_initial_position(args.input_file, args.output_file)


def run_inject_backbone_only_command(argv):
    parser = argparse.ArgumentParser(description="Inject Upside backbone-only nodes into a stage .up file.")
    parser.add_argument("up_file")
    parser.add_argument("protein_source")
    parser.add_argument("upside_home")
    parser.add_argument("rama_library")
    parser.add_argument("rama_sheet_mixing")
    parser.add_argument("hbond_energy")
    parser.add_argument("reference_state_rama")
    args = parser.parse_args(argv)
    require_existing_file(args.up_file)
    require_existing_file(args.protein_source)
    require_existing_file(args.rama_library)
    require_existing_file(args.rama_sheet_mixing)
    require_existing_file(args.hbond_energy)
    require_existing_file(args.reference_state_rama)
    upside_config = load_upside_config(args.upside_home)
    residue_ids, residue_to_ncac = collect_runtime_ncac_mapping(args.up_file)
    residue_names = parse_backbone_source_residue_names(args.protein_source)
    ref_n_atom = write_backbone_nodes(
        upside_config=upside_config,
        up_file=args.up_file,
        residue_ids=residue_ids,
        residue_to_ncac=residue_to_ncac,
        residue_names=residue_names,
        rama_library=args.rama_library,
        rama_sheet_mixing=args.rama_sheet_mixing,
        hbond_energy=args.hbond_energy,
        reference_state_rama=args.reference_state_rama,
    )
    print(
        f"Injected backbone-only Upside nodes into {args.up_file}: "
        f"n_res={len(residue_ids)} ref_backbone_atoms={ref_n_atom}"
    )


def run_validate_backbone_only_command(argv):
    parser = argparse.ArgumentParser(description="Validate backbone-only stage files and table artifacts.")
    parser.add_argument("up_file", nargs="?", default=None, type=Path)
    parser.add_argument("--table-csv", type=Path, default=None)
    parser.add_argument("--table-meta", type=Path, default=None)
    parser.add_argument("--cross-table-csv", type=Path, default=None)
    parser.add_argument("--cross-table-meta", type=Path, default=None)
    parser.add_argument("--cross-artifact", type=Path, default=None)
    args = parser.parse_args(argv)
    if (
        args.up_file is None
        and args.table_csv is None
        and args.table_meta is None
        and args.cross_table_csv is None
        and args.cross_table_meta is None
        and args.cross_artifact is None
    ):
        raise ValueError("Provide at least one artifact to validate")
    up_summary = validate_up_file(require_existing_file(args.up_file)) if args.up_file is not None else None
    table_summary = validate_table_csv(require_existing_file(args.table_csv)) if args.table_csv is not None else None
    meta_summary = validate_table_meta(require_existing_file(args.table_meta)) if args.table_meta is not None else None
    cross_table_summary = (
        validate_cross_table_csv(require_existing_file(args.cross_table_csv))
        if args.cross_table_csv is not None
        else None
    )
    cross_meta_summary = (
        validate_cross_table_meta(require_existing_file(args.cross_table_meta))
        if args.cross_table_meta is not None
        else None
    )
    cross_artifact_summary = (
        validate_cross_artifact(require_existing_file(args.cross_artifact))
        if args.cross_artifact is not None
        else None
    )
    if up_summary is not None:
        print(f"OK: {args.up_file}")
        print(f"  backbone_only_stage n_atom={up_summary['n_atom']} n_res={up_summary['n_res']}")
    else:
        print("OK: table artifacts")
    if table_summary is not None:
        print(f"  depth_table rows={table_summary['n_rows']}")
    if meta_summary is not None:
        tp = meta_summary["term_presence"]
        print(
            "  membrane_terms "
            f"backbone_cb_energy={int(bool(tp.get('backbone_cb_energy', 0)))} "
            f"hbond_hb_energy={int(bool(tp.get('hbond_hb_energy', 0)))} "
            f"sidechain_term={int(bool(tp.get('sidechain_term', 0)))} "
            f"sheet_term={int(bool(tp.get('sheet_term', 0)))}"
        )
    if cross_table_summary is not None:
        print(
            "  cross_table "
            f"rows={cross_table_summary['n_rows']} "
            f"protein_classes={cross_table_summary['n_protein_class']} "
            f"env_types={cross_table_summary['n_env_type']}"
        )
    if cross_meta_summary is not None:
        rmse = cross_meta_summary["fit_rmse_summary"]
        print(
            "  cross_meta "
            f"env_types={cross_meta_summary['n_env_type']} "
            f"calibration_mode={cross_meta_summary['calibration_mode']} "
            f"descriptor_name={cross_meta_summary['descriptor_name']} "
            f"fit_rmse_min={float(rmse['min']):.6f} "
            f"fit_rmse_median={float(rmse['median']):.6f} "
            f"fit_rmse_max={float(rmse['max']):.6f}"
        )
    if cross_artifact_summary is not None:
        print(
            "  cross_artifact "
            f"schema={cross_artifact_summary['schema']} "
            f"calibration_mode={cross_artifact_summary.get('calibration_mode', 'n/a')} "
            f"descriptor_name={cross_artifact_summary.get('descriptor_name', 'n/a')} "
            f"dry_schema={cross_artifact_summary.get('dry_schema', 'n/a')} "
            f"dry_types={cross_artifact_summary.get('n_dry_type', 'n/a')} "
            f"env_types={cross_artifact_summary['n_env_type']} "
            f"n_radial={cross_artifact_summary['n_radial']}"
        )


def run_validate_stage_potentials_command(argv):
    parser = argparse.ArgumentParser(
        description="Validate that a stage file does not contain forbidden membrane runtime nodes."
    )
    parser.add_argument("up_file", type=Path)
    args = parser.parse_args(argv)
    summary = validate_stage_potential_nodes(require_existing_file(args.up_file))
    print(f"OK: {args.up_file}")
    print(
        "  stage_potentials "
        f"n_nodes={summary['n_potential_nodes']} "
        "forbidden_membrane_nodes=0"
    )


def run_build_depth_table_command(argv):
    parser = argparse.ArgumentParser(description="Build a depth-dependent dry-MARTINI bead interaction table.")
    parser.add_argument("--bilayer-pdb", type=Path, default=Path("pdb/bilayer.MARTINI.pdb"))
    parser.add_argument("--lipid-itp", type=Path, default=Path("ff_dry/dry_martini_v2.1_lipids.itp"))
    parser.add_argument("--membrane-h5", type=Path, default=Path("../../parameters/ff_2.1/membrane.h5"))
    parser.add_argument("--output-csv", type=Path, default=Path("outputs/depth_interaction_table.csv"))
    parser.add_argument("--output-json", type=Path, default=Path("outputs/depth_interaction_table.meta.json"))
    parser.add_argument("--lipid-resnames", nargs="+", default=["DOPC", "DOP"])
    args = parser.parse_args(argv)
    bilayer_pdb = require_existing_file(args.bilayer_pdb)
    lipid_itp = require_existing_file(args.lipid_itp)
    membrane_h5 = require_existing_file(args.membrane_h5)
    atom_to_type = parse_lipid_itp_atom_types(lipid_itp, args.lipid_resnames)
    z_center, bead_depth_stats, bead_signed_samples = parse_bilayer_depths(bilayer_pdb, args.lipid_resnames)
    cb_z, hb_z, cb_energy, hb_energy, cb_names, term_presence = load_membrane_channels(membrane_h5)
    bead_sample_stats = summarize_depth_samples(bead_signed_samples, cb_z, hb_z, cb_energy, hb_energy)
    type_signed_samples = pool_depth_samples_by_type(atom_to_type, bead_signed_samples)
    type_sample_stats = summarize_depth_samples(type_signed_samples, cb_z, hb_z, cb_energy, hb_energy)
    symmetry_report = audit_membrane_symmetry(cb_z, hb_z, cb_energy, hb_energy, cb_names)
    rows = []
    for bead_name, stats in bead_depth_stats.items():
        sample_stats = bead_sample_stats[bead_name]
        rows.append(
            {
                "bead_name": bead_name,
                "bead_type": atom_to_type.get(bead_name, "UNKNOWN"),
                "mean_depth": stats["mean_depth"],
                "std_depth": stats["std_depth"],
                "count": stats["count"],
                "cb_backbone_burial0_mean": float(sample_stats["cb_backbone_burial0_mean"]),
                "cb_backbone_burial1_mean": float(sample_stats["cb_backbone_burial1_mean"]),
                "cb_backbone_mean": float(sample_stats["cb_backbone_mean"]),
                "hb_donor_unbound": float(sample_stats["hb_donor_unbound"]),
                "hb_donor_bound": float(sample_stats["hb_donor_bound"]),
                "hb_acceptor_unbound": float(sample_stats["hb_acceptor_unbound"]),
                "hb_acceptor_bound": float(sample_stats["hb_acceptor_bound"]),
                "hb_mean": float(sample_stats["hb_mean"]),
            }
        )
    rows.sort(key=lambda row: row["mean_depth"], reverse=True)
    fieldnames = list(REQUIRED_DEPTH_TABLE_COLUMNS)
    write_csv_rows(args.output_csv, fieldnames, rows)
    metadata = {
        "bilayer_pdb": str(bilayer_pdb),
        "lipid_itp": str(lipid_itp),
        "membrane_h5": str(membrane_h5),
        "lipid_resnames": [x.upper() for x in args.lipid_resnames],
        "z_center": z_center,
        "n_bead_names": len(rows),
        "cb_residue_names": cb_names,
        "depth_definition": "|z - z_center| (no PO4 anchor)",
        "aggregation_mode": "empirical_depth_samples",
        "term_presence": term_presence,
        "bead_sample_stats": bead_sample_stats,
        "type_sample_stats": type_sample_stats,
        "bead_signed_depth_samples": {
            key: [float(x) for x in np.asarray(values, dtype=np.float64).tolist()]
            for key, values in bead_signed_samples.items()
        },
        "type_signed_depth_samples": {
            key: [float(x) for x in np.asarray(values, dtype=np.float64).tolist()]
            for key, values in type_signed_samples.items()
        },
        "symmetry_report": symmetry_report,
    }
    write_summary(args.output_json, metadata)
    print(f"Wrote depth interaction table: {args.output_csv}")
    print(f"Wrote metadata: {args.output_json}")
    print(f"Bilayer z_center = {z_center:.6f} A")
    for key, value in term_presence.items():
        print(f"  {key}={int(bool(value))}")


def run_build_backbone_cross_table_command(argv):
    parser = argparse.ArgumentParser(description="Build the backbone cross interaction table artifact.")
    parser.add_argument("--ff-itp", type=Path, default=SCRIPT_DIR / "ff_dry" / "dry_martini_v2.1.itp")
    parser.add_argument("--depth-table-csv", type=Path, default=SCRIPT_DIR / "outputs" / "depth_interaction_table.csv")
    parser.add_argument("--depth-table-meta", type=Path, default=None)
    parser.add_argument("--output-csv", type=Path, default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.csv")
    parser.add_argument("--output-json", type=Path, default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.meta.json")
    parser.add_argument("--output-h5", type=Path, default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.h5")
    parser.add_argument("--radial-min", type=float, default=2.0)
    parser.add_argument("--radial-max", type=float, default=12.0)
    parser.add_argument("--sample-step", type=float, default=0.1)
    parser.add_argument("--n-radial", type=int, default=12)
    parser.add_argument("--ridge", type=float, default=1e-4)
    args = parser.parse_args(argv)
    ff_itp = require_existing_file(args.ff_itp)
    depth_table_csv = require_existing_file(args.depth_table_csv)
    depth_table_meta = require_existing_file(args.depth_table_meta) if args.depth_table_meta is not None else None
    if args.radial_min <= 0.0 or args.radial_max <= args.radial_min:
        raise ValueError("Require 0 < radial_min < radial_max")
    if args.sample_step <= 0.0 or args.n_radial <= 0 or args.ridge < 0.0:
        raise ValueError("Invalid radial basis settings")
    repo_root = Path(__file__).resolve().parents[2]
    atomtypes = parse_dry_atomtypes(ff_itp)
    param_table = read_martini3_nonbond_params(str(ff_itp))
    depth_by_type, depth_target_info = load_depth_targets(depth_table_csv, depth_table_meta)
    descriptor_by_role = {
        role: {
            env_type: compute_role_pair_descriptor(ROLE_PROXIES[role], env_type, param_table)
            for env_type in atomtypes
        }
        for role in BACKBONE_CLASSES
    }
    anchor_models_by_role = build_role_anchor_models(depth_by_type, param_table)
    radial_grid = np.arange(args.radial_min, args.radial_max + 0.5 * args.sample_step, args.sample_step)
    if radial_grid[-1] > args.radial_max:
        radial_grid[-1] = args.radial_max
    rows, table_tensor, pair_summaries = build_cross_type_records(
        atomtypes=atomtypes,
        param_table=param_table,
        anchor_models_by_role=anchor_models_by_role,
        radial_grid=radial_grid,
    )
    overlap_report = build_overlap_report(depth_by_type, pair_summaries)
    fit_rmses = np.asarray([item["fit_rmse"] for item in pair_summaries.values()], dtype=np.float64)
    legacy_vs_refined = depth_target_info.get("legacy_vs_refined", {})
    anchor_targets_by_role = {}
    for role in BACKBONE_CLASSES:
        anchor_model = anchor_models_by_role[role]
        role_targets = {}
        for anchor_point in anchor_model["raw_anchor_points"]:
            label = anchor_point["env_type"]
            role_targets[label] = {
                "proxy_type": anchor_model["proxy_type"],
                "descriptor": float(anchor_point["descriptor"]),
                "sigma_nm": float(anchor_point["sigma_nm"]),
                "refined_well_depth": float(anchor_point["target_value"]),
                "target_field": anchor_model["target_field"],
                "target_description": anchor_model["target_description"],
            }
            if label in legacy_vs_refined:
                legacy_key = anchor_model["target_field"]
                role_targets[label][f"legacy_{legacy_key}"] = float(
                    legacy_vs_refined[label].get(legacy_key.replace("mean_", "legacy_mean_"), legacy_vs_refined[label]["legacy_mean_hb_attraction"])
                )
                role_targets[label][f"refined_{legacy_key}"] = float(
                    legacy_vs_refined[label].get(legacy_key.replace("mean_", "refined_mean_"), role_targets[label]["refined_well_depth"])
                )
                role_targets[label][f"delta_{legacy_key}"] = (
                    role_targets[label][f"refined_{legacy_key}"] - role_targets[label][f"legacy_{legacy_key}"]
                )
        anchor_targets_by_role[role] = role_targets
    metadata = {
        "artifact_schema": SPLINE_CROSS_SCHEMA,
        "sources": {
            "ff_itp": repo_relative(ff_itp, repo_root),
            "depth_table_csv": repo_relative(depth_table_csv, repo_root),
            "depth_table_meta": repo_relative(depth_table_meta, repo_root) if depth_table_meta is not None else None,
        },
        "protein_classes": list(BACKBONE_CLASSES),
        "environment_types": atomtypes,
        "radial_grid": {
            "min": float(args.radial_min),
            "max": float(args.radial_max),
            "step": float(args.sample_step),
            "n_samples": int(radial_grid.shape[0]),
            "representation": "direct_spline_table",
            "legacy_rbm_settings": {
                "n_radial_basis": int(args.n_radial),
                "ridge": float(args.ridge),
            },
        },
        "calibration": {
            "calibration_mode": CROSS_CALIBRATION_MODE,
            "descriptor_name": ROLE_DESCRIPTOR_NAME,
            "rmin_name": ROLE_RMIN_NAME,
            "role_target_fields": ROLE_TARGET_FIELDS,
            "role_target_descriptions": ROLE_TARGET_DESCRIPTIONS,
            "anchor_models_by_role": {
                role: {
                    "proxy_type": anchor_models_by_role[role]["proxy_type"],
                    "target_field": anchor_models_by_role[role]["target_field"],
                    "target_description": anchor_models_by_role[role]["target_description"],
                    "descriptor_name": anchor_models_by_role[role]["descriptor_name"],
                    "rmin_name": anchor_models_by_role[role]["rmin_name"],
                    "merged_anchor_labels": anchor_models_by_role[role]["labels"],
                    "merged_anchor_descriptors": [float(x) for x in anchor_models_by_role[role]["x"]],
                    "merged_anchor_depths": [float(y) for y in anchor_models_by_role[role]["y"]],
                    "raw_anchor_points": anchor_models_by_role[role]["raw_anchor_points"],
                }
                for role in BACKBONE_CLASSES
            },
            "covered_depth_types": {key: depth_by_type[key] for key in DOPC_COVERED_TYPES if key in depth_by_type},
            "anchor_targets_by_role": anchor_targets_by_role,
            "legacy_vs_refined": legacy_vs_refined,
        },
        "role_proxies": ROLE_PROXIES,
        "role_scaling": CROSS_CALIBRATION_MODE,
        "type_descriptor_by_role": descriptor_by_role,
        "pair_summary": pair_summaries,
        "fit_rmse_summary": {
            "min": float(np.min(fit_rmses)),
            "median": float(np.median(fit_rmses)),
            "max": float(np.max(fit_rmses)),
        },
        "overlap_report": overlap_report,
        "symmetry_report": depth_target_info.get("symmetry_report"),
        "depth_target_source": depth_target_info.get("source", "csv_only"),
    }
    write_csv_rows(args.output_csv, REQUIRED_CROSS_TABLE_COLUMNS, rows)
    write_summary(args.output_json, metadata)
    write_cross_artifact(
        args.output_h5,
        BACKBONE_CLASSES,
        atomtypes,
        radial_grid,
        table_tensor,
        metadata,
    )
    print(f"Wrote spline cross table CSV: {args.output_csv}")
    print(f"Wrote spline cross table metadata: {args.output_json}")
    print(f"Wrote spline cross table artifact: {args.output_h5}")


def run_inject_backbone_cross_command(argv):
    parser = argparse.ArgumentParser(description="Inject a MARTINI backbone cross-potential node into prepared .up files.")
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--up-files", nargs="*", default=[])
    parser.add_argument("--overwrite-existing", type=int, choices=[0, 1], default=1)
    args = parser.parse_args(argv)
    artifact_path = require_existing_file(args.artifact)
    up_paths = collect_cross_up_paths(args.manifest, args.up_files)
    artifact = load_cross_artifact(artifact_path)
    if artifact["protein_classes"] != list(BACKBONE_CLASSES):
        raise ValueError(
            f"{artifact_path}: expected protein classes {list(BACKBONE_CLASSES)}, "
            f"got {artifact['protein_classes']}"
        )
    for up_path in up_paths:
        protein_indices, protein_class_indices, env_indices, env_class_indices = collect_cross_indices(up_path, artifact)
        inject_cross_node(
            up_path,
            artifact,
            protein_indices,
            protein_class_indices,
            env_indices,
            env_class_indices,
            bool(args.overwrite_existing),
        )
        node_name = RBM_CROSS_NODE_NAME if artifact["schema"] == RBM_CROSS_SCHEMA else SPLINE_CROSS_NODE_NAME
        print(
            f"Injected {node_name} into {up_path} "
            f"(protein={protein_indices.size}, env={env_indices.size})"
        )


COMMANDS = {
    "handoff": run_handoff_command,
    "inject-backbone-only": run_inject_backbone_only_command,
    "validate-backbone-only": run_validate_backbone_only_command,
    "validate-stage-potentials": run_validate_stage_potentials_command,
    "build-depth-table": run_build_depth_table_command,
    "build-backbone-cross-table": run_build_backbone_cross_table_command,
    "inject-backbone-cross": run_inject_backbone_cross_command,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in COMMANDS:
        COMMANDS[argv[0]](argv[1:])
        return
    main_prepare(argv)


if __name__ == "__main__":
    main()
