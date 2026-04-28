#!/usr/bin/env python3

import argparse
import contextlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from martini_prepare_system_lib import (
    center_of_mass,
    build_sc_martini_h5,
    collect_bb_map,
    convert_stage,
    compute_lipid_residue_indices,
    coords,
    estimate_salt_pairs,
    extract_protein_cg_atoms,
    infer_effective_ion_volume_fraction_from_template,
    infer_protein_charge_from_cg,
    inject_stage7_sc_table_nodes,
    lipid_resname,
    parse_pdb,
    place_ions,
    remove_overlapping_lipids,
    set_initial_position,
    set_box_from_lipid_xy,
    set_coords,
    tile_and_crop_bilayer_lipids,
    validate_hybrid_mapping,
    validate_backbone_reference_frame,
    write_hybrid_mapping_h5,
    write_pdb,
    DEFAULT_SC_TABLE_JSON,
)


PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"


def parse_prepare_args(argv=None):
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

    parser.add_argument("--bilayer-pdb", default=str(REPO_ROOT / "parameters" / "dryMARTINI" / "DOPC.pdb"))
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
        "--protein-placement-mode",
        choices=["embed", "outside-top", "outside-bottom"],
        default="embed",
        help="How to place the protein relative to the bilayer during mixed-system preparation.",
    )
    parser.add_argument(
        "--protein-orientation-mode",
        choices=["input", "lay-flat"],
        default="input",
        help="Whether to keep the input protein orientation or rotate it so its thinnest axis becomes the bilayer normal.",
    )
    parser.add_argument(
        "--protein-surface-gap",
        type=float,
        default=6.0,
        help="Requested minimum gap in Angstrom between the protein surface and the lipid surface for outside-of-bilayer placement.",
    )
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
        else (WORKFLOW_DIR / "pdb" / f"{args.pdb_id}.MARTINI.pdb")
    )
    runtime_itp = (
        Path(args.runtime_itp_output).expanduser().resolve()
        if args.runtime_itp_output
        else (WORKFLOW_DIR / "pdb" / f"{args.pdb_id}_proa.itp")
    )
    return runtime_pdb, runtime_itp


def copy_if_different(src: Path, dst: Path):
    src_resolved = src.expanduser().resolve()
    dst_resolved = dst.expanduser().resolve()
    if src_resolved == dst_resolved:
        return
    shutil.copy2(src_resolved, dst_resolved)


def canonicalize_axis_sign(axis):
    axis = np.asarray(axis, dtype=float)
    dominant = int(np.argmax(np.abs(axis)))
    if axis[dominant] < 0.0:
        axis = -axis
    return axis


def lay_flat_rotation_basis(protein_xyz):
    centered = protein_xyz - center_of_mass(protein_xyz)
    if centered.shape[0] < 3:
        return np.eye(3, dtype=float)

    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]

    axis_x = canonicalize_axis_sign(eigvecs[:, order[0]])
    axis_z = canonicalize_axis_sign(eigvecs[:, order[-1]])
    axis_y = np.cross(axis_z, axis_x)
    axis_y_norm = float(np.linalg.norm(axis_y))
    if axis_y_norm <= 1.0e-8:
        axis_y = canonicalize_axis_sign(eigvecs[:, order[1]])
    else:
        axis_y = axis_y / axis_y_norm

    axis_z = np.cross(axis_x, axis_y)
    axis_z = axis_z / np.linalg.norm(axis_z)
    basis = np.column_stack((axis_x, axis_y, axis_z))
    if np.linalg.det(basis) < 0.0:
        basis[:, 1] *= -1.0
    return basis


def orient_protein_xyz(protein_xyz, orientation_mode):
    protein_center = center_of_mass(protein_xyz)
    centered = protein_xyz - protein_center
    if orientation_mode == "input":
        return centered
    if orientation_mode == "lay-flat":
        return centered @ lay_flat_rotation_basis(protein_xyz)
    raise ValueError(f"Unsupported protein orientation mode: {orientation_mode}")


def place_protein_xyz(protein_xyz_centered, bilayer_xyz, bilayer_center, placement_mode, surface_gap):
    if placement_mode == "embed":
        return protein_xyz_centered + bilayer_center

    translation = np.array([float(bilayer_center[0]), float(bilayer_center[1]), 0.0], dtype=float)
    if placement_mode == "outside-top":
        translation[2] = float(bilayer_xyz[:, 2].max()) + float(surface_gap) - float(
            protein_xyz_centered[:, 2].min()
        )
        return protein_xyz_centered + translation
    if placement_mode == "outside-bottom":
        translation[2] = float(bilayer_xyz[:, 2].min()) - float(surface_gap) - float(
            protein_xyz_centered[:, 2].max()
        )
        return protein_xyz_centered + translation

    raise ValueError(f"Unsupported protein placement mode: {placement_mode}")


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
    ff_dir = Path(
        os.environ.get("UPSIDE_MARTINI_FF_DIR", str(REPO_ROOT / "parameters" / "dryMARTINI"))
    ).expanduser()
    if not ff_dir.is_absolute():
        ff_dir = (REPO_ROOT / ff_dir).resolve()
    else:
        ff_dir = ff_dir.resolve()
    ff_file = ff_dir / "dry_martini_v2.1.itp"
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
    if not args.protein_cg_pdb:
        raise ValueError("--protein-cg-pdb is required for mode=both")
    if not args.bilayer_pdb:
        raise ValueError("--bilayer-pdb is required for mode=both")

    protein_cg_pdb = Path(args.protein_cg_pdb).expanduser().resolve()
    bilayer_pdb = Path(args.bilayer_pdb).expanduser().resolve()
    if not protein_cg_pdb.exists():
        raise FileNotFoundError(f"Protein CG PDB not found: {protein_cg_pdb}")
    if not bilayer_pdb.exists():
        raise FileNotFoundError(f"Bilayer PDB not found: {bilayer_pdb}")

    protein_atoms_raw, _ = parse_pdb(protein_cg_pdb)
    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    protein_atoms = extract_protein_cg_atoms(protein_atoms_raw)
    bilayer_lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not protein_atoms:
        raise ValueError("No atoms found in protein CG PDB.")
    if not bilayer_lipid_atoms:
        raise ValueError("No lipid residues found in bilayer template.")

    protein_xyz_raw = coords(protein_atoms)
    bilayer_xyz = coords(bilayer_lipid_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    bmin = bilayer_xyz.min(axis=0)
    bmax = bilayer_xyz.max(axis=0)
    base_side = max(float(bmax[0] - bmin[0]), float(bmax[1] - bmin[1]))

    if args.protein_placement_mode == "embed" and args.protein_orientation_mode == "input":
        protein_center = center_of_mass(protein_xyz_raw)
        shift = bilayer_center - protein_center
        protein_xyz = protein_xyz_raw + shift
        set_coords(protein_atoms, protein_xyz)

        pmin = protein_xyz.min(axis=0)
        pmax = protein_xyz.max(axis=0)
        pspan = pmax - pmin
        pcenter_xy = 0.5 * (pmin[:2] + pmax[:2])

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
    else:
        protein_xyz_oriented = orient_protein_xyz(protein_xyz_raw, args.protein_orientation_mode)
        pmin = protein_xyz_oriented.min(axis=0)
        pmax = protein_xyz_oriented.max(axis=0)
        pspan = pmax - pmin

        min_required_xy = 3.0 * pspan[:2] + 2.0 * float(args.box_padding_xy)
        target_side = max(base_side * float(max(1.0, args.xy_scale)), float(np.max(min_required_xy)))
        target_xy_min = np.array(
            [bilayer_center[0] - 0.5 * target_side, bilayer_center[1] - 0.5 * target_side],
            dtype=float,
        )
        target_xy_max = np.array(
            [bilayer_center[0] + 0.5 * target_side, bilayer_center[1] + 0.5 * target_side],
            dtype=float,
        )
        bilayer_lipids = tile_and_crop_bilayer_lipids(
            bilayer_atoms=bilayer_atoms,
            bilayer_box=bilayer_box,
            target_xy_min=target_xy_min,
            target_xy_max=target_xy_max,
        )

        bilayer_lipid_atoms_cropped = [a for a in bilayer_lipids if lipid_resname(a["resname"])]
        if not bilayer_lipid_atoms_cropped:
            raise ValueError("No lipid residues remained after bilayer tiling/cropping.")
        bilayer_xyz_cropped = coords(bilayer_lipid_atoms_cropped)
        protein_xyz = place_protein_xyz(
            protein_xyz_centered=protein_xyz_oriented,
            bilayer_xyz=bilayer_xyz_cropped,
            bilayer_center=center_of_mass(bilayer_xyz_cropped),
            placement_mode=args.protein_placement_mode,
            surface_gap=args.protein_surface_gap,
        )
        set_coords(protein_atoms, protein_xyz)
    lipid_residues, keep_nonlipid = compute_lipid_residue_indices(bilayer_lipids)
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_lipids,
        protein_atoms=protein_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=float(args.protein_lipid_cutoff),
    )

    packed_atoms = protein_atoms + bilayer_kept
    min_box_z_target = float(3.0 * pspan[2])
    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        pad_z=float(args.box_padding_z),
        force_square_xy=True,
        min_box_z=min_box_z_target,
        center_lipid_in_z=True,
    )

    effective_vol_frac = infer_effective_ion_volume_fraction_from_template(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        salt_molar=float(args.salt_molar),
    )
    protein_charge = (
        int(args.protein_net_charge)
        if args.protein_net_charge is not None
        else int(infer_protein_charge_from_cg(protein_atoms))
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
    protein_xyz_final = coords(protein_atoms)
    bilayer_kept_lipid_atoms = [a for a in bilayer_kept if lipid_resname(a["resname"])]
    bilayer_xyz_final = coords(bilayer_kept_lipid_atoms) if bilayer_kept_lipid_atoms else coords(bilayer_kept)

    if args.protein_itp:
        protein_itp = Path(args.protein_itp).expanduser().resolve()
        if not protein_itp.exists():
            raise FileNotFoundError(f"Protein ITP not found: {protein_itp}")
        copy_if_different(protein_itp, runtime_itp)

    mapping_summary = {}
    if args.hybrid_mapping_output:
        if not args.protein_aa_pdb:
            raise ValueError(
                "--protein-aa-pdb is required when --hybrid-mapping-output is provided."
            )
        protein_aa_pdb = Path(args.protein_aa_pdb).expanduser().resolve()
        if not protein_aa_pdb.exists():
            raise FileNotFoundError(f"Protein AA PDB not found: {protein_aa_pdb}")

        protein_aa_atoms, _ = parse_pdb(protein_aa_pdb)
        if not protein_aa_atoms:
            raise ValueError(f"No atoms found in protein AA PDB: {protein_aa_pdb}")

        frame_diag = validate_backbone_reference_frame(
            protein_aa_atoms,
            protein_atoms,
            min_matched_residues=max(1, int(args.bb_aa_min_matched_residues)),
            max_rigid_rmsd=float(args.bb_aa_max_rigid_rmsd),
            context=f"{protein_aa_pdb.name} -> {protein_cg_pdb.name}",
        )
        bb_entries = collect_bb_map(protein_aa_atoms, protein_atoms)

        mapping_h5 = Path(args.hybrid_mapping_output).expanduser().resolve()
        mapping_h5.parent.mkdir(parents=True, exist_ok=True)
        env_atom_indices = list(range(len(protein_atoms), len(all_atoms)))
        write_hybrid_mapping_h5(
            mapping_h5,
            bb_entries=bb_entries,
            total_martini_atoms=len(all_atoms),
            env_atom_indices=env_atom_indices,
            n_protein_atoms=len(protein_atoms),
        )

        mapping_summary["protein_aa_pdb"] = str(protein_aa_pdb)
        mapping_summary["mapping_h5"] = str(mapping_h5)
        mapping_summary["backbone_frame_check"] = frame_diag
        mapping_summary["bb_map_entries"] = int(len(bb_entries))

        if args.hybrid_bb_map_json_output:
            mapping_json = Path(args.hybrid_bb_map_json_output).expanduser().resolve()
            write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})
            mapping_summary["mapping_json"] = str(mapping_json)

    summary = {
        "mode": "both",
        "input_protein_cg_pdb": str(protein_cg_pdb),
        "input_bilayer_pdb": str(bilayer_pdb),
        "runtime_pdb": str(runtime_pdb),
        "runtime_itp": str(runtime_itp) if args.protein_itp else None,
        "protein_placement_mode": str(args.protein_placement_mode),
        "protein_orientation_mode": str(args.protein_orientation_mode),
        "protein_surface_gap_requested_angstrom": float(args.protein_surface_gap),
        "xy_scale": float(args.xy_scale),
        "base_xy_side_angstrom": float(base_side),
        "target_xy_side_angstrom": float(target_side),
        "box_angstrom": [float(v) for v in box_lengths],
        "protein_charge_used": int(protein_charge),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "protein_atoms": int(len(protein_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "lipid_residues_removed": int(removed_lipids),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "protein_span_angstrom": [float(v) for v in pspan],
        "protein_final_z_min": float(protein_xyz_final[:, 2].min()),
        "protein_final_z_max": float(protein_xyz_final[:, 2].max()),
        "lipid_final_z_min": float(bilayer_xyz_final[:, 2].min()),
        "lipid_final_z_max": float(bilayer_xyz_final[:, 2].max()),
        "protein_top_clearance_angstrom": float(protein_xyz_final[:, 2].min() - bilayer_xyz_final[:, 2].max()),
        "protein_bottom_clearance_angstrom": float(bilayer_xyz_final[:, 2].min() - protein_xyz_final[:, 2].max()),
    }
    summary.update(mapping_summary)
    return summary


def run_stage_conversion(args, runtime_pdb: Path, runtime_itp: Path):
    prev_pdb = os.environ.get("UPSIDE_RUNTIME_PDB_FILE")
    prev_itp = os.environ.get("UPSIDE_RUNTIME_ITP_FILE")

    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(runtime_pdb)
    if args.mode in {"protein", "both"}:
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


def main():
    args = parse_prepare_args()
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
        if args.mode in {"protein", "both"}:
            if not runtime_itp.exists():
                raise FileNotFoundError(
                    f"Protein ITP required for mode={args.mode} stage conversion: {runtime_itp}. "
                    "Provide --protein-itp when running with --prepare-structure 1."
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


def run_prepare_command(argv):
    args = parse_prepare_args(argv)
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
        if args.mode in {"protein", "both"}:
            if not runtime_itp.exists():
                raise FileNotFoundError(
                    f"Protein ITP required for mode={args.mode} stage conversion: {runtime_itp}. "
                    "Provide --protein-itp when running with --prepare-structure 1."
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


def run_build_sc_martini_h5_command(argv):
    parser = argparse.ArgumentParser(
        description="Build a native-unit martini.h5 SC table from assembled SC-training JSON."
    )
    parser.add_argument(
        "--sc-table-json",
        default=str(DEFAULT_SC_TABLE_JSON),
        help="Path to assembled sc_table.json from SC-training.",
    )
    parser.add_argument(
        "--output-h5",
        default="parameters/ff_2.1/martini.h5",
        help="Output HDF5 file path.",
    )
    args = parser.parse_args(argv)
    build_sc_martini_h5(Path(args.sc_table_json), Path(args.output_h5))


def run_validate_hybrid_mapping_command(argv):
    parser = argparse.ArgumentParser(
        description="Validate hybrid mapping HDF5 schema and index consistency."
    )
    parser.add_argument("mapping_h5", type=Path, help="Path to hybrid mapping HDF5 file.")
    parser.add_argument("--n-atom", type=int, default=None, help="Optional expected n_atom value.")
    args = parser.parse_args(argv)
    validate_hybrid_mapping(args.mapping_h5, n_atom=args.n_atom)


def run_set_initial_position_command(argv):
    parser = argparse.ArgumentParser(
        description="Copy final coordinates from one stage file into the next stage input."
    )
    parser.add_argument("input_file", help="Source stage file.")
    parser.add_argument("output_file", help="Target stage file to update.")
    args = parser.parse_args(argv)
    set_initial_position(args.input_file, args.output_file)


def run_inject_stage7_sc_command(argv):
    parser = argparse.ArgumentParser(
        description="Inject stage-7 dry-MARTINI SC table coupling into a prepared .up file."
    )
    parser.add_argument("up_file", help="Target stage-7 .up file to modify in place.")
    parser.add_argument("martini_h5", help="Native-unit martini.h5 table library.")
    parser.add_argument("upside_home", help="UPSIDE_HOME used to locate upside_config.py.")
    parser.add_argument("rama_library", help="Upside rama.dat path.")
    parser.add_argument("rama_sheet_mixing", help="Upside sheet-mixing path.")
    parser.add_argument("hbond_energy", help="Upside hbond.h5 path.")
    parser.add_argument("reference_state_rama", help="Upside rama_reference.pkl path.")
    parser.add_argument(
        "--protein-itp",
        help="Protein ITP used to recover residue names when /input/sequence is absent or mismatched.",
    )
    args = parser.parse_args(argv)
    inject_stage7_sc_table_nodes(
        up_file=Path(args.up_file),
        martini_h5=Path(args.martini_h5),
        upside_home=Path(args.upside_home),
        rama_library=Path(args.rama_library),
        rama_sheet_mixing=Path(args.rama_sheet_mixing),
        hbond_energy=Path(args.hbond_energy),
        reference_state_rama=Path(args.reference_state_rama),
        protein_itp=Path(args.protein_itp) if args.protein_itp else None,
    )


PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX",
}
LIPID_RESIDUES = {"DOP", "DOPC"}
SC_LIBRARY_REQUIRED_DATASETS = [
    "grid_nm",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_kj_mol",
    "rotamer_angular_energy_kj_mol",
    "rotamer_angular_profile",
]


@contextlib.contextmanager
def temporary_env(updates):
    previous = {key: os.environ.get(key) for key in updates}
    try:
        for key, value in updates.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = str(value)
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def env_default(name, default):
    return os.environ.get(name, default)


def env_int(name, default):
    return int(os.environ.get(name, str(default)))


def env_float(name, default):
    return float(os.environ.get(name, str(default)))


def workflow_path(value, base=WORKFLOW_DIR):
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def generate_random_seed():
    seed = random.SystemRandom().randint(1, 2**32 - 1)
    return seed if seed != 0 else 1


def split_comment(line: str):
    if ";" in line:
        body, comment = line.split(";", 1)
        return body.rstrip(), ";" + comment
    return line.rstrip(), ""


def run_checked(cmd, cwd=None, log_file=None):
    print(" ".join(str(x) for x in cmd))
    if log_file is None:
        subprocess.run([str(x) for x in cmd], cwd=str(cwd) if cwd else None, check=True)
        return

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, [str(x) for x in cmd])


def parse_top_molecule_counts(top_file: Path, moleculetype_name: str):
    copy_count = 0
    in_molecules = False
    for raw in top_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            in_molecules = stripped[1:-1].strip().lower() == "molecules"
            continue
        if not in_molecules:
            continue
        parts = stripped.split()
        if len(parts) >= 2 and parts[0] == moleculetype_name:
            copy_count += int(parts[1])
    return max(copy_count, 1)


def resolve_itp_from_top(top_file: Path):
    top_file = top_file.resolve()
    candidates = []
    include_re = re.compile(r'#include\s+"([^"]+\.itp)"')
    for line in top_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = include_re.search(line)
        if not match:
            continue
        path = top_file.parent / match.group(1)
        if path.exists():
            candidates.append(path)

    def has_atoms(path):
        text = path.read_text(encoding="utf-8", errors="ignore").lower()
        return "[ atoms ]" in text or "[atoms]" in text

    for candidate in candidates:
        if has_atoms(candidate):
            return candidate
    for candidate in sorted(top_file.parent.glob("*.itp")):
        if has_atoms(candidate):
            return candidate
    raise FileNotFoundError(f"Could not resolve protein ITP from {top_file}")


def materialize_runtime_itp_from_top(top_file: Path, source_itp: Path, out_itp: Path):
    source_lines = source_itp.read_text(encoding="utf-8", errors="ignore").splitlines()
    moleculetype_name = None
    current_section = None
    atom_count = 0
    residue_count = 0
    for raw in source_lines:
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
            continue
        if not stripped or stripped.startswith(";"):
            continue
        if current_section == "moleculetype" and moleculetype_name is None:
            moleculetype_name = stripped.split()[0]
        elif current_section == "atoms":
            parts = stripped.split()
            if len(parts) >= 6:
                atom_count = max(atom_count, int(parts[0]))
                residue_count = max(residue_count, int(parts[2]))

    if not moleculetype_name or atom_count <= 0 or residue_count <= 0:
        raise ValueError(f"Could not resolve molecule shape from {source_itp}")

    copy_count = parse_top_molecule_counts(top_file, moleculetype_name)
    if copy_count == 1:
        out_itp.write_text("\n".join(source_lines) + "\n", encoding="utf-8")
        return

    section_index_offsets = {
        "bonds": 2,
        "constraints": 2,
        "pairs": 2,
        "angles": 3,
        "dihedrals": 4,
        "position_restraints": 1,
    }

    def transform_line(raw: str, section: str | None, atom_offset: int, residue_offset: int):
        body, comment = split_comment(raw)
        stripped = body.strip()
        if not stripped or stripped.startswith(";") or (stripped.startswith("[") and stripped.endswith("]")):
            return raw
        parts = body.split()
        if section == "atoms" and len(parts) >= 6:
            parts[0] = str(int(parts[0]) + atom_offset)
            parts[2] = str(int(parts[2]) + residue_offset)
            parts[5] = str(int(parts[5]) + atom_offset)
        elif section in section_index_offsets:
            for idx in range(min(section_index_offsets[section], len(parts))):
                parts[idx] = str(int(parts[idx]) + atom_offset)
        elif section == "exclusions":
            for idx in range(len(parts)):
                parts[idx] = str(int(parts[idx]) + atom_offset)
        transformed = " ".join(parts)
        return f"{transformed} {comment}" if comment else transformed

    rendered = []
    for copy_index in range(copy_count):
        atom_offset = copy_index * atom_count
        residue_offset = copy_index * residue_count
        current_section = None
        if copy_index > 0:
            rendered.append("")
        for raw in source_lines:
            stripped = raw.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                current_section = stripped[1:-1].strip().lower()
                rendered.append(raw)
                continue
            rendered.append(transform_line(raw, current_section, atom_offset, residue_offset))
    out_itp.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def prepare_workflow_protein_inputs(args):
    mass_ff_file = workflow_path(args.mass_ff_file).resolve()
    if args.martinize_enable:
        martinize_dir = workflow_path(args.run_dir) / "martinize"
        martinize_dir.mkdir(parents=True, exist_ok=True)
        cg_pdb = martinize_dir / f"{args.runtime_pdb_id}.MARTINI.pdb"
        top_file = martinize_dir / f"{args.runtime_pdb_id}.top"
        cmd = [
            sys.executable,
            str(workflow_path(args.martinize_script).resolve()),
            "-f",
            str(workflow_path(args.protein_aa_pdb).resolve()),
            "-x",
            str(cg_pdb.resolve()),
            "-o",
            str(top_file.resolve()),
            "-ff",
            args.martinize_ff,
            "-name",
            args.martinize_molname,
        ]
        print("Running MARTINI martinize command:")
        run_checked(cmd, cwd=martinize_dir)
        source_itp = resolve_itp_from_top(top_file)
        runtime_itp = martinize_dir / f"{args.runtime_pdb_id}.runtime.itp"
        materialize_runtime_itp_from_top(top_file, source_itp, runtime_itp)
        assert_protein_itp_mass_compatibility(runtime_itp)
        return cg_pdb, runtime_itp

    if not args.protein_cg_pdb or not args.protein_itp:
        raise ValueError("MARTINIZE_ENABLE=0 requires protein CG PDB and protein ITP.")
    protein_cg = workflow_path(args.protein_cg_pdb).resolve()
    protein_itp = workflow_path(args.protein_itp).resolve()
    if not protein_cg.exists():
        raise FileNotFoundError(protein_cg)
    if not protein_itp.exists():
        raise FileNotFoundError(protein_itp)
    if not mass_ff_file.exists():
        raise FileNotFoundError(mass_ff_file)
    assert_protein_itp_mass_compatibility(protein_itp)
    return protein_cg, protein_itp


def pdb_min_protein_lipid_distance(pdb_file: Path):
    protein_xyz = []
    lipid_xyz = []
    with pdb_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:21].strip().upper()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if resname in PROTEIN_RESIDUES:
                protein_xyz.append(xyz)
            elif resname in LIPID_RESIDUES:
                lipid_xyz.append(xyz)
    if not protein_xyz or not lipid_xyz:
        return float("nan")
    protein = np.asarray(protein_xyz, dtype=float)
    lipid = np.asarray(lipid_xyz, dtype=float)
    min_d2 = float("inf")
    for i in range(0, lipid.shape[0], 2000):
        block = lipid[i:i + 2000]
        d = block[:, None, :] - protein[None, :, :]
        min_d2 = min(min_d2, float(np.einsum("ijk,ijk->ij", d, d).min()))
    return min_d2 ** 0.5


def prepare_workflow_hybrid_artifacts(args, protein_cg: Path, protein_itp: Path):
    print("=== Stage 0: Hybrid Packing + Mapping Export ===")
    cutoff = float(args.protein_lipid_cutoff)
    while True:
        print(f"Hybrid packing attempt with protein-lipid cutoff: {cutoff:.3f} A")
        run_prepare_command(
            [
                "--mode", "both",
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.hybrid_packed_pdb),
                "--runtime-itp-output", str(args.runtime_itp_file),
                "--prepare-structure", "1",
                "--protein-aa-pdb", str(workflow_path(args.protein_aa_pdb).resolve()),
                "--protein-cg-pdb", str(protein_cg.resolve()),
                "--protein-itp", str(protein_itp.resolve()),
                "--hybrid-mapping-output", str(args.hybrid_mapping_file),
                "--hybrid-bb-map-json-output", str(args.hybrid_prep_dir / "hybrid_bb_map.json"),
                "--bilayer-pdb", str(workflow_path(args.bilayer_pdb).resolve()),
                "--salt-molar", str(args.salt_molar),
                "--protein-lipid-cutoff", f"{cutoff:.6g}",
                "--ion-cutoff", str(args.ion_cutoff),
                "--xy-scale", str(args.xy_scale),
                "--box-padding-xy", str(args.box_padding_xy),
                "--box-padding-z", str(args.box_padding_z),
                "--protein-placement-mode", args.protein_placement_mode,
                "--protein-orientation-mode", args.protein_orientation_mode,
                "--protein-surface-gap", str(args.protein_surface_gap),
                "--bb-aa-min-matched-residues", str(args.bb_aa_min_matched_residues),
                "--bb-aa-max-rigid-rmsd", str(args.bb_aa_max_rigid_rmsd),
                "--seed", str(args.prep_seed),
                "--summary-json", str(args.hybrid_prep_dir / "hybrid_prep_summary.json"),
            ]
        )
        if not args.hybrid_packed_pdb.exists():
            raise FileNotFoundError(args.hybrid_packed_pdb)
        min_gap = pdb_min_protein_lipid_distance(args.hybrid_packed_pdb)
        if np.isfinite(min_gap) and min_gap >= float(args.protein_lipid_min_gap):
            print(
                f"Hybrid packing accepted: min protein-lipid distance {min_gap:.6f} A "
                f"(target >= {args.protein_lipid_min_gap} A)"
            )
            break
        if cutoff >= float(args.protein_lipid_cutoff_max):
            raise RuntimeError(
                "Hybrid packing still overpacked near protein.\n"
                f"  Observed min protein-lipid distance: {min_gap:.6f} A\n"
                f"  Target min distance: {args.protein_lipid_min_gap} A\n"
                f"  Reached cutoff limit: {args.protein_lipid_cutoff_max} A"
            )
        cutoff += float(args.protein_lipid_cutoff_step)
        print(f"Hybrid packing too tight near protein. Retrying with cutoff {cutoff:.3f} A.")

    if not args.hybrid_mapping_file.exists():
        raise FileNotFoundError(args.hybrid_mapping_file)
    if args.hybrid_validate:
        validate_hybrid_mapping(args.hybrid_mapping_file)
    copy_if_different(args.hybrid_packed_pdb, args.runtime_pdb_file)
    copy_if_different(protein_itp, args.runtime_itp_file)
    print(f"Runtime MARTINI PDB: {args.runtime_pdb_file}")
    print(f"Runtime protein ITP: {args.runtime_itp_file}")


def replace_dataset(parent, name, data):
    if name in parent:
        del parent[name]
    parent.create_dataset(name, data=data)


def h5_as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def pad_bytes(value, dtype):
    n = np.dtype(dtype).itemsize
    raw = str(value).encode("ascii", errors="ignore")[:n]
    return raw.ljust(n, b" ")


def inject_hybrid_mapping(up_file: Path, mapping_file: Path):
    import h5py

    groups = ["hybrid_control", "hybrid_bb_map", "hybrid_env_topology"]
    optional_groups = ("chain_break",)
    component_mass = {"N": 14.0, "CA": 12.0, "C": 12.0, "O": 16.0}
    component_names_default = ("N", "CA", "C", "O")
    component_scale = 72.0 / float(sum(component_mass.values()))

    with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
        src_inp = src["/input"]
        dst_inp = dst.require_group("input")
        base_n_atom = int(dst_inp["pos"].shape[0])
        src_mem = src_inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32)
        if base_n_atom != int(src_mem.shape[0]):
            raise ValueError(
                f"Hybrid mapping n_atom mismatch for {up_file}: up has {base_n_atom}, "
                f"mapping has {src_mem.shape[0]}"
            )
        for group in groups:
            if group not in src_inp:
                raise ValueError(f"Missing mapping group in {mapping_file}: /input/{group}")
            if group in dst_inp:
                del dst_inp[group]
            src.copy(src_inp[group], dst_inp, name=group)
        for group in optional_groups:
            if group not in src_inp:
                continue
            if group in dst_inp:
                del dst_inp[group]
            src.copy(src_inp[group], dst_inp, name=group)

        bb_grp = dst_inp["hybrid_bb_map"]
        env_grp = dst_inp["hybrid_env_topology"]
        pot_grp = dst_inp["potential"]["martini_potential"]

        bb_residue = bb_grp["bb_residue_index"][:].astype(np.int32)
        bb_ref_idx = (
            bb_grp["reference_atom_indices"][:].astype(np.int32)
            if "reference_atom_indices" in bb_grp
            else bb_grp["atom_indices"][:].astype(np.int32)
        )
        n_bb = int(bb_ref_idx.shape[0])
        if bb_ref_idx.ndim != 2 or bb_ref_idx.shape[1] != 4:
            raise ValueError("hybrid_bb_map reference/active indices must have shape (n_bb,4)")
        if "reference_atom_coords" in bb_grp:
            bb_ref_xyz = bb_grp["reference_atom_coords"][:].astype(np.float32)
        else:
            bb_ref_xyz = np.zeros((n_bb, 4, 3), dtype=np.float32)
        if bb_ref_xyz.shape != (n_bb, 4, 3):
            raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")

        comp_names = list(component_names_default)
        if "reference_atom_names" in bb_grp and bb_grp["reference_atom_names"][:].shape == (4,):
            comp_names = [
                h5_as_text(x).strip().upper() or component_names_default[i]
                for i, x in enumerate(bb_grp["reference_atom_names"][:])
            ]

        ref_values = []
        ref_lookup = {}
        for raw_ref_idx in bb_ref_idx.reshape(-1):
            raw_ref_idx = int(raw_ref_idx)
            if raw_ref_idx < 0 or raw_ref_idx in ref_lookup:
                continue
            ref_lookup[raw_ref_idx] = len(ref_values)
            ref_values.append(raw_ref_idx)
        n_ref_index = len(ref_values)
        ref_offset = base_n_atom

        if n_ref_index > 0:
            pos_dtype = dst_inp["pos"].dtype
            mom_dtype = dst_inp["mom"].dtype
            vel_dtype = dst_inp["vel"].dtype
            mass_dtype = dst_inp["mass"].dtype
            charge_dtype = dst_inp["charges"].dtype
            type_dtype = dst_inp["type"].dtype
            name_dtype = dst_inp["atom_names"].dtype if "atom_names" in dst_inp else type_dtype
            role_dtype = dst_inp["atom_roles"].dtype if "atom_roles" in dst_inp else type_dtype
            residue_dtype = dst_inp["residue_ids"].dtype
            molecule_dtype = dst_inp["molecule_ids"].dtype
            pclass_dtype = dst_inp["particle_class"].dtype

            ref_pos = np.zeros((n_ref_index, 3, 1), dtype=pos_dtype)
            ref_mom = np.zeros((n_ref_index, 3, 1), dtype=mom_dtype)
            ref_vel = np.zeros((n_ref_index, 3), dtype=vel_dtype)
            ref_mass = np.ones((n_ref_index,), dtype=mass_dtype)
            ref_charge = np.zeros((n_ref_index,), dtype=charge_dtype)
            ref_type = np.full((n_ref_index,), pad_bytes("AA", type_dtype), dtype=type_dtype)
            ref_name = np.full((n_ref_index,), pad_bytes("AA", name_dtype), dtype=name_dtype)
            ref_role = np.full((n_ref_index,), pad_bytes("AA", role_dtype), dtype=role_dtype)
            ref_residue = np.full((n_ref_index,), -1, dtype=residue_dtype)
            ref_molecule = np.full((n_ref_index,), 0, dtype=molecule_dtype)
            ref_pclass = np.full((n_ref_index,), pad_bytes("PROTEINAA", pclass_dtype), dtype=pclass_dtype)

            for k in range(n_bb):
                resid = int(bb_residue[k])
                for d in range(4):
                    ref_idx = int(bb_ref_idx[k, d])
                    if ref_idx < 0:
                        continue
                    compact_idx = ref_lookup[ref_idx]
                    ref_pos[compact_idx, :, 0] = bb_ref_xyz[k, d, :]
                    cname = comp_names[d] if d < len(comp_names) else component_names_default[d]
                    ref_name[compact_idx] = pad_bytes(cname, name_dtype)
                    ref_role[compact_idx] = pad_bytes(cname, role_dtype)
                    ref_residue[compact_idx] = resid
                    ref_mass[compact_idx] = mass_dtype.type(
                        (component_scale * float(component_mass.get(cname, 12.0))) / 12.0
                    )

            for name, appendix in [
                ("pos", ref_pos),
                ("mom", ref_mom),
                ("vel", ref_vel),
                ("mass", ref_mass),
                ("charges", ref_charge),
                ("type", ref_type),
                ("atom_names", ref_name),
                ("atom_roles", ref_role),
                ("residue_ids", ref_residue),
                ("molecule_ids", ref_molecule),
                ("particle_class", ref_pclass),
            ]:
                if name not in dst_inp:
                    continue
                base = dst_inp[name][:]
                if base.shape[0] != base_n_atom:
                    raise ValueError(f"{up_file}: /input/{name} length changed before augmentation")
                replace_dataset(dst_inp, name, np.concatenate([base, appendix.astype(base.dtype, copy=False)], axis=0))

            pot_atom_indices = pot_grp["atom_indices"][:]
            pot_ai_append = np.arange(ref_offset, ref_offset + n_ref_index, dtype=pot_atom_indices.dtype)
            replace_dataset(pot_grp, "atom_indices", np.concatenate([pot_atom_indices, pot_ai_append], axis=0))
            pot_charges = pot_grp["charges"][:]
            replace_dataset(
                pot_grp,
                "charges",
                np.concatenate([pot_charges, np.zeros((n_ref_index,), dtype=pot_charges.dtype)], axis=0),
            )

        n_atom_aug = int(dst_inp["pos"].shape[0])
        bb_runtime_idx = np.full((n_bb, 4), -1, dtype=np.int32)
        bb_runtime_mask = np.zeros((n_bb, 4), dtype=np.int8)
        bb_runtime_w = np.zeros((n_bb, 4), dtype=np.float32)
        for k in range(n_bb):
            raw_w = np.zeros((4,), dtype=np.float32)
            for d in range(4):
                ref_idx = int(bb_ref_idx[k, d])
                if ref_idx < 0:
                    continue
                bb_runtime_idx[k, d] = ref_offset + ref_lookup[ref_idx]
                bb_runtime_mask[k, d] = 1
                cname = comp_names[d] if d < len(comp_names) else component_names_default[d]
                raw_w[d] = np.float32(component_mass.get(cname, 12.0))
            wsum = float(raw_w.sum())
            if wsum > 0.0:
                bb_runtime_w[k, :] = raw_w / wsum

        replace_dataset(bb_grp, "atom_indices", bb_runtime_idx)
        replace_dataset(bb_grp, "atom_mask", bb_runtime_mask)
        replace_dataset(bb_grp, "weights", bb_runtime_w)
        bb_grp.attrs["atom_index_space"] = np.bytes_("stage_runtime")
        bb_grp.attrs["reference_index_space"] = np.bytes_("protein_aa_pdb_0based")
        bb_grp.attrs["reference_index_offset"] = np.int32(ref_offset)
        bb_grp.attrs["reference_index_count"] = np.int32(n_ref_index)

        membership = np.full((n_atom_aug,), -1, dtype=np.int32)
        membership[:base_n_atom] = src_mem
        if n_ref_index > 0:
            membership[ref_offset:ref_offset + n_ref_index] = 0
        replace_dataset(env_grp, "protein_membership", membership)
        if env_grp["protein_membership"].shape[0] != dst_inp["pos"].shape[0]:
            raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch")


def set_stage_label(up_file: Path, stage_label: str):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("stage_parameters")
        grp.attrs["enable"] = np.int8(1)
        grp.attrs["current_stage"] = np.bytes_(stage_label)


def set_hybrid_control_mode(up_file: Path, activation_stage: str, preprod_mode="rigid"):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("hybrid_control")
        grp.attrs["enable"] = np.int8(1)
        grp.attrs["activation_stage"] = np.bytes_(activation_stage)
        grp.attrs["preprod_protein_mode"] = np.bytes_(preprod_mode)


def set_hybrid_production_controls(up_file: Path, args):
    import h5py

    scale = float(args.protein_env_interface_scale)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(f"protein_env_interface_scale must be finite and > 0, got {scale!r}")
    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("hybrid_control")
        grp.attrs["production_nonprotein_hard_sphere"] = np.int8(0)
        grp.attrs["protein_env_interface_scale"] = np.float32(scale)
        grp.attrs["sc_env_lj_force_cap"] = np.float32(args.sc_env_lj_force_cap)
        grp.attrs["sc_env_coul_force_cap"] = np.float32(args.sc_env_coul_force_cap)
        grp.attrs["sc_env_relax_steps"] = np.int32(args.sc_env_relax_steps)
        grp.attrs["sc_env_backbone_hold_steps"] = np.int32(args.sc_env_backbone_hold_steps)
        grp.attrs["sc_env_po4_z_hold_steps"] = np.int32(args.sc_env_po4_z_hold_steps)
        grp.attrs["sc_env_po4_z_clamp_enabled"] = np.int8(1 if args.sc_env_po4_z_clamp_enable else 0)
        grp.attrs["nonprotein_hs_force_cap"] = np.float32(args.nonprotein_hs_force_cap)


def set_barostat_type(up_file: Path, barostat_type: int):
    import tables as tb

    with tb.open_file(up_file, "r+") as t:
        if "/input/barostat" in t:
            t.root.input.barostat._v_attrs.type = int(barostat_type)


def set_production_backbone_fix_rigid(up_file: Path):
    import h5py

    required_roles = ("BB", "N", "CA", "C", "O")
    with h5py.File(up_file, "r+") as h5:
        inp = h5.require_group("input")
        membership = inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32, copy=False)
        n_atom = int(inp["pos"].shape[0])
        if membership.shape[0] != n_atom:
            raise ValueError(f"{up_file}: protein_membership length mismatch")
        role_source = "atom_roles" if "atom_roles" in inp else "atom_names"
        role_values = [h5_as_text(x).strip().upper() for x in inp[role_source][:]]
        role_counts = {role: 0 for role in required_roles}
        selected = []
        for atom_idx, role in enumerate(role_values):
            if membership[atom_idx] < 0 or role not in role_counts:
                continue
            selected.append(atom_idx)
            role_counts[role] += 1
        missing = [role for role, count in role_counts.items() if count == 0]
        if missing:
            raise ValueError(f"{up_file}: missing protein backbone roles: {','.join(missing)}")
        bb_count = role_counts["BB"]
        mismatched = [role for role in required_roles[1:] if role_counts[role] != bb_count]
        if mismatched:
            raise ValueError(f"{up_file}: inconsistent backbone role counts: {role_counts}")
        existing = np.zeros((0,), dtype=np.int32)
        if "fix_rigid" in inp and "atom_indices" in inp["fix_rigid"]:
            existing = inp["fix_rigid"]["atom_indices"][:].astype(np.int32, copy=False)
        merged = np.unique(np.concatenate([existing, np.asarray(selected, dtype=np.int32)]))
        if "fix_rigid" in inp:
            del inp["fix_rigid"]
        grp = inp.require_group("fix_rigid")
        grp.attrs["enable"] = np.int8(1)
        grp.attrs["selection"] = np.bytes_("protein_backbone_roles_bb_n_ca_c_o")
        grp.attrs["selection_source"] = np.bytes_(role_source)
        grp.attrs["selection_count"] = np.int32(merged.shape[0])
        for role, count in role_counts.items():
            grp.attrs[f"count_{role}"] = np.int32(count)
        grp.create_dataset("atom_indices", data=merged.astype(np.int32, copy=False))


def ensure_sc_martini_library(args):
    import h5py

    library = args.sc_martini_library
    need_build = True
    if library.exists():
        with h5py.File(library, "r") as h5:
            missing = [name for name in SC_LIBRARY_REQUIRED_DATASETS if name not in h5]
        need_build = bool(missing)
        if need_build:
            print(f"NOTICE: {library} is missing required SC datasets; rebuilding it.")
    if need_build:
        if not args.sc_martini_table_json.exists():
            raise FileNotFoundError(args.sc_martini_table_json)
        print(f"Building {library} from {args.sc_martini_table_json}")
        build_sc_martini_h5(args.sc_martini_table_json, library)
    with h5py.File(library, "r") as h5:
        missing = [name for name in SC_LIBRARY_REQUIRED_DATASETS if name not in h5]
    if missing:
        raise ValueError(f"{library} missing required datasets after build: {','.join(missing)}")


def assert_hybrid_stage_active(up_file: Path, expected_stage: str, expected_activation: str):
    import h5py

    with h5py.File(up_file, "r") as h5:
        inp = h5["/input"]
        stage = h5_as_text(inp["stage_parameters"].attrs.get("current_stage", b"")).strip()
        hy = inp["hybrid_control"].attrs
        enable = int(hy.get("enable", 0))
        activation = h5_as_text(hy.get("activation_stage", b"")).strip()
        if stage != expected_stage or enable != 1 or activation != expected_activation:
            raise ValueError(
                f"{up_file}: expected stage={expected_stage}, activation={expected_activation}, enable=1; "
                f"got stage={stage}, activation={activation}, enable={enable}"
            )
        if "martini_sc_table_1body" not in inp["potential"] and expected_stage == "production":
            raise ValueError(f"{up_file}: missing production martini_sc_table_1body node")
        env_membership = inp["hybrid_env_topology"]["protein_membership"][:]
        if not np.any(env_membership < 0):
            raise ValueError(f"{up_file}: no non-protein environment targets found")
    print(f"Hybrid activation verified: stage={expected_stage}, activation_stage={expected_activation}, enable=1")


def stage_npt_targets(stage_label: str, args):
    if stage_label in {"6.0", "6.1"}:
        return 0.0, 0.0
    if stage_label in {"6.2", "6.3", "6.4", "6.5", "6.6"}:
        return float(args.bar_1_to_eup_per_a3), float(args.bar_1_to_eup_per_a3)
    return 0.0, 0.0


def stage_conversion_env(args, stage_label: str, prepare_stage: str, npt_enable: int, lipidhead_fc: float):
    target_pxy, target_pz = stage_npt_targets(stage_label, args)
    return {
        "UPSIDE_HOME": str(args.upside_home),
        "UPSIDE_MARTINI_FF_DIR": str(args.martini_ff_dir),
        "UPSIDE_MARTINI_ENERGY_CONVERSION": str(args.martini_energy_conversion),
        "UPSIDE_MARTINI_LENGTH_CONVERSION": str(args.martini_length_conversion),
        "UPSIDE_SIMULATION_STAGE": prepare_stage,
        "UPSIDE_NPT_ENABLE": str(int(npt_enable)),
        "UPSIDE_NPT_TARGET_PXY": str(target_pxy),
        "UPSIDE_NPT_TARGET_PZ": str(target_pz),
        "UPSIDE_NPT_TAU": str(args.npt_tau),
        "UPSIDE_NPT_COMPRESSIBILITY": str(args.compressibility_3e4_bar_inv_to_a3_per_eup),
        "UPSIDE_NPT_COMPRESSIBILITY_XY": str(args.compressibility_3e4_bar_inv_to_a3_per_eup),
        "UPSIDE_NPT_COMPRESSIBILITY_Z": "0.0",
        "UPSIDE_NPT_INTERVAL": str(args.npt_interval),
        "UPSIDE_NPT_SEMI": "1",
        "UPSIDE_NPT_DEBUG": "1",
        "UPSIDE_EWALD_ENABLE": str(args.ewald_enable),
        "UPSIDE_EWALD_ALPHA": str(args.ewald_alpha),
        "UPSIDE_EWALD_KMAX": str(args.ewald_kmax),
        "UPSIDE_BILAYER_LIPIDHEAD_FC": str(lipidhead_fc),
    }


def prepare_stage_file(args, target_file: Path, prepare_stage: str, npt_enable: int, barostat_type: int, lipidhead_fc: float, stage_label: str, protein_itp: Path):
    with temporary_env(stage_conversion_env(args, stage_label, prepare_stage, npt_enable, lipidhead_fc)):
        run_prepare_command(
            [
                "--mode", args.universal_prep_mode,
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.runtime_pdb_file),
                "--runtime-itp-output", str(args.runtime_itp_file),
                "--prepare-structure", "0",
                "--stage", prepare_stage,
                "--run-dir", str(args.run_dir),
                "--summary-json", str(args.hybrid_prep_dir / f"stage_{prepare_stage}.summary.json"),
            ]
        )
    prepared_tmp = args.run_dir / "test.input.up"
    if not prepared_tmp.exists():
        raise FileNotFoundError(prepared_tmp)
    target_file.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(prepared_tmp), str(target_file))
    inject_hybrid_mapping(target_file, args.hybrid_mapping_file)
    set_hybrid_control_mode(target_file, args.hybrid_preprod_activation_stage)
    set_stage_label(target_file, stage_label)
    if stage_label == "production":
        ensure_sc_martini_library(args)
        set_hybrid_control_mode(target_file, "production")
        set_hybrid_production_controls(target_file, args)
        inject_stage7_sc_table_nodes(
            up_file=target_file,
            martini_h5=args.sc_martini_library,
            upside_home=args.upside_home,
            rama_library=args.upside_rama_library,
            rama_sheet_mixing=args.upside_rama_sheet_mixing,
            hbond_energy=args.upside_hbond_energy,
            reference_state_rama=args.upside_reference_state_rama,
            protein_itp=protein_itp,
        )
        if args.prod_70_backbone_fix_rigid_enable:
            set_production_backbone_fix_rigid(target_file)
        assert_hybrid_stage_active(target_file, "production", "production")
    if npt_enable:
        set_barostat_type(target_file, barostat_type)


def handoff_initial_position(args, input_file: Path, output_file: Path, mode="default"):
    refresh = "1" if mode == "production_hybrid" else "0"
    with temporary_env(
        {
            "UPSIDE_SET_INITIAL_STRICT_COPY": str(args.strict_stage_handoff),
            "UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS": refresh,
            "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION": "0",
        }
    ):
        set_initial_position(input_file, output_file)


def run_minimization_stage(args, stage_label: str, up_file: Path, max_iter: int):
    print(f"=== Stage {stage_label}: Minimization ===")
    cmd = [
        args.upside_executable,
        up_file,
        "--duration", "0",
        "--frame-interval", "1",
        "--temperature", args.temperature,
        "--time-step", args.min_time_step,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", args.seed,
        "--integrator", "v",
        "--disable-recentering",
        "--minimize",
        "--min-max-iter", max_iter,
        "--min-energy-tol", "1e-6",
        "--min-force-tol", "1e-3",
        "--min-step", "0.01",
    ]
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.log")


def run_md_stage(args, stage_label: str, input_file: Path, output_file: Path, nsteps: int, dt: float, frame_steps: int):
    effective_frame_steps = int(frame_steps)
    if effective_frame_steps >= int(nsteps):
        effective_frame_steps = max(1, int(nsteps) // 10)
        print(f"NOTICE: frame_steps ({frame_steps}) >= nsteps ({nsteps}); using frame_steps={effective_frame_steps}")
    frame_interval = f"{effective_frame_steps * float(dt):.10g}"
    if input_file.resolve() != output_file.resolve():
        shutil.copy2(input_file, output_file)
        handoff_initial_position(args, input_file, output_file)
    print(f"=== Stage {stage_label}: MD ===")
    cmd = [
        args.upside_executable,
        output_file,
        "--duration-steps", nsteps,
        "--frame-interval", frame_interval,
        "--temperature", args.temperature,
        "--time-step", dt,
        "--thermostat-timescale", args.thermostat_timescale,
        "--thermostat-interval", args.thermostat_interval,
        "--seed", args.seed,
        "--integrator", "v",
        "--disable-recentering",
    ]
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.log")


def extract_stage_vtf(args, stage_label: str, stage_file: Path, mode: str):
    vtf_file = args.run_dir / f"{args.pdb_id}.stage_{stage_label}.vtf"
    print(f"=== Stage {stage_label}: VTF Extraction (mode {mode}) ===")
    cmd = [
        sys.executable,
        args.extract_vtf_script,
        stage_file,
        vtf_file,
        stage_file,
        args.runtime_pdb_id,
        "--mode",
        mode,
        "--split-segments",
    ]
    run_checked(cmd)


def infer_stage70_label_from_file(stage_file: Path, pdb_id: str):
    match = re.fullmatch(rf"{re.escape(pdb_id)}\.stage_(7\.\d+)\.up", stage_file.name)
    return match.group(1) if match else ""


def infer_next_stage70_label(source_file: Path, checkpoint_dir: Path, pdb_id: str):
    pattern = re.compile(rf"^{re.escape(pdb_id)}\.stage_7\.(\d+)\.up$")
    indices = []
    match = pattern.fullmatch(source_file.name)
    if match:
        indices.append(int(match.group(1)))
    if checkpoint_dir.is_dir():
        for path in checkpoint_dir.glob(f"{pdb_id}.stage_7*.up"):
            match = pattern.fullmatch(path.name)
            if path.is_file() and match:
                indices.append(int(match.group(1)))
    return f"7.{max(indices) + 1}" if indices else "7.1"


def resolve_previous_stage70_from_run_dir(previous_run_dir: Path, pdb_id: str):
    checkpoint_dir = previous_run_dir / "checkpoints"
    search_dir = checkpoint_dir if checkpoint_dir.is_dir() else previous_run_dir
    pattern = re.compile(rf"^{re.escape(pdb_id)}\.stage_7\.(\d+)\.up$")
    candidates = []
    if search_dir.is_dir():
        for path in search_dir.glob(f"{pdb_id}.stage_7*.up"):
            match = pattern.fullmatch(path.name)
            if path.is_file() and match:
                candidates.append((int(match.group(1)), path.stat().st_mtime_ns, str(path), path))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    return candidates[0][3]


def autodetect_previous_stage70(args):
    if args.continue_stage_70_from:
        return args.continue_stage_70_from
    if args.previous_stage7_file:
        return args.previous_stage7_file
    if args.previous_run_dir:
        return resolve_previous_stage70_from_run_dir(args.previous_run_dir, args.pdb_id)
    if not args.auto_continue_from_previous_run:
        return None
    if not args.auto_continue_glob:
        return None
    pattern = re.compile(rf"^{re.escape(args.pdb_id)}\.stage_7\.(\d+)\.up$")
    candidates = []
    for path in (WORKFLOW_DIR / "outputs").glob(args.auto_continue_glob):
        match = pattern.fullmatch(path.name)
        if path.is_file() and match:
            candidates.append((int(match.group(1)), path.stat().st_mtime_ns, str(path), path))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    return candidates[0][3]


def resolve_continuation_outputs(args):
    source = autodetect_previous_stage70(args)
    if source is None:
        return None, None, None
    source = source.resolve()
    label = args.continue_stage_70_label
    if not label and args.continue_stage_70_output:
        label = infer_stage70_label_from_file(args.continue_stage_70_output, args.pdb_id)
        if not label:
            raise ValueError(f"CONTINUE_STAGE_70_OUTPUT must be named {args.pdb_id}.stage_7.N.up")
    if not label:
        label = infer_next_stage70_label(source, args.checkpoint_dir, args.pdb_id)
    if not re.fullmatch(r"7\.\d+", label):
        raise ValueError(f"Continuation stage label must be numeric stage_7.N, got {label}")
    output = args.continue_stage_70_output or (args.checkpoint_dir / f"{args.pdb_id}.stage_{label}.up")
    return source, output, label


def run_stage70_continuation(args, source_file: Path, output_file: Path, stage_label: str):
    if not source_file.exists():
        raise FileNotFoundError(source_file)
    assert_hybrid_stage_active(source_file, "production", "production")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if source_file.resolve() != output_file.resolve():
        shutil.copy2(source_file, output_file)
    handoff_initial_position(args, source_file, output_file, "production_hybrid")
    assert_hybrid_stage_active(output_file, "production", "production")
    run_md_stage(args, stage_label, output_file, output_file, args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
    extract_stage_vtf(args, stage_label, output_file, "2")


def normalize_hybrid_workflow_args(args):
    args.upside_home = workflow_path(args.upside_home, REPO_ROOT).resolve()
    args.run_dir = workflow_path(args.run_dir).resolve()
    args.checkpoint_dir = workflow_path(args.checkpoint_dir).resolve() if args.checkpoint_dir else args.run_dir / "checkpoints"
    args.log_dir = workflow_path(args.log_dir).resolve() if args.log_dir else args.run_dir / "logs"
    args.hybrid_prep_dir = workflow_path(args.hybrid_prep_dir).resolve() if args.hybrid_prep_dir else args.run_dir / "hybrid_prep"
    args.runtime_pdb_file = workflow_path(args.runtime_pdb_file).resolve() if args.runtime_pdb_file else args.hybrid_prep_dir / f"{args.runtime_pdb_id}.MARTINI.pdb"
    args.runtime_itp_file = workflow_path(args.runtime_itp_file).resolve() if args.runtime_itp_file else args.hybrid_prep_dir / f"{args.runtime_pdb_id}_proa.itp"
    args.hybrid_mapping_file = workflow_path(args.hybrid_mapping_file).resolve() if args.hybrid_mapping_file else args.hybrid_prep_dir / "hybrid_mapping.h5"
    args.hybrid_packed_pdb = workflow_path(args.hybrid_packed_pdb).resolve() if args.hybrid_packed_pdb else args.hybrid_prep_dir / "hybrid_packed.MARTINI.pdb"
    args.upside_executable = workflow_path(args.upside_executable, args.upside_home).resolve() if args.upside_executable else args.upside_home / "obj" / "upside"
    args.martini_ff_dir = workflow_path(args.martini_ff_dir, args.upside_home).resolve()
    args.mass_ff_file = workflow_path(args.mass_ff_file).resolve() if args.mass_ff_file else args.martini_ff_dir / "dry_martini_v2.1.itp"
    args.sc_martini_library = workflow_path(args.sc_martini_library, args.upside_home).resolve()
    args.sc_martini_table_json = workflow_path(args.sc_martini_table_json, args.upside_home).resolve()
    args.upside_rama_library = workflow_path(args.upside_rama_library, args.upside_home).resolve()
    args.upside_rama_sheet_mixing = workflow_path(args.upside_rama_sheet_mixing, args.upside_home).resolve()
    args.upside_hbond_energy = workflow_path(args.upside_hbond_energy, args.upside_home).resolve()
    args.upside_reference_state_rama = workflow_path(args.upside_reference_state_rama, args.upside_home).resolve()
    args.martinize_script = workflow_path(args.martinize_script).resolve()
    args.extract_vtf_script = workflow_path(args.extract_vtf_script).resolve()
    args.continue_stage_70_from = workflow_path(args.continue_stage_70_from).resolve() if args.continue_stage_70_from else None
    args.continue_stage_70_output = workflow_path(args.continue_stage_70_output).resolve() if args.continue_stage_70_output else None
    args.previous_stage7_file = workflow_path(args.previous_stage7_file).resolve() if args.previous_stage7_file else None
    args.previous_run_dir = workflow_path(args.previous_run_dir).resolve() if args.previous_run_dir else None
    if args.prep_seed is None:
        args.prep_seed = generate_random_seed()
    if args.seed is None:
        args.seed = generate_random_seed()
        if args.seed == args.prep_seed:
            args.seed = generate_random_seed()
    for path in [args.run_dir, args.checkpoint_dir, args.log_dir, args.hybrid_prep_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return args


def run_hybrid_workflow_command(argv):
    parser = argparse.ArgumentParser(description="Run the hybrid 1RKL dry-MARTINI workflow.")
    parser.add_argument("--pdb-id", default=env_default("PDB_ID", "1rkl"))
    parser.add_argument("--runtime-pdb-id", default=env_default("RUNTIME_PDB_ID", None))
    parser.add_argument("--upside-home", default=env_default("UPSIDE_HOME", str(REPO_ROOT)))
    parser.add_argument("--run-dir", default=env_default("RUN_DIR", "outputs/martini_test_1rkl_hybrid"))
    parser.add_argument("--checkpoint-dir", default=env_default("CHECKPOINT_DIR", None))
    parser.add_argument("--log-dir", default=env_default("LOG_DIR", None))
    parser.add_argument("--hybrid-prep-dir", default=env_default("HYBRID_PREP_DIR", None))
    parser.add_argument("--protein-aa-pdb", default=env_default("PROTEIN_AA_PDB", None))
    parser.add_argument("--protein-cg-pdb", default=env_default("PROTEIN_CG_PDB", ""))
    parser.add_argument("--protein-itp", default=env_default("PROTEIN_ITP", ""))
    parser.add_argument("--bilayer-pdb", default=env_default("BILAYER_PDB", None))
    parser.add_argument("--runtime-pdb-file", default=env_default("RUNTIME_PDB_FILE", None))
    parser.add_argument("--runtime-itp-file", default=env_default("RUNTIME_ITP_FILE", None))
    parser.add_argument("--hybrid-mapping-file", default=env_default("HYBRID_MAPPING_FILE", None))
    parser.add_argument("--hybrid-packed-pdb", default=env_default("HYBRID_PACKED_PDB", None))
    parser.add_argument("--universal-prep-mode", default=env_default("UNIVERSAL_PREP_MODE", "both"))
    parser.add_argument("--hybrid-validate", type=int, default=env_int("HYBRID_VALIDATE", 1))
    parser.add_argument("--hybrid-preprod-activation-stage", default=env_default("HYBRID_PREPROD_ACTIVATION_STAGE", "__hybrid_disabled__"))
    parser.add_argument("--martinize-enable", type=int, default=env_int("MARTINIZE_ENABLE", 1))
    parser.add_argument("--martinize-ff", default=env_default("MARTINIZE_FF", "martini22"))
    parser.add_argument("--martinize-molname", default=env_default("MARTINIZE_MOLNAME", "PROA"))
    parser.add_argument("--martinize-script", default=env_default("MARTINIZE_SCRIPT", str(PY_DIR / "martini_martinize.py")))
    parser.add_argument("--extract-vtf-script", default=env_default("EXTRACT_VTF_SCRIPT", str(PY_DIR / "martini_extract_vtf.py")))
    parser.add_argument("--salt-molar", type=float, default=env_float("SALT_MOLAR", 0.15))
    parser.add_argument("--protein-lipid-cutoff", type=float, default=env_float("PROTEIN_LIPID_CUTOFF", 4.5))
    parser.add_argument("--ion-cutoff", type=float, default=env_float("ION_CUTOFF", 4.0))
    parser.add_argument("--xy-scale", type=float, default=env_float("XY_SCALE", 1.0))
    parser.add_argument("--box-padding-xy", type=float, default=env_float("BOX_PADDING_XY", 0.0))
    parser.add_argument("--box-padding-z", type=float, default=env_float("BOX_PADDING_Z", 20.0))
    parser.add_argument("--protein-placement-mode", choices=["embed", "outside-top", "outside-bottom"], default=env_default("PROTEIN_PLACEMENT_MODE", "embed"))
    parser.add_argument("--protein-orientation-mode", choices=["input", "lay-flat"], default=env_default("PROTEIN_ORIENTATION_MODE", "input"))
    parser.add_argument("--protein-surface-gap", type=float, default=env_float("PROTEIN_SURFACE_GAP", 6.0))
    parser.add_argument("--protein-lipid-min-gap", type=float, default=env_float("PROTEIN_LIPID_MIN_GAP", 4.5))
    parser.add_argument("--protein-lipid-cutoff-step", type=float, default=env_float("PROTEIN_LIPID_CUTOFF_STEP", 0.5))
    parser.add_argument("--protein-lipid-cutoff-max", type=float, default=env_float("PROTEIN_LIPID_CUTOFF_MAX", 8.0))
    parser.add_argument("--temperature", type=float, default=env_float("TEMPERATURE", 0.8647))
    parser.add_argument("--thermostat-timescale", type=float, default=env_float("THERMOSTAT_TIMESCALE", 5.0))
    parser.add_argument("--thermostat-interval", type=int, default=env_int("THERMOSTAT_INTERVAL", -1))
    parser.add_argument("--strict-stage-handoff", type=int, default=env_int("STRICT_STAGE_HANDOFF", 1))
    parser.add_argument("--min-60-max-iter", type=int, default=env_int("MIN_60_MAX_ITER", 500))
    parser.add_argument("--min-61-max-iter", type=int, default=env_int("MIN_61_MAX_ITER", 500))
    parser.add_argument("--eq-62-nsteps", type=int, default=env_int("EQ_62_NSTEPS", 500))
    parser.add_argument("--eq-63-nsteps", type=int, default=env_int("EQ_63_NSTEPS", 500))
    parser.add_argument("--eq-64-nsteps", type=int, default=env_int("EQ_64_NSTEPS", 500))
    parser.add_argument("--eq-65-nsteps", type=int, default=env_int("EQ_65_NSTEPS", 500))
    parser.add_argument("--eq-66-nsteps", type=int, default=env_int("EQ_66_NSTEPS", 500))
    parser.add_argument("--prod-70-nsteps", type=int, default=env_int("PROD_70_NSTEPS", 10000))
    parser.add_argument("--eq-time-step", type=float, default=env_float("EQ_TIME_STEP", 0.010))
    parser.add_argument("--prod-time-step", type=float, default=env_float("PROD_TIME_STEP", 0.002))
    parser.add_argument("--min-time-step", type=float, default=env_float("MIN_TIME_STEP", 0.010))
    parser.add_argument("--eq-frame-steps", type=int, default=env_int("EQ_FRAME_STEPS", 1000))
    parser.add_argument("--prod-frame-steps", type=int, default=env_int("PROD_FRAME_STEPS", 50))
    parser.add_argument("--prod-70-npt-enable", type=int, default=env_int("PROD_70_NPT_ENABLE", 0))
    parser.add_argument("--prod-70-backbone-fix-rigid-enable", type=int, default=env_int("PROD_70_BACKBONE_FIX_RIGID_ENABLE", 0))
    parser.add_argument("--prep-seed", default=os.environ.get("PREP_SEED"))
    parser.add_argument("--seed", default=os.environ.get("SEED"))
    parser.add_argument("--continue-stage-70-from", default=env_default("CONTINUE_STAGE_70_FROM", ""))
    parser.add_argument("--continue-stage-70-output", default=env_default("CONTINUE_STAGE_70_OUTPUT", ""))
    parser.add_argument("--continue-stage-70-label", default=env_default("CONTINUE_STAGE_70_LABEL", ""))
    parser.add_argument("--previous-run-dir", default=env_default("PREVIOUS_RUN_DIR", ""))
    parser.add_argument("--previous-stage7-file", default=env_default("PREVIOUS_STAGE7_FILE", ""))
    parser.add_argument("--auto-continue-from-previous-run", type=int, default=env_int("AUTO_CONTINUE_FROM_PREVIOUS_RUN", 0))
    parser.add_argument("--auto-continue-glob", default=env_default("AUTO_CONTINUE_GLOB", ""))
    parser.add_argument("--upside-executable", default=env_default("UPSIDE_EXECUTABLE", None))
    parser.add_argument("--martini-ff-dir", default=env_default("UPSIDE_MARTINI_FF_DIR", str(REPO_ROOT / "parameters" / "dryMARTINI")))
    parser.add_argument("--mass-ff-file", default=env_default("MASS_FF_FILE", ""))
    parser.add_argument("--sc-martini-library", default=env_default("SC_MARTINI_LIBRARY", "parameters/ff_2.1/martini.h5"))
    parser.add_argument("--sc-martini-table-json", default=env_default("SC_MARTINI_TABLE_JSON", "SC-training/runs/default/results/assembled/sc_table.json"))
    parser.add_argument("--upside-rama-library", default=env_default("UPSIDE_RAMA_LIBRARY", "parameters/common/rama.dat"))
    parser.add_argument("--upside-rama-sheet-mixing", default=env_default("UPSIDE_RAMA_SHEET_MIXING", "parameters/ff_2.1/sheet"))
    parser.add_argument("--upside-hbond-energy", default=env_default("UPSIDE_HBOND_ENERGY", "parameters/ff_2.1/hbond.h5"))
    parser.add_argument("--upside-reference-state-rama", default=env_default("UPSIDE_REFERENCE_STATE_RAMA", "parameters/common/rama_reference.pkl"))
    parser.add_argument("--sc-env-lj-force-cap", type=float, default=25.0)
    parser.add_argument("--sc-env-coul-force-cap", type=float, default=25.0)
    parser.add_argument("--nonprotein-hs-force-cap", type=float, default=100.0)
    parser.add_argument("--sc-env-po4-z-clamp-enable", type=int, default=1)
    parser.add_argument("--sc-env-relax-steps", type=int, default=150)
    parser.add_argument("--sc-env-backbone-hold-steps", type=int, default=200)
    parser.add_argument("--sc-env-po4-z-hold-steps", type=int, default=150)
    parser.add_argument("--bb-aa-max-rigid-rmsd", type=float, default=1.5)
    parser.add_argument("--bb-aa-min-matched-residues", type=int, default=8)
    parser.add_argument("--npt-tau", type=float, default=4.0)
    parser.add_argument("--npt-interval", type=int, default=10)
    parser.add_argument("--ewald-enable", type=int, default=1)
    parser.add_argument("--ewald-alpha", type=float, default=0.2)
    parser.add_argument("--ewald-kmax", type=int, default=5)
    parser.add_argument("--prod-70-barostat-type", type=int, default=1)
    parser.add_argument("--martini-energy-conversion", type=float, default=2.914952774272)
    parser.add_argument("--martini-length-conversion", type=float, default=10.0)
    parser.add_argument("--bar-1-to-eup-per-a3", type=float, default=0.000020659477)
    parser.add_argument("--compressibility-3e4-bar-inv-to-a3-per-eup", type=float, default=14.521180763676)
    parser.add_argument("--protein-env-interface-scale", type=float, default=1.0)
    args = parser.parse_args(argv)
    if args.runtime_pdb_id is None:
        args.runtime_pdb_id = f"{args.pdb_id}_hybrid"
    if args.protein_aa_pdb is None:
        args.protein_aa_pdb = f"pdb/{args.pdb_id}.pdb"
    if args.bilayer_pdb is None:
        args.bilayer_pdb = str(Path(args.upside_home) / "parameters" / "dryMARTINI" / "DOPC.pdb")
    args.martinize_enable = bool(args.martinize_enable)
    args.hybrid_validate = bool(args.hybrid_validate)
    args.prod_70_backbone_fix_rigid_enable = bool(args.prod_70_backbone_fix_rigid_enable)
    args.auto_continue_from_previous_run = bool(args.auto_continue_from_previous_run)
    args.prep_seed = int(args.prep_seed) if args.prep_seed not in (None, "") else None
    args.seed = int(args.seed) if args.seed not in (None, "") else None
    args = normalize_hybrid_workflow_args(args)

    if not args.upside_executable.exists():
        raise FileNotFoundError(args.upside_executable)
    for required in [
        workflow_path(args.protein_aa_pdb).resolve(),
        workflow_path(args.bilayer_pdb).resolve(),
        args.upside_rama_library,
        args.upside_hbond_energy,
        args.upside_reference_state_rama,
        args.martinize_script,
        args.extract_vtf_script,
    ]:
        if not required.exists():
            raise FileNotFoundError(required)

    source, output, label = resolve_continuation_outputs(args)
    print("=== Hybrid Dry MARTINI Workflow ===")
    print(f"Protein ID: {args.pdb_id}")
    print(f"Runtime PDB ID: {args.runtime_pdb_id}")
    print(f"Preparation seed: {args.prep_seed}")
    print(f"Simulation seed: {args.seed}")
    print(f"Run directory: {args.run_dir}")
    if source:
        print("Continuation mode: production only")
        print(f"Continuation source: {source}")
        print(f"Continuation output: {output}")
        run_stage70_continuation(args, source, output, label)
        return

    protein_cg, protein_itp = prepare_workflow_protein_inputs(args)
    prepare_workflow_hybrid_artifacts(args, protein_cg, protein_itp)

    files = {
        "prepared_60": args.checkpoint_dir / f"{args.pdb_id}.stage_6.0.prepared.up",
        "stage_60": args.checkpoint_dir / f"{args.pdb_id}.stage_6.0.up",
        "prepared_61": args.checkpoint_dir / f"{args.pdb_id}.stage_6.1.prepared.up",
        "stage_61": args.checkpoint_dir / f"{args.pdb_id}.stage_6.1.up",
        "stage_62": args.checkpoint_dir / f"{args.pdb_id}.stage_6.2.up",
        "prepared_63": args.checkpoint_dir / f"{args.pdb_id}.stage_6.3.prepared.up",
        "stage_63": args.checkpoint_dir / f"{args.pdb_id}.stage_6.3.up",
        "prepared_64": args.checkpoint_dir / f"{args.pdb_id}.stage_6.4.prepared.up",
        "stage_64": args.checkpoint_dir / f"{args.pdb_id}.stage_6.4.up",
        "prepared_65": args.checkpoint_dir / f"{args.pdb_id}.stage_6.5.prepared.up",
        "stage_65": args.checkpoint_dir / f"{args.pdb_id}.stage_6.5.up",
        "prepared_66": args.checkpoint_dir / f"{args.pdb_id}.stage_6.6.prepared.up",
        "stage_66": args.checkpoint_dir / f"{args.pdb_id}.stage_6.6.up",
        "prepared_70": args.checkpoint_dir / f"{args.pdb_id}.stage_7.0.prepared.up",
        "stage_70": args.checkpoint_dir / f"{args.pdb_id}.stage_7.0.up",
    }

    prepare_stage_file(args, files["prepared_60"], "minimization", 1, 0, 0, "minimization", protein_itp)
    shutil.copy2(files["prepared_60"], files["stage_60"])
    run_minimization_stage(args, "6.0", files["stage_60"], args.min_60_max_iter)
    extract_stage_vtf(args, "6.0", files["stage_60"], "1")

    prepare_stage_file(args, files["prepared_61"], "npt_prod", 1, 0, 0, "minimization", protein_itp)
    shutil.copy2(files["prepared_61"], files["stage_61"])
    handoff_initial_position(args, files["stage_60"], files["stage_61"])
    run_minimization_stage(args, "6.1", files["stage_61"], args.min_61_max_iter)
    extract_stage_vtf(args, "6.1", files["stage_61"], "1")

    prepare_stage_file(args, files["stage_62"], "npt_equil", 1, 0, 200, "minimization", protein_itp)
    handoff_initial_position(args, files["stage_61"], files["stage_62"])
    run_md_stage(args, "6.2", files["stage_62"], files["stage_62"], args.eq_62_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.2", files["stage_62"], "1")

    prepare_stage_file(args, files["prepared_63"], "npt_equil_reduced", 1, 0, 100, "minimization", protein_itp)
    shutil.copy2(files["prepared_63"], files["stage_63"])
    handoff_initial_position(args, files["stage_62"], files["stage_63"])
    run_md_stage(args, "6.3", files["stage_63"], files["stage_63"], args.eq_63_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.3", files["stage_63"], "1")

    prepare_stage_file(args, files["prepared_64"], "npt_prod", 1, 0, 50, "minimization", protein_itp)
    shutil.copy2(files["prepared_64"], files["stage_64"])
    handoff_initial_position(args, files["stage_63"], files["stage_64"])
    run_md_stage(args, "6.4", files["stage_64"], files["stage_64"], args.eq_64_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.4", files["stage_64"], "1")

    prepare_stage_file(args, files["prepared_65"], "npt_prod", 1, 0, 20, "minimization", protein_itp)
    shutil.copy2(files["prepared_65"], files["stage_65"])
    handoff_initial_position(args, files["stage_64"], files["stage_65"])
    run_md_stage(args, "6.5", files["stage_65"], files["stage_65"], args.eq_65_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.5", files["stage_65"], "1")

    prepare_stage_file(args, files["prepared_66"], "npt_prod", 1, 0, 10, "minimization", protein_itp)
    shutil.copy2(files["prepared_66"], files["stage_66"])
    handoff_initial_position(args, files["stage_65"], files["stage_66"])
    run_md_stage(args, "6.6", files["stage_66"], files["stage_66"], args.eq_66_nsteps, args.eq_time_step, args.eq_frame_steps)
    extract_stage_vtf(args, "6.6", files["stage_66"], "1")

    prepare_stage_file(args, files["prepared_70"], "npt_prod", args.prod_70_npt_enable, args.prod_70_barostat_type, 0, "production", protein_itp)
    shutil.copy2(files["prepared_70"], files["stage_70"])
    handoff_initial_position(args, files["stage_66"], files["stage_70"], "production_hybrid")
    assert_hybrid_stage_active(files["stage_70"], "production", "production")
    run_md_stage(args, "7.0", files["stage_70"], files["stage_70"], args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
    extract_stage_vtf(args, "7.0", files["stage_70"], "2")
    print("=== Workflow Complete ===")


if __name__ == "__main__":
    command_handlers = {
        "build-sc-martini-h5": run_build_sc_martini_h5_command,
        "inject-stage7-sc": run_inject_stage7_sc_command,
        "run-hybrid-workflow": run_hybrid_workflow_command,
        "set-initial-position": run_set_initial_position_command,
        "validate-hybrid-mapping": run_validate_hybrid_mapping_command,
    }
    argv = sys.argv[1:]
    if argv and argv[0] in command_handlers:
        command_handlers[argv[0]](argv[1:])
    else:
        run_prepare_command(argv)
