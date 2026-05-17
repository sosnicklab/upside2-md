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

from martini_build_tables import build_martini_tables
from martini_prepare_system_lib import (
    build_backbone_with_virtual_bb,
    center_of_mass,
    convert_stage,
    compute_lipid_residue_indices,
    coords,
    estimate_salt_pairs,
    extract_protein_backbone_atoms_from_aa,
    infer_effective_ion_volume_fraction_from_template,
    infer_protein_charge_from_residues,
    inject_particles_table,
    inject_cg_lipid_nodes,
    inject_stage7_sc_table_nodes,
    lipid_resname,
    map_backbone_types_from_martinize_fallback,
    parse_pdb,
    place_ions,
    remove_overlapping_lipids,
    set_initial_position,
    set_box_from_lipid_xy,
    set_coords,
    tile_and_crop_bilayer_lipids,
    validate_hybrid_mapping,
    write_hybrid_mapping_h5,
    write_pdb,
    write_stage_debug_outputs,
)


PY_DIR = Path(__file__).resolve().parent
REPO_ROOT = PY_DIR.parent
WORKFLOW_DIR = REPO_ROOT / "example" / "16.MARTINI"

DEFAULT_SC_ENV_LJ_FORCE_CAP = 25.0
DEFAULT_SC_ENV_COUL_FORCE_CAP = 25.0
DEFAULT_NONPROTEIN_HS_FORCE_CAP = 100.0
DEFAULT_SC_ENV_PO4_Z_CLAMP_ENABLE = 1
DEFAULT_SC_ENV_RELAX_STEPS = 150
DEFAULT_SC_ENV_BACKBONE_HOLD_STEPS = 200
DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS = 150
DEFAULT_NPT_TAU = 4.0
DEFAULT_NPT_INTERVAL = 10
DEFAULT_PROD_70_BAROSTAT_TYPE = 1
DEFAULT_MARTINI_ENERGY_CONVERSION = 2.914952774272
DEFAULT_MARTINI_LENGTH_CONVERSION = 10.0
DEFAULT_BAR_1_TO_EUP_PER_A3 = 0.000020659477
DEFAULT_COMPRESSIBILITY_3E4_BAR_INV_TO_A3_PER_EUP = 14.521180763676
DEFAULT_PROTEIN_ENV_INTERFACE_SCALE = 1.0


def parse_prepare_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Unified preparation script for bilayer-only, protein-only, or mixed "
            "protein+bilayer systems. Can optionally convert prepared structure "
            "to UPSIDE input."
        )
    )
    parser.add_argument("--mode", choices=["both"], required=True)
    parser.add_argument("--pdb-id", required=True, help="Runtime PDB id for stage conversion")
    parser.add_argument("--runtime-pdb-output", default=None)
    parser.add_argument("--prepare-structure", type=int, default=1, choices=[0, 1])
    parser.add_argument("--stage", default=None, help="stage name for UPSIDE input conversion")
    parser.add_argument("--run-dir", default="outputs/martini_test")

    parser.add_argument("--bilayer-pdb", default=str(REPO_ROOT / "parameters" / "dryMARTINI" / "DOPC.pdb"))
    parser.add_argument("--protein-aa-pdb", default=None)
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
    parser.add_argument("--summary-json", default=None)
    return parser.parse_args(argv)


def runtime_paths(args):
    runtime_pdb = (
        Path(args.runtime_pdb_output).expanduser().resolve()
        if args.runtime_pdb_output
        else (WORKFLOW_DIR / "pdb" / f"{args.pdb_id}.MARTINI.pdb")
    )
    return runtime_pdb


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


def prepare_mixed_structure(args, runtime_pdb):
    if not args.protein_aa_pdb:
        raise ValueError("--protein-aa-pdb is required for mode=both")
    if not args.bilayer_pdb:
        raise ValueError("--bilayer-pdb is required for mode=both")

    protein_aa_pdb = Path(args.protein_aa_pdb).expanduser().resolve()
    bilayer_pdb = Path(args.bilayer_pdb).expanduser().resolve()
    if not protein_aa_pdb.exists():
        raise FileNotFoundError(f"Protein AA PDB not found: {protein_aa_pdb}")
    if not bilayer_pdb.exists():
        raise FileNotFoundError(f"Bilayer PDB not found: {bilayer_pdb}")

    protein_aa_atoms, _ = parse_pdb(protein_aa_pdb)
    if not protein_aa_atoms:
        raise ValueError(f"No atoms found in protein AA PDB: {protein_aa_pdb}")
    protein_backbone_atoms, _ = extract_protein_backbone_atoms_from_aa(protein_aa_atoms)
    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    bilayer_lipid_atoms = [a for a in bilayer_atoms if lipid_resname(a["resname"])]
    if not bilayer_lipid_atoms:
        raise ValueError("No lipid residues found in bilayer template.")

    protein_xyz_raw = coords(protein_backbone_atoms)
    bilayer_xyz = coords(bilayer_lipid_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    bmin = bilayer_xyz.min(axis=0)
    bmax = bilayer_xyz.max(axis=0)
    base_side = max(float(bmax[0] - bmin[0]), float(bmax[1] - bmin[1]))

    if args.protein_placement_mode == "embed" and args.protein_orientation_mode == "input":
        protein_center = center_of_mass(protein_xyz_raw)
        shift = bilayer_center - protein_center
        protein_xyz = protein_xyz_raw + shift
        set_coords(protein_backbone_atoms, protein_xyz)

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
        set_coords(protein_backbone_atoms, protein_xyz)
    lipid_residues, keep_nonlipid = compute_lipid_residue_indices(bilayer_lipids)
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_lipids,
        protein_atoms=protein_backbone_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=float(args.protein_lipid_cutoff),
    )

    packed_atoms = protein_backbone_atoms + bilayer_kept
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
        else int(infer_protein_charge_from_residues(protein_backbone_atoms))
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

    bb_type_by_residue = map_backbone_types_from_martinize_fallback(protein_backbone_atoms)
    protein_atoms, bb_entries = build_backbone_with_virtual_bb(protein_backbone_atoms, bb_type_by_residue)
    all_atoms = protein_atoms + bilayer_kept + ion_atoms
    write_pdb(runtime_pdb, all_atoms, box_lengths)
    protein_xyz_final = coords(protein_backbone_atoms)
    bilayer_kept_lipid_atoms = [a for a in bilayer_kept if lipid_resname(a["resname"])]
    bilayer_xyz_final = coords(bilayer_kept_lipid_atoms) if bilayer_kept_lipid_atoms else coords(bilayer_kept)

    mapping_summary = {}
    if args.hybrid_mapping_output:
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
        mapping_summary["bb_fallback_typed_residues"] = int(len(bb_type_by_residue))
        mapping_summary["bb_map_entries"] = int(len(bb_entries))

        if args.hybrid_bb_map_json_output:
            mapping_json = Path(args.hybrid_bb_map_json_output).expanduser().resolve()
            write_summary(mapping_json, {"bb_entries": bb_entries, "count": len(bb_entries)})
            mapping_summary["mapping_json"] = str(mapping_json)

    summary = {
        "mode": "both",
        "input_protein_aa_pdb": str(protein_aa_pdb),
        "input_bilayer_pdb": str(bilayer_pdb),
        "runtime_pdb": str(runtime_pdb),
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


def run_stage_conversion(args, runtime_pdb: Path):
    prev_pdb = os.environ.get("UPSIDE_RUNTIME_PDB_FILE")

    os.environ["UPSIDE_RUNTIME_PDB_FILE"] = str(runtime_pdb)
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

        os.environ.pop("UPSIDE_RUNTIME_ITP_FILE", None)


def write_summary(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def run_prepare_command(argv):
    args = parse_prepare_args(argv)
    runtime_pdb = runtime_paths(args)
    runtime_pdb.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "mode": args.mode,
        "pdb_id": args.pdb_id,
        "prepare_structure": bool(args.prepare_structure),
        "stage": args.stage,
        "run_dir": args.run_dir,
    }

    if args.prepare_structure:
        summary.update(prepare_mixed_structure(args, runtime_pdb))
    else:
        if not runtime_pdb.exists():
            raise FileNotFoundError(
                f"Runtime PDB not found for stage conversion: {runtime_pdb}. "
                "Run with --prepare-structure 1 first."
            )

    if args.stage:
        run_stage_conversion(args, runtime_pdb)
        summary["upside_input"] = str(Path(args.run_dir).expanduser().resolve() / "test.input.up")

    if args.summary_json:
        summary_path = Path(args.summary_json).expanduser().resolve()
    else:
        summary_path = runtime_pdb.with_suffix(".prep_summary.json")
    write_summary(summary_path, summary)
    print(f"Preparation summary written to: {summary_path}")


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


def prepare_workflow_hybrid_artifacts(args):
    print("=== Stage 0: Hybrid Packing + Mapping Export ===")
    cutoff = float(args.protein_lipid_cutoff)
    while True:
        print(f"Hybrid packing attempt with protein-lipid cutoff: {cutoff:.3f} A")
        run_prepare_command(
            [
                "--mode", "both",
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.hybrid_packed_pdb),
                "--prepare-structure", "1",
                "--protein-aa-pdb", str(workflow_path(args.protein_aa_pdb).resolve()),
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
    print(f"Runtime MARTINI PDB: {args.runtime_pdb_file}")


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

    with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
        src_inp = src["/input"]
        dst_inp = dst.require_group("input")
        base_n_atom = int(dst_inp["pos"].shape[0])
        src_mem = src_inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32)

        # If CG lipid coarse-graining reduced the atom count, remap the membership
        old_to_new = None
        if base_n_atom != int(src_mem.shape[0]):
            remap_grp = dst_inp.get("hybrid_remap")
            if remap_grp is not None and "old_to_new" in remap_grp:
                old_to_new = remap_grp["old_to_new"][:].astype(np.int32)
                n_old = int(remap_grp.attrs["n_old"])
                if n_old == src_mem.shape[0]:
                    print(f"  Remapping hybrid membership from {n_old} → {base_n_atom} atoms "
                          f"(CG lipid coarse-graining)")
                    new_membership = np.full(base_n_atom, -1, dtype=np.int32)
                    for old_idx in range(n_old):
                        new_idx = old_to_new[old_idx]
                        new_membership[new_idx] = src_mem[old_idx]
                    src_mem = new_membership
                else:
                    raise ValueError(
                        f"CG remap n_old={n_old} doesn't match mapping "
                        f"n_atom={src_mem.shape[0]}"
                    )
            else:
                raise ValueError(
                    f"Hybrid mapping n_atom mismatch for {up_file}: up has {base_n_atom}, "
                    f"mapping has {src_mem.shape[0]}, and no CG lipid remap data found"
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

        # If we remapped membership, fix up the stored membership in the output
        if old_to_new is not None and "protein_membership" in dst_inp["hybrid_env_topology"]:
            del dst_inp["hybrid_env_topology"]["protein_membership"]
            dst_inp["hybrid_env_topology"].create_dataset(
                "protein_membership", data=src_mem.astype(np.int32))

        bb_grp = dst_inp["hybrid_bb_map"]
        env_grp = dst_inp["hybrid_env_topology"]
        if bb_grp["atom_indices"].shape[0] != bb_grp["bb_residue_index"].shape[0]:
            raise ValueError("hybrid_bb_map atom_indices/bb_residue_index size mismatch")
        if bb_grp["atom_indices"].shape[1] != 4:
            raise ValueError("hybrid_bb_map/atom_indices must have shape (n_bb,4)")

        # Remap BB atom_indices through old_to_new if CG lipids changed atom ordering
        if old_to_new is not None:
            bb_indices = bb_grp["atom_indices"][:].astype(np.int32)
            bb_remapped = old_to_new[bb_indices]
            del dst_inp["hybrid_bb_map"]["atom_indices"]
            dst_inp["hybrid_bb_map"].create_dataset("atom_indices", data=bb_remapped.astype(np.int32))

        bb_grp.attrs["atom_index_space"] = np.bytes_("stage_runtime")
        bb_grp.attrs["reference_index_space"] = np.bytes_("stage_runtime")
        bb_grp.attrs["reference_index_offset"] = np.int32(0)
        bb_grp.attrs["reference_index_count"] = np.int32(0)
        membership = env_grp["protein_membership"][:].astype(np.int32, copy=False)
        if membership.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch")


def set_stage_label(up_file: Path, stage_label: str):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        grp = h5.require_group("input").require_group("stage_parameters")
        grp.attrs["enable"] = np.int8(1)
        grp.attrs["current_stage"] = np.bytes_(stage_label)


def set_hybrid_control_mode(up_file: Path, activation_stage: str, preprod_mode="rigid_body"):
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


def set_debug_rigid_protein(up_file: Path):
    import h5py

    with h5py.File(up_file, "r+") as h5:
        inp = h5["/input"]
        atom_indices = np.zeros(0, dtype=np.int32)
        if "particle_class" in inp:
            classes = np.asarray([
                x.decode("ascii") if isinstance(x, (bytes, np.bytes_)) else str(x)
                for x in inp["particle_class"][:]
            ], dtype=object)
            atom_indices = np.where(np.char.upper(classes.astype(str)) == "PROTEIN")[0].astype(np.int32)
        if atom_indices.size == 0 and "hybrid_env_topology/protein_membership" in inp:
            membership = np.asarray(inp["hybrid_env_topology/protein_membership"][:], dtype=np.int32)
            atom_indices = np.where(membership >= 0)[0].astype(np.int32)

        if atom_indices.size == 0:
            raise ValueError(f"{up_file}: DEBUG_RIGID_PROTEIN requested but no protein atoms were found")

        if "fix_rigid" in inp:
            del inp["fix_rigid"]
        grp = inp.create_group("fix_rigid")
        grp.attrs["enable"] = np.int8(1)
        grp.create_dataset("atom_indices", data=atom_indices, dtype=np.int32)
        print(f"DEBUG_RIGID_PROTEIN: fixed {atom_indices.size} protein atoms in {up_file}")


def set_barostat_type(up_file: Path, barostat_type: int):
    import tables as tb

    with tb.open_file(up_file, "r+") as t:
        if "/input/barostat" in t:
            t.root.input.barostat._v_attrs.type = int(barostat_type)


def ensure_sc_martini_library(args):
    import h5py

    library = args.sc_martini_library
    if not library.exists():
        raise FileNotFoundError(library)
    with h5py.File(library, "r") as h5:
        sc_grp = h5["sc_table"] if "sc_table" in h5 else h5
        missing = [name for name in SC_LIBRARY_REQUIRED_DATASETS if name not in sc_grp]
    if missing:
        raise ValueError(f"{library} missing required SC datasets: {','.join(missing)}")


def assert_hybrid_stage_active(
    up_file: Path,
    expected_stage: str,
    expected_activation: str,
    require_interface_nodes: bool = False,
):
    import h5py

    with h5py.File(up_file, "r") as h5:
        inp = h5["/input"]
        pot = inp["potential"]
        stage = h5_as_text(inp["stage_parameters"].attrs.get("current_stage", b"")).strip()
        hy = inp["hybrid_control"].attrs
        enable = int(hy.get("enable", 0))
        activation = h5_as_text(hy.get("activation_stage", b"")).strip()
        if stage != expected_stage or enable != 1 or activation != expected_activation:
            raise ValueError(
                f"{up_file}: expected stage={expected_stage}, activation={expected_activation}, enable=1; "
                f"got stage={stage}, activation={activation}, enable={enable}"
            )
        if require_interface_nodes:
            missing = [
                node
                for node in ("martini_potential", "martini_sc_table_1body")
                if node not in pot
            ]
            atom_names = []
            if "atom_names" in inp:
                atom_names = [
                    x.decode("ascii") if isinstance(x, (bytes, np.bytes_)) else str(x)
                    for x in inp["atom_names"][:]
                ]
            if any(name.strip().upper() == "CGL" for name in atom_names) and "cg_lipid_sc" not in pot:
                missing.append("cg_lipid_sc")
            if missing:
                raise ValueError(
                    f"{up_file}: missing required hybrid interface node(s): {', '.join(missing)}"
                )
        env_membership = inp["hybrid_env_topology"]["protein_membership"][:]
        if not np.any(env_membership < 0):
            raise ValueError(f"{up_file}: no non-protein environment targets found")
        if "type" in inp:
            atom_types = np.asarray([
                h5_as_text(x).strip().upper() for x in inp["type"][:]
            ], dtype=object)
            lipid_internal = np.isin(atom_types, ["CGL", "CGLD"])
            if np.any(lipid_internal & (env_membership >= 0)):
                bad = np.where(lipid_internal & (env_membership >= 0))[0][:8]
                raise ValueError(
                    f"{up_file}: CGL/CGLD sites must not have protein membership; "
                    f"bad indices={bad.tolist()}"
                )
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
        "UPSIDE_BILAYER_LIPIDHEAD_FC": str(lipidhead_fc),
        "UPSIDE_WRITE_DEBUG_PDB": "0",
    }


def inject_hybrid_interface_nodes(args, target_file: Path, current_stage: str, activation_stage: str):
    ensure_sc_martini_library(args)
    set_hybrid_control_mode(target_file, activation_stage)
    set_hybrid_production_controls(target_file, args)
    inject_stage7_sc_table_nodes(
        up_file=target_file,
        martini_h5=args.sc_martini_library,
        upside_home=args.upside_home,
        rama_library=args.upside_rama_library,
        rama_sheet_mixing=args.upside_rama_sheet_mixing,
        hbond_energy=args.upside_hbond_energy,
        reference_state_rama=args.upside_reference_state_rama,
    )
    # SC placement node must exist before CG lipid node injection
    # so that cg_lipid_sc (CG↔SC) can find it.
    inject_cg_lipid_nodes(up_file=target_file, martini_h5=args.martini_table_h5)
    assert_hybrid_stage_active(
        target_file,
        current_stage,
        activation_stage,
        require_interface_nodes=True,
    )


def prepare_stage_file(args, target_file: Path, prepare_stage: str, npt_enable: int, barostat_type: int, lipidhead_fc: float, stage_label: str):
    with temporary_env(stage_conversion_env(args, stage_label, prepare_stage, npt_enable, lipidhead_fc)):
        run_prepare_command(
            [
                "--mode", args.universal_prep_mode,
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.runtime_pdb_file),
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
    inject_particles_table(up_file=target_file, martini_h5=args.martini_table_h5)
    if stage_label == "production":
        inject_hybrid_interface_nodes(args, target_file, "production", "production")
    else:
        inject_hybrid_interface_nodes(
            args,
            target_file,
            stage_label,
            args.hybrid_preprod_activation_stage,
        )
    if npt_enable:
        set_barostat_type(target_file, barostat_type)
    write_stage_debug_outputs(target_file, debug_dir=args.run_dir / "debug")


def write_stage70_initial_debug_files(args):
    print("=== Stage 7.0 Initial Debug PDB Export ===")
    with temporary_env({
        **stage_conversion_env(args, "7.0", "npt_prod", args.prod_70_npt_enable, 0),
        "UPSIDE_INITIAL_DEBUG_ONLY": "1",
        "UPSIDE_INITIAL_DEBUG_PREFIX": f"{args.pdb_id}.stage_7.0",
    }):
        run_prepare_command(
            [
                "--mode", args.universal_prep_mode,
                "--pdb-id", args.runtime_pdb_id,
                "--runtime-pdb-output", str(args.runtime_pdb_file),
                "--prepare-structure", "0",
                "--stage", "npt_prod",
                "--run-dir", str(args.run_dir),
                "--summary-json", str(args.hybrid_prep_dir / "stage_7.0.initial_debug.summary.json"),
            ]
        )


def handoff_initial_position(args, input_file: Path, output_file: Path, mode="default", previous_dt=None):
    preserve_transition = "1" if mode == "production_restart" and previous_dt is not None else "0"
    with temporary_env(
        {
            "UPSIDE_SET_INITIAL_STRICT_COPY": str(args.strict_stage_handoff),
            "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION": "0",
            "UPSIDE_SET_INITIAL_PRESERVE_HYBRID_TRANSITION": preserve_transition,
            "UPSIDE_SET_INITIAL_TIME_STEP": str(previous_dt or 0),
        }
    ):
        set_initial_position(input_file, output_file)


def input_momentum_restart_valid(up_file: Path):
    import h5py

    with h5py.File(up_file, "r") as h5:
        if "/input/mom" not in h5:
            return False
        return int(h5["/input/mom"].attrs.get("restart_valid", 0)) != 0


def output_restart_state_valid(up_file: Path):
    import h5py

    with h5py.File(up_file, "r") as h5:
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            return False
        return int(h5["/output/mom"].attrs.get("restart_final_state_valid", 0)) != 0


def mark_output_restart_state(up_file: Path, nsteps: int, dt: float):
    import h5py

    expected_time = float(nsteps) * float(dt)
    with h5py.File(up_file, "r+") as h5:
        if "/output/time" not in h5 or h5["/output/time"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no output/time; cannot mark restart state")
        if "/output/mom" not in h5 or h5["/output/mom"].shape[0] == 0:
            raise RuntimeError(f"{up_file} has no output/mom; cannot mark restart state")
        final_time = float(h5["/output/time"][-1])
        tolerance = max(1e-6, abs(float(dt)) * 0.51)
        if abs(final_time - expected_time) > tolerance:
            raise RuntimeError(
                f"{up_file} final output time {final_time:.10g} does not match expected {expected_time:.10g}"
            )
        h5["/output/mom"].attrs["restart_final_state_valid"] = np.int8(1)
        h5["/output/mom"].attrs["restart_duration_steps"] = np.int64(nsteps)
        h5["/output/mom"].attrs["restart_time_step"] = float(dt)
        h5["/output/mom"].attrs["restart_final_time"] = final_time


def run_minimization_stage(args, stage_label: str, up_file: Path, max_iter: int):
    if int(max_iter) <= 0:
        print(f"=== Stage {stage_label}: Minimization skipped (max_iter <= 0) ===")
        return
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
        "--record-momentum",
    ]
    if input_momentum_restart_valid(output_file):
        cmd.append("--restart-using-momentum")
    write_stage_debug_outputs(output_file, debug_dir=args.run_dir / "debug")
    run_checked(cmd, log_file=args.log_dir / f"stage_{stage_label}.log")
    mark_output_restart_state(output_file, nsteps, dt)


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


def infer_next_stage70_label_from_source(source_file: Path, pdb_id: str):
    pattern = re.compile(rf"^{re.escape(pdb_id)}\.stage_7\.(\d+)\.up$")
    match = pattern.fullmatch(source_file.name)
    if not match:
        return "7.1"
    return f"7.{int(match.group(1)) + 1}"


def resolve_continuation_outputs(args):
    if args.previous_run_dir or args.previous_stage7_file:
        raise ValueError(
            "Legacy continuation discovery flags are no longer supported. "
            "Use --continue-stage-70-from with an explicit source stage file."
        )
    if args.auto_continue_from_previous_run or args.auto_continue_glob:
        raise ValueError(
            "Auto continuation discovery is no longer supported. "
            "Use --continue-stage-70-from with an explicit source stage file."
        )
    if args.continue_stage_70_from is None:
        return None, None, None
    source = args.continue_stage_70_from.resolve()
    label = args.continue_stage_70_label
    if not label and args.continue_stage_70_output:
        label = infer_stage70_label_from_file(args.continue_stage_70_output, args.pdb_id)
        if not label:
            raise ValueError(f"CONTINUE_STAGE_70_OUTPUT must be named {args.pdb_id}.stage_7.N.up")
    if not label:
        label = infer_next_stage70_label_from_source(source, args.pdb_id)
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
    if not output_restart_state_valid(source_file):
        raise ValueError(
            "Production continuation requires a source stage generated by the current workflow with "
            "a validated final restart state. Regenerate the previous production stage with the current workflow."
        )
    handoff_initial_position(args, source_file, output_file, "production_restart", args.prod_time_step)
    if not input_momentum_restart_valid(output_file):
        raise ValueError(
            "Production continuation requires restart-valid /output/mom in the source stage. "
            "Regenerate the previous production stage with the current workflow so momentum is recorded."
        )
    if args.debug_rigid_protein:
        set_debug_rigid_protein(output_file)
    assert_hybrid_stage_active(output_file, "production", "production")
    run_md_stage(args, stage_label, output_file, output_file, args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
    extract_stage_vtf(args, stage_label, output_file, "2")


def normalize_hybrid_workflow_args(args):
    args.upside_home = workflow_path(args.upside_home, REPO_ROOT).resolve()
    args.run_dir = workflow_path(args.run_dir).resolve()
    args.checkpoint_dir = args.run_dir / "checkpoints"
    args.log_dir = args.run_dir / "logs"
    args.hybrid_prep_dir = args.run_dir / "hybrid_prep"
    args.runtime_pdb_file = args.hybrid_prep_dir / f"{args.runtime_pdb_id}.MARTINI.pdb"
    args.hybrid_mapping_file = args.hybrid_prep_dir / "hybrid_mapping.h5"
    args.hybrid_packed_pdb = args.hybrid_prep_dir / "hybrid_packed.MARTINI.pdb"
    args.upside_executable = args.upside_home / "obj" / "upside"
    args.martini_ff_dir = args.upside_home / "parameters" / "dryMARTINI"
    args.mass_ff_file = args.martini_ff_dir / "dry_martini_v2.1.itp"
    args.martini_table_h5 = args.run_dir / "martini.h5"
    args.sc_martini_library = args.martini_table_h5
    args.upside_rama_library = args.upside_home / "parameters" / "common" / "rama.dat"
    args.upside_rama_sheet_mixing = args.upside_home / "parameters" / "ff_2.1" / "sheet"
    args.upside_hbond_energy = args.upside_home / "parameters" / "ff_2.1" / "hbond.h5"
    args.upside_reference_state_rama = args.upside_home / "parameters" / "common" / "rama_reference.pkl"
    args.universal_prep_mode = "both"
    args.hybrid_validate = True
    args.hybrid_preprod_activation_stage = "minimization"
    args.sc_env_lj_force_cap = DEFAULT_SC_ENV_LJ_FORCE_CAP
    args.sc_env_coul_force_cap = DEFAULT_SC_ENV_COUL_FORCE_CAP
    args.nonprotein_hs_force_cap = DEFAULT_NONPROTEIN_HS_FORCE_CAP
    args.sc_env_po4_z_clamp_enable = DEFAULT_SC_ENV_PO4_Z_CLAMP_ENABLE
    args.sc_env_relax_steps = DEFAULT_SC_ENV_RELAX_STEPS
    args.sc_env_backbone_hold_steps = DEFAULT_SC_ENV_BACKBONE_HOLD_STEPS
    args.sc_env_po4_z_hold_steps = DEFAULT_SC_ENV_PO4_Z_HOLD_STEPS
    args.npt_tau = DEFAULT_NPT_TAU
    args.npt_interval = DEFAULT_NPT_INTERVAL
    args.prod_70_barostat_type = DEFAULT_PROD_70_BAROSTAT_TYPE
    args.martini_energy_conversion = DEFAULT_MARTINI_ENERGY_CONVERSION
    args.martini_length_conversion = DEFAULT_MARTINI_LENGTH_CONVERSION
    args.bar_1_to_eup_per_a3 = DEFAULT_BAR_1_TO_EUP_PER_A3
    args.compressibility_3e4_bar_inv_to_a3_per_eup = DEFAULT_COMPRESSIBILITY_3E4_BAR_INV_TO_A3_PER_EUP
    args.protein_env_interface_scale = DEFAULT_PROTEIN_ENV_INTERFACE_SCALE
    args.debug_rigid_protein = bool(args.debug_rigid_protein)
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


def _decode_stage_strings(values):
    out = []
    for value in np.asarray(values):
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(value).strip())
    return np.asarray(out, dtype=object)


def _detect_has_bonded_environment_particles(up_file: Path) -> bool:
    """Return True when the stage has real bonded dry-MARTINI environment pairs."""
    import h5py

    with h5py.File(up_file, "r") as h5:
        bond_path = "/input/potential/dist_spring/id"
        if bond_path not in h5:
            print("  No dist_spring bonds found in stage input")
            print("  Extended pre-7.0 equilibrium stages needed: False")
            return False

        bonds = np.asarray(h5[bond_path][:], dtype=np.int64)
        if bonds.ndim != 2 or bonds.shape[1] < 2 or bonds.size == 0:
            print("  Empty dist_spring bond list in stage input")
            print("  Extended pre-7.0 equilibrium stages needed: False")
            return False

        inp = h5["/input"]
        atom_types = (
            _decode_stage_strings(inp["type"][:])
            if "type" in inp
            else np.asarray([""] * int(inp["pos"].shape[0]), dtype=object)
        )
        env_mask = None
        membership_path = "/input/hybrid_env_topology/protein_membership"
        if membership_path in h5:
            membership = np.asarray(h5[membership_path][:], dtype=np.int32)
            env_mask = membership < 0
        elif "particle_class" in inp:
            classes = np.char.upper(_decode_stage_strings(inp["particle_class"][:]).astype(str))
            env_mask = (classes != "PROTEIN") & (classes != "PROTEINAA")
        else:
            env_mask = np.ones(atom_types.shape[0], dtype=bool)

        bonded_env_pairs = []
        ignored_orientation_bonds = 0
        for raw_i, raw_j in bonds[:, :2]:
            i = int(raw_i)
            j = int(raw_j)
            if i < 0 or j < 0 or i >= atom_types.shape[0] or j >= atom_types.shape[0]:
                continue
            pair_types = {str(atom_types[i]).strip().upper(), str(atom_types[j]).strip().upper()}
            if pair_types == {"CGL", "CGLD"}:
                ignored_orientation_bonds += 1
                continue
            if bool(env_mask[i]) and bool(env_mask[j]):
                bonded_env_pairs.append((i, j, tuple(sorted(pair_types))))

    print(f"  Ignored synthetic CGL-CGLD orientation bonds: {ignored_orientation_bonds}")
    print(f"  Real bonded dry-MARTINI environment pairs: {len(bonded_env_pairs)}")
    if bonded_env_pairs:
        preview = ", ".join(
            f"{i}:{j}({'/'.join(pair_types)})" for i, j, pair_types in bonded_env_pairs[:5]
        )
        print(f"  Example bonded environment pairs: {preview}")
    print(f"  Extended pre-7.0 equilibrium stages needed: {bool(bonded_env_pairs)}")
    return bool(bonded_env_pairs)


def run_hybrid_workflow_command(argv):
    parser = argparse.ArgumentParser(description="Run the hybrid 1RKL dry-MARTINI workflow.")
    parser.add_argument("--pdb-id", default=env_default("PDB_ID", "1rkl"))
    parser.add_argument("--runtime-pdb-id", default=env_default("RUNTIME_PDB_ID", None))
    parser.add_argument("--upside-home", default=env_default("UPSIDE_HOME", str(REPO_ROOT)))
    parser.add_argument("--run-dir", default=env_default("RUN_DIR", "outputs/martini_test_1rkl_hybrid"))
    parser.add_argument("--protein-aa-pdb", default=env_default("PROTEIN_AA_PDB", None))
    parser.add_argument("--bilayer-pdb", default=env_default("BILAYER_PDB", None))
    parser.add_argument("--extract-vtf-script", default=env_default("EXTRACT_VTF_SCRIPT", str(PY_DIR / "martini_extract_vtf.py")))
    parser.add_argument("--salt-molar", type=float, default=env_float("SALT_MOLAR", 0.15))
    parser.add_argument("--protein-lipid-cutoff", type=float, default=env_float("PROTEIN_LIPID_CUTOFF", 4.5))
    parser.add_argument("--ion-cutoff", type=float, default=env_float("ION_CUTOFF", 10.0))
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
    parser.add_argument("--min-61-max-iter", type=int, default=env_int("MIN_61_MAX_ITER", 0))
    parser.add_argument("--eq-60-nsteps", type=int, default=env_int("EQ_60_NSTEPS", 500))
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
    parser.add_argument("--prep-seed", default=os.environ.get("PREP_SEED"))
    parser.add_argument("--seed", default=os.environ.get("SEED"))
    parser.add_argument("--continue-stage-70-from", default=env_default("CONTINUE_STAGE_70_FROM", ""))
    parser.add_argument("--continue-stage-70-output", default=env_default("CONTINUE_STAGE_70_OUTPUT", ""))
    parser.add_argument("--continue-stage-70-label", default=env_default("CONTINUE_STAGE_70_LABEL", ""))
    parser.add_argument("--previous-run-dir", default=env_default("PREVIOUS_RUN_DIR", ""))
    parser.add_argument("--previous-stage7-file", default=env_default("PREVIOUS_STAGE7_FILE", ""))
    parser.add_argument("--auto-continue-from-previous-run", type=int, default=env_int("AUTO_CONTINUE_FROM_PREVIOUS_RUN", 0))
    parser.add_argument("--auto-continue-glob", default=env_default("AUTO_CONTINUE_GLOB", ""))
    parser.add_argument("--initial-debug-only", type=int, default=env_int("INITIAL_DEBUG_ONLY", 0))
    parser.add_argument("--debug-rigid-protein", type=int, default=env_int("DEBUG_RIGID_PROTEIN", 0))
    args = parser.parse_args(argv)
    if args.runtime_pdb_id is None:
        args.runtime_pdb_id = f"{args.pdb_id}_hybrid"
    if args.protein_aa_pdb is None:
        args.protein_aa_pdb = f"pdb/{args.pdb_id}.pdb"
    if args.bilayer_pdb is None:
        args.bilayer_pdb = str(Path(args.upside_home) / "parameters" / "dryMARTINI" / "DOPC.pdb")
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

    prepare_workflow_hybrid_artifacts(args)
    write_stage70_initial_debug_files(args)

    if int(args.initial_debug_only):
        print("=== Initial Debug Only Complete ===")
        return

    # Extract the canonical dry-MARTINI DOPC bead reference for CG lipid table fitting.
    # The system bilayer PDB still controls placement; the potential uses the force-field
    # reference lipid so table coefficients do not depend on one packed snapshot conformation.
    cg_lipid_config = None
    try:
        from martini_prepare_system_lib import _DOPC_ATOM_NAMES

        table_ref_pdb = (args.upside_home / "parameters" / "dryMARTINI" / "DOPC.pdb").resolve()
        dopc_atoms = []
        with open(table_ref_pdb) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                resname = line[17:21].strip().upper()
                if resname not in ("DOPC", "DOP"):
                    continue
                aname = line[12:16].strip().upper()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                dopc_atoms.append({"name": aname, "x": x, "y": y, "z": z})

        n_per_lipid = 14
        if len(dopc_atoms) >= n_per_lipid:
            first = dopc_atoms[:n_per_lipid]
            com = np.mean([[a["x"], a["y"], a["z"]] for a in first], axis=0)
            ref_bead_positions = np.array(
                [[a["x"] - com[0], a["y"] - com[1], a["z"] - com[2]] for a in first]
            )
            ref_bead_positions_nm = ref_bead_positions * 0.1

            from martini_build_tables import _DOPC_BEAD_TYPES

            cg_lipid_config = {
                "ref_bead_positions_nm": ref_bead_positions_nm,
                "bead_types": list(_DOPC_BEAD_TYPES),
            }
            print(f"Extracted CG lipid table reference from {table_ref_pdb.name}: "
                  f"{len(dopc_atoms)} DOPC atoms ({len(dopc_atoms) // n_per_lipid} lipids)")
        else:
            print(f"Fewer than {n_per_lipid} DOPC atoms in {table_ref_pdb}, "
                  f"skipping CG lipid table build")
    except Exception as e:
        print(f"Could not extract CG lipid reference positions: {e}")

    skip_build = os.environ.get("UPSIDE_SKIP_MARTINI_BUILD", "")
    if skip_build and args.martini_table_h5.exists():
        print(f"UPSIDE_SKIP_MARTINI_BUILD set, reusing existing {args.martini_table_h5}")
    else:
        build_martini_tables(
            output_path=args.martini_table_h5,
            dry_ff_path=args.mass_ff_file,
            martinize_path=args.upside_home / "py" / "martinize.py",
            sidechain_lib_path=args.upside_home / "parameters" / "ff_2.1" / "sidechain.h5",
            forcefield_name="martini22",
            cg_lipid_config=cg_lipid_config,
        )

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

    print("=== Stage 6.0: rigid-protein NPT box relaxation ===")
    prepare_stage_file(args, files["prepared_60"], "npt_equil", 1, 0, 0, "minimization")
    shutil.copy2(files["prepared_60"], files["stage_60"])
    run_md_stage(
        args,
        "6.0",
        files["stage_60"],
        files["stage_60"],
        args.eq_60_nsteps,
        args.eq_time_step,
        args.eq_frame_steps,
    )
    extract_stage_vtf(args, "6.0", files["stage_60"], "1")
    needs_pre70 = _detect_has_bonded_environment_particles(files["stage_60"])

    if needs_pre70:
        print("=== Bonded dry-MARTINI environment detected -> running extended pre-7.0 equilibrium ===")

        prepare_stage_file(args, files["prepared_61"], "npt_prod", 1, 0, 0, "minimization")
        shutil.copy2(files["prepared_61"], files["stage_61"])
        handoff_initial_position(args, files["stage_60"], files["stage_61"])
        run_minimization_stage(args, "6.1", files["stage_61"], args.min_61_max_iter)
        extract_stage_vtf(args, "6.1", files["stage_61"], "1")

        prepare_stage_file(args, files["stage_62"], "npt_equil", 1, 0, 200, "minimization")
        handoff_initial_position(args, files["stage_61"], files["stage_62"])
        run_md_stage(args, "6.2", files["stage_62"], files["stage_62"], args.eq_62_nsteps, args.eq_time_step, args.eq_frame_steps)
        extract_stage_vtf(args, "6.2", files["stage_62"], "1")

        prepare_stage_file(args, files["prepared_63"], "npt_equil_reduced", 1, 0, 100, "minimization")
        shutil.copy2(files["prepared_63"], files["stage_63"])
        handoff_initial_position(args, files["stage_62"], files["stage_63"])
        run_md_stage(args, "6.3", files["stage_63"], files["stage_63"], args.eq_63_nsteps, args.eq_time_step, args.eq_frame_steps)
        extract_stage_vtf(args, "6.3", files["stage_63"], "1")

        prepare_stage_file(args, files["prepared_64"], "npt_prod", 1, 0, 50, "minimization")
        shutil.copy2(files["prepared_64"], files["stage_64"])
        handoff_initial_position(args, files["stage_63"], files["stage_64"])
        run_md_stage(args, "6.4", files["stage_64"], files["stage_64"], args.eq_64_nsteps, args.eq_time_step, args.eq_frame_steps)
        extract_stage_vtf(args, "6.4", files["stage_64"], "1")

        prepare_stage_file(args, files["prepared_65"], "npt_prod", 1, 0, 20, "minimization")
        shutil.copy2(files["prepared_65"], files["stage_65"])
        handoff_initial_position(args, files["stage_64"], files["stage_65"])
        run_md_stage(args, "6.5", files["stage_65"], files["stage_65"], args.eq_65_nsteps, args.eq_time_step, args.eq_frame_steps)
        extract_stage_vtf(args, "6.5", files["stage_65"], "1")

        prepare_stage_file(args, files["prepared_66"], "npt_prod", 1, 0, 10, "minimization")
        shutil.copy2(files["prepared_66"], files["stage_66"])
        handoff_initial_position(args, files["stage_65"], files["stage_66"])
        run_md_stage(args, "6.6", files["stage_66"], files["stage_66"], args.eq_66_nsteps, args.eq_time_step, args.eq_frame_steps)
        extract_stage_vtf(args, "6.6", files["stage_66"], "1")

        prepare_stage_file(args, files["prepared_70"], "npt_prod", args.prod_70_npt_enable, args.prod_70_barostat_type, 0, "production")
        shutil.copy2(files["prepared_70"], files["stage_70"])
        handoff_initial_position(args, files["stage_66"], files["stage_70"], "production_hybrid")
        if args.debug_rigid_protein:
            set_debug_rigid_protein(files["stage_70"])
        assert_hybrid_stage_active(files["stage_70"], "production", "production")
        run_md_stage(args, "7.0", files["stage_70"], files["stage_70"], args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
        extract_stage_vtf(args, "7.0", files["stage_70"], "2")
    else:
        print("=== No bonded dry-MARTINI environment pairs -> handoff from stage 6.0 to stage 7.0 ===")
        prepare_stage_file(args, files["prepared_70"], "npt_prod", args.prod_70_npt_enable, args.prod_70_barostat_type, 0, "production")
        shutil.copy2(files["prepared_70"], files["stage_70"])
        handoff_initial_position(args, files["stage_60"], files["stage_70"], "production_hybrid")
        if args.debug_rigid_protein:
            set_debug_rigid_protein(files["stage_70"])
        assert_hybrid_stage_active(files["stage_70"], "production", "production")
        run_md_stage(args, "7.0", files["stage_70"], files["stage_70"], args.prod_70_nsteps, args.prod_time_step, args.prod_frame_steps)
        extract_stage_vtf(args, "7.0", files["stage_70"], "2")

    print("=== Workflow Complete ===")


if __name__ == "__main__":
    argv = sys.argv[1:]
    if not argv or argv[0] != "run-hybrid-workflow":
        raise SystemExit(
            "Unsupported command. Use:\n"
            "  martini_prepare_system.py run-hybrid-workflow [options]"
        )
    run_hybrid_workflow_command(argv[1:])
