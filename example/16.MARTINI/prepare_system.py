#!/usr/bin/env python3

import argparse
import json
import os
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

from prepare_system_lib import (
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
    parse_pdb,
    place_ions,
    remove_overlapping_lipids,
    set_box_from_lipid_xy,
    tile_and_crop_bilayer_lipids,
    validate_backbone_reference_frame,
    write_hybrid_mapping_h5,
    write_pdb,
)


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args():
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
    return parser.parse_args()


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
    shift = bilayer_center - protein_center
    protein_xyz = protein_xyz + shift
    for atom, c in zip(protein_atoms, protein_xyz):
        atom["x"], atom["y"], atom["z"] = float(c[0]), float(c[1]), float(c[2])
    if uses_martini_protein:
        protein_aa_xyz = coords(protein_aa_atoms)
        protein_aa_xyz = protein_aa_xyz + shift
        for atom, c in zip(protein_aa_atoms, protein_aa_xyz):
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


def main():
    args = parse_args()
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


if __name__ == "__main__":
    main()
