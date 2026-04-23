#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np

from martini_prepare_system_lib import (
    DEFAULT_WORKFLOW_BILAYER_PDB,
    DEFAULT_SC_TABLE_JSON,
    build_sc_martini_h5,
    build_sidechain_centroid_proxy_atoms,
    center_of_mass,
    collect_aa_backbone_map,
    compute_lipid_residue_indices,
    convert_stage,
    coords,
    estimate_salt_pairs,
    extract_backbone_sequence,
    extract_protein_aa_atoms,
    extract_protein_aa_backbone_atoms,
    infer_effective_ion_volume_fraction_from_template,
    infer_protein_charge_from_aa,
    inject_stage7_sc_table_nodes,
    lipid_resname,
    parse_pdb,
    place_ions,
    remove_overlapping_lipids,
    set_box_from_lipid_xy,
    set_coords,
    set_initial_position,
    tile_and_crop_bilayer_lipids,
    validate_hybrid_mapping,
    write_backbone_metadata_h5,
    write_pdb,
)


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)


def abs_path(path):
    return os.path.abspath(os.path.expanduser(path))


def runtime_pdb_path(args):
    if args.runtime_pdb_output:
        return abs_path(args.runtime_pdb_output)
    return os.path.join(REPO_ROOT, "example", "16.MARTINI", "pdb", "%s.MARTINI.pdb" % args.pdb_id)


def parse_prepare_args(argv):
    parser = argparse.ArgumentParser(
        description="Prepare the default 1RKL mixed AA-backbone + dry-MARTINI workflow input."
    )
    parser.add_argument("--pdb-id", required=True, help="Runtime PDB id for stage conversion")
    parser.add_argument("--runtime-pdb-output", default=None)
    parser.add_argument("--prepare-structure", type=int, default=1, choices=[0, 1])
    parser.add_argument("--stage", default=None, help="Stage name for UPSIDE input conversion")
    parser.add_argument("--run-dir", default="outputs/martini_test")
    parser.add_argument(
        "--bilayer-pdb",
        default=DEFAULT_WORKFLOW_BILAYER_PDB,
    )
    parser.add_argument("--protein-aa-pdb", required=True)
    parser.add_argument("--hybrid-mapping-output", default=None)
    parser.add_argument("--xy-scale", type=float, default=1.0)
    parser.add_argument("--box-padding-xy", type=float, default=0.0)
    parser.add_argument("--box-padding-z", type=float, default=20.0)
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--protein-lipid-cutoff", type=float, default=3.0)
    parser.add_argument("--protein-net-charge", type=int, default=None)
    parser.add_argument("--stage-lipidhead-fc", type=float, default=0.0)
    parser.add_argument("--npt-enable", type=int, default=0)
    parser.add_argument("--npt-target-pxy", type=float, default=0.0)
    parser.add_argument("--npt-target-pz", type=float, default=0.0)
    return parser.parse_args(argv)


def prepare_mixed_structure(args, runtime_pdb):
    protein_aa_pdb = abs_path(args.protein_aa_pdb)
    bilayer_pdb = abs_path(args.bilayer_pdb)

    if not os.path.exists(protein_aa_pdb):
        raise IOError("Protein AA PDB not found: %s" % protein_aa_pdb)
    if not os.path.exists(bilayer_pdb):
        raise IOError("Bilayer PDB not found: %s" % bilayer_pdb)
    if args.xy_scale < 1.0:
        raise ValueError("--xy-scale must be >= 1.0")

    protein_aa_atoms_raw, _ = parse_pdb(protein_aa_pdb)
    bilayer_atoms, bilayer_box = parse_pdb(bilayer_pdb)
    protein_aa_atoms_full = extract_protein_aa_atoms(protein_aa_atoms_raw)
    protein_atoms = extract_protein_aa_backbone_atoms(protein_aa_atoms_full)
    protein_env_atoms = [atom.copy() for atom in protein_atoms]
    protein_env_atoms.extend(build_sidechain_centroid_proxy_atoms(protein_aa_atoms_full))
    bilayer_lipid_atoms = [atom for atom in bilayer_atoms if lipid_resname(atom["resname"])]

    if not protein_atoms:
        raise ValueError("No AA backbone atoms found in protein AA PDB")
    if not bilayer_lipid_atoms:
        raise ValueError("No lipid residues found in bilayer template")

    protein_xyz = coords(protein_atoms)
    protein_env_xyz = coords(protein_env_atoms)
    bilayer_xyz = coords(bilayer_lipid_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    bilayer_span = bilayer_xyz.max(axis=0) - bilayer_xyz.min(axis=0)
    base_side = max(float(bilayer_span[0]), float(bilayer_span[1]))

    shift = bilayer_center - center_of_mass(protein_xyz)
    protein_xyz = protein_xyz + shift
    protein_env_xyz = protein_env_xyz + shift
    set_coords(protein_atoms, protein_xyz)
    set_coords(protein_env_atoms, protein_env_xyz)

    env_span = protein_env_xyz.max(axis=0) - protein_env_xyz.min(axis=0)
    env_center_xy = 0.5 * (protein_env_xyz.min(axis=0)[:2] + protein_env_xyz.max(axis=0)[:2])
    min_required_xy = 3.0 * env_span[:2] + 2.0 * float(args.box_padding_xy)
    target_side = max(base_side * float(args.xy_scale), float(np.max(min_required_xy)))
    target_xy_min = np.array(
        [env_center_xy[0] - 0.5 * target_side, env_center_xy[1] - 0.5 * target_side],
        dtype=float,
    )
    target_xy_max = np.array(
        [env_center_xy[0] + 0.5 * target_side, env_center_xy[1] + 0.5 * target_side],
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
        protein_atoms=protein_env_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=float(args.protein_lipid_cutoff),
    )
    _ = removed_lipids

    packed_atoms = protein_atoms + bilayer_kept
    box_lengths = set_box_from_lipid_xy(
        all_atoms=packed_atoms,
        lipid_atoms=bilayer_kept,
        pad_z=float(args.box_padding_z),
        force_square_xy=True,
        min_box_z=float(3.0 * env_span[2]),
        center_lipid_in_z=True,
    )

    effective_vol_frac = infer_effective_ion_volume_fraction_from_template(
        bilayer_atoms=bilayer_atoms,
        bilayer_box=bilayer_box,
        salt_molar=float(args.salt_molar),
    )
    if args.protein_net_charge is not None:
        protein_charge = int(args.protein_net_charge)
    else:
        protein_charge = int(infer_protein_charge_from_aa(protein_aa_atoms_full))

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

    if args.hybrid_mapping_output:
        mapping_h5 = abs_path(args.hybrid_mapping_output)
        mapping_dir = os.path.dirname(mapping_h5)
        if mapping_dir:
            os.makedirs(mapping_dir, exist_ok=True)
        bb_entries = collect_aa_backbone_map(protein_atoms)
        write_backbone_metadata_h5(
            path=mapping_h5,
            bb_entries=bb_entries,
            total_atoms=len(all_atoms),
            env_atom_indices=list(range(len(protein_atoms), len(all_atoms))),
            n_protein_atoms=len(protein_atoms),
            sequence=extract_backbone_sequence(protein_atoms),
        )
        validate_hybrid_mapping(mapping_h5, n_atom=len(all_atoms))


def run_stage_conversion(args, runtime_pdb):
    convert_stage(
        pdb_id=args.pdb_id,
        stage=args.stage,
        run_dir=args.run_dir,
        runtime_pdb_file=runtime_pdb,
        stage_lipidhead_fc=args.stage_lipidhead_fc,
        npt_enable=args.npt_enable,
        npt_target_pxy=args.npt_target_pxy,
        npt_target_pz=args.npt_target_pz,
    )


def run_prepare(argv):
    args = parse_prepare_args(argv)
    runtime_pdb = runtime_pdb_path(args)
    runtime_dir = os.path.dirname(runtime_pdb)
    if runtime_dir:
        os.makedirs(runtime_dir, exist_ok=True)

    if args.prepare_structure:
        prepare_mixed_structure(args, runtime_pdb)
    elif not os.path.exists(runtime_pdb):
        raise IOError("Runtime PDB not found for stage conversion: %s" % runtime_pdb)

    if args.stage:
        run_stage_conversion(args, runtime_pdb)


def run_build_sc_martini_h5(argv):
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
    build_sc_martini_h5(args.sc_table_json, args.output_h5)


def run_set_initial_position(argv):
    parser = argparse.ArgumentParser(
        description="Copy final coordinates from one stage file into the next stage input."
    )
    parser.add_argument("input_file", help="Source stage file.")
    parser.add_argument("output_file", help="Target stage file to update.")
    args = parser.parse_args(argv)
    set_initial_position(args.input_file, args.output_file)


def run_inject_stage7_sc(argv):
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
    args = parser.parse_args(argv)
    inject_stage7_sc_table_nodes(
        up_file=args.up_file,
        martini_h5=args.martini_h5,
        upside_home=args.upside_home,
        rama_library=args.rama_library,
        rama_sheet_mixing=args.rama_sheet_mixing,
        hbond_energy=args.hbond_energy,
        reference_state_rama=args.reference_state_rama,
    )


def main():
    argv = sys.argv[1:]
    if argv and argv[0] == "build-sc-martini-h5":
        run_build_sc_martini_h5(argv[1:])
    elif argv and argv[0] == "inject-stage7-sc":
        run_inject_stage7_sc(argv[1:])
    elif argv and argv[0] == "set-initial-position":
        run_set_initial_position(argv[1:])
    else:
        run_prepare(argv)


if __name__ == "__main__":
    main()
