#!/usr/bin/env python3

import argparse
import importlib.util
import os
from pathlib import Path

import h5py
import numpy as np
import tables as tb
import _pickle as cPickle


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


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Inject Upside backbone-only nodes into a stage .up file using runtime "
            "N/CA/C carriers. This removes rotamer/placement sidechain nodes."
        )
    )
    p.add_argument("up_file")
    p.add_argument("protein_source")
    p.add_argument("upside_home")
    p.add_argument("rama_library")
    p.add_argument("rama_sheet_mixing")
    p.add_argument("hbond_energy")
    p.add_argument("reference_state_rama")
    return p.parse_args()


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
    name = name.strip().upper()
    return aliases.get(name, name)


def parse_itp_residue_names(path):
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
            if role != "BB":
                continue
            if resnr in seen:
                continue
            seen.add(resnr)
            resnames.append(normalize_resname(parts[3]))
    return resnames


def parse_pdb_residue_names(path):
    aa_res = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX",
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


def parse_residue_names(path):
    suffix = Path(path).suffix.lower()
    if suffix == ".pdb":
        return parse_pdb_residue_names(path)
    return parse_itp_residue_names(path)


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


def require_file(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    return p


def load_upside_config(upside_home):
    upside_config_py = Path(upside_home) / "py" / "upside_config.py"
    require_file(upside_config_py)
    spec = importlib.util.spec_from_file_location("upside_config_runtime", str(upside_config_py))
    uc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(uc)
    return uc


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
            rid = int(resid)
            if rid not in seen:
                residue_ids.append(rid)
                seen.add(rid)
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
            if rid not in residue_to_ncac:
                residue_to_ncac[rid] = (n_rt, ca_rt, c_rt)

    return residue_ids, residue_to_ncac


def write_backbone_nodes(
    uc,
    up_file,
    residue_ids,
    residue_to_ncac,
    itp_resnames,
    rama_library,
    rama_sheet_mixing,
    hbond_energy,
    reference_state_rama,
):
    if len(itp_resnames) != len(residue_ids):
        raise ValueError(
            f"ITP/BB residue count mismatch: ITP has {len(itp_resnames)} residues, bb map has {len(residue_ids)}"
        )

    ref_n_atom = 3 * len(residue_ids)
    atom_map = np.full((ref_n_atom,), -1, dtype=np.int64)
    for i, rid in enumerate(residue_ids):
        if rid not in residue_to_ncac:
            raise ValueError(f"Missing N/CA/C runtime mapping for residue {rid}")
        n_idx, ca_idx, c_idx = residue_to_ncac[rid]
        atom_map[3 * i + 0] = n_idx
        atom_map[3 * i + 1] = ca_idx
        atom_map[3 * i + 2] = c_idx

    with h5py.File(up_file, "r+") as up:
        inp = up["/input"]
        pot = inp["potential"]

        # Remove sidechain/placement nodes and any previous backbone node set.
        for node_name in SIDECHAIN_NODES + BACKBONE_NODES:
            if node_name in pot:
                del pot[node_name]

        if "sequence" in inp:
            del inp["sequence"]
        seq_data = np.asarray([np.bytes_(x) for x in itp_resnames], dtype="S3")
        inp.create_dataset("sequence", data=seq_data)

    fasta_seq = np.array(itp_resnames)
    spring_args = type(
        "SpringArgs",
        (),
        {"bond_stiffness": 48.0, "angle_stiffness": 175.0, "omega_stiffness": 30.0},
    )()

    with tb.open_file(up_file, mode="a") as tf:
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
            rama_library,
            rama_sheet_mixing,
            secstr_bias="",
            mode="mixture",
            param_deriv=False,
        )
        uc.write_affine_alignment(len(fasta_seq))

        ref_state_cor = np.log(cPickle.load(open(reference_state_rama, "rb"), encoding="latin1"))
        ref_state_cor -= ref_state_cor.mean()
        grp = tf.create_group(tf.root.input.potential, "rama_map_pot_ref")
        grp._v_attrs.arguments = np.array([b"rama_coord"])
        grp._v_attrs.log_pot = 0
        uc.create_array(grp, "residue_id", obj=np.arange(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_map_id", obj=np.zeros(len(fasta_seq), dtype=np.int32))
        uc.create_array(grp, "rama_pot", obj=ref_state_cor[None])

        uc.write_infer_H_O(fasta_seq, np.array([], dtype=np.int32))
        uc.write_count_hbond(fasta_seq, False)
        uc.write_short_hbond(fasta_seq, hbond_energy)
        uc.write_rama_coord2()
        uc.write_backbone_pair(fasta_seq)

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

        for grp_name, dset_name in remap_datasets:
            ds = pot[grp_name][dset_name]
            ds[...] = remap_atom_index_array(ds[:], atom_map).astype(ds.dtype, copy=False)

    return ref_n_atom


def main():
    args = parse_args()
    require_file(args.up_file)
    require_file(args.protein_source)
    require_file(args.rama_library)
    require_file(args.rama_sheet_mixing)
    require_file(args.hbond_energy)
    require_file(args.reference_state_rama)

    uc = load_upside_config(args.upside_home)
    residue_ids, residue_to_ncac = collect_runtime_ncac_mapping(args.up_file)
    itp_resnames = parse_residue_names(args.protein_source)
    ref_n_atom = write_backbone_nodes(
        uc=uc,
        up_file=args.up_file,
        residue_ids=residue_ids,
        residue_to_ncac=residue_to_ncac,
        itp_resnames=itp_resnames,
        rama_library=args.rama_library,
        rama_sheet_mixing=args.rama_sheet_mixing,
        hbond_energy=args.hbond_energy,
        reference_state_rama=args.reference_state_rama,
    )
    print(
        f"Injected backbone-only Upside nodes into {args.up_file}: "
        f"n_res={len(residue_ids)} ref_backbone_atoms={ref_n_atom}"
    )


if __name__ == "__main__":
    main()
