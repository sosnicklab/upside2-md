#!/usr/bin/env python3

import argparse
import importlib.util
import types
import warnings
from pathlib import Path
import _pickle as cPickle

import h5py
import numpy as np
import tables as tb


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
LEGACY_STAGE7_NODES = [
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
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


def parse_args():
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
    return parser.parse_args()


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
            if low == "[ atoms ]" or low == "[atoms]":
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


def build_selected_residue_rows(sequence, restype_to_index):
    cb_index = []
    residue_table_index = []
    skipped = []
    for seq_idx, resname in enumerate(sequence):
        norm = normalize_resname(resname)
        table_idx = restype_to_index.get(norm)
        if table_idx is None:
            skipped.append((seq_idx, norm))
            continue
        cb_index.append(seq_idx)
        residue_table_index.append(table_idx)
    if not cb_index:
        raise ValueError("No stage-7 residues matched martini.h5 restype_order")
    return (
        np.asarray(cb_index, dtype=np.int32),
        np.asarray(residue_table_index, dtype=np.int32),
        skipped,
    )


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

    return (
        np.asarray(env_atom_index, dtype=np.int32),
        np.asarray(env_target_index, dtype=np.int32),
    )


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
    spring_args = types.SimpleNamespace(
        bond_stiffness=48.0,
        angle_stiffness=175.0,
        omega_stiffness=30.0,
    )

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


def main():
    args = parse_args()
    up_file = Path(args.up_file).expanduser().resolve()
    martini_h5 = Path(args.martini_h5).expanduser().resolve()
    upside_home = Path(args.upside_home).expanduser().resolve()
    rama_library = Path(args.rama_library).expanduser().resolve()
    rama_sheet_mixing = Path(args.rama_sheet_mixing).expanduser().resolve()
    hbond_energy = Path(args.hbond_energy).expanduser().resolve()
    reference_state_rama = Path(args.reference_state_rama).expanduser().resolve()
    protein_itp = Path(args.protein_itp).expanduser().resolve() if args.protein_itp else None

    if not up_file.exists():
        raise SystemExit(f"ERROR: stage file not found: {up_file}")
    if not martini_h5.exists():
        raise SystemExit(f"ERROR: martini.h5 not found: {martini_h5}")
    for path in [rama_library, rama_sheet_mixing, hbond_energy, reference_state_rama]:
        if not path.exists():
            raise SystemExit(f"ERROR: required Upside input not found: {path}")

    with h5py.File(martini_h5, "r") as sc_lib:
        restype_order = decode_string_array(sc_lib["restype_order"])
        target_order = decode_string_array(sc_lib["target_order"])
        grid_nm = sc_lib["grid_nm"][:].astype(np.float32)
        energy_kj_mol = sc_lib["energy_kj_mol"][:].astype(np.float32)

    restype_to_index = {name: i for i, name in enumerate(restype_order)}
    target_to_index = {name: i for i, name in enumerate(target_order)}

    with h5py.File(up_file, "r+") as up:
        inp = up["input"]
        pot = inp["potential"]
        martini_potential = pot["martini_potential"]

        residue_ids, affine_atoms = build_affine_atoms(inp)
        sequence = resolve_sequence(inp, len(residue_ids), protein_itp)

        cb_index, residue_table_index, skipped = build_selected_residue_rows(sequence, restype_to_index)
        env_atom_index, env_target_index = build_env_rows(inp, target_to_index)

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

        g_cb = recreate_group(pot, "placement_fixed_point_only_CB")
        g_cb.attrs["arguments"] = np.asarray([np.bytes_("affine_alignment")])
        g_cb.create_dataset("affine_residue", data=np.arange(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset("layer_index", data=np.zeros(len(sequence), dtype=np.int32), dtype=np.int32)
        g_cb.create_dataset("placement_data", data=CB_PLACEMENT, dtype=np.float32)

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
        g_sc = recreate_group(pot, "martini_sc_table_potential")
        g_sc.attrs["arguments"] = np.asarray([np.bytes_("pos"), np.bytes_("placement_fixed_point_only_CB")])
        g_sc.attrs["energy_conversion_kj_per_eup"] = np.float32(
            martini_potential.attrs["energy_conversion_kj_per_eup"]
        )
        g_sc.attrs["length_conversion_angstrom_per_nm"] = np.float32(
            martini_potential.attrs["length_conversion_angstrom_per_nm"]
        )
        g_sc.attrs["x_len"] = np.float32(martini_potential.attrs["x_len"])
        g_sc.attrs["y_len"] = np.float32(martini_potential.attrs["y_len"])
        g_sc.attrs["z_len"] = np.float32(martini_potential.attrs["z_len"])
        g_sc.create_dataset("cb_index", data=cb_index, dtype=np.int32)
        g_sc.create_dataset("residue_table_index", data=residue_table_index, dtype=np.int32)
        g_sc.create_dataset("env_atom_index", data=env_atom_index, dtype=np.int32)
        g_sc.create_dataset("env_target_index", data=env_target_index, dtype=np.int32)
        g_sc.create_dataset("grid_nm", data=grid_nm, dtype=np.float32)
        g_sc.create_dataset("energy_kj_mol", data=energy_kj_mol, dtype=np.float32)
        g_sc.create_dataset(
            "restype_order",
            data=np.asarray([np.bytes_(x) for x in restype_order], dtype="S4"),
        )
        g_sc.create_dataset(
            "target_order",
            data=np.asarray([np.bytes_(x) for x in target_order], dtype="S8"),
        )

        print(
            f"Injected martini_sc_table_potential into {up_file}: "
            f"n_residues={len(cb_index)} skipped={len(skipped)} n_env={len(env_atom_index)} "
            f"n_restypes={len(restype_order)} n_targets={len(target_order)}"
        )


if __name__ == "__main__":
    main()
