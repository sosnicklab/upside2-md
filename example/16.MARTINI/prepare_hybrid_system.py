#!/usr/bin/env python3

import argparse
import json
import shlex
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

NA_AVOGADRO = 6.02214076e23


@dataclass
class Config:
    protein_pdb: Path
    bilayer_pdb: Path
    output_dir: Path
    protein_cg_pdb: Optional[Path]
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
        "--martinize-cmd",
        default=None,
        help=(
            "Optional martinize command template used when --protein-cg-pdb is not given. "
            "Use placeholders {input} and {output}."
        ),
    )
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--protein-lipid-cutoff", type=float, default=5.0)
    parser.add_argument("--ion-cutoff", type=float, default=4.0)
    parser.add_argument("--box-padding-xy", type=float, default=15.0)
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
                "resname": line[17:20].strip(),
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
            chain = atom["chain"] if atom["chain"] else "A"
            segid = atom["segid"][:4].ljust(4) if atom["segid"] else "    "
            f.write(
                f"{atom['record']:<6}{idx:5d} {atom['name'][:4]:>4} {atom['resname'][:3]:>3} {chain:1}"
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
        if atom["resname"].upper() == "DOPC":
            key = (atom["chain"], atom["resseq"], atom["icode"])
            lipid_residues[key].append(idx)
        else:
            keep_nonlipid.append(idx)
    return lipid_residues, keep_nonlipid


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


def estimate_salt_pairs(box_lengths, salt_molar):
    volume_a3 = float(box_lengths[0] * box_lengths[1] * box_lengths[2])
    volume_l = volume_a3 * 1e-27
    pairs = int(round(salt_molar * NA_AVOGADRO * volume_l))
    return max(0, pairs)


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


def collect_bb_map(protein_aa_atoms, protein_cg_atoms):
    aa_by_res = defaultdict(dict)
    for idx, atom in enumerate(protein_aa_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"])
        aa_by_res[key][atom["name"].strip().upper()] = idx

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

        atom_map = aa_by_res[found_key]
        names = ["N", "CA", "C", "O"]
        idxs = [atom_map.get(n, -1) for n in names]
        mask = [1 if i >= 0 else 0 for i in idxs]
        n_present = sum(mask)
        if n_present == 0:
            continue
        weights = [1.0 / n_present if m else 0.0 for m in mask]
        bb_entries.append(
            {
                "bb_resseq": atom["resseq"],
                "bb_chain": atom["chain"],
                "bb_icode": atom["icode"],
                "bb_atom_index": cg_idx,
                "atom_indices": idxs,
                "atom_mask": mask,
                "weights": weights,
            }
        )
    return bb_entries


def collect_sc_map(protein_aa_atoms, protein_cg_atoms):
    aa_by_res = defaultdict(dict)
    for idx, atom in enumerate(protein_aa_atoms):
        key = (atom["chain"], atom["resseq"], atom["icode"])
        aa_by_res[key][atom["name"].strip().upper()] = idx

    sc_entries = []
    for cg_idx, atom in enumerate(protein_cg_atoms):
        atom_name = atom["name"].upper()
        if atom_name == "BB":
            continue
        key = (atom["chain"], atom["resseq"], atom["icode"])
        if key not in aa_by_res:
            cands = [k for k in aa_by_res.keys() if k[1] == atom["resseq"]]
            if len(cands) == 1:
                key = cands[0]
            else:
                continue
        atom_map = aa_by_res[key]
        target_names = ["N", "CA", "C", "O"]
        target_indices = [atom_map.get(n, -1) for n in target_names]
        valid = [i for i in target_indices if i >= 0]
        if not valid:
            continue
        n_valid = len(valid)
        target_weights = [1.0 / n_valid if idx >= 0 else 0.0 for idx in target_indices]
        sc_entries.append(
            {
                "proxy_atom_index": cg_idx,
                "residue_index": atom["resseq"],
                "rotamer_id": 0,
                "rotamer_prob": 1.0,
                "proxy_type": 0,
                "local_pos": [0.0, 0.0, 0.0],
                "proj_target_indices": target_indices,
                "proj_weights": target_weights,
                "protein_id": 0,
            }
        )
    return sc_entries


def write_hybrid_mapping_h5(
    path: Path,
    bb_entries,
    sc_entries,
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
        ctrl.attrs["coupling_align_enable"] = np.int8(0)
        ctrl.attrs["coupling_align_debug"] = np.int8(0)
        ctrl.attrs["coupling_align_interval"] = np.int32(100)
        ctrl.attrs["schema_version"] = np.int32(1)

        bb_grp = inp.create_group("hybrid_bb_map")
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

        sc_grp = inp.create_group("hybrid_sc_map")
        n_sc = len(sc_entries)
        sc_grp.create_dataset(
            "residue_index",
            data=np.array([x["residue_index"] for x in sc_entries], dtype=np.int32),
        )
        sc_grp.create_dataset(
            "rotamer_offset",
            data=np.arange(n_sc + 1, dtype=np.int32),
        )
        sc_grp.create_dataset(
            "rotamer_id",
            data=np.array([x["rotamer_id"] for x in sc_entries], dtype=np.int32),
        )
        sc_grp.create_dataset(
            "proxy_type",
            data=np.array([x["proxy_type"] for x in sc_entries], dtype=np.int32),
        )
        sc_grp.create_dataset(
            "proxy_atom_index",
            data=np.array([x["proxy_atom_index"] for x in sc_entries], dtype=np.int32),
        )
        sc_grp.create_dataset(
            "rotamer_probability",
            data=np.array([x["rotamer_prob"] for x in sc_entries], dtype=np.float32),
        )
        sc_grp.create_dataset(
            "local_pos",
            data=np.array([x["local_pos"] for x in sc_entries], dtype=np.float32).reshape(n_sc, 3),
        )
        sc_grp.create_dataset(
            "proj_target_indices",
            data=np.array([x["proj_target_indices"] for x in sc_entries], dtype=np.int32).reshape(n_sc, 4),
        )
        sc_grp.create_dataset(
            "proj_weights",
            data=np.array([x["proj_weights"] for x in sc_entries], dtype=np.float32).reshape(n_sc, 4),
        )
        sc_grp.create_dataset(
            "protein_id",
            data=np.array([x["protein_id"] for x in sc_entries], dtype=np.int32),
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

    bilayer_xyz = coords(bilayer_atoms)
    protein_xyz = coords(protein_cg_atoms)
    bilayer_center = center_of_mass(bilayer_xyz)
    protein_center = center_of_mass(protein_xyz)
    translated = protein_xyz + (bilayer_center - protein_center)
    set_coords(protein_cg_atoms, translated)

    lipid_residues, keep_nonlipid = compute_lipid_residue_indices(bilayer_atoms)
    bilayer_kept, removed_lipids = remove_overlapping_lipids(
        bilayer_atoms=bilayer_atoms,
        protein_atoms=protein_cg_atoms,
        lipid_residues=lipid_residues,
        keep_nonlipid=keep_nonlipid,
        cutoff=config.protein_lipid_cutoff,
    )

    packed_atoms = protein_cg_atoms + bilayer_kept
    box_lengths = resize_and_shift(
        packed_atoms,
        pad_xy=config.box_padding_xy,
        pad_z=config.box_padding_z,
    )
    if bilayer_box is not None:
        box_lengths = np.maximum(box_lengths, np.array(bilayer_box, dtype=float))

    protein_charge = (
        config.protein_net_charge
        if config.protein_net_charge is not None
        else infer_protein_charge_from_cg(protein_cg_atoms)
    )
    salt_pairs = estimate_salt_pairs(box_lengths, config.salt_molar)
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
    set_coords(all_atoms, np.mod(coords(all_atoms), box_lengths))

    packed_pdb = config.output_dir / "hybrid_packed.MARTINI.pdb"
    write_pdb(packed_pdb, all_atoms, box_lengths)

    bb_entries = collect_bb_map(protein_aa_atoms, protein_cg_atoms)
    sc_entries = collect_sc_map(protein_aa_atoms, protein_cg_atoms)
    mapping_h5 = config.output_dir / "hybrid_mapping.h5"
    env_atom_indices = list(range(len(protein_cg_atoms), len(all_atoms)))
    write_hybrid_mapping_h5(
        mapping_h5,
        bb_entries=bb_entries,
        sc_entries=sc_entries,
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
        "protein_charge_used": int(protein_charge),
        "salt_molar": float(config.salt_molar),
        "salt_pairs_target": int(salt_pairs),
        "na_added": int(n_na),
        "cl_added": int(n_cl),
        "protein_atoms_cg": int(len(protein_cg_atoms)),
        "bilayer_atoms_kept": int(len(bilayer_kept)),
        "lipid_residues_removed": int(removed_lipids),
        "ion_atoms_added": int(len(ion_atoms)),
        "total_atoms": int(len(all_atoms)),
        "bb_map_entries": int(len(bb_entries)),
        "sc_map_entries": int(len(sc_entries)),
    }
    write_summary(config.output_dir / "hybrid_prep_summary.json", summary)

    print(f"Packed system written to: {packed_pdb}")
    print(f"Hybrid mapping HDF5 written to: {mapping_h5}")
    print(f"Summary written to: {config.output_dir / 'hybrid_prep_summary.json'}")


if __name__ == "__main__":
    main()
