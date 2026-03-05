#!/usr/bin/env python3

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


SIDECHAIN_DATASET_CANDIDATES = (
    "sc_energy",
    "sidechain_energy",
    "sidechain",
)

SHEET_DATASET_CANDIDATES = (
    "sheet_energy",
    "sheet",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a depth-dependent dry-MARTINI bead interaction table using "
            "Upside membrane channels (cb_energy for backbone, hb_energy for hbond)."
        )
    )
    parser.add_argument(
        "--bilayer-pdb",
        type=Path,
        default=Path("pdb/bilayer.MARTINI.pdb"),
        help="Bilayer MARTINI PDB file.",
    )
    parser.add_argument(
        "--lipid-itp",
        type=Path,
        default=Path("ff_dry/dry_martini_v2.1_lipids.itp"),
        help="Lipid ITP used to map bead name -> MARTINI type.",
    )
    parser.add_argument(
        "--membrane-h5",
        type=Path,
        default=Path("../../parameters/ff_2.1/membrane.h5"),
        help="Upside membrane potential file.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("outputs/depth_interaction_table.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("outputs/depth_interaction_table.meta.json"),
        help="Output JSON metadata path.",
    )
    parser.add_argument(
        "--lipid-resnames",
        nargs="+",
        default=["DOPC", "DOP"],
        help="Lipid residue names to include from PDB.",
    )
    return parser.parse_args()


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def parse_itp_atom_types(path: Path, allowed_molecules):
    allowed_molecules = {x.upper() for x in allowed_molecules}
    atom_to_type = {}

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    current_section = ""
    current_molecule = None
    i = 0
    while i < len(lines):
        raw = lines[i]
        line = raw.split(";", 1)[0].strip()
        i += 1
        if not line:
            continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line.strip("[]").strip().lower()
            if current_section == "moleculetype":
                current_molecule = None
            continue

        if current_section == "moleculetype":
            if current_molecule is not None:
                continue
            current_molecule = line.split()[0].upper()
            continue

        if current_section != "atoms":
            continue
        if current_molecule not in allowed_molecules:
            continue

        parts = line.split()
        if len(parts) < 5:
            continue
        bead_type = parts[1].strip()
        bead_name = parts[4].strip().upper()
        atom_to_type.setdefault(bead_name, bead_type)

    if not atom_to_type:
        raise ValueError(
            f"No atoms parsed from {path} for molecules {sorted(allowed_molecules)}"
        )
    return atom_to_type


def parse_bilayer_depths(path: Path, lipid_resnames):
    lipid_resnames = {x.upper() for x in lipid_resnames}
    bead_z = defaultdict(list)
    all_lipid_z = []

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if resname not in lipid_resnames:
                continue
            atom_name = line[12:16].strip().upper()
            try:
                z = float(line[46:54])
            except ValueError:
                continue
            bead_z[atom_name].append(z)
            all_lipid_z.append(z)

    if not all_lipid_z:
        raise ValueError(f"No lipid atoms found in {path} for resnames {sorted(lipid_resnames)}")

    z_center = float(np.mean(np.asarray(all_lipid_z, dtype=np.float64)))

    stats = {}
    for bead_name, z_values in bead_z.items():
        z_arr = np.asarray(z_values, dtype=np.float64)
        depths = np.abs(z_arr - z_center)
        stats[bead_name] = {
            "count": int(depths.size),
            "mean_depth": float(np.mean(depths)),
            "std_depth": float(np.std(depths)),
        }
    return z_center, stats


def load_membrane_channels(path: Path):
    with h5py.File(path, "r") as h5:
        if "cb_energy" not in h5:
            raise ValueError(f"Missing cb_energy in {path}")
        if "hb_energy" not in h5:
            raise ValueError(f"Missing hb_energy in {path}")

        cb_energy = h5["cb_energy"][:].astype(np.float64)
        hb_energy = h5["hb_energy"][:].astype(np.float64)
        cb_z_min = float(h5.attrs["cb_z_min"])
        cb_z_max = float(h5.attrs["cb_z_max"])
        hb_z_min = float(h5.attrs["hb_z_min"])
        hb_z_max = float(h5.attrs["hb_z_max"])
        cb_names = [x.decode("ascii") if isinstance(x, bytes) else str(x) for x in h5["names"][:]]

        sidechain_present = any(name in h5 for name in SIDECHAIN_DATASET_CANDIDATES)
        sheet_present = any(name in h5 for name in SHEET_DATASET_CANDIDATES)

    if cb_energy.ndim != 3:
        raise ValueError(f"cb_energy must be rank-3, got shape {cb_energy.shape}")
    if hb_energy.ndim != 3:
        raise ValueError(f"hb_energy must be rank-3, got shape {hb_energy.shape}")
    if hb_energy.shape[0] != 2 or hb_energy.shape[1] != 2:
        raise ValueError(f"hb_energy expected shape [2,2,n], got {hb_energy.shape}")

    cb_z = np.linspace(cb_z_min, cb_z_max, cb_energy.shape[2], dtype=np.float64)
    hb_z = np.linspace(hb_z_min, hb_z_max, hb_energy.shape[2], dtype=np.float64)

    term_presence = {
        "backbone_cb_energy": True,
        "hbond_hb_energy": True,
        "sidechain_term": sidechain_present,
        "sheet_term": sheet_present,
    }
    return cb_z, hb_z, cb_energy, hb_energy, cb_names, term_presence


def interp_symmetric(z_grid, y_values, depth):
    # Depth is defined as |z - z_center|. Evaluate both +depth and -depth and average.
    d = abs(float(depth))
    yp = np.interp(d, z_grid, y_values, left=float(y_values[0]), right=float(y_values[-1]))
    yn = np.interp(-d, z_grid, y_values, left=float(y_values[0]), right=float(y_values[-1]))
    return 0.5 * (yp + yn)


def evaluate_at_depth(depth, cb_z, hb_z, cb_energy, hb_energy):
    cb_lvl0 = np.array([interp_symmetric(cb_z, cb_energy[i, 0, :], depth) for i in range(cb_energy.shape[0])])
    cb_lvl1 = np.array([interp_symmetric(cb_z, cb_energy[i, 1, :], depth) for i in range(cb_energy.shape[0])])

    hb_donor_unbound = interp_symmetric(hb_z, hb_energy[0, 0, :], depth)
    hb_donor_bound = interp_symmetric(hb_z, hb_energy[0, 1, :], depth)
    hb_acceptor_unbound = interp_symmetric(hb_z, hb_energy[1, 0, :], depth)
    hb_acceptor_bound = interp_symmetric(hb_z, hb_energy[1, 1, :], depth)

    hb_all = np.array(
        [hb_donor_unbound, hb_donor_bound, hb_acceptor_unbound, hb_acceptor_bound],
        dtype=np.float64,
    )

    return {
        "cb_backbone_burial0_mean": float(np.mean(cb_lvl0)),
        "cb_backbone_burial1_mean": float(np.mean(cb_lvl1)),
        "cb_backbone_mean": float(np.mean(np.concatenate([cb_lvl0, cb_lvl1]))),
        "hb_donor_unbound": float(hb_donor_unbound),
        "hb_donor_bound": float(hb_donor_bound),
        "hb_acceptor_unbound": float(hb_acceptor_unbound),
        "hb_acceptor_bound": float(hb_acceptor_bound),
        "hb_mean": float(np.mean(hb_all)),
    }


def main():
    args = parse_args()
    bilayer_pdb = require_file(args.bilayer_pdb.resolve())
    lipid_itp = require_file(args.lipid_itp.resolve())
    membrane_h5 = require_file(args.membrane_h5.resolve())

    atom_to_type = parse_itp_atom_types(lipid_itp, args.lipid_resnames)
    z_center, bead_depth_stats = parse_bilayer_depths(bilayer_pdb, args.lipid_resnames)
    cb_z, hb_z, cb_energy, hb_energy, cb_names, term_presence = load_membrane_channels(membrane_h5)

    rows = []
    for bead_name, stats in bead_depth_stats.items():
        depth = stats["mean_depth"]
        values = evaluate_at_depth(depth, cb_z, hb_z, cb_energy, hb_energy)
        rows.append(
            {
                "bead_name": bead_name,
                "bead_type": atom_to_type.get(bead_name, "UNKNOWN"),
                "mean_depth": stats["mean_depth"],
                "std_depth": stats["std_depth"],
                "count": stats["count"],
                **values,
            }
        )

    rows.sort(key=lambda x: x["mean_depth"], reverse=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
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
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "bilayer_pdb": str(bilayer_pdb),
        "lipid_itp": str(lipid_itp),
        "membrane_h5": str(membrane_h5),
        "lipid_resnames": [x.upper() for x in args.lipid_resnames],
        "z_center": z_center,
        "n_bead_names": len(rows),
        "cb_residue_names": cb_names,
        "depth_definition": "|z - z_center| (no PO4 anchor)",
        "term_presence": term_presence,
    }
    with args.output_json.open("w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)

    print(f"Wrote depth interaction table: {args.output_csv}")
    print(f"Wrote metadata: {args.output_json}")
    print(f"Bilayer z_center = {z_center:.6f} A")
    print("membrane.h5 term presence:")
    for key, value in term_presence.items():
        print(f"  {key}={int(bool(value))}")


if __name__ == "__main__":
    main()
