#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


REQUIRED_BACKBONE_NODES = (
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
)

FORBIDDEN_SIDECHAIN_NODES = (
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
    "placement_point_vector_only",
    "affine_alignment",
)

REQUIRED_TABLE_COLUMNS = (
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


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Validate backbone-only production stage (.up) and optional depth "
            "interaction-table artifacts."
        )
    )
    p.add_argument(
        "up_file",
        type=Path,
        nargs="?",
        default=None,
        help="Optional stage .up file to validate",
    )
    p.add_argument("--table-csv", type=Path, default=None, help="Optional depth table CSV")
    p.add_argument("--table-meta", type=Path, default=None, help="Optional depth table metadata JSON")
    return p.parse_args()


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path


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
            ds = pot[group_name][dataset_name][:]
            arr = np.asarray(ds)
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

        seq = inp["sequence"][:]
        if seq.ndim != 1:
            raise ValueError("/input/sequence must be rank-1")
        if seq.shape[0] == 0:
            raise ValueError("/input/sequence is empty")

    return {
        "n_atom": n_atom,
        "n_res": int(seq.shape[0]),
    }


def validate_table_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header row")
        missing = [c for c in REQUIRED_TABLE_COLUMNS if c not in reader.fieldnames]
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

            for key in REQUIRED_TABLE_COLUMNS:
                if key in ("bead_name", "bead_type"):
                    continue
                try:
                    value = float(row[key])
                except ValueError as exc:
                    raise ValueError(f"CSV row {i} field {key} is not numeric") from exc
                if not np.isfinite(value):
                    raise ValueError(f"CSV row {i} field {key} is not finite")

    return {
        "n_rows": len(rows),
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
    tp = data["term_presence"]
    for key in ("backbone_cb_energy", "hbond_hb_energy"):
        if key not in tp or not bool(tp[key]):
            raise ValueError(f"Metadata term_presence indicates missing required term: {key}")

    return {
        "term_presence": tp,
    }


def main():
    args = parse_args()
    if args.up_file is None and args.table_csv is None and args.table_meta is None:
        raise ValueError("Provide at least one artifact: up_file and/or --table-csv/--table-meta")

    up_summary = None
    up_file = None
    if args.up_file is not None:
        up_file = require_file(args.up_file.resolve())
        up_summary = validate_up_file(up_file)

    table_summary = None
    if args.table_csv is not None:
        table_csv = require_file(args.table_csv.resolve())
        table_summary = validate_table_csv(table_csv)

    meta_summary = None
    if args.table_meta is not None:
        table_meta = require_file(args.table_meta.resolve())
        meta_summary = validate_table_meta(table_meta)

    if up_summary is not None:
        print(f"OK: {up_file}")
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


if __name__ == "__main__":
    main()
