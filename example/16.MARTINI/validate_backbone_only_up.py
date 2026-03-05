#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


REQUIRED_BACKBONE_NODES = (
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
)

FORBIDDEN_SIDECHAIN_NODES = (
    "rotamer",
    "placement_fixed_scalar",
    "placement_fixed_point_vector_only",
    "placement_point_vector_only",
)

BACKBONE_CLASSES = ("N", "CA", "C", "O")
DOPC_COVERED_TYPES = ("Q0", "Qa", "Na", "C1", "C3")

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
    "r_angstrom",
    "target_potential",
    "fitted_potential",
    "r_min_angstrom",
    "well_depth",
    "base_depth",
    "role_scale",
    "sigma_nm",
    "epsilon_p2",
    "epsilon_proxy",
    "fit_rmse",
)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Validate backbone-only production stage (.up), depth interaction "
            "table artifacts, and radial backbone cross-table artifacts."
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
    p.add_argument("--cross-table-csv", type=Path, default=None, help="Optional radial cross table CSV")
    p.add_argument("--cross-table-meta", type=Path, default=None, help="Optional radial cross table metadata JSON")
    p.add_argument("--cross-artifact", type=Path, default=None, help="Optional radial cross table HDF5 artifact")
    return p.parse_args()


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def decode_bytes_array(arr):
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


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

        if "martini_rbm_cross_potential" in pot:
            cross = pot["martini_rbm_cross_potential"]
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
                        f"Missing cross-potential dataset: /input/potential/martini_rbm_cross_potential/{name}"
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
            if np.any(~np.isfinite(radial_centers)) or np.any(~np.isfinite(radial_widths)) or np.any(radial_widths <= 0):
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
            if "protein_classes" in cross:
                protein_classes = decode_bytes_array(cross["protein_classes"][:]).tolist()
                if len(protein_classes) != weights.shape[0]:
                    raise ValueError("Cross-potential protein_classes length mismatch")
            if "environment_classes" in cross:
                env_classes = decode_bytes_array(cross["environment_classes"][:]).tolist()
                if len(env_classes) != weights.shape[1]:
                    raise ValueError("Cross-potential environment_classes length mismatch")

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
        missing = [c for c in REQUIRED_DEPTH_TABLE_COLUMNS if c not in reader.fieldnames]
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
                try:
                    value = float(row[key])
                except ValueError as exc:
                    raise ValueError(f"CSV row {i} field {key} is not numeric") from exc
                if not np.isfinite(value):
                    raise ValueError(f"CSV row {i} field {key} is not finite")

    return {
        "n_rows": len(rows),
    }


def validate_cross_table_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("Cross-table CSV has no header row")
        missing = [c for c in REQUIRED_CROSS_TABLE_COLUMNS if c not in reader.fieldnames]
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
            if not protein_class or not env_type:
                raise ValueError(f"Cross-table CSV row {i} has empty class/type")
            protein_classes.add(protein_class)
            env_types.add(env_type)

            key = (protein_class, env_type, row["r_angstrom"])
            if key in seen:
                raise ValueError(f"Cross-table CSV duplicate row key: {key}")
            seen.add(key)

            for field in REQUIRED_CROSS_TABLE_COLUMNS:
                if field in ("protein_class", "env_type"):
                    continue
                try:
                    value = float(row[field])
                except ValueError as exc:
                    raise ValueError(f"Cross-table CSV row {i} field {field} is not numeric") from exc
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
    tp = data["term_presence"]
    for key in ("backbone_cb_energy", "hbond_hb_energy"):
        if key not in tp or not bool(tp[key]):
            raise ValueError(f"Metadata term_presence indicates missing required term: {key}")

    return {
        "term_presence": tp,
    }


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
    for key in DOPC_COVERED_TYPES:
        if key not in calibration["covered_depth_types"]:
            raise ValueError(f"Cross-table metadata missing covered type: {key}")

    overlap = data["overlap_report"]
    if "covered_types" not in overlap or "pearson_r" not in overlap:
        raise ValueError("Cross-table metadata overlap_report is incomplete")

    fit_summary = data["fit_rmse_summary"]
    for key in ("min", "median", "max"):
        if key not in fit_summary:
            raise ValueError(f"Cross-table metadata fit_rmse_summary missing {key}")
        if not np.isfinite(float(fit_summary[key])):
            raise ValueError(f"Cross-table metadata fit_rmse_summary {key} is not finite")

    return {
        "n_env_type": len(env_types),
        "fit_rmse_summary": fit_summary,
    }


def validate_cross_artifact(path: Path):
    with h5py.File(path, "r") as h5:
        required = (
            "/classes/protein",
            "/classes/environment",
            "/rbm/cross/radial_centers",
            "/rbm/cross/radial_widths",
            "/rbm/cross/weights",
        )
        for key in required:
            if key not in h5:
                raise ValueError(f"Cross artifact missing {key}")

        protein_classes = decode_bytes_array(h5["/classes/protein"][:]).tolist()
        env_classes = decode_bytes_array(h5["/classes/environment"][:]).tolist()
        centers = h5["/rbm/cross/radial_centers"][:]
        widths = h5["/rbm/cross/radial_widths"][:]
        weights = h5["/rbm/cross/weights"][:]

        schema = str(h5["/meta"].attrs.get("schema", ""))
        if schema != "martini_backbone_cross_v1":
            raise ValueError(f"Unexpected cross artifact schema: {schema}")
        if protein_classes != list(BACKBONE_CLASSES):
            raise ValueError(f"Cross artifact protein classes mismatch: {protein_classes}")
        if len(env_classes) != 38:
            raise ValueError(f"Cross artifact expected 38 environment classes, got {len(env_classes)}")
        if centers.ndim != 1 or widths.ndim != 1:
            raise ValueError("Cross artifact radial centers/widths must be rank-1")
        if weights.shape != (len(BACKBONE_CLASSES), len(env_classes), centers.shape[0]):
            raise ValueError(
                f"Cross artifact weights shape mismatch: {weights.shape}"
            )
        if not np.all(np.isfinite(centers)) or not np.all(np.isfinite(widths)) or not np.all(np.isfinite(weights)):
            raise ValueError("Cross artifact contains non-finite values")
        if np.any(widths <= 0.0):
            raise ValueError("Cross artifact radial_widths must be positive")

    return {
        "n_env_type": len(env_classes),
        "n_radial": int(centers.shape[0]),
    }


def main():
    args = parse_args()
    if (
        args.up_file is None
        and args.table_csv is None
        and args.table_meta is None
        and args.cross_table_csv is None
        and args.cross_table_meta is None
        and args.cross_artifact is None
    ):
        raise ValueError(
            "Provide at least one artifact: up_file, depth table artifacts, and/or cross-table artifacts"
        )

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

    cross_table_summary = None
    if args.cross_table_csv is not None:
        cross_table_csv = require_file(args.cross_table_csv.resolve())
        cross_table_summary = validate_cross_table_csv(cross_table_csv)

    cross_meta_summary = None
    if args.cross_table_meta is not None:
        cross_table_meta = require_file(args.cross_table_meta.resolve())
        cross_meta_summary = validate_cross_table_meta(cross_table_meta)

    cross_artifact_summary = None
    if args.cross_artifact is not None:
        cross_artifact = require_file(args.cross_artifact.resolve())
        cross_artifact_summary = validate_cross_artifact(cross_artifact)

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
            f"fit_rmse_min={float(rmse['min']):.6f} "
            f"fit_rmse_median={float(rmse['median']):.6f} "
            f"fit_rmse_max={float(rmse['max']):.6f}"
        )
    if cross_artifact_summary is not None:
        print(
            "  cross_artifact "
            f"env_types={cross_artifact_summary['n_env_type']} "
            f"n_radial={cross_artifact_summary['n_radial']}"
        )


if __name__ == "__main__":
    main()
