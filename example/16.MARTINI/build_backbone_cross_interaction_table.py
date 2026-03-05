#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from prepare_system_lib import read_martini3_nonbond_params


BACKBONE_CLASSES = ("N", "CA", "C", "O")
SCRIPT_DIR = Path(__file__).resolve().parent
ROLE_PROXIES = {
    "N": "Nd",
    "CA": "C3",
    "C": "P4",
    "O": "Na",
}
BASE_PROXY = "P2"
DOPC_COVERED_TYPES = ("Q0", "Qa", "Na", "C1", "C3")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a radial dry-MARTINI vs Upside-backbone interaction table "
            "using the DOPC dry-type trend from the existing depth table."
        )
    )
    parser.add_argument(
        "--ff-itp",
        type=Path,
        default=SCRIPT_DIR / "ff_dry" / "dry_martini_v2.1.itp",
        help="Dry MARTINI parameter file with [ atomtypes ] and [ nonbond_params ].",
    )
    parser.add_argument(
        "--depth-table-csv",
        type=Path,
        default=SCRIPT_DIR / "outputs" / "depth_interaction_table.csv",
        help="Existing DOPC depth interaction table used for calibration.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.csv",
        help="Output sampled table CSV.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.meta.json",
        help="Output metadata JSON.",
    )
    parser.add_argument(
        "--output-h5",
        type=Path,
        default=SCRIPT_DIR / "outputs" / "backbone_cross_interaction_table.h5",
        help="Output HDF5 artifact for martini_rbm_cross_potential.",
    )
    parser.add_argument(
        "--radial-min",
        type=float,
        default=2.0,
        help="Minimum sampled distance in Angstrom.",
    )
    parser.add_argument(
        "--radial-max",
        type=float,
        default=12.0,
        help="Maximum sampled distance in Angstrom.",
    )
    parser.add_argument(
        "--sample-step",
        type=float,
        default=0.1,
        help="Sample spacing for the exported radial table in Angstrom.",
    )
    parser.add_argument(
        "--n-radial",
        type=int,
        default=12,
        help="Number of Gaussian radial basis centers in the runtime artifact.",
    )
    parser.add_argument(
        "--ridge",
        type=float,
        default=1e-4,
        help="Ridge regularization used when fitting Gaussian-basis weights.",
    )
    return parser.parse_args()


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def repo_relative(path: Path, repo_root: Path):
    try:
        return str(path.resolve().relative_to(repo_root))
    except ValueError:
        return str(path.resolve())


def parse_atomtypes(path: Path):
    atomtypes = []
    in_atomtypes = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("[") and line.endswith("]"):
            in_atomtypes = line.strip("[]").strip().lower() == "atomtypes"
            continue
        if not in_atomtypes or not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 6 and parts[3] == "A":
            atomtypes.append(parts[0])

    if not atomtypes:
        raise ValueError(f"No atom types parsed from {path}")
    return atomtypes


def aggregate_depth_table(path: Path):
    per_type = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            bead_type = row["bead_type"].strip()
            entry = per_type.setdefault(
                bead_type,
                {
                    "count": 0,
                    "depth_sum": 0.0,
                    "hb_sum": 0.0,
                    "cb_sum": 0.0,
                },
            )
            entry["count"] += 1
            entry["depth_sum"] += float(row["mean_depth"])
            entry["hb_sum"] += -float(row["hb_mean"])
            entry["cb_sum"] += float(row["cb_backbone_mean"])

    if not per_type:
        raise ValueError(f"No rows parsed from {path}")

    out = {}
    for bead_type, entry in per_type.items():
        count = float(entry["count"])
        out[bead_type] = {
            "n_rows": int(entry["count"]),
            "mean_depth": entry["depth_sum"] / count,
            "mean_hb_attraction": entry["hb_sum"] / count,
            "mean_cb_value": entry["cb_sum"] / count,
        }
    return out


def pair_param(param_table, type1, type2):
    if (type1, type2) not in param_table:
        raise KeyError(f"Missing nonbonded parameter pair ({type1}, {type2})")
    return param_table[(type1, type2)]


def build_anchor_model(depth_by_type, param_table):
    required = ("Q0", "Qa", "Na", "C1", "C3")
    missing = [key for key in required if key not in depth_by_type]
    if missing:
        raise ValueError(f"Depth table is missing required DOPC dry types: {missing}")

    q_depth = 0.5 * (
        depth_by_type["Q0"]["mean_hb_attraction"] + depth_by_type["Qa"]["mean_hb_attraction"]
    )
    anchor_points = [
        (pair_param(param_table, BASE_PROXY, "Q0")[1], q_depth, "Q"),
        (pair_param(param_table, BASE_PROXY, "Na")[1], depth_by_type["Na"]["mean_hb_attraction"], "N"),
        (pair_param(param_table, BASE_PROXY, "C1")[1], depth_by_type["C1"]["mean_hb_attraction"], "C_light"),
        (pair_param(param_table, BASE_PROXY, "C3")[1], depth_by_type["C3"]["mean_hb_attraction"], "C_deep"),
    ]
    x = np.asarray([item[0] for item in anchor_points], dtype=np.float64)
    y = np.asarray([item[1] for item in anchor_points], dtype=np.float64)
    labels = [item[2] for item in anchor_points]
    if np.any(np.diff(x) <= 0):
        raise ValueError("Anchor epsilon values must be strictly increasing")
    return x, y, labels


def pchip_slopes(x, y):
    n = x.shape[0]
    h = np.diff(x)
    delta = np.diff(y) / h
    m = np.zeros(n, dtype=np.float64)

    if n == 2:
        m[:] = delta[0]
        return m

    for k in range(1, n - 1):
        if delta[k - 1] == 0.0 or delta[k] == 0.0 or np.sign(delta[k - 1]) != np.sign(delta[k]):
            m[k] = 0.0
        else:
            w1 = 2.0 * h[k] + h[k - 1]
            w2 = h[k] + 2.0 * h[k - 1]
            m[k] = (w1 + w2) / ((w1 / delta[k - 1]) + (w2 / delta[k]))

    m0 = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    if np.sign(m0) != np.sign(delta[0]):
        m0 = 0.0
    elif np.sign(delta[0]) != np.sign(delta[1]) and abs(m0) > 3.0 * abs(delta[0]):
        m0 = 3.0 * delta[0]
    m[0] = m0

    mn = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    if np.sign(mn) != np.sign(delta[-1]):
        mn = 0.0
    elif np.sign(delta[-1]) != np.sign(delta[-2]) and abs(mn) > 3.0 * abs(delta[-1]):
        mn = 3.0 * delta[-1]
    m[-1] = mn

    return m


def pchip_interpolate(x, y, xq):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xq = np.asarray(xq, dtype=np.float64)
    m = pchip_slopes(x, y)

    idx = np.searchsorted(x, xq, side="right") - 1
    idx = np.clip(idx, 0, x.shape[0] - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    m0 = m[idx]
    m1 = m[idx + 1]
    h = x1 - x0
    t = (xq - x0) / h

    h00 = (2.0 * t**3) - (3.0 * t**2) + 1.0
    h10 = (t**3) - (2.0 * t**2) + t
    h01 = (-2.0 * t**3) + (3.0 * t**2)
    h11 = (t**3) - (t**2)
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def interpolate_depth(epsilon_value, anchor_x, anchor_y):
    epsilon_value = float(epsilon_value)
    if epsilon_value <= anchor_x[0]:
        return float(anchor_y[0])
    if epsilon_value <= anchor_x[-1]:
        return float(pchip_interpolate(anchor_x, anchor_y, np.array([epsilon_value], dtype=np.float64))[0])

    slope = (anchor_y[-1] - anchor_y[-2]) / (anchor_x[-1] - anchor_x[-2])
    value = anchor_y[-1] + slope * (epsilon_value - anchor_x[-1])
    cap = 1.25 * anchor_y[-1]
    return float(min(value, cap))


def merge_control_points(x, y):
    order = np.argsort(x)
    x_sorted = np.asarray(x, dtype=np.float64)[order]
    y_sorted = np.asarray(y, dtype=np.float64)[order]

    merged_x = []
    merged_y = []
    i = 0
    while i < x_sorted.shape[0]:
        j = i + 1
        values = [y_sorted[i]]
        while j < x_sorted.shape[0] and abs(x_sorted[j] - x_sorted[i]) < 1e-8:
            values.append(y_sorted[j])
            j += 1
        merged_x.append(float(x_sorted[i]))
        merged_y.append(float(np.mean(np.asarray(values, dtype=np.float64))))
        i = j

    if len(merged_x) < 2:
        raise ValueError("Need at least two unique control points")
    return np.asarray(merged_x, dtype=np.float64), np.asarray(merged_y, dtype=np.float64)


def radial_control_points(r_min, depth, radial_min, radial_max):
    x = np.asarray(
        [
            max(radial_min, 0.72 * r_min),
            max(radial_min + 0.2, 0.90 * r_min),
            r_min,
            min(radial_max, 1.35 * r_min),
            radial_max,
        ],
        dtype=np.float64,
    )
    y = np.asarray(
        [
            3.0 * depth,
            1.2 * depth,
            -depth,
            -0.30 * depth,
            0.0,
        ],
        dtype=np.float64,
    )
    return merge_control_points(x, y)


def gaussian_basis(r, centers, widths):
    delta = (r[:, None] - centers[None, :]) / widths[None, :]
    return np.exp(-0.5 * delta * delta)


def fit_gaussian_weights(target_values, centers, widths, radial_grid, ridge):
    basis = gaussian_basis(radial_grid, centers, widths)
    gram = basis.T @ basis
    rhs = basis.T @ target_values
    reg = ridge * np.eye(centers.shape[0], dtype=np.float64)
    weights = np.linalg.solve(gram + reg, rhs)
    fitted = basis @ weights
    rmse = float(np.sqrt(np.mean((fitted - target_values) ** 2)))
    return weights, fitted, rmse


def decode_ring_scale(atom_type):
    return 0.75 if atom_type.startswith("S") else 1.0


def build_type_records(atomtypes, param_table, anchor_x, anchor_y, radial_grid, centers, widths, ridge):
    rows = []
    weight_tensor = np.zeros((len(BACKBONE_CLASSES), len(atomtypes), len(centers)), dtype=np.float32)
    pair_summaries = {}

    for env_index, env_type in enumerate(atomtypes):
        sigma_base, epsilon_base = pair_param(param_table, BASE_PROXY, env_type)
        base_depth = interpolate_depth(epsilon_base, anchor_x, anchor_y) * decode_ring_scale(env_type)

        for role_index, role in enumerate(BACKBONE_CLASSES):
            proxy = ROLE_PROXIES[role]
            sigma_role = sigma_base
            epsilon_role = epsilon_base
            if (proxy, env_type) in param_table:
                sigma_role, epsilon_role = param_table[(proxy, env_type)]

            role_scale = np.clip(epsilon_role / epsilon_base, 0.6, 1.4) if epsilon_base > 0.0 else 1.0
            well_depth = float(base_depth * role_scale)
            r_min = float((2.0 ** (1.0 / 6.0)) * sigma_role * 10.0)
            ctrl_x, ctrl_y = radial_control_points(r_min, well_depth, float(radial_grid[0]), float(radial_grid[-1]))
            target = pchip_interpolate(ctrl_x, ctrl_y, radial_grid)
            weights, fitted, rmse = fit_gaussian_weights(target, centers, widths, radial_grid, ridge)
            weight_tensor[role_index, env_index, :] = weights.astype(np.float32)

            pair_summaries[f"{role}:{env_type}"] = {
                "r_min_angstrom": r_min,
                "well_depth": well_depth,
                "base_depth": float(base_depth),
                "role_scale": float(role_scale),
                "sigma_nm": float(sigma_role),
                "epsilon_p2": float(epsilon_base),
                "epsilon_proxy": float(epsilon_role),
                "fit_rmse": rmse,
                "control_r": [float(x) for x in ctrl_x],
                "control_v": [float(y) for y in ctrl_y],
            }

            for r_value, target_value, fitted_value in zip(radial_grid, target, fitted):
                rows.append(
                    {
                        "protein_class": role,
                        "env_type": env_type,
                        "r_angstrom": float(r_value),
                        "target_potential": float(target_value),
                        "fitted_potential": float(fitted_value),
                        "r_min_angstrom": r_min,
                        "well_depth": well_depth,
                        "base_depth": float(base_depth),
                        "role_scale": float(role_scale),
                        "sigma_nm": float(sigma_role),
                        "epsilon_p2": float(epsilon_base),
                        "epsilon_proxy": float(epsilon_role),
                        "fit_rmse": rmse,
                    }
                )

    return rows, weight_tensor, pair_summaries


def build_overlap_report(depth_by_type, pair_summaries):
    overlap_rows = []
    predicted = []
    observed = []
    for env_type in DOPC_COVERED_TYPES:
        role_depths = [pair_summaries[f"{role}:{env_type}"]["well_depth"] for role in BACKBONE_CLASSES]
        pred_value = float(np.mean(np.asarray(role_depths, dtype=np.float64)))
        obs_value = float(depth_by_type[env_type]["mean_hb_attraction"])
        overlap_rows.append(
            {
                "env_type": env_type,
                "observed_hb_attraction": obs_value,
                "predicted_role_mean_well_depth": pred_value,
                "mean_depth": float(depth_by_type[env_type]["mean_depth"]),
                "n_rows": int(depth_by_type[env_type]["n_rows"]),
            }
        )
        predicted.append(pred_value)
        observed.append(obs_value)

    predicted_arr = np.asarray(predicted, dtype=np.float64)
    observed_arr = np.asarray(observed, dtype=np.float64)
    corr = float(np.corrcoef(predicted_arr, observed_arr)[0, 1]) if predicted_arr.size >= 2 else 1.0
    return {
        "covered_types": overlap_rows,
        "pearson_r": corr,
    }


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_metadata(path: Path, metadata):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def write_artifact(path: Path, protein_classes, env_types, centers, widths, weights, metadata):
    path.parent.mkdir(parents=True, exist_ok=True)
    str_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["schema"] = "martini_backbone_cross_v1"
        meta.attrs["radial_min"] = np.float32(metadata["radial_grid"]["min"])
        meta.attrs["radial_max"] = np.float32(metadata["radial_grid"]["max"])
        meta.attrs["sample_step"] = np.float32(metadata["radial_grid"]["step"])
        meta.attrs["n_env_type"] = np.int32(len(env_types))
        meta.attrs["n_protein_class"] = np.int32(len(protein_classes))

        sources = h5.create_group("sources")
        sources.attrs["ff_itp"] = metadata["sources"]["ff_itp"]
        sources.attrs["depth_table_csv"] = metadata["sources"]["depth_table_csv"]

        classes = h5.create_group("classes")
        classes.create_dataset("protein", data=np.asarray(protein_classes, dtype=str_dtype))
        classes.create_dataset("environment", data=np.asarray(env_types, dtype=str_dtype))

        rbm = h5.create_group("rbm")
        cross = rbm.create_group("cross")
        cross.attrs["cutoff"] = np.float32(metadata["radial_grid"]["max"])
        cross.create_dataset("radial_centers", data=np.asarray(centers, dtype=np.float32))
        cross.create_dataset("radial_widths", data=np.asarray(widths, dtype=np.float32))
        cross.create_dataset("weights", data=np.asarray(weights, dtype=np.float32))

        views = h5.create_group("views")
        views.create_dataset("protein_classes", data=np.asarray(protein_classes, dtype=str_dtype))
        views.create_dataset("environment_classes", data=np.asarray(env_types, dtype=str_dtype))


def main():
    args = parse_args()
    ff_itp = require_file(args.ff_itp.resolve())
    depth_table_csv = require_file(args.depth_table_csv.resolve())
    if args.radial_min <= 0.0 or args.radial_max <= args.radial_min:
        raise ValueError("Require 0 < radial_min < radial_max")
    if args.sample_step <= 0.0:
        raise ValueError("--sample-step must be positive")
    if args.n_radial <= 0:
        raise ValueError("--n-radial must be positive")
    if args.ridge < 0.0:
        raise ValueError("--ridge must be non-negative")

    repo_root = Path(__file__).resolve().parents[2]
    atomtypes = parse_atomtypes(ff_itp)
    param_table = read_martini3_nonbond_params(str(ff_itp))
    depth_by_type = aggregate_depth_table(depth_table_csv)
    anchor_x, anchor_y, anchor_labels = build_anchor_model(depth_by_type, param_table)

    radial_grid = np.arange(args.radial_min, args.radial_max + 0.5 * args.sample_step, args.sample_step)
    if radial_grid[-1] > args.radial_max:
        radial_grid[-1] = args.radial_max
    centers = np.linspace(args.radial_min, args.radial_max, args.n_radial, dtype=np.float64)
    if args.n_radial == 1:
        widths = np.asarray([max(1.0, args.radial_max - args.radial_min)], dtype=np.float64)
    else:
        widths = np.full(args.n_radial, centers[1] - centers[0], dtype=np.float64)

    rows, weight_tensor, pair_summaries = build_type_records(
        atomtypes=atomtypes,
        param_table=param_table,
        anchor_x=anchor_x,
        anchor_y=anchor_y,
        radial_grid=radial_grid,
        centers=centers,
        widths=widths,
        ridge=args.ridge,
    )
    overlap_report = build_overlap_report(depth_by_type, pair_summaries)

    fit_rmses = np.asarray([item["fit_rmse"] for item in pair_summaries.values()], dtype=np.float64)
    metadata = {
        "sources": {
            "ff_itp": repo_relative(ff_itp, repo_root),
            "depth_table_csv": repo_relative(depth_table_csv, repo_root),
        },
        "protein_classes": list(BACKBONE_CLASSES),
        "environment_types": atomtypes,
        "radial_grid": {
            "min": float(args.radial_min),
            "max": float(args.radial_max),
            "step": float(args.sample_step),
            "n_samples": int(radial_grid.shape[0]),
            "n_radial_basis": int(args.n_radial),
            "ridge": float(args.ridge),
        },
        "calibration": {
            "anchor_labels": anchor_labels,
            "anchor_epsilons": [float(x) for x in anchor_x],
            "anchor_depths": [float(y) for y in anchor_y],
            "covered_depth_types": {
                key: depth_by_type[key]
                for key in DOPC_COVERED_TYPES
                if key in depth_by_type
            },
        },
        "role_proxies": ROLE_PROXIES,
        "pair_summary": pair_summaries,
        "fit_rmse_summary": {
            "min": float(np.min(fit_rmses)),
            "median": float(np.median(fit_rmses)),
            "max": float(np.max(fit_rmses)),
        },
        "overlap_report": overlap_report,
    }

    write_csv(args.output_csv, rows)
    write_metadata(args.output_json, metadata)
    write_artifact(
        path=args.output_h5,
        protein_classes=BACKBONE_CLASSES,
        env_types=atomtypes,
        centers=centers,
        widths=widths,
        weights=weight_tensor,
        metadata=metadata,
    )

    print(f"Wrote radial cross table CSV: {args.output_csv}")
    print(f"Wrote radial cross table metadata: {args.output_json}")
    print(f"Wrote radial cross table artifact: {args.output_h5}")
    print(f"Covered dry types from depth table: {', '.join(DOPC_COVERED_TYPES)}")
    print(
        "Fit RMSE summary: "
        f"min={metadata['fit_rmse_summary']['min']:.6f} "
        f"median={metadata['fit_rmse_summary']['median']:.6f} "
        f"max={metadata['fit_rmse_summary']['max']:.6f}"
    )


if __name__ == "__main__":
    main()
