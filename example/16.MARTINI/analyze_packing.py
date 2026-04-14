#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np


TAIL_BONDS: tuple[tuple[str, str, str], ...] = (
    ("GL1", "C1A", "sn1_gl1_c1a"),
    ("C1A", "C2A", "sn1_c1a_c2a"),
    ("C2A", "D3A", "sn1_c2a_d3a"),
    ("D3A", "C4A", "sn1_d3a_c4a"),
    ("C4A", "C5A", "sn1_c4a_c5a"),
    ("GL2", "C1B", "sn2_gl2_c1b"),
    ("C1B", "C2B", "sn2_c1b_c2b"),
    ("C2B", "D3B", "sn2_c2b_d3b"),
    ("D3B", "C4B", "sn2_d3b_c4b"),
    ("C4B", "C5B", "sn2_c4b_c5b"),
)


def _normalize_output_positions(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr)
    if out.ndim == 4:
        if out.shape[-1] == 3:
            return np.asarray(out[:, 0, :, :], dtype=np.float64)
        if out.shape[1] == 1 and out.shape[2] == 3:
            return np.asarray(np.transpose(out[:, :, :, 0], (0, 2, 1)), dtype=np.float64)
    if out.ndim == 3:
        if out.shape[-1] == 3:
            return np.asarray(out, dtype=np.float64)
        if out.shape[1] == 3:
            return np.asarray(np.transpose(out, (0, 2, 1)), dtype=np.float64)
    raise ValueError(f"Unexpected output/pos shape: {arr.shape}")


def _normalize_time(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr)
    while out.ndim > 1:
        out = out[:, 0]
    return np.asarray(out, dtype=np.float64).reshape(-1)


def _normalize_box(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr)
    if out.ndim == 3:
        if out.shape[-1] == 3:
            return np.asarray(out[:, 0, :], dtype=np.float64)
        if out.shape[1] == 3:
            return np.asarray(np.transpose(out[:, :, 0], (0, 1)), dtype=np.float64)
    if out.ndim == 2 and out.shape[-1] == 3:
        return np.asarray(out, dtype=np.float64)
    raise ValueError(f"Unexpected output/box shape: {arr.shape}")


def _decode_str_array(dataset: h5py.Dataset) -> np.ndarray:
    arr = dataset[:]
    out: list[str] = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


def _load_box_frames(h5: h5py.File, n_frame: int) -> np.ndarray:
    if "output/box" in h5:
        return _normalize_box(h5["output/box"][:])
    grp = h5["input/potential/martini_potential"]
    box = np.array([grp.attrs["x_len"], grp.attrs["y_len"], grp.attrs["z_len"]], dtype=np.float64)
    return np.repeat(box[None, :], n_frame, axis=0)


def _minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def _select_lipid_molecules(atom_names: np.ndarray, molecule_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    po4_mask = atom_names == "PO4"
    if not np.any(po4_mask):
        raise ValueError("Could not identify lipid PO4 beads from atom_names")
    lipid_molecule_ids = np.unique(molecule_ids[po4_mask])
    lipid_mask = np.isin(molecule_ids, lipid_molecule_ids)
    po4_indices = np.where(lipid_mask & po4_mask)[0]
    if po4_indices.size == 0:
        raise ValueError("Could not infer lipid molecules from PO4 beads")
    return lipid_molecule_ids, po4_indices


def _build_tail_bond_index_map(
    atom_names: np.ndarray,
    molecule_ids: np.ndarray,
    lipid_molecule_ids: np.ndarray,
) -> dict[str, np.ndarray]:
    label_to_pairs: dict[str, list[tuple[int, int]]] = {label: [] for _, _, label in TAIL_BONDS}
    for molecule_id in lipid_molecule_ids:
        atom_idx = np.where(molecule_ids == molecule_id)[0]
        name_to_index = {str(atom_names[idx]): int(idx) for idx in atom_idx}
        for atom_a, atom_b, label in TAIL_BONDS:
            idx_a = name_to_index.get(atom_a)
            idx_b = name_to_index.get(atom_b)
            if idx_a is None or idx_b is None:
                continue
            label_to_pairs[label].append((idx_a, idx_b))
    out: dict[str, np.ndarray] = {}
    for label, pairs in label_to_pairs.items():
        if pairs:
            out[label] = np.asarray(pairs, dtype=np.int32)
    if not out:
        raise ValueError("Could not build any lipid tail-bond pairs for order analysis")
    return out


def _mean_and_std(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def analyze_stage7_packing(
    stage_file: Path,
    *,
    burn_in_fraction: float,
    burn_in_frames: int | None,
    target_apl: float | None,
    target_thickness: float | None,
    target_order: float | None,
) -> tuple[dict[str, Any], list[dict[str, float]]]:
    with h5py.File(stage_file, "r") as h5:
        if "output/pos" not in h5:
            raise ValueError(f"Missing output/pos in {stage_file}")
        pos = _normalize_output_positions(h5["output/pos"][:])
        n_frame_total, n_atom, _ = pos.shape
        if n_frame_total < 2:
            raise ValueError(f"Need at least 2 frames for packing analysis, found {n_frame_total}")

        if "output/time" in h5:
            times = _normalize_time(h5["output/time"][:])
        else:
            times = np.arange(n_frame_total, dtype=np.float64)
        if times.shape[0] != n_frame_total:
            raise ValueError("output/time length does not match output/pos frames")

        box = _load_box_frames(h5, n_frame_total)
        atom_names = _decode_str_array(h5["input/atom_names"])
        molecule_ids = np.asarray(h5["input/molecule_ids"][:], dtype=np.int32)

        lipid_molecule_ids, po4_indices = _select_lipid_molecules(atom_names, molecule_ids)
        bond_pairs = _build_tail_bond_index_map(atom_names, molecule_ids, lipid_molecule_ids)

        if burn_in_frames is None:
            burn_start = max(1, int(math.floor(burn_in_fraction * n_frame_total)))
        else:
            burn_start = int(burn_in_frames)
        burn_start = max(0, min(burn_start, n_frame_total - 1))

        pos = pos[burn_start:]
        times = times[burn_start:]
        box = box[burn_start:]
        n_frame_used = pos.shape[0]

        time_series: list[dict[str, float]] = []
        area_xy_values: list[float] = []
        apl_upper_values: list[float] = []
        apl_lower_values: list[float] = []
        apl_mean_values: list[float] = []
        thickness_values: list[float] = []
        upper_count_values: list[int] = []
        lower_count_values: list[int] = []
        order_values: list[float] = []
        order_by_label_values: dict[str, list[float]] = {label: [] for label in bond_pairs}

        for frame in range(n_frame_used):
            box_frame = box[frame]
            area_xy = float(box_frame[0] * box_frame[1])
            po4_z = pos[frame, po4_indices, 2]
            midplane_z = float(np.median(po4_z))
            upper_mask = po4_z >= midplane_z
            lower_mask = ~upper_mask
            n_upper = int(np.count_nonzero(upper_mask))
            n_lower = int(np.count_nonzero(lower_mask))
            if n_upper == 0 or n_lower == 0:
                raise ValueError("Leaflet split failed: one leaflet is empty")

            upper_z_mean = float(np.mean(po4_z[upper_mask]))
            lower_z_mean = float(np.mean(po4_z[lower_mask]))
            thickness = upper_z_mean - lower_z_mean
            apl_upper = area_xy / float(n_upper)
            apl_lower = area_xy / float(n_lower)
            apl_mean = 0.5 * (apl_upper + apl_lower)

            frame_order_means: list[float] = []
            for label, pairs in bond_pairs.items():
                disp = pos[frame, pairs[:, 1], :] - pos[frame, pairs[:, 0], :]
                disp = _minimum_image(disp, box_frame[None, :])
                norms = np.linalg.norm(disp, axis=1)
                valid = norms > 1e-12
                if not np.any(valid):
                    continue
                cos2 = (disp[valid, 2] / norms[valid]) ** 2
                s = 0.5 * (3.0 * cos2 - 1.0)
                label_mean = float(np.mean(s))
                order_by_label_values[label].append(label_mean)
                frame_order_means.append(label_mean)

            order_mean = float(np.mean(frame_order_means)) if frame_order_means else float("nan")

            time_series.append(
                {
                    "frame_index": float(burn_start + frame),
                    "time": float(times[frame]),
                    "box_x": float(box_frame[0]),
                    "box_y": float(box_frame[1]),
                    "box_z": float(box_frame[2]),
                    "area_xy": area_xy,
                    "n_leaflet_upper": float(n_upper),
                    "n_leaflet_lower": float(n_lower),
                    "apl_upper": apl_upper,
                    "apl_lower": apl_lower,
                    "apl_mean": apl_mean,
                    "po4_midplane_z": midplane_z,
                    "po4_upper_z_mean": upper_z_mean,
                    "po4_lower_z_mean": lower_z_mean,
                    "po4_thickness": thickness,
                    "tail_order_mean": order_mean,
                }
            )

            area_xy_values.append(area_xy)
            apl_upper_values.append(apl_upper)
            apl_lower_values.append(apl_lower)
            apl_mean_values.append(apl_mean)
            thickness_values.append(thickness)
            upper_count_values.append(n_upper)
            lower_count_values.append(n_lower)
            if math.isfinite(order_mean):
                order_values.append(order_mean)

        summary: dict[str, Any] = {
            "stage_file": str(stage_file),
            "n_atom": int(n_atom),
            "n_frames_total": int(n_frame_total),
            "n_frames_used": int(n_frame_used),
            "burn_in_frames": int(burn_start),
            "n_lipid_molecules": int(lipid_molecule_ids.shape[0]),
            "n_po4": int(po4_indices.shape[0]),
            "leaflet_upper_count_mean": float(np.mean(np.asarray(upper_count_values, dtype=np.float64))),
            "leaflet_lower_count_mean": float(np.mean(np.asarray(lower_count_values, dtype=np.float64))),
            "area_xy_angstrom2": _mean_and_std(area_xy_values),
            "apl_angstrom2": {
                "upper": _mean_and_std(apl_upper_values),
                "lower": _mean_and_std(apl_lower_values),
                "mean": _mean_and_std(apl_mean_values),
            },
            "po4_thickness_angstrom": _mean_and_std(thickness_values),
            "tail_order_parameter": _mean_and_std(order_values) if order_values else None,
            "tail_order_parameter_by_bond": {
                label: _mean_and_std(values)
                for label, values in sorted(order_by_label_values.items())
                if values
            },
        }
        if target_apl is not None:
            summary["target_apl_angstrom2"] = float(target_apl)
            summary["target_apl_abs_error_angstrom2"] = float(
                abs(summary["apl_angstrom2"]["mean"]["mean"] - float(target_apl))
            )
        if target_thickness is not None:
            summary["target_po4_thickness_angstrom"] = float(target_thickness)
            summary["target_po4_thickness_abs_error_angstrom"] = float(
                abs(summary["po4_thickness_angstrom"]["mean"] - float(target_thickness))
            )
        if target_order is not None and summary["tail_order_parameter"] is not None:
            summary["target_tail_order_parameter"] = float(target_order)
            summary["target_tail_order_abs_error"] = float(
                abs(summary["tail_order_parameter"]["mean"] - float(target_order))
            )
        return summary, time_series


def write_timeseries_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise ValueError("No time-series rows available to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze membrane packing metrics from a hybrid stage_7.0.up file."
    )
    parser.add_argument("stage_file", type=Path, help="Path to stage_7.0.up")
    parser.add_argument("--json-out", type=Path, default=None, help="Optional summary JSON output path")
    parser.add_argument("--csv-out", type=Path, default=None, help="Optional per-frame CSV output path")
    parser.add_argument(
        "--burn-in-fraction",
        type=float,
        default=0.20,
        help="Fraction of frames to discard as burn-in when --burn-in-frames is not provided",
    )
    parser.add_argument(
        "--burn-in-frames",
        type=int,
        default=None,
        help="Explicit number of frames to discard as burn-in",
    )
    parser.add_argument("--target-apl", type=float, default=None, help="Optional target APL in Angstrom^2")
    parser.add_argument(
        "--target-thickness",
        type=float,
        default=None,
        help="Optional target PO4 leaflet thickness in Angstrom",
    )
    parser.add_argument(
        "--target-order",
        type=float,
        default=None,
        help="Optional target tail orientational order parameter",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary, time_series = analyze_stage7_packing(
        args.stage_file,
        burn_in_fraction=float(args.burn_in_fraction),
        burn_in_frames=args.burn_in_frames,
        target_apl=args.target_apl,
        target_thickness=args.target_thickness,
        target_order=args.target_order,
    )

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    if args.csv_out is not None:
        args.csv_out.parent.mkdir(parents=True, exist_ok=True)
        write_timeseries_csv(args.csv_out, time_series)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
