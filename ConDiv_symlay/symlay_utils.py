#!/usr/bin/env python3
"""Shared helpers for the constrained symmetric-layer ConDiv workflow."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import tables as tb


def _locate_project_root() -> Path:
    env_root = os.environ.get("CONDIV_PROJECT_ROOT")
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if (root / "py" / "run_upside.py").exists():
            return root

    starts = [Path(__file__).resolve(), Path.cwd().resolve()]
    seen: set[Path] = set()
    for start in starts:
        for candidate in [start.parent, *start.parents]:
            if candidate in seen:
                continue
            seen.add(candidate)
            if (candidate / "py" / "run_upside.py").exists() and (candidate / "parameters").exists():
                return candidate
    raise RuntimeError(
        "Could not locate project root containing py/run_upside.py. "
        "Set CONDIV_PROJECT_ROOT to /Users/yinhan/Documents/upside2-md."
    )


PROJECT_ROOT = _locate_project_root()
PY_DIR = PROJECT_ROOT / "py"
OBJ_DIR = PROJECT_ROOT / "obj"
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))
if str(OBJ_DIR) not in sys.path:
    sys.path.insert(0, str(OBJ_DIR))

import upside_engine as ue


LAYER_TEMPLATE_SCHEMA = "condiv_symlay_layer_template_v1"
LAYER_MANIFEST_SCHEMA = "condiv_symlay_layer_manifest_v1"
DEFAULT_SUPPORT_MARGIN = 0.5
DEFAULT_DENSE_GRID_SIZE = 801


def repo_relative(path: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(resolved)


def load_json(path: Path) -> dict:
    return json.loads(path.expanduser().resolve().read_text(encoding="utf-8"))


def load_layer_template(path: Path) -> dict:
    data = load_json(path)
    if data.get("schema", "") != LAYER_TEMPLATE_SCHEMA:
        raise ValueError(f"Unexpected layer template schema in {path}: {data.get('schema', '')!r}")
    half_slots = data.get("half_slots", [])
    half_type_sequence = data.get("half_type_sequence", [])
    full_type_sequence = data.get("full_type_sequence", [])
    if len(half_slots) != len(half_type_sequence):
        raise ValueError("half_slots length must match half_type_sequence length")
    slot_types = [str(slot["type"]).strip() for slot in half_slots]
    if slot_types != [str(x).strip() for x in half_type_sequence]:
        raise ValueError("half_type_sequence must match half_slots types exactly")
    mirrored = slot_types + slot_types[::-1]
    if [str(x).strip() for x in full_type_sequence] != mirrored:
        raise ValueError("full_type_sequence must equal half_type_sequence mirrored")
    return data


def parse_first_molecule_atom_types(itp_path: Path, molecule_name: str) -> Dict[str, str]:
    lines = itp_path.expanduser().resolve().read_text(encoding="utf-8", errors="ignore").splitlines()
    current_section = ""
    current_molecule = None
    atom_to_type: Dict[str, str] = {}
    target = molecule_name.strip().upper()
    capture_atoms = False

    for raw in lines:
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[]").strip().lower()
            if capture_atoms and current_section == "atoms" and section != "atoms":
                break
            current_section = section
            if section == "moleculetype":
                current_molecule = None
                capture_atoms = False
            continue

        if current_section == "moleculetype":
            if current_molecule is None:
                current_molecule = line.split()[0].upper()
                capture_atoms = current_molecule == target and not atom_to_type
            continue

        if current_section != "atoms" or not capture_atoms:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        atom_to_type[parts[4].strip().upper()] = parts[1].strip()

    if not atom_to_type:
        raise ValueError(f"Could not parse first {molecule_name} atom block from {itp_path}")
    return atom_to_type


def parse_bilayer_signed_depths(pdb_path: Path, allowed_resnames: Iterable[str]) -> Tuple[float, Dict[str, np.ndarray]]:
    allowed = {str(x).strip().upper() for x in allowed_resnames}
    bead_z: Dict[str, List[float]] = {}
    all_lipid_z: List[float] = []
    with pdb_path.expanduser().resolve().open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:21].strip().upper()
            if resname not in allowed:
                continue
            atom_name = line[12:16].strip().upper()
            try:
                z = float(line[46:54])
            except ValueError:
                continue
            bead_z.setdefault(atom_name, []).append(z)
            all_lipid_z.append(z)
    if not all_lipid_z:
        raise ValueError(f"No lipid atoms found in {pdb_path} for resnames {sorted(allowed)}")
    z_center = float(np.mean(np.asarray(all_lipid_z, dtype=np.float64)))
    signed = {
        bead_name: np.asarray(values, dtype=np.float64) - z_center for bead_name, values in bead_z.items()
    }
    return z_center, signed


def build_layer_manifest(template: dict, lipid_itp: Path, bilayer_pdb: Path) -> dict:
    atom_to_type = parse_first_molecule_atom_types(lipid_itp, template["lipid_resname"])
    z_center, bead_signed = parse_bilayer_signed_depths(bilayer_pdb, template["bilayer_resnames"])

    half_slots = []
    for slot_index, slot in enumerate(template["half_slots"]):
        slot_id = str(slot["slot_id"]).strip()
        slot_type = str(slot["type"]).strip()
        members = [str(x).strip().upper() for x in slot["members"]]
        if not members:
            raise ValueError(f"Slot {slot_id} has no members")
        expected_types = {atom_to_type.get(name, "") for name in members}
        if expected_types != {slot_type}:
            raise ValueError(
                f"Slot {slot_id} expected type {slot_type} but parsed member types were {sorted(expected_types)}"
            )
        missing = [name for name in members if name not in bead_signed]
        if missing:
            raise ValueError(f"Slot {slot_id} is missing bead samples for {missing}")
        signed = np.concatenate([np.asarray(bead_signed[name], dtype=np.float64) for name in members], axis=0)
        upper = signed[signed > 0.0]
        lower = signed[signed < 0.0]
        if upper.size == 0 or lower.size == 0:
            raise ValueError(
                f"Slot {slot_id} requires both leaflets; got upper={upper.size} lower={lower.size}"
            )
        half_slots.append(
            {
                "slot_index": int(slot_index),
                "slot_id": slot_id,
                "type": slot_type,
                "members": members,
                "count_total": int(signed.size),
                "count_upper": int(upper.size),
                "count_lower": int(lower.size),
                "mean_abs_depth": float(np.mean(np.abs(signed))),
                "std_abs_depth": float(np.std(np.abs(signed))),
                "upper_mean_signed_depth": float(np.mean(upper)),
                "lower_mean_signed_depth": float(np.mean(lower)),
            }
        )

    sorted_half = sorted(half_slots, key=lambda row: row["mean_abs_depth"])
    positive_projection_depths = [float(row["mean_abs_depth"]) for row in sorted_half]
    positive_projection_slot_ids = [str(row["slot_id"]) for row in sorted_half]

    return {
        "schema": LAYER_MANIFEST_SCHEMA,
        "template_schema": template["schema"],
        "source": {
            "template": repo_relative(Path(template.get("_source_path", "layer_template.json"))),
            "lipid_itp": repo_relative(lipid_itp),
            "bilayer_pdb": repo_relative(bilayer_pdb),
        },
        "lipid_resname": str(template["lipid_resname"]).strip(),
        "bilayer_resnames": [str(x).strip().upper() for x in template["bilayer_resnames"]],
        "z_center": float(z_center),
        "half_type_sequence": [str(x).strip() for x in template["half_type_sequence"]],
        "full_type_sequence": [str(x).strip() for x in template["full_type_sequence"]],
        "half_slots": half_slots,
        "positive_projection_depths": positive_projection_depths,
        "positive_projection_slot_ids": positive_projection_slot_ids,
    }


def load_layer_manifest(path: Path) -> dict:
    data = load_json(path)
    if data.get("schema", "") != LAYER_MANIFEST_SCHEMA:
        raise ValueError(f"Unexpected layer manifest schema in {path}: {data.get('schema', '')!r}")
    if "positive_projection_depths" not in data or not data["positive_projection_depths"]:
        raise ValueError("Layer manifest missing positive_projection_depths")
    return data


def write_layer_manifest_csv(path: Path, manifest: dict) -> None:
    import csv

    fieldnames = [
        "slot_index",
        "slot_id",
        "type",
        "members",
        "count_total",
        "count_upper",
        "count_lower",
        "mean_abs_depth",
        "std_abs_depth",
        "upper_mean_signed_depth",
        "lower_mean_signed_depth",
    ]
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in manifest["half_slots"]:
            out = dict(row)
            out["members"] = ",".join(row["members"])
            writer.writerow(out)


def write_layer_manifest_plot(path: Path, manifest: dict) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = path.parent / ".mpl-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    half_slots = manifest["half_slots"]
    y = np.arange(len(half_slots), dtype=np.float64)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), constrained_layout=True)

    ax = axes[0]
    ax.axvline(0.0, color="0.5", linestyle="--", linewidth=1.0)
    for idx, row in enumerate(half_slots):
        upper = float(row["upper_mean_signed_depth"])
        lower = float(row["lower_mean_signed_depth"])
        spread = float(row["std_abs_depth"])
        label = f"{row['slot_id']} ({row['type']})"
        ax.hlines(y[idx], lower, upper, color="#4c78a8", alpha=0.35, linewidth=2.0)
        ax.errorbar([upper], [y[idx]], xerr=[spread], fmt="o", color="#4c78a8", capsize=3)
        ax.errorbar([lower], [y[idx]], xerr=[spread], fmt="s", color="#f58518", capsize=3)
        ax.text(upper + 0.3, y[idx] + 0.05, label, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([row["type"] for row in half_slots])
    ax.set_xlabel("Signed depth z - z_center (A)")
    ax.set_title("Measured slot depths by leaflet")
    ax.grid(alpha=0.2, linewidth=0.6)

    ax = axes[1]
    full_types = manifest["full_type_sequence"]
    full_depths = manifest["positive_projection_depths"][::-1] + manifest["positive_projection_depths"]
    full_signed = [float(x) for x in full_depths[: len(full_types)]]
    full_signed = [abs(depth) for depth in full_signed[: len(manifest["positive_projection_depths"])]]
    mirrored = [-depth for depth in full_signed[::-1]] + full_signed
    ax.axhline(0.0, color="0.5", linestyle="--", linewidth=1.0)
    ax.scatter(np.arange(len(mirrored), dtype=np.float64), mirrored, c="#54a24b", s=35)
    ax.set_xticks(np.arange(len(full_types), dtype=np.float64))
    ax.set_xticklabels(full_types, rotation=45, ha="right")
    ax.set_ylabel("Representative signed depth (A)")
    ax.set_title("Canonical mirrored type sequence")
    ax.grid(alpha=0.2, linewidth=0.6)

    fig.suptitle("ConDiv_symlay DOPC topology-slot manifest", fontsize=13)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def load_membrane_arrays(path: Path) -> dict:
    with tb.open_file(str(path.expanduser().resolve()), "r") as t:
        root = t.root
        return {
            "cb": root.cb_energy[:].astype(np.float64),
            "icb": root.icb_energy[:].astype(np.float64),
            "hb": root.hb_energy[:].astype(np.float64),
            "ihb": root.ihb_energy[:].astype(np.float64),
            "cb_z_min": float(root._v_attrs.cb_z_min),
            "cb_z_max": float(root._v_attrs.cb_z_max),
            "hb_z_min": float(root._v_attrs.hb_z_min),
            "hb_z_max": float(root._v_attrs.hb_z_max),
        }


def z_to_spline_x(z_values: np.ndarray, z_min: float, z_max: float, n_node: int) -> np.ndarray:
    if n_node < 3:
        raise ValueError("Need at least 3 spline coefficients")
    span = float(z_max) - float(z_min)
    if span <= 0.0:
        raise ValueError(f"Invalid support span: z_min={z_min}, z_max={z_max}")
    x = (np.asarray(z_values, dtype=np.float64) - float(z_min)) * float(n_node - 2) / span
    return np.clip(x, 0.0, float(n_node - 2))


def spline_internal_z(z_min: float, z_max: float, n_node: int) -> np.ndarray:
    x = np.arange(1, n_node - 1, dtype=np.float64)
    return float(z_min) + x * (float(z_max) - float(z_min)) / float(n_node - 2)


def sample_spline_curve(coeff: np.ndarray, z_min: float, z_max: float, z_values: np.ndarray) -> np.ndarray:
    x = z_to_spline_x(z_values, z_min, z_max, int(coeff.size))
    return np.asarray(ue.clamped_spline_value(np.asarray(coeff, dtype=np.float64), x), dtype=np.float64)


def symmetric_support(original_min: float, original_max: float, slot_depths: Sequence[float], margin: float) -> Tuple[float, float]:
    max_slot = max(float(x) for x in slot_depths)
    support_max = float(math.ceil(max(abs(float(original_min)), abs(float(original_max)), max_slot + float(margin))))
    return -support_max, support_max


def _anchor_voronoi_means(
    positive_depths: np.ndarray,
    positive_values: np.ndarray,
    anchors: np.ndarray,
) -> np.ndarray:
    boundaries = np.empty(anchors.size + 1, dtype=np.float64)
    boundaries[0] = 0.0
    boundaries[-1] = anchors[-1]
    if anchors.size > 1:
        boundaries[1:-1] = 0.5 * (anchors[:-1] + anchors[1:])
    anchor_values = np.zeros_like(anchors)
    for idx in range(anchors.size):
        left = boundaries[idx]
        right = boundaries[idx + 1]
        if idx == anchors.size - 1:
            mask = (positive_depths >= left) & (positive_depths <= right)
        else:
            mask = (positive_depths >= left) & (positive_depths < right)
        if not np.any(mask):
            nearest = int(np.argmin(np.abs(positive_depths - anchors[idx])))
            anchor_values[idx] = float(positive_values[nearest])
        else:
            anchor_values[idx] = float(np.mean(positive_values[mask]))
    return anchor_values


def project_curve_to_symmetric_layers(
    coeff: np.ndarray,
    source_support: Tuple[float, float],
    target_support: Tuple[float, float],
    positive_slot_depths: Sequence[float],
    dense_grid_size: int = DEFAULT_DENSE_GRID_SIZE,
) -> Tuple[np.ndarray, dict]:
    coeff = np.asarray(coeff, dtype=np.float64)
    if dense_grid_size < 101:
        raise ValueError("dense_grid_size must be at least 101")
    src_min, src_max = (float(source_support[0]), float(source_support[1]))
    dst_min, dst_max = (float(target_support[0]), float(target_support[1]))
    if abs(dst_min + dst_max) > 1e-6:
        raise ValueError(f"Target support must be symmetric; got ({dst_min}, {dst_max})")

    dense_grid = np.linspace(dst_min, dst_max, int(dense_grid_size), dtype=np.float64)
    dense_pos = dense_grid[dense_grid >= 0.0]
    sampled_pos = sample_spline_curve(coeff, src_min, src_max, dense_pos)
    sampled_neg = sample_spline_curve(coeff, src_min, src_max, -dense_pos)
    even_pos = 0.5 * (sampled_pos + sampled_neg)

    slot_depths = sorted(
        {
            float(depth)
            for depth in positive_slot_depths
            if float(depth) > 0.0 and float(depth) < dst_max - 1e-9
        }
    )
    anchors = np.asarray([0.0, *slot_depths, dst_max], dtype=np.float64)
    anchor_values = _anchor_voronoi_means(dense_pos, even_pos, anchors)

    recon_pos = np.interp(dense_pos, anchors, anchor_values)
    recon_full = np.interp(np.abs(dense_grid), anchors, anchor_values)
    internal_z = spline_internal_z(dst_min, dst_max, int(coeff.size))
    target_internal = np.interp(np.abs(internal_z), anchors, anchor_values)
    new_coeff = np.asarray(ue.clamped_spline_solve(target_internal), dtype=np.float64)
    solved_dense = sample_spline_curve(new_coeff, dst_min, dst_max, dense_grid)

    residual = solved_dense - recon_full
    return new_coeff, {
        "dense_grid_size": int(dense_grid_size),
        "n_anchor": int(anchors.size),
        "slot_depths": slot_depths,
        "anchor_depths": anchors.tolist(),
        "anchor_values": anchor_values.tolist(),
        "projection_residual_max_abs": float(np.max(np.abs(residual))),
        "projection_residual_rms": float(np.sqrt(np.mean(np.square(residual)))),
        "symmetry_residual_max_abs": float(
            np.max(np.abs(solved_dense - solved_dense[::-1]))
        ),
    }


def project_membrane_arrays(
    arrays: dict,
    manifest: dict,
    dense_grid_size: int = DEFAULT_DENSE_GRID_SIZE,
    support_margin: float = DEFAULT_SUPPORT_MARGIN,
) -> Tuple[dict, dict]:
    slot_depths = [float(x) for x in manifest["positive_projection_depths"]]
    cb_support = symmetric_support(arrays["cb_z_min"], arrays["cb_z_max"], slot_depths, support_margin)
    hb_support = symmetric_support(arrays["hb_z_min"], arrays["hb_z_max"], slot_depths, support_margin)
    projected = {
        "cb": np.zeros_like(arrays["cb"], dtype=np.float64),
        "icb": np.zeros_like(arrays["icb"], dtype=np.float64),
        "hb": np.zeros_like(arrays["hb"], dtype=np.float64),
        "ihb": np.zeros_like(arrays["ihb"], dtype=np.float64),
        "cb_z_min": float(cb_support[0]),
        "cb_z_max": float(cb_support[1]),
        "hb_z_min": float(hb_support[0]),
        "hb_z_max": float(hb_support[1]),
    }
    summary = {"cb": [], "icb": [], "hb": [], "ihb": []}

    for i in range(arrays["cb"].shape[0]):
        for j in range(arrays["cb"].shape[1]):
            coeff, stats = project_curve_to_symmetric_layers(
                arrays["cb"][i, j],
                (arrays["cb_z_min"], arrays["cb_z_max"]),
                cb_support,
                slot_depths,
                dense_grid_size=dense_grid_size,
            )
            projected["cb"][i, j] = coeff
            summary["cb"].append(stats)

    for i in range(arrays["icb"].shape[0]):
        coeff, stats = project_curve_to_symmetric_layers(
            arrays["icb"][i],
            (arrays["cb_z_min"], arrays["cb_z_max"]),
            cb_support,
            slot_depths,
            dense_grid_size=dense_grid_size,
        )
        projected["icb"][i] = coeff
        summary["icb"].append(stats)

    for i in range(arrays["hb"].shape[0]):
        for j in range(arrays["hb"].shape[1]):
            coeff, stats = project_curve_to_symmetric_layers(
                arrays["hb"][i, j],
                (arrays["hb_z_min"], arrays["hb_z_max"]),
                hb_support,
                slot_depths,
                dense_grid_size=dense_grid_size,
            )
            projected["hb"][i, j] = coeff
            summary["hb"].append(stats)

    for i in range(arrays["ihb"].shape[0]):
        for j in range(arrays["ihb"].shape[1]):
            coeff, stats = project_curve_to_symmetric_layers(
                arrays["ihb"][i, j],
                (arrays["hb_z_min"], arrays["hb_z_max"]),
                hb_support,
                slot_depths,
                dense_grid_size=dense_grid_size,
            )
            projected["ihb"][i, j] = coeff
            summary["ihb"].append(stats)

    return projected, summary


def write_projected_membrane_file(source_path: Path, dest_path: Path, projected: dict) -> None:
    dest_path = dest_path.expanduser().resolve()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.expanduser().resolve() != dest_path:
        import shutil

        shutil.copyfile(str(source_path.expanduser().resolve()), str(dest_path))
    with tb.open_file(str(dest_path), "a") as t:
        t.root.cb_energy[:] = np.asarray(projected["cb"], dtype=np.float32)
        t.root.icb_energy[:] = np.asarray(projected["icb"], dtype=np.float32)
        t.root.hb_energy[:] = np.asarray(projected["hb"], dtype=np.float32)
        t.root.ihb_energy[:] = np.asarray(projected["ihb"], dtype=np.float32)
        t.root._v_attrs.cb_z_min = float(projected["cb_z_min"])
        t.root._v_attrs.cb_z_max = float(projected["cb_z_max"])
        t.root._v_attrs.hb_z_min = float(projected["hb_z_min"])
        t.root._v_attrs.hb_z_max = float(projected["hb_z_max"])


def _summarize_stat_rows(rows: Sequence[dict]) -> dict:
    max_proj = max((float(row["projection_residual_max_abs"]) for row in rows), default=0.0)
    max_sym = max((float(row["symmetry_residual_max_abs"]) for row in rows), default=0.0)
    rms_vals = np.asarray([float(row["projection_residual_rms"]) for row in rows], dtype=np.float64)
    return {
        "n_curve": int(len(rows)),
        "max_projection_residual": float(max_proj),
        "median_projection_rms": float(np.median(rms_vals)) if rms_vals.size else 0.0,
        "max_symmetry_residual": float(max_sym),
    }


def analyze_membrane_constraints(
    membrane_path: Path,
    manifest: dict,
    dense_grid_size: int = DEFAULT_DENSE_GRID_SIZE,
    support_margin: float = DEFAULT_SUPPORT_MARGIN,
) -> dict:
    arrays = load_membrane_arrays(membrane_path)
    projected, summary = project_membrane_arrays(
        arrays,
        manifest,
        dense_grid_size=dense_grid_size,
        support_margin=support_margin,
    )
    support_checks = {
        "cb_is_symmetric": bool(abs(arrays["cb_z_min"] + arrays["cb_z_max"]) <= 1e-6),
        "hb_is_symmetric": bool(abs(arrays["hb_z_min"] + arrays["hb_z_max"]) <= 1e-6),
        "cb_covers_slots": bool(arrays["cb_z_max"] >= max(manifest["positive_projection_depths"])),
        "hb_covers_slots": bool(arrays["hb_z_max"] >= max(manifest["positive_projection_depths"])),
    }
    support_checks["cb_matches_required_support"] = bool(
        abs(arrays["cb_z_min"] - projected["cb_z_min"]) <= 1e-6
        and abs(arrays["cb_z_max"] - projected["cb_z_max"]) <= 1e-6
    )
    support_checks["hb_matches_required_support"] = bool(
        abs(arrays["hb_z_min"] - projected["hb_z_min"]) <= 1e-6
        and abs(arrays["hb_z_max"] - projected["hb_z_max"]) <= 1e-6
    )
    return {
        "schema": "condiv_symlay_constraint_report_v1",
        "membrane_path": repo_relative(membrane_path),
        "support": {
            "cb": [float(arrays["cb_z_min"]), float(arrays["cb_z_max"])],
            "hb": [float(arrays["hb_z_min"]), float(arrays["hb_z_max"])],
            "required_cb": [float(projected["cb_z_min"]), float(projected["cb_z_max"])],
            "required_hb": [float(projected["hb_z_min"]), float(projected["hb_z_max"])],
        },
        "support_checks": support_checks,
        "projection_summary": {
            "cb": _summarize_stat_rows(summary["cb"]),
            "icb": _summarize_stat_rows(summary["icb"]),
            "hb": _summarize_stat_rows(summary["hb"]),
            "ihb": _summarize_stat_rows(summary["ihb"]),
        },
    }
