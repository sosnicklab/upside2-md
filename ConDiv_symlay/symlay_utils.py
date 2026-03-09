#!/usr/bin/env python3
"""Shared helpers for the constrained symmetric-layer ConDiv workflow."""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from functools import lru_cache
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
LAYER_MANIFEST_SCHEMA = "condiv_symlay_layer_manifest_v2"
LAYER_MANIFEST_SCHEMA_V1 = "condiv_symlay_layer_manifest_v1"
DEFAULT_SUPPORT_MARGIN = 0.5
DEFAULT_DENSE_GRID_SIZE = 801
Q0_SURFACE_OFFSET_ANGSTROM = 4.7
PAIR_TEACHER_SCHEMA = "condiv_symlay_pair_teacher_v2"
PAIR_SCORE_MODE = "explicit_shell_projected_radial_pair_table_v1"
PAIR_ROLE_PROXIES = {"N": "Qd", "CA": "C1", "C": "C2", "O": "Qa"}
PAIR_PROXY_FOR_TARGET = {
    "cb_mean": ("C1", "C2"),
    "hb_donor": ("Qd",),
    "hb_acceptor": ("Qa",),
}
PAIR_SLOT_WEIGHT_MODEL = "simplified_q0_qa_1_else_2_v1"
PAIR_COVERED_TYPES = ("Q0", "Qa", "Na", "C1", "C3")
PAIR_DESCRIPTOR_NAME = "role_proxy_pair_epsilon_v1"
PAIR_RMIN_NAME = "role_proxy_pair_sigma_v1"
PAIR_RADIAL_MIN = 2.0
PAIR_RADIAL_MAX = 12.0
PAIR_RADIAL_STEP = 0.1
PAIR_RADIAL_SHAPE_MODEL = "descriptor_interpolated_radial_table_v1"
PAIR_KERNEL_MODEL = "bilayer_shell_histogram_v1"


def repo_relative(path: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(resolved)


def load_json(path: Path) -> dict:
    return json.loads(path.expanduser().resolve().read_text(encoding="utf-8"))


def _default_slot_weight(slot_type: str) -> float:
    return 1.0 if str(slot_type).strip().upper() in {"Q0", "QA"} else 2.0


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
                "actual_member_count": int(len(members)),
                "slot_weight": float(_default_slot_weight(slot_type)),
                "count_total": int(signed.size),
                "count_upper": int(upper.size),
                "count_lower": int(lower.size),
                "mean_abs_depth": float(np.mean(np.abs(signed))),
                "std_abs_depth": float(np.std(np.abs(signed))),
                "upper_mean_signed_depth": float(np.mean(upper)),
                "lower_mean_signed_depth": float(np.mean(lower)),
            }
        )

    q0_slots = [row for row in half_slots if row["type"] == "Q0"]
    if len(q0_slots) != 1:
        raise ValueError(f"Expected exactly one Q0 slot in half bilayer, found {len(q0_slots)}")
    q0_slot = q0_slots[0]
    upper_surface_z = float(q0_slot["upper_mean_signed_depth"]) + Q0_SURFACE_OFFSET_ANGSTROM
    lower_surface_z = float(q0_slot["lower_mean_signed_depth"]) - Q0_SURFACE_OFFSET_ANGSTROM

    for row in half_slots:
        upper_surface_depth = float(upper_surface_z - float(row["upper_mean_signed_depth"]))
        lower_surface_depth = float(float(row["lower_mean_signed_depth"]) - lower_surface_z)
        row["upper_surface_depth"] = upper_surface_depth
        row["lower_surface_depth"] = lower_surface_depth
        row["mean_surface_depth"] = 0.5 * (upper_surface_depth + lower_surface_depth)

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
        "surface_reference": {
            "definition": "Q0 signed z +/- 4.7A",
            "q0_surface_offset_angstrom": float(Q0_SURFACE_OFFSET_ANGSTROM),
            "upper_surface_z": upper_surface_z,
            "lower_surface_z": lower_surface_z,
        },
        "topology_weight_model": PAIR_SLOT_WEIGHT_MODEL,
        "half_slots": half_slots,
        "positive_projection_depths": positive_projection_depths,
        "positive_projection_slot_ids": positive_projection_slot_ids,
    }


def load_layer_manifest(path: Path) -> dict:
    data = load_json(path)
    schema = data.get("schema", "")
    if schema not in {LAYER_MANIFEST_SCHEMA, LAYER_MANIFEST_SCHEMA_V1}:
        raise ValueError(f"Unexpected layer manifest schema in {path}: {data.get('schema', '')!r}")
    if "positive_projection_depths" not in data or not data["positive_projection_depths"]:
        raise ValueError("Layer manifest missing positive_projection_depths")
    if "surface_reference" not in data:
        q0_slots = [row for row in data["half_slots"] if row["type"] == "Q0"]
        if len(q0_slots) != 1:
            raise ValueError(f"{path}: cannot infer surface reference without exactly one Q0 slot")
        q0_slot = q0_slots[0]
        data["surface_reference"] = {
            "definition": "Q0 signed z +/- 4.7A",
            "q0_surface_offset_angstrom": float(Q0_SURFACE_OFFSET_ANGSTROM),
            "upper_surface_z": float(q0_slot["upper_mean_signed_depth"]) + Q0_SURFACE_OFFSET_ANGSTROM,
            "lower_surface_z": float(q0_slot["lower_mean_signed_depth"]) - Q0_SURFACE_OFFSET_ANGSTROM,
        }
    for row in data["half_slots"]:
        row.setdefault("actual_member_count", int(len(row.get("members", []))))
        row.setdefault("slot_weight", float(_default_slot_weight(row["type"])))
        if "mean_surface_depth" not in row:
            upper_surface_z = float(data["surface_reference"]["upper_surface_z"])
            lower_surface_z = float(data["surface_reference"]["lower_surface_z"])
            upper_surface_depth = float(upper_surface_z - float(row["upper_mean_signed_depth"]))
            lower_surface_depth = float(float(row["lower_mean_signed_depth"]) - lower_surface_z)
            row["upper_surface_depth"] = upper_surface_depth
            row["lower_surface_depth"] = lower_surface_depth
            row["mean_surface_depth"] = 0.5 * (upper_surface_depth + lower_surface_depth)
    data.setdefault("topology_weight_model", PAIR_SLOT_WEIGHT_MODEL)
    return data


def write_layer_manifest_csv(path: Path, manifest: dict) -> None:
    import csv

    fieldnames = [
        "slot_index",
        "slot_id",
        "type",
        "members",
        "actual_member_count",
        "slot_weight",
        "count_total",
        "count_upper",
        "count_lower",
        "mean_abs_depth",
        "mean_surface_depth",
        "std_abs_depth",
        "upper_mean_signed_depth",
        "lower_mean_signed_depth",
        "upper_surface_depth",
        "lower_surface_depth",
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


_SPLINE_INFLUENCE_CACHE: dict[tuple, np.ndarray] = {}


def spline_influence_matrix(sample_z: Sequence[float], z_min: float, z_max: float, n_node: int) -> np.ndarray:
    rounded_depths = tuple(round(float(x), 6) for x in sample_z)
    key = (rounded_depths, round(float(z_min), 6), round(float(z_max), 6), int(n_node))
    cached = _SPLINE_INFLUENCE_CACHE.get(key)
    if cached is not None:
        return cached

    z = np.asarray(sample_z, dtype=np.float64)
    x = z_to_spline_x(z, z_min, z_max, int(n_node))
    basis = np.zeros((z.shape[0], int(n_node)), dtype=np.float64)
    eye = np.eye(int(n_node), dtype=np.float64)
    for idx in range(int(n_node)):
        basis[:, idx] = np.asarray(ue.clamped_spline_value(eye[idx], x), dtype=np.float64)
    _SPLINE_INFLUENCE_CACHE[key] = basis
    return basis


def weighted_affine_trend_residual(values: np.ndarray, scores: np.ndarray, weights: np.ndarray) -> dict:
    y = np.asarray(values, dtype=np.float64)
    s = np.asarray(scores, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    x = np.column_stack([np.ones_like(s), s])
    wx = w[:, None] * x
    gram = x.T @ wx
    rhs = x.T @ (w * y)
    beta = np.linalg.solve(gram + 1e-8 * np.eye(gram.shape[0], dtype=np.float64), rhs)
    fitted = x @ beta
    residual = y - fitted
    grad_values = w * residual
    return {
        "loss": 0.5 * float(np.sum(w * residual * residual)),
        "beta_intercept": float(beta[0]),
        "beta_scale": float(beta[1]),
        "fitted": fitted,
        "residual": residual,
        "grad_values": grad_values,
        "weighted_rms_residual": float(np.sqrt(np.sum(w * residual * residual) / max(np.sum(w), 1e-8))),
        "max_abs_residual": float(np.max(np.abs(residual))) if residual.size else 0.0,
    }


@lru_cache(maxsize=1)
def _load_martini_lib_module():
    module_path = PROJECT_ROOT / "example" / "16.MARTINI" / "lib.py"
    spec = importlib.util.spec_from_file_location("upside_martini_lib", str(module_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load MARTINI helper module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_dry_martini_param_table(ff_itp: Path):
    module = _load_martini_lib_module()
    table = module.read_martini3_nonbond_params(str(ff_itp.expanduser().resolve()))
    if not table:
        raise ValueError(f"Failed to read dry MARTINI pair parameters from {ff_itp}")
    return table


def pair_param(param_table, type1: str, type2: str) -> Tuple[float, float]:
    key1 = (str(type1).strip(), str(type2).strip())
    key2 = (key1[1], key1[0])
    if key1 in param_table:
        sigma, epsilon = param_table[key1]
        return float(sigma), float(epsilon)
    if key2 in param_table:
        sigma, epsilon = param_table[key2]
        return float(sigma), float(epsilon)
    raise KeyError(f"Missing dry MARTINI pair parameters for ({type1}, {type2})")


def pair_role_descriptor_and_sigma(param_table, proxies: Sequence[str], env_type: str) -> Tuple[float, float]:
    sigmas = []
    descriptors = []
    for proxy in proxies:
        sigma_nm, epsilon = pair_param(param_table, proxy, env_type)
        sigmas.append(float(sigma_nm))
        descriptors.append(float(epsilon))
    return float(np.mean(np.asarray(sigmas, dtype=np.float64))), float(
        np.mean(np.asarray(descriptors, dtype=np.float64))
    )


def _resolve_manifest_source_path(raw_path: str) -> Path:
    path = Path(str(raw_path)).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (PROJECT_ROOT / path).resolve()


def _build_slot_pooling(manifest: dict) -> dict:
    half_slots = manifest["half_slots"]
    slot_ids = [str(row["slot_id"]).strip() for row in half_slots]
    slot_types = [str(row["type"]).strip() for row in half_slots]
    slot_depths = np.asarray([float(row["mean_abs_depth"]) for row in half_slots], dtype=np.float64)
    slot_surface_depths = np.asarray([float(row["mean_surface_depth"]) for row in half_slots], dtype=np.float64)
    covered_types = list(PAIR_COVERED_TYPES)
    slot_pooling = np.zeros((len(covered_types), len(half_slots)), dtype=np.float64)
    covered_weights = np.zeros(len(covered_types), dtype=np.float64)

    for type_index, covered_type in enumerate(covered_types):
        for slot_index, row in enumerate(half_slots):
            if str(row["type"]).strip() != covered_type:
                continue
            member_count = int(row.get("actual_member_count", len(row.get("members", [])) or 1))
            weight = float(row.get("slot_weight", _default_slot_weight(row["type"]))) * float(max(member_count, 1))
            slot_pooling[type_index, slot_index] = weight
            covered_weights[type_index] += weight
        if covered_weights[type_index] <= 0.0:
            raise ValueError(f"Layer manifest is missing required covered type {covered_type}")
        slot_pooling[type_index, :] /= covered_weights[type_index]

    return {
        "covered_types": covered_types,
        "slot_ids": slot_ids,
        "slot_types": slot_types,
        "slot_depths": slot_depths,
        "slot_surface_depths": slot_surface_depths,
        "slot_pooling": slot_pooling,
        "covered_weights": covered_weights,
    }


def _sample_reference_role_slot_values(membrane_arrays: dict, slot_depths: np.ndarray) -> dict:
    cb_values = np.zeros((membrane_arrays["cb"].shape[0], slot_depths.shape[0]), dtype=np.float64)
    for aa_index in range(membrane_arrays["cb"].shape[0]):
        burial0 = sample_spline_curve(
            membrane_arrays["cb"][aa_index, 0],
            membrane_arrays["cb_z_min"],
            membrane_arrays["cb_z_max"],
            slot_depths,
        )
        burial1 = sample_spline_curve(
            membrane_arrays["cb"][aa_index, 1],
            membrane_arrays["cb_z_min"],
            membrane_arrays["cb_z_max"],
            slot_depths,
        )
        cb_values[aa_index, :] = 0.5 * (burial0 + burial1)

    donor_unpaired = sample_spline_curve(
        membrane_arrays["hb"][0, 0],
        membrane_arrays["hb_z_min"],
        membrane_arrays["hb_z_max"],
        slot_depths,
    )
    donor_paired = sample_spline_curve(
        membrane_arrays["hb"][0, 1],
        membrane_arrays["hb_z_min"],
        membrane_arrays["hb_z_max"],
        slot_depths,
    )
    acceptor_unpaired = sample_spline_curve(
        membrane_arrays["hb"][1, 0],
        membrane_arrays["hb_z_min"],
        membrane_arrays["hb_z_max"],
        slot_depths,
    )
    acceptor_paired = sample_spline_curve(
        membrane_arrays["hb"][1, 1],
        membrane_arrays["hb_z_min"],
        membrane_arrays["hb_z_max"],
        slot_depths,
    )

    return {
        "cb_mean": np.mean(cb_values, axis=0),
        "hb_donor": 0.5 * (donor_unpaired + donor_paired),
        "hb_acceptor": 0.5 * (acceptor_unpaired + acceptor_paired),
    }


def _parse_bilayer_bead_records(manifest: dict) -> list[dict]:
    source = manifest.get("source", {})
    if "lipid_itp" not in source or "bilayer_pdb" not in source:
        raise ValueError("Layer manifest must record source lipid_itp and bilayer_pdb for pair-teacher construction")
    lipid_itp = _resolve_manifest_source_path(source["lipid_itp"])
    bilayer_pdb = _resolve_manifest_source_path(source["bilayer_pdb"])
    atom_to_type = parse_first_molecule_atom_types(lipid_itp, manifest["lipid_resname"])
    allowed_resnames = {str(x).strip().upper() for x in manifest["bilayer_resnames"]}
    records = []
    with bilayer_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:21].strip().upper()
            if resname not in allowed_resnames:
                continue
            atom_name = line[12:16].strip().upper()
            bead_type = atom_to_type.get(atom_name)
            if bead_type is None:
                continue
            try:
                coord = np.asarray(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    dtype=np.float64,
                )
            except ValueError:
                continue
            records.append({"atom_name": atom_name, "bead_type": str(bead_type).strip(), "coord": coord})
    if not records:
        raise ValueError(f"No bilayer bead records parsed from {bilayer_pdb}")
    return records


def _build_shell_kernel(manifest: dict, env_types: Sequence[str], radial_grid: np.ndarray) -> np.ndarray:
    records = _parse_bilayer_bead_records(manifest)
    coord_by_type: Dict[str, List[np.ndarray]] = {}
    for record in records:
        coord_by_type.setdefault(record["bead_type"], []).append(record["coord"])
    coord_by_type = {
        bead_type: np.asarray(coords, dtype=np.float64) for bead_type, coords in coord_by_type.items() if coords
    }

    if radial_grid.ndim != 1 or radial_grid.size < 2:
        raise ValueError("radial_grid must be rank-1 with at least two samples")
    radial_step = float(radial_grid[1] - radial_grid[0])
    if radial_step <= 0.0:
        raise ValueError("radial_grid must be strictly increasing")
    edges = np.empty(radial_grid.shape[0] + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (radial_grid[:-1] + radial_grid[1:])
    edges[0] = max(0.0, float(radial_grid[0]) - 0.5 * radial_step)
    edges[-1] = float(radial_grid[-1]) + 0.5 * radial_step

    env_index = {str(env_type): idx for idx, env_type in enumerate(env_types)}
    kernel = np.zeros((len(PAIR_COVERED_TYPES), len(env_types), radial_grid.shape[0]), dtype=np.float64)
    for covered_index, covered_type in enumerate(PAIR_COVERED_TYPES):
        centers = coord_by_type.get(covered_type)
        if centers is None or centers.size == 0:
            raise ValueError(f"Bilayer template is missing coordinates for covered type {covered_type}")
        accum = np.zeros((len(env_types), radial_grid.shape[0]), dtype=np.float64)
        for center in centers:
            for env_type, env_coords in coord_by_type.items():
                if env_type not in env_index:
                    continue
                delta = env_coords - center[None, :]
                dist = np.sqrt(np.sum(delta * delta, axis=1))
                if env_type == covered_type:
                    dist = dist[dist > 1e-6]
                hist, _ = np.histogram(dist, bins=edges)
                accum[env_index[env_type], :] += hist.astype(np.float64)
        kernel[covered_index, :, :] = accum / float(centers.shape[0])
    return kernel


def _merge_anchor_points_with_labels(x: np.ndarray, y: np.ndarray, labels: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, list[str]]:
    order = np.argsort(np.asarray(x, dtype=np.float64))
    x_sorted = np.asarray(x, dtype=np.float64)[order]
    y_sorted = np.asarray(y, dtype=np.float64)[order]
    labels_sorted = [str(labels[idx]) for idx in order.tolist()]
    merged_x: List[float] = []
    merged_y: List[float] = []
    merged_labels: List[str] = []
    cursor = 0
    while cursor < x_sorted.shape[0]:
        next_cursor = cursor + 1
        y_values = [float(y_sorted[cursor])]
        label_values = [labels_sorted[cursor]]
        while next_cursor < x_sorted.shape[0] and abs(float(x_sorted[next_cursor]) - float(x_sorted[cursor])) < 1e-8:
            y_values.append(float(y_sorted[next_cursor]))
            label_values.append(labels_sorted[next_cursor])
            next_cursor += 1
        merged_x.append(float(x_sorted[cursor]))
        merged_y.append(float(np.mean(np.asarray(y_values, dtype=np.float64))))
        merged_labels.append("+".join(label_values))
        cursor = next_cursor
    if len(merged_x) < 2:
        raise ValueError("Need at least two unique dry-MARTINI anchor descriptors")
    return (
        np.asarray(merged_x, dtype=np.float64),
        np.asarray(merged_y, dtype=np.float64),
        merged_labels,
    )


def _pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n_point = x.shape[0]
    h = np.diff(x)
    delta = np.diff(y) / h
    slopes = np.zeros(n_point, dtype=np.float64)
    if n_point == 2:
        slopes[:] = delta[0]
        return slopes
    for index in range(1, n_point - 1):
        if delta[index - 1] == 0.0 or delta[index] == 0.0 or np.sign(delta[index - 1]) != np.sign(delta[index]):
            slopes[index] = 0.0
        else:
            w1 = 2.0 * h[index] + h[index - 1]
            w2 = h[index] + 2.0 * h[index - 1]
            slopes[index] = (w1 + w2) / ((w1 / delta[index - 1]) + (w2 / delta[index]))
    slope0 = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    if np.sign(slope0) != np.sign(delta[0]):
        slope0 = 0.0
    elif np.sign(delta[0]) != np.sign(delta[1]) and abs(slope0) > 3.0 * abs(delta[0]):
        slope0 = 3.0 * delta[0]
    slopes[0] = slope0
    slopen = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    if np.sign(slopen) != np.sign(delta[-1]):
        slopen = 0.0
    elif np.sign(delta[-1]) != np.sign(delta[-2]) and abs(slopen) > 3.0 * abs(delta[-1]):
        slopen = 3.0 * delta[-1]
    slopes[-1] = slopen
    return slopes


def _pchip_interpolate(x: np.ndarray, y: np.ndarray, x_query: np.ndarray) -> np.ndarray:
    slopes = _pchip_slopes(x, y)
    interval = np.searchsorted(x, x_query, side="right") - 1
    interval = np.clip(interval, 0, x.shape[0] - 2)
    x0 = x[interval]
    x1 = x[interval + 1]
    y0 = y[interval]
    y1 = y[interval + 1]
    m0 = slopes[interval]
    m1 = slopes[interval + 1]
    h = x1 - x0
    t = (x_query - x0) / h
    h00 = (2.0 * t**3) - (3.0 * t**2) + 1.0
    h10 = (t**3) - (2.0 * t**2) + t
    h01 = (-2.0 * t**3) + (3.0 * t**2)
    h11 = (t**3) - (t**2)
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def _interpolate_anchor_value(descriptor: float, anchor_x: np.ndarray, anchor_y: np.ndarray) -> float:
    descriptor = float(descriptor)
    if descriptor <= float(anchor_x[0]):
        return float(anchor_y[0])
    if descriptor <= float(anchor_x[-1]):
        return float(_pchip_interpolate(anchor_x, anchor_y, np.asarray([descriptor], dtype=np.float64))[0])
    slope = float((anchor_y[-1] - anchor_y[-2]) / (anchor_x[-1] - anchor_x[-2]))
    return float(anchor_y[-1] + slope * (descriptor - anchor_x[-1]))


def _merge_control_points(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x)
    x_sorted = np.asarray(x, dtype=np.float64)[order]
    y_sorted = np.asarray(y, dtype=np.float64)[order]
    merged_x: List[float] = []
    merged_y: List[float] = []
    cursor = 0
    while cursor < x_sorted.shape[0]:
        next_cursor = cursor + 1
        y_values = [float(y_sorted[cursor])]
        while next_cursor < x_sorted.shape[0] and abs(float(x_sorted[next_cursor]) - float(x_sorted[cursor])) < 1e-8:
            y_values.append(float(y_sorted[next_cursor]))
            next_cursor += 1
        merged_x.append(float(x_sorted[cursor]))
        merged_y.append(float(np.mean(np.asarray(y_values, dtype=np.float64))))
        cursor = next_cursor
    if len(merged_x) < 2:
        raise ValueError("Need at least two unique radial control points")
    return np.asarray(merged_x, dtype=np.float64), np.asarray(merged_y, dtype=np.float64)


def _radial_control_points(r_min_angstrom: float, attraction_depth: float) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(
        [
            max(PAIR_RADIAL_MIN, 0.72 * r_min_angstrom),
            max(PAIR_RADIAL_MIN + 0.2, 0.90 * r_min_angstrom),
            r_min_angstrom,
            min(PAIR_RADIAL_MAX, 1.35 * r_min_angstrom),
            PAIR_RADIAL_MAX,
        ],
        dtype=np.float64,
    )
    y = np.asarray(
        [
            3.0 * attraction_depth,
            1.2 * attraction_depth,
            -attraction_depth,
            -0.30 * attraction_depth,
            0.0,
        ],
        dtype=np.float64,
    )
    return _merge_control_points(x, y)


def _build_role_pair_table(
    param_table,
    role_name: str,
    covered_types: Sequence[str],
    covered_values: np.ndarray,
    env_types: Sequence[str],
    radial_grid: np.ndarray,
) -> Tuple[np.ndarray, dict]:
    proxy_types = PAIR_PROXY_FOR_TARGET[role_name]
    anchor_descriptors = []
    anchor_labels = []
    for covered_type in covered_types:
        _sigma_nm, descriptor = pair_role_descriptor_and_sigma(param_table, proxy_types, covered_type)
        anchor_descriptors.append(float(descriptor))
        anchor_labels.append(str(covered_type))
    anchor_x, anchor_y, merged_labels = _merge_anchor_points_with_labels(
        np.asarray(anchor_descriptors, dtype=np.float64),
        np.asarray(covered_values, dtype=np.float64),
        anchor_labels,
    )

    table = np.zeros((len(env_types), radial_grid.shape[0]), dtype=np.float64)
    pair_summary = {}
    for env_index, env_type in enumerate(env_types):
        sigma_nm, descriptor = pair_role_descriptor_and_sigma(param_table, proxy_types, env_type)
        r_min_angstrom = float((2.0 ** (1.0 / 6.0)) * sigma_nm * 10.0)
        attraction_depth = _interpolate_anchor_value(float(descriptor), anchor_x, anchor_y)
        ctrl_x, ctrl_y = _radial_control_points(r_min_angstrom, attraction_depth)
        table[env_index, :] = _pchip_interpolate(ctrl_x, ctrl_y, radial_grid)
        pair_summary[str(env_type)] = {
            "descriptor_value": float(descriptor),
            "sigma_nm": float(sigma_nm),
            "r_min_angstrom": r_min_angstrom,
            "attraction_depth": float(attraction_depth),
            "control_r": [float(x) for x in ctrl_x],
            "control_v": [float(y) for y in ctrl_y],
        }

    return table, {
        "descriptor_name": PAIR_DESCRIPTOR_NAME,
        "rmin_name": PAIR_RMIN_NAME,
        "anchor_descriptor": [float(x) for x in anchor_x],
        "anchor_value": [float(y) for y in anchor_y],
        "anchor_label": merged_labels,
        "pair_summary": pair_summary,
    }


def _project_table_to_covered_scores(table_values: np.ndarray, shell_kernel: np.ndarray) -> np.ndarray:
    if shell_kernel.shape[1:] != table_values.shape:
        raise ValueError(
            f"Shell kernel shape {shell_kernel.shape} is incompatible with table values {table_values.shape}"
        )
    return np.tensordot(shell_kernel, table_values, axes=([1, 2], [0, 1])).astype(np.float64)


def build_pair_teacher(manifest: dict, ff_itp: Path, membrane_file: Path) -> dict:
    param_table = load_dry_martini_param_table(ff_itp)
    membrane_arrays = load_membrane_arrays(Path(membrane_file))
    slot_info = _build_slot_pooling(manifest)
    slot_pooling = slot_info["slot_pooling"]
    reference_slot_values = _sample_reference_role_slot_values(membrane_arrays, slot_info["slot_depths"])

    env_types = sorted({str(atom_type) for pair in param_table.keys() for atom_type in pair})
    radial_grid = np.arange(PAIR_RADIAL_MIN, PAIR_RADIAL_MAX + 0.5 * PAIR_RADIAL_STEP, PAIR_RADIAL_STEP, dtype=np.float64)
    if radial_grid[-1] > PAIR_RADIAL_MAX:
        radial_grid[-1] = PAIR_RADIAL_MAX
    shell_kernel = _build_shell_kernel(manifest, env_types, radial_grid)

    role_tables = {}
    role_calibration = {}
    covered_reference = {}
    covered_scores = {}
    for role_name, slot_values in reference_slot_values.items():
        covered_values = slot_pooling @ np.asarray(slot_values, dtype=np.float64)
        table, calibration = _build_role_pair_table(
            param_table,
            role_name,
            slot_info["covered_types"],
            covered_values,
            env_types,
            radial_grid,
        )
        role_tables[role_name] = table
        role_calibration[role_name] = calibration
        covered_reference[role_name] = covered_values.tolist()
        covered_scores[role_name] = _project_table_to_covered_scores(table, shell_kernel).tolist()

    return {
        "schema": PAIR_TEACHER_SCHEMA,
        "score_mode": PAIR_SCORE_MODE,
        "ff_itp": repo_relative(ff_itp),
        "reference_membrane_file": repo_relative(Path(membrane_file)),
        "slot_weight_model": manifest.get("topology_weight_model", PAIR_SLOT_WEIGHT_MODEL),
        "covered_types": slot_info["covered_types"],
        "covered_weights": slot_info["covered_weights"].tolist(),
        "slot_ids": slot_info["slot_ids"],
        "slot_types": slot_info["slot_types"],
        "slot_depths": slot_info["slot_depths"].tolist(),
        "slot_surface_depths": slot_info["slot_surface_depths"].tolist(),
        "slot_pooling": slot_pooling.tolist(),
        "env_types": env_types,
        "radial_grid": radial_grid.tolist(),
        "radial_shape_model": PAIR_RADIAL_SHAPE_MODEL,
        "kernel_model": PAIR_KERNEL_MODEL,
        "shell_kernel": shell_kernel.tolist(),
        "role_proxies": {key: list(value) for key, value in PAIR_PROXY_FOR_TARGET.items()},
        "role_calibration": role_calibration,
        "reference_covered_values": covered_reference,
        "covered_scores": covered_scores,
        "table_values": {key: value.tolist() for key, value in role_tables.items()},
        "sources": manifest.get("source", {}),
    }


def compute_pair_regularization_gradient(
    param,
    init_membrane_file: str,
    manifest: dict,
    pair_teacher: dict,
    cb_weight: float,
    donor_weight: float,
    acceptor_weight: float,
):
    init_arrays = load_membrane_arrays(Path(init_membrane_file))
    slot_depths = np.asarray(pair_teacher["slot_depths"], dtype=np.float64)
    slot_pooling = np.asarray(pair_teacher["slot_pooling"], dtype=np.float64)
    covered_weights = np.asarray(pair_teacher["covered_weights"], dtype=np.float64)
    cb_scores = np.asarray(pair_teacher["covered_scores"]["cb_mean"], dtype=np.float64)
    donor_scores = np.asarray(pair_teacher["covered_scores"]["hb_donor"], dtype=np.float64)
    acceptor_scores = np.asarray(pair_teacher["covered_scores"]["hb_acceptor"], dtype=np.float64)

    cb_basis = spline_influence_matrix(slot_depths, init_arrays["cb_z_min"], init_arrays["cb_z_max"], param.cb.shape[-1])
    hb_basis = spline_influence_matrix(slot_depths, init_arrays["hb_z_min"], init_arrays["hb_z_max"], param.hb.shape[-1])

    grad = type(param)(
        np.zeros_like(param.cb, dtype=np.float32),
        np.zeros_like(param.icb, dtype=np.float32),
        np.zeros_like(param.hb, dtype=np.float32),
        np.zeros_like(param.ihb, dtype=np.float32),
    )
    stats = {
        "pair_loss_cb": 0.0,
        "pair_loss_hb_donor": 0.0,
        "pair_loss_hb_acceptor": 0.0,
        "pair_loss_total": 0.0,
        "pair_cb_residual_rms_median": 0.0,
        "pair_cb_residual_max": 0.0,
        "pair_hb_donor_residual_rms": 0.0,
        "pair_hb_acceptor_residual_rms": 0.0,
    }

    cb_rms = []
    cb_max = []
    if float(cb_weight) > 0.0:
        for aa in range(param.cb.shape[0]):
            slot_values = 0.5 * (
                cb_basis @ np.asarray(param.cb[aa, 0], dtype=np.float64)
                + cb_basis @ np.asarray(param.cb[aa, 1], dtype=np.float64)
            )
            y = slot_pooling @ slot_values
            trend = weighted_affine_trend_residual(y, cb_scores, covered_weights)
            grad_type = float(cb_weight) * trend["grad_values"]
            grad_slot = slot_pooling.T @ grad_type
            grad_curve = 0.5 * (cb_basis.T @ grad_slot)
            grad.cb[aa, 0] += grad_curve.astype(np.float32)
            grad.cb[aa, 1] += grad_curve.astype(np.float32)
            stats["pair_loss_cb"] += float(cb_weight) * trend["loss"]
            cb_rms.append(float(trend["weighted_rms_residual"]))
            cb_max.append(float(trend["max_abs_residual"]))

    def _accumulate_hb(tp_index: int, score_vector: np.ndarray, field_name: str, field_rms_name: str, weight_value: float):
        if float(weight_value) <= 0.0:
            return
        slot_values = 0.5 * (
            hb_basis @ np.asarray(param.hb[tp_index, 0], dtype=np.float64)
            + hb_basis @ np.asarray(param.hb[tp_index, 1], dtype=np.float64)
        )
        y = slot_pooling @ slot_values
        trend = weighted_affine_trend_residual(y, score_vector, covered_weights)
        grad_type = float(weight_value) * trend["grad_values"]
        grad_slot = slot_pooling.T @ grad_type
        grad_curve = 0.5 * (hb_basis.T @ grad_slot)
        grad.hb[tp_index, 0] += grad_curve.astype(np.float32)
        grad.hb[tp_index, 1] += grad_curve.astype(np.float32)
        stats[field_name] += float(weight_value) * trend["loss"]
        stats[field_rms_name] = float(trend["weighted_rms_residual"])
        stats[f"{field_name}_max_abs_residual"] = float(trend["max_abs_residual"])

    _accumulate_hb(0, donor_scores, "pair_loss_hb_donor", "pair_hb_donor_residual_rms", donor_weight)
    _accumulate_hb(1, acceptor_scores, "pair_loss_hb_acceptor", "pair_hb_acceptor_residual_rms", acceptor_weight)

    stats["pair_loss_total"] = (
        float(stats["pair_loss_cb"]) + float(stats["pair_loss_hb_donor"]) + float(stats["pair_loss_hb_acceptor"])
    )
    if cb_rms:
        stats["pair_cb_residual_rms_median"] = float(np.median(np.asarray(cb_rms, dtype=np.float64)))
        stats["pair_cb_residual_max"] = float(np.max(np.asarray(cb_max, dtype=np.float64)))

    return grad, stats


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
