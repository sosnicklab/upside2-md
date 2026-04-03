#!/usr/bin/env python3
"""SC training workflow for dry-MARTINI sidechain-type tables."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import shlex
import subprocess as sp
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import h5py
import numpy as np

SCHEMA_MANIFEST = "sc_training_manifest_v3"
SCHEMA_TASK = "sc_training_task_result_v3"
SCHEMA_TABLE = "sc_training_table_v3"
SCHEMA_SLURM_ROUND = "sc_training_slurm_round_v1"
SCHEMA_SIDECHAINS = "sc_training_sidechains_v1"

COULOMB_K_DRY_KJ_NM = 138.935458 / 15.0
ANGSTROM_TO_NM = 0.1
CANONICAL_RESIDUES = (
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
)

POSITIVE_TYPES = {"Qda", "Qd", "SQda", "SQd"}
NEGATIVE_TYPES = {"Qa", "SQa"}

CANONICAL_CB_POSITION_ANG = (0.0, 0.94375626, 1.2068012)
_cb_norm = math.sqrt(sum(x * x for x in CANONICAL_CB_POSITION_ANG))
CANONICAL_CB_VECTOR_UNIT = tuple(x / _cb_norm for x in CANONICAL_CB_POSITION_ANG)


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _workflow_dir() -> Path:
    return Path(__file__).resolve().parent


def _data_dir() -> Path:
    return _workflow_dir() / "data"


def _project_root() -> Path:
    return _workflow_dir().parent


def _default_sidechain_library() -> Path:
    return _project_root() / "parameters" / "ff_2.1" / "sidechain.h5"


def _manifest_path(base_dir: Path) -> Path:
    return base_dir / "training_manifest.json"


def _task_result_dir(base_dir: Path) -> Path:
    return base_dir / "results" / "tasks"


def _assembled_dir(base_dir: Path) -> Path:
    return base_dir / "results" / "assembled"


def _slurm_dir(base_dir: Path) -> Path:
    return base_dir / "slurm"


def _benchmark_dir(base_dir: Path) -> Path:
    return base_dir / "benchmark"


def _current_run_record() -> Path:
    return _workflow_dir() / ".sc_training_current_run_dir"


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _write_text(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def _record_current_run_dir(base_dir: Path) -> None:
    _current_run_record().write_text(str(base_dir.resolve()) + "\n", encoding="utf-8")


def _infer_type_charge(bead_type: str) -> float:
    bead_type = bead_type.strip()
    if bead_type in POSITIVE_TYPES:
        return 1.0
    if bead_type in NEGATIVE_TYPES:
        return -1.0
    return 0.0


def _resolve_base_dir(base_dir_arg: str | None) -> Path:
    if base_dir_arg:
        return Path(base_dir_arg).expanduser().resolve()
    record = _current_run_record()
    if record.exists():
        return Path(record.read_text(encoding="utf-8").strip()).expanduser().resolve()
    return (_workflow_dir() / "runs" / "default").resolve()


def _load_martini_forcefield(martinize_path: Path, forcefield_name: str) -> Dict[str, List[str]]:
    if martinize_path.suffix.lower() == ".json":
        payload = _load_json(martinize_path)
        residues = payload.get("residues", payload)
        if not isinstance(residues, dict):
            raise RuntimeError(f"Invalid sidechain JSON payload in {martinize_path}")
        if "forcefield_name" in payload and forcefield_name and payload["forcefield_name"] != forcefield_name:
            raise RuntimeError(
                f"Sidechain JSON forcefield mismatch: requested '{forcefield_name}', "
                f"file provides '{payload['forcefield_name']}'"
            )
        residue_map: Dict[str, List[str]] = {}
        for residue in CANONICAL_RESIDUES:
            raw = residues.get(residue, [])
            residue_map[residue] = [str(tok).strip() for tok in raw if str(tok).strip()]
        return residue_map

    spec = importlib.util.spec_from_file_location("sc_training_martinize_runtime", martinize_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load martinize module from {martinize_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, forcefield_name):
        raise RuntimeError(f"Forcefield '{forcefield_name}' not found in {martinize_path}")
    ff = getattr(module, forcefield_name)()
    residue_map: Dict[str, List[str]] = {}
    for residue in CANONICAL_RESIDUES:
        raw = ff.sidechains.get(residue, [])
        if not raw:
            residue_map[residue] = []
            continue
        bead_tokens = [str(tok).strip() for tok in raw[0] if str(tok).strip()]
        residue_map[residue] = [tok for tok in bead_tokens if tok != "D"]
    return residue_map


def _parse_dry_forcefield(ff_path: Path) -> Tuple[List[str], Dict[Tuple[str, str], Dict[str, float]]]:
    macros: Dict[str, Tuple[float, float]] = {}
    atomtypes: List[str] = []
    pair_params: Dict[Tuple[str, str], Dict[str, float]] = {}
    section = ""

    with ff_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            stripped = raw.split(";", 1)[0].strip()
            if not stripped:
                continue
            if stripped.startswith("#define"):
                parts = stripped.split()
                if len(parts) == 4:
                    try:
                        macros[parts[1]] = (float(parts[2]), float(parts[3]))
                    except ValueError:
                        pass
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                section = stripped[1:-1].strip().lower()
                continue

            parts = stripped.split()
            if section == "atomtypes":
                atomtypes.append(parts[0])
                continue
            if section != "nonbond_params":
                continue
            if len(parts) < 4 or parts[2] != "1":
                continue

            type_i = parts[0]
            type_j = parts[1]
            sigma_nm: float
            epsilon_kj: float
            if len(parts) == 4:
                macro = parts[3]
                if macro not in macros:
                    raise RuntimeError(f"Unknown dry-MARTINI macro '{macro}' in {ff_path}")
                sigma_nm, epsilon_kj = macros[macro]
            else:
                sigma_nm = float(parts[3])
                epsilon_kj = float(parts[4])

            payload = {
                "sigma_nm": sigma_nm,
                "epsilon_kj_mol": epsilon_kj,
            }
            pair_params[(type_i, type_j)] = payload
            pair_params[(type_j, type_i)] = payload

    return atomtypes, pair_params


def _parse_target_spec(spec: str) -> Dict[str, Any]:
    parts = [part.strip() for part in spec.split(":")]
    if len(parts) == 1:
        label = parts[0]
        bead_type = label
        charge = _infer_type_charge(bead_type)
    elif len(parts) == 2:
        label, bead_type = parts
        charge = _infer_type_charge(bead_type)
    elif len(parts) == 3:
        label, bead_type, charge_raw = parts
        charge = float(charge_raw)
    else:
        raise ValueError(f"Invalid target spec '{spec}'. Use label[:bead_type[:charge]].")

    return {
        "label": label,
        "bead_type": bead_type,
        "charge": charge,
    }


def _default_target_atomtypes(atomtypes: List[str]) -> List[str]:
    return list(atomtypes)


def _build_targets(specs: Iterable[str], atomtypes: List[str]) -> List[Dict[str, Any]]:
    if specs:
        return [_parse_target_spec(spec) for spec in specs]
    return [
        {
            "label": atomtype,
            "bead_type": atomtype,
            "charge": _infer_type_charge(atomtype),
        }
        for atomtype in _default_target_atomtypes(atomtypes)
    ]


def _linspace(start: float, stop: float, count: int) -> List[float]:
    if count < 2:
        return [start]
    step = (stop - start) / float(count - 1)
    return [start + step * i for i in range(count)]


def _fibonacci_sphere(count: int) -> List[List[float]]:
    if count < 1:
        raise ValueError("direction count must be positive")
    if count == 1:
        return [[0.0, 0.0, 1.0]]

    directions: List[List[float]] = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for idx in range(count):
        y = 1.0 - (2.0 * idx) / float(count - 1)
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden_angle * idx
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        directions.append([x, y, z])
    return directions


def _dot3(a: Iterable[float], b: Iterable[float]) -> float:
    ax, ay, az = a
    bx, by, bz = b
    return float(ax) * float(bx) + float(ay) * float(by) + float(az) * float(bz)


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _mean_over_prefix_axes(values: Any) -> List[float]:
    shape = getattr(values, "shape", ())
    if not shape:
        return [float(values)]
    if len(shape) == 1:
        return [float(x) for x in values[:]]

    total = 1
    for dim in shape[:-1]:
        total *= int(dim)
    sums = [0.0] * int(shape[-1])
    for prefix in _iter_index_prefix(tuple(int(x) for x in shape[:-1])):
        row = values[prefix]
        for i, val in enumerate(row):
            sums[i] += float(val)
    return [val / float(total) for val in sums]


def _iter_index_prefix(shape: Tuple[int, ...]) -> Iterable[Tuple[int, ...]]:
    if not shape:
        yield ()
        return
    first, rest = shape[0], shape[1:]
    for i in range(first):
        for tail in _iter_index_prefix(rest):
            yield (i,) + tail


def _load_sidechain_orientation_library(sidechain_lib_path: Path) -> Dict[str, Dict[str, Any]]:
    if not sidechain_lib_path.exists():
        raise RuntimeError(f"Sidechain orientation library not found: {sidechain_lib_path}")

    with h5py.File(sidechain_lib_path, "r") as h5:
        restype_order = [x.decode("ascii") if isinstance(x, bytes) else str(x) for x in h5["restype_order"][:]]
        start_stop_bead = h5["rotamer_start_stop_bead"][:]
        rotamer_center_fixed = h5["rotamer_center_fixed"][:, :6]

        prob_weights = None
        if "rotamer_prob_fixed" in h5:
            prob_weights = [float(x) for x in h5["rotamer_prob_fixed"][:].reshape(-1)]
        elif "rotamer_prob" in h5:
            prob_weights = _mean_over_prefix_axes(h5["rotamer_prob"])

        residue_info: Dict[str, Dict[str, Any]] = {}
        for residue_index, residue_name in enumerate(restype_order):
            start, stop, n_bead = [int(x) for x in start_stop_bead[residue_index]]
            if stop <= start or n_bead <= 0:
                residue_info[residue_name] = {
                    "center_nm": [],
                    "vector_unit": [],
                    "weight": [],
                }
                continue

            n_row = stop - start
            if n_row % n_bead != 0:
                raise RuntimeError(
                    f"Sidechain orientation rows are not divisible by n_bead for residue {residue_name}: "
                    f"rows={n_row} n_bead={n_bead}"
                )
            n_rot = n_row // n_bead
            block = rotamer_center_fixed[start:stop].reshape(n_rot, n_bead, 6)
            centers_nm = [
                [float(val) * ANGSTROM_TO_NM for val in row]
                for row in block[:, :, :3].mean(axis=1)
            ]
            vectors = []
            for row in block[:, 0, 3:6]:
                norm = math.sqrt(sum(float(x) * float(x) for x in row))
                if norm <= 0.0:
                    raise RuntimeError(f"Invalid zero sidechain vector in {sidechain_lib_path} for residue {residue_name}")
                vectors.append([float(x) / norm for x in row])

            if prob_weights is None:
                weights = [1.0 / float(n_rot)] * n_rot
            else:
                row_weights = prob_weights[start:stop:n_bead]
                if len(row_weights) != n_rot:
                    raise RuntimeError(
                        f"Rotamer-weight shape mismatch in {sidechain_lib_path} for residue {residue_name}: "
                        f"expected {n_rot}, got {len(row_weights)}"
                    )
                total = sum(max(0.0, float(x)) for x in row_weights)
                if total <= 0.0:
                    weights = [1.0 / float(n_rot)] * n_rot
                else:
                    weights = [max(0.0, float(x)) / total for x in row_weights]

            residue_info[residue_name] = {
                "center_nm": centers_nm,
                "vector_unit": vectors,
                "weight": weights,
                "n_rotamer": n_rot,
                "n_bead": n_bead,
            }

    return residue_info


def _accumulate_on_cos_grid(
    cos_theta: float,
    value: float,
    cos_grid: List[float],
    value_sum: List[float],
    weight_sum: List[float],
) -> None:
    if len(cos_grid) == 1:
        value_sum[0] += value
        weight_sum[0] += 1.0
        return

    step = cos_grid[1] - cos_grid[0]
    coord = (cos_theta - cos_grid[0]) / step
    if coord <= 0.0:
        value_sum[0] += value
        weight_sum[0] += 1.0
        return
    if coord >= len(cos_grid) - 1:
        value_sum[-1] += value
        weight_sum[-1] += 1.0
        return

    lo = int(math.floor(coord))
    hi = lo + 1
    hi_weight = coord - float(lo)
    lo_weight = 1.0 - hi_weight
    value_sum[lo] += lo_weight * value
    weight_sum[lo] += lo_weight
    value_sum[hi] += hi_weight * value
    weight_sum[hi] += hi_weight


def _factorize_one_sided_orientation(
    sampled_energy_grid: List[List[float]],
    cos_theta_grid: List[float],
) -> Tuple[List[float], List[float], List[float], List[List[float]], float]:
    sampled = np.asarray(sampled_energy_grid, dtype=np.float64)
    if sampled.ndim != 2:
        raise RuntimeError(f"Expected a 2D sampled energy grid, got shape {sampled.shape}")

    radial = sampled.mean(axis=0)
    residual = sampled - radial[None, :]

    if np.allclose(residual, 0.0):
        angular_profile = np.zeros(sampled.shape[0], dtype=np.float64)
        angular_radial = np.zeros(sampled.shape[1], dtype=np.float64)
        fitted = radial[None, :].repeat(sampled.shape[0], axis=0)
        rms_error = 0.0
    else:
        u, s, vh = np.linalg.svd(residual, full_matrices=False)
        angular_profile = u[:, 0].copy()
        angular_radial = (s[0] * vh[0, :]).copy()

        profile_mean = float(angular_profile.mean())
        if abs(profile_mean) > 1.0e-12:
            radial = radial + profile_mean * angular_radial
            angular_profile = angular_profile - profile_mean

        if float(np.dot(angular_profile, np.asarray(cos_theta_grid, dtype=np.float64))) < 0.0:
            angular_profile *= -1.0
            angular_radial *= -1.0

        max_abs = float(np.max(np.abs(angular_profile)))
        if max_abs > 0.0:
            angular_profile /= max_abs
            angular_radial *= max_abs

        fitted = radial[None, :] + angular_profile[:, None] * angular_radial[None, :]
        rms_error = float(np.sqrt(np.mean((sampled - fitted) ** 2)))

    return (
        [float(x) for x in radial],
        [float(x) for x in angular_profile],
        [float(x) for x in angular_radial],
        [[float(x) for x in row] for row in fitted.tolist()],
        rms_error,
    )


def _task_code(residue: str, target_label: str) -> str:
    safe_label = "".join(ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in target_label)
    return f"{residue}__{safe_label}"


def _build_manifest(
    base_dir: Path,
    martinize_path: Path,
    dry_ff_path: Path,
    sidechain_lib_path: Path,
    forcefield_name: str,
    target_specs: List[str],
    r_min_nm: float,
    r_max_nm: float,
    r_count: int,
    direction_count: int,
    cos_theta_count: int,
    benchmark_script: Path,
) -> Dict[str, Any]:
    residue_map = _load_martini_forcefield(martinize_path, forcefield_name)
    atomtypes, pair_params = _parse_dry_forcefield(dry_ff_path)
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    targets = _build_targets(target_specs, atomtypes)
    directions = _fibonacci_sphere(direction_count)
    cos_theta_grid = _linspace(-1.0, 1.0, cos_theta_count)

    tasks: List[Dict[str, Any]] = []
    for residue in sorted(residue_map):
        bead_types = residue_map[residue]
        if not bead_types:
            continue
        orientation = orientation_map.get(residue)
        if not orientation or not orientation["center_nm"]:
            raise RuntimeError(f"Missing orientation geometry for residue {residue} in {sidechain_lib_path}")
        bead_charges = [_infer_type_charge(bead_type) for bead_type in bead_types]
        for target in targets:
            missing = [bead for bead in bead_types if (bead, target["bead_type"]) not in pair_params]
            if missing:
                raise RuntimeError(
                    f"Missing dry-MARTINI nonbond parameter for residue {residue} "
                    f"target {target['bead_type']} via beads {sorted(set(missing))}"
                )
            tasks.append(
                {
                    "task_id": len(tasks),
                    "code": _task_code(residue, target["label"]),
                    "residue": residue,
                    "sidechain_bead_types": bead_types,
                    "sidechain_bead_charges": bead_charges,
                    "sidechain_effective_center_nm": orientation["center_nm"],
                    "sidechain_effective_vector_unit": orientation["vector_unit"],
                    "sidechain_rotamer_weight": orientation["weight"],
                    "target": dict(target),
                }
            )

    return {
        "schema": SCHEMA_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "workflow_dir": str(_workflow_dir()),
        "project_root": str(_project_root()),
        "martinize_path": str(martinize_path),
        "sidechain_source_path": str(martinize_path),
        "dry_forcefield_path": str(dry_ff_path),
        "sidechain_library_path": str(sidechain_lib_path),
        "forcefield_name": forcefield_name,
        "benchmark_script": str(benchmark_script),
        "grid": {
            "r_min_nm": r_min_nm,
            "r_max_nm": r_max_nm,
            "r_count": r_count,
            "direction_count": direction_count,
            "direction_vectors_unit": directions,
            "cos_theta_count": cos_theta_count,
            "cos_theta_grid": cos_theta_grid,
            "cb_anchor_nm": [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG],
            "cb_vector_unit": list(CANONICAL_CB_VECTOR_UNIT),
            "sample_count_per_task": r_count * direction_count,
        },
        "assumptions": {
            "effective_model": "sum_beadwise_rotamer_weighted_oriented_cb_table",
            "geometry_note": (
                "Training samples target dry-MARTINI particles in the original Upside residue-local CB frame. "
                "Each residue uses the sidechain library's rotamer-resolved effective centers, projected into a "
                "single dry-MARTINI table over radial distance and cos(theta) with respect to the canonical CA->CB vector."
            ),
            "target_selection": {
                "mode": "default_all_dry_martini_atomtypes" if not target_specs else "explicit_user_targets",
                "default_target_count": len(targets),
            },
            "units": {
                "distance": "nm",
                "energy": "kJ/mol",
                "charge": "e",
            },
            "coulomb_constant_kj_mol_nm_e2": COULOMB_K_DRY_KJ_NM,
            "runtime_unit_policy": {
                "native_units_only": True,
                "note": (
                    "Training outputs remain in native dry-MARTINI units. Simulation-side "
                    "conversion must be supplied separately as runtime parameters."
                ),
            },
            "orientation_contract": {
                "anchor": "CB",
                "vector": "canonical_CA_to_CB_unit_vector",
                "n1_source": "backbone_defined_CA_to_CB_vector",
                "n12_source": "sampled_CB_to_target_direction",
                "angular_coordinate": "cos_theta",
                "note": (
                    "The one-sided orientation term uses the same backbone-defined CB vector in training "
                    "and runtime. The sampled dry-MARTINI target position contributes only the CB-to-target "
                    "direction used for n12."
                ),
                "sidechain_library_units": "angstrom",
                "training_table_units": "nm / kJ/mol / e",
            },
        },
        "targets": targets,
        "tasks": tasks,
    }


def _pair_energy(r_nm: float, sigma_nm: float, epsilon_kj_mol: float, qi: float, qj: float) -> Tuple[float, float, float]:
    sr = sigma_nm / r_nm
    sr2 = sr * sr
    sr6 = sr2 * sr2 * sr2
    sr12 = sr6 * sr6
    lj = 4.0 * epsilon_kj_mol * (sr12 - sr6)
    coulomb = COULOMB_K_DRY_KJ_NM * qi * qj / r_nm if qi and qj else 0.0
    return lj + coulomb, lj, coulomb


def _run_task(base_dir: Path, manifest: Dict[str, Any], task: Dict[str, Any]) -> Path:
    dry_ff_path = Path(manifest["dry_forcefield_path"]).expanduser().resolve()
    _atomtypes, pair_params = _parse_dry_forcefield(dry_ff_path)
    grid = manifest["grid"]
    r_values = _linspace(float(grid["r_min_nm"]), float(grid["r_max_nm"]), int(grid["r_count"]))
    direction_vectors = [list(vec) for vec in grid["direction_vectors_unit"]]
    cos_theta_grid = [float(x) for x in grid["cos_theta_grid"]]
    cb_anchor_nm = [float(x) for x in grid["cb_anchor_nm"]]
    cb_vector_unit = [float(x) for x in grid["cb_vector_unit"]]
    target = task["target"]
    target_type = str(target["bead_type"])
    target_charge = float(target["charge"])
    rotamer_centers_nm = [list(row) for row in task["sidechain_effective_center_nm"]]
    rotamer_weights = [float(x) for x in task["sidechain_rotamer_weight"]]
    rotamer_vectors = [list(row) for row in task["sidechain_effective_vector_unit"]]
    n_rotamer = len(rotamer_centers_nm)

    angular_energy = [[0.0 for _ in r_values] for _ in cos_theta_grid]
    angular_lj = [[0.0 for _ in r_values] for _ in cos_theta_grid]
    angular_coulomb = [[0.0 for _ in r_values] for _ in cos_theta_grid]
    rotamer_angular_energy = [
        [[0.0 for _ in r_values] for _ in cos_theta_grid]
        for _ in range(n_rotamer)
    ]
    position_samples_nm: List[List[float]] = []
    sample_cos_theta: List[float] = []
    sample_energy: List[float] = []
    sample_lj: List[float] = []
    sample_coulomb: List[float] = []
    component_details: List[Dict[str, Any]] = []

    for residue_bead, bead_charge in zip(task["sidechain_bead_types"], task["sidechain_bead_charges"]):
        params = pair_params[(residue_bead, target_type)]
        component_details.append(
            {
                "bead_type": residue_bead,
                "bead_charge": bead_charge,
                "sigma_nm": params["sigma_nm"],
                "epsilon_kj_mol": params["epsilon_kj_mol"],
            }
        )

    for ir, r_nm in enumerate(r_values):
        energy_sum = [0.0 for _ in cos_theta_grid]
        energy_weight = [0.0 for _ in cos_theta_grid]
        lj_sum = [0.0 for _ in cos_theta_grid]
        lj_weight = [0.0 for _ in cos_theta_grid]
        coul_sum = [0.0 for _ in cos_theta_grid]
        coul_weight = [0.0 for _ in cos_theta_grid]
        rotamer_energy_sum = [[0.0 for _ in cos_theta_grid] for _ in range(n_rotamer)]
        rotamer_energy_weight = [[0.0 for _ in cos_theta_grid] for _ in range(n_rotamer)]

        for direction in direction_vectors:
            target_position_nm = [
                cb_anchor_nm[0] + r_nm * float(direction[0]),
                cb_anchor_nm[1] + r_nm * float(direction[1]),
                cb_anchor_nm[2] + r_nm * float(direction[2]),
            ]
            # n1 is the residue-local backbone-defined CB vector. The sampled target
            # direction contributes only n12 for the one-sided angular coordinate.
            cos_theta = _clamp(-_dot3(direction, cb_vector_unit), -1.0, 1.0)
            total_energy = 0.0
            total_lj = 0.0
            total_coulomb = 0.0

            for irot, (center_nm, rot_weight) in enumerate(zip(rotamer_centers_nm, rotamer_weights)):
                dx = target_position_nm[0] - float(center_nm[0])
                dy = target_position_nm[1] - float(center_nm[1])
                dz = target_position_nm[2] - float(center_nm[2])
                dist_nm = max(1.0e-6, math.sqrt(dx * dx + dy * dy + dz * dz))

                rot_energy = 0.0
                rot_lj = 0.0
                rot_coulomb = 0.0
                for residue_bead, bead_charge in zip(task["sidechain_bead_types"], task["sidechain_bead_charges"]):
                    params = pair_params[(residue_bead, target_type)]
                    pair_total, pair_lj, pair_coul = _pair_energy(
                        r_nm=dist_nm,
                        sigma_nm=float(params["sigma_nm"]),
                        epsilon_kj_mol=float(params["epsilon_kj_mol"]),
                        qi=float(bead_charge),
                        qj=target_charge,
                    )
                    rot_energy += pair_total
                    rot_lj += pair_lj
                    rot_coulomb += pair_coul

                _accumulate_on_cos_grid(
                    cos_theta,
                    rot_energy,
                    cos_theta_grid,
                    rotamer_energy_sum[irot],
                    rotamer_energy_weight[irot],
                )
                total_energy += rot_weight * rot_energy
                total_lj += rot_weight * rot_lj
                total_coulomb += rot_weight * rot_coulomb

            position_samples_nm.append([r_nm * float(direction[0]), r_nm * float(direction[1]), r_nm * float(direction[2])])
            sample_cos_theta.append(cos_theta)
            sample_energy.append(total_energy)
            sample_lj.append(total_lj)
            sample_coulomb.append(total_coulomb)
            _accumulate_on_cos_grid(cos_theta, total_energy, cos_theta_grid, energy_sum, energy_weight)
            _accumulate_on_cos_grid(cos_theta, total_lj, cos_theta_grid, lj_sum, lj_weight)
            _accumulate_on_cos_grid(cos_theta, total_coulomb, cos_theta_grid, coul_sum, coul_weight)

        for ia in range(len(cos_theta_grid)):
            if energy_weight[ia] <= 0.0 or lj_weight[ia] <= 0.0 or coul_weight[ia] <= 0.0:
                raise RuntimeError(
                    f"Cos(theta) bin received no samples for residue {task['residue']} target {target_type} "
                    f"at r={r_nm:.4f} nm; increase direction-count or reduce cos-theta-count"
                )
            angular_energy[ia][ir] = energy_sum[ia] / energy_weight[ia]
            angular_lj[ia][ir] = lj_sum[ia] / lj_weight[ia]
            angular_coulomb[ia][ir] = coul_sum[ia] / coul_weight[ia]
            for irot in range(n_rotamer):
                if rotamer_energy_weight[irot][ia] <= 0.0:
                    raise RuntimeError(
                        f"Cos(theta) bin received no rotamer samples for residue {task['residue']} target {target_type} "
                        f"rotamer={irot} at r={r_nm:.4f} nm"
                    )
                rotamer_angular_energy[irot][ia][ir] = (
                    rotamer_energy_sum[irot][ia] / rotamer_energy_weight[irot][ia]
                )

    radial_energy, angular_profile, angular_radial_energy, fitted_energy, factorization_rms_error = (
        _factorize_one_sided_orientation(angular_energy, cos_theta_grid)
    )
    rotamer_radial_energy = []
    rotamer_angular_profile = []
    rotamer_angular_radial_energy = []
    rotamer_fitted_energy = []
    rotamer_factorization_rms_error = []
    for rot_grid in rotamer_angular_energy:
        rot_radial, rot_profile, rot_angular, rot_fitted, rot_rms = _factorize_one_sided_orientation(
            rot_grid, cos_theta_grid
        )
        rotamer_radial_energy.append(rot_radial)
        rotamer_angular_profile.append(rot_profile)
        rotamer_angular_radial_energy.append(rot_angular)
        rotamer_fitted_energy.append(rot_fitted)
        rotamer_factorization_rms_error.append(rot_rms)

    payload = {
        "schema": SCHEMA_TASK,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "task_id": int(task["task_id"]),
        "code": task["code"],
        "residue": task["residue"],
        "target": target,
        "sidechain_bead_types": list(task["sidechain_bead_types"]),
        "sidechain_bead_charges": list(task["sidechain_bead_charges"]),
        "sidechain_effective_center_nm": rotamer_centers_nm,
        "sidechain_effective_vector_unit": rotamer_vectors,
        "sidechain_rotamer_weight": rotamer_weights,
        "cb_anchor_nm": cb_anchor_nm,
        "cb_vector_unit": cb_vector_unit,
        "component_pairs": component_details,
        "grid_nm": r_values,
        "cos_theta_grid": cos_theta_grid,
        "direction_vectors_unit": direction_vectors,
        "position_samples_nm": position_samples_nm,
        "sample_cos_theta": sample_cos_theta,
        "sampled_grid_energy_kj_mol": angular_energy,
        "sampled_grid_lj_kj_mol": angular_lj,
        "sampled_grid_coulomb_kj_mol": angular_coulomb,
        "energy_kj_mol": fitted_energy,
        "radial_energy_kj_mol": radial_energy,
        "angular_energy_kj_mol": angular_radial_energy,
        "angular_profile": angular_profile,
        "fit_rms_error_kj_mol": factorization_rms_error,
        "rotamer_count": n_rotamer,
        "rotamer_probability_fixed": rotamer_weights,
        "rotamer_sampled_grid_energy_kj_mol": rotamer_angular_energy,
        "rotamer_energy_kj_mol": rotamer_fitted_energy,
        "rotamer_radial_energy_kj_mol": rotamer_radial_energy,
        "rotamer_angular_energy_kj_mol": rotamer_angular_radial_energy,
        "rotamer_angular_profile": rotamer_angular_profile,
        "rotamer_fit_rms_error_kj_mol": rotamer_factorization_rms_error,
        "sample_energy_kj_mol": sample_energy,
        "sample_lj_kj_mol": sample_lj,
        "sample_coulomb_kj_mol": sample_coulomb,
        "residue_effective_charge": float(sum(task["sidechain_bead_charges"])),
        "assumption": manifest["assumptions"]["effective_model"],
    }
    out_path = _task_result_dir(base_dir) / f"{task['code']}.json"
    _write_json(out_path, payload)
    return out_path


def _load_manifest_for_base(base_dir: Path) -> Dict[str, Any]:
    manifest_path = _manifest_path(base_dir)
    if not manifest_path.exists():
        raise RuntimeError(f"Training manifest not found: {manifest_path}")
    return _load_json(manifest_path)


def _find_task(manifest: Dict[str, Any], task_id: int) -> Dict[str, Any]:
    for task in manifest["tasks"]:
        if int(task["task_id"]) == int(task_id):
            return task
    raise RuntimeError(f"Task id {task_id} not found in manifest")


def assemble_results(base_dir: Path) -> Path:
    manifest = _load_manifest_for_base(base_dir)
    assembled_entries: List[Dict[str, Any]] = []
    by_residue: Dict[str, Dict[str, Any]] = {}

    for task in manifest["tasks"]:
        result_path = _task_result_dir(base_dir) / f"{task['code']}.json"
        if not result_path.exists():
            raise RuntimeError(f"Missing task result: {result_path}")
        result = _load_json(result_path)
        assembled_entries.append(result)
        by_residue.setdefault(result["residue"], {})[result["target"]["label"]] = {
            "target": result["target"],
            "component_pairs": result["component_pairs"],
            "grid_nm": result["grid_nm"],
            "cos_theta_grid": result["cos_theta_grid"],
            "direction_vectors_unit": result["direction_vectors_unit"],
            "position_samples_nm": result["position_samples_nm"],
            "sample_cos_theta": result["sample_cos_theta"],
            "energy_kj_mol": result["energy_kj_mol"],
            "sampled_grid_energy_kj_mol": result["sampled_grid_energy_kj_mol"],
            "sampled_grid_lj_kj_mol": result["sampled_grid_lj_kj_mol"],
            "sampled_grid_coulomb_kj_mol": result["sampled_grid_coulomb_kj_mol"],
            "radial_energy_kj_mol": result["radial_energy_kj_mol"],
            "angular_energy_kj_mol": result["angular_energy_kj_mol"],
            "angular_profile": result["angular_profile"],
            "fit_rms_error_kj_mol": result["fit_rms_error_kj_mol"],
            "rotamer_count": result["rotamer_count"],
            "rotamer_probability_fixed": result["rotamer_probability_fixed"],
            "rotamer_sampled_grid_energy_kj_mol": result["rotamer_sampled_grid_energy_kj_mol"],
            "rotamer_energy_kj_mol": result["rotamer_energy_kj_mol"],
            "rotamer_radial_energy_kj_mol": result["rotamer_radial_energy_kj_mol"],
            "rotamer_angular_energy_kj_mol": result["rotamer_angular_energy_kj_mol"],
            "rotamer_angular_profile": result["rotamer_angular_profile"],
            "rotamer_fit_rms_error_kj_mol": result["rotamer_fit_rms_error_kj_mol"],
            "sample_energy_kj_mol": result["sample_energy_kj_mol"],
            "sample_lj_kj_mol": result["sample_lj_kj_mol"],
            "sample_coulomb_kj_mol": result["sample_coulomb_kj_mol"],
            "sidechain_effective_center_nm": result["sidechain_effective_center_nm"],
            "sidechain_effective_vector_unit": result["sidechain_effective_vector_unit"],
            "sidechain_rotamer_weight": result["sidechain_rotamer_weight"],
            "cb_anchor_nm": result["cb_anchor_nm"],
            "cb_vector_unit": result["cb_vector_unit"],
            "residue_effective_charge": result["residue_effective_charge"],
        }

    payload = {
        "schema": SCHEMA_TABLE,
        "created_at_utc": _now_utc(),
        "manifest_path": str(_manifest_path(base_dir)),
        "assumptions": manifest["assumptions"],
        "grid": manifest["grid"],
        "forcefield_name": manifest["forcefield_name"],
        "entries": assembled_entries,
        "tables_by_residue": by_residue,
    }
    out_path = _assembled_dir(base_dir) / "sc_table.json"
    _write_json(out_path, payload)

    summary = {
        "schema": f"{SCHEMA_TABLE}_summary",
        "created_at_utc": _now_utc(),
        "manifest_path": str(_manifest_path(base_dir)),
        "n_tasks": len(assembled_entries),
        "n_residues": len(by_residue),
        "n_targets": len(manifest["targets"]),
        "residues": sorted(by_residue),
        "targets": [target["label"] for target in manifest["targets"]],
    }
    _write_json(_assembled_dir(base_dir) / "sc_table_summary.json", summary)
    return out_path


def _phase_env(phase: str, key: str) -> str:
    return os.environ.get(f"SC_TRAIN_{phase}_{key}", os.environ.get(f"SC_TRAIN_SBATCH_{key}", "")).strip()


def _optional_sbatch_directives(phase: str) -> List[str]:
    directives: List[str] = []
    for key, flag in (
        ("PARTITION", "partition"),
        ("ACCOUNT", "account"),
        ("QOS", "qos"),
        ("CONSTRAINT", "constraint"),
        ("MEM", "mem"),
    ):
        value = _phase_env(phase, key)
        if value:
            directives.append(f"#SBATCH --{flag}={value}")
    return directives


def _run_sbatch(cmd: List[str]) -> str:
    proc = sp.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or f"sbatch failed: {' '.join(cmd)}")
    return proc.stdout.strip().split(";", 1)[0]


def _array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _slurm_dir(base_dir) / "train-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=sc-train-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('SC_TRAIN_TRAIN_WALLTIME', '04:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(os.environ.get('SC_TRAIN_CPUS_PER_TASK', '1'))}",
        f"#SBATCH --array=0-{n_tasks - 1}",
    ]
    lines.extend(_optional_sbatch_directives("TRAIN"))
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"SCRIPT_DIR={shlex.quote(str(_workflow_dir()))}",
            f"ROUND_MANIFEST={shlex.quote(str(round_manifest_path))}",
            'PYTHON_BIN="${SC_TRAIN_PYTHON:-python3}"',
            "",
            "if [ -f /etc/profile.d/modules.sh ]; then",
            "  source /etc/profile.d/modules.sh",
            "fi",
            "",
            "if command -v module >/dev/null 2>&1; then",
            "  module load python/3.11.9 || true",
            "  module load cmake || true",
            "  module load openmpi || true",
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -n "${SC_TRAIN_PROJECT_ROOT:-}" ] && [ -f "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate" ]; then',
            '  source "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate"',
            'fi',
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-array-task --round-manifest "$ROUND_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID"',
        ]
    )
    return "\n".join(lines) + "\n"


def _collector_script_content(round_manifest_path: Path, base_dir: Path) -> str:
    output_path = _slurm_dir(base_dir) / "collect-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=sc-collect-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('SC_TRAIN_COLLECT_WALLTIME', '01:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
    ]
    lines.extend(_optional_sbatch_directives("COLLECT"))
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"SCRIPT_DIR={shlex.quote(str(_workflow_dir()))}",
            f"ROUND_MANIFEST={shlex.quote(str(round_manifest_path))}",
            'PYTHON_BIN="${SC_TRAIN_PYTHON:-python3}"',
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -n "${SC_TRAIN_PROJECT_ROOT:-}" ] && [ -f "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate" ]; then',
            '  source "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate"',
            'fi',
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" assemble-results --round-manifest "$ROUND_MANIFEST"',
        ]
    )
    return "\n".join(lines) + "\n"


def _benchmark_script_content(base_dir: Path) -> str:
    output_path = _slurm_dir(base_dir) / "benchmark-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=sc-bench-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('SC_TRAIN_BENCH_WALLTIME', '24:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
    ]
    lines.extend(_optional_sbatch_directives("BENCH"))
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"SCRIPT_DIR={shlex.quote(str(_workflow_dir()))}",
            f"BASE_DIR={shlex.quote(str(base_dir))}",
            'PYTHON_BIN="${SC_TRAIN_PYTHON:-python3}"',
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -n "${SC_TRAIN_PROJECT_ROOT:-}" ] && [ -f "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate" ]; then',
            '  source "${SC_TRAIN_PROJECT_ROOT}/.venv/bin/activate"',
            'fi',
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-benchmark --base-dir "$BASE_DIR" --execute',
        ]
    )
    return "\n".join(lines) + "\n"


def submit_slurm(base_dir: Path, run_benchmark: bool, no_submit: bool) -> int:
    manifest = _load_manifest_for_base(base_dir)
    round_manifest = {
        "schema": SCHEMA_SLURM_ROUND,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "n_tasks": len(manifest["tasks"]),
        "run_benchmark": bool(run_benchmark),
    }
    round_manifest_path = _slurm_dir(base_dir) / "round_manifest.json"
    _write_json(round_manifest_path, round_manifest)

    array_script = _slurm_dir(base_dir) / "train_array.sbatch"
    collect_script = _slurm_dir(base_dir) / "collect_results.sbatch"
    benchmark_script = _slurm_dir(base_dir) / "run_benchmark.sbatch"
    _write_text(array_script, _array_script_content(round_manifest_path, base_dir, len(manifest["tasks"])), executable=True)
    _write_text(collect_script, _collector_script_content(round_manifest_path, base_dir), executable=True)
    if run_benchmark:
        _write_text(benchmark_script, _benchmark_script_content(base_dir), executable=True)

    print("Round manifest:", round_manifest_path)
    print("Array script:", array_script)
    print("Collector script:", collect_script)
    if run_benchmark:
        print("Benchmark script:", benchmark_script)

    if no_submit:
        print("Staged Slurm scripts only; no jobs submitted.")
        return 0

    if not shutil_which("sbatch"):
        raise RuntimeError("sbatch not available in PATH")

    train_job = _run_sbatch(["sbatch", "--parsable", str(array_script)])
    collect_job = _run_sbatch(["sbatch", "--parsable", f"--dependency=afterok:{train_job}", str(collect_script)])
    record: Dict[str, Any] = {
        "created_at_utc": _now_utc(),
        "train_job_id": train_job,
        "collect_job_id": collect_job,
    }
    print("Submitted training array job:", train_job)
    print("Submitted collector job:", collect_job)
    if run_benchmark:
        bench_job = _run_sbatch(["sbatch", "--parsable", f"--dependency=afterok:{collect_job}", str(benchmark_script)])
        record["benchmark_job_id"] = bench_job
        print("Submitted benchmark job:", bench_job)

    _write_json(_slurm_dir(base_dir) / "round_submission.json", record)
    return 0


def shutil_which(binary: str) -> str | None:
    for path in os.environ.get("PATH", "").split(os.pathsep):
        if not path:
            continue
        candidate = Path(path) / binary
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def run_benchmark(base_dir: Path, execute: bool) -> int:
    manifest = _load_manifest_for_base(base_dir)
    assembled_path = _assembled_dir(base_dir) / "sc_table.json"
    if not assembled_path.exists():
        raise RuntimeError(f"Assembled SC table not found: {assembled_path}")

    benchmark_script = Path(manifest["benchmark_script"]).expanduser().resolve()
    if not benchmark_script.exists():
        raise RuntimeError(f"Benchmark script not found: {benchmark_script}")

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = _benchmark_dir(base_dir) / stamp
    run_dir.mkdir(parents=True, exist_ok=True)
    table_copy = run_dir / "sc_table.json"
    table_copy.write_text(assembled_path.read_text(encoding="utf-8"), encoding="utf-8")
    env_file = run_dir / "benchmark_env.sh"
    env_file.write_text(
        "\n".join(
            [
                f"export SC_TRAINING_BASE_DIR={shlex.quote(str(base_dir))}",
                f"export SC_TRAINING_TABLE_PATH={shlex.quote(str(table_copy))}",
                f"export SC_TRAINING_WORKFLOW_DIR={shlex.quote(str(_workflow_dir()))}",
                f"export SC_TRAINING_NOTE={shlex.quote('First-pass table staged; runtime integration still required.')}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    command = (
        f"source {shlex.quote(str(env_file))} && "
        f"cd {shlex.quote(str(benchmark_script.parent))} && "
        f"bash {shlex.quote(str(benchmark_script.name))}"
    )
    command_path = run_dir / "benchmark_command.sh"
    _write_text(command_path, "#!/bin/bash\nset -euo pipefail\n" + command + "\n", executable=True)

    print("Benchmark staging directory:", run_dir)
    print("Benchmark command:", command_path)

    if not execute:
        return 0

    log_path = run_dir / "benchmark.log"
    proc = sp.run(
        ["bash", str(command_path)],
        cwd=str(benchmark_script.parent),
        capture_output=True,
        text=True,
        check=False,
    )
    log_path.write_text(
        proc.stdout + ("\n" if proc.stdout and not proc.stdout.endswith("\n") else "") + proc.stderr,
        encoding="utf-8",
    )
    print("Benchmark log:", log_path)
    if proc.returncode != 0:
        raise RuntimeError(f"Benchmark script failed with exit code {proc.returncode}")
    return 0


def cmd_init_run(args: argparse.Namespace) -> int:
    base_dir = _resolve_base_dir(args.base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)
    manifest = _build_manifest(
        base_dir=base_dir,
        martinize_path=Path(args.martinize).expanduser().resolve(),
        dry_ff_path=Path(args.dry_forcefield).expanduser().resolve(),
        sidechain_lib_path=Path(args.sidechain_library).expanduser().resolve(),
        forcefield_name=args.forcefield,
        target_specs=list(args.target or []),
        r_min_nm=float(args.r_min_nm),
        r_max_nm=float(args.r_max_nm),
        r_count=int(args.r_count),
        direction_count=int(args.direction_count),
        cos_theta_count=int(args.cos_theta_count),
        benchmark_script=Path(args.benchmark_script).expanduser().resolve(),
    )
    _write_json(_manifest_path(base_dir), manifest)
    _record_current_run_dir(base_dir)
    print("Initialized SC-training run:", base_dir)
    print("Manifest:", _manifest_path(base_dir))
    print("Tasks:", len(manifest["tasks"]))
    return 0


def cmd_run_local(args: argparse.Namespace) -> int:
    base_dir = _resolve_base_dir(args.base_dir)
    manifest = _load_manifest_for_base(base_dir)
    _record_current_run_dir(base_dir)
    for task in manifest["tasks"]:
        result_path = _task_result_dir(base_dir) / f"{task['code']}.json"
        if result_path.exists() and not args.force:
            continue
        _run_task(base_dir, manifest, task)
    assembled = assemble_results(base_dir)
    print("Assembled table:", assembled)
    return 0


def cmd_run_array_task(args: argparse.Namespace) -> int:
    round_manifest = _load_json(Path(args.round_manifest).expanduser().resolve())
    base_dir = Path(round_manifest["base_dir"]).expanduser().resolve()
    manifest = _load_manifest_for_base(base_dir)
    task = _find_task(manifest, int(args.task_id))
    result_path = _run_task(base_dir, manifest, task)
    print("Completed task:", result_path)
    return 0


def cmd_assemble(args: argparse.Namespace) -> int:
    if args.round_manifest:
        round_manifest = _load_json(Path(args.round_manifest).expanduser().resolve())
        base_dir = Path(round_manifest["base_dir"]).expanduser().resolve()
    else:
        base_dir = _resolve_base_dir(args.base_dir)
    assembled = assemble_results(base_dir)
    print("Assembled table:", assembled)
    return 0


def cmd_submit_slurm(args: argparse.Namespace) -> int:
    base_dir = _resolve_base_dir(args.base_dir)
    _record_current_run_dir(base_dir)
    return submit_slurm(base_dir=base_dir, run_benchmark=bool(args.run_benchmark), no_submit=bool(args.no_submit))


def cmd_run_benchmark(args: argparse.Namespace) -> int:
    base_dir = _resolve_base_dir(args.base_dir)
    return run_benchmark(base_dir=base_dir, execute=bool(args.execute))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="SC training workflow")
    sub = parser.add_subparsers(dest="cmd", required=True)

    default_martinize = _project_root() / "example" / "16.MARTINI" / "martinize.py"
    default_dry_ff = _project_root() / "example" / "16.MARTINI" / "ff_dry" / "dry_martini_v2.1.itp"
    default_sidechain_library = _default_sidechain_library()
    default_benchmark = _project_root() / "example" / "16.MARTINI" / "run_sim_1rkl.sh"

    p_init = sub.add_parser("init-run", help="Create a training manifest in a run directory")
    p_init.add_argument("--base-dir", help="Run directory")
    p_init.add_argument(
        "--martinize",
        default=str(default_martinize),
        help="Path to martinize.py or a sidechain-definition JSON",
    )
    p_init.add_argument("--dry-forcefield", default=str(default_dry_ff), help="Path to dry MARTINI .itp")
    p_init.add_argument(
        "--sidechain-library",
        default=str(default_sidechain_library),
        help="Path to the original Upside sidechain.h5 library used for residue-local CB-frame geometry",
    )
    p_init.add_argument("--forcefield", default="martini22", help="Forcefield class name inside martinize.py")
    p_init.add_argument("--target", action="append", help="Target spec label[:bead_type[:charge]]; repeatable")
    p_init.add_argument("--r-min-nm", type=float, default=0.25, help="Minimum training radius in nm")
    p_init.add_argument("--r-max-nm", type=float, default=1.20, help="Maximum training radius in nm")
    p_init.add_argument("--r-count", type=int, default=96, help="Number of radial samples")
    p_init.add_argument(
        "--direction-count",
        type=int,
        default=24,
        help="Number of spherical directions sampled around the sidechain at each radius",
    )
    p_init.add_argument(
        "--cos-theta-count",
        type=int,
        default=13,
        help="Number of uniform cos(theta) bins in the assembled orientation-aware table",
    )
    p_init.add_argument("--benchmark-script", default=str(default_benchmark), help="Canonical benchmark workflow script")
    p_init.set_defaults(func=cmd_init_run)

    p_local = sub.add_parser("run-local", help="Run all tasks locally and assemble results")
    p_local.add_argument("--base-dir", help="Run directory")
    p_local.add_argument("--force", action="store_true", help="Recompute task outputs even if they already exist")
    p_local.set_defaults(func=cmd_run_local)

    p_array = sub.add_parser("run-array-task", help="Run a single staged Slurm array task")
    p_array.add_argument("--round-manifest", required=True, help="Generated Slurm round manifest JSON")
    p_array.add_argument("--task-id", required=True, type=int, help="Task index")
    p_array.set_defaults(func=cmd_run_array_task)

    p_assemble = sub.add_parser("assemble-results", help="Assemble completed task outputs")
    p_assemble.add_argument("--base-dir", help="Run directory")
    p_assemble.add_argument("--round-manifest", help="Optional generated Slurm round manifest JSON")
    p_assemble.set_defaults(func=cmd_assemble)

    p_submit = sub.add_parser("submit-slurm", help="Stage and optionally submit Slurm jobs")
    p_submit.add_argument("--base-dir", help="Run directory")
    p_submit.add_argument("--run-benchmark", action="store_true", help="Submit a dependent benchmark job after assembly")
    p_submit.add_argument("--no-submit", action="store_true", help="Generate Slurm scripts without calling sbatch")
    p_submit.set_defaults(func=cmd_submit_slurm)

    p_bench = sub.add_parser("run-benchmark", help="Stage or execute the canonical benchmark workflow")
    p_bench.add_argument("--base-dir", help="Run directory")
    p_bench.add_argument("--execute", action="store_true", help="Execute the benchmark instead of only staging the command")
    p_bench.set_defaults(func=cmd_run_benchmark)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
