#!/usr/bin/env python3
"""Build dry-MARTINI spline tables from ITP-derived parameters."""

from __future__ import annotations

import math
import os
import shutil
import importlib.util
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Set, Tuple

import h5py
import numpy as np

from martini_cg_lipid_params import derive_dopc_cg_params
from martini_itp_reader import (
    infer_charge_from_atomtype,
    load_martini_forcefield,
    parse_dry_forcefield,
    parse_itp_atomtype_masses,
)

COULOMB_K_DRY_KJ_NM = 138.935458 / 15.0
ENERGY_CONVERSION_KJ_PER_EUP = 2.914952774272
LENGTH_CONVERSION_A_PER_NM = 10.0
ANGSTROM_TO_NM = 0.1
DEFAULT_PRODUCTION_TEMP_UPSIDE = 0.8647
DEFAULT_PRODUCTION_KBT_KJ_MOL = DEFAULT_PRODUCTION_TEMP_UPSIDE * ENERGY_CONVERSION_KJ_PER_EUP
DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE = float(
    os.environ.get("UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE", "25.0")
)
DEFAULT_PMF_AVERAGE_TEMP_UPSIDE = DEFAULT_PRODUCTION_TEMP_UPSIDE
PARTICLES_GRID_N = 1000
PARTICLES_R_MIN_A = 0.0
PARTICLES_R_MAX_A = 12.0
DRY_MARTINI_NONBONDED_CUTOFF_NM = PARTICLES_R_MAX_A * ANGSTROM_TO_NM
NUMERICAL_DISTANCE_GUARD_NM = 1.0e-6

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)

CANONICAL_CB_POSITION_ANG = (0.0, 0.94375626, 1.2068012)
_cb_norm = math.sqrt(sum(x * x for x in CANONICAL_CB_POSITION_ANG))
CANONICAL_CB_VECTOR_UNIT = tuple(x / _cb_norm for x in CANONICAL_CB_POSITION_ANG)

SCHEMA_PARTICLES = "martini_particles_combined"
SCHEMA_SC = "martini_sc_combined"

# Lazy-loaded from ITP at first CG lipid table build.
_CURRENT_CG_BONDS: list | None = None
_CURRENT_CG_ANGLES: list | None = None


def _positive_int_env(name: str, default: int) -> int:
    text = os.environ.get(name, "").strip()
    if not text:
        return int(default)
    try:
        value = int(text)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {text!r}") from exc
    return max(1, value)


def _float_env(name: str, default: float) -> float:
    text = os.environ.get(name, "").strip()
    if not text:
        return float(default)
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f"{name} must be a float, got {text!r}") from exc


def _table_worker_count(task_count: int = 0) -> int:
    explicit = os.environ.get("UPSIDE_MARTINI_TABLE_WORKERS", "").strip()
    if explicit:
        return min(_positive_int_env("UPSIDE_MARTINI_TABLE_WORKERS", 1), max(1, task_count or 1))

    for name in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        text = os.environ.get(name, "").strip()
        if text:
            return min(_positive_int_env(name, 1), max(1, task_count or 1))

    cpu_count = os.cpu_count() or 1
    default = max(1, min(cpu_count - 1 if cpu_count > 1 else 1, 8))
    return min(default, max(1, task_count or 1))


def _parallel_map_ordered(label: str, func: Callable[[Any], Any], tasks: list) -> list:
    workers = _table_worker_count(len(tasks))
    if workers <= 1 or len(tasks) <= 1:
        return [func(task) for task in tasks]

    context_name = os.environ.get("UPSIDE_MARTINI_MP_CONTEXT", "fork").strip() or "fork"
    try:
        context = mp.get_context(context_name)
    except ValueError:
        context = mp.get_context()
    print(f"  Parallel {label}: {len(tasks)} tasks on {workers} worker(s)")
    chunksize = max(1, min(128, (len(tasks) + workers * 8 - 1) // (workers * 8)))
    try:
        with context.Pool(processes=workers) as pool:
            return list(pool.imap(func, tasks, chunksize=chunksize))
    except (OSError, PermissionError) as exc:
        print(f"  Parallel {label}: process workers unavailable ({exc}); using threads")
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(func, tasks))


def _bead_frame_angles(count: int) -> np.ndarray:
    count = max(1, int(count))
    return np.linspace(0.0, 2.0 * np.pi, count, endpoint=False, dtype=np.float64)


def _boltzmann_free_energy_kj_mol(
    values_kj_mol: list[float] | np.ndarray,
    temperature: float,
    weights: list[float] | np.ndarray | None = None,
) -> float:
    values = np.asarray(values_kj_mol, dtype=np.float64)
    if values.size == 0:
        return 0.0
    if weights is None:
        weights_arr = np.ones(values.shape, dtype=np.float64)
    else:
        weights_arr = np.asarray(weights, dtype=np.float64)
    positive = weights_arr > 0.0
    if not np.any(positive):
        return float(np.mean(values))
    values = values[positive]
    weights_arr = weights_arr[positive]
    if temperature <= 0.0:
        return float(np.average(values, weights=weights_arr))

    kbt = float(temperature) * ENERGY_CONVERSION_KJ_PER_EUP
    emin = float(np.min(values))
    z = float(np.sum(weights_arr * np.exp(-(values - emin) / kbt)) / np.sum(weights_arr))
    return float(emin - kbt * math.log(max(z, 1e-300)))


def _bead_frame_count(kind: str, default: int = 1) -> int:
    default = int(default)
    specific = os.environ.get(f"UPSIDE_MARTINI_{kind.upper()}_BEAD_FRAME_COUNT", "").strip()
    if specific:
        return _positive_int_env(f"UPSIDE_MARTINI_{kind.upper()}_BEAD_FRAME_COUNT", default)
    shared = os.environ.get("UPSIDE_MARTINI_BEAD_FRAME_COUNT", "").strip()
    if shared:
        return _positive_int_env("UPSIDE_MARTINI_BEAD_FRAME_COUNT", default)
    return default


def _truthy_env(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _write_common_table_contract_attrs(
    group: h5py.Group,
    table_family: str,
    source_object: str,
    target_object: str,
    projection_ensemble: str,
    runtime_representation: str,
    correction_layer: str = "none",
) -> None:
    group.attrs["table_generation_contract"] = "explicit_dry_martini_constituent_projection"
    group.attrs["table_family"] = table_family
    group.attrs["source_object"] = source_object
    group.attrs["target_object"] = target_object
    group.attrs["projection_ensemble"] = projection_ensemble
    group.attrs["base_energy_source"] = "dry_martini_pair_params_lj_coulomb"
    group.attrs["base_energy_helper"] = "_compute_pair_energy_and_gradient"
    group.attrs["correction_layer"] = correction_layer
    group.attrs["unit_contract"] = "native_dry_martini_nm_kjmol_e_to_upside_runtime_attrs"
    group.attrs["runtime_representation"] = runtime_representation
    group.attrs["nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    group.attrs["numerical_distance_guard_nm"] = np.float32(NUMERICAL_DISTANCE_GUARD_NM)


def _cg_lipid_pair_radial_support(
    ref_bead_positions_nm: np.ndarray,
    knot_spacing_ang: float,
) -> tuple[float, int]:
    ref = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if ref.ndim == 3:
        ref = ref.reshape(-1, 3)
    max_radius_nm = float(np.max(np.linalg.norm(ref, axis=1)))
    cutoff_nm = DRY_MARTINI_NONBONDED_CUTOFF_NM + 2.0 * max_radius_nm
    n_knot_radial = int(math.ceil(cutoff_nm * LENGTH_CONVERSION_A_PER_NM / knot_spacing_ang)) + 2
    r_max_nm = float((n_knot_radial - 2) * knot_spacing_ang / LENGTH_CONVERSION_A_PER_NM)
    return r_max_nm, n_knot_radial


def _cg_lipid_target_radial_support(
    ref_bead_positions_nm: np.ndarray,
    knot_spacing_ang: float,
) -> tuple[float, int]:
    ref = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if ref.ndim == 3:
        ref = ref.reshape(-1, 3)
    max_radius_nm = float(np.max(np.linalg.norm(ref, axis=1)))
    cutoff_nm = DRY_MARTINI_NONBONDED_CUTOFF_NM + max_radius_nm
    n_knot_radial = int(math.ceil(cutoff_nm * LENGTH_CONVERSION_A_PER_NM / knot_spacing_ang)) + 2
    r_max_nm = float((n_knot_radial - 2) * knot_spacing_ang / LENGTH_CONVERSION_A_PER_NM)
    return r_max_nm, n_knot_radial


def _ensure_cg_bonds_angles(lipids_itp_path: Path):
    global _CURRENT_CG_BONDS, _CURRENT_CG_ANGLES
    if _CURRENT_CG_BONDS is not None:
        return
    from martini_itp_reader import parse_dopc_from_itp
    dopc = parse_dopc_from_itp(lipids_itp_path)
    _CURRENT_CG_BONDS = list(dopc["bonds"])
    _CURRENT_CG_ANGLES = list(dopc["angles"])


def compute_dopc_bonded_energy(positions: np.ndarray, lipids_itp_path: Path) -> float:
    _ensure_cg_bonds_angles(lipids_itp_path)
    positions = np.asarray(positions, dtype=np.float64)
    energy = 0.0
    for i, j, r0, k in _CURRENT_CG_BONDS:
        dr = positions[i] - positions[j]
        r = float(np.linalg.norm(dr))
        energy += 0.5 * float(k) * (r - float(r0)) ** 2
    for i, j, k_idx, theta0_deg, k_ang in _CURRENT_CG_ANGLES:
        r_ij = positions[i] - positions[j]
        r_kj = positions[k_idx] - positions[j]
        d_ij = float(np.linalg.norm(r_ij))
        d_kj = float(np.linalg.norm(r_kj))
        if d_ij < 1e-8 or d_kj < 1e-8:
            continue
        cos_theta = float(np.dot(r_ij, r_kj)) / (d_ij * d_kj)
        cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
        cos_theta0 = math.cos(math.radians(float(theta0_deg)))
        energy += 0.5 * float(k_ang) * (cos_theta - cos_theta0) ** 2
    return float(energy)


def compute_dopc_intramolecular_nonbonded_energy(
    positions: np.ndarray,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
) -> float:
    _ensure_cg_bonds_angles(lipids_itp_path)
    positions = np.asarray(positions, dtype=np.float64)
    bond_pairs = {
        (min(int(i), int(j)), max(int(i), int(j)))
        for i, j, _r0, _k in (_CURRENT_CG_BONDS or [])
    }
    energy = 0.0
    for i in range(len(positions) - 1):
        for j in range(i + 1, len(positions)):
            if (i, j) in bond_pairs:
                continue
            params = pair_params.get((bead_types[i], bead_types[j]))
            if params is None:
                params = pair_params.get((bead_types[j], bead_types[i]))
            if params is None:
                raise RuntimeError(
                    f"Missing pair params for DOPC intramolecular pair ({bead_types[i]}, {bead_types[j]})"
                )
            dr = positions[j] - positions[i]
            r = float(np.linalg.norm(dr))
            if r > DRY_MARTINI_NONBONDED_CUTOFF_NM:
                continue
            eff_r = max(r, float(dist_min_nm))
            sigma = float(params["sigma_nm"])
            epsilon = float(params["epsilon_kj_mol"])
            sr = sigma / eff_r
            sr2 = sr * sr
            sr6 = sr2 * sr2 * sr2
            energy += 4.0 * epsilon * (sr6 * sr6 - sr6)
            qi = float(bead_charges[i])
            qj = float(bead_charges[j])
            if qi != 0.0 and qj != 0.0:
                energy += COULOMB_K_DRY_KJ_NM * qi * qj / eff_r
    return float(energy)


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _dot3(a: Iterable[float], b: Iterable[float]) -> float:
    return float(sum(x * y for x, y in zip(a, b)))


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
        directions.append([math.cos(theta) * radius, y, math.sin(theta) * radius])
    return directions


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


def _accumulate_sample_on_cos_grid(
    cos_theta: float,
    value: float,
    sample_grid: List[List[float]],
    sample_weight_grid: List[List[float]],
    cos_grid: List[float],
    sample_weight: float = 1.0,
) -> None:
    if len(cos_grid) == 1:
        sample_grid[0].append(value)
        sample_weight_grid[0].append(sample_weight)
        return
    step = cos_grid[1] - cos_grid[0]
    coord = (cos_theta - cos_grid[0]) / step
    if coord <= 0.0:
        sample_grid[0].append(value)
        sample_weight_grid[0].append(sample_weight)
        return
    if coord >= len(cos_grid) - 1:
        sample_grid[-1].append(value)
        sample_weight_grid[-1].append(sample_weight)
        return
    lo = int(math.floor(coord))
    hi = lo + 1
    hi_weight = coord - float(lo)
    lo_weight = 1.0 - hi_weight
    sample_grid[lo].append(value)
    sample_weight_grid[lo].append(sample_weight * lo_weight)
    sample_grid[hi].append(value)
    sample_weight_grid[hi].append(sample_weight * hi_weight)


def _factorize_one_sided_orientation(
    sampled_energy_grid: List[List[float]],
    cos_theta_grid: List[float],
) -> Tuple[List[float], List[float], List[float], float]:
    sampled = np.asarray(sampled_energy_grid, dtype=np.float64)
    radial = sampled.mean(axis=0)
    residual = sampled - radial[None, :]

    if np.allclose(residual, 0.0):
        angular_profile = np.zeros(sampled.shape[0], dtype=np.float64)
        angular_radial = np.zeros(sampled.shape[1], dtype=np.float64)
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
        reconstruction = radial[None, :] + angular_profile[:, None] * angular_radial[None, :]
        floor = sampled.min(axis=0)
        radial += np.maximum(floor[None, :] - reconstruction, 0.0).max(axis=0)
        core_mask = sampled.max(axis=0) > 1.0e6
        if np.any(core_mask):
            radial[core_mask] = sampled.min(axis=0)[core_mask]
            angular_radial[core_mask] = 0.0
        sample_min = sampled.min(axis=0)
        sample_max = sampled.max(axis=0)
        radial = np.clip(radial, sample_min, sample_max)
        angular_bound = np.maximum(np.minimum(radial - sample_min, sample_max - radial), 0.0)
        angular_radial = np.clip(angular_radial, -angular_bound, angular_bound)
        reconstruction = radial[None, :] + angular_profile[:, None] * angular_radial[None, :]
        rms_error = float(np.sqrt(np.mean((sampled - reconstruction) ** 2)))

    return (
        [float(x) for x in radial],
        [float(x) for x in angular_profile],
        [float(x) for x in angular_radial],
        rms_error,
    )


def _load_sidechain_orientation_library(
    sidechain_lib_path: Path,
) -> Dict[str, Dict[str, Any]]:
    if not sidechain_lib_path.exists():
        raise RuntimeError(f"Sidechain orientation library not found: {sidechain_lib_path}")
    with h5py.File(sidechain_lib_path, "r") as h5:
        restype_order = [
            x.decode("ascii") if isinstance(x, bytes) else str(x)
            for x in h5["restype_order"][:]
        ]
        start_stop_bead = h5["rotamer_start_stop_bead"][:]
        rotamer_center_fixed = h5["rotamer_center_fixed"][:, :6]
        prob_weights = None
        if "rotamer_prob_fixed" in h5:
            prob_weights = [float(x) for x in h5["rotamer_prob_fixed"][:].reshape(-1)]

        def _mean_over_prefix_axes(values):
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
            return [v / float(total) for v in sums]

        def _iter_index_prefix(shape):
            if not shape:
                yield ()
                return
            first, rest = shape[0], shape[1:]
            for i in range(first):
                for tail in _iter_index_prefix(rest):
                    yield (i,) + tail

        if prob_weights is None and "rotamer_prob" in h5:
            prob_weights = _mean_over_prefix_axes(h5["rotamer_prob"])

        residue_info: Dict[str, Dict[str, Any]] = {}
        for residue_index, residue_name in enumerate(restype_order):
            start, stop, n_bead = [int(x) for x in start_stop_bead[residue_index]]
            if stop <= start or n_bead <= 0:
                residue_info[residue_name] = {
                    "center_nm": [], "vector_unit": [], "weight": [],
                }
                continue
            n_row = stop - start
            if n_row % n_bead != 0:
                raise RuntimeError(
                    f"Sidechain orientation rows not divisible by n_bead for {residue_name}"
                )
            n_rot = n_row // n_bead
            block = rotamer_center_fixed[start:stop].reshape(n_rot, n_bead, 6)
            centers_nm = [
                [float(val) * ANGSTROM_TO_NM for val in row]
                for row in block[:, :, :3].mean(axis=1)
            ]
            # Per-bead positions in nm for each rotamer: (n_rot, n_bead, 3)
            bead_positions_nm = np.asarray(block[:, :, :3], dtype=np.float64) * ANGSTROM_TO_NM
            vectors = []
            for row in block[:, 0, 3:6]:
                norm = math.sqrt(sum(float(x) * float(x) for x in row))
                if norm <= 0.0:
                    raise RuntimeError(
                        f"Invalid zero sidechain vector for residue {residue_name}"
                    )
                vectors.append([float(x) / norm for x in row])
            if prob_weights is None:
                weights = [1.0 / float(n_rot)] * n_rot
            else:
                row_weights = prob_weights[start:stop:n_bead]
                total = sum(max(0.0, float(x)) for x in row_weights)
                weights = [max(0.0, float(x)) / total for x in row_weights] if total > 0 else [1.0 / float(n_rot)] * n_rot
            residue_info[residue_name] = {
                "center_nm": centers_nm,
                "vector_unit": vectors,
                "weight": weights,
                "n_rotamer": n_rot,
                "n_bead": n_bead,
                "bead_positions_nm": bead_positions_nm,
            }
    return residue_info


def _load_martini_sidechain_offsets_nm(
    martinize_path: Path,
    forcefield_name: str,
) -> Dict[str, np.ndarray]:
    """Build MARTINI sidechain bead offsets centered on the CG sidechain point.

    The Upside rotamer library provides one center/vector per sidechain rotamer.
    MARTINI can represent that sidechain by multiple beads; this helper derives
    their local offsets from the MARTINI sidechain bonded geometry.
    """
    martinize_path = Path(martinize_path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location(
        "martini_sidechain_geometry_runtime", martinize_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load martinize module from {martinize_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, forcefield_name):
        raise RuntimeError(f"Forcefield '{forcefield_name}' not found in {martinize_path}")
    ff = getattr(module, forcefield_name)()

    offsets: Dict[str, np.ndarray] = {}
    for residue in CANONICAL_RESIDUES:
        raw = ff.sidechains.get(residue, [])
        bead_types = [str(tok).strip() for tok in raw[0]] if raw else []
        bead_types = [tok for tok in bead_types if tok and tok != "D"]
        n_bead = len(bead_types)
        if n_bead == 0:
            offsets[residue] = np.zeros((0, 3), dtype=np.float64)
            continue
        if n_bead == 1:
            offsets[residue] = np.zeros((1, 3), dtype=np.float64)
            continue

        bonds = list(raw[1]) if len(raw) > 1 else []
        connectivity = getattr(ff, "connectivity", {}).get(residue, [[]])
        conn_bonds = list(connectivity[0]) if connectivity else []

        def bond_length(a: int, b: int, default: float) -> float:
            key = {int(a), int(b)}
            for idx, pair in enumerate(conn_bonds):
                if {int(pair[0]), int(pair[1])} == key and idx < len(bonds):
                    return float(bonds[idx][0])
            return float(default)

        bb_sc = bond_length(0, 1, float(bonds[0][0]) if bonds else 0.32)
        sc_sc = bond_length(1, 2, float(bonds[1][0]) if len(bonds) > 1 else 0.27)

        if n_bead == 2:
            pos = np.array(
                [
                    [0.0, 0.0, bb_sc],
                    [0.0, 0.0, bb_sc + sc_sc],
                ],
                dtype=np.float64,
            )
        elif n_bead == 3:
            half = 0.5 * sc_sc
            axial = math.sqrt(max(sc_sc * sc_sc - half * half, 0.0))
            pos = np.array(
                [
                    [0.0, 0.0, bb_sc],
                    [half, 0.0, bb_sc + axial],
                    [-half, 0.0, bb_sc + axial],
                ],
                dtype=np.float64,
            )
        elif n_bead == 4:
            half = 0.5 * sc_sc
            axial = math.sqrt(max(sc_sc * sc_sc - half * half, 0.0))
            pos = np.array(
                [
                    [0.0, 0.0, bb_sc],
                    [half, 0.0, bb_sc + axial],
                    [-half, 0.0, bb_sc + axial],
                    [0.0, 0.0, bb_sc + 2.0 * axial],
                ],
                dtype=np.float64,
            )
        else:
            raise RuntimeError(
                f"Unsupported MARTINI sidechain bead count for {residue}: {n_bead}"
            )
        offsets[residue] = pos - pos.mean(axis=0, keepdims=True)
    return offsets


def _expand_rotamer_sidechain_positions(
    orientation: Dict[str, Any],
    residue: str,
    offsets_nm: np.ndarray,
) -> List[np.ndarray]:
    centers = [np.asarray(c, dtype=np.float64) for c in orientation["center_nm"]]
    vectors = [np.asarray(v, dtype=np.float64) for v in orientation["vector_unit"]]
    out = []
    for irot, (center, vector) in enumerate(zip(centers, vectors)):
        if center.shape != (3,):
            raise RuntimeError(f"Invalid rotamer center for {residue} rotamer {irot}")
        if vector.shape != (3,) or float(np.linalg.norm(vector)) <= 1.0e-12:
            raise RuntimeError(f"Invalid sidechain vector for {residue} rotamer {irot}")
        rot = _rotation_to_align_z_np(vector)
        out.append(center[None, :] + (rot @ offsets_nm.T).T)
    return out


def _build_particles_group(
    h5: h5py.File,
    atomtypes: List[str],
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    active_atom_types: Set[str],
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
) -> None:
    active_types = [t for t in atomtypes if t in active_atom_types]
    if not active_types:
        active_types = list(atomtypes)

    charges = np.asarray([infer_charge_from_atomtype(t) for t in atomtypes], dtype=np.float32)
    type_to_charge = {t: infer_charge_from_atomtype(t) for t in atomtypes}

    # Find unique (eps, sig) pairs for active types (pair with self and cross)
    unique_eps_sig: Set[Tuple[float, float]] = set()
    seen_types = set()
    for ti in active_types:
        for tj in active_types:
            key = (ti, tj)
            if key in pair_params and key not in seen_types:
                seen_types.add(key)
                p = pair_params[key]
                eps = p["epsilon_kj_mol"] / energy_conv
                sig = p["sigma_nm"] * length_conv
                unique_eps_sig.add((eps, sig))

    eps_sig_list = sorted(unique_eps_sig)
    active_charges = sorted(set(type_to_charge[t] for t in active_types))
    qq_set = sorted(set(qi * qj for qi in active_charges for qj in active_charges))

    # Build combined LJ+Coulomb grids for each (eps, sig, qq) triple
    triples: List[Tuple[float, float, float]] = []
    grids: List[np.ndarray] = []
    coulomb_k_eup = COULOMB_K_DRY_KJ_NM * length_conv / energy_conv

    for eps, sig in eps_sig_list:
        for qq in qq_set:
            grid = np.zeros(PARTICLES_GRID_N, dtype=np.float64)
            for i in range(PARTICLES_GRID_N):
                r = PARTICLES_R_MIN_A + i * (PARTICLES_R_MAX_A - PARTICLES_R_MIN_A) / (PARTICLES_GRID_N - 1)
                r = max(r, 0.1 * sig)
                r2 = r * r
                r3 = r2 * r
                r6 = r3 * r3
                sig2 = sig * sig
                sig6 = sig2 * sig2 * sig2
                sig12 = sig6 * sig6
                lj = 4.0 * eps * (sig12 / (r6 * r6) - sig6 / r6)
                coul = coulomb_k_eup * qq / max(r, 1.0e-6) if abs(qq) > 1e-10 else 0.0
                grid[i] = lj + coul
            triples.append((eps, sig, qq))
            grids.append(grid)

    n_triples = len(triples)
    g = h5.create_group("particles")
    g.attrs["schema"] = SCHEMA_PARTICLES
    g.attrs["n_points"] = PARTICLES_GRID_N
    g.attrs["r_min_ang"] = np.float32(PARTICLES_R_MIN_A)
    g.attrs["r_max_ang"] = np.float32(PARTICLES_R_MAX_A)
    g.attrs["coulomb_k_eup"] = np.float32(coulomb_k_eup)

    eps_arr = np.asarray([t[0] for t in triples], dtype=np.float64)
    sig_arr = np.asarray([t[1] for t in triples], dtype=np.float64)
    qq_arr = np.asarray([t[2] for t in triples], dtype=np.float64)
    combined = np.zeros((n_triples, PARTICLES_GRID_N), dtype=np.float64)
    for idx, grid in enumerate(grids):
        combined[idx, :] = grid

    g.create_dataset("unique_eps_eup", data=eps_arr, dtype=np.float64)
    g.create_dataset("unique_sig_ang", data=sig_arr, dtype=np.float64)
    g.create_dataset("unique_charge_product", data=qq_arr, dtype=np.float64)
    g.create_dataset("combined_energy_grids", data=combined, dtype=np.float64)
    g.create_dataset(
        "type_order",
        data=np.asarray([np.bytes_(x) for x in atomtypes], dtype="S8"),
    )
    g.create_dataset("type_charge", data=charges, dtype=np.float32)

    print(
        f"  particles: {n_triples} combined (eps,sig,qq) triples, "
        f"{len(active_types)} active types, {PARTICLES_GRID_N} radial points"
    )


def _run_sc_task(
    residue: str,
    sidechain_bead_types: List[str],
    sidechain_bead_charges: List[float],
    target_type: str,
    target_charge: float,
    rotamer_bead_positions_nm: List[List[List[float]]],
    rotamer_weights: List[float],
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    r_values: List[float],
    direction_vectors: List[List[float]],
    cos_theta_grid: List[float],
    cb_anchor_nm: List[float],
    cb_vector_unit: List[float],
    sidechain_bead_frame_angles: List[float],
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    temperature: float = 0.0,
) -> Dict[str, Any]:
    n_rotamer = len(rotamer_bead_positions_nm)
    n_angle = len(cos_theta_grid)
    n_radial = len(r_values)
    rotamer_positions = [
        np.asarray(pos, dtype=np.float64) for pos in rotamer_bead_positions_nm
    ]
    n_bead = len(sidechain_bead_types)
    for irot, pos in enumerate(rotamer_positions):
        if pos.shape != (n_bead, 3):
            raise RuntimeError(
                f"{residue} rotamer {irot} has bead geometry shape {pos.shape}, "
                f"expected ({n_bead}, 3)"
            )

    angular_energy = [[0.0 for _ in r_values] for _ in cos_theta_grid]
    rotamer_angular_energy = [
        [[0.0 for _ in r_values] for _ in cos_theta_grid]
        for _ in range(n_rotamer)
    ]

    cb_anchor_arr = np.asarray(cb_anchor_nm, dtype=np.float64)
    cb_vector_arr = np.asarray(cb_vector_unit, dtype=np.float64)
    cb_vector_arr /= max(float(np.linalg.norm(cb_vector_arr)), 1e-12)
    bead_frame_angles = [float(x) for x in sidechain_bead_frame_angles] or [0.0]

    for ir, r_nm in enumerate(r_values):
        energy_samples = [[] for _ in cos_theta_grid]
        energy_sample_weights = [[] for _ in cos_theta_grid]
        rotamer_energy_samples = [[[] for _ in cos_theta_grid] for _ in range(n_rotamer)]
        rotamer_energy_sample_weights = [[[] for _ in cos_theta_grid] for _ in range(n_rotamer)]

        for direction in direction_vectors:
            target_pos_nm = [
                cb_anchor_nm[0] + r_nm * direction[0],
                cb_anchor_nm[1] + r_nm * direction[1],
                cb_anchor_nm[2] + r_nm * direction[2],
            ]
            cos_theta = _clamp(-_dot3(direction, cb_vector_unit), -1.0, 1.0)
            for irot, (sc_positions, rot_weight) in enumerate(
                zip(rotamer_positions, rotamer_weights)
            ):
                target = np.asarray(target_pos_nm, dtype=np.float64).reshape(1, 3)
                for sc_frame_angle in bead_frame_angles:
                    framed_sc_positions = _rotate_points_about_axis_np(
                        sc_positions, cb_vector_arr, sc_frame_angle, cb_anchor_arr
                    )
                    rot_energy, _, _ = _compute_pair_energy_and_gradient(
                        framed_sc_positions,
                        target,
                        sidechain_bead_types,
                        [target_type],
                        sidechain_bead_charges,
                        [target_charge],
                        pair_params,
                        dist_min_nm=dist_min_nm,
                    )
                    _accumulate_sample_on_cos_grid(
                        cos_theta,
                        float(rot_energy),
                        rotamer_energy_samples[irot],
                        rotamer_energy_sample_weights[irot],
                        cos_theta_grid,
                    )
                    _accumulate_sample_on_cos_grid(
                        cos_theta,
                        float(rot_energy),
                        energy_samples,
                        energy_sample_weights,
                        cos_theta_grid,
                        sample_weight=float(rot_weight),
                    )

        for ia in range(n_angle):
            if not energy_samples[ia]:
                raise RuntimeError(
                    f"Cos(theta) bin empty for {residue} target {target_type} at r={r_nm:.4f} nm"
                )
            angular_energy[ia][ir] = _boltzmann_free_energy_kj_mol(
                energy_samples[ia],
                temperature,
                energy_sample_weights[ia],
            )
            for irot in range(n_rotamer):
                if not rotamer_energy_samples[irot][ia]:
                    raise RuntimeError(
                        f"Rotamer cos(theta) bin empty for {residue} target {target_type} "
                        f"rotamer={irot} at r={r_nm:.4f} nm"
                    )
                rotamer_angular_energy[irot][ia][ir] = _boltzmann_free_energy_kj_mol(
                    rotamer_energy_samples[irot][ia],
                    temperature,
                    rotamer_energy_sample_weights[irot][ia],
                )

    radial_energy, angular_profile, angular_radial_energy, rms_error = (
        _factorize_one_sided_orientation(angular_energy, cos_theta_grid)
    )

    rotamer_radial_energy = []
    rotamer_angular_profile = []
    rotamer_angular_radial_energy = []
    rotamer_rms_error = []
    for irot in range(n_rotamer):
        rr, rp, ra, rrm = _factorize_one_sided_orientation(
            rotamer_angular_energy[irot], cos_theta_grid
        )
        rotamer_radial_energy.append([float(v) for v in rr])
        rotamer_angular_profile.append(rp)
        rotamer_angular_radial_energy.append([float(v) for v in ra])
        rotamer_rms_error.append(rrm)

    return {
        "residue": residue,
        "target_label": target_type,
        "grid_nm": r_values,
        "cos_theta_grid": cos_theta_grid,
        "radial_energy_kj_mol": radial_energy,
        "angular_energy_kj_mol": angular_radial_energy,
        "angular_profile": angular_profile,
        "rotamer_count": n_rotamer,
        "rotamer_probability_fixed": rotamer_weights,
        "rotamer_radial_energy_kj_mol": rotamer_radial_energy,
        "rotamer_angular_energy_kj_mol": rotamer_angular_radial_energy,
        "rotamer_angular_profile": rotamer_angular_profile,
        "rotamer_full_energy_kj_mol": rotamer_angular_energy,
        "factorization_rms_error": rms_error,
        "sidechain_bead_frame_count": len(bead_frame_angles),
        "azimuthal_average": (
            "energy_expectation"
            if temperature <= 0.0
            else (
                "tempered_boltzmann_free_energy"
                if abs(float(temperature) - DEFAULT_PRODUCTION_TEMP_UPSIDE) > 1e-8
                else "boltzmann_free_energy"
            )
        ),
    }


def _run_sc_task_from_dict(task: dict) -> Dict[str, Any]:
    return _run_sc_task(**task)


def _build_sc_table_group(
    h5: h5py.File,
    residue_map: Dict[str, List[str]],
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    sidechain_lib_path: Path,
    active_residue_names: List[str],
    active_target_types: List[str],
    martini_sidechain_offsets_nm: Dict[str, np.ndarray] | None = None,
    r_min_nm: float = 0.25,
    r_max_nm: float = 1.20,
    r_count: int = 96,
    direction_count: int = 24,
    cos_theta_count: int = 13,
    average_temperature: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
) -> None:
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    r_values = _linspace(r_min_nm, r_max_nm, r_count)
    direction_vectors = _fibonacci_sphere(direction_count)
    cos_theta_grid = _linspace(-1.0, 1.0, cos_theta_count)
    cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
    cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)
    sidechain_bead_frame_angles = [float(x) for x in _bead_frame_angles(_bead_frame_count("SC", 1))]

    # Determine residue-task and target-task lists
    residue_tasks = []
    for residue in active_residue_names:
        if residue not in residue_map:
            continue
        bead_types = residue_map[residue]
        if not bead_types:
            continue
        orientation = orientation_map.get(residue)
        if not orientation or len(orientation.get("bead_positions_nm", [])) == 0:
            raise RuntimeError(
                f"Missing orientation geometry for residue {residue} in {sidechain_lib_path}"
            )
        offsets_nm = (
            martini_sidechain_offsets_nm.get(residue)
            if martini_sidechain_offsets_nm is not None
            else None
        )
        if offsets_nm is not None:
            bead_positions_nm = _expand_rotamer_sidechain_positions(
                orientation, residue, np.asarray(offsets_nm, dtype=np.float64)
            )
        else:
            bead_positions_nm = orientation["bead_positions_nm"]
        bead_charges = [infer_charge_from_atomtype(bt) for bt in bead_types]
        residue_tasks.append((residue, bead_types, bead_charges, orientation, bead_positions_nm))

    if not residue_tasks:
        print("  sc_table: no active residues with sidechains, skipping")
        return

    # Collect results
    by_residue: Dict[str, Dict[str, Dict[str, Any]]] = {}
    sc_tasks = []
    for residue, bead_types, bead_charges, orientation, bead_positions_nm in residue_tasks:
        for target_type in active_target_types:
            missing = [bt for bt in bead_types if (bt, target_type) not in pair_params]
            if missing:
                raise RuntimeError(
                    f"Missing nonbond param for residue {residue} target {target_type}"
                )
            sc_tasks.append(
                {
                    "residue": residue,
                    "sidechain_bead_types": list(bead_types),
                    "sidechain_bead_charges": list(bead_charges),
                    "target_type": target_type,
                    "target_charge": infer_charge_from_atomtype(target_type),
                    "rotamer_bead_positions_nm": bead_positions_nm,
                    "rotamer_weights": orientation["weight"],
                    "pair_params": pair_params,
                    "r_values": r_values,
                    "direction_vectors": direction_vectors,
                    "cos_theta_grid": cos_theta_grid,
                    "cb_anchor_nm": cb_anchor_nm,
                    "cb_vector_unit": cb_vector_unit,
                    "sidechain_bead_frame_angles": sidechain_bead_frame_angles,
                    "dist_min_nm": NUMERICAL_DISTANCE_GUARD_NM,
                    "temperature": float(average_temperature),
                }
            )

    print(
        f"  sc_table direction-vector sampling: {direction_count} target directions, "
        f"{len(sidechain_bead_frame_angles)} SC bead-frame sample(s)"
    )
    for result in _parallel_map_ordered("SC-particle table", _run_sc_task_from_dict, sc_tasks):
        residue = result["residue"]
        target_type = result["target_label"]
        if residue not in by_residue:
            by_residue[residue] = {}
        by_residue[residue][target_type] = result

    residues = sorted(by_residue)
    targets = list(active_target_types)
    first = by_residue[residues[0]][targets[0]]
    ref_grid_nm = np.asarray(first["grid_nm"], dtype=np.float32)
    ref_cos_grid = np.asarray(first["cos_theta_grid"], dtype=np.float32)
    n_r, n_cos = ref_grid_nm.size, ref_cos_grid.size
    n_t = len(targets)
    max_rot = max(int(by_residue[r][targets[0]]["rotamer_count"]) for r in residues)

    rad_e = np.zeros((len(residues), n_t, n_r), dtype=np.float32)
    ang_e = np.zeros((len(residues), n_t, n_r), dtype=np.float32)
    ang_p = np.zeros((len(residues), n_t, n_cos), dtype=np.float32)
    rc = np.zeros((len(residues),), dtype=np.float32)
    rpf = np.zeros((len(residues), max_rot), dtype=np.float32)
    r_rad = np.zeros((len(residues), max_rot, n_t, n_r), dtype=np.float32)
    r_ang = np.zeros((len(residues), max_rot, n_t, n_r), dtype=np.float32)
    r_prof = np.zeros((len(residues), max_rot, n_t, n_cos), dtype=np.float32)
    r_full = np.zeros((len(residues), max_rot, n_t, n_r, n_cos), dtype=np.float32)

    for ri, residue in enumerate(residues):
        for ti, target in enumerate(targets):
            e = by_residue[residue][target]
            n_rot = int(e["rotamer_count"])
            if ti == 0:
                rc[ri] = float(n_rot)
                rpf[ri, :n_rot] = np.asarray(e["rotamer_probability_fixed"], dtype=np.float32)
            rad_e[ri, ti, :] = np.asarray(e["radial_energy_kj_mol"], dtype=np.float32)
            ang_e[ri, ti, :] = np.asarray(e["angular_energy_kj_mol"], dtype=np.float32)
            ang_p[ri, ti, :] = np.asarray(e["angular_profile"], dtype=np.float32)
            r_rad[ri, :n_rot, ti, :] = np.asarray(e["rotamer_radial_energy_kj_mol"], dtype=np.float32)
            r_ang[ri, :n_rot, ti, :] = np.asarray(e["rotamer_angular_energy_kj_mol"], dtype=np.float32)
            r_prof[ri, :n_rot, ti, :] = np.asarray(e["rotamer_angular_profile"], dtype=np.float32)
            r_full[ri, :n_rot, ti, :, :] = np.asarray(
                e["rotamer_full_energy_kj_mol"], dtype=np.float32
            ).transpose(0, 2, 1)

    g = h5.create_group("sc_table")
    _write_common_table_contract_attrs(
        g,
        table_family="SC-particle",
        source_object="sidechain_rotamer_bead_ensemble",
        target_object="dry_martini_particle_type",
        projection_ensemble="sidechain_rotamers_x_sidechain_bead_frame_x_target_direction_grid",
        runtime_representation="full_radial_angular_grid",
    )
    g.attrs["schema"] = SCHEMA_SC
    g.attrs["forcefield_name"] = "martini22"
    g.attrs["sample_dist_min_nm"] = np.float32(NUMERICAL_DISTANCE_GUARD_NM)
    g.attrs["sidechain_bead_frame_count"] = np.int32(len(sidechain_bead_frame_angles))
    g.attrs["orientation_sampling"] = "target_direction_vector_grid"
    g.attrs["azimuthal_average"] = first.get("azimuthal_average", "")
    g.attrs["azimuthal_average_temperature_upside"] = np.float32(average_temperature)
    g.attrs["excluded_area_source"] = "none_direct_dry_martini_sc_particle_table"
    g.attrs["attractive_control_source"] = "retained_direct_dry_martini_sc_particle_table"
    g.attrs["isotropic_background_source"] = "none_direct_dry_martini_sc_particle_table"
    g.attrs["relaxation"] = "rigid_rotated_geometry"
    g.attrs["runtime_representation"] = "full_radial_angular_grid"
    g.attrs["spline_control_quantity"] = "direct_rotamer_free_energy_kj_mol"
    g.create_dataset("restype_order", data=np.asarray([np.bytes_(x) for x in residues], dtype="S4"))
    g.create_dataset("target_order", data=np.asarray([np.bytes_(x) for x in targets], dtype="S8"))
    g.create_dataset("grid_nm", data=ref_grid_nm, dtype=np.float32)
    g.create_dataset("cos_theta_grid", data=ref_cos_grid, dtype=np.float32)
    g.create_dataset("rotamer_count", data=rc, dtype=np.float32)
    g.create_dataset("rotamer_probability_fixed", data=rpf, dtype=np.float32)
    for name, arr in [
        ("radial_energy_kj_mol", rad_e),
        ("angular_energy_kj_mol", ang_e),
        ("angular_profile", ang_p),
        ("rotamer_radial_energy_kj_mol", r_rad),
        ("rotamer_angular_energy_kj_mol", r_ang),
        ("rotamer_angular_profile", r_prof),
        ("rotamer_full_energy_kj_mol", r_full),
    ]:
        g.create_dataset(name, data=arr, dtype=np.float32)

    print(
        f"  sc_table: {len(residues)} residues x {len(targets)} targets, "
        f"{n_cos} angular x {n_r} radial points, max {max_rot} rotamers"
    )


# The C++ deBoor_value_and_deriv in spline.h implements uniform cubic B-spline
# evaluation with integer knot spacing. It does NOT use a custom knot vector.
# We replicate that algorithm here for fitting control points that the C++ can
# consume directly.


def _deBoor_uniform_basis_weights(t: float, n_control: int) -> np.ndarray:
    """Return (n_control,) array of basis weights at parameter t.

    Replicates C++ deBoor_value_and_deriv / uniform_deBoor_algorithm:
      x_bin = int(t), excess = t - x_bin
      Uses 4 consecutive control points starting at x_bin-1.
      The spline value is the dot product of weights with control points.

    For integer t (excess=0): weights = [1/6, 2/3, 1/6, 0] at positions
    [x_bin-1, x_bin, x_bin+1, x_bin+2].

    Note: The C++ code reads 4 consecutive floats starting at x_bin-1. At the
    upper boundary (x_bin+2 >= n_control), it reads into the next array in
    memory (Ang1 to Ang2). The 4th coefficient weight is zero for integer t, so
    this is harmless. We allow x_bin up to n_control-1 (one-past-end read).
    """
    w = np.zeros(n_control, dtype=np.float64)
    x_bin = int(t)
    if x_bin < 1 or x_bin >= n_control:
        return w
    excess = t - x_bin

    for di in range(4):
        ci = x_bin - 1 + di
        if ci < 0 or ci >= n_control:
            continue

        c_arr = np.zeros(4)
        c_arr[di] = 1.0

        a11 = (excess + 2.0) / 3.0
        a12 = (excess + 1.0) / 3.0
        a13 = excess / 3.0

        c11 = (1.0 - a11) * c_arr[0] + a11 * c_arr[1]
        c12 = (1.0 - a12) * c_arr[1] + a12 * c_arr[2]
        c13 = (1.0 - a13) * c_arr[2] + a13 * c_arr[3]

        a22 = (excess + 1.0) / 2.0
        a23 = excess / 2.0

        c22 = (1.0 - a22) * c11 + a22 * c12
        c23 = (1.0 - a23) * c12 + a23 * c13

        c33 = (1.0 - excess) * c22 + excess * c23

        w[ci] = float(c33)

    return w


def _deBoor_uniform_basis_matrix(t_samples: np.ndarray, n_control: int) -> np.ndarray:
    """Compute (n_samples, n_control) basis matrix for uniform deBoor evaluation."""
    basis = np.zeros((len(t_samples), n_control), dtype=np.float64)
    for si, t in enumerate(t_samples):
        basis[si, :] = _deBoor_uniform_basis_weights(float(t), n_control)
    return basis


def _fit_angular_bspline(
    t_samples: np.ndarray,
    ang_values: np.ndarray,
    n_control: int = 15,
    smooth: float = 0.01,
) -> np.ndarray:
    """Fit control points for uniform deBoor angular B-spline.

    Uses the exact same evaluation algorithm as C++ spline.h:deBoor_value_and_deriv,
    then solves with Tikhonov regularization on second differences.

    The basis matrix has full column rank (every control point contributes to
    at least one sample) so the fit is well-conditioned even with tiny smooth.
    """
    basis = _deBoor_uniform_basis_matrix(t_samples, n_control)

    D2 = np.zeros((n_control - 2, n_control), dtype=np.float64)
    for i in range(n_control - 2):
        D2[i, i] = 1.0
        D2[i, i + 1] = -2.0
        D2[i, i + 2] = 1.0

    A = basis.T @ basis + float(smooth) * D2.T @ D2
    b = basis.T @ ang_values
    control = np.linalg.solve(A, b)
    return control


def _second_difference_matrix(n_control: int) -> np.ndarray:
    d2 = np.zeros((n_control - 2, n_control), dtype=np.float64)
    for i in range(n_control - 2):
        d2[i, i] = 1.0
        d2[i, i + 1] = -2.0
        d2[i, i + 2] = 1.0
    return d2


def _cubic_bspline_basis_one(t: float, knots: np.ndarray, i: int) -> float:
    """Evaluate cubic B-spline basis function N_{i,3}(t) using Cox-de Boor recurrence."""
    n0 = 1.0 if knots[i] <= t < knots[i + 1] else 0.0
    if n0 == 0.0 and (t < knots[i] or t >= knots[i + 4]):
        return 0.0

    def _n(i: int, k: int) -> float:
        if k == 0:
            return 1.0 if knots[i] <= t < knots[i + 1] else 0.0
        c1 = 0.0
        denom1 = knots[i + k] - knots[i]
        if denom1 > 1e-15:
            c1 = (t - knots[i]) / denom1 * _n(i, k - 1)
        c2 = 0.0
        denom2 = knots[i + k + 1] - knots[i + 1]
        if denom2 > 1e-15:
            c2 = (knots[i + k + 1] - t) / denom2 * _n(i + 1, k - 1)
        return c1 + c2

    return _n(i, 3)


def _cubic_bspline_basis_values(t_samples: np.ndarray, knot_vector: np.ndarray) -> np.ndarray:
    """Compute (n_samples, n_control) matrix of cubic B-spline basis values."""
    n_control = len(knot_vector) - 4
    n_samples = len(t_samples)
    basis = np.zeros((n_samples, n_control), dtype=np.float64)
    for si, t in enumerate(t_samples):
        for ci in range(n_control):
            basis[si, ci] = _cubic_bspline_basis_one(float(t), knot_vector, ci)
    return basis


def _fit_radial_bspline(
    t_samples: np.ndarray,
    y_samples: np.ndarray,
    knot_vector: np.ndarray,
    smooth: float = 0.0,
) -> np.ndarray:
    """Fit control points for Upside's clamped deBoor radial evaluator."""
    n_control = len(knot_vector) - 4
    basis = np.zeros((len(t_samples), n_control), dtype=np.float64)
    for si, t in enumerate(t_samples):
        t = float(t)
        if t <= 1.0:
            basis[si, 0:3] = (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
        elif t >= n_control - 2:
            basis[si, n_control - 3:n_control] = (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
        else:
            basis[si, :] = _deBoor_uniform_basis_weights(t, n_control)

    if smooth > 0.0:
        D2 = np.zeros((n_control - 2, n_control), dtype=np.float64)
        for i in range(n_control - 2):
            D2[i, i] = 1.0
            D2[i, i + 1] = -2.0
            D2[i, i + 2] = 1.0
        A = basis.T @ basis + float(smooth) * D2.T @ D2
        b = basis.T @ y_samples
        control = np.linalg.solve(A, b)
    else:
        control, _, _, _ = np.linalg.lstsq(basis, y_samples, rcond=None)

    return control


def _radial_deboor_basis_matrix(t_samples: np.ndarray, n_control: int) -> np.ndarray:
    """Compute the clamped radial spline basis used by Upside's runtime evaluator."""
    basis = np.zeros((len(t_samples), n_control), dtype=np.float64)
    for si, t in enumerate(t_samples):
        t = float(t)
        if t <= 1.0:
            basis[si, 0:3] = (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
        elif t >= n_control - 2:
            basis[si, n_control - 3:n_control] = (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
        else:
            basis[si, :] = _deBoor_uniform_basis_weights(t, n_control)
    return basis


def _fit_radial_angular_tensor_bspline(
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    values_kj_mol: np.ndarray,
    n_knot_radial: int,
    n_knot_angular: int,
    knot_spacing_ang: float,
    energy_conversion: float = ENERGY_CONVERSION_KJ_PER_EUP,
    smooth: float = 0.01,
) -> np.ndarray:
    """Fit tensor-product controls for E(r, cos theta) in Upside units."""
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    values = np.asarray(values_kj_mol, dtype=np.float64)
    if values.shape != (r_values_nm.size, cos_theta_grid.size):
        raise ValueError(
            "radial/angular grid shape mismatch: "
            f"{values.shape} vs {(r_values_nm.size, cos_theta_grid.size)}"
        )

    rad_knot_vector = np.zeros(n_knot_radial + 4, dtype=np.float64)
    rad_knot_vector[4:-4] = np.arange(1, n_knot_radial - 3, dtype=np.float64)
    rad_knot_vector[-4:] = rad_knot_vector[-5]
    t_radial = r_values_nm * LENGTH_CONVERSION_A_PER_NM / float(knot_spacing_ang)
    t_angular = (cos_theta_grid + 1.0) * (float(n_knot_angular - 3) * 0.5) + 1.0

    radial_basis = _radial_deboor_basis_matrix(t_radial, n_knot_radial)
    radial_rhs = values / float(energy_conversion)
    if smooth > 0.0:
        radial_d2 = _second_difference_matrix(n_knot_radial)
        radial_op = radial_basis.T @ radial_basis + float(smooth) * radial_d2.T @ radial_d2
        radial_controls = np.linalg.solve(radial_op, radial_basis.T @ radial_rhs)
    else:
        radial_controls, _, _, _ = np.linalg.lstsq(radial_basis, radial_rhs, rcond=None)

    angular_basis = _deBoor_uniform_basis_matrix(t_angular, n_knot_angular)
    angular_rhs = radial_controls.T
    if smooth > 0.0:
        angular_d2 = _second_difference_matrix(n_knot_angular)
        angular_op = angular_basis.T @ angular_basis + float(smooth) * angular_d2.T @ angular_d2
        controls = np.linalg.solve(angular_op, angular_basis.T @ angular_rhs).T
    else:
        controls, _, _, _ = np.linalg.lstsq(angular_basis, angular_rhs, rcond=None)
        controls = controls.T
    return controls


def _evaluate_radial_angular_tensor_controls(
    control_tensor: np.ndarray,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    knot_spacing_ang: float,
) -> np.ndarray:
    control_tensor = np.asarray(control_tensor, dtype=np.float64)
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    if control_tensor.ndim != 2:
        raise ValueError("radial/angular control tensor must have shape (n_radial, n_angular)")
    n_knot_radial, n_knot_angular = control_tensor.shape
    t_radial = r_values_nm * LENGTH_CONVERSION_A_PER_NM / float(knot_spacing_ang)
    t_angular = (cos_theta_grid + 1.0) * (float(n_knot_angular - 3) * 0.5) + 1.0
    radial_basis = _radial_deboor_basis_matrix(t_radial, n_knot_radial)
    angular_basis = _deBoor_uniform_basis_matrix(t_angular, n_knot_angular)
    radial_eval = radial_basis @ control_tensor
    return radial_eval @ angular_basis.T


def _fit_radial_angular_angular_tensor_bspline(
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    values_kj_mol: np.ndarray,
    n_knot_radial: int,
    n_knot_angular: int,
    knot_spacing_ang: float,
    energy_conversion: float = ENERGY_CONVERSION_KJ_PER_EUP,
    smooth: float = 0.01,
) -> np.ndarray:
    """Fit tensor-product controls for E(r, cos theta1, cos theta2) in Upside units."""
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    values = np.asarray(values_kj_mol, dtype=np.float64)
    expected = (r_values_nm.size, cos_theta_grid.size, cos_theta_grid.size)
    if values.shape != expected:
        raise ValueError(f"radial/angular/angular grid shape mismatch: {values.shape} vs {expected}")

    rad_knot_vector = np.zeros(n_knot_radial + 4, dtype=np.float64)
    rad_knot_vector[4:-4] = np.arange(1, n_knot_radial - 3, dtype=np.float64)
    rad_knot_vector[-4:] = rad_knot_vector[-5]
    t_radial = r_values_nm * LENGTH_CONVERSION_A_PER_NM / float(knot_spacing_ang)
    t_angular = (cos_theta_grid + 1.0) * (float(n_knot_angular - 3) * 0.5) + 1.0

    radial_basis = _radial_deboor_basis_matrix(t_radial, n_knot_radial)
    radial_rhs = values.reshape(len(t_radial), -1) / float(energy_conversion)
    if smooth > 0.0:
        radial_d2 = _second_difference_matrix(n_knot_radial)
        radial_op = radial_basis.T @ radial_basis + float(smooth) * radial_d2.T @ radial_d2
        radial_controls = np.linalg.solve(radial_op, radial_basis.T @ radial_rhs)
    else:
        radial_controls, _, _, _ = np.linalg.lstsq(radial_basis, radial_rhs, rcond=None)
    radial_controls = radial_controls.reshape(
        n_knot_radial,
        cos_theta_grid.size,
        cos_theta_grid.size,
    )

    angular_basis = _deBoor_uniform_basis_matrix(t_angular, n_knot_angular)
    if smooth > 0.0:
        angular_d2 = _second_difference_matrix(n_knot_angular)
        angular_op = angular_basis.T @ angular_basis + float(smooth) * angular_d2.T @ angular_d2
        solve_angular = lambda rhs: np.linalg.solve(angular_op, angular_basis.T @ rhs)
    else:
        solve_angular = lambda rhs: np.linalg.lstsq(angular_basis, rhs, rcond=None)[0]

    tmp = np.zeros((n_knot_radial, n_knot_angular, cos_theta_grid.size), dtype=np.float64)
    for ir in range(n_knot_radial):
        tmp[ir, :, :] = solve_angular(radial_controls[ir, :, :])

    controls = np.zeros((n_knot_radial, n_knot_angular, n_knot_angular), dtype=np.float64)
    for ir in range(n_knot_radial):
        controls[ir, :, :] = solve_angular(tmp[ir, :, :].T).T
    return controls


def _evaluate_radial_angular_angular_tensor_controls(
    control_tensor: np.ndarray,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    knot_spacing_ang: float,
) -> np.ndarray:
    """Evaluate a stored SC full-tensor control lattice on the sampled table grid."""
    controls = np.asarray(control_tensor, dtype=np.float64)
    if controls.ndim != 3:
        raise ValueError(f"expected 3D SC control tensor, got shape {controls.shape}")

    n_knot_radial, n_knot_ang1, n_knot_ang2 = controls.shape
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    t_radial = r_values_nm * LENGTH_CONVERSION_A_PER_NM / float(knot_spacing_ang)
    t_ang1 = (cos_theta_grid + 1.0) * (float(n_knot_ang1 - 3) * 0.5) + 1.0
    t_ang2 = (cos_theta_grid + 1.0) * (float(n_knot_ang2 - 3) * 0.5) + 1.0

    radial_basis = _radial_deboor_basis_matrix(t_radial, n_knot_radial)
    angular_basis1 = _deBoor_uniform_basis_matrix(t_ang1, n_knot_ang1)
    angular_basis2 = _deBoor_uniform_basis_matrix(t_ang2, n_knot_ang2)
    return np.einsum(
        "ri,aj,bk,ijk->rab",
        radial_basis,
        angular_basis1,
        angular_basis2,
        controls,
        optimize=True,
    )


def _fit_angular_angular_tensor_bspline(
    cos_theta_grid: np.ndarray,
    values_kj_mol: np.ndarray,
    n_knot_angular: int,
    energy_conversion: float = ENERGY_CONVERSION_KJ_PER_EUP,
    smooth: float = 0.01,
) -> np.ndarray:
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    values = np.asarray(values_kj_mol, dtype=np.float64)
    expected = (cos_theta_grid.size, cos_theta_grid.size)
    if values.shape != expected:
        raise ValueError(f"angular/angular grid shape mismatch: {values.shape} vs {expected}")

    t_angular = (cos_theta_grid + 1.0) * (float(n_knot_angular - 3) * 0.5) + 1.0
    tmp = np.zeros((n_knot_angular, cos_theta_grid.size), dtype=np.float64)
    for ia2 in range(cos_theta_grid.size):
        tmp[:, ia2] = _fit_angular_bspline(
            t_angular,
            values[:, ia2] / float(energy_conversion),
            n_control=n_knot_angular,
            smooth=smooth,
        )

    controls = np.zeros((n_knot_angular, n_knot_angular), dtype=np.float64)
    for ia1 in range(n_knot_angular):
        controls[ia1, :] = _fit_angular_bspline(
            t_angular,
            tmp[ia1, :],
            n_control=n_knot_angular,
            smooth=smooth,
        )
    return controls


_CG_DERIVED_NUMERIC_ATTRS = (
    "contact_nm",
    "contact_ang",
    "max_sigma_nm",
    "orientation_length_ang",
    "orientation_mass_g_mol",
    "orientation_bond_fc_eup_a2",
    "orientation_projected_bond_fc_eup_a2",
    "orientation_carrier_bond_fc_factor",
    "transverse_inertia_g_mol_a2",
    "head_tail_span_ang",
    "tail_projection_ang",
    "max_axis_radius_ang",
    "max_perp_radius_ang",
    "energy_conversion_kj_per_eup",
    "length_conversion_ang_per_nm",
)

_CG_DERIVED_STRING_ATTRS = (
    "mass_source",
    "contact_source",
    "orientation_length_source",
    "orientation_mass_source",
    "orientation_bond_fc_source",
    "orientation_projected_bond_fc_source",
)


def _write_h5_atomically(output_path: Path, writer: Callable[[h5py.File], None]) -> None:
    """Write an HDF5 file through a sibling temp file, then atomically replace."""
    output_path = Path(output_path).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(f".{output_path.name}.tmp.{os.getpid()}")
    try:
        if tmp_path.exists():
            tmp_path.unlink()
        with h5py.File(tmp_path, "w") as h5:
            writer(h5)
        os.replace(tmp_path, output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _write_cg_derived_attrs(group, derived_params: dict) -> None:
    group.attrs["derivation_schema"] = "dry_martini_dopc_derived"
    for attr_name in _CG_DERIVED_NUMERIC_ATTRS:
        if attr_name in derived_params:
            group.attrs[attr_name] = np.float32(derived_params[attr_name])
    for attr_name in _CG_DERIVED_STRING_ATTRS:
        if attr_name in derived_params:
            group.attrs[attr_name] = derived_params[attr_name]

def _compute_pair_energy_and_gradient(
    pos1: np.ndarray,
    pos2: np.ndarray,
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    cutoff_nm: float | None = DRY_MARTINI_NONBONDED_CUTOFF_NM,
) -> tuple:
    """Compute non-bonded (LJ+Coulomb) energy and per-bead gradients for a pair of lipids.

    pos1, pos2: (14, 3) float64 arrays in nm.

    `dist_min_nm` is a numerical singularity guard, not a model parameter.
    Production table builds pass `NUMERICAL_DISTANCE_GUARD_NM`.

    Returns (total_energy_kj_mol, grad1, grad2) where grad1, grad2 are (14, 3).
    """
    n1, n2 = pos1.shape[0], pos2.shape[0]
    total_energy = 0.0
    grad1 = np.zeros_like(pos1)
    grad2 = np.zeros_like(pos2)

    for i in range(n1):
        for j in range(n2):
            dx = pos2[j] - pos1[i]
            dist_sq = float(np.dot(dx, dx))
            dist = float(np.sqrt(dist_sq))
            if cutoff_nm is not None and dist > float(cutoff_nm):
                continue
            eff_dist = max(dist, dist_min_nm)

            key = (bead_types1[i], bead_types2[j])
            params = pair_params.get(key)
            if params is None:
                params = pair_params.get((bead_types2[j], bead_types1[i]))
            if params is None:
                raise RuntimeError(
                    f"Missing pair params for ({bead_types1[i]}, {bead_types2[j]})"
                )

            sigma_nm = params["sigma_nm"]
            epsilon_kj = params["epsilon_kj_mol"]

            sr = sigma_nm / eff_dist
            sr2 = sr * sr
            sr6 = sr2 * sr2 * sr2
            lj = 4.0 * epsilon_kj * (sr6 * sr6 - sr6)
            total_energy += lj

            if dist > 1e-10:
                dlj_dr = 4.0 * epsilon_kj * (-12.0 * sr6 * sr6 / eff_dist + 6.0 * sr6 / eff_dist)
                unit = dx / dist
                grad1[i] -= dlj_dr * unit
                grad2[j] += dlj_dr * unit

            q1 = bead_charges1[i]
            q2 = bead_charges2[j]
            if q1 and q2 and dist > 1e-10:
                coul_eff = max(dist, dist_min_nm)
                coul = COULOMB_K_DRY_KJ_NM * q1 * q2 / coul_eff
                total_energy += coul
                dcoul_dr = -COULOMB_K_DRY_KJ_NM * q1 * q2 / (coul_eff * coul_eff)
                unit = dx / dist
                grad1[i] -= dcoul_dr * unit
                grad2[j] += dcoul_dr * unit

    return total_energy, grad1, grad2


def _compute_cg_pair_energy(
    r_nm: float,
    dir1: np.ndarray,
    dir2: np.ndarray,
    frame_angle1: float,
    frame_angle2: float,
    ref_bead1: np.ndarray,
    ref_bead2: np.ndarray,
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
) -> float:
    """Compute LJ + Coulomb energy between two CG lipid particles at given geometry.

    Each CG lipid's beads are placed at COM_i + R(dir_i) * ref_bead_i[k].
    """
    R1 = _rotation_to_align_z_np(dir1) @ _rotation_about_axis_np(
        np.array([0.0, 0.0, 1.0], dtype=np.float64), frame_angle1
    )
    R2 = _rotation_to_align_z_np(dir2) @ _rotation_about_axis_np(
        np.array([0.0, 0.0, 1.0], dtype=np.float64), frame_angle2
    )

    com1 = np.zeros(3)
    com2 = np.array([r_nm, 0.0, 0.0])

    pos1 = com1[None, :] + (R1 @ ref_bead1.T).T
    pos2 = com2[None, :] + (R2 @ ref_bead2.T).T

    total_energy, _, _ = _compute_pair_energy_and_gradient(
        pos1, pos2, bead_types1, bead_types2,
        bead_charges1, bead_charges2, pair_params,
        dist_min_nm=dist_min_nm,
        cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
    )
    return total_energy


_PAIR_RELAX_TAIL_INDICES = np.asarray([4, 5, 6, 7, 8, 9, 10, 11, 12, 13], dtype=np.int32)
_PAIR_RELAX_TAIL_CHAINS = (
    np.asarray([4, 5, 6, 7, 8], dtype=np.int32),
    np.asarray([9, 10, 11, 12, 13], dtype=np.int32),
)


def _compute_dopc_internal_energy(
    positions_nm: np.ndarray,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
) -> float:
    return (
        compute_dopc_bonded_energy(positions_nm, lipids_itp_path)
        + compute_dopc_intramolecular_nonbonded_energy(
            positions_nm,
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            lipids_itp_path=lipids_itp_path,
        )
    )


def _recenter_bead_positions_to_target_com(
    positions_nm: np.ndarray,
    target_com_nm: np.ndarray,
) -> np.ndarray:
    shift = np.asarray(target_com_nm, dtype=np.float64) - np.mean(positions_nm, axis=0)
    return np.asarray(positions_nm, dtype=np.float64) + shift[None, :]


def _orient_canonical_lipid_reference(
    ref_bead_nm: np.ndarray,
    axis: np.ndarray,
    target_com_nm: np.ndarray,
    frame_angle: float = 0.0,
) -> np.ndarray:
    ref = np.asarray(ref_bead_nm, dtype=np.float64)
    rot = _rotation_to_align_z_np(np.asarray(axis, dtype=np.float64)) @ _rotation_about_axis_np(
        np.array([0.0, 0.0, 1.0], dtype=np.float64),
        float(frame_angle),
    )
    oriented = (rot @ ref.T).T
    return _recenter_bead_positions_to_target_com(oriented, target_com_nm)


def _tail_axial_compaction_weights(ref_bead_nm: np.ndarray) -> np.ndarray:
    ref_bead_nm = np.asarray(ref_bead_nm, dtype=np.float64)
    weights = np.zeros(ref_bead_nm.shape[0], dtype=np.float64)
    for chain in _PAIR_RELAX_TAIL_CHAINS:
        z = ref_bead_nm[chain, 2]
        z0 = float(z[0])
        z1 = float(np.max(z))
        span = max(z1 - z0, 1.0e-6)
        weights[chain] = np.clip((z - z0) / span, 0.0, 1.0)
    return weights


def _tail_radial_unit_vectors(
    ref_bead_nm: np.ndarray,
    axis: np.ndarray,
) -> tuple[np.ndarray, float]:
    positions = np.asarray(ref_bead_nm, dtype=np.float64)
    unit_axis = np.asarray(axis, dtype=np.float64)
    unit_axis /= max(float(np.linalg.norm(unit_axis)), 1.0e-12)
    origin = np.asarray(positions[0], dtype=np.float64)
    rel = positions - origin[None, :]
    radial = rel - np.outer(rel @ unit_axis, unit_axis)
    radial_norm = np.linalg.norm(radial, axis=1)
    radial_unit = np.zeros_like(radial)
    valid = radial_norm > 1.0e-12
    radial_unit[valid] = radial[valid] / radial_norm[valid, None]
    max_radial = float(np.max(radial_norm[_PAIR_RELAX_TAIL_INDICES]))
    return radial_unit, max_radial


def _apply_tail_axial_compaction(
    ref_bead_nm: np.ndarray,
    axis: np.ndarray,
    compaction_nm: np.ndarray | float,
    tail_weights: np.ndarray,
) -> np.ndarray:
    positions = np.asarray(ref_bead_nm, dtype=np.float64).copy()
    delta = np.asarray(compaction_nm, dtype=np.float64)
    if delta.ndim == 0:
        delta = np.full(len(_PAIR_RELAX_TAIL_CHAINS), float(delta), dtype=np.float64)
    if delta.shape != (len(_PAIR_RELAX_TAIL_CHAINS),):
        raise ValueError("tail axial compaction expects one scalar per tail chain")
    axis = np.asarray(axis, dtype=np.float64)
    for chain_idx, chain in enumerate(_PAIR_RELAX_TAIL_CHAINS):
        if delta[chain_idx] <= 0.0:
            continue
        positions[chain] -= (
            float(delta[chain_idx])
            * tail_weights[chain, None]
            * axis[None, :]
        )
    return positions


def _apply_tail_axial_shift(
    ref_bead_nm: np.ndarray,
    axis: np.ndarray,
    shift_nm: np.ndarray | float,
    tail_weights: np.ndarray,
) -> np.ndarray:
    positions = np.asarray(ref_bead_nm, dtype=np.float64).copy()
    delta = np.asarray(shift_nm, dtype=np.float64)
    if delta.ndim == 0:
        delta = np.full(len(_PAIR_RELAX_TAIL_CHAINS), float(delta), dtype=np.float64)
    if delta.shape != (len(_PAIR_RELAX_TAIL_CHAINS),):
        raise ValueError("tail axial shift expects one scalar per tail chain")
    axis = np.asarray(axis, dtype=np.float64)
    for chain_idx, chain in enumerate(_PAIR_RELAX_TAIL_CHAINS):
        if abs(float(delta[chain_idx])) <= 1.0e-12:
            continue
        positions[chain] -= (
            float(delta[chain_idx])
            * tail_weights[chain, None]
            * axis[None, :]
        )
    return positions


def _apply_tail_compaction_modes(
    ref_bead_nm: np.ndarray,
    axis: np.ndarray,
    axial_compaction_nm: np.ndarray | float,
    radial_compaction_nm: float,
    tail_weights: np.ndarray,
    radial_unit: np.ndarray,
) -> np.ndarray:
    positions = _apply_tail_axial_compaction(
        ref_bead_nm,
        axis,
        axial_compaction_nm,
        tail_weights,
    )
    radial_delta = max(float(radial_compaction_nm), 0.0)
    if radial_delta <= 0.0:
        return positions
    positions = np.asarray(positions, dtype=np.float64).copy()
    positions[_PAIR_RELAX_TAIL_INDICES] -= (
        radial_delta
        * tail_weights[_PAIR_RELAX_TAIL_INDICES, None]
        * np.asarray(radial_unit[_PAIR_RELAX_TAIL_INDICES], dtype=np.float64)
    )
    return positions


def _pair_conditioned_tail_relaxation_effective_energy(
    r_nm: float,
    dir1: np.ndarray,
    dir2: np.ndarray,
    frame_angle1: float,
    frame_angle2: float,
    ref_bead1: np.ndarray,
    ref_bead2: np.ndarray,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    step_schedule_nm: tuple[float, ...] = (0.20, 0.10, 0.05, 0.025),
    max_compaction_nm: float = 0.9,
) -> tuple[float, float]:
    dir1 = np.asarray(dir1, dtype=np.float64)
    dir2 = np.asarray(dir2, dtype=np.float64)
    ref_bead1 = np.asarray(ref_bead1, dtype=np.float64)
    ref_bead2 = np.asarray(ref_bead2, dtype=np.float64)
    r_axis = np.array([float(r_nm), 0.0, 0.0], dtype=np.float64)
    target_com1 = np.zeros(3, dtype=np.float64)
    target_com2 = r_axis.copy()

    pos1 = _orient_canonical_lipid_reference(ref_bead1, dir1, target_com1, frame_angle1)
    pos2 = _orient_canonical_lipid_reference(ref_bead2, dir2, target_com2, frame_angle2)
    tail_weights1 = _tail_axial_compaction_weights(ref_bead1)
    tail_weights2 = _tail_axial_compaction_weights(ref_bead2)
    radial_unit1, max_radial1 = _tail_radial_unit_vectors(pos1, dir1)
    radial_unit2, max_radial2 = _tail_radial_unit_vectors(pos2, dir2)

    intra0 = (
        _compute_dopc_internal_energy(pos1, bead_types, bead_charges, pair_params, lipids_itp_path)
        + _compute_dopc_internal_energy(
            pos2 - r_axis[None, :],
            bead_types,
            bead_charges,
            pair_params,
            lipids_itp_path,
        )
    )
    inter0, _, _ = _compute_pair_energy_and_gradient(
        pos1,
        pos2,
        bead_types,
        bead_types,
        bead_charges,
        bead_charges,
        pair_params,
        dist_min_nm=dist_min_nm,
        cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
    )
    rigid_total = float(intra0 + inter0)
    state_quantum_nm = min(float(step) for step in step_schedule_nm if float(step) > 0.0)
    lipid1_state_cache: dict[tuple[tuple[int, ...], int], tuple[np.ndarray, float]] = {}
    lipid2_state_cache: dict[tuple[tuple[int, ...], int], tuple[np.ndarray, float]] = {}
    pair_energy_cache: dict[
        tuple[tuple[tuple[int, ...], int], tuple[tuple[int, ...], int]],
        float,
    ] = {}

    def state_key(
        comp_nm: np.ndarray,
        radial_nm: float,
    ) -> tuple[tuple[int, ...], int]:
        comp = np.asarray(comp_nm, dtype=np.float64)
        comp_key = tuple(int(x) for x in np.rint(comp / state_quantum_nm))
        radial_key = int(np.rint(float(radial_nm) / state_quantum_nm))
        return comp_key, radial_key

    def local_state(
        cache: dict[tuple[tuple[int, ...], int], tuple[np.ndarray, float]],
        pos: np.ndarray,
        axis: np.ndarray,
        target_com: np.ndarray,
        tail_weights: np.ndarray,
        radial_unit: np.ndarray,
        comp_nm: np.ndarray,
        radial_nm: float,
        intra_shift_nm: np.ndarray | None = None,
    ) -> tuple[tuple[np.ndarray, float], tuple[tuple[int, ...], int]]:
        key = state_key(comp_nm, radial_nm)
        cached = cache.get(key)
        if cached is not None:
            return cached, key
        local = _recenter_bead_positions_to_target_com(
            _apply_tail_compaction_modes(
                pos,
                axis,
                comp_nm,
                radial_nm,
                tail_weights,
                radial_unit,
            ),
            target_com,
        )
        intra_positions = local
        if intra_shift_nm is not None:
            intra_positions = local - np.asarray(intra_shift_nm, dtype=np.float64)[None, :]
        intra = _compute_dopc_internal_energy(
            intra_positions,
            bead_types,
            bead_charges,
            pair_params,
            lipids_itp_path,
        )
        cached = (local, float(intra))
        cache[key] = cached
        return cached, key

    def total_energy(
        comp1_nm: np.ndarray,
        comp2_nm: np.ndarray,
        radial1_nm: float,
        radial2_nm: float,
    ) -> float:
        (local1, intra1), key1 = local_state(
            lipid1_state_cache,
            pos1,
            dir1,
            target_com1,
            tail_weights1,
            radial_unit1,
            comp1_nm,
            radial1_nm,
        )
        (local2, intra2), key2 = local_state(
            lipid2_state_cache,
            pos2,
            dir2,
            target_com2,
            tail_weights2,
            radial_unit2,
            comp2_nm,
            radial2_nm,
            intra_shift_nm=r_axis,
        )
        pair_key = (key1, key2)
        e_inter = pair_energy_cache.get(pair_key)
        if e_inter is None:
            e_inter, _, _ = _compute_pair_energy_and_gradient(
                local1,
                local2,
                bead_types,
                bead_types,
                bead_charges,
                bead_charges,
                pair_params,
                dist_min_nm=dist_min_nm,
                cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
            )
            e_inter = float(e_inter)
            pair_energy_cache[pair_key] = e_inter
        return float(intra1 + intra2 + e_inter)

    n_chain = len(_PAIR_RELAX_TAIL_CHAINS)
    cur_comp1 = np.zeros(n_chain, dtype=np.float64)
    cur_comp2 = np.zeros(n_chain, dtype=np.float64)
    cur_radial1 = 0.0
    cur_radial2 = 0.0
    cur_total = rigid_total
    for step_nm in step_schedule_nm:
        step = float(step_nm)
        if step <= 0.0:
            continue
        improved = True
        while improved:
            improved = False
            for lipid_id in (0, 1):
                for chain_idx in range(n_chain):
                    best_total = cur_total
                    best_comp1 = cur_comp1.copy()
                    best_comp2 = cur_comp2.copy()
                    best_radial1 = float(cur_radial1)
                    best_radial2 = float(cur_radial2)
                    for sign in (-1.0, 1.0):
                        trial_comp1 = cur_comp1.copy()
                        trial_comp2 = cur_comp2.copy()
                        trial_radial1 = float(cur_radial1)
                        trial_radial2 = float(cur_radial2)
                        if lipid_id == 0:
                            trial_comp1[chain_idx] = float(
                                np.clip(cur_comp1[chain_idx] + sign * step, 0.0, max_compaction_nm)
                            )
                        else:
                            trial_comp2[chain_idx] = float(
                                np.clip(cur_comp2[chain_idx] + sign * step, 0.0, max_compaction_nm)
                            )
                        if np.array_equal(trial_comp1, cur_comp1) and np.array_equal(trial_comp2, cur_comp2):
                            continue
                        e_trial = total_energy(
                            trial_comp1,
                            trial_comp2,
                            trial_radial1,
                            trial_radial2,
                        )
                        if e_trial < best_total - 1.0e-9:
                            best_total = e_trial
                            best_comp1 = trial_comp1
                            best_comp2 = trial_comp2
                            best_radial1 = trial_radial1
                            best_radial2 = trial_radial2
                    if best_total < cur_total - 1.0e-9:
                        cur_comp1 = best_comp1.copy()
                        cur_comp2 = best_comp2.copy()
                        cur_radial1 = float(best_radial1)
                        cur_radial2 = float(best_radial2)
                        cur_total = float(best_total)
                        improved = True
                best_total = cur_total
                best_comp1 = cur_comp1.copy()
                best_comp2 = cur_comp2.copy()
                best_radial1 = float(cur_radial1)
                best_radial2 = float(cur_radial2)
                radial_limit = max_radial1 if lipid_id == 0 else max_radial2
                for sign in (-1.0, 1.0):
                    trial_comp1 = cur_comp1.copy()
                    trial_comp2 = cur_comp2.copy()
                    trial_radial1 = float(cur_radial1)
                    trial_radial2 = float(cur_radial2)
                    if lipid_id == 0:
                        trial_radial1 = float(np.clip(cur_radial1 + sign * step, 0.0, radial_limit))
                    else:
                        trial_radial2 = float(np.clip(cur_radial2 + sign * step, 0.0, radial_limit))
                    if trial_radial1 == cur_radial1 and trial_radial2 == cur_radial2:
                        continue
                    e_trial = total_energy(
                        trial_comp1,
                        trial_comp2,
                        trial_radial1,
                        trial_radial2,
                    )
                    if e_trial < best_total - 1.0e-9:
                        best_total = e_trial
                        best_comp1 = trial_comp1
                        best_comp2 = trial_comp2
                        best_radial1 = trial_radial1
                        best_radial2 = trial_radial2
                if best_total < cur_total - 1.0e-9:
                    cur_comp1 = best_comp1.copy()
                    cur_comp2 = best_comp2.copy()
                    cur_radial1 = float(best_radial1)
                    cur_radial2 = float(best_radial2)
                    cur_total = float(best_total)
                    improved = True

    return rigid_total, float(min(cur_total, rigid_total))


def _pair_conditioned_tail_relaxation_correction_grid(
    ref_bead_positions_nm: np.ndarray,
    ref_bead_positions_nm2: np.ndarray | None,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    average_temperature: float,
    azimuthal_count: int = 1,
    bead_frame_count: int = 1,
    face_cos_min: float = 0.5,
    radial_cutoff_nm: float = 2.5,
    mask_mode: str = "binary",
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
) -> dict:
    refs1 = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    refs2 = refs1 if ref_bead_positions_nm2 is None else _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm2)
    r_values = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta = np.asarray(cos_theta_grid, dtype=np.float64)
    face_weight = _build_cross_leaflet_face_weight(
        r_values,
        cos_theta,
        float(face_cos_min),
        float(radial_cutoff_nm),
        mask_mode=mask_mode,
    )
    phi_values = np.linspace(0.0, 2.0 * np.pi, max(int(azimuthal_count), 1), endpoint=False)
    n12_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    dirs1 = _directions_with_dot_np(-n12_axis, cos_theta, phi_values)
    dirs2 = _directions_with_dot_np(n12_axis, cos_theta, phi_values)
    bead_frame_angles = _bead_frame_angles(max(int(bead_frame_count), 1))
    correction_grid = np.zeros((r_values.size, cos_theta.size, cos_theta.size), dtype=np.float64)
    rigid_grid = np.full_like(correction_grid, np.nan)
    relaxed_grid = np.full_like(correction_grid, np.nan)
    temp = float(average_temperature)
    if temp <= 0.0:
        raise ValueError("pair-conditioned tail relaxation correction requires positive average temperature")

    tasks = []
    for ir, r_nm in enumerate(r_values):
        for ia1 in range(cos_theta.size):
            for ia2 in range(cos_theta.size):
                weight = float(face_weight[ir, ia1, ia2])
                if weight <= 0.0:
                    continue
                tasks.append(
                    {
                        "ir": ir,
                        "ia1": ia1,
                        "ia2": ia2,
                        "r_nm": float(r_nm),
                        "weight": weight,
                        "dirs1": np.asarray(dirs1[ia1], dtype=np.float64),
                        "dirs2": np.asarray(dirs2[ia2], dtype=np.float64),
                        "bead_frame_angles": list(bead_frame_angles),
                        "refs1": np.asarray(refs1, dtype=np.float64),
                        "refs2": np.asarray(refs2, dtype=np.float64),
                        "bead_types": list(bead_types),
                        "bead_charges": list(bead_charges),
                        "pair_params": pair_params,
                        "lipids_itp_path": Path(lipids_itp_path).expanduser().resolve(),
                        "dist_min_nm": float(dist_min_nm),
                        "temperature": temp,
                    }
                )

    for ir, ia1, ia2, rigid_fe, relaxed_fe, weighted_corr in _parallel_map_ordered(
        "Pair-conditioned tail relaxation grid",
        _run_pair_conditioned_tail_relaxation_grid_task,
        tasks,
    ):
        rigid_grid[ir, ia1, ia2] = rigid_fe
        relaxed_grid[ir, ia1, ia2] = relaxed_fe
        correction_grid[ir, ia1, ia2] = weighted_corr

    return {
        "correction_grid_kj_mol": correction_grid,
        "rigid_grid_kj_mol": rigid_grid,
        "relaxed_grid_kj_mol": relaxed_grid,
        "face_mask": np.asarray(face_weight, dtype=np.float32),
        "face_cos_min": float(face_cos_min),
        "radial_cutoff_nm": float(radial_cutoff_nm),
        "mask_mode": str(mask_mode),
    }


def _run_pair_conditioned_tail_relaxation_grid_task(
    task: dict,
) -> tuple[int, int, int, float, float, float]:
    _ensure_cg_bonds_angles(Path(task["lipids_itp_path"]).expanduser().resolve())
    rigid_samples = []
    relaxed_samples = []
    for dir1 in np.asarray(task["dirs1"], dtype=np.float64):
        for dir2 in np.asarray(task["dirs2"], dtype=np.float64):
            for frame_angle1 in task["bead_frame_angles"]:
                for frame_angle2 in task["bead_frame_angles"]:
                    for ref1 in np.asarray(task["refs1"], dtype=np.float64):
                        for ref2 in np.asarray(task["refs2"], dtype=np.float64):
                            rigid, relaxed = _pair_conditioned_tail_relaxation_effective_energy(
                                float(task["r_nm"]),
                                dir1,
                                dir2,
                                float(frame_angle1),
                                float(frame_angle2),
                                ref1,
                                ref2,
                                bead_types=task["bead_types"],
                                bead_charges=task["bead_charges"],
                                pair_params=task["pair_params"],
                                lipids_itp_path=task["lipids_itp_path"],
                                dist_min_nm=float(task["dist_min_nm"]),
                            )
                            rigid_samples.append(rigid)
                            relaxed_samples.append(relaxed)
    rigid_fe = _boltzmann_free_energy_kj_mol(
        np.asarray(rigid_samples, dtype=np.float64),
        float(task["temperature"]),
    )
    relaxed_fe = _boltzmann_free_energy_kj_mol(
        np.asarray(relaxed_samples, dtype=np.float64),
        float(task["temperature"]),
    )
    weighted_corr = min(0.0, float(relaxed_fe - rigid_fe)) * float(task["weight"])
    return (
        int(task["ir"]),
        int(task["ia1"]),
        int(task["ia2"]),
        float(rigid_fe),
        float(relaxed_fe),
        float(weighted_corr),
    )


def _apply_pair_conditioned_tail_relaxation_to_cg_result(
    result_cg: dict,
    ref_bead_positions_nm: np.ndarray,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    face_cos_min: float = 0.5,
    radial_cutoff_nm: float = 2.5,
    mask_mode: str = "binary",
) -> tuple[dict, dict]:
    energy_grid = np.asarray(result_cg["energy_grid_raw"], dtype=np.float64)
    if "reference_energy_eup" not in result_cg:
        raise RuntimeError("pair-conditioned tail relaxation requires reference_energy_eup in result_cg")
    reference_energy_kj_mol = float(result_cg["reference_energy_eup"]) * ENERGY_CONVERSION_KJ_PER_EUP
    correction = _pair_conditioned_tail_relaxation_correction_grid(
        np.asarray(ref_bead_positions_nm, dtype=np.float64),
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        lipids_itp_path=lipids_itp_path,
        ref_bead_positions_nm2=None,
        r_values_nm=np.asarray(result_cg["r_values_nm"], dtype=np.float64),
        cos_theta_grid=np.asarray(result_cg["cos_theta_grid"], dtype=np.float64),
        average_temperature=float(
            result_cg.get(
                "azimuthal_average_temperature_upside",
                DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
            )
        ),
        azimuthal_count=int(result_cg.get("azimuthal_count", 1)),
        bead_frame_count=int(result_cg.get("bead_frame_count", 1)),
        face_cos_min=face_cos_min,
        radial_cutoff_nm=radial_cutoff_nm,
        mask_mode=mask_mode,
        dist_min_nm=float(result_cg.get("sample_dist_min_nm", NUMERICAL_DISTANCE_GUARD_NM)),
    )
    corrected_grid = np.asarray(energy_grid + correction["correction_grid_kj_mol"], dtype=np.float64)
    kbt_kj_mol = float(DEFAULT_PRODUCTION_TEMP_UPSIDE) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt_kj_mol <= 0.0:
        raise RuntimeError("pair-conditioned tail relaxation requires positive control temperature")
    reduced_energy = (corrected_grid - reference_energy_kj_mol) / kbt_kj_mol
    control_grid = np.log1p(np.maximum(reduced_energy, 0.0))
    tensor = _fit_radial_angular_angular_tensor_bspline(
        np.asarray(result_cg["r_values_nm"], dtype=np.float64),
        np.asarray(result_cg["cos_theta_grid"], dtype=np.float64),
        control_grid,
        n_knot_radial=int(result_cg["n_radial"]),
        n_knot_angular=int(result_cg["n_angular"]),
        knot_spacing_ang=float(result_cg["knot_spacing_ang"]),
        energy_conversion=1.0,
        smooth=float(result_cg["fit_smooth"]),
    )
    control_min = float(np.min(control_grid))
    control_max = float(np.max(control_grid))
    updated = dict(result_cg)
    updated["energy_grid_raw"] = corrected_grid
    updated["interaction_param"] = np.clip(tensor, control_min, control_max)
    updated["pair_relaxation_correction_source"] = (
        "source_conditioned_two_lipid_collective_tail_axial_compaction"
    )
    updated["pair_relaxation_energy_basis"] = "total_effective_energy_kj_mol"
    updated["pair_relaxation_face_cos_min"] = float(correction["face_cos_min"])
    updated["pair_relaxation_radial_cutoff_nm"] = float(correction["radial_cutoff_nm"])
    updated["pair_relaxation_mask_mode"] = str(correction.get("mask_mode", mask_mode))
    updated["pair_relaxation_correction_grid_kj_mol"] = correction["correction_grid_kj_mol"]
    updated["pair_relaxation_rigid_grid_kj_mol"] = correction["rigid_grid_kj_mol"]
    updated["pair_relaxation_relaxed_grid_kj_mol"] = correction["relaxed_grid_kj_mol"]
    updated["pair_relaxation_face_mask"] = correction["face_mask"]
    return updated, correction


def apply_pair_conditioned_tail_relaxation_correction_to_dopc_h5(
    table_h5_path: Path,
    dry_ff_path: Path,
    lipids_itp_path: Path,
    face_cos_min: float = 0.5,
    radial_cutoff_nm: float = 2.5,
    mask_mode: str = "binary",
) -> dict:
    table_h5_path = Path(table_h5_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    lipids_itp_path = Path(lipids_itp_path).expanduser().resolve()
    _atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    from martini_itp_reader import parse_dopc_from_itp
    dopc = parse_dopc_from_itp(lipids_itp_path)

    with h5py.File(table_h5_path, "r+") as h5:
        cglt = h5["cg_lipid_table"]
        pair_grp = cglt["cg_lipid_pair"]
        ref_nm = np.asarray(cglt["ref_bead_positions_nm"][:], dtype=np.float64)
        energy_grid = np.asarray(pair_grp["energy_grid_raw_kj_mol"][:], dtype=np.float64)
        n_radial, n_angle1, n_angle2 = energy_grid.shape
        if n_angle1 != n_angle2:
            raise RuntimeError("CGL pair grid must be square in angular dimensions")
        result_cg = {
            "energy_grid_raw": energy_grid,
            "reference_energy_eup": float(pair_grp["reference_energy_eup"][0, 0]),
            "r_values_nm": np.linspace(
                float(pair_grp.attrs["fit_r_min_nm"]),
                float(pair_grp.attrs["fit_r_max_nm"]),
                int(n_radial),
                dtype=np.float64,
            ),
            "cos_theta_grid": np.linspace(-1.0, 1.0, int(n_angle1), dtype=np.float64),
            "n_radial": int(pair_grp.attrs["n_radial"]),
            "n_angular": int(pair_grp.attrs["n_angular"]),
            "knot_spacing_ang": float(pair_grp.attrs["knot_spacing_ang"]),
            "fit_smooth": float(pair_grp.attrs["fit_smooth"]),
            "sample_dist_min_nm": float(pair_grp.attrs.get("sample_dist_min_nm", NUMERICAL_DISTANCE_GUARD_NM)),
            "azimuthal_average_temperature_upside": float(
                pair_grp.attrs.get("azimuthal_average_temperature_upside", DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE)
            ),
        }
        updated, correction = _apply_pair_conditioned_tail_relaxation_to_cg_result(
            result_cg,
            ref_nm,
            bead_types=list(dopc["bead_types"]),
            bead_charges=list(dopc["bead_charges"]),
            pair_params=pair_params,
            lipids_itp_path=lipids_itp_path,
            face_cos_min=face_cos_min,
            radial_cutoff_nm=radial_cutoff_nm,
            mask_mode=mask_mode,
        )
        pair_grp["energy_grid_raw_kj_mol"][...] = np.asarray(updated["energy_grid_raw"], dtype=np.float32)
        pair_grp["interaction_param"][...] = np.asarray(updated["interaction_param"], dtype=np.float32).reshape(1, 1, -1)
        pair_grp["reference_energy_eup"][...] = np.asarray([[updated["reference_energy_eup"]]], dtype=np.float32)
        pair_grp.attrs["pair_relaxation_correction_source"] = updated["pair_relaxation_correction_source"]
        pair_grp.attrs["pair_relaxation_energy_basis"] = updated.get(
            "pair_relaxation_energy_basis",
            "intermolecular_energy_kj_mol",
        )
        pair_grp.attrs["pair_relaxation_face_cos_min"] = np.float32(correction["face_cos_min"])
        pair_grp.attrs["pair_relaxation_radial_cutoff_nm"] = np.float32(correction["radial_cutoff_nm"])
        pair_grp.attrs["pair_relaxation_mask_mode"] = str(correction.get("mask_mode", mask_mode))
        for dset_name, values in (
            ("pair_relaxation_correction_grid_kj_mol", updated["pair_relaxation_correction_grid_kj_mol"]),
            ("pair_relaxation_rigid_grid_kj_mol", updated["pair_relaxation_rigid_grid_kj_mol"]),
            ("pair_relaxation_relaxed_grid_kj_mol", updated["pair_relaxation_relaxed_grid_kj_mol"]),
            ("pair_relaxation_face_mask", updated["pair_relaxation_face_mask"]),
        ):
            if dset_name in pair_grp:
                del pair_grp[dset_name]
            pair_grp.create_dataset(dset_name, data=np.asarray(values))
    print(f"Applied pair-conditioned tail-relaxation correction to {table_h5_path}")
    return correction


def _build_cross_leaflet_face_mask(
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    face_cos_min: float,
    radial_cutoff_nm: float,
) -> np.ndarray:
    return _build_cross_leaflet_face_weight(
        r_values_nm,
        cos_theta_grid,
        face_cos_min,
        radial_cutoff_nm,
        mask_mode="binary",
    ) > 0.5


def _build_cross_leaflet_face_weight(
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    face_cos_min: float,
    radial_cutoff_nm: float,
    mask_mode: str = "binary",
) -> np.ndarray:
    r_values = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta = np.asarray(cos_theta_grid, dtype=np.float64)
    face_cos = float(face_cos_min)
    radial_cutoff = float(radial_cutoff_nm)
    if mask_mode == "binary":
        return (
            (r_values[:, None, None] <= radial_cutoff)
            & (cos_theta[None, :, None] <= -face_cos)
            & (cos_theta[None, None, :] <= -face_cos)
        ).astype(np.float64)
    if mask_mode != "smooth":
        raise ValueError(f"Unknown cross-leaflet face mask mode {mask_mode!r}")
    denom = max(1.0 - face_cos, 1.0e-6)
    ang1 = np.clip(((-cos_theta) - face_cos) / denom, 0.0, 1.0)
    ang2 = np.clip(((-cos_theta) - face_cos) / denom, 0.0, 1.0)
    radial = np.clip(1.0 - np.maximum(r_values, 0.0) / max(radial_cutoff, 1.0e-6), 0.0, 1.0)
    return radial[:, None, None] * ang1[None, :, None] * ang2[None, None, :]


def _fit_compaction_result_for_cg_result(
    result_cg: dict,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    face_cos_min: float | None = None,
    radial_cutoff_nm: float = 2.5,
    pair_relax_state_correction: bool = False,
    mask_mode: str = "binary",
    correction_center_mode: str = "base",
    pair_state_model: str = "bilinear",
) -> dict:
    from martini_itp_reader import parse_dopc_from_itp

    compaction_pool_conformers = _positive_int_env("UPSIDE_CGL_COMPACTION_POOL_CONFORMERS", 32)
    compaction_burnin_steps = _positive_int_env("UPSIDE_CGL_COMPACTION_BURNIN_STEPS", 20000)
    compaction_steps_per_conf = _positive_int_env("UPSIDE_CGL_COMPACTION_STEPS_PER_CONFORMER", 500)
    compaction_representatives = _positive_int_env("UPSIDE_CGL_COMPACTION_STATE_REPRESENTATIVES", 2)
    compaction_self_bins = _positive_int_env("UPSIDE_CGL_COMPACTION_SELF_BINS", 12)
    compaction_fit_smooth = _float_env(
        "UPSIDE_CGL_COMPACTION_FIT_SMOOTH",
        float(result_cg.get("fit_smooth", 0.1)),
    )
    dopc = parse_dopc_from_itp(lipids_itp_path)
    compaction_refs_nm = sample_isolated_dopc_bonded_conformers(
        dopc,
        lipids_itp_path=lipids_itp_path,
        pair_params=pair_params,
        conformer_count=compaction_pool_conformers,
        temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
        seed=1777,
        mc_burnin_steps=compaction_burnin_steps,
        mc_steps_per_conformer=compaction_steps_per_conf,
    )
    compaction_values_ang = _dopc_tail_extension_series_ang(compaction_refs_nm)
    compaction_states = _select_compaction_state_representatives(
        compaction_refs_nm,
        compaction_values_ang,
        representative_count=compaction_representatives,
    )
    compaction_self = _fit_compaction_self_pmf(
        compaction_values_ang,
        temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
        n_bin=compaction_self_bins,
        smooth=0.01,
    )
    compaction_result = _fit_pair_compaction_result_for_states(
        result_cg=result_cg,
        compaction_states=compaction_states,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        lipids_itp_path=lipids_itp_path,
        face_cos_min=face_cos_min,
        radial_cutoff_nm=radial_cutoff_nm,
        pair_relax_state_correction=pair_relax_state_correction,
        mask_mode=mask_mode,
        correction_center_mode=correction_center_mode,
        pair_state_model=pair_state_model,
        fit_smooth=compaction_fit_smooth,
    )
    compaction_result["self"] = compaction_self
    return compaction_result


def _fit_pair_compaction_result_for_states(
    result_cg: dict,
    compaction_states: dict,
    bead_types: list[str],
    bead_charges: list[float],
    pair_params: dict,
    lipids_itp_path: Path,
    face_cos_min: float | None = None,
    radial_cutoff_nm: float = 2.5,
    pair_relax_state_correction: bool = False,
    mask_mode: str = "binary",
    correction_center_mode: str = "base",
    pair_state_model: str = "bilinear",
    fit_smooth: float | None = None,
) -> dict:
    compaction_fit_smooth = _float_env(
        "UPSIDE_CGL_COMPACTION_FIT_SMOOTH",
        float(result_cg.get("fit_smooth", 0.1) if fit_smooth is None else fit_smooth),
    )
    compressed_center_ang = compaction_states.get("compressed_center_ang")
    has_compressed_state = (
        "compressed_refs_nm" in compaction_states
        and compressed_center_ang is not None
        and np.isfinite(float(compressed_center_ang))
        and float(compressed_center_ang) > float(compaction_states["compact_center_ang"]) + 1.0e-6
    )

    def sample_state_grid(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
        return _sample_cg_pair_energy_grid(
            refs_i_nm,
            refs_j_nm,
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            r_values_nm=result_cg["r_values_nm"],
            cos_theta_grid=result_cg["cos_theta_grid"],
            azimuthal_count=int(result_cg["azimuthal_count"]),
            bead_frame_count=int(result_cg["bead_frame_count"]),
            dist_min_nm=1.0e-6,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )

    def sample_symmetric_state_grid(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
        forward = sample_state_grid(refs_i_nm, refs_j_nm)
        reverse = sample_state_grid(refs_j_nm, refs_i_nm)
        return 0.5 * (forward + np.swapaxes(reverse, 1, 2))

    state_grid_ee = sample_state_grid(
        compaction_states["extended_refs_nm"],
        compaction_states["extended_refs_nm"],
    )
    state_grid_cc = sample_state_grid(
        compaction_states["compact_refs_nm"],
        compaction_states["compact_refs_nm"],
    )
    state_grid_ec = sample_symmetric_state_grid(
        compaction_states["extended_refs_nm"],
        compaction_states["compact_refs_nm"],
    )
    state_grid_xx = None
    state_grid_ex = None
    state_grid_cx = None
    if has_compressed_state:
        state_grid_xx = sample_state_grid(
            compaction_states["compressed_refs_nm"],
            compaction_states["compressed_refs_nm"],
        )
        state_grid_ex = sample_symmetric_state_grid(
            compaction_states["extended_refs_nm"],
            compaction_states["compressed_refs_nm"],
        )
        state_grid_cx = sample_symmetric_state_grid(
            compaction_states["compact_refs_nm"],
            compaction_states["compressed_refs_nm"],
        )
    if pair_relax_state_correction:
        pr_face_cos_min = _float_env("UPSIDE_CGL_PAIR_RELAX_FACE_COS_MIN", 0.5)
        pr_radial_cutoff_nm = _float_env("UPSIDE_CGL_PAIR_RELAX_RADIAL_CUTOFF_NM", 2.5)
        def sample_relax_correction(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray | None) -> dict:
            return _pair_conditioned_tail_relaxation_correction_grid(
                refs_i_nm,
                refs_j_nm,
                bead_types=bead_types,
                bead_charges=bead_charges,
                pair_params=pair_params,
                lipids_itp_path=lipids_itp_path,
                r_values_nm=result_cg["r_values_nm"],
                cos_theta_grid=result_cg["cos_theta_grid"],
                average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                azimuthal_count=int(result_cg.get("azimuthal_count", 1)),
                bead_frame_count=int(result_cg.get("bead_frame_count", 1)),
                face_cos_min=pr_face_cos_min,
                radial_cutoff_nm=pr_radial_cutoff_nm,
                dist_min_nm=1.0e-6,
            )

        def symmetric_relax_correction(refs_i_nm: np.ndarray, refs_j_nm: np.ndarray) -> np.ndarray:
            forward = sample_relax_correction(refs_i_nm, refs_j_nm)
            reverse = sample_relax_correction(refs_j_nm, refs_i_nm)
            return 0.5 * (
                np.asarray(forward["correction_grid_kj_mol"], dtype=np.float64)
                + np.swapaxes(np.asarray(reverse["correction_grid_kj_mol"], dtype=np.float64), 1, 2)
            )

        state_grid_ee = state_grid_ee + np.asarray(
            sample_relax_correction(compaction_states["extended_refs_nm"], None)["correction_grid_kj_mol"],
            dtype=np.float64,
        )
        state_grid_cc = state_grid_cc + np.asarray(
            sample_relax_correction(compaction_states["compact_refs_nm"], None)["correction_grid_kj_mol"],
            dtype=np.float64,
        )
        state_grid_ec = state_grid_ec + symmetric_relax_correction(
            compaction_states["extended_refs_nm"],
            compaction_states["compact_refs_nm"],
        )
        if has_compressed_state:
            state_grid_xx = state_grid_xx + np.asarray(
                sample_relax_correction(compaction_states["compressed_refs_nm"], None)["correction_grid_kj_mol"],
                dtype=np.float64,
            )
            state_grid_ex = state_grid_ex + symmetric_relax_correction(
                compaction_states["extended_refs_nm"],
                compaction_states["compressed_refs_nm"],
            )
            state_grid_cx = state_grid_cx + symmetric_relax_correction(
                compaction_states["compact_refs_nm"],
                compaction_states["compressed_refs_nm"],
            )
    compact_probability = float(compaction_states["compact_probability"])
    state_model = str(pair_state_model).strip().lower()
    if has_compressed_state and state_model != "bilinear":
        raise ValueError(
            "Compressed-state pair compaction requires the bilinear pair-state model"
        )
    if has_compressed_state:
        extended_probability = float(
            compaction_states.get("extended_probability", max(0.0, 1.0 - compact_probability))
        )
        compact_middle_probability = float(
            compaction_states.get("compact_middle_probability", compact_probability)
        )
        compressed_probability = float(compaction_states.get("compressed_probability", 0.0))
        total_probability = (
            extended_probability + compact_middle_probability + compressed_probability
        )
        if total_probability <= 1.0e-12:
            extended_probability = max(0.0, 1.0 - compact_probability)
            compact_middle_probability = compact_probability
            compressed_probability = 0.0
            total_probability = extended_probability + compact_middle_probability
        extended_probability /= total_probability
        compact_middle_probability /= total_probability
        compressed_probability /= total_probability
        state_grid_avg = (
            (extended_probability * extended_probability) * state_grid_ee
            + (2.0 * extended_probability * compact_middle_probability) * state_grid_ec
            + (compact_middle_probability * compact_middle_probability) * state_grid_cc
            + (2.0 * extended_probability * compressed_probability) * state_grid_ex
            + (2.0 * compact_middle_probability * compressed_probability) * state_grid_cx
            + (compressed_probability * compressed_probability) * state_grid_xx
        )
    else:
        extended_probability = 1.0 - compact_probability
        state_grid_avg = (
            (extended_probability * extended_probability) * state_grid_ee
            + (2.0 * compact_probability * extended_probability) * state_grid_ec
            + (compact_probability * compact_probability) * state_grid_cc
        )
    base_grid_kj_mol = np.asarray(result_cg["energy_grid_raw"], dtype=np.float64)
    reference_energy_kj_mol = float(result_cg["reference_energy_eup"]) * ENERGY_CONVERSION_KJ_PER_EUP
    control_kbt_kj_mol = DEFAULT_PRODUCTION_TEMP_UPSIDE * ENERGY_CONVERSION_KJ_PER_EUP
    if control_kbt_kj_mol <= 0.0:
        raise ValueError("CGL compaction correction requires positive control temperature")

    def to_log1p_control(grid_kj_mol: np.ndarray) -> np.ndarray:
        reduced = (np.asarray(grid_kj_mol, dtype=np.float64) - reference_energy_kj_mol) / control_kbt_kj_mol
        return np.log1p(np.maximum(reduced, 0.0))

    base_control_grid = to_log1p_control(base_grid_kj_mol)
    state_avg_control_grid = to_log1p_control(state_grid_avg)
    center_mode = str(correction_center_mode).strip().lower()
    if center_mode == "base":
        center_control_grid = base_control_grid
    elif center_mode == "state_average":
        center_control_grid = state_avg_control_grid
    else:
        raise ValueError(f"Unsupported compaction correction center mode: {correction_center_mode!r}")
    correction_ee_control = to_log1p_control(state_grid_ee) - center_control_grid
    correction_ec_control = to_log1p_control(state_grid_ec) - center_control_grid
    correction_cc_control = to_log1p_control(state_grid_cc) - center_control_grid
    correction_ex_control = None
    correction_cx_control = None
    correction_xx_control = None
    if has_compressed_state:
        correction_ex_control = to_log1p_control(state_grid_ex) - center_control_grid
        correction_cx_control = to_log1p_control(state_grid_cx) - center_control_grid
        correction_xx_control = to_log1p_control(state_grid_xx) - center_control_grid
    if state_model == "additive":
        correction_cc_control = 2.0 * correction_ec_control - correction_ee_control
    elif state_model != "bilinear":
        raise ValueError(f"Unsupported compaction pair-state model: {pair_state_model!r}")
    compaction_face_weight = None
    if face_cos_min is not None:
        compaction_face_weight = _build_cross_leaflet_face_weight(
            result_cg["r_values_nm"],
            result_cg["cos_theta_grid"],
            float(face_cos_min),
            float(radial_cutoff_nm),
            mask_mode=mask_mode,
        )
        if mask_mode == "binary":
            compaction_face_mask = compaction_face_weight > 0.5
            correction_ee_control = np.where(compaction_face_mask, correction_ee_control, 0.0)
            correction_ec_control = np.where(compaction_face_mask, correction_ec_control, 0.0)
            correction_cc_control = np.where(compaction_face_mask, correction_cc_control, 0.0)
            if correction_ex_control is not None:
                correction_ex_control = np.where(compaction_face_mask, correction_ex_control, 0.0)
                correction_cx_control = np.where(compaction_face_mask, correction_cx_control, 0.0)
                correction_xx_control = np.where(compaction_face_mask, correction_xx_control, 0.0)
        else:
            correction_ee_control = correction_ee_control * compaction_face_weight
            correction_ec_control = correction_ec_control * compaction_face_weight
            correction_cc_control = correction_cc_control * compaction_face_weight
            if correction_ex_control is not None:
                correction_ex_control = correction_ex_control * compaction_face_weight
                correction_cx_control = correction_cx_control * compaction_face_weight
                correction_xx_control = correction_xx_control * compaction_face_weight

    def fit_correction_tensor(correction_grid_control: np.ndarray) -> np.ndarray:
        tensor = _fit_radial_angular_angular_tensor_bspline(
            result_cg["r_values_nm"],
            result_cg["cos_theta_grid"],
            correction_grid_control,
            n_knot_radial=int(result_cg["n_radial"]),
            n_knot_angular=int(result_cg["n_angular"]),
            knot_spacing_ang=float(result_cg["knot_spacing_ang"]),
            energy_conversion=1.0,
            smooth=compaction_fit_smooth,
        )
        corr_min = float(np.min(correction_grid_control))
        corr_max = float(np.max(correction_grid_control))
        return np.clip(tensor, corr_min, corr_max)

    compaction_result = {
        "compact_center_ang": float(compaction_states["compact_center_ang"]),
        "extended_center_ang": float(compaction_states["extended_center_ang"]),
        "compact_probability": compact_probability,
        "pool_size": int(compaction_states["pool_size"]),
        "compact_pool_size": int(compaction_states["compact_pool_size"]),
        "extended_pool_size": int(compaction_states["extended_pool_size"]),
        "representative_compact_count": int(compaction_states["compact_refs_nm"].shape[0]),
        "representative_extended_count": int(compaction_states["extended_refs_nm"].shape[0]),
        "compaction_min_ang": float(compaction_states["compaction_min_ang"]),
        "compaction_max_ang": float(compaction_states["compaction_max_ang"]),
        "compaction_p05_ang": float(compaction_states["compaction_p05_ang"]),
        "compaction_p95_ang": float(compaction_states["compaction_p95_ang"]),
        "delta_extended_extended": fit_correction_tensor(correction_ee_control).reshape(-1),
        "delta_extended_compact": fit_correction_tensor(correction_ec_control).reshape(-1),
        "delta_compact_compact": fit_correction_tensor(correction_cc_control).reshape(-1),
        "grid_extended_extended_kj_mol": state_grid_ee,
        "grid_extended_compact_kj_mol": state_grid_ec,
        "grid_compact_compact_kj_mol": state_grid_cc,
        "grid_average_kj_mol": state_grid_avg,
        "correction_center_mode": center_mode,
        "pair_state_model": state_model,
        "fit_smooth": float(compaction_fit_smooth),
    }
    if has_compressed_state:
        compaction_result["delta_extended_compressed"] = fit_correction_tensor(
            correction_ex_control
        ).reshape(-1)
        compaction_result["delta_compact_compressed"] = fit_correction_tensor(
            correction_cx_control
        ).reshape(-1)
        compaction_result["delta_compressed_compressed"] = fit_correction_tensor(
            correction_xx_control
        ).reshape(-1)
        compaction_result["grid_extended_compressed_kj_mol"] = state_grid_ex
        compaction_result["grid_compact_compressed_kj_mol"] = state_grid_cx
        compaction_result["grid_compressed_compressed_kj_mol"] = state_grid_xx
        compaction_result["representative_compressed_count"] = int(
            compaction_states["compressed_refs_nm"].shape[0]
        )
        compaction_result["compressed_pool_size"] = int(
            compaction_states.get("compressed_pool_size", compaction_states["compressed_refs_nm"].shape[0])
        )
        compaction_result["compressed_center_ang"] = float(compressed_center_ang)
    if compaction_face_weight is not None:
        compaction_result["face_mask"] = np.asarray(compaction_face_weight, dtype=np.float32)
        compaction_result["face_cos_min"] = float(face_cos_min)
        compaction_result["radial_cutoff_nm"] = float(radial_cutoff_nm)
        compaction_result["mask_mode"] = str(mask_mode)
    return compaction_result


def _log1p_reduced_control_from_energy_grid(
    energy_grid_kj_mol: np.ndarray,
    reference_energy_eup: float,
    temperature_upside: float,
) -> np.ndarray:
    energy_grid_kj_mol = np.asarray(energy_grid_kj_mol, dtype=np.float64)
    reference_energy_kj_mol = float(reference_energy_eup) * ENERGY_CONVERSION_KJ_PER_EUP
    kbt_kj_mol = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt_kj_mol <= 0.0:
        raise ValueError("log1p reduced control conversion requires positive temperature")
    reduced_energy = (energy_grid_kj_mol - reference_energy_kj_mol) / kbt_kj_mol
    return np.log1p(np.maximum(reduced_energy, 0.0))


def _fit_single_cgl_state_delta_full_tensor(
    base_energy_grid_kj_mol: np.ndarray,
    extended_energy_grid_kj_mol: np.ndarray,
    compact_energy_grid_kj_mol: np.ndarray,
    compressed_energy_grid_kj_mol: np.ndarray | None,
    reference_energy_eup: float,
    temperature_upside: float,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    n_knot_radial: int,
    n_knot_angular: int,
    knot_spacing_ang: float,
    smooth: float,
    base_interaction_param: np.ndarray | None = None,
) -> dict:
    if base_interaction_param is not None:
        base_control = _evaluate_radial_angular_angular_tensor_controls(
            np.asarray(base_interaction_param, dtype=np.float64).reshape(
                int(n_knot_radial),
                int(n_knot_angular),
                int(n_knot_angular),
            ),
            r_values_nm,
            cos_theta_grid,
            knot_spacing_ang,
        )
    else:
        base_control = _log1p_reduced_control_from_energy_grid(
            base_energy_grid_kj_mol,
            reference_energy_eup,
            temperature_upside,
        )
    extended_control = _log1p_reduced_control_from_energy_grid(
        extended_energy_grid_kj_mol,
        reference_energy_eup,
        temperature_upside,
    )
    compact_control = _log1p_reduced_control_from_energy_grid(
        compact_energy_grid_kj_mol,
        reference_energy_eup,
        temperature_upside,
    )
    compressed_control = None
    if compressed_energy_grid_kj_mol is not None:
        compressed_control = _log1p_reduced_control_from_energy_grid(
            compressed_energy_grid_kj_mol,
            reference_energy_eup,
            temperature_upside,
        )
    delta_extended_grid = extended_control - base_control
    delta_compact_grid = compact_control - base_control

    def fit_delta(delta_grid: np.ndarray) -> np.ndarray:
        tensor = _fit_radial_angular_angular_tensor_bspline(
            r_values_nm,
            cos_theta_grid,
            delta_grid,
            n_knot_radial=n_knot_radial,
            n_knot_angular=n_knot_angular,
            knot_spacing_ang=knot_spacing_ang,
            energy_conversion=1.0,
            smooth=float(smooth),
        )
        delta_min = float(np.min(delta_grid))
        delta_max = float(np.max(delta_grid))
        return np.clip(tensor, delta_min, delta_max)

    payload = {
        "delta_extended": fit_delta(delta_extended_grid).reshape(-1),
        "delta_compact": fit_delta(delta_compact_grid).reshape(-1),
        "grid_extended_kj_mol": np.asarray(extended_energy_grid_kj_mol, dtype=np.float32),
        "grid_compact_kj_mol": np.asarray(compact_energy_grid_kj_mol, dtype=np.float32),
    }
    if compressed_control is not None and compressed_energy_grid_kj_mol is not None:
        payload["delta_compressed"] = fit_delta(compressed_control - base_control).reshape(-1)
        payload["grid_compressed_kj_mol"] = np.asarray(
            compressed_energy_grid_kj_mol,
            dtype=np.float32,
        )
    return payload


def _build_single_cgl_state_delta_radial_angular(
    base_energy_grid_kj_mol: np.ndarray,
    extended_energy_grid_kj_mol: np.ndarray,
    compact_energy_grid_kj_mol: np.ndarray,
    compressed_energy_grid_kj_mol: np.ndarray | None,
    reference_energy_eup: float,
    temperature_upside: float,
    r_values_nm: np.ndarray | None = None,
    cos_theta_grid: np.ndarray | None = None,
    n_knot_radial: int | None = None,
    n_knot_angular: int | None = None,
    knot_spacing_ang: float | None = None,
    smooth: float = 0.01,
    base_interaction_param: np.ndarray | None = None,
) -> dict:
    if (
        base_interaction_param is not None
        and r_values_nm is not None
        and cos_theta_grid is not None
        and n_knot_radial is not None
        and n_knot_angular is not None
        and knot_spacing_ang is not None
    ):
        base_control = _evaluate_radial_angular_tensor_controls(
            np.asarray(base_interaction_param, dtype=np.float64).reshape(
                int(n_knot_radial),
                int(n_knot_angular),
            ),
            np.asarray(r_values_nm, dtype=np.float64),
            np.asarray(cos_theta_grid, dtype=np.float64),
            float(knot_spacing_ang),
        )
    else:
        base_control = _log1p_reduced_control_from_energy_grid(
            base_energy_grid_kj_mol,
            reference_energy_eup,
            temperature_upside,
        )
    extended_control = _log1p_reduced_control_from_energy_grid(
        extended_energy_grid_kj_mol,
        reference_energy_eup,
        temperature_upside,
    )
    compact_control = _log1p_reduced_control_from_energy_grid(
        compact_energy_grid_kj_mol,
        reference_energy_eup,
        temperature_upside,
    )
    compressed_control = None
    if compressed_energy_grid_kj_mol is not None:
        compressed_control = _log1p_reduced_control_from_energy_grid(
            compressed_energy_grid_kj_mol,
            reference_energy_eup,
            temperature_upside,
        )
    delta_extended_grid = extended_control - base_control
    delta_compact_grid = compact_control - base_control
    if (
        r_values_nm is not None
        and cos_theta_grid is not None
        and n_knot_radial is not None
        and n_knot_angular is not None
        and knot_spacing_ang is not None
    ):
        delta_extended = _fit_radial_angular_tensor_bspline(
            np.asarray(r_values_nm, dtype=np.float64),
            np.asarray(cos_theta_grid, dtype=np.float64),
            np.asarray(delta_extended_grid, dtype=np.float64),
            int(n_knot_radial),
            int(n_knot_angular),
            float(knot_spacing_ang),
            energy_conversion=1.0,
            smooth=float(smooth),
        ).reshape(-1)
        delta_compact = _fit_radial_angular_tensor_bspline(
            np.asarray(r_values_nm, dtype=np.float64),
            np.asarray(cos_theta_grid, dtype=np.float64),
            np.asarray(delta_compact_grid, dtype=np.float64),
            int(n_knot_radial),
            int(n_knot_angular),
            float(knot_spacing_ang),
            energy_conversion=1.0,
            smooth=float(smooth),
        ).reshape(-1)
        delta_compressed = None
        if compressed_control is not None:
            delta_compressed = _fit_radial_angular_tensor_bspline(
                np.asarray(r_values_nm, dtype=np.float64),
                np.asarray(cos_theta_grid, dtype=np.float64),
                np.asarray(compressed_control - base_control, dtype=np.float64),
                int(n_knot_radial),
                int(n_knot_angular),
                float(knot_spacing_ang),
                energy_conversion=1.0,
                smooth=float(smooth),
            ).reshape(-1)
    else:
        delta_extended = delta_extended_grid.reshape(-1)
        delta_compact = delta_compact_grid.reshape(-1)
        delta_compressed = None if compressed_control is None else (compressed_control - base_control).reshape(-1)
    payload = {
        "delta_extended": np.asarray(delta_extended, dtype=np.float32),
        "delta_compact": np.asarray(delta_compact, dtype=np.float32),
        "grid_extended_kj_mol": np.asarray(extended_energy_grid_kj_mol, dtype=np.float32),
        "grid_compact_kj_mol": np.asarray(compact_energy_grid_kj_mol, dtype=np.float32),
    }
    if delta_compressed is not None and compressed_energy_grid_kj_mol is not None:
        payload["delta_compressed"] = np.asarray(delta_compressed, dtype=np.float32)
        payload["grid_compressed_kj_mol"] = np.asarray(
            compressed_energy_grid_kj_mol,
            dtype=np.float32,
        )
    return payload


def apply_compaction_correction_to_dopc_h5(
    table_h5_path: Path,
    dry_ff_path: Path,
    lipids_itp_path: Path,
    face_cos_min: float | None = 0.5,
    radial_cutoff_nm: float = 2.5,
    pair_relax_state_correction: bool = False,
    mask_mode: str = "binary",
    correction_center_mode: str = "base",
    pair_state_model: str = "bilinear",
) -> dict:
    table_h5_path = Path(table_h5_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    lipids_itp_path = Path(lipids_itp_path).expanduser().resolve()
    _atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    from martini_itp_reader import parse_dopc_from_itp

    dopc = parse_dopc_from_itp(lipids_itp_path)
    with h5py.File(table_h5_path, "r+") as h5:
        cglt = h5["cg_lipid_table"]
        pair_grp = cglt["cg_lipid_pair"]
        energy_grid = np.asarray(pair_grp["energy_grid_raw_kj_mol"][:], dtype=np.float64)
        n_radial, n_angle1, n_angle2 = energy_grid.shape
        if n_angle1 != n_angle2:
            raise RuntimeError("CGL pair grid must be square in angular dimensions")
        result_cg = {
            "energy_grid_raw": energy_grid,
            "reference_energy_eup": float(pair_grp["reference_energy_eup"][0, 0]),
            "r_values_nm": np.linspace(
                float(pair_grp.attrs["fit_r_min_nm"]),
                float(pair_grp.attrs["fit_r_max_nm"]),
                int(n_radial),
                dtype=np.float64,
            ),
            "cos_theta_grid": np.linspace(-1.0, 1.0, int(n_angle1), dtype=np.float64),
            "n_radial": int(pair_grp.attrs["n_radial"]),
            "n_angular": int(pair_grp.attrs["n_angular"]),
            "knot_spacing_ang": float(pair_grp.attrs["knot_spacing_ang"]),
            "fit_smooth": float(pair_grp.attrs["fit_smooth"]),
            "azimuthal_count": int(pair_grp.attrs["azimuthal_count"]),
            "bead_frame_count": int(pair_grp.attrs["cgl_bead_frame_count"]),
        }
        compaction_result = _fit_compaction_result_for_cg_result(
            result_cg,
            bead_types=list(dopc["bead_types"]),
            bead_charges=list(dopc["bead_charges"]),
            pair_params=pair_params,
            lipids_itp_path=lipids_itp_path,
            face_cos_min=face_cos_min,
            radial_cutoff_nm=radial_cutoff_nm,
            pair_relax_state_correction=pair_relax_state_correction,
            mask_mode=mask_mode,
            correction_center_mode=correction_center_mode,
            pair_state_model=pair_state_model,
        )
        pair_grp.attrs["compaction_pair_correction"] = (
            "isolated_source_two_state_log1p_control_delta_relative_to_base_pair_table"
        )
        pair_grp.attrs["compact_state_center_ang"] = np.float32(compaction_result["compact_center_ang"])
        pair_grp.attrs["extended_state_center_ang"] = np.float32(compaction_result["extended_center_ang"])
        if "cg_lipid_compaction" in cglt:
            del cglt["cg_lipid_compaction"]
        comp_grp = cglt.create_group("cg_lipid_compaction")
        compaction_tau = float(os.environ.get("CG_LIPID_COMPACTION_THERMOSTAT_TIMESCALE", "5.0"))
        _write_common_table_contract_attrs(
            comp_grp,
            table_family="CGL-Compaction",
            source_object="isolated_dopc_axial_extension_coordinate",
            target_object="dynamic_cgl_hidden_axial_extension_coordinate",
            projection_ensemble="isolated_dopc_mc_ensemble",
            runtime_representation=(
                "clamped_bspline_self_pmf_plus_two_state_log1p_control_delta_relative_to_base_pair_table"
            ),
            correction_layer="isolated_source_compaction_state",
        )
        comp_grp.attrs["schema"] = "cg_lipid_compaction_v1"
        comp_grp.attrs["coordinate_name"] = "head_to_tail_midpoint_length_ang"
        comp_grp.attrs["coordinate_reference"] = "head_to_tail_midpoint_distance"
        comp_grp.attrs["boltzmann_temperature_upside"] = np.float32(DEFAULT_PRODUCTION_TEMP_UPSIDE)
        comp_grp.attrs["thermostat_timescale"] = np.float32(compaction_tau)
        comp_grp.attrs["mass_up"] = np.float32(
            compaction_result["self"]["effective_stiffness_eup_a2"] * compaction_tau * compaction_tau
        )
        comp_grp.attrs["effective_stiffness_eup_a2"] = np.float32(
            compaction_result["self"]["effective_stiffness_eup_a2"]
        )
        comp_grp.attrs["self_coord_min_ang"] = np.float32(compaction_result["self"]["coord_min_ang"])
        comp_grp.attrs["self_coord_max_ang"] = np.float32(compaction_result["self"]["coord_max_ang"])
        comp_grp.attrs["self_coord_spacing_ang"] = np.float32(
            compaction_result["self"]["coord_spacing_ang"]
        )
        comp_grp.attrs["self_n_knot"] = np.int32(compaction_result["self"]["n_knot"])
        comp_grp.attrs["compact_state_center_ang"] = np.float32(compaction_result["compact_center_ang"])
        comp_grp.attrs["extended_state_center_ang"] = np.float32(compaction_result["extended_center_ang"])
        comp_grp.attrs["compact_state_probability"] = np.float32(compaction_result["compact_probability"])
        comp_grp.attrs["correction_center_mode"] = compaction_result["correction_center_mode"]
        comp_grp.attrs["pair_state_model"] = compaction_result["pair_state_model"]
        if "face_cos_min" in compaction_result:
            comp_grp.attrs["face_cos_min"] = np.float32(compaction_result["face_cos_min"])
            comp_grp.attrs["radial_cutoff_nm"] = np.float32(compaction_result["radial_cutoff_nm"])
            comp_grp.attrs["mask_mode"] = compaction_result["mask_mode"]
        comp_grp.attrs["compaction_pool_size"] = np.int32(compaction_result["pool_size"])
        comp_grp.attrs["compact_pool_size"] = np.int32(compaction_result["compact_pool_size"])
        comp_grp.attrs["extended_pool_size"] = np.int32(compaction_result["extended_pool_size"])
        comp_grp.attrs["representative_compact_count"] = np.int32(compaction_result["representative_compact_count"])
        comp_grp.attrs["representative_extended_count"] = np.int32(compaction_result["representative_extended_count"])
        comp_grp.attrs["compaction_min_ang"] = np.float32(compaction_result["compaction_min_ang"])
        comp_grp.attrs["compaction_max_ang"] = np.float32(compaction_result["compaction_max_ang"])
        comp_grp.attrs["compaction_p05_ang"] = np.float32(compaction_result["compaction_p05_ang"])
        comp_grp.attrs["compaction_p95_ang"] = np.float32(compaction_result["compaction_p95_ang"])
        for name, values in (
            ("self_coeff", compaction_result["self"]["self_coeff_eup"]),
            ("pmf_centers_ang", compaction_result["self"]["pmf_centers_ang"]),
            ("pmf_values_kj_mol", compaction_result["self"]["pmf_values_kj_mol"]),
            ("delta_extended_extended", compaction_result["delta_extended_extended"]),
            ("delta_extended_compact", compaction_result["delta_extended_compact"]),
            ("delta_compact_compact", compaction_result["delta_compact_compact"]),
            ("grid_extended_extended_kj_mol", compaction_result["grid_extended_extended_kj_mol"]),
            ("grid_extended_compact_kj_mol", compaction_result["grid_extended_compact_kj_mol"]),
            ("grid_compact_compact_kj_mol", compaction_result["grid_compact_compact_kj_mol"]),
            ("grid_average_kj_mol", compaction_result["grid_average_kj_mol"]),
        ):
            comp_grp.create_dataset(name, data=np.asarray(values, dtype=np.float32))
        for name in (
            "delta_extended_compressed",
            "delta_compact_compressed",
            "delta_compressed_compressed",
            "grid_extended_compressed_kj_mol",
            "grid_compact_compressed_kj_mol",
            "grid_compressed_compressed_kj_mol",
        ):
            if name in compaction_result:
                comp_grp.create_dataset(name, data=np.asarray(compaction_result[name], dtype=np.float32))
        if "face_mask" in compaction_result:
            comp_grp.create_dataset("face_mask", data=np.asarray(compaction_result["face_mask"], dtype=np.float32))
    print(f"Applied compaction correction to {table_h5_path}")
    return compaction_result


def _lipid_pair_param_matrices(
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n1 = len(bead_types1)
    n2 = len(bead_types2)
    sigma = np.zeros((n1, n2), dtype=np.float64)
    epsilon = np.zeros((n1, n2), dtype=np.float64)
    charge_product = np.zeros((n1, n2), dtype=np.float64)
    for i, bt1 in enumerate(bead_types1):
        for j, bt2 in enumerate(bead_types2):
            params = pair_params.get((bt1, bt2))
            if params is None:
                params = pair_params.get((bt2, bt1))
            if params is None:
                raise RuntimeError(f"Missing pair params for ({bt1}, {bt2})")
            sigma[i, j] = float(params["sigma_nm"])
            epsilon[i, j] = float(params["epsilon_kj_mol"])
            charge_product[i, j] = float(bead_charges1[i]) * float(bead_charges2[j])
    return sigma, epsilon, charge_product


def _cg_lipid_configurations(
    dirs: np.ndarray,
    bead_frame_angles: np.ndarray,
    ref_nm: np.ndarray,
) -> np.ndarray:
    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    configs = []
    refs = np.asarray(ref_nm, dtype=np.float64)
    for dir_vec in np.asarray(dirs, dtype=np.float64):
        r_base = _rotation_to_align_z_np(np.asarray(dir_vec, dtype=np.float64))
        for frame_angle in bead_frame_angles:
            rot = r_base @ _rotation_about_axis_np(z_axis, float(frame_angle))
            for ref_conf in refs:
                configs.append((rot @ ref_conf.T).T)
    return np.asarray(configs, dtype=np.float64)


def _compute_cg_pair_energy_samples_vectorized(
    r_nm: float,
    dirs1: np.ndarray,
    dirs2: np.ndarray,
    bead_frame_angles: np.ndarray,
    ref_nm1: np.ndarray,
    bead_types1: list,
    bead_charges1: list,
    pair_params: dict,
    dist_min_nm: float,
    ref_nm2: np.ndarray | None = None,
    bead_types2: list | None = None,
    bead_charges2: list | None = None,
    chunk_size: int = 16,
) -> np.ndarray:
    if ref_nm2 is None:
        ref_nm2 = ref_nm1
    if bead_types2 is None:
        bead_types2 = bead_types1
    if bead_charges2 is None:
        bead_charges2 = bead_charges1

    configs1 = _cg_lipid_configurations(dirs1, bead_frame_angles, ref_nm1)
    configs2 = _cg_lipid_configurations(dirs2, bead_frame_angles, ref_nm2)
    configs2 = configs2 + np.array([float(r_nm), 0.0, 0.0], dtype=np.float64)[None, None, :]

    sigma, epsilon, charge_product = _lipid_pair_param_matrices(
        bead_types1,
        bead_types2,
        bead_charges1,
        bead_charges2,
        pair_params,
    )
    sigma = sigma[None, None, :, :]
    epsilon = epsilon[None, None, :, :]
    charge_product = charge_product[None, None, :, :]
    charged = charge_product != 0.0

    out = []
    chunk_size = max(1, int(chunk_size))
    for start in range(0, configs1.shape[0], chunk_size):
        p1 = configs1[start : start + chunk_size]
        delta = configs2[None, :, None, :, :] - p1[:, None, :, None, :]
        dist = np.sqrt(np.sum(delta * delta, axis=-1))
        active = dist <= DRY_MARTINI_NONBONDED_CUTOFF_NM
        eff_dist = np.maximum(dist, float(dist_min_nm))
        sr = sigma / eff_dist
        sr2 = sr * sr
        sr6 = sr2 * sr2 * sr2
        energy = 4.0 * epsilon * (sr6 * sr6 - sr6)
        coulomb = COULOMB_K_DRY_KJ_NM * charge_product / eff_dist
        energy = np.where(active, energy, 0.0)
        energy = energy + np.where(active & charged & (dist > 1.0e-10), coulomb, 0.0)
        out.append(np.sum(energy, axis=(2, 3)).reshape(-1))
    return np.concatenate(out)


def _direction_with_dot_np(axis: np.ndarray, dot_value: float, phi: float) -> np.ndarray:
    """Return a unit vector v with dot(v, axis)=dot_value."""
    axis = np.asarray(axis, dtype=np.float64)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1e-12:
        raise ValueError("axis must be non-zero")
    axis = axis / axis_norm

    if abs(axis[0]) < 0.9:
        tangent1 = np.cross(axis, np.array([1.0, 0.0, 0.0], dtype=np.float64))
    else:
        tangent1 = np.cross(axis, np.array([0.0, 1.0, 0.0], dtype=np.float64))
    tangent1 /= np.linalg.norm(tangent1)
    tangent2 = np.cross(axis, tangent1)

    c = float(np.clip(dot_value, -1.0, 1.0))
    s = math.sqrt(max(0.0, 1.0 - c * c))
    return c * axis + s * (math.cos(phi) * tangent1 + math.sin(phi) * tangent2)


def _directions_with_dot_np(axis: np.ndarray, dot_grid: np.ndarray, phi_values: np.ndarray) -> np.ndarray:
    dirs = np.zeros((len(dot_grid), len(phi_values), 3), dtype=np.float64)
    for ia, dot_value in enumerate(dot_grid):
        for ip, phi in enumerate(phi_values):
            dirs[ia, ip] = _direction_with_dot_np(axis, float(dot_value), float(phi))
    return dirs


def _rotation_about_axis_np(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=np.float64)
    norm = float(np.linalg.norm(axis))
    if norm <= 1.0e-12:
        raise ValueError("rotation axis must be non-zero")
    axis = axis / norm
    x, y, z = axis
    c = math.cos(float(angle))
    s = math.sin(float(angle))
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=np.float64,
    )


def _rotate_points_about_axis_np(
    points: np.ndarray,
    axis: np.ndarray,
    angle: float,
    origin: np.ndarray,
) -> np.ndarray:
    rot = _rotation_about_axis_np(axis, angle)
    pts = np.asarray(points, dtype=np.float64)
    org = np.asarray(origin, dtype=np.float64)
    return org[None, :] + (rot @ (pts - org[None, :]).T).T


def _canonicalize_lipid_reference_to_z(ref_bead_positions_nm: np.ndarray) -> np.ndarray:
    ref = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    direction = ((ref[8] + ref[13]) * 0.5) - ref[0]
    direction /= max(float(np.linalg.norm(direction)), 1e-12)
    z_axis = direction
    if abs(z_axis[0]) < 0.99:
        x_axis = np.cross([1.0, 0.0, 0.0], z_axis)
    else:
        x_axis = np.cross([0.0, 1.0, 0.0], z_axis)
    x_axis /= max(float(np.linalg.norm(x_axis)), 1e-12)
    y_axis = np.cross(z_axis, x_axis)
    rot_local_to_ref = np.array([x_axis, y_axis, z_axis], dtype=np.float64).T
    return (rot_local_to_ref.T @ ref.T).T


def _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm: np.ndarray) -> np.ndarray:
    refs = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if refs.shape == (14, 3):
        refs = refs[None, :, :]
    if refs.ndim != 3 or refs.shape[1:] != (14, 3):
        raise ValueError(f"lipid reference ensemble must be (14, 3) or (n, 14, 3), got {refs.shape}")
    return np.asarray([_canonicalize_lipid_reference_to_z(ref) for ref in refs], dtype=np.float64)


def _reference_summary_geometry(ref_bead_positions_nm: np.ndarray) -> np.ndarray:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    return np.mean(refs, axis=0)


def _dopc_tail_extension_ang(ref_bead_positions_nm: np.ndarray) -> float:
    ref = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    head = ref[0]
    tail_mid = 0.5 * (ref[8] + ref[13])
    extension = float(np.linalg.norm(tail_mid - head))
    if extension <= 1.0e-12:
        raise ValueError("DOPC axial extension coordinate requires a finite head-tail axis")
    return extension * LENGTH_CONVERSION_A_PER_NM


def _dopc_tail_extension_series_ang(ref_bead_positions_nm: np.ndarray) -> np.ndarray:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    return np.asarray([_dopc_tail_extension_ang(ref) for ref in refs], dtype=np.float64)


def _retarget_lipid_reference_tail_extension_ang(
    ref_bead_positions_nm: np.ndarray,
    target_extension_ang: float,
) -> np.ndarray:
    ref = _canonicalize_lipid_reference_to_z(ref_bead_positions_nm)
    target = float(target_extension_ang)
    if not np.isfinite(target) or target <= 0.0:
        raise ValueError(f"Target tail extension must be positive and finite, got {target_extension_ang!r}")
    current = _dopc_tail_extension_ang(ref)
    delta_ang = target - current
    if abs(delta_ang) <= 1.0e-12:
        return ref.copy()
    axis = ((ref[8] + ref[13]) * 0.5) - ref[0]
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1.0e-12:
        raise ValueError("Cannot retarget lipid tail extension for a degenerate head-tail axis")
    tail_weights = _tail_axial_compaction_weights(ref)
    retargeted = _apply_tail_axial_shift(
        ref,
        axis / axis_norm,
        -(delta_ang * ANGSTROM_TO_NM),
        tail_weights,
    )
    return _canonicalize_lipid_reference_to_z(retargeted)


def _select_reference_ensemble_representatives(
    ref_bead_positions_nm: np.ndarray,
    representative_count: int,
) -> dict:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    count = max(1, int(representative_count))
    if refs.shape[0] <= count:
        compaction = _dopc_tail_extension_series_ang(refs)
        return {
            "representative_refs_nm": refs.copy(),
            "representative_indices": np.arange(refs.shape[0], dtype=np.int32),
            "representative_compaction_ang": compaction.astype(np.float64),
            "pool_size": int(refs.shape[0]),
        }

    compaction = _dopc_tail_extension_series_ang(refs)
    order = np.argsort(compaction)
    # Use equal-probability bin centers so the reduced ensemble represents the
    # bulk shell distribution instead of forcing rare compaction extrema into
    # every fit.
    rank_targets = ((np.arange(count, dtype=np.float64) + 0.5) * float(order.size) / float(count)) - 0.5
    selected_ranks: list[int] = []
    used_ranks: set[int] = set()
    for target in rank_targets:
        base_rank = int(round(float(target)))
        for delta in range(order.size):
            candidate_ranks = (base_rank,) if delta == 0 else (base_rank - delta, base_rank + delta)
            chosen = None
            for cand in candidate_ranks:
                if 0 <= cand < order.size and cand not in used_ranks:
                    chosen = cand
                    break
            if chosen is not None:
                used_ranks.add(chosen)
                selected_ranks.append(chosen)
                break
    if len(selected_ranks) != count:
        raise RuntimeError(
            f"Failed to select {count} representative reference conformers from pool of {refs.shape[0]}"
        )
    representative_indices = order[np.asarray(selected_ranks, dtype=np.int32)]
    representative_indices = representative_indices[
        np.argsort(compaction[representative_indices], kind="stable")
    ]
    return {
        "representative_refs_nm": refs[representative_indices].copy(),
        "representative_indices": representative_indices.astype(np.int32),
        "representative_compaction_ang": compaction[representative_indices].astype(np.float64),
        "pool_size": int(refs.shape[0]),
    }


def _select_compaction_state_representatives(
    ref_bead_positions_nm: np.ndarray,
    compaction_ang: np.ndarray,
    representative_count: int,
) -> dict:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    compaction = np.asarray(compaction_ang, dtype=np.float64)
    if refs.shape[0] != compaction.size:
        raise ValueError("Compaction coordinate count must match reference conformer count")
    if compaction.size < 4:
        raise ValueError("Compaction-state redesign requires at least four isolated conformers")

    q25 = float(np.quantile(compaction, 0.25))
    q75 = float(np.quantile(compaction, 0.75))
    # For the axial head-to-tail extension coordinate, bilayer-like compact
    # tails correspond to larger extension, not smaller.
    compact_mask = compaction >= q75
    extended_mask = compaction <= q25
    if int(np.count_nonzero(compact_mask)) < 1 or int(np.count_nonzero(extended_mask)) < 1:
        raise ValueError("Failed to identify compact and extended isolated DOPC states")

    compact_center = float(np.mean(compaction[compact_mask]))
    extended_center = float(np.mean(compaction[extended_mask]))
    if not (compact_center > extended_center):
        raise ValueError(
            f"Expected compact center > extended center, got {compact_center:.4f} <= {extended_center:.4f}"
        )

    def select_subset(mask: np.ndarray, center: float) -> np.ndarray:
        idx = np.where(mask)[0]
        order = idx[np.argsort(np.abs(compaction[idx] - center))]
        count = max(1, min(int(representative_count), order.size))
        return refs[order[:count]]

    compact_refs = select_subset(compact_mask, compact_center)
    extended_refs = select_subset(extended_mask, extended_center)
    compact_probability = float(np.mean(compact_mask))
    return {
        "compact_refs_nm": compact_refs,
        "extended_refs_nm": extended_refs,
        "compact_center_ang": compact_center,
        "extended_center_ang": extended_center,
        "compact_probability": compact_probability,
        "pool_size": int(refs.shape[0]),
        "compact_pool_size": int(np.count_nonzero(compact_mask)),
        "extended_pool_size": int(np.count_nonzero(extended_mask)),
        "compaction_min_ang": float(np.min(compaction)),
        "compaction_max_ang": float(np.max(compaction)),
        "compaction_p05_ang": float(np.quantile(compaction, 0.05)),
        "compaction_p95_ang": float(np.quantile(compaction, 0.95)),
    }


def _select_compaction_state_representatives_by_center(
    ref_bead_positions_nm: np.ndarray,
    compact_center_ang: float,
    extended_center_ang: float,
    representative_count: int,
    compact_probability: float | None = None,
) -> dict:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    compaction = _dopc_tail_extension_series_ang(refs)
    if refs.shape[0] < 2:
        raise ValueError("Center-matched compaction-state redesign requires at least two conformers")

    compact_center = float(compact_center_ang)
    extended_center = float(extended_center_ang)
    if not (compact_center > extended_center):
        raise ValueError(
            f"Expected compact center > extended center, got {compact_center:.4f} <= {extended_center:.4f}"
        )

    midpoint = 0.5 * (compact_center + extended_center)
    compact_mask = compaction >= midpoint
    extended_mask = compaction <= midpoint
    if int(np.count_nonzero(compact_mask)) < 1 or int(np.count_nonzero(extended_mask)) < 1:
        median = float(np.median(compaction))
        compact_mask = compaction >= median
        extended_mask = compaction <= median
    if int(np.count_nonzero(compact_mask)) < 1 or int(np.count_nonzero(extended_mask)) < 1:
        raise ValueError(
            "Failed to identify both compact and extended conformers in the reference ensemble"
        )

    def select_subset(mask: np.ndarray, center: float) -> tuple[np.ndarray, np.ndarray]:
        idx = np.where(mask)[0]
        order = idx[np.argsort(np.abs(compaction[idx] - center), kind="stable")]
        count = max(1, min(int(representative_count), order.size))
        selected = order[:count]
        return refs[selected], compaction[selected]

    compact_refs, compact_selected = select_subset(compact_mask, compact_center)
    extended_refs, extended_selected = select_subset(extended_mask, extended_center)
    if compact_probability is None:
        compact_probability = float(np.mean(compaction >= midpoint))

    return {
        "compact_refs_nm": compact_refs,
        "extended_refs_nm": extended_refs,
        "compact_center_ang": compact_center,
        "extended_center_ang": extended_center,
        "compact_probability": float(compact_probability),
        "pool_size": int(refs.shape[0]),
        "compact_pool_size": int(np.count_nonzero(compact_mask)),
        "extended_pool_size": int(np.count_nonzero(extended_mask)),
        "compaction_min_ang": float(np.min(compaction)),
        "compaction_max_ang": float(np.max(compaction)),
        "compaction_p05_ang": float(np.quantile(compaction, 0.05)),
        "compaction_p95_ang": float(np.quantile(compaction, 0.95)),
        "compact_selected_compaction_ang": np.asarray(compact_selected, dtype=np.float64),
        "extended_selected_compaction_ang": np.asarray(extended_selected, dtype=np.float64),
    }


def _augment_compaction_states_with_compressed_branch(
    ref_bead_positions_nm: np.ndarray,
    compaction_states: dict,
    representative_count: int,
    supplemental_refs_nm: np.ndarray | None = None,
    sample_compaction_ang: np.ndarray | None = None,
) -> dict:
    compact_center = float(compaction_states["compact_center_ang"])
    sample_values = None
    if sample_compaction_ang is not None:
        sample_values = np.asarray(sample_compaction_ang, dtype=np.float64)
        sample_values = sample_values[np.isfinite(sample_values)]
        sample_upper = sample_values[sample_values > compact_center + 1.0e-6]
        if sample_upper.size >= max(4, int(representative_count)):
            compressed_center = float(np.mean(sample_upper))
            compact_refs = np.asarray(compaction_states["compact_refs_nm"], dtype=np.float64)
            count = max(1, min(int(representative_count), compact_refs.shape[0]))
            selected = np.arange(count, dtype=np.int32)
            augmented = dict(compaction_states)
            augmented["compressed_refs_nm"] = np.asarray(
                [
                    _retarget_lipid_reference_tail_extension_ang(ref, compressed_center)
                    for ref in compact_refs[selected]
                ],
                dtype=np.float64,
            )
            augmented["compressed_center_ang"] = compressed_center
            augmented["compressed_pool_size"] = int(sample_upper.size)
            augmented["compressed_selected_compaction_ang"] = np.asarray(
                [compressed_center] * count,
                dtype=np.float64,
            )
            augmented["compressed_center_source"] = "target_compaction_high_population_mean"
            return augmented

    pool_refs = []
    pool_compaction = []

    def append_upper_tail(refs_nm: np.ndarray | None):
        if refs_nm is None:
            return
        refs = _canonicalize_lipid_reference_ensemble_to_z(refs_nm)
        compaction = _dopc_tail_extension_series_ang(refs)
        if refs.shape[0] != compaction.size:
            raise ValueError("Compressed compaction-state support requires matched reference conformers")
        above_compact = compaction > compact_center + 1.0e-6
        if not np.any(above_compact):
            return
        pool_refs.append(refs[above_compact])
        pool_compaction.append(compaction[above_compact])

    append_upper_tail(ref_bead_positions_nm)
    append_upper_tail(supplemental_refs_nm)
    if not pool_compaction:
        return dict(compaction_states)

    refs = np.concatenate(pool_refs, axis=0)
    compaction = np.concatenate(pool_compaction, axis=0)
    compressed_center = float(np.median(compaction))
    order = np.argsort(np.abs(compaction - compressed_center), kind="stable")
    count = max(1, min(int(representative_count), order.size))
    selected = order[:count]
    augmented = dict(compaction_states)
    augmented["compressed_refs_nm"] = refs[selected]
    augmented["compressed_center_ang"] = compressed_center
    augmented["compressed_pool_size"] = int(compaction.size)
    augmented["compressed_selected_compaction_ang"] = np.asarray(
        compaction[selected],
        dtype=np.float64,
    )
    return augmented


def _fit_compaction_self_pmf(
    compaction_ang: np.ndarray,
    temperature_upside: float,
    n_bin: int = 12,
    smooth: float = 0.01,
) -> dict:
    values = np.asarray(compaction_ang, dtype=np.float64)
    if values.ndim != 1 or values.size < 4:
        raise ValueError("Compaction self PMF requires at least four samples")
    kbt_kj_mol = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt_kj_mol <= 0.0:
        raise ValueError("Compaction self PMF temperature must be positive")

    lo = float(np.quantile(values, 0.05))
    hi = float(np.quantile(values, 0.95))
    span = max(hi - lo, 1.0e-3)
    margin = max(0.15 * span, 0.1)
    coord_min = float(min(np.min(values), lo - margin))
    coord_max = float(max(np.max(values), hi + margin))
    n_center = max(8, int(n_bin))
    centers = np.linspace(coord_min, coord_max, n_center, dtype=np.float64)
    spacing = float(centers[1] - centers[0]) if centers.size > 1 else 0.25
    edges = _grid_edges_from_centers(centers, lo=centers[0] - 0.5 * spacing, hi=centers[-1] + 0.5 * spacing)
    hist, _ = np.histogram(values, bins=edges)
    prob = hist.astype(np.float64) / max(1.0, float(np.sum(hist)))
    floor = 1.0 / max(float(values.size) * 1000.0, 1.0)
    pmf_kj_mol = -kbt_kj_mol * np.log(np.maximum(prob, floor))
    pmf_kj_mol -= float(np.min(pmf_kj_mol))

    n_knot = int(centers.size) + 2
    knot_vector = np.zeros(n_knot + 4, dtype=np.float64)
    knot_vector[4:-4] = np.arange(1, n_knot - 3, dtype=np.float64)
    knot_vector[-4:] = knot_vector[-5]
    coord = 1.0 + (centers - centers[0]) / max(spacing, 1.0e-6)
    coeff = _fit_radial_bspline(
        coord,
        pmf_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        knot_vector,
        smooth=float(smooth),
    )
    coeff = np.clip(
        coeff,
        float(np.min(pmf_kj_mol)) / ENERGY_CONVERSION_KJ_PER_EUP,
        float(np.max(pmf_kj_mol)) / ENERGY_CONVERSION_KJ_PER_EUP,
    )
    min_idx = int(np.argmin(pmf_kj_mol))
    if 0 < min_idx < pmf_kj_mol.size - 1:
        curvature_kj_mol_a2 = (
            float(pmf_kj_mol[min_idx - 1])
            - 2.0 * float(pmf_kj_mol[min_idx])
            + float(pmf_kj_mol[min_idx + 1])
        ) / max(spacing * spacing, 1.0e-12)
    else:
        curvature_kj_mol_a2 = kbt_kj_mol / max(spacing * spacing, 1.0e-12)
    effective_stiffness = max(
        curvature_kj_mol_a2 / ENERGY_CONVERSION_KJ_PER_EUP,
        1.0e-3,
    )
    return {
        "self_coeff_eup": coeff,
        "coord_min_ang": float(centers[0]),
        "coord_max_ang": float(centers[-1]),
        "coord_spacing_ang": float(spacing),
        "n_knot": int(n_knot),
        "pmf_centers_ang": centers,
        "pmf_values_kj_mol": pmf_kj_mol,
        "effective_stiffness_eup_a2": float(effective_stiffness),
        "sample_count": int(values.size),
    }


def _single_cgl_state_reparameterized_tables(
    delta_extended: np.ndarray,
    delta_compact: np.ndarray,
    grid_extended_kj_mol: np.ndarray | None,
    grid_compact_kj_mol: np.ndarray | None,
    source_extended_center_ang: float,
    source_compact_center_ang: float,
    target_extended_center_ang: float,
    target_compact_center_ang: float,
) -> dict:
    denom = float(source_compact_center_ang) - float(source_extended_center_ang)
    if abs(denom) <= 1.0e-12:
        raise ValueError(
            "Single-CGL endpoint reparameterization requires distinct source centers, got "
            f"{source_extended_center_ang:.6f} / {source_compact_center_ang:.6f}"
        )

    def interp(coord_ang: float) -> tuple[float, float]:
        compact_weight = (float(coord_ang) - float(source_extended_center_ang)) / denom
        return 1.0 - compact_weight, compact_weight

    w_ext_ext, w_ext_comp = interp(target_extended_center_ang)
    w_cmp_ext, w_cmp_comp = interp(target_compact_center_ang)
    delta_extended = np.asarray(delta_extended, dtype=np.float64)
    delta_compact = np.asarray(delta_compact, dtype=np.float64)
    payload = {
        "delta_extended": np.asarray(
            w_ext_ext * delta_extended + w_ext_comp * delta_compact,
            dtype=np.float32,
        ),
        "delta_compact": np.asarray(
            w_cmp_ext * delta_extended + w_cmp_comp * delta_compact,
            dtype=np.float32,
        ),
    }
    if grid_extended_kj_mol is not None and grid_compact_kj_mol is not None:
        grid_extended = np.asarray(grid_extended_kj_mol, dtype=np.float64)
        grid_compact = np.asarray(grid_compact_kj_mol, dtype=np.float64)
        payload["grid_extended_kj_mol"] = np.asarray(
            w_ext_ext * grid_extended + w_ext_comp * grid_compact,
            dtype=np.float32,
        )
        payload["grid_compact_kj_mol"] = np.asarray(
            w_cmp_ext * grid_extended + w_cmp_comp * grid_compact,
            dtype=np.float32,
        )
    return payload


def _pair_compaction_reparameterized_tables(
    delta_extended_extended: np.ndarray,
    delta_extended_compact: np.ndarray,
    delta_compact_compact: np.ndarray,
    grid_extended_extended_kj_mol: np.ndarray | None,
    grid_extended_compact_kj_mol: np.ndarray | None,
    grid_compact_compact_kj_mol: np.ndarray | None,
    grid_average_kj_mol: np.ndarray | None,
    source_extended_center_ang: float,
    source_compact_center_ang: float,
    target_extended_center_ang: float,
    target_compact_center_ang: float,
    delta_extended_compressed: np.ndarray | None = None,
    delta_compact_compressed: np.ndarray | None = None,
    delta_compressed_compressed: np.ndarray | None = None,
    grid_extended_compressed_kj_mol: np.ndarray | None = None,
    grid_compact_compressed_kj_mol: np.ndarray | None = None,
    grid_compressed_compressed_kj_mol: np.ndarray | None = None,
    source_compressed_center_ang: float | None = None,
    target_compressed_center_ang: float | None = None,
) -> dict:
    def single_mix_weights(
        coord_ang: float,
        extended_center_ang: float,
        compact_center_ang: float,
        compressed_center_ang: float | None = None,
    ) -> tuple[float, float, float]:
        denom_lo = float(compact_center_ang) - float(extended_center_ang)
        if abs(denom_lo) <= 1.0e-12:
            raise ValueError(
                "Pair-compaction reparameterization requires distinct extended/compact centers, got "
                f"{extended_center_ang:.6f} / {compact_center_ang:.6f}"
            )
        if (
            compressed_center_ang is None
            or not np.isfinite(float(compressed_center_ang))
            or float(compressed_center_ang) <= float(compact_center_ang) + 1.0e-6
        ):
            s = (float(coord_ang) - float(extended_center_ang)) / denom_lo
            return 1.0 - s, s, 0.0
        if float(coord_ang) <= float(compact_center_ang):
            s = (float(coord_ang) - float(extended_center_ang)) / denom_lo
            return 1.0 - s, s, 0.0
        denom_hi = float(compressed_center_ang) - float(compact_center_ang)
        if abs(denom_hi) <= 1.0e-12:
            raise ValueError(
                "Pair-compaction reparameterization requires distinct compact/compressed centers, got "
                f"{compact_center_ang:.6f} / {compressed_center_ang:.6f}"
            )
        t = (float(coord_ang) - float(compact_center_ang)) / denom_hi
        return 0.0, 1.0 - t, t

    has_source_compressed = (
        delta_extended_compressed is not None
        and delta_compact_compressed is not None
        and delta_compressed_compressed is not None
        and source_compressed_center_ang is not None
        and np.isfinite(float(source_compressed_center_ang))
        and float(source_compressed_center_ang) > float(source_compact_center_ang) + 1.0e-6
    )
    has_any_source_compressed_payload = any(
        value is not None
        for value in (
            delta_extended_compressed,
            delta_compact_compressed,
            delta_compressed_compressed,
            grid_extended_compressed_kj_mol,
            grid_compact_compressed_kj_mol,
            grid_compressed_compressed_kj_mol,
        )
    )
    if has_any_source_compressed_payload and not has_source_compressed:
        raise ValueError(
            "Pair-compaction reparameterization requires a complete compressed-state source payload"
        )
    has_target_compressed = (
        target_compressed_center_ang is not None
        and np.isfinite(float(target_compressed_center_ang))
        and float(target_compressed_center_ang) > float(target_compact_center_ang) + 1.0e-6
    )
    if has_target_compressed and not has_source_compressed:
        raise ValueError(
            "Pair-compaction reparameterization cannot target a compressed state without compressed source tensors"
        )

    delta_ee = np.asarray(delta_extended_extended, dtype=np.float64)
    delta_ec = np.asarray(delta_extended_compact, dtype=np.float64)
    delta_cc = np.asarray(delta_compact_compact, dtype=np.float64)
    delta_ex = (
        np.asarray(delta_extended_compressed, dtype=np.float64)
        if delta_extended_compressed is not None
        else None
    )
    delta_cx = (
        np.asarray(delta_compact_compressed, dtype=np.float64)
        if delta_compact_compressed is not None
        else None
    )
    delta_xx = (
        np.asarray(delta_compressed_compressed, dtype=np.float64)
        if delta_compressed_compressed is not None
        else None
    )

    def mix(
        coord_i_ang: float,
        coord_j_ang: float,
        ee: np.ndarray,
        ec: np.ndarray,
        cc: np.ndarray,
        ex: np.ndarray | None = None,
        cx: np.ndarray | None = None,
        xx: np.ndarray | None = None,
    ) -> np.ndarray:
        wi_ext, wi_comp, wi_x = single_mix_weights(
            coord_i_ang,
            source_extended_center_ang,
            source_compact_center_ang,
            source_compressed_center_ang if has_source_compressed else None,
        )
        wj_ext, wj_comp, wj_x = single_mix_weights(
            coord_j_ang,
            source_extended_center_ang,
            source_compact_center_ang,
            source_compressed_center_ang if has_source_compressed else None,
        )
        out = (
            (wi_ext * wj_ext) * ee
            + (wi_ext * wj_comp + wi_comp * wj_ext) * ec
            + (wi_comp * wj_comp) * cc
        )
        if has_source_compressed:
            if ex is None or cx is None or xx is None:
                raise ValueError("Compressed-state pair mix requires EX/CX/XX tensors")
            out = out + (
                (wi_ext * wj_x + wi_x * wj_ext) * ex
                + (wi_comp * wj_x + wi_x * wj_comp) * cx
                + (wi_x * wj_x) * xx
            )
        return out

    payload = {
        "delta_extended_extended": np.asarray(
            mix(
                target_extended_center_ang,
                target_extended_center_ang,
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        ),
        "delta_extended_compact": np.asarray(
            mix(
                target_extended_center_ang,
                target_compact_center_ang,
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        ),
        "delta_compact_compact": np.asarray(
            mix(
                target_compact_center_ang,
                target_compact_center_ang,
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        ),
    }
    if has_target_compressed:
        payload["delta_extended_compressed"] = np.asarray(
            mix(
                target_extended_center_ang,
                float(target_compressed_center_ang),
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        )
        payload["delta_compact_compressed"] = np.asarray(
            mix(
                target_compact_center_ang,
                float(target_compressed_center_ang),
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        )
        payload["delta_compressed_compressed"] = np.asarray(
            mix(
                float(target_compressed_center_ang),
                float(target_compressed_center_ang),
                delta_ee,
                delta_ec,
                delta_cc,
                delta_ex,
                delta_cx,
                delta_xx,
            ),
            dtype=np.float32,
        )

    if (
        grid_extended_extended_kj_mol is not None
        and grid_extended_compact_kj_mol is not None
        and grid_compact_compact_kj_mol is not None
    ):
        grid_ee = np.asarray(grid_extended_extended_kj_mol, dtype=np.float64)
        grid_ec = np.asarray(grid_extended_compact_kj_mol, dtype=np.float64)
        grid_cc = np.asarray(grid_compact_compact_kj_mol, dtype=np.float64)
        grid_ex = (
            np.asarray(grid_extended_compressed_kj_mol, dtype=np.float64)
            if grid_extended_compressed_kj_mol is not None
            else None
        )
        grid_cx = (
            np.asarray(grid_compact_compressed_kj_mol, dtype=np.float64)
            if grid_compact_compressed_kj_mol is not None
            else None
        )
        grid_xx = (
            np.asarray(grid_compressed_compressed_kj_mol, dtype=np.float64)
            if grid_compressed_compressed_kj_mol is not None
            else None
        )
        payload["grid_extended_extended_kj_mol"] = np.asarray(
            mix(
                target_extended_center_ang,
                target_extended_center_ang,
                grid_ee,
                grid_ec,
                grid_cc,
                grid_ex,
                grid_cx,
                grid_xx,
            ),
            dtype=np.float32,
        )
        payload["grid_extended_compact_kj_mol"] = np.asarray(
            mix(
                target_extended_center_ang,
                target_compact_center_ang,
                grid_ee,
                grid_ec,
                grid_cc,
                grid_ex,
                grid_cx,
                grid_xx,
            ),
            dtype=np.float32,
        )
        payload["grid_compact_compact_kj_mol"] = np.asarray(
            mix(
                target_compact_center_ang,
                target_compact_center_ang,
                grid_ee,
                grid_ec,
                grid_cc,
                grid_ex,
                grid_cx,
                grid_xx,
            ),
            dtype=np.float32,
        )
        if has_target_compressed and grid_ex is not None and grid_cx is not None and grid_xx is not None:
            payload["grid_extended_compressed_kj_mol"] = np.asarray(
                mix(
                    target_extended_center_ang,
                    float(target_compressed_center_ang),
                    grid_ee,
                    grid_ec,
                    grid_cc,
                    grid_ex,
                    grid_cx,
                    grid_xx,
                ),
                dtype=np.float32,
            )
            payload["grid_compact_compressed_kj_mol"] = np.asarray(
                mix(
                    target_compact_center_ang,
                    float(target_compressed_center_ang),
                    grid_ee,
                    grid_ec,
                    grid_cc,
                    grid_ex,
                    grid_cx,
                    grid_xx,
                ),
                dtype=np.float32,
            )
            payload["grid_compressed_compressed_kj_mol"] = np.asarray(
                mix(
                    float(target_compressed_center_ang),
                    float(target_compressed_center_ang),
                    grid_ee,
                    grid_ec,
                    grid_cc,
                    grid_ex,
                    grid_cx,
                    grid_xx,
                ),
                dtype=np.float32,
            )
        if grid_average_kj_mol is not None:
            payload["grid_average_kj_mol"] = np.asarray(grid_average_kj_mol, dtype=np.float32)
    return payload


def _unit_perpendicular(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=np.float64)
    trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if abs(float(np.dot(vector, trial))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    perp = trial - vector * float(np.dot(vector, trial))
    norm = float(np.linalg.norm(perp))
    if norm <= 1.0e-12:
        raise ValueError("Cannot construct perpendicular vector for ITP DOPC reference")
    return perp / norm


def _place_bonded_bead(
    positions: np.ndarray,
    atom_i: int,
    atom_j: int,
    atom_k: int,
    bond_length_nm: float,
    angle_deg: float,
    bend_sign: float = 1.0,
) -> None:
    p_i = positions[int(atom_i)]
    p_j = positions[int(atom_j)]
    u = p_i - p_j
    norm = float(np.linalg.norm(u))
    if norm <= 1.0e-12:
        raise ValueError("Cannot place ITP DOPC reference bead from coincident parent atoms")
    u /= norm
    theta = math.radians(float(angle_deg))
    perp = _unit_perpendicular(u) * (1.0 if bend_sign >= 0.0 else -1.0)
    direction = math.cos(theta) * u + math.sin(theta) * perp
    positions[int(atom_k)] = p_j + float(bond_length_nm) * direction


def _build_dopc_reference_from_itp_topology(dopc: dict) -> np.ndarray:
    """Construct an isolated DOPC reference from ITP bonded geometry only.

    The reference is a deterministic tree embedding of the DOPC topology.  It
    avoids using a bilayer-template PDB conformation while still preserving the
    dry-MARTINI bond lengths and available bond-angle preferences.
    """
    atom_names = list(dopc["atom_names"])
    if len(atom_names) != 14:
        raise ValueError(f"DOPC topology must have 14 beads, got {len(atom_names)}")
    bonds = {(min(i, j), max(i, j)): (float(r0), float(k)) for i, j, r0, k in dopc["bonds"]}
    angles = {(i, j, k): (float(theta), float(force_k)) for i, j, k, theta, force_k in dopc["angles"]}

    def bond(i: int, j: int) -> float:
        key = (min(i, j), max(i, j))
        if key not in bonds:
            raise ValueError(f"Missing DOPC bond {atom_names[i]}-{atom_names[j]}")
        return bonds[key][0]

    def angle(i: int, j: int, k: int, default: float = 180.0) -> float:
        if (i, j, k) in angles:
            return angles[(i, j, k)][0]
        if (k, j, i) in angles:
            return angles[(k, j, i)][0]
        return float(default)

    pos = np.full((14, 3), np.nan, dtype=np.float64)
    # Head/glycerol scaffold.  The exact global orientation is irrelevant; the
    # table builder canonicalizes the final head-to-tail axis.
    pos[2] = np.array([0.0, 0.0, 0.0], dtype=np.float64)  # GL1
    pos[1] = np.array([0.0, 0.0, -bond(1, 2)], dtype=np.float64)  # PO4
    pos[0] = pos[1] + np.array([0.0, 0.0, -bond(0, 1)], dtype=np.float64)  # NC3
    _place_bonded_bead(pos, 1, 2, 3, bond(2, 3), angle(1, 2, 3), bend_sign=1.0)  # GL2
    _place_bonded_bead(pos, 1, 2, 4, bond(2, 4), angle(1, 2, 4), bend_sign=-1.0)  # C1A

    # A tail: GL1-C1A-C2A-D3A-C4A-C5A.
    chain_a = [2, 4, 5, 6, 7, 8]
    for prev, cur, nxt in zip(chain_a[:-2], chain_a[1:-1], chain_a[2:]):
        _place_bonded_bead(pos, prev, cur, nxt, bond(cur, nxt), angle(prev, cur, nxt), bend_sign=-1.0)

    # B tail starts from GL2.  The ITP has no GL1-GL2-C1B angle, so place the
    # branch roughly parallel to the A tail before following explicit angles.
    tail_axis = pos[4] - pos[2]
    tail_axis /= float(np.linalg.norm(tail_axis))
    pos[9] = pos[3] + bond(3, 9) * tail_axis
    chain_b = [3, 9, 10, 11, 12, 13]
    for prev, cur, nxt in zip(chain_b[:-2], chain_b[1:-1], chain_b[2:]):
        _place_bonded_bead(pos, prev, cur, nxt, bond(cur, nxt), angle(prev, cur, nxt), bend_sign=1.0)

    if not np.all(np.isfinite(pos)):
        raise ValueError("Failed to construct finite ITP-derived DOPC reference")
    pos -= np.mean(pos, axis=0)
    return pos


def sample_isolated_dopc_bonded_conformers(
    dopc: dict,
    lipids_itp_path: Path,
    pair_params: dict,
    conformer_count: int = 2,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    seed: int = 1729,
    mc_burnin_steps: int = 5000,
    mc_steps_per_conformer: int = 2000,
    proposal_sigma_nm: float = 0.025,
) -> np.ndarray:
    """Sample isolated DOPC conformers from bonded and intramolecular nonbonded terms.

    This is intentionally protein-free and bilayer-free. It broadens the
    direct constituent projection with dry-MARTINI bonded flexibility plus the
    lipid's own intramolecular nonbonded interactions, while avoiding bilayer
    PMFs, IBI targets, force matching, or template conformers.
    """
    conformer_count = max(1, int(conformer_count))
    base = _build_dopc_reference_from_itp_topology(dopc)

    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt <= 0.0:
        raise ValueError("isolated DOPC conformer temperature must be positive")
    burnin_steps = max(0, int(mc_burnin_steps))
    steps_per_conformer = max(1, int(mc_steps_per_conformer))
    sigma = float(proposal_sigma_nm)
    if sigma <= 0.0:
        raise ValueError("isolated DOPC conformer proposal sigma must be positive")

    rng = np.random.default_rng(int(seed))
    current = np.asarray(base, dtype=np.float64).copy()
    current -= current.mean(axis=0, keepdims=True)
    current_energy = (
        compute_dopc_bonded_energy(current, lipids_itp_path)
        + compute_dopc_intramolecular_nonbonded_energy(
            current,
            bead_types=list(dopc["bead_types"]),
            bead_charges=list(dopc["bead_charges"]),
            pair_params=pair_params,
            lipids_itp_path=lipids_itp_path,
        )
    )
    accepted = 0
    attempted = 0

    def step_once() -> None:
        nonlocal accepted, attempted, current, current_energy
        attempted += 1
        trial = current.copy()
        bead = int(rng.integers(0, trial.shape[0]))
        trial[bead] += rng.normal(0.0, sigma, size=3)
        trial -= trial.mean(axis=0, keepdims=True)
        trial_energy = (
            compute_dopc_bonded_energy(trial, lipids_itp_path)
            + compute_dopc_intramolecular_nonbonded_energy(
                trial,
                bead_types=list(dopc["bead_types"]),
                bead_charges=list(dopc["bead_charges"]),
                pair_params=pair_params,
                lipids_itp_path=lipids_itp_path,
            )
        )
        delta = trial_energy - current_energy
        if delta <= 0.0 or float(rng.random()) < math.exp(-delta / kbt):
            current = trial
            current_energy = trial_energy
            accepted += 1

    for _ in range(burnin_steps):
        step_once()

    refs = []
    while len(refs) < conformer_count:
        for _ in range(steps_per_conformer):
            step_once()
        refs.append(current.copy())

    acceptance = accepted / max(1, attempted)
    print(
        "  Sampled isolated DOPC intramolecular conformers: "
        f"{len(refs)} conformer(s), burn-in={burnin_steps}, "
        f"acceptance={acceptance:.3f}, T={temperature_upside:.4f} Upside, "
        f"seed={int(seed)}"
    )
    return _canonicalize_lipid_reference_ensemble_to_z(np.asarray(refs, dtype=np.float64))


def _derive_dopc_cg_params_from_reference_ensemble(
    ref_bead_positions_nm: np.ndarray,
    bead_types: list,
    pair_params: dict,
    bead_masses_g_mol: list | None,
    bonds: list | None,
) -> dict:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    per_conf = [
        derive_dopc_cg_params(
            ref_bead_positions_nm=ref,
            bead_types=bead_types,
            pair_params=pair_params,
            bead_masses_g_mol=bead_masses_g_mol,
            bonds=bonds,
            energy_conversion_kj_per_eup=ENERGY_CONVERSION_KJ_PER_EUP,
            length_conversion_ang_per_nm=LENGTH_CONVERSION_A_PER_NM,
        )
        for ref in refs
    ]
    out = dict(per_conf[0])
    mean_attrs = (
        "orientation_length_ang",
        "orientation_bond_fc_eup_a2",
        "orientation_projected_bond_fc_eup_a2",
        "transverse_inertia_g_mol_a2",
        "head_tail_span_ang",
        "tail_projection_ang",
    )
    for name in mean_attrs:
        out[name] = float(np.mean([float(p[name]) for p in per_conf]))
    max_attrs = (
        "max_axis_radius_ang",
        "max_perp_radius_ang",
    )
    for name in max_attrs:
        out[name] = float(np.max([float(p[name]) for p in per_conf]))
    out["orientation_mass_g_mol"] = (
        float(out["transverse_inertia_g_mol_a2"])
        / max(float(out["orientation_length_ang"]) ** 2, 1.0e-12)
    )
    out["conformer_count"] = int(refs.shape[0])
    out["orientation_length_source"] = "ensemble_mean_DOPC_COM_to_tail_midpoint_projection"
    out["orientation_mass_source"] = "ensemble_mean_DOPC_transverse_rotational_inertia"
    out["orientation_bond_fc_source"] = "ensemble_mean_2x_projected_DOPC_bond_stiffness"
    out["orientation_projected_bond_fc_source"] = "ensemble_mean_sum_DOPC_bond_k_projected_on_orientation_axis"
    return out


def _cg_lipid_conformer_pair_indices(n_conformer: int) -> list[tuple[int, int]]:
    n_conformer = int(n_conformer)
    if n_conformer < 1:
        raise ValueError("n_conformer must be positive")
    return [(i, j) for i in range(n_conformer) for j in range(n_conformer)]


def _decode_h5_string_array(values) -> list[str]:
    return [
        x.decode("ascii", errors="ignore").strip() if isinstance(x, bytes) else str(x).strip()
        for x in values
    ]


def _read_output_box_lengths_ang(h5: h5py.File) -> np.ndarray | None:
    if "/output/box" in h5 and h5["/output/box"].shape[0] > 0:
        box = np.asarray(h5["/output/box"][-1], dtype=np.float64).reshape(-1)
        if box.size >= 3 and np.all(box[:3] > 0.0):
            return box[:3]
    for path in (
        "/input/potential/martini_potential",
        "/input/potential/cg_lipid_pair",
        "/input/potential/dist_spring",
    ):
        if path in h5:
            attrs = h5[path].attrs
            if all(name in attrs for name in ("x_len", "y_len", "z_len")):
                box = np.asarray([attrs["x_len"], attrs["y_len"], attrs["z_len"]], dtype=np.float64)
                if np.all(box > 0.0):
                    return box
    return None


def _extract_dopc_conformer_pool_from_upside_h5(
    up_file: Path,
    frame_start_fraction: float = 0.5,
    min_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_distance_ang: float = 0.0,
) -> np.ndarray:
    """Extract a pooled set of centered 14-bead DOPC conformers from one Upside trajectory."""
    up_file = Path(up_file).expanduser().resolve()
    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        atom_names = _decode_h5_string_array(h5["/input/atom_names"][:])
        particle_class = (
            _decode_h5_string_array(h5["/input/particle_class"][:])
            if "/input/particle_class" in h5
            else [""] * len(atom_names)
        )
        molecule_ids = np.asarray(h5["/input/molecule_ids"][:], dtype=np.int64)
        pos = h5["/output/pos"]
        n_frame = int(pos.shape[0])
        box_ang = _read_output_box_lengths_ang(h5)

        lipid_molecules = []
        for mol_id in sorted(set(int(x) for x in molecule_ids)):
            idx = np.where(molecule_ids == mol_id)[0]
            if idx.size != 14:
                continue
            if any(particle_class[i] and particle_class[i] != "LIPID" for i in idx):
                continue
            names = [atom_names[i] for i in idx]
            if names[:4] != ["NC3", "PO4", "GL1", "GL2"]:
                continue
            lipid_molecules.append(idx)
        if not lipid_molecules:
            raise RuntimeError(f"No 14-bead DOPC lipid molecules found in {up_file}")

        nonlipid_index = np.asarray(
            [i for i, cls in enumerate(particle_class) if cls and cls != "LIPID"],
            dtype=np.int64,
        )
        first_frame = int(round(float(frame_start_fraction) * float(max(0, n_frame - 1))))
        frame_candidates = list(range(first_frame, n_frame))
        if not frame_candidates:
            frame_candidates = [n_frame - 1]

        if (
            min_nonlipid_xy_distance_ang > 0.0
            and max_nonlipid_xy_distance_ang > 0.0
            and max_nonlipid_xy_distance_ang < min_nonlipid_xy_distance_ang
        ):
            raise RuntimeError(
                "max_nonlipid_xy_distance_ang must be >= min_nonlipid_xy_distance_ang "
                f"for {up_file}"
            )

        conformers = []
        filtered_near_protein = 0
        filtered_far_from_protein = 0
        filtered_far_from_protein_3d = 0
        for frame in frame_candidates:
            frame_pos = np.asarray(pos[frame, 0, :, :], dtype=np.float64)
            nonlipid_pos = frame_pos[nonlipid_index, :] if nonlipid_index.size else np.zeros((0, 3), dtype=np.float64)
            for lipid in lipid_molecules:
                beads_ang = np.asarray(frame_pos[lipid, :], dtype=np.float64)
                delta = beads_ang - beads_ang[0][None, :]
                if box_ang is not None:
                    delta -= box_ang[None, :] * np.round(delta / box_ang[None, :])
                beads_ang = beads_ang[0][None, :] + delta
                center_ang = np.mean(beads_ang, axis=0)
                delta_nonlipid = None
                if min_nonlipid_xy_distance_ang > 0.0 and nonlipid_pos.size:
                    delta_nonlipid = _minimum_image_ang(nonlipid_pos - center_ang[None, :], box_ang)
                    min_nonlipid_xy = float(
                        np.sqrt(np.min(np.sum(delta_nonlipid[:, :2] * delta_nonlipid[:, :2], axis=1)))
                    )
                    if min_nonlipid_xy < float(min_nonlipid_xy_distance_ang):
                        filtered_near_protein += 1
                        continue
                elif nonlipid_pos.size:
                    delta_nonlipid = _minimum_image_ang(nonlipid_pos - center_ang[None, :], box_ang)
                    min_nonlipid_xy = float(
                        np.sqrt(np.min(np.sum(delta_nonlipid[:, :2] * delta_nonlipid[:, :2], axis=1)))
                    )
                else:
                    min_nonlipid_xy = float("inf")
                    min_nonlipid = float("inf")
                if nonlipid_pos.size:
                    if delta_nonlipid is None:
                        delta_nonlipid = _minimum_image_ang(nonlipid_pos - center_ang[None, :], box_ang)
                    min_nonlipid = float(np.sqrt(np.min(np.sum(delta_nonlipid * delta_nonlipid, axis=1))))
                else:
                    min_nonlipid = float("inf")
                if max_nonlipid_xy_distance_ang > 0.0 and min_nonlipid_xy > float(max_nonlipid_xy_distance_ang):
                    filtered_far_from_protein += 1
                    continue
                if max_nonlipid_distance_ang > 0.0 and min_nonlipid > float(max_nonlipid_distance_ang):
                    filtered_far_from_protein_3d += 1
                    continue
                beads_nm = (beads_ang - center_ang[None, :]) * ANGSTROM_TO_NM
                conformers.append(beads_nm)
    if not conformers:
        raise RuntimeError(
            f"No DOPC conformers remained after filtering {up_file} with "
            f"min_nonlipid_xy_distance_ang={min_nonlipid_xy_distance_ang:.3f}, "
            f"max_nonlipid_xy_distance_ang={max_nonlipid_xy_distance_ang:.3f}, "
            f"max_nonlipid_distance_ang={max_nonlipid_distance_ang:.3f}"
        )
    pool = np.asarray(conformers, dtype=np.float64)
    pool = _canonicalize_lipid_reference_ensemble_to_z(pool)
    kept = int(pool.shape[0])
    total = len(frame_candidates) * len(lipid_molecules)
    if (
        min_nonlipid_xy_distance_ang > 0.0
        or max_nonlipid_xy_distance_ang > 0.0
        or max_nonlipid_distance_ang > 0.0
    ):
        filter_bits = []
        if min_nonlipid_xy_distance_ang > 0.0:
            filter_bits.append(
                f"excluding {filtered_near_protein} near-protein lipids with XY cutoff "
                f"{float(min_nonlipid_xy_distance_ang):.1f} A"
            )
        if max_nonlipid_xy_distance_ang > 0.0:
            filter_bits.append(
                f"excluding {filtered_far_from_protein} far-from-protein lipids with XY cutoff "
                f"{float(max_nonlipid_xy_distance_ang):.1f} A"
            )
        if max_nonlipid_distance_ang > 0.0:
            filter_bits.append(
                f"excluding {filtered_far_from_protein_3d} far-from-protein lipids with 3D cutoff "
                f"{float(max_nonlipid_distance_ang):.1f} A"
            )
        print(
            f"  Loaded DOPC conformer pool from {up_file}: kept {kept} of {total} "
            f"frame-lipid samples after {'; '.join(filter_bits)}"
        )
    else:
        print(
            f"  Loaded DOPC conformer pool from {up_file}: kept {kept} of {total} "
            f"frame-lipid samples"
        )
    return pool


def load_dopc_conformers_from_upside_h5(
    up_file: Path,
    max_conformers: int = 4,
    frame_start_fraction: float = 0.5,
    min_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_distance_ang: float = 0.0,
) -> np.ndarray:
    """Extract centered DOPC conformers from one full-resolution Upside trajectory."""
    return load_dopc_conformers_from_upside_h5_pool(
        [up_file],
        max_conformers=max_conformers,
        frame_start_fraction=frame_start_fraction,
        min_nonlipid_xy_distance_ang=min_nonlipid_xy_distance_ang,
        max_nonlipid_xy_distance_ang=max_nonlipid_xy_distance_ang,
        max_nonlipid_distance_ang=max_nonlipid_distance_ang,
    )


def _balanced_reference_pool_counts(pool_sizes: Iterable[int], total_count: int) -> list[int]:
    pool_sizes = [max(0, int(size)) for size in pool_sizes]
    total = min(max(0, int(total_count)), sum(pool_sizes))
    counts = [0] * len(pool_sizes)
    remaining = total
    while remaining > 0:
        active = [i for i, size in enumerate(pool_sizes) if counts[i] < size]
        if not active:
            break
        share, extra = divmod(remaining, len(active))
        progressed = False
        for rank, idx in enumerate(active):
            capacity = pool_sizes[idx] - counts[idx]
            add = min(capacity, share + (1 if rank < extra else 0))
            if add > 0:
                counts[idx] += add
                remaining -= add
                progressed = True
        if not progressed:
            break
    return counts


def load_dopc_conformers_from_upside_h5_pool(
    up_files: Iterable[Path],
    max_conformers: int = 4,
    frame_start_fraction: float = 0.5,
    min_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_xy_distance_ang: float = 0.0,
    max_nonlipid_distance_ang: float = 0.0,
) -> np.ndarray:
    """Extract centered DOPC conformers from one or more full-resolution Upside trajectories."""
    resolved_paths = [Path(path).expanduser().resolve() for path in up_files]
    if not resolved_paths:
        raise RuntimeError("No Upside trajectory paths were provided for DOPC conformer extraction")
    max_conformers = max(1, int(max_conformers))
    pools = [
        _extract_dopc_conformer_pool_from_upside_h5(
            path,
            frame_start_fraction=frame_start_fraction,
            min_nonlipid_xy_distance_ang=min_nonlipid_xy_distance_ang,
            max_nonlipid_xy_distance_ang=max_nonlipid_xy_distance_ang,
            max_nonlipid_distance_ang=max_nonlipid_distance_ang,
        )
        for path in resolved_paths
    ]
    pooled = np.concatenate(pools, axis=0)
    n_pool = int(pooled.shape[0])
    n_sample = min(max_conformers, n_pool)
    if n_sample < max_conformers:
        print(
            f"  Requested {max_conformers} DOPC conformers but only "
            f"{n_pool} pooled frame-lipid samples are available"
        )
    if len(pools) == 1:
        selection = _select_reference_ensemble_representatives(
            pooled,
            representative_count=n_sample,
        )
        refs = np.asarray(selection["representative_refs_nm"], dtype=np.float64)
        source_count_desc = ""
    else:
        per_source_counts = _balanced_reference_pool_counts(
            [pool.shape[0] for pool in pools],
            n_sample,
        )
        refs_by_source = []
        source_bits = []
        for path, pool, count in zip(resolved_paths, pools, per_source_counts):
            if count <= 0:
                continue
            selection = _select_reference_ensemble_representatives(
                pool,
                representative_count=count,
            )
            refs_by_source.append(np.asarray(selection["representative_refs_nm"], dtype=np.float64))
            source_bits.append(f"{path.stem}:{count}")
        if not refs_by_source:
            raise RuntimeError("No DOPC representatives were selected from the pooled reference trajectories")
        refs = np.concatenate(refs_by_source, axis=0)
        order = np.argsort(_dopc_tail_extension_series_ang(refs))
        refs = refs[order]
        source_count_desc = f" (balanced per-source counts: {', '.join(source_bits)})"
    source_desc = str(resolved_paths[0]) if len(resolved_paths) == 1 else f"{len(resolved_paths)} trajectories"
    compaction = _dopc_tail_extension_series_ang(refs)
    if (
        min_nonlipid_xy_distance_ang > 0.0
        or max_nonlipid_xy_distance_ang > 0.0
        or max_nonlipid_distance_ang > 0.0
    ):
        filter_desc = []
        if min_nonlipid_xy_distance_ang > 0.0:
            filter_desc.append(f"XY >= {float(min_nonlipid_xy_distance_ang):.1f} A")
        if max_nonlipid_xy_distance_ang > 0.0:
            filter_desc.append(f"XY <= {float(max_nonlipid_xy_distance_ang):.1f} A")
        if max_nonlipid_distance_ang > 0.0:
            filter_desc.append(f"3D <= {float(max_nonlipid_distance_ang):.1f} A")
        print(
            f"  Loaded {refs.shape[0]} DOPC conformers from {source_desc} "
            f"using representative pooled samples after protein-distance filtering "
            f"({' and '.join(filter_desc)}){source_count_desc}"
        )
    else:
        print(
            f"  Loaded {refs.shape[0]} DOPC conformers from {source_desc} "
            f"using representative pooled samples{source_count_desc}"
        )
    print(
        f"  Bilayer reference compaction summary: mean={float(np.mean(compaction)):.3f} A, "
        f"p25={float(np.quantile(compaction, 0.25)):.3f} A, "
        f"median={float(np.median(compaction)):.3f} A, "
        f"p75={float(np.quantile(compaction, 0.75)):.3f} A"
    )
    return refs


def _dopc_lipid_molecules_from_upside_h5(h5: h5py.File, up_file: Path) -> list[np.ndarray]:
    atom_names = _decode_h5_string_array(h5["/input/atom_names"][:])
    particle_class = (
        _decode_h5_string_array(h5["/input/particle_class"][:])
        if "/input/particle_class" in h5
        else [""] * len(atom_names)
    )
    molecule_ids = np.asarray(h5["/input/molecule_ids"][:], dtype=np.int64)
    lipid_molecules = []
    for mol_id in sorted(set(int(x) for x in molecule_ids)):
        idx = np.where(molecule_ids == mol_id)[0]
        if idx.size != 14:
            continue
        if any(particle_class[i] and particle_class[i] != "LIPID" for i in idx):
            continue
        names = [atom_names[i] for i in idx]
        if names[:4] != ["NC3", "PO4", "GL1", "GL2"]:
            continue
        lipid_molecules.append(idx)
    if not lipid_molecules:
        raise RuntimeError(f"No 14-bead DOPC lipid molecules found in {up_file}")
    return lipid_molecules


def _extract_cgl_centers_and_orientations_from_frame(
    frame_pos_ang: np.ndarray,
    lipid_molecules: list[np.ndarray],
    box_ang: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    centers = []
    directions = []
    for lipid in lipid_molecules:
        beads_ang = np.asarray(frame_pos_ang[lipid, :], dtype=np.float64)
        delta = beads_ang - beads_ang[0][None, :]
        if box_ang is not None:
            delta -= box_ang[None, :] * np.round(delta / box_ang[None, :])
        beads_ang = beads_ang[0][None, :] + delta
        center = np.mean(beads_ang, axis=0)
        direction = ((beads_ang[8] + beads_ang[13]) * 0.5) - beads_ang[0]
        norm = float(np.linalg.norm(direction))
        if norm <= 1.0e-12:
            continue
        centers.append(center)
        directions.append(direction / norm)
    if not centers:
        raise RuntimeError("No valid DOPC centers/orientations extracted from frame")
    return np.asarray(centers, dtype=np.float64), np.asarray(directions, dtype=np.float64)


def _cgl_particle_indices_from_upside_h5(h5: h5py.File) -> np.ndarray:
    if "/input/type" not in h5:
        return np.zeros(0, dtype=np.int64)
    atom_types = _decode_h5_string_array(h5["/input/type"][:])
    return np.asarray([i for i, name in enumerate(atom_types) if name == "CGL"], dtype=np.int64)


def _extract_dynamic_cgl_centers_and_orientations_from_frame(
    h5: h5py.File,
    frame: int,
    cgl_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    centers = np.asarray(h5["/output/pos"][frame, 0, cgl_indices, :], dtype=np.float64)
    if "/output/cgl_orientation" in h5:
        directions = np.asarray(h5["/output/cgl_orientation"][frame, 0, :, :], dtype=np.float64)
    elif "/input/potential/cgl_orientation_state/direction" in h5:
        directions = np.asarray(h5["/input/potential/cgl_orientation_state/direction"][:], dtype=np.float64)
    else:
        raise RuntimeError("CGL trajectory lacks /output/cgl_orientation and input orientation state")
    if directions.shape != centers.shape:
        raise RuntimeError(
            f"CGL orientation shape {directions.shape} does not match center shape {centers.shape}"
        )
    norm = np.linalg.norm(directions, axis=1)
    valid = norm > 1.0e-12
    if not np.any(valid):
        raise RuntimeError("No valid CGL orientations extracted from frame")
    return centers[valid], directions[valid] / norm[valid, None]


def _grid_edges_from_centers(values: np.ndarray, lo: float | None = None, hi: float | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    edges = np.empty(values.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    if values.size > 1:
        edges[0] = values[0] - 0.5 * float(values[1] - values[0])
        edges[-1] = values[-1] + 0.5 * float(values[-1] - values[-2])
    else:
        edges[0] = values[0] - 0.5
        edges[-1] = values[0] + 0.5
    if lo is not None:
        edges[0] = float(lo)
    if hi is not None:
        edges[-1] = float(hi)
    return edges


def _center_orientation_pair_counts_from_upside_h5(
    up_file: Path,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    frame_start_fraction: float = 0.5,
) -> dict:
    up_file = Path(up_file).expanduser().resolve()
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    if r_values_nm.ndim != 1 or r_values_nm.size < 2:
        raise ValueError("r_values_nm must be a 1D grid with at least two points")
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    edges_nm = _grid_edges_from_centers(r_values_nm, lo=0.0)
    angle_edges = _grid_edges_from_centers(cos_theta_grid, lo=-1.0, hi=1.0)

    counts = np.zeros((r_values_nm.size, cos_theta_grid.size, cos_theta_grid.size), dtype=np.float64)
    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        cgl_indices = _cgl_particle_indices_from_upside_h5(h5)
        if cgl_indices.size:
            source_kind = "dynamic_vector_cgl"
            lipid_molecules = []
        else:
            source_kind = "full_dopc"
            lipid_molecules = _dopc_lipid_molecules_from_upside_h5(h5, up_file)
        pos = h5["/output/pos"]
        n_frame = int(pos.shape[0])
        first_frame = int(round(float(frame_start_fraction) * float(max(0, n_frame - 1))))
        frame_indices = list(range(first_frame, n_frame))
        if not frame_indices:
            frame_indices = [n_frame - 1]
        box_ang = _read_output_box_lengths_ang(h5)
        for frame in frame_indices:
            if source_kind == "dynamic_vector_cgl":
                centers_ang, directions = _extract_dynamic_cgl_centers_and_orientations_from_frame(
                    h5,
                    frame,
                    cgl_indices,
                )
            else:
                centers_ang, directions = _extract_cgl_centers_and_orientations_from_frame(
                    np.asarray(pos[frame, 0, :, :], dtype=np.float64),
                    lipid_molecules,
                    box_ang,
                )
            for i in range(centers_ang.shape[0] - 1):
                delta = centers_ang[i + 1 :, :] - centers_ang[i][None, :]
                if box_ang is not None:
                    delta -= box_ang[None, :] * np.round(delta / box_ang[None, :])
                r_ang = np.linalg.norm(delta, axis=1)
                valid = r_ang > 1.0e-12
                if not np.any(valid):
                    continue
                delta = delta[valid]
                r_ang = r_ang[valid]
                n12 = delta / r_ang[:, None]
                dir1 = directions[i][None, :]
                dir2 = directions[i + 1 :, :][valid]
                a1 = -np.sum(dir1 * n12, axis=1)
                a2 = np.sum(dir2 * n12, axis=1)
                r_nm = r_ang * ANGSTROM_TO_NM
                samples = np.column_stack([r_nm, a1, a2])
                reverse_samples = np.column_stack([r_nm, a2, a1])
                hist, _ = np.histogramdd(
                    np.vstack([samples, reverse_samples]),
                    bins=(edges_nm, angle_edges, angle_edges),
                )
                counts += hist.astype(np.float64)
    print(
        f"  Loaded center/orientation pair counts from {up_file}: "
        f"{int(np.sum(counts))} ordered pair samples, source={source_kind}, "
        f"{int(np.count_nonzero(counts > 0.0))}/{counts.size} occupied tensor bins"
    )
    return {
        "counts": counts,
        "occupied": counts > 0.0,
        "source": source_kind,
        "frame_start_fraction": float(frame_start_fraction),
        "pair_sample_count": int(np.sum(counts)),
        "occupied_radial_bins": int(np.count_nonzero(np.any(counts > 0.0, axis=(1, 2)))),
        "occupied_tensor_bins": int(np.count_nonzero(counts > 0.0)),
    }


def _fluid_bilayer_pair_pmf_from_upside_h5(
    up_file: Path,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    temperature_upside: float,
    frame_start_fraction: float = 0.5,
) -> dict:
    dist = _center_orientation_pair_counts_from_upside_h5(
        up_file,
        r_values_nm,
        cos_theta_grid,
        frame_start_fraction=frame_start_fraction,
    )
    counts = np.asarray(dist["counts"], dtype=np.float64)
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    edges_nm = _grid_edges_from_centers(r_values_nm, lo=0.0)

    shell_volume = (4.0 * math.pi / 3.0) * (edges_nm[1:] ** 3 - edges_nm[:-1] ** 3)
    density = counts / np.maximum(shell_volume[:, None, None], 1.0e-12)
    finite = density > 0.0
    pmf = np.full_like(counts, np.nan, dtype=np.float64)
    if np.any(finite):
        kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
        pmf[finite] = -kbt * np.log(density[finite])
        pmf[finite] -= float(np.min(pmf[finite]))
    print(
        f"  Loaded fluid-bilayer orientation-resolved PMF from {up_file}: "
        f"{int(np.sum(counts))} ordered lipid-pair samples, "
        f"{int(np.count_nonzero(finite))}/{counts.size} occupied tensor bins"
    )
    return {
        "pmf_kj_mol": pmf,
        "counts": counts,
        "occupied": finite,
        "source": "full_dopc_fluid_bilayer_oriented_pair_pmf",
        "frame_start_fraction": float(dist["frame_start_fraction"]),
        "pair_sample_count": int(dist["pair_sample_count"]),
        "occupied_radial_bins": int(np.count_nonzero(np.any(finite, axis=(1, 2)))),
        "occupied_tensor_bins": int(np.count_nonzero(finite)),
    }


def _relative_entropy_pair_ibi_correction_from_upside_h5(
    target_up_file: Path,
    model_up_file: Path,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    temperature_upside: float,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    step_size: float = 1.0,
) -> dict:
    step_size = float(step_size)
    if not (0.0 < step_size <= 1.0):
        raise ValueError(f"IBI step_size must be in (0, 1], got {step_size}")
    target = _center_orientation_pair_counts_from_upside_h5(
        target_up_file,
        r_values_nm,
        cos_theta_grid,
        frame_start_fraction=target_frame_start_fraction,
    )
    model = _center_orientation_pair_counts_from_upside_h5(
        model_up_file,
        r_values_nm,
        cos_theta_grid,
        frame_start_fraction=model_frame_start_fraction,
    )
    target_counts = np.asarray(target["counts"], dtype=np.float64)
    model_counts = np.asarray(model["counts"], dtype=np.float64)
    target_total = max(float(np.sum(target_counts)), 1.0)
    model_total = max(float(np.sum(model_counts)), 1.0)
    target_probability = target_counts / target_total
    model_probability = model_counts / model_total
    sampled = (target_counts > 0.0) & (model_counts > 0.0)
    raw_correction = np.full_like(target_counts, np.nan, dtype=np.float64)
    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    raw_correction[sampled] = kbt * np.log(model_probability[sampled] / target_probability[sampled])
    correction = raw_correction * step_size
    raw_finite = raw_correction[np.isfinite(raw_correction)]
    finite = correction[np.isfinite(correction)]
    if raw_finite.size:
        raw_corr_min = float(np.min(raw_finite))
        raw_corr_max = float(np.max(raw_finite))
        raw_corr_mean = float(np.mean(raw_finite))
    else:
        raw_corr_min = raw_corr_max = raw_corr_mean = 0.0
    if finite.size:
        corr_min = float(np.min(finite))
        corr_max = float(np.max(finite))
        corr_mean = float(np.mean(finite))
    else:
        corr_min = corr_max = corr_mean = 0.0
    print(
        f"  IBI CGL-CGL correction from target={target_up_file} and model={model_up_file}: "
        f"{int(np.count_nonzero(sampled))}/{sampled.size} sampled tensor bins, "
        f"raw deltaU range {raw_corr_min:.3f}..{raw_corr_max:.3f} kJ/mol, "
        f"applied step={step_size:.3f} range {corr_min:.3f}..{corr_max:.3f} kJ/mol"
    )
    return {
        "correction_kj_mol": correction,
        "raw_correction_kj_mol": raw_correction,
        "sampled": sampled,
        "target_counts": target_counts,
        "model_counts": model_counts,
        "target_pair_sample_count": int(target["pair_sample_count"]),
        "model_pair_sample_count": int(model["pair_sample_count"]),
        "target_occupied_tensor_bins": int(target["occupied_tensor_bins"]),
        "model_occupied_tensor_bins": int(model["occupied_tensor_bins"]),
        "sampled_tensor_bins": int(np.count_nonzero(sampled)),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "step_size": step_size,
        "raw_correction_min_kj_mol": raw_corr_min,
        "raw_correction_max_kj_mol": raw_corr_max,
        "raw_correction_mean_kj_mol": raw_corr_mean,
        "correction_min_kj_mol": corr_min,
        "correction_max_kj_mol": corr_max,
        "correction_mean_kj_mol": corr_mean,
    }


def _cgl_density_frames_from_upside_h5(
    up_file: Path,
    frame_start_fraction: float = 0.5,
    max_frames: int = 400,
) -> list[tuple[np.ndarray, np.ndarray | None]]:
    up_file = Path(up_file).expanduser().resolve()
    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        pos = h5["/output/pos"]
        frames = _sampled_frame_indices(
            int(pos.shape[0]),
            frame_start_fraction=frame_start_fraction,
            max_frames=max_frames,
        )

        cgl_indices = _cgl_particle_indices_from_upside_h5(h5)
        if cgl_indices.size:
            source_kind = "dynamic_vector_cgl"
            lipid_molecules = []
        else:
            source_kind = "full_dopc"
            lipid_molecules = _dopc_lipid_molecules_from_upside_h5(h5, up_file)
        box_ang = _read_output_box_lengths_ang(h5)
        out = []
        for frame in frames:
            if source_kind == "dynamic_vector_cgl":
                centers_ang, _ = _extract_dynamic_cgl_centers_and_orientations_from_frame(
                    h5,
                    frame,
                    cgl_indices,
                )
            else:
                centers_ang, _ = _extract_cgl_centers_and_orientations_from_frame(
                    np.asarray(pos[frame, 0, :, :], dtype=np.float64),
                    lipid_molecules,
                    box_ang,
                )
            out.append((centers_ang.astype(np.float64), None if box_ang is None else box_ang.astype(np.float64)))
    print(f"  Loaded {len(out)} local-density frames from {up_file}, source={source_kind}")
    return out


def _minimum_image_ang(delta: np.ndarray, box_ang: np.ndarray | None) -> np.ndarray:
    out = np.asarray(delta, dtype=np.float64).copy()
    if box_ang is None:
        return out
    box = np.asarray(box_ang, dtype=np.float64)
    valid = box > 0.0
    out[..., valid] -= box[valid] * np.round(out[..., valid] / box[valid])
    return out


def _derive_density_kernel_cutoff_ang(
    target_frames: list[tuple[np.ndarray, np.ndarray | None]],
    hist_max_ang: float = 24.0,
) -> tuple[float, dict]:
    distances = []
    for centers_ang, box_ang in target_frames:
        if centers_ang.shape[0] < 4:
            continue
        mid_z = float(np.median(centers_ang[:, 2]))
        for mask in (centers_ang[:, 2] <= mid_z, centers_ang[:, 2] > mid_z):
            xy = centers_ang[mask, :2]
            if xy.shape[0] < 2:
                continue
            box_xy = None if box_ang is None else box_ang[:2]
            for i in range(xy.shape[0] - 1):
                delta = xy[i + 1 :, :] - xy[i][None, :]
                if box_xy is not None:
                    delta -= box_xy[None, :] * np.round(delta / box_xy[None, :])
                r = np.sqrt(np.sum(delta * delta, axis=1))
                distances.append(r)
    if not distances:
        raise RuntimeError("Cannot derive density kernel cutoff: no same-leaflet distances")
    distances = np.concatenate(distances)
    distances = distances[(distances > 1.0e-6) & (distances < float(hist_max_ang))]
    if distances.size < 100:
        raise RuntimeError("Cannot derive density kernel cutoff: too few same-leaflet distances")
    edges = np.linspace(0.0, float(hist_max_ang), 97)
    counts, edges = np.histogram(distances, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak_window = (centers >= 4.0) & (centers <= 12.0)
    if not np.any(peak_window):
        raise RuntimeError("Cannot derive density kernel cutoff: no first-shell peak window")
    peak_candidates = np.where(peak_window)[0]
    peak_idx = int(peak_candidates[int(np.argmax(counts[peak_candidates]))])
    search = np.arange(peak_idx + 1, len(centers))
    search = search[centers[search] <= 18.0]
    min_idx = None
    for idx in search:
        if idx <= 0 or idx >= len(counts) - 1:
            continue
        if counts[idx] <= counts[idx - 1] and counts[idx] <= counts[idx + 1]:
            min_idx = int(idx)
            break
    if min_idx is None and search.size:
        min_idx = int(search[int(np.argmin(counts[search]))])
    if min_idx is None:
        min_idx = int(np.argmin(np.abs(centers - min(18.0, 1.45 * centers[peak_idx]))))
    cutoff = float(centers[min_idx])
    cutoff = max(8.0, min(18.0, cutoff))
    return cutoff, {
        "rdf_hist_centers_ang": centers,
        "rdf_hist_counts": counts.astype(np.float64),
        "rdf_first_peak_ang": float(centers[peak_idx]),
        "rdf_first_minimum_ang": cutoff,
        "rdf_distance_sample_count": int(distances.size),
    }


def _density_kernel_values(r_ang: np.ndarray, cutoff_ang: float) -> np.ndarray:
    x = np.asarray(r_ang, dtype=np.float64) / float(cutoff_ang)
    out = np.zeros_like(x, dtype=np.float64)
    inside = x < 1.0
    y = 1.0 - x[inside] * x[inside]
    out[inside] = y * y
    return out


def _local_density_samples(
    frames: list[tuple[np.ndarray, np.ndarray | None]],
    cutoff_ang: float,
) -> np.ndarray:
    samples = []
    for centers_ang, box_ang in frames:
        n = centers_ang.shape[0]
        rho = np.zeros(n, dtype=np.float64)
        mid_z = float(np.median(centers_ang[:, 2]))
        leaflet_masks = (
            centers_ang[:, 2] <= mid_z,
            centers_ang[:, 2] > mid_z,
        )
        for mask in leaflet_masks:
            members = np.where(mask)[0]
            if members.size < 2:
                continue
            leaflet_centers = centers_ang[members, :]
            for local_i in range(members.size - 1):
                delta = _minimum_image_ang(
                    leaflet_centers[local_i + 1 :, :] - leaflet_centers[local_i][None, :],
                    box_ang,
                )
                r = np.sqrt(np.sum(delta * delta, axis=1))
                valid = (r > 1.0e-12) & (r < float(cutoff_ang))
                if not np.any(valid):
                    continue
                w = _density_kernel_values(r[valid], cutoff_ang)
                global_i = int(members[local_i])
                global_j = members[local_i + 1 :][valid]
                rho[global_i] += float(np.sum(w))
                np.add.at(rho, global_j, w)
        samples.append(rho)
    if not samples:
        raise RuntimeError("No local-density samples were computed")
    return np.concatenate(samples)


def _build_relative_entropy_density_table(
    target_up_file: Path,
    model_up_file: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    n_kernel: int = 24,
    n_embedding: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    target_frames = _cgl_density_frames_from_upside_h5(
        target_up_file,
        frame_start_fraction=target_frame_start_fraction,
        max_frames=max_frames,
    )
    model_frames = _cgl_density_frames_from_upside_h5(
        model_up_file,
        frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
    )
    cutoff_ang, cutoff_meta = _derive_density_kernel_cutoff_ang(target_frames)

    target_rho = _local_density_samples(target_frames, cutoff_ang)
    model_rho = _local_density_samples(model_frames, cutoff_ang)
    rho_lo = float(min(np.min(target_rho), np.min(model_rho)))
    rho_hi = float(max(np.max(target_rho), np.max(model_rho)))
    pad = max(0.05 * (rho_hi - rho_lo), 0.05)
    rho_min = max(0.0, rho_lo - pad)
    rho_max = rho_hi + pad
    hist_edges = np.linspace(rho_min, rho_max, int(n_hist) + 1, dtype=np.float64)
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    target_counts, _ = np.histogram(target_rho, bins=hist_edges)
    model_counts, _ = np.histogram(model_rho, bins=hist_edges)
    target_prob = (target_counts.astype(np.float64) + float(pseudocount))
    model_prob = (model_counts.astype(np.float64) + float(pseudocount))
    target_prob /= float(np.sum(target_prob))
    model_prob /= float(np.sum(model_prob))
    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    correction_kj_mol = kbt * np.log(model_prob / target_prob)
    correction_kj_mol -= float(np.average(correction_kj_mol, weights=target_prob))

    kernel_spacing = float(cutoff_ang) / float(int(n_kernel) - 2)
    kernel_knot_vector = np.zeros(int(n_kernel) + 4, dtype=np.float64)
    kernel_knot_vector[4:-4] = np.arange(1, int(n_kernel) - 3, dtype=np.float64)
    kernel_knot_vector[-4:] = kernel_knot_vector[-5]
    kernel_r = np.linspace(0.0, float(cutoff_ang), 160, dtype=np.float64)
    kernel_coeff = _fit_radial_bspline(
        kernel_r / kernel_spacing,
        _density_kernel_values(kernel_r, cutoff_ang),
        kernel_knot_vector,
        smooth=0.0,
    )

    rho_spacing = (rho_max - rho_min) / float(int(n_embedding) - 3)
    emb_knot_vector = np.zeros(int(n_embedding) + 4, dtype=np.float64)
    emb_knot_vector[4:-4] = np.arange(1, int(n_embedding) - 3, dtype=np.float64)
    emb_knot_vector[-4:] = emb_knot_vector[-5]
    emb_t = 1.0 + (hist_centers - rho_min) / rho_spacing
    embedding_coeff = _fit_radial_bspline(
        emb_t,
        correction_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        emb_knot_vector,
        smooth=float(smooth),
    )

    finite = correction_kj_mol[np.isfinite(correction_kj_mol)]
    print(
        f"  CGL density correction: cutoff={cutoff_ang:.3f} A, "
        f"rho={rho_min:.3f}..{rho_max:.3f}, "
        f"deltaU={float(np.min(finite)):.3f}..{float(np.max(finite)):.3f} kJ/mol"
    )
    return {
        "kernel_coeff": kernel_coeff.astype(np.float32),
        "embedding_coeff": embedding_coeff.astype(np.float32),
        "kernel_cutoff_ang": float(cutoff_ang),
        "kernel_knot_spacing_ang": float(kernel_spacing),
        "kernel_n_knot": int(n_kernel),
        "embedding_rho_min": float(rho_min),
        "embedding_rho_max": float(rho_max),
        "embedding_rho_spacing": float(rho_spacing),
        "embedding_n_knot": int(n_embedding),
        "target_rho": target_rho.astype(np.float32),
        "model_rho": model_rho.astype(np.float32),
        "rho_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "model_counts": model_counts.astype(np.float32),
        "correction_kj_mol": correction_kj_mol.astype(np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "target_frame_count": int(len(target_frames)),
        "model_frame_count": int(len(model_frames)),
        "target_sample_count": int(target_rho.size),
        "model_sample_count": int(model_rho.size),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "target_upside_h5_path": str(Path(target_up_file).expanduser().resolve()),
        "model_upside_h5_path": str(Path(model_up_file).expanduser().resolve()),
        **cutoff_meta,
    }


def add_cg_lipid_density_table_to_h5(
    base_h5_path: Path,
    output_path: Path,
    target_upside_h5_path: Path,
    model_upside_h5_path: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
) -> None:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    density = _build_relative_entropy_density_table(
        target_upside_h5_path,
        model_upside_h5_path,
        temperature_upside=temperature_upside,
        target_frame_start_fraction=target_frame_start_fraction,
        model_frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
    )

    def _writer(h5: h5py.File) -> None:
        with h5py.File(base_h5_path, "r") as src:
            for key, value in src.attrs.items():
                h5.attrs[key] = value
            for key in src.keys():
                src.copy(key, h5)
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_density" in cg_grp:
            del cg_grp["cg_lipid_density"]
        dens = cg_grp.create_group("cg_lipid_density")
        dens.attrs["schema"] = "cg_lipid_density_relative_entropy_v1"
        dens.attrs["source"] = "relative_entropy_local_density_target_full_dopc_model_one_particle_cgl"
        dens.attrs["runtime_representation"] = "sum_i_F_of_local_density_spline_kernel"
        dens.attrs["correction_layer"] = "conservative_many_body_local_density"
        dens.attrs["one_particle_cgl_preserved"] = np.int32(1)
        dens.attrs["orientation_potential"] = np.int32(0)
        dens.attrs["force_applies_to"] = "cgl_centers_only"
        dens.attrs["kernel_source"] = "compact_quadratic_kernel_with_cutoff_from_full_dopc_same_leaflet_center_rdf_first_minimum"
        dens.attrs["embedding_update_rule"] = "F_density=kBT*ln(P_model(rho)/P_target(rho))"
        dens.attrs["boltzmann_temperature_upside"] = np.float32(temperature_upside)
        dens.attrs["energy_conversion_kj_per_eup"] = np.float32(ENERGY_CONVERSION_KJ_PER_EUP)
        dens.attrs["length_conversion_ang_per_nm"] = np.float32(LENGTH_CONVERSION_A_PER_NM)
        for attr_name in (
            "kernel_cutoff_ang",
            "kernel_knot_spacing_ang",
            "kernel_n_knot",
            "embedding_rho_min",
            "embedding_rho_max",
            "embedding_rho_spacing",
            "embedding_n_knot",
            "pseudocount",
            "smooth",
            "target_frame_count",
            "model_frame_count",
            "target_sample_count",
            "model_sample_count",
            "target_frame_start_fraction",
            "model_frame_start_fraction",
            "target_upside_h5_path",
            "model_upside_h5_path",
            "rdf_first_peak_ang",
            "rdf_first_minimum_ang",
            "rdf_distance_sample_count",
        ):
            value = density[attr_name]
            if isinstance(value, (int, np.integer)):
                dens.attrs[attr_name] = np.int32(value)
            elif isinstance(value, (float, np.floating)):
                dens.attrs[attr_name] = np.float32(value)
            else:
                dens.attrs[attr_name] = value
        dens.create_dataset("kernel_coeff", data=density["kernel_coeff"])
        dens.create_dataset("embedding_coeff", data=density["embedding_coeff"])
        dens.create_dataset("target_rho", data=density["target_rho"])
        dens.create_dataset("model_rho", data=density["model_rho"])
        dens.create_dataset("rho_hist_centers", data=density["rho_hist_centers"])
        dens.create_dataset("target_counts", data=density["target_counts"])
        dens.create_dataset("model_counts", data=density["model_counts"])
        dens.create_dataset("correction_kj_mol", data=density["correction_kj_mol"])
        dens.create_dataset("rdf_hist_centers_ang", data=density["rdf_hist_centers_ang"].astype(np.float32))
        dens.create_dataset("rdf_hist_counts", data=density["rdf_hist_counts"].astype(np.float32))

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} with cg_lipid_density")


def _cgl_contact_frames_from_upside_h5(
    up_file: Path,
    frame_start_fraction: float = 0.5,
    max_frames: int = 400,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray | None]]:
    up_file = Path(up_file).expanduser().resolve()
    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        pos = h5["/output/pos"]
        n_frame = int(pos.shape[0])
        first_frame = int(round(float(frame_start_fraction) * float(max(0, n_frame - 1))))
        frames = list(range(first_frame, n_frame))
        if not frames:
            frames = [n_frame - 1]
        max_frames = max(1, int(max_frames))
        if len(frames) > max_frames:
            pick = np.linspace(0, len(frames) - 1, max_frames, dtype=int)
            frames = [frames[int(i)] for i in pick]

        cgl_indices = _cgl_particle_indices_from_upside_h5(h5)
        if cgl_indices.size:
            source_kind = "dynamic_vector_cgl"
            lipid_molecules = []
        else:
            source_kind = "full_dopc"
            lipid_molecules = _dopc_lipid_molecules_from_upside_h5(h5, up_file)
        box_ang = _read_output_box_lengths_ang(h5)
        out = []
        for frame in frames:
            if source_kind == "dynamic_vector_cgl":
                centers_ang, directions = _extract_dynamic_cgl_centers_and_orientations_from_frame(
                    h5,
                    frame,
                    cgl_indices,
                )
            else:
                centers_ang, directions = _extract_cgl_centers_and_orientations_from_frame(
                    np.asarray(pos[frame, 0, :, :], dtype=np.float64),
                    lipid_molecules,
                    box_ang,
                )
            out.append(
                (
                    centers_ang.astype(np.float64),
                    directions.astype(np.float64),
                    None if box_ang is None else box_ang.astype(np.float64),
                )
            )
    print(f"  Loaded {len(out)} cross-leaflet contact frames from {up_file}, source={source_kind}")
    return out


def _sampled_frame_indices(
    n_frame: int,
    frame_start_fraction: float,
    max_frames: int,
) -> list[int]:
    first_frame = int(round(float(frame_start_fraction) * float(max(0, n_frame - 1))))
    frames = list(range(first_frame, n_frame))
    if not frames:
        frames = [n_frame - 1]
    max_frames = max(1, int(max_frames))
    if len(frames) > max_frames:
        pick = np.linspace(0, len(frames) - 1, max_frames, dtype=int)
        frames = [frames[int(i)] for i in pick]
    return frames


def _pair_cross_leaflet_face_weight(
    delta_ang: np.ndarray,
    dir1: np.ndarray,
    dir2: np.ndarray,
    radial_cutoff_ang: float,
    face_cos_min: float,
    smooth_weight: bool,
) -> np.ndarray:
    delta = np.asarray(delta_ang, dtype=np.float64)
    r = np.sqrt(np.sum(delta * delta, axis=1))
    valid = (r > 1.0e-12) & (r <= float(radial_cutoff_ang))
    out = np.zeros(delta.shape[0], dtype=np.float64)
    if not np.any(valid):
        return out
    unit = delta[valid] / r[valid, None]
    dir2_valid = np.asarray(dir2, dtype=np.float64)[valid]
    face1 = np.sum(unit * np.asarray(dir1, dtype=np.float64)[None, :], axis=1)
    face2 = -np.sum(unit * dir2_valid, axis=1)
    if not smooth_weight:
        out[valid] = ((face1 >= float(face_cos_min)) & (face2 >= float(face_cos_min))).astype(np.float64)
        return out
    ang_denom = max(1.0 - float(face_cos_min), 1.0e-6)
    w1 = np.clip((face1 - float(face_cos_min)) / ang_denom, 0.0, 1.0)
    w2 = np.clip((face2 - float(face_cos_min)) / ang_denom, 0.0, 1.0)
    wr = np.clip(1.0 - r[valid] / max(float(radial_cutoff_ang), 1.0e-6), 0.0, 1.0)
    out[valid] = wr * w1 * w2
    return out


def _cross_leaflet_contact_samples(
    frames: list[tuple[np.ndarray, np.ndarray, np.ndarray | None]],
    radial_cutoff_ang: float,
    face_cos_min: float,
    smooth_weight: bool,
) -> np.ndarray:
    samples = []
    for centers_ang, directions, box_ang in frames:
        n = int(centers_ang.shape[0])
        if n < 2:
            continue
        survival = np.ones(n, dtype=np.float64)
        leaflet_side = centers_ang[:, 2] >= float(np.median(centers_ang[:, 2]))
        for ai in range(n - 1):
            delta = _minimum_image_ang(
                centers_ang[ai + 1 :, :] - centers_ang[ai][None, :],
                box_ang,
            )
            cross = leaflet_side[ai + 1 :] != leaflet_side[ai]
            if not np.any(cross):
                continue
            local_idx = np.where(cross)[0]
            weight = _pair_cross_leaflet_face_weight(
                delta[local_idx],
                directions[ai],
                directions[ai + 1 :][local_idx],
                radial_cutoff_ang=radial_cutoff_ang,
                face_cos_min=face_cos_min,
                smooth_weight=smooth_weight,
            )
            if not np.any(weight > 0.0):
                continue
            surv = np.maximum(0.0, 1.0 - weight)
            survival[ai] *= float(np.prod(surv))
            survival[ai + 1 + local_idx] *= surv
        samples.append(1.0 - survival)
    if not samples:
        raise RuntimeError("No cross-leaflet contact samples were computed")
    return np.concatenate(samples)


def _build_relative_entropy_contact_embedding_table(
    target_up_file: Path,
    model_up_file: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = False,
    n_embedding: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    target_frames = _cgl_contact_frames_from_upside_h5(
        target_up_file,
        frame_start_fraction=target_frame_start_fraction,
        max_frames=max_frames,
    )
    model_frames = _cgl_contact_frames_from_upside_h5(
        model_up_file,
        frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
    )
    target_contact = _cross_leaflet_contact_samples(
        target_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
    )
    model_contact = _cross_leaflet_contact_samples(
        model_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
    )
    coord_lo = float(min(np.min(target_contact), np.min(model_contact)))
    coord_hi = float(max(np.max(target_contact), np.max(model_contact)))
    pad = max(0.05 * (coord_hi - coord_lo), 0.02)
    coord_min = max(0.0, coord_lo - pad)
    coord_max = min(1.0, coord_hi + pad)
    hist_edges = np.linspace(coord_min, coord_max, int(n_hist) + 1, dtype=np.float64)
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    target_counts, _ = np.histogram(target_contact, bins=hist_edges)
    model_counts, _ = np.histogram(model_contact, bins=hist_edges)
    target_prob = target_counts.astype(np.float64) + float(pseudocount)
    model_prob = model_counts.astype(np.float64) + float(pseudocount)
    target_prob /= float(np.sum(target_prob))
    model_prob /= float(np.sum(model_prob))
    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    correction_kj_mol = kbt * np.log(model_prob / target_prob)
    correction_kj_mol -= float(np.average(correction_kj_mol, weights=target_prob))

    coord_spacing = max(coord_max - coord_min, 1.0e-3) / float(int(n_embedding) - 3)
    emb_knot_vector = np.zeros(int(n_embedding) + 4, dtype=np.float64)
    emb_knot_vector[4:-4] = np.arange(1, int(n_embedding) - 3, dtype=np.float64)
    emb_knot_vector[-4:] = emb_knot_vector[-5]
    emb_t = 1.0 + (hist_centers - coord_min) / coord_spacing
    embedding_coeff = _fit_radial_bspline(
        emb_t,
        correction_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        emb_knot_vector,
        smooth=float(smooth),
    )

    finite = correction_kj_mol[np.isfinite(correction_kj_mol)]
    print(
        f"  CGL contact embedding: cutoff={radial_cutoff_ang:.3f} A, "
        f"face_cos_min={face_cos_min:.3f}, coord={coord_min:.3f}..{coord_max:.3f}, "
        f"deltaU={float(np.min(finite)):.3f}..{float(np.max(finite)):.3f} kJ/mol"
    )
    return {
        "embedding_coeff": embedding_coeff.astype(np.float32),
        "embedding_coord_min": float(coord_min),
        "embedding_coord_max": float(coord_max),
        "embedding_coord_spacing": float(coord_spacing),
        "embedding_n_knot": int(n_embedding),
        "contact_radial_cutoff_ang": float(radial_cutoff_ang),
        "contact_face_cos_min": float(face_cos_min),
        "contact_smooth_weight": int(bool(smooth_weight)),
        "target_contact": target_contact.astype(np.float32),
        "model_contact": model_contact.astype(np.float32),
        "contact_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "model_counts": model_counts.astype(np.float32),
        "correction_kj_mol": correction_kj_mol.astype(np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "target_frame_count": int(len(target_frames)),
        "model_frame_count": int(len(model_frames)),
        "target_sample_count": int(target_contact.size),
        "model_sample_count": int(model_contact.size),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "target_upside_h5_path": str(Path(target_up_file).expanduser().resolve()),
        "model_upside_h5_path": str(Path(model_up_file).expanduser().resolve()),
    }


def add_cg_lipid_contact_embedding_table_to_h5(
    base_h5_path: Path,
    output_path: Path,
    target_upside_h5_path: Path,
    model_upside_h5_path: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = False,
) -> None:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    embedding = _build_relative_entropy_contact_embedding_table(
        target_upside_h5_path,
        model_upside_h5_path,
        temperature_upside=temperature_upside,
        target_frame_start_fraction=target_frame_start_fraction,
        model_frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
    )

    def _writer(h5: h5py.File) -> None:
        with h5py.File(base_h5_path, "r") as src:
            for key, value in src.attrs.items():
                h5.attrs[key] = value
            for key in src.keys():
                src.copy(key, h5)
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_contact_embedding" in cg_grp:
            del cg_grp["cg_lipid_contact_embedding"]
        embed = cg_grp.create_group("cg_lipid_contact_embedding")
        embed.attrs["schema"] = "cg_lipid_contact_embedding_relative_entropy_v1"
        embed.attrs["source"] = "relative_entropy_cross_leaflet_contact_target_full_dopc_model_one_particle_cgl"
        embed.attrs["runtime_representation"] = "sum_i_F_of_cross_leaflet_contact_field"
        embed.attrs["correction_layer"] = "conservative_many_body_cross_leaflet_contact_embedding"
        embed.attrs["one_particle_cgl_preserved"] = np.int32(1)
        embed.attrs["boltzmann_temperature_upside"] = np.float32(temperature_upside)
        embed.attrs["energy_conversion_kj_per_eup"] = np.float32(ENERGY_CONVERSION_KJ_PER_EUP)
        embed.attrs["length_conversion_ang_per_nm"] = np.float32(LENGTH_CONVERSION_A_PER_NM)
        for attr_name in (
            "embedding_coord_min",
            "embedding_coord_max",
            "embedding_coord_spacing",
            "embedding_n_knot",
            "contact_radial_cutoff_ang",
            "contact_face_cos_min",
            "contact_smooth_weight",
            "pseudocount",
            "smooth",
            "target_frame_count",
            "model_frame_count",
            "target_sample_count",
            "model_sample_count",
            "target_frame_start_fraction",
            "model_frame_start_fraction",
            "target_upside_h5_path",
            "model_upside_h5_path",
        ):
            value = embedding[attr_name]
            if isinstance(value, (int, np.integer)):
                embed.attrs[attr_name] = np.int32(value)
            elif isinstance(value, (float, np.floating)):
                embed.attrs[attr_name] = np.float32(value)
            else:
                embed.attrs[attr_name] = value
        embed.create_dataset("embedding_coeff", data=embedding["embedding_coeff"])
        embed.create_dataset("target_contact", data=embedding["target_contact"])
        embed.create_dataset("model_contact", data=embedding["model_contact"])
        embed.create_dataset("contact_hist_centers", data=embedding["contact_hist_centers"])
        embed.create_dataset("target_counts", data=embedding["target_counts"])
        embed.create_dataset("model_counts", data=embedding["model_counts"])
        embed.create_dataset("correction_kj_mol", data=embedding["correction_kj_mol"])

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} with cg_lipid_contact_embedding")


def _cross_leaflet_gap_distance_for_frame(
    centers_ang: np.ndarray,
    directions: np.ndarray,
    box_ang: np.ndarray | None,
    radial_cutoff_ang: float,
    face_cos_min: float,
    smooth_weight: bool,
    fallback_gap_ang: float | None = None,
) -> np.ndarray:
    if fallback_gap_ang is None:
        fallback_gap_ang = float(radial_cutoff_ang)
    centers_ang = np.asarray(centers_ang, dtype=np.float64)
    directions = np.asarray(directions, dtype=np.float64)
    n = int(centers_ang.shape[0])
    if n < 1:
        return np.zeros(0, dtype=np.float64)
    numerator = np.zeros(n, dtype=np.float64)
    denominator = np.zeros(n, dtype=np.float64)
    leaflet_side = centers_ang[:, 2] >= float(np.median(centers_ang[:, 2]))
    for ai in range(n - 1):
        delta = _minimum_image_ang(
            centers_ang[ai + 1 :, :] - centers_ang[ai][None, :],
            box_ang,
        )
        cross = leaflet_side[ai + 1 :] != leaflet_side[ai]
        if not np.any(cross):
            continue
        local_idx = np.where(cross)[0]
        local_delta = delta[local_idx]
        weight = _pair_cross_leaflet_face_weight(
            local_delta,
            directions[ai],
            directions[ai + 1 :][local_idx],
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
        )
        positive = weight > 0.0
        if not np.any(positive):
            continue
        local_idx = local_idx[positive]
        local_delta = local_delta[positive]
        weight = weight[positive]
        distance = np.sqrt(np.sum(local_delta * local_delta, axis=1))
        numerator[ai] += float(np.sum(weight * distance))
        denominator[ai] += float(np.sum(weight))
        global_j = ai + 1 + local_idx
        np.add.at(numerator, global_j, weight * distance)
        np.add.at(denominator, global_j, weight)
    gap = np.full(n, float(fallback_gap_ang), dtype=np.float64)
    valid = denominator > 1.0e-12
    gap[valid] = numerator[valid] / denominator[valid]
    return gap


def _cross_leaflet_gap_distance_samples(
    frames: list[tuple[np.ndarray, np.ndarray, np.ndarray | None]],
    radial_cutoff_ang: float,
    face_cos_min: float,
    smooth_weight: bool,
    fallback_gap_ang: float | None = None,
) -> np.ndarray:
    samples = []
    for centers_ang, directions, box_ang in frames:
        gap = _cross_leaflet_gap_distance_for_frame(
            centers_ang,
            directions,
            box_ang,
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
            fallback_gap_ang=fallback_gap_ang,
        )
        if gap.size < 1:
            continue
        samples.append(gap)
    if not samples:
        raise RuntimeError("No cross-leaflet gap-distance samples were computed")
    return np.concatenate(samples)


def _full_dopc_gap_and_compaction_samples_from_upside_h5(
    up_file: Path,
    frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
) -> dict:
    up_file = Path(up_file).expanduser().resolve()
    if fallback_gap_ang is None:
        fallback_gap_ang = float(radial_cutoff_ang)
    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        pos = h5["/output/pos"]
        frames = _sampled_frame_indices(
            int(pos.shape[0]),
            frame_start_fraction=frame_start_fraction,
            max_frames=max_frames,
        )
        lipid_molecules = _dopc_lipid_molecules_from_upside_h5(h5, up_file)
        box_ang = _read_output_box_lengths_ang(h5)
        gap_samples = []
        compaction_samples = []
        for frame in frames:
            centers = []
            directions = []
            compaction = []
            for lipid in lipid_molecules:
                beads_ang = np.asarray(pos[frame, 0, lipid, :], dtype=np.float64)
                delta = beads_ang - beads_ang[0][None, :]
                if box_ang is not None:
                    delta -= box_ang[None, :] * np.round(delta / box_ang[None, :])
                beads_ang = beads_ang[0][None, :] + delta
                direction = ((beads_ang[8] + beads_ang[13]) * 0.5) - beads_ang[0]
                norm = float(np.linalg.norm(direction))
                if norm <= 1.0e-12:
                    continue
                centers.append(np.mean(beads_ang, axis=0))
                directions.append(direction / norm)
                beads_nm = (beads_ang - np.mean(beads_ang, axis=0)[None, :]) * ANGSTROM_TO_NM
                compaction.append(_dopc_tail_extension_ang(beads_nm))
            if not centers:
                continue
            centers_ang = np.asarray(centers, dtype=np.float64)
            dirs = np.asarray(directions, dtype=np.float64)
            gap = _cross_leaflet_gap_distance_for_frame(
                centers_ang,
                dirs,
                box_ang,
                radial_cutoff_ang=radial_cutoff_ang,
                face_cos_min=face_cos_min,
                smooth_weight=smooth_weight,
                fallback_gap_ang=fallback_gap_ang,
            )
            if gap.size != len(compaction):
                raise RuntimeError(
                    f"Gap/compaction sample count mismatch in {up_file} frame {frame}: "
                    f"{gap.size} vs {len(compaction)}"
                )
            gap_samples.append(gap.astype(np.float64))
            compaction_samples.append(np.asarray(compaction, dtype=np.float64))
    if not gap_samples:
        raise RuntimeError(f"No full-DOPC gap/compaction samples could be extracted from {up_file}")
    print(f"  Loaded {len(gap_samples)} aligned full-DOPC gap/compaction frames from {up_file}")
    return {
        "gap_ang": np.concatenate(gap_samples),
        "compaction_ang": np.concatenate(compaction_samples),
        "frame_count": int(len(gap_samples)),
        "sample_count": int(sum(len(x) for x in gap_samples)),
    }


def _select_cgl_compaction_reference_dataset_name(
    cg_grp: h5py.Group,
    reference_dataset_name: str | None = None,
) -> str:
    ref_dataset = str(reference_dataset_name or "").strip()
    if ref_dataset:
        if ref_dataset not in cg_grp:
            raise RuntimeError(f"cg_lipid_table lacks {ref_dataset}")
        return ref_dataset
    if "interface_ref_bead_positions_nm" in cg_grp:
        return "interface_ref_bead_positions_nm"
    if "ref_bead_positions_nm" in cg_grp:
        return "ref_bead_positions_nm"
    raise RuntimeError("cg_lipid_table lacks a compaction reference ensemble")


def _normalize_compaction_coordinate_values(
    compaction_ang: np.ndarray,
    compact_center_ang: float,
    extended_center_ang: float,
    clip: bool = True,
) -> np.ndarray:
    values = np.asarray(compaction_ang, dtype=np.float64)
    denom = float(compact_center_ang) - float(extended_center_ang)
    if abs(denom) <= 1.0e-12:
        raise ValueError(
            f"Compaction-state centers must be distinct, got "
            f"compact={compact_center_ang:.6f} A, extended={extended_center_ang:.6f} A"
        )
    coord = (values - float(extended_center_ang)) / denom
    if clip:
        coord = np.clip(coord, 0.0, 1.0)
    return coord


def _reference_compaction_state_metadata_from_ensemble(
    ref_ensemble_nm: np.ndarray,
    representative_count: int,
) -> dict:
    refs = _canonicalize_lipid_reference_ensemble_to_z(ref_ensemble_nm)
    compaction_ang = _dopc_tail_extension_series_ang(refs)
    if compaction_ang.size < 4:
        raise ValueError("Reference-ensemble compaction metadata requires at least four conformers")
    seed_states = _select_compaction_state_representatives(
        refs,
        compaction_ang,
        representative_count=representative_count,
    )
    states = _select_compaction_state_representatives_by_center(
        refs,
        compact_center_ang=float(seed_states["compact_center_ang"]),
        extended_center_ang=float(seed_states["extended_center_ang"]),
        representative_count=representative_count,
    )
    states["compaction_ang"] = np.asarray(compaction_ang, dtype=np.float64)
    return states


def _load_reference_compaction_state_metadata_from_h5(
    compaction_reference_h5_path: Path,
    reference_dataset_name: str | None = None,
    representative_count: int | None = None,
) -> dict:
    compaction_reference_h5_path = Path(compaction_reference_h5_path).expanduser().resolve()
    count = max(
        1,
        int(representative_count or _positive_int_env("UPSIDE_CGL_COMPACTION_STATE_REPRESENTATIVES", 2)),
    )
    with h5py.File(compaction_reference_h5_path, "r") as h5:
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{compaction_reference_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        ref_dataset = _select_cgl_compaction_reference_dataset_name(cg_grp, reference_dataset_name)
        ref_ensemble_nm = np.asarray(cg_grp[ref_dataset][:], dtype=np.float64)
    states = _reference_compaction_state_metadata_from_ensemble(
        ref_ensemble_nm,
        representative_count=count,
    )
    states["reference_dataset"] = ref_dataset
    return states


def _load_physical_compaction_state_centers_from_h5(
    compaction_reference_h5_path: Path,
) -> tuple[float, float]:
    compaction_reference_h5_path = Path(compaction_reference_h5_path).expanduser().resolve()
    try:
        states = _load_reference_compaction_state_metadata_from_h5(
            compaction_reference_h5_path,
            representative_count=1,
        )
        compact = float(states["compact_center_ang"])
        extended = float(states["extended_center_ang"])
    except (RuntimeError, ValueError, KeyError):
        with h5py.File(compaction_reference_h5_path, "r") as h5:
            if "cg_lipid_table/cg_lipid_compaction" not in h5:
                raise RuntimeError(
                    f"{compaction_reference_h5_path} lacks cg_lipid_table/cg_lipid_compaction"
                )
            comp_grp = h5["cg_lipid_table/cg_lipid_compaction"]
            compact = float(
                comp_grp.attrs.get(
                    "reference_compact_center_ang",
                    comp_grp.attrs.get("compact_state_center_ang", 0.0),
                )
            )
            extended = float(
                comp_grp.attrs.get(
                    "reference_extended_center_ang",
                    comp_grp.attrs.get("extended_state_center_ang", 0.0),
                )
            )
    if compact <= extended:
        compact, extended = extended, compact
    if compact - extended <= 1.0:
        raise RuntimeError(
            f"{compaction_reference_h5_path} does not contain physical compaction-state centers: "
            f"compact={compact:.6f} A, extended={extended:.6f} A"
        )
    return compact, extended


def _cg_lipid_compaction_metadata_needs_refresh(cg_grp: h5py.Group) -> bool:
    if "cg_lipid_compaction" not in cg_grp:
        return False
    comp_grp = cg_grp["cg_lipid_compaction"]
    required_pair_dsets = [
        "delta_extended_extended",
        "delta_extended_compact",
        "delta_compact_compact",
    ]
    compressed_state_center = float(comp_grp.attrs.get("compressed_state_center_ang", np.nan))
    if np.isfinite(compressed_state_center):
        required_pair_dsets.extend(
            (
                "delta_extended_compressed",
                "delta_compact_compressed",
                "delta_compressed_compressed",
            )
        )
    if "self_coeff" not in comp_grp or any(name not in comp_grp for name in required_pair_dsets):
        return True
    self_coeff = np.asarray(comp_grp["self_coeff"][:], dtype=np.float64)
    if self_coeff.size < 1 or not np.isfinite(self_coeff).all() or float(np.max(np.abs(self_coeff))) <= 1.0e-6:
        return True
    compact = float(comp_grp.attrs.get("reference_compact_center_ang", 0.0))
    extended = float(comp_grp.attrs.get("reference_extended_center_ang", 0.0))
    if compact - extended <= 1.0:
        return True
    pair_compact = float(comp_grp.attrs.get("pair_reference_compact_center_ang", np.nan))
    pair_extended = float(comp_grp.attrs.get("pair_reference_extended_center_ang", np.nan))
    if (
        not np.isfinite(pair_compact)
        or not np.isfinite(pair_extended)
        or abs(pair_compact - compact) > 1.0e-3
        or abs(pair_extended - extended) > 1.0e-3
    ):
        return True
    if "pair_compaction_state_source" not in comp_grp.attrs:
        return True
    if np.isfinite(compressed_state_center):
        pair_compressed = float(comp_grp.attrs.get("pair_reference_compressed_center_ang", np.nan))
        ref_compressed = float(comp_grp.attrs.get("reference_compressed_center_ang", np.nan))
        if not np.isfinite(pair_compressed):
            return True
        if np.isfinite(ref_compressed) and abs(pair_compressed - ref_compressed) > 1.0e-3:
            return True
    self_compaction_sample_source = _decode_h5_scalar_attr(
        comp_grp.attrs.get("self_compaction_sample_source", "")
    ).strip()
    if np.isfinite(compressed_state_center) and "target_compaction_ang" in comp_grp:
        if self_compaction_sample_source != "target_compaction_ang_projection":
            return True
    try:
        ref_dataset = _select_cgl_compaction_reference_dataset_name(cg_grp)
        ref_ensemble_nm = np.asarray(cg_grp[ref_dataset][:], dtype=np.float64)
        ref_compaction_ang = _dopc_tail_extension_series_ang(ref_ensemble_nm)
        ref_q = _normalize_compaction_coordinate_values(ref_compaction_ang, compact, extended)
    except (RuntimeError, ValueError):
        return True
    endpoint_fraction = float(np.mean((ref_q <= 1.0e-6) | (ref_q >= 1.0 - 1.0e-6)))
    return endpoint_fraction >= 0.85


def _single_cgl_compaction_group_needs_refresh(
    cg_grp: h5py.Group,
    group_name: str,
) -> bool:
    if group_name not in cg_grp:
        return False
    grp = cg_grp[group_name]
    if "delta_extended" not in grp or "delta_compact" not in grp:
        return True
    compact = float(grp.attrs.get("compact_state_center_ang", np.nan))
    extended = float(grp.attrs.get("extended_state_center_ang", np.nan))
    if abs(compact - 1.0) > 1.0e-3 or abs(extended - 0.0) > 1.0e-3:
        return True
    ref_compact = float(grp.attrs.get("reference_compact_center_ang", np.nan))
    ref_extended = float(grp.attrs.get("reference_extended_center_ang", np.nan))
    if not np.isfinite(ref_compact) or not np.isfinite(ref_extended):
        return True
    if "compaction_state_source" not in grp.attrs:
        return True
    if "cg_lipid_compaction" not in cg_grp:
        return True

    comp_grp = cg_grp["cg_lipid_compaction"]
    expected_ref_compact = float(comp_grp.attrs.get("reference_compact_center_ang", np.nan))
    expected_ref_extended = float(comp_grp.attrs.get("reference_extended_center_ang", np.nan))
    if not np.isfinite(expected_ref_compact) or not np.isfinite(expected_ref_extended):
        return True
    if abs(ref_compact - expected_ref_compact) > 1.0e-3:
        return True
    if abs(ref_extended - expected_ref_extended) > 1.0e-3:
        return True
    expected_runtime_compressed = float(comp_grp.attrs.get("compressed_state_center_ang", np.nan))
    group_runtime_compressed = float(grp.attrs.get("compressed_state_center_ang", np.nan))
    expected_ref_compressed = float(comp_grp.attrs.get("reference_compressed_center_ang", np.nan))
    group_ref_compressed = float(grp.attrs.get("reference_compressed_center_ang", np.nan))
    if np.isfinite(expected_runtime_compressed):
        if "delta_compressed" not in grp or "grid_compressed_kj_mol" not in grp:
            return True
        if not np.isfinite(group_runtime_compressed) or abs(group_runtime_compressed - expected_runtime_compressed) > 1.0e-3:
            return True
        if (
            np.isfinite(expected_ref_compressed)
            and (
                not np.isfinite(group_ref_compressed)
                or abs(group_ref_compressed - expected_ref_compressed) > 1.0e-3
            )
        ):
            return True
    elif "delta_compressed" in grp or np.isfinite(group_runtime_compressed):
        return True

    expected_prob = float(comp_grp.attrs.get("compact_state_probability", np.nan))
    group_prob = float(grp.attrs.get("compact_state_probability", np.nan))
    if np.isfinite(expected_prob):
        if not np.isfinite(group_prob) or abs(group_prob - expected_prob) > 1.0e-6:
            return True

    expected_source = _decode_h5_scalar_attr(comp_grp.attrs.get("compaction_state_source", "")).strip()
    group_source = _decode_h5_scalar_attr(grp.attrs.get("compaction_state_source", "")).strip()
    if expected_source and group_source != expected_source:
        return True

    return False


def _fit_gap_compact_probability_response(
    gap_samples_ang: np.ndarray,
    compact_probability_samples: np.ndarray,
    n_response: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
    coord_range_source: np.ndarray | None = None,
) -> dict:
    gap_samples_ang = np.asarray(gap_samples_ang, dtype=np.float64)
    compact_probability_samples = np.asarray(compact_probability_samples, dtype=np.float64)
    if gap_samples_ang.shape != compact_probability_samples.shape:
        raise ValueError("gap and compact-probability samples must have matching shapes")
    if gap_samples_ang.size < 1:
        raise ValueError("gap response fit requires at least one sample")
    coord_source = gap_samples_ang if coord_range_source is None else np.asarray(coord_range_source, dtype=np.float64)
    coord_lo = float(np.min(coord_source))
    coord_hi = float(np.max(coord_source))
    pad = max(0.05 * (coord_hi - coord_lo), 0.25)
    coord_min = max(0.0, coord_lo - pad)
    coord_max = coord_hi + pad
    hist_edges = np.linspace(coord_min, coord_max, int(n_hist) + 1, dtype=np.float64)
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    target_counts, _ = np.histogram(gap_samples_ang, bins=hist_edges)
    target_q_sum, _ = np.histogram(
        gap_samples_ang,
        bins=hist_edges,
        weights=np.clip(compact_probability_samples, 0.0, 1.0),
    )
    prior_mean = float(np.mean(np.clip(compact_probability_samples, 0.0, 1.0)))
    compact_probability = (target_q_sum + float(pseudocount) * prior_mean) / (
        target_counts.astype(np.float64) + float(pseudocount)
    )
    sampled = target_counts > 0
    if np.any(sampled):
        first = int(np.argmax(sampled))
        last = int(len(sampled) - 1 - np.argmax(sampled[::-1]))
        compact_probability[:first] = compact_probability[first]
        compact_probability[last + 1 :] = compact_probability[last]
        compact_probability[first : last + 1] = np.minimum.accumulate(
            compact_probability[first : last + 1]
        )
    else:
        compact_probability = np.minimum.accumulate(compact_probability)
    compact_probability = np.clip(compact_probability, 1.0e-3, 1.0 - 1.0e-3)

    coord_spacing = max(coord_max - coord_min, 1.0e-3) / float(int(n_response) - 3)
    knot_vector = np.zeros(int(n_response) + 4, dtype=np.float64)
    knot_vector[4:-4] = np.arange(1, int(n_response) - 3, dtype=np.float64)
    knot_vector[-4:] = knot_vector[-5]
    coord_t = 1.0 + (hist_centers - coord_min) / coord_spacing
    response_coeff = _fit_radial_bspline(
        coord_t,
        compact_probability,
        knot_vector,
        smooth=float(smooth),
    )
    return {
        "response_coeff": response_coeff.astype(np.float32),
        "response_coord_min": float(coord_min),
        "response_coord_max": float(coord_max),
        "response_coord_spacing": float(coord_spacing),
        "response_n_knot": int(n_response),
        "gap_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "compact_probability_hist": compact_probability.astype(np.float32),
        "prior_mean_probability": float(prior_mean),
    }


def _build_relative_entropy_gap_embedding_table(
    target_up_file: Path,
    model_up_file: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
    n_embedding: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    target_frames = _cgl_contact_frames_from_upside_h5(
        target_up_file,
        frame_start_fraction=target_frame_start_fraction,
        max_frames=max_frames,
    )
    model_frames = _cgl_contact_frames_from_upside_h5(
        model_up_file,
        frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
    )
    if fallback_gap_ang is None:
        fallback_gap_ang = float(radial_cutoff_ang)
    target_gap = _cross_leaflet_gap_distance_samples(
        target_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )
    model_gap = _cross_leaflet_gap_distance_samples(
        model_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )
    coord_lo = float(min(np.min(target_gap), np.min(model_gap)))
    coord_hi = float(max(np.max(target_gap), np.max(model_gap)))
    pad = max(0.05 * (coord_hi - coord_lo), 0.25)
    coord_min = max(0.0, coord_lo - pad)
    coord_max = coord_hi + pad
    hist_edges = np.linspace(coord_min, coord_max, int(n_hist) + 1, dtype=np.float64)
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    target_counts, _ = np.histogram(target_gap, bins=hist_edges)
    model_counts, _ = np.histogram(model_gap, bins=hist_edges)
    target_prob = target_counts.astype(np.float64) + float(pseudocount)
    model_prob = model_counts.astype(np.float64) + float(pseudocount)
    target_prob /= float(np.sum(target_prob))
    model_prob /= float(np.sum(model_prob))
    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    correction_kj_mol = kbt * np.log(model_prob / target_prob)
    correction_kj_mol -= float(np.average(correction_kj_mol, weights=target_prob))

    coord_spacing = max(coord_max - coord_min, 1.0e-3) / float(int(n_embedding) - 3)
    emb_knot_vector = np.zeros(int(n_embedding) + 4, dtype=np.float64)
    emb_knot_vector[4:-4] = np.arange(1, int(n_embedding) - 3, dtype=np.float64)
    emb_knot_vector[-4:] = emb_knot_vector[-5]
    emb_t = 1.0 + (hist_centers - coord_min) / coord_spacing
    embedding_coeff = _fit_radial_bspline(
        emb_t,
        correction_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        emb_knot_vector,
        smooth=float(smooth),
    )

    finite = correction_kj_mol[np.isfinite(correction_kj_mol)]
    print(
        f"  CGL gap embedding: cutoff={radial_cutoff_ang:.3f} A, "
        f"face_cos_min={face_cos_min:.3f}, gap={coord_min:.3f}..{coord_max:.3f} A, "
        f"deltaU={float(np.min(finite)):.3f}..{float(np.max(finite)):.3f} kJ/mol"
    )
    return {
        "embedding_coeff": embedding_coeff.astype(np.float32),
        "embedding_coord_min": float(coord_min),
        "embedding_coord_max": float(coord_max),
        "embedding_coord_spacing": float(coord_spacing),
        "embedding_n_knot": int(n_embedding),
        "gap_radial_cutoff_ang": float(radial_cutoff_ang),
        "gap_face_cos_min": float(face_cos_min),
        "gap_smooth_weight": int(bool(smooth_weight)),
        "gap_fallback_ang": float(fallback_gap_ang),
        "target_gap": target_gap.astype(np.float32),
        "model_gap": model_gap.astype(np.float32),
        "gap_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "model_counts": model_counts.astype(np.float32),
        "correction_kj_mol": correction_kj_mol.astype(np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "target_frame_count": int(len(target_frames)),
        "model_frame_count": int(len(model_frames)),
        "target_sample_count": int(target_gap.size),
        "model_sample_count": int(model_gap.size),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "target_upside_h5_path": str(Path(target_up_file).expanduser().resolve()),
        "model_upside_h5_path": str(Path(model_up_file).expanduser().resolve()),
    }


def add_cg_lipid_gap_embedding_table_to_h5(
    base_h5_path: Path,
    output_path: Path,
    target_upside_h5_path: Path,
    model_upside_h5_path: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
) -> None:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    embedding = _build_relative_entropy_gap_embedding_table(
        target_upside_h5_path,
        model_upside_h5_path,
        temperature_upside=temperature_upside,
        target_frame_start_fraction=target_frame_start_fraction,
        model_frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )

    def _writer(h5: h5py.File) -> None:
        with h5py.File(base_h5_path, "r") as src:
            for key, value in src.attrs.items():
                h5.attrs[key] = value
            for key in src.keys():
                src.copy(key, h5)
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_gap_embedding" in cg_grp:
            del cg_grp["cg_lipid_gap_embedding"]
        embed = cg_grp.create_group("cg_lipid_gap_embedding")
        embed.attrs["schema"] = "cg_lipid_gap_embedding_relative_entropy_v1"
        embed.attrs["source"] = "relative_entropy_cross_leaflet_gap_target_full_dopc_model_one_particle_cgl"
        embed.attrs["runtime_representation"] = "sum_i_F_of_cross_leaflet_weighted_gap_distance"
        embed.attrs["correction_layer"] = "conservative_many_body_cross_leaflet_gap_embedding"
        embed.attrs["one_particle_cgl_preserved"] = np.int32(1)
        embed.attrs["boltzmann_temperature_upside"] = np.float32(temperature_upside)
        embed.attrs["energy_conversion_kj_per_eup"] = np.float32(ENERGY_CONVERSION_KJ_PER_EUP)
        embed.attrs["length_conversion_ang_per_nm"] = np.float32(LENGTH_CONVERSION_A_PER_NM)
        for attr_name in (
            "embedding_coord_min",
            "embedding_coord_max",
            "embedding_coord_spacing",
            "embedding_n_knot",
            "gap_radial_cutoff_ang",
            "gap_face_cos_min",
            "gap_smooth_weight",
            "gap_fallback_ang",
            "pseudocount",
            "smooth",
            "target_frame_count",
            "model_frame_count",
            "target_sample_count",
            "model_sample_count",
            "target_frame_start_fraction",
            "model_frame_start_fraction",
            "target_upside_h5_path",
            "model_upside_h5_path",
        ):
            value = embedding[attr_name]
            if isinstance(value, (int, np.integer)):
                embed.attrs[attr_name] = np.int32(value)
            elif isinstance(value, (float, np.floating)):
                embed.attrs[attr_name] = np.float32(value)
            else:
                embed.attrs[attr_name] = value
        embed.create_dataset("embedding_coeff", data=embedding["embedding_coeff"])
        embed.create_dataset("target_gap", data=embedding["target_gap"])
        embed.create_dataset("model_gap", data=embedding["model_gap"])
        embed.create_dataset("gap_hist_centers", data=embedding["gap_hist_centers"])
        embed.create_dataset("target_counts", data=embedding["target_counts"])
        embed.create_dataset("model_counts", data=embedding["model_counts"])
        embed.create_dataset("correction_kj_mol", data=embedding["correction_kj_mol"])

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} with cg_lipid_gap_embedding")


def _build_target_gap_pmf_embedding_table(
    target_up_file: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
    n_embedding: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    target_frames = _cgl_contact_frames_from_upside_h5(
        target_up_file,
        frame_start_fraction=target_frame_start_fraction,
        max_frames=max_frames,
    )
    if fallback_gap_ang is None:
        fallback_gap_ang = float(radial_cutoff_ang)
    target_gap = _cross_leaflet_gap_distance_samples(
        target_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )
    coord_lo = float(np.min(target_gap))
    coord_hi = float(np.max(target_gap))
    pad = max(0.05 * (coord_hi - coord_lo), 0.25)
    coord_min = max(0.0, coord_lo - pad)
    coord_max = coord_hi + pad
    hist_edges = np.linspace(coord_min, coord_max, int(n_hist) + 1, dtype=np.float64)
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    target_counts, _ = np.histogram(target_gap, bins=hist_edges)
    target_prob = target_counts.astype(np.float64) + float(pseudocount)
    target_prob /= float(np.sum(target_prob))
    kbt = float(temperature_upside) * ENERGY_CONVERSION_KJ_PER_EUP
    pmf_kj_mol = -kbt * np.log(np.maximum(target_prob, 1.0e-12))
    pmf_kj_mol -= float(np.average(pmf_kj_mol, weights=target_prob))

    coord_spacing = max(coord_max - coord_min, 1.0e-3) / float(int(n_embedding) - 3)
    emb_knot_vector = np.zeros(int(n_embedding) + 4, dtype=np.float64)
    emb_knot_vector[4:-4] = np.arange(1, int(n_embedding) - 3, dtype=np.float64)
    emb_knot_vector[-4:] = emb_knot_vector[-5]
    emb_t = 1.0 + (hist_centers - coord_min) / coord_spacing
    embedding_coeff = _fit_radial_bspline(
        emb_t,
        pmf_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        emb_knot_vector,
        smooth=float(smooth),
    )

    print(
        f"  CGL target gap PMF embedding: cutoff={radial_cutoff_ang:.3f} A, "
        f"face_cos_min={face_cos_min:.3f}, gap={coord_min:.3f}..{coord_max:.3f} A, "
        f"pmf={float(np.min(pmf_kj_mol)):.3f}..{float(np.max(pmf_kj_mol)):.3f} kJ/mol"
    )
    return {
        "embedding_coeff": embedding_coeff.astype(np.float32),
        "embedding_coord_min": float(coord_min),
        "embedding_coord_max": float(coord_max),
        "embedding_coord_spacing": float(coord_spacing),
        "embedding_n_knot": int(n_embedding),
        "gap_radial_cutoff_ang": float(radial_cutoff_ang),
        "gap_face_cos_min": float(face_cos_min),
        "gap_smooth_weight": int(bool(smooth_weight)),
        "gap_fallback_ang": float(fallback_gap_ang),
        "target_gap": target_gap.astype(np.float32),
        "gap_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "pmf_kj_mol": pmf_kj_mol.astype(np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "target_frame_count": int(len(target_frames)),
        "target_sample_count": int(target_gap.size),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "target_upside_h5_path": str(Path(target_up_file).expanduser().resolve()),
    }


def add_cg_lipid_target_gap_pmf_embedding_table_to_h5(
    base_h5_path: Path,
    output_path: Path,
    target_upside_h5_path: Path,
    temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    target_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
) -> None:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    embedding = _build_target_gap_pmf_embedding_table(
        target_upside_h5_path,
        temperature_upside=temperature_upside,
        target_frame_start_fraction=target_frame_start_fraction,
        max_frames=max_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )

    def _writer(h5: h5py.File) -> None:
        with h5py.File(base_h5_path, "r") as src:
            for key, value in src.attrs.items():
                h5.attrs[key] = value
            for key in src.keys():
                src.copy(key, h5)
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_gap_embedding" in cg_grp:
            del cg_grp["cg_lipid_gap_embedding"]
        embed = cg_grp.create_group("cg_lipid_gap_embedding")
        embed.attrs["schema"] = "cg_lipid_gap_embedding_target_pmf_v1"
        embed.attrs["source"] = "target_pmf_cross_leaflet_gap_full_dopc"
        embed.attrs["runtime_representation"] = "sum_i_F_of_cross_leaflet_weighted_gap_distance"
        embed.attrs["correction_layer"] = "conservative_many_body_cross_leaflet_gap_target_pmf"
        embed.attrs["one_particle_cgl_preserved"] = np.int32(1)
        embed.attrs["boltzmann_temperature_upside"] = np.float32(temperature_upside)
        embed.attrs["energy_conversion_kj_per_eup"] = np.float32(ENERGY_CONVERSION_KJ_PER_EUP)
        embed.attrs["length_conversion_ang_per_nm"] = np.float32(LENGTH_CONVERSION_A_PER_NM)
        for attr_name in (
            "embedding_coord_min",
            "embedding_coord_max",
            "embedding_coord_spacing",
            "embedding_n_knot",
            "gap_radial_cutoff_ang",
            "gap_face_cos_min",
            "gap_smooth_weight",
            "gap_fallback_ang",
            "pseudocount",
            "smooth",
            "target_frame_count",
            "target_sample_count",
            "target_frame_start_fraction",
            "target_upside_h5_path",
        ):
            value = embedding[attr_name]
            if isinstance(value, (int, np.integer)):
                embed.attrs[attr_name] = np.int32(value)
            elif isinstance(value, (float, np.floating)):
                embed.attrs[attr_name] = np.float32(value)
            else:
                embed.attrs[attr_name] = value
        embed.create_dataset("embedding_coeff", data=embedding["embedding_coeff"])
        embed.create_dataset("target_gap", data=embedding["target_gap"])
        embed.create_dataset("gap_hist_centers", data=embedding["gap_hist_centers"])
        embed.create_dataset("target_counts", data=embedding["target_counts"])
        embed.create_dataset("pmf_kj_mol", data=embedding["pmf_kj_mol"])

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} with target gap PMF embedding")


def _build_gap_compaction_response_table(
    target_up_file: Path,
    model_up_file: Path,
    target_compaction_reference_h5_path: Path | None = None,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
    n_response: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    if fallback_gap_ang is None:
        fallback_gap_ang = float(radial_cutoff_ang)
    model_frames = _cgl_contact_frames_from_upside_h5(
        model_up_file,
        frame_start_fraction=model_frame_start_fraction,
        max_frames=max_frames,
    )
    model_gap = _cross_leaflet_gap_distance_samples(
        model_frames,
        radial_cutoff_ang=radial_cutoff_ang,
        face_cos_min=face_cos_min,
        smooth_weight=smooth_weight,
        fallback_gap_ang=fallback_gap_ang,
    )
    if target_compaction_reference_h5_path is None:
        target_frames = _cgl_contact_frames_from_upside_h5(
            target_up_file,
            frame_start_fraction=target_frame_start_fraction,
            max_frames=max_frames,
        )
        target_gap = _cross_leaflet_gap_distance_samples(
            target_frames,
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
            fallback_gap_ang=fallback_gap_ang,
        )
        coord_lo = float(min(np.min(target_gap), np.min(model_gap)))
        coord_hi = float(max(np.max(target_gap), np.max(model_gap)))
        pad = max(0.05 * (coord_hi - coord_lo), 0.25)
        coord_min = max(0.0, coord_lo - pad)
        coord_max = coord_hi + pad
        hist_edges = np.linspace(coord_min, coord_max, int(n_hist) + 1, dtype=np.float64)
        hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
        target_counts, _ = np.histogram(target_gap, bins=hist_edges)
        model_counts, _ = np.histogram(model_gap, bins=hist_edges)
        target_prob = target_counts.astype(np.float64) + float(pseudocount)
        model_prob = model_counts.astype(np.float64) + float(pseudocount)
        target_prob /= float(np.sum(target_prob))
        model_prob /= float(np.sum(model_prob))
        compact_probability = target_prob / np.maximum(target_prob + model_prob, 1.0e-12)
        sampled = (target_counts + model_counts) > 0
        if np.any(sampled):
            first = int(np.argmax(sampled))
            last = int(len(sampled) - 1 - np.argmax(sampled[::-1]))
            compact_probability[:first] = compact_probability[first]
            compact_probability[last + 1 :] = compact_probability[last]
            compact_probability[first : last + 1] = np.minimum.accumulate(
                compact_probability[first : last + 1]
            )
        else:
            compact_probability = np.minimum.accumulate(compact_probability)
        compact_probability = np.clip(compact_probability, 1.0e-3, 1.0 - 1.0e-3)

        coord_spacing = max(coord_max - coord_min, 1.0e-3) / float(int(n_response) - 3)
        knot_vector = np.zeros(int(n_response) + 4, dtype=np.float64)
        knot_vector[4:-4] = np.arange(1, int(n_response) - 3, dtype=np.float64)
        knot_vector[-4:] = knot_vector[-5]
        coord_t = 1.0 + (hist_centers - coord_min) / coord_spacing
        response_coeff = _fit_radial_bspline(
            coord_t,
            compact_probability,
            knot_vector,
            smooth=float(smooth),
        )
        response_mode = "gap_compact_probability_classifier"
        prior_mean_probability = float(np.mean(compact_probability))
        target_compaction_ang = None
        target_compact_probability = None
        reference_compact_center_ang = None
        reference_extended_center_ang = None
        target_frame_count = int(len(target_frames))
        target_sample_count = int(target_gap.size)
    else:
        target_data = _full_dopc_gap_and_compaction_samples_from_upside_h5(
            target_up_file,
            frame_start_fraction=target_frame_start_fraction,
            max_frames=max_frames,
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
            fallback_gap_ang=fallback_gap_ang,
        )
        target_gap = np.asarray(target_data["gap_ang"], dtype=np.float64)
        target_compaction_ang = np.asarray(target_data["compaction_ang"], dtype=np.float64)
        reference_compact_center_ang, reference_extended_center_ang = _load_physical_compaction_state_centers_from_h5(
            target_compaction_reference_h5_path
        )
        target_compact_probability = np.clip(
            (target_compaction_ang - reference_extended_center_ang)
            / max(reference_compact_center_ang - reference_extended_center_ang, 1.0e-12),
            0.0,
            1.0,
        )
        fit = _fit_gap_compact_probability_response(
            target_gap,
            target_compact_probability,
            n_response=n_response,
            n_hist=n_hist,
            pseudocount=pseudocount,
            smooth=smooth,
            coord_range_source=np.concatenate([target_gap, model_gap]),
        )
        response_coeff = np.asarray(fit["response_coeff"], dtype=np.float64)
        coord_min = float(fit["response_coord_min"])
        coord_max = float(fit["response_coord_max"])
        coord_spacing = float(fit["response_coord_spacing"])
        hist_centers = np.asarray(fit["gap_hist_centers"], dtype=np.float64)
        target_counts = np.asarray(fit["target_counts"], dtype=np.float64)
        compact_probability = np.asarray(fit["compact_probability_hist"], dtype=np.float64)
        hist_edges = _grid_edges_from_centers(hist_centers, lo=coord_min, hi=coord_max)
        model_counts, _ = np.histogram(model_gap, bins=hist_edges)
        response_mode = "gap_target_compact_state_projection"
        prior_mean_probability = float(fit["prior_mean_probability"])
        target_frame_count = int(target_data["frame_count"])
        target_sample_count = int(target_data["sample_count"])

    print(
        f"  CGL gap compaction response: mode={response_mode}, cutoff={radial_cutoff_ang:.3f} A, "
        f"face_cos_min={face_cos_min:.3f}, gap={coord_min:.3f}..{coord_max:.3f} A, "
        f"compact_prob={float(np.min(compact_probability)):.3f}..{float(np.max(compact_probability)):.3f}"
    )
    out = {
        "response_coeff": np.asarray(response_coeff, dtype=np.float32),
        "response_coord_min": float(coord_min),
        "response_coord_max": float(coord_max),
        "response_coord_spacing": float(coord_spacing),
        "response_n_knot": int(n_response),
        "gap_radial_cutoff_ang": float(radial_cutoff_ang),
        "gap_face_cos_min": float(face_cos_min),
        "gap_smooth_weight": int(bool(smooth_weight)),
        "gap_fallback_ang": float(fallback_gap_ang),
        "target_gap": target_gap.astype(np.float32),
        "model_gap": model_gap.astype(np.float32),
        "gap_hist_centers": hist_centers.astype(np.float32),
        "target_counts": target_counts.astype(np.float32),
        "model_counts": model_counts.astype(np.float32),
        "compact_probability_hist": compact_probability.astype(np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "prior_mean_probability": float(prior_mean_probability),
        "implicit_response_mode": response_mode,
        "target_frame_count": int(target_frame_count),
        "model_frame_count": int(len(model_frames)),
        "target_sample_count": int(target_sample_count),
        "model_sample_count": int(model_gap.size),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "target_upside_h5_path": str(Path(target_up_file).expanduser().resolve()),
        "model_upside_h5_path": str(Path(model_up_file).expanduser().resolve()),
    }
    if target_compaction_reference_h5_path is not None:
        out["target_compaction_reference_h5_path"] = str(
            Path(target_compaction_reference_h5_path).expanduser().resolve()
        )
        out["target_compaction_ang"] = target_compaction_ang.astype(np.float32)
        out["target_compact_probability"] = target_compact_probability.astype(np.float32)
        out["reference_compact_center_ang"] = float(reference_compact_center_ang)
        out["reference_extended_center_ang"] = float(reference_extended_center_ang)
    return out


def _build_gap_compaction_response_from_stored_h5(
    base_h5_path: Path,
    target_compaction_reference_h5_path: Path | None = None,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
    n_response: int = 32,
    n_hist: int = 96,
    pseudocount: float = 0.5,
    smooth: float = 1.0e-3,
) -> dict:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    fallback_gap = float(fallback_gap_ang if fallback_gap_ang is not None else radial_cutoff_ang)
    with h5py.File(base_h5_path, "r") as h5:
        if "cg_lipid_table/cg_lipid_compaction" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table/cg_lipid_compaction")
        comp_grp = h5["cg_lipid_table/cg_lipid_compaction"]
        target_gap = np.asarray(comp_grp["target_gap"][:], dtype=np.float64)
        model_gap = np.asarray(comp_grp["model_gap"][:], dtype=np.float64)
        target_frame_count = int(comp_grp.attrs.get("target_frame_count", 0))
        model_frame_count = int(comp_grp.attrs.get("model_frame_count", 0))
        target_sample_count = int(comp_grp.attrs.get("target_sample_count", target_gap.size))
        model_sample_count = int(comp_grp.attrs.get("model_sample_count", model_gap.size))
        target_frame_start_fraction = float(comp_grp.attrs.get("target_frame_start_fraction", 0.5))
        model_frame_start_fraction = float(comp_grp.attrs.get("model_frame_start_fraction", 0.5))
        target_upside_h5_path = comp_grp.attrs.get("target_upside_h5_path", "")
        model_upside_h5_path = comp_grp.attrs.get("model_upside_h5_path", "")
        if target_compaction_reference_h5_path is None:
            raise RuntimeError(
                "stored-sample gap-response rebuild requires target_compaction_reference_h5_path"
            )
        if "target_compaction_ang" not in comp_grp:
            raise RuntimeError(f"{base_h5_path} lacks stored target_compaction_ang samples")
        target_compaction_ang = np.asarray(comp_grp["target_compaction_ang"][:], dtype=np.float64)

    reference_compact_center_ang, reference_extended_center_ang = _load_physical_compaction_state_centers_from_h5(
        target_compaction_reference_h5_path
    )
    target_compact_probability = _normalize_compaction_coordinate_values(
        target_compaction_ang,
        compact_center_ang=reference_compact_center_ang,
        extended_center_ang=reference_extended_center_ang,
    )
    fit = _fit_gap_compact_probability_response(
        target_gap,
        target_compact_probability,
        n_response=n_response,
        n_hist=n_hist,
        pseudocount=pseudocount,
        smooth=smooth,
        coord_range_source=np.concatenate([target_gap, model_gap]),
    )
    hist_centers = np.asarray(fit["gap_hist_centers"], dtype=np.float64)
    hist_edges = _grid_edges_from_centers(
        hist_centers,
        lo=float(fit["response_coord_min"]),
        hi=float(fit["response_coord_max"]),
    )
    model_counts, _ = np.histogram(model_gap, bins=hist_edges)
    return {
        "response_coeff": np.asarray(fit["response_coeff"], dtype=np.float32),
        "response_coord_min": float(fit["response_coord_min"]),
        "response_coord_max": float(fit["response_coord_max"]),
        "response_coord_spacing": float(fit["response_coord_spacing"]),
        "response_n_knot": int(n_response),
        "gap_radial_cutoff_ang": float(radial_cutoff_ang),
        "gap_face_cos_min": float(face_cos_min),
        "gap_smooth_weight": int(bool(smooth_weight)),
        "gap_fallback_ang": fallback_gap,
        "target_gap": target_gap.astype(np.float32),
        "model_gap": model_gap.astype(np.float32),
        "gap_hist_centers": hist_centers.astype(np.float32),
        "target_counts": np.asarray(fit["target_counts"], dtype=np.float32),
        "model_counts": np.asarray(model_counts, dtype=np.float32),
        "compact_probability_hist": np.asarray(fit["compact_probability_hist"], dtype=np.float32),
        "pseudocount": float(pseudocount),
        "smooth": float(smooth),
        "prior_mean_probability": float(fit["prior_mean_probability"]),
        "implicit_response_mode": "gap_target_compact_state_projection_stored_samples",
        "target_frame_count": int(target_frame_count),
        "model_frame_count": int(model_frame_count),
        "target_sample_count": int(target_sample_count),
        "model_sample_count": int(model_sample_count),
        "target_frame_start_fraction": float(target_frame_start_fraction),
        "model_frame_start_fraction": float(model_frame_start_fraction),
        "target_upside_h5_path": target_upside_h5_path,
        "model_upside_h5_path": model_upside_h5_path,
        "target_compaction_reference_h5_path": str(
            Path(target_compaction_reference_h5_path).expanduser().resolve()
        ),
        "target_compaction_ang": target_compaction_ang.astype(np.float32),
        "target_compact_probability": target_compact_probability.astype(np.float32),
        "reference_compact_center_ang": float(reference_compact_center_ang),
        "reference_extended_center_ang": float(reference_extended_center_ang),
    }


def add_cg_lipid_gap_compaction_response_to_h5(
    base_h5_path: Path,
    output_path: Path,
    target_upside_h5_path: Path,
    model_upside_h5_path: Path,
    target_compaction_reference_h5_path: Path | None = None,
    target_frame_start_fraction: float = 0.5,
    model_frame_start_fraction: float = 0.5,
    max_frames: int = 400,
    radial_cutoff_ang: float = 25.0,
    face_cos_min: float = 0.5,
    smooth_weight: bool = True,
    fallback_gap_ang: float | None = None,
    single_cgl_compaction: bool = True,
    dry_ff_path: Path | None = None,
    lipids_itp_path: Path | None = None,
    martinize_path: Path | None = None,
    sidechain_lib_path: Path | None = None,
    forcefield_name: str = "martini22",
) -> None:
    base_h5_path = Path(base_h5_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    try:
        response = _build_gap_compaction_response_table(
            target_upside_h5_path,
            model_upside_h5_path,
            target_compaction_reference_h5_path=target_compaction_reference_h5_path,
            target_frame_start_fraction=target_frame_start_fraction,
            model_frame_start_fraction=model_frame_start_fraction,
            max_frames=max_frames,
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
            fallback_gap_ang=fallback_gap_ang,
        )
    except FileNotFoundError as exc:
        print(
            "  Gap-response rebuild fallback: "
            f"using stored H5 samples because an original response source is unavailable ({exc})"
        )
        response = _build_gap_compaction_response_from_stored_h5(
            base_h5_path=base_h5_path,
            target_compaction_reference_h5_path=target_compaction_reference_h5_path,
            radial_cutoff_ang=radial_cutoff_ang,
            face_cos_min=face_cos_min,
            smooth_weight=smooth_weight,
            fallback_gap_ang=fallback_gap_ang,
        )
    single_cgl_corrections = None
    if single_cgl_compaction:
        with h5py.File(base_h5_path, "r") as src:
            cg_grp = src.get("cg_lipid_table")
            need_sc = bool(
                cg_grp is not None
                and _single_cgl_compaction_group_needs_refresh(cg_grp, "cg_lipid_sc")
            )
            need_target = bool(
                cg_grp is not None
                and _single_cgl_compaction_group_needs_refresh(cg_grp, "cg_lipid_target")
            )
            need_compaction_metadata = bool(
                cg_grp is not None and _cg_lipid_compaction_metadata_needs_refresh(cg_grp)
            )
            if need_compaction_metadata:
                need_sc = bool(cg_grp is not None and "cg_lipid_sc" in cg_grp)
                need_target = bool(cg_grp is not None and "cg_lipid_target" in cg_grp)
        if need_sc or need_target or need_compaction_metadata:
            single_cgl_corrections = _build_single_cgl_compaction_corrections_from_base_h5(
                base_h5_path=base_h5_path,
                dry_ff_path=dry_ff_path,
                lipids_itp_path=lipids_itp_path,
                martinize_path=martinize_path,
                sidechain_lib_path=sidechain_lib_path,
                forcefield_name=forcefield_name,
                include_sc=need_sc,
                include_target=need_target,
            )

    def _writer(h5: h5py.File) -> None:
        with h5py.File(base_h5_path, "r") as src:
            for key, value in src.attrs.items():
                h5.attrs[key] = value
            for key in src.keys():
                src.copy(key, h5)
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_compaction" not in cg_grp or "cg_lipid_pair" not in cg_grp:
            raise RuntimeError(f"{base_h5_path} lacks compaction-aware CGL pair tables")
        comp_grp = cg_grp["cg_lipid_compaction"]
        if "gap_response_coeff" in comp_grp:
            del comp_grp["gap_response_coeff"]
        for dataset_name in (
            "gap_hist_centers",
            "target_gap",
            "model_gap",
            "target_counts",
            "model_counts",
            "compact_probability_hist",
            "target_compaction_ang",
            "target_compact_probability",
        ):
            if dataset_name in comp_grp:
                del comp_grp[dataset_name]
        comp_grp.attrs["implicit_response_mode"] = response["implicit_response_mode"]
        comp_grp.attrs["gap_response_coord_min_ang"] = np.float32(response["response_coord_min"])
        comp_grp.attrs["gap_response_coord_max_ang"] = np.float32(response["response_coord_max"])
        comp_grp.attrs["gap_response_coord_spacing_ang"] = np.float32(response["response_coord_spacing"])
        comp_grp.attrs["gap_response_n_knot"] = np.int32(response["response_n_knot"])
        comp_grp.attrs["gap_response_radial_cutoff_ang"] = np.float32(response["gap_radial_cutoff_ang"])
        comp_grp.attrs["gap_response_face_cos_min"] = np.float32(response["gap_face_cos_min"])
        comp_grp.attrs["gap_response_smooth_weight"] = np.int32(response["gap_smooth_weight"])
        comp_grp.attrs["gap_response_fallback_ang"] = np.float32(response["gap_fallback_ang"])
        for attr_name in (
            "pseudocount",
            "smooth",
            "target_frame_count",
            "model_frame_count",
            "target_sample_count",
            "model_sample_count",
            "target_frame_start_fraction",
            "model_frame_start_fraction",
            "target_upside_h5_path",
            "model_upside_h5_path",
            "prior_mean_probability",
        ):
            value = response[attr_name]
            if isinstance(value, (int, np.integer)):
                comp_grp.attrs[attr_name] = np.int32(value)
            elif isinstance(value, (float, np.floating)):
                comp_grp.attrs[attr_name] = np.float32(value)
            else:
                comp_grp.attrs[attr_name] = value
        for attr_name in (
            "target_compaction_reference_h5_path",
            "reference_compact_center_ang",
            "reference_extended_center_ang",
        ):
            if attr_name not in response:
                continue
            value = response[attr_name]
            if isinstance(value, (float, np.floating)):
                comp_grp.attrs[attr_name] = np.float32(value)
            else:
                comp_grp.attrs[attr_name] = value
        comp_grp.create_dataset("gap_response_coeff", data=response["response_coeff"])
        comp_grp.create_dataset("gap_hist_centers", data=response["gap_hist_centers"])
        comp_grp.create_dataset("target_gap", data=response["target_gap"])
        comp_grp.create_dataset("model_gap", data=response["model_gap"])
        comp_grp.create_dataset("target_counts", data=response["target_counts"])
        comp_grp.create_dataset("model_counts", data=response["model_counts"])
        comp_grp.create_dataset("compact_probability_hist", data=response["compact_probability_hist"])
        if "target_compaction_ang" in response:
            comp_grp.create_dataset("target_compaction_ang", data=response["target_compaction_ang"])
            comp_grp.create_dataset(
                "target_compact_probability",
                data=response["target_compact_probability"],
            )
        if single_cgl_corrections is not None:
            _apply_single_cgl_compaction_corrections_to_h5(cg_grp, single_cgl_corrections)

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} with gap-driven compaction response")


def _default_dry_martini_repo_paths() -> dict[str, Path]:
    repo_root = Path(__file__).resolve().parent.parent
    return {
        "dry_ff_path": repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1.itp",
        "lipids_itp_path": repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1_lipids.itp",
        "martinize_path": repo_root / "py" / "martinize.py",
        "sidechain_lib_path": repo_root / "parameters" / "ff_2.1" / "sidechain.h5",
    }


def _decode_h5_string_array(values: np.ndarray) -> list[str]:
    out = []
    for value in np.asarray(values):
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(value).strip())
    return out


def _decode_h5_scalar_attr(value) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip()
    return str(value).strip()


def _build_single_cgl_compaction_corrections_from_base_h5(
    base_h5_path: Path,
    dry_ff_path: Path | None = None,
    lipids_itp_path: Path | None = None,
    martinize_path: Path | None = None,
    sidechain_lib_path: Path | None = None,
    forcefield_name: str = "martini22",
    include_sc: bool = True,
    include_target: bool = True,
    include_pair_compaction: bool = True,
    reference_dataset_name: str | None = None,
) -> dict:
    defaults = _default_dry_martini_repo_paths()
    dry_ff_path = Path(dry_ff_path or defaults["dry_ff_path"]).expanduser().resolve()
    lipids_itp_path = Path(lipids_itp_path or defaults["lipids_itp_path"]).expanduser().resolve()
    martinize_path = Path(martinize_path or defaults["martinize_path"]).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path or defaults["sidechain_lib_path"]).expanduser().resolve()
    base_h5_path = Path(base_h5_path).expanduser().resolve()

    from martini_itp_reader import parse_dopc_from_itp

    _, pair_params = parse_dry_forcefield(dry_ff_path)
    dopc = parse_dopc_from_itp(lipids_itp_path)
    compaction_representatives = _positive_int_env(
        "UPSIDE_CGL_COMPACTION_STATE_REPRESENTATIVES", 2
    )
    compaction_self_bins = _positive_int_env(
        "UPSIDE_CGL_COMPACTION_SELF_BINS", 12
    )
    compaction_state_source = (
        os.environ.get("UPSIDE_CGL_COMPACTION_STATE_SOURCE", "auto").strip().lower()
        or "auto"
    )
    compaction_pool_conformers = _positive_int_env(
        "UPSIDE_CGL_COMPACTION_POOL_CONFORMERS", 32
    )
    compaction_burnin_steps = _positive_int_env(
        "UPSIDE_CGL_COMPACTION_BURNIN_STEPS", 20000
    )
    compaction_steps_per_conf = _positive_int_env(
        "UPSIDE_CGL_COMPACTION_STEPS_PER_CONFORMER", 500
    )
    isolated_compaction_refs_nm = None

    def ensure_isolated_compaction_refs_nm() -> np.ndarray:
        nonlocal isolated_compaction_refs_nm
        if isolated_compaction_refs_nm is None:
            isolated_compaction_refs_nm = sample_isolated_dopc_bonded_conformers(
                dopc,
                lipids_itp_path=lipids_itp_path,
                pair_params=pair_params,
                conformer_count=compaction_pool_conformers,
                temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
                seed=1777,
                mc_burnin_steps=compaction_burnin_steps,
                mc_steps_per_conformer=compaction_steps_per_conf,
            )
        return isolated_compaction_refs_nm

    with h5py.File(base_h5_path, "r") as h5:
        if "cg_lipid_table" not in h5:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_table")
        cg_grp = h5["cg_lipid_table"]
        if "cg_lipid_compaction" not in cg_grp:
            raise RuntimeError(f"{base_h5_path} lacks cg_lipid_compaction")
        ref_dataset = _select_cgl_compaction_reference_dataset_name(cg_grp, reference_dataset_name)
        ref_ensemble_nm = np.asarray(cg_grp[ref_dataset][:], dtype=np.float64)
        bead_charges = list(np.asarray(cg_grp["bead_charges"][:], dtype=np.float64))
        comp_grp = cg_grp["cg_lipid_compaction"]
        compaction_tau = float(comp_grp.attrs.get("thermostat_timescale", 5.0))
        stored_compact_center_ang = float(comp_grp.attrs.get("reference_compact_center_ang", 0.0))
        stored_extended_center_ang = float(comp_grp.attrs.get("reference_extended_center_ang", 0.0))
        if stored_compact_center_ang <= stored_extended_center_ang:
            stored_compact_center_ang, stored_extended_center_ang = (
                stored_extended_center_ang,
                stored_compact_center_ang,
            )
        stored_compact_probability = float(comp_grp.attrs.get("compact_state_probability", 0.5))
        compact_center_ang = 0.0
        extended_center_ang = 0.0
        compact_probability = 0.5
        compaction_coord_samples = None
        compaction_states = None
        compaction_state_source_used = None
        compressed_center_ang = None
        compressed_state_coord = None
        clip_compaction_state_values = False
        if compaction_state_source in (
            "reference_ensemble",
            "reference",
            "auto",
        ):
            try:
                compaction_states = _reference_compaction_state_metadata_from_ensemble(
                    ref_ensemble_nm,
                    representative_count=compaction_representatives,
                )
                compact_center_ang = float(compaction_states["compact_center_ang"])
                extended_center_ang = float(compaction_states["extended_center_ang"])
                compact_probability = float(compaction_states["compact_probability"])
                compaction_coord_samples = np.asarray(compaction_states["compaction_ang"], dtype=np.float64)
                compaction_state_source_used = f"{ref_dataset}_quantile_center_matched"
                clip_compaction_state_values = True
            except ValueError as exc:
                if compaction_state_source in ("reference_ensemble", "reference"):
                    raise
                print(
                    "  Single-CGL compaction retrofit: "
                    f"falling back from {ref_dataset} center-matched states ({exc})"
                )
        if compaction_states is None and compaction_state_source in ("stored_contract", "stored", "base_h5_contract", "auto"):
            if stored_compact_center_ang - stored_extended_center_ang > 1.0:
                compaction_states = _select_compaction_state_representatives_by_center(
                    ref_ensemble_nm,
                    compact_center_ang=stored_compact_center_ang,
                    extended_center_ang=stored_extended_center_ang,
                    representative_count=compaction_representatives,
                    compact_probability=stored_compact_probability,
                )
                compact_center_ang = float(stored_compact_center_ang)
                extended_center_ang = float(stored_extended_center_ang)
                compact_probability = float(stored_compact_probability)
                compaction_coord_samples = _dopc_tail_extension_series_ang(ref_ensemble_nm)
                compaction_state_source_used = f"{ref_dataset}_stored_contract_matched"
                clip_compaction_state_values = False
            elif compaction_state_source != "auto":
                raise RuntimeError(
                    "Stored-contract compaction-state repair requires valid "
                    "reference_compact_center_ang / reference_extended_center_ang metadata"
                )
        if compaction_states is None:
            compaction_refs_nm = ensure_isolated_compaction_refs_nm()
            compaction_values_ang = _dopc_tail_extension_series_ang(compaction_refs_nm)
            compaction_states = _select_compaction_state_representatives(
                compaction_refs_nm,
                compaction_values_ang,
                representative_count=compaction_representatives,
            )
            compact_center_ang = float(compaction_states["compact_center_ang"])
            extended_center_ang = float(compaction_states["extended_center_ang"])
            compact_probability = float(compaction_states["compact_probability"])
            compaction_coord_samples = np.asarray(compaction_values_ang, dtype=np.float64)
            compaction_state_source_used = "isolated_dopc_mc_center_quantile"
            clip_compaction_state_values = True
        if compaction_coord_samples is None:
            raise RuntimeError("Failed to determine compaction-coordinate samples for CGL metadata repair")
        ref_compaction_values_ang = _dopc_tail_extension_series_ang(ref_ensemble_nm)
        compressed_supplemental_refs_nm = None
        sample_compaction_support_ang = None
        if "target_compaction_ang" in comp_grp:
            sample_compaction_support_ang = np.asarray(
                comp_grp["target_compaction_ang"][:],
                dtype=np.float64,
            )
        if int(np.count_nonzero(ref_compaction_values_ang > compact_center_ang + 1.0e-6)) < max(2, compaction_representatives):
            extra_ref_pools = []
            for dataset_name in ("sc_interface_ref_bead_positions_nm", "ref_bead_positions_nm"):
                if dataset_name != ref_dataset and dataset_name in cg_grp:
                    extra_ref_pools.append(np.asarray(cg_grp[dataset_name][:], dtype=np.float64))
            extra_ref_pools.append(np.asarray(ensure_isolated_compaction_refs_nm(), dtype=np.float64))
            if extra_ref_pools:
                compressed_supplemental_refs_nm = np.concatenate(extra_ref_pools, axis=0)
        compaction_states = _augment_compaction_states_with_compressed_branch(
            ref_ensemble_nm,
            compaction_states,
            representative_count=compaction_representatives,
            supplemental_refs_nm=compressed_supplemental_refs_nm,
            sample_compaction_ang=sample_compaction_support_ang,
        )
        compressed_center_ang = compaction_states.get("compressed_center_ang")
        if compressed_center_ang is not None and float(compressed_center_ang) > compact_center_ang + 1.0e-6:
            compressed_state_coord = float(
                _normalize_compaction_coordinate_values(
                    np.asarray([compressed_center_ang], dtype=np.float64),
                    compact_center_ang=compact_center_ang,
                    extended_center_ang=extended_center_ang,
                    clip=False,
                )[0]
            )
            clip_compaction_state_values = False
        correction_compact_probability = float(compact_probability)
        if 1.0e-6 < stored_compact_probability < 1.0 - 1.0e-6:
            correction_compact_probability = float(stored_compact_probability)
        self_compaction_coord_samples = np.asarray(compaction_coord_samples, dtype=np.float64)
        self_compaction_sample_source = "compaction_state_samples"
        if compressed_state_coord is not None and "target_compaction_ang" in comp_grp:
            self_compaction_coord_samples = np.asarray(
                comp_grp["target_compaction_ang"][:],
                dtype=np.float64,
            )
            self_compaction_sample_source = "target_compaction_ang_projection"
        compaction_state_values = _normalize_compaction_coordinate_values(
            self_compaction_coord_samples,
            compact_center_ang=compact_center_ang,
            extended_center_ang=extended_center_ang,
            clip=clip_compaction_state_values,
        )
        if compressed_center_ang is not None and float(compressed_center_ang) > compact_center_ang + 1.0e-6:
            probability_compaction_values_ang = np.asarray(self_compaction_coord_samples, dtype=np.float64)
            ext_cmp_midpoint = 0.5 * (extended_center_ang + compact_center_ang)
            cmp_x_midpoint = 0.5 * (compact_center_ang + float(compressed_center_ang))
            extended_mask = probability_compaction_values_ang < ext_cmp_midpoint
            compressed_mask = probability_compaction_values_ang >= cmp_x_midpoint
            compact_middle_mask = (~extended_mask) & (~compressed_mask)
            compaction_states["extended_probability"] = float(np.mean(extended_mask))
            compaction_states["compact_middle_probability"] = float(np.mean(compact_middle_mask))
            compaction_states["compressed_probability"] = float(np.mean(compressed_mask))
        compaction_self = _fit_compaction_self_pmf(
            compaction_state_values,
            temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
            n_bin=compaction_self_bins,
            smooth=0.01,
        )
        corrections = {
            "compact_state_center_ang": 1.0,
            "extended_state_center_ang": 0.0,
            "compact_state_probability": correction_compact_probability,
            "reference_compact_center_ang": compact_center_ang,
            "reference_extended_center_ang": extended_center_ang,
            "compaction_state_source": compaction_state_source_used,
            "compaction_reference_dataset": ref_dataset,
            "self_compaction_sample_source": self_compaction_sample_source,
            "cg_lipid_compaction": {
                "self_coeff": np.asarray(compaction_self["self_coeff_eup"], dtype=np.float32),
                "pmf_centers_ang": np.asarray(compaction_self["pmf_centers_ang"], dtype=np.float32),
                "pmf_values_kj_mol": np.asarray(compaction_self["pmf_values_kj_mol"], dtype=np.float32),
                "coord_min_ang": float(compaction_self["coord_min_ang"]),
                "coord_max_ang": float(compaction_self["coord_max_ang"]),
                "coord_spacing_ang": float(compaction_self["coord_spacing_ang"]),
                "n_knot": int(compaction_self["n_knot"]),
                "effective_stiffness_eup_a2": float(compaction_self["effective_stiffness_eup_a2"]),
                "mass_up": float(
                    compaction_self["effective_stiffness_eup_a2"] * compaction_tau * compaction_tau
                ),
            },
        }
        if compressed_state_coord is not None and compressed_center_ang is not None:
            corrections["compressed_state_center_ang"] = float(compressed_state_coord)
            corrections["reference_compressed_center_ang"] = float(compressed_center_ang)

        pair_grp = cg_grp.get("cg_lipid_pair")
        if include_pair_compaction and pair_grp is not None:
            source_pair_compact_ang = float(
                comp_grp.attrs.get(
                    "pair_reference_compact_center_ang",
                    comp_grp.attrs.get("reference_compact_center_ang", np.nan),
                )
            )
            source_pair_extended_ang = float(
                comp_grp.attrs.get(
                    "pair_reference_extended_center_ang",
                    comp_grp.attrs.get("reference_extended_center_ang", np.nan),
                )
            )
            source_pair_compressed_ang = float(
                comp_grp.attrs.get(
                    "pair_reference_compressed_center_ang",
                    comp_grp.attrs.get("reference_compressed_center_ang", np.nan),
                )
            )
            has_source_pair_compressed_payload = all(
                dataset_name in comp_grp
                for dataset_name in (
                    "delta_extended_compressed",
                    "delta_compact_compressed",
                    "delta_compressed_compressed",
                )
            )
            if (
                np.isfinite(source_pair_compact_ang)
                and np.isfinite(source_pair_extended_ang)
                and source_pair_compact_ang - source_pair_extended_ang > 1.0
                and all(
                    dataset_name in comp_grp
                    for dataset_name in (
                        "delta_extended_extended",
                        "delta_extended_compact",
                        "delta_compact_compact",
                    )
                )
                and (
                    compressed_state_coord is None
                    or (
                        has_source_pair_compressed_payload
                        and np.isfinite(source_pair_compressed_ang)
                        and source_pair_compressed_ang > source_pair_compact_ang + 1.0e-6
                    )
                )
            ):
                pair_compaction = _pair_compaction_reparameterized_tables(
                    delta_extended_extended=np.asarray(
                        comp_grp["delta_extended_extended"][:],
                        dtype=np.float64,
                    ),
                    delta_extended_compact=np.asarray(
                        comp_grp["delta_extended_compact"][:],
                        dtype=np.float64,
                    ),
                    delta_compact_compact=np.asarray(
                        comp_grp["delta_compact_compact"][:],
                        dtype=np.float64,
                    ),
                    grid_extended_extended_kj_mol=np.asarray(
                        comp_grp["grid_extended_extended_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_extended_extended_kj_mol" in comp_grp else None,
                    grid_extended_compact_kj_mol=np.asarray(
                        comp_grp["grid_extended_compact_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_extended_compact_kj_mol" in comp_grp else None,
                    grid_compact_compact_kj_mol=np.asarray(
                        comp_grp["grid_compact_compact_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_compact_compact_kj_mol" in comp_grp else None,
                    grid_average_kj_mol=np.asarray(
                        comp_grp["grid_average_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_average_kj_mol" in comp_grp else None,
                    source_extended_center_ang=source_pair_extended_ang,
                    source_compact_center_ang=source_pair_compact_ang,
                    target_extended_center_ang=extended_center_ang,
                    target_compact_center_ang=compact_center_ang,
                    delta_extended_compressed=np.asarray(
                        comp_grp["delta_extended_compressed"][:],
                        dtype=np.float64,
                    ) if "delta_extended_compressed" in comp_grp else None,
                    delta_compact_compressed=np.asarray(
                        comp_grp["delta_compact_compressed"][:],
                        dtype=np.float64,
                    ) if "delta_compact_compressed" in comp_grp else None,
                    delta_compressed_compressed=np.asarray(
                        comp_grp["delta_compressed_compressed"][:],
                        dtype=np.float64,
                    ) if "delta_compressed_compressed" in comp_grp else None,
                    grid_extended_compressed_kj_mol=np.asarray(
                        comp_grp["grid_extended_compressed_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_extended_compressed_kj_mol" in comp_grp else None,
                    grid_compact_compressed_kj_mol=np.asarray(
                        comp_grp["grid_compact_compressed_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_compact_compressed_kj_mol" in comp_grp else None,
                    grid_compressed_compressed_kj_mol=np.asarray(
                        comp_grp["grid_compressed_compressed_kj_mol"][:],
                        dtype=np.float64,
                    ) if "grid_compressed_compressed_kj_mol" in comp_grp else None,
                    source_compressed_center_ang=float(
                        source_pair_compressed_ang
                    ),
                    target_compressed_center_ang=float(compressed_center_ang)
                    if compressed_state_coord is not None and compressed_center_ang is not None
                    else None,
                )
            else:
                pair_relax_state_correction = "pair_relax" in _decode_h5_scalar_attr(
                    comp_grp.attrs.get("source", "")
                )
                face_cos_min = (
                    float(comp_grp.attrs["face_cos_min"])
                    if "face_cos_min" in comp_grp.attrs
                    else None
                )
                radial_cutoff_nm = float(comp_grp.attrs.get("radial_cutoff_nm", 2.5))
                mask_mode = _decode_h5_scalar_attr(comp_grp.attrs.get("mask_mode", "binary")) or "binary"
                correction_center_mode = (
                    _decode_h5_scalar_attr(comp_grp.attrs.get("correction_center_mode", "base")) or "base"
                )
                pair_state_model = (
                    _decode_h5_scalar_attr(comp_grp.attrs.get("pair_state_model", "bilinear")) or "bilinear"
                )
                pair_energy_grid = np.asarray(pair_grp["energy_grid_raw_kj_mol"][:], dtype=np.float64)
                n_radial, n_angle1, n_angle2 = pair_energy_grid.shape
                if n_angle1 != n_angle2:
                    raise RuntimeError("CGL pair retrofit requires a square angular grid")
                pair_result_cg = {
                    "energy_grid_raw": pair_energy_grid,
                    "reference_energy_eup": float(pair_grp["reference_energy_eup"][0, 0]),
                    "r_values_nm": np.linspace(
                        float(pair_grp.attrs["fit_r_min_nm"]),
                        float(pair_grp.attrs["fit_r_max_nm"]),
                        int(n_radial),
                        dtype=np.float64,
                    ),
                    "cos_theta_grid": np.linspace(-1.0, 1.0, int(n_angle1), dtype=np.float64),
                    "n_radial": int(pair_grp.attrs["n_radial"]),
                    "n_angular": int(pair_grp.attrs["n_angular"]),
                    "knot_spacing_ang": float(pair_grp.attrs["knot_spacing_ang"]),
                    "fit_smooth": float(pair_grp.attrs.get("fit_smooth", 0.1)),
                    "azimuthal_count": int(pair_grp.attrs.get("azimuthal_count", 1)),
                    "bead_frame_count": int(pair_grp.attrs.get("cgl_bead_frame_count", 1)),
                }
                pair_compaction = _fit_pair_compaction_result_for_states(
                    result_cg=pair_result_cg,
                    compaction_states=compaction_states,
                    bead_types=list(dopc["bead_types"]),
                    bead_charges=list(dopc["bead_charges"]),
                    pair_params=pair_params,
                    lipids_itp_path=lipids_itp_path,
                    face_cos_min=face_cos_min,
                    radial_cutoff_nm=radial_cutoff_nm,
                    pair_relax_state_correction=pair_relax_state_correction,
                    mask_mode=mask_mode,
                    correction_center_mode=correction_center_mode,
                    pair_state_model=pair_state_model,
                )
            corrections["cg_lipid_compaction"].update(
                {
                    "delta_extended_extended": np.asarray(
                        pair_compaction["delta_extended_extended"], dtype=np.float32
                    ),
                    "delta_extended_compact": np.asarray(
                        pair_compaction["delta_extended_compact"], dtype=np.float32
                    ),
                    "delta_compact_compact": np.asarray(
                        pair_compaction["delta_compact_compact"], dtype=np.float32
                    ),
                    "pair_reference_compact_center_ang": float(compact_center_ang),
                    "pair_reference_extended_center_ang": float(extended_center_ang),
                    "pair_compaction_state_source": str(compaction_state_source_used),
                }
            )
            if "delta_extended_compressed" in pair_compaction:
                corrections["cg_lipid_compaction"]["delta_extended_compressed"] = np.asarray(
                    pair_compaction["delta_extended_compressed"],
                    dtype=np.float32,
                )
                corrections["cg_lipid_compaction"]["delta_compact_compressed"] = np.asarray(
                    pair_compaction["delta_compact_compressed"],
                    dtype=np.float32,
                )
                corrections["cg_lipid_compaction"]["delta_compressed_compressed"] = np.asarray(
                    pair_compaction["delta_compressed_compressed"],
                    dtype=np.float32,
                )
                corrections["cg_lipid_compaction"]["pair_reference_compressed_center_ang"] = float(
                    compressed_center_ang
                )
            for dataset_name in (
                "grid_extended_extended_kj_mol",
                "grid_extended_compact_kj_mol",
                "grid_compact_compact_kj_mol",
                "grid_extended_compressed_kj_mol",
                "grid_compact_compressed_kj_mol",
                "grid_compressed_compressed_kj_mol",
                "grid_average_kj_mol",
            ):
                if dataset_name in pair_compaction:
                    corrections["cg_lipid_compaction"][dataset_name] = np.asarray(
                        pair_compaction[dataset_name],
                        dtype=np.float32,
                    )
            if "correction_center_mode" in pair_compaction:
                corrections["cg_lipid_compaction"]["correction_center_mode"] = str(
                    pair_compaction["correction_center_mode"]
                )
            if "pair_state_model" in pair_compaction:
                corrections["cg_lipid_compaction"]["pair_state_model"] = str(
                    pair_compaction["pair_state_model"]
                )
            if "face_mask" in pair_compaction:
                corrections["cg_lipid_compaction"]["face_mask"] = np.asarray(
                    pair_compaction["face_mask"], dtype=np.float32
                )
                corrections["cg_lipid_compaction"]["face_cos_min"] = float(
                    pair_compaction["face_cos_min"]
                )
                corrections["cg_lipid_compaction"]["radial_cutoff_nm"] = float(
                    pair_compaction["radial_cutoff_nm"]
                )
                corrections["cg_lipid_compaction"]["mask_mode"] = str(
                    pair_compaction["mask_mode"]
                )

        sc_grp = cg_grp.get("cg_lipid_sc")
        if include_sc and sc_grp is not None:
            source_sc_compact_ang = float(sc_grp.attrs.get("reference_compact_center_ang", np.nan))
            source_sc_extended_ang = float(sc_grp.attrs.get("reference_extended_center_ang", np.nan))
            if (
                np.isfinite(source_sc_compact_ang)
                and np.isfinite(source_sc_extended_ang)
                and "delta_extended" in sc_grp
                and "delta_compact" in sc_grp
            ):
                corrections["cg_lipid_sc"] = {
                    "delta_extended": np.asarray(sc_grp["delta_extended"][:], dtype=np.float32),
                    "delta_compact": np.asarray(sc_grp["delta_compact"][:], dtype=np.float32),
                }
                if "grid_extended_kj_mol" in sc_grp:
                    corrections["cg_lipid_sc"]["grid_extended_kj_mol"] = np.asarray(
                        sc_grp["grid_extended_kj_mol"][:],
                        dtype=np.float32,
                    )
                if "grid_compact_kj_mol" in sc_grp:
                    corrections["cg_lipid_sc"]["grid_compact_kj_mol"] = np.asarray(
                        sc_grp["grid_compact_kj_mol"][:],
                        dtype=np.float32,
                    )
                if compressed_state_coord is not None:
                    restype_order = _decode_h5_string_array(sc_grp["restype_order"][:])
                    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
                    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
                    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
                        martinize_path, forcefield_name
                    )
                    cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
                    cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)
                    sc_n_modes = int(sc_grp.attrs.get("n_modes", 0))
                    if sc_n_modes != 0:
                        raise RuntimeError("single-CGL compaction retrofit only supports full-tensor cg_lipid_sc")
                    sc_n_radial = int(sc_grp.attrs["n_radial"])
                    sc_n_angular = int(sc_grp.attrs["n_angular"])
                    sc_r_count = sc_n_radial
                    sc_cos_theta_count = max(sc_n_angular - 2, 1)
                    sc_temperature = float(
                        sc_grp.attrs.get("boltzmann_temperature_upside", DEFAULT_PRODUCTION_TEMP_UPSIDE)
                    )
                    sc_average_temperature = float(
                        sc_grp.attrs.get(
                            "azimuthal_average_temperature_upside",
                            DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                        )
                    )
                    sc_angular_smooth = float(sc_grp.attrs.get("angular_smooth", 0.1))
                    sc_tasks = []
                    for ri, residue in enumerate(restype_order):
                        if residue not in residue_map or residue not in orientation_map:
                            raise RuntimeError(
                                f"{base_h5_path} cg_lipid_sc residue {residue!r} is missing orientation/source data"
                            )
                        sc_bead_types = residue_map[residue]
                        sc_bead_charges = [infer_charge_from_atomtype(bt) for bt in sc_bead_types]
                        orientation = orientation_map[residue]
                        sc_positions_by_rotamer = _expand_rotamer_sidechain_positions(
                            orientation,
                            residue,
                            np.asarray(martini_sidechain_offsets_nm[residue], dtype=np.float64),
                        )
                        sc_tasks.append(
                            {
                                "ri": ri,
                                "residue": residue,
                                "ref_bead_positions_nm": ref_ensemble_nm,
                                "cg_bead_types": list(dopc["bead_types"]),
                                "cg_bead_charges": list(dopc["bead_charges"]),
                                "target_type": "CGL",
                                "target_charge": 0.0,
                                "rotamer_bead_positions_nm": sc_positions_by_rotamer,
                                "rotamer_weights": orientation["weight"],
                                "sc_bead_types": list(sc_bead_types),
                                "sc_bead_charges": list(sc_bead_charges),
                                "pair_params": pair_params,
                                "cb_anchor_nm": cb_anchor_nm,
                                "cb_vector_unit": cb_vector_unit,
                                "r_min_nm": float(sc_grp.attrs.get("fit_r_min_nm", 0.30)),
                                "r_max_nm": float(sc_grp.attrs["fit_r_max_nm"]),
                                "r_count": sc_r_count,
                                "cos_theta_count": sc_cos_theta_count,
                                "azimuthal_count": int(sc_grp.attrs.get("azimuthal_count", 4)),
                                "sidechain_bead_frame_count": int(
                                    sc_grp.attrs.get("sidechain_bead_frame_count", 1)
                                ),
                                "cg_bead_frame_count": int(sc_grp.attrs.get("cgl_bead_frame_count", 1)),
                                "n_modes": sc_n_modes,
                                "n_knot_radial": sc_n_radial,
                                "n_knot_angular": sc_n_angular,
                                "angular_smooth": sc_angular_smooth,
                                "knot_spacing_ang": float(sc_grp.attrs["knot_spacing_ang"]),
                                "temperature": sc_temperature,
                                "average_temperature": sc_average_temperature,
                            }
                        )
                    compressed_sc_results = {
                        ri: result_sc
                        for ri, _residue, result_sc in _parallel_map_ordered(
                            "CG-SC compressed-state retrofit table",
                            _fit_cg_lipid_sc_quadspline_from_dict,
                            [
                                {
                                    **task,
                                    "ref_bead_positions_nm": np.asarray(
                                        compaction_states["compressed_refs_nm"],
                                        dtype=np.float64,
                                    ),
                                }
                                for task in sc_tasks
                            ],
                        )
                    }
                    n_sc_types = len(restype_order)
                    sc_n_param = sc_n_radial * sc_n_angular * sc_n_angular
                    sc_delta_compressed = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
                    sc_grid_compressed = np.zeros(
                        (n_sc_types, sc_r_count, sc_cos_theta_count, sc_cos_theta_count),
                        dtype=np.float32,
                    )
                    for ri in range(n_sc_types):
                        base_interaction_param = np.asarray(
                            sc_grp["interaction_param"][ri, 0, :],
                            dtype=np.float64,
                        )
                        base_reference_energy = float(
                            np.asarray(sc_grp["reference_energy_eup"][ri, 0], dtype=np.float64)
                        )
                        r_values_nm = np.linspace(
                            float(sc_grp.attrs.get("fit_r_min_nm", 0.30)),
                            float(sc_grp.attrs["fit_r_max_nm"]),
                            sc_r_count,
                            dtype=np.float64,
                        )
                        cos_theta_grid = np.linspace(-1.0, 1.0, sc_cos_theta_count, dtype=np.float64)
                        delta = _fit_single_cgl_state_delta_full_tensor(
                            base_energy_grid_kj_mol=np.asarray(
                                compressed_sc_results[ri]["energy_grid_raw"],
                                dtype=np.float64,
                            ),
                            extended_energy_grid_kj_mol=np.asarray(
                                compressed_sc_results[ri]["energy_grid_raw"],
                                dtype=np.float64,
                            ),
                            compact_energy_grid_kj_mol=np.asarray(
                                compressed_sc_results[ri]["energy_grid_raw"],
                                dtype=np.float64,
                            ),
                            compressed_energy_grid_kj_mol=np.asarray(
                                compressed_sc_results[ri]["energy_grid_raw"],
                                dtype=np.float64,
                            ),
                            reference_energy_eup=base_reference_energy,
                            base_interaction_param=base_interaction_param,
                            temperature_upside=float(sc_temperature),
                            r_values_nm=r_values_nm,
                            cos_theta_grid=cos_theta_grid,
                            n_knot_radial=sc_n_radial,
                            n_knot_angular=sc_n_angular,
                            knot_spacing_ang=float(sc_grp.attrs["knot_spacing_ang"]),
                            smooth=sc_angular_smooth,
                        )
                        sc_delta_compressed[ri, 0, :] = np.asarray(
                            delta["delta_compressed"],
                            dtype=np.float32,
                        )
                        sc_grid_compressed[ri] = np.asarray(
                            delta["grid_compressed_kj_mol"],
                            dtype=np.float32,
                        )
                    corrections["cg_lipid_sc"]["delta_compressed"] = sc_delta_compressed
                    corrections["cg_lipid_sc"]["grid_compressed_kj_mol"] = sc_grid_compressed
            else:
                restype_order = _decode_h5_string_array(sc_grp["restype_order"][:])
                orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
                residue_map = load_martini_forcefield(martinize_path, forcefield_name)
                martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
                    martinize_path, forcefield_name
                )
                cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
                cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)
                sc_n_modes = int(sc_grp.attrs.get("n_modes", 0))
                if sc_n_modes != 0:
                    raise RuntimeError("single-CGL compaction retrofit only supports full-tensor cg_lipid_sc")
                sc_n_radial = int(sc_grp.attrs["n_radial"])
                sc_n_angular = int(sc_grp.attrs["n_angular"])
                sc_r_count = sc_n_radial
                sc_cos_theta_count = max(sc_n_angular - 2, 1)
                sc_tasks = []
                sc_temperature = float(
                    sc_grp.attrs.get("boltzmann_temperature_upside", DEFAULT_PRODUCTION_TEMP_UPSIDE)
                )
                sc_average_temperature = float(
                    sc_grp.attrs.get(
                        "azimuthal_average_temperature_upside",
                        DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                    )
                )
                sc_angular_smooth = float(sc_grp.attrs.get("angular_smooth", 0.1))
                for ri, residue in enumerate(restype_order):
                    if residue not in residue_map or residue not in orientation_map:
                        raise RuntimeError(
                            f"{base_h5_path} cg_lipid_sc residue {residue!r} is missing orientation/source data"
                        )
                    sc_bead_types = residue_map[residue]
                    sc_bead_charges = [infer_charge_from_atomtype(bt) for bt in sc_bead_types]
                    orientation = orientation_map[residue]
                    sc_positions_by_rotamer = _expand_rotamer_sidechain_positions(
                        orientation,
                        residue,
                        np.asarray(martini_sidechain_offsets_nm[residue], dtype=np.float64),
                    )
                    sc_tasks.append(
                        {
                            "ri": ri,
                            "residue": residue,
                            "ref_bead_positions_nm": ref_ensemble_nm,
                            "cg_bead_types": list(dopc["bead_types"]),
                            "cg_bead_charges": list(dopc["bead_charges"]),
                            "target_type": "CGL",
                            "target_charge": 0.0,
                            "rotamer_bead_positions_nm": sc_positions_by_rotamer,
                            "rotamer_weights": orientation["weight"],
                            "sc_bead_types": list(sc_bead_types),
                            "sc_bead_charges": list(sc_bead_charges),
                            "pair_params": pair_params,
                            "cb_anchor_nm": cb_anchor_nm,
                            "cb_vector_unit": cb_vector_unit,
                            "r_min_nm": float(sc_grp.attrs.get("fit_r_min_nm", 0.30)),
                            "r_max_nm": float(sc_grp.attrs["fit_r_max_nm"]),
                            "r_count": sc_r_count,
                            "cos_theta_count": sc_cos_theta_count,
                            "azimuthal_count": int(sc_grp.attrs.get("azimuthal_count", 4)),
                            "sidechain_bead_frame_count": int(
                                sc_grp.attrs.get("sidechain_bead_frame_count", 1)
                            ),
                            "cg_bead_frame_count": int(sc_grp.attrs.get("cgl_bead_frame_count", 1)),
                            "n_modes": sc_n_modes,
                            "n_knot_radial": sc_n_radial,
                            "n_knot_angular": sc_n_angular,
                            "angular_smooth": sc_angular_smooth,
                            "knot_spacing_ang": float(sc_grp.attrs["knot_spacing_ang"]),
                            "temperature": sc_temperature,
                            "average_temperature": sc_average_temperature,
                        }
                    )

                base_sc_results = {
                    ri: result_sc
                    for ri, _residue, result_sc in _parallel_map_ordered(
                        "CG-SC base retrofit table",
                        _fit_cg_lipid_sc_quadspline_from_dict,
                        sc_tasks,
                    )
                }
                extended_sc_results = {
                    ri: result_sc
                    for ri, _residue, result_sc in _parallel_map_ordered(
                        "CG-SC extended-state retrofit table",
                        _fit_cg_lipid_sc_quadspline_from_dict,
                        [
                            {
                                **task,
                                "ref_bead_positions_nm": np.asarray(
                                    compaction_states["extended_refs_nm"],
                                    dtype=np.float64,
                                ),
                            }
                            for task in sc_tasks
                        ],
                    )
                }
                compact_sc_results = {
                    ri: result_sc
                    for ri, _residue, result_sc in _parallel_map_ordered(
                        "CG-SC compacted-state retrofit table",
                        _fit_cg_lipid_sc_quadspline_from_dict,
                        [
                            {
                                **task,
                                "ref_bead_positions_nm": np.asarray(
                                    compaction_states["compact_refs_nm"],
                                    dtype=np.float64,
                                ),
                            }
                            for task in sc_tasks
                        ],
                    )
                }
                compressed_sc_results = None
                if compressed_state_coord is not None:
                    compressed_sc_results = {
                        ri: result_sc
                        for ri, _residue, result_sc in _parallel_map_ordered(
                            "CG-SC compressed-state retrofit table",
                            _fit_cg_lipid_sc_quadspline_from_dict,
                            [
                                {
                                    **task,
                                    "ref_bead_positions_nm": np.asarray(
                                        compaction_states["compressed_refs_nm"],
                                        dtype=np.float64,
                                    ),
                                }
                                for task in sc_tasks
                            ],
                        )
                    }
                n_sc_types = len(restype_order)
                sc_n_param = sc_n_radial * sc_n_angular * sc_n_angular
                sc_delta_extended = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
                sc_delta_compact = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
                sc_delta_compressed = (
                    np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
                    if compressed_sc_results is not None
                    else None
                )
                sc_grid_extended = np.zeros(
                    (n_sc_types, sc_r_count, sc_cos_theta_count, sc_cos_theta_count),
                    dtype=np.float32,
                )
                sc_grid_compact = np.zeros_like(sc_grid_extended)
                sc_grid_compressed = (
                    np.zeros_like(sc_grid_extended)
                    if compressed_sc_results is not None
                    else None
                )
                for ri in range(n_sc_types):
                    base_sc = base_sc_results[ri]
                    delta = _fit_single_cgl_state_delta_full_tensor(
                        base_energy_grid_kj_mol=np.asarray(base_sc["energy_grid_raw"], dtype=np.float64),
                        extended_energy_grid_kj_mol=np.asarray(
                            extended_sc_results[ri]["energy_grid_raw"],
                            dtype=np.float64,
                        ),
                        compact_energy_grid_kj_mol=np.asarray(
                            compact_sc_results[ri]["energy_grid_raw"],
                            dtype=np.float64,
                        ),
                        compressed_energy_grid_kj_mol=(
                            np.asarray(
                                compressed_sc_results[ri]["energy_grid_raw"],
                                dtype=np.float64,
                            )
                            if compressed_sc_results is not None
                            else None
                        ),
                        reference_energy_eup=float(base_sc["reference_energy_eup"]),
                        base_interaction_param=np.asarray(base_sc["interaction_param"], dtype=np.float64),
                        temperature_upside=float(sc_temperature),
                        r_values_nm=np.asarray(base_sc["r_values_nm"], dtype=np.float64),
                        cos_theta_grid=np.asarray(base_sc["cos_theta_grid"], dtype=np.float64),
                        n_knot_radial=sc_n_radial,
                        n_knot_angular=sc_n_angular,
                        knot_spacing_ang=float(sc_grp.attrs["knot_spacing_ang"]),
                        smooth=sc_angular_smooth,
                    )
                    sc_delta_extended[ri, 0, :] = np.asarray(delta["delta_extended"], dtype=np.float32)
                    sc_delta_compact[ri, 0, :] = np.asarray(delta["delta_compact"], dtype=np.float32)
                    sc_grid_extended[ri] = np.asarray(delta["grid_extended_kj_mol"], dtype=np.float32)
                    sc_grid_compact[ri] = np.asarray(delta["grid_compact_kj_mol"], dtype=np.float32)
                    if sc_delta_compressed is not None and sc_grid_compressed is not None:
                        sc_delta_compressed[ri, 0, :] = np.asarray(
                            delta["delta_compressed"],
                            dtype=np.float32,
                        )
                        sc_grid_compressed[ri] = np.asarray(
                            delta["grid_compressed_kj_mol"],
                            dtype=np.float32,
                        )
                corrections["cg_lipid_sc"] = {
                    "delta_extended": sc_delta_extended,
                    "delta_compact": sc_delta_compact,
                    "grid_extended_kj_mol": sc_grid_extended,
                    "grid_compact_kj_mol": sc_grid_compact,
                }
                if sc_delta_compressed is not None and sc_grid_compressed is not None:
                    corrections["cg_lipid_sc"]["delta_compressed"] = sc_delta_compressed
                    corrections["cg_lipid_sc"]["grid_compressed_kj_mol"] = sc_grid_compressed

        target_grp = cg_grp.get("cg_lipid_target")
        if include_target and target_grp is not None:
            source_target_compact_ang = float(target_grp.attrs.get("reference_compact_center_ang", np.nan))
            source_target_extended_ang = float(target_grp.attrs.get("reference_extended_center_ang", np.nan))
            if (
                np.isfinite(source_target_compact_ang)
                and np.isfinite(source_target_extended_ang)
                and "delta_extended" in target_grp
                and "delta_compact" in target_grp
            ):
                corrections["cg_lipid_target"] = {
                    "delta_extended": np.asarray(target_grp["delta_extended"][:], dtype=np.float32),
                    "delta_compact": np.asarray(target_grp["delta_compact"][:], dtype=np.float32),
                }
                if "grid_extended_kj_mol" in target_grp:
                    corrections["cg_lipid_target"]["grid_extended_kj_mol"] = np.asarray(
                        target_grp["grid_extended_kj_mol"][:],
                        dtype=np.float32,
                    )
                if "grid_compact_kj_mol" in target_grp:
                    corrections["cg_lipid_target"]["grid_compact_kj_mol"] = np.asarray(
                        target_grp["grid_compact_kj_mol"][:],
                        dtype=np.float32,
                    )
                if compressed_state_coord is not None:
                    target_types = _decode_h5_string_array(target_grp["target_order"][:])
                    n_types = len(target_types)
                    n_knot_radial = int(target_grp.attrs["n_radial"])
                    n_knot_angular = int(target_grp.attrs["n_angular"])
                    knot_spacing_ang = float(target_grp.attrs["knot_spacing_ang"])
                    target_temperature = float(
                        target_grp.attrs.get("boltzmann_temperature_upside", DEFAULT_PRODUCTION_TEMP_UPSIDE)
                    )
                    target_average_temperature = float(
                        target_grp.attrs.get(
                            "azimuthal_average_temperature_upside",
                            DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                        )
                    )
                    r_max_ang = float((n_knot_radial - 2) * knot_spacing_ang)
                    r_sample_ang = np.linspace(float(knot_spacing_ang), r_max_ang, n_knot_radial)
                    r_sample_nm = r_sample_ang / LENGTH_CONVERSION_A_PER_NM
                    t_angular_sample = np.linspace(1.0, float(n_knot_angular - 2), n_knot_angular)
                    cos_theta_grid = (t_angular_sample - 1.0) / (0.5 * float(n_knot_angular - 3)) - 1.0
                    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
                    phi_values = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)
                    bead_frame_angles = _bead_frame_angles(
                        int(target_grp.attrs.get("cgl_bead_frame_count", 1))
                    )
                    target_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
                    orientation_dirs = _directions_with_dot_np(target_axis, cos_theta_grid, phi_values)
                    target_tasks = [
                        {
                            "ti": ti,
                            "target_type": tgt_type,
                            "bead_types": list(dopc["bead_types"]),
                            "bead_charges": bead_charges,
                            "pair_params": pair_params,
                            "ref_nm": ref_ensemble_nm,
                            "r_sample_nm": r_sample_nm,
                            "cos_theta_grid": cos_theta_grid,
                            "orientation_dirs": orientation_dirs,
                            "bead_frame_angles": bead_frame_angles,
                            "temperature": target_temperature,
                            "average_temperature": target_average_temperature,
                        }
                        for ti, tgt_type in enumerate(target_types)
                    ]
                    compressed_target_results = {
                        ti: np.asarray(energy_grid, dtype=np.float64)
                        for ti, _tgt_type, _control_flat, _ref_kj_mol, energy_grid in _parallel_map_ordered(
                            "CGL-particle compressed-state retrofit table",
                            _run_cgl_target_type_task,
                            [
                                {
                                    **task,
                                    "ref_nm": np.asarray(
                                        compaction_states["compressed_refs_nm"],
                                        dtype=np.float64,
                                    ),
                                }
                                for task in target_tasks
                            ],
                        )
                    }
                    target_ref_energy = np.asarray(
                        target_grp["reference_energy_eup"][:],
                        dtype=np.float64,
                    ).reshape(-1)
                    n_param = n_knot_radial * n_knot_angular
                    target_delta_compressed = np.zeros((1, n_types, n_param), dtype=np.float32)
                    target_grid_compressed = np.zeros(
                        (n_types, len(r_sample_nm), len(cos_theta_grid)),
                        dtype=np.float32,
                    )
                    for ti in range(n_types):
                        base_interaction_param = np.asarray(
                            target_grp["interaction_param"][0, ti, :],
                            dtype=np.float64,
                        )
                        delta = _build_single_cgl_state_delta_radial_angular(
                            base_energy_grid_kj_mol=compressed_target_results[ti],
                            extended_energy_grid_kj_mol=compressed_target_results[ti],
                            compact_energy_grid_kj_mol=compressed_target_results[ti],
                            compressed_energy_grid_kj_mol=compressed_target_results[ti],
                            reference_energy_eup=float(target_ref_energy[ti]),
                            temperature_upside=float(target_temperature),
                            r_values_nm=r_sample_nm,
                            cos_theta_grid=cos_theta_grid,
                            n_knot_radial=n_knot_radial,
                            n_knot_angular=n_knot_angular,
                            knot_spacing_ang=knot_spacing_ang,
                            base_interaction_param=base_interaction_param,
                        )
                        target_delta_compressed[0, ti, :] = np.asarray(
                            delta["delta_compressed"],
                            dtype=np.float32,
                        )
                        target_grid_compressed[ti] = np.asarray(
                            delta["grid_compressed_kj_mol"],
                            dtype=np.float32,
                        )
                    corrections["cg_lipid_target"]["delta_compressed"] = target_delta_compressed
                    corrections["cg_lipid_target"]["grid_compressed_kj_mol"] = target_grid_compressed
            else:
                target_types = _decode_h5_string_array(target_grp["target_order"][:])
                n_types = len(target_types)
                n_knot_radial = int(target_grp.attrs["n_radial"])
                n_knot_angular = int(target_grp.attrs["n_angular"])
                knot_spacing_ang = float(target_grp.attrs["knot_spacing_ang"])
                target_temperature = float(
                    target_grp.attrs.get("boltzmann_temperature_upside", DEFAULT_PRODUCTION_TEMP_UPSIDE)
                )
                target_average_temperature = float(
                    target_grp.attrs.get(
                        "azimuthal_average_temperature_upside",
                        DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                    )
                )
                r_max_ang = float((n_knot_radial - 2) * knot_spacing_ang)
                r_sample_ang = np.linspace(float(knot_spacing_ang), r_max_ang, n_knot_radial)
                r_sample_nm = r_sample_ang / LENGTH_CONVERSION_A_PER_NM
                t_angular_sample = np.linspace(1.0, float(n_knot_angular - 2), n_knot_angular)
                cos_theta_grid = (t_angular_sample - 1.0) / (0.5 * float(n_knot_angular - 3)) - 1.0
                cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
                phi_values = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)
                bead_frame_angles = _bead_frame_angles(
                    int(target_grp.attrs.get("cgl_bead_frame_count", 1))
                )
                target_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
                orientation_dirs = _directions_with_dot_np(target_axis, cos_theta_grid, phi_values)
                target_tasks = [
                    {
                        "ti": ti,
                        "target_type": tgt_type,
                        "bead_types": list(dopc["bead_types"]),
                        "bead_charges": bead_charges,
                        "pair_params": pair_params,
                        "ref_nm": ref_ensemble_nm,
                        "r_sample_nm": r_sample_nm,
                        "cos_theta_grid": cos_theta_grid,
                        "orientation_dirs": orientation_dirs,
                        "bead_frame_angles": bead_frame_angles,
                        "temperature": target_temperature,
                        "average_temperature": target_average_temperature,
                    }
                    for ti, tgt_type in enumerate(target_types)
                ]
                base_target_results = {
                    ti: np.asarray(energy_grid, dtype=np.float64)
                    for ti, _tgt_type, _control_flat, _ref_kj_mol, energy_grid in _parallel_map_ordered(
                        "CGL-particle base retrofit table",
                        _run_cgl_target_type_task,
                        target_tasks,
                    )
                }
                extended_target_results = {
                    ti: np.asarray(energy_grid, dtype=np.float64)
                    for ti, _tgt_type, _control_flat, _ref_kj_mol, energy_grid in _parallel_map_ordered(
                        "CGL-particle extended-state retrofit table",
                        _run_cgl_target_type_task,
                        [
                            {
                                **task,
                                "ref_nm": np.asarray(compaction_states["extended_refs_nm"], dtype=np.float64),
                            }
                            for task in target_tasks
                        ],
                    )
                }
                compact_target_results = {
                    ti: np.asarray(energy_grid, dtype=np.float64)
                    for ti, _tgt_type, _control_flat, _ref_kj_mol, energy_grid in _parallel_map_ordered(
                        "CGL-particle compacted-state retrofit table",
                        _run_cgl_target_type_task,
                        [
                            {
                                **task,
                                "ref_nm": np.asarray(compaction_states["compact_refs_nm"], dtype=np.float64),
                            }
                            for task in target_tasks
                        ],
                    )
                }
                compressed_target_results = None
                if compressed_state_coord is not None:
                    compressed_target_results = {
                        ti: np.asarray(energy_grid, dtype=np.float64)
                        for ti, _tgt_type, _control_flat, _ref_kj_mol, energy_grid in _parallel_map_ordered(
                            "CGL-particle compressed-state retrofit table",
                            _run_cgl_target_type_task,
                            [
                                {
                                    **task,
                                    "ref_nm": np.asarray(
                                        compaction_states["compressed_refs_nm"],
                                        dtype=np.float64,
                                    ),
                                }
                                for task in target_tasks
                            ],
                        )
                    }
                target_ref_energy = np.asarray(
                    target_grp["reference_energy_eup"][:],
                    dtype=np.float64,
                ).reshape(-1)
                n_param = n_knot_radial * n_knot_angular
                target_delta_extended = np.zeros((1, n_types, n_param), dtype=np.float32)
                target_delta_compact = np.zeros((1, n_types, n_param), dtype=np.float32)
                target_delta_compressed = (
                    np.zeros((1, n_types, n_param), dtype=np.float32)
                    if compressed_target_results is not None
                    else None
                )
                target_grid_extended = np.zeros(
                    (n_types, len(r_sample_nm), len(cos_theta_grid)),
                    dtype=np.float32,
                )
                target_grid_compact = np.zeros_like(target_grid_extended)
                target_grid_compressed = (
                    np.zeros_like(target_grid_extended)
                    if compressed_target_results is not None
                    else None
                )
                for ti in range(n_types):
                    delta = _build_single_cgl_state_delta_radial_angular(
                        base_energy_grid_kj_mol=base_target_results[ti],
                        extended_energy_grid_kj_mol=extended_target_results[ti],
                        compact_energy_grid_kj_mol=compact_target_results[ti],
                        compressed_energy_grid_kj_mol=(
                            compressed_target_results[ti]
                            if compressed_target_results is not None
                            else None
                        ),
                        reference_energy_eup=float(target_ref_energy[ti]),
                        temperature_upside=float(target_temperature),
                        r_values_nm=r_sample_nm,
                        cos_theta_grid=cos_theta_grid,
                        n_knot_radial=n_knot_radial,
                        n_knot_angular=n_knot_angular,
                        knot_spacing_ang=knot_spacing_ang,
                    )
                    target_delta_extended[0, ti, :] = np.asarray(
                        delta["delta_extended"],
                        dtype=np.float32,
                    )
                    target_delta_compact[0, ti, :] = np.asarray(
                        delta["delta_compact"],
                        dtype=np.float32,
                    )
                    target_grid_extended[ti] = np.asarray(
                        delta["grid_extended_kj_mol"],
                        dtype=np.float32,
                    )
                    target_grid_compact[ti] = np.asarray(
                        delta["grid_compact_kj_mol"],
                        dtype=np.float32,
                    )
                    if target_delta_compressed is not None and target_grid_compressed is not None:
                        target_delta_compressed[0, ti, :] = np.asarray(
                            delta["delta_compressed"],
                            dtype=np.float32,
                        )
                        target_grid_compressed[ti] = np.asarray(
                            delta["grid_compressed_kj_mol"],
                            dtype=np.float32,
                        )
                corrections["cg_lipid_target"] = {
                    "delta_extended": target_delta_extended,
                    "delta_compact": target_delta_compact,
                    "grid_extended_kj_mol": target_grid_extended,
                    "grid_compact_kj_mol": target_grid_compact,
                }
                if target_delta_compressed is not None and target_grid_compressed is not None:
                    corrections["cg_lipid_target"]["delta_compressed"] = target_delta_compressed
                    corrections["cg_lipid_target"]["grid_compressed_kj_mol"] = target_grid_compressed

    print(
        "  Single-CGL compaction retrofit: "
        f"ref={ref_dataset}, pair={'delta_compact_compact' in corrections['cg_lipid_compaction']}, "
        f"SC={'cg_lipid_sc' in corrections}, "
        f"target={'cg_lipid_target' in corrections}, "
        f"state source={compaction_state_source_used}, "
        f"reference centers=({extended_center_ang:.3f}, {compact_center_ang:.3f}) A"
        + (
            f", compressed={float(compressed_center_ang):.3f} A"
            if compressed_center_ang is not None and compressed_state_coord is not None
            else ""
        )
    )
    return corrections


def _apply_single_cgl_compaction_corrections_to_h5(cg_grp: h5py.Group, corrections: dict) -> None:
    if "cg_lipid_compaction" in cg_grp and "cg_lipid_compaction" in corrections:
        comp_grp = cg_grp["cg_lipid_compaction"]
        payload = corrections["cg_lipid_compaction"]
        for dataset_name in (
            "self_coeff",
            "pmf_centers_ang",
            "pmf_values_kj_mol",
        ):
            if dataset_name in comp_grp:
                del comp_grp[dataset_name]
        pair_payload_names = (
            "delta_extended_extended",
            "delta_extended_compact",
            "delta_compact_compact",
            "delta_extended_compressed",
            "delta_compact_compressed",
            "delta_compressed_compressed",
            "grid_extended_extended_kj_mol",
            "grid_extended_compact_kj_mol",
            "grid_compact_compact_kj_mol",
            "grid_extended_compressed_kj_mol",
            "grid_compact_compressed_kj_mol",
            "grid_compressed_compressed_kj_mol",
            "grid_average_kj_mol",
            "face_mask",
        )
        if any(dataset_name in payload for dataset_name in pair_payload_names):
            for dataset_name in pair_payload_names:
                if dataset_name in comp_grp:
                    del comp_grp[dataset_name]
        comp_grp.attrs["compact_state_center_ang"] = np.float32(corrections["compact_state_center_ang"])
        comp_grp.attrs["extended_state_center_ang"] = np.float32(corrections["extended_state_center_ang"])
        comp_grp.attrs["compact_state_probability"] = np.float32(corrections["compact_state_probability"])
        comp_grp.attrs["reference_compact_center_ang"] = np.float32(
            corrections["reference_compact_center_ang"]
        )
        comp_grp.attrs["reference_extended_center_ang"] = np.float32(
            corrections["reference_extended_center_ang"]
        )
        if "compressed_state_center_ang" in corrections:
            comp_grp.attrs["compressed_state_center_ang"] = np.float32(
                corrections["compressed_state_center_ang"]
            )
        if "reference_compressed_center_ang" in corrections:
            comp_grp.attrs["reference_compressed_center_ang"] = np.float32(
                corrections["reference_compressed_center_ang"]
            )
        comp_grp.attrs["self_coord_min_ang"] = np.float32(payload["coord_min_ang"])
        comp_grp.attrs["self_coord_max_ang"] = np.float32(payload["coord_max_ang"])
        comp_grp.attrs["self_coord_spacing_ang"] = np.float32(payload["coord_spacing_ang"])
        comp_grp.attrs["self_n_knot"] = np.int32(payload["n_knot"])
        comp_grp.attrs["effective_stiffness_eup_a2"] = np.float32(
            payload["effective_stiffness_eup_a2"]
        )
        comp_grp.attrs["mass_up"] = np.float32(payload["mass_up"])
        if "correction_center_mode" in payload:
            comp_grp.attrs["correction_center_mode"] = str(payload["correction_center_mode"])
        if "pair_state_model" in payload:
            comp_grp.attrs["pair_state_model"] = str(payload["pair_state_model"])
        if "pair_reference_compact_center_ang" in payload:
            comp_grp.attrs["pair_reference_compact_center_ang"] = np.float32(
                payload["pair_reference_compact_center_ang"]
            )
        if "pair_reference_extended_center_ang" in payload:
            comp_grp.attrs["pair_reference_extended_center_ang"] = np.float32(
                payload["pair_reference_extended_center_ang"]
            )
        if "pair_reference_compressed_center_ang" in payload:
            comp_grp.attrs["pair_reference_compressed_center_ang"] = np.float32(
                payload["pair_reference_compressed_center_ang"]
            )
        if "pair_compaction_state_source" in payload:
            comp_grp.attrs["pair_compaction_state_source"] = np.bytes_(
                str(payload["pair_compaction_state_source"])
            )
        if "face_cos_min" in payload:
            comp_grp.attrs["face_cos_min"] = np.float32(payload["face_cos_min"])
            comp_grp.attrs["radial_cutoff_nm"] = np.float32(payload["radial_cutoff_nm"])
            comp_grp.attrs["mask_mode"] = str(payload["mask_mode"])
        if "compaction_state_source" in corrections:
            comp_grp.attrs["compaction_state_source"] = np.bytes_(str(corrections["compaction_state_source"]))
        if "compaction_reference_dataset" in corrections:
            comp_grp.attrs["compaction_reference_dataset"] = np.bytes_(
                str(corrections["compaction_reference_dataset"])
            )
        if "self_compaction_sample_source" in corrections:
            comp_grp.attrs["self_compaction_sample_source"] = np.bytes_(
                str(corrections["self_compaction_sample_source"])
            )
        comp_grp.create_dataset("self_coeff", data=np.asarray(payload["self_coeff"], dtype=np.float32))
        comp_grp.create_dataset(
            "pmf_centers_ang",
            data=np.asarray(payload["pmf_centers_ang"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "pmf_values_kj_mol",
            data=np.asarray(payload["pmf_values_kj_mol"], dtype=np.float32),
        )
        for dataset_name in (
            "delta_extended_extended",
            "delta_extended_compact",
            "delta_compact_compact",
            "delta_extended_compressed",
            "delta_compact_compressed",
            "delta_compressed_compressed",
            "grid_extended_extended_kj_mol",
            "grid_extended_compact_kj_mol",
            "grid_compact_compact_kj_mol",
            "grid_extended_compressed_kj_mol",
            "grid_compact_compressed_kj_mol",
            "grid_compressed_compressed_kj_mol",
            "grid_average_kj_mol",
        ):
            if dataset_name in payload:
                comp_grp.create_dataset(
                    dataset_name,
                    data=np.asarray(payload[dataset_name], dtype=np.float32),
                )
        if "face_mask" in payload:
            comp_grp.create_dataset(
                "face_mask",
                data=np.asarray(payload["face_mask"], dtype=np.float32),
            )

    for group_name in ("cg_lipid_sc", "cg_lipid_target"):
        if group_name not in corrections:
            continue
        grp = cg_grp[group_name]
        for dataset_name in (
            "delta_extended",
            "delta_compact",
            "delta_compressed",
            "grid_extended_kj_mol",
            "grid_compact_kj_mol",
            "grid_compressed_kj_mol",
        ):
            if dataset_name in grp:
                del grp[dataset_name]
        grp.attrs["correction_layer"] = "source_derived_cgl_compaction_state"
        grp.attrs["compact_state_center_ang"] = np.float32(corrections["compact_state_center_ang"])
        grp.attrs["extended_state_center_ang"] = np.float32(corrections["extended_state_center_ang"])
        grp.attrs["compact_state_probability"] = np.float32(corrections["compact_state_probability"])
        grp.attrs["reference_compact_center_ang"] = np.float32(
            corrections["reference_compact_center_ang"]
        )
        grp.attrs["reference_extended_center_ang"] = np.float32(
            corrections["reference_extended_center_ang"]
        )
        if "compressed_state_center_ang" in corrections:
            grp.attrs["compressed_state_center_ang"] = np.float32(
                corrections["compressed_state_center_ang"]
            )
        if "reference_compressed_center_ang" in corrections:
            grp.attrs["reference_compressed_center_ang"] = np.float32(
                corrections["reference_compressed_center_ang"]
            )
        if "compaction_state_source" in corrections:
            grp.attrs["compaction_state_source"] = np.bytes_(str(corrections["compaction_state_source"]))
        payload = corrections[group_name]
        grp.create_dataset("delta_extended", data=np.asarray(payload["delta_extended"], dtype=np.float32))
        grp.create_dataset("delta_compact", data=np.asarray(payload["delta_compact"], dtype=np.float32))
        if "delta_compressed" in payload:
            grp.create_dataset(
                "delta_compressed",
                data=np.asarray(payload["delta_compressed"], dtype=np.float32),
            )
        grp.create_dataset(
            "grid_extended_kj_mol",
            data=np.asarray(payload["grid_extended_kj_mol"], dtype=np.float32),
        )
        grp.create_dataset(
            "grid_compact_kj_mol",
            data=np.asarray(payload["grid_compact_kj_mol"], dtype=np.float32),
        )
        if "grid_compressed_kj_mol" in payload:
            grp.create_dataset(
                "grid_compressed_kj_mol",
                data=np.asarray(payload["grid_compressed_kj_mol"], dtype=np.float32),
            )


def _load_cg_lipid_pair_energy_grid_from_h5(
    table_h5_path: Path,
    expected_shape: tuple[int, int, int],
    r_values_nm: np.ndarray,
) -> np.ndarray:
    table_h5_path = Path(table_h5_path).expanduser().resolve()
    with h5py.File(table_h5_path, "r") as h5:
        if "/cg_lipid_table/cg_lipid_pair" in h5:
            pair_grp = h5["/cg_lipid_table/cg_lipid_pair"]
        elif "/input/potential/cg_lipid_pair" in h5:
            pair_grp = h5["/input/potential/cg_lipid_pair"]
        else:
            raise RuntimeError(f"{table_h5_path} lacks a CGL pair table")
        if "energy_grid_raw_kj_mol" not in pair_grp:
            raise RuntimeError(f"{table_h5_path} lacks cg_lipid_pair/energy_grid_raw_kj_mol")
        energy_grid = np.asarray(pair_grp["energy_grid_raw_kj_mol"][:], dtype=np.float64)
        if energy_grid.shape != expected_shape:
            raise RuntimeError(
                f"{table_h5_path} energy grid shape {energy_grid.shape} "
                f"does not match expected {expected_shape}"
            )
        fit_r_min_nm = float(pair_grp.attrs.get("fit_r_min_nm", r_values_nm[0]))
        fit_r_max_nm = float(pair_grp.attrs.get("fit_r_max_nm", r_values_nm[-1]))
        if (
            abs(fit_r_min_nm - float(r_values_nm[0])) > 1.0e-6
            or abs(fit_r_max_nm - float(r_values_nm[-1])) > 1.0e-6
        ):
            raise RuntimeError(
                f"{table_h5_path} radial support {fit_r_min_nm:.6f}-{fit_r_max_nm:.6f} nm "
                f"does not match expected {float(r_values_nm[0]):.6f}-{float(r_values_nm[-1]):.6f} nm"
            )
    print(f"  CGL-CGL: loaded previous pair energy grid from {table_h5_path}")
    return energy_grid


def _run_cg_pair_tensor_task(task: dict) -> tuple[int, int, int, float]:
    ir = int(task["ir"])
    ia1 = int(task["ia1"])
    ia2 = int(task["ia2"])
    r_nm = float(task["r_nm"])
    dirs1 = task["dirs1"]
    dirs2 = task["dirs2"]
    bead_frame_angles = task["bead_frame_angles"]
    ref_nm1 = np.asarray(task.get("ref_nm1", task.get("ref_nm")), dtype=np.float64)
    if ref_nm1.ndim == 2:
        ref_nm1 = ref_nm1[None, :, :]
    ref_nm2 = np.asarray(task.get("ref_nm2", ref_nm1), dtype=np.float64)
    if ref_nm2.ndim == 2:
        ref_nm2 = ref_nm2[None, :, :]
    conformer_pairs = task.get("conformer_pairs")
    if conformer_pairs is None:
        conformer_pairs = []
    temperature = float(task.get("temperature", 0.0))
    all_pairs = _cg_lipid_conformer_pair_indices(ref_nm1.shape[0])
    if conformer_pairs:
        sample_values = []
        for dir1 in dirs1:
            for dir2 in dirs2:
                for frame_angle1 in bead_frame_angles:
                    for frame_angle2 in bead_frame_angles:
                        for i_conf, j_conf in conformer_pairs:
                            sample_values.append(
                                _compute_cg_pair_energy(
                                    r_nm,
                                    np.asarray(dir1, dtype=np.float64),
                                    np.asarray(dir2, dtype=np.float64),
                                    float(frame_angle1),
                                    float(frame_angle2),
                                    ref_nm1[int(i_conf)],
                                    ref_nm2[int(j_conf)],
                                    task.get("bead_types1", task["bead_types"]),
                                    task.get("bead_types2", task["bead_types"]),
                                    task.get("bead_charges1", task["bead_charges"]),
                                    task.get("bead_charges2", task["bead_charges"]),
                                    task["pair_params"],
                                    dist_min_nm=float(task["dist_min_nm"]),
                                )
                            )
        values = np.asarray(sample_values, dtype=np.float64)
    else:
        values = _compute_cg_pair_energy_samples_vectorized(
            r_nm,
            dirs1,
            dirs2,
            bead_frame_angles,
            ref_nm1,
            task.get("bead_types1", task["bead_types"]),
            task.get("bead_charges1", task["bead_charges"]),
            task["pair_params"],
            dist_min_nm=float(task["dist_min_nm"]),
            ref_nm2=ref_nm2,
            bead_types2=task.get("bead_types2", task["bead_types"]),
            bead_charges2=task.get("bead_charges2", task["bead_charges"]),
        )
    return ir, ia1, ia2, _boltzmann_free_energy_kj_mol(values, temperature)


def _sample_cg_pair_energy_grid(
    ref_bead_positions_nm1: np.ndarray,
    ref_bead_positions_nm2: np.ndarray,
    bead_types: list,
    bead_charges: list,
    pair_params: dict,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    azimuthal_count: int,
    bead_frame_count: int,
    dist_min_nm: float,
    average_temperature: float,
) -> np.ndarray:
    r_values = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta = np.asarray(cos_theta_grid, dtype=np.float64)
    phi_values = np.linspace(0.0, 2.0 * np.pi, int(azimuthal_count), endpoint=False, dtype=np.float64)
    bead_frame_angles = _bead_frame_angles(bead_frame_count)
    n12_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    dirs1 = _directions_with_dot_np(-n12_axis, cos_theta, phi_values)
    dirs2 = _directions_with_dot_np(n12_axis, cos_theta, phi_values)

    ref1 = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm1)
    ref2 = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm2)
    tasks = []
    for ir, r_nm in enumerate(r_values):
        for ia1 in range(cos_theta.size):
            for ia2 in range(cos_theta.size):
                tasks.append(
                    {
                        "ir": ir,
                        "ia1": ia1,
                        "ia2": ia2,
                        "r_nm": float(r_nm),
                        "dirs1": dirs1[ia1],
                        "dirs2": dirs2[ia2],
                        "bead_frame_angles": bead_frame_angles,
                        "ref_nm1": ref1,
                        "ref_nm2": ref2,
                        "bead_types": list(bead_types),
                        "bead_charges": list(bead_charges),
                        "pair_params": pair_params,
                        "dist_min_nm": float(dist_min_nm),
                        "temperature": float(average_temperature),
                    }
                )
    energy_grid = np.zeros((r_values.size, cos_theta.size, cos_theta.size), dtype=np.float64)
    for ir, ia1, ia2, energy in _parallel_map_ordered("CG-CG compaction table", _run_cg_pair_tensor_task, tasks):
        energy_grid[ir, ia1, ia2] = energy
    return energy_grid


def _extract_dopc_beads_centers_orientations_from_frame(
    frame_pos_ang: np.ndarray,
    lipid_molecules: list[np.ndarray],
    box_ang: np.ndarray | None,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
    bead_frames = []
    centers = []
    directions = []
    for lipid in lipid_molecules:
        beads_ang = np.asarray(frame_pos_ang[lipid, :], dtype=np.float64)
        delta = beads_ang - beads_ang[0][None, :]
        if box_ang is not None:
            delta -= box_ang[None, :] * np.round(delta / box_ang[None, :])
        beads_ang = beads_ang[0][None, :] + delta
        center = np.mean(beads_ang, axis=0)
        direction = ((beads_ang[8] + beads_ang[13]) * 0.5) - beads_ang[0]
        norm = float(np.linalg.norm(direction))
        if norm <= 1.0e-12:
            continue
        bead_frames.append(beads_ang)
        centers.append(center)
        directions.append(direction / norm)
    if not bead_frames:
        raise RuntimeError("No valid DOPC bead frames extracted from frame")
    return bead_frames, np.asarray(centers, dtype=np.float64), np.asarray(directions, dtype=np.float64)


def _least_squares_generalized_force_energy_grid(
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    base_energy_grid_kj_mol: np.ndarray,
    sampled: np.ndarray,
    counts: np.ndarray,
    mean_dudr: np.ndarray,
    mean_duda1: np.ndarray,
    mean_duda2: np.ndarray,
    min_count: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    var_index = np.full(sampled.shape, -1, dtype=np.int32)
    sampled_indices = np.argwhere(sampled)
    for idx, (ir, ia1, ia2) in enumerate(sampled_indices):
        var_index[int(ir), int(ia1), int(ia2)] = int(idx)
    n_var = int(sampled_indices.shape[0])
    if n_var == 0:
        return (
            np.asarray(base_energy_grid_kj_mol, dtype=np.float64).copy(),
            np.zeros_like(sampled, dtype=bool),
            {
                "lsq_variable_count": 0,
                "lsq_equation_count": 0,
                "lsq_residual_rms_kj_mol": 0.0,
                "lsq_solver": "none_no_sampled_bins",
            },
        )

    rows = []
    cols = []
    data = []
    rhs = []
    graph = [[] for _ in range(n_var)]

    def add_difference_equation(lo_key, hi_key, delta_coord: float, deriv_value: float, weight: float) -> None:
        row = len(rhs)
        lo = int(var_index[lo_key])
        hi = int(var_index[hi_key])
        rows.extend([row, row])
        cols.extend([lo, hi])
        data.extend([-weight, weight])
        rhs.append(weight * float(deriv_value) * float(delta_coord))
        graph[lo].append(hi)
        graph[hi].append(lo)

    count_scale = math.sqrt(float(max(1, int(min_count))))
    n_r, n_a1, n_a2 = sampled.shape
    for ir in range(n_r):
        for ia1 in range(n_a1):
            for ia2 in range(n_a2):
                if not sampled[ir, ia1, ia2]:
                    continue
                if ir + 1 < n_r and sampled[ir + 1, ia1, ia2]:
                    deriv = 0.5 * (float(mean_dudr[ir, ia1, ia2]) + float(mean_dudr[ir + 1, ia1, ia2]))
                    if np.isfinite(deriv):
                        weight = math.sqrt(min(float(counts[ir, ia1, ia2]), float(counts[ir + 1, ia1, ia2]))) / count_scale
                        add_difference_equation(
                            (ir, ia1, ia2),
                            (ir + 1, ia1, ia2),
                            float(r_values_nm[ir + 1] - r_values_nm[ir]),
                            deriv,
                            weight,
                        )
                if ia1 + 1 < n_a1 and sampled[ir, ia1 + 1, ia2]:
                    deriv = 0.5 * (float(mean_duda1[ir, ia1, ia2]) + float(mean_duda1[ir, ia1 + 1, ia2]))
                    if np.isfinite(deriv):
                        weight = math.sqrt(min(float(counts[ir, ia1, ia2]), float(counts[ir, ia1 + 1, ia2]))) / count_scale
                        add_difference_equation(
                            (ir, ia1, ia2),
                            (ir, ia1 + 1, ia2),
                            float(cos_theta_grid[ia1 + 1] - cos_theta_grid[ia1]),
                            deriv,
                            weight,
                        )
                if ia2 + 1 < n_a2 and sampled[ir, ia1, ia2 + 1]:
                    deriv = 0.5 * (float(mean_duda2[ir, ia1, ia2]) + float(mean_duda2[ir, ia1, ia2 + 1]))
                    if np.isfinite(deriv):
                        weight = math.sqrt(min(float(counts[ir, ia1, ia2]), float(counts[ir, ia1, ia2 + 1]))) / count_scale
                        add_difference_equation(
                            (ir, ia1, ia2),
                            (ir, ia1, ia2 + 1),
                            float(cos_theta_grid[ia2 + 1] - cos_theta_grid[ia2]),
                            deriv,
                            weight,
                        )

    seen = np.zeros(n_var, dtype=bool)
    for start in range(n_var):
        if seen[start]:
            continue
        stack = [start]
        component = []
        seen[start] = True
        while stack:
            node = stack.pop()
            component.append(node)
            for nb in graph[node]:
                if not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        anchor = max(
            component,
            key=lambda idx: (
                int(sampled_indices[idx, 0]),
                float(counts[tuple(sampled_indices[idx])]),
            ),
        )
        key = tuple(int(x) for x in sampled_indices[anchor])
        weight = math.sqrt(float(counts[key])) / count_scale
        row = len(rhs)
        rows.append(row)
        cols.append(anchor)
        data.append(weight)
        rhs.append(weight * float(base_energy_grid_kj_mol[key]))

    if not rhs:
        return (
            np.asarray(base_energy_grid_kj_mol, dtype=np.float64).copy(),
            np.zeros_like(sampled, dtype=bool),
            {
                "lsq_variable_count": n_var,
                "lsq_equation_count": 0,
                "lsq_residual_rms_kj_mol": 0.0,
                "lsq_solver": "none_no_equations",
            },
        )

    b = np.asarray(rhs, dtype=np.float64)
    solver = "scipy.sparse.linalg.lsqr"
    try:
        from scipy import sparse
        from scipy.sparse import linalg as sparse_linalg

        matrix = sparse.coo_matrix(
            (np.asarray(data, dtype=np.float64), (np.asarray(rows), np.asarray(cols))),
            shape=(len(rhs), n_var),
        ).tocsr()
        solution = sparse_linalg.lsqr(matrix, b, atol=1.0e-8, btol=1.0e-8, iter_lim=10000)[0]
        residual = matrix @ solution - b
    except Exception:
        solver = "numpy.linalg.lstsq"
        matrix = np.zeros((len(rhs), n_var), dtype=np.float64)
        for row, col, value in zip(rows, cols, data):
            matrix[int(row), int(col)] += float(value)
        solution = np.linalg.lstsq(matrix, b, rcond=None)[0]
        residual = matrix @ solution - b

    energy = np.asarray(base_energy_grid_kj_mol, dtype=np.float64).copy()
    updated = np.zeros_like(sampled, dtype=bool)
    for idx, key in enumerate(sampled_indices):
        ikey = tuple(int(x) for x in key)
        energy[ikey] = float(solution[idx])
        updated[ikey] = True
    rms = float(np.sqrt(np.mean(np.asarray(residual, dtype=np.float64) ** 2))) if residual.size else 0.0
    return energy, updated, {
        "lsq_variable_count": n_var,
        "lsq_equation_count": int(len(rhs)),
        "lsq_residual_rms_kj_mol": rms,
        "lsq_solver": solver,
    }


def _force_matched_pair_energy_from_upside_h5(
    up_file: Path,
    r_values_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    base_energy_grid_kj_mol: np.ndarray,
    bead_types: list,
    bead_charges: list,
    pair_params: dict,
    core_min_nm: float,
    frame_start_fraction: float = 0.5,
    max_frames: int = 100,
    min_count: int = 4,
    mode: str = "radial",
) -> dict:
    up_file = Path(up_file).expanduser().resolve()
    mode = str(mode).strip().lower()
    if mode not in {"radial", "generalized"}:
        raise ValueError(f"force_match_mode must be 'radial' or 'generalized', got {mode!r}")
    r_values_nm = np.asarray(r_values_nm, dtype=np.float64)
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    r_edges = _grid_edges_from_centers(r_values_nm, lo=0.0)
    angle_edges = _grid_edges_from_centers(cos_theta_grid, lo=-1.0, hi=1.0)
    counts = np.zeros(base_energy_grid_kj_mol.shape, dtype=np.float64)
    dudr_sum = np.zeros_like(counts)
    duda1_sum = np.zeros_like(counts)
    duda2_sum = np.zeros_like(counts)
    angular_counts = np.zeros_like(counts)

    with h5py.File(up_file, "r") as h5:
        if "/output/pos" not in h5:
            raise RuntimeError(f"{up_file} is missing /output/pos")
        lipid_molecules = _dopc_lipid_molecules_from_upside_h5(h5, up_file)
        pos = h5["/output/pos"]
        n_frame = int(pos.shape[0])
        first_frame = int(round(float(frame_start_fraction) * float(max(0, n_frame - 1))))
        frame_indices = list(range(first_frame, n_frame))
        if not frame_indices:
            frame_indices = [n_frame - 1]
        max_frames = max(1, int(max_frames))
        if len(frame_indices) > max_frames:
            select = np.linspace(0, len(frame_indices) - 1, max_frames)
            frame_indices = [frame_indices[int(round(float(x)))] for x in select]
        box_ang = _read_output_box_lengths_ang(h5)

        for frame in frame_indices:
            bead_frames_ang, centers_ang, directions = _extract_dopc_beads_centers_orientations_from_frame(
                np.asarray(pos[frame, 0, :, :], dtype=np.float64),
                lipid_molecules,
                box_ang,
            )
            for i in range(len(bead_frames_ang) - 1):
                for j in range(i + 1, len(bead_frames_ang)):
                    delta_ang = centers_ang[j] - centers_ang[i]
                    if box_ang is not None:
                        delta_ang -= box_ang * np.round(delta_ang / box_ang)
                    r_ang = float(np.linalg.norm(delta_ang))
                    if r_ang <= 1.0e-12:
                        continue
                    r_nm = r_ang * ANGSTROM_TO_NM
                    ir = int(np.searchsorted(r_edges, r_nm, side="right") - 1)
                    if ir < 0 or ir >= r_values_nm.size:
                        continue
                    n12 = delta_ang / r_ang
                    a1 = float(-np.dot(directions[i], n12))
                    a2 = float(np.dot(directions[j], n12))
                    ia1 = int(np.searchsorted(angle_edges, a1, side="right") - 1)
                    ia2 = int(np.searchsorted(angle_edges, a2, side="right") - 1)
                    if ia1 < 0 or ia1 >= cos_theta_grid.size or ia2 < 0 or ia2 >= cos_theta_grid.size:
                        continue

                    shifted_j_ang = bead_frames_ang[j] + (centers_ang[i] + delta_ang - centers_ang[j])[None, :]
                    _, grad_i, grad_j = _compute_pair_energy_and_gradient(
                        bead_frames_ang[i] * ANGSTROM_TO_NM,
                        shifted_j_ang * ANGSTROM_TO_NM,
                        bead_types,
                        bead_types,
                        bead_charges,
                        bead_charges,
                        pair_params,
                        dist_min_nm=NUMERICAL_DISTANCE_GUARD_NM,
                    )
                    force_i = -np.sum(grad_i, axis=0)
                    force_j = -np.sum(grad_j, axis=0)
                    dudr = float(np.dot(force_i, n12))

                    counts[ir, ia1, ia2] += 1.0
                    dudr_sum[ir, ia1, ia2] += dudr

                    if mode == "generalized":
                        center_i_nm = centers_ang[i] * ANGSTROM_TO_NM
                        center_j_shift_nm = (centers_ang[i] + delta_ang) * ANGSTROM_TO_NM
                        beads_i_nm = bead_frames_ang[i] * ANGSTROM_TO_NM
                        shifted_j_nm = shifted_j_ang * ANGSTROM_TO_NM
                        torque_i = np.sum(np.cross(beads_i_nm - center_i_nm[None, :], -grad_i), axis=0)
                        torque_j = np.sum(np.cross(shifted_j_nm - center_j_shift_nm[None, :], -grad_j), axis=0)
                        axis1 = np.cross(directions[i], n12)
                        axis2 = np.cross(n12, directions[j])
                        axis1_norm2 = float(np.dot(axis1, axis1))
                        axis2_norm2 = float(np.dot(axis2, axis2))
                        if axis1_norm2 > 1.0e-8 and axis2_norm2 > 1.0e-8:
                            duda1 = float(np.dot(torque_i, axis1) / axis1_norm2)
                            duda2 = float(np.dot(torque_j, axis2) / axis2_norm2)
                            angular_counts[ir, ia1, ia2] += 1.0
                            duda1_sum[ir, ia1, ia2] += duda1
                            duda2_sum[ir, ia1, ia2] += duda2

                    e_rev = -n12
                    dudr_rev = float(np.dot(force_j, e_rev))
                    counts[ir, ia2, ia1] += 1.0
                    dudr_sum[ir, ia2, ia1] += dudr_rev
                    if mode == "generalized":
                        axis1_rev = np.cross(directions[j], e_rev)
                        axis2_rev = np.cross(e_rev, directions[i])
                        axis1_rev_norm2 = float(np.dot(axis1_rev, axis1_rev))
                        axis2_rev_norm2 = float(np.dot(axis2_rev, axis2_rev))
                        if axis1_rev_norm2 > 1.0e-8 and axis2_rev_norm2 > 1.0e-8:
                            duda1_rev = float(np.dot(torque_j, axis1_rev) / axis1_rev_norm2)
                            duda2_rev = float(np.dot(torque_i, axis2_rev) / axis2_rev_norm2)
                            angular_counts[ir, ia2, ia1] += 1.0
                            duda1_sum[ir, ia2, ia1] += duda1_rev
                            duda2_sum[ir, ia2, ia1] += duda2_rev

    mean_dudr = np.full_like(counts, np.nan, dtype=np.float64)
    sampled = counts >= float(max(1, int(min_count)))
    mean_dudr[sampled] = dudr_sum[sampled] / counts[sampled]
    mean_duda1 = np.full_like(counts, np.nan, dtype=np.float64)
    mean_duda2 = np.full_like(counts, np.nan, dtype=np.float64)
    angular_sampled = angular_counts >= float(max(1, int(min_count)))
    mean_duda1[angular_sampled] = duda1_sum[angular_sampled] / angular_counts[angular_sampled]
    mean_duda2[angular_sampled] = duda2_sum[angular_sampled] / angular_counts[angular_sampled]
    energy = np.asarray(base_energy_grid_kj_mol, dtype=np.float64).copy()
    updated = np.zeros_like(sampled, dtype=bool)
    core_mask = r_values_nm[:, None, None] >= float(core_min_nm)
    sampled &= core_mask
    lsq_info = {
        "lsq_variable_count": 0,
        "lsq_equation_count": 0,
        "lsq_residual_rms_kj_mol": 0.0,
        "lsq_solver": "radial_segment_integration",
    }

    if mode == "generalized":
        generalized_sampled = sampled & angular_sampled
        energy, updated, lsq_info = _least_squares_generalized_force_energy_grid(
            r_values_nm,
            cos_theta_grid,
            base_energy_grid_kj_mol,
            generalized_sampled,
            counts,
            mean_dudr,
            mean_duda1,
            mean_duda2,
            min_count=min_count,
        )
        sampled = generalized_sampled
    else:
        for ia1 in range(cos_theta_grid.size):
            for ia2 in range(cos_theta_grid.size):
                valid = np.where(sampled[:, ia1, ia2])[0]
                if valid.size < 2:
                    continue
                breaks = np.where(np.diff(valid) > 1)[0]
                segment_starts = np.r_[0, breaks + 1]
                segment_ends = np.r_[breaks, valid.size - 1]
                for start_idx, end_idx in zip(segment_starts, segment_ends):
                    segment = valid[start_idx : end_idx + 1]
                    if segment.size < 2:
                        continue
                    high = int(segment[-1])
                    energy[high, ia1, ia2] = base_energy_grid_kj_mol[high, ia1, ia2]
                    updated[high, ia1, ia2] = True
                    for pos_idx in range(segment.size - 2, -1, -1):
                        lo = int(segment[pos_idx])
                        hi = int(segment[pos_idx + 1])
                        dr = float(r_values_nm[hi] - r_values_nm[lo])
                        avg_dudr = 0.5 * (float(mean_dudr[lo, ia1, ia2]) + float(mean_dudr[hi, ia1, ia2]))
                        energy[lo, ia1, ia2] = energy[hi, ia1, ia2] - avg_dudr * dr
                        updated[lo, ia1, ia2] = True

    print(
        f"  CGL-CGL {mode} force matching from {up_file}: "
        f"{int(np.sum(counts))} ordered pair-force samples, "
        f"{int(np.count_nonzero(updated))}/{updated.size} tensor bins updated"
    )
    return {
        "energy_grid_kj_mol": energy,
        "counts": counts,
        "mean_dudr_kj_mol_nm": mean_dudr,
        "mean_duda1_kj_mol": mean_duda1,
        "mean_duda2_kj_mol": mean_duda2,
        "angular_counts": angular_counts,
        "updated": updated,
        "mode": mode,
        "frame_count": int(len(frame_indices)),
        "pair_force_sample_count": int(np.sum(counts)),
        "updated_tensor_bins": int(np.count_nonzero(updated)),
        "sampled_tensor_bins": int(np.count_nonzero(sampled)),
        "frame_start_fraction": float(frame_start_fraction),
        "max_frames": int(max_frames),
        "min_count": int(min_count),
        **lsq_info,
    }


def _fit_cg_lipid_quadspline(
    ref_bead_positions_nm: np.ndarray,
    bead_types: list,
    bead_charges: list,
    pair_params: dict,
    r_min_nm: float = 0.25,
    r_max_nm: float = 0.70,
    r_count: int = 16,
    cos_theta_count: int = 7,
    azimuthal_count: int = 2,
    bead_frame_count: int = 1,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    knot_spacing_ang: float = 1.4,
    n_modes: int = 4,
    n_knot_radial: int = 14,
    n_knot_angular: int = 9,
    cg_smooth: float = 0.01,
    temperature: float = 0.0,
    average_temperature: float | None = None,
    fluid_pair_pmf_upside_h5_path: Path | None = None,
    fluid_pair_pmf_frame_start_fraction: float = 0.5,
    fluid_pair_pmf_core_min_nm: float = 0.0,
    force_match_upside_h5_path: Path | None = None,
    force_match_frame_start_fraction: float = 0.5,
    force_match_max_frames: int = 100,
    force_match_min_count: int = 4,
    force_match_mode: str = "radial",
    ibi_base_pair_h5_path: Path | None = None,
    ibi_target_upside_h5_path: Path | None = None,
    ibi_model_upside_h5_path: Path | None = None,
    ibi_target_frame_start_fraction: float = 0.5,
    ibi_model_frame_start_fraction: float = 0.5,
    ibi_step_size: float = 1.0,
) -> dict:
    """Fit full tensor-product B-spline parameters for CG lipid-CG lipid interactions.

    The runtime table is stored as a log1p-reduced energy control using
    ``temperature`` as the numerical transform scale.  ``average_temperature``
    controls unresolved azimuthal/bead-frame averaging; values <= 0 use direct
    energy expectation instead of a pair PMF.
    """
    if average_temperature is None:
        average_temperature = temperature
    average_temperature = float(average_temperature)
    r_values = np.asarray(_linspace(r_min_nm, r_max_nm, r_count), dtype=np.float64)
    cos_theta_grid = np.asarray(_linspace(-1.0, 1.0, cos_theta_count), dtype=np.float64)
    n_angle = cos_theta_count
    n_radial = r_count

    phi_values = np.linspace(0.0, 2.0 * np.pi, azimuthal_count, endpoint=False)
    n12_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    dirs1 = _directions_with_dot_np(-n12_axis, cos_theta_grid, phi_values)
    dirs2 = _directions_with_dot_np(n12_axis, cos_theta_grid, phi_values)

    ref_nm = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    conformer_pairs = _cg_lipid_conformer_pair_indices(ref_nm.shape[0])

    bead_frame_angles = _bead_frame_angles(bead_frame_count)
    if ibi_base_pair_h5_path is not None:
        energy_grid = _load_cg_lipid_pair_energy_grid_from_h5(
            ibi_base_pair_h5_path,
            (n_radial, n_angle, n_angle),
            r_values,
        )
    else:
        energy_grid = np.zeros((n_radial, n_angle, n_angle), dtype=np.float64)
        total_samples = (
            n_radial * n_angle * n_angle
            * azimuthal_count * azimuthal_count
            * len(bead_frame_angles) * len(bead_frame_angles)
            * len(conformer_pairs)
        )
        print(f"  Sampling CG-CG energy: {n_radial} radial x {n_angle}^2 angular "
              f"x {azimuthal_count}^2 azimuthal x {len(bead_frame_angles)}^2 bead-frame "
              f"x {len(conformer_pairs)} conformer-pairs = {total_samples} samples, "
              "direct rotated ensemble geometry, full tensor")

        tasks = []
        for ir, r_nm in enumerate(r_values):
            for ia1 in range(n_angle):
                for ia2 in range(n_angle):
                    tasks.append(
                        {
                            "ir": ir,
                            "ia1": ia1,
                            "ia2": ia2,
                            "r_nm": float(r_nm),
                            "dirs1": dirs1[ia1],
                            "dirs2": dirs2[ia2],
                            "bead_frame_angles": bead_frame_angles,
                            "ref_nm": ref_nm,
                            "bead_types": list(bead_types),
                            "bead_charges": list(bead_charges),
                            "pair_params": pair_params,
                            "dist_min_nm": float(dist_min_nm),
                            "temperature": average_temperature,
                            "conformer_pairs": conformer_pairs,
                        }
                    )
        for ir, ia1, ia2, energy in _parallel_map_ordered("CG-CG table", _run_cg_pair_tensor_task, tasks):
            energy_grid[ir, ia1, ia2] = energy

    fluid_pmf = None
    if fluid_pair_pmf_upside_h5_path is not None and ibi_base_pair_h5_path is None:
        fluid_pmf = _fluid_bilayer_pair_pmf_from_upside_h5(
            fluid_pair_pmf_upside_h5_path,
            r_values,
            cos_theta_grid,
            temperature_upside=temperature,
            frame_start_fraction=fluid_pair_pmf_frame_start_fraction,
        )
        occupied = np.asarray(fluid_pmf["occupied"], dtype=bool)
        pmf = np.asarray(fluid_pmf["pmf_kj_mol"], dtype=np.float64)
        if float(fluid_pair_pmf_core_min_nm) > 0.0:
            occupied &= r_values[:, None, None] >= float(fluid_pair_pmf_core_min_nm)
        energy_grid[occupied] = pmf[occupied]
        print(
            "  CGL-CGL: using fluid-bilayer orientation-resolved PMF for occupied bins "
            "outside the direct dry-MARTINI overlap core"
        )

    force_match = None
    if force_match_upside_h5_path is not None:
        force_match = _force_matched_pair_energy_from_upside_h5(
            force_match_upside_h5_path,
            r_values,
            cos_theta_grid,
            energy_grid,
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            core_min_nm=fluid_pair_pmf_core_min_nm,
            frame_start_fraction=force_match_frame_start_fraction,
            max_frames=force_match_max_frames,
            min_count=force_match_min_count,
            mode=force_match_mode,
        )
        energy_grid = np.asarray(force_match["energy_grid_kj_mol"], dtype=np.float64)

    ibi_correction = None
    if ibi_model_upside_h5_path is not None:
        if ibi_target_upside_h5_path is None:
            if fluid_pair_pmf_upside_h5_path is None:
                raise ValueError("IBI correction requires a target trajectory or fluid_pair_pmf_upside_h5_path")
            ibi_target_upside_h5_path = fluid_pair_pmf_upside_h5_path
        ibi_correction = _relative_entropy_pair_ibi_correction_from_upside_h5(
            ibi_target_upside_h5_path,
            ibi_model_upside_h5_path,
            r_values,
            cos_theta_grid,
            temperature_upside=temperature,
            target_frame_start_fraction=ibi_target_frame_start_fraction,
            model_frame_start_fraction=ibi_model_frame_start_fraction,
            step_size=ibi_step_size,
        )
        ibi_mask = np.asarray(ibi_correction["sampled"], dtype=bool)
        if float(fluid_pair_pmf_core_min_nm) > 0.0:
            ibi_mask &= r_values[:, None, None] >= float(fluid_pair_pmf_core_min_nm)
        delta_u = np.asarray(ibi_correction["correction_kj_mol"], dtype=np.float64)
        energy_grid[ibi_mask] += delta_u[ibi_mask]
        print(
            "  CGL-CGL: applied relative-entropy/IBI correction in "
            f"{int(np.count_nonzero(ibi_mask))}/{ibi_mask.size} sampled non-core bins "
            f"with step={float(ibi_correction['step_size']):.3f}"
        )

    reference_energy_kj_mol = float(np.min(energy_grid))
    kbt_kj_mol = float(temperature) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt_kj_mol <= 0.0:
        raise RuntimeError("CGL-CGL log1p reduced-PMF table requires positive temperature")
    reduced_energy = (energy_grid - reference_energy_kj_mol) / kbt_kj_mol
    control_grid = np.log1p(np.maximum(reduced_energy, 0.0))

    radial_mean = energy_grid.mean(axis=(1, 2))
    attractive_background = np.zeros_like(radial_mean)

    knot_spacing = float(knot_spacing_ang)
    tensor = _fit_radial_angular_angular_tensor_bspline(
        r_values,
        cos_theta_grid,
        control_grid,
        n_knot_radial=n_knot_radial,
        n_knot_angular=n_knot_angular,
        knot_spacing_ang=knot_spacing,
        energy_conversion=1.0,
        smooth=cg_smooth,
    )
    r_min_ang = float(r_min_nm) * LENGTH_CONVERSION_A_PER_NM
    n_core = max(0, min(n_knot_radial - 1, int(math.ceil((r_min_ang - 1e-6) / knot_spacing))))
    short_range_core_kj_mol = float(np.max(energy_grid[0]))
    short_range_core_eup = _fit_angular_angular_tensor_bspline(
        cos_theta_grid,
        control_grid[0],
        n_knot_angular=n_knot_angular,
        energy_conversion=1.0,
        smooth=cg_smooth,
    )
    if n_core > 0:
        tensor[:n_core, :, :] = np.maximum(tensor[:n_core, :, :], short_range_core_eup[None, :, :])
    control_min = float(np.min(control_grid))
    control_max = float(np.max(control_grid))
    tensor = np.clip(tensor, control_min, control_max)
    attractive_control_count = int(np.count_nonzero(energy_grid < 0.0))
    rms_error = 0.0
    interaction_param = tensor.reshape(-1)

    return {
        "tensor_knots": tensor,
        "interaction_param": interaction_param,
        "rms_error": rms_error,
        "energy_grid_raw": energy_grid,
        "reference_energy_eup": reference_energy_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
        "attractive_radial_background_kj_mol": attractive_background,
        "azimuthal_average": (
            "tempered_boltzmann_free_energy"
            if average_temperature > 0.0 and abs(average_temperature - float(temperature)) > 1e-8
            else ("energy_expectation" if average_temperature <= 0.0 else "boltzmann_free_energy")
        ),
        "azimuthal_average_temperature_upside": float(average_temperature),
        "isotropic_background_source": (
            "full_dopc_projected_generalized_force_matching"
            if force_match is not None and force_match["mode"] == "generalized"
            else "full_dopc_projected_radial_force_matching"
            if force_match is not None
            else (
            "fluid_bilayer_orientation_resolved_pair_pmf"
            if fluid_pmf is not None
            else "none_full_resolved_dry_martini_pair_table"
            )
        ),
        "excluded_area_source": (
            "direct_dry_martini_overlap_core_for_force_matched_table"
            if force_match is not None
            else (
            "direct_dry_martini_overlap_core_for_fluid_pmf_table"
            if fluid_pmf is not None
            else "none_full_resolved_dry_martini_pair_table"
            )
        ),
        "excluded_area_nonnegative_rows": 0,
        "attractive_control_source": (
            "full_dopc_projected_generalized_force_matching"
            if force_match is not None and force_match["mode"] == "generalized"
            else "full_dopc_projected_radial_force_matching"
            if force_match is not None
            else (
            "fluid_bilayer_orientation_resolved_pair_pmf"
            if fluid_pmf is not None
            else "retained_full_resolved_dry_martini_pair_table"
            )
        ),
        "attractive_control_count": attractive_control_count,
        "core_boundary_source": (
            "force_matched_generalized_projection_plus_direct_dry_martini_overlap_core"
            if force_match is not None and force_match["mode"] == "generalized"
            else "force_matched_radial_projection_plus_direct_dry_martini_overlap_core"
            if force_match is not None
            else (
            "fluid_bilayer_radial_pmf_plus_direct_dry_martini_overlap_core"
            if fluid_pmf is not None
            else "angular_resolved_first_sampled_dry_martini_energy"
            )
        ),
        "core_boundary_row": int(n_core),
        "unresolved_core_rows": int(n_core),
        "unresolved_core_energy_kj_mol": short_range_core_kj_mol,
        "r_values_nm": r_values,
        "cos_theta_grid": cos_theta_grid,
        "knot_spacing_ang": knot_spacing,
        "cutoff_ang": float((n_knot_radial - 2) * knot_spacing),
        "n_modes": 0,
        "n_radial": int(n_knot_radial),
        "n_angular": int(n_knot_angular),
        "fit_smooth": float(cg_smooth),
        "azimuthal_count": int(azimuthal_count),
        "bead_frame_count": int(len(bead_frame_angles)),
        "conformer_count": int(ref_nm.shape[0]),
        "conformer_pair_count": int(len(conformer_pairs)),
        "schema": "cg_lipid_pair_full",
        "sample_dist_min_nm": float(dist_min_nm),
        "ibi_base_pair_h5_path": (
            "" if ibi_base_pair_h5_path is None else str(Path(ibi_base_pair_h5_path).expanduser().resolve())
        ),
        "fluid_pmf_pair_sample_count": 0 if fluid_pmf is None else int(fluid_pmf["pair_sample_count"]),
        "fluid_pmf_occupied_radial_bins": 0 if fluid_pmf is None else int(fluid_pmf["occupied_radial_bins"]),
        "fluid_pmf_occupied_tensor_bins": 0 if fluid_pmf is None else int(fluid_pmf["occupied_tensor_bins"]),
        "fluid_pmf_frame_start_fraction": (
            0.0 if fluid_pmf is None else float(fluid_pmf["frame_start_fraction"])
        ),
        "fluid_pmf_core_min_nm": float(fluid_pair_pmf_core_min_nm) if fluid_pmf is not None else 0.0,
        "force_match_enabled": int(force_match is not None),
        "force_match_frame_count": 0 if force_match is None else int(force_match["frame_count"]),
        "force_match_pair_force_sample_count": 0 if force_match is None else int(force_match["pair_force_sample_count"]),
        "force_match_sampled_tensor_bins": 0 if force_match is None else int(force_match["sampled_tensor_bins"]),
        "force_match_updated_tensor_bins": 0 if force_match is None else int(force_match["updated_tensor_bins"]),
        "force_match_frame_start_fraction": 0.0 if force_match is None else float(force_match["frame_start_fraction"]),
        "force_match_max_frames": 0 if force_match is None else int(force_match["max_frames"]),
        "force_match_min_count": 0 if force_match is None else int(force_match["min_count"]),
        "force_match_mode": "" if force_match is None else str(force_match["mode"]),
        "force_match_lsq_variable_count": 0 if force_match is None else int(force_match["lsq_variable_count"]),
        "force_match_lsq_equation_count": 0 if force_match is None else int(force_match["lsq_equation_count"]),
        "force_match_lsq_residual_rms_kj_mol": (
            0.0 if force_match is None else float(force_match["lsq_residual_rms_kj_mol"])
        ),
        "force_match_lsq_solver": "" if force_match is None else str(force_match["lsq_solver"]),
        "force_match_counts": (
            np.zeros_like(energy_grid, dtype=np.float64)
            if force_match is None
            else np.asarray(force_match["counts"], dtype=np.float64)
        ),
        "force_match_mean_dudr_kj_mol_nm": (
            np.full_like(energy_grid, np.nan, dtype=np.float64)
            if force_match is None
            else np.asarray(force_match["mean_dudr_kj_mol_nm"], dtype=np.float64)
        ),
        "force_match_mean_duda1_kj_mol": (
            np.full_like(energy_grid, np.nan, dtype=np.float64)
            if force_match is None
            else np.asarray(force_match["mean_duda1_kj_mol"], dtype=np.float64)
        ),
        "force_match_mean_duda2_kj_mol": (
            np.full_like(energy_grid, np.nan, dtype=np.float64)
            if force_match is None
            else np.asarray(force_match["mean_duda2_kj_mol"], dtype=np.float64)
        ),
        "force_match_angular_counts": (
            np.zeros_like(energy_grid, dtype=np.float64)
            if force_match is None
            else np.asarray(force_match["angular_counts"], dtype=np.float64)
        ),
        "force_match_updated": (
            np.zeros_like(energy_grid, dtype=np.int8)
            if force_match is None
            else np.asarray(force_match["updated"], dtype=np.int8)
        ),
        "ibi_enabled": int(ibi_correction is not None),
        "ibi_target_pair_sample_count": 0 if ibi_correction is None else int(ibi_correction["target_pair_sample_count"]),
        "ibi_model_pair_sample_count": 0 if ibi_correction is None else int(ibi_correction["model_pair_sample_count"]),
        "ibi_target_occupied_tensor_bins": 0 if ibi_correction is None else int(ibi_correction["target_occupied_tensor_bins"]),
        "ibi_model_occupied_tensor_bins": 0 if ibi_correction is None else int(ibi_correction["model_occupied_tensor_bins"]),
        "ibi_sampled_tensor_bins": 0 if ibi_correction is None else int(ibi_correction["sampled_tensor_bins"]),
        "ibi_target_frame_start_fraction": 0.0 if ibi_correction is None else float(ibi_correction["target_frame_start_fraction"]),
        "ibi_model_frame_start_fraction": 0.0 if ibi_correction is None else float(ibi_correction["model_frame_start_fraction"]),
        "ibi_step_size": 0.0 if ibi_correction is None else float(ibi_correction["step_size"]),
        "ibi_raw_correction_min_kj_mol": (
            0.0 if ibi_correction is None else float(ibi_correction["raw_correction_min_kj_mol"])
        ),
        "ibi_raw_correction_max_kj_mol": (
            0.0 if ibi_correction is None else float(ibi_correction["raw_correction_max_kj_mol"])
        ),
        "ibi_raw_correction_mean_kj_mol": (
            0.0 if ibi_correction is None else float(ibi_correction["raw_correction_mean_kj_mol"])
        ),
        "ibi_correction_min_kj_mol": 0.0 if ibi_correction is None else float(ibi_correction["correction_min_kj_mol"]),
        "ibi_correction_max_kj_mol": 0.0 if ibi_correction is None else float(ibi_correction["correction_max_kj_mol"]),
        "ibi_correction_mean_kj_mol": 0.0 if ibi_correction is None else float(ibi_correction["correction_mean_kj_mol"]),
        "ibi_correction_grid_kj_mol": (
            np.full_like(energy_grid, np.nan, dtype=np.float64)
            if ibi_correction is None
            else np.asarray(ibi_correction["correction_kj_mol"], dtype=np.float64)
        ),
        "ibi_target_counts": (
            np.zeros_like(energy_grid, dtype=np.float64)
            if ibi_correction is None
            else np.asarray(ibi_correction["target_counts"], dtype=np.float64)
        ),
        "ibi_model_counts": (
            np.zeros_like(energy_grid, dtype=np.float64)
            if ibi_correction is None
            else np.asarray(ibi_correction["model_counts"], dtype=np.float64)
        ),
        "energy_transform": (
            "log1p_reduced_energy_expectation"
            if average_temperature <= 0.0
            else (
                "log1p_reduced_tempered_pmf"
                if abs(average_temperature - float(temperature)) > 1e-8
                else "log1p_reduced_pmf"
            )
        ),
        "spline_control_quantity": (
            "log1p_reduced_energy_expectation"
            if average_temperature <= 0.0
            else (
                "log1p_reduced_tempered_free_energy"
                if abs(average_temperature - float(temperature)) > 1e-8
                else "log1p_reduced_free_energy"
            )
        ),
    }


def _fit_cg_lipid_sc_quadspline(
    ref_bead_positions_nm: np.ndarray,
    cg_bead_types: list,
    cg_bead_charges: list,
    target_type: str,
    target_charge: float,
    rotamer_bead_positions_nm: list,
    rotamer_weights: list,
    sc_bead_types: list,
    sc_bead_charges: list,
    pair_params: dict,
    cb_anchor_nm: list,
    cb_vector_unit: list,
    r_min_nm: float = 0.30,
    r_max_nm: float = 0.70,
    r_count: int = 16,
    cos_theta_count: int = 9,
    azimuthal_count: int = 4,
    sidechain_bead_frame_count: int = 1,
    cg_bead_frame_count: int = 1,
    n_modes: int = 4,
    n_knot_radial: int = 14,
    n_knot_angular: int = 15,
    angular_smooth: float = 0.01,
    knot_spacing_ang: float = 1.4,
    temperature: float = 0.0,
    average_temperature: float | None = None,
) -> dict:
    """Fit full multimode B-spline params for sidechain-CG lipid interactions.

    The parameter layout matches _fit_cg_lipid_quadspline:
      V0(r) + sum_m Ang1_m(a_sc) * Ang2_m(a_cg) * Vm(r)
    where source1 is the sidechain vector and source2 is the CG lipid vector.
    """
    r_values = np.asarray(_linspace(r_min_nm, r_max_nm, r_count), dtype=np.float64)
    cos_theta_grid = np.asarray(_linspace(-1.0, 1.0, cos_theta_count), dtype=np.float64)
    n_angle = cos_theta_count
    n_radial = r_count
    if average_temperature is None:
        average_temperature = temperature
    average_temperature = float(average_temperature)

    phi_values = np.linspace(0.0, 2.0 * np.pi, azimuthal_count, endpoint=False)
    sidechain_bead_frame_angles = _bead_frame_angles(sidechain_bead_frame_count)
    cg_bead_frame_angles = _bead_frame_angles(cg_bead_frame_count)

    ref_nm = _canonicalize_lipid_reference_ensemble_to_z(ref_bead_positions_nm)
    cb_vec = np.asarray(cb_vector_unit, dtype=np.float64)
    cb_vec /= max(float(np.linalg.norm(cb_vec)), 1e-12)
    sc_to_cg_dirs = _directions_with_dot_np(-cb_vec, cos_theta_grid, phi_values)
    n_rotamer = len(rotamer_bead_positions_nm)

    rotamer_sc_positions = []
    n_sc_bead = len(sc_bead_types)
    for irot in range(n_rotamer):
        pos = np.asarray(rotamer_bead_positions_nm[irot], dtype=np.float64)
        if pos.shape != (n_sc_bead, 3):
            raise RuntimeError(
                f"CG-SC {target_type} rotamer {irot} has bead geometry shape {pos.shape}, "
                f"expected ({n_sc_bead}, 3)"
            )
        rotamer_sc_positions.append(pos)

    print(f"  Sampling CG-SC energy for target={target_type}: "
          f"{n_radial} radial x {n_angle}^2 angular x {azimuthal_count}^2 azimuthal "
          f"x {len(sidechain_bead_frame_angles)} SC bead-frame x {len(cg_bead_frame_angles)} CGL bead-frame, "
          f"x {ref_nm.shape[0]} conformers, direct rotated ensemble geometry, modes={n_modes}")

    # Energy grid: (n_radial, n_angle_sc, n_angle_cg), matching runtime source order.
    energy_grid = np.zeros((n_radial, n_angle, n_angle), dtype=np.float64)

    for ir, r_nm in enumerate(r_values):
        for ia_sc in range(n_angle):
            for ia_cg in range(n_angle):
                sample_energies = []
                sample_weights = []
                for ip_sc in range(azimuthal_count):
                    dir_to_cg = sc_to_cg_dirs[ia_sc, ip_sc]
                    cg_com = r_nm * dir_to_cg
                    for ip_cg in range(azimuthal_count):
                        dir_cg = _direction_with_dot_np(
                            dir_to_cg, cos_theta_grid[ia_cg], phi_values[ip_cg]
                        )

                        R_cg_base = _rotation_to_align_z_np(dir_cg)

                        for irot in range(n_rotamer):
                            w = rotamer_weights[irot]
                            if w <= 0.0:
                                continue

                            sc_positions = rotamer_sc_positions[irot]
                            for sc_frame_angle in sidechain_bead_frame_angles:
                                framed_sc_positions = _rotate_points_about_axis_np(
                                    sc_positions,
                                    cb_vec,
                                    float(sc_frame_angle),
                                    np.asarray(cb_anchor_nm, dtype=np.float64),
                                )
                                for cg_frame_angle in cg_bead_frame_angles:
                                    R_cg = R_cg_base @ _rotation_about_axis_np(
                                        np.array([0.0, 0.0, 1.0], dtype=np.float64),
                                        float(cg_frame_angle),
                                    )
                                    for ref_conf in ref_nm:
                                        cg_positions = cg_com[None, :] + (R_cg @ ref_conf.T).T

                                        pair_energy, _, _ = _compute_pair_energy_and_gradient(
                                            cg_positions, framed_sc_positions, cg_bead_types, sc_bead_types,
                                            cg_bead_charges, sc_bead_charges, pair_params,
                                            dist_min_nm=NUMERICAL_DISTANCE_GUARD_NM,
                                        )
                                        sample_energies.append(float(pair_energy))
                                        sample_weights.append(float(w))

                energy_grid[ir, ia_sc, ia_cg] = _boltzmann_free_energy_kj_mol(
                    sample_energies,
                    average_temperature,
                    sample_weights,
                )

    if int(n_modes) == 0:
        reference_energy_kj_mol = float(np.min(energy_grid))
        kbt_kj_mol = float(temperature) * ENERGY_CONVERSION_KJ_PER_EUP
        if kbt_kj_mol <= 0.0:
            raise RuntimeError("SC-CGL log1p reduced table requires positive transform temperature")
        reduced_energy = (energy_grid - reference_energy_kj_mol) / kbt_kj_mol
        control_grid = np.log1p(np.maximum(reduced_energy, 0.0))
        tensor = _fit_radial_angular_angular_tensor_bspline(
            r_values,
            cos_theta_grid,
            control_grid,
            n_knot_radial=n_knot_radial,
            n_knot_angular=n_knot_angular,
            knot_spacing_ang=knot_spacing_ang,
            energy_conversion=1.0,
            smooth=angular_smooth,
        )
        r_min_ang = float(r_min_nm) * LENGTH_CONVERSION_A_PER_NM
        n_core = max(0, min(n_knot_radial - 1, int(math.ceil((r_min_ang - 1e-6) / float(knot_spacing_ang)))))
        core_control = _fit_angular_angular_tensor_bspline(
            cos_theta_grid,
            control_grid[0],
            n_knot_angular=n_knot_angular,
            energy_conversion=1.0,
            smooth=angular_smooth,
        )
        if n_core > 0:
            tensor[:n_core, :, :] = np.maximum(tensor[:n_core, :, :], core_control[None, :, :])
        control_min = float(np.min(control_grid))
        control_max = float(np.max(control_grid))
        tensor = np.clip(tensor, control_min, control_max)
        recon = control_grid
        rms_error = 0.0
        short_range_core_kj_mol = float(np.max(energy_grid[0]))
        return {
            "tensor_knots": tensor,
            "interaction_param": tensor.reshape(-1),
            "rms_error": rms_error,
            "v_radial_raw": energy_grid.mean(axis=(1, 2)),
            "energy_grid_raw": energy_grid,
            "reference_energy_eup": reference_energy_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP,
            "short_range_core_kj_mol": short_range_core_kj_mol,
            "raw_short_range_core_kj_mol": short_range_core_kj_mol,
            "short_range_core_rows": int(n_core),
            "short_range_core_source": "angular_resolved_first_sampled_dry_martini_energy",
            "excluded_area_source": "none_full_resolved_dry_martini_sc_cgl_table",
            "excluded_area_nonnegative_rows": 0,
            "isotropic_background_source": "none_full_resolved_dry_martini_sc_cgl_table",
            "attractive_control_source": "retained_full_resolved_dry_martini_sc_cgl_table",
            "ang1_raw": np.zeros((0, n_angle), dtype=np.float64),
            "ang2_raw": np.zeros((0, n_angle), dtype=np.float64),
            "v_angular_raw": np.zeros((0, n_radial), dtype=np.float64),
            "r_values_nm": r_values,
            "cos_theta_grid": cos_theta_grid,
            "knot_spacing_ang": float(knot_spacing_ang),
            "cutoff_ang": float((n_knot_radial - 2) * float(knot_spacing_ang)),
            "taper_width_ang": float(knot_spacing_ang),
            "n_modes": 0,
            "n_radial": int(n_knot_radial),
            "n_angular": int(n_knot_angular),
            "angular_smooth": float(angular_smooth),
            "azimuthal_count": int(azimuthal_count),
            "azimuthal_average": (
                "energy_expectation"
                if average_temperature <= 0.0
                else (
                    "tempered_boltzmann_free_energy"
                    if abs(average_temperature - float(temperature)) > 1e-8
                    else "boltzmann_free_energy"
                )
            ),
            "azimuthal_average_temperature_upside": float(average_temperature),
            "sidechain_bead_frame_count": int(len(sidechain_bead_frame_angles)),
            "cg_bead_frame_count": int(len(cg_bead_frame_angles)),
            "conformer_count": int(ref_nm.shape[0]),
            "energy_transform": (
                "log1p_reduced_energy_expectation"
                if average_temperature <= 0.0
                else (
                    "log1p_reduced_tempered_pmf"
                    if abs(average_temperature - float(temperature)) > 1e-8
                    else "log1p_reduced_pmf"
                )
            ),
            "spline_control_quantity": (
                "log1p_reduced_energy_expectation"
                if average_temperature <= 0.0
                else (
                    "log1p_reduced_tempered_free_energy"
                    if abs(average_temperature - float(temperature)) > 1e-8
                    else "log1p_reduced_free_energy"
                )
            ),
        }

    radial_mean = energy_grid.mean(axis=(1, 2))
    attractive_background = np.zeros_like(radial_mean)
    v_radial = energy_grid.mean(axis=(1, 2))
    residual_fit = energy_grid - v_radial[:, None, None]

    mode_ang1 = []
    mode_ang2 = []
    mode_radial = []
    residual_remaining = residual_fit.copy()
    for _ in range(int(n_modes)):
        u_all = np.zeros((n_radial, n_angle), dtype=np.float64)
        v_all = np.zeros((n_radial, n_angle), dtype=np.float64)
        s_all = np.zeros(n_radial, dtype=np.float64)
        for ir in range(n_radial):
            mat = residual_remaining[ir]
            if np.allclose(mat, 0.0):
                continue
            u, s, vh = np.linalg.svd(mat, full_matrices=False)
            if len(s) == 0:
                continue
            u_all[ir] = u[:, 0]
            v_all[ir] = vh[0, :]
            s_all[ir] = s[0]

        ref_idx = int(np.argmax(np.abs(s_all)))
        u_ref = u_all[ref_idx].copy()
        for ir in range(n_radial):
            if np.dot(u_all[ir], u_ref) < 0.0:
                u_all[ir] *= -1.0
                v_all[ir] *= -1.0

        weights = np.abs(s_all)
        if float(weights.sum()) > 1e-15:
            ang1 = np.average(u_all, axis=0, weights=weights)
            ang2 = np.average(v_all, axis=0, weights=weights)
        else:
            ang1 = np.zeros(n_angle)
            ang2 = np.zeros(n_angle)

        max_abs = max(float(np.max(np.abs(ang1))), float(np.max(np.abs(ang2))), 1e-15)
        ang1 /= max_abs
        ang2 /= max_abs

        basis = np.outer(ang1, ang2)
        denom = float(np.dot(basis.ravel(), basis.ravel()))
        vm = np.zeros(n_radial, dtype=np.float64)
        if denom > 1e-15:
            for ir in range(n_radial):
                vm[ir] = float(np.dot(residual_remaining[ir].ravel(), basis.ravel()) / denom)
            residual_remaining -= basis[None, :, :] * vm[:, None, None]

        mode_ang1.append(ang1)
        mode_ang2.append(ang2)
        mode_radial.append(vm)

    recon = v_radial[:, None, None]
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        recon = recon + np.outer(ang1, ang2)[None, :, :] * vm[:, None, None]
    target_grid = v_radial[:, None, None] + residual_fit
    rms_error = float(np.sqrt(np.mean((target_grid - recon) ** 2)))

    # Fit B-splines
    inv_conv = 1.0 / ENERGY_CONVERSION_KJ_PER_EUP

    inv_dtheta = (n_knot_angular - 3) / 2.0
    t_angular = (cos_theta_grid + 1.0) * inv_dtheta + 1.0

    knot_spacing = float(knot_spacing_ang)
    t_radial_ang = r_values * 10.0 / knot_spacing
    rad_knot_vector = np.zeros(n_knot_radial + 4, dtype=np.float64)
    rad_knot_vector[4:-4] = np.arange(1, n_knot_radial - 3, dtype=np.float64)
    rad_knot_vector[-4:] = rad_knot_vector[-5]

    v0_knots = _fit_radial_bspline(t_radial_ang, v_radial, rad_knot_vector, smooth=0.01) * inv_conv
    n_unconstrained = max(1, min(n_knot_radial - 1, int(math.ceil(float(r_values[0]) * 10.0 / knot_spacing)) + 1))
    raw_short_range_core_kj_mol = float(np.max(energy_grid[0]))
    short_range_core_kj_mol = raw_short_range_core_kj_mol
    short_range_core_eup = raw_short_range_core_kj_mol * inv_conv
    control_min = float(np.min(energy_grid)) * inv_conv
    control_max = float(np.max(energy_grid)) * inv_conv
    v0_knots[:n_unconstrained] = short_range_core_eup
    param_parts = [v0_knots]
    ang1_knots_all = []
    ang2_knots_all = []
    vm_knots_all = []
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        ang1_knots = _fit_angular_bspline(t_angular, ang1, n_knot_angular, smooth=angular_smooth)
        ang2_knots = _fit_angular_bspline(t_angular, ang2, n_knot_angular, smooth=angular_smooth)
        vm_knots = _fit_radial_bspline(t_radial_ang, vm, rad_knot_vector, smooth=0.01) * inv_conv
        ang1_knots = np.clip(ang1_knots, -1.0, 1.0)
        ang2_knots = np.clip(ang2_knots, -1.0, 1.0)
        vm_knots[:n_unconstrained] = 0.0
        ang1_knots_all.append(ang1_knots)
        ang2_knots_all.append(ang2_knots)
        vm_knots_all.append(vm_knots)
        param_parts.extend([ang1_knots, ang2_knots, vm_knots])

    v0_knots[:] = np.clip(v0_knots, control_min, control_max)
    for vm_knots in vm_knots_all:
        vm_knots[:] = np.clip(vm_knots, control_min, control_max)

    attractive_control_count = int(np.count_nonzero(v0_knots < 0.0))
    for vm_knots in vm_knots_all:
        attractive_control_count += int(np.count_nonzero(vm_knots < 0.0))
    param_parts = [v0_knots]
    for ang1_knots, ang2_knots, vm_knots in zip(ang1_knots_all, ang2_knots_all, vm_knots_all):
        param_parts.extend([ang1_knots, ang2_knots, vm_knots])

    interaction_param = np.concatenate(param_parts)

    taper_width_ang = knot_spacing
    cutoff_ang = float(r_values[-1] * 10.0 + taper_width_ang)

    return {
        "v0_knots": v0_knots,
        "ang1_knots": np.asarray(ang1_knots_all),
        "ang2_knots": np.asarray(ang2_knots_all),
        "v_mode_knots": np.asarray(vm_knots_all),
        "interaction_param": interaction_param,
        "rms_error": rms_error,
        "v_radial_raw": v_radial,
        "short_range_core_kj_mol": short_range_core_kj_mol,
        "raw_short_range_core_kj_mol": raw_short_range_core_kj_mol,
        "short_range_core_rows": int(n_unconstrained),
        "short_range_core_source": "max_first_sampled_dry_martini_energy_expectation",
        "excluded_area_source": "none_full_resolved_dry_martini_sc_cgl_table",
        "excluded_area_nonnegative_rows": 0,
        "isotropic_background_source": (
            "none_full_resolved_dry_martini_sc_cgl_table"
        ),
        "attractive_control_source": (
            "retained_full_resolved_dry_martini_sc_cgl_table"
        ),
        "ang1_raw": np.asarray(mode_ang1),
        "ang2_raw": np.asarray(mode_ang2),
        "v_angular_raw": np.asarray(mode_radial),
        "r_values_nm": r_values,
        "cos_theta_grid": cos_theta_grid,
        "knot_spacing_ang": knot_spacing,
        "cutoff_ang": cutoff_ang,
        "taper_width_ang": taper_width_ang,
        "n_modes": int(len(mode_radial)),
        "n_radial": int(n_knot_radial),
        "n_angular": int(n_knot_angular),
        "angular_smooth": float(angular_smooth),
        "azimuthal_count": int(azimuthal_count),
        "azimuthal_average": (
            "energy_expectation"
            if average_temperature <= 0.0
            else (
                "tempered_boltzmann_free_energy"
                if abs(average_temperature - float(temperature)) > 1e-8
                else "boltzmann_free_energy"
            )
        ),
            "azimuthal_average_temperature_upside": float(average_temperature),
            "sidechain_bead_frame_count": int(len(sidechain_bead_frame_angles)),
            "cg_bead_frame_count": int(len(cg_bead_frame_angles)),
        }


def _fit_cg_lipid_sc_quadspline_from_dict(task: dict) -> tuple[int, str, dict]:
    task = dict(task)
    ri = int(task.pop("ri"))
    residue = str(task.pop("residue"))
    return ri, residue, _fit_cg_lipid_sc_quadspline(**task)


def _rotation_to_align_z_np(dir_vec: np.ndarray) -> np.ndarray:
    """Build rotation matrix that maps z-axis to dir_vec."""
    z_axis = dir_vec / np.linalg.norm(dir_vec)
    if abs(z_axis[0]) < 0.99:
        x_axis = np.cross([1.0, 0.0, 0.0], z_axis)
    else:
        x_axis = np.cross([0.0, 1.0, 0.0], z_axis)
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    return np.array([x_axis, y_axis, z_axis]).T


def _run_cgl_effective_lj_task(task: dict) -> tuple[str, dict]:
    target_type = str(task["target_type"])
    bead_types = task["bead_types"]
    pair_params = task["pair_params"]
    ref_nm = np.asarray(task["ref_nm"], dtype=np.float64)
    if ref_nm.ndim == 2:
        ref_nm = ref_nm[None, :, :]
    dir_array = task["dir_array"]
    bead_frame_angles = task["bead_frame_angles"]
    r_values = task["r_values"]
    r6 = r_values ** 6
    r12 = r_values ** 12

    def _resolve_self_params(bt: str) -> dict:
        p = pair_params.get((bt, bt))
        if p is None:
            raise RuntimeError(f"No self-interaction LJ params for type {bt}")
        return p

    bead_params = []
    for bt in bead_types:
        params = pair_params.get((bt, target_type)) or pair_params.get((target_type, bt))
        if params is None:
            p_self = _resolve_self_params(bt)
            p_tgt = _resolve_self_params(target_type)
            params = {
                "sigma_nm": (p_self["sigma_nm"] + p_tgt["sigma_nm"]) / 2.0,
                "epsilon_kj_mol": math.sqrt(p_self["epsilon_kj_mol"] * p_tgt["epsilon_kj_mol"]),
            }
        bead_params.append(params)

    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    avg_energy = np.zeros(len(r_values), dtype=np.float64)
    for ir, r in enumerate(r_values):
        energy_sum = 0.0
        sample_count = 0
        for dir_vec in dir_array:
            target_pos = float(r) * dir_vec
            for frame_angle in bead_frame_angles:
                rot = _rotation_about_axis_np(z_axis, float(frame_angle))
                for ref_conf in ref_nm:
                    framed_ref = (rot @ ref_conf.T).T
                    total_lj = 0.0
                    for bead_pos, params in zip(framed_ref, bead_params):
                        d = target_pos - bead_pos
                        dist = float(math.sqrt(float(np.dot(d, d))))
                        if dist < 0.001:
                            dist = 0.001
                        sig = params["sigma_nm"]
                        eps = params["epsilon_kj_mol"]
                        sr = sig / dist
                        sr2 = sr * sr
                        sr6 = sr2 * sr2 * sr2
                        total_lj += 4.0 * eps * (sr6 * sr6 - sr6)
                    energy_sum += total_lj
                    sample_count += 1
        avg_energy[ir] = energy_sum / max(sample_count, 1)

    attractive_mask = avg_energy < 0.0
    if attractive_mask.sum() >= 3:
        y = avg_energy[attractive_mask] * r12[attractive_mask]
        x = r6[attractive_mask]
        fit = np.polyfit(x, y, 1)
        B = -float(fit[0])
        A = float(fit[1])
        if A > 0.0 and B > 0.0:
            sigma_eff = (A / B) ** (1.0 / 6.0)
            epsilon_eff = B * B / (4.0 * A)
        else:
            imin = int(np.argmin(avg_energy))
            sigma_eff = r_values[imin] / (2.0 ** (1.0 / 6.0))
            epsilon_eff = max(-float(avg_energy[imin]), 0.01)
    else:
        imin = int(np.argmin(avg_energy))
        sigma_eff = r_values[imin] / (2.0 ** (1.0 / 6.0))
        epsilon_eff = 0.01

    uncapped_sigma_eff = max(0.1, min(float(sigma_eff), 5.0))
    return target_type, {
        "sigma_nm": uncapped_sigma_eff,
        "epsilon_kj_mol": max(0.01, min(float(epsilon_eff), 100.0)),
        "uncapped_sigma_nm": uncapped_sigma_eff,
        "orientation_count": len(dir_array),
        "bead_frame_count": int(len(bead_frame_angles)),
        "conformer_count": int(ref_nm.shape[0]),
    }


def _compute_cgl_effective_lj_params(
    ref_bead_positions_nm: np.ndarray,
    bead_types: list,
    pair_params: dict,
    r_min_nm: float = 0.2,
    r_max_nm: float = 3.0,
    n_radial: int = 100,
    n_orientations: int = 200,
    bead_frame_count: int = 1,
) -> dict:
    """Compute effective LJ (sigma_nm, epsilon_kj_mol) for CGL syntype with each target type.

    Orientation-averages the total LJ interaction between all 14 DOPC beads and a
    point particle of each target type, then fits effective LJ(12,6) parameters.

    Returns dict: target_type to dict(sigma_nm=float, epsilon_kj_mol=float)
    """
    ref_nm = _canonicalize_lipid_reference_ensemble_to_z(np.asarray(ref_bead_positions_nm, dtype=np.float64))
    n_beads = ref_nm.shape[1]
    if n_beads != len(bead_types):
        raise ValueError(f"Bead count mismatch: {n_beads} positions vs {len(bead_types)} types")

    all_types = sorted(set(k[0] for k in pair_params.keys()) | set(k[1] for k in pair_params.keys()))

    directions = _fibonacci_sphere(n_orientations)
    dir_array = np.asarray(directions, dtype=np.float64)
    bead_frame_angles = _bead_frame_angles(bead_frame_count)

    r_values = np.linspace(r_min_nm, r_max_nm, n_radial)
    tasks = [
        {
            "target_type": target_type,
            "bead_types": list(bead_types),
            "pair_params": pair_params,
            "ref_nm": ref_nm,
            "dir_array": dir_array,
            "bead_frame_angles": bead_frame_angles,
            "r_values": r_values,
        }
        for target_type in all_types
    ]
    return dict(_parallel_map_ordered("CGL effective LJ metadata", _run_cgl_effective_lj_task, tasks))


def _build_cg_lipid_tables(
    h5: h5py.File,
    pair_params: dict,
    sidechain_lib_path: Path,
    martinize_path: Path,
    forcefield_name: str,
    active_residue_names: list,
    ref_bead_positions_nm: np.ndarray | None = None,
    bead_types: list | None = None,
    bead_charges: list | None = None,
    bead_masses_g_mol: dict | None = None,
    lipids_itp_path: Path | None = None,
    r_min_nm: float = 0.30,
    r_max_nm: float = 0.70,
    r_count: int = 24,
    cos_theta_count: int = 13,
    azimuthal_count: int = 4,
    fluid_pair_pmf_upside_h5_path: Path | None = None,
    fluid_pair_pmf_frame_start_fraction: float = 0.5,
    force_match_upside_h5_path: Path | None = None,
    force_match_frame_start_fraction: float = 0.5,
    force_match_max_frames: int = 100,
    force_match_min_count: int = 4,
    force_match_mode: str = "radial",
    ibi_base_pair_h5_path: Path | None = None,
    ibi_target_upside_h5_path: Path | None = None,
    ibi_model_upside_h5_path: Path | None = None,
    ibi_target_frame_start_fraction: float = 0.5,
    ibi_model_frame_start_fraction: float = 0.5,
) -> None:
    """Build CG lipid pair and CG lipid-SC quadspline tables and store in HDF5."""
    if ref_bead_positions_nm is None:
        print("  cg_lipid_table: no reference bead positions provided, skipping")
        return

    if bead_types is None or bead_charges is None:
        raise ValueError(
            "bead_types and bead_charges must be provided; "
            "parse them from the ITP via martini_itp_reader.parse_dopc_from_itp()"
        )

    bead_charges = [float(q) for q in bead_charges]

    if lipids_itp_path is not None:
        _ensure_cg_bonds_angles(lipids_itp_path)

    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if ref_nm.shape == (14, 3):
        ref_ensemble_nm = _canonicalize_lipid_reference_ensemble_to_z(ref_nm)
    elif ref_nm.ndim == 3 and ref_nm.shape[1:] == (14, 3):
        ref_ensemble_nm = _canonicalize_lipid_reference_ensemble_to_z(ref_nm)
    else:
        raise ValueError(f"ref_bead_positions_nm must be (14, 3) or (n, 14, 3), got {ref_nm.shape}")

    print("\n=== CG Lipid Table Building ===")
    bead_mass_values = None
    if bead_masses_g_mol is not None:
        bead_mass_values = [bead_masses_g_mol[bt] for bt in bead_types]
    derived_params = _derive_dopc_cg_params_from_reference_ensemble(
        ref_bead_positions_nm=ref_ensemble_nm,
        bead_types=bead_types,
        pair_params=pair_params,
        bead_masses_g_mol=bead_mass_values,
        bonds=_CURRENT_CG_BONDS,
    )
    contact_nm = float(derived_params["contact_nm"])
    print(
        "  DOPC-derived CGL params: "
        f"contact={derived_params['contact_ang']:.3f} A, "
        f"orientation_length={derived_params['orientation_length_ang']:.3f} A, "
        f"orientation_mass={derived_params['orientation_mass_g_mol']:.3f} g/mol, "
        f"orientation_bond_fc={derived_params['orientation_bond_fc_eup_a2']:.3f} E_up/A^2"
    )

    # Resolution control: "coarse", "medium", or "fine" (default).
    _res = os.environ.get("UPSIDE_CG_LIPID_RESOLUTION", "fine").strip().lower()
    cg_knot_spacing_ang = 0.35
    cg_r_max_nm, cg_n_radial = _cg_lipid_pair_radial_support(ref_ensemble_nm, cg_knot_spacing_ang)
    sc_knot_spacing_ang = 1.4
    sc_fit_r_max_nm, sc_n_radial = _cg_lipid_target_radial_support(ref_ensemble_nm, sc_knot_spacing_ang)
    if _res == "coarse":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 2, 8, 5, 2
    elif _res == "medium":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 4, 12, 7, 2
    else:
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 4, 16, 9, 4
    cg_r_count = max(_cg_r, cg_n_radial)
    cg_cos_theta_count = min(cos_theta_count, _cg_ct)
    sc_r_count = max(min(r_count, _sc_r), sc_n_radial)
    sc_cos_theta_count = min(cos_theta_count, _sc_ct)
    sc_azimuthal_count = min(azimuthal_count, _sc_az)
    cg_bead_frame_count = _bead_frame_count("CGL", 8)
    sc_bead_frame_count = _bead_frame_count("SC", 1)
    print(f"  CG lipid resolution: {_res} "
          f"(CG: {cg_r_count}r x {cg_cos_theta_count}^2 theta x {_cg_az}^2 phi x {cg_bead_frame_count}^2 bead-frame, "
          f"SC: {sc_r_count}r x {sc_cos_theta_count}^2 theta x {sc_azimuthal_count}^2 phi "
          f"x {sc_bead_frame_count} SC bead-frame x {cg_bead_frame_count} CGL bead-frame)")
    print(
        "  CGL-CGL radial support: "
        f"r_max={cg_r_max_nm:.3f} nm, n_radial={cg_n_radial}, "
        "source=2*max_DOPC_bead_radius_plus_dry_MARTINI_cutoff"
    )
    print(
        "  CG-SC radial support: "
        f"r_max={sc_fit_r_max_nm:.3f} nm, n_radial={sc_n_radial}, "
        "source=max_DOPC_bead_radius_plus_dry_MARTINI_cutoff"
    )

    print(f"  CGL tables use direct rotated ensemble geometry ({ref_ensemble_nm.shape[0]} conformer(s))")

    # CG-CG directional spline from direct rotated DOPC bead geometries.
    cg_n_angular = min(cg_cos_theta_count + 2, 15)
    cg_fit_smooth = _float_env("UPSIDE_CGL_PAIR_FIT_SMOOTH", 0.1)
    result_cg = _fit_cg_lipid_quadspline(
        ref_bead_positions_nm=ref_ensemble_nm,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        r_min_nm=0.035,
        r_max_nm=cg_r_max_nm,
        r_count=cg_r_count,
        cos_theta_count=cg_cos_theta_count,
        azimuthal_count=_cg_az,
        bead_frame_count=cg_bead_frame_count,
        dist_min_nm=1e-6,
        knot_spacing_ang=cg_knot_spacing_ang,
        n_modes=4,
        n_knot_radial=cg_n_radial,
        n_knot_angular=cg_n_angular,
        cg_smooth=cg_fit_smooth,
        temperature=DEFAULT_PRODUCTION_TEMP_UPSIDE,
        average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        fluid_pair_pmf_upside_h5_path=fluid_pair_pmf_upside_h5_path,
        fluid_pair_pmf_frame_start_fraction=fluid_pair_pmf_frame_start_fraction,
        fluid_pair_pmf_core_min_nm=contact_nm,
        force_match_upside_h5_path=force_match_upside_h5_path,
        force_match_frame_start_fraction=force_match_frame_start_fraction,
        force_match_max_frames=force_match_max_frames,
        force_match_min_count=force_match_min_count,
        force_match_mode=force_match_mode,
        ibi_base_pair_h5_path=ibi_base_pair_h5_path,
        ibi_target_upside_h5_path=ibi_target_upside_h5_path,
        ibi_model_upside_h5_path=ibi_model_upside_h5_path,
        ibi_target_frame_start_fraction=ibi_target_frame_start_fraction,
        ibi_model_frame_start_fraction=ibi_model_frame_start_fraction,
    )
    print(
        "  CG-CG: full tensor table "
        f"{result_cg['n_radial']}r x {result_cg['n_angular']}^2 angular, "
        f"max|E| = {float(np.max(np.abs(result_cg['energy_grid_raw']))):.3f} kJ/mol"
    )

    pair_relaxation_result = None
    if lipids_itp_path is not None and _truthy_env("UPSIDE_CGL_PAIR_RELAX_BASE"):
        pair_relaxation_face_cos_min = _float_env("UPSIDE_CGL_PAIR_RELAX_FACE_COS_MIN", 0.5)
        pair_relaxation_radial_cutoff_nm = _float_env("UPSIDE_CGL_PAIR_RELAX_RADIAL_CUTOFF_NM", 2.5)
        result_cg, pair_relaxation_result = _apply_pair_conditioned_tail_relaxation_to_cg_result(
            result_cg,
            ref_ensemble_nm,
            bead_types=list(bead_types),
            bead_charges=list(bead_charges),
            pair_params=pair_params,
            lipids_itp_path=lipids_itp_path,
            face_cos_min=pair_relaxation_face_cos_min,
            radial_cutoff_nm=pair_relaxation_radial_cutoff_nm,
        )
        print(
            "  CG-CG pair-relaxation base correction: "
            f"face cos>={pair_relaxation_face_cos_min:.3f}, "
            f"r<={pair_relaxation_radial_cutoff_nm:.3f} nm, "
            f"min delta={float(np.min(pair_relaxation_result['correction_grid_kj_mol'])):.3f} kJ/mol"
        )

    compaction_states = None
    compaction_result = None
    if lipids_itp_path is not None:
        from martini_itp_reader import parse_dopc_from_itp

        compaction_pool_conformers = _positive_int_env(
            "UPSIDE_CGL_COMPACTION_POOL_CONFORMERS", 32
        )
        compaction_burnin_steps = _positive_int_env(
            "UPSIDE_CGL_COMPACTION_BURNIN_STEPS", 20000
        )
        compaction_steps_per_conf = _positive_int_env(
            "UPSIDE_CGL_COMPACTION_STEPS_PER_CONFORMER", 500
        )
        compaction_representatives = _positive_int_env(
            "UPSIDE_CGL_COMPACTION_STATE_REPRESENTATIVES", 2
        )
        compaction_self_bins = _positive_int_env(
            "UPSIDE_CGL_COMPACTION_SELF_BINS", 12
        )
        dopc = parse_dopc_from_itp(lipids_itp_path)
        compaction_refs_nm = sample_isolated_dopc_bonded_conformers(
            dopc,
            lipids_itp_path=lipids_itp_path,
            pair_params=pair_params,
            conformer_count=compaction_pool_conformers,
            temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
            seed=1777,
            mc_burnin_steps=compaction_burnin_steps,
            mc_steps_per_conformer=compaction_steps_per_conf,
        )
        compaction_values_ang = _dopc_tail_extension_series_ang(compaction_refs_nm)
        compaction_states = _select_compaction_state_representatives(
            compaction_refs_nm,
            compaction_values_ang,
            representative_count=compaction_representatives,
        )
        compaction_self = _fit_compaction_self_pmf(
            compaction_values_ang,
            temperature_upside=DEFAULT_PRODUCTION_TEMP_UPSIDE,
            n_bin=compaction_self_bins,
            smooth=0.01,
        )
        state_grid_ee = _sample_cg_pair_energy_grid(
            compaction_states["extended_refs_nm"],
            compaction_states["extended_refs_nm"],
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            r_values_nm=result_cg["r_values_nm"],
            cos_theta_grid=result_cg["cos_theta_grid"],
            azimuthal_count=_cg_az,
            bead_frame_count=cg_bead_frame_count,
            dist_min_nm=1.0e-6,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )
        state_grid_cc = _sample_cg_pair_energy_grid(
            compaction_states["compact_refs_nm"],
            compaction_states["compact_refs_nm"],
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            r_values_nm=result_cg["r_values_nm"],
            cos_theta_grid=result_cg["cos_theta_grid"],
            azimuthal_count=_cg_az,
            bead_frame_count=cg_bead_frame_count,
            dist_min_nm=1.0e-6,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )
        state_grid_ec_forward = _sample_cg_pair_energy_grid(
            compaction_states["extended_refs_nm"],
            compaction_states["compact_refs_nm"],
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            r_values_nm=result_cg["r_values_nm"],
            cos_theta_grid=result_cg["cos_theta_grid"],
            azimuthal_count=_cg_az,
            bead_frame_count=cg_bead_frame_count,
            dist_min_nm=1.0e-6,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )
        state_grid_ec_reverse = _sample_cg_pair_energy_grid(
            compaction_states["compact_refs_nm"],
            compaction_states["extended_refs_nm"],
            bead_types=bead_types,
            bead_charges=bead_charges,
            pair_params=pair_params,
            r_values_nm=result_cg["r_values_nm"],
            cos_theta_grid=result_cg["cos_theta_grid"],
            azimuthal_count=_cg_az,
            bead_frame_count=cg_bead_frame_count,
            dist_min_nm=1.0e-6,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )
        state_grid_ec = 0.5 * (
            state_grid_ec_forward + np.swapaxes(state_grid_ec_reverse, 1, 2)
        )
        compact_probability = float(compaction_states["compact_probability"])
        extended_probability = 1.0 - compact_probability
        state_grid_avg = (
            (extended_probability * extended_probability) * state_grid_ee
            + (2.0 * compact_probability * extended_probability) * state_grid_ec
            + (compact_probability * compact_probability) * state_grid_cc
        )
        base_grid_kj_mol = np.asarray(result_cg["energy_grid_raw"], dtype=np.float64)
        reference_energy_kj_mol = (
            float(result_cg["reference_energy_eup"]) * ENERGY_CONVERSION_KJ_PER_EUP
        )
        control_kbt_kj_mol = DEFAULT_PRODUCTION_TEMP_UPSIDE * ENERGY_CONVERSION_KJ_PER_EUP
        if control_kbt_kj_mol <= 0.0:
            raise ValueError("CGL compaction correction requires positive control temperature")

        def to_log1p_control(grid_kj_mol: np.ndarray) -> np.ndarray:
            reduced = (np.asarray(grid_kj_mol, dtype=np.float64) - reference_energy_kj_mol) / control_kbt_kj_mol
            return np.log1p(np.maximum(reduced, 0.0))

        base_control_grid = to_log1p_control(base_grid_kj_mol)
        correction_ee_control = to_log1p_control(state_grid_ee) - base_control_grid
        correction_ec_control = to_log1p_control(state_grid_ec) - base_control_grid
        correction_cc_control = to_log1p_control(state_grid_cc) - base_control_grid
        compaction_face_mask = None
        compaction_face_cos_min = os.environ.get("UPSIDE_CGL_COMPACTION_FACE_COS_MIN", "").strip()
        if compaction_face_cos_min:
            face_cos_min = float(compaction_face_cos_min)
            radial_cutoff_nm = _float_env("UPSIDE_CGL_COMPACTION_FACE_RADIAL_CUTOFF_NM", 2.5)
            compaction_face_mask = _build_cross_leaflet_face_mask(
                result_cg["r_values_nm"],
                result_cg["cos_theta_grid"],
                face_cos_min,
                radial_cutoff_nm,
            )
            correction_ee_control = np.where(compaction_face_mask, correction_ee_control, 0.0)
            correction_ec_control = np.where(compaction_face_mask, correction_ec_control, 0.0)
            correction_cc_control = np.where(compaction_face_mask, correction_cc_control, 0.0)

        def fit_correction_tensor(correction_grid_control: np.ndarray) -> np.ndarray:
            compaction_fit_smooth = _float_env("UPSIDE_CGL_COMPACTION_FIT_SMOOTH", cg_fit_smooth)
            tensor = _fit_radial_angular_angular_tensor_bspline(
                result_cg["r_values_nm"],
                result_cg["cos_theta_grid"],
                correction_grid_control,
                n_knot_radial=result_cg["n_radial"],
                n_knot_angular=result_cg["n_angular"],
                knot_spacing_ang=result_cg["knot_spacing_ang"],
                energy_conversion=1.0,
                smooth=compaction_fit_smooth,
            )
            corr_min = float(np.min(correction_grid_control))
            corr_max = float(np.max(correction_grid_control))
            return np.clip(tensor, corr_min, corr_max)

        compaction_result = {
            "self": compaction_self,
            "compact_center_ang": float(compaction_states["compact_center_ang"]),
            "extended_center_ang": float(compaction_states["extended_center_ang"]),
            "compact_probability": compact_probability,
            "pool_size": int(compaction_states["pool_size"]),
            "compact_pool_size": int(compaction_states["compact_pool_size"]),
            "extended_pool_size": int(compaction_states["extended_pool_size"]),
            "representative_compact_count": int(compaction_states["compact_refs_nm"].shape[0]),
            "representative_extended_count": int(compaction_states["extended_refs_nm"].shape[0]),
            "compaction_min_ang": float(compaction_states["compaction_min_ang"]),
            "compaction_max_ang": float(compaction_states["compaction_max_ang"]),
            "compaction_p05_ang": float(compaction_states["compaction_p05_ang"]),
            "compaction_p95_ang": float(compaction_states["compaction_p95_ang"]),
            "delta_extended_extended": fit_correction_tensor(correction_ee_control).reshape(-1),
            "delta_extended_compact": fit_correction_tensor(correction_ec_control).reshape(-1),
            "delta_compact_compact": fit_correction_tensor(correction_cc_control).reshape(-1),
            "grid_extended_extended_kj_mol": state_grid_ee,
            "grid_extended_compact_kj_mol": state_grid_ec,
            "grid_compact_compact_kj_mol": state_grid_cc,
            "grid_average_kj_mol": state_grid_avg,
            "fit_smooth": float(compaction_fit_smooth),
        }
        if compaction_face_mask is not None:
            compaction_result["face_mask"] = compaction_face_mask.astype(np.int8)
            compaction_result["face_cos_min"] = float(face_cos_min)
            compaction_result["radial_cutoff_nm"] = float(radial_cutoff_nm)
        print(
            "  CGL compaction state: "
            f"pool={compaction_result['pool_size']}, compact p={compact_probability:.3f}, "
            f"compact center={compaction_result['compact_center_ang']:.3f} A, "
            f"extended center={compaction_result['extended_center_ang']:.3f} A"
        )

    # CG-SC quadspline
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
        martinize_path, forcefield_name
    )

    residues = [r for r in active_residue_names
                if r in residue_map and residue_map[r] and r in orientation_map]
    if _truthy_env("UPSIDE_CG_LIPID_SKIP_SC"):
        print("  cg_lipid_sc: skipping SC-CGL table generation for CGL-only screening")
        residues = []
    sc_n_modes = 0
    sc_n_angular = min(sc_cos_theta_count + 2, 15)
    sc_angular_smooth = 0.1
    sc_taper_width_ang = sc_knot_spacing_ang
    sc_cutoff_ang = float(sc_fit_r_max_nm * 10.0 + sc_taper_width_ang)
    if sc_n_modes == 0:
        sc_n_param = sc_n_radial * sc_n_angular * sc_n_angular
    else:
        sc_n_param = sc_n_radial + sc_n_modes * (2 * sc_n_angular + sc_n_radial)
    sc_rms_error = np.zeros(len(residues), dtype=np.float32)
    sc_short_range_core = np.zeros(len(residues), dtype=np.float32)
    sc_short_range_core_rows = np.zeros(len(residues), dtype=np.int32)
    sc_reference_energy = np.zeros(len(residues), dtype=np.float32)
    sc_compaction_delta_extended = None
    sc_compaction_delta_compact = None
    sc_compaction_extended_grid = None
    sc_compaction_compact_grid = None
    first_sc_result = None

    if not residues:
        print("  cg_lipid_sc: no active residues with sidechains, skipping")
        n_sc_types = 0
        interaction_param_sc = np.zeros((0, 1, sc_n_param), dtype=np.float64)
        sc_residue_names = []
    else:
        cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
        cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)

        n_sc_types = len(residues)
        interaction_param_sc = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float64)
        sc_residue_names = []
        base_sc_results = {}

        sc_fit_tasks = []
        for ri, residue in enumerate(residues):
            sc_bead_types = residue_map[residue]
            sc_bead_charges = [infer_charge_from_atomtype(bt) for bt in sc_bead_types]
            orientation = orientation_map[residue]
            sc_positions_by_rotamer = _expand_rotamer_sidechain_positions(
                orientation,
                residue,
                np.asarray(martini_sidechain_offsets_nm[residue], dtype=np.float64),
            )

            sc_fit_tasks.append(
                {
                    "ri": ri,
                    "residue": residue,
                    "ref_bead_positions_nm": ref_ensemble_nm,
                    "cg_bead_types": list(bead_types),
                    "cg_bead_charges": list(bead_charges),
                    "target_type": "CGL",
                    "target_charge": 0.0,
                    "rotamer_bead_positions_nm": sc_positions_by_rotamer,
                    "rotamer_weights": orientation["weight"],
                    "sc_bead_types": list(sc_bead_types),
                    "sc_bead_charges": list(sc_bead_charges),
                    "pair_params": pair_params,
                    "cb_anchor_nm": cb_anchor_nm,
                    "cb_vector_unit": cb_vector_unit,
                    "r_min_nm": r_min_nm,
                    "r_max_nm": sc_fit_r_max_nm,
                    "r_count": sc_r_count,
                    "cos_theta_count": sc_cos_theta_count,
                    "azimuthal_count": sc_azimuthal_count,
                    "sidechain_bead_frame_count": sc_bead_frame_count,
                    "cg_bead_frame_count": cg_bead_frame_count,
                    "n_modes": sc_n_modes,
                    "n_knot_radial": sc_n_radial,
                    "n_knot_angular": sc_n_angular,
                    "angular_smooth": sc_angular_smooth,
                    "knot_spacing_ang": sc_knot_spacing_ang,
                    "temperature": DEFAULT_PRODUCTION_TEMP_UPSIDE,
                    "average_temperature": DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
                }
            )

        for ri, residue, result_sc in _parallel_map_ordered(
            "CG-SC table", _fit_cg_lipid_sc_quadspline_from_dict, sc_fit_tasks
        ):
            if first_sc_result is None:
                first_sc_result = result_sc
            interaction_param_sc[ri, 0, :] = result_sc["interaction_param"]
            sc_rms_error[ri] = np.float32(result_sc["rms_error"])
            sc_short_range_core[ri] = np.float32(result_sc["short_range_core_kj_mol"])
            sc_short_range_core_rows[ri] = np.int32(result_sc["short_range_core_rows"])
            sc_reference_energy[ri] = np.float32(result_sc.get("reference_energy_eup", 0.0))
            sc_cutoff_ang = min(sc_cutoff_ang, float(result_sc["cutoff_ang"]))
            sc_residue_names.append(residue)
            base_sc_results[ri] = result_sc
            print(f"  CG-SC({residue}): RMS error = {result_sc['rms_error']:.4f} kJ/mol, "
                  f"modes = {result_sc['n_modes']}")

        if compaction_states is not None and base_sc_results:
            sc_compaction_delta_extended = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
            sc_compaction_delta_compact = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float32)
            sc_compaction_extended_grid = np.zeros(
                (n_sc_types, sc_r_count, sc_cos_theta_count, sc_cos_theta_count),
                dtype=np.float32,
            )
            sc_compaction_compact_grid = np.zeros_like(sc_compaction_extended_grid, dtype=np.float32)
            sc_extended_tasks = []
            sc_compact_tasks = []
            for task in sc_fit_tasks:
                sc_extended_tasks.append({
                    **task,
                    "ref_bead_positions_nm": np.asarray(
                        compaction_states["extended_refs_nm"],
                        dtype=np.float64,
                    ),
                })
                sc_compact_tasks.append({
                    **task,
                    "ref_bead_positions_nm": np.asarray(
                        compaction_states["compact_refs_nm"],
                        dtype=np.float64,
                    ),
                })

            extended_sc_results = {
                ri: result_sc
                for ri, _residue, result_sc in _parallel_map_ordered(
                    "CG-SC extended-state table",
                    _fit_cg_lipid_sc_quadspline_from_dict,
                    sc_extended_tasks,
                )
            }
            compact_sc_results = {
                ri: result_sc
                for ri, _residue, result_sc in _parallel_map_ordered(
                    "CG-SC compacted-state table",
                    _fit_cg_lipid_sc_quadspline_from_dict,
                    sc_compact_tasks,
                )
            }
            for ri in range(n_sc_types):
                base_sc = base_sc_results[ri]
                delta = _fit_single_cgl_state_delta_full_tensor(
                    base_energy_grid_kj_mol=np.asarray(base_sc["energy_grid_raw"], dtype=np.float64),
                    extended_energy_grid_kj_mol=np.asarray(
                        extended_sc_results[ri]["energy_grid_raw"],
                        dtype=np.float64,
                    ),
                    compact_energy_grid_kj_mol=np.asarray(
                        compact_sc_results[ri]["energy_grid_raw"],
                        dtype=np.float64,
                    ),
                    compressed_energy_grid_kj_mol=None,
                    reference_energy_eup=float(base_sc["reference_energy_eup"]),
                    base_interaction_param=np.asarray(base_sc["interaction_param"], dtype=np.float64),
                    temperature_upside=float(DEFAULT_PRODUCTION_TEMP_UPSIDE),
                    r_values_nm=np.asarray(base_sc["r_values_nm"], dtype=np.float64),
                    cos_theta_grid=np.asarray(base_sc["cos_theta_grid"], dtype=np.float64),
                    n_knot_radial=int(base_sc["n_radial"]),
                    n_knot_angular=int(base_sc["n_angular"]),
                    knot_spacing_ang=float(base_sc["knot_spacing_ang"]),
                    smooth=float(sc_angular_smooth),
                )
                sc_compaction_delta_extended[ri, 0, :] = np.asarray(
                    delta["delta_extended"],
                    dtype=np.float32,
                )
                sc_compaction_delta_compact[ri, 0, :] = np.asarray(
                    delta["delta_compact"],
                    dtype=np.float32,
                )
                sc_compaction_extended_grid[ri] = np.asarray(
                    delta["grid_extended_kj_mol"],
                    dtype=np.float32,
                )
                sc_compaction_compact_grid[ri] = np.asarray(
                    delta["grid_compact_kj_mol"],
                    dtype=np.float32,
                )

    # Store in HDF5
    cg_grp = h5.create_group("cg_lipid_table")
    cg_grp.attrs["bead_charge_source"] = "dry_martini_v2.1_lipids.itp:DOPC_atoms"
    cg_grp.attrs["lipid_net_charge"] = np.float32(sum(bead_charges))
    cg_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_grp.attrs["conformer_count"] = np.int32(ref_ensemble_nm.shape[0])
    cg_grp.attrs["conformer_source"] = "dynamic_vector_cgl_reference_ensemble"
    cg_grp.create_dataset("bead_charges", data=np.asarray(bead_charges, dtype=np.float32))
    cg_grp.create_dataset("ref_bead_positions_nm", data=ref_ensemble_nm.astype(np.float32))
    _write_cg_derived_attrs(cg_grp, derived_params)

    # CG-CG pair
    cg_pair_grp = cg_grp.create_group("cg_lipid_pair")
    if compaction_result is not None:
        cg_pair_correction = "source_derived_compaction_state"
    elif int(result_cg.get("force_match_enabled", 0)):
        cg_pair_correction = (
            "projected_generalized_force_matching"
            if str(result_cg.get("force_match_mode", "")) == "generalized"
            else "projected_radial_force_matching"
        )
    elif int(result_cg.get("ibi_enabled", 0)):
        cg_pair_correction = (
            "iterative_ibi_sampled_bins"
            if result_cg.get("ibi_base_pair_h5_path")
            else "fluid_pmf_ibi_sampled_bins"
        )
    elif int(result_cg.get("fluid_pmf_pair_sample_count", 0)) > 0:
        cg_pair_correction = "fluid_pmf_sampled_bins"
    else:
        cg_pair_correction = "none"
    _write_common_table_contract_attrs(
        cg_pair_grp,
        table_family="CGL-CGL",
        source_object="dopc_cgl_constituent_bead_ensemble",
        target_object="dopc_cgl_constituent_bead_ensemble",
        projection_ensemble="cgl_conformer_pairs_x_two_cgl_orientations_x_bead_frames",
        runtime_representation="log1p_reduced_tensor_bspline",
        correction_layer=cg_pair_correction,
    )
    pair_param = result_cg["interaction_param"].astype(np.float32)
    cg_pair_grp.create_dataset(
        "interaction_param",
        data=pair_param.reshape(1, 1, pair_param.size),
    )
    cg_pair_grp.create_dataset(
        "reference_energy_eup",
        data=np.asarray([[result_cg["reference_energy_eup"]]], dtype=np.float32),
    )
    cg_pair_grp.attrs["n_cg_types"] = 1
    cg_pair_grp.attrs["rms_error_kj_mol"] = np.float32(result_cg["rms_error"])
    cg_pair_grp.attrs["schema"] = result_cg["schema"]
    cg_pair_grp.attrs["radial_mode"] = "full_tensor"
    cg_pair_grp.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
    cg_pair_grp.attrs["energy_transform"] = result_cg["energy_transform"]
    cg_pair_grp.attrs["log1p_reduced_transform"] = np.int32(1)
    cg_pair_grp.attrs["boltzmann_temperature_upside"] = np.float32(DEFAULT_PRODUCTION_TEMP_UPSIDE)
    cg_pair_grp.attrs["spline_control_quantity"] = result_cg["spline_control_quantity"]
    cg_pair_grp.attrs["sample_dist_min_nm"] = np.float32(result_cg["sample_dist_min_nm"])
    cg_pair_grp.attrs["sample_dist_min_source"] = "numerical_zero_guard_only"
    cg_pair_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_pair_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_pair_grp.attrs["fit_r_min_nm"] = np.float32(result_cg["r_values_nm"][0])
    cg_pair_grp.attrs["fit_r_max_nm"] = np.float32(result_cg["r_values_nm"][-1])
    cg_pair_grp.attrs["radial_support_source"] = "2*max_dopc_bead_radius_plus_dry_martini_cutoff"
    cg_pair_grp.attrs["n_modes"] = np.int32(result_cg["n_modes"])
    cg_pair_grp.attrs["n_radial"] = np.int32(result_cg["n_radial"])
    cg_pair_grp.attrs["n_angular"] = np.int32(result_cg["n_angular"])
    cg_pair_grp.attrs["fit_smooth"] = np.float32(result_cg["fit_smooth"])
    cg_pair_grp.attrs["azimuthal_count"] = np.int32(result_cg["azimuthal_count"])
    cg_pair_grp.attrs["cgl_bead_frame_count"] = np.int32(result_cg["bead_frame_count"])
    cg_pair_grp.attrs["conformer_count"] = np.int32(result_cg["conformer_count"])
    cg_pair_grp.attrs["conformer_pair_count"] = np.int32(result_cg["conformer_pair_count"])
    cg_pair_grp.attrs["orientation_sampling"] = "both_cgl_direction_vectors"
    cg_pair_grp.attrs["knot_spacing_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["cutoff_ang"] = np.float32(result_cg["cutoff_ang"])
    cg_pair_grp.attrs["taper_width_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["azimuthal_average"] = result_cg["azimuthal_average"]
    cg_pair_grp.attrs["azimuthal_average_temperature_upside"] = np.float32(
        result_cg["azimuthal_average_temperature_upside"]
    )
    cg_pair_grp.attrs["isotropic_background_source"] = result_cg["isotropic_background_source"]
    cg_pair_grp.attrs["isotropic_background_min_kj_mol"] = np.float32(
        float(np.min(result_cg["attractive_radial_background_kj_mol"]))
    )
    cg_pair_grp.attrs["excluded_area_source"] = result_cg["excluded_area_source"]
    cg_pair_grp.attrs["excluded_area_nonnegative_rows"] = np.int32(result_cg["excluded_area_nonnegative_rows"])
    cg_pair_grp.attrs["attractive_control_source"] = result_cg["attractive_control_source"]
    cg_pair_grp.attrs["attractive_control_count"] = np.int32(result_cg["attractive_control_count"])
    cg_pair_grp.attrs["unresolved_core_source"] = result_cg["core_boundary_source"]
    cg_pair_grp.attrs["unresolved_core_boundary_row"] = np.int32(result_cg["core_boundary_row"])
    cg_pair_grp.attrs["unresolved_core_rows"] = np.int32(result_cg["unresolved_core_rows"])
    cg_pair_grp.attrs["unresolved_core_energy_kj_mol"] = np.float32(result_cg["unresolved_core_energy_kj_mol"])
    cg_pair_grp.attrs["fluid_pmf_pair_sample_count"] = np.int64(result_cg.get("fluid_pmf_pair_sample_count", 0))
    cg_pair_grp.attrs["fluid_pmf_occupied_radial_bins"] = np.int32(result_cg.get("fluid_pmf_occupied_radial_bins", 0))
    cg_pair_grp.attrs["fluid_pmf_occupied_tensor_bins"] = np.int32(result_cg.get("fluid_pmf_occupied_tensor_bins", 0))
    cg_pair_grp.attrs["fluid_pmf_frame_start_fraction"] = np.float32(
        result_cg.get("fluid_pmf_frame_start_fraction", 0.0)
    )
    cg_pair_grp.attrs["fluid_pmf_core_min_nm"] = np.float32(result_cg.get("fluid_pmf_core_min_nm", 0.0))
    cg_pair_grp.attrs["force_match_enabled"] = np.int32(result_cg.get("force_match_enabled", 0))
    cg_pair_grp.attrs["force_match_mode"] = result_cg.get("force_match_mode", "")
    cg_pair_grp.attrs["force_match_projection"] = (
        "radial_and_angular_generalized_forces_on_cgl_pair_coordinates"
        if result_cg.get("force_match_mode", "") == "generalized"
        else "radial_generalized_force_on_cgl_pair_coordinate"
    )
    cg_pair_grp.attrs["force_match_boundary_condition"] = (
        "connected_component_high_r_anchor_keeps_base_energy_constant"
        if result_cg.get("force_match_mode", "") == "generalized"
        else "segment_high_r_keeps_base_energy_constant"
    )
    cg_pair_grp.attrs["force_match_unsampled_bin_policy"] = "retain_existing_table_value"
    cg_pair_grp.attrs["force_match_frame_count"] = np.int32(result_cg.get("force_match_frame_count", 0))
    cg_pair_grp.attrs["force_match_pair_force_sample_count"] = np.int64(
        result_cg.get("force_match_pair_force_sample_count", 0)
    )
    cg_pair_grp.attrs["force_match_sampled_tensor_bins"] = np.int32(
        result_cg.get("force_match_sampled_tensor_bins", 0)
    )
    cg_pair_grp.attrs["force_match_updated_tensor_bins"] = np.int32(
        result_cg.get("force_match_updated_tensor_bins", 0)
    )
    cg_pair_grp.attrs["force_match_frame_start_fraction"] = np.float32(
        result_cg.get("force_match_frame_start_fraction", 0.0)
    )
    cg_pair_grp.attrs["force_match_max_frames"] = np.int32(result_cg.get("force_match_max_frames", 0))
    cg_pair_grp.attrs["force_match_min_count"] = np.int32(result_cg.get("force_match_min_count", 0))
    cg_pair_grp.attrs["force_match_lsq_variable_count"] = np.int32(
        result_cg.get("force_match_lsq_variable_count", 0)
    )
    cg_pair_grp.attrs["force_match_lsq_equation_count"] = np.int32(
        result_cg.get("force_match_lsq_equation_count", 0)
    )
    cg_pair_grp.attrs["force_match_lsq_residual_rms_kj_mol"] = np.float32(
        result_cg.get("force_match_lsq_residual_rms_kj_mol", 0.0)
    )
    cg_pair_grp.attrs["force_match_lsq_solver"] = result_cg.get("force_match_lsq_solver", "")
    cg_pair_grp.attrs["ibi_enabled"] = np.int32(result_cg.get("ibi_enabled", 0))
    cg_pair_grp.attrs["ibi_update_rule"] = "U_new=U_old+kBT*ln(P_model/P_target)_sampled_bins_only"
    cg_pair_grp.attrs["ibi_unsampled_bin_policy"] = "retain_existing_table_value"
    cg_pair_grp.attrs["ibi_base_pair_h5_path"] = result_cg.get("ibi_base_pair_h5_path", "")
    cg_pair_grp.attrs["ibi_target_pair_sample_count"] = np.int64(result_cg.get("ibi_target_pair_sample_count", 0))
    cg_pair_grp.attrs["ibi_model_pair_sample_count"] = np.int64(result_cg.get("ibi_model_pair_sample_count", 0))
    cg_pair_grp.attrs["ibi_target_occupied_tensor_bins"] = np.int32(result_cg.get("ibi_target_occupied_tensor_bins", 0))
    cg_pair_grp.attrs["ibi_model_occupied_tensor_bins"] = np.int32(result_cg.get("ibi_model_occupied_tensor_bins", 0))
    cg_pair_grp.attrs["ibi_sampled_tensor_bins"] = np.int32(result_cg.get("ibi_sampled_tensor_bins", 0))
    cg_pair_grp.attrs["ibi_target_frame_start_fraction"] = np.float32(
        result_cg.get("ibi_target_frame_start_fraction", 0.0)
    )
    cg_pair_grp.attrs["ibi_model_frame_start_fraction"] = np.float32(
        result_cg.get("ibi_model_frame_start_fraction", 0.0)
    )
    cg_pair_grp.attrs["ibi_step_size"] = np.float32(result_cg.get("ibi_step_size", 0.0))
    cg_pair_grp.attrs["ibi_raw_correction_min_kj_mol"] = np.float32(
        result_cg.get("ibi_raw_correction_min_kj_mol", 0.0)
    )
    cg_pair_grp.attrs["ibi_raw_correction_max_kj_mol"] = np.float32(
        result_cg.get("ibi_raw_correction_max_kj_mol", 0.0)
    )
    cg_pair_grp.attrs["ibi_raw_correction_mean_kj_mol"] = np.float32(
        result_cg.get("ibi_raw_correction_mean_kj_mol", 0.0)
    )
    cg_pair_grp.attrs["ibi_correction_min_kj_mol"] = np.float32(result_cg.get("ibi_correction_min_kj_mol", 0.0))
    cg_pair_grp.attrs["ibi_correction_max_kj_mol"] = np.float32(result_cg.get("ibi_correction_max_kj_mol", 0.0))
    cg_pair_grp.attrs["ibi_correction_mean_kj_mol"] = np.float32(result_cg.get("ibi_correction_mean_kj_mol", 0.0))
    if compaction_result is not None:
        cg_pair_grp.attrs["compaction_pair_correction"] = (
            "isolated_source_two_state_log1p_control_delta_relative_to_base_pair_table"
        )
        cg_pair_grp.attrs["compact_state_center_ang"] = np.float32(
            compaction_result["compact_center_ang"]
        )
        cg_pair_grp.attrs["extended_state_center_ang"] = np.float32(
            compaction_result["extended_center_ang"]
        )
    if "pair_relaxation_correction_source" in result_cg:
        cg_pair_grp.attrs["pair_relaxation_correction_source"] = result_cg["pair_relaxation_correction_source"]
        cg_pair_grp.attrs["pair_relaxation_energy_basis"] = result_cg.get(
            "pair_relaxation_energy_basis",
            "intermolecular_energy_kj_mol",
        )
        cg_pair_grp.attrs["pair_relaxation_face_cos_min"] = np.float32(
            result_cg["pair_relaxation_face_cos_min"]
        )
        cg_pair_grp.attrs["pair_relaxation_radial_cutoff_nm"] = np.float32(
            result_cg["pair_relaxation_radial_cutoff_nm"]
        )
    _write_cg_derived_attrs(cg_pair_grp, derived_params)
    cg_pair_grp.create_dataset("energy_grid_raw_kj_mol", data=result_cg["energy_grid_raw"].astype(np.float32))
    if "pair_relaxation_correction_grid_kj_mol" in result_cg:
        cg_pair_grp.create_dataset(
            "pair_relaxation_correction_grid_kj_mol",
            data=np.asarray(result_cg["pair_relaxation_correction_grid_kj_mol"], dtype=np.float32),
        )
        cg_pair_grp.create_dataset(
            "pair_relaxation_rigid_grid_kj_mol",
            data=np.asarray(result_cg["pair_relaxation_rigid_grid_kj_mol"], dtype=np.float32),
        )
        cg_pair_grp.create_dataset(
            "pair_relaxation_relaxed_grid_kj_mol",
            data=np.asarray(result_cg["pair_relaxation_relaxed_grid_kj_mol"], dtype=np.float32),
        )
        cg_pair_grp.create_dataset(
            "pair_relaxation_face_mask",
            data=np.asarray(result_cg["pair_relaxation_face_mask"], dtype=np.int8),
        )
    if int(result_cg.get("force_match_enabled", 0)):
        cg_pair_grp.create_dataset(
            "force_match_counts",
            data=result_cg["force_match_counts"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "force_match_mean_dudr_kj_mol_nm",
            data=result_cg["force_match_mean_dudr_kj_mol_nm"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "force_match_mean_duda1_kj_mol",
            data=result_cg["force_match_mean_duda1_kj_mol"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "force_match_mean_duda2_kj_mol",
            data=result_cg["force_match_mean_duda2_kj_mol"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "force_match_angular_counts",
            data=result_cg["force_match_angular_counts"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "force_match_updated",
            data=result_cg["force_match_updated"].astype(np.int8),
        )
    if int(result_cg.get("ibi_enabled", 0)):
        cg_pair_grp.create_dataset(
            "ibi_correction_grid_kj_mol",
            data=result_cg["ibi_correction_grid_kj_mol"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "ibi_target_counts",
            data=result_cg["ibi_target_counts"].astype(np.float32),
        )
        cg_pair_grp.create_dataset(
            "ibi_model_counts",
            data=result_cg["ibi_model_counts"].astype(np.float32),
        )
    cg_pair_grp.create_dataset(
        "attractive_radial_background_kj_mol",
        data=result_cg["attractive_radial_background_kj_mol"].astype(np.float32),
    )

    if compaction_result is not None:
        comp_grp = cg_grp.create_group("cg_lipid_compaction")
        compaction_tau = float(os.environ.get("CG_LIPID_COMPACTION_THERMOSTAT_TIMESCALE", "5.0"))
        _write_common_table_contract_attrs(
            comp_grp,
            table_family="CGL-Compaction",
            source_object="isolated_dopc_axial_extension_coordinate",
            target_object="dynamic_cgl_hidden_axial_extension_coordinate",
            projection_ensemble="isolated_dopc_mc_ensemble",
            runtime_representation=(
                "clamped_bspline_self_pmf_plus_two_state_log1p_control_delta_relative_to_base_pair_table"
            ),
            correction_layer="isolated_source_compaction_state",
        )
        comp_grp.attrs["schema"] = "cg_lipid_compaction_v1"
        comp_grp.attrs["coordinate_name"] = "head_to_tail_midpoint_length_ang"
        comp_grp.attrs["coordinate_reference"] = "head_to_tail_midpoint_distance"
        comp_grp.attrs["boltzmann_temperature_upside"] = np.float32(DEFAULT_PRODUCTION_TEMP_UPSIDE)
        comp_grp.attrs["thermostat_timescale"] = np.float32(compaction_tau)
        comp_grp.attrs["mass_up"] = np.float32(
            compaction_result["self"]["effective_stiffness_eup_a2"] * compaction_tau * compaction_tau
        )
        comp_grp.attrs["effective_stiffness_eup_a2"] = np.float32(
            compaction_result["self"]["effective_stiffness_eup_a2"]
        )
        comp_grp.attrs["self_coord_min_ang"] = np.float32(
            compaction_result["self"]["coord_min_ang"]
        )
        comp_grp.attrs["self_coord_max_ang"] = np.float32(
            compaction_result["self"]["coord_max_ang"]
        )
        comp_grp.attrs["self_coord_spacing_ang"] = np.float32(
            compaction_result["self"]["coord_spacing_ang"]
        )
        comp_grp.attrs["self_n_knot"] = np.int32(compaction_result["self"]["n_knot"])
        comp_grp.attrs["compact_state_center_ang"] = np.float32(
            compaction_result["compact_center_ang"]
        )
        comp_grp.attrs["extended_state_center_ang"] = np.float32(
            compaction_result["extended_center_ang"]
        )
        comp_grp.attrs["compact_state_probability"] = np.float32(
            compaction_result["compact_probability"]
        )
        if "face_cos_min" in compaction_result:
            comp_grp.attrs["face_cos_min"] = np.float32(compaction_result["face_cos_min"])
            comp_grp.attrs["radial_cutoff_nm"] = np.float32(compaction_result["radial_cutoff_nm"])
        comp_grp.attrs["compaction_pool_size"] = np.int32(compaction_result["pool_size"])
        comp_grp.attrs["compact_pool_size"] = np.int32(compaction_result["compact_pool_size"])
        comp_grp.attrs["extended_pool_size"] = np.int32(compaction_result["extended_pool_size"])
        comp_grp.attrs["representative_compact_count"] = np.int32(
            compaction_result["representative_compact_count"]
        )
        comp_grp.attrs["representative_extended_count"] = np.int32(
            compaction_result["representative_extended_count"]
        )
        if "representative_compressed_count" in compaction_result:
            comp_grp.attrs["representative_compressed_count"] = np.int32(
                compaction_result["representative_compressed_count"]
            )
        if "compressed_pool_size" in compaction_result:
            comp_grp.attrs["compressed_pool_size"] = np.int32(
                compaction_result["compressed_pool_size"]
            )
        comp_grp.attrs["compaction_min_ang"] = np.float32(compaction_result["compaction_min_ang"])
        comp_grp.attrs["compaction_max_ang"] = np.float32(compaction_result["compaction_max_ang"])
        comp_grp.attrs["compaction_p05_ang"] = np.float32(compaction_result["compaction_p05_ang"])
        comp_grp.attrs["compaction_p95_ang"] = np.float32(compaction_result["compaction_p95_ang"])
        comp_grp.create_dataset(
            "self_coeff",
            data=np.asarray(compaction_result["self"]["self_coeff_eup"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "pmf_centers_ang",
            data=np.asarray(compaction_result["self"]["pmf_centers_ang"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "pmf_values_kj_mol",
            data=np.asarray(compaction_result["self"]["pmf_values_kj_mol"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "delta_extended_extended",
            data=np.asarray(compaction_result["delta_extended_extended"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "delta_extended_compact",
            data=np.asarray(compaction_result["delta_extended_compact"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "delta_compact_compact",
            data=np.asarray(compaction_result["delta_compact_compact"], dtype=np.float32),
        )
        for dataset_name in (
            "delta_extended_compressed",
            "delta_compact_compressed",
            "delta_compressed_compressed",
            "grid_extended_compressed_kj_mol",
            "grid_compact_compressed_kj_mol",
            "grid_compressed_compressed_kj_mol",
        ):
            if dataset_name in compaction_result:
                comp_grp.create_dataset(
                    dataset_name,
                    data=np.asarray(compaction_result[dataset_name], dtype=np.float32),
                )
        comp_grp.create_dataset(
            "grid_extended_extended_kj_mol",
            data=np.asarray(compaction_result["grid_extended_extended_kj_mol"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_extended_compact_kj_mol",
            data=np.asarray(compaction_result["grid_extended_compact_kj_mol"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_compact_compact_kj_mol",
            data=np.asarray(compaction_result["grid_compact_compact_kj_mol"], dtype=np.float32),
        )
        comp_grp.create_dataset(
            "grid_average_kj_mol",
            data=np.asarray(compaction_result["grid_average_kj_mol"], dtype=np.float32),
        )
        if "face_mask" in compaction_result:
            comp_grp.create_dataset(
                "face_mask",
                data=np.asarray(compaction_result["face_mask"], dtype=np.int8),
            )

    # CG-SC
    cg_sc_grp = cg_grp.create_group("cg_lipid_sc")
    _write_common_table_contract_attrs(
        cg_sc_grp,
        table_family="CGL-SC",
        source_object="sidechain_rotamer_bead_ensemble",
        target_object="dopc_cgl_constituent_bead_ensemble",
        projection_ensemble="sidechain_rotamers_x_cgl_conformers_x_two_orientations_x_bead_frames",
        runtime_representation="log1p_reduced_tensor_bspline",
        correction_layer="source_derived_cgl_compaction_state" if sc_compaction_delta_extended is not None else "none",
    )
    cg_sc_grp.create_dataset(
        "interaction_param",
        data=interaction_param_sc.astype(np.float32),
    )
    cg_sc_grp.create_dataset(
        "restype_order",
        data=np.asarray([np.bytes_(x) for x in sc_residue_names], dtype="S4"),
    )
    cg_sc_grp.create_dataset(
        "reference_energy_eup",
        data=sc_reference_energy[:n_sc_types, None].astype(np.float32),
    )
    cg_sc_grp.create_dataset("rms_error_kj_mol", data=sc_rms_error[:n_sc_types])
    cg_sc_grp.attrs["n_sc_types"] = n_sc_types
    cg_sc_grp.attrs["n_cg_types"] = 1
    cg_sc_grp.attrs["schema"] = (
        "cg_lipid_sc_full_tensor" if sc_n_modes == 0 else "cg_lipid_sc_quadspline"
    )
    cg_sc_grp.attrs["radial_mode"] = "full_tensor" if sc_n_modes == 0 else "full_multimode"
    cg_sc_grp.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
    cg_sc_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_sc_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_sc_grp.attrs["radial_support_source"] = "max_dopc_bead_radius_plus_dry_martini_cutoff"
    cg_sc_grp.attrs["sample_dist_min_nm"] = np.float32(NUMERICAL_DISTANCE_GUARD_NM)
    cg_sc_grp.attrs["sample_dist_min_source"] = "numerical_zero_guard_only"
    cg_sc_grp.attrs["fit_r_min_nm"] = np.float32(r_min_nm)
    cg_sc_grp.attrs["fit_r_max_nm"] = np.float32(sc_fit_r_max_nm)
    cg_sc_grp.attrs["n_modes"] = np.int32(sc_n_modes)
    cg_sc_grp.attrs["n_radial"] = np.int32(sc_n_radial)
    cg_sc_grp.attrs["n_angular"] = np.int32(sc_n_angular)
    cg_sc_grp.attrs["angular_smooth"] = np.float32(sc_angular_smooth if n_sc_types else 0.0)
    cg_sc_grp.attrs["azimuthal_count"] = np.int32(sc_azimuthal_count if n_sc_types else 0)
    cg_sc_grp.attrs["sidechain_bead_frame_count"] = np.int32(sc_bead_frame_count if n_sc_types else 0)
    cg_sc_grp.attrs["cgl_bead_frame_count"] = np.int32(cg_bead_frame_count if n_sc_types else 0)
    cg_sc_grp.attrs["conformer_count"] = np.int32(ref_ensemble_nm.shape[0])
    cg_sc_grp.attrs["orientation_sampling"] = "sidechain_and_cgl_direction_vectors"
    cg_sc_grp.attrs["knot_spacing_ang"] = np.float32(sc_knot_spacing_ang)
    cg_sc_grp.attrs["cutoff_ang"] = np.float32(sc_cutoff_ang)
    cg_sc_grp.attrs["taper_width_ang"] = np.float32(sc_taper_width_ang)
    cg_sc_grp.attrs["azimuthal_average"] = (
        first_sc_result["azimuthal_average"] if n_sc_types and first_sc_result is not None else ""
    )
    cg_sc_grp.attrs["azimuthal_average_temperature_upside"] = np.float32(
        first_sc_result["azimuthal_average_temperature_upside"]
        if n_sc_types and first_sc_result is not None
        else 0.0
    )
    if int(sc_n_modes) == 0:
        cg_sc_grp.attrs["energy_transform"] = (
            first_sc_result["energy_transform"] if n_sc_types and first_sc_result is not None else ""
        )
        cg_sc_grp.attrs["spline_control_quantity"] = (
            first_sc_result["spline_control_quantity"] if n_sc_types and first_sc_result is not None else ""
        )
        cg_sc_grp.attrs["log1p_reduced_transform"] = np.int32(1)
        cg_sc_grp.attrs["boltzmann_temperature_upside"] = np.float32(DEFAULT_PRODUCTION_TEMP_UPSIDE)
    cg_sc_grp.attrs["short_range_core_source"] = (
        first_sc_result["short_range_core_source"] if n_sc_types and first_sc_result is not None else ""
    )
    cg_sc_grp.attrs["excluded_area_source"] = (
        str(first_sc_result["excluded_area_source"])
        if n_sc_types and first_sc_result
        else ""
    )
    cg_sc_grp.attrs["isotropic_background_source"] = (
        str(first_sc_result["isotropic_background_source"])
        if n_sc_types and first_sc_result
        else ""
    )
    cg_sc_grp.attrs["attractive_control_source"] = (
        str(first_sc_result["attractive_control_source"])
        if n_sc_types and first_sc_result
        else ""
    )
    cg_sc_grp.attrs["excluded_area_nonnegative_rows"] = np.int32(
        int(first_sc_result["excluded_area_nonnegative_rows"])
        if n_sc_types and first_sc_result
        else 0
    )
    cg_sc_grp.create_dataset("short_range_core_kj_mol", data=sc_short_range_core[:n_sc_types])
    cg_sc_grp.create_dataset("short_range_core_rows", data=sc_short_range_core_rows[:n_sc_types])
    if sc_compaction_delta_extended is not None:
        cg_sc_grp.attrs["compact_state_center_ang"] = np.float32(
            compaction_states["compact_center_ang"]
        )
        cg_sc_grp.attrs["extended_state_center_ang"] = np.float32(
            compaction_states["extended_center_ang"]
        )
        cg_sc_grp.attrs["compact_state_probability"] = np.float32(
            compaction_states["compact_probability"]
        )
        cg_sc_grp.create_dataset(
            "delta_extended",
            data=sc_compaction_delta_extended.astype(np.float32),
        )
        cg_sc_grp.create_dataset(
            "delta_compact",
            data=sc_compaction_delta_compact.astype(np.float32),
        )
        cg_sc_grp.create_dataset(
            "grid_extended_kj_mol",
            data=sc_compaction_extended_grid.astype(np.float32),
        )
        cg_sc_grp.create_dataset(
            "grid_compact_kj_mol",
            data=sc_compaction_compact_grid.astype(np.float32),
        )
    _write_cg_derived_attrs(cg_sc_grp, derived_params)

    # Store the SC bead type names covered by this quadspline so that
    # convert_stage() can zero the corresponding MartiniPotential entries.
    if sc_residue_names:
        sc_bead_types_set: set = set()
        for r in sc_residue_names:
            sc_bead_types_set.update(str(bt) for bt in residue_map.get(r, []))
        sc_bead_types = sorted(sc_bead_types_set)
        cg_sc_grp.create_dataset(
            "sc_bead_types",
            data=np.asarray([np.bytes_(x) for x in sc_bead_types], dtype="S8"),
        )

    print(
        f"  Stored: CG-CG (1x1x{pair_param.size}), "
        f"CG-SC ({n_sc_types}x1x{interaction_param_sc.shape[-1]}) in {h5.filename}"
    )

    # Effective LJ parameters are retained only as target-type metadata for
    # directional CGL-target table construction. Generic MartiniPotential CGL
    # pairs are excluded during stage conversion.
    print("  Computing CGL target-type metadata...")
    effective_lj = _compute_cgl_effective_lj_params(
        ref_bead_positions_nm=ref_ensemble_nm,
        bead_types=bead_types,
        pair_params=pair_params,
        bead_frame_count=cg_bead_frame_count,
    )

    # Store effective LJ parameters so convert_stage() can read them back
    # instead of recomputing target metadata with different lipid conformations.
    eff_grp = cg_grp.create_group("effective_lj")
    eff_types = sorted(effective_lj.keys())
    eff_sigmas = np.array([effective_lj[t]["sigma_nm"] for t in eff_types], dtype=np.float32)
    eff_epsilons = np.array([effective_lj[t]["epsilon_kj_mol"] for t in eff_types], dtype=np.float32)
    eff_uncapped_sigmas = np.array(
        [effective_lj[t].get("uncapped_sigma_nm", effective_lj[t]["sigma_nm"]) for t in eff_types],
        dtype=np.float32,
    )
    eff_types_enc = np.array([np.bytes_(t) for t in eff_types], dtype="S8")
    eff_grp.create_dataset("target_types", data=eff_types_enc)
    eff_grp.create_dataset("sigma_nm", data=eff_sigmas)
    eff_grp.create_dataset("epsilon_kj_mol", data=eff_epsilons)
    eff_grp.create_dataset("uncapped_sigma_nm", data=eff_uncapped_sigmas)
    eff_grp.attrs["source"] = "orientation_average_metadata_not_runtime"
    eff_grp.attrs["cgl_bead_frame_count"] = np.int32(cg_bead_frame_count)
    eff_grp.attrs["conformer_count"] = np.int32(ref_ensemble_nm.shape[0])
    eff_grp.attrs["orientation_sampling"] = "fibonacci_cgl_direction_vectors"
    _write_cg_derived_attrs(eff_grp, derived_params)

    # Build directional B-spline tables for CGL against all non-CGL target types.
    # After this, MartiniPotential CGL-X can be omitted for all X.
    _build_cgl_target_table(
        h5, cg_grp, effective_lj,
        ref_bead_positions_nm=ref_ensemble_nm,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        derived_params=derived_params,
        temperature=DEFAULT_PRODUCTION_TEMP_UPSIDE,
        average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        compaction_states=compaction_states,
    )
    print()


def _sample_cgl_target_energy_grid(
    ref_nm: np.ndarray,
    tgt_type: str,
    bead_types: list,
    bead_charges: list,
    pair_params: dict,
    r_sample_nm: np.ndarray,
    cos_theta_grid: np.ndarray,
    orientation_dirs: np.ndarray,
    bead_frame_angles: list[float],
    average_temperature: float,
) -> np.ndarray:
    ref_nm = np.asarray(ref_nm, dtype=np.float64)
    if ref_nm.ndim == 2:
        ref_nm = ref_nm[None, :, :]
    energy_grid = np.zeros((len(r_sample_nm), cos_theta_grid.size), dtype=np.float64)
    target_charge = infer_charge_from_atomtype(tgt_type)

    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    target_pos = np.zeros((1, 3), dtype=np.float64)
    for ir, r_nm in enumerate(r_sample_nm):
        target_pos[0, :] = (float(r_nm), 0.0, 0.0)
        for ia in range(cos_theta_grid.size):
            sample_energies = []
            for direction in orientation_dirs[ia]:
                rot_base = _rotation_to_align_z_np(direction)
                for frame_angle in bead_frame_angles:
                    rot = rot_base @ _rotation_about_axis_np(z_axis, float(frame_angle))
                    for ref_conf in ref_nm:
                        cg_positions = (rot @ ref_conf.T).T
                        e, _, _ = _compute_pair_energy_and_gradient(
                            cg_positions,
                            target_pos,
                            bead_types,
                            [tgt_type],
                            bead_charges,
                            [target_charge],
                            pair_params,
                            dist_min_nm=1e-6,
                        )
                        sample_energies.append(float(e))
            energy_grid[ir, ia] = _boltzmann_free_energy_kj_mol(
                sample_energies,
                average_temperature,
            )
    return energy_grid


def _run_cgl_target_type_task(task: dict) -> tuple[int, str, np.ndarray, float, np.ndarray]:
    ti = int(task["ti"])
    tgt_type = str(task["target_type"])
    bead_types = task["bead_types"]
    bead_charges = task["bead_charges"]
    pair_params = task["pair_params"]
    r_sample_nm = task["r_sample_nm"]
    cos_theta_grid = task["cos_theta_grid"]
    orientation_dirs = task["orientation_dirs"]
    bead_frame_angles = task["bead_frame_angles"]
    temperature = float(task.get("temperature", 0.0))
    average_temperature = float(task.get("average_temperature", temperature))
    energy_grid = _sample_cgl_target_energy_grid(
        ref_nm=np.asarray(task["ref_nm"], dtype=np.float64),
        tgt_type=tgt_type,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        r_sample_nm=r_sample_nm,
        cos_theta_grid=cos_theta_grid,
        orientation_dirs=orientation_dirs,
        bead_frame_angles=bead_frame_angles,
        average_temperature=average_temperature,
    )

    reference_energy_kj_mol = float(np.min(energy_grid))
    kbt_kj_mol = float(temperature) * ENERGY_CONVERSION_KJ_PER_EUP
    if kbt_kj_mol <= 0.0:
        raise RuntimeError("cg_lipid_target Boltzmann-weight PMF table requires positive temperature")
    reduced_energy = (energy_grid - reference_energy_kj_mol) / kbt_kj_mol
    control = np.log1p(np.maximum(reduced_energy, 0.0))
    return ti, tgt_type, control.reshape(-1), reference_energy_kj_mol, energy_grid


def _build_cgl_target_table(
    h5: h5py.File,
    cg_grp: h5py.Group,
    effective_lj: dict,
    ref_bead_positions_nm: np.ndarray | None = None,
    bead_types: list | None = None,
    bead_charges: list | None = None,
    pair_params: dict | None = None,
    derived_params: dict | None = None,
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
    n_knot_radial: int = 14,
    n_knot_angular: int = 321,
    knot_spacing_ang: float = 0.35,
    temperature: float = 0.0,
    average_temperature: float | None = None,
    compaction_states: dict | None = None,
) -> None:
    """Build directional tensor B-spline tables for CGL-point targets."""
    target_types = sorted(t for t in effective_lj if t != "CGL")
    if not target_types:
        print("  cg_lipid_target: no target types, skipping")
        return
    if average_temperature is None:
        average_temperature = temperature
    average_temperature = float(average_temperature)

    # Sample densely for an accurate B-spline fit.  The interaction samples
    # explicit DOPC bead-vs-target energies and only uses a near-zero numerical
    # guard to avoid division by zero for exact bead/target overlaps.
    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64) if ref_bead_positions_nm is not None else None
    explicit_source = (
        ref_nm is not None
        and (ref_nm.shape == (14, 3) or (ref_nm.ndim == 3 and ref_nm.shape[1:] == (14, 3)))
        and bead_types is not None
        and bead_charges is not None
        and pair_params is not None
    )
    if not explicit_source:
        raise RuntimeError("cg_lipid_target requires explicit DOPC bead geometry and dry-MARTINI parameters")
    ref_nm = _canonicalize_lipid_reference_ensemble_to_z(ref_nm)
    max_radius_ang = float(np.max(np.linalg.norm(ref_nm.reshape(-1, 3), axis=1)) * length_conv)
    extended_radial_count = (
        int(math.ceil((max_radius_ang + DRY_MARTINI_NONBONDED_CUTOFF_NM * length_conv) / knot_spacing_ang))
        + 2
    )
    n_knot_radial = max(int(n_knot_radial), extended_radial_count)
    n_types = len(target_types)
    n_param = n_knot_radial * n_knot_angular
    interaction_param = np.zeros((1, n_types, n_param), dtype=np.float64)
    reference_energy_eup = np.zeros((1, n_types), dtype=np.float64)
    r_min_ang = float(knot_spacing_ang) + 0.1
    r_max_ang = float((n_knot_radial - 2) * knot_spacing_ang)
    r_sample_ang = np.linspace(float(knot_spacing_ang), r_max_ang, n_knot_radial)
    r_sample_nm = r_sample_ang / length_conv
    t_angular_sample = np.linspace(1.0, float(n_knot_angular - 2), n_knot_angular)
    cos_theta_grid = (t_angular_sample - 1.0) / (0.5 * float(n_knot_angular - 3)) - 1.0
    cos_theta_grid = np.asarray(cos_theta_grid, dtype=np.float64)
    base_energy_grid = np.zeros((n_types, len(r_sample_nm), len(cos_theta_grid)), dtype=np.float64)
    target_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    phi_values = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)
    bead_frame_angles = _bead_frame_angles(_bead_frame_count("CGL", 1))
    orientation_dirs = _directions_with_dot_np(target_axis, cos_theta_grid, phi_values)
    target_fit_smooth = 0.01

    target_tasks = [
        {
            "ti": ti,
            "target_type": tgt_type,
            "bead_types": list(bead_types),
            "bead_charges": list(bead_charges),
            "pair_params": pair_params,
            "ref_nm": ref_nm,
            "r_sample_nm": r_sample_nm,
            "cos_theta_grid": cos_theta_grid,
            "orientation_dirs": orientation_dirs,
            "bead_frame_angles": bead_frame_angles,
            "n_knot_radial": n_knot_radial,
            "n_knot_angular": n_knot_angular,
            "knot_spacing_ang": knot_spacing_ang,
            "energy_conv": energy_conv,
            "r_min_ang": r_min_ang,
            "temperature": float(temperature),
            "average_temperature": float(average_temperature),
        }
        for ti, tgt_type in enumerate(target_types)
    ]
    for ti, _tgt_type, control_flat, reference_energy_kj_mol, energy_grid in _parallel_map_ordered(
        "CGL-particle target table", _run_cgl_target_type_task, target_tasks
    ):
        control_grid = np.asarray(control_flat, dtype=np.float64).reshape(
            len(r_sample_nm),
            len(cos_theta_grid),
        )
        interaction_param[0, ti, :] = _fit_radial_angular_tensor_bspline(
            r_sample_nm,
            cos_theta_grid,
            control_grid,
            n_knot_radial=n_knot_radial,
            n_knot_angular=n_knot_angular,
            knot_spacing_ang=knot_spacing_ang,
            energy_conversion=1.0,
            smooth=target_fit_smooth,
        ).reshape(-1)
        reference_energy_eup[0, ti] = float(reference_energy_kj_mol) / float(energy_conv)
        base_energy_grid[ti] = np.asarray(energy_grid, dtype=np.float64)

    target_compaction_delta_extended = None
    target_compaction_delta_compact = None
    target_compaction_extended_grid = None
    target_compaction_compact_grid = None
    if compaction_states is not None:
        extended_tasks = [
            {
                **task,
                "ref_nm": np.asarray(compaction_states["extended_refs_nm"], dtype=np.float64),
            }
            for task in target_tasks
        ]
        compact_tasks = [
            {
                **task,
                "ref_nm": np.asarray(compaction_states["compact_refs_nm"], dtype=np.float64),
            }
            for task in target_tasks
        ]
        target_compaction_delta_extended = np.zeros((1, n_types, n_param), dtype=np.float32)
        target_compaction_delta_compact = np.zeros((1, n_types, n_param), dtype=np.float32)
        target_compaction_extended_grid = np.zeros_like(base_energy_grid, dtype=np.float32)
        target_compaction_compact_grid = np.zeros_like(base_energy_grid, dtype=np.float32)

        extended_results = {
            ti: np.asarray(energy_grid, dtype=np.float64)
            for ti, _tgt_type, _control_flat, _reference_energy_kj_mol, energy_grid in _parallel_map_ordered(
                "CGL-particle target extended-state table",
                _run_cgl_target_type_task,
                extended_tasks,
            )
        }
        compact_results = {
            ti: np.asarray(energy_grid, dtype=np.float64)
            for ti, _tgt_type, _control_flat, _reference_energy_kj_mol, energy_grid in _parallel_map_ordered(
                "CGL-particle target compacted-state table",
                _run_cgl_target_type_task,
                compact_tasks,
            )
        }
        for ti in range(n_types):
            delta = _build_single_cgl_state_delta_radial_angular(
                base_energy_grid_kj_mol=base_energy_grid[ti],
                extended_energy_grid_kj_mol=extended_results[ti],
                compact_energy_grid_kj_mol=compact_results[ti],
                compressed_energy_grid_kj_mol=None,
                reference_energy_eup=float(reference_energy_eup[0, ti]),
                temperature_upside=float(temperature),
                r_values_nm=r_sample_nm,
                cos_theta_grid=cos_theta_grid,
                n_knot_radial=n_knot_radial,
                n_knot_angular=n_knot_angular,
                knot_spacing_ang=knot_spacing_ang,
                smooth=target_fit_smooth,
            )
            target_compaction_delta_extended[0, ti, :] = np.asarray(
                delta["delta_extended"], dtype=np.float32
            )
            target_compaction_delta_compact[0, ti, :] = np.asarray(
                delta["delta_compact"], dtype=np.float32
            )
            target_compaction_extended_grid[ti] = np.asarray(
                delta["grid_extended_kj_mol"], dtype=np.float32
            )
            target_compaction_compact_grid[ti] = np.asarray(
                delta["grid_compact_kj_mol"], dtype=np.float32
            )

    target_grp = cg_grp.create_group("cg_lipid_target")
    _write_common_table_contract_attrs(
        target_grp,
        table_family="CGL-particle",
        source_object="dopc_cgl_constituent_bead_ensemble",
        target_object="dry_martini_particle_type",
        projection_ensemble="cgl_conformers_x_cgl_orientations_x_bead_frames",
        runtime_representation="log1p_reduced_radial_angular_bspline",
        correction_layer="source_derived_cgl_compaction_state" if compaction_states is not None else "none",
    )
    target_grp.create_dataset(
        "interaction_param",
        data=interaction_param.astype(np.float32),
    )
    target_grp.create_dataset(
        "reference_energy_eup",
        data=reference_energy_eup.astype(np.float32),
    )
    target_grp.create_dataset(
        "target_order",
        data=np.asarray([np.bytes_(x) for x in target_types], dtype="S8"),
    )
    target_grp.attrs["n_target_types"] = np.int32(n_types)
    target_grp.attrs["n_cg_types"] = np.int32(1)
    target_grp.attrs["schema"] = "cg_lipid_target"
    target_grp.attrs["n_modes"] = np.int32(0)
    target_grp.attrs["n_radial"] = np.int32(n_knot_radial)
    target_grp.attrs["n_angular"] = np.int32(n_knot_angular)
    target_grp.attrs["azimuthal_count"] = np.int32(len(phi_values))
    target_grp.attrs["radial_sample_count"] = np.int32(len(r_sample_ang))
    target_grp.attrs["angular_sample_count"] = np.int32(len(cos_theta_grid))
    target_grp.attrs["cgl_bead_frame_count"] = np.int32(len(bead_frame_angles))
    target_grp.attrs["conformer_count"] = np.int32(ref_nm.shape[0])
    target_grp.attrs["orientation_sampling"] = "cgl_direction_vector"
    target_grp.attrs["knot_spacing_ang"] = np.float32(knot_spacing_ang)
    cutoff_ang = float((n_knot_radial - 2) * knot_spacing_ang)
    target_grp.attrs["cutoff_ang"] = np.float32(cutoff_ang)
    target_grp.attrs["taper_width_ang"] = np.float32(knot_spacing_ang)
    target_grp.attrs["azimuthal_average"] = (
        "energy_expectation"
        if average_temperature <= 0.0
        else (
            "tempered_boltzmann_free_energy"
            if abs(average_temperature - float(temperature)) > 1e-8
            else "boltzmann_free_energy"
        )
    )
    target_grp.attrs["azimuthal_average_temperature_upside"] = np.float32(average_temperature)
    target_grp.attrs["energy_transform"] = (
        "log1p_reduced_energy_expectation"
        if average_temperature <= 0.0
        else (
            "log1p_reduced_tempered_pmf"
            if abs(average_temperature - float(temperature)) > 1e-8
            else "log1p_reduced_pmf"
        )
    )
    target_grp.attrs["log1p_reduced_transform"] = np.int32(1)
    target_grp.attrs["boltzmann_weight_transform"] = np.int32(0)
    target_grp.attrs["boltzmann_temperature_upside"] = np.float32(temperature)
    target_grp.attrs["spline_control_quantity"] = (
        "log1p_reduced_energy_expectation"
        if average_temperature <= 0.0
        else (
            "log1p_reduced_tempered_free_energy"
            if abs(average_temperature - float(temperature)) > 1e-8
            else "log1p_reduced_free_energy"
        )
    )
    target_grp.attrs["unresolved_core_source"] = "angular_resolved_first_sampled_dry_martini_energy"
    target_grp.attrs["excluded_area_source"] = "none_full_resolved_dry_martini_cgl_target_table"
    target_grp.attrs["excluded_area_nonnegative_rows"] = np.int32(0)
    target_grp.attrs["source"] = "explicit_dopc_directional"
    target_grp.attrs["sample_dist_min_nm"] = np.float32(1e-6)
    target_grp.attrs["sample_dist_min_source"] = "numerical_zero_guard_only"
    target_grp.attrs["isotropic_background_source"] = "none_full_resolved_dry_martini_cgl_target_table"
    target_grp.attrs["attractive_control_source"] = "retained_full_resolved_dry_martini_cgl_target_table"
    target_grp.attrs["relaxation"] = "ensemble_rotated_geometry"
    target_grp.attrs["angle_convention"] = "ang=n_cgl_dot_n12"
    target_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    target_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    target_grp.attrs["radial_support_source"] = "max_dopc_bead_radius_plus_dry_martini_cutoff"
    if compaction_states is not None:
        target_grp.attrs["compact_state_center_ang"] = np.float32(
            compaction_states["compact_center_ang"]
        )
        target_grp.attrs["extended_state_center_ang"] = np.float32(
            compaction_states["extended_center_ang"]
        )
        target_grp.attrs["compact_state_probability"] = np.float32(
            compaction_states["compact_probability"]
        )
        target_grp.create_dataset(
            "delta_extended",
            data=target_compaction_delta_extended.astype(np.float32),
        )
        target_grp.create_dataset(
            "delta_compact",
            data=target_compaction_delta_compact.astype(np.float32),
        )
        target_grp.create_dataset(
            "grid_extended_kj_mol",
            data=target_compaction_extended_grid.astype(np.float32),
        )
        target_grp.create_dataset(
            "grid_compact_kj_mol",
            data=target_compaction_compact_grid.astype(np.float32),
        )
    if derived_params:
        _write_cg_derived_attrs(target_grp, derived_params)

    print(f"  cg_lipid_target: {n_types} target types, "
          f"{n_knot_radial} radial x {n_knot_angular} angular knots, "
          f"cutoff={cutoff_ang:.1f} A, source={target_grp.attrs['source']}")


def build_martini_tables(
    output_path: Path,
    dry_ff_path: Path,
    martinize_path: Path,
    sidechain_lib_path: Path,
    forcefield_name: str = "martini22",
    active_residue_names: Iterable[str] | None = None,
    active_atom_types: Set[str] | None = None,
    r_count: int = 24,
    direction_count: int = 16,
    cos_theta_count: int = 13,
    cg_lipid_config: dict | None = None,
) -> None:
    output_path = Path(output_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    martinize_path = Path(martinize_path).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
        martinize_path, forcefield_name
    )
    active_residue_names = list(active_residue_names or CANONICAL_RESIDUES)
    active_atom_types = set(active_atom_types or atomtypes)

    cg_lipid_config = dict(cg_lipid_config or {})
    bead_types = cg_lipid_config.get("bead_types")
    bead_charges = cg_lipid_config.get("bead_charges")
    lipids_itp_path = cg_lipid_config.get("lipids_itp_path")
    if lipids_itp_path is None:
        candidate = dry_ff_path.parent / "dry_martini_v2.1_lipids.itp"
        lipids_itp_path = candidate if candidate.exists() else None

    if bead_types is not None and bead_charges is None:
        bead_charges = [infer_charge_from_atomtype(bt) for bt in bead_types]

    with h5py.File(output_path, "w") as h5:
        h5.attrs["schema"] = "martini_combined"
        _build_particles_group(h5, atomtypes, pair_params, active_atom_types)
        _build_sc_table_group(
            h5,
            residue_map,
            pair_params,
            sidechain_lib_path,
            active_residue_names=active_residue_names,
            active_target_types=sorted(active_atom_types),
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )
        if cg_lipid_config:
            _build_cg_lipid_tables(
                h5,
                pair_params=pair_params,
                sidechain_lib_path=sidechain_lib_path,
                martinize_path=martinize_path,
                forcefield_name=forcefield_name,
                active_residue_names=active_residue_names,
                ref_bead_positions_nm=cg_lipid_config.get("ref_bead_positions_nm"),
                bead_types=bead_types,
                bead_charges=bead_charges,
                bead_masses_g_mol=parse_itp_atomtype_masses(dry_ff_path),
                lipids_itp_path=Path(lipids_itp_path) if lipids_itp_path else None,
                r_count=r_count,
                cos_theta_count=cos_theta_count,
                azimuthal_count=direction_count,
            )
    print(f"Built {output_path}")


def build_particle_h5(
    output_path: Path,
    dry_ff_path: Path,
) -> None:
    """Generate particle.h5 with all particle types from the ITP."""
    output_path = Path(output_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)

    def _writer(h5: h5py.File) -> None:
        h5.attrs["schema"] = SCHEMA_PARTICLES
        _build_particles_group(h5, atomtypes, pair_params, set(atomtypes))

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} ({len(atomtypes)} types)")


def build_sidechain_h5(
    output_path: Path,
    dry_ff_path: Path,
    martinize_path: Path,
    sidechain_lib_path: Path,
    forcefield_name: str = "martini22",
) -> None:
    """Generate sidechain.h5 with all residues and target types."""
    output_path = Path(output_path).expanduser().resolve()
    martinize_path = Path(martinize_path).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path).expanduser().resolve()
    atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
        martinize_path, forcefield_name
    )

    def _writer(h5: h5py.File) -> None:
        h5.attrs["schema"] = SCHEMA_SC
        _build_sc_table_group(
            h5, residue_map, pair_params, sidechain_lib_path,
            active_residue_names=list(CANONICAL_RESIDUES),
            active_target_types=atomtypes,
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path}")


def build_dopc_h5(
    output_path: Path,
    dry_ff_path: Path,
    lipids_itp_path: Path,
    martinize_path: Path,
    sidechain_lib_path: Path,
    dopc_pdb_path: Path,
    forcefield_name: str = "martini22",
    conformer_upside_h5_path: Path | None = None,
    conformer_upside_h5_paths: list[Path] | None = None,
    conformer_count: int = 4,
    conformer_frame_start_fraction: float = 0.5,
    conformer_min_nonlipid_xy_distance_ang: float = 0.0,
    conformer_max_nonlipid_xy_distance_ang: float = 0.0,
    conformer_max_nonlipid_distance_ang: float = 0.0,
    isolated_conformer_count: int = 8,
    isolated_conformer_pool_count: int = 32,
    isolated_conformer_temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    isolated_conformer_seed: int = 1729,
    isolated_conformer_burnin_steps: int = 5000,
    isolated_conformer_mc_steps: int = 2000,
    isolated_conformer_proposal_sigma_nm: float = 0.025,
    fluid_pair_pmf_upside_h5_path: Path | None = None,
    fluid_pair_pmf_frame_start_fraction: float = 0.5,
    force_match_upside_h5_path: Path | None = None,
    force_match_frame_start_fraction: float = 0.5,
    force_match_max_frames: int = 100,
    force_match_min_count: int = 4,
    force_match_mode: str = "radial",
    ibi_base_pair_h5_path: Path | None = None,
    ibi_target_upside_h5_path: Path | None = None,
    ibi_model_upside_h5_path: Path | None = None,
    ibi_target_frame_start_fraction: float = 0.5,
    ibi_model_frame_start_fraction: float = 0.5,
) -> None:
    """Generate dopc.h5 with DOPC CG lipid tables."""
    output_path = Path(output_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    lipids_itp_path = Path(lipids_itp_path).expanduser().resolve()
    martinize_path = Path(martinize_path).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path).expanduser().resolve()
    dopc_pdb_path = Path(dopc_pdb_path).expanduser().resolve()
    from martini_itp_reader import parse_dopc_from_itp

    atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    atomtype_masses = parse_itp_atomtype_masses(dry_ff_path)
    dopc = parse_dopc_from_itp(lipids_itp_path)
    reference_upside_paths = []
    if conformer_upside_h5_path is not None:
        reference_upside_paths.append(Path(conformer_upside_h5_path).expanduser().resolve())
    if conformer_upside_h5_paths is not None:
        reference_upside_paths.extend(
            Path(path).expanduser().resolve() for path in conformer_upside_h5_paths
        )

    ref_bead_positions_nm, reference_source, representative_metadata = _resolve_dopc_reference_ensemble(
        dopc=dopc,
        pair_params=pair_params,
        lipids_itp_path=lipids_itp_path,
        conformer_upside_h5_path=conformer_upside_h5_path,
        conformer_upside_h5_paths=conformer_upside_h5_paths,
        conformer_count=conformer_count,
        conformer_frame_start_fraction=conformer_frame_start_fraction,
        conformer_min_nonlipid_xy_distance_ang=conformer_min_nonlipid_xy_distance_ang,
        conformer_max_nonlipid_xy_distance_ang=conformer_max_nonlipid_xy_distance_ang,
        conformer_max_nonlipid_distance_ang=conformer_max_nonlipid_distance_ang,
        isolated_conformer_count=isolated_conformer_count,
        isolated_conformer_pool_count=isolated_conformer_pool_count,
        isolated_conformer_temperature_upside=isolated_conformer_temperature_upside,
        isolated_conformer_seed=isolated_conformer_seed,
        isolated_conformer_burnin_steps=isolated_conformer_burnin_steps,
        isolated_conformer_mc_steps=isolated_conformer_mc_steps,
        isolated_conformer_proposal_sigma_nm=isolated_conformer_proposal_sigma_nm,
    )

    def _writer(h5: h5py.File) -> None:
        h5.attrs["dopc_reference_source"] = reference_source
        h5.attrs["dopc_reference_conformer_count"] = np.int32(
            int(np.asarray(ref_bead_positions_nm).reshape(-1, 14, 3).shape[0])
        )
        if reference_upside_paths:
            h5.attrs["dopc_reference_up_file_count"] = np.int32(len(reference_upside_paths))
            h5.attrs["dopc_reference_up_files"] = np.asarray(
                [np.bytes_(str(path)) for path in reference_upside_paths],
                dtype="S512",
            )
            if conformer_min_nonlipid_xy_distance_ang > 0.0:
                h5.attrs["dopc_reference_min_nonlipid_xy_distance_ang"] = np.float32(
                    conformer_min_nonlipid_xy_distance_ang
                )
            if conformer_max_nonlipid_xy_distance_ang > 0.0:
                h5.attrs["dopc_reference_max_nonlipid_xy_distance_ang"] = np.float32(
                    conformer_max_nonlipid_xy_distance_ang
                )
            if conformer_max_nonlipid_distance_ang > 0.0:
                h5.attrs["dopc_reference_max_nonlipid_distance_ang"] = np.float32(
                    conformer_max_nonlipid_distance_ang
                )
        else:
            h5.attrs["dopc_isolated_conformer_temperature_upside"] = np.float32(
                isolated_conformer_temperature_upside
            )
            h5.attrs["dopc_isolated_conformer_seed"] = np.int32(isolated_conformer_seed)
            h5.attrs["dopc_isolated_conformer_burnin_steps"] = np.int32(
                isolated_conformer_burnin_steps
            )
            h5.attrs["dopc_isolated_conformer_mc_steps"] = np.int32(isolated_conformer_mc_steps)
            h5.attrs["dopc_isolated_conformer_proposal_sigma_nm"] = np.float32(
                isolated_conformer_proposal_sigma_nm
            )
            h5.attrs["dopc_isolated_conformer_pool_count"] = np.int32(
                int(max(int(isolated_conformer_pool_count), int(isolated_conformer_count)))
            )
            if representative_metadata is not None:
                h5.attrs["dopc_isolated_conformer_representative_count"] = np.int32(
                    int(np.asarray(representative_metadata["representative_refs_nm"]).shape[0])
                )
                h5.attrs["dopc_isolated_conformer_representative_compaction_min_ang"] = np.float32(
                    float(np.min(representative_metadata["representative_compaction_ang"]))
                )
                h5.attrs["dopc_isolated_conformer_representative_compaction_max_ang"] = np.float32(
                    float(np.max(representative_metadata["representative_compaction_ang"]))
                )
        _build_cg_lipid_tables(
            h5,
            pair_params=pair_params,
            sidechain_lib_path=sidechain_lib_path,
            martinize_path=martinize_path,
            forcefield_name=forcefield_name,
            active_residue_names=list(CANONICAL_RESIDUES),
            ref_bead_positions_nm=ref_bead_positions_nm,
            bead_types=dopc["bead_types"],
            bead_charges=dopc["bead_charges"],
            bead_masses_g_mol=atomtype_masses,
            lipids_itp_path=lipids_itp_path,
            fluid_pair_pmf_upside_h5_path=fluid_pair_pmf_upside_h5_path,
            fluid_pair_pmf_frame_start_fraction=fluid_pair_pmf_frame_start_fraction,
            force_match_upside_h5_path=force_match_upside_h5_path,
            force_match_frame_start_fraction=force_match_frame_start_fraction,
            force_match_max_frames=force_match_max_frames,
            force_match_min_count=force_match_min_count,
            force_match_mode=force_match_mode,
            ibi_base_pair_h5_path=ibi_base_pair_h5_path,
            ibi_target_upside_h5_path=ibi_target_upside_h5_path,
            ibi_model_upside_h5_path=ibi_model_upside_h5_path,
            ibi_target_frame_start_fraction=ibi_target_frame_start_fraction,
            ibi_model_frame_start_fraction=ibi_model_frame_start_fraction,
        )

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path}")


def _resolve_dopc_reference_ensemble(
    dopc: dict,
    pair_params: dict,
    lipids_itp_path: Path,
    conformer_upside_h5_path: Path | None = None,
    conformer_upside_h5_paths: list[Path] | None = None,
    conformer_count: int = 4,
    conformer_frame_start_fraction: float = 0.5,
    conformer_min_nonlipid_xy_distance_ang: float = 0.0,
    conformer_max_nonlipid_xy_distance_ang: float = 0.0,
    conformer_max_nonlipid_distance_ang: float = 0.0,
    isolated_conformer_count: int = 8,
    isolated_conformer_pool_count: int = 32,
    isolated_conformer_temperature_upside: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
    isolated_conformer_seed: int = 1729,
    isolated_conformer_burnin_steps: int = 5000,
    isolated_conformer_mc_steps: int = 2000,
    isolated_conformer_proposal_sigma_nm: float = 0.025,
    append_isolated_conformers: bool = False,
) -> tuple[np.ndarray, str, dict | None]:
    reference_upside_paths = []
    if conformer_upside_h5_path is not None:
        reference_upside_paths.append(Path(conformer_upside_h5_path).expanduser().resolve())
    if conformer_upside_h5_paths is not None:
        reference_upside_paths.extend(
            Path(path).expanduser().resolve() for path in conformer_upside_h5_paths
        )
    if reference_upside_paths:
        ref_bead_positions_nm = load_dopc_conformers_from_upside_h5_pool(
            reference_upside_paths,
            max_conformers=conformer_count,
            frame_start_fraction=conformer_frame_start_fraction,
            min_nonlipid_xy_distance_ang=conformer_min_nonlipid_xy_distance_ang,
            max_nonlipid_xy_distance_ang=conformer_max_nonlipid_xy_distance_ang,
            max_nonlipid_distance_ang=conformer_max_nonlipid_distance_ang,
        )
        uses_xy_shell = conformer_max_nonlipid_xy_distance_ang > 0.0
        uses_3d_shell = conformer_max_nonlipid_distance_ang > 0.0
        uses_bulk_filter = conformer_min_nonlipid_xy_distance_ang > 0.0
        if (
            len(reference_upside_paths) == 1
            and not uses_bulk_filter
            and not uses_xy_shell
            and not uses_3d_shell
        ):
            reference_source = "full_resolution_upside_conformer_trajectory"
        elif (
            len(reference_upside_paths) == 1
            and uses_bulk_filter
            and not uses_xy_shell
            and not uses_3d_shell
        ):
            reference_source = "full_resolution_upside_bulk_lipid_filtered_conformer_trajectory"
        elif (
            len(reference_upside_paths) == 1
            and not uses_bulk_filter
            and uses_xy_shell
            and not uses_3d_shell
        ):
            reference_source = "full_resolution_upside_interface_lipid_conformer_trajectory"
        elif (
            len(reference_upside_paths) == 1
            and not uses_bulk_filter
            and not uses_xy_shell
            and uses_3d_shell
        ):
            reference_source = "full_resolution_upside_local_3d_lipid_conformer_trajectory"
        elif not uses_bulk_filter and not uses_xy_shell and not uses_3d_shell:
            reference_source = "pooled_full_resolution_upside_conformer_trajectory"
        elif uses_bulk_filter and not uses_xy_shell and not uses_3d_shell:
            reference_source = "pooled_full_resolution_upside_bulk_lipid_filtered_conformer_trajectory"
        elif not uses_bulk_filter and uses_xy_shell and not uses_3d_shell:
            reference_source = "pooled_full_resolution_upside_interface_lipid_conformer_trajectory"
        elif not uses_bulk_filter and not uses_xy_shell and uses_3d_shell:
            reference_source = "pooled_full_resolution_upside_local_3d_lipid_conformer_trajectory"
        else:
            reference_source = "pooled_full_resolution_upside_shell_filtered_conformer_trajectory"
        traj_metadata = {
            "source_kind": "upside_trajectory_pool",
            "reference_paths": [str(path) for path in reference_upside_paths],
            "min_nonlipid_xy_distance_ang": float(conformer_min_nonlipid_xy_distance_ang),
            "max_nonlipid_xy_distance_ang": float(conformer_max_nonlipid_xy_distance_ang),
            "max_nonlipid_distance_ang": float(conformer_max_nonlipid_distance_ang),
        }
        if not append_isolated_conformers:
            return (
                np.asarray(ref_bead_positions_nm, dtype=np.float64),
                reference_source,
                traj_metadata,
            )

        pool_count = max(int(isolated_conformer_pool_count), int(isolated_conformer_count))
        pool_refs_nm = sample_isolated_dopc_bonded_conformers(
            dopc,
            lipids_itp_path=lipids_itp_path,
            pair_params=pair_params,
            conformer_count=pool_count,
            temperature_upside=isolated_conformer_temperature_upside,
            seed=isolated_conformer_seed,
            mc_burnin_steps=isolated_conformer_burnin_steps,
            mc_steps_per_conformer=isolated_conformer_mc_steps,
            proposal_sigma_nm=isolated_conformer_proposal_sigma_nm,
        )
        representative_metadata = _select_reference_ensemble_representatives(
            pool_refs_nm,
            representative_count=isolated_conformer_count,
        )
        mixed_refs_nm = np.concatenate(
            [
                np.asarray(ref_bead_positions_nm, dtype=np.float64),
                np.asarray(representative_metadata["representative_refs_nm"], dtype=np.float64),
            ],
            axis=0,
        )
        return (
            np.asarray(mixed_refs_nm, dtype=np.float64),
            f"{reference_source}_plus_isolated_dopc_representative_ensemble",
            {
                **traj_metadata,
                "source_kind": "upside_trajectory_plus_isolated_representative",
                "traj_reference_source": reference_source,
                "traj_conformer_count": int(np.asarray(ref_bead_positions_nm).shape[0]),
                "isolated_conformer_count": int(
                    np.asarray(representative_metadata["representative_refs_nm"]).shape[0]
                ),
                "isolated_reference_pool_count": int(pool_count),
                "isolated_reference_seed": int(isolated_conformer_seed),
            },
        )

    pool_count = max(int(isolated_conformer_pool_count), int(isolated_conformer_count))
    pool_refs_nm = sample_isolated_dopc_bonded_conformers(
        dopc,
        lipids_itp_path=lipids_itp_path,
        pair_params=pair_params,
        conformer_count=pool_count,
        temperature_upside=isolated_conformer_temperature_upside,
        seed=isolated_conformer_seed,
        mc_burnin_steps=isolated_conformer_burnin_steps,
        mc_steps_per_conformer=isolated_conformer_mc_steps,
        proposal_sigma_nm=isolated_conformer_proposal_sigma_nm,
    )
    representative_metadata = _select_reference_ensemble_representatives(
        pool_refs_nm,
        representative_count=isolated_conformer_count,
    )
    return (
        np.asarray(representative_metadata["representative_refs_nm"], dtype=np.float64),
        "isolated_dopc_itp_bonded_plus_intramolecular_nonbonded_mc_representative_ensemble",
        {
            **representative_metadata,
            "source_kind": "isolated_representative",
        },
    )
