#!/usr/bin/env python3
"""Build dry-MARTINI spline tables from ITP-derived parameters."""

from __future__ import annotations

import math
import os
import importlib.util
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
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
PARTICLES_GRID_N = 1000
PARTICLES_R_MIN_A = 0.0
PARTICLES_R_MAX_A = 12.0
DRY_MARTINI_NONBONDED_CUTOFF_NM = PARTICLES_R_MAX_A * ANGSTROM_TO_NM

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)

CANONICAL_CB_POSITION_ANG = (0.0, 0.94375626, 1.2068012)
_cb_norm = math.sqrt(sum(x * x for x in CANONICAL_CB_POSITION_ANG))
CANONICAL_CB_VECTOR_UNIT = tuple(x / _cb_norm for x in CANONICAL_CB_POSITION_ANG)

SCHEMA_PARTICLES = "martini_particles_combined_v1"
SCHEMA_SC = "martini_sc_combined_v1"

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
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            return list(pool.map(func, tasks))
    except (OSError, PermissionError) as exc:
        print(f"  Parallel {label}: process workers unavailable ({exc}); using threads")
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(func, tasks))


def _bead_frame_angles(count: int) -> np.ndarray:
    count = max(1, int(count))
    return np.linspace(0.0, 2.0 * np.pi, count, endpoint=False, dtype=np.float64)


def _bead_frame_count(kind: str, default: int = 1) -> int:
    default = int(default)
    specific = os.environ.get(f"UPSIDE_MARTINI_{kind.upper()}_BEAD_FRAME_COUNT", "").strip()
    if specific:
        return _positive_int_env(f"UPSIDE_MARTINI_{kind.upper()}_BEAD_FRAME_COUNT", default)
    shared = os.environ.get("UPSIDE_MARTINI_BEAD_FRAME_COUNT", "").strip()
    if shared:
        return _positive_int_env("UPSIDE_MARTINI_BEAD_FRAME_COUNT", default)
    return default


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
    return _compute_lipid_bonded_energy(positions)


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
    rel_relax_steps: int = 0,
    sc_restraint_k: float = 5000.0,
    dist_min_nm: float = 0.10,
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
        energy_sum = [0.0 for _ in cos_theta_grid]
        energy_weight = [0.0 for _ in cos_theta_grid]
        rotamer_energy_sum = [[0.0 for _ in cos_theta_grid] for _ in range(n_rotamer)]
        rotamer_energy_weight = [[0.0 for _ in cos_theta_grid] for _ in range(n_rotamer)]

        for direction in direction_vectors:
            target_pos_nm = [
                cb_anchor_nm[0] + r_nm * direction[0],
                cb_anchor_nm[1] + r_nm * direction[1],
                cb_anchor_nm[2] + r_nm * direction[2],
            ]
            cos_theta = _clamp(-_dot3(direction, cb_vector_unit), -1.0, 1.0)
            total_energy = 0.0

            for irot, (sc_positions, rot_weight) in enumerate(
                zip(rotamer_positions, rotamer_weights)
            ):
                rot_energy_sum = 0.0
                target = np.asarray(target_pos_nm, dtype=np.float64).reshape(1, 3)
                for sc_frame_angle in bead_frame_angles:
                    if rel_relax_steps > 0:
                        body_rot = (_rotation_about_axis_np(np.array([0.0, 0.0, 1.0]), float(sc_frame_angle))
                                     @ sc_positions.T).T
                        sc_body, rot_energy = _relax_sc_beads(
                            init_body_positions=body_rot,
                            ref_body_positions=np.asarray(sc_positions, dtype=np.float64),
                            sc_bead_types=list(sidechain_bead_types),
                            sc_bead_charges=list(sidechain_bead_charges),
                            k_restraint=sc_restraint_k,
                            pair_params=pair_params,
                            direction_lab=cb_vector_arr,
                            cb_anchor_lab=cb_anchor_arr,
                            partner_lab_positions=target,
                            partner_bead_types=[target_type],
                            partner_bead_charges=[target_charge],
                            dist_min_nm=float(dist_min_nm),
                            rel_relax_steps=rel_relax_steps,
                        )
                        # Add WCA excluded area per-bead-pair using relaxed positions.
                        sc_lab = _dopc_body_to_lab(sc_body, cb_vector_arr, cb_anchor_arr)
                        for bead_pos, bead_type in zip(sc_lab, sidechain_bead_types):
                            delta = target[0] - bead_pos
                            d_nm = max(1e-6, float(np.linalg.norm(delta)))
                            params = pair_params[(bead_type, target_type)]
                            sigma_nm = float(params["sigma_nm"])
                            contact_nm = (2.0 ** (1.0 / 6.0)) * sigma_nm
                            if d_nm < contact_nm:
                                sr = sigma_nm / max(d_nm, 1e-6)
                                sr6 = sr ** 6
                                rot_energy += 4.0 * DEFAULT_PRODUCTION_KBT_KJ_MOL * (sr6 * sr6 - sr6) + DEFAULT_PRODUCTION_KBT_KJ_MOL
                    else:
                        framed_sc_positions = _rotate_points_about_axis_np(
                            sc_positions, cb_vector_arr, sc_frame_angle, cb_anchor_arr
                        )
                        rot_energy = 0.0
                        for bead_pos, bead_type, bead_charge in zip(
                            framed_sc_positions, sidechain_bead_types, sidechain_bead_charges
                        ):
                            delta = target[0] - bead_pos
                            dist_nm = max(1.0e-6, float(np.linalg.norm(delta)))
                            params = pair_params[(bead_type, target_type)]
                            sigma_nm = params["sigma_nm"]
                            epsilon_kj = params["epsilon_kj_mol"]
                            sr = sigma_nm / dist_nm
                            sr2 = sr * sr
                            sr6 = sr2 * sr2 * sr2
                            lj = 4.0 * epsilon_kj * (sr6 * sr6 - sr6)
                            coul = (
                                COULOMB_K_DRY_KJ_NM * bead_charge * target_charge / dist_nm
                                if bead_charge and target_charge
                                else 0.0
                            )
                            rot_energy += lj + coul
                            # Add WCA excluded area per-bead-pair.
                            contact_nm = (2.0 ** (1.0 / 6.0)) * sigma_nm
                            if dist_nm < contact_nm:
                                sr_wca = sigma_nm / max(dist_nm, 1e-6)
                                sr6_wca = sr_wca ** 6
                                rot_energy += 4.0 * DEFAULT_PRODUCTION_KBT_KJ_MOL * (sr6_wca * sr6_wca - sr6_wca) + DEFAULT_PRODUCTION_KBT_KJ_MOL
                    rot_energy_sum += rot_energy
                rot_energy = rot_energy_sum / float(len(bead_frame_angles))

                _accumulate_on_cos_grid(
                    cos_theta, rot_energy, cos_theta_grid,
                    rotamer_energy_sum[irot], rotamer_energy_weight[irot],
                )
                total_energy += rot_weight * rot_energy

            _accumulate_on_cos_grid(
                cos_theta, total_energy, cos_theta_grid,
                energy_sum, energy_weight,
            )

        for ia in range(n_angle):
            if energy_weight[ia] <= 0.0:
                raise RuntimeError(
                    f"Cos(theta) bin empty for {residue} target {target_type} at r={r_nm:.4f} nm"
                )
            angular_energy[ia][ir] = energy_sum[ia] / energy_weight[ia]
            for irot in range(n_rotamer):
                if rotamer_energy_weight[irot][ia] <= 0.0:
                    raise RuntimeError(
                        f"Rotamer cos(theta) bin empty for {residue} target {target_type} "
                        f"rotamer={irot} at r={r_nm:.4f} nm"
                    )
                rotamer_angular_energy[irot][ia][ir] = (
                    rotamer_energy_sum[irot][ia] / rotamer_energy_weight[irot][ia]
                )

    radial_energy, angular_profile, angular_radial_energy, rms_error = (
        _factorize_one_sided_orientation(angular_energy, cos_theta_grid)
    )

    # Attractive removal + background subtraction: radial and angular coupling
    # must be nonnegative. Negative values encode non-transferable two-body
    # attraction that produces unphysical many-body effects.
    radial_energy = [max(0.0, float(v)) for v in radial_energy]
    angular_radial_energy = [max(0.0, float(v)) for v in angular_radial_energy]

    rotamer_radial_energy = []
    rotamer_angular_profile = []
    rotamer_angular_radial_energy = []
    rotamer_rms_error = []
    for irot in range(n_rotamer):
        rr, rp, ra, rrm = _factorize_one_sided_orientation(
            rotamer_angular_energy[irot], cos_theta_grid
        )
        rotamer_radial_energy.append([max(0.0, float(v)) for v in rr])
        rotamer_angular_profile.append(rp)
        rotamer_angular_radial_energy.append([max(0.0, float(v)) for v in ra])
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
        "factorization_rms_error": rms_error,
        "sidechain_bead_frame_count": len(bead_frame_angles),
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
    fit_relax_steps: int = 0,
    sc_restraint_k: float = 5000.0,
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
                    "rel_relax_steps": int(fit_relax_steps),
                    "sc_restraint_k": float(sc_restraint_k),
                    "dist_min_nm": 0.10,
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

    g = h5.create_group("sc_table")
    g.attrs["schema"] = SCHEMA_SC
    g.attrs["forcefield_name"] = "martini22"
    g.attrs["fit_relax_steps"] = np.int32(fit_relax_steps)
    g.attrs["sidechain_bead_frame_count"] = np.int32(len(sidechain_bead_frame_angles))
    g.attrs["orientation_sampling"] = "target_direction_vector_grid"
    g.attrs["excluded_area_source"] = "wca_per_bead_pair_contact_kbt"
    g.attrs["attractive_control_source"] = "nontransferable_attraction_removed"
    g.attrs["isotropic_background_source"] = "attractive_radial_removed"
    g.attrs["relaxation"] = "lab_frame_position_restraints" if fit_relax_steps > 0 else "rigid_rotated_geometry"
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

    radial_controls = np.zeros((n_knot_radial, cos_theta_grid.size), dtype=np.float64)
    for ia in range(cos_theta_grid.size):
        radial_controls[:, ia] = _fit_radial_bspline(
            t_radial,
            values[:, ia] / float(energy_conversion),
            rad_knot_vector,
            smooth=smooth,
        )

    controls = np.zeros((n_knot_radial, n_knot_angular), dtype=np.float64)
    for ir in range(n_knot_radial):
        controls[ir, :] = _fit_angular_bspline(
            t_angular,
            radial_controls[ir, :],
            n_control=n_knot_angular,
            smooth=smooth,
        )
    return controls


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

    radial_controls = np.zeros(
        (n_knot_radial, cos_theta_grid.size, cos_theta_grid.size),
        dtype=np.float64,
    )
    for ia1 in range(cos_theta_grid.size):
        for ia2 in range(cos_theta_grid.size):
            radial_controls[:, ia1, ia2] = _fit_radial_bspline(
                t_radial,
                values[:, ia1, ia2] / float(energy_conversion),
                rad_knot_vector,
                smooth=smooth,
            )

    tmp = np.zeros((n_knot_radial, n_knot_angular, cos_theta_grid.size), dtype=np.float64)
    for ir in range(n_knot_radial):
        for ia2 in range(cos_theta_grid.size):
            tmp[ir, :, ia2] = _fit_angular_bspline(
                t_angular,
                radial_controls[ir, :, ia2],
                n_control=n_knot_angular,
                smooth=smooth,
            )

    controls = np.zeros((n_knot_radial, n_knot_angular, n_knot_angular), dtype=np.float64)
    for ir in range(n_knot_radial):
        for ia1 in range(n_knot_angular):
            controls[ir, ia1, :] = _fit_angular_bspline(
                t_angular,
                tmp[ir, ia1, :],
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
    "transverse_inertia_g_mol_a2",
    "head_tail_span_ang",
    "tail_projection_ang",
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
    group.attrs["derivation_schema"] = "dry_martini_dopc_derived_v1"
    for attr_name in _CG_DERIVED_NUMERIC_ATTRS:
        if attr_name in derived_params:
            group.attrs[attr_name] = np.float32(derived_params[attr_name])
    for attr_name in _CG_DERIVED_STRING_ATTRS:
        if attr_name in derived_params:
            group.attrs[attr_name] = derived_params[attr_name]

def _compute_lipid_bonded_energy(positions: np.ndarray) -> float:
    """Compute harmonic bond + cosine-based angle energy for ONE DOPC lipid.

    positions: (14, 3) float64 array of bead positions in nm.
    Returns energy in kJ/mol.
    """
    energy = 0.0
    for i, j, r0, k in _CURRENT_CG_BONDS:
        dr = positions[i] - positions[j]
        r = float(np.sqrt(np.dot(dr, dr)))
        energy += 0.5 * k * (r - r0) ** 2
    # Angles: V(theta) = 0.5 * k * (cos(theta) - cos(theta0))^2
    for i, j, k_idx, theta0_deg, k_ang in _CURRENT_CG_ANGLES:
        r_ij = positions[i] - positions[j]
        r_kj = positions[k_idx] - positions[j]
        d_ij = float(np.sqrt(np.dot(r_ij, r_ij)))
        d_kj = float(np.sqrt(np.dot(r_kj, r_kj)))
        if d_ij < 1e-8 or d_kj < 1e-8:
            continue
        cos_theta = float(np.dot(r_ij, r_kj)) / (d_ij * d_kj)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        cos_theta0 = math.cos(math.radians(theta0_deg))
        energy += 0.5 * float(k_ang) * (cos_theta - cos_theta0) ** 2
    return energy


def _compute_lipid_bonded_energy_and_gradient(
    positions: np.ndarray,
) -> tuple[float, np.ndarray]:
    """Compute harmonic bond + cosine angle energy and analytic gradient.

    positions: (14, 3) float64 array in any frame (rotation/translation invariant).
    Returns (energy_kj_mol, gradient).
    """
    n = positions.shape[0]
    grad = np.zeros_like(positions, dtype=np.float64)
    energy = 0.0

    for i, j, r0, k in _CURRENT_CG_BONDS:
        dr = positions[i] - positions[j]
        r = float(np.sqrt(np.dot(dr, dr)))
        if r < 1e-12:
            continue
        de = float(k) * (r - r0) / r
        grad[i] += de * dr
        grad[j] -= de * dr
        energy += 0.5 * float(k) * (r - r0) ** 2

    for i, j, k_idx, theta0_deg, k_ang in _CURRENT_CG_ANGLES:
        r_ij = positions[i] - positions[j]
        r_kj = positions[k_idx] - positions[j]
        d_ij = float(np.sqrt(np.dot(r_ij, r_ij)))
        d_kj = float(np.sqrt(np.dot(r_kj, r_kj)))
        if d_ij < 1e-8 or d_kj < 1e-8:
            continue
        cos_theta = float(np.dot(r_ij, r_kj)) / (d_ij * d_kj)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        cos_theta0 = math.cos(math.radians(theta0_deg))
        dcos = cos_theta - cos_theta0
        energy += 0.5 * float(k_ang) * dcos * dcos

        coeff = float(k_ang) * dcos / (d_ij * d_kj)
        grad_i = coeff * (r_kj - cos_theta * (d_kj / d_ij) * r_ij)
        grad_k = coeff * (r_ij - cos_theta * (d_ij / d_kj) * r_kj)
        grad[i] += grad_i
        grad[k_idx] += grad_k
        grad[j] -= (grad_i + grad_k)

    return energy, grad


def _project_com_direction_constraints(
    positions: np.ndarray,
    masses: np.ndarray,
    head_idx: int,
    tail_indices: list,
    target_dir: np.ndarray | None = None,
    max_iter: int = 5,
    fixed_bead_idx: int | None = None,
) -> np.ndarray:
    """Project positions to satisfy COM=0 and direction=target_dir constraints.

    positions: (N, 3) bead positions (body frame).
    masses: (N,) bead masses.
    head_idx: index of head bead (NC3 for DOPC, CB for SC).
    tail_indices: indices of tail beads (C5A/C5B for DOPC, all non-CB SC beads).
    target_dir: direction to align to (default z-axis).
    fixed_bead_idx: if set, pin this bead to origin instead of subtracting COM.
        Use for SC beads where CB (index 0) is the anchor point.

    Applies position constraints and Rodrigues rotation iteratively.
    """
    if target_dir is None:
        target_dir = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    else:
        target_dir = np.asarray(target_dir, dtype=np.float64)
        target_dir = target_dir / float(np.linalg.norm(target_dir))

    pos = np.asarray(positions, dtype=np.float64).copy()
    for _ in range(max_iter):
        if fixed_bead_idx is not None:
            # Pin the anchor bead to origin
            anchor_pos = pos[fixed_bead_idx].copy()
            pos = pos - anchor_pos[None, :]
        else:
            # COM constraint
            weighted_center = np.average(pos, axis=0, weights=masses)
            pos = pos - weighted_center[None, :]

        # Direction constraint: rotate (tail_mid - head) to target_dir
        head = pos[head_idx]
        tail_indices_arr = np.asarray(tail_indices, dtype=np.intp)
        tail_mid = np.mean(pos[tail_indices_arr], axis=0)
        current_dir = tail_mid - head
        current_norm = float(np.linalg.norm(current_dir))
        if current_norm < 1e-12:
            break
        current_dir = current_dir / current_norm

        v = np.cross(current_dir, target_dir)
        s = float(np.linalg.norm(v))
        c = float(np.dot(current_dir, target_dir))
        if s < 1e-12:
            break
        v = v / s
        # Rodrigues rotation about origin
        pos = (
            pos * c
            + np.cross(v[None, :], pos) * s
            + v[None, :] * np.sum(v[None, :] * pos, axis=1, keepdims=True) * (1.0 - c)
        )
    return pos


def _clamp_bead_radial(
    lab_positions: np.ndarray,
    ref_body: np.ndarray,
    head_idx: int = 0,
    tail_indices: tuple = (8, 13),
) -> None:
    """Clamp each bead's radial distance from the direction axis to reference.

    In a dense bilayer, neighboring lipids occupy the surrounding cylindrical
    volume. Beads cannot expand radially outward beyond their reference
    envelope — this is steric exclusion by neighbors, not a fitted restraint.
    Beads may contract inward but not expand outward.

    This is applied in-place to lab_positions.
    """
    tail_arr = np.asarray(tail_indices, dtype=np.intp)
    head = lab_positions[head_idx]
    tail_mid = np.mean(lab_positions[tail_arr], axis=0)
    direction = tail_mid - head
    dir_norm = float(np.linalg.norm(direction))
    if dir_norm < 1e-12:
        return
    direction /= dir_norm

    com = np.mean(lab_positions, axis=0)
    delta = lab_positions - com
    axial = np.dot(delta, direction)
    perp = delta - np.outer(axial, direction)
    r_perp = np.sqrt(np.sum(perp ** 2, axis=1))
    ref_r = np.sqrt(np.asarray(ref_body, dtype=np.float64)[:, 0] ** 2
                    + np.asarray(ref_body, dtype=np.float64)[:, 1] ** 2)

    exceed = (r_perp > ref_r) & (r_perp > 1e-12)
    if np.any(exceed):
        scale = ref_r[exceed] / r_perp[exceed]
        perp[exceed] *= scale[:, None]
        lab_positions[:] = com + np.outer(axial, direction) + perp


def _compute_position_restraint_energy_and_gradient(
    positions: np.ndarray,
    ref_positions: np.ndarray,
    k_restraint: float,
) -> tuple[float, np.ndarray]:
    """Harmonic position restraint: E = 0.5 * k * sum (pos - ref)^2."""
    delta = np.asarray(positions, dtype=np.float64) - np.asarray(ref_positions, dtype=np.float64)
    energy = 0.5 * float(k_restraint) * float(np.sum(delta * delta))
    gradient = float(k_restraint) * delta
    return energy, gradient


def _dopc_body_to_lab(
    body_positions: np.ndarray,
    direction_lab: np.ndarray,
    com_lab: np.ndarray,
) -> np.ndarray:
    """Transform DOPC bead positions from body frame (COM=0, dir=z) to lab frame."""
    R = _rotation_to_align_z_np(np.asarray(direction_lab, dtype=np.float64))
    return com_lab[None, :] + (R @ np.asarray(body_positions, dtype=np.float64).T).T


def _dopc_lab_to_body(
    lab_positions: np.ndarray,
    direction_lab: np.ndarray,
    com_lab: np.ndarray,
) -> np.ndarray:
    """Transform DOPC bead positions from lab frame to body frame (COM=0, dir=z)."""
    R = _rotation_to_align_z_np(np.asarray(direction_lab, dtype=np.float64))
    return ((np.asarray(lab_positions, dtype=np.float64) - com_lab[None, :]) @ R)


def _relax_dopc_beads(
    init_body_positions: np.ndarray,
    ref_body_positions: np.ndarray,
    bead_types: list,
    bead_charges: list,
    bead_masses: np.ndarray,
    pair_params: dict,
    direction_lab: np.ndarray,
    com_lab: np.ndarray,
    partner_lab_positions: np.ndarray,
    partner_bead_types: list,
    partner_bead_charges: list,
    dist_min_nm: float = 0.20,
    rel_relax_steps: int = 50,
    step_size_nm: float = 0.0005,
    grad_tolerance: float = 1e-4,
) -> tuple[np.ndarray, float]:
    """Relax DOPC bead positions with COM and direction constraints.

    The DOPC beads relax against fixed partner positions (point particle,
    SC beads, or CGLipid beads). The DOPC COM and direction vector are held
    fixed at their lab-frame values. Internal energy includes bonded terms
    (bonds + angles) and nonbonded with the partner.

    Returns (relaxed_body_positions, final_energy_kj_mol).
    """
    if rel_relax_steps <= 0:
        lab_pos = _dopc_body_to_lab(init_body_positions, direction_lab, com_lab)
        energy, _, _ = _compute_pair_energy_and_gradient(
            lab_pos, partner_lab_positions,
            bead_types, partner_bead_types,
            bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=dist_min_nm,
        )
        return init_body_positions.copy(), energy

    body_pos = np.asarray(init_body_positions, dtype=np.float64).copy()
    step = float(step_size_nm)
    head_idx = 0
    tail_indices = [8, 13]

    for _ in range(rel_relax_steps):
        lab_pos = _dopc_body_to_lab(body_pos, direction_lab, com_lab)

        # Nonbonded energy + gradient with partner
        nb_energy, grad_lab, _ = _compute_pair_energy_and_gradient(
            lab_pos, partner_lab_positions,
            bead_types, partner_bead_types,
            bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=dist_min_nm,
        )

        # Transform nonbonded gradient to body frame
        R = _rotation_to_align_z_np(np.asarray(direction_lab, dtype=np.float64))
        grad_body_nb = (grad_lab @ R)

        # Bonded energy + gradient (finite differences)
        bonded_energy, grad_body_bonded = _bonded_gradient_body(
            body_pos, bead_masses, head_idx, tail_indices
        )

        total_grad = grad_body_nb + grad_body_bonded
        max_grad = float(np.max(np.abs(total_grad)))
        if max_grad < grad_tolerance:
            break

        # Backtracking line search
        accepted = False
        for _ in range(3):
            new_body = body_pos - step * total_grad
            new_body = _project_com_direction_constraints(
                new_body, bead_masses, head_idx, tail_indices
            )
            new_lab = _dopc_body_to_lab(new_body, direction_lab, com_lab)
            new_nb, _, _ = _compute_pair_energy_and_gradient(
                new_lab, partner_lab_positions,
                bead_types, partner_bead_types,
                bead_charges, partner_bead_charges,
                pair_params, dist_min_nm=dist_min_nm,
            )
            new_bonded = _compute_lipid_bonded_energy(new_body)
            if new_nb + new_bonded <= nb_energy + bonded_energy + 1e-10:
                body_pos = new_body
                step = min(step * 1.1, float(step_size_nm) * 2.0)
                accepted = True
                break
            step *= 0.5
        if not accepted:
            body_pos = _project_com_direction_constraints(
                body_pos - step * total_grad, bead_masses, head_idx, tail_indices
            )
            step = float(step_size_nm)

    # Final energy from relaxed positions
    lab_pos = _dopc_body_to_lab(body_pos, direction_lab, com_lab)
    final_energy, _, _ = _compute_pair_energy_and_gradient(
        lab_pos, partner_lab_positions,
        bead_types, partner_bead_types,
        bead_charges, partner_bead_charges,
        pair_params, dist_min_nm=dist_min_nm,
    )
    return body_pos, final_energy


def _relax_dopc_lab(
    init_body_positions: np.ndarray,
    ref_body_positions: np.ndarray,
    bead_types: list,
    bead_charges: list,
    bead_masses: np.ndarray,
    pair_params: dict,
    direction_lab: np.ndarray,
    com_lab: np.ndarray,
    partner_lab_positions: np.ndarray,
    partner_bead_types: list,
    partner_bead_charges: list,
    dist_min_nm: float = 0.10,
    rel_relax_steps: int = 50,
    step_size_nm: float = 0.0005,
    grad_tolerance: float = 1e-4,
    com_restraint_k: float = 500.0,
    dir_restraint_k: float = 500.0,
    plane_constraint: bool = False,
) -> tuple[np.ndarray, float]:
    """Relax a single DOPC molecule in the lab frame with soft restraints.

    Unlike _relax_dopc_beads which enforces hard COM and direction constraints
    via iterative projection, this function uses harmonic restraints on the
    COM and head-to-tail vector. This allows physical bead deformation while
    keeping the molecule near the sampled geometry.

    Returns (relaxed_body_positions, final_energy_kj_mol).
    """
    head_idx = 0
    tail_indices = [8, 13]
    tail_arr = np.asarray(tail_indices, dtype=np.intp)
    n_beads = len(bead_types)

    # Compute reference span vector in lab frame from the body-frame geometry.
    ref_head_body = np.asarray(ref_body_positions[head_idx], dtype=np.float64)
    ref_tail_mid_body = np.mean(np.asarray(ref_body_positions, dtype=np.float64)[tail_arr], axis=0)
    ref_span = float(np.linalg.norm(ref_tail_mid_body - ref_head_body))
    ref_span_vec_lab = np.asarray(direction_lab, dtype=np.float64) * ref_span

    # Place beads at initial positions in lab frame.
    pos = _dopc_body_to_lab(init_body_positions, direction_lab, com_lab)

    if rel_relax_steps <= 0:
        energy, _, _ = _compute_pair_energy_and_gradient(
            pos, partner_lab_positions,
            bead_types, partner_bead_types,
            bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=dist_min_nm,
        )
        return init_body_positions.copy(), energy

    ref_bonded = _compute_lipid_bonded_energy(init_body_positions)
    step = float(step_size_nm)
    inv_n = 1.0 / float(n_beads)
    inv_n_tails = 1.0 / float(len(tail_indices))

    for _ in range(rel_relax_steps):
        nb_energy, grad_nb, _ = _compute_pair_energy_and_gradient(
            pos, partner_lab_positions,
            bead_types, partner_bead_types,
            bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=0.01,
        )

        bond_energy, grad_bond = _compute_lipid_bonded_energy_and_gradient(pos)

        # COM restraint
        current_com = np.mean(pos, axis=0)
        com_delta = current_com - np.asarray(com_lab, dtype=np.float64)
        com_energy = 0.5 * float(com_restraint_k) * float(np.dot(com_delta, com_delta))
        com_grad_per_dim = float(com_restraint_k) * com_delta * inv_n

        # Head-to-tail vector restraint
        head_pos = pos[head_idx]
        tail_mid = np.mean(pos[tail_arr], axis=0)
        current_span_vec = tail_mid - head_pos
        span_delta = current_span_vec - ref_span_vec_lab
        dir_energy = 0.5 * float(dir_restraint_k) * float(np.dot(span_delta, span_delta))
        dir_coeff = float(dir_restraint_k) * span_delta

        total_grad = grad_nb + grad_bond
        total_grad += np.full_like(pos, com_grad_per_dim)
        total_grad[head_idx] -= dir_coeff
        for ti in tail_indices:
            total_grad[ti] += dir_coeff * inv_n_tails

        total_energy = nb_energy + bond_energy + com_energy + dir_energy

        # Gradient clipping: max 0.005 nm per bead per step.
        disp = step * total_grad
        max_disp = float(np.max(np.sqrt(np.sum(disp ** 2, axis=1))))
        if max_disp > 0.005:
            scale = 0.005 / max_disp
            total_grad *= scale

        max_grad = float(np.max(np.abs(total_grad)))
        if max_grad < grad_tolerance:
            break

        # Backtracking line search in lab frame.
        accepted = False
        for _ in range(3):
            new_pos = pos - step * total_grad
            if plane_constraint:
                _clamp_bead_radial(new_pos, ref_body_positions, head_idx=0, tail_indices=(8, 13))
            new_nb, _, _ = _compute_pair_energy_and_gradient(
                new_pos, partner_lab_positions,
                bead_types, partner_bead_types,
                bead_charges, partner_bead_charges,
                pair_params, dist_min_nm=0.01,
            )
            new_bond = _compute_lipid_bonded_energy(new_pos)
            new_com = np.mean(new_pos, axis=0)
            new_com_delta = new_com - np.asarray(com_lab, dtype=np.float64)
            new_com_e = 0.5 * float(com_restraint_k) * float(np.dot(new_com_delta, new_com_delta))
            new_head = new_pos[head_idx]
            new_tail_mid = np.mean(new_pos[tail_arr], axis=0)
            new_span_delta = (new_tail_mid - new_head) - ref_span_vec_lab
            new_dir_e = 0.5 * float(dir_restraint_k) * float(np.dot(new_span_delta, new_span_delta))
            if new_nb + new_bond + new_com_e + new_dir_e <= total_energy + 1e-10:
                pos = new_pos
                step = min(step * 1.1, float(step_size_nm) * 2.0)
                accepted = True
                break
            step *= 0.5
        if not accepted:
            pos -= step * total_grad
            if plane_constraint:
                _clamp_bead_radial(pos, ref_body_positions, head_idx=0, tail_indices=(8, 13))
            step = float(step_size_nm)

    # Final energy: physical LJ + bonded deformation cost.
    final_nb, _, _ = _compute_pair_energy_and_gradient(
        pos, partner_lab_positions,
        bead_types, partner_bead_types,
        bead_charges, partner_bead_charges,
        pair_params, dist_min_nm=0.10,
        cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
    )
    final_body = _dopc_lab_to_body(pos, direction_lab, com_lab)
    final_bonded = _compute_lipid_bonded_energy(pos)
    deformation_cost = final_bonded - ref_bonded
    return final_body, final_nb + deformation_cost


def _bonded_gradient_body(
    body_positions: np.ndarray,
    masses: np.ndarray,
    head_idx: int,
    tail_indices: list,
    eps: float = 1e-6,
    fixed_bead_idx: int | None = None,
) -> tuple[float, np.ndarray]:
    """Compute bonded energy gradient via central finite differences.

    The gradient is computed in the unconstrained space with constraint
    projection applied to each finite-difference sample so that the
    bonded gradient does not fight the geometric constraints.
    """
    n_beads = body_positions.shape[0]
    grad = np.zeros_like(body_positions)
    for i in range(n_beads):
        for d in range(3):
            pos_plus = body_positions.copy()
            pos_minus = body_positions.copy()
            pos_plus[i, d] += eps
            pos_minus[i, d] -= eps
            pos_plus = _project_com_direction_constraints(
                pos_plus, masses, head_idx, tail_indices,
                fixed_bead_idx=fixed_bead_idx,
            )
            pos_minus = _project_com_direction_constraints(
                pos_minus, masses, head_idx, tail_indices,
                fixed_bead_idx=fixed_bead_idx,
            )
            e_plus = _compute_lipid_bonded_energy(pos_plus)
            e_minus = _compute_lipid_bonded_energy(pos_minus)
            grad[i, d] = (e_plus - e_minus) / (2.0 * eps)
    energy = _compute_lipid_bonded_energy(body_positions)
    return energy, grad


def _relax_dopc_pair(
    init_body1: np.ndarray,
    init_body2: np.ndarray,
    ref_body: np.ndarray,
    bead_types: list,
    bead_charges: list,
    bead_masses: np.ndarray,
    pair_params: dict,
    dir1_lab: np.ndarray,
    dir2_lab: np.ndarray,
    r_nm: float,
    dist_min_nm: float = 0.20,
    rel_relax_steps: int = 50,
    step_size_nm: float = 0.0005,
    grad_tolerance: float = 1e-4,
    soft_core_alpha: float = 0.0,
    plane_constraint: bool = False,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Relax two DOPC molecules in lab frame with bonded deformation cost.

    Both lipids are free to translate, rotate, and deform internally during
    relaxation. The bonded deformation penalty (relative to reference bead
    geometry) resists excessive distortion. This matches the original
    _relax_lipid_pair approach and produces a physically softer effective
    pair potential than rigid COM/direction constraints.

    When soft_core_alpha > 0, uses soft-core LJ during relaxation to
    prevent bead explosion from steric clashes. Final energy is always
    evaluated with the physical LJ (alpha=0) and includes the bonded
    deformation cost.

    Returns (relaxed_body1, relaxed_body2, final_energy_kj_mol) where
    body positions are in the body frame (COM=0, dir=z).
    """
    com1_lab = np.zeros(3, dtype=np.float64)
    com2_lab = np.array([r_nm, 0.0, 0.0], dtype=np.float64)

    # Place beads at rigid reference positions in lab frame.
    pos1 = _dopc_body_to_lab(init_body1, dir1_lab, com1_lab)
    pos2 = _dopc_body_to_lab(init_body2, dir2_lab, com2_lab)

    if rel_relax_steps <= 0:
        energy, _, _ = _compute_pair_energy_and_gradient(
            pos1, pos2, bead_types, bead_types,
            bead_charges, bead_charges, pair_params,
            dist_min_nm=dist_min_nm,
        )
        return init_body1.copy(), init_body2.copy(), energy

    p1 = pos1.copy()
    p2 = pos2.copy()
    step = float(step_size_nm)
    ref_bond1 = _compute_lipid_bonded_energy(init_body1)
    ref_bond2 = _compute_lipid_bonded_energy(init_body2)
    head_idx = 0
    tail_indices = [8, 13]

    for _ in range(rel_relax_steps):
        nb_energy, grad1_nb, grad2_nb = _compute_pair_energy_and_gradient(
            p1, p2, bead_types, bead_types,
            bead_charges, bead_charges, pair_params,
            dist_min_nm=0.01,
            soft_core_alpha=soft_core_alpha,
        )

        # Bonded energy + gradient in lab frame (analytic, rotation/translation invariant).
        bond1, grad1_b = _compute_lipid_bonded_energy_and_gradient(p1)
        bond2, grad2_b = _compute_lipid_bonded_energy_and_gradient(p2)

        grad1 = grad1_nb + grad1_b
        grad2 = grad2_nb + grad2_b

        # Gradient clipping: max 0.005 nm per bead per step (matching original).
        d1_norms = np.sqrt(np.sum((step * grad1) ** 2, axis=1))
        d2_norms = np.sqrt(np.sum((step * grad2) ** 2, axis=1))
        max_disp = max(float(np.max(d1_norms)), float(np.max(d2_norms)), 1e-15)
        if max_disp > 0.005:
            scale = 0.005 / max_disp
            grad1 *= scale
            grad2 *= scale

        max_grad = max(float(np.max(np.abs(grad1))), float(np.max(np.abs(grad2))))
        if max_grad < grad_tolerance:
            break

        # Line search in lab frame.
        accepted = False
        for _ in range(3):
            new_p1 = p1 - step * grad1
            new_p2 = p2 - step * grad2
            if plane_constraint:
                _clamp_bead_radial(new_p1, ref_body, head_idx=0, tail_indices=(8, 13))
                _clamp_bead_radial(new_p2, ref_body, head_idx=0, tail_indices=(8, 13))
            new_nb, _, _ = _compute_pair_energy_and_gradient(
                new_p1, new_p2, bead_types, bead_types,
                bead_charges, bead_charges, pair_params,
                dist_min_nm=0.01,
                soft_core_alpha=soft_core_alpha,
            )
            new_b1 = _compute_lipid_bonded_energy(new_p1)
            new_b2 = _compute_lipid_bonded_energy(new_p2)
            if new_nb + new_b1 + new_b2 <= nb_energy + bond1 + bond2 + 1e-10:
                p1, p2 = new_p1, new_p2
                step = min(step * 1.1, float(step_size_nm) * 2.0)
                accepted = True
                break
            step *= 0.5
        if not accepted:
            p1 -= step * grad1
            p2 -= step * grad2
            if plane_constraint:
                _clamp_bead_radial(p1, ref_body, head_idx=0, tail_indices=(8, 13))
                _clamp_bead_radial(p2, ref_body, head_idx=0, tail_indices=(8, 13))
            step = float(step_size_nm)

    # Final energy: physical LJ (no soft-core) + bonded deformation cost.
    final_nb, _, _ = _compute_pair_energy_and_gradient(
        p1, p2, bead_types, bead_types,
        bead_charges, bead_charges, pair_params,
        dist_min_nm=0.10,
        soft_core_alpha=0.0,
        cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
    )
    final_body1 = _dopc_lab_to_body(p1, dir1_lab, com1_lab)
    final_body2 = _dopc_lab_to_body(p2, dir2_lab, com2_lab)
    final_bond1 = _compute_lipid_bonded_energy(p1)
    final_bond2 = _compute_lipid_bonded_energy(p2)
    deformation_cost = (final_bond1 - ref_bond1) + (final_bond2 - ref_bond2)
    return final_body1, final_body2, final_nb + deformation_cost


def _relax_sc_beads(
    init_body_positions: np.ndarray,
    ref_body_positions: np.ndarray,
    sc_bead_types: list,
    sc_bead_charges: list,
    k_restraint: float,
    pair_params: dict,
    direction_lab: np.ndarray,
    cb_anchor_lab: np.ndarray,
    partner_lab_positions: np.ndarray,
    partner_bead_types: list,
    partner_bead_charges: list,
    dist_min_nm: float = 0.10,
    rel_relax_steps: int = 50,
    step_size_nm: float = 0.0005,
    grad_tolerance: float = 1e-4,
) -> tuple[np.ndarray, float]:
    """Relax SC beads in lab frame with harmonic position restraints.

    SC beads relax in the lab frame against fixed partner positions. The
    internal energy is a harmonic position restraint to the reference
    lab-frame positions (derived from the body-frame reference and the
    fixed CB anchor + direction). No hard geometric constraints are applied.

    Returns (relaxed_body_positions, final_energy_kj_mol) where final energy
    includes the nonbonded interaction plus the deformation cost (restraint
    energy relative to the reference).
    """
    n_beads = len(sc_bead_types)
    dir_lab = np.asarray(direction_lab, dtype=np.float64)
    anchor_lab = np.asarray(cb_anchor_lab, dtype=np.float64)
    ref_lab = _dopc_body_to_lab(np.asarray(ref_body_positions, dtype=np.float64), dir_lab, anchor_lab)
    init_lab = _dopc_body_to_lab(np.asarray(init_body_positions, dtype=np.float64), dir_lab, anchor_lab)

    if rel_relax_steps <= 0 or n_beads <= 1:
        energy, _, _ = _compute_pair_energy_and_gradient(
            init_lab, partner_lab_positions,
            sc_bead_types, partner_bead_types,
            sc_bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=dist_min_nm,
        )
        return init_body_positions.copy(), energy

    pos = init_lab.copy()
    step = float(step_size_nm)

    for _ in range(rel_relax_steps):
        nb_energy, grad_nb, _ = _compute_pair_energy_and_gradient(
            pos, partner_lab_positions,
            sc_bead_types, partner_bead_types,
            sc_bead_charges, partner_bead_charges,
            pair_params, dist_min_nm=0.01,
        )

        restraint_energy, grad_restraint = _compute_position_restraint_energy_and_gradient(
            pos, ref_lab, k_restraint
        )

        total_grad = grad_nb + grad_restraint

        # Gradient clipping: max 0.005 nm per bead per step.
        disp = step * total_grad
        max_disp = float(np.max(np.sqrt(np.sum(disp ** 2, axis=1))))
        if max_disp > 0.005:
            total_grad *= 0.005 / max_disp

        max_grad = float(np.max(np.abs(total_grad)))
        if max_grad < grad_tolerance:
            break

        total_energy = nb_energy + restraint_energy

        accepted = False
        for _ in range(3):
            new_pos = pos - step * total_grad
            new_nb, _, _ = _compute_pair_energy_and_gradient(
                new_pos, partner_lab_positions,
                sc_bead_types, partner_bead_types,
                sc_bead_charges, partner_bead_charges,
                pair_params, dist_min_nm=0.01,
            )
            new_re, _ = _compute_position_restraint_energy_and_gradient(
                new_pos, ref_lab, k_restraint
            )
            if new_nb + new_re <= total_energy + 1e-10:
                pos = new_pos
                step = min(step * 1.1, float(step_size_nm) * 2.0)
                accepted = True
                break
            step *= 0.5
        if not accepted:
            pos -= step * total_grad
            step = float(step_size_nm)

    final_nb, _, _ = _compute_pair_energy_and_gradient(
        pos, partner_lab_positions,
        sc_bead_types, partner_bead_types,
        sc_bead_charges, partner_bead_charges,
        pair_params, dist_min_nm=0.10,
    )
    final_re, _ = _compute_position_restraint_energy_and_gradient(pos, ref_lab, k_restraint)
    final_body = _dopc_lab_to_body(pos, dir_lab, anchor_lab)
    return final_body, final_nb + final_re


def _compute_pair_energy_and_gradient(
    pos1: np.ndarray,
    pos2: np.ndarray,
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
    dist_min_nm: float = 0.20,
    soft_core_alpha: float = 0.0,
    cutoff_nm: float | None = DRY_MARTINI_NONBONDED_CUTOFF_NM,
) -> tuple:
    """Compute non-bonded (LJ+Coulomb) energy and per-bead gradients for a pair of lipids.

    pos1, pos2: (14, 3) float64 arrays in nm.

    If soft_core_alpha > 0, uses GROMACS-style soft-core LJ:
      r_eff = (r^6 + alpha * sigma^6)^(1/6)
    This smoothly caps the repulsive energy at short distances, preventing
    catastrophic LJ energies when beads overlap.

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

            if soft_core_alpha > 0.0 and dist < 2.0 * sigma_nm:
                # Soft-core LJ: r_eff^6 = r^6 + alpha * sigma^6
                sc_offset = float(soft_core_alpha) * (sigma_nm ** 6)
                r6_eff = dist_sq ** 3 + sc_offset
                r_eff = r6_eff ** (1.0 / 6.0)
                inv_r6 = 1.0 / r6_eff
                inv_r12 = inv_r6 * inv_r6
                lj = 4.0 * epsilon_kj * (sigma_nm ** 12 * inv_r12 - sigma_nm ** 6 * inv_r6)
                total_energy += lj

                if dist > 1e-10:
                    # d(r_eff)/d(dist) = dist^5 / r_eff^5
                    dr_eff_dr = (dist ** 5) / (r_eff ** 5)
                    dlj_dr_eff = 4.0 * epsilon_kj * (
                        -12.0 * sigma_nm ** 12 / (r_eff ** 13)
                        + 6.0 * sigma_nm ** 6 / (r_eff ** 7)
                    )
                    dlj_dr = dlj_dr_eff * dr_eff_dr
                    factor = dlj_dr / dist
                    gx = factor * dx
                    grad1[i] -= gx
                    grad2[j] += gx
            else:
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
    dist_min_nm: float = 0.20,
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
        dist_min_nm=dist_min_nm, soft_core_alpha=0.0,
        cutoff_nm=DRY_MARTINI_NONBONDED_CUTOFF_NM,
    )
    return total_energy


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


def _run_cg_pair_tensor_task(task: dict) -> tuple[int, int, int, float]:
    ir = int(task["ir"])
    ia1 = int(task["ia1"])
    ia2 = int(task["ia2"])
    r_nm = float(task["r_nm"])
    dirs1 = task["dirs1"]
    dirs2 = task["dirs2"]
    bead_frame_angles = task["bead_frame_angles"]
    ref_nm = task["ref_nm"]
    rel_relax_steps = int(task.get("rel_relax_steps", 0))
    soft_core_alpha = float(task.get("soft_core_alpha", 0.0))
    plane_constraint = bool(task.get("plane_constraint", True))
    sample_values = []
    for dir1 in dirs1:
        for dir2 in dirs2:
            for frame_angle1 in bead_frame_angles:
                for frame_angle2 in bead_frame_angles:
                    if rel_relax_steps > 0:
                        z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
                        rot1 = _rotation_about_axis_np(z_axis, float(frame_angle1))
                        rot2 = _rotation_about_axis_np(z_axis, float(frame_angle2))
                        init_body1 = (rot1 @ ref_nm.T).T
                        init_body2 = (rot2 @ ref_nm.T).T
                        _, _, energy = _relax_dopc_pair(
                            init_body1=init_body1,
                            init_body2=init_body2,
                            ref_body=ref_nm,
                            bead_types=task["bead_types"],
                            bead_charges=task["bead_charges"],
                            bead_masses=task.get("bead_masses"),
                            pair_params=task["pair_params"],
                            dir1_lab=np.asarray(dir1, dtype=np.float64),
                            dir2_lab=np.asarray(dir2, dtype=np.float64),
                            r_nm=r_nm,
                            dist_min_nm=float(task["dist_min_nm"]),
                            rel_relax_steps=rel_relax_steps,
                            soft_core_alpha=soft_core_alpha,
                            plane_constraint=plane_constraint,
                        )
                        sample_values.append(energy)
                    else:
                        sample_values.append(
                            _compute_cg_pair_energy(
                                r_nm,
                                np.asarray(dir1, dtype=np.float64),
                                np.asarray(dir2, dtype=np.float64),
                                float(frame_angle1),
                                float(frame_angle2),
                                ref_nm,
                                ref_nm,
                                task["bead_types"],
                                task["bead_types"],
                                task["bead_charges"],
                                task["bead_charges"],
                                task["pair_params"],
                                dist_min_nm=float(task["dist_min_nm"]),
                            )
                        )
    return ir, ia1, ia2, float(np.mean(sample_values))


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
    dist_min_nm: float = 0.25,
    knot_spacing_ang: float = 1.4,
    excluded_area_contact_nm: float | None = None,
    n_modes: int = 4,
    n_knot_radial: int = 14,
    n_knot_angular: int = 15,
    cg_smooth: float = 0.01,
    rel_relax_steps: int = 0,
    bead_masses: list | None = None,
    relax_soft_core_alpha: float = 0.0,
    temperature: float = 0.0,
    plane_constraint: bool = False,
) -> dict:
    """Fit full tensor-product B-spline parameters for CG lipid-CG lipid interactions.

    When temperature > 0, the azimuthal average is computed as a Boltzmann-weighted
    free energy (potential of mean force) at the given temperature in T_up units.
    This correctly suppresses steric-clash configurations that have negligible
    Boltzmann weight, rather than letting one clash dominate the arithmetic mean.
    """
    r_values = np.asarray(_linspace(r_min_nm, r_max_nm, r_count), dtype=np.float64)
    cos_theta_grid = np.asarray(_linspace(-1.0, 1.0, cos_theta_count), dtype=np.float64)
    n_angle = cos_theta_count
    n_radial = r_count

    phi_values = np.linspace(0.0, 2.0 * np.pi, azimuthal_count, endpoint=False)
    n12_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    dirs1 = _directions_with_dot_np(-n12_axis, cos_theta_grid, phi_values)
    dirs2 = _directions_with_dot_np(n12_axis, cos_theta_grid, phi_values)

    energy_grid = np.zeros((n_radial, n_angle, n_angle), dtype=np.float64)
    ref_nm = _canonicalize_lipid_reference_to_z(ref_bead_positions_nm)

    bead_frame_angles = _bead_frame_angles(bead_frame_count)
    total_samples = (
        n_radial * n_angle * n_angle
        * azimuthal_count * azimuthal_count
        * len(bead_frame_angles) * len(bead_frame_angles)
    )
    print(f"  Sampling CG-CG energy: {n_radial} radial x {n_angle}^2 angular "
          f"x {azimuthal_count}^2 azimuthal x {len(bead_frame_angles)}^2 bead-frame "
          f"= {total_samples} samples, direct rotated geometry, full tensor")

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
                        "rel_relax_steps": int(rel_relax_steps),
                        "soft_core_alpha": float(relax_soft_core_alpha),
                        "temperature": float(temperature),
                        "plane_constraint": bool(plane_constraint),
                        "bead_masses": (
                            np.asarray(bead_masses, dtype=np.float64)
                            if bead_masses is not None
                            else None
                        ),
                    }
                )
    for ir, ia1, ia2, energy in _parallel_map_ordered("CG-CG table", _run_cg_pair_tensor_task, tasks):
        r_nm_val = float(r_values[ir])
        # WCA excluded area: add repulsive energy below contact distance.
        # This ensures the sampled energies reflect physical volume exclusion.
        if excluded_area_contact_nm is not None and r_nm_val < float(excluded_area_contact_nm):
            sigma_nm = float(excluded_area_contact_nm) / (2.0 ** (1.0 / 6.0))
            sr = sigma_nm / max(r_nm_val, 1e-6)
            sr6 = sr ** 6
            energy += 4.0 * DEFAULT_PRODUCTION_KBT_KJ_MOL * (sr6 * sr6 - sr6) + DEFAULT_PRODUCTION_KBT_KJ_MOL
        energy_grid[ir, ia1, ia2] = energy

    # Isotropic background: subtract mean attractive contribution so the
    # table captures only orientation-dependent effects.
    radial_mean = energy_grid.mean(axis=(1, 2))
    attractive_background = np.minimum(radial_mean, 0.0)
    if np.any(attractive_background < 0.0):
        energy_grid = energy_grid - attractive_background[:, None, None]

    knot_spacing = float(knot_spacing_ang)
    tensor = _fit_radial_angular_angular_tensor_bspline(
        r_values,
        cos_theta_grid,
        energy_grid,
        n_knot_radial=n_knot_radial,
        n_knot_angular=n_knot_angular,
        knot_spacing_ang=knot_spacing,
        smooth=cg_smooth,
    )
    r_min_ang = float(r_min_nm) * LENGTH_CONVERSION_A_PER_NM
    n_core = max(0, min(n_knot_radial - 1, int(math.ceil((r_min_ang - 1e-6) / knot_spacing))))
    short_range_core_kj_mol = float(np.max(energy_grid[0]))
    short_range_core_eup = short_range_core_kj_mol / ENERGY_CONVERSION_KJ_PER_EUP
    if n_core > 0:
        tensor[:n_core, :, :] = np.maximum(tensor[:n_core, :, :], short_range_core_eup)
    control_min = float(np.min(energy_grid)) / ENERGY_CONVERSION_KJ_PER_EUP
    control_max = float(np.max(energy_grid)) / ENERGY_CONVERSION_KJ_PER_EUP
    tensor = np.clip(tensor, control_min, control_max)
    excluded_area_rows = 0
    if excluded_area_contact_nm is not None and excluded_area_contact_nm > 0.0:
        excluded_area_rows = max(
            0,
            min(
                n_knot_radial,
                int(math.ceil(float(excluded_area_contact_nm) * LENGTH_CONVERSION_A_PER_NM / knot_spacing)),
            ),
        )
        if excluded_area_rows > 0:
            tensor[:excluded_area_rows, :, :] = np.maximum(tensor[:excluded_area_rows, :, :], 0.0)
    attractive_control_count = int(np.count_nonzero(tensor < 0.0))
    if attractive_control_count:
        tensor = np.maximum(tensor, 0.0)
    rms_error = 0.0
    interaction_param = tensor.reshape(-1)

    return {
        "tensor_knots": tensor,
        "interaction_param": interaction_param,
        "rms_error": rms_error,
        "energy_grid_raw": energy_grid,
        "attractive_radial_background_kj_mol": attractive_background,
        "azimuthal_average": "energy_expectation" if temperature <= 0.0 else "boltzmann_free_energy",
        "isotropic_background_source": "attractive_radial_angular_mean_subtracted" if np.any(attractive_background < 0.0) else "none_full_resolved_dry_martini_pair_table",
        "excluded_area_source": "wca_dopc_contact_kbt" if excluded_area_rows > 0 else "",
        "excluded_area_nonnegative_rows": excluded_area_rows,
        "attractive_control_source": "nontransferable_many_neighbor_cgl_cgl_attraction_removed" if attractive_control_count > 0 else "retained_full_resolved_dry_martini_pair_table",
        "attractive_control_count": attractive_control_count,
        "core_boundary_source": "max_first_sampled_dry_martini_energy_expectation",
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
        "azimuthal_count": int(azimuthal_count),
        "bead_frame_count": int(len(bead_frame_angles)),
        "schema": "cg_lipid_pair_full_v1",
        "rel_relax_steps": int(rel_relax_steps),
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
    knot_spacing_ang: float = 1.4,
    excluded_area_contact_nm: float | None = None,
    rel_relax_steps: int = 0,
    sc_restraint_k: float = 5000.0,
    cg_bead_masses: np.ndarray | None = None,
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

    phi_values = np.linspace(0.0, 2.0 * np.pi, azimuthal_count, endpoint=False)
    sidechain_bead_frame_angles = _bead_frame_angles(sidechain_bead_frame_count)
    cg_bead_frame_angles = _bead_frame_angles(cg_bead_frame_count)

    ref_nm = _canonicalize_lipid_reference_to_z(ref_bead_positions_nm)
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
          f"direct rotated geometry, modes={n_modes}")

    # Energy grid: (n_radial, n_angle_sc, n_angle_cg), matching runtime source order.
    energy_grid = np.zeros((n_radial, n_angle, n_angle), dtype=np.float64)

    for ir, r_nm in enumerate(r_values):
        for ia_sc in range(n_angle):
            for ia_cg in range(n_angle):
                energy_sum = 0.0
                weight_sum = 0.0
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
                            pair_energy_sum = 0.0
                            pair_sample_count = 0
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
                                    cg_positions = cg_com[None, :] + (R_cg @ ref_nm.T).T

                                    if rel_relax_steps > 0:
                                        # Lab-frame SC relaxation with position restraints.
                                        sc_body = _dopc_lab_to_body(
                                            framed_sc_positions, cb_vec,
                                            np.asarray(cb_anchor_nm, dtype=np.float64))
                                        sc_body, _ = _relax_sc_beads(
                                            init_body_positions=sc_body,
                                            ref_body_positions=sc_body.copy(),
                                            sc_bead_types=sc_bead_types,
                                            sc_bead_charges=sc_bead_charges,
                                            k_restraint=sc_restraint_k,
                                            pair_params=pair_params,
                                            direction_lab=cb_vec,
                                            cb_anchor_lab=np.asarray(cb_anchor_nm, dtype=np.float64),
                                            partner_lab_positions=cg_positions,
                                            partner_bead_types=cg_bead_types,
                                            partner_bead_charges=cg_bead_charges,
                                            dist_min_nm=0.10,
                                            rel_relax_steps=rel_relax_steps,
                                        )
                                        framed_sc_positions = _dopc_body_to_lab(
                                            sc_body, cb_vec,
                                            np.asarray(cb_anchor_nm, dtype=np.float64))

                                        # Lab-frame DOPC relaxation with soft COM/dir restraints.
                                        cg_body = _dopc_lab_to_body(cg_positions, dir_cg, cg_com)
                                        cg_body, _ = _relax_dopc_lab(
                                            init_body_positions=cg_body,
                                            ref_body_positions=ref_nm,
                                            bead_types=cg_bead_types,
                                            bead_charges=cg_bead_charges,
                                            bead_masses=cg_bead_masses,
                                            pair_params=pair_params,
                                            direction_lab=dir_cg,
                                            com_lab=cg_com,
                                            partner_lab_positions=framed_sc_positions,
                                            partner_bead_types=sc_bead_types,
                                            partner_bead_charges=sc_bead_charges,
                                            dist_min_nm=0.10,
                                            rel_relax_steps=rel_relax_steps,
                                            plane_constraint=True,
                                        )
                                        cg_positions = _dopc_body_to_lab(cg_body, dir_cg, cg_com)

                                    pair_energy, _, _ = _compute_pair_energy_and_gradient(
                                        cg_positions, framed_sc_positions, cg_bead_types, sc_bead_types,
                                        cg_bead_charges, sc_bead_charges, pair_params,
                                        dist_min_nm=0.10, soft_core_alpha=0.0,
                                    )
                                    pair_energy_sum += pair_energy
                                    pair_sample_count += 1
                            pair_energy = pair_energy_sum / max(pair_sample_count, 1)

                            energy_sum += w * pair_energy
                            weight_sum += w

                energy_grid[ir, ia_sc, ia_cg] = energy_sum / max(weight_sum, 1e-15)

    # WCA excluded area: add repulsive energy below contact distance.
    if excluded_area_contact_nm is not None:
        for ir, r_nm in enumerate(r_values):
            if float(r_nm) < float(excluded_area_contact_nm):
                sigma_nm = float(excluded_area_contact_nm) / (2.0 ** (1.0 / 6.0))
                sr = sigma_nm / max(float(r_nm), 1e-6)
                sr6 = sr ** 6
                wca = 4.0 * DEFAULT_PRODUCTION_KBT_KJ_MOL * (sr6 * sr6 - sr6) + DEFAULT_PRODUCTION_KBT_KJ_MOL
                energy_grid[ir, :, :] += wca

    # Isotropic background: subtract mean attractive contribution.
    radial_mean = energy_grid.mean(axis=(1, 2))
    attractive_background = np.minimum(radial_mean, 0.0)
    if np.any(attractive_background < 0.0):
        energy_grid = energy_grid - attractive_background[:, None, None]
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

    n_knot_angular = 15
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
    v0_knots[:n_unconstrained] = short_range_core_eup
    param_parts = [v0_knots]
    ang1_knots_all = []
    ang2_knots_all = []
    vm_knots_all = []
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        ang1_knots = _fit_angular_bspline(t_angular, ang1, n_knot_angular, smooth=0.01)
        ang2_knots = _fit_angular_bspline(t_angular, ang2, n_knot_angular, smooth=0.01)
        vm_knots = _fit_radial_bspline(t_radial_ang, vm, rad_knot_vector, smooth=0.01) * inv_conv
        ang1_knots = np.clip(ang1_knots, -1.0, 1.0)
        ang2_knots = np.clip(ang2_knots, -1.0, 1.0)
        vm_knots[:n_unconstrained] = 0.0
        ang1_knots_all.append(ang1_knots)
        ang2_knots_all.append(ang2_knots)
        vm_knots_all.append(vm_knots)
        param_parts.extend([ang1_knots, ang2_knots, vm_knots])

    excluded_area_rows = 0
    if excluded_area_contact_nm is not None and excluded_area_contact_nm > 0.0:
        excluded_area_rows = max(
            0,
            min(
                n_knot_radial,
                int(math.ceil(float(excluded_area_contact_nm) * LENGTH_CONVERSION_A_PER_NM / knot_spacing)) + 1,
            ),
        )
        if excluded_area_rows > 0:
            v0_knots[:excluded_area_rows] = np.maximum(v0_knots[:excluded_area_rows], 0.0)
            for vm_knots in vm_knots_all:
                vm_knots[:excluded_area_rows] = 0.0
            short_range_core_kj_mol = max(0.0, raw_short_range_core_kj_mol)

    # Attractive removal: clamp all negative control points to zero.
    # Negative B-spline controls produce unphysical many-body attraction
    # that is non-transferable from the two-body training.
    attractive_control_count = 0
    attractive_control_count += int(np.count_nonzero(v0_knots < 0.0))
    v0_knots = np.maximum(v0_knots, 0.0)
    for vm_knots in vm_knots_all:
        attractive_control_count += int(np.count_nonzero(vm_knots < 0.0))
        vm_knots[:] = np.maximum(vm_knots, 0.0)
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
        "excluded_area_source": (
            "wca_dopc_contact_kbt"
            if excluded_area_rows > 0
            else ""
        ),
        "excluded_area_nonnegative_rows": int(excluded_area_rows),
        "isotropic_background_source": (
            "attractive_radial_angular_mean_subtracted"
            if np.any(attractive_background < 0.0)
            else "none_full_resolved"
        ),
        "attractive_control_source": (
            "nontransferable_attraction_removed"
            if attractive_control_count > 0
            else "retained_full_resolved"
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
        "azimuthal_count": int(azimuthal_count),
        "sidechain_bead_frame_count": int(len(sidechain_bead_frame_angles)),
        "cg_bead_frame_count": int(len(cg_bead_frame_angles)),
        "rel_relax_steps": int(rel_relax_steps),
    }


def _fit_cg_lipid_sc_quadspline_from_dict(task: dict) -> tuple[int, str, dict]:
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
    ref_nm = task["ref_nm"]
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
                framed_ref = (rot @ ref_nm.T).T
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
    ref_nm = _canonicalize_lipid_reference_to_z(np.asarray(ref_bead_positions_nm, dtype=np.float64))
    n_beads = ref_nm.shape[0]
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
    plane_constraint: bool = False,
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
    if ref_nm.shape != (14, 3):
        raise ValueError(f"ref_bead_positions_nm must be (14, 3), got {ref_nm.shape}")

    print("\n=== CG Lipid Table Building ===")
    bead_mass_values = None
    if bead_masses_g_mol is not None:
        bead_mass_values = [bead_masses_g_mol[bt] for bt in bead_types]
    derived_params = derive_dopc_cg_params(
        ref_bead_positions_nm=ref_nm,
        bead_types=bead_types,
        pair_params=pair_params,
        bead_masses_g_mol=bead_mass_values,
        bonds=_CURRENT_CG_BONDS,
        energy_conversion_kj_per_eup=ENERGY_CONVERSION_KJ_PER_EUP,
        length_conversion_ang_per_nm=LENGTH_CONVERSION_A_PER_NM,
    )
    contact_nm = float(derived_params["contact_nm"])
    sc_fit_r_max_nm = min(float(r_max_nm), contact_nm)
    print(
        "  DOPC-derived CGL params: "
        f"contact={derived_params['contact_ang']:.3f} A, "
        f"orientation_length={derived_params['orientation_length_ang']:.3f} A, "
        f"orientation_mass={derived_params['orientation_mass_g_mol']:.3f} g/mol, "
        f"orientation_bond_fc={derived_params['orientation_bond_fc_eup_a2']:.3f} E_up/A^2"
    )

    # Resolution control: "coarse", "medium", or "fine" (default).  The CGL-CGL
    # excluded-area projection is kept at the full grid because reduced angular
    # sampling admits unphysical same-leaflet overlap in bilayer-only validation.
    _res = os.environ.get("UPSIDE_CG_LIPID_RESOLUTION", "fine").strip().lower()
    if _res == "coarse":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 2, 8, 5, 2
    elif _res == "medium":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 2, 12, 7, 2
    else:
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 2, 16, 9, 4
    cg_bead_frame_count = _bead_frame_count("CGL", 1)
    sc_bead_frame_count = _bead_frame_count("SC", 1)
    print(f"  CG lipid resolution: {_res} "
          f"(CG: {_cg_r}r x {_cg_ct}^2 theta x {_cg_az}^2 phi x {cg_bead_frame_count}^2 bead-frame, "
          f"SC: {_sc_r}r x {_sc_ct}^2 theta x {_sc_az}^2 phi "
          f"x {sc_bead_frame_count} SC bead-frame x {cg_bead_frame_count} CGL bead-frame)")

    fit_relax_steps = _positive_int_env("UPSIDE_MARTINI_FIT_RELAX_STEPS", 50)
    sc_restraint_k = float(os.environ.get("UPSIDE_MARTINI_SC_RESTRAINT_K", "5000.0"))
    if fit_relax_steps > 0:
        print(f"  Hidden-bead relaxation: {fit_relax_steps} steps, "
              f"SC restraint k={sc_restraint_k:.1f} kJ/mol/nm^2")

    # CG-CG directional spline from relaxed DOPC bead geometries when
    # fit_relax_steps > 0, otherwise direct rigid-geometry projection.
    result_cg = _fit_cg_lipid_quadspline(
        ref_bead_positions_nm=ref_nm,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        r_min_nm=0.50,
        r_max_nm=1.68,
        r_count=min(r_count, _cg_r),
        cos_theta_count=min(cos_theta_count, _cg_ct),
        azimuthal_count=_cg_az,
        bead_frame_count=cg_bead_frame_count,
        dist_min_nm=0.25,
        knot_spacing_ang=1.4,
        excluded_area_contact_nm=contact_nm,
        n_modes=4,
        n_knot_radial=14,
        n_knot_angular=15,
        cg_smooth=0.01,
        rel_relax_steps=fit_relax_steps,
        bead_masses=bead_mass_values,
        relax_soft_core_alpha=0.0,
        plane_constraint=plane_constraint,
    )
    print(
        "  CG-CG: full tensor table "
        f"{result_cg['n_radial']}r x {result_cg['n_angular']}^2 angular, "
        f"max|E| = {float(np.max(np.abs(result_cg['energy_grid_raw']))):.3f} kJ/mol"
    )

    # CG-SC quadspline
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
        martinize_path, forcefield_name
    )

    residues = [r for r in active_residue_names
                if r in residue_map and residue_map[r] and r in orientation_map]
    sc_n_modes = int(os.environ.get("UPSIDE_CG_LIPID_SC_N_MODES", "4"))
    sc_n_radial = 14
    sc_n_angular = 15
    sc_knot_spacing_ang = 1.4
    sc_taper_width_ang = sc_knot_spacing_ang
    sc_cutoff_ang = float(sc_fit_r_max_nm * 10.0 + sc_taper_width_ang)
    sc_n_param = sc_n_radial + sc_n_modes * (2 * sc_n_angular + sc_n_radial)
    sc_rms_error = np.zeros(len(residues), dtype=np.float32)
    sc_short_range_core = np.zeros(len(residues), dtype=np.float32)
    sc_short_range_core_rows = np.zeros(len(residues), dtype=np.int32)

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

        sc_fit_tasks = []
        first_sc_result = None
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
                    "ref_bead_positions_nm": ref_nm,
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
                    "r_count": min(r_count, _sc_r),
                    "cos_theta_count": min(cos_theta_count, _sc_ct),
                    "azimuthal_count": min(azimuthal_count, _sc_az),
                    "sidechain_bead_frame_count": sc_bead_frame_count,
                    "cg_bead_frame_count": cg_bead_frame_count,
                    "n_knot_radial": sc_n_radial,
                    "knot_spacing_ang": sc_knot_spacing_ang,
                    "excluded_area_contact_nm": contact_nm,
                    "rel_relax_steps": int(fit_relax_steps),
                    "sc_restraint_k": float(sc_restraint_k),
                    "cg_bead_masses": (
                        np.asarray(bead_mass_values, dtype=np.float64)
                        if bead_mass_values is not None
                        else None
                    ),
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
            sc_cutoff_ang = min(sc_cutoff_ang, float(result_sc["cutoff_ang"]))
            sc_residue_names.append(residue)
            print(f"  CG-SC({residue}): RMS error = {result_sc['rms_error']:.4f} kJ/mol, "
                  f"modes = {result_sc['n_modes']}")

    # Store in HDF5
    cg_grp = h5.create_group("cg_lipid_table")
    cg_grp.attrs["bead_charge_source"] = "dry_martini_v2.1_lipids.itp:DOPC_atoms"
    cg_grp.attrs["lipid_net_charge"] = np.float32(sum(bead_charges))
    cg_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_grp.create_dataset("bead_charges", data=np.asarray(bead_charges, dtype=np.float32))
    _write_cg_derived_attrs(cg_grp, derived_params)

    # CG-CG pair
    cg_pair_grp = cg_grp.create_group("cg_lipid_pair")
    pair_param = result_cg["interaction_param"].astype(np.float32)
    cg_pair_grp.create_dataset(
        "interaction_param",
        data=pair_param.reshape(1, 1, pair_param.size),
    )
    cg_pair_grp.attrs["n_cg_types"] = 1
    cg_pair_grp.attrs["rms_error_kj_mol"] = np.float32(result_cg["rms_error"])
    cg_pair_grp.attrs["schema"] = result_cg["schema"]
    cg_pair_grp.attrs["radial_mode"] = "full_tensor"
    cg_pair_grp.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
    cg_pair_grp.attrs["fit_relax_steps"] = np.int32(result_cg.get("rel_relax_steps", 0))
    cg_pair_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_pair_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_pair_grp.attrs["fit_r_min_nm"] = np.float32(result_cg["r_values_nm"][0])
    cg_pair_grp.attrs["fit_r_max_nm"] = np.float32(1.68)
    cg_pair_grp.attrs["n_modes"] = np.int32(result_cg["n_modes"])
    cg_pair_grp.attrs["n_radial"] = np.int32(result_cg["n_radial"])
    cg_pair_grp.attrs["n_angular"] = np.int32(result_cg["n_angular"])
    cg_pair_grp.attrs["azimuthal_count"] = np.int32(result_cg["azimuthal_count"])
    cg_pair_grp.attrs["cgl_bead_frame_count"] = np.int32(result_cg["bead_frame_count"])
    cg_pair_grp.attrs["orientation_sampling"] = "both_cgl_direction_vectors"
    cg_pair_grp.attrs["knot_spacing_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["cutoff_ang"] = np.float32(result_cg["cutoff_ang"])
    cg_pair_grp.attrs["taper_width_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["azimuthal_average"] = result_cg["azimuthal_average"]
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
    _write_cg_derived_attrs(cg_pair_grp, derived_params)
    cg_pair_grp.create_dataset("energy_grid_raw_kj_mol", data=result_cg["energy_grid_raw"].astype(np.float32))
    cg_pair_grp.create_dataset(
        "attractive_radial_background_kj_mol",
        data=result_cg["attractive_radial_background_kj_mol"].astype(np.float32),
    )

    # CG-SC
    cg_sc_grp = cg_grp.create_group("cg_lipid_sc")
    cg_sc_grp.create_dataset(
        "interaction_param",
        data=interaction_param_sc.astype(np.float32),
    )
    cg_sc_grp.create_dataset(
        "restype_order",
        data=np.asarray([np.bytes_(x) for x in sc_residue_names], dtype="S4"),
    )
    cg_sc_grp.create_dataset("rms_error_kj_mol", data=sc_rms_error[:n_sc_types])
    cg_sc_grp.attrs["n_sc_types"] = n_sc_types
    cg_sc_grp.attrs["n_cg_types"] = 1
    cg_sc_grp.attrs["schema"] = "cg_lipid_sc_quadspline_v2"
    cg_sc_grp.attrs["radial_mode"] = "full_multimode"
    cg_sc_grp.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
    cg_sc_grp.attrs["fit_relax_steps"] = np.int32(fit_relax_steps if n_sc_types else 0)
    cg_sc_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    cg_sc_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
    cg_sc_grp.attrs["fit_r_min_nm"] = np.float32(r_min_nm)
    cg_sc_grp.attrs["fit_r_max_nm"] = np.float32(sc_fit_r_max_nm)
    cg_sc_grp.attrs["n_modes"] = np.int32(sc_n_modes)
    cg_sc_grp.attrs["n_radial"] = np.int32(sc_n_radial)
    cg_sc_grp.attrs["n_angular"] = np.int32(sc_n_angular)
    cg_sc_grp.attrs["azimuthal_count"] = np.int32(min(azimuthal_count, _sc_az) if n_sc_types else 0)
    cg_sc_grp.attrs["sidechain_bead_frame_count"] = np.int32(sc_bead_frame_count if n_sc_types else 0)
    cg_sc_grp.attrs["cgl_bead_frame_count"] = np.int32(cg_bead_frame_count if n_sc_types else 0)
    cg_sc_grp.attrs["orientation_sampling"] = "sidechain_and_cgl_direction_vectors"
    cg_sc_grp.attrs["knot_spacing_ang"] = np.float32(sc_knot_spacing_ang)
    cg_sc_grp.attrs["cutoff_ang"] = np.float32(sc_cutoff_ang)
    cg_sc_grp.attrs["taper_width_ang"] = np.float32(sc_taper_width_ang)
    cg_sc_grp.attrs["azimuthal_average"] = "energy_expectation"
    cg_sc_grp.attrs["short_range_core_source"] = "max_first_sampled_dry_martini_energy_expectation"
    cg_sc_grp.attrs["excluded_area_source"] = "wca_dopc_contact_kbt" if n_sc_types else ""
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
        int(math.ceil(contact_nm * LENGTH_CONVERSION_A_PER_NM / sc_knot_spacing_ang)) + 1
        if n_sc_types
        else 0
    )
    cg_sc_grp.create_dataset("short_range_core_kj_mol", data=sc_short_range_core[:n_sc_types])
    cg_sc_grp.create_dataset("short_range_core_rows", data=sc_short_range_core_rows[:n_sc_types])
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
        ref_bead_positions_nm=ref_nm,
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
    eff_grp.attrs["orientation_sampling"] = "fibonacci_cgl_direction_vectors"
    _write_cg_derived_attrs(eff_grp, derived_params)

    # Build directional B-spline tables for CGL against all non-CGL target types.
    # After this, MartiniPotential CGL-X can be omitted for all X.
    _build_cgl_target_table(
        h5, cg_grp, effective_lj,
        ref_bead_positions_nm=ref_nm,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        derived_params=derived_params,
        rel_relax_steps=fit_relax_steps,
        cg_bead_masses=np.asarray(bead_mass_values, dtype=np.float64) if bead_mass_values is not None else None,
    )
    print()


def _run_cgl_target_type_task(task: dict) -> tuple[int, str, np.ndarray, int]:
    ti = int(task["ti"])
    tgt_type = str(task["target_type"])
    bead_types = task["bead_types"]
    bead_charges = task["bead_charges"]
    pair_params = task["pair_params"]
    ref_nm = task["ref_nm"]
    r_sample_nm = task["r_sample_nm"]
    cos_theta_grid = task["cos_theta_grid"]
    orientation_dirs = task["orientation_dirs"]
    bead_frame_angles = task["bead_frame_angles"]
    n_knot_radial = int(task["n_knot_radial"])
    n_knot_angular = int(task["n_knot_angular"])
    knot_spacing_ang = float(task["knot_spacing_ang"])
    energy_conv = float(task["energy_conv"])
    r_min_ang = float(task["r_min_ang"])
    contact_ang = float(task["contact_ang"])

    energy_grid = np.zeros((len(r_sample_nm), cos_theta_grid.size), dtype=np.float64)
    target_charge = infer_charge_from_atomtype(tgt_type)

    z_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    target_pos = np.zeros((1, 3), dtype=np.float64)
    rel_relax_steps = int(task.get("rel_relax_steps", 0))
    cg_bead_masses = task.get("cg_bead_masses")
    for ir, r_nm in enumerate(r_sample_nm):
        target_pos[0, :] = (float(r_nm), 0.0, 0.0)
        for ia in range(cos_theta_grid.size):
            e_sum = 0.0
            sample_count = 0
            for direction in orientation_dirs[ia]:
                rot_base = _rotation_to_align_z_np(direction)
                for frame_angle in bead_frame_angles:
                    if rel_relax_steps > 0:
                        init_body = (_rotation_about_axis_np(z_axis, float(frame_angle)) @ ref_nm.T).T
                        _, energy = _relax_dopc_lab(
                            init_body_positions=init_body,
                            ref_body_positions=ref_nm,
                            bead_types=bead_types,
                            bead_charges=bead_charges,
                            bead_masses=cg_bead_masses,
                            pair_params=pair_params,
                            direction_lab=np.asarray(direction, dtype=np.float64),
                            com_lab=np.zeros(3, dtype=np.float64),
                            partner_lab_positions=target_pos,
                            partner_bead_types=[tgt_type],
                            partner_bead_charges=[target_charge],
                            dist_min_nm=0.10,
                            rel_relax_steps=rel_relax_steps,
                            plane_constraint=True,
                        )
                        e_sum += energy
                    else:
                        rot = rot_base @ _rotation_about_axis_np(z_axis, float(frame_angle))
                        cg_positions = (rot @ ref_nm.T).T
                        e, _, _ = _compute_pair_energy_and_gradient(
                            cg_positions,
                            target_pos,
                            bead_types,
                            [tgt_type],
                            bead_charges,
                            [target_charge],
                            pair_params,
                            dist_min_nm=0.10,
                            soft_core_alpha=0.0,
                        )
                        e_sum += e
                    sample_count += 1
            energy_grid[ir, ia] = e_sum / max(sample_count, 1)

    # Isotropic background: subtract mean attractive contribution.
    radial_mean = energy_grid.mean(axis=1)
    attractive_bg = np.minimum(radial_mean, 0.0)
    if np.any(attractive_bg < 0.0):
        energy_grid = energy_grid - attractive_bg[:, None]

    control = _fit_radial_angular_tensor_bspline(
        r_sample_nm,
        cos_theta_grid,
        energy_grid,
        n_knot_radial=n_knot_radial,
        n_knot_angular=n_knot_angular,
        knot_spacing_ang=knot_spacing_ang,
        energy_conversion=energy_conv,
        smooth=0.01,
    )
    n_core = max(1, min(n_knot_radial - 1, int(math.ceil(r_min_ang / knot_spacing_ang))))
    boundary_row = min(n_core, n_knot_radial - 1)
    control[:n_core, :] = control[boundary_row:boundary_row + 1, :]
    nonnegative_rows = 0
    if contact_ang > 0.0:
        nonnegative_rows = max(
            0,
            min(n_knot_radial, int(math.ceil(contact_ang / knot_spacing_ang)) + 1),
        )
        if nonnegative_rows:
            control[:nonnegative_rows, :] = np.maximum(control[:nonnegative_rows, :], 0.0)
    # Attractive removal: clamp all negative B-spline controls to zero.
    control = np.maximum(control, 0.0)
    return ti, tgt_type, control.reshape(-1), nonnegative_rows


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
    n_knot_angular: int = 15,
    knot_spacing_ang: float = 1.4,
    rel_relax_steps: int = 0,
    cg_bead_masses: np.ndarray | None = None,
) -> None:
    """Build directional tensor B-spline tables for CGL-point targets."""
    target_types = sorted(t for t in effective_lj if t != "CGL")
    if not target_types:
        print("  cg_lipid_target: no target types, skipping")
        return

    n_types = len(target_types)
    n_param = n_knot_radial * n_knot_angular
    interaction_param = np.zeros((1, n_types, n_param), dtype=np.float64)

    # Sample densely for an accurate B-spline fit.  The interaction samples
    # explicit DOPC bead-vs-target energies and only uses a bead-scale distance
    # floor to avoid evaluating below the resolved dry-MARTINI core.
    r_min_ang = float(knot_spacing_ang) + 0.1
    r_max_ang = float((n_knot_radial - 2) * knot_spacing_ang)
    n_sample = 80
    r_sample_ang = np.linspace(r_min_ang, r_max_ang, n_sample)
    r_sample_nm = r_sample_ang / length_conv
    cos_theta_grid = np.asarray(_linspace(-1.0, 1.0, 9), dtype=np.float64)
    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64) if ref_bead_positions_nm is not None else None
    explicit_source = (
        ref_nm is not None
        and ref_nm.shape == (14, 3)
        and bead_types is not None
        and bead_charges is not None
        and pair_params is not None
    )
    if not explicit_source:
        raise RuntimeError("cg_lipid_target requires explicit DOPC bead geometry and dry-MARTINI parameters")
    ref_nm = _canonicalize_lipid_reference_to_z(ref_nm)
    target_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    phi_values = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)
    bead_frame_angles = _bead_frame_angles(_bead_frame_count("CGL", 1))
    orientation_dirs = _directions_with_dot_np(target_axis, cos_theta_grid, phi_values)

    contact_ang = float(derived_params.get("contact_ang", 0.0)) if derived_params else 0.0
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
            "contact_ang": contact_ang,
            "rel_relax_steps": int(rel_relax_steps),
            "cg_bead_masses": cg_bead_masses,
        }
        for ti, tgt_type in enumerate(target_types)
    ]
    nonnegative_rows = 0
    for ti, _tgt_type, control_flat, rows in _parallel_map_ordered(
        "CGL-particle target table", _run_cgl_target_type_task, target_tasks
    ):
        interaction_param[0, ti, :] = control_flat
        nonnegative_rows = max(nonnegative_rows, int(rows))

    target_grp = cg_grp.create_group("cg_lipid_target")
    target_grp.create_dataset(
        "interaction_param",
        data=interaction_param.astype(np.float32),
    )
    target_grp.create_dataset(
        "target_order",
        data=np.asarray([np.bytes_(x) for x in target_types], dtype="S8"),
    )
    target_grp.attrs["n_target_types"] = np.int32(n_types)
    target_grp.attrs["n_cg_types"] = np.int32(1)
    target_grp.attrs["schema"] = "cg_lipid_target_v1"
    target_grp.attrs["n_modes"] = np.int32(0)
    target_grp.attrs["n_radial"] = np.int32(n_knot_radial)
    target_grp.attrs["n_angular"] = np.int32(n_knot_angular)
    target_grp.attrs["azimuthal_count"] = np.int32(len(phi_values))
    target_grp.attrs["cgl_bead_frame_count"] = np.int32(len(bead_frame_angles))
    target_grp.attrs["orientation_sampling"] = "cgl_direction_vector"
    target_grp.attrs["knot_spacing_ang"] = np.float32(knot_spacing_ang)
    cutoff_ang = float((n_knot_radial - 2) * knot_spacing_ang)
    target_grp.attrs["cutoff_ang"] = np.float32(cutoff_ang)
    target_grp.attrs["taper_width_ang"] = np.float32(knot_spacing_ang)
    target_grp.attrs["azimuthal_average"] = "energy_expectation"
    target_grp.attrs["unresolved_core_source"] = "first_resolved_dry_martini_energy_expectation"
    target_grp.attrs["excluded_area_source"] = "dopc_contact_nonnegative_controls"
    target_grp.attrs["excluded_area_nonnegative_rows"] = np.int32(nonnegative_rows)
    target_grp.attrs["source"] = "explicit_dopc_directional"
    target_grp.attrs["fit_relax_steps"] = np.int32(rel_relax_steps)
    target_grp.attrs["isotropic_background_source"] = "attractive_radial_mean_subtracted"
    target_grp.attrs["attractive_control_source"] = "nontransferable_attraction_removed"
    target_grp.attrs["relaxation"] = "lab_frame_soft_restraints" if rel_relax_steps > 0 else "rigid_rotated_geometry"
    target_grp.attrs["angle_convention"] = "ang=n_cgl_dot_n12"
    target_grp.attrs["bead_nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    target_grp.attrs["bead_nonbonded_cutoff_source"] = "generic_martini_potential_cutoff"
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

    sc_fit_relax_steps = _positive_int_env("UPSIDE_MARTINI_FIT_RELAX_STEPS", 50)
    sc_restraint_k = float(os.environ.get("UPSIDE_MARTINI_SC_RESTRAINT_K", "5000.0"))

    with h5py.File(output_path, "w") as h5:
        h5.attrs["schema"] = "martini_combined_v1"
        _build_particles_group(h5, atomtypes, pair_params, active_atom_types)
        _build_sc_table_group(
            h5,
            residue_map,
            pair_params,
            sidechain_lib_path,
            active_residue_names=active_residue_names,
            active_target_types=sorted(active_atom_types),
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            fit_relax_steps=sc_fit_relax_steps,
            sc_restraint_k=sc_restraint_k,
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

    sc_fit_relax_steps = _positive_int_env("UPSIDE_MARTINI_FIT_RELAX_STEPS", 50)
    sc_restraint_k = float(os.environ.get("UPSIDE_MARTINI_SC_RESTRAINT_K", "5000.0"))

    def _writer(h5: h5py.File) -> None:
        h5.attrs["schema"] = SCHEMA_SC
        _build_sc_table_group(
            h5, residue_map, pair_params, sidechain_lib_path,
            active_residue_names=list(CANONICAL_RESIDUES),
            active_target_types=atomtypes,
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            fit_relax_steps=sc_fit_relax_steps,
            sc_restraint_k=sc_restraint_k,
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

    with open(dopc_pdb_path) as f:
        dopc_atoms = []
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:21].strip().upper()
            if resname not in ("DOPC", "DOP"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dopc_atoms.append([x, y, z])

    n_per_lipid = 14
    first = dopc_atoms[:n_per_lipid]
    com = np.mean(first, axis=0)
    ref_bead_positions_nm = np.array(
        [[a[0] - com[0], a[1] - com[1], a[2] - com[2]] for a in first]
    ) * 0.1

    def _writer(h5: h5py.File) -> None:
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
        )

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path}")


def build_interlipid_h5(output_path: Path) -> None:
    """Generate interlipid.h5 as an empty placeholder."""
    output_path = Path(output_path).expanduser().resolve()
    def _writer(h5: h5py.File) -> None:
        g = h5.create_group("cross_lipid")
        g.attrs["schema"] = "cross_lipid_v1"
        g.attrs["n_lipid_types"] = 1
        g.create_dataset("interaction_param", data=np.zeros((0,), dtype=np.float64))

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path} (empty cross-lipid placeholder)")
