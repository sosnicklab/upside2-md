#!/usr/bin/env python3
"""Build a per-run martini.h5 with combined LJ+Coulomb spline tables.

Generates particle-particle combined energy grids and sidechain-environment
orientation-aware combined energy tables, trimmed to only the types and
residues actually present in the current system.
"""

from __future__ import annotations

import importlib.util
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

COULOMB_K_DRY_KJ_NM = 138.935458 / 15.0
ENERGY_CONVERSION_KJ_PER_EUP = 2.914952774272
LENGTH_CONVERSION_A_PER_NM = 10.0
ANGSTROM_TO_NM = 0.1

PARTICLES_GRID_N = 1000
PARTICLES_R_MIN_A = 0.0
PARTICLES_R_MAX_A = 12.0

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)

POSITIVE_TYPES = {"Qda", "Qd", "SQda", "SQd"}
NEGATIVE_TYPES = {"Qa", "SQa"}

CANONICAL_CB_POSITION_ANG = (0.0, 0.94375626, 1.2068012)
_cb_norm = math.sqrt(sum(x * x for x in CANONICAL_CB_POSITION_ANG))
CANONICAL_CB_VECTOR_UNIT = tuple(x / _cb_norm for x in CANONICAL_CB_POSITION_ANG)

SCHEMA_PARTICLES = "martini_particles_combined_v1"
SCHEMA_SC = "martini_sc_combined_v1"


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _infer_type_charge(bead_type: str) -> float:
    bt = bead_type.strip()
    if bt in POSITIVE_TYPES:
        return 1.0
    if bt in NEGATIVE_TYPES:
        return -1.0
    return 0.0


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
        rms_error = float(np.sqrt(np.mean((sampled - (radial[None, :] + angular_profile[:, None] * angular_radial[None, :])) ** 2)))

    return (
        [float(x) for x in radial],
        [float(x) for x in angular_profile],
        [float(x) for x in angular_radial],
        rms_error,
    )


# ---------------------------------------------------------------------------
# Forcefield parsing
# ---------------------------------------------------------------------------

def _parse_dry_forcefield(
    ff_path: Path,
) -> Tuple[List[str], Dict[Tuple[str, str], Dict[str, float]]]:
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
            type_i, type_j = parts[0], parts[1]
            if len(parts) == 4:
                macro = parts[3]
                if macro not in macros:
                    raise RuntimeError(f"Unknown dry-MARTINI macro '{macro}' in {ff_path}")
                sigma_nm, epsilon_kj = macros[macro]
            else:
                sigma_nm = float(parts[3])
                epsilon_kj = float(parts[4])
            payload = {"sigma_nm": sigma_nm, "epsilon_kj_mol": epsilon_kj}
            pair_params[(type_i, type_j)] = payload
            pair_params[(type_j, type_i)] = payload
    return atomtypes, pair_params


def _load_martini_forcefield(
    martinize_path: Path, forcefield_name: str
) -> Dict[str, List[str]]:
    spec = importlib.util.spec_from_file_location(
        "sc_training_martinize_runtime", martinize_path
    )
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


# ---------------------------------------------------------------------------
# Particle-particle combined grids
# ---------------------------------------------------------------------------

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

    charges = np.asarray([_infer_type_charge(t) for t in atomtypes], dtype=np.float32)
    type_to_charge = {t: _infer_type_charge(t) for t in atomtypes}

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


# ---------------------------------------------------------------------------
# Sidechain-environment combined tables
# ---------------------------------------------------------------------------

def _run_sc_task(
    residue: str,
    sidechain_bead_types: List[str],
    sidechain_bead_charges: List[float],
    target_type: str,
    target_charge: float,
    rotamer_centers_nm: List[List[float]],
    rotamer_weights: List[float],
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    r_values: List[float],
    direction_vectors: List[List[float]],
    cos_theta_grid: List[float],
    cb_anchor_nm: List[float],
    cb_vector_unit: List[float],
) -> Dict[str, Any]:
    n_rotamer = len(rotamer_centers_nm)
    n_angle = len(cos_theta_grid)
    n_radial = len(r_values)

    angular_energy = [[0.0 for _ in r_values] for _ in cos_theta_grid]
    rotamer_angular_energy = [
        [[0.0 for _ in r_values] for _ in cos_theta_grid]
        for _ in range(n_rotamer)
    ]

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

            for irot, (center_nm, rot_weight) in enumerate(
                zip(rotamer_centers_nm, rotamer_weights)
            ):
                dx = target_pos_nm[0] - center_nm[0]
                dy = target_pos_nm[1] - center_nm[1]
                dz = target_pos_nm[2] - center_nm[2]
                dist_nm = max(1.0e-6, math.sqrt(dx * dx + dy * dy + dz * dz))

                rot_energy = 0.0
                for bead_type, bead_charge in zip(sidechain_bead_types, sidechain_bead_charges):
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

    rotamer_radial_energy = []
    rotamer_angular_profile = []
    rotamer_angular_radial_energy = []
    rotamer_rms_error = []
    for irot in range(n_rotamer):
        rr, rp, ra, rrm = _factorize_one_sided_orientation(
            rotamer_angular_energy[irot], cos_theta_grid
        )
        rotamer_radial_energy.append(rr)
        rotamer_angular_profile.append(rp)
        rotamer_angular_radial_energy.append(ra)
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
    }


def _build_sc_table_group(
    h5: h5py.File,
    residue_map: Dict[str, List[str]],
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    sidechain_lib_path: Path,
    active_residue_names: List[str],
    active_target_types: List[str],
    r_min_nm: float = 0.25,
    r_max_nm: float = 1.20,
    r_count: int = 96,
    direction_count: int = 24,
    cos_theta_count: int = 13,
) -> None:
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    r_values = _linspace(r_min_nm, r_max_nm, r_count)
    direction_vectors = _fibonacci_sphere(direction_count)
    cos_theta_grid = _linspace(-1.0, 1.0, cos_theta_count)
    cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
    cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)

    # Determine residue-task and target-task lists
    residue_tasks = []
    for residue in active_residue_names:
        if residue not in residue_map:
            continue
        bead_types = residue_map[residue]
        if not bead_types:
            continue
        orientation = orientation_map.get(residue)
        if not orientation or not orientation["center_nm"]:
            raise RuntimeError(
                f"Missing orientation geometry for residue {residue} in {sidechain_lib_path}"
            )
        bead_charges = [_infer_type_charge(bt) for bt in bead_types]
        residue_tasks.append((residue, bead_types, bead_charges, orientation))

    if not residue_tasks:
        print("  sc_table: no active residues with sidechains, skipping")
        return

    # Collect results
    by_residue: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for residue, bead_types, bead_charges, orientation in residue_tasks:
        by_residue[residue] = {}
        for target_type in active_target_types:
            missing = [bt for bt in bead_types if (bt, target_type) not in pair_params]
            if missing:
                raise RuntimeError(
                    f"Missing nonbond param for residue {residue} target {target_type}"
                )
            result = _run_sc_task(
                residue=residue,
                sidechain_bead_types=bead_types,
                sidechain_bead_charges=bead_charges,
                target_type=target_type,
                target_charge=_infer_type_charge(target_type),
                rotamer_centers_nm=orientation["center_nm"],
                rotamer_weights=orientation["weight"],
                pair_params=pair_params,
                r_values=r_values,
                direction_vectors=direction_vectors,
                cos_theta_grid=cos_theta_grid,
                cb_anchor_nm=cb_anchor_nm,
                cb_vector_unit=cb_vector_unit,
            )
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


# ---------------------------------------------------------------------------
# Uniform deBoor B-spline utilities (matching C++ spline.h)
# ---------------------------------------------------------------------------
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
    memory (Ang1→Ang2). The 4th coefficient weight is zero for integer t, so
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


# ---------------------------------------------------------------------------
# CG lipid energy sampling
# ---------------------------------------------------------------------------

# DOPC bead types in ITP order (must match _DOPC_ATOM_NAMES in martini_prepare_system_lib.py)
_DOPC_BEAD_TYPES = [
    "Q0", "Qa", "Na", "Na",
    "C1", "C1", "C3", "C1", "C1",
    "C1", "C1", "C3", "C1", "C1",
]

# DOPC bonded topology from dry_martini_v2.1_lipids.itp (0-based indices)
_DOPC_BONDS: list = [
    (0, 1, 0.450, 1250.0),   # Q0-Qa     mb_np
    (1, 2, 0.450, 1250.0),   # Qa-Na     mb_pg1
    (2, 3, 0.370, 1250.0),   # Na-Na     mb_gg
    (2, 4, 0.480, 1250.0),   # Na-C1     mb_cc
    (4, 5, 0.480, 1250.0),   # C1-C1     mb_cc
    (5, 6, 0.480, 1250.0),   # C1-C3     mb_cc
    (6, 7, 0.480, 1250.0),   # C3-C1     mb_cc
    (7, 8, 0.480, 1250.0),   # C1-C1     mb_cc
    (3, 9, 0.480, 1250.0),   # Na-C1     mb_cc
    (9, 10, 0.480, 1250.0),  # C1-C1     mb_cc
    (10, 11, 0.480, 1250.0), # C1-C3     mb_cc
    (11, 12, 0.480, 1250.0), # C3-C1     mb_cc
    (12, 13, 0.480, 1250.0), # C1-C1     mb_cc
]

# Cosine-based angle potentials: V(θ) = 0.5 * k * (cos(θ) - cos(θ0))²
_DOPC_ANGLES: list = [
    (1, 2, 3, 120.0, 25.0),    # Qa-Na-Na      ma_pgg
    (1, 2, 4, 180.0, 25.0),    # Qa-Na-C1      ma_pgc
    (2, 4, 5, 180.0, 35.0),    # Na-C1-C1      ma_gcc
    (4, 5, 6, 180.0, 35.0),    # C1-C1-C3      ma_ccc
    (5, 6, 7, 120.0, 45.0),    # C1-C3-C1      ma_cdc
    (6, 7, 8, 180.0, 35.0),    # C3-C1-C1      ma_ccc
    (3, 9, 10, 180.0, 35.0),   # Na-C1-C1      ma_gcc
    (9, 10, 11, 180.0, 35.0),  # C1-C1-C3      ma_ccc
    (10, 11, 12, 120.0, 45.0), # C1-C3-C1      ma_cdc
    (11, 12, 13, 180.0, 35.0), # C3-C1-C1      ma_ccc
]


def _compute_lipid_bonded_energy(positions: np.ndarray) -> float:
    """Compute harmonic bond + cosine-based angle energy for ONE DOPC lipid.

    positions: (14, 3) float64 array of bead positions in nm.
    Returns energy in kJ/mol.
    """
    energy = 0.0
    # Bonds: V(r) = 0.5 * k * (r - r0)²
    for i, j, r0, k in _DOPC_BONDS:
        dr = positions[i] - positions[j]
        r = float(np.sqrt(np.dot(dr, dr)))
        energy += 0.5 * k * (r - r0) ** 2
    # Angles: V(θ) = 0.5 * k * (cos(θ) - cos(θ0))²
    for i, j, k_idx, theta0_deg, k_ang in _DOPC_ANGLES:
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


def _compute_lipid_bonded_energy_and_gradient(positions: np.ndarray):
    """Compute harmonic bond + cosine-based angle energy and analytical gradient.

    positions: (14, 3) float64 array of bead positions in nm.
    Returns (energy_kj_mol, grad) where grad is (14, 3).
    """
    n_beads = positions.shape[0]
    energy = 0.0
    grad = np.zeros_like(positions)

    # Bonds: V(r) = 0.5 * k * (r - r0)²
    for i, j, r0, k in _DOPC_BONDS:
        dr = positions[i] - positions[j]
        r = float(np.sqrt(np.dot(dr, dr)))
        if r < 1e-10:
            continue
        de_dr = float(k) * (r - r0)
        energy += 0.5 * float(k) * (r - r0) ** 2
        g = (de_dr / r) * dr
        grad[i] += g
        grad[j] -= g

    # Angles: V(θ) = 0.5 * k * (cos(θ) - cos(θ0))²
    for i, j, k_idx, theta0_deg, k_ang in _DOPC_ANGLES:
        r_ij = positions[i] - positions[j]
        r_kj = positions[k_idx] - positions[j]
        d_ij = float(np.sqrt(np.dot(r_ij, r_ij)))
        d_kj = float(np.sqrt(np.dot(r_kj, r_kj)))
        if d_ij < 1e-8 or d_kj < 1e-8:
            continue
        cos_theta = float(np.dot(r_ij, r_kj)) / (d_ij * d_kj)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        cos_theta0 = math.cos(math.radians(theta0_deg))
        dc = cos_theta - cos_theta0
        energy += 0.5 * float(k_ang) * dc * dc

        de_dc = float(k_ang) * dc
        r_ij_unit = r_ij / d_ij
        r_kj_unit = r_kj / d_kj
        dc_di = (r_kj_unit - cos_theta * r_ij_unit) / d_ij
        dc_dk = (r_ij_unit - cos_theta * r_kj_unit) / d_kj
        dc_dj = -(dc_di + dc_dk)
        grad[i] += de_dc * dc_di
        grad[j] += de_dc * dc_dj
        grad[k_idx] += de_dc * dc_dk

    return energy, grad


def _relax_lipid_pair(
    pos1: np.ndarray,
    pos2: np.ndarray,
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
    n_steps: int = 200,
) -> tuple:
    """Bond-constrained relaxation with real LJ (no distance floor during descent).

    Uses dist_min=0.01 during relaxation (purely numerical safety) so the
    energy surface retains genuine gradients everywhere, allowing the descent
    to push overlapping beads apart. Gradient clipping (max_disp=0.005 nm/step)
    keeps the enormous LJ forces bounded. No line search — fixed-step descent
    always moves toward lower energy.

    Final energy is evaluated with dist_min=0.10 to bound values for SVD.

    Returns (relaxed_pos1, relaxed_pos2, final_nonbonded_energy).
    """
    p1 = pos1.copy()
    p2 = pos2.copy()
    step_size = 0.0005

    for _ in range(n_steps):
        nb_energy, grad1_nb, grad2_nb = _compute_pair_energy_and_gradient(
            p1, p2, bead_types1, bead_types2,
            bead_charges1, bead_charges2, pair_params,
            dist_min_nm=0.01, soft_core_alpha=0.0,
        )
        bond1, grad1_b = _compute_lipid_bonded_energy_and_gradient(p1)
        bond2, grad2_b = _compute_lipid_bonded_energy_and_gradient(p2)

        grad1 = grad1_nb + grad1_b
        grad2 = grad2_nb + grad2_b

        # Gradient clipping: max 0.005 nm per bead per step
        disp1 = step_size * grad1
        disp2 = step_size * grad2
        d1_norms = np.sqrt(np.sum(disp1 ** 2, axis=1))
        d2_norms = np.sqrt(np.sum(disp2 ** 2, axis=1))
        max_grad_disp = max(float(np.max(d1_norms)), float(np.max(d2_norms)), 1e-15)
        if max_grad_disp > 0.005:
            scale = 0.005 / max_grad_disp
            grad1 *= scale
            grad2 *= scale

        p1 -= step_size * grad1
        p2 -= step_size * grad2

    # Evaluate final energy with distance floor for bounded SVD input
    final_nb, _, _ = _compute_pair_energy_and_gradient(
        p1, p2, bead_types1, bead_types2,
        bead_charges1, bead_charges2, pair_params,
        dist_min_nm=0.10, soft_core_alpha=0.0,
    )
    return p1, p2, final_nb


def _relax_lipid_against_points(
    cg_pos: np.ndarray,
    target_positions: np.ndarray,
    cg_bead_types: list,
    target_bead_types: list,
    cg_bead_charges: list,
    target_bead_charges: list,
    pair_params: dict,
    n_steps: int = 200,
) -> tuple:
    """Relax CG lipid beads against fixed target beads with real LJ.

    CG beads are bond-connected (can stretch/bend). Target beads are fixed
    (e.g., sidechain rotamer centers). Same fixed-step descent + gradient
    clipping as _relax_lipid_pair. Final energy uses dist_min=0.10 for SVD.

    Returns (relaxed_cg_pos, final_nb_energy).
    """
    p_cg = cg_pos.copy()
    step_size = 0.0005

    for _ in range(n_steps):
        nb_energy, grad_cg, grad_target = _compute_pair_energy_and_gradient(
            p_cg, target_positions, cg_bead_types, target_bead_types,
            cg_bead_charges, target_bead_charges, pair_params,
            dist_min_nm=0.01, soft_core_alpha=0.0,
        )
        bond_energy, grad_bond = _compute_lipid_bonded_energy_and_gradient(p_cg)

        grad_cg = grad_cg + grad_bond  # target grad is discarded (fixed)

        disp = step_size * grad_cg
        norms = np.sqrt(np.sum(disp ** 2, axis=1))
        max_grad_disp = max(float(np.max(norms)), 1e-15)
        if max_grad_disp > 0.005:
            grad_cg *= 0.005 / max_grad_disp

        p_cg -= step_size * grad_cg

    final_nb, _, _ = _compute_pair_energy_and_gradient(
        p_cg, target_positions, cg_bead_types, target_bead_types,
        cg_bead_charges, target_bead_charges, pair_params,
        dist_min_nm=0.10, soft_core_alpha=0.0,
    )
    return p_cg, final_nb


def _compute_cg_pair_energy(
    r_nm: float,
    dir1: np.ndarray,
    dir2: np.ndarray,
    ref_bead1: np.ndarray,
    ref_bead2: np.ndarray,
    bead_types1: list,
    bead_types2: list,
    bead_charges1: list,
    bead_charges2: list,
    pair_params: dict,
    relax_steps: int = 0,
    dist_min_nm: float = 0.20,
) -> float:
    """Compute LJ + Coulomb energy between two CG lipid particles at given geometry.

    Each CG lipid's beads are placed at COM_i + R(dir_i) * ref_bead_i[k].
    If relax_steps > 0, runs bond-constrained minimization with the real LJ
    before computing energy, allowing bonds to stretch and absorb LJ repulsion
    rather than diverging.
    """
    def _rotation_to_align_z(dir_vec):
        z_axis = dir_vec / np.linalg.norm(dir_vec)
        if abs(z_axis[0]) < 0.99:
            x_axis = np.cross([1.0, 0.0, 0.0], z_axis)
        else:
            x_axis = np.cross([0.0, 1.0, 0.0], z_axis)
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)
        return np.array([x_axis, y_axis, z_axis]).T

    R1 = _rotation_to_align_z(dir1)
    R2 = _rotation_to_align_z(dir2)

    com1 = np.zeros(3)
    com2 = np.array([r_nm, 0.0, 0.0])

    pos1 = com1[None, :] + (R1 @ ref_bead1.T).T
    pos2 = com2[None, :] + (R2 @ ref_bead2.T).T

    if relax_steps > 0:
        _, _, total_energy = _relax_lipid_pair(
            pos1, pos2, bead_types1, bead_types2,
            bead_charges1, bead_charges2, pair_params,
            n_steps=relax_steps,
        )
        return total_energy

    total_energy, _, _ = _compute_pair_energy_and_gradient(
        pos1, pos2, bead_types1, bead_types2,
        bead_charges1, bead_charges2, pair_params,
        dist_min_nm=dist_min_nm, soft_core_alpha=0.0,
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
    relax_steps: int = 0,
    dist_min_nm: float = 0.25,
    knot_spacing_ang: float = 1.4,
    residual_cap_kj_mol: float | None = 250.0,
    energy_cap_kj_mol: float | None = 500.0,
    n_modes: int = 4,
    n_knot_radial: int = 14,
) -> dict:
    """Fit full multi-mode B-spline parameters for CG lipid ↔ CG lipid interactions.

    Samples the energy landscape E(r, cosθ1, cosθ2) and decomposes via:
      E(r, cosθ1, cosθ2) ≈ V0(r) + Σm Ang1_m(cosθ1) * Ang2_m(cosθ2) * Vm(r)

    The CGL-CGL isotropic core is not supplied by martini_potential for this
    schema; V0(r) is therefore the true angular mean, not a residual baseline.
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

    total_samples = n_radial * n_angle * n_angle * azimuthal_count * azimuthal_count
    print(f"  Sampling CG↔CG energy: {n_radial} radial × {n_angle}² angular "
          f"× {azimuthal_count}² azimuthal = {total_samples} samples, "
          f"relax={relax_steps}, modes={n_modes}")

    for ir, r_nm in enumerate(r_values):
        for ia1 in range(n_angle):
            for ia2 in range(n_angle):
                energy_sum = 0.0
                for ip1 in range(azimuthal_count):
                    dir1 = dirs1[ia1, ip1]
                    for ip2 in range(azimuthal_count):
                        dir2 = dirs2[ia2, ip2]
                        energy_sum += _compute_cg_pair_energy(
                            r_nm, dir1, dir2, ref_nm, ref_nm,
                            bead_types, bead_types,
                            bead_charges, bead_charges,
                            pair_params,
                            relax_steps=relax_steps,
                            dist_min_nm=dist_min_nm,
                        )
                energy_grid[ir, ia1, ia2] = energy_sum / (azimuthal_count * azimuthal_count)

    if energy_cap_kj_mol is not None and energy_cap_kj_mol > 0.0:
        fit_grid = np.clip(
            energy_grid,
            -float(energy_cap_kj_mol),
            float(energy_cap_kj_mol),
        )
    else:
        fit_grid = energy_grid

    v_radial = np.mean(fit_grid, axis=(1, 2))
    residual = fit_grid - v_radial[:, None, None]
    if residual_cap_kj_mol is not None and residual_cap_kj_mol > 0.0:
        residual_fit = np.clip(residual, -float(residual_cap_kj_mol), float(residual_cap_kj_mol))
    else:
        residual_fit = residual

    mode_ang1 = []
    mode_ang2 = []
    mode_radial = []
    residual_remaining = residual_fit.copy()
    for mode_idx in range(int(n_modes)):
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
        v_ref = v_all[ref_idx].copy()
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

        ang_sym = 0.5 * (ang1 + ang2)
        ang1 = ang_sym.copy()
        ang2 = ang_sym.copy()
        max_abs = max(float(np.max(np.abs(ang1))), float(np.max(np.abs(ang2))), 1e-15)
        ang1 /= max_abs
        ang2 /= max_abs

        basis = np.outer(ang1, ang2)
        denom = float(np.dot(basis.ravel(), basis.ravel()))
        vm = np.zeros(n_radial, dtype=np.float64)
        if denom > 1e-15:
            for ir in range(n_radial):
                vm[ir] = float(np.dot(residual_remaining[ir].ravel(), basis.ravel()) / denom)
            if residual_cap_kj_mol is not None and residual_cap_kj_mol > 0.0:
                cap = float(residual_cap_kj_mol)
                vm = np.clip(vm, -cap, cap)
            residual_remaining -= basis[None, :, :] * vm[:, None, None]

        mode_ang1.append(ang1)
        mode_ang2.append(ang2)
        mode_radial.append(vm)

    recon = v_radial[:, None, None]
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        recon = recon + np.outer(ang1, ang2)[None, :, :] * vm[:, None, None]
    target_grid = v_radial[:, None, None] + residual_fit
    rms_error = float(np.sqrt(np.mean((target_grid - recon) ** 2)))

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
    param_parts = [v0_knots]
    ang1_knots_all = []
    ang2_knots_all = []
    vm_knots_all = []
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        ang1_knots = _fit_angular_bspline(t_angular, ang1, n_knot_angular, smooth=0.01)
        ang2_knots = _fit_angular_bspline(t_angular, ang2, n_knot_angular, smooth=0.01)
        vm_knots = _fit_radial_bspline(t_radial_ang, vm, rad_knot_vector, smooth=0.01) * inv_conv
        ang1_knots_all.append(ang1_knots)
        ang2_knots_all.append(ang2_knots)
        vm_knots_all.append(vm_knots)
        param_parts.extend([ang1_knots, ang2_knots, vm_knots])

    interaction_param = np.concatenate(param_parts)

    return {
        "v0_knots": v0_knots,
        "ang1_knots": np.asarray(ang1_knots_all),
        "ang2_knots": np.asarray(ang2_knots_all),
        "v_mode_knots": np.asarray(vm_knots_all),
        "interaction_param": interaction_param,
        "rms_error": rms_error,
        "v_radial_raw": v_radial,
        "residual_cap_kj_mol": residual_cap_kj_mol,
        "energy_cap_kj_mol": energy_cap_kj_mol,
        "ang1_raw": np.asarray(mode_ang1),
        "ang2_raw": np.asarray(mode_ang2),
        "v_angular_raw": np.asarray(mode_radial),
        "r_values_nm": r_values,
        "cos_theta_grid": cos_theta_grid,
        "knot_spacing_ang": knot_spacing,
        "cutoff_ang": float((n_knot_radial - 2) * knot_spacing),
        "n_modes": int(len(mode_radial)),
        "n_radial": int(n_knot_radial),
        "n_angular": int(n_knot_angular),
    }


def _fit_cg_lipid_sc_quadspline(
    ref_bead_positions_nm: np.ndarray,
    cg_bead_types: list,
    cg_bead_charges: list,
    target_type: str,
    target_charge: float,
    rotamer_centers_nm: list,
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
    relax_steps: int = 120,
    n_modes: int = 4,
    n_knot_radial: int = 14,
    knot_spacing_ang: float = 1.4,
    residual_cap_kj_mol: float | None = 250.0,
    energy_cap_kj_mol: float | None = 500.0,
) -> dict:
    """Fit full multimode B-spline params for sidechain ↔ CG lipid interactions.

    The parameter layout matches _fit_cg_lipid_quadspline:
      V0(r) + sum_m Ang1_m(a_sc) * Ang2_m(a_cg) * Vm(r)
    where source1 is the sidechain vector and source2 is the CG lipid vector.
    """
    r_values = np.asarray(_linspace(r_min_nm, r_max_nm, r_count), dtype=np.float64)
    cos_theta_grid = np.asarray(_linspace(-1.0, 1.0, cos_theta_count), dtype=np.float64)
    n_angle = cos_theta_count
    n_radial = r_count

    phi_values = np.linspace(0.0, 2.0 * np.pi, azimuthal_count, endpoint=False)

    ref_nm = _canonicalize_lipid_reference_to_z(ref_bead_positions_nm)
    cb_vec = np.asarray(cb_vector_unit, dtype=np.float64)
    cb_vec /= max(float(np.linalg.norm(cb_vec)), 1e-12)
    sc_to_cg_dirs = _directions_with_dot_np(-cb_vec, cos_theta_grid, phi_values)
    n_rotamer = len(rotamer_centers_nm)

    # Pre-build SC bead position arrays for each rotamer
    rotamer_sc_positions = []
    for irot in range(n_rotamer):
        center = np.asarray(rotamer_centers_nm[irot], dtype=np.float64)
        if center.ndim == 1:
            rotamer_sc_positions.append(center.reshape(1, 3))
        else:
            rotamer_sc_positions.append(center)

    print(f"  Sampling CG↔SC energy for target={target_type}: "
          f"{n_radial} radial × {n_angle}² angular × {azimuthal_count} azimuthal, "
          f"relax={relax_steps}, modes={n_modes}")

    # Energy grid: (n_radial, n_angle_sc, n_angle_cg), matching runtime source order.
    energy_grid = np.zeros((n_radial, n_angle, n_angle), dtype=np.float64)

    for ir, r_nm in enumerate(r_values):
        for ia_sc in range(n_angle):
            for ia_cg in range(n_angle):
                energy_sum = 0.0
                weight_sum = 0.0
                for ip in range(azimuthal_count):
                    dir_to_cg = sc_to_cg_dirs[ia_sc, ip]
                    cg_com = r_nm * dir_to_cg
                    dir_cg = _direction_with_dot_np(dir_to_cg, cos_theta_grid[ia_cg], 0.0)

                    R_cg = _rotation_to_align_z_np(dir_cg)
                    cg_positions = cg_com[None, :] + (R_cg @ ref_nm.T).T

                    for irot in range(n_rotamer):
                        w = rotamer_weights[irot]
                        if w <= 0.0:
                            continue

                        sc_positions = rotamer_sc_positions[irot]

                        if relax_steps > 0:
                            _, pair_energy = _relax_lipid_against_points(
                                cg_positions, sc_positions,
                                cg_bead_types, sc_bead_types,
                                cg_bead_charges, sc_bead_charges,
                                pair_params, n_steps=relax_steps,
                            )
                        else:
                            pair_energy, _, _ = _compute_pair_energy_and_gradient(
                                cg_positions, sc_positions, cg_bead_types, sc_bead_types,
                                cg_bead_charges, sc_bead_charges, pair_params,
                                dist_min_nm=0.20, soft_core_alpha=0.0,
                            )

                        energy_sum += w * pair_energy
                        weight_sum += w

                energy_grid[ir, ia_sc, ia_cg] = energy_sum / max(weight_sum, 1e-15)

    if energy_cap_kj_mol is not None and energy_cap_kj_mol > 0.0:
        fit_grid = np.clip(energy_grid, -float(energy_cap_kj_mol), float(energy_cap_kj_mol))
    else:
        fit_grid = energy_grid

    v_radial = np.mean(fit_grid, axis=(1, 2))
    residual = fit_grid - v_radial[:, None, None]
    if residual_cap_kj_mol is not None and residual_cap_kj_mol > 0.0:
        residual_fit = np.clip(residual, -float(residual_cap_kj_mol), float(residual_cap_kj_mol))
    else:
        residual_fit = residual

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
            if residual_cap_kj_mol is not None and residual_cap_kj_mol > 0.0:
                vm = np.clip(vm, -float(residual_cap_kj_mol), float(residual_cap_kj_mol))
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
    param_parts = [v0_knots]
    ang1_knots_all = []
    ang2_knots_all = []
    vm_knots_all = []
    for ang1, ang2, vm in zip(mode_ang1, mode_ang2, mode_radial):
        ang1_knots = _fit_angular_bspline(t_angular, ang1, n_knot_angular, smooth=0.01)
        ang2_knots = _fit_angular_bspline(t_angular, ang2, n_knot_angular, smooth=0.01)
        vm_knots = _fit_radial_bspline(t_radial_ang, vm, rad_knot_vector, smooth=0.01) * inv_conv
        ang1_knots_all.append(ang1_knots)
        ang2_knots_all.append(ang2_knots)
        vm_knots_all.append(vm_knots)
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
        "residual_cap_kj_mol": residual_cap_kj_mol,
        "energy_cap_kj_mol": energy_cap_kj_mol,
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
    }


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


def _cgl_effective_lj_sigma_cap_nm() -> float:
    raw = os.environ.get("UPSIDE_CG_LIPID_MAX_EFFECTIVE_SIGMA_NM", "0.9")
    try:
        cap = float(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid UPSIDE_CG_LIPID_MAX_EFFECTIVE_SIGMA_NM={raw!r}") from exc
    if cap <= 0.0:
        return float("inf")
    return cap


def _compute_cgl_effective_lj_params(
    ref_bead_positions_nm: np.ndarray,
    bead_types: list,
    pair_params: dict,
    r_min_nm: float = 0.2,
    r_max_nm: float = 3.0,
    n_radial: int = 100,
    n_orientations: int = 200,
) -> dict:
    """Compute effective LJ (sigma_nm, epsilon_kj_mol) for CGL syntype with each target type.

    Orientation-averages the total LJ interaction between all 14 DOPC beads and a
    point particle of each target type, then fits effective LJ(12,6) parameters.

    Returns dict: target_type → dict(sigma_nm=float, epsilon_kj_mol=float)
    """
    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    n_beads = ref_nm.shape[0]
    if n_beads != len(bead_types):
        raise ValueError(f"Bead count mismatch: {n_beads} positions vs {len(bead_types)} types")

    all_types = sorted(set(k[0] for k in pair_params.keys()) | set(k[1] for k in pair_params.keys()))

    directions = _fibonacci_sphere(n_orientations)
    dir_array = np.asarray(directions, dtype=np.float64)

    r_values = np.linspace(r_min_nm, r_max_nm, n_radial)
    r6 = r_values ** 6
    r12 = r_values ** 12

    # Cache self-interaction params for combination-rule fallback
    self_params_cache: dict = {}

    def _resolve_self_params(bt: str) -> dict:
        if bt not in self_params_cache:
            p = pair_params.get((bt, bt))
            if p is None:
                raise RuntimeError(f"No self-interaction LJ params for type {bt}")
            self_params_cache[bt] = p
        return self_params_cache[bt]

    sigma_cap_nm = _cgl_effective_lj_sigma_cap_nm()
    result: dict = {}
    for target_type in all_types:
        bead_params = []
        for bt in bead_types:
            key = (bt, target_type)
            params = pair_params.get(key)
            if params is None:
                key = (target_type, bt)
                params = pair_params.get(key)
            if params is None:
                p_self = _resolve_self_params(bt)
                p_tgt = _resolve_self_params(target_type)
                sig = (p_self["sigma_nm"] + p_tgt["sigma_nm"]) / 2.0
                eps = math.sqrt(p_self["epsilon_kj_mol"] * p_tgt["epsilon_kj_mol"])
                params = {"sigma_nm": sig, "epsilon_kj_mol": eps}
            bead_params.append(params)

        avg_energy = np.zeros(n_radial)
        for ir, r in enumerate(r_values):
            energy_sum = 0.0
            for dir_vec in dir_array:
                target_pos = r * dir_vec
                total_lj = 0.0
                for b in range(n_beads):
                    d = target_pos - ref_nm[b]
                    dist = float(math.sqrt(float(np.dot(d, d))))
                    if dist < 0.001:
                        dist = 0.001
                    sig = bead_params[b]["sigma_nm"]
                    eps = bead_params[b]["epsilon_kj_mol"]
                    sr = sig / dist
                    sr2 = sr * sr
                    sr6 = sr2 * sr2 * sr2
                    total_lj += 4.0 * eps * (sr6 * sr6 - sr6)
                energy_sum += total_lj
            avg_energy[ir] = energy_sum / n_orientations

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

        uncapped_sigma_eff = max(0.1, min(sigma_eff, 5.0))
        sigma_eff = min(uncapped_sigma_eff, sigma_cap_nm)
        epsilon_eff = max(0.01, min(epsilon_eff, 100.0))
        result[target_type] = {
            "sigma_nm": sigma_eff,
            "epsilon_kj_mol": epsilon_eff,
            "uncapped_sigma_nm": uncapped_sigma_eff,
            "sigma_cap_nm": sigma_cap_nm,
        }

    return result


def _compute_cgl_self_lj(
    ref_bead_positions_nm: np.ndarray,
    bead_types: list,
    pair_params: dict,
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
    n_orientations: int = 100,
) -> tuple:
    """Fit effective LJ(12,6) for CGL↔CGL by orientation-averaging two 14-bead lipids.

    Unlike _compute_cgl_effective_lj_params which computes CGL↔point_particle,
    this computes the full CGL↔CGL interaction with both lipids as extended objects.
    Returns (eps_eup, sig_ang) for the effective LJ.
    """
    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    n_beads = ref_nm.shape[0]

    # Use a modest number of orientation pairs for efficiency
    directions = _fibonacci_sphere(min(n_orientations, 100))
    dir_array = np.asarray(directions, dtype=np.float64)
    n_dir = len(dir_array)

    # Sample COM-COM distances; use fewer points for speed
    r_values = np.linspace(0.5, 2.0, 40)
    r6 = r_values ** 6
    r12 = r_values ** 12

    def _rotation_to_align_z_np_local(dir_vec):
        z_axis = dir_vec / np.linalg.norm(dir_vec)
        if abs(z_axis[0]) < 0.99:
            x_axis = np.cross([1.0, 0.0, 0.0], z_axis)
        else:
            x_axis = np.cross([0.0, 1.0, 0.0], z_axis)
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)
        return np.array([x_axis, y_axis, z_axis]).T

    avg_energy = np.zeros(len(r_values))
    # Use shuffled indices for random pairing (avoids O(n²) orientation loop)
    rng = np.random.default_rng(42)

    for ir, r_nm in enumerate(r_values):
        energy_sum = 0.0
        count = 0
        # Randomly pair orientations for efficiency
        idx1 = np.arange(n_dir)
        idx2 = rng.permutation(n_dir)
        for i in range(n_dir):
            dir1 = dir_array[idx1[i]]
            dir2 = dir_array[idx2[i]]
            R1 = _rotation_to_align_z_np_local(dir1)
            pos1 = (R1 @ ref_nm.T).T  # (14, 3) — CG lipid 1 at origin

            R2 = _rotation_to_align_z_np_local(dir2)
            offset2 = np.array([r_nm, 0.0, 0.0])
            pos2 = offset2[None, :] + (R2 @ ref_nm.T).T

            total_lj = 0.0
            for b1 in range(n_beads):
                for b2 in range(n_beads):
                    dx = pos2[b2] - pos1[b1]
                    dist = float(np.sqrt(np.dot(dx, dx)))
                    if dist < 0.30:
                        dist = 0.30  # soft core
                    bt1 = bead_types[b1]
                    bt2 = bead_types[b2]
                    key = (bt1, bt2)
                    params = pair_params.get(key)
                    if params is None:
                        params = pair_params.get((bt2, bt1))
                    if params is None:
                        raise RuntimeError(f"Missing pair params for ({bt1}, {bt2})")
                    sig = params["sigma_nm"]
                    eps = params["epsilon_kj_mol"]
                    sr = sig / dist
                    sr6 = sr ** 6
                    total_lj += 4.0 * eps * (sr6 * sr6 - sr6)

            energy_sum += total_lj
            count += 1

        avg_energy[ir] = energy_sum / max(count, 1)

    # Fit LJ(12,6): E(r) = 4*eps*((sig/r)^12 - (sig/r)^6)
    # → E*r^12 = 4*eps*sig^12 - 4*eps*sig^6 * r^6
    # → y = A - B*x where y=E*r^12, x=r^6, A=4*eps*sig^12, B=4*eps*sig^6
    attractive_mask = avg_energy < 0.0
    if attractive_mask.sum() >= 3:
        y = avg_energy[attractive_mask] * r12[attractive_mask]
        x = r6[attractive_mask]
        fit = np.polyfit(x, y, 1)
        B_fit = -float(fit[0])
        A_fit = float(fit[1])
        if A_fit > 0.0 and B_fit > 0.0:
            sigma_nm = (A_fit / B_fit) ** (1.0 / 6.0)
            epsilon_kj = B_fit * B_fit / (4.0 * A_fit)
        else:
            imin = int(np.argmin(avg_energy))
            sigma_nm = r_values[imin] / (2.0 ** (1.0 / 6.0))
            epsilon_kj = max(-float(avg_energy[imin]), 0.01)
    else:
        imin = int(np.argmin(avg_energy))
        sigma_nm = r_values[imin] / (2.0 ** (1.0 / 6.0))
        epsilon_kj = max(-float(avg_energy[imin]), 0.01)

    sigma_nm = max(0.2, min(sigma_nm, 2.0))
    epsilon_kj = max(0.01, min(epsilon_kj, 50.0))

    # Physical override: the orientation-averaged LJ fit with soft-core produces
    # an artificially small sigma (~0.4 nm) because bead-bead overlap is
    # suppressed.  Use a physically motivated sigma based on the lipid's
    # spatial extent: max bead-COM distance ≈ 0.5 nm, so two lipids cannot
    # approach closer than ~1.0 nm without bead overlap.
    # sigma = 1.0 nm (10 Å) provides a repulsive core below ~1.0 nm.
    # epsilon = 2.9 kJ/mol (1.0 E_up) keeps the isotropic attraction weak;
    # the quadspline captures orientation-dependent interactions.
    sigma_nm = 1.0
    epsilon_kj = ENERGY_CONVERSION_KJ_PER_EUP  # 1.0 E_up ≈ 2.915 kJ/mol

    eps_eup = epsilon_kj / energy_conv
    sig_ang = sigma_nm * length_conv

    return eps_eup, sig_ang


def _extend_particles_group_for_cgl(
    h5: h5py.File,
    pair_params: dict,
    effective_lj: dict,
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
    ref_bead_positions_nm: np.ndarray = None,
    bead_types: list = None,
) -> None:
    """Add CGL syntype to existing particles group with combined energy grids.

    Reads the current particles group, appends CGL to type_order and type_charge,
    adds new (eps,sig,qq) triples for CGL↔all_types, and regenerates all combined
    energy grid arrays.
    """
    pg = h5["particles"]
    type_order_raw = pg["type_order"][:]
    type_order = [
        x.decode("ascii") if isinstance(x, bytes) else str(x) for x in type_order_raw
    ]
    if "CGL" in type_order:
        return

    existing_charges = pg["type_charge"][:].astype(np.float64)
    existing_eps = pg["unique_eps_eup"][:].astype(np.float64)
    existing_sig = pg["unique_sig_ang"][:].astype(np.float64)
    existing_qq = pg["unique_charge_product"][:].astype(np.float64)
    existing_grids = pg["combined_energy_grids"][:].astype(np.float64)

    r_min = float(pg.attrs["r_min_ang"])
    r_max = float(pg.attrs["r_max_ang"])
    n_pts = int(pg.attrs["n_points"])

    coulomb_k_eup = COULOMB_K_DRY_KJ_NM * length_conv / energy_conv

    # Build new triples and grids for CGL↔each_existing_type
    new_eps_list = list(existing_eps)
    new_sig_list = list(existing_sig)
    new_qq_list = list(existing_qq)
    new_grids_list = [existing_grids[i] for i in range(len(existing_eps))]

    existing_triple_set = set(
        (float(existing_eps[i]), float(existing_sig[i]), float(existing_qq[i]))
        for i in range(len(existing_eps))
    )

    cgl_charge = 0.0
    for ti_idx, ti_name in enumerate(type_order):
        tgt_charge = float(existing_charges[ti_idx])
        qq = cgl_charge * tgt_charge
        eff = effective_lj.get(ti_name)
        if eff is None:
            raise RuntimeError(f"No effective LJ params for CGL↔{ti_name}")

        eps_eup = float(eff["epsilon_kj_mol"]) / energy_conv
        sig_ang = float(eff["sigma_nm"]) * length_conv

        triple_key = (eps_eup, sig_ang, qq)
        if triple_key in existing_triple_set:
            continue
        existing_triple_set.add(triple_key)

        grid = np.zeros(n_pts, dtype=np.float64)
        for i in range(n_pts):
            r = r_min + i * (r_max - r_min) / (n_pts - 1)
            r = max(r, 0.1 * sig_ang)
            r2 = r * r
            r6 = r2 * r2 * r2
            sig2 = sig_ang * sig_ang
            sig6 = sig2 * sig2 * sig2
            sig12 = sig6 * sig6
            if eps_eup != 0.0 and sig_ang != 0.0:
                lj = 4.0 * eps_eup * (sig12 / (r6 * r6) - sig6 / r6)
            else:
                lj = 0.0
            coul = coulomb_k_eup * qq / max(r, 1.0e-6) if abs(qq) > 1e-10 else 0.0
            grid[i] = lj + coul

        new_eps_list.append(eps_eup)
        new_sig_list.append(sig_ang)
        new_qq_list.append(qq)
        new_grids_list.append(grid)

    # Add CGL↔CGL self-interaction grid.
    # The effective LJ from orientation-averaging single-target interactions gives
    # sigma ~ 1.7 nm (the lipid diameter), which produces enormous repulsion at
    # typical bilayer COM distances (~0.8 nm). For CGL↔CGL, compute the
    # orientation-averaged pair energy directly and fit a proper effective LJ.
    cgl_self_eps, cgl_self_sig = _compute_cgl_self_lj(
        ref_bead_positions_nm=ref_bead_positions_nm,
        bead_types=bead_types,
        pair_params=pair_params,
        energy_conv=energy_conv,
        length_conv=length_conv,
    )
    cgl_self_qq = 0.0
    cgl_self_qq = 0.0

    triple_key = (cgl_self_eps, cgl_self_sig, cgl_self_qq)
    if triple_key not in existing_triple_set:
        existing_triple_set.add(triple_key)
        grid = np.zeros(n_pts, dtype=np.float64)
        for i in range(n_pts):
            r = r_min + i * (r_max - r_min) / (n_pts - 1)
            r = max(r, 0.1 * cgl_self_sig)
            r2 = r * r
            r6 = r2 * r2 * r2
            sig2 = cgl_self_sig * cgl_self_sig
            sig6 = sig2 * sig2 * sig2
            sig12 = sig6 * sig6
            lj = 4.0 * cgl_self_eps * (sig12 / (r6 * r6) - sig6 / r6)
            grid[i] = lj
        new_eps_list.append(cgl_self_eps)
        new_sig_list.append(cgl_self_sig)
        new_qq_list.append(cgl_self_qq)
        new_grids_list.append(grid)

    # Update type_order and type_charge
    new_type_order = type_order + ["CGL"]
    new_type_charges = np.concatenate([existing_charges, np.array([cgl_charge], dtype=np.float64)])

    # Rebuild combined arrays
    n_triples = len(new_eps_list)
    all_grids = np.zeros((n_triples, n_pts), dtype=np.float64)
    for idx, grid in enumerate(new_grids_list):
        all_grids[idx, :] = grid

    # Delete and recreate datasets
    for ds_name in ["type_order", "type_charge", "unique_eps_eup", "unique_sig_ang",
                     "unique_charge_product", "combined_energy_grids"]:
        del pg[ds_name]

    pg.create_dataset("type_order",
                       data=np.asarray([np.bytes_(x) for x in new_type_order], dtype="S8"))
    pg.create_dataset("type_charge", data=new_type_charges.astype(np.float32))
    pg.create_dataset("unique_eps_eup", data=np.asarray(new_eps_list, dtype=np.float64))
    pg.create_dataset("unique_sig_ang", data=np.asarray(new_sig_list, dtype=np.float64))
    pg.create_dataset("unique_charge_product", data=np.asarray(new_qq_list, dtype=np.float64))
    pg.create_dataset("combined_energy_grids", data=all_grids, dtype=np.float64)

    print(f"  particles: extended with CGL → {n_triples} combined triples, "
          f"{len(new_type_order)} types")


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
    r_min_nm: float = 0.30,
    r_max_nm: float = 0.70,
    r_count: int = 24,
    cos_theta_count: int = 13,
    azimuthal_count: int = 4,
) -> None:
    """Build CG lipid pair and CG lipid ↔ SC quadspline tables and store in HDF5."""
    if ref_bead_positions_nm is None:
        print("  cg_lipid_table: no reference bead positions provided, skipping")
        return

    if bead_types is None:
        bead_types = _DOPC_BEAD_TYPES
    if bead_charges is None:
        bead_charges = [_infer_type_charge(bt) for bt in bead_types]

    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if ref_nm.shape != (14, 3):
        raise ValueError(f"ref_bead_positions_nm must be (14, 3), got {ref_nm.shape}")

    print("\n=== CG Lipid Table Building ===")

    # Resolution control: "coarse", "medium", or "fine" (default)
    # Coarse/medium give reasonable accuracy for testing at much lower cost.
    _res = os.environ.get("UPSIDE_CG_LIPID_RESOLUTION", "fine").strip().lower()
    if _res == "coarse":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 8, 5, 2, 8, 5, 2
    elif _res == "medium":
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 12, 7, 2, 12, 7, 2
    else:
        _cg_r, _cg_ct, _cg_az, _sc_r, _sc_ct, _sc_az = 16, 7, 2, 16, 9, 4
    print(f"  CG lipid resolution: {_res} "
          f"(CG: {_cg_r}r×{_cg_ct}²θ×{_cg_az}²φ, "
          f"SC: {_sc_r}r×{_sc_ct}²θ×{_sc_az}²φ)")

    # CG ↔ CG quadspline. Keep the particles table as the isotropic core and
    # fit only the residual directional correction; moderate bead relaxation
    # to soften unphysical rigid-body overlaps while preserving repulsive core.
    relax_steps = 50
    result_cg = _fit_cg_lipid_quadspline(
        ref_bead_positions_nm=ref_nm,
        bead_types=bead_types,
        bead_charges=bead_charges,
        pair_params=pair_params,
        r_min_nm=0.70,
        r_max_nm=1.68,
        r_count=min(r_count, _cg_r),
        cos_theta_count=min(cos_theta_count, _cg_ct),
        azimuthal_count=_cg_az,
        relax_steps=relax_steps,
        dist_min_nm=0.25,
        knot_spacing_ang=1.4,
        residual_cap_kj_mol=250.0,
        n_modes=4,
        n_knot_radial=14,
    )
    print(f"  CG↔CG: RMS error = {result_cg['rms_error']:.4f} kJ/mol, "
          f"modes = {result_cg['n_modes']}, "
          f"max|V0| = {float(np.max(np.abs(result_cg['v_radial_raw']))):.3f} kJ/mol")
    cg_pair_scale = float(os.environ.get("UPSIDE_CG_LIPID_PAIR_SCALE", "1.0"))
    print(f"  CG↔CG table scale kappa = {cg_pair_scale:.6g}")

    # CG ↔ SC quadspline
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    residue_map = _load_martini_forcefield(martinize_path, forcefield_name)

    residues = [r for r in active_residue_names
                if r in residue_map and residue_map[r] and r in orientation_map]
    sc_n_modes = int(os.environ.get("UPSIDE_CG_LIPID_SC_N_MODES", "4"))
    sc_n_radial = 14
    sc_n_angular = 15
    sc_knot_spacing_ang = 1.4
    sc_taper_width_ang = sc_knot_spacing_ang
    sc_cutoff_ang = float(r_max_nm * 10.0 + sc_taper_width_ang)
    sc_n_param = sc_n_radial + sc_n_modes * (2 * sc_n_angular + sc_n_radial)
    sc_rms_error = np.zeros(len(residues), dtype=np.float32)

    if not residues:
        print("  cg_lipid_sc: no active residues with sidechains, skipping")
        n_sc_types = 0
        interaction_param_sc = np.zeros((0, 1, sc_n_param), dtype=np.float64)
        sc_residue_names = []
    else:
        cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
        cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)
        relax_steps = 120

        n_sc_types = len(residues)
        interaction_param_sc = np.zeros((n_sc_types, 1, sc_n_param), dtype=np.float64)
        sc_residue_names = []

        for ri, residue in enumerate(residues):
            sc_bead_types = residue_map[residue]
            sc_bead_charges = [_infer_type_charge(bt) for bt in sc_bead_types]
            orientation = orientation_map[residue]

            result_sc = _fit_cg_lipid_sc_quadspline(
                ref_bead_positions_nm=ref_nm,
                cg_bead_types=bead_types,
                cg_bead_charges=bead_charges,
                target_type="CGL",
                target_charge=0.0,
                rotamer_centers_nm=orientation["center_nm"],
                rotamer_weights=orientation["weight"],
                sc_bead_types=sc_bead_types,
                sc_bead_charges=sc_bead_charges,
                pair_params=pair_params,
                cb_anchor_nm=cb_anchor_nm,
                cb_vector_unit=cb_vector_unit,
                r_min_nm=r_min_nm,
                r_max_nm=r_max_nm,
                r_count=min(r_count, _sc_r),
                cos_theta_count=min(cos_theta_count, _sc_ct),
                azimuthal_count=min(azimuthal_count, _sc_az),
                relax_steps=relax_steps,
                n_modes=sc_n_modes,
                n_knot_radial=sc_n_radial,
                knot_spacing_ang=sc_knot_spacing_ang,
            )
            interaction_param_sc[ri, 0, :] = result_sc["interaction_param"]
            sc_rms_error[ri] = np.float32(result_sc["rms_error"])
            sc_cutoff_ang = min(sc_cutoff_ang, float(result_sc["cutoff_ang"]))
            sc_residue_names.append(residue)
            print(f"  CG↔SC({residue}): RMS error = {result_sc['rms_error']:.4f} kJ/mol, "
                  f"modes = {result_sc['n_modes']}")

    # Store in HDF5
    cg_grp = h5.create_group("cg_lipid_table")

    # CG ↔ CG pair
    cg_pair_grp = cg_grp.create_group("cg_lipid_pair")
    pair_param = result_cg["interaction_param"].astype(np.float64, copy=True)
    nr = int(result_cg["n_radial"])
    na = int(result_cg["n_angular"])
    pair_param[:nr] *= cg_pair_scale
    mode_stride = 2 * na + nr
    for mode_idx in range(int(result_cg["n_modes"])):
        radial_start = nr + mode_idx * mode_stride + 2 * na
        pair_param[radial_start:radial_start + nr] *= cg_pair_scale
    pair_param = pair_param.astype(np.float32)
    cg_pair_grp.create_dataset(
        "interaction_param",
        data=pair_param.reshape(1, 1, pair_param.size),
    )
    cg_pair_grp.attrs["n_cg_types"] = 1
    cg_pair_grp.attrs["rms_error_kj_mol"] = np.float32(result_cg["rms_error"])
    cg_pair_grp.attrs["schema"] = "cg_lipid_quadspline_v3"
    cg_pair_grp.attrs["radial_mode"] = "full_multimode"
    cg_pair_grp.attrs["angle_convention"] = "ang1=-n1_dot_n12;ang2=n2_dot_n12"
    cg_pair_grp.attrs["fit_relax_steps"] = np.int32(0)
    cg_pair_grp.attrs["fit_r_min_nm"] = np.float32(0.70)
    cg_pair_grp.attrs["fit_r_max_nm"] = np.float32(1.68)
    cg_pair_grp.attrs["n_modes"] = np.int32(result_cg["n_modes"])
    cg_pair_grp.attrs["n_radial"] = np.int32(result_cg["n_radial"])
    cg_pair_grp.attrs["n_angular"] = np.int32(result_cg["n_angular"])
    cg_pair_grp.attrs["pair_scale"] = np.float32(cg_pair_scale)
    cg_pair_grp.attrs["knot_spacing_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["cutoff_ang"] = np.float32(result_cg["cutoff_ang"])
    cg_pair_grp.attrs["taper_width_ang"] = np.float32(result_cg["knot_spacing_ang"])
    cg_pair_grp.attrs["residual_cap_kj_mol"] = np.float32(result_cg["residual_cap_kj_mol"])
    cg_pair_grp.attrs["energy_cap_kj_mol"] = np.float32(result_cg["energy_cap_kj_mol"])
    cg_pair_grp.create_dataset("v0_raw_kj_mol", data=result_cg["v_radial_raw"].astype(np.float32))
    cg_pair_grp.create_dataset("mode_radial_raw_kj_mol", data=result_cg["v_angular_raw"].astype(np.float32))

    # CG ↔ SC
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
    cg_sc_grp.attrs["fit_relax_steps"] = np.int32(120 if n_sc_types else 0)
    cg_sc_grp.attrs["fit_r_min_nm"] = np.float32(r_min_nm)
    cg_sc_grp.attrs["fit_r_max_nm"] = np.float32(r_max_nm)
    cg_sc_grp.attrs["n_modes"] = np.int32(sc_n_modes)
    cg_sc_grp.attrs["n_radial"] = np.int32(sc_n_radial)
    cg_sc_grp.attrs["n_angular"] = np.int32(sc_n_angular)
    cg_sc_grp.attrs["knot_spacing_ang"] = np.float32(sc_knot_spacing_ang)
    cg_sc_grp.attrs["cutoff_ang"] = np.float32(sc_cutoff_ang)
    cg_sc_grp.attrs["taper_width_ang"] = np.float32(sc_taper_width_ang)
    cg_sc_grp.attrs["residual_cap_kj_mol"] = np.float32(250.0)
    cg_sc_grp.attrs["energy_cap_kj_mol"] = np.float32(500.0)

    print(
        f"  Stored: CG↔CG (1×1×{pair_param.size}), "
        f"CG↔SC ({n_sc_types}×1×{interaction_param_sc.shape[-1]}) in {h5.filename}"
    )

    # Extend particles group with CGL syntype for isotropic LJ interactions
    # (CG lipid ↔ ions/water/BB/env beads)
    print("  Computing effective CGL LJ parameters for particles table extension...")
    effective_lj = _compute_cgl_effective_lj_params(
        ref_bead_positions_nm=ref_nm,
        bead_types=bead_types,
        pair_params=pair_params,
    )
    _extend_particles_group_for_cgl(h5, pair_params, effective_lj,
                                     ref_bead_positions_nm=ref_nm,
                                     bead_types=bead_types)

    # Store effective LJ parameters so convert_stage() can read them back
    # instead of recomputing with potentially different lipid conformations.
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
    eff_grp.attrs["max_effective_sigma_nm"] = np.float32(_cgl_effective_lj_sigma_cap_nm())
# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def build_martini_tables(
    output_path: Path,
    dry_ff_path: Path,
    martinize_path: Path,
    sidechain_lib_path: Path,
    forcefield_name: str = "martini22",
    active_residue_names: List[str] | None = None,
    active_atom_types: Set[str] | None = None,
    r_min_nm: float = 0.25,
    r_max_nm: float = 1.20,
    r_count: int = 96,
    direction_count: int = 24,
    cos_theta_count: int = 13,
    cg_lipid_config: dict | None = None,
) -> Path:
    output_path = Path(output_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    martinize_path = Path(martinize_path).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path).expanduser().resolve()

    atomtypes, pair_params = _parse_dry_forcefield(dry_ff_path)

    if active_atom_types is None:
        active_atom_types = set(atomtypes)
    if active_residue_names is None:
        active_residue_names = list(CANONICAL_RESIDUES)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(output_path, "w") as h5:
        h5.attrs["schema"] = "martini_combined_v1"
        h5.attrs["created_for_run"] = str(output_path.parent.name)

        print(f"Building martini tables -> {output_path}")

        _build_particles_group(
            h5, atomtypes, pair_params, active_atom_types,
        )

        residue_map = _load_martini_forcefield(martinize_path, forcefield_name)

        # Determine active target types from active atom types (non-protein env types)
        sc_active_targets = sorted(
            t for t in active_atom_types if t in atomtypes
        )
        if sc_active_targets:
            _build_sc_table_group(
                h5, residue_map, pair_params, sidechain_lib_path,
                active_residue_names=list(active_residue_names),
                active_target_types=sc_active_targets,
                r_min_nm=r_min_nm, r_max_nm=r_max_nm, r_count=r_count,
                direction_count=direction_count, cos_theta_count=cos_theta_count,
            )
        else:
            print("  sc_table: no active target types, skipping")

        # Build CG lipid tables if configuration is provided
        if cg_lipid_config is not None:
            _build_cg_lipid_tables(
                h5,
                pair_params=pair_params,
                sidechain_lib_path=sidechain_lib_path,
                martinize_path=martinize_path,
                forcefield_name=forcefield_name,
                active_residue_names=list(active_residue_names),
                ref_bead_positions_nm=cg_lipid_config.get("ref_bead_positions_nm"),
                bead_types=cg_lipid_config.get("bead_types"),
                bead_charges=cg_lipid_config.get("bead_charges"),
                r_min_nm=r_min_nm,
                r_max_nm=min(r_max_nm, 0.70),  # quadspline cutoff ~7 Å
                r_count=r_count,
                cos_theta_count=cos_theta_count,
            )

    print(f"Built {output_path}")
    return output_path
