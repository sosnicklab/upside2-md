#!/usr/bin/env python3
"""Build a per-run martini.h5 with combined LJ+Coulomb spline tables.

Generates particle-particle combined energy grids and sidechain-environment
orientation-aware combined energy tables, trimmed to only the types and
residues actually present in the current system.
"""

from __future__ import annotations

import importlib.util
import math
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

    print(f"Built {output_path}")
    return output_path
