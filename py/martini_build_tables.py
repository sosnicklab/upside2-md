#!/usr/bin/env python
"""Build dry-MARTINI spline tables from ITP-derived parameters."""

import math
import os
import importlib.util
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Set, Tuple

import h5py
import numpy as np

from martini_itp_reader import (
    infer_charge_from_atomtype,
    load_martini_forcefield,
    parse_dry_forcefield,
)

COULOMB_K_DRY_KJ_NM = 138.935458 / 15.0
ENERGY_CONVERSION_KJ_PER_EUP = 2.914952774272
LENGTH_CONVERSION_A_PER_NM = 10.0
ANGSTROM_TO_NM = 0.1
DEFAULT_PRODUCTION_TEMP_UPSIDE = 0.8647
DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE = float(
    os.environ.get("UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE", "25.0")
)
PARTICLES_GRID_N = 1000
# The grid starts above zero because the dry-MARTINI pair potential diverges at r = 0 and the table must
# BE that potential, not a variant of it. 0.3 A is far inside anything reachable: the repulsive core there
# is ~5e17 E_up/A for the bulkiest bead pair, three orders above any other force in the system. Measured
# consequence of the previous choice (r_min = 0 with the radius floored at 0.1*sigma, which made the
# potential CONSTANT and therefore force-free below 0.47-0.60 A): pairs reached 0.0355 A and felt nothing,
# ~1440 times per 9 h run, while the true LJ force there is 5.8e29 E_up/A. See findings 92.
PARTICLES_R_MIN_A = 0.3
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

    # MARTINI 2 nonbonded contract: BOTH terms must reach the cutoff smoothly.
    #
    # A bare truncation leaves a step in the energy at r_c, and for a charged pair that step is
    # k*qq/r_c = 2.65 E_up at r_c = 12 A -- about 3.5 kT. Every time such a pair crosses the cutoff the
    # integrator gains or loses that much, so energy is not conserved and the error accumulates in
    # whichever degrees of freedom keep crossing. With hundreds of mobile ions in an implicit-solvent box
    # that is a continuous stochastic energy source, and it is what let a glpG micelle replica run away
    # (findings 79-80). This is the correction flagged as outstanding in findings 75.
    #
    # Coulomb -> reaction field with a conducting boundary (eps_rf -> infinity), the MARTINI convention:
    #     E(r) = k*qq * [ 1/r + k_rf*r^2 - c_rf ],  k_rf = 1/(2 r_c^3),  c_rf = 3/(2 r_c)
    # which gives E(r_c) = 0 AND dE/dr(r_c) = 0 exactly. eps_r = 15 is already folded into coulomb_k_eup.
    # LJ -> potential-shift (GROMACS `vdw-modifier = Potential-shift`): subtract E_LJ(r_c). The residual
    # force at r_c is ~1e-2 E_up/A, three orders below the Coulomb step being removed.
    r_c = float(PARTICLES_R_MAX_A)
    k_rf = 1.0 / (2.0 * r_c ** 3)
    c_rf = 3.0 / (2.0 * r_c)

    for eps, sig in eps_sig_list:
        sig6 = sig ** 6
        sig12 = sig6 * sig6
        lj_shift = 4.0 * eps * (sig12 / r_c ** 12 - sig6 / r_c ** 6)
        for qq in qq_set:
            r = np.linspace(PARTICLES_R_MIN_A, PARTICLES_R_MAX_A, PARTICLES_GRID_N)
            r6 = r ** 6
            grid = 4.0 * eps * (sig12 / (r6 * r6) - sig6 / r6) - lj_shift
            if abs(qq) > 1e-10:
                grid = grid + coulomb_k_eup * qq * (1.0 / r + k_rf * r * r - c_rf)
            if not np.isfinite(grid).all():
                raise ValueError(
                    "combined grid for (eps=%g, sig=%g, qq=%g) is not finite over [%g, %g] A"
                    % (eps, sig, qq, PARTICLES_R_MIN_A, PARTICLES_R_MAX_A))
            triples.append((eps, sig, qq))
            grids.append(grid)

    # The spline table is a representation of the dry-MARTINI pair potential, never a variant of it, so
    # assert the equivalence rather than assume it: every grid point must equal the published functional
    # form -- potential-shifted LJ plus reaction-field Coulomb -- evaluated at that radius.
    for (eps, sig, qq), grid in zip(triples, grids):
        r = np.linspace(PARTICLES_R_MIN_A, PARTICLES_R_MAX_A, PARTICLES_GRID_N)
        reference = 4.0 * eps * ((sig / r) ** 12 - (sig / r) ** 6) \
            - 4.0 * eps * ((sig / r_c) ** 12 - (sig / r_c) ** 6)
        if abs(qq) > 1e-10:
            reference = reference + coulomb_k_eup * qq * (1.0 / r + k_rf * r * r - c_rf)
        # Scale the tolerance by the energy scale rather than by |reference| alone: the potential-shifted
        # LJ crosses zero just inside the cutoff, so a pure relative measure there amplifies float64
        # round-off (the two algebraically identical spellings of sigma^12/r^12 differ in the last bits)
        # into an apparent 1e-12 error on a quantity that is itself ~1e-15 E_up.
        scale = np.maximum(np.abs(reference), max(abs(eps), 1.0))
        deviation = np.abs(grid - reference) / scale
        if deviation.max() > 1e-10:
            raise ValueError(
                "tabulated grid for (eps=%g, sig=%g, qq=%g) departs from the analytic dry-MARTINI form "
                "by %.3e (scaled) at r = %.4f A" % (eps, sig, qq, deviation.max(),
                                                    float(r[int(np.argmax(deviation))])))

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
    # Convert native dry-MARTINI units (nm, kJ/mol) to Upside runtime units (Angstrom, E_up)
    # once here at build time, so the runtime h5/config stores Upside-unit values directly.
    g.create_dataset("grid_ang", data=ref_grid_nm * LENGTH_CONVERSION_A_PER_NM, dtype=np.float32)
    g.create_dataset("cos_theta_grid", data=ref_cos_grid, dtype=np.float32)
    g.create_dataset("rotamer_count", data=rc, dtype=np.float32)
    g.create_dataset("rotamer_probability_fixed", data=rpf, dtype=np.float32)
    for name, arr in [
        ("radial_energy_eup", rad_e),
        ("angular_energy_eup", ang_e),
        ("rotamer_radial_energy_eup", r_rad),
        ("rotamer_angular_energy_eup", r_ang),
        ("rotamer_full_energy_eup", r_full),
    ]:
        g.create_dataset(name, data=arr / ENERGY_CONVERSION_KJ_PER_EUP, dtype=np.float32)
    # angular_profile / rotamer_angular_profile are dimensionless cos-profile multipliers.
    for name, arr in [
        ("angular_profile", ang_p),
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

            # Same MARTINI 2 nonbonded contract as the particle-particle grids: LJ potential-shift and
            # reaction-field Coulomb, so energy and force both reach the cutoff smoothly. Truncating
            # bare leaves a step of k*q1q2/r_c at r_c -- ~3.5 kT for a charged pair at 1.2 nm -- which
            # breaks energy conservation every time a pair crosses (findings 75, 79-80).
            sr = sigma_nm / eff_dist
            sr2 = sr * sr
            sr6 = sr2 * sr2 * sr2
            lj = 4.0 * epsilon_kj * (sr6 * sr6 - sr6)
            if cutoff_nm is not None:
                src6 = (sigma_nm / float(cutoff_nm)) ** 6
                lj -= 4.0 * epsilon_kj * (src6 * src6 - src6)
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
                kq = COULOMB_K_DRY_KJ_NM * q1 * q2
                if cutoff_nm is not None:
                    r_c = float(cutoff_nm)
                    k_rf = 1.0 / (2.0 * r_c ** 3)
                    c_rf = 3.0 / (2.0 * r_c)
                    coul = kq * (1.0 / coul_eff + k_rf * coul_eff * coul_eff - c_rf)
                    dcoul_dr = kq * (-1.0 / (coul_eff * coul_eff) + 2.0 * k_rf * coul_eff)
                else:
                    coul = kq / coul_eff
                    dcoul_dr = -kq / (coul_eff * coul_eff)
                total_energy += coul
                unit = dx / dist
                grad1[i] -= dcoul_dr * unit
                grad2[j] += dcoul_dr * unit

    return total_energy, grad1, grad2


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


def build_martini_h5(
    output_path: Path,
    dry_ff_path: Path,
    martinize_path: Path,
    sidechain_lib_path: Path,
    forcefield_name: str = "martini22",
) -> None:
    """Generate martini.h5 with both the particle and sidechain groups."""
    output_path = Path(output_path).expanduser().resolve()
    dry_ff_path = Path(dry_ff_path).expanduser().resolve()
    martinize_path = Path(martinize_path).expanduser().resolve()
    sidechain_lib_path = Path(sidechain_lib_path).expanduser().resolve()
    atomtypes, pair_params = parse_dry_forcefield(dry_ff_path)
    residue_map = load_martini_forcefield(martinize_path, forcefield_name)
    martini_sidechain_offsets_nm = _load_martini_sidechain_offsets_nm(
        martinize_path, forcefield_name
    )

    def _writer(h5: h5py.File) -> None:
        h5.attrs["schema"] = "martini_ff_combined"
        h5.attrs["schema_particles"] = SCHEMA_PARTICLES
        h5.attrs["schema_sc"] = SCHEMA_SC
        _build_particles_group(h5, atomtypes, pair_params, set(atomtypes))
        _build_sc_table_group(
            h5, residue_map, pair_params, sidechain_lib_path,
            active_residue_names=list(CANONICAL_RESIDUES),
            active_target_types=atomtypes,
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path}")
