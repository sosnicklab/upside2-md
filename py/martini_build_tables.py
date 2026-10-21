#!/usr/bin/env python
"""Build dry-MARTINI spline tables from ITP-derived parameters.

    python py/martini_build_tables.py [--upside-home <repo>] [--force]

writes parameters/ff_2.1/martini.h5 (the /particles and /sc_table groups). The number of parallel
table-build workers comes from UPSIDE_MARTINI_TABLE_WORKERS, then the Slurm CPU allocation, then the
local CPU count.
"""

import argparse
import math
import os
import importlib.util
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Set, Tuple

import h5py
import numpy as np

from martini_itp_reader import (
    infer_charge_from_atomtype,
    load_martini_forcefield,
    parse_dry_forcefield,
    parse_itp_file,
    sc_target_label,
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
# The SC-particle table is an azimuthal average, so it must sample the azimuth densely enough to
# converge. 24 target directions and one side-chain bead frame left about two directions per cos node
# and one-node ridges of up to 1e9 E_up that gave single-evaluation force spikes of 33,000 E_up/A in
# glpG (findings 4.10, 4.13). At 37 theta nodes the polar bins need 2400 directions; 2400 x 24 against
# 9600 x 24 differs by 0.01 E_up at the 95th percentile.
SC_TABLE_DIRECTION_COUNT = 2400
SC_TABLE_BEAD_FRAME_COUNT = 24
# The runtime evaluates the SC-particle table as a clamped cubic spline (zero end slopes, as the particle
# tables) in the CB-target distance r and the angle theta between the CB vector and the CB-to-target
# direction. Both clamped conditions hold exactly for the physical table: the azimuthal average is even in
# theta about 0 and pi, and the radial grid runs past the reach of every side-chain bead (the largest
# bead-anchor distance plus the 1.2 nm pair cutoff), so its last nodes are exactly zero. In cos(theta) the
# end slopes are not zero and a clamped spline would be wrong there. 37 theta nodes (5 degrees) represent
# the table to 0.017 E_up at the 95th percentile (0.58 E_up at worst, beside steep walls) against a 73-node
# build; 13 nodes left 0.19 E_up and overshoots of 28 E_up.
SC_TABLE_R_STEP_NM = 0.01
SC_TABLE_THETA_COUNT = 37

CANONICAL_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)

# N, CA, C, CB of the reference residue. The rotamer bead centres of sidechain.h5 are given in the
# affine_alignment frame, whose origin is the centroid of N/CA/C, and the runtime places CB in that frame
# (py/martini_prepare_system_lib.py CB_PLACEMENT), so the table's CB anchor is the CB relative to that
# centroid. The CB vector is CB - CA, a difference, so the origin cancels there.
_REFERENCE_N_CA_C_CB_ANG = np.array([
    [-1.19280531, -0.83127186, 0.0],
    [ 0.0,         0.0,        0.0],
    [ 1.25222632, -0.87268266, 0.0],
    [ 0.0,         0.94375626, 1.2068012],
])
CANONICAL_CB_POSITION_ANG = tuple(_REFERENCE_N_CA_C_CB_ANG[3] - _REFERENCE_N_CA_C_CB_ANG[:3].mean(axis=0))
_cb_direction = _REFERENCE_N_CA_C_CB_ANG[3] - _REFERENCE_N_CA_C_CB_ANG[1]
CANONICAL_CB_VECTOR_UNIT = tuple(_cb_direction / np.linalg.norm(_cb_direction))

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
    values_kj_mol: np.ndarray,
    weights: np.ndarray,
    sampled: np.ndarray,
    temperature: float,
) -> np.ndarray:
    """Weighted Boltzmann free energy over the last axis, of the samples a bin was given (`sampled`).
    It runs over the positive-weight samples, and is their weighted mean at temperature 0; a bin whose
    samples all have weight 0 takes their plain mean."""
    values, weights, sampled = np.broadcast_arrays(
        np.asarray(values_kj_mol, dtype=np.float64), np.asarray(weights, dtype=np.float64), sampled
    )
    positive = weights > 0.0
    weights = np.where(positive, weights, 0.0)
    weight_sum = weights.sum(axis=-1)
    norm = np.where(weight_sum > 0.0, weight_sum, 1.0)
    if temperature <= 0.0:
        free = (weights * np.where(positive, values, 0.0)).sum(axis=-1) / norm
    else:
        kbt = float(temperature) * ENERGY_CONVERSION_KJ_PER_EUP
        emin = np.where(positive, values, np.inf).min(axis=-1)
        emin = np.where(np.isfinite(emin), emin, 0.0)
        shifted = np.where(positive, values - emin[..., None], 0.0)
        z = (weights * np.exp(-shifted / kbt)).sum(axis=-1) / norm
        free = emin - kbt * np.log(np.maximum(z, 1e-300))
    plain_mean = np.where(sampled, values, 0.0).sum(axis=-1) / np.maximum(sampled.sum(axis=-1), 1)
    return np.where(weight_sum > 0.0, free, plain_mean)


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
    group.attrs["base_energy_helper"] = "_pair_energy_kj_mol"
    group.attrs["correction_layer"] = correction_layer
    group.attrs["unit_contract"] = "native_dry_martini_nm_kjmol_e_to_upside_runtime_attrs"
    group.attrs["runtime_representation"] = runtime_representation
    group.attrs["nonbonded_cutoff_nm"] = np.float32(DRY_MARTINI_NONBONDED_CUTOFF_NM)
    group.attrs["numerical_distance_guard_nm"] = np.float32(NUMERICAL_DISTANCE_GUARD_NM)


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


def _grid_bin_weights(
    coord_values: np.ndarray,
    grid: List[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Linear binning of each sample's coordinate onto the uniform grid, as (weights, sampled), each
    (n_sample, n_grid). A sample at or beyond an end node goes whole to it; an interior one is split
    between its two nodes, and one on a node still reaches the next node, with weight 0."""
    coord_values = np.asarray(coord_values, dtype=np.float64)
    n_sample, n_grid = coord_values.size, len(grid)
    weights = np.zeros((n_sample, n_grid), dtype=np.float64)
    sampled = np.zeros((n_sample, n_grid), dtype=bool)
    if n_grid == 1:
        weights[:, 0] = 1.0
        sampled[:, 0] = True
        return weights, sampled
    coord = (coord_values - grid[0]) / (grid[1] - grid[0])
    rows = np.arange(n_sample)
    for end, at_end in ((0, coord <= 0.0), (n_grid - 1, coord >= n_grid - 1)):
        weights[rows[at_end], end] = 1.0
        sampled[rows[at_end], end] = True
    inner = (coord > 0.0) & (coord < n_grid - 1)
    lo = np.floor(coord[inner]).astype(np.int64)
    hi_weight = coord[inner] - lo
    weights[rows[inner], lo] = 1.0 - hi_weight
    weights[rows[inner], lo + 1] = hi_weight
    sampled[rows[inner], lo] = True
    sampled[rows[inner], lo + 1] = True
    return weights, sampled


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
    theta_grid: List[float],
    cb_anchor_nm: List[float],
    cb_vector_unit: List[float],
    sidechain_bead_frame_angles: List[float],
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    temperature: float = 0.0,
) -> Dict[str, Any]:
    n_bead = len(sidechain_bead_types)
    rotamer_positions = [np.asarray(pos, dtype=np.float64) for pos in rotamer_bead_positions_nm]
    for irot, pos in enumerate(rotamer_positions):
        if pos.shape != (n_bead, 3):
            raise RuntimeError(
                f"{residue} rotamer {irot} has bead geometry shape {pos.shape}, "
                f"expected ({n_bead}, 3)"
            )
    n_rotamer = len(rotamer_positions)
    n_angle = len(theta_grid)
    n_radial = len(r_values)

    bead_params = []
    for bead_type, bead_charge in zip(sidechain_bead_types, sidechain_bead_charges):
        params = pair_params.get((bead_type, target_type))
        if params is None:
            params = pair_params.get((target_type, bead_type))
        if params is None:
            raise RuntimeError(f"Missing pair params for ({bead_type}, {target_type})")
        charge_product = bead_charge * target_charge if bead_charge and target_charge else 0.0
        bead_params.append((params["sigma_nm"], params["epsilon_kj_mol"], charge_product))

    cb_anchor_arr = np.asarray(cb_anchor_nm, dtype=np.float64)
    cb_vector_arr = np.asarray(cb_vector_unit, dtype=np.float64)
    cb_vector_arr /= max(float(np.linalg.norm(cb_vector_arr)), 1e-12)
    bead_frame_angles = [float(x) for x in sidechain_bead_frame_angles] or [0.0]
    # every rotamer's beads in every bead frame: (n_rotamer, n_frame, n_bead, 3)
    framed_positions = np.stack([
        np.stack([
            _rotate_points_about_axis_np(pos, cb_vector_arr, angle, cb_anchor_arr)
            for angle in bead_frame_angles
        ])
        for pos in rotamer_positions
    ])

    # theta is the angle between the CB vector and the CB-to-target direction: 0 on the side-chain side
    directions = np.asarray(direction_vectors, dtype=np.float64)
    theta = np.arccos(np.clip(directions @ cb_vector_arr, -1.0, 1.0))
    bin_weight, bin_sampled = _grid_bin_weights(theta, theta_grid)
    if not bin_sampled.any(axis=0).all():
        raise RuntimeError(f"Theta bin empty for {residue} target {target_type}")
    # A bin's samples run over (bead frame, target direction), each direction carrying its bin weight
    # in every frame.
    sample_weight = np.tile(bin_weight.T, (1, len(bead_frame_angles)))
    sample_sampled = np.tile(bin_sampled.T, (1, len(bead_frame_angles)))

    rotamer_energy = np.zeros((n_rotamer, n_radial, n_angle), dtype=np.float64)
    for ir, r_nm in enumerate(r_values):
        target_pos_nm = cb_anchor_arr[None, :] + r_nm * directions
        dist_nm = np.linalg.norm(
            framed_positions[:, :, :, None, :] - target_pos_nm[None, None, None, :, :], axis=-1
        )
        energy = sum(
            _pair_energy_kj_mol(dist_nm[:, :, ib], *bead_params[ib], dist_min_nm=dist_min_nm)
            for ib in range(n_bead)
        )
        samples = energy.reshape(n_rotamer, -1)
        rotamer_energy[:, ir, :] = _boltzmann_free_energy_kj_mol(
            samples[:, None, :], sample_weight[None], sample_sampled[None], temperature
        )

    return {
        "residue": residue,
        "target_label": sc_target_label(target_type, target_charge),
        "rotamer_count": n_rotamer,
        "rotamer_probability_fixed": rotamer_weights,
        "rotamer_full_energy_kj_mol": rotamer_energy,
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
    active_targets: List[Tuple[str, float]],
    martini_sidechain_offsets_nm: Dict[str, np.ndarray] | None = None,
    r_step_nm: float = SC_TABLE_R_STEP_NM,
    direction_count: int = SC_TABLE_DIRECTION_COUNT,
    bead_frame_count: int = SC_TABLE_BEAD_FRAME_COUNT,
    theta_count: int = SC_TABLE_THETA_COUNT,
    average_temperature: float = DEFAULT_PRODUCTION_TEMP_UPSIDE,
) -> None:
    """Write /sc_table for every (residue, environment target), a target being a bead type with its
    topology charge (`active_targets` holds (type, charge) pairs)."""
    orientation_map = _load_sidechain_orientation_library(sidechain_lib_path)
    direction_vectors = _fibonacci_sphere(direction_count)
    theta_grid = _linspace(0.0, math.pi, theta_count)
    cb_anchor_nm = [x * ANGSTROM_TO_NM for x in CANONICAL_CB_POSITION_ANG]
    cb_vector_unit = list(CANONICAL_CB_VECTOR_UNIT)
    sidechain_bead_frame_angles = [float(x) for x in _bead_frame_angles(bead_frame_count)]

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

    # Radial grid from the CB anchor out past the reach of every bead: a target farther than the largest
    # bead-anchor distance plus the pair cutoff interacts with no bead, so the last two nodes, and with
    # them the table's value and slope at its end, are exactly zero.
    bead_reach_nm = max(
        float(np.linalg.norm(np.asarray(pos, dtype=np.float64) - np.asarray(cb_anchor_nm), axis=-1).max())
        for _, _, _, _, positions in residue_tasks for pos in positions
    )
    n_r = int(math.ceil((DRY_MARTINI_NONBONDED_CUTOFF_NM + bead_reach_nm) / r_step_nm)) + 2
    r_values = [r_step_nm * i for i in range(n_r)]

    # Collect results
    by_residue: Dict[str, Dict[str, Dict[str, Any]]] = {}
    sc_tasks = []
    for residue, bead_types, bead_charges, orientation, bead_positions_nm in residue_tasks:
        for target_type, target_charge in active_targets:
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
                    "target_charge": float(target_charge),
                    "rotamer_bead_positions_nm": bead_positions_nm,
                    "rotamer_weights": orientation["weight"],
                    "pair_params": pair_params,
                    "r_values": r_values,
                    "direction_vectors": direction_vectors,
                    "theta_grid": theta_grid,
                    "cb_anchor_nm": cb_anchor_nm,
                    "cb_vector_unit": cb_vector_unit,
                    "sidechain_bead_frame_angles": sidechain_bead_frame_angles,
                    "dist_min_nm": NUMERICAL_DISTANCE_GUARD_NM,
                    "temperature": float(average_temperature),
                }
            )

    print(
        f"  sc_table direction-vector sampling: {direction_count} target directions, "
        f"{len(sidechain_bead_frame_angles)} SC bead-frame sample(s), {n_r} radial x {theta_count} theta nodes"
    )
    for result in _parallel_map_ordered("SC-particle table", _run_sc_task_from_dict, sc_tasks):
        by_residue.setdefault(result["residue"], {})[result["target_label"]] = result

    residues = sorted(by_residue)
    targets = [sc_target_label(t, q) for t, q in active_targets]
    n_t = len(targets)
    max_rot = max(int(by_residue[r][targets[0]]["rotamer_count"]) for r in residues)

    rc = np.zeros((len(residues),), dtype=np.float32)
    rpf = np.zeros((len(residues), max_rot), dtype=np.float32)
    r_full = np.zeros((len(residues), max_rot, n_t, n_r, theta_count), dtype=np.float64)
    for ri, residue in enumerate(residues):
        for ti, target in enumerate(targets):
            e = by_residue[residue][target]
            n_rot = int(e["rotamer_count"])
            if ti == 0:
                rc[ri] = float(n_rot)
                rpf[ri, :n_rot] = np.asarray(e["rotamer_probability_fixed"], dtype=np.float32)
            r_full[ri, :n_rot, ti] = e["rotamer_full_energy_kj_mol"]
    tail = np.abs(r_full[:, :, :, -2:, :]).max()
    if tail != 0.0:
        raise RuntimeError(f"sc_table does not vanish past the bead reach: |E| = {tail:g} kJ/mol at the last nodes")

    g = h5.create_group("sc_table")
    _write_common_table_contract_attrs(
        g,
        table_family="SC-particle",
        source_object="sidechain_rotamer_bead_ensemble",
        target_object="dry_martini_particle_type_and_charge",
        projection_ensemble="sidechain_rotamers_x_sidechain_bead_frame_x_target_direction_grid",
        runtime_representation="clamped_cubic_spline_r_theta",
    )
    g.attrs["schema"] = SCHEMA_SC
    g.attrs["forcefield_name"] = "martini22"
    g.attrs["sample_dist_min_nm"] = np.float32(NUMERICAL_DISTANCE_GUARD_NM)
    g.attrs["sidechain_bead_frame_count"] = np.int32(len(sidechain_bead_frame_angles))
    g.attrs["target_direction_count"] = np.int32(direction_count)
    g.attrs["orientation_sampling"] = "target_direction_vector_grid"
    g.attrs["azimuthal_average"] = by_residue[residues[0]][targets[0]]["azimuthal_average"]
    g.attrs["azimuthal_average_temperature_upside"] = np.float32(average_temperature)
    g.attrs["excluded_area_source"] = "none_direct_dry_martini_sc_particle_table"
    g.attrs["attractive_control_source"] = "retained_direct_dry_martini_sc_particle_table"
    g.attrs["isotropic_background_source"] = "none_direct_dry_martini_sc_particle_table"
    g.attrs["relaxation"] = "rigid_rotated_geometry"
    g.attrs["spline_control_quantity"] = "direct_rotamer_free_energy_kj_mol"
    g.attrs["bead_reach_nm"] = np.float32(bead_reach_nm)
    g.create_dataset("restype_order", data=np.asarray([np.bytes_(x) for x in residues], dtype="S4"))
    g.create_dataset("target_order", data=np.asarray([np.bytes_(x) for x in targets], dtype="S8"))
    g.create_dataset("target_type", data=np.asarray([np.bytes_(t) for t, _ in active_targets], dtype="S8"))
    g.create_dataset("target_charge", data=np.asarray([q for _, q in active_targets], dtype=np.float32))
    # Convert native dry-MARTINI units (nm, kJ/mol) to Upside runtime units (Angstrom, E_up)
    # once here at build time, so the runtime h5/config stores Upside-unit values directly.
    g.create_dataset("grid_ang", data=np.asarray(r_values) * LENGTH_CONVERSION_A_PER_NM, dtype=np.float32)
    g.create_dataset("theta_grid", data=np.asarray(theta_grid), dtype=np.float32)
    g.create_dataset("rotamer_count", data=rc, dtype=np.float32)
    g.create_dataset("rotamer_probability_fixed", data=rpf, dtype=np.float32)
    g.create_dataset("rotamer_full_energy_eup", data=r_full / ENERGY_CONVERSION_KJ_PER_EUP, dtype=np.float32)

    print(
        f"  sc_table: {len(residues)} residues x {len(targets)} targets, "
        f"{theta_count} theta x {n_r} radial points, max {max_rot} rotamers"
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


def _pair_energy_kj_mol(
    dist_nm: np.ndarray,
    sigma_nm: float,
    epsilon_kj_mol: float,
    charge_product: float,
    dist_min_nm: float = NUMERICAL_DISTANCE_GUARD_NM,
    cutoff_nm: float = DRY_MARTINI_NONBONDED_CUTOFF_NM,
) -> np.ndarray:
    """Dry-MARTINI pair energy (kJ/mol) at an array of distances (nm).

    The same MARTINI 2 nonbonded contract as the particle-particle grids: LJ potential-shift and
    reaction-field Coulomb, so energy and force both reach the cutoff smoothly. Truncating bare leaves a
    step of k*q1q2/r_c at r_c -- ~3.5 kT for a charged pair at 1.2 nm -- which breaks energy
    conservation every time a pair crosses (findings 75, 79-80).

    `dist_min_nm` is a numerical singularity guard, not a model parameter. Production table builds
    pass `NUMERICAL_DISTANCE_GUARD_NM`.
    """
    dist_nm = np.asarray(dist_nm, dtype=np.float64)
    eff_dist = np.maximum(dist_nm, dist_min_nm)
    sr = sigma_nm / eff_dist
    sr2 = sr * sr
    sr6 = sr2 * sr2 * sr2
    src6 = (sigma_nm / float(cutoff_nm)) ** 6
    energy = 4.0 * epsilon_kj_mol * (sr6 * sr6 - sr6) - 4.0 * epsilon_kj_mol * (src6 * src6 - src6)
    if charge_product:
        kq = COULOMB_K_DRY_KJ_NM * charge_product
        r_c = float(cutoff_nm)
        k_rf = 1.0 / (2.0 * r_c ** 3)
        c_rf = 3.0 / (2.0 * r_c)
        coulomb = kq * (1.0 / eff_dist + k_rf * eff_dist * eff_dist - c_rf)
        energy = energy + np.where(dist_nm > 1e-10, coulomb, 0.0)
    return np.where(dist_nm > float(cutoff_nm), 0.0, energy)


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


def _environment_targets(dry_ff_path: Path) -> List[Tuple[str, float]]:
    """Every (bead type, charge) of the environment molecule topologies beside the force field file: the
    side-chain table is built per target with its own charge, because a bead type does not fix the charge
    (the PC choline bead is Q0 with charge +1, where the type-name convention gives 0)."""
    targets = set()
    for itp in sorted(Path(dry_ff_path).parent.glob("*.itp")):
        if itp.resolve() == Path(dry_ff_path).resolve():
            continue
        for molecule in parse_itp_file(itp)["molecules"].values():
            for atom in molecule["atoms"]:
                targets.add((str(atom["type"]), float(atom["charge"])))
    if not targets:
        raise RuntimeError(f"No environment topologies found beside {dry_ff_path}")
    return sorted(targets)


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
    environment_targets = _environment_targets(dry_ff_path)

    def _writer(h5: h5py.File) -> None:
        h5.attrs["schema"] = "martini_ff_combined"
        h5.attrs["schema_particles"] = SCHEMA_PARTICLES
        h5.attrs["schema_sc"] = SCHEMA_SC
        _build_particles_group(h5, atomtypes, pair_params, set(atomtypes))
        _build_sc_table_group(
            h5, residue_map, pair_params, sidechain_lib_path,
            active_residue_names=list(CANONICAL_RESIDUES),
            active_targets=environment_targets,
            martini_sidechain_offsets_nm=martini_sidechain_offsets_nm,
            average_temperature=DEFAULT_TEMPERED_AVERAGE_TEMP_UPSIDE,
        )

    _write_h5_atomically(output_path, _writer)
    print(f"Built {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build parameters/ff_2.1/martini.h5.")
    parser.add_argument("--upside-home", type=Path, default=Path(__file__).resolve().parent.parent,
                        help="repository root (default: the one holding this script)")
    parser.add_argument("--force", action="store_true", help="rebuild an existing martini.h5")
    args = parser.parse_args()
    root = args.upside_home.expanduser().resolve()
    output = root / "parameters" / "ff_2.1" / "martini.h5"
    if output.exists() and not args.force:
        print(f"{output} exists; pass --force to rebuild it")
    else:
        build_martini_h5(
            output_path=output,
            dry_ff_path=root / "example" / "16.MARTINI" / "dryMARTINI_itp" / "dry_martini_v2.1.itp",
            martinize_path=root / "py" / "martinize.py",
            sidechain_lib_path=root / "parameters" / "ff_2.1" / "sidechain.h5",
        )
