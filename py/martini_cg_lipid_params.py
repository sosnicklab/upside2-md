#!/usr/bin/env python3
"""Shared dry-MARTINI DOPC-derived parameters for single-particle CG lipids."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np

ENERGY_CONVERSION_KJ_PER_EUP = 2.914952774272
LENGTH_CONVERSION_A_PER_NM = 10.0

DOPC_ATOM_NAMES = (
    "NC3", "PO4", "GL1", "GL2",
    "C1A", "C2A", "D3A", "C4A", "C5A",
    "C1B", "C2B", "D3B", "C4B", "C5B",
)

DOPC_BEAD_TYPES = (
    "Q0", "Qa", "Na", "Na",
    "C1", "C1", "C3", "C1", "C1",
    "C1", "C1", "C3", "C1", "C1",
)

DOPC_BEAD_CHARGES = (
    1.0, -1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
)

# DOPC bonded topology from dry_martini_v2.1_lipids.itp (0-based indices).
DOPC_BONDS = (
    (0, 1, 0.450, 1250.0),
    (1, 2, 0.450, 1250.0),
    (2, 3, 0.370, 1250.0),
    (2, 4, 0.480, 1250.0),
    (4, 5, 0.480, 1250.0),
    (5, 6, 0.480, 1250.0),
    (6, 7, 0.480, 1250.0),
    (7, 8, 0.480, 1250.0),
    (3, 9, 0.480, 1250.0),
    (9, 10, 0.480, 1250.0),
    (10, 11, 0.480, 1250.0),
    (11, 12, 0.480, 1250.0),
    (12, 13, 0.480, 1250.0),
)


def _pair_param(pair_params: dict, type_i: str, type_j: str) -> dict | None:
    return pair_params.get((type_i, type_j)) or pair_params.get((type_j, type_i))


def dopc_max_sigma_nm(bead_types: Iterable[str], pair_params: dict) -> float:
    sigmas = []
    types = list(bead_types)
    for ti in types:
        for tj in types:
            params = _pair_param(pair_params, ti, tj)
            if params is not None:
                sigmas.append(float(params["sigma_nm"]))
    if not sigmas:
        raise ValueError("No dry-MARTINI nonbonded sigmas found for DOPC bead types")
    return max(sigmas)


def parse_dry_martini_masses(ff_file: str | Path) -> dict[str, float]:
    """Read dry-MARTINI atomtype masses from an ITP file."""
    ff_path = Path(ff_file).expanduser().resolve()
    masses: dict[str, float] = {}
    in_atomtypes = False
    with ff_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("["):
                in_atomtypes = line.replace(" ", "") == "[atomtypes]"
                continue
            if not in_atomtypes:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                masses[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return masses


def derive_dopc_cg_params(
    ref_bead_positions_nm: np.ndarray,
    bead_types: Iterable[str],
    pair_params: dict,
    bead_masses_g_mol: Iterable[float] | None = None,
    bonds: Iterable[tuple[int, int, float, float]] = DOPC_BONDS,
    energy_conversion_kj_per_eup: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conversion_ang_per_nm: float = LENGTH_CONVERSION_A_PER_NM,
) -> dict:
    """Derive CGL geometry and contact parameters from a DOPC reference lipid.

    The returned values are in Upside stage units where noted.  This keeps the
    single-particle lipid model tied to the dry-MARTINI bead geometry and ITP
    parameters instead of standalone stabilization constants.
    """
    ref_nm = np.asarray(ref_bead_positions_nm, dtype=np.float64)
    if ref_nm.shape != (14, 3):
        raise ValueError(f"ref_bead_positions_nm must have shape (14,3), got {ref_nm.shape}")

    bead_types = list(bead_types)
    if len(bead_types) != 14:
        raise ValueError(f"DOPC bead type list must have 14 entries, got {len(bead_types)}")

    if energy_conversion_kj_per_eup <= 0.0:
        raise ValueError("energy_conversion_kj_per_eup must be positive")
    if length_conversion_ang_per_nm <= 0.0:
        raise ValueError("length_conversion_ang_per_nm must be positive")

    ref_ang = ref_nm * length_conversion_ang_per_nm
    com_ang = np.mean(ref_ang, axis=0)
    rel_ang = ref_ang - com_ang
    head_ang = ref_ang[0]
    tail_mid_ang = 0.5 * (ref_ang[8] + ref_ang[13])
    axis = tail_mid_ang - head_ang
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1.0e-8:
        raise ValueError("Cannot derive DOPC orientation axis from coincident head/tail sites")
    unit_axis = axis / axis_norm

    projection_ang = rel_ang @ unit_axis
    r2_perp = np.sum(rel_ang * rel_ang, axis=1) - projection_ang * projection_ang
    r2_perp = np.maximum(r2_perp, 0.0)

    if bead_masses_g_mol is None:
        masses = np.full(14, 72.0, dtype=np.float64)
        mass_source = "uniform_72_g_mol_fallback"
    else:
        masses = np.asarray(list(bead_masses_g_mol), dtype=np.float64)
        if masses.shape != (14,) or not np.all(masses > 0.0):
            raise ValueError("DOPC bead masses must be 14 positive values")
        mass_source = "dry_martini_atomtypes"

    orientation_length_ang = abs(float(np.dot(tail_mid_ang - com_ang, unit_axis)))
    if orientation_length_ang <= 1.0e-6:
        orientation_length_ang = 0.5 * axis_norm
    transverse_inertia = float(np.sum(masses * r2_perp))
    orientation_mass_g = transverse_inertia / max(orientation_length_ang * orientation_length_ang, 1.0e-12)

    projected_k = 0.0
    for i, j, _r0_nm, k_kj_mol_nm2 in bonds:
        dr_nm = ref_nm[int(i)] - ref_nm[int(j)]
        norm_nm = float(np.linalg.norm(dr_nm))
        if norm_nm <= 1.0e-12:
            continue
        cos_axis = float(abs(np.dot(dr_nm / norm_nm, unit_axis)))
        projected_k += float(k_kj_mol_nm2) * cos_axis * cos_axis
    if projected_k <= 0.0:
        projected_k = float(np.median([b[3] for b in bonds]))
    bond_fc_eup_a2 = projected_k / (
        energy_conversion_kj_per_eup
        * length_conversion_ang_per_nm
        * length_conversion_ang_per_nm
    )

    max_sigma_nm = dopc_max_sigma_nm(bead_types, pair_params)
    contact_nm = (2.0 ** (1.0 / 6.0)) * max_sigma_nm

    return {
        "contact_nm": float(contact_nm),
        "contact_ang": float(contact_nm * length_conversion_ang_per_nm),
        "max_sigma_nm": float(max_sigma_nm),
        "orientation_length_ang": float(orientation_length_ang),
        "orientation_mass_g_mol": float(orientation_mass_g),
        "orientation_bond_fc_eup_a2": float(bond_fc_eup_a2),
        "transverse_inertia_g_mol_a2": float(transverse_inertia),
        "head_tail_span_ang": float(axis_norm),
        "tail_projection_ang": float(orientation_length_ang),
        "max_perp_radius_ang": float(np.sqrt(np.max(r2_perp))),
        "mass_source": mass_source,
        "contact_source": "2^(1/6)*max_dopc_pair_sigma",
        "orientation_length_source": "DOPC_COM_to_tail_midpoint_projection",
        "orientation_mass_source": "DOPC_transverse_rotational_inertia",
        "orientation_bond_fc_source": "sum_DOPC_bond_k_projected_on_orientation_axis",
        "energy_conversion_kj_per_eup": float(energy_conversion_kj_per_eup),
        "length_conversion_ang_per_nm": float(length_conversion_ang_per_nm),
    }
