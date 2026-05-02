#!/usr/bin/env python3
"""Build parameters/dryMARTINI/particles.h5 from the dry-MARTINI forcefield.

Pre-computes LJ energy grids per unique (epsilon, sigma) pair and a Coulomb
reference grid (qq = 1.0 e^2) on a 1000-point uniform radial grid [0, 12] A.
The C++ MartiniPotential constructor loads these grids, fits cubic splines,
and assigns them to PairParams — removing the last direct force-field formula
evaluation from the simulation code.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np

COULOMB_K_DRY_KJ_NM = 138.935458 / 15.0
ENERGY_CONVERSION_KJ_PER_EUP = 2.914952774272
LENGTH_CONVERSION_A_PER_NM = 10.0

GRID_N_POINTS = 1000
GRID_R_MIN_A = 0.0
GRID_R_MAX_A = 12.0

POSITIVE_TYPES = {"Qda", "Qd", "SQda", "SQd"}
NEGATIVE_TYPES = {"Qa", "SQa"}

SCHEMA = "martini_particles_v1"


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

            type_i = parts[0]
            type_j = parts[1]
            if len(parts) == 4:
                macro = parts[3]
                if macro not in macros:
                    raise RuntimeError(
                        f"Unknown dry-MARTINI macro '{macro}' in {ff_path}"
                    )
                sigma_nm, epsilon_kj = macros[macro]
            else:
                sigma_nm = float(parts[3])
                epsilon_kj = float(parts[4])

            payload = {"sigma_nm": sigma_nm, "epsilon_kj_mol": epsilon_kj}
            pair_params[(type_i, type_j)] = payload
            pair_params[(type_j, type_i)] = payload

    return atomtypes, pair_params


def _infer_type_charge(bead_type: str) -> float:
    bead_type = bead_type.strip()
    if bead_type in POSITIVE_TYPES:
        return 1.0
    if bead_type in NEGATIVE_TYPES:
        return -1.0
    return 0.0


def compute_lj_grid(
    sigma_nm: float,
    epsilon_kj_mol: float,
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
) -> np.ndarray:
    sig_a = sigma_nm * length_conv
    eps_eup = epsilon_kj_mol / energy_conv
    grid = np.zeros(GRID_N_POINTS, dtype=np.float64)
    for i in range(GRID_N_POINTS):
        r = GRID_R_MIN_A + i * (GRID_R_MAX_A - GRID_R_MIN_A) / (GRID_N_POINTS - 1)
        if r < 0.1 * sig_a:
            r = 0.1 * sig_a
        r2 = r * r
        r3 = r2 * r
        r6 = r3 * r3
        sig2 = sig_a * sig_a
        sig6 = sig2 * sig2 * sig2
        sig12 = sig6 * sig6
        inv_r6 = sig6 / r6
        inv_r12 = sig12 / (r6 * r6)
        grid[i] = 4.0 * eps_eup * (inv_r12 - inv_r6)
    return grid


def compute_coulomb_ref_grid(coulomb_k_eup: float) -> np.ndarray:
    grid = np.zeros(GRID_N_POINTS, dtype=np.float64)
    for i in range(GRID_N_POINTS):
        r = GRID_R_MIN_A + i * (GRID_R_MAX_A - GRID_R_MIN_A) / (GRID_N_POINTS - 1)
        if r == 0.0:
            r = 1.0e-6
        grid[i] = coulomb_k_eup / r
    return grid


def _extract_unique_lj_params(
    pair_params: Dict[Tuple[str, str], Dict[str, float]],
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
) -> List[Tuple[float, float]]:
    seen: set = set()
    result: List[Tuple[float, float]] = []
    for (_ti, _tj), p in sorted(pair_params.items()):
        eps = p["epsilon_kj_mol"] / energy_conv
        sig = p["sigma_nm"] * length_conv
        key = (eps, sig)
        if key not in seen:
            seen.add(key)
            result.append(key)
    return result


def build_particles_h5(
    output_path: Path,
    dry_ff_path: Path,
    energy_conv: float = ENERGY_CONVERSION_KJ_PER_EUP,
    length_conv: float = LENGTH_CONVERSION_A_PER_NM,
) -> None:
    atomtypes, pair_params = _parse_dry_forcefield(dry_ff_path)

    unique_lj = _extract_unique_lj_params(pair_params, energy_conv, length_conv)
    n_lj = len(unique_lj)

    coulomb_k_eup = COULOMB_K_DRY_KJ_NM * length_conv / energy_conv

    eps_arr = np.asarray([p[0] for p in unique_lj], dtype=np.float64)
    sig_arr = np.asarray([p[1] for p in unique_lj], dtype=np.float64)
    lj_grids = np.zeros((n_lj, GRID_N_POINTS), dtype=np.float64)
    for idx, (eps_eup, sig_a) in enumerate(unique_lj):
        sigma_nm = sig_a / length_conv
        epsilon_kj = eps_eup * energy_conv
        lj_grids[idx, :] = compute_lj_grid(sigma_nm, epsilon_kj, energy_conv, length_conv)

    coul_ref = compute_coulomb_ref_grid(coulomb_k_eup)

    charges = np.asarray([_infer_type_charge(t) for t in atomtypes], dtype=np.float32)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as h5:
        h5.attrs["schema"] = SCHEMA
        h5.attrs["lj_n_unique"] = n_lj
        h5.attrs["lj_n_points"] = GRID_N_POINTS
        h5.attrs["lj_r_min_ang"] = np.float32(GRID_R_MIN_A)
        h5.attrs["lj_r_max_ang"] = np.float32(GRID_R_MAX_A)
        h5.attrs["coul_n_points"] = GRID_N_POINTS
        h5.attrs["coul_r_min_ang"] = np.float32(GRID_R_MIN_A)
        h5.attrs["coul_r_max_ang"] = np.float32(GRID_R_MAX_A)
        h5.attrs["coulomb_k_eup"] = np.float32(coulomb_k_eup)

        h5.create_dataset("lj_unique_eps_eup", data=eps_arr)
        h5.create_dataset("lj_unique_sig_ang", data=sig_arr)
        h5.create_dataset("lj_energy_grids", data=lj_grids)
        h5.create_dataset("coulomb_reference_grid", data=coul_ref)
        h5.create_dataset(
            "type_order",
            data=np.asarray([np.bytes_(x) for x in atomtypes], dtype="S8"),
        )
        h5.create_dataset("type_charge", data=charges)

    print(
        f"Built {output_path}: {n_lj} unique LJ pairs, "
        f"{len(atomtypes)} types, {GRID_N_POINTS} radial points"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build particles.h5 from dry-MARTINI forcefield"
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build", help="Build particles.h5")
    p_build.add_argument(
        "--dry-forcefield",
        required=True,
        type=Path,
        help="Path to dry_martini_v2.1.itp",
    )
    p_build.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: parameters/dryMARTINI/particles.h5)",
    )
    p_build.add_argument(
        "--energy-conversion",
        type=float,
        default=ENERGY_CONVERSION_KJ_PER_EUP,
        help="kJ/mol per EUP",
    )
    p_build.add_argument(
        "--length-conversion",
        type=float,
        default=LENGTH_CONVERSION_A_PER_NM,
        help="Angstrom per nm",
    )

    args = parser.parse_args()

    if args.cmd == "build":
        output = args.output
        if output is None:
            script_dir = Path(__file__).resolve().parent
            output = script_dir.parent / "parameters" / "dryMARTINI" / "particles.h5"
        else:
            output = output.expanduser().resolve()

        build_particles_h5(
            output_path=output,
            dry_ff_path=Path(args.dry_forcefield).expanduser().resolve(),
            energy_conv=args.energy_conversion,
            length_conv=args.length_conversion,
        )
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
