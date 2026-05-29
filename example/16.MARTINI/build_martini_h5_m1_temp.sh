#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ -f "${PROJECT_ROOT}/.venv/bin/activate" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
else
    echo "ERROR: missing Python environment: ${PROJECT_ROOT}/.venv" >&2
    exit 1
fi

source "${PROJECT_ROOT}/source.sh"

set -u

export UPSIDE_HOME="${PROJECT_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
export PATH="${PROJECT_ROOT}/obj:$PATH"

if [ -z "${UPSIDE_MARTINI_TABLE_WORKERS+x}" ]; then
    if command -v sysctl >/dev/null 2>&1; then
        UPSIDE_MARTINI_TABLE_WORKERS="$(sysctl -n hw.logicalcpu)"
    else
        UPSIDE_MARTINI_TABLE_WORKERS="$(getconf _NPROCESSORS_ONLN)"
    fi
    export UPSIDE_MARTINI_TABLE_WORKERS
fi

UPSIDE_MARTINI_FIT_RELAX_STEPS=50
export UPSIDE_MARTINI_FIT_RELAX_STEPS

DOPC_H5="${PROJECT_ROOT}/parameters/dryMARTINI/dopc.h5"

echo "=== Rebuilding only cg_lipid_pair in ${DOPC_H5} ==="
echo "Using ${UPSIDE_MARTINI_TABLE_WORKERS} MARTINI table worker(s), ${UPSIDE_MARTINI_FIT_RELAX_STEPS} fit relax step(s)"

python3 - "$@" << 'PYEOF'
import os, sys, math, json
from pathlib import Path
import h5py, numpy as np

repo_root = Path(os.environ["UPSIDE_HOME"])
sys.path.insert(0, str(repo_root / "py"))

from martini_build_tables import (
    _fit_cg_lipid_quadspline,
    _bead_frame_count,
    _positive_int_env,
    LENGTH_CONVERSION_A_PER_NM,
    ENERGY_CONVERSION_KJ_PER_EUP,
    DRY_MARTINI_NONBONDED_CUTOFF_NM,
    _ensure_cg_bonds_angles,
    _canonicalize_lipid_reference_to_z,
    derive_dopc_cg_params,
    _write_cg_derived_attrs,
)
import martini_build_tables as _mbt
from martini_itp_reader import parse_dry_forcefield, parse_dopc_from_itp

dry_ff_path = repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1.itp"
lipids_itp_path = repo_root / "parameters" / "dryMARTINI" / "dry_martini_v2.1_lipids.itp"
dopc_pdb_path = repo_root / "parameters" / "dryMARTINI" / "DOPC.pdb"
dopc_h5_path = Path(os.environ.get("DOPC_H5", str(repo_root / "parameters" / "dryMARTINI" / "dopc.h5")))

# Parse inputs (same as build_dopc_h5)
_, pair_params = parse_dry_forcefield(dry_ff_path)
from martini_itp_reader import parse_itp_atomtype_masses
atomtype_masses = parse_itp_atomtype_masses(dry_ff_path)

dopc = parse_dopc_from_itp(lipids_itp_path)
bead_types = dopc["bead_types"]
bead_charges = [float(q) for q in dopc["bead_charges"]]

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

ref_nm = _canonicalize_lipid_reference_to_z(ref_bead_positions_nm)
bead_mass_values = [atomtype_masses[bt] for bt in bead_types]

_ensure_cg_bonds_angles(lipids_itp_path)

derived_params = derive_dopc_cg_params(
    ref_bead_positions_nm=ref_nm,
    bead_types=bead_types,
    pair_params=pair_params,
    bead_masses_g_mol=bead_mass_values,
    bonds=_mbt._CURRENT_CG_BONDS,
    energy_conversion_kj_per_eup=ENERGY_CONVERSION_KJ_PER_EUP,
    length_conversion_ang_per_nm=LENGTH_CONVERSION_A_PER_NM,
)
contact_nm = float(derived_params["contact_nm"])

fit_relax_steps = _positive_int_env("UPSIDE_MARTINI_FIT_RELAX_STEPS", 50)
cg_bead_frame_count = _bead_frame_count("CGL", 1)

_cg_r, _cg_ct, _cg_az = 16, 7, 2
r_min_nm = 0.50
r_max_nm = 1.68

print(f"  contact_nm={contact_nm:.4f}, fit_relax_steps={fit_relax_steps}")
print("  Rebuilding CG-CG full tensor table...")

result_cg = _fit_cg_lipid_quadspline(
    ref_bead_positions_nm=ref_nm,
    bead_types=bead_types,
    bead_charges=bead_charges,
    pair_params=pair_params,
    r_min_nm=r_min_nm,
    r_max_nm=r_max_nm,
    r_count=min(24, _cg_r),
    cos_theta_count=min(13, _cg_ct),
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
    plane_constraint=False,
)

print(
    f"  CG-CG result: {result_cg['n_radial']}r x {result_cg['n_angular']}^2 angular, "
    f"max|E| = {float(np.max(np.abs(result_cg['energy_grid_raw']))):.3f} kJ/mol, "
    f"excluded_area_nonnegative_rows = {result_cg['excluded_area_nonnegative_rows']}, "
    f"attractive_control_count = {result_cg['attractive_control_count']}"
)

# Replace cg_lipid_pair group in-place
print(f"  Writing to {dopc_h5_path} ...")
with h5py.File(dopc_h5_path, "r+") as h5:
    cg_grp = h5["cg_lipid_table"]
    if "cg_lipid_pair" in cg_grp:
        del cg_grp["cg_lipid_pair"]

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

print("  Done. cg_lipid_pair rebuilt successfully.")
PYEOF
