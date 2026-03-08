#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../source.sh"
source "${SCRIPT_DIR}/../../.venv/bin/activate"
set -euo pipefail
cd "${SCRIPT_DIR}"

# Hybrid 1RKL workflow:
# 0) Hybrid preparation (packed MARTINI + hybrid mapping export)
# 1) Stage input generation (dry MARTINI)
# 2) Inject hybrid mapping + stage-control metadata into each stage .up file
# 2.5) Validate the published MARTINI force-field artifact and inject the production cross term
# 3) Run 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0
# 4) Extract key VTF trajectories

# =============================================================================
# USER CONFIGURATION
# =============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        PDB_ID=*)
            PDB_ID="${1#*=}"
            shift
            ;;
        *)
            echo "Unknown parameter $1"
            exit 1
            ;;
    esac
done

PDB_ID="${PDB_ID:-1rkl}"
RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_hybrid}"

INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_hybrid}"
CHECKPOINT_DIR="${CHECKPOINT_DIR:-${RUN_DIR}/checkpoints}"
LOG_DIR="${LOG_DIR:-${RUN_DIR}/logs}"
HYBRID_PREP_DIR="${HYBRID_PREP_DIR:-${RUN_DIR}/hybrid_prep}"

PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
BILAYER_PDB="${BILAYER_PDB:-pdb/bilayer.MARTINI.pdb}"
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${SCRIPT_DIR}/prepare_system.py}"
UNIVERSAL_PREP_MODE="${UNIVERSAL_PREP_MODE:-both}"

RUNTIME_PDB_FILE="${RUNTIME_PDB_FILE:-${HYBRID_PREP_DIR}/${RUNTIME_PDB_ID}.MARTINI.pdb}"
HYBRID_MAPPING_FILE="${HYBRID_MAPPING_FILE:-${HYBRID_PREP_DIR}/hybrid_mapping.h5}"
HYBRID_PACKED_PDB="${HYBRID_PACKED_PDB:-${HYBRID_PREP_DIR}/hybrid_packed.MARTINI.pdb}"
HYBRID_VALIDATE="${HYBRID_VALIDATE:-1}"

PREPARED_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.prepared.up"
STAGE_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
PREPARED_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.prepared.up"
STAGE_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
STAGE_62_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.2.up"
PREPARED_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.prepared.up"
STAGE_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.up"
PREPARED_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.prepared.up"
STAGE_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.up"
PREPARED_65_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.prepared.up"
STAGE_65_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.up"
PREPARED_66_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.prepared.up"
STAGE_66_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up"
PREPARED_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.prepared.up"
STAGE_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up"

BACKBONE_CROSS_ENABLE="${BACKBONE_CROSS_ENABLE:-1}"
BACKBONE_CROSS_REBUILD="${BACKBONE_CROSS_REBUILD:-0}"
DEPTH_TABLE_CSV="${DEPTH_TABLE_CSV:-${RUN_DIR}/depth_interaction_table.csv}"
DEPTH_TABLE_META="${DEPTH_TABLE_META:-${RUN_DIR}/depth_interaction_table.meta.json}"
BACKBONE_CROSS_TABLE_CSV="${BACKBONE_CROSS_TABLE_CSV:-${RUN_DIR}/backbone_cross_interaction_table.csv}"
BACKBONE_CROSS_TABLE_META="${BACKBONE_CROSS_TABLE_META:-${RUN_DIR}/backbone_cross_interaction_table.meta.json}"
BACKBONE_CROSS_TABLE_H5="${BACKBONE_CROSS_TABLE_H5:-${UPSIDE_HOME}/parameters/ff_2.1/martini.h5}"
export UPSIDE_MARTINI_FORCEFIELD_H5="${UPSIDE_MARTINI_FORCEFIELD_H5:-${BACKBONE_CROSS_TABLE_H5}}"
BACKBONE_CROSS_BUILD_SCRIPT="${BACKBONE_CROSS_BUILD_SCRIPT:-${SCRIPT_DIR}/build_martini_parameters.py}"

# Stage-0 hybrid structure prep controls
SALT_MOLAR="${SALT_MOLAR:-0.15}"
# More conservative default exclusion to prevent lipids packed too close to protein surface.
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-4.5}"
ION_CUTOFF="${ION_CUTOFF:-4.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PREP_SEED="${PREP_SEED:-2026}"
PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-4.5}"
PROTEIN_LIPID_CUTOFF_STEP="${PROTEIN_LIPID_CUTOFF_STEP:-0.5}"
PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-8.0}"
BB_AA_MIN_MATCHED_RESIDUES="${BB_AA_MIN_MATCHED_RESIDUES:-8}"
BB_AA_MAX_RIGID_RMSD="${BB_AA_MAX_RIGID_RMSD:-1.5}"

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"
STRICT_STAGE_HANDOFF="${STRICT_STAGE_HANDOFF:-1}"
HYBRID_CONTROL_ENABLE="${HYBRID_CONTROL_ENABLE:-1}"
HYBRID_ACTIVATION_STAGE="${HYBRID_ACTIVATION_STAGE:-production}"
HYBRID_PREPROD_PROTEIN_MODE="${HYBRID_PREPROD_PROTEIN_MODE:-rigid}"
HYBRID_PREP_RUNTIME_MODE="${HYBRID_PREP_RUNTIME_MODE:-dry_martini_prep}"
HYBRID_ACTIVE_RUNTIME_MODE="${HYBRID_ACTIVE_RUNTIME_MODE:-aa_backbone_explicit_lipid}"
HYBRID_PREPROD_LIPID_HEADGROUP_ROLES="${HYBRID_PREPROD_LIPID_HEADGROUP_ROLES:-PO4}"
PREPROD_FIX_RIGID_ENABLE="${PREPROD_FIX_RIGID_ENABLE:-1}"
PREPROD_RIGID_ASSERT_ENABLE="${PREPROD_RIGID_ASSERT_ENABLE:-1}"
PREPROD_RIGID_TOL="${PREPROD_RIGID_TOL:-1e-6}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-5000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
# Production runs with injected Upside all-atom backbone nodes; use a smaller
# stable timestep than MARTINI-only stages.
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-5000}"

COMP_3E4_BAR_INV_TO_A3_PER_EUP="${COMP_3E4_BAR_INV_TO_A3_PER_EUP:-14.521180763676}"
# Unit conversion from AGENTS.md:
#   1 bar = 0.000020659477 E_up / Angstrom^3
BAR_1_TO_EUP_PER_A3="${BAR_1_TO_EUP_PER_A3:-0.000020659477}"
NPT_REF_P_ZERO="${NPT_REF_P_ZERO:-0.0}"
NPT_REF_P_ONE_BAR="${NPT_REF_P_ONE_BAR:-$BAR_1_TO_EUP_PER_A3}"
export UPSIDE_NPT_TARGET_PXY="${UPSIDE_NPT_TARGET_PXY:-$NPT_REF_P_ZERO}"
export UPSIDE_NPT_TARGET_PZ="${UPSIDE_NPT_TARGET_PZ:-$NPT_REF_P_ZERO}"
export UPSIDE_NPT_TAU="${UPSIDE_NPT_TAU:-4.0}"
export UPSIDE_NPT_COMPRESSIBILITY="${UPSIDE_NPT_COMPRESSIBILITY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_XY="${UPSIDE_NPT_COMPRESSIBILITY_XY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_Z="${UPSIDE_NPT_COMPRESSIBILITY_Z:-0.0}"
export UPSIDE_NPT_INTERVAL="${UPSIDE_NPT_INTERVAL:-10}"
export UPSIDE_NPT_SEMI="${UPSIDE_NPT_SEMI:-1}"
export UPSIDE_NPT_DEBUG="${UPSIDE_NPT_DEBUG:-1}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"
PROD_70_BAROSTAT_TYPE="${PROD_70_BAROSTAT_TYPE:-1}"

export UPSIDE_OVERWRITE_SPLINES="${UPSIDE_OVERWRITE_SPLINES:-1}"
export UPSIDE_EWALD_ENABLE="${UPSIDE_EWALD_ENABLE:-1}"
export UPSIDE_EWALD_ALPHA="${UPSIDE_EWALD_ALPHA:-0.2}"
export UPSIDE_EWALD_KMAX="${UPSIDE_EWALD_KMAX:-5}"
export UPSIDE_MARTINI_FF_DIR="${UPSIDE_MARTINI_FF_DIR:-ff_dry}"
UPSIDE_RAMA_LIBRARY="${UPSIDE_RAMA_LIBRARY:-${UPSIDE_HOME}/parameters/common/rama.dat}"
UPSIDE_RAMA_SHEET_MIXING="${UPSIDE_RAMA_SHEET_MIXING:-${UPSIDE_HOME}/parameters/ff_2.1/sheet}"
UPSIDE_HBOND_ENERGY="${UPSIDE_HBOND_ENERGY:-${UPSIDE_HOME}/parameters/ff_2.1/hbond.h5}"
UPSIDE_REFERENCE_STATE_RAMA="${UPSIDE_REFERENCE_STATE_RAMA:-${UPSIDE_HOME}/parameters/common/rama_reference.pkl}"

# =============================================================================
# VALIDATION
# =============================================================================
if [ -z "${UPSIDE_HOME:-}" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi

if [ ! -f "${PROTEIN_AA_PDB}" ]; then
    echo "ERROR: protein AA PDB not found: ${PROTEIN_AA_PDB}"
    exit 1
fi
if [ ! -f "${BILAYER_PDB}" ]; then
    echo "ERROR: bilayer PDB not found: ${BILAYER_PDB}"
    exit 1
fi
if [ ! -f "${UPSIDE_RAMA_LIBRARY}" ]; then
    echo "ERROR: Upside rama library not found: ${UPSIDE_RAMA_LIBRARY}"
    exit 1
fi
if [ ! -f "${UPSIDE_HBOND_ENERGY}" ]; then
    echo "ERROR: Upside hbond energy file not found: ${UPSIDE_HBOND_ENERGY}"
    exit 1
fi
if [ ! -f "${UPSIDE_REFERENCE_STATE_RAMA}" ]; then
    echo "ERROR: Upside reference-state rama file not found: ${UPSIDE_REFERENCE_STATE_RAMA}"
    exit 1
fi
if [ ! -f "${UNIVERSAL_PREP_SCRIPT}" ]; then
    echo "ERROR: universal prep script not found: ${UNIVERSAL_PREP_SCRIPT}"
    exit 1
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$HYBRID_PREP_DIR"

set_barostat_type() {
    local up_file="$1"
    local barostat_type="$2"
    python3 - "$up_file" "$barostat_type" << 'PY'
import sys
import tables as tb
up_file = sys.argv[1]
barostat_type = int(sys.argv[2])
with tb.open_file(up_file, 'r+') as t:
    if '/input/barostat' in t:
        t.root.input.barostat._v_attrs.type = barostat_type
PY
}

set_stage_label() {
    local up_file="$1"
    local stage_label="$2"
    python3 - "$up_file" "$stage_label" << 'PY'
import sys
import h5py
import numpy as np
up_file = sys.argv[1]
stage_label = sys.argv[2]
with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    grp = inp.require_group("stage_parameters")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["current_stage"] = np.bytes_(stage_label)
PY
}

inject_hybrid_mapping() {
    local up_file="$1"
    local mapping_file="$2"
    python3 - "$up_file" "$mapping_file" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
mapping_file = sys.argv[2]
groups = [
    "hybrid_bb_map",
    "hybrid_env_topology",
]
optional_groups = [
    "hybrid_control",
]

BB_COMPONENT_MASS = {
    "N": 14.0,
    "CA": 12.0,
    "C": 12.0,
    "O": 16.0,
}
BB_COMPONENT_SUM = float(sum(BB_COMPONENT_MASS.values()))
DEFAULT_COMPONENT_NAMES = ("N", "CA", "C", "O")

def as_text(x):
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", errors="ignore")
    return str(x)

def pad_bytes(value, dtype):
    n = np.dtype(dtype).itemsize
    raw = str(value).encode("ascii", errors="ignore")[:n]
    return raw.ljust(n, b" ")

def replace_dataset(parent, name, data):
    if name in parent:
        del parent[name]
    parent.create_dataset(name, data=data)

with h5py.File(mapping_file, "r") as src, h5py.File(up_file, "r+") as dst:
    src_inp = src["/input"]
    dst_inp = dst.require_group("input")

    base_n_atom = int(dst["/input/pos"].shape[0])
    src_mem = src["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    mem_len = int(src_mem.shape[0])
    if base_n_atom != mem_len:
        raise ValueError(
            f"Hybrid mapping n_atom mismatch for {up_file}: up has {base_n_atom}, mapping has {mem_len}"
        )

    for g in groups:
        if g not in src_inp:
            raise ValueError(f"Missing mapping group in {mapping_file}: /input/{g}")
        if g in dst_inp:
            del dst_inp[g]
        src.copy(src_inp[g], dst_inp, name=g)

    for g in optional_groups:
        if g not in src_inp:
            continue
        if g in dst_inp:
            del dst_inp[g]
        src.copy(src_inp[g], dst_inp, name=g)

    bb_grp = dst_inp["hybrid_bb_map"]
    env_grp = dst_inp["hybrid_env_topology"]
    pot_grp = dst_inp["potential"]["martini_potential"]

    bb_residue = bb_grp["bb_residue_index"][:].astype(np.int32)
    bb_ref_idx = (
        bb_grp["reference_atom_indices"][:].astype(np.int32)
        if "reference_atom_indices" in bb_grp
        else bb_grp["atom_indices"][:].astype(np.int32)
    )
    n_bb = int(bb_ref_idx.shape[0])
    if bb_ref_idx.ndim != 2 or bb_ref_idx.shape[1] != 4:
        raise ValueError("hybrid_bb_map reference/active indices must have shape (n_bb,4)")

    if "reference_atom_coords" in bb_grp:
        bb_ref_xyz = bb_grp["reference_atom_coords"][:].astype(np.float32)
        if bb_ref_xyz.shape != (n_bb, 4, 3):
            raise ValueError("hybrid_bb_map/reference_atom_coords must have shape (n_bb,4,3)")
    else:
        bb_ref_xyz = np.zeros((n_bb, 4, 3), dtype=np.float32)

    comp_names = list(DEFAULT_COMPONENT_NAMES)
    if "reference_atom_names" in bb_grp:
        names_raw = bb_grp["reference_atom_names"][:]
        if names_raw.shape == (4,):
            comp_names = [as_text(x).strip().upper() or DEFAULT_COMPONENT_NAMES[i] for i, x in enumerate(names_raw)]

    valid_ref = bb_ref_idx >= 0
    max_ref_idx = int(np.max(bb_ref_idx[valid_ref])) if np.any(valid_ref) else -1
    n_ref_index = max_ref_idx + 1
    ref_offset = base_n_atom

    if n_ref_index > 0:
        pos_dtype = dst_inp["pos"].dtype
        mom_dtype = dst_inp["mom"].dtype
        vel_dtype = dst_inp["vel"].dtype
        mass_dtype = dst_inp["mass"].dtype
        charge_dtype = dst_inp["charges"].dtype
        type_dtype = dst_inp["type"].dtype
        name_dtype = dst_inp["atom_names"].dtype if "atom_names" in dst_inp else type_dtype
        role_dtype = dst_inp["atom_roles"].dtype if "atom_roles" in dst_inp else type_dtype
        residue_dtype = dst_inp["residue_ids"].dtype
        molecule_dtype = dst_inp["molecule_ids"].dtype
        pclass_dtype = dst_inp["particle_class"].dtype

        ref_pos = np.zeros((n_ref_index, 3, 1), dtype=pos_dtype)
        ref_mom = np.zeros((n_ref_index, 3, 1), dtype=mom_dtype)
        ref_vel = np.zeros((n_ref_index, 3), dtype=vel_dtype)
        ref_mass = np.ones((n_ref_index,), dtype=mass_dtype)
        ref_charge = np.zeros((n_ref_index,), dtype=charge_dtype)
        ref_type = np.full((n_ref_index,), pad_bytes("AA", type_dtype), dtype=type_dtype)
        ref_name = np.full((n_ref_index,), pad_bytes("AA", name_dtype), dtype=name_dtype)
        ref_role = np.full((n_ref_index,), pad_bytes("AA", role_dtype), dtype=role_dtype)
        ref_residue = np.full((n_ref_index,), -1, dtype=residue_dtype)
        ref_molecule = np.full((n_ref_index,), 0, dtype=molecule_dtype)
        ref_pclass = np.full((n_ref_index,), pad_bytes("PROTEINAA", pclass_dtype), dtype=pclass_dtype)

        for k in range(n_bb):
            resid = int(bb_residue[k])
            for d in range(4):
                ref_idx = int(bb_ref_idx[k, d])
                if ref_idx < 0:
                    continue
                ref_pos[ref_idx, :, 0] = bb_ref_xyz[k, d, :]
                cname = comp_names[d] if d < len(comp_names) else DEFAULT_COMPONENT_NAMES[d]
                ref_name[ref_idx] = pad_bytes(cname, name_dtype)
                ref_role[ref_idx] = pad_bytes(cname, role_dtype)
                ref_residue[ref_idx] = resid
                comp_mass = float(BB_COMPONENT_MASS.get(cname, 12.0))
                ref_mass[ref_idx] = mass_dtype.type(comp_mass / BB_COMPONENT_SUM)

        def append_input_dataset(name, appendix):
            if name not in dst_inp:
                return
            base = dst_inp[name][:]
            if base.shape[0] != base_n_atom:
                raise ValueError(f"{up_file}: /input/{name} length changed unexpectedly before augmentation")
            merged = np.concatenate([base, appendix.astype(base.dtype, copy=False)], axis=0)
            replace_dataset(dst_inp, name, merged)

        append_input_dataset("pos", ref_pos)
        append_input_dataset("mom", ref_mom)
        append_input_dataset("vel", ref_vel)
        append_input_dataset("mass", ref_mass)
        append_input_dataset("charges", ref_charge)
        append_input_dataset("type", ref_type)
        append_input_dataset("atom_names", ref_name)
        append_input_dataset("atom_roles", ref_role)
        append_input_dataset("residue_ids", ref_residue)
        append_input_dataset("molecule_ids", ref_molecule)
        append_input_dataset("particle_class", ref_pclass)

        pot_atom_indices = pot_grp["atom_indices"][:]
        if pot_atom_indices.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: /input/potential/martini_potential/atom_indices length mismatch")
        pot_ai_append = np.arange(ref_offset, ref_offset + n_ref_index, dtype=pot_atom_indices.dtype)
        replace_dataset(pot_grp, "atom_indices", np.concatenate([pot_atom_indices, pot_ai_append], axis=0))

        pot_charges = pot_grp["charges"][:]
        if pot_charges.shape[0] != base_n_atom:
            raise ValueError(f"{up_file}: /input/potential/martini_potential/charges length mismatch")
        pot_q_append = np.zeros((n_ref_index,), dtype=pot_charges.dtype)
        replace_dataset(pot_grp, "charges", np.concatenate([pot_charges, pot_q_append], axis=0))

    n_atom_aug = int(dst_inp["pos"].shape[0])
    if n_atom_aug != base_n_atom + max(0, n_ref_index):
        raise ValueError(f"{up_file}: unexpected n_atom after reference index augmentation")

    bb_runtime_idx = np.full((n_bb, 4), -1, dtype=np.int32)
    bb_runtime_mask = np.zeros((n_bb, 4), dtype=np.int8)
    bb_runtime_w = np.zeros((n_bb, 4), dtype=np.float32)
    for k in range(n_bb):
        raw_w = np.zeros((4,), dtype=np.float32)
        for d in range(4):
            ref_idx = int(bb_ref_idx[k, d])
            if ref_idx < 0:
                continue
            run_idx = ref_offset + ref_idx
            bb_runtime_idx[k, d] = run_idx
            bb_runtime_mask[k, d] = 1
            cname = comp_names[d] if d < len(comp_names) else DEFAULT_COMPONENT_NAMES[d]
            raw_w[d] = np.float32(BB_COMPONENT_MASS.get(cname, 12.0))
        wsum = float(raw_w.sum())
        if wsum > 0.0:
            bb_runtime_w[k, :] = raw_w / wsum

    replace_dataset(bb_grp, "atom_indices", bb_runtime_idx)
    replace_dataset(bb_grp, "atom_mask", bb_runtime_mask)
    replace_dataset(bb_grp, "weights", bb_runtime_w)
    bb_grp.attrs["atom_index_space"] = np.bytes_("stage_runtime")
    bb_grp.attrs["reference_index_space"] = np.bytes_("protein_aa_pdb_0based")
    bb_grp.attrs["reference_index_offset"] = np.int32(ref_offset)
    bb_grp.attrs["reference_index_count"] = np.int32(max(0, n_ref_index))

    membership = np.full((n_atom_aug,), -1, dtype=np.int32)
    membership[:base_n_atom] = src_mem
    if n_ref_index > 0:
        membership[ref_offset:ref_offset + n_ref_index] = 0
    replace_dataset(env_grp, "protein_membership", membership)

    if env_grp["protein_membership"].shape[0] != dst_inp["pos"].shape[0]:
        raise ValueError(f"{up_file}: hybrid_env_topology/protein_membership length mismatch after augmentation")
PY
}

set_hybrid_control() {
    local up_file="$1"
    local enable="$2"
    local activation_stage="$3"
    local preprod_mode="$4"
    local prep_runtime_mode="$5"
    local active_runtime_mode="$6"
    local preprod_headgroup_roles="$7"
    python3 - "$up_file" "$enable" "$activation_stage" "$preprod_mode" "$prep_runtime_mode" "$active_runtime_mode" "$preprod_headgroup_roles" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
enable = int(sys.argv[2])
activation_stage = sys.argv[3]
preprod_mode = sys.argv[4]
prep_runtime_mode = sys.argv[5]
active_runtime_mode = sys.argv[6]
preprod_headgroup_roles = sys.argv[7]

with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    grp = inp.require_group("hybrid_control")
    grp.attrs["enable"] = np.int8(enable)
    grp.attrs["activation_stage"] = np.bytes_(activation_stage)
    grp.attrs["preprod_protein_mode"] = np.bytes_(preprod_mode)
    grp.attrs["prep_runtime_mode"] = np.bytes_(prep_runtime_mode)
    grp.attrs["active_runtime_mode"] = np.bytes_(active_runtime_mode)
    grp.attrs["preprod_lipid_headgroup_roles"] = np.bytes_(preprod_headgroup_roles)
PY
}

set_fix_rigid_from_membership() {
    local up_file="$1"
    local enable="$2"
    python3 - "$up_file" "$enable" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
enable = int(sys.argv[2])

with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    if enable:
        if "/input/hybrid_env_topology/protein_membership" not in h5:
            raise ValueError(f"{up_file}: missing /input/hybrid_env_topology/protein_membership")
        membership = h5["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
        protein_idx = np.where(membership >= 0)[0].astype(np.int32)
        grp = inp.require_group("fix_rigid")
        grp.attrs["enable"] = np.int8(1)
        if "atom_indices" in grp:
            del grp["atom_indices"]
        grp.create_dataset("atom_indices", data=protein_idx)
    else:
        grp = inp.require_group("fix_rigid")
        grp.attrs["enable"] = np.int8(0)
        if "atom_indices" in grp:
            del grp["atom_indices"]
PY
}

assert_protein_rigid_stage() {
    local stage_label="$1"
    local up_file="$2"
    local tol="${3:-$PREPROD_RIGID_TOL}"
    if [ "${PREPROD_FIX_RIGID_ENABLE}" != "1" ] || [ "${PREPROD_RIGID_ASSERT_ENABLE}" != "1" ]; then
        return
    fi
    python3 - "$stage_label" "$up_file" "$tol" << 'PY'
import sys
import h5py
import numpy as np

stage_label = sys.argv[1]
up_file = sys.argv[2]
tol = float(sys.argv[3])

with h5py.File(up_file, "r") as h5:
    if "/input/fix_rigid" not in h5:
        raise SystemExit(f"ERROR: Stage {stage_label} missing /input/fix_rigid for rigid-protein validation")
    fix_grp = h5["/input/fix_rigid"]
    if int(fix_grp.attrs.get("enable", 0)) == 0:
        raise SystemExit(f"ERROR: Stage {stage_label} has /input/fix_rigid disabled during rigid-protein validation")
    if "/input/hybrid_env_topology/protein_membership" not in h5:
        raise SystemExit(f"ERROR: Stage {stage_label} missing protein_membership for rigid validation")
    if "/output/pos" not in h5 or h5["/output/pos"].shape[0] == 0:
        raise SystemExit(f"ERROR: Stage {stage_label} has no output positions for rigid validation")

    membership = h5["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    protein_idx = np.where(membership >= 0)[0]
    if protein_idx.size == 0:
        raise SystemExit(f"ERROR: Stage {stage_label} has no protein atoms in protein_membership")

    start = h5["/input/pos"][:, :, 0][protein_idx].astype(np.float64)
    final = h5["/output/pos"][-1, 0, :, :][protein_idx].astype(np.float64)
    disp = np.linalg.norm(final - start, axis=1)
    max_disp = float(np.max(disp))
    mean_disp = float(np.mean(disp))
    print(
        f"Rigid-protein check stage {stage_label}: "
        f"n_protein={protein_idx.size} max_disp={max_disp:.6g} mean_disp={mean_disp:.6g} tol={tol:.6g}"
    )
    if not np.isfinite(max_disp) or max_disp > tol:
        raise SystemExit(
            f"ERROR: Stage {stage_label} moved protein atoms despite rigid hold "
            f"(max_disp={max_disp:.6g}, tol={tol:.6g})"
        )
PY
}

assert_production_handoff_protein_particles() {
    local source_file="$1"
    local target_file="$2"
    local tol="${3:-$PREPROD_RIGID_TOL}"
    if [ "${PREPROD_RIGID_ASSERT_ENABLE}" != "1" ]; then
        return
    fi
    python3 - "$source_file" "$target_file" "$tol" << 'PY'
import sys
import h5py
import numpy as np

source_file = sys.argv[1]
target_file = sys.argv[2]
tol = float(sys.argv[3])

with h5py.File(source_file, "r") as src, h5py.File(target_file, "r") as dst:
    if "/output/pos" not in src or src["/output/pos"].shape[0] == 0:
        raise SystemExit(f"ERROR: handoff source has no output positions: {source_file}")
    if "/input/hybrid_env_topology/protein_membership" not in src:
        raise SystemExit(f"ERROR: handoff source missing protein_membership: {source_file}")
    if "/input/hybrid_env_topology/protein_membership" not in dst:
        raise SystemExit(f"ERROR: handoff target missing protein_membership: {target_file}")

    src_mem = src["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    dst_mem = dst["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
    compare_n = min(src_mem.shape[0], dst_mem.shape[0])
    if compare_n <= 0:
        raise SystemExit("ERROR: no overlapping atom range for production handoff validation")

    carrier_offset = compare_n
    if "/input/hybrid_bb_map" in dst:
        carrier_offset = int(dst["/input/hybrid_bb_map"].attrs.get("reference_index_offset", compare_n))
    carrier_offset = max(0, min(compare_n, carrier_offset))

    protein_idx = np.where(src_mem[:carrier_offset] >= 0)[0]
    if protein_idx.size == 0:
        protein_idx = np.where(src_mem[:compare_n] >= 0)[0]
    if protein_idx.size == 0:
        raise SystemExit("ERROR: no protein particles found for production handoff validation")

    src_pos = src["/output/pos"][-1, 0, :, :][protein_idx].astype(np.float64)
    dst_pos = dst["/input/pos"][:, :, 0][protein_idx].astype(np.float64)
    disp = np.linalg.norm(dst_pos - src_pos, axis=1)
    max_disp = float(np.max(disp))
    mean_disp = float(np.mean(disp))
    print(
        "Production handoff protein-particle check: "
        f"n_protein={protein_idx.size} max_disp={max_disp:.6g} mean_disp={mean_disp:.6g} tol={tol:.6g}"
    )
    if not np.isfinite(max_disp) or max_disp > tol:
        raise SystemExit(
            f"ERROR: production handoff changed original protein particles "
            f"(max_disp={max_disp:.6g}, tol={tol:.6g})"
        )
PY
}

inject_backbone_only_production_nodes() {
    local up_file="$1"
    local protein_source="$2"
    local upside_home="$3"
    local rama_library="$4"
    local rama_sheet_mixing="$5"
    local hbond_energy="$6"
    local reference_state_rama="$7"
    if [ ! -f "${UNIVERSAL_PREP_SCRIPT}" ]; then
        echo "ERROR: prep script not found: ${UNIVERSAL_PREP_SCRIPT}"
        exit 1
    fi
    python3 "${UNIVERSAL_PREP_SCRIPT}" inject-backbone-only \
        "${up_file}" \
        "${protein_source}" \
        "${upside_home}" \
        "${rama_library}" \
        "${rama_sheet_mixing}" \
        "${hbond_energy}" \
        "${reference_state_rama}"
}

prepare_backbone_cross_artifacts() {
    if [ "${BACKBONE_CROSS_ENABLE}" != "1" ]; then
        return
    fi

    if [ "${BACKBONE_CROSS_REBUILD}" = "1" ]; then
        echo "ERROR: BACKBONE_CROSS_REBUILD is no longer supported inside run_sim_1rkl.sh"
        echo "Run the standalone parameter builder first:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi

    if [ ! -f "${BACKBONE_CROSS_BUILD_SCRIPT}" ]; then
        echo "ERROR: backbone cross builder script not found: ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi

    if [ ! -f "${BACKBONE_CROSS_TABLE_H5}" ]; then
        echo "ERROR: prebuilt backbone cross artifact not found: ${BACKBONE_CROSS_TABLE_H5}"
        echo "Run the standalone parameter builder first:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi

    local cached_schema
    local cached_calibration_mode
    local cached_dry_schema
    local cached_dry_n_type
    IFS='|' read -r cached_schema cached_calibration_mode cached_dry_schema cached_dry_n_type <<<"$(
        python3 - "${BACKBONE_CROSS_TABLE_H5}" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "r") as h5:
    if "/meta" not in h5:
        print("|||")
    else:
        meta = h5["/meta"]
        dry = h5["/dry/nonbond"] if "/dry/nonbond" in h5 else None
        print(
            f"{str(meta.attrs.get('schema', ''))}|"
            f"{str(meta.attrs.get('calibration_mode', ''))}|"
            f"{str(dry.attrs.get('schema', '')) if dry is not None else ''}|"
            f"{int(dry.attrs.get('n_type', 0)) if dry is not None else 0}"
        )
PY
    )"

    if [ "${cached_schema}" != "martini_backbone_spline_cross_v1" ]; then
        echo "ERROR: backbone cross artifact schema is '${cached_schema:-unknown}', expected martini_backbone_spline_cross_v1"
        echo "Rebuild the published artifact with:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi
    if [ "${cached_calibration_mode}" != "role_specific_proxy_pair_rows_v1" ]; then
        echo "ERROR: backbone cross artifact calibration_mode is '${cached_calibration_mode:-unknown}', expected role_specific_proxy_pair_rows_v1"
        echo "Rebuild and republish the updated artifact with:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi
    if [ "${cached_dry_schema}" != "martini_dry_nonbond_v1" ]; then
        echo "ERROR: MARTINI dry nonbond artifact schema is '${cached_dry_schema:-unknown}', expected martini_dry_nonbond_v1"
        echo "Rebuild and republish the updated artifact with:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi
    if [ "${cached_dry_n_type}" != "38" ]; then
        echo "ERROR: MARTINI dry nonbond artifact type count is '${cached_dry_n_type:-unknown}', expected 38"
        echo "Rebuild and republish the updated artifact with:"
        echo "  python3 ${BACKBONE_CROSS_BUILD_SCRIPT}"
        exit 1
    fi

    echo "=== Stage 0.5: Use Prebuilt MARTINI Force-Field Artifact ==="
    echo "Using prebuilt MARTINI force-field artifact: ${BACKBONE_CROSS_TABLE_H5}"
    echo "Cross artifact calibration: ${cached_calibration_mode}"
    echo "Dry nonbond artifact schema: ${cached_dry_schema} (${cached_dry_n_type} types)"

    if [ "${HYBRID_VALIDATE}" = "1" ]; then
        local validator_cmd=(
            python3 "${UNIVERSAL_PREP_SCRIPT}" validate-backbone-only
            --cross-artifact "${BACKBONE_CROSS_TABLE_H5}"
        )
        if [ -f "${BACKBONE_CROSS_TABLE_CSV}" ]; then
            validator_cmd+=(
                --cross-table-csv "${BACKBONE_CROSS_TABLE_CSV}"
            )
        fi
        if [ -f "${BACKBONE_CROSS_TABLE_META}" ]; then
            validator_cmd+=(
                --cross-table-meta "${BACKBONE_CROSS_TABLE_META}"
            )
        fi
        if [ -f "${DEPTH_TABLE_CSV}" ]; then
            validator_cmd+=(
                --table-csv "${DEPTH_TABLE_CSV}"
            )
        fi
        if [ -f "${DEPTH_TABLE_META}" ]; then
            validator_cmd+=(
                --table-meta "${DEPTH_TABLE_META}"
            )
        fi
        "${validator_cmd[@]}"
    fi
}

inject_backbone_cross_production_node() {
    local up_file="$1"
    if [ "${BACKBONE_CROSS_ENABLE}" != "1" ]; then
        return
    fi

    python3 "${UNIVERSAL_PREP_SCRIPT}" inject-backbone-cross \
        --artifact "${BACKBONE_CROSS_TABLE_H5}" \
        --up-files "${up_file}"
}

prepare_hybrid_artifacts() {
    echo "=== Stage 0: Hybrid Packing + Mapping Export ==="

    local packing_cutoff="${PROTEIN_LIPID_CUTOFF}"
    local packing_min_gap="nan"

    while true; do
        echo "Hybrid packing attempt with protein-lipid cutoff: ${packing_cutoff} Å"

        python3 "${UNIVERSAL_PREP_SCRIPT}" \
            --mode "both" \
            --pdb-id "${RUNTIME_PDB_ID}" \
            --runtime-pdb-output "${HYBRID_PACKED_PDB}" \
            --prepare-structure 1 \
            --protein-aa-pdb "${PROTEIN_AA_PDB}" \
            --hybrid-mapping-output "${HYBRID_MAPPING_FILE}" \
            --hybrid-bb-map-json-output "${HYBRID_PREP_DIR}/hybrid_bb_map.json" \
            --bilayer-pdb "${BILAYER_PDB}" \
            --salt-molar "${SALT_MOLAR}" \
            --protein-lipid-cutoff "${packing_cutoff}" \
            --ion-cutoff "${ION_CUTOFF}" \
            --box-padding-xy "${BOX_PADDING_XY}" \
            --box-padding-z "${BOX_PADDING_Z}" \
            --bb-aa-min-matched-residues "${BB_AA_MIN_MATCHED_RESIDUES}" \
            --bb-aa-max-rigid-rmsd "${BB_AA_MAX_RIGID_RMSD}" \
            --seed "${PREP_SEED}" \
            --summary-json "${HYBRID_PREP_DIR}/hybrid_prep_summary.json"

        if [ ! -f "${HYBRID_PACKED_PDB}" ]; then
            echo "ERROR: hybrid packed PDB not found: ${HYBRID_PACKED_PDB}"
            exit 1
        fi

        packing_min_gap="$(python3 - "${HYBRID_PREP_DIR}/hybrid_prep_summary.json" << 'PY'
import json
import math
import sys

with open(sys.argv[1], "r", encoding="utf-8") as fh:
    summary = json.load(fh)
value = float(summary.get("min_protein_lipid_distance", float("nan")))
if not math.isfinite(value):
    raise SystemExit("nan")
print(f"{value:.6f}")
PY
)"

        if python3 - "${packing_min_gap}" "${PROTEIN_LIPID_MIN_GAP}" << 'PY'
import sys
observed = float(sys.argv[1])
target = float(sys.argv[2])
raise SystemExit(0 if observed >= target else 1)
PY
        then
            echo "Hybrid packing accepted: min protein-lipid distance ${packing_min_gap} Å (target >= ${PROTEIN_LIPID_MIN_GAP} Å)"
            break
        fi

        if python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_MAX}" << 'PY'
import sys
cutoff = float(sys.argv[1])
limit = float(sys.argv[2])
raise SystemExit(0 if cutoff >= limit else 1)
PY
        then
            echo "ERROR: hybrid packing still overpacked near protein."
            echo "  Observed min protein-lipid distance: ${packing_min_gap} Å"
            echo "  Target min distance: ${PROTEIN_LIPID_MIN_GAP} Å"
            echo "  Reached cutoff limit: ${PROTEIN_LIPID_CUTOFF_MAX} Å"
            echo "  Try increasing PROTEIN_LIPID_CUTOFF_MAX and/or BOX_PADDING_XY."
            exit 1
        fi

        packing_cutoff="$(python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_STEP}" << 'PY'
import sys
cutoff = float(sys.argv[1])
step = float(sys.argv[2])
print(f"{cutoff + step:.3f}")
PY
)"
        echo "Hybrid packing too tight near protein (min ${packing_min_gap} Å). Retrying with cutoff ${packing_cutoff} Å."
    done

    if [ ! -f "${HYBRID_MAPPING_FILE}" ]; then
        echo "ERROR: hybrid mapping file not found: ${HYBRID_MAPPING_FILE}"
        exit 1
    fi

    if [ "$(python3 - "${HYBRID_PACKED_PDB}" "${RUNTIME_PDB_FILE}" << 'PY'
import os
import sys
print(int(os.path.abspath(sys.argv[1]) == os.path.abspath(sys.argv[2])))
PY
    )" = "0" ]; then
        cp -f "${HYBRID_PACKED_PDB}" "${RUNTIME_PDB_FILE}"
    fi
    echo "Runtime MARTINI PDB: ${RUNTIME_PDB_FILE}"
}

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local barostat_type="$4"
    local lipidhead_fc="${5:-0}"
    local stage_label="${6:-minimization}"

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"
    export UPSIDE_BILAYER_LIPIDHEAD_FC="$lipidhead_fc"

    python3 "${UNIVERSAL_PREP_SCRIPT}" \
        --mode "${UNIVERSAL_PREP_MODE}" \
        --pdb-id "${RUNTIME_PDB_ID}" \
        --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
        --prepare-structure 0 \
        --stage "$prepare_stage" \
        --run-dir "$RUN_DIR" \
        --summary-json "${HYBRID_PREP_DIR}/stage_${prepare_stage}.summary.json"

    local prepared_tmp="${RUN_DIR}/test.input.up"
    if [ ! -f "$prepared_tmp" ]; then
        echo "ERROR: preparation failed for stage ${prepare_stage}: $prepared_tmp not found"
        exit 1
    fi

    mv -f "$prepared_tmp" "$target_file"
    inject_hybrid_mapping "$target_file" "${HYBRID_MAPPING_FILE}"
    set_stage_label "$target_file" "$stage_label"
    set_hybrid_control \
        "$target_file" \
        "$HYBRID_CONTROL_ENABLE" \
        "$HYBRID_ACTIVATION_STAGE" \
        "$HYBRID_PREPROD_PROTEIN_MODE" \
        "$HYBRID_PREP_RUNTIME_MODE" \
        "$HYBRID_ACTIVE_RUNTIME_MODE" \
        "$HYBRID_PREPROD_LIPID_HEADGROUP_ROLES"
    if [ "$stage_label" = "production" ]; then
        set_fix_rigid_from_membership "$target_file" "0"
    elif [ "${PREPROD_FIX_RIGID_ENABLE}" = "1" ]; then
        set_fix_rigid_from_membership "$target_file" "1"
    fi
    if [ "$stage_label" = "production" ]; then
        inject_backbone_only_production_nodes \
            "$target_file" \
            "${PROTEIN_AA_PDB}" \
            "${UPSIDE_HOME}" \
            "${UPSIDE_RAMA_LIBRARY}" \
            "${UPSIDE_RAMA_SHEET_MIXING}" \
            "${UPSIDE_HBOND_ENERGY}" \
            "${UPSIDE_REFERENCE_STATE_RAMA}"
        inject_backbone_cross_production_node "$target_file"
        if [ "${HYBRID_VALIDATE}" = "1" ]; then
            local validator_cmd=(
                python3 "${UNIVERSAL_PREP_SCRIPT}" validate-backbone-only
                "$target_file"
            )
            if [ "${BACKBONE_CROSS_ENABLE}" = "1" ]; then
                validator_cmd+=(
                    --cross-artifact "${BACKBONE_CROSS_TABLE_H5}"
                )
                if [ -f "${BACKBONE_CROSS_TABLE_CSV}" ]; then
                    validator_cmd+=(
                        --cross-table-csv "${BACKBONE_CROSS_TABLE_CSV}"
                    )
                fi
                if [ -f "${BACKBONE_CROSS_TABLE_META}" ]; then
                    validator_cmd+=(
                        --cross-table-meta "${BACKBONE_CROSS_TABLE_META}"
                    )
                fi
            fi
            "${validator_cmd[@]}"
        fi
    fi

    if [ "${HYBRID_VALIDATE}" = "1" ]; then
        python3 "${UNIVERSAL_PREP_SCRIPT}" validate-stage-potentials "$target_file"
    fi

    if [ "$npt_enable" = "1" ]; then
        set_barostat_type "$target_file" "$barostat_type"
    fi
}

set_stage_npt_targets() {
    local stage_label="$1"
    case "$stage_label" in
        6.0|6.1)
            # CHARMM-GUI step6.0/6.1 minimization.mdp: ref_p = 0.0 0.0
            export UPSIDE_NPT_TARGET_PXY="$NPT_REF_P_ZERO"
            export UPSIDE_NPT_TARGET_PZ="$NPT_REF_P_ZERO"
            ;;
        6.2|6.3|6.4|6.5|6.6)
            # CHARMM-GUI step6.2-6.6 equilibration.mdp: ref_p defaults to 1 bar.
            # semi-isotropic coupling keeps z compressibility at 0.0, so z scaling remains disabled.
            export UPSIDE_NPT_TARGET_PXY="$NPT_REF_P_ONE_BAR"
            export UPSIDE_NPT_TARGET_PZ="$NPT_REF_P_ONE_BAR"
            ;;
        *)
            return
            ;;
    esac
    echo "NPT targets for stage ${stage_label}: Pxy=${UPSIDE_NPT_TARGET_PXY}, Pz=${UPSIDE_NPT_TARGET_PZ} (E_up/Angstrom^3)"
}

run_minimization_stage() {
    local stage_label="$1"
    local up_file="$2"
    local max_iter="$3"

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    echo "=== Stage ${stage_label}: Minimization ==="
    echo "Input/Output: $up_file"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$up_file"
        "--duration" "0"
        "--frame-interval" "1"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$MIN_TIME_STEP"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "v"
        "--disable-recentering"
        "--minimize"
        "--min-max-iter" "$max_iter"
        "--min-energy-tol" "1e-6"
        "--min-force-tol" "1e-3"
        "--min-step" "0.01"
    )

    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Stage ${stage_label} failed"
        exit 1
    fi
}

run_md_stage() {
    local stage_label="$1"
    local input_file="$2"
    local output_file="$3"
    local nsteps="$4"
    local dt="$5"
    local frame_steps="$6"

    local effective_frame_steps="$frame_steps"

    if (( effective_frame_steps >= nsteps )); then
        effective_frame_steps=$(( nsteps / 10 ))
        if (( effective_frame_steps < 1 )); then
            effective_frame_steps=1
        fi
        echo "NOTICE: frame_steps (${frame_steps}) >= nsteps (${nsteps}); using frame_steps=${effective_frame_steps}"
    fi

    # Frame interval remains time-based in CLI. Duration uses step-count override.
    local frame_interval
    frame_interval="$(awk -v n="$effective_frame_steps" -v dt="$dt" 'BEGIN{printf "%.10g", n*dt}')"

    if [ "$input_file" != "$output_file" ]; then
        cp -f "$input_file" "$output_file"
        handoff_initial_position "$input_file" "$output_file"
    fi

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    echo "=== Stage ${stage_label}: MD ==="
    echo "Input:  $input_file"
    echo "Output: $output_file"
    echo "nsteps=${nsteps}, dt=${dt}, duration(steps)=${nsteps}, frame_steps=${effective_frame_steps}, frame_interval(time)=${frame_interval}"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$output_file"
        "--duration-steps" "$nsteps"
        "--frame-interval" "$frame_interval"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$dt"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "v"
        "--disable-recentering"
    )

    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Stage ${stage_label} failed"
        exit 1
    fi
}

handoff_initial_position() {
    local input_file="$1"
    local output_file="$2"
    local mode="${3:-default}"

    local refresh_backbone="0"
    if [ "$mode" = "production_backbone" ]; then
        # Refresh AA carrier coordinates onto
        # current BB anchors to avoid a first-step force/energy spike.
        refresh_backbone="1"
    fi

    UPSIDE_SET_INITIAL_STRICT_COPY="$STRICT_STAGE_HANDOFF" \
    UPSIDE_SET_INITIAL_REFRESH_BACKBONE_CARRIERS="$refresh_backbone" \
        python3 "${UNIVERSAL_PREP_SCRIPT}" handoff "$input_file" "$output_file"
}

extract_stage_vtf() {
    local stage_label="$1"
    local stage_file="$2"
    local mode="$3"
    local vtf_file="${RUN_DIR}/${PDB_ID}.stage_${stage_label}.vtf"

    echo "=== Stage ${stage_label}: VTF Extraction (mode ${mode}) ==="
    echo "Input:  $stage_file"
    echo "Output: $vtf_file"
    python3 "${SCRIPT_DIR}/extract_martini_vtf.py" "$stage_file" "$vtf_file" "$stage_file" "$RUNTIME_PDB_ID" --mode "$mode"
}

echo "=== Dual Dry-MARTINI / Upside 1RKL Workflow ==="
echo "Protein ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Universal prep: ${UNIVERSAL_PREP_SCRIPT} (mode=${UNIVERSAL_PREP_MODE})"
echo "Hybrid prep: $HYBRID_PREP_DIR"
echo "Simulation stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
echo "Pre-production mode: explicit fix_rigid protein hold + dry-MARTINI environment"
echo "Production mode: dry-MARTINI bilayer + Upside backbone + backbone cross table"
echo

prepare_hybrid_artifacts
prepare_backbone_cross_artifacts

# 6.0: soft-core minimization (pre-production / hybrid inactive)
set_stage_npt_targets "6.0"
prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "0" "minimization"
cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
assert_protein_rigid_stage "6.0" "$STAGE_60_FILE"
extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

# 6.1: hard minimization (pre-production / hybrid inactive)
set_stage_npt_targets "6.1"
prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0" "0" "minimization"
cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
handoff_initial_position "$STAGE_60_FILE" "$STAGE_61_FILE"
run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
assert_protein_rigid_stage "6.1" "$STAGE_61_FILE"
extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

# 6.2: soft equilibration with rigid protein hold
set_stage_npt_targets "6.2"
prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "0" "200" "minimization"
handoff_initial_position "$STAGE_61_FILE" "$STAGE_62_FILE"
run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
assert_protein_rigid_stage "6.2" "$STAGE_62_FILE"
extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

# 6.3: reduced softening equilibration with rigid protein hold
set_stage_npt_targets "6.3"
prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "0" "100" "minimization"
cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
handoff_initial_position "$STAGE_62_FILE" "$STAGE_63_FILE"
run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
assert_protein_rigid_stage "6.3" "$STAGE_63_FILE"
extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

# 6.4-6.6: hard equilibration with restraint ramp and rigid protein hold
set_stage_npt_targets "6.4"
prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "0" "50" "minimization"
cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
handoff_initial_position "$STAGE_63_FILE" "$STAGE_64_FILE"
run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
assert_protein_rigid_stage "6.4" "$STAGE_64_FILE"
extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

set_stage_npt_targets "6.5"
prepare_stage_file "$PREPARED_65_FILE" "npt_prod" "1" "0" "20" "minimization"
cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
handoff_initial_position "$STAGE_64_FILE" "$STAGE_65_FILE"
run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
assert_protein_rigid_stage "6.5" "$STAGE_65_FILE"
extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

set_stage_npt_targets "6.6"
prepare_stage_file "$PREPARED_66_FILE" "npt_prod" "1" "0" "10" "minimization"
cp -f "$PREPARED_66_FILE" "$STAGE_66_FILE"
handoff_initial_position "$STAGE_65_FILE" "$STAGE_66_FILE"
run_md_stage "6.6" "$STAGE_66_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
assert_protein_rigid_stage "6.6" "$STAGE_66_FILE"
extract_stage_vtf "6.6" "$STAGE_66_FILE" "1"

# 7.0: production (dry bilayer + Upside backbone + cross table)
prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE" "0" "production"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
handoff_initial_position "$STAGE_66_FILE" "$STAGE_70_FILE" "production_backbone"
assert_production_handoff_protein_particles "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"

echo
echo "=== Workflow Complete ==="
echo "Hybrid prep:"
echo "  Packed PDB: ${HYBRID_PACKED_PDB}"
echo "  Mapping:    ${HYBRID_MAPPING_FILE}"
if [ "${BACKBONE_CROSS_ENABLE}" = "1" ]; then
    echo "Force-field artifacts:"
    echo "  MARTINI H5: ${BACKBONE_CROSS_TABLE_H5}"
    if [ -f "${DEPTH_TABLE_CSV}" ]; then
        echo "  Depth CSV:  ${DEPTH_TABLE_CSV}"
    fi
    if [ -f "${BACKBONE_CROSS_TABLE_CSV}" ]; then
        echo "  Cross CSV:  ${BACKBONE_CROSS_TABLE_CSV}"
    fi
fi
echo "Checkpoints:"
echo "  6.0: $STAGE_60_FILE"
echo "  6.1: $STAGE_61_FILE"
echo "  6.2: $STAGE_62_FILE"
echo "  6.3: $STAGE_63_FILE"
echo "  6.4: $STAGE_64_FILE"
echo "  6.5: $STAGE_65_FILE"
echo "  6.6: $STAGE_66_FILE"
echo "  7.0: $STAGE_70_FILE"
echo "VTF:"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.0.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.1.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.2.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.3.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.4.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.5.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_6.6.vtf"
echo "  ${RUN_DIR}/${PDB_ID}.stage_7.0.vtf"
