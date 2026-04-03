#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../source.sh"
source "${SCRIPT_DIR}/../../.venv/bin/activate"
set -euo pipefail
cd "${SCRIPT_DIR}"

# Production-only runner for 1RKL hybrid workflow.
# Assumes stages 6.0-6.6 and hybrid prep artifacts already exist.

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
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${SCRIPT_DIR}/prepare_system.py}"
UNIVERSAL_PREP_MODE="${UNIVERSAL_PREP_MODE:-both}"

RUNTIME_PDB_FILE="${RUNTIME_PDB_FILE:-${HYBRID_PREP_DIR}/${RUNTIME_PDB_ID}.MARTINI.pdb}"
RUNTIME_ITP_FILE="${RUNTIME_ITP_FILE:-${HYBRID_PREP_DIR}/${RUNTIME_PDB_ID}_proa.itp}"
HYBRID_MAPPING_FILE="${HYBRID_MAPPING_FILE:-${HYBRID_PREP_DIR}/hybrid_mapping.h5}"

STAGE_66_FILE="${STAGE_66_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up}"
PREPARED_70_FILE="${PREPARED_70_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.prepared.up}"
STAGE_70_FILE="${STAGE_70_FILE:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up}"

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-5000}"
# Production runs with injected Upside all-atom backbone nodes; use a smaller
# stable timestep than MARTINI-only stages.
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
# Leave the stage-7 protein backbone non-rigid by default in the hybrid
# workflow. Set this to 1 to re-enable the production fix_rigid mask over
# dry-MARTINI BB beads and AA backbone carrier roles BB/N/CA/C/O.
PROD_70_BACKBONE_FIX_RIGID_ENABLE="${PROD_70_BACKBONE_FIX_RIGID_ENABLE:-0}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-5000}"

# AGENTS.md: 1 E_up = 2.914952774272 kJ/mol
export UPSIDE_MARTINI_ENERGY_CONVERSION="${UPSIDE_MARTINI_ENERGY_CONVERSION:-2.914952774272}"
# Dry-MARTINI native lengths are in nm; simulation inputs use Angstrom.
export UPSIDE_MARTINI_LENGTH_CONVERSION="${UPSIDE_MARTINI_LENGTH_CONVERSION:-10}"

COMP_3E4_BAR_INV_TO_A3_PER_EUP="${COMP_3E4_BAR_INV_TO_A3_PER_EUP:-14.521180763676}"
export UPSIDE_NPT_TARGET_PXY="${UPSIDE_NPT_TARGET_PXY:-0.0}"
export UPSIDE_NPT_TARGET_PZ="${UPSIDE_NPT_TARGET_PZ:-0.0}"
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
SC_MARTINI_LIBRARY="${SC_MARTINI_LIBRARY:-${UPSIDE_HOME}/parameters/ff_2.1/martini.h5}"
SC_MARTINI_TABLE_JSON_DEFAULT="${SC_MARTINI_TABLE_JSON_DEFAULT:-${UPSIDE_HOME}/SC-training/runs/default/results/assembled/sc_table.json}"
SC_MARTINI_TABLE_JSON="${SC_MARTINI_TABLE_JSON:-${SC_MARTINI_TABLE_JSON_DEFAULT}}"
SC_MARTINI_BUILD_SCRIPT="${SC_MARTINI_BUILD_SCRIPT:-${SCRIPT_DIR}/build_sc_martini_h5.py}"
SC_TABLE_STAGE7_INJECTOR="${SC_TABLE_STAGE7_INJECTOR:-${SCRIPT_DIR}/inject_sc_table_stage7.py}"
UPSIDE_RAMA_LIBRARY="${UPSIDE_RAMA_LIBRARY:-${UPSIDE_HOME}/parameters/common/rama.dat}"
UPSIDE_RAMA_SHEET_MIXING="${UPSIDE_RAMA_SHEET_MIXING:-${UPSIDE_HOME}/parameters/ff_2.1/sheet}"
UPSIDE_HBOND_ENERGY="${UPSIDE_HBOND_ENERGY:-${UPSIDE_HOME}/parameters/ff_2.1/hbond.h5}"
UPSIDE_REFERENCE_STATE_RAMA="${UPSIDE_REFERENCE_STATE_RAMA:-${UPSIDE_HOME}/parameters/common/rama_reference.pkl}"
SC_ENV_LJ_FORCE_CAP="${SC_ENV_LJ_FORCE_CAP:-25.0}"
SC_ENV_COUL_FORCE_CAP="${SC_ENV_COUL_FORCE_CAP:-25.0}"
SC_ENV_RELAX_STEPS="${SC_ENV_RELAX_STEPS:-150}"
SC_ENV_BACKBONE_HOLD_STEPS="${SC_ENV_BACKBONE_HOLD_STEPS:-200}"
SC_ENV_PO4_Z_HOLD_STEPS="${SC_ENV_PO4_Z_HOLD_STEPS:-$SC_ENV_RELAX_STEPS}"
SC_ENV_RELAX_DT="${SC_ENV_RELAX_DT:-0.002}"
SC_ENV_RESTRAINT_K="${SC_ENV_RESTRAINT_K:-5.0}"
SC_ENV_MAX_DISPLACEMENT="${SC_ENV_MAX_DISPLACEMENT:-2.0}"
SC_ENV_PO4_Z_CLAMP_ENABLE="${SC_ENV_PO4_Z_CLAMP_ENABLE:-1}"
SC_ENV_PO4_Z_CLAMP_MODE="${SC_ENV_PO4_Z_CLAMP_MODE:-initial}"
PRODUCTION_NONPROTEIN_HARD_SPHERE="${PRODUCTION_NONPROTEIN_HARD_SPHERE:-0}"
NONPROTEIN_HS_FORCE_CAP="${NONPROTEIN_HS_FORCE_CAP:-100.0}"
NONPROTEIN_HS_POTENTIAL_CAP="${NONPROTEIN_HS_POTENTIAL_CAP:-5000.0}"
INTEGRATION_RMSD_ALIGN_ENABLE="${INTEGRATION_RMSD_ALIGN_ENABLE:-1}"

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
if [ ! -f "${SC_MARTINI_BUILD_SCRIPT}" ]; then
    echo "ERROR: SC MARTINI builder script not found: ${SC_MARTINI_BUILD_SCRIPT}"
    exit 1
fi
if [ ! -f "${SC_TABLE_STAGE7_INJECTOR}" ]; then
    echo "ERROR: stage-7 SC injector script not found: ${SC_TABLE_STAGE7_INJECTOR}"
    exit 1
fi

if [ ! -f "${RUNTIME_PDB_FILE}" ]; then
    echo "ERROR: runtime MARTINI PDB not found: ${RUNTIME_PDB_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid runtime assets."
    exit 1
fi
if [ ! -f "${RUNTIME_ITP_FILE}" ]; then
    echo "ERROR: runtime protein ITP not found: ${RUNTIME_ITP_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid runtime assets."
    exit 1
fi
if [ ! -f "${HYBRID_MAPPING_FILE}" ]; then
    echo "ERROR: hybrid mapping file not found: ${HYBRID_MAPPING_FILE}"
    echo "Run full run_sim_1rkl.sh once to generate hybrid prep outputs."
    exit 1
fi
if [ ! -f "${STAGE_66_FILE}" ]; then
    echo "ERROR: stage 6.6 checkpoint not found: ${STAGE_66_FILE}"
    echo "This script assumes stages 6.0-6.6 have already completed."
    exit 1
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

ensure_sc_martini_library() {
    local need_build="0"
    if [ -f "${SC_MARTINI_LIBRARY}" ]; then
        if python3 - "${SC_MARTINI_LIBRARY}" << 'PY'
import sys
import h5py

required = [
    "grid_nm",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_kj_mol",
    "rotamer_angular_energy_kj_mol",
    "rotamer_angular_profile",
]

path = sys.argv[1]
with h5py.File(path, "r") as h5:
    missing = [name for name in required if name not in h5]
    if missing:
        raise SystemExit(1)
PY
        then
            return
        fi
        echo "NOTICE: ${SC_MARTINI_LIBRARY} is missing rotamer-resolved SC datasets; rebuilding it."
        need_build="1"
    else
        need_build="1"
    fi

    if [ ! -f "${SC_MARTINI_TABLE_JSON}" ]; then
        echo "ERROR: SC MARTINI table JSON not found: ${SC_MARTINI_TABLE_JSON}"
        echo "The stage-7 SC injector now requires the rotamer-resolved martini.h5 schema."
        exit 1
    fi
    if [ "${need_build}" = "1" ]; then
        echo "Building ${SC_MARTINI_LIBRARY} from ${SC_MARTINI_TABLE_JSON}"
    fi
    python3 "${SC_MARTINI_BUILD_SCRIPT}" \
        --sc-table-json "${SC_MARTINI_TABLE_JSON}" \
        --output-h5 "${SC_MARTINI_LIBRARY}"
    if ! python3 - "${SC_MARTINI_LIBRARY}" << 'PY'
import sys
import h5py

required = [
    "grid_nm",
    "cos_theta_grid",
    "rotamer_count",
    "rotamer_probability_fixed",
    "rotamer_radial_energy_kj_mol",
    "rotamer_angular_energy_kj_mol",
    "rotamer_angular_profile",
]

path = sys.argv[1]
with h5py.File(path, "r") as h5:
    missing = [name for name in required if name not in h5]
    if missing:
        raise SystemExit(
            "rebuilt SC MARTINI library is still missing required datasets: "
            + ",".join(missing)
        )
PY
    then
        echo "ERROR: rebuilt SC MARTINI library failed schema validation: ${SC_MARTINI_LIBRARY}"
        exit 1
    fi
}

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

set_hybrid_control_mode() {
    local up_file="$1"
    local activation_stage="$2"
    local preprod_mode="${3:-rigid}"
    python3 - "$up_file" "$activation_stage" "$preprod_mode" << 'PY'
import sys
import h5py
import numpy as np
up_file = sys.argv[1]
activation_stage = sys.argv[2]
preprod_mode = sys.argv[3]
with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["activation_stage"] = np.bytes_(activation_stage)
    grp.attrs["preprod_protein_mode"] = np.bytes_(preprod_mode)
PY
}

set_hybrid_production_controls() {
    local up_file="$1"
    local nonprotein_hard_sphere="$2"
    local rmsd_align_enable="$3"
    python3 - "$up_file" "$nonprotein_hard_sphere" "$rmsd_align_enable" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
nonprotein_hard_sphere = int(sys.argv[2])
rmsd_align_enable = int(sys.argv[3])

with h5py.File(up_file, "r+") as h5:
    grp = h5.require_group("input").require_group("hybrid_control")
    grp.attrs["production_nonprotein_hard_sphere"] = np.int8(1 if nonprotein_hard_sphere else 0)
    grp.attrs["integration_rmsd_align_enable"] = np.int8(1 if rmsd_align_enable else 0)
PY
}

set_production_backbone_fix_rigid() {
    local up_file="$1"
    python3 - "$up_file" << 'PY'
import sys
import h5py
import numpy as np

up_file = sys.argv[1]
required_roles = ("BB", "N", "CA", "C", "O")

def as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)

with h5py.File(up_file, "r+") as h5:
    inp = h5.require_group("input")
    if "hybrid_env_topology" not in inp:
        raise SystemExit(f"ERROR: missing /input/hybrid_env_topology in {up_file}")
    if "pos" not in inp:
        raise SystemExit(f"ERROR: missing /input/pos in {up_file}")

    membership = inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32, copy=False)
    n_atom = int(inp["pos"].shape[0])
    if membership.shape[0] != n_atom:
        raise SystemExit(
            f"ERROR: {up_file} protein_membership length {membership.shape[0]} does not match n_atom {n_atom}"
        )

    role_source = ""
    if "atom_roles" in inp:
        role_source = "atom_roles"
    elif "atom_names" in inp:
        role_source = "atom_names"
    else:
        raise SystemExit(f"ERROR: missing /input/atom_roles or /input/atom_names in {up_file}")

    role_values = [as_text(x).strip().upper() for x in inp[role_source][:]]
    if len(role_values) != n_atom:
        raise SystemExit(
            f"ERROR: {up_file} {role_source} length {len(role_values)} does not match n_atom {n_atom}"
        )

    role_counts = {role: 0 for role in required_roles}
    selected = []
    for atom_idx, role in enumerate(role_values):
        if membership[atom_idx] < 0 or role not in role_counts:
            continue
        selected.append(atom_idx)
        role_counts[role] += 1

    missing_roles = [role for role, count in role_counts.items() if count == 0]
    if missing_roles:
        raise SystemExit(
            f"ERROR: {up_file} missing protein backbone roles for fix_rigid selection: "
            f"{','.join(missing_roles)}"
        )

    bb_count = role_counts["BB"]
    mismatched_roles = [role for role in required_roles[1:] if role_counts[role] != bb_count]
    if mismatched_roles:
        detail = ", ".join(f"{role}={role_counts[role]}" for role in required_roles)
        raise SystemExit(
            f"ERROR: {up_file} protein backbone role counts are inconsistent for fix_rigid selection ({detail})"
        )

    existing = np.zeros((0,), dtype=np.int32)
    if "fix_rigid" in inp and "atom_indices" in inp["fix_rigid"]:
        existing = inp["fix_rigid"]["atom_indices"][:].astype(np.int32, copy=False)

    merged = np.unique(np.concatenate([existing, np.asarray(selected, dtype=np.int32)]))

    if "fix_rigid" in inp:
        del inp["fix_rigid"]
    grp = inp.require_group("fix_rigid")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["selection"] = np.bytes_("protein_backbone_roles_bb_n_ca_c_o")
    grp.attrs["selection_source"] = np.bytes_(role_source)
    grp.attrs["selection_count"] = np.int32(merged.shape[0])
    for role, count in role_counts.items():
        grp.attrs[f"count_{role}"] = np.int32(count)
    grp.create_dataset("atom_indices", data=merged.astype(np.int32, copy=False))
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
    "hybrid_control",
    "hybrid_bb_map",
    "hybrid_env_topology",
]

BB_COMPONENT_MASS = {
    "N": 14.0,
    "CA": 12.0,
    "C": 12.0,
    "O": 16.0,
}
BB_BEAD_MASS = 72.0
BB_COMPONENT_SUM = float(sum(BB_COMPONENT_MASS.values()))
BB_COMPONENT_SCALE = (BB_BEAD_MASS / BB_COMPONENT_SUM) if BB_COMPONENT_SUM > 0.0 else 1.0
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
                ref_mass[ref_idx] = mass_dtype.type((BB_COMPONENT_SCALE * comp_mass) / 12.0)

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

inject_stage7_sc_table_nodes() {
    local up_file="$1"
    local sc_library="$2"
    python3 "${SC_TABLE_STAGE7_INJECTOR}" \
        "$up_file" \
        "$sc_library" \
        "${UPSIDE_HOME}" \
        "${UPSIDE_RAMA_LIBRARY}" \
        "${UPSIDE_RAMA_SHEET_MIXING}" \
        "${UPSIDE_HBOND_ENERGY}" \
        "${UPSIDE_REFERENCE_STATE_RAMA}" \
        --protein-itp "${RUNTIME_ITP_FILE}"
}

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local barostat_type="$4"
    local lipidhead_fc="${5:-0}"
    local stage_label="${6:-production}"

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"
    export UPSIDE_BILAYER_LIPIDHEAD_FC="$lipidhead_fc"

    python3 "${UNIVERSAL_PREP_SCRIPT}" \
        --mode "${UNIVERSAL_PREP_MODE}" \
        --pdb-id "${RUNTIME_PDB_ID}" \
        --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
        --runtime-itp-output "${RUNTIME_ITP_FILE}" \
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
    if [ "$stage_label" = "production" ]; then
        ensure_sc_martini_library
        set_hybrid_control_mode "$target_file" "production" "rigid"
        set_hybrid_production_controls \
            "$target_file" \
            "$PRODUCTION_NONPROTEIN_HARD_SPHERE" \
            "$INTEGRATION_RMSD_ALIGN_ENABLE"
        inject_stage7_sc_table_nodes \
            "$target_file" \
            "${SC_MARTINI_LIBRARY}"
        if [ "${PROD_70_BACKBONE_FIX_RIGID_ENABLE}" = "1" ]; then
            set_production_backbone_fix_rigid "$target_file"
        fi
    fi

    if [ "$npt_enable" = "1" ]; then
        set_barostat_type "$target_file" "$barostat_type"
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

    local frame_interval
    frame_interval="$(awk -v n="$effective_frame_steps" -v dt="$dt" 'BEGIN{printf "%.10g", n*dt}')"

    if [ "$input_file" != "$output_file" ]; then
        cp -f "$input_file" "$output_file"
        python3 set_initial_position.py "$input_file" "$output_file"
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

extract_stage_vtf() {
    local stage_label="$1"
    local stage_file="$2"
    local mode="$3"
    local vtf_file="${RUN_DIR}/${PDB_ID}.stage_${stage_label}.vtf"

    echo "=== Stage ${stage_label}: VTF Extraction (mode ${mode}) ==="
    echo "Input:  $stage_file"
    echo "Output: $vtf_file"
    python3 extract_martini_vtf.py "$stage_file" "$vtf_file" "$stage_file" "$RUNTIME_PDB_ID" --mode "$mode"
}

echo "=== Hybrid 1RKL Production-Only Runner ==="
echo "Protein ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Universal prep: ${UNIVERSAL_PREP_SCRIPT} (mode=${UNIVERSAL_PREP_MODE})"
echo "Assumed completed checkpoint: ${STAGE_66_FILE}"
echo

prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE" "0" "production"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
python3 set_initial_position.py "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"

echo
echo "=== Production-Only Run Complete ==="
echo "Checkpoint:"
echo "  7.0: $STAGE_70_FILE"
echo "VTF:"
echo "  ${RUN_DIR}/${PDB_ID}.stage_7.0.vtf"
