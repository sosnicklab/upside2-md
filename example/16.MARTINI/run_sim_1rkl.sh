#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CALLER_UPSIDE_HOME="${UPSIDE_HOME:-}"
if [ "${UPSIDE_SKIP_SOURCE_SH:-0}" != "1" ]; then
    source "${PROJECT_ROOT}/source.sh"
    if [ -n "${CALLER_UPSIDE_HOME}" ]; then
        export UPSIDE_HOME="${CALLER_UPSIDE_HOME}"
    fi
elif [ -n "${CALLER_UPSIDE_HOME}" ]; then
    export UPSIDE_HOME="${CALLER_UPSIDE_HOME}"
else
    export UPSIDE_HOME="${PROJECT_ROOT}"
fi
if [ -d "${PROJECT_ROOT}/.venv/bin" ]; then
    export PATH="${PROJECT_ROOT}/.venv/bin:$PATH"
fi
export PATH="${PROJECT_ROOT}/obj:$PATH"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
set -euo pipefail
cd "${SCRIPT_DIR}"
PYTHON_WORKFLOW_DIR="${PROJECT_ROOT}/py"

generate_random_seed() {
    local seed=""
    if [ -r /dev/urandom ] && command -v od >/dev/null 2>&1; then
        seed="$(od -An -N4 -tu4 /dev/urandom 2>/dev/null || true)"
        seed="${seed//[[:space:]]/}"
    fi
    if [ -z "$seed" ] || [ "$seed" = "0" ]; then
        seed="$(( ((RANDOM & 32767) << 17) ^ ((RANDOM & 32767) << 2) ^ (RANDOM & 3) ))"
    fi
    if [ -z "$seed" ] || [ "$seed" = "0" ]; then
        seed="1"
    fi
    printf '%s\n' "$seed"
}

if [ -z "${PREP_SEED:-}" ]; then
    PREP_SEED="$(generate_random_seed)"
fi
if [ -z "${SEED:-}" ]; then
    SEED="$(generate_random_seed)"
    if [ "$SEED" = "$PREP_SEED" ]; then
        SEED="$(generate_random_seed)"
    fi
fi

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
RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_aabb}"

INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_aabb}"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"
LOG_DIR="${RUN_DIR}/logs"
PREP_DIR="${RUN_DIR}/prep"

PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
BILAYER_PDB="${BILAYER_PDB:-pdb/DOPC.pdb}"
UNIVERSAL_PREP_SCRIPT="${PYTHON_WORKFLOW_DIR}/martini_prepare_system.py"
EXTRACT_VTF_SCRIPT="${PYTHON_WORKFLOW_DIR}/martini_extract_vtf.py"

RUNTIME_PDB_FILE="${PREP_DIR}/${RUNTIME_PDB_ID}.MARTINI.pdb"
BACKBONE_METADATA_FILE="${PREP_DIR}/backbone_metadata.h5"

PREPARED_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.prepared.up"
STAGE_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
PREPARED_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.prepared.up"
STAGE_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
PREPARED_62_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.2.prepared.up"
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
CONTINUE_STAGE_70_FROM="${CONTINUE_STAGE_70_FROM:-}"
CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.continue.up}"
CONTINUE_STAGE_70_LABEL="${CONTINUE_STAGE_70_LABEL:-7.0_continue}"

SALT_MOLAR="${SALT_MOLAR:-0.15}"
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-4.5}"
XY_SCALE="${XY_SCALE:-1.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"

readonly ION_CUTOFF="4.0"
readonly PROTEIN_LIPID_MIN_GAP="4.5"
readonly PROTEIN_LIPID_CUTOFF_STEP="0.5"
readonly PROTEIN_LIPID_CUTOFF_MAX="8.0"

TEMPERATURE="${TEMPERATURE:-0.8647}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-10000}"

readonly MIN_TIME_STEP="0.010"
readonly EQ_TIME_STEP="0.010"
readonly PROD_TIME_STEP="0.002"
readonly EQ_FRAME_STEPS="1000"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-50}"

readonly BAR_1_TO_EUP_PER_A3="0.000020659477"
readonly NPT_REF_P_ZERO="0.0"
readonly NPT_REF_P_ONE_BAR="${BAR_1_TO_EUP_PER_A3}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"
STAGE_NPT_TARGET_PXY="${NPT_REF_P_ZERO}"
STAGE_NPT_TARGET_PZ="${NPT_REF_P_ZERO}"

SC_MARTINI_LIBRARY="${UPSIDE_HOME}/parameters/ff_2.1/martini.h5"
SC_MARTINI_TABLE_JSON="${UPSIDE_HOME}/SC-training/runs/default/results/assembled/sc_table.json"
UPSIDE_RAMA_LIBRARY="${UPSIDE_HOME}/parameters/common/rama.dat"
UPSIDE_RAMA_SHEET_MIXING="${UPSIDE_HOME}/parameters/ff_2.1/sheet"
UPSIDE_HBOND_ENERGY="${UPSIDE_HOME}/parameters/ff_2.1/hbond.h5"
UPSIDE_REFERENCE_STATE_RAMA="${UPSIDE_HOME}/parameters/common/rama_reference.pkl"

if [ -z "${UPSIDE_HOME:-}" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi
for required_file in \
    "${PROTEIN_AA_PDB}" \
    "${BILAYER_PDB}" \
    "${UNIVERSAL_PREP_SCRIPT}" \
    "${EXTRACT_VTF_SCRIPT}" \
    "${UPSIDE_RAMA_LIBRARY}" \
    "${UPSIDE_HBOND_ENERGY}" \
    "${UPSIDE_REFERENCE_STATE_RAMA}"
do
    if [ ! -f "${required_file}" ]; then
        echo "ERROR: required file not found: ${required_file}"
        exit 1
    fi
done

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$PREP_DIR"

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

with h5py.File(sys.argv[1], "r") as h5:
    missing = [name for name in required if name not in h5]
    if missing:
        raise SystemExit(1)
PY
        then
            return
        fi
        echo "NOTICE: ${SC_MARTINI_LIBRARY} is missing rotamer-resolved datasets; rebuilding it."
        need_build="1"
    else
        need_build="1"
    fi

    if [ ! -f "${SC_MARTINI_TABLE_JSON}" ]; then
        echo "ERROR: SC MARTINI table JSON not found: ${SC_MARTINI_TABLE_JSON}"
        exit 1
    fi
    if [ "${need_build}" = "1" ]; then
        python3 "${UNIVERSAL_PREP_SCRIPT}" build-sc-martini-h5 \
            --sc-table-json "${SC_MARTINI_TABLE_JSON}" \
            --output-h5 "${SC_MARTINI_LIBRARY}"
    fi
}

set_stage_label() {
    local up_file="$1"
    local stage_label="$2"
    python3 - "$up_file" "$stage_label" << 'PY'
import sys
import h5py
import numpy as np

with h5py.File(sys.argv[1], "r+") as h5:
    grp = h5.require_group("input").require_group("stage_parameters")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["current_stage"] = np.bytes_(sys.argv[2])
PY
}

inject_runtime_metadata() {
    local up_file="$1"
    local metadata_file="$2"
    python3 - "$up_file" "$metadata_file" << 'PY'
import sys
import h5py

up_file = sys.argv[1]
metadata_file = sys.argv[2]
groups = ("hybrid_bb_map", "hybrid_env_topology", "sequence")

with h5py.File(metadata_file, "r") as src, h5py.File(up_file, "r+") as dst:
    src_inp = src["/input"]
    dst_inp = dst.require_group("input")
    for name in groups:
        if name not in src_inp:
            raise ValueError(f"Missing metadata entry: /input/{name}")
        if name in dst_inp:
            del dst_inp[name]
        src.copy(src_inp[name], dst_inp, name=name)
PY
}

promote_stage_from_previous() {
    local source_file="$1"
    local target_file="$2"
    local stage_label="$3"
    local npt_enable="$4"
    local lipidhead_fc="$5"
    local target_pxy="$6"
    local target_pz="$7"
    cp -f "$source_file" "$target_file"
    python3 - "$target_file" "$stage_label" "$npt_enable" "$lipidhead_fc" "$target_pxy" "$target_pz" << 'PY'
import re
import sys
import h5py
import numpy as np

target_file = sys.argv[1]
stage_label = sys.argv[2]
npt_enable = int(sys.argv[3])
lipidhead_fc = float(sys.argv[4])
target_pxy = float(sys.argv[5])
target_pz = float(sys.argv[6])

energy_conversion = 2.914952774272
length_conversion = 10.0
restraint_spring = np.float32(0.0)
if lipidhead_fc > 0.0:
    restraint_spring = np.float32(
        lipidhead_fc / (energy_conversion * length_conversion * length_conversion)
    )

def as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)

with h5py.File(target_file, "r+") as h5:
    inp = h5["/input"]
    if "/output/pos" in h5 and h5["/output/pos"].shape[0] > 0:
        last_pos = h5["/output/pos"][-1, 0, :, :]
    else:
        last_pos = inp["pos"][:, :, 0]

    last_box = None
    if "/output/box" in h5 and h5["/output/box"].size > 0:
        last_box = h5["/output/box"][-1]
        if len(last_box.shape) == 2 and last_box.shape[1] == 3:
            last_box = last_box[0]
    elif "/input/potential/martini_potential" in h5:
        pot_grp = h5["/input/potential/martini_potential"]
        if all(k in pot_grp.attrs for k in ("x_len", "y_len", "z_len")):
            last_box = np.array(
                [pot_grp.attrs["x_len"], pot_grp.attrs["y_len"], pot_grp.attrs["z_len"]],
                dtype=np.float32,
            )

    for key in list(h5.keys()):
        if key == "output" or re.fullmatch(r"output_previous_\\d+", str(key)):
            del h5[key]

    if "pos" in inp:
        del inp["pos"]
    inp.create_dataset("pos", data=last_pos[:, :, np.newaxis])

    stage_grp = inp.require_group("stage_parameters")
    stage_grp.attrs["enable"] = np.int8(1)
    stage_grp.attrs["current_stage"] = np.bytes_(stage_label)

    pot = inp["potential"]
    if "restraint_position" in pot:
        del pot["restraint_position"]
    if restraint_spring > 0.0:
        membership = inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32)
        atom_names = [as_text(x).strip().upper() for x in inp["atom_names"][:]]
        restraint_indices = [
            atom_idx
            for atom_idx, (membership_value, atom_name) in enumerate(zip(membership, atom_names))
            if membership_value < 0 and atom_name == "PO4"
        ]
        if not restraint_indices:
            raise ValueError("Cannot create lipid-head restraints: no environment PO4 atoms found")
        restraint_indices = np.asarray(restraint_indices, dtype=np.int32)
        ref_pos = last_pos[restraint_indices].astype(np.float32)
        spring_const_xyz = np.zeros((len(restraint_indices), 3), dtype=np.float32)
        spring_const_xyz[:, 2] = restraint_spring
        grp = pot.create_group("restraint_position")
        grp.attrs["arguments"] = np.asarray([np.bytes_("pos")])
        grp.attrs["initialized"] = np.uint8(1)
        grp.create_dataset("restraint_indices", data=restraint_indices, dtype=np.int32)
        grp.create_dataset("ref_pos", data=ref_pos, dtype=np.float32)
        grp.create_dataset("spring_const_xyz", data=spring_const_xyz, dtype=np.float32)
        grp.create_dataset(
            "spring_const",
            data=np.full((len(restraint_indices),), restraint_spring, dtype=np.float32),
            dtype=np.float32,
        )

    if "fix_rigid" in inp:
        del inp["fix_rigid"]
    if "barostat" in inp:
        del inp["barostat"]

    if npt_enable:
        grp = inp.create_group("barostat")
        grp.attrs["enable"] = np.int8(npt_enable)
        grp.attrs["target_p_xy"] = target_pxy
        grp.attrs["target_p_z"] = target_pz
        grp.attrs["tau_p"] = 4.0
        grp.attrs["compressibility_xy"] = 14.521180763676
        grp.attrs["compressibility_z"] = 0.0
        grp.attrs["interval"] = 10
        grp.attrs["semi_isotropic"] = 1
        grp.attrs["debug"] = 1

    if last_box is not None and "martini_potential" in pot:
        martini_potential = pot["martini_potential"]
        martini_potential.attrs["x_len"] = float(last_box[0])
        martini_potential.attrs["y_len"] = float(last_box[1])
        martini_potential.attrs["z_len"] = float(last_box[2])
PY
}

set_preproduction_spatial_holds() {
    local up_file="$1"
    python3 - "$up_file" << 'PY'
import sys
import h5py
import numpy as np

required_roles = ("N", "CA", "C", "O")

def as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)

with h5py.File(sys.argv[1], "r+") as h5:
    inp = h5["/input"]
    membership = inp["hybrid_env_topology"]["protein_membership"][:].astype(np.int32)
    roles = [as_text(value).strip().upper() for value in inp["atom_roles"][:]]
    atom_names = [as_text(value).strip().upper() for value in inp["atom_names"][:]]

    fixed_atoms = []
    counts = {role: 0 for role in required_roles}
    for atom_idx, role in enumerate(roles):
        if membership[atom_idx] < 0 or role not in counts:
            continue
        fixed_atoms.append(atom_idx)
        counts[role] += 1

    missing = [role for role, count in counts.items() if count == 0]
    if missing:
        raise ValueError(f"Missing AA backbone roles for preproduction hold: {missing}")

    ca_count = counts["CA"]
    bad = [role for role in required_roles if counts[role] != ca_count]
    if bad:
        detail = ", ".join(f"{role}={counts[role]}" for role in required_roles)
        raise ValueError(f"Inconsistent AA backbone role counts for preproduction hold ({detail})")

    z_fixed_atoms = []
    for atom_idx, atom_name in enumerate(atom_names):
        if membership[atom_idx] >= 0:
            continue
        if atom_name == "PO4":
            z_fixed_atoms.append(atom_idx)
    if not z_fixed_atoms:
        raise ValueError("Requested PO4 Z hold but found no environment PO4 atoms")

    if "fix_rigid" in inp:
        del inp["fix_rigid"]
    grp = inp.create_group("fix_rigid")
    grp.attrs["enable"] = np.int8(1)
    grp.attrs["mode"] = np.bytes_("rigid_body")
    grp.attrs["selection"] = np.bytes_("aa_backbone_rigid_body_hold")
    grp.attrs["selection_count"] = np.int32(len(fixed_atoms))
    grp.create_dataset("atom_indices", data=np.asarray(fixed_atoms, dtype=np.int32))

    if z_fixed_atoms:
        grp.attrs["z_selection"] = np.bytes_("environment_po4_z_hold")
        grp.attrs["z_selection_count"] = np.int32(len(z_fixed_atoms))
        grp.create_dataset("z_atom_indices", data=np.asarray(z_fixed_atoms, dtype=np.int32))
PY
}

validate_production_stage_file() {
    local up_file="$1"
    python3 - "$up_file" << 'PY'
import sys
import h5py
import numpy as np

def as_text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return str(value)

with h5py.File(sys.argv[1], "r") as h5:
    inp = h5["/input"]
    stage = as_text(inp["stage_parameters"].attrs.get("current_stage", b"")).strip().lower()
    if stage != "production":
        raise ValueError(f"Expected production stage file, found {stage!r}")
    if "sequence" not in inp:
        raise ValueError("Missing /input/sequence")
    if "hybrid_bb_map" not in inp:
        raise ValueError("Missing /input/hybrid_bb_map")
    if "hybrid_env_topology" not in inp:
        raise ValueError("Missing /input/hybrid_env_topology")
    membership = inp["hybrid_env_topology"]["protein_membership"][:]
    if membership.shape[0] != inp["pos"].shape[0]:
        raise ValueError("protein_membership length mismatch")
PY
}

inject_sidechain_nodes() {
    local up_file="$1"
    ensure_sc_martini_library
    python3 "${UNIVERSAL_PREP_SCRIPT}" inject-stage7-sc \
        "$up_file" \
        "${SC_MARTINI_LIBRARY}" \
        "${UPSIDE_HOME}" \
        "${UPSIDE_RAMA_LIBRARY}" \
        "${UPSIDE_RAMA_SHEET_MIXING}" \
        "${UPSIDE_HBOND_ENERGY}" \
        "${UPSIDE_REFERENCE_STATE_RAMA}"
}

prepare_backbone_artifacts() {
    echo "=== Stage 0: AA Backbone Packing Export ==="
    local packing_cutoff="${PROTEIN_LIPID_CUTOFF}"
    local packing_min_gap="nan"

    while true; do
        echo "Packing attempt with protein-lipid cutoff: ${packing_cutoff} Å"
        python3 "${UNIVERSAL_PREP_SCRIPT}" \
            --pdb-id "${RUNTIME_PDB_ID}" \
            --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
            --prepare-structure 1 \
            --protein-aa-pdb "${PROTEIN_AA_PDB}" \
            --hybrid-mapping-output "${BACKBONE_METADATA_FILE}" \
            --bilayer-pdb "${BILAYER_PDB}" \
            --salt-molar "${SALT_MOLAR}" \
            --protein-lipid-cutoff "${packing_cutoff}" \
            --ion-cutoff "${ION_CUTOFF}" \
            --xy-scale "${XY_SCALE}" \
            --box-padding-xy "${BOX_PADDING_XY}" \
            --box-padding-z "${BOX_PADDING_Z}" \
            --seed "${PREP_SEED}"

        packing_min_gap="$(python3 - "${RUNTIME_PDB_FILE}" << 'PY'
import sys
import numpy as np

protein_names = {"N", "CA", "C", "O"}
lipid_residues = {"DOP", "DOPC"}
protein_xyz = []
lipid_xyz = []
with open(sys.argv[1], "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip().upper()
        resname = line[17:21].strip().upper()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom_name in protein_names:
            protein_xyz.append(xyz)
        elif resname in lipid_residues:
            lipid_xyz.append(xyz)
if not protein_xyz or not lipid_xyz:
    raise SystemExit("nan")
protein_xyz = np.asarray(protein_xyz, dtype=float)
lipid_xyz = np.asarray(lipid_xyz, dtype=float)
min_d2 = float("inf")
chunk = 2000
for i in range(0, lipid_xyz.shape[0], chunk):
    block = lipid_xyz[i:i + chunk]
    diff = block[:, None, :] - protein_xyz[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", diff, diff)
    min_d2 = min(min_d2, float(d2.min()))
print(f"{min_d2 ** 0.5:.6f}")
PY
)"

        if python3 - "${packing_min_gap}" "${PROTEIN_LIPID_MIN_GAP}" << 'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= float(sys.argv[2]) else 1)
PY
        then
            echo "Packing accepted: min protein-lipid distance ${packing_min_gap} Å"
            break
        fi

        if python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_MAX}" << 'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= float(sys.argv[2]) else 1)
PY
        then
            echo "ERROR: packing remained too tight near the AA backbone."
            echo "  Observed min protein-lipid distance: ${packing_min_gap} Å"
            echo "  Target min distance: ${PROTEIN_LIPID_MIN_GAP} Å"
            echo "  Reached cutoff limit: ${PROTEIN_LIPID_CUTOFF_MAX} Å"
            exit 1
        fi

        packing_cutoff="$(python3 - "${packing_cutoff}" "${PROTEIN_LIPID_CUTOFF_STEP}" << 'PY'
import sys
print(f"{float(sys.argv[1]) + float(sys.argv[2]):.3f}")
PY
)"
        echo "Packing too tight near the AA backbone (min ${packing_min_gap} Å). Retrying with cutoff ${packing_cutoff} Å."
    done
}

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local lipidhead_fc="${4:-0}"
    local stage_label="${5:-minimization}"
    local target_pxy="${6:-${NPT_REF_P_ZERO}}"
    local target_pz="${7:-${NPT_REF_P_ZERO}}"

    python3 "${UNIVERSAL_PREP_SCRIPT}" \
        --pdb-id "${RUNTIME_PDB_ID}" \
        --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
        --prepare-structure 0 \
        --stage "$prepare_stage" \
        --run-dir "$RUN_DIR" \
        --protein-aa-pdb "${PROTEIN_AA_PDB}" \
        --stage-lipidhead-fc "${lipidhead_fc}" \
        --npt-enable "${npt_enable}" \
        --npt-target-pxy "${target_pxy}" \
        --npt-target-pz "${target_pz}"

    local prepared_tmp="${RUN_DIR}/test.input.up"
    if [ ! -f "$prepared_tmp" ]; then
        echo "ERROR: preparation failed for stage ${prepare_stage}: ${prepared_tmp} not found"
        exit 1
    fi

    mv -f "$prepared_tmp" "$target_file"
    inject_runtime_metadata "$target_file" "${BACKBONE_METADATA_FILE}"
    set_stage_label "$target_file" "$stage_label"
    inject_sidechain_nodes "$target_file"

    if [ "$stage_label" != "production" ]; then
        set_preproduction_spatial_holds "$target_file"
    fi
}

set_stage_npt_targets() {
    local stage_label="$1"
    case "$stage_label" in
        6.0|6.1)
            STAGE_NPT_TARGET_PXY="$NPT_REF_P_ZERO"
            STAGE_NPT_TARGET_PZ="$NPT_REF_P_ZERO"
            ;;
        6.2|6.3|6.4|6.5|6.6)
            STAGE_NPT_TARGET_PXY="$NPT_REF_P_ONE_BAR"
            STAGE_NPT_TARGET_PZ="$NPT_REF_P_ONE_BAR"
            ;;
        *)
            STAGE_NPT_TARGET_PXY="$NPT_REF_P_ZERO"
            STAGE_NPT_TARGET_PZ="$NPT_REF_P_ZERO"
            return
            ;;
    esac
}

run_minimization_stage() {
    local stage_label="$1"
    local up_file="$2"
    local max_iter="$3"
    local log_file="${LOG_DIR}/stage_${stage_label}.log"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$up_file"
        "--duration" "0"
        "--frame-interval" "1"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$MIN_TIME_STEP"
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
        echo "ERROR: Stage ${stage_label} minimization failed"
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
    fi

    local frame_interval
    frame_interval="$(awk -v n="$effective_frame_steps" -v dt="$dt" 'BEGIN{printf "%.10g", n*dt}')"

    if [ "$input_file" != "$output_file" ]; then
        cp -f "$input_file" "$output_file"
        handoff_initial_position "$input_file" "$output_file"
    fi

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$output_file"
        "--duration-steps" "$nsteps"
        "--frame-interval" "$frame_interval"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$dt"
        "--seed" "$SEED"
        "--integrator" "v"
        "--disable-recentering"
    )

    if ! "${cmd[@]}" 2>&1 | tee "$log_file"; then
        echo "ERROR: Stage ${stage_label} MD failed"
        exit 1
    fi
}

handoff_initial_position() {
    local input_file="$1"
    local output_file="$2"
    UPSIDE_SET_INITIAL_STRICT_COPY="1" \
    UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS="0" \
    UPSIDE_SET_INITIAL_RECENTER_PRODUCTION="0" \
        python3 "${UNIVERSAL_PREP_SCRIPT}" set-initial-position "$input_file" "$output_file"
}

list_stage_vtf_outputs() {
    local stage_label="$1"
    local stage_file="$2"
    python3 - "$RUN_DIR" "$PDB_ID" "$stage_label" "$stage_file" << 'PY'
import os
import re
import sys
import h5py

run_dir, pdb_id, stage_label, stage_file = sys.argv[1:5]
prefix = os.path.join(run_dir, f"{pdb_id}.stage_{stage_label}")
pattern = re.compile(r"^output_previous_(\d+)$")
with h5py.File(stage_file, "r") as f:
    groups = []
    for key in f.keys():
        match = pattern.fullmatch(str(key))
        if match and isinstance(f[key], h5py.Group):
            groups.append((int(match.group(1)), str(key)))
    groups.sort()
    output_groups = [name for _, name in groups]
    if "output" in f and isinstance(f["output"], h5py.Group):
        output_groups.append("output")
if len(output_groups) <= 1:
    print(f"{prefix}.vtf")
else:
    for segment_index, _ in enumerate(output_groups):
        print(f"{prefix}.segment_{segment_index}.vtf")
PY
}

print_stage_vtf_outputs() {
    local stage_label="$1"
    local stage_file="$2"
    while IFS= read -r vtf_file; do
        [ -n "$vtf_file" ] || continue
        echo "  ${vtf_file}"
    done < <(list_stage_vtf_outputs "$stage_label" "$stage_file")
}

extract_stage_vtf() {
    local stage_label="$1"
    local stage_file="$2"
    local mode="$3"
    local vtf_file="${RUN_DIR}/${PDB_ID}.stage_${stage_label}.vtf"
    python3 "${EXTRACT_VTF_SCRIPT}" "$stage_file" "$vtf_file" "$stage_file" "$RUNTIME_PDB_ID" --mode "$mode" --split-segments
}

run_stage70_continuation() {
    local source_file="$1"
    local output_file="$2"
    local stage_label="${3:-7.0_continue}"

    if [ ! -f "$source_file" ]; then
        echo "ERROR: continuation source not found: ${source_file}"
        exit 1
    fi

    validate_production_stage_file "$source_file"
    mkdir -p "$(dirname "$output_file")"

    if [ "$(python3 - "$source_file" "$output_file" << 'PY'
import os
import sys
print(int(os.path.abspath(sys.argv[1]) == os.path.abspath(sys.argv[2])))
PY
)" = "0" ]; then
        cp -f "$source_file" "$output_file"
    fi

    handoff_initial_position "$source_file" "$output_file"
    validate_production_stage_file "$output_file"
    run_md_stage "$stage_label" "$output_file" "$output_file" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
    extract_stage_vtf "$stage_label" "$output_file" "2"
}

echo "=== AA-Backbone 1RKL Workflow ==="
echo "Protein ID: ${PDB_ID}"
echo "Runtime PDB ID: ${RUNTIME_PDB_ID}"
echo "Preparation seed: ${PREP_SEED}"
echo "Simulation seed: ${SEED}"
echo "Prep script: ${UNIVERSAL_PREP_SCRIPT}"
echo "VTF extractor: ${EXTRACT_VTF_SCRIPT}"
echo "Prep dir: ${PREP_DIR}"
if [ -n "${CONTINUE_STAGE_70_FROM}" ]; then
    echo "Continuation mode: production stage 7.0 only"
    echo "Continuation source: ${CONTINUE_STAGE_70_FROM}"
    echo "Continuation output: ${CONTINUE_STAGE_70_OUTPUT}"
else
    echo "Simulation stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
fi
echo

if [ -n "${CONTINUE_STAGE_70_FROM}" ]; then
    run_stage70_continuation "${CONTINUE_STAGE_70_FROM}" "${CONTINUE_STAGE_70_OUTPUT}" "${CONTINUE_STAGE_70_LABEL}"
else
    prepare_backbone_artifacts

    set_stage_npt_targets "6.0"
    prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "minimization" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
    run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
    extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

    set_stage_npt_targets "6.1"
    promote_stage_from_previous "$STAGE_60_FILE" "$PREPARED_61_FILE" "minimization" "1" "0" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_61_FILE"
    cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
    run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
    extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

    set_stage_npt_targets "6.2"
    promote_stage_from_previous "$STAGE_61_FILE" "$PREPARED_62_FILE" "minimization" "1" "200" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_62_FILE"
    cp -f "$PREPARED_62_FILE" "$STAGE_62_FILE"
    run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

    set_stage_npt_targets "6.3"
    promote_stage_from_previous "$STAGE_62_FILE" "$PREPARED_63_FILE" "minimization" "1" "100" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_63_FILE"
    cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
    run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

    set_stage_npt_targets "6.4"
    promote_stage_from_previous "$STAGE_63_FILE" "$PREPARED_64_FILE" "minimization" "1" "50" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_64_FILE"
    cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
    run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

    set_stage_npt_targets "6.5"
    promote_stage_from_previous "$STAGE_64_FILE" "$PREPARED_65_FILE" "minimization" "1" "20" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_65_FILE"
    cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
    run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

    set_stage_npt_targets "6.6"
    promote_stage_from_previous "$STAGE_65_FILE" "$PREPARED_66_FILE" "minimization" "1" "10" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    set_preproduction_spatial_holds "$PREPARED_66_FILE"
    cp -f "$PREPARED_66_FILE" "$STAGE_66_FILE"
    run_md_stage "6.6" "$STAGE_66_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
    extract_stage_vtf "6.6" "$STAGE_66_FILE" "1"

    promote_stage_from_previous "$STAGE_66_FILE" "$PREPARED_70_FILE" "production" "$PROD_70_NPT_ENABLE" "0" "$STAGE_NPT_TARGET_PXY" "$STAGE_NPT_TARGET_PZ"
    cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
    validate_production_stage_file "$STAGE_70_FILE"
    run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
    extract_stage_vtf "7.0" "$STAGE_70_FILE" "2"
fi

echo
echo "=== Workflow Complete ==="
if [ -n "${CONTINUE_STAGE_70_FROM}" ]; then
    echo "Continuation source: ${CONTINUE_STAGE_70_FROM}"
    echo "Continuation checkpoint: ${CONTINUE_STAGE_70_OUTPUT}"
    echo "Continuation VTF:"
    print_stage_vtf_outputs "${CONTINUE_STAGE_70_LABEL}" "${CONTINUE_STAGE_70_OUTPUT}"
else
    echo "Prep:"
    echo "  Packed PDB: ${RUNTIME_PDB_FILE}"
    echo "  Metadata:   ${BACKBONE_METADATA_FILE}"
    echo "Checkpoints:"
    echo "  6.0: ${STAGE_60_FILE}"
    echo "  6.1: ${STAGE_61_FILE}"
    echo "  6.2: ${STAGE_62_FILE}"
    echo "  6.3: ${STAGE_63_FILE}"
    echo "  6.4: ${STAGE_64_FILE}"
    echo "  6.5: ${STAGE_65_FILE}"
    echo "  6.6: ${STAGE_66_FILE}"
    echo "  7.0: ${STAGE_70_FILE}"
    echo "VTF:"
    print_stage_vtf_outputs "6.0" "$STAGE_60_FILE"
    print_stage_vtf_outputs "6.1" "$STAGE_61_FILE"
    print_stage_vtf_outputs "6.2" "$STAGE_62_FILE"
    print_stage_vtf_outputs "6.3" "$STAGE_63_FILE"
    print_stage_vtf_outputs "6.4" "$STAGE_64_FILE"
    print_stage_vtf_outputs "6.5" "$STAGE_65_FILE"
    print_stage_vtf_outputs "6.6" "$STAGE_66_FILE"
    print_stage_vtf_outputs "7.0" "$STAGE_70_FILE"
fi
