#!/bin/bash
set -eo pipefail

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

set -u

if [ -d "${PROJECT_ROOT}/.venv/bin" ]; then
    export PATH="${PROJECT_ROOT}/.venv/bin:$PATH"
fi
export PATH="${PROJECT_ROOT}/obj:$PATH"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"

cd "${SCRIPT_DIR}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        PDB_ID=*)
            PDB_ID="${1#*=}"
            shift
            ;;
        *)
            echo "Unknown parameter $1" >&2
            exit 1
            ;;
    esac
done

PDB_ID="${PDB_ID:-1rkl}"
TEMPERATURE="${TEMPERATURE:-0.8}"

# Full-resolution dynamics use one g-JF step for every physical degree of freedom. Production maps the
# requested factor-four-corrected 40 ps clock onto the bare mobility of each MARTINI particle. A protein
# carrier receives the additive FDT friction of the DOPC beads inside the 12 A interaction range. Molecular
# DOPC COM diffusion remains an independently measured validation observable.
export UPSIDE_PROTEIN_TIME_PS_PER_STEP="${UPSIDE_PROTEIN_TIME_PS_PER_STEP:-40.0}"
export UPSIDE_MARTINI_TIME_FACTOR="${UPSIDE_MARTINI_TIME_FACTOR:-4.0}"
export UPSIDE_DOPC_TARGET_DIFFUSION_UM2_S="${UPSIDE_DOPC_TARGET_DIFFUSION_UM2_S:-11.5}"
export UPSIDE_DOPC_REFERENCE_TEMPERATURE_UP="${TEMPERATURE}"
export UPSIDE_DRY_MARTINI_RELAXATION_PS="${UPSIDE_DRY_MARTINI_RELAXATION_PS:-4.0}"

RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_hybrid_full}"
RUN_DIR="${RUN_DIR:-outputs/martini_${PDB_ID}_hybrid_full}"
PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
BILAYER_PDB="${BILAYER_PDB:-${UPSIDE_HOME}/parameters/dryMARTINI/DOPC.pdb}"
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${PROJECT_ROOT}/py/martini_prepare_system.py}"
EXTRACT_VTF_SCRIPT="${EXTRACT_VTF_SCRIPT:-${PROJECT_ROOT}/py/martini_extract_vtf.py}"

SALT_MOLAR="${SALT_MOLAR:-0.15}"
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-0.0}"
ION_CUTOFF="${ION_CUTOFF:-10.0}"
XY_SCALE="${XY_SCALE:-1.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-embed}"
PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-input}"
PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-0.0}"
PROTEIN_LIPID_CUTOFF_STEP="${PROTEIN_LIPID_CUTOFF_STEP:-0.5}"
PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-8.0}"

THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-5.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
STRICT_STAGE_HANDOFF="${STRICT_STAGE_HANDOFF:-1}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-0}"
MIN_70_MAX_ITER="${MIN_70_MAX_ITER:-500}"
EQ_60_NSTEPS="${EQ_60_NSTEPS:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_BURNIN_NSTEPS="${PROD_70_BURNIN_NSTEPS:-40000}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-50000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.009}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.009}"
EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-50}"
if [ "${EQ_TIME_STEP}" != "${PROD_TIME_STEP}" ]; then
    echo "ERROR: equilibration and production must use the same MARTINI timestep" >&2
    exit 1
fi
export UPSIDE_MARTINI_TIME_STEP_UP="${PROD_TIME_STEP}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"

PREP_SEED="${PREP_SEED:-}"
SEED="${SEED:-}"
CONTINUE_STAGE_70_FROM="${CONTINUE_STAGE_70_FROM:-}"
CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-}"
CONTINUE_STAGE_70_LABEL="${CONTINUE_STAGE_70_LABEL:-}"
AUTO_CONTINUE_STAGE_70="${AUTO_CONTINUE_STAGE_70:-0}"

if [ ! -f "${UNIVERSAL_PREP_SCRIPT}" ]; then
    echo "ERROR: preparation script not found: ${UNIVERSAL_PREP_SCRIPT}" >&2
    exit 1
fi

find_latest_stage7() {
    local best_file=""
    local best_idx=-1
    local path=""
    local file=""
    local idx=""
    shopt -s nullglob
    for path in "${RUN_DIR}/checkpoints/${PDB_ID}.stage_7."*.up; do
        [ -f "$path" ] || continue
        file="$(basename "$path")"
        if [[ "$file" =~ ^${PDB_ID}\.stage_7\.([0-9]+)\.up$ ]]; then
            idx="${BASH_REMATCH[1]}"
            if [ "$idx" -gt "$best_idx" ]; then
                best_idx="$idx"
                best_file="$path"
            fi
        fi
    done
    shopt -u nullglob
    if [ -n "$best_file" ]; then
        printf '%s\n' "$best_file"
    fi
}

ensure_martini_parameter_files() {
    local martini_h5="${UPSIDE_HOME}/parameters/ff_2.1/martini.h5"
    if [ ! -f "${martini_h5}" ]; then
        echo "MARTINI force-field file missing. Generating..."
        python3 "${PROJECT_ROOT}/py/martini_gen_params.py" --upside-home "${UPSIDE_HOME}"
    fi
}

if [ -z "${CONTINUE_STAGE_70_FROM}" ] && [ "${AUTO_CONTINUE_STAGE_70}" = "1" ]; then
    CONTINUE_STAGE_70_FROM="$(find_latest_stage7)"
fi

if [ -n "${CONTINUE_STAGE_70_FROM}" ]; then
    echo "Detected continuation source: ${CONTINUE_STAGE_70_FROM}"
fi

ensure_martini_parameter_files

workflow_args=(
    run-hybrid-workflow
    --pdb-id "${PDB_ID}"
    --runtime-pdb-id "${RUNTIME_PDB_ID}"
    --upside-home "${UPSIDE_HOME}"
    --run-dir "${RUN_DIR}"
    --protein-aa-pdb "${PROTEIN_AA_PDB}"
    --bilayer-pdb "${BILAYER_PDB}"
    --extract-vtf-script "${EXTRACT_VTF_SCRIPT}"
    --salt-molar "${SALT_MOLAR}"
)
workflow_args+=(
    --protein-lipid-cutoff "${PROTEIN_LIPID_CUTOFF}"
    --ion-cutoff "${ION_CUTOFF}"
    --xy-scale "${XY_SCALE}"
    --box-padding-xy "${BOX_PADDING_XY}"
    --box-padding-z "${BOX_PADDING_Z}"
    --protein-placement-mode "${PROTEIN_PLACEMENT_MODE}"
    --protein-orientation-mode "${PROTEIN_ORIENTATION_MODE}"
    --protein-surface-gap "${PROTEIN_SURFACE_GAP}"
    --protein-lipid-min-gap "${PROTEIN_LIPID_MIN_GAP}"
    --protein-lipid-cutoff-step "${PROTEIN_LIPID_CUTOFF_STEP}"
    --protein-lipid-cutoff-max "${PROTEIN_LIPID_CUTOFF_MAX}"
    --temperature "${TEMPERATURE}"
    --thermostat-timescale "${THERMOSTAT_TIMESCALE}"
    --thermostat-interval "${THERMOSTAT_INTERVAL}"
    --strict-stage-handoff "${STRICT_STAGE_HANDOFF}"
    --min-60-max-iter "${MIN_60_MAX_ITER}"
    --min-61-max-iter "${MIN_61_MAX_ITER}"
    --min-70-max-iter "${MIN_70_MAX_ITER}"
    --eq-60-nsteps "${EQ_60_NSTEPS}"
    --eq-62-nsteps "${EQ_62_NSTEPS}"
    --eq-63-nsteps "${EQ_63_NSTEPS}"
    --eq-64-nsteps "${EQ_64_NSTEPS}"
    --eq-65-nsteps "${EQ_65_NSTEPS}"
    --eq-66-nsteps "${EQ_66_NSTEPS}"
    --prod-70-burnin-nsteps "${PROD_70_BURNIN_NSTEPS}"
    --prod-70-nsteps "${PROD_70_NSTEPS}"
    --eq-time-step "${EQ_TIME_STEP}"
    --prod-time-step "${PROD_TIME_STEP}"
    --eq-frame-steps "${EQ_FRAME_STEPS}"
    --prod-frame-steps "${PROD_FRAME_STEPS}"
    --prod-70-npt-enable "${PROD_70_NPT_ENABLE}"
    --prep-seed "${PREP_SEED}"
    --seed "${SEED}"
    --continue-stage-70-from "${CONTINUE_STAGE_70_FROM}"
    --continue-stage-70-output "${CONTINUE_STAGE_70_OUTPUT}"
    --continue-stage-70-label "${CONTINUE_STAGE_70_LABEL}"
)

python3 "${UNIVERSAL_PREP_SCRIPT}" "${workflow_args[@]}"
