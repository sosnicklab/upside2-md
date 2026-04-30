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
RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_hybrid}"
RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_hybrid}"

PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
PROTEIN_CG_PDB="${PROTEIN_CG_PDB:-}"
PROTEIN_ITP="${PROTEIN_ITP:-}"
BILAYER_PDB="${BILAYER_PDB:-${UPSIDE_HOME}/parameters/dryMARTINI/DOPC.pdb}"
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${PROJECT_ROOT}/py/martini_prepare_system.py}"

MARTINIZE_ENABLE="${MARTINIZE_ENABLE:-1}"
MARTINIZE_FF="${MARTINIZE_FF:-martini22}"
MARTINIZE_MOLNAME="${MARTINIZE_MOLNAME:-PROA}"
MARTINIZE_SCRIPT="${MARTINIZE_SCRIPT:-${PROJECT_ROOT}/py/martini_martinize.py}"
EXTRACT_VTF_SCRIPT="${EXTRACT_VTF_SCRIPT:-${PROJECT_ROOT}/py/martini_extract_vtf.py}"

SALT_MOLAR="${SALT_MOLAR:-0.15}"
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-4.5}"
ION_CUTOFF="${ION_CUTOFF:-4.0}"
XY_SCALE="${XY_SCALE:-1.0}"
BOX_PADDING_XY="${BOX_PADDING_XY:-0.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-embed}"
PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-input}"
PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-4.5}"
PROTEIN_LIPID_CUTOFF_STEP="${PROTEIN_LIPID_CUTOFF_STEP:-0.5}"
PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-8.0}"

TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-5.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
STRICT_STAGE_HANDOFF="${STRICT_STAGE_HANDOFF:-1}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-10000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"
EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-50}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"
PROD_70_BACKBONE_FIX_RIGID_ENABLE="${PROD_70_BACKBONE_FIX_RIGID_ENABLE:-0}"

PREP_SEED="${PREP_SEED:-}"
SEED="${SEED:-}"
CONTINUE_STAGE_70_FROM="${CONTINUE_STAGE_70_FROM:-}"
CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-}"
CONTINUE_STAGE_70_LABEL="${CONTINUE_STAGE_70_LABEL:-}"
PREVIOUS_RUN_DIR="${PREVIOUS_RUN_DIR:-}"
PREVIOUS_STAGE7_FILE="${PREVIOUS_STAGE7_FILE:-}"
AUTO_CONTINUE_FROM_PREVIOUS_RUN="${AUTO_CONTINUE_FROM_PREVIOUS_RUN:-0}"
AUTO_CONTINUE_GLOB="${AUTO_CONTINUE_GLOB:-}"

if [ ! -f "${UNIVERSAL_PREP_SCRIPT}" ]; then
    echo "ERROR: preparation script not found: ${UNIVERSAL_PREP_SCRIPT}" >&2
    exit 1
fi

find_latest_stage7_in_pattern() {
    local pattern="$1"
    local best_file=""
    local best_idx=-1
    local path=""
    local file=""
    local idx=""
    shopt -s nullglob
    for path in $pattern; do
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

if [ -z "${CONTINUE_STAGE_70_FROM}" ] && [ -n "${PREVIOUS_STAGE7_FILE}" ] && [ -f "${PREVIOUS_STAGE7_FILE}" ]; then
    CONTINUE_STAGE_70_FROM="${PREVIOUS_STAGE7_FILE}"
fi

if [ -z "${CONTINUE_STAGE_70_FROM}" ] && [ -n "${PREVIOUS_RUN_DIR}" ]; then
    CONTINUE_STAGE_70_FROM="$(find_latest_stage7_in_pattern "${PREVIOUS_RUN_DIR}/checkpoints/${PDB_ID}.stage_7.*.up")"
fi

if [ -z "${CONTINUE_STAGE_70_FROM}" ]; then
    CONTINUE_STAGE_70_FROM="$(find_latest_stage7_in_pattern "${RUN_DIR}/checkpoints/${PDB_ID}.stage_7.*.up")"
fi

if [ -z "${CONTINUE_STAGE_70_FROM}" ] && [ "${AUTO_CONTINUE_FROM_PREVIOUS_RUN}" = "1" ] && [ -n "${AUTO_CONTINUE_GLOB}" ]; then
    CONTINUE_STAGE_70_FROM="$(find_latest_stage7_in_pattern "${SCRIPT_DIR}/outputs/${AUTO_CONTINUE_GLOB}")"
fi

if [ -n "${CONTINUE_STAGE_70_FROM}" ]; then
    echo "Detected continuation source: ${CONTINUE_STAGE_70_FROM}"
fi

python3 "${UNIVERSAL_PREP_SCRIPT}" run-hybrid-workflow \
    --pdb-id "${PDB_ID}" \
    --runtime-pdb-id "${RUNTIME_PDB_ID}" \
    --upside-home "${UPSIDE_HOME}" \
    --run-dir "${RUN_DIR}" \
    --protein-aa-pdb "${PROTEIN_AA_PDB}" \
    --protein-cg-pdb "${PROTEIN_CG_PDB}" \
    --protein-itp "${PROTEIN_ITP}" \
    --bilayer-pdb "${BILAYER_PDB}" \
    --martinize-enable "${MARTINIZE_ENABLE}" \
    --martinize-ff "${MARTINIZE_FF}" \
    --martinize-molname "${MARTINIZE_MOLNAME}" \
    --martinize-script "${MARTINIZE_SCRIPT}" \
    --extract-vtf-script "${EXTRACT_VTF_SCRIPT}" \
    --salt-molar "${SALT_MOLAR}" \
    --protein-lipid-cutoff "${PROTEIN_LIPID_CUTOFF}" \
    --ion-cutoff "${ION_CUTOFF}" \
    --xy-scale "${XY_SCALE}" \
    --box-padding-xy "${BOX_PADDING_XY}" \
    --box-padding-z "${BOX_PADDING_Z}" \
    --protein-placement-mode "${PROTEIN_PLACEMENT_MODE}" \
    --protein-orientation-mode "${PROTEIN_ORIENTATION_MODE}" \
    --protein-surface-gap "${PROTEIN_SURFACE_GAP}" \
    --protein-lipid-min-gap "${PROTEIN_LIPID_MIN_GAP}" \
    --protein-lipid-cutoff-step "${PROTEIN_LIPID_CUTOFF_STEP}" \
    --protein-lipid-cutoff-max "${PROTEIN_LIPID_CUTOFF_MAX}" \
    --temperature "${TEMPERATURE}" \
    --thermostat-timescale "${THERMOSTAT_TIMESCALE}" \
    --thermostat-interval "${THERMOSTAT_INTERVAL}" \
    --strict-stage-handoff "${STRICT_STAGE_HANDOFF}" \
    --min-60-max-iter "${MIN_60_MAX_ITER}" \
    --min-61-max-iter "${MIN_61_MAX_ITER}" \
    --eq-62-nsteps "${EQ_62_NSTEPS}" \
    --eq-63-nsteps "${EQ_63_NSTEPS}" \
    --eq-64-nsteps "${EQ_64_NSTEPS}" \
    --eq-65-nsteps "${EQ_65_NSTEPS}" \
    --eq-66-nsteps "${EQ_66_NSTEPS}" \
    --prod-70-nsteps "${PROD_70_NSTEPS}" \
    --eq-time-step "${EQ_TIME_STEP}" \
    --prod-time-step "${PROD_TIME_STEP}" \
    --min-time-step "${MIN_TIME_STEP}" \
    --eq-frame-steps "${EQ_FRAME_STEPS}" \
    --prod-frame-steps "${PROD_FRAME_STEPS}" \
    --prod-70-npt-enable "${PROD_70_NPT_ENABLE}" \
    --prod-70-backbone-fix-rigid-enable "${PROD_70_BACKBONE_FIX_RIGID_ENABLE}" \
    --prep-seed "${PREP_SEED}" \
    --seed "${SEED}" \
    --continue-stage-70-from "${CONTINUE_STAGE_70_FROM}" \
    --continue-stage-70-output "${CONTINUE_STAGE_70_OUTPUT}" \
    --continue-stage-70-label "${CONTINUE_STAGE_70_LABEL}" \
    --previous-run-dir "" \
    --previous-stage7-file "" \
    --auto-continue-from-previous-run "0" \
    --auto-continue-glob ""
