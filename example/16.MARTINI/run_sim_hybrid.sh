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
LIPID_RESOLUTION="${LIPID_RESOLUTION:-coarse}"
if [ "${LIPID_RESOLUTION}" != "coarse" ] && [ "${LIPID_RESOLUTION}" != "full" ]; then
    echo "ERROR: LIPID_RESOLUTION must be coarse or full" >&2
    exit 1
fi

coarse_hybrid_args=()
if [ "${LIPID_RESOLUTION}" = "coarse" ]; then
    export UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE="${UPSIDE_MARTINI_TEMPERED_AVERAGE_TEMP_UPSIDE:-25.0}"
    export EXPLICIT_IONS="${EXPLICIT_IONS:-1}"
    export CG_LIPID_MASS_SCALE="${CG_LIPID_MASS_SCALE:-0.012}"
    export CG_LIPID_ROTATIONAL_THERMOSTAT_TIMESCALE="${CG_LIPID_ROTATIONAL_THERMOSTAT_TIMESCALE:-0.008}"
    export UPSIDE_CGL_COMPACTION_IMPLICIT_RESPONSE="${UPSIDE_CGL_COMPACTION_IMPLICIT_RESPONSE:-0}"
    export CGL_GLE_MEMORY_TAUS="${CGL_GLE_MEMORY_TAUS:-0.2,2.0}"
    export CGL_GLE_COUPLINGS="${CGL_GLE_COUPLINGS:-0.30375,0.2205}"
    export CGL_GLE_TEMPERATURE_GRID="${CGL_GLE_TEMPERATURE_GRID:-0.7,0.8,0.8647,0.9,1.0,1.1,1.2}"
    export CGL_GLE_COUPLING_SCALES="${CGL_GLE_COUPLING_SCALES:-1,1;1,1;1,1;1.013,1.057;1.05,1.22;1.20,1.39;1.50,1.50}"
    export CGL_GLE_MEMORY_TAU_SCALES="${CGL_GLE_MEMORY_TAU_SCALES:-0.33,0.33;0.50,0.50;1,1;0.85,0.95;1.10,1.33;1.35,1.58;1.70,1.70}"
    export CGL_GLE_REPLACE_MARKOVIAN="${CGL_GLE_REPLACE_MARKOVIAN:-1}"
    export PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-1}"
    export CGL_VTF_DISPLAY_MODE="${CGL_VTF_DISPLAY_MODE:-rod}"
    coarse_hybrid_args=(
        --explicit-ions "${EXPLICIT_IONS}"
        --cg-lipid-mass-scale "${CG_LIPID_MASS_SCALE}"
        --cgl-gle-memory-taus "${CGL_GLE_MEMORY_TAUS}"
        --cgl-gle-couplings "${CGL_GLE_COUPLINGS}"
        --cgl-gle-replace-markovian "${CGL_GLE_REPLACE_MARKOVIAN}"
    )
    if [ -n "${CGL_GLE_TEMPERATURE_GRID}" ] || [ -n "${CGL_GLE_COUPLING_SCALES}" ] || [ -n "${CGL_GLE_MEMORY_TAU_SCALES}" ]; then
        coarse_hybrid_args+=(
            --cgl-gle-temperature-grid "${CGL_GLE_TEMPERATURE_GRID}"
            --cgl-gle-coupling-scales "${CGL_GLE_COUPLING_SCALES}"
            --cgl-gle-memory-tau-scales "${CGL_GLE_MEMORY_TAU_SCALES}"
        )
    fi
fi

if [ "${LIPID_RESOLUTION}" = "full" ]; then
    default_suffix="_hybrid_full"
else
    default_suffix="_hybrid"
fi

RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}${default_suffix}}"
RUN_DIR="${RUN_DIR:-outputs/martini_${PDB_ID}${default_suffix}}"
PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/${PDB_ID}.pdb}"
BILAYER_PDB="${BILAYER_PDB:-${UPSIDE_HOME}/parameters/dryMARTINI/DOPC.pdb}"
UNIVERSAL_PREP_SCRIPT="${UNIVERSAL_PREP_SCRIPT:-${PROJECT_ROOT}/py/martini_prepare_system.py}"
EXTRACT_VTF_SCRIPT="${EXTRACT_VTF_SCRIPT:-${PROJECT_ROOT}/py/martini_extract_vtf.py}"
CGL_VTF_DISPLAY_MODE="${CGL_VTF_DISPLAY_MODE:-rod}"

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

TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-5.0}"
CG_LIPID_THERMOSTAT_TIMESCALE="${CG_LIPID_THERMOSTAT_TIMESCALE:-0.0}"
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
PROD_70_NSTEPS="${PROD_70_NSTEPS:-10000}"
STAGE_70_BURNIN_PROTEIN_RESTRAINT_SPRING="${STAGE_70_BURNIN_PROTEIN_RESTRAINT_SPRING:-10.0}"
STAGE_70_RELEASE_SC_ENV_BACKBONE_HOLD_STEPS="${STAGE_70_RELEASE_SC_ENV_BACKBONE_HOLD_STEPS:-0}"
STAGE_70_RELEASE_SC_ENV_PO4_Z_HOLD_STEPS="${STAGE_70_RELEASE_SC_ENV_PO4_Z_HOLD_STEPS:-0}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-50}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"

PREP_SEED="${PREP_SEED:-}"
SEED="${SEED:-}"
CONTINUE_STAGE_70_FROM="${CONTINUE_STAGE_70_FROM:-}"
CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-}"
CONTINUE_STAGE_70_LABEL="${CONTINUE_STAGE_70_LABEL:-}"

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
    local martini_ff_dir="${UPSIDE_HOME}/parameters/dryMARTINI"
    local param_file=""
    local required_params=(
        "${martini_ff_dir}/particle.h5"
        "${martini_ff_dir}/sidechain.h5"
    )
    if [ "${LIPID_RESOLUTION}" = "coarse" ]; then
        required_params+=("${martini_ff_dir}/dopc.h5")
    fi
    for param_file in "${required_params[@]}"; do
        if [ ! -f "${param_file}" ]; then
            echo "One or more MARTINI parameter files missing. Generating..."
            python3 "${PROJECT_ROOT}/py/martini_gen_params.py" --upside-home "${UPSIDE_HOME}"
            break
        fi
    done
}

if [ -z "${CONTINUE_STAGE_70_FROM}" ]; then
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
    --cg-lipid-vtf-display "${CGL_VTF_DISPLAY_MODE}"
    --salt-molar "${SALT_MOLAR}"
)
if [ "${LIPID_RESOLUTION}" = "coarse" ]; then
    workflow_args+=("${coarse_hybrid_args[@]}")
fi
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
    --cg-lipid-thermostat-timescale "${CG_LIPID_THERMOSTAT_TIMESCALE}"
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
    --stage-70-burnin-protein-restraint-spring "${STAGE_70_BURNIN_PROTEIN_RESTRAINT_SPRING}"
    --stage-70-release-sc-env-backbone-hold-steps "${STAGE_70_RELEASE_SC_ENV_BACKBONE_HOLD_STEPS}"
    --stage-70-release-sc-env-po4-z-hold-steps "${STAGE_70_RELEASE_SC_ENV_PO4_Z_HOLD_STEPS}"
    --eq-time-step "${EQ_TIME_STEP}"
    --prod-time-step "${PROD_TIME_STEP}"
    --eq-frame-steps "${EQ_FRAME_STEPS}"
    --prod-frame-steps "${PROD_FRAME_STEPS}"
    --prod-70-npt-enable "${PROD_70_NPT_ENABLE}"
    --lipid-resolution "${LIPID_RESOLUTION}"
    --prep-seed "${PREP_SEED}"
    --seed "${SEED}"
    --continue-stage-70-from "${CONTINUE_STAGE_70_FROM}"
    --continue-stage-70-output "${CONTINUE_STAGE_70_OUTPUT}"
    --continue-stage-70-label "${CONTINUE_STAGE_70_LABEL}"
)

python3 "${UNIVERSAL_PREP_SCRIPT}" "${workflow_args[@]}"
