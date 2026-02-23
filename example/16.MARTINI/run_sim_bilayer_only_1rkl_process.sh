#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../source.sh"
source "${SCRIPT_DIR}/../../.venv/bin/activate"
set -euo pipefail
cd "${SCRIPT_DIR}"

# Bilayer-only workflow with stage/process parity to run_sim_1rkl.sh:
# 0) Bilayer-only preparation (tile/crop bilayer + ion placement, no protein)
# 1) Stage input generation (dry MARTINI)
# 2) Run 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0
# 3) Extract per-stage VTF trajectories

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

PDB_ID="${PDB_ID:-bilayer}"
RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID}_bilayer_only}"

INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="${RUN_DIR:-outputs/martini_test_bilayer_only_1rkl_process}"
CHECKPOINT_DIR="${CHECKPOINT_DIR:-${RUN_DIR}/checkpoints}"
LOG_DIR="${LOG_DIR:-${RUN_DIR}/logs}"

BILAYER_PDB="${BILAYER_PDB:-pdb/bilayer.MARTINI.pdb}"
RUNTIME_PDB_FILE="${SCRIPT_DIR}/pdb/${RUNTIME_PDB_ID}.MARTINI.pdb"
BILAYER_PREP_SCRIPT="${BILAYER_PREP_SCRIPT:-${SCRIPT_DIR}/prepare_system.py}"
BILAYER_PREP_ENABLE="${BILAYER_PREP_ENABLE:-1}"
BILAYER_PREP_DIR="${BILAYER_PREP_DIR:-${RUN_DIR}/bilayer_prep}"

# Unified prep controls
BILAYER_PREP_MODE="${BILAYER_PREP_MODE:-bilayer}"
# Bilayer prep controls (same prep-process lineage as run_sim_1rkl stage-0 packing)
BILAYER_XY_SCALE="${BILAYER_XY_SCALE:-2.0}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"
ION_CUTOFF="${ION_CUTOFF:-4.0}"
BOX_PADDING_Z="${BOX_PADDING_Z:-20.0}"
PREP_SEED="${PREP_SEED:-2026}"
PROTEIN_CG_PDB="${PROTEIN_CG_PDB:-}"
PROTEIN_ITP="${PROTEIN_ITP:-}"
PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-3.0}"
PROTEIN_NET_CHARGE="${PROTEIN_NET_CHARGE:-}"

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

# Simulation controls
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"

MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-500}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-500}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-500}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-500}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-500}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-500}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-500}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-100000}"

EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.002}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-100000}"

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

if [ ! -f "${BILAYER_PDB}" ]; then
    echo "ERROR: bilayer PDB not found: ${BILAYER_PDB}"
    exit 1
fi

if [ "${BILAYER_PREP_ENABLE}" = "1" ]; then
    if [ ! -f "${BILAYER_PREP_SCRIPT}" ]; then
        echo "ERROR: bilayer prep script not found: ${BILAYER_PREP_SCRIPT}"
        exit 1
    fi
fi

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR" "$BILAYER_PREP_DIR"

prepare_runtime_pdb() {
    if [ "${BILAYER_PREP_ENABLE}" = "1" ]; then
        echo "=== Stage 0: Bilayer Expansion Prep ==="
        local prep_cmd=(
            python3 "${BILAYER_PREP_SCRIPT}"
            --mode "${BILAYER_PREP_MODE}" \
            --pdb-id "${RUNTIME_PDB_ID}" \
            --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
            --prepare-structure 1 \
            --bilayer-pdb "${BILAYER_PDB}" \
            --xy-scale "${BILAYER_XY_SCALE}" \
            --box-padding-z "${BOX_PADDING_Z}" \
            --salt-molar "${SALT_MOLAR}" \
            --ion-cutoff "${ION_CUTOFF}" \
            --seed "${PREP_SEED}" \
            --summary-json "${BILAYER_PREP_DIR}/structure_prep_summary.json"
        )
        if [ -n "${PROTEIN_CG_PDB}" ]; then
            prep_cmd+=(--protein-cg-pdb "${PROTEIN_CG_PDB}")
        fi
        if [ -n "${PROTEIN_ITP}" ]; then
            prep_cmd+=(--protein-itp "${PROTEIN_ITP}")
        fi
        if [ -n "${PROTEIN_LIPID_CUTOFF}" ]; then
            prep_cmd+=(--protein-lipid-cutoff "${PROTEIN_LIPID_CUTOFF}")
        fi
        if [ -n "${PROTEIN_NET_CHARGE}" ]; then
            prep_cmd+=(--protein-net-charge "${PROTEIN_NET_CHARGE}")
        fi
        "${prep_cmd[@]}"

        if [ ! -f "${RUNTIME_PDB_FILE}" ]; then
            echo "ERROR: runtime PDB not found after prep: ${RUNTIME_PDB_FILE}"
            exit 1
        fi

        echo "Runtime PDB prepared from enlarged bilayer: ${RUNTIME_PDB_FILE}"
    else
        cp -f "${BILAYER_PDB}" "${RUNTIME_PDB_FILE}"
        echo "Runtime PDB copied directly from template bilayer: ${RUNTIME_PDB_FILE}"
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

prepare_stage_file() {
    local target_file="$1"
    local prepare_stage="$2"
    local npt_enable="$3"
    local barostat_type="$4"
    local lipidhead_fc="${5:-0}"
    local stage_label="${6:-minimization}"

    # Retain this parameter for call-site parity with run_sim_1rkl.sh.
    : "${stage_label}"

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"
    export UPSIDE_BILAYER_LIPIDHEAD_FC="$lipidhead_fc"

    python3 "${BILAYER_PREP_SCRIPT}" \
        --mode "${BILAYER_PREP_MODE}" \
        --pdb-id "${RUNTIME_PDB_ID}" \
        --runtime-pdb-output "${RUNTIME_PDB_FILE}" \
        --prepare-structure 0 \
        --stage "$prepare_stage" \
        --run-dir "$RUN_DIR" \
        --summary-json "${BILAYER_PREP_DIR}/stage_${prepare_stage}.summary.json"

    local prepared_tmp="${RUN_DIR}/test.input.up"
    if [ ! -f "$prepared_tmp" ]; then
        echo "ERROR: preparation failed for stage ${prepare_stage}: $prepared_tmp not found"
        exit 1
    fi

    mv -f "$prepared_tmp" "$target_file"

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
            # Semi-isotropic coupling keeps z compressibility at 0.0.
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

echo "=== Bilayer-Only Dry MARTINI Workflow (run_sim_1rkl process) ==="
echo "Bilayer input: ${BILAYER_PDB}"
echo "Prep script: ${BILAYER_PREP_SCRIPT} (mode=${BILAYER_PREP_MODE})"
echo "Bilayer prep: ${BILAYER_PREP_DIR} (enabled=${BILAYER_PREP_ENABLE}, xy_scale=${BILAYER_XY_SCALE})"
echo "PDB ID: $PDB_ID"
echo "Runtime PDB ID: $RUNTIME_PDB_ID"
echo "Simulation stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
echo

prepare_runtime_pdb

# 6.0: soft-core minimization
set_stage_npt_targets "6.0"
prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0" "0" "minimization"
cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"
extract_stage_vtf "6.0" "$STAGE_60_FILE" "1"

# 6.1: hard minimization
set_stage_npt_targets "6.1"
prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0" "0" "minimization"
cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
python3 set_initial_position.py "$STAGE_60_FILE" "$STAGE_61_FILE"
run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"
extract_stage_vtf "6.1" "$STAGE_61_FILE" "1"

# 6.2: soft equilibration
set_stage_npt_targets "6.2"
prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "0" "200" "minimization"
python3 set_initial_position.py "$STAGE_61_FILE" "$STAGE_62_FILE"
run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.2" "$STAGE_62_FILE" "1"

# 6.3: reduced softening equilibration
set_stage_npt_targets "6.3"
prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "0" "100" "minimization"
cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
python3 set_initial_position.py "$STAGE_62_FILE" "$STAGE_63_FILE"
run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.3" "$STAGE_63_FILE" "1"

# 6.4-6.6: hard equilibration with restraint ramp
set_stage_npt_targets "6.4"
prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "0" "50" "minimization"
cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
python3 set_initial_position.py "$STAGE_63_FILE" "$STAGE_64_FILE"
run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.4" "$STAGE_64_FILE" "1"

set_stage_npt_targets "6.5"
prepare_stage_file "$PREPARED_65_FILE" "npt_prod" "1" "0" "20" "minimization"
cp -f "$PREPARED_65_FILE" "$STAGE_65_FILE"
python3 set_initial_position.py "$STAGE_64_FILE" "$STAGE_65_FILE"
run_md_stage "6.5" "$STAGE_65_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.5" "$STAGE_65_FILE" "1"

set_stage_npt_targets "6.6"
prepare_stage_file "$PREPARED_66_FILE" "npt_prod" "1" "0" "10" "minimization"
cp -f "$PREPARED_66_FILE" "$STAGE_66_FILE"
python3 set_initial_position.py "$STAGE_65_FILE" "$STAGE_66_FILE"
run_md_stage "6.6" "$STAGE_66_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
extract_stage_vtf "6.6" "$STAGE_66_FILE" "1"

# 7.0: production
prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE" "0" "production"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
python3 set_initial_position.py "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"
extract_stage_vtf "7.0" "$STAGE_70_FILE" "1"

echo
echo "=== Workflow Complete ==="
echo "Runtime PDB:"
echo "  ${RUNTIME_PDB_FILE}"
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
