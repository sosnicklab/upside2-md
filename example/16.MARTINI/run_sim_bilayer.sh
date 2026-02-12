#!/bin/bash
source ../../source.sh
source ../../.venv/bin/activate

# Dry MARTINI bilayer workflow aligned to CHARMM-GUI Gromacs stages:
# 6.0 soft-core minimization
# 6.1 hard minimization
# 6.2-6.6 hard equilibration (semi-isotropic Berendsen)
# 7.0 hard production (configurable: Parrinello-Rahman by default)

set -euo pipefail

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

INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"
LOG_DIR="${RUN_DIR}/logs"

PREPARED_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.prepared.up"
STAGE_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
PREPARED_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.prepared.up"
STAGE_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
STAGE_62_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.2.up"
PREPARED_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.prepared.up"
STAGE_63_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.up"
PREPARED_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.prepared.up"
STAGE_64_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.up"
STAGE_65_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.up"
STAGE_66_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up"
PREPARED_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.prepared.up"
STAGE_70_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up"

# Temperature mapping: 303.15 K / 350.588235 K_per_T_up ~= 0.8647
TEMPERATURE="${TEMPERATURE:-0.8647}"
THERMOSTAT_TIMESCALE="${THERMOSTAT_TIMESCALE:-4.0}"
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED="${SEED:-7090685331}"

# Gromacs nsteps from dry MARTINI mdp files
MIN_60_MAX_ITER="${MIN_60_MAX_ITER:-5000}"
MIN_61_MAX_ITER="${MIN_61_MAX_ITER:-5000}"
EQ_62_NSTEPS="${EQ_62_NSTEPS:-50000}"
EQ_63_NSTEPS="${EQ_63_NSTEPS:-50000}"
EQ_64_NSTEPS="${EQ_64_NSTEPS:-50000}"
EQ_65_NSTEPS="${EQ_65_NSTEPS:-50000}"
EQ_66_NSTEPS="${EQ_66_NSTEPS:-50000}"
PROD_70_NSTEPS="${PROD_70_NSTEPS:-1500000}"

# dt from dry MARTINI mdp
EQ_TIME_STEP="${EQ_TIME_STEP:-0.010}"
PROD_TIME_STEP="${PROD_TIME_STEP:-0.020}"
MIN_TIME_STEP="${MIN_TIME_STEP:-0.010}"

# output cadence in integration steps (mdp uses 1000 for equil, 5000 for production)
EQ_FRAME_STEPS="${EQ_FRAME_STEPS:-1000}"
PROD_FRAME_STEPS="${PROD_FRAME_STEPS:-5000}"

# Unit conversions from AGENTS.md / CLAUDE.md
BAR_TO_EUP="${BAR_TO_EUP:-0.000020659477}"
COMP_3E4_BAR_INV_TO_A3_PER_EUP="${COMP_3E4_BAR_INV_TO_A3_PER_EUP:-14.521180763676}"

# NPT defaults for dry MARTINI implicit bilayer:
# semi-isotropic, tensionless membrane (Pxy=0), fixed normal axis (beta_z=0).
# With beta_z=0, the z reference pressure is inert; keep ref Pz=0 to match paper/Gromacs mdp.
export UPSIDE_NPT_TARGET_PXY="${UPSIDE_NPT_TARGET_PXY:-0.0}"
export UPSIDE_NPT_TARGET_PZ="${UPSIDE_NPT_TARGET_PZ:-0.0}"
export UPSIDE_NPT_TAU="${UPSIDE_NPT_TAU:-4.0}"
export UPSIDE_NPT_COMPRESSIBILITY="${UPSIDE_NPT_COMPRESSIBILITY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_XY="${UPSIDE_NPT_COMPRESSIBILITY_XY:-$COMP_3E4_BAR_INV_TO_A3_PER_EUP}"
export UPSIDE_NPT_COMPRESSIBILITY_Z="${UPSIDE_NPT_COMPRESSIBILITY_Z:-0.0}"
export UPSIDE_NPT_INTERVAL="${UPSIDE_NPT_INTERVAL:-10}"
export UPSIDE_NPT_SEMI="${UPSIDE_NPT_SEMI:-1}"
export UPSIDE_NPT_DEBUG="${UPSIDE_NPT_DEBUG:-1}"
PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-1}"
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

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

duration_from_nsteps() {
    local nsteps="$1"
    local dt="$2"
    awk -v n="$nsteps" -v d="$dt" 'BEGIN { printf "%.8f", n*d }'
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

    export UPSIDE_SIMULATION_STAGE="$prepare_stage"
    export UPSIDE_NPT_ENABLE="$npt_enable"

    python3 prepare_martini.py "$PDB_ID" --stage "$prepare_stage" "$RUN_DIR"

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

    local duration
    duration="$(duration_from_nsteps "$nsteps" "$dt")"
    local frame_interval
    frame_interval="$(duration_from_nsteps "$frame_steps" "$dt")"

    if [ "$input_file" != "$output_file" ]; then
        cp -f "$input_file" "$output_file"
        python3 set_initial_position.py "$input_file" "$output_file"
    fi

    local log_file="${LOG_DIR}/stage_${stage_label}.log"
    echo "=== Stage ${stage_label}: MD ==="
    echo "Input:  $input_file"
    echo "Output: $output_file"
    echo "nsteps=${nsteps}, dt=${dt}, duration=${duration}"

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$output_file"
        "--duration" "$duration"
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

echo "=== Dry MARTINI Bilayer Workflow ==="
echo "PDB ID: $PDB_ID"
echo "Sequence: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
echo "WARNING: Gromacs BILAYER_LIPIDHEAD_FC restraint ramp (200/100/50/20/10) is not yet implemented in prepare_martini.py; this script matches stage order, timesteps, and pressure-coupling behavior."
echo

# 6.0: soft-core minimization, NPT on, Berendsen
prepare_stage_file "$PREPARED_60_FILE" "minimization" "1" "0"
cp -f "$PREPARED_60_FILE" "$STAGE_60_FILE"
run_minimization_stage "6.0" "$STAGE_60_FILE" "$MIN_60_MAX_ITER"

# 6.1: hard minimization, NPT on, Berendsen
prepare_stage_file "$PREPARED_61_FILE" "npt_prod" "1" "0"
cp -f "$PREPARED_61_FILE" "$STAGE_61_FILE"
python3 set_initial_position.py "$STAGE_60_FILE" "$STAGE_61_FILE"
run_minimization_stage "6.1" "$STAGE_61_FILE" "$MIN_61_MAX_ITER"

# 6.2-6.6: hard equilibration, NPT on, Berendsen
# 6.2: soft potential (smoother start, compensates for missing lipid-head restraint ramp)
prepare_stage_file "$STAGE_62_FILE" "npt_equil" "1" "0"
python3 set_initial_position.py "$STAGE_61_FILE" "$STAGE_62_FILE"
run_md_stage "6.2" "$STAGE_62_FILE" "$STAGE_62_FILE" "$EQ_62_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"

# 6.3: reduced softening
prepare_stage_file "$PREPARED_63_FILE" "npt_equil_reduced" "1" "0"
cp -f "$PREPARED_63_FILE" "$STAGE_63_FILE"
python3 set_initial_position.py "$STAGE_62_FILE" "$STAGE_63_FILE"
run_md_stage "6.3" "$STAGE_63_FILE" "$STAGE_63_FILE" "$EQ_63_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"

# 6.4-6.6: hard potential with Berendsen
prepare_stage_file "$PREPARED_64_FILE" "npt_prod" "1" "0"
cp -f "$PREPARED_64_FILE" "$STAGE_64_FILE"
python3 set_initial_position.py "$STAGE_63_FILE" "$STAGE_64_FILE"
run_md_stage "6.4" "$STAGE_64_FILE" "$STAGE_64_FILE" "$EQ_64_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
run_md_stage "6.5" "$STAGE_64_FILE" "$STAGE_65_FILE" "$EQ_65_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"
run_md_stage "6.6" "$STAGE_65_FILE" "$STAGE_66_FILE" "$EQ_66_NSTEPS" "$EQ_TIME_STEP" "$EQ_FRAME_STEPS"

# 7.0: hard production (default Parrinello-Rahman, set PROD_70_NPT_ENABLE=0 for NVT)
prepare_stage_file "$PREPARED_70_FILE" "npt_prod" "$PROD_70_NPT_ENABLE" "$PROD_70_BAROSTAT_TYPE"
cp -f "$PREPARED_70_FILE" "$STAGE_70_FILE"
python3 set_initial_position.py "$STAGE_66_FILE" "$STAGE_70_FILE"
run_md_stage "7.0" "$STAGE_70_FILE" "$STAGE_70_FILE" "$PROD_70_NSTEPS" "$PROD_TIME_STEP" "$PROD_FRAME_STEPS"

# VTF outputs for key checkpoints
VTF_61_FILE="${RUN_DIR}/${PDB_ID}.stage_6.1.vtf"
VTF_66_FILE="${RUN_DIR}/${PDB_ID}.stage_6.6.vtf"
VTF_70_FILE="${RUN_DIR}/${PDB_ID}.stage_7.0.vtf"

python3 extract_martini_vtf.py "$STAGE_61_FILE" "$VTF_61_FILE" "$STAGE_61_FILE" "$PDB_ID"
python3 extract_martini_vtf.py "$STAGE_66_FILE" "$VTF_66_FILE" "$STAGE_66_FILE" "$PDB_ID"
python3 extract_martini_vtf.py "$STAGE_70_FILE" "$VTF_70_FILE" "$STAGE_70_FILE" "$PDB_ID"

echo
echo "=== Workflow Complete ==="
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
echo "  $VTF_61_FILE"
echo "  $VTF_66_FILE"
echo "  $VTF_70_FILE"
