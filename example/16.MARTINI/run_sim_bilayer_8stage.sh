#!/bin/bash

# MARTINI 3.0 8-Stage Bilayer Equilibration Workflow
# Implements CHARMM-GUI MARTINI protocol (Steps 6.0-7.0)
#
# Stage 6.0: Minimization with soft-core potentials
# Stage 6.1: Minimization without soft-core
# Stage 6.2: MD dt=0.002, k_lipid=200, Berendsen
# Stage 6.3: MD dt=0.005, k_lipid=100, Berendsen
# Stage 6.4: MD dt=0.010, k_lipid=50, Berendsen
# Stage 6.5: MD dt=0.015, k_lipid=25, Berendsen
# Stage 6.6: MD dt=0.020, k_lipid=10, Berendsen
# Stage 7.0: MD Production dt=0.020, k=0, Parrinello-Rahman

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PDB_ID="${PDB_ID:-bilayer}"

# Directories
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_8stage"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"
LOG_DIR="${RUN_DIR}/logs"

# Create directories
mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

# Simulation parameters (UPSIDE units)
TEMPERATURE=0.8  # ~303K in UPSIDE units
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED=12345

# Stage durations (in MD steps)
MIN_STEPS_60=500
MIN_STEPS_61=500
MD_STEPS_62=2000
MD_STEPS_63=2000
MD_STEPS_64=2000
MD_STEPS_65=2000
MD_STEPS_66=2000
MD_STEPS_70=5000

FRAME_INTERVAL=20

# NPT settings (UPSIDE units: 1 bar = 0.000020659 E_up/Angstrom^3)
export UPSIDE_NPT_TARGET_PXY=0.000020659
export UPSIDE_NPT_TARGET_PZ=0.000020659
export UPSIDE_NPT_TAU=1.0
export UPSIDE_NPT_COMPRESSIBILITY=3e-4
export UPSIDE_NPT_INTERVAL=10
export UPSIDE_NPT_SEMI=1
export UPSIDE_NPT_DEBUG=1

# =============================================================================
# VALIDATION
# =============================================================================
if [ -z "$UPSIDE_HOME" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set!"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi

echo "=== MARTINI 3.0 8-Stage Bilayer Equilibration Workflow ==="
echo "PDB ID: $PDB_ID"
echo "Stages: 6.0 -> 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5 -> 6.6 -> 7.0"
echo

# =============================================================================
# STAGE 0: PREPARE INPUT FILES
# =============================================================================
echo "=== Stage 0: Preparing Input Files ==="
source ../../source.sh
source ../../.venv/bin/activate

export UPSIDE_OVERWRITE_SPLINES=1
export UPSIDE_NPT_ENABLE=1
export UPSIDE_BAROSTAT_TYPE=0  # Berendsen for equilibration stages

python3 prepare_martini.py "$PDB_ID"

PREPARED_FILE="${RUN_DIR}/test.input.up"
INPUT_FILE="${INPUTS_DIR}/${PDB_ID}.up"

if [ ! -f "$PREPARED_FILE" ]; then
    echo "ERROR: Input preparation failed - file not found: $PREPARED_FILE"
    exit 1
fi

cp -f "$PREPARED_FILE" "$INPUT_FILE"
echo "Input prepared: $INPUT_FILE"
echo

# =============================================================================
# STAGE 6.0: MINIMIZATION WITH SOFT-CORE
# =============================================================================
echo "=== Stage 6.0: Minimization with Soft-Core Potentials ==="

STAGE_60_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
cp -f "$INPUT_FILE" "$STAGE_60_FILE"

# Enable soft-core potentials
export UPSIDE_SOFTEN_LJ=1
export UPSIDE_LJ_ALPHA=0.2
export UPSIDE_SOFTEN_COULOMB=1
export UPSIDE_SLATER_ALPHA=2.0

CMD_60=(
    "$UPSIDE_EXECUTABLE"
    "$STAGE_60_FILE"
    "--duration" "0"
    "--frame-interval" "1"
    "--temperature" "$TEMPERATURE"
    "--time-step" "0.1"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "$MIN_STEPS_60"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)

echo "Running minimization with soft-core..."
if "${CMD_60[@]}" 2>&1 | tee "${LOG_DIR}/stage_6.0.log"; then
    echo "Stage 6.0 completed"
else
    echo "ERROR: Stage 6.0 failed!"
    exit 1
fi
echo

# =============================================================================
# STAGE 6.1: MINIMIZATION WITHOUT SOFT-CORE
# =============================================================================
echo "=== Stage 6.1: Minimization without Soft-Core ==="

STAGE_61_FILE="${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
cp -f "$STAGE_60_FILE" "$STAGE_61_FILE"
python3 set_initial_position.py "$STAGE_60_FILE" "$STAGE_61_FILE"

# Disable soft-core potentials
export UPSIDE_SOFTEN_LJ=0
export UPSIDE_SOFTEN_COULOMB=0

CMD_61=(
    "$UPSIDE_EXECUTABLE"
    "$STAGE_61_FILE"
    "--duration" "0"
    "--frame-interval" "1"
    "--temperature" "$TEMPERATURE"
    "--time-step" "0.1"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "$MIN_STEPS_61"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)

echo "Running minimization without soft-core..."
if "${CMD_61[@]}" 2>&1 | tee "${LOG_DIR}/stage_6.1.log"; then
    echo "Stage 6.1 completed"
else
    echo "ERROR: Stage 6.1 failed!"
    exit 1
fi
echo

# =============================================================================
# HELPER FUNCTION: RUN MD STAGE
# =============================================================================
run_md_stage() {
    local stage_name=$1
    local prev_file=$2
    local timestep=$3
    local steps=$4
    local barostat_type=$5

    local stage_file="${CHECKPOINT_DIR}/${PDB_ID}.stage_${stage_name}.up"

    echo "=== Stage ${stage_name}: MD Equilibration ==="
    echo "  Timestep: ${timestep}"
    echo "  Steps: ${steps}"
    echo "  Barostat: $([ $barostat_type -eq 1 ] && echo 'Parrinello-Rahman' || echo 'Berendsen')"

    # Copy and update initial positions
    cp -f "$prev_file" "$stage_file"
    python3 set_initial_position.py "$prev_file" "$stage_file"

    # Update barostat type in HDF5 file
    python3 - <<END
import h5py
with h5py.File("$stage_file", 'r+') as f:
    if '/input/barostat' in f:
        f['/input/barostat'].attrs['type'] = $barostat_type
END

    local cmd=(
        "$UPSIDE_EXECUTABLE"
        "$stage_file"
        "--duration" "$steps"
        "--frame-interval" "$FRAME_INTERVAL"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$timestep"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "vel_verlet"
        "--disable-recentering"
    )

    if "${cmd[@]}" 2>&1 | tee "${LOG_DIR}/stage_${stage_name}.log"; then
        echo "Stage ${stage_name} completed"
    else
        echo "ERROR: Stage ${stage_name} failed!"
        exit 1
    fi
    echo

    echo "$stage_file"
}

# =============================================================================
# STAGES 6.2-6.6: MD EQUILIBRATION WITH BERENDSEN
# =============================================================================
STAGE_62_FILE=$(run_md_stage "6.2" "$STAGE_61_FILE" "0.002" "$MD_STEPS_62" "0")
STAGE_63_FILE=$(run_md_stage "6.3" "$STAGE_62_FILE" "0.005" "$MD_STEPS_63" "0")
STAGE_64_FILE=$(run_md_stage "6.4" "$STAGE_63_FILE" "0.010" "$MD_STEPS_64" "0")
STAGE_65_FILE=$(run_md_stage "6.5" "$STAGE_64_FILE" "0.015" "$MD_STEPS_65" "0")
STAGE_66_FILE=$(run_md_stage "6.6" "$STAGE_65_FILE" "0.020" "$MD_STEPS_66" "0")

# =============================================================================
# STAGE 7.0: PRODUCTION WITH PARRINELLO-RAHMAN
# =============================================================================
STAGE_70_FILE=$(run_md_stage "7.0" "$STAGE_66_FILE" "0.020" "$MD_STEPS_70" "1")

# =============================================================================
# SUMMARY
# =============================================================================
echo "=== 8-Stage Workflow Complete ==="
echo "Checkpoints:"
echo "  Stage 6.0: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.0.up"
echo "  Stage 6.1: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.1.up"
echo "  Stage 6.2: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.2.up"
echo "  Stage 6.3: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.3.up"
echo "  Stage 6.4: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.4.up"
echo "  Stage 6.5: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.5.up"
echo "  Stage 6.6: ${CHECKPOINT_DIR}/${PDB_ID}.stage_6.6.up"
echo "  Stage 7.0: ${CHECKPOINT_DIR}/${PDB_ID}.stage_7.0.up (Production)"
echo
echo "Done."
