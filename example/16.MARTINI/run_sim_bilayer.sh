#!/bin/bash

# MARTINI 3.0 Complete Multi-Stage Simulation Workflow
# Stages: Prepare -> Minimization -> NPT Equilibration -> NVT Production -> VTF Generation
# Each stage uses the output checkpoint from the previous stage as its starting point.
#
# ENVIRONMENT VARIABLES:
#   PDB_ID                - PDB identifier (default: bilayer)
#   UPSIDE_NPT_ENABLE     - Enable NPT barostat (1=NPT for equilibration)
#   UPSIDE_NPT_TARGET_PXY - Target lateral pressure (reduced units)
#   UPSIDE_NPT_TARGET_PZ  - Target normal pressure (reduced units)
#   UPSIDE_NPT_TAU        - Barostat time constant
#   UPSIDE_NPT_INTERVAL   - Steps between barostat applications

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PDB_ID="${PDB_ID:-bilayer}"

# Directories
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"

# Filenames
INPUT_FILE="${INPUTS_DIR}/${PDB_ID}.up"
MINIMIZED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.minimized.up"
NPT_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt.up"
NVT_FILE="${CHECKPOINT_DIR}/${PDB_ID}.nvt.up"
VTF_FILE="${RUN_DIR}/${PDB_ID}.vtf"
LOG_DIR="${RUN_DIR}/logs"

# Simulation parameters
TEMPERATURE=0.8
TIME_STEP=0.1
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED=12345

# Stage durations (in MD steps)
MIN_STEPS="${MIN_STEPS:-500}"
NPT_STEPS="${NPT_STEPS:-2000}"
NVT_STEPS="${NVT_STEPS:-5000}"
FRAME_INTERVAL="${FRAME_INTERVAL:-20}"

# Softening parameters for minimization
export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-1}
export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.2}
export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-1}
export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-2.0}

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

# Create directories
mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

echo "=== MARTINI 3.0 Multi-Stage Workflow ==="
echo "PDB ID: $PDB_ID"
echo "Stages: Prepare -> Minimization -> NPT -> NVT -> VTF"
echo "  MIN_STEPS:  $MIN_STEPS"
echo "  NPT_STEPS:  $NPT_STEPS"
echo "  NVT_STEPS:  $NVT_STEPS"
echo

# =============================================================================
# STAGE 1: PREPARE INPUT FILES
# =============================================================================
echo "=== Stage 1: Preparing Input Files ==="
source ../../source.sh
source ../../.venv/bin/activate

# Prepare with softened potentials enabled (for minimization compatibility)
export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}
python3 prepare_martini.py "$PDB_ID"

PREPARED_FILE="${RUN_DIR}/test.input.up"
if [ ! -f "$PREPARED_FILE" ]; then
    echo "ERROR: Input preparation failed - file not found: $PREPARED_FILE"
    exit 1
fi

# Copy to normalized location
cp -f "$PREPARED_FILE" "$INPUT_FILE"
echo "Input prepared: $INPUT_FILE ($(du -h "$INPUT_FILE" | cut -f1))"
echo

# =============================================================================
# STAGE 2: ENERGY MINIMIZATION
# =============================================================================
echo "=== Stage 2: Energy Minimization ==="
echo "Running $MIN_STEPS steps with softened potentials"
echo "Input:  $INPUT_FILE"
echo "Output: $MINIMIZED_FILE"

CMD_MIN=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$MIN_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "1000"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)

START_TIME=$(date +%s)
if "${CMD_MIN[@]}" 2>&1 | tee "${LOG_DIR}/minimization.log"; then
    END_TIME=$(date +%s)
    echo "Minimization completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: Minimization failed!"
    exit 1
fi

# Save minimized checkpoint (input file was modified in-place)
cp -f "$INPUT_FILE" "$MINIMIZED_FILE"
echo "Minimized checkpoint: $MINIMIZED_FILE"
echo

# =============================================================================
# STAGE 3: NPT EQUILIBRATION (pressure coupling enabled)
# =============================================================================
echo "=== Stage 3: NPT Equilibration ==="
echo "Running $NPT_STEPS steps with Berendsen barostat"
echo "Input:  $MINIMIZED_FILE"
echo "Output: $NPT_FILE"

# Enable NPT for this stage
export UPSIDE_NPT_ENABLE=1
export UPSIDE_NPT_TARGET_PXY=${UPSIDE_NPT_TARGET_PXY:-1.0}
export UPSIDE_NPT_TARGET_PZ=${UPSIDE_NPT_TARGET_PZ:-1.0}
export UPSIDE_NPT_TAU=${UPSIDE_NPT_TAU:-1.0}
export UPSIDE_NPT_INTERVAL=${UPSIDE_NPT_INTERVAL:-10}

echo "NPT settings: Pxy=${UPSIDE_NPT_TARGET_PXY} Pz=${UPSIDE_NPT_TARGET_PZ} tau=${UPSIDE_NPT_TAU} interval=${UPSIDE_NPT_INTERVAL}"

# Work from minimized checkpoint
cp -f "$MINIMIZED_FILE" "$NPT_FILE"

CMD_NPT=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_FILE"
    "--duration" "$NPT_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--continue"
)

START_TIME=$(date +%s)
if "${CMD_NPT[@]}" 2>&1 | tee "${LOG_DIR}/npt_equilibration.log"; then
    END_TIME=$(date +%s)
    echo "NPT equilibration completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT equilibration failed!"
    exit 1
fi

echo "NPT checkpoint: $NPT_FILE"
echo

# =============================================================================
# STAGE 4: NVT PRODUCTION (fixed volume)
# =============================================================================
echo "=== Stage 4: NVT Production ==="
echo "Running $NVT_STEPS steps with fixed volume (NPT disabled)"
echo "Input:  $NPT_FILE"
echo "Output: $NVT_FILE"

# Disable NPT for production
export UPSIDE_NPT_ENABLE=0

# Work from NPT checkpoint
cp -f "$NPT_FILE" "$NVT_FILE"

CMD_NVT=(
    "$UPSIDE_EXECUTABLE"
    "$NVT_FILE"
    "--duration" "$NVT_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--continue"
)

START_TIME=$(date +%s)
if "${CMD_NVT[@]}" 2>&1 | tee "${LOG_DIR}/nvt_production.log"; then
    END_TIME=$(date +%s)
    echo "NVT production completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NVT production failed!"
    exit 1
fi

echo "NVT production checkpoint: $NVT_FILE"
echo

# =============================================================================
# STAGE 5: VTF GENERATION
# =============================================================================
echo "=== Stage 5: Generating VTF ==="
echo "Extracting trajectory from: $NVT_FILE"
echo "Output VTF: $VTF_FILE"

if python3 extract_martini_vtf.py "$NVT_FILE" "$VTF_FILE" "$NVT_FILE" "$PDB_ID"; then
    VTF_SIZE=$(du -h "$VTF_FILE" | cut -f1)
    echo "VTF file generated: $VTF_FILE ($VTF_SIZE)"
else
    echo "ERROR: VTF generation failed!"
    exit 1
fi

echo
echo "=== Workflow Complete ==="
echo "Checkpoints:"
echo "  Prepared:   $INPUT_FILE"
echo "  Minimized:  $MINIMIZED_FILE"
echo "  NPT:        $NPT_FILE"
echo "  NVT:        $NVT_FILE"
echo "Trajectory:   $VTF_FILE"
echo
echo "To visualize: vmd $VTF_FILE"
echo "Done."
