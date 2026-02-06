#!/bin/bash
source ../../source.sh
source ../../.venv/bin/activate

# MARTINI 3.0 Water Simulation Workflow
# Stages: Prepare -> Minimization -> NPT Equilibration -> NPT Production -> VTF Generation

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PDB_ID="water"

# Directories
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_water_test"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"

# Filenames - per-stage .up files
PREPARED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.prepared.up"
MINIMIZED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.minimized.up"
NPT_EQUIL_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_equil.up"
NPT_PROD_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_prod.up"
VTF_FILE="${RUN_DIR}/${PDB_ID}.vtf"
LOG_DIR="${RUN_DIR}/logs"

# Simulation parameters
TEMPERATURE=0.8
TIME_STEP=0.1
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL="-1"
SEED=12345

# Stage durations (in MD steps)
MIN_STEPS="500"
NPT_EQUIL_STEPS="2000"
NPT_PROD_STEPS="5000"
FRAME_INTERVAL="20"

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

echo "=== MARTINI 3.0 Water Simulation Workflow ===="
echo "PDB ID: $PDB_ID"
echo "Stages: Prepare -> Minimization -> NPT Equilibration -> NPT Production -> VTF"
echo "  NPT Equilibration: $NPT_EQUIL_STEPS steps"
echo "  NPT Production:     $NPT_PROD_STEPS steps"
echo

# =============================================================================
# STAGE 1: PREPARE INPUT FILE
# =============================================================================
echo "=== Stage 1: Preparing Input File (${PREPARED_FILE}) ===="
source ../../source.sh
source ../../.venv/bin/activate

# Prepare with NPT enabled (isotropic pressure coupling for water)
export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}
export UPSIDE_NPT_ENABLE=${UPSIDE_NPT_ENABLE:-1}
export UPSIDE_NPT_SEMI=${UPSIDE_NPT_SEMI:-0}  # Isotropic pressure coupling
export UPSIDE_NPT_COMPRESSIBILITY=${UPSIDE_NPT_COMPRESSIBILITY:-2.1782}  # Water compressibility (4.5e-5 bar⁻¹)
export UPSIDE_SIMULATION_STAGE="minimization"
python3 prepare_martini.py "$PDB_ID" --stage "minimization" "$RUN_DIR"

# Move prepared file to checkpoint directory
PREPARED_TMP="${RUN_DIR}/test.input.up"
if [ -f "$PREPARED_TMP" ]; then
    mv -f "$PREPARED_TMP" "$PREPARED_FILE"
    echo "Input prepared: $PREPARED_FILE ($(du -h "$PREPARED_FILE" | cut -f1))"
else
    echo "ERROR: Input preparation failed - file not found: $PREPARED_TMP"
    exit 1
fi

echo

# =============================================================================
# STAGE 2: ENERGY MINIMIZATION (${MINIMIZED_FILE})
# =============================================================================
echo "=== Stage 2: Energy Minimization (${MINIMIZED_FILE}) ===="
echo "Running gradient descent minimization"
echo "Input:  $PREPARED_FILE"
echo "Output: $MINIMIZED_FILE"

# Copy prepared file to minimization file
cp -f "$PREPARED_FILE" "$MINIMIZED_FILE"

CMD_MIN=(
    "$UPSIDE_EXECUTABLE"
    "$MINIMIZED_FILE"
    "--duration" "0"
    "--frame-interval" "1"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "$MIN_STEPS"
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

echo "Minimized checkpoint: $MINIMIZED_FILE"
echo

# =============================================================================
# STAGE 3: NPT EQUILIBRATION (Berendsen barostat)
# =============================================================================
echo "=== Stage 3: NPT Equilibration (${NPT_EQUIL_FILE}) ===="
echo "Running $NPT_EQUIL_STEPS steps"
echo "Input:  $MINIMIZED_FILE"
echo "Output: $NPT_EQUIL_FILE"

# Generate NPT equilibration file
export UPSIDE_SIMULATION_STAGE="npt_equil"
export UPSIDE_NPT_SEMI=${UPSIDE_NPT_SEMI:-0}  # Isotropic pressure coupling
export UPSIDE_NPT_COMPRESSIBILITY=${UPSIDE_NPT_COMPRESSIBILITY:-2.1782}  # Water compressibility (4.5e-5 bar⁻¹)
python3 prepare_martini.py "$PDB_ID" --stage "npt_equil" "$RUN_DIR"

# Move prepared file to checkpoint directory
NPT_EQUIL_TMP="${RUN_DIR}/test.input.up"
if [ -f "$NPT_EQUIL_TMP" ]; then
    mv -f "$NPT_EQUIL_TMP" "$NPT_EQUIL_FILE"
else
    echo "ERROR: NPT equilibration preparation failed - file not found: $NPT_EQUIL_TMP"
    exit 1
fi

# Set initial position in NPT file to last frame from minimization
python3 set_initial_position.py "$MINIMIZED_FILE" "$NPT_EQUIL_FILE"

CMD_NPT_EQUIL=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_EQUIL_FILE"
    "--duration" "$NPT_EQUIL_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
)

START_TIME=$(date +%s)
if "${CMD_NPT_EQUIL[@]}" 2>&1 | tee "${LOG_DIR}/npt_equilibration.log"; then
    END_TIME=$(date +%s)
    echo "NPT equilibration completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT equilibration failed!"
    exit 1
fi

echo "NPT equilibration checkpoint: $NPT_EQUIL_FILE"
echo

# =============================================================================
# STAGE 4: NPT PRODUCTION (Parrinello-Rahman barostat)
# =============================================================================
echo "=== Stage 4: NPT Production (${NPT_PROD_FILE}) ===="
echo "Running $NPT_PROD_STEPS steps with Parrinello-Rahman barostat"
echo "Input:  $NPT_EQUIL_FILE"
echo "Output: $NPT_PROD_FILE"

# Generate NPT production file
export UPSIDE_SIMULATION_STAGE="npt_prod"
export UPSIDE_NPT_SEMI=${UPSIDE_NPT_SEMI:-0}  # Isotropic pressure coupling
export UPSIDE_NPT_COMPRESSIBILITY=${UPSIDE_NPT_COMPRESSIBILITY:-2.1782}  # Water compressibility (4.5e-5 bar⁻¹)
python3 prepare_martini.py "$PDB_ID" --stage "npt_prod" "$RUN_DIR"

# Move prepared file to checkpoint directory
NPT_PROD_TMP="${RUN_DIR}/test.input.up"
if [ -f "$NPT_PROD_TMP" ]; then
    mv -f "$NPT_PROD_TMP" "$NPT_PROD_FILE"
else
    echo "ERROR: NPT production preparation failed - file not found: $NPT_PROD_TMP"
    exit 1
fi

# Set initial position in production file to last frame from equilibration
python3 set_initial_position.py "$NPT_EQUIL_FILE" "$NPT_PROD_FILE"

CMD_NPT_PROD=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_PROD_FILE"
    "--duration" "$NPT_PROD_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
)

START_TIME=$(date +%s)
if "${CMD_NPT_PROD[@]}" 2>&1 | tee "${LOG_DIR}/npt_production.log"; then
    END_TIME=$(date +%s)
    echo "NPT production completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT production failed!"
    exit 1
fi

echo "NPT production checkpoint: $NPT_PROD_FILE"
echo

# =============================================================================
# STAGE 5: VTF GENERATION (ALL STAGES)
# =============================================================================
echo "=== Stage 5: Generating VTF Files for All Stages ===="

# Minimization VTF
MIN_VTF_FILE="${RUN_DIR}/${PDB_ID}.minimized.vtf"
echo "Extracting minimization trajectory from: $MINIMIZED_FILE"
echo "Output VTF: $MIN_VTF_FILE"
if python3 extract_martini_vtf.py "$MINIMIZED_FILE" "$MIN_VTF_FILE" "$MINIMIZED_FILE" "$PDB_ID"; then
    MIN_VTF_SIZE=$(du -h "$MIN_VTF_FILE" | cut -f1)
    echo "Minimization VTF generated: $MIN_VTF_FILE ($MIN_VTF_SIZE)"
else
    echo "ERROR: Minimization VTF generation failed!"
    exit 1
fi

# NPT Equilibration VTF
NPT_EQUIL_VTF_FILE="${RUN_DIR}/${PDB_ID}.npt_equil.vtf"
echo "Extracting NPT equilibration trajectory from: $NPT_EQUIL_FILE"
echo "Output VTF: $NPT_EQUIL_VTF_FILE"
if python3 extract_martini_vtf.py "$NPT_EQUIL_FILE" "$NPT_EQUIL_VTF_FILE" "$NPT_EQUIL_FILE" "$PDB_ID"; then
    NPT_EQUIL_VTF_SIZE=$(du -h "$NPT_EQUIL_VTF_FILE" | cut -f1)
    echo "NPT equilibration VTF generated: $NPT_EQUIL_VTF_FILE ($NPT_EQUIL_VTF_SIZE)"
else
    echo "ERROR: NPT equilibration VTF generation failed!"
    exit 1
fi

# NPT Production VTF
NPT_PROD_VTF_FILE="${RUN_DIR}/${PDB_ID}.npt_prod.vtf"
echo "Extracting NPT production trajectory from: $NPT_PROD_FILE"
echo "Output VTF: $NPT_PROD_VTF_FILE"
if python3 extract_martini_vtf.py "$NPT_PROD_FILE" "$NPT_PROD_VTF_FILE" "$NPT_PROD_FILE" "$PDB_ID"; then
    NPT_PROD_VTF_SIZE=$(du -h "$NPT_PROD_VTF_FILE" | cut -f1)
    echo "NPT production VTF generated: $NPT_PROD_VTF_FILE ($NPT_PROD_VTF_SIZE)"
else
    echo "ERROR: NPT production VTF generation failed!"
    exit 1
fi

echo
echo "=== Workflow Complete ===="
echo "Checkpoints:"
echo "  Prepared:          $PREPARED_FILE"
echo "  Minimized:         $MINIMIZED_FILE"
echo "  NPT Equilibration: $NPT_EQUIL_FILE"
echo "  NPT Production:    $NPT_PROD_FILE"
echo "Trajectories:"
echo "  Minimization:      $MIN_VTF_FILE"
echo "  NPT Equilibration: $NPT_EQUIL_VTF_FILE"
echo "  NPT Production:    $NPT_PROD_VTF_FILE"
echo
echo "To visualize: vmd <vtf_file>"
echo "Done."
