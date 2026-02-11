#!/bin/bash
source ../../source.sh
source ../../.venv/bin/activate

# MARTINI 3.0 Bilayer Simulation Workflow (CHARMM-GUI Protocol)
# Stages: Prepare -> Minimization -> NPT Equilibration (Berendsen) -> NPT Production (Parrinello-Rahman) -> VTF Generation
# Each stage uses a separate .up file generated at the start of the stage
#
# ENVIRONMENT VARIABLES:
#   PDB_ID                - PDB identifier (default: bilayer)
#   UPSIDE_NPT_TARGET_PXY - Target lateral pressure (default: 0.000020659 = 1 bar)
#   UPSIDE_NPT_TARGET_PZ  - Target normal pressure (default: 0.000020659 = 1 bar)
#   UPSIDE_NPT_TAU        - Barostat time constant (default: 1.0)
#   UPSIDE_NPT_INTERVAL   - Steps between barostat applications (default: 10)
#   UPSIDE_BAROSTAT_TYPE  - Barostat type: 0=Berendsen, 1=Parrinello-Rahman

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION
# =============================================================================
# Extract PDB_ID from command line arguments if provided
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

# Default to bilayer if not specified
PDB_ID="${PDB_ID:-bilayer}"

# Directories
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"

# Filenames - per-stage .up files
PREPARED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.prepared.up"
MINIMIZED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.minimized.up"
NPT_EQUIL_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_equil.up"
NPT_EQUIL_REDUCED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_equil_reduced.up"
NPT_PROD_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_prod.up"
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
NVT_EQUIL_STEPS="${NVT_EQUIL_STEPS:-1000}"
NPT_EQUIL_STEPS="${NPT_EQUIL_STEPS:-2000}"
NPT_PROD_STEPS="${NPT_PROD_STEPS:-5000}"
FRAME_INTERVAL="${FRAME_INTERVAL:-20}"

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

echo "=== MARTINI 3.0 Bilayer Workflow (CHARMM-GUI Protocol) ==="
echo "PDB ID: $PDB_ID"
echo "Stages: Prepare -> Minimization -> NPT Equilibration (Softened, Berendsen) -> NPT Equilibration (Reduced Softening, Berendsen) -> NPT Production (Hard, Parrinello-Rahman) -> VTF"
echo "  NPT Equilibration 1: $NPT_EQUIL_STEPS steps (softened potentials, lj_alpha=0.2, slater_alpha=2.0)"
echo "  NPT Equilibration 2: $NPT_EQUIL_STEPS steps (reduced softening, lj_alpha=0.05, slater_alpha=0.5)"
echo "  NPT Production:     $NPT_PROD_STEPS steps (hard particles)"
echo

# =============================================================================
# STAGE 1: PREPARE INPUT FILE
# =============================================================================
echo "=== Stage 1: Preparing Input File (${PREPARED_FILE}) ==="
source ../../source.sh
source ../../.venv/bin/activate

# Prepare with NPT enabled (so barostat configuration is stored in input file)
export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}
export UPSIDE_NPT_ENABLE=${UPSIDE_NPT_ENABLE:-1}
export UPSIDE_EWALD_ENABLE=${UPSIDE_EWALD_ENABLE:-1}
export UPSIDE_EWALD_ALPHA=${UPSIDE_EWALD_ALPHA:-0.2}
export UPSIDE_EWALD_KMAX=${UPSIDE_EWALD_KMAX:-5}
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
echo "=== Stage 2: Energy Minimization (${MINIMIZED_FILE}) ==="
echo "Running gradient descent minimization with softened potentials"
echo "Input:  $PREPARED_FILE"
echo "Output: $MINIMIZED_FILE"

# Copy prepared file to minimization file
cp -f "$PREPARED_FILE" "$MINIMIZED_FILE"

CMD_MIN=(
    "$UPSIDE_EXECUTABLE"
    "$MINIMIZED_FILE"
    # Required arguments even for minimization (executable parses all args before checking --minimize)
    "--duration" "0"          # 0 duration means no MD steps after minimization
    "--frame-interval" "1"    # Required by executable
    "--temperature" "$TEMPERATURE"  # Required by executable
    "--time-step" "$TIME_STEP"      # Required by executable
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"  # Required by executable
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"    # Required by executable
    "--seed" "$SEED"          # Required by executable
    "--integrator" "vel_verlet"    # Required by executable
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

echo "Minimized checkpoint: $MINIMIZED_FILE"
echo

# =============================================================================
# STAGE 3: NPT EQUILIBRATION WITH SOFTENED POTENTIALS (Berendsen barostat)
# =============================================================================
echo "=== Stage 3: NPT Equilibration (Softened Potentials, Berendsen) (${NPT_EQUIL_FILE}) ==="
echo "Running $NPT_EQUIL_STEPS steps with softened potentials to relax system"
echo "Input:  $MINIMIZED_FILE"
echo "Output: $NPT_EQUIL_FILE"

# Generate NPT equilibration file with softened potentials
export UPSIDE_SIMULATION_STAGE="npt_equil"
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
# STAGE 4: NPT PRODUCTION (Reduced Softening, Berendsen barostat)
# =============================================================================
echo "=== Stage 4: NPT Equilibration (Reduced Softening, Berendsen) (${NPT_EQUIL_REDUCED_FILE}) ==="
echo "Running $NPT_EQUIL_STEPS steps with reduced softening potentials"
echo "Input:  $NPT_EQUIL_FILE"
echo "Output: $NPT_EQUIL_REDUCED_FILE"

# Generate NPT equilibration file with reduced softening potentials
export UPSIDE_SIMULATION_STAGE="npt_equil_reduced"
python3 prepare_martini.py "$PDB_ID" --stage "npt_equil_reduced" "$RUN_DIR"

# Move prepared file to checkpoint directory
NPT_EQUIL_REDUCED_TMP="${RUN_DIR}/test.input.up"
if [ -f "$NPT_EQUIL_REDUCED_TMP" ]; then
    mv -f "$NPT_EQUIL_REDUCED_TMP" "$NPT_EQUIL_REDUCED_FILE"
else
    echo "ERROR: NPT equilibration (reduced) preparation failed - file not found: $NPT_EQUIL_REDUCED_TMP"
    exit 1
fi

# Set initial position in reduced softening file to last frame from previous equilibration
python3 set_initial_position.py "$NPT_EQUIL_FILE" "$NPT_EQUIL_REDUCED_FILE"

CMD_NPT_EQUIL_REDUCED=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_EQUIL_REDUCED_FILE"
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
if "${CMD_NPT_EQUIL_REDUCED[@]}" 2>&1 | tee "${LOG_DIR}/npt_equilibration_reduced.log"; then
    END_TIME=$(date +%s)
    echo "NPT equilibration (reduced softening) completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT equilibration (reduced softening) failed!"
    exit 1
fi

echo "NPT equilibration (reduced softening) checkpoint: $NPT_EQUIL_REDUCED_FILE"
echo

# =============================================================================
# STAGE 5: NPT PRODUCTION (Hard, Parrinello-Rahman barostat)
# =============================================================================
echo "=== Stage 5: NPT Production (Hard, Parrinello-Rahman) (${NPT_PROD_FILE}) ==="
echo "Running $NPT_PROD_STEPS steps with Parrinello-Rahman barostat"
echo "Input:  $NPT_EQUIL_REDUCED_FILE"
echo "Output: $NPT_PROD_FILE"

# Generate NPT production file with hard potentials
export UPSIDE_SIMULATION_STAGE="npt_prod"
python3 prepare_martini.py "$PDB_ID" --stage "npt_prod" "$RUN_DIR"

# Move prepared file to checkpoint directory
NPT_PROD_TMP="${RUN_DIR}/test.input.up"
if [ -f "$NPT_PROD_TMP" ]; then
    mv -f "$NPT_PROD_TMP" "$NPT_PROD_FILE"
else
    echo "ERROR: NPT production preparation failed - file not found: $NPT_PROD_TMP"
    exit 1
fi

# Set initial position in production file to last frame from reduced softening equilibration
python3 set_initial_position.py "$NPT_EQUIL_REDUCED_FILE" "$NPT_PROD_FILE"

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
# STAGE 6: VTF GENERATION (ALL STAGES)
# =============================================================================
echo "=== Stage 6: Generating VTF Files for All Stages ==="

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
echo "=== Workflow Complete ==="
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
