#!/bin/bash

# MARTINI 3.0 Optimized Simulation Script
# Complete workflow: prepare input, optimize interaction table, run simulation, generate VTF
#
# POTENTIAL MODE TOGGLE:
# - Set POTENTIAL_MODE="soft" (default) for softened potentials (two-stage simulation)
# - Set POTENTIAL_MODE="regular" for regular potentials only (single-stage simulation)
# - Can be overridden via environment: POTENTIAL_MODE=regular ./run_sim.sh
#
# PME CONFIGURATION:
# - Set USE_PME="1" (default) to enable Particle Mesh Ewald for long-range Coulomb interactions
# - Set USE_PME="0" to use standard short-range Coulomb cutoff
# - Can be overridden via environment: USE_PME=0 ./run_sim.sh
# - PME provides O(N log N) scaling for large systems vs O(N²) for standard cutoff

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION - Set your PDB ID and potential mode here
# =============================================================================
# Change this value to use a different PDB file (e.g., "bilayer", "protein", etc.)
# The script will look for pdb/{PDB_ID}.MARTINI.pdb
PDB_ID="bilayer"

# POTENTIAL MODE TOGGLE
# Set to "soft" to use softened potentials (default)
# Set to "regular" to use regular potentials only
# Can be overridden via environment variable: POTENTIAL_MODE=regular ./run_sim.sh
POTENTIAL_MODE="${POTENTIAL_MODE:-soft}"

# PME CONFIGURATION
# Set to "1" to enable Particle Mesh Ewald for long-range Coulomb interactions
# Set to "0" to use standard short-range Coulomb cutoff
# Can be overridden via environment variable: USE_PME=1 ./run_sim.sh
# PME provides O(N log N) scaling for large systems vs O(N²) for standard cutoff
# Uses improved 4th-order B-splines and optimized FFT for better accuracy
USE_PME="${USE_PME:-1}"

# THERMOSTAT INTERVAL (optional; set default if unset)
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"

# NVT ENSEMBLE CONFIGURATION
# Disable NPT barostat for NVT ensemble (constant volume)
export UPSIDE_NPT_ENABLE=0

# ADVANCED PME TUNING (optional overrides)
# These can be set via environment variables for fine-tuning:
# UPSIDE_PME_ALPHA=0.2      # Ewald screening parameter (0.1-0.5 range)
# UPSIDE_PME_RCUT=10.0      # Real space cutoff in Angstroms
# UPSIDE_PME_NX=64          # Grid size X (must be power of 2: 16,32,64,128...)
# UPSIDE_PME_NY=64          # Grid size Y (must be power of 2: 16,32,64,128...)
# UPSIDE_PME_NZ=64          # Grid size Z (must be power of 2: 16,32,64,128...)
# UPSIDE_PME_ORDER=4        # B-spline order (3 or 4 recommended)
# =============================================================================

# Validate potential mode
if [ "$POTENTIAL_MODE" != "soft" ] && [ "$POTENTIAL_MODE" != "regular" ]; then
    echo "ERROR: Invalid POTENTIAL_MODE: $POTENTIAL_MODE"
    echo "Valid options: 'soft' or 'regular'"
    exit 1
fi

# Validate PME setting
if [ "$USE_PME" != "0" ] && [ "$USE_PME" != "1" ]; then
    echo "ERROR: Invalid USE_PME: $USE_PME"
    echo "Valid options: '0' (disable) or '1' (enable)"
    exit 1
fi

# Validate PME grid sizes (must be powers of 2 for optimal FFT performance)
if [ "$USE_PME" = "1" ]; then
    validate_power_of_2() {
        local val=$1
        local name=$2
        if [ "$val" -gt 0 ] && [ $((val & (val - 1))) -eq 0 ]; then
            return 0  # Valid power of 2
        else
            echo "WARNING: $name=$val is not a power of 2. FFT performance may be suboptimal."
            echo "Recommended values: 16, 32, 64, 128, 256..."
            return 1
        fi
    }
    
    # Check grid sizes if they're set
    [ -n "$UPSIDE_PME_NX" ] && validate_power_of_2 "$UPSIDE_PME_NX" "UPSIDE_PME_NX"
    [ -n "$UPSIDE_PME_NY" ] && validate_power_of_2 "$UPSIDE_PME_NY" "UPSIDE_PME_NY"
    [ -n "$UPSIDE_PME_NZ" ] && validate_power_of_2 "$UPSIDE_PME_NZ" "UPSIDE_PME_NZ"
fi


# Command line argument overrides the configuration above
if [ $# -eq 1 ]; then
    PDB_ID="$1"
    echo "Command line argument overrides configuration: using PDB ID: $PDB_ID"
elif [ $# -gt 1 ]; then
    echo "Usage: $0 [PDB_ID]"
    echo "  PDB_ID: The PDB identifier (optional, overrides configuration)"
    echo "  Configuration PDB ID: $PDB_ID"
    echo "  Potential mode: $POTENTIAL_MODE"
    exit 1
else
    echo "Using configured PDB ID: $PDB_ID"
fi

echo "Potential mode: $POTENTIAL_MODE"
echo "PME enabled: $USE_PME"
if [ "$USE_PME" = "1" ]; then
    echo "PME configuration:"
    echo "  Grid size: ${UPSIDE_PME_NX:-64}x${UPSIDE_PME_NY:-64}x${UPSIDE_PME_NZ:-64}"
    echo "  Alpha: ${UPSIDE_PME_ALPHA:-0.2}"
    echo "  Real space cutoff: ${UPSIDE_PME_RCUT:-10.0} Å"
    echo "  B-spline order: ${UPSIDE_PME_ORDER:-4}"
    echo "  Features: Improved FFT, 4th-order B-splines, spline tables"
fi
echo "Thermostat interval: $THERMOSTAT_INTERVAL"
echo "Ensemble: NVT (constant volume, no barostat)"

# Configuration
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
INPUT_FILE="${INPUTS_DIR}/${PDB_ID}.up"   # normalized input path
OUTPUT_FILE="${RUN_DIR}/test.run.optimized.up"
LOG_FILE="${RUN_DIR}/test.run.optimized.log"
VTF_FILE="${RUN_DIR}/${PDB_ID}.vtf"   # VTF in run directory

# Softening run mode: soft_only or soft_then_regular (can override via env SOFT_RUN_MODE)
SOFT_RUN_MODE="${SOFT_RUN_MODE:-soft_then_regular}"

# Check if UPSIDE_HOME is set
if [ -z "$UPSIDE_HOME" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set!"
    echo "Please set it to your UPSIDE installation directory."
    exit 1
fi

UPSIDE_HOME="$UPSIDE_HOME"
UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"

# Check if UPSIDE executable exists
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    echo "Please build UPSIDE first."
    exit 1
fi

echo "=== MARTINI 3.0 Complete Optimized Workflow ==="
echo "PDB ID: $PDB_ID"
echo "Normalized input: $INPUT_FILE"
echo "Simulation output: $OUTPUT_FILE"
echo "Log file: $LOG_FILE"
echo "VTF file: $VTF_FILE"
echo "UPSIDE executable: $UPSIDE_EXECUTABLE"
echo

# Ensure directories
mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR"

# GROMACS MARTINI Simulation Stages Parameters
# Following CHARMM-GUI GROMACS MARTINI protocol exactly

# Stage 1: Minimization (steps 6.0 & 6.1)
MIN_STEPS_1=5000    # Step 6.0: 5000 steps (soft potentials)
MIN_STEPS_2=5000    # Step 6.1: 5000 steps (regular potentials)

# Stage 2: Equilibration with restraints (steps 6.2-6.6)
# Step 6.2: Strong restraints (1000 fc, 200 bilayer fc)
EQ_STEPS_1=500     # 500,000 steps / 1000 = 500 steps (for testing)
EQ_DT_1=0.002
EQ_RESTRAINT_1=1000
EQ_BILAYER_RESTRAINT_1=200

# Step 6.3: Medium restraints (500 fc, 100 bilayer fc)
EQ_STEPS_2=200     # 200,000 steps / 1000 = 200 steps (for testing)
EQ_DT_2=0.005
EQ_RESTRAINT_2=500
EQ_BILAYER_RESTRAINT_2=100

# Step 6.4: Weak restraints (250 fc, 50 bilayer fc)
EQ_STEPS_3=100     # 100,000 steps / 1000 = 100 steps (for testing)
EQ_DT_3=0.010
EQ_RESTRAINT_3=250
EQ_BILAYER_RESTRAINT_3=50

# Step 6.5: Very weak restraints (100 fc, 20 bilayer fc)
EQ_STEPS_4=50      # 50,000 steps / 1000 = 50 steps (for testing)
EQ_DT_4=0.015
EQ_RESTRAINT_4=100
EQ_BILAYER_RESTRAINT_4=20

# Step 6.6: Final equilibration (50 fc, 10 bilayer fc)
EQ_STEPS_5=50      # 50,000 steps / 1000 = 50 steps (for testing)
EQ_DT_5=0.020
EQ_RESTRAINT_5=50
EQ_BILAYER_RESTRAINT_5=10

# Stage 3: Production (step 7) - Set to 20 for testing (original: 1,500,000 steps)
PRODUCTION_STEPS=20
PRODUCTION_DT=0.020

# Common parameters
FRAME_INTERVAL=20
TEMPERATURE=0.8  # ~303.15 K in UPSIDE units
THERMOSTAT_TIMESCALE=0.135
SEED=12345

# Integrators: NVT ensemble (constant volume, no barostat)
INTEGRATOR_MIN="${INTEGRATOR_MIN:-nvtc}"
INTEGRATOR_EQ="${INTEGRATOR_EQ:-nvtc}"
INTEGRATOR_PROD="${INTEGRATOR_PROD:-nvtc}"
MAX_FORCE="${MAX_FORCE:-}"

echo "GROMACS MARTINI Simulation Stages:"
echo "  Stage 1: Minimization"
echo "    - Step 6.0: ${MIN_STEPS_1} steps (SOFT potentials)"
echo "    - Step 6.1: ${MIN_STEPS_2} steps (REGULAR potentials)"
echo "  Stage 2: Equilibration with restraints (5 stages)"
echo "    - Step 6.2: Strong restraints (${EQ_RESTRAINT_1} fc, ${EQ_BILAYER_RESTRAINT_1} bilayer fc, ${EQ_STEPS_1} steps, dt=${EQ_DT_1})"
echo "    - Step 6.3: Medium restraints (${EQ_RESTRAINT_2} fc, ${EQ_BILAYER_RESTRAINT_2} bilayer fc, ${EQ_STEPS_2} steps, dt=${EQ_DT_2})"
echo "    - Step 6.4: Weak restraints (${EQ_RESTRAINT_3} fc, ${EQ_BILAYER_RESTRAINT_3} bilayer fc, ${EQ_STEPS_3} steps, dt=${EQ_DT_3})"
echo "    - Step 6.5: Very weak restraints (${EQ_RESTRAINT_4} fc, ${EQ_BILAYER_RESTRAINT_4} bilayer fc, ${EQ_STEPS_4} steps, dt=${EQ_DT_4})"
echo "    - Step 6.6: Final equilibration (${EQ_RESTRAINT_5} fc, ${EQ_BILAYER_RESTRAINT_5} bilayer fc, ${EQ_STEPS_5} steps, dt=${EQ_DT_5})"
echo "  Stage 3: Production (Step 7)"
echo "    - Production: ${PRODUCTION_STEPS} steps (REGULAR potentials, no restraints, dt=${PRODUCTION_DT})"
echo "  Temperature: $TEMPERATURE UPSIDE units (~303.15 K)"
echo "  Thermostat timescale: $THERMOSTAT_TIMESCALE"
echo "  Integrators: min=$INTEGRATOR_MIN, eq=$INTEGRATOR_EQ, prod=$INTEGRATOR_PROD (NVT ensemble)"
echo

# Step 1: Prepare input files (produce example/16.MARTINI/outputs/martini_test/test.input.up)
echo "=== Step 1: Preparing Input Files ==="
echo "Running prepare_martini.py with PDB ID: $PDB_ID"

# Set softening options based on potential mode
if [ "$POTENTIAL_MODE" = "soft" ]; then
    echo "Using softened potentials"
    export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-1}
    export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.2}
    export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-1}
    export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-2.0}
else
    echo "Using regular potentials"
    export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-0}
    export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.0}
    export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-0}
    export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-0.0}
fi
# Ensure force cap default consistent with stabilized run
# export UPSIDE_FORCE_CAP=${UPSIDE_FORCE_CAP:-50}

export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}
# NVT ensemble: disable NPT barostat for constant volume simulation
export UPSIDE_NPT_ENABLE=0

# PME configuration - Optimized settings for improved accuracy and performance
export UPSIDE_USE_PME=${USE_PME}
export UPSIDE_PME_ALPHA=${UPSIDE_PME_ALPHA:-0.2}        # Ewald screening parameter (optimized)
export UPSIDE_PME_RCUT=${UPSIDE_PME_RCUT:-10.0}         # Real space cutoff in Angstroms
export UPSIDE_PME_NX=${UPSIDE_PME_NX:-64}               # Grid size X (power of 2 for optimal FFT)
export UPSIDE_PME_NY=${UPSIDE_PME_NY:-64}               # Grid size Y (power of 2 for optimal FFT)
export UPSIDE_PME_NZ=${UPSIDE_PME_NZ:-64}               # Grid size Z (power of 2 for optimal FFT)
export UPSIDE_PME_ORDER=${UPSIDE_PME_ORDER:-4}          # B-spline interpolation order (4th order for accuracy)

echo "Potential options (env): UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ} UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA} UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB} UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA} UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES}"
echo "PME options (env): UPSIDE_USE_PME=${UPSIDE_USE_PME} UPSIDE_PME_ALPHA=${UPSIDE_PME_ALPHA} UPSIDE_PME_RCUT=${UPSIDE_PME_RCUT} UPSIDE_PME_NX=${UPSIDE_PME_NX} UPSIDE_PME_NY=${UPSIDE_PME_NY} UPSIDE_PME_NZ=${UPSIDE_PME_NZ} UPSIDE_PME_ORDER=${UPSIDE_PME_ORDER}"
source ../../source.sh
source ../../.venv/bin/activate

python3 prepare_martini.py "$PDB_ID"

PREPARED_FILE_PATH="${RUN_DIR}/test.input.up"
if [ ! -f "$PREPARED_FILE_PATH" ]; then
    echo "ERROR: Input preparation failed - file not found: $PREPARED_FILE_PATH"
    exit 1
fi

# Normalize to inputs/pdb_id.up
cp -f "$PREPARED_FILE_PATH" "$INPUT_FILE"
INPUT_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
echo "Input normalized to: $INPUT_FILE ($INPUT_SIZE)"
echo

# Step 2: Optimize the interaction table (overwrite inputs/pdb_id.up)
echo "=== Step 2: Optimizing Interaction Table ==="
TMP_OPT_FILE="${INPUT_FILE}.tmp"
python3 optimize_interaction_table.py "$INPUT_FILE" "$TMP_OPT_FILE"

if [ ! -f "$TMP_OPT_FILE" ]; then
    echo "ERROR: Optimization failed - optimized file not created!"
    exit 1
fi

mv -f "$TMP_OPT_FILE" "$INPUT_FILE"
OPTIMIZED_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
echo "Optimization complete!"
echo "  Optimized input: $INPUT_FILE ($OPTIMIZED_SIZE)"
echo

# Step 3: GROMACS MARTINI Simulation Stages
echo "=== Step 3: Running GROMACS MARTINI Simulation Stages ==="

# Stage 1: Minimization (steps 6.0 & 6.1)
echo "-- Stage 1: Minimization (Steps 6.0 & 6.1) --"

# Step 6.0: First minimization with SOFT potentials (GROMACS free-energy = yes)
echo "Step 6.0: First minimization (${MIN_STEPS_1} steps) with SOFT potentials"
export UPSIDE_SOFTEN_LJ=1
export UPSIDE_LJ_ALPHA=0.2
export UPSIDE_SOFTEN_COULOMB=1
export UPSIDE_SLATER_ALPHA=2.0

CMD_MIN1=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "0"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "0.001"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_MIN"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "$MIN_STEPS_1"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)
[ -n "$MAX_FORCE" ] && CMD_MIN1+=("--max-force" "$MAX_FORCE")
echo "Command (minimization 6.0): ${CMD_MIN1[*]}"
echo "Using SOFT potentials: LJ_ALPHA=${UPSIDE_LJ_ALPHA}, SLATER_ALPHA=${UPSIDE_SLATER_ALPHA}"
START_TIME_MIN1=$(date +%s)
if "${CMD_MIN1[@]}" 2>&1 | tee "$LOG_FILE"; then
    END_TIME_MIN1=$(date +%s)
    echo "Minimization 6.0 completed in $((END_TIME_MIN1 - START_TIME_MIN1)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Minimization 6.0 failed!"
    exit 1
fi

# Step 6.1: Second minimization with REGULAR potentials (GROMACS free-energy = no)
echo "Step 6.1: Second minimization (${MIN_STEPS_2} steps) with REGULAR potentials"
export UPSIDE_SOFTEN_LJ=0
export UPSIDE_LJ_ALPHA=0.0
export UPSIDE_SOFTEN_COULOMB=0
export UPSIDE_SLATER_ALPHA=0.0

CMD_MIN2=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "0"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "0.001"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_MIN"
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "$MIN_STEPS_2"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)
[ -n "$MAX_FORCE" ] && CMD_MIN2+=("--max-force" "$MAX_FORCE")
echo "Command (minimization 6.1): ${CMD_MIN2[*]}"
echo "Using REGULAR potentials: LJ_ALPHA=${UPSIDE_LJ_ALPHA}, SLATER_ALPHA=${UPSIDE_SLATER_ALPHA}"
START_TIME_MIN2=$(date +%s)
if "${CMD_MIN2[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_MIN2=$(date +%s)
    echo "Minimization 6.1 completed in $((END_TIME_MIN2 - START_TIME_MIN2)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Minimization 6.1 failed!"
    exit 1
fi

# Stage 2: Equilibration with restraints (steps 6.2-6.6)
echo "-- Stage 2: Equilibration with Restraints (Steps 6.2-6.6) --"
echo "Switching to regular potentials for equilibration stages"

# Switch to regular potentials for equilibration
export UPSIDE_SOFTEN_LJ=0
export UPSIDE_LJ_ALPHA=0.0
export UPSIDE_SOFTEN_COULOMB=0
export UPSIDE_SLATER_ALPHA=0.0
echo "Using regular potentials: LJ_ALPHA=${UPSIDE_LJ_ALPHA}, SLATER_ALPHA=${UPSIDE_SLATER_ALPHA}"

# Step 6.2: Strong restraints
echo "Step 6.2: Strong restraints (${EQ_RESTRAINT_1} fc, ${EQ_BILAYER_RESTRAINT_1} bilayer fc, ${EQ_STEPS_1} steps, dt=${EQ_DT_1})"
CMD_EQ1=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$EQ_STEPS_1"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$EQ_DT_1"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_EQ"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_EQ1+=("--max-force" "$MAX_FORCE")
echo "Command (equilibration 6.2): ${CMD_EQ1[*]}"
START_TIME_EQ1=$(date +%s)
if "${CMD_EQ1[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_EQ1=$(date +%s)
    echo "Equilibration 6.2 completed in $((END_TIME_EQ1 - START_TIME_EQ1)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Equilibration 6.2 failed!"
    exit 1
fi

# Step 6.3: Medium restraints
echo "Step 6.3: Medium restraints (${EQ_RESTRAINT_2} fc, ${EQ_BILAYER_RESTRAINT_2} bilayer fc, ${EQ_STEPS_2} steps, dt=${EQ_DT_2})"
CMD_EQ2=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$EQ_STEPS_2"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$EQ_DT_2"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_EQ"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_EQ2+=("--max-force" "$MAX_FORCE")
echo "Command (equilibration 6.3): ${CMD_EQ2[*]}"
START_TIME_EQ2=$(date +%s)
if "${CMD_EQ2[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_EQ2=$(date +%s)
    echo "Equilibration 6.3 completed in $((END_TIME_EQ2 - START_TIME_EQ2)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Equilibration 6.3 failed!"
    exit 1
fi

# Step 6.4: Weak restraints
echo "Step 6.4: Weak restraints (${EQ_RESTRAINT_3} fc, ${EQ_BILAYER_RESTRAINT_3} bilayer fc, ${EQ_STEPS_3} steps, dt=${EQ_DT_3})"
CMD_EQ3=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$EQ_STEPS_3"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$EQ_DT_3"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_EQ"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_EQ3+=("--max-force" "$MAX_FORCE")
echo "Command (equilibration 6.4): ${CMD_EQ3[*]}"
START_TIME_EQ3=$(date +%s)
if "${CMD_EQ3[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_EQ3=$(date +%s)
    echo "Equilibration 6.4 completed in $((END_TIME_EQ3 - START_TIME_EQ3)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Equilibration 6.4 failed!"
    exit 1
fi

# Step 6.5: Very weak restraints
echo "Step 6.5: Very weak restraints (${EQ_RESTRAINT_4} fc, ${EQ_BILAYER_RESTRAINT_4} bilayer fc, ${EQ_STEPS_4} steps, dt=${EQ_DT_4})"
CMD_EQ4=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$EQ_STEPS_4"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$EQ_DT_4"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_EQ"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_EQ4+=("--max-force" "$MAX_FORCE")
echo "Command (equilibration 6.5): ${CMD_EQ4[*]}"
START_TIME_EQ4=$(date +%s)
if "${CMD_EQ4[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_EQ4=$(date +%s)
    echo "Equilibration 6.5 completed in $((END_TIME_EQ4 - START_TIME_EQ4)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Equilibration 6.5 failed!"
    exit 1
fi

# Step 6.6: Final equilibration
echo "Step 6.6: Final equilibration (${EQ_RESTRAINT_5} fc, ${EQ_BILAYER_RESTRAINT_5} bilayer fc, ${EQ_STEPS_5} steps, dt=${EQ_DT_5})"
CMD_EQ5=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$EQ_STEPS_5"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$EQ_DT_5"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_EQ"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_EQ5+=("--max-force" "$MAX_FORCE")
echo "Command (equilibration 6.6): ${CMD_EQ5[*]}"
START_TIME_EQ5=$(date +%s)
if "${CMD_EQ5[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_EQ5=$(date +%s)
    echo "Equilibration 6.6 completed in $((END_TIME_EQ5 - START_TIME_EQ5)) seconds"
    # Copy final structure from output to input for next stage
    echo "Preparing structure for next stage..."
    TMP_FILE="${INPUT_FILE}.tmp"
    python3 ../../py/generate_restart_config.py "$INPUT_FILE" "$TMP_FILE"
    mv "$TMP_FILE" "$INPUT_FILE"
else
    echo "ERROR: Equilibration 6.6 failed!"
    exit 1
fi

# Stage 3: Production (step 7)
echo "-- Stage 3: Production (Step 7) --"
echo "Production run (${PRODUCTION_STEPS} steps, dt=${PRODUCTION_DT}) with regular potentials"
echo "Using regular potentials: LJ_ALPHA=${UPSIDE_LJ_ALPHA}, SLATER_ALPHA=${UPSIDE_SLATER_ALPHA}"
CMD_PROD=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$PRODUCTION_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$PRODUCTION_DT"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_PROD"
    "--disable-recentering"
)
[ -n "$MAX_FORCE" ] && CMD_PROD+=("--max-force" "$MAX_FORCE")
echo "Command (production): ${CMD_PROD[*]}"
START_TIME_PROD=$(date +%s)
if "${CMD_PROD[@]}" 2>&1 | tee -a "$LOG_FILE"; then
    END_TIME_PROD=$(date +%s)
    echo "Production completed in $((END_TIME_PROD - START_TIME_PROD)) seconds"
else
    echo "ERROR: Production failed!"
    exit 1
fi

echo
echo "=== Simulation Complete! ==="

# Step 4: Generate VTF at RUN_DIR/pdb_id.vtf from the unified .up file
echo
echo "=== Step 4: Generating VTF ==="
echo "Generating VTF file: $VTF_FILE"
if python3 extract_martini_vtf.py "$INPUT_FILE" "$VTF_FILE" "$INPUT_FILE" "$PDB_ID"; then
    VTF_SIZE=$(du -h "$VTF_FILE" | cut -f1)
    echo "VTF file generated successfully: $VTF_FILE ($VTF_SIZE)"
else
    echo "ERROR: VTF file generation failed!"
    exit 1
fi

echo
echo "=== Complete GROMACS MARTINI Workflow Summary ==="
echo "Unified input/output file: $INPUT_FILE ($(du -h "$INPUT_FILE" | cut -f1))"
echo "  GROMACS MARTINI simulation stages completed:"
echo "    Step 6.0: Minimization with SOFT potentials (free-energy = yes)"
echo "    Step 6.1: Minimization with REGULAR potentials (free-energy = no)"
echo "    Steps 6.2-6.6: Equilibration with REGULAR potentials + restraints"
echo "    Step 7: Production with REGULAR potentials (no restraints)"
echo "  Results saved under: /output in $INPUT_FILE"
echo "VTF: $VTF_FILE ($VTF_SIZE)"
echo
echo "Done."


