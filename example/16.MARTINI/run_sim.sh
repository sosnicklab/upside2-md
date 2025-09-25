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

# Simulation parameters (from original run_martini.py)
DURATION=5000
FRAME_INTERVAL=20
TEMPERATURE=0.8
TIME_STEP=0.01
THERMOSTAT_TIMESCALE=0.135
#THERMOSTAT_TIMESCALE=5
#THERMOSTAT_INTERVAL=-1
SEED=12345
# Integrators: allow separate choices for softened/min/min and regular stages
INTEGRATOR_SOFT="${INTEGRATOR_SOFT:-nvtc}"
INTEGRATOR_MIN="${INTEGRATOR_MIN:-nvtc}"
INTEGRATOR_REG="${INTEGRATOR_REG:-nvtc}"
MAX_FORCE="${MAX_FORCE:-50}"

echo "Simulation parameters:"
echo "  Duration: $DURATION steps"
echo "  Frame interval: $FRAME_INTERVAL steps"
echo "  Temperature: $TEMPERATURE UPSIDE units (~280 K)"
echo "  Time step: $TIME_STEP"
echo "  Thermostat timescale: $THERMOSTAT_TIMESCALE"
echo "  Integrators: soft=$INTEGRATOR_SOFT, min=$INTEGRATOR_MIN, reg=$INTEGRATOR_REG"
echo

# Step 1: Prepare input files (produce example/16.MARTINI/outputs/martini_test/test.input.up)
echo "=== Step 1: Preparing Input Files ==="
echo "Running prepare_martini.py with PDB ID: $PDB_ID"

# Set softening options based on potential mode
if [ "$POTENTIAL_MODE" = "soft" ]; then
    echo "Using softened potentials"
    export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-1}
    export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.01}
    export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-1}
    export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-1.0}
else
    echo "Using regular potentials"
    export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-0}
    export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.0}
    export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-0}
    export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-0.0}
fi

export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}

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

# Step 3: Simulation run
if [ "$POTENTIAL_MODE" = "regular" ]; then
    echo "=== Step 3: Running Regular Potential Simulation ==="
    echo "Running for $DURATION steps with regular potentials"
    
    # Single stage: regular potential for full duration
    CMD_REG=(
        "$UPSIDE_EXECUTABLE"
        "$INPUT_FILE"
        "--duration" "$DURATION"
        "--frame-interval" "$FRAME_INTERVAL"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$TIME_STEP"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "$INTEGRATOR"
        "--max-force" "$MAX_FORCE"
        "--disable-recentering"
    )
    echo "Command (regular): ${CMD_REG[*]}"
    START_TIME=$(date +%s)
    if "${CMD_REG[@]}" 2>&1 | tee "$LOG_FILE"; then
        END_TIME=$(date +%s)
        echo "Regular simulation completed in $((END_TIME - START_TIME)) seconds"
    else
        echo "ERROR: Regular simulation failed!"
        exit 1
    fi
    
else
    # Soft potential mode (original two-stage logic)
    if [ "$SOFT_RUN_MODE" = "soft_only" ]; then
        echo "=== Step 3: Running Softened-Only Simulation ==="
    else
        echo "=== Step 3: Running Two-Stage Simulation (softened -> regular) ==="
    fi

    SOFT_STEPS=100
    REMAINING_STEPS=$(( DURATION > SOFT_STEPS ? DURATION - SOFT_STEPS : 0 ))

# Stage 3.1: softened run (writes to INPUT_FILE:/output)
echo "-- Stage 1: Softened potential for $SOFT_STEPS steps --"
CMD_SOFT=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$SOFT_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR_SOFT"
    "--max-force" "$MAX_FORCE"
    "--disable-recentering"
)
echo "Command (soft): ${CMD_SOFT[*]}"
START_TIME_SOFT=$(date +%s)
if "${CMD_SOFT[@]}" 2>&1 | tee "$LOG_FILE"; then
    END_TIME_SOFT=$(date +%s)
    echo "Softened run completed in $((END_TIME_SOFT - START_TIME_SOFT)) seconds"
else
    echo "ERROR: Softened stage failed!"
    exit 1
fi

# Stage 3.2: If requested, prepare restart from softened end state and disable softening in INPUT_FILE
if [ "$SOFT_RUN_MODE" = "soft_then_regular" ] && [ "$REMAINING_STEPS" -gt 0 ]; then
python3 - "$INPUT_FILE" << 'PYEOF'
import sys, h5py
path = sys.argv[1]
with h5py.File(path, 'r+') as f:
    # Move /output -> /output_soft (keep softened logs)
    if 'output_soft' in f:
        del f['output_soft']
    if 'output' in f:
        f.move('output', 'output_soft')
    else:
        raise SystemExit('No /output found after softened run')
    # Copy last frame positions to /input/pos (shape: [n_atom, 3, 1])
    pos = f['/output_soft/pos']
    last = pos.shape[0]-1
    last_pos = pos[last]
    # Handle possible systems axis: last_pos could be [1, n_atom, 3]
    if last_pos.ndim == 3 and last_pos.shape[0] == 1:
        last_pos = last_pos[0]
    ipos = f['/input/pos']
    if ipos.shape[0] != last_pos.shape[0] or ipos.shape[1] != last_pos.shape[1]:
        raise SystemExit('Shape mismatch between output pos and input pos')
    ipos[...] = last_pos[:, :, None]
    # Disable softening attributes for regular run
    pot = f['/input/potential/martini_potential']
    if 'coulomb_soften' in pot.attrs:
        pot.attrs.modify('coulomb_soften', 0)
    else:
        pot.attrs['coulomb_soften'] = 0
    if 'lj_soften' in pot.attrs:
        pot.attrs.modify('lj_soften', 0)
    else:
        pot.attrs['lj_soften'] = 0
PYEOF

# Stage 3.2.5: Energy minimization between soft and regular
    echo "-- Stage 2: Minimization (pre-regular) --"
    CMD_MIN=(
        "$UPSIDE_EXECUTABLE"
        "$INPUT_FILE"
        "--duration" "0"                # run only minimization then exit
        "--frame-interval" "$FRAME_INTERVAL"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$TIME_STEP"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "$INTEGRATOR_MIN"
        "--max-force" "${MAX_FORCE_MIN:-20}"
        "--minimize"                  # enable minimization
        "--min-max-iter" "1000"
        "--min-energy-tol" "1e-6"
        "--min-force-tol" "1e-3"
        "--min-step" "0.1"
        "--disable-recentering"
    )
    echo "Command (min): ${CMD_MIN[*]}"
    if ! "${CMD_MIN[@]}" 2>&1 | tee -a "$LOG_FILE"; then
        echo "ERROR: Minimization stage failed!"
        exit 1
    fi

# Stage 3.3: regular run for remaining steps (writes to INPUT_FILE:/output)
    echo "-- Stage 3: Regular potential for $REMAINING_STEPS steps --"
    CMD_REG=(
        "$UPSIDE_EXECUTABLE"
        "$INPUT_FILE"
        "--duration" "$REMAINING_STEPS"
        "--frame-interval" "$FRAME_INTERVAL"
        "--temperature" "$TEMPERATURE"
        "--time-step" "$TIME_STEP"
        "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
        "--thermostat-interval" "$THERMOSTAT_INTERVAL"
        "--seed" "$SEED"
        "--integrator" "$INTEGRATOR_REG"
        "--max-force" "$MAX_FORCE"
        "--disable-recentering"
    )
    echo "Command (regular): ${CMD_REG[*]}"
    START_TIME_REG=$(date +%s)
    if "${CMD_REG[@]}" 2>&1 | tee -a "$LOG_FILE"; then
        END_TIME_REG=$(date +%s)
        echo "Regular run completed in $((END_TIME_REG - START_TIME_REG)) seconds"
    else
        echo "ERROR: Regular stage failed!"
        exit 1
    fi
elif [ "$SOFT_RUN_MODE" = "soft_only" ]; then
    echo "Skipping regular stage: SOFT_RUN_MODE=soft_only"
else
    echo "Skipping regular stage: REMAINING_STEPS=0"
fi

fi  # End of POTENTIAL_MODE check

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
echo "=== Complete Workflow Summary ==="
echo "Unified input/output file: $INPUT_FILE ($(du -h "$INPUT_FILE" | cut -f1))"
echo "  Softened stage saved under: /output_soft in $INPUT_FILE"
echo "  Regular stage saved under: /output in $INPUT_FILE"
echo "VTF: $VTF_FILE ($VTF_SIZE)"
echo
echo "Done."


