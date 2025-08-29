#!/bin/bash

# MARTINI 3.0 Optimized Simulation Script
# Complete workflow: prepare input, optimize interaction table, run simulation, generate VTF

set -e  # Exit on any error

# Configuration
PDB_ID="1rkl"
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
INPUT_FILE="${INPUTS_DIR}/${PDB_ID}.up"   # normalized input path
OUTPUT_FILE="${RUN_DIR}/test.run.optimized.up"
LOG_FILE="${RUN_DIR}/test.run.optimized.log"
VTF_FILE="${RUN_DIR}/${PDB_ID}.vtf"   # VTF in run directory

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
DURATION=200
FRAME_INTERVAL=20
TEMPERATURE=0.8
TIME_STEP=0.01
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL=1
PRESSURE=0.000020652136
BAROSTAT_TIMESCALE=0.001
BAROSTAT_INTERVAL=1.0
SEED=12345
INTEGRATOR="npt"

echo "Simulation parameters:"
echo "  Duration: $DURATION steps"
echo "  Frame interval: $FRAME_INTERVAL steps"
echo "  Temperature: $TEMPERATURE UPSIDE units (~280 K)"
echo "  Time step: $TIME_STEP"
echo "  Thermostat timescale: $THERMOSTAT_TIMESCALE"
echo "  Pressure: $PRESSURE UPSIDE units (1 bar = 0.000020652136 Å⁻³)"
echo "  Integrator: $INTEGRATOR"
echo

# Step 1: Prepare input files (produce example/16.MARTINI/outputs/martini_test/test.input.up)
echo "=== Step 1: Preparing Input Files ==="
echo "Running prepare_martini.py..."
source ../../source.sh

python prepare_martini.py

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
python optimize_interaction_table.py "$INPUT_FILE" "$TMP_OPT_FILE"

if [ ! -f "$TMP_OPT_FILE" ]; then
    echo "ERROR: Optimization failed - optimized file not created!"
    exit 1
fi

mv -f "$TMP_OPT_FILE" "$INPUT_FILE"
OPTIMIZED_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
echo "Optimization complete!"
echo "  Optimized input: $INPUT_FILE ($OPTIMIZED_SIZE)"
echo

# Step 3: Run simulation with normalized/optimized input
echo "=== Step 3: Running Simulation ==="
CMD=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    "--duration" "$DURATION"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--pressure" "$PRESSURE"
    "--barostat-timescale" "$BAROSTAT_TIMESCALE"
    "--barostat-interval" "$BAROSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "$INTEGRATOR"
    "--mass" "72.0"
    "--disable-initial-thermalization"
    "--restart-using-momentum"
    "--output" "$OUTPUT_FILE"
)

echo "Starting simulation..."
echo "Command: ${CMD[*]}"
echo

START_TIME=$(date +%s)
if "${CMD[@]}" 2>&1 | tee "$LOG_FILE"; then
    echo
    echo "=== Simulation Complete! ==="
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    echo "Duration: ${RUNTIME} seconds"
    echo "Results saved to: $OUTPUT_FILE"
    echo "Log saved to: $LOG_FILE"
    if [ -f "$OUTPUT_FILE" ]; then
        OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "Output file size: $OUTPUT_SIZE"
    fi
else
    echo "ERROR: Simulation failed!"
    exit 1
fi

# Step 4: Generate VTF at RUN_DIR/pdb_id.vtf
echo
echo "=== Step 4: Generating VTF ==="
echo "Generating VTF file: $VTF_FILE"
if python extract_martini_vtf.py "$OUTPUT_FILE" "$VTF_FILE" "$INPUT_FILE"; then
    VTF_SIZE=$(du -h "$VTF_FILE" | cut -f1)
    echo "VTF file generated successfully: $VTF_FILE ($VTF_SIZE)"
else
    echo "ERROR: VTF file generation failed!"
    exit 1
fi

echo
echo "=== Complete Workflow Summary ==="
echo "Input (normalized): $INPUT_FILE ($OPTIMIZED_SIZE)"
echo "Simulation output: $OUTPUT_FILE${OUTPUT_SIZE:+ ($OUTPUT_SIZE)}"
echo "VTF: $VTF_FILE ($VTF_SIZE)"
echo
echo "Done."


