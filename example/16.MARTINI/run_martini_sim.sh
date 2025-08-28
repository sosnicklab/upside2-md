#!/bin/bash

# MARTINI 3.0 Protein-Lipid System Simulation Script
# Runs the actual UPSIDE simulation using the prepared input file.

set -e  # Exit on any error

# Configuration
RUN_DIR="outputs/martini_test"
INPUT_FILE="${RUN_DIR}/test.input.up"
OUTPUT_FILE="${RUN_DIR}/test.run.up"
LOG_FILE="${RUN_DIR}/test.run.log"

# Check if UPSIDE_HOME is set
if [ -z "$UPSIDE_HOME" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set!"
    echo "Please set it to your UPSIDE installation directory."
    exit 1
fi

UPSIDE_HOME="$UPSIDE_HOME"
UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "Please run prepare_martini.py first."
    exit 1
fi

# Check if UPSIDE executable exists
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    echo "Please build UPSIDE first."
    exit 1
fi

echo "=== MARTINI 3.0 Protein-Lipid System Simulation ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Log file: $LOG_FILE"
echo "UPSIDE executable: $UPSIDE_EXECUTABLE"
echo

# Simulation parameters (from original run_martini.py)
DURATION=200
FRAME_INTERVAL=20
TEMPERATURE=0.8
TIME_STEP=0.01
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL=1
PRESSURE=0.0
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
echo "  Pressure: $PRESSURE"
echo "  Integrator: $INTEGRATOR"
echo

# Build UPSIDE command (correct format from original script)
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
    "--disable-initial-thermalization"
    "--restart-using-momentum"
)

echo "Starting simulation..."
echo "Command: ${CMD[*]}"
echo

# Run simulation
START_TIME=$(date +%s)

# Run the command and capture both stdout and stderr
if "${CMD[@]}" 2>&1 | tee "$LOG_FILE"; then
    echo
    echo "=== Simulation Complete! ==="
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "Duration: ${DURATION} seconds"
    echo "Output file: $OUTPUT_FILE"
    echo "Log file: $LOG_FILE"
    
    # Check if output file was created
    if [ -f "$OUTPUT_FILE" ]; then
        FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "Output file size: $FILE_SIZE"
    else
        echo "WARNING: Output file was not created!"
    fi
else
    echo
    echo "ERROR: Simulation failed!"
    exit 1
fi
