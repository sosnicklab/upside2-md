#!/bin/bash

# Define current directory variables
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

# Path math
PROJECT_ROOT="$WORK_DIR/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

# Source the user's venv
source "$WORK_DIR/venv/bin/activate"

# Export PYTHONPATH (Critical for finding 'upside_engine' and 'rotamer_parameter_estimation')
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# --- Configuration ---
# Point to the subdirectory (test_00) where init wrote the file
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"
steps=50
# ---------------------

echo "Running Zero Descent with: python3 (from venv)"
echo "Target Checkpoint: $checkpoint"
echo "Steps: $steps"

# Run the zero descent script
python3 -u ConDiv_zero.py "$checkpoint" "$steps" 2>&1 | tee -a "$WORK_DIR/run_zero.output"
