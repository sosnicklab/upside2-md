#!/bin/bash

# Define current directory variables
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

# Path math
PROJECT_ROOT="$WORK_DIR/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

# Export PYTHONPATH
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# Setup Run Parameters
mode=restart

# --- FIX: Point to the subdirectory (test_00) where init wrote the file ---
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"
# ------------------------------------------------------------------------

step=40

# Run using the VENV python
echo "Running ConDiv with: $VENV_PYTHON"
echo "Checkpoint: $checkpoint"

$VENV_PYTHON ConDiv.py $mode $checkpoint $step | tee -a "$WORK_DIR/run.output"