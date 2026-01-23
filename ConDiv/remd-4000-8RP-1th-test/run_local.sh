#!/bin/bash

# Define current directory variables
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

# Path math
PROJECT_ROOT="$WORK_DIR/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

# Source the user's venv as requested
source "$WORK_DIR/venv/bin/activate"

# Export PYTHONPATH
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# Setup Run Parameters
mode=restart

# --- FIX: Point to the subdirectory (test_00) where init wrote the file ---
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"
latest_dir=$(ls -d "$WORK_DIR"/test_00/epoch_*_minibatch_* 2>/dev/null | sort | tail -n 1)
if [ ! -z "$latest_dir" ] && [ -f "$latest_dir/checkpoint.pkl" ]; then
    echo "Found newer checkpoint in: $latest_dir"
    checkpoint="$latest_dir/checkpoint.pkl"
fi
# ------------------------------------------------------------------------

step=40

# Run using the VENV python
echo "Running ConDiv with: python3 (from venv)"
echo "Checkpoint: $checkpoint"

python3 -u ConDiv.py $mode $checkpoint $step 2>&1 | tee -a "$WORK_DIR/run.output"