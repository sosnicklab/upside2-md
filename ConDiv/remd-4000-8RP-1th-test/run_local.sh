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

# Export PYTHONPATH
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# Setup Run Parameters
mode=restart
step=40

# --- AUTO-RESUME LOGIC --------------------------------------------------
# 1. Default to the initial checkpoint (Start from 0)
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"

# 2. Look for the latest "epoch_XX_minibatch_XX" folder
#    'sort' ensures epoch_01 comes after epoch_00, and 'tail -n 1' grabs the last one.
latest_dir=$(ls -d "$WORK_DIR/test_00"/epoch_*_minibatch_* 2>/dev/null | sort | tail -n 1)

# 3. If a later folder exists AND contains a checkpoint, use it instead.
if [ ! -z "$latest_dir" ]; then
    if [ -f "$latest_dir/checkpoint.pkl" ]; then
        echo "--> Found Resume Point: $latest_dir"
        checkpoint="$latest_dir/checkpoint.pkl"
    else
        echo "--> Latest folder empty or corrupt. Falling back to previous or initial."
    fi
fi
# ------------------------------------------------------------------------

# Run using the VENV python
echo "Running ConDiv with: python3 (from venv)"
echo "Loading Checkpoint: $checkpoint"

# NOTE: Added quotes around "$checkpoint" to handle paths safely
python3 -u ConDiv.py $mode "$checkpoint" $step 2>&1 | tee -a "$WORK_DIR/run.output"