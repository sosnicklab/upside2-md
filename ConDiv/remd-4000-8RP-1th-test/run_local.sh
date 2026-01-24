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
step=76  # Run exactly 40 batches, then STOP.

# --- RESUME LOGIC (Fixes Rollback) --------------------------------------
# 1. Default to the initial checkpoint (Epoch 0)
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"

# 2. Get all epoch folders, sorted by Time (newest first)
#    ls -t orders by modification time, so we check the most recent attempts first.
all_dirs=$(ls -td "$WORK_DIR/test_00"/epoch_*_minibatch_* 2>/dev/null)

# 3. Find the first folder that actually contains a checkpoint file
for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
        echo "--> Resuming from: $d"
        checkpoint="$d/checkpoint.pkl"
        break  # We found the latest valid save. Stop searching.
    fi
done
# ------------------------------------------------------------------------

# Run using the VENV python
echo "Running ConDiv with: python3 (from venv)"
echo "Loading Checkpoint: $checkpoint"

python3 -u ConDiv.py $mode "$checkpoint" $step 2>&1 | tee -a "$WORK_DIR/run.output"