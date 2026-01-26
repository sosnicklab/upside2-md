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
step=152  # Run exactly 76 steps and then QUIT.

# --- ROBUST RESUME LOGIC (NO ROLLBACK) ----------------------------------
# 1. Start with initial as fallback
checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"

# 2. Get ALL epoch folders, sorted reverse (epoch_01 before epoch_00)
#    We use 'sort -r' to guarantee the correct alphabetical order.
all_dirs=$(ls -d "$WORK_DIR/test_00"/epoch_*_minibatch_* 2>/dev/null | sort -r)

# 3. Find the most recent VALID checkpoint
#    This loop prevents rollback by skipping empty/corrupt folders at the end.
for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
        echo "--> Found Valid Resume Point: $d"
        checkpoint="$d/checkpoint.pkl"
        break  # Stop searching immediately once found
    fi
done
# ------------------------------------------------------------------------

# Run using the VENV python
echo "Running ConDiv with: python3 (from venv)"
echo "Loading Checkpoint: $checkpoint"

# Run once. Stop sharp.
python3 -u ConDiv.py $mode "$checkpoint" $step 2>&1 | tee -a "$WORK_DIR/run.output"