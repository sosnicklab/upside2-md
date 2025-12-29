#!/bin/bash
set -e
WORK_DIR=$(pwd)
# Source the user's venv as requested
source "$WORK_DIR/venv/bin/activate"

PROJECT_ROOT="$WORK_DIR/../.."
export PYTHONPATH="$WORK_DIR:$PROJECT_ROOT/obj:$PROJECT_ROOT/py:$PROJECT_ROOT/src:$PYTHONPATH"

# Clean previous run
rm -rf "$WORK_DIR/test_debug"

# 1. Re-initialize to generate the smaller dataset
echo "Initializing Debug Run..."
python3 debug_ConDiv.py initialize init_param upside_input pdb_list "$WORK_DIR/test_debug"

# 2. Run the debug loop
checkpoint="$WORK_DIR/test_debug/initial_checkpoint.pkl"
echo "Running Debug Loop..."
python3 debug_ConDiv.py restart "$checkpoint" 5
