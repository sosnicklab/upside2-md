#!/bin/bash
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"
PROJECT_ROOT="$WORK_DIR/../.."
export PYTHONPATH="$WORK_DIR:$PROJECT_ROOT/obj:$PROJECT_ROOT/py:$PROJECT_ROOT/src:$PYTHONPATH"

rm -rf "$WORK_DIR/test_debug"

# Use the debug script
# 1. Re-initialize to generate the smaller dataset
echo "Initializing Debug Run..."
$VENV_PYTHON debug_ConDiv.py initialize init_param upside_input pdb_list $WORK_DIR/test_debug

# 2. Run the debug loop
checkpoint="$WORK_DIR/test_debug/initial_checkpoint.pkl"
echo "Running Debug Loop..."
$VENV_PYTHON debug_ConDiv.py restart $checkpoint 5
