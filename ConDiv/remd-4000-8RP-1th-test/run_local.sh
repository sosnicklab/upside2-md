#!/bin/bash

# Define current directory variables
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

# Path math: Assuming we are in projectroot/ConDiv/folder
# We need to reach projectroot/obj and projectroot/py
PROJECT_ROOT="$WORK_DIR/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

# Check if venv python exists
if [ ! -f "$VENV_PYTHON" ]; then
    echo "Error: Virtual environment not found at $VENV_PYTHON"
    echo "Please run setup_venv.sh first."
    exit 1
fi

# Export PYTHONPATH
# 1. WORK_DIR: for ConDiv.py itself
# 2. UPSIDE_LIB: for upside_engine.so
# 3. UPSIDE_PY: for rotamer_parameter_estimation.py (projectroot/py)
# 4. PROJECT_ROOT/src: (Optional) if you saved source files there instead of py
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# Setup Run Parameters
mode=restart
checkpoint="$WORK_DIR/initial_checkpoint.pkl"
step=40

# Run using the VENV python
echo "Running ConDiv with: $VENV_PYTHON"
echo "PYTHONPATH: $PYTHONPATH"

$VENV_PYTHON ConDiv.py $mode $checkpoint $step | tee -a "$WORK_DIR/run.output"