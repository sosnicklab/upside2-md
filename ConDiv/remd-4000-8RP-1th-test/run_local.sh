#!/bin/bash

# Define current directory variables
WORK_DIR=$(pwd)
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

# Assuming your compiled engine is in 'obj' relative to this folder,
# or typically one level up if ConDiv is in a subfolder.
# Adjust "../obj" if your object folder is elsewhere.
UPSIDE_LIB="$WORK_DIR/../obj" 

# Check if venv python exists
if [ ! -f "$VENV_PYTHON" ]; then
    echo "Error: Virtual environment not found at $VENV_PYTHON"
    echo "Please run setup_venv.sh first."
    exit 1
fi

# Export PYTHONPATH so Python finds the C++ extension and local src
# We add WORK_DIR (for ConDiv imports) and UPSIDE_LIB (for the engine)
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$PYTHONPATH"

# Setup Run Parameters
# If 'test_00' logic is required, adjust path variables below. 
# Here we assume we are running IN the folder containing the checkpoint.
mode=restart
checkpoint="$WORK_DIR/initial_checkpoint.pkl"
step=40

# Run using the VENV python
echo "Running ConDiv with: $VENV_PYTHON"
echo "PYTHONPATH: $PYTHONPATH"

$VENV_PYTHON ConDiv.py $mode $checkpoint $step | tee -a "$WORK_DIR/run.output"