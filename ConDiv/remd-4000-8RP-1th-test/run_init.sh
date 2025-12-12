#!/bin/bash

# Define paths
X=$(pwd)
Y=test_00

# Ensure venv matches your setup
VENV_PYTHON="$X/venv/bin/python3"
PROJECT_ROOT="$X/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

# Export PYTHONPATH
export PYTHONPATH="$X:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

# Initialization Parameters
mode=initialize
initial_param=init_param       # Folder with initial parameters (Must exist)
upside_input_dir=upside_input  # Folder with protein inputs (Must exist)
pdb_list=pdb_list              # File with list of PDBs (Must exist)

echo "Running Initialization..."
$VENV_PYTHON ConDiv.py $mode $initial_param $upside_input_dir $pdb_list $X/$Y
