#!/bin/bash

# Get the directory where this script is running
CURRENT_DIR=$(pwd)
VENV_DIR="$CURRENT_DIR/venv"

echo "Setting up Python virtual environment in $VENV_DIR..."

# 1. Create venv if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    # Try to use python3.10 if available, otherwise fallback to system python3
    if command -v python3.10 &> /dev/null; then
        echo "Using Python 3.10..."
        python3.10 -m venv "$VENV_DIR"
    else
        echo "Using default Python 3..."
        python3 -m venv "$VENV_DIR"
    fi
    echo "Virtual environment created."
else
    echo "Virtual environment already exists."
fi

# 2. Activate venv
source "$VENV_DIR/bin/activate"

# 3. Upgrade pip and build tools
pip install --upgrade pip setuptools wheel

# 4. Set compilers and paths for Apple Silicon (M1/M2)
# This fixes the missing header errors (lzo1x.h, etc.)
export HDF5_DIR=$(brew --prefix hdf5)
export LZO_DIR=$(brew --prefix lzo)
export BLOSC_DIR=$(brew --prefix c-blosc)
export BZIP2_DIR=$(brew --prefix bzip2)

# Flags to help the compiler find the headers
export CFLAGS="-I$HDF5_DIR/include -I$LZO_DIR/include -I$BLOSC_DIR/include -I$BZIP2_DIR/include"
export LDFLAGS="-L$HDF5_DIR/lib -L$LZO_DIR/lib -L$BLOSC_DIR/lib -L$BZIP2_DIR/lib"

# 5. Install dependencies
echo "Installing dependencies..."
# Install Cython and Numpy first as they are build deps for tables
pip install cython numpy
# Install tables explicitly with the flags
HDF5_DIR=$HDF5_DIR pip install tables
# Install remaining deps
# prody is required by py/PDB_to_initial_structure.py (used by input conversion)
pip install scipy torch ipython prody

# 6. Clean up
find . -name "*.pyc" -delete

echo ""
echo "----------------------------------------------------------------"
echo "Setup complete!"
echo "To run your training, use:"
echo "  $VENV_DIR/bin/python3 ConDiv.py ..."
echo "----------------------------------------------------------------"
