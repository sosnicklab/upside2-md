#!/bin/bash
rm -rf .venv

echo "Creating virtual environment in .venv..."
python3.11 -m venv .venv

source .venv/bin/activate
pip install --upgrade pip setuptools wheel

echo "Installing core scientific packages..."
# jax is pymbar's solver backend
pip install numpy scipy h5py tables matplotlib mdtraj pymbar jax pandas ProDy biopython scikit-learn

echo "Installing optional packages..."
pip install pyhdx==0.4.3 "hdxms-datasets<0.2"

echo "Testing core imports..."
python -c "import numpy, scipy, h5py, tables, matplotlib, mdtraj, pymbar, jax, pandas, prody, Bio, sklearn"
echo "Setup complete! Run 'source .venv/bin/activate' to start working."
