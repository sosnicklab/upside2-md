#!/bin/bash
rm -rf .venv

echo "Creating virtual environment in .venv..."
python3.11 -m venv .venv

source .venv/bin/activate
pip install --upgrade pip setuptools wheel

echo "Installing core scientific packages..."
pip install h5py tables matplotlib mdtraj pymbar pandas ProDy scikit-learn jax

echo "Installing optional packages..."
pip install colorcet pyhdx==0.4.3 "hdxms-datasets<0.2"

echo "Testing core imports..."
python -c "import h5py, tables, matplotlib, mdtraj, pymbar, pandas, prody, sklearn, jax"
echo "Setup complete! Run 'source .venv/bin/activate' to start working."