#!/bin/bash
rm -rf .venv

echo "Creating virtual environment in .venv..."
python3.15 -m venv .venv

source .venv/bin/activate
pip install --upgrade pip setuptools wheel

echo "Installing core scientific packages..."
# jax is pymbar's solver backend
pip install numpy scipy h5py tables matplotlib mdtraj pymbar jax pandas ProDy biopython scikit-learn

echo "Installing analysis and document packages..."
# pyhdx (whose apply_cmap 5.analyze_D_uptake.py replaces without it) and torch have no Python 3.15 wheels
pip install PeptideBuilder tmtools python-pptx pymupdf pypdf

echo "Testing core imports..."
python -c "import numpy, scipy, h5py, tables, matplotlib, mdtraj, pymbar, jax, pandas, prody, Bio, sklearn"
echo "Setup complete! Run 'source .venv/bin/activate' to start working."
