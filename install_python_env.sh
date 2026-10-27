#!/bin/bash
rm -rf .venv

echo "Creating virtual environment in .venv..."
python3.14 -m venv .venv

source .venv/bin/activate
pip install --upgrade pip setuptools wheel

echo "Installing core scientific packages..."
# jax is pymbar's solver backend
pip install numpy scipy h5py tables matplotlib mdtraj pymbar jax pandas ProDy biopython scikit-learn

echo "Installing analysis and document packages..."
# pydssp brings torch; pyhdx.plot needs proplot, which pins matplotlib 3.4, so 5.analyze_D_uptake.py uses its own apply_cmap
pip install PeptideBuilder tmtools python-pptx pymupdf pypdf pydssp pyhdx==0.4.3 "hdxms-datasets<0.2" colorcet

echo "Testing imports..."
python -c "import numpy, scipy, h5py, tables, matplotlib, mdtraj, pymbar, jax, pandas, prody, Bio, sklearn"
python -c "import PeptideBuilder, tmtools, pptx, pymupdf, pypdf, torch, pydssp, pyhdx, hdxms_datasets, colorcet"
echo "Setup complete! Run 'source .venv/bin/activate' to start working."
