#!/bin/bash

echo "Setting up universal Python environment for UPSIDE2..."

# Create virtual environment in the root directory
python3 -m venv .venv

# Activate virtual environment and install packages
source .venv/bin/activate
pip install --upgrade pip 
pip install --upgrade h5py tables matplotlib mdtraj pymbar jax pandas
deactivate

echo "Python environment setup completed!"
echo "Virtual environment created at: $(pwd)/.venv"
echo "To activate the environment, run: source .venv/bin/activate"

# Update source.sh to include virtual environment activation
if [ -f "source.sh" ]; then
    echo "" >> source.sh
    echo "# Activate Python virtual environment if it exists" >> source.sh
    echo "if [ -d \"\$UPSIDE_HOME/.venv\" ]; then" >> source.sh
    echo "    source \$UPSIDE_HOME/.venv/bin/activate" >> source.sh
    echo "    echo \"Python virtual environment activated\"" >> source.sh
    echo "else" >> source.sh
    echo "    echo \"Warning: Python virtual environment not found at \$UPSIDE_HOME/.venv\"" >> source.sh
    echo "    echo \"Run 'bash setup_python_env.sh' to create the environment\"" >> source.sh
    echo "fi" >> source.sh
    
    echo "Updated source.sh with virtual environment activation"
else
    echo "Warning: source.sh not found - virtual environment activation not added"
fi
