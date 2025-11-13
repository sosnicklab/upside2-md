# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

## Common Commands

### Installation & Build
```bash
# Install dependencies and compile C++ core
./install.sh

# Alternative: Manual build
source source.sh  # Sets up environment variables
rm -rf obj/* && cd obj && cmake ../src/ && make
```

### Simulation Workflow
```bash
# Convert PDB to initial structure files
upside/py/PDB_to_initial_structure.py input.pdb output_basename

# Create simulation configuration (.up file)
upside/py/upside_config.py --output simulation.up --fasta output_basename.fasta ...

# Run constant temperature simulation
upside/obj/upside --duration 1e7 --temperature 0.85 simulation.up

# Run replica exchange simulation
upside/obj/upside --duration 1e7 --temperature 0.5,0.53,0.56 --swap-set 0-1 --swap-set 1-2 config_0.up config_1.up config_2.up
```

### Analysis
```bash
# List contents of .up file
upside/py/attr_overview.py simulation.up

# Convert to VMD format
upside/py/extract_vtf simulation.up simulation.vtf

# Load in Python with MDTraj
import sys
sys.path.append('upside/py')
import mdtraj_upside as mu
traj = mu.load_upside_traj('simulation.up')
```

## Code Architecture

### C++ Core (src/)
The C++ core handles the simulation engine and numerical computations:

- **main.cpp**: Entry point for the `upside` executable, parses command line arguments using TCLAP.
- **DerivEngine**: Core simulation engine that manages differentiable computation nodes (DerivComputation) to calculate potentials and derivatives.
- **MonteCarloSampler**: Implements Monte Carlo moves (like pivot moves) to enhance sampling.
- **Thermostat**: Controls temperature in NVT simulations.
- **H5Support**: Handles HDF5 input/output for simulation configurations and trajectories.
- **OpenMP**: Used for parallelism (replica exchange, force calculations).

Key files:
- src/main.cpp
- src/deriv_engine.h/.cpp
- src/monte_carlo_sampler.h/.cpp
- src/thermostat.h/.cpp
- src/h5_support.h/.cpp

### Python Layer (py/)
The Python layer provides user-friendly interfaces for configuration and analysis:

- **upside_config.py**: Main configuration script for creating .up simulation files.
- **PDB_to_initial_structure.py**: Processes PDB files into Upside-compatible initial structures.
- **run_upside.py**: High-level script for running simulations.
- **mdtraj_upside.py**: Integrates with MDTraj for trajectory analysis and visualization.
- **upside_engine.py**: Python wrapper around the C++ core library.

Key files:
- py/upside_config.py
- py/PDB_to_initial_structure.py
- py/run_upside.py
- py/mdtraj_upside.py

### Examples
The `example/` directory contains tutorials for various use cases:
- Getting started with basic simulations
- Replica exchange simulation
- Trajectory analysis
- Membrane simulations
- Pulling simulations
- MARTINI force field integration

## Directory Structure
```
/
├── src/               # C++ core source code
├── py/                # Python scripts
├── example/           # Tutorials and examples
├── parameters/        # Force field parameters
├── cmake/             # CMake modules
├── obj/               # Build output directory
├── README.md          # Main documentation
├── install.sh         # Installation script
└── source.sh          # Environment setup script
```

## Dependencies
- **C++**: CMake 2.8+, C++11 compiler, HDF5 1.8+, Eigen 3.0+
- **Python**: Python 3.7+, Numpy, Scipy, PyTables, Prody, MDTraj (optional)
- for compling this project, run install_M1.sh at project root
- for running any script in this project, source source.sh and .venv/bin/activate at project root