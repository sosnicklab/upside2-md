
# Code Structure

## Root Directory
- **src/**: C++ core simulation engine
- **py/**: Python scripts for configuration and analysis
- **example/**: Tutorials and examples for various use cases
- **parameters/**: Force field parameters
- **cmake/**: CMake modules
- **obj/**: Build output directory (contains executable and shared library)
- **inputs/**: Input files for simulations
- **outputs/**: Output files from simulations

## C++ Core Files
- **main.cpp**: Entry point for upside executable
- **deriv_engine.h/.cpp**: Core simulation engine managing differentiable computation nodes
- **monte_carlo_sampler.h/.cpp**: Implements Monte Carlo moves (pivot moves)
- **thermostat.h/.cpp**: Controls temperature in NVT simulations
- **h5_support.h/.cpp**: Handles HDF5 input/output
- **membrane_potential.cpp**: Membrane interactions
- **martini.cpp**: MARTINI force field integration

## Python Scripts
- **upside_config.py**: Main configuration script for creating .up simulation files
- **PDB_to_initial_structure.py**: Processes PDB files into Upside-compatible structures
- **run_upside.py**: High-level script for running simulations
- **mdtraj_upside.py**: Integrates with MDTraj for trajectory analysis
- **upside_engine.py**: Python wrapper around C++ core library
- **upside_martini.py**: MARTINI force field integration
