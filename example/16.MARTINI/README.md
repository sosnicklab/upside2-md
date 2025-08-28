# MARTINI 3.0 Optimized Simulation Workflow

This directory contains a clean, optimized workflow for running MARTINI 3.0 coarse-grained molecular dynamics simulations with UPSIDE.

## Essential Files

### Core Workflow Scripts
- **`prepare_martini.py`** - Prepares the initial HDF5 input file from PDB structure
- **`optimize_interaction_table.py`** - Optimizes the interaction table for memory efficiency
- **`run_martini_optimized.sh`** - Main simulation script with optimization
- **`extract_martini_vtf.py`** - Extracts trajectory to VTF format for visualization

### Input/Output Directories
- **`pdb/`** - Contains input PDB structures and topology files
- **`outputs/`** - Contains simulation results and optimized files
- **`ff3.00/`** - MARTINI 3.0 force field parameters

## Quick Start

1. **Prepare the system:**
   ```bash
   source ../../source.sh
   python prepare_martini.py
   ```

2. **Run optimized simulation (includes VTF generation):**
   ```bash
   ./run_martini_optimized.sh
   ```

3. **Manual trajectory extraction (if needed):**
   ```bash
   python extract_martini_vtf.py outputs/martini_test/test.run.optimized.up trajectory.vtf
   python extract_martini_vtf.py outputs/martini_test/test.run.optimized.up trajectory.pdb
   ```

## Features

- **Memory optimization**: ~87.5% reduction in interaction table size
- **Complete MARTINI 3.0 support**: Bonds, angles, dihedrals, and non-bonded interactions
- **Proper unit conversions**: Correctly derived from energy conversion factor (2.914952774272 kJ/mol → E_up)
- **Additional exclusions**: Proper handling of ITP file exclusions
- **Protein-lipid systems**: Full support for complex biomolecular systems
- **Automatic VTF/PDB generation**: Visualization files created automatically after simulation
- **Comprehensive workflow**: Single command completes preparation, optimization, simulation, and visualization

## Backup Files

All development, debugging, and analysis scripts have been moved to the `bak/` directory for reference.

## Requirements

- UPSIDE installation with MARTINI support
- Python 3.x with h5py, numpy, tables
- MARTINI 3.0 force field files
