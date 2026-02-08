# Project: 16.MARTINI (Upside Molecular Dynamics Simulation)

## Overview
This project implements a MARTINI force field based molecular dynamics simulation workflow for lipid bilayers and water systems using the Upside MD package.

## Key Commands

### Simulation Setup and Execution
- `run_sim_bilayer.sh`: Run bilayer simulation
- `run_sim_bilayer_8stage.sh`: Run 8-stage bilayer equilibration
- `run_sim_water.sh`: Run water simulation
- `run_sim.sh`: Generic simulation runner
- `test_run.sh`: Test run script

### Preparation and Analysis
- `prepare_martini.py`: Generate Upside simulation input files (.up) for MARTINI force field
- `set_initial_position.py`: Set initial positions in simulation files
- `extract_martini_vtf.py`: Extract trajectory files in VTF format
- `scan_water_diffusion.py`: Analyze water diffusion from simulation data
- `plot_diffusion_results.py`: Plot diffusion analysis results

### Optimization
- `optimize_interaction_table.py`: Optimize MARTINI interaction tables

## Environment Setup
- Python scripts require Upside MD package and HDF5 library
- Shell scripts should be run from project root

## Example Workflow
1. Run water simulation: `./run_sim_water.sh`
2. Analyze diffusion: `python scan_water_diffusion.py`
3. Plot results: `python plot_diffusion_results.py`
