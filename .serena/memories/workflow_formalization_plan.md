# Workflow Formalization Plan for run_sim_bilayer.sh

## Project Goal
Formalize the workflow logic of `run_sim_bilayer.sh` to ensure each simulation stage uses a separate `.up` file with appropriate parameters (including softened/hard potentials) and coordinates passed from previous stages.

## Current Workflow Analysis
The current `run_sim_bilayer.sh` implements:
1. **Stage 1: Prepare Input Files** - Generates initial .up file with MARTINI 3.0 parameters using `prepare_martini.py`
2. **Stage 2: Energy Minimization** - Runs gradient descent minimization with softened potentials
3. **Stage 3: NPT Equilibration (Softened Potentials, Berendsen)** - Equilibrates system with softened potentials
4. **Stage 4: NPT Equilibration (Reduced Softening, Berendsen)** - Reduces softening parameters
5. **Stage 5: NPT Production (Hard, Parrinello-Rahman)** - Runs production simulation with hard potentials
6. **Stage 6: VTF Generation** - Extracts trajectories for visualization

## Key Issues to Address
1. The workflow uses environment variables for parameter passing, which is not explicit
2. Parameters are not clearly separated between stages
3. The current implementation modifies files in-place rather than generating new .up files for each stage
4. Stage-specific parameters (like softening) are not properly documented or validated

## Architecture Decisions
1. **Per-Stage .up Files**: Each stage will have its own input .up file generated at the start of the stage
2. **Parameter Inheritance**: Each stage's .up file will inherit coordinates from previous stage's output
3. **Explicit Parameterization**: All stage-specific parameters will be clearly defined in the shell script
4. **HDF5 Attribute Modification**: Python will be used to modify HDF5 attributes directly when generating stage files

## Implementation Phases

### Phase 1: Modify prepare_martini.py
- Add support for generating stage-specific .up files
- Add parameters for softening (lj_alpha, slater_alpha)
- Add support for barostat configuration (Berendsen vs Parrinello-Rahman)
- Improve error handling and validation

### Phase 2: Modify extract_martini_vtf.py
- Ensure VTF extraction works with per-stage .up files
- Improve box dimension handling from simulation output
- Add support for extracting trajectories from specific stages

### Phase 3: Modify src/main.cpp
- Ensure stage-specific parameters are properly read from .up files
- Add support for Berendsen and Parrinello-Rahman barostats
- Improve minimization and equilibration logic

### Phase 4: Modify src/martini.cpp
- Ensure MartiniPotential reads softening parameters from HDF5 attributes
- Improve spline generation for softened potentials
- Add support for stage-specific potential parameters

### Phase 5: Modify src/box.cpp
- Improve NPT barostat implementation
- Add support for Berendsen and Parrinello-Rahman barostats
- Improve box dimension update logic

### Phase 6: Update run_sim_bilayer.sh
- Refactor to use per-stage .up files
- Add explicit parameterization for each stage
- Improve error handling and validation
- Add stage-specific comments and documentation

## Expected Changes
- Each stage will have its own .up file (e.g., bilayer.minimized.up, bilayer.npt_equil.up, etc.)
- Parameters like softening will be clearly defined for each stage
- Coordinates will be passed from previous stage's output
- The workflow will be more modular and maintainable
