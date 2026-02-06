# Task Plan

## Project Goal
Formalize the workflow logic of `run_sim_bilayer.sh` to ensure each simulation stage uses a separate `.up` file with appropriate parameters (including softened/hard potentials) and coordinates passed from previous stages.

## Architecture & Key Decisions
- **Per-Stage .up Files**: Each stage will have its own input .up file generated at the start of the stage
- **Parameter Inheritance**: Each stage's .up file will inherit coordinates from previous stage's output
- **Explicit Parameterization**: All stage-specific parameters will be clearly defined in the shell script
- **HDF5 Attribute Modification**: Python will be used to modify HDF5 attributes directly when generating stage files

## Execution Phases
- [x] **Phase 1: Modify prepare_martini.py**
  - [x] Add support for generating stage-specific .up files
  - [x] Add parameters for softening (lj_alpha, slater_alpha)
  - [x] Add support for barostat configuration (Berendsen vs Parrinello-Rahman)
  - [x] Improve error handling and validation
- [x] **Phase 2: Modify run_sim_bilayer.sh**
  - [x] Refactor to use per-stage .up files
  - [x] Add explicit parameterization for each stage
  - [x] Improve error handling and validation
  - [x] Add stage-specific comments and documentation
- [x] **Phase 3: Verify Workflow**
  - [x] Run simulation to ensure all stages complete
  - [x] Check per-stage .up files are created correctly
  - [x] Verify softening parameters and barostat types are correctly set

## Known Errors / Blockers
- **HDF5 File Handling**: Ensure that the HDF5 library is available in the Python environment
- **Permissions**: Ensure the script has execute permissions
- **Stage Parameter Passing**: Ensure that stage-specific parameters are correctly passed to prepare_martini.py
