# Task Plan

## Project Goal
Fix the issue where softened particles weren't switching to hard particles after minimization in `run_sim_bilayer.sh`, and address the resulting system instability.

## Architecture & Key Decisions
- **Problem 1**: Potential parameters were stored in HDF5 file attributes and not being updated during simulation, causing softened particles to remain active in subsequent stages.
- **Solution**: Added Python code to directly modify HDF5 file attributes after each stage.
- **Problem 2**: Sudden switch from fully softened to hard particles caused extreme potential values and system instability.
- **Solution**: Implemented progressive softening reduction with two NPT equilibration stages.

## Execution Phases
- [x] **Phase 1: Identify Root Cause**
  - [x] Examine run_sim_bilayer.sh script
  - [x] Check how prepare_martini.py stores parameters in HDF5 file
- [x] **Phase 2: Fix Softening Switching**
  - [x] Add Python code to modify HDF5 attributes
  - [x] Test the fix
- [x] **Phase 3: Address System Instability**
  - [x] Implement progressive softening reduction (two NPT equilibration stages)
  - [x] Update script with new stages
- [x] **Phase 4: Verify Fix**
  - [x] Run simulation to ensure stable potential values
  - [x] Check logs for softening parameter updates

## Known Errors / Blockers
- **System Instability After Initial Fix**: Caused by sudden switch from fully softened to hard particles. Fixed by progressive softening reduction.
