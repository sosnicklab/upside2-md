# Task Plan

## Project Goal
Run and debug the MARTINI 3.0 8-stage bilayer equilibration workflow at `example/16.MARTINI/run_sim_bilayer_8stage.sh`.

## Architecture & Key Decisions
- **Problem 1**: `prepare_martini.py` had a hardcoded output directory (`outputs/martini_test/`) that didn't match the shell script's expected directory (`outputs/martini_8stage/`).
- **Solution**: Modified `prepare_martini.py` to accept an optional output directory argument and updated the shell script to pass the correct directory.

## Execution Phases
- [x] **Phase 1: Initial Setup**
  - [x] Check project activation and onboarding status
  - [x] Read the `run_sim_bilayer_8stage.sh` script
  - [x] Check project structure and dependencies
  - [x] Verify C++ core was compiled
- [x] **Phase 2: Fix Output Directory Issue**
  - [x] Modify `prepare_martini.py` to accept output directory parameter
  - [x] Update `run_sim_bilayer_8stage.sh` to pass output directory
- [x] **Phase 3: Run Workflow**
  - [x] Run stages 0-6.1 (successfully completed)
  - [x] Troubleshoot stage 6.2 (completed after long run)
  - [x] Run stage 6.3 manually (completed)
  - [x] Run stages 6.4-6.6 manually (completed)
  - [x] Run stage 7.0 (Production with Parrinello-Rahman, completed)
- [x] **Phase 4: Update Logs**
  - [x] Update stage 6.4-7.0 logs with actual output

## Known Errors / Blockers
- **File name too long**: Occurred when processing stage 6.3. Fixed by running the stage manually.
