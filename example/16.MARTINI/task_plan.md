# Task Plan

## Project Goal
- Run `run_sim_1ubq.sh` successfully and debug failures until the workflow can start the simulation.

## Architecture & Key Decisions
- Use the existing shell workflow and upstream Upside scripts rather than rewriting the run pipeline.
- Follow project environment requirements by sourcing `.venv/bin/activate` and `source.sh` before execution.
- Apply minimal, local fixes needed to unblock startup while preserving current simulation setup.

## Execution Phases
- [x] Phase 1: Reproduce the failure by running `run_sim_1ubq.sh` and capture the first blocking error.
- [x] Phase 2: Identify root cause from script/config/runtime outputs.
- [x] Phase 3: Implement minimal fix and re-run.
- [x] Phase 4: Verify the simulation starts (initial run enters startup stage without immediate fatal error).
- [x] Phase 5: Document outcomes and any remaining blockers.

## Known Errors / Blockers
- Previous blocker fixed: `set_initial_position.py` crashed when `/output/pos` existed but had zero frames.
- Non-blocking runtime warnings remain in logs from HDF5 close diagnostics.

## Revised Decisions
- None yet.
