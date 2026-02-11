# Task: Add MARTINI backbone rigid hold

## 1. Project Goal
Implement a rigid hold for a group of particles (MARTINI protein backbone beads) in `src/martini.cpp`, expose a calling interface in `src/main.cpp`, and apply the rigid hold in `example/16.MARTINI/run_sim_1ubq.sh` before production (including minimization and soft-particle stages).

## 2. Architecture & Key Decisions
* **Scope**: Add a MARTINI-specific rigid-hold feature without changing existing simulation architecture.
* **Interface**: Expose a minimal CLI/config entry point in `src/main.cpp` to activate the MARTINI backbone rigid hold, using explicit BB indices from input when available.
* **Workflow**: Enable the hold in the example script prior to production, covering minimization and soft-particle stages.

## 3. Execution Phases
- [x] Phase 1: Inspect current MARTINI handling in `src/martini.cpp` and CLI/config wiring in `src/main.cpp`.
- [x] Phase 2: Implement backbone rigid hold and plumb the interface.
- [x] Phase 3: Update `example/16.MARTINI/run_sim_1ubq.sh` to enable the hold for minimization and soft-particle runs.

## 4. Known Errors / Blockers
- None yet.
