# Task: Audit and tighten rigid-hold energy coupling behavior in C++

## 1. Project Goal
Verify rigid-hold behavior in `src/*` and address the likely energy-coupling path where held atoms can indirectly drive non-held atoms through the NPT pressure path.

## 2. Architecture & Key Decisions
- Scope:
  - Inspect rigid-hold enforcement in `src/martini.cpp`, `src/deriv_engine.cpp`, `src/main.cpp`, and `src/box.cpp`.
  - Keep interface and file format behavior unchanged.
- Key decisions:
  - Keep existing rigid-hold semantics (held atoms still interact unless explicitly decoupled).
  - Patch only the NPT pressure estimator so constrained DOFs do not feed barostat scaling.
- Revised Decisions:
  - None yet.

## 3. Execution Phases
- [x] Phase 1: Trace rigid-hold registration and enforcement in minimization/MD.
- [x] Phase 2: Trace force, integration, and barostat call order for coupling paths.
- [x] Phase 3: Implement minimal `src/box.cpp` fix to exclude constrained DOFs from pressure estimation.
- [x] Phase 4: Build and run focused verification.
- [x] Phase 5: Record results in `progress.md`.

## 4. Known Errors / Blockers
- None currently.
