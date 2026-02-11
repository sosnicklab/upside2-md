# Task: Debug NPT + Ewald explosion in MARTINI bilayer example

## 1. Objective
Identify and fix the source of pressure/particle blow-up in `example/16.MARTINI/run_sim_bilayer.sh` between commits `9829d3a70f0e6cee8e7f63a454179b4de87718e4` and `b819b310fcf400e9a815439bfe2b02d0d88cec93`, using correct physical unit conversions (no parameter fudging).

## 2. Architecture & Key Decisions
* **Scope of change**: Keep Ewald implementation and NPT logic intact; only adjust code where it is incorrect or mismatched to physical units.
* **Debug focus**: Compare behavior across the two commits, inspect Ewald real/reciprocal/self terms, pressure tensor contributions, and unit conversions.
* **Evidence-driven**: Use logs, diffs, and derived unit checks to pinpoint the regression instead of tuning parameters.

## 3. Execution Phases
- [x] Phase 1: Collect evidence (diffs, logs, relevant config inputs).
- [x] Phase 2: Trace unit conversions and pressure/force contributions for NPT+Ewald.
- [x] Phase 3: Implement minimal fix with unit-consistent logic.
- [x] Phase 4: Validate with the bilayer workflow and document results.

## 4. Known Errors / Blockers
- None yet.
