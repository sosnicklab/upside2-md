# Task Plan

## Project Goal
- Remove the NPT box-dimension debug output from the C++ runtime so simulation logs no longer emit messages such as:
  - `[NPT] t 0.500 box 111.93 111.93 110.20`

## Architecture & Key Decisions
- Keep the NPT barostat behavior unchanged.
- Remove only the box-dimension debug `printf` calls in `src/box.cpp`.
- Keep the non-finite-pressure warning path, since that is still useful runtime diagnostics and is not just box-dimension spam.
- Do not change the HDF5/debug wiring unless it becomes dead or required for correctness; this cleanup should be minimal-impact.

## Execution Phases
- [x] Phase 1: Locate all NPT box-dimension debug prints in the C++ runtime.
- [x] Phase 2: Remove the box-dimension debug output with minimal source changes.
- [x] Phase 3: Rebuild and verify the cleanup.

## Known Errors / Blockers
- No blocker identified.

## Review
- Removed the two NPT box-dimension debug messages from `src/box.cpp`:
  - the barostat registration message no longer prints the initial box dimensions,
  - the per-update `[NPT] t ... box ...` message was removed entirely.
- Kept the non-finite-pressure warning path unchanged.
- Verification passed:
  - `rg -n "\\[NPT\\].*box|box %.2f %.2f %.2f|box %.2f x %.2f x %.2f" src/box.cpp src/main.cpp`
  - `source .venv/bin/activate && source source.sh && cmake --build obj -j4`
