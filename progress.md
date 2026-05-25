# Progress Log

## 2026-05-24 MARTINI Branch Cleanup
- Actions taken:
  - Restored low-risk shared-file changes toward master style and removed unused
    MARTINI placeholders.
  - Split the MARTINI C++ implementation into focused subsystem files and added
    a private internal header for cross-file declarations.
  - Removed generated dry-MARTINI HDF5 artifacts from source control scope and
    added an ignore rule for regenerated copies.
  - Added a public MARTINI table-builder facade used by the CG-lipid example and
    made the minimal test parse DOPC topology from the ITP instead of duplicating
    bead metadata.
  - Scoped the minimal CG-lipid test to GLY so it does not fit unused sidechain
    tables.
- Files modified:
  - `.gitignore`
  - `example/16.MARTINI/readme.md`
  - `example/16.MARTINI/test_cg_lipid/run_test.py`
  - `parameters/dryMARTINI/*.h5` generated artifacts removed
  - `plan.md`
  - `progress.md`
  - `findings.md`
  - `py/martini_build_tables.py`
  - `py/run_upside.py`
  - `py/upside_martini.py` removed
  - `src/CMakeLists_M1.txt`
  - `src/CMakeLists_Other.txt`
  - `src/bead_interaction.h`
  - `src/extra_particles.h` removed
  - `src/martini*.cpp`
  - `src/martini_internal.h`
  - `src/rotamer.cpp`
  - `src/state_logger.h`
- Verification:
  - `python3 -m py_compile` passed for MARTINI Python modules and examples.
  - `cmake --build obj` passed.
  - `git diff --check` passed.
  - `git diff -- py/martinize.py` is empty.
  - `example/16.MARTINI/test_cg_lipid/run_test.py` passed: generated
    `martini.h5`, converted the stage, injected CG-lipid nodes, and completed
    minimization with final potential `-183.983780`.
  - Full production workflow parity was not rerun because it is expensive; this
    is documented in `plan.md`.
  - Wider refactors from Claude's proposal that need their own parity baseline
    are documented as remaining risk in `plan.md`.
