# 2026-05-24 MARTINI Branch Cleanup

## Project Goal
- Rewrite the MARTINI additions on `martini-dev` so they are concise, closer to
  master style, and easier to maintain.
- Preserve the current physical model and keep `example/16.MARTINI` workflows at
  metric parity with the current branch.
- Treat generated MARTINI parameter files as rebuildable artifacts.

## Architecture & Key Decisions
- Scope is MARTINI-only: `src/box.*`, `src/martini*`, shared C++ files touched
  for MARTINI support, `py/martini_*`, and `example/16.MARTINI`.
- `py/martinize.py` is vendored third-party code. Do not edit, format, delete,
  annotate, or otherwise touch it.
- Restore master parity in shared files unless the MARTINI workflow needs the
  change.
- Keep all SC-env, BB-env, CGL-CGL, CGL-SC, and CGL-target interactions active.
- Use metric parity, not byte-identical trajectories, as workflow acceptance.

## Execution Phases
- [x] Phase 0: Capture baseline checks and short-run metrics where practical.
- [x] Phase 1: Remove unused placeholder code and restore accidental shared
      changes.
- [x] Phase 2: Split `src/martini.cpp` by subsystem and update build files.
- [x] Phase 3: Deduplicate C++ integration/helpers without changing runtime
      semantics.
- [x] Phase 4: Clean MARTINI Python modules and example wrappers.
- [x] Phase 5: Run build/tests and record parity results.

## Known Errors / Blockers
- Full 1RKL/1AFO production workflow parity was not rerun in this cleanup pass;
  the focused CG-lipid workflow was used as the executable MARTINI check.
- Larger structural refactors from Claude's proposal, such as splitting
  `py/martini_prepare_system_lib.py` and refactoring shared integrator internals,
  should only be done with a dedicated parity baseline because they carry a
  wider workflow risk than the verified cleanup pass here.

## Review
- `src/martini.cpp` was split into subsystem files with shared private
  declarations in `src/martini_internal.h`; build files were updated.
- Restored shared files toward master style where MARTINI did not need the
  change, removed unused MARTINI placeholders, and left `py/martinize.py`
  untouched.
- Deleted generated dry-MARTINI HDF5 artifacts from version control scope and
  ignored regenerated copies.
- Verification passed:
  - `python3 -m py_compile` for MARTINI Python modules and example scripts.
  - `cmake --build obj`.
  - `git diff --check`.
  - `example/16.MARTINI/test_cg_lipid/run_test.py` completed table build,
    stage conversion, table injection, and minimization.
