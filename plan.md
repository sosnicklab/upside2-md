## 2026-05-24 MARTINI Branch Cleanup

### Project Goal
- Straighten the current MARTINI C++/Python workflow logic for the active files:
  `src/martini*`, `src/box*`, branch diffs in `src/main.*`, `py/martini_*.py`,
  and the 1RKL/1AFO example launchers.
- Make 1RKL and 1AFO launchers share the same workflow semantics within each
  lipid-resolution class: stage order, defaults, restart rules, output layout,
  and production log format.
- Remove debug-only user-facing controls and debug progress output.
- Normalize all branch-only MARTINI Python/C++ code against the master checkout
  at `/Users/yinhan/Documents/upside2-md-master`: style, naming, ownership
  boundaries, and duplication.
- Keep all physical MARTINI/Upside interface interactions active; do not disable,
  scale, or retune SC-env, BB-env, CGL-CGL, CGL-SC, CGL-target, CGL-ion excluded
  volume, or leaflet-orientation physics.

### Architecture & Key Decisions
- Use one neutral shared launcher for the MARTINI hybrid workflow and keep the
  named 1RKL/1AFO scripts as thin compatibility wrappers.
- Treat schema names and version attributes as compatibility guards, not debug
  marks; preserve them unless a real schema migration is needed.
- Delete debug-only controls from production workflow code instead of merely
  hiding them: no `DEBUG_RIGID_PROTEIN`, `INITIAL_DEBUG_ONLY`,
  `UPSIDE_WRITE_DEBUG_PDB`, or C++ potential-component dump flag.
- Progress logs should report only `protein_potential` and `total_potential`
  for hybrid MARTINI runs. CGL component potentials remain active internally
  and continue contributing to `DerivEngine::potential`.
- Keep MARTINI minimization and `--duration-steps` support because the workflow
  relies on them for production handoff and exact step-count semantics.
- Keep generated outputs unstaged and avoid any git state-changing commands.
- Prefer removing branch-only wrapper layers and duplicated readers/calculators
  when an existing standard-library or project helper path is available.

### Execution Phases
- [x] Phase 1: Remove debug-only workflow controls and output.
- [x] Phase 2: Audit branch-vs-master diff and identify duplicate/new-code
      cleanup targets.
- [x] Phase 3: Refactor MARTINI Python scripts for style and logic.
- [x] Phase 4: Refactor MARTINI C++ sources and main/box deltas for style and
      ownership.
- [x] Phase 5: Re-run syntax/build/log-format verification and document review.

### Known Errors / Blockers
- None.

### Review
- Moved the generic MARTINI topology parser from stage conversion into
  `py/martini_itp_reader.py`, and removed the duplicated parser from
  `py/martini_prepare_system_lib.py`.
- Removed section banners, debug wording, verbose generated-looking diagnostics,
  and non-ASCII progress text from the MARTINI Python/C++ files touched in this
  cleanup.
- Restored `src/main.h` to master parity, kept master-style OpenMP scheduling in
  `src/main.cpp`, and preserved the required hybrid progress line:
  `protein_potential ... total_potential ...` without CGL component debug terms.
- Verification passed: Python compile, shell syntax, workflow help, clean-shell
  launcher bootstrap check, C++ build, `git diff --check`, ASCII/style scan, and
  removed-debug-control scan. Remaining C++ build warnings are pre-existing.
