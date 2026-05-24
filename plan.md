## 2026-05-24 MARTINI Workflow Cleanup

### Project Goal
- Straighten the current MARTINI C++/Python workflow logic for the active files:
  `src/martini*`, `src/box*`, branch diffs in `src/main.*`, `py/martini_*.py`,
  and the 1RKL/1AFO example launchers.
- Make 1RKL and 1AFO launchers share the same workflow semantics within each
  lipid-resolution class: stage order, defaults, restart rules, output layout,
  and production log format.
- Remove debug-only user-facing controls and debug progress output.
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

### Execution Phases
- [x] Phase 1: Update task notes and keep progress concise.
- [x] Phase 2: Add shared launcher and convert 1RKL/1AFO scripts to wrappers.
- [x] Phase 3: Remove debug-only Python workflow paths and duplicated stage-7
      orchestration.
- [x] Phase 4: Clean C++ public declarations and progress logging.
- [x] Phase 5: Run syntax/build/log-format verification and document results.

### Known Errors / Blockers
- None.

### Review
- Shared launch semantics now live in `example/16.MARTINI/run_sim_hybrid.sh`;
  1RKL/1AFO coarse/full scripts only set system-specific defaults.
- Python workflow no longer exposes or invokes debug-only initial-structure or
  stage-debug output paths. Stage 7 production handoff uses one helper path for
  direct and extended equilibration flows.
- MARTINI C++ declarations moved to `src/martini.h`; `src/main.h` is back to the
  public `upside_main` entry point. The hybrid progress line reports
  `protein_potential` and `total_potential` only.
- Verification passed: Python compile, shell syntax, debug-log scan, and
  `cmake --build obj`.
