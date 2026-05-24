## 2026-05-24 MARTINI Python Runner

### Project Goal
- Update `example/16.MARTINI/run.py` so it follows the compact style of the
  other example `run.py` scripts.
- Make it runnable from a local Mac shell or from Slurm.
- Accept an OPM-style input PDB and drive the existing MARTINI hybrid workflow.
- Add a lipid representation toggle for `coarse` CG-lipid and `full`
  single-particle/full bead mode.
- Expose lipid composition configuration, currently restricted to 100% DOPC.

### Architecture & Key Decisions
- Keep `run.py` as a thin orchestrator. Reuse `py/martini_prepare_system.py
  run-hybrid-workflow` and do not duplicate stage-generation or simulation
  logic.
- Use environment variables only for the lower-level workflow settings already
  consumed by the MARTINI scripts.
- For Slurm mode, write a small self-contained submission script that sets the
  project environment explicitly and sets `UPSIDE_SKIP_SOURCE_SH=1`.
- Do not alter MARTINI force-field physics or interaction tables.

- [x] Phase 1: Inspect existing `run.py`, example runner style, and workflow
      arguments.
- [x] Phase 2: Rewrite `example/16.MARTINI/run.py` as the shared runner.
- [x] Phase 3: Verify syntax/help and document behavior.

### Known Errors / Blockers
- None.

### Review
- Rewrote `example/16.MARTINI/run.py` in the editable style used by the other
  example runners.
- The runner now accepts an OPM-oriented protein PDB path, validates the current
  100% DOPC lipid composition setting, and maps `lipid_model` to the existing
  workflow's `--lipid-resolution`.
- Local mode sets `UPSIDE_HOME`, `PATH`, and `PYTHONPATH` and generates missing
  MARTINI parameter files before running the workflow.
- Slurm mode writes and submits a self-contained batch script with module loads,
  repo `.venv` activation, explicit environment exports, and
  `UPSIDE_SKIP_SOURCE_SH=1`.
- Verification passed: `python3 -m py_compile`, style/ASCII scan, and
  `git diff --check`.
