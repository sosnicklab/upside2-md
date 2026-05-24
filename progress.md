# Progress Log

## 2026-05-24 MARTINI Python Runner
- Actions taken:
  - Replaced `example/16.MARTINI/run.py` with a compact editable example runner.
  - Added explicit OPM PDB input, DOPC-only lipid composition validation, and a
    `lipid_model` toggle that maps to the existing workflow's coarse/full lipid
    resolution.
  - Added local Mac execution setup using the repo `.venv` when available.
  - Added Slurm submission path that writes a self-contained batch script with
    module setup, `.venv` activation, explicit `UPSIDE_HOME/PATH/PYTHONPATH`,
    and `UPSIDE_SKIP_SOURCE_SH=1`.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `example/16.MARTINI/run.py`
- Verification:
  - `python3 -m py_compile example/16.MARTINI/run.py` passed.
  - ASCII/style scan found no debug markers or non-ASCII text in `run.py`.
  - `git diff --check` passed for the touched files.
  - Full execution was not run because it would launch the MARTINI simulation.
