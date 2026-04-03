# Task Plan

## Project Goal
- Add example-style entrypoint files for `example/16.MARTINI`:
  - `readme.md`
  - `run.py`
- Match the lightweight style of existing examples such as `example/02.ReplicaExchangeSimulation`, while keeping the real workflow implemented in `run_sim_1rkl.sh`.

## Architecture & Key Decisions
- Keep `run_sim_1rkl.sh` as the actual MARTINI workflow implementation.
- Add `run.py` as a thin Python wrapper that:
  - exposes the main workflow knobs near the top of the file,
  - passes them to `run_sim_1rkl.sh` through environment variables,
  - invokes the shell workflow with `subprocess`.
- Add `readme.md` in the same simple style as the other example directories:
  - short description,
  - minimal run commands,
  - short note about outputs and visualization.
- This task intentionally adds `run.py` even though the directory was previously reduced to four Python files, because the user explicitly requested a standard example-style Python entrypoint.

## Execution Phases
- [x] Phase 1: Inspect existing example `readme.md` / `run.py` style and current MARTINI workflow entrypoints.
- [x] Phase 2: Implement `example/16.MARTINI/run.py` as a wrapper around `run_sim_1rkl.sh`.
- [x] Phase 3: Implement `example/16.MARTINI/readme.md` with matching example documentation style.
- [x] Phase 4: Verify the new files parse cleanly and document the result.

## Known Errors / Blockers
- No remaining blocker for this task.

## Review
- Added [run.py](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/run.py) as a thin wrapper around [run_sim_1rkl.sh](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/run_sim_1rkl.sh):
  - exposes the main workflow settings at the top of the file,
  - passes them as environment overrides,
  - launches the shell workflow with `subprocess`.
- Added [readme.md](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/readme.md) in the lightweight style used by the other example directories.
- Verification passed:
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile example/16.MARTINI/run.py example/16.MARTINI/prepare_system.py example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/extract_martini_vtf.py example/16.MARTINI/martinize.py`
  - `source .venv/bin/activate && source source.sh && bash -n example/16.MARTINI/run_sim_1rkl.sh`
