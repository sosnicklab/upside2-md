# Task Plan

## Project Goal
- Reduce Python files under `example/16.MARTINI` to exactly four:
  - `martinize.py`
  - one simulation-input generation script
  - one VTF conversion script
  - one helper module
- Keep `run_sim_1rkl.sh` working by routing all stage-preparation, validation, SC-table build, stage-7 injection, and handoff utilities through the single input-generation script plus the helper module.

## Architecture & Key Decisions
- Keep `martinize.py` unchanged as the explicit MARTINI coarse-graining script.
- Keep `extract_martini_vtf.py` as the explicit trajectory/VTF export script.
- Use `prepare_system.py` as the only other workflow-facing Python entry point.
- Merge the logic from `build_sc_martini_h5.py`, `inject_sc_table_stage7.py`, `set_initial_position.py`, and `validate_hybrid_mapping.py` into `prepare_system.py` plus `prepare_system_lib.py`.
- Update `run_sim_1rkl.sh` to call `prepare_system.py` subcommands instead of separate Python entry points.
- Remove obsolete Python scripts after the merged path is verified.

## Execution Phases
- [x] Phase 1: Recreate task tracking and record the corrected four-file requirement.
- [x] Phase 2: Merge extra Python entry-point logic into `prepare_system_lib.py` and expose it through `prepare_system.py`.
- [x] Phase 3: Rewire `run_sim_1rkl.sh` to the consolidated Python interface.
- [x] Phase 4: Delete obsolete Python scripts and caches, then verify the final tree and workflow parsing.

## Known Errors / Blockers
- No remaining blockers for this cleanup task.

## Review
- Final retained Python files under `example/16.MARTINI`:
  - `martinize.py`
  - `prepare_system.py`
  - `extract_martini_vtf.py`
  - `prepare_system_lib.py`
- Removed obsolete Python entry points:
  - `build_sc_martini_h5.py`
  - `inject_sc_table_stage7.py`
  - `set_initial_position.py`
  - `validate_hybrid_mapping.py`
- `prepare_system.py` now owns the extra workflow entry points through subcommands:
  - `build-sc-martini-h5`
  - `inject-stage7-sc`
  - `set-initial-position`
  - `validate-hybrid-mapping`
- Verification passed:
  - `find example/16.MARTINI -maxdepth 1 -name '*.py' | sort`
  - `rg -n "build_sc_martini_h5\.py|inject_sc_table_stage7\.py|set_initial_position\.py|validate_hybrid_mapping\.py" example/16.MARTINI/run_sim_1rkl.sh example/16.MARTINI/*.py`
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile example/16.MARTINI/prepare_system.py example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/extract_martini_vtf.py example/16.MARTINI/martinize.py`
  - `source .venv/bin/activate && source source.sh && bash -n example/16.MARTINI/run_sim_1rkl.sh`
