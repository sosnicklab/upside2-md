# Task: Enforce exact 6.6->7.0 handoff without production recenter

## 1. Project Goal
Follow user-required semantics: `6.6 -> 7.0` should copy coordinates exactly (no recenter), with stage changes limited to interaction table/stage settings (NPT->NVT behavior) for `example/16.MARTINI/outputs/martini_test_1rkl_rigid_dry`.

## 2. Architecture & Key Decisions
- Scope:
  - Remove production recenter mode usage from rigid-dry workflow script.
  - Reapply strict handoff for current `6.6 -> 7.0` run output.
  - Validate transition exactness and initial 7.0 potential sanity.
- Key decisions:
  - Use strict copy for all stage handoffs in `run_sim_1rkl_rigid_dry.sh`.
  - Keep `set_initial_position.py` granular toggles available but do not use production recenter in this workflow.
- Revised Decisions:
  - None yet.

## 3. Execution Phases
- [x] Phase 1: Confirm user-required handoff semantics (no recenter for `6.6 -> 7.0`).
- [x] Phase 2: Remove production-center handoff mode from rigid-dry workflow.
- [x] Phase 3: Reapply strict `6.6 -> 7.0` handoff in current output set.
- [x] Phase 4: Validate exact transition and initial `7.0` potential sanity.
- [x] Phase 5: Validate syntax and log outcomes in `progress.md` + `lessons.md`.

## 4. Known Errors / Blockers
- None currently.

## 5. Follow-Up Task (2026-02-23): Prepare-System-Only Preflight

### Goal
- Keep `prepare_system.py` as the only workflow prep entrypoint and add strict AA-backbone (`N/CA/C/O`) vs MARTINI-`BB` coordinate-frame preflight checks before hybrid mapping export.

### Execution
- [x] Confirm `run_sim_1rkl.sh` and `run_sim_1rkl_rigid_dry.sh` prep entrypoint is `prepare_system.py`.
- [x] Remove runtime dependency on deleted `prepare_hybrid_system.py`.
- [x] Add reusable backbone-frame diagnostics + validation (`matched residues`, rigid-fit RMSD, translation-only RMSD, centroid shift, rotation angle).
- [x] Enforce preflight in `prepare_system.py --mode both` mapping path with fail-fast thresholds.
- [x] Validate successful prep on `1rkl` assets and fail-fast behavior on mismatched AA/CG inputs.

### Review
- `prepare_system.py` now imports shared helpers from `prepare_system_lib.py` (non-entrypoint utility module), so workflows continue to execute only through `prepare_system.py`.
- Preflight passes for `1rkl.AA.pdb -> 1rkl_hybrid.MARTINI.pdb` and records diagnostics in summary JSON.
- Preflight correctly fails on inconsistent AA/CG inputs with explicit error diagnostics.

## 6. Follow-Up Task (2026-02-23): Port Rigid-Dry Handoff/Prep Fixes to run_sim_1rkl.sh

### Goal
- Apply the same validated stage handoff behavior from `run_sim_1rkl_rigid_dry.sh` to `run_sim_1rkl.sh` for transitions including `6.0->6.1` and `6.6->7.0`, while keeping hybrid-production behavior unchanged.

### Execution
- [x] Compare stage handoff implementation between `run_sim_1rkl.sh` and `run_sim_1rkl_rigid_dry.sh`.
- [x] Add strict handoff control (`STRICT_STAGE_HANDOFF`) and `handoff_initial_position()` wrapper in `run_sim_1rkl.sh`.
- [x] Replace direct `set_initial_position.py` transition calls in `run_sim_1rkl.sh` with the wrapper.
- [x] Confirm preparation path remains `prepare_system.py` (no `prepare_hybrid_system.py` usage).
- [x] Validate script syntax and preparation smoke.

### Review
- `run_sim_1rkl.sh` now uses strict-copy-capable handoff wrapper for all stage transitions, including `6.0->6.1` and `6.6->7.0`.
- Stage-0 preparation path for `run_sim_1rkl.sh` remains universal prep via `prepare_system.py`, and prep smoke confirms AA/BB preflight diagnostics are active.
- `run_sim_1rkl.sh` now exposes and forwards AA/BB preflight thresholds (`BB_AA_MIN_MATCHED_RESIDUES`, `BB_AA_MAX_RIGID_RMSD`) into stage-0 `prepare_system.py`.

## 7. Follow-Up Task (2026-02-23): Single Conversion Entrypoint Cleanup

### Goal
- Enforce a single conversion script entrypoint (`prepare_system.py`) across workflows, while keeping conversion helper logic in library modules.

### Execution
- [x] Remove workflow-level calls to `prepare_martini.py` and switch them to `prepare_system.py`.
- [x] Refactor `prepare_system.py` conversion path to call `prepare_martini` as an imported library function (no script subprocess hop).
- [x] Keep `prepare_system_lib.py` as utility-only module and retain `prepare_system.py` as workflow entrypoint.
- [x] Validate updated workflows and stage conversion smoke paths.

### Review
- Workflows now call only `prepare_system.py` for stage conversion (including bilayer and production-only helper workflows).
- `prepare_system.py` no longer shells out to `prepare_martini.py`; it imports and invokes conversion in-process.
- `prepare_martini.py` is now library-only at runtime (`__main__` exits with a pointer to `prepare_system.py`), enforcing a single conversion script entrypoint.
- Conversion behavior validated for both `bilayer` and `both` stage-generation smoke tests.

## 8. Follow-Up Task (2026-02-23): Remove Legacy prepare_martini.py

### Goal
- Delete `prepare_martini.py` entirely after proving it is no longer used, while preserving conversion behavior through `prepare_system.py`.

### Execution
- [x] Confirm remaining code-level dependency graph for `prepare_martini.py`.
- [x] Move conversion implementation into `prepare_system_convert_lib.py`.
- [x] Rewire `prepare_system.py` to import conversion from `prepare_system_convert_lib.py`.
- [x] Delete `prepare_martini.py`.
- [x] Re-validate conversion paths and workflow syntax checks.

### Review
- Legacy conversion script removed.
- Single operational conversion entrypoint remains `prepare_system.py`.
- Conversion helper split is now:
  - `prepare_system_lib.py` (all packing/mapping + stage conversion internals)

## 9. Follow-Up Task (2026-02-23): Collapse to One Lib File

### Goal
- Enforce exactly one entry script and one lib file for preparation/conversion.

### Execution
- [x] Merge stage-conversion helper code into `prepare_system_lib.py`.
- [x] Rewire `prepare_system.py` to use `prepare_system_lib.convert_stage(...)`.
- [x] Delete `prepare_system_convert_lib.py`.
- [x] Re-validate workflows and stage conversion smoke tests.

### Review
- Final file layout is now:
  - `prepare_system.py` (entry script)
  - `prepare_system_lib.py` (single helper lib)
