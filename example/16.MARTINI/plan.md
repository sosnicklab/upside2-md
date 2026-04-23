# Workflow Update: 1AFO Default Auto-Resume

## Project Goal
- Make `run_sim_1afo.sh` default to continuing from the most recent previous stage-7 checkpoint of a prior `1afo` AABB run.
- If no prior `1afo` stage-7 artifact exists, fall back to the scratch workflow.
- Preserve the current `run_sim_1afo_outlipid.sh` behavior without letting the new base-script auto-resume accidentally switch an outlipid run onto an AABB checkpoint.

## Architecture & Key Decisions
- Keep `run_sim_1afo.sh` as a thin wrapper over `run_sim_1rkl.sh`, but add the same explicit-override and auto-detect continuation logic used by the outlipid wrapper.
- Search only `outputs/martini_test_1afo_aabb*/checkpoints/` for prior AABB stage-7 artifacts.
- Add a dedicated opt-out flag in `run_sim_1afo.sh` so callers can disable AABB auto-resume when needed.
- Have `run_sim_1afo_outlipid.sh` set that opt-out before delegating, so outlipid continuation remains scoped to outlipid artifacts only.

## Execution Phases
- [x] Phase 1: Confirm the current `1afo` wrappers and identify the safe auto-resume integration point.
- [x] Phase 2: Patch `run_sim_1afo.sh` for default auto-resume with scratch fallback, and isolate `run_sim_1afo_outlipid.sh` from that new base behavior.
- [x] Phase 3: Verify shell syntax plus both scratch and continuation behavior in isolated harnesses.

## Known Errors / Blockers
- No blocker.

## Review
- Updated `run_sim_1afo.sh` so it now behaves like the outlipid wrapper with respect to continuation:
  - explicit precedence remains:
    - `CONTINUE_STAGE_70_FROM`
    - `PREVIOUS_STAGE7_FILE`
    - `PREVIOUS_RUN_DIR`
  - otherwise it auto-detects the newest prior AABB stage-7 artifact under `outputs/martini_test_1afo_aabb*/checkpoints/`,
  - if one is found it continues into `outputs/martini_test_1afo_aabb_continue`,
  - if none is found it falls back to the scratch default `outputs/martini_test_1afo_aabb`.
- Added an internal guard so `run_sim_1afo_outlipid.sh` does not inherit the new AABB auto-resume path:
  - it now exports `DISABLE_1AFO_AABB_AUTO_CONTINUE=1` before delegating to `run_sim_1afo.sh`,
  - so outlipid continuation remains driven only by outlipid stage-7 artifacts.
- Verification:
  - `bash -n run_sim_1afo.sh`
  - `bash -n run_sim_1afo_outlipid.sh`
  - isolated AABB harness with no prior stage-7 file:
    - `CONTINUE_STAGE_70_FROM=`
    - `RUN_DIR=outputs/martini_test_1afo_aabb`
  - isolated AABB harness with prior `1afo.stage_7.0.up`:
    - `CONTINUE_STAGE_70_FROM=<detected AABB stage-7 file>`
    - `RUN_DIR=outputs/martini_test_1afo_aabb_continue`
  - isolated outlipid harness with only an AABB prior stage-7 file:
    - `CONTINUE_STAGE_70_FROM=` and `DISABLE_1AFO_AABB_AUTO_CONTINUE=1`
  - isolated outlipid harness with an outlipid prior stage-7 file:
    - `CONTINUE_STAGE_70_FROM=<detected outlipid stage-7 file>`
    - `RUN_DIR=outputs/martini_test_1afo_outlipid_continue`
