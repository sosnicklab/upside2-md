# Workflow Update: Default Outlipid Auto-Resume

## Project Goal
- Update `run_sim_1rkl_outlipid.sh` so it resumes from the last available stage-7 frame of a previous outlipid run by default.
- If no previous run artifact is found, fall back to the normal scratch workflow.
- Preserve the existing explicit continuation overrides and keep the change scoped to the wrapper script.

## Architecture & Key Decisions
- Implement auto-resume only in `run_sim_1rkl_outlipid.sh`; do not change the shared base workflow unless a concrete incompatibility appears.
- Keep the current explicit controls first:
  - `CONTINUE_STAGE_70_FROM`
  - `PREVIOUS_STAGE7_FILE`
  - `PREVIOUS_RUN_DIR`
- Only when none of those are set, auto-detect the most recent prior outlipid production checkpoint under `example/16.MARTINI/outputs/`.
- Search only outlipid-specific output trees to avoid resuming from unrelated embedded or hybrid runs.

## Execution Phases
- [x] Phase 1: Inspect the current wrapper continuation logic and the local output naming/layout used for previous outlipid runs.
- [x] Phase 2: Patch the wrapper to auto-detect the latest previous stage-7 artifact, with scratch fallback when none exists.
- [x] Phase 3: Verify shell syntax and exercise the detection result against the current local `outputs/` tree.

## Known Errors / Blockers
- No blocker.

## Review
- Updated `run_sim_1rkl_outlipid.sh` to derive the real workflow directory from `BASE_WORKFLOW_SCRIPT`, then auto-detect the newest prior outlipid stage-7 file under `outputs/` when no explicit continuation source is provided.
- Preserved the explicit continuation precedence:
  - `CONTINUE_STAGE_70_FROM`
  - `PREVIOUS_STAGE7_FILE`
  - `PREVIOUS_RUN_DIR`
- Added an escape hatch for intentional scratch starts:
  - `AUTO_CONTINUE_FROM_PREVIOUS_RUN=0`
- Verification:
  - `bash -n run_sim_1rkl_outlipid.sh`
  - isolated `/tmp` harness with no prior outlipid stage-7 file:
    - `CONTINUE_STAGE_70_FROM=` and `RUN_DIR=outputs/martini_test_1rkl_outlipid`
  - isolated `/tmp` harness with a prior outlipid `checkpoints/1rkl.stage_7.0.up`:
    - `CONTINUE_STAGE_70_FROM=<detected stage-7 file>`
    - `RUN_DIR=outputs/martini_test_1rkl_outlipid_continue`
  - current local repo outputs:
    - no outlipid stage-7 artifact exists yet, so the wrapper will currently fall back to scratch in this workspace.
