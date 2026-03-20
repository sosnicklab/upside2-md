# Task Plan

## Project Goal
Make the hybrid stage-7 workflow leave the protein backbone non-rigid again by default in both the main runner and the production-only test runner.

## Architecture & Key Decisions
- Keep this as a workflow-default change unless inspection shows the stage-7 scripts no longer route through the existing rigid-mask helper.
- Change both runner defaults together so `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` stay aligned.
- Preserve the existing `set_production_backbone_fix_rigid()` helper and stage-7 wiring; only the default enable flag should flip back off.

## Execution Phases
- [x] Re-scope task tracking for the stage-7 non-rigid protein-backbone workflow change
- [x] Update the hybrid runner defaults so production stage 7 no longer enables the backbone rigid mask by default
- [x] Verify the edited workflow scripts and document the result

## Known Errors / Blockers
- None currently.

## Review
- Restored the stage-7 workflow to leave the protein backbone non-rigid by default in:
  - `example/16.MARTINI/run_sim_1rkl.sh`
  - `example/16.MARTINI/test_prod_run_sim_1rkl.sh`
- Kept the existing `set_production_backbone_fix_rigid()` helper and stage-7 conditional wiring unchanged; only the default flag was flipped back to `0`.
- Verified both edited scripts with `bash -n`.
