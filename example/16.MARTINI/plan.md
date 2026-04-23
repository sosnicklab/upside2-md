# Verification: 1AFO Outlipid Unaffected By AABB Slurm Fix

## Project Goal
- Verify that the recent Slurm/bootstrap modification to `run_sim_1afo.sh` does not change the intended behavior of `run_sim_1afo_outlipid.sh`.
- Confirm that the outlipid wrapper still:
  - uses outlipid defaults,
  - resumes only from outlipid stage-7 artifacts,
  - ignores AABB auto-resume,
  - and delegates under Slurm with the correct project-root environment.

## Architecture & Key Decisions
- Treat this as a verification task first; only patch code if the check exposes a regression.
- Validate behavior in isolated harnesses so the real repo outputs are not disturbed.
- Explicitly test both:
  - an AABB prior stage-7 artifact present without an outlipid artifact,
  - an outlipid prior stage-7 artifact present.

## Execution Phases
- [x] Phase 1: Inspect the current `1afo` AABB and outlipid wrappers after the Slurm fix.
- [x] Phase 2: Run isolated Slurm-style harnesses to verify that `run_sim_1afo_outlipid.sh` still uses outlipid-only continuation behavior and correct environment bootstrap.
- [x] Phase 3: Record the verification result in the local trackers and summarize whether the outlipid wrapper is unaffected.

## Known Errors / Blockers
- No blocker.

## Review
- `run_sim_1afo_outlipid.sh` remains unaffected by the `run_sim_1afo.sh` Slurm bootstrap change.
- Verified properties:
  - it still exports `DISABLE_1AFO_AABB_AUTO_CONTINUE=1` before delegation,
  - with only an AABB prior stage-7 artifact present, it does not resume and keeps:
    - `CONTINUE_STAGE_70_FROM=`
    - `RUN_DIR=outputs/martini_test_1afo_outlipid`
    - `RUNTIME_PDB_ID=1afo_outlipid`
  - with an outlipid prior stage-7 artifact present, it still resumes from that outlipid file and uses:
    - `RUN_DIR=outputs/martini_test_1afo_outlipid_continue`
  - under a simulated `sbatch` context with a wrong inherited `UPSIDE_HOME`, it still delegates with:
    - `UPSIDE_HOME=<resolved project root>`
    - `UPSIDE_SKIP_SOURCE_SH=1`
    - `PYTHONPATH=<resolved project root>/py`
    - `PATH` containing `<resolved project root>/obj`
- Verification:
  - `bash -n run_sim_1afo_outlipid.sh`
  - isolated Slurm-style harness with only `outputs/martini_test_1afo_aabb/checkpoints/1afo.stage_7.0.up`
  - isolated Slurm-style harness with `outputs/martini_test_1afo_outlipid/checkpoints/1afo.stage_7.0.up`
