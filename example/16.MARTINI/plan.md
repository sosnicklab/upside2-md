# Workflow Addition: 1AFO Scripts

## Project Goal
- Create `run_sim_1afo.sh` and `run_sim_1afo_outlipid.sh` based on the current `1rkl` workflows.
- The new scripts should target the local `1AFO` input structure while keeping the existing workflow behavior and interfaces as intact as possible.
- Keep changes minimal and scoped to the new scripts plus task tracking.

## Architecture & Key Decisions
- Reuse the current `1rkl` workflows through thin wrappers instead of copying the full implementation.
- Keep the user-facing workflow names and output directories in lowercase `1afo`.
- Point the default input structure to the actual repo file `pdb/1AFO.pdb` because that is the file present in this workspace.
- Preserve the current outlipid auto-resume pattern, but specialize it to `1afo` run directories and stage filenames.

## Execution Phases
- [x] Phase 1: Inspect the current `1rkl` workflow scripts and confirm the local `1AFO` input artifact.
- [x] Phase 2: Add `run_sim_1afo.sh` and `run_sim_1afo_outlipid.sh` as minimal wrappers with `1afo` defaults.
- [x] Phase 3: Verify shell syntax and the key environment/argument wiring for both new scripts.

## Known Errors / Blockers
- No blocker.

## Review
- Added `run_sim_1afo.sh` as a thin wrapper over `run_sim_1rkl.sh`:
  - defaults `PDB_ID=1afo`,
  - defaults `RUN_DIR=outputs/martini_test_1afo_aabb`,
  - defaults `RUNTIME_PDB_ID=1afo_aabb`,
  - defaults `PROTEIN_AA_PDB=pdb/1AFO.pdb`.
- Added `run_sim_1afo_outlipid.sh` as a `1afo`-specialized version of the current outlipid wrapper:
  - delegates to `run_sim_1afo.sh`,
  - uses `1afo` output names and checkpoint filenames,
  - auto-detects previous `1afo.stage_7.0.up` / `1afo.stage_7.0.continue.up` files under `outputs/martini_test_1afo_outlipid*/checkpoints/`,
  - falls back to scratch when no previous outlipid stage-7 file exists.
- Verification:
  - `bash -n run_sim_1afo.sh`
  - `bash -n run_sim_1afo_outlipid.sh`
  - isolated wrapper harness for `run_sim_1afo.sh` confirmed:
    - `arg1=PDB_ID=1afo`
    - `RUN_DIR=outputs/martini_test_1afo_aabb`
    - `RUNTIME_PDB_ID=1afo_aabb`
    - `PROTEIN_AA_PDB=pdb/1AFO.pdb`
  - isolated outlipid harness confirmed:
    - with a prior `1afo.stage_7.0.up`, continuation is auto-selected and `RUN_DIR=outputs/martini_test_1afo_outlipid_continue`,
    - without a prior stage-7 file, `CONTINUE_STAGE_70_FROM` stays empty and `RUN_DIR=outputs/martini_test_1afo_outlipid`.
