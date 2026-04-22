# Workflow-16 Local Bilayer Asset Update

## Project Goal
- Move `DOPC.pdb` from `parameters/dryMARTINI/` into this example's local `pdb/` directory.
- Update `run_sim_1rkl.sh` and `run_sim_1rkl_outlipid.sh` so workflow-16 uses the example-local bilayer file by default.
- Update `readme.md` to state that the bilayer structure should be generated from CHARMM-GUI.

## Architecture & Key Decisions
- Keep the change local to `example/16.MARTINI/`; do not change shared prep code or unrelated examples.
- Make `pdb/DOPC.pdb` the default bilayer path so the example is self-contained.
- Update both shell entry points explicitly, even though the wrapper inherits the base workflow, so the local bilayer default is obvious at both user-facing launch surfaces.

## Execution Phases
- [x] Phase 1: Re-read `plan.md`, `findings.md`, `progress.md`, and inspect current bilayer-path usage in the workflow scripts and README.
- [x] Phase 2: Patch the local trackers, scripts, and README for the example-local bilayer asset and CHARMM-GUI documentation.
- [x] Phase 3: Move `DOPC.pdb` into `example/16.MARTINI/pdb/` and verify the updated references.

## Known Errors / Blockers
- Moving the source file out of `parameters/dryMARTINI/` requires a write outside the example directory, so the file move may need escalated permission.

## Review
- Moved the default bilayer asset into `example/16.MARTINI/pdb/DOPC.pdb`.
- Updated both user-facing shell entry points to default `BILAYER_PDB` to `pdb/DOPC.pdb`.
- Updated `readme.md` to state that the bilayer structure should be generated from CHARMM-GUI and stored locally in the example.
- Verification:
  - `ls -la pdb/DOPC.pdb`
  - confirmed `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/DOPC.pdb` no longer exists after the move
  - `bash -n run_sim_1rkl.sh`
  - `bash -n run_sim_1rkl_outlipid.sh`
