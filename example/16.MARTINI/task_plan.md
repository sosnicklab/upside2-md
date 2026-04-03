# Task Plan

## Project Goal
- Fix the MARTINI workflow mass force-field path handling so an absolute `MASS_FF_FILE` does not get prefixed with `example/16.MARTINI`.

## Architecture & Key Decisions
- Root cause is in `prepare_protein_inputs()` inside `run_sim_1rkl.sh`:
  - `MASS_FF_FILE` now defaults to an absolute path under `${UPSIDE_HOME}/parameters/dryMARTINI`,
  - but the script was still constructing `mass_ff_path` as `"${SCRIPT_DIR}/${MASS_FF_FILE}"`,
  - producing an invalid doubled path for absolute inputs.
- Fix only the shell path resolution:
  - resolve `MASS_FF_FILE` through `abspath(expanduser(...))`,
  - do not change the MARTINI workflow contract otherwise.

## Execution Phases
- [x] Phase 1: Audit `MASS_FF_FILE` path handling in the shell workflow.
- [x] Phase 2: Fix absolute-vs-relative mass force-field path resolution.
- [x] Phase 3: Verify the updated script behavior and syntax.

## Known Errors / Blockers
- No blocker identified.

## Review
- Fixed [run_sim_1rkl.sh](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/run_sim_1rkl.sh) so `prepare_protein_inputs()` resolves `MASS_FF_FILE` safely whether it is absolute or relative.
- Root cause was the old `"${SCRIPT_DIR}/${MASS_FF_FILE}"` join against an already-absolute default path.
- Verification passed:
  - `python3 -c 'import os; p=\"/Users/yinhan/Documents/upside2-md-martini/parameters/dryMARTINI/dry_martini_v2.1.itp\"; print(os.path.abspath(os.path.expanduser(p)))'`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
