# Task Plan

## Project Goal
Change the dry-MARTINI Upside hybrid workflow in `example/16.MARTINI/run_sim_1rkl.sh` so that:
- SC-env, SC-BB, and env-BB LJ/Coulomb force-cap removal ramps over the first `100` production steps
- bilayer `PO4` z-coordinate holding remains active for the first `150` production steps
- the existing `200`-step hybrid interaction startup window otherwise stays unchanged

## Architecture & Key Decisions
- Limit the change to workflow defaults in `example/16.MARTINI/run_sim_1rkl.sh` unless inspection shows the runtime cannot already represent separate windows.
- Preserve the current hybrid-control attribute interface (`sc_env_relax_steps`, `sc_env_po4_z_hold_steps`) instead of introducing new controls.
- Verify the runtime path consumes those attrs independently before deciding whether any C++ edit is necessary.
- Keep the change minimal: update only the requested schedule and the comments that describe it.

## Execution Phases
- [x] Inspect the workflow script and runtime handling for SC relax steps and PO4 z-hold steps
- [x] Refresh task tracking files for this schedule-change task
- [x] Patch `example/16.MARTINI/run_sim_1rkl.sh` to use `100` relax steps and `150` PO4 z-hold steps
- [x] Verify the edited workflow by checking the script text and shell syntax

## Known Errors / Blockers
- None currently. The runtime already appears to parse and enforce separate `sc_env_relax_steps` and `sc_env_po4_z_hold_steps` values.

## Review
- Updated `example/16.MARTINI/run_sim_1rkl.sh` so `SC_ENV_RELAX_STEPS` now defaults to `100` while `SC_ENV_PO4_Z_HOLD_STEPS` is explicitly kept at `150`.
- Kept `SC_ENV_BACKBONE_HOLD_STEPS` at `200`, preserving the existing startup interaction window outside the requested ramp change.
- Verified the edited control block and ran `bash -n example/16.MARTINI/run_sim_1rkl.sh` successfully.
