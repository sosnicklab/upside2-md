# Task Plan

## Project Goal
Audit and fix the dry-MARTINI Upside hybrid integration workflow so that:
- sidechain rotamer numbering is preserved exactly through `example/16.MARTINI/run_sim_1rkl.sh` and `src/martini.cpp`
- hybrid sidechain structure sourcing uses MARTINI 2.2 sidechain structure rather than all-atom structure where the current workflow is inconsistent

## Architecture & Key Decisions
- Treat `parameters/ff_2.1/sidechain.h5` and `py/upside_config.py` as the source of truth for rotamer state numbering.
- Compare script-time node injection and runtime decoding directly rather than inferring from prior notes.
- Make the minimal change necessary to preserve numbering and switch structure sourcing without altering unrelated hybrid behavior.
- Verify with concrete artifacts or targeted checks after editing.

## Execution Phases
- [x] Refresh task tracking files and read existing task context for the MARTINI hybrid workflow
- [x] Inspect `run_sim_1rkl.sh` and `src/martini.cpp` for SC numbering and sidechain structure sourcing
- [x] Implement minimal fixes for exact SC numbering and MARTINI 2.2 sidechain structure usage
- [x] Verify with targeted checks and document findings

## Review
- Confirmed that stage-7 rotamer injection already mirrors Upside's encoded `id_seq` layout; added explicit validation so numbering mismatches fail early.
- Changed SC map generation to follow the MARTINI CG sidechain proxy structure directly instead of gating on all-atom residue lookup.
- Added runtime validation in `src/martini.cpp` so placement-group numbering and SC row assignments must preserve exact 0-based rotamer ids.
- Verified Python syntax, shell syntax, C++ build, and a targeted MARTINI-asset replay for SC map generation and rotamer numbering.

## Known Errors / Blockers
- Need to determine whether numbering drift occurs in stage-file injection, runtime placement-group decoding, or both.
- Need to identify the exact current all-atom sidechain structure dependency before replacing it with MARTINI 2.2 structure.
