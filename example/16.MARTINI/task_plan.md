# Task Plan

## Project Goal
- Fix `scan_lipid_diffusion.py` so Slurm workflow generation covers the intended tau and temperature sweep (cross-checked against `scan_water_diffusion.py`).

## Architecture & Key Decisions
- Keep the existing script-generation workflow and Slurm array model unchanged.
- Make the tau/T sweep builders in lipid mirror water defaults so behavior is consistent across both scan drivers.
- Keep the change minimal: only sweep-range logic and argument-parsing structure, while preserving bilayer-specific runtime behavior.

## Execution Phases
- [x] Phase 1: Confirm current lipid vs water sweep mismatch.
- [x] Phase 2: Implement lipid sweep/CLI updates.
- [x] Phase 3: Generate scripts and verify folder/task counts and Slurm array bounds.
- [x] Phase 4: Record outcomes and any residual blockers.

## Known Errors / Blockers
- No active blockers.

## Revised Decisions
- Removed the intermediate broad CLI tau/T parameterization and kept a minimal parity change with `scan_water_diffusion.py` as requested.
