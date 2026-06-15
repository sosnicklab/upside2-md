# Project Goal

Ensure the 1RKL and 1AFO CGL hybrid workflows use the same CGL cutoff rules and verify that the CGL lipid model is faster than the full-resolution lipid model for both systems.

# Architecture & Key Decisions

- Do not disable, scale, or bypass SC-env, BB-env, SC-CGL, CGL-target, or CGL-CGL interactions.
- CGL cutoff metadata should be model/table driven, not protein/workflow specific.
- 1RKL and 1AFO may differ in protein size and particle counts, but the CGL-CGL and CGL-target cutoff parameters must come from the same `dopc.h5` metadata and runtime defaults.
- Performance acceptance is system-local: compare 1RKL CGL against 1RKL full lipid, and 1AFO CGL against 1AFO full lipid.

# Execution Phases

- [x] Phase 1: Compare CGL cutoff metadata and runtime attrs across active 1RKL and 1AFO coarse outputs.
- [x] Phase 2: Compare active coarse/full stage-7 timings for 1RKL and 1AFO.
- [x] Phase 3: Fix any workflow-specific cutoff divergence.
- [x] Phase 4: Rebuild/regenerate only if needed, then verify both systems.
- [x] Phase 5: Update project notes with the cutoff and performance result.

# Known Errors / Blockers

- None currently.

# Review

- Active 1RKL and 1AFO coarse outputs use the same CGL cutoff metadata:
  CGL-CGL `cutoff_ang=41.3`, CGL-target `cutoff_ang=26.6`, SC-CGL `cutoff_ang=26.6`,
  bead nonbonded cutoff `1.2 nm`, DOPC axis radius `14.55848 A`, and DOPC perpendicular radius `4.114344 A`.
- The 1RKL and 1AFO wrapper scripts do not set separate CGL cutoffs; both delegate to `run_sim_hybrid.sh`.
- Current active stage-7 timing comparisons pass the performance requirement:
  1RKL CGL `1707.99 us/systems/step` versus full lipid `3328.29 us/systems/step`;
  1AFO CGL `1923.47 us/systems/step` versus full lipid `2522.99 us/systems/step`.
- Current active CGL orientation checks are clean:
  1RKL aligned-z min/p05/mean `0.770/0.849/0.949`, no bad-parallel or flipped rods;
  1AFO aligned-z min/p05/mean `0.832/0.864/0.961`, no bad-parallel or flipped rods.
- No code change was needed because there was no cutoff divergence to fix.
