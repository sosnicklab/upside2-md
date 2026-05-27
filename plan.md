# 2026-05-26 MARTINI 1AFO/1RKL Stage-7 Diagnostics

## Project Goal
- Diagnose the bent helix in `martini_1afo_hybrid_full/1afo.stage_7.0.vtf`,
  especially whether the production-stage initial protein equals the input PDB
  geometry or has already been changed before production output.
- Diagnose abnormal lipid orientations/z coordinates in
  `martini_1afo_hybrid/1afo.stage_7.0.vtf`.
- Diagnose secondary-structure loss in both `martini_1rkl_hybrid` and
  `martini_1rkl_hybrid_full`.
- Fix confirmed bugs physically; do not disable or scale away SC-env, BB-env,
  lipid-lipid, CGL-SC, or CGL-target interactions.

## Architecture & Key Decisions
- Inspect actual HDF5 checkpoints and VTF files under `example/16.MARTINI/outputs`.
  Use copied files for any command that might rewrite HDF5.
- Separate raw stage handoff coordinates, production-frame coordinates, and VTF
  display/extraction artifacts before changing force-field code.
- Compare full and CGL lipid modes against the same protein reference and stage
  boundaries. A valid fix must preserve charged termini for every chain
  fragment and preserve additive backbone/sidechain force ownership.

## Execution Phases
- [x] Phase 1: Inventory stage-7 HDF5/VTF/log artifacts for 1AFO and 1RKL.
- [x] Phase 2: Quantify protein geometry at raw input, minimized/burn-in
      handoff, and VTF frame 0 where available.
- [x] Phase 3: Quantify lipid orientation vectors and z distributions in the
      1AFO coarse output and distinguish runtime versus VTF display issues.
- [x] Phase 4: Identify confirmed code or force-field bug and patch only that.
- [x] Phase 5: Run focused verification and document remaining limits.

## Known Errors / Blockers
- Existing trajectories under `example/16.MARTINI/outputs` remain stale. The
  fixed SC-particle behavior requires regenerating `parameters/dryMARTINI/sidechain.h5`
  and then rerunning the affected workflows.

## Review
- Prepared stage-7 protein coordinates match the reference/PDB geometry after
  rigid alignment (`~2-3e-6 A` CA RMSD). Production frame 0 is exactly the
  promoted stage-7 input, so visible bends are burn-in results, not import or
  VTF extraction artifacts.
- Full-resolution lipid damage is explained by the old SC-particle table
  factorization, which created a runtime minimum around `-4.7e11 kJ/mol`.
  A temporary rebuilt table with the patched factorization has minimum
  `-17.3 kJ/mol`.
- CGL lipid tilt is a runtime CGLD orientation drift. A copied 1AFO CGL
  burn-in with the derived leaflet-normal term kept all 182 CGL vectors upright
  over 4000 steps (`min |n_z|=0.765`, zero below `0.5`).
- Verified: Python compile, C++ build, `git diff --check`, termini charges in
  all four current 1AFO/1RKL stage files.
