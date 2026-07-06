# Project Goal

Fix the dryMARTINI hybrid interface so the coarse model resolves the `1afo`
collapse without introducing the coarse-only bent-protein artifact seen in
`1rkl`, and while preserving the `1rkl` helical secondary structure across the
bilayer leaflet interface with one generic physical force field for both
systems.

# Architecture & Key Decisions

- Keep the physical model intact.
  Do not disable SC-env or BB-env interactions, do not add restraints, and do
  not introduce any orientational twisting parameter.
- Treat the current failure as a CGL compaction-state force-field defect.
  The coarse and full-resolution `1afo` workflows start from the same packed
  structure, but only the coarse stage-7 run collapses, so the root cause is
  in the coarse compaction-state mapping and correction tables rather than in
  static lipid occupancy.
- Repair the reference-state contract at the H5 layer.
  The current endpoint-remapped contract widened the physical compact/extended
  centers enough to suppress the compression-side response that the old
  pair/target correction tables assumed.
  The delivered generic fix keeps the validated current SC correction and the
  existing pair/target compaction tensors, but rebuilds the single-CGL self
  PMF around a bounded reference-ensemble state contract.
  The compact/extended centers now come from the live DOPC reference ensemble,
  the runtime prior remains the calibrated stored probability, and the self
  PMF stays inside the endpoint-interpolation family instead of opening a
  runaway high-`q` compact basin.
- Keep the prepared-file compaction state in the runtime hidden coordinate.
  Fresh hybrid stage-7 preparation still seeds
  `cgl_compaction_state/value` from physical tail-extension Angstroms even
  though the runtime node expects the normalized compaction coordinate defined
  by `compact_state_center_ang` / `extended_state_center_ang`.
  Once the repaired self PMF extends outside `[0, 1]`, the old injector
  heuristic silently skips the physical-to-runtime mapping and the C++ runtime
  clamps every lipid to the upper hidden-state bound on startup.
  The injector must affine-map from physical reference centers to runtime
  state centers regardless of the self-range bounds.
- Keep the force field generic.
  `1afo` and `1rkl` must share the same H5 tables and runtime logic.
  Do not solve the regression with system-specific settings, alternate
  workflows, or any protein-specific parameter branch.
- Preserve the runtime hidden-state design.
  The coarse runtime should still use the existing normalized endpoint tables
  (`extended=0`, `compact=1`) for SC/target/pair interpolation, but the self
  PMF and restart range must not be hard-clipped to that interval when the
  repaired physical state contract needs a compression-side response.
  The fix is to derive the physical compact/extended centers from the same
  reference ensemble being fit and then normalize the self PMF on that bounded
  coordinate. That preserves compression response physically without
  extrapolating the shared SC/target/pair tables beyond their trained compact
  endpoint.
- Validate against the full-resolution morphology, not just coarse stability.
  For `1rkl`, the acceptance target is not only “no vertical collapse”; the
  coarse stage-7 output also must not show the spurious helix bending that is
  absent from the matching full-resolution stage-7 output.
  The acceptance target also includes local secondary-structure stability at
  the bilayer leaflet interface rather than only global chain shape metrics.

# Execution Phases

- [completed] Compare coarse and full-resolution `1rkl` stage-7 outputs with a
  local secondary-structure diagnostic focused on the bilayer leaflet
  interface rather than only a global bend metric.
- [completed] Reproduce the `1rkl` leaflet-interface instability from the
  saved restart state and isolate the remaining fresh-start failure to the
  prepared-file compaction-state coordinate contract.
- [completed] Implement the minimal physical generic fix without adding any
  protein-specific branch or orientational restraint.
  The generator now prefers reference-ensemble compaction centers in `auto`,
  preserves the calibrated runtime compact prior, and keeps stored-contract
  fallback only as a fallback path.
- [completed] Rebuild the installed `parameters/dryMARTINI/dopc.h5`, reinject
  the affected stage-7 files with the corrected runtime compaction-state
  seeding, and validate the same installed H5 on both `1rkl` and `1afo`.
- [completed] Update `plan.md`, `findings.md`, and `progress.md` with the
  final diagnosis, implementation, and validation evidence.

# Known Errors / Blockers

- None for the delivered bounded-compaction path.

# Review

- Installed H5:
  `parameters/dryMARTINI/dopc.h5` now carries the bounded reference-ensemble
  compaction contract with
  `reference_extended_center_ang=12.13485 A`,
  `reference_compact_center_ang=26.96317 A`,
  `compact_state_probability=0.3023485`,
  and `self_coord_min/max=-0.06337/1.08474`.
- `1rkl` validation from fresh reinjected stage-7 start:
  CA RMSD to the full-resolution control improved from `3.55 A` on the bad
  coarse source to `2.81 A`, the local interface line-fit score improved from
  `0.1350` to `0.0940`, and late H-bonds improved from `23.26` to `29.41`
  (full-resolution control `29.50`).
- `1afo` validation from fresh reinjected stage-7 start:
  total CA RMSD to the full-resolution control improved from `4.25 A` to
  `2.91 A`, chain A RMSD improved from `2.51 A` to `1.74 A`, chain B RMSD
  improved from `3.70 A` to `2.71 A`, and late H-bonds improved from `73.65`
  to `80.21` (full-resolution control `84.17`).
