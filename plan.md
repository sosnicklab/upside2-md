# Project Goal

Fix the dryMARTINI hybrid interface so the coarse model resolves the `1afo`
collapse without introducing the coarse-only bent-protein artifact seen in
`1rkl`, while keeping one generic physical force field for both systems.

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
  The delivered generic fix keeps the validated current SC correction, restores
  the old pair/target compaction tables, and repairs the self PMF around that
  same old physical state contract so both extension and compression response
  remain active.
- Keep the force field generic.
  `1afo` and `1rkl` must share the same H5 tables and runtime logic.
  Do not solve the regression with system-specific settings, alternate
  workflows, or any protein-specific parameter branch.
- Preserve the runtime hidden-state design.
  The coarse runtime should still use the existing normalized endpoint tables
  (`extended=0`, `compact=1`) for SC/target/pair interpolation, but the self
  PMF and restart range must not be hard-clipped to that interval when the
  repaired physical state contract needs a small compression-side overshoot.
- Validate against the full-resolution morphology, not just coarse stability.
  For `1rkl`, the acceptance target is not only “no vertical collapse”; the
  coarse stage-7 output also must not show the spurious helix bending that is
  absent from the matching full-resolution stage-7 output.

# Execution Phases

- [completed] Compare coarse and full-resolution `1rkl` stage-7 outputs and
  quantify the bent-protein artifact with a shape-sensitive metric rather than
  a coarse tilt-only check.
- [completed] Reproduce the `1rkl` bend from the saved restart state and
  isolate whether it comes from the current generic compaction-table repair or
  from a broader protein-facing force-field design issue.
- [completed] Identify the responsible generic force-field term and implement
  the minimal physical fix without adding any protein-specific branch.
- [completed] Rebuild the affected H5 tables and validate the same generic fix
  on both `1rkl` and `1afo`.
- [completed] Update `plan.md`, `findings.md`, and `progress.md` with the final
  diagnosis, implementation, and validation evidence.

# Known Errors / Blockers

- None for the delivered H5 update.

# Review

- The current generic H5 fix removes most of the coarse-only `1rkl` bending by
  restoring the old compression-side compaction contract in the pair/target
  tables while keeping a nonzero repaired self PMF.
- The same H5 remains generic across `1rkl` and `1afo`; no protein-specific
  parameters or workflow branches were introduced.
