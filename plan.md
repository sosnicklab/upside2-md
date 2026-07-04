# Project Goal

Fix the dryMARTINI hybrid interface so the coarse `1afo` two-helix system keeps
a physical lipid-mediated separation without regressing the accepted `1rkl`
tilt behavior.

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
  Use the stored DOPC reference ensemble inside
  `parameters/dryMARTINI/dopc.h5` to define physical compact and extended tail
  centers, then propagate the same contract consistently into:
  `cg_lipid_compaction`, `cg_lipid_sc`, and `cg_lipid_target`.
- Keep the first repair physically minimal.
  The user-corrected root cause is the missing compression-side correction in
  the protein-facing force field, so rebuild `cg_lipid_compaction/self_coeff`
  plus the `cg_lipid_sc` and `cg_lipid_target` compact-vs-extended tables
  under the repaired centers first.
  Leave the existing CGL-CGL pair tensor in place unless the controlled restart
  validation still shows the `1afo` collapse.
- Preserve the runtime hidden-state design.
  The coarse runtime should use a normalized `q in [0, 1]` compaction
  coordinate, with physical tail-extension centers retained as explicit H5
  metadata for initialization and diagnostics.
- Validate from the same saved `1afo` stage-7 source state.
  Compare controlled short continuations before and after the H5 repair so the
  force-field effect is isolated from stochastic repacking noise.

# Execution Phases

- [completed] Verify that coarse and full-resolution `1afo` start from the
  same packed structure and that the regression is coarse-force-field specific.
- [completed] Confirm the compaction-state defect in the installed
  `parameters/dryMARTINI/dopc.h5`.
  The prepared coarse `1afo` state was seeded near the compact endpoint, the
  live `cg_lipid_compaction/self_coeff` had been zero, and the stored physical
  state centers were stale.
- [completed] Patch the table-generation path so compaction metadata can be
  rebuilt from the stored DOPC reference ensemble and so stale single-CGL
  groups are detectable.
- [completed] Overwrite `parameters/dryMARTINI/dopc.h5` with the minimal
  consistent repair needed for the user-corrected root cause:
  update the self compaction PMF plus the SC/target compact-vs-extended
  correction tables under the repaired physical state centers, while keeping
  the runtime model unchanged.
- [completed] Reinject a controlled `1afo` stage-7 diagnostic from the same saved
  state and verify that compaction no longer saturates at the compact endpoint
  and that helix separation improves.
- [completed] Run a short `1rkl` diagnostic to confirm the accepted tilt behavior
  is preserved.
- [in-progress] Update `plan.md`, `findings.md`, and `progress.md` with the final
  repaired diagnosis, implementation, and validation evidence.

# Known Errors / Blockers

- None at the current validation scope.

# Review

- Full-resolution and coarse `1afo` use the same packed hybrid structure, so
  the coarse collapse is not explained by missing lipids in the inter-helix
  corridor alone.
- The current installed H5 is internally inconsistent:
  `cg_lipid_compaction` now uses repaired physical centers
  (`12.13 A` extended, `26.96 A` compact), while `cg_lipid_sc` and
  `cg_lipid_target` still use stale centers (`13.38 A`, `19.62 A`).
- A controlled repaired-state continuation showed that fixing only the hidden
  state initialization and self PMF reduces compact-endpoint saturation, but it
  does not by itself recover the `1afo` helix separation.
  The remaining mismatch is in the stale compact-vs-extended pair/SC/target
  correction tables.
- The delivered repair uses the actual runtime state model instead of an
  expensive re-fit.
  The C++ runtime combines SC/target corrections linearly in the hidden state
  and combines the CGL-CGL endpoint tensors bilinearly, so the stale endpoint
  tables can be remapped exactly from the old physical centers in
  `parameters/dryMARTINI/dopc.h5.bak` to the repaired centers in the live H5.
- Validation from the saved stage-7 restart state now moves in the right
  direction for `1afo`.
  Over a 10,000-step continuation from the same saved final stage-7 frame,
  the chain-center XY separation increased from `8.41 A` to `8.79 A` and the
  mean CGL compaction dropped from `0.9999` to `0.9178` instead of staying
  pinned at the compact endpoint.
- A 10,000-step `1rkl` continuation stayed in the accepted tilted regime.
  With a simple backbone principal-axis diagnostic, the final tilt changed only
  modestly (`31.15 -> 28.88 deg` using the full-chain metric) and did not show
  a return to the obvious vertical-collapse failure mode.
