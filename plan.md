# Project Goal

Fix the dryMARTINI hybrid force field so protein-facing CGL-target interactions
gain the missing single-CGL tail-compaction relief without changing the legacy
`cg_lipid_target` base tensor that previously kept the protein structure stable.

# Architecture & Key Decisions

- Keep the conservative physics unchanged apart from restoring the missing
  target-side single-CGL hidden-state overlay. Do not change timestep,
  temperature, masses, damping, or add restraints.
- Preserve the accepted July 15, 2026 legacy `cg_lipid_target` base table
  exactly. The repair must add the target `delta_*` payload on top of that base
  tensor instead of rebuilding the target interaction from a different source
  ensemble.
- Use the existing single-CGL retrofit machinery in
  [py/martini_build_tables.py](/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py)
  to refresh overlays in-place on an existing DOPC H5. The new public entry
  point should copy an existing H5 and only rewrite the single-CGL compaction
  payload plus associated metadata.
- Keep the SC-CGL and CGL-CGL contracts unchanged. The repaired target overlay
  must use the same isolated-DOPC single-lipid reference contract already used
  by `cg_lipid_sc` and `cg_lipid_compaction`.
- Validate with artifact invariants first, then runtime checks:
  1AFO stage-6 wrapped-z geometry must pass; 1RKL must retain protein structure
  by RMSD/Rg/H-bond metrics, not tilt/depth alone.

# Execution Phases

- [x] Identify why the previous target fix destabilized the protein and compare
  the repaired H5 against the last stable target-base artifact.
- [x] Add a source-level helper that refreshes single-CGL overlays on top of an
  existing DOPC H5 without rebuilding the target base tensor.
- [x] Rebuild `parameters/dryMARTINI/dopc.h5` from the stable pre-fix backup,
  backing up the bad rebuilt file first.
- [x] Verify the repaired H5 invariants:
  `cg_lipid_target/interaction_param` and `reference_energy_eup` must remain
  exactly equal to the stable backup, while target `delta_*` datasets and
  hidden-state attrs are added.
- [x] Validate 1AFO stage 6.0 directly against the failure in
  `example/16.MARTINI/1afo.out`.
- [x] Validate 1RKL structural stability with a production continuation segment
  using the repaired H5.
- [x] Update `progress.md` and `findings.md` with the corrected diagnosis and
  the new validation criteria.

# Known Errors / Blockers

- A full fresh 1AFO stage-7.0 production rerun has not been completed yet under
  the repaired H5. Thursday, July 16, 2026 validation confirms the original
  stage-6 failure is fixed, but the remaining long-production 1AFO tilt behavior
  still needs a dedicated rerun before calling the interface fully closed.
- The configured master-repository reference path exists, but it does not
  contain the current dryMARTINI builder file layout needed for a direct
  line-by-line comparison. The artifact/runtime contract in this workspace was
  therefore used as the implementation baseline.

## Revised Decisions

- Do not rebuild `cg_lipid_target` from the direct builder path when applying
  the single-CGL target fix unless the legacy target source contract can be
  reproduced exactly. The previous repair failed because it replaced the
  accepted 12-conformer target base tensor with a 2-conformer tensor from
  `ref_bead_positions_nm`.
- Treat the stable backup
  `parameters/dryMARTINI/dopc.h5.pre_target_cgl_compaction_fix_20260715-231157.bak`
  as the accepted target-base contract for this repair.

# Review

- Repaired artifact:
  `parameters/dryMARTINI/dopc.h5` now preserves the legacy target base tensor
  exactly and adds
  `cg_lipid_target/{delta_extended,delta_compact,delta_compressed}` plus the
  explicit single-CGL target metadata.
- 1AFO validation on Thursday, July 16, 2026:
  fresh stage 6.0 rerun under
  `example/16.MARTINI/outputs/martini_1afo_hybrid_stage60_overlay_fix_20260716/`
  completed without wrapped-z CGL rejection and produced a valid
  `1afo.stage_6.0.up`.
- 1RKL validation on Thursday, July 16, 2026:
  continuation checkpoint
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/checkpoints/1rkl.stage_7.1.up`
  remained structurally stable over the sampled production segment, with
  protein RMSD `3.31 A`, CA RMSD `3.49 A`, and Rg `11.73 -> 11.82 A`.
