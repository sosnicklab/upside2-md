# Findings

## Root Cause

- On Thursday, July 16, 2026, the failed production artifacts were
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf` and
  `example/16.MARTINI/1afo.out`.
- The original target bug was real: protein-facing `cg_lipid_target` lacked the
  explicit single-CGL `delta_*` payload even though SC-CGL already used
  `cgl_compaction_state`.
- The previous repair was still wrong because it changed the target base
  physics. The rebuilt `parameters/dryMARTINI/dopc.h5` replaced the stable
  legacy target tensor with a different source contract:
  new target attrs were
  `reference_dataset_name=ref_bead_positions_nm`, `conformer_count=2`,
  `head_tail_span_ang=26.228395`, `n_radial=89`;
  the stable pre-fix target attrs were
  `reference_dataset_name=isolated_dopc_intramolecular_mc`,
  `conformer_count=12`, `head_tail_span_ang=19.912987`, `n_radial=91`.
- That contract change is large enough to explain both user-reported regressions
  from Thursday, July 16, 2026:
  1RKL lost structural integrity, and 1AFO failed during stage 6.0 with
  wrapped-z CGL outliers.
- The physically correct target repair is therefore:
  preserve the accepted legacy `cg_lipid_target` base tensor exactly, and add
  the missing single-CGL target overlay on top of it using the same isolated
  DOPC compaction-state contract already used by `cg_lipid_sc` and
  `cg_lipid_compaction`.

## Implementation Result

- Added
  [add_single_cgl_compaction_to_dopc_h5](/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py:9970)
  to copy an existing DOPC H5 and refresh only the single-CGL overlays.
- Rebuilt `parameters/dryMARTINI/dopc.h5` from
  `parameters/dryMARTINI/dopc.h5.pre_target_cgl_compaction_fix_20260715-231157.bak`
  while backing up the bad rebuilt file as
  `parameters/dryMARTINI/dopc.h5.bad_target_rebuild_20260716-082037.bak`.
- Verified the repaired H5 invariants:
  `cg_lipid_target/interaction_param` is exactly equal to the stable backup,
  `reference_energy_eup` is exactly equal to the stable backup,
  and the repaired target group now carries
  `delta_extended`, `delta_compact`, `delta_compressed`,
  `single_cgl_state_reference_dataset=isolated_dopc_intramolecular_mc`,
  `compaction_state_source=isolated_dopc_mc_stored_pair_contract_matched`,
  and the stored compressed-state metadata.

## Validation Result

- 1AFO stage 6.0 rerun on Thursday, July 16, 2026 under
  `example/16.MARTINI/outputs/martini_1afo_hybrid_stage60_overlay_fix_20260716/`
  completed cleanly.
  `cgl_wrapped_z_outliers(...)` returned `[]`, so the exact wrapped-z geometry
  rejection from `example/16.MARTINI/1afo.out` is no longer present.
- 1RKL repaired-force-field continuation on Thursday, July 16, 2026 under
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/checkpoints/1rkl.stage_7.1.up`
  remained structurally stable over the sampled production segment:
  protein RMSD `3.31 A`, CA RMSD `3.49 A`, protein Rg `11.73 -> 11.82 A`,
  Rg range `11.27 .. 12.37 A`, mean CGL compaction `0.479 -> 0.473`.
- The runtime log for that 1RKL segment kept H-bonds around `11-12` and did
  not show the protein blowup seen with the previous target rebuild.

## Lessons

- When adding a hidden-state conservative correction to a legacy force-field
  artifact, first prove whether the accepted behavior depends on the existing
  base tensor. If it does, overlay the correction; do not rebuild the base
  interaction from a different reference contract.
- For this interface, validating only tilt/depth is insufficient. After any
  target-CGL change, also check:
  1AFO stage-6 wrapped-z geometry,
  1RKL protein RMSD/Rg,
  and a secondary-structure proxy such as H-bond retention.
- If a user reports that a “fix” destabilized the protein, treat that as a
  correction to the validation standard and update the success criteria before
  making more force-field changes.
