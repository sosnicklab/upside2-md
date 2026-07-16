# Progress

## Current Phase

- Added
  [add_single_cgl_compaction_to_dopc_h5](/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py:9970)
  so a DOPC H5 can be copied and have only its single-CGL overlay payloads
  refreshed, preserving the existing SC-CGL and CGL-target base tensors.
- Backed up the bad rebuilt target file to
  `parameters/dryMARTINI/dopc.h5.bad_target_rebuild_20260716-082037.bak`.
- Rebuilt `parameters/dryMARTINI/dopc.h5` from
  `parameters/dryMARTINI/dopc.h5.pre_target_cgl_compaction_fix_20260715-231157.bak`
  with `include_sc=False, include_target=True`.
- Verified the repaired target artifact on Thursday, July 16, 2026:
  `cg_lipid_target/interaction_param` and `reference_energy_eup` are exactly
  identical to the stable backup, while the target group now carries
  `delta_extended`, `delta_compact`, `delta_compressed`, and the isolated-DOPC
  single-CGL metadata.
- Ran a fresh 1AFO stage-6 rerun into
  `example/16.MARTINI/outputs/martini_1afo_hybrid_stage60_overlay_fix_20260716/`.
  Stage 6.0 completed successfully, `cgl_wrapped_z_outliers(...)` returned
  `[]`, and the run produced
  `checkpoints/1afo.stage_6.0.up` plus `1afo.stage_6.0.vtf`.
- Ran a repaired-force-field 1RKL production continuation into
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/checkpoints/1rkl.stage_7.1.up`
  and sampled 200 output frames before manual stop for analysis.

## Validation Snapshot

- 1AFO stage 6.0:
  no wrapped-z CGL outliers; the exact failure reported in
  `example/16.MARTINI/1afo.out` did not recur.
- 1RKL continuation:
  protein RMSD `3.31 A`, CA RMSD `3.49 A`, protein Rg `11.73 -> 11.82 A`,
  Rg range `11.27 .. 12.37 A`, mean CGL compaction `0.479 -> 0.473`.
- 1RKL runtime log over the sampled production segment kept H-bonds in the
  `~11-12` range and did not show the runaway protein-potential behavior seen
  with the previous target rebuild.

## Next

- Run a full fresh 1AFO stage-7.0 production validation if the remaining
  chain-tilt behavior needs to be rechecked under the repaired overlay.
