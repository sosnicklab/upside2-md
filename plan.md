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
- Keep the bounded reference-ensemble H5 contract.
  The installed `parameters/dryMARTINI/dopc.h5` is still the right generic
  self-PMF repair: it preserves the calibrated compact prior and the shared
  pair/target tensors while restoring a bounded compression-support contract at
  the table-generation layer.
- Treat the remaining stock `1rkl` bend as a compressed-regime pair-table
  defect, not as a stale injector/runtime bug.
  The direct stock-surface replay with the new self-PMF source keeps the
  hidden state in the physical compressed regime during burn-in, but the
  production run still degrades once `q` settles around `1.5..1.85`.
  That means the hidden-state source mismatch is real but secondary:
  after it is fixed, the remaining failure is the old two-state
  `cg_lipid_pair` compaction model extrapolating past its compact endpoint.
- Keep the force field generic.
  `1afo` and `1rkl` must share the same H5 tables and runtime logic.
  Do not solve the regression with system-specific settings, alternate
  workflows, or any protein-specific parameter branch.
- Rebuild one consistent compaction contract across self, SC, target, and pair
  tables.
  The live `cg_lipid_compaction` group now uses the new bounded reference
  centers, but `cg_lipid_sc` and `cg_lipid_target` still advertise the older
  physical compact / extended centers and provenance.
  The next fix must restore one generic compact-side response across those
  groups rather than changing runtime rules per interaction family.
  The stock-stage replay evidence now shows that simple two-state endpoint
  remaps are not sufficient: they can make the compact-side response stronger,
  but they either destabilize `1afo`, destroy `1rkl` backbone H-bonds, or both.
- Revised decision: extend the single-CGL compaction model with an explicit
  compressed-tail branch above the current compact endpoint.
  Source inspection showed that the retrofit builder deliberately clips the
  self-state coordinate back into the two-endpoint interval when the centers
  come from the same reference ensemble, and the SC/target runtime only mixes
  `delta_extended` and `delta_compact` linearly in `q`.
  The next candidate must therefore preserve the unclipped reference-ensemble
  high-`q` support in the self PMF and add an explicit compressed-state table
  for SC/target so the runtime no longer has to extrapolate the compact-side
  slope indefinitely.
- Revised decision: the compressed-tail branch must use a representative
  compressed sampled-state center, not the single most extreme high-`q`
  reference conformer or the median of three sparse static outliers.
  The first three-branch candidate improved `1rkl` hbonds and compaction, but
  deeper inspection showed that the installed reference ensembles only contain
  three samples above the compact center and the isolated-DOPC conformer pool
  contributes none there.
  The next candidate must therefore derive the compressed branch from the
  stored full-DOPC sampled compaction distribution (`target_compaction_ang`)
  and synthesize matching representative geometries from the compact refs with
  the same signed axial tail mode used by the physical relaxation model.
- Validate against the full-resolution morphology, not just coarse stability.
  For `1rkl`, the acceptance target is not only “no vertical collapse”; the
  coarse stage-7 output also must not show the spurious helix bending that is
  absent from the matching full-resolution stage-7 output.
  The acceptance target also includes local secondary-structure stability at
  the bilayer leaflet interface rather than only global chain shape metrics.
- Revised decision: extend the CGL-CGL pair compaction model to the same
  three-state basis as SC/target.
  After the self PMF is rebuilt from the stored full-DOPC compaction samples,
  the live `1rkl` production replay enters a real compressed regime instead of
  pinning at the self boundary.
  At that point the remaining regression comes from the pair correction path,
  which still only has the old `EE / EC / CC` tensor family and therefore
  extrapolates beyond support when `q > compact_state_center`.
  The next candidate must add explicit `EX / CX / XX` pair tensors plus the
  matching runtime interpolation, using the same generic compressed-state
  reference ensemble as the SC/target branch.
- Revised decision: the compressed pair branch must include the physical
  pair-conditioned tail-relaxation response, not only the rigid pair geometry.
  A fast diagnostic candidate that preserved the old validated `EE / EC / CC`
  tensors and added only rigid-state `EX / CX / XX` compressed tables improved
  early burn-in but failed `1rkl` late production:
  the final stage drifted back to worse-than-stock bend and H-bond metrics.
  The next candidate must therefore carry the compressed pair states through
  the same relaxation-corrected force path as the old states rather than
  treating compressed support as a rigid-only overlay.
- Revised decision: the sampled-center `32.599 A` contract can be propagated
  analytically to self, SC, target, and the existing compressed pair payload,
  but that analytical pair-center remap is not sufficient for final
  acceptance.
  Reparameterizing the exact-v8 compressed tables to the new sampled center
  materially improves `1rkl` H-bonds and bend relative to bad stock, but the
  final production segment still over-extends the helix relative to the
  full-resolution control.
  The next candidate must therefore keep the updated sampled-center
  self/SC/target contract and rebuild the compressed `EX / CX / XX` pair
  branch exactly at that center rather than only reweighting the old
  `29.17 A` compressed pair tensors.
- Accepted runtime fix: terminal-state compaction interpolation must saturate,
  not extrapolate, once the hidden state leaves the representative interval.
  Source inspection plus stock-surface replay showed that the shared
  `compute_single_compaction_mix()` routine kept extending the linear
  compact-to-compressed segment for `q > compressed_state_center`, which gives
  negative compact weights and overweighted compressed tensors in the same SC,
  target, and pair paths.
  The accepted generic implementation clamps the three-state mix only at the
  physical terminal states while preserving the same H5 tables and parameters
  for both systems.
- Accepted H5 choice: keep the `v13` sampled high-population compressed center
  (`30.347 A`) with the runtime saturation fix and reject the higher-center
  `v11` H5.
  With the runtime bug fixed, `v11` still weakens `1rkl` late production:
  `hbonds ~= 28.88`, `bend RMS ~= 4.08 A`, `Rg ~= 12.81 A`.
  The `v13` H5 remains the better generic balance across the two validated
  stock stage-7 workflows.

# Execution Phases

- [completed] Recheck the actual stock `1rkl` stage-7 artifacts, not only the
  temporary reinjected restart surface.
  The stock prepared file already uses the live bounded compaction-state
  contract, but the saved production trajectory still loses H-bonds and bends.
- [completed] Revert the rejected runtime-only clamp experiment.
- [completed] Isolate the stale-table / remap hypotheses with stock-surface
  stage-7 replay against multiple H5 candidates.
- [completed] Extend the single-CGL self/SC/target compaction response with a
  physical compressed-tail branch and a full-DOPC self-PMF source.
- [completed] Extend the CGL-CGL pair compaction correction to the same
  three-state compressed basis, including the compressed-state pair-relaxation
  response, and validate it on the stock stage-7 surface.
- [completed] Rebuild `parameters/dryMARTINI/dopc.h5` with the validated
  three-state self/SC/target/pair contract, then rerun the stock-surface
  `1rkl` and `1afo` stage-7 validations from the real workflow inputs.
- [completed] Update the tracker files with the validated stock-workflow
  evidence and any remaining blockers.

# Known Errors / Blockers

- No open implementation blocker remains for the accepted generic fix path.
- Rejected final-stage alternatives:
  full runtime clamp, pair-only runtime clamp, two-state endpoint remaps,
  rigid-only compressed pair overlays, and the higher-center `v11` H5 all
  regress at least one stock stage-7 workflow.

# Review

- Accepted source fix:
  `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
  now saturates the shared three-state compaction mix at the terminal extended
  and compressed branches instead of extrapolating beyond them.
- Accepted default H5:
  `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  now matches
  `/private/tmp/dopc_20260706_sc_target_pair_compressed_branch_candidate_v13_targetmean_reparam.h5`
  with the previous file backed up under
  `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/backup/install_20260706_123918_satmix_runtime_fix/dopc.h5`.
- Stock stage-7 validation used the real prepared-file workflow surface with
  fixed stock seeds:
  `1rkl` seed `3411623971`,
  `1afo` seed `1771339876`.
- Accepted validation artifacts:
  `/private/tmp/cgl_compression_branch_validation_1rkl_v13_targetmean_reparam_satmix_r5`
  and
  `/private/tmp/cgl_compression_branch_validation_1afo_v13_targetmean_reparam_satmix_r1`.
- Summary metrics from those accepted runs:
  `1rkl` late half:
  `hbonds ~= 30.43`,
  `Rg ~= 12.39 A`,
  `bend RMS ~= 3.91 A`,
  `end-to-end ~= 40.65 A`.
  `1afo` late half:
  `hbonds ~= 78.07`,
  `Rg ~= 15.13 A`,
  `bend RMS ~= 7.35 A`,
  `end-to-end ~= 45.96 A`.
