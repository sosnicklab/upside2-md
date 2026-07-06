# Progress Log

## Current Task: Accepted Generic Compression Fix

- Actions taken:
  - Patched the shared three-state compaction mixer in
    `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
    so terminal extended/compressed branches saturate instead of linearly
    extrapolating beyond support.
  - Rebuilt the C++ core and replayed the real stock stage-7 workflow surface
    with fixed stock seeds for both benchmark systems using the accepted
    `v13` compressed-branch H5.
  - Replayed and rejected the higher-center `v11` H5 on the same patched
    runtime because it regressed `1rkl` late production.
  - Backed up the old installed default H5 and promoted
    `/private/tmp/dopc_20260706_sc_target_pair_compressed_branch_candidate_v13_targetmean_reparam.h5`
    to
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - Accepted `1rkl` stock-surface stage-7 replay:
    `/private/tmp/cgl_compression_branch_validation_1rkl_v13_targetmean_reparam_satmix_r5`
    with late-half
    `hbonds ~= 30.43`,
    `Rg ~= 12.39 A`,
    `bend RMS ~= 3.91 A`,
    `end-to-end ~= 40.65 A`,
    and no compressed-branch runtime extrapolation.
  - Accepted `1afo` stock-surface stage-7 replay:
    `/private/tmp/cgl_compression_branch_validation_1afo_v13_targetmean_reparam_satmix_r1`
    with late-half
    `hbonds ~= 78.07`,
    `Rg ~= 15.13 A`,
    `bend RMS ~= 7.35 A`,
    `end-to-end ~= 45.96 A`.
  - Rejected `1rkl` follow-up on the higher-center `v11` H5:
    `/private/tmp/cgl_compression_branch_validation_1rkl_v11_targetdist_reparam_satmix_r1`
    fell to
    `hbonds ~= 28.88`,
    `bend RMS ~= 4.08 A`,
    `Rg ~= 12.81 A`,
    so `v11` is not the installed default.

- Failures and fixes:
  - The actual root-cause bug was not “missing lipid” or a missing table alone.
    The force application itself was still unphysical because the shared
    compressed-state interpolation produced negative weights once
    `q > compressed_state_center`.
    Fixing that runtime path was necessary before any H5 choice could be
    judged fairly.

## Current Task: Reopened Stock-Workflow `1rkl` Bend Fix

- Actions taken:
  - Reused the exact-v8 compressed pair candidate as a source H5 and ran the
    shared single-CGL builder against it so the new sampled compressed center
    (`32.599 A`) could propagate analytically across self, SC, target, and the
    existing compressed pair payload without another full exact pair solve.
  - Verified that the resulting candidate writes a consistent compressed-state
    contract into all three groups:
    `cg_lipid_compaction`,
    `cg_lipid_sc`,
    and `cg_lipid_target` now all advertise
    `reference_compressed_center_ang ~= 32.599 A`.
  - Replayed `1rkl` on the stock stage-7 surface with that sampled-center
    reparameterized candidate and monitored the full `40k` burn-in plus the
    `10k` production segment against the bad stock failure mode.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - The sampled-center analytical remap is real:
    the rebuilt H5 at
    `/private/tmp/dopc_20260706_sc_target_pair_compressed_branch_candidate_v11_targetdist_reparam.h5`
    carries
    `reference_compressed_center_ang ~= 32.599 A`,
    `pair_reference_compressed_center_ang ~= 32.599 A`,
    and matching SC/target compressed attrs.
  - On fixed-seed stock stage-7 replay for `1rkl`, the burn-in phase stays
    materially stronger than bad stock all the way to the end:
    H-bonds remain mostly in the `29..34` range while `Rg` stays near
    `12.7..12.8` through the whole `40k`-step burn-in.
  - The late saved-state metrics improve over bad stock:
    `hbonds ~= 31.03` vs `29.47`,
    `bend RMS ~= 3.39 A` vs `3.80 A`,
    `end-to-end ~= 45.75 A` vs `41.96 A`.

- Failures and fixes:
  - The sampled-center analytical remap is still not the final answer.
    The `1rkl` production segment expands to `CA Rg ~= 13.37 A`, above both
    the bad stock coarse run and the full-resolution control, and the aligned
    CA overlay is still too extended.
    The next iteration is an exact rebuild of the compressed
    `EX / CX / XX` pair branch on top of the improved sampled-center
    self/SC/target contract.

- Actions taken:
  - Audited the compressed-state source after the rejected exact pair
    candidate and quantified the full mismatch between the static H5
    references and the sampled full-DOPC compaction distribution.
  - Verified that the attempted isolated-DOPC supplement is a dead end for
    this bug: it contributes zero conformers above the current compact center,
    so it cannot move the compressed branch off the rejected `29.17 A`
    support point.
  - Confirmed that the stored
    `cg_lipid_compaction/target_compaction_ang` samples are the actual high-q
    physical source the self PMF already follows, with support extending to
    `39.52 A`.
  - Patched `py/martini_build_tables.py` so the compressed branch now derives
    its center from the sampled full-DOPC upper quartile and synthesizes the
    compressed representative geometries from the compact refs using a signed
    axial tail mode.
  - Patched the three-state probability split so it is computed from the
    sampled compaction distribution instead of the tiny static reference pool.
  - Added temporary helpers to rebuild a full single-CGL candidate H5 and to
    evaluate stage-7 outputs from `.up` files with late-frame H-bond, `CA`
    `Rg`, bend RMS, end-to-end, and CGL-compaction statistics.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/tmp_build_single_cgl_candidate_full.py`
  - `/Users/yinhan/Documents/upside2-md/tmp_build_pair_candidate_exact_compressed.py`
  - `/Users/yinhan/Documents/upside2-md/tmp_eval_stage70.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - The new compressed-state source is active in the builder:
    `_augment_compaction_states_with_compressed_branch(...)` now reports
    `compressed_center_source=target_compaction_upper_quartile_mean` and
    `compressed_center_ang ~= 32.599 A`.
  - The synthesized compressed representative ref lands exactly at that target
    extension (`32.599 A`) rather than inheriting the old sparse-reference
    `29.17 A` median.
  - The new H5 evaluator reproduces the known late-stage stock/full
    comparison direction for `1rkl`:
    coarse has lower H-bonds and shorter end-to-end distance than the
    full-resolution control.

- Failures and fixes:
  - The first attempt to improve compressed support by appending isolated-DOPC
    conformers was rejected before replay because source inspection showed it
    would rebuild the same old compressed branch.
    The task is now focused on the sampled full-DOPC compaction source and the
    signed axial tail deformation path instead of another sparse-reference
    supplement.

- Actions taken:
  - Replayed `1rkl` on the stock stage-7 surface after fixing the self-PMF
    source and adding compressed-state support to SC/target.
  - Verified that the new self-PMF source removes the old self-boundary trap:
    the live compaction state now occupies the physical compressed regime
    during burn-in and early production (`mean q ~= 1.56`, `p95 ~= 1.84`).
  - Compared that improved hidden-state behavior against the production
    degradation and confirmed the remaining defect is downstream of the
    self/SC/target branch.
  - Traced the remaining path in `src/martini_cg_lipid.cpp`,
    `py/martini_build_tables.py`, and `py/martini_prepare_system_lib.py` and
    confirmed that `cg_lipid_pair` still implements only the old two-state
    `EE / EC / CC` correction family.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - The self-PMF source fix is real and active:
    the prepared file and the live output now both visit compressed `q`
    values well above the compact endpoint.
  - The remaining failure is also real:
    despite that corrected hidden-state occupancy, the stock-surface `1rkl`
    production replay still loses H-bonds and bends.
  - The current source confirms the mismatch:
    pair compaction has no `EX / CX / XX` branch even though the hidden state
    can now reach that regime.

- Failures and fixes:
  - The previous “self/SC/target compressed branch” fix is necessary but not
    sufficient.
    The task is now focused on extending the pair-compaction model family
    itself rather than reshuffling the same two-state endpoint tables.

- Actions taken:
  - Rechecked the actual stock artifact
    `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf`
    after the bounded-H5 install instead of relying on the temporary replay
    surface.
  - Confirmed the stock prepared file already uses the live bounded
    compaction-state contract from `parameters/dryMARTINI/dopc.h5`:
    `q min/max/mean ~= 0.3681 / 1.0847 / 0.8149`.
  - Compared the stock prepared state, the saved stock production state, and
    the earlier temporary restart validation state; they are not the same
    workflow surface, so the earlier replay result was not sufficient for
    acceptance.
  - Inspected the stock burn-in and production logs and confirmed that the
    real saved production run still loses H-bonds and expands slightly while
    bending.
  - Traced the runtime correction paths in `src/martini_cg_lipid.cpp` and
    verified that `cg_lipid_pair`, `cg_lipid_sc`, and `cg_lipid_target` all
    use unclamped compaction interpolation coordinates.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`

- Verification:
  - The installed H5 remains correctly wired into the stock prepared file, so
    the remaining failure is not a stale-table installation problem.
  - Source inspection confirmed the remaining runtime issue:
    the SC and target corrections use raw `q` and `1-q` weights, and the pair
    correction normalizes `q` to `si/sj` without clamping before bilinear
    mixing.

- Failures and fixes:
  - The earlier temporary replay acceptance was invalid for the user-facing
    deliverable because the stock `stage_7.0.vtf` path still regressed.
    The task is reopened around the runtime interpolation rule rather than the
    H5 generation path.

## Current Task: Two-State Compression Hypotheses Rejected

- Actions taken:
  - Reverted the rejected runtime pair-clamp experiment from
    `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`.
  - Patched `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
    so stale SC/target groups are detectable when their reference centers and
    compaction-state provenance no longer match the active compaction contract.
  - Built multiple stock-stage candidate H5 files under `/private/tmp/` and
    replayed stage 7.0 from the saved stock prepared files with the stock
    seeds for both `1rkl` and `1afo`.
  - Candidates tested:
    full endpoint reparameterization of self+SC+target+pair;
    SC+target-only endpoint reparameterization;
    SC-only endpoint reparameterization;
    target-only endpoint reparameterization;
    and a self-only “selected representative center” contract.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/src/martini_cg_lipid.cpp`
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - Every candidate was checked on the stock stage-7 replay surface, not only
    on a reinjected temporary restart.
  - The common replay outcome is consistent:
    stronger compact-side weighting can shrink `1rkl` globally but also
    destroys the local H-bond network and degrades `1afo`.
  - The best “moderate” self-only center shift seeded
    `1rkl` with `q min/max/mean ~= 0.257 / 1.150 / 0.904`, but even that still
    drove the late burn-in H-bonds down into the single digits.

- Failures and fixes:
  - No tested two-state contract repair satisfies both systems.
    The remaining defect now looks like a model-family limitation rather than
    a stale H5-install or one-line remap bug.

- Actions taken:
  - Fixed the compressed-state representative-selection bug in
    `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py` so the
    supplemented upper-tail compressed center is taken from the representative
    combined pool (`29.17 A`) instead of the sparse reference-ensemble
    maximum (`31.16 A`).
  - Built and validated a fast diagnostic pair candidate that preserved the
    old validated `EE / EC / CC` pair tensors and added only rigid-state
    `EX / CX / XX` compressed pair tensors:
    `/private/tmp/dopc_20260706_sc_target_pair_compressed_branch_candidate_v7_centerfix_rigidx_preserve.h5`.
  - Replayed stock stage 7.0 for `1rkl` with the saved stock seed
    `3411623971` from the real prepared-file surface.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/tmp_build_pair_candidate.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - The corrected compressed-center selection is real:
    the builder now chooses `compressed_center_ang ~= 29.17 A`
    (`compressed_state_coord ~= 1.1488`) rather than the old singleton
    `31.16 A` outlier.
  - The diagnostic H5 writes the expected compressed pair datasets and attrs:
    `delta_extended_compressed`, `delta_compact_compressed`,
    `delta_compressed_compressed`, and
    `pair_reference_compressed_center_ang`.
  - The rigid-only compressed pair branch is rejected by the real stock
    `1rkl` stage-7 replay:
    late-half metrics degrade to
    `hbonds ~= 26.60`, `bend RMS ~= 3.32 A`,
    `end-to-end ~= 34.98 A`,
    worse than the bad stock coarse source and farther from the full control.

- Failures and fixes:
  - A compressed rigid pair overlay is not sufficient.
    The next iteration must extend the pair-conditioned tail-relaxation
    response to the compressed `EX / CX / XX` branch instead of only
    correcting the rigid compressed geometry.
