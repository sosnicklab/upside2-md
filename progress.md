# Progress Log

## Current Task: Fully Physical CGL-SC Replacement Attempt

- Actions taken:
  - Implemented a pooled-reference SC rebuild path in
    `py/martini_build_tables.py` and `py/martini_gen_params.py` so a single
    `cg_lipid_sc` retrain can combine shell-filtered full-resolution DOPC
    conformers with representative isolated-DOPC conformers under one generic
    fitting rule.
  - Verified the new code path with a one-row `SER` probe overlay:
    shell-filtered + isolated reference loading, SC overlay fit, single-CGL
    compaction retrofit, and provenance attrs all completed cleanly.
  - Built a first full pooled-reference candidate
    `/private/tmp/dopc_shell15xy3d18_plusisolated_scfull_20260702.h5`
    using `4` shell conformers plus `4` isolated conformers while preserving
    the current shell-derived target tables.
  - Verified the candidate runtime path on a fresh
    `1rkl stage_6.0 -> stage_7.0` replay preparation:
    `cg_lipid_pair`, `cg_lipid_target`, and `cg_lipid_rotamer_sc` all use
    `cgl_compaction_state` with no retained implicit pair-response attrs once
    stage preparation completes.
  - Ran the matched `1rkl` replay gate for the first pooled candidate and
    rejected it relative to the installed hydromix baseline.
  - Added fail-loud validation for overlay rebuild completeness in
    `py/martini_build_tables.py` and `py/martini_prepare_system_lib.py` after
    an incomplete pooled build silently dropped the explicit SC compaction path.
  - Rebuilt the second pooled candidate cleanly as
    `/private/tmp/dopc_shell15xy3d18_plusiso2x4_scfull_retry_20260702.h5`,
    verified the explicit physical runtime path again on a fresh replay prep,
    and ran the matched `1rkl` replay gate.
  - Rejected the second pooled candidate as well.
    Its production replay lost a large fraction of late hbonds and worsened CA
    RMSD materially relative to the hydromix baseline.
  - Pivoted the active replacement route to shell-family completion:
    rebuild the remaining nonhydrophobic SC rows from the same shell-filtered
    reference ensemble already used by the accepted target and hydrophobic SC
    rows.
  - Built a matched reinjection replay harness that reuses a successful
    prepared `stage_7.0` template, swaps in a candidate H5 with
    `inject_cg_lipid_nodes()`, and runs the same shared `stage_6.0` handoff
    and production seed for cleaner candidate A/B tests.
  - Used that harness to test a targeted nonhydrophobic-only pooled candidate
    built from the physical `4 shell + 4 isolated` SC fit while preserving the
    installed hydrophobic shell rows.
  - Rejected that targeted candidate as well:
    despite high hbonds and near-baseline CGL transport, it still drove the
    protein too vertical and too expanded in late production.
  - Re-ran the original full pooled `4 shell + 4 isolated` candidate on the
    same matched reinjection replay harness and rejected it there too.
    The fair A/B confirmed the main regression is protein expansion and
    vertical alignment in production, not a transport-timescale artifact from
    the earlier replay setup.
  - Pivoted the next physical iteration toward the shell reference ensemble
    itself:
    rebuild the nonhydrophobic shell rows from a broader full-trajectory shell
    conformer pool rather than only the late-half shell around already
    verticalized proteins.
  - Built the broadened full-trajectory nonhydrophobic-shell candidate
    `/private/tmp/dopc_shell15xy3d18_start0_nonhydrophobic_candidate_20260702.h5`
    and ran the matched `1rkl` reinjection replay gate with the shared
    stage-6 handoff and production seed.
  - Rejected that broadened nonhydrophobic-shell candidate.
    It stayed compact and well hydrogen-bonded and preserved the CGL
    transport rate, but it became much more vertical by every tested
    orientation surrogate.
  - Tested the broadened hydrophobic-only shell candidate and rejected it
    quickly after the matched `1rkl` replay showed early expansion and hbond
    loss.
  - Built the broadened all-row shell candidate
    `/private/tmp/dopc_shell15xy3d18_start0_allshell_candidate_20260702.h5`
    and ran the matched `1rkl` reinjection replay gate.
  - Rejected the broadened all-row shell candidate as well.
    It kept hbonds high, compacted the protein strongly, and preserved the
    CGL lateral transport timescale, but it drove the protein almost perfectly
    vertical by every tested axis surrogate.
  - Inspected the shell-reference reduction rule itself and found a concrete
    upstream bias:
    trajectory-derived DOPC representatives are chosen by exact compaction-rank
    endpoints, so every reduced shell ensemble always includes the absolute
    minimum and maximum compaction conformers.
  - Started one more shell-plus-isolated rebuild before that issue was fully
    understood:
    `/private/tmp/dopc_shell15xy3d18_start0_plusiso4_scfull_20260702.h5`
    with the broadened shell pool plus `4` isolated representatives.
    That long build is still in flight, but it uses the endpoint-biased
    selector and is no longer the cleanest next branch.
  - Replanned the next physical iteration around the representative-selection
    fix instead of continuing to stack more candidates on the same biased shell
    reduction rule.
  - Patched `py/martini_build_tables.py` so reduced DOPC reference ensembles
    are chosen from equal-probability compaction bins rather than from the
    exact minimum and maximum compaction ranks.
  - Verified the patch directly on the pooled shell reference:
    the broadened shell reduced set moved from `4.215-30.316 A` with an
    `8/4` (`1rkl/1afo`) split to `12.896-26.126 A` with a `6/6` split.
  - Restarted the broadened shell-plus-isolated all-row rebuild as
    `/private/tmp/dopc_shell15xy3d18_start0_plusiso4_scfull_midquant_20260702.h5`
    so the next candidate isolates the representative-selection fix without
    adding a new interaction model.
  - Built and executed a 1-step `1rkl` stage-7 smoke harness from the shared
    stage-6 source to verify that the matched replay path still works before
    spending a full candidate build on it.
  - The exact midquant shell-plus-isolated rebuild completed the expensive base
    SC overlay, but the follow-on explicit-compaction retrofit remained too
    expensive to finish cleanly in this turn.
    The incomplete file was kept only as a source of the new `interaction_param`
    table, not as a physically valid final artifact.
  - Created a clearly marked screening-only proxy by combining that new
    midquant shell-plus-isolated `interaction_param` table with the closest
    completed shell-plus-isolated delta datasets, then ran a matched `1rkl`
    replay to decide whether the selector fix was directionally promising.
  - Rejected the shell-plus-isolated midquant branch on that screen.
    The production replay was already expanding by `5-6` time units, so there
    is no reason to spend another long exact-retrofit pass on this family
    before changing the reference ensemble itself.
  - Audited the next pooled full-bilayer branch and found another upstream
    bias:
    without source balancing, the reduced `4`-conformer pooled trajectory
    reference was split `3/1` across `1rkl/1afo` by raw frame count.
  - Patched `py/martini_build_tables.py` so pooled trajectory references are
    balanced per source trajectory before concatenation.
    The preferred full-bilayer reduced set is now `2/2` across
    `1rkl/1afo`, with compactions
    `17.143, 17.800, 23.063, 23.497 A`.
  - Started the exact source-balanced full-bilayer plus isolated rebuild as
    `/private/tmp/dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_20260702.h5`.
  - Diagnosed the first unreadable temp H5 from that branch as a workflow
    collision, not a physics failure:
    multiple whole-build `martini_gen_params.py` processes were writing the
    same output path at once.
  - Killed the stale competing builds and restarted one clean exact rebuild on
    the same source-balanced branch with a single output writer.
  - Revalidated the replay-metric path directly from the stage-7 H5 files.
    The hydromix baseline reproduces
    `late hbonds ~= 30.62`,
    `late CA RMSD ~= 1.76 A`,
    and the accepted CGL transport values
    `Dxy_half ~= 0.421`,
    `Dxy_quarter ~= 0.409`.
    The rejected all-shell replay reproduces the expected near-vertical
    orientation and faster transport.
  - Finished the exact source-balanced full-bilayer candidate build and ran
    the matched `1rkl` replay gate.
    The branch stayed compact and preserved the CGL transport timescale, but
    it was still too vertical to accept.
  - Built a screening proxy on the same new source-balanced base SC table that
    reused the installed SC compaction deltas.
    That proxy stayed close to baseline and did not show the exact branch's
    verticalization, isolating the remaining problem to the exact single-CGL
    SC compaction-delta retrofit rather than to the rebuilt base SC table.
  - Inspected the retrofit path and found the source mismatch directly:
    the base SC table was fit on the pooled overlay ensemble, but the explicit
    SC compact/extended delta endpoints were still being rebuilt from a
    separate isolated-DOPC MC pool.
  - Patched `py/martini_build_tables.py` so the single-CGL SC/target retrofit
    now defaults to selecting compact and extended representative conformers
    from the same overlay reference ensemble, matched to the stored physical
    compaction centers, with the isolated-MC route retained only as fallback.
  - Verified the patch with `python3 -m py_compile py/martini_build_tables.py`.
  - Built the accepted SC-only retrofit artifact
    `/private/tmp/dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_scstateoverlay_scintf_nologin8_20260702.h5`
    from the source-balanced full-bilayer base H5 and confirmed that only the
    SC compact/extended delta datasets changed.
  - Re-ran the matched `1rkl` replay on that artifact and recovered a near-
    baseline orientation while preserving compactness and transport:
    `late hbonds 29.08`,
    `late Rg 12.50 A`,
    `late CA RMSD 1.86 A`,
    `axis_trim1 0.740`,
    `pca_full 0.888`,
    `Dxy_half / quarter 0.439 / 0.443`.
  - Re-ran the same artifact on matched `1afo` replay.
    The run completed cleanly through minimization, 40k burn-in, and 10k
    production with seed `20260703`, then measured
    `late hbonds 83.96`,
    `late Rg 14.83 A`,
    `late CA RMSD 2.68 A`,
    `overall pca_full 0.996`,
    `Dxy_half / quarter 0.459 / 0.456`.
  - Compared the accepted candidate directly against the installed
    `parameters/dryMARTINI/dopc.h5` and verified that every bilayer-only CGL
    dataset and attr is bitwise identical.
    The accepted promotion step changes only
    `cg_lipid_sc/delta_extended` and `cg_lipid_sc/delta_compact`, so the
    standalone bilayer timescale is preserved exactly.
  - Promoted the accepted candidate into
    `parameters/dryMARTINI/dopc.h5` and refreshed
    `parameters/dryMARTINI/dopc.h5.bak` from the previous default artifact.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_gen_params.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system_lib.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - The one-row pooled-reference probe completed end to end and produced the
    expected provenance:
    shell target preserved,
    `dopc_interface_sc_reference_source =
    pooled_full_resolution_upside_shell_filtered_conformer_trajectory_plus_isolated_dopc_representative_ensemble`,
    and the isolated component counts/seed were recorded.
  - The first full pooled candidate replayed `1rkl stage_6.0 -> stage_7.0`
    successfully with the explicit compaction-state runtime path confirmed in
    the prepared file after stage preparation finished.
  - Matched replay metrics for the first full pooled candidate versus the
    hydromix baseline:
    `late hbonds 29.67 vs 30.62`,
    `late Rg 12.63 vs 12.76 A`,
    `late CA RMSD 2.14 vs 1.76 A`,
    `axis_trim1 0.792 vs 0.730`,
    `Dxy_step_half 0.433 vs 0.421`,
    `Dxy_step_quarter 0.440 vs 0.409`.
    Interpretation:
    the candidate stayed physically clean at runtime but regressed the `1rkl`
    replay gate, so it was not carried forward to `1afo`.
  - Matched replay metrics for the second pooled candidate versus the
    hydromix baseline:
    `late hbonds 19.74 vs 30.62`,
    `late Rg 12.38 vs 12.76 A`,
    `late CA RMSD 2.63 vs 1.76 A`,
    `axis_trim1 0.720 vs 0.730`,
    `Dxy_step_half 0.437 vs 0.421`,
    `Dxy_step_quarter 0.439 vs 0.409`.
    Interpretation:
    more isolated bias makes the replay substantially less stable, so the
    pooled isolated route is not the active path forward.
  - Matched replay metrics for the nonhydrophobic-only pooled candidate versus
    the hydromix baseline on the shared seed-controlled replay harness:
    `late hbonds 34.26 vs 30.62`,
    `late Rg 13.30 vs 12.76 A`,
    `late CA RMSD 2.71 vs 1.76 A`,
    `axis_trim1 0.819 vs 0.730`,
    `Dxy_step_half 0.408 vs 0.421`,
    `Dxy_step_quarter 0.411 vs 0.409`.
    Interpretation:
    keeping only the nonhydrophobic rows pooled is still not enough; the
    protein remains overexpanded and too aligned with the bilayer normal.
  - Matched replay metrics for the full pooled `4 shell + 4 isolated`
    candidate on the same seed-controlled replay harness:
    `late hbonds 28.40 vs 30.62`,
    `late Rg 13.48 vs 12.76 A`,
    `late CA RMSD 2.26 vs 1.76 A`,
    `axis_trim1 0.809 vs 0.730`,
    `Dxy_step_half 0.416 vs 0.421`,
    `Dxy_step_quarter 0.415 vs 0.409`.
    Interpretation:
    the full pooled candidate also fails on a clean replay A/B, so the next
    hypothesis is not “more isolated regularization” but “less biased shell
    reference sampling.”
  - Matched replay metrics for the broadened all-row shell candidate versus
    the hydromix baseline:
    `late hbonds 34.95 vs 30.62`,
    `late Rg 11.68 vs 12.84 A`,
    `Dxy_step_half 0.464 vs 0.421`,
    `Dxy_step_quarter 0.461 vs 0.409`,
    but orientation surrogates become far more vertical:
    `axis_ee_3_22_half 0.987 vs 0.749`,
    `axis_pca_trim1_half 0.979 vs 0.840`,
    `axis_pca_3_22_half 0.989 vs 0.839`.
    Interpretation:
    broadening the shell pool alone does not fix the generic drift; it
    amplifies the orientational bias while keeping timescale acceptable.
  - Reference-pool audit for the shell selector:
    with `12` representatives, the current reduction picks exact compaction
    percentiles `0, 9.1, ..., 90.9, 100`.
    For the broadened shell pool that means the reduced ensemble always
    includes `4.215 A` and `30.316 A`, the absolute pool extrema.
    That gives a concrete physical reason to patch the selection rule before
    more expensive rebuilds.
  - Post-patch reference-pool audit:
    the same broadened shell pool now reduces to
    `12.896, 15.148, 16.600, 17.738, 18.732, 19.715, 20.602, 21.479, 22.375, 23.406, 24.561, 26.126 A`
    with an even `6/6` source split across `1rkl` and `1afo`.
    The late-half shell pool shows the same behavior.
  - Stage-7 harness smoke test:
    a copied prepared template plus shared `stage_6.0.up` ran cleanly through
    handoff, 1-step production, and VTF extraction, so the next candidate can
    be validated immediately after the H5 build finishes.
  - Shell-plus-isolated midquant screen result:
    the matched `1rkl` replay started near the baseline (`hbonds 33.8`,
    `Rg 12.7 A`) but by steps `2550-2900` had degraded to roughly
    `hbonds 22-24` and `Rg 13.7-13.9 A`.
    That is enough to reject the branch directionally even before finishing the
    full production window.
  - Source-balanced full-bilayer reference audit:
    the preferred pooled full-bilayer reduced set is now balanced `2/2`
    across `1rkl/1afo`,
    instead of the previous raw-count-driven `3/1` split.
  - Replay-metric calibration:
    on the hydromix `1rkl` baseline H5,
    direct analysis reproduces
    `late hbonds 30.624`,
    `late CA Rg 12.843 A`,
    `late CA RMSD 1.756 A`,
    `Dxy_half 0.421`,
    `Dxy_quarter 0.409`.
    On the rejected all-shell replay, the same script reproduces
    `late hbonds 34.95`,
    `late CA Rg 11.67 A`,
    `late CA RMSD 1.52 A`,
    and near-vertical PCA surrogates, confirming that the current gate is
    measuring the same failure mode as the earlier reviews.
  - Exact source-balanced full-bilayer `1rkl` replay:
    `late hbonds 32.96`,
    `late Rg 12.48 A`,
    `late CA RMSD 1.83 A`,
    `axis_trim1_half 0.847`,
    `pca_full_half 0.929`,
    `Dxy_half 0.428`,
    `Dxy_quarter 0.430`.
    Interpretation:
    the branch is structurally stable and keeps the timescale close, but it is
    still too vertical to accept as the generic fix.
  - Proxy source-balanced `1rkl` screen with installed SC compaction deltas:
    `late hbonds 25.87`,
    `late Rg 12.77 A`,
    `late CA RMSD 1.30 A`,
    `axis_trim1_all 0.708`,
    `Dxy_all 0.421`.
    Interpretation:
    the source-balanced base SC table is not the remaining failure mode; the
    exact explicit SC compaction-delta layer is.

- Result:
  - A generic pooled-reference SC retrain path now exists in source and is
    reproducible from the CLI.
  - Both tested pooled candidates are rejected on `1rkl`, despite remaining
    physically clean at runtime.
  - The matched reinjection replay harness now makes further candidate tests
    cheaper and cleaner.
  - The shell-family completion route and the targeted nonhydrophobic pooled
    route are both rejected on `1rkl`.
  - The full pooled `4 shell + 4 isolated` candidate is also rejected on that
    same fair replay harness.
  - The broader full-trajectory nonhydrophobic-shell retry is also rejected:
    it keeps `1rkl` compact and strongly hydrogen-bonded and preserves the
    CGL timescale, but it drives the protein more vertical than baseline.
  - The broadened hydrophobic-only and broadened all-row shell retries are now
    rejected too.
  - The representative-selector patch is implemented and verified on the shell
    reference pool itself.
  - The shell-plus-isolated family remains rejected even with the corrected
    selector.
  - The source-balanced pooled full-bilayer branch is the current best
    physical base family.
  - The remaining active method is to keep that base table and replace only
    the exact SC compaction-delta retrofit with overlay-center-matched state
    selection before re-gating `1rkl` and then `1afo`.
