# Findings

## Lessons

- Do not assign causality from a coarse packed-geometry artifact without
  checking the matching full-resolution control workflow first.
  For this project, if full-resolution and coarse runs start from the same
  packed `1afo` structure but only coarse collapses, the root cause is not the
  static corridor occupancy by itself.
- When the master baseline does not contain the subsystem being cleaned, define
  scope from the accepted workflows and the prepared-file contract rather than
  from the raw presence of flags in the current code.
  For this project, dead-branch cleanup has to be justified by workflow
  reachability, not by guesswork about why an old toggle exists.
- Do not infer the accepted dryMARTINI correction contract from stale prepared
  files alone when the source has changed later in the same turn.
  For this project, confirm the live contract with a fresh reinjection or a
  controlled continuation from the same saved stage before deciding which path
  is physically responsible.
- When a cleanup touches a force-correction path that has a known physical role,
  do not treat the branch as dead based only on code reachability or smoke-run
  success.
  First verify that the intended correction is still active in a realistic
  stage-7 output or prepared file and that the known behavioral regression does
  not reappear.
- When the tracked force-field H5 is unchanged but a saved trajectory regresses,
  isolate the runtime surface before touching table generation again.
  A byte-for-byte `dopc.h5` match to `HEAD` means the next diagnostic should be
  injector/runtime A/B from the same saved state, not a blind retrain or H5
  rewrite.
- When a cleanup or refactor request names concrete workflows for verification,
  make those workflows part of the plan before calling the task done.
  For this project, workflow validation is part of the deliverable, not an
  optional follow-up.
- Before using trajectory arrays for transport metrics, confirm the exact H5
  shape and slice the system axis explicitly.
  A mistaken read of `output/pos` as `(frame, atom, xyz)` instead of
  `(frame, system, atom, xyz)` can fabricate a severe diffusion regression.
- Do not reach first for script-level thermostat, GLE, or workflow-parameter
  tuning when a hybrid bilayer timescale shifts after a physics change.
  First isolate which interaction family changed the effective drag and fix
  that physical coupling path before touching transport knobs.
- Do not claim a force-field retrain is delivered until the installed file path
  used by the default workflow has been updated or the exception is stated
  explicitly.
  In this project the acceptance target is usually
  `parameters/dryMARTINI/dopc.h5`, not just a generated overlay under
  `example/16.MARTINI/outputs/...`.
- When comparing stochastic hybrid workflow reruns, fix or reuse the same
  packing and simulation seeds before treating metric deltas as an H5-only
  effect.
  Fresh random starts are good for smoke testing the installed path, but they
  are not a strict A/B against an earlier validated run.
- When a candidate runtime fix spans table generation, injection, and C++
  evaluation, do not trust the source diff alone.
  Inspect a prepared `.up` file and verify that the new node arguments and
  datasets are actually present before interpreting any rerun as validation of
  the intended physics.
- When a hidden state changes coordinate systems across table-generation
  stages, inspect both the stored values and the allowed domain in the prepared
  `.up` file.
  The explicit CGL compaction runtime looked wired correctly at the node level,
  but the prepared state was still seeded from physical tail-extension values
  into a normalized `0..1` overlay coordinate, silently clamping every lipid at
  startup.
- When a validated H5 artifact becomes the default workflow input, verify that
  the current generator entrypoint can recreate its provenance-critical options.
  Helper functions and local one-off scripts are not enough; if the CLI cannot
  drive shell filtering, residue-subset SC rebuilds, or SC row merges, the
  implementation is still incomplete.
- When an exact retrain regresses but a proxy that reuses one correction layer
  does not, isolate the failure to that correction layer before rebuilding the
  whole table family again.
  For the CGL-SC issue, the exact/proxy replay split is what showed that the
  remaining vertical bias was in the explicit SC compaction deltas, not in the
  new source-balanced base `interaction_param` table.
- When a compaction correction model is already linear or bilinear in its
  hidden state, do not jump straight to an expensive representative-state
  re-fit after the physical endpoint centers move.
  First check whether the stale tables can be reparameterized exactly from the
  old physical centers to the new ones inside the existing model family.

## External / Technical Findings

- 2026-07-04: The stale CGL compaction tables can be repaired exactly by
  endpoint reparameterization; a fresh expensive re-fit is not required for
  this fix.
  - `parameters/dryMARTINI/dopc.h5.bak` preserves the old physical endpoint
    contract:
    `extended=13.381275 A`, `compact=19.615337 A`,
    `compact_state_probability=0.3023485`.
  - The repaired live
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    already carried the corrected self-PMF contract:
    `reference_extended_center_ang=12.134849 A`,
    `reference_compact_center_ang=26.963173 A`,
    `compact_state_probability=0.5833333`.
  - The C++ runtime in `src/martini_cg_lipid.cpp` combines
    `cg_lipid_sc/delta_extended` and `delta_compact` linearly in the hidden
    compaction coordinate and combines the
    `cg_lipid_compaction/delta_extended_extended`,
    `delta_extended_compact`, and `delta_compact_compact` tensors bilinearly.
  - Therefore the stale endpoint tables can be remapped analytically from the
    old physical centers to the new ones without changing the runtime model:
    the SC/target tables by affine endpoint evaluation and the CGL-CGL pair
    tables by bilinear endpoint evaluation.
  - The delivered repair now writes one consistent contract across
    `cg_lipid_compaction`, `cg_lipid_sc`, and `cg_lipid_target`, with
    `extended_state_center_ang=0`, `compact_state_center_ang=1`,
    repaired physical `reference_*` attrs, and explicit provenance markers
    (`endpoint_reparameterized_from_dopc_h5_bak` for the remapped tables).
  - Validation from the saved stage-7 restart state moved the coarse `1afo`
    system in the right direction:
    the original saved `stage_7.0` source ends at `8.41 A` chain-center XY
    separation with mean compaction `0.9999`, while the repaired 10,000-step
    continuation from that same final frame ends at `8.79 A` separation with
    mean compaction `0.9178`.
  - A 10,000-step repaired `1rkl` continuation stayed in the accepted tilted
    regime by a simple backbone principal-axis diagnostic
    (`31.15 -> 28.88 deg` on the full-chain metric), so the fix did not
    reproduce the obvious old vertical-collapse mode.

- 2026-07-04: The current `1afo` collapse is a coarse compaction-state
  force-field defect, not a missing-lipid packing defect.
  - The saved coarse
    `example/16.MARTINI/outputs/martini_1afo_hybrid/checkpoints/1afo.stage_7.0.up`
    does show real inward drift:
    the chain-center XY separation shrinks from about `11.50 A` to `8.41 A`
    across the saved frames.
  - The matching coarse and full-resolution `1afo` workflows start from the
    same packed hybrid structure:
    `example/16.MARTINI/outputs/martini_1afo_hybrid/hybrid_prep/hybrid_prep_summary.json`
    and
    `example/16.MARTINI/outputs/martini_1afo_hybrid_full/hybrid_prep/hybrid_prep_summary.json`
    report the same lipid removal and box geometry.
    The full-resolution stage-7 trajectory remains stable, so static corridor
    occupancy is not the root cause by itself.
  - The installed
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    contained a broken compaction-state contract:
    `cg_lipid_compaction/self_coeff` was zero and its physical state-center
    metadata was stale relative to the stored DOPC reference ensemble.
  - The coarse prepared `1afo` stage-7 file therefore mapped physical tail
    extensions of about `17.6..31.1 A` onto a clipped hidden coordinate with
    mean `q=0.985`, seeding almost every lipid at the compact endpoint before
    dynamics.
  - Repairing only the hidden-state initialization and self PMF is necessary
    but not sufficient.
    A controlled repaired-state continuation reduces the compact-endpoint
    fraction from about `0.989` to about `0.176`, yet the two helices still
    collapse inward.
  - The remaining defect is broader than SC/target metadata alone.
    The retrofit path had repaired the hidden-state self PMF, but it had not
    regenerated the compact-vs-extended CGL-CGL pair tensor stored inside
    `cg_lipid_compaction`, even though the prepared `cg_lipid_pair` node reads
    its compact-state deltas directly from that group.
  - `cg_lipid_sc` and `cg_lipid_target` still use stale compact/extended
    centers (`13.381 / 19.615 A`) while the repaired `cg_lipid_compaction`
    group now uses the correct physical centers (`12.135 / 26.963 A`), and the
    pair tensor lacked explicit repaired-state provenance markers as well.
    That means the compression-side and extension-side pair/SC/target
    corrections were not being applied under one consistent state definition.

- 2026-07-04: The user-reported Jul 4 `stage_7.0` outputs are pre-repair
  artifacts; they do not reflect the current injector/runtime contract.
  - The master repository does not contain matching MARTINI workflow files or
    MARTINI-specific `py/` / `src/` files at the corresponding paths, so the
    cleanup scope was defined by tracing the four accepted workflows rather
    than by pairing files line-by-line against master.
  - The reported files were written at
    `2026-07-04 02:07` through `02:41` local time.
    Inspection of those saved `stage_7.0.prepared.up` files shows the older
    implicit SC/target contract:
    `cg_lipid_rotamer_sc.arguments = ["placement_fixed_point_vector_only", "compose_vector6d"]`,
    `cg_lipid_target.arguments = ["compose_vector6d", "pos"]`,
    plus `implicit_compaction_gap_response = 1` and `gap_response_coeff`.
  - Direct reinjection with the current
    `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system_lib.py`
    now produces the repaired committed contract:
    `cg_lipid_pair.arguments = ["compose_vector6d", "cgl_compaction_state"]`,
    `cg_lipid_rotamer_sc.arguments = ["placement_fixed_point_vector_only", "compose_vector6d", "cgl_compaction_state"]`,
    and
    `cg_lipid_target.arguments = ["compose_vector6d", "pos", "cgl_compaction_state"]`.
  - Controlled 2,000-step stage-7.1 continuations from the same saved Jul 4
    stage-7 sources do not show the repaired explicit path pushing the proteins
    toward a more vertical state than the implicit-source controls.
    `1rkl` source / implicit / repaired-explicit tilts:
    `46.56 / 45.27 / 46.64 deg`.
    `1afo` chain-A source / implicit / repaired-explicit tilts:
    `18.17 / 19.93 / 19.83 deg`.
    `1afo` chain-B source / implicit / repaired-explicit tilts:
    `16.35 / 14.81 / 15.17 deg`.
  - The installed
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    matches `HEAD` exactly on the inspected pair, compaction, SC, and target
    correction datasets.
  - Interpretation:
    the regression the user observed was tied to pre-repair stage-7 artifacts,
    not to drift in the tracked rigid-DOPC correction tables.
- 2026-07-03: The accepted generic fix is the overlay-center-matched SC
  compaction retrofit on the source-balanced full-bilayer base H5.
  - Accepted H5:
    `/private/tmp/dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_scstateoverlay_scintf_nologin8_20260702.h5`
  - Matched `1rkl` replay:
    `late hbonds 29.08`,
    `late Rg 12.50 A`,
    `late CA RMSD 1.86 A`,
    `axis_trim1_half 0.740`,
    `pca_full_half 0.888`,
    `Dxy_half 0.439`,
    `Dxy_quarter 0.443`.
  - Matched `1afo` replay:
    `late hbonds 83.96`,
    `late Rg 14.83 A`,
    `late CA RMSD 2.68 A`,
    `overall pca_full_half 0.996`,
    `overall pca_trim1_half 0.997`,
    `Dxy_half 0.459`,
    `Dxy_quarter 0.456`.
  - Interpretation:
    the fix removes the remaining source mismatch in the SC state correction,
    stays generic across both proteins, and keeps the hybrid CGL clock close to
    the accepted replay baselines.
- 2026-07-03: The accepted SC retrofit does not change the bilayer-only CGL
  model at all.
  - Direct H5 comparison between the accepted candidate and the installed
    `parameters/dryMARTINI/dopc.h5` shows that
    `cg_lipid_pair/interaction_param`,
    all `cg_lipid_compaction/*` datasets used by the CGL-only model,
    and the CGL orientation/compaction attrs are bitwise identical.
  - The only changed datasets are
    `cg_lipid_sc/delta_extended` and `cg_lipid_sc/delta_compact`.
    The target tables are also unchanged.
  - Interpretation:
    the accepted promotion step preserves the standalone bilayer timescale and
    bilayer mechanics exactly; only the protein side-chain response to CGL
    compaction changes.

- 2026-07-02: The source-balanced full-bilayer branch is viable at the base SC
  table level; the remaining replay regression comes from the exact SC
  compaction-delta retrofit.
  - Exact source-balanced candidate:
    `/private/tmp/dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_20260702.h5`
  - Matched `1rkl` exact replay:
    `late hbonds 32.96`,
    `late Rg 12.48 A`,
    `late CA RMSD 1.83 A`,
    `axis_trim1_half 0.847`,
    `pca_full_half 0.929`,
    `Dxy_half 0.428`,
    `Dxy_quarter 0.430`.
  - Screening proxy on the same new base SC table but with the installed SC
    compaction deltas:
    `/private/tmp/dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_proxybasecomp_20260702.h5`
    with
    `late hbonds 25.87`,
    `late Rg 12.77 A`,
    `late CA RMSD 1.30 A`,
    `axis_trim1_all 0.708`,
    `Dxy_all 0.421`.
  - Interpretation:
    the new full-bilayer/source-balanced base SC `interaction_param` family is
    not the remaining problem.
    The failure sits in how the single-CGL SC compaction deltas are being
    rebuilt after that base table is fit.
- 2026-07-02: The previous exact SC compaction retrofit was mixing sources in
  a biased way.
  - Before the patch:
    `_build_single_cgl_compaction_corrections_from_base_h5()` rebuilt compact
    and extended SC endpoint ensembles from
    `sample_isolated_dopc_bonded_conformers()` plus quartile selection even
    when the base H5 already stored a pooled overlay conformer ensemble and
    physical compaction centers.
  - Code change:
    the retrofit now defaults to selecting compact and extended
    representatives directly from the same overlay reference ensemble, matched
    to `reference_extended_center_ang` and `reference_compact_center_ang`.
    The isolated-MC path remains only as a fallback when those overlay centers
    are unavailable.
  - Interpretation:
    this keeps the explicit SC state correction on the same physical conformer
    family as the rebuilt base SC table and removes the last obvious source
    mismatch in the `1rkl`-failing branch.
- 2026-07-02: The template-based replay harness is the robust path for fresh
  exact stage-7 reinjection tests.
  - A fresh `prepare_stage_file` replay path can fail with
    `cgl_orientation_state rotational_thermostat_timescale must be positive`
    if the prepared file does not carry the rotational thermostat metadata.
  - Reusing a known-good `stage_7.0.prepared.up` template and reinjecting the
    candidate CGL nodes with `inject_cg_lipid_nodes()` avoids that workflow
    issue and keeps the replay comparison low-noise.
  - Interpretation:
    this is a workflow guardrail, not a physics signal.

- 2026-07-02: A generic pooled full-bilayer reference ensemble must be
  balanced across source trajectories, not weighted by raw pooled frame count.
  - Before the patch:
    with `4` pooled full-bilayer representatives from
    `1rkl.stage_7.0 + 1afo.stage_7.0`,
    the reducer picked `3` conformers from `1rkl` and only `1` from `1afo`
    because the `1rkl` trajectory contributed more frame-lipid samples.
    The selected compactions were
    `15.342, 19.215, 21.992, 24.916 A`.
  - Code change:
    `load_dopc_conformers_from_upside_h5_pool()` now allocates representative
    counts per source trajectory with `_balanced_reference_pool_counts()` and
    reduces each source independently before concatenating the selected
    conformers.
  - After the patch:
    the same `4` representatives split evenly as `2/2`
    (`1rkl/1afo`),
    with selected compactions
    `17.143, 17.800, 23.063, 23.497 A`.
  - Interpretation:
    a multi-protein reference ensemble should represent proteins, not just raw
    sample volume.
    Source balancing is a generic physical data-selection rule, not a
    per-protein force tweak.
- 2026-07-02: The unreadable temp H5 on the first source-balanced full-bilayer
  attempt was caused by concurrent writers, not by a new physical-model bug.
  - Observed failure:
    the first exact source-balanced build left only
    `/private/tmp/.dopc_fulltraj_start0_plusiso4x4_scfull_midquant_srcbal_20260702.h5.tmp.*`,
    and that temp file could not be opened by `h5py`.
  - Root cause:
    multiple `martini_gen_params.py` processes were simultaneously targeting
    the same output path.
    Once those stale builds were terminated, a single clean rebuild could be
    launched on the same branch.
  - Interpretation:
    that corruption does not count against the physical validity of the
    source-balanced full-bilayer method itself.
- 2026-07-02: The trajectory-derived shell representative selector is itself a
  likely physical bias source.
  - Current behavior:
    `_select_reference_ensemble_representatives()` chooses representative
    conformers by exact compaction rank endpoints.
    With `12` representatives from the pooled shell-filtered
    `1rkl + 1afo` interface ensemble, the selected compactions are the exact
    `0, 9.1, ..., 90.9, 100` percentiles rather than probability-centered
    quantiles.
  - For the broadened shell pool (`frame_start_fraction=0.0`,
    `XY <= 15 A`, `3D <= 18 A`), that means the selector always injects the
    rarest compaction extremes into the reduced ensemble:
    `4.215 A` at the absolute minimum and `30.316 A` at the absolute maximum.
    The late-half shell pool does the same with `7.114 A` and `29.832 A`.
  - The pooled selector also inherits the raw trajectory count imbalance.
    On the broadened shell pool the chosen `12` representatives split
    `8` from `1rkl` and `4` from `1afo`;
    on the late-half shell pool they split `7` from `1rkl` and `5` from
    `1afo`.
  - Interpretation:
    every shell-heavy retrain so far is being fit against reduced shell
    ensembles that overemphasize rare interface deformations.
    Replacing endpoint selection with probability-centered representatives is a
    physical upstream fix and should be tested before more shell-table rebuilds.
- 2026-07-02: The probability-centered replacement selector immediately removes
  the worst shell-pool artifacts before any new force-field build.
  - Code change:
    `_select_reference_ensemble_representatives()` now targets equal-probability
    bin centers instead of the exact minimum and maximum compaction ranks.
  - On the same broadened shell pool (`frame_start_fraction=0.0`,
    `XY <= 15 A`, `3D <= 18 A`), the reduced `12`-conformer shell set changes
    from
    `4.215, ..., 30.316 A`
    to
    `12.896, ..., 26.126 A`
    and the source split changes from `8/4` (`1rkl/1afo`) to `6/6`.
  - On the late-half shell pool, the reduced set similarly changes from
    `7.114, ..., 29.832 A`
    to
    `12.669, ..., 25.972 A`,
    again with a `6/6` source split.
  - Interpretation:
    the fix is not cosmetic.
    It removes forced shell-compaction outliers and also makes the pooled
    multi-protein reference ensemble much less trajectory-count biased.
- 2026-07-02: The probability-centered selector fix does not by itself rescue
  the shell-plus-isolated SC family.
  - Screening artifact:
    `/private/tmp/dopc_shell15xy3d18_start0_plusiso4_scfull_midquant_screenproxy_20260702.h5`
    paired the new midquant shell-plus-isolated `interaction_param` table with
    the closest available completed shell-plus-isolated compaction-delta set so
    the branch could be screened before spending another long exact-retrofit
    pass.
  - Matched `1rkl` screening result:
    after minimization and burn-in, the production replay was already clearly
    regressing by `5.8` time units.
    The logged protein `Rg` climbed from `12.7 A` at step `0` to
    `13.9 A` by step `2900`, while hbonds fell from `33.8` to about
    `22-24`.
  - Interpretation:
    fixing the representative-endpoint bias is worthwhile, but it does not
    remove the deeper problem that shell-proximal reference conformers still
    encode a bad generic protein-coupling tendency.
    The next physical source ensemble should remove the shell conditioning
    itself rather than keep iterating inside that family.
- 2026-07-02: Broadening the hydrophobic shell rows alone does not produce a
  valid replay.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_start0_hydrophobic_candidate_20260702.h5`
    built by rebuilding only the hydrophobic rows from the broadened
    full-trajectory shell-filtered ensemble.
  - Replay result on matched `1rkl`:
    production quickly became obviously bad, with `Rg` around `13.3-13.5 A`
    and materially reduced hbonds.
  - Interpretation:
    the orientation problem is not fixed by swapping just the hydrophobic shell
    rows, and that branch is not worth full completion.
- 2026-07-02: The broadened all-row shell retrain is also rejected on matched
  `1rkl`.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_start0_allshell_candidate_20260702.h5`
    built from the broadened full-trajectory shell-filtered ensemble for all
    SC rows.
  - Matched replay result versus hydromix baseline on `1rkl`:
    `late hbonds 34.95 vs 30.62`,
    `late Rg 11.68 vs 12.84 A`,
    `Dxy_half 0.464 vs 0.421`,
    `Dxy_quarter 0.461 vs 0.409`,
    but orientation surrogates become much more vertical:
    `axis_ee_3_22_half 0.987 vs 0.749`,
    `axis_pca_trim1_half 0.979 vs 0.840`,
    `axis_pca_3_22_half 0.989 vs 0.839`.
  - Interpretation:
    broadening the shell ensemble helps compactness and keeps the CGL transport
    rate in range, but the shell-only retrain still overcouples the protein to
    the bilayer normal and is not an acceptable generic fix.
- 2026-07-02: Broadening only the nonhydrophobic shell rows to the full
  shell-filtered trajectory does not solve the vertical-drift problem.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_start0_nonhydrophobic_candidate_20260702.h5`
    built by overwriting only
    `ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR`
    on the installed base with rows retrained from the broadened
    `frame_start_fraction=0.0` shell-filtered ensemble.
  - Matched replay result on `1rkl`:
    exact late hbond mean rises strongly
    (`35.81` vs `30.62`),
    production-log late `Rg` stays compact
    (`~12.07 A` vs `~12.84 A`),
    and the CGL transport rate remains acceptable
    (`Dxy_half 0.436` vs `0.421`,
    `Dxy_quarter 0.437` vs `0.409`).
  - But the branch is clearly more vertical, independent of the exact axis
    window:
    full-CA end-to-end `|cos(theta_z)|` rises from `0.841` to `0.942`,
    the previously useful core-window end-to-end surrogate (`res 3..22`)
    rises from `0.749` to `0.932`,
    and full-CA PCA surrogates rise similarly
    (`0.839 -> 0.939` for full-CA PCA,
    `0.840 -> 0.936` for trimmed PCA).
  - Interpretation:
    broadening the shell reference only for the nonhydrophobic rows tightens
    the interface and preserves the timescale, but it does not remove the
    bad orientational bias.
    The remaining dominant bias likely sits in the hydrophobic shell rows that
    were left unchanged from the late-half shell build.
- 2026-07-02: The more isolated-biased pooled-reference retry fails `1rkl`
  more severely than the first pooled candidate and should not be promoted.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_plusiso2x4_scfull_retry_20260702.h5`
    with `2` shell-filtered conformers plus `4` isolated conformers, target
    tables preserved from the current shell-derived base.
  - Runtime path:
    the repaired build completed fully and the replay-prepared stage used the
    explicit local compaction-state path for pair, target, and SC nodes with
    no restored implicit pair-response attrs.
  - Matched replay result versus hydromix baseline on `1rkl`:
    `late hbonds 19.74 vs 30.62`,
    `late Rg 12.38 vs 12.76 A`,
    `late CA RMSD 2.63 vs 1.76 A`,
    `axis_trim1 0.720 vs 0.730`,
    `Dxy_step_half 0.437 vs 0.421`,
    `Dxy_step_quarter 0.439 vs 0.409`.
  - Interpretation:
    increasing the isolated component does not solve the generic vertical-drift
    issue cleanly; it instead strips too much interface support and destabilizes
    the protein replay while leaving lateral transport only slightly faster.
    The pooled isolated route is no longer the preferred physical replacement
    path.
- 2026-07-02: A matched reinjection replay harness now exists for cheap,
  seed-controlled `stage_6.0 -> stage_7.0` candidate gates.
  - Method:
    reuse the successful prepared `stage_7.0` template from a working replay,
    reinject `cg_lipid_pair`, `cg_lipid_target`, and `cg_lipid_rotamer_sc`
    directly with `inject_cg_lipid_nodes()`, hand off the same shared
    `stage_6.0.up` coordinates, and run the same minimization, burn-in, and
    production stages under the hydromix replay seed.
  - Benefit:
    this removes ambiguity from earlier candidate comparisons that mixed
    different replay setups or seeds, so remaining physical candidates can be
    screened on the same trajectory handoff and production random stream.
- 2026-07-02: Rebuilding only the nine nonhydrophobic rows from the pooled
  physical `4 shell + 4 isolated` SC fit does not fix the vertical-bias issue.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_plusiso4x4_nonhydrophobic_candidate_20260702.h5`
    built by overwriting only
    `ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR`
    on the installed hydrophobic-shell base.
  - Matched replay result versus hydromix baseline on `1rkl` with the shared
    `stage_6.0` handoff and seed `20260702`:
    `late hbonds 34.26 vs 30.62`,
    `late Rg 13.30 vs 12.76 A`,
    `late CA RMSD 2.71 vs 1.76 A`,
    `axis_trim1 0.819 vs 0.730`,
    `Dxy_step_half 0.408 vs 0.421`,
    `Dxy_step_quarter 0.411 vs 0.409`.
  - Interpretation:
    softening only the nonhydrophobic rows keeps the protein intact and
    preserves the CGL timescale, but it still leaves too much of the shell-row
    vertical bias in place and overexpands the protein in production.
    The hydrophobic shell rows remain implicated in the bad orientation trend.
- 2026-07-02: A strict same-seed replay rejects the full pooled `4 shell + 4
  isolated` SC candidate as well.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_plusisolated_scfull_20260702.h5`
    tested on the shared reinjection replay harness with the hydromix
    production seed `20260702`.
  - Matched replay result versus hydromix baseline on `1rkl`:
    `late hbonds 28.40 vs 30.62`,
    `late Rg 13.48 vs 12.76 A`,
    `late CA RMSD 2.26 vs 1.76 A`,
    `axis_trim1 0.809 vs 0.730`,
    `Dxy_step_half 0.416 vs 0.421`,
    `Dxy_step_quarter 0.415 vs 0.409`.
  - Interpretation:
    the earlier pooled full-row candidate does not survive a fair replay A/B.
    Its main failure mode is not the CGL timescale but protein expansion and
    lingering vertical-bias coupling in production.
- 2026-07-02: The next physically motivated hypothesis is that the shell
  reference ensemble itself is biased because current shell retrains only use
  the late half of trajectories after the proteins have already drifted toward
  the bilayer normal.
  - Practical implication:
    instead of adding more isolated-DOPC regularization, try rebuilding the
    nonhydrophobic shell rows from a broader shell-filtered ensemble that
    includes earlier trajectory frames by lowering
    `--dopc-target-overlay-reference-frame-start-fraction`.
- 2026-07-02: A single generic SC retrain can now use one pooled physical
  reference ensemble that combines shell-filtered bilayer conformers with
  representative isolated-DOPC conformers.
  - Code change:
    `_resolve_dopc_reference_ensemble()` and
    `build_dopc_target_overlay_h5()` now support appending representative
    isolated conformers to a trajectory-derived overlay reference via the new
    CLI flag
    `--dopc-target-overlay-reference-include-isolated`.
  - Provenance:
    rebuilt overlays record
    `dopc_interface_sc_reference_source =
    pooled_full_resolution_upside_shell_filtered_conformer_trajectory_plus_isolated_dopc_representative_ensemble`
    along with the isolated component count, pool count, and seed.
  - Probe validation:
    a one-row `SER` pooled-reference overlay completed cleanly through the
    single-CGL compaction retrofit, confirming the new path is not just a
    source diff but an executable rebuild route.
- 2026-07-02: The first full pooled-reference SC candidate stayed physical at
  runtime but still failed the matched `1rkl` replay gate.
  - Candidate:
    `/private/tmp/dopc_shell15xy3d18_plusisolated_scfull_20260702.h5`
    with `4` shell-filtered conformers plus `4` isolated conformers, target
    tables preserved from the current shell-derived base.
  - Runtime-path verification:
    after stage preparation completes, the fresh replay-prepared
    `1rkl.stage_7.0.prepared.up` uses
    `cg_lipid_pair.arguments = ["compose_vector6d", "cgl_compaction_state"]`,
    `cg_lipid_target.arguments = ["compose_vector6d", "pos", "cgl_compaction_state"]`,
    and
    `cg_lipid_rotamer_sc.arguments = ["placement_fixed_point_vector_only", "compose_vector6d", "cgl_compaction_state"]`
    with no retained implicit pair-response attrs.
  - Matched replay result versus hydromix baseline on `1rkl`:
    `late hbonds 29.67 vs 30.62`,
    `late Rg 12.63 vs 12.76 A`,
    `late CA RMSD 2.14 vs 1.76 A`,
    `axis_trim1 0.792 vs 0.730`,
    `Dxy_step_half 0.433 vs 0.421`,
    `Dxy_step_quarter 0.440 vs 0.409`.
  - Conclusion:
    the `4+4` pooled SC retrain is still too shell-biased to replace the
    installed mixed-source artifact.
- 2026-07-02: Do not inspect a replay-prepared `.up` file for final node
  wiring until stage preparation/injection has actually finished.
  - During the active `prepare_stage_file()` write, a transient read of
    `1rkl.stage_7.0.prepared.up` showed a half-written `cg_lipid_pair` with
    implicit-response attrs that disappeared once injection completed.
  - Rule:
    wait for the stage-preparation completion line before using a prepared-file
    node audit as evidence about the runtime physics path.
- 2026-07-02: Rebuilding a new overlay from a mixed base H5 must scrub stale
  interface provenance attrs instead of inheriting the old mix story.
  - Real metadata bug:
    `build_dopc_target_overlay_h5()` copied top-level attrs from the base H5
    and then rebuilt SC or target tables, but a full SC rebuild from the
    installed hydrophobic mix still inherited stale attrs such as
    `dopc_interface_sc_reference_source` and
    `dopc_interface_sc_replaced_residue_names`.
  - Code fix:
    `py/martini_build_tables.py` now deletes stale
    `dopc_interface_reference_*`,
    `dopc_interface_sc_*`,
    `interface_reference_*`, and
    `sc_interface_reference_*` attrs before writing overlay provenance, and it
    now writes the SC reference attrs even when SC and target are rebuilt
    together.
- 2026-07-02: The least-ad-hoc route to a completely physical replacement is
  not a fresh 18-row retrain from scratch; it is to finish the current shell
  table family by retraining only the nine nonhydrophobic SC rows from the
  same shell-filtered reference ensemble already used by the installed target
  tables and preserved hydrophobic SC rows.
  - Physical logic:
    the installed `parameters/dryMARTINI/dopc.h5` already carries shell-derived
    target tables and shell-derived hydrophobic SC rows.
    Rebuilding `ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR` from the same
    pooled `1rkl + 1afo` shell reference would make every SC row shell-derived
    without introducing a new mixed source.
  - Practical blocker:
    both the full 18-row shell rebuild and the cheaper nine-row selective shell
    rebuild remained expensive in this environment.
    The selective build reached the intended
    `Parallel CG-SC overlay table: 9 tasks on 8 worker(s)` path and could be
    moved outside the sandbox to avoid the thread fallback, but it still did
    not finish within this turn, so no validated all-shell H5 was promoted.
- 2026-07-02: The current installed hydrophobic-shell SC mix was physically
  usable but not reproducible from the repo until the row-merge path was wired.
  - Real implementation gap:
    `py/martini_build_tables.py` already contained
    `_overwrite_sc_rows_from_source_group()`, but nothing called it.
    `py/martini_gen_params.py` also did not expose the shell-filter distance
    options or residue-subset SC rebuild option needed to recreate the installed
    `dopc.h5`.
  - Code fix:
    `build_dopc_target_overlay_h5()` now accepts merge-only
    `sc_row_source_h5_path`, so a final SC row-mixed overlay can be built
    without forcing a target or SC rebuild in the last step.
    `martini_gen_params.py` now exposes:
    `--dopc-target-overlay-reference-max-protein-xy-distance-ang`,
    `--dopc-target-overlay-reference-max-protein-distance-ang`,
    `--dopc-target-overlay-rebuild-sc-residue`, and
    `--dopc-target-overlay-sc-row-source-h5`.
  - Cheap end-to-end validation:
    a merge-only CLI run built `/tmp/dopc_selfmerge_20260702.h5` successfully,
    and checked SC/target datasets matched the installed
    `parameters/dryMARTINI/dopc.h5` exactly (`max_abs_diff = 0` on inspected
    tables).
    A selective `SER` subset CLI run also reached the shell-lattice SC rebuild
    stage and then the single-CGL compaction-correction stage before being
    stopped for runtime cost.
  - Runtime validation:
    the live default hybrid checkpoints under
    `example/16.MARTINI/outputs/martini_1rkl_hybrid` and
    `.../martini_1afo_hybrid` are using the installed shell-mix table
    signature (`cg_lipid_rotamer_sc n_radial = 25`,
    `cg_lipid_target n_radial = 91`) with explicit
    `cgl_compaction_state` arguments for SC and target nodes.
    Their late-half summaries remain acceptable:
    `1rkl hbonds ~= 31.51`, `Rg ~= 12.89 A`, `Dxy ~= 0.287 / 0.284`;
    `1afo hbonds ~= 80.43`, `Rg ~= 15.63 A`,
    `Dxy ~= 0.311 / 0.322`.
- 2026-07-01: Removing the nonphysical pair-side gap-response shortcut does
  not slow the bilayer; it shifts the current matched-seed `1rkl` workflow
  faster than the retained explicit branch.
  - Code result:
    `example/16.MARTINI/run_sim_hybrid.sh` now defaults
    `UPSIDE_CGL_COMPACTION_IMPLICIT_RESPONSE=0`, and
    `py/martini_prepare_system_lib.py` no longer copies gap-response attrs or
    `gap_response_coeff` into explicit `cg_lipid_pair` nodes.
  - Prepared-file verification:
    reinjecting a stage-7 file with the new code gives
    `cg_lipid_pair.arguments = ["compose_vector6d", "cgl_compaction_state"]`
    with no `implicit_compaction_response`,
    no `implicit_compaction_gap_response`,
    no `gap_response_coord_*` attrs, and no `gap_response_coeff` dataset.
  - Validation run:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/timescale_seeded_20260701_1rkl_pairexplicit_match`
    finishes cleanly and keeps the protein stable enough
    (`late hbonds ~= 33.25`, `late Rg ~= 12.88 A`).
  - Transport result:
    on the same short-lag CGL estimator used for the rerun comparison,
    the physical pair-explicit run is faster than the retained current
    explicit branch:
    `Dxy_step_half / quarter ~= 0.435 / 0.439`
    versus `0.304 / 0.302` for the matched current explicit run and
    `0.700 / 0.677` for the old too-fast implicit-response control.
  - Interpretation:
    the user’s slowdown concern is not supported.
    The remaining problem is now the opposite direction:
    once the nonphysical runtime shortcut is removed, the installed explicit
    pair-state model is not yet calibrated to preserve the accepted bilayer
    clock.
- 2026-07-01: The current branch is acceptable on the project's diffusion gate,
  but it is not fully "physical" under the strictest interpretation of the
  project's modeling rules.
  - Good:
    the active prepared stage has no membrane lateral potential and no
    separate external alignment field.
    The explicit `cg_lipid_target` and `cg_lipid_rotamer_sc` paths depend only
    on local CGL-protein geometry plus each lipid's own compaction state.
  - Real caveat:
    `cg_lipid_pair` still runs through the implicit gap-response layer.
    In `src/martini_cg_lipid.cpp`, that layer first classifies leaflets by the
    median `z` coordinate (`classify_leaflets_by_median_z`) and then builds a
    cross-leaflet gap field from those frame-dependent labels.
    That means the pair-side compaction response is not a purely local,
    bilayer-agnostic Hamiltonian.
  - Real caveat:
    the production workflow still relies on calibrated auxiliary dynamics:
    `input/cgl_gle`, `cgl_orientation_state.rotational_thermostat_timescale`,
    and `cgl_compaction_state.thermostat_timescale`.
    These are coarse-grained transport / hidden-state controls, not equilibrium
    interface potentials, but they are still phenomenological rather than
    fully first-principles dynamics.
  - Minor caveat:
    the prepared SC explicit-state metadata still mixes normalized compaction
    coordinates (`cgl_compaction_state` on `0..1`) with physical state-center
    attrs copied from the SC table (`13.38 / 19.62 A`).
    The current SC runtime does not use those attrs, so this is not an active
    force bug, but it leaves the model self-contradictory.
- 2026-07-01: The apparent severe CGL slowdown was mostly a measurement bug,
  not a real failure of the current explicit branch.
  - The diffusion analysis must read `output/pos` as
    `(n_frame, n_system, n_atom, 3)`.
    The earlier much smaller `Dxy` values came from mishandling that shape.
  - The retained project target in
    `example/16.MARTINI/cg_lipid_potentials.tex` is about
    `0.25 A^2/tu`.
  - Recomputed current-branch reruns are close to that target:
    - matched-seed explicit rerun
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/timescale_seeded_20260701_1rkl_explicit_match`
      gives `Dxy_half=0.261`, `Dxy_quarter=0.211`
    - earlier explicit fresh rerun
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/timescale_fresh_20260701_1rkl_current40`
      gives `Dxy_half=0.322`, `Dxy_quarter=0.259`
  - The old implicit-response path is much faster than the accepted clock.
    The current-table implicit control at
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/timescale_fresh_20260701_1rkl_implicit40d`
    gives `Dxy_half=1.462`, `Dxy_quarter=1.364`.
    Its late helix-axis `|cos(theta_z)|` is also more vertical
    (`~0.961` half-mean) than the current explicit rerun
    (`~0.886` half-mean).
  - A matched-seed SC-implicit-only split is also faster than the retained
    explicit branch and slightly weaker structurally.
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/timescale_seeded_20260701_1rkl_scimplicit_match`
    gives `Dxy_half=0.300`, `Dxy_quarter=0.291`,
    late hbond sum about `30.71` versus `31.65` for the matched explicit run,
    and late `Rg` about `13.26 A` versus `12.82 A`.
  - Conclusion:
    the current explicit compaction branch is not too slow relative to the
    project timescale target.
    Reverting SC or target to the older implicit-response default would move
    the CGL clock in the wrong direction.
- 2026-06-30: The current uncommitted branch combines two pieces of the same
  intended fix.
  - `py/martini_build_tables.py` adds support for training interface-shell
    DOPC reference ensembles by filtering full-resolution lipids with a
    maximum protein XY distance.
  - `src/martini_cg_lipid.cpp` adds optional
    `cgl_compaction_state` support to `cg_lipid_target`.
- 2026-06-30: The missing wiring was in
  `py/martini_prepare_system_lib.py`.
  - Archived hybrid prepared files still showed
    `cg_lipid_target.arguments = ["compose_vector6d", "pos"]`.
  - No `cgl_compaction_state` node existed in the archived
    `martini_1rkl_hybrid` or `martini_1afo_hybrid` stage-7 checkpoints.
  - Result:
    the previous reruns did not test the new explicit expandable-tail path,
    even though the H5 carried target compaction correction datasets.
- 2026-06-30: The previous agent already trained an interface-target overlay
  artifact at
  `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_interfacefull15a_t12_overlay/dopc_target_interfacefull15a_t12_overlay.h5`.
  - It preserves pair/SC tables and rebuilds only the target-facing interface
    tables.
  - Its interface reference source is the pooled full-resolution
    `1rkl` + `1afo` stage-7 lipid ensemble with `XY <= 15 A` to the protein.
- 2026-06-30: The prior validation directory for that overlay,
  `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_interfacefull15a_t12/`,
  contains burn-in output only and no finished production log.
  It is not a valid acceptance run for the user request.
- 2026-06-30: Wiring `cg_lipid_target` to an explicit `cgl_compaction_state`
  requires both sides of the interface fix.
  - The injector must create `cgl_compaction_state`,
    `cg_lipid_compaction_self`, and target arguments
    `["compose_vector6d", "pos", "cgl_compaction_state"]`.
  - The C++ target node must accept the optional third argument, use it in
    `compute_value`, and propagate sensitivity back to the compaction state in
    `propagate_deriv`.
- 2026-06-30: Fresh validation after completing that wiring gives a partial
  improvement, not a full acceptance.
  - New `1rkl` run:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_interfacefull15a_t12_explicit`
    reduces late chain-A axis alignment from about `0.926` to `0.837`.
  - New `1afo` run:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_target_interfacefull15a_t12_explicit`
    does not resolve the issue consistently; chain B remains almost fully
    aligned with the CGL axis field late in the trajectory.
- 2026-06-30: The fully retrained explicit-SC overlay exposed a second
  explicit-state runtime bug in preparation.
  - The overlay compaction self potential is defined on a normalized
    `0..1` hidden coordinate.
  - Fresh prepared files were still writing
    `cgl_compaction_state/value ~= 17.6..31.1` from
    `compose_vector6d/initial_compaction_ang`.
  - The C++ state node clamps to its declared domain on load, so the intended
    startup state was being collapsed onto the compact boundary.
- 2026-06-30: Mapping the injected startup state into the overlay coordinate
  system fixes that inconsistency materially.
  - New prepared files now show
    `cgl_compaction_state/value ~= 0.676..1.000` for both `1rkl` and `1afo`
    under the full retrained overlay.
  - Fresh reruns:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_targetqfix_fullsc_explicit_initmap`
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_targetqfix_fullsc_explicit_initmap`
    recover most of the lost `1rkl` gain:
    late chain-A `|cos(theta_z)|` improves from about `0.792` to `0.709`,
    late CA RMSD from about `2.67 A` to `2.11 A`.
  - `1afo` does not become a clean win under the corrected runtime:
    chain-A late RMSD worsens from about `1.92 A` to `3.13 A`, and chain-B
    late `|cos(theta_z)|` drifts from about `0.968` to `0.973`.
  - Conclusion:
    the remaining problem is no longer the explicit-state wiring.
    It is now concentrated in the retrained SC correction content or its
    reference ensemble, especially for the `1afo` interface case.
- 2026-07-01: A true local 3D protein shell improves the over-alignment bug,
  but only the tighter `3D <= 18 A` branch is worthwhile.
  - Full retrained local shell:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_targetqfix_fullsc_shell3d18`
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_targetqfix_fullsc_shell3d18`
    improve alignment on both systems at once.
    `1rkl` late `|cos(theta_z)|` drops to about `0.674`, and `1afo` chain B
    drops to about `0.916`.
  - Loosening the shell to `3D <= 20 A` is not an improvement.
    `1rkl` worsens to about `0.691 / 2.68 A`, and `1afo` chain B worsens to
    about `0.950 / 2.65 A`, even though the branch runs softer.
  - Conclusion:
    the local-shell direction is correct, but `3D <= 20 A` should be rejected.
- 2026-07-01: The remaining `1afo` chain-B tradeoff is now isolated to the
  retrained SC layer rather than the target tables.
  - Target-only local shell:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_targetonly_shell3d18`
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_targetonly_shell3d18`
    keep the `3D <= 18 A` target tables while preserving the older SC tables.
  - Result:
    `1rkl` improves further (`|cos(theta_z)|` about `0.601`,
    late RMSD about `2.19 A`) and `1afo` chain A improves strongly
    (`|cos(theta_z)|` about `0.857`, late RMSD about `2.14 A`),
    but `1afo` chain B returns to near-vertical
    (`|cos(theta_z)|` about `0.983`, late RMSD about `2.63 A`).
  - Conclusion:
    the local target shell should be retained.
    The remaining chain-B alignment gain in `fullsc_shell3d18` is coming from
    the retrained SC tables, and so is a substantial part of the remaining
    structural penalty.
- 2026-07-01: Residue-selective SC row support is implemented, but the current
  SC fitting path is still the practical bottleneck.
  - `py/martini_build_tables.py` now has helper support to preserve or
    overwrite selected `cg_lipid_sc` residue rows when composing mixed H5
    overlays.
  - A hydrophobic-first selective candidate
    (`CYS, ILE, LEU, MET, PHE, PRO, TRP, TYR, VAL`) compiles cleanly at the
    code level, but both the direct selective rebuild and a reduced
    selected-residue-only fit path remained long-running enough that this turn
    did not finish a fresh candidate H5 for reruns.
- 2026-07-01: The current SC-compaction retrofit path is numerically
  inconsistent with the runtime tables it produces.
  - The new retrofit builder in
    `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
    reconstructs a reusable SC base grid from existing `interaction_param`
    controls (`_reconstruct_sc_base_result_from_table`) and then fits
    `delta_extended` / `delta_compact` against newly sampled state grids.
  - A direct reconstruction check on
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    shows that evaluating `interaction_param + delta_*` does not reproduce the
    stored SC state grids for several residues, especially `LYS`, `GLU`,
    `VAL`, `ASP`, `ILE`, and `THR`.
    Compact-state mismatches are especially large
    (for example `LYS` compact RMS about `5214 kJ/mol`,
    max error about `165287 kJ/mol`).
  - The runtime only consumes the fitted `interaction_param` and `delta_*`
    controls, not the stored raw `grid_*_kj_mol` datasets.
    Conclusion:
    the SC explicit-state issue is more likely in table construction /
    retrofit consistency than in the current `cg_lipid_rotamer_sc`
    interpolation loop itself.
- 2026-07-01: The prepared SC explicit-state metadata still mixes normalized
  and physical compaction coordinate conventions.
  - In the current committed `dopc.h5`, the compaction self potential is on
    `0..1` with `cg_lipid_compaction` centers `0.0 / 1.0`.
  - The SC table still carries physical state-center attrs
    `13.381 / 19.615 A`.
  - `copy_single_cgl_compaction()` currently prefers the SC-table attrs when
    populating `cg_lipid_rotamer_sc`, so the prepared node carries physical
    center metadata even when the actual hidden state is normalized.
  - The present SC runtime does not use those attrs, so this is not the main
    numerical fault, but it leaves the prepared model self-contradictory and
    mirrors the same kind of scale mismatch that previously broke
    `cg_lipid_target`.
- 2026-07-01: The SC selective-rebuild path needed one more code fix to be
  practically usable.
  - Rebuilding selected SC rows against a different reference ensemble was
    silently changing the SC spline lattice (`n_radial`, `fit_r_max_nm`),
    which made row-level splicing back into the current shell overlay invalid.
  - `py/martini_build_tables.py` now preserves the base overlay SC support
    when a selective rebuild is sourced from an existing H5, so subset builds
    can be merged safely with `_overwrite_sc_rows_from_source_group()`.
- 2026-07-01: The best current compromise is a hydrophobic-only shell SC mix.
  - Built subset source:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_interfacefull15a_t12_overlay/dopc_target_interfacefull15a_t12_overlay.shell3d18_sc_nonhydrophobic_base_subset.h5`
    containing isolated-ensemble replacements for
    `ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR` on the preserved shell-18
    lattice.
  - Final mixed candidate:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_interfacefull15a_t12_overlay/dopc_target_interfacefull15a_t12_overlay.shell3d18_sc_hydrophobic_mix.h5`
    keeps local-shell rows for `CYS, ILE, LEU, MET, PHE, PRO, TRP, TYR, VAL`
    and overwrites the non-hydrophobic rows from the subset source.
  - Validation:
    `1rkl` late `|cos(theta_z)|` improves from the archived `0.799` to
    `0.722` and late CA RMSD from `2.69 A` to `1.54 A`.
    `1afo` chain A improves strongly (`0.971 / 3.44 A -> 0.839 / 1.64 A`).
    `1afo` chain B reduces verticalization relative to the archived hybrid
    (`0.965 -> 0.938`) while keeping RMSD moderate (`1.51 A -> 1.65 A`) and
    greatly improving over the broad shell branch RMSD (`2.21 A -> 1.65 A`).
- 2026-07-01: The validated candidate had not yet been installed into the
  default workflow path.
  - `example/16.MARTINI/run_sim_hybrid.sh` defaults to
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    unless `DOPC_H5` is set.
  - The fresh validation runs were performed with
    `DOPC_H5=/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_interfacefull15a_t12_overlay/dopc_target_interfacefull15a_t12_overlay.shell3d18_sc_hydrophobic_mix.h5`.
  - The installed
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    still had its older `2026-06-29` timestamp before promotion, which is why
    the user did not see an updated `.h5` in the default location.
- 2026-07-01: The validated candidate is now installed at the default path.
  - Backup of the pre-promotion file:
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/backup/install_20260701_132601/dopc.h5`
  - Installed file:
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  - Verification:
    the installed file matches the validated candidate
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_interfacefull15a_t12_overlay/dopc_target_interfacefull15a_t12_overlay.shell3d18_sc_hydrophobic_mix.h5`
    byte-for-byte (`cmp`) with identical SHA-256
    `e8d09c35ed37a9f65bcb491c32a7dce9a200a049944c5761205d5563c0c318c8`.
- 2026-07-01: The default stock workflows now consume the installed candidate
  cleanly with no `DOPC_H5` override.
  - Fresh stock-path validation directories:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260701_1rkl_cleanstart_stockdefault_hydromix`
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260701_1afo_cleanstart_stockdefault_hydromix`
  - The preparation logs for both runs explicitly use the force-field
    directory `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI` and
    inject `cg_lipid_pair`, `cg_lipid_rotamer_sc`, and `cg_lipid_target` into
    the prepared stage-7 files.
  - Late metrics from these fresh runs differ from the earlier
    override-backed hydrophobic-mix validations because the stochastic seeds
    differ.
    Stock-path rerun late metrics:
    `1rkl A 0.732 / 2.72 A`,
    `1afo A 0.904 / 1.93 A`,
    `1afo B 0.958 / 2.33 A`.
    Earlier override-backed hydrophobic-mix metrics:
    `1rkl A 0.722 / 1.54 A`,
    `1afo A 0.839 / 1.64 A`,
    `1afo B 0.938 / 1.65 A`.
  - Conclusion:
    the installation and stock-path injection are verified; exact behavioral
    parity would require rerunning with matched seeds or a reused packed
    starting structure.
