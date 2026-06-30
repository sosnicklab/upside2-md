# Progress Log

## Current Task: Universal Coarse-Hybrid Protein Stability

- Actions taken:
  - Reconstructed protein-only DSSP and RMSD directly from the cited hybrid
    checkpoints.
  - Verified that `1rkl stage_6.0` is structurally stable and that the
    secondary-structure loss begins only after release into unrestrained
    `stage_7.0` production.
  - Replayed stage-7 production from the promoted burn-in input to separate
    activation-path defects from conservative bilayer or transport defects.
  - Tested three activation-fix families:
    - short restrained release after burn-in: rejected
    - global post-burn-in hybrid transition reset / hold: rejected as a
      universal default because it hurts `1afo`
    - generic opt-in post-burn-in release-hold support: kept only as
      diagnostic instrumentation, not as a retained production fix
  - Removed the active `1rkl`-only launcher override after the user required
    universal simulation settings across proteins.
  - Localized the remaining coarse defect to the protein-lipid interface:
    coarse `1rkl` loses its native tilted transmembrane pose and frays mainly
    at the helix caps, while coarse `1afo` and the full-lipid `1rkl`
    reference remain materially more stable.
  - Confirmed that the shared CGL target force field is underresolved at build
    time: the active `dopc.h5` was built from only two isolated DOPC
    conformers, and direct ensemble probes show large seed sensitivity in the
    protein-facing `cg_lipid_target` support and PMFs.
  - Implemented a shared table-builder redesign to replace the old
    two-conformer default with a deterministic representative ensemble:
    sample a larger isolated-DOPC pool and select compaction-stratified
    representatives for the default CGL base tables.
  - Started a full consistent `dopc.h5` rebuild so pair, SC, and target
    tables all use the same patched source path.
  - Aborted that first full rebuild before any write phase after confirming it
    was spending the wall clock in the already-accepted `cg_lipid_pair` tensor
    rather than on the protein-facing defect surface.
  - Rebuilt only `cg_lipid_target` on top of the accepted `dopc.h5` base using
    the new `8`-representative / `32`-pool isolated-DOPC source path.
  - Injected that rebuilt target node into a short coarse `1rkl stage_7.0`
    handoff smoke and confirmed that the runtime consumed the new
    `84 x 321` target tensor without error.
  - Ran the matching short coarse `1afo stage_7.0` handoff smoke with the same
    shared settings and confirmed the rebuilt target node also injects and runs
    cleanly there.
  - Restored `parameters/dryMARTINI/dopc.h5` to the last clean consistent
    baseline after the two diagnostic smokes and removed the extra temporary
    backup H5 files created during this iteration.
  - Confirmed that a full representative `cg_lipid_sc` rebuild is still too
    expensive to use as the next retained branch.
  - Built a cheaper representative-target overlay H5 that preserves the
    accepted pair/SC/compaction tables from the clean baseline and rebuilds
    only `effective_lj` plus `cg_lipid_target` from the `32 -> 8`
    isolated-DOPC representative ensemble.
  - Added that overlay builder to source as
    `build_dopc_target_overlay_h5(...)` plus CLI support in
    `py/martini_gen_params.py`.
  - Identified an internal inconsistency in the first representative-SC overlay
    path before trusting it for `1rkl` validation:
    the rebuilt interface PMFs still inherited compact/extended correction
    tensors from the old baseline `dopc.h5`, so the branch did not actually
    replace the stale tail-compaction response on the protein-facing SC/target
    surface.
  - Patched the overlay builder to regenerate interface compaction corrections
    from the rebuilt representative ensemble instead of copying baseline
    target/SC deltas forward, and removed the now-dead baseline-delta copy
    helper from `py/martini_build_tables.py`.
  - Aborted the first SC-overlay H5 generation after confirming it was still
    using the inconsistent correction path, then restarted a fresh build on
    the corrected compaction-consistent branch.
  - Validated the overlay branch in longer same-settings coarse continuations
    for both `1rkl` and `1afo`, then measured protein-only and bilayer metrics
    directly from the resulting `stage_7.1` trajectories.
  - Resolved the clean-start workflow path for retained-branch acceptance:
    use a fresh `RUN_DIR` so `run_sim_hybrid.sh` cannot auto-detect an old
    `stage_7.*` continuation source.
  - Ran a true clean-start retained-branch `1afo` workflow in
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_cleanstart_target_rep8/`
    with the representative-target overlay H5 and explicit seeds.
  - Measured the finished clean-start `1afo` production directly from the
    stage file and confirmed a correct bilayer plus normal production kinetics,
    while noting that the protein keeps strong secondary structure but explores
    a looser aligned pose than the archived coarse baseline.
  - Searched the live repo, master repo, and backup repos for the requested
    `elk` target and found no matching input, output, or launcher artifact.
  - Removed the aborted SC-representative overlay H5 and the aborted
    direct-from-prepared production probe outputs after they proved
    non-comparable.
  - Implemented a new bilayer-derived CGL interface-reference path in the
    overlay builder:
    - pool one or more full-resolution Upside trajectories
    - exclude lipids within a configurable XY distance of non-lipid particles
    - preserve the trajectory/filter provenance in the overlay H5 metadata
  - Identified and fixed a new source bug in that pooled-bilayer path:
    the first implementation selected retained conformers by evenly spaced
    pooled frame index instead of by the compaction coordinate, which biased
    the retained ensemble toward an over-compact subset.
    The builder now reuses the deterministic compaction-stratified
    representative selector for the pooled bilayer ensemble too.
  - Built and validated two bilayer-derived target-only overlay branches for
    clean-start `1rkl`, and rejected both.
    - first `10 A` XY-filtered pooled-bilayer branch before the selector fix
    - corrected representative pooled-bilayer branch after the selector fix
  - Aborted a full bilayer-derived `cg_lipid_target + cg_lipid_sc` overlay
    build after confirming that the SC worker pool was too slow to use as the
    first practical falsification branch.
  - Patched the overlay builder so it can preserve the current best target
    term while rebuilding only `cg_lipid_sc` from a bilayer-derived reference
    ensemble.
    The new source path keeps the base `effective_lj` plus
    `cg_lipid_target`, writes an SC-specific reference dataset, and applies an
    SC-only compaction retrofit.
  - Built the new SC-only bilayer-derived overlay at
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_sc_bulkbilayer10a_t12_rep_overlay/dopc_target_rep8_sc_bulkbilayer10a_t12_rep_overlay.h5`
    from the pooled `1rkl` + `1afo` full-bilayer reference ensemble with the
    same `10 A` XY protein-distance filter and `12` retained representative
    conformers.
  - Ran a fresh clean-start `1rkl` workflow on that SC-only overlay and
    re-extracted the retained acceptance metrics directly from the finished
    stage file.
  - Cross-checked the same overlay on a fresh clean-start `1afo` workflow to
    verify that the `1rkl` improvement did not come from a benchmark-specific
    regression.
  - Reconstructed the exact clean-start metric definitions directly from the
    stage files:
    - kinetic ratio from the final `stage_7.0.log` summary line
    - orientation flips from sign reversals relative to the initial CGL
      production orientation
    - `1rkl` helix-core occupancy from simplified DSSP on residue slice `3:28`
      (0-based), matching the retained baseline

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_hybrid.sh`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_1rkl.sh`
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_gen_params.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`

- New output artifacts:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_ss_probe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_release_probe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_transitionreset_probe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_transitionhold_probe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_transitionhold_helperprobe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage7_transitionhold_helperprobe/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage7_currentreplay/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage7_transitionhold1000/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage7_transitionhold1000/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_cli/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_bulkbilayer10a_t12_overlay/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_bulkbilayer10a_t12_rep_overlay/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_sc_bulkbilayer10a_t12_rep_overlay/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage70_target_rep8_long/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage70_target_rep8_long/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_cleanstart_target_rep8/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_bulkbilayer10a_t12/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_bulkbilayer10a_t12_rep/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_rep8_sc_bulkbilayer10a_t12_rep/`
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_target_rep8_sc_bulkbilayer10a_t12_rep/`

- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py`
  - `python3 -m py_compile py/martini_build_tables.py py/martini_gen_params.py`
  - `bash -n example/16.MARTINI/run_sim_hybrid.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `1rkl` cited artifact diagnosis:
    - `stage_6.0` DSSP is unchanged from start to finish
    - cited `stage_7.0` input remains close to stage 6, but unrestrained
      production drifts to about `0.317 nm` protein RMSD in the archived run
  - retained `1rkl` helper replay with the opt-in `2000/1500` release hold:
    - final DSSP match to the promoted burn-in input: `0.839`
    - final protein RMSD to the promoted burn-in input: `0.227 nm`
    - coarse bilayer late-half separation remains about `16.15 A`
    - coarse bilayer late-half `6x6` occupancy CV is about `0.147`
  - `1afo` cross-check:
    - the same release hold is not safe as a repo-wide default
    - a matched baseline replay without the hold keeps `1afo` at
      `DSSP_match≈0.986`, while the hold degrades that benchmark
  - coarse-vs-reference interpretation:
    - coarse `1rkl` loses a large part of its initial tilt during stage 7,
      while the full-lipid `1rkl` reference largely preserves it
    - the helical loss is concentrated near the membrane-entry and
      membrane-exit caps rather than in the buried helix core
  - representative-ensemble force-field probe:
    - four nearby seeds with the old two-conformer target design give
      `n_radial=76..84` and coarse non-core target-PMF variability of about
      `2.66e5 kJ/mol` RMS for `N0`
    - the redesigned representative-ensemble path drops that same non-core
      `N0` variability to about `1.24e3 kJ/mol` in a direct `32 -> 8`
      representative probe
  - short `1rkl` target-rebuild smoke:
    - output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage70_target_rep8_smoke/`
    - injected `cg_lipid_target` node uses `n_radial=84`,
      `n_angular=321`, `cutoff=28.7 A`
    - minimization plus `500` burn-in / `200` production steps completed
      successfully under the shared coarse settings
    - production log shows temperature mean `0.8647` and hbond count
      `31.0 -> 32.7` over the short smoke
  - short `1afo` target-rebuild smoke:
    - output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage70_target_rep8_smoke/`
    - minimization plus `500` burn-in / `200` production steps completed
      successfully under the same shared coarse settings
    - production log shows temperature mean `0.8647` and hbond count in the
      `90.9..93.0` range over the short smoke
  - H5 cleanup:
    - restored live `dopc.h5` to a consistent state with
      `pair/sc/target conformer_count=2/2/2`
    - removed:
      `dopc.h5.bak_pre_refensemblefix_20260629`,
      `dopc.h5.bak_pre_target_rep8_20260629`,
      and `dopc.h5.bak_target_bsplinefix_20260629`
  - representative-target overlay builder:
    - `python3 -m py_compile py/martini_build_tables.py py/martini_gen_params.py`
    - `python3 py/martini_gen_params.py --upside-home /Users/yinhan/Documents/upside2-md --dopc-target-overlay-base-h5 /Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5 --dopc-target-overlay-output /Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_cli/dopc_target_rep8_cli.h5`
    - resulting CLI-built overlay keeps `pair/sc/target conformer_count=2/2/8`
      and preserves target compaction deltas
  - representative-SC overlay correction-path patch:
    - `python3 -m py_compile py/martini_build_tables.py py/martini_gen_params.py`
    - confirmed that the corrected builder now regenerates interface
      compaction corrections after rebuilding the representative SC/target PMFs
      instead of copying the old baseline `delta_extended` / `delta_compact`
      tensors forward
    - the new corrected overlay build is running at
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_sc_compact_overlay/`
      and must complete before the next `1rkl` validation step
  - longer representative-target benchmark validation:
    - `1rkl` continuation output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1rkl_stage70_target_rep8_long/`
    - `1afo` continuation output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_stage70_target_rep8_long/`
    - matched `5..10` time-window metrics:
      `1rkl` target-overlay vs baseline coarse:
      `RMSD_t10=2.10 A vs 1.92 A`, `hbond_t10=32.94 vs 26.15`,
      `sep_5_10=16.17 A vs 16.81 A`, `|n_z|=0.9943 vs 0.9946`,
      `occ_cv=0.160 vs 0.160`, zero flips, and same-leaflet `nn_p05=7.07 A vs 6.53 A`
      `1afo` target-overlay vs baseline coarse:
      `RMSD_t10=2.59 A vs 2.99 A`, `hbond_t10=84.11 vs 80.76`,
      `sep_5_10=15.76 A vs 16.47 A`, `|n_z|=0.9943 vs 0.9944`,
      `occ_cv=0.341 vs 0.278`, zero flips, and same-leaflet `nn_p05=7.04 A vs 6.62 A`
    - continuation caveat:
      those longer comparison runs came from 200-step smoke outputs with
      manually marked restart-valid final states, so they are suitable for
      protein/bilayer comparison but not for transport acceptance
  - clean-start `1afo` retained-branch validation:
    - output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_cleanstart_target_rep8/`
    - workflow:
      full fresh `stage_6.0 -> stage_7.0` handoff under the shared coarse
      settings with the representative-target overlay H5, not a continuation
    - production kinetic ratio:
      `avg_kinetic_energy/1.5kT≈1.053`
    - bilayer metrics:
      `sep_initial≈16.55 A`, `sep_late≈16.70 A`, `sep_max≈16.92 A`,
      `|n_z|_late≈0.99335`, `occ_cv_late≈0.332`,
      `nn_mean/p05≈7.61/6.74 A`, zero flips
    - protein metrics:
      `DSSP_final≈0.986`, `DSSP_late≈0.970`,
      `hbond_final/late≈77.4/81.5`
      aligned CA RMSD to the stage reference is looser than the archived
      coarse baseline (`final≈5.14 A` vs `≈3.52 A`), but the secondary
      structure and hbond signal remain strong
  - unresolved second benchmark:
    - no local `elk` / `1elk` target exists in the live repo, master repo,
      or backups, so that clean-start validation is blocked pending the exact
      user-supplied identifier
  - bilayer-derived target-only overlay builder:
    - `python3 -m py_compile py/martini_build_tables.py py/martini_gen_params.py`
    - first pooled-bilayer target-only overlay:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_bulkbilayer10a_t12_overlay/dopc_target_bulkbilayer10a_t12_overlay.h5`
      built from two full-resolution trajectories with a `10 A` non-lipid XY
      cutoff and `12` retained conformers
    - selector-fix verification:
      the same pooled-bilayer source now retains a representative conformer
      set with `mean tail extension≈20.95 A` and `median≈21.18 A`, instead of
      the earlier biased compact subset at `mean≈18.11 A`
    - corrected representative target-only overlay:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_bulkbilayer10a_t12_rep_overlay/dopc_target_bulkbilayer10a_t12_rep_overlay.h5`
  - bilayer-derived `1rkl` clean-start validation:
    - first pooled-bilayer branch output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_bulkbilayer10a_t12/`
    - corrected representative branch output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_bulkbilayer10a_t12_rep/`
    - both branches keep an acceptable bilayer
      (`sep_late≈16.34..16.37 A`, `|n_z|≈0.9935`, zero flips),
      but both are rejected as protein fixes:
      first branch:
      `avg_kinetic/1.5kT≈1.479`, `RMSD_final≈3.27 A`,
      `RMSD_late≈2.64 A`, `hbond_late≈28.41`
      corrected representative branch:
      `avg_kinetic/1.5kT≈1.615`, `RMSD_final≈3.34 A`,
      `RMSD_late≈2.78 A`, `helix_late≈0.622`,
      `helix_core_late≈0.752`, `hbond_late≈22.90`
  - SC-only bilayer-derived overlay build:
    - direct source build:
      `build_dopc_target_overlay_h5(..., rebuild_sc=True, rebuild_target=False, ...)`
      on top of
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_cli/dopc_target_rep8_cli.h5`
    - resulting overlay metadata:
      preserves representative target/effective tables
      (`conformer_count=8`), adds
      `sc_interface_ref_bead_positions_nm`, and rebuilds
      `cg_lipid_sc` from the pooled bilayer reference
      (`conformer_count=12`) with an SC-only compaction retrofit
  - fresh `1rkl` validation on the SC-only bilayer overlay:
    - output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1rkl_cleanstart_target_rep8_sc_bulkbilayer10a_t12_rep/`
    - production kinetic ratio:
      `avg_kinetic_energy/1.5kT≈1.262`
    - bilayer metrics:
      `sep_late≈16.11 A`, `sep_max≈16.25 A`,
      `|n_z|_late≈0.99333`, `occ_cv_late≈0.226`,
      `nn_mean/p05≈7.44/6.44 A`, zero flips
    - protein metrics:
      `CA_RMSD_final/late≈2.79/2.05 A`,
      `helix_full_final/late≈0.581/0.684`,
      `helix_core_final/late≈0.72/0.845`,
      `hbond_final/late≈32.16/31.36`
    - endpoint note:
      the last frame has weaker helix occupancy than the late-half mean, but
      the last-20-frame averages stay close to the retained baseline while RMSD
      and hbond metrics improve
  - fresh `1afo` sanity check on the same SC-only overlay:
    - output:
      `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/validate_20260630_1afo_cleanstart_target_rep8_sc_bulkbilayer10a_t12_rep/`
    - production kinetic ratio:
      `avg_kinetic_energy/1.5kT≈1.081`
    - bilayer metrics:
      `sep_late≈16.40 A`, `sep_max≈16.62 A`,
      `|n_z|_late≈0.99317`, `occ_cv_late≈0.322`,
      `nn_mean/p05≈7.50/6.59 A`, zero flips
    - protein metrics:
      `CA_RMSD_final/late≈5.14/4.42 A`,
      `DSSP_final/late≈0.972/0.962`,
      `hbond_final/late≈77.74/80.93`
    - interpretation:
      slight DSSP / late-hbond softening versus the retained target-only
      baseline, but no clear regression in bilayer behavior, kinetics, or
      overall protein stability

- Retained conclusion:
  - the `1rkl` stage-7 issue is not a bilayer-only transport defect
  - the measured release hold is useful only as a probe and is not a valid
    retained fix under universal settings
  - the next root-cause surface is the shared coarse protein-lipid force-field
    design, especially the terms that control helix-cap packing and tilt
  - the shared force-field issue now has a concrete build-time root cause:
    `cg_lipid_target` was both underresolved by a two-conformer isolated DOPC
    base ensemble and, before the latest source patch, numerically distorted by
    writing sampled PMF values directly into spline-control storage
  - the current best retained branch is the representative-target overlay:
    it improves the protein-facing behavior of both `1rkl` and `1afo` while
    keeping leaflet separation, orientation, and same-leaflet spacing
    acceptable under the shared coarse settings
  - the direct bilayer-derived target-only retraining branch is not retained.
    Even after fixing the pooled-bilayer representative-selector bug, the
    `1rkl` clean-start branch slows production kinetics and weakens late
    protein secondary-structure / hbond metrics relative to the current
    representative-target baseline.
  - the new best retained branch is the SC-only bilayer-derived overlay.
    Preserving the representative target term while rebuilding only
    `cg_lipid_sc` from the pooled bilayer ensemble improves `1rkl` on the main
    protein and transport metrics without materially regressing `1afo`.
  - the remaining blocker is promotion strategy rather than root-cause
    diagnosis: decide whether to promote the SC-only overlay into the live
    `dopc.h5` workflow default or keep it as an explicit overlay artifact until
    the user requests the promotion step.
