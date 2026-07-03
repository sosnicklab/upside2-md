# Project Goal

Replace the current hydrophobic-row-mixed CGL/DOPC artifact with a single
physically derived force-field artifact that fixes the rigid-DOPC
protein-alignment issue generically, without per-protein tuning, runtime
orientation bias, or other nonphysical coupling shortcuts, while preserving
accepted hybrid structural behavior and CGL transport timescale.

# Architecture & Key Decisions

- Keep the runtime model local and generic.
  No per-protein switches, no bilayer-aware twisting term, and no SC row
  splicing from different artifacts if a fully physical candidate is accepted.
- Prefer a single pooled physical reference ensemble over residue-source splicing.
  If a pure shell-derived all-row SC fit regresses validation, the next
  acceptable route is to rebuild SC from one uniform reference ensemble that
  combines trajectory-derived DOPC conformers with generic isolated DOPC
  conformers under the same fitting rule for every residue.
- Treat trajectory-derived DOPC representatives as probability-mass summaries,
  not extrema summaries.
  A physical reduced shell ensemble should represent common compaction states,
  so the representative selector must not force the absolute minimum and
  maximum shell conformers into every overlay fit.
- Balance pooled representatives across source trajectories.
  A generic multi-protein reference ensemble should not let the longest
  trajectory dominate the reduced conformer set purely by raw frame-lipid
  count; representative counts should be split evenly across source
  trajectories, capped only by each source pool size.
- Keep explicit SC compaction deltas on the same conformer family as the base
  SC table.
  When the base H5 already stores a physical overlay ensemble and physical
  compaction centers, the single-CGL SC retrofit should derive its compact and
  extended endpoint ensembles from that overlay set rather than from a
  separate isolated-DOPC MC pool.
- Accept the source-balanced full-bilayer base with overlay-center-matched SC
  endpoints as the final generic fix.
  The validated replacement keeps the pooled full-bilayer base SC table,
  rebuilds only the SC compact/extended delta endpoints from the same overlay
  reference ensemble around the stored physical centers, and leaves every
  bilayer-only CGL dataset identical to the installed baseline H5.
- Prefer completing the existing shell-derived table family over continuing the
  pooled isolated route once matched replay shows the pooled candidates regress.
  The installed target tables and hydrophobic SC rows are already shell-derived;
  the cleanest next fully physical candidate is to rebuild the remaining
  nonhydrophobic SC rows from the same shell-filtered reference ensemble.
- Stay inside the active dryMARTINI interface scope:
  `py/martini_build_tables.py`,
  `py/martini_gen_params.py`,
  `py/martini_prepare_system_lib.py`,
  `src/martini_cg_lipid.cpp`,
  `parameters/dryMARTINI/dopc.h5`,
  and directly related validation scripts only when needed.
- Prefer a single-source retrained artifact over code changes.
  The current remaining nonphysical piece is in the installed H5 content, not
  in the explicit SC/target runtime path.
- Preserve existing physical-model contracts.
  Do not disable hybrid interface interactions, do not reintroduce pair-side
  implicit gap response, and do not add orientational or leaflet-aware bias to
  the CGL force field.
- Accept a candidate only if both systems remain valid.
  `1rkl` and `1afo` must stay structurally acceptable, show no obvious forced
  vertical drift, and keep CGL lateral transport in the accepted range.

# Execution Phases

- [x] Review prior shell candidates, their validation runs, and installed-artifact
      provenance to select the cleanest single-source physical candidate.
- [x] Reject or accept the fully shell-derived all-row SC candidate against
      matched replay controls.
- [x] Implement pooled-reference overlay support so a single CG-SC rebuild can
      use shell-filtered and isolated DOPC conformers without residue-source
      splicing.
- [x] Build and gate the pooled-reference SC candidates on the preserved
      shell-target base.
      Both tested pooled candidates remained runtime-physical but failed the
      matched `1rkl` replay gate, so that route is not the active promotion
      path.
- [x] Build the shell-family completion candidate by rebuilding the remaining
      nonhydrophobic SC rows from the same shell-filtered reference ensemble as
      the accepted target and hydrophobic rows.
- [x] Test whether the shell reference ensemble itself is biased by using a
      broader full-trajectory shell conformer pool rather than only the late
      trajectory half when rebuilding the nonhydrophobic SC rows.
- [x] Probe whether the remaining vertical-bias signal is dominated by the
      hydrophobic shell rows by rebuilding those rows from the broadened
      full-trajectory shell ensemble before committing to a broadened all-row
      shell retrain.
- [x] Validate the broadened all-row shell retrain on matched `1rkl` replay.
      It preserved compactness and timescale but drove the protein nearly
      vertical, so that branch is rejected.
- [x] Patch the trajectory-derived DOPC representative selector so shell
      overlays use probability-centered conformers instead of forcing shell
      compaction extrema into every reduced ensemble.
- [x] Rebuild and screen the shell-plus-isolated candidate with the updated
      representative-selection rule to determine whether the selector bias was
      the dominant failure mode.
- [x] Rebuild the source-balanced pooled full-bilayer plus isolated candidate
      and use exact/proxy matched `1rkl` replays to isolate whether the
      remaining regression sits in the base SC table or in the explicit SC
      compaction-delta layer.
- [x] Replace the isolated-MC single-CGL SC compaction retrofit with an
      overlay-center-matched state selection, then re-gate the same
      source-balanced candidate on matched `1rkl`.
- [x] If `1rkl` passes, validate the same candidate on matched `1afo` replay
      for orientation, RMSD/Rg/hbonds, and CGL `Dxy`.
- [x] Summarize whether the accepted fix is fully physical, what evidence supports
      it, and any remaining caveats or blockers.
- [x] Promote the accepted H5 into `parameters/dryMARTINI/dopc.h5` with a
      backup of the previous default artifact.

# Known Errors / Blockers

- No active technical blocker remains for the accepted fix.
  The source-balanced full-bilayer base plus overlay-center-matched SC
  compaction retrofit is now validated on both replay systems.
- The rejected branches remain rejected.
  Pure shell-only completion, pooled shell-plus-isolated retrains, and the
  broadened shell families all produced either protein expansion or excessive
  vertical alignment on the matched replay harness.
- A fresh standalone 288-lipid rerun is not required to validate the last SC
  retrofit step.
  The accepted H5 is bitwise identical to the installed `dopc.h5` for every
  bilayer-only dataset and attribute, so the bilayer-only CGL clock and
  leaflet mechanics are exactly unchanged by the accepted patch.

# Review

- Accepted status:
  the generic physical fix is the source-balanced full-bilayer base H5 with
  SC compact/extended delta endpoints rebuilt from the same overlay reference
  ensemble around the stored physical compaction centers.
  No per-protein switch, no extra twisting/orientational term, and no bilayer-
  specific runtime bias were introduced.
- Validation summary:
  on matched `1rkl`, the accepted candidate gives
  `late hbonds 29.08` vs `30.62`,
  `late Rg 12.50 A` vs `12.84 A`,
  `late CA RMSD 1.86 A` vs `1.76 A`,
  `axis_trim1 0.740` vs `0.730`,
  `pca_full 0.888` vs `0.848`,
  `Dxy_half / quarter 0.439 / 0.443` vs `0.421 / 0.409`.
  On matched `1afo`, it gives
  `late hbonds 83.96` vs `80.11`,
  `late Rg 14.83 A` vs `14.88 A`,
  `late CA RMSD 2.68 A` vs `2.94 A`,
  `overall pca_full 0.996` vs `0.998`,
  `Dxy_half / quarter 0.459 / 0.456` vs `0.435 / 0.440`.
- Bilayer-timescale conclusion:
  the accepted H5 changes only `cg_lipid_sc/delta_extended` and
  `cg_lipid_sc/delta_compact`.
  All bilayer-only CGL tables and attrs are bitwise identical to the installed
  `parameters/dryMARTINI/dopc.h5`, so the standalone CGL bilayer timescale and
  mechanics are exactly preserved by this acceptance step.
- Delivery status:
  the accepted H5 is now installed at `parameters/dryMARTINI/dopc.h5`, with
  the previous default copied to `parameters/dryMARTINI/dopc.h5.bak`.
