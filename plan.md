# Project Goal

Fix the remaining `1rkl` coarse-hybrid protein defect at the shared
force-field / implementation level without introducing protein-specific
simulation settings, non-physical bilayer restraints, or bilayer-state terms
inside the CGL force field.

# Architecture & Key Decisions

- Treat `stage_6.0` as the structural baseline.
  The cited `1rkl` artifact keeps native secondary structure through stage 6,
  so the defect starts later in the hybrid workflow.
- Separate stage-7 activation defects from equilibrium force-field defects.
  Replays must start from the promoted stage-7 burn-in input so the diagnosis
  targets the actual release into unrestrained production.
- Reject benchmark-specific workflow tweaks as retained fixes.
  A stage-7 release hold that helps `1rkl` but harms `1afo` is useful as a
  diagnostic probe only, not as a production method.
- Keep universal simulation settings across proteins.
  Any retained fix must live in shared force-field construction or shared
  runtime code and must not depend on per-protein launcher overrides.
- Focus on the coarse protein-lipid interface before revisiting transport.
  Current evidence points to `1rkl` coarse helix-cap / tilt sensitivity rather
  than a bilayer-only CGL transport failure or a generic Upside core collapse.
- Treat the default isolated-DOPC reference ensemble as part of the force
  field, not as a harmless build shortcut.
  The live `dopc.h5` was built from only two isolated DOPC conformers, and
  direct probes show that this makes `cg_lipid_target` support and PMFs
  seed-sensitive in exactly the cap-sensitive regime that differentiates
  `1rkl` from `1afo`.
- Retain the representative-ensemble `cg_lipid_target` overlay as the current
  practical branch.
  The validated near-term fix preserves the accepted pair/SC/compaction tables
  from the clean baseline `dopc.h5` and rebuilds only `effective_lj` plus
  `cg_lipid_target` from the `32 -> 8` isolated-DOPC representative ensemble.
- Treat the remaining `1rkl` gap as an SC-CGL interfacial design problem until
  disproven by clean-start validation.
  The target-overlay branch fixed bilayer behavior and improved protein-facing
  backbone terms, but the remaining cap / tilt mismatch is still most
  consistent with the preserved baseline `cg_lipid_sc` surface.
- Treat clean-start workflow validation as a separate acceptance gate.
  The earlier continuation-based `1rkl` and `1afo` replays were sufficient to
  compare protein and bilayer behavior, but promotion now requires a fresh
  workflow run that goes through the normal stage-6 to stage-7 handoff in a
  new run directory with the retained overlay H5.
- Retain the bilayer-derived SC-only overlay as the current best candidate.
  The latest source fix preserves the representative isolated-DOPC target term
  and rebuilds only `cg_lipid_sc` plus its SC-specific compaction corrections
  from the pooled bilayer-derived bulk-DOPC ensemble.

# Execution Phases

- [x] Quantify the `1rkl stage_7.0` defect and confirm it starts after stage 6.
- [x] Reject stage-7 activation-hold fixes as universal settings.
- [x] Map the failing coarse force-field terms against the localized `1rkl`
      helix-cap and tilt defect.
- [x] Implement and validate the representative-target overlay as the first
      shared force-field correction.
- [x] Run clean-start validation for `1afo` with the representative-target
      overlay H5 and remeasure protein and bilayer behavior.
- [x] Rebuild the representative overlay with a source-derived `cg_lipid_sc`
      surface and validate clean-start `1rkl`.
- [x] Replace the isolated-lipid interface reference ensemble with a
      bilayer-derived bulk-DOPC ensemble for the protein-facing CGL tables,
      then validate clean-start `1rkl`.
- [x] Preserve the current best target term and rebuild only `cg_lipid_sc`
      from the bilayer-derived bulk-DOPC ensemble, then validate clean-start
      `1rkl`.
- [ ] If the bilayer-derived interface branch still misses `1rkl`, isolate the
      exact sidechain classes / force-field terms responsible and redesign that
      shared interface term directly.
- [x] Re-run bilayer and `1afo` sanity checks on the retained `1rkl` fix.

# Known Errors / Blockers

- The measured post-burn-in release hold remains in source as optional
  instrumentation, but it is not a retained fix and must stay disabled by
  default.
- The live `parameters/dryMARTINI/dopc.h5` remains on the last clean
  internally consistent baseline (`pair/sc/target=2/2/2` conformers).
  The new best branch still lives as an overlay H5 plus a source-level overlay
  builder, not as the promoted default file.
- The representative isolated-DOPC SC overlay does not close the `1rkl` gap.
  It improves late helix occupancy but shifts the whole helix upward relative
  to the membrane midplane, worsens final aligned CA RMSD, and slightly
  enlarges the leaflet spacing.
- The runtime geometry audit has not yet found a simple target-angle or
  cross-leaflet sign-convention bug.
  The C++ and Python implementations agree on the core target and
  cross-leaflet face-angle conventions checked so far.
- The next root-cause branch is the training/reference ensemble itself:
  the user-identified defect is bilayer-compacted DOPC tails missing from the
  protein-facing PMF construction, so the next redesign should derive the
  interface ensemble from bulk bilayer lipids rather than from isolated DOPC
  conformers.
- The target-only bilayer-derived interface rebuild is rejected.
  Two clean-start `1rkl` runs on pooled-bilayer target-only overlays kept the
  bilayer ordered but worsened production kinetics and overall protein
  acceptance metrics relative to the retained representative-target baseline.
- The SC-only bilayer-derived overlay looks retained as a force-field fix, but
  it still has a mild endpoint-interpretation caveat for `1rkl`.
  Late-half `1rkl` RMSD and hbond metrics improve materially, the production
  kinetic ratio improves to `≈1.262`, and the bilayer stays correct, but the
  final production frame has lower helix occupancy than the late-half mean.
  The last-20-frame averages stay close to the retained baseline, so this
  currently reads as an endpoint fluctuation rather than a clear instability.
- Fresh `1afo` sanity validation on the SC-only overlay stays within the
  accepted envelope rather than showing a benchmark-breaking regression.
  Relative to the retained target-only baseline, `1afo` keeps a correct
  bilayer and nearly unchanged production kinetics, with only a small decrease
  in DSSP / late hbond signal that does not suggest a secondary-structure
  collapse.

# Review

- Confirmed facts to carry forward:
  - `1rkl stage_6.0` remains structurally stable.
  - The large drift appears only after release into unrestrained
    `stage_7.0` production.
  - The previously tested `1rkl`-only release hold is not acceptable because
    it degrades `1afo`.
  - The most likely universal defect surface is the coarse protein-lipid
    interface that sets helix-cap packing and tilt preferences.
  - The current protein-facing CGL force field is underresolved at build time:
    `cg_lipid_target` was built from a two-conformer isolated DOPC base
    ensemble, and the old target table also misused sampled PMF values as
    spline coefficients instead of fitted spline controls.
  - The current best retained branch is the representative-target overlay:
    keep the accepted baseline pair/SC/compaction tables, rebuild
    `effective_lj` plus `cg_lipid_target` from the `32 -> 8` isolated-DOPC
    representative ensemble, and route workflows to that H5 via `DOPC_H5`.
  - The SC-overlay clean-start `1rkl` run did not validate as a retained fix.
    Rebuilding `cg_lipid_sc` from the same isolated-DOPC representative
    ensemble raises helix occupancy, but it also lifts the helix relative to
    the bilayer center and worsens the final global pose. This means the
    remaining defect is deeper than simply rebuilding SC-CGL from a larger
    isolated reference ensemble.
  - The next critical falsification step is a bilayer-derived interface
    rebuild that samples bulk DOPC conformers from full-resolution bilayer
    trajectories while excluding lipids near embedded proteins. If that closes
    the `1rkl` gap, then the residual issue was the reference-ensemble design
    itself rather than the runtime convention.
  - That broad bilayer-derived target-only rebuild did not close the gap.
    The next targeted falsification step is therefore SC-only:
    keep the retained representative-target overlay fixed and replace only the
    `cg_lipid_sc` term with a bilayer-derived bulk-ensemble retrain.
  - That SC-only bilayer-derived retrain is now the current best branch.
    The new overlay at
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_dopc_target_rep8_sc_bulkbilayer10a_t12_rep_overlay/dopc_target_rep8_sc_bulkbilayer10a_t12_rep_overlay.h5`
    preserves the representative target term
    (`target/effective conformer_count=8`, generic interface source unchanged)
    while rebuilding only `cg_lipid_sc` from the pooled bilayer ensemble
    (`sc conformer_count=12`, SC-specific reference dataset and compaction
    retrofit). Clean-start `1rkl` improves to
    `avg_kinetic/1.5kT≈1.262`, `RMSD_final/late≈2.79/2.05 A`,
    `hbond_final/late≈32.16/31.36`, with a correct bilayer
    (`sep_late≈16.11 A`, `|n_z|≈0.9933`, `nn_p05≈6.44 A`, `occ_cv≈0.226`,
    zero flips). Fresh `1afo` on the same overlay keeps a correct bilayer and
    healthy kinetics (`avg_kinetic/1.5kT≈1.081`) with only minor DSSP/hbond
    softening relative to the retained target-only baseline.
  - Clean-start promotion requires a new run directory. The generic launcher
    auto-detects a prior `stage_7.*` checkpoint when `RUN_DIR` is reused, so
    validation must use a fresh output path rather than the old continuation
    directories.
  - Fresh `1afo` validation on the representative-target overlay completed
    through the normal stage-6 to stage-7 workflow at
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/codex_1afo_cleanstart_target_rep8/`.
    The production branch keeps a correct bilayer
    (`sep_late≈16.70 A`, `|n_z|≈0.9933`, `nn_p05≈6.74 A`, zero flips) and a
    normal kinetic ratio (`avg_kinetic/1.5kT≈1.053`), while protein secondary
    structure remains high (`DSSP_final≈0.986`, `DSSP_late≈0.970`).
    Relative to the older coarse baseline, the clean-start overlay run has
    better hbonds but a larger final aligned CA RMSD (`≈5.14 A` vs `≈3.52 A`),
    so the remaining protein question is whether that extra drift is an
    acceptable pose/tilt excursion or needs another force-field iteration.
