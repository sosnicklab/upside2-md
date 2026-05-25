# 2026-05-25 1AFO Coarse First-Frame Bend

## Project Goal
- Debug why `example/16.MARTINI/outputs/martini_1afo_hybrid` has one bent
  helix at the very first production frame while
  `martini_1afo_hybrid_full` has both helices relatively straight.
- Determine whether the issue is in stage-7 minimization, burn-in promotion,
  VTF/output extraction, or CGL interface physics.
- Fix only confirmed bugs with a physical change: no disabled, scaled-away, or
  bypassed SC-env, BB-env, lipid-lipid, CGL-SC, or CGL-target interactions.

## Architecture & Key Decisions
- Compare exact existing HDF5 outputs first; do not infer from previous copied
  probes.
- Use copied files for any runtime diagnostics that could write to HDF5.
- Treat the single-particle lipid as a coarse-grained version of the explicit
  DOPC model. Differences must come from valid coarse-graining, not missing
  force paths or workflow bookkeeping.
- N/CA/C/O carrier atoms must continue to receive additive force contributions
  from sidechain/rotamer coupling and projected BB-proxy environment coupling.
  Only true virtual BB proxy gradients are projected through the BB map.
- Stage handoff must be physically conservative: if minimization or burn-in is
  performed, the promoted `/input` coordinates/momenta/box must correspond to
  the intended state and the first saved production frame must match that
  state.

## Execution Phases
- [x] Phase 1: Quantify 1AFO coarse/full helix geometry at stage-6 output,
      stage-7 prepared input, stage-7 minimized/burn-in handoff, and stage-7
      first output frame.
- [x] Phase 2: Compare CGL target, CGL-SC, BB proxy, terminal charge, and
      chain-break metadata between coarse and full stage-7 files.
- [x] Phase 3: Inspect minimization and burn-in workflow code for stale
      coordinate promotion, wrong output frame source, or coarse-only node
      injection differences.
- [x] Phase 4: Implement the smallest physical bug fix and update MARTINI
      documentation if model behavior changes.
- [x] Phase 5: Verify with focused copied-file regeneration/minimization and
      geometry metrics before reporting completion.

## Known Errors / Blockers
- The old `martini_1afo_hybrid` stage-7 output was stale relative to the
  current C++ force-routing binary. Its stage-7 minimization log started from
  `162.96 E_up`, while the current binary evaluates the same prepared file at
  `2867.25 E_up` and minimizes to a different physical basin.
- The ignored generated files under
  `example/16.MARTINI/outputs/martini_1afo_hybrid` were regenerated for
  stage 7 from the current binary. Source code did not need another physics
  change beyond the already-applied direct-vs-projected carrier routing fix.

## Review
- The minimizer itself was not the source of the first-frame bend: copied
  minimization-only geometry changed 1AFO fragment bends only from
  `40.50/42.89 deg` to `40.07/42.19 deg`.
- The first production frame comes after a 40k-step stage-7 burn-in promotion,
  not immediately after minimization. With the current binary, that promoted
  coarse state has fragment bends `28.98/21.38 deg` and hbond `82.37`.
- Updated generated output verification:
  `martini_1afo_hybrid/checkpoints/1afo.stage_7.0.up` now has first-frame
  fragment bends `28.98/21.38 deg` and final hbond `84.70`; the VTF was
  regenerated from that checkpoint.
