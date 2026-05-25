# 2026-05-25 1AFO Coarse-Lipid Helix Bend

## Project Goal
- Diagnose why one helix still bends in the 1AFO single-particle lipid
  workflow while both helices remain straight in the full-resolution lipid
  workflow.
- Determine whether the observed artifact is from stale outputs, incomplete
  CGL-SC rotamer coupling, or another physical mismatch in the coarse lipid
  model.
- Fix any confirmed bug without disabling, weakening, or bypassing hybrid
  protein/environment interactions.

## Architecture & Key Decisions
- Keep SC-env, BB-env, CGL-CGL, CGL-SC, and CGL-target interactions active.
- Treat full DOPC and single-particle CGL as alternate representations of the
  same dry-MARTINI lipid environment.
- The current 1AFO coarse outputs already contain `cg_lipid_rotamer_sc`; the
  remaining bend is not from stale missing CGL-SC coupling.
- Use copied HDF5 files for runtime checks; do not run diagnostics directly on
  original output checkpoints.
- The 1AFO kink is introduced during stage-7 burn-in, not during the stage-7
  handoff minimization.
- Charged terminal BB bead types (`Qd`, `Qa`) must stay charged in the generic
  dry-MARTINI force field, but the single-particle CGL target projection should
  treat those terminal targets like ions: excluded-volume-only against CGL.  The
  directional explicit-DOPC charged-target PMF is not transferable as an
  additive attraction to the CGL projection after lipid headgroup, solvent, and
  ion degrees of freedom have been integrated out.

## Execution Phases
- [x] Phase 1: Check whether the current 1AFO coarse outputs were regenerated
      with `cg_lipid_rotamer_sc` or are stale.
- [x] Phase 2: Quantify 1AFO helix geometry and component energies in coarse
      versus full outputs.
- [x] Phase 3: Isolate any remaining implementation or model mismatch and make
      the smallest physical fix.
- [x] Phase 4: Run focused verification and document residual risks.

## Review
- Implemented `cg_lipid_target_charged_bb_excluded_volume` for charged
  terminal BB proxy targets only.  The fix keeps charged termini in the
  dry-MARTINI protein model while preventing reuse of the non-transferable
  attractive explicit-DOPC charged-target PMF as an additive CGL attraction.
- Verification on a regenerated copied 1AFO stage-7 burn-in kept both helices
  in the straight/full-lipid range: final half-chain angles were `32.34 deg`
  and `36.88 deg`, versus `37.33 deg` and `85.72 deg` in the old coarse output.

## Known Errors / Blockers
- None yet.
