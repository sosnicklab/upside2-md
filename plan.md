# 2026-05-27 MARTINI Full-Lipid Secondary Structure and CGL Table Fix

## Project Goal
- Fix bilayer orientation in CG-lipid simulations. Lipids around protein are "messy."
- All fixes must be physical — no parameter twisting, no additional orientational potentials.
- Everything must derive from ITP force field files.

## Architecture & Key Decisions
- **Root cause identified (2026-05-28)**: The CG-CG B-spline had severe angular underdetermination.
  15×15=225 angular control points fitted to 7×7=49 angular data points per radial distance
  (4.6× underdetermination). With Tikhonov regularization λ=0.01, the B-spline developed
  spurious oscillations between sample points, creating unphysical angular energy features.
- **Fix**: Reduce `n_knot_angular` from 15 to `min(cos_theta_count + 2, 15)` = 9, and increase
  `smooth` from 0.01 to 0.1. This matches knot resolution to data resolution and provides
  adequate regularization.
- The same underdetermination pattern exists in the SC-CGL table (9²=81 data → 15²=225 controls,
  2.8× ratio) but is less severe. Not yet addressed.
- No additional orientation potentials. No parameter twisting. B-spline coefficients are
  fully derived from ITP-defined bead interactions.

## Execution Phases
- [x] Phase 1: Remove `cg_lipid_leaflet_orientation` code, injection, and docs.
- [x] Phase 2: Audit current regenerated `sidechain.h5` and full-resolution stage files.
- [x] Phase 3: Patch the confirmed full-resolution force-field bug physically.
- [x] Phase 4: Add CG-CG excluded-area nonnegativity (contact_nm_nonnegativity).
- [x] Phase 5: Diagnose root cause of messy bilayer orientation → angular B-spline underdetermination.
- [x] Phase 6: Fix CG-CG B-spline: n_knot_angular=9, cg_smooth=0.1. Bilayer test PASSES.
- [ ] Phase 7: Re-run full protein+lipid workflow to verify orientation improvement.

## Known Errors / Blockers
- Full workflow re-run needed to confirm the B-spline fix resolves orientation in protein+lipid system.
- SC-CGL table may have similar (less severe) underdetermination issue; not yet addressed.
- r[5] (8.9 Å) has an extreme raw energy value (1971.7 kJ/mol) that dominates the B-spline fit.
  May need investigation of relaxation convergence at this distance.

## Files Changed
| File | Change |
|---|---|
| `py/martini_build_tables.py` | Added `n_knot_angular`, `cg_smooth` params to `_fit_cg_lipid_quadspline`; call site passes `n_knot_angular=min(_cg_ct+2,15)`, `cg_smooth=0.1` |
| `example/16.MARTINI/build_martini_h5_m1_temp.sh` | Updated to pass new B-spline params |
| `example/16.MARTINI/cg_lipid_potentials.tex` | Updated N_theta 15→9, added regularization rationale |
| `parameters/dryMARTINI/dopc.h5` | Regenerated with 9×9 angular knots, smooth=0.1 |
