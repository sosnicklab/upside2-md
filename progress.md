# Progress

## Current Phase

- Restored every modified C++ and force-field builder file exactly to
  `52f637e`; removed the inactive explicit-state/sieve implementation.
- Applied the narrow transport delta: standard Verlet remains one step,
  coarse equilibration and production use `dt=0.004`; the committed GLE
  parameters are unchanged.
- Added one preparation invariant required by the time contract: protein
  carriers use Upside mass `1`; dryMARTINI environment particles retain
  native mass divided by 12.
- Static checks pass: Python compilation, shell syntax, mass classification,
  unchanged physical runtime diff, and installed DOPC H5 hash
  `ba7049317e69a371812c69da8876edaa4ce81c21bef4c75a82a761176862b46c`.
- Git history located the requested earlier calibration at `9fb81b5`:
  mass `0.012`, rotational tau `0.008`, and the same exact-step GLE
  schedule. Commit `52f637e` later changed only the mass default to `0.3`
  while adding the required idempotent-scaling fix.
- Fresh one-step A/B rejected direct restoration of `0.012` on the later
  pair-relaxed table: raw `D=0.910 A^2/T_up`, leaflet separation
  `16.4 -> 24.9 A`, maximum leaflet reassignment `0.337`, and degraded
  orientation.
- Accepted `mass_scale=0.05` after estimator and literature reconciliation.
  The fresh 404-CGL run records the exact one-step clock, remains finite and
  ordered. Multi-origin x56 diffusion converges to `11.47/11.45 um^2/s`
  over `10-40/20-60 T_up`, matching the published `11.5 um^2/s` DOPC value
  at 303 K. A matched `0.075` early probe is rejected as too fast.
- Fresh short 1rkl and 1afo stage-7.0 runs pass the initial hybrid gate. Both
  record exactly `8 T_up` for 2,000 steps, keep protein mass `1`, preserve the
  initial secondary structure and protein-bilayer angle, retain closed ordered
  bilayers, and have finite coordinates and energies.
- 1rkl final hbond/DSSP/angle are `35.8/0.903/36.7 deg`; 1afo final
  hbond/DSSP are `89.1/1.000`, with chain angles remaining within their initial
  `15-18 deg` range.
- Sequential stages 7.1-7.3 also pass. Over 320 ns, 1RKL has mean/final DSSP
  `0.925/0.903`, mean angle `35.9 deg` versus `34.4 deg` initially, and finite
  stage-7.3 dynamics. 1AFO has mean/final DSSP `0.989/1.000`, retains both
  initial shallow angles, and remains finite.
- Rewrote `cg_lipid_potentials.tex` as a cohesive Methods--Results manuscript
  using the installed force-field provenance and fresh validation. `latexmk`
  completes without warnings; all seven rendered pages pass visual inspection.
- Final syntax, diff, hash, and build checks pass. C++ and force-field builder
  files are unchanged from `52f637e`.

## Next

- No required implementation work remains. Longer independent replicas would
  narrow transport uncertainty but are not needed to resolve the reported bug.
