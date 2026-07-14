# Progress

## Current Phase

- Confirmed that no workflow, preparation, Upside, or persistent `screen`
  process remains active.
- Verified all four requested wrappers through the real auto-continuation
  selector without running dynamics. Coarse 1RKL, coarse 1AFO, and full 1RKL
  select their active `stage_7.0.up`; full 1AFO finds no active checkpoint
  because its unstable stage 7.0 is intentionally archived and rejected.
- Stopped the two mistakenly restarted coarse workflows before they reached
  production. Their previously completed coarse `stage_7.0.up` checkpoints
  remain intact.
- Located the missing completed full-resolution `stage_7.0.up` checkpoints in
  `/Users/yinhan/Documents/upside2-md-bak`. They have restart-valid momentum,
  but current validation correctly rejects them: protein carrier mass is `6`
  instead of `1`, and their production timestep is `0.002` instead of `0.004`.
- Enabled latest-stage auto-selection by default in the four requested system
  wrappers. The shared workflow remains opt-in.
- Previously started four persistent, named `screen` loops. Both coarse
  workflows resumed from stage 7.0, completed stage 7.1, and automatically
  advanced to stage 7.2 without rerunning preparation; those generated stages
  were later rejected because the old continuation path minimized them.
- Previously started both full workflows fresh because the archived checkpoints violate
  current mass/timestep invariants. Full 1AFO became unstable abruptly at
  burn-in step 30,500 (`protein_potential > 10000`, sustained hbond loss) and
  was stopped at step 31,500. Full 1RKL and both coarse loops were subsequently
  stopped as well.
- Validation across the auto-generated coarse segments found non-exact restart
  boundaries (`0.45-0.83 A` position errors and large momentum errors) and
  strong orientation drift. Root cause: `run_stage70_continuation` minimized
  every copied production checkpoint, which replaced the promoted position and
  invalidated momentum before MD. All workflows are now stopped, and the
  continuation minimization has been removed.
- Isolated 500-step continuations now preserve position, momentum, and
  compaction exactly for both proteins; orientation differs only by one
  float32 normalization ULP (`1.19e-7`). No minimization is invoked and source
  hashes remain unchanged.
- The same test reveals that the retained fresh 1RKL stage 7.0 had already
  rotated from `34 deg` to roughly `4 deg`. This source checkpoint is rejected;
  executable/source parity must be established before a new baseline run.

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
- Expanded `cg_lipid_potentials.tex` into a reconstruction-level paper. The
  Methods now specify the source Hamiltonian, all conformer ensembles,
  tempered projection, spline equations and resolutions, physical pair/SC tail
  minimizations, compaction PMF, runtime derivatives, and exact GLE process.
- `latexmk` completes without warnings; all 10 final PDF pages pass visual
  inspection. The manuscript numbers match installed H5 metadata and accepted
  trajectories.
- Final syntax, diff, hash, and build checks pass. C++ and force-field builder
  files are unchanged from `52f637e`.

## Next

- Verify and rebuild the Upside executable from the restored source, rerun fresh
  coarse stage 7.0, then continue only if the protein-angle gate passes.
