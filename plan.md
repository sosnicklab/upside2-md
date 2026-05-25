# 2026-05-25 1RKL Bend and MARTINI Equivalence Audit

## Project Goal
- Debug the new strange bend in the 1RKL single-particle lipid result.
- Ensure charged terminal BB typing is physically consistent: every protein
  strand/fragment has an N-terminal `Qd` and C-terminal `Qa`, including broken
  proteins such as 1AFO, in both full-resolution and single-particle lipid
  workflows.
- Audit lipid-lipid and lipid-sidechain interaction ownership/equivalence
  between full-resolution lipid and single-particle CGL workflows.
- Fix confirmed bugs without disabling, scaling away, or bypassing hybrid
  protein/environment interactions.

## Architecture & Key Decisions
- Keep all SC-env, BB-env, lipid-lipid, CGL-CGL, CGL-SC, CGL-target, and
  generic dry-MARTINI interactions active.
- Treat single-particle lipid as a coarse-grained representation of the same
  dry-MARTINI DOPC model, not a different protein or lipid force field.
- First audit exact current HDF5 outputs and generated nodes before editing
  model code again.
- Use copied HDF5 files for runtime checks; do not run diagnostics directly on
  original output checkpoints.
- Terminal charge assignment must be residue-fragment based and must not depend
  on whether lipids are explicit DOPC or CGL.
- User correction: N/CA/C/O carrier atoms must still receive the intended force
  contributions from both sidechain/rotamer coupling and backbone/environment
  coupling.  Do not fix full/coarse mismatch by removing their environment
  force path.
- Revised direction: verify the C++ derivative accumulation path.  Carrier
  atoms need additive force accumulation from sidechain/rotamer and projected
  backbone/environment paths; direct carrier forces, when present, must be
  added to the carrier atom and not re-projected as BB-proxy forces.
- Confirmed C++ bug: mapped N/CA/C/O carriers were treated as projectable BB
  proxy sites because they have `ROLE_BB` and a BB map index.  Direct carrier
  forces must be added to the carrier atom; only virtual BB proxy forces should
  be projected through the BB map.

## Execution Phases
- [x] Phase 1: Identify the exact 1RKL single-particle bend timing and compare
      geometry against full-resolution lipid outputs.
- [x] Phase 2: Audit terminal BB type/charge assignment for 1RKL and broken
      1AFO in both workflows.
- [x] Phase 3: Audit lipid-lipid and lipid-SC ownership/equivalence between
      full-resolution and single-particle workflows.
- [x] Phase 4: Isolate and implement the smallest physical bug fix.
- [x] Phase 5: Run focused verification and update the MARTINI documentation.

## Known Errors / Blockers
- Existing `example/16.MARTINI/outputs` artifacts are stale with respect to the
  C++ projection fix until the workflows are rerun.
- A copied 1RKL run with every carrier added as an independent CGL target failed
  (`hbond=6.5`, kinetic ratio `2.096`; after projection correction, 10k-step
  hbond still fell to `9.6`), so independent CGL target copies are rejected.

## Review
- Corrected after user feedback: N/CA/C/O carrier atoms remain force
  accumulation sites for sidechain/rotamer and backbone/environment
  contributions. The fix is not to remove their force path; it is to prevent
  direct carrier gradients from being re-projected as virtual BB-proxy
  gradients.
- Final accepted implementation keeps CGL targets on BB proxies and non-protein
  point targets, keeps charged BB endpoints per fragment, and routes true
  virtual BB-proxy gradients through the BB map only.
- Focused validation passed: C++ build, Python compile, 1RKL 40k-equivalent
  copied probe (`hbond=31.43`, kinetic ratio `0.989`), and 1AFO fragment target
  audit (`72` BB proxies plus `98` ions, no missing BB proxy targets).
