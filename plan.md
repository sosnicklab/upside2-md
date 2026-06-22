# Project Goal

Resolve the persisted CGL z-coordinate/geometric bilayer issue visible in the existing 1AFO and 1RKL coarse VTF outputs, including CGL particles squeezed outside the bilayer in 1RKL.

# Architecture & Key Decisions

- Diagnose from the named outputs before changing code: `example/16.MARTINI/outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf` and `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf`.
- Preserve all hybrid SC-env and BB-env interactions; do not disable or scale physical interactions to hide geometry problems.
- Distinguish physical CGL center positions from VTF display rods, but treat bad rod geometry as a real output bug if it misrepresents the model.
- Keep fixes scoped to CGL preparation/extraction/workflow code unless output evidence points elsewhere.
- Do not use git write operations. All edits remain unstaged.

# Execution Phases

- [x] Inspect named VTF/H5/log outputs and determine whether they were generated before or after the latest salt/z-prep fix.
- [x] Quantify CGL center z, display rod z, and outlier/squeezed CGL behavior for 1AFO and 1RKL.
- [x] Compare output-time CGL geometry against available full-resolution references and prepared/start geometry.
- [x] Identify the runtime or extraction root cause.
- [x] Implement a valid physical/workflow fix that preserves calibrated CGL timescale.
- [x] Verify with targeted geometry checks, syntax tests, and no mobility/mass/restraint changes.

# Known Errors / Blockers

- The existing named output files under `example/16.MARTINI/outputs/` still
  reflect the old run until the full workflows are regenerated. The patched
  default workflows were verified in shortened fresh `/tmp` runs.

# Current Findings

- The named outputs were regenerated today and include explicit ions in
  `hybrid_prep_summary.json`: 1RKL has `93` ion atoms and 1AFO has `98`.
- 1AFO physical CGL centers are geometrically coherent at the final frame: no
  CGL center is more than `8 A` from its orientation-sign leaflet median.
- 1RKL has two physical CGL-center outliers in the final H5 frame: CGL ordinal
  `34` at raw `z=46.67 A` in the upper leaflet and CGL ordinal `55` at raw
  `z=-40.19 A` in the lower leaflet. These are real center positions, not only
  VTF display rods.
- Those 1RKL CGLs are already outliers at stage-7 frame 0 and are created
  during stage-6 relaxation from close cross-leaflet CGL neighbors. Preserving
  raw full-DOPC geometric COM z exactly makes the one-particle CGL initial
  state too compressed for the current conservative table.
- The VTF rod display path was replacing/renaming the original CGL center atom
  with a synthetic PO4 endpoint. That makes CGL z positions look shifted from
  DOPC COM even when the physical center is coherent.
- For 1RKL, two physical CGL centers are ejected by stage-6 minimization:
  `6.0.prepared` has no outliers, but stage-6 frame 0 already has CGL ordinals
  `34` and `55` at about `+/-42.8 A`. Stage-6 MD and stage 7 inherit those
  minimized coordinates.
- Raising the initial leaflet z separation is not robust: z-separation overrides
  of `15 A` and table-derived `18.8 A` can create larger minimizer ejections in
  one or both systems. The safer direction is to keep COM-derived CGL z
  positions and stabilize only the stage-6 minimization handoff.
- Invalidated finding: changing the CGL mass scale or adding z restraints can
  suppress the observed z motion but destroys the calibrated CGL timescale.
  Those are not acceptable fixes.
- The remaining 1RKL outlier is seed/basin-sensitive. Seeds `11`, `22`, and
  `2027` pass with calibrated `CG_LIPID_MASS_SCALE=0.012`, while `2026`
  produces a wrapped CGL z outlier after stage 6. The valid workflow fix is to
  reject bad stage-6 coarse CGL geometry and retry the same physical workflow
  with incremented seeds, rather than changing CGL mobility or adding forces.

# Revised Decisions

- Preserve the physical CGL center atom in VTF rod mode and append synthetic rod
  atoms instead of moving the center atom to a PO4 endpoint.
- Preserve the calibrated `CG_LIPID_MASS_SCALE=0.012` default.
- Do not add CGL z restraints or other mobility-limiting terms.
- Add bounded coarse-CGL stage-6 geometry rejection/retry. The retry regenerates
  stage-0 hybrid packing and reruns stage 6 with incremented seeds. It preserves
  `CG_LIPID_MASS_SCALE=0.012`, GLE parameters, force tables, and no CGL z
  restraints.

# Review

- Implemented VTF rod-mode output so the original CGL atom remains at the
  physical CGL center; synthetic PO4/C1A rod atoms are appended for display.
- Restored default CGL z-separation conditioning to `0.0`, preserving source
  DOPC COM-derived CGL z placement unless explicitly overridden.
- Reverted invalid dynamic changes: restored `CG_LIPID_MASS_SCALE=0.012` and
  removed temporary CGL z restraints.
- Added bounded coarse stage-6 CGL geometry retry after detecting wrapped-z
  outliers. Bad seed `2026` is rejected and retried with `2027`.
- Removed no diagnostic scripts in this pass because no matching unneeded
  diagnostic scripts remain in `example/16.MARTINI`.
- Previous shortened validation using mass/restraint changes is invalid and
  must not be used as evidence.
- Valid verification passed: 1RKL with bad seeds `2026/2026` rejected one
  stage-6 wrapped-z outlier and retried with `2027/2027`; final CGL
  `mass0=1.008`, wrapped `n_out=0`, max wrapped `|z|=18.54 A`. 1AFO with
  `2026/2026` required no retry; final CGL `mass0=1.008`, wrapped `n_out=0`,
  max wrapped `|z|=16.23 A`.
