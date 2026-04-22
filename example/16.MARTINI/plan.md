# 1RKL AABB Stage-7 Oxygen Failure

## Project Goal
- Diagnose the reported problems in `outputs/martini_test_1rkl_aabb`.
- Check whether bilayer drift is a real runtime error and whether it explains the stage-7 oxygen blow-up.
- Fix the smallest root cause in the active `1rkl` AA-backbone workflow so stage `7.0` keeps backbone oxygens physically attached when the backbone is released.

## Architecture & Key Decisions
- Treat the reported bilayer shift as a hypothesis, not the assumed root cause.
- Prefer direct artifact measurements from the reported run over VTF-only visual impressions.
- Keep the current AA-backbone workflow structure unless the diagnostics prove it is fundamentally unsound.
- If stage-7 oxygens are being ignored by the native Upside 3-site backbone nodes, add explicit local geometry constraints for the runtime `O` atoms instead of changing the broader workflow design.

## Execution Phases
- [x] Phase 1: Re-read local task trackers and inspect the reported `martini_test_1rkl_aabb` artifacts plus the active workflow script.
- [x] Phase 2: Quantify bilayer/protein motion from the HDF5 outputs and verify whether bilayer drift is the actual failure.
- [x] Phase 3: Audit the stage-7 backbone node injection and identify why explicit oxygens destabilize after the rigid hold is removed.
- [x] Phase 4: Implement the smallest code fix for stage-7 explicit oxygen geometry.
- [x] Phase 5: Run targeted verification on regenerated stage-7 artifacts and record the outcome.

## Known Errors / Blockers
- Current diagnosis from `outputs/martini_test_1rkl_aabb`:
  - preproduction bilayer drift relative to the fixed protein is small (`< 0.15 Å` in XY, `0 Å` in Z for PO4 COM),
  - the stage-7 failure is the explicit backbone oxygen coordinates: `C-O` distances grow from about `1.23 Å` to tens of Å.
- Root cause confirmed:
  - `inject_backbone_nodes(...)` writes the native Upside 3-site backbone machinery for `N/CA/C`,
  - explicit runtime `O` atoms are present in the stage file but receive no direct restoring bond/angle geometry once `fix_rigid` is removed.

## Review
- Kept the current AA-backbone workflow design and patched only the stage-file injector in `py/martini_prepare_system_lib.py`.
- Added explicit local geometry for runtime carbonyl oxygens after the native 3-site backbone nodes are injected:
  - one `C-O` bond per residue,
  - one `CA-C-O` angle per residue,
  - one `O-C-N(next)` angle for each nonterminal residue.
- Verification on a regenerated production handoff file showed:
  - original `outputs/martini_test_1rkl_aabb/checkpoints/1rkl.stage_7.0.up`:
    - `C-O` max = `102.05 Å`,
    - frame-10 `C-O` max = `3.26 Å`,
    - relative `PO4`/protein drift = `4.77 Å` XY, `4.16 Å` Z,
    - absolute `PO4` drift = `0.31 Å` XY,
  - fixed smoke artifact `outputs/martini_test_1rkl_aabb_oxygen_fix_smoke/checkpoints/1rkl.stage_7.0.handoff.up` after `1000` MD steps:
    - `C-O` range remained `0.96–1.58 Å`,
    - last-frame `C-O` mean = `1.20 Å`,
    - relative `PO4`/protein drift = `0.46 Å` XY, `0.26 Å` Z,
    - absolute `PO4` drift = `0.04 Å` XY.
