# Progress Log

## Current Task: Generic `1rkl` / `1afo` CGL Compaction-State Repair

- Actions taken:
  - Rechecked the `1afo` regression against the matching full-resolution
    workflow and confirmed they start from the same packed structure.
  - Measured the saved coarse `stage_7.0` collapse and verified the
    full-resolution control remains stable from the same starting pack.
  - Inspected the installed `parameters/dryMARTINI/dopc.h5` and found that the
    original `cg_lipid_compaction/self_coeff` was zero and that the stored
    physical compact/extended tail centers were stale.
  - Patched `py/martini_build_tables.py` so the compaction reference centers,
    normalized hidden-state coordinate, and self PMF can be rebuilt from the
    stored DOPC reference ensemble, and so stale `cg_lipid_sc` and
    `cg_lipid_target` compaction groups are detectable.
  - Profiled the exact pair-relaxation rebuild path and confirmed it is not a
    practical first-line repair locally:
    a representative pair-relax evaluation costs about `0.028 s`, implying a
    multi-hour full re-fit.
  - Verified in `src/martini_cg_lipid.cpp` that the runtime uses a linear
    SC/target compaction mix and a bilinear CGL-CGL compaction mix.
  - Replaced the slow rebuild script with an exact endpoint-remap repair that
    uses `parameters/dryMARTINI/dopc.h5.bak` as the old physical state
    contract and analytically rewrites the stale SC/target/pair endpoint tables
    onto the repaired physical centers already stored in the live H5.
  - Quantified the user-reported coarse-only `1rkl` bending against the
    matching full-resolution control with a backbone line-RMS metric and
    reproduced the bend from the saved stage-7 start state.
  - Ran same-start-state A/B replays that isolated the endpoint-remapped
    pair tensor as the main new bend source, with a smaller contribution from
    the remapped target tables.
  - Iterated through generic H5 variants and found that the best cross-system
    fix keeps the current SC correction, restores the old pair/target
    compaction tables, and repairs the self PMF around the old physical
    compaction contract from `dopc.h5.bak`.
  - Overwrote `parameters/dryMARTINI/dopc.h5` with that delivered artifact.
  - Reinjected the delivered tables into copies of the saved
    `1afo.stage_7.0.up` and `1rkl.stage_7.0.up` files and ran 10,000-step
    production continuations from those saved stage-7 restart states.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  - `/Users/yinhan/Documents/upside2-md/findings.md`

- Verification:
  - `python -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py`
    succeeded under the project environment.
  - The delivered live H5 now reports the restored old physical state contract
    with a repaired nonzero self PMF:
    `reference_extended_center_ang=13.381275 A`,
    `reference_compact_center_ang=19.615337 A`,
    `compact_state_probability=0.3023485`,
    `self_coord_min_ang=-0.0874`,
    `self_coord_max_ang=1.1418`.
  - `1rkl` validation:
    the bad coarse source ends at backbone line-RMS `4.6494`, while the
    delivered live-H5 continuation from the same saved stage-7 start ends at
    `3.8695`; the full-resolution control ends at `3.6087`.
  - `1afo` validation:
    the delivered live-H5 continuation stays away from the earlier collapsed
    state, with chain-center XY separation
    `11.5248 -> 10.6894` across 10,000 steps.

- Failures and fixes:
  - A remap-only repair of the initial compaction coordinate did not hold:
    the state still relaxed back to the compact endpoint, which exposed the
    missing self-PMF defect.
  - Repairing only `cg_lipid_compaction` removed the hidden-state saturation
    but did not recover `1afo` separation, which exposed the remaining stale
    compact-vs-extended pair/SC/target correction tables.
  - The exact representative-state pair-relaxation rebuild was too expensive
    to use as the main delivery path, so the fix was changed to an exact
    analytical endpoint remap inside the runtime’s existing linear/bilinear
    state model.
  - The analytical endpoint remap was not enough:
    it helped `1afo`, but it also introduced a new coarse-only `1rkl`
    bending mode because the widened physical state contract suppressed the
    compression side assumed by the old pair/target tables.
  - Allowing the full unclipped old-contract hidden-state range was also too
    aggressive; it drove the compaction coordinate out to about `2.85` and
    made the bend worse again.
