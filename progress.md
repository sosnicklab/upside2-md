# Progress Log

## Current Task: `1afo` CGL Compaction-State Repair

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
  - Overwrote `parameters/dryMARTINI/dopc.h5` with that repaired artifact.
  - Reinjected the repaired tables into copies of the saved
    `1afo.stage_7.0.up` and `1rkl.stage_7.0.up` files and ran 10,000-step
    production continuations from those saved stage-7 restart states.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`
  - `/Users/yinhan/Documents/upside2-md/rebuild_dopc_h5_once.py`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5.bak`

- Verification:
  - `python -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py`
    succeeded under the project environment.
  - The repaired live H5 now reports a consistent state contract:
    `cg_lipid_compaction`, `cg_lipid_sc`, and `cg_lipid_target` all use
    runtime centers `0.0/1.0`, physical reference centers
    `12.134849 / 26.963173 A`, and updated compact-state probability
    `0.5833333`.
  - The remapped tables differ materially from the backup old-state tables,
    especially on the compact side:
    `cg_lipid_sc/delta_compact` RMS shift `~1.28`,
    `cg_lipid_target/delta_compact` RMS shift `~3.66`,
    `cg_lipid_compaction/delta_compact_compact` RMS shift `~2.15`.
  - `1afo` repaired continuation:
    from the same saved final stage-7 frame, chain-center XY separation moved
    `8.41 -> 8.79 A` and mean CGL compaction moved `0.9999 -> 0.9178`.
  - `1rkl` repaired continuation:
    10,000 steps completed stably; a simple full-chain backbone principal-axis
    metric changed only `31.15 -> 28.88 deg`, which does not indicate a return
    to the obvious vertical-collapse failure mode.

- Failures and fixes:
  - A remap-only repair of the initial compaction coordinate did not hold:
    the state still relaxed back to the compact endpoint, which exposed the
    missing self-PMF defect.
  - Repairing only `cg_lipid_compaction` removed the hidden-state saturation
    but did not recover `1afo` separation, which exposed the remaining stale
    compact-vs-extended pair/SC/target correction tables.
  - Rebuilding the full retrofit from stdin failed because Python
    multiprocessing cannot reload `<stdin>` as a module, so the work was moved
    into the real file `rebuild_dopc_h5_once.py`.
  - The exact representative-state pair-relaxation rebuild was too expensive
    to use as the main delivery path, so the fix was changed to an exact
    analytical endpoint remap inside the runtime’s existing linear/bilinear
    state model.
  - A partial old-control `1afo` continuation was started for an extra A/B but
    interrupted once the repaired restart already showed the expected outward
    separation and reduced compaction saturation.
