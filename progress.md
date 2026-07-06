# Progress Log

## Current Task: Generic `1rkl` / `1afo` CGL Compaction-State Repair

- Actions taken:
  - Patched `py/martini_build_tables.py` so the single-CGL compaction retrofit
    now prefers reference-ensemble physical centers in `auto`, keeps the
    stored calibrated compact prior, and fits the self PMF on the bounded
    center-derived hidden coordinate instead of the unstable fully-unclipped
    stored-contract coordinate.
  - Patched `_apply_single_cgl_compaction_corrections_to_h5()` so a partial
    self-PMF update no longer deletes untouched pair-compaction tensors from
    `cg_lipid_compaction`.
  - Overwrote `parameters/dryMARTINI/dopc.h5` from the live backup with the
    bounded reference-ensemble self-PMF correction while preserving the
    existing pair-compaction tensors.
  - Reinjected fresh copies of
    `example/16.MARTINI/outputs/martini_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`
    and
    `example/16.MARTINI/outputs/martini_1afo_hybrid/checkpoints/1afo.stage_7.0.up`
    from the installed H5 and confirmed the seeded compaction coordinate is
    bounded (`q mean ~= 0.815`, `q max ~= 1.085`) instead of clamping into the
    runaway compact basin.
  - Ran fixed-seed 10,000-step validations from those reinjected restart
    copies:
    `/private/tmp/1rkl_stage7_reference_bounded_10k.up` and
    `/private/tmp/1afo_stage7_reference_bounded_10k.up`.
  - Measured the same acceptance diagnostics against the bad coarse sources
    and the full-resolution controls.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
  - `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5.pre_reference_ensemble_bounded_compaction_fix.bak`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`

- Verification:
  - `python -m py_compile py/martini_build_tables.py`
    succeeded under the project environment after the generator/apply-path
    changes.
  - Installed H5 verification:
    `parameters/dryMARTINI/dopc.h5` now has the bounded compaction contract
    with
    `reference_extended_center_ang=12.13485 A`,
    `reference_compact_center_ang=26.96317 A`,
    `compact_state_probability=0.3023485`,
    `self_coord_min/max=-0.06337/1.08474`,
    and it preserves the existing pair-compaction datasets
    `delta_extended_extended`, `delta_extended_compact`,
    `delta_compact_compact`.
  - `1rkl` validation from fresh reinjected stage-7 restart:
    CA RMSD to full-resolution improved from `3.55 A` on the bad coarse
    source to `2.81 A`;
    local interface score improved from `0.1350` to `0.0940`;
    late H-bond mean improved from `23.26` to `29.41`
    (full-resolution control `29.50`);
    late compaction mean moved from `1.0578` to `0.9533`.
  - `1afo` validation from fresh reinjected stage-7 restart:
    total CA RMSD to full-resolution improved from `4.25 A` to `2.91 A`;
    chain A RMSD improved from `2.51 A` to `1.74 A`;
    chain B RMSD improved from `3.70 A` to `2.71 A`;
    late H-bond mean improved from `73.65` to `80.21`
    (full-resolution control `84.17`);
    late compaction mean moved from `1.0622` to `0.9503`.

- Failures and fixes:
  - Fully unclipping the stored-contract hidden coordinate was not valid:
    it reopened a runaway compact basin (`q` up to about `2.85`) and destroyed
    the protein H-bond network even though it improved the coarse bend metric.
  - Rebuilding from an older backup that lacked the pair-compaction tensors
    produced a structurally stale H5 that the injector rejected.
    The fix was to preserve the live pair-compaction payload and patch the
    H5 apply path so a self-only update does not delete those untouched
    datasets.
