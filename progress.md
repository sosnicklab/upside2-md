# Progress Log

## Current Task: CGL Leaflet-Distance Source Investigation

- Actions taken:
  - Reviewed existing `plan.md`, `findings.md`, and `progress.md`.
  - Reframed the active plan around source-code causes of excessive CGL
    leaflet separation and bilayer-only validation.
  - Quantified retained 1RKL/1AFO hybrid outputs. Physical CGL center leaflet
    separation is about 30 A; VTF synthetic rod endpoints inflate apparent
    head/head separation to about 57 A.
  - Patched the shared C++ PBC minimum-image helper to use rounded box-image
    reduction instead of a single add/subtract.
  - Added `prepare --mode bilayer` so a CGL-only bilayer can be generated
    without faking a protein.
  - Generated a 72-CGL bilayer-only input under
    `example/16.MARTINI/outputs/codex_cgl_bilayer_pbc`.
  - Fixed the CGL table validator default to match the installed `Tavg=25.0`
    H5 contract.
  - Patched the generic MARTINI potential to accept zero optimized pair rows
    as a no-op, which is required for pure CGL-only spline-node systems.
  - Ran a short CGL-only bilayer validation and a center-separation energy
    scan. The validation still fails geometry because the current installed
    one-particle CGL pair table favors a wider CGL-center separation than the
    source DOPC COM geometry.
- Files modified:
  - `src/box.h`
  - `src/martini_potential.cpp`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Test results:
  - Python compile passed for dryMARTINI Python files.
  - `make -C obj upside -j2` passed.
  - Bilayer-only prepare/conversion passed.
  - Bilayer-only CGL spline-node injection passed.
  - Short bilayer-only CGL MD ran, but geometry failed: source DOPC CGL COM
    separation `12.77 A` expanded toward `33 A`.
  - Energy scan confirmed the current pair table is lower in total energy near
    `24 A` than near the source COM separation.

## Previous Task Summary: dryMARTINI Interface Refactor

- Actions taken:
  - Compared dryMARTINI-related Python, C++, parameter, and MARTINI example
    paths against `/Users/yinhan/Documents/upside2-md-master`.
  - Confirmed the source/interface dryMARTINI code is additive relative to
    master; ignored generated outputs, caches, PDFs, checkpoints, and backup
    H5 files for rewrite scope.
  - Read human-written master Python and C++ style baselines, excluding
    `example/00.AnalysisScripts`.
  - Reworked the hybrid workflow Python command into named, linear helpers for
    config parsing, validation, stage files, stage 6, pre-7.0, and stage 7.
  - Centralized repeated hybrid C++ stage-fixing logic in one helper.
  - Moved MARTINI parameter-file detection/generation in the shell wrapper into
    a named function.
  - Removed ignored/versioned H5 development backups from
    `parameters/dryMARTINI` and removed the stale temp H5 file under
    `example/16.MARTINI/outputs`.
  - Removed inactive runtime flags from hybrid control, stage parameters,
    fixed-rigid setup, barostat setup, and CGL GLE setup. These paths now use
    group/config presence and stage activation directly.
  - Removed the dryMARTINI workflow multi-step integrator invocation. The
    workflow now invokes the standard Verlet integrator with no `--inner-step`.
    Master comparison confirmed the general C++ multi-step integrator is not
    dryMARTINI-only, so it remains in place.
  - Fixed the reported `example/16.MARTINI/*.out` failures. All four logs
    failed after stage 6.0 because the driver still used `inner_step=3` for
    `--integrator v`, producing output time `15` where the workflow expected
    `5`. `src/main.cpp` now uses `inner_step=1` unless `--integrator mv` is
    explicitly selected.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `findings.md`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `src/box.cpp`
  - `src/box.h`
  - `src/main.cpp`
  - `src/martini_cg_lipid.cpp`
  - `src/martini_fix_rigid.cpp`
  - `src/martini_hybrid.cpp`
  - `src/martini_internal.h`
  - `src/martini_potential.cpp`
  - `src/martini_stage_params.cpp`
  - `example/16.MARTINI/run_sim_hybrid.sh`
- Verification:
  - Python compile passed for the dryMARTINI Python files.
  - Shell syntax checks passed for MARTINI build/run scripts.
  - `martini_prepare_system.py run-hybrid-workflow --help` rendered
    successfully.
  - `make -C obj upside -j2` passed.
  - `git diff --check` passed for touched source/script files.
  - H5 filename audit now finds only canonical force-field H5 names and
    workflow `hybrid_mapping.h5` outputs.
  - Inactive-flag scan finds no hybrid interaction disable flags. Remaining
    `enable` names are NPT setup controls, `--disable-recentering`, predicate
    names, and stale-table provenance metadata (`force_match_enabled`,
    `ibi_enabled`).
  - Command-construction check shows `--integrator v`, no `--inner-step`, and
    direct `dt` timing (`frame_interval=0.01` for 5 frames at `dt=0.002`).
  - Direct binary smoke test on a copied checkpoint wrote final time `0.1` for
    10 steps at `dt=0.01` with `--integrator v`.
  - Short isolated workflow checks completed for `1rkl`, `1rkl_full`, `1afo`,
    and `1afo_full` under `example/16.MARTINI/outputs/codex_verify_*`.
