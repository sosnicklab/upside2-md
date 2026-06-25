# Project Goal

Find and fix the source-code cause of excessive leaflet separation in the
dryMARTINI Upside hybrid CGL bilayer, without tuning parameters, disabling
hybrid interactions, adding orientational/z restraints, or baking bilayer
target information into the CGL force field.

# Architecture & Key Decisions

- Scope is defined by repository diffs against the master checkout, limited to
  dryMARTINI-related files: `py/martini_*`, `src/martini*`, `src/box.*`,
  `example/16.MARTINI`, and `parameters/dryMARTINI`.
- Diff result: the dryMARTINI interface source is additive relative to master.
  Practical rewrite scope is source/interface code, not generated outputs,
  caches, PDFs, checkpoint files, or backup H5 artifacts.
- Current source/script scope:
  - Python: `py/martini_build_tables.py`, `py/martini_cg_lipid_params.py`,
    `py/martini_extract_vtf.py`, `py/martini_gen_params.py`,
    `py/martini_itp_reader.py`, `py/martini_prepare_system.py`,
    `py/martini_prepare_system_lib.py`, `py/martinize.py`
  - C++: `src/box.cpp`, `src/box.h`, `src/martini.cpp`, `src/martini.h`,
    `src/martini_cg_lipid.cpp`, `src/martini_fix_rigid.cpp`,
    `src/martini_hybrid.cpp`, `src/martini_internal.h`,
    `src/martini_masses.cpp`, `src/martini_potential.cpp`,
    `src/martini_stage_params.cpp`
  - MD scripts: `example/16.MARTINI/build_martini_h5_m1.sh`,
    `example/16.MARTINI/build_martini_h5_slurm.sh`,
    `example/16.MARTINI/run.py`,
    `example/16.MARTINI/run_sim_*.sh`
- Ignore `example/00.AnalysisScripts` completely as a style baseline and edit
  target.
- Preserve all dryMARTINI physical-model contracts: no disabling SC-env,
  BB-env, CGL, or table-based interactions; no arbitrary force-field caps or
  bilayer-derived CGL parameter hacks.
- Prefer direct, clean implementations for active development files
  (`example/16.MARTINI/*`, `src/box.*`, `src/martini*`,
  `py/martini_*.py`) over wrappers, compatibility layers, or disabled legacy
  blocks.
- Refined rewrite target: keep validated physics/table kernels behaviorally
  intact in this pass and clean the dryMARTINI interface layer first:
  configuration parsing, shell handoff, workflow validation, continuation,
  stage preparation, retry policy, and production execution should read as one
  linear workflow.
- Keep git read-only. File edits remain unstaged and uncommitted.
- New integrator cleanup decision: `/Users/yinhan/Documents/upside2-md-master`
  already contains Upside's general multi-step integrator (`mv`),
  `--inner-step`, and `build_integrator_levels`. Do not remove the general
  engine integrator as dryMARTINI-only code. The dryMARTINI workflow should
  simply avoid that path by invoking the standard Verlet integrator (`v`)
  without passing `--inner-step`.
- Revised decision after log triage: the standard Verlet driver path must use
  an active `inner_step` of `1`; only explicit `--integrator mv` should read
  and apply `--inner-step`. Otherwise `--integrator v` still advances and logs
  time as if it were three inner steps, which breaks dryMARTINI restart timing.
- Current investigation must validate a bilayer-only CGL model before judging
  protein-hybrid outputs. Required gates are: two coherent leaflets, CGL centers
  evenly distributed in x-y, CGL orientations consistent with leaflet normals,
  leaflet separation comparable to the initialized physical bilayer, and
  lateral diffusion compared to DOPC after the `14 * 4` coarse-graining factor.
- Parameter twisting is out of scope. Valid fixes must come from setup,
  coordinate persistence, unit conversion, wrapping, force/table application,
  integrator/time accounting, or output/extraction bugs.

# Execution Phases

- [x] Preserve prior refactor/timing results as background context.
- [x] Quantify the reported retained VTF outputs and separate physical CGL
      centers from visualization rods.
- [x] Trace CGL-only and hybrid setup/runtime code for source-level causes of
      excessive leaflet distance.
- [x] Implement the smallest source-level fix, if a root cause is confirmed.
- [x] Resolve CGL table metadata validation for bilayer-only spline-node
      injection.
- [x] Validate with a bilayer-only CGL run before interpreting hybrid protein
      systems.
- [x] Document evidence, remaining risks, and verification results.

# Known Errors / Blockers

- Current user-reported failure: retained hybrid outputs
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf` and
  `example/16.MARTINI/outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf` appear
  to have a larger-than-expected CGL leaflet distance.
- Bilayer-only validation now runs but fails the geometry target. With the
  current installed one-particle CGL pair table, the DOPC template CGL COM
  separation starts at `12.77 A`, while a short CGL-only run expands toward
  `33 A`. Energy scans show the current pair table favors wider center
  separation, with much lower total pair energy near `24 A` than near the
  source COM separation.

# Review

- Diff scope showed the dryMARTINI source/interface files are additive relative
  to `/Users/yinhan/Documents/upside2-md-master`; generated outputs, caches,
  PDFs, checkpoints, and backup H5 files were excluded from rewrite scope.
- Master style baseline was taken from human-written `py/extract_vtf.py`,
  `py/upside_config.py`, `src/spring.cpp`, and `src/deriv_engine.h`, not from
  `example/00.AnalysisScripts`.
- Refactored `py/martini_prepare_system.py` so environment config, CLI parsing,
  validation, stage file naming, stage-6 retry, optional pre-7.0 equilibration,
  and stage-7 production are named, linear helpers instead of one patched
  command body.
- Refactored `src/martini_hybrid.cpp` to centralize repeated dynamic fixed
  atom / rigid group / z-hold stage logic while preserving the prior
  registration and refresh enabled-state behavior.
- Refactored `example/16.MARTINI/run_sim_hybrid.sh` to move MARTINI parameter
  library detection/generation into a named function before workflow handoff.
- Verification passed: Python compile for dryMARTINI Python files, shell
  syntax for MARTINI scripts, `martini_prepare_system.py run-hybrid-workflow
  --help`, `make -C obj upside -j2`, and `git diff --check` for touched files.
- Follow-up audit removed ignored/versioned H5 development backups from
  `parameters/dryMARTINI` and the stale temp H5 file from
  `example/16.MARTINI/outputs`. Remaining dryMARTINI H5 names are canonical:
  `dopc.h5`, `particle.h5`, `sidechain.h5`, and workflow
  `hybrid_mapping.h5` outputs.
- Removed inactive runtime flags from the dryMARTINI interface. Hybrid, stage,
  fixed-rigid, barostat, and CGL GLE configuration now uses group/config
  presence and stage activation directly rather than `enable`/`enabled`
  attributes. The only remaining `enabled` strings in scoped dryMARTINI files
  are stale force-field provenance metadata names used to reject invalid
  bilayer-trained tables.
- Final verification passed: dryMARTINI Python compile, MARTINI shell syntax,
  workflow help rendering, inactive-flag scan, canonical H5 filename audit,
  `git diff --check`, and `make -C obj upside -j2`.
- Integrator cleanup removed the dryMARTINI workflow `mv`/`--inner-step` path.
  `run_md_stage` now constructs `--integrator v` with no `--inner-step`; frame
  and restart timing use the requested `dt` directly. Master comparison showed
  the general C++ multi-step integrator is upstream master behavior, so the
  multi-step APIs remain in place.
- Follow-up log triage found that plain `--integrator v` still inherited an
  active `inner_step=3` in `src/main.cpp`, so stage 6.0 output time became
  `15` instead of `5`. The driver now initializes `inner_step=1` and only
  applies the `--inner-step` value when `--integrator mv` is explicitly used.
  Short isolated checks for `1rkl`, `1rkl_full`, `1afo`, and `1afo_full`
  completed through `=== Workflow Complete ===`.
