# Progress Log

## 2026-07-11 Stage-7 Continuation, CGL Timescale, and 1rkl Orientation

- Inspected the failed original `stage_7.3` checkpoints.  Both `1rkl` and
  `1afo` exploded because CGL masses were multiplied by
  `CG_LIPID_MASS_SCALE` on every continuation:
  `84 -> 1.008 -> 0.012096 -> 0.000145 -> 1.74e-6`.
- Patched `py/martini_prepare_system.py` so CGL mass scaling is idempotent:
  the unscaled CGL mass is stored on `/input/mass`, reused on continuation,
  and already repeatedly scaled checkpoints are rejected with a clear error.
- Patched `example/16.MARTINI/run_sim_hybrid.sh` coarse defaults:
  `CG_LIPID_MASS_SCALE=0.3`, `PROD_70_BURNIN_NSTEPS=0`, and
  `STAGE_70_BURNIN_PROTEIN_RESTRAINT_SPRING=0.0`.
- Updated the active CGL-target angular table in
  `parameters/dryMARTINI/dopc.h5` by reversing the target angular axis and
  added matching source-generator provenance in `py/martini_build_tables.py`.
- Root-caused the fresh `1rkl` vertical orientation: the old coarse restrained
  stage-7 burn-in shifted the bilayer relative to the fixed protein by about
  `3-6 A` before unrestrained production.  With the burn-in disabled, the
  protein/bilayer depth offset remains near zero and the helix stays tilted.

## Validation Outputs

- Patched `1rkl` chain:
  `example/16.MARTINI/outputs/martini_1rkl_noburnin_probe_20260711`
  completed `stage_7.0` through `stage_7.3`.
- Patched `1afo` continuation chain:
  `example/16.MARTINI/outputs/martini_1afo_finalfix_cont_20260711`
  completed `stage_7.1` through `stage_7.3` from a valid `stage_7.0`.
- `1rkl` CA principal-axis angle to bilayer normal:
  `35.26 -> 30.27`, `30.18 -> 33.62`, `33.55 -> 32.67`,
  `32.62 -> 31.97` degrees across stages `7.0` to `7.3`.
- `1rkl.stage_7.3`: finite coordinates, fixed box
  `110.452 x 110.452 x 97.551 A`, hbond sum `33.7`, CGL mass `25.2`,
  CGL xy step mean/p95/max `0.0367/0.0780/0.1912 A`.
- `1afo.stage_7.3`: finite coordinates, fixed box
  `91.292 x 91.292 x 149.907 A`, hbond sum `86.1`, CGL mass `25.2`,
  CGL xy step mean/p95/max `0.0404/0.0849/0.1655 A`.
- Full-lipid DOPC COM reference from local full-lipid trajectories:
  `1rkl` `0.0340/0.0682/0.1331 A` and `1afo`
  `0.0339/0.0679/0.1331 A` mean/p95/max per saved frame.  The patched CGL
  clock is therefore not slow; it is slightly faster than full-lipid COM but
  close.
- Hard-error scans of the patched validation `.out` logs found no tracebacks,
  fatal Python errors, runtime errors, overflows, or explosion messages.
- Verification commands passed:
  `python -m py_compile py/martini_build_tables.py py/martini_prepare_system.py`
  and `bash -n example/16.MARTINI/run_sim_hybrid.sh`.

## Current Task: DOPC Tail-Compression Bug Fix

- Installed and validated a corrected DOPC force-field artifact:
  `parameters/dryMARTINI/dopc.h5`.
- Backed up the rejected previous artifact to
  `parameters/dryMARTINI/backup/dopc.before_pair_baseonly_20260710.h5`.
- Patched the table builder, validator, and stage injector so CGL-CGL runtime
  uses the pair-relaxed base spline only:
  `cg_lipid_pair` consumes `compose_vector6d`, with no runtime pair-q delta or
  q-grid payloads.
- Kept the explicit hidden tail-compression state active for the q self PMF and
  SC-CGL.  CGL-target remains an active base spline with no q overlay.
- Verified `parameters/dryMARTINI/dopc.h5` passes
  `_validate_cg_lipid_table_schema`; `_cg_lipid_compaction_metadata_needs_refresh`
  returns `False`; installed `cg_lipid_compaction` has
  `pair_runtime_compaction_overlay=disabled_pair_relaxed_base_only` and no
  runtime pair payloads.
- Ran `python -m py_compile py/martini_build_tables.py
  py/martini_prepare_system_lib.py py/martini_prepare_system.py`.

## Validation Outputs

- Fresh `1rkl` validation:
  `example/16.MARTINI/outputs/martini_1rkl_hybrid_pairbaseonly_full_20260710`
- Fresh `1afo` validation:
  `example/16.MARTINI/outputs/martini_1afo_hybrid_pairbaseonly_full_20260710`

Both systems completed `stage_7.0`, `stage_7.1`, and `stage_7.2`.

## Acceptance Results

- `1rkl` log summaries:
  - `stage_7.0`: final hbond `32.6`, tail hbond mean `30.58`, final Rg `13.1`.
  - `stage_7.1`: final hbond `36.8`, tail hbond mean `37.07`, final Rg `13.4`.
  - `stage_7.2`: final hbond `32.5`, tail hbond mean `35.80`, final Rg `13.7`.
- `1afo` log summaries:
  - `stage_7.0`: final hbond `87.5`, tail hbond mean `83.67`, final Rg `15.3`.
  - `stage_7.1`: final hbond `77.5`, tail hbond mean `76.06`, final Rg `15.0`.
  - `stage_7.2`: final hbond `77.9`, tail hbond mean `80.47`, final Rg `14.8`.
- All prepared/output H5 files retained
  `pair_runtime_compaction_overlay=disabled_pair_relaxed_base_only` and no
  CGL-CGL runtime q payloads.
- Final CGL z ranges remained finite:
  - `1rkl.stage_7.2`: `[-7.18, 12.26] A`
  - `1afo.stage_7.2`: `[-9.10, 12.04] A`
- Protein-proximal CGL q is no longer suppressed:
  - `1rkl.stage_7.2`: q within 20 A of protein `0.458`, far-field q `0.415`.
  - `1afo.stage_7.2`: q within 20 A of protein `0.530`, far-field q `0.478`.
- DSSP-style checks retained secondary structure:
  - `1rkl`: 7.0-start to 7.2-final match `0.968`, final H/E/C `24/0/7`.
  - `1afo`: 7.0-start to 7.2-final match `0.986`, final H/E/C `55/0/17`.

## Conclusion

The bug is fixed for the acceptance surface tested here.  The rejected
runtime pair-q overlay failure did not recur, and the original local
protein-contact q-suppression symptom is corrected in both benchmark systems.

## 2026-07-10 Workflow `.out` Failure

- Inspected `example/16.MARTINI/1rkl.out`, `1afo.out`, `1rkl_full.out`, and
  `1afo_full.out`.
- All four logs completed stage-7.0 MD and failed afterward during VTF
  extraction.  The common error was child Python startup failure:
  `Fatal Python error: init_sys_streams` with `OSError: [Errno 9] Bad file
  descriptor`.
- Root cause: `run_checked(..., log_file=None)` let helper subprocesses inherit
  workflow stdio.  Background/redirected workflow launches can leave inherited
  descriptors invalid for a child Python interpreter.
- Patched `py/martini_prepare_system.py` so helper subprocesses use
  `stdin=subprocess.DEVNULL`, pipe stdout/stderr, and stream captured output
  back to the workflow log.
- Verified the patched no-log path by running `martini_extract_vtf.py` through
  `run_checked` with closed stdin and redirected output; it completed and wrote
  `/tmp/1rkl.stage_7.0.run_checked.vtf`.
- Regenerated the four missing final VTF files from the completed checkpoints:
  `outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf`,
  `outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf`,
  `outputs/martini_1rkl_hybrid_full/1rkl.stage_7.0.vtf`, and
  `outputs/martini_1afo_hybrid_full/1afo.stage_7.0.vtf`.
