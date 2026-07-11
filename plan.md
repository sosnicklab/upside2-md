# Project Goal

Maintain the validated dryMARTINI DOPC hybrid workflow and keep the example
`1rkl`/`1afo` workflows runnable from `example/16.MARTINI`.

# Architecture & Key Decisions

- Keep all hybrid interface interactions active.  Do not disable SC-env,
  BB-env, CGL-target, or CGL-CGL physics to work around failures.
- The accepted DOPC runtime contract is pair-relaxed base-only for CGL-CGL:
  `cg_lipid_pair` consumes only `compose_vector6d`; `cgl_compaction_state`
  remains active for the self PMF and SC-CGL; CGL-target remains base-only with
  no q overlay.
- `parameters/dryMARTINI/dopc.h5` is the shared force-field artifact.  Back up
  old H5 files and overwrite this path; do not create versioned force-field H5
  names.
- Coarse CGL stage-7 production should start from the handed-off/minimized
  state without the restrained burn-in by default.  The restrained coarse
  burn-in can translate the bilayer relative to a fixed protein and create an
  unphysical depth offset before unrestrained production.
- Workflow helper subprocesses must be robust to background/batch invocation
  where inherited stdio descriptors may be closed or invalid.

# Current Task: Stage-7 CGL Continuation, Timescale, and 1rkl Orientation

- [completed] Inspect failed `stage_7.3` outputs for `1rkl` and `1afo`.
- [completed] Identify the explosion root cause: CGL masses were multiplied by
  `CG_LIPID_MASS_SCALE` again on every continuation.
- [completed] Patch continuation metadata so CGL mass scaling is idempotent and
  already repeatedly scaled checkpoints are rejected.
- [completed] Regenerate `stage_7.1` through `stage_7.3` from valid `stage_7.0`
  checkpoints and verify both systems no longer explode.
- [completed] Inspect all `example/16.MARTINI/**/*.out` logs and connect failures
  to either invalid old checkpoints or remaining code defects.
- [completed] Validate CGL transport timescale against full-lipid DOPC COM motion.
- [completed] Diagnose why `1rkl` rotates toward the bilayer normal during
  stable stage-7 production.
- [completed] Patch the minimal orientation/release issue without disabling any
  hybrid interactions.
- [completed] Update the active CGL-target table convention in
  `parameters/dryMARTINI/dopc.h5`, re-run the necessary `1rkl`/`1afo`
  validation stages, and document the result in `progress.md` and
  `findings.md`.

# Known Errors / Blockers

- Existing original `stage_7.1+` checkpoints in `martini_1rkl_hybrid` and
  `martini_1afo_hybrid` are invalid because their CGL masses were already
  repeatedly scaled.  Continue from a valid `stage_7.0` or rerun fresh.
- No current blocker.  The patched coarse workflow keeps `1rkl` tilted through
  `stage_7.3` and both benchmark continuations stay finite.
