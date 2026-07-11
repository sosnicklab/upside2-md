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
- Stage-7 production removes the burn-in protein restraint, minimizes the
  released Hamiltonian, then starts unrestrained production.
- Workflow helper subprocesses must be robust to background/batch invocation
  where inherited stdio descriptors may be closed or invalid.

# Current Task: `*.out` Workflow Failure

- [completed] Scan `example/16.MARTINI/*.out` logs.
- [completed] Identify first failure: all four workflows complete stage-7.0 MD
  and fail during VTF extraction, when `martini_extract_vtf.py` starts with
  invalid inherited stdio and raises `Fatal Python error: init_sys_streams`.
- [completed] Patch the workflow subprocess launcher so helper commands get
  explicit valid stdin/stdout/stderr plumbing.
- [completed] Re-run a representative failed workflow or extraction path under
  redirected/background-style output.
- [completed] Update `progress.md` and `findings.md` with the result.

# Known Errors / Blockers

- No blocker.  The failure was workflow subprocess stdio handling, not
  simulation physics or corrupted checkpoints.
