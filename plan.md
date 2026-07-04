# Project Goal

Repair the dryMARTINI SC-CGL rigid-DOPC regression by restoring the current
workflow-used correction path to the last committed injector/runtime behavior,
while keeping the interface code clean and avoiding dead flags or test-only
branches.

# Architecture & Key Decisions

- Treat the last committed repository state (`HEAD`, commit `3df4df4`) as the
  behavioral reference for this regression task.
- Define scope from the concrete correction path used by the MARTINI hybrid
  workflows, from direct `HEAD` vs worktree reinjection comparisons, and from
  the fresh Jul 4 stage-7 artifacts the user reported.
- Restore only the SC/target compaction semantics that the committed workflow
  code actually uses.
  Do not reintroduce unused pair mean-field or implicit-response toggles.
- Keep the pair correction split intact:
  `cg_lipid_pair` remains on the committed explicit compaction-state path,
  and `cg_lipid_rotamer_sc` / `cg_lipid_target` must match the committed
  runtime contract in current source.
- Treat the reported Jul 4 `stage_7.0` files as pre-repair artifacts.
  They are valid regression examples, but they do not reflect the current
  injector after the repair landed.
- Verify with controlled reruns from the reported stage-7 sources, not only
  source diffs or prepared-file inspection.

# Execution Phases

- [completed] Compare the reported Jul 4 stage-7 outputs, current prepared
  files, and `HEAD` injector/runtime behavior to isolate the regression scope.
- [completed] Restore the committed SC/target compaction path in Python
  injection and the C++ runtime without reintroducing unrelated dead workflow
  flags.
- [completed] Rebuild and rerun targeted MARTINI stage-7.1 continuations from
  the reported `1rkl` and `1afo` stage-7 sources to confirm the repaired path
  does not push the proteins toward the bad vertical state.
- [in-progress] Update `progress.md` and `findings.md` with the corrected
  lesson, repaired scope, and verification results.

# Known Errors / Blockers

- Full fresh reruns of
  `run_sim_1afo.sh`,
  `run_sim_1rkl.sh`,
  `run_sim_1afo_full.sh`,
  and `run_sim_1rkl_full.sh`
  have not been repeated after this targeted regression repair.
  The current verification is limited to controlled stage-7.1 continuations
  from the user-reported Jul 4 stage-7 sources.

# Review

- The repaired injector/runtime path now matches `HEAD` for the workflow-used
  SC/target correction contract.
- Controlled 2,000-step stage-7.1 continuations from the reported Jul 4
  `1rkl` and `1afo` stage-7 sources do not show the repaired path driving the
  proteins toward a more vertical state than the corresponding implicit-source
  controls.
