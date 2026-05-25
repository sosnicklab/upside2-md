# 2026-05-24 MARTINI Secondary-Structure Divergence

## Project Goal
- Find why protein secondary-structure results differ between the full-scale
  lipid workflow and the single-particle lipid workflow under
  `example/16.MARTINI`.
- Identify actual bugs, fix them with minimal code changes, and verify the
  workflows use consistent protein/environment physics.

## Architecture & Key Decisions
- Treat the single-particle lipid model as a CG representation of the same
  dry-MARTINI lipid environment, not a separate protein force field.
- Keep SC-env, BB-env, CGL-CGL, CGL-SC, and CGL-target interactions active.
- Compare workflow scripts, generated inputs, and existing outputs before
  editing runtime code.
- Do not change shared physics parameters unless diagnostics show a direct bug.
- Root-cause decision: coarse CGL-SC was applied as a standalone CB potential,
  while full-resolution DOPC-SC entered the sidechain rotamer one-body energy.
  The physical fix is to feed the dry-MARTINI-derived CGL-SC spline into the
  rotamer solver as a one-body term, preserving active CGL-SC forces through
  rotamer probabilities instead of bypassing sidechain conformation selection.
- Runtime node decision: the rotamer-coupled CGL-SC term is named
  `cg_lipid_rotamer_sc` to avoid the deriv-engine prefix collision with the
  legacy standalone `cg_lipid_sc` node type.

## Execution Phases
- [x] Phase 1: Inventory relevant full-scale and single-particle workflows and
      output artifacts.
- [x] Phase 2: Compare generated HDF5 inputs and logs for force-field/node
      differences that can affect protein secondary structure.
- [x] Phase 3: Isolate root cause and implement the smallest correct fix.
- [x] Phase 4: Run focused verification and document residual risks.

## Known Errors / Blockers
- None yet.

## Review
- Root cause: full lipid mode exposed explicit DOPC beads to
  `martini_sc_table_1body`, but single-particle lipid mode only exposed ions
  there; its CGL-SC table was evaluated outside the rotamer one-body solver.
- Fix: inject `cg_lipid_rotamer_sc` as a rotamer one-body coordinate node using
  the existing dry-MARTINI CGL-SC spline, append it to `rotamer` arguments, and
  refresh generated `cg_lipid_*` nodes when reinjecting a prepared file.
- Verification: Python syntax check passed, C++ build passed, copied 1RKL
  injection produced 117/117 matched rotamer rows and 282 CGL targets, and a
  one-step `obj/upside` smoke test ran with nonzero `rotamer_1body_energy2`
  entries from the new CGL-SC one-body channel.
