# Findings

## 2026-04-14 (Current Workflow Target)
- `hybrid-interface-sweep/` must be a bilayer-only calibration workflow.
- The requested fix is to use softened potentials for both LJ and Coulomb in production-stage `7.0`.
- The requested change belongs only in `hybrid-interface-sweep/`.
  - `example/16.MARTINI/` is not part of this correction.

## 2026-04-14 (Scientific Direction)
- The current scalar `interaction_scale` workflow is the wrong control surface for this task.
- The user’s concern is that plain coefficient scaling leaves the LJ and Coulomb singular shapes intact near `r -> 0`.
- The restored sweep should therefore operate on the potential shape controls already supported by the engine:
  - `lj_soften`
  - `lj_soften_alpha`
  - `coulomb_soften`
  - `slater_alpha`

## 2026-04-14 (Stage-7 Rewrite Surface)
- The bilayer prep path already writes the required MARTINI node for production-stage patching:
  - `/input/potential/martini_potential`
- Rewriting only the staged `7.0` attrs is sufficient to change the production nonbonded forces seen by the bilayer-only diffusion measurement.
- The corrected workflow should record the applied production softening settings in task results so the analysis tree remains self-describing.

## 2026-04-14 (Sweep Surface And Compatibility)
- The workflow should sweep two non-negative controls:
  - `lj_alpha`
  - `slater_alpha`
- The earlier softening workflow grid remains a reasonable default restoration target:
  - `lj_alpha = 0.0, 0.025, 0.05, 0.10, 0.20`
  - `slater_alpha = 0.0, 0.25, 0.50, 1.00, 2.00`
- Scalar-manifest directories must be rejected explicitly by a schema bump so stale runs are not mixed into the restored softening workflow.

## 2026-04-14 (Analysis Surface)
- The existing bilayer-only diffusion analysis remains the right observable family:
  - `PO4` lateral diffusion relative to bilayer COM,
  - physical-unit diffusion conversion via the configured `ps/step` assumption,
  - reciprocal-diffusion viscosity proxy,
  - optional target-diffusion ranking.
- Only the grouping axis and recommendation surface need to change:
  - from `interaction_scale`
  - back to `(lj_alpha, slater_alpha)`.

## 2026-04-14 (Runtime Execution Surface)
- The workflow was still failing locally even after the semantic restoration because the preparation helpers launched bare `python3` subprocesses.
- Since the wrappers and generated Slurm scripts already choose a specific interpreter, the workflow subprocesses must inherit that same interpreter.
- Using `sys.executable` inside `workflow.py` keeps:
  - local wrapper runs,
  - direct workflow invocations,
  - and Slurm array tasks
  on one consistent Python environment with the required packages.

## 2026-04-14 (Downloaded Softening Analysis Rerun)
- The local repo copy currently contains the downloaded `hybrid-interface-sweep/analysis/` tree only:
  - task-level analysis JSON files,
  - assembled CSV/JSON summaries,
  - Slurm logs and manifests.
- It does not contain the corresponding full sweep task tree locally, so this rerun is a validation/review of the downloaded analysis artifacts rather than a fresh extraction from raw `stage_7.0.up` files.
- The downloaded softening analysis is numerically complete:
  - `75 / 75` analysis task JSON files succeeded,
  - `0` failed tasks,
  - all `25` conditions have full `3 / 3` replicate coverage.
- Internal consistency is strong across the full sweep:
  - every task uses `160` post-burn-in frames,
  - every task uses `102` `PO4` beads,
  - all diffusion values are finite,
  - task-level `PO4` fit quality spans `R^2 = 0.9780 -> 0.9927`.
- The saved recommendation from the downloaded artifacts is supported by the condition table:
  - best mean diffusion: `lj_alpha = 0.1`, `slater_alpha = 0.5`, `D = 0.8621 um^2/s`,
  - but that branch is somewhat noisier than the nearby alternatives (`CV = 0.082`).
- Two practical neighboring branches remain important for interpretation:
  - lower-perturbation high-fluidity branch: `(lj_alpha, slater_alpha) = (0.0, 1.0)` with `D = 0.8584 um^2/s`, `CV = 0.047`,
  - highest-stability near-top branch: `(0.1, 2.0)` with `D = 0.8510 um^2/s`, `CV = 0.013`, `min R^2 = 0.9888`.
- Relative to the same provisional `40 ps/step` target proxy from `hybrid_timescale.md`:
  - target at workflow temperature `T = 0.8647` remains about `2.892 um^2/s`,
  - the downloaded softening sweep reaches only `0.733 -> 0.862 um^2/s`,
  - so no tested condition demonstrates a target match.

## Lessons
- When the user points out that scaling leaves a singular pair-potential shape unchanged, restore the shape-control surface instead of trying to force the scalar workflow to fit.
- When a workflow semantics change invalidates old run trees, bump the manifest/result schemas immediately so reuse fails loudly.
- When a workflow script is launched from a managed environment, subprocess Python calls must use `sys.executable` rather than bare `python3`.
- When the user provides only a downloaded analysis tree, state explicitly whether the rerun is a fresh recomputation from raw checkpoints or a validation/review of the downloaded artifacts.
