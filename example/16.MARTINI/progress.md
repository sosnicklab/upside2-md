# Progress Log

## 2026-02-11
- Initialized required workflow files: `task_plan.md`, `findings.md`, `progress.md`.
- Prepared to reproduce failure for `run_sim_1ubq.sh`.
- Ran `bash run_sim_1ubq.sh`; initial failure was environment invocation from wrong path (`.venv/bin/activate` not found when pre-sourcing in current dir). Switched to running script directly in-place since script already sources `../../source.sh` and `../../.venv/bin/activate`.
- Reproduced workflow blocker at stage handoff:
  `set_initial_position.py` failed with `IndexError: Index (-1) out of range for empty dimension` when reading `f['/output/pos'][-1, ...]`.
- Modified `set_initial_position.py`:
  - Added guard to use `/output/pos` only when it exists and has at least one frame.
  - Added fallback to `/input/pos` for empty-output files (e.g., minimize-only files).
- Re-ran workflow with short verification run:
  `NPT_EQUIL_STEPS=20 NPT_PROD_STEPS=20 FRAME_INTERVAL=5 bash run_sim_1ubq.sh`.
- Verification results:
  - `outputs/martini_test/logs/npt_equilibration.log` shows completion with checkpoint write.
  - `outputs/martini_test/logs/npt_production.log` shows simulation startup (`0 / 20 ...` progress line).
  - Checkpoints produced: `1ubq.prepared.up`, `1ubq.minimized.up`, `1ubq.npt_equil.up`, `1ubq.npt_equil_reduced.up`, `1ubq.npt_prod.up`.
