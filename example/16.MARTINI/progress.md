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

## 2026-02-13
- Investigated scan sweep mismatch between `scan_lipid_diffusion.py` and `scan_water_diffusion.py`.
- Confirmed lipid previously used a narrow hardcoded grid (`T=0.75..1.00`, `tau=[1,2,4,8]`) while water used a broad scan grid.
- Updated `scan_lipid_diffusion.py` with minimal parity changes:
  - Added argparse-based entrypoint style matching water scan script.
  - Changed `build_temperature_range()` to `np.arange(0.600, 1.001, 0.02)`.
  - Changed `build_tau_range()` to merged low/high/special tau logic used by water scan.
  - Kept bilayer-specific runtime script generation and diffusion-analysis logic unchanged.
- Verification run:
  - Command: `source ../../.venv/bin/activate && source ../../source.sh && python3 scan_lipid_diffusion.py --base-dir lipid_diffusion_check`
  - Result: `Generated 651 simulation scripts.`
  - Verified `lipid_diffusion_check/tasks.txt` has 651 lines and `lipid_diffusion_check/run_scan.slurm` contains `#SBATCH --array=1-651`.
