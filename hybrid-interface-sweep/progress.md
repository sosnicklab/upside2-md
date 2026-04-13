# Progress Log

## 2026-04-13 (Initialization)
- Created a dedicated root-level workflow folder for sweeping the hybrid `PROTEIN_ENV_INTERFACE_SCALE`.
- Recorded the implementation goal, design choices, and execution phases in the local tracker files before code changes.

## 2026-04-13 (Implementation)
- Added `hybrid-interface-sweep/workflow.py` with:
  - manifest generation,
  - task execution through `example/16.MARTINI/run_sim_1rkl.sh`,
  - task-local `RUN_DIR` allocation,
  - result aggregation,
  - Slurm array and collector staging.
- Added:
  - `hybrid-interface-sweep/run_local.sh`
  - `hybrid-interface-sweep/submit_remote_round.sh`
  - `hybrid-interface-sweep/README.md`

## 2026-04-13 (Verification)
- Ran `python3 -m py_compile hybrid-interface-sweep/workflow.py`.
- Ran `bash -n hybrid-interface-sweep/run_local.sh`.
- Ran `bash -n hybrid-interface-sweep/submit_remote_round.sh`.
- Staged a no-submit Slurm run under `/tmp/hybrid_interface_sweep_smoke`.
- Verified the staged Slurm output includes:
  - `slurm/round_manifest.json`
  - `slurm/train_array.sbatch`
  - `slurm/collect_results.sbatch`
- Ran a reduced one-task local smoke run with captured hybrid env overrides:
  - `MIN_60_MAX_ITER=1`
  - `MIN_61_MAX_ITER=1`
  - `EQ_62_NSTEPS=1`
  - `EQ_63_NSTEPS=1`
  - `EQ_64_NSTEPS=1`
  - `EQ_65_NSTEPS=1`
  - `EQ_66_NSTEPS=1`
  - `PROD_70_NSTEPS=1`
  - `EQ_FRAME_STEPS=1`
  - `PROD_FRAME_STEPS=1`
- Verified the smoke task:
  - completed successfully,
  - wrote `results/tasks/scale0p85_r01.json`,
  - produced task-local `run/checkpoints/1rkl.stage_7.0.up`,
  - assembled `1` successful task and `1` completed scale under `assembled/summary.json`.
