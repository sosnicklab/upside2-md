## Project Goal

Fix the remote training workflow for `ConDiv_symlay` so each Upside simulation runs as a separate Slurm job instead of as a single multi-CPU job, and add a Slurm-executable driver that performs the FF update after all simulations finish and submits the next round of training.

## Architecture & Key Decisions

- Keep the existing training workflow structure and touch only the orchestration layer.
- Replace the current "one Slurm allocation runs all simulations" behavior with one generated Slurm array script per minibatch round.
- Each array element runs exactly one protein simulation, so scheduling is per simulation while the round still uses one `sbatch` file.
- Keep the fixed Slurm update driver script that runs after the whole array finishes, performs the FF update from completed outputs, and submits the next training round.
- Use `afterany` on the array job for the update job so one failed simulation does not deadlock a whole minibatch; the Python finalizer will decide whether enough successful outputs exist to continue.
- Preserve current training/config interfaces where possible; prefer generating wrapper scripts over changing simulation internals.

## Execution Phases

- [x] Inspect existing `ConDiv_symlay` orchestration, Slurm generation, and training loop entry points.
- [x] Identify why the current remote run collapses into a single job despite requesting many CPUs.
- [x] Refactor the staged workflow to emit one simulation array `sbatch` script per round.
- [x] Update FF update / next-round submission flow to depend on the array job.
- [x] Verify array staging behavior and update docs.

## Known Errors / Blockers

- Current remote training round appears to request many CPUs inside one Slurm allocation, so only one simulation is running instead of independent jobs.
- Cluster-specific details may be encoded in existing scripts; those must be preserved unless they cause the orchestration bug.

## Review

- Revised the separate-job workflow to use one Slurm array script per round:
  - `ConDiv_symlay/submit_remote_round.sh` still submits the next round.
  - `ConDiv_symlay/slurm_round.py` now stages one `simulate_array.sbatch` file with one array task per protein and submits one array job instead of one `sbatch` per protein.
  - `ConDiv_symlay/run_remote_update.sh` remains the single dependent FF-update job.
- Added `run-array-task` handling in `ConDiv_symlay/slurm_round.py` so each array task resolves its protein target from the round manifest and runs the matching worker spec.
- Kept the reusable minibatch finalization split in `ConDiv_symlay/ConDiv_mem.py`; only the Slurm orchestration changed.
- Updated `ConDiv_symlay/README.md` to describe the array workflow and the single generated `simulate_array.sbatch` file.
- Verification completed:
  - `python3 -m py_compile ConDiv_symlay/ConDiv_mem.py ConDiv_symlay/slurm_round.py ConDiv_symlay/training_control.py`
  - `bash -n ConDiv_symlay/submit_remote_round.sh ConDiv_symlay/run_remote_update.sh`
  - local dry-run staging with `python3 ConDiv_symlay/slurm_round.py submit-round --base-dir <tmpdir> --run-steps 2 --no-submit`
  - verified generated `simulate_array.sbatch` contains `#SBATCH --array=0-14` and dispatches by `SLURM_ARRAY_TASK_ID`
- Remaining limitation:
  - I still could not run real `sbatch` submissions or a full FF-update finalize on this machine, because there is no live Slurm cluster here and the bundled sample round does not include persisted `*.divergence.pkl` outputs.
