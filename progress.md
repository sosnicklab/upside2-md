## Progress Log

- 2026-03-11: Initialized task tracking files for the `ConDiv_symlay` Slurm workflow redesign.
- 2026-03-11: Inspected `ConDiv_symlay/run_remote.sh`, `ConDiv_symlay/training_control.py`, and `ConDiv_symlay/ConDiv_mem.py`; confirmed the current remote workflow is a single Slurm allocation with internal worker fan-out.
- 2026-03-11: Chose a staged design: generate one `sbatch` file per protein simulation, submit a dependent FF-update Slurm job with `afterany`, and let that update job submit the next training round.
- 2026-03-11: Refactored `ConDiv_symlay/ConDiv_mem.py` into reusable minibatch helpers so staged Slurm rounds and the legacy in-process flow share the same FF-update logic.
- 2026-03-11: Added `ConDiv_symlay/slurm_round.py`, `ConDiv_symlay/submit_remote_round.sh`, and `ConDiv_symlay/run_remote_update.sh` to stage per-target job files, submit them, finalize the round, and submit the next round.
- 2026-03-11: Updated `ConDiv_symlay/run_init.sh` and `ConDiv_symlay/README.md` to point to the separate-job Slurm workflow and document the legacy single-allocation mode separately.
- 2026-03-11: Verified Python syntax, shell syntax, and local round staging with a temporary checkpoint copy; could not verify live Slurm submission or full FF finalization without cluster access and generated divergence files.
- 2026-03-11: User corrected the Slurm design: the simulation side must be a single array `sbatch` file per round, with one array task per protein simulation. Refactor in progress.
- 2026-03-11: Refactored `ConDiv_symlay/slurm_round.py` to stage one `simulate_array.sbatch` file per round, submit one array job, and dispatch tasks by `SLURM_ARRAY_TASK_ID` through the round manifest.
- 2026-03-11: Re-verified Python syntax, shell syntax, and local no-submit staging; confirmed the generated simulation script is a single job-array sbatch file.
