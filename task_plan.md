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

## Follow-up Goal

Fix the `ConDiv_symlay` worker post-processing crash so completed simulation outputs can still produce `*.divergence.pkl`, allowing the update job to finalize `epoch_00_minibatch_00` and submit later rounds.

## Revised Decisions

- Treat the current failure as a worker-side post-processing crash after the `upside` simulation completes, not as a Slurm dependency/orchestration bug.
- Keep the existing round/update submission flow unchanged; fix the worker so `epoch_00_minibatch_00` can be rerun from `initial_checkpoint.pkl`.
- Isolate the native `upside_engine` contrast evaluation in a short-lived helper subprocess that persists its result before process teardown, so a native cleanup crash cannot take down the main worker before `*.divergence.pkl` is written.
- Some downloaded targets also crash when a large contrast batch is evaluated in one helper call, even after teardown is bypassed. Handle that by recursively bisecting failing frame batches and concatenating the per-batch derivative tensors.
- Do not force users to rerun completed simulations when only post-processing failed. The finalize/update path should recover missing `*.divergence.pkl` files from existing `.run.*.h5` outputs before it attempts the FF update.
- Treat empty inner-gradient blocks from recovered divergences as zero gradients with the same shape as the checkpoint parameter tensors, so the pair regularizer and optimizer update remain shape-consistent.

## Follow-up Phases

- [x] Confirm where the worker writes `*.divergence.pkl` and which native step runs immediately before that write.
- [x] Implement the isolated divergence-evaluation helper in `ConDiv_symlay/ConDiv_mem.py`.
- [x] Verify syntax and CLI plumbing for the new helper entrypoint.

## Follow-up Blockers

- The native `double free or corruption (!prev)` crash happens on the cluster after simulation completion, so local verification is limited to syntax and code-path checks unless a representative run directory is available locally.

## Follow-up Review

- Added an internal `compute-divergence` helper mode in `ConDiv_symlay/ConDiv_mem.py` that computes the native `upside_engine` contrast in a short-lived subprocess, flushes the result pickle, and exits via `os._exit(0)` to avoid crashing the main worker during native teardown.
- Updated `py/upside_engine.py` so the helper subprocess can explicitly skip `free_deriv_engine`, and extended the helper path to recursively bisect failing frame batches and concatenate their results.
- Switched `main_worker()` to call the isolated helper for the contrast step, added update-time recovery for missing divergence files, and normalized empty inner-gradient blocks against checkpoint parameter shapes while keeping the existing RMSD/COM calculations, divergence format, checkpoint flow, and Slurm orchestration unchanged.
- Verification completed:
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile ConDiv_symlay/ConDiv_mem.py py/upside_engine.py`
  - replayed the downloaded `ConDiv_symlay/test_dimer3/epoch_00_minibatch_00` outputs through `ConDiv_mem._compute_divergence_isolated(...)` for all 15 targets; all returned 300-frame derivative tensors successfully
  - rebuilt all 15 missing divergence files in the downloaded `epoch_00_minibatch_00` round and confirmed `finalize_minibatch_from_outputs()` completed successfully against that recovered data
  - created a fresh temporary copy of `ConDiv_symlay/test_dimer3`, deleted all recovered divergence files, patched the copied checkpoint/manifest to local paths, and ran `python3 ConDiv_symlay/slurm_round.py finalize-round --round-manifest <temp manifest>` successfully end-to-end; it recovered all 15 divergence files, wrote `epoch_00_minibatch_00/checkpoint.pkl`, and updated training progress/status outputs

## Current Correction

- The user confirmed `ConDiv_symlay/slurm-46899622.out` is from the current code, so the previous follow-up review is not sufficient evidence that the bug is fixed on the active code path.
- Local reproduction also showed the current helper path can still crash child Python processes during divergence recovery, which matches the user-visible macOS popup failures even when the parent process retries and eventually succeeds.

## Current Goal

Eliminate the native membrane-gradient corruption that crashes `ConDiv_symlay` divergence recovery, while preserving bitwise-identical simulation results before and after the fix.

## Current Decisions

- The remaining crash is a native bug in membrane parameter-derivative evaluation, not a simulation-integrator bug and not primarily a Python orchestration bug.
- The required native fix is to stop writing right-side spline derivative contributions past the end of the derivative buffer in `src/membrane_potential.cpp`.
- Keep the Python-side helper hardening (`divergence_helper.py`, single-threaded helper env, lighter `Upside` metadata reads) because it reduces unrelated teardown noise, but it is not the root-cause fix.
- Restore the normal helper batch limit to `64` after the native fix; the analysis path should now succeed without relying on artificially small batches.
- Preserve simulation behavior exactly. Any change in `src/membrane_potential.cpp` is acceptable only if seeded old-vs-new simulation outputs remain identical.

## Current Phases

- [x] Reproduce the current helper crash path directly from the downloaded `epoch_00_minibatch_00` artifacts.
- [x] Identify the concrete native root cause in membrane derivative evaluation.
- [x] Patch the native out-of-bounds writes and rebuild the affected binaries.
- [x] Verify stress-target recovery and end-to-end `finalize-round` on a fresh local copy of the downloaded round.
- [x] Verify seeded old-vs-new simulation parity for the `src/membrane_potential.cpp` change.

## Current Review

- Root cause: `MembraneCBPotential::get_param_deriv()` and `MembraneSurfCBPotential::get_param_deriv()` in `src/membrane_potential.cpp` both allocated `deriv(np)` and then wrote right-side spline contributions at `shift + np`, corrupting memory during membrane-gradient extraction.
- Native fix applied: remove the extra `+ np` offset in those two derivative-write paths so both left- and right-side spline terms index the allocated `np`-length buffer correctly.
- Analysis verification:
  - previously unstable targets `2oar`, `3ukm`, and `5vb2` now complete divergence recovery with zero helper retries at `CONDIV_DIVERGENCE_BATCH_LIMIT=64`
  - a fresh `python3 ConDiv_symlay/slurm_round.py finalize-round --round-manifest <patched temp manifest>` replay on a clean temp copy of `ConDiv_symlay/test_dimer3` rebuilt all 15 divergence files, wrote `checkpoint.pkl`, and reported no recovery failures
- Simulation-parity verification:
  - built a temporary pre-patch `upside` binary and ran old/new binaries on identical stripped copies of `2oar.run.1.h5` with the same seed
  - compared all shared `/output` datasets and found exact matches for every dataset, including `pos`, `potential`, `kinetic`, and `hbond`
