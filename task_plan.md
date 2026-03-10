Project Goal
- Fix `ConDiv_symlay` Slurm path/run-directory resolution so `sbatch run_remote.sh` uses the intended checkout and matching initialized run directory instead of a stale inherited path.

Architecture & Key Decisions
- Keep `SCRIPT_DIR` / `PROJECT_ROOT` resolution based on the submitting checkout as before.
- Add a workflow-local run-directory record file, `ConDiv_symlay/.condiv_current_run_dir`, written by `run_init.sh` after successful initialization.
- Make zero-argument `run_remote.sh` prefer that recorded run directory over inherited `BASE_DIR` values, because the intended workflow is `./run_init.sh` then `sbatch run_remote.sh`.
- Preserve explicit custom-run support through `run_init.sh`: if the user initializes a custom `BASE_DIR`, the recorded run directory updates to that path and later zero-argument Slurm restarts follow it.

Execution Phases
- [x] Inspect the wrapper path/bootstrap logic against the reported Slurm log.
- [x] Patch `run_init.sh` and `run_remote.sh` to record and reuse the active run directory from the current checkout.
- [x] Verify that a stale inherited `BASE_DIR` no longer overrides the recorded current run directory.

Known Errors / Blockers
- None.

Review
- Patched [run_init.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_init.sh) to normalize the resolved run directory and record it in `ConDiv_symlay/.condiv_current_run_dir` after a successful init.
- Patched [run_remote.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_remote.sh) so the recorded run directory now wins for the normal zero-argument restart path.
- Updated [README.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/README.md) to document that behavior.
- Verification:
  - `bash -n ConDiv_symlay/run_init.sh ConDiv_symlay/run_remote.sh`
  - local wrapper-resolution probe with `BASE_DIR=/definitely/wrong/path` and a recorded run dir pointing at `/tmp/condiv_symlay_inner_baseline_smoke_fast`
  - the wrapper printed:
    - `base dir: /tmp/condiv_symlay_inner_baseline_smoke_fast`
    - `run dir record: .../ConDiv_symlay/.condiv_current_run_dir`
    - `checkpoint: /tmp/condiv_symlay_inner_baseline_smoke_fast/initial_checkpoint.pkl`
