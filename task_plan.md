Project Goal
- Fix the `ConDiv_symlay` wrapper scripts so bare `./run_init.sh` and `sbatch run_remote.sh` keep the training run under the current `ConDiv_symlay` folder and do not get redirected by stale environment from other workflows.

Architecture & Key Decisions
- Stop trusting the generic inherited `BASE_DIR` in `ConDiv_symlay` wrappers.
- Keep the default training run under the same `ConDiv_symlay/` directory as the wrapper scripts, using `ConDiv_symlay/test_${PROFILE}`.
- Keep the workflow-local run-directory record file, `ConDiv_symlay/.condiv_current_run_dir`, written by `run_init.sh` after successful initialization.
- Make `run_local.sh` and `run_remote.sh` prefer the recorded run directory, with a discovery fallback for pre-patch initialized runs under the current checkout.
- Preserve the existing Slurm spool-copy path resolution in `run_remote.sh`.

Execution Phases
- [x] Inspect current `ConDiv_symlay` wrappers against the corrected scope.
- [x] Patch `ConDiv_symlay` wrappers to ignore stale generic `BASE_DIR` and keep the run under the current `ConDiv_symlay/` folder.
- [x] Verify locally and update tracking files, findings, progress, lessons, and README.

Known Errors / Blockers
- None.

Review
- Patched [run_init.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_init.sh), [run_local.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_local.sh), [run_remote.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_remote.sh), and [run_validate_rounds.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_validate_rounds.sh) so the training run stays under the current `ConDiv_symlay/` directory and no longer trusts a stale inherited `BASE_DIR`.
- Kept `ConDiv_symlay/.condiv_current_run_dir` recording in `run_init.sh`, with local/Slurm restart scripts preferring that recorded path and `run_remote.sh` still resolving the real checkout from `SLURM_SUBMIT_DIR` / `CONDIV_PROJECT_ROOT`.
- Updated [README.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/README.md) to remove the `BASE_DIR` override guidance for the normal path and document the current-checkout run behavior.
- Verification:
  - `bash -n ConDiv_symlay/run_init.sh ConDiv_symlay/run_local.sh ConDiv_symlay/run_remote.sh ConDiv_symlay/run_validate_rounds.sh`
  - `BASE_DIR=/definitely/wrong/path bash ConDiv_symlay/run_init.sh` initialized under `/Users/yinhan/Documents/upside2-md/ConDiv_symlay/test_dimer3`
  - `BASE_DIR=/definitely/wrong/path bash ConDiv_symlay/run_remote.sh 0` printed:
    - `base dir: /Users/yinhan/Documents/upside2-md/ConDiv_symlay/test_dimer3`
    - `run dir record: /Users/yinhan/Documents/upside2-md/ConDiv_symlay/.condiv_current_run_dir`
    - `checkpoint: /Users/yinhan/Documents/upside2-md/ConDiv_symlay/test_dimer3/initial_checkpoint.pkl`
