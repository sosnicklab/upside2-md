Project Goal
- Hard-code `ConDiv_symlay` Slurm training to use all 48 CPUs in a fixed `6 workers x 8 CPUs` layout while preserving the Train-style one-worker-per-protein launch model.

Architecture & Key Decisions
- Match the `Train/ConDiv.py` worker model instead of the temporary replica-slot model:
  - one Python worker per protein
  - each worker runs a local `upside` replica-exchange bundle
  - one optimizer update after the entire minibatch completes
- Under Slurm, dispatch each protein worker through `srun --ntasks=1 --cpus-per-task=<omp_threads>` instead of running workers as plain local subprocesses.
- Hard-code the Slurm CPU layout in `run_remote.sh`:
  - `#SBATCH --ntasks-per-node=6`
  - `#SBATCH --cpus-per-task=8`
  - `CONDIV_N_REPLICA=8`
  - `CONDIV_OMP_THREADS=8`
  - `CONDIV_MAX_PARALLEL_WORKERS=6`
- Fail fast if the live Slurm allocation differs from the required `48` CPUs / `6` task slots, so the workflow does not silently drift to another CPU layout.
- Preserve the existing `ConDiv_symlay` run-directory fixes:
  - workflow-local default run dir
  - `.condiv_current_run_dir`
  - Slurm spool-copy path resolution

Execution Phases
- [x] Patch `ConDiv_symlay` Slurm wrapper to a fixed 48-CPU layout.
- [x] Update docs and task tracking for the hard-coded CPU model.
- [x] Verify syntax and resolved fixed settings.

Known Errors / Blockers
- None.

Review
- Kept the `ConDiv_symlay` Slurm worker model aligned with the reference `/Users/yinhan/Documents/Train` workflow:
  - one `srun --ntasks=1 --cpus-per-task=<omp_threads>` worker step per protein
  - each worker launches one local `upside` replica bundle
  - FF update still happens after the full minibatch
- Changed [run_remote.sh](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/run_remote.sh) so the Slurm CPU layout is now hard-coded to:
  - `#SBATCH --ntasks-per-node=6`
  - `#SBATCH --cpus-per-task=8`
  - `CONDIV_N_REPLICA=8`
  - `CONDIV_OMP_THREADS=8`
  - `CONDIV_MAX_PARALLEL_WORKERS=6`
  - with fail-fast checks if the allocation is not exactly `48` CPUs / `6` task slots
- Updated [README.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/README.md) to document the hard-coded `6 x 8 = 48` layout.
- Verification:
  - `bash -n ConDiv_symlay/run_remote.sh`
  - `python3 -m py_compile ConDiv_symlay/ConDiv_mem.py`
  - fake-Slurm import probe under the project venv/source bootstrap returned:
    - `n_replica=8`
    - `omp_threads=8`
    - `max_parallel_workers=0` in the Python default config, with the wrapper now overriding that to a fixed `6`
  - direct `_run_worker_subprocess(...)` command capture returned:
    - `['srun', '--exclusive', '--nodes=1', '--ntasks=1', '--cpus-per-task=8', '--slurmd-debug=0', '--output=/tmp/1orq.output_worker', 'python3', 'worker.py', 'worker', '1orq']`
