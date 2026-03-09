Project Goal
- Fix the `ConDiv_symlay` training failure where worker-generated `*.run.N.h5` configuration files cannot be opened reliably by `upside` during epoch execution.

Architecture & Key Decisions
- Treat this as a workflow-runtime bug in config generation or HDF5 file handling, not as a cluster-environment issue until the code path is ruled out.
- Trace the worker call graph from minibatch launch through config generation, HDF5 writes, and `upside` launch.
- Prefer a minimal fix at the file-generation or file-lifetime boundary instead of broad workflow changes.
- Verify by reproducing the failing path locally as far as possible and by checking the resulting HDF5 files explicitly.

Execution Phases
- [x] Locate the worker/config-generation path that produces `*.run.N.h5`.
- [x] Inspect the failing HDF5 files and logs to determine whether the files are missing, truncated, or left in an invalid state.
- [x] Identify the root cause in code and implement the minimal fix.
- [x] Verify the fix by reproducing the relevant generation/open path or by running the smallest viable training smoke test.
- [x] Record review notes and verification status.

Known Errors / Blockers
- User reports `upside` errors opening generated `*.run.N.h5` files inside `epoch_00_minibatch_00`.

Review
- Root cause identified in the Slurm launch model: `ConDiv_symlay` launched a Python worker through `srun`, then launched `upside` through a second nested `srun` with `--ntasks=n_replica`. `obj/upside` is not MPI-linked, so that path spawned multiple independent `upside` processes against the same HDF5 config bundle.
- Implemented the minimal runtime fix in `ConDiv_symlay/ConDiv_mem.py`: keep `srun` only at the Python-worker boundary and launch `upside` locally inside the worker step.
- Reduced transient config storage for each target by eliminating the redundant rewritten `base.h5` copy. The template file is now moved into `run.0.h5` after cloning the remaining replicas.
- Local verification: syntax check passed for `ConDiv_symlay/ConDiv_mem.py`; a local config-generation smoke test produced valid `run.0.h5` and `run.1.h5`, validated against the membrane model, and removed the temporary `base.h5` as intended.
