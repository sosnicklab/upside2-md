## Lessons

- When a user reports a scheduling/orchestration failure, inspect how work maps onto cluster job units before assuming CPU counts imply parallel execution.
- When a user specifies a Slurm array workflow, do not translate that into many separate `sbatch` files; keep one array script and map one simulation to one array task.
- When introducing new cluster workflow entrypoints, add a preflight check for required files in the remote checkout so missing-sync failures are explicit.
- When a user says an error artifact comes from the current code, do not trust prior notes or earlier local reruns over that statement; re-reproduce the active code path directly before claiming the bug is fixed.
