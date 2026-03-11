## Lessons

- When a user reports a scheduling/orchestration failure, inspect how work maps onto cluster job units before assuming CPU counts imply parallel execution.
- When a user specifies a Slurm array workflow, do not translate that into many separate `sbatch` files; keep one array script and map one simulation to one array task.
