## Findings

- `ConDiv_symlay/run_remote.sh` is intentionally hard-coded to a single 48-task Slurm allocation and then dispatches workers via `srun`, so the cluster only sees one job ID even though multiple worker steps run inside it.
- `ConDiv_symlay/ConDiv_mem.py` combines three concerns in `run_minibatch()`: staging the minibatch forcefield inputs, running each target worker, and aggregating divergences into the FF update. The new Slurm workflow needs those concerns split so simulation jobs and the update job can be separate Slurm jobs.
- The existing minibatch logic already tolerates worker failures as long as at least one target succeeds. The new dependency chain should therefore use `afterany` rather than `afterok` for the update job.
- The staged workflow must preserve `CONDIV_SYMLAY_LAYER_MANIFEST` during both initial round staging and update-job resubmission; otherwise the temporary membrane file can skip the symmetry projection step.
- User correction: the simulation side must be one Slurm array `sbatch` file per round, not one `sbatch` file per protein.
- The simplest compliant array design is one generated `simulate_array.sbatch` file per round plus the round manifest; `SLURM_ARRAY_TASK_ID` indexes into the manifest task order and resolves the per-target worker spec.
