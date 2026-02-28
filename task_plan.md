# MARTINI Relaxation Slurm Adaptation Plan

## Project Goal
Make `example/16.MARTINI/run_relax_6x_rigid_dry.sh` usable in a Slurm workflow by submitting all selected MARTINI relaxation runs as Slurm jobs so multiple systems can run in parallel.

## Architecture & Key Decisions
1. Keep `example/16.MARTINI/run_relax_6x_rigid_dry.sh` as the single-system stage runner. It already encapsulates the required sequential stage handoff for one system.
2. Add Slurm support in `example/16.MARTINI/batch_relax_rigid_dry.py`, because that is the layer that already iterates over multiple PDB systems.
3. Submit one Slurm job per PDB system. Stages `6.0` to `6.6` remain sequential inside each job; parallelism is across systems.
4. Preserve the existing local execution path so current non-Slurm usage does not regress.
5. Generate per-system Slurm wrapper scripts in each run directory instead of relying on large `sbatch --wrap` strings. This keeps quoting predictable and leaves an auditable job payload on disk.
6. Support optional waiting/polling so the batch driver can either submit and return immediately or wait for all submitted jobs before writing the relaxed manifest.

## Execution Phases
- [x] Inspect current rigid dry MARTINI workflow and identify the orchestration layer to change.
- [x] Add Slurm runner CLI options and submission helpers to `batch_relax_rigid_dry.py`.
- [x] Preserve existing direct local execution behavior.
- [x] Add submission metadata and completion handling for Slurm jobs.
- [x] Validate Python syntax and CLI help.
- [x] Validate shell syntax for the single-system runner remains clean.
- [x] Record verification results and operating assumptions.

## Known Errors / Blockers
1. Slurm commands may be unavailable inside the current sandbox, so full end-to-end job submission may not be executable here.
2. Cluster-specific defaults such as account, partition, and memory cannot be hardcoded safely; the new interface must leave them configurable.

## Review Criteria
1. Existing local `batch_relax_rigid_dry.py` behavior still works without Slurm flags.
2. Slurm mode submits one job per selected system and allows those jobs to execute in parallel.
3. Each submitted job has a reproducible wrapper script and log path under its run directory.
4. The relaxed manifest is produced correctly after local runs, and after Slurm runs when waiting is enabled.

## Review (This Run)
1. Updated `example/16.MARTINI/batch_relax_rigid_dry.py` to support `--runner slurm` while keeping `--runner local` as the default path.
2. Added Slurm submission controls for account, partition, walltime, CPUs, memory, job-name prefix, extra `sbatch` args, wait/poll behavior, and dry-run wrapper generation.
3. Slurm mode now writes a per-system wrapper script (`submit_relax.slurm.sh`) and submission metadata (`batch_relax.log`) under each run directory, then submits one job per selected PDB system.
4. When `--slurm-wait 1` is used, the batch driver polls `squeue` until all submitted jobs leave the queue and then validates that each expected `stage_6.6.up` output exists before adding it to the relaxed manifest.
5. Validation completed in sandbox:
   - `python3 -m py_compile example/16.MARTINI/batch_relax_rigid_dry.py` passed.
   - `python3 example/16.MARTINI/batch_relax_rigid_dry.py --help` passed.
   - `bash -n example/16.MARTINI/run_relax_6x_rigid_dry.sh` passed.
   - Slurm dry-run on real manifest input passed:
     `python3 example/16.MARTINI/batch_relax_rigid_dry.py --manifest train-data/training_manifest.json --runner slurm --dry-run 1 --slurm-wait 0 --limit 1 --output-root /tmp/upside2-md-slurm-dryrun`.
   - Dry-run generated `/private/tmp/upside2-md-slurm-dryrun/1ors/submit_relax.slurm.sh` with only the workflow-specific exported variables and the expected `sbatch` command in `/private/tmp/upside2-md-slurm-dryrun/1ors/batch_relax.log`.
   - Follow-up wrapper correction applied: generated Slurm scripts now initialize the module system when needed and run `module load cmake` plus `module load openapi` before the workflow payload.
6. Remaining verification gap:
   - Real `sbatch`/`squeue` execution was not attempted in this sandbox, so live Slurm submission and wait-path completion must be exercised on the target cluster.
