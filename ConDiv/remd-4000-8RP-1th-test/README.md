# ConDiv Training Workflow

This directory contains scripts for initializing and running ConDiv training locally or on a Slurm cluster.

## Prerequisites

- Working directory: this folder (`remd-4000-8RP-1th-test`)
- Python virtual environment at `venv/`
- Input assets present:
  - `init_param/`
  - `upside_input/`
  - `pdb_list`

## 1. Initialize (login node)

Run once before first training submission:

```bash
./run_init.sh
```

This generates `test_00/initial_checkpoint.pkl` and setup files.

## 2. Run One Training Round On Slurm

Load required modules first:

```bash
module load gcc/12.2.0
module load cmake
module load openmpi
```

Submit:

```bash
sbatch run_remote.sh
```

Default restart iterations per submission: `152`.

### Override iteration count

Use positional argument:

```bash
sbatch run_remote.sh 200
```

Or environment variable:

```bash
sbatch --export=RUN_STEPS=200 run_remote.sh
```

## 3. Continue Next Round

After the job finishes, submit again:

```bash
sbatch run_remote.sh
```

`run_remote.sh` automatically resumes from the latest valid checkpoint in:

- `test_00/epoch_*_minibatch_*/checkpoint.pkl`

If none is found, it falls back to:

- `test_00/initial_checkpoint.pkl`

## Logs

- Training log (appended): `run.output`
- Slurm stdout/stderr: `slurm-condiv-train-<jobid>.out`

## Local Run (optional)

For local execution instead of Slurm:

```bash
./run_local.sh
```
