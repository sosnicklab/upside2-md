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

## Build New Dataset From RCSB (X-ray + NMR)

Use `prepare_rcsb_upside_input.py` to query/filter RCSB, download structures, and generate Upside-format inputs.

Example:

```bash
source venv/bin/activate
python3 prepare_rcsb_upside_input.py \
  --output-root rcsb_dataset \
  --filter-profile contrastive \
  --max-candidates-per-class 1000 \
  --max-keep-per-class 250
```

Generated outputs:

- `rcsb_dataset/xray/upside_input/` (accepted X-ray structures as `.fasta/.initial.pkl/.chi`)
- `rcsb_dataset/nmr/upside_input/` (accepted NMR structures as `.fasta/.initial.pkl/.chi`)
- `rcsb_dataset/lists/pdb_list_xray`
- `rcsb_dataset/lists/pdb_list_nmr`
- `rcsb_dataset/lists/pdb_list_combined`
- `rcsb_dataset/lists/structure_index_xray.csv` (includes `deposit_year`)
- `rcsb_dataset/lists/structure_index_nmr.csv` (includes `deposit_year`)
- `rcsb_dataset/lists/structure_index_combined.csv` (includes `deposit_year`)
- `rcsb_dataset/manifest.csv` (keep/skip reasons)

Default filters are aligned to current project data:

- `contrastive` profile: sequence length 50 to 100 residues (`sidechain`: 50 to 500)
- X-ray resolution <= 2.2 A
- pairwise sequence similarity <= 30% (SequenceMatcher-based approximation)
- globularity outlier removal on log(Nres)-log(Rg) with RANSAC-style fit
- excludes non-protein / nucleic-acid-containing entries
- excludes AlphaFold/ModelArchive-style identifiers (`AF_*`, `MA_*`)
- excludes non-single-chain proteins
- excludes heme-containing structures (e.g., `HEM`, `HEA`, `HEC`, `HEO`)
- excludes structures containing non-standard amino-acid residues
- excludes membrane-related proteins (keyword/title based metadata filter)
- applies post-conversion chain-break filter (`max_sep < 2.0`)

Conversion details:

- Uses the same initial conversion tool as `example/01.GettingStarted/0.run.sh`:
  - `py/PDB_to_initial_structure.py ... --record-chain-breaks`
- Then writes training-compatible `*.initial.pkl` from generated `*.initial.npy`.
