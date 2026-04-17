# hybrid-interface-sweep

This folder now runs an interface-scale RMSF calibration for `1rkl`.

It compares two methods:

- Reference: `example/08.MembraneSimulation/3.fixed_curvature.run.py` physics for `1rkl`
- Swept method: `example/16.MARTINI/run_sim_1rkl.sh`

The sweep variable is `PROTEIN_ENV_INTERFACE_SCALE`. The analysis now:

- derives a fixed embedded-residue region from the reference membrane geometry,
- ignores solution residues for scoring,
- fits a trend line over embedded-region RMSD amplitude difference versus `interface_scale`,
- still writes embedded-region RMSF profiles as a secondary diagnostic.

## Scope

- Each `reference` task runs one independent example-08-style fixed-curvature replicate for `1rkl`.
- Each `hybrid` task runs one example-16 hybrid replicate for one `interface_scale`.
- The workflow supports local execution and Slurm array staging/submission.
- The workflow also supports a post-run embedded-region analysis phase over completed outputs.

## Default Sweep

- targeted confirmation sweep:
  - `interface_scale = 0.55, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725, 0.75, 0.775, 0.80, 0.85`
- `reference_replicates = 4`
- `hybrid_replicates = 5`
- `pdb_id = 1rkl`
- RMSF burn-in fraction: `0.20`
- trend-line sample count: `201`
- longer default trajectories for smoother RMSF estimates:
  - reference `REFERENCE_DURATION = 200001`
  - hybrid `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1000`
  - hybrid `PROD_70_NSTEPS = 50000`
  - hybrid `EQ_FRAME_STEPS = 250`
  - hybrid `PROD_FRAME_STEPS = 100`

The default output directory is `hybrid-interface-sweep/runs/default`.

Old softening or diffusion manifests are intentionally incompatible with this workflow. Reinitialize with `HYBRID_SWEEP_FORCE_INIT=1` or use a fresh `HYBRID_SWEEP_BASE_DIR`.

## Run Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_HYBRID_REPLICATES`
- `HYBRID_SWEEP_REFERENCE_REPLICATES`
- `HYBRID_SWEEP_REPLICATES`
  - fallback that applies the same count to both task families
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_BURN_IN_FRACTION`
- `HYBRID_SWEEP_TRENDLINE_SAMPLES`
- `HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN`
- `HYBRID_SWEEP_FORCE_INIT=1`
- `HYBRID_SWEEP_MAX_TASKS`
- `HYBRID_SWEEP_START_TASK`
- `HYBRID_SWEEP_OVERWRITE=1`

## Submit To Slurm

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/submit_remote_round.sh
```

This stages and submits a Slurm array with one array element per task from `sweep_manifest.json`, followed by a dependent collector job.

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_INTERFACE_SCALES`
- `HYBRID_SWEEP_HYBRID_REPLICATES`
- `HYBRID_SWEEP_REFERENCE_REPLICATES`
- `HYBRID_SWEEP_REPLICATES`
- `HYBRID_SWEEP_PDB_ID`
- `HYBRID_SWEEP_SEED`
- `HYBRID_SWEEP_BURN_IN_FRACTION`
- `HYBRID_SWEEP_TRENDLINE_SAMPLES`
- `HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN`
- `HYBRID_SWEEP_FORCE_INIT=1`
- `HYBRID_SWEEP_NO_SUBMIT=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

Optional Slurm resource environment variables:

- `HYBRID_SWEEP_TRAIN_WALLTIME`
- `HYBRID_SWEEP_COLLECT_WALLTIME`
- `HYBRID_SWEEP_CPUS_PER_TASK`
- `HYBRID_SWEEP_SBATCH_PARTITION`
- `HYBRID_SWEEP_SBATCH_ACCOUNT`
- `HYBRID_SWEEP_SBATCH_QOS`
- `HYBRID_SWEEP_SBATCH_CONSTRAINT`
- `HYBRID_SWEEP_SBATCH_MEM`

## Runtime Overrides

### Reference Method

The example-08 reference side is configured from these environment variables at `init-run` time:

- `REFERENCE_TEMPERATURE`
- `REFERENCE_DURATION`
- `REFERENCE_FRAME_INTERVAL`
- `REFERENCE_MEMBRANE_THICKNESS`
- `REFERENCE_USE_CURVATURE`
- `REFERENCE_CURVATURE_RADIUS`
- `REFERENCE_CURVATURE_SIGN`

Defaults match `example/08.MembraneSimulation/3.fixed_curvature.run.py` for `1rkl`.
The workflow lengthens the default reference duration to `200001` steps to reduce RMSF noise while keeping the same reference physics.

### Hybrid Method

The hybrid side replays a whitelist of `example/16.MARTINI/run_sim_1rkl.sh` overrides captured at `init-run` time.

The workflow now seeds longer default hybrid durations into that captured environment unless you override them explicitly:

- `EQ_62_NSTEPS = 1000`
- `EQ_63_NSTEPS = 1000`
- `EQ_64_NSTEPS = 1000`
- `EQ_65_NSTEPS = 1000`
- `EQ_66_NSTEPS = 1000`
- `PROD_70_NSTEPS = 50000`
- `EQ_FRAME_STEPS = 250`
- `PROD_FRAME_STEPS = 100`

Examples:

- `TEMPERATURE`
- `THERMOSTAT_TIMESCALE`
- `MIN_60_MAX_ITER`
- `EQ_62_NSTEPS` through `EQ_66_NSTEPS`
- `PROD_70_NSTEPS`
- `EQ_FRAME_STEPS`
- `PROD_FRAME_STEPS`
- `PROD_70_NPT_ENABLE`
- `PROD_70_BAROSTAT_TYPE`
- `PROTEIN_LIPID_CUTOFF`
- `SALT_MOLAR`
- `SC_ENV_LJ_FORCE_CAP`
- `SC_ENV_COUL_FORCE_CAP`

To capture additional hybrid environment keys, set:

- `HYBRID_SWEEP_EXTRA_ENV_KEYS=KEY1,KEY2,...`

To capture additional reference-side environment keys, set:

- `HYBRID_SWEEP_EXTRA_REFERENCE_KEYS=KEY1,KEY2,...`

## Analysis Outputs

The post-run analysis derives the membrane-embedded residue region from the example-08 reference trajectories and then compares only that selected region.

The analysis uses:

- rigid 3D alignment over the selected embedded-region backbone carrier atoms
- embedded-region residue RMSF from the mapped CA-like backbone atom
- embedded-region framewise RMSD amplitude as the primary calibration metric
- condition-level embedded-region RMSD amplitude difference versus `interface_scale`
- a fitted trend line over `condition_embedded_region_rmsd_delta_vs_reference_angstrom`

The embedded-region selector is reference-driven:

- compute residue membrane occupancy from the reference membrane geometry after burn-in,
- keep residues whose occupancy is at least `HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN` (default `0.50`),
- apply that same residue subset to all reference and hybrid task analyses.

The analysis writes under `BASE_DIR/analysis/`:

- `analysis_manifest.json`
- `results/tasks/*.json`
- `assembled/task_results.csv`
- `assembled/residue_rmsf_profiles.csv`
- `assembled/reference_profile.csv`
- `assembled/condition_profiles.csv`
- `assembled/condition_summary.csv`
- `assembled/trendline_points.csv`
- `assembled/interface_scale_vs_rmsf_difference.png`
- `assembled/interface_scale_vs_rmsf_difference.svg`
- `assembled/best_interface_scale_rmsf_vs_reference.png`
- `assembled/best_interface_scale_rmsf_vs_reference.svg`
- `assembled/failed_tasks.csv`
- `assembled/recommendation_summary.json`
- `assembled/summary.json`

The main plot-ready files are:

- observed condition points: `assembled/condition_summary.csv`
- fitted trend line: `assembled/trendline_points.csv`
- rendered figure: `assembled/interface_scale_vs_rmsf_difference.png`
- vector figure: `assembled/interface_scale_vs_rmsf_difference.svg`
- best-scale embedded-region RMSF overlay: `assembled/best_interface_scale_rmsf_vs_reference.png`
- best-scale embedded-region RMSF overlay vector: `assembled/best_interface_scale_rmsf_vs_reference.svg`

The recommended scale is saved in:

- `assembled/recommendation_summary.json`

That recommendation file also records which sampled scale was used for the RMSF overlay plot:

- exact best sampled scale when that is the chosen target,
- otherwise the nearest stable sampled scale to the trend-line recommendation.

Important:

- old `analysis_manifest.json` / `analysis/results/tasks/*.json` files from the full-protein RMSF version are not compatible with this embedded-region RMSD analysis;
- rerun `init-analysis` with `HYBRID_SWEEP_FORCE_ANALYSIS_INIT=1` and recompute analysis tasks before assembling.

## Run Analysis Locally

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/run_analysis_local.sh
```

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_FORCE_ANALYSIS_INIT=1`
- `HYBRID_SWEEP_ANALYSIS_MAX_TASKS`
- `HYBRID_SWEEP_ANALYSIS_START_TASK`
- `HYBRID_SWEEP_ANALYSIS_NO_ASSEMBLE=1`
- `HYBRID_SWEEP_OVERWRITE=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

## Submit Analysis To Slurm

```bash
source .venv/bin/activate
source source.sh
./hybrid-interface-sweep/submit_analysis.sh
```

This follows the same pattern as the run phase:

- `init-analysis` discovers successful reference and hybrid outputs,
- `submit-analysis-slurm` stages an analysis Slurm array,
- a dependent collector job assembles the RMSF tables and trend-line outputs.

Useful wrapper environment variables:

- `HYBRID_SWEEP_BASE_DIR`
- `HYBRID_SWEEP_FORCE_ANALYSIS_INIT=1`
- `HYBRID_SWEEP_NO_SUBMIT=1`
- `HYBRID_SWEEP_OVERWRITE=1`
- `HYBRID_SWEEP_PYTHON`
- `HYBRID_SWEEP_HDF5_MODULE`

Optional analysis Slurm resource environment variables:

- `HYBRID_SWEEP_ANALYSIS_WALLTIME`
- `HYBRID_SWEEP_ANALYSIS_COLLECT_WALLTIME`
- `HYBRID_SWEEP_ANALYSIS_CPUS_PER_TASK`
- `HYBRID_SWEEP_SBATCH_PARTITION`
- `HYBRID_SWEEP_SBATCH_ACCOUNT`
- `HYBRID_SWEEP_SBATCH_QOS`
- `HYBRID_SWEEP_SBATCH_CONSTRAINT`
- `HYBRID_SWEEP_SBATCH_MEM`
