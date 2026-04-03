# SC-training

This folder contains the dry-MARTINI sidechain-type to dry-MARTINI particle-type training workflow used to build the stage-7 `martini.h5` sidechain library.

## Scope

- The workflow is independent from the legacy SC back-mapping / rotamer-placement production path.
- It derives per-residue/per-target training curves from dry-MARTINI forcefield inputs already present in this repository.
- It can run locally or stage/submit Slurm array jobs.
- By default it builds the complete table over all `38` dry-MARTINI particle types from `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`, giving `18 x 38 = 684` residue-target tasks for the current non-empty canonical sidechains.
- The canonical benchmark workflow target is `example/16.MARTINI/run_sim_1rkl.sh`.
- The default geometry/data sources are:
  - dry-MARTINI sidechain bead types from `example/16.MARTINI/martinize.py` (`martini22`);
  - dry-MARTINI nonbond parameters from `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`;
  - original Upside residue-local sidechain geometry from `parameters/ff_2.1/sidechain.h5`.

## Current Contract

The current workflow follows the user-specified one-sided Upside orientation form:

- each residue sidechain is represented as the set of dry-MARTINI bead types defined by `example/16.MARTINI/martinize.py` for `martini22`;
- target dry-MARTINI particles are sampled in the original residue-local `CB` frame recovered from the Upside sidechain library;
- the angular coordinate is the one-sided Upside sidechain vector term `-n_1 . n_12`, where:
  - `n_1` is the backbone-defined canonical `CA -> CB` unit vector shared with Upside runtime;
  - `n_12` is the sampled `CB ->` dry-MARTINI target direction;
- the assembled forcefield is stored in factorized form:
  - `V_hybrid(r,theta) = V_radial(r) + ang_1(-n_1 . n_12) * V_angular(r)`;
- task JSON outputs retain the sampled directional surface and the fitted factorized form so fit quality can be audited.

This contract is written into the manifest and assembled outputs so training and runtime use the same orientation model.

## Unit Contract

The workflow emits native dry-MARTINI units only:

- distance: `nm`
- energy: `kJ/mol`
- charge: `e`
- Coulomb constant: `138.935458 / 15 = 9.262363866666668 kJ/mol/nm/e^2`

No runtime conversion values are baked into the training outputs. Simulation-side code must supply its own unit-conversion parameters when these tables are consumed.

## Local Usage

From the project root:

```bash
source .venv/bin/activate
source source.sh
./SC-training/run_local.sh
```

Optional environment variables:

- `SC_TRAIN_BASE_DIR`: run directory, default `SC-training/runs/default`
- `RUN_BENCHMARK=1`: run the canonical `run_sim_1rkl.sh` benchmark after local training assembly
- `SC_TRAIN_PYTHON`: Python interpreter to use if you do not want the default `python3`

## Slurm Usage

From the project root on a Slurm system:

```bash
source .venv/bin/activate
source source.sh
./SC-training/submit_remote_round.sh
```

Optional environment variables:

- `SC_TRAIN_BASE_DIR`: run directory, default `SC-training/runs/default`
- `SC_TRAIN_RUN_BENCHMARK=1`: submit a dependent benchmark job after result assembly
- `SC_TRAIN_NO_SUBMIT=1`: stage Slurm scripts without calling `sbatch`
- `SC_TRAIN_PYTHON`: Python interpreter to use if you do not want the default `python3`
- `SC_TRAIN_SBATCH_PARTITION`
- `SC_TRAIN_SBATCH_ACCOUNT`
- `SC_TRAIN_SBATCH_QOS`
- `SC_TRAIN_SBATCH_CONSTRAINT`
- `SC_TRAIN_SBATCH_MEM`
- `SC_TRAIN_TRAIN_WALLTIME`
- `SC_TRAIN_COLLECT_WALLTIME`
- `SC_TRAIN_BENCH_WALLTIME`

To stage scripts without submitting:

```bash
python3 SC-training/workflow.py submit-slurm --base-dir SC-training/runs/default --no-submit
```

## Main Commands

```bash
python3 SC-training/workflow.py init-run --base-dir SC-training/runs/default
python3 SC-training/workflow.py run-local --base-dir SC-training/runs/default
python3 SC-training/workflow.py assemble-results --base-dir SC-training/runs/default
python3 SC-training/workflow.py submit-slurm --base-dir SC-training/runs/default
python3 SC-training/workflow.py run-benchmark --base-dir SC-training/runs/default --execute
```

## Outputs

Each run directory contains:

- `training_manifest.json`: residue/target task manifest and training assumptions
- `results/tasks/*.json`: one result per residue-target task
- `results/assembled/sc_table.json`: assembled table
- `results/assembled/sc_table_summary.json`: simple summary
- `slurm/`: generated Slurm manifests and scripts
- `benchmark/`: staged benchmark commands and optional benchmark run logs

## Benchmark Note

The workflow already targets `example/16.MARTINI/run_sim_1rkl.sh` as the benchmark entrypoint because that is the real production workflow. Runtime consumption of the assembled SC table is still a separate C++ integration step.

If you upload only the `SC-training/` folder, training and Slurm staging will still work. The benchmark path remains optional and still requires access to the full Upside/MARTINI workflow repository.

The benchmark / simulation side now expects explicit unit-conversion parameters in the environment:

- `UPSIDE_MARTINI_ENERGY_CONVERSION`
- `UPSIDE_MARTINI_LENGTH_CONVERSION`
