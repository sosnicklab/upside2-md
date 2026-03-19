# Analysis Workflow

This directory contains the maintained HDX/HX-MS analysis workflow for this repository.

The active workflow scripts are:

- `analysis.sh`
- `0.run_HXMS.py`
- `1.config.py`
- `2.traj_ana.sh`
- `3.get_protaction_states.sh`
- `4.calc_D_uptake.py`
- `5.analyze_D_uptake.py`

The numbered scripts are the workflow steps. `analysis.sh` is the recommended driver for running the full workflow in order.

## Expected Directory Layout

Run the workflow from a simulation directory that contains:

- `inputs/`
- `outputs/`
- `pdb/`
- `results/` will be created if needed

Typical required inputs are:

- `inputs/<pdb_id>.fasta`
- `inputs/<pdb_id>.chain_breaks`
- `inputs/<pdb_id>.initial.npy` if you use native initialization in `1.config.py`
- `pdb/<pdb_id>_rawuptake_HXMS.csv` if you want the HX-MS comparison workflow
- `outputs/<sim_id>/<pdb_id>.run.<replica>.up` for Upside trajectory analysis

The shell-based steps require the Upside environment (`UPSIDE_HOME`) and the Python environment from the project root.

## Recommended Full Workflow

Use `analysis.sh` unless you need to run individual steps by hand.

The driver does the following:

- sources `${PROJECT_ROOT}/source.sh`
- activates `${PROJECT_ROOT}/.venv`
- runs the numbered workflow in order
- can run locally or submit itself through Slurm
- can skip the experiment-dependent branch when no experimental HX-MS data are available

The default driver order is:

1. `0.run_HXMS.py`
2. `1.config.py`
3. `2.traj_ana.sh`
4. `3.get_protaction_states.sh`
5. `4.calc_D_uptake.py` with `analysis_mode=uptake`
6. `4.calc_D_uptake.py` with `analysis_mode=stability`
7. `4.calc_D_uptake.py` with `analysis_mode=pca`
8. `5.analyze_D_uptake.py` with `analysis_mode=uptake`
9. `5.analyze_D_uptake.py` with `analysis_mode=dg_summary`

If `SKIP_EXPERIMENT_DATA="true"`, the driver skips:

- `0.run_HXMS.py`
- `5.analyze_D_uptake.py` with `analysis_mode=uptake`

## Run From This Repository

Example from a checked-in workflow directory:

```bash
cd example/GlpG_test
bash ../../analysis_scripts/analysis.sh
```

Before running, edit the `#CHECKME` settings at the top of `analysis.sh`:

- `RUNNER`
- `WORK_DIR`
- `PDB_ID`
- `SIM_ID`
- `N_REP`
- `START_FRAME`
- `SKIP_EXPERIMENT_DATA`
- `HXMS_METHOD`
- `PROTEIN_STATE`
- `EXP_DATA_FILE`
- Slurm settings if `RUNNER="slurm"`

## Copy Workflow Into A Simulation Directory

If you want the workflow scripts to live directly inside a simulation directory, copy the contents of `analysis_scripts/` into that directory, then point the copied workflow back to the main project root.

Example:

```bash
SIM_DIR=/path/to/my_simulation
cp -r ./analysis_scripts/* "$SIM_DIR"/
cd "$SIM_DIR"
```

After copying:

1. Open `analysis.sh`.
2. Change the `PROJECT_ROOT=` assignment so it points to the repository root that contains `source.sh`, `.venv`, `py/`, and `obj/`.
3. Update the `#CHECKME` settings for the target system.
4. Run the workflow from the simulation directory.

Local run:

```bash
bash ./analysis.sh
```

Slurm run:

```bash
bash ./analysis.sh
```

For Slurm, first set `RUNNER="slurm"` in `analysis.sh`.

Notes for copied workflows:

- Keep `WORK_DIR="./"` if you want the copied scripts to operate on the current simulation directory.
- Keep the `helpers/` directory beside the copied top-level scripts.
- The shell steps still depend on `${PROJECT_ROOT}/source.sh` to define `UPSIDE_HOME`.

## What Each Script Does

### `analysis.sh`

This is the full-workflow driver.

Use it when you want one place to set the main run parameters and execute the active workflow end to end. It is the only script that manages:

- environment activation
- local versus Slurm execution
- skipping the experiment-dependent branch
- fixed default step ordering

### `0.run_HXMS.py`

This preprocesses experimental HX-MS uptake data from:

- `pdb/<pdb_id>_rawuptake_HXMS.csv`

It writes peptide metadata and experiment-derived uptake arrays under:

- `outputs/<sim_id>/`

Use it when you have raw experimental HX-MS data and want to compare experiment against simulation in step 5.

Supported subworkflows:

- `HXMS_SUBWORKFLOW=all`
- `HXMS_SUBWORKFLOW=normalized`
- `HXMS_SUBWORKFLOW=stretched`

Main outputs include:

- `<pdb_id>_pep_ids.csv`
- normalized `%D` arrays
- stretched-exponential fit arrays
- HX-MS preprocessing PDFs

### `1.config.py`

This builds the Upside HDX configuration file for the target system.

Default usage:

- `workflow_action=config`

It reads the model inputs in `inputs/` and writes:

- `inputs/<pdb_id>-HDX.up`

This step is Upside-specific. It requires `UPSIDE_HOME` and the Upside parameter files from the main project.

Special mode:

- `workflow_action=select_h5_frames`

In this mode, the script filters existing trajectory `.up` files by RMSD and writes `selected.*.up` files.

### `2.traj_ana.sh`

This extracts per-replica observables from the Upside trajectories in:

- `outputs/<sim_id>/`

For each replica it writes the arrays used by downstream steps, including:

- `results/<pdb_id>_<sim_id>_<replica>_Energy.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Rg.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Rmsd.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Hbond.npy`
- `results/<pdb_id>_<sim_id>_<replica>_T.npy`

It also writes a `.vtf` trajectory for visualization.

This step calls:

- `helpers/get_info_from_upside_traj.py`
- `${UPSIDE_HOME}/py/extract_vtf.py`

### `3.get_protaction_states.sh`

This extracts protection-state arrays from the Upside trajectories.

For each replica it writes:

- `results/<pdb_id>_<sim_id>_<replica>_PS.npy`

It also writes:

- `results/<pdb_id>.resid`

This step calls `${UPSIDE_HOME}/py/get_protection_state.py`, so it is also Upside-specific.

### `4.calc_D_uptake.py`

This is the main simulation-side analysis step. It reads the arrays produced by steps 2 and 3 and runs one of three branches:

- `analysis_mode=uptake`
- `analysis_mode=stability`
- `analysis_mode=pca`

`analysis_mode=uptake`:

- computes simulated deuterium uptake
- writes whole-protein and peptide-level uptake outputs
- can run without experimental data, although experiment-derived peptide/time definitions are still needed for the full peptide-comparison branch

`analysis_mode=stability`:

- computes residue-level stability and protection free-energy style summaries
- writes convergence plots
- writes residue `dG` and related outputs
- uses optional overlays from `pdb/<pdb_id>_HXMS.csv`, `pdb/<pdb_id>_NMR.csv`, and `pdb/<pdb_id>_NMR_MS.csv` when present

`analysis_mode=pca`:

- runs trajectory PCA and related RMSD/H-bond diagnostics

### `5.analyze_D_uptake.py`

This is the comparison and summary step. It supports:

- `analysis_mode=uptake`
- `analysis_mode=dg_summary`
- `analysis_mode=compare_hxms`

`analysis_mode=uptake`:

- compares experimental HX-MS peptide uptake against simulated uptake
- writes comparison plots and cooperativity-style summaries

`analysis_mode=dg_summary`:

- summarizes residue-level `dG` outputs across one or more proteins

`analysis_mode=compare_hxms`:

- compares HX-MS peptide uptake across proteins using `pdb_ids` and `protein_states`

## Manual Run Order

If you do not want to use `analysis.sh`, activate the environment first:

```bash
source "${PROJECT_ROOT}/source.sh"
source "${PROJECT_ROOT}/.venv/bin/activate"
```

Then run the steps from the workflow directory.

With experimental HX-MS data:

```bash
pdb_id=MY_PROTEIN sim_id=MY_SIM python 0.run_HXMS.py
workflow_action=config pdb_id=MY_PROTEIN python 1.config.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 bash 2.traj_ana.sh
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 bash 3.get_protaction_states.sh
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=uptake python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=stability python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=pca python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=uptake python 5.analyze_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=dg_summary python 5.analyze_D_uptake.py
```

Without experimental HX-MS data:

```bash
workflow_action=config pdb_id=MY_PROTEIN python 1.config.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 bash 2.traj_ana.sh
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 bash 3.get_protaction_states.sh
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=uptake python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=stability python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=pca python 4.calc_D_uptake.py
pdb_id=MY_PROTEIN sim_id=MY_SIM n_rep=48 analysis_mode=dg_summary python 5.analyze_D_uptake.py
```

## Using Results From Other MD Software

The full numbered workflow is Upside-specific through steps 1 to 3.

You can still reuse the downstream analysis for trajectories from other MD software, but you must first produce the same files that steps 2 and 3 normally write into `results/`.

Minimum files needed by step 4:

- `results/<pdb_id>_<sim_id>_<replica>_Energy.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Rg.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Rmsd.npy`
- `results/<pdb_id>_<sim_id>_<replica>_Hbond.npy`
- `results/<pdb_id>_<sim_id>_<replica>_T.npy`
- `results/<pdb_id>_<sim_id>_<replica>_PS.npy`
- `results/<pdb_id>.resid`
- `inputs/<pdb_id>.fasta`

Recommended adaptation points:

- Adapt `helpers/get_info_from_upside_traj.py` if you want a step-2-like script that reads your trajectory format and writes the observable arrays expected by step 4.
- Replace the protection-state generation from step 3 if `${UPSIDE_HOME}/py/get_protection_state.py` does not apply to your trajectories. The important requirement is to generate compatible `*_PS.npy` arrays and `results/<pdb_id>.resid`.

Once those files exist, you can start from:

- `4.calc_D_uptake.py`
- `5.analyze_D_uptake.py`

If you also want experiment-versus-simulation comparison, you still need the experiment-derived HX-MS outputs from step 0 or equivalent files with the same content.

## Main Output Families

The exact outputs depend on the chosen modes, but the main result families are:

- experiment-derived peptide arrays in `outputs/<sim_id>/`
- per-replica observable arrays in `results/`
- per-replica protection states in `results/`
- whole-protein and peptide uptake outputs in `results/`
- residue `dG` and stability plots in `results/`
- PCA plots in `results/`
- comparison and summary plots in `results/`
