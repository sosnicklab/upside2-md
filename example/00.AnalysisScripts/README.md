# HX-MS Analysis Workflow

The numbered scripts are run in order by `analysis.sh`. Edit the `# checkme` settings at the top of `analysis.sh` before running.

## Quick Start

```bash
source ../../source.sh
bash analysis.sh
```

## Workflow Order

1. `0.run_HXMS.py` — preprocess experimental HX-MS uptake data (skipped if `skip_experiment_data=true`)
2. `1.config.py` — build Upside HDX configuration file
3. `2.traj_ana.sh` — extract per-replica observables from Upside trajectories
4. `3.get_protaction_states.sh` — extract protection states from Upside trajectories
5. `4.calc_D_uptake.py analysis_mode=uptake` — compute simulated deuterium uptake
6. `4.calc_D_uptake.py analysis_mode=stability` — residue-level stability summaries
7. `5.analyze_D_uptake.py analysis_mode=uptake` — compare experiment to simulation (skipped if `skip_experiment_data=true`)
8. `5.analyze_D_uptake.py analysis_mode=dg_summary` — summarize residue dG across proteins
9. `6.generate_hx_plots.py` — legacy HX overview plots (optional, auto-detect)

## Manual Run

```bash
source ../../source.sh
source ../../.venv/bin/activate

pdb_id=GlpG sim_id=Mt_500 n_rep=1 protein_state=pd6 \
    HXMS_method=stretch_exp exp_data_file=GlpG_rawuptake_HXMS.csv \
    python 0.run_HXMS.py

workflow_action=config pdb_id=GlpG python 1.config.py

pdb_id=GlpG sim_id=Mt_500 n_rep=1 bash 2.traj_ana.sh
pdb_id=GlpG sim_id=Mt_500 n_rep=1 bash 3.get_protaction_states.sh

pdb_id=GlpG sim_id=Mt_500 n_rep=1 analysis_mode=uptake \
    protein_state=pd6 exp_data_file=GlpG_rawuptake_HXMS.csv \
    python 4.calc_D_uptake.py

pdb_id=GlpG sim_id=Mt_500 n_rep=1 analysis_mode=stability \
    python 4.calc_D_uptake.py

pdb_id=GlpG sim_id=Mt_500 n_rep=1 analysis_mode=uptake \
    protein_state=pd6 HXMS_method=stretch_exp \
    python 5.analyze_D_uptake.py

pdb_id=GlpG sim_id=Mt_500 n_rep=1 analysis_mode=dg_summary \
    python 5.analyze_D_uptake.py

hx_plot_prefix=glpg hx_plot_state=pd6 python 6.generate_hx_plots.py
```

## Directory Layout

```
inputs/    — fasta, chain_breaks, initial.npy, HDX .up config
outputs/   — simulation trajectories (.up files per replica)
pdb/       — raw uptake CSV, optional NMR/MS overlays
results/   — created automatically, all output goes here
helpers/   — utility modules used by the workflow steps
```

## Slurm

Set `runner=slurm` in `analysis.sh` and adjust `slurm_*` settings. The driver submits itself via `sbatch`.

## Outputs

- per-replica arrays in `results/`: Energy, Rmsd, Rg, Hbond, T, PS
- whole-protein and peptide uptake in `results/`: `*_percentD*.csv`, `*_percentD_feats.csv`
- uptake plots: `D_uptake*.pdf`, `D_uptake*.png`
- dG and stability plots, COF summaries
- experiment-simulation comparison plots (if exp data available)

## REMD and Metadynamics analysis

Two standalone helpers analyze enhanced-sampling runs (independent of the HX-MS pipeline):

**Temperature-REMD** — `helpers/analyze_remd.py`. Upside REMD writes one `.up` per slot at a fixed
temperature (`output/temperature`, in Upside units so kT equals that value in E_up). The tool reads all
replicas, tabulates per-temperature observables, and reconstructs `F(CV)` at a target temperature by
MBAR-reweighting the pooled samples (falls back to a per-temperature histogram if `pymbar` is absent).

```bash
python helpers/analyze_remd.py 'outputs/sim/prot.run.*.up' results/prot --cv rg --target-kT 0.8647
python helpers/analyze_remd.py 'outputs/sim/prot.run.*.up' results/prot --cv dist --cv-atoms 12 87
```

Outputs `*_remd_obs.csv` (T, ⟨E⟩, ⟨CV⟩), `*_remd_fes.npy` (cv, F), `*_remd.png`.

**Metadynamics** — `helpers/analyze_metadynamics.py`. Reconstructs the free-energy surface from the
metadynamics node's logged bias. The last frame of `output/metadynamics_bias` is the converged bias
`V(s)`; the FES is `-(γ/(γ-1))·V(s)` for well-tempered runs (γ = `bias_factor`), else `-V(s)`.

```bash
python helpers/analyze_metadynamics.py outputs/sim/prot.run.0.up results/prot --target-kT 0.8647
```

Outputs `*_fes.npy` (cv, F in E_up), `*_cv.npy` (CV time series), `*_metad.png`.

## Using Other MD Software

Steps 1–3 are Upside-specific. Steps 4–6 only need `.npy` arrays in `results/`. To use trajectories from another MD engine, produce these files with matching names and shapes:

| File | Shape | Content |
|---|---|---|
| `<pdb>_<sim>_<rep>_Energy.npy` | `(n_frames,)` or `(n_frames, 1)` | per-frame potential energy |
| `<pdb>_<sim>_<rep>_Rmsd.npy` | `(n_frames,)` | per-frame CA RMSD |
| `<pdb>_<sim>_<rep>_Rg.npy` | `(n_frames,)` | per-frame radius of gyration |
| `<pdb>_<sim>_<rep>_Hbond.npy` | `(n_frames,)` | per-frame total hbond count |
| `<pdb>_<sim>_<rep>_T.npy` | scalar | simulation temperature in Kelvin |
| `<pdb>_<sim>_<rep>_PS.npy` | `(n_frames, n_nh)` | per-frame binary protection state |
| `<pdb>.resid` | text, one int per line | NH donor residue indices (non-prolines) |
| `<pdb>.fasta` | FASTA | protein sequence (in `inputs/`) |

Rmsd, Rg should be in Angstroms. PS.npy columns must match the order of residues in `.resid`. Hbond.npy can be zeros — step 4 does not use it. Energy.npy can be constant (zero) if no per-frame energies are available; MBAR weights become uniform and temperature scanning still works through `k_chem(T)`.
