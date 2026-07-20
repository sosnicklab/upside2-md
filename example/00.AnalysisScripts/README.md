# HX-MS Analysis Workflow

The numbered scripts are run in order by `analysis.sh`. Edit the `# checkme` settings at the top of `analysis.sh` before running.

## Quick Start

```bash
source ../../source.sh
bash analysis.sh
```

The workflow accepts ordinary Upside trajectories by default. Full dry-MARTINI
hybrid trajectories are projected onto the ordinary protein HDX trajectory
contract before the numbered analysis steps run; the HDX estimators themselves
are not duplicated.

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

## Full dry-MARTINI hybrid trajectories

Set `trajectory_mode=martini_hybrid`, the directory containing the full hybrid
H5 files, and a filename pattern containing `{replica}`. Step 2 calls
`py/martini_hdx_project.py` to create ordinary protein-only trajectory views in
`outputs/<sim_id>/`. The views contain mapped N/CA/C coordinates but retain the
full protein-plus-bilayer potential and the original Upside temperature for
MBAR. Steps 2--6 then run through the same code used by `example/04.HDX` and the
ordinary workflow.

```bash
work_dir=/path/to/hdx_work \
pdb_id=1rkl sim_id=martini_hybrid n_rep=1 start_frame=100 \
trajectory_mode=martini_hybrid \
hybrid_source_dir=/path/to/martini_1rkl_hybrid_full/checkpoints \
hybrid_source_pattern='1rkl.stage_7.{replica}.up' \
legacy_T_range=0.8647 T_targets=0.8647 target_T=0.8647 \
skip_experiment_data=true \
bash analysis.sh
```

For a single-temperature trajectory, set every analysis target to that sampled
Upside temperature. Do not temperature-reweight one hybrid replica to an
unsampled condition. `example/16.MARTINI/run_hdx_analysis.sh` reads the
temperature from replica 0 and supplies these settings automatically.

The work directory must contain the normal `inputs/<pdb>.fasta`,
`inputs/<pdb>.initial.npy`, `inputs/<pdb>.chain_breaks`, and `pdb/` directory.
`1.config.py` builds the same heavy-atom-coverage `inputs/<pdb>-HDX.up` used for
ordinary Upside HDX analysis.

Step 3 always saves the unmodified protein result as
`<prefix>_PS_protein.npy`. With no additional input, downstream
`<prefix>_PS.npy` is identical. A calibrated per-frame water-accessibility array
may be supplied with `water_accessibility_dir` and a
`water_accessibility_pattern` containing `{replica}`. The array must match PS
shape and contain values in `[0,1]`; the combined exchange probability is

```text
(1 - protein_protection) * water_accessibility
```

Dry MARTINI has no water, so the workflow does not invent this calibration. A
global between-phosphate-planes Boolean is not a quantitative substitute: it
saturates stable transmembrane amides. Leave the accessibility pattern empty
unless the array has been calibrated against explicit-water simulation or an
experimental reference.

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
| `<pdb>_<sim>_<rep>_T.npy` | scalar | simulation `kT` in Upside units, for example `0.8647` |
| `<pdb>_<sim>_<rep>_PS.npy` | `(n_frames, n_nh)` | per-frame protection probability in `[0,1]` |
| `<pdb>.resid` | text, one int per line | NH donor residue indices (non-prolines) |
| `<pdb>.fasta` | FASTA | protein sequence (in `inputs/`) |

Rmsd and Rg should be in Angstroms. `T.npy` is not Kelvin: step 4 uses it
directly in `beta=1/T`, then performs a separate Kelvin conversion for intrinsic
exchange chemistry. PS columns must match the order of residues in `.resid`.
Hbond.npy can be zeros because step 4 does not use it. Energy.npy can be constant
only when no per-frame energies are available; for a hybrid trajectory it must
be the full coupled-system potential, not the protein-only analysis energy.

### Timescale interpretation

The uptake model is EX2-like and uses equilibrium protection probabilities:
`k_obs = k_chem * (1 - <PS>)`. Simulated elapsed time is not used as the HX
labeling time. A mismatch between lipid lateral diffusion and the empirical
protein clock therefore does not directly rescale the calculated uptake curve.
It does affect decorrelation and convergence. Quantitative results require
replica overlap or independent-seed agreement, block convergence of every
residue's protection probability, timestep convergence, and experimental
peptide-level validation.
