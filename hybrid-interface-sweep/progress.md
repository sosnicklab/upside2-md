# Progress Log

## 2026-04-13 (Workflow Redesign)
- Reopened `hybrid-interface-sweep/` after the user corrected the scientific goal:
  - the existing workflow answered a hybrid protein/bilayer interface-scale question,
  - but the requested workflow needed to be bilayer-only and softening-focused.
- Rewrote the local tracker files to make bilayer-only production softening, bilayer diffusion, and a viscosity proxy the source of truth for the task.

## 2026-04-13 (Implementation)
- Replaced `hybrid-interface-sweep/workflow.py` so the sweep now:
  - prepares bilayer-only systems through `py/martini_prepare_system.py --mode bilayer`,
  - runs a bilayer `6.0 -> 7.0` stage ladder,
  - sweeps `(lj_alpha, slater_alpha, replicate)` instead of `PROTEIN_ENV_INTERFACE_SCALE`,
  - rewrites the staged production `martini_potential` attrs for task-specific softening,
  - assembles bilayer-only task summaries,
  - discovers successful `stage_7.0.up` files for analysis,
  - computes `PO4` lateral diffusion relative to bilayer COM,
  - converts diffusion into physical units using the configured `ps/step` assumption,
  - reports a reciprocal-diffusion viscosity proxy,
  - writes `recommendation_summary.json` from the assembled analysis.
- Updated wrapper scripts:
  - `hybrid-interface-sweep/run_local.sh`
  - `hybrid-interface-sweep/submit_remote_round.sh`
- Rewrote `hybrid-interface-sweep/README.md` around the corrected bilayer-only parameter surface and outputs.

## 2026-04-13 (Verification)
- Ran static verification:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
- Ran a reduced one-task bilayer smoke sweep under `/tmp/hybrid_interface_softening_smoke` with:
  - `lj_alpha = 0.05`
  - `slater_alpha = 0.5`
  - `replicates = 1`
  - `MIN_60_MAX_ITER = 1`
  - `MIN_61_MAX_ITER = 1`
  - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
  - `PROD_70_NSTEPS = 40`
  - `EQ_FRAME_STEPS = 1`
  - `PROD_FRAME_STEPS = 1`
- Verified the sweep smoke outputs:
  - `results/tasks/lj0p05_coul0p5_r01.json` reports `success = true`
  - `assembled/summary.json` reports one completed task and one completed condition
  - staged `bilayer.stage_7.0.up` carries:
    - `lj_soften = 1`
    - `lj_soften_alpha = 0.05`
    - `coulomb_soften = 1`
    - `slater_alpha = 0.5`
- Ran reduced local analysis on the same smoke tree and verified:
  - `analysis/results/tasks/lj0p05_coul0p5_r01.json`
  - `analysis/assembled/task_results.csv`
  - `analysis/assembled/condition_summary.csv`
  - `analysis/assembled/recommendation_summary.json`
  - `analysis/assembled/summary.json`
- Verified no-submit Slurm staging:
  - `slurm/round_manifest.json`
  - `slurm/run_array.sbatch`
  - `slurm/collect_results.sbatch`
  - `analysis/slurm/round_manifest.json`
  - `analysis/slurm/analyze_array.sbatch`
  - `analysis/slurm/collect_analysis.sbatch`
