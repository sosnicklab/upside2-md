# Hybrid Interface Softening Sweep

## Project Goal
- Replace the current mis-scoped `hybrid-interface-sweep` with a bilayer-only dry-MARTINI calibration workflow.
- Sweep production softening settings for Lennard-Jones and Coulomb interactions, measure bilayer fluidity, and use that calibration to choose softening intended for later protein/bilayer interface use.
- Keep the implementation self-contained under `hybrid-interface-sweep/`.

## Architecture & Key Decisions
- Do not use the full hybrid `example/16.MARTINI/run_sim_1rkl.sh` path for this workflow:
  - that path measures the wrong Hamiltonian,
  - and its protein-relative observables do not answer the calibration question.
- Use the existing bilayer preparation surface in `py/martini_prepare_system.py --mode bilayer`.
- Keep the membrane stage ladder close to workflow `16` so the calibration still uses the same dry-MARTINI preparation and equilibration style:
  - `6.0` / `6.1` minimization,
  - `6.2` / `6.3` softened equilibration,
  - `6.4` / `6.5` / `6.6` hard equilibration,
  - `7.0` production.
- Apply the sweep only to the production-stage bilayer Hamiltonian by rewriting `input/potential/martini_potential` attrs in the staged `7.0` file:
  - `lj_soften`,
  - `lj_soften_alpha`,
  - `coulomb_soften`,
  - `slater_alpha`.
- Treat `lj_alpha <= 0` as LJ softening disabled and `slater_alpha <= 0` as Coulomb softening disabled.
- Analyze bilayer-only observables:
  - lipid `PO4` lateral diffusion relative to bilayer COM as the main fluidity signal,
  - reciprocal diffusion as a simple viscosity proxy because the repo does not already expose a direct membrane-viscosity observable.
- Report converted physical units using the existing calibration assumption:
  - `40 ps` per integrator step,
  - production `dt` taken from the staged run settings.
- Allow an optional target diffusion in physical units so the workflow can rank candidate softening settings automatically when the user provides the desired match point.

## Execution Phases
- [x] Phase 1: Rewrite the workflow tracker files to reflect the bilayer-only softening goal and the user correction.
- [x] Phase 2: Replace the workflow runner so it stages bilayer-only tasks over `(lj_alpha, slater_alpha, replicate)` instead of `interface_scale`.
- [x] Phase 3: Replace the analysis path with bilayer-only diffusion / viscosity-proxy summaries and recommendation output.
- [x] Phase 4: Update wrapper scripts and README for the new parameter surface and outputs.
- [x] Phase 5: Run static verification and at least a reduced smoke path, then record review notes.

## Known Errors / Blockers
- There is no direct viscosity estimator already wired into this repo for the dry-MARTINI bilayer workflow.
- The replacement workflow should therefore expose:
  - raw lateral diffusion,
  - converted physical diffusion units,
  - and a reciprocal-diffusion viscosity proxy,
  while keeping the recommendation logic explicit about that limitation.

## Review
- Implemented:
  - bilayer-only stage orchestration in `workflow.py`,
  - production-stage softening attr injection for `lj_soften` / `lj_soften_alpha` / `coulomb_soften` / `slater_alpha`,
  - bilayer-only `PO4` diffusion analysis with physical-unit conversion and viscosity proxy,
  - recommendation assembly output,
  - updated wrapper args and README.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced local smoke run under `/tmp/hybrid_interface_softening_smoke` with:
    - `lj_alpha = 0.05`
    - `slater_alpha = 0.5`
    - `replicates = 1`
    - `MIN_60_MAX_ITER = 1`
    - `MIN_61_MAX_ITER = 1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 40`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - verified the smoke task wrote:
    - `results/tasks/lj0p05_coul0p5_r01.json`
    - `assembled/task_results.csv`
    - `assembled/condition_summary.csv`
    - `assembled/summary.json`
  - verified staged `bilayer.stage_7.0.up` attrs:
    - `lj_soften = 1`
    - `lj_soften_alpha = 0.05`
    - `coulomb_soften = 1`
    - `slater_alpha = 0.5`
  - reduced analysis smoke run wrote:
    - `analysis/results/tasks/lj0p05_coul0p5_r01.json`
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/recommendation_summary.json`
    - `analysis/assembled/summary.json`
  - verified no-submit Slurm staging wrote:
    - `slurm/round_manifest.json`
    - `slurm/run_array.sbatch`
    - `slurm/collect_results.sbatch`
    - `analysis/slurm/round_manifest.json`
    - `analysis/slurm/analyze_array.sbatch`
    - `analysis/slurm/collect_analysis.sbatch`
