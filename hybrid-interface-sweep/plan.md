# Hybrid Interface Softening Sweep

## Project Goal
- Correct `hybrid-interface-sweep/` so it calibrates a bilayer-only scalar softening factor `f` for later interface use.
- Sweep one simple interaction-strength factor rather than softened-potential shape parameters.
- Keep all changes inside `hybrid-interface-sweep/`; do not modify `example/16.MARTINI/`.

## Architecture & Key Decisions
- Use the existing bilayer preparation surface in `py/martini_prepare_system.py --mode bilayer`.
- Keep the bilayer stage ladder already implemented in `hybrid-interface-sweep/workflow.py`:
  - `6.0` / `6.1` minimization,
  - `6.2` / `6.3` softened equilibration from the existing bilayer prep path,
  - `6.4` / `6.5` / `6.6` hard equilibration,
  - `7.0` production.
- Apply the sweep only to the staged production `7.0` Hamiltonian inside the sweep workflow by rewriting `input/potential/martini_potential/coefficients`:
  - `epsilon *= f` for LJ,
  - `q_i *= sqrt(f)` and `q_j *= sqrt(f)` for Coulomb,
  - this is equivalent to multiplying evaluated LJ and Coulomb pair energies and forces by `f`.
- Disable production-stage `lj_soften` / `coulomb_soften` attrs in the rewritten `7.0` file so the sweep does not mix two different control semantics.
- Treat `f` as a true softening factor:
  - finite,
  - within `[0, 1]`,
  - grouped and reported as a single scalar sweep axis.
- Extend the default sweep below `0.6` because the current downloaded grid never reached the provisional `40 ps/step` target proxy.
- Keep bilayer-only observables unchanged:
  - lipid `PO4` lateral diffusion relative to bilayer COM,
  - reciprocal diffusion as an explicitly labeled viscosity proxy,
  - optional target diffusion ranking in physical units.
- Bump workflow manifest/result schemas so stale `(lj_alpha, slater_alpha)` manifests fail loudly instead of being reused accidentally.

## Execution Phases
- [x] Phase 1: Update the tracker files and design notes to the scalar-factor semantics and the “do not touch `example/16.MARTINI/`” correction.
- [x] Phase 2: Rewrite `workflow.py` to sweep a single `interaction_scale` and patch staged `7.0` coefficients instead of production softening attrs.
- [x] Phase 3: Update local/Slurm wrappers, assembled summaries, analysis grouping, and README to the single-factor surface.
- [x] Phase 4: Run static verification and a reduced smoke path on the corrected workflow, then record review notes.

## Known Errors / Blockers
- Old run directories created with the prior `(lj_alpha, slater_alpha)` schema are intentionally incompatible with the corrected workflow.
- The workflow still reports a diffusion-derived viscosity proxy rather than a direct membrane viscosity observable because the repo does not already expose one for this bilayer path.

## Review
- Implemented:
  - single-axis `interaction_scale` sweep manifest/task generation,
  - staged production `7.0` MARTINI coefficient rewrite with `epsilon *= f` and `q *= sqrt(f)`,
  - production-stage `lj_soften` / `coulomb_soften` disablement,
  - assembled run/analysis summaries grouped by `interaction_scale`,
  - updated local/Slurm wrappers and README for `HYBRID_SWEEP_INTERACTION_SCALES`,
  - schema bumps so stale `(lj_alpha, slater_alpha)` manifests are rejected.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced local smoke run under `/tmp/hybrid_interface_scale_smoke_v2` with:
    - `interaction_scale = 0.85`
    - `replicates = 1`
    - `MIN_60_MAX_ITER = 1`
    - `MIN_61_MAX_ITER = 1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 40`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - direct staged HDF5 verification on `bilayer.stage_7.0.prepared.up` versus `bilayer.stage_7.0.up`:
    - `epsilon` ratio = `0.85`
    - charge ratios = `0.921954 = sqrt(0.85)`
    - `lj_soften = 0`
    - `coulomb_soften = 0`
  - reduced local analysis run wrote:
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/recommendation_summary.json`
    - `analysis/assembled/summary.json`
  - no-submit Slurm staging wrote:
    - `slurm/round_manifest.json`
    - `slurm/run_array.sbatch`
    - `slurm/collect_results.sbatch`
    - `analysis/slurm/round_manifest.json`
    - `analysis/slurm/analyze_array.sbatch`
    - `analysis/slurm/collect_analysis.sbatch`
- Downloaded scalar-factor analysis review:
  - assembled analysis is complete for `5` scales with `3 / 3` replicates each,
  - best tested in-grid condition is `interaction_scale = 0.6`,
  - the tested range still does not reach the provisional `40 ps/step` target proxy near `2.892 um^2/s`.
- Default sweep extension:
  - widened the default scale list to `1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1`,
  - so a fresh default run now probes materially lower interaction strength without requiring manual override.
  - verified by `init-run` on `/tmp/hybrid_interface_scale_defaults_check_v3`, which produced `36` tasks across the widened scale list.
