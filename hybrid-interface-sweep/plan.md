# Hybrid Interface Softening Sweep Restoration

## Project Goal
- Restore `hybrid-interface-sweep/` so it sweeps softened nonbonded potentials for both LJ and Coulomb in bilayer-only production runs.
- Remove the current scalar `interaction_scale` rewrite from the workflow.
- Keep all changes inside `hybrid-interface-sweep/`; do not modify `example/16.MARTINI/`.

## Architecture & Key Decisions
- Keep the existing bilayer preparation surface in `py/martini_prepare_system.py --mode bilayer`.
- Keep the bilayer stage ladder already implemented in `hybrid-interface-sweep/workflow.py`:
  - `6.0` / `6.1` minimization,
  - `6.2` / `6.3` softened equilibration from the existing bilayer prep path,
  - `6.4` / `6.5` / `6.6` hard equilibration,
  - `7.0` production.
- Apply the sweep only to the staged production `7.0` Hamiltonian inside the sweep workflow by rewriting `input/potential/martini_potential` attrs:
  - `lj_soften = 1` when `lj_alpha > 0`,
  - `lj_soften_alpha = lj_alpha`,
  - `coulomb_soften = 1` when `slater_alpha > 0`,
  - `slater_alpha = slater_alpha`.
- Treat `lj_alpha` and `slater_alpha` as independent non-negative sweep axes; default to the previously working grid:
  - `lj_alpha = 0.0, 0.025, 0.05, 0.10, 0.20`
  - `slater_alpha = 0.0, 0.25, 0.50, 1.00, 2.00`
- Bump manifest/result schemas so stale scalar-factor manifests fail loudly instead of being reused accidentally.
- Restore assembled run/analysis grouping to `(lj_alpha, slater_alpha)` conditions.

## Execution Phases
- [x] Phase 1: Rewrite the local tracker files to the restored softening semantics and record the user correction.
- [x] Phase 2: Patch `workflow.py` task schemas, CLI surface, and stage-7 mutation logic from scalar scaling back to softened-potential attrs.
- [x] Phase 3: Patch wrappers and README so local/Slurm entrypoints expose `HYBRID_SWEEP_LJ_ALPHAS` and `HYBRID_SWEEP_SLATER_ALPHAS`.
- [x] Phase 4: Run static verification and a reduced smoke path that proves the staged `7.0` file carries the requested softening attrs.
- [x] Phase 5: Record review notes, compatibility caveats, and verification results.

## Known Errors / Blockers
- Base directories initialized with the scalar `interaction_scale` manifest are intentionally incompatible with the restored workflow and must be reinitialized.
- Any previously downloaded scalar-factor analysis is not valid for selecting softened-potential parameters.

## Review
- Implemented:
  - restored two-axis `(lj_alpha, slater_alpha, replicate)` manifest/task generation,
  - restored production-stage `7.0` softening patching via `lj_soften`, `lj_soften_alpha`, `coulomb_soften`, and `slater_alpha`,
  - restored run/analysis assembly grouping and recommendation selection on `(lj_alpha, slater_alpha)`,
  - restored wrapper and README env surfaces to `HYBRID_SWEEP_LJ_ALPHAS` and `HYBRID_SWEEP_SLATER_ALPHAS`,
  - bumped manifest/result schemas to `v3` so stale scalar-factor trees fail loudly.
- Runtime fix:
  - `workflow.py` now launches preparation subprocesses with `sys.executable` instead of bare `python3`, so the prep path stays on the same interpreter selected by the wrapper or Slurm entrypoint.
- Verification completed:
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced local smoke run under `/tmp/hybrid_interface_softening_restore_smoke` with:
    - `lj_alpha = 0.05`
    - `slater_alpha = 0.5`
    - `replicates = 1`
    - `MIN_60_MAX_ITER = 1`
    - `MIN_61_MAX_ITER = 1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 40`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - direct staged HDF5 verification on `bilayer.stage_7.0.up` confirmed:
    - `lj_soften = 1`
    - `lj_soften_alpha = 0.05`
    - `coulomb_soften = 1`
    - `slater_alpha = 0.5`
  - reduced local analysis run wrote:
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/recommendation_summary.json`
    - `analysis/assembled/summary.json`
