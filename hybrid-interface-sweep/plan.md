# Hybrid Interface RMSF Calibration Sweep

## Project Goal
- Replace `hybrid-interface-sweep/` with an interface-scale calibration workflow for `1rkl`.
- Use the implicit-membrane `example/08.MembraneSimulation` method as the RMSF reference.
- Use the hybrid dry-MARTINI `example/16.MARTINI` method as the swept method.
- Sweep the Upside/dry-MARTINI interface scaling factor and choose it from an RMSF-error trend line.
- Keep local and Slurm execution paths working, including a post-run analysis phase.

## Architecture & Key Decisions
- Revert the sweep surface back to `PROTEIN_ENV_INTERFACE_SCALE`; the softening workflow is no longer the right control for this task.
- Represent the sweep as two task families inside one manifest:
  - `reference` tasks: run `1rkl` with the example-08 fixed-curvature implicit-membrane method.
  - `hybrid` tasks: run `1rkl` with the example-16 hybrid membrane workflow at one `interface_scale`.
- Do not call the example-08 script directly.
  - `3.fixed_curvature.run.py` shells out to `sbatch` and bundles replicas only for convenience.
  - Because `exchange = False`, independent per-replicate runs are equivalent and fit the sweep/Slurm model better.
- Keep the example-08 reference physics aligned to `3.fixed_curvature.run.py`:
  - `pdb_id = 1rkl`
  - `membrane_thickness = 24.8`
  - `use_curvature = True`
  - `curvature_radius = 120.0`
  - `curvature_sign = 1`
  - `temperature = 0.80`
  - `integrator = mv`
  - `--disable-recentering`
- Keep the hybrid production runner aligned to `example/16.MARTINI/run_sim_1rkl.sh`.
- Preserve caller-owned repo paths for the hybrid runner.
  - Current `source.sh` points `UPSIDE_HOME` at a non-existent checkout.
  - The workflow therefore needs the example-16 shell script to honor an already-set `UPSIDE_HOME`.
- RMSF comparison surface:
  - compare the shared `1rkl` protein backbone only,
  - align frames rigidly in 3D using all mapped backbone carrier atoms,
  - compute per-residue RMSF from the mapped CA-like backbone atom for each residue,
  - aggregate the example-08 reference profile across reference replicates,
  - score each hybrid replicate and condition by RMSF profile mismatch to that reference.
- Analysis outputs must be plot-ready:
  - per-task RMSF summaries,
  - per-residue RMSF profiles,
  - condition-level RMSF error table versus `interface_scale`,
  - fitted trend-line samples,
  - recommendation JSON with both best sampled scale and trend-line-selected scale.
- Fit a simple polynomial trend line over condition-level RMSF RMSE.
  - Default to quadratic when at least three scales are available.
  - Choose the recommended scale from the minimum fitted error within the sampled scale range.

## Execution Phases
- [x] Phase 1: Rewrite tracker files for the new interface-scale RMSF calibration workflow.
- [x] Phase 2: Replace `workflow.py` so it stages reference and hybrid tasks, runs them locally or via Slurm, and records the right metadata.
- [x] Phase 3: Update the example-16 shell runner just enough to preserve caller-supplied repo paths.
- [x] Phase 4: Update wrapper scripts and README for the new run and analysis surfaces.
- [x] Phase 5: Run static verification and reduced smoke tests for both methods plus the assembled trend-line analysis.
- [x] Phase 6: Record review notes, verification results, and caveats.

## Known Errors / Blockers
- `source.sh` currently exports `UPSIDE_HOME=/Users/yinhan/Documents/upside2-md-martini`, but that path does not exist on this machine.
- No pre-existing example-08 output artifact is stored in the repo, so the workflow must generate its own reference trajectories before hybrid RMSF comparison.

## Review
- Implemented:
  - replaced the old softening workflow with a mixed-task `reference + hybrid` interface-scale RMSF sweep,
  - added example-08-style reference task generation and direct Upside execution for `1rkl`,
  - restored hybrid task execution through `example/16.MARTINI/run_sim_1rkl.sh` with per-task `PROTEIN_ENV_INTERFACE_SCALE`,
  - added RMSF analysis that writes residue profiles, condition summaries, fitted trend-line samples, and a recommendation JSON,
  - updated wrapper scripts and README for the new interface-scale run surface,
  - patched `example/16.MARTINI/run_sim_1rkl.sh` so caller-provided `UPSIDE_HOME` survives its local `source.sh`.
- Runtime hardening:
  - all four shell wrappers now prefer repo-local `.venv/bin/python3` directly before falling back to `VIRTUAL_ENV` or system `python3`.
- Verification completed:
  - `source .venv/bin/activate && source source.sh && .venv/bin/python -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - reduced local smoke sweep under `/tmp/hybrid_interface_rmsf_smoke` with:
    - `interface_scale = 0.85`
    - `reference_replicates = 1`
    - `hybrid_replicates = 1`
    - `REFERENCE_DURATION = 20`
    - `REFERENCE_FRAME_INTERVAL = 1`
    - `MIN_60_MAX_ITER = 1`
    - `MIN_61_MAX_ITER = 1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 20`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - reduced analysis run on the same smoke tree wrote:
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/residue_rmsf_profiles.csv`
    - `analysis/assembled/reference_profile.csv`
    - `analysis/assembled/condition_profiles.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/trendline_points.csv`
    - `analysis/assembled/recommendation_summary.json`
    - `analysis/assembled/summary.json`
  - verified hybrid production control on the smoke tree:
    - `input/hybrid_control/protein_env_interface_scale = 0.85`
  - verified no-submit Slurm staging on the same smoke tree:
    - `slurm/round_manifest.json`
    - `slurm/run_array.sbatch`
    - `slurm/collect_results.sbatch`
    - `analysis/slurm/round_manifest.json`
    - `analysis/slurm/analyze_array.sbatch`
    - `analysis/slurm/collect_analysis.sbatch`
