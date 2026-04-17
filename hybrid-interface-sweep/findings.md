# Findings

## 2026-04-17 (Per-Scale Overlay Presentation)
- The fitted `interface_scale_vs_rmsf_difference.png` plot is too weak to be the main presentation artifact for this workflow.
- For this calibration, the user wants to inspect direct reference-vs-hybrid RMSF overlays and choose the scale visually.
- The analysis should therefore render one overlay for every sampled stable scale, not only one “best” overlay.
- The workflow defaults were shifted to a high-end exploratory range because the current overlay at `0.85` still sits above the reference amplitude:
  - `0.825, 0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.00, 1.025, 1.05, 1.075, 1.10, 1.125, 1.15, 1.175, 1.20`
- For this exploratory pass, one long hybrid replicate per scale is a better use of compute than many noisy replicates across the same narrow range.

## 2026-04-17 (Lesson: Match The Primary Output To The Real Decision)
- When the user needs to convince others with a figure, a fitted aggregate trend is not automatically the right artifact.
- The primary output should match the actual decision rule.
- In this workflow, that means direct per-scale RMSF overlays are primary and derived trend plots are secondary.

## 2026-04-17 (Lesson: Leave Range Headroom When The Boundary Point Is Still On One Side)
- When the highest sampled scale is still clearly above the reference amplitude, stop treating that edge as the practical cap.
- Extend the exploratory range beyond the current boundary instead of overcommitting to denser sampling below it.

## 2026-04-17 (Embedded-Region Analysis Target)
- The current full-protein per-residue RMSF fit is not aligned with the intended scientific question.
- The calibration should ignore solution-exposed residues and focus only on the membrane-embedded region.
- Detailed per-residue RMSF shape inside that region is secondary.
- The primary fit metric is now the overall fluctuation amplitude of the embedded region:
  - implemented as mean framewise embedded-region RMSD to the embedded-region mean structure.

## 2026-04-17 (Reference-Driven Embedded Selector)
- The cleanest shared selector is to derive the membrane-embedded residue subset from the reference method only, then apply that same residue subset to every hybrid task.
- The selector now uses the stored membrane geometry in the reference stage file:
  - `input/potential/cb_membrane_potential` or `input/potential/cb_surf_membrane_potential`
  - `half_thickness`
  - `use_curvature`
  - `curvature_radius`
  - `curvature_sign`
- For curved reference membranes, the depth coordinate should follow the engine definition rather than plain `z`:
  - `depth = curvature_sign * (|xyz - center| - curvature_radius)`
- Residues are selected by reference occupancy after burn-in, with default threshold:
  - `HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN = 0.50`

## 2026-04-17 (Embedded-Region RMSD Metric Verification)
- A synthetic end-to-end analysis under `/private/tmp/hybrid_interface_embedded_metric_smoke` verified the new behavior.
- The derived embedded subset was exactly the intended membrane region:
  - residues `2, 3, 4, 5`
- In that synthetic test:
  - one hybrid condition matched the reference embedded-region fluctuation amplitude while making the solution residues much noisier,
  - another hybrid condition increased the embedded-region fluctuation amplitude.
- The assembled analysis correctly chose the amplitude-matched condition:
  - `best_sampled_interface_scale = 0.60`
  - `metric = condition_embedded_region_rmsd_delta_vs_reference_angstrom`
- This confirms the new fit is driven by embedded-region RMSD amplitude and is insensitive to solution-region noise.

## 2026-04-17 (Downloaded `default/` Bundle Compatibility)
- The downloaded `hybrid-interface-sweep/default/` tree is sufficient for a local rerun of the current analysis even without the raw `tasks/` directory.
- Why it works:
  - `default/sweep_manifest.json` is present,
  - `default/analysis/analysis_manifest.json` is already `hybrid_interface_rmsf_sweep_analysis_manifest_v2`,
  - `default/analysis/results/tasks/*.json` are already `hybrid_interface_rmsf_sweep_analysis_result_v2`.
- Consequence:
  - local `assemble-analysis` can be rerun directly from the downloaded analysis-task JSONs,
  - the missing raw `stage_7.0.up` files do not block final assembly in this specific bundle.

## 2026-04-17 (Downloaded `default/` Bundle Result)
- Local assembly of the downloaded `default/` bundle completed successfully under:
  - `hybrid-interface-sweep/default/analysis/assembled/`
- Coverage after assembly:
  - `59` analysis tasks total,
  - `44` successful task analyses,
  - `15` analysis errors,
  - `23` additional hybrid tasks filtered by the protein-stability gate,
  - `17` stable hybrid tasks retained,
  - `9` hybrid conditions used for fitting out of `11` sampled scales.
- Embedded reference region:
  - residues `7 -> 24`
  - `18` residues total
- Best sampled scale by the current metric:
  - `interface_scale = 0.675`
  - `condition_embedded_region_rmsd_delta_vs_reference_angstrom = 0.1441`
  - but only `1 / 4` stable replicate remained.
- Trendline result:
  - fitted recommendation `interface_scale = 0.841`
  - nearest sampled stable overlay scale `= 0.85`
  - fit quality is very weak:
    - `R^2 = 0.053`
- Interpretation:
  - the assembled analysis is valid and current,
  - but the fitted optimum is not persuasive because replicate coverage is sparse and the trend fit is poor.

## 2026-04-14 (Current Workflow Target)
- `hybrid-interface-sweep/` is no longer a bilayer-only diffusion or softening workflow.
- The new target is an interface-scale calibration workflow for `1rkl`.
- The desired calibration observable is protein backbone RMSF:
  - reference method: `example/08.MembraneSimulation`
  - swept method: `example/16.MARTINI`

## 2026-04-14 (Reference Surface: Example 08)
- The relevant `1rkl` example-08 workflow is `example/08.MembraneSimulation/3.fixed_curvature.run.py`.
- Its reference settings are materially different from the generic membrane example:
  - `pdb_id = 1rkl`
  - `membrane_thickness = 24.8`
  - `use_curvature = True`
  - `curvature_radius = 120.0`
  - `curvature_sign = 1`
  - `temperature = 0.80`
  - `exchange = False`
- Because `exchange = False`, running each reference replicate independently is equivalent to the bundled replica launch used in the example script.
- The example-08 config surface for `1rkl` exposes the same 31-residue backbone triplet structure needed for RMSF comparison:
  - `input/potential/affine_alignment/atoms` has shape `(31, 3)`.

## 2026-04-14 (Hybrid Surface: Example 16)
- The relevant hybrid control already exists in `example/16.MARTINI/run_sim_1rkl.sh` and in the engine:
  - env var `PROTEIN_ENV_INTERFACE_SCALE`
  - HDF5 attr `input/hybrid_control/protein_env_interface_scale`
- The hybrid stage-7 output for `1rkl` already exposes a directly comparable 31-residue backbone map:
  - `input/hybrid_bb_map/aa_atom_index` has shape `(31, 3)`.
- This supports a shared analysis surface:
  - rigid alignment over mapped backbone carrier atoms,
  - per-residue RMSF from the mapped CA-like atom.

## 2026-04-14 (Runtime Path Hazard)
- `source.sh` in this checkout currently exports:
  - `UPSIDE_HOME=/Users/yinhan/Documents/upside2-md-martini`
- That path does not exist locally.
- `example/16.MARTINI/run_sim_1rkl.sh` currently re-sources `source.sh`, so it will clobber any caller-provided `UPSIDE_HOME` unless the script is patched to preserve it.
- Consequence:
  - the new sweep should set repo-local runtime paths itself,
  - and the example-16 runner must honor an already-set `UPSIDE_HOME`.

## 2026-04-14 (Wrapper Python Resolution)
- Relying on shell activation state alone is not robust enough for this workflow.
- In the first smoke attempt, the wrapper fell back to a system `python3` that imported `h5py` from Homebrew but had no `numpy` available in that interpreter.
- The shell wrappers therefore need to prefer repo-local `.venv/bin/python3` directly before falling back to `VIRTUAL_ENV` or system `python3`.

## 2026-04-14 (Reduced RMSF Smoke Verification)
- A reduced end-to-end smoke run under `/tmp/hybrid_interface_rmsf_smoke` completed successfully with:
  - `reference_replicates = 1`
  - `hybrid_replicates = 1`
  - `interface_scale = 0.85`
  - shortened reference and hybrid stage lengths for verification only.
- The hybrid task reached and saved:
  - `tasks/scale0p85_r01/run/checkpoints/1rkl.stage_7.0.up`
- The staged hybrid control in that file carries the requested scale:
  - `input/hybrid_control/protein_env_interface_scale = 0.85`
- The reduced analysis also completed successfully and wrote the full assembled RMSF output set:
  - `task_results.csv`
  - `residue_rmsf_profiles.csv`
  - `reference_profile.csv`
  - `condition_profiles.csv`
  - `condition_summary.csv`
  - `trendline_points.csv`
  - `recommendation_summary.json`
  - `summary.json`
- With only one sampled scale in the smoke tree, the trend-line fit correctly degrades to a degree-0 constant fit and recommends that same sampled scale.

## 2026-04-15 (Protein-Only RMSF Audit)
- The hybrid `stage_7.0.up` file does dump the full mixed system:
  - smoke `output/pos` shape is `(20, 1, 4323, 3)`.
- The RMSF analysis does not use all `4323` particles.
- On the verified smoke artifact, it selects the residue map from:
  - `input/potential/affine_alignment/atoms`
  - shape `(31, 3)`
  - total selected atoms `= 93`
- The selected atoms are protein-only backbone carrier atoms:
  - `selected_particle_classes = ['PROTEINAA']`
  - `selected_atom_roles = ['C', 'CA', 'N']`
- The analysis now validates this explicitly and will fail if the selected RMSF atoms are not protein-only.

## 2026-04-15 (Slurm Analysis Failure: Kabsch SVD Non-Convergence)
- The remote analysis failure
  - `numpy.linalg.LinAlgError: SVD did not converge`
  from `_align_points_kabsch(...)` indicates the selected protein-backbone frame data can become non-finite for at least some tasks.
- On a `3 x 3` covariance matrix, SVD non-convergence is far more likely to come from `NaN`/`Inf` contamination than from a genuinely hard linear-algebra problem.
- The RMSF analysis is now hardened to:
  - keep the protein-only selection check,
  - drop non-finite selected backbone frames after burn-in,
  - require at least five finite frames to remain,
  - retry the 3x3 SVD with a tiny diagonal jitter if the first call fails.
- The per-task analysis payload now records:
  - `n_frames_dropped_nonfinite`
- A synthetic HDF5 with one NaN protein frame now succeeds and reports:
  - `n_frames_dropped_nonfinite = 1`
  - `n_frames_used = 5`

## 2026-04-15 (Protein Stability Filter For Destroyed Hybrid Trajectories)
- Some hybrid `1rkl` runs can destroy the protein badly enough that RMSF matching is no longer meaningful.
- The analysis now treats stability as a per-hybrid-trajectory gate relative to the reference ensemble.
- The current filter compares each hybrid task against reference baselines using four ratios:
  - mean residue RMSF,
  - max residue RMSF,
  - mean CA radius of gyration,
  - mean CA span.
- Default filter limits are:
  - `stability_mean_rmsf_ratio_max = 3.0`
  - `stability_max_rmsf_ratio_max = 3.0`
  - `stability_ca_rg_ratio_max = 1.75`
  - `stability_ca_span_ratio_max = 1.75`
- Unstable hybrid tasks remain visible in:
  - `task_results.csv`
  - `failed_tasks.csv`
  with explicit ratio values and filter reasons.
- Unstable hybrid tasks are excluded from:
  - `condition_summary.csv` metrics,
  - `condition_profiles.csv`,
  - trend-line fitting,
  - recommended scale selection.
- The shortened local smoke hybrid run is an example destroyed trajectory under this filter:
  - `mean_rmsf_ratio = 13.09`
  - `max_rmsf_ratio = 52.65`
  - `ca_rg_ratio = 3.02`
  - `ca_span_ratio = 5.37`

## 2026-04-15 (Existing Run Tree Compatibility)
- Existing completed run trees may have been initialized before the stability-filter keys were added to `analysis_settings`.
- Re-analysis therefore must not assume those manifest keys already exist.
- The analysis layer now backfills missing settings with current defaults when:
  - creating `analysis_manifest.json`,
  - running analysis tasks,
  - assembling analysis outputs.

## 2026-04-15 (No-Stable Analysis State)
- A legitimate outcome of the new filter is that zero hybrid trajectories remain stable enough for fitting.
- That case should not crash analysis or leave stale old recommendation files behind.
- The assembled outputs now represent this state explicitly:
  - `recommendation_summary.json` gets `analysis_status = no_stable_hybrid_trajectories`,
  - recommendation and trend-line scalars become `null`,
  - `trendline_points.csv` is written with headers only,
  - `summary.json` records `trendline_available = 0`.

## 2026-04-15 (Downloaded Task Bundle Shape)
- The downloaded `hybrid-interface-sweep/tasks/` bundle is not the raw run tree.
- It contains per-task analysis result JSON files directly:
  - `4` successful reference task results,
  - `14` successful hybrid task results,
  - `4` hybrid analysis failures.
- Because the raw remote `stage_7.0.up` files are not present locally, this rerun is an assembled analysis from downloaded task-result JSONs rather than a fresh frame extraction from HDF5 checkpoints.
- The bundle also uses mixed analysis payload versions:
  - most successful task JSONs are older results with RMSF-only fields,
  - only newer results include geometry fields like `max_rmsf_angstrom`, `ca_radius_of_gyration_*`, and `ca_span_*`.

## 2026-04-15 (Legacy Task-Result Compatibility)
- The current assembled-analysis code must tolerate downloaded task-result bundles that mix old and new analysis payload schemas.
- For legacy successful payloads, the assembler can safely derive:
  - `max_rmsf_angstrom` from the saved `residue_rmsf_angstrom` profile.
- Geometry-based stability fields cannot be reconstructed from legacy RMSF-only task JSONs, so those checks must be skipped when the required values are absent.
- With this compatibility layer, the downloaded bundle assembled successfully and still filtered obvious destroyed trajectories by the available RMSF ratios.

## 2026-04-15 (Current Stable Sweep Region)
- The latest assembled downloaded task bundle supports a practical stable default sweep region of about:
  - `interface_scale = 0.40 -> 1.00`
- In that bundle:
  - `0.40` had full `3 / 3` stable replicate coverage,
  - `0.55`, `0.70`, `0.85`, and `1.00` all retained at least one stable replicate,
  - `0.25` is closer to the instability boundary.
- A denser default grid inside `0.40 -> 1.00` is therefore a better next calibration surface than continuing to sample the lower-`scale` destabilized branch sparsely.

## 2026-04-15 (Longer RMSF Calibration Defaults)
- To make the scaling-factor trend easier to defend, the workflow defaults should reduce trajectory-noise rather than relying only on more sampled scales.
- The workflow now lengthens the default runs used for new sweeps:
  - reference `REFERENCE_DURATION = 200001`,
  - hybrid `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1000`,
  - hybrid `PROD_70_NSTEPS = 50000`,
  - hybrid `EQ_FRAME_STEPS = 250`,
  - hybrid `PROD_FRAME_STEPS = 100`.

## 2026-04-15 (Rendered Plot Artifact)
- Plot-ready CSVs are useful for notebooks, but they are not sufficient when the goal is to persuade others quickly that the scale calibration is working.
- The assembled analysis should therefore render a figure directly.
- The workflow now writes:
  - `analysis/assembled/interface_scale_vs_rmsf_difference.png`
  - `analysis/assembled/interface_scale_vs_rmsf_difference.svg`
- The rendered figure shows:
  - stable condition points,
  - fit trend line,
  - best sampled point,
  - stable/completed replicate counts per scale.
- Matplotlib on this machine needs writable local cache directories.
  - Setting `MPLCONFIGDIR` and `XDG_CACHE_HOME` inside the plotting step avoids noisy font/cache warnings in analysis logs.

## 2026-04-16 (Downloaded Long-Sweep Run Bundle)
- The latest downloaded bundle includes:
  - `results/tasks/*.json`
  - `analysis/results/tasks/*.json`
  - `slurm/`
  - `analysis/slurm/`
- It does not include the original `sweep_manifest.json`, but that manifest is reconstructable from the downloaded task-result JSONs.
- Downloaded run coverage:
  - `43` run-task JSONs total,
  - `42` successful runs,
  - `1` failed run.
- Downloaded analysis coverage:
  - `42` analysis-task JSONs total,
  - `37` successful analyses,
  - `5` analysis failures.

## 2026-04-16 (Latest Run Failure Root Cause)
- The only run-stage failure in the downloaded long sweep is:
  - `scale0p75_r02`
- It is not an analysis bug.
- The stage-7 dynamics itself blew up:
  - `Rg` rose into the `1e5 -> 1e6 A` range,
  - `martini_potential` rose into about `1e15`,
  - the C++ engine then segfaulted inside spline / interaction evaluation during `Stage 7.0`.
- This failed task should be treated as a genuine unstable trajectory, not retried in analysis.

## 2026-04-16 (Latest Analysis Failure Root Cause)
- The five analysis-task failures are all the same failure mode:
  - `Too few finite backbone frames remain after filtering for RMSF analysis: 0`
- Affected tasks:
  - `scale0p5_r02`
  - `scale0p55_r03`
  - `scale0p6_r01`
  - `scale0p9_r02`
  - `scale0p95_r03`
- Interpretation:
  - these runs completed far enough to write `stage_7.0.up`,
  - but every post-burn-in selected protein-backbone frame was non-finite,
  - so the analysis correctly rejected them as unusable.

## 2026-04-16 (Latest Assembled Validity)
- The latest downloaded long sweep is still usable for partial calibration analysis after filtering failed and unstable tasks.
- Assembled counts after local rerun:
  - `37` successful task analyses,
  - `9` additional successful hybrid analyses filtered as unstable,
  - `24` stable hybrid analyses retained,
  - `13` interface-scale conditions still represented in the fitted table.
- Important limitation:
  - stable replicate coverage is uneven,
  - several scales rely on only `1` or `2` stable replicates,
  - the quadratic trend fit is weak:
    - `R^2 = 0.187`
- Consequence:
  - the sweep is valid enough to inspect and plot,
  - but the trend line is not strong enough to treat the fitted minimum as a definitive calibration result.

## 2026-04-16 (Latest Practical Ranking)
- Raw best sampled condition from the latest assembled bundle:
  - `interface_scale = 0.60`
  - `condition_profile_rmse = 0.350 A`
  - but only `1 / 3` stable replicates remained.
- Trendline minimum from the latest bundle:
  - `interface_scale = 0.583`
  - fitted RMSE `= 0.488 A`
  - fit quality `R^2 = 0.187`
- More defensible branches by stable coverage:
  - `0.75`:
    - `RMSE = 0.364 A`
    - `2 / 2` stable completed replicates
    - one separate run crashed in stage 7
  - `0.65`:
    - `RMSE = 0.421 A`
    - full `3 / 3` stable replicate coverage
  - `0.80`:
    - `RMSE = 0.429 A`
    - full `3 / 3` stable replicate coverage

## 2026-04-16 (Targeted Confirmation Sweep Defaults)
- The broad `0.40 -> 1.00` default sweep is no longer the best default calibration surface.
- The latest assembled long sweep supports a narrower confirmation window centered on the current best-supported region:
  - `interface_scale = 0.55, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725, 0.75, 0.775, 0.80, 0.85`
- Replicate count matters more than further widening the scale range at this stage.
- New default hybrid replicate count:
  - `hybrid_replicates = 5`
- A fresh `init-run` now produces:
  - `11` sampled scales,
  - `55` hybrid tasks total,
  - plus the unchanged reference task family.

## 2026-04-16 (Best-Scale RMSF Overlay Plot)
- The assembled analysis should not stop at the scale-vs-RMSE trend plot when the goal is to convince others that the interface scaling is physically meaningful.
- The analysis now renders a second presentation plot:
  - `analysis/assembled/best_interface_scale_rmsf_vs_reference.png`
  - `analysis/assembled/best_interface_scale_rmsf_vs_reference.svg`
- This plot overlays:
  - the reference per-residue RMSF mean and spread,
  - the selected hybrid condition per-residue RMSF mean and spread.
- Selection rule:
  - if the fitted recommended scale is also a sampled stable scale, use it directly,
  - otherwise use the nearest stable sampled scale to the fitted recommendation,
  - if no fitted recommendation is available, fall back to the best sampled stable scale.
- On the current assembled dataset:
  - `trendline_recommended_interface_scale = 0.583`
  - `profile_comparison_interface_scale = 0.60`
  - `profile_comparison_selection_basis = nearest_stable_sample_to_trendline_recommendation`

## 2026-04-14 (Analysis Surface)
- There is no existing RMSF analysis helper in the repo for this calibration.
- The new analysis should therefore compute and store:
  - per-task residue RMSF profiles,
  - aggregate reference RMSF profile,
  - per-hybrid-task RMSF mismatch against that reference,
  - condition-level summaries over `interface_scale`,
  - fit samples for a trend line used to choose the recommended scale.
- Plot-ready outputs should be explicit CSV files rather than requiring notebook-only postprocessing.

## 2026-04-14 (Current Workflow Target)
- `hybrid-interface-sweep/` must be a bilayer-only calibration workflow.
- The requested fix is to use softened potentials for both LJ and Coulomb in production-stage `7.0`.
- The requested change belongs only in `hybrid-interface-sweep/`.
  - `example/16.MARTINI/` is not part of this correction.

## 2026-04-14 (Scientific Direction)
- The current scalar `interaction_scale` workflow is the wrong control surface for this task.
- The user’s concern is that plain coefficient scaling leaves the LJ and Coulomb singular shapes intact near `r -> 0`.
- The restored sweep should therefore operate on the potential shape controls already supported by the engine:
  - `lj_soften`
  - `lj_soften_alpha`
  - `coulomb_soften`
  - `slater_alpha`

## 2026-04-14 (Stage-7 Rewrite Surface)
- The bilayer prep path already writes the required MARTINI node for production-stage patching:
  - `/input/potential/martini_potential`
- Rewriting only the staged `7.0` attrs is sufficient to change the production nonbonded forces seen by the bilayer-only diffusion measurement.
- The corrected workflow should record the applied production softening settings in task results so the analysis tree remains self-describing.

## 2026-04-14 (Sweep Surface And Compatibility)
- The workflow should sweep two non-negative controls:
  - `lj_alpha`
  - `slater_alpha`
- The earlier softening workflow grid remains a reasonable default restoration target:
  - `lj_alpha = 0.0, 0.025, 0.05, 0.10, 0.20`
  - `slater_alpha = 0.0, 0.25, 0.50, 1.00, 2.00`
- Scalar-manifest directories must be rejected explicitly by a schema bump so stale runs are not mixed into the restored softening workflow.

## 2026-04-14 (Analysis Surface)
- The existing bilayer-only diffusion analysis remains the right observable family:
  - `PO4` lateral diffusion relative to bilayer COM,
  - physical-unit diffusion conversion via the configured `ps/step` assumption,
  - reciprocal-diffusion viscosity proxy,
  - optional target-diffusion ranking.
- Only the grouping axis and recommendation surface need to change:
  - from `interaction_scale`
  - back to `(lj_alpha, slater_alpha)`.

## 2026-04-14 (Runtime Execution Surface)
- The workflow was still failing locally even after the semantic restoration because the preparation helpers launched bare `python3` subprocesses.
- Since the wrappers and generated Slurm scripts already choose a specific interpreter, the workflow subprocesses must inherit that same interpreter.
- Using `sys.executable` inside `workflow.py` keeps:
  - local wrapper runs,
  - direct workflow invocations,
  - and Slurm array tasks
  on one consistent Python environment with the required packages.

## 2026-04-14 (Downloaded Softening Analysis Rerun)
- The local repo copy currently contains the downloaded `hybrid-interface-sweep/analysis/` tree only:
  - task-level analysis JSON files,
  - assembled CSV/JSON summaries,
  - Slurm logs and manifests.
- It does not contain the corresponding full sweep task tree locally, so this rerun is a validation/review of the downloaded analysis artifacts rather than a fresh extraction from raw `stage_7.0.up` files.
- The downloaded softening analysis is numerically complete:
  - `75 / 75` analysis task JSON files succeeded,
  - `0` failed tasks,
  - all `25` conditions have full `3 / 3` replicate coverage.
- Internal consistency is strong across the full sweep:
  - every task uses `160` post-burn-in frames,
  - every task uses `102` `PO4` beads,
  - all diffusion values are finite,
  - task-level `PO4` fit quality spans `R^2 = 0.9780 -> 0.9927`.
- The saved recommendation from the downloaded artifacts is supported by the condition table:
  - best mean diffusion: `lj_alpha = 0.1`, `slater_alpha = 0.5`, `D = 0.8621 um^2/s`,
  - but that branch is somewhat noisier than the nearby alternatives (`CV = 0.082`).
- Two practical neighboring branches remain important for interpretation:
  - lower-perturbation high-fluidity branch: `(lj_alpha, slater_alpha) = (0.0, 1.0)` with `D = 0.8584 um^2/s`, `CV = 0.047`,
  - highest-stability near-top branch: `(0.1, 2.0)` with `D = 0.8510 um^2/s`, `CV = 0.013`, `min R^2 = 0.9888`.
- Relative to the same provisional `40 ps/step` target proxy from `hybrid_timescale.md`:
  - target at workflow temperature `T = 0.8647` remains about `2.892 um^2/s`,
  - the downloaded softening sweep reaches only `0.733 -> 0.862 um^2/s`,
  - so no tested condition demonstrates a target match.

## Lessons
- When the user says the solution region should be ignored, the fit metric must be redefined around a membrane-region selector instead of merely relabeling the old full-protein RMSF outputs.
- When the user says per-residue shape does not matter, do not keep fitting profile RMSE; use a scalar fluctuation-amplitude metric instead.
- When the user points out that scaling leaves a singular pair-potential shape unchanged, restore the shape-control surface instead of trying to force the scalar workflow to fit.
- When a workflow semantics change invalidates old run trees, bump the manifest/result schemas immediately so reuse fails loudly.
- When a workflow script is launched from a managed environment, subprocess Python calls must use `sys.executable` rather than bare `python3`.
- When the user provides only a downloaded analysis tree, state explicitly whether the rerun is a fresh recomputation from raw checkpoints or a validation/review of the downloaded artifacts.
- When the user redirects the scientific target entirely, rewrite the workflow around the new observable instead of trying to preserve the previous analysis pipeline.
- When the user asks for a convincing calibration plot, add a direct reference-vs-selected-condition RMSF overlay instead of assuming the trend plot alone is enough.
- When an example script re-sources environment setup, preserve caller-supplied runtime paths if the sweep depends on running that script from the current checkout.
- When wrappers are meant to use a repo-local virtualenv, prefer the concrete `.venv/bin/python3` path before trusting the ambient shell state.
- When a workflow output contains both target and environment particles, record and validate the exact atom-class subset used by the analysis instead of leaving that assumption implicit.
- When a cluster-only `SVD did not converge` failure appears in a small rigid-alignment solve, assume non-finite frame data first and harden the analysis to filter bad frames instead of treating it as a pure numerical-method issue.
- When the user says some trajectories destroy the protein, exclude those trajectories from condition averages and trend fitting instead of letting them bias the calibration surface.
- When adding new analysis settings to a persisted manifest, backfill defaults for older run trees so re-analysis still works on already-completed sweeps.
- When a filtered-data workflow leaves no valid samples for fitting, emit an explicit no-fit result state rather than aborting after partially overwriting assembled outputs.
- When the user provides a downloaded `tasks/` bundle, inspect whether it is raw task directories or already-computed task-result JSONs before assuming the local rerun can read HDF5 checkpoints.
- When downloaded task-result JSONs mix legacy and current schemas, normalize them during assembly instead of requiring the user to regenerate the whole analysis remotely.
- When a user wants a scale-calibration result to be persuasive rather than merely inspectable, emit a rendered PNG/SVG figure directly from assembled analysis instead of stopping at CSV tables.
- When Matplotlib runs in constrained environments, point its config and cache directories at writable analysis-local paths so plotting does not clutter Slurm logs with cache warnings.
- When a downloaded run bundle is missing the original manifest but includes per-task result JSONs, reconstruct the manifest locally instead of blocking on a rerun.
- When the fitted trend line has very low `R^2`, report that explicitly and prefer discussing best-sampled and best-covered branches separately rather than overselling the fitted optimum.
