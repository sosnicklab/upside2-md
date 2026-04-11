# Findings

## 2026-04-06 (Repo Reuse For Bilayer Diffusion Workflow)
- `py/martini_prepare_system.py` already supports a bilayer-only `prepare -> stage conversion` path.
- The generated `.up` files already contain the metadata needed for direct MSD analysis:
  - `input/atom_names`
  - `input/residue_ids`
  - `input/molecule_ids`
  - `input/particle_class`
  - box dimensions via `output/box` or `input/potential/martini_potential` attrs.
- The bundled `parameters/dryMARTINI/DOPC.pdb` template currently contains:
  - `72` DOPC lipids,
  - `8` `NA`,
  - `8` `CL`.
- Existing reusable Slurm patterns are available in:
  - `SC-training/workflow.py`
  - `SC-training/run_local.sh`
  - `SC-training/submit_remote_round.sh`
  - `/Users/yinhan/Documents/upside2-md/ConDiv_symlay/slurm_round.py`
  - `/Users/yinhan/Documents/upside2-md/ConDiv_symlay/submit_remote_round.sh`

## 2026-04-06 (Workflow Verification Findings)
- In the current bilayer prep path, lipid beads are not labeled as `LIPID` in `input/particle_class`; the DOPC beads appear as `OTHER`, while ions are labeled as `ION`.
- The analyzer therefore needs a fallback lipid selection path based on `PO4`-containing molecules rather than assuming `particle_class == LIPID`.
- Leaflet assignment is more reliable when taken from the early-stage bilayer geometry (`stage_6.0.up`) than from the production-stage first frame, because reduced smoke runs can distort the production reference geometry.
- Wrapped `z`-coordinate splitting is required for leaflet assignment; a raw mean-`z` split is not safe when one leaflet sits across the periodic boundary.
- The deliberately tiny smoke setting used for workflow verification can show large, nonphysical energies. It is adequate for plumbing verification only and should not be used to assess fluidity values.

## 2026-04-06 (Slurm Wrapper Portability Findings)
- A shell wrapper that stores `PYTHON_BIN=python3` before module loading / activation is not stable on the cluster; the eventual `python3` on `PATH` can become the module interpreter instead of the project venv.
- For these workflows, wrapper scripts should resolve `PYTHON_BIN` after venv activation and prefer `${VIRTUAL_ENV}/bin/python3` when available.
- The workflow runtime environment should override `UPSIDE_HOME` to the current repository root rather than inheriting an ambient value from the caller, because `martini_prepare_system_lib.py` reads `os.environ['UPSIDE_HOME']` directly.
- When child commands fail, the workflow should report the log path directly; otherwise the user only sees the command line and has to hunt for the actual error.

## 2026-04-06 (Cluster Python 3.9 Compatibility Finding)
- `py/martini_prepare_system_lib.py` is imported during `init-run`, and that module contained PEP 604 annotations (`Path | None`, `int | None`) without deferred annotation evaluation.
- On Python 3.9 this does not fail at parse time; it fails at import-time annotation evaluation with `TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'`.
- The minimal compatibility fix is `from __future__ import annotations` at the top of `py/martini_prepare_system_lib.py`.

## 2026-04-06 (Checkout Root Preservation Finding)
- `bilayer-lateral-diffusion/workflow.py` originally used `Path(__file__).resolve()` to compute its repo root.
- On the cluster, that can canonicalize `/home/...` into `/project/...` and make the workflow run against a different checkout than the one the user invoked.
- For this workflow, `Path(os.path.abspath(__file__))` is the safer choice because it preserves the launched checkout path without symlink canonicalization.

## 2026-04-06 (Cluster HDF5 Module Finding)
- On the target cluster, `obj/upside` requires the HDF5 runtime module to be loaded explicitly; otherwise Slurm array tasks fail with `libhdf5.so.310: cannot open shared object file`.
- The correct site-specific fix is to load `hdf5/1.14.3` alongside the existing `python`, `cmake`, and `openmpi` module loads.
- This requirement applies to both the wrapper path and the generated Slurm scripts, not just the submit step.

## 2026-04-08 (Stage-7-Only Post-Analysis Findings)
- The existing `_analyze_production()` implementation can analyze a `stage_7.0.up` checkpoint without a stage-6 reference by using the stage-7 file as its own leaflet reference.
- That makes a separate post hoc analysis workflow viable after users delete `6.*` checkpoint files to save space.
- Writing post hoc outputs under `stage7-analysis/` avoids overwriting the original `results/` and `assembled/` trees.
- Packing results for download should include the manifest plus analysis outputs by default; the stage-7 checkpoint files should remain optional because they can dominate archive size.
- Verified directly on one existing `stage_7.0.up` file that `reference_stage_file=None` causes the analyzer to use the stage-7 file itself as `leaflet_reference_file`.

## 2026-04-08 (Downloaded Report Calibration Findings)
- For the completed `stage7-analysis/` tree, all structural counts matched:
  - `525` manifest stage-7 files,
  - `525` successful tasks,
  - `525` task CSV rows,
  - `175` condition rows,
  - `525` per-task JSON files.
- Under the user-provided calibration of `40 ps` per Upside integration step and the workflow step size of `0.01`, one workflow time unit corresponds to `4 ns`.
- That implies:
  - `nm^2/ns = diffusion_mean_nm2_per_time / 4`
  - `um^2/s = diffusion_mean_nm2_per_time * 250`
- Across this scan, mass scaling is the dominant diffusion-control knob:
  - geometric-mean diffusion changes by about `141x` from mass `1.0` to `0.1`,
  - temperature contributes about `2.3x` from `0.7` to `1.1`,
  - damping contributes only about `1.27x` across `4` to `64`.
- The fastest robust conditions in the current grid all use `mass_scale = 0.1`; mass scales `>= 0.25` remain much slower under the `40 ps/timestep` calibration.
- The single fastest mean condition, `(T=0.7, tau=64, mass=0.1)`, is also relatively noisy across replicates (`CV ~ 0.45`), so candidate ranking should consider variability, not just the mean diffusion.
- The assembled `thickness_angstrom_mean` values span both positive and negative signs in this dataset, so thickness should not be used as a validity gate or ranking metric without clarifying its sign convention first.

## 2026-04-08 (Fixed-Temperature Physical Calibration Findings)
- Using the repository unit contract from `AGENTS.md`, the scanned temperatures map to:
  - `0.7 -> 245.4 K`
  - `0.8 -> 280.5 K`
  - `0.9 -> 315.5 K`
  - `1.0 -> 350.6 K`
  - `1.1 -> 385.6 K`
- Using the DOPC supported-bilayer diffusion values reported by Bag et al. 2014 (`298 K: 2.46`, `303 K: 2.77`, `308 K: 3.07`, `313 K: 3.40 um^2/s`), an Arrhenius fit to the rounded table values gives an inferred `E_A ~ 16.66 kJ/mol`.
- Comparing the scan to that fixed-temperature target shows that the bilayer remains too slow at every scanned temperature even at `mass_scale = 0.1`.
- Fitting `log D` versus `log mass_scale` at fixed temperature gives a strong power-law relation with exponent near `-2.1` to `-2.25`, while damping contributes only a weak secondary effect.
- A stable temperature-specific correction therefore looks like:
  - `T_up = 0.7`: `tau ~ 16`, `mass ~ 0.068`
  - `T_up = 0.8`: `tau ~ 4`, `mass ~ 0.042`
  - `T_up = 0.9`: `tau ~ 6`, `mass ~ 0.032`
  - `T_up = 1.0`: `tau ~ 8`, `mass ~ 0.030`
  - `T_up = 1.1`: `tau ~ 6`, `mass ~ 0.030`
- These are extrapolated recommendations, not directly validated scan points, because all of the required masses are below the current scan minimum `0.1`.
- The physical target source constrains `298-313 K` directly; `0.7`, `0.8`, `1.0`, and `1.1` therefore rely on Arrhenius extrapolation outside the source measurement window.

## 2026-04-08 (Thermostat Einstein-Relation Bug)
- The current MARTINI runtime is internally inconsistent:
  - `src/deriv_engine.cpp` updates positions with the loaded per-atom masses from `/input/mass`,
  - but `src/thermostat.cpp` previously applied the same OU noise amplitude to every atom, as if every atom had unit mass.
- For an OU thermostat acting on momentum, the correct stationary distribution is `p ~ N(0, mT)`, so the exact noise term must scale as `sqrt(mass) * sqrt(T * (1 - exp(-2 dt / tau)))`.
- Without that `sqrt(mass)` factor, the code enforces `⟨p^2⟩ = T` instead of `mT`, which implies a free-particle diffusion scaling `D ~ T tau / m^2` rather than the Einstein-consistent `D ~ T tau / m`.
- That matches the scan symptom that the fitted diffusion dependence on mass was near `m^-2`, which is a strong sign the thermostat bug was contaminating the calibration.
- The weak / non-monotonic damping dependence in the bilayer scan is therefore not safe to interpret physically until the scan is rerun with the corrected thermostat.
- Kinetic-energy logging also needed a mass-aware correction because the old logger reported `0.5 <p^2>` rather than `0.5 <p^2 / m>` when masses were active.

## 2026-04-09 (Post-Fix Rerun Findings)
- The refreshed post-fix `stage7-analysis/` tree is internally clean:
  - `525` stage-7 files,
  - `525` successful tasks,
  - `175` condition rows,
  - no task with `fit_r2_mean < 0.95`.
- After fixing the thermostat, the previous strong mass trend disappeared:
  - stable-window temperature effect is about `2.87x`,
  - stable-window mass effect is about `1.36x`,
  - stable-window damping effect is about `1.05x`.
- That means the old `mass ~ m^-2` behavior was indeed a thermostat artifact rather than a trustworthy physical calibration signal.
- The `tau = 64` branch still shows unusually high replicate variance and should not be used as the main calibration anchor even after the thermostat fix.
- In the corrected scan, the closest stable scanned points to the physical DOPC target all sit at `mass = 0.1`, but they only reach roughly `15%`, `9%`, `6%`, `4%`, and `3%` of target diffusion for `T_up = 0.7, 0.8, 0.9, 1.0, 1.1`.
- Therefore the current post-fix grid does not resolve a reliable bilayer mass correction. Any literal extrapolation below the tested range is too severe to save as a calibrated recommendation.

## 2026-04-09 (Next Rerun Strategy After Thermostat Consistency Audit)
- Once workflow `16` is aligned back to the standard example thermostat default, the cleanest next bilayer rerun is a mass-only scan at fixed `tau = 5.0`.
- The refreshed post-fix scan does not justify reopening damping as a main knob yet:
  - inside the stable window, damping changes diffusion only weakly,
  - the fastest stable point at every tested temperature already sat at the lowest tested mass `0.1`.
- Therefore the next useful probe is to push mass lower on a logarithmic grid while holding both the thermostat convention and the temperature grid fixed:
  - `damping = 5.0`
  - `temperature = 0.70, 0.80, 0.90, 1.00, 1.10`
  - `mass_scale = 0.1, 0.05, 0.02, 0.01, 0.005, 0.002`
  - `replicates = 3`
- If that focused rerun still leaves the bilayer far too slow, then mass-only retiming is not enough and the remaining mismatch is not a simple thermostat/mass issue.

## 2026-04-10 (Preset-Aware Stage-7 Analysis Wrapper)
- The focused `tau5_mass_only` rerun introduced a second canonical run directory under `bilayer-lateral-diffusion/runs/`.
- A stage-7 analysis wrapper that hardcodes `runs/default` is therefore no longer sufficient for preset-based workflows.
- The clean wrapper-level fix is to let `submit_stage7_analysis.sh` derive its default base directory from a named preset such as `BILAYER_DIFF_RUN_PRESET=tau5_mass_only`, while still letting `BILAYER_DIFF_BASE_DIR` override everything.
- A matching preset wrapper (`submit_tau5_mass_only_stage7_analysis.sh`) keeps the run and post-analysis entrypoints symmetric and avoids repeated manual `BASE_DIR` overrides.

## 2026-04-10 (Focused Tau-5 Mass-Only Rerun Findings)
- The refreshed focused rerun is internally clean:
  - `90` stage-7 files,
  - `90` task rows,
  - `30` condition rows,
  - no task below `R^2 < 0.95`.
- With damping fixed at `tau = 5`, lowering mass from `0.1` to `0.002` increases geometric-mean diffusion by only about `1.61x`, while temperature still changes it by about `3.20x`.
- The best robust point in the current grid is:
  - `T_up = 1.1`,
  - `tau = 5`,
  - `mass = 0.005`,
  - `D = 0.365 um^2/s`,
  - `CV = 0.060`.
- The very lowest tested mass `0.002` is already noisy at hotter temperatures:
  - `T_up = 0.9` has `CV ~ 0.79`,
  - `T_up = 1.1` has `CV ~ 0.36`.
- Even with the deeper mass scan, the nearest in-grid points remain far below the physical DOPC target at all temperatures, so the current rerun still does not resolve a reportable mass correction.

## 2026-04-10 (Recommended Joint Follow-Up Grid)
- The focused mass-only rerun says two things at once:
  - lower mass does help,
  - but by itself it is not enough.
- The earlier broad post-fix damping scan says:
  - moderate damping changes are weak,
  - very large `tau` branches become noisy and are poor calibration anchors.
- Therefore the next useful scan is not another full broad grid and not another mass-only grid. It is a focused joint scan:
  - keep temperature fixed to the same Upside temperature for each condition,
  - keep `tau = 5` in the grid as the aligned baseline,
  - test a moderate weaker-damping window `5, 8, 12, 16, 24`,
  - focus mass on `0.02, 0.01, 0.005, 0.002, 0.001`.
- This grid is targeted at the current uncertainty:
  - can lower mass and moderately weaker damping combine constructively,
  - and where does the stability/noise boundary appear before reaching the clearly noisy extreme-`tau` regime.

## 2026-04-10 (Refreshed Joint Mass+Damping Reanalysis)
- The refreshed downloaded joint scan contains the expected full grid:
  - `375` stage-7 files,
  - `375` task rows,
  - `125` condition rows,
  - damping values `5, 8, 12, 16, 24`,
  - mass scales `0.02, 0.01, 0.005, 0.002, 0.001`.
- The raw assembled dataset is not globally valid as a calibration surface:
  - `62` task rows have `NaN` diffusion outputs,
  - `28` conditions produce non-finite condition-level diffusion,
  - `46` conditions exceed `CV > 0.25`,
  - `1` task falls below `R^2 < 0.95`.
- The first-pass adjacent-mass filter was too weak because it only rejected the immediately offending point. That allowed a lighter-mass branch to remain “stable” even after a large jump had already revealed a new dynamical regime.
- Treating the first large adjacent-mass jump as the start of an unstable lighter-mass branch fixes that failure mode. With a `4x` jump threshold, the robust subset drops to `49 / 125` conditions and no longer promotes the obvious `tau = 24` blow-up branch as a calibration candidate.
- Under the corrected robust subset:
  - temperature changes diffusion by about `2.95x`,
  - damping changes diffusion by about `1.86x`,
  - mass changes diffusion by about `1.32x`.
- The trustworthy part of the joint scan is still too slow:
  - best robust point is `T = 1.1`, `tau = 5`, `mass = 0.01`, `D = 0.435 um^2/s`,
  - no robust condition reaches `1.0 um^2/s`.
- The full grid does reach target-like or larger diffusion, but only on excluded branches:
  - the onset appears in `tau = 16-24` with low mass,
  - those branches are characterized by NaNs, large replicate variance, or abrupt mass-to-mass jumps,
  - they are instability-frontier probes, not reusable calibration points.
- Therefore the current joint scan does not resolve a robust per-temperature mass/damping correction. The operational conclusion is still:
  - `tau = 5` remains the safest damping anchor,
  - hotter temperatures and lighter masses help modestly,
  - a trustworthy timescale match is not yet demonstrated.

## 2026-04-10 (Boundary-Focused Follow-Up Grid)
- The refreshed joint analysis narrows the useful search region:
  - the trustworthy region still centers on `tau = 5-8`,
  - the instability frontier starts somewhere before the old `tau = 24` branch and is already visible by `tau = 16`,
  - `mass = 0.001` is no longer worth scanning in the default grid because it repeatedly produced NaNs or pathological diffusion.
- The next efficient scan is therefore a boundary-focused grid rather than another wide factorial pass:
  - `tau = 5, 8, 10, 12, 14, 16`
  - `mass = 0.02, 0.015, 0.01, 0.0075, 0.005, 0.003, 0.002`
  - `temperature = 0.70, 0.80, 0.90, 1.00, 1.10`
  - `replicates = 3`
- This grid preserves the robust anchors, refines the onset region on the damping axis, and samples the only mass window that still looks scientifically useful.

## 2026-04-11 (Boundary-Scan Reanalysis)
- The refreshed boundary scan is structurally better than the previous broad joint scan:
  - `630` stage-7 files,
  - `630` successful tasks,
  - `210` condition rows,
  - only `16` task-level `NaN` diffusion values,
  - only `9` non-finite condition rows.
- The robust subset grows substantially on this grid:
  - `125 / 210` conditions pass the current filters,
  - compared with `49 / 125` in the previous broad joint scan.
- But the scientific conclusion does not improve:
  - robust-subset damping effect is only about `1.33x`,
  - robust-subset mass effect is only about `1.30x`,
  - best robust point is `T = 1.1`, `tau = 5`, `mass = 0.005`, `D = 0.394 um^2/s`,
  - no robust condition reaches `1.0 um^2/s`.
- The boundary-focused scan therefore strengthens a sharper conclusion:
  - `tau = 5` remains the best production damping anchor,
  - `tau = 8` is usable but not better,
  - `tau = 10-16` primarily maps the approach to instability rather than exposing a better robust calibration region.
- The instability frontier is now clearer:
  - it appears repeatedly at `tau = 12-16`,
  - especially with `mass = 0.003` or `0.002`,
  - the excluded fast branches show either abrupt mass-to-mass jumps, high CV, or non-finite outputs.
- Therefore the current evidence no longer supports “try a bit more damping” as the next main idea. The useful interpretation is now:
  - damping retuning is largely exhausted as a clean acceleration mechanism in this regime,
  - the timescale mismatch is not being solved by the currently tested mass range either.

## 2026-04-11 (Integrator/Time-Step Diagnostic Design)
- In this codebase the integrator choice is global per run, not per subsystem:
  - `main.cpp` selects one integrator mode for the full simulation,
  - the current bilayer workflow and workflow `16` both use `integrator = v`,
  - a dryMARTINI-only stability probe therefore has to be approximated as a bilayer-only workflow, not as “different integrators for protein and bilayer” inside one coupled run.
- `mv` is not the obvious next dryMARTINI stability lever here because the useful bilayer-force `integrator_level` split is not established in the present workflow surface.
- The clean diagnostic comparison is therefore:
  - `integrator = v` versus `nvtc`,
  - on the already-identified unstable boundary points,
  - while shrinking `dt` before introducing any force clipping.
- `nvtc` should default to `max-force = 0.0` for this diagnostic. Otherwise a successful run could be stabilized by clipping rather than by the stepper itself, which would blur the interpretation.
- For wrapper design, placing each `(integrator, dt)` pair in its own run directory is simpler than mixing them inside one manifest:
  - Slurm array cardinality stays local to one scientific grid,
  - stage-7 discovery remains unchanged,
  - post-analysis can reuse the existing `submit_stage7_analysis.sh` machinery with only a thin loop wrapper.

## Lessons
- When the user corrects a physical scaling assumption, revise the calibration strategy before coding against the old assumption.
- For MARTINI trajectory analysis in this repo, verify the actual HDF5 metadata on a real stage file before hard-coding selection logic.
- When a cluster traceback shows an import-time type-annotation failure, check for missing `from __future__ import annotations` before assuming the runtime logic is wrong.
- When a workflow is intended to run from a specific checkout on shared filesystems, avoid `Path(...).resolve()` for repo-root discovery unless canonicalization is explicitly desired.
- When the user provides a cluster-specific module dependency, propagate it into every execution surface that launches the binary: submit wrapper, local wrapper, and generated Slurm scripts.
- When the user trims old checkpoint stages after production is complete, treat post-analysis as a separate workflow with its own manifest and outputs rather than trying to reuse the original simulation task pipeline.
- When the user specifies that temperature must remain fixed, do not absorb a timescale mismatch by moving temperature; use the repository’s unit-conversion contract first, then solve the remaining mismatch with the intended dynamical knobs.
- When the user wants the protein thermostat preserved, analyze the fixed-thermostat slice first and only recommend mass changes unless the runtime itself proves a damping change is necessary.
- When the user questions a physical scaling law, inspect the exact numerical update rule in the runtime before trusting empirical scan fits; a missing mass factor can masquerade as a parameter trend.
- When a post-fix rerun qualitatively changes the trend, update every downstream summary artifact that encoded the old interpretation; do not leave scripts or markdown with stale recommendation text.
- When borrowing a physical timestep calibration from other examples, audit the full dynamical default set first, especially `--thermostat-timescale` and `--time-step`; matching only one knob is not enough.
- When a workflow already supports the needed scientific parameters in its Python CLI, prefer exposing them through the wrapper and a focused preset script instead of changing the global scientific defaults for every run.
- When an analysis script encodes a prior scan’s default parameter window, resolve requested parameter subsets against the values actually present in the current dataset instead of failing on a new focused rerun.
- When regenerating saved report text after a focused rerun, make the markdown describe the actual scanned parameter ranges rather than reusing wording from an earlier broader grid.
- When a focused rerun shows one knob is too weak by itself, the next preset should test a targeted combination around the current robust region instead of jumping back to a full coarse factorial scan.
- When one analysis script introduces a robustness filter, every downstream calibration script must use the same filter; otherwise the repo will mix “fast” points with “trustworthy” points and save contradictory recommendations.
- When the current scan already identifies a clear instability frontier, the next preset should bracket that frontier with finer spacing and should drop parameter values that are already conclusively pathological.
- When a boundary-focused scan confirms that the robust optimum stays on the baseline branch, treat that as evidence against the broader hypothesis; do not keep reopening the same knob range without a new mechanism.
- When verifying a new preset wrapper, read the live manifest schema before fabricating checkpoint paths; this workflow keys task directories by `code`, not by an assumed `task_name`.
