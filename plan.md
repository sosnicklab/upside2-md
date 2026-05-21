## 2026-05-21 1RKL/1AFO Production Energy Drift Debug

### Project Goal
- Debug why `example/16.MARTINI/run_sim_1rkl.sh` output logs in `example/16.MARTINI/1rkl.*.log` and `example/16.MARTINI/1afo.*.log` still show `protein_potential` and `total_potential` accumulating downward instead of fluctuating around equilibrium, even though the protein remains structurally stable.
- Resolve the root cause without disabling, scaling, or empirically retuning physical interactions.
- If interaction data must change, patch the existing generated `.h5` artifacts under `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/` rather than regenerating the full table workflow.
- Compare MARTINI production timestep semantics against original Upside example workflows so the protein/Upside part does not silently run on a different timescale.
- Remove redundant explicit Upside CLI timestep/integrator settings from the MARTINI workflow when master/default behavior already covers them.

### Architecture & Key Decisions
- Treat the drift as a runtime/integrator, restart, or force-accounting issue until component metrics prove a table defect.
- Preserve SC-env, BB-env, CGL-SC, CGL-target, CGL-CGL, generic dry-MARTINI, ion, and bonded interactions.
- Prefer a minimal source fix in the active C++/Python execution path over parameter changes.
- Validate on copied or existing generated artifacts from `martini_test_1rkl_hybrid` to avoid expensive `.h5` regeneration.
- Finding: the drift is dominated by `protein_potential`; `cg_lipid_pair` and `cg_lipid_sc` stay comparatively bounded in the reported logs.
- Finding: the commented Predescu three-stage coefficient normalization in `DerivEngine::integration_cycle(VecArray,float,float,IntegratorType)` matches Predescu et al. 2012. In that paper, the stage coefficients sum to `q=3` and one outer application advances `q*dt`; `dt` is the average time per force evaluation, not the outer-step transit time.
- Finding: the coefficient convention is acceptable for existing Upside workflows because their timestep choices and trained potential behavior are calibrated around the existing engine semantics. The hybrid MARTINI workflow is different: it introduces physical MARTINI masses, stiff CGL-CGLD orientation carriers, tabulated dry-MARTINI interface forces, and a rigid stage-6 to flexible production release.
- Decision: do not change `DerivEngine::integration_cycle`; keep shared Upside/master behavior unchanged.
- Finding: stage 6 keeps the protein rigid, so its internal hbond count remains `4.5`; production releases the protein and hbond count rises into the 20s while `protein_potential` lowers. The stage-7 handoff needs a physical full-Hamiltonian relaxation before MD.
- Finding: the C++ minimization wrapper always switches to `"minimization"` stage; for a production-stage file this minimizes a different Hamiltonian than production uses and can produce nonsensical objectives.
- Revised Decision: production-stage minimization must preserve the production stage. Preproduction files can still switch to `"minimization"` as before.
- User constraint: shared Upside workflows matching `/Users/yinhan/Documents/upside2-md-master` must keep the existing default C++ integrator and minimization behavior.
- Finding: original Upside workflows such as `example/02.ReplicaExchangeSimulation/run.py` and `example/08.MembraneSimulation/0.normal.run.py` do not pass `--time-step`, `--integrator`, or `--inner-step`. Master uses default `--integrator v`, `--time-step 0.009`, and schedules public time with `inner_step=3`.
- Finding: master does use `--time-step` through `py/run_upside.py` and uses `--integrator mv`/`--inner-step` only for explicit multi-step workflows such as curvature and multistep examples. Explicit `--integrator v` is redundant for MARTINI because `v` is already the default.
- Finding: current `src/main.cpp` had drifted from master by setting `inner_step=1` for `v`; that would change ordinary Upside workflow scheduling and must be restored.
- Revised Decision: restore master-compatible `v` public-time scheduling (`inner_step=3`) and set MARTINI production back to `PROD_TIME_STEP=0.002`. Stability should come from physical masses plus production-stage minimization, not a hidden threefold timescale slowdown.
- Revised Decision: keep default minimization stage switching unchanged; add opt-in `--minimize-preserve-stage` and use it only for the stage-7 MARTINI production handoff minimization.
- New Question: using a separate MARTINI integrator weakens claims inherited from original Upside. Test whether removing the mass-aware MARTINI runtime path and using the shared `v` integrator preserves stability after the production-Hamiltonian handoff minimization.
- Finding: deleting `/input/mass` on a copied stage-7 artifact is unstable with the shared `v` integrator: total potential falls to `-24239.58 E_up`, Rg grows to `14.9 A`, and kinetic ratio rises to `5.70`.
- Finding: scaling pair energies by mass is not physical because it changes the scalar potential and configurational distribution; a pair has two generally different masses, so per-particle mass-weighted pair forces would not be a single conservative pair potential.
- Finding: copied 1RKL timing validation with `PROD_TIME_STEP=0.002` and restored `inner_step=3` scheduling wrote `output/time` from `0.0` to `60.0`, frame spacing `0.3`, and `restart_public_time_step=0.006`; protein Rg remained about `12.9 A`.
- Final Decision: retain physical masses, use the shared default Verlet integrator without passing `--integrator v`, restore master-compatible public-time scheduling, keep MARTINI MD `--time-step` settings, and keep the production-Hamiltonian handoff minimization.
- Follow-up Result: removed redundant MARTINI `--integrator v`, removed minimization-only `--time-step`/`MIN_TIME_STEP`, and kept MD `--time-step` for the intentional MARTINI equilibration/production timesteps.

- [x] Phase 1: Inspect current logs, stage files, diffs, and exact reported output artifacts.
- [x] Phase 2: Quantify which energy components are drifting and whether drift correlates with restart boundaries, thermostat/barostat behavior, or integration step count.
- [x] Phase 3: Trace the responsible source path and identify the smallest physical correction.
- [x] Phase 4: Patch source or existing HDF5 artifacts as appropriate.
- [x] Phase 5: Run focused validation on the existing output artifacts and compare energy behavior before/after.
- [x] Phase 6: Test shared-integrator/unit-mass alternative and choose final workflow path.
- [x] Phase 7: Compare original Upside workflow timestep/integrator settings and correct MARTINI production settings if needed.
- [x] Phase 8: Remove redundant explicit CLI flags and update `cg_lipid_potentials.tex` with the final timestep convention.

### Current Investigation
- [x] Phase 9: Parse the current 1RKL and 1AFO logs and compare component trends, restart boundaries, hbond/Rg stability, and output timing.
- [x] Phase 10: Inspect the exact generated HDF5 outputs for per-node potential components to determine whether `protein_potential` is a real physical relaxation or a logging/accounting accumulation bug.
- [x] Phase 11: Implement the smallest physical correction. Do not weaken, scale, disable, or retune interactions.
- [x] Phase 12: Validate on copied exact artifacts and update documentation.

### Known Errors / Blockers
- Existing generated checkpoint outputs were produced with the old attractive CGL-ion target node.  The source fix applies when the workflow prepares/injects new stage files; patching old checkpoints in place would make their stored trajectories inconsistent with their input potential metadata.

### Review
- The Predescu coefficient pattern in the shared overload is consistent with the paper; no shared `DerivEngine` implementation change is made.
- MARTINI-specific changes:
  - Use shared default Verlet behavior; no separate MARTINI integrator remains and the workflow no longer passes redundant `--integrator v`.
  - Retain physical `/input/mass`; do not scale energies by mass.
  - Restore master-compatible `v` public-time scheduling (`inner_step=3` for time accounting).
  - Keep MARTINI production input timestep at `0.002`; frame/restart public time uses `3*dt`.
  - Keep MD `--time-step` explicit for MARTINI equilibration/production, but omit `--inner-step` and omit minimization-only `--time-step`.
  - Add `--minimize-preserve-stage` and use it for stage-7 MARTINI production handoff minimization.
  - `example/16.MARTINI/run_sim_1rkl.sh` exposes `MIN_70_MAX_ITER`, default `500`.
- Copied-artifact validation:
  - Original stage-7.0 log: total potential `164.6 -> -1068.85 E_up`, kinetic ratio `1.222`.
  - Unit-mass test failed: total potential `-505.9 -> -24239.58 E_up`, kinetic ratio `5.70`.
  - Shared `v` with physical masses and `dt=0.002`: copied timing validation wrote final output time `60.0`, frame spacing `0.3`, restart public timestep `0.006`, kinetic ratio `0.972`, and protein Rg about `12.9 A`.
- Latest root cause:
  - The old `protein_potential` log bucket included `cg_lipid_target`; actual Upside protein internal energy is not the multi-thousand-`E_up` sink.
  - The real `total_potential` drift is dominated by `cg_lipid_target`, especially mobile ions.  In 1RKL stage 7.3, split diagnostics gave about `-1558 E_up` from BB targets and `-3401 E_up` from ion targets; another copied continuation drove the ion part to about `-4170 E_up`.
  - The fix splits mobile ions out of the CGL-target node and derives an excluded-volume-only ion table from the current `cg_lipid_target` table by clipping negative controls to zero, while leaving BB/protein targets unchanged.
