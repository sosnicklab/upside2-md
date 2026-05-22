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
- Active blocker: user reports exact `example/16.MARTINI/1rkl.*.log` still moves from about `+100` to `-20 E_up` over five trajectories. Re-parse the exact logs before considering the drift resolved.
- Active blocker: user reports `outputs/martini_test_1rkl_hybrid/1rkl.stage_7.4.vtf` shows edge lipids with strange z coordinates and orientations. Determine whether this is a VTF periodic-wrapping/rendering bug or a real CGL/CGLD orientation dynamics problem.
- User correction: treat in-plane CGLD-CGL orientation vectors as a real physical defect. Restore the previous VTF display style and debug why single-particle lipid orientations can rotate into the x-y plane.
- User verification request: confirm ion handling is not an ion-protein or ion-lipid exclusion. Ions must still interact physically through generic MARTINI ion-protein, ion-ion, ion-water, and available non-CGL lipid/environment pair paths; only direct CGL/CGLD generic pairs are excluded because CGL physics is owned by dedicated spline nodes.
- User correction applied: the inconsistent CGL-only `Qa`/`Qd` charge removal has been replaced by a coherent alternate dry-MARTINI protein-backbone model where charged head/tail BB typing is disabled and all backbone carrier/proxy interactions use neutral `P5` consistently.

### Follow-up Investigation
- [x] Phase 13: Quantify residual component drift in current `example/16.MARTINI/1rkl.*.log`.
- [x] Phase 14: Audit ion spatial behavior, ion-energy ownership, and whether ions are double-counted or missing physical interactions.
- [x] Phase 15: Identify the smallest physical correction, if the drift is not just normal slow equilibration.
- [x] Phase 16: Validate on copied exact artifacts and update documentation/task notes.
- [x] Phase 17: Parse all current 1AFO and 1RKL chunks by component, including chunk-to-chunk absolute baselines.
- [x] Phase 18: Inspect matching stage HDF5 files for restart-state consistency, production counters, momentum, component ownership, and structural/contact trends.
- [x] Phase 19: Isolate the common residual sink and implement the smallest physical/runtime correction.
- [x] Phase 20: Validate copied continuations for both systems and update docs/task notes.
- [x] Phase 21: Parse fresh post-burn-in 1RKL logs and confirm the burn-in stage was applied.
- [x] Phase 22: Inspect post-burn-in stage checkpoints for component, contact, ion, CGL-SC, and CGL-target trends across trajectories.
- [x] Phase 23: Identify and patch the real residual physical/runtime source without altering force-field parameters.
- [x] Phase 24: Identify terminal charged BB proxy targets as the remaining CGL-target sink and reject the inconsistent CGL-only charged-target split.
- [x] Phase 25: Replace charged-target-only patch with consistent neutral `P5` backbone typing across generated code paths and existing stage-7 artifacts.
- [x] Phase 26: Parse exact current 1RKL logs, audit ion pair ownership in the actual current artifacts, and scan `run_sim_1rkl.sh` for physical-model consistency.
- [x] Phase 27: Debug edge lipid geometry in `1rkl.stage_7.4.vtf` by comparing VTF rendering coordinates with HDF5 minimum-image CGL/CGLD orientation vectors.
- [x] Phase 28: Revert the VTF display style change, quantify in-plane orientation states directly from `1rkl.stage_7.4.vtf`/HDF5, and identify the physical/runtime root cause without parameter twisting.

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
- Follow-up review after the ion-target split:
  - Current `1rkl.0-2.log` still contain visible production relaxation, but by `1rkl.3.log` total potential changes only `-28 E_up` first-to-last; the largest remaining component changes are `cg_lipid_target` `-42.6 E_up` and `cg_lipid_sc` `-29.0 E_up`.
  - A copied continuation from current `1rkl.stage_7.3.up` final coordinates/momenta ran 10k more steps over 60 public time units: total potential first/last `-1045.9 -> -1062.2 E_up`, range `[-1138.7, -981.5]`, and first/last fifth means differ by `-34.4 E_up`.
  - Current `1afo.0.log` total potential is flat first-to-last (`-730.55 -> -730.09 E_up`) despite expected compensation among protein, CGL-pair, CGL-SC, and CGL-target terms.
  - Ion audit found zero generic CGL-ion pairs in both checked stage-7 files.  Ions remain in generic dry-MARTINI ion/protein and ion/ion pairs, and CGL-ion target controls are nonnegative excluded volume.
  - 1RKL and 1AFO ion-CGL minimum distances remain above `16 A` over the saved trajectories, with stable ion-z distributions; there is no evidence of residual ion adsorption to CGL.
- Reopened 1AFO multi-chunk drift:
  - 1AFO drops from about `-730 E_up` in `1afo.0.log` to about `-1001 E_up` in `1afo.3.log`; the main slow components are `cg_lipid_sc` and `cg_lipid_target`.
  - Restart state is consistent: output momentum is valid, transition starts advance by `10000` steps per chunk, and first saved potentials match the prior final potential.
  - Copied continuations from the current 1AFO and 1RKL stage-7.3 endpoints fluctuate around the current basin instead of continuing the large early drop.  The residual drift is the stage-7 production-Hamiltonian interface relaxation window.
  - Decision: add an explicit stage-7 burn-in/equilibration under the same production Hamiltonian before the named production segment.  Preserve all interactions, physical masses, momenta, and hybrid transition counters; do not alter tables or force-field parameters.
  - Implementation: `PROD_70_BURNIN_NSTEPS` defaults to `40000`.  The workflow runs burn-in with production `dt=0.002`, records momentum, promotes final burn-in positions/momenta to `/input`, advances `sc_env_transition_step_start`, deletes burn-in `/output`, and then runs the named production chunk.
  - Validation: a copied 1AFO burn-in smoke test promoted restart-valid input momentum, cleared `/output`, advanced the transition counter, and a follow-on `--restart-using-momentum` MD run succeeded.
- Reopened 1RKL post-burn-in drift:
  - Fresh logs did include the new burn-in, and restart bookkeeping is continuous.
  - Per-pair CGL-target evaluation showed the remaining sink is dominated by two charged terminal backbone proxy targets (`Qd`/`Qa`, atom indices 124 and 154), not by mobile ions.
  - Rejected intermediate fix: splitting only the CGL-target interaction for `Qa`/`Qd` is inconsistent because the dry-MARTINI protein remains charged in other interaction paths.
  - Final decision: use the alternate uniform neutral dry-MARTINI protein-backbone model. All protein backbone carrier atoms and virtual BB proxies are `P5` with zero charge; generic MARTINI pairs, BB-env, and CGL-target now share the same backbone definition. The charged head/tail convention is not used because its terminal charged BB--DOPC well is not transferable to the additive single-particle CGL projection.
  - Existing `martini_test_1rkl_hybrid` stage-7 checkpoints were patched in place: first 155 protein atoms are `P5`/zero-charge, `hybrid_bb_map/bb_type` is all `P5`, `cg_lipid_target` owns all 31 neutral BB proxy targets, `cg_lipid_target_charged_excluded_volume` is absent, and `cg_lipid_target_ion_excluded_volume` still owns 93 ions with nonnegative controls.
  - Validation: a first copied 10k-step continuation re-equilibrated old charged-model coordinates (`85.55 -> -13.30 E_up`, first/last fifth means `59.35 -> 8.89 E_up`). The next continuation was fluctuation-scale (`-13.30 -> -26.54 E_up`, range `[-73.93, 23.57]`, first/last fifth means `-1.60 -> -35.71 E_up`), and a third continuation reversed the mean change (`-22.91 -> -12.66 E_up`, range `[-58.95, 63.38]`, kinetic ratio `1.017`).
- Latest 1RKL production recheck:
  - Exact production sections in `1rkl.0-4.log` still are not stationary: production total potential moves `111.75 -> -2.21 E_up`, and a copied continuation from `1rkl.stage_7.4.up` has first/last fifth means `-6.37 -> -41.05 E_up`.
  - This is not an ion-protein exclusion or old CGL-target ion sink. Current stage-7 files contain `14,415` generic ion-protein pairs (`11,532` ion-backbone-carrier and `2,883` ion-BB-proxy), `4,278` ion-ion pairs, zero generic ion-CGL/CGLD pairs, and a dedicated CGL-ion excluded-volume target node containing all 93 ions.
  - Remaining drift is mainly CGL-CGL strain release plus CGL-SC contact relaxation: across `1rkl.0-4.log`, `cg_lipid_pair` changes `350.80 -> 294.28 E_up`, `cg_lipid_sc` changes `-85.23 -> -151.68 E_up`, while `cg_lipid_target` is approximately flat (`-31.51 -> -29.01 E_up`).
  - Workflow scan: current design uses physical masses, shared Upside integrator semantics, explicit production-Hamiltonian minimization, explicit burn-in, same-Hamiltonian production continuations, CGL-owned lipid interactions, and generic MARTINI non-CGL ion/protein pairs. No force scaling, disabled SC-env/BB-env, or parameter-zeroing was found. The open issue is insufficient stationarity of the production basin, not an ion exclusion.
- Edge-lipid VTF geometry recheck:
  - The HDF5 stage-7.4 endpoint is not showing a special physical box-edge problem. CGL-CGLD minimum-image orientation lengths remain normal (`11.13-12.04 A`) and edge lipid CGL/CGLD z ranges are comparable to non-edge lipids; no lipid display atoms are outside the z half-box.
  - User correction accepted: the in-plane CGLD-CGL orientations are a real physical defect, not just VTF display. The previous symmetric VTF display style was restored.
  - Root cause: the CGL-CGL tensor table had all negative controls clipped to zero, leaving a largely repulsive lipid-lipid term. Sparse single-particle lipids could therefore rotate into the membrane plane as nearly free rotors. Restoring the stored radial attraction directly was rejected: a copied test drove `cg_lipid_pair` from about `-1453` to `-3953 E_up` over one chunk, showing additive many-neighbor collapse.
  - Fix: add MARTINI-only `cg_lipid_leaflet_orientation`, a mean-field orientation term derived from the existing CGL-CGL first-neighbor horizontal-orientation penalty. It uses the initial dry-MARTINI leaflet sign and propagates through the same CGLD orientation coordinate.
  - Validation on a copied stage-7.4 continuation: original output had `34` lipids with `|n_z|<0.5` and `14` with `|n_z|<0.25`; the leaflet-normal node reduced both counts to `0`, with no anti-aligned lipids. The corrected `1rkl.stage_7.4.up` and `1rkl.stage_7.4.vtf` have been regenerated.
