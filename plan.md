
## 2026-05-20 1RKL Residual Protein Instability

### Project Goal
- Debug the newly reported residual instability in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf`.
- Determine whether the remaining failure is lipid collapse, protein-lipid over-attraction, sidechain/environment routing, backbone environment routing, stage minimization/handoff, or stale generated tables.
- Keep SC-env, BB-env, CGL-SC, CGL-target, CGL-CGL, generic dry-MARTINI, ion, and bonded interactions active and physical.
- Validate candidate fixes by modifying copied `.h5/.up` artifacts from the reported output before considering a full workflow rerun.

### Architecture & Key Decisions
- Treat the current named output as authoritative, even if previous direct-HDF5 validation looked improved.
- Reopened Decision: the user reports `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` still looks unchanged and loses secondary structure. The exact named output must be audited before accepting any validation from a separate directory.
- Reopened Decision: protein stability must be evaluated with explicit secondary-structure proxies, including hbond retention and CA alignment/RMSD, not only lipid spacing and protein Rg.
- Compare scalar logs with coordinate-derived metrics: CA RMSD/Rg, hbond sum, CGL-CGL nearest distances, CGL-protein distances, CGL-SC contacts, lipid orientation, and table attrs.
- If the protein fails while lipid spacing remains physical, focus on CGL-SC/CGL-target/protein-env physical table semantics or force routing, not CGL-CGL.
- If the failure is from stale tables or missing stage-6 minimization, fix rebuild/guard/workflow behavior rather than changing physical interactions.
- Revised Decision: stage minimization must promote its final minimized coordinates and box back into `/input` and clear restart momentum before the following MD stage. Running a minimizer that leaves its result only under `/output` is a no-op for the subsequent stage.
- Revised Decision: the minimization promotion fix is necessary but not sufficient. Direct copied-HDF5 validation improves CGL-protein contact geometry, but CA RMSD still reaches about `5.2 Å` and CGL-protein minimum distance reaches `2.6 Å`, so the remaining failure path is protein-lipid interface semantics rather than only stage handoff.
- Revised Decision: DOPC CG table construction must use the actual lipid molecule charges from `dry_martini_v2.1_lipids.itp`, not bead-type inferred charges. `NC3` has type `Q0` but charge `+1`; treating `Q0` as neutral makes projected DOPC net negative and creates nonphysical attraction to charged BB/ion targets.
- Revised Decision: the remaining CGL-target failure is spline overshoot, not an interaction to disable. Charged target rows (`Qd/Qa`) in the fitted directional tensor contain deep negative controls while the short-range physical samples are repulsive; constrain only the dry-MARTINI repulsive rows to nonnegative controls, matching the existing excluded-area treatment for CGL-CGL.
- Revised Decision: the current exact output uses fresh charge/target tables, but CGL-SC sub-contact controls remain attractive. Because a CGL-SC separation inside the DOPC-derived contact shell represents unresolved sidechain overlap with the hidden DOPC molecule, those contact-region controls must be non-attractive while resolved outside-contact CGL-SC attraction remains active.

### Execution Phases
- [x] Phase 1: Re-audit the exact reported 1RKL output, generated `martini.h5`, checkpoint HDF5, logs, timestamps, and VTF-equivalent geometry.
- [x] Phase 2: Determine whether the exact output is stale or whether the current code still fails when applied in-place to copied `.h5` artifacts.
- [x] Phase 3: If current code still fails, identify the remaining physical interaction path driving secondary-structure loss.
- [x] Phase 4: Implement the minimal physical fix without disabling, scaling, excluding, capping, or restraining required interactions.
- [x] Phase 5: Quick-validate by direct copied-HDF5 modification from the reported output, including hbond/CA RMSD checks.
- [x] Phase 6: Run focused regression checks and update model documentation/logs.

### Known Errors / Blockers
- The reported exact `martini_test_1rkl_hybrid/martini.h5` is now intentionally rejected by the stale guard because it lacks CGL-SC DOPC-contact controls.
- The full workflow was not rerun; validation used direct copied-HDF5 modification in `martini_test_1rkl_hybrid_sc_contactcheck` as requested.

### Review
- Root causes:
  - Stage minimization did not promote minimized coordinates and box back into `/input`, so subsequent MD started from the unrelaxed input.
  - DOPC CG table construction inferred bead charges from bead type, missing the lipid-specific NC3 `+1` charge.
  - CGL-target fitted spline controls had attractive short-range overshoot in the DOPC excluded-area domain for charged target rows.
  - Fresh exact-output validation exposed a remaining CGL-SC table error: spline controls inside the DOPC contact shell stayed attractive, making sub-contact sidechain/lipid overlap an artificial energy sink.
- Implemented physical fixes:
  - Minimized stages now promote final `/output/pos` and `/output/box` to `/input`, clear restart momentum, and keep minimization logs separate.
  - DOPC CG table fitting uses explicit DOPC molecule charges and records `bead_charge_source` plus `lipid_net_charge`.
  - CGL-target rows whose radial support contributes inside DOPC contact are constrained to nonnegative controls and guarded by stale-table validation.
  - CGL-SC rows whose radial support lies inside the DOPC contact shell now have nonnegative isotropic controls and zero angular-mode amplitudes, so unresolved overlap is excluded while outside-contact CGL-SC attraction remains active.
- Direct copied-HDF5 1RKL validation:
  - Using the promoted stage-7 input and patched target table, final CA RMSD improved from `5.203 Å` to `3.902 Å`; Rg stayed `12.861 Å` instead of compacting to `11.866 Å`.
  - Final CGL-protein min/p05 improved from `2.620/9.248 Å` to `3.722/8.385 Å`.
  - Final CGL-CGL nearest-neighbor min/p05 was `5.676/5.884 Å`; lipid orientation p05 remained `0.796` with one transient near-parallel outlier.
  - CGL-SC contact-control validation in `martini_test_1rkl_hybrid_sc_contactcheck` improved final hbond retention from `1.22` to `17.96`, final CA RMSD from `5.86 Å` to `3.66 Å`, and final CGL-SC energy from `-1319 E_up` to `-73.6 E_up`; CGL-CGL final min/p05 stayed physical at `5.781/5.982 Å`.
- Verification:
  - `python3 -m py_compile py/martini_cg_lipid_params.py py/martini_build_tables.py py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed.
  - Stale guard rejects the reported exact `martini.h5` for missing CGL-SC DOPC-contact controls and accepts the patched contact-control validation `martini.h5`.
  - Generated `example/16.MARTINI/outputs/martini_test_1rkl_hybrid_sc_contactcheck/1rkl.stage_7.0.vtf` for VMD inspection.
  - `git diff --check` passed.

## 2026-05-19 1RKL/1AFO Protein Stability Regression

### Project Goal
- Debug the new unstable protein outputs in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid` and `example/16.MARTINI/outputs/martini_test_1afo_hybrid`.
- Identify whether the instability is from protein core terms, protein-env interfaces, CGL-target/CGL-SC/CGL-CGL tables, restart/handoff state, or table/schema mismatch.
- Keep all protein-lipid, protein-env, lipid-lipid, SC-env, and BB-env interactions physical and active; do not zero, exclude, scale, restrain, or cap interactions to hide the issue.
- Update `example/16.MARTINI/cg_lipid_potentials.tex` if the physical model changes.

### Architecture & Key Decisions
- Start from the named outputs, not assumptions from bilayer-only validation.
- Compare both proteins to find common failure signatures: timing, component energies, protein geometry, lipid morphology, and injected table schemas.
- Treat existing output `martini.h5` files as suspect until their attrs match the current schema guard.
- If the failure is from a stale table, fix workflow reuse/guard behavior rather than altering physics.
- If the failure is from a force-routing or handoff bug, patch the implementation while preserving full interface inclusion.
- Validate the debug on 1RKL by directly patching copied HDF5 artifacts from `martini_test_1rkl_hybrid`, avoiding a full workflow rerun.
- Revised Decision: CGL-CGL sampled rows keep charged dry-MARTINI bead relaxation, but unresolved inner spline rows must not inherit a weak boundary row. Fill only those unresolved rows from the maximum first sampled dry-MARTINI energy expectation and keep the DOPC-contact WCA excluded-area projection; this is a force-field-derived finite core, not a chosen cap.
- Revised Decision: run a physical stage-6.0 minimization before rigid-protein NPT dynamics by default. The packed bilayer starts with finite residual stress; relaxing it under the full active protein-env, protein-lipid, lipid-lipid, ion, and bonded interactions prevents the hot handoff that drives later CGL-SC attraction. This is an integration/handoff correction, not an interaction exclusion, scale, restraint, or empirical cap.

### Execution Phases
- [x] Phase 1: Audit output files, logs, stage HDF5 schemas, and table attrs for both proteins.
- [x] Phase 2: Quantify protein and bilayer geometry over production frames.
- [x] Phase 3: Identify the shared unstable interaction path or stale-artifact mismatch.
- [x] Phase 4: Implement the minimal physical fix.
- [x] Phase 5: Validate with focused reruns or direct HDF5 tests, plus bilayer-only regression if lipid terms change.
- [x] Phase 6: Update docs and logs.

### Known Errors / Blockers
- Existing `martini_test_1rkl_hybrid/martini.h5` and `martini_test_1afo_hybrid/martini.h5` are stale and are intentionally rejected by the schema guard. Rebuild them to get the new CGL-CGL unresolved-core attrs.
- Full workflows can be expensive; direct copied-HDF5 validation was used for 1RKL as requested.

### Review
- Root cause: both reported outputs show lipid collapse before protein failure. In 1RKL stage 7, CGL-CGL nearest-neighbor p05 falls from `4.313 Å` to `0.215 Å`, CGL-protein minimum reaches `0.267 Å`, kinetic energy reaches `67.6`, and CA RMSD reaches `8.81 Å`.
- Implemented physical fixes:
  - CGL-CGL sampled rows still use charged dry-MARTINI relaxation, but unresolved inner spline rows are floored from the maximum first sampled dry-MARTINI energy expectation rather than a weak boundary row or chosen cap.
  - Stage 6.0 now runs a default `500`-iteration minimization under the full active Hamiltonian before rigid-protein NPT MD.
  - The stale-table guard rejects old `martini.h5` files with `unresolved_core_source=first_resolved_dry_martini_energy_expectation`.
- Direct 1RKL copied-HDF5 validation:
  - Source-equivalent CGL-CGL table plus minimized stage-6 handoff kept final stage-7 CGL-CGL nearest-neighbor min/p05 at `5.968/6.200 Å` versus bad `0.067/0.215 Å`.
  - Final CGL-protein min/p05 was `5.197/9.729 Å` versus bad `0.267/1.328 Å`.
  - Final CA RMSD was `4.23 Å` and Rg `12.61 Å` versus bad CA RMSD `8.81 Å` and Rg `8.94 Å`.
- Bilayer-only source validation passed:
  - `bash example/16.MARTINI/test_cg_bilayer/run_test.sh --resolution coarse --stage60-npt-contract --steps 1000 --dt 0.01 --frame-steps 50`;
  - final aligned-z min/p05/mean `0.664/0.838/0.936`, `bad_parallel=0`, `bad_flip=0`, leaflet crossings `0/0`;
  - same-leaflet nearest-neighbor min/p05 `4.708/5.650 Å`, CGL-CGLD length RMSD `0.261 Å`.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system.py py/martini_prepare_system_lib.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed.
  - `git diff --check` passed.

## 2026-05-17 1RKL Protein/Bilayer Stability Regression

### Project Goal
- Debug and resolve the reported `1rkl.stage_7.0.vtf` protein unfolding after the latest CG-lipid orientation model changes.
- Determine separately whether the bilayer is stable and whether protein unfolding is driven by CGL-CGL, CGL-SC, CGL-target, generic Martini pairs, or stage handoff/setup.
- Preserve physical hybrid interface interactions; do not disable SC-env, BB-env, CGL-SC, or CGL-target terms, and do not use empirical scales, orientation pins, or restraints to hide the issue.
- Keep `example/16.MARTINI/cg_lipid_potentials.tex` synchronized if the model changes.

### Architecture & Key Decisions
- Diagnose from the current `example/16.MARTINI/outputs/martini_test_1rkl_hybrid` artifacts before changing model code.
- Compare protein geometry, lipid morphology, and component potentials over `stage_6.0` and `stage_7.0`.
- If a CGL-target or CGL-SC term is responsible, fix the physical table/injection/evaluator semantics rather than turning the interaction off.
- Re-run focused 1RKL validation after any fix, and re-run bilayer-only validation if CGL-CGL or CGLD mechanics are touched.
- Revised Decision: the full workflow is too slow for this debug cycle; validate by rebuilding only the affected CGL-CGL spline table from existing 1RKL packed coordinates, patching copied HDF5 artifacts, and running short MD from those patched HDF5 files.
- Revised Decision: relaxation during explicit-bead CGL table fitting is physical only if it charges bonded deformation relative to the unrelaxed bead reference; free bead relaxation makes overlaps artificially cheap and destabilizes the bilayer.
- Revised Decision: after NPT handoff, every PBC-aware potential node must receive the handed-off box, not only `martini_potential`.
- Revised Decision: standalone CGL-CGL bilayers are stable at the stage-6 timestep, so the remaining 1RKL collapse must be addressed in CGL-target table construction rather than by weakening protein/interface interactions.
- Revised Decision: CGL-target fitting needs the same unresolved-overlap treatment as CGL-CGL and CGL-SC, but not by numeric caps or pairwise PMF; use dry-MARTINI energy-expectation sampling and fill unresolved radial knots from the first resolved energy-expectation boundary.
- Revised Decision: CG-lipid potential tables should be fit from the canonical dry-MARTINI DOPC reference geometry; packed system coordinates remain placement/input geometry and should not change the force-field table coefficients.
- Revised Decision: preproduction/minimization to production handoff should transfer coordinates and box only, then rethermalize velocities at the production temperature. Momentum reuse is reserved for exact production-to-production continuations with matching restart semantics.
- Revised Decision: hidden `CGLD` orientation sites are lipid internal coordinates, not protein. Hybrid membership remapping must initialize unmapped runtime atoms as non-protein and assert that both `CGL` and `CGLD` have negative protein membership.
- Revised Decision: do not exclude `BB` or any other protein/environment target from `cg_lipid_target`. The last working protein-env behavior depended on complete interface inclusion; any instability must be fixed by correcting physical table semantics, force routing, or handoff state, not by removing interaction paths.
- Revised Decision: `cg_lipid_target` forces on active virtual `BB` sites must use the same BB-proxy-to-carrier gradient projection as `martini_potential`; otherwise BB-target interactions are present in energy but their protein-side gradient is lost or misapplied.
- Revised Decision: backbone carrier atoms (`N/CA/C/O`) are not independent dry-MARTINI targets. `cg_lipid_target` must include the physical `BB` proxy site and project its force to carriers, but it must not also target each carrier atom as a separate CGL interaction site.
- Revised Decision: the unresolved CGL-CGL excluded-volume core must come from the underlying dry-MARTINI energy expectation. The prior fixed `5000 kJ/mol` core cap was stable but not acceptable model physics.
- Revised Decision: generated stages must reject stale `martini.h5` tables that lack the CGL-CGL energy-expectation unresolved-core attrs, because silently reusing the latest bad output table would reproduce the protein/bilayer collapse.
- Revised Decision: arbitrary energy caps are not acceptable model physics. Replace the fixed CGL-CGL hard cap and sampled-energy clipping with dry-MARTINI-derived energy-expectation table values; unresolved radial knots must be filled from the first resolved energy-expectation boundary rather than a chosen numeric cap.
- Revised Decision: the dry-MARTINI pair PMF over unresolved lipid azimuths caused many-body lipid aggregation in the latest output. Treat it as a rejected model for additive pair potentials; re-evaluate using force-field energy expectation plus physically derived unresolved-core boundary rather than pairwise PMF.
- Revised Decision: verify SC-env and BB-env interactions are present and nonzero in the current generated output before relying on protein stability as evidence.
- Revised Decision: the first energy-expectation boundary-row fill is also insufficient for CGL-CGL, because a flat attractive unresolved core lets same-leaflet CGL centers collapse below dry-MARTINI contact in bilayer-only validation. The unresolved core must be sampled directly at shorter CGL separations from the explicit DOPC force field with charged relaxation, not filled by a constant boundary row.
- Revised Decision: for CGL-CGL separations below the DOPC-derived contact distance, unresolved azimuths should use the dry-MARTINI steric envelope rather than a plain energy expectation. This is a projection of the excluded volume occupied by hidden bead azimuths, not an empirical cap.
- Revised Decision: the CGL-CGL tensor fit must enforce nonnegative spline controls on rows that contribute to sub-contact separations. Negative sub-contact spline coefficients are fit overshoot, not dry-MARTINI physics, because the steric-envelope samples represent excluded volume there.
- Revised Decision: the steric-envelope fit still allowed rod overlap in bilayer-only validation, so the simpler physical CGL-CGL projection is rigid canonical DOPC bead sampling. The runtime CG lipid has only center/vector DOFs, so relaxing hidden beads during CGL-CGL table fitting adds an internal deformation channel that the simulation cannot represent.
- Revised Decision: rigid canonical DOPC sampling creates unresolved-core energies too stiff for the current spline/timestep. The accepted finite-core representation must be derived from the maximum relaxed dry-MARTINI sampled energy below the DOPC contact distance, recorded in `martini.h5`, instead of a hand-picked cap.
- Revised Decision: a flat maximum-energy core across all rows inside DOPC contact is also rejected. It is force-field-derived but physically overbroad for the spline support: it repels normal contact configurations, injects heat, and destroys the bilayer. Preserve the sampled dry-MARTINI radial/angular shape and only use force-field-derived short-range values where the spline is genuinely unresolved.
- Revised Decision: CGL-SC fixed fitting caps are not acceptable for final model documentation. Remove fixed sampled-energy/residual clipping from CGL-SC and record sidechain short-range rows from sampled dry-MARTINI energy expectations when SC rows exist.
- Revised Decision: audit all protein-lipid and lipid-lipid terms before finalizing: classify each as explicit dry-MARTINI-derived physics, standard Upside protein physics, or documented numerical table-fitting regularization.
- Revised Decision: a divergent WCA excluded-area projection sampled down to `1.4 Å` is numerically invalid for the spline table and was rejected. Use the DOPC-contact WCA only over the physically resolved compressed CGL-center range starting at `5.0 Å`.
- Revised Decision: reduced CGL-CGL angular grids are not acceptable because they admit same-leaflet overlap. CGL-CGL fitting always uses the full grid (`16r x 7^2 theta x 2^2 phi`); reduced resolution may only apply to sidechain table construction.
- Revised Decision: within the DOPC-derived excluded-area domain, negative CGL-CGL spline controls are fitting overshoot, not physical attraction. Clamp only those excluded-area rows to nonnegative values and record the row count in `martini.h5`.

### Execution Phases
- [x] Phase 1: Audit current 1RKL output geometry, energies, and injected node schemas.
- [x] Phase 2: Identify the unstable interaction path and root cause.
- [x] Phase 3: Patch the minimal physical model or injection bug.
- [x] Phase 4: Re-audit the freshly reported `martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` instability.
- [x] Phase 5: Patch the root cause and revalidate without disabling or scaling physical interactions.
- [x] Phase 6: Document final results and update model notes if needed.

### Known Errors / Blockers
- Existing `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/martini.h5` is stale and is intentionally rejected for reuse; a fresh workflow must rebuild it to get the current CGL-CGL/CGL-target table schema.
- Full 1RKL workflow was not rerun in this pass because the user requested bilayer-only validation or direct HDF5 modification where possible.

### Review
- Root cause: the reported protein unfolding followed lipid collapse. The current bad artifact reaches same-leaflet CGL XY nearest distance near `0.91 Å` by production frame 50 and then the protein Rg grows to `25.9 Å`.
- Fixed table construction without disabling physics:
  - explicit-bead lipid relaxation now charges dry-MARTINI bonded deformation relative to the unrelaxed reference;
  - CGL-target tables now use dry-MARTINI energy-expectation samples and fill unresolved radial knots from the first resolved energy-expectation boundary;
  - CG-lipid potential tables use the canonical dry-MARTINI `DOPC.pdb` reference geometry, not a packed trajectory snapshot;
  - NPT handoff updates box attributes on every PBC-aware potential node.
- Validation on copied artifacts under `example/16.MARTINI/outputs/martini_test_1rkl_hybrid_pairpatch_check`:
  - earlier capped and PMF validation is superseded by the energy-expectation table path below;
  - 1RKL stage 6 short run kept 3D CGL-CGL nearest distances physical: final min/p05/mean `6.945/7.182/7.737 Å`;
  - 1RKL stage 7 short run for `5000` steps kept protein Rg in `12.8 -> 11.4 Å` instead of unfolding;
  - stage 7 final 3D CGL-CGL nearest min/p05/mean was `6.739/6.996/7.577 Å`;
  - stage 7 lipid vectors did not become parallel to the bilayer: `abs(n_z)` min/p05/mean at 10 ps was `0.616/0.758/0.922`, with `0` particles below `0.35`.
- Bilayer-only validation with the full-resolution table passed:
  - `bash example/16.MARTINI/test_cg_bilayer/run_test.sh --resolution full --stage60-npt-contract --steps 1000 --frame-steps 25`;
  - final aligned orientation min/p05/mean `0.833/0.885/0.966`;
  - same-leaflet nearest-neighbor min/p05 `6.424/6.731 Å`;
  - CGL-CGLD length min/max/RMSD `10.815/11.716/0.237 Å`.
- Additional bilayer-only validation at the stage-6 timestep passed:
  - `bash example/16.MARTINI/test_cg_bilayer/run_test.sh --resolution full --skip-build --skip-stage --stage60-npt-contract --steps 500 --dt 0.01 --frame-steps 50`;
  - final aligned orientation min/p05/mean `0.836/0.870/0.956`;
  - same-leaflet nearest-neighbor min/p05 `6.007/6.274 Å`.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system.py py/martini_prepare_system_lib.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - `cmake --build obj --target upside` passed;
  - `git diff --check` passed.
- Final corrective review after the user reiterated not to exclude anything:
  - No physical SC-env, BB-env, CGL-SC, CGL-target, or CGL-CGL interaction path was disabled or scaled away.
  - `cg_lipid_target` includes the physical virtual `BB` proxy and environment targets; only duplicate `N/CA/C/O` carrier rows are omitted because they are not separate dry-MARTINI target sites.
  - CGL-target gradients on virtual `BB` now route through the hybrid BB proxy projection, preserving protein-side forces.
  - The fixed `5000 kJ/mol` CGL-CGL core cap and the additive pairwise PMF replacement have both been rejected; CGL-CGL and CGL-target tables now use dry-MARTINI energy-expectation values and first-resolved energy-expectation boundary rows for unresolved radial knots.
  - Previous 50k-step focused hybrid validation of the cap-based prototype kept same-leaflet 3D CGL-CGL nearest distances near `7.6-7.9 Å`, final protein Rg near `11.2 Å`, final hbond sum near `19`, and final aligned CA RMSD near `4.35 Å`; this stability result is not accepted as final evidence because the cap was nonphysical.
  - Previous PMF bilayer validation is rejected because the latest output showed many-body lipid aggregation and overlapping branches.
  - Current stale output `martini.h5` is rejected by the new guard with a rebuild message for the missing CGL-CGL/CGL-target energy-expectation unresolved-core attrs.
  - Fresh bilayer-only validation on the current CGL-CGL model passed for `5000` steps at `dt=0.002` using `bash example/16.MARTINI/test_cg_bilayer/run_test.sh --resolution coarse --steps 5000 --frame-steps 100`; the harness now reports `CG: 16r x 7^2 theta x 2^2 phi` for CGL-CGL even under coarse mode.
  - Final bilayer metrics: aligned-z min/p05/mean `0.431/0.695/0.898`, `bad_parallel=0`, `bad_flip=0`, leaflet crossings `0/0`, same-leaflet nearest-neighbor min/p05 `4.937/5.769 Å`, and CGL-CGLD length RMSD `0.275 Å`.
  - Fresh table schema audit passed with `fit_r_min_nm=0.5`, `fit_relax_steps=50`, `excluded_area_source=wca_dopc_contact_kbt`, `excluded_area_epsilon_kj_mol=2.5205595`, `excluded_area_nonnegative_rows=5`, and no energy-cap attrs.
  - Final cleanup removed inactive cap branches from table generation; Python compile and `git diff --check` passed afterward.

## 2026-05-16 Physical CG-Lipid Orientation Model

### Project Goal
- Fix CG-lipid orientation so bilayer particles remain physical vectors driven by dry-MARTINI-derived pair interactions.
- Keep SC-env, BB-env, CGL-SC, CGL-CGL, CGL-target interactions active; do not use orientation pins, z-restraints, empirical scales, or disabled physics.
- Validate with a bilayer-only simulation that checks particle distribution and orientation.
- Keep `example/16.MARTINI/cg_lipid_potentials.tex` synchronized with the implemented model.

### Architecture & Key Decisions
- Treat `CGLD` as an internal orientation coordinate only. It may carry mass and the DOPC-derived CGL-CGLD bond, but it must not be included in generic nonbonded pairs.
- Exclude both `CGL` and `CGLD` from generic `martini_potential` pair generation. Dedicated spline nodes own every CGL interaction.
- Replace the current low-rank CGL-CGL angular approximation with a full cubic B-spline table over `(r, a1, a2)`, where `a1=-n1.rhat` and `a2=n2.rhat`.
- Replace the radial-only CGL-X table with a directional CGL-target table over `(r, a)`, where `a=n_CGL.rhat`; build it from explicit DOPC bead-target dry-MARTINI energies.
- Remove the optional orientation spring production/debug path from active injection.
- Extend the bilayer-only harness to fail on bad leaflet distribution, near-parallel bilayer lipids, broken CGL-CGLD lengths, and unphysical same-leaflet overlaps.

### Execution Phases
- [x] Phase 1: Patch table generation schemas for full CGL-CGL and directional CGL-target tables.
- [x] Phase 2: Patch runtime C++ nodes and injection logic for full/directional tables and CGL/CGLD pair exclusion.
- [x] Phase 3: Extend bilayer-only diagnostics and pass/fail criteria.
- [x] Phase 4: Update `cg_lipid_potentials.tex` to the new model.
- [x] Phase 5: Run syntax/build/table/bilayer validation and document results.

### Known Errors / Blockers
- None remaining for the bilayer-only validation path.

### Review
- `py/martini_build_tables.py` now writes `cg_lipid_pair_full_v1` as a full tensor spline over `(r, a1, a2)` with `3150` coefficients per CGL-CGL type pair, and writes directional `cg_lipid_target_v1` tables over `(r, a)` with `source=explicit_dopc_directional`.
- Removed obsolete effective CGL generic-particle table extension, pair scaling, orientation spring, and radial `cg_lipid_isotropic` runtime path.
- `py/martini_prepare_system_lib.py` now excludes `CGL` and `CGLD` from generic `martini_potential` pair generation and injects CGL interactions through dedicated spline nodes.
- `src/martini_cg_lipid.cpp` evaluates the full CGL-CGL tensor and CGL-target tensor, propagating orientation derivatives through the CGLD-position-derived vector.
- `example/16.MARTINI/test_cg_bilayer/run_test.py` now fails on near-parallel bilayer lipids, leaflet crossings, broken CGL-CGLD lengths, and unphysical same-leaflet overlaps.
- Bilayer-only validation passed for `1000` steps with the stage-6.0 NPT contract:
  - final aligned orientation `z` min/p05/mean `0.550/0.819/0.952`;
  - `bad_parallel=0`, `bad_flip=0`, leaflet crossings `0/0`;
  - same-leaflet nearest-neighbor min/p05 `6.484/6.811 Å`;
  - CGL-CGLD length reference `11.139 Å`, final min/max/RMSD `10.677/11.548/0.175 Å`.
- HDF5 audit of the regenerated bilayer stage found `0` generic Martini pairs and `0` pairs touching `CGL` or `CGLD`; injected `cg_lipid_pair` carries `schema=cg_lipid_pair_full_v1`, `n_modes=0`, and parameter shape `(1, 1, 3150)`.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - `cmake --build obj --target upside` passed;
  - `bash example/16.MARTINI/test_cg_bilayer/run_test.sh --resolution coarse --stage60-npt-contract --steps 1000 --frame-steps 25` passed;
  - `git diff --check` passed.

## 2026-05-15 VTF Lipid Endpoint Coloring

### Project Goal
- Make CG lipid bars in VMD visibly distinguish hydrophilic and hydrophobic sides under Coloring Method `Name`, `Type`, `Element`, and `ResName`.
- Keep the displayed lipid rod length derived from the starting DOPC PDB geometry without making rods visually overlong.

### Architecture & Key Decisions
- Current VTF uses `name LIPH type HYDROPHILIC resname UNK` for hydrophilic endpoints and `name LIPT type HYDROPHOBIC resname DOPC` for hydrophobic endpoints, so `ResName` and `Element` coloring are not robustly distinct.
- VTF atom records can carry `atomicnumber`; use nitrogen-like hydrophilic endpoints and carbon-like hydrophobic endpoints for `Element` coloring.
- Drop VMD `ResType` support and do not emit companion `.vmd` files.
- Represent each lipid as two same-colored half-rods meeting at the CGL center, using actual DOPC MARTINI bead categories for VMD Martini coloring:
  - hydrophilic half: `PO4/Qa/LIPH/15`;
  - hydrophobic half: `C1A/C1/LIPT/6`.
- Use the DOPC-derived `orientation_length_ang` from the starting PDB as the total visual rod span instead of the full NC3-to-tail endpoint span.

### Execution Phases
- [x] Phase 1: Patch VTF endpoint metadata and atom-record atomic numbers.
- [x] Phase 2: Remove companion VMD script generation and `ResType` support.
- [x] Phase 3: Split lipid rods into hydrophilic and hydrophobic VTF half-rods.
- [x] Phase 4: Regenerate and verify the 1RKL VTF side metadata and shortened span.
- [x] Phase 5: Document results.

### Known Errors / Blockers
- None.

### Review
- Patched `py/martini_extract_vtf.py` so lipid display endpoints carry distinct MARTINI `Name`, `Type`, `ResName`, and `atomicnumber` metadata:
  - hydrophilic endpoint: `name=PO4`, `type=Qa`, `resname=LIPH`, `atomicnumber=15`;
  - hydrophobic endpoint: `name=C1A`, `type=C1`, `resname=LIPT`, `atomicnumber=6`.
- Added `atomicnumber` to all emitted VTF atom records, inferred from atom names/types for non-lipid particles.
- Removed same-stem `.vmd` companion script generation and dropped `ResType` support.
- Split each visible lipid into two side-colored half-rods so VMD can color both halves directly from VTF atom metadata.
- Replaced the full NC3-to-tail display span with the DOPC-derived orientation length from the starting PDB. Regenerated `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` reports `282` lipids with total displayed rod span `11.631003 Å`.
- Verification:
  - `python3 -m py_compile py/martini_extract_vtf.py` passed;
  - regenerated the latest 1RKL VTF;
  - VTF atom audit found `564` hydrophilic `PO4/Qa/LIPH/15` side atoms and `564` hydrophobic `C1A/C1/LIPT/6` side atoms;
  - VTF bond audit found `282` hydrophilic half-bonds and `282` hydrophobic half-bonds.

## 2026-05-15 Restore Preproduction Interface Physics

### Project Goal
- Fix `example/16.MARTINI/run_sim_1rkl.sh` after the latest result in `outputs/martini_test_1rkl_hybrid` shows `stage_6.0` creates protein-lipid clashes before production.
- Keep SC-env, BB-env, and CGL-SC interface interactions active; do not disable, zero, or tune physical interactions to hide the failure.
- Preserve the restored rigid-protein `6.0` NPT relaxation stage and the direct `6.0 -> 7.0` handoff for current CG-lipid systems.

### Architecture & Key Decisions
- Root cause: `stage_6.0.up` has `cg_lipid_pair` and `martini_potential`, but lacks `martini_sc_table_1body` and `cg_lipid_sc`; `stage_7.0.up` restores them after `6.0` already moved lipids into the rigid protein.
- Use the same hybrid interface injection path for preproduction and production stages, while keeping stage labels/activation modes distinct:
  - preproduction: `current_stage=minimization`, `activation_stage=minimization`, rigid protein;
  - production: `current_stage=production`, `activation_stage=production`.
- Add assertions that hybrid stages with protein and CG lipids fail if required interface nodes are missing.

### Execution Phases
- [x] Phase 1: Refactor hybrid interface node injection so preproduction stages get `martini_sc_table_1body` and `cg_lipid_sc`.
- [x] Phase 2: Add required-node assertions for hybrid preproduction and production stages.
- [x] Phase 3: Verify regenerated `6.0` no longer creates the observed CGL-protein clash.
- [x] Phase 4: Document results and residual risks.

### Known Errors / Blockers
- None remaining for the reproduced `6.0` handoff failure.

### Review
- Added shared hybrid interface injection in `py/martini_prepare_system.py` so preproduction and production stages both carry `martini_sc_table_1body`, `cg_lipid_sc`, `cg_lipid_pair`, and `martini_potential`.
- Kept preproduction as `current_stage=minimization`, `activation_stage=minimization`, and rigid-protein mode; production remains `current_stage=production`, `activation_stage=production`.
- Added required-node assertions so hybrid stages with CG lipids fail before MD if physical interface nodes are missing.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed;
  - regenerated `/private/tmp/upside_preprod_fix_full/checkpoints/1rkl.stage_6.0.up` contains all required interface nodes;
  - `6.0` CGL-protein minimum distance is `7.608 Å` at input and `5.357 Å` after 500 steps, instead of the latest bad `2.275 Å`;
  - `6.0` 3D CGL-CGL nearest distance improves from `6.233 Å` to `7.126 Å`;
  - 10000-step production smoke ends at Rg `13.2 Å` with observed range `12.8-14.1 Å`, not the bad `17.9 Å` drift;
  - production `cg_lipid_sc` starts at `18.88 E_up` and stays in `[-14.70, 18.88] E_up`, not the bad `1152.51 E_up` startup spike;
  - `git diff --check` passed.

## 2026-05-15 CG-Lipid Physical Parameter Rationalization

### Project Goal
- Remove current CG-lipid "twisting" parameters from the 1RKL MARTINI workflow where they can be derived from dry-MARTINI ITP/PDB data.
- Keep the single-particle CGL model, active SC-env/BB-env interactions, and spline-table runtime path intact.
- Preserve SVD/table resolution choices as numerical approximation controls.

### Architecture & Key Decisions
- Add a shared DOPC-derived geometry/force-field helper used by table building and stage conversion.
- Derive CGL contact spacing from dry-MARTINI nonbonded sigma (`2^(1/6) * max_sigma`) instead of fixed `7.0 Å`.
- Treat CGLD as an internal orientation-coordinate carrier:
  - derive its length from the DOPC COM-to-tail midpoint projection;
  - derive its mass from DOPC transverse rotational inertia;
  - derive its bond force constant from projected DOPC bonded stiffness.
- Remove the default `cg_lipid_orientation_spring`; do not pin lipid direction to the initial structure in normal runs.
- Replace the global CGL-vs-non-CGL sigma cap with explicit orientation-averaged CGL-X radial spline tables.
- Keep taper width equal to one radial knot interval and keep SVD/mode count as resolution controls.
- Add stale-table guards so old `martini.h5` files with capped CGL LJ metadata are rejected for fresh stage injection.

### Execution Phases
- [x] Phase 1: Add DOPC derivation helper and route table/stage code through it.
- [x] Phase 2: Patch CGL orientation mechanics, initial spacing, and CGL-X isotropic table generation.
- [x] Phase 3: Remove default orientation spring injection and add stale-table compatibility guard.
- [x] Phase 4: Run syntax/build/table/workflow smoke checks.
- [x] Phase 5: Document review and remaining risks.

### Known Errors / Blockers
- Existing generated `martini.h5` files may lack the new derivation attrs and should be rebuilt rather than silently reused.

### Review
- Added `py/martini_cg_lipid_params.py` as the shared DOPC-derived parameter helper for atom names, bead types, ITP masses, contact distance, CGLD geometry/mass, and projected orientation-bond stiffness.
- Patched table generation so CG-CG and CG-SC fitting use the derived DOPC contact radius, CGL-X isotropic tables are explicit DOPC orientation averages, and generated tables carry `dry_martini_dopc_derived_v1` attrs.
- Patched stage conversion so CGLD length, mass, orientation bond stiffness, and same-leaflet conditioning default to derived physical values; environment overrides remain explicit overrides.
- Removed default orientation-spring injection by setting `UPSIDE_CG_LIPID_ORIENT_SPRING_K` default to `0.0`.
- Added stale-table validation that rejects old capped CGL effective-LJ tables before reuse or injection.
- Verification:
  - `python3 -m py_compile py/martini_cg_lipid_params.py py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py` passed;
  - DOPC helper sanity check gave `contact_ang=6.959265`, `orientation_length_ang=11.139272`, `orientation_mass_g_mol=77.048875`, `orientation_bond_fc_eup_a2=39.435978`, `max_sigma_nm=0.620000`;
  - coarse table smoke built `/private/tmp/cg_lipid_derived_smoke.h5` with `cg_lipid_isotropic.source=explicit_dopc_orientation_average`, no `max_effective_sigma_nm`, and derived attrs;
  - stale-table guard rejected a synthetic capped legacy table;
  - `INITIAL_DEBUG_ONLY=1 RUN_DIR=outputs/cg_lipid_derived_initial_debug PREP_SEED=12345 SEED=12345 bash example/16.MARTINI/run_sim_1rkl.sh` completed, adding CGLD sites with derived defaults and conditioning same-leaflet spacing to `6.959 Å`;
  - `git diff --check` passed.

## 2026-05-07 1RKL CG-SC Range Stabilization

### Project Goal
- Stop the hybrid 1RKL protein collapse caused by nonphysical long-range CG lipid-sidechain attraction.
- Keep the directional single-particle lipid model and protein/environment interface interactions enabled.
- Make production logs distinguish protein stability from CG lipid energetic terms.

### Architecture & Key Decisions
- The DOPC-sidechain table is fitted only through `fit_r_max_nm=0.7` (`7.0 Å`), so runtime use must not extend the attractive fitted table to the old `16.8 Å` cutoff.
- Use a compact physical support of `fit_r_max_nm * 10 + taper_width_ang`, which is `8.4 Å` with the current `1.4 Å` taper.
- Keep the DOPC-DOPC `cg_lipid_pair` cutoff at its existing support; this fix is specific to `cg_lipid_sc`.
- Add a stale-table guard during HDF5 stage injection so old `martini.h5` files with `cg_lipid_sc.cutoff_ang=16.8` are capped without requiring users to rebuild tables before debugging.
- Update logging so `cg_lipid_pair` and `cg_lipid_sc` no longer inflate the reported protein potential.

### Execution Phases
- [x] Phase 1: Patch CG lipid-sidechain table generation and stale-table injection cutoff handling.
- [x] Phase 2: Add debug summary visibility for the active `cg_lipid_sc` cutoff and CGL-protein radial density.
- [x] Phase 3: Patch runtime log component accounting for protein, CGL-CGL, and CGL-SC potentials.
- [x] Phase 4: Run syntax/build and focused stale-table/short-run validation.

### Known Errors / Blockers
- Existing 1RKL output uses `cg_lipid_sc.cutoff_ang=16.8` even though the fitted SC radial support is only `7.0 Å`, causing lipid particles to over-attract the protein and collapse it over time.
- Current `prot_potential` logging includes `cg_lipid_pair` and `cg_lipid_sc`, obscuring whether the protein itself is destabilizing.

### Review
- Patched `py/martini_build_tables.py` so newly generated `cg_lipid_sc` tables store cutoff as fitted support plus taper (`8.4 Å` for current 1RKL/DOPC-sidechain settings), while leaving `cg_lipid_pair` unchanged.
- Patched `py/martini_prepare_system_lib.py` so stale `martini.h5` files with `cg_lipid_sc.cutoff_ang=16.8` inject into stages with the capped `8.4 Å` runtime cutoff.
- Added debug summary fields for active `cg_lipid_sc` attrs and CGL-near-protein distance counts.
- Patched `src/main.cpp` hybrid progress logging to report `protein_potential`, `cg_lipid_pair`, and `cg_lipid_sc` separately.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py` passed;
  - `cmake --build obj --target upside` passed with pre-existing warnings only;
  - stale-table injection test converted old table attrs `cutoff_ang=16.8`, `fit_r_max_nm=0.7`, `taper_width_ang=1.4` into injected `cutoff_ang=8.4`;
  - generated stage debug JSON reports the injected `cg_lipid_sc` cutoff and protein-near CGL counts;
  - 100-step smoke printed separated terms and started from `cg_lipid_sc=-3.65 E_up`;
  - 10000-step temp validation kept protein Rg near `12.3-12.8 Å`, ending at `12.7 Å`; `cg_lipid_sc` ended near `-100.27 E_up` instead of the previous multi-thousand `E_up` collapse driver;
  - `git diff --check` passed.

## 2026-05-12 Restore 6.0 NPT Relaxation and Two-Endpoint Lipid VTF

### Project Goal
- Restore an active rigid-protein `6.0` NPT relaxation stage before production for the current CG-lipid workflows.
- Preserve the longer pre-production route for systems with real bonded dry-MARTINI environment particles.
- Replace the visible one-sided lipid VTF ghost representation with hydrophilic/hydrophobic endpoint particles connected by a bond.

### Architecture & Key Decisions
- `6.0` becomes a real MD/NPT stage in the shared workflow path, using the existing non-production rigid-body protein machinery and semi-isotropic barostat handoff.
- Add an explicit `EQ_60_NSTEPS` / `--eq-60-nsteps` control; default it to `500` to match the current short equilibration stage scale.
- Keep the dormant `6.1-6.6` path, but trigger it on any real bonded dry-MARTINI environment particles, not only lipids. Synthetic `CGL-CGLD` orientation bonds do not count.
- Store CG lipid display endpoint offsets derived from the original DOPC reference head/tail span so VTF extraction can render a centered hydrophilic-hydrophobic vector at the physical reference length.
- Keep production-stage NPT behavior unchanged.

### Execution Phases
- [x] Phase 1: Restore mandatory `6.0` NPT scheduling, configuration surface, and bonded-environment routing.
- [x] Phase 2: Persist lipid display endpoint metadata and replace VTF ghost rendering with two typed endpoint particles plus a bond.
- [x] Phase 3: Extend focused regression coverage for `6.0`, bonded-path detection, VTF topology, and bilayer shrink reporting.
- [x] Phase 4: Run syntax, targeted workflow, VTF, and bilayer comparison verification.

### Known Errors / Blockers
- No archived bilayer-only numerical shrink artifact from `3a98be1a4b4ebbbdf645fd1db1dcb84efa86af1e` exists in-tree, so the comparison was made against that commit's recovered stage/barostat contract plus a fresh local bilayer-only probe.

### Review
- Fresh `1rkl` workflow now runs mandatory `6.0` NPT MD, keeps the protein rigid during that stage, and hands both coordinates and box into `7.0`.
- The extended `6.1-6.6` path remains available only when real bonded dry-MARTINI environment pairs are present; synthetic `CGL-CGLD` orientation bonds are ignored.
- VTF extraction now emits typed hydrophilic/hydrophobic lipid endpoints with one bond per lipid and a mean first-frame displayed span of `24.270467 Å` in the reduced `1rkl` smoke artifact.
- Bilayer-only `stage60`-contract smoke used the recovered zero-target semi-isotropic Berendsen settings and showed effectively flat short-horizon XY drift (`+0.008194 Å` over 500 steps, `z` unchanged), so there was no meaningful box shrink to claim from that isolated probe.

## 2026-05-12 Direct Bilayer-Only NPT Comparison vs Dry-MARTINI

### Project Goal
- Build a bilayer-only bead-resolved dry-MARTINI NPT probe under `/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI`.
- Run it under settings matched to the current single-particle bilayer probe and compare `Lx/Ly` relaxation directly.

### Architecture & Key Decisions
- Add only the minimal dry-MARTINI lipid-only stage-conversion support needed to prepare a DOPC-only input.
- Keep the comparison probe separate from the hybrid protein workflow so it does not alter production orchestration.
- Reuse the same zero-target semi-isotropic Berendsen NPT contract, timestep, frame spacing, temperature, and step count used by the current single-particle probe.
- Compare box trajectories numerically from HDF5 output, not from console logs alone.

### Execution Phases
- [x] Phase 1: Add the dry-MARTINI lipid-only preparation probe in the external repo.
- [x] Phase 2: Run matched dry-MARTINI and single-particle bilayer NPT jobs.
- [x] Phase 3: Extract XY box trajectories and decide whether the single-particle model matches dry-MARTINI relaxation.

### Known Errors / Blockers
- None for this comparison pass.

### Review
- Added `/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/test_dry_bilayer/` plus the minimal lipid-only converter guard needed by the external dry-MARTINI repo.
- Matched `500`-step NPT probes used the same `50.091999 Å` initial XY side length, `dt=0.002`, frame spacing `25`, temperature `0.8647`, and zero-target semi-isotropic Berendsen contract.
- The dry-MARTINI bead-resolved bilayer ended at `50.071121 Å` in XY (`-0.020878 Å` side change, `-0.083341%` XY area).
- The current single-particle bilayer ended at `50.100193 Å` in XY (`+0.008194 Å` side change, `+0.032718%` XY area).
- The signs differ over the matched probe, so the current single-particle bilayer NPT box response is not consistent with the bead-resolved dry-MARTINI reference on this test.

## 2026-05-07 1RKL VTF Ghost Mapping and CGL-Ion Startup Stability

### Project Goal
- Fix the strange first-frame VTF visualization in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid`.
- Remove the artificial single-particle lipid/ion startup clash that destabilizes `run_sim_1rkl.sh`.
- Add diagnostics that make CGL-ion/protein proximity visible before running expensive pair table generation and dynamics.

### Architecture & Key Decisions
- Treat the current HDF5 stage input as the source of truth: frame 0 CGL coordinates are slab-like and ions are not truly bonded.
- Fix VTF CG lipid ghost atoms by remapping original HDF5 CGL indices to output atom indices before adding CGH atoms and bonds.
- Keep the directional CGL-CGL spline table and protein SC/BB environment interactions enabled.
- Cap only the isotropic effective CGL LJ used against non-CGL dry-MARTINI particle types; the current orientation-averaged lipid-vs-point fit creates a nonphysical `~1.8 nm` sigma for ions.
- Increase the default 1RKL ion placement cutoff to match the capped CGL effective radius.

### Execution Phases
- [x] Phase 1: Patch VTF CGH atom/bond/position remapping.
- [x] Phase 2: Cap non-CGL effective CGL LJ parameters and annotate table attrs.
- [x] Phase 3: Add CGL-ion/protein debug metrics to generated JSON summaries.
- [x] Phase 4: Update 1RKL default ion cutoff and validate scripts.
- [x] Phase 5: Re-extract current VTF and verify no CGH bonds originate from ions/protein.
- [x] Phase 6: Regenerate/debug-check a 1RKL initial input and verify the CGL-ion LJ spike is removed.
- [x] Phase 7: Add an optional rigid-protein debug mode for production lipid-stability isolation.

### Known Errors / Blockers
- Current output VTF uses original HDF5 atom indices for CGL ghost bonds after mode-2 remapping, so CGH atoms appear attached to NA/CL/protein atoms.
- Current CGL-ion effective LJ has sigma near `17.9 Å`; ions placed at `~6-8 Å` produce an initial ion-lipid energy spike of roughly `145k E_up`.

### Review
- Fixed `py/martini_extract_vtf.py` so CGH direction-marker atoms use original HDF5 CGL/CGLD vectors but attach to remapped output CGL atom indices.
- Capped non-CGL effective CGL LJ sigma at `UPSIDE_CG_LIPID_MAX_EFFECTIVE_SIGMA_NM=0.9` by default in table generation and in stage conversion, including stale `martini.h5` reuse.
- Changed default `ION_CUTOFF` for the 1RKL workflows and direct `martini_prepare_system.py` use from `4.0 Å` to `10.0 Å`.
- Added `cgl_partner_stats` to debug summaries, including nearest CGL-ion/protein distances and generated `martini_potential` LJ totals/maxima when pair tables are present.
- Added `DEBUG_RIGID_PROTEIN=1` / `--debug-rigid-protein 1`, which injects `/input/fix_rigid/atom_indices` for protein atoms in production stage files without disabling protein-environment interactions.
- Verification:
  - `python3 -m py_compile py/martini_extract_vtf.py py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py` passed;
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh example/16.MARTINI/run_sim_1rkl_outlipid.sh example/16.MARTINI/run_sim_1afo.sh example/16.MARTINI/run_sim_1afo_outlipid.sh` passed;
  - re-extracted current `1rkl.stage_7.0.up` to `/tmp/1rkl.stage_7.0.fixed.vtf`; all `282` CGH bonds are `CGL-CGH` and `0` CGH bonds originate from ions/protein;
  - current unstable stage debug summary now reports the old root cause: CGL-ion LJ max `97404.794 E_up`, CGL-ion LJ sum `145929.948 E_up`, and min CGL-ion distance `6.366 Å`;
  - with the sigma cap applied to the same old geometry, the estimated CGL-ion LJ sum drops to `22.963 E_up` and max pair to `22.167 E_up`;
  - `INITIAL_DEBUG_ONLY=1 RUN_DIR=outputs/debug_1rkl_ljfix_initial PREP_SEED=12345 SEED=12345 bash example/16.MARTINI/run_sim_1rkl.sh` completed and generated debug PDB/JSON without `martini.h5`;
  - new initial-debug JSON reports min CGL-ion distance `21.028 Å`, no warnings, and no ion bonds;
  - direct rigid debug injection fixed `155` protein atoms in a copied stage file;
  - `git diff --check` passed.

## 2026-05-07 DOPC-DOPC vs DOPC-Sidechain Directional Table Audit

### Project Goal
- Verify whether the DOPC-sidechain coarse-grained table calculation uses the same directional method as the single-particle DOPC-DOPC interaction.
- Identify any mismatch in angular sign convention, radial spline construction, orientation relaxation/sampling, unit conversion, or runtime evaluation.

### Architecture & Key Decisions
- Treat DOPC-DOPC as the reference implementation for directional single-particle lipid interactions.
- Compare generation code and runtime code, not just HDF5 schema names.
- Do not modify physics until the mismatch, if any, is proven from source and/or generated table attrs.
- Revised Decision: DOPC-sidechain must use the same full multimode directional table layout and runtime evaluator style as DOPC-DOPC, with SC-specific rotamer averaging retained.

### Execution Phases
- [x] Phase 1: Locate DOPC-DOPC and DOPC-sidechain table generation paths.
- [x] Phase 2: Compare angular sign convention, radial basis, mode schema, relaxation, and sampling.
- [x] Phase 3: Compare runtime evaluator paths and injected HDF5 attrs.
- [x] Phase 4: Report whether they are equivalent or patch the mismatch.

### Known Errors / Blockers
- Source audit found a real mismatch: DOPC-DOPC uses `cg_lipid_quadspline_v3` full multimode params, but DOPC-sidechain still uses fixed 54-parameter `cg_lipid_sc_quadspline_v1`.

### Review
- DOPC-DOPC reference path:
  - `_fit_cg_lipid_quadspline()` builds `cg_lipid_quadspline_v3`;
  - table layout is full multimode: `V0(r) + sum_m Ang1_m(a1) * Ang2_m(a2) * Vm(r)`;
  - runtime uses `eval_multimode_pair()` with `ang1=-n1_dot_n12;ang2=n2_dot_n12`.
- Original DOPC-sidechain path was not equivalent:
  - `_fit_cg_lipid_sc_quadspline()` emitted fixed 54-parameter `cg_lipid_sc_quadspline_v1`;
  - runtime used the older single-mode `eval_quadspline()` path;
  - radial knots/default cutoff differed from the DOPC-DOPC full multimode table.
- Patched DOPC-sidechain to use the same full multimode parameter layout and runtime evaluator style:
  - `cg_lipid_sc_quadspline_v2`;
  - full multimode variable-length params;
  - same angle convention and source order;
  - same 15 angular knots, 14 radial knots, 1.4 Å radial knot spacing, 16.8 Å cutoff, and taper attrs;
  - retained SC-specific rotamer averaging and lipid-against-fixed-SC relaxation.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py` passed;
  - `cmake --build obj --target upside` passed;
  - direct coarse table-generation check produced `cg_lipid_pair` shape `(1,1,278)` with `schema=cg_lipid_quadspline_v3`, `n_modes=6`, `n_radial=14`, `n_angular=15`;
  - the same check produced `cg_lipid_sc` shape `(1,1,190)` with `schema=cg_lipid_sc_quadspline_v2`, `radial_mode=full_multimode`, `angle_convention=ang1=-n1_dot_n12;ang2=n2_dot_n12`, `n_modes=4`, `n_radial=14`, `n_angular=15`;
  - `git diff --check` passed.


## 2026-05-06 Initial Structure Debug PDB Export

### Project Goal
- Generate visualization-ready PDB diagnostics for every generated MARTINI input so incorrect initial structures can be debugged before simulation.
- Prove whether apparent ion bonds are real HDF5 topology or visualization artifacts.
- Make hidden CGLD orientation sites inspectable without confusing normal bilayer visualization.
- Generate a dedicated single-particle lipid bilayer PDB containing only physical CGL lipid particles.

### Architecture & Key Decisions
- Add a reusable stage-input debug exporter in `py/martini_prepare_system_lib.py`.
- Emit debug files under the simulation output directory `debug/`, keyed by the `.up` filename stem.
- Write two PDBs: all particles including CGLD, and a visible model excluding CGLD.
- Write a third bilayer-only PDB containing only physical `CGL` particles.
- Do not write PDB `CONECT`; actual topology is reported in a TSV from `/input/potential/dist_spring/id`.
- Fix CGLD `molecule_ids` to match parent CGL particles so molecule grouping is not misleading.

### Execution Phases
- [x] Phase 1: Implement debug PDB/TSV/JSON export from generated `.up` inputs.
- [x] Phase 2: Wire debug export into simple stage conversion and finalized hybrid checkpoint generation.
- [x] Phase 3: Fix CGLD molecule metadata.
- [x] Phase 4: Run syntax and small bilayer input generation verification.
- [x] Phase 5: Inspect generated debug reports and document whether ions are truly bonded.
- [x] Phase 6: Add and verify a dedicated single-particle lipid bilayer PDB export.
- [x] Phase 7: Add a fast `run_sim_1rkl.sh` initial-debug mode that writes PDBs after packing/CG conversion and exits before MARTINI table generation, pair-list generation, and dynamics.
- [x] Phase 8: Export debug PDBs for the actual `run_sim_1rkl.sh` stage files passed to simulation during full runs, and make the default apply to all example/16.MARTINI shell workflows.
- [x] Phase 9: Emit the stage-7 single-particle lipid PDB in normal workflow runs immediately after packing, not only under `INITIAL_DEBUG_ONLY=1`.

### Known Errors / Blockers
- Existing `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` was generated before the CGLD molecule-id fix; its debug summary now reports `282` CGLD molecule-id mismatches. Regenerate the workflow input to clear this stale metadata.
- User correction: the first implementation verified the bilayer-only test and prepared-stage export, but did not prove that `run_sim_1rkl.sh` emits a bilayer PDB for the exact stage-7 simulation input.
- User correction: full stage conversion is too slow for visual initial-structure debugging because it computes MARTINI pair interactions; the workflow needs a PDB-only path that skips `martini.h5` generation/use and pair-list generation.
- User correction: normal `run_sim_1rkl.sh` and inherited workflows still need to generate the single-particle lipid PDB without requiring users to set `INITIAL_DEBUG_ONLY=1`.

### Review
- Added `write_stage_debug_outputs()`:
  - writes `<stem>.input_debug.visible.pdb`, `<stem>.input_debug.all.pdb`, `<stem>.input_debug.bilayer.pdb`, `<stem>.input_debug.bonds.tsv`, and `<stem>.input_debug.summary.json`;
  - omits PDB `CONECT` so viewers do not invent misleading bond topology from the debug PDB;
  - records true `dist_spring` topology, CGLD vector lengths, lipid leaflet stats, particle counts, and warnings.
- Wired debug export into `convert_stage()`, finalized hybrid `prepare_stage_file()`, and the bilayer-only test workflow.
- Normal hybrid workflow runs now call the same fast stage-7 initial-debug export immediately after packing, before MARTINI table building.
- Added `INITIAL_DEBUG_ONLY=1` / `--initial-debug-only` for `run_sim_1rkl.sh` and inherited workflows (`1rkl_outlipid`, `1afo`, `1afo_outlipid`):
  - performs real workflow packing and CG lipid conversion;
  - writes stage-7 debug PDB/TSV/JSON files;
  - exits before MARTINI table generation, HDF5 pair-list generation, and dynamics.
- Full runs also export debug PDBs immediately before each MD stage from the exact `.up` file passed to `upside`.
- Fixed generated `molecule_ids` so each CGLD orientation site shares the parent CGL lipid molecule.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system_lib.py py/martini_prepare_system.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh example/16.MARTINI/run_sim_1rkl_outlipid.sh example/16.MARTINI/run_sim_1afo.sh example/16.MARTINI/run_sim_1afo_outlipid.sh` passed;
  - `python example/16.MARTINI/test_cg_bilayer/run_test.py --skip-build --skip-run` generated debug PDB/TSV/JSON files;
  - bilayer-only debug summary reports `ion_bond_count=0`, `bond_count=72`, `cgl_cgld_bond_count=72`, `molecule_id_count=72`, and no warnings;
  - bilayer-only `<stem>.input_debug.bilayer.pdb` contains exactly `72` physical `CGL` atom records, no `CGLD`, and no ions;
  - `INITIAL_DEBUG_ONLY=1 RUN_DIR=outputs/debug_1rkl_initial_only bash example/16.MARTINI/run_sim_1rkl.sh` generated 1RKL stage-7 debug files without `martini.h5`;
  - generated 1RKL bilayer PDB contains `282` physical `CGL` atom records, no `CGLD`, no ions;
  - generated 1RKL summary reports `ion_bond_count=0`, `bond_count=282`, all bonds are CGL-CGLD, and no warnings.
  - normal `run_sim_1rkl.sh` with `RUN_DIR=outputs/debug_1rkl_normal_generates_pdb` generated `debug/1rkl.stage_7.0.input_debug.bilayer.pdb` before table building; the file contains `282` CGL atom records, no `CGLD`, no ions.

## 2026-05-06 Full Directional CG Lipid Bilayer Stabilization

### Project Goal
- Stabilize the single-particle DOPC bilayer over time without permanent leaflet restraints.
- Replace the unstable combination of large isotropic CGL-CGL LJ plus one residual directional mode with a full directional CGL-CGL spline table.
- Preserve the single physical CGL particle plus hidden rotatable CGLD orientation carrier.

### Architecture & Key Decisions
- Use `cg_lipid_pair` as the only physical CGL-CGL interaction when CG lipid tables are present.
- Keep CGL interactions with water, ions, protein dry-MARTINI sites, SC-env, and BB-env intact.
- Fit a full multi-mode directional table:
  `V(r,a1,a2) = V0(r) + sum_m A_m(a1) * B_m(a2) * Vm(r)`.
- Default to four residual SVD modes, explicit tail-tail and side-by-side cohesion modes, clamped radial splines, and the existing angular spline basis.
- Keep old static `direction` fallback, but require regenerated CG lipid table artifacts for the new CGL-CGL pair schema.
- Add initial CGL geometry conditioning only as preparation, not as a production-time restraint.

### Execution Phases
- [x] Phase 1: Patch table generation for full multi-mode CGL-CGL fitting and validation diagnostics.
- [x] Phase 2: Patch C++ `cg_lipid_pair` runtime to consume full multi-mode params and preserve derivative propagation.
- [x] Phase 3: Patch stage injection so CGL-CGL `martini_potential` coefficients are zeroed when `cg_lipid_pair` is active.
- [x] Phase 4: Add CGL-only initial bilayer conditioning and stronger bilayer-only diagnostics.
- [x] Phase 5: Run syntax/build/table/bilayer verification and document final behavior.

### Known Errors / Blockers
- Longer and larger bilayer validation remains a calibration problem, but the 72-DOPC bilayer-only test no longer shows monotonic thickness spreading through `50000` steps.

### Review
- Implemented `cg_lipid_quadspline_v3`:
  - CGL-CGL pair params are now variable-length full multimode splines;
  - table contains `V0(r)`, four residual SVD modes, one explicit tail-tail cohesion mode, and one explicit same-leaflet side-by-side cohesion mode;
  - `kappa` scaling applies only to energy-valued radial channels, not dimensionless angular splines.
- Removed the unstable CGL-CGL isotropic core:
  - CGL-CGL rows in `martini_potential` are zero when CG lipid pair tables are active;
  - CGL interactions with other dry-MARTINI targets remain table-driven.
- Added preparation-only CGL spacing conditioning:
  - same-leaflet XY nearest-neighbor min/p05 improved from `5.347/5.524 Å` to `7.000/7.000 Å`;
  - CGLD orientation sites move with their CGL partners.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - `cmake --build obj --target upside` passed;
  - regenerated bilayer-only table stores `schema=cg_lipid_quadspline_v3`, `n_modes=6`, `n_radial=14`, and `1×1×278` CGL-CGL params;
  - 100-step bilayer-only smoke stayed finite with initial potential `-39.46 E_up`;
  - 50000-step bilayer-only run stayed finite and slab-like: `z_std` changed from `6.512 Å` to `5.968 Å` at time `100`, with unit-normalized dynamic vectors.

## 2026-05-06 Dynamic CG Lipid Orientation Sites

### Project Goal
- Allow single-particle DOPC lipid vectors to rotate during simulation instead of remaining fixed at their initial directions.
- Preserve the one physical CGL COM particle model while adding a hidden orientation carrier that receives torque through the existing directional potential.

### Architecture & Key Decisions
- Add one hidden noninteracting orientation site per CGL lipid (`CGLD`) during stage conversion.
- Compute each runtime lipid direction as the normalized minimum-image displacement from CGL COM to its orientation site when `compose_vector6d/orientation_index` exists.
- Propagate direction derivatives back to both CGL COM and orientation site through the normalized-vector Jacobian.
- Keep existing static `compose_vector6d/direction` behavior as a fallback for old stage files.
- Constrain the CGL-CGLD distance with a harmonic bond using the existing bond infrastructure.

### Execution Phases
- [x] Phase 1: Inspect existing bond and pair table infrastructure for a minimal dynamic-orientation implementation.
- [x] Phase 2: Patch stage conversion to append CGLD sites, masses, zero nonbonded interactions, and CGL-CGLD bonds.
- [x] Phase 3: Patch `compose_vector6d` runtime to derive and backpropagate dynamic directions.
- [x] Phase 4: Update bilayer-only diagnostics to report direction rotation.
- [x] Phase 5: Run syntax/build and bilayer-only smoke verification, then record results.

### Known Errors / Blockers
- The dynamic orientation carrier now receives torque, but the bilayer-only smoke still spreads in `z` over `5000` steps; remaining instability is a force-model/calibration issue rather than a fixed-vector implementation bug.
- The dynamic direction minimum-image box in `compose_vector6d` is read from stage attrs; if future NPT runs enable large box changes, this path should be upgraded to consume live box updates.

### Review
- Added hidden `CGLD` orientation sites during CG lipid conversion:
  - one site per CGL particle;
  - zero charge and zero nonbonded interactions;
  - harmonic CGL-CGLD bond through existing `dist_spring`;
  - orientation length and bond force constant configurable via `UPSIDE_CG_LIPID_ORIENTATION_LENGTH` and `UPSIDE_CG_LIPID_ORIENTATION_BOND_FC`.
- Updated `compose_vector6d` so runtime lipid directions are normalized CGL-to-CGLD displacements when `orientation_index` exists, with direction derivatives propagated back to both particles.
- Updated bilayer-only and VTF diagnostics to derive displayed CGL directions from orientation sites when present.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system_lib.py py/martini_extract_vtf.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - `cmake --build obj --target upside` passed;
  - 100-step bilayer-only dynamic-vector smoke stayed finite and reported direction rotation up to `1.122 deg`;
  - 5000-step bilayer-only dynamic-vector run stayed finite, directions remained normalized, and direction rotation reached `149.876 deg` max / `65.677 deg` mean;
  - 5000-step `z_std` improved from the previous static-vector `13.384 Å` to `10.390 Å`, but the bilayer is not yet physically stable.

## 2026-05-06 CG Lipid Bilayer-Only Isolation

### Project Goal
- Verify whether the single-particle DOPC instability starts from bad initial placement or from the CG lipid force model itself.
- Add a smaller bilayer-only test for single-particle DOPC lipids and use it to assess short-horizon stability without protein/interface terms.

### Architecture & Key Decisions
- First audit the latest `run_sim_1rkl.sh` generated stage files: CGL positions, direction vectors, bilayer leaflet separation, XY nearest-neighbor spacing, box attrs, and whether CG lipid nodes use the intended table attrs.
- Build the bilayer-only example using existing MARTINI preparation/conversion utilities rather than inventing a separate file format.
- Keep the base CGL-CGL LJ and `cg_lipid_pair` active in the isolated test; do not disable the CG lipid directional correction just to make the test pass.
- Use short smoke simulations for local diagnosis, then document whether a longer run remains required for physical stability.

### Execution Phases
- [x] Phase 1: Inspect current 1RKL output and verify initial single-particle lipid placement.
- [x] Phase 2: Add or adapt a bilayer-only single-particle DOPC test workflow.
- [x] Phase 3: Run the bilayer-only workflow and short stability smoke test.
- [x] Phase 4: Analyze bilayer-only trajectory geometry/energy and identify the next root cause.
- [x] Phase 5: Record findings, verification, and residual blockers.

### Known Errors / Blockers
- The latest inspected 1RKL artifact has bad CGL placement: upper-leaflet CGL z reaches below the midplane and same-leaflet XY nearest-neighbor distance drops to about `2.245 Å`.
- A bilayer-only NVT test stays finite after `5000` steps with the corrected table, but the single-particle model spreads in `z`; long-lived bilayer morphology still needs a physical stabilization target beyond this numerical smoke test.

### Review
- Added `example/16.MARTINI/test_cg_bilayer/run_test.py` and `run_test.sh` for a stripped 72-DOPC bilayer-only single-particle workflow.
- Fixed lipid-only stage conversion in `py/martini_prepare_system_lib.py` so DOPC-only MARTINI PDBs do not require fake protein backbone atoms.
- Fixed CG lipid table generation in `py/martini_build_tables.py`:
  - angular basis knots remain dimensionless;
  - radial B-spline fitting now uses the same clamped deBoor basis as runtime;
  - `V_angular(r)` is fit by least-squares projection and clipped to the residual cap instead of averaging residual/product ratios.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed;
  - bilayer-only table rebuild produced `radial maxabs=0`, `V_angular maxabs=19.5169 E_up`;
  - 100-step DOPC-only smoke stayed finite with max CGL displacement `0.100 Å`;
  - 5000-step DOPC-only run stayed finite, with `z_std` growing from `6.512 Å` to `13.384 Å`.

## 2026-05-05 CG Lipid Quadspline Bilayer Stabilization

### Project Goal
- Stabilize the single-particle directional CG lipid workflow used by `example/16.MARTINI/run_sim_1rkl.sh`.
- Make the generated CG lipid spline table and runtime evaluator match the intended directional potential:
  `V = kappa * (V_radial(r12) + Ang1(-n1*n12) * Ang2(n2*n12) * V_angular(r12))`.
- Preserve the active dry-MARTINI environment interactions; do not disable SC-env, BB-env, or the base MARTINI lipid interactions.

### Architecture & Key Decisions
- Keep the dry-MARTINI CGL-CGL LJ interaction in `martini_potential` as the isotropic excluded-volume core.
- Convert `cg_lipid_pair` into a residual directional correction by zeroing its fitted isotropic radial component and fitting only the orientation-dependent angular residual.
- Do not relax lipid internal beads while fitting CG-CG tables; relaxation can erase the repulsive wall and overfit an already represented isotropic term.
- Fit the CG-CG residual outside the CGL-CGL core over `7-14 Å` with a CG-specific radial knot spacing, rather than training a huge angular correction on bead-overlap configurations.
- Bound the CG-CG directional residual during SVD fitting so explicit-bead overlap outliers are left to the isotropic CGL-CGL core instead of becoming enormous angular spline coefficients.
- Train angular samples in the same pair-axis sign convention used at runtime: lipid 1 uses `-n1 dot n12`, lipid 2 uses `n2 dot n12`.
- Replace the generic CG lipid `InteractionGraph` runtime path with CG-specific evaluators that use minimum-image periodic distances and the explicit sign convention above.
- Apply a finite-range taper to the CG residual quadspline correction so the right-clamped radial spline cannot become a constant all-to-all attraction beyond the fitted range.
- Treat CG lipid directions as static vector descriptors for this stabilization pass. Direction derivatives are accumulated where the upstream node can consume them, but `compose_vector6d` still drops lipid-direction derivatives.

### Execution Phases
- [x] Phase 1: Record current workflow findings and add a stabilization plan.
- [x] Phase 2: Patch CG lipid table generation for residual CG-CG fitting and runtime-aligned angular sampling.
- [x] Phase 3: Patch runtime CG lipid pair/SC evaluators for minimum-image PBC and explicit angular signs.
- [x] Phase 4: Regenerate/check the focused table artifacts and run syntax/build/smoke verification.
- [x] Phase 5: Document review results and residual physical risks.

### Known Errors / Blockers
- Full bilayer stability is a physical validation problem and may require a longer workflow run than the local smoke tests can finish in one turn.
- Dynamic lipid orientation/rotational DOFs are intentionally deferred; this patch stabilizes the current static-direction model rather than adding a new integrator state.

### Review
- Fixed table generation:
  - CG-CG fitting now uses `relax_steps=0`;
  - CG-CG angular samples use the runtime pair-axis convention `Ang1(-n1 dot n12)` and `Ang2(n2 dot n12)`;
  - `cg_lipid_pair` stores a zero isotropic radial channel and a bounded directional residual over `7-14 Å`;
  - generated tables record schema, angle convention, radial mode, fit range, residual cap, knot spacing, and cutoff attrs.
- Fixed runtime evaluation:
  - `cg_lipid_pair` and `cg_lipid_sc` now use CG-specific evaluators with minimum-image PBC, explicit angular signs, and finite-range tapering;
  - NPT box updates now reach CG lipid potential nodes through the `PotentialNode` virtual box updater;
  - stage injection copies box dimensions and CG pair spline-range attrs into CG lipid nodes.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py` passed;
  - `cmake --build obj --target upside` passed;
  - focused CG-CG fit on the test DOPC reference produced zero radial knots, finite angular knots (`max |V_angular knot| ~= 92.46 E_up`), and `cutoff_ang=14.0`.

## 2026-04-30 Hybrid AA-Direct Forensic Regression Fix

### Project Goal
- Restore the pre-regression hybrid physics while keeping the AA-direct architecture (no martinize runtime path), with stable stage-6 rigid behavior and valid stage-7 mode-2 extraction.

### Architecture & Key Decisions
- Keep `hybrid_bb_map/bb_atom_index=-1` as a valid AA-direct mapping sentinel; drive BB sites from mapped N/CA/C/O carriers instead of aliasing CA coordinates.
- Keep SC-env and BB-env interactions active; fix force projection and site construction rather than damping/disable workarounds.
- Enforce pre-production rigidity by geometric rigid projection (best-fit rigid transform each step), not only force/momentum filtering.
- Restore martinize-equivalent charged termini BB typing (`Qd`/`Qa`) in extracted backbone mapping and preserve explicit dry-MARTINI-to-Upside unit conversion contract.

### Execution Phases
- [x] Phase 1: Audit and remove CA-alias/fallback regressions in C++ hybrid runtime.
- [x] Phase 2: Rewire BB virtual-site position/force projection for AA carriers without proxy-coordinate overwrite.
- [x] Phase 3: Restore BB mapping parity details (termini typing and BB charge assignment) in Python prep.
- [x] Phase 4: Verify compile, mode-2 extraction, rigid stage-6 geometry behavior, and stage-7 release behavior.

### Known Errors / Blockers
- Reduced-horizon workflow smoke runs can still show large MARTINI energies because minimization/equilibration are intentionally truncated; use full schedule for physical stability assessment.

### Review
- Fixed C++ regression points:
  - removed runtime `bb<-ca` aliasing path and preserved sentinel semantics for missing BB proxy indices;
  - prevented carrier CA atoms from being treated as fixed proxies in production constraint routing;
  - enabled COM-mapped BB site evaluation for AA carriers and conservative BB force projection back to N/CA/C/O carriers (including CA contribution);
  - stopped proxy refresh paths from overwriting physical AA carrier coordinates in AA-direct mode.
- Fixed rigid-body enforcement:
  - rigid groups now apply a weighted best-fit rigid transform each step and project atom coordinates onto that rigid manifold, while preserving net translation/rotation.
- Fixed Python mapping parity:
  - BB mapping now applies charged termini fragment overrides (`Qd`/`Qa`);
  - AA-backbone MARTINI BB charges are assigned from BB type (`Qd/Qa/...`) instead of hardcoded zero.
- Verification summary:
  - `cmake --build obj --target upside` passed;
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py` passed;
  - mode-2 extraction from stage-7 artifact now succeeds with valid AA-direct carrier fallback;
  - short stage-6 run keeps protein pair distances rigid (mean pair-distance drift ~`3e-05 Å` on sampled set);
  - short stage-7 run shows rigidity released (pair-distance drift becomes non-zero);
  - reduced full workflow run completes through stage `7.0` and mode-2 VTF extraction without backbone-map errors.

## 2026-04-28 Hybrid MARTINI Workflow Refactor

### Project Goal
- Decouple `example/16.MARTINI/run_sim_1rkl.sh` orchestration from low-level MARTINI constants, stabilization heuristics, and HDF5 mutation details.
- Make `py/martini_prepare_system.py` own the typed workflow defaults and stage metadata writes through argparse.
- Preserve the current hybrid MARTINI physics: SC-env and BB-env interface interactions must remain active in production.

### Architecture & Key Decisions
- Rewrite `run_sim_1rkl.sh` as a thin launcher that only bootstraps the environment and forwards high-level workflow controls.
- Add a `run-hybrid-workflow` Python command that owns stage prep, stage execution, continuation naming, HDF5 stage mutation, SC library validation/build, and VTF extraction dispatch.
- Move the named force caps, relaxation ramps, backbone-frame heuristics, NPT internals, Ewald alpha, and unit constants into Python argparse defaults.
- Remove dead shell-only low-level variables rather than preserving inert compatibility surfaces.
- Keep production non-protein hard-sphere replacement disabled by default and do not add any switch that disables SC-env or BB-env hybrid interface interactions.

### Execution Phases
- [x] Phase 1: Add Python workflow config/orchestration and argparse defaults.
- [x] Phase 2: Replace the base bash workflow with a high-level launcher.
- [x] Phase 3: Update associated wrappers to expose only high-level defaults.
- [x] Phase 4: Verify syntax, parser surface, code search, and reduced workflow behavior where feasible.
- [x] Phase 5: Record review results and remaining blockers.

### Known Errors / Blockers
- None identified at task start.

### Review
- Implemented:
  - `py/martini_prepare_system.py run-hybrid-workflow` now owns hybrid packing, stage-file conversion, HDF5 metadata injection, continuation naming, stage execution, SC-library validation/build, and VTF extraction.
  - `run_sim_1rkl.sh` is now a high-level launcher with orchestration controls only; low-level force caps, ramps, NPT/Ewald internals, and unit constants are argparse defaults in Python.
  - `run_sim_1rkl_outlipid.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_outlipid.sh` are thin wrappers that set high-level system/placement/run defaults.
  - production stage files now write migrated hybrid-control attrs from Python and explicitly keep `production_nonprotein_hard_sphere=0`.
- Physics integrity audit:
  - reduced full workflow generated a production stage containing both `martini_potential` and `martini_sc_table_1body`;
  - production `hybrid_control` had `activation_stage=production`, `sc_env_lj_force_cap=25`, `sc_env_coul_force_cap=25`, `sc_env_relax_steps=150`, `sc_env_backbone_hold_steps=200`, `sc_env_po4_z_hold_steps=150`, `sc_env_po4_z_clamp_enabled=1`, and `nonprotein_hs_force_cap=100`;
  - production environment membership retained `4025` non-protein targets.
- Verification:
  - `bash -n` passed for all four workflow shell entrypoints.
  - `.venv/bin/python -m py_compile py/martini_prepare_system.py` passed.
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --help` exposed the migrated argparse surface.
  - Search found no migrated low-level variable names in the bash wrappers or `run.py`.
  - Reduced continuation smoke wrote `/tmp/hybrid_refactor_continue/checkpoints/1rkl.stage_7.1.up` and VTF.
  - Reduced fresh workflow wrote `/tmp/hybrid_refactor_full/checkpoints/1rkl.stage_7.0.up` and VTF after completing stages `6.0` through `7.0`.

## 2026-04-28 Hybrid MARTINI First-Frame Geometry Regression

### Project Goal
- Restore stable 1RKL initial geometry after the workflow refactor; the first VTF frame must retain the proper helical backbone while preserving all active SC-env and BB-env hybrid interface interactions.

### Architecture & Key Decisions
- Treat the regression as a stage-file geometry/mapping issue unless runtime evidence proves otherwise.
- Keep the fix in `py/martini_prepare_system.py`; do not reintroduce low-level shell controls or wrapper compatibility paths.
- Preserve direct stage-runtime BB mapping and production SC/BB environment potentials.
- Compact appended AA reference atoms to the actual referenced N/CA/C/O atoms instead of preserving sparse raw AA PDB index gaps as runtime particles.
- Stage-7 CB placement must consume `hybrid_bb_map/atom_indices` after runtime remapping, not recompute runtime indices from raw `reference_atom_indices`.

### Execution Phases
- [x] Phase 1: Inspect generated stage-7 mapping metadata and identify the sparse-reference atom expansion.
- [x] Phase 2: Patch hybrid mapping injection to append only compact referenced AA atoms.
- [x] Phase 3: Verify reduced 1RKL workflow output and first-frame mapping integrity.
- [x] Phase 4: Record results and remaining risk.

### Known Errors / Blockers
- Reduced one-step workflows remain invalid as stability benchmarks because intentionally collapsed 6.x relaxation can leave enormous MARTINI energies before stage 7. Use them only for wiring and geometry checks.

### Review
- Fixed:
  - `py/martini_prepare_system.py` now appends only the compact set of referenced AA backbone atoms and remaps runtime BB carriers through that compact lookup.
  - `py/martini_prepare_system_lib.py` stage-7 CB placement now consumes the already-remapped `hybrid_bb_map/atom_indices` directly.
- Verification:
  - reduced fresh run completed through `/tmp/hybrid_refactor_fix_full2/checkpoints/1rkl.stage_7.0.up` and `/tmp/hybrid_refactor_fix_full2/1rkl.stage_7.0.vtf`;
  - production stage has `n_atom=4214`, `reference_index_count=124`, runtime BB indices `4090..4213`, and all runtime BB indices in bounds;
  - appended AA role counts are exactly `31` each for `N`, `CA`, `C`, and `O`;
  - production HDF5 still contains both `martini_potential` and `martini_sc_table_1body`, with `activation_stage=production` and `production_nonprotein_hard_sphere=0`;
  - stage-7 input backbone distances are normal before dynamics: N-CA mean `1.471 A`, CA-C mean `1.539 A`, C-O mean `1.234 A`, CA-CA next mean `3.906 A`.

## 2026-04-28 Hybrid MARTINI Stage-7 Explosion Regression

### Project Goal
- Stop the current stage-7 1RKL production blow-up while keeping SC-env and BB-env interface potentials active.

### Architecture & Key Decisions
- Treat the user-reported trajectory as authoritative: the workflow is still wrong if Rg jumps from `12.6 A` to `17862.9 A` by step `50`.
- Compare the runtime force carrier model against the last stable direct-Upside hybrid semantics before making more cleanup edits.
- Do not disable `martini_potential`, `martini_sc_table_1body`, SC-env, or BB-env interactions.
- Do not add debug bypass flags.

### Execution Phases
- [x] Phase 1: Re-open regression notes and capture the user correction.
- [x] Phase 2: Inspect stage-7 force carrier and pair semantics against a known stable artifact.
- [x] Phase 3: Patch the wrong carrier/integration path directly.
- [x] Phase 4: Replay a short production window and verify Rg and energies do not explode.

### Known Errors / Blockers
- User-reported production run explodes by step `50` after the compact mapping fix:
  - step `0`: `Rg 12.6 A`, total `-12776.13`;
  - step `50`: `Rg 17862.9 A`, total `1.111e+12`;
  - step `100`: `Rg 34374319.3 A`, total `2.324e+20`.

### Review
- Fixed two runtime issues:
  - `src/martini.cpp` now uses the remapped `hybrid_bb_map/atom_indices` as the runtime N/CA/C/O carrier indices for backbone O refresh instead of recomputing `reference_index_offset + raw_reference_index`;
  - `martini_sc_table_1body` and the legacy SC-table node now apply the configured SC-env force cap to table point/vector gradients whenever the cap is positive, keeping SC-env active but preventing launch-force startup events.
- Verification:
  - rebuilt `obj/upside` successfully;
  - replayed the saved failing handoff for `1000` steps from `/tmp/hybrid_stage7_replay_cap1000.up`;
  - replay stayed finite with CA Rg `12.5-12.7 A` through the full replay;
  - output HDF5 is all finite, with CA Rg `12.759 A` at frame `0` and `12.723 A` at the last saved frame;
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed;
  - `bash -n` passed for all four MARTINI workflow entrypoint scripts.

## 2026-04-28 Hybrid MARTINI Production Restart Parity

### Project Goal
- Fix stage-7 production restart behavior so split runs (`7.0 -> 7.1 -> ...`) preserve the same protein stability as an uninterrupted production run.

### Architecture & Key Decisions
- Restart must preserve physical state, not just positions: momenta/velocities, box, hybrid startup counters, and dynamic hold/ramp state must be equivalent.
- Continue using active SC-env and BB-env interactions; do not add restart-only debug bypasses.
- Compare split-vs-continuous production using the saved stage-7 handoff before declaring the bug fixed.

### Execution Phases
- [x] Phase 1: Document the restart regression and inspect continuation workflow.
- [x] Phase 2: Identify missing restart state by comparing uninterrupted and split runs.
- [x] Phase 3: Patch continuation state transfer.
- [x] Phase 4: Verify split-vs-continuous stability and log results.

### Known Errors / Blockers
- No active blocker remains for production restart state transfer.

### Review
- Fixed:
  - production MD runs now record momentum so stage files contain restartable `output/mom`;
  - handoff now copies the last saved `output/mom` to `input/mom` and marks it restart-valid;
  - continuation now passes `--restart-using-momentum` only when copied momentum is valid;
  - production continuation now preserves the hybrid SC-env transition counter with `hybrid_control.sc_env_transition_step_start`;
  - `src/martini.cpp` initializes `sc_env_transition_step` from that restart attr instead of always restarting the ramp/hold from zero.
- Verification:
  - direct continuous `1000`-step replay from the saved stage-7 handoff stayed stable with Rg around `12.7-12.9 A`;
  - split replay (`500 + 500`) stayed stable after restart with Rg around `12.8-13.0 A`;
  - actual workflow continuation wrote `/tmp/restart_workflow_check/checkpoints/1rkl.stage_7.1.up`, invoked `--record-momentum --restart-using-momentum`, and kept Rg finite around `12.9 A`;
  - generated continuation file has `input/mom.restart_valid=1`, `output/mom`, and `sc_env_transition_step_start=450`;
  - `cmake --build obj --target upside`, Python compile checks, and shell syntax checks passed.

## 2026-04-29 Hybrid MARTINI Restart Instability Still Reproduces

### Project Goal
- Resolve the still-observed production restart instability when starting from a previous production stage.

### Architecture & Key Decisions
- Treat the prior restart fix as incomplete.
- Audit whether continuation is using the true final integrator state or merely the last logged frame.
- If exact restart state is unavailable, fail clearly instead of silently rethermalizing or restarting from stale logged state.
- Production-to-production continuation must preserve saved coordinates and momenta exactly; hybrid reference carriers must not be rebuilt during restart handoff.
- MD output must include the true final state so a later continuation does not restart from the last pre-final logging interval.
- Preserve active SC-env and BB-env interactions.

### Execution Phases
- [x] Phase 1: Record unresolved restart instability correction.
- [x] Phase 2: Audit restart source state actually available in current production files.
- [x] Phase 3: Fix continuation to use true final state or fail loudly when unavailable.
- [x] Phase 4: Verify restart from a real production stage generated by the workflow.

### Known Errors / Blockers
- Existing production files created without `--record-momentum` cannot be exact restart sources; continuation now fails before MD and tells the user to regenerate the previous production stage with the current workflow.
- Existing production files created before the validated final-state marker also cannot be accepted as exact restart sources, even if they contain `output/mom`.

### Review
- Fixed:
  - production `7.x -> 7.x` continuation now uses `production_restart` handoff, which copies saved final coordinates and momentum exactly and does not refresh hybrid carriers;
  - continuation now requires both a validated final-state marker on `output/mom` and `input/mom.restart_valid=1` after handoff;
  - `upside` now appends a final MD output sample after the integration loop so the next stage restarts from the true final state;
  - workflow MD stages mark `output/mom.restart_final_state_valid=1` only after verifying the final output time matches the requested duration;
  - hybrid transition carry-forward now advances by the full final output time instead of the last periodic pre-final frame.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed;
  - `bash -n` passed for all four MARTINI workflow shell entrypoints;
  - `cmake --build obj --target upside` passed with only existing warnings;
  - a patched 500-step production source wrote final `output/time=1.0` and matching `output/mom`;
  - workflow continuation from that source used `--restart-using-momentum`, kept Rg finite from `12.6 A` to `12.3 A`, and produced finite positions/momenta;
  - continuation input exactly matched source final state: position max diff `0.0`, momentum max diff `0.0`, `input/mom.restart_valid=1`, `sc_env_transition_step_start=500`;
  - workflow-generated `stage_7.1` was restarted into `stage_7.2`, used `--restart-using-momentum`, stayed stable from `Rg 12.2 A` to `12.1 A`, and produced `sc_env_transition_step_start=1000`;
  - `stage_7.2` input exactly matched `stage_7.1` final state: position max diff `0.0`, momentum max diff `0.0`;
  - continuation from a source with `output/mom` but no validated final-state marker failed before MD with the expected final restart-state error;
  - continuation from an older source without `output/mom` failed before MD with the expected restart-valid momentum error.

## 2026-04-24 dry-MARTINI Runtime Acceleration Import

### Project Goal
- Import the dry-MARTINI runtime acceleration logic from `/Users/yinhan/Documents/upside2-md_temp/src/martini.cpp` into this repo's `src/martini.cpp`.
- Keep the scope limited to performance paths used by the active dry-MARTINI runtime, not unrelated feature changes from the temp checkout.

### Architecture & Key Decisions
- Only transplant acceleration mechanisms that directly reduce hot-loop work in the current runtime:
  - cached Verlet-style pairlists for `MartiniPotential`,
  - cached per-pair parameter metadata so the main nonbonded loop avoids repeated coefficient unpacking / spline-map lookup,
  - cached active contact lists for `martini_sc_table_1body`, which is the active production SC/environment node in this checkout.
- Do not import unrelated temp-file changes such as fix-rigid extensions, stage-parameter cleanup, or removal of legacy compatibility nodes unless they are required for the acceleration paths to compile.
- Preserve backward compatibility in `src/martini.cpp`:
  - keep existing node names and interfaces intact,
  - keep legacy `martini_sc_table_potential` available.

### Execution Phases
- [x] Confirm the exact acceleration deltas in the temp file and map them onto the active runtime surfaces in this checkout.
- [x] Patch `src/martini.cpp` with the selected acceleration paths only.
- [x] Build and verify the modified MARTINI runtime.
- [x] Document the imported optimizations and verification results.

### Known Errors / Blockers
- The checked-in `obj/` build cache is stale from the old moved repo path (`/Users/yinhan/Documents/upside2-md-martini`), so verification in this pass used a fresh out-of-tree build directory under `/tmp` instead of the existing `obj/`.

### Review
- Confirmed the temp checkout acceleration scope:
  - `MartiniPotential` adds a cached Verlet-style pairlist plus cached per-pair parameter metadata/pointer indirection;
  - the active SC/environment runtime in this checkout is `martini_sc_table_1body`, and the temp checkout accelerates that path with a cached active-contact list.
- Implemented in `src/martini.cpp`:
  - `MartiniPotential` now compacts pair coefficients into a unique parameter table, caches direct spline pointers per unique parameter row, and rebuilds an active pairlist only when atoms move beyond half the configured skin or the box changes;
  - `MartiniScTableOneBody` now caches row/environment contacts within `cutoff + cache_buffer` and reuses them across both value and derivative passes until the cache is invalidated by motion or box changes;
  - the NPT box updater now also updates `martini_sc_table_1body`, so the active-contact cache remains valid under box scaling for the live production node.
- Verification:
  - configured a fresh build tree with `cmake -S src -B /tmp/upside2-md-build-20260424`;
  - built `upside`, `upside_calculation`, and `upside_engine` successfully with `cmake --build /tmp/upside2-md-build-20260424 --target upside upside_calculation upside_engine`.
- Residual warnings:
  - build still reports a pre-existing `%p` formatting warning in `martini_masses::get_mass(...)` and unrelated existing warnings in other translation units; no new build errors remain from the acceleration import.

## 2026-04-14 MARTINI Workflow Python Runtime Repair

### Project Goal
- Make `example/16.MARTINI/run_sim_1rkl.sh` use the repository `.venv` reliably so `py/martini_prepare_system_lib.py` can import `h5py`.
- Repair `install_python_env.sh` so a repo-local `.venv` still works after the repository has been moved on disk.

### Architecture & Key Decisions
- Treat the reported `ModuleNotFoundError: h5py` as an interpreter-selection bug, not a missing dependency bug:
  - the current `.venv/bin/python` imports `h5py` successfully;
  - the failure comes from `run_sim_1rkl.sh` resolving `python3` from Homebrew `python3.14` after sourcing a stale `.venv/bin/activate`.
- Fix the workflow wrapper directly by prepending the repo-local `.venv/bin` to `PATH` instead of depending on `activate` side effects.
- Fix the installer in place by repairing stale `.venv/bin/activate*` files and stale console-script shebangs under `.venv/bin` to the current repo path, avoiding an unnecessary full environment rebuild.

### Execution Phases
- [x] Confirm the interpreter/path mismatch and inspect the current `.venv` metadata.
- [x] Patch `example/16.MARTINI/run_sim_1rkl.sh` and `install_python_env.sh` with the minimal robust fix.
- [x] Verify the repaired runtime path and document the result.

### Known Errors / Blockers
- None so far. The local `.venv` already contains `h5py`, so verification does not depend on a fresh network install.

### Review
- Root cause verified:
  - `.venv/bin/python` imports `h5py`;
  - Homebrew `python3.14` does not;
  - the old `.venv/bin/activate` still pointed at `/Users/yinhan/Documents/upside2-md-martini/.venv`, so `run_sim_1rkl.sh` was picking the wrong interpreter.
- Implemented:
  - `example/16.MARTINI/run_sim_1rkl.sh` now prepends the repo-local `.venv/bin` to `PATH` instead of sourcing `activate`;
  - `install_python_env.sh` now repairs stale `activate*` files and stale `.venv/bin` Python-script shebangs in place.
- Verification:
  - `bash -n install_python_env.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -lc 'PROJECT_ROOT="$PWD"; source "$PROJECT_ROOT/source.sh"; export PATH="$PROJECT_ROOT/.venv/bin:$PATH"; python3 -c "import sys, h5py; print(sys.executable); print(h5py.__version__)"'`
  - `env PIP_NO_INDEX=1 bash install_python_env.sh`
  - `bash -lc 'source .venv/bin/activate && python3 -c "import sys, h5py; print(sys.executable); print(h5py.__version__)"'`
  - `bash -lc 'source .venv/bin/activate && pip --version'`

## 2026-04-14 Python Environment Installer Portability

### Project Goal
- Make `install_python_env.sh` succeed on both Linux and local macOS hosts while preserving the existing Python-package intent for UPSIDE2.

### Architecture & Key Decisions
- Keep Python `3.11` as the required interpreter version for the repository environment.
- Interpreter selection will stay configuration-driven:
  - honor `PYTHON_BIN` when provided;
  - otherwise discover a compatible Python `3.11` interpreter from common Linux/macOS locations instead of assuming a single command name.
- macOS setup will export Homebrew-backed dependency hints for packages that may need native HDF5 / compression libraries during pip installation.
- Core scientific packages remain required.
- Optional analysis packages remain best-effort so a platform-specific failure in `pyhdx`-related extras does not abort the whole environment bootstrap.

### Execution Phases
- [x] Audit the current installer and local Darwin environment.
- [x] Patch `install_python_env.sh` for cross-platform interpreter and package installation behavior.
- [x] Verify shell syntax and installer logic, then document the result.

### Known Errors / Blockers
- No active code blocker remains for this installer task.
- Verification completed in the current Darwin arm64 shell with `PIP_NO_INDEX=1`, which exercised interpreter discovery, macOS build-hint setup, venv reuse, and the post-install import checks without requiring fresh downloads.

## 2026-04-01 Stage-7 SC Table Integration

### Project Goal
- Consume the completed external `/Users/yinhan/Documents/SC-training` training results by building a native-unit `martini.h5`.
- Replace the legacy stage-7 sidechain back-mapping / rotamer production path in `example/16.MARTINI/run_sim_1rkl.sh` with a direct dry-MARTINI sidechain table interaction.
- Keep the new sidechain force field inactive for stages `6.0` through `6.6` and active only in production stage `7.0`.

### Architecture & Key Decisions
- `martini.h5` will be built from the assembled `sc_table.json` output and stored as a native-unit library (`nm`, `kJ/mol`) with a simple table schema:
  - residue order,
  - target dry-MARTINI type order,
  - uniform radial grid,
  - radial energy table.
- Simulation-side conversion into Upside units remains parameterized:
  - the stage-7 runtime node will read explicit `energy_conversion_kj_per_eup` and `length_conversion_angstrom_per_nm` attrs;
  - no numeric unit conversion will be baked into training outputs.
- Revised decision after user correction on 2026-04-02:
  - restore the reduced simulation-mass path (`ff_mass / 12`) in preparation and stage-file augmentation;
  - keep the weighted BB mapping revert in place;
  - retain the scientific-notation logging improvement only as a readability fix.
- Revised decision after user correction:
  - the stage-7 dry-MARTINI coupling must be integrated directly with Upside carriers, not evaluated on parallel protein MARTINI proxy coordinates and projected afterward.
- Sidechain/environment coupling will remain a separate stage-7 potential node, but it must be directly Upside-carried:
  - use deterministic Upside sidechain-bead coordinates from `affine_alignment -> placement_fixed_point_only_CB`;
  - act only between selected protein residues and non-protein dry-MARTINI particles;
  - add equal-and-opposite force to the dry particle and the Upside sidechain-bead coordinate so the backbone receives the feedback through the coord-node derivative chain;
  - all legacy protein MARTINI `SC` proxy nonbonded interactions in `martini_potential` must be disabled in stage 7.
- Backbone/environment coupling must also be directly Upside-carried:
  - treat the Upside `CA` carrier as the dry-MARTINI protein `BB` interaction site;
  - keep dry-MARTINI `BB` to surrounding non-protein dry particles active, but apply the resulting force directly to `CA` and the dry particle rather than to a proxy `BB` bead followed by projection;
  - revised after regression report: the old weighted `BB -> N/CA/C/O` active map has been restored for now; the attempted CA-only active BB remap has been reverted pending a deeper force-path audit;
  - revised after the 2026-04-02 production-drift audit:
    - direct `BB -> CA` and direct `CB` stage-7 forces must still honor the existing startup protein-feedback ramp (`sc_env_backbone_hold_steps`) instead of feeding full force back onto the protein from step `0`;
    - stage-7 must not evaluate any protein-internal MARTINI proxy-proxy terms, because those legacy bonded/nonbonded terms only add bookkeeping energy once the protein is carried directly by Upside;
  - do not evaluate protein internal SC-backbone dry-MARTINI interactions, because those are already handled on the Upside side.
- The workflow must remove the legacy production SC path together with its dead control surface:
  - no stage-7 `rotamer`, `placement_fixed_scalar`, or `placement_fixed_point_vector_only`;
  - no stage-7 `hybrid_sc_map` requirement;
  - no stage-7 SC-relaxation / proxy-control attributes that only served the old probabilistic proxy path.
- Revised decision after the 2026-04-02 cleanup pass:
  - remove the dead prep-only `hybrid_sc_map` export/validation path as well; the active workflow only needs `hybrid_control`, `hybrid_bb_map`, and `hybrid_env_topology`;
  - `src/box.cpp` and `src/box.h` were audited during this pass and left unchanged because they no longer contain the retired stage-7 proxy/rotamer subsystem, only active box/NPT infrastructure.
- `ALA` and `GLY` remain excluded from the new SC table because the trained library covers the current `18` non-empty canonical sidechain types; those residues still participate through the existing backbone/environment MARTINI path only.

### Execution Phases
- [x] Phase A: Build `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- [x] Phase B: Implement a new stage-7 SC table runtime node and wire its box-scaling behavior.
- [x] Phase C: Replace the production-side script injection in `run_sim_1rkl.sh` with stage-7 CB + SC-table setup and remove the legacy SC rotamer path.
- [x] Phase D: Update the production helper workflow to match the new stage-7 path.
- [x] Phase E: Remove the remaining parallel protein-MARTINI force path so stage 7 uses direct Upside-carried BB/SC interactions only.
- [x] Phase F: Verify `martini.h5` generation, build integrity, stage-7 HDF5 injection structure, and workflow/runtime behavior for the direct-Upside design.

### Known Errors / Blockers
- No code-path blocker remains for stage-7 SC-table injection in this environment:
  - the repaired repository `.venv` provides both `h5py` and `tables`;
  - the stage-7 injector now falls back to the protein ITP when fresh prepared files omit `/input/sequence`.
- Reduced-step full-workflow runs are still not a stability benchmark for the hybrid physics:
  - the intentionally shortened `6.4-6.6` relaxation settings used for fast verification can drive the pre-production MARTINI state to enormous energies before stage 7;
  - that does not block the stage-7 integration itself, but any physical assessment still needs the intended relaxation horizon / benchmark settings.

### Review
- 2026-04-02 legacy-code removal verification:
  - physically removed the dead probabilistic SC / rotamer / placement subsystem from `src/martini.cpp` instead of leaving it compiled behind runtime gates;
  - removed the dead prep-only `hybrid_sc_map` generation and validation path from `example/16.MARTINI/prepare_system_lib.py`, `example/16.MARTINI/prepare_system.py`, and `example/16.MARTINI/validate_hybrid_mapping.py`;
  - audited `src/box.cpp` and `src/box.h` and left them intact because the remaining code there is still-active barostat/box plumbing rather than disabled hybrid legacy logic;
  - verification:
    - `cmake --build obj` passed,
    - `python -m py_compile example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/prepare_system.py example/16.MARTINI/validate_hybrid_mapping.py` passed under `.venv`,
    - `bash -n example/16.MARTINI/run_sim_1rkl.sh` and `bash -n example/16.MARTINI/test_prod_run_sim_1rkl.sh` passed,
    - a fresh prep run in `/tmp/hybrid_mapping_bb_only` produced a mapping file whose `/input` group contains only `hybrid_bb_map`, `hybrid_control`, and `hybrid_env_topology`,
    - a shortened active-workflow run in `/tmp/legacy_cleanup_short` accepted that mapping schema and completed the full current `run_sim_1rkl.sh` ladder through stage `7.0`.
- 2026-04-02 production-drift fix verification:
  - restored the startup protein-feedback ramp on the direct stage-7 force paths in `src/martini.cpp`:
    - direct `BB`/environment MARTINI forces now scale only the protein-side feedback by `compute_sc_backbone_feedback_mix(...)`,
    - `martini_sc_table_potential` now does the same for `CB` feedback while leaving the environment-side force unchanged;
  - removed the remaining active protein-internal MARTINI proxy-proxy path during stage 7 by making `allow_protein_pair_by_rule(...)` reject all protein-proxy pairs in the active hybrid runtime;
  - targeted replay verification from the same saved production handoff (`example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`) now stays negative through the user’s failure window:
    - archived current run at step `2050`: `martini_potential 26387.78`, `total 26216.96`;
    - patched replay at step `2050`: `martini_potential -25985.45`, `total -26147.28`.
- Implemented `example/16.MARTINI/build_sc_martini_h5.py` and generated `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- Implemented `example/16.MARTINI/inject_sc_table_stage7.py` so stage 7 now:
  - preserves the existing Upside backbone-node augmentation,
  - deletes the legacy SC production path,
  - injects deterministic CB coordinates and the new `martini_sc_table_potential`.
- Updated both `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` to use the new stage-7-only SC table path and auto-build the library if only the assembled training JSON is available.
- Implemented the new `martini_sc_table_potential` runtime in `src/martini.cpp`:
  - CB-to-env radial table evaluation,
  - equal-and-opposite force accumulation,
  - additive backbone feedback through the CB coord node,
  - participation in NPT box scaling.
- Removed the remaining parallel protein proxy path from the active stage-7 runtime/workflow:
  - `martini_potential` now skips protein `SC` proxy MARTINI pairs entirely during active hybrid production;
  - protein `BB` proxy MARTINI pairs now use the Upside `CA` carrier coordinate as the interaction site and accumulate force directly onto `CA` plus the dry-MARTINI partner;
  - the active `hybrid_bb_map` is now CA-only in both exported mapping files and injected runtime stage files (`atom_mask=[0,1,0,0]`, `weights=[0,1,0,0]`), and runtime parsing now rejects any remaining active `N/C/O` carrier rows;
  - stage-file injection no longer copies `hybrid_sc_map` into runtime stages, the old shell-level SC control helper is gone from the active workflow path, and the legacy stage-7 rotamer augmenter now hard-fails if called.
- Verification completed in this environment:
  - `python3 -m py_compile` passed for the new Python scripts;
  - `bash -n` passed for both workflow shell scripts;
  - `cmake --build obj` passed;
  - `h5ls -r parameters/ff_2.1/martini.h5` confirmed the generated library shape.
- Regression revert verification completed after the user reported a broken minimization stage:
  - reverted the attempted native-mass switch and the CA-only active BB remap;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run again reaches stage `6.0` with a normal minimization startup:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - this replaces the broken regression signature reported by the user:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/30400945802214465601536.00/...`
- Native-mass-only restoration verification completed on 2026-04-02:
  - restored the native mass path while leaving the weighted 4-carrier BB map in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run still starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - fresh `1rkl.stage_6.0.up` in `/tmp/martini_mass_restore` again shows native dry-MARTINI mass values in `/input/mass` (`72.0` and `45.0` for the physical dry particles).
- Reduced-mass restoration verification completed after the user reported box blow-up with original masses:
  - restored reduced simulation masses (`ff_mass / 12`) in preparation and stage-file augmentation while leaving the scientific-notation logging change in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run in `/tmp/martini_mass_reduced` again starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - direct HDF5 inspection of `1rkl.stage_6.0.up` confirms the reduced mass values are back in `/input/mass`:
    - `4082` particles at `6.0`,
    - `8` particles at `3.75`.
- Logging/readability follow-up on 2026-04-02:
  - the huge `martini_potential` / `total` values seen in current stage-7 logs are real runtime values, not a `printf` type mismatch or summation bug;
  - `src/main.cpp` now formats large-magnitude energies in scientific notation for readability only.
- Additional end-to-end verification completed after the `.venv` repair:
  - the focused production helper `example/16.MARTINI/test_prod_run_sim_1rkl.sh` now runs end to end from a real `6.6 -> 7.0` handoff, injects the new SC table successfully, and completes a short stage-7 production replay;
  - direct HDF5 inspection confirmed stage gating:
    - `6.6` contains no `martini_sc_table_potential`, no `placement_fixed_point_only_CB`, no `hybrid_sc_map`, and keeps hybrid activation disabled;
    - `7.0.prepared` contains `martini_sc_table_potential` and `placement_fixed_point_only_CB`, remains free of legacy `rotamer` nodes, and activates the hybrid path only for production;
  - the real top-level workflow `example/16.MARTINI/run_sim_1rkl.sh` completes end to end under shortened verification settings with the direct-Upside stage-7 path.

## 2026-03-31 SC-training Workflow Addendum

### Project Goal
- Establish a new `SC-training/` workflow that derives dry-MARTINI sidechain-type to dry-MARTINI particle-type training data and assembled tables independently of the legacy SC back-mapping runtime path.
- Make the workflow runnable both locally and on Slurm.
- Use `example/16.MARTINI/run_sim_1rkl.sh` as the canonical benchmark workflow target once the trained tables are consumed by the hybrid runtime.

### Architecture & Key Decisions
- The new workflow will live under a dedicated top-level folder `SC-training/`.
- The workflow will treat the current SC back-mapping / rotamer-placement production path as legacy design context rather than a dependency.
- The replacement SC/dry-MARTINI interaction must be a two-way force term: each SC/dry pair contributes equal-and-opposite force to the dry particle and to the effective protein-side SC site, and the protein-side contribution must be projected back onto the Upside carriers and added to the existing backbone/sidechain forces before integration.
- Existing backbone/environment MARTINI coupling remains in scope and separate from the new SC table:
  - the effective protein backbone dry-MARTINI `BB`/`CA` anchor should still interact with surrounding non-protein dry-MARTINI particles;
  - the new SC table should not reintroduce generic protein SC-backbone MARTINI interactions, because protein internal SC-backbone coupling is already handled on the Upside side.
- Residue-sidechain dry-MARTINI bead definitions will be sourced from `example/16.MARTINI/martinize.py` using the `martini22` forcefield definitions already used by the example workflow.
- Dry-MARTINI pair parameters will be sourced from `example/16.MARTINI/ff_dry/dry_martini_v2.1*.itp`.
- The unit contract must be explicit and parameterized:
  - `SC-training/` artifacts remain in native dry-MARTINI units (`nm`, `kJ/mol`, `e`, dielectric `15`);
  - simulation-side conversion into Upside units must come from explicit parameters / attrs, not baked conversion numbers in training outputs or MARTINI runtime code;
  - the current MARTINI workflow expects `UPSIDE_MARTINI_ENERGY_CONVERSION` and `UPSIDE_MARTINI_LENGTH_CONVERSION` to be provided explicitly when generating simulation inputs.
- The initial training workflow will emphasize reproducible task generation, batch execution, and result assembly; any simplifying assumptions used for the first pass must be explicit in workflow outputs and README docs.
- Slurm support should follow the proven ConDiv workflow shape: explicit path resolution, environment setup, generated manifests/specs, array-friendly execution, and a submission wrapper rather than an ad hoc one-off `sbatch` command.

### Execution Phases
- [x] Phase A: Inspect and reuse the relevant MARTINI residue/type definitions and Slurm workflow patterns.
- [x] Phase B: Create `SC-training/` scripts for manifest generation, per-task training execution, and final result assembly.
- [x] Phase C: Add local and Slurm entrypoints that stage and run training jobs from the new folder.
- [x] Phase D: Add minimal documentation describing assumptions, outputs, and benchmark expectations.
- [x] Phase E: Verify script integrity and basic dry-run behavior.

### Known Errors / Blockers
- Runtime consumption of the new SC training table is not yet implemented in `src/martini.cpp`, so benchmark orchestration may need to stop at workflow staging/logging until that integration exists.
- That runtime integration still needs to enforce the new mechanical requirements explicitly:
  - SC/dry interactions must remain equal-and-opposite;
  - SC feedback must be added to existing Upside carrier gradients rather than replacing them;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling must remain active while protein internal SC-backbone coupling stays excluded from the new table.
- The first-pass training workflow may need an explicit simplifying geometry assumption for residue-level sidechain reduction; that assumption must be documented, not hidden.
- `SC-training/` is now self-contained for training and Slurm staging, but the optional benchmark hook still targets `example/16.MARTINI/run_sim_1rkl.sh` and therefore still requires the full repository when benchmark execution is requested.

### Review
- Implemented a new `SC-training/` workflow with:
  - manifest generation from `example/16.MARTINI/martinize.py` (`martini22`) and `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`,
  - per-task residue/target training jobs,
  - assembled output tables,
  - local runner,
  - Slurm staging/submission wrappers,
  - benchmark staging against `example/16.MARTINI/run_sim_1rkl.sh`.
- Smoke verification succeeded in `/tmp/sc_training_smoke`:
  - manifest generation completed,
  - full local run completed,
  - assembled summary reported `684` tasks across `18` residues and `38` target particle types,
  - Slurm array, collector, and benchmark scripts were generated successfully with `--no-submit`.
- Follow-up workflow correction restored the default target set to the full dry-MARTINI atomtype list from the bundled forcefield after a temporary over-narrowing to non-ring types:
  - default table cardinality is `18 x 38 = 684` residue-target tasks;
  - manifest/task outputs now keep the full target set while retaining the explicit spherical surrounding-position samples around the sidechain center.
- Follow-up architecture review confirmed additional force-path requirements for the future runtime hookup:
  - the new SC/dry term must remain a two-way force and its protein-side feedback must be added on top of existing Upside forces;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling remains required and should stay distinct from the new SC table;
  - dry-MARTINI unit handling is now defined as native-unit training plus parameterized simulation-side conversion, not a baked training/runtime numeric contract.
- Portability follow-up verification succeeded from a standalone folder copy under `/private/tmp`:
  - `SC-training/data/` now bundles the required dry-MARTINI nonbond parameter file and martini22 sidechain residue definitions;
  - `workflow.py init-run`, `workflow.py run-local`, `run_local.sh`, and stage-only `submit_remote_round.sh` all work from the copied folder without the parent repository;
  - generated Slurm array/collector scripts now reference the copied `SC-training/` folder directly instead of requiring `source.sh` or repo-root paths for training.
- Slurm-array verification succeeded from a fresh staged run:
  - `submit_remote_round.sh` now reaches `workflow.py submit-slurm` reliably even when repo-root `source.sh` would otherwise trip over an unset `MY_PYTHON`;
  - generated training script uses `#SBATCH --array=0-683` for the current manifest and dispatches one residue-target task per array index via `run-array-task --task-id "$SLURM_ARRAY_TASK_ID"`;
  - staged round manifest count and training manifest count matched exactly (`684` / `684`), with a separate collector script for post-array assembly.
- Batch-spooled launcher verification succeeded after the user-reported Midway failure:
  - `submit_remote_round.sh` and `run_local.sh` no longer trust `dirname "$0"` as the workflow root when executed from a Slurm spool copy;
  - the wrappers now locate the real `SC-training` directory by checking candidate directories for `workflow.py`, including `SLURM_SUBMIT_DIR` and `SLURM_SUBMIT_DIR/SC-training`;
  - simulated `sbatch submit_remote_round.sh` from inside `SC-training/` and simulated `sbatch SC-training/submit_remote_round.sh` from the repo root both staged cleanly and produced matching `684`-task manifests.
- Unit-contract follow-up verification succeeded:
  - training manifests now record native-unit policy only and no longer emit baked Upside conversion numbers;
  - `prepare_system_lib.py` now requires explicit conversion env vars and writes those MARTINI unit-conversion attrs into the simulation input;
  - `src/martini.cpp` now derives Coulomb scaling from those attrs instead of a literal hardcoded runtime constant.

### Phase 1: Geometric Mapping (Defining the Training Ground)
Before you can train the splines, you must define the physical relationship between the two models.
1. **Define the Upside Virtual Bead:** For each of the 20 standard amino acids, identify the Upside definition for the virtual side-chain bead position ($\vec{r}_{SC}$) and its normal vector ($\vec{n}$) relative to the backbone.
2. **Define the MARTINI Target:** For each amino acid, map the corresponding multi-bead Dry-MARTINI representation (e.g., mapping a Tryptophan side chain to its 4 constituent MARTINI beads). 

### Phase 2: Generating the 2D Energy Landscapes
You need to calculate the "ground truth" thermodynamic interactions using the Dry-MARTINI force field. Because Dry-MARTINI is an implicit solvent model, this can be done in a vacuum.
1. **Set Up the Grid:** Place your multi-bead MARTINI side chain at the origin. Place a single target Dry-MARTINI particle (e.g., a standard $C1$ lipid bead or $P4$ polar bead) in the space around it.
2. **Scan the Space:** Systematically move the single target bead across a 2D grid defined by:
   * **Distance ($r_{12}$):** The distance from the target bead to the Upside virtual center ($\vec{r}_{SC}$).
   * **Angle ($\theta_1$):** The angle between the Upside normal vector ($\vec{n}$) and the vector connecting the two beads ($\vec{n}_{12}$).
3. **Calculate Energy:** At every grid point, calculate the sum of the standard Dry-MARTINI Lennard-Jones and Coulomb interactions between the single target bead and the multi-bead side chain.
4. **Repeat:** Do this for all 20 amino acids against all relevant Dry-MARTINI particle types.

### Phase 3: Mathematical Spline Fitting
You must compress those 2D energy landscapes into the modified, factorized Upside equation.
1. **Set Up the Equation:** Use the simplified hybrid Hamiltonian that drops the second angle term:
   $$V_{hybrid} = V_{radial}(r_{12}) + ang_1(-n_1 \cdot n_{12}) \cdot V_{angular}(r_{12})$$
2. **Initialize Splines:** Set up $V_{radial}$, $ang_1$, and $V_{angular}$ as 1D natural cubic splines with adjustable control points (knots).
3. **Optimize:** Use a least-squares fitting algorithm or gradient descent to adjust the spline knots until the output of the $V_{hybrid}$ equation closely matches the 2D energy grid you generated in Phase 2.

### Phase 4: Engine Integration (C++ Implementation)
Append the new physics into the Upside MD engine's core loop.
1. **Introduce Real Particles:** Modify the Upside state to accept independent, mass-bearing Dry-MARTINI particles with standard $x, y, z$ coordinates and momenta.
2. **Calculate the Force:** In the force-evaluation loop, evaluate $V_{hybrid}$ and calculate its spatial gradient (force) and angular derivative (torque).
3. **Push the Real Bead:** Apply the negative spatial gradient directly to the Dry-MARTINI particle's momentum.
4. **Push the Ghost Bead:** Apply the exact inverse spatial vector to the Upside virtual side chain ($\vec{r}_{SC}$). Apply the angular derivative as a torque to the virtual normal vector ($\vec{n}$).
5. **Distribute to the Backbone:** Feed the virtual force and torque into Upside's existing chain-rule matrices to distribute the physical forces onto the $N$, $C\alpha$, and $C$ backbone atoms.

### Phase 5: Validation and Thermodynamic Tuning
Run test simulations to ensure the hybrid physics do not break the model.
1. **Check Momentum Conservation:** Monitor the total momentum of the simulation box. If the sum of forces on the MARTINI beads and the Upside backbone atoms does not equal exactly zero at every time step, there is a bug in your chain-rule implementation.
2. **Tune the Scaling Factor ($\kappa$):** Introduce a global scaling multiplier to the hybrid potential: $V_{total} = \kappa(V_{hybrid})$. Run a folded protein in a box of your Dry-MARTINI particles. If the protein denatures because the MARTINI interactions are overwhelmingly strong compared to Upside's internal folding potentials, lower $\kappa$ until the thermodynamic balance is restored.

What specific biological system or environment are you aiming to simulate first once this hybrid engine is built?

## 2026-04-13 Interface-Only Hybrid Scaling Calibration

### Project Goal
- Calibrate one shared dry-MARTINI interaction scale from bilayer-only DOPC diffusion runs.
- Apply that scale only to Upside/dry-MARTINI interface interactions in active hybrid production.
- Leave bilayer-bilayer interactions unchanged in the hybrid workflow.

### Architecture & Key Decisions
- Use one shared scalar `protein_env_interface_scale` rather than separate LJ and Coulomb knobs.
- In active hybrid stage `7.0`, apply the scale only to cross-interface protein/environment interactions:
  - direct `martini_potential` protein-environment pairs,
  - `martini_sc_table_1body`,
  - legacy `martini_sc_table_potential` for backward compatibility.
- Preserve current hybrid startup logic such as the protein-feedback ramp; the new factor is an additional interface-strength control, not a replacement for startup stabilization.
- Calibrate the scalar in the bilayer-only workflow by rewriting dry-MARTINI pair coefficients before simulation:
  - `epsilon *= pair_scale`
  - `q_i *= sqrt(pair_scale)`
  - `q_j *= sqrt(pair_scale)`
- Keep bilayer calibration separate from hybrid runtime application:
  - bilayer-only runs may scale all MARTINI pairs to match target fluidity,
  - hybrid production must apply the chosen factor only to protein-environment interface terms.
- Correct the active-stage assumption from earlier notes:
  - the current production SC path in this checkout is `martini_sc_table_1body`,
  - `martini_sc_table_potential` remains only as a compatibility surface.

### Execution Phases
- [x] Phase A: Extend hybrid-control schema and runtime state with `protein_env_interface_scale`.
- [x] Phase B: Apply interface-only scaling in `src/martini.cpp` for direct pair terms and both SC-table paths.
- [x] Phase C: Expose the production control in `example/16.MARTINI/run_sim_1rkl.sh`.
- [x] Phase D: Add bilayer-only pair-scale calibration support plus a fixed sweep wrapper and cross-run report.
- [x] Phase E: Verify build/script integrity and record review notes, progress, and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the implementation must respect two different semantics:
  - bilayer calibration runs scale all MARTINI interactions inside the bilayer-only workflow,
  - hybrid production must scale only protein-environment interface terms.

### Review
- Implemented the new hybrid control attr `protein_env_interface_scale` in:
  - `src/martini.cpp`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/run_sim_1rkl.sh`
- Runtime behavior now matches the intended split:
  - active direct protein-environment `martini_potential` pairs scale only the cross-interface LJ+Coulomb interaction,
  - bilayer-bilayer and other non-interface pairs remain unchanged in the hybrid workflow,
  - `martini_sc_table_1body` and legacy `martini_sc_table_potential` both respond to the same scale.
- Added bilayer calibration support in:
  - `bilayer-lateral-diffusion/workflow.py`
  - `bilayer-lateral-diffusion/submit_interface_scale_calibration_round.sh`
  - `bilayer-lateral-diffusion/report_interface_scale_calibration.py`
- Verification completed:
  - `python3 -m py_compile py/martini_prepare_system_lib.py bilayer-lateral-diffusion/workflow.py bilayer-lateral-diffusion/report_interface_scale_calibration.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`
  - `bash -n bilayer-lateral-diffusion/run_local.sh`
  - `bash -n bilayer-lateral-diffusion/submit_interface_scale_calibration_round.sh`
  - `cmake --build obj --target upside`
  - no-submit wrapper staging for `/tmp/bilayer_interface_scale_calibration/pair_scale_0p85`
  - reduced local bilayer smoke runs at `pair_scale=1.0` and `0.85`
  - direct HDF5 inspection of staged `stage_6.0.prepared.up` files confirmed:
    - `epsilon` ratio = `0.85`
    - charge ratios = `sqrt(0.85) = 0.921954`
  - `report_interface_scale_calibration.py` executed successfully on a temporary two-run smoke tree.

## 2026-04-13 Root-Level Hybrid Interface Sweep Workflow

### Project Goal
- Add a new project-root workflow folder that sweeps `PROTEIN_ENV_INTERFACE_SCALE` over the real hybrid `1RKL` production path.
- Reuse `example/16.MARTINI/run_sim_1rkl.sh` as the single simulation entrypoint.
- Submit the sweep to Slurm using the same array/collector pattern already proven in `bilayer-lateral-diffusion/`.

### Architecture & Key Decisions
- Create a dedicated root-level folder `hybrid-interface-sweep/`.
- Keep the sweep workflow thin:
  - generate a manifest of `(interface_scale, replicate)` tasks,
  - run one full `run_sim_1rkl.sh` instance per task in its own task-local `RUN_DIR`,
  - collect only execution/status metadata rather than adding a new scientific analyzer in v1.
- Do not duplicate hybrid simulation logic in Python; invoke `example/16.MARTINI/run_sim_1rkl.sh` directly with environment overrides.
- Follow the `bilayer-lateral-diffusion` Slurm submission shape:
  - `workflow.py`
  - `run_local.sh`
  - `submit_remote_round.sh`
  - generated array and collector `.sbatch` scripts
- Preserve reproducibility by capturing a whitelist of relevant `run_sim_1rkl.sh` environment overrides into the manifest at `init-run` time.

### Execution Phases
- [x] Phase F: Add root and local tracker entries for the new hybrid sweep workflow.
- [x] Phase G: Create `hybrid-interface-sweep/workflow.py` with manifest, task runner, assembly, and Slurm submission.
- [x] Phase H: Add `run_local.sh`, `submit_remote_round.sh`, `README.md`, and local tracker files.
- [x] Phase I: Run static checks plus a reduced smoke run and no-submit Slurm staging.
- [x] Phase J: Record review notes and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the new workflow must reuse the real hybrid shell entrypoint rather than silently diverging from `example/16.MARTINI/run_sim_1rkl.sh`.

### Review
- Added a new root-level workflow folder:
  - `hybrid-interface-sweep/workflow.py`
  - `hybrid-interface-sweep/run_local.sh`
  - `hybrid-interface-sweep/submit_remote_round.sh`
  - `hybrid-interface-sweep/README.md`
  - local `plan.md`, `progress.md`, and `findings.md`
- The new workflow:
  - generates a manifest of `(interface_scale, replicate)` tasks,
  - runs the real `example/16.MARTINI/run_sim_1rkl.sh` entrypoint per task,
  - assigns each task its own task-local `RUN_DIR`,
  - captures a whitelist of relevant hybrid env overrides into the manifest,
  - stages Slurm array and collector scripts using the same pattern as `bilayer-lateral-diffusion`.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - no-submit Slurm staging under `/tmp/hybrid_interface_sweep_smoke`
  - reduced local smoke run with:
    - `interface_scale = 0.85`
    - `replicates = 1`
    - `MIN_60_MAX_ITER=1`
    - `MIN_61_MAX_ITER=1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 1`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - smoke result verification:
    - task result JSON recorded `success = true`
    - task-local `stage_7.0.up` exists at `/tmp/hybrid_interface_sweep_smoke/tasks/scale0p85_r01/run/checkpoints/1rkl.stage_7.0.up`
    - assembled summary recorded `1` successful task and `1` completed scale.

## 2026-04-13 Hybrid Sweep Post-Run Analysis

### Project Goal
- Add a post-run analysis path for `hybrid-interface-sweep/` that measures hybrid diffusion/fluidity signals from completed `stage_7.0.up` files.
- Make that analysis runnable locally and on Slurm.
- Produce downloadable per-task and assembled analysis artifacts for later offline review.

### Architecture & Key Decisions
- Reuse the sweep manifest/task tree and analyze existing `stage_7.0.up` outputs post hoc rather than mixing analysis into the simulation task runner.
- Measure two signals per completed sweep task:
  - protein lateral COM MSD relative to bilayer COM, using protein `CA` carrier atoms from `hybrid_bb_map`,
  - lipid `PO4` lateral MSD as a bilayer fluidity guardrail.
- Use the same basic MSD treatment as the bilayer workflow:
  - XY unwrapping,
  - bilayer COM drift removal,
  - burn-in removal,
  - linear `MSD = 4Dt + b` fit over a fixed lag window.
- Add the analysis as workflow subcommands plus a thin Slurm wrapper under `hybrid-interface-sweep/`, following the same array/collector pattern as the existing run workflow.

### Execution Phases
- [x] Phase K: Add root and local tracker entries for the analysis extension.
- [x] Phase L: Extend `hybrid-interface-sweep/workflow.py` with manifest discovery, per-task analysis, assembly, and Slurm submission for analysis.
- [x] Phase M: Add analysis wrapper/docs in `hybrid-interface-sweep/`.
- [x] Phase N: Run static checks plus reduced analysis smoke verification.
- [x] Phase O: Record review notes and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the analysis must select protein carriers and lipid molecules from the actual stage-file schema used by `run_sim_1rkl.sh`, not from assumptions copied from other workflows.

### Review
- Extended `hybrid-interface-sweep/workflow.py` with a full post-run analysis path:
  - `init-analysis`
  - `run-analysis-local`
  - `run-analysis-task`
  - `assemble-analysis`
  - `submit-analysis-slurm`
- Added analysis wrappers and docs:
  - `hybrid-interface-sweep/run_analysis_local.sh`
  - `hybrid-interface-sweep/submit_analysis.sh`
  - updated `hybrid-interface-sweep/README.md`
- Analysis outputs now measure:
  - protein lateral COM diffusion relative to the bilayer COM,
  - lipid `PO4` lateral diffusion as a bilayer guardrail,
  - per-task MSD and linear-diffusion fit metadata in `analysis/results/tasks/*.json`,
  - assembled CSV and summary JSON outputs under `analysis/assembled/`.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced real sweep run in `/tmp/hybrid_interface_sweep_analysis_smoke` with `PROD_70_NSTEPS=40`
  - local analysis run on that smoke output
  - no-submit analysis Slurm staging under `/tmp/hybrid_interface_sweep_analysis_smoke/analysis/slurm`
  - generated artifacts verified:
    - `analysis/analysis_manifest.json`
    - `analysis/results/tasks/scale0p85_r01.json`
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/summary.json`

## 2026-04-29 C++ Core Cleanup For 1RKL Default Path

### Project Goal
- Refactor and clean `src/main.cpp`, `src/box.h`, `src/box.cpp`, and `src/martini.cpp` so only logic required by the default `example/16.MARTINI/run_sim_1rkl.sh` pipeline remains.

### Architecture & Key Decisions
- Scope is strict default-path execution only; optional workflow/debug branches are removed.
- Keep SC-env and BB-env interaction potentials active; do not disable or zero these force paths.
- Keep barostat/Ewald and hybrid runtime features only where needed by the active 1RKL pipeline.
- Use `/Users/yinhan/Documents/upside2-md-master/src/main.cpp` as reference for removing temporary CLI/debug deviations in `src/main.cpp`.

### Execution Phases
- [x] Phase 1: Clean `src/main.cpp` against master and keep only required workflow deltas.
- [x] Phase 2: Remove dead/debug scaffolding from `src/box.h` and `src/box.cpp`.
- [x] Phase 3: Remove dead wrappers/placeholders/debug branches from `src/martini.cpp` while preserving active physics.
- [x] Phase 4: Build and run reduced workflow verification for default 1RKL path.
- [x] Phase 5: Document results in `progress.md` and review section.

### Known Errors / Blockers
- None at start.

### Review
- Implemented cleanup in the four requested files:
  - `src/main.cpp`: removed split-potential/Rg logging scaffolding, removed temporary CLI branches (`nvtc`, `max-force`, `martini-hold-backbone`), kept only default-path flags required by `run_sim_1rkl.sh` (`duration-steps`, minimization controls, momentum restart controls).
  - `src/box.h` and `src/box.cpp`: removed barostat debug toggles/prints and equilibrium-freeze scaffolding; preserved active NPT scaling and Ewald reciprocal computation paths.
  - `src/martini.cpp`: removed dead placeholder and compatibility scaffolding (`martini_masses::martini_integration_cycle`, `ConjugateGradientMinimizer` node, `minimize_structure_with_regular_potential` stub), removed SC-energy-dump plumbing and nonprotein hard-sphere experimental branch, and removed debug/status prints.
- Physics integrity preserved:
  - SC-env and BB-env force paths remain active.
  - `martini_potential` and `martini_sc_table_1body` both remain present in generated production stage files.
- Verification:
  - `cmake --build obj --target upside` succeeded.
  - Reduced default-workflow run succeeded end-to-end:
    - `RUN_DIR=/tmp/cpp_cleanup_smoke` with minimized stage step counts.
    - Stage outputs through `1rkl.stage_7.0.up` and VTF extraction completed.
  - Stage-7 production file check confirmed:
    - `/input/potential/martini_potential` exists;
    - `/input/potential/martini_sc_table_1body` exists.

## 2026-04-29 Python Workflow Cleanup (Baseline + Explicit 7.x Continuation)

### Project Goal
- Refactor and clean `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, and `py/martini_extract_vtf.py` to remove dead/legacy/debug-hoarded logic.
- Keep only functionality required by `example/16.MARTINI/run_sim_1rkl.sh` fresh execution and explicit production continuation (`--continue-stage-70-from` -> `7.1`, `7.2`, ...).
- Preserve SC-env and BB-env interface interactions fully active.

### Architecture & Key Decisions
- Keep `run-hybrid-workflow` as the only intended high-level entrypoint in `martini_prepare_system.py`.
- Remove legacy command surfaces and wrapper flows not needed by the active pipeline.
- Keep continuation path only for explicit source file input; remove auto-discovery branches.
- Flatten defensive/legacy branching in VTF extraction to a direct path used by workflow outputs.

### Execution Phases
- [x] Phase 1: Prune `martini_prepare_system.py` CLI/command surfaces and continuation discovery branches.
- [x] Phase 2: Remove dead standalone/legacy conversion layers from `martini_prepare_system_lib.py` while preserving required functions.
- [x] Phase 3: Simplify `martini_extract_vtf.py` to workflow-required VTF behavior.
- [x] Phase 4: Run compile/syntax checks and validate retained continuation + physics-critical attrs remain intact.
- [x] Phase 5: Record review results and residual risks.

### Known Errors / Blockers
- None at task start.

### Review
- Implemented:
  - `py/martini_prepare_system.py` now supports only the `run-hybrid-workflow` command surface; legacy subcommands and fallback execution paths were removed.
  - `run-hybrid-workflow` parser now matches flags passed by `run_sim_1rkl.sh`; internal physics and stage defaults moved to explicit in-code constants.
  - Continuation now requires explicit `--continue-stage-70-from`; auto-discovery (`previous-run`/`auto-continue`) code paths were removed and now fail fast if used.
  - `py/martini_prepare_system_lib.py` removed dead standalone prep CLI/`main` orchestration and unused wrapper functions (`create_production_input`, `main_always_fixed`, martinize wrapper path).
  - `py/martini_extract_vtf.py` was simplified to VTF-only extraction and direct trajectory-group handling; defensive NaN frame substitution and unused PDB/output-group branches were removed.
- Physics integrity preserved:
  - production hybrid-control writes still enforce active SC-env configuration (`sc_env_*` attrs plus positive interface scale) and keep `production_nonprotein_hard_sphere=0`;
  - production stage injection still requires and validates active `martini_sc_table_1body` and `martini_potential`.
- Verification:
  - `.venv/bin/python -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py`
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --help`
  - `.venv/bin/python py/martini_extract_vtf.py --help`
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --previous-run-dir foo` (verified explicit-source continuation enforcement with clear error).

## 2026-04-29 Hybrid MARTINI Architecture Simplification (AA->Backbone-Only, No CG Protein)

### Project Goal
- Refactor the hybrid MARTINI workflow to remove the AA->CG martinize pipeline entirely and run hybrid staging from direct AA input while keeping only protein backbone carriers (`N/CA/C/O`) in runtime stage files.
- Enforce protein rigid-body behavior (free translation/rotation, fixed internal geometry) during stage 6.x and release that constraint for stage 7.x production.
- Keep SC-env and BB-env interface interactions active in both stages.

### Architecture & Key Decisions
- Delete `py/martini_martinize.py` and remove all martinize execution/config from workflow scripts.
- Runtime protein representation in both 6.x and 7.x is backbone-only (`N/CA/C/O`), with no CG protein particles.
- BB-env is evaluated from per-residue virtual BB COM built from `N/CA/C/O` carriers and force-distributed back to those carriers.
- Stage 6.x hybrid activation must be on; stage 7.x remains hybrid-active with rigid constraint released.
- Use `pydssp` for residue secondary-structure assignment and inline the MARTINI backbone type mapping logic needed for BB typing.

### Execution Phases
- [ ] Phase 1: Remove martinize surfaces and add DSSP dependency.
- [ ] Phase 2: Refactor Python prep to AA-backbone-only runtime packing/mapping and stage artifact generation.
- [ ] Phase 3: Implement C++ virtual-BB interface force path and rigid-body stage-6 constraint path.
- [ ] Phase 4: Clean dead code paths and remove deleted script references.
- [ ] Phase 5: Run build/syntax/workflow verification and document review/results.

### Known Errors / Blockers
- None at start.

## 2026-04-29 AA-Direct Hybrid MARTINI Simplification

### Project Goal
- Remove the AA->CG martinize conversion path from the 1RKL hybrid workflow.
- Use AA backbone atoms directly for stage 6.x/7.x hybrid interaction with dry-MARTINI bilayer.
- Keep stage 6.x protein rigid as a free rigid body (translation/rotation allowed), then release that constraint in stage 7.x.

### Architecture & Key Decisions
- Delete `py/martini_martinize.py`; keep only DSSP-driven AA-backbone -> dry-MARTINI BB type mapping inside `py/martini_prepare_system_lib.py`.
- Make `run_sim_1rkl.sh` call `run-hybrid-workflow` with AA input only (`--protein-aa-pdb`) and no martinize/ITP path.
- Keep hybrid stage control as `activation_stage=minimization` with `preprod_protein_mode=rigid_body` for stage 6.x and switch to `activation_stage=production` at stage 7.x.
- Implement runtime rigid-body group projection in `src/martini.cpp` through `martini_fix_rigid::set_dynamic_rigid_groups(...)` instead of absolute position pinning.

### Execution Phases
- [x] Phase 1: Remove martinize CLI and file dependencies from workflow scripts.
- [x] Phase 2: Keep only AA backbone BB mapping path in prep code and remove dead CG mapper helpers.
- [x] Phase 3: Ensure stage 6.x rigid-body hold and stage 7.x release in C++ hybrid transitions.
- [x] Phase 4: Validate syntax/build and check for stale references.

### Known Errors / Blockers
- End-to-end smoke run is blocked in current local `.venv` by missing optional dependency `tables` (`ModuleNotFoundError`).

### Review
- Deleted `py/martini_martinize.py`.
- Refactored prep/workflow to AA-only input and backbone-only mapping.
- C++ rigid-body dynamic group path compiles and is wired into stage transition control.
- Static verification passed (`py_compile`, `cmake --build obj --target upside`).

## 2026-04-30 Rigid Invariance and Handoff Geometry Integrity

### Project Goal
- Enforce mathematically rigid protein backbone geometry during minimization and stage 6.x preproduction.
- Guarantee stage 6.6 -> 7.0 handoff copies AA backbone coordinates exactly (no per-residue coordinate rewrites).

### Architecture & Key Decisions
- Keep stage-6 rigid behavior as true rigid-body projection in C++ minimization path, including trial line-search states.
- Remove coordinate refresh paths that can overwrite AA backbone carriers during handoff/runtime.
- Keep BB-env and SC-env interface forces active; only remove coordinate mutation paths.

### Execution Phases
- [x] Phase 1: Patch minimization to rigid-project trial/accepted/final states.
- [x] Phase 2: Remove Python handoff hybrid-carrier refresh from active path.
- [x] Phase 3: Remove runtime hybrid coordinate overwrite call before MARTINI pair force evaluation.
- [x] Phase 4: Build and run geometry invariance checks.

### Known Errors / Blockers
- None in this patch round.

### Review
- `src/martini.cpp`: minimization now rigid-projects before first compute, each line-search trial, and final state flush.
- `py/martini_prepare_system.py` / `py/martini_prepare_system_lib.py`: removed `UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS` path and `refresh_hybrid_reference_carriers`.
- `src/martini.cpp`: removed hybrid pre-force coordinate refresh call, eliminating runtime coordinate overwrites.
- Verification:
  - `cmake --build obj --target upside` passed.
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py` passed.
  - Stage 6.0 minimization invariance check: pair-distance drift mean `1.2e-07 A`, max `1.9e-06 A`.
  - Stage 6.6 -> 7.0 handoff check: protein atom coordinate diff mean/max `0.0 A`.

## 2026-04-30 Stage-7 Helix Start Forensics

### Project Goal
- Identify where AA-backbone helix geometry is lost before stage `7.0` and align workflow defaults with the last known-good no-preprod-minimization behavior.

### Architecture & Key Decisions
- Validate geometry from checkpoint files directly (not only VMD rendering).
- Keep stage `6.0/6.1` minimization opt-in rather than default, matching the user’s last-working expectation.

### Execution Phases
- [x] Phase 1: Compare AA backbone vs stage `6.0/6.6/7.0` checkpoint geometry.
- [x] Phase 2: Validate stage-7 VTF first-frame helix assignment independently.
- [x] Phase 3: Disable stage `6.0/6.1` minimization by default and keep explicit opt-in.
- [x] Phase 4: Re-run reduced fresh workflow and verify stage-7 start helix.

### Known Errors / Blockers
- None in this patch round.

### Review
- Forensic result on current artifacts: stage `7.0` start backbone is still helix-consistent with AA reference; no structural corruption observed in checkpoint geometry.
- Changed defaults so stage `6.0/6.1` minimization is skipped unless user sets positive iteration counts.
- Fresh reduced run with new defaults confirmed stage-7 input and VTF first frame both retain helix signature (`23` helical residues by DSSP on 31-residue backbone).
- Added strict stage-7-start vs PDB comparison on a fresh longer run:
  - `/tmp/helix_longer_check/checkpoints/1rkl.stage_7.0.up` matches PDB backbone geometry up to rigid transform within micro-Angstrom numerical noise.
  - Existing older run artifacts under `example/16.MARTINI/outputs/martini_test_1rkl_hybrid` still show millangstrom-level drift and should not be used as baseline for patched behavior.

## 2026-04-30 BB Type + Interface Logic Parity (AA Runtime)

### Project Goal
- Keep AA-backbone runtime mode while restoring last-working BB typing policy and BB-env/SC-env interaction semantics.

### Architecture & Key Decisions
- Remove DSSP dependency from BB type assignment; use martinize fallback with fixed coil mapping plus termini/chain-break `Qd/Qa` overrides.
- Preserve AA virtual-BB COM position evaluation from `hybrid_bb_map` carriers (`N/CA/C/O`) and keep gradient projection back to carriers.
- Align interaction control flow with last-working behavior by removing the newly-added cross-interface role filter in `martini_potential`.

### Execution Phases
- [x] Phase 1: Replace DSSP-based BB typing with martinize fallback typing.
- [x] Phase 2: Wire fallback BB typing into both hybrid mapping export and stage conversion.
- [x] Phase 3: Restore BB-env/SC-env interaction control flow parity in C++ while keeping AA COM projection.
- [x] Phase 4: Verify Python syntax and C++ build.

### Known Errors / Blockers
- None in this patch round.

### Review
- `py/martini_prepare_system_lib.py` now computes BB types with martinize fallback (`ss="C"` table behavior) and applies termini/chain-break charged BB rules; DSSP usage removed from this path.
- `py/martini_prepare_system.py` now consumes the fallback mapper and reports `bb_fallback_typed_residues` in prep summaries.
- `src/martini.cpp` now uses last-working BB-env/SC-env pair-routing flow (removed added cross-interface role gate) while retaining AA virtual-BB COM site projection for carrier atoms.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed.
  - `cmake --build obj --target upside` passed.

## 2026-04-30 BB Proxy COM + Force Distributor Restoration

### Project Goal
- Keep current AA-runtime framework while restoring last-working BB semantics:
  - BB is a dedicated virtual proxy (not CA);
  - BB position is COM of N/CA/C/O;
  - BB-env force is distributed back to N/CA/C/O carriers.

### Architecture & Key Decisions
- Preserve explicit virtual `BB` atoms in runtime protein coordinates and keep `ROLE_BB` assignment only for `BB` atom names.
- Keep mapped carrier atoms (`N/CA/C/O`) excluded from direct protein-env MARTINI pair coupling to avoid double-counting alongside BB-env projection.
- Ensure BB gradient projection is not startup-only; projection must stay active through production stage interaction flow whenever a BB map row exists.

### Execution Phases
- [x] Phase 1: Confirm last-committed BB semantics in `src/martini.cpp` and runtime mapping schema.
- [x] Phase 2: Restore virtual-BB export path in Python prep and allow BB particle types in stage conversion.
- [x] Phase 3: Patch C++ pair loop so BB gradients are always projected to carriers in active hybrid stages.
- [x] Phase 4: Rebuild and run syntax checks.

### Known Errors / Blockers
- No long-horizon stage-7 replay run executed in this patch round yet.

### Review
- `src/martini.cpp` now routes mapped BB gradients to carrier atoms (`N/CA/C/O`) independent of startup ramp state.
- Runtime keeps BB role separate from CA role and refreshes BB coordinates from mapped COM in active hybrid.
- Build and syntax checks passed (`cmake --build obj --target upside`, `python -m py_compile`).
