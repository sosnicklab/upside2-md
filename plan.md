# Project Goal

Build and maintain a physically defensible single-vector DOPC coarse-grained lipid (CGL) force field for the UPSIDE/dry-MARTINI hybrid workflow. The active cleanup task is to keep the hybrid C++/Python implementation aligned with the current tempered-PMF/direct-geometry schema by removing unused compatibility paths, stale versioned metadata, scaffold code, and dead controls without changing the physical interaction model.

# Architecture & Key Decisions

- Runtime CGL orientation must come from the vector-particle spline interactions, especially CGL-CGL and SC-CGL. Do not add a standalone CGL orientation potential.
- Table generation must evaluate direct dry-MARTINI nonbonded energies from rotated full-resolution bead geometries over the sampled runtime direction-vector coordinates.
- CGL-CGL tables must retain the resolved dry-MARTINI lipid-lipid attraction in the spline. Do not subtract it into an absent radial background or replace it with a separate orientation correction.
- Revised Decision: CGL-CGL unresolved azimuthal/bead-frame samples should use a tempered direct PMF as the next two-body model. After fixing runtime log-control injection, a 50k bilayer-only validation with the production-temperature Boltzmann pair PMF still collapsed same-leaflet packing (`p05 ~3.86 A`) and produced a flip, so the pair PMF remains too favorable for independently chosen hidden axial rotations in a dense single-vector bilayer. Direct energy expectation was too stiff and launched the bilayer apart. A tempered PMF keeps direct dry-MARTINI energies and the three runtime coordinates while avoiding capping, normalization, hidden relaxation, or an added orientation potential.
- CGL axial bead-frame rotation is an unresolved coordinate in the single-vector lipid. CGL-CGL and CGL-SC table generation should sample multiple rotations around the CGL axis; CGL-CGL uses a tempered two-body PMF for dense-bilayer transferability.
- Revised Decision: SC-CGL now uses the extended-support full tensor and invertible transformed control, analogous to the successful CGL-CGL/CGL-target representation. This is still a two-body direct-geometry spline table and does not introduce normalization, capping, hidden relaxation, interaction scaling, or a standalone CGL orientation potential.
- Phase 2C implementation result: the accepted hybrid surface uses direct rotated dry-MARTINI geometry for CGL-CGL, SC-CGL, CGL-particle, and SC-particle. Active table builds use only a near-zero numerical singularity guard (`1e-6 nm`) rather than physical-distance capping.
- Revised Decision: SC-particle must also move away from the separable radial plus angular-profile factorization when validating full-resolution lipids. Fresh full 1RKL exposed a factorization artifact in `rotamer_angular_energy_kj_mol` down to `-149236.78 kJ/mol`, which is not a physical dry-MARTINI well. The next implementation target is a direct full radial-by-angular SC-particle table evaluated by the runtime node, with no force/energy caps, interaction scaling, or hidden relaxation.
- Revised Decision: Active dry-MARTINI table builds must not floor pair distances at `0.10 nm` or similar values. Use only a near-zero numerical singularity guard (`1e-6 nm`) so the generated spline represents the physical dry-MARTINI core instead of a capped core.
- CGL-CGL needs denser unresolved sampling than the initial `2^2` azimuthal / fixed axial table. The `4`-frame axial PMF produced the large energy reduction, and `4^2` azimuthal sampling improved orientation metrics but did not stabilize the bilayer; the next transferable two-body step is denser axial bead-frame PMF sampling.
- CGL-CGL unresolved short-range rows should preserve angular dependence from the first sampled radial shell instead of copying one worst-case angular clash into every orientation.
- CGL-CGL radial support must cover the extended resolved lipid bead geometry, not only the CGL COM contact distance. A COM cutoff of `16.8 A` truncates physically relevant tail-tail contacts for head-to-tail CGL vectors and leaves the bilayer smoke dominated by cross-leaflet repulsion.
- CGL-particle radial support must also cover the extended resolved lipid bead geometry, not only the CGL COM contact distance. A target can interact with any resolved DOPC bead represented by the CGL vector.
- Coarse initial CGL bilayers must not inherit wrapped/full-lipid fixture artifacts. Preparation conditions CGL/CGLD leaflet membership from the CGL direction sign when available and sets the default CGL leaflet COM separation to the DOPC-derived display-tail zero-gap scale while moving orientation sites with their parent CGL particles.
- Coarse initial same-leaflet CGL spacing should use the box area per leaflet lipid when available, not only the bead contact distance. The contact distance is a lower bound for excluded volume, not a physically packed bilayer lattice spacing.
- Coarse initial CGL conditioning in hybrid systems must also preserve protein/ion target clearance. Lipid COM moves that improve CGL-CGL packing must not create CGL-target overlaps after hybrid packing.
- Hybrid CGL initial conditioning must clear targets against the swept resolved DOPC bead geometry, not only against the CGL COM. The CGL-particle table integrates unresolved axial bead-frame rotations, so initial placement must avoid target overlaps with that swept two-body geometry.
- Hybrid CGL initial conditioning should use the full DOPC-derived contact distance for swept-bead CGL-target clearance. The earlier `0.65 * contact` target left charged ions/proxies inside the physical direct CGL-target hard core, producing large but legitimate table energies rather than a table-generation defect.
- The priority model remains a transferable two-body CGL method, because the same direct rotated-geometry idea must extend from CGL-CGL to CGL-SC where many-neighbor normalization has no clear analogue.
- CGL-CGL many-neighbor normalization is a fallback only if the two-body table semantics are exhausted. If used, it must be explicit table metadata derived from a calibration rule, not an environment-only tuning knob.
- Spline fits must be well-conditioned. Angular B-spline control counts should track the angular data resolution instead of greatly exceeding it.
- CGL-particle target splines need enough angular resolution to localize the charged headgroup hard core near `ang=-1`. Coarse angular controls smear the physical headgroup singularity into nearby ion orientations and create unphysical `1e9 E_up` excluded-volume values at ~4 A swept-bead clearance.
- CGL-particle target splines must also use enough radial and angular fit resolution to prevent hard-core sampled rows from bleeding into non-overlap head-side geometries. A 1RKL probe found raw direct dry-MARTINI energy `-5.4 kJ/mol` at `r~19 A, ang~-0.936`, while the old fitted table evaluated the same geometry near `1.7e8 E_up`.
- Revised Decision: CGL-particle target PMFs should use an invertible `log1p((E - E_ref) / kBT)` spline control instead of raw Boltzmann weights. The float32 Boltzmann-weight representation underflows in hard-core regions and creates an artificial flat ~71 `E_up` plateau with weak restoring force. The log1p reduced-free-energy control remains the same two-body physical PMF surface, avoids the old hard-core bleed-through, and does not introduce capping, normalization, hidden relaxation, or an orientation potential.
- CGL direct-geometry spline evaluators must not flatten the short-range radial coordinate. The existing clamped radial B-spline returns zero derivative below one knot spacing, which acts as a hard-core force cap even when the table contains a physical repulsive core. CGL radial evaluation should linearly extrapolate the first in-domain spline slope below the first knot while still tapering normally at the cutoff.
- CGL-CGL short-range validation must check the fitted/runtime PMF, not only raw table metadata. A physical direct PMF can still fail if the spline representation smears or flattens the hard core enough for minimization to find collapsed same-leaflet states. The next acceptable fixes are better two-body table representation choices, such as denser radial support or an invertible transformed control, not capping, normalization, hidden relaxation, or an added orientation potential.
- Revised Decision: CGL-CGL direct sampling must include the physical overlap core. The previous table started at `5 A`, where some unresolved axial/azimuthal PMF sectors are already attractive, then extrapolated everything below that radius. The CGL-CGL table should start at the first radial knot inside the core and use a near-zero numerical bead-distance guard only to avoid division by zero in exact overlaps.
- Revised Decision: CGL-CGL should fit an invertible `log1p((E - E_ref) / kBT)` reduced-PMF control, matching the CGL-target representation. Raw-energy fitting of the short-core table smeared true extreme-angle hard cores into ordinary long-range angular sectors; the transform keeps the same physical two-body PMF while making the spline numerically representable.
- Revised Decision: Runtime CGL-CGL injection must copy the pair log-control metadata and `reference_energy_eup` into the `.up` node. The generated `dopc.h5` may contain the physical inverse-transform table, but if `inject_cg_lipid_nodes()` omits `log1p_reduced_transform`, `boltzmann_temperature_upside`, or `reference_energy_eup`, Upside interprets spline controls as energies and silently weakens the hard core by orders of magnitude.
- CGL-target charged-BB and ion interactions must use the same direct physical CGL-target table as other targets. A runtime nonnegative excluded-volume projection is not a physical potential; it removes attractive/canceling parts of the resolved DOPC electrostatics and leaves one-sided headgroup repulsion.
- Extended CGL-CGL radial tensors need explicit smoothing verification because the larger COM range can amplify gradients. The `fit_smooth=0.5` trial worsened hidden orientation-site stretching, so the current best-known two-body setting is restored to `fit_smooth=0.1` while the remaining CGL-CGLD length issue is diagnosed separately.
- CGLD is a numerical orientation-coordinate carrier for a unit vector, not a physical dry-MARTINI bead. Its harmonic CGL-CGLD bond must be stiff enough to keep the hidden carrier near the derived vector length under two-body table torques. Keep the dry-MARTINI projected bonded stiffness as metadata, but use an explicit `2.0x` carrier-constraint stiffness in runtime preparation unless a true holonomic constraint is added.
- Dry-MARTINI training artifacts remain in native units. Runtime conversion into Upside units happens in the table-building/simulation path.
- All simulation interactions computed by Upside must be spline-table based, except existing Ewald behavior.
- Full workflow verification must use regenerated handoff artifacts, not stale output directories or pre-minimization prepared files.
- Hybrid protein-lipid packing clearance defaults should be derived from the active dry-MARTINI DOPC nonbonded contact distance, not from a stale hardcoded Angstrom value. This is a physical dry-MARTINI exclusion distance used during preparation, not an interaction scaling or cap.
- Full-resolution hybrid protein BB proxy types must come from structure-derived secondary structure when no protein ITP is supplied. The previous coil-only martinize fallback typed all 1AFO helical BB proxies mostly as `P5`, making full-lipid BB-env interactions physically wrong for a transmembrane helix while keeping the interactions enabled.
- Stage-6 protein position restraints must select protein atoms (`protein_membership >= 0`), not environment atoms.
- Stage-7 production-Hamiltonian burn-in may use protein position restraints as an equilibration protocol, provided SC-env and BB-env interactions remain present at full scale and the restraints are removed before the recorded production segment. This is not an interaction scale, cap, or exclusion; it prevents the equilibration handoff from using a distorted protein endpoint as the production starting structure.
- Experimental Decision: Phase 10 deliberately applies production-temperature PMF hidden-state averaging, without tempering, to CGL-CGL, SC-CGL, CGL-particle, and SC-particle. This is the user's preferred one-method model if it validates. Runtime interactions remain full-strength and spline/table based; no SC-env/BB-env interaction may be disabled, capped, or scaled down to force a pass.
- Revised Experimental Decision: Phase 11 uses the same `tau=10.0` tempered PMF hidden-state reduction for CGL-CGL, SC-CGL, CGL-particle, and SC-particle. This is a uniform configurational PMF over unresolved rigid orientations, not a guarantee of correct real-time dynamics. Validation must therefore check structure and dynamics-adjacent observables rather than assuming the PMF preserves kinetics.

# Execution Phases

- [x] Phase 16: Add exact pairlist acceleration for CGL-specific runtime nodes.
  - [x] Implement cached pairlists for CGL-CGL, SC-CGL, CGL-rotamer-SC, and CGL-target loops using each table's existing cutoff plus a rebuild buffer.
  - [x] Keep force-field semantics unchanged: no interaction removal inside the table cutoff, no shorter physical cutoff, no cap/scale/orientation potential.
  - [x] Build C++ and run focused smoke/short stability tests comparing runtime behavior.
  - [x] Benchmark coarse 1RKL against the current full-lipid timing enough to confirm the bottleneck is reduced.

- [x] Phase 15: Clean up `cg_lipid_potentials.tex` to match the cleaned hybrid implementation.
  - [x] Inspect the current TeX for patch-history language, stale mixed-method descriptions, removed attrs/files, and contradictions with regenerated H5 metadata.
  - [x] Rewrite the method flow as one coherent derivation for the current direct-geometry tempered-PMF tables.
  - [x] Verify the TeX compiles and document any remaining caveats.

- [x] Phase 14: Clean up hybrid-interface C++ and Python code.
  - [x] Inventory active hybrid-interface source files and generated-schema entry points.
  - [x] Identify unused scaffolding, obsolete compatibility branches, stale experiment/version metadata, and dead helper code in `src/martini*`, `src/box.*`, and `py/martini_*.py`.
  - [x] Remove only code proven unused by call-site search or superseded active schema; keep runtime physics and spline-table semantics unchanged.
  - [x] Run Python compile checks, C++ build checks, regenerated-H5 metadata audits, and at least a lightweight workflow preparation/simulation smoke test.
  - [x] Document remaining cleanup risks or intentionally retained legacy interfaces.

- [x] Phase 13: Diagnose reported CGL bilayer visual gap and messy orientation.
  - [x] Identify the relevant current CGL output directories under `example/16.MARTINI/outputs`.
  - [x] Quantify physical CGL leaflet geometry: COM separation, centerline tail gap, resolved bead-envelope gap, orientation distributions, flips, and crossings.
  - [x] Compare physical HDF5 geometry against VTF display geometry to determine whether the gap is real or visualization-driven.
  - [x] Decide whether a physical preparation/table fix or visualization/export adjustment is warranted.

- [x] Phase 12: Physical-integrity audit of the uniform tempered-PMF implementation.
  - [x] Audit source code for twist/torsion coordinates, force caps, arbitrary interaction scaling, and standalone CGL orientation potentials.
  - [x] Audit installed H5 table metadata and datasets for cap/scaling/orientation-potential markers.
  - [x] Audit generated runtime `.up` files from the focused Phase 11 run for the same physical-integrity markers.
  - [x] Report whether the active implementation satisfies the physical-design constraints.

- [x] Phase 11: Universal tempered-PMF force-field experiment with dynamics caveat.
  - [x] Revise generator defaults so CGL-CGL, SC-CGL, CGL-particle, and SC-particle all use `tau=10.0` tempered PMF.
  - [x] Rebuild `parameters/dryMARTINI/dopc.h5` and `parameters/dryMARTINI/sidechain.h5`.
  - [x] Audit rebuilt H5 metadata for all four interaction classes.
  - [x] Run focused coarse validation first, then decide whether to run the remaining workflows based on bilayer/protein stability.
  - [x] Quantify dynamics-sensitive observables available from the generated trajectories, including CGL lateral MSD and orientation autocorrelation when the output cadence is sufficient.
  - [x] Decide whether the dynamics mismatch can be handled by reporting/calibrating an effective time scale, or whether the force-field reduction itself must change.
  - [x] Report structural stability and the limitations of interpreting dynamics from tempered-PMF trajectories.

- [x] Phase 10: Universal production-PMF force-field experiment.
  - [x] Route CGL-CGL, SC-CGL, CGL-particle, and SC-particle table generation through production-temperature PMF hidden-state averaging.
  - [x] Rebuild `parameters/dryMARTINI/dopc.h5` and `parameters/dryMARTINI/sidechain.h5`.
  - [x] Audit rebuilt H5 metadata for all four interaction classes.
  - [ ] Run `example/16.MARTINI/run_sim_1afo.sh`.
  - [x] Run `example/16.MARTINI/run_sim_1rkl.sh`.
  - [ ] Run `example/16.MARTINI/run_sim_1afo_full.sh`.
  - [ ] Run `example/16.MARTINI/run_sim_1rkl_full.sh`.
  - [x] Quantify bilayer stability, CGL orientation, and protein Rg/hbond for the first completed coarse workflow.
  - [x] Report whether the one-method production-PMF model works across all pairs.

- [x] Phase 9: Universal tempered-PMF force-field experiment.
  - [x] Route CGL-CGL, SC-CGL, CGL-particle, and SC-particle table generation through `tau=10.0` hidden-state PMF averaging.
  - [x] Rebuild `parameters/dryMARTINI/dopc.h5` and `parameters/dryMARTINI/sidechain.h5`.
  - [x] Audit rebuilt H5 metadata for all four interaction classes.
  - [x] Run `example/16.MARTINI/run_sim_1afo.sh`.
  - [x] Run `example/16.MARTINI/run_sim_1rkl.sh`.
  - [x] Run `example/16.MARTINI/run_sim_1afo_full.sh`.
  - [x] Run `example/16.MARTINI/run_sim_1rkl_full.sh`.
  - [x] Quantify bilayer stability, CGL orientation, protein Rg/hbond, and secondary-structure stability.
  - [x] Iterate only with physically justified table/representation fixes if validation fails.
  - [x] Report whether one tempered-PMF method works across all pairs.

- [x] Phase 8: Clean up `cg_lipid_potentials.tex` after the CGL-CGL/SC-CGL averaging discussion.
  - [x] Rewrite unresolved-coordinate equations as one coherent method.
  - [x] Make CGL-CGL tempered PMF and SC-CGL weighted energy expectation explicit.
  - [x] Verify the TeX compiles.

- [x] Phase 7: Fix 1AFO full-lipid stage-7 secondary-structure distortion.
  - [x] Reproduce the reported VTF distortion against HDF5 coordinates and stage logs.
  - [x] Trace onset to stage-7 production burn-in, not VTF export or stage-6 preparation.
  - [x] Implement structure-derived MARTINI BB proxy typing and mapping provenance.
  - [x] Fix the protein position-restraint selector bug.
  - [x] Run compile/syntax checks.
  - [x] Diagnose remaining stage-7 burn-in force source after structure-derived BB typing.
  - [x] Implement the next physically justified fix.
  - [x] Regenerate and validate a fresh 1AFO full-lipid workflow.
  - [x] Replace the active output only after validation passes.

# Known Errors / Blockers

- None for Phase 14. The initial rebuild check exposed two stale cleanup surfaces and both were fixed: `example/16.MARTINI/build_martini_h5_m1.sh` still passed removed fit-relax controls, and `py/martini_prepare_system.py` still wrote scalar protein position restraints while the cleaned runtime expects `spring_const_xyz`.
- Phase 11 dynamics caveat: the focused 1RKL coarse structural diagnostic passes after a corrected production restart, but CGL lateral COM dynamics are faster than the active full-resolution DOPC COM reference at longer short-time lags. At lag `15.0` time units, CGL MSD is `1.400 A^2` and full DOPC COM MSD is `0.470 A^2`, implying an approximate effective-time scale factor of `0.34` for that observable on this short trajectory. This is a calibration issue, not something the tempered PMF itself can guarantee.
- Phase 13 result: the current active coarse VTF gap is not a physical bilayer void. Direction-sign leaflet assignment gives no CGL flips or leaflet crossings in active `martini_1afo_hybrid` or `martini_1rkl_hybrid`; final tail-centerline gaps are negative (`-4.897 A` for 1AFO, `-3.807 A` for 1RKL), so the CGL tail centerlines overlap across the midplane. The visible gap comes from sparse vector rendering and viewer treatment of the synthetic CGL display radius, not from missing central attraction. The orientation looks somewhat noisy because CGL rods tilt thermally; final aligned-z min/p05/mean is `0.478/0.846/0.947` for 1AFO and `0.711/0.837/0.943` for 1RKL.
- Phase 12 physical-integrity audit passed. Installed `dopc.h5` and `sidechain.h5`, plus the Phase 11 prepared/production `.up` files, have no CGL twist/orientation-potential node, no active force caps, no arbitrary interaction scaling, no hidden relaxation, no excluded-area/nonnegative projection, and no active duplicate-row normalization. A stale unread `nonprotein_hs_force_cap=100` metadata default was found and cleaned to `0.0`.
- Phase 10 result: the universal production-temperature PMF experiment fails on the first coarse validation workflow, `run_sim_1rkl.sh`. The workflow completed, but the final bilayer/protein metrics are unacceptable: CGL aligned-z min/p05/mean `-0.063/0.448/0.821`, `bad_parallel=7`, `bad_flip=1`, leaflet crossings `4/4`, same-leaflet NN min/p05 `2.285/2.886 A`, protein hbond first/final/min/last20 `28.74/6.93/0.56/6.43`, and protein Rg first/final/last20 `12.56/10.40/10.37 A`. The remaining three workflows are intentionally not launched until the model decision is revised, because the one-method production-PMF model already fails the required acceptance criteria.
- Phase 9 result: the universal `tau=10.0` tempered-PMF experiment fails the coarse CGL-orientation criterion. All four workflows completed, and protein/full-lipid metrics stayed stable, but final coarse CGL orientations contained flipped rods: `1AFO` had `bad_parallel=2`, `bad_flip=2`; `1RKL` had `bad_parallel=4`, `bad_flip=4`. Same-leaflet spacing and leaflet crossing checks did not indicate bilayer collapse, so the failure is specifically orientation transferability of the one-method table scheme. No physical fix was applied because the direct conclusion is that applying the same tempered-PMF reduction to all four pair classes does not work as a production model.
- None for Phase 8. The methods TeX now distinguishes hidden rigid-body orientations from internal DOPC conformations and derives weighted expectation and PMF reductions; the Phase 11 update supersedes the earlier mixed-method text with uniform `tau=10.0` tempered PMF.
- None for Phase 7. The structure-typed but unrestrained-burn-in validation `outputs/phase7_1afo_full_secondary_types` failed and was not promoted. The accepted validation `outputs/phase7_1afo_full_burnin_restraint` uses full-strength SC-env/BB-env interactions, applies protein position restraints only during stage-7 burn-in, removes those restraints before production, and has replaced the active `outputs/martini_1afo_hybrid_full` output. The previous active output is backed up as `outputs/martini_1afo_hybrid_full.pre_phase7_burnin_restraint_backup_20260611_152618`.

# Review

- Phase 16 implementation keeps the physical table cutoffs unchanged and adds only exact neighbor-cache pruning outside each table cutoff plus a `pairlist_buffer_ang` rebuild buffer. Focused checks passed: C++ rebuild, bilayer-only smoke validation, full 1RKL workflow smoke, and a 5000-step direct 1RKL stage-7 stability run.

- Phase 15 TeX cleanup replaced the patch-history methods text with a single current-method description covering the CGL coordinate, rigid resolved dry-MARTINI energy, uniform `tau=10.0` tempered-PMF hidden-orientation reduction, log1p spline transform, CGL-CGL, SC-CGL, CGL-particle, SC-particle, runtime ownership, numerical evaluation, validation scope, and units. `pdflatex -interaction=nonstopmode -halt-on-error cg_lipid_potentials.tex` passed from the repo environment; LaTeX reported only overfull-box warnings from long HDF5 schema strings.
- Phase 14 cleanup passed: Python compile for touched Martini modules, `bash -n` for MARTINI workflow wrappers, full C++ `make -C obj -j4`, regenerated `parameters/dryMARTINI/{particle,sidechain,dopc}.h5`, stale metadata audit for cap/scale/relax/schema markers, `example/16.MARTINI/test_cg_bilayer/run_test.sh`, and a short isolated 1RKL hybrid smoke in `outputs/phase14_1rkl_smoke2` through stage 6.0, stage 7.0 burn-in, stage 7.0 production, and VTF extraction. Runtime stage files contain SC-CGL full tensor `(18, 1, 2541)` and CGL-target `(1, 38, 25038)` tables with no removed force-cap/interface-scale/fit-relax attrs.
- Phase 8 TeX cleanup passed two `pdflatex -interaction=nonstopmode -halt-on-error cg_lipid_potentials.tex` runs from the repo-root environment. LaTeX reported only small pre-existing-style overfull boxes, not errors.
- Phase 11 uniform tempered-PMF implementation passed Python compile checks, H5 metadata validation for CGL-CGL/SC-CGL/CGL-particle/SC-particle, and a `pdflatex -interaction=nonstopmode -halt-on-error cg_lipid_potentials.tex` check. Focused 1RKL coarse structural validation passed, but CGL lateral MSD remains faster than the full-resolution DOPC COM reference at longer short-time lags, so lipid dynamics need effective-time calibration or a separate dissipative/friction model.
- Phase 12 audit checks passed after setting the stale unread `nonprotein_hs_force_cap` default and Phase 11 checkpoint attributes to zero. `py_compile` passed for the touched Python files. The audit verified CGL-CGL, SC-CGL, CGL-particle, and SC-particle tables use `tau=10.0` tempered PMF with `fit_relax_steps=0`, `sample_dist_min_nm=1e-6`, no cap/scale attrs, and no twist attrs; runtime `martini_potential/force_cap=0`, `protein_env_interface_scale=1`, SC-env force caps are zero, and no CGL orientation/twist potential nodes are present.
- Phase 13 diagnostic used stored `compose_vector6d/direction` signs for leaflet assignment. Median-z leaflet splitting can falsely report flips in wrapped/tiled CGL outputs; with the stored signs the active 1AFO and 1RKL coarse outputs have zero flips, zero crossings, overlapping tail centerlines, and acceptable same-leaflet spacing. No force-field or preparation change was made.
- Active `martini_1afo_hybrid_full` stage-7 metrics after the fix: hbond first/final/min/last20 `86.74/72.97/62.61/69.77`, Rg first/final/last20 `15.81/15.60/15.55 A`, within-chain final CA consecutive min/p50/max `3.23/3.77/4.31 A`, within-chain final CA(i)-CA(i+4) min/p50/max `4.65/6.29/12.30 A`.
- Production restraint audit: `/input/potential/restraint_position` is absent from the active stage-7 production file.
- Mapping provenance: `61` structure-geometry residues, `7` coil fallback residues, and `4` terminal charge overrides.

- [x] Phase 6: Document current CGL-CGL training/sampling path.
  - [x] Trace active generator and method documentation.
  - [x] Verify installed `parameters/dryMARTINI/dopc.h5` metadata.
  - [x] Answer the user's CGL-CGL rigidity and sampling questions without changing force-field files.

- [x] Phase 5: Resolve current visible `martini_1afo_hybrid` CGL gap / horizontal particle.
  - [x] Parse the active VTF and identify the exact displayed CGL rod(s) lying horizontally in the mid-gap.
  - [x] Map displayed VTF rod indices back to HDF5 CGL/CGLD indices and compare raw, centralized, and minimum-image direction vectors.
  - [x] Determine whether the artifact is caused by VTF centralization/wrapping, stale VTF output, or actual HDF5 trajectory geometry.
  - [x] Implement the smallest physical or visualization-only fix consistent with the cause.
  - [x] Regenerate the affected VTF/workflow output and verify no horizontal mid-gap CGL remains.

- [x] Phase 4: Follow-up visualization and 1RKL full-lipid helix audit.
  - [x] Compare CGL physical COM/vector geometry against VTF display head/tail construction.
  - [x] Inspect `martini_extract_vtf.py` display offsets/radii and determine whether the remaining visible gap is caused by rendering representation.
  - [x] Quantify 1RKL full-lipid protein hbond/Rg and local helix geometry over the reported/current trajectory.
  - [x] Decide whether code changes are needed; if so, keep them to visualization or physically justified preparation fixes only.
  - [x] Re-run targeted validation after any change.

- [x] Phase 3: Reproduce and fix current reported artifacts.
  - [x] Identify whether the reported artifacts are present in current outputs or are stale visualization/analysis results.
  - [x] Audit `example/16.MARTINI/run_sim_*.sh` and generated stage files for differences between CGL and full-resolution lipid workflows.
  - [x] Quantify bilayer gap, CGL orientation, leaflet separation, trapped CGL events, and protein secondary-structure / hbond retention from trajectories.
  - [x] Trace root cause to physical table generation, runtime injection, initial packing, workflow staging, or visualization/export.
  - [x] Implement the smallest physical fix and regenerate only required artifacts.
  - [x] Re-run targeted validations, then the affected `run_sim_*.sh` workflows.

- [x] Phase 1: Establish a stable CGL-only bilayer path using dynamic CGL vectors and production spline tables.
- [x] Phase 2A: Repair angular spline conditioning in production table builds.
  - [x] Restore CGL-CGL production call-site conditioning: use angular controls tied to sampled `cos_theta_count` and stronger smoothing.
  - [x] Apply the same conditioning principle to SC-CGL fitting and metadata.
  - [x] Enforce direct rotated-geometry sampling for CGL-CGL and SC-CGL production tables; hidden-bead relaxation is not part of the accepted model.
  - [x] Keep the change limited to fit conditioning; do not alter the physical sampled energy surface.
- [x] Phase 2B: Rebuild/check `dopc.h5` and verify focused CGL-CGL/SC-CGL table metadata and fit behavior.
  - [x] Reduced temporary PHE smoke build passed with direct rotated geometry, `n_angular=5`, smoothing `0.1`, and expected SC parameter shape.
  - [x] Production `parameters/dryMARTINI/dopc.h5` rebuilt with direct geometry metadata: CGL-CGL `n_angular=9`, `fit_smooth=0.1`, `fit_relax_steps=0`; SC-CGL `n_angular=11`, `angular_smooth=0.1`, `fit_relax_steps=0`.
  - [x] Restored CGL-CGL attraction retention, Boltzmann PMF averaging over unresolved azimuthal samples, and angular-resolved unresolved-core rows.
  - [x] Replaced the fixed CGL bead-frame table build with transferable two-body unresolved axial PMF averaging for CGL-CGL and CGL-SC.
  - [x] Increase CGL-CGL unresolved azimuthal sampling to `4^2`, rebuild/check installed tables, and rerun the focused bilayer smoke.
  - [x] Increase CGL axial bead-frame sampling to `8` for CGL-CGL/CGL-SC table generation, rebuild/check installed CGL-CGL, and rerun the focused bilayer smoke.
  - [x] Extend CGL-CGL radial support from resolved DOPC bead radius plus dry-MARTINI cutoff, rebuild/check installed CGL-CGL, and rerun the focused bilayer smoke.
  - [x] Add CGL initial leaflet-z de-overlap conditioning, then rerun the focused bilayer smoke.
  - [x] Use area-per-lipid same-leaflet XY spacing as the default coarse initial conditioning target, then rerun the focused bilayer smoke.
  - [x] Restore CGL-CGL tensor smoothing to the best-known `0.1` setting after the `0.5` trial worsened CGL-CGLD length, rebuild/check installed CGL-CGL, and rerun the focused bilayer smoke.
  - [x] Treat CGLD as a numerical unit-vector carrier by using explicit `2.0x` carrier-constraint stiffness over the projected DOPC bonded stiffness, then rerun the focused bilayer smoke.
  - [x] Keep CGL-CGL many-neighbor normalization deferred; focused validation passed with the transferable two-body table path.
- [x] Phase 2C: Extend the accepted two-body direct-geometry method to the full hybrid interaction surface and run full-system workflows.
  - [x] Audit/rebuild `dopc.h5`, `sidechain.h5`, and `particle.h5` so SC-CGL, CGL-particle, and SC-particle use physical direct dry-MARTINI geometry with unresolved bead-frame/azimuthal PMF averaging where applicable.
  - [x] Confirm no CGL-CGL/SC-CGL/CGL-particle/SC-particle table uses pair scaling, many-neighbor normalization, hidden-bead relaxation, physical-distance capping, or a standalone CGL orientation potential. CGL-CGL uses an explicit tempered two-body PMF over unresolved axial/azimuthal coordinates.
  - [x] Add or identify validation metrics for bilayer stability, CGL orientation, protein structural stability, and secondary-structure stability using existing Upside-native observables first.
  - [x] `example/16.MARTINI/run_sim_1rkl.sh`
  - [x] `example/16.MARTINI/run_sim_1afo.sh`
  - [x] `example/16.MARTINI/run_sim_1rkl_full.sh`
  - [x] `example/16.MARTINI/run_sim_1afo_full.sh`
  - [x] Implement extended-support full-tensor SC-CGL with transformed control and runtime inverse transform.
  - [x] Gate the new SC-CGL with H5 metadata audit, CGL-only bilayer smoke, and short 1RKL/1AFO CGL stage-6 diagnostics before launching full stage-7 workflows.
  - [x] Replace separable SC-particle factorization with a direct full radial-by-angular runtime table and rebuild `sidechain.h5`.
  - [x] Remove SC-env runtime force caps from the physical validation path.
  - [x] Remove active `0.10 nm` table-build pair-distance floors from SC-CGL and SC-particle, rebuild affected H5 artifacts, and rerun the acceptance tests.
  - [x] Re-run fresh full-resolution 1RKL/1AFO workflows after the SC-particle full-tensor fix and distance-floor audit.
- [x] Phase 2D: Rewrite `example/16.MARTINI/cg_lipid_potentials.tex` as a coherent JCTC methods section for the accepted physical table method.
  - [x] Remove obsolete descriptions of hidden-bead relaxation, distance floors, force caps, WCA/excluded-area projections, and staged debugging choices.
  - [x] Describe the accepted direct dry-MARTINI geometry, unresolved-coordinate averaging, log1p spline transforms, table ownership, and validation matrix.
  - [x] Verify the TeX no longer contradicts the installed H5 metadata or physical-model constraints.
- [x] Phase 2E: Re-audit CGL-CGL, SC-CGL, SC-particle, and CGL-particle physical integrity after the method rewrite.
  - [x] Confirm installed H5 tables have no twist coordinate, capping attributes, excluded-area projections, hidden relaxation, or interaction scaling.
  - [x] Confirm generated no-floor stage-7 `.up` files use unity protein-environment scale, zero SC-env force caps, zero generic Martini force cap, and no CGL orientation/twist nodes.

# Known Errors / Blockers

- Phase 3 reported artifacts are resolved in regenerated outputs. Root causes were preparation/staging defaults, not disabled interactions: median-COM CGL leaflet classification misclassified wrapped/tiled lipids near the bilayer center; the previous CGL z conditioning added a display-tail contact gap; full-lipid 1AFO used a stale protein-lipid packing clearance below the dry-MARTINI DOPC contact distance.
- Phase 5 current `martini_1afo_hybrid` artifact is resolved in the active output. The visible horizontal particle was VTF rod ordinal `104`, residue `130`, and mapped to HDF5 CGL/CGLD indices `464/736`; the defective output had final aligned-z `0.244` and was born during dynamics from a prepared cross-leaflet lateral near-overlap (`2.542 A`). CGL preparation now conditions opposite-leaflet lateral CGL spacing using the DOPC-derived `max_perp_radius_ang` envelope. Active regenerated `martini_1afo_hybrid` final frame: aligned-z min/p05/mean `0.700/0.881/0.954`, bad mid-gap HDF5 CGLs `[]`, bad mid-gap VTF rods `[]`, ordinal `104` aligned-z `0.890`, hbond first/final/min/last20 `85.87/85.02/74.68/80.78`, Rg first/final/last20 `15.26/15.38/15.47 A`.
- Phase 4 CGL gap follow-up points to visualization, not a physical bilayer void. Active `martini_1rkl_hybrid` final frame has CGL centerline tail overlap `-2.487 A` and resolved bead-envelope overlap `-10.715 A`; active `martini_1afo_hybrid` final frame has centerline tail overlap `-0.922 A` and bead-envelope overlap `-9.150 A`. The VTF exporter now emits CGL display radii from `max_perp_radius_ang=4.114 A`.
- Phase 4 1RKL full-lipid helix check does not show global de-folding. DSSP helix count is `14 -> 18` from production frame 0 to final, hbond first/final/min/last20 is `21.30/23.85/12.55/23.94`, and Rg first/final/last20 is `13.14/13.08/13.08 A`. The largest final local C-alpha i-to-i+4 stretch is around resseq `25-29` (`+3.34 A`), consistent with a small local C-terminal deformation rather than a workflow-wide loss of secondary structure.
- Fresh Phase 3 CGL validation passes. `1RKL` CGL final: tail gap `-0.513 A`, aligned-z min/p05/mean `0.599/0.873/0.955`, no bad-parallel CGLs, no flips, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.254/6.927 A`, protein last20 hbond/Rg `29.88/11.74`. `1AFO` CGL final: tail gap `-0.754 A`, aligned-z min/p05/mean `0.797/0.886/0.962`, no bad-parallel CGLs, no flips, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.331/6.909 A`, protein last20 hbond/Rg `74.29/15.45`.
- Fresh Phase 3 1AFO full-lipid validation improves the disrupted protein output. Old output production hbond first/final/min/last20 was `33.35/31.27/23.89/32.69` with Rg last20 `15.99`; regenerated output is `35.67/53.69/29.24/50.07` with Rg last20 `15.11`.
- Focused 50k CGL-only bilayer validation passes with the installed two-body path after runtime log-control injection, tempered CGL-CGL PMF, and no `0.10 nm` table-build distance floor: aligned-z min/p05/mean `0.953/0.966/0.988`, no flips/crossings, same-leaflet NN min/p05 `6.470/6.471 A`, CGL-CGLD RMS length deviation `0.105 A`.
- Fresh no-floor CGL-mode workflow validation passes. `1RKL` CGL: aligned-z min/p05/mean/finalmean `-0.374/0.860/0.942/0.939`, no leaflet crossings, final NN p05 lower/upper `6.890/6.755 A`, production last20 hbonds/Rg/total `28.54/12.13/-5653.6`, kinetic ratio `0.989`. `1AFO` CGL: aligned-z min/p05/mean/finalmean `0.798/0.913/0.971/0.971`, no leaflet crossings, final NN p05 lower/upper `7.458/6.967 A`, production last20 hbonds/Rg/total `70.19/14.66/-3218.9`, kinetic ratio `0.932`.
- Fresh no-floor full-resolution lipid workflow validation passes. `1RKL` full: DOPC leaflet separation `13.160 -> 13.055 A`, head-tail `|dz|` final p05/mean `8.095/13.249 A`, production last20 hbonds/Rg/total `25.08/13.44/-23906.9`, kinetic ratio `0.994`. `1AFO` full: DOPC leaflet separation `11.783 -> 11.613 A`, head-tail `|dz|` final p05/mean `7.975/13.386 A`, production last20 hbonds/Rg/total `54.24/14.30/-15685.6`, kinetic ratio `0.986`.

# Review

Phase 2A code-level implementation passed syntax checks, wrapper syntax checks, `git diff --check`, and a reduced temporary table build. Phase 2B focused CGL-only validation passes without CGL-CGL normalization. Phase 2C now passes the requested fresh validation matrix on 1RKL and 1AFO with both CGL and full-resolution lipid models using physical direct-geometry tables. Phase 2D rewrote the TeX method section around the accepted model and passed a local `pdflatex` compile check. Phase 2E re-audited the four requested interaction classes and found no active twist parameter, capping, arbitrary interaction scaling, or added CGL orientation potential in the installed tables or no-floor stage-7 files. Phase 3 reproduced the user-reported artifacts in current outputs, fixed the preparation/staging causes with physical defaults, and passed regenerated default-length 1RKL/1AFO CGL plus 1AFO full-lipid validations. Phase 4 fixed the CGL VTF display envelope and found no evidence that the current 1RKL full-lipid local helix deformation requires a physical-model change. Phase 5 corrected the remaining active `martini_1afo_hybrid` horizontal CGL by adding physical initial cross-leaflet CGL lateral de-overlap and replacing the active output with a validated regenerated run.
