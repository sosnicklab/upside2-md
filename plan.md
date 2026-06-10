# Project Goal

Build a physically defensible single-vector DOPC coarse-grained lipid (CGL) force field for the UPSIDE/dry-MARTINI hybrid workflow. The current active task is to debug the user-observed CGL bilayer gap / horizontally trapped CGL artifact in `martini_1rkl_hybrid` and `martini_1afo_hybrid`, plus secondary-structure disruption in `martini_1afo_hybrid_full`, without disabling SC-env/BB-env interactions, adding ad-hoc CGL orientation potentials, capping forces, arbitrary scaling, or twisting parameters.

# Architecture & Key Decisions

- Runtime CGL orientation must come from the vector-particle spline interactions, especially CGL-CGL and SC-CGL. Do not add a standalone CGL orientation potential.
- Table generation must evaluate direct dry-MARTINI nonbonded energies from rotated full-resolution bead geometries over the sampled runtime direction-vector coordinates.
- CGL-CGL tables must retain the resolved dry-MARTINI lipid-lipid attraction in the spline. Do not subtract it into an absent radial background or replace it with a separate orientation correction.
- Revised Decision: CGL-CGL unresolved azimuthal/bead-frame samples should use a tempered direct PMF as the next two-body model. After fixing runtime log-control injection, a 50k bilayer-only validation with the production-temperature Boltzmann pair PMF still collapsed same-leaflet packing (`p05 ~3.86 A`) and produced a flip, so the pair PMF remains too favorable for independently chosen hidden axial rotations in a dense single-vector bilayer. Direct energy expectation was too stiff and launched the bilayer apart. A tempered PMF keeps direct dry-MARTINI energies and the three runtime coordinates while avoiding capping, normalization, hidden relaxation, or an added orientation potential. SC-CGL/CGL-target remain production-temperature Boltzmann PMFs unless their own validation shows the same dense-phase representability failure.
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

# Execution Phases

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
- Fresh Phase 3 CGL validation passes. `1RKL` CGL final: tail gap `-0.513 A`, aligned-z min/p05/mean `0.599/0.873/0.955`, no bad-parallel CGLs, no flips, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.254/6.927 A`, protein last20 hbond/Rg `29.88/11.74`. `1AFO` CGL final: tail gap `-0.754 A`, aligned-z min/p05/mean `0.797/0.886/0.962`, no bad-parallel CGLs, no flips, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.331/6.909 A`, protein last20 hbond/Rg `74.29/15.45`.
- Fresh Phase 3 1AFO full-lipid validation improves the disrupted protein output. Old output production hbond first/final/min/last20 was `33.35/31.27/23.89/32.69` with Rg last20 `15.99`; regenerated output is `35.67/53.69/29.24/50.07` with Rg last20 `15.11`.
- Focused 50k CGL-only bilayer validation passes with the installed two-body path after runtime log-control injection, tempered CGL-CGL PMF, and no `0.10 nm` table-build distance floor: aligned-z min/p05/mean `0.953/0.966/0.988`, no flips/crossings, same-leaflet NN min/p05 `6.470/6.471 A`, CGL-CGLD RMS length deviation `0.105 A`.
- Fresh no-floor CGL-mode workflow validation passes. `1RKL` CGL: aligned-z min/p05/mean/finalmean `-0.374/0.860/0.942/0.939`, no leaflet crossings, final NN p05 lower/upper `6.890/6.755 A`, production last20 hbonds/Rg/total `28.54/12.13/-5653.6`, kinetic ratio `0.989`. `1AFO` CGL: aligned-z min/p05/mean/finalmean `0.798/0.913/0.971/0.971`, no leaflet crossings, final NN p05 lower/upper `7.458/6.967 A`, production last20 hbonds/Rg/total `70.19/14.66/-3218.9`, kinetic ratio `0.932`.
- Fresh no-floor full-resolution lipid workflow validation passes. `1RKL` full: DOPC leaflet separation `13.160 -> 13.055 A`, head-tail `|dz|` final p05/mean `8.095/13.249 A`, production last20 hbonds/Rg/total `25.08/13.44/-23906.9`, kinetic ratio `0.994`. `1AFO` full: DOPC leaflet separation `11.783 -> 11.613 A`, head-tail `|dz|` final p05/mean `7.975/13.386 A`, production last20 hbonds/Rg/total `54.24/14.30/-15685.6`, kinetic ratio `0.986`.

# Review

Phase 2A code-level implementation passed syntax checks, wrapper syntax checks, `git diff --check`, and a reduced temporary table build. Phase 2B focused CGL-only validation passes without CGL-CGL normalization. Phase 2C now passes the requested fresh validation matrix on 1RKL and 1AFO with both CGL and full-resolution lipid models using physical direct-geometry tables. Phase 2D rewrote the TeX method section around the accepted model and passed a local `pdflatex` compile check. Phase 2E re-audited the four requested interaction classes and found no active twist parameter, capping, arbitrary interaction scaling, or added CGL orientation potential in the installed tables or no-floor stage-7 files. Phase 3 reproduced the user-reported artifacts in current outputs, fixed the preparation/staging causes with physical defaults, and passed regenerated default-length 1RKL/1AFO CGL plus 1AFO full-lipid validations.
