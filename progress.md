# Progress Log

## 2026-06-10 Phase 5 Active 1AFO Horizontal CGL Fix
- Actions taken:
  - Reproduced the user's current visible `martini_1afo_hybrid` artifact in both VTF and HDF5. The horizontal rod is VTF ordinal `104`, residue `130`, mapping to HDF5 CGL/CGLD indices `464/736`.
  - Corrected the previous display-only diagnosis: the active output had an actual horizontal CGL, with final aligned-z `0.244` and bad mid-gap HDF5/VTF index `[104]`.
  - Traced the cause to the prepared coarse initial condition: the bad CGL started upright, but its nearest opposite-leaflet CGL was only `2.542 A` away laterally; the prepared 1AFO bilayer had cross-leaflet nearest-XY min/p05 `0.796/1.819 A`.
  - Added initial-condition-only cross-leaflet CGL lateral de-overlap using the DOPC-derived `max_perp_radius_ang` envelope. This does not add a runtime potential, cap, scaling, or orientation restraint.
  - Ran a short 1AFO CGL probe and a default-length 1AFO CGL validation under `outputs/phase5_1afo_crossxy_*`; both removed the bad mid-gap horizontal CGL.
  - Preserved the defective active output as `example/16.MARTINI/outputs/martini_1afo_hybrid.pre_crossxy_backup_20260610_154545`.
  - Replaced active `example/16.MARTINI/outputs/martini_1afo_hybrid` with the validated default-length regenerated run.
- Files modified:
  - `py/martini_prepare_system_lib.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python -m py_compile py/martini_prepare_system_lib.py py/martini_prepare_system.py py/martini_extract_vtf.py` passed before workflow validation.
  - Default validation `outputs/phase5_1afo_crossxy_final`: final aligned-z min/p05/mean `0.700/0.881/0.954`, bad mid-gap HDF5 CGLs `[]`, bad mid-gap VTF rods `[]`.
  - Active `outputs/martini_1afo_hybrid`: final aligned-z min/p05/mean `0.700/0.881/0.954`, bad mid-gap HDF5 CGLs `[]`, bad mid-gap VTF rods `[]`, ordinal `104` aligned-z `0.890`.
  - Active protein observables stayed stable: hbond first/final/min/last20 `85.87/85.02/74.68/80.78`, Rg first/final/last20 `15.26/15.38/15.47 A`.

## 2026-06-10 Phase 4 Visualization and 1RKL Full-Lipid Follow-up
- Actions taken:
  - Rechecked the active `outputs/martini_1rkl_hybrid` and `outputs/martini_1afo_hybrid` CGL trajectories that the VTFs are generated from.
  - Confirmed the apparent remaining CGL bilayer gap is a display issue: centerline tail endpoints already overlap and the resolved DOPC bead envelope overlaps substantially across the midplane.
  - Patched `py/martini_extract_vtf.py` so synthetic CGL display atoms emit a VTF `radius` equal to the canonical DOPC max perpendicular bead radius (`4.114 A`).
  - Regenerated `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf` and `example/16.MARTINI/outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf`.
  - Audited `outputs/martini_1rkl_hybrid_full/checkpoints/1rkl.stage_7.0.up` for the reported small alpha-helix de-folding.
- Files modified:
  - `py/martini_extract_vtf.py`
  - `plan.md`
  - `progress.md`
  - `findings.md`
- Verification:
  - `python -m py_compile py/martini_extract_vtf.py` passed.
  - Regenerated VTFs contain `radius 4.114` on CGL display atoms.
  - CGL final-frame metrics: 1RKL centerline tail gap `-2.487 A`, bead-envelope gap `-10.715 A`; 1AFO centerline tail gap `-0.922 A`, bead-envelope gap `-9.150 A`; neither has bad-parallel CGLs, flips, central COM lipids, or leaflet crossings.
  - 1RKL full-lipid protein metrics: DSSP helix count `14 -> 18` from production frame 0 to final, hbond first/final/min/last20 `21.30/23.85/12.55/23.94`, Rg first/final/last20 `13.14/13.08/13.08 A`; largest final local CA i-to-i+4 stretch is resseq `25-29` at `+3.34 A`.

## 2026-06-10 Phase 3 Reported Artifact Reproduction
- Actions taken:
  - Started a fresh debug pass for the user-reported CGL bilayer gap / trapped horizontal CGL and 1AFO full-lipid secondary-structure disruption.
  - Updated `plan.md` to make the contradiction between prior pass metrics and current observed artifacts the active blocker.
  - Reproduced the CGL artifact directly in `.up` files: current 1RKL CGL stage-7 has five bad CGL vectors, including three near-horizontal midplane lipids and two upper-leaflet classified lipids crossed into the lower leaflet; current 1AFO CGL has two crossed CGLs.
  - Found the 1RKL CGL issue is already present in stage-6 prepared input for four tiled corner lipids because median-COM leaflet classification misclassifies lower-leaflet lipids whose COMs sit near the bilayer center.
  - Reproduced the 1AFO full-lipid protein disruption in stage-7 burn-in: the hbond score starts near 92 after stage-7 minimization but drops to roughly 33 before production output.
  - Patched CGL preparation to classify leaflets by CGL direction when both signs are present and to stop default z-expansion of CGL COM planes beyond the resolved DOPC geometry.
  - Patched hybrid workflow defaults so protein-lipid packing clearance is derived from dry-MARTINI DOPC max contact distance unless explicitly overridden.
  - First fresh CGL rerun showed the no-z-expansion variant starts from too-overlapped CGL COM planes and leaves horizontal CGL vectors after minimization/production. Revised the default z conditioning to the derived zero-gap display-tail spacing `2 * tail_projection_ang`, replacing the old gap-producing `2 * tail_projection_ang + contact`.
  - Zero-gap short CGL probes passed the CGL artifact checks:
    - 1RKL final `aligned_z min/p05/mean=0.701/0.849/0.955`, `bad_parallel=0`, `bad_flip=0`, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.818/7.190 A`.
    - 1AFO final `aligned_z min/p05/mean=0.790/0.876/0.959`, `bad_parallel=0`, `bad_flip=0`, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.994/7.206 A`.
  - Default-length regenerated CGL runs passed the reported artifact checks:
    - 1RKL final `tail_gap=-0.513 A`, `aligned_z min/p05/mean=0.599/0.873/0.955`, `bad_parallel=0`, `bad_flip=0`, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.254/6.927 A`, protein last20 hbond/Rg `29.88/11.74`.
    - 1AFO final `tail_gap=-0.754 A`, `aligned_z min/p05/mean=0.797/0.886/0.962`, `bad_parallel=0`, `bad_flip=0`, no central COM lipids, no leaflet crossings, same-leaflet NN min/p05 `6.331/6.909 A`, protein last20 hbond/Rg `74.29/15.45`.
  - Regenerated 1AFO full-lipid run improved secondary-structure retention relative to the old output: hbond first/final/min/last20 changed from `33.35/31.27/23.89/32.69` to `35.67/53.69/29.24/50.07`; Rg last20 changed from `15.99` to `15.11`.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/run_sim_hybrid.sh`
- Verification:
  - `python -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed.
  - `bash -n example/16.MARTINI/run_sim_hybrid.sh` passed.
  - Fresh default-length workflows completed: `run_sim_1rkl.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_full.sh` under `example/16.MARTINI/outputs/phase3_*`.

## 2026-06-10 Phase 2E Physical Integrity Re-audit
- Actions taken:
  - Re-audited installed `dopc.h5` and `sidechain.h5` metadata for CGL-CGL, SC-CGL, CGL-particle, and SC-particle.
  - Searched active table-generation, preparation, and runtime code for twist/orientation nodes, cap attributes, force caps, excluded-area projections, normalization, and interface scales.
  - Audited the four no-floor stage-7 `.up` files from the accepted 1RKL/1AFO CGL/full validation matrix.
- Verification:
  - All four table classes passed assertions for `fit_relax_steps=0`, `sample_dist_min_nm=1e-6`, numerical-zero-guard metadata, no excluded-area rows/projections, no cap attributes, and direct retained dry-MARTINI source metadata.
  - CGL-CGL and SC-CGL use full tensors with `n_modes=0`; unresolved axial/bead-frame sampling is table-generation quadrature, not a runtime twist coordinate.
  - Stage-7 files have `protein_env_interface_scale=1.0`, `sc_env_lj_force_cap=0.0`, `sc_env_coul_force_cap=0.0`, generic Martini `force_cap=0.0`, and no CGL leaflet/orientation/twist nodes.
  - Caveat: CGLD still has the documented numerical carrier bond and its `2x` stiffness metadata. That is separate from the four audited nonbonded table interactions and is not an added CGL orientation potential.

## 2026-06-10 Phase 2D JCTC Method Rewrite
- Actions taken:
  - Rewrote `example/16.MARTINI/cg_lipid_potentials.tex` as a coherent methods section for the accepted direct dry-MARTINI geometry model.
  - Removed obsolete method blocks describing hidden-bead relaxation, finite distance floors, WCA/excluded-area projections, force caps, separable SC-CGL/SC-particle assumptions, and staged debugging choices.
  - Documented CGL-CGL, SC-CGL, CGL-particle, and SC-particle table construction, unresolved-coordinate averaging, log1p inverse transforms, runtime ownership, units, and validation criteria.
- Files modified:
  - `example/16.MARTINI/cg_lipid_potentials.tex`
  - `plan.md`
- Verification:
  - H5 metadata was rechecked against `parameters/dryMARTINI/dopc.h5` and `parameters/dryMARTINI/sidechain.h5`.
  - Stale-term scan found only explicit negative statements such as no twisting coordinate, no normalization, no nonnegative projection, and no force cap.
  - `pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/upside2_tex_check example/16.MARTINI/cg_lipid_potentials.tex` passed and produced a 6-page PDF.

## 2026-06-10 Phase 2C Direct-Geometry Hybrid Acceptance
- Actions taken:
  - Rebuilt `sidechain.h5` and `dopc.h5` after removing active `0.10 nm` dry-MARTINI pair-distance floors from SC-particle and SC-CGL table generation. Active tables now use only a `1e-6 nm` numerical singularity guard.
  - Added the CGL-target `log1p((E - E_ref) / kBT)` reduced-PMF spline control and runtime inverse transform so charged target hard cores remain physical without float32 Boltzmann-weight underflow.
  - Added the same runtime inverse-transform support for CGL-CGL pair tables and fixed `inject_cg_lipid_nodes()` so generated stage `.up` files copy `log1p_reduced_transform`, `boltzmann_temperature_upside`, `energy_transform`, `spline_control_quantity`, and `reference_energy_eup`.
  - Added radial low-end extrapolation in the CGL spline evaluators to remove the zero-derivative plateau below the first radial knot.
  - Replaced the 4-mode separable SC-CGL path with an extended-support full tensor and transformed control.
  - Replaced separable SC-particle runtime use with a direct full radial-by-angular rotamer table and removed SC-env runtime force caps from the physical validation path.
  - Ran the requested fresh validation matrix after the no-floor rebuild: 1RKL/1AFO with CGL lipids and 1RKL/1AFO with full-resolution lipids.
- Files modified:
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `src/martini_cg_lipid.cpp`
  - `src/martini_internal.h`
  - `src/martini_potential.cpp`
  - `parameters/dryMARTINI/sidechain.h5`
  - `parameters/dryMARTINI/dopc.h5`
- Verification:
  - `python -m py_compile py/martini_build_tables.py py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed.
  - `make -C obj -j2` passed.
  - Production H5 metadata audit confirmed direct geometry, `fit_relax_steps=0`, no pair scaling, no many-neighbor normalization, no hidden relaxation, no standalone CGL orientation potential, `sample_dist_min_nm=1e-6`, CGL-target `log1p_reduced_transform=1`, and CGL-CGL `log1p_reduced_transform=1`.
  - Runtime stage-7 `.up` audit confirmed CGL-CGL pair transform metadata and `pair_interaction/reference_energy_eup` are present after injection.
  - Focused CGL-only 50k bilayer validation passed: aligned-z min/p05/mean `0.953/0.966/0.988`, no flips/crossings, same-leaflet NN min/p05 `6.470/6.471 A`, CGL-CGLD RMS length deviation `0.105 A`.
  - 1RKL CGL completed: aligned-z min/p05/mean/finalmean `-0.374/0.860/0.942/0.939`, no leaflet crossings, final NN p05 lower/upper `6.890/6.755 A`; protein production last20 hbonds/Rg/total `28.54/12.13/-5653.6`; kinetic ratio `0.989`.
  - 1AFO CGL completed: aligned-z min/p05/mean/finalmean `0.798/0.913/0.971/0.971`, no leaflet crossings, final NN p05 lower/upper `7.458/6.967 A`; protein production last20 hbonds/Rg/total `70.19/14.66/-3218.9`; kinetic ratio `0.932`.
  - 1RKL full-resolution completed: DOPC leaflet separation `13.160 -> 13.055 A`, head-tail `|dz|` final p05/mean `8.095/13.249 A`; protein production last20 hbonds/Rg/total `25.08/13.44/-23906.9`; kinetic ratio `0.994`.
  - 1AFO full-resolution completed: DOPC leaflet separation `11.783 -> 11.613 A`, head-tail `|dz|` final p05/mean `7.975/13.386 A`; protein production last20 hbonds/Rg/total `54.24/14.30/-15685.6`; kinetic ratio `0.986`.

## 2026-06-09 CGL-CGL Table Conditioning and Bilayer Smoke
- Actions taken:
  - Restored production table-fit conditioning in `py/martini_build_tables.py`: CGL-CGL angular controls now follow the actual sampled angular grid (`7 -> 9`), SC-CGL uses `9 -> 11`, and both use stronger angular smoothing.
  - Enforced direct rotated-geometry CGL table generation by making `UPSIDE_MARTINI_FIT_RELAX_STEPS=0` valid/default and rejecting hidden-bead relaxation for CGL-CGL/SC-CGL.
  - Rebuilt `parameters/dryMARTINI/dopc.h5` atomically with direct-geometry metadata.
  - Fixed CGL-CGL table physics inconsistencies: retained resolved dry-MARTINI attractions, implemented Boltzmann free-energy averaging for unresolved CGL-CGL azimuthal samples, and changed unresolved short-range rows to angular-resolved first-shell controls.
  - Extended CGL-CGL radial support to `42 A`, increased unresolved CGL bead-frame sampling to `8`, restored extended-table smoothing to `0.1` after a failed `0.5` trial, and kept normalization deferred.
  - Added initial CGL leaflet z de-overlap and area-derived same-leaflet XY conditioning for coarse bilayer preparation.
  - Treated CGLD as a numerical unit-vector carrier by using explicit `2.0x` carrier stiffness over the projected DOPC bonded stiffness while retaining projected stiffness metadata.
- Files modified:
  - `py/martini_build_tables.py`
  - `py/martini_cg_lipid_params.py`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/build_martini_h5_m1.sh`
  - `example/16.MARTINI/build_martini_h5_slurm.sh`
  - `example/16.MARTINI/build_martini_h5_m1_temp.sh`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
  - `parameters/dryMARTINI/dopc.h5`
  - `plan.md`
  - `progress.md`
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_gen_params.py` passed.
  - `bash -n` passed for the M1, Slurm, and temporary M1 MARTINI H5 rebuild scripts.
  - `git diff --check` passed before the production rebuild.
  - Reduced temporary PHE table build passed with `fit_relax_steps=0`, expected CGL-CGL/SC-CGL metadata, and expected SC parameter shape.
  - Installed `dopc.h5` CGL-CGL metadata now has `fit_smooth=0.1`, `fit_r_max_nm=4.2`, `n_radial=32`, `cutoff_ang=42.0`, `cgl_bead_frame_count=8`, and `azimuthal_average=boltzmann_free_energy`.
  - Focused installed-table 200-step CGL-only bilayer smoke passed: no flips or leaflet crossings, aligned-z min/p05/mean `0.527/0.786/0.930`, same-leaflet NN min/p05 `4.369/5.230 A`, CGL-CGLD length min/max/rmsdev `10.679/13.475/0.648 A`.
- Current blocker:
  - Full protein+lipid workflows remain unverified after the focused CGL-only pass.

## 2026-05-29 Stage 6.0 Minimization Overshoot Fix
- Problem: 1rkl hybrid system's stage 6.0 minimization overshoots from ~350k E_up to ~-90k E_up
  because the protein (in rigid groups) drifts into the lipid beads, finding deep LJ attractive
  wells. 1afo works fine because it has fewer lipid beads relative to its larger protein.
- Fix: Inject harmonic position restraints (k=10.0 E_up/A^2) on all protein atoms before stage
  6.0 minimization. The rigid body solver converts per-atom restraint forces to COM + orientation
  restraints, keeping the protein near its initial position during minimization. Restraints are
  removed before MD.
- Verified:
  - 1rkl: min 409k→+5.5k (no overshoot), MD Rg=12.8 stable, stage 7.0 Rg=12.7 stable
  - 1afo: min 791k→+33.6k (no overshoot), MD Rg=15.9 stable, stage 7.0 Rg=15.8 stable
- Files modified:
  - `py/martini_prepare_system.py` — added `inject_protein_position_restraints`,
    `remove_protein_position_restraints`; integrated into stage 6.0 workflow

## 2026-05-28 CG-CG B-Spline Angular Regularization Fix
- Actions taken:
  - Diagnosed root cause of messy bilayer orientation: CG-CG B-spline had severe
    angular underdetermination (15×15=225 controls fitting 7×7=49 data points per
    radial distance). With Tikhonov λ=0.01, the fit created spurious angular
    oscillations between sample points, producing unphysical orientational energy
    features not present in the underlying force field.
  - Added `n_knot_angular` and `cg_smooth` parameters to `_fit_cg_lipid_quadspline`.
  - Changed call site in `_build_cg_lipid_tables` to use
    `n_knot_angular=min(_cg_ct + 2, 15) = 9` and `cg_smooth=0.1`.
  - Rebuilt `dopc.h5` cg_lipid_pair group: 14×9² = 1134 params (was 3150).
  - Bilayer test PASSES: validation=PASS, no flips, no crossings.
  - Angular fine-grid std reduced 62-70% across all radial distances.
  - Updated `cg_lipid_potentials.tex` with new N_θ=9 and regularization rationale.
- Files modified:
  - `py/martini_build_tables.py` — added params, updated call site
  - `example/16.MARTINI/build_martini_h5_m1_temp.sh` — updated call for rebuild
  - `example/16.MARTINI/cg_lipid_potentials.tex` — N_θ 15→9, added regularization section
  - `parameters/dryMARTINI/dopc.h5` — regenerated cg_lipid_pair
  - `plan.md` — updated with diagnosis and fix
  - `findings.md` — added B-spline underdetermination lesson
  - `progress.md` — this entry
- Remaining:
  - Full protein+lipid workflow re-run needed to verify orientation improvement.
  - SC-CGL table has similar (less severe) underdetermination; not yet addressed.

## 2026-05-27 Full-Lipid/CGL Correction
- Actions taken:
  - Removed the standalone `cg_lipid_leaflet_orientation` runtime potential and
    injection path.
  - Changed CGL--CGL table generation so `dopc.h5` retains the resolved
    dry-MARTINI lipid--lipid attraction instead of subtracting/clipping the
    attractive background and adding a separate orientation correction.
  - Reworked SC-particle factorization for full-resolution lipid mode: signed
    SVD is retained in the resolved range, but rows are floored against sampled
    dry-MARTINI minima and unresolved hard-core rows are finite radial barriers.
  - Updated `cg_lipid_potentials.tex` and lessons to reflect the user's
    correction.
  - Removed hidden-bead relaxation from CGL--CGL and SC--CGL table generation;
    tables now use direct rotated full-resolution bead geometries at sampled
    direction vectors.
  - Made MARTINI `.h5` builders write through sibling temp files and atomically
    replace completed outputs.
  - Regenerated `parameters/dryMARTINI/dopc.h5` and `sidechain.h5`.
- Verification:
  - Python compile and C++ build passed after removing the CGL orientation node.
  - Temporary rebuilt SC table:
    reconstruction min/max `-17.31/1.04e6 kJ/mol`, no huge negative wells.
  - CGL--CGL smoke table retained attractive controls in `dopc.h5` path:
    sampled raw min/max `-30.67/1.59e6 kJ/mol`, fitted controls
    `-10.52/5.45e5 E_up`, `attractive_control_source=retained_full_resolved_dry_martini_pair_table`.
  - Final `dopc.h5` check: CGL--CGL `fit_relax_steps=0`,
    `isotropic_background_source=none_full_resolved_dry_martini_pair_table`,
    `excluded_area_nonnegative_rows=0`, retained attractive controls
    (`count=434`), and SC--CGL `fit_relax_steps=0`.
  - Final `sidechain.h5` check: SC-particle reconstruction min/max
    `-15.68/9.62e5 kJ/mol`; no huge negative lipid-sidechain well remains.
  - Current stage files keep termini charged: 1AFO fragments are
    `(0,36) Qd->Qa` and `(36,72) Qd->Qa`; 1RKL is `(0,31) Qd->Qa`, in both
    CGLipid and full-resolution outputs.
  - `python3 -m py_compile ...`, `cmake --build obj`, and `git diff --check`
    passed.
- Failure and fix:
  - A copied trajectory test from `stage_7.0.prepared.up` was invalid because
    that file is pre-minimization and has severe full-lipid clashes; the run was
    stopped after NaNs. Use regenerated workflow handoff, not prepared files,
    for trajectory validation.
  - A sandboxed DOPC rebuild fell back to threads and was too slow; rerunning
    outside the sandbox used process workers and completed the regenerated
    `dopc.h5`.

## 2026-05-26 Stage-7 1AFO/1RKL Debug
- Actions taken:
  - Compared stage-7 prepared inputs, promoted production inputs, and output
    frame 0 for `martini_1afo_hybrid`, `martini_1afo_hybrid_full`,
    `martini_1rkl_hybrid`, and `martini_1rkl_hybrid_full`.
  - Identified a force-field bug: an unphysical negative SC-particle
    factorization well affecting full-resolution lipid burn-in.
  - An intermediate `cg_lipid_leaflet_orientation` attempt was removed on
    2026-05-27; CGL orientation behavior now belongs in the `dopc.h5`
    CGL--CGL spline table.
  - Added thread fallback for table builds when process pools are disallowed.
  - Updated `example/16.MARTINI/cg_lipid_potentials.tex`.
- Files modified:
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system_lib.py`
  - `src/martini_cg_lipid.cpp`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
  - `findings.md`
  - `progress.md`
- Test results:
  - Existing stage-7 prepared protein coordinates match reference mapping
    (`RMSD=0.000 A`); frame-0 bends come from burn-in-promoted input state.
  - Temporary rebuilt SC-particle table had runtime minimum `-17.3 kJ/mol`,
    versus `-4.7e11 kJ/mol` from the stale table.

## 2026-05-26 MARTINI H5 Rebuild Scripts
- Actions taken:
  - Added local M1 rebuild script for all dry-MARTINI `.h5` files.
  - Added Slurm rebuild script for all dry-MARTINI `.h5` files.
  - Both scripts call `py/martini_gen_params.py --force --upside-home <repo>`
    so outputs are written under `parameters/dryMARTINI`.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `example/16.MARTINI/build_martini_h5_m1.sh`
  - `example/16.MARTINI/build_martini_h5_slurm.sh`
- Test results:
  - `bash -n example/16.MARTINI/build_martini_h5_m1.sh` passed.
  - `bash -n example/16.MARTINI/build_martini_h5_slurm.sh` passed.
  - Both scripts are executable.
- Not run:
  - Full table regeneration, because it is long-running and rewrites
    production parameter files.

## 2026-05-26 Direction-Vector MARTINI Table Builds
- Actions taken:
  - Audited SC-particle, CGL-particle, SC-CGL, and CGL-CGL table construction
    for direction-vector sampling completeness.
  - Reopened the previous around-vector wording after the user clarified that
    direction vectors are the intended potential inputs.
  - Changed optional around-vector bead-frame quadrature to default to one
    sample and renamed the exposed controls/metadata to bead-frame terminology.
  - Updated `example/16.MARTINI/cg_lipid_potentials.tex` so CGL-CGL, SC-CGL,
    CGL-particle, and SC-particle sections describe direction-vector sampling
    as the physical spline input.
  - Added multiprocessing across independent table slices, using
    `UPSIDE_MARTINI_TABLE_WORKERS` first, then Slurm CPU allocation variables,
    then local CPU count.
  - Confirmed `py/martini_gen_params.py --help` exposes bead-frame controls
    rather than the previous around-vector control wording.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `py/martini_build_tables.py`
  - `py/martini_gen_params.py`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
- Test results:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_gen_params.py py/martini_prepare_system_lib.py py/martini_prepare_system.py example/16.MARTINI/test_cg_bilayer/run_test.py` passed.
  - Focused CGL-CGL smoke generated a `3 x 3 x 3` raw grid with
    `azimuthal_count=2` and `bead_frame_count=2`.
  - Focused CGL-particle smoke generated a C1/Qa target table with
    `orientation_sampling=cgl_direction_vector` and `cgl_bead_frame_count=2`.
  - Focused SC-particle smoke generated a PHE x C1 table with
    `orientation_sampling=target_direction_vector_grid` and
    `sidechain_bead_frame_count=2`.
  - Focused PHE SC-CGL smoke generated parameters with
    `sidechain_bead_frame_count=2` and `cg_bead_frame_count=2`.
  - `git diff --check` passed.
- Failures and fixes:
  - The local sandbox blocks process-pool semaphore queries, so smoke tests
    printed the intended worker count and fell back to one worker. The code
    keeps multiprocessing enabled for normal M1/Slurm environments.

## 2026-05-25 Orientation-Resolved MARTINI CG Tables
- Actions taken:
  - Re-scoped the active task around the corrected physical requirement:
    CG lipid and CG sidechain spline tables must be generated by rotating
    resolved full-resolution bead models over sampled orientations.
  - Documented that runtime CG lipid geometry must not be derived from one
    packed lipid conformation.
  - Changed coarse runtime CGLD geometry to derive from canonical
    `parameters/dryMARTINI/DOPC.pdb` rather than the first packed lipid.
  - Removed the extra `cg_lipid_leaflet_orientation` potential; CGL
    orientational forces are owned by the orientation-dependent spline tables.
  - Changed SC-particle and SC-CGL table fitting to expand each rotamer
    center/vector into resolved MARTINI sidechain bead positions from MARTINI
    bonded geometry.
  - Corrected the model notes: SC-particle is shared in both lipid modes for
    non-CGL environment particles; full-resolution mode also uses it for
    explicit lipid beads, while CGLipid mode replaces only explicit-lipid
    SC-particle contacts with SC-CGL. There is no dry-MARTINI SC-SC table.
  - Updated the CGL-only bilayer test harness to use installed production
    `particle.h5` plus `dopc.h5` by default, while keeping `--rebuild-tables`
    for explicit local table refits.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `example/16.MARTINI/test_cg_bilayer/run_test.py`
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system_lib.py`
  - `src/martini_cg_lipid.cpp`
  - `src/main.cpp`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
- Test results:
  - `python3 -m py_compile example/16.MARTINI/test_cg_bilayer/run_test.py py/martini_build_tables.py py/martini_prepare_system_lib.py py/martini_prepare_system.py` passed.
  - Canonical DOPC runtime geometry check reproduced current table metadata:
    `orientation_length_ang=11.139272`, `orientation_mass_g_mol=77.048875`,
    `orientation_bond_fc_eup_a2=39.435978`; canonical display offsets are
    `head=-14.558480 A`, `tail=11.139272 A`.
  - Focused PHE SC-env table smoke test generated a resolved 3-bead PHE table.
  - Focused PHE SC-CGL smoke test exercised the two-azimuth resolved-sidechain
    fitting path.
  - CGL-only bilayer validation with installed production tables passed:
    200-step default smoke and 2000-step NVT run. The 2000-step run had no
    flips, no leaflet crossings, same-leaflet nearest-neighbor min/p05
    `6.304/6.483 A`, and CGL-CGLD length min/max/rmsdev
    `10.732/11.515/0.161 A`.
  - Existing stale 1AFO coarse output fails the new compose/table geometry
    guard, as expected, until regenerated.
  - `cmake --build obj` passed.
- Failures and fixes:
  - Initial full-sidechain smoke test failed because the Upside rotamer library
    has one placement bead per rotamer. The fix now derives MARTINI SC bead
    offsets from `martinize.py` and places them at each rotamer center/vector.
  - The old CGL-only bilayer harness imported removed helper symbols and tried
    to refit tables by default. It now reads DOPC bead types from the ITP parser,
    skips obsolete debug output calls, and defaults to installed production
    tables for runtime validation.
