# Progress Log

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
