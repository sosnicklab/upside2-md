## 2026-05-20 1RKL Replica-2 Shape and Box Drift Debug

### Project Goal
- Debug `example/16.MARTINI/run_sim_1rkl.sh` using the exact artifacts in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/`.
- Explain why `1rkl.stage_7.0.vtf` has residual secondary-structure drift, why `1rkl.stage_7.2.vtf` loses protein shape, why the bilayer becomes elliptical in XY, and why box `total_potential` keeps decreasing.
- Validate fixes by modifying copied `.h5/.up` artifacts from the reported output before any long table-generation workflow rerun.
- Preserve all physical interface terms: SC-env, BB-env, CGL-SC, CGL-target, CGL-CGL, generic dry-MARTINI, ion, and bonded interactions remain active.
- Keep `example/16.MARTINI/cg_lipid_potentials.tex` synchronized with any physical model change.

### Architecture & Key Decisions
- Start from the exact named output and compare all available stage-7 replicas, not only trajectory 0.
- Treat monotonic box-potential lowering and XY ellipse formation as likely barostat/box-coupling or pressure-unit symptoms until metrics prove otherwise.
- Separate protein secondary-structure loss from membrane morphology: measure hbond retention, CA RMSD/Rg, bilayer XY covariance/aspect ratio, box vectors, lipid nearest distances, and component energies.
- Finding: stage 7 uses a fixed square box after NPT handoff, so the XY ellipse is coordinate/membrane morphology drift, not anisotropic box-vector drift.
- Finding: the total-potential drift is dominated by `cg_lipid_pair`; the number of attractive CGL-CGL pairs grows from about 2,794 at stage-7 start to over 6,200 by the end of `7.2`.
- Finding: copied-HDF5 rigid-protein relaxation still drives `cg_lipid_pair` from about `-90.9k` to `-116.7k E_up` in 10k steps, so production release timing is not the root cause.
- Revised Decision: explicit DOPC-derived CGL table sampling must use the same dry-MARTINI bead-level nonbonded cutoff as the generic MARTINI runtime. The current DOPC projection sums bead LJ/Coulomb tails beyond the generic `12 A` cutoff, creating a many-neighbor attraction that the underlying bead model would not apply.
- Finding: the bead-cutoff-only copied-HDF5 patch still descends rapidly, so the remaining failure is the isotropic cohesive background in isolated DOPC-DOPC sampling being reused as a many-neighbor additive pair potential.
- Revised Decision: CGL-CGL fitting should retain dry-MARTINI excluded volume and orientation-dependent residuals, but subtract the attractive radial angular mean from each sampled row. Positive/repulsive row means remain unchanged. This removes only the non-transferable isotropic cohesive background; it does not disable CGL-CGL interactions.
- Finding: radial-background subtraction alone still leaves attractive spline controls that act as a many-neighbor energy sink. A copied-HDF5 validation with all attractive CGL-CGL controls removed keeps `cg_lipid_pair` near a small positive value and preserves protein Rg near `12.1-12.4 A` over 10k steps.
- Revised Decision: CGL-CGL owns excluded volume and orientation-dependent repulsive compatibility for the single-particle membrane; attractive isolated-DOPC pair controls are not transferable to a many-neighbor additive CGL model and must be removed from the CGL-CGL spline controls. Protein-lipid and target interactions remain active and can still be attractive where physically resolved.
- If the instability comes from barostat or handoff semantics, patch those mechanics without altering force-field parameters.
- If an interaction table is still physically wrong, patch table construction and guard stale artifacts; do not introduce empirical scales, caps, exclusions, or restraints.

### Execution Phases
- [x] Phase 1: Inventory exact output artifacts, stage files, logs, and current source diffs.
- [x] Phase 2: Quantify per-replica protein, bilayer, box, and component-energy drift for stage 7.
- [x] Phase 3: Identify the root cause path and decide whether it is table semantics, force routing, barostat, or restart/handoff.
- [x] Phase 4: Patch copied HDF5 artifacts from `martini_test_1rkl_hybrid` for focused validation.
- [x] Phase 5: Implement minimal source/doc changes and update `cg_lipid_potentials.tex` if the physical model changes.
- [x] Phase 6: Run focused validation and document final results.

### Known Errors / Blockers
- The existing reported `martini_test_1rkl_hybrid/martini.h5` is now intentionally stale and rejected by the schema guard. Rebuild or regenerate tables before rerunning the workflow normally.

### Review
- Root cause: the exact output uses a CGL-CGL table that carries isolated DOPC-DOPC attractive controls into an additive many-neighbor single-particle membrane. The membrane lowers energy by creating more attractive CGL-CGL neighbors; CGL-CGL interacting pairs grow from about `2,794` at stage-7 start to over `6,200` by the end of `7.2`, producing XY morphology drift and later protein deformation.
- Implemented physical corrections:
  - Explicit DOPC-derived CG-lipid table sampling now records and uses the generic dry-MARTINI bead-level nonbonded cutoff (`1.2 nm`).
  - CGL-CGL table fitting subtracts the attractive radial angular mean from sampled rows and removes remaining attractive CGL-CGL spline controls, leaving CGL-CGL as excluded-volume/orientation-compatibility physics for the one-particle membrane.
  - Stale-table validation rejects old `martini.h5` files that lack the cutoff/background/attractive-control metadata.
  - `cg_lipid_potentials.tex` documents the cutoff, background subtraction, and non-transferable CGL-CGL attraction removal.
- Copied-HDF5 validation:
  - Rigid-protein validation and bead-cutoff-only validation both failed; they still showed rapid CGL-CGL energy descent.
  - The accepted copied-HDF5 validation in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid_repcheck` ran `10k` total stage-7 steps from the reported starting HDF5 with patched CGL-CGL controls.
  - Final CGL nearest-neighbor min/p05/mean was `6.278/6.729/7.821 Å`; XY aspect stayed about `1.08`.
  - Final protein Rg was `12.47 Å`, hbond sum `26.75`, and aligned CA RMSD `4.16 Å`.
- Verification:
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py` passed.
  - Stale guard rejects the reported old `martini.h5`.
  - `git diff --check` passed.
