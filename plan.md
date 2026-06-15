# Project Goal

Improve the performance of the dry-MARTINI/Upside hybrid interface for CGL particles by reducing unnecessary CGL pairlist work while preserving the validated physical table potentials and hybrid interface interactions.

# Architecture & Key Decisions

- Do not disable, scale, or bypass SC-env, BB-env, SC-CGL, CGL-target, or CGL-CGL interactions.
- CGL tables remain direct dry-MARTINI-derived spline tables. The table cutoff and taper are the physical force-field boundary.
- The performance change should target runtime candidate-pair construction, not the table energy surface.
- Because a CGL particle is a vector representing an extended DOPC body, the broad spherical support can be replaced only by a conservative directional support test: keep every pair that could place any represented DOPC bead pair within the dry-MARTINI cutoff.
- CGL-specific code under `src/martini_cg_lipid.cpp` is in the clean-slate set, so keep the implementation direct and remove any unused scaffolding introduced during the change.
- Apply the first directional filter to CGL-CGL and CGL-target pairs. Keep SC-CGL on its current spherical support until a side-chain body extent is stored explicitly enough to make a directional filter conservative for every rotamer row.

# Execution Phases

- [x] Phase 1: Inspect current CGL cutoff and pairlist paths.
- [x] Phase 2: Identify how dry-MARTINI cutoff and CGL geometry metadata are written into H5.
- [x] Phase 3: Implement conservative directional CGL pairlist filtering with minimal runtime scope.
- [x] Phase 4: Verify with build/tests and, if practical, compare pair counts or runtime behavior.
- [x] Phase 5: Document results and any residual risks.

# Known Errors / Blockers

- Stale hybrid `.up` files created before this change do not contain
  `compose_vector6d/max_axis_radius_ang`; reinjection correctly rejects them.
  Verification requiring CGL-target nodes must use a freshly generated stage.

# Review

- Implemented directional CGL pairlist filtering for CGL-CGL and CGL-target by
  testing the CGL capsule envelope against the dry-MARTINI 1.2 nm bead cutoff
  plus the pairlist buffer. The spline table cutoff/taper and force-field
  surface are unchanged.
- Added `max_axis_radius_ang` to DOPC-derived CGL metadata and propagated
  `bead_nonbonded_cutoff_nm`, `length_conversion_ang_per_nm`,
  `max_axis_radius_ang`, and `max_perp_radius_ang` into runtime CGL-CGL and
  CGL-target nodes.
- Pairlist cache invalidation now accounts for CGL vector rotation by treating
  endpoint rotation as body motion scaled by the CGL axial radius.
- Verification passed: Python compile, C++ build, regenerated MARTINI H5 files,
  rebuilt `dopc.h5` metadata audit, focused CGL bilayer smoke, fresh short 1RKL
  hybrid smoke through stage 7.0, and candidate-count audit on the fresh 1RKL
  stage.
