# AA-Backbone 1RKL Workflow

## Project Goal
- Replace the `1rkl` dry-MARTINI protein presentation in `example/16.MARTINI/run_sim_1rkl.sh` with a runtime system that contains:
  - all-atom protein backbone atoms (`N/CA/C/O`),
  - DOPC lipid environment,
  - ions and water/environment particles,
  - probabilistic sidechain-lipid interactions.
- Hold the AA backbone rigid through stages `6.0` to `6.6`.
- Release the AA backbone by default at stage `7.0`, while keeping an optional production rigid-backbone override.
- Remove the retired CG/hybrid workflow-16 code instead of leaving blocked or dormant branches behind.

## Architecture & Key Decisions
- Keep the shared HDF5 metadata group names `hybrid_bb_map` and `hybrid_env_topology` because active runtime/export helpers already consume them.
- Remove the workflow-16 dry-MARTINI protein generation path:
  - no martinize call,
  - no MARTINI protein CG PDB input,
  - no MARTINI protein ITP dependency for stage conversion.
- Represent protein atoms in stage files directly as runtime AA backbone atoms:
  - use a dedicated placeholder MARTINI type for zero nonbonded interaction,
  - assign physical per-atom masses from `N/CA/C/O`,
  - mark them as `PROTEINAA`.
- Activate the probabilistic sidechain machinery on all prepared stages `6.0` through `7.0`, not only production.
- Use explicit stage-file spatial holds before production:
  - absolute `fix_rigid` hold for the AA backbone,
  - `z_atom_indices` hold for environment `PO4` atoms.
- Retire `PROTEIN_ENV_INTERFACE_SCALE` for this workflow rather than remapping it.

## Execution Phases
- [x] Phase 1: Replace the local plan/findings entries for the AA-backbone-only task and record the user correction about deleting retired code.
- [x] Phase 2: Refactor prep metadata generation in `py/martini_prepare_system.py` / `py/martini_prepare_system_lib.py` to emit AA-backbone runtime metadata and no CG protein path.
- [x] Phase 3: Refactor `convert_stage(...)` so AA backbone atoms are accepted directly and receive zero nonbonded placeholder typing plus physical masses.
- [x] Phase 4: Simplify `example/16.MARTINI/run_sim_1rkl.sh` to the single AA-backbone workflow, inject sidechain nodes for stages `6.0` through `7.0`, and remove retired controls.
- [x] Phase 5: Update continuation and VTF export helpers for the AA-backbone schema.
- [x] Phase 6: Run syntax checks and reduced smoke tests, then record the review and verification results.
- [x] Phase 7: Audit the current prep output against the last committed bilayer coordinates and restore committed initial bilayer `z` placement if the AA-backbone prep changed it.
- [x] Phase 8: Re-verify prep artifacts and reduced workflow behavior after the bilayer-baseline fix.

## Known Errors / Blockers
- No blocker is known yet.
- The reduced smoke workflow still shows extremely large MARTINI energies in several preproduction stages, but the same reduced run completes through `7.0`; this task changed representation/plumbing, not the underlying packing calibration.

## Review
- Implemented the AA-backbone-only `1rkl` workflow and removed the retired workflow-16 dry-MARTINI protein prep/control surface from the active scripts.
- Verified:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - in-memory Python syntax checks for `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, and `py/martini_extract_vtf.py`
  - direct prep smoke to `outputs/aabb_smoke`
  - direct stage `6.0` conversion smoke
  - reduced end-to-end AA-backbone workflow run to `outputs/martini_test_1rkl_aabb_smoke`
  - reduced stage-7 continuation smoke to `outputs/martini_test_1rkl_aabb_continue_smoke`
- Artifact checks on the smoke outputs confirmed:
  - `stage_6.0` contains `fix_rigid` over 124 AA backbone atoms,
  - `stage_7.0` omits `fix_rigid` by default,
  - `particle_class` uses `PROTEINAA`,
  - runtime AA backbone atom roles are `N/CA/C/O`,
  - sidechain nodes are present on both `6.0` and `7.0`.
- Follow-up bugfix verification confirmed the preproduction spatial holds:
  - rebuilt `obj/upside` after adding user-defined `z_atom_indices` support in `/input/fix_rigid`,
  - reduced workflow run `outputs/martini_test_1rkl_aabb_hold_smoke` completed through `7.0`,
  - `stage_6.2.up` contains `124` absolute protein holds plus `282` `PO4` Z holds,
  - over a 5-step `6.2` MD segment:
    - protein max displacement = `0.0 Å`,
    - `PO4` max `|Δz|` = `0.0 Å`,
  - `stage_7.0.up` still omits `fix_rigid` / `z_atom_indices` by default.
- Follow-up prep-baseline audit confirmed the initial bilayer is back on the committed `z` frame:
  - the regression source was the prep path using the backbone-only atom cloud for XY sizing, lipid overlap removal, and `min_box_z_target`,
  - the fix uses a prep-only envelope of AA backbone atoms plus one sidechain centroid proxy per residue while keeping the runtime system backbone-only,
  - direct prep comparison against `outputs/martini_test_1rkl_hybrid/hybrid_prep/hybrid_packed.MARTINI.pdb` now shows:
    - identical lipid count: `280`,
    - identical XY box: `111.159 x 111.159 Å`,
    - restored box `z`: `109.665 Å` versus committed `110.196 Å`,
    - restored `PO4` mean `z`: `53.656 Å` versus committed `54.052 Å`,
    - unchanged `PO4` thickness within `0.006 Å`,
  - reduced workflow run `outputs/martini_test_1rkl_aabb_bilayerfix_smoke` completed through `7.0`,
  - its `stage_6.2.up` still contains `124` absolute protein holds and `280` `PO4` Z holds with zero observed displacement in the 5-step segment.
