# Workflow-16 Preproduction Rigid-Body Hold

## Project Goal
- Replace the AA-backbone workflow's absolute preproduction protein pinning with a rigid-body hold that preserves internal geometry while still allowing whole-protein translation and rotation.
- Keep the existing stage ladder and public shell interface unchanged.
- Preserve the current hard `PO4` `z` hold through preproduction and the current unfixed default production stage.

## Architecture & Key Decisions
- The observed production-stage "bilayer z shift" in `outputs/martini_test_1rkl_aabb` is not a true bulk bilayer drift bug:
  - bilayer/environment/total `z` drift stays small,
  - the large motion is protein-relative drift after stage `7.0` releases the preproduction hold.
- The requested workflow semantics are now:
  - protein rigid during minimization and preproduction,
  - protein not locked to absolute space,
  - rigid-body translation and rotation both allowed,
  - environment `PO4` atoms still fixed in `z` during preproduction.
- Do not add a new stage or new shell controls.
- The fix belongs in the shared `fix_rigid` runtime path:
  - extend `fix_rigid` with a new internal `rigid_body` mode,
  - keep existing missing-mode behavior as the current absolute pin to avoid regressions.
- Workflow `16` should opt into the new mode only for its preproduction protein hold.

## Execution Phases
- [x] Phase 1: Replace the local trackers for the rigid-body preproduction task and record the new user correction.
- [x] Phase 2: Extend the `fix_rigid` engine path to support rigid-body minimization and MD while preserving the legacy absolute-pin mode.
- [x] Phase 3: Update `run_sim_1rkl.sh` to write `fix_rigid.mode = rigid_body` for preproduction protein holds while leaving `PO4` `z` holds unchanged.
- [x] Phase 4: Rebuild the binary, run targeted syntax/workflow verification, and measure whether the AA-backbone production drift is reduced without reintroducing bilayer drift.

## Known Errors / Blockers
- No blocker yet.

## Review
- Implemented the new `fix_rigid` mode in the engine:
  - legacy missing-mode behavior still maps to the old absolute pin,
  - `mode = rigid_body` now keeps the selected atoms on a best-fit rigid transform,
  - MD projection rebuilds rigid-body-consistent momenta from translational and rotational motion,
  - minimization projection keeps the selected atoms on the rigid-body manifold while still allowing net translation/rotation descent.
- Updated the AA-backbone workflow preproduction writer:
  - `example/16.MARTINI/run_sim_1rkl.sh` now writes `fix_rigid.mode = rigid_body`,
  - the `PO4` `z_atom_indices` hold remains unchanged,
  - production stage `7.0` still removes `fix_rigid` entirely.
- Verification summary:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `PYTHONPYCACHEPREFIX=/tmp/upside_pycache python -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py example/16.MARTINI/run.py`
  - `cmake --build obj -j4`
  - reduced workflow smoke:
    - `RUN_DIR=example/16.MARTINI/outputs/martini_test_1rkl_rigidbody_smoke`
    - `MIN_60_MAX_ITER=5 MIN_61_MAX_ITER=5`
    - `EQ_62_NSTEPS=10 EQ_63_NSTEPS=10 EQ_64_NSTEPS=10 EQ_65_NSTEPS=10 EQ_66_NSTEPS=10`
    - `PROD_70_NSTEPS=200 PROD_FRAME_STEPS=10`
    - completed stages `6.0` through `7.0`
  - direct artifact checks on the reduced smoke outputs confirmed:
    - `stage_6.2.up` and `stage_6.6.up` contain `fix_rigid.mode = rigid_body`,
    - `stage_6.2.up` and `stage_6.6.up` preserve protein internal distances to about `1e-6 Å`,
    - `stage_6.2.up` and `stage_6.6.up` keep `PO4` `z` holds exact (`0.0 Å` drift),
    - preproduction protein COM motion is no longer exactly zero,
    - `stage_7.0.up` contains no `fix_rigid`.
  - direct production continuation from the prepared stage-7 handoff:
    - `obj/upside .../1rkl.stage_7.0.long.up --duration-steps 1000 --frame-interval 0.1 ...`
    - over the resulting 20 output frames:
      - bilayer `z` drift max = `0.066 Å`,
      - protein `z` drift max = `0.052 Å`,
      - bilayer-minus-protein relative `z` change max = `0.118 Å`.
