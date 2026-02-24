# Task: Create `run_sim_1rkl_new.sh` with Hybrid Active at Production

## 1. Project Goal
Add a new workflow `example/16.MARTINI/run_sim_1rkl_new.sh`, based on `run_sim_1rkl_rigid_dry.sh`, that activates Upside hybrid at stage `7.0` and enables all designed hybrid production components.

## 2. Architecture & Key Decisions
- Scope:
  - Clone `run_sim_1rkl_rigid_dry.sh` into a new workflow script.
  - Keep the same stage topology/process (`6.0 -> 7.0`).
  - Enable hybrid activation at production stage.
  - Enable production hybrid controls and rotamer/backbone augmentation.
- Key decisions:
  - Use `HYBRID_ACTIVATION_STAGE=production` as default in the new workflow.
  - Keep preproduction rigid behavior via existing `preprod_protein_mode=rigid`.
  - Use production handoff carrier refresh for `6.6 -> 7.0`.
  - Use `PROD_TIME_STEP=0.002` and stage `7.0` VTF extraction mode `2`.
- Revised Decisions:
  - Add explicit production hybrid assertion before running stage `7.0`.

## 3. Execution Phases
- [x] Phase 1: Inspect rigid-dry workflow and identify hybrid-disable points.
- [x] Phase 2: Create new workflow script and switch activation defaults.
- [x] Phase 3: Wire production hybrid stack (SC controls + rotamer nodes + handoff mode).
- [x] Phase 4: Validate script syntax and run short smoke workflow.
- [x] Phase 5: Log results in `progress.md`.

## 4. Known Errors / Blockers
- None currently.
