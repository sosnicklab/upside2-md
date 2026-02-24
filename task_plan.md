# Task: Fix AA Backbone Corner Artifact in `run_sim_1rkl_rigid_dry.sh` VTF Output

## 1. Project Goal
Verify all-atom (AA) protein-backbone mapping remains properly aligned to dry-MARTINI backbone in production stage, and fix VTF extraction so mapped AA atoms are not rendered at simulation-box corners.

## 2. Architecture & Key Decisions
- Scope:
  - Inspect production-stage hybrid mapping datasets in `.up` files and validate mapped AA backbone coordinates versus MARTINI BB anchors.
  - Inspect `extract_martini_vtf.py` mode-2 coordinate handling for hybrid outputs.
  - Patch extraction/mapping handling as needed, without changing simulation physics.
- Key decisions:
  - Prefer fixing post-processing (`extract_martini_vtf.py`) if mapping in `.up` is already correct.
  - If mapping data itself is inconsistent, patch upstream mapping/injection logic in workflow/prep path.
- Revised Decisions:
  - Keep `mode 1` output semantics as requested: include all MARTINI particles and AA backbone.
  - Exclude raw `PROTEINAA` carrier particles from visualization and reconstruct AA backbone from `hybrid_bb_map` in extraction.

## 3. Execution Phases
- [x] Phase 1: Reproduce and inspect production `.up`/`.vtf` artifact for `run_sim_1rkl_rigid_dry.sh`.
- [x] Phase 2: Validate AA-backbone mapping integrity in hybrid datasets.
- [x] Phase 3: Patch VTF extraction (or mapping path) to correctly handle mapped AA coordinates.
- [x] Phase 4: Re-run focused extraction/workflow checks to confirm corner artifact is gone.
- [x] Phase 5: Log outcomes in `progress.md`.

## 4. Known Errors / Blockers
- None currently.
