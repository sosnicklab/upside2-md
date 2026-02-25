# Task: Enforce RMSD Alignment for Per-Step AA Coordinate Updates

## 1. Project Goal
Ensure that all-atom coordinates updated during hybrid production integration are RMSD-aligned to the previous-step all-atom coordinates before they are applied, to remove rigid-body rotation drift.

## 2. Architecture & Key Decisions
- Scope:
  - Inspect hybrid AA update path in C++ runtime (`src/martini.cpp` and related integration code).
  - Verify whether per-step RMSD alignment already exists and whether it is applied at the correct point.
  - Patch logic so updated AA coordinates are aligned against previous-step AA coordinates before commit.
  - Keep behavior gated by existing hybrid control flags when appropriate.
- Key decisions:
  - Use previous-step AA backbone anchors as the alignment reference.
  - Apply alignment in the integration update path (not only during output/extraction).
  - Preserve existing workflow interfaces and configuration keys.
- Revised Decisions:
  - Invoke RMSD alignment before each force-evaluation/update substep in all integration-cycle variants.

## 3. Execution Phases
- [x] Phase 1: Locate current AA update + alignment code path in runtime.
- [x] Phase 2: Implement/adjust per-step RMSD alignment-before-apply behavior.
- [x] Phase 3: Build and run focused validation to confirm behavior.
- [x] Phase 4: Record outcomes in `progress.md`.

## 4. Known Errors / Blockers
- None currently.
