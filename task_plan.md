# Task Plan

## Project Goal
- Find the source of the bug causing the all-atom protein backbone to break during the production stage of `example/16.MARTINI/run_sim_1rkl.sh`, and implement a minimal fix.

## Architecture & Key Decisions
- Start from the provided reproduction path in `example/16.MARTINI/` and trace into `src/` only as needed.
- Prefer the smallest fix that restores correct behavior without unrelated architectural changes.
- Verification must include at least one focused reproduction or invariant check tied to backbone integrity during production.
- Production-stage diagnosis is based on the saved stage-7 artifacts already present under `example/16.MARTINI/outputs/`.
- The active hybrid workflow currently injects four AA backbone carriers (`N/CA/C/O`), but the injected Upside backbone force field only constrains `N/CA/C`. The `O` carrier must therefore be treated as a derived coordinate, not an unconstrained dynamic degree of freedom.
- The production VTF export should prefer the actual runtime AA carrier coordinates when they exist; translation-only reconstruction from per-residue `BB` offsets is not sufficiently faithful.

## Execution Phases
- [x] Review existing project findings relevant to MARTINI / hybrid workflow
- [x] Inspect the example workflow and identify the production-stage configuration path
- [x] Trace the code path in `src/` that affects backbone handling in this workflow
- [x] Identify the root cause and document the decision
- [x] Implement the fix with minimal surface area
- [x] Verify behavior with focused checks
- [x] Record review notes and residual risks

## Known Errors / Blockers
- Current root cause:
  - `augment_production_rotamer_nodes()` injects `N/CA/C/O` carriers, but the generated production backbone nodes only constrain `N/CA/C`.
  - Saved stage-7 output confirms unconstrained `O` drift: `C-O` carrier distance grows from about `1.23 Å` at frame 0 to a mean of about `17.9 Å` and a max of about `30.3 Å` by the last saved frame in `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`.
  - `extract_martini_vtf.py` mode 2 compounds the issue by ignoring actual runtime carriers and reconstructing per-residue AA backbone coordinates from reference geometry plus a pure translation from the current `BB` position.

## Review
- Implemented:
  - `src/martini.cpp` now restores each active backbone `O` carrier from the live `N/CA/C` frame using the stored `hybrid_bb_map/reference_atom_coords` before MARTINI `BB` refresh is computed.
  - `example/16.MARTINI/extract_martini_vtf.py` now prefers actual runtime AA carrier coordinates from `hybrid_bb_map/atom_indices` when those indices resolve to runtime roles `N/CA/C/O`; the old translation-only reconstruction remains fallback-only.
- Verification:
  - rebuilt with `cmake --build obj`;
  - strict-copy replay from a broken saved stage-7 state (`/tmp/1rkl.stage_7.0.o_fix_replay_strict.up`) started from `C-O` mean/max `17.90/30.30 Å` in `/input/pos`;
  - after one production step, the first saved output frame had `C-O` mean/max `1.32/1.57 Å`;
  - exporter check on that replay file reported `use_runtime_carriers=True` and `reconstruct_backbone_aa(...)` matched the actual runtime carriers exactly (`0.0 Å` mean/max mismatch on the checked frame).
- Residual risk:
  - the runtime fix reconstructs `O` coordinates from the local `N/CA/C` frame before force evaluation and logging, but does not introduce a new explicit `O`-specific force term or momentum constraint. This resolves the observed backbone breakage and BB-refresh contamination with minimal impact, but the `O` carrier still behaves as a derived coordinate rather than a fully independent dynamical atom.
