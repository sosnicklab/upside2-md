# Project Goal
Identify the modernization differences between the original Theano-based workflow and the PyTorch port that are causing gradient explosions, then propose precise corrections to restore stable gradients without artificial clipping.

## Architecture & Key Decisions
- Compare original vs modernized code paths file-by-file and function-by-function to trace gradient flow and data shapes.
- Preserve existing architecture unless a change is required for correctness; record any changes under Revised Decisions.
- Treat numerical precision (float32 vs float64) and parameter packing/unpacking as first-class suspects.

### Revised Decisions
- None yet.

## Execution Phases
- [x] Phase 1: Establish baseline by diffing the four target file pairs and logging key behavioral differences.
- [x] Phase 2: Trace gradient flow from `compute_divergence` through `backprop_deriv` to `get_d_obj_torch`, including tensor dtype and shape checks.
- [x] Phase 3: Identify mismatches (e.g., regularization scaling, metadata usage, kwargs changes) and map them to gradient instability.
- [x] Phase 4: Propose specific code fixes with file/line references and rationale.

## Known Errors / Blockers
- Original reference paths for `rotamer_parameter_estimation.py`, `run_upside.py`, and `upside_engine.py` live under `/Users/yinhan/Documents/upside2-md/py/` rather than the repo root.
