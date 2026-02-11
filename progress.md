# Progress Log

2026-02-11
- Initialized required workflow files: `task_plan.md`, `findings.md`, `progress.md`.
- Diffed key file pairs and located original file paths under `/Users/yinhan/Documents/upside2-md/py/`.
- Inspected `ConDiv.py` vs `ConDiv_original.py`, `rotamer_parameter_estimation.py`, `run_upside.py`, and `upside_engine.py`.
- Traced gradient flow from `compute_divergence` to `backprop_deriv` to `get_d_obj_torch` and noted mismatches.
- Rewrote `get_d_obj_torch()` to compute gradients w.r.t. packed `lparam` and restore the original regularization.
- Updated `compute_divergence()` to align sheet finite-difference handling with original semantics and per-restype maps.
- Removed gradient clipping and global norm rescaling from `backprop_deriv()`/`run_minibatch()`.
