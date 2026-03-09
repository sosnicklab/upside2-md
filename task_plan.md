Project Goal
- Make `ConDiv_symlay` use one fixed DOPC membrane thickness for all training targets instead of reading per-target `.thickness` files.

Architecture & Key Decisions
- Preserve the current training and worker interfaces as much as possible; patch only the thickness source.
- Add one workflow-level config value for membrane thickness, exposed through an environment variable so the fixed DOPC value is explicit and easy to override.
- Default the fixed thickness to `30.2 A`, which is the median of the current `ConDiv_symlay/upside_input2/*.thickness` values. This removes per-target variation while staying close to existing practice.
- Remove the hard dependency on per-target `.thickness` files during training-set construction.
- Keep restart compatibility with old cached checkpoints by ignoring serialized target-specific thickness values at worker-dispatch time and using the current global workflow thickness instead.

Execution Phases
- [x] Audit all remaining membrane-thickness inputs and choose the minimal global-config path.
- [x] Patch `ConDiv_symlay/ConDiv_mem.py` to use one global membrane thickness everywhere.
- [x] Update README / tracking files to describe the fixed-thickness DOPC model.
- [x] Run focused verification:
  - syntax check
  - local config-generation smoke test
  - confirm worker args no longer depend on target-specific thickness
- [x] Record review notes and residual risks.

Known Errors / Blockers
- Current code still reads `base + ".thickness"` for every target and passes that per-target value into `membrane_thickness`.
- Existing cached `cd_training.pkl` files may still contain varying `Target.thickness` values, so worker dispatch must ignore them to make the change effective on restart.

Review
- Added one global workflow thickness knob, `CONDIV_DOPC_THICKNESS`, in [ConDiv_mem.py]. The default is `30.2 A`, chosen as the median of the current `ConDiv_symlay/upside_input2/*.thickness` files so the workflow stops varying thickness by target while staying close to prior practice.
- `ConDiv_symlay/ConDiv_mem.py` no longer reads `base + ".thickness"` during training-set construction. New `Target` rows are populated with the global thickness instead.
- Worker dispatch now passes `state["membrane_thickness"]` instead of `target.thickness`, so old cached `cd_training.pkl` files with varied serialized thicknesses are neutralized on restart.
- The initialized run state now stores `membrane_thickness`, and `_apply_runtime_config()` preserves that stored value unless `CONDIV_DOPC_THICKNESS` is explicitly set in the current environment.
- Updated `ConDiv_symlay/README.md` to document the fixed-thickness DOPC model and the new override variable.
- Local verification:
  - `python3 -m py_compile ConDiv_symlay/ConDiv_mem.py` passed.
  - One-target initialization smoke test with `CONDIV_DOPC_THICKNESS=30.2` stored `30.2` in both `cd_training.pkl` and `condiv_init.pkl`.
  - A second one-target smoke test used target `1orq`, whose legacy file is `29.4`, and still stored `30.2`, confirming the workflow now ignores per-target `.thickness` files.
- Residual risk:
  - Existing code still carries the legacy `Target.thickness` field for checkpoint compatibility, but the active worker path no longer trusts it.
