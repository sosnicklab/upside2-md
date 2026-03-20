# Progress Log

## 2026-03-20
- Initialized task tracking files for the MARTINI hybrid workflow backbone-break debugging task.
- Inspected `example/16.MARTINI/run_sim_1rkl.sh` and the stage-7 hybrid runtime code in `src/martini.cpp`.
- Used the saved production artifact `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` to validate the failure mode directly from HDF5 output instead of relying only on logs or VTF appearance.
- Identified the simulation-side defect:
  - stage-7 injects active `N/CA/C/O` carriers, but the generated Upside backbone force field only restrains `N/CA/C`;
  - the `O` carriers are therefore unconstrained and drift far from `C` during production while still contributing to `hybrid_bb_map` BB refresh weights.
- Identified the export-side defect:
  - `extract_martini_vtf.py` reconstructs AA backbone coordinates from reference geometry plus BB translation instead of using the actual runtime carrier coordinates, which can make the exported stage-7 backbone look much worse than the real simulation state.
- Patched `src/martini.cpp`:
  - added runtime loading of `hybrid_bb_map/reference_atom_coords`;
  - restore each active backbone `O` carrier from the current `N/CA/C` frame before `refresh_bb_positions_if_active(...)` computes the MARTINI `BB` proxy position.
- Patched `example/16.MARTINI/extract_martini_vtf.py`:
  - when stage files already contain runtime AA carriers under `hybrid_bb_map/atom_indices` with roles `N/CA/C/O`, export those coordinates directly instead of using the older translation-only reconstruction.
- Rebuilt successfully with `cmake --build obj`.
- Verified runtime repair on a focused strict-copy replay:
  - seeded `/tmp/1rkl.stage_7.0.o_fix_replay_strict.up` from the broken last saved production frame with hybrid carrier refresh disabled at handoff;
  - confirmed pre-run `/input/pos` `C-O` mean/max `17.90/30.30 Å`;
  - ran one production step with the rebuilt `obj/upside`;
  - confirmed first saved `/output/pos` frame `C-O` mean/max `1.32/1.57 Å`.
- Verified exporter repair:
  - `build_backbone_projection_map(...)` reports `use_runtime_carriers=True` on the replay file;
  - exported AA backbone coordinates match actual runtime carriers exactly on the checked frame (`0.0 Å` mean/max mismatch).
