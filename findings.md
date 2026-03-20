# Findings

## External / Technical Findings
- Saved baseline production artifact:
  - file: `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`
  - stage-7 output contains actual protein AA carrier atoms with roles `N/CA/C/O` on runtime indices `4090..4322`.
- Actual carrier geometry from that artifact:
  - restrained `N-CA` and `CA-C` bonds remain near expected values over saved frames.
  - `C-O` is not part of the injected `Distance3D` / `Spring_bond` backbone restraint set and drifts badly in the saved output.
  - measured `C-O` carrier distance changes from about `1.234 Å` at frame 0 to a mean of about `17.90 Å` and a max of about `30.30 Å` at the last saved frame.
- Source-path finding:
  - `example/16.MARTINI/run_sim_1rkl.sh` injects four active BB carriers with weights `[N, CA, C, O]`.
  - `augment_production_rotamer_nodes()` builds production backbone nodes using only `N/CA/C` (`ref_n_atom = 3 * len(residue_ids)`), so `O` receives no direct backbone restraints.
  - `src/martini.cpp::refresh_bb_positions_if_active(...)` still uses the four-carrier weighted COM, so drifting `O` contaminates the refreshed MARTINI `BB` proxy.
- Export-path finding:
  - `example/16.MARTINI/extract_martini_vtf.py::reconstruct_backbone_aa(...)` does not use the actual runtime AA carriers when they exist.
  - It backmaps each residue by adding a pure translation from the live `BB` coordinate to the stored reference `N/CA/C/O` geometry, with no residue rotation.
  - In the saved stage-7 artifact, that reconstructed backbone diverges strongly from the actual runtime carriers by the last saved frame:
    - mean carrier-position mismatch about `7.12 Å`,
    - max carrier-position mismatch about `20.75 Å`,
    - mean peptide `C(i)-N(i+1)` bond in the reconstruction about `7.78 Å` versus about `1.37 Å` in the actual runtime carriers.
- Implemented fix:
  - `src/martini.cpp` now reads `hybrid_bb_map/reference_atom_coords` into runtime state and reconstructs the active `O` carrier from the current `N/CA/C` local frame before BB COM refresh.
  - `example/16.MARTINI/extract_martini_vtf.py` now exports actual runtime carrier coordinates when `hybrid_bb_map/atom_indices` resolve to runtime roles `N/CA/C/O`.
- Focused verification:
  - rebuilt successfully with `cmake --build obj`;
  - strict-copy replay from the broken last saved frame confirmed runtime repair:
    - pre-run `/input/pos` `C-O` mean/max `17.90/30.30 Å`,
    - first saved `/output/pos` frame after replay `C-O` mean/max `1.32/1.57 Å`;
  - exporter verification on that replay file:
    - `use_runtime_carriers=True`,
    - exported AA backbone coordinates matched actual runtime carriers exactly on the checked frame (`0.0 Å` mean/max mismatch).

## Lessons
- None yet.
