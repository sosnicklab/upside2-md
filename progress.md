# Progress Log

## 2026-06-15 CGL Directional Pairlist Support

- Actions taken:
  - Added DOPC-derived `max_axis_radius_ang` metadata.
  - Propagated CGL body geometry and dry-MARTINI bead cutoff attrs into runtime CGL-CGL and CGL-target nodes.
  - Replaced broad spherical CGL-CGL and CGL-target candidate filters with conservative capsule-envelope filters using the dry-MARTINI 1.2 nm bead cutoff plus pairlist buffer.
  - Kept SC-CGL pairlists on their existing spherical support until side-chain body extents are explicit enough for a conservative directional filter.
  - Rebuilt installed MARTINI parameter files, including `parameters/dryMARTINI/dopc.h5`.
- Files modified:
  - `src/martini_cg_lipid.cpp`
  - `py/martini_cg_lipid_params.py`
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system_lib.py`
  - `parameters/dryMARTINI/dopc.h5`
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - `python3 -m py_compile py/martini_cg_lipid_params.py py/martini_build_tables.py py/martini_prepare_system_lib.py` passed.
  - `make -C obj -j4` passed.
  - `example/16.MARTINI/build_martini_h5_m1.sh` completed successfully.
  - Rebuilt `dopc.h5` records `max_axis_radius_ang=14.55848`, `max_perp_radius_ang=4.114344`, and `bead_nonbonded_cutoff_nm=1.2` on CGL table groups.
  - `example/16.MARTINI/test_cg_bilayer/run_test.sh` passed with `validation=PASS`.
  - Fresh short `1rkl` coarse hybrid smoke under `outputs/phase_cgl_directional_1rkl_smoke` completed through stage 7.0.
  - Candidate audit on the fresh 1RKL stage: CGL-CGL candidates reduced from `14678` spherical pairs to `4579` capsule pairs; CGL-target candidates reduced from `1238` to `679`.
- Residual notes:
  - Old hybrid `.up` files without `compose_vector6d/max_axis_radius_ang` are intentionally rejected during reinjection with rebuilt CGL tables. Regenerate stages instead of mixing old geometry metadata with new tables.
