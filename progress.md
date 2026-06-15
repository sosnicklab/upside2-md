# Progress Log

## 2026-06-15 Cross-System CGL Cutoff and Speed Check

- Actions taken:
  - Compared active 1RKL and 1AFO coarse HDF5 cutoff metadata for CGL-CGL, CGL-target, and SC-CGL runtime nodes.
  - Checked the 1RKL and 1AFO wrapper scripts for workflow-specific CGL cutoff overrides.
  - Compared active stage-7 coarse/full timings for both systems.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - 1RKL and 1AFO both use CGL-CGL `cutoff_ang=41.3`, CGL-target `cutoff_ang=26.6`, SC-CGL `cutoff_ang=26.6`, bead nonbonded cutoff `1.2 nm`, DOPC axis radius `14.55848 A`, and DOPC perpendicular radius `4.114344 A`.
  - 1RKL timing: CGL `1707.99 us/systems/step`; full lipid `3328.29 us/systems/step`.
  - 1AFO timing: CGL `1923.47 us/systems/step`; full lipid `2522.99 us/systems/step`.
  - Final CGL orientation checks have no bad-parallel or flipped rods for either active system.
- Residual notes:
  - No code change was required because the active outputs already use the same cutoff metadata and satisfy the CGL-faster-than-full requirement for both systems.

## 2026-06-15 Active 1RKL CGL Orientation and Timing Check

- Actions taken:
  - Audited the reported active `1rkl.stage_7.0.vtf` issue from the HDF5 CGL/CGLD vectors rather than the display file alone.
  - Confirmed the old active stage-7 trajectory had real CGL orientation outliers, already present at the stage-7 production handoff.
  - Regenerated a default 1RKL coarse hybrid trajectory with the current directional CGL pairlist and strict stage handoff.
  - Replaced `example/16.MARTINI/outputs/martini_1rkl_hybrid` with the validated regenerated output.
  - Preserved the previous active output as `example/16.MARTINI/outputs/martini_1rkl_hybrid.pre_directional_pairlist_backup_20260615_140424`.
  - Compared current CGL and full-resolution timings for 1RKL hybrid and same-size DOPC-only bilayer tests.
- Files and outputs modified:
  - `example/16.MARTINI/outputs/martini_1rkl_hybrid`
  - `example/16.MARTINI/outputs/martini_1rkl_hybrid.pre_directional_pairlist_backup_20260615_140424`
  - `plan.md`
  - `findings.md`
  - `progress.md`
- Verification:
  - Active replacement `checkpoints/1rkl.stage_7.0.up` final CGL metrics: aligned-z min/p05/mean `0.770/0.849/0.949`, `bad_parallel=0`, `bad_flip=0`, leaflet crossings `0/0`.
  - Active replacement `logs/stage_7.0.log`: coarse CGL hybrid `1707.99 us/systems/step`.
  - Existing active full-resolution comparator `logs/stage_7.0.log`: `3328.29 us/systems/step`.
  - Same-size DOPC-only bilayer timings measured from command output: CGL `350.29 us/systems/step`; full-resolution lipid `797.32 us/systems/step`.
- Residual notes:
  - No additional performance code change was made in this phase because the current measured CGL paths are already faster than the comparable full-resolution paths.

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
