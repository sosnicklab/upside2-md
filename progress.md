# Progress Log

## 2026-05-20 (1RKL Replica-2 Shape and Box Drift Debug)
- Actions taken:
  - Reopened the 1RKL stability debug for the exact output path after the user reported residual trajectory-0 secondary-structure drift, trajectory-2 protein deformation, XY bilayer ellipse formation, and monotonic box `total_potential` lowering.
  - Replaced the completed old plan with an active plan focused on per-replica metrics, box/barostat behavior, and direct copied-HDF5 validation.
  - Audited stage-7 HDF5/log artifacts and found stage 7 uses a fixed square box, while CGL-CGL contact count and `cg_lipid_pair` energy drift strongly over `7.0 -> 7.2`.
  - Quantified geometry: stage-7 XY aspect ratio grows from about `1.09` to `1.31`, CGL nearest-neighbor p05 decreases from about `6.50 Å` to `5.68 Å`, and protein hbond retention collapses during `7.1`.
  - Selected copied-HDF5 validation: rigid-protein membrane relaxation under the same Hamiltonian, then production release from the relaxed membrane.
  - Found rigid-protein relaxation and bead-cutoff-only CGL-CGL patches still fail: both continue strong `cg_lipid_pair` energy descent.
  - Implemented CG-lipid table changes so explicit DOPC projections use the generic dry-MARTINI bead cutoff, subtract attractive CGL-CGL radial background, and remove remaining attractive CGL-CGL controls.
  - Updated stale-table validation and `cg_lipid_potentials.tex`.
  - Validated the accepted CGL-CGL control patch by directly modifying copied HDF5 files from the reported output.
- Files modified:
  - `py/martini_build_tables.py`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
  - `plan.md`
  - `progress.md`
- Verification:
  - Copied-HDF5 accepted validation: `example/16.MARTINI/outputs/martini_test_1rkl_hybrid_repcheck/checkpoints/1rkl.stage_7.0.up` plus continuation `1rkl.stage_7.1.up`, `10k` total steps.
  - Final metrics: hbond sum `26.75`, CA RMSD `4.16 Å`, Rg `12.47 Å`, XY aspect `1.08`, CGL nearest-neighbor min/p05/mean `6.278/6.729/7.821 Å`.
  - `python3 -m py_compile py/martini_build_tables.py py/martini_prepare_system_lib.py` passed.
  - Stale guard rejects the reported old `martini.h5`.
  - `git diff --check` passed.
