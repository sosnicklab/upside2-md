# Progress Log

## 2026-05-25 1AFO Coarse-Lipid Helix Bend
- Actions taken:
  - Confirmed current 1AFO coarse outputs already include
    `cg_lipid_rotamer_sc`; the remaining bend was not stale missing CGL-SC
    coupling.
  - Localized the kink to stage-7 burn-in: stage-6 output and stage-7 prepared
    input were not bent, and handoff minimization alone did not bend either
    helix.
  - Tested charged-terminal CGL controls on copied HDF5 artifacts and found the
    coarse chain-1 bend disappears when `Qd`/`Qa` BB proxy CGL controls are
    excluded-volume-only.
  - Implemented a physical split for charged BB proxy targets:
    `cg_lipid_target_charged_bb_excluded_volume` keeps nonnegative CGL controls
    while ordinary BB/protein targets keep the full CGL-target table.
  - Updated `example/16.MARTINI/cg_lipid_potentials.tex` to document why the
    explicit-DOPC charged-terminal PMF is not transferable to the additive CGL
    projection.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
- Test results:
  - `python3 -m py_compile py/martini_prepare_system_lib.py py/martini_prepare_system.py`
    passed.
  - Regenerated CGL nodes on `/private/tmp/1afo_charged_bb_split_inject.up`:
    68 ordinary CGL targets, 4 charged-BB excluded-volume targets, and 98 ion
    excluded-volume targets; charged-BB and ion control minima are `0.0`.
  - One-step `obj/upside` smoke test on the regenerated copied HDF5 passed.
  - Full generated-path 1AFO stage-7 handoff plus 40k-step burn-in on
    `/private/tmp/1afo_generated_split_stage7_burnin.up` passed; final helix
    half-chain angles were `32.34 deg` and `36.88 deg`, compared with the old
    coarse output's `37.33 deg` and `85.72 deg`.
- Failures and fixes:
  - Initial copied injection check used the wrong example-local DOPC HDF5 path;
    reran with root `parameters/dryMARTINI/dopc.h5`.
