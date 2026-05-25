# Progress Log

## 2026-05-25 1RKL Additive Carrier Force Routing
- Actions taken:
  - Audited current 1RKL and 1AFO full/coarse HDF5 outputs for terminal BB
    endpoint types, endpoint charges, generic MARTINI pair ownership, and CGL
    node ownership.
  - Rejected the earlier BB-proxy-only carrier exclusion after the user
    correction: N/CA/C/O carriers must still accumulate sidechain/rotamer and
    backbone-environment force contributions.
  - Found the implementation bug in the C++ projection path: mapped N/CA/C/O
    carriers have `ROLE_BB` plus a BB map index, so direct carrier gradients
    were being routed through BB-proxy projection instead of being added
    directly.
  - Updated the projection checks so only true virtual BB proxies are projected;
    direct carrier gradients remain direct.
  - Kept CGL targets proxy-based rather than adding independent N/CA/C/O carrier
    target copies, because copied 1RKL probes showed independent carrier
    targets collapse the secondary structure.
  - Fixed fragment metadata and validation so every residue fragment, including
    broken 1AFO strands, is checked for `Qd/+1` and `Qa/-1` BB endpoints in
    both full and coarse workflows.
  - Updated `example/16.MARTINI/cg_lipid_potentials.tex` to document additive
    carrier force accumulation and direct-vs-projected gradient routing.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `src/martini_hybrid.cpp`
  - `src/martini_potential.cpp`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
- Test results:
  - `cmake --build obj` passed.
  - `python3 -m py_compile py/martini_prepare_system_lib.py py/martini_prepare_system.py`
    passed.
  - Copied 1RKL proxy-target/additive-force probe passed stage validation,
    minimized from `657.02` to `60.93 E_up`, and completed 40k-equivalent MD
    with final hbond `31.43` and kinetic ratio `0.989`; the stale coarse output
    ends at hbond `23.71`.
  - Rejected independent carrier CGL targets: 40k probe final hbond `6.5`;
    after projection correction, 10k probe still fell to hbond `9.6`.
  - Current 1AFO CGL injection on a copied file produced `182 CGL x 170`
    target particles (`72` BB proxies plus `98` ions), retained both fragment
    endpoint pairs (`0-35` and `36-71`) as `Qd/+1` and `Qa/-1`, and missed no
    BB proxy target.
- Failures and fixes:
  - A temporary 1AFO injection using the small `test_cg_lipid` table skipped
    `cg_lipid_rotamer_sc` because that fixture lacks CG-SC table entries; the
    endpoint and target-count validation was rerun without treating that
    fixture limitation as a production missing-node failure.

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
