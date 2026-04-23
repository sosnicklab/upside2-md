# Workflow Check: 1RKL Outside-Bilayer Script Validity

## Project Goal
- Determine whether `run_sim_1rkl_outlipid.sh` is still valid after the recent Upside engine changes and the updates made to `run_sim_1rkl.sh`.
- Preserve the out-of-bilayer workflow intent:
  - start `1rkl` outside the bilayer,
  - keep the script consistent with current helper interfaces and runtime attrs,
  - avoid changing shared architecture unless a concrete incompatibility is found.

## Architecture & Key Decisions
- Treat this as a workflow compatibility audit first, not an automatic refactor.
- Compare `run_sim_1rkl_outlipid.sh` directly against the current `run_sim_1rkl.sh` and any shared helper/runtime interfaces it depends on.
- If the outlipid script is broken, fix only the concrete incompatibilities required to restore validity.
- Validation should prefer:
  - shell syntax checks,
  - config-generation / dry-run style checks,
  - targeted runtime or parser verification only where needed.

## Execution Phases
- [x] Phase 1: Compare `run_sim_1rkl_outlipid.sh` with `run_sim_1rkl.sh` and identify changed assumptions.
- [x] Phase 2: Validate the outlipid script against current helpers/runtime, and patch only confirmed incompatibilities if needed.
- [x] Phase 3: Record results in the workflow trackers and summarize whether the script is still valid.

## Known Errors / Blockers
- `run_sim_1rkl_outlipid.sh` is not currently valid for its stated purpose because the current prep path does not consume its outside-placement overrides.

## Review
- Comparison result:
  - `run_sim_1rkl_outlipid.sh` is now a thin wrapper that delegates execution to the updated `run_sim_1rkl.sh`,
  - this means it does inherit the recent shared workflow changes, including the new preproduction `fix_rigid.mode = rigid_body` path and the current continuation logic.
- Confirmed incompatibilities:
  - the wrapper exports `PROTEIN_PLACEMENT_MODE`, `PROTEIN_ORIENTATION_MODE`, and `PROTEIN_SURFACE_GAP`,
  - the current prep implementation in `py/martini_prepare_system.py` / `py/martini_prepare_system_lib.py` does not reference those settings at all,
  - the wrapper also exports `PROTEIN_LIPID_MIN_GAP` and `PROTEIN_LIPID_CUTOFF_MAX`, but the current base script hardcodes those values as readonly constants, so the wrapper cannot override them.
- Targeted validation:
  - `bash -n run_sim_1rkl.sh`
  - `bash -n run_sim_1rkl_outlipid.sh`
  - direct stage-0 prep check from repo root with the wrapper's current packing parameters:
    - `python3 py/martini_prepare_system.py --pdb-id 1rkl_outlipid_check ... --xy-scale 1.35 --box-padding-xy 20.0 --box-padding-z 50.0 --protein-lipid-cutoff 5.0`
  - measured on `example/16.MARTINI/outputs/outlipid_validity_check/prep/1rkl_outlipid_check.MARTINI.pdb`:
    - protein `z_min = 55.143 Å`
    - upper-leaflet `PO4 z_max = 91.635 Å`
    - clearance above upper leaflet = `-36.492 Å`
- Conclusion:
  - the generated structure is still embedded in the bilayer rather than above it,
  - so `run_sim_1rkl_outlipid.sh` is not currently valid as an "initially outside of bilayer" workflow.
