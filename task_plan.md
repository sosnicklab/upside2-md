# Task Plan

## Project Goal
- Find the source of the delayed stage-7 hybrid instability in `example/16.MARTINI/run_sim_1rkl.sh` that remains stable through about 5000 steps but explodes by about 10000 steps, and implement a minimal fix in `src/`.

## Architecture & Key Decisions
- Start from the provided reproduction path in `example/16.MARTINI/` and trace into `src/` only as needed.
- Prefer the smallest fix that restores correct behavior without unrelated architectural changes.
- Verification must reach the reported delayed-failure horizon or an equivalent targeted replay from the saved stage-7 artifact.
- Production-stage diagnosis is based on the saved stage-7 artifacts already present under `example/16.MARTINI/outputs/`.
- The previously fixed unconstrained `O` carrier bug is now treated as a resolved prerequisite issue, not sufficient explanation for the reported 10000-step explosion.
- The active-`BB` virtual-site momentum defect is also treated as a real but non-dominant issue because the user reran the full workflow and reproduced the same 5000 -> 10000 blow-up.
- Focus next on hybrid algorithm and implementation details in `src/` that can inject energy gradually, especially coordinate refresh, state handoff, alignment semantics, and timescale coupling between the dry-MARTINI and Upside subsystems.

## Execution Phases
- [x] Review existing project findings relevant to MARTINI / hybrid workflow
- [x] Inspect the example workflow and identify the production-stage configuration path
- [x] Confirm the delayed instability from the saved stage-7 production log
- [x] Inspect saved stage-7 artifacts and existing project diagnostics to localize the instability class
- [x] Trace the code path in `src/` that can produce the delayed energy blow-up
- [x] Identify and retire the non-dominant delayed-instability candidates already disproved by the user's rerun
- [x] Implement the current minimal source fix with minimal surface area
- [ ] Verify behavior at the delayed instability horizon
- [x] Record review notes and residual risks

## Known Errors / Blockers
- Confirmed delayed-failure reproduction:
  - `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/logs/stage_7.0.log` is stable at step `5000` and blows up by step `10000`:
    - `5000 / 100000 ... potential   595.84, martini_potential -26362.14, total -25766.30`
    - `10000 / 100000 ... potential 14036.12, martini_potential   516.82, total 14552.94`
- Resolved prerequisite issue:
  - stage-7 injects active `N/CA/C/O` carriers, but the generated Upside backbone force field only restrained `N/CA/C`; the unconstrained `O` carrier drift was fixed by reconstructing `O` from the live `N/CA/C` frame before BB refresh.
- Revised root-cause focus:
  - the active-`BB` virtual-site momentum defect was real and is fixed, but it was not sufficient to remove the user's observed full-workflow explosion;
  - the current source-level bug is the stage-7 `integration_rmsd_align` path mutating real protein coordinates and momenta during production even though the original hybrid design only intended rigid alignment for coupling/reference coordinates;
  - that raw-state rewrite mixes an external rigid-body projection into the live integrator state every cycle, which is another non-Hamiltonian energy-injection path and matches the class of delayed hybrid instability under investigation.
- Remaining blocker:
  - I have not yet completed a full 10000-step production replay with the current alignment fix, so long-horizon validation remains open.
  - The exact production replay was run and verified through its `5000`-step checkpoint (`time 10.0`) with changed energies, then stopped because reaching the full `10000`-step horizon is hour-scale in this environment.

## Review
- Prior completed fix:
  - `src/martini.cpp` now restores each active backbone `O` carrier from the live `N/CA/C` frame using the stored `hybrid_bb_map/reference_atom_coords` before MARTINI `BB` refresh is computed.
  - `example/16.MARTINI/extract_martini_vtf.py` now prefers actual runtime AA carrier coordinates from `hybrid_bb_map/atom_indices` when those indices resolve to runtime roles `N/CA/C/O`; the old translation-only reconstruction remains fallback-only.
- Prior verification:
  - rebuilt with `cmake --build obj`;
  - strict-copy replay from a broken saved stage-7 state (`/tmp/1rkl.stage_7.0.o_fix_replay_strict.up`) started from `C-O` mean/max `17.90/30.30 Å` in `/input/pos`;
  - after one production step, the first saved output frame had `C-O` mean/max `1.32/1.57 Å`;
  - exporter check on that replay file reported `use_runtime_carriers=True` and `reconstruct_backbone_aa(...)` matched the actual runtime carriers exactly (`0.0 Å` mean/max mismatch on the checked frame).
- Current review status:
  - implemented a new `src/` fix for active hybrid `BB` proxy dynamics:
    - active `BB` proxy atoms are now installed into the dynamic fixed-atom mask whenever hybrid production is active, so they are treated as virtual coordinates rather than independent thermostatted/integrated particles;
    - initial momentum zeroing now runs after startup thermalization as well, so those virtual `BB` sites do not begin stage 7 with ghost momentum.
  - implemented the current `src/` fix for the reopened long-horizon instability:
    - `src/martini.cpp::align_active_protein_coordinates(...)` no longer rotates or translates the real integrated protein coordinates / momenta during stage-7 production;
    - the path is now reduced to BB refresh plus RMSD / reference bookkeeping, so rigid alignment remains coupling-side logic instead of mutating the live trajectory state.
  - verification completed:
    - rebuilt with `cmake --build obj`;
    - short recorded-momentum replay from the exact `6.6 -> 7.0` workflow handoff confirmed the pre-fix bug signature had existed: frame-0 active `BB` momentum norm was nonzero on all 31 proxies (mean `4.55e-2`, max `8.74e-2`);
    - the same replay after the fix confirmed all active `BB` proxy momenta are exactly zero on every checked saved frame (`5/5` frames, max `0.0`);
    - active `BB` positions still track the carrier-weighted COM after the fix (max mismatch about `2.15e-6 Å` across the checked frames);
    - early production energies remained well-behaved in the 5-step per-frame replay (`martini_potential -23613.40 -> -23801.97` from step `0 -> 4`);
    - a 5-step A/B replay from the exact stage-7 handoff showed that `integration_rmsd_align_enable=1` and `integration_rmsd_align_enable=0` now produce the same saved raw trajectory to numerical noise only (position / potential max absolute difference about `3.4e-5`);
    - the exact fixed-code production replay from the real workflow handoff was carried through the critical midpoint and no longer matches the broken baseline at step `5000`:
      - broken baseline `5000`: `potential 595.84`, `martini_potential -26362.14`, `total -25766.30`;
      - current replay `5000` / `time 10.0`: `potential 768.92`, `martini_potential -26414.40`, `total -25645.48`.
  - residual risk:
    - the full reported 10000-step horizon has not been rerun to completion yet, so this remains a source-level fix with strong short-horizon equivalence checks plus a verified healthy midpoint, but not full end-to-end long-run confirmation;
    - reconstructed backbone `O` carriers still show independent nonzero momentum in short saved replays, which is a separate virtual-site correctness issue that likely needs proper chain-rule force redistribution rather than another quick fixed-atom patch.
