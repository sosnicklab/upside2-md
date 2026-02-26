# Findings

## 2026-02-17 (Runtime Semantics Audit)
- `output/potential` is logged in `../../src/main.cpp` via `sys->engine.potential`, and `../../src/deriv_engine.cpp` computes this as the sum over all active `potential_term` nodes. It is not a protein-only energy stream.
- Pre-production rigid hold in hybrid mode is configured via `/input/hybrid_control` (`preprod_protein_mode=rigid`) and enforced dynamically in C++ (`../../src/martini.cpp`) using `protein_membership` and `atom_roles`; `/input/fix_rigid` does not need to exist in stage files for this path.
- In pre-production stage files (6.0-6.6), observed potential nodes are MARTINI-side nodes (`martini_potential`, `dist_spring`, `angle_spring`, `dihedral_spring`, plus `restraint_position` where present), consistent with total-system potential variation while protein coordinates remain fixed.

## 2026-02-17 (CHARMM-GUI MDP Pressure Coupling Inputs)
- Source files inspected: `/Users/yinhan/Downloads/charmm-gui-7090685331/gromacs/step6.0_minimization.mdp` through `step6.6_equilibration.mdp` and `step7_production.mdp`.
- Pressure coupling settings observed:
  - Stages `6.0-6.1`: `Pcoupl=berendsen`, `Pcoupltype=semiisotropic`, `tau-p=4.0`, `compressibility=3e-4 0.0`, `ref_p=0.0 0.0`.
  - Stages `6.2-6.6`: `Pcoupl=berendsen`, `Pcoupltype=semiisotropic`, `tau-p=4.0`, `compressibility=3e-4 0.0`; `ref_p` not explicit (GROMACS default `1.0 bar` behavior used for equilibration intent).
  - Stage `7`: `pcoupl=no`.
- Unit conversion used (from AGENTS/CLAUDE project guidance): `1 bar = 0.000020659477 E_up/Angstrom^3`.

## 2026-02-24 (Hybrid Production Blow-Up Isolation)
- Controlled probes on `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.runtime_debug_hybrid_nosc.up` show that runaway expansion occurs whenever `hybrid_active=1`, even when:
  - `n_sc=0` (SC probabilistic path effectively off),
  - `production_nonprotein_hard_sphere=0`,
  - `integration_rmsd_align_enable=0`,
  - `exclude_intra_protein_martini=0`.
- Equivalent probes with hybrid inactive (`activation_stage=__hybrid_disabled__` or `enable=0`) remain stable over the same settings.
- Runtime diagnostics from the long failing run (`1rkl.stage_7.0.up`) indicate `diagnostics/sc_env_energy_total` remains O(10^2) while total MARTINI potential explodes to O(10^17), so SC env energy is not the dominant diverging term.
- BB coupling-specific probe: zeroing `/input/hybrid_bb_map/atom_mask` while keeping hybrid active removes the runaway behavior, isolating instability to active BB coupling logic.
- The unstable code path is centered in `../../src/martini.cpp`:
  - `refresh_bb_positions_if_active(...)` (active-stage BB coordinate overwrite from mapped carriers),
  - `project_bb_gradient_if_active(...)` (active-stage BB gradient redistribution/zeroing).

## 2026-02-24 (BB/AA Carrier Mass Consistency)
- BB bead mass was verified as `72/12 = 6.0` in runtime inputs, while injected AA backbone carriers were previously set to unscaled masses (`14/12`, `12/12`, `12/12`, `16/12`) summing to `54/12 = 4.5`.
- To enforce carrier mass sum equal to BB mass, required scaling is `72/54`.
- Effective carrier masses after applying `72/54` are:
  - `N: 72/54 * 14/12 = 1.5556`
  - `CA: 72/54 * 12/12 = 1.3333`
  - `C: 72/54 * 12/12 = 1.3333`
  - `O: 72/54 * 16/12 = 1.7778`
- Implemented in `example/16.MARTINI/run_sim_1rkl.sh` by scaling injected `ref_mass` values in `inject_hybrid_mapping()`.

## 2026-02-26 (Pre-production Pressure Coupling Debug)
- Stage-6.3 runtime confirms physical pressure target setup is correct in Upside units:
  - `target_p_xy = target_p_z = 2.0659477e-05 E_up/Angstrom^3` (`1 bar`),
  - `compressibility_xy = 14.521180763676 Angstrom^3/E_up` (`3e-4 bar^-1`),
  - `compressibility_z = 0.0` (semi-isotropic, no Z scaling).
- Reproduced observed expansion regime from the same stage-6.3 checkpoint:
  - short replay (`60` steps): `box_xy 111.08 -> 111.21` with positive instantaneous `Pxy ~ 1.47e-2 .. 3.06e-2 E_up/Angstrom^3`.
  - longer replay (`2000` steps): `box_xy 111.08 -> 116.54` at step `1000`, matching user-visible transient over-expansion.
- Root runtime defect found in NPT implementation:
  - `simulation_box::npt::update_node_boxes(...)` was a no-op, so node-local box caches used by minimum-image bonded terms (`dist_spring`, `angle_spring`, `dihedral_spring`) were stale during intra-stage barostat scaling.
- Implemented fix:
  - Added barostat callback registration in `src/box.h`/`src/box.cpp`.
  - Registered MARTINI-side callback in `src/martini.cpp` to scale node-local box lengths every NPT update.
- Important follow-up from validation:
  - An attempted `delta_t * interval` barostat-time aggregation was tested and reverted after it caused severe stage-6.3 over-expansion (`111.08 -> 136.79` at step `1000` in replay).
  - Current patch keeps existing coupling strength and fixes only node-box propagation correctness.

## 2026-02-26 (User-Directed Stage-6.3-Only Verification)
- Ran stage-6.3 only from existing stage-6.2 chain (no stage-6.0/6.1/6.4+ rerun):
  - reset `checkpoints/1rkl.stage_6.3.up` from `checkpoints/1rkl.stage_6.3.prepared.up`,
  - strict handoff from `checkpoints/1rkl.stage_6.2.up`,
  - reran `5000` steps and regenerated `outputs/martini_test_1rkl_hybrid/1rkl.stage_6.3.vtf`.
- Observed stage-6.3 box evolution after periodic-nonbonded fix:
  - initial: `111.086 x 111.086 x 110.196`,
  - at `1000` steps: `111.95 x 111.95 x 110.20`,
  - final frame: `112.001 x 112.001 x 110.196`,
  - max over saved frames: `112.020 x 112.020 x 110.196`.
- Comparison against pre-fix stage-6.3 behavior from same workflow:
  - previous `1000`-step box sample was `~116.50 x 116.50 x 110.20`,
  - current rerun reduces that early expansion by about `4.55 Å` in XY.
