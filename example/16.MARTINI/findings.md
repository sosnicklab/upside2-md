# Findings

## 2026-03-05 (Pre-Production Rigid Control Regression Root Cause)
- Current generated stage files under `outputs/martini_test_1rkl_hybrid/checkpoints/` do not contain `/input/hybrid_control`; spot checks on `1rkl.stage_6.6.up` and `1rkl.stage_7.0.up` confirmed the group is absent.
- The existing C++ pre-production rigid-mask logic in `../../src/martini.cpp` only activates when `/input/hybrid_control` exists with:
  - `enable=1`
  - `preprod_protein_mode=rigid`
  - `prep_runtime_mode=dry_martini_prep`
- `run_sim_1rkl.sh` had been reduced to injecting only `hybrid_bb_map` and `hybrid_env_topology`, so the workflow no longer emitted the control metadata required to freeze protein before production.
- Restoring rigidity does not require a C++ change; it is sufficient for the workflow to synthesize `/input/hybrid_control` into each prepared stage file and set production activation explicitly (`activation_stage=production`, `active_runtime_mode=aa_backbone_explicit_lipid`).

## 2026-03-05 (Backbone-Only Workflow Schema Audit)
- Workflow-level legacy hook sweep over `example/16.MARTINI` (`*.sh`, `*.py`) confirms removal of:
  - stage-activation/sc-control hook helpers (`set_hybrid_activation_stage`, `set_hybrid_sc_controls`, `assert_hybrid_stage_active`),
  - legacy production handoff env var (`UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS`),
  - legacy SC/hard-sphere/RMSD runtime knob wiring in shell workflows.
- Active mapping workflow now uses only:
  - `/input/hybrid_bb_map`
  - `/input/hybrid_env_topology`
- Mapping writer (`prepare_system_lib.py`) no longer emits legacy control/sidechain mapping groups, and mapping validator has been rewritten to this backbone-only schema.

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

## 2026-03-05 (Membrane Channel Audit + Depth-Table Assumptions)
- `parameters/ff_2.1/membrane.h5` datasets present:
  - `cb_energy` shape `(20, 2, 18)` with residue names `ALA..VAL`,
  - `hb_energy` shape `(2, 2, 18)`,
  - `icb_energy`, `ihb_energy`, `burial_nodes`, `names`.
- No explicit sidechain/sheet membrane channels were found in this file:
  - missing candidates: `sc_energy`, `sidechain_energy`, `sheet_energy`.
- File attributes confirm z-grids used for interpolation:
  - `cb_z_min=-17.0`, `cb_z_max=15.0`
  - `hb_z_min=-33.0`, `hb_z_max=15.0`
- Bilayer depth extraction from `example/16.MARTINI/pdb/bilayer.MARTINI.pdb`:
  - lipid resnames present in template: `DOP` (plus ions `NA`, `CL`),
  - per-bead depth computed as `|z-z_center|` with `z_center=42.484580 A`.
- Produced interaction-table artifact (`outputs/depth_interaction_table.csv`) has 14 DOPC bead names with mapped MARTINI bead types from `dry_martini_v2.1_lipids.itp`.
- Metadata artifact confirms current-term availability for this plan:
  - `backbone_cb_energy=1`
  - `hbond_hb_energy=1`
  - `sidechain_term=0`
  - `sheet_term=0`
