# Findings

## 2026-02-28 (SC Probability Weighting Re-Check)
- The capped-to-uncapped force ramp added in `../../src/martini.cpp` does not replace or bypass the existing SC probability weighting path.
- Live SC row probabilities are still refreshed from the rotamer node (`refresh_sc_row_probabilities_from_rotamer(...)`) and normalized per proxy into `sc_row_prob_norm`.
- SC-env and allowed same-residue SC-BB interactions still use those normalized probabilities as priors in the Boltzmann mixture (`z_shift += prior * exp(-shifted)`, `w = prior * exp(-shifted) / z_shift`); the new `sc_force_uncap_mix` only blends capped and uncapped pair-force vectors inside `eval_pair_force(...)`.

## 2026-02-28 (Production SC Force-Cap Transition Wiring)
- Before this change, `../../src/martini.cpp` parsed `sc_env_lj_force_cap`, `sc_env_coul_force_cap`, and `sc_env_relax_steps`, but the active probabilistic SC coupling path always applied fixed force caps; `sc_env_relax_steps` did not affect SC-env or same-residue SC-BB pair evaluation.
- The production SC coupling path now uses a per-activation runtime counter to compute an uncapped-force mix factor from `0.0 -> 1.0` over `sc_env_relax_steps` active production evaluations.
- The ramp is applied only to SC-env and allowed same-residue SC-BB LJ/Coulomb pair forces, by blending each capped pair-force vector with its uncapped counterpart; hard-sphere and non-SC pair paths are unchanged.

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

## 2026-03-02 (SC Cap-Removal Schedule Tightening)
- The active SC force-mix function in `../../src/martini.cpp` already reaches full uncapped force when `sc_env_transition_step >= sc_env_relax_steps - 1`; no runtime algorithm change was required for the user's request.
- Reducing the default `sc_env_relax_steps` from `200` to `150` is sufficient to leave at least the last `50` steps of the existing `200`-step SC-env/SC-BB interaction window on regular LJ/Coulomb forces.
- Updated defaults in both `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh`, and aligned the runtime fallback default in `../../src/martini.cpp` to avoid mismatched behavior when the HDF5 control attribute is absent.

## 2026-03-02 (Rigid-Hold / PO4 Clamp Audit After 150-Step Ramp Change)
- The production-stage `sc_env_transition_step` counter only drives the capped-to-uncapped force mix in `../../src/martini.cpp`; it does not gate protein rigid masks or PO4 z-clamp duration.
- Protein rigid hold is currently pre-production only: when hybrid becomes active (`current_stage == activation_stage`), `../../src/martini.cpp` clears both dynamic fixed-atom and dynamic z-fixed-atom masks instead of keeping backbone/protein coordinates frozen.
- The production PO4 z-clamp is not a true coordinate constraint and is not limited to 150 steps:
  - it stores each environment `PO4` atom's initial `z` coordinate after hybrid activation,
  - uses that stored `z` only when constructing SC-env interaction geometry,
  - zeroes only the `z` component of the SC-env gradient on those `PO4` atoms,
  - and remains enabled for the full active hybrid stage while `sc_env_po4_z_clamp_enabled=1`.

## 2026-03-02 (Active-Stage Backbone / PO4 Hold Implementation)
- Added explicit hybrid-control attributes for active-stage startup holds:
  - `sc_env_backbone_hold_steps` default `200`,
  - `sc_env_po4_z_hold_steps` default `150`.
- The user clarified that the startup hold must not freeze the active Upside backbone carriers or the refreshed MARTINI `BB` coordinates in space.
- Final semantics:
  - during the first `200` active SC transition steps, protein still evolves through the normal Upside `N/CA/C/O -> BB refresh -> RMSD align` path,
  - but BB and probabilistic SC coupling gradients are prevented from feeding back onto protein coordinates during that window,
  - so startup coupling is one-way onto environment rather than a hard coordinate freeze on protein.
- The PO4 hold atom set is derived from non-protein `atom_roles == "PO4"` and is enforced as an actual integrator/barostat z-fixed mask, not just an SC-interaction-space clamp.
- Step counting semantics:
  - `sc_env_transition_step` still increments from the SC coupling path during active production evaluations,
  - hold-mask refresh is now applied at the end of each MD step from `../../src/main.cpp`, so the configured hold remains active for the full current step and drops only before the next one.
- Active RMSD alignment remains enabled during the startup hold so new Upside backbone coordinates continue to define the refreshed active `BB` positions.
- Validation from a 2-step production smoke (`example/16.MARTINI/test_prod_run_sim_1rkl.sh`):
  - protein atoms (`protein_membership >= 0`): `2.064e-3 Å` max displacement between saved frames 0 and 1,
  - active `BB` proxy atoms: `3.586e-4 Å` max displacement,
  - active `N/CA/C/O` carrier targets from `hybrid_bb_map`: `2.064e-3 Å` max displacement,
  - `PO4` z displacement: `0.0 Å`,
  - `PO4` x/y displacement: nonzero (`1.268e-4 Å` max).

## 2026-03-03 (Production Startup Energy Build-Up Follow-Up)
- The startup cap-removal logic previously only reached the probabilistic SC row-evaluation path; deterministic protein-env BB pairs were still evaluated with `eval_pair_force(..., 0.f, 0.f, 1.f, ...)`, so they bypassed the production startup force-cap transition entirely.
- Advancing `sc_env_transition_step` only when `use_probabilistic_sc` is true was also too narrow for the intended startup-window semantics; the active production transition counter should advance for the full hybrid-active startup window.
- A cap-magnitude experiment did not move the early production trend:
  - replay with `SC_ENV_LJ_FORCE_CAP=50`, `SC_ENV_COUL_FORCE_CAP=50` matched the `25/25` replay through step `450`,
  - indicating the observed startup energy build-up was not being limited by the `25`-force cap threshold itself.
- A force-mix shape experiment also did not change the replay over the first `450` steps:
  - replacing the linear uncapping schedule with a quadratic ease-in produced identical logged potentials in the targeted `500`-step replay,
  - indicating the dominant startup issue was not the SC cap-release curve by itself.
- The effective lever was the protein feedback hold:
  - changing BB/SC gradient feedback onto protein from a hard `0` over the first `200` steps to a gradual `0 -> 1` ramp over `sc_env_backbone_hold_steps` strongly reduced the startup MARTINI potential.
- Targeted production-only replay results (`test_prod_run_sim_1rkl.sh`, `500` steps, same seed and starting checkpoint):
  - previous cap-only replay after deterministic BB-env capping:
    - steps `50/100/150/200/300/400/450`: MARTINI potential `4070/3877/4149/4449/5468/6788/7784`
  - new feedback-ramp replay:
    - steps `50/100/150/200/300/400/450`: MARTINI potential `3777/3366/3115/3046/2992/2767/2764`
- Interpretation:
  - the dominant startup energy accumulation was coming from keeping protein-environment back-reaction fully off while the environment was already responding to the refreshed Upside-driven backbone;
  - gradually restoring protein feedback over the same existing 200-step window dissipates that overlap energy much more effectively than only adjusting the startup SC force-cap parameters.
- Full-horizon confirmation from the `5000`-step replay (`outputs/martini_test_1rkl_hybrid_phase17_full/logs/stage_7.0.log`):
  - new MARTINI potential at `500/1000/1500/2000/2500/3000/3500/4000/4500`:
    `2962/2607/2295/2187/2370/1936/1973/1884/1812`;
  - user-reported original MARTINI potential:
    `7117/23902/54017/94657/135397/183071/237411/340015/453706`;
  - the new run remains in the `~1.8e3-3.0e3` range through the full stage instead of entering the previous monotonic runaway regime.
