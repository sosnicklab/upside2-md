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

## 2026-03-03 (PO4 Z-Hold Window Runtime Re-Check)
- Code path still reads as intended:
  - `active_sc_env_po4_z_hold_enabled(...)` uses `sc_env_transition_step < sc_env_po4_z_hold_steps`,
  - default `sc_env_po4_z_hold_steps` is `150`,
  - `refresh_transition_holds_for_engine(...)` is called at the MD step boundary from `../../src/main.cpp`.
- Runtime validation contradicts the intended `150`-step hold:
  - targeted production-only replay: `152` steps, `frame_steps=1`, run dir `outputs/martini_test_1rkl_hybrid_po4hold`;
  - checkpoint inspected: `outputs/martini_test_1rkl_hybrid_po4hold/checkpoints/1rkl.stage_7.0.up`;
  - non-protein `PO4` atom count: `280`.
- Measured `PO4` z displacement from saved frames:
  - frames `0..75`: max absolute `z` displacement `0.0 Å`,
  - first nonzero frame: `76`, max absolute `z` displacement `9.16e-05 Å`,
  - frame `150`: max absolute `z` displacement `1.0627e-01 Å`,
  - frame `151`: max absolute `z` displacement `1.0877e-01 Å`.
- Measured `PO4` xy displacement confirms these atoms are otherwise moving:
  - frame `150`: max in-plane displacement `2.2726e-01 Å`.
- Inference from the observed release at frame `76`:
  - the effective hold duration is about `75` MD steps, not `150`;
  - this is consistent with `sc_env_transition_step` being consumed roughly twice per MD step in the active production integrator path, though that still needs a code-level root-cause fix.

## 2026-03-03 (PO4 Z-Hold Counter Fix Validation)
- Root cause was the transition counter location, not the hold predicate:
  - `sc_env_transition_step` was still incrementing inside active MARTINI force evaluation,
  - active force evaluation happens more than once per MD step,
  - so the `150`-step startup windows were being consumed in force calls rather than completed MD steps.
- Runtime fix:
  - removed the active-stage `sc_env_transition_step += 1` from MARTINI force evaluation,
  - moved the increment into `refresh_transition_holds_for_engine(...)`, which is called from `../../src/main.cpp` at the MD step boundary.
- Post-fix targeted validation:
  - replay: `152` production steps with `frame_steps=1`,
  - run dir: `outputs/martini_test_1rkl_hybrid_po4hold_fix`,
  - checkpoint inspected: `outputs/martini_test_1rkl_hybrid_po4hold_fix/checkpoints/1rkl.stage_7.0.up`,
  - non-protein `PO4` atom count: `280`.
- Measured `PO4` displacement after the fix:
  - max absolute `z` displacement through frame `150`: `0.0 Å`,
  - max absolute `z` displacement at frame `149`: `0.0 Å`,
  - max absolute `z` displacement at frame `150`: `0.0 Å`,
  - first nonzero `z` displacement at frame `151`: `9.155e-05 Å`,
  - max in-plane displacement at frame `150`: `2.2728e-01 Å`.
- Conclusion:
  - bilayer `PO4` `z` is now held exactly constant for the first `150` MD steps of the `200`-step SC-env-BB startup window,
  - release begins on the following saved step/frame (`151`), while `x/y` motion remains unconstrained throughout the hold.

## 2026-03-04 (Pre-7.0 Workflow Parity Audit vs Rigid Dry)
- The stage `6.0-6.6` driver blocks in `run_sim_1rkl.sh` and `run_sim_1rkl_rigid_dry.sh` are identical for:
  - stage ordering,
  - `prepare_stage_file(...)` arguments,
  - NPT target selection,
  - minimization/MD invocation,
  - stage-to-stage handoff sequence,
  - pre-production VTF extraction mode.
- Existing generated checkpoints confirm no physical drift before stage `7.0`:
  - `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.0.up` matches `outputs/martini_test_1rkl_rigid_dry/checkpoints/1rkl.stage_6.0.up` in final coordinates and box exactly;
  - same exact match for stages `6.2` and `6.6`.
- The only pre-7.0 file-level differences found were in hybrid-control metadata and inert AA carrier masses:
  - hybrid stage files carried `activation_stage="production"` inherited from mapping injection, while rigid-dry explicitly overwrote pre-production stages to `activation_stage="__hybrid_disabled__"`;
  - both workflows already had `preprod_protein_mode="rigid"`;
  - injected `PROTEINAA` carrier masses differ (`124` atoms; hybrid uses scaled masses summing to MARTINI `BB` mass, rigid-dry still uses the older unscaled masses), but all non-`PROTEINAA` masses are identical and the pre-7.0 coordinates/boxes are unchanged.
- Practical implication:
  - pre-7.0 physical state is already the same between the two workflows;
  - the hybrid script still benefited from an explicit pre-production activation override so stage-file semantics match the rigid workflow instead of depending on the current `"minimization"` stage-label mismatch to keep hybrid inactive.

## 2026-03-04 (Stage-7 Initial Energy Mismatch Root Cause)
- The stage-7 initial-energy mismatch is not caused by pre-production drift or by different `6.6 -> 7.0` handoff coordinates:
  - the hybrid and rigid-dry stage `6.0`, `6.2`, and `6.6` checkpoints already match exactly in coordinates and box;
  - the stage-7 prepared files also have identical `/input/pos` on the shared MARTINI coordinate array before production-specific runtime changes act.
- The hybrid stage-7 file carries a different potential graph from the rigid-dry stage-7 file:
  - hybrid has additional production Upside nodes under `/input/potential`: `affine_alignment`, `backbone_pairs`, `hbond_energy`, `infer_H_O`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, `protein_hbond`, `rama_coord`, `rama_map_pot`, `rama_map_pot_ref`, `rotamer`;
  - rigid-dry has only `angle_spring`, `dihedral_spring`, `dist_spring`, and `martini_potential`.
- Reconstructed handoff test from the same hybrid `stage_6.6` checkpoint:
  - copied the hybrid `stage_7.0.prepared.up`,
  - applied the same `set_initial_position.py` handoff used by `run_sim_1rkl.sh` (`STRICT_COPY=1`, `REFRESH_HYBRID_CARRIERS=1`, no recenter),
  - then ran two 1-step probes from those same coordinates.
- Probe A: hybrid stage-7 file unchanged (`activation_stage="production"`, `hybrid_active=1`)
  - initial energy: `Upside/MARTINI/Total = 663.56 / 4268.90 / 4932.46`
  - matches the user-observed hybrid log exactly.
- Probe B: same handed-off hybrid stage-7 file, but only `activation_stage` changed to `__hybrid_disabled__` (`hybrid_active=0`)
  - initial energy: `Upside/MARTINI/Total = 663.56 / -19534.69 / -18871.13`
  - MARTINI term matches the rigid-dry log exactly.
- Conclusion:
  - the `+663.56` `potential` term comes from the extra Upside production nodes present only in the hybrid stage-7 file;
  - the MARTINI sign flip (`-19534.69 -> +4268.90`) is caused entirely by the active hybrid runtime path, not by coordinates.
- The active hybrid runtime path changes stage-7 MARTINI energetics in ways that remove the large negative cohesive terms present in the rigid-dry workflow:
  - most intra-protein MARTINI pairs are filtered out when hybrid is active (`BB-BB` off, `SC-SC` off, `BB-SC` only same residue);
  - non-protein/non-protein MARTINI nonbonded is switched to hard-sphere-like repulsion while hybrid is active.

## 2026-03-19 (Stage-7 Backbone Rigid Hold via Workflow Fix-Rigid Mask)
- `example/16.MARTINI/run_sim_1rkl.sh` can keep the stage-7 backbone rigid without disabling hybrid production by writing `/input/fix_rigid` directly into the prepared production `.up` file after hybrid node augmentation.
- The smallest reliable stage-7 selection is driven by prepared-file metadata, not hardcoded indices:
  - require `hybrid_env_topology/protein_membership >= 0` to restrict the mask to protein atoms,
  - require atom roles/names `BB`, `N`, `CA`, `C`, and `O` so both dry-MARTINI backbone beads and injected all-atom backbone carriers are included.
- Existing `/input/fix_rigid/atom_indices` should be preserved and unioned with the new backbone mask rather than overwritten silently.
- A fast-fail validation is worthwhile in the workflow helper:
  - missing any one of `BB/N/CA/C/O` should abort stage preparation,
  - per-role counts should match the `BB` count so partial carrier injection does not produce a silently incomplete rigid mask.

## 2026-03-19 (Initial Diagnostic: Rigid Backbone + Active Hybrid Stage-7 Failure)
- This diagnosis was incomplete. Later Phase 22 probes show that rigid backbone + active hybrid is stable when `production_nonprotein_hard_sphere=0`; the rigid mask itself is not the root cause.
- Using the archived stable phase-17 prepared stage-7 file plus the archived stage-6.6 handoff, adding only the new stage-7 backbone `fix_rigid` mask reproduces the user's accumulating-energy trace exactly:
  - initial `Upside/MARTINI/Total = 663.56 / 4268.90 / 4932.46`,
  - step `500`: `663.56 / 3471.02 / 4134.58`,
  - step `1000`: `663.56 / 5938.75 / 6602.31`.
- In that reproduced run, protein observables stay flat (`potential`, `Rg`, hbonds), so the pump is not coming from changing Upside backbone structure or from handoff drift; it is the active stage-7 MARTINI kernel acting on a protein that is no longer allowed to respond.
- Rigid protein alone is stable:
  - archived `run_sim_1rkl_rigid_dry.sh` stage-7 output remains monotone and bounded (`martini_potential -19534.69 -> -22112.04` over `0 -> 4500`).
- Disabling active integration RMSD alignment does not materially change the reproduced rigid-hybrid trace at the first sampled checkpoint (`3464.62` vs `3471.02` at step `500`), so alignment is not the dominant cause.
- Practical conclusion:
  - active stage-7 hybrid production assumes protein `BB/N/CA/C/O` can absorb the BB/SC coupling it applies;
  - if those backbone DOFs are fixed, the physically consistent workflow is to keep the explicit backbone `fix_rigid` mask but disable stage-7 hybrid activation, not to leave `hybrid_active=1`.
- Validated fixed control combination from the same archived handoff:
  - set `activation_stage="__hybrid_disabled__"` and `preprod_protein_mode="free"` while keeping the backbone `fix_rigid` mask,
  - resulting energy trend is stable and decreasing: initial `663.56 / -19534.69 / -18871.13`, step `500` `663.56 / -22979.09 / -22315.53`, step `1000` `663.56 / -24867.03 / -24203.48`.

## 2026-03-20 (Behavior Lesson: Do Not Treat Hybrid Disablement as a Fix)
- User correction pattern:
  - this repository is explicitly about hybrid simulation, so a result that becomes stable only after disabling hybrid is not a valid fix.
- Working rule:
  - when a workaround turns off the defining system feature, record it only as diagnostic evidence and immediately broaden the investigation into `src/` and the system design before calling the issue fixed.
- Applied here:
  - the rigid-backbone stage-7 `activation_stage="__hybrid_disabled__"` path remains useful to isolate the failing branch,
  - but the real task is now to identify which active-hybrid runtime path or Hamiltonian change in `src/` causes the energy pump and repair that while keeping hybrid active.

## 2026-03-20 (Phase 22 Runtime Audit: Fixed Atoms vs Active Hybrid)
- The current fixed-atom implementation does not make protein coordinates globally immutable once active hybrid is on:
  - `../../src/martini.cpp::apply_fix_rigid_md(...)` only zeros derivatives and momenta for fixed atoms.
  - `../../src/deriv_engine.cpp` skips fixed atoms inside the integrator update loops, but `../../src/martini.cpp::align_active_protein_coordinates(...)` still rotates/translates all protein atoms before force evaluation, and `refresh_bb_positions_if_active(...)` still overwrites `BB` coordinates from mapped carriers.
- Therefore the earlier "one-way forcing into fully frozen atoms" explanation is incomplete:
  - fixed atoms are protected against integrator updates and barostat scaling,
  - but they are not protected against active-hybrid coordinate rewrites that happen outside the integrator.
- The active-hybrid Hamiltonian also changes much more than the earlier rigid-backbone explanation captured:
  - `allow_intra_protein_pair_if_active(...)` removes most intra-protein MARTINI bonded/nonbonded terms when hybrid is active.
  - `production_nonprotein_hard_sphere=1` replaces all non-protein/non-protein MARTINI nonbonded interactions with WCA-like repulsion in production.
- On the exact same archived rigid stage-7 handoff, keeping `hybrid_active=1` but setting `exclude_intra_protein_martini=0` and `production_nonprotein_hard_sphere=0` changes the initial MARTINI energy from `+4268.90` to `-23653.72`.
- Working conclusion for Phase 22:
  - the rigid-backbone failure cannot be blamed only on projected forces landing on fixed atoms;
  - active-hybrid coordinate rewrites and the active-hybrid Hamiltonian edits both need to be audited as potential root causes.

## 2026-03-20 (Phase 22 Root Cause: Non-Protein Hard-Sphere Branch)
- The rigid stage-7 control matrix on the exact same archived `6.6 -> 7.0` handoff isolates `production_nonprotein_hard_sphere` as the dominant failing branch:
  - baseline rigid active hybrid (`exclude_intra=1`, `nonprotein_hs=1`): `4268.90 -> 3471.02 -> 5938.75` at `0/500/1000`;
  - rigid active hybrid with only `exclude_intra=0`: `4228.57 -> 3448.09 -> 6024.10`;
  - rigid active hybrid with only `nonprotein_hs=0`: `-23613.40 -> -23939.15 -> -21299.20`;
  - rigid active hybrid with both toggles off: `-23653.72 -> -23962.09 -> -21213.85`.
- Interpretation:
  - disabling intra-protein MARTINI exclusions is not the important change for this failure mode;
  - replacing all non-protein/non-protein MARTINI LJ+Coulomb interactions with repulsive-only WCA in active hybrid production removes the environment's cohesive energy and drives the observed MARTINI energy pump.
- The same conclusion survives outside the rigid diagnostic:
  - with hybrid still active, no rigid mask, and `nonprotein_hs=0`, the non-rigid stage-7 probe stayed stable through step `1000` with `Rg ~12.9 A`, `martini_potential -24897.98` at step `500`, and `-25198.87` at step `1000`.
- Implementation rule:
  - `production_nonprotein_hard_sphere` should default to `0` in the runtime and generated hybrid inputs.
  - workflow stage preparation should keep hybrid active and write `production_nonprotein_hard_sphere=0` explicitly instead of disabling hybrid when the backbone rigid mask is used.
- Verification:
  - after the runtime default change, deleting the `production_nonprotein_hard_sphere` attribute entirely from a production file still produced `nonprotein_hs=0` in the runtime parse log and the expected negative initial MARTINI energy (`-23613.40`).
