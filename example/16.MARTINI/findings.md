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

## 2026-04-02 (SC-Table Orientation Contract Audit)
- `SC-training` first-pass semantics are explicitly isotropic:
  - `SC-training/README.md` says sampled target positions lie on spherical shells around an effective sidechain center and that sidechain geometry is not reconstructed, so directional samples remain isotropic within a shell.
  - `SC-training/workflow.py` writes the same assumption into the manifest as `effective_model = "sum_beadwise_colocated_spherical_shells"`.
- The assembled training output contains directional sample metadata, but the runtime library intentionally discards it:
  - `example/16.MARTINI/build_sc_martini_h5.py` reads only `grid_nm` and `energy_kj_mol` from each residue/target entry and writes only those radial arrays into `martini.h5`.
  - therefore the runtime SC-table library contains no sidechain orientation vector, normal vector, or anisotropic angular dependence.
- The runtime stage-7 SC-table path also performs a purely radial lookup:
  - `example/16.MARTINI/inject_sc_table_stage7.py` creates `placement_fixed_point_only_CB` from `affine_alignment` plus a fixed `CB_PLACEMENT`.
  - `../../src/martini.cpp::MartiniScTablePotential` computes `dr = cbp - envp`, reduces it to `dist = |dr|`, and evaluates the spline only in that scalar radius.
- Conclusion:
  - there is no trained sidechain-normal orientation term to preserve in the current first-pass SC table;
  - runtime usage is consistent with that isotropic training contract.
- Important caveat:
  - the affine frame still matters geometrically because it places the `CB` proxy in space from `N/CA/C`, including an out-of-plane component;
  - however, that affects only where the isotropic center sits, not the shape of the trained forcefield.
  - the training workflow names the origin an "effective sidechain center", while runtime currently uses `CB` as that anchor; that is a center-choice assumption, not an orientation inconsistency.

## 2026-04-02 (User Correction: MARTINI SC Table Must Follow Upside Orientation Contract)
- The isotropic audit above is not the final design target; it only describes the incomplete baseline currently in the repository.
- The user clarified the intended contract:
  - original Upside sidechain potentials include vector terms for sidechain orientation,
  - the MARTINI sidechain training/export/runtime path must do the same,
  - the only model change versus original Upside is the interaction source: dry-MARTINI sidechain beads against dry-MARTINI particle types, yielding one table per amino-acid residue type and dry-MARTINI target type.
- The original Upside orientation source is already present locally and matches the user-provided ConDiv reference:
  - `parameters/ff_2.1/sidechain.h5` and `/Users/yinhan/Documents/upside2-md-ConDiv/ConDiv/remd-4000-8RP-1th-test/init_param/sidechain.h5` have identical `rotamer_center_fixed` contents.
  - `py/upside_config.py::write_environment(...)` constructs `placement_fixed_point_vector_only_CB` from the canonical `CB` position and a `CA -> CB` unit vector by default.
  - `src/environment.cpp` consumes that 6D `CB` placement through `dp = dot(displace_unitvec, rvec1)`, confirming the original environment-side sidechain model is orientation-aware through a single sidechain vector.
- The implementation target for the MARTINI SC table is therefore:
  - train in the original residue-local `CB` frame,
  - store a `distance x cos(theta)` table in `martini.h5`,
  - inject/use `placement_fixed_point_vector_only_CB` at stage 7,
  - evaluate the runtime table against both `CB` point and `CB` vector, not radius alone.

## 2026-04-02 (Behavior Lesson: Do Not Freeze a Baseline Audit Into the Final Contract)
- User correction pattern:
  - an audit of the current implementation is not permission to preserve that implementation when the scientific intent is different.
- Working rule:
  - when the user points to an original reference workflow or forcefield contract, treat the current code as a baseline to compare against, then implement the reference contract rather than optimizing around the baseline shortcut.

## 2026-04-02 (User Correction: Required MARTINI SC Functional Form)
- The required one-sided sidechain/dry-MARTINI potential is not an arbitrary 2D surface.
- The user-specified contract is:
  - `V_hybrid(r,theta) = V_radial(r) + ang_1(-n_1 . n_12) * V_angular(r)`.
- This matches the original Upside SC-SC interaction structure after removing the second orientation-dependent factor for the dry-MARTINI particle, which has no orientation state of its own.
- Implementation consequence:
  - training may still sample the full directional energy surface,
  - but the stored/exported/runtime forcefield must be factorized into a radial baseline term, a one-dimensional angular profile `ang_1`, and a radial angular-amplitude term `V_angular(r)`.

## 2026-04-02 (Original Upside SC-SC Factorization Source)
- The original SC-SC training form is explicitly factorized in `/Users/yinhan/Documents/upside2-md-ConDiv/py/rotamer_parameter_estimation.py::quadspline_energy(...)`:
  - spline parameters are split into `dp1`, `dp2`, `uni`, and `direc`,
  - the returned energy tensor is `uni + dp1 * dp2 * direc`.
- Interpretation for the MARTINI hybrid workflow:
  - `uni` is the radial baseline term,
  - `dp1` and `dp2` are the two partner-orientation factors,
  - `direc` is the radial directional-amplitude term.
- Therefore the user-specified dry-MARTINI reduction is consistent with original Upside:
  - because the dry-MARTINI target particle has no orientation state, the second partner factor is removed,
  - leaving the one-sided form `V_radial(r) + ang_1(-n_1 . n_12) * V_angular(r)`.

## 2026-04-02 (User Clarification: n1 Must Be Backbone-Defined)
- The intended `n_1` vector is the sidechain orientation with respect to the protein backbone, not a direction inferred from the placed dry-MARTINI target particle.
- The current implementation already matches that contract:
  - `SC-training/workflow.py` computes the angular coordinate as `-dot(direction, cb_vector_unit)`, where `cb_vector_unit` is the canonical backbone-defined `CA -> CB` vector and `direction` is only the sampled `CB -> target` direction.
  - `src/martini.cpp::MartiniScTableOneBody` uses the same runtime-side `CB` vector from `placement_fixed_point_vector_only_CB` and evaluates `angular_coord = -dot(displace_unitvec, cbv)`.
- Practical consequence:
  - training and runtime share the same definition of `n_1`,
  - the dry-MARTINI particle contributes only `n_12`,
  - no separate target-defined orientation state is needed at runtime.

## 2026-04-02 (SC/Dry-MARTINI Probabilistic Weighting Audit)
- Original Upside static sidechain weighting is carried by the rotamer solver, not by a post-hoc average:
  - `/Users/yinhan/Documents/upside2-md-ConDiv/py/upside_config.py::write_rotamer_placement(...)` writes fixed rotamer one-body energies through `placement_fixed_scalar` from `rotamer_prob_fixed`.
  - `../../src/rotamer.cpp` consumes one-body inputs per rotamer state and combines them with pair interactions through the standard belief-propagation `rotamer` node.
- In this repo, the local sidechain library does not contain `rotamer_prob_fixed`; it contains only `rotamer_prob`.
  - Averaging `parameters/ff_2.1/sidechain.h5::rotamer_prob` over its two backbone-context axes already yields per-residue rotamer probabilities whose residue-local sums are `~1.0`, so it is a valid fixed-prior fallback for static placement/training.
- Current MARTINI SC path matches that Upside weighting contract:
  - `SC-training/workflow.py` derives `sidechain_rotamer_weight` from the same residue-local mean `rotamer_prob` fallback and uses it when assembling the residue-averaged training target.
  - `example/16.MARTINI/inject_sc_table_stage7.py` builds `placement_fixed_scalar` from the same residue-local mean `rotamer_prob` fallback when `rotamer_prob_fixed` is absent.
  - `example/16.MARTINI/inject_sc_table_stage7.py` wires the dry-MARTINI table into `rotamer(arguments=[placement_fixed_point_vector_only, placement_fixed_scalar, martini_sc_table_1body])`.
  - `../../src/martini.cpp::MartiniScTableOneBody` therefore supplies per-rotamer one-body energies to the existing Upside `rotamer` solver instead of applying an external deterministic average.
- Numerical verification on a smoke-built temporary library:
  - derived fixed priors from `sidechain.h5::rotamer_prob`,
  - exported `martini.h5::rotamer_probability_fixed`,
  - recovered probabilities from injected `placement_fixed_scalar` via `exp(-E)` and residue-local normalization,
  - all agreed to floating-point tolerance, with worst-case max-abs mismatch below `4e-7`.
- Operational verification:
  - temporary `martini.h5` rebuild and stage-7 injection smoke both succeeded;
  - zero-duration `obj/upside` runtime probe on the injected file also succeeded after restoring missing legacy unit-conversion attrs on the copied probe stage file.
- Conclusion:
  - yes, sidechain/dry-MARTINI interactions need probabilistic weighting if the goal is parity with Upside sidechain treatment;
  - the current implementation now uses that same weight setting by routing the hybrid SC term through the standard rotamer machinery rather than bypassing it.

## 2026-04-02 (Stage-7 Crash Root Cause: Stale SC Artifacts)
- The reported `inject_sc_table_stage7.py` crash on missing `rotamer_count` was caused by stale default training artifacts, not by a bad HDF5 read path in the injector.
- Verified root cause:
  - `SC-training/runs/default/results/assembled/sc_table.json` was still `sc_training_table_v2` and lacked per-rotamer keys like `rotamer_count`.
  - `parameters/ff_2.1/martini.h5` therefore also lacked `rotamer_count`, `rotamer_probability_fixed`, and the per-rotamer radial/angular datasets.
  - `example/16.MARTINI/run_sim_1rkl.sh::ensure_sc_martini_library()` and the production-only test script only checked file existence, so they allowed the stale `martini.h5` to pass through to stage-7 injection.
- Fix applied:
  - both workflow scripts now default `SC_MARTINI_TABLE_JSON` to `SC-training/runs/default/results/assembled/sc_table.json`,
  - validate the required rotamer-resolved `martini.h5` datasets before production injection,
  - rebuild `martini.h5` automatically when the file is missing or stale.
  - `example/16.MARTINI/inject_sc_table_stage7.py` now surfaces a direct schema/rebuild error instead of failing later with a raw missing-dataset exception.
- Artifact repair:
  - reran the full default `SC-training` workflow (`init-run`, then `run-local --force`) so the assembled table now has schema `sc_training_table_v3`.
  - rebuilt `parameters/ff_2.1/martini.h5` from that refreshed assembled table; it now includes the required rotamer-resolved datasets.
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh` passed.
  - `bash -n example/16.MARTINI/test_prod_run_sim_1rkl.sh` passed.
  - `python3 -m py_compile example/16.MARTINI/inject_sc_table_stage7.py` passed.
  - stage-7 injection smoke using the rebuilt repo `martini.h5` succeeded.

## 2026-04-02 (Behavior Lesson: Validate Artifact Schema, Not Just Timestamps or Existence)
- User correction pattern:
  - claiming an artifact is "updated" based on file presence or modification time is not good enough when runtime behavior depends on specific datasets/schema.
- Working rule:
  - for generated forcefield artifacts, verify the required schema directly before concluding they are current or before wiring them into runtime.
  - if a workflow depends on a generated HDF5 artifact, the guard should validate required datasets, not only file existence.

## 2026-04-02 (SC-Training Orientation Vector / No-Backbone Audit)
- The SC-training angular coordinate uses the same default vector construction as Upside's environment-side `CB` placement:
  - `SC-training/workflow.py` defines `CANONICAL_CB_VECTOR_UNIT` from the canonical `CB` position `(0.0, 0.94375626, 1.2068012)` and uses it in `cos_theta = -dot(direction, cb_vector_unit)`.
  - `py/upside_config.py::write_environment(...)` and the ConDiv reference both construct `placement_fixed_point_vector_only_CB` with the same canonical `CB` point and the default `CA -> CB` unit vector `(ref_pos[3] - ref_pos[1]) / |ref_pos[3] - ref_pos[1]|`.
  - Direct numeric comparison gives exact agreement (`max_abs_diff = 0.0`) between the training vector and the default Upside vector.
- Important interpretation:
  - `SC-training/workflow.py` also loads per-rotamer vectors from `sidechain.h5::rotamer_center_fixed[:,3:6]`, but those are carried as metadata and are not used to define the one-sided `n_1` orientation coordinate.
  - The actual `n_1` used in the training surface is backbone-defined, matching the runtime/UpSide contract.
- Backbone particles are not included in SC-training pairs:
  - the training bead list comes from `_load_martini_forcefield(...)`, which reads martinize `ff.sidechains` entries for each residue rather than any backbone bead definition.
  - Current residue sidechain bead sets are sidechain-only, and a full check over the default residue map found no `BB` bead in any residue.
  - The refreshed default training manifest uses sidechain bead types `{AC1, AC2, C3, C5, N0, P1, P4, P5, Qa, Qd, SC4, SC5, SNd, SP1}`; `BB` is absent.

## 2026-04-02 (Hybrid BB Force Transfer Audit)
- The active-stage `BB-BB` exclusion is already stronger than the user requirement:
  - `src/martini.cpp::allow_protein_pair_by_rule(...)` currently returns `false` for any intra-protein MARTINI proxy pair while hybrid is active.
  - That gate feeds `allow_intra_protein_pair_if_active(...)` / `allow_multibody_term_if_active(...)`, so it suppresses intra-protein proxy nonbonded pairs, bond terms, and angle terms, not only `BB-BB`.
- The actual force-transfer bug was in the BB protein-environment pair path:
  - `src/martini.cpp::MartiniPotential::compute_value(...)` was evaluating active-stage protein `BB` interactions at `CA` by replacing `p1/p2` with `direct_ca_atom_for_bb_proxy(...)`.
  - It also wrote the protein-side BB force directly onto `CA`, so BB coupling bypassed the hybrid BB carrier map and did not use the same live BB proxy geometry that is refreshed from the backbone carriers.
- The corrected runtime contract is now:
  - refresh live `BB` proxy coordinates before MARTINI pair evaluation,
  - evaluate `BB` environment interactions at the BB proxy coordinate itself,
  - project the resulting protein-side BB gradient back through `hybrid_bb_map` weights to the mapped backbone carriers (`N/CA/C/O` when present),
  - let those projected BB forces sum naturally with SC table point/vector forces on the shared backbone DOFs.

## 2026-04-03 (Should BB Be Merged Into SC Training?)
- `example/16.MARTINI/martinize.py::martini22` defines backbone typing by secondary structure plus a small residue-specific override table:
  - default `BB` types by secondary structure are `N0, Nda, N0, Nd, Na, Nda, Nda, P5, P5` for `F,E,H,1,2,3,T,S,C` respectively.
  - only `ALA`, `PRO`, and `HYP` override that default pattern in `bbtyp`.
  - the actual runtime lookup is `bbGetBead(residue, ss)`, which first checks residue-specific overrides and otherwise falls back to the default state-dependent bead.
- Exact count for the forcefield used here (`martini22`):
  - unique backbone bead types used by `bbGetBead(...)`: `C5, N0, Na, Nd, Nda, P4, P5`.
  - therefore `n = 7` backbone dry-MARTINI bead types.
  - across the 20 canonical amino acids, those assignments collapse to only 3 residue-pattern families:
    - default pattern for 18 residues,
    - `ALA` pattern,
    - `PRO` pattern.
- Consequence for training design:
  - a `20 x n = 140` residue-specific BB table family is not how dry-MARTINI backbone nonbonded identity is defined in `martinize.py`.
  - if the goal is to stay faithful to the FF being used in simulation, BB training should be keyed by the BB bead type actually assigned by MARTINI, not by amino-acid identity on top of that.
  - with the current 38 dry-MARTINI target particle types in the assembled hybrid library, that means:
    - bead-type-faithful BB replacement would be `7 x 38 = 266` BB-target tables,
    - forcing a residue-specific `20 x 7 x 38 = 5320` table family would be a much larger and scientifically weaker parameterization unless there is a deliberate non-MARTINI reason to add amino-acid-specific effective BB terms.
- Recommendation:
  - do not fold BB into the existing SC-training workflow as one merged SC+BB table contract.
  - if learned BB terms are desired, add a separate BB-training workflow keyed by MARTINI BB bead type (and its own local backbone frame if orientation is later introduced), while keeping SC training residue/rotamer-based.

## 2026-04-03 (VMD NewCartoon Failure on Exported VTF)
- Root cause of the VMD error
  - `example/16.MARTINI/extract_martini_vtf.py` was writing VTF atom lines as only `atom <i> name <aname>`.
  - Without `resid/resname/chain/segid`, VMD assigned fallback residue metadata `X 0 X` to every atom when building the temporary PDB used for `NewCartoon`/Stride.
  - That caused Stride to see one huge fake residue and fail with:
    - `too many atoms in residue X 0 X`
    - followed by missing Stride output / cartoon generation failure.
- Secondary metadata issue found during the audit
  - mode-1/mode-2 backmapped protein backbone atoms were being labeled with residue name `PRO` for every residue.
  - The prepared stage file already carries `/input/sequence`, and for the tested stage file its length matches the unique `hybrid_bb_map/bb_residue_index` order, so the real residue names can be restored during export.
- Exporter fix
  - VTF atom records now include:
    - `name`, `resid`, `resname`, `segid`, and `chain`.
  - Backmapped protein backbone residue names are now derived from `/input/sequence` via the ordered `hybrid_bb_map/bb_residue_index` mapping instead of being hard-coded to `PRO`.
  - PDB export was updated in parallel to include the same chain-aware residue metadata.
- Verified output on the stage-7 prepared file:
  - VTF protein lines now look like:
    - `atom ... name N resid 4 resname ASP segid sA chain A`
    - `atom ... name CA resid 4 resname ASP segid sA chain A`
  - so VMD should now see standard backbone residues rather than the fallback anonymous residue.

## 2026-04-03 (Workflow Cleanup Dependency Audit)
- `example/16.MARTINI/run_sim_1rkl.sh` direct local-script dependencies are:
  - `prepare_system.py`
  - `martinize.py`
  - `build_sc_martini_h5.py`
  - `inject_sc_table_stage7.py`
  - `validate_hybrid_mapping.py`
  - `set_initial_position.py`
  - `extract_martini_vtf.py`
- The only extra Python module dependency within that kept set is:
  - `prepare_system.py -> prepare_system_lib.py`
- Workflow data/assets that must remain checked in for that path are:
  - `ff_dry/`
  - `pdb/`
- Not required as checked-in workflow dependencies:
  - alternate run scripts,
  - test/scan/plot scripts,
  - debug spline text files,
  - cached `__pycache__` / tool-state directories,
  - generated `inputs/` and `outputs/` directories.

## 2026-04-03 (Behavior Lesson: Confirm Exact Retention Contract Before Cleanup)
- User correction pattern:
  - when a cleanup request gives a structural rule like "only allowed files", do not infer a broader or narrower file budget; restate and implement the exact allowed set.
- Working rule:
  - for file-count or retention constraints, mirror the user's exact allowance literally and update the task tracker before deleting or consolidating files.
  - do not remove tracking files like `task_plan.md` unless the user explicitly includes them in the deletion scope.

## 2026-04-03 (Four-File Python Consolidation Result)
- The allowed Python layout under `example/16.MARTINI` is now:
  - `martinize.py`
  - `prepare_system.py`
  - `extract_martini_vtf.py`
  - `prepare_system_lib.py`
- Consolidation design that preserved the workflow:
  - `prepare_system.py` remains the simulation-input generation entry point and now dispatches the former auxiliary CLIs via subcommands.
  - `prepare_system_lib.py` now contains the logic previously split across:
    - `build_sc_martini_h5.py`
    - `inject_sc_table_stage7.py`
    - `set_initial_position.py`
    - `validate_hybrid_mapping.py`
  - `run_sim_1rkl.sh` now calls:
    - `prepare_system.py build-sc-martini-h5`
    - `prepare_system.py inject-stage7-sc`
    - `prepare_system.py set-initial-position`
    - `prepare_system.py validate-hybrid-mapping`
- Verification outcome:
  - only four `.py` files remain at top level in `example/16.MARTINI`;
  - no retained script still references the deleted Python filenames;
  - retained scripts parse cleanly.

## 2026-04-03 (Legacy Spline/Force Debug File Removal)
- Runtime source of the legacy debug files was isolated to `src/martini.cpp`:
  - `DihedralSpring` wrote `dihedral_splines.txt`;
  - `DistSpring` wrote `bond_splines.txt`;
  - `AngleSpring` wrote `angle_splines.txt`;
  - `MartiniPotential` wrote `all_splines.txt` and `force_debug.txt`.
- Those file outputs were controlled by HDF5 attrs that are no longer needed:
  - `debug_mode`
  - `force_debug_mode`
  - `overwrite_spline_tables`
- Stage-file generation was still enabling those attrs in `example/16.MARTINI/prepare_system_lib.py`, and `example/16.MARTINI/run_sim_1rkl.sh` still exported `UPSIDE_OVERWRITE_SPLINES`; both hooks are now removed.
- After cleanup, source-code references to those debug filenames remain only in tracker/history markdown, not in the runtime or preparation code.
