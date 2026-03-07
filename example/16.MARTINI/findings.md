# Findings

## 2026-03-06 (Primary-Box PBC Containment in Current `run_sim_1rkl.sh`)
- Fresh current-code rerun root cause split across three boundaries:
  - stage conversion copied boxed runtime-PDB coordinates verbatim into `.up`, including atoms exactly on the upper box face (`coord == L`), which are outside the half-open primary box convention `[0, L)`
  - stage handoff copied the previous stage's sampled coordinates verbatim into the next stage input and refreshed AA backbone carriers without wrapping afterward
  - `src/main.cpp` logged raw `engine.pos->output` into `/output/pos`, so sampled trajectories could drift slightly negative/above-box under PBC and later show very large out-of-box values once the known `6.4+` instability kicked in
- Fresh pre-fix proof from `example/16.MARTINI/outputs/test_pbc_probe_current/`:
  - `stage_6.0.up` had `3` out-of-box atoms already in `/input/pos`
  - `stage_6.0.up` last sampled `/output/pos` frame had `12` out-of-box atoms
  - `stage_6.2.up` last sampled frame had `19` out-of-box atoms
  - `stage_6.4.up` last sampled frame had `3997` out-of-box atoms while the separate MARTINI instability exploded coordinates
- Fix implemented:
  - added `wrap_xyz_array_to_box(...)` in `lib.py`
  - applied it in `lib.py::convert_stage()` before writing boxed `.up` input coordinates
  - applied it in `prepare_system.py::set_initial_position()` after handoff/backbone refresh and before writing `/input/pos`
  - changed the `pos` logger in `src/main.cpp` to wrap sampled coordinates against the live NPT box when present, otherwise against the static MARTINI box attrs
- Build note:
  - `cmake --build obj -j 4` initially failed because generated link rules still referenced `/opt/homebrew/opt/hdf5/lib/libhdf5.320.0.0.dylib`
  - `cmake -S src -B obj` refreshed the build files to the installed HDF5 version and the rebuild then succeeded
- Fresh post-fix proof from `example/16.MARTINI/outputs/test_pbc_wrapped/`:
  - `stage_6.0.up`, `stage_6.2.up`, `stage_6.4.up`, `stage_6.6.up`, and `stage_7.0.up` each had `0` out-of-box atoms in both `/input/pos` and the last sampled `/output/pos`
  - exported `1rkl.stage_7.0.vtf` had `0` out-of-box atoms in both frames
  - the separate late-stage physical instability still exists, but it no longer escapes the sampled/exported primary-box coordinates in this workflow

## 2026-03-06 (Three-Script Workflow Consolidation)
- The active workflow dependency surface had already shrunk to six Python helper entrypoints:
  - stage handoff
  - backbone-only node injection
  - backbone-only validation
  - depth-table build
  - backbone cross-table build
  - backbone cross-node injection
- Those helpers were only invoked from `run_sim_1rkl.sh` and `run_sim_bilayer.sh`; the remaining Python files in `example/16.MARTINI` were not part of the live workflow path.
- Cleanest consolidation path:
  - rename `prepare_system_lib.py` to `lib.py`
  - keep structural prep/stage conversion in `lib.py`
  - expose the six live helper behaviors as `prepare_system.py` subcommands
  - delete the standalone helper scripts
- Reduced rerun proof from `example/16.MARTINI/outputs/test_three_script_cleanup/`:
  - `prepare_system.py build-depth-table` wrote the depth CSV/JSON successfully
  - `prepare_system.py build-backbone-cross-table` rebuilt `parameters/ff_2.1/martini.h5`
  - `prepare_system.py validate-backbone-only` passed on both table artifacts and the prepared production `.up`
  - `prepare_system.py inject-backbone-only` and `prepare_system.py inject-backbone-cross` both executed successfully on stage `7.0`
  - `prepare_system.py handoff` was exercised repeatedly by the workflow and preserved rigid-protein checks through stages `6.0-6.6`
  - the workflow completed through stage `7.0` with only the three remaining Python files in the directory
- One stale cleanup bug appeared during the first reduced rerun:
  - `run_sim_1rkl.sh` still checked deleted `${DEPTH_TABLE_BUILDER}` / `${BACKBONE_CROSS_*}` variables in its validation preamble
  - removing that dead variable check fixed the run immediately; no further deleted-script references remained in the active shell/workflow files

## 2026-03-06 (AA-Only Workflow Cleanup Root Cause + Proof)
- The remaining MARTINI-protein leak was not in `run_sim_1rkl.sh` anymore; it was in `prepare_system_lib.convert_stage()`, which treated any existing `pdb/{pdb_id}_proa.itp` as proof of a mixed protein-lipid system even when the runtime MARTINI PDB contained only lipids and ions.
- That stale file-existence check caused stage conversion to:
  - print `=== Mixed Protein-Lipid System Detected ===`
  - load `pdb/1rkl_hybrid_proa.itp`
  - add protein topology/connectivity assumptions back into an AA-only workflow
- Fix implemented:
  - `prepare_system_lib.py` now scans the runtime PDB itself for MARTINI protein atoms and only loads protein topology/connectivity when those atoms are actually present
  - `prepare_system.py` now exports `UPSIDE_RUNTIME_ITP_FILE` only when the runtime PDB contains MARTINI protein atoms and only requires a protein ITP in that case
- Reduced end-to-end rerun (`example/16.MARTINI/outputs/test_aa_only_workflow/`) proved the corrected behavior:
  - stage converter repeatedly printed `=== Environment-Only MARTINI System Detected ===`
  - stage converter printed `Ignoring protein topology file because runtime PDB has no MARTINI protein atoms: .../pdb/1rkl_hybrid_proa.itp`
  - packed runtime MARTINI PDB `hybrid_prep/1rkl_hybrid.MARTINI.pdb` contained zero protein MARTINI atoms
  - production checkpoint `checkpoints/1rkl.stage_7.0.up` had `particle_class_counts {'ION': 125, 'OTHER': 4424, 'PROTEINAA ': 124}`, so the only protein particles present at runtime were the injected AA carriers
  - exported stage-7 VTF ended with `atom 4549 name N`, `atom 4550 name CA`, `atom 4551 name C`, `atom 4552 name O`, and continued in the same `N/CA/C/O` pattern through the final protein atoms, with no `BB` or `SC*` protein atoms present

## 2026-03-06 (PDB Placement + Recenter Root Cause)
- The user-visible “protein sticking out” in `1rkl.stage_7.0.vtf` came from two distinct frame bugs:
  - prep/stage generation bugs:
    - `prepare_system.py` / `prepare_system_lib.py` aligned the protein into the bilayer template frame instead of shifting the bilayer into the protein PDB frame
    - `prepare_system_lib.convert_stage()` then re-centered MARTINI positions again when writing `.up` stage files
    - later `inject_hybrid_mapping()` appended `PROTEINAA` carriers from the mapping/reference frame, so protein carriers and MARTINI lipids ended up in different coordinate frames even in fresh stage files
  - export bug:
    - `extract_martini_vtf.py` re-centered mode-2 output by periodic protein COM, so the VTF no longer showed the raw stage coordinates
- Fix implemented:
  - prep now shifts the bilayer template into the protein PDB frame instead of shifting the protein alone
  - stage conversion now preserves the boxed runtime PDB coordinates and no longer recenters/wraps them before writing `.up`
  - `set_initial_position.py` no longer has production BB recenter logic
  - `extract_martini_vtf.py` no longer recenters output and now writes the actual stage input frame first
- Proof from the fresh reduced run `example/16.MARTINI/outputs/test_pdb_position_preserved2/`:
  - packed-system prep summary recorded the applied bilayer alignment shift:
    - `[-25.227, -25.873, -42.056]` Å
  - stage `6.0` input geometry is now aligned:
    - protein center `[59.972, 58.817, 57.250]`
    - lipid center `[56.454, 58.318, 57.726]`
    - offset `[3.518, 0.500, -0.475]` Å
  - regenerated `1rkl.stage_7.0.vtf` now matches the raw stage file exactly:
    - frame `0` vs `stage_7.0.up` input `max_abs_diff = 5e-4 Å`
    - frame `1` vs `stage_7.0.up` first output `max_abs_diff = 5e-4 Å`
  - the previously referenced artifact `outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` was regenerated with the new exporter and its frame `0` also now matches its stage input within `5e-4 Å`
- Separate remaining blocker:
  - later stages can still generate very large MARTINI energies, and in the reduced rerun the stage `7.0` input lipid coordinates were already badly distorted by the time production started
  - that is no longer a preparation/export-frame bug; it is the separate environment-side stability issue already observed in this workflow

## 2026-03-06 (C++ MARTINI-Protein Runtime Removal + Rebuilt Proof)
- The remaining C++ MARTINI-protein presentation was split into two pieces:
  - `src/main.cpp` still exposed a MARTINI `BB` rigid-hold CLI hook/registration path
  - `src/martini.cpp` still contained the dead `martini_hybrid` runtime layer with protein-membership state, BB refresh/projection, and old hybrid coupling machinery
- Fix implemented:
  - removed the `--martini-hold-backbone` CLI path and the `register_fix_rigid_backbone_for_engine(...)` registration from `src/main.cpp`
  - removed `register_fix_rigid_backbone_for_engine(...)` and the entire dead `martini_hybrid` namespace from `src/martini.cpp`
  - simplified the active MARTINI pair kernel so it no longer consults hybrid protein state or MARTINI-protein proxy logic; the only remaining protein/environment coupling path is the explicit backbone cross node
- Post-fix source sweep:
  - `rg` over `src/` found no remaining matches for `martini_hybrid::`, `HybridRuntimeState`, `protein_membership`, `register_fix_rigid_backbone_for_engine`, or `martini-hold-backbone`
  - the only remaining protein-aware MARTINI code in `src/martini.cpp` is `martini_rbm_cross_potential`, which reads explicit `protein_atom_indices` and `protein_class_index` from injected production nodes and couples AA backbone carriers to environment atoms
- Rebuild + rerun proof:
  - `cmake --build obj -j 4` passed after the C++ cleanup
  - reduced rerun `example/16.MARTINI/outputs/test_aa_only_cpp_cleanup/` completed through stage `7.0`
  - pre-production stages `6.0-6.6` again printed `Rigid-protein check ... max_disp=0`
  - production handoff again printed `Production handoff protein-particle check: n_protein=124 max_disp=0`
  - packed runtime MARTINI PDB residue counts were exactly `{'DOPC': 4424, 'NA': 64, 'CL': 61}` with no protein residue present
  - stage `7.0` particle classes were `{'OTHER': 4424, 'ION': 125, 'PROTEINAA': 124}`
  - the tail of `1rkl.stage_7.0.vtf` remained `N/CA/C/O` only

## 2026-03-05 (Pre-Production Rigid Control Regression Root Cause)
- Before the workflow patch on 2026-03-05, generated stage files under `outputs/martini_test_1rkl_hybrid/checkpoints/` did not contain `/input/hybrid_control`; spot checks on `1rkl.stage_6.6.up` and `1rkl.stage_7.0.up` confirmed the group was absent at that point.
- The existing C++ pre-production rigid-mask logic in `../../src/martini.cpp` only activates when `/input/hybrid_control` exists with:
  - `enable=1`
  - `preprod_protein_mode=rigid`
  - `prep_runtime_mode=dry_martini_prep`
- `run_sim_1rkl.sh` had been reduced to injecting only `hybrid_bb_map` and `hybrid_env_topology`, so the workflow no longer emitted the control metadata required to freeze protein before production.
- Restoring rigidity does not require a C++ change; it is sufficient for the workflow to synthesize `/input/hybrid_control` into each prepared stage file and set production activation explicitly (`activation_stage=production`, `active_runtime_mode=aa_backbone_explicit_lipid`).

## 2026-03-06 (Why Hybrid Control Still Did Not Freeze Pre-Production Protein)
- After `/input/hybrid_control` was restored in current stage files, stage `6.6` still showed severe internal protein deformation relative to stage `6.0`; minimum-image pairwise backbone/carrier distances changed by tens of Angstroms, so this was not just translation or wrapping noise.
- The active runtime binary still contains hybrid-control strings, but the current source/runtime path no longer registers hybrid state from `../../src/main.cpp`:
  - `../../src/martini.cpp` still defines `martini_hybrid::register_hybrid_for_engine(...)`
  - `../../src/main.cpp` no longer calls it
  - therefore `/input/hybrid_control` attrs in `run_sim_1rkl.sh` stage files are inert for pre-production rigidity in the current workflow
- The active and sufficient rigid mechanism is `/input/fix_rigid`, which is still registered in `../../src/main.cpp` and is honored by both:
  - `martini_fix_rigid::apply_fix_rigid_{minimization,md}` in `../../src/martini.cpp`
  - fixed-atom skipping during semi-isotropic barostat scaling in `../../src/box.cpp`
- A short MD probe on a temp pre-production stage file with explicit `/input/fix_rigid` enabled kept all 298 protein atoms exactly stationary (`max_abs_displacement = 0.0` after 20 steps), confirming this is the correct workflow-level enforcement path.

## 2026-03-06 (Direct C++ Audit: What Rigid-Hold Code Exists)
- Yes, rigid-hold code exists in `src/*.cpp`.
- Active/runtime-wired path:
  - `../../src/main.cpp` registers `/input/fix_rigid` via `martini_fix_rigid::register_fix_rigid_for_engine(...)`
  - `../../src/main.cpp` applies `martini_fix_rigid::apply_fix_rigid_md(...)` before integration, after integration, and again after barostat updates
  - `../../src/martini.cpp` implements the `martini_fix_rigid` registry and the actual force/momentum zeroing
  - `../../src/box.cpp` skips fixed atoms during semi-isotropic scaling, so fixed atoms are not deformed by the barostat
  - `../../src/deriv_engine.cpp` builds fixed-atom masks from the same registry
- Present but currently inactive path:
  - `../../src/martini.cpp` still implements hybrid-driven pre-production rigidity through `martini_hybrid::register_hybrid_for_engine(...)` and `martini_hybrid::update_stage_for_engine(...)`
  - `../../src/main.cpp` no longer declares or calls `register_hybrid_for_engine(...)`
  - result: `/input/hybrid_control` pre-production rigid attrs are not consumed by the current runtime unless that registration hook is reconnected
- Bottom line:
  - C++ rigid-hold support is present
  - the live mechanism is `fix_rigid`
  - the hybrid-metadata rigid mechanism is dead code in the current main-program wiring

## 2026-03-06 (Workflow-Proven Pre-Production Rigidity + Production Handoff)
- `example/16.MARTINI/run_sim_1rkl.sh` now contains two runtime assertions:
  - `assert_protein_rigid_stage(...)`
    - checks that protein coordinates in `/input/pos` and the last `/output/pos` frame are identical after each stage `6.0-6.6`
  - `assert_production_handoff_protein_particles(...)`
    - checks that the original protein particles in `stage_7.0.up` match the last frame of `stage_6.6.up`
    - excludes intentionally refreshed appended AA carrier atoms from this handoff comparison
- Reduced end-to-end workflow run succeeded with the assertions enabled:
  - run directory: `example/16.MARTINI/outputs/test_fix_rigid_verify2/`
  - overrides used to keep runtime short:
    - `MIN_60_MAX_ITER=5`
    - `MIN_61_MAX_ITER=5`
    - `EQ_62_NSTEPS=20`
    - `EQ_63_NSTEPS=20`
    - `EQ_64_NSTEPS=20`
    - `EQ_65_NSTEPS=20`
    - `EQ_66_NSTEPS=20`
    - `PROD_70_NSTEPS=1`
    - `EQ_FRAME_STEPS=10`
    - `PROD_FRAME_STEPS=1`
    - `BACKBONE_CROSS_ENABLE=0`
- Observed script-side proof lines:
  - `Rigid-protein check stage 6.0 ... max_disp=0`
  - `Rigid-protein check stage 6.1 ... max_disp=0`
  - `Rigid-protein check stage 6.2 ... max_disp=0`
  - `Rigid-protein check stage 6.3 ... max_disp=0`
  - `Rigid-protein check stage 6.4 ... max_disp=0`
  - `Rigid-protein check stage 6.5 ... max_disp=0`
  - `Rigid-protein check stage 6.6 ... max_disp=0`
  - `Production handoff protein-particle check ... max_disp=0`
- Interpretation:
  - protein atoms did not move during any pre-production stage in the tested workflow run
  - the protein state handed into production matched the rigid stage-6.6 endpoint for the original protein particles

## 2026-03-06 (Fresh Confirmation on Current Script State)
- A second reduced rerun against the current `run_sim_1rkl.sh` state also passed with the same assertions:
  - run directory: `example/16.MARTINI/outputs/test_fix_rigid_verify3/`
  - same reduced settings as the previous proof run (`MIN_60_MAX_ITER=5`, `MIN_61_MAX_ITER=5`, `EQ_62-66_NSTEPS=20`, `PROD_70_NSTEPS=1`, `BACKBONE_CROSS_ENABLE=0`)
- Observed directly from the fresh rerun:
  - stages `6.0-6.6` each printed `Rigid-protein check ... max_disp=0`
  - the production transition printed `Production handoff protein-particle check: n_protein=65 max_disp=0`
  - the workflow then completed through stage `7.0`
- This removes ambiguity about stale artifacts: the currently checked-in workflow script itself now reproduces the rigid pre-production behavior on a fresh run.

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

## 2026-03-06 (Fresh Independent Verification of Current Rigid-Protein Workflow)
- A fresh reduced rerun of the current checked-in `example/16.MARTINI/run_sim_1rkl.sh` completed through stage `7.0` in:
  - `example/16.MARTINI/outputs/test_double_check_rigid_20260306/`
- Fresh script-side proof lines from that run:
  - `Rigid-protein check stage 6.0: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.1: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.2: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.3: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.4: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.5: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.6: n_protein=124 max_disp=0`
  - `Production handoff protein-particle check: n_protein=124 max_disp=0`
- Independent post-run HDF5 inspection on the generated checkpoints confirmed the same result without relying on shell logging:
  - stages `6.0-6.6` each contain `/input/fix_rigid` with `enable=1`
  - stages `6.0-6.6` each have `124` protein particles by `/input/hybrid_env_topology/protein_membership`
  - for each pre-production stage, protein `input` coordinates equal the stage last `output` coordinates exactly (`max_disp=0`)
  - stage `7.0` has `/input/fix_rigid enable=0`, so production is not accidentally still frozen
  - the original protein particles in `stage_7.0.up` exactly match the `stage_6.6.up` last-frame protein coordinates (`max_disp=0`)
- Interpretation:
  - the current workflow still enforces rigid protein hold during pre-production by the active `fix_rigid` path
  - the protein coordinates remain unchanged right up to the production handoff
  - the separate large MARTINI-energy instability remains present, but it is not a protein-drift or handoff-mutation regression

## 2026-03-06 (Refined Spline Cross-Force-Field Build)
- The dry-table descriptor `mean epsilon to all dry types` cleanly separates the five requested anchor types, unlike `epsilon(P2,type)` which leaves `Q0` and `Qa` degenerate:
  - `Q0 = 1.0526315789473684`
  - `Qa = 1.142105263157895`
  - `Na = 1.8000000000000005`
  - `C3 = 2.544736842105263`
  - `C1 = 2.744736842105263`
- Fresh refined artifact build in `example/16.MARTINI/outputs/test_refined_cross_20260306/` produced:
  - depth metadata with `aggregation_mode=empirical_depth_samples`,
  - per-bead and per-type sample summaries,
  - raw signed depth samples,
  - symmetry audit metadata,
  - spline artifact schema `martini_backbone_spline_cross_v1`.
- Symmetry audit on the current `parameters/ff_2.1/membrane.h5` shows the bilayer model is not symmetric:
  - worst raw mismatch is `cb:ARG:burial0`
  - `max_abs_diff = 5.329984533383209`
  - worst sampled depth = `13.1 A`
  - values at the worst point are `-0.40438` vs `4.92560`
- The generated spline artifact contains identical `N/CA/C/O` tables by construction and was confirmed by direct HDF5 comparison (`allclose` exact within stored values).
- Mass-normalization checks for the new convention:
  - standard dry bead `72 / 72 = 1.0`
  - ring bead `45 / 72 = 0.625`
  - injected carrier masses are `N=14/54`, `CA=12/54`, `C=12/54`, `O=16/54`, summing to `1.0`
- Migration detail discovered during validation:
  - prepared production `.up` files may already contain legacy `martini_rbm_cross_potential`,
  - when injecting the new spline node, the injector must remove the legacy node too; otherwise both cross-force terms coexist in the stage file.
- Runtime validation result:
  - a copied `1rkl.stage_7.0.prepared.up` accepted the new spline node, passed HDF5 validation, and completed a one-step `obj/upside` run,
  - the startup MARTINI energy remains extremely large, so schema/runtime integration is working but physical calibration remains poor.
- Canonical artifact publication result:
  - `parameters/ff_2.1/martini.h5` has now been replaced with the refined spline artifact
  - validated schema: `martini_backbone_spline_cross_v1`
  - direct SHA1 comparison confirms it is byte-identical to `outputs/test_refined_cross_20260306/backbone_cross_interaction_table.h5`
- Workflow-separation result:
  - `run_sim_1rkl.sh` no longer contains the `build-depth-table` or `build-backbone-cross-table` calls,
  - the simulation workflow now requires a prebuilt published artifact and points users to `example/16.MARTINI/build_martini_parameters.py` when the artifact is missing or outdated,
  - the standalone builder script successfully reproduces the spline artifact and publishes a valid `martini.h5` outside the simulation workflow.

## 2026-03-06 (User Feedback: Backbone `O` Is the Worst Case)
- User-reported behavior: both the current spline-table version and the previous version perform especially badly for backbone `O` among `N/CA/C/O`.
- This is consistent with the present builder design:
  - `example/16.MARTINI/prepare_system.py` currently writes identical environment tables for all four backbone classes,
  - the current anchor target is a single scalar derived from the averaged hydrogen-bond channel (`mean_hb_attraction`), not a role-specific donor/acceptor decomposition.
- Technical implication:
  - the current “all backbone classes identical” simplification is likely too aggressive,
  - `O` should be treated separately in the next refinement, and the next model should probably use acceptor-specific membrane information for `O` instead of the shared averaged target used now.

## 2026-03-06 (Role-Specific Backbone Calibration Applied)
- The published builder now uses role-specific anchor targets:
  - `N -> mean_hb_donor_attraction`
  - `O -> mean_hb_acceptor_attraction`
  - `CA/C -> mean_cb_value`
- For the `Q0` anchor in the rebuilt metadata:
  - `N` target field `mean_hb_donor_attraction = 0.16819868543005373`
  - `CA` target field `mean_cb_value = 0.8677337670194787`
  - `C` target field `mean_cb_value = 0.8677337670194787`
  - `O` target field `mean_hb_acceptor_attraction = 0.44794905782531336`
- Direct HDF5 checks on the rebuilt/published artifact show the role split took effect:
  - schema remains `martini_backbone_spline_cross_v1`
  - `max |N-O|` for the `Q0` spline table is `2.141817092895508`
  - `max |CA-C|` for the `Q0` spline table is `0.0`
  - `Q0` minima changed from shared to role-specific:
    - `N = -1.287758231163025`
    - `O = -3.4295754432678223`
- The rebuilt artifact was published to `parameters/ff_2.1/martini.h5` and confirmed byte-identical to `example/16.MARTINI/outputs/test_o_role_refine_20260306/backbone_cross_interaction_table.h5` (`SHA1 8c9f0542661c9ca7104e4d9423f30185d160e07b`).
- Runtime integration remained valid:
  - injected copied production stage passed `validate-backbone-only`
  - one-step `obj/upside` startup/run completed
  - the very large MARTINI energy remained essentially unchanged, so the change corrected class differentiation but not the main physical instability.

## 2026-03-06 (Production Rg Regression: Legacy Selector Was Reading Environment Beads)
- The active runtime `Rg` path in `../../src/main.cpp` had regressed to a legacy selector:
  - `collect_rg_backbone_indices_from_sequence()` returned the first `3*n_res` particle indices
  - `collect_rg_backbone_indices()` used that selector directly with no `hybrid_bb_map` override
- In the current AA-only MARTINI workflow this is invalid because protein atoms are appended after the environment:
  - in `outputs/test_pbc_wrapped/checkpoints/1rkl.stage_7.0.up`, `n_res = 31`, so the legacy selector used indices `0..92`
  - those indices have `protein_membership < 0` for every selected particle, so they are all non-protein environment beads
  - the explicit mapped protein-backbone atoms live in `/input/hybrid_bb_map/{atom_indices,atom_mask}` with `124` valid indices spanning `4549..4672`
- Direct checkpoint measurement confirms the reporting error:
  - same `stage_7.0.up` input coordinates give `Rg = 57.26 A` with the legacy `0..92` selector
  - the mapped protein-backbone selector gives `Rg = 12.70 A`
- Conclusion:
  - the inflated first-step production `Rg` was a reporting-selection regression, not evidence that the rigid protein geometry had changed
