# Task Plan

## 2026-03-06 Reset: Restore Production Rg to Protein Backbone Only

### Project Goal (Current)
- Verify how production-stage `Rg` is calculated in the active `example/16.MARTINI/run_sim_1rkl.sh` workflow.
- Fix the runtime status-line `Rg` so it reflects the protein backbone only in the AA-only MARTINI workflow.
- Preserve the existing stage dynamics and PBC/output behavior; only correct the reporting path.

### Architecture & Key Decisions (Current)
- Treat this as a runtime reporting regression first, not as a physical deformation regression, because the user reports the protein remains rigid while `Rg` changes sharply.
- Validate the active runtime selector against actual `stage_7.0` checkpoint data before changing the code.
- Prefer the workflow's explicit backbone mapping (`/input/hybrid_bb_map`) over positional assumptions such as "first `3*n_res` atoms".
- Keep a conservative fallback for non-hybrid inputs when the backbone mapping is absent.

### Execution Phases (Current Run)
- [x] Phase A: Trace the current runtime `Rg` path and compare it against generated production-stage checkpoint data.
- [x] Phase B: Implement the minimal fix so runtime `Rg` uses mapped protein-backbone particles.
- [x] Phase C: Rebuild and verify the first production-stage `Rg` printout on a reduced run.

### Known Errors / Blockers (Current Run)
- None yet.

### Review (Current Run)
- Root cause was a runtime reporting regression in `src/main.cpp`: `collect_rg_backbone_indices()` had fallen back to legacy `3*n_res` indexing.
- In the active AA-only MARTINI workflow that selector pointed at particles `0..92`, which are environment beads, not protein atoms.
- The generated production checkpoint `example/16.MARTINI/outputs/test_pbc_wrapped/checkpoints/1rkl.stage_7.0.up` contains the explicit protein-backbone mapping under `/input/hybrid_bb_map`:
  - valid mapped backbone atoms: `124`
  - mapped index range: `4549..4672`
- `src/main.cpp` now prefers `/input/hybrid_bb_map/{atom_indices,atom_mask}` for `Rg` selection and falls back to the legacy sequence-based selector only when that mapping is absent.
- Verification passed after rebuild:
  - direct HDF5 check on `stage_7.0.up`: legacy selector `Rg = 57.26 A`, mapped backbone selector `Rg = 12.70 A`
  - one-step runtime probe on a copied production checkpoint printed `Rg  12.7 A` at step `0`

## 2026-03-06 Reset: Keep MARTINI Workflow Coordinates Inside Primary PBC Box

### Project Goal (Current)
- Verify whether the current `example/16.MARTINI/run_sim_1rkl.sh` workflow keeps all stage coordinates inside the active periodic box.
- Fix the narrowest layer that lets particles appear outside the primary box during preparation, stage handoff, simulation output, or export.
- Preserve the existing physical semantics: no protein-only recentering, no architecture drift, and no changes beyond what is needed for PBC containment.

### Architecture & Key Decisions (Current)
- Treat this as a coordinate-frame/boxing bug first, not a force-field bug.
- Validate on fresh artifacts from the current workflow code, not older mixed-protein outputs that predate the AA-only cleanup.
- Prefer a single canonical wrapping policy at workflow boundaries:
  - prepared stage inputs should start inside the primary box
  - stage handoff should keep positions in the current box
  - exported trajectories should not make in-box coordinates look out-of-box
- If runtime C++ output still leaks out of `[0, L)`, fix the minimal live boundary where wrapped coordinates are expected instead of reintroducing protein-centric recentering.

### Execution Phases (Current Run)
- [x] Phase A: Generate fresh current-code artifacts and identify the first boundary where coordinates leave `[0, L)`.
- [x] Phase B: Implement the minimal PBC containment fix in the responsible workflow/runtime layer.
- [x] Phase C: Re-run a reduced workflow and prove all particles remain inside the primary box for the checked stages/exports.

### Known Errors / Blockers (Current Run)
- Older `example/16.MARTINI/outputs/test_fix_rigid_verify*` artifacts still reflect pre-AA-only workflow states and cannot be used as the final correctness signal for this task.
- The separate late-stage MARTINI instability is still present in stages `6.4+` and remains out of scope for this PBC containment fix.

### Review (Current Run)
- Added a shared coordinate wrapper in `lib.py` and applied it when stage conversion writes boxed coordinates into `.up` files.
- Updated `prepare_system.py` handoff logic to wrap the carried last-frame coordinates after any backbone-carrier refresh and before writing the next stage input.
- Updated `src/main.cpp` position logging so sampled `output/pos` frames are wrapped against the live NPT box or the static MARTINI box before being written to HDF5.
- Reconfigured the local build directory with `cmake -S src -B obj` because the existing generated link rules referenced a missing HDF5 soname (`libhdf5.320.0.0.dylib`), then rebuilt successfully with `cmake --build obj -j 4`.
- Fresh reduced rerun in `example/16.MARTINI/outputs/test_pbc_wrapped/` confirmed zero out-of-box atoms in both `/input/pos` and the last sampled `/output/pos` for stages `6.0`, `6.2`, `6.4`, `6.6`, and `7.0`.
- Export verification also passed: `example/16.MARTINI/outputs/test_pbc_wrapped/1rkl.stage_7.0.vtf` had `0` out-of-box atoms in both frames.

## 2026-03-06 Reset: Three-Script Workflow Cleanup

### Project Goal (Current)
- Reduce `example/16.MARTINI` to exactly three Python files:
  - `prepare_system.py`
  - `extract_martini_vtf.py`
  - `lib.py`
- Keep the existing workflow behavior, including the recent PDB-placement and AA-only protein fixes.
- Remove direct shell/script dependencies on helper executables such as `set_initial_position.py`, `inject_backbone_only_nodes.py`, and the table-builder/validator scripts.

### Architecture & Key Decisions (Current)
- Rename `prepare_system_lib.py` to `lib.py` and keep the existing prep/stage-conversion logic there.
- Move the active helper script entrypoints into `prepare_system.py` as subcommands, so the workflow has one operational Python command surface plus the VTF exporter.
- Update `run_sim_1rkl.sh` and `run_sim_bilayer.sh` to call only `prepare_system.py` and `extract_martini_vtf.py`.
- Delete obsolete Python scripts after the shell workflows no longer reference them.
- Validate with syntax checks plus a workflow smoke test on the updated script surface, not just import checks.

### Execution Phases (Current Run)
- [x] Phase A: Rename the helper library to `lib.py` and add the required helper subcommands to `prepare_system.py`.
- [x] Phase B: Update shell workflows to use only `prepare_system.py` and `extract_martini_vtf.py`.
- [x] Phase C: Remove the extra Python files under `example/16.MARTINI`.
- [x] Phase D: Run syntax/workflow validation and record the result.

### Known Errors / Blockers (Current Run)
- The pre-existing large MARTINI-energy instability in later stages remains out of scope for this cleanup.

### Review (Current Run)
- `example/16.MARTINI` now contains exactly three Python files: `prepare_system.py`, `extract_martini_vtf.py`, and `lib.py`.
- `prepare_system_lib.py` was renamed to `lib.py`, and the active helper entrypoints were folded into `prepare_system.py` as subcommands:
  - `handoff`
  - `inject-backbone-only`
  - `validate-backbone-only`
  - `build-depth-table`
  - `build-backbone-cross-table`
  - `inject-backbone-cross`
- `run_sim_1rkl.sh` and `run_sim_bilayer.sh` now call only `prepare_system.py` and `extract_martini_vtf.py`.
- The obsolete helper scripts were deleted from `example/16.MARTINI`.
- Validation passed:
  - `python3 -m py_compile example/16.MARTINI/prepare_system.py example/16.MARTINI/lib.py example/16.MARTINI/extract_martini_vtf.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh example/16.MARTINI/run_sim_bilayer.sh`
  - reduced end-to-end rerun in `outputs/test_three_script_cleanup/` completed through stage `7.0` using the new subcommands, including table build, validation, backbone-node injection, handoff, cross-node injection, and VTF export
  - the reduced rerun still showed the pre-existing large MARTINI energies in late pre-production, but the three-script workflow itself completed successfully

## 2026-03-06 Reset: Preserve PDB Bilayer Placement

### Project Goal (Current)
- Preserve the protein-bilayer placement implied by the input protein PDB when generating the initial hybrid system.
- If the assembled system is translated to center it in the box, apply that translation to protein and bilayer together as one rigid body.
- Remove production-stage recentering from the example Python workflow/export path so stage `7.0` VTF/PDB output shows the actual initial structure.

### Architecture & Key Decisions (Current)
- The prep path must no longer translate the protein to the bilayer template frame by itself.
- Instead, the bilayer template is aligned into the protein PDB frame before tiling/cropping, and later box placement shifts the assembled system rigidly.
- `set_initial_position.py` must not perform any production BB recentering; handoff is strict copy plus optional carrier refresh only.
- `extract_martini_vtf.py` must not apply protein-centric visualization recentering or box-centering for production output; exported coordinates should match the stage file.
- The example runtime continues to invoke `obj/upside` with `--disable-recentering`, so no workflow-stage C++ recentering is active.

### Execution Phases (Current Run)
- [x] Phase A: Remove prep-time protein-only translation and preserve the PDB placement by moving the bilayer frame instead.
- [x] Phase B: Remove production-stage Python recentering from stage handoff and VTF export.
- [x] Phase C: Run targeted validation and compare source PDB placement against packed/stage-7/exported coordinates.

### Known Errors / Blockers (Current Run)
- Initial placement/export bug is fixed, but later pre-production environment coordinates can still become badly distorted because of the separate large MARTINI-energy instability in stages `6.4-7.0`.

### Review (Current Run)
- `prepare_system.py` and `prepare_system_lib.py` no longer translate the protein into the bilayer template frame by itself; they align the bilayer template into the protein PDB frame and then let later box placement move the assembled system rigidly.
- `prepare_system_lib.convert_stage()` no longer recenters/wraps MARTINI coordinates before writing `.up` stage files, so MARTINI particles and injected `PROTEINAA` carriers stay in the same coordinate frame.
- `set_initial_position.py` no longer contains production BB recentering; stage handoff is now strict copy plus optional carrier refresh only.
- `extract_martini_vtf.py` no longer applies protein-centric visualization recentering and now prepends the actual stage input frame before output frames.
- Validation on `outputs/test_pdb_position_preserved2/` shows the initial packed/stage geometry is now aligned:
  - stage `6.0` input protein-lipid center offset is only `[3.518, 0.500, -0.475]` Å
  - regenerated `stage_7.0.vtf` frame `0` matches `stage_7.0.up` input coordinates within `5e-4` Å
  - the previously referenced `outputs/martini_test_1rkl_hybrid/1rkl.stage_7.0.vtf` was regenerated and its frame `0` also now matches the stage input within `5e-4` Å

## 2026-03-06 Reset: AA-Only Protein Workflow Cleanup

### Project Goal (Current)
- Remove all workflow logic that converts protein structure into a MARTINI representation.
- Keep the protein in pure Upside all-atom backbone form only (`N`, `CA`, `C`, `O`) through preparation, stage conversion, runtime, and export.
- Keep MARTINI representation only for the environment (lipids, ions, solvent/beads as applicable).
- Remove C++ runtime codepaths that still model or special-case protein as MARTINI particles.

### Architecture & Key Decisions (Current)
- The workflow must not call martinize, resolve protein MARTINI ITPs, or infer protein topology from MARTINI protein assets.
- Hybrid mapping remains only as AA-backbone carrier/reference metadata plus environment topology metadata.
- Stage conversion must treat the runtime MARTINI PDB as environment-only when no protein MARTINI topology is explicitly provided.
- Production export must use the carried Upside backbone coordinates directly and must not reconstruct protein atoms from MARTINI `BB` beads.
- The active C++ runtime must not register MARTINI-protein backbone hold hooks, hybrid protein-membership state, or MARTINI intra-protein filtering logic.
- Protein/environment coupling must come only from explicit AA-backbone production nodes and not from MARTINI-protein proxy logic.

### Execution Phases (Current Run)
- [x] Phase A: Remove any residual protein-MARTINI auto-detection from stage conversion and production-node injection.
- [x] Phase B: Verify workflow preparation emits environment-only runtime MARTINI assets plus AA-only protein reference/carrier metadata.
- [x] Phase C: Run a reduced end-to-end workflow and confirm no MARTINI protein topology/conversion path is used.
- [x] Phase D: Verify stage `7.0` exported protein contains only `N/CA/C/O` coordinates sourced from Upside carriers.
- [x] Phase E: Remove remaining C++ MARTINI-protein runtime hooks and dead hybrid registration paths.
- [x] Phase F: Rebuild, rerun a reduced end-to-end workflow, and record results in progress/findings logs.

### Known Errors / Blockers (Current Run)
- Reduced validation run still shows extremely large MARTINI energies later in equilibration/production. This task only removed protein MARTINI conversion paths and verified AA-only protein representation; it did not resolve the separate stability/energy issue.

### Review (Current Run)
- `prepare_system_lib.convert_stage()` now detects MARTINI protein presence from the runtime PDB instead of blindly keying off `*_proa.itp` file existence.
- `prepare_system.py` now exports `UPSIDE_RUNTIME_ITP_FILE` only when the runtime PDB actually contains MARTINI protein atoms and only requires a protein ITP in that case.
- Reduced workflow rerun `outputs/test_aa_only_workflow/` completed through stage `7.0` without using MARTINI protein topology:
  - converter printed `Environment-Only MARTINI System Detected`
  - converter printed `Ignoring protein topology file because runtime PDB has no MARTINI protein atoms`
  - packed runtime MARTINI PDB contained zero protein atoms
  - stage `7.0` runtime state contained `124` `PROTEINAA` particles and no MARTINI protein particles
  - stage `7.0` VTF export ended with only `N/CA/C/O` protein atoms, with no `BB` or `SC*` protein proxies
- C++ runtime cleanup is now also complete:
  - `src/main.cpp` no longer exposes or registers the MARTINI `BB` rigid-hold hook
  - `src/martini.cpp` no longer contains the dead `martini_hybrid` runtime layer or MARTINI-protein BB refresh/projection logic
  - rebuilt binary `obj/upside` passed a fresh reduced rerun in `outputs/test_aa_only_cpp_cleanup/`
  - the packed runtime MARTINI PDB in that rerun contained only environment residues (`DOPC`, `NA`, `CL`)
  - stage `7.0` runtime state again contained only `PROTEINAA` protein particles, and the exported VTF again ended with `N/CA/C/O` atoms only

## 2026-03-05 Reset: Backbone-Only + Depth Table

### Project Goal (Current)
- Remove legacy dry-MARTINI/Upside hybrid physical-interaction execution paths.
- Keep simulation target as dry-MARTINI explicit bilayer plus Upside all-atom backbone only (`N,CA,C,O`).
- Remove sidechain back-mapping/coupling from active production workflow.
- Build a depth-dependent interaction-table artifact from membrane potential data (no PO4 anchoring).

### Architecture & Key Decisions (Current)
- `cb_energy` in `parameters/ff_2.1/membrane.h5` is used as the backbone channel.
- `hb_energy` in `parameters/ff_2.1/membrane.h5` is used as the hydrogen-bond channel.
- Sidechain term is deferred (no dedicated sidechain membrane channel currently in `membrane.h5`).
- Sheet term is omitted for now (no explicit sheet channel in `membrane.h5`).
- Depth coordinate is `|z - z_center|`, where `z_center` is bilayer center from lipid geometry; PO4 is not used as an anchor.
- This run produces and validates interaction-table artifacts only; no runtime force-law injection of the new table yet.

### Execution Phases (Current Run)
- [x] Phase A: Disable core runtime hybrid hooks (`main.cpp`, `main.h`, `deriv_engine.cpp`) so legacy hybrid interaction paths are inactive.
- [x] Phase B: Add backbone-only production node injection helper (no rotamer/sidechain mapping nodes).
- [x] Phase C: Update example production workflow to use backbone-only injection and remove sidechain back-mapping activation.
- [x] Phase D: Add depth interaction-table builder for DOPC bead types from bilayer geometry + membrane potential channels (`cb`, `hb`).
- [x] Phase E: Add validation checks for backbone-only stage files and interaction-table outputs.
- [x] Phase F: Run targeted validation (syntax/build/smoke checks) and document results in progress/findings logs.
- [x] Phase G: Remove remaining legacy hybrid workflow hooks from shell/Python orchestration and keep backbone-only mapping schema (`hybrid_bb_map` + `hybrid_env_topology`) only.

### Review (2026-03-05 Current Run)
- Runtime hybrid registration/alignment hooks were removed from `src/main.h`, `src/main.cpp`, and `src/deriv_engine.cpp`.
- Backbone-only stage injector was implemented in `inject_backbone_only_nodes.py` and now removes sidechain potential nodes.
- Production workflow scripts (`run_sim_1rkl.sh`, `run_sim_1rkl_new.sh`) now call backbone-only injection; rotamer/placement augmentation blocks were removed.
- Rigid/relax/prod test workflow scripts were cleaned to remove legacy activation/sc-control hooks and legacy handoff env toggles.
- `inject_hybrid_mapping()` in active workflows now copies only backbone/environment mapping groups (`hybrid_bb_map`, `hybrid_env_topology`).
- Mapping export/validation utilities were updated to backbone-only schema expectations (no workflow dependence on legacy control/SC mapping groups).
- Depth table artifacts were generated successfully:
  - `outputs/depth_interaction_table.csv`
  - `outputs/depth_interaction_table.meta.json`
- `membrane.h5` term check confirms:
  - backbone (`cb_energy`) present
  - hbond (`hb_energy`) present
  - sidechain term absent
  - sheet term absent

## Project Goal
- Define and implement a hybrid dry-MARTINI + Upside workflow where protein dynamics are handled by Upside, environment is handled by dry MARTINI, and coupling starts only at production stage.
- Add a dedicated rigid dry-MARTINI workflow variant where the protein stays fixed for all stages (including stage 7.0), production-stage Upside coupling is never activated, and trajectory export avoids SC/backbone back-mapping.

## Architecture & Key Decisions
- Data preparation includes packing OPM protein `pdb/1rkl.pdb` into MARTINI DOPC bilayer `pdb/bilayer.MARTINI.pdb`, with MARTINI 2.2 coarse-graining and 0.15 M NaCl ion placement.
- Protein backbone state is represented by Upside atoms (`N, CA, C, O`) as primary DOFs.
- MARTINI protein backbone beads (`BB`) are computed from COM of mapped `N, CA, C, O` atoms each cycle; no per-step rebuilding of atoms from BB.
- Dry-MARTINI forces on protein `BB`/mapped sidechain proxies are projected onto protein atoms/DOFs in Upside; Upside then applies its own force field update.
- Hybrid force exchange is disabled before production stage; pre-production keeps protein rigid to relax environment.
- In production-stage hybrid mode, MARTINI intra-protein interactions are selectively filtered to avoid double-counting protein internal physics already modeled by Upside: no BB-BB, no SC-SC, and BB-SC only for same residue.
- Rotamer mapping is used for protein-environment coupling and for allowed same-residue protein `SC-BB` pairs in production so sidechain proxy forces are probability-weighted by live rotamer marginals.
- Hybrid `.up` schema will add dedicated groups under `/input`: `hybrid_control`, `hybrid_bb_map`, `hybrid_sc_map`, and `hybrid_env_topology`.
- Add optional rigid-body drift control for Upside protein motion in production: RMSD alignment or equivalent COM+rotation removal mode, enabled by configuration and evaluated for physical side effects.

## Execution Phases
- [x] Phase 1: Finalize hybrid data schema for `.up` and packed-system metadata.
- [x] Phase 2: Implement data-preparation script (CG conversion, bilayer packing, ions, mapping export).
- [x] Phase 3: Implement Upside input parsing for hybrid mapping blocks and stage-gated runtime flags.
- [x] Phase 4: Implement hybrid coupling loop in Upside (MARTINI->Upside force projection, BB COM refresh, production-only activation).
- [x] Phase 5: Implement sidechain rotamer-environment mapping path with intra-protein MARTINI exclusions.
- [x] Phase 6: Implement/benchmark rigid-body drift mitigation mode (frame-to-frame alignment option) and diagnostics.
- [x] Phase 7: Validate correctness, stability, and energy/force accounting; document usage and limitations.
- [x] Phase 8: Audit pre-production rigidity and potential-output semantics in `run_sim_1rkl.sh` workflow artifacts.
- [x] Phase 9: Implement split potential logging (`Upside-only` vs `MARTINI`) and enforce pre-production NPT coupling to match CHARMM-GUI stage settings.
- [x] Phase 10: Add and validate a local `.up` input-visualization script (`visualize_up_input.py`) with direct executable invocation and sane defaults.
- [x] Phase 11: Remove production-stage MARTINI fallback paths in hybrid SC force projection and fix stage handoff/production stability regressions.
- [x] Phase 12: Add a separate rigid-protein workflow script derived from `run_sim_1rkl.sh` that disables production-stage Upside/SC paths and keeps dry-MARTINI rigid protein interactions through stage 7.0.

## Known Errors / Blockers
- Requires concrete integration points in Upside C++ step loop and force modules.
- Rigid-body alignment scope confirmed: apply alignment to coupling coordinates only; preserve raw integrated coordinates in saved trajectories.
- Phase 4 assumes valid BB mapping indices into current position array (`bb_atom_index` preferred; residue/name inference fallback may fail on nonstandard topologies).
- Phase 5 probabilistic sidechain coupling is implemented in C++; remaining work is physical calibration/stability validation for production trajectories.
- Phase 6 alignment uses BB-frame rigid transform between consecutive steps (not full Kabsch over all protein atoms); validate stability impact in Phase 7.
- Phase 7 completed smoke/integrity validation; full long-horizon physical validation remains to be benchmarked on production trajectories.
- Hybrid mapping export now stores BB component targets in protein-AA PDB index space; stage-file injection must convert this reference index space into runtime coordinate indices before simulation.
- Legacy static protein CG assets (`pdb/1rkl.MARTINI.pdb`, `pdb/1rkl_proa.itp`) are not reliable for dry-MARTINI compatibility; workflow mitigation is runtime MARTINI2 generation via `martinize.py -ff martini22`.
- Stage 6.0 HDF5 close-ID runtime failure is resolved (hybrid-control parsing bug in C++ fixed); remaining workflow issue is downstream VTF extraction mismatch when using `extract_martini_vtf.py` with a non-runtime structure source.
- Current `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up` has no `rotamer` or `placement*_point_vector_only` nodes, so live probabilistic rotamer weighting is inactive there and hybrid sidechain rows remain static (`rotamer_id=0` fallback).
- After removing SC placement/probability fallback, production files lacking `rotamer` + `placement*_point_vector_only` nodes will fail fast in hybrid mode and must be regenerated with those nodes present.
- Current stage-generation workflow (`prepare_martini.py` in `run_sim_1rkl.sh` / `test_prod_run_sim_1rkl.sh`) does not emit Upside `affine_alignment`/`placement_*`/`rotamer` nodes, so production-stage strict SC coupling fails unless scripts augment stage-7 files with these nodes.
- Stage-7 script augmentation now injects required `affine_alignment`/`placement_fixed_*`/`rotamer` nodes with `ff_2.1/sidechain.h5` compatibility; previous production SIGBUS (degenerate affine triplets) is fixed. Remaining blocker is production-stage physical instability (rapid potential blow-up/NaN) rather than a runtime memory crash.
- Fixed 2026-02-18: `../../src/martini.cpp` SC rigid mapping now restores target centroid when mapping centroid-relative `local` coordinates (`mapped = R*local + tgt_centroid`), removing the previously confirmed translation error.
- Fixed 2026-02-18 follow-up: SC row-expansion now caps per-residue base proxy rows by placement-point capacity per rotamer (minimum placement group size), preventing over-expanded proxy rows from residues whose MARTINI proxy count exceeds placement representation.
- Remaining blocker (post-fixes): long-horizon production stability/physical calibration still needs benchmarking; 20-step production probes are now finite and avoid previous `1e7-1e12` spike regimes, but extended trajectories still require validation.
- Fixed 2026-02-19: production SC probabilistic coupling no longer uses linear `sum_r p_r E_r` aggregation; runtime now applies Boltzmann-consistent `-log sum_r p_r exp(-E_r)` mixing (with stabilized numerics), removing reproduced `1e5-1e7` short-horizon MARTINI potential spikes from low-probability clash rotamers.
- Fixed 2026-02-20: SC inner relaxation loop now integrates with force direction (`-dE/dx`) instead of gradient-ascent (`+dE/dx`), addressing reproduced production energy blow-up from unstable SC relaxation dynamics.
- Fixed 2026-02-20: stage handoff (`set_initial_position.py`) now recenters production-stage protein BB to box center using periodic circular-center estimation and global coordinate shift+wrap, preventing protein corner placement/PBC split at stage 7 start.
- Fixed 2026-02-20: production mode-2 trajectory export (`extract_martini_vtf.py`) now recenters by periodic protein COM (residue `PRO`) before wrapping, preventing visualization-only PBC splitting that made centered proteins appear at box corners.
- `output/potential` is a total-engine scalar (sum over active potential nodes), not a protein-only all-atom Upside-backbone channel; this can vary during pre-production even with rigid protein because environment terms evolve.
- Remaining blocker: full long-horizon validation with exact production SC settings (`sc_env_relax_steps=200`) is still running/needs benchmarking; short diagnostics are improved but not yet a complete 5k-step confirmation.
- New blocker (2026-02-24): production-stage instability is now isolated to the hybrid-active BB coupling path in `../../src/martini.cpp` (`refresh_bb_positions_if_active` + `project_bb_gradient_if_active`), not the non-protein hard-sphere toggle or SC environment loop. Hybrid-active runs with BB mapping enabled still show runaway expansion, while equivalent runs with active stage disabled or BB map masks zeroed remain stable.

## Revised Decisions
- Hybrid coupling starts only at production stage; pre-production protein remains rigid.
- Pre-production rigid hold keeps the entire protein fully fixed in space (`x,y,z`) while lipid `PO4` headgroup beads are constrained only along `z` (`x,y` remain free); other lipid beads remain unconstrained.
- MARTINI intra-protein interactions in production hybrid mode follow selective rules:
  - no BB-BB (bonded or nonbonded),
  - no SC-SC (bonded or nonbonded),
  - BB-SC allowed only for same residue.
- BB MARTINI positions are derived each cycle from COM of Upside `N,CA,C,O`; no inverse backbone reconstruction from BB.
- Deterministic single-rotamer sidechain projection in C++ is removed; sidechain coupling uses per-rotamer probabilistic mapping rows from `hybrid_sc_map`.
- Sidechain probabilistic projection/force evaluation is implemented in Upside C++ runtime and is active only during production-stage hybrid coupling; preparation scripts remain unchanged for pre-production rigid setup.
- If enabled, frame-to-frame RMSD alignment is used only in coupling frame construction and does not overwrite primary trajectory coordinates.
- Packed-system hybrid mapping artifacts must keep all BB/SC projection indices in the same runtime coordinate space as `protein_membership`; validator/runtime now enforce protein-local bounds.
- Hybrid sidechain placement in production will be sourced from Upside rotamer placement coordinates (`placement_fixed_point_vector_only` / `placement_point_vector_only`) per rotamer state when available, with runtime fallback to `hybrid_sc_map` projected rows if placement nodes are absent.
- Probabilistic sidechain weights in production will use live Upside rotamer marginals from C++ runtime rotamer state when available, with fallback to static `hybrid_sc_map/rotamer_probability` only when rotamer marginals are unavailable.
- Hybrid runtime may add ordering dependencies so rotamer/placement sidechain states are computed before MARTINI coupling each step; this is a scheduling change only (no Python preparation changes).
- Sidechain-environment force transfer in hybrid production prefers explicit `hybrid_sc_map` projection targets/weights generated from martinize bonded topology, with residue-based MARTINI `BB` target lookup retained only as fallback.
- Production-stage MARTINI protein interaction filtering is runtime-enforced in C++ (based on `protein_membership`, `atom_roles`, and `residue_ids`) so non-production stages remain unaffected while protein is rigid.
- Production-stage MARTINI nonbonded interactions for non-protein particles are hard-sphere-like (repulsive-only) while protein-involving MARTINI interactions remain dry-MARTINI for hybrid coupling.
- Production probabilistic sidechain coupling includes allowed same-residue protein `SC-BB` interactions (not only protein-environment), with per-rotamer force/potential weighting from live `rotamer` marginals.
- Production probabilistic sidechain MARTINI coupling (SC-environment and allowed same-residue SC-BB) must use Boltzmann-consistent free-energy mixing across rotamer states (`-log sum_r p_r exp(-E_r)`) instead of linear energy averaging (`sum_r p_r E_r`) to avoid instability from low-probability clash rotamers.
- Preparation exports per-`BB` backbone reference metadata (`N,CA,C,O` atom indices/coordinates plus comment strings) in `hybrid_bb_map` using protein-AA PDB index space; stage-file injection maps these references into runtime coordinate carriers for active BB COM refresh and BB-force projection.
- Runtime hybrid scheduling enforces `martini_potential` execution after all coordinate-node updates each integration sub-step, so MARTINI coupling always consumes the latest backbone-derived coordinates for that sub-step.
- `run_sim_1rkl.sh` defaults to MARTINI2.2 protein coarse-graining via local `martinize.py` (`-ff martini22`) and performs a preflight protein-ITP-vs-dry-FF mass compatibility check to prevent silent MARTINI3/dry-MARTINI mixing.
- Hybrid rotamer/placement node binding is discovered from the built `DerivEngine` graph at registration time (not by pre-scanning HDF5 `/input/potential`), to avoid HDF5 handle/path edge cases in hybrid-control parsing.
- Bilayer preparation for hybrid packing must tile/crop lipid coordinates in XY so lipid edges extend to at least one protein-length margin beyond protein edges on each side, and the XY box edges are set exactly to lipid edges.
- Ion placement for target salinity uses a template-calibrated effective ion-accessible volume fraction (derived from bilayer template ions/box) to avoid full-box overestimation in dry-MARTINI packed systems.
- Hybrid bilayer packing enforces square XY dimensions in the generated system by expanding the shorter side to the longer side after minimum protein-margin constraints are applied.
- Hybrid packed PDB output canonicalizes lipid residue naming to `DOPC` (accepting `DOP`/`DOPC` as inputs) and preserves 4-character residue-name formatting.
- Hybrid box sizing enforces a minimum `Z` length of `3x` the protein `Z` span in generated packed systems.
- Hybrid packing places the lipid bilayer midplane at the box `Z` center in generated PDBs.
- Pre-production rigid hold (protein + configured lipid headgroup roles, default `PO4`) is fixed in space, not only force-constrained.
- Pre-production fixed-in-space enforcement uses an integration-level fixed-atom mask (skip momentum/position updates for fixed atoms) plus force/momentum zeroing; direct per-step coordinate rewrites are avoided to prevent early-step NaN instability.
- Lipid-overlap removal for protein insertion uses a tighter default cutoff (`3.0 Å`) to avoid creating an artificially large annular cavity around the protein.
- NPT barostat updates in hybrid pre-production rigid mode remain enabled; scaling is applied only to unconstrained DOFs (protein fixed atoms are not scaled, z-fixed headgroups remain z-fixed).
- `run_sim_1rkl.sh` converts requested MD steps to Upside CLI time units (`duration`, `frame-interval`) to avoid unintended 100x step inflation.
- VTF export is stage-wise: modes are `mode 1` (all MARTINI particles) for non-production stages and `mode 2` (non-protein MARTINI + protein all-atom backbone backmapped from `hybrid_bb_map`) for production stage.
- Upside CLI supports `--duration-steps` for step-count-driven stage lengths, and `run_sim_1rkl.sh` uses it so stage step counts are independent of `dt`.
- Sidechain-to-backbone force transfer mapping will prefer explicit `hybrid_sc_map` projection targets/weights (prepared from martinize.py protein bonded topology/force constants) with BB-residue lookup retained only as runtime fallback.
- Minimization mode no longer forces a post-minimization stage switch to `production`; runtime now restores the pre-minimization stage label to prevent production-only hybrid/UpSide activation in pre-production stage files.
- Production-stage hard-sphere behavior is implemented in C++ MARTINI runtime as repulsive-only WCA-like nonbonded interactions for non-protein/non-protein pairs; this is controlled by `hybrid_control` and active only when hybrid is active.
- Sidechain placement/probability fallback is removed in production hybrid mode: SC mapping requires live Upside `rotamer` and `placement*_point_vector_only` nodes, and static `hybrid_sc_map` fallback is no longer used for SC position/probability evaluation.
- SC rotamer projection in production uses a residue-level rigid transform between reference and current rotamer placement frames so mapped SC geometry remains rigid across rotamer states.
- Production-stage workflow scripts will explicitly augment stage-7 `.up` files with `affine_alignment`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, and `rotamer` nodes (built from sidechain library + current MARTINI protein mapping) before running production.
- Pre-production rigid hold verification should use runtime hybrid metadata (`hybrid_control` + `protein_membership`) and trajectory displacement checks; presence of `/input/fix_rigid` datasets is not required because pre-production fixed masks are generated dynamically in C++ runtime.
- Potential logging will be split:
  - `output/potential`: non-MARTINI (Upside-side) potential contribution only.
  - `output/martini_potential`: MARTINI potential contribution, preserving visibility of the previously observed scalar behavior for MARTINI stages.
- Pre-production NPT coupling should follow CHARMM-GUI stage intent:
  - stages `6.0-6.1`: `ref_p = 0.0 bar`, `tau_p = 4.0`, `compressibility = 3e-4 bar^-1 (xy), 0.0 (z)`.
  - stages `6.2-6.6`: `ref_p = 1.0 bar` (semi-isotropic; `z` compressibility remains `0.0`), `tau_p = 4.0`, same compressibility.
  - convert pressure to Upside units using `1 bar = 0.000020659477 E_up/Angstrom^3`.
- Hybrid pre-production rigid hold and NPT are both enabled by scaling only unconstrained DOFs during barostat box/coordinate updates (fully fixed atoms remain fixed; z-fixed atoms keep fixed z).
- Production probabilistic SC-environment coupling will cap LJ and Coulomb force magnitudes (configurable via `hybrid_control` attributes) to prevent force spikes.
- Production probabilistic SC-environment coupling will run an inner SC-only relaxation loop (`200` steps by default) with fixed environment/backbone and SC restrained near mapped rotamer positions; only the final relaxed-step SC->BB force is projected back to backbone.
- Inner SC relaxation integration treats `dE/dx` as gradient and advances positions using force direction (`x += dt * (-dE/dx)`), with restraint contribution in the same force convention.
- Production stage file preparation (`run_sim_1rkl.sh`, `test_prod_run_sim_1rkl.sh`) will explicitly set SC coupling control attrs in `/input/hybrid_control` (`sc_env_lj_force_cap`, `sc_env_coul_force_cap`, `sc_env_relax_steps`, `sc_env_relax_dt`, `sc_env_restraint_k`, `sc_env_max_displacement`) instead of relying on runtime defaults.
- Production SC rigid mapping in C++ must treat `sc_local_pos` as centroid-relative coordinates and re-add target centroid after rotation (`mapped = R*local + tgt_centroid`) to avoid frame-translation artifacts.
- Production SC row handling in C++ must enforce per-residue proxy capacity from placement-point support and select the safest proxy subset each step (by environment-distance metric) when proxy rows exceed representable capacity.
- For residues whose reference placement groups have fewer than 3 points (underdetermined rigid frame), SC local offsets are refreshed each step from current proxy coordinates.
- Underdetermined-group SC target displacement relative to current proxy is capped by `sc_env_max_displacement` during mapped-position construction to prevent large rotamer-shift jumps.
- BB force transfer no longer allows BB self-target projection: runtime BB COM refresh and BB gradient projection skip `atom_indices == bb_atom_index`, and preparation emits mass-fraction component mappings for four BB components with nearest non-BB protein fallback only when direct all-atom component indices are unavailable in the active runtime index space.
- Hybrid stage-file injection will materialize protein-AA backbone coordinate carriers (indexed by preserved PDB backbone indices via a fixed runtime offset), convert `hybrid_bb_map` active targets from PDB index space to these runtime carriers, and use those carriers for cycle-wise BB recomputation and BB-force projection.
- Production stage preparation must generate Upside backbone force-field nodes directly inside the hybrid stage `.up` (single-file workflow) by invoking `upside_config` writer functions in-process; avoid generating/consuming a separate backbone reference `.up` template.
- Production hybrid SC force transfer must target injected AA backbone carrier atoms (`N/CA/C/O`) in runtime index space; projection to MARTINI `BB`/`SC` atoms is disallowed.
- Stage handoff (`set_initial_position.py`) must realign injected AA carrier coordinates to current stage BB proxy coordinates before production starts to avoid BB/carrier frame discontinuities.
- Stage handoff (`set_initial_position.py`) must also recenter production-stage protein BB to box center (periodic circular-center on BB coordinates) and wrap coordinates back into the box before production integration.
- Production hybrid integration applies in-place RMSD alignment of protein coordinates (BB-based Kabsch/Horn rigid transform) once at the start of each C++ integration cycle before force evaluation/update, so Upside backbone updates and hybrid coupling operate in a rotation-stabilized frame.
- Production-stage timestep for hybrid runs with injected Upside all-atom backbone nodes must default to `0.002` (not MARTINI-only `0.020`) to maintain integration stability.
- During preparation, all-atom backbone reference coordinates (`N/CA/C/O`) are RMSD-aligned (Kabsch on residue backbone COMs) to MARTINI `BB` positions before writing `hybrid_bb_map/reference_atom_coords`.
- New workflow variant requirement: keep `hybrid_control` enabled for dynamic rigid-mask enforcement but set `activation_stage` to a non-used token so hybrid protein coupling is never active in any stage; this preserves rigid-protein dry-MARTINI interactions without production Upside force exchange.
- New workflow variant requirement: do not inject production rotamer/backbone nodes and use VTF extraction mode `1` at stage 7.0 to disable SC/backbone back-mapping output.
- Stage-file injection in `run_sim_1rkl.sh` must scale AA backbone carrier masses by `72/54` so each `N/CA/C/O` carrier set sums to MARTINI `BB` mass `72/12`, giving per-carrier masses `N=1.56`, `CA=1.33`, `C=1.33`, `O=1.78` (Upside mass units).

## Active Task (Radial Backbone Cross Table)

### Goal
- Implement a true distance-dependent dry-MARTINI vs Upside-backbone interaction-table workflow that uses the DOPC dry-type trend from `outputs/depth_interaction_table.csv` to span all dry-MARTINI types and emits a runtime-ready artifact for `martini_rbm_cross_potential`.

### Architecture & Key Decisions (Active Task)
- Keep the existing depth-table generator unchanged; use its DOPC bead-type aggregates as calibration inputs for the new cross table.
- Runtime form is the existing Gaussian radial-basis node `martini_rbm_cross_potential`; no new C++ potential node is added in this run.
- Backbone classes are the explicit injected AA roles `N`, `CA`, `C`, `O`.
- Environment classes are dry-MARTINI `/input/type` values and must cover all 38 dry atom types from `ff_dry/dry_martini_v2.1.itp`.
- The target radial curve is arbitrary smooth spline data, not LJ; dry-MARTINI `sigma/epsilon` are used only to span trend, well position, and role-specific scaling.
- DOPC calibration uses bead-type aggregates from the existing depth table and sets both ranking and scale for extrapolation to uncovered dry types.

### Execution Phases (Active Task)
- [x] Phase H1: Add radial cross-table generator and artifact writer.
- [x] Phase H2: Add standalone injector for prepared `.up` files using the generated artifact.
- [x] Phase H3: Extend validation tooling for the new CSV/JSON/HDF5 artifacts.
- [x] Phase H4: Run syntax and generation/injection verification and record results.

### Known Errors / Blockers (Active Task)
- None yet.

### Review (Active Task)
- Added new tooling:
  - `build_backbone_cross_interaction_table.py`
  - `inject_backbone_cross_potential.py`
  - extended `validate_backbone_only_up.py` for cross-table CSV/JSON/HDF5 validation
- Generated artifacts:
  - `outputs/backbone_cross_interaction_table.csv`
  - `outputs/backbone_cross_interaction_table.meta.json`
  - `outputs/backbone_cross_interaction_table.h5`
- Verification results:
  - `python3 -m py_compile build_backbone_cross_interaction_table.py inject_backbone_cross_potential.py validate_backbone_only_up.py` passed.
  - `python3 build_backbone_cross_interaction_table.py` passed and produced a `4 x 38 x 12` runtime weight tensor.
  - `python3 validate_backbone_only_up.py --cross-table-csv ... --cross-table-meta ... --cross-artifact ...` passed.
  - Injector check on a temporary copy of `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.prepared.up` passed and wrote a `martini_rbm_cross_potential` node with:
    - `protein_atom_indices=124`
    - `env_atom_indices=4025`
    - `weights_shape=(4, 38, 12)`

## Active Task (Dual Dry-MARTINI / Upside Runtime Wiring)

### Goal
- Update `run_sim_1rkl.sh` so the production-stage workflow actually runs the intended mixed simulation:
  - dry-MARTINI potentials for bilayer/environment
  - Upside all-atom-backbone potentials for protein
  - `martini_rbm_cross_potential` from the new backbone cross-table artifact for protein-environment coupling

### Architecture & Key Decisions
- Keep pre-production stages `6.0-6.6` unchanged: dry-MARTINI preparation plus hybrid mapping only, with no injected Upside-backbone or cross-potential nodes.
- Build the depth-table and radial backbone cross-table artifacts once per workflow run, immediately after hybrid packing/mapping generation and before stage file preparation.
- Inject both production-only node sets into stage `7.0` prepared files:
  - backbone-only Upside nodes from `inject_backbone_only_nodes.py`
  - cross potential from `inject_backbone_cross_potential.py`
- Rely on existing `aa_backbone_explicit_lipid` runtime behavior so `martini_potential` continues to cover bilayer/environment while suppressing protein-internal MARTINI interactions during production.
- Validate only production-stage prepared files against the backbone/cross-table validator, because earlier stages intentionally lack the production-only nodes.

### Execution Phases
- [x] Phase I1: Finish `run_sim_1rkl.sh` call-site wiring for backbone cross-table artifact generation and production-node injection.
- [x] Phase I2: Update workflow summary text/comments so stage `7.0` semantics are explicit.
- [x] Phase I3: Run shell syntax and targeted artifact/injection validation, then record results.

### Known Errors / Blockers
- None yet.

### Review
- `run_sim_1rkl.sh` now calls `prepare_backbone_cross_artifacts` immediately after hybrid packing/mapping preparation, so depth-table and radial cross-table artifacts are built before any stage files are prepared.
- Production-stage preparation (`stage_label=production`) now injects both backbone-only Upside nodes and `martini_rbm_cross_potential`; validation passes cross-table flags only when `BACKBONE_CROSS_ENABLE=1`.
- Workflow messaging now makes the mixed runtime explicit:
  - banner text updated to `Dual Dry-MARTINI / Upside`
  - stage `7.0` comment updated to `dry bilayer + Upside backbone + cross table`
  - completion summary now reports the generated depth/cross artifact paths
- Verification completed:
  - `bash -n run_sim_1rkl.sh` passed.
  - `build_depth_interaction_table.py` and `build_backbone_cross_interaction_table.py` succeeded writing artifacts under `outputs/martini_test_1rkl_hybrid/`.
  - `validate_backbone_only_up.py --table-* --cross-*` passed on those artifacts.
  - Production-path smoke test passed on a temporary copy of `1rkl.stage_7.0.prepared.up` after applying the same injector order as the script (`inject_backbone_only_nodes.py` then `inject_backbone_cross_potential.py`).

## Active Task (Cached MARTINI Hybrid Parameter File)

### Goal
- Move the runtime dry-MARTINI/Upside cross-interaction artifact into the shared parameter tree as `parameters/ff_2.1/martini.h5` so normal simulation runs reuse it instead of rebuilding it every time.

### Architecture & Key Decisions
- Treat `parameters/ff_2.1/martini.h5` as the default runtime artifact consumed by `inject_backbone_cross_potential.py`.
- Keep the CSV/JSON outputs as optional build products in the run directory for inspection and regeneration only; they are no longer required for every run.
- Default workflow behavior is cache-first:
  - if `parameters/ff_2.1/martini.h5` exists, reuse it
  - rebuild only when the file is missing or an explicit rebuild flag is set
- Production-stage validation must tolerate cached-only runs by passing `--cross-artifact` always and only adding `--cross-table-csv/--cross-table-meta` when those files exist.

### Execution Phases
- [x] Phase J1: Update `run_sim_1rkl.sh` defaults and artifact-preparation logic for cached `parameters/ff_2.1/martini.h5`.
- [x] Phase J2: Generate/install `parameters/ff_2.1/martini.h5`.
- [x] Phase J3: Re-run syntax and targeted artifact/injection validation, then record results.

### Known Errors / Blockers
- None yet.

### Review
- `run_sim_1rkl.sh` now defaults `BACKBONE_CROSS_TABLE_H5` to `parameters/ff_2.1/martini.h5` via `${UPSIDE_HOME}/parameters/ff_2.1/martini.h5`.
- Added cache-first behavior:
  - if `martini.h5` exists and `BACKBONE_CROSS_REBUILD!=1`, the workflow reuses it and skips rebuilding depth/cross artifacts
  - if the file is missing or `BACKBONE_CROSS_REBUILD=1`, the workflow rebuilds and writes the HDF5 artifact back to `parameters/ff_2.1/martini.h5`
- Production-stage validation now always uses `--cross-artifact` and only passes `--cross-table-csv/--cross-table-meta` when those optional files exist, so cached-only runs do not fail due to missing run-local reports.
- Installed shared runtime artifact:
  - `parameters/ff_2.1/martini.h5`
- Verification completed:
  - `bash -n run_sim_1rkl.sh` passed.
  - `build_backbone_cross_interaction_table.py ... --output-h5 parameters/ff_2.1/martini.h5` passed.
  - `validate_backbone_only_up.py --cross-artifact parameters/ff_2.1/martini.h5` passed.
  - Production-path smoke test on a temp `.up` copy passed using the shared artifact path:
    - `inject_backbone_only_nodes.py`
    - `inject_backbone_cross_potential.py --artifact parameters/ff_2.1/martini.h5`
    - `validate_backbone_only_up.py ... --cross-artifact parameters/ff_2.1/martini.h5`

## Active Task (Stage-7 Runtime Compatibility Fixes)

### Goal
- Fix the actual stage-7.0 runtime startup errors in the backbone-only + cross-potential production path so `run_sim_1rkl.sh` can launch MD without failing on missing node dependencies or cross-node shape mismatches.

### Architecture & Key Decisions
- `affine_alignment` is a required backbone dependency for `backbone_pairs`; it must be kept/generated in backbone-only production stage files and must not be classified as a forbidden sidechain node.
- `martini_rbm_cross_potential` must use dense local class indices inside each injected `.up` file, because the C++ loader infers class counts from the max class index present in that file rather than from the artifact-wide class list.
- Validation must check runtime-relevant graph/schema constraints, not only artifact presence:
  - require `affine_alignment` alongside the other backbone nodes
  - verify cross-node `weights.shape == (max(protein_class_index)+1, max(env_class_index)+1, n_radial)`

### Execution Phases
- [x] Phase K1: Restore `affine_alignment` in backbone-only node injection and align validation with that dependency.
- [x] Phase K2: Dense-reindex injected cross-potential classes and update validation for runtime-compatible `weights` shapes.
- [x] Phase K3: Re-run temp stage-7.0 injection plus one-step `upside` startup verification.

### Known Errors / Blockers
- One-step startup now succeeds, but the temp production smoke test still reports an extremely large initial `martini_potential`; that is a physical/calibration issue, not the startup/schema failure reported here.

### Review
- `inject_backbone_only_nodes.py` now treats `affine_alignment` as part of the required backbone node set and explicitly writes it before `backbone_pairs`.
- `validate_backbone_only_up.py` now requires `affine_alignment` and validates injected `martini_rbm_cross_potential` node dimensions against the actual class indices stored in the stage file.
- `inject_backbone_cross_potential.py` now compresses used protein/environment classes to dense local indices and slices the weight tensor accordingly before writing the node into a `.up` file.
- Verification completed on a temp copy of `1rkl.stage_7.0.prepared.up`:
  - backbone-only injection + cross injection passed validation
  - injected cross node now uses `weights.shape == (4, 6, 12)` for the six dry types present in this system
  - one-step `obj/upside` startup passed the previous `backbone_pairs -> affine_alignment` and `weights shape` failures and entered MD

## Active Task (Pre-Production Rigid Control Restoration)

### Goal
- Restore the explicit stage-file control metadata in `run_sim_1rkl.sh` so stages `6.0-6.6` actually hold the protein rigid before the production stage while keeping stage `7.0` activation semantics unchanged.

### Architecture & Key Decisions
- Keep the current backbone-only mapping schema (`hybrid_bb_map` + `hybrid_env_topology`) as the required mapping payload.
- Treat `/input/hybrid_control` as runtime control metadata that the workflow synthesizes directly into each prepared stage file instead of requiring it from the mapping export.
- Preserve any optional `hybrid_control` group from the mapping file if present, then overwrite the workflow-critical attrs so behavior is deterministic:
  - `enable=1`
  - `activation_stage=production`
  - `preprod_protein_mode=rigid`
  - `prep_runtime_mode=dry_martini_prep`
  - `active_runtime_mode=aa_backbone_explicit_lipid`
  - `preprod_lipid_headgroup_roles=PO4`
- Limit the fix to `run_sim_1rkl.sh`; do not change C++ runtime behavior for this request.

### Execution Phases
- [x] Phase L1: Patch `run_sim_1rkl.sh` to write the required `hybrid_control` attrs into every prepared stage file.
- [x] Phase L2: Run syntax plus temp-HDF5 verification to confirm the generated control block is present and carries the intended values.
- [x] Phase L3: Record the root cause and verification results in `findings.md` and `progress.md`.

### Known Errors / Blockers
- Existing generated stage files under `outputs/martini_test_1rkl_hybrid/checkpoints/` were produced before this fix and currently lack `/input/hybrid_control`; they must be regenerated to benefit from the restored rigid-preproduction behavior.

### Review
- `run_sim_1rkl.sh` now restores `/input/hybrid_control` directly into every prepared stage file after mapping injection and stage-label assignment.
- The workflow writes the control attrs needed by the existing C++ rigid-mask logic:
  - `enable=1`
  - `activation_stage=production`
  - `preprod_protein_mode=rigid`
  - `prep_runtime_mode=dry_martini_prep`
  - `active_runtime_mode=aa_backbone_explicit_lipid`
  - `preprod_lipid_headgroup_roles=PO4`
- `inject_hybrid_mapping()` also preserves an optional mapping-provided `hybrid_control` group if one appears in future artifacts, but the workflow still overwrites the required control attrs so stage semantics stay deterministic.
- Verification completed:
  - `bash -n run_sim_1rkl.sh` passed.
  - A temp copy of `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.6.up` was updated via the new `set_hybrid_control()` helper from the script itself and then inspected with `h5py`; the resulting attrs matched the intended pre-production rigid / production activation settings while `current_stage` remained `minimization`.

## Active Task (Explicit Pre-Production `fix_rigid` Enforcement)

### Goal
- Ensure `run_sim_1rkl.sh` actually keeps the protein backbone rigid through stages `6.0-6.6` in the currently active runtime, so stage `7.0` starts from an undeformed protein.

### Architecture & Key Decisions
- Root cause: the workflow was writing `/input/hybrid_control`, but the currently active runtime path no longer registers hybrid state from `src/main.cpp`, so those attrs were not consumed during stage `6.x` execution.
- Use the already active `/input/fix_rigid` path instead of re-enabling dormant hybrid registration:
  - pre-production stage files (`stage_label != production`) write `/input/fix_rigid` with `atom_indices = where(protein_membership >= 0)`
  - production stage files write `/input/fix_rigid.enable = 0` and remove `atom_indices`
- Keep the `hybrid_control` attrs in place for schema/runtime intent, but do not rely on them for pre-production rigidity in this workflow.
- Limit the fix to `example/16.MARTINI/run_sim_1rkl.sh`; do not change C++ code for this request.

### Execution Phases
- [x] Phase M1: Add explicit pre-production `fix_rigid` stage-file injection in `run_sim_1rkl.sh`.
- [x] Phase M2: Validate generated `fix_rigid` metadata for both pre-production and production stage semantics.
- [x] Phase M3: Run a short pre-production MD probe and confirm zero protein displacement under the current binary.

### Known Errors / Blockers
- Existing generated stage files under `outputs/martini_test_1rkl_hybrid/checkpoints/` were produced before this fix and must be regenerated to pick up explicit `fix_rigid` enforcement.

### Review
- `run_sim_1rkl.sh` now writes `/input/fix_rigid` directly from `hybrid_env_topology/protein_membership` for pre-production stages and disables it for production.
- This path uses the existing active runtime hooks (`register_fix_rigid_for_engine`, `apply_fix_rigid_md`, and barostat fixed-atom skipping in `src/box.cpp`) instead of the inactive hybrid-control registration path.
- Verification completed:
  - `bash -n run_sim_1rkl.sh` passed.
  - temp pre-production stage-file metadata probe on `1rkl.stage_6.4.prepared.up` showed:
    - `/input/fix_rigid.enable = 1`
    - `atom_indices.shape[0] = 298`
  - temp production stage-file metadata probe on `1rkl.stage_7.0.prepared.up` showed:
    - `/input/fix_rigid.enable = 0`
    - no `atom_indices` dataset
  - short MD probe on a temp copy of `1rkl.stage_6.4.prepared.up` with `fix_rigid` enabled:
    - `duration-steps = 20`
    - protein displacement from `/input/pos` to last `/output/pos` frame was exactly zero (`max_abs_displacement = 0.0` over 298 protein atoms)

## Active Task (C++ Rigid-Path Audit)

### Goal
- Confirm whether any rigid-hold code still exists in `src/*.cpp`, and distinguish active runtime hooks from dead/unregistered code paths.

### Architecture & Key Decisions
- Limit this task to source inspection; do not change C++ behavior in this audit.
- Report the runtime split explicitly:
  - live `fix_rigid` path
  - dead hybrid-driven pre-production rigid path
- Use exact source references so the next code change can target the correct call site instead of guessing from workflow metadata.

## 2026-03-06 (Refine Dry-MARTINI <-> Upside Cross Force Field)

### Project Goal (Current)
- Replace the current approximate RBM cross-force-field builder with a refined dry-MARTINI-to-Upside workflow that:
  - computes `Q0/Qa/Na/C1/C3` anchor targets from empirical DOPC depth occupancy
  - audits but does not yet act on membrane asymmetry
  - spans those anchors to all dry-MARTINI types through a dry-table trend instead of per-type fitting
  - writes a true distance-dependent spline interaction for protein/environment pairs
  - normalizes dry-MARTINI masses by dividing by `72`, with matching BB-carrier normalization

### Architecture & Key Decisions (Current)
- Keep the active shell workflow surface (`prepare_system.py`, `run_sim_1rkl.sh`) and extend it in place.
- The refined build path uses empirical signed `z-z_center` samples from the DOPC template, not just per-type mean depths.
- Membrane asymmetry is measured from raw `membrane.h5` curves and recorded in metadata, but the parameter build still uses a symmetrized evaluation for this pass.
- Type spanning uses one global trend from dry-MARTINI interaction-table strength to Upside well depth; no per-type force fitting and no role-specific scaling.
- Runtime cross interaction moves from `martini_rbm_cross_potential` to a new spline-table node; backward compatibility for legacy artifacts may remain in code, but the workflow should emit the new spline schema by default.
- Distance dependence is represented as direct radial spline tables over a uniform distance grid.
- The next refinement will no longer force identical backbone-class tables:
  - `N` will use donor-channel anchor targets,
  - `O` will use acceptor-channel anchor targets,
  - `CA/C` will use burial-channel anchor targets.
- Dry-MARTINI environment masses use `mass/72`; injected AA backbone carrier masses are normalized so each `N+CA+C+O` carrier set sums to `1.0`.
- Parameter calculation and publication to `parameters/ff_2.1/martini.h5` must live in one standalone script under `example/16.MARTINI/` and must not be performed inside the simulation workflow.

### Execution Phases (Current Run)
- [x] Phase R1: Extend the depth-table builder with empirical sample pooling and symmetry-audit metadata.
- [x] Phase R2: Replace the Python cross-table builder/output schema with global-trend spline tables and normalized carrier masses.
- [x] Phase R3: Add the new C++ spline cross-potential node and injector/runtime wiring.
- [x] Phase R4: Normalize dry-MARTINI stage masses to `mass/72` and align injected carrier masses to the same convention.
- [x] Phase R5: Run targeted validation: Python syntax, C++ build, artifact validation, injection/startup smoke test, and mass sanity checks.

### Known Errors / Blockers (Current Run)
- The existing production physical-instability issue may still remain after the force-field representation change; this run only replaces the builder/runtime form and mass convention, not full physical calibration.
- New user finding: the force field performs especially poorly for backbone `O`. This likely invalidates the current simplifying decision that all four backbone classes (`N/CA/C/O`) should share identical spline tables; the next refinement should split `O` back out, and likely reintroduce at least partial donor/acceptor role specificity.
- New user requirement: parameter generation must be decoupled from `run_sim_1rkl.sh`; the workflow should consume a prebuilt `parameters/ff_2.1/martini.h5` instead of rebuilding it on demand.

### Review (Current Run)
- `example/16.MARTINI/prepare_system.py` now:
  - builds depth rows from empirical signed depth samples instead of `V(mean_depth)`,
  - records `type_sample_stats`, raw sample lists, and a raw membrane `V(+d)` vs `V(-d)` asymmetry audit,
  - spans all 38 dry types from separate `Q0/Qa/Na/C1/C3` anchors using the dry-table descriptor `mean_epsilon_to_all_dry_types`,
  - writes spline-table artifacts with schema `martini_backbone_spline_cross_v1`,
  - validates and injects both legacy RBM and new spline artifacts, removing legacy cross nodes during overwrite.
- `src/martini.cpp` now registers `martini_spline_cross_potential`, which loads direct radial tables, fits clamped 1D splines at startup, and evaluates protein-environment energy/forces from those spline tables.
- `example/16.MARTINI/lib.py` now normalizes dry-MARTINI environment masses with `mass/72`.
- `example/16.MARTINI/run_sim_1rkl.sh` now:
  - normalizes injected AA carrier masses as `14/54`, `12/54`, `12/54`, `16/54`,
  - consumes a prebuilt `parameters/ff_2.1/martini.h5` instead of building parameters inline,
  - errors out with a pointer to the standalone builder script if the artifact is missing or has the wrong schema.
- `example/16.MARTINI/build_martini_parameters.py` now provides the standalone parameter-build path:
  - builds the empirical depth table,
  - builds the spline cross artifact,
  - validates the intermediate CSV/JSON/HDF5 artifacts,
  - publishes the final runtime artifact to `parameters/ff_2.1/martini.h5` or another requested output path.
- Validation completed in `example/16.MARTINI/outputs/test_refined_cross_20260306`:
  - `python3 -m py_compile example/16.MARTINI/prepare_system.py example/16.MARTINI/lib.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `cmake --build obj -j 4`
  - rebuilt `depth_interaction_table.*` and `backbone_cross_interaction_table.*`
  - `validate-backbone-only` passed on the rebuilt CSV/JSON/HDF5 artifacts
  - injected `martini_spline_cross_potential` into a copied `1rkl.stage_7.0.prepared.up` and `validate-backbone-only` passed on the resulting stage file
  - one-step `obj/upside` startup/run completed on the copied stage file with the new spline node active
- Canonical runtime artifact update completed:
  - copied `example/16.MARTINI/outputs/test_refined_cross_20260306/backbone_cross_interaction_table.h5`
    to `parameters/ff_2.1/martini.h5`
  - destination validates as `martini_backbone_spline_cross_v1`
  - destination SHA1 matches the rebuilt refined artifact exactly
- Standalone-builder validation completed:
  - `python3 -m py_compile example/16.MARTINI/build_martini_parameters.py`
  - `python3 example/16.MARTINI/build_martini_parameters.py --output-dir example/16.MARTINI/outputs/test_parameter_builder_script_20260306 --publish-output example/16.MARTINI/outputs/test_parameter_builder_script_20260306/martini.published.h5`
  - published test artifact validates as `martini_backbone_spline_cross_v1` with `38` environment classes
- Role-specific backbone refinement completed:
  - rebuilt and published `parameters/ff_2.1/martini.h5` from `example/16.MARTINI/outputs/test_o_role_refine_20260306/backbone_cross_interaction_table.h5`
  - builder now calibrates `N` from donor-channel anchors, `O` from acceptor-channel anchors, and `CA/C` from burial-channel anchors
  - injected the updated published artifact into a copied `1rkl.stage_7.0.prepared.up`, `validate-backbone-only` passed, and one-step `obj/upside` startup/run still completed
- Verified diagnostics:
  - artifact schema: `martini_backbone_spline_cross_v1`
  - `O` no longer shares the same spline table as `N/CA/C`
  - `CA` and `C` remain intentionally identical
  - example published-artifact differences:
    - `max |N-O|` for `Q0` table = `2.1418`
    - `max |CA-C|` for `Q0` table = `0.0`
    - `Q0` well minima: `N=-1.2878`, `O=-3.4296`
  - anchor descriptors remain distinct: `Q0=1.0526`, `Qa=1.1421`, `Na=1.8000`, `C3=2.5447`, `C1=2.7447`
  - symmetry audit worst case from the current membrane model: `cb:ARG:burial0`, `max_abs_diff=5.32998` at `13.1 A`
  - mass normalization sanity: `72 -> 1.0`, `45 -> 0.625`, carrier mass sum `N+CA+C+O = 1.0`

### Execution Phases
- [x] Phase N1: Search `src/*.cpp` / `src/*.h` for rigid/fix/hybrid stage-control hooks.
- [x] Phase N2: Inspect exact call sites in `src/main.cpp`, `src/martini.cpp`, `src/box.cpp`, and `src/deriv_engine.cpp`.
- [x] Phase N3: Record which rigid mechanisms are active vs inactive.

### Review
- Live rigid-hold code exists in C++:
  - `src/martini.cpp`: `martini_fix_rigid` registry, `/input/fix_rigid` reader, and force/momentum zeroing
  - `src/main.cpp`: registration of `fix_rigid` and repeated `apply_fix_rigid_md(...)` calls inside MD
  - `src/box.cpp`: fixed-atom skipping during barostat scaling and pressure estimation
  - `src/deriv_engine.cpp`: fixed-atom mask helpers
- The hybrid-driven pre-production rigid path also exists in `src/martini.cpp`, but it is not currently wired from `src/main.cpp`:
  - `martini_hybrid::register_hybrid_for_engine(...)` exists
  - `martini_hybrid::update_stage_for_engine(...)` exists
  - `src/main.cpp` does not declare or call `register_hybrid_for_engine(...)`
- Conclusion:
  - yes, there is hold-rigid code in `src/*.cpp`
  - the active runtime-enforced mechanism is `fix_rigid`
  - the hybrid-metadata-driven rigid path is currently dead unless `register_hybrid_for_engine(...)` is reconnected in `src/main.cpp`

## Active Task (Workflow-Level Rigid + Handoff Guarantees)

### Goal
- Make `example/16.MARTINI/run_sim_1rkl.sh` prove during execution that the protein stays rigid in every pre-production stage and that the rigid protein state is handed unchanged into stage `7.0`.

### Architecture & Key Decisions
- Keep `/input/fix_rigid` as the runtime enforcement path.
- Add workflow-level assertions instead of relying on offline inspection:
  - `assert_protein_rigid_stage(...)` checks `/input/pos` vs last `/output/pos` for protein atoms after each stage `6.0-6.6`
  - `assert_production_handoff_protein_particles(...)` checks `6.6` last-frame vs `7.0` input for original protein particles only, excluding intentionally refreshed appended AA carrier atoms
- Reuse the same tolerance knob for both checks:
  - `PREPROD_RIGID_ASSERT_ENABLE=1`
  - `PREPROD_RIGID_TOL=1e-6`

### Execution Phases
- [x] Phase O1: Add per-stage rigid assertions to `run_sim_1rkl.sh`.
- [x] Phase O2: Add a production-handoff assertion that ignores refreshed appended carrier atoms and checks original protein particles only.
- [x] Phase O3: Run reduced end-to-end workflow tests and verify all rigid/handoff assertions pass.

### Review
- `run_sim_1rkl.sh` now contains:
  - `assert_protein_rigid_stage(...)`
  - `assert_production_handoff_protein_particles(...)`
- Reduced end-to-end workflow run passed:
  - command env:
    - `RUN_DIR=outputs/test_fix_rigid_verify2`
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
- Assertions observed during the real script run:
  - `Rigid-protein check stage 6.0 ... max_disp=0`
  - `Rigid-protein check stage 6.1 ... max_disp=0`
  - `Rigid-protein check stage 6.2 ... max_disp=0`
  - `Rigid-protein check stage 6.3 ... max_disp=0`
  - `Rigid-protein check stage 6.4 ... max_disp=0`
  - `Rigid-protein check stage 6.5 ... max_disp=0`
  - `Rigid-protein check stage 6.6 ... max_disp=0`
  - `Production handoff protein-particle check ... max_disp=0`
- Fresh rerun against the current script state also passed with the same proof lines:
  - `RUN_DIR=outputs/test_fix_rigid_verify3`
  - stages `6.0-6.6` again reported `max_disp=0`
  - the `6.6 -> 7.0` handoff again reported `max_disp=0`
  - the workflow completed through stage `7.0`
- Supporting reduced-run artifacts:
  - `example/16.MARTINI/outputs/test_fix_rigid_verify/`
  - `example/16.MARTINI/outputs/test_fix_rigid_verify2/`
  - `example/16.MARTINI/outputs/test_fix_rigid_verify3/`

## Active Task (Production-Start Protein Representation)

### Goal
- Make the production-stage workflow display and validate the actual folded AA backbone carriers at stage start, instead of a residue-wise BB reconstruction that can make the protein appear broken.

### Architecture & Key Decisions
- Treat the user-visible production-start representation as the acceptance criterion.
- Prefer actual appended `PROTEINAA` carrier coordinates when present in a stage file.
- Keep BB-based reconstruction only as a fallback for older files that do not contain appended AA carriers.
- Verify both:
  - the stage `7.0` `.up` input contains the expected carrier coordinates
  - the stage `7.0` VTF/PDB export path uses those carrier coordinates instead of reconstructing from BB anchors

### Execution Phases
- [ ] Phase P1: Audit `extract_martini_vtf.py` mode-2 protein path against production-stage `.up` carrier coordinates.
- [ ] Phase P2: Patch export logic to prefer actual appended AA carriers when available.
- [ ] Phase P3: Rerun a reduced workflow and verify the exported production-start protein matches the stage `7.0` carrier coordinates.

### Known Errors / Blockers
- None yet.

## 2026-03-06 (Fresh Double-Check: Pre-Production Rigid Hold + Unchanged Production Handoff)

### Goal
- Re-verify the current `example/16.MARTINI/run_sim_1rkl.sh` state on a fresh reduced rerun.
- Prove two things on newly generated artifacts:
  - the protein remains rigid through stages `6.0-6.6`
  - the original protein coordinates handed into stage `7.0` are unchanged from the `6.6` last frame

### Architecture & Key Decisions
- Do not patch the workflow unless the live rerun or direct HDF5 inspection finds a regression.
- Use both:
  - the script's built-in assertions (`assert_protein_rigid_stage`, `assert_production_handoff_protein_particles`)
  - an independent post-run HDF5 check over the generated checkpoint files

### Execution Phases
- [x] Phase Q1: Re-read the live workflow code paths for `fix_rigid`, stage handoff, and production preparation.
- [x] Phase Q2: Run a fresh reduced `6.0 -> 7.0` workflow from the current script.
- [x] Phase Q3: Independently inspect the generated HDF5 stage files for rigid hold and exact `6.6 -> 7.0` protein-coordinate equality.

### Known Errors / Blockers
- The separate late-stage MARTINI environment instability still produces very large MARTINI energies, but it did not break the rigid-protein guarantee in this verification run.

### Review
- Fresh reduced rerun completed through stage `7.0`:
  - run directory: `example/16.MARTINI/outputs/test_double_check_rigid_20260306/`
  - overrides:
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
- Built-in workflow proof lines from the fresh run:
  - `Rigid-protein check stage 6.0: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.1: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.2: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.3: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.4: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.5: n_protein=124 max_disp=0`
  - `Rigid-protein check stage 6.6: n_protein=124 max_disp=0`
  - `Production handoff protein-particle check: n_protein=124 max_disp=0`
- Independent HDF5 post-check on the same artifacts confirmed:
  - stages `6.0-6.6` each have `/input/fix_rigid enable=1`
  - stages `6.0-6.6` each have `n_protein=124` and exact `input -> last output` protein-coordinate equality (`max_disp=0`)
  - stage `7.0` has `production_fix_rigid=0` by design
  - the original protein particles in `stage_7.0.up` match the `stage_6.6.up` last output exactly (`max_disp=0`)
- Conclusion:
  - current `run_sim_1rkl.sh` still holds the protein rigid before production
  - the protein coordinates handed into production are unchanged
  - no workflow code change was required
