
## 2026-04-28 Hybrid MARTINI Workflow Refactor

### Project Goal
- Decouple `example/16.MARTINI/run_sim_1rkl.sh` orchestration from low-level MARTINI constants, stabilization heuristics, and HDF5 mutation details.
- Make `py/martini_prepare_system.py` own the typed workflow defaults and stage metadata writes through argparse.
- Preserve the current hybrid MARTINI physics: SC-env and BB-env interface interactions must remain active in production.

### Architecture & Key Decisions
- Rewrite `run_sim_1rkl.sh` as a thin launcher that only bootstraps the environment and forwards high-level workflow controls.
- Add a `run-hybrid-workflow` Python command that owns stage prep, stage execution, continuation naming, HDF5 stage mutation, SC library validation/build, and VTF extraction dispatch.
- Move the named force caps, relaxation ramps, backbone-frame heuristics, NPT internals, Ewald alpha, and unit constants into Python argparse defaults.
- Remove dead shell-only low-level variables rather than preserving inert compatibility surfaces.
- Keep production non-protein hard-sphere replacement disabled by default and do not add any switch that disables SC-env or BB-env hybrid interface interactions.

### Execution Phases
- [x] Phase 1: Add Python workflow config/orchestration and argparse defaults.
- [x] Phase 2: Replace the base bash workflow with a high-level launcher.
- [x] Phase 3: Update associated wrappers to expose only high-level defaults.
- [x] Phase 4: Verify syntax, parser surface, code search, and reduced workflow behavior where feasible.
- [x] Phase 5: Record review results and remaining blockers.

### Known Errors / Blockers
- None identified at task start.

### Review
- Implemented:
  - `py/martini_prepare_system.py run-hybrid-workflow` now owns hybrid packing, stage-file conversion, HDF5 metadata injection, continuation naming, stage execution, SC-library validation/build, and VTF extraction.
  - `run_sim_1rkl.sh` is now a high-level launcher with orchestration controls only; low-level force caps, ramps, NPT/Ewald internals, and unit constants are argparse defaults in Python.
  - `run_sim_1rkl_outlipid.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_outlipid.sh` are thin wrappers that set high-level system/placement/run defaults.
  - production stage files now write migrated hybrid-control attrs from Python and explicitly keep `production_nonprotein_hard_sphere=0`.
- Physics integrity audit:
  - reduced full workflow generated a production stage containing both `martini_potential` and `martini_sc_table_1body`;
  - production `hybrid_control` had `activation_stage=production`, `sc_env_lj_force_cap=25`, `sc_env_coul_force_cap=25`, `sc_env_relax_steps=150`, `sc_env_backbone_hold_steps=200`, `sc_env_po4_z_hold_steps=150`, `sc_env_po4_z_clamp_enabled=1`, and `nonprotein_hs_force_cap=100`;
  - production environment membership retained `4025` non-protein targets.
- Verification:
  - `bash -n` passed for all four workflow shell entrypoints.
  - `.venv/bin/python -m py_compile py/martini_prepare_system.py` passed.
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --help` exposed the migrated argparse surface.
  - Search found no migrated low-level variable names in the bash wrappers or `run.py`.
  - Reduced continuation smoke wrote `/tmp/hybrid_refactor_continue/checkpoints/1rkl.stage_7.1.up` and VTF.
  - Reduced fresh workflow wrote `/tmp/hybrid_refactor_full/checkpoints/1rkl.stage_7.0.up` and VTF after completing stages `6.0` through `7.0`.

## 2026-04-28 Hybrid MARTINI First-Frame Geometry Regression

### Project Goal
- Restore stable 1RKL initial geometry after the workflow refactor; the first VTF frame must retain the proper helical backbone while preserving all active SC-env and BB-env hybrid interface interactions.

### Architecture & Key Decisions
- Treat the regression as a stage-file geometry/mapping issue unless runtime evidence proves otherwise.
- Keep the fix in `py/martini_prepare_system.py`; do not reintroduce low-level shell controls or wrapper compatibility paths.
- Preserve direct stage-runtime BB mapping and production SC/BB environment potentials.
- Compact appended AA reference atoms to the actual referenced N/CA/C/O atoms instead of preserving sparse raw AA PDB index gaps as runtime particles.
- Stage-7 CB placement must consume `hybrid_bb_map/atom_indices` after runtime remapping, not recompute runtime indices from raw `reference_atom_indices`.

### Execution Phases
- [x] Phase 1: Inspect generated stage-7 mapping metadata and identify the sparse-reference atom expansion.
- [x] Phase 2: Patch hybrid mapping injection to append only compact referenced AA atoms.
- [x] Phase 3: Verify reduced 1RKL workflow output and first-frame mapping integrity.
- [x] Phase 4: Record results and remaining risk.

### Known Errors / Blockers
- Reduced one-step workflows remain invalid as stability benchmarks because intentionally collapsed 6.x relaxation can leave enormous MARTINI energies before stage 7. Use them only for wiring and geometry checks.

### Review
- Fixed:
  - `py/martini_prepare_system.py` now appends only the compact set of referenced AA backbone atoms and remaps runtime BB carriers through that compact lookup.
  - `py/martini_prepare_system_lib.py` stage-7 CB placement now consumes the already-remapped `hybrid_bb_map/atom_indices` directly.
- Verification:
  - reduced fresh run completed through `/tmp/hybrid_refactor_fix_full2/checkpoints/1rkl.stage_7.0.up` and `/tmp/hybrid_refactor_fix_full2/1rkl.stage_7.0.vtf`;
  - production stage has `n_atom=4214`, `reference_index_count=124`, runtime BB indices `4090..4213`, and all runtime BB indices in bounds;
  - appended AA role counts are exactly `31` each for `N`, `CA`, `C`, and `O`;
  - production HDF5 still contains both `martini_potential` and `martini_sc_table_1body`, with `activation_stage=production` and `production_nonprotein_hard_sphere=0`;
  - stage-7 input backbone distances are normal before dynamics: N-CA mean `1.471 A`, CA-C mean `1.539 A`, C-O mean `1.234 A`, CA-CA next mean `3.906 A`.

## 2026-04-28 Hybrid MARTINI Stage-7 Explosion Regression

### Project Goal
- Stop the current stage-7 1RKL production blow-up while keeping SC-env and BB-env interface potentials active.

### Architecture & Key Decisions
- Treat the user-reported trajectory as authoritative: the workflow is still wrong if Rg jumps from `12.6 A` to `17862.9 A` by step `50`.
- Compare the runtime force carrier model against the last stable direct-Upside hybrid semantics before making more cleanup edits.
- Do not disable `martini_potential`, `martini_sc_table_1body`, SC-env, or BB-env interactions.
- Do not add debug bypass flags.

### Execution Phases
- [x] Phase 1: Re-open regression notes and capture the user correction.
- [x] Phase 2: Inspect stage-7 force carrier and pair semantics against a known stable artifact.
- [x] Phase 3: Patch the wrong carrier/integration path directly.
- [x] Phase 4: Replay a short production window and verify Rg and energies do not explode.

### Known Errors / Blockers
- User-reported production run explodes by step `50` after the compact mapping fix:
  - step `0`: `Rg 12.6 A`, total `-12776.13`;
  - step `50`: `Rg 17862.9 A`, total `1.111e+12`;
  - step `100`: `Rg 34374319.3 A`, total `2.324e+20`.

### Review
- Fixed two runtime issues:
  - `src/martini.cpp` now uses the remapped `hybrid_bb_map/atom_indices` as the runtime N/CA/C/O carrier indices for backbone O refresh instead of recomputing `reference_index_offset + raw_reference_index`;
  - `martini_sc_table_1body` and the legacy SC-table node now apply the configured SC-env force cap to table point/vector gradients whenever the cap is positive, keeping SC-env active but preventing launch-force startup events.
- Verification:
  - rebuilt `obj/upside` successfully;
  - replayed the saved failing handoff for `1000` steps from `/tmp/hybrid_stage7_replay_cap1000.up`;
  - replay stayed finite with CA Rg `12.5-12.7 A` through the full replay;
  - output HDF5 is all finite, with CA Rg `12.759 A` at frame `0` and `12.723 A` at the last saved frame;
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed;
  - `bash -n` passed for all four MARTINI workflow entrypoint scripts.

## 2026-04-28 Hybrid MARTINI Production Restart Parity

### Project Goal
- Fix stage-7 production restart behavior so split runs (`7.0 -> 7.1 -> ...`) preserve the same protein stability as an uninterrupted production run.

### Architecture & Key Decisions
- Restart must preserve physical state, not just positions: momenta/velocities, box, hybrid startup counters, and dynamic hold/ramp state must be equivalent.
- Continue using active SC-env and BB-env interactions; do not add restart-only debug bypasses.
- Compare split-vs-continuous production using the saved stage-7 handoff before declaring the bug fixed.

### Execution Phases
- [x] Phase 1: Document the restart regression and inspect continuation workflow.
- [x] Phase 2: Identify missing restart state by comparing uninterrupted and split runs.
- [x] Phase 3: Patch continuation state transfer.
- [x] Phase 4: Verify split-vs-continuous stability and log results.

### Known Errors / Blockers
- No active blocker remains for production restart state transfer.

### Review
- Fixed:
  - production MD runs now record momentum so stage files contain restartable `output/mom`;
  - handoff now copies the last saved `output/mom` to `input/mom` and marks it restart-valid;
  - continuation now passes `--restart-using-momentum` only when copied momentum is valid;
  - production continuation now preserves the hybrid SC-env transition counter with `hybrid_control.sc_env_transition_step_start`;
  - `src/martini.cpp` initializes `sc_env_transition_step` from that restart attr instead of always restarting the ramp/hold from zero.
- Verification:
  - direct continuous `1000`-step replay from the saved stage-7 handoff stayed stable with Rg around `12.7-12.9 A`;
  - split replay (`500 + 500`) stayed stable after restart with Rg around `12.8-13.0 A`;
  - actual workflow continuation wrote `/tmp/restart_workflow_check/checkpoints/1rkl.stage_7.1.up`, invoked `--record-momentum --restart-using-momentum`, and kept Rg finite around `12.9 A`;
  - generated continuation file has `input/mom.restart_valid=1`, `output/mom`, and `sc_env_transition_step_start=450`;
  - `cmake --build obj --target upside`, Python compile checks, and shell syntax checks passed.

## 2026-04-29 Hybrid MARTINI Restart Instability Still Reproduces

### Project Goal
- Resolve the still-observed production restart instability when starting from a previous production stage.

### Architecture & Key Decisions
- Treat the prior restart fix as incomplete.
- Audit whether continuation is using the true final integrator state or merely the last logged frame.
- If exact restart state is unavailable, fail clearly instead of silently rethermalizing or restarting from stale logged state.
- Production-to-production continuation must preserve saved coordinates and momenta exactly; hybrid reference carriers must not be rebuilt during restart handoff.
- MD output must include the true final state so a later continuation does not restart from the last pre-final logging interval.
- Preserve active SC-env and BB-env interactions.

### Execution Phases
- [x] Phase 1: Record unresolved restart instability correction.
- [x] Phase 2: Audit restart source state actually available in current production files.
- [x] Phase 3: Fix continuation to use true final state or fail loudly when unavailable.
- [x] Phase 4: Verify restart from a real production stage generated by the workflow.

### Known Errors / Blockers
- Existing production files created without `--record-momentum` cannot be exact restart sources; continuation now fails before MD and tells the user to regenerate the previous production stage with the current workflow.
- Existing production files created before the validated final-state marker also cannot be accepted as exact restart sources, even if they contain `output/mom`.

### Review
- Fixed:
  - production `7.x -> 7.x` continuation now uses `production_restart` handoff, which copies saved final coordinates and momentum exactly and does not refresh hybrid carriers;
  - continuation now requires both a validated final-state marker on `output/mom` and `input/mom.restart_valid=1` after handoff;
  - `upside` now appends a final MD output sample after the integration loop so the next stage restarts from the true final state;
  - workflow MD stages mark `output/mom.restart_final_state_valid=1` only after verifying the final output time matches the requested duration;
  - hybrid transition carry-forward now advances by the full final output time instead of the last periodic pre-final frame.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py` passed;
  - `bash -n` passed for all four MARTINI workflow shell entrypoints;
  - `cmake --build obj --target upside` passed with only existing warnings;
  - a patched 500-step production source wrote final `output/time=1.0` and matching `output/mom`;
  - workflow continuation from that source used `--restart-using-momentum`, kept Rg finite from `12.6 A` to `12.3 A`, and produced finite positions/momenta;
  - continuation input exactly matched source final state: position max diff `0.0`, momentum max diff `0.0`, `input/mom.restart_valid=1`, `sc_env_transition_step_start=500`;
  - workflow-generated `stage_7.1` was restarted into `stage_7.2`, used `--restart-using-momentum`, stayed stable from `Rg 12.2 A` to `12.1 A`, and produced `sc_env_transition_step_start=1000`;
  - `stage_7.2` input exactly matched `stage_7.1` final state: position max diff `0.0`, momentum max diff `0.0`;
  - continuation from a source with `output/mom` but no validated final-state marker failed before MD with the expected final restart-state error;
  - continuation from an older source without `output/mom` failed before MD with the expected restart-valid momentum error.

## 2026-04-24 dry-MARTINI Runtime Acceleration Import

### Project Goal
- Import the dry-MARTINI runtime acceleration logic from `/Users/yinhan/Documents/upside2-md_temp/src/martini.cpp` into this repo's `src/martini.cpp`.
- Keep the scope limited to performance paths used by the active dry-MARTINI runtime, not unrelated feature changes from the temp checkout.

### Architecture & Key Decisions
- Only transplant acceleration mechanisms that directly reduce hot-loop work in the current runtime:
  - cached Verlet-style pairlists for `MartiniPotential`,
  - cached per-pair parameter metadata so the main nonbonded loop avoids repeated coefficient unpacking / spline-map lookup,
  - cached active contact lists for `martini_sc_table_1body`, which is the active production SC/environment node in this checkout.
- Do not import unrelated temp-file changes such as fix-rigid extensions, stage-parameter cleanup, or removal of legacy compatibility nodes unless they are required for the acceleration paths to compile.
- Preserve backward compatibility in `src/martini.cpp`:
  - keep existing node names and interfaces intact,
  - keep legacy `martini_sc_table_potential` available.

### Execution Phases
- [x] Confirm the exact acceleration deltas in the temp file and map them onto the active runtime surfaces in this checkout.
- [x] Patch `src/martini.cpp` with the selected acceleration paths only.
- [x] Build and verify the modified MARTINI runtime.
- [x] Document the imported optimizations and verification results.

### Known Errors / Blockers
- The checked-in `obj/` build cache is stale from the old moved repo path (`/Users/yinhan/Documents/upside2-md-martini`), so verification in this pass used a fresh out-of-tree build directory under `/tmp` instead of the existing `obj/`.

### Review
- Confirmed the temp checkout acceleration scope:
  - `MartiniPotential` adds a cached Verlet-style pairlist plus cached per-pair parameter metadata/pointer indirection;
  - the active SC/environment runtime in this checkout is `martini_sc_table_1body`, and the temp checkout accelerates that path with a cached active-contact list.
- Implemented in `src/martini.cpp`:
  - `MartiniPotential` now compacts pair coefficients into a unique parameter table, caches direct spline pointers per unique parameter row, and rebuilds an active pairlist only when atoms move beyond half the configured skin or the box changes;
  - `MartiniScTableOneBody` now caches row/environment contacts within `cutoff + cache_buffer` and reuses them across both value and derivative passes until the cache is invalidated by motion or box changes;
  - the NPT box updater now also updates `martini_sc_table_1body`, so the active-contact cache remains valid under box scaling for the live production node.
- Verification:
  - configured a fresh build tree with `cmake -S src -B /tmp/upside2-md-build-20260424`;
  - built `upside`, `upside_calculation`, and `upside_engine` successfully with `cmake --build /tmp/upside2-md-build-20260424 --target upside upside_calculation upside_engine`.
- Residual warnings:
  - build still reports a pre-existing `%p` formatting warning in `martini_masses::get_mass(...)` and unrelated existing warnings in other translation units; no new build errors remain from the acceleration import.

## 2026-04-14 MARTINI Workflow Python Runtime Repair

### Project Goal
- Make `example/16.MARTINI/run_sim_1rkl.sh` use the repository `.venv` reliably so `py/martini_prepare_system_lib.py` can import `h5py`.
- Repair `install_python_env.sh` so a repo-local `.venv` still works after the repository has been moved on disk.

### Architecture & Key Decisions
- Treat the reported `ModuleNotFoundError: h5py` as an interpreter-selection bug, not a missing dependency bug:
  - the current `.venv/bin/python` imports `h5py` successfully;
  - the failure comes from `run_sim_1rkl.sh` resolving `python3` from Homebrew `python3.14` after sourcing a stale `.venv/bin/activate`.
- Fix the workflow wrapper directly by prepending the repo-local `.venv/bin` to `PATH` instead of depending on `activate` side effects.
- Fix the installer in place by repairing stale `.venv/bin/activate*` files and stale console-script shebangs under `.venv/bin` to the current repo path, avoiding an unnecessary full environment rebuild.

### Execution Phases
- [x] Confirm the interpreter/path mismatch and inspect the current `.venv` metadata.
- [x] Patch `example/16.MARTINI/run_sim_1rkl.sh` and `install_python_env.sh` with the minimal robust fix.
- [x] Verify the repaired runtime path and document the result.

### Known Errors / Blockers
- None so far. The local `.venv` already contains `h5py`, so verification does not depend on a fresh network install.

### Review
- Root cause verified:
  - `.venv/bin/python` imports `h5py`;
  - Homebrew `python3.14` does not;
  - the old `.venv/bin/activate` still pointed at `/Users/yinhan/Documents/upside2-md-martini/.venv`, so `run_sim_1rkl.sh` was picking the wrong interpreter.
- Implemented:
  - `example/16.MARTINI/run_sim_1rkl.sh` now prepends the repo-local `.venv/bin` to `PATH` instead of sourcing `activate`;
  - `install_python_env.sh` now repairs stale `activate*` files and stale `.venv/bin` Python-script shebangs in place.
- Verification:
  - `bash -n install_python_env.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -lc 'PROJECT_ROOT="$PWD"; source "$PROJECT_ROOT/source.sh"; export PATH="$PROJECT_ROOT/.venv/bin:$PATH"; python3 -c "import sys, h5py; print(sys.executable); print(h5py.__version__)"'`
  - `env PIP_NO_INDEX=1 bash install_python_env.sh`
  - `bash -lc 'source .venv/bin/activate && python3 -c "import sys, h5py; print(sys.executable); print(h5py.__version__)"'`
  - `bash -lc 'source .venv/bin/activate && pip --version'`

## 2026-04-14 Python Environment Installer Portability

### Project Goal
- Make `install_python_env.sh` succeed on both Linux and local macOS hosts while preserving the existing Python-package intent for UPSIDE2.

### Architecture & Key Decisions
- Keep Python `3.11` as the required interpreter version for the repository environment.
- Interpreter selection will stay configuration-driven:
  - honor `PYTHON_BIN` when provided;
  - otherwise discover a compatible Python `3.11` interpreter from common Linux/macOS locations instead of assuming a single command name.
- macOS setup will export Homebrew-backed dependency hints for packages that may need native HDF5 / compression libraries during pip installation.
- Core scientific packages remain required.
- Optional analysis packages remain best-effort so a platform-specific failure in `pyhdx`-related extras does not abort the whole environment bootstrap.

### Execution Phases
- [x] Audit the current installer and local Darwin environment.
- [x] Patch `install_python_env.sh` for cross-platform interpreter and package installation behavior.
- [x] Verify shell syntax and installer logic, then document the result.

### Known Errors / Blockers
- No active code blocker remains for this installer task.
- Verification completed in the current Darwin arm64 shell with `PIP_NO_INDEX=1`, which exercised interpreter discovery, macOS build-hint setup, venv reuse, and the post-install import checks without requiring fresh downloads.

## 2026-04-01 Stage-7 SC Table Integration

### Project Goal
- Consume the completed external `/Users/yinhan/Documents/SC-training` training results by building a native-unit `martini.h5`.
- Replace the legacy stage-7 sidechain back-mapping / rotamer production path in `example/16.MARTINI/run_sim_1rkl.sh` with a direct dry-MARTINI sidechain table interaction.
- Keep the new sidechain force field inactive for stages `6.0` through `6.6` and active only in production stage `7.0`.

### Architecture & Key Decisions
- `martini.h5` will be built from the assembled `sc_table.json` output and stored as a native-unit library (`nm`, `kJ/mol`) with a simple table schema:
  - residue order,
  - target dry-MARTINI type order,
  - uniform radial grid,
  - radial energy table.
- Simulation-side conversion into Upside units remains parameterized:
  - the stage-7 runtime node will read explicit `energy_conversion_kj_per_eup` and `length_conversion_angstrom_per_nm` attrs;
  - no numeric unit conversion will be baked into training outputs.
- Revised decision after user correction on 2026-04-02:
  - restore the reduced simulation-mass path (`ff_mass / 12`) in preparation and stage-file augmentation;
  - keep the weighted BB mapping revert in place;
  - retain the scientific-notation logging improvement only as a readability fix.
- Revised decision after user correction:
  - the stage-7 dry-MARTINI coupling must be integrated directly with Upside carriers, not evaluated on parallel protein MARTINI proxy coordinates and projected afterward.
- Sidechain/environment coupling will remain a separate stage-7 potential node, but it must be directly Upside-carried:
  - use deterministic Upside sidechain-bead coordinates from `affine_alignment -> placement_fixed_point_only_CB`;
  - act only between selected protein residues and non-protein dry-MARTINI particles;
  - add equal-and-opposite force to the dry particle and the Upside sidechain-bead coordinate so the backbone receives the feedback through the coord-node derivative chain;
  - all legacy protein MARTINI `SC` proxy nonbonded interactions in `martini_potential` must be disabled in stage 7.
- Backbone/environment coupling must also be directly Upside-carried:
  - treat the Upside `CA` carrier as the dry-MARTINI protein `BB` interaction site;
  - keep dry-MARTINI `BB` to surrounding non-protein dry particles active, but apply the resulting force directly to `CA` and the dry particle rather than to a proxy `BB` bead followed by projection;
  - revised after regression report: the old weighted `BB -> N/CA/C/O` active map has been restored for now; the attempted CA-only active BB remap has been reverted pending a deeper force-path audit;
  - revised after the 2026-04-02 production-drift audit:
    - direct `BB -> CA` and direct `CB` stage-7 forces must still honor the existing startup protein-feedback ramp (`sc_env_backbone_hold_steps`) instead of feeding full force back onto the protein from step `0`;
    - stage-7 must not evaluate any protein-internal MARTINI proxy-proxy terms, because those legacy bonded/nonbonded terms only add bookkeeping energy once the protein is carried directly by Upside;
  - do not evaluate protein internal SC-backbone dry-MARTINI interactions, because those are already handled on the Upside side.
- The workflow must remove the legacy production SC path together with its dead control surface:
  - no stage-7 `rotamer`, `placement_fixed_scalar`, or `placement_fixed_point_vector_only`;
  - no stage-7 `hybrid_sc_map` requirement;
  - no stage-7 SC-relaxation / proxy-control attributes that only served the old probabilistic proxy path.
- Revised decision after the 2026-04-02 cleanup pass:
  - remove the dead prep-only `hybrid_sc_map` export/validation path as well; the active workflow only needs `hybrid_control`, `hybrid_bb_map`, and `hybrid_env_topology`;
  - `src/box.cpp` and `src/box.h` were audited during this pass and left unchanged because they no longer contain the retired stage-7 proxy/rotamer subsystem, only active box/NPT infrastructure.
- `ALA` and `GLY` remain excluded from the new SC table because the trained library covers the current `18` non-empty canonical sidechain types; those residues still participate through the existing backbone/environment MARTINI path only.

### Execution Phases
- [x] Phase A: Build `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- [x] Phase B: Implement a new stage-7 SC table runtime node and wire its box-scaling behavior.
- [x] Phase C: Replace the production-side script injection in `run_sim_1rkl.sh` with stage-7 CB + SC-table setup and remove the legacy SC rotamer path.
- [x] Phase D: Update the production helper workflow to match the new stage-7 path.
- [x] Phase E: Remove the remaining parallel protein-MARTINI force path so stage 7 uses direct Upside-carried BB/SC interactions only.
- [x] Phase F: Verify `martini.h5` generation, build integrity, stage-7 HDF5 injection structure, and workflow/runtime behavior for the direct-Upside design.

### Known Errors / Blockers
- No code-path blocker remains for stage-7 SC-table injection in this environment:
  - the repaired repository `.venv` provides both `h5py` and `tables`;
  - the stage-7 injector now falls back to the protein ITP when fresh prepared files omit `/input/sequence`.
- Reduced-step full-workflow runs are still not a stability benchmark for the hybrid physics:
  - the intentionally shortened `6.4-6.6` relaxation settings used for fast verification can drive the pre-production MARTINI state to enormous energies before stage 7;
  - that does not block the stage-7 integration itself, but any physical assessment still needs the intended relaxation horizon / benchmark settings.

### Review
- 2026-04-02 legacy-code removal verification:
  - physically removed the dead probabilistic SC / rotamer / placement subsystem from `src/martini.cpp` instead of leaving it compiled behind runtime gates;
  - removed the dead prep-only `hybrid_sc_map` generation and validation path from `example/16.MARTINI/prepare_system_lib.py`, `example/16.MARTINI/prepare_system.py`, and `example/16.MARTINI/validate_hybrid_mapping.py`;
  - audited `src/box.cpp` and `src/box.h` and left them intact because the remaining code there is still-active barostat/box plumbing rather than disabled hybrid legacy logic;
  - verification:
    - `cmake --build obj` passed,
    - `python -m py_compile example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/prepare_system.py example/16.MARTINI/validate_hybrid_mapping.py` passed under `.venv`,
    - `bash -n example/16.MARTINI/run_sim_1rkl.sh` and `bash -n example/16.MARTINI/test_prod_run_sim_1rkl.sh` passed,
    - a fresh prep run in `/tmp/hybrid_mapping_bb_only` produced a mapping file whose `/input` group contains only `hybrid_bb_map`, `hybrid_control`, and `hybrid_env_topology`,
    - a shortened active-workflow run in `/tmp/legacy_cleanup_short` accepted that mapping schema and completed the full current `run_sim_1rkl.sh` ladder through stage `7.0`.
- 2026-04-02 production-drift fix verification:
  - restored the startup protein-feedback ramp on the direct stage-7 force paths in `src/martini.cpp`:
    - direct `BB`/environment MARTINI forces now scale only the protein-side feedback by `compute_sc_backbone_feedback_mix(...)`,
    - `martini_sc_table_potential` now does the same for `CB` feedback while leaving the environment-side force unchanged;
  - removed the remaining active protein-internal MARTINI proxy-proxy path during stage 7 by making `allow_protein_pair_by_rule(...)` reject all protein-proxy pairs in the active hybrid runtime;
  - targeted replay verification from the same saved production handoff (`example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`) now stays negative through the user’s failure window:
    - archived current run at step `2050`: `martini_potential 26387.78`, `total 26216.96`;
    - patched replay at step `2050`: `martini_potential -25985.45`, `total -26147.28`.
- Implemented `example/16.MARTINI/build_sc_martini_h5.py` and generated `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- Implemented `example/16.MARTINI/inject_sc_table_stage7.py` so stage 7 now:
  - preserves the existing Upside backbone-node augmentation,
  - deletes the legacy SC production path,
  - injects deterministic CB coordinates and the new `martini_sc_table_potential`.
- Updated both `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` to use the new stage-7-only SC table path and auto-build the library if only the assembled training JSON is available.
- Implemented the new `martini_sc_table_potential` runtime in `src/martini.cpp`:
  - CB-to-env radial table evaluation,
  - equal-and-opposite force accumulation,
  - additive backbone feedback through the CB coord node,
  - participation in NPT box scaling.
- Removed the remaining parallel protein proxy path from the active stage-7 runtime/workflow:
  - `martini_potential` now skips protein `SC` proxy MARTINI pairs entirely during active hybrid production;
  - protein `BB` proxy MARTINI pairs now use the Upside `CA` carrier coordinate as the interaction site and accumulate force directly onto `CA` plus the dry-MARTINI partner;
  - the active `hybrid_bb_map` is now CA-only in both exported mapping files and injected runtime stage files (`atom_mask=[0,1,0,0]`, `weights=[0,1,0,0]`), and runtime parsing now rejects any remaining active `N/C/O` carrier rows;
  - stage-file injection no longer copies `hybrid_sc_map` into runtime stages, the old shell-level SC control helper is gone from the active workflow path, and the legacy stage-7 rotamer augmenter now hard-fails if called.
- Verification completed in this environment:
  - `python3 -m py_compile` passed for the new Python scripts;
  - `bash -n` passed for both workflow shell scripts;
  - `cmake --build obj` passed;
  - `h5ls -r parameters/ff_2.1/martini.h5` confirmed the generated library shape.
- Regression revert verification completed after the user reported a broken minimization stage:
  - reverted the attempted native-mass switch and the CA-only active BB remap;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run again reaches stage `6.0` with a normal minimization startup:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - this replaces the broken regression signature reported by the user:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/30400945802214465601536.00/...`
- Native-mass-only restoration verification completed on 2026-04-02:
  - restored the native mass path while leaving the weighted 4-carrier BB map in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run still starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - fresh `1rkl.stage_6.0.up` in `/tmp/martini_mass_restore` again shows native dry-MARTINI mass values in `/input/mass` (`72.0` and `45.0` for the physical dry particles).
- Reduced-mass restoration verification completed after the user reported box blow-up with original masses:
  - restored reduced simulation masses (`ff_mass / 12`) in preparation and stage-file augmentation while leaving the scientific-notation logging change in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run in `/tmp/martini_mass_reduced` again starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - direct HDF5 inspection of `1rkl.stage_6.0.up` confirms the reduced mass values are back in `/input/mass`:
    - `4082` particles at `6.0`,
    - `8` particles at `3.75`.
- Logging/readability follow-up on 2026-04-02:
  - the huge `martini_potential` / `total` values seen in current stage-7 logs are real runtime values, not a `printf` type mismatch or summation bug;
  - `src/main.cpp` now formats large-magnitude energies in scientific notation for readability only.
- Additional end-to-end verification completed after the `.venv` repair:
  - the focused production helper `example/16.MARTINI/test_prod_run_sim_1rkl.sh` now runs end to end from a real `6.6 -> 7.0` handoff, injects the new SC table successfully, and completes a short stage-7 production replay;
  - direct HDF5 inspection confirmed stage gating:
    - `6.6` contains no `martini_sc_table_potential`, no `placement_fixed_point_only_CB`, no `hybrid_sc_map`, and keeps hybrid activation disabled;
    - `7.0.prepared` contains `martini_sc_table_potential` and `placement_fixed_point_only_CB`, remains free of legacy `rotamer` nodes, and activates the hybrid path only for production;
  - the real top-level workflow `example/16.MARTINI/run_sim_1rkl.sh` completes end to end under shortened verification settings with the direct-Upside stage-7 path.

## 2026-03-31 SC-training Workflow Addendum

### Project Goal
- Establish a new `SC-training/` workflow that derives dry-MARTINI sidechain-type to dry-MARTINI particle-type training data and assembled tables independently of the legacy SC back-mapping runtime path.
- Make the workflow runnable both locally and on Slurm.
- Use `example/16.MARTINI/run_sim_1rkl.sh` as the canonical benchmark workflow target once the trained tables are consumed by the hybrid runtime.

### Architecture & Key Decisions
- The new workflow will live under a dedicated top-level folder `SC-training/`.
- The workflow will treat the current SC back-mapping / rotamer-placement production path as legacy design context rather than a dependency.
- The replacement SC/dry-MARTINI interaction must be a two-way force term: each SC/dry pair contributes equal-and-opposite force to the dry particle and to the effective protein-side SC site, and the protein-side contribution must be projected back onto the Upside carriers and added to the existing backbone/sidechain forces before integration.
- Existing backbone/environment MARTINI coupling remains in scope and separate from the new SC table:
  - the effective protein backbone dry-MARTINI `BB`/`CA` anchor should still interact with surrounding non-protein dry-MARTINI particles;
  - the new SC table should not reintroduce generic protein SC-backbone MARTINI interactions, because protein internal SC-backbone coupling is already handled on the Upside side.
- Residue-sidechain dry-MARTINI bead definitions will be sourced from `example/16.MARTINI/martinize.py` using the `martini22` forcefield definitions already used by the example workflow.
- Dry-MARTINI pair parameters will be sourced from `example/16.MARTINI/ff_dry/dry_martini_v2.1*.itp`.
- The unit contract must be explicit and parameterized:
  - `SC-training/` artifacts remain in native dry-MARTINI units (`nm`, `kJ/mol`, `e`, dielectric `15`);
  - simulation-side conversion into Upside units must come from explicit parameters / attrs, not baked conversion numbers in training outputs or MARTINI runtime code;
  - the current MARTINI workflow expects `UPSIDE_MARTINI_ENERGY_CONVERSION` and `UPSIDE_MARTINI_LENGTH_CONVERSION` to be provided explicitly when generating simulation inputs.
- The initial training workflow will emphasize reproducible task generation, batch execution, and result assembly; any simplifying assumptions used for the first pass must be explicit in workflow outputs and README docs.
- Slurm support should follow the proven ConDiv workflow shape: explicit path resolution, environment setup, generated manifests/specs, array-friendly execution, and a submission wrapper rather than an ad hoc one-off `sbatch` command.

### Execution Phases
- [x] Phase A: Inspect and reuse the relevant MARTINI residue/type definitions and Slurm workflow patterns.
- [x] Phase B: Create `SC-training/` scripts for manifest generation, per-task training execution, and final result assembly.
- [x] Phase C: Add local and Slurm entrypoints that stage and run training jobs from the new folder.
- [x] Phase D: Add minimal documentation describing assumptions, outputs, and benchmark expectations.
- [x] Phase E: Verify script integrity and basic dry-run behavior.

### Known Errors / Blockers
- Runtime consumption of the new SC training table is not yet implemented in `src/martini.cpp`, so benchmark orchestration may need to stop at workflow staging/logging until that integration exists.
- That runtime integration still needs to enforce the new mechanical requirements explicitly:
  - SC/dry interactions must remain equal-and-opposite;
  - SC feedback must be added to existing Upside carrier gradients rather than replacing them;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling must remain active while protein internal SC-backbone coupling stays excluded from the new table.
- The first-pass training workflow may need an explicit simplifying geometry assumption for residue-level sidechain reduction; that assumption must be documented, not hidden.
- `SC-training/` is now self-contained for training and Slurm staging, but the optional benchmark hook still targets `example/16.MARTINI/run_sim_1rkl.sh` and therefore still requires the full repository when benchmark execution is requested.

### Review
- Implemented a new `SC-training/` workflow with:
  - manifest generation from `example/16.MARTINI/martinize.py` (`martini22`) and `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`,
  - per-task residue/target training jobs,
  - assembled output tables,
  - local runner,
  - Slurm staging/submission wrappers,
  - benchmark staging against `example/16.MARTINI/run_sim_1rkl.sh`.
- Smoke verification succeeded in `/tmp/sc_training_smoke`:
  - manifest generation completed,
  - full local run completed,
  - assembled summary reported `684` tasks across `18` residues and `38` target particle types,
  - Slurm array, collector, and benchmark scripts were generated successfully with `--no-submit`.
- Follow-up workflow correction restored the default target set to the full dry-MARTINI atomtype list from the bundled forcefield after a temporary over-narrowing to non-ring types:
  - default table cardinality is `18 x 38 = 684` residue-target tasks;
  - manifest/task outputs now keep the full target set while retaining the explicit spherical surrounding-position samples around the sidechain center.
- Follow-up architecture review confirmed additional force-path requirements for the future runtime hookup:
  - the new SC/dry term must remain a two-way force and its protein-side feedback must be added on top of existing Upside forces;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling remains required and should stay distinct from the new SC table;
  - dry-MARTINI unit handling is now defined as native-unit training plus parameterized simulation-side conversion, not a baked training/runtime numeric contract.
- Portability follow-up verification succeeded from a standalone folder copy under `/private/tmp`:
  - `SC-training/data/` now bundles the required dry-MARTINI nonbond parameter file and martini22 sidechain residue definitions;
  - `workflow.py init-run`, `workflow.py run-local`, `run_local.sh`, and stage-only `submit_remote_round.sh` all work from the copied folder without the parent repository;
  - generated Slurm array/collector scripts now reference the copied `SC-training/` folder directly instead of requiring `source.sh` or repo-root paths for training.
- Slurm-array verification succeeded from a fresh staged run:
  - `submit_remote_round.sh` now reaches `workflow.py submit-slurm` reliably even when repo-root `source.sh` would otherwise trip over an unset `MY_PYTHON`;
  - generated training script uses `#SBATCH --array=0-683` for the current manifest and dispatches one residue-target task per array index via `run-array-task --task-id "$SLURM_ARRAY_TASK_ID"`;
  - staged round manifest count and training manifest count matched exactly (`684` / `684`), with a separate collector script for post-array assembly.
- Batch-spooled launcher verification succeeded after the user-reported Midway failure:
  - `submit_remote_round.sh` and `run_local.sh` no longer trust `dirname "$0"` as the workflow root when executed from a Slurm spool copy;
  - the wrappers now locate the real `SC-training` directory by checking candidate directories for `workflow.py`, including `SLURM_SUBMIT_DIR` and `SLURM_SUBMIT_DIR/SC-training`;
  - simulated `sbatch submit_remote_round.sh` from inside `SC-training/` and simulated `sbatch SC-training/submit_remote_round.sh` from the repo root both staged cleanly and produced matching `684`-task manifests.
- Unit-contract follow-up verification succeeded:
  - training manifests now record native-unit policy only and no longer emit baked Upside conversion numbers;
  - `prepare_system_lib.py` now requires explicit conversion env vars and writes those MARTINI unit-conversion attrs into the simulation input;
  - `src/martini.cpp` now derives Coulomb scaling from those attrs instead of a literal hardcoded runtime constant.

### Phase 1: Geometric Mapping (Defining the Training Ground)
Before you can train the splines, you must define the physical relationship between the two models.
1. **Define the Upside Virtual Bead:** For each of the 20 standard amino acids, identify the Upside definition for the virtual side-chain bead position ($\vec{r}_{SC}$) and its normal vector ($\vec{n}$) relative to the backbone.
2. **Define the MARTINI Target:** For each amino acid, map the corresponding multi-bead Dry-MARTINI representation (e.g., mapping a Tryptophan side chain to its 4 constituent MARTINI beads). 

### Phase 2: Generating the 2D Energy Landscapes
You need to calculate the "ground truth" thermodynamic interactions using the Dry-MARTINI force field. Because Dry-MARTINI is an implicit solvent model, this can be done in a vacuum.
1. **Set Up the Grid:** Place your multi-bead MARTINI side chain at the origin. Place a single target Dry-MARTINI particle (e.g., a standard $C1$ lipid bead or $P4$ polar bead) in the space around it.
2. **Scan the Space:** Systematically move the single target bead across a 2D grid defined by:
   * **Distance ($r_{12}$):** The distance from the target bead to the Upside virtual center ($\vec{r}_{SC}$).
   * **Angle ($\theta_1$):** The angle between the Upside normal vector ($\vec{n}$) and the vector connecting the two beads ($\vec{n}_{12}$).
3. **Calculate Energy:** At every grid point, calculate the sum of the standard Dry-MARTINI Lennard-Jones and Coulomb interactions between the single target bead and the multi-bead side chain.
4. **Repeat:** Do this for all 20 amino acids against all relevant Dry-MARTINI particle types.

### Phase 3: Mathematical Spline Fitting
You must compress those 2D energy landscapes into the modified, factorized Upside equation.
1. **Set Up the Equation:** Use the simplified hybrid Hamiltonian that drops the second angle term:
   $$V_{hybrid} = V_{radial}(r_{12}) + ang_1(-n_1 \cdot n_{12}) \cdot V_{angular}(r_{12})$$
2. **Initialize Splines:** Set up $V_{radial}$, $ang_1$, and $V_{angular}$ as 1D natural cubic splines with adjustable control points (knots).
3. **Optimize:** Use a least-squares fitting algorithm or gradient descent to adjust the spline knots until the output of the $V_{hybrid}$ equation closely matches the 2D energy grid you generated in Phase 2.

### Phase 4: Engine Integration (C++ Implementation)
Append the new physics into the Upside MD engine's core loop.
1. **Introduce Real Particles:** Modify the Upside state to accept independent, mass-bearing Dry-MARTINI particles with standard $x, y, z$ coordinates and momenta.
2. **Calculate the Force:** In the force-evaluation loop, evaluate $V_{hybrid}$ and calculate its spatial gradient (force) and angular derivative (torque).
3. **Push the Real Bead:** Apply the negative spatial gradient directly to the Dry-MARTINI particle's momentum.
4. **Push the Ghost Bead:** Apply the exact inverse spatial vector to the Upside virtual side chain ($\vec{r}_{SC}$). Apply the angular derivative as a torque to the virtual normal vector ($\vec{n}$).
5. **Distribute to the Backbone:** Feed the virtual force and torque into Upside's existing chain-rule matrices to distribute the physical forces onto the $N$, $C\alpha$, and $C$ backbone atoms.

### Phase 5: Validation and Thermodynamic Tuning
Run test simulations to ensure the hybrid physics do not break the model.
1. **Check Momentum Conservation:** Monitor the total momentum of the simulation box. If the sum of forces on the MARTINI beads and the Upside backbone atoms does not equal exactly zero at every time step, there is a bug in your chain-rule implementation.
2. **Tune the Scaling Factor ($\kappa$):** Introduce a global scaling multiplier to the hybrid potential: $V_{total} = \kappa(V_{hybrid})$. Run a folded protein in a box of your Dry-MARTINI particles. If the protein denatures because the MARTINI interactions are overwhelmingly strong compared to Upside's internal folding potentials, lower $\kappa$ until the thermodynamic balance is restored.

What specific biological system or environment are you aiming to simulate first once this hybrid engine is built?

## 2026-04-13 Interface-Only Hybrid Scaling Calibration

### Project Goal
- Calibrate one shared dry-MARTINI interaction scale from bilayer-only DOPC diffusion runs.
- Apply that scale only to Upside/dry-MARTINI interface interactions in active hybrid production.
- Leave bilayer-bilayer interactions unchanged in the hybrid workflow.

### Architecture & Key Decisions
- Use one shared scalar `protein_env_interface_scale` rather than separate LJ and Coulomb knobs.
- In active hybrid stage `7.0`, apply the scale only to cross-interface protein/environment interactions:
  - direct `martini_potential` protein-environment pairs,
  - `martini_sc_table_1body`,
  - legacy `martini_sc_table_potential` for backward compatibility.
- Preserve current hybrid startup logic such as the protein-feedback ramp; the new factor is an additional interface-strength control, not a replacement for startup stabilization.
- Calibrate the scalar in the bilayer-only workflow by rewriting dry-MARTINI pair coefficients before simulation:
  - `epsilon *= pair_scale`
  - `q_i *= sqrt(pair_scale)`
  - `q_j *= sqrt(pair_scale)`
- Keep bilayer calibration separate from hybrid runtime application:
  - bilayer-only runs may scale all MARTINI pairs to match target fluidity,
  - hybrid production must apply the chosen factor only to protein-environment interface terms.
- Correct the active-stage assumption from earlier notes:
  - the current production SC path in this checkout is `martini_sc_table_1body`,
  - `martini_sc_table_potential` remains only as a compatibility surface.

### Execution Phases
- [x] Phase A: Extend hybrid-control schema and runtime state with `protein_env_interface_scale`.
- [x] Phase B: Apply interface-only scaling in `src/martini.cpp` for direct pair terms and both SC-table paths.
- [x] Phase C: Expose the production control in `example/16.MARTINI/run_sim_1rkl.sh`.
- [x] Phase D: Add bilayer-only pair-scale calibration support plus a fixed sweep wrapper and cross-run report.
- [x] Phase E: Verify build/script integrity and record review notes, progress, and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the implementation must respect two different semantics:
  - bilayer calibration runs scale all MARTINI interactions inside the bilayer-only workflow,
  - hybrid production must scale only protein-environment interface terms.

### Review
- Implemented the new hybrid control attr `protein_env_interface_scale` in:
  - `src/martini.cpp`
  - `py/martini_prepare_system_lib.py`
  - `example/16.MARTINI/run_sim_1rkl.sh`
- Runtime behavior now matches the intended split:
  - active direct protein-environment `martini_potential` pairs scale only the cross-interface LJ+Coulomb interaction,
  - bilayer-bilayer and other non-interface pairs remain unchanged in the hybrid workflow,
  - `martini_sc_table_1body` and legacy `martini_sc_table_potential` both respond to the same scale.
- Added bilayer calibration support in:
  - `bilayer-lateral-diffusion/workflow.py`
  - `bilayer-lateral-diffusion/submit_interface_scale_calibration_round.sh`
  - `bilayer-lateral-diffusion/report_interface_scale_calibration.py`
- Verification completed:
  - `python3 -m py_compile py/martini_prepare_system_lib.py bilayer-lateral-diffusion/workflow.py bilayer-lateral-diffusion/report_interface_scale_calibration.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -n bilayer-lateral-diffusion/submit_remote_round.sh`
  - `bash -n bilayer-lateral-diffusion/run_local.sh`
  - `bash -n bilayer-lateral-diffusion/submit_interface_scale_calibration_round.sh`
  - `cmake --build obj --target upside`
  - no-submit wrapper staging for `/tmp/bilayer_interface_scale_calibration/pair_scale_0p85`
  - reduced local bilayer smoke runs at `pair_scale=1.0` and `0.85`
  - direct HDF5 inspection of staged `stage_6.0.prepared.up` files confirmed:
    - `epsilon` ratio = `0.85`
    - charge ratios = `sqrt(0.85) = 0.921954`
  - `report_interface_scale_calibration.py` executed successfully on a temporary two-run smoke tree.

## 2026-04-13 Root-Level Hybrid Interface Sweep Workflow

### Project Goal
- Add a new project-root workflow folder that sweeps `PROTEIN_ENV_INTERFACE_SCALE` over the real hybrid `1RKL` production path.
- Reuse `example/16.MARTINI/run_sim_1rkl.sh` as the single simulation entrypoint.
- Submit the sweep to Slurm using the same array/collector pattern already proven in `bilayer-lateral-diffusion/`.

### Architecture & Key Decisions
- Create a dedicated root-level folder `hybrid-interface-sweep/`.
- Keep the sweep workflow thin:
  - generate a manifest of `(interface_scale, replicate)` tasks,
  - run one full `run_sim_1rkl.sh` instance per task in its own task-local `RUN_DIR`,
  - collect only execution/status metadata rather than adding a new scientific analyzer in v1.
- Do not duplicate hybrid simulation logic in Python; invoke `example/16.MARTINI/run_sim_1rkl.sh` directly with environment overrides.
- Follow the `bilayer-lateral-diffusion` Slurm submission shape:
  - `workflow.py`
  - `run_local.sh`
  - `submit_remote_round.sh`
  - generated array and collector `.sbatch` scripts
- Preserve reproducibility by capturing a whitelist of relevant `run_sim_1rkl.sh` environment overrides into the manifest at `init-run` time.

### Execution Phases
- [x] Phase F: Add root and local tracker entries for the new hybrid sweep workflow.
- [x] Phase G: Create `hybrid-interface-sweep/workflow.py` with manifest, task runner, assembly, and Slurm submission.
- [x] Phase H: Add `run_local.sh`, `submit_remote_round.sh`, `README.md`, and local tracker files.
- [x] Phase I: Run static checks plus a reduced smoke run and no-submit Slurm staging.
- [x] Phase J: Record review notes and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the new workflow must reuse the real hybrid shell entrypoint rather than silently diverging from `example/16.MARTINI/run_sim_1rkl.sh`.

### Review
- Added a new root-level workflow folder:
  - `hybrid-interface-sweep/workflow.py`
  - `hybrid-interface-sweep/run_local.sh`
  - `hybrid-interface-sweep/submit_remote_round.sh`
  - `hybrid-interface-sweep/README.md`
  - local `plan.md`, `progress.md`, and `findings.md`
- The new workflow:
  - generates a manifest of `(interface_scale, replicate)` tasks,
  - runs the real `example/16.MARTINI/run_sim_1rkl.sh` entrypoint per task,
  - assigns each task its own task-local `RUN_DIR`,
  - captures a whitelist of relevant hybrid env overrides into the manifest,
  - stages Slurm array and collector scripts using the same pattern as `bilayer-lateral-diffusion`.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - no-submit Slurm staging under `/tmp/hybrid_interface_sweep_smoke`
  - reduced local smoke run with:
    - `interface_scale = 0.85`
    - `replicates = 1`
    - `MIN_60_MAX_ITER=1`
    - `MIN_61_MAX_ITER=1`
    - `EQ_62_NSTEPS ... EQ_66_NSTEPS = 1`
    - `PROD_70_NSTEPS = 1`
    - `EQ_FRAME_STEPS = 1`
    - `PROD_FRAME_STEPS = 1`
  - smoke result verification:
    - task result JSON recorded `success = true`
    - task-local `stage_7.0.up` exists at `/tmp/hybrid_interface_sweep_smoke/tasks/scale0p85_r01/run/checkpoints/1rkl.stage_7.0.up`
    - assembled summary recorded `1` successful task and `1` completed scale.

## 2026-04-13 Hybrid Sweep Post-Run Analysis

### Project Goal
- Add a post-run analysis path for `hybrid-interface-sweep/` that measures hybrid diffusion/fluidity signals from completed `stage_7.0.up` files.
- Make that analysis runnable locally and on Slurm.
- Produce downloadable per-task and assembled analysis artifacts for later offline review.

### Architecture & Key Decisions
- Reuse the sweep manifest/task tree and analyze existing `stage_7.0.up` outputs post hoc rather than mixing analysis into the simulation task runner.
- Measure two signals per completed sweep task:
  - protein lateral COM MSD relative to bilayer COM, using protein `CA` carrier atoms from `hybrid_bb_map`,
  - lipid `PO4` lateral MSD as a bilayer fluidity guardrail.
- Use the same basic MSD treatment as the bilayer workflow:
  - XY unwrapping,
  - bilayer COM drift removal,
  - burn-in removal,
  - linear `MSD = 4Dt + b` fit over a fixed lag window.
- Add the analysis as workflow subcommands plus a thin Slurm wrapper under `hybrid-interface-sweep/`, following the same array/collector pattern as the existing run workflow.

### Execution Phases
- [x] Phase K: Add root and local tracker entries for the analysis extension.
- [x] Phase L: Extend `hybrid-interface-sweep/workflow.py` with manifest discovery, per-task analysis, assembly, and Slurm submission for analysis.
- [x] Phase M: Add analysis wrapper/docs in `hybrid-interface-sweep/`.
- [x] Phase N: Run static checks plus reduced analysis smoke verification.
- [x] Phase O: Record review notes and lessons.

### Known Errors / Blockers
- No blocker is known yet, but the analysis must select protein carriers and lipid molecules from the actual stage-file schema used by `run_sim_1rkl.sh`, not from assumptions copied from other workflows.

### Review
- Extended `hybrid-interface-sweep/workflow.py` with a full post-run analysis path:
  - `init-analysis`
  - `run-analysis-local`
  - `run-analysis-task`
  - `assemble-analysis`
  - `submit-analysis-slurm`
- Added analysis wrappers and docs:
  - `hybrid-interface-sweep/run_analysis_local.sh`
  - `hybrid-interface-sweep/submit_analysis.sh`
  - updated `hybrid-interface-sweep/README.md`
- Analysis outputs now measure:
  - protein lateral COM diffusion relative to the bilayer COM,
  - lipid `PO4` lateral diffusion as a bilayer guardrail,
  - per-task MSD and linear-diffusion fit metadata in `analysis/results/tasks/*.json`,
  - assembled CSV and summary JSON outputs under `analysis/assembled/`.
- Verification completed:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced real sweep run in `/tmp/hybrid_interface_sweep_analysis_smoke` with `PROD_70_NSTEPS=40`
  - local analysis run on that smoke output
  - no-submit analysis Slurm staging under `/tmp/hybrid_interface_sweep_analysis_smoke/analysis/slurm`
  - generated artifacts verified:
    - `analysis/analysis_manifest.json`
    - `analysis/results/tasks/scale0p85_r01.json`
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/summary.json`

## 2026-04-29 C++ Core Cleanup For 1RKL Default Path

### Project Goal
- Refactor and clean `src/main.cpp`, `src/box.h`, `src/box.cpp`, and `src/martini.cpp` so only logic required by the default `example/16.MARTINI/run_sim_1rkl.sh` pipeline remains.

### Architecture & Key Decisions
- Scope is strict default-path execution only; optional workflow/debug branches are removed.
- Keep SC-env and BB-env interaction potentials active; do not disable or zero these force paths.
- Keep barostat/Ewald and hybrid runtime features only where needed by the active 1RKL pipeline.
- Use `/Users/yinhan/Documents/upside2-md-master/src/main.cpp` as reference for removing temporary CLI/debug deviations in `src/main.cpp`.

### Execution Phases
- [x] Phase 1: Clean `src/main.cpp` against master and keep only required workflow deltas.
- [x] Phase 2: Remove dead/debug scaffolding from `src/box.h` and `src/box.cpp`.
- [x] Phase 3: Remove dead wrappers/placeholders/debug branches from `src/martini.cpp` while preserving active physics.
- [x] Phase 4: Build and run reduced workflow verification for default 1RKL path.
- [x] Phase 5: Document results in `progress.md` and review section.

### Known Errors / Blockers
- None at start.

### Review
- Implemented cleanup in the four requested files:
  - `src/main.cpp`: removed split-potential/Rg logging scaffolding, removed temporary CLI branches (`nvtc`, `max-force`, `martini-hold-backbone`), kept only default-path flags required by `run_sim_1rkl.sh` (`duration-steps`, minimization controls, momentum restart controls).
  - `src/box.h` and `src/box.cpp`: removed barostat debug toggles/prints and equilibrium-freeze scaffolding; preserved active NPT scaling and Ewald reciprocal computation paths.
  - `src/martini.cpp`: removed dead placeholder and compatibility scaffolding (`martini_masses::martini_integration_cycle`, `ConjugateGradientMinimizer` node, `minimize_structure_with_regular_potential` stub), removed SC-energy-dump plumbing and nonprotein hard-sphere experimental branch, and removed debug/status prints.
- Physics integrity preserved:
  - SC-env and BB-env force paths remain active.
  - `martini_potential` and `martini_sc_table_1body` both remain present in generated production stage files.
- Verification:
  - `cmake --build obj --target upside` succeeded.
  - Reduced default-workflow run succeeded end-to-end:
    - `RUN_DIR=/tmp/cpp_cleanup_smoke` with minimized stage step counts.
    - Stage outputs through `1rkl.stage_7.0.up` and VTF extraction completed.
  - Stage-7 production file check confirmed:
    - `/input/potential/martini_potential` exists;
    - `/input/potential/martini_sc_table_1body` exists.

## 2026-04-29 Python Workflow Cleanup (Baseline + Explicit 7.x Continuation)

### Project Goal
- Refactor and clean `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, and `py/martini_extract_vtf.py` to remove dead/legacy/debug-hoarded logic.
- Keep only functionality required by `example/16.MARTINI/run_sim_1rkl.sh` fresh execution and explicit production continuation (`--continue-stage-70-from` -> `7.1`, `7.2`, ...).
- Preserve SC-env and BB-env interface interactions fully active.

### Architecture & Key Decisions
- Keep `run-hybrid-workflow` as the only intended high-level entrypoint in `martini_prepare_system.py`.
- Remove legacy command surfaces and wrapper flows not needed by the active pipeline.
- Keep continuation path only for explicit source file input; remove auto-discovery branches.
- Flatten defensive/legacy branching in VTF extraction to a direct path used by workflow outputs.

### Execution Phases
- [x] Phase 1: Prune `martini_prepare_system.py` CLI/command surfaces and continuation discovery branches.
- [x] Phase 2: Remove dead standalone/legacy conversion layers from `martini_prepare_system_lib.py` while preserving required functions.
- [x] Phase 3: Simplify `martini_extract_vtf.py` to workflow-required VTF behavior.
- [x] Phase 4: Run compile/syntax checks and validate retained continuation + physics-critical attrs remain intact.
- [x] Phase 5: Record review results and residual risks.

### Known Errors / Blockers
- None at task start.

### Review
- Implemented:
  - `py/martini_prepare_system.py` now supports only the `run-hybrid-workflow` command surface; legacy subcommands and fallback execution paths were removed.
  - `run-hybrid-workflow` parser now matches flags passed by `run_sim_1rkl.sh`; internal physics and stage defaults moved to explicit in-code constants.
  - Continuation now requires explicit `--continue-stage-70-from`; auto-discovery (`previous-run`/`auto-continue`) code paths were removed and now fail fast if used.
  - `py/martini_prepare_system_lib.py` removed dead standalone prep CLI/`main` orchestration and unused wrapper functions (`create_production_input`, `main_always_fixed`, martinize wrapper path).
  - `py/martini_extract_vtf.py` was simplified to VTF-only extraction and direct trajectory-group handling; defensive NaN frame substitution and unused PDB/output-group branches were removed.
- Physics integrity preserved:
  - production hybrid-control writes still enforce active SC-env configuration (`sc_env_*` attrs plus positive interface scale) and keep `production_nonprotein_hard_sphere=0`;
  - production stage injection still requires and validates active `martini_sc_table_1body` and `martini_potential`.
- Verification:
  - `.venv/bin/python -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py`
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --help`
  - `.venv/bin/python py/martini_extract_vtf.py --help`
  - `.venv/bin/python py/martini_prepare_system.py run-hybrid-workflow --previous-run-dir foo` (verified explicit-source continuation enforcement with clear error).
