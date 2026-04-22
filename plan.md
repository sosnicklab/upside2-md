
## 2026-04-22 MARTINI Python Legacy-Style Cleanup

### Project Goal
- Restyle `py/martini_extract_vtf.py`, `py/martini_prepare_system_lib.py`, and `py/martini_prepare_system.py` so they read like the older `py/` scripts in this repo and stop sounding AI-written, while preserving the active Example 16 workflow interface.

### Architecture & Key Decisions
- Preserve the current CLI and import surface used by `example/16.MARTINI/run_sim_1rkl.sh`.
- Rewrite toward the existing repo-local Python idiom:
  - plain string / `os.path` path handling instead of `Path`,
  - no type hints or `from __future__`,
  - flatter helper layout and simpler command dispatch,
  - shorter comments and direct error/status messages.
- Keep behavior stable unless a change is required to preserve current semantics after the style cleanup.

### Execution Phases
- [x] Update the root task-tracking files for this cleanup and record the implementation/verification path.
- [x] Rewrite `py/martini_prepare_system.py` toward the older repo script style without changing its command names or flags.
- [x] Rewrite `py/martini_extract_vtf.py` toward the older repo script style without changing its current CLI contract.
- [x] Rewrite `py/martini_prepare_system_lib.py` toward the older repo script style while preserving helper names and call signatures used by the workflow.
- [x] Run syntax checks, CLI help smoke checks, and a reduced real workflow verification.
- [x] Record the review notes, verification results, and any style-specific findings.

### Known Errors / Blockers
- None so far.

### Review
- Rewrote [py/martini_prepare_system.py](/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py:1) into the flatter repo style:
  - plain `os.path` handling,
  - no `Path` objects or type hints,
  - direct command dispatch instead of a command-handler table,
  - unchanged command names and flags.
- Reworked [py/martini_extract_vtf.py](/Users/yinhan/Documents/upside2-md/py/martini_extract_vtf.py:1) to match the rest of `py/` more closely:
  - string-based path handling,
  - simpler extractor status output,
  - reduced comment noise,
  - same positional args plus `--mode` and `--split-segments`.
- Simplified [py/martini_prepare_system_lib.py](/Users/yinhan/Documents/upside2-md/py/martini_prepare_system_lib.py:1) without changing the active MARTINI workflow surface:
  - removed `from __future__`, type annotations, and `Path`-typed interfaces,
  - introduced string-based path helpers,
  - shortened overly explanatory diagnostics and banner logging,
  - preserved the helper names used by `martini_prepare_system.py` and the Example 16 wrapper.
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py`
  - `bash -lc 'source .venv/bin/activate && source source.sh && python3 py/martini_prepare_system.py -h'`
  - `bash -lc 'source .venv/bin/activate && source source.sh && python3 py/martini_prepare_system.py build-sc-martini-h5 -h'`
  - `bash -lc 'source .venv/bin/activate && source source.sh && python3 py/martini_prepare_system.py inject-stage7-sc -h'`
  - `bash -lc 'source .venv/bin/activate && source source.sh && python3 py/martini_prepare_system.py set-initial-position -h'`
  - `bash -lc 'source .venv/bin/activate && source source.sh && python3 py/martini_extract_vtf.py -h'`
  - reduced real workflow run:
    - `bash -lc 'source .venv/bin/activate && source source.sh && RUN_DIR=/tmp/martini_legacy_style_verify MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1 EQ_62_NSTEPS=1 EQ_63_NSTEPS=1 EQ_64_NSTEPS=1 EQ_65_NSTEPS=1 EQ_66_NSTEPS=1 PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl.sh'`
- Observed result:
  - the help surfaces remained available under the repo environment;
  - the reduced workflow completed through fresh stage `7.0`;
  - fresh VTF extraction still ran for stages `6.0` through `7.0` on the rewritten script path.

## 2026-04-22 True Stage-to-Stage Handoff

### Project Goal
- Make every stage in `example/16.MARTINI/run_sim_1rkl.sh` start from the previous stage result, not from a freshly regenerated stage template with only coordinates copied forward.

### Architecture & Key Decisions
- Keep fresh generation only for `6.0`, since it has no predecessor.
- For `6.1` through `7.0`, derive each prepared file from the previous completed stage file:
  - copy the previous stage file,
  - move the last output coordinates into `/input/pos`,
  - remove stale `/output` and `/output_previous_*` groups,
  - patch only the stage control deltas (`current_stage`, barostat attrs, lipid-head restraint node, rigid-hold selection).
- Create the lipid-head restraint node during promotion from the carried-forward coordinates instead of relying on a fresh stage template.
- Preserve the current shell-visible checkpoint pattern (`*.prepared.up` then runnable `*.up`) while changing the source of truth for stage initialization.

### Execution Phases
- [x] Audit the actual differences between adjacent prepared stage files and determine whether the stage graph is otherwise identical.
- [x] Patch the workflow so every stage after `6.0` is promoted from the previous completed stage rather than regenerated from scratch.
- [x] Verify the reduced default workflow and confirm the promoted prepared files match the previous stage outputs.

### Known Errors / Blockers
- None so far.

### Review
- Adjacent stage audit on fresh reduced artifacts showed that the hybrid/runtime graph stayed constant across the ladder:
  - `atom_names`, `atom_roles`, `sequence`, `hybrid_bb_map`, `hybrid_env_topology`, and `martini_sc_table_1body` matched between neighboring stages;
  - the real per-stage deltas were only:
    - last-step coordinates / box,
    - `stage_parameters.current_stage`,
    - barostat targets,
    - presence/value of `restraint_position`,
    - `fix_rigid` selection state.
- Patched [example/16.MARTINI/run_sim_1rkl.sh](/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_1rkl.sh:270) to add generic `promote_stage_from_previous(...)` and use it for:
  - `6.1`
  - `6.2`
  - `6.3`
  - `6.4`
  - `6.5`
  - `6.6`
  - `7.0`
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - reduced fresh workflow run:
    - `source .venv/bin/activate && source source.sh && RUN_DIR=/tmp/cleanup_1rkl_true_stage_handoff MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1 EQ_62_NSTEPS=1 EQ_63_NSTEPS=1 EQ_64_NSTEPS=1 EQ_65_NSTEPS=1 EQ_66_NSTEPS=1 PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl.sh`
  - direct HDF5 promotion checks:
    - `6.1.prepared /input/pos` matches `6.0 /output/pos[-1]`,
    - `6.2.prepared /input/pos` matches `6.1 /output/pos[-1]`,
    - `7.0.prepared /input/pos` matches `6.6 /output/pos[-1]`,
    - promoted prepared files no longer contain stale `/output` groups.

## 2026-04-22 Production Promotion Cleanup

### Project Goal
- Remove redundant stage-7 reinjection from the default `example/16.MARTINI/run_sim_1rkl.sh` path by promoting the existing `6.6` prepared file into production instead of regenerating and reinjecting a fresh `7.0` file.

### Architecture & Key Decisions
- Keep the existing AA-backbone/hybrid-node injection for stages prepared from scratch (`6.0` to `6.6`), because those files still originate from `convert_stage(...)`.
- Stop rebuilding production from scratch:
  - derive `7.0.prepared` from `6.6.prepared`,
  - switch `current_stage` to `production`,
  - remove preproduction-only controls (`restraint_position`, `fix_rigid`, `barostat`) before the normal `6.6 -> 7.0` coordinate handoff.
- Preserve the existing optional production controls:
  - if `PROD_70_NPT_ENABLE=1`, recreate `/input/barostat` from the current environment settings;
  - if `PROD_70_BACKBONE_FIX_RIGID_ENABLE=1`, continue applying the production `fix_rigid` selection after promotion.

### Execution Phases
- [x] Confirm the exact schema differences between fresh `6.6.prepared` and `7.0.prepared`.
- [x] Patch the workflow to promote `6.6.prepared` into `7.0.prepared` instead of regenerating/injecting production from scratch.
- [x] Verify the default reduced workflow still reaches `7.0` and that stage-7 preparation no longer depends on a fresh injection pass.

### Known Errors / Blockers
- None so far.

### Review
- Confirmed on the fresh reduced baseline artifacts that `6.6.prepared` and `7.0.prepared` already shared the same hybrid payload:
  - identical `atom_names`, `atom_roles`, `sequence`, `hybrid_bb_map`, `hybrid_env_topology`, `placement_fixed_point_vector_only_CB`, and `martini_sc_table_1body` datasets;
  - the real differences were only:
    - `stage_parameters.current_stage`,
    - presence of `/input/potential/restraint_position`,
    - presence of `/input/barostat`,
    - presence of `/input/fix_rigid`.
- Patched `example/16.MARTINI/run_sim_1rkl.sh` to:
  - add `promote_preproduction_to_production(...)`,
  - build `7.0.prepared` from `6.6.prepared`,
  - remove the preproduction-only controls during promotion,
  - recreate production barostat only when `PROD_70_NPT_ENABLE=1`,
  - keep the optional production `fix_rigid` override after the normal `6.6 -> 7.0` handoff.
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - reduced fresh workflow run:
    - `source .venv/bin/activate && source source.sh && RUN_DIR=/tmp/cleanup_1rkl_production_promotion MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1 EQ_62_NSTEPS=1 EQ_63_NSTEPS=1 EQ_64_NSTEPS=1 EQ_65_NSTEPS=1 EQ_66_NSTEPS=1 PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl.sh`
- Observed result:
  - the reduced run completed through `7.0`;
  - the workflow log shows the last fresh injection at `6.6.prepared`, followed by direct production handoff without a new `inject-stage7-sc` pass;
  - fresh promoted `7.0.prepared` has `current_stage=production` and no `barostat`, `fix_rigid`, or `restraint_position`, while retaining the same hybrid payload datasets as `6.6.prepared`.

## 2026-04-21 1RKL Python Workflow Cleanup

### Project Goal
- Remove dead prep/extraction code from `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, and `py/martini_extract_vtf.py`, keeping only the fresh default `example/16.MARTINI/run_sim_1rkl.sh` path.

### Architecture & Key Decisions
- Keep only the mixed AA-backbone + DOPC + ion preparation path used by the default wrapper.
- Keep only the Python helper commands still needed by the default wrapper:
  - `build-sc-martini-h5`
  - `inject-stage7-sc`
  - `set-initial-position`
- Inline hybrid-mapping validation into mixed-system preparation instead of keeping a separate public subcommand.
- Keep VTF extraction modes `1` and `2` plus split-segment handling, but drop the extra generic CLI/output surface.

### Execution Phases
- [x] Audit the current Python entrypoints and helper surface against the default `run_sim_1rkl.sh` path.
- [x] Patch `py/martini_prepare_system.py` to remove unused prep modes, deprecated flags, duplicate wrapper logic, and dead subcommands.
- [x] Patch `py/martini_prepare_system_lib.py` to remove unused helpers and default-flow-inactive branches from stage conversion.
- [x] Patch `py/martini_extract_vtf.py` and `example/16.MARTINI/run_sim_1rkl.sh` to match the narrowed VTF/prep interface.
- [x] Verify the default workflow end to end on a reduced fresh run.

### Known Errors / Blockers
- None so far.

### Review
- Static verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_extract_vtf.py`
  - repo sweep confirmed removal of the deleted prep/extractor options and dead helper names.
- Fresh workflow verification:
  - ran a reduced fresh `example/16.MARTINI/run_sim_1rkl.sh` in `/tmp/cleanup_1rkl_python_verify` with one-step minimization/equilibration/production settings;
  - the run completed through stage `7.0`.
- Observed behavior on the narrowed Python path:
  - mixed-system preparation still wrote the runtime PDB and backbone metadata HDF5;
  - stage conversion still produced fresh `6.0` to `7.0` checkpoints;
  - VTF extraction still produced mode-`1` files for `6.x` and mode-`2` file for `7.0`;
  - stage-7 SC injection still produced `martini_sc_table_1body` on the fresh production checkpoint.

## 2026-04-21 1RKL-Only Runtime Cleanup

### Project Goal
- Remove runtime and generator code that is not needed by a freshly generated default `example/16.MARTINI/run_sim_1rkl.sh` workflow.
- Restrict compatibility to the current `1rkl` generation path rather than preserving old checked-in `.up` artifacts.

### Architecture & Key Decisions
- Keep the current stage-7 rotamer-weighted path because the active injector still emits:
  - `rotamer`
  - `placement_fixed_point_vector_only`
  - `placement_fixed_scalar`
  - `martini_sc_table_1body`
- Remove the unused radial-only `martini_sc_table_potential` runtime because the current injector does not emit it.
- Keep active Berendsen NPT, Ewald, fix-rigid, hybrid state, and current spring/restraint nodes.
- Remove dead stage-parameter force-constant plumbing while preserving only stage-label tracking (`current_stage`) for hybrid activation/minimization handoff.

### Execution Phases
- [x] Audit the fresh `1rkl` generation/runtime path and separate active nodes from dead compatibility code.
- [x] Patch `src/martini.cpp`, `src/main.cpp`, `src/main.h`, `src/box.cpp`, and `src/box.h` to remove dead hooks and legacy box/barostat fallbacks.
- [x] Patch `py/martini_prepare_system_lib.py` to stop writing dead stage-parameter payloads and the legacy barostat `compressibility` attr.
- [x] Rebuild and run targeted verification for the current `1rkl` path.
- [x] Record verification results and final review notes.

### Known Errors / Blockers
- None so far.

### Review
- Build verification:
  - `source .venv/bin/activate && source source.sh && cmake --build obj`
  - rebuild passed; only pre-existing warnings remained.
- Generator verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py`
  - `python3 -m py_compile py/martini_prepare_system_lib.py`
- Fresh workflow verification:
  - ran a reduced fresh `example/16.MARTINI/run_sim_1rkl.sh` in `/tmp/cleanup_1rkl_verify` with one-step minimization/equilibration/production settings;
  - the run completed through fresh stage `7.0`.
- Fresh stage-file schema verification:
  - `/input/stage_parameters` now contains only attrs `enable` and `current_stage` plus PyTables metadata; no dead bond/angle child groups remain;
  - `/input/barostat` now contains only `compressibility_xy` and `compressibility_z`; the legacy scalar `compressibility` attr is absent;
  - fresh `/input/potential` still contains the required active nodes for the current workflow, including:
    - pre-production `martini_potential`, `dist_spring`, `angle_spring`, `restraint_position`, `rotamer`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, `placement_fixed_point_vector_only_CB`, `martini_sc_table_1body`;
    - production `martini_potential`, `dist_spring`, `angle_spring`, `rotamer`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, `placement_fixed_point_vector_only_CB`, `martini_sc_table_1body`;
  - fresh stage files do not contain `martini_sc_table_potential`.
- Repo sweep verification:
  - repo-wide search found no remaining code references to `wrap_positions`, `martini_integration_cycle`, `clear_masses_for_engine`, `apply_stage_bond_params`, `apply_stage_angle_params`, `clear_stage_params_for_engine`, or `martini_sc_table_potential`.

## 2026-04-21 Example 16 Dead Path Removal Audit

### Project Goal
- Determine whether the Parrinello-Rahman barostat path in `src/box.cpp` and the Conjugate Gradient minimizer path in `src/martini.cpp` are actually used by `example/16.MARTINI/run_sim_1rkl.sh`.
- If either path is not used by that workflow, remove it completely together with any stale Example 16 configuration surface that only exists to select it.

### Architecture & Key Decisions
- Audit the real workflow path from `example/16.MARTINI/run_sim_1rkl.sh` through `src/main.cpp` before editing any runtime code.
- Treat "used by the workflow" as the code path that can be selected by the current Example 16 wrapper and its default stage preparation/invocation logic.
- Remove dead selection/configuration knobs together with dead implementations so the workflow surface matches the surviving runtime behavior.

### Execution Phases
- [x] Trace the Example 16 workflow to the runtime entry points for NPT and minimization.
- [x] Confirm whether the Parrinello-Rahman barostat branch and CG minimizer branch are actually reachable from that workflow.
- [x] Remove the unused runtime/code-path pieces and stale workflow configuration knobs.
- [x] Verify script syntax and rebuild the runtime.

### Known Errors / Blockers
- None so far.

### Review
- Confirmed workflow reachability:
  - `example/16.MARTINI/run_sim_1rkl.sh` minimization uses `--minimize`, which `src/main.cpp` routes to `martini_run_minimization(...)`;
  - the `ConjugateGradientMinimizer` node and `minimize_structure_with_regular_potential(...)` helper in `src/martini.cpp` were dead.
  - Example 16 pre-production NPT stages write Berendsen barostat type `0`;
  - stage `7.0` keeps NPT disabled by default, so the default workflow does not use Parrinello-Rahman at all.
- Removed:
  - the dead CG minimizer node/helper code from `src/martini.cpp`;
  - the Parrinello-Rahman branch and related state/type plumbing from `src/box.cpp` and `src/box.h`;
  - the dead Example 16 `barostat_type` selection surface from:
    - `example/16.MARTINI/run_sim_1rkl.sh`
    - `example/16.MARTINI/run.py`
    - `py/martini_prepare_system_lib.py`
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `python3 -m py_compile example/16.MARTINI/run.py py/martini_prepare_system_lib.py`
  - `cmake --build obj`

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
