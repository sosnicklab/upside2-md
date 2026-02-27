# RBM Cross-Potential + Training/Data Prep Implementation Plan

## Project Goal
Implement an end-to-end RBM cross-potential workflow for MARTINI hybrid systems:
1. Add a new trainable C++ potential node for protein-environment cross interactions.
2. Add a training script that learns RBM cross weights from prepared `.up` inputs and writes an artifact.
3. Add a download/preparation script that fetches inputs and prepares training-ready systems.

## Scope For This Cycle (MVP)
1. Implement `martini_rbm_cross_potential` in `src/martini.cpp` with force/energy and trainable parameter API.
2. Register the node in MARTINI node registration paths.
3. Add a Python training utility under `example/16.MARTINI/` to train and export RBM parameters.
4. Add a Python download+preparation utility under `example/16.MARTINI/` to build training manifests.
5. Verify build/syntax and preserve default runtime behavior when RBM node is absent.
6. Add a batch relaxation runner for rigid dry-MARTINI preprocessing that outputs a training-ready relaxed manifest.

## Architecture & Key Decisions
1. Keep RBM node additive and opt-in via HDF5 node presence only; no behavior changes for existing configs.
2. Use Gaussian radial basis expansion over protein-env distances with class-conditioned weights.
3. Provide `get_param`, `set_param`, and `get_param_deriv` for optimizer compatibility.
4. Gate RBM evaluation to active AA-backbone explicit-lipid mode when hybrid state exists.
5. Keep preparation/training flow script-based in `example/16.MARTINI/` to avoid broad core changes.

## `src/` Minimal-Change Rule (Mandatory)
1. For files under `src/`, only touch logic required for the RBM cross-potential node and registration.
2. Compare touched `src` files against `/Users/yinhan/Documents/upside2-md-master/src`.
3. Confirm no functional changes to master-equivalent paths unless explicitly required for RBM activation.

## Execution Phases

### Phase 1: Schema + Planning
- [x] Define RBM node HDF5 schema (attrs + datasets) and parameter flattening order.
- [x] Define training artifact schema and input manifest format.
- [x] Define integration points with existing `prepare_system.py` workflow.

### Phase 2: C++ RBM Node
- [x] Implement `martini_rbm_cross_potential` in `src/martini.cpp`.
- [x] Implement force and potential computation with cutoff and radial basis features.
- [x] Implement parameter interfaces (`get_param`, `set_param`, `get_param_deriv`).
- [x] Add node registration in both static and explicit registrar paths.

### Phase 3: Training Script
- [x] Add `example/16.MARTINI/train_rbm_cross_potential.py`.
- [x] Load `.up` inputs and derive protein/env class assignments.
- [x] Train RBM cross weights (positive vs negative contrastive batches).
- [x] Save artifact (including weights and metadata) and optional `.up` injection.

### Phase 4: Download + Preparation Script
- [x] Add `example/16.MARTINI/download_and_prepare_training_inputs.py`.
- [x] Download source structures from OPM/RCSB fallback.
- [x] Invoke existing preparation pipeline per system and collect outputs.
- [x] Write reproducible manifest for training script consumption.

### Phase 5: Verification
- [x] Run targeted C++ build check.
- [x] Run Python syntax checks on new/modified scripts.
- [x] Dry-run training script help/argument parsing.
- [x] Diff touched `src/` files against master clone and confirm limited scope.

### Phase 6: Batch Relaxation Runner
- [x] Add `example/16.MARTINI/batch_relax_rigid_dry.py`.
- [x] Consume `training_manifest.json` and execute per-system rigid dry workflow.
- [x] Capture relaxed stage outputs (`stage_6.6.up` by default) into a new manifest for training.
- [x] Verify CLI help + py_compile and provide command examples.

### Phase 7: Port Bugfixes From `upside2-md-martini`
- [x] Diff commit range `dd111dbffd9ab4842860ea6628e4a966bda3a049 -> 12b944d1bc743ea1104757cb94d2b49277b8295a`.
- [x] Apply source-code bugfix hunks to this branch (`src/box.h`, `src/box.cpp`, `src/martini.cpp`).
- [x] Run build verification after port.

## Known Errors / Blockers
1. Network availability may limit end-to-end download testing inside sandbox.
2. No existing golden RBM artifact in repo for numerical regression.

## Review Criteria
1. Existing MARTINI behavior remains unchanged when RBM node is not present.
2. `src/` edits remain minimal and localized.
3. New scripts are runnable with clear CLI contracts and deterministic outputs.
4. Batch relaxation can resume/reuse completed runs and emits a complete relaxed manifest.

## Review (This Run)
1. Added `martini_rbm_cross_potential` node in `src/martini.cpp` with:
   - HDF5 schema loading for atom indices/classes, radial basis, and class-conditioned weights.
   - Gradient-based force/energy computation over protein-env pairs within cutoff.
   - Trainable parameter API (`get_param`, `set_param`, `get_param_deriv`).
   - Registration in both `RegisterNodeType` and explicit `MartiniNodeRegistrar`.
2. Added `example/16.MARTINI/train_rbm_cross_potential.py`:
   - Consumes `.up` files (direct list or manifest).
   - Builds protein/environment class vocabularies from existing input metadata.
   - Trains RBM cross weights with positive-vs-shuffled negative feature contrast.
   - Writes artifact HDF5 and optional in-place injection of `/input/potential/martini_rbm_cross_potential`.
3. Added `example/16.MARTINI/download_and_prepare_training_inputs.py`:
   - Reads PDB ID list.
   - Downloads AA PDB from OPM with RCSB fallback.
   - Calls `prepare_system.py` per ID and records output `.up` files.
   - Writes reproducible manifest (`systems` + `failed`) for training.
4. Verification evidence:
   - `cmake --build /tmp/upside2-md-build-rbm -j4` succeeded.
   - `python3 -m py_compile` passed for both new scripts.
   - `--help` dry-runs passed for both scripts.
   - Smoke training/injection test passed on synthetic `.up` (`/tmp/rbm_smoke.input.up` -> `/tmp/rbm_smoke_artifact.h5`).
5. Master-parity caveat for `src/`:
   - Only touched `src` file is `src/martini.cpp`.
   - Requested reference file `/Users/yinhan/Documents/upside2-md-master/src/martini.cpp` is absent, so direct file-level parity diff for this file is not possible.
6. Downloader robustness updates after real-list failures:
   - `example/16.MARTINI/download_and_prepare_training_inputs.py` now supports martinizable AA sanitization (`--aa-clean-mode largest_complete_chain`) that keeps the largest contiguous chain segment with complete residue heavy atoms before `martinize.py`.
   - Added residue-level heavy-atom completeness checks (including histidine variants and ILE `CD`/`CD1` alias handling) to avoid known `martinize.py` failures from missing side-chain/backbone atoms.
   - Added martinize cache validation (`top_was_generated_from`) so stale `.top`/`.itp` outputs are not reused when AA input changes.
7. Post-fix verification:
   - Re-ran problematic entries with the user command style and confirmed success:
     - `4ryo` -> `train-data/runs/4ryo/test.input.up`
     - `4dve` -> `train-data/runs/4dve/test.input.up`
   - Additional subset re-run (first 20 IDs from user list) completed with `Prepared systems: 20` and `Failed systems: 0`, including previously failing IDs (`4x5m`, `4p79`).
8. Batch relaxation runner implementation:
   - Added `example/16.MARTINI/batch_relax_rigid_dry.py` to:
     - read `training_manifest.json`,
     - run a standalone stage-6 relax workflow per system with downloaded/generated AA/CG/ITP inputs,
     - force `MARTINIZE_ENABLE=0` to avoid remartinization drift,
     - write relaxed manifest entries keyed by `stage_6.6.up`,
     - set manifest `up_file` to the relaxed stage path (keeping original as `prepared_up_file`) so `train_rbm_cross_potential.py` consumes relaxed inputs directly.
   - Supports resume via `--skip-existing 1` and per-system logs (`batch_relax.log`).
9. Batch runner verification:
   - `python3 -m py_compile example/16.MARTINI/batch_relax_rigid_dry.py` passed.
   - CLI help check passed.
   - Smoke run passed on one system (`1ors`) with short stage counts and produced:
     - `/Users/yinhan/Documents/upside2-md/train-data/relax_batch_smoke/1ors/checkpoints/1ors.stage_6.6.up`
     - `/Users/yinhan/Documents/upside2-md/train-data/relax_batch_smoke/relaxed_training_manifest.json`
   - Resume check passed: rerun reports `SKIPPED_EXISTING` for completed system.
10. Standalone stage-6 relax script per user correction:
   - Added `example/16.MARTINI/run_relax_6x_rigid_dry.sh` that does not depend on `run_sim_1rkl_rigid_dry.sh`.
   - Script explicitly runs only stages `6.0` -> `6.6` and extracts VTF files for each stage.
   - Updated `batch_relax_rigid_dry.py` default `--run-script` to this new standalone script and locked `--target-stage` to `6.6`.
   - Verified on smoke run (`1ors`) that checkpoints and VTF outputs include `stage_6.0..6.6` and no `stage_7.0` output.
11. Ported bugfixes from `/Users/yinhan/Documents/upside2-md-martini` commit range:
   - Source diff basis: `dd111dbffd9ab4842860ea6628e4a966bda3a049 -> 12b944d1bc743ea1104757cb94d2b49277b8295a`.
   - Applied `src/box.h` updates:
     - Added `NodeBoxUpdater` callback type and `register_node_box_updater(...)` API.
   - Applied `src/box.cpp` updates:
     - Added atomic callback storage and implementation of `register_node_box_updater(...)`.
     - Replaced `update_node_boxes(...)` placeholder with callback dispatch.
   - Applied `src/martini.cpp` updates:
     - `MartiniPotential` now stores/reads box dimensions and applies minimum-image to pair displacement in force evaluation.
     - Added explicit NPT box-scaling propagation callback (`update_martini_node_boxes`) for `martini_potential`, `dist_spring`, `angle_spring`, `dihedral_spring`.
     - Registered callback via `simulation_box::npt::register_node_box_updater(...)` during node registrar initialization.
   - Verification:
     - `cmake --build /tmp/upside2-md-build-rbm -j4` succeeded after port.
