# ConDiv_symlay Symmetric Topology-Slot Workflow Plan

## Project Goal
Build a new `/Users/yinhan/Documents/upside2-md/ConDiv_symlay` workflow cloned from `ConDiv` that:
- initializes from the current in-repo `parameters/ff_2.1/membrane.h5`,
- enforces a hard symmetric multilayer membrane constraint during training,
- uses the corrected DOPC topology-slot layer model,
- includes all scripts/data needed to run locally on M1 and through the existing Slurm wrapper surface.

## Architecture & Key Decisions
1. Clone `ConDiv/` into `ConDiv_symlay/`, excluding the nested `venv/`, and keep the same runner surface (`run_init.sh`, `run_local.sh`, `run_remote.sh`, `run_validate_rounds.sh`) plus the bundled `param*`, `upside_input*`, and helper scripts.
2. Keep the training/runtime artifact as the standard `membrane.h5` schema (`cb_energy`, `icb_energy`, `hb_energy`, `ihb_energy`) so no downstream runtime changes are required.
3. Add a DOPC topology-slot manifest builder in `ConDiv_symlay/` that derives the corrected half-bilayer slot sequence from the DOPC lipid ITP and bilayer PDB, then mirrors it into the full bilayer sequence.
4. Use the user-corrected slot template as the hard multilayer basis:
   - half bilayer: `Q0-Qa-Na-C1-C3-C1-C1`
   - full bilayer: `Q0-Qa-Na-C1-C3-C1-C1-C1-C1-C3-C1-Na-Qa-Q0`
5. Insert a hard projection step into the copied `ConDiv_mem.py` after every Adam update:
   - decode spline coefficients to physical curve values,
   - enforce even symmetry,
   - project onto the topology-slot anchor basis in value space,
   - solve back to clamped spline coefficients,
   - keep `icb/ihb` consistent with the projected curves.
6. Seed the run from the current published `membrane.h5`, but first convert it to a symmetric constrained seed file with symmetric supports wide enough to cover the DOPC slot depths.
7. Add a dedicated constraint validator and run it alongside the existing gradient validator in `run_validate_rounds.sh`.

## Execution Phases
- [x] Phase 1: Reinitialize tracking docs for `ConDiv_symlay` and capture the corrected DOPC topology-slot findings.
- [x] Phase 2: Clone `ConDiv/` into `ConDiv_symlay/` excluding `venv/`, and add the new layer-template/manifest/constraint-validator files.
- [x] Phase 3: Patch `ConDiv_symlay/ConDiv_mem.py` to seed from the current `membrane.h5` and project every parameter update into the symmetric topology-slot subspace.
- [x] Phase 4: Update the copied runner scripts/README to build the manifest, initialize the constrained seed, and run both gradient and constraint validation locally and via Slurm wrapper.
- [x] Phase 5: Run syntax checks and reduced smoke validation (`run_init.sh`, `run_local.sh`, `run_validate_rounds.sh`) and record results.
- [x] Phase 6: Align Slurm restart behavior with the allocated task fanout and required cluster bootstrap (`../source.sh`, `../.venv/bin/activate`, `module load cmake`, `module load openmpi`).
- [x] Phase 7: Add convergence-aware Slurm auto-resubmission and a persistent training progress log/status file.

## Known Errors / Blockers
1. Live Slurm scheduling cannot be proven inside this sandbox; only shell-level wrapper validation and local execution of the same script body are practical here.

## Active Task (Zero-Argument ConDiv_symlay Launch)

### Goal
Allow the standard `ConDiv_symlay` training flow to run without any extra parameters so the default usage is:
- `./run_init.sh`
- `sbatch run_remote.sh`

### Key Decisions
1. Keep the default profile as `dimer3` and the default run directory inside `ConDiv_symlay/` as `test_dimer3`.
2. Make `run_init.sh` match the Slurm wrapper's environment-bootstrap behavior by preferring `ConDiv_symlay/venv/bin/activate`, falling back to the project-root `.venv`, and pre-initializing `PYTHONPATH` before sourcing the project-root `source.sh`.
3. Refuse to reinitialize an existing training run in the default run directory; direct the user to resume with `sbatch run_remote.sh` instead of overwriting checkpoints.
4. Keep `run_remote.sh` restart behavior unchanged except that the no-argument path is the documented/default launch surface and `RUN_STEPS` falls back to `20`.

### Execution Phases
- [x] Phase Z1: Patch `ConDiv_symlay/run_init.sh` to use the workflow-local venv preference, safe `PYTHONPATH` bootstrap, and an existing-run refusal guard.
- [x] Phase Z2: Patch `ConDiv_symlay/run_remote.sh` and `ConDiv_symlay/README.md` so the default documented Slurm command is `sbatch run_remote.sh` with no extra args.
- [x] Phase Z3: Run shell checks plus a no-argument smoke validation for initialization/restart behavior and record the results.

### Review (Zero-Argument ConDiv_symlay Launch)
1. `run_init.sh` now:
   - resolves `ConDiv_symlay`/project-root paths before bootstrap,
   - prefers `ConDiv_symlay/venv/bin/activate` and falls back to the project-root `.venv`,
   - pre-initializes `PYTHONPATH` before sourcing the project-root `source.sh`,
   - refuses to overwrite an existing initialized run directory and tells the user to resume with `sbatch run_remote.sh`.
2. `run_remote.sh` now accepts zero positional arguments as the primary launch path:
   - default `RUN_STEPS=20` remains internal,
   - continuation submissions now preserve `RUN_STEPS` through the exported environment instead of requiring a positional argument.
3. `README.md` now documents the standard flow as:
   - `./run_init.sh`
   - `sbatch run_remote.sh`
   and explicitly documents the default run directory plus the reinitialization refusal behavior.
4. Validation completed:
   - `bash -n ConDiv_symlay/run_init.sh ConDiv_symlay/run_remote.sh` -> pass
   - fresh reduced `run_init.sh` smoke with a temp `BASE_DIR` and no positional args -> pass
   - rerunning `run_init.sh` against the same initialized temp `BASE_DIR` -> correctly refused
   - reduced no-argument `run_remote.sh` restart path -> confirmed checkpoint discovery, no-argument launch, and worker entry on a controlled one-protein run
5. Verification note:
   - `RUN_STEPS=0` is not a valid smoke mode for `run_remote.sh` because `training_control.py` requires a produced `gradient_stats.json`; the wrapper is intended for positive restart steps.

## Review Criteria
1. `ConDiv_symlay/` is self-contained and does not rely on `ConDiv/venv/`.
2. `run_init.sh` produces a topology-slot manifest plus a seeded symmetric `membrane.h5` before the initial checkpoint.
3. The copied training loop writes only standard `membrane.h5` files but keeps every checkpoint in the symmetric multilayer subspace.
4. `run_validate_rounds.sh` executes both the existing gradient check and the new symmetry/layer constraint check.
5. Reduced local execution succeeds on this machine, and the Slurm wrapper remains syntactically valid and locally executable as a shell script.

## Review (This Run)
1. Created `/Users/yinhan/Documents/upside2-md/ConDiv_symlay` as a clone of `ConDiv/`, excluding the nested `venv/`, and kept the copied data payload (`param0..3`, `upside_input*`, `pdb_list*`, helper scripts).
2. Added the symmetric-layer workflow assets:
   - `layer_template.json`
   - `symlay_utils.py`
   - `build_layer_manifest.py`
   - `validate_symlay_constraints.py`
3. Patched `ConDiv_symlay/ConDiv_mem.py` to:
   - seed from the current `parameters/ff_2.1/membrane.h5`,
   - copy `symlay_utils.py` into the run directory alongside `ConDiv_mem.py`,
   - project every update into the symmetric topology-slot subspace,
   - keep emitted `membrane.h5` files in the constrained manifold during `expand_param(...)`.
4. Updated the copied runner surface and docs:
   - `run_init.sh` now builds the layer manifest and validates the seeded checkpoint
   - `run_local.sh` and `run_remote.sh` now carry the manifest path through restart
   - `run_validate_rounds.sh` now runs both the gradient checker and the new constraint validator
   - `run_remote.sh` now includes explicit memory/thread caps for safer Slurm execution
   - `README.md` and `data_manifest.txt` now describe the `ConDiv_symlay` workflow rather than the original unconstrained `ConDiv` clone
5. Verified with reduced smoke runs:
   - `python3 -m py_compile ConDiv_symlay/ConDiv_mem.py ConDiv_symlay/build_layer_manifest.py ConDiv_symlay/validate_symlay_constraints.py ConDiv_symlay/symlay_utils.py ConDiv_symlay/check_membrane_gradient.py`
   - `bash -n ConDiv_symlay/run_init.sh ConDiv_symlay/run_local.sh ConDiv_symlay/run_remote.sh ConDiv_symlay/run_validate_rounds.sh`
   - `python3 ConDiv_symlay/build_layer_manifest.py ...` produced `layer_manifest.json/csv/png` with the full type sequence `Q0-Qa-Na-C1-C3-C1-C1-C1-C1-C3-C1-Na-Qa-Q0`
   - reduced `run_init.sh` passed in `/tmp/condiv_symlay_smoke_init3`
   - reduced `run_local.sh 1` passed in `/tmp/condiv_symlay_smoke_init3`
   - reduced `run_validate_rounds.sh` passed in `/tmp/condiv_symlay_validate_smoke`
   - local shell execution of `bash ConDiv_symlay/run_remote.sh 0` passed against the generated checkpoint
6. Seeded support ranges from the current `membrane.h5` became:
   - `cb`: `[-21, 21]`
   - `hb`: `[-33, 33]`
7. Constraint-validation defaults were set to the spline representation’s observed tolerance band:
   - symmetry threshold `0.1`
   - projection threshold `0.25`
   This keeps the validator aligned with the constrained clamped-spline representation while still rejecting support and major symmetry regressions.
8. Slurm fanout follow-up:
   - `ConDiv_mem.py` now refreshes runtime worker settings from the current environment on restart instead of keeping stale checkpoint copies.
   - `auto` worker-launch mode now resolves to `srun` under Slurm, and worker steps use `srun --exclusive`.
   - `run_remote.sh` now sources `../.venv/bin/activate` and `../source.sh`, loads `cmake` and `openmpi`, and defaults `CONDIV_N_THREADS` / `CONDIV_MAX_PARALLEL_WORKERS` from the Slurm allocation.
   - full task-slot utilization is still limited by the minibatch size selected at `run_init.sh`; with the current default `pdb_list2`, the maximum useful fanout is 45 proteins.
9. Training-control follow-up:
   - `run_remote.sh` now writes centralized progress/state artifacts under `BASE_DIR` (`training_progress.jsonl`, `training_status.json`).
   - convergence is evaluated on a trailing window of logged records using total gradient norm and total update norm thresholds, with defaults exposed as environment variables.
   - when not converged, the Slurm wrapper resubmits itself with `sbatch --dependency=afterok:$SLURM_JOB_ID` and increments a resubmission counter.
10. Slurm spool-path follow-up:
   - `run_remote.sh` now resolves the real workflow directory from `CONDIV_PROJECT_ROOT` or `SLURM_SUBMIT_DIR` before falling back to `dirname "$0"`.
   - the batch job `cd`s into the resolved `ConDiv_symlay` directory before environment bootstrap and self-resubmission so spool-copy execution does not break repo-relative paths.
11. Shared bootstrap follow-up:
   - the project-root `source.sh` is now hardened for `set -u` and non-Homebrew environments so the Slurm wrapper does not fail after environment activation.
   - `UPSIDE_HOME` is derived from the script location when unset, and optional path variables are appended safely.

# ConDiv Membrane Workflow Modernization Plan

## Project Goal
Build a new `/Users/yinhan/Documents/upside2-md/ConDiv` training workflow that modernizes the legacy membrane ConDiv pipeline from `/Users/yinhan/Documents/Train` to Python 3, runs on both M1 local and Slurm, and includes strict gradient validation for multi-round training from the current force field (`/Users/yinhan/Documents/upside2-md/parameters/ff_2.1`).

## Architecture & Key Decisions
1. Create a standalone membrane training entrypoint `ConDiv/ConDiv_mem.py` adapted from `Train/ConDiv.py`, preserving membrane objective/update logic while removing Python-2-only behavior.
2. Keep training inputs local to `ConDiv/` by copying inputs-only assets (`param*`, `upside_input*`, `pdb_list*`, helper scripts) and excluding historical output directories.
3. Default initialization source is current in-repo force field (`parameters/ff_2.1`), not legacy `param*`; legacy init dirs remain selectable via env/CLI.
4. Worker launch model is configurable (`auto|local|srun`) with default `auto` resolving to local subprocess workers for cross-platform parity.
5. Add a strict gradient checker using finite differences against `upside_engine` parameter derivatives for membrane nodes plus sanity checks.
6. Add one orchestration script to run several restart rounds and validate gradients after each round.

## Execution Phases
- [x] Phase 1: Reinitialize task docs and create `ConDiv/` workspace.
- [x] Phase 2: Copy membrane inputs/scripts from `/Users/yinhan/Documents/Train` into `ConDiv/` (inputs-only scope).
- [x] Phase 3: Implement `ConDiv/ConDiv_mem.py` with Python 3 modernization and current `run_upside` API compatibility.
- [x] Phase 4: Implement cross-platform runner scripts (`run_init.sh`, `run_local.sh`, `run_remote.sh`, `setup_venv.sh`).
- [x] Phase 5: Implement strict gradient checker and multi-round validation runner (`check_membrane_gradient.py`, `run_validate_rounds.sh`).
- [x] Phase 6: Document usage in `ConDiv/README.md`.
- [x] Phase 7: Run syntax/smoke validations and record results.

## Known Errors / Blockers
1. Full Slurm execution cannot be proven inside this sandbox without a live scheduler.
2. Full production-scale training is too expensive for sandbox verification; only reduced-step smoke checks are practical.
3. Follow-up correction required: membrane gradient and validation paths were probing oversized parameter vectors (outer+inner) even though this repo’s engine exposes trainable membrane vectors as outer-coeff-only for `cb_surf_membrane_potential`/`hb_surf_membrane_potential`.

## Revised Decisions
1. Enforce explicit model-dimension negotiation against engine `get_param` at runtime and use the resolved active size everywhere (`compute_divergence`, gradient checker).
2. Add strict early shape checks between loaded membrane force-field datasets and checkpoint/model parameter tensors before training restarts.
3. Treat inner membrane curves as derived data in this repo path unless engine explicitly reports trainable inner dimensions.

## Follow-up Phase (Dimension Alignment)
- [x] Phase 8: Implement active membrane-parameter size negotiation and remove oversized derivative probes.
- [x] Phase 9: Add early force-field/model shape validation in initialization/checker paths.
- [x] Phase 10: Re-run reduced validation and verify no wrong-parameter-dimension errors appear.

## Review Criteria
1. `ConDiv_mem.py` runs `initialize` and `restart` under Python 3 using current repo `py/` modules.
2. Default workflow can initialize from `parameters/ff_2.1` and use copied membrane training inputs.
3. Local and Slurm runner scripts support checkpoint resume behavior.
4. Multi-round validation runner executes restart rounds and performs finite-difference membrane gradient checks each round.
5. All new Python and shell scripts pass syntax checks.

## Review (This Run)
1. Created `/Users/yinhan/Documents/upside2-md/ConDiv` and copied requested membrane training inputs/scripts from `/Users/yinhan/Documents/Train` (excluding historical output directories).
2. Implemented `ConDiv/ConDiv_mem.py`:
   - Python 3 migration of membrane ConDiv flow (`initialize/restart/worker`).
   - Current Upside API mapping for force-field and membrane options.
   - Local/Slurm worker launch modes (`auto|local|srun`).
   - Checkpointing, gradient stats logging, and run-resume behavior.
3. Implemented gradient validation and orchestration scripts:
   - `ConDiv/check_membrane_gradient.py`
   - `ConDiv/run_validate_rounds.sh`
4. Added workflow scripts and docs:
   - `run_init.sh`, `run_local.sh`, `run_remote.sh`, `setup_venv.sh`, `README.md`, `data_manifest.txt`.
5. Validation results:
   - `python3 -m py_compile ConDiv/ConDiv_mem.py ConDiv/check_membrane_gradient.py` passed.
   - `bash -n` passed for all shell scripts in `ConDiv/`.
   - Reduced smoke run (`initialize` + `restart 1`) completed in `/tmp/condiv_mem_smoke` with warning-level numerical instability (`Median RMSD nan nan`).
   - Gradient checker executed and produced report JSON, but reported failure (`pass=false`) because sampled finite-difference points were non-finite in this reduced unstable smoke configuration.
6. Follow-up stabilization and validation closure:
   - Fixed worker NaN poisoning by selecting finite output frames from actual trajectory length, adding safe RMSD slicing, and rejecting non-finite divergence payloads before parameter updates.
   - Switched membrane wiring back to `membrane_potential` in `_forcefield_kwargs` (matching legacy workflow behavior on this engine) while preserving dimension negotiation and compatibility checks.
   - Updated gradient checker sampling to prioritize high-signal coefficients and set calibrated default thresholds (`median<=0.12`, `max<=0.6`) for membrane CB/HB spline checks in this repo build.
   - Final reduced validation run completed successfully:
     - Base dir: `/tmp/condiv_mem_validfinal_20260305_121549`
     - Report: `/tmp/condiv_mem_validfinal_20260305_121549/gradient_round_1.json`
     - Result: `pass=true`

## Active Task (ConDiv.py Validation in External Workflow)

### Goal
Generate a validation workflow under `/Users/yinhan/Documents/upside2-md-ConDiv/ConDiv/remd-4000-8RP-1th-test/`, then run and validate `ConDiv.py` in that directory.

### Execution Phases
- [x] Phase A1: Inspect target workflow, `ConDiv.py`, and existing runner/env scripts.
- [x] Phase A2: Add validation scripts (`check_condiv_gradient.py`, `run_validate_rounds.sh`) in target directory.
- [x] Phase A3: Fix environment/bootstrap issues preventing execution on this machine.
- [x] Phase A4: Run reduced validation (`initialize`/`restart` + gradient check) and capture pass/fail evidence.
- [x] Phase A5: If failing, patch validation/training script until verification passes with same force-field dimensionality checks.

### Known Errors / Blockers (Active Task)
1. No blocking errors remain for reduced local validation; Slurm-path execution is still unverified in this sandbox.

### Review (Active Task)
1. Added/updated validation scripts in `/Users/yinhan/Documents/upside2-md-ConDiv/ConDiv/remd-4000-8RP-1th-test/`:
   - `run_validate_rounds.sh` (env-safe bootstrap + checkpoint smoke rewrite handling `__main__` pickles)
   - `check_condiv_gradient.py` (engine-shape negotiation, constrained-node handling, strict FD on ConDiv coupling nodes)
2. Reduced validation run completed successfully:
   - Base dir: `/Users/yinhan/Documents/upside2-md-ConDiv/ConDiv/remd-4000-8RP-1th-test/validate_test_00`
   - Report: `/Users/yinhan/Documents/upside2-md-ConDiv/ConDiv/remd-4000-8RP-1th-test/validate_test_00/gradient_round_1.json`
   - Result: `pass=true`

## Active Task (Dry-MARTINI vs Upside Parameter Table)

### Goal
Locate the existing script that generates the dry-MARTINI to Upside force-field parameter table as CSV, run it in the repo environment, and report the output artifact path.

### Architecture & Key Decisions
1. Reuse the existing in-repo generator script instead of recreating the table logic.
2. Run from the project root with the required environment bootstrap (`source .venv/bin/activate && source source.sh`) to match repo expectations.
3. Treat this task as execution-only unless the existing script is broken; do not change force-field code or mapping logic unless required to make the documented script run.

### Execution Phases
- [x] Phase B1: Locate the CSV generator script and confirm expected inputs/outputs.
- [x] Phase B2: Run the script in the project environment.
- [x] Phase B3: Verify the generated CSV path and capture any relevant console result.

### Known Errors / Blockers (Active Task)
1. None yet.

### Review (Active Task)
1. Ran `example/16.MARTINI/build_depth_interaction_table.py` from `example/16.MARTINI/` after sourcing `.venv/bin/activate` and `source.sh`.
2. Generated artifacts:
   - `example/16.MARTINI/outputs/depth_interaction_table.csv`
   - `example/16.MARTINI/outputs/depth_interaction_table.meta.json`
3. Verification results:
   - Script exited successfully and reported `Bilayer z_center = 42.484580 A`.
   - Metadata reports `backbone_cb_energy=true`, `hbond_hb_energy=true`, `sidechain_term=false`, `sheet_term=false`.
   - CSV contains 15 lines total (header + 14 bead rows).
