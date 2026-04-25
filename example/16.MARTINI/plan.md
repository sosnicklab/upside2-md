# Out-of-Bilayer 1RKL Workflow

## Project Goal
- Create `example/16.MARTINI/run_sim_1rkl_outlipid.sh` based on `run_sim_1rkl.sh`.
- Start the `1rkl` protein outside the bilayer, laid on top of the upper leaflet instead of being prepared with the embedded placement path.
- Increase the simulation box so the outside-start system has room if the protein unfolds.
- Generate a reviewable prepared PDB and keep the workflow directly submitable with `sbatch`.
- Add a continuation mode so the out-of-bilayer wrapper can resume production from a previous `stage_7.0.up`.
- Ensure VTF extraction writes one file per continued trajectory segment instead of flattening archived `output_previous_*` groups into one output.
- Ensure the workflow seeds are generated per run by default instead of using fixed literal numbers.

## Architecture & Key Decisions
- Keep `example/16.MARTINI/run_sim_1rkl.sh` backward-compatible by adding optional preparation controls only.
- Implement the new starting geometry in `py/martini_prepare_system.py`, because that is the actual stage-0 mixed-system builder used by the workflow.
- Make the new placement explicit:
  - preserve existing centered/embed behavior as the default,
  - add an outside-top placement mode for the new workflow,
  - add a deterministic "lay flat" orientation so the protein is placed on the bilayer surface rather than only shifted in `z`.
- Size the out-of-bilayer system with larger padding defaults in the new wrapper script instead of changing the base workflow defaults.
- Use the stage-0 packed MARTINI PDB as the review artifact so the reviewed structure matches the actual simulation input geometry.
- Implement continuation in `run_sim_1rkl.sh` as an optional stage-7-only resume path so the wrapper can reuse existing production restart logic without duplicating the execution command.
- Default continuation output to a new `stage_7.0.continue.up` file to avoid mutating the previous run artifact unless the user explicitly asks for in-place continuation.
- Handle multi-segment continuation at VTF extraction time:
  - keep the existing `${stage}.vtf` name for single-segment stage files,
  - emit `${stage}.segment_<n>.vtf` files when a stage file contains `output_previous_*` groups plus the current `output`.
- Keep seed behavior explicit:
  - `PREP_SEED` and `SEED` should be generated at runtime when unset,
  - explicit user-provided env vars must still override the generated values for reproducibility.

## Execution Phases
- [x] Phase 1: Update the prep script with optional outside-top placement and flat orientation controls while preserving current defaults.
- [x] Phase 2: Wire the new controls into `run_sim_1rkl.sh` through optional environment variables only.
- [x] Phase 3: Add `run_sim_1rkl_outlipid.sh` with Slurm headers and larger-box defaults.
- [x] Phase 4: Generate a review PDB for the new starting structure and inspect the placement summary.
- [x] Phase 5: Run syntax / compile verification and record the final review notes.
- [x] Phase 6: Add stage-7 continuation support to `run_sim_1rkl.sh` and expose it through the out-of-bilayer wrapper.
- [x] Phase 7: Verify the continuation path and record the updated review notes.
- [x] Phase 8: Fix Slurm spool-copy path resolution in the out-of-bilayer wrapper and verify it with a spool-style smoke test.
- [x] Phase 9: Align the Slurm wrapper environment with a self-contained cluster setup and document the rule in `AGENTS.md`.
- [x] Phase 10: Make VTF extraction segment-aware for continued stage files and verify the multi-segment naming path.
- [x] Phase 11: Replace fixed workflow seed defaults with runtime-generated seeds and verify both the base workflow and the out-of-bilayer wrapper.
- [x] Phase 12: Remove the dead `martini_integration_stage(...)` helper from `src/martini.cpp` and verify the active MARTINI mass-aware path still builds.

## Known Errors / Blockers
- The review PDB will be MARTINI resolution because the actual packed simulation input is a MARTINI structure.
- No open blocker from the implementation pass.

## Review
- Implemented:
  - optional mixed-system placement/orientation controls in `py/martini_prepare_system.py`
  - opt-in env wiring in `example/16.MARTINI/run_sim_1rkl.sh`
  - new Slurm wrapper `example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - optional `CONTINUE_STAGE_70_FROM` / `CONTINUE_STAGE_70_OUTPUT` production-resume path in `example/16.MARTINI/run_sim_1rkl.sh`
  - wrapper aliases `PREVIOUS_RUN_DIR` and `PREVIOUS_STAGE7_FILE` in `example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - segment-aware VTF extraction in `py/martini_extract_vtf.py` and `example/16.MARTINI/run_sim_1rkl.sh`
  - runtime-generated default `PREP_SEED` and `SEED` in `example/16.MARTINI/run_sim_1rkl.sh`, with wrapper pass-through documentation in `example/16.MARTINI/run_sim_1rkl_outlipid.sh`
- Review artifact:
  - `example/16.MARTINI/outputs/martini_test_1rkl_outlipid/hybrid_prep/hybrid_packed.MARTINI.pdb`
- Continuation artifact:
  - `example/16.MARTINI/outputs/martini_test_1rkl_outlipid_continue_smoke/checkpoints/1rkl.stage_7.0.continue.up`
- Observed prepared geometry:
  - placement = `outside-top`
  - orientation = `lay-flat`
  - requested protein-top clearance = `6.0 Å`
  - final box = `177.472 x 177.472 x 160.621 Å`
  - protein span after flat layout = `43.177 x 16.204 x 11.812 Å`
  - removed lipid residues = `0`
- Observed continuation behavior:
  - wrapper accepts `PREVIOUS_RUN_DIR` and maps it to the previous `checkpoints/1rkl.stage_7.0.up`
  - default continuation output is a new `stage_7.0.continue.up` file in the current `RUN_DIR`
  - one-step production continuation completed successfully and wrote `1rkl.stage_7.0_continue.vtf`
  - continued stage files with multiple archived trajectory groups now extract to per-segment files such as `stage_7.0_continue.segment_0.vtf` and `stage_7.0_continue.segment_1.vtf`
- Observed seed behavior:
  - fixed defaults `PREP_SEED=2026` and `SEED=7090685331` were removed
  - two clean-environment base-workflow runs generated different seed pairs:
    - base A: `PREP_SEED=230538232`, `SEED=2946250751`
    - base B: `PREP_SEED=702128154`, `SEED=3561662407`
  - two clean-environment wrapper runs generated different seed pairs:
    - wrapper A: `PREP_SEED=1002878124`, `SEED=479116090`
    - wrapper B: `PREP_SEED=3764234752`, `SEED=606192356`
  - explicit `PREP_SEED` / `SEED` env vars still override the generated defaults
- Slurm constraint:
  - `run_sim_1rkl_outlipid.sh` now uses `#SBATCH --time=36:00:00` to match the user-provided hard cluster limit
  - `run_sim_1rkl_outlipid.sh` now resolves the base workflow from `SLURM_SUBMIT_DIR`, `PWD`, or `BASE_WORKFLOW_SCRIPT` instead of relying only on `BASH_SOURCE[0]`
  - `run_sim_1rkl_outlipid.sh` now loads the cluster modules and activates `PROJECT_ROOT/.venv` before running the workflow
  - `run_sim_1rkl.sh` now supports `UPSIDE_SKIP_SOURCE_SH=1` so Slurm wrappers can bypass the local Apple-Silicon bootstrap
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - no-write Python syntax check for `py/martini_extract_vtf.py` via `compile(...)`
  - `python` source-compilation check for `py/martini_prepare_system.py`
  - direct stage-0 preparation run using the new outside-top + lay-flat options, writing summary JSON and packed PDB under `outputs/martini_test_1rkl_outlipid/hybrid_prep/`
  - continuation smoke test:
    - `PREVIOUS_RUN_DIR=/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_test_1rkl_hybrid`
    - `RUN_DIR=/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_test_1rkl_outlipid_continue_smoke`
    - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - spool-style wrapper smoke test from `/tmp` with `SLURM_SUBMIT_DIR=/Users/yinhan/Documents/upside2-md/example/16.MARTINI`
  - wrapper-managed environment smoke test:
    - `PREVIOUS_RUN_DIR=/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_test_1rkl_hybrid`
    - `RUN_DIR=/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_test_1rkl_outlipid_envfix_smoke`
    - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - clean-environment base-workflow seed smoke tests:
    - `env -i ... CONTINUE_STAGE_70_FROM=... RUN_DIR=...seedcheck_base_a PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl.sh`
    - `env -i ... CONTINUE_STAGE_70_FROM=... RUN_DIR=...seedcheck_base_b PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl.sh`
    - observed different generated seed pairs and successful continuation in both runs
  - clean-environment wrapper seed smoke tests:
    - `env -i ... PREVIOUS_RUN_DIR=... RUN_DIR=...seedcheck_wrapper_a PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl_outlipid.sh`
    - `env -i ... PREVIOUS_RUN_DIR=... RUN_DIR=...seedcheck_wrapper_b PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash example/16.MARTINI/run_sim_1rkl_outlipid.sh`
    - observed different generated seed pairs and successful continuation in both runs
  - single-segment extractor smoke test:
    - `python py/martini_extract_vtf.py ... --split-segments` on `outputs/martini_test_1rkl_outlipid_envfix_smoke/checkpoints/1rkl.stage_7.0.continue.up`
    - observed one output file: `/tmp/segment_single_check.vtf`
  - synthetic multi-segment extractor smoke test:
    - created `/tmp/1rkl.stage_7.0.segment_test2.up` with both `output_previous_0` and `output`
    - ran `python py/martini_extract_vtf.py ... --split-segments`
    - observed output files: `/tmp/segment_multi_check.segment_0.vtf` and `/tmp/segment_multi_check.segment_1.vtf`
- Integrator audit (2026-04-20):
  - `run_sim_1rkl.sh` launches all minimization and MD stages with `--integrator v`, so the workflow stays on the standard Verlet CLI path.
  - The runtime still becomes MARTINI mass-aware because `src/main.cpp` always calls `martini_masses::load_masses_for_engine(...)`, and the generated workflow `.up` files contain `/input/mass`.
  - The actual mass-aware position updates happen inside the inline `martini_masses::has_masses(...)` branches in `src/deriv_engine.cpp`.
  - The custom helper `martini_masses::martini_integration_stage(...)` in `src/martini.cpp` is not wired into the execution path and has no call sites.
- Dead-helper cleanup (2026-04-20):
  - removed `martini_masses::martini_integration_stage(...)` from `src/martini.cpp`
  - confirmed no remaining source references under `src/`
  - verified a clean out-of-tree configure/build from `src/CMakeLists.txt` in `/tmp/upside2-md-stagecheck`

## 2026-04-24 (Workflow Port: 1AFO + Temp-Repo Wrapper Sync)

### Project Goal
- Add `run_sim_1afo.sh` and `run_sim_1afo_outlipid.sh` to this repo and align the existing `1rkl` outlipid wrapper with the temp-repo workflow behavior, while keeping the current repo on its hybrid prep/runtime path.

### Architecture & Key Decisions
- Do not copy the temp repo's full `run_sim_1rkl.sh`.
  - This repo's base workflow is the hybrid workflow and already owns the correct prep/stage-injection path.
- Implement `1afo` as thin wrappers around this repo's `run_sim_1rkl.sh`.
  - Override only the PDB id, input PDB path, run-directory defaults, outlipid placement defaults, and continuation autodetection rules.
- Patch the shared hybrid prep/runtime path for multi-chain `1AFO`.
  - The wrappers alone are insufficient because the current hybrid metadata still collapses chain-local residue numbers and does not propagate chain-break metadata into stage-7 topology reconstruction.
- Keep continuation compatibility with this repo's current base workflow.
  - Superseded on 2026-04-25: wrapper autodetection now uses numeric `stage_7.<n>.up` names only.

### Execution Phases
- [x] Phase 1: Patch shared hybrid metadata/runtime code for unique multi-chain backbone residue indexing and chain-break propagation.
- [x] Phase 2: Patch VTF export so appended AA-backbone atoms preserve chain ids and skip cross-chain peptide links.
- [x] Phase 3: Add `run_sim_1afo.sh` and `run_sim_1afo_outlipid.sh` with Slurm-safe bootstrap and wrapper-level continuation autodetection.
- [x] Phase 4: Update `run_sim_1rkl_outlipid.sh` to the same wrapper autodetection/bootstrap pattern adapted to this repo.
- [x] Phase 5: Run shell syntax checks and a reduced local `1afo` smoke path, then document review notes.

### Known Errors / Blockers
- No open blocker from the implementation pass.

### Review
- Implemented:
  - multi-chain-safe hybrid metadata in `py/martini_prepare_system_lib.py`
  - chain-aware AA-backbone VTF export in `py/martini_extract_vtf.py`
  - optional `chain_break` propagation during hybrid stage-file injection in `example/16.MARTINI/run_sim_1rkl.sh`
  - multi-chain runtime ITP materialization from martinize `.top` in `example/16.MARTINI/run_sim_1rkl.sh`
  - new hybrid wrappers `example/16.MARTINI/run_sim_1afo.sh` and `example/16.MARTINI/run_sim_1afo_outlipid.sh`
  - updated outlipid wrapper autodetection/bootstrap in `example/16.MARTINI/run_sim_1rkl_outlipid.sh`
- Verified:
  - `bash -n` on all four workflow entrypoints
  - `.venv/bin/python3 -m py_compile py/martini_prepare_system_lib.py py/martini_extract_vtf.py`
  - reduced `1afo` end-to-end hybrid smoke:
    - `RUN_DIR=outputs/1afo_port_smoke2`
    - all minimization / equilibration / production lengths reduced to `1`
    - completed through `outputs/1afo_port_smoke2/checkpoints/1afo.stage_7.0.up`
  - stage-7 artifact audit:
    - `n_bb = 72`
    - `unique_bb_residue_index = 72`
    - `bb_chain_ids = ['A', 'B']`
    - `chain_first_residue = [36]`
    - `chain_counts = [1, 1]`
    - no `Distance3D` bond between sequence indices `35` and `36`
    - `example/16.MARTINI/outputs/1afo_port_smoke2/1afo.stage_7.0.vtf` contains both `segid sA chain A` and `segid sB chain B`

## 2026-04-25 (Numeric Continuation Trajectory Naming)

### Project Goal
- Replace ugly `_continue` checkpoint and VTF names with numeric continuation stage names.
- Ensure repeated continuations naturally produce the next numbered trajectory instead of names like `_continue_continue`.

### Architecture & Key Decisions
- Keep the existing stage `7.0` production checkpoint as the first production trajectory segment.
- Name resumed production segments as `stage_7.1`, `stage_7.2`, etc.
- Determine the next continuation index from existing checkpoints in the target `RUN_DIR` and the source checkpoint name, so reruns can continue from either a previous run directory or an explicit previous stage file.
- Keep explicit `CONTINUE_STAGE_70_OUTPUT` as an override for compatibility, but require it to use numeric `stage_7.<n>.up` naming.
- This supersedes the earlier `.continue` default naming decision.
- After user correction, all current bash workflow autodetection should use numeric `stage_7.<n>.up` names only; legacy `.continue` matching is no longer kept in bash wrappers.

### Execution Phases
- [x] Phase 1: Inspect continuation output selection and VTF naming in the workflow scripts.
- [x] Phase 2: Patch continuation naming to choose numeric `stage_7.<n>` checkpoint and VTF outputs by default.
- [x] Phase 3: Update wrapper autodetection if it still prefers legacy `.continue` names.
- [x] Phase 4: Run syntax checks and a reduced continuation smoke test.
- [x] Phase 5: Remove legacy `.continue` file matching from all bash workflows and verify `1afo` continuation output.

### Known Errors / Blockers
- No blocker identified yet.

### Review
- Implemented in `run_sim_1rkl.sh`:
  - removed the default `stage_7.0.continue.up` checkpoint name,
  - removed the default `7.0_continue` VTF label,
  - added numeric continuation output resolution for `stage_7.1`, `stage_7.2`, etc.,
  - kept explicit `CONTINUE_STAGE_70_OUTPUT` and `CONTINUE_STAGE_70_LABEL` overrides only for numeric stage names.
- Wrapper review:
  - `run_sim_1rkl_outlipid.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_outlipid.sh` now auto-detect only numeric `stage_7.<n>.up` continuation checkpoints.
  - `run_sim_1afo.sh` delegates to the same numeric continuation resolver and was verified directly.
- Verification:
  - `bash -n` passed for `run_sim_1rkl.sh`, `run_sim_1rkl_outlipid.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_outlipid.sh`.
  - reduced continuation from `1rkl.stage_7.0.up` wrote `outputs/numeric_continue_smoke/checkpoints/1rkl.stage_7.1.up` and `outputs/numeric_continue_smoke/1rkl.stage_7.1.vtf`.
  - reduced continuation from `1rkl.stage_7.1.up` wrote `outputs/numeric_continue_smoke/checkpoints/1rkl.stage_7.2.up` and `outputs/numeric_continue_smoke/1rkl.stage_7.2.vtf`.
  - reduced continuation through `run_sim_1afo.sh` from `1afo.stage_7.0.up` wrote `outputs/numeric_continue_1afo_smoke/checkpoints/1afo.stage_7.1.up` and `outputs/numeric_continue_1afo_smoke/1afo.stage_7.1.vtf`.
  - code search found no remaining `.continue` filename matching or hard-coded default `stage_7.0.continue` / `7.0_continue` names in the workflow entrypoints.

## 2026-04-25 (1RKL Outlipid Continuation Autodetect Fix)

### Project Goal
- Fix `run_sim_1rkl_outlipid.sh` continuation from a previous trajectory.
- Preserve numeric-only continuation checkpoint naming.

### Architecture & Key Decisions
- Resolve previous outlipid checkpoints from the real workflow directory derived from `BASE_WORKFLOW_SCRIPT`, not from `SCRIPT_DIR`.
- Keep explicit `PREVIOUS_RUN_DIR`, `PREVIOUS_STAGE7_FILE`, and `CONTINUE_STAGE_70_FROM` behavior unchanged.
- Keep autodetection numeric-only: `1rkl.stage_7.<n>.up`.

### Execution Phases
- [x] Phase 1: Reproduce/confirm the wrapper continuation failure path.
- [x] Phase 2: Patch `run_sim_1rkl_outlipid.sh` to use the real workflow directory for autodetection.
- [x] Phase 3: Verify local continuation and Slurm-spool-style autodetection.

### Known Errors / Blockers
- No blocker identified yet.

### Review
- Implemented in `run_sim_1rkl_outlipid.sh`:
  - added `WORKFLOW_DIR` from the resolved `BASE_WORKFLOW_SCRIPT`,
  - changed auto-continuation search from `${SCRIPT_DIR}/outputs` to `${WORKFLOW_DIR}/outputs`.
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - explicit `PREVIOUS_RUN_DIR` continuation through the wrapper wrote:
    - `outputs/martini_test_1rkl_outlipid_continue_fix_smoke/checkpoints/1rkl.stage_7.1.up`
    - `outputs/martini_test_1rkl_outlipid_continue_fix_smoke/1rkl.stage_7.1.vtf`
  - Slurm-spool-style wrapper run from `/tmp` with `SLURM_SUBMIT_DIR` set auto-detected the prior outlipid numeric checkpoint and wrote:
    - `outputs/martini_test_1rkl_outlipid_continue_fix_smoke/checkpoints/1rkl.stage_7.2.up`
    - `outputs/martini_test_1rkl_outlipid_continue_fix_smoke/1rkl.stage_7.2.vtf`

## 2026-04-25 (1AFO Outlipid Stage Prep OOM)

### Project Goal
- Fix `run_sim_1afo_outlipid.sh` freezing locally and being OOM-killed on Slurm during stage-file preparation.
- Keep the out-of-bilayer start geometry usable without requiring excessive Slurm memory.

### Architecture & Key Decisions
- Diagnose the stage-prep memory source before changing resource requests.
- Keep the shared base workflow behavior intact unless the root cause is shared.
- Prefer reducing pathological wrapper-generated system size or avoiding unnecessary dense allocations over raising memory.

### Execution Phases
- [x] Phase 1: Inspect `1afo` outlipid wrapper defaults, generated prep artifacts, and the killed stage-prep code path.
- [x] Phase 2: Identify the memory blow-up and patch the smallest responsible code path.
- [x] Phase 3: Verify `1afo` outlipid reduced workflow reaches stage `7.0` or the formerly killed stage-prep call without freezing/OOM.

### Known Errors / Blockers
- No open blocker from the implementation pass.

### Review
- Root cause:
  - the shared MARTINI stage writer built a full Python list of all non-bonded pairs and a full per-pair `float32[4]` coefficient table before writing HDF5,
  - the `1afo_outlipid` box has `17966` MARTINI atoms, which produces `161363448` non-bonded pairs after exclusions,
  - that dense Python-side allocation caused the local freeze and Slurm OOM during the stage-prep conversion call.
- Implemented in `py/martini_prepare_system_lib.py`:
  - stream non-bonded `pairs` directly into an HDF5 extendable array,
  - store unique coefficient rows once in `coefficients`,
  - store one integer `coefficient_indices` entry per pair,
  - mark the MARTINI group with `optimized_format = 1`, matching the existing optimized runtime loader in `src/martini.cpp`.
- Verification:
  - no-write source compile of `py/martini_prepare_system_lib.py` passed,
  - reduced `1afo` base workflow completed through `outputs/optimized_format_1afo_smoke/checkpoints/1afo.stage_7.0.up`,
  - reduced `1afo_outlipid` wrapper completed through `outputs/oomfix_1afo_outlipid_wrapper_smoke/checkpoints/1afo.stage_7.0.up`,
  - the formerly killed prep path wrote `161363448` non-bonded pairs with `32` unique coefficient rows,
  - runtime loaded the optimized table with `MARTINI: Using optimized interaction table format`.
