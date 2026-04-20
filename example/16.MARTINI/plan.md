# Out-of-Bilayer 1RKL Workflow

## Project Goal
- Create `example/16.MARTINI/run_sim_1rkl_outlipid.sh` based on `run_sim_1rkl.sh`.
- Start the `1rkl` protein outside the bilayer, laid on top of the upper leaflet instead of being prepared with the embedded placement path.
- Increase the simulation box so the outside-start system has room if the protein unfolds.
- Generate a reviewable prepared PDB and keep the workflow directly submitable with `sbatch`.
- Add a continuation mode so the out-of-bilayer wrapper can resume production from a previous `stage_7.0.up`.
- Ensure VTF extraction writes one file per continued trajectory segment instead of flattening archived `output_previous_*` groups into one output.

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
  - single-segment extractor smoke test:
    - `python py/martini_extract_vtf.py ... --split-segments` on `outputs/martini_test_1rkl_outlipid_envfix_smoke/checkpoints/1rkl.stage_7.0.continue.up`
    - observed one output file: `/tmp/segment_single_check.vtf`
  - synthetic multi-segment extractor smoke test:
    - created `/tmp/1rkl.stage_7.0.segment_test2.up` with both `output_previous_0` and `output`
    - ran `python py/martini_extract_vtf.py ... --split-segments`
    - observed output files: `/tmp/segment_multi_check.segment_0.vtf` and `/tmp/segment_multi_check.segment_1.vtf`
