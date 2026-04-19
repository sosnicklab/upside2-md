# Out-of-Bilayer 1RKL Workflow

## Project Goal
- Create `example/16.MARTINI/run_sim_1rkl_outlipid.sh` based on `run_sim_1rkl.sh`.
- Start the `1rkl` protein outside the bilayer, laid on top of the upper leaflet instead of being prepared with the embedded placement path.
- Increase the simulation box so the outside-start system has room if the protein unfolds.
- Generate a reviewable prepared PDB and keep the workflow directly submitable with `sbatch`.

## Architecture & Key Decisions
- Keep `example/16.MARTINI/run_sim_1rkl.sh` backward-compatible by adding optional preparation controls only.
- Implement the new starting geometry in `py/martini_prepare_system.py`, because that is the actual stage-0 mixed-system builder used by the workflow.
- Make the new placement explicit:
  - preserve existing centered/embed behavior as the default,
  - add an outside-top placement mode for the new workflow,
  - add a deterministic "lay flat" orientation so the protein is placed on the bilayer surface rather than only shifted in `z`.
- Size the out-of-bilayer system with larger padding defaults in the new wrapper script instead of changing the base workflow defaults.
- Use the stage-0 packed MARTINI PDB as the review artifact so the reviewed structure matches the actual simulation input geometry.

## Execution Phases
- [x] Phase 1: Update the prep script with optional outside-top placement and flat orientation controls while preserving current defaults.
- [x] Phase 2: Wire the new controls into `run_sim_1rkl.sh` through optional environment variables only.
- [x] Phase 3: Add `run_sim_1rkl_outlipid.sh` with Slurm headers and larger-box defaults.
- [x] Phase 4: Generate a review PDB for the new starting structure and inspect the placement summary.
- [x] Phase 5: Run syntax / compile verification and record the final review notes.

## Known Errors / Blockers
- The review PDB will be MARTINI resolution because the actual packed simulation input is a MARTINI structure.
- No open blocker from the implementation pass.

## Review
- Implemented:
  - optional mixed-system placement/orientation controls in `py/martini_prepare_system.py`
  - opt-in env wiring in `example/16.MARTINI/run_sim_1rkl.sh`
  - new Slurm wrapper `example/16.MARTINI/run_sim_1rkl_outlipid.sh`
- Review artifact:
  - `example/16.MARTINI/outputs/martini_test_1rkl_outlipid/hybrid_prep/hybrid_packed.MARTINI.pdb`
- Observed prepared geometry:
  - placement = `outside-top`
  - orientation = `lay-flat`
  - requested protein-top clearance = `6.0 Å`
  - final box = `177.472 x 177.472 x 160.621 Å`
  - protein span after flat layout = `43.177 x 16.204 x 11.812 Å`
  - removed lipid residues = `0`
- Verification:
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `bash -n example/16.MARTINI/run_sim_1rkl_outlipid.sh`
  - `python` source-compilation check for `py/martini_prepare_system.py`
  - direct stage-0 preparation run using the new outside-top + lay-flat options, writing summary JSON and packed PDB under `outputs/martini_test_1rkl_outlipid/hybrid_prep/`
