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
