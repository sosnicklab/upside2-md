# 2026-05-26 MARTINI H5 Rebuild Scripts

## Project Goal
- Add two scripts under `example/16.MARTINI` that regenerate every
  dry-MARTINI `.h5` file required by the example workflows.
- The local script must use the M1/local environment bootstrap and run table
  generation in parallel.
- The Slurm script must use the cluster module/venv bootstrap and run table
  generation in parallel using allocated CPUs.

## Architecture & Key Decisions
- Use the existing single source of truth, `py/martini_gen_params.py --force`,
  so both scripts generate the same files under `parameters/dryMARTINI`:
  `particle.h5`, `sidechain.h5`, `dopc.h5`, and `interlipid.h5`.
- Local M1 setup follows the project rule: activate `.venv`, then source
  `source.sh`, then call the Python generator with an explicit worker count.
- Slurm setup avoids the local `source.sh` bootstrap and instead loads modules,
  activates the repo `.venv`, exports `UPSIDE_HOME`, `PATH`, and `PYTHONPATH`,
  and sets `UPSIDE_SKIP_SOURCE_SH=1`.
- Worker count remains overrideable with `UPSIDE_MARTINI_TABLE_WORKERS`; when
  unset, the local script uses local logical CPU count and the Slurm script
  uses `SLURM_CPUS_PER_TASK`.

## Execution Phases
- [x] Phase 1: Add local M1 rebuild script.
- [x] Phase 2: Add Slurm rebuild script.
- [x] Phase 3: Verify shell syntax and record usage.

## Known Errors / Blockers
- Full table regeneration is intentionally not run during script verification
  because it is long-running and rewrites production parameter files.

## Review
- Added `example/16.MARTINI/build_martini_h5_m1.sh` for local M1 rebuilds.
  It activates `.venv`, sources `source.sh`, sets an explicit local worker
  count from logical CPUs when not overridden, and calls
  `py/martini_gen_params.py --force`.
- Added `example/16.MARTINI/build_martini_h5_slurm.sh` for Slurm rebuilds.
  It uses the module/venv bootstrap, sets `UPSIDE_SKIP_SOURCE_SH=1`, uses
  `SLURM_CPUS_PER_TASK` when `UPSIDE_MARTINI_TABLE_WORKERS` is unset, and
  calls the same generator.
- Both scripts forward extra command-line arguments to `martini_gen_params.py`,
  so optional table controls such as `--bead-frame-count` remain available.
