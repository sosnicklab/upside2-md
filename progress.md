# Progress Log

## 2026-05-24 MARTINI Workflow Cleanup
- Actions taken:
  - Replaced stale drift-debug task notes with the current cleanup plan.
  - Added `run_sim_hybrid.sh` as the shared MARTINI launcher.
  - Converted 1RKL/1AFO coarse/full scripts into thin wrappers with system and
    resolution defaults.
  - Removed Python debug-only workflow switches and debug PDB/JSON/TSV writers
    from production preparation paths.
  - Consolidated stage-7 production handoff into a single helper used by both
    direct and extended equilibration flows.
  - Added `src/martini.h`, reduced `src/main.h` to the C entry point, and routed
    C++ MARTINI users through the new header.
  - Removed the C++ potential-component dump flag and shortened hybrid progress
    output to `protein_potential` plus `total_potential`.
  - Fixed `run_sim_hybrid.sh` so legacy `source.sh` is sourced before enabling
    `set -u`, avoiding unbound `PYTHONPATH` failures from clean shells.
  - Consolidated MARTINI ITP topology parsing into `py/martini_itp_reader.py`
    and removed the duplicated parser from `py/martini_prepare_system_lib.py`.
  - Cleaned branch-only MARTINI Python/C++ style: removed banner sections,
    debug/generator-looking wording, verbose conversion printouts, and
    non-ASCII progress/comment text while leaving physical interactions intact.
  - Restored `src/main.h` to master parity and kept `src/main.cpp` closer to
    master output/threading style.
- Files modified:
  - `plan.md`
  - `progress.md`
  - `example/16.MARTINI/run_sim_hybrid.sh`
  - `example/16.MARTINI/run_sim_1rkl.sh`
  - `example/16.MARTINI/run_sim_1afo.sh`
  - `example/16.MARTINI/run_sim_1rkl_full.sh`
  - `example/16.MARTINI/run_sim_1afo_full.sh`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `py/martini_itp_reader.py`
  - `py/martini_build_tables.py`
  - `py/martini_extract_vtf.py`
  - `py/martini_gen_params.py`
  - `src/martini.h`
  - `src/main.h`
  - `src/main.cpp`
  - `src/box.cpp`
  - `src/martini.cpp`
  - `src/martini_cg_lipid.cpp`
  - `src/deriv_engine.cpp`
  - `src/thermostat.cpp`
  - `src/martini_hybrid_runtime.h` removed after merging its declaration into
    `src/martini.h`.
- Verification:
  - `bash -n` passed for all MARTINI run scripts.
  - Reproduced the launcher from a clean shell with `PYTHONPATH`,
    `CPLUS_INCLUDE_PATH`, `LIBRARY_PATH`, and `LD_LIBRARY_PATH` unset; it now
    reaches argument parsing instead of failing in `source.sh`.
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_gen_params.py py/martini_build_tables.py py/martini_cg_lipid_params.py py/martini_extract_vtf.py py/martini_itp_reader.py` passed.
  - `python3 py/martini_prepare_system.py run-hybrid-workflow --help` passed
    and exposes no removed debug/restart compatibility flags.
  - Final scan found no active debug workflow controls or component dump format
    strings in the edited files.
  - Final ASCII/style scan found no non-ASCII text in edited MARTINI Python/C++
    files.
  - `git diff --check` passed for the edited files.
  - `cmake --build obj` passed; remaining warnings are pre-existing compiler
    warnings in unrelated code paths.
