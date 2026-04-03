# Task Plan

## Project Goal
- Remove legacy spline/debug text file generation from the MARTINI runtime:
  - `all_splines.txt`
  - `angle_splines.txt`
  - `bond_splines.txt`
  - `dihedral_splines.txt`
  - `force_debug.txt`
- Keep the cleaned `example/16.MARTINI` Python layout intact while ensuring the runtime and preparation path no longer emit or configure those debug files.

## Architecture & Key Decisions
- Remove debug file emission at the source in `src/martini.cpp` rather than relying on workflow cleanup.
- Remove preparation-side attributes that only existed to drive spline/force debug file output.
- Keep runtime physics unchanged apart from deleting unused file-debug behavior.

## Execution Phases
- [x] Phase 1: Audit runtime and preparation code paths that still generate or configure the legacy spline/debug text files.
- [x] Phase 2: Remove the runtime file-writing logic and the matching preparation-side debug attributes.
- [x] Phase 3: Verify no references remain and rebuild the affected target.

## Known Errors / Blockers
- No remaining blocker for this cleanup task.

## Review
- Removed runtime file emission from [martini.cpp](/Users/yinhan/Documents/upside2-md-martini/src/martini.cpp):
  - deleted `all_splines.txt`, `bond_splines.txt`, `angle_splines.txt`, `dihedral_splines.txt`, and `force_debug.txt` write paths;
  - removed the associated `debug_mode`, `force_debug_mode`, `overwrite_spline_tables`, and file-stream state from the MARTINI runtime classes.
- Removed preparation-side hooks from [prepare_system_lib.py](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/prepare_system_lib.py):
  - no longer writes `debug_mode`, `force_debug_mode`, or `overwrite_spline_tables` attrs into generated stage files.
- Removed the now-dead workflow env export from [run_sim_1rkl.sh](/Users/yinhan/Documents/upside2-md-martini/example/16.MARTINI/run_sim_1rkl.sh).
- Verification passed:
  - `rg -n "all_splines\.txt|angle_splines\.txt|bond_splines\.txt|dihedral_splines\.txt|force_debug\.txt|force_debug_mode|overwrite_spline_tables|UPSIDE_OVERWRITE_SPLINES|debug_mode = 1|_v_attrs\.debug_mode" -S src example/16.MARTINI`
  - `source .venv/bin/activate && source source.sh && cmake --build obj -j4`
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile example/16.MARTINI/prepare_system.py example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/extract_martini_vtf.py example/16.MARTINI/martinize.py`
  - `source .venv/bin/activate && source source.sh && bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `find . \( -name 'all_splines.txt' -o -name 'angle_splines.txt' -o -name 'bond_splines.txt' -o -name 'dihedral_splines.txt' -o -name 'force_debug.txt' \) | sort`
- Read-only merge audit for `origin/master -> martini`:
  - current branch is clean (`git status --short --branch` showed `## martini...origin/martini`);
  - merge base with `origin/master` is `94761377598955310efb435cfb3347984ce4a12e`;
  - `origin/master` has 10 commits not in `martini`;
  - overlap since the merge base is small: only `install_python_env.sh` and `py/upside_engine.py` changed on both sides.
- Read-only `git merge-tree` predicts two real content conflicts:
  - `install_python_env.sh`
  - `py/upside_engine.py`
- Suggested resolution strategy:
  - for `install_python_env.sh`, prefer the `origin/master` refactor because it already defaults to Python `3.11` and adds stricter environment setup checks;
  - for `py/upside_engine.py`, prefer the `origin/master` `load_upside_library()` helper over the branch’s simpler `platform.system()` branch selection.
- Manual merge implementation strategy:
  - because write-side `git merge` is prohibited by repo policy, reproduce the merge result directly in the working tree;
  - copy all 25 files changed on `origin/master` since merge base `94761377598955310efb435cfb3347984ce4a12e` into the current tree from `origin/master`;
  - for overlapping files (`install_python_env.sh`, `py/upside_engine.py`), keep the `origin/master` version exactly as requested;
  - verify those 25 files match `origin/master` after the sync, leaving the rest of the branch-specific MARTINI work untouched.
- Manual merge result:
  - synced all 25 `origin/master`-changed files into the working tree;
  - the two overlap/conflict files now use the `origin/master` versions:
    - [install_python_env.sh](/Users/yinhan/Documents/upside2-md-martini/install_python_env.sh)
    - [upside_engine.py](/Users/yinhan/Documents/upside2-md-martini/py/upside_engine.py)
  - direct blob/hash verification confirmed every master-changed file now matches `origin/master` exactly in the working tree.
