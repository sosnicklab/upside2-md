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
