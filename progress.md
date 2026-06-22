# Progress Log

## Current Task: Coarse Salt/CGL Z Fix And Diagnostic Cleanup

- Actions taken:
  - Initialized a focused task plan for the salt-ion, CGL z-coordinate, and
    diagnostic cleanup issues.
  - Restored the coarse workflow default to include explicit ions.
  - Changed default CGL initial z conditioning so CGL particles remain at the
    source DOPC geometric COM z unless an explicit z-separation override is
    requested.
  - Removed untracked diagnostic scripts not referenced by the four requested
    workflows or shared preparation path.
- Files modified:
  - `example/16.MARTINI/run_sim_hybrid.sh`
  - `py/martini_prepare_system_lib.py`
  - `plan.md`
  - `progress.md`
- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py
    py/martini_prepare_system_lib.py` passed after sourcing `.venv` and
    `source.sh`.
  - `bash -n` passed for `run_sim_hybrid.sh`, `run_sim_1rkl.sh`,
    `run_sim_1rkl_full.sh`, `run_sim_1afo.sh`, and `run_sim_1afo_full.sh`.
  - Targeted 1RKL preparation under `/tmp/upside2_md_verify_salt_cgl` produced
    `48` NA atoms and `45` CL atoms with explicit ions enabled; direct CGL
    collapse produced `276` CGL particles for `276` DOPC lipids, with max
    absolute CGL z minus source DOPC COM z `1.4e-14 A`.
  - Targeted 1AFO preparation under `/tmp/upside2_md_verify_salt_cgl_1afo`
    produced `47` NA atoms and `51` CL atoms with explicit ions enabled;
    direct CGL collapse produced `176` CGL particles for `176` DOPC lipids,
    with max absolute CGL z minus source DOPC COM z `2.8e-14 A`.
  - A full stage conversion was intentionally stopped because CGL
    target-clearance conditioning is expensive for a quick verification; the
    lower-level prep/collapse check covers the reported salt and initial z
    placement bugs.
