# Progress Log

## Current Task: dryMARTINI Workflow-Path Cleanup

- Actions taken:
  - Traced the exact dryMARTINI path used by
    `run_sim_1afo.sh`,
    `run_sim_1rkl.sh`,
    `run_sim_1afo_full.sh`,
    and `run_sim_1rkl_full.sh`
    through `example/16.MARTINI/run_sim_hybrid.sh`,
    `py/martini_prepare_system.py`,
    `py/martini_prepare_system_lib.py`,
    and `src/martini_cg_lipid.cpp`.
  - Confirmed the master repository does not carry matching MARTINI workflow or
    MARTINI-specific `py/` / `src/` files, so the cleanup scope was defined by
    reachable workflow behavior rather than by a one-to-one master file list.
  - Identified that the user-reported
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_1rkl_hybrid/*stage_7.0*`
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/outputs/martini_1afo_hybrid/*stage_7.0*`
    files were generated on Jul 4 before the current injector repair.
    Those prepared files still carry the earlier implicit SC/target contract.
  - Restored the current workflow-used injector/runtime path to the last
    committed behavior for SC and target compaction correction:
    explicit `cgl_compaction_state` arguments with `delta_extended` /
    `delta_compact` datasets, while keeping `cg_lipid_pair` on the committed
    explicit compaction-state path.
  - Built controlled continuation sources from the reported stage-7 outputs and
    reinjected the repaired CGL nodes for direct A/B verification against the
    stale implicit-source artifacts.

- Files modified:
  - `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_hybrid.sh`
  - `/Users/yinhan/Documents/upside2-md/py/martini_build_tables.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_gen_params.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system.py`
  - `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system_lib.py`
  - `/Users/yinhan/Documents/upside2-md/plan.md`
  - `/Users/yinhan/Documents/upside2-md/findings.md`
  - `/Users/yinhan/Documents/upside2-md/progress.md`

- Verification:
  - `python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py py/martini_build_tables.py` passed.
  - `make -C obj -j4` passed after restoring the runtime behavior.
  - Direct `git show HEAD:parameters/dryMARTINI/dopc.h5` comparison confirms
    the installed
    `/Users/yinhan/Documents/upside2-md/parameters/dryMARTINI/dopc.h5`
    is byte-for-byte identical to `HEAD` on the inspected compaction datasets.
  - Reinjection of the current
    `/Users/yinhan/Documents/upside2-md/py/martini_prepare_system_lib.py`
    into a reported stage-7 template now yields:
    `cg_lipid_pair.arguments = ["compose_vector6d", "cgl_compaction_state"]`,
    `cg_lipid_rotamer_sc.arguments = ["placement_fixed_point_vector_only", "compose_vector6d", "cgl_compaction_state"]`,
    and
    `cg_lipid_target.arguments = ["compose_vector6d", "pos", "cgl_compaction_state"]`.
  - Controlled 2,000-step stage-7.1 continuations from the user-reported Jul 4
    stage-7 sources completed for both proteins:
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/example/16.MARTINI/outputs/diag_1rkl_implicit_cont`,
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/example/16.MARTINI/outputs/diag_1rkl_explicit_cont`,
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/example/16.MARTINI/outputs/diag_1afo_implicit_cont`,
    and
    `/Users/yinhan/Documents/upside2-md/example/16.MARTINI/example/16.MARTINI/outputs/diag_1afo_explicit_cont`.
  - Measured end-to-end / tilt from the final continuation frames:
    `1rkl` source `34.71 A / 46.56 deg`,
    implicit-control continuation `33.57 A / 45.27 deg`,
    repaired explicit continuation `34.34 A / 46.64 deg`.
    `1afo` source chain-A/B `39.17 A / 18.17 deg` and `39.10 A / 16.35 deg`;
    implicit-control continuation `41.61 A / 19.93 deg` and `37.44 A / 14.81 deg`;
    repaired explicit continuation `41.96 A / 19.83 deg` and `38.16 A / 15.17 deg`.

- Failures and fixes:
  - The Jul 4 `stage_7.0` artifacts the user pointed to were generated before
    the current injector repair, so they could not be used as proof of the
    current source behavior without a controlled rerun.
  - The tracked correction tables in `parameters/dryMARTINI/dopc.h5` did not
    regress; the repair work was isolated to the injector/runtime path and
    verified from the same saved stage-7 sources.
