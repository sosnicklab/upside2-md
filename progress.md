# Progress Log

## 2026-05-24 MARTINI Secondary-Structure Divergence
- Actions taken:
  - Compared full-resolution and single-particle lipid stage-7 inputs/logs.
  - Found that coarse mode kept CGL-SC active, but evaluated it as a standalone
    CB potential rather than through the sidechain rotamer one-body solver.
  - Added a rotamer-coupled CGL-SC node, `cg_lipid_rotamer_sc`, and made CG
    lipid node injection refresh generated `cg_lipid_*` nodes to avoid stale
    standalone SC terms.
  - Updated the CG lipid derivation document to describe the runtime ownership
    change.
- Files modified:
  - `plan.md`
  - `findings.md`
  - `progress.md`
  - `py/martini_prepare_system.py`
  - `py/martini_prepare_system_lib.py`
  - `src/martini_cg_lipid.cpp`
  - `example/16.MARTINI/cg_lipid_potentials.tex`
- Test results:
  - `source .venv/bin/activate && source source.sh && python3 -m py_compile py/martini_prepare_system.py py/martini_prepare_system_lib.py`
    passed.
  - `source .venv/bin/activate && source source.sh && cmake --build obj`
    passed.
  - Copied 1RKL stage-7 injection generated `cg_lipid_rotamer_sc`, removed the
    old `cg_lipid_sc`, matched 117/117 rotamer rows, and connected 282 CGL
    particles.
  - One-step `obj/upside` smoke test on `/private/tmp/cgl_sc_1body_test.up`
    passed; output contained nonzero `rotamer_1body_energy2` CGL-SC entries.
- Failures and fixes:
  - Initial smoke test failed because `cg_lipid_sc_1body` collided with the
    legacy registered prefix `cg_lipid_sc`; renamed the new node to
    `cg_lipid_rotamer_sc`.
