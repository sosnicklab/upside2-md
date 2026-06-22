# Project Goal

Fix the coarse hybrid 1AFO/1RKL workflows so their outputs retain salt ions
and place CGL particles on the same lipid-center z convention as the matching
full-resolution workflows. Remove diagnostic scripts that are not required by
`run_sim_1afo_full.sh`, `run_sim_1afo.sh`, `run_sim_1rkl_full.sh`, or
`run_sim_1rkl.sh`.

# Architecture & Key Decisions

- Preserve all hybrid SC-env and BB-env interactions; do not disable or scale
  physical interactions to avoid workflow problems.
- Keep edits scoped to the MARTINI example workflow and preparation code unless
  evidence points elsewhere.
- Treat the full-resolution scripts as the reference for salt retention and
  lipid COM z placement.
- Do not use git write operations. All edits remain unstaged.

# Execution Phases

- [x] Inspect current workflow scripts, shared wrappers, and local status.
- [x] Identify why coarse 1AFO/1RKL lose salt ions while full workflows retain
  them.
- [x] Identify why coarse CGL z coordinates diverge from full DOPC COM z
  placement.
- [x] Implement the smallest workflow/preparation fix.
- [x] Remove unused diagnostic scripts not needed by the four requested
  workflows.
- [x] Verify syntax and run targeted output/preparation checks.

# Known Errors / Blockers

- None.

# Review

- Root cause for missing salts: the coarse branch of
  `run_sim_hybrid.sh` defaulted `EXPLICIT_IONS=0`, while full-resolution runs
  retained the preparer default of explicit ions enabled.
- Root cause for CGL z mismatch: CGL particles were initially placed at DOPC
  geometric COMs, then an automatic leaflet z-separation conditioning step
  shifted the upper and lower leaflets apart. Full-resolution DOPC does not
  apply this shift.
- Fix applied: coarse hybrid workflows now default `EXPLICIT_IONS=1`, and CGL
  z-separation conditioning defaults to `0.0` so source DOPC COM z placement is
  preserved unless explicitly overridden by `UPSIDE_CG_LIPID_MIN_LEAFLET_Z_SEP`.
- Cleanup applied: removed untracked diagnostic scripts that are not referenced
  by `run_sim_1afo_full.sh`, `run_sim_1afo.sh`, `run_sim_1rkl_full.sh`,
  `run_sim_1rkl.sh`, their shared wrapper, or `py/`.
- Verification: Python and shell syntax checks passed. Targeted coarse prep
  produced salts for both reported systems: 1RKL `48` NA / `45` CL and 1AFO
  `47` NA / `51` CL. CGL z matched source DOPC COM z with max absolute
  difference `1.4e-14 A` for 1RKL and `2.8e-14 A` for 1AFO.
