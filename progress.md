# Progress log

## 2026-09-01: current state

**Poster delivered 2026-08-19** (ver5, OneDrive). glpG wildtype (79HIS) ΔG-vs-residue figure complete: 22 836 frames/replica, 123/203 amides resolved at 298 K, Spearman 0.676 vs implicit membrane (pooled 4 repeats, 170 amides).

**midway2 REMD (48890657–60, all four variants):** block 1 running.
- 79HIS: 74%, 4 rollbacks
- 79HIS_S115T: 88%, 5 rollbacks; **active NaN in replica 27** — watch
- 79ALA: 45%, 10 rollbacks (high)
- 79ALA_S115T: 96%, 5 rollbacks

**midway3 REMD block 2 (56105310–13):** still PENDING (Priority ~130k).

**NP footprinting (six K190-proximal orientations, block 3 running):** post-fix block-2 measurement (18 246 frames, 3.2% adsorbed-and-compact) shows none of Carlson et al.'s five target lysines contacted (K190=0.000, K525=0.000, K541=0.033). Simulation CONTRADICTS paper's central claim (K190 most protected). No footprint npz exists for current block-3 data; run `np_footprint.py` after block 3 completes. Rg reaches 230.9 Å on run.3 — protein over-unfolds.

**Before any cluster HDX analysis:** re-upload `calc_hdx_ht.py` and `4.calc_D_uptake.py` from repo; cluster copy drifts from repo.

## Key completed milestones

| Date | Milestone |
|------|-----------|
| 2026-08-10 | findings-88 fix (BB force path) deployed; all seeds rebuilt |
| 2026-08-13 | MBAR reference subtraction fixed (findings 91); dG plots delivered |
| 2026-08-15 | BB proxy reworked onto `infer_H_O`; CB placement bug fixed (findings 102) |
| 2026-08-16 | Cluster rigid-protein bug found and fixed (current_stage must be "production") |
| 2026-08-17 | Membrane accessibility term wired into HDX pipeline (findings 113) |
| 2026-08-18 | NP campaign rebuilt with corrected CB + r_min_ang; relaunched |
| 2026-08-19 | Poster finalised; NP block-2 re-measured; all claims corrected vs Carlson et al. |
