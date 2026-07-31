# Progress log

## Current phase — production runs live (2026-07-31)

The dry-MARTINI hybrid prep was fully redesigned (CHARMM-GUI bilayer + martinize CG envelope + OPM
orientation + xy-only MC-barostat NPT + solvent-volume ions); see `plan.md` and
`example/16.MARTINI/readme.md`. Both HDX and NP-footprinting production are now running.

### glpG-DDM HDX REMD — 4 variants, RUNNING, healthy
Jobs 52832557 (79HIS), 52832591 (79HIS_S115T), 52832592 (79ALA), 52832593 (79ALA_S115T) on midway3
caslake, in `/home/yinhanw/project/glpG_DDM_REMD/`, calling the beagle3 binary. 48-replica T-REMD,
ladder `linspace(√0.70,√0.90,48)²`, self-resubmitting. Verified: all folded (hbonds 140–243), finite,
0 NaN/abort, block 1 in progress. Rare protein_potential transients (≤0.004% of frames, recover) are
negligible for MBAR. Next: HDX ΔG pipeline once trajectories accumulate.

### NP footprinting (1AO6 on 5 nm MPA-AuNP) — 6 faces, RUNNING, healthy
Job 52852032 in `/home/yinhanw/project/NP-1AO6/prod/`. Re-run after the dt=0.009 run went NaN
(backbone tore when a segment crept into the gold LJ core). Fixes: dt 0.009→0.005, graded soft-core
LJ/coulomb equilibration, 0.15 M salt over solvent volume (verified 0.150 M all 6 faces), gold NP
centered, and **box sized from the gold-centered reach** (fixed a 180-0-0 PBC wrap-split: box was
sized from the frame bbox but centered on the gold, so the longest-z orientation straddled the −z
edge). Calibration passed clean on all 6 (Rg 25.9, protein_potential −200 to −238, KE/1.5kT ≈ 1.01–1.03,
no NaN); now in production chunk 1. 6 independent trajectories, single T, no REMD.

## Done this phase
- `py/martini_prepare_system{,_lib}.py`: CHARMM-GUI parser, martinize CG, OPM Kabsch orient,
  solvent-volume ions, MC-barostat (xy-only, tensionless) stage config. `src/box.cpp/.h`: MC-only
  barostat (Berendsen/PR + virial removed), rebuilt clean.
- `example/16.MARTINI/readme.md`: rewritten as one coherent end-to-end HDX guide
  (prepare seed → T-REMD production → MBAR analysis).
- `scratchpad/NP-footprinting/np_hybrid.py`: reach-based box sizing; 6 faces rebuilt + uploaded.
- Stray files moved to scratchpad / deleted per user request.

## Notes carried forward
- Git strictly read-only; all edits uncommitted. Master (`upside2-md-master`) parity preserved.
- midway3 default shell is zsh (1-indexed arrays) — use explicit index mapping in upload loops.
- beagle3 binary needs `module load hdf5/1.14.3` at runtime (libhdf5.so.310).
