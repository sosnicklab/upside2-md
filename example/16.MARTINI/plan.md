# Plan — membrane-aware HDX protection criterion (item 3) + workflow speed

## Goal
Codex HDX-trust checklist, item 3 only (user-selected scope): implement the
"open AND water-accessible" protection criterion for the Upside/dry-MARTINI
hybrid and wire it into the hybrid so per-residue open probabilities can be
computed. Verify on existing trajectories. (Cluster campaigns + experimental
benchmark = out of scope, user drives.)

Separately: investigate why example/16.MARTINI workflows are slow (delegated to
a background profiling agent).

## Key facts (verified against real data)
- Hybrid `.up`: `/output/pos` = (nf,1,n_atom,3) in Angstrom; trajectory lives in
  the run file (not loadable via mdtraj_upside, which assumes 3*n_res atoms).
- `get_protection_state.py` CANNOT run on the hybrid: it reads
  `sigmoid_coupling_environment`/`hbbb_coverage`/`environment_coverage_hb`/`surface`
  nodes that the hybrid config does NOT have, and feeds backbone-only coords.
- Hybrid config DOES have: `protein_hbond` (col 6 = backbone H-bond score, first
  n_donor rows = donors), `infer_H_O.donors` (residue + id; id[:,1] = amide N atom),
  `particle_class` (PROTEIN/LIPID/ION), `atom_names` (incl. lipid `PO4`).
- Midplane: PO4 z is bimodal → midplane = mean(PO4_z); leaflet planes = mean of
  PO4 above/below midplane. (median is wrong — lands in one leaflet.)

## Criterion (manuscript §4.8)
PS = 1 (protected) iff  backbone-H-bonded (score > hbond_cut)  OR  lipid-buried
(amide N between the two phosphate planes). Exchange-competent (PS=0) iff
open AND water-accessible (N beyond the phosphate plane, in the aqueous region).
Reduces to H-bond-only when no lipids present. 1 = protected (matches PS.npy
convention → drop-in for HXMS steps 4-6).

## Steps
- [x] Investigate hybrid trajectory/config structure + signal availability
- [x] Write `py/martini_protection_state.py` (hybrid-aware extractor)
- [x] Verify on clean stage-6.6 1RKL + 1AFO trajectories:
      1RKL (fully embedded helix) open 65.5% -> 0.0% (all rescued);
      1AFO (mixed) 25.4% -> 13.4% (~12% rescued, ~13% still exchange-competent).
      Mechanics verified; NOT dynamic convergence (only near-rigid eq frames exist).
- [x] Report + hand back speed findings

## Production stage-7.0 blow-up — DIAGNOSED (corrected per user)
Soft->hard table progression across 6.0-6.6 is BY DESIGN (gradual minimization),
NOT a bug. Confirmed: spline tables (combined_energy_grids) harden 6.0/6.2(63.7)
-> 6.3(689) -> 6.4/6.5/6.6/7.0(12583 @r=3A, ALL byte-identical). So the ladder
completes at 6.4 and there is NO hand-off cliff at 6.6->7.0 (same table). The
high energy under the harder table is expected and meant to be minimized away.

ACTUAL root cause = hard-core + g-JF BEAD-EJECTION when the raw core first engages.
Stage 6.4 is the first full-hard-table stage; its MD frame 0 is clean (68 A) but by
frame 1 lipid tail beads launch to ~8.4e5 A (runaway over the first frame interval).
Stages 6.2-6.6 run run_md_stage (MD), NOT run_minimization_stage - only 6.0/6.1
minimize. So the hard table engages at 6.4 with no minimization against it first;
contacts benign under 6.3's softer table sit on the steep r^-12 wall, and the g-JF
step gives a runaway displacement that ejects beads. Ejected beads PBC-wrap onto
neighbors and accumulate through 6.5/6.6, so stage-7.0 (40000-step MD) inherits a
corrupted lipid system -> initial E 9.2e13, stage-7.0 min line-search fails ->
production diverges/SIGABRT. Same instability class as the production divergence.
Lead: the current g-JF (martini_brownian.cpp) has NO displacement cap; the earlier
implicit Brownian integrator had a "displacement-cap safety net" (memory 2026-07-16)
that would have caught exactly this - the manuscript removed it ("no displacement
cap"). Tension with the "no arbitrary capping" rule.

Fix directions (physics-preserving; NOT implemented - needs user go-ahead):
 (a) Run a minimization against the hard table at the hardening step (6.4) before
     MD, so hard-core contacts relax before dynamics (currently 6.4 goes straight
     to MD).
 (b) Finer hardness ladder near the soft->full-hard jump (6.3->6.4 is 689->12583).
 (c) Make the g-JF step robust to hard-core force spikes when the raw core engages
     (revisit the removed displacement-cap safety net vs the no-cap rule).

## Script cleanup (py/martini_*.py) — RESULT: nothing to delete
7 files. 6 are LIVE/library (prepare_system, prepare_system_lib, itp_reader,
build_tables, gen_params, extract_vtf) - all referenced by run.py/run_sim_hybrid.sh
(git-grep confirmed). The 7th (martini_protection_state.py) is the item-3 tool just
created this turn (untracked, intentionally standalone like get_protection_state.py).
No orphaned/dead scripts exist; deleting any would break the workflow.
