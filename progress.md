# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`. Entries before
10-24 are condensed by line of work; the full log to 10-27 is in git history (commit `c526abc`).

---

## Summary to 2026-10-06

* **09-18 to 09-24, trainer.** The ConDiv port then in use trained FF1's form: spline burial, no
  backbone desolvation, no DSE (findings 9t, 9u). It was replaced by Kleinmann's port of Peng's FF2
  trainer, restored to the SI. From ff2.1, 8 of 9 groups held a fixed point (9v). AWH on 40 glycine
  dipeptides finished, and ff14SB agrees with ff99SB-ILDN to 0.045 nats (9r).
* **09-27 to 09-30, round 2 (basin offsets).** ff_3.0 was released and withdrawn on 09-30 (findings
  1.14-1.15). Its glpG chains ran on seeds left with the FF1-form ff_3.0's coverage tables (3.11).
* **10-01 to 10-04, round 3 (`ff30_glyhb`).** AWH glycine map, glycine H-bond offsets in
  `HBondEnergy` (12-entry configs bitwise unchanged), side-chain rate / 10, the selection panel and
  release by checkpoint selection. The panels lost helix and the offsets went left-handed (1.17);
  stopped 10-04.
* **10-02, cleanup.** `training/` was cut to five files and made site-neutral. `findings.md` was
  condensed, and CLAUDE.md's compile and Slurm sections were rewritten (user).
* **10-02 to 10-05, rotamer-BP validation figure for Tobin.** Finished; it is in
  `/beagle3/.../bp_validation/static/`. The deck is not sent.
* **10-04/05, round 4 designed.** BioEmu-fitted glycine library (`rama31.dat`; push probe d -0.014
  [-0.041, +0.012]), non-glycine maps kept NDRD (1.19); ff30_bio and ff30_gdepth started on caslake
  with automatic validation; the poly-Gly reference started (Phase 12).
* **10-06, restart at dt 0.009** (user). dt 0.015 destroyed free replicas, and both epoch-0
  panels lost folding and helix (1.21-1.22). The dt 0.009 runs and the matched control
  ff21_ctrl_dt009 started on midway2. TM4 moved to 12 seeds. Those counts were later found invalid
  (3.11).

## Summary 2026-10-07 to 2026-10-23

* **Round-4 ff3.0 training, 10-07 to 10-12 (plan.md Phase 11).** Two arms were added 10-07 (user),
  every run starting from force fields byte-identical to its dt 0.009 twin's: SI rate (`ff30_bio_si`,
  `ff30_gdepth_si`, every group but rot at the SI's learning rate) and frozen H-bond (`hb`, `dhb`,
  `hbg`, `sheet` at rate 0: `ff30_bio_fz`, `ff30_gdepth_fz`, control `ff21_ctrl_fz`). How Jumper and
  Peng trained the H-bond energy is findings 1.23. Three runs converged at step 76 and were released:
  `ff30_gdepth_dt009` as ff_3.0_gdepth on 10-09 (1.26; copied to `parameters/ff_3.0/`), `ff30_bio_si`
  as ff_3.0_bio_si and `ff30_bio_dt009` as ff_3.0_bio on 10-10 (1.31, 1.32). Epoch-end panels: 1.21.
  The user stopped the frozen arms (10-09), dropped the BioEmu-map line for the basin-depth path (10-11,
  "as it is expandable to all residues") and stopped `ff30_gdepth_si` before its gate (10-12). Each was
  deleted with its derived products (manifests in `checks/disk_reclaim_20261011/`; lesson 10.25), which
  leaves ff_3.0_gdepth as ff3.0.
* **The local glpG TM4 test, 10-07 to 10-11 (findings 1.24-1.25).** `patch_glpg.py` never wrote the
  two coverage tables, so every TM4 count before 10-07 ran on the retired FF1-form ff_3.0's (3.11);
  fixed and deployed with the user's approval. User correction: TM4 is judged only in the hybrid,
  with the in-training force field (10.19). On fixed inputs each epoch end and release ran 12 or 24
  seeds in a local queue (`tm4_local/`, `tm4_queue.sh`, `tm4_compare.py`; sets and tables in
  `tm4_local/runs_cov/`); from 10-09 the primary test is last-block DSSP alpha of 135-151 (user;
  10.20). No in-training checkpoint was resolved from its start, and by resampling a true 0.10 / 0.15
  difference is detected 18% / 31% of the time at 12 seeds, 34% / 58% at 24. Retired 10-11 (user:
  "that metric no longer matters"). 10-15: ff2.1's TM4 unfolding in the user's figure is the old
  table (lesson 10.27).
* **Operations, 10-08 to 10-13.** The midway2 node exclusion law (user, 10-08):
  `/project/trsosnic/yinhan/slurm/midway2.args` is the one list, with `README.md`, `update_pending.sh`
  and `sbatch_wrapper` (`~/bin/sbatch`); remote_jobs.md §0d and CLAUDE.md carry it. midway2-0027,
  0088, 0103, 0116 and 0151 were excluded on their failure records. The Mac Studio watched the jobs by
  hourly cron with a cluster `WATCH_STATUS.md`; a MacBook Pro session kept its records in separate files
  and merged them (10.22). 10-11: `/project/trsosnic` at 98%, 928 G freed at the user's word
  (`checks/disk_reclaim_20261011/`).
* **Poly-Gly reference (plan.md Phase 12), 10-07 to 10-11.** Four 7.5 nm dodecahedron starts
  (~28,950 atoms) in self-chaining production; stopped by the user at ~97 ns per replica once the
  table fix had repaired TM4's glycine flips. Checkpoints kept in `polygly/`.
* **glpG TM1 and the side-chain/lipid table (plan.md Phases 13, 15), 10-09 to 10-10.** TM helices
  across 9 sets x 12 seeds (findings 4.6) and the TM1 diagnosis (4.7). Upside's membrane H-bond term
  (option A, the user's choice) holds TM1 but tears TM4 on ff2.1 (4.8); the TM1 kicks and both tears
  start at side-chain/lipid table spikes from under-sampling (4.10). The converged table (600 x 12)
  holds TM1 with no training: last-block alpha 0.931 against d9_03's 0.685, p 0.002 (4.13). In
  production 10-10 (user): `py/martini_build_tables.py` samples 600 x 12 by default,
  `parameters/ff_2.1/martini.h5` rebuilt, `patch_glpg.py` copies its table into patched seeds, queued
  glpG seeds re-patched (`checks/sc_table_20261010/`). On it ff2.1 holds TM4 (0.944 against its
  old-table 0.610, p 0.002) and no tested candidate improves on ff2.1 (4.15). The membrane H-bond term
  was not adopted (user): with the table fixed its gain is below the spread of two table-only sets, and
  it imposes the bilayer geometry from outside. `py/martini_extract_vtf.py`'s stationary-atom bug fixed (3.8).
* **ff3.1, per-map basin depths (plan.md Phase 14), 10-09 to 10-13.** Research and probe (findings
  1.28-1.29); the user chose per-map depths under an L2 prior. Trainer `training/ff31/` (1.30) ran
  `ff31_depth_all` and `ff31_depth_data` from ff2.1 to step 76, watched by
  `training/ff31/steric_check.py`, added after Jumper's question (1.33). Review (1.34): `all` is worse
  than ff3.0, `data` ties it, and two trainings of the ff2.1 control differ by 0.14 folded.
  `ff_3.1_depth_data` was released 10-13 and validated on the Peng set and glpG (1.35). Retired 10-26
  by the merged ff3.1.
* **ff3.2 (plan.md Phase 17), 10-15 to 10-18.** ff31_depth_data's trainer started from ff3.0
  (`$P/training/ff32_depth_data`); stopped 10-18 in step 69 for the angular fix and rerun (Phase 21).
* **glpG HDX (plan.md Phase 16), 10-10 to 10-23.** Literature (findings 5.5): TM cores exchange
  slowly but finitely (PF ~10^5-10^7), and our binary `acc` makes lipid-shielded PF infinite by
  construction. The local TM4 readout on the fixed table (10-11) and the ff3.0 and ff3.1 REMD readouts
  (10-14, 10-17) came before the 10-17 placement fix (5.5); the ff2.1 fixed-table REMD was cancelled
  10-14 (user) before its HDX, and a 10-13 request to switch the readout to ff3.0 was withdrawn
  (lesson 10.26). The PI's H145A TM5/TM6 question answered from the hybrid and implicit runs (4.16;
  `~/Downloads/glpG_H145A_TM5_TM6/`, `~/Downloads/glpG_WT_implicit_ff2.1_TM5_TM6/`). On the
  placement-fixed ff3.0 chains (10-22) TM4's core 141-144 is censored at 298 K in every variant and
  block, and block 3 against block 2 differs by RMS 0.65-0.91 kcal/mol against jackknife sigma ~0.1
  (5.6; figures in `~/Downloads/glpG_HDX_ff3.0_blocks_20261022/`, lesson 10.36). 10-23: HDX labels are
  PDB - 1, missed in every TM4-vs-HXMS comparison since 10-11 (10.38); the recheck with A61's
  accessibility shift undone is in 5.5-5.6. The user then closed the TM4 HDX question, to be read on
  Upside 3's glpG run, and dropped HDX on ff3.1's repaired chains.
* **The hybrid side-chain placement defect (plan.md Phase 18), 10-17 to 10-18.** Every hybrid seed
  placed its side-chain beads at `rcf[li[li]]` (findings 4.17). Fixed in
  `py/martini_prepare_system_lib.py`; at the user's word the 8 running glpG chains were cancelled and
  all 10 rerun from repaired seeds (`checks/placefix_20261017/`); block 1 holds TM4 under ff3.0 and
  ff3.1 (1.35). 10-18: 1.1 T freed on `/project`, the defective chains' trajectories included
  (`checks/disk_reclaim_20261018/`).
* **Deterministic against probabilistic side chains (plan.md Phases 19-21, 25-26), 10-17 to 10-21.**
  A second panel of 89 mdCATH domains (`ff3_selection/`): ff_2.1_nobp folds 0.264 against
  ff_2.1_ctrl38's 0.365, p 7e-5 (findings 1.21). The protocol audit (Phase 20) found the comparison
  clean but at matched length, not at no-BP's fixed point, and found the side-chain pair term's
  angular gradient misattributed in master's engine (3.12). Fixed in `src/interaction_graph.h` and installed in `$P/obj` (Phase 21);
  no-BP and its control continued from step 38, ff3.1 and ff3.2 rerun from step 0 as `_angfix`
  copies. At the fixed points no-BP folds -0.111 [-0.160, -0.065], p 1e-5 (Phase 25, 1.36). On TM
  score (Phase 26, user): mdCATH -0.034 (p 7e-8), 44 domains -0.026 (p 0.007), Peng natives -0.029
  (p 0.18), de novo +0.010. Lesson 10.32.
* **Released ff2.1 and the retraining loss (plan.md Phases 27-28), 10-21.** Every force field this
  workflow trains from ff2.1 folds 0.13-0.23 less than released ff2.1 on the 44-domain panel, mostly
  from epoch 0, the loss spread over side-chain table, backbone burial and H-bond; at dt 0.009 released
  ff2.1 folds as at 0.015 (findings 1.37). Released ff2.1 was submitted on all three tests
  (`checks/ff21_released_20261021/`).
* **ff3.1.1 and ff3.2.1 (plan.md Phases 29-30), 10-21 to 10-23.** Mode `pooled` in
  `training/ff31/rama_basin.py`: one shared offset per central residue and basin plus per-map
  deviations, width 0.5 by cross-validation. At the user's word only the depths train
  (`ff311_depths_only`, `ff321_depths_only`, from NDRD's and ff3.0's maps). 10-23: glycine alpha_L's
  gap closed 0.517 to 0.179 and still closing, every other map agrees to the runs' noise (findings
  1.34); both continued 76 -> 114 under an agreement rule fixed before the data (Phase 30).
* **Large proteins de novo (plan.md Phase 22), 10-19.** On mdCATH the folded-fraction drop at 151-250
  residues is a metric artifact; local de novo runs of 155-179-residue proteins and the control end
  in compact non-native globules, and vectorial growth equals bulk (findings 11e, upside3.md 8.3.2).
  Study data in `scratchpad/large_protein_20261019/`.
* **The `environment_coverage_hb` typing defect, 10-19** (findings 3.13, upside3.md 4.1). Fixed in
  `src/environment.cpp` and `py/upside_config.py`, old configs bitwise equal; upside3.md written.
* **Upside 3 (plan.md Phases 23-24, upside3.md), 10-19 to 10-22.** The user made martini-dev Upside 3,
  with fixes there only (CLAUDE.md reworded). Planning: findings 11f (reports in
  `scratchpad/upside3_planning_20261019/`). The audit: a harness in `tests/` (`fixtures.py`,
  `node_fd.py`, `regress.py`, `config_diff.py`, `registry.py`, `index_oracle.py`, `line_integral.py`,
  `invariance.py`), pass 1 and pass 2 of every unit by reviewers (`pass2/*/report.md`), every fix
  reproduced before and verified after; ledger A1-A233 in findings 3.14. Among them: REMD and the
  thermostat (A3-A28, A143-A149), the dry-MARTINI tables and barostat with `martini.h5` rebuilt
  (A74-A76), the side-chain table rebuilt from ideal chi-row conformers (A192), the hybrid preparation
  (A70-A73, A83, A183-A191), torch out of the glue and trainers (A204), the sequence-distance exclusion
  key (A231). Bilayers: CHARMM-GUI step5 starts interdigitated and GROMACS shows the same early box
  jump, so preps start from GROMACS-equilibrated membranes (findings 8); glpG's solvation gate (8b,
  A102, A126). The closing criteria (upside3.md 4.3) were met 10-22. The Upside 3 tree
  `/beagle3/trsosnic/yinhan/upside3` was built on midway2 and named in CLAUDE.md; Peng's 32 arms and
  the glpG prep and REMD were submitted as the validation of the debug, through
  `example/16.MARTINI/run.py` rewritten to submit stages as jobs and REMD in blocks (A232-A233).
  Lessons 10.29-10.31, 10.33, 10.35.
* **The Boltz-2 HDX sampler (diffusion_dynamics.md Phases 0-1), 10-21 to 10-23.** User: a modified
  Boltz to replace Upside for HDX sampling, built on Boltz, never BioEmu (lesson 10.34; BioEmu2 code
  allowed). Installed on `/beagle3` (Boltz at `b1ebfc4`, torch 2.14.1+cu126), run on midway3 `gpu`.
  Phase 0: `step_scale` 1.0 causes the bond-length scatter and the `md` tag widens the ensemble; Phase 1
  uses `md` at 1.0. Our clone fixes Boltz's sample chunking, and the wrappers fail a skipped batch. At
  10^3 samples every measured amide but one lies above the bound RT ln(1000) = 3.75-4.09 kcal/mol
  (findings 11b); a Ub psWT set (like for like with the experiment) and HEEH_rd4_0097 at 10^4 (rank
  order) were added.
  Cluster lessons in 10.14.
* **Records and machines.** Handoffs for cold starts on another computer on 10-07, 10-14 and 10-18
  (remote_jobs.md "Resume here"); the user's pushes merged on request with local edits unstaged;
  scratchpad cleanup 10-15, 224 GB (`scratchpad/CLEANUP_20261015_manifest.txt`); the plotting venv
  rebuilt per `manage.sh` 10-23 (Python 3.15.0, matplotlib 3.11.2).

## 2026-10-24

* **Status and the Boltz Phase 1 readout** (user: "check status", "proceed to next step for Boltz
  job"). Logs clean. HEEH at 10^4 samples, 278 K: halves agree (Fisher p >= 0.125), 4 of 20 measured
  amides never open where experiment has 3.2-3.8, Spearman 0.04 against experiment (Upside FF2 0.32);
  Ub psWT below the censored experiment at 1 of 24 (Upside 7). Fixed the name clash that had crashed
  `example/00.AnalysisScripts/helpers/plot_ref_style.py` on every input since 10-21.
* **Boltz-2 at 10^5 samples** (user: "Submit the jobs for next steps"). HEEH's hairpin opens 25-35x
  less than in experiment, so EHEE and Ub psWT went to 10^5 (`phase1_seeds.sbatch`, seeds 13-112, an
  hourly submitter for the gpu QOS's 12-job limit). Boltz's 5 files per sample filled the group's
  `/beagle3` inode quota by 12:20 and all 20 jobs failed (lesson 10.39). Fix (user: "Fix the problem
  for me"): `pack_predictions.sh` tars each predictions directory once the tar lists every file,
  `check_samples.py` reads the tars, each seed is packed as it ends; resubmitted. GNU tar's "file
  changed as we read it" exit on complete GPFS files now passes on to `tar --compare` and the file
  count. Tools in `boltz_hdx/` (scratchpad and cluster).
* **BioEmu2 HDX side test set up** (user chose "Lift it, side test only" for the 10-21 rule). The
  exception is in findings 10.34, the memory and diffusion_dynamics.md sections 1 and 7. midway3
  `/project/trsosnic/yinhan/bioemu2_hdx`, venv on uv's CPython 3.13.17 (the `python/3.13.11` module
  lacks `_ctypes`); one V100 job for EHEE and Ub psWT at 10^3 and HEEH at 10^4.
* **ff3.1.1 / ff3.2.1.** ~33 min a step, KE/1.5kT at most 1.018, glycine alpha_L's gap 0.179 at round
  4. The collector's status line (user: "fix it") counts flagged steps from the newest steric file:
  ff3.1.1 has six above the threshold, not the two first reported, all minibatch 1 or 14 (findings
  1.33).
* **Watch started 23:45** (user: "continue monitoring midway2 and midway3"). The first pass found
  inode-exhaustion casualties the 13:10 check missed: Upside 3's REMD block 1 ended with its
  `block_count` write refused, so block 2 never queued (48 replica files intact), and four Peng chunks
  failed in 9 s. Both resubmitted at the user's word. New read-only probe `~/watch_probe/watch_mdw2.sh`.

## 2026-10-25

* **Watch passes overnight.** The EHEE and Ub psWT 10^5 sets and the first BioEmu2 job finished,
  every seed `ok`. ff3.1.1 reached step 114. New flagged steric steps 109 (ff3.1.1) and 90 (ff3.2.1),
  both minibatch 14, whose restrained runs read as high (findings 1.33). The trainer's round line
  counts "|c| > 0.1" over every offset against the trainable depths only, so it can exceed its total.
* **Boltz-2 and BioEmu2 readouts** (user: "Work on boltz and bioemu2 job"). `pool_seeds.py` pools each
  protein's 10 jobs; `samples_to_traj.py` and `dg_phase1.py` are importable, one readout for every set,
  old outputs reproduced byte for byte. Boltz-2 at 10^5: EHEE ranked (Spearman 0.74, bootstrap
  0.51-0.86; Upside FF2 0.61), Ub psWT not (-0.21; Upside 0.43); the openings are rare large
  excursions, not broken chains. BioEmu2: HEEH Spearman 0.24 (Upside 0.32), its unfolded state sampled
  at ~4e-4 but the hairpin overprotected; EHEE over-open at 9 of 15 amides at 10^3.
* **Controlled 10^5 tests** (user: "propose and apply fix"; plan approved). From the raw flags,
  protection is burial-dominated: H-bonds break in 1-95% of samples at the measured amides, burial
  keeps them protected, so only rare excursions open them; seeds and BioEmu2 batches are exchangeable
  (dispersion ~1.0). Pre-registered and submitted on midway3: Boltz-2 Ub psWT with EHEE's query-only
  MSA (`inputs/ubqps_nomsa.yaml`; nine jobs held for a fold gate), BioEmu2 EHEE and Ub psWT at 10^5.
  `dg_phase1.py` gained a block bootstrap: HEEH -0.06 to 0.17 (Boltz-2), 0.14 to 0.35 (BioEmu2).

## 2026-10-26

* **00:50-01:20, the two questions** (user: "is ff3.1.1 and ff3.2.1 currently equivalent?", "check if
  current upside 3 result is consistent with upside 2 result"). Over rounds 0-6 (`dump_rama_state.py`,
  `compare_pooled_runs.py`): agreement (a) is met at round 6, (b) falls short only by a constant the
  data cannot move (findings 1.34). Upside 3's glpG block 1 against Upside 2 with the same code: KE,
  C-N and five helices match, TM4 is lost in prep stage 7.0 (3.15). 44 Peng scoring jobs submitted for
  the engine comparison.
* **09:40-10:55, merged ff3.1 released** (user: "merge into ff3.1 now and run validation. Remove any
  existing ff3.1, ff3.2 ..."; plan approved: retire the released sets and results, glpG on Upside 3
  only; plan.md Phase 30). Merge written and checked in `scratchpad/ff31_merge/`,
  released to four trees, the old ff3.1/ff3.2 sets and results moved aside
  (`checks/ff31_merged_20261026/release.sh`). Soluble validation submitted: 2 panels, 2 TM jobs, 32
  Peng arms, a collector. A234: the hybrid prep takes `--force-field` (for ff_2.1 the old and new
  injection are identical, for ff_3.1 only `rama_pot` differs), deployed to the Upside 3 tree; glpG
  with ff3.1 submitted there. Peng engine comparison: natives agree, the 6 finished de novo arms are
  lower on Upside 3. Handoff for "continue monitoring": remote_jobs.md "Resume here", "Handoff" and
  "The watch" rewritten, CLAUDE.md names the phrase, bundle `checks/handoff_20261026/` on the cluster.
* **12:25 on, the all-residue Ramachandran PDF** (user: every amino acid, ff2.1 against ff3.1, in
  `gly_rama_ff321_nonhelix_zero.png`'s form). `~/Downloads/rama_all_ff21_ff31.pdf`, 9 pages, one
  central residue per row (20 + cis-proline): the pooled coil map, ff2.1 (NDRD's X|right|ALL) and ff3.1
  (that plus the mean change of X's 39 pool maps), on the non-helical zero. Page 1 ranks each basin's
  depth change (blue deeper, red shallower) and overall the share of probability that changes basin
  (GLY 10.0%, GLN 8.1%, LYS 6.0%); page 2 the native training sites per basin (410 training proteins,
  46 of 456 held out, 44,330 sites), a basin with no trainable offset hollow. An untrained basin moves
  only with its residue's shared shift (alpha_L dG equal to pPII''s to 0.001), whose sign follows the
  summed free-minus-native gap of the untrained basins. Against the user's earlier glycine figures the
  arrays are pixel-identical (ff2.1 to `gly_rama_bias_nonhelix_zero.png`, the retired ff3.1 to
  `gly_rama_ff31_nonhelix_zero.png`), only the labels differ, and the merged ff3.1's alpha_L is -2.58
  against the retired one's -2.95. Open: on the old 37-map set glycine reads -2.568 / -1.504 / -0.575
  at the reference cells against findings 1.34's -2.574 / -1.510 / -0.579, not traced
  (`scratchpad/ff31_merge/` is on the other Mac). Scratchpad scripts: `ff31_training_counts_remote.py`,
  `ff31_training_gap_remote.py`, `rama_all_pooled_maps.py`, `plot_rama_all_pooled.py`.
* **13:55-14:25, MacBook Pro, Upside 3 glpG with ff3.0** (user: "upside2 did unfold TM4 ... before a
  bug fix", "do submit ff3.0"). The Upside 2 ff2.1 TM4 references all predate the 10-17 placement fix:
  findings 3.15 corrected, lesson 10.40. A235: the hybrid preparation's side-chain check compares only
  the three `sidechain.h5` datasets the side-chain/lipid table is built from, so ff3.0 runs (ff_2.1,
  ff_3.0, ff_3.1 pass; a changed chi row or bead count is refused); deployed to the Upside 3 tree and
  glpG with ff3.0 submitted (`upside3_validation/glpg_ff30`).
* **14:50-22:30, MacBook Pro, glpG in a DDM micelle (plan.md Phase 31)** (user: "charmm-gui fails to
  build micelle this size. Let's build it locally", "you should implement the missing dihedral").
  Micelle sizes from Oliver 2013 SAXS and glpG's measured belt (N 140 empty, 150-200 around glpG).
  A236: `dihedral_spring` (GROMACS type 1), the reader's dihedral key bug, ITPs under FLEXIBLE, the
  official DDM ITP; equal to GROMACS on 600 dihedrals, node_fd clean, regress 0 of 4260 arrays differ.
  GROMACS self-assembly of a frozen CG glpG with 200 DDM, 200 ns (`scratchpad/ddm_micelle/`): all bound,
  128 in the belt and 72 outside, maltose heads within 0.6 nm of 56% of the belt core's exposed beads
  (0% in POPE/POPG), so the complex is not used yet. 140 DDM alone form the SAXS micelle, so the head
  coat comes from the frozen Martini glpG's contact or the excess detergent (findings 7.5).
* **15:45-16:10, MacBook Pro, Python environment** (user: "transit everything to python 3.15 if
  possible"). Homebrew had removed python@3.11, the `.venv`'s interpreter; the dead venv is kept as
  `.venv.bak_py311_20261026` and `.venv` was rebuilt on Homebrew Python 3.15.1, with regress, node_fd,
  the GROMACS dihedral check and index_oracle as on 3.13. Replaced by 3.14 on 10-27.

## 2026-10-27

* **08:50-09:38, MacBook Pro, Boltz-2 MSA control and BioEmu2 EHEE at 10^5** (diffusion_dynamics.md
  Phase 1 and section 7). The flag code's scratch engine, rebuilt at `6d295b0`, reproduces the earlier
  flags byte for byte under Python 3.15; 40 job npz pooled in ledger order. MSA depth is ruled out on
  Ub psWT (censored Spearman -0.06, bootstrap -0.25 to 0.09). BioEmu2 opens EHEE through a 0.18%
  unfolded population, all 15 amides at 3.5-4.7 kcal/mol (censored Spearman 0.46; Upside FF2 0.61).
* **12:00-13:55, Python 3.14 on both Macs** (user: "transit to 3.14 then"; on the MacBook Pro, "Git
  pull and merge. Also check if we need to fix the venv", the read-only git rule lifted for that
  request). Python 3.15 has no wheels for torch (so pydssp), pyhdx 0.4.3 or scikit-image; 3.14 and 3.13
  install the full set. `.venv` on both machines is Homebrew Python 3.14.10, built by
  `install_python_env.sh`; old venvs kept as `.venv.bak_py311_20261027` (Mac Studio) and
  `.venv.bak_py315_20261027` (MacBook Pro). `pyhdx.plot` imports under no version (its proplot extra
  pins matplotlib 3.4), so 5.analyze_D_uptake.py keeps its own `apply_cmap`. On the Mac Studio
  `source.sh` was regenerated from `source_arm` and the engine rebuilt against `eigen@3` (backups
  `source.sh.bak_pre_py315_20261027`, `obj/*.bak_pre_rebuild_20261027`); the MacBook Pro fast-forwarded
  096486cb -> 3005a3fd. Verified on both: the 28 fixtures equal but for paths in `/input/args` and
  `/input/adv_args`, regress 0 of 4064 arrays differ; on the Mac Studio mdtraj DSSP, index oracle and
  node_fd output identical (the 1rkl hybrid fixtures, MacBook only, not rerun); on the MacBook Pro the
  BioEmu2 and Boltz HDX readouts byte for byte. Present under every stack, settled 10-29
  (below): index oracle flags on 6 datasets (`hb_surf_membrane_potential` and
  `hb_membrane_potential` index lengths in 1rkl_surface and 1rkl_channel_m0, `rama_map_pot/rama_pot`
  rows 1, 4, 8 in chig_rama); node_fd one bad entry in each of 6 checks over 1dfn_cavity, 1qhj_memb,
  1rkl_channel, 1ubq_rotamer, chig_hbbb. node_fd's probes come from one RNG across all configs, so its
  logs compare only over the same config list in the same order.
* **12:32-12:55, MacBook Pro, BioEmu2 Ub psWT at 10^5** (diffusion_dynamics.md section 7). 100,000
  samples pooled; seed 1700's one physicality-rejected frame has an ordinary backbone and is kept.
  BioEmu2 never leaves ubiquitin's native basin (CA-RMSD at most 4.17 A), so experiment's least stable
  amides V70, W45, I13 and L67 never open (6-12 expected each) and the order is lost (censored Spearman
  0.06, bootstrap -0.07 to 0.43; Upside FF2 0.43). With HEEH and EHEE, BioEmu2 at its defaults does not
  rank these proteins' HDX as well as Upside FF2; section 7's 10^5 item closed.
* **15:20-15:55, MacBook Pro, Gate 1** (diffusion_dynamics.md Phase 1). Phase 1's three small Boltz
  sample sets copied from midway3, so `run_flags.sh` runs end to end here (88 flag files and 24 dG
  outputs reproduce byte for byte). Peng 2022 SI Tables S4-S7, each amide's urea slope against the
  global m-value (`flags/urea_classes.py`): EHEE's measured amides subglobal, HEEH's helices global and
  its hairpin subglobal, Ub psWT's amides below the bound local, the four Boltz-2 never opens among
  them. Gate 1's HDX half recorded as "wrong states"; the risk that mdCATH lacks the openings now holds
  for global and subglobal ones only.
* **16:05-19:10, MacBook Pro, Phase 2: an add-on to frozen Boltz-2** (user: "You are not allowed to
  retrain boltz ... But you can train small models as addons"; diffusion_dynamics.md Phase 2). A global
  add-on cannot help: Boltz-2's ensembles mixed with experiment's own dG_global keep Spearman 0.74,
  0.04 and -0.19, so the errors sit inside the folded ensemble. Plan approved: a learned tilt potential
  sampled through Boltz-2's guidance and Feynman-Kac resampling, trained on Ferrari et al.'s multiplexed
  HDX-MS, gated by E1, a no-training test of whether steered Boltz-2 can estimate its own rare openings;
  the 10-21 Phases 2-3 dropped. E1 code in `scratchpad/boltz_hdx/addon/` (`sampler.py`, `opening.py`,
  `e1_rare.py`, `e1.sbatch`, `e1_readout.py`); parity passes in a CPU sandbox with the cluster's Boltz
  tree, a strides bug in `opening.py` fixed, the GPU pilot submitted on midway3. The first potential
  acted on Boltz's O atoms, not the estimator's inferred geometry, and the 10^5 flags show V70, I13, W45
  and L67 buried far past the open threshold in every sample. Lesson 10.41.
* **19:10-22:30, MacBook Pro, watch and handoff** (user: "Record this session's progress and I will
  quit and work on another computer"). ff3.1's mdCATH TM is level with released ff2.1 and above the
  control; ff3.1 Peng BBA de novo complete (final frames finite, longest C-N 1.68 A, KE/1.5kT
  0.999-1.008) and scored with BBL native (3 of 32). The Mac Studio session's submissions recorded
  (8 Peng seed replicates; ff3.0 glpG `n_blocks` 1 -> 3). The midway3 socket master ended about 21:40.
  Handoff: watch cron deleted, owner none, the watch prompt rewritten for the live jobs, the bundle
  given the MacBook Pro's records, working-tree patch, memory and git-ignored Phase 2 tools.
* **10-28 10:45-11:00, ff3.2's priority and the beta' fix** (user: "ff3.2 should be kept at a higher priority than
  ff3.1. Check if we have the beta issue addressed in ff3.2"). Both nativefix chains, and the debug's, had died
  10:41 in a Slurm controller outage. The step-0 gate, never run since submission, reads FAIL as written (glycine
  beta' 0.656 against 0.75; old target 0.356); `scratchpad/gate_ceiling.py` shows a 20 deg spring gives at most
  0.726 there, and the replica sits near that ceiling in every basin but non-glycine beta' (n 18). Links 49635669
  (ff3.2) and 49635327 (ff3.1, nice 10) held, then released 10:58 on "proceed with your best decision"; the debug
  chains restarted 10:58 on their own and were put at nice 10 behind ff3.2. Files: plan.md Phase 33, findings.md
  1.34, remote_jobs.md.
* **10-28 11:10-11:50, MacBook Pro, training terms for HDX, large proteins and IDRs** (after the PINN paper was
  judged no help; user: "the highest priority is HDX. The one following it is large protein performance, then by IDR
  sequences"; "check if you can do analysis with [midway2 results]. if not, write it for later plan"). Survey:
  ubiquitin's Peng-benchmark
  REMD is usable; no Upside EHEE/HEEH run, no RNase HI REMD, no IDP run, no kept ConDiv step trajectories, and the
  benchmark homeodomain ladder sits above its transition. M1 on the Mac (`scratchpad/hdx_term_20261028/`): the stock
  estimator on ff2.1 and ff3.0 native arms, MBAR to 273 and 300 K. ff2.1 reproduces Peng's FF2 (Spearman 0.90); 14-15
  of 24 measured amides resolve at 273 K (bound 5.6); ff2.1 opens the globally exchanging amides within the folded
  state; ff3.0's ubiquitin melts about 0.06 T_up lower, with 3 ladder round trips; the definition shift is 0.12
  kcal/mol. Files: findings.md 5.7, upside3.md 8.2 (design and steps H2-H5, then M2-M6 and T1) and 11 item 9, plan.md Phase 23
  item 5. No cluster job submitted; three arms copied to the Mac read-only. Then "record the plannings. We will run
  it later": M2, M3 (Peng's FF2 protocol from SI pp. 12-13, about 570 node-hours) and M6 (now H2, H3 and B0) specified in upside3.md
  8.2, the item added to remote_jobs.md "Resume here" step 7, `fetch.sh` added, and the handoff bundle's records
  and M1 tools refreshed (old records in `records.bak_20261028_1150/`).
* **10-28 11:05-13:00, MacBook Pro, disordered regions** (user: "for upside3, investigate the possibility of
  enhancing intrinsic disordered region simulation accuracy"; then "record everything in upside3.md. We will
  implement this later"). Code: training never sees a disordered chain; the DSE target is ConDiv's SARW (backbone
  excluded volume, Rama only). Literature by sub-agent (`scratchpad/idp_20261028/lit/review.md`): no Upside-IDP
  comparison with experiment exists; a 16-entry SAXS panel. User corrections: Upside's electrostatics are implicit in
  the trained terms (findings 10.42), and no explicit electrostatic term is added (double counting). B0 launched on
  the Mac 11:56 (ff2.1 x 6 and SARW x 3, then ff3.1 x 6; 30,000 tu from random coil); alpha-synuclein pilot: ff2.1
  35.6 and 31.0 A against 35.5 measured, SARW 40.1, Rg tau 340-1,300 tu. Files: upside3.md 8.4 (with B0, 11
  item 10, 12), plan.md Phase 34, findings.md 10.42.
* **10-28 16:10-19:50, MacBook Pro, B0 read** (stage A at 16:10, stage B at 19:37).
  * Results, ff2.1 / ff3.1:
    * Offset 0.886 / 0.850, with bootstrap intervals that exclude 1.
    * r of log(Rg/N^0.6) 0.46 / 0.49.
    * ACTR and A1-LCD spend long stretches as globules.
    * Cv78/Cv14 is 1.056 / 1.076 against 0.753.
    * 11 / 10 of 16 entries are unresolved.
    * ff3.1 differs from ff2.1 only for CoRNID (31 against 48 A) and Hst5.
    * The SARW offset is 1.20.
  * Correction: ACTR helix 1 seeds range 0.00-0.73. Stage A's "close to measured" came from the seed mean alone, and
    4.2 now records helix 1 as unresolved.
  * Fix: `stat_ineff.py` failed on HDF5 file locks, because it opened stage B files the driver was still writing.
    It now covers only the models in `results.txt`'s seeds line.
  * New: `read_b0.py`, which applies 4.2's rules, including per-seed helicity and globule fractions.
  * Recorded in upside3.md 8.4, 11 item 10 and 12, and plan.md Phase 34.
  * Then, on "record this ... I will work on it later in another session":
    * Added 4.2 steps B1 (ACTR at 298 K) and B2 (sampling).
    * Added a "Resume here" for a cold session.
    * Wrote B0's answers into I2 and I3.
    * Flagged that ACTR's SAXS source is unchecked. The literature sub-agent once sent the user's email in a Crossref
  User-Agent header, against its instructions.
* **10-28 23:40-10-29 00:20, Mac Studio, the Upside 3 harness on HEAD** (user: "read upside3 design, run any local test
  that is needed and you can run right now"). Code unchanged since A236; 28 fixtures regenerated, equal to 10-27's
  but for paths. Regress 0 of 4064 arrays against both 10-27 records. Invariance remd, thermostat, mc (pivot
  acceptance 0.7438) and chig equipartition PASS. Index oracle: the 10-27 flags were the oracle's, stale since A155
  and A168; fixed as A237 (`tests/index_oracle.py`), 0 mismatches on 28, three injected errors caught. node_fd's 6
  flagged directions line-integrate to noise (1ubq_rotamer's at width 5, below float32 at the default). The
  1rkl_channel one comes from the `surface` score, a count with no derivative: the channel energy steps with no
  force (up to 1.4 E_up per step along CB placement), recorded for the user (findings 3.14, upside3.md 11 item 15).
  Not run: the hybrid fixtures (CHARMM-GUI inputs and baseline h53 on the MacBook Pro only); 4.2 B1 and M2 (their
  scripts and data are on the MacBook Pro and the cluster). Files: `tests/index_oracle.py`, findings.md, plan.md
  Phase 24 item 2c, upside3.md 11.
* **10-29 08:06-09:45, Mac Studio, status pass and the owed items** (user: "check remote job status", then "fix them
  all"; both sockets reopened, two Duo pushes). Found `ff311_glypair` hung since 22:41 10-28 in step 34: a
  `compute_divergence` Pool child for 3g0m died of a glibc double free in `RotamerSidechain`'s destructor and
  `pool.map` waited for its lost task (remote_jobs.md §1 line 2, §8). Boltz E1 diagnosis 61755382 failed its rule 1
  (no reproduction at the same seed), so it was not read further. Then:
  - 3g0m's step files copied aside md5-checked (`checks/ff31_debug_20261027/hang_3g0m_step34/`), the hung srun step
    cancelled 09:20; the link ended FAILED and insurance 49636476 waits to rerun step 34.
  - ff3.2's epoch-1 panel `n32_01` submitted (49640339).
  - Boltz reproducibility job 61786565 (six repeats of one setting on one GPU: current, pre-trace and deterministic
    mode), its reading fixed first in diffusion_dynamics.md Phase 2 E1; `addon/e1_repro.sbatch`, `repro_compare.py`
    (dry-run on the existing npz). Read 10:15: no pair reproduces, deterministic mode included, parting at step 0;
    the repeat test must be statistical (`log_z` s.d. 1.6 over 8 runs), the user's choice. On the user's "go with
    your recommendation", the diagnosis's rules 2-3 read statistically (rule fixed first; `addon/diag_rules23.py`):
    no shrinking from 256 to 1,024 particles, pointing to a defect in the weights at low power; step 0 carries 26
    nats in every run, and the loss is in the potential increments, not the guidance ratio.
  - Round 2 read for ff3.2 and `_mb2` (plan.md Phases 32-33; `gly_LR` reproduced round 1's -0.567 and -1.065 first).
  - glpG on Upside 3 read: ff2.1 blocks 2-3 and block 1 of ff3.1 and ff3.0 (findings 3.15): TM4 0.44, 0.26, 0.14
    under ff2.1, 0.42 under ff3.1, 0.85 under ff3.0 (Upside 2 0.96).
  Files: remote_jobs.md, plan.md, findings.md, diffusion_dynamics.md, progress.md; scratchpad/boltz_hdx/addon.
* **10-29, Mac Studio, `upside3.md` restructured** (user: "read through upside3.md, do some research and formalize the
  design ... make it logical and ordered, and fill in the gap between steps"; "and run any local test that are
  needed"). Plan approved (`~/.claude/plans/mighty-prancing-candle.md`), then one rewrite: scope and definitions,
  standards for every step, a roadmap of stages S1-S7 and track G with gates, the audit, units and the training tree,
  a temperature section, a base-force-field section, the training changes (labels H1-H5, L1-L6 with D1-D4, B0-B2,
  I1-I3; GPU ladder V0-V3), validation and release, GPU, one decisions table, files. Research: two record-reading
  sub-agents, a literature sub-agent (Peng's S14-S16 calibration, HDX fitting precedents, Ferrari's per-domain data,
  how force fields accumulate changes) and an independent structural review. Local test: the channel potential in
  dynamics on 1rkl_channel (force-less energy 1.43 E_up over 9 tu against an s.d. of 1.66; xy translation
  invariant). Checks: every old number kept except stale code line numbers; every section reference resolves; no
  em dash or listed word. Cross-references updated in plan.md, findings.md, remote_jobs.md, progress.md, CLAUDE.md
  and two test headers. Files: `upside3.md`, findings.md (3.14 note, 5.7, 11b), plan.md Phase 23 item 6.
  Then the user's GPU ruling ("the same results as the current ffs, without retraining"): `upside3.md` 10 rewritten
  to it (the current force fields from the same files, to the audited CPU engine's results; energies and forces
  only, so V4 and the parameter-derivative check dropped; G3 and G4 no longer wait on S5 or the side-chain choice),
  with the roadmap, 1.2, 7.1, 9.1 and decision 12; plan.md Phase 23.
* **10-29 12:20-12:55, Mac Studio, status pass and the Boltz toy test** (user: "yes. and check status"). Status:
  nothing failed. ff3.2 at step 62 and `_mb2` at step 61, round 3 read for both (plan.md Phases 32-33; a forked
  sub-agent ran the readout, `gly_LR` reproduced round 2's -0.480 and -0.913 first); ff3.1's Peng arms 19 of 32
  complete (four since 08:10, checked finite, C-N 1.33-1.35 A, KE/1.5kT 0.985-1.026), 13 scored; Upside 3's cspA de
  novo complete (30 of 32); glpG ff2.1 block 4 about 43% through on block 3's pace; queued jobs unchanged. Found
  Upside 3's Peng rescoring not run since 10-26 (8 finished de novo arms, 16 jobs), proposed to the user. Toy test
  (rule fixed first in diffusion_dynamics.md Phase 2 E1): `addon/toy.py` runs `sampler.py` on Boltz's own
  `AtomDiffusion` (copied to `addon/sandbox/`) with a Gaussian-ubiquitin denoiser, parity equal; truth from 10^6
  plain samples, 592 steered runs on 18 cores in 11 min. All four settings pass (mean Z-hat / Z within 2 standard
  errors of 1 at 1,024 and 4,096 particles), so the weights and `log_z` are right and I61's shortfall points to the
  estimator's variance with Boltz-2's network. A toy run changes with the batch size (rounding, amplified by
  resampling), like the GPU's. Files: remote_jobs.md, plan.md, diffusion_dynamics.md, progress.md;
  scratchpad/boltz_hdx/addon (`toy.py`, `toy_readout.py`, `sandbox/`, `runs/toy/`).
* **10-29 13:05-13:25, Mac Studio, the rescoring and resampling at low noise only** (user: "run both for me"). The
  eight unscored Upside 3 de novo arms and their Upside 2 controls submitted by `submit_scores.sh` (49643806-49643821).
  `sampler.py` gained `resampling_t_max` (`e1_rare.py`, `toy.py`: `--resample-below`), inert by default (three toy
  runs bitwise equal); rule fixed first in diffusion_dynamics.md Phase 2 E1, limit sigma_data 16 A. Toy at the limit:
  all four settings pass (576 runs). Then 61805263 on a V100: I61 kappa 2, 256 particles, seeds 1 and 2; its readout
  `addon/lownoise_read.py` (checked against `e1_readout.py` on a tuning run). A slip: an unsplit ssh variable in zsh
  skipped the beagle3 backup while rsync went through; the backups were put in place from the md5-equal local
  copies. Files: remote_jobs.md, diffusion_dynamics.md, progress.md; scratchpad/boltz_hdx/addon (`sampler.py`,
  `e1_rare.py`, `toy.py`, `lownoise_read.py`, `e1_lownoise_ubqps.txt`, `runs/toy_low/`).
* **10-29 15:20-, MacBook Pro, Upside 3 on the GPU (plan.md Phase 35; user: "start implementing GPU for upside3").**
  Plan approved (`~/.claude/plans/swift-roaming-floyd.md`). Measured the GPUs: midway2 `gpu2` K80 (driver 460, CUDA
  11.2), midway3 `gpu` V100/RTX6000/A100 (V100 driver 535); one CUDA 11.5 binary from midway2 runs on both (G0 probe:
  Threefry bits equal, fixed-point sums bitwise and repeatable). Development tree `/beagle3/trsosnic/yinhan/upside3_gpu`.
  Built G1 (`device.h`, `device_engine`, flattening, fixed-point accumulation, `--gpu` loop, C API, `tests/device_v1.py`)
  and the device nodes of every ff2.1 node but rotamer (G2, G3). CPU bitwise unchanged through each step
  (`tests/regress.py`, 0 of 4260 arrays) after two corrections: pair potentials moved to headers keep internal
  linkage, and `hbond_energy`'s scores stay inline on the CPU (findings 11g). V1 passes on the host backend for
  bonded, backbone and no-rotamer configs of 1ubq, chig, 1dfn; on the K80 the bonded and backbone configs pass but
  for one discontinuity of the CPU model (`hbond_energy`'s exact-zero test). Found two CPU defects (findings 11g):
  that test, and `hb_environment_coverage`'s unwritten derivative lanes (fixed). Files: src (`device*.{h,cu}`,
  `edge_params.h`, `hbond_interaction.h`, `environment_interaction.h`, `eig4.h`, `backbone_steric.h`, the ported
  nodes' `make_device`, `main.cpp`, CMake), `py/upside_engine.py`, `tests/device_v1.py`; scratchpad/gpu_work.
* **10-29 16:15-17:05, Mac Studio, bottom-up side-chain pair tables investigated and recorded** (user: "for SC-SC
  interactions in upside3, I plan to do it in a different way ... Investigate and record this in upside3"; then
  "you need to make corresponding change in hybrid spline table as well"). `upside3.md` 7.4 (new; Steps renumbered
  7.5): the model today, routes A and B, the F0-F5 diagnostic, a geometric screen of the one-bead side chain from the
  hybrid's conformer set (ARG, LYS and TRP first: distal groups move 1.6-2.2 A RMS within one state; planar groups turn
  up to 71 degrees), the literature (two sub-agents; UNRES's Amber/TIP3P PMFs failed with one site where the screen
  points), the hybrid's corresponding change, open questions and BU1-BU6. Found on the way: `rotamer_center_fixed` is
  a trained output of Jumper's fit (its `program_text`), 2.7 A on average from master's `get_beads` points for ARG;
  up.md 2.1 corrected. Files: upside3.md (7.4, 7.5, status, 1.3, 1.4, 8.3.4, 8.5, 11 item 16, 12), up.md 2.1, plan.md
  Phase 36, progress.md; scratchpad/sc_bottomup_20261029 (`sc_geometry_screen.py`, `bead_def_check.py`, outputs,
  `lit_bottomup.md`, `lit_data_mapping.md`).
* **10-29 17:10-18:25, Mac Studio, the user's three rulings on 7.4** (user: "Try to avoid adding a second bead ...";
  "The correction could be simply a spline table of same dimention"; "Find the best ff to use for me"). The F ladder
  rewritten on one bead per state: forms are sums of products of 1-D tables (a full six-coordinate grid would be 136 GB
  over 3,741 state pairs); a bead frame (F4) completes the pair; more states (F5) carry the per-neighbour conformer
  choice, a gap computable from finer states' tables. The correction is a zero-start table of the CG table's
  dimension, linear in its parameters (a product of zero tables has no gradient). Force field: two sub-agents (cut off
  once by the usage limit, resumed), DES-Amber with TIP4P-D chosen, ff19SB/OPC as the pilot's second force field, with
  a salt-bridge check rule fixed before the pilot. findings 10.43 (the lesson), memory `sc-one-bead-cg-tables`. Files:
  upside3.md (7.4, 11), plan.md Phase 36, findings.md 10.43, progress.md; scratchpad/sc_bottomup_20261029
  (`ff_association.md`, `ff_proteins.md`).
* **10-29 21:00-10-30 00:10, MacBook Pro, GPU validation, then the hybrid on the device** (user: "fix the issue.
  after the validation is completed, proceed to implement hybrid GPU support"). A238 done; with it the K80's V1
  passes every node of 9 ff2.1 fixtures but chig_rama's total potential, whose one rounding-scale move came out 0 by
  coincidence (V1 now takes the RMS of 4). V2 (14-replica 1ubq REMD, K80 and V100, with and without rotamer) agrees
  within 2.9 SE; the GPUs are 1.2-3.7x slower per system step than 4 CPU cores. The hybrid: `device_martini.cu`
  (`martini_hybrid_position`, `martini_potential` on the CPU's filtered pair table, `martini_sc_table_1body`,
  `dist_spring`, `angle_spring`), `DeviceDynamics` with masses, fixed and z-fixed atoms, the g-JF step and the PO4
  hold, `--gpu` for `/input/brownian`. CPU bitwise (0 of 4260) after bisecting a sharing that was not: the side-chain
  table and springs keep their CPU code (findings 11g). V1 passes on the host backend for 1rkl 7.0 and 7.1; a paired
  g-JF run differs by 4 ulps after one round. Build g5a installed; V1 and hybrid V2 submitted on the K80 and V100
  (remote_jobs.md line 5). Files: src (`device_martini.cu`, `martini_interaction.h`, `device_engine.{h,cu}`,
  `device_nodes.h`, `device_rotamer.cu`, `spline.h`, `box.h`, `martini_potential.cpp`, `martini_hybrid.cpp`,
  `main.cpp`, CMake), `tests/device_v1.py`, `tests/device_dyn.py`; plan.md, findings.md, upside3.md, remote_jobs.md;
  scratchpad/gpu_work (`regress_hybrid.sh`, `variant.py`, `hyb_dyn/pair.sh`).
* **10-29 22:40-23:45, Mac Studio, the user's two points on 7.4** (the 136 GB and sites; "in simulation they are usually
  inside protein").
  * **Size.** The 136 GB is the dense six-coordinate grid in float32, file and memory; every proposed form is under
    11 MB.
  * **Speed.** ff2.1 ubiquitin and CA II, profiled with `sample`:
    * side-chain bead pairs take 25% of the step, BP 32-34%, the two coverage nodes 27-30%;
    * an isotropic site-pair kernel costs 3.0 ns against today's 9.8-10.5 ns;
    * martini22 sites give 1.6x and 2.8x the bead pairs, so an estimated +1% and +21% step (wrong; measured 10-30).
  * **Shape.** The wall screen: with a centroid bead, ARG, TYR, TRP, PHE and LYS still move their contact wall
    2.1-2.8 A over direction, which products cannot follow. Sites for those five are recommended to the user as a
    question.
  * **Burial.** 30 native PDBs: 58% of nonlocal contacts both buried, 15% both exposed; nonpolar pairs 80% buried,
    charged pairs 54% exposed. Two limits mixed by burial are proposed. A third force-field review (`ff_buried.md`)
    found the Amber line shares its dry apolar pairs, so DES-Amber stays, with buried checks added to the pilot.
  * findings 10.44. Files: upside3.md 7.4 and 11-12, findings.md, progress.md, plan.md Phase 36;
    scratchpad/sc_bottomup_20261029 (`sc_wall_screen.py`, `contact_burial.py`, `speed/`, `ff_buried.md`).
* **10-30 00:00-07:45, Mac Studio, the user's rule on sites and implicit water** ("if more than 1 bead can make it
  more accurate and faster, then yes"; "check if implicit water ff exist and can work here").
  * **Sites in the engine.** `build_sites.py` libraries with more than one bead per state for ARG, LYS, TRP, TYR and
    PHE. Co-located copies reproduce ff2.1 to 3e-5 in energy and force; martini22 sites pass `node_fd.py`. Step +21%
    (ubiquitin) and +39% (CA II), +14% and +33% with isotropic site tables. The 10-29 estimate was wrong (findings
    10.45).
  * **One-bead forms** timed by probe engines (a copy of `src/` with the kernel replaced, never the repo): F2, F3 and
    F4 on the pairs that involve the five +2-4%, +6%, and at least +17-21%.
  * **Shape** (steric proxy, 10 pairs, held-out half): F2 and F3 gain at most 0.013 over F0. A frame on the bead
    (F4) cuts F3's residual 3.2-5 times and leaves 1.2-4.8 times less than sites. Unregularised F4 put false wells at
    overlapped poses on two pairs; a cross-validated smoothness penalty fixed one, bounded angular tables (F4b) both.
    One bead stays by the user's rule.
  * **Implicit water** (sub-agent review, `ff_implicit.md`): none can be the reference; W_water stays explicit; PB or
    an accurate GB is a candidate for U_dry's electrostatics after a buried-charge check.
  * Files: upside3.md (7.4, 8.3.4, 11, 12), plan.md Phase 36, findings.md 10.45, progress.md, memory
    `sc-one-bead-cg-tables`; scratchpad/sc_bottomup_20261029 (`sites/`, `ff_implicit.md`).
* **10-30 00:10-01:55, MacBook Pro, the hybrid's V2 by seeds.** V1 of build g5a passes every node of the 11 fixtures,
  the hybrid's included, on the K80 and the V100. The hybrid's one-seed V2 reached 4.7 (K80) and 5.3 (V100) block
  standard errors, and its CPU control (seeds 7 and 8) 7.6, so the block error understates the hybrid's (lipid D's
  overlapping windows; slow bath relaxation). `tests/device_dyn.py` rewritten as `run` (one arm at one seed, means to
  an npz, trajectories in a temporary directory) and `compare` (Welch's t over seeds, Bonferroni); job scripts
  `v2run.sbatch` (arrays) and `v2cmp.sbatch` replace `v2.sbatch` and `v2c.sbatch`. Build g5b installed. V2 submitted
  as 8 seeds each on broadwl, the K80 and the V100 (remote_jobs.md line 5). Files: `tests/device_dyn.py`; plan.md,
  findings.md, upside3.md, remote_jobs.md, progress.md; scratchpad/gpu_work/g0.
* **10-30 01:55-07:40, MacBook Pro, G5 speed and the hybrid's GPU lipid excess.** Speed (plan.md item 9):
  `device.h`'s `lanes_for`, `AccumLocal`, `slots_for` and a warp-shuffle `block_scan`, used by the pair-list
  rebuild, `martini_potential`, `martini_sc_table_1body`, BP's incidence lists and the capacity-sized kernels; BP's
  passes fused (builds g6a, g6b, g6c). Each is bitwise equal to g5b on the K80; per system step 1ubq without rotamer
  360 to 245 us, full 1ubq 867 to 526 us, the hybrid 14.5 to 6.5 ms; batch scaling on the K80 and V100 in findings
  11g. V2 of the hybrid by seeds now fails: 40 GPU seeds run the lipids 0.35% hot (4.5 SE) against 24 CPU and
  host-backend seeds, with or without FMA. Ruled out by measurement: races, g-JF coefficients, wrong forces during a
  GPU run (`audit.py`), damping and noise differences (`onestep.sbatch`). The no-exchange arms are running. Files:
  `src/device.h`, `device_interaction_graph.h`, `device_martini.cu`, `device_rotamer.cu`, `device_affine.cu`,
  `device_environment.cu`, `device_hbond.cu`; plan.md, findings.md, upside3.md, remote_jobs.md; scratchpad/gpu_work
  (`build_hostbackend.sh`, `bitwise_gpu.sh`, `g0/` scripts, `configs/1rkl70_dyn10.up`).
* **10-30 07:40-10:45, MacBook Pro, the hybrid's GPU lipid excess found and fixed (A239).** The user's `git pull`:
  origin's one commit fast-forwarded, findings.md, plan.md, progress.md and upside3.md merged three-way with the local
  edits (two appended-entry conflicts kept both sides in date order). Without exchange the K80 stayed +0.42% hot, so
  the cause was in plain dynamics. The g-JF position update's force term, below half an ulp of a large coordinate,
  rounded away where the compiler kept the source order (nvcc; clang for z); gcc on broadwl had reassociated it.
  Measured: the reproducer gives the excess's size, and the K80 with the step summed first matches the CPU (+0.1 SE).
  Both engines now use x + b dt/m (p + dt/2 f + beta/2); regress moves only the hybrid dynamics (12 of 4260 arrays).
  V2 on build g7b running (CPU, K80, V100). Files: `src/martini_brownian.cpp`, `src/device_engine.cu`; findings.md
  (A239, 11g), plan.md, upside3.md, remote_jobs.md; scratchpad/gpu_work (`deadzone_repro.py`, `g0/` diagnostics,
  `dev12_home`, `dev12_record.npz`).
* **10-30 10:45-12:20, MacBook Pro, the hybrid's V2 passes on the fixed build; g7b installed.** g7b against the CPU
  (REMD of 1rkl 7.0, 16 CPU, 16 K80, 8 V100 seeds): Bonferroni 1 on both GPUs, lipid KE within 1.8 SE in every
  replica. The old host backend without exchange sat with the CPU (+1.0 SE; 10 of 16 seeds, 6 hit the 3 h limit on
  shared nodes). g7b installed as the GPU tree's `obj/` (byte-equal to `obj_g7b`); `obj_g6c_nofma`, `obj_g7a` and the
  diagnostic scripts (`audit.py`, `audit.sbatch`, `onestep.sbatch`, `noswap.py`, `nsrun.sbatch`) deleted on the
  cluster and the Mac. Files: findings.md (A239, 11g), plan.md, upside3.md, remote_jobs.md.
* **10-30 12:20-14:30, MacBook Pro, the rest of the validation ladder submitted (user: "run them").** ff3.1 configs of
  chig, 1ubq and 1dfn (the same 31 nodes as ff2.1): V1 passes on the host backend. V1 of 14 configs on the K80,
  V100, RTX6000 and A100; the protein V2 by seeds (1ubq without and with rotamer, and on ff3.1); `tests/device_dyn.py`
  now also compares the backbone angles. V3 at a tenth of the Peng duration (user's choice): a harness from
  `bp_validation` (`v3_run.py`, `v3.sbatch`, `v3_analyse.py`, `v3cmp.sbatch`), smoke-tested at 400 tu, then 32 CPU and
  8 V100 arms; a detached driver runs the analysis when they end. Files: `tests/device_dyn.py`; plan.md (item 10),
  remote_jobs.md; scratchpad/gpu_work (`ff31_configs.py`, `v3_driver.sh`, `g0/` V3 scripts, `configs/*_ff31.up`).
* **10-30 17:06-17:20, MacBook Pro, status pass** (user: "Check remote job status. How is ff3.2 training"; not the
  owner, no cron). Nothing of lines 1-3 failed since 14:30 10-29; 14 newly finished arms read finite, KE/1.5kT
  0.999-1.009. ff3.2 (`ff32_nativefix`) at its target 76 (11:58 10-30): round 4 read (glycine's data pull gone, its
  -ln(P_aL/P_aR) -0.425, margin +0.065), epoch-1 panel `n32_01` read (0.414 against released ff2.1's 0.615, below
  `a32_01`); `_mb2` at 76; `ff311_glypair`'s step 34 rerun clean; released ff2.1's Peng DONE (native 0.552, de novo
  0.427). Files: remote_jobs.md (header, Owed, §1, table), plan.md Phase 33 items 6-7, progress.md; cluster
  `WATCH_STATUS.md` and bundle `records/`, `scratchpad/handoff_20261030/round4_read.sh`.
* **10-30 17:30-17:45, MacBook Pro, ff3.2 extended** (user: "What do you recommend to do with ff3.2?", then "Go ahead
  and continue training ff3.2"). `convergence_gate.py` at step 76 by hand: NOT CONVERGED, `dhb` alone (p 0.0045, cut
  0.0050). Extended to 95 by the same chain (49670425; driver md5 `54383004`, step 76's checkpoint, `.chain_starts` 0
  11 65); `submit_new.sh` (dry-run first) submitted `n32_02` 49670431 and `n32_03` 49670432. Files: plan.md Phase 33
  items 7-8, remote_jobs.md, progress.md; cluster `checks/ff31_fix_20261028/` (`gate_step76_ff32_nativefix.txt`,
  `submitted.log`).
* **10-30 18:50-19:55, MacBook Pro, V1 on every GPU type and the protein V2 by seeds pass** (user: "check status").
  V1 on g7b: every quantity of the 11 fixtures and the three ff3.1 configs passes on the K80, V100, RTX6000 and A100,
  with the same worst error per config on all four. Protein V2 (1ubq REMD, 16 CPU, 16 K80, 8 V100 seeds; ff2.1 without
  and with rotamer, ff3.1): Bonferroni 0.087 to 1; the per-seed means over replicas within 2.0 SE. KE/1.5kT is about
  1.005 on the CPU for the protein too (open). V3 running: CPU arms 5-41% at 18:52, 3 of 8 V100 arms 10%, 5 queued.
  Files: findings.md 11g, plan.md (item 10), upside3.md (status), remote_jobs.md.

## 2026-10-31

* **13:54-14:10, MacBook Pro, status pass** (user: "Check remote job status, especially ff3.2 training ... the basin
  difference between alpha L and R for GLY"; not the owner, no cron). ff3.2 reached 95 at 03:10; the step-95 gate by
  hand reads CONVERGED (dhb p 0.027); round 5 read (glycine -ln(P_aL/P_aR) -0.439, ff3.0 -0.598, ff3.1 -0.825;
  rounds 1-4 reproduced exactly); `n32_02` 0.436 and `n32_03` 0.431 read; `n32_04` submitted by `submit_new.sh`. Every
  ConDiv chain at its target (ff3.1 nativefix 76, both debug arms 76; `ff311_glypair` NODE_FAIL in step 65 recovered
  by its insurance link). Files: remote_jobs.md, plan.md Phase 33 items 6-8; cluster `WATCH_STATUS.md`, bundle
  `records/`, `checks/ff31_fix_20261028/gate_step95_ff32_nativefix.txt`.
* **14:15-16:00, MacBook Pro, what is wrong in ff3.2's folding** (user: "Check what is wrong in ff3.2 that fails to fold
  protein"). `n32_04` folds 0.447 (-0.169 against released ff2.1), so ff3.2 stops by Phase 33's rule. Read-only
  diagnosis on the kept panel and Peng frames: tertiary swelling with secondary structure intact; ff3.2 lowers
  unfolded-state energy relative to native on all 16 Peng proteins, about half through the DSE-driven scalars
  (backbone burial scale, H-bond energies); A193's SARW fix missing from ff3.2's trainer; the 3g0m double free
  reproduced once on the login node, memcheck of a small case clean. One 1-core job (valgrind, 49687255). Files:
  findings.md 1.37, plan.md Phase 33 items 7-9, remote_jobs.md; cluster `checks/ff32_fold_20261031/` (`fail_mode.py`,
  `decomp.py`, `decomp2.py`, `dir_split.py`, `crystal_view.py`, `grad_profile.py`, `peng_profile.py`, `peng_blocks.py`,
  `peng_full.py`, `repro_pd.py`, their `.txt` outputs).
* **16:30-17:05, MacBook Pro, ff3.1 deployed and validating; the fix of the all-term loss started** (user: "1. Deploy
  ff3.1, run validations (Peng, glpG, and other tests we ran) 2. Fix the folding issue and continue working on ff3.2";
  asked which ff3.1, the user chose the fixed-target retrain). Before that, read for the user: glycine alpha_L - alpha_R
  of the released ff3.1 -1.064 and of `ff31_nativefix` by round -0.904 to -0.781 (settled), and the temperature reading
  of ff3.2 (6% hotter for global stability, not term by term; findings 1.37). Plan.md Phase 37: `ff31_nativefix` step 76
  extracted and checked, the merge's 5 Peng chunks and collector cancelled, released to the four trees and the Mac
  repo, the merge's results moved aside (`.bak_pre_nativefix_20261031`), and its tests submitted (2 panels, 2 TM jobs,
  32 Peng arms, collector, glpG prep and REMD on Upside 3). Phase 38: two test arms built from ff3.2's build with one
  trainer line each (DSE off; the SI's SARW), checkpoints checked against ff3.2's, SARW change checked statically,
  both chains submitted with their own panel and unfolded-shift hand-off. The watch probe's ff3.1 paths updated (md5
  `cf834d4b`). Files: plan.md Phases 33 (item 6), 37, 38; findings.md 1.37; remote_jobs.md; progress.md;
  `parameters/ff_3.1/rama.dat` (+ backup); cluster `checks/ff31_nativefix_release_20261031/`,
  `checks/ff32_fix_20261031/`, `~/watch_probe/watch_mdw2.sh`.
* **10-31 08:30-22:55, MacBook Pro, V3 passes.** The 32 CPU arms ended by 08:30, the 8 V100 arms by 22:35, each in
  one chunk; the driver's `v3cmp.sbatch` (49703078): Bonferroni 0.12 to 1 per protein and force field, 0.93 overall,
  T_mid within 3 K. The validation ladder of item 10 is complete; the rest of 10.2's ladder is item 11. Files:
  findings.md 11g, plan.md (items 10-11), upside3.md (status), remote_jobs.md.
* **10-31 23:10-23:50, MacBook Pro, GPU work handed off** (user: "record everything in md file and resume when I says
  'resume GPU jobs'. I will restart on another computer"). Commit check for the user: no scaffolding, the Mac's src
  checksum-equal to g7b's source, no CPU source changed since the regress baseline; `parameters/ff_3.1/rama.dat` is
  the other session's ff3.1 release. Readiness told the user (plan.md Phase 35 "Resume here"). Bundle
  `/project/trsosnic/yinhan/checks/handoff_gpu_20261031/`: working-tree patch on `cce51ce0` (checked: applied to HEAD
  it reproduces the tree), untracked files, records, `scratchpad/gpu_work` tools, the regression fixtures. Files:
  CLAUDE.md (the "resume GPU jobs" trigger), plan.md Phase 35, remote_jobs.md, upside3.md, progress.md.
* **10-31 23:25-23:55, MacBook Pro, status pass.** ff3.1's panels read (`paired_npz.py`): 44 domains fold 0.631
  against released ff2.1's 0.603, mdCATH 0.515 against 0.488, level with the merge; its glpG preparation keeps TM4
  (stage 7.0 DSSP mean 0.95, the merge 0.36). Owed glpG blocks read: ff2.1 block 4 TM4 0.05, ff3.0 block 2 0.82.
  cspA de novo requeued after NODE_FAIL on midway2-0189. Phase 38 arms at steps 9 and 6 of 19. Files: findings.md
  1.34, 3.15; plan.md Phase 37; remote_jobs.md; cluster `glpg_prep_tm4.txt`, `glpg_u3_block4.txt`,
  `glpg_ff30_block2.txt`.
* **11-01 00:05-00:35, Mac Studio, "resume GPU jobs".** Set up from the bundle `checks/handoff_gpu_20261031/`: the
  committed code (HEAD `7f31f72`) equals the GPU tree's `src py tests cmake`; `scratchpad/gpu_work` and
  `upside3_audit` restored; `dev12_home` built and regressed, 0 of 4260 arrays against `dev12_record.npz`. User chose
  (a), then "run both" (a and the hybrid's CPU validation). (a): the GPU tree's `ff_3.1/rama.dat` replaced by the
  release (old kept as `.bak_pre_nativefix_20261031`), the three ff3.1 configs rewritten (only the Ramachandran
  tables differ), V1 passes on the host backend and the V100 (61986517). The hybrid: glpG's two running blocks on
  Upside 3 healthy (no error line), read at their ends 11-02. Files: plan.md Phase 35, remote_jobs.md line 5,
  findings.md 11g.
* **11-01 01:30 and 08:50 CST, Mac Studio, status passes.** ff3.1's mdCATH TM read (mean TM 0.731 against released
  ff2.1's 0.727, +0.004 [-0.003, +0.012]). ff3.1's 32 Peng chunks resubmitted; Phase 38's `ff32_dse0` (DSE target off)
  read: `z32_00` folds 0.575, +0.127 over `n32_00`, -0.040 [-0.076, -0.007] against released ff2.1, unfolded shift
  +3.8 E_up (was -11.9); `ff32_sarwsi` in step 17. Files: findings.md 1.34, 1.37; plan.md Phases 37, 38;
  remote_jobs.md; cluster `WATCH_STATUS.md`.
* **11-01 09:10-11:40, Mac Studio, hybrid glpG HDX lone points** (user: the ff3.0 dG figures on the 11-02 deck's
  slides 17-18 show finite points inside the off-scale TM bands, the implicit `79*.png` do not). Diagnosed (findings
  5.8): each lone point rests on 1-40 exchange frames, seen at one plotted T because MBAR on the coupled energy pools
  ~5 rungs per T (implicit: ~20); part of them are flickers of the 7 A tail-contact test at the bilayer midplane.
  Six shielding rules scored on per-frame geometry of all 112 replicas (49726840), wrongly reported as clearing
  none (a proxy metric; lesson 10.46). User chose "A61 fix + censor": the four sets regenerated with the repo's
  accessibility code (49727925), dG by the repo pipeline on the Mac Studio (reproduces the cluster's censored flags
  and bounds exactly), `plot_ref_style.py` draws censored amides off scale.
* **11-01 14:07-14:30, MacBook Pro, status pass** (midway2 socket reopened, one Duo push). Phase 38 step 1 read: DSE
  off (`z32_00`) folds 0.575 against `n32_00`'s 0.449 and released ff2.1's 0.615, unfolded shift +3.80 E_up; the SI's
  SARW (`s32_00`) 0.446, -9.93. Step 2 put to the user. ff3.1's mdCATH TM level with released ff2.1 (+0.004). Upside 3
  ubiquitin de novo COMPLETE (final frames clean); the last two Upside 3 de novo arms and their controls queued for
  scoring (49736365-49736368). Second NODE_FAIL requeue (alpha3d native, midway2-0156). Files: findings.md 1.37,
  plan.md Phase 38, remote_jobs.md.
* **11-01 14:20-15:40, Mac Studio, HDX lone points continued** (user: TM4 and TM5 still show finite points). The membrane term got a core
  test: inside the ester planes the contact radius is the core amide N-tail g(r)'s first minimum, 8.12 A, measured
  like the 7.00 A; `first_minimum` now searches the first well. Checked against the geometry arrays (1.00000),
  regenerated on the cluster (49737262), dG on the Mac Studio. TM5 clears in all four variants, TM4 in 79HIS; three
  single-event TM4 points remain (findings 5.8). Slides 17-18 rebuilt from `figs/hdx_ff30_core_20261101/`,
  validated against the previous build and the 09-14 original, only those slides, notes and images changed. Files:
  `py/martini_hdx_membrane_accessibility.py`, `example/16.MARTINI/readme.md`, `example/00.AnalysisScripts/README.md`,
  `example/00.AnalysisScripts/helpers/plot_ref_style.py`; findings.md 5.2, 5.3, 5.3b, 5.8, 10.46; plan.md Known
  Errors; remote_jobs.md section 1 line 2.
* **11-01 14:30-15:00, MacBook Pro, ff3.2 continued** (user: "Proceed to continue work on ff3.2. Double check it is
  training all parameters, not holding anything fixed"). Parameter check (`param_moves.py`, `ckpt_rates.py`,
  `params_check.txt`): every group with an engine gradient moved in epoch 0; held only what the engine gives no
  derivative for (backbone burial centre, sharpness and H-bond weight; environment coverage geometry; H-bond entries
  4-11; bead positions) and the basin offsets without data. `rg_check.py`: without the DSE unfolded chains are 2.2%
  smaller than ff2.1's, with it 3.9% larger. New `after_training.sbatch` for `ff32_dse0` (panel, shift, Rg check, next
  epoch, gate from 76; gate tested on step 19); training on from step 19 (49737476). Files: plan.md Phase 38,
  findings.md 1.37, remote_jobs.md; cluster `checks/ff32_fix_20261031/`, `ff32_dse0/after_training.sbatch`.
* **11-02 00:00-00:20, Mac Studio, status pass.** The MacBook Pro's 14:44 records (Phase 38 read, the user's ff3.2
  decision, its lesson) merged three-way into this Mac's (base `7f31f72`; its lesson renumbered 10.47). midway2 faulty
  since 21:19 11-01 (Slurm "Unexpected message received"; from 23:57 its login node opens no session): ff3.2's link
  49737476 failed in epoch 1 (insurance queued), 28 of ff3.1's 32 Peng chains broke at resubmission; restart script
  `checks/ff31_nativefix_release_20261031/restart_broken.sh` written, not run. Read through midway3: ff3.0's glpG
  block 3 (done, TM4 0.85 at T 0.70) and `compare_u3_u2.py` over all 16 de novo pairs (-0.030 [-0.049, -0.014],
  p 0.004). Files: findings.md 1.37, 3.15, 10.47; plan.md Phases 24, 38; upside3.md; remote_jobs.md; progress.md.
* **11-02 03:10-07:35, Mac Studio, watching for midway2.** The midway2 socket master died 03:11; both login nodes take
  TCP on port 22 but send no SSH banner (tested from midway3), so no reconnect was tried and a banner watch runs
  through midway3. ff3.1's glpG block 1 on Upside 3 ended 07:19 and was read on the midway3 login: TM4 DSSP 0.74 at
  T 0.70 (ff3.0 0.85, the merge 0.42), falling within the block; KE, C-N and TM1 as the other runs. Files:
  findings.md 3.15; plan.md Phase 37; remote_jobs.md; progress.md.
* **11-02 12:15-12:55, MacBook Pro, status pass; ff3.2's end of training** (user: "check remote job status, ff3.2
  training status", then "resume broke jobs ... we do not need convergence, as stated by Peng and Jumper. figure out a
  way to determine end of training", then "check ff3.2 basin depth difference between GLY alpha L and R, and whether
  it keeps protein folded"). midway2's Slurm came back about 09:43 and wrote the jobs lost in the fault as TIMEOUT;
  ff3.2's insurance link resumed epoch 1 at minibatch 13. The broken Peng chains and the seed scoring were restarted
  by another session on this Mac (12:18); checked here: 30 of 30 resumed from their checkpoints with the law. Peng's
  SI and Jumper's thesis and paper re-read (findings 1.23): neither tests convergence or picks a checkpoint on held-out
  proteins. The convergence gate is replaced by FF2's 76 steps: `ff32_dse0/after_training.sbatch` now stops at 76
  (md5 `11bd82e9`; run under stubs at steps 38 and 76; `.bak_pre_budget_20261102`), and step 76 is the result unless
  an earlier epoch folds better at p < 0.05. Epoch 0 read: glycine's pooled alpha_L - alpha_R -1.064 to -0.877 (ff3.0
  -0.836); `z32_00` keeps helical and left-handed glycines level with released ff2.1. Files: plan.md (end of training,
  Phase 38 items 6-7), findings.md 1.23 and 10.48, remote_jobs.md, progress.md.
* **11-02 12:10-13:15, MacBook Pro, a second session: midway2's return and Boltz E1's particle count** (user: "check
  remote job status, GPU jobs status", then "Restart 30 Peng chains" and "Score 8 seed replicates", then "run the GPU
  jobs. We need to continue working on trying to make Boltz work for HDX" and "Particle scaling, full"). Opened the
  midway2 socket (one Duo push); restarted ff3.1's 30 broken Peng chains (49743748-49743779) and submitted the seed
  scoring (49743788-49743795; the four alpha3d jobs done 12:34). Boltz: the 256-particle traces show the weights
  collapsing onto a few particles, so the low-noise setting goes out at 1,024, 4,096 and 16,384 particles on midway3
  V100s (62073874-62073880, readout 62073881), its rule fixed first in diffusion_dynamics.md Phase 2 E1; no code
  change, `scale_read.py` checked against the earlier runs. Files: diffusion_dynamics.md, remote_jobs.md, plan.md
  Phase 37, progress.md; `$B/addon/scale/`, `scale_read.py`; WATCH_STATUS.md.
* **11-02 13:15-13:45, MacBook Pro, the second session: Upside 3 GPU validation, the soluble remainder** (user: "is the
  upside3 GPU validataion completed?", then "continue the validation" and "Soluble remainder"). Readings fixed first in
  plan.md Phase 35 item 11. `tests/device_v1.py --frames`, `tests/device_dyn.py nve` (the logged kinetic energy is per
  atom: total = potential + n_atom x kinetic), `g0/v1x_prep.py`, `g0/v3_dist.py` on a refactored `v3_analyse.py`,
  `g0/v2seq.sbatch`; host-backend tests pass on the Mac. Submitted V1 at scale (16 Peng proteins, BP 1e-3 and 1e-6,
  10 frames each), V2 drift (1ubq full and ff3.1, 16 CPU and 8 V100 seeds) and V3 distributions (from item 10's data).
  Files: plan.md Phase 35, remote_jobs.md line 5, progress.md; tests/device_v1.py, tests/device_dyn.py; GPU tree
  `checks/g0/`.
* **11-02 15:10-17:00, MacBook Pro, the second session: readings, A240, the handoff** (user: "estimate completion
  time for upside3 GPU validation and implementation", then "fix the bug you find, then dump everything into md
  files. I will start on another computer"). V3 distributions and V2 energy drift pass (Bonferroni 1, and 0.65 and
  1); both engines drift about 90 E_up per 1000 tu without the thermostat (Known Errors). V1 at scale failed at
  BBA's and NTL9's hottest native frame: weak H-bond occupancies kept only whole steps of 2^-24, from `1-hb` and
  `1-expf(-s)`. Fixed as A240 with `log1pf` and `expm1f` on both engines (`expm1f` alone changed nothing, which
  located the per-pair subtraction); reproducer, node_fd, regress (835 arrays, all reachable from protein_hbond) and
  host-backend V1 pass; GPU build g8a made, its V100 V1 queued (62096441, 62096444). Corrected plan.md: glpG's TM4
  REMD needs no rigid group or barostat. Read the Peng seed replicates (top7's engine gap beyond its seed spread).
  Wrote the handoff bundle `checks/handoff_gpu_20261102/`. Files: src/hbond.cpp, src/hbond_interaction.h,
  src/Float4.h, src/device_hbond.cu; plan.md Phases 24 and 35 and Known Errors, findings.md 3.14 (A240), 3.15, 11g,
  remote_jobs.md, progress.md; `scratchpad/gpu_work/` (`a240_repro.py`, `a240_regress_check.py`,
  `g0/v1home.sbatch`, `configs/v1x/`, `dev13_record.npz`, `dev14_record.npz`).
* **11-02 12:55-13:55, MacBook Pro, ff3.2's remaining loss** (user: "find the source of the remaining loss and try to
  resolve it"). Diagnosed on `z32_00` (DSE off, -0.040 folded against released ff2.1) with the 10-31 scripts and a new
  `decomp3.py`: the loss is tertiary (long-range contacts -0.03, Rg +1%), through a directed change of the side-chain
  hydrophobe coverage (-0.71 E_up of the loose-minus-native gap) and burial (-0.34) that favours the crystal over the
  thermal native (3.79 E_up) and lowers compact non-native states that no replica presents to the native target.
  The user chose the fix: the restrained replica held in the native basin (a one-sided wall at 1.2 d0 on native CA
  contacts, the panel's definition) instead of the 1.0 A elastic network. `ConDiv_basin.py` written and smoke-tested
  on 3g0m (wall exact; restrained RMSD 0.92 to 1.33 A at Q 0.96); arm `ff32_basin` built and submitted, one epoch.
  Six one-file reset panels (14:20): the side-chain file (+0.054) and the burial (+0.041) each bring `z32_00` level
  with released ff2.1, the backbone burial's reset costs -0.096, as the decomposition predicted. On the user's
  "double check we are not forcing ff3.2 to produce the result we want", every ff3.2 edit was audited and the arm's
  decision rule fixed before its results (plan.md Phase 38 item 10). Files: plan.md Phase 38 items 8-10,
  findings.md 1.37, remote_jobs.md, progress.md; cluster `checks/ff32_rest_20261102/`.
* **11-02 14:25-16:50, MacBook Pro, ff3.2's completion estimate and the handoff** (user: "estimate completion time
  for ff3.2", then "dump everything into md files. I will start on another computer"). Step times read from both
  logs: `ff32_dse0`'s link on midway2-[0316-0327] takes 4,000-4,600 s a step since the fault, while `ff32_basin` on
  other nodes takes 1,700-2,000 s, so the slowdown follows that link's nodes (cause not found). Estimate: step 76
  between 15:00 11-03 and 15:00 11-04, the full validation about 3 days after; with the basin restraint adopted,
  1-2 days later. Handoff: records updated, the glycine depth scripts made path-independent (they reproduce round 1's
  -0.877 / -0.640) and copied with the records, the working-tree patch and this Mac's memory into the cluster's
  handoff bundle. Files: remote_jobs.md (header, §1 line 2, "Resume here" step 4), plan.md Phase 38 "Resume here",
  progress.md; bundle `checks/handoff_20261026/` (README addendum, `records/`, `scratchpad/handoff_20261102/`).
* **11-02 17:15-17:30, Mac Studio, A240's GPU check and install** (user: "continue upside3 GPU jobs monitor"). Opened
  both sockets. V1 on build g8a (V100) read by plan.md Phase 35 item 11's rule: 62096441 passes every quantity of the
  32 Peng configs at both BP tolerances (BBA and NTL9, which failed `hbond_energy` on g7b, now worst 8.0e-6 and
  1.4e-5), 62096444 every quantity of the 14 standard configs (1dfn's `hbond_energy` worst 1.5e-4 gone). No job of the
  GPU tree on either cluster, so `install_gpu.sh g8a` (17:22; `obj` = `obj_g8a`, g7b kept as `obj_g7b`). Item 11's
  soluble part is complete. Refreshed this Mac's `scratchpad/gpu_work` from the 11-02 bundle (two replaced scripts
  kept as `.bak_pre_mbp_20261102`). Files: plan.md Phase 35 ("Resume here", item 11), findings.md 3.14 (A240) and
  11g, remote_jobs.md (header, "Owed by hand", job table, line 5), progress.md.
* **11-02 17:18-20:40, Mac Studio, ff3.2's epoch 1 and the throttled nodes** (user: "continue ff3.2 monitor"). Opened the midway2 socket.
  Epoch 1 of `ff32_dse0` ended 16:47 (24 of 24 proteins every step, KE/1.5kT 0.994-1.014, H-bond margin +0.126) and
  its after-training job queued `z32_01`, its shift and epoch 2 (49748334). Read by Phase 38 item 6: `z32_01` folds
  0.578, level with `z32_00`, -0.038 [-0.072, -0.005] under released ff2.1; unfolded shift +5.61 E_up; round 2's
  glycine alpha_L - alpha_R -0.791; the free replicas' Rg at T 1.1 0.957 of epoch 0's, so the unfolded chain keeps
  compacting, against item 6's condition (reported). The 70-min steps are the nodes' clock: on midway2-[0316-0327],
  where epoch 2 ran again, every busy core sits at 1,199 MHz against 2,500-2,900 on normal nodes. At the user's word
  (18:10) the 12 nodes went into the node law (README row, `--test-only` checked, `update_pending.sh`) and epoch 2's
  link 49748334 was cancelled 18:12 after step 39's checkpoint; its insurance 49748335 resumed there 19:11 on 22
  normal nodes (cores 2.0-2.8 GHz, checked), ran step 40 and lost a node 19:49; 49748605 resumed at step 40 the same
  minute (step 41 1,681 s). No ff3.2 job writes to `/beagle3`, whose inode quota is spent (20:35 check).
  `ff32_basin` at step 11 of 19 by 20:12, clean. Split the work with the GPU session on this Mac (it owns the watch, cron
  `1b3ab805`; this session keeps the ff3.2 entries). Files: remote_jobs.md (header, "Resume here" step 4, §0d, §1
  line 2, job table), plan.md Phase 38 ("Resume here", item 6), findings.md 1.37, progress.md; cluster
  `slurm/midway2.args` and `README.md`; `scratchpad/handoff_20261102/` fetched from the bundle (round files in
  `gly/rounds/`, the clock probe `node_clock.sh` added).
* **11-02 17:35-17:56, Mac Studio, the watch and the plan for (c)** (user: "plan that, while still monitoring").
  Took the jobs: watch cron `1b3ab805` (hourly at :23), owner the Mac Studio's session `upside2-md-a0`; the other
  Mac Studio session `upside2-md-fe` reads ff3.2 by hand and owns its entries (agreed by message). Surveyed the
  cluster runners (soluble REMD, the mdCATH panel, glpG's Upside 3 runs) and found on the host backend that one
  `--gpu` process cannot hold different proteins (`max_n_edge`, `rama_map_pot`; findings 11g). Wrote plan.md Phase
  35 item 12. Files: plan.md Phase 35, findings.md 11g, remote_jobs.md (owner lines), progress.md; cluster
  WATCH_STATUS.md.
* **11-02 17:58-18:22, Mac Studio, item 12 steps 1-2** (user approved steps 1-2 with glpG on ff3.1, mdCATH on ff2.1
  and ff3.1, the largest system's capacity). Device batches of different proteins: `DevicePairlist` takes each
  system's `max_n_edge` (slabs of the largest, each system held to its own; `backbone_pairs` and the rotamer's
  residue pairs the same), `rama_map_pot` uploads each system's maps (equal ones once); `device_v1.py --batch`.
  Host backend: one-protein V1 identical to HEAD's, mixed batches pass every quantity, trajectories bitwise, regress
  0 of 4260. The hybrid `run.py` gains `remd_engine` (its stale per-block comment corrected). Build g9a, `home_g9a`,
  glpG's V1 frames; 7 midway3 jobs submitted 18:20 (bitwise, mixed V1, glpG V1, two glpG probes through `run.py`).
  Files: src/device_interaction_graph.h, src/device_nodes.h, src/device_rotamer.cu, src/device_affine.cu,
  src/backbone_steric.cpp, src/rama_map_pot.cpp, tests/device_v1.py, example/16.MARTINI/run.py; plan.md Phase 35
  item 12, remote_jobs.md line 5, progress.md; `scratchpad/gpu_work/g0/` (`v1g_prep.py`, `glpg_probe.sh`).
* **11-02 20:26, Mac Studio, watch pass and item 12 steps 1-2 read.** g9a passes on a V100: bitwise to g8a (1ubq,
  1rkl), mixed V1 (ff2.1, ff3.1), glpG's V1 at both tolerances. The V100 glpG probe ran 48 replicas at 850 tu an hour
  (14 times the CPU block), 4.4 of 16 GiB, then failed writing `block_count`: `/beagle3/trsosnic`'s group inode soft
  quota (1,276,989 of 1,268,000) had its grace expire between 19:30 and 19:46, which also cost Boltz E1's two
  16,384-particle arms their results and failed its readout. `/project/trsosnic` is 140 G over its block soft quota
  with 6 days of grace. Install of g9a held; the deletions are the user's. Files: remote_jobs.md (header, Owed by
  hand, lines 4-5, job table, disk section), plan.md item 12, findings.md 10.39, progress.md; cluster WATCH_STATUS.md.
* **11-02 21:51-21:56, Mac Studio, disk cleanup** (user: "do a cleanup of the disk, check what is no longer used").
  `/beagle3` files from 1,276,989 to 1,263,519, under the 1,268,000 soft quota: the records' safe-to-delete
  discards and the read V100 probe replicas deleted; `chi1_perplexity`, `pylibs` and the GPU builds before g7b packed
  into compared tars in `/project/trsosnic/yinhan/checks/disk_reclaim_20261102/`. File creation works again; every
  job started in the block swept (ff3.1's collector and scorings wrote to `/project`). g9a installed; the A100 probe
  resubmitted (62120973). `/project`: at the user's choice the retired ff3.1 glpG chains' replicas (180 G) deleted
  22:40, 23 G under its block soft quota.
  Then (user: "some ff3.0 trajectories can be cleaned", choosing the failed and superseded Peng sets) the replica
  files of `ff3_benchmark/runs/*_ff_3.0` and `runs_superseded/ff_3.0_20260930-022510` (44 G, 896 files) deleted 22:56.
  11-03 00:24-00:38 (user: "clean up old boltz and bioemu2 files, if no longer in use"): Boltz Phase 0-1 runs'
  data (23.8 G, logs kept) and BioEmu2's samples, env and weights (11 G) deleted; Boltz Phase 2 untouched.
  Files: remote_jobs.md (header, Owed by hand, line 5, job table, disk section, 1b), plan.md item 12, findings.md
  (pylibs pointer), progress.md.
* **11-03 00:37-03:10, Mac Studio (`upside2-md-fe`), ff3.0/ff3.1 trajectories and `ff32_basin`'s end** (user: "We can
  delete some ff3.0 and ff3.1 trajectories to free some space"). Agreed with the watch session, which had no ff3.x
  reclaim in flight; surveyed every ff3.0/ff3.1 set (size, files, what still reads it) and put four to the user, who
  chose the retired ff3.1 sets and the ff3.0 HDX frame copies. A checked manifest (2,300 `*.run.<r>.up` files, no
  hardlink, none written in 24 h) and `reclaim_ff3x_20261103.sh` deleted 486 G 00:50-00:57: `/project` 3,468 to
  3,094 GiB (480 G under its soft quota), `/beagle3` 1,392 files; the HDX results, outputs and the chains' replicas
  verified kept. `ff32_basin` COMPLETED 00:48 (19 of 19 steps clean); its Rg check reads 0.943 at T 1.1. Its
  mechanism job first failed on the half-written panel and was requeued after it (with `tm_panel`); all readouts
  ended by 02:56 and were read 03:05 by plan.md Phase 38 item 10's rule: (a) and (b) both met (folded 0.630 against
  `z32_00`'s 0.575, TM level, crystal and loose gaps about halved, 3g0m's 3-7 A dip gone); reported to the user, who
  decides. Files: remote_jobs.md (header, Owed by hand, line 2, job table, disk section, 1b), plan.md (Phase 38
  Resume here, item 10), findings.md 1.37 "The fix", progress.md; cluster `checks/disk_reclaim_20261102/` (script,
  manifest, `reclaim.log`), `checks/ff32_rest_20261102/mech_w32_00.sbatch`.
* **11-03 07:30-07:50, Mac Studio (`upside2-md-fe`), ff3.2 continues under the native-basin restraint** (user: "Ok.
  Proceed", to the recommendation to adopt `ff32_basin` and let `ff32_dse0` finish as the comparison). plan.md Phase
  38 item 11 written first. `ff32_basin`'s `after_training.sbatch` replaced by `ff32_dse0`'s per-epoch form with tags
  `w32_EE` and a new `mech.sbatch` (crystal and loose gaps per epoch), the one-epoch version backed up; checked by a
  stub run of every branch and `--test-only`, no test job. Training keeps no frames, so the proposed hot-replica
  native-contact check was replaced by the per-epoch gaps. First link 49757705 submitted 07:39, PD on priority
  (broadwl full) behind `ff32_dse0`'s epoch-3 link, which with `z32_02` waits since 05:50; `ff32_dse0`'s step 57
  reached 05:50, Rg at T 1.1 0.923, round 3's glycine alpha_L - alpha_R -0.758. Files: plan.md (Phase 38 Resume here,
  items 7 and 11), remote_jobs.md (header, Resume here, line 2, job table), progress.md; cluster
  `checks/ff32_rest_20261102/` (`after_training_ff32_basin.sbatch`, `mech.sbatch`, `stub_epochs.sh`,
  `stub_epochs.txt`, `submitted.log`), `$P/training/ff32_basin/after_training.sbatch`.
  08:05 (user: "Decide for me"): `ff32_dse0`'s epoch-3 link set to nice 100, so ff3.2's link starts first.
* **11-03 07:35-08:15, Mac Studio (`upside2-md-a0`), `/beagle3` file quota** (user: "Delete some files to clear the
  file limit", "Maybe tar some of them then delete"). `du --inodes` over `/beagle3/trsosnic`: ours are 209,049 of the
  group's 1,254,101 (soft quota 1,268,000); `anaconda/` (a removed account) 507,205, `schwartznw/` 339,897, `yiheng/`
  180,347. Three unused-by-any-job candidates put to the user, who chose Boltz's two: `reclaim_boltz_phase2.sh`
  deleted `boltz_hdx/cache/mols/` (45,228, the extracted copy of `mols.tar`, compared first) and packed `boltz_hdx/.venv`
  (45,422) into `boltz_venv.tar` on `/project`, compared and counted, then deleted it. `df -i` 1,254,102 to 1,163,475;
  both are unpacked before any E1 rerun. Files: remote_jobs.md (Owed by hand, line 4, disk section),
  diffusion_dynamics.md (Status), progress.md; cluster `checks/disk_reclaim_20261102/` (script, manifests, compares,
  `reclaim.log`, `boltz_venv.tar`).
* **11-03 11:50-12:45, Mac Studio, the most-probable-state side chain** (user, the PI's request: "after we have 50%,
  30%, 20% probabilities, we only accept 50% state ... We will need to retrain the ff"; "decide them for me, and go
  ahead and implement and deploy it. Do not make the core engine change git commitable"). plan.md Phase 39, findings
  1.38. Patched `$P`'s `rotamer.cpp` (attribute `most_probable_state`) in the git-ignored
  `scratchpad/bpmax_20261103/`; checked on the Mac (BP bitwise, max-mode energy rebuilt independently, forces and
  parameter derivatives by central differences, dynamics: KE/1.5kT 0.8-2.4% above BP, NVE drift 8-40 times BP's);
  built the separate tree `/project/trsosnic/yinhan/upside2-md-bpmax` on midway2 (BP bitwise `$P/obj` there too),
  the run `ff21_bpmax_dt009` from the control's inputs, and submitted setup and the first link to step 76
  (remote_jobs.md §1 line 1). Files: plan.md, findings.md, remote_jobs.md, progress.md; repo `src/` and `py/` untouched.
* **11-03 12:50-13:05, Mac Studio, the third arm's program state checked** (user: "Double check you train it from
  ff2.1, and the program should be the same state as the other tests"). init_param md5-equal to released `ff_2.1`
  and to both arms'; trainer, Python, venv and inputs the same. The engine was not: both earlier arms ran steps 0-37
  before the angular-gradient fix (10-18). The first link (49766483, started 12:33) and its successor were cancelled
  in step 0, a second build `obj_preangfix` made and shown bitwise to the earlier arms' steps 0-37 binary in BP mode,
  and the run restarted as link 49766810 to step 38 with an automatic switch to the fixed build at 38 (stub-tested).
  Files: plan.md, remote_jobs.md, progress.md; cluster run README, `setup/{at_target,switch_engine}.sh`.
* **11-03 15:00-16:35, Mac Studio, HDX m-values for Boltz-2 and BioEmu2** (user, relaying the PI: "in Peng 2022
  paper, Figure 5, it is not only delta G we care about. We also care about the slopes"). Peng's Figs 5-7 and 10 and
  Table S14 score each amide's dG and its m-value (slope against urea); our readouts scored dG only. The 0 M m-value
  is s times the protected NHs lost on opening, so its correlations need no calibration. Experiment's m (eq. S18 on
  Tables S5-S7) reproduces Table S14's FF2 m-value correlations to within 0.03. Boltz-2 opens every measured amide
  locally (0.9-6.1 NHs); HEEH m Spearman -0.18 (FF2 0.475); EHEE and Ub psWT unresolved at 10^5. BioEmu2 HEEH 0.465.
  The high-urea half needs s from each model's own global unfolding, which Boltz-2 does not sample. 16:00-16:20
  (user: "yes. proceed"): BioEmu2's HEEH high-urea half, s from its own unfolded samples. At 0 M the whole Table S14
  row matches FF2; at 2.13 M it is not determined, since 4-8 unfolded samples carry it and one holds 6% of the weight.
  Then the MacBook Pro uploaded BioEmu2's EHEE and Ub psWT 10^5 readout files to `$E/readout/` (16:25); pulled
  md5-equal, `.resid` identical, the 10-27 dG readout reproduced. EHEE opens only globally (14-23 NHs; m Spearman
  0.18, FF2 0.40), and its global stability is 3.5-3.7 kcal/mol against 8.57, so its 4 M row cannot be compared
  with FF2's. Ub psWT opens 6 of 23, all locally, too few to score. Files: diffusion_dynamics.md (Status, Phase 1,
  section 7), remote_jobs.md (BioEmu2 paragraph, 11-03 cleanup line), findings.md 10.50, progress.md; scratchpad
  `boltz_hdx/flags/urea/` (`parse_si.py`, `check_s14.py`, `mvalue.py`, `high_urea.py`, their logs) and
  `flags/reference/*_urea.csv`.
* **11-03 16:38-16:52, Mac Studio, Boltz E1's particle-count readout** (user: "yes", to reading the finished runs
  after asking whether frozen Boltz can be made to sample HDX-correct frames). The eight finished runs pulled from
  `$B/addon/runs/` and read with `scale_read.py` here (its readout job had failed with the lost 16,384 arms). The 256
  runs reproduce; at 1,024 and 4,096 particles no seed reaches the -11.1 line or calibration (Feynman-Kac open
  fraction 2e-8 to 1.3e-5 against 3.8e-4), and the mean does not rise from 1,024 to 4,096. Formally undecided without
  the 16,384 arms. Files: diffusion_dynamics.md (Status, Phase 2 E1), remote_jobs.md (header, Resume here, line 4
  detail, job table), progress.md; `scratchpad/boltz_hdx/addon/runs/scale_readout.txt`, also at `$B/addon/runs/`.
* **11-03 17:38-17:55, Mac Studio, ff3.2 progress check** (user: "check ff3.2 training progress"). `ff32_basin`
  (ff3.2) link 49757705 R from 11:28, step 31 at 17:33, 30 min a step; `check_step.py` on steps 19-31: 24 of 24
  proteins, KE/1.5kT at most 1.015, H-bond margin +0.098 to +0.106 (under +0.10 at 24-26 and 31; user notified).
  `ff32_dse0` link 49754482 step 69 at 17:37, margin rising to +0.159. `z32_02` read: 0.599, level with released ff2.1
  (-0.016 [-0.059, +0.026]); shift +6.38 E_up. Files: remote_jobs.md (header, Resume here, line 2, job table),
  plan.md (Phase 38 Resume here step 3, items 6 and 11), progress.md.
* **11-03 17:10-18:15, Mac Studio, E1 closed; Phase 3 and its toy test T0** (user: "proceed with your recommendation.
  but don't give up on the big goal, try another potential solution"). E1 closed in diffusion_dynamics.md and
  remote_jobs.md (its finished jobs' row deleted). Phase 3 written with its rules before the data: replica exchange on
  frozen Boltz-2 whose Metropolis move is Boltz-2's own renoise-denoise, a ladder along Peng's protected-NH count, MBAR.
  `addon/remc.py` (the move as Boltz's sampling loop from any step, bitwise to `AtomDiffusion.sample` from step 0),
  `toy_remc.py`, `t0_read.py`, `remc_boltz.py` (T1's driver, not yet run). T0 on the toy: parity and invariance pass;
  the ladder failed at 400 sweeps (top rungs not equilibrated, rung means still rising) and passed 8 of 8 at 3,000
  sweeps read over the last 1,500 (P(d >= 7.5) 4.15e-4 against 3.70e-4; ln Z within 1 se). T1 (about 3.5 V100-hours,
  node-local venv so no `/beagle3` files) proposed. Files: diffusion_dynamics.md (Status, Decision 6, Phase 2 E1,
  Phase 3), remote_jobs.md, progress.md; `scratchpad/boltz_hdx/addon/` (`runs/t0/`).
* **11-03 18:15-18:40, Mac Studio, Phase 3's T1 submitted** (user: "proceed"). `remc_boltz.py` rewritten for the
  real model (the sampler's conditioning captured from Boltz2.forward, parity and self-repeat checks, CA-RMSD to
  `1UBQ_psWT_ref.pdb` checked against scipy's alignment), `t1.sbatch` (venv and CCD molecules extracted to the node's
  `/tmp`, so no `/beagle3` files; the processed input remade by `boltz predict`), `remc.run` simplified (toy still
  bitwise). Deployed md5-equal, `--test-only` passed, submitted as 62199777 on midway3's V100 queue. `t1_read.py` written
  and tested on synthetic chains; that test found rule 2's reference error zero for a never-opened donor, corrected
  before any data (diffusion_dynamics.md Phase 3 T1). Files: diffusion_dynamics.md, remote_jobs.md (§1 line 4, owed
  by hand, job table), progress.md; `scratchpad/boltz_hdx/addon/` (`remc_boltz.py`, `t1.sbatch`, `t1_read.py`).
* **11-03 18:00-18:50, Mac Studio, GPU: soluble proteins closed; item 12 steps 3-5 started** (user: "close what is
  validated, and proceed to next validations"). Closed in plan.md Phase 35 (Resume here, items 6 and 11) and
  upside3.md 10.2-10.3: soluble proteins on ff2.1 and ff3.1 pass V0-V3 but the mdCATH panel; G2's speed target is
  recorded as not met (one V100 = 14 broadwl cores on proteins). Step 3, glpG V2: `checks/configs/glpg.stage_7.0.up`
  promoted from `glpg_ff31`'s checkpoint (start bitwise the CPU REMD's), 16 CPU and 8 V100 seeds submitted. Step 4:
  `g0/glpg_tm4.sh` (new) submitted one 48-replica V100 block through `run.py`, read later over the CPU block's t
  0-2150. Step 5: `checks/mdc/mdc_run.py` and `mdc.sbatch` (new; panel.py's configs and schedule, 4 length batches,
  either engine), tested on the midway2 login node; the CPU reference's 8 batches and the V100 probe of batch 4
  submitted. Files: plan.md, upside3.md, remote_jobs.md (line 5), progress.md; `scratchpad/gpu_work/g0/glpg_tm4.sh`,
  `scratchpad/gpu_work/mdc/`.
* **11-03 22:05-22:20, Mac Studio, GPU status** (user: "check status"). Step 3, glpG V2 by seeds, passes (Bonferroni
  0.94 over 48, largest 2.6 SE); the mdCATH CPU reference is complete (712 runs to t 9990, finite, KE/1.5kT 1.006);
  step 4's V100 block at 2,110 tu; step 5's probe still on Priority (start estimated 00:00 11-04). Files: plan.md
  (Phase 35 Resume here, item 12 steps 3 and 5), findings.md 11g, upside3.md 10.2, remote_jobs.md (header, line 5,
  job table), progress.md.
* **11-03 23:10-23:55, Mac Studio, GPU: step 4 read, step 5's probe diagnosed and resubmitted** (watcher: 62200634
  COMPLETED, 62201376 FAILED). Step 4: glpG's V100 TM4 block (2,790 tu) read over the CPU block's t 0-2150 by
  `g0/glpg_tm4_read.py` (glpg_u3_block1.py cut to a time window); T 0.70's TM4 DSSP 0.653 against 0.741 (the CPU
  halves' spread), so every rung of both blocks was read (`g0/glpg_tm4_rungs.py`): the difference changes sign along
  the ladder and the ladder mean agrees (0.706 against 0.689, the same fall by quarter); KE and C-N equal. Step 5:
  the V100 probe's log was lost with its temporary directory; on the Mac host backend batch 4's 88 systems run, 22 of
  them for 12 tu match the CPU engine (KE/1.5kT 2.1-4.0 on both before the anneal), and the footprint (32.1 GB
  against the CPU engine's 12.0 GB) shows about 20 GB of device arrays, from every interaction graph's full-pair
  capacity times the systems. `mdc_run.py` now prints upside's log on a failure; the probe resubmitted on an A100
  (62225254). Files: plan.md (Phase 35 Resume here, item 12 steps 4 and 5), findings.md 11g, upside3.md 10.2 and
  10.3, remote_jobs.md (header, line 5, job table, §1), progress.md; `scratchpad/gpu_work/g0/glpg_tm4_read.py`,
  `glpg_tm4_rungs.py`, `scratchpad/gpu_work/mdc/mdc_run.py`.
* **11-04 00:20-00:45, Mac Studio, GPU: probe read, host loop parallelized (g9b), GPU panel submitted.** The A100
  probe 62225254: 19,755 MiB (confirms the V100's memory limit), 68.5 us per system step, host work 19.1% (over
  item 12's 10%). `src/main.cpp`: one `for_each_system` (the CPU branch's trapped OpenMP loop) now runs both the
  CPU branch and the `--gpu` host work. Mac: CPU regress 0 of 4260 against `dev14_record.npz`; `mc_bitwise.sh` (new)
  with 4 threads bitwise equal to g9a on both engines (mixed batches, frames and moves every tu); a mixed batch of
  pre-release and released ff3.1 configs is refused by name, correctly. Batch 3's device memory 10.4 GB (Mac). g9b
  built on midway2, installed 00:36; `g0/bitwise.sbatch` now fails its job on a difference; check 62227671 and the
  panel's 8 jobs 62227672-79 behind it. Files: src/main.cpp, plan.md (Phase 35 Resume here, item 12 step 5),
  findings.md 11g, remote_jobs.md (header, line 5, job table, §1), progress.md; `scratchpad/gpu_work/mc_bitwise.sh`,
  `g0/bitwise.sbatch`, `dev16_home/`, `dev16_record.npz`.
* **11-04 08:30-08:40, Mac Studio, GPU: g9b's check and the mdCATH GPU panel read** (watcher: 62227671-79
  COMPLETED). g9b is bitwise g9a on a V100 (frames, swaps and pivot moves every tu). The panel's 712 runs pass the CPU
  reference's checks (t 9990, finite, KE/1.5kT 1.006, C-N at most 1.85 A); 1.40-1.49 times the 28-core broadwl
  node's speed on V100s (batches 1-3) and 2.04 on the A100 (batch 4). Folded fraction GPU - CPU: ff2.1 -0.002
  [-0.023, +0.019] p 0.91, ff3.1 -0.003 [-0.030, +0.023] p 0.92. TM-score comparison (step 5's test) submitted on
  broadwl. Files: plan.md (Phase 35 Resume here, item 12 step 5), findings.md 11g, upside3.md 10.2,
  remote_jobs.md (header, line 5, job table, §1), progress.md.
* **11-04 10:39-10:50, Mac Studio, GPU: step 5 passes, item 12 closed** (user: "git push all local changes";
  watcher: 49803361-62 COMPLETED). The 7 modified files committed as `be8f652` and pushed to `origin/martini-dev`.
  TM-score of the GPU panel against the CPU reference: mean TM GPU - CPU -0.000 [-0.005, +0.005] p 0.86 (ff2.1),
  +0.000 [-0.005, +0.006] p 0.91 (ff3.1), so the engines agree by step 5's rule; one of 28 subgroup rows at p 0.046
  (ff2.1 mixed stratum, TM > 0.5, GPU higher). Batches of different proteins closed; step 6's report given.
  Files: plan.md (Phase 35 Resume here, items 11 and 12), upside3.md (opening, 10.2, 10.3 G5), findings.md 11g,
  remote_jobs.md (header, checklist, job table, §1 line 5), progress.md.
* **11-04 10:40-11:30, MacBook Pro: status pass, Boltz T1 fixed, GPU step 5 read, closures and disk surveyed**
  (user: "check remote job status", then "fix the failed jobs, check which task branch can be closed, and check what
  can be removed"). Since 09:26 11-03 every job COMPLETED but Boltz T1 62199777, which failed its parity check: the
  module's `num_sampling_steps` is Boltz2's default 5, so `remc` ran a 5-step schedule against the sampler's 200.
  Fixed; the rerun 62262317 passed parity and failed on a strided `atom_idx` index; fixed and the post-check path run
  end to end on the Mac with Boltz stubbed; rerun 62263853. Upside 3 GPU item 12 step 5 passes (mean TM GPU - CPU
  ff2.1 -0.000 [-0.005, +0.005], ff3.1 +0.000 [-0.005, +0.006]). Read beside, not recorded in their lines (owned by
  the Mac Studio): ff3.2 `w32_01` 0.633, `w32_02` 0.674 (+0.059 [+0.015, +0.103] over released ff2.1), `z32_03`
  0.606; bpmax at step 38 on the fixed engine; ff3.1 Peng 23 of 32 COMPLETE. Then (user: record the analysed results
  in the 1102 report folder, then delete them; Mac Studio stopped, no watch): readouts of seven finished datasets
  copied to `GroupMeetings/1102/results_archive_20261104/` (README) and `checks/archive_1102_20261104/`, then 559 GB
  of their trajectories deleted by manifest (`checks/disk_reclaim_20261104/`; `/project` 989 G free, `/beagle3` 1.2 T,
  files 1,163,220). Owner none. Files: `scratchpad/boltz_hdx/addon/` `remc_boltz.py`, `remc.py` (and the cluster's,
  with backups), diffusion_dynamics.md (Phase 3 T1), plan.md (Phase 35), remote_jobs.md (header, Owed, Handoff, the
  watch, line 4, line 5, job table, Disk), progress.md. Then `anaconda/` (user: check its use, a replacement with
  fewer files in a public directory, a tar backup first): no current use found; replacement
  `/project/trsosnic/shared/anaconda_py39` (anaconda's 426 packages at their versions, 166k files, 420 of 429 import
  on both clusters, HDX parity) and backup `shared/archive/beagle3_anaconda_20261104.tar` (637,439 entries); then
  deleted at the user's word, 14:03-14:42: `/beagle3` files 1,163,216 to 656,016 of 1,268,000. Then ff3.0 (user: transitional, delete its trajectories once analysed and on the slides): the owed TM4 over
  blocks 1-5 of its four Upside 2 chains read (frays in the 79ALA backgrounds, holds in the 79HIS ones; findings 1.35),
  readouts archived, 260.3 GB of ff3.0 trajectories deleted (`checks/disk_reclaim_20261104_ff30/`). Then ff3.1 (user: also transitional; all ff3.0/ff3.1 go once ff3.2 is validated): ff3.2 at step 70 of 76;
  129.9 GB of analysed ff3.1 data deleted (`checks/disk_reclaim_20261104_ff31/`), readouts archived. `a32` deleted
  (43.0 GB). ff3.2 epochs 0-2 read (plan.md Phase 38 item 11): glycine alpha_L - alpha_R -0.774 at round 3 (ff3.1
  -0.781); folds 0.674 against released ff2.1's 0.615; the loose gap and hot-chain compaction grow. Diagnosed
  (findings 1.37): both are states the training contrast never contains; compaction is backbone burial (-1.95 of
  -2.23 E_up), the loose gap the directed side-chain packing change; fixes listed in plan.md Phase 38 item 11. The
  user chose no change: epoch 3's read and the epoch choice run detached (`checks/ff32_epoch3_20261104/`, 49816524),
  restart points archived (`checks/ff32_checkpoints_20261104/`).
* **11-04 16:20-16:40, Mac Studio: git pull and a status pass** (user: "git pull and merge with local", then "check
  remote status, find what is completed, what is closed, and if we need to submit any appending jobs"). Pulled
  `986af0b` (fast-forward) and re-applied this computer's uncommitted record edits over it, merging the MacBook Pro's
  Phase 35 and line 5 entries by hand (step 6 closed here; the panel's trajectories deleted there). Status, read only:
  every job since 00:00 11-04 COMPLETED but `tar_anaconda` (stopped in its compare) and T1's 62262317 (rerun);
  running are bpmax link 49806998 (engine switched at step 38, workers on `obj_angfix/upside`, step 49, KE/1.5kT
  0.998-1.037), `ff32_basin` at step 72, ff3.1 Peng's 9 de novo arms (23 of 32 scored) and Boltz T1 62263853 (parity
  passed, sweep 50 of 100). Every running line has its successor queued; nothing to submit. Then (user: "decide for
  me what needs to be run"): Phase 32 judged by its rule, failed (rounds 3-4 averaged, aR +0.036 and aL +0.397
  against ff3.0's -0.107 / +0.517; `ff311_glypair` still pulled, t +2.7), not rerun; H2 moved to ff3.2's ubiquitin
  arms, since ff3.1's native arm's replicas went in the 11-04 reclaim; what follows each ff3.2 verdict fixed in
  Phase 38 item 11; this computer's Boltz `remc.py` / `remc_boltz.py` pulled md5-equal from the cluster (backups
  `.bak_pre_nsteps_20261104`). Files: remote_jobs.md (header, Resume steps 4 and 7, job table, lines 1-2),
  plan.md (Phases 23, 32, 38, 39), upside3.md 8.2, progress.md; cluster `WATCH_STATUS.md` and the bundle's
  `records/`.
* **11-04 16:45-17:05, Mac Studio: the watch taken, with auto-start** (user: "monitor and auto start jobs"). Owner
  the Mac Studio, session `upside2-md-b9`, cron `16048789` hourly at :23 (expires 11-11), whose prompt also submits
  the jobs decided in advance (remote_jobs.md "The watch"). Prepared and dry-run (`plan` mode, `sbatch --test-only`,
  nothing queued, test copies removed): `checks/ff32_validate_20261104/` (ff3.2's validation: stages the chosen
  epoch as `ff_3.2_<epoch>`, mdCATH panel and TM, Peng's 32 arms and collector, glpG on Upside 3) and
  `checks/bpmax_fixedpoint_20261104/` (bpmax's 44-domain panel in its tree and both TM comparisons). Files:
  remote_jobs.md (header, Handoff, The watch, job table), plan.md (Phases 38, 39), progress.md; cluster
  `WATCH_STATUS.md`.
* **11-04 18:25-18:40, Mac Studio: watch pass; Boltz T1 read, fails rule 2.** T1 62263853 COMPLETED 17:37 (parity
  equal); `t1.npz` and its log pulled to `scratchpad/boltz_hdx/addon/runs/boltz_remc_t1_ubqps_v100_62263853/`, read
  by `t1_read.py`: rule 2 fails (two donors at z +4.1 and +3.8, each opened once or twice in the 1,280-sample
  reference, exact p 0.06 and 0.02), rules 3-4 pass, all deviations toward more protection; T2 not submitted, the
  next run is the user's. ff3.2's `ff32_basin` stopped at step 76 (budget); its panel `w32_03`, shift, gaps and the
  epoch read queued themselves. ff3.1 NTL9 de novo COMPLETE (final frames finite, KE/1.5kT 1.003-1.007). Files:
  diffusion_dynamics.md Phase 3 T1, remote_jobs.md (header, owed items, line 4, job table), progress.md; cluster
  `WATCH_STATUS.md` and the bundle's `records/`.
* **11-04 19:25-19:45, Mac Studio: ff3.2 epoch 3 validated by the watch; Boltz T1b submitted.** `summary.txt` 19:13:
  `w32_03` chosen, LEVEL with released ff2.1 (+0.029 [-0.008, +0.066] p 0.065; `w32_02` not better at p < 0.05), so
  `checks/ff32_validate_20261104/submit.sh` ran (`plan`, then `apply` 19:26): `ff_3.2_w32_03` staged md5-checked in
  both trees, the mdCATH panel and TM, Peng's 32 arms and collector, glpG prep and REMD block 1. Then (user: "proceed
  with your recommendation on boltz") T1b: T1 with run seed 2 and an 8x same-trunk reference, rule written before
  its data (diffusion_dynamics.md Phase 3 T1b); `remc_boltz.py` observes 1,024 replicas per call (one call at T1's
  sizes; on T1's final frames chunks give equal flags and C-N, CA-RMSD within 2e-5 A, flags equal T1's record);
  `--test-only`, submitted 62293143. Files: plan.md Phase 38 item 11, diffusion_dynamics.md Phase 3, remote_jobs.md
  (header, owed items, The watch, line 4, job table), progress.md; cluster `WATCH_STATUS.md`, `addon/remc_boltz.py`,
  `addon/t1b.sbatch`; the bundle's `records/`.
* **11-04 23:48-00:05, Mac Studio: the watch found the Peng scorer cannot score ff3.1's de novo arms.** The
  collector's pass 49819033 submitted NTL9 de novo's score job 49829587, which printed `skipped` in 14 s: the scorer
  calibrates the residue offset on the native arm's first frame, and ff3.1's native replicas were deleted 11-04.
  A read-only probe showed the native arm's `input/<stem>.up` gives the same offsets on all 16 proteins, and equals
  frame 0 within 2.4e-7 A where replicas exist. The fix, reading the frame from `input/pos`, is proposed and awaits
  the user, since two running collectors use that scorer. glpG's ff3.2 prep completed (KE/1.5kT 0.998). Files:
  remote_jobs.md (header, owed items, line 2, job table), findings.md 10.51, progress.md; cluster `WATCH_STATUS.md`,
  `~/watch_probe/calib_input_probe.py`.
* **11-05 01:50-02:55, Mac Studio: ff3.2's mdCATH panel and TM read; Boltz T1b started.** `mdc_w32_03` (356 runs)
  folds 0.566 of 89 domains against released ff2.1's 0.488 (+0.078, p 9.7e-07) and ff3.1's 0.515 (+0.051, p
  5.6e-05); mean TM 0.739 against 0.727 (+0.012, p 1.1e-04) and 0.731 (+0.008, p 0.023), TM > 0.5 level. glpG's
  REMD block 1 started 01:55. T1b started 01:01, parity equal. ff3.1's BBL de novo COMPLETE 03:09 (final frames
  finite, KE/1.5kT 1.003-1.006), unscorable until the scorer fix. Files: plan.md Phase 38 item 11, remote_jobs.md (header,
  owed items, line 4, job table), progress.md; cluster `WATCH_STATUS.md`.
* **11-05 06:48-07:10, Mac Studio: Boltz T1b read, fails rule 2; nothing submitted.** 62293143 COMPLETED 06:03.
  `t1_read.py`: parity equal; residue 2's protected fraction +4.1 se against the limit 3.4 (0.9989 against 0.9968 in
  the 8,448-sample reference, counts 23 and 27 open); mean N +2.4 se, CA-RMSD +1.6 se, C-N level, no trend. T1's
  failing rows now read +0.7 and -1.1. 6 of the 7 rows over 2 se have the chains more protected. The excess sits
  after moves at 0.92-1.90 A; across one move the small sigma add protection and 7.67-15.95 A take it back, each
  within 2 se, and the mixed chain cannot isolate a sigma. ff3.2's 32 Peng arms all resubmitted at 05:57. Files:
  diffusion_dynamics.md Phase 3 T1b, remote_jobs.md (header, owed items, The watch, line 4, job table),
  progress.md; Mac `addon/runs/boltz_remc_t1b_ubqps_v100_62293143/` (`readout.txt`, `beside.py`, `beside.txt`);
  cluster `WATCH_STATUS.md`.
* **11-05 07:48-08:00, Mac Studio: the most-probable-state side chain (`ff21_bpmax_dt009`) converged at step 76;
  its evaluation submitted.** Link 49806998 reached 76 at 07:08 (24 of 24 proteins, no error line, step 76's
  KE/1.5kT at most 1.031). Gate 49846050 CONVERGED (smallest p 0.033, `dhb`); its mdCATH panel `mdc_bpmax_03`
  49846115 queued itself, configs in the mode. `checks/bpmax_fixedpoint_20261104/submit.sh plan` passed
  `--test-only` on all three jobs; `apply` at 07:49 submitted the 44-domain panel 49847544 and TM jobs 49847545 and
  49847546. ff3.2's collector at 07:30: 0 of 32 scored (first chunks only). glpG ff3.2 REMD block 1 at ff3.1's pace,
  no error lines. Files: plan.md Phase 39 items 4c and 5, remote_jobs.md (header, The watch, line 1, line 2's ff3.2
  row, job table), progress.md; cluster `WATCH_STATUS.md`.
* **11-05 09:48-09:55, Mac Studio: the most-probable-state side chain's 44-domain panel read.** `p44_bpmax_03`
  49847544 COMPLETED 09:14 (176 runs, KE/1.5kT 1.009-1.035, no non-finite output), TM 49847545 at 09:18: mean TM
  0.652, level with the BP control (0.650, p 0.51), above the fixed-state arm (0.623, +0.029, p 0.0023), below
  released ff2.1 (0.681, -0.028, p 3.6e-06). The mdCATH panel and its TM run on. Files: plan.md Phase 39 item 5,
  remote_jobs.md (header, line 1, job table), progress.md; cluster `WATCH_STATUS.md`.
* **11-05 10:50-11:55, Mac Studio: ff3.2's glpG preparation read; Boltz T1c submitted; ff3.2's validation checked.**
  User: "is ff3.2 keeping protein folded and glpG TM4 stable?", "decide for me the next move for Boltz. also make sure
  all ff3.2 validation jobs are submitted". glpG prep (`glpg_u3_prep_tm4.py`): TM4 DSSP through stage 7.0 mean 0.94,
  last frame 0.82 (ff3.1 0.95, 0.82); TM1 last frame 0.74 (ff3.1 1.00). REMD block 1 writes `/output` only at its
  end; its log at T 0.70, t 300-600: 94.7 backbone H-bonds, Rg 19.7 A (ff3.1 94.0, 19.9; ff2.1 without TM4 87.4).
  Boltz: the one-move test proposed at 06:55 was dropped before any run (T1b's per-move s.d. of N 1.26-2.28 gives
  0.020-0.036 se at 4,096 pairs against effects of 0.02-0.05); T1c, T1b with only the moves at 7.67 and 15.95 A,
  rule written first, `--test-only` passed, 62369218. ff3.2's validation: all 32 Peng arms queued, collector
  pending, glpG REMD running, mdCATH panel and TM read, staged in both trees; nothing missing. Files:
  diffusion_dynamics.md Phase 3 T1c, remote_jobs.md (header, owed items, The watch, line 2's ff3.2, line 4, job
  table), progress.md; `addon/t1c.sbatch`, `addon/submitted.log`; cluster `checks/ff32_validate_20261104/glpg_prep_tm4.txt`.
* **11-05 12:00-12:15, Mac Studio: ff3.2 is the target; ff3.0 and ff3.1 transitional** (user: "ff3.1 and ff3.0 no
  longer matters, if the ff3.2 validation is positive. they are transitional"). Recorded in plan.md "Architecture and
  Key Decisions". The ff3.1 scorer fix is deferred to the case where ff3.2's validation fails. ff3.2's Peng collector
  now compares against released ff2.1 alone (`checks/ff32_validate_20261104/collect.sbatch` `55615d2a`, backup
  `.bak_pre_ff21only_20261105`; `bash -n` and `--test-only` pass; the two-column `compare_peng.py` checked on ff2.1
  against ctrl38), since its paired summary keeps only proteins scored in every column. Files: plan.md,
  remote_jobs.md (owed items, line 2's ff3.2 row), progress.md.
* **11-05 14:40-14:50, Mac Studio: the most-probable-state side chain's mdCATH panel read; Phase 39 done.**
  `mdc_bpmax_03` 49846115 COMPLETED 13:55 (356 runs, KE/1.5kT 1.009-1.032, no non-finite output), TM 49847546 at
  14:23: mean TM 0.696, level with the BP control (0.702, p 0.25), above the fixed-state arm (0.667, p 1.5e-06), below
  released ff2.1 (0.727, p 7.1e-10), as on the 44 domains. ff3.2's collector now runs the ff2.1-only file (49854603).
  Files: plan.md Phase 39 item 5, remote_jobs.md (header, line 1, run table, job table), progress.md; cluster
  `WATCH_STATUS.md`.
* **11-05 16:40-19:50, Mac Studio: watch passes; ff3.2's first Peng score; T1c runs.** All 36 Peng chunks that ended
  16:13-16:29 resubmitted themselves with the node law. The 11:55 pass missed ff3.1 gpW de novo's `[bench] COMPLETE`
  (11:12); its final frames were read 16:45 (finite, KE/1.5kT 1.000-1.007 over 14 replicas). ff3.2 BBA native
  COMPLETE 19:30 (final frames finite, KE/1.5kT 1.002-1.007), scored 19:34: mean TM 0.360, median RMSD 7.70 A,
  against released ff2.1's 0.416 and 2.41 A; one arm, the paired verdict waits for `peng_ff32.txt`. Boltz T1c
  62369218 R from 16:57, `[check]` equal. Files: remote_jobs.md (header, line 2's ff3.1 text, line 4, job table),
  progress.md; cluster `WATCH_STATUS.md`.
* **11-05 19:55-20:40, Mac Studio: why ff3.2 loses BBA** (user: "check why ff3.2 fails so bad on BBA"). Read-only
  scripts on the BBA native arms (`checks/ff32_bba_20261105/`): at T 0.65 ff3.2 turns BBA into one helix over
  residues 1-26 (GLY8 alpha_R 0.86); the helix-minus-native energy gap moves -3.37 E_up, mostly the backbone burial
  (-2.49; its scale -0.301 to -0.136, halved only in the DSE-off runs `z32`/`w32`), the pair table (-1.44) and the
  glycine maps (-0.92). Other native arms level with ff2.1 but for weaker helix gains in WW domain and NTL9. Files:
  findings.md 1.39, progress.md.
* **11-06 00:40-01:15, Mac Studio: watch pass; Boltz T1c passes, T2a stops the large-move ladder.** Both sockets had
  expired; each reopened once (Duo). T1c 62369218 COMPLETED 22:38, `t1.npz` md5 equal, `t1_read.py`: rules 1-4 pass.
  T2a's rule was written into diffusion_dynamics.md before its read: the large moves are near-independent (N across a
  move r +0.029 and -0.011, limit 0.2), so a tilted ladder on them gains nothing over plain draws reweighted; no T2
  built or submitted, the options go to the user. midway2: ff3.1 NuG2 and hyp de novo and ff3.2 BBL native COMPLETE,
  final frames finite (KE/1.5kT 1.001-1.007); ff3.1 proteinG de novo resubmitted 49863650. Files: diffusion_dynamics.md
  Phase 3 T1c and T2a, remote_jobs.md (owed items, the watch, line 4, job table), progress.md; new
  `scratchpad/boltz_hdx/addon/t2a_read.py`; cluster `WATCH_STATUS.md`.
* **11-06 01:20-01:45, Mac Studio: does ff3.2's BBA failure reach mdCATH** (user: "proceed", to the offered check).
  `checks/ff32_bba_20261105/mdc_ss.py` on all 89 mdCATH domains: ff3.2 is 1% of residues more helical than the
  all-atom reference (ff2.1 0.6%), in alpha and mixed domains; the 26 beta-rich domains keep their helix level and
  fold better (+0.083, p 0.005). `bba_onset.py`: BBA's native left at 8,400 tu, near the panel's 10,000 tu, so the
  panel tests a BBA-pace takeover only at its end; Peng's long arms decide. The BBA failure is not a temperature
  shift (folded fraction falls on cooling). Files: findings.md 1.39, progress.md.
* **11-06 01:50-05:50, Mac Studio: Boltz T1d built, checked and submitted** (user: "Proceed with your recommendation",
  option (b) of T2a). Rule written into diffusion_dynamics.md Phase 3 T1d before any code ran. `remc.py` (`run`'s
  `level`, `denoise_from`'s `end`, `noised_level`), `remc_boltz.py --noised` with a split parity check, `toy_remc.py
  --noised --sigmas`, `t0_read.py --invariance`; backups `.bak_pre_noised_20261106`. Checks: T0's toy ladder rerun
  bitwise equal over 450 sweeps; the noised chain on the toy passes T0's rule 2 (-2.4, -2.3, -2.2 se, all one side);
  `--test-only` passed. T1d 62419935 submitted 05:42 (V100, about 5 h). Watch cron `5c4918e1` replaces `8935beb8`
  (reads T1d; T2 only to `--test-only`). Peng: 32 chunks resubmitted 02:44-02:59; ff3.1 cspA de novo and ff3.2 BBA
  de novo COMPLETE (final frames finite; cspA's 101-frame last chunk reads 0.988-1.006 from the restart transient,
  its previous chunk 1.003-1.006); ff3.2 BBL native scored 0.428 against ff2.1's 0.457. Files: diffusion_dynamics.md,
  remote_jobs.md, progress.md; `scratchpad/boltz_hdx/addon/` code and `runs/t1d_toy/`; cluster `WATCH_STATUS.md`.
* **11-06 06:35-11:40, Mac Studio: watch passes; Boltz T1d fails; Peng arms complete.** T1d 62419935 COMPLETED 11:09
  (wall 17,123 s), `t1.npz` md5 `d54b7d33` equal; `t1_read.py`: rule 1 passes (both `[check]` lines equal), rules 2
  and 3 fail (residue 73 +3.7 against 3.4; mean N 58.711 against 58.609, +3.9 se, toward protection, nearly twice
  T1b's +0.056), rule 4 passes (+0.7 se). Neither move pulls the raised chains back (change in N per move +0.003 at
  1.90 A, -0.003 at 3.95 A), so both carry the bias. Memory r +0.33 / +0.16; by T1d's rule nothing is built or
  submitted, and closing Phase 3 goes to the user. Peng: ff3.1 protein G (06:14) and protein L (09:21) de novo, ff3.2
  WW domain native and de novo (09:46, 09:55) and NTL9 native (11:23) COMPLETE, final frames finite over 14 replicas
  (KE/1.5kT 1.002-1.009); ff3.2 BBA de novo scored 0.378 against ff2.1's 0.411 (median RMSD 6.69 against 2.47 A);
  ff3.1's eight de novo score jobs skipped as expected (findings 10.51). Files: diffusion_dynamics.md Phase 3 T1d,
  remote_jobs.md, progress.md; `scratchpad/boltz_hdx/addon/runs/boltz_remc_t1d_ubqps_v100_62419935/`.
* **11-06 11:45-12:55, Mac Studio: the user's three asks (11:50).** (1) "hybrid preparation can be on CPU. Also modify
  the scripts in examples so they can run on GPU on midway2 or midway3. After that you can mark upside3 GPU as
  complete": `engine` switch in examples 01, 02, 03, 04, 05, 10, 11, 12, 13, 15 (06-09 use refused nodes),
  `16.MARTINI/run.py` `cluster` (its GPU block had depended across clusters), upside3.md 10.5; stub runs: `cpu`
  command lines unchanged, `gpu` adds only `--gpu` and the GPU job options; Phase 35 marked complete. (2) "move to
  next Boltz job and submit it for me": T2b, the N ladder on T1d's move (rule in diffusion_dynamics.md Phase 3 T2b;
  `t2_read.py` with a histogram MBAR, tested on synthetic ladders; `t2b.sbatch`, `--test-only` passed); its toy
  check (i) fails T0's rule 3 (tails 0.55-0.74 of the truth), so nothing is submitted and two diagnostic toy ladders
  run. (3) "fix the ff3.2, train it and re-do the validations": the interim Peng read (49882604; the scorer's method on
  finished frames, matched to released ff2.1) is level, dTM -0.004 [-0.017, +0.009], losses in BBA, WW and cspA de
  novo; `ff32_pcd` (persistent free chains, plan.md Phase 38 item 12, `ConDiv_pcd.py` `abb27c42`) built, checked
  (static test, compare_init, stub of after_training) and started (49883975, R 12:26; its first step writes the free
  frames). Lesson findings 10.52. Files: plan.md, diffusion_dynamics.md, upside3.md, findings.md, remote_jobs.md,
  progress.md, the examples; checks dirs `ff32_interim_20261106`, `ff32_pcd_20261106`; `scratchpad/boltz_hdx/addon/`
  (`t2_read.py`, `t2b.sbatch`, `t2_toy_diag.py`, `t0_read.py --stationary`, `runs/t2_toy/`).
* **11-06 13:10-13:30, Mac Studio: more proteins a step for ff3.2's fix.** The user relayed Jumper's suggestion (first
  read as more epochs: `epoch_drift.py` shows ff32_basin's side-chain table and sheet still drifting at 76 steps, and
  the convergence gate blind to it; then clarified as more proteins each round). `noise_scale.py`: ff32_basin's
  gradient noise scale 127-2,221 proteins by group in epoch 0 and beyond 456 after, so 24 a step is 1-16% signal
  (findings 1.40). The user chose 57 a step, 80 steps: `ff32_pcd` stopped after 1 step; `ff32_pcd57`
  (`make_pcd57.py`, `match_order.py`, `train_chain57.sbatch`, `after_training_ff32_pcd57.sbatch` reading epochs 1, 3,
  5, 7, 9; stub run of every branch) started 13:24 (49885329, 798 CPUs). plan.md Phase 38 item 12 revised before any
  data. Watch cron `13bea105` replaced `5c4918e1` at 13:10.
* **11-06 16:35-19:20, Mac Studio: watch passes; midway2-0261 excluded.** ff3.2 glpG REMD block 1 read late (TM4 0.977
  at T 0.70; the sacct-window miss, findings 10.53). ff32_pcd57 epoch 1's free replicas checked to start from epoch 0's
  final frames (`checks/ff32_pcd_20261106/pcd_start_check.py`, 57 of 57). Three score jobs died in a second on
  midway2-0261: its group-name lookup fails, so `module` init breaks `set -e` scripts (reproduced inside a running job
  with `srun --overlap`). Excluded at the user's word, `slurm/midway2.args` and `README.md` (backups
  `*.bak_pre_0261_20261106`), pending jobs updated, `--test-only` refuses it (remote_jobs.md §0d).
* **11-06 21:15-21:30, Mac Studio: Peng for the most-probable-state side chain** (user: "run the remaining validation
  for ff21_bpmax_dt009"; plan.md Phase 39 item 6). The other arms' Peng releases re-extracted md5-equal from their
  `epoch_01_minibatch_18`, so this arm's own was released as `ff_2.1_bpmax38` into the bpmax tree only. New in
  `checks/bpmax_peng_20261106/`: `release.sh`, `bench_bpmax.sbatch` (`bench.sbatch` with the tree's env), `submit.sh`,
  `collect.sbatch`, `static_check.py` (flag set in the config, patched binary, 100 tu finite). `--test-only` on all
  33, then 32 arms and the collector submitted 21:20, all R; every replica config carries `most_probable_state = 1`.
  Watch cron `d08a8402` replaced `7cb0bd8c` to cover them.
* **11-06 22:45-23:10, Mac Studio: ff32_pcd57 epoch 1 read** (plan.md Phase 38 item 12). Mechanism (a) passes
  (new `checks/ff32_pcd_20261106/mech_a.py`, read-only); the 44 fold 0.684 against released ff2.1's 0.615 at dt 0.015
  (`ff21_released`; `p44_ff21_released` is the dt 0.009 panel and is not the reference for these); the backbone burial
  scale fell to ff32_basin's step-76 value by step 16. Glycine rounds pulled to `scratchpad/handoff_20261102/gly/rounds/p57/`.
  Boltz T2b seed 5 COMPLETED 22:52; the read waits for seed 6.
* **11-07 00:47-01:05, Mac Studio: Boltz T2b read, fails.** Both seeds pulled md5-equal, `t2_read.py`: stationary,
  parity equal, calibration fails on I61 (0.47 of plain, limit 0.5), support fails (V70, I13, W45, L67 unopened); by
  the rule Phase 3 closes at plain resolution, left to the user (diffusion_dynamics.md Phase 3 T2b). 25 ff3.2 Peng
  chunks resubmitted at midnight; gpW native's `[bench] COMPLETE` of 22:50 read late, finite (findings 10.53 widened:
  a pass carries the time its sacct ran).
* **11-07 07:47-08:10, Mac Studio: ff32_pcd57 epoch 3 read** (plan.md Phase 38 item 12). The 44 fold 0.687 against
  released ff2.1's 0.615 and `w32_03`'s 0.645 (+0.042 [+0.001, +0.086]); shift +13.94, gap -1.150; glycine settled at
  -0.77; the backbone burial scale fell to -0.053 at step 32, reported to the user. Glycine rounds 3-4 of both runs
  pulled md5-equal to `scratchpad/handoff_20261102/gly/rounds/{p57,basin}/`. ff3.2 Peng: 11 more arms scored against
  released ff2.1 (20 of 32, 7 higher, 13 lower).
* **11-07 15:47-16:05, Mac Studio: ff32_pcd57 epoch 5 read** (plan.md Phase 38 item 12). The 44 fold 0.715, +0.099
  over released ff2.1 and +0.070 over `w32_03`; shift +14.58, gap -1.313; glycine -0.80 to -0.83; burial scale -0.035
  at step 48. Glycine rounds 5-6 pulled md5-equal to `scratchpad/handoff_20261102/gly/rounds/p57/`.
* **11-07 16:00-16:35, Mac Studio: Boltz T2c proposed, checked and submitted** (user: "Propose next boltz job for me
  and run it"; diffusion_dynamics.md Phase 3 T2c, Decision 8). One ladder per amide on its opening excess (I61, V70,
  I13, W45, L67), rule written first. `remc_boltz.py` gained `--targets` and `--lambdas`; new `t2c.sbatch`,
  `t2c_read.py`, `runs/t2c_check/` (excess check, synthetic jobs, rule 0's false-alarm rate). Check (ii) found rule 0's
  error wrong for four ladders (t with 3 dof) before any data; fixed. Four V100 jobs 62551232-62551235 queued.
* **11-07 22:00-22:15, Mac Studio: ff32_pcd57 against every flag** (user: "check if ff3.2 pass all the flags we set
  before", then "I meant the current training ff3.2"). `check_step.py` at each epoch end and the steric check, both
  read-only, never run on this run before: the H-bond margin was under +0.10 in epochs 1-4 and recovered to +0.139;
  STERIC FLAG on free non-glycine beta' + pPII' (0.0043 against 0.0042). Recorded in plan.md item 12; lesson
  findings 10.54.
* **11-07 23:50-00:05, Mac Studio: ff32_pcd57 epoch 7 read** (plan.md Phase 38 item 12). The 44 fold 0.701, +0.086
  over released ff2.1 and +0.056 over `w32_03`, level with epoch 5; shift +13.48, gap -1.660; glycine -0.82 to -0.84;
  the backbone burial scale stopped falling (-0.050, -0.038). Glycine rounds 7-8 pulled md5-equal to
  `scratchpad/handoff_20261102/gly/rounds/p57/`. Most-probable-state BBL native COMPLETE, finite.
* **11-08 05:50-06:15, Mac Studio: ff32_pcd57's step-80 build** (plan.md Phase 38 item 12). Built
  `checks/ff32_p57_validate_20261108/` from `ff32_epoch3_20261104` and `ff32_validate_20261104` (names, source check,
  `mdc_w32_03` as the TM panel's beside reference, glpG directory); rule (b)'s consequence written before the data
  (BELOW `w32_03` stops, the user decides). Dry run with `p57_07` standing in: `choose.py` reproduces epoch 7's
  numbers, `submit.sh plan` passes `--test-only` on every job. The backbone burial scale near zero in epoch 9
  (-0.0006 to -0.024). Most-probable-state BBA and BBL native scored, level with the BP control; its 30 other arms
  rolled into their fourth chunks.
* **11-08 06:50-07:05, Mac Studio: Boltz T2c read; ff32_pcd57 at step 80** (diffusion_dynamics.md Phase 3 T2c; plan.md
  Phase 38 item 12). T2c's four jobs COMPLETED 06:05-06:31, pulled md5-equal to `addon/runs/`, read by `t2c_read.py`:
  calibration fails (I61's ladders not stationary, +3.7 se against 3), so no amide is read; nothing submitted, the
  next step is the user's. ff32_pcd57's last link COMPLETED 06:49 at step 80; epoch 9's panel, shift and gaps queued.
* **11-08 07:10-08:10, Mac Studio: Boltz T2d proposed and submitted** (user: "Propose next boltz job and submit it
  for me"; diffusion_dynamics.md Decision 9 and Phase 3 T2d). Diagnosed T2c (`addon/runs/t2c_diag/diag.py`): the
  target ladders open once, at sweeps 617-1766, and I61 drifts toward protection at all 14 rungs; per-trunk plain
  counts are not overdispersed, so not trunk sampling. T2d asks the decade question of exact plain samples: 5 x 10^5
  more at seeds 113-612 through `e1_rare.py --unguided`, read by `t2d_read.py` (pooling by Fisher, Poisson bounds
  against experiment / 10), rule written first; checks: code md5-equal, reader right on five known sets,
  `--test-only` of 20 jobs. 12 submitted (62607649-62607660); 8 blocks owed to the `gpu` QOS's 12-job limit.
* **11-08 08:45-09:05, Mac Studio: ff32_pcd57 epoch 9 read, its choice and full validation** (plan.md Phase 38 item
  12). Epoch 9 folds 0.706 (+0.091 over released ff2.1, +0.062 over `w32_03`); `read.sh`: rule (b) passes, step 80
  chosen, above ff2.1; `submit.sh` plan (35 `--test-only`) then apply: `ff_3.2_p57_09` staged, mdCATH panel and TM,
  Peng's 32 arms and collector, glpG on Upside 3. Glycine rounds 9-10 pulled md5-equal. Nine Peng chunks (8 ff3.2
  de novo, ff3.1 ubiquitin de novo) resubmitted themselves 07:16-07:32; ff3.2 NTL9 de novo scored level with ff2.1.
* **11-08 09:10-09:25, Mac Studio: w32_03's validation cleared** (user: "you can clear the previous version of
  ff3.2's validation"). Cancelled its 8 de novo Peng chunks and collector (24 of 32 scored); deleted its Peng runs
  (30 GB), mdCATH runs (47 GB, after checking all 356 TM series that `ff_3.2_p57_09`'s TM reads are cached), glpG
  outputs (9.3 GB) and the staged `ff_3.2_w32_03` copies; readouts, jsons, logs, its 44-domain panel and params kept.
  upside3.md H2 now names `ff_3.2_p57_09`'s ubiquitin arms.
* **11-08 14:40-17:45, Mac Studio: watch through a /beagle3 stall; ff3.2's mdCATH read.** /beagle3 hung on both
  clusters about 14:30-15:15 (midway3 sessions returned nothing; recorded, no Duo relaunch). Boltz T2d's 62607650 went
  NODE_FAIL in it and Slurm requeued it under its id, so its run directory and 15 seeds are kept; the 145 seed files
  then were pulled md5-equal and sane. 24 bpmax Peng chunks resubmitted at their budget; 8 of 32 bpmax arms scored.
  `ff_3.2_p57_09` on mdCATH: mean TM 0.745 against released ff2.1's 0.727 (+0.018 [+0.012, +0.025]), above on every
  subset; `w32_03` +0.006 beside (plan.md Phase 38 item 12). ff3.1's last Peng arm finished (unscorable until the fix).
* **11-08 18:39-18:50, Mac Studio: Boltz T2d's /tmp fault found and fixed** (diffusion_dynamics.md Phase 3 T2d;
  findings.md 10.55). Three T2d jobs COMPLETED; four FAILED 127 at their next seed's start, each minutes after a
  neighbour on its node ended: midway3's epilog deletes every /tmp entry of the user at any job's end, and the venv
  was extracted to /tmp. `t2d.sbatch` now uses `$SLURM_TMPDIR` (`65da611b`); `--test-only` passed; the four went back
  into their run directories with seeds 413-487 (62625768-62625774). Seeds written before the stop are whole.
* **11-09 05:40-06:05, Mac Studio: Boltz T2d read** (diffusion_dynamics.md Phase 3 T2d). All 500 seeds pulled
  md5-equal, pooled (5 x 10^5, each seed 1,000 finite unsteered samples) and flagged by Phase 1's pipeline;
  `t2d_read.py`: not pooled with the 10^5 (D52, p 1.4e-4), read alone; V70 too stable by more than a decade, I13,
  W45 and L67 unresolved at a decade. Nothing submitted. Files: remote_jobs.md, diffusion_dynamics.md.
* **11-09 08:25-08:50, Mac Studio: Boltz E2a proposed and submitted** (user: "propose and run next BOltz job";
  diffusion_dynamics.md Phase 2 E2a, Decision 10). Ferrari et al.'s HDX-NMR set (Zenodo, public) gives 11 new
  site-resolved domains of the add-on's class; plain Boltz-2 at 5 x 10^4 a domain tests whether its own order holds
  there, the add-on's premise, against Upside FF2's 0.45. Rule fixed before data; inputs, MSAs and `nmr.sbatch`
  (`t2d.sbatch` with the input named) built; `--test-only` on 22; 12 submitted on RTX 6000 (V100s full), 10 owed.
  Failure and fix: a cut to 16 GB (from T2d's 5.7 GB peak) ran out of memory on midway3-0283, whose job tmp is in
  RAM; back to 40 GB and resubmitted before any seed (findings 10.55). Watch cron `0afd5100` replaces `87e5ce00`. Files: diffusion_dynamics.md, remote_jobs.md; scratchpad/boltz_hdx/addon
  (nmr_domains.py, prep_nmr.py, nmr.sbatch, inputs/nmr_domains.tsv), flags/reference/*_experiment.csv.
* **11-09 14:12-15:55, MacBook Pro: status pass, Boltz's faulted GPU, the backbone term's loss and the DSE diagnosed**
  (not the owner, no cron; plan.md Phase 38 item 13, findings 1.41). Status: every midway2 job since 07:39 COMPLETED,
  9 Peng arms COMPLETE with finite final frames; ff3.2's BBA native scored 0.388 mean TM (released ff2.1 0.416,
  `w32_03` 0.360). Boltz E2a 62687511 ran on midway3-0286's faulted GPU (`2F:00.0`, "System is not in ready state");
  cancelled, its block rerun as 62715074, 0286 excluded from every E2a job, its directory kept as `*.faulted_gpu`.
  Diagnosis, read-only (`checks/dse_diag_20261109/`): without the DSE the contrast alone takes the backbone solvation
  scale to zero in every epoch, persistent chains or not; the DSE path is the SI's; its size follows the unfolded
  frames' residual H-bonds, not their Rg. Its confirming run `dse_visits` (rule fixed first; smoke test passed)
  submitted at the user's word; no training until it is read. Failure and fix: an early reading of the term's sign
  called the free chains more buried; the engine's formula shows the term counts exposed backbone, and the record
  states that. Files: findings.md 1.41, plan.md Phase 38 item 13, remote_jobs.md, progress.md.
* **11-09 16:45-17:15, Mac Studio: the DSE run read** (user: "continue working on debug ff3.2"; plan.md Phase 38 item
  13, findings 1.41; not the owner, no cron). Before the verdict: the worker reproduces `ff32_pcd57`'s `main_worker`
  line for line, the reader implements the item's rule as written, and the 72 protein-visits are physical (Rg ladder
  from the native's to the SARW's, persisted frames finite, bond spread equal to the trainer's own frames;
  `checks/dse_diag_20261109/visits_sanity.py`). `dse_visits_read.txt`: verdict (c). h does not fall with persistence
  (-0.016 [-0.054, +0.021]); the colder bracketing replica holds most residual H-bonds (0.315 against 0.075), but the
  warmer alone balances the bb_env scale at 0.173 [0.042, 0.282], under 0.3. Nothing submitted. Then (user: "we have
  peng's training at ~/Documents") `Train(1).zip`'s 18 trainer copies read: all membrane-potential training with
  no DSE path, so the lambda of Peng's soluble run stays unknown. Files: findings.md 1.41, plan.md Phase 38 item
  13, remote_jobs.md, progress.md.
* **11-09 17:20-20:05, Mac Studio: λ's provenance, Peng's answer, ff3.3 built and submitted** (plan.md Phase 38
  item 14, findings 1.41, 9v). The DSE is Peng's (the paper's main text and SI p. 4, λ 0.3 "empirically set");
  Kleinmann's port has λ 0 in every commit, among local test settings, with Peng's original absent. Peng (via the
  user): ff2.1 used 0.3 after a training that did not start from scratch; raise it 0.1 -> 0.2 -> 0.3. The user chose
  the thirds schedule, a start from ff3.2 and the name ff3.3. Built `ff33_ramp57` (ff32_pcd57's step-80 checkpoint
  moved, driver `dcc9896f`); static test, stub run and `--test-only` passed. Submitted the first link 50022816 and
  the unfolded-state H-bond test on ff3.2 (50022814, 50022815). Files: plan.md, findings.md, remote_jobs.md,
  progress.md; cluster `checks/ff33_ramp_20261109/`.
* **11-10 13:35-16:10, MacBook Pro: status, Boltz E2a, and why our trainer pushes ff2.1 off** (not the owner, no
  cron; plan.md Phase 38 items 15-16, findings 1.42, 9v). Both sockets opened with one-shot copies of the expect
  scripts (the repo's re-send a rejected password). Status: every job since 15:00 11-09 COMPLETED; ff3.3's panels
  fall with lambda (`r33_11` 0.627, `r33_13` 0.547); ff3.2's glpG TM4 read 0.773; 7 Peng arms COMPLETE. At the
  user's word the 10 owed Boltz E2a blocks went in (62844278-62844287). Then, at the user's "assume DSE is not the
  problem" and "there has to be something else": the resume is clean; the backbone burial scale alone carries ff3.3's
  loss (7 reset panels); ff_2.0 (ff2.1 without its undocumented burial edit) folds level with ff2.1 and feels the same
  push at the control's step (fixed-point test), so the edit is not the cause; the same scale carries the control's
  drop from ff2.1 (-0.301 to -0.404 in epoch 0). Peng's FF2 trainer is not on GitHub or anywhere readable on the
  clusters (sweep of 18,085 files). Submitted the replica-resolved visit (50045068, readout 50045069) to find which
  sampling choice sets the scale's balance. Failures and fixes: a projection readout failed its own self-test
  (replaced before any ff_2.0 data, the rule revised in the plan); the replica worker's first smoke test lost its
  output to a cleanup-order bug (fixed, smoke passed). Files: plan.md, findings.md, remote_jobs.md, progress.md;
  cluster `checks/ff33_cause_20261110/` (README), `checks/peng_trainer_search_20261110/`.
* **11-10 16:20-16:55: the replica-resolved visit read** (not the owner, no cron; plan.md Phase 38 item 16, findings
  1.42). Reopened midway2's socket (one Duo push). Reviewed `replica_worker.py` against the control's `main_worker`
  (same temperatures, anneal, swaps, reweighting, DSE selection). Fixed `rep_read.py` before the readout ran: a protein
  with no unfolded frames under a DSE variant had been dropped from the NSE sum too (backup `.bak_nanprot_20261110`;
  functional check on the finished pickles, values not read). 50045068 (48 of 48) and 50045069 COMPLETED 16:42.
  Verdict by the rule: no sampling choice tested is the difference; the scale balances at +0.018 [-0.034, +0.071].
  Exploratory `rep_explore.txt`: the NSE grows steeply with the free replica's temperature. Read Peng's `Train(1).zip`
  copies and the port: restraint form the same (struck); Peng's 2022 membrane copies keep the native out of the swaps,
  ours, the port's and his 2019 soluble trainer do not (a weak candidate, untested); no archive of Peng's holds a DSE
  path (six searched by content); the scale's rate is the port's 0.025 a step, Peng's unknown. A first pair of
  zsh wait loops never ran their ssh (no word splitting of `$S`); replaced by a bash loop. Files: plan.md,
  findings.md, remote_jobs.md, progress.md; cluster `checks/ff33_cause_20261110/` (README, `rep_explore.py`,
  `rep_explore.txt`), `WATCH_STATUS.md` record line.
* **11-10 17:00-18:15: three candidates for the scale's balance tested, none is it** (plan.md Phase 38 item 16,
  findings 1.42). Searched Peng's archives by content: no DSE path anywhere. The native-out-of-swaps visit (50045830):
  +0.070 [+0.012, +0.138], not the difference. Code review of the whole path (engine derivative, config, anneal,
  exchange, integrator levels, contrast, Adam) found one undocumented departure from Peng's worker: it starts the
  hottest free replica and the SARW from a random coil, ours from the native. Its visit (50046185): +0.053 [-0.001,
  +0.120], not the difference. Built `ff21_coilstart` (the control's epoch 0 with that start; driver diff, moved
  checkpoint, static test, `--test-only`), link 50046188 R from 18:11, panel `cs_00` by item 16 (b)'s rule about
  08:00-10:00 11-11. Failures and fixes: shell quoting of inline scripts (moved to piped files); the static test's
  comparison helper (string arrays, the command-line record). Files: plan.md, findings.md, remote_jobs.md, progress.md;
  cluster `checks/coilstart_20261110/` (README), `checks/ff33_cause_20261110/`.
* **11-10 18:15-18:25: the control at Jumper's learning rate** (plan.md Phase 38 item 16 (c)). Jumper's FF1 factor is
  0.25 (`ConDiv_original.py` 539-546), which gives the SI's FF2 rates; ours and the port's is 0.5. Built
  `ff21_jumperlr` (driver one constant from the control's; checkpoint moved with both Adam rate fields halved, the
  rest checked equal), link 50046810 PD, panels `jl_00`, `jl_01` by the item's rule. A first mover run failed its own
  check (the copied solver shared the source's object); fixed with a deep copy before the change.
* **11-10 18:30-21:20: watch, the SI on held parameters, and a second code review.** Watcher on `ff21_coilstart` and
  `ff21_jumperlr` (both healthy, step 3 at 20:15); a 2cnj srun launch failure was the trainer's designed relaunch,
  so that line left the watcher's error patterns. The SI and main text train all parameters together, with no frozen
  group or stage (findings 1.42). Second review of the training code (plan.md item 17): no new defect; the known
  ledger gaps A29 and A30 are master's; the time step is the one protocol difference not yet measured on the balance
  (a dt 0.015 visit proposed, not started). Files: plan.md, findings.md, progress.md; scratch copies of the cluster
  code in the session scratchpad `review2/`.
* **11-10 21:20-21:40: the time step's visit** (plan.md Phase 38 item 18, user: "yes. proceed to test it"). Rule
  fixed in plan.md first; `checks/dt015_20261110/` built (worker = the dt 0.009 visit's with `time_step` 0.015 and a
  C-N bond record; readout paired against the dt 0.009 visit), checked by diff, bond indexing on natives and a dry
  run of the readout (reproduces +0.018, change 0); one fix before submission (the verdict reused a fresh bootstrap).
  Submitted array 50050586 and readout 50050587; added to the watcher. Files: plan.md, remote_jobs.md, progress.md;
  cluster `checks/dt015_20261110/` (README), `checks/coilstart_20261110/watch_state.sh`.
* **11-10 22:45: the time step's visit read** (plan.md Phase 38 item 18). 48 of 48 ran 21:39-22:38; by the rule
  it moves the balance toward 0.3 but is not the whole difference: +0.088 [+0.036, +0.151] against +0.018, change
  +0.069 [+0.031, +0.110], through the NSE. Health: KE/1.5kT 1.013-1.022 (1.006 at 0.009), no replica destroyed,
  C-N bonds to 2.15 A. Files: plan.md, findings.md, remote_jobs.md, progress.md.
* **11-10 23:00: ff3.3 epoch 15 read** (plan.md Phase 38 item 14). `r33_15` (lambda 0.2) folds 0.608: level with
  released ff2.1, +0.061 over `r33_13`, -0.098 under `p57_09`; shift -1.56. The scale stayed near -0.24 through
  epochs 14-15 and falls again under lambda 0.3 (-0.276 at step 131). Files: plan.md, remote_jobs.md, progress.md.
* **11-10 23:00-23:25: folded fraction dropped; ff3.3's validation armed** (user: "Then we do not care about it at
  all ... run its validation after the training", after Peng's paper was found to use TM-score, RMSD, H-bond
  cooperativity and melting-curve dG, no folded fraction). Plan item 14's 44-domain gate and epoch choice removed;
  the validated checkpoint is step 160. Built `checks/ff33_validate_20261110/` from ff3.2's validation (no gate,
  h added, `mdc_p57_09` beside); dry run passed `--test-only` on 38 jobs; hourly waiter 50052248 submits it at
  step 160. midway2 login2's /project writes hung intermittently (a copy and two log writes arrived empty; redone and
  md5-checked). Lesson findings 10.56; memory `group-folding-metrics` updated. Files: plan.md, findings.md,
  remote_jobs.md, progress.md.
* **11-10 23:30-00:00: does the folded-fraction drop matter in Peng's metrics?** (plan.md Phase 38 item 19; user: "the
  only metric we care about are Peng's and Jumper's metrices ... drop it"). Rule fixed first: TM-score and Ca-RMSD on
  the five 44-domain panels with trajectories, `r33_13` vs `p57_09` deciding. Built `checks/fold_vs_tm_20261110/`
  (`--test-only` passed), submitted 50052608; added to the watcher with ff3.3's waiter. The two ff2.1 runs keep
  going until it reads. Files: plan.md, remote_jobs.md, progress.md.
* **11-11 00:00-00:15: the drop is real in Peng's metrics** (plan.md Phase 38 item 19, user: the drop "might be a real
  defect in training"). From records: ff2.1 continued by our trainer at lambda 0.3 (`ff_2.1_ctrl38`) loses on Peng's
  benchmark, native -0.066 and de novo -0.128 TM against released ff2.1 (which reproduces his FF2), below FF1; and
  TM on mdCATH and the 44 domains. `fold_vs_tm` 50052608 cancelled as redundant (4 min run, partial outputs not
  results). The two ff2.1 checks are now read in TM; their references `ff21_released_tr`, `c9_00_tr` submitted
  (50052737-8). Files: plan.md, findings.md, remote_jobs.md, progress.md.
* **11-11 00:50-01:00: TM references ran.** `ff21_released_tr`, `c9_00_tr` COMPLETED 00:51 (176 of 176, parameter files
  md5-equal to the originals'). Submitted `tm_refs` 50053390 (the control's epoch-0 drop in TM and RMSD; caches both
  for `cs_00`, `jl_00`, `jl_01`); watcher restarted with it. Files: plan.md, remote_jobs.md, progress.md.
* **11-11 01:10: `tm_refs` read** (50053390, 9 min). The control's epoch-0 drop in Peng's and Jumper's metrics: mean
  TM -0.038 [-0.047, -0.030] (42 of 44 lower), Ca-RMSD +0.45 A (40 of 44 higher). References cached for `cs_00`,
  `jl_00`, `jl_01`. Files: plan.md, findings.md, remote_jobs.md, progress.md.
* **11-11 01:15-02:00: constrained DSE designed and recorded** (plan.md Phase 38 item 20; user: "what can we do to
  improve it?", "record this in md file, expect to run this on another computer tomorrow afternoon"). Wrote
  `constrained_dse.md` (design, rule, build and checks; copy `checks/dsecon_20261111/DESIGN.md`); measured ff2.1's h
  on the 48 training proteins (`h_target.py`: 0.145 [0.116, 0.179]); found the trainer's divergence files keep only
  the combined contrast, so the worker must record the three H-bond sums. CLAUDE.md gains a "run the constrained DSE
  test" entry. Files: constrained_dse.md, plan.md, CLAUDE.md, remote_jobs.md, progress.md.
* **11-11 09:45-10:25: status pass, TM read submitted, constrained DSE amended, built and submitted** (user: "show me
  remote job status", "Proceed with step 1", "double check that it is not picking results", "Yes, and go ahead and
  implement and test repaired training"). Every job since 22:00 COMPLETED; `cs_00` 0.418 and `jl_00` 0.464 folded
  (not read); `tm_bc` 50061363 reads both in TM (item 16 (b), (c)(1)), `r33_17` beside. The constrained DSE's rule
  amended before any data (binding criterion, controller check, exhaustive branches, screen, follow-up rule); built
  in `checks/dsecon_20261111/` (driver `3150a435`), every check passed after one fix (h summed in float64, as the
  target); first link 50061538 submitted. Files: constrained_dse.md (+ cluster `DESIGN.md`), plan.md items 20,
  remote_jobs.md, progress.md.
* **11-11 11:45-12:20: the constrained DSE's follow-up built and gated** (user: Peng asked for his code, no guarantee;
  "Do our planned debug assuming he cannot find it"). Searched midway for Peng's FF2 training output first (none
  readable; findings 9v). Built `ff32_dsecon57` (ff3.2 from step 80 under the constraint) in
  `checks/dsecon_k32_20261111/` as ff3.3 was built; checks passed; target visit 50063663 (114 proteins, released ff2.1)
  with readout 50063672; gate 50063731 starts the run if the screen allows; the validation at 160 generated from
  ff3.3's and dry-run (38 jobs). Faults found and fixed before they cost anything: a `r33_19` rename that would have
  compared the follow-up with itself (caught by the diff), the first gate dying under `set -eo pipefail`
  (findings 10.57). `dc_read.py` prints the TM class. Files: constrained_dse.md, plan.md item 20, findings 9v and
  10.57, remote_jobs.md, memory.
* **11-11 12:20-13:00: engine parity, test A** (user: "use the upside version at that time to test it", then "run the
  test that cost least amount of time"; plan.md item 21). Found Peng's FF2-era engines (`upside-v1.96` 2020-09,
  `upside-ff2.0` 2021-08) and Python 2 env; wrote the ConDiv configuration with each toolchain for 1hyp, 3ffy, 1kng;
  evaluated each with each engine. The backbone burial and H-bond derivatives that set lambda* are Peng's to 1e-7;
  the only differences are our two deliberate engine fixes (confirmed with our pre-patch builds) and Peng's 2020 CB
  placement. Files: `checks/engine_parity_20261111/`, plan.md item 21, findings 9v.
* **11-11 14:15-14:40: test B built and submitted; this morning's results read** (user, after the trading session: "do
  it now, record in md file so I can resume it from another computer"). Test B (plan.md item 21, rule fixed before
  the data): `replica_worker_peng.py` generated from item 16's worker with Peng's 2021 writer (his Python 2
  `run_upside.upside_config`), binary (his and our `run_upside` differ by a print) and library; smoke 50066492 ->
  array 50066493 -> readout 50066494. Read: the follow-up's target h 0.1394 [0.1189, 0.1630]; item 15's reset tests in
  TM (the backbone burial scale the one resolved block, +0.014 of ff3.3's -0.024). Submitted `tm_jl01` for item 16
  (c) part (2). Files: `checks/engine_parity_20261111/testB/`, `checks/fold_vs_tm_20261110/tm_jl01.sbatch`, plan.md
  items 15, 20, 21, remote_jobs.md.
* **11-11 15:25-15:35: monitoring pass** (user: "connect to remote and monitor corresponding jobs"). midway2's socket
  reopened (one Duo push); no job failed since 14:00. Test B's smoke passed (its KE/1.5kT 0.57 is a 300 tu run's
  warm-up; the array's 8000 tu visits read 0.998-1.016 on Peng's binary), array at 19 of 48. Read: `jl_01` in TM, level
  with the control's step 19, so Jumper's rate sets only the speed (plan.md item 16 (c), findings 1.42); ff3.3's h
  -0.109 [-0.197, -0.028] under ff3.2's and level with ff2.1's (item 14 rule (a) met). `ff21_dsecon` at step 9, lambda
  0.032 next. A read-only 10-min poller on the Mac reports changes. Files: plan.md items 14, 16, 21; findings 1.42;
  remote_jobs.md.
* **11-11 16:20: test B read** (plan.md item 21, findings 9v): on Peng's 2021 Upside the backbone burial scale balances at
  lambda* +0.047 [+0.002, +0.100], under 0.3, as ours (+0.018): the engine is not the difference. Health equal on both
  engines; the coinciding restrained-RMSD summaries checked per protein (they differ, median 0.013 A).
* **11-11 16:28: interim read of the constrained DSE screen at step 10** (user: "check the fix we proposed for DSE, is
  it working?"; plan.md item 20). Controller as designed (lambda 0 -> 0.038, h just above target, binding met); the
  scale goes the NSE's way (-0.301 -> -0.19, control -0.40); free replica Ca-RMSD -0.21 [-0.30, -0.11] A against the
  control on the same proteins and steps. Not the verdict (step 19). Files: `checks/dsecon_20261111/interim/`.
* **11-11 16:35-16:56: ff3.4 proposed, named and submitted** (user: "propose a ff3.4, with the fixed training, start
  from ff3.2", then "submit ff3.4 training now ... if the training repair validation passes, we will use ff3.4 as the
  official version"). ff3.4 = `ff32_dsecon57` (plan.md Phase 38 item 22): first link 50067557 and validation waiter
  50067558 submitted by hand as the gate would, gate cancelled; staged name renamed to `ff_3.4_dc32_19` (generator
  fixed for a header comment), dry run 38 of 38 `--test-only`, identical to the morning's but for the name. Files:
  `checks/dsecon_k32_20261111/` (README, `submitted.log`), plan.md item 22, remote_jobs.md, memory
  `project_ff34_constrained_dse`.
* **11-11 17:50-18:17: the constrained DSE rejected; early stopping adopted; new ff3.4 submitted** (user: "you are
  essentially just twisting lambda", then "Peng stopped early by picking a model with best result", then "stop ff3.4,
  and start a new ff3.4 with the updated training workflow, or apply it in in-training ff3.3 ... lambda=0.3, with best
  ff performance possible"). Lesson findings 10.58. Constrained ff3.4 cancelled 18:07. Early-stopping rule (the
  one-standard-error rule on the 44-domain panel in TM-score and Ca-RMSD) written before any data,
  `checks/early_stop_20261111/`. ff3.3: panels r33_16 and r33_18 and the selection es_r33 submitted 18:11. New ff3.4
  `ff34_lam03_57` = ff3.2's run at the SI's lambda 0.3 (driver one line from ff3.2's, checkpoint equal but paths,
  stub of every branch, --test-only), first link 50067712 submitted 18:16. Files: plan.md item 23, remote_jobs.md,
  memory `project_ff34_constrained_dse`.
* **11-11 19:51: ff3.3's early stopping read** (`checks/early_stop_20261111/es_r33.txt`; the rule's check passed):
  r33_19 chosen; validation TM rises 0.663 -> 0.680 over the lambda-0.3 epochs 16-19 and is level with released ff2.1
  (-0.003 [-0.011, +0.006]), under ff3.2 (-0.018). The stopping rule never fired. plan.md item 23.
* **11-11 21:41: ff3.3's mdCATH TM read**: r33_19 0.721 against released ff2.1's 0.727, -0.006 [-0.013, +0.000]
  (not wholly under 0); against ff3.2 -0.025. Agrees with the validation panel's ranking (findings 1.43). plan.md item 14.
* **11-11 22:11: the constrained screen read as a measurement** (`dc_read.txt`): partial; dc_00 TM -0.012 [-0.020,
  -0.004] against released ff2.1, +0.026 against the control, at a ninth of the control's DSE push; no weight in
  0.05-0.3 holds ff2.1 in place (findings 1.43). Its link finished; no job left on it.
* **11-11 22:25-22:38: ff3.3 extended; the panel check submitted** (user: "maybe extend ff3.3 a bit, as the score is
  still rising", "also do your check on 44-domain"). `ff33_ext57` (ff3.3 from step 160, its 0.3 stage to epoch 23;
  build, stub and --test-only passed) first link 50082451; `es_proxy` 50082449 (the panel against mdCATH and Peng,
  rule fixed before the comparison). Files `checks/ff33_ext_20261111/`, `checks/early_stop_20261111/`.
* **11-11 22:50: the 44-domain panel checked against the test suites** (`proxy_check.txt`): it predicts mdCATH
  (Spearman +0.98) but not Peng's native arms (every ff3.x 0.009-0.016 under ff2.1 there, level or above on the
  panel), the loss carried by small beta-containing proteins the panel lacks. findings 1.43, plan.md item 23.
* **11-11 23:49-23:57: ff3.3's whole curve; ff3.4 cancelled and deleted.** Epochs 11 and 15 scored (0.687, 0.681):
  each raise of the DSE weight costs TM and the stage recovers it, each ending near released ff2.1; figure
  `checks/early_stop_20261111/ff33_curve.png`. At the user's word ff3.4 (`ff34_lam03_57`) cancelled at step 9 and
  its run directory and epoch-0 panel deleted; scripts, logs and text readouts kept. plan.md item 23.
* **11-12 00:12: the rejected constrained ff3.4 deleted** (user: "delete ff3.4"): `ff32_dsecon57`'s run directory
  (39 GB) and its target visit's files; scripts, logs and readouts kept. `/project` 824 GB available after both ff3.4
  deletions (719 GB before them; GPFS showed the freed space a minute late).
* **11-12 08:20: ff3.3 extension epochs 20-21 read early**: 0.671 and 0.676 TM (released ff2.1 0.682), RMSD 2.33 and
  2.27 A; epoch 19 (0.680, 2.19 A) still the best at lambda 0.3. Link 2 running to step 192. plan.md item 23.
* **11-12 12:21: ff3.3 extension epoch 22 read early**: 0.668 TM (released ff2.1 -0.014 [-0.022, -0.006]); epochs
  20-22 all under epoch 19's 0.680. Epoch 23 and the selection at step 192 this evening. plan.md item 23.
* **11-12 16:05: ff3.3's stopping rule fired**: epoch 23 0.669 TM, -2.4 SE under epoch 19 after epoch 22's -2.7. The
  prepared second extension (`ff33_ext2_57`, epochs 24-27) not submitted, by the rule fixed before the read; the user
  decides. plan.md item 23, findings 1.43.
* **11-12 17:15: ff3.3's first Peng arm scored**: BBA native 0.366 TM against released ff2.1's 0.416 (RMSD median
  7.6 A against 2.4), as `w32_03`; one of 32, no reading before DONE. plan.md item 14.
* **11-12 17:47: which parts of ff3.3 fail BBA** (user: "Look at the ff and check what part of ff fail it"):
  findings 1.39's decomposition on ff3.3's BBA arm; the push toward the all-helix state (-4.12 E_up) is spread over
  the side-chain file (-1.87), the Rama maps (-1.56: glycine -0.90, others -0.66) and side-chain burial (-1.02); the
  backbone burial term is no longer in it. findings 1.44, plan.md item 24.
* **11-12 18:11: BBA swap test submitted** (user: "Run the test"): five BBA native arms on staged force fields,
  ff3.3 with one part from ff2.1 (glycine maps, sidechain.h5, environment.h5, whole map) and ff2.1 with ff3.3's map;
  rule fixed before they ran. `checks/ff33_bba_20261112/`, plan.md item 24.
* **11-12 23:52: BBL decomposition read** (job `bbl33`): ff3.3 keeps BBL's helices and moves its loop; the energy
  shift toward non-native (-5.36 E_up) is in the same three files as BBA's: sidechain.h5 -2.62 (hydrophobe coverage
  -2.17), rama.dat -1.75, environment.h5 -1.04. findings 1.44, plan.md item 24.
* **11-13 00:49: BBL swap arms submitted** (user: "submit three BBL arms"): ff3.3 with ff2.1's sidechain.h5, map, or
  environment.h5 on Peng's BBL native arm; rule (TM > 0.5 fraction, midpoint of ff2.1 and ff3.3) fixed before they ran.
  `checks/ff33_bba_20261112/` README step 3, plan.md item 24 step 3.
* **11-13 08:40: BBA swap test and ff3.3's TM4 read** (job `reads1113`, `swap_read.py`): ff2.1's map restores BBA
  in ff3.3 (folded 0.520 against ff2.1 0.848, ff3.3 0.088), ff3.3's map loses half of it in ff2.1 (0.424); the
  glycine-map, sidechain.h5 and environment.h5 swaps do not restore. ff3.3's glpG TM4 0.555 at T 0.70 (ff3.2 0.773,
  ff2.1 0.44). ff3.3's BBA de novo 0.361 (ff2.1 0.411). findings 1.44, lesson 10.61; plan.md item 24 step 2.
* **11-13 08:43: one more BBA swap arm, the extension deleted** (user: "Submit the run, and delete the extension").
  `ff_bba_33nongly21` (ff3.3 with ff2.1's map except ff3.3's central-glycine entries; only coil/dimer_pot differs
  between the libraries) submitted 50153802 with its rule fixed first; `ff33_ext2_57`'s run directory (296 MB, never
  run) deleted, `ff33_ext57`'s checkpoint md5 unchanged. plan.md item 24 step 2b, item 23.
* **11-13 10:45: `ff33_dsemap57` built, not submitted** (user: map in the DSE, start at ff3.3's lambda-0 point,
  published side-chain rate, no phi/psi springs, SARW temperature fixed): trainer edits by `make_dsemap.py` on
  ff3.3's driver and ff3.2's rama_basin; static tests and the after-training stub pass; first link --test-only.
  findings 10.63-10.64; plan.md item 25.
* **11-13 11:00: `ff33_dsemap57` running** (first link 50156016, steps 80-96 at lambda 0): the schedule now has two
  lambda-0 epochs before ff3.3's ramp (user), to step 176; driver replaced in place before any job, checks re-passed.
* **11-13 13:48: ff3.1, ff3.2 and ff3.3 validations cancelled** (user: their training is flawed): 32 jobs, the 26
  running ff3.3 Peng arms, ff3.2's last 3 de novo arms and the three collectors; none resubmitted. Peng stopped at
  6 (ff3.3), 26 (ff3.2) and 23 (ff3.1) of 32 scored; runs, scores and read-outs kept. Kept running: `ff33_dsemap57`,
  the swap arms, bpmax. remote_jobs.md, plan.md Phase 37 item 4 and Phase 38 items 12, 14.
