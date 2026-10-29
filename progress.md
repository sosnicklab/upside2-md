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
  in compact non-native globules, and vectorial growth equals bulk (findings 11e, upside3.md 5.2).
  Study data in `scratchpad/large_protein_20261019/`.
* **The `environment_coverage_hb` typing defect, 10-19** (findings 3.13, upside3.md 5.4). Fixed in
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
  A102, A126). The closing criteria (upside3.md 2.2 item 5) were met 10-22. The Upside 3 tree
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
  BioEmu2 and Boltz HDX readouts byte for byte. Present under every stack, not
  investigated: index oracle flags on 6 datasets (`hb_surf_membrane_potential` and
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
  kcal/mol. Files: findings.md 5.7, upside3.md 4.1 (design and steps M2-M6, T1) and 7 item 5, plan.md Phase 23
  item 5. No cluster job submitted; three arms copied to the Mac read-only. Then "record the plannings. We will run
  it later": M2, M3 (Peng's FF2 protocol from SI pp. 12-13, about 570 node-hours) and M6 specified in upside3.md
  4.1, the item added to remote_jobs.md "Resume here" step 7, `fetch.sh` added, and the handoff bundle's records
  and M1 tools refreshed (old records in `records.bak_20261028_1150/`).
* **10-28 11:05-13:00, MacBook Pro, disordered regions** (user: "for upside3, investigate the possibility of
  enhancing intrinsic disordered region simulation accuracy"; then "record everything in upside3.md. We will
  implement this later"). Code: training never sees a disordered chain; the DSE target is ConDiv's SARW (backbone
  excluded volume, Rama only). Literature by sub-agent (`scratchpad/idp_20261028/lit/review.md`): no Upside-IDP
  comparison with experiment exists; a 16-entry SAXS panel. User corrections: Upside's electrostatics are implicit in
  the trained terms (findings 10.42), and no explicit electrostatic term is added (double counting). B0 launched on
  the Mac 11:56 (ff2.1 x 6 and SARW x 3, then ff3.1 x 6; 30,000 tu from random coil); alpha-synuclein pilot: ff2.1
  35.6 and 31.0 A against 35.5 measured, SARW 40.1, Rg tau 340-1,300 tu. Files: upside3.md 4.2 (with 4.1 M6, 7
  item 6, 8), plan.md Phase 34, findings.md 10.42.
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
  * Recorded in upside3.md 4.2, 7 item 6 and 8, and plan.md Phase 34.
  * Then, on "record this ... I will work on it later in another session":
    * Added 4.2 steps B1 (ACTR at 298 K) and B2 (sampling).
    * Added a "Resume here" for a cold session.
    * Wrote B0's answers into I2 and I3.
    * Flagged that ACTR's SAXS source is unchecked. The literature sub-agent once sent the user's email in a Crossref
  User-Agent header, against its instructions.
