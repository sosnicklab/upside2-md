# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`. The full log up
to 10-07 12:00 is in git history.

---

## Summary to 2026-10-06

* **09-18 to 09-24, trainer.** The ConDiv port then in use trained FF1's form: spline burial, no
  backbone desolvation, no DSE (findings 9t, 9u). It was replaced by Kleinmann's port of Peng's FF2
  trainer, restored to the SI. From ff2.1, 8 of 9 groups held a fixed point (9v). AWH on 40 glycine
  dipeptides finished, and ff14SB agrees with ff99SB-ILDN to 0.045 nats (9r).
* **09-27 to 09-30, round 2 (basin offsets).** ff_3.0 was released 09-30 02:25 and withdrawn at
  20:34 (findings 1.14-1.15). Its glpG chains ran on seeds left with the FF1-form ff_3.0's coverage
  tables (3.11).
* **10-01 to 10-04, round 3 (`ff30_glyhb`).** AWH glycine map, glycine H-bond offsets in
  `HBondEnergy` (12-entry configs bitwise unchanged), side-chain rate / 10, the selection panel and
  release by checkpoint selection. The panels lost helix and the offsets went left-handed (1.17);
  stopped 10-04 21:58.
* **10-02, cleanup.** `training/` was cut to five files and made site-neutral. `findings.md` was
  condensed, and CLAUDE.md's compile and Slurm sections were rewritten (user).
* **10-02 to 10-05, rotamer-BP validation figure for Tobin.** Finished; it is in
  `/beagle3/.../bp_validation/static/`. The deck is not sent.
* **10-04/05, round 4 designed.**
  - Library: BioEmu-fitted glycine library (`rama31.dat`); push probe d -0.014 [-0.041, +0.012];
    non-glycine maps kept NDRD (1.19).
  - Runs: ff30_bio and ff30_gdepth started on caslake; automatic validation deployed.
  - Poly-Gly reference started (Phase 12).
* **10-06, restart at dt 0.009** (user). dt 0.015 destroyed free replicas, and both epoch-0
  panels lost folding and helix (1.21-1.22). The dt 0.009 runs and the matched control
  ff21_ctrl_dt009 started on midway2. TM4 moved to 12 seeds. Those counts were later found invalid
  (3.11).

## 2026-10-07

* **Frozen arm** (user). `ff30_bio_fz`, `ff30_gdepth_fz` and the control `ff21_ctrl_fz` run with
  `hb`, `dhb`, `hbg` and `sheet` at learning rate 0. Initial force fields byte-identical to their
  twins'. ff30_bio_fz started 10:05; at step 0 `hbond.h5` and `sheet` were md5-identical to
  ff2.1's.
* **How Jumper and Peng trained the H-bond** (user question; findings 1.23). Read Jumper's thesis,
  both 2018 papers, Peng's 2022 SI and every trainer.
  - No soluble trainer froze or staged `hb`. It trained jointly at a small rate, and both authors
    stopped early. Without DSE, H-bonds got stronger; Peng's DSE pulls them weaker.
  - Our `hb` rate is twice the SI's.
  - ff_2.0 to ff_2.1 kept `hbond.h5`, `sheet` and `bb_env.dat` byte-identical.
  - Corrected 9e and 1.17: ff2.1's H-bond energies came from FF2's ConDiv.
* **TM4 test defect found and fixed** (user approved; findings 3.11).
  - Defect: `patch_glpg.py` never wrote the two coverage tables. The live seeds carry the FF1-form
    ff_3.0's (rewritten so by the 09-30 deploy), so every TM4 count to date is invalid. At the seed
    frame the rotamer energy is -19.3 as tested against -121.7 with bio_start's own tables.
  - Fixed script deployed: local, `$P/training/` and `tm4_local/scripts/`, through login1 (login2's
    `/project` hangs). The pristine-seed gate passes as before. Notes for other computers are in
    remote_jobs.md and the cluster README.
  - Local inputs re-patched; pre-fix runs moved to `runs_precov_20261007/`.
  - The other Claude session (the watch) coordinated throughout.
* **No remote job needed cancelling.** Trainings, panels and Peng arms build their inputs with
  `upside_config` from the force field's own files. Only patched glpG inputs were affected.
* **TM4 reruns on fixed inputs.** bio_start ran 11:24-~13:10, then b9_01 to ~15:00. c9_00 was
  extracted 11:59 and b9_00 fetched.
  - The driver was replaced twice (a kill tested on a dummy first; running sets untouched). The
    last, `rerun_cov3_20261007.sh`, stops after b9_01 for the 16:00 shutdown. The rest of the queue
    is in remote_jobs.md "Handoff".
* **User correction: TM4 is judged only in the hybrid** (findings 10.19, memory
  `tm4-judged-in-hybrid-only`).
* **Handoff** (user: this computer shuts down at 16:00).
  - plan.md rewritten around the TM4 repair job (530 -> ~190 lines). remote_jobs.md header and
    "Resume here" became one handoff section, and the §1 TM4 block was replaced.
  - Pre-fix TM4 claims corrected in plan.md, remote_jobs.md and findings.md (1.14-1.22, 3.11, the new
    12c); architecture.md and GLY_sym.md moved from round 3 to round 4.
  - findings 1.13-1.22 and 9c-9w condensed from drafts by four agents, each checked and spliced
    (findings.md 4,704 -> 4,325 lines; the four files together 6,817 -> 5,645). Every number in
    each draft traces to its original. One substantive fix: the backbone term's engine derivative
    is `scale` only (9u now agrees with 9v).
* **13:00-13:20, SI-rate arm and the frozen runs** (user). `ff30_bio_si` and `ff30_gdepth_si` were
  submitted: the two designs with every group but rot at the SI's learning rate, otherwise their
  dt 0.009 twins. Built by `checks/si_init_20261007/` scripts, with initial force fields
  byte-identical to the twins'. The user left the frozen runs to my judgement:
  - ff30_bio_fz kept (running; with b9 and bs, the bio design gets three H-bond rates);
  - ff30_gdepth_fz and ff21_ctrl_fz cancelled before they started (queue start 10-08 10:44 and
    23:59).
  - `submit_new.sh` maps `bs_EE` and `ds_EE` (backup `.bak_pre_si_20261007`).
* **bio_start on fixed inputs (12 seeds, 11:24-13:15):** unwound 5, flipped 4, every seed with
  KE/1.5kT 1.001-1.011. The same force field on pre-fix inputs gave 4 and 1 (Fisher p 1.00, 0.32).
  b9_01 started 13:19:30 after a waiter deadlock (findings 10.14) and ends about 15:10. The watch
  records results in findings 1.24.
* **13:35 ff21_ctrl_fz resubmitted** (user): with c9 it is the cleanest test of the hypothesis behind
  the freeze, H-bond drift destabilizing TM4, with no glycine change. That hypothesis rested on
  pre-fix TM4 counts and an unresolved b9_00 difference (p 0.41), and is untested on fixed inputs.
* **14:00-14:10 poly-Gly production** (user approved). The collapse was read with gmx gyrate,
  polystat and mindist -pi: compact within 4 ns, reopening, image distance >= 5.2 nm. Built four
  7.5 nm dodecahedron starts (~28,950 atoms) and submitted the self-chaining `prod.sbatch` as
  49204428. The collapse ran at 14 ns/day, not the planned 55; 1 us may take about a month.
* **b9_01 on fixed inputs (12 seeds, 13:19-15:14):** unwound 6, flipped 5, against bio_start's 5 and
  4 (Fisher p 1.00). Mean last-block TM4 0.83 against 0.90, Mann-Whitney one-sided p 0.23.
  Unresolved; slightly toward unwound. The local queue ended for the shutdown; "continue jobs" on
  the next computer starts at remote_jobs.md "Resume here" (pointer in CLAUDE.md).
* **16:20-16:35, Mac Studio takes over ("continue jobs").** midway2 socket live (landed on
  login1). Link 49194446 (ff30_bio_dt009) ended normally at step 53, 16:02; successor queued.
  New steps all healthy (24 of 24, KE/1.5kT <= 1.013). Margins: b9 +0.015 (step 53), d9 +0.110
  (step 29, nearing the +0.10 hold), c9 +0.127 (step 26). Nothing else started.
  - Local TM4 test brought up to date: fixed `patch_glpg.py` installed, pre-fix `runs/` and
    `patched/` moved aside, eight queue inputs re-patched, `c9_00` and `runs_cov/` fetched.
  - Parity: the re-patched bio_start input equals the MacBook Pro's in all 150 `/input` datasets,
    and a 200 tu run reproduces its seed-1 log frame for frame.
  - New driver `scripts/tm4_queue.sh` runs `tm4_queue.txt` one set at a time; ff21_released started
    16:29. Watch cron `cad6ca59` hourly at :23.
* **17:46-18:40, watch pass.** ff30_bio_si started 16:34; its rates are confirmed from the checkpoint
  and the first step's changes. ff30_gdepth_si started 17:20. Panel c9_00: folded 0.390, below
  b9_00 and d9_00 (findings 1.22). TM4 ff21_released (Mac Studio, 16:29-18:26): unwound 9, flipped
  5 against bio_start's 5 and 4, p 0.21, not resolved (findings 1.24). New
  `scripts/tm4_compare.py` reproduces the b9_01 table exactly; `events_vs_tm4.py` reads logs
  beside the runs (both copied to the cluster).
* **19:47-20:35, watch passes.** All steps 24 of 24, KE/1.5kT at most 1.016. ff30_gdepth_dt009's
  margin fell to +0.097 at step 33: hold, user told. TM4 c9_00 (Mac Studio, 18:26-20:25): unwound
  7, flipped 1 against ff21_released's 9 and 5, not resolved (findings 1.24); set and table copied to
  `tm4_local/runs_cov/`. gdepth_start started 20:25.
* **21:45-21:55, frozen arm completed (user).** ff30_gdepth_fz resubmitted unchanged as 49206070
  after checking its inputs against ff30_gdepth_dt009's (plan.md Phase 11). ff21_ctrl_dt009
  finished at its target (21:25); c9_01 extracted, patched, queued for TM4, panel 49206076.
  ff21_ctrl_fz started 21:16; ff30_bio_fz resumed from step 4 at 21:46. ff30_gdepth_dt009's first
  link timed out in step 35 (remote_jobs.md §8); its successor resumes from step 34.
* **22:10-22:55, watch passes.** ff21_ctrl_fz freeze check passed at step 0. TM4 gdepth_start
  (20:25-22:25): unwound 7, flipped 2, not resolved from ff2.1 or bio_start (findings 1.24); copied
  to `tm4_local/runs_cov/`. polygly production started ~22:41; minimisation stopped at 5000 steps
  with Fmax 223-1411 on protein atoms, as the collapse's did; equilibration running.
* **23:47, watch pass.** Panel c9_01: folded 0.405, helix -0.029 (c9_00 0.390; findings 1.22).
  polygly equilibrated 23:03 at 33.7 ns/day per replica; image distance 3.5-4.9 nm over the first
  1.1 ns. Frozen runs hold +0.1919; ff30_bio_si +0.125 at step 13.
* **00:10-00:35 10-08, handoff prepared** (user: the jobs move to the work computer ~10:00).
  remote_jobs.md "Handoff" rewritten for the Mac Studio to work computer move (last-pass steps,
  queue state, setup with an exact parity check), header, §1 and the stored watch prompt updated
  (ff30_gdepth_fz added). Watch probes saved on midway2 as `~/watch_probe/stock.sh` and `steps.sh`.
  ff30_bio_dt009 link 2 started 00:12 (step 54); ff21_ctrl_fz's link 49204264 ended in NODE_FAIL at
  00:12 in step 2, its successor queued.
* **00:26, TM4 d9_00:** unwound 5, flipped 2 against gdepth_start's 7 and 2, not resolved (findings
  1.24); copied to `tm4_local/runs_cov/`. b9_00 started 00:26.
* **02:23-03:00, TM4 b9_00 and three epoch ends.** b9_00: unwound 4, flipped 4 against bio_start's
  5 and 4, not resolved (findings 1.24); copied to `tm4_local/runs_cov/`. fp_e00 started 02:23.
  `submit_new.sh` submitted panels b9_02, d9_01, bs_00 (49207052-54); the three are extracted to
  `checks/r4_epochs/`, patched and queued (d9_01 first).
* **04:19-05:00, TM4 fp_e00; bz_00.** fp_e00: unwound 5, flipped 4 against ff21_released's 9 and 5,
  not resolved (findings 1.24); copied to `tm4_local/runs_cov/`. c9_01 started 04:19. bz_00
  extracted (hbond.h5, sheet identical to bio_start's), panel 49207340, patched and queued.
* **06:15-07:00, TM4 c9_01; ff30_bio_si hold.** c9_01: unwound 10, flipped 6 against ff2.1's 9 and
  5 and c9_00's 7 and 1, not resolved (findings 1.24); copied to `tm4_local/runs_cov/` with a
  second table against c9_00. d9_01 started 06:15. ff30_bio_si's margin fell below +0.10 at step 26
  (+0.097); user told, it trains on.
* **08:12-08:50, TM4 d9_01; ff30_bio_dt009 node failure; ds_00.** d9_01: unwound 6, flipped 4
  against gdepth_start's 7 and 2 and d9_00's 5 and 2, not resolved (findings 1.24); copied to
  `tm4_local/runs_cov/`. bs_00 started 08:12. ff30_bio_dt009's link 49194448 ended in NODE_FAIL at
  08:32 in step 64 (batch host midway2-0027, which failed ff21_ctrl_fz's link at 00:12); successor
  49206580 resumes from step 63. ds_00 extracted, panel 49207993, patched and queued.
* **09:10-09:25, 24-seed frozen comparisons; the Mac Studio keeps the jobs.** User approved 24 seeds
  for bz_00 against b9_00 and cz_00 against c9_00, and then kept the jobs on the Mac Studio with
  another computer standing by. `tm4_compare.py` reads seeds 1-24 from `runs_cov/` and `runs/`
  with a per-set n (reproduces the earlier tables; mock 24-seed test as expected); `tm4_queue.sh`
  takes `<tag> 13 24` lines (stub test), installed by rename under the running driver, which exits
  after bs_00 when a waiter starts the new one (both scripts copied to `tm4_local/scripts/`, backups
  `.bak_pre_24seed_20261008`). remote_jobs.md "Handoff" rewritten as one owner at a time, with a
  live `tm4_local/WATCH_STATUS.md` on the cluster that each watch pass rewrites.
* **09:25-09:40, midway2 node exclusion law** (user; plan approved). `/project/trsosnic/yinhan/slurm/
  midway2.args` is the one list (midway2-0027 added for two NODE_FAILs as batch host), with
  `README.md` (each node's record), `update_pending.sh` and `sbatch_wrapper` (linked as
  `~/bin/sbatch`). 13 `slurm.args` symlinked to it; `submit_remd.sh` and `bench.sbatch` read it;
  `train_chain.sbatch` and training `README.md` comments point at it (backups
  `.bak_pre_law_20261008`). Tests with `--test-only` and stubs all as expected; all 13 pending jobs
  updated, dependencies intact. remote_jobs.md §0d and CLAUDE.md "Default Cluster" carry the rule.
* **10:48-11:05, TM4 bs_00.** bs_00: unwound 6, flipped 5 against bio_start's 5 and 4 (p 1.00)
  and its twin b9_00's 4 and 4 (p 0.68, 1.00), not resolved (findings 1.24); copied to
  `tm4_local/runs_cov/` with both tables (md5 checked). The waiter started the seed-range driver at
  10:08 and bz_00 runs (ends ~12:05). check_step on 7 new steps: 24 of 24, KE/1.5kT at most 1.014;
  no new hold. The permission check refused `submit_new.sh`; it had nothing to submit (every epoch
  end has its panel). User (11:05): the watch keeps submitting panels on its own.
* **11:35-11:50, bs_01 panel and TM4; ff30_gdepth_si node failure.** `submit_new.sh` (user
  approved 11:05 and 11:40) submitted bs_01's panel 49208587, the first panel under the exclusion
  law (ExcNodeList carries 0027). bs_01 extracted (rama.dat identical to bs_00's, hbond.h5 moved),
  patched and queued behind `b9_00 13 24`. ff30_gdepth_si's link 49204256 ended NODE_FAIL at 11:07
  in step 21 on midway2-0027 (third failure there today; it started before the law); successor
  49205173 resumes from step 20, estimate 10-09 20:34. Recorded in remote_jobs.md, the law's
  README and `WATCH_STATUS.md`.
* **12:05-12:20, TM4 bz_00 seeds 1-12; `tm4_compare.py` reads a first half.** The script refused
  bz_00 while its seeds 13-24 ran; it now takes an optional third argument N (12 or 24) that reads
  seeds 1-N of each set (link dir `SET_vs_REF_N/`; backup `.bak_pre_nseeds_20261008`, md5
  `46984192...`; reproduces b9_00's table with and without N; bz_00 without N still refused).
  bz_00: unwound 3, flipped 0 against 5 / 4 (bio_start) and 4 / 4 (b9_00); flips one-sided p 0.047,
  two-sided 0.093, unwound not resolved; s8 unwinds to 0.45 with GLY143 flipped in blocks 1-3
  (findings 1.24). Seeds 1-12 and both tables copied to `tm4_local/runs_cov/`, the script to
  `tm4_local/scripts/` (old one backed up there).
* **14:00-14:20, TM4 bz_00 at 24 seeds; bz_01.** bz_00's seeds 13-24 unwound 5 and flipped 4 of
  12, bio_start's counts; at 24 seeds 8 and 4 against bio_start's 5 and 4 of 12 (one-sided p 0.45,
  0.24), last-block TM4 0.897 against 0.899. The first half's zero flips (p 0.047) did not repeat
  (findings 1.24). Seeds 13-24 and `tm4_compare_cov_bz_00_24.txt` copied to `tm4_local/runs_cov/`
  (48 files, md5 checked). bz_01 (step 37, 13:58): panel 49210456 (carries the law), extracted
  (hbond.h5, sheet, rama.dat identical to bz_00's), patched, queued behind bs_01; its patched input
  has bio_start's H-bond energy and Rama map exactly. b9_00 seeds 13-24 run until ~16:00.
* **15:56-16:10, TM4 b9_00 and bz_00 at 24 seeds each.** b9_00's seeds 13-24 unwound 3 and flipped
  1; at 24 seeds 7 and 5 against bio_start's 5 and 4 of 12 (two-sided p 0.48, 0.44). bz_00 against
  b9_00, 24 seeds each: 8 against 7 unwound, 4 against 5 flipped, p 1.00; TM4 0.897 against 0.914.
  TM1 is higher in bz_00 (0.891 against 0.813, added Mann-Whitney p 0.003, both halves the same
  way), and bz_00 against bio_start in TM1 is p 0.25 (findings 1.24). Tables
  `tm4_compare_cov_bz_00_24_vs_b9_00.txt` and `tm4_compare_cov_b9_00_24.txt` and b9_00's seeds
  13-24 copied to `tm4_local/runs_cov/` (md5 checked). bs_01 runs from 15:56.
* **17:47-18:10, d9_02; TM4 bs_01.** d9_02 (ff30_gdepth_dt009 step 56, depth round 3, dL - dR
  +0.602): panel 49212678 from `submit_new.sh` (carries the law), extracted, patched, queued last.
  bs_01: unwound 5, flipped 2 of 12, against bio_start's 5 and 4 (p 1.00, 0.64) and b9_01's 6 and 5
  (p 1.00, 0.37); last-block TM4 0.879 against b9_01's 0.833. s9 unwinds to 0.51 with GLY149
  flipped, TM1 0.32 and a total-potential jump of 18081 (findings 1.24). Seeds and tables copied to
  `tm4_local/runs_cov/`. bz_01 runs from 17:53.
* **18:47-20:10, TM4 bz_01.** bz_01 (ff30_bio_fz step 37, H-bond and sheet frozen): unwound 6,
  flipped 2 of 12, against bio_start's 5 and 4 (p 1.00, 0.64) and b9_01's 6 and 5 (p 1.00, 0.37);
  last-block TM4 0.846 against 0.899 and 0.833, block means falling as b9_01's do; s12 ends at
  0.22, the lowest on fixed inputs. The freeze does not stop the epoch-1 lean. A dataset-by-dataset
  comparison of the patched inputs shows bz_00 and bz_01 differ from bio_start in the rotamer pair
  interactions as well as the two coverage tables; findings 1.24 had listed only the coverage
  tables, corrected. Seeds and tables copied to `tm4_local/runs_cov/`. b9_02 runs from 19:49.
* **20:47-21:50, TM4 b9_02 and three panels.** b9_02 (ff30_bio_dt009 step 56): unwound 8, flipped 4
  of 12, against bio_start's 5 and 4 (p 0.41, 1.00) and b9_01's 6 and 5 (p 0.68, 1.00); last-block
  TM4 0.830 against 0.899 and 0.833. Its loss comes in the last block (block means 0.99, 0.98,
  0.93, 0.83), and two of the eight unwound end at 0.889 (findings 1.24). Panels b9_02, d9_01 and
  bs_00 (34 domains): folded 0.448, 0.457, 0.457, each dominated by gdepth_start in helix; b9_00
  does not dominate bs_00 (findings 1.21). Seeds and tables copied to `tm4_local/runs_cov/` (md5
  checked). ff30_gdepth_fz's first link started 21:34; its successor 49214121 carries the law, which
  closes the law's end-to-end check (remote_jobs.md §0d). ds_00 runs from 21:45.
* **21:47-22:55, panels bz_00, ds_00, bs_01; ff30_gdepth_fz freeze check.** Folded fraction on 36
  domains: bz_00 0.450, bs_01 0.430 and ds_00 0.412 against their port-rate twins' 0.478, 0.465
  and 0.477; d9_00 dominates ds_00 in helix, b9_01 dominates bs_01 in helical glycine, and b9_00
  does not dominate bz_00 (findings 1.21). ff30_gdepth_fz steps 0-1 clean; its freeze check passed
  (hbond.h5 and sheet identical to ff2.1's, `checks/fz_init_20261007/dz_step00`). ff21_ctrl_fz's
  successor 49205852 started 22:28 from step 2 as designed; its next link carries the law.
* **22:47-23:59, TM4 ds_00 and panel bz_01.** ds_00 (ff30_gdepth_si step 18): unwound 9, flipped 3
  of 12, against gdepth_start's 7 and 2 (p 0.67, 1.00) and d9_00's 5 and 2 (p 0.21, 1.00);
  last-block TM4 0.775, the lowest mean of any set (against d9_00's 0.900, added Mann-Whitney
  p 0.05). s7's jump of 9732 at t 3776 heats the protein to 4.3 and the lipids to 3.0 times their
  kinetic energy; its run KE/1.5kT is 1.080, below the 1.2 flag (findings 1.24). Panel bz_01:
  folded 0.461 against b9_01's 0.465 and bz_00's 0.450, so the frozen run's gap closes at epoch 1
  (findings 1.21). Seeds and tables copied to `tm4_local/runs_cov/`. c9_00 seeds 13-24 run from
  23:41.
* **01:47-02:10 10-09, TM4 c9_00 at 24 seeds; bz_02.** c9_00's seeds 13-24 unwound 9 and flipped 8
  of 12, against seeds 1-12's 7 and 1 (same input, binary and flags; p 0.009 between the halves).
  At 24 seeds against ff21_released: 16 and 9 against 9 and 5 of 12 (p 0.72, 1.00), last-block TM4
  0.843 against 0.827. Its 12-seed lean toward fewer flips was chance, and c9_01's flip excess over
  c9_00 (one-sided p 0.034 at 12 seeds) is p 0.50 against the 24; findings 1.24 corrected in the
  c9_00, c9_01, d9_00, d9_01 and fp_e00 paragraphs and the resolution note. bz_02 (ff30_bio_fz
  step 56) extracted (frozen files identical to bz_01's), panel 49214736 submitted, patched and
  queued after d9_02 (from 01:37). New steps clean (d9 to 66 at +0.073; bz to 56; dz to 6; cz to 7).
* **02:15-03:15 10-09, TM4 secondary structure by DSSP (user).** The TM4 test scored helix from
  phi/psi boxes alone. Added DSSP (mdtraj, on N, CA, C and Upside's own carbonyl O from
  `infer_H_O`) to `tm4_local.py` and `tm4_compare.py`: alpha-helix and any-helix fractions of
  135-151 per block, alpha per residue, Mann-Whitney on each (backups `*.bak_pre_dssp_20261009`).
  Old and new scripts give identical old lines on four comparisons; all 28 tables regenerated with
  no earlier line changed (pre-DSSP copies in `tables_pre_dssp_20261009/`). DSSP alpha is 0.10-0.22
  below the box in every set: where they disagree TM4 has turned partly into pi-helix (i->i+4 O...N
  3.9-5.3 A) or frayed at 147-151. The order of sets is nearly the same (Spearman 0.88) and nothing
  resolves (findings 1.25, lesson 10.20). Uploaded scripts, README and tables to `tm4_local/`.
* **03:09-03:16 10-09, watch pass; panel d9_02.** ff30_bio_fz's link 49202419 COMPLETED 02:55
  after its 54 steps (5-58); successor 49206071 pending. Panel d9_02 COMPLETED 02:04: folded
  0.439, helix -0.030 (41 domains), dominated by gdepth_start in helix; 0.439 against d9_01's
  0.457 on 39 shared domains (findings 1.21). ff21_ctrl_fz's step 10 relaunched 12 workers (5
  twice) and all 24 ran. New steps clean (d9 to 67 at +0.080; bz to 58; dz to 9; cz to 9).
  Successor start estimates slipped: 49206580 to 17:06, 49205039 and 49205173 to 10-10 06:10.
* **03:20-03:55 10-09, DSSP primary (user); TM4 d9_02.** `tm4_compare.py` now leads with the primary
  test, last-block DSSP alpha of 135-151 by a two-sided Mann-Whitney with a bootstrap interval; the
  dihedral counts are secondary (backups `*.bak_pre_dssp_primary_20261009`). All 30 tables carry it.
  No in-training checkpoint is resolved from its start (findings 1.25 table). d9_02 (01:37-03:35):
  alpha 0.765 against gdepth_start's 0.657 (p 0.58) and d9_01's 0.615 (p 0.21), levelling off over
  the last two blocks, no energy jump. Power by resampling: at 12 seeds a true 0.10 / 0.15
  difference is detected 18% / 31% of the time, at 24 seeds 34% / 58%. bz_02 runs from 03:35.
* **04:04-04:15 10-09, watch pass; ff21_ctrl_fz link ending on midway2-0088.** Step 11 of
  ff21_ctrl_fz (from 03:42) lost 8 srun launches; 2r2y, 3jyz and 4qbo failed all three on
  midway2-0088 ("Invalid job credential"), so the step raises and link 49205852 ends; successor
  49214161 resumes from step 10. 0088's second such event (first: ff30_bio_fz 10-07 14:40), the
  documented record to exclude it; asked the user. Steps 10 of cz and dz and 68 of d9 clean (d9
  margin +0.087). Every pending job's start estimate is now 10-10 03:33-03:53.
* **04:12 10-09, midway2-0088 excluded (user).** Added to `/project/trsosnic/yinhan/slurm/midway2.args`
  with its record in README.md (backups `*.bak_pre_0088_20261009`); `--test-only` accepted the list,
  the wrapper refused `-w midway2-0088`, and `update_pending.sh` set all 9 pending jobs to it. Every
  midway2 `slurm.args` symlink reads it. ff21_ctrl_fz's link 49205852 FAILED 04:11 as expected (3 of
  24 workers failed); successor 49214161 resumes from step 10.
* **05:47 10-09, watch pass; TM4 bz_02.** Steps 70 of d9 (margin +0.091) and 12-13 of dz clean
  (24 of 24, KE/1.5kT at most 1.016); every pending start estimate 10-10 03:52. bz_02 finished
  05:33: DSSP alpha 0.856 against its twin b9_02's 0.694 (+0.162 [-0.010, +0.330], p 0.046), the
  first resolved primary test, only just; against bio_start 0.769, p 0.30. Tables
  `tm4_compare_cov_bz_02.txt`, `_vs_b9_02.txt`; findings 1.24-1.25, remote_jobs §1 and Handoff,
  plan.md. The TM4 queue is empty until dz_00.
* **08:47 10-09, watch pass; dz_00 and two links started.** ff30_bio_dt009's successor 49206580
  (08:23, `step 64 of 76`) and ff30_gdepth_si's 49205173 (08:40, `step 21 of 76`) started as
  designed, and their successors carry the law. ff30_gdepth_fz reached its epoch-0 end (step 18,
  depth round 1 dL - dR +0.567). dz_00 extracted (hbond.h5 and sheet md5-identical to gdepth_start's),
  panel 49216103 submitted by submit_new.sh, patched, and TM4 started 08:50. d9 steps 72-74 clean
  (margin +0.082); bz_02's panel runs (116 of 176).
* **09:47 10-09, watch pass; ff30_gdepth_dt009 converged, two NODE_FAILs.** ff30_gdepth_dt009
  reached step 76 (08:49); gate 49216104 CONVERGED (steps 58-76, every p > 0.005, dhb 0.045) and
  released step 75 as ff_3.0_gdepth (round-4 library; margin +0.082, under the hold); 32 benchmark
  arms and 4 glpG chains queued, the first arms' logs sane. d9_03 extracted (byte-identical to the
  release), panel 49216387, patched, queued behind dz_00. At 09:17:42 ff30_bio_dt009's and
  ff30_gdepth_si's links NODE_FAILed with other users' jobs on excluded nodes (the controller-link
  signature); midway2-0103 and 0116 NOT_RESPONDING after; successors resume from steps 64 and 21;
  asked the user about excluding the two. Panels: bz_02 0.424 (b9_02 0.448), dz_00 0.432 (d9_00
  0.477, dominating). polygly link 1 COMPLETED at -maxh (~48 ns per replica, mindist -pi >= 2.35
  nm). findings 1.21 and new 1.26; remote_jobs §1, §8 and Handoff.
* **10:30 10-09, midway2-0103 and 0116 excluded (user).** Added to `/project/trsosnic/yinhan/slurm/midway2.args`
  with a README.md row (backups `*.bak_pre_0103_20261009`); `--test-only` accepted the list, the
  wrapper refused both nodes as it refuses 0088 and accepted 0110, and `update_pending.sh` set all
  19 pending jobs to it. No running job sits on either node.
* **10:40 10-09, TM4 residues along ff30_gdepth_dt009 (user question).** New
  `scratchpad/ff3_local_test/scripts/tm4_residue_ss.py` (tm4_local.py's own reader and DSSP) tables
  the last-block DSSP code per TM4 residue; `tm4_residue_ss_d9.txt` for gdepth_start, d9_00-02 and
  ff2.1. The weak end moves (N-terminal at d9_00, C-terminal 148-151 at d9_01-02); ff2.1's
  midplane pi loss is absent. findings 1.25 "Where it is lost".
* **10:42 10-09, ff_3.0 locally (user).** Copied ff30_gdepth_dt009's release (six files) into
  `parameters/ff_3.0/` (untracked; the two Sep 10 `.bak` files there kept); md5 equal to
  `$P/parameters/ff_3.0_gdepth` and byte-equal to `ff/d9_03`. Recorded the rama/sheet/hbond paths
  a run must pass (plan.md Phase 11, findings 1.26, remote_jobs 1b).
* **10:47 10-09, watch pass; TM4 dz_00.** dz_00 finished 10:46: DSSP alpha 0.791 against
  gdepth_start's 0.657 (p 0.47) and its twin d9_00's 0.801 (p 0.80), not resolved; strong N-terminal
  half, weak C-terminal, the reverse of d9_00; s11's jump 23501 the largest on fixed inputs.
  d9_03 started 10:46. Steps 55 of bs (+0.058) and 21-22 of dz clean; all 32 benchmark arms run.
* **10:58 10-09, bs_02.** ff30_bio_si step 56 clean (margin +0.055, 24 of 24, KE/1.5kT at most
  1.014); panel 49216513 submitted by submit_new.sh with the new law; extracted (rama.dat as
  bs_01's), patched, queued behind d9_03.
* **11:52-15:50 10-09, MacBook Pro session: glpG TM1** (user). Its records were kept in separate files
  while the Mac Studio owned the shared ones, and merged here the same day (findings 10.22). A read-only status
  check at 11:52: bio runs at steps 57-64 of 76, the 32 ff_3.0_gdepth benchmark arms running, 4 glpG
  chains queued.
  - Glycine map figures, the three maps on one zero (findings 1.27). Local glpG pair, ff2.1 and
    ff3.0 seed 1, bitwise equal to the cluster's seed-1 runs (1.26); figures
    `~/Downloads/glpg_tm4_ss_ff21_vs_ff30.png` (slide-13 format, phi/psi) and
    `glpg_tm4_dssp_ff21_vs_ff30.png` (DSSP alpha), the six TM helices marked from our PDB (11).
  - TM helices across 9 sets x 12 seeds (4.6) and the TM1 diagnosis (4.7); the user chose option A,
    Upside's membrane H-bond term (plan.md Phase 13). Built and verified 13:45-13:52; 24 runs
    started 13:48. The ff3.0 twelve were stopped at t 60 for the 16:00 shutdown; the ff2.1 twelve
    finished 15:39: helicity up in every TM helix (TM1 0.79 -> 0.90, TM4 primary 0.609 -> 0.788), but
    5 seeds jump > 3000 E_up and two tear TM4's backbone (C-N 9.5 A) (4.8).
  - Transient backbone excursions found in the baseline sets and localised by a dense-frame replay
    (4.9). The VTF's stationary atom is residue 210's unused O slot; `py/martini_extract_vtf.py`
    fixed (user's commit 14a27cc3) and the VTF regenerated (3.8). The DSSP companion figure differs
    from slide 13 because of its criterion (10.20).
* **12:47 10-09, watch pass; TM4 d9_03, the released ff_3.0.** d9_03 finished 12:43: DSSP alpha
  0.726 against gdepth_start's 0.657 (+0.069 [-0.165, +0.300], p 0.64) and d9_02's 0.765 (p 1.00),
  not resolved; against ff2.1's 0.610 (added) p 0.25; 5 unwound, 4 flipped, no energy jump; weakest
  at 148-151 and 140-142, no pi-helix (`tm4_residue_ss_d9_03.txt`). 27 files and the ff2.1 table in
  runs_cov, md5 verified. bs_02 started 12:43. Steps 25-26 of dz (frozen; 2j6b returned) and 58 of
  bs (+0.050) clean. findings 1.24-1.26, remote_jobs §1 and Handoff, WATCH_STATUS.md.
* **14:47 10-09, watch pass; TM4 bs_02.** bs_02 finished 14:38: DSSP alpha 0.684 against
  bio_start's 0.769 (p 0.58) and its twin b9_02's 0.694 (p 0.93), not resolved; strong N-terminal,
  weak C-terminal, the reverse of b9_02; 5 unwound, 2 flipped, no energy jump. The TM4 queue is
  empty. Steps 28-29 of dz (frozen) and 60 of bs (+0.043) clean. findings 1.24-1.25, remote_jobs §1
  and Handoff, WATCH_STATUS.md.
* **18:47-20:47 10-09, watch passes; dz_01 extracted and TM4 dz_01.** ff30_gdepth_fz's epoch-1 end
  (step 37, 18:36): extracted (hbond.h5 and sheet = gdepth_start's, rama.dat its round-2 library,
  dL - dR +0.574), panel 49218034 by submit_new.sh (carries the law), patched, TM4 18:49-20:45:
  DSSP alpha 0.780 against gdepth_start's 0.657 (p 0.58) and its twin d9_01's 0.615 (+0.165
  [-0.056, +0.385], p 0.069), not resolved; strong at 144-149 where d9_01 is weakest; 5 unwound,
  3 flipped (all GLY143). dz steps 34-41 (frozen) and bs steps 63-67 (+0.037 to +0.042) clean.
  findings 1.24-1.25, remote_jobs §1 and Handoff, WATCH_STATUS.md.
* **21:20 10-09, session closed for an upgrade (user).** The watch cron `cad6ca59` ends with it;
  nothing local runs (TM4 queue empty after dz_01). remote_jobs "Handoff" has the restart steps
  (one catch-up pass, then CronCreate with the prompt in "The watch"); "Resume here" step 4 lists the
  current owed checks. WATCH_STATUS.md carries the restart note.
* **21:25-23:00 10-09, Mac Studio after the restart: watch, md merge, TM1 kick found.** One watch
  pass (21:27; ff30_gdepth_fz step 42 clean; 25 of 32 ff_3.0_gdepth benchmark chunks carry the
  pre-10:30 exclusion list, asked); cron `45cc9e4a`. The MacBook Pro's records merged into plan.md
  Phase 13, findings 1.26-1.27, 3.8, 4.6-4.9, 10.14 and 10.20-10.22, and its three files deleted.
  TM1 (user: resolve it locally in ff3.0): option A on ff3.0, 12 seeds, running from 21:58
  (`scratchpad/glpg_tm1/glpg_tm1_hbmem/runs`); a read-only kick probe (bitwise against the stored
  runs) found the TM1-middle kicks to be the side-chain/lipid table's under-sampling (findings 4.10);
  a converged table is being built (`kick_ff30/sc_table_600x12`).
* **22:40-23:05 10-09, ff3.1 research and the route-1 probe (user).** Literature and trainer study
  (findings 1.28): the parameter count is no obstacle (Upside already trains 31,285), the data per
  depth and the native-placement target are. User chose route 1 (per-map depths under an L2 prior,
  no new simulations; plan.md Phase 14). Probe on c9's recorded basin populations
  (`scratchpad/ff31_probe/`, findings 1.29): cross-validated width 0.3 nat, per-map detail transfers
  less well than pooled offsets, and the transferable push is a global alpha_R deepening.
* **23:05-23:50 10-09, frozen arms stopped; ff3.1 built and submitted (user).** Cancelled
  ff30_gdepth_fz's link 49206070 (after its successor 49214121), the successors 49206071 and 49214161,
  and the dz_01 panel 49218034 (user: the dt 0.009 training already keeps TM4 stable). ff3.1 trainer
  `training/ff31/` (per-map basin depths; findings 1.30), tested offline (reproduces round 4's
  round 1 exactly), deployed as `ff31_depth_all` (49218340) and `ff31_depth_data` (49218341), both
  from ff2.1 at NDRD, target 76, no gate. `check_step.py` and `submit_new.sh` extended (backups
  `.bak_pre_ff31_20261009`); watch cron `cfa1b8c7`. remote_jobs.md, WATCH_STATUS.md.
* **00:09-00:20 10-10, TM1: option A on ff3.0 finished.** TM1 last-block alpha 0.818 against d9_03's
  0.685 (p 0.073, not resolved), no jump above 3000 (findings 4.11). The probe build gained a start
  evaluation (`UPSIDE_KICK_PROBE_START`) and was rebuilt in place while the fixed-table seeds ran on
  it; all 12 survived, but a live binary is rebuilt only into a fresh path from now on. Replays of
  ff2.1 + term s1 (to t 1700) and s2 (to t 2960) started, frame for frame with the MacBook Pro's logs,
  to read the tears of findings 4.8.
* **00:25-00:35 10-10, TM1 literature (user question) and the combined test.** Literature agent
  (findings 4.12): no experiment shows GlpG's TM1 opening; the core costs 6.3-9 E_up per broken
  H-bond against the hybrid's 1.96; option A supplies about the low end, A and B must not be
  combined. Fixed table + term set started 00:32 (`run_probe_set.sh`, `runs_sc600x12_hbmem/`).
* **01:40 10-10, ff2.1 + term's tears explained.** The probe replays (bit for bit with the MacBook
  Pro's runs) show both TM4 tears start with a side-chain/lipid table spike on MET142's CA (77,407 and
  153,463 E_up/A), not the term (findings 4.8, 4.10). Every kick site so far is a MET (34, 35, 142).
* **02:07 10-10, TM1 resolved by the table fix.** The fixed-table set finished: TM1 last-block
  alpha 0.931 against d9_03's 0.685 (p 0.002), TM4 0.924 against 0.726 (p 0.10), no evaluation above
  1000 E_up/A, one C-N excursion frame against 33 (findings 4.13). Production (martini.h5, the
  builder, the cluster) waits for the user. The fixed-table + term set runs to ~04:10.
* **02:47 10-10, watch pass; ff30_bio_si converged.** Step 75 at 02:18; gate 49218437 converged
  (steps 58-76, lowest dhb p 0.0125) and released ff_3.0_bio_si (= bs_03, md5 in both trees); 32
  benchmark arms (10 R, first logs clean) and 4 glpG chains queued; bs_03 panel 49218518 (submit_new).
  bs_03 extracted, patched (matches the gate's own patch), TM4 from 02:50; the TM1 combined set's seeds
  9-12 paused for it and resumed when it finishes. remote_jobs.md §1, WATCH_STATUS.md, plan.md Phase 11.
* **04:48 10-10, TM4 bs_03 (= ff_3.0_bio_si) finished.** DSSP alpha 0.739 against bio_start's 0.769
  (p 0.89), bs_02's 0.684 (p 0.71) and d9_03's 0.726 (p 1.00, added): not resolved; unwound 4,
  flipped 2 of 12. s1 carries a transient protein-potential jump to +200,071 at t 1490 (recovered,
  TM4 0.99). Findings 1.24 table, 1.25, new 1.31; tables `tm4_compare_cov_bs_03*.txt`, runs to
  `runs_cov/` (local hardlinks, cluster rsync). The TM1 seeds 9-12 resumed at 04:46.
  remote_jobs.md §1, WATCH_STATUS.md.
* **06:05 10-10, TM1: the term on the fixed table (Phase 13 step 5 done).** The fixed-table + term
  set finished 05:49: TM1 0.964 and TM4 0.959 against d9_03's 0.685 and 0.726 (p < 0.001, 0.009),
  both resolved; against the table fix alone TM1 +0.033 (p 0.049), TM4 p 0.21; no C-N excursion, no
  jump; three evaluations above 1000 E_up/A, protein side-chain contacts (ASN185, SER119), not the
  table or the term. Findings 4.14 (new), plan.md Phase 13. TM1 data and README_mac_studio.txt
  uploaded to `checks/tm1_hbmem_20261009/mac_studio_20261010/` (5.7 GB, sizes verified 06:31).
* **08:30-08:58 10-10, the side-chain/lipid table fix in production (plan.md Phase 15; findings
  4.15).** User: apply it so future runs use it, everywhere (re-patch the queued glpG chains, rerun
  TM4 references). `py/martini_build_tables.py`: the SC task as arrays, sampling 600 x 12 by default,
  env overrides removed; at 24 x 1 it reproduces the old table (9 cells one float32 step off), at
  600 x 12 the tested table to float32 rounding; deterministic. `parameters/ff_2.1/martini.h5`
  rebuilt (`2e612261...`, old `.bak_pre_sc600x12_20261010`). `patch_glpg.py` (`97ca3294...`) copies
  martini.h5's table into patched seeds. Cluster: martini.h5 in both trees, builder in /beagle3, patch
  tool in `$P/training/` and `tm4_local/scripts/` (`checks/sc_table_20261010/deploy.sh`); the 8
  queued glpG seeds re-patched and checked (`repatch.sh`; first check script failed on PyTables'
  vlen strings, rewritten with h5py). TM4: six `_sc` reference sets queued from 08:50 (d9_03_sc,
  bs_03_sc, ff21_bioT1_6_sc, gdepth_start_sc, ff21_released_sc, d9_01_sc). up.md 2.6; watch prompt
  and cron (`aa713659`, replacing `cfa1b8c7`); remote_jobs.md, WATCH_STATUS.md.
* **10:50 10-10, watch pass; first fixed-table results.** d9_03_sc (08:50-10:45) reproduces Phase
  13's table-only set within noise (TM1 0.882 against 0.931, p 0.19; TM4 0.884 against 0.924, p 0.24;
  no jump), not frame for frame (float32-rounding table difference, probe build); table
  `glpg_tm1/tm1_compare_d9_03_sc_vs_sc600x12.txt`, runs to `runs_cov/`. ff_3.0_gdepth's 79HIS glpG
  chains started 10:16 on the fixed table (replicas 0, 13, 27 checked). ff_3.0_bio_si proteinL de
  novo NODE_FAILed 09:52 on 0151 and was requeued by Slurm. findings 4.15, plan.md Phase 15,
  remote_jobs.md §1, WATCH_STATUS.md.
* **11:20 10-10, the membrane H-bond term not adopted (user).** With the table fixed its gain is
  below the spread between two table-only sets; it imposes the bilayer geometry from outside. findings
  4.14-4.15, plan.md Phase 13/15 and Known Errors, WATCH_STATUS.md. ff30_bio_dt009's successor
  49215780 R from 10:44 with the owed first lines (remote_jobs.md §1).
* **11:43 10-10, watch pass; d9_03 panel.** d9_03 (= ff_3.0_gdepth) folds 0.464, helix -0.030 on
  41 domains, dominated by gdepth_start in helix, beta and helical glycine (findings 1.21).
  ff30_gdepth_si's link 49216034 R from 11:10 (first line as owed); ff30_bio_dt009 step 65 +0.012; all
  4 ff_3.0_gdepth glpG chains R on the fixed table. remote_jobs.md §1, WATCH_STATUS.md.
* **12:43 10-10, watch pass; bs_03_sc TM4 and bs_02 panel.** bs_03_sc against d9_03_sc 0.951 against
  0.884 (p 0.10), no jump (findings 4.15); bs_02 panel 0.429, helix -0.030, dominated by gdepth_start
  in helix, not by its twin b9_02 (findings 1.21). ff31_depth_all's first link R from 12:39 (first
  lines as owed); polygly link 2 R from 11:55; all 8 glpG validation chains R on the fixed table.
  ff30_bio_dt009 steps 66-67, ff30_gdepth_si step 22 clean. remote_jobs.md §1, WATCH_STATUS.md.
* **13:43 10-10, watch pass; bs_03 panel.** bs_03 (= ff_3.0_bio_si) folds 0.467, helix -0.028 on 41
  domains, dominated by gdepth_start in helix; with d9_03 0.467 against 0.464, helix -0.027 against
  -0.033 (findings 1.21). ff31_depth_all step 0 clean (offsets 0 at round 0); ff30_bio_dt009 steps
  68-69 and ff30_gdepth_si step 23 clean; 20 bio_si benchmark chunks rolled over. remote_jobs.md §1,
  WATCH_STATUS.md.
* **14:43 10-10, watch pass; bs_03_sc against its start.** bs_03_sc 0.951 against ff21_bioT1_6_sc's
  0.939 (p 0.77), no jump (findings 4.15); gdepth_start_sc running. ff30_bio_dt009 steps 70-71,
  ff30_gdepth_si step 24 and ff31_depth_all step 1 clean; ff_3.0_gdepth BBL native benchmark at
  target. findings 4.15, plan.md Phase 15, remote_jobs.md §1, WATCH_STATUS.md.
* **16:43 10-10, watch pass; ff30_bio_dt009 gate CONVERGED, d9_03_sc against its start.**
  ff30_bio_dt009 reached step 76 at 16:37; its gate passed every group (lowest dhb p 0.13) and
  released b9_03 as ff_3.0_bio, md5 verified; 32 benchmark arms and 4 glpG chains queued with the
  law (findings 1.32). b9_03 extracted, its panel submitted, patched and queued as b9_03_sc.
  d9_03_sc against gdepth_start_sc 0.884 against 0.947, p 0.13 (findings 4.15). findings 1.32 and
  4.15, plan.md Phase 11, remote_jobs.md §1 and "The TM4 queue", WATCH_STATUS.md.
* **18:43 10-10, watch pass; ff2.1 on the fixed table.** ff21_released_sc holds TM4: 0.944 against
  its old-table 0.610 (p 0.002), unwound 0 against 9 of 12; every start reads alike and neither
  tested candidate improves on ff2.1 (findings 4.15). Tables `tm4_compare_cov_ff21_released_sc*.txt`,
  `_gdepth_start_sc.txt`, `_d9_03_sc_vs_ff21_released_sc.txt`, `_bs_03_sc_vs_ff21_released_sc.txt`.
  ff_3.0_bio's 32 benchmark arms R with clean first logs. findings 4.15, plan.md Phase 11,
  remote_jobs.md §1, WATCH_STATUS.md.
* **20:43 10-10, watch pass; d9_01_sc.** d9_01_sc against gdepth_start_sc 0.925 against 0.947,
  p 0.45, not resolved, 0 unwound or flipped (findings 4.15); tables `tm4_compare_cov_d9_01_sc.txt`
  and `_d9_03_sc_vs_d9_01_sc.txt`, set and tables copied to runs_cov, md5 verified. b9_03_sc
  running. ff30_gdepth_si step 30 and ff31_depth_all step 7 clean. findings 4.15, remote_jobs.md §1
  and "The TM4 queue", WATCH_STATUS.md.
* **21:38 10-10, the ff3.1 steric check (user, after Jumper).** New `training/ff31/steric_check.py`
  (deployed `checks/ff31_steric/`, md5 e35e3d98): each step's library against a hard-sphere
  Ramachandran map, side-chain and backbone clash classes, per coil map and per training residue's
  engine map, plus the non-glycine phi > 0 populations of the free ensembles. Verified: CB within
  0.08 A of 1AO6's, alanine's allowed share 21%, round 0 x1.0000; on the predicted first round no map
  above x1.17. ff31_depth_all step 8 (+0.167): no flag. Watch step 3 and the cron prompt carry it
  (cron `c6b0479e` replaces `aa713659`). findings 1.33, plan.md Phase 14 step 5, remote_jobs.md,
  WATCH_STATUS.md.
* **22:46 10-10, watch pass; b9_03_sc, and ff30_gdepth_si's node failure.** b9_03_sc (= ff_3.0_bio)
  0.950 against its start ff21_bioT1_6_sc's 0.939 (p 0.36) and its twin bs_03_sc's 0.951 (p 0.30),
  not resolved; against d9_03_sc 0.884, p just under 0.05 with the bootstrap touching zero (added);
  set and 4 tables copied to runs_cov, md5 verified (findings 4.15). ff30_gdepth_si's link 49216034
  NODE_FAILed 20:52:43 in step 32 with another user's job on midway2-0151; step 31 (+0.127) is
  complete and successor 49219651 is PD. ff31_depth_all steps 9-10 clean, steric check no flag; the
  watch's steric command now names its step (without it the script reads the newest step).
  findings 4.15, remote_jobs.md §1 and watch step 3, WATCH_STATUS.md.
* **10-10, literature: TM-helix HDX in a bilayer.** Three literature agents plus the lab's H145A HXMS
  file. TM cores exchange slowly but finitely (PF ~10^5-10^7); exchange needs both an opening and
  water/OH- access, and the bilayer suppresses both. Our binary `acc` makes lipid-shielded PF infinite by
  construction. Corrected the TM4 140-144 uptake quoted in findings 5.1 and 7.1 (34% at 24 h, not 50%).
  findings 5.5.
* **10-11, HDX follow-up (records only, no code changed).** Recorded the uptake PF = 1000 cap, the
  "non-exchanging" wording and the interface-boundary check as pending the user's decision (findings 5.5,
  plan.md Known Errors). TM4: crystal SASA and a tail-sized probe show TM4 protein-buried with ~5 amides
  within 7 A of lipid-reachable surface, matching the hybrid's shielded set (open item closed in
  remote_jobs.md); the HXMS TM4 peptide gives dG ~6.1-7.4 kcal/mol, so the implicit's off-scale TM4 is
  consistent with experiment and the hybrid's 3.5-4.9 is too labile. Corrected findings 5.3b's model ranking.
* **23:46 10-10, watch pass; quiet.** ff31_depth_all step 11 (+0.158) clean, steric check no flag;
  ff_3.0_bio_si's second benchmark chunks resubmitted at their wall. **00:00 10-11:** session
  upside2-md-87 asked for a glpG HDX analysis on the fixed-table REMD chains (reported user-approved,
  low priority); recorded as deferred to block 2 (~10-13) in remote_jobs.md §1. The WATCH_STATUS.md
  upload stalled at ~00:40 with `/project` unresponsive on midway2-login1 (`echo` and `hostname`
  answer); the local copy (md5 c33ad9e6) is the one to upload. If the stalled write lands, the
  remote file reads e3ee8d8f (same text, one path missing). remote_jobs.md §1, WATCH_STATUS.md.
* **00:46 10-11, watch pass.** ff31_depth_data started 00:06 (24 nodes; first lines as owed,
  successor 49221772). Reads of the cluster's WATCH_STATUS.md still hang on midway2-login1 (a GPFS
  wait on that one file; the rest of `/project` answered) until ~00:56; it then cleared and the file
  was rewritten at 01:00, md5 verified. remote_jobs.md §1.
* **01:46-04:46 10-11, watch passes.** ff31_depth_data's first steps 0-4 and ff31_depth_all's
  12-15 clean, steric check no flag (largest residue x1.046). ff30_gdepth_si's successor resumed
  step 32 at 03:27 (30 min steps). ff_3.0_bio's 4 glpG chains started 03:27, their replicas on
  martini.h5's table. **b9_03 panel** (= ff_3.0_bio): 0.452, helix -0.029, dominated by
  gdepth_start in helix; of the three released candidates bs_03 leads (0.467) and dominates d9_03
  in helix and beta (findings 1.21). remote_jobs.md §1, WATCH_STATUS.md.
* **05:46-06:46 10-11, watch passes; two epoch ends.** ff31_depth_all step 18: depth round 1 as
  predicted (held-out gap left 0.990, 66 depths beyond 0.1, findings 1.30); the steric check on the
  first moved library has no flag (coil maps at most x1.10, findings 1.33). ff30_gdepth_si step 37
  (+0.111), depth round 2. ds_01 and fa_00 extracted (rama.dat md5-equal to the runs' libraries),
  panels submitted, glpG inputs patched; their local TM4 sets held at the user's word (06:48, "run
  another session's simulation first"), and upside2-md-87 reports nothing to run locally. findings
  1.30, 1.33, remote_jobs.md §1 and "The TM4 queue", WATCH_STATUS.md.
* **07:09-09:20 10-11, local TM4 HDX on the fixed table (plan.md Phase 16; user: "run this
  session's job first").** ff21_released_sc, 12 seeds, T 0.85, 4000 tu, frame interval 2, on the
  Mac Studio (`scratchpad/tm4_hdx_local/`), ahead of the TM4 queue (released 09:05). Analysis on the
  stale set's ff2.1 HDX topology (md5 166c9d12; H-bond parameters equal to the input's) with the
  repo's HDX steps, through a scratch projection because the local runs log at `basic` (no
  `/output/hbond`). Result: TM4 137-148 lose protein protection in at most 0.06% of frames at
  298 K, against 10-25% in the surviving pre-fix 48-rung set; 142-144 never exchange (dG > 5.8);
  DSSP alpha 0.94-0.95, none unwound. The 10-10 T 0.80 sets agree. Findings 5.5 rewritten,
  including the stale set's provenance (files dated 08-27, not 09-04). A stray `results/` written
  to the repo root by the first test was removed. 17:20: the user chose a remote REMD over a local
  old-table control (the Mac Studio is busy); an ff2.1 79HIS fixed-table chain with HDX after block
  2 was requested of upside2-md-98 (plan.md Phase 16). Held on a 98%-full `/project` until the
  user approved a 928 G reclaim; submitted 17:41 for 2 blocks by upside2-md-98, seed checks in
  `checks/ff21sc_remd_20261011/`.
* **11:46 10-11, watch pass; ds_01_sc and two panels.** ds_01_sc 0.894 against gdepth_start_sc
  0.947 (p 0.19) and d9_01_sc 0.925 (p 0.62), not resolved; set and tables in runs_cov, md5
  verified (findings 4.15). Panels ds_01 0.420 and fa_00 0.471 (ff3.1's first end), each dominated
  by gdepth_start in helix (findings 1.21). ff30_gdepth_si step 47 margin +0.100008, on the hold
  line. fa_00_sc running. findings 1.21, 4.15, remote_jobs.md §1, WATCH_STATUS.md.
* **12:46 10-11, watch pass; ff30_gdepth_si on hold.** Its H-bond margin fell below +0.10 at step 48
  (+0.099957; step 49 +0.099370), read at full precision from the extracted hbond.h5; the user was
  told and the run trains on. ff31_depth_all 25 and ff31_depth_data 13-14 clean, no steric flag.
  remote_jobs.md §1, WATCH_STATUS.md.
* **13:46 10-11, watch pass; fa_00_sc.** ff3.1's first end in glpG TM4: 0.926 against
  ff21_released_sc 0.944 (p 0.60), not resolved, 0 unwound or flipped; lower at 135-145, higher at
  149-151 (findings 4.15); set and table in runs_cov, md5 verified. Local queue empty. Training
  steps clean, no steric flag. findings 4.15, remote_jobs.md §1, WATCH_STATUS.md.
* **15:46 10-11, watch pass; ff31_depth_data's first epoch end.** Depth round 1 as predicted
  (held-out gap left 0.994, 3 of 1,442 depths beyond 0.1; findings 1.30); the steric check on its
  library has no flag (findings 1.33). fd_00 extracted (rama.dat = round 1's), panel 49223540
  running, patched, TM4 fd_00_sc running from 15:48. findings 1.30, 1.33, remote_jobs.md §1,
  WATCH_STATUS.md.
* **17:15-17:45 10-11: disk and a new REMD chain.** Checking headroom for session upside2-md-87's
  request found `/project/trsosnic` at 98% (88 G free), with the 12 glpG chains needing ~410 G more.
  At the user's word 928 G was freed (panel `.up` trajectories 397 G, NP-1AO6 run data 489 G, the
  four ff_3.0_cancelled glpG directories 40 G; manifests in `checks/disk_reclaim_20261011/`); the NP
  and panel scripts were kept, two 1.5 KB NP analysis scripts inside a deleted directory were lost.
  `panel.py select` reproduces its table from the npz alone. Then, user-approved for 2 blocks, glpG
  REMD 49223635 (79HIS, ff2.1 on the fixed table) from a seed checked against the pristine ff2.1
  seed, the live seed and the Mac's input. The watch gains a disk step (2b). remote_jobs.md §1,
  "The watch", WATCH_STATUS.md.
* **17:46 10-11, watch pass; fd_00_sc and two panels.** fd_00_sc 0.914 against ff21_released_sc
  0.944 (p 0.23) and fa_00_sc 0.926 (p 0.53), not resolved, one seed with GLY149 flipped; the two
  ff3.1 variants fray TM4 at opposite ends (findings 4.15). Panels fd_00 0.479 and ds_02 0.406, both
  dominated by gdepth_start in helix (findings 1.21). ds_02_sc running; d9_02_sc, its fixed-table
  twin, patched and queued. The ff2.1_sc REMD chain calibrated at 639.6 ms/step.
* **18:46-19:46 10-11, watch passes; ds_02_sc.** ds_02_sc 0.889 against gdepth_start_sc 0.947
  (p 0.47), not resolved, weakest at 135-139 (findings 4.15). ff_3.0_gdepth's 79HIS glpG chain
  finished block 1 whole and started block 2. ff31_depth_all margin +0.103 at step 33, near the
  hold line; training steps clean, no steric flag. remote_jobs.md §1, WATCH_STATUS.md.
* **21:46 10-11, watch pass; d9_02_sc.** d9_02_sc 0.916 against gdepth_start_sc 0.947 (p 0.17) and
  ds_02_sc against it 0.889 (p 0.67), not resolved; across ff30_gdepth_dt009's epoch ends TM4's
  C-terminal end frays step by step, the SI-rate run's does not (findings 4.15). Local queue empty.
  ff_3.0_gdepth's 79HIS_S115T glpG chain in block 2. TM4 queue paragraph in remote_jobs.md rewritten
  as one block. findings 4.15, remote_jobs.md §1, WATCH_STATUS.md.
* **22:46 10-11, watch pass; polygly link 2.** Link 2 ended at its -maxh at 94.4-94.7 ns per replica,
  32.4-32.6 ns/day; mindist -pi over the appended trajectories 2.31-3.28 nm (gmx 2024.4 needs
  `module load gcc/10.1.0` on the login node, now in the §1 row). Link 3 runs. All four ff_3.0_gdepth
  glpG chains and ff_3.0_bio_si's 79ALA_S115T are in block 2. Training steps clean, no steric flag.
  remote_jobs.md §1, WATCH_STATUS.md.
* **23:05 10-11, polygly stopped (user).** The user judged poly-Gly no longer needed, since the table
  fix repaired TM4's glycine flips (findings 4.15); agreed. STOP written, successor 49223804 and link 3
  49219705 cancelled, ~97 ns per replica checkpointed and kept (`polygly/` is only 213 MB, so nothing
  was deleted). plan.md Phase 12 closed; remote_jobs.md §1 row moved to "Finished runs"; the watch
  prompt text no longer names polygly.
* **23:10-23:25 10-11, the local TM4 glycine-flip test retired (user).** "We have a working ff and we
  are trying to expand it so that metric no longer matters"; then, 23:15, "I meant the test for tm4 GLY
  flips, not for tm4 hdx". Nothing was running or queued locally. Watch steps 4-5 (extract, patch,
  queue, compare) replaced by a retirement note; the watch is steps 1-5 now. TM4 HDX continues (the
  ff2.1_sc REMD chain, the HDX requests). Watch cron `eae86453` (23:20, expires 10-18) replaces
  `c6b0479e` with the prompt in remote_jobs.md. plan.md Phase 14 step 4, remote_jobs.md "The watch",
  "The TM4 queue" and §1, WATCH_STATUS.md.
* **23:30-23:40 10-11, frozen-H-bond arms deleted (user).** Re-checked: no job runs from ff30_bio_fz,
  ff30_gdepth_fz or ff21_ctrl_fz (cancelled 10-09 23:09). Their step directories (6.1 G) went first;
  then, at the user's "they are wrong and no longer needed", the whole run directories and their three
  first-link logs (manifests in `checks/disk_reclaim_20261011/`). Their derived bz/dz products are kept
  pending the user. Lesson 10.25. remote_jobs.md, findings 10.25, WATCH_STATUS.md.
* **23:50 10-11, the BioEmu-map line dropped (user).** "We decided to go with basin depth training
  path (as it is expandable to all residues)." The three BioEmu runs had finished; their 57 running
  validation jobs (64 benchmark arms' chunks, 8 glpG chains) were cancelled after STOP files, none
  resubmitted. Deleted with manifests: step directories (39 G, mostly ff30_bio's half-written step),
  the 8 glpG chain directories (92 G), the benchmark replica files (40 G on /beagle3); setup, released
  parameters, seeds, logs, extracted ends and panel npz kept. /project 1,144 G free. plan.md decision
  "The basin-depth path", findings 1.32, remote_jobs.md ("Resume here" rewritten, §1 rows moved to
  "Finished runs"), memory, WATCH_STATUS.md.
* **23:45 10-11, the dropped lines deleted with their inputs (user).** "Delete them as well. Also no
  need to leave input file for this case." Cluster: the BioEmu run directories and logs, released
  ff_3.0_bio and ff_3.0_bio_si (both trees), the 8 glpG chains' seeds and logs, the 64 benchmark run
  directories and logs, and both lines' extracted ends, panels and TM4 sets (851 targets, 12,133
  files, 48 G apparent; manifest `checks/disk_reclaim_20261011/manifest_bio_and_derived_20261012.txt`).
  Mac: 871 local files (45 G). Kept: comparison tables, `select_last_*` tables, the BioEmu library
  itself (repo `rama31.dat`, `checks/gly_bioemu_map/`). Panel selection checked without bio_start;
  `submit_new.sh` lists only the live lines (it then submitted fa_01's panel). remote_jobs.md's
  current-state header and "Handoff" rewritten without the retired lines. plan.md, remote_jobs.md,
  WATCH_STATUS.md.
* **23:50 10-11, watch pass; ff31_depth_all's epoch-1 end.** Depth round 2: 235 of 5,040 depths beyond
  0.1, held-out gap left 0.994 (findings 1.30); steric check on round 2's library no flag (findings
  1.33); fa_01 panel running. ff31_depth_data step 29 +0.110, ff30_gdepth_si step 70 +0.1004. Watch
  cron `56a4d6f9` replaces `eae86453` with a prompt for the basin-depth path only. Tonight's deletion
  times corrected to 23:45 10-11.
* **00:25 10-12, ff30_gdepth_si stopped and deleted (user).** "Since ff3.0 is working, it is the only
  branch we need (and its expansion to all residues)." Successor 49222183 and link 49219651 were
  cancelled at step 72, before the target and gate, so nothing was released. Deleted with its inputs
  (cluster 54.5 G, mostly a half-written step; Mac 8 G): run directory, first-link log, ds_00-02
  extracted ends, panels and TM4 sets. `submit_new.sh` lists d9, c9, fa, fd only. Watch cron
  `c618e779` replaces `56a4d6f9`. plan.md, findings 1.32, remote_jobs.md, memory, WATCH_STATUS.md.
* **01:51 10-12, watch pass; panel fa_01.** fa_01 (ff31_depth_all step 37, round 2's library) folded
  0.444, helix -0.033 on 39 domains, dominated in helix by gdepth_start and by fa_00 (findings 1.21).
  ff31_depth_all link 2 49219752 resumed at step 39 after link 1's wall (00:39), step 39 +0.1052;
  ff31_depth_data step 31 +0.1068; steric no flag on all four new steps; KE/1.5kT at most 1.012.
  WATCH_STATUS.md rewritten short; stale frozen-arm and ds-set lines in remote_jobs.md §1 corrected.
* **06:46 10-12, watch pass; ff31_depth_data's epoch-1 end.** Depth round 2 (step 37): 26 of 1,442
  depths beyond 0.1, step rms 0.0185, held-out gap left 0.9959 (findings 1.30); steric check on round
  2's library no flag (findings 1.33); fd_01 panel 49224221 submitted. ff31_depth_all step 46 +0.108.
* **07:46 10-12, watch pass; the first steric flag.** ff31_depth_data step 39: free non-glycine beta' +
  pPII' 0.0050 above steps 0-18's 0.0047 (restrained 0.0045); the same minibatch read 0.0047 at step 1;
  clash classes unchanged. User notified; nothing changed (findings 1.33).
* **09:46 10-12, watch pass; panel fd_01.** fd_01 (ff31_depth_data step 37, round 2's library) folded
  0.472, helix -0.028 on 43 domains; gdepth_start dominates it in helix, fd_00 does not (findings 1.21).
  Steric: no new flag (fd steps 40-42 at most 0.0035). ff31_depth_all step 50 +0.1032.
* **12:46 10-12, watch pass; ff31_depth_all on margin hold.** Step 53 margin +0.0982, below +0.10
  (step 52 +0.1008); user notified, the run trains on. ff31_depth_data link 1 TIMEOUT 12:07 at its wall,
  link 2 49221772 resumed at step 45 (+0.1118).
* **14:46 10-12, watch pass; ff31_depth_all's epoch-2 end.** Depth round 3 (step 56): 345 of 5,040
  depths beyond 0.1, held-out gap left 0.9992 (transfer 1.0%, 0.6%, 0.1% over rounds 1-3; findings
  1.30); steric on round 3's library no flag, largest map x1.22 (findings 1.33); margin +0.0950, on
  hold since step 53; fa_02 panel submitted. ff31_depth_data step 49 +0.1091.
* **15:46 10-12, watch pass; panel fa_02.** fa_02 (ff31_depth_all step 56, round 3's library) folded
  0.404, helix -0.028 on 40 domains; fa_00-02 fold 0.471, 0.444, 0.404 (findings 1.21). Steric no new
  flag; ff31_depth_all step 57 +0.0979 (hold), ff31_depth_data step 50 +0.1101.
* **18:46 10-12, watch pass; ff31_depth_data's epoch-2 end.** Depth round 3 (step 56): 39 of 1,442
  depths beyond 0.1, held-out gap left 0.9983 (transfer 0.6%, 0.4%, 0.2%; findings 1.30); steric on
  round 3's library no new flag (findings 1.33); fd_02 panel running. ff31_depth_all step 62 +0.0986.
  New in the queue: the MacBook session's ff21_nobp_dt009 and the ff_2.1_ctrl38 benchmark (34 jobs,
  user's request, recorded in the run's README); a §1 row points to it; this watch leaves them alone.
* **19:46 10-12, watch pass; panel fd_02.** fd_02 (ff31_depth_data step 56, round 3) folded 0.449,
  helix -0.029 on 40 domains; fd_00-02 0.479, 0.472, 0.449, fa_02 0.404 on 39 shared (findings 1.21).
  Step 58 (step 39's minibatch, round 3) reads 0.0047, at the steric threshold: step 39 was scatter
  (findings 1.33). ff31_depth_all step 64 +0.0969, ff31_depth_data step 59 +0.1174.
* **02:46 10-13, watch pass; glpG REMD ff2.1_sc block 2.** Block 1 49223635 COMPLETED whole 02:40;
  block 2 49299095 R (`block_count` 2); STOP owed at block 2's last chunk (~12:30 10-14). ff31_depth_all
  step 73 +0.0754 (hold), ff31_depth_data step 73 +0.1430; steric no new flag.
* **03:46 10-13, watch pass; both ff3.1 runs at target 76.** ff31_depth_data 03:37 (margin +0.141,
  round 4: 58 depths beyond 0.1, gap left 0.9994), ff31_depth_all 03:53 (margin +0.070, round 4: 418
  beyond 0.1, gap left 1.0001); both stopped for the user's review, no gate (findings 1.30); steric on
  both round-4 libraries no flag (findings 1.33); last panels fd_03 49311326 and fa_03 49313290
  submitted. NTL9 de novo benchmark COMPLETE. WATCH_STATUS.md rewritten for the review.
* **05:50 10-13, watch pass; ff3.1's last panels.** fa_03 (ff31_depth_all step 75) folded 0.406, helix
  -0.030; fd_03 (ff31_depth_data step 75) 0.452, -0.025; on 37 shared domains fd_03 leads; gdepth_start
  (0.603, -0.018) dominates both in helix (findings 1.21, plan.md Phase 14 step 4). 79ALA_S115T's glpG
  chain finished block 2 whole 05:37, block 3 R. Logs clean.
* **10:20-11:00 10-13, user: no HDX on ff3.0, midway2-0151 excluded, ff3.1 reviewed.** HDX on
  ff_3.0_gdepth's chains declined (remote_jobs.md §1 spec kept for the ff2.1_sc chain). midway2-0151
  added to `slurm/midway2.args` with its README row (two NODE_FAILs 10-10; backups `*.bak_pre_0151_20261013`;
  the wrapper refuses `-w midway2-0151`; no job pending). ff3.1 review (findings 1.34, plan.md Phase 14
  step 4): matched-step panels and paired bootstraps (`checks/ff31_review/`): `all` worse than ff3.0,
  `data` ties it; two trainings of the ff2.1 control differ by 0.14 folded, so single-run gaps of
  ~0.05 are not attributable.
* **12:15 10-13, user: validate ff3.1 on the Peng set and glpG TM4.** ff31_depth_data's step 75 released
  as `ff_3.1_depth_data` by `validate_ff.sh` from launcher 49350700 (`ff31_depth_data/validate_fd03.sbatch`):
  both trees md5 verified and equal to the panel's fd_03; 32 Peng arms 49350797-49350828 and 4 glpG chains
  49350829-49350832, all with the law including 0151. Checked: patch round trip 3.6e-15; each glpG seed
  differs from its ff_3.0_gdepth twin only in the five force-field datasets (fixed table carried); a
  benchmark input's hbond equals the release; first frames sane. plan.md Phase 14 step 6, remote_jobs.md §1.
* **16:00-16:40 10-13, user: answer the PI on H145A TM5/TM6 from ff3.0 (findings 4.16).** Read the
  ff_3.0_gdepth glpG chains' rungs 0 (T 0.70) and 14 (T 0.80) and the ff2.1_sc chain's (backbone npz
  on midway2 `checks/tm5_tm6_contacts_20261013/`, login-node read, nothing submitted). Contacts, packing
  and DSSP with the established definitions; reference from the seeds, after a replica file's
  `/input/pos` proved to be the current chunk start (findings 6.2). H145A TM5 stays helical; its start
  contacts with TM6 go within ~300 tu at both rungs, TM5 packs on TM2 instead; T 0.70 and T 0.80 maps
  differ within block-to-block noise. Figures, one per plot, `~/Downloads/glpG_H145A_TM5_TM6/`;
  scripts `scratchpad/tm5_tm6/`. Files modified: findings.md (4.16, 6.2), progress.md.
* **16:45-17:20 10-13, user: the same analysis on the implicit bilayer (findings 4.16, 4.1).** No ff3.0
  or H145A implicit run exists; read the four ff2.1 WT implicit-membrane REMD runs
  (`implicit_79HIS_src{,2,3,4}`, rungs 0 and 24) through the same `analyze.py` (implicit loader added;
  hybrid values unchanged by the refactor; O placement checked equal to the engine's). TM5 fully
  helical; TM5-TM6 Q 0.20 / 0.16, TM5 onto TM2 as in the hybrid; T 0.70 vs 0.80 maps within run-to-run
  noise; the core drifts 4.8-6.1 A, so the 4.1 deviation does not need the dry-MARTINI coupling.
  Figures one per plot: `~/Downloads/glpG_WT_implicit_ff2.1_TM5_TM6/` (hybrid set regenerated by the
  same script). Files modified: findings.md (4.1, 4.16), progress.md.
* **20:25-21:00 10-13, glpG TM4 HDX: switch to ff_3.0_gdepth requested, then withdrawn (plan.md
  Phase 16).** Read the REMD driver logs of ff2.1_sc and ff_3.0_gdepth 79HIS/79ALA, blocks 1-3: on the
  fixed table ff2.1 keeps more backbone H-bonds than ff_3.0_gdepth (T 0.80-0.87, block 2: 181-186
  against 158-172), which with 4.15 contradicts the reason given for the switch; the user was told
  after the choice (lesson 10.26). The relayed "ff3.0 only" (stop ff2.1_sc, HDX on ff_3.0_gdepth block
  2) was not applied: upside2-md-98 asked the user directly, who kept the current plan. ff2.1_sc runs
  to block 2's end (~13:05 10-14) and gets its HDX; no ff3.0 HDX. Files: plan.md, findings.md (5.5
  H-bond table, 10.26), progress.md.
* **01:00-01:20 10-14, user: "you can cancel ff2.1 job. also check if TM4 secondary structure is stable
  (helical) under ff3.1".** `scancel 49299095` after re-reading the queue (Upside traps SIGTERM, KillWait
  30 s); all 28 ff2.1_sc replica files read back (30 rotated groups, finite), 18 G kept, HDX not run.
  TM4 under ff_3.1_depth_data: new read-only `checks/tm4_ff31_20261014/tm4_ff31.py` (login node, 21 s,
  rotated groups only, locking off); TM4 134-151 0.97-1.00 per group at T 0.70 over block 1's first
  6-7 chunks, as ff_3.0_gdepth on the same chunks, glycines never at phi > 0. Files modified: findings.md
  (1.35), plan.md (Phase 14 step 6, Phase 16), remote_jobs.md (§1, Finished runs, checklist), progress.md.
* **01:16 10-14, glpG TM4 HDX: Phase 16 closed on the local readout.** The user cancelled the ff2.1_sc
  REMD chain directly with upside2-md-98 (`scancel 49299095` 01:11:22, block 2 at chunk step 6,930 of
  10,494; sacct CANCELLED, data kept); its HDX and the pre-HDX checks were not run. With the ff3.0 HDX
  declined, no REMD test of TM4 remains. Files: plan.md (Phase 16 header and items, summary lines),
  findings.md (5.5 limits), progress.md.
* **01:20-01:40 10-14, glpG TM4 HDX: ff3.0 readout requested, phase handed to upside2-md-98.** Data
  check (driver logs): ff_3.0_gdepth block 2 whole in all four chains, ~5,100-5,400 frames per
  replica (300 per chunk), so a 298 K censored bound of ~5.1 kcal/mol; ff_3.1_depth_data in block 1
  (equilibration) until ~10-16. Local readouts compared at T 0.80: ff3.0 135, 137-147 censored, 136
  and 148-151 finite (148-149 mostly one seed); ff2.1 135-138 finite, 139-149 censored; implicit ff2.1
  never opens. Told the user TM4 HDX is not resolved until the REMD readout. Relayed the ff3.0 HDX
  request (msg 974b0511) and, at the user's request, handed Phase 16 and the deferred fixes to
  upside2-md-98 (msg 404b8295). Files: plan.md (Phase 16), progress.md.
* **01:20-02:00 10-14, TM4 HDX on ff_3.0_gdepth block 2 (user, relayed 01:19 and confirmed here: "Yes, 79ALA then
  79HIS"; upside2-md-87 handed over its Phase 16 work and closed).** Pre-checks (`checks/hdx_ff30_20261014/prechecks.py`):
  KE/1.5kT 0.998-1.002 per slot, swap acceptance 0.36-0.53, no ladder round trip in block 2 (walkers stay near the
  temperature their energy fits; reseed continuity and swap-set alternation verified). New launcher `hdx_block2.sbatch`
  (hdx_cluster.sbatch's steps; block 2 copied out of the live file by `extract_block.py`, tested on one group; ff_3.0_gdepth
  HDX topology; the repo's `martini_hdx_project.py`, as the cluster copy had drifted; dG on PS and PS_protein; `readout.py`).
  The user asked for midway2; no midway2 environment has mdtraj, pymbar and matplotlib, so 60668572-60668573 run on
  caslake from 01:59. Files modified: plan.md (Phase 16), remote_jobs.md (§1), progress.md.
* **02:04-02:15 10-14: TM4 HDX readout, ff_3.0_gdepth block 2.** Both jobs COMPLETED in 5 min; every replica
  joined with all block-2 frames (5,723 / 5,421, none dropped). TM4's core censored at every rung (above 5.85
  kcal/mol at 298 K; 79ALA's 137 and 139 at the bound), 142-144 against the pre-fix 4.1-4.4; only 130-134 and
  149-151 resolve. Figures (SciencePlots, censored amides at the bound) and readouts in
  `~/Downloads/glpG_TM4_HDX_ff3.0_gdepth_block2/`. Files modified: findings.md (5.5), plan.md (Phase 16),
  remote_jobs.md (§1, Finished runs), progress.md.
* **02:10-02:20 10-14, user: "Do the analysis on ff3.1 as well".** ff_3.1_depth_data is in block 1 (equilibration)
  until ~23:30 10-14, so the readout waits for its block 2 (~late morning 10-16); launcher and pre-check now take
  the force field and block (pre-check verified identical on ff3.0's 79ALA). Owed in remote_jobs.md and
  WATCH_STATUS. Files modified: plan.md (Phase 16), remote_jobs.md, progress.md.
* **02:46-15:52 10-14, hourly watch passes (Mac Studio), every log clean.** ff_3.1_depth_data's 32 arms ran
  their first two 10.5 h chunks and resubmitted themselves (third chunks 49400794-49400941); BBA native COMPLETE
  11:17 at its target (KE/1.5kT 0.999-1.007). ff_3.0_gdepth cspA de novo COMPLETE 06:32 (KE/1.5kT 1.001-1.006),
  ubiquitin de novo in its next chunk 49401175; its 79HIS and 79ALA_S115T chains finished block 3 whole 15:35 and
  15:44 and run block 4. Files modified: remote_jobs.md (§1 snapshot and rows).
* **16:00-16:10 10-14, handoff (user: "record everything in md file so another session can take over with
  \"continue jobs\", and it may be on another computer. I want to close this session").** Owner set to none in
  remote_jobs.md and the cluster's WATCH_STATUS.md; the next "continue jobs" session takes the jobs without waiting
  for "take over". Tools moved off this session's scratchpad to the cluster: `~/watch_probe/logs.sh` (log scan),
  `~/watch_probe/final_frames.py` (a finished arm's final frames), `checks/hdx_ff30_20261014/plot_tm4.py`. The ff3.1
  HDX recipe written out in remote_jobs.md §1, with the block-to-group mapping measured on ff_3.0_gdepth 79HIS
  (output_previous_0 the seed's; blocks 1-3 = groups 1-17, 18-35, 36-54) and the exact sbatch line; plan.md's
  claim that `extract_block.py` checks the calibration group corrected. Watch cron `c618e779` deleted.
  Files modified: remote_jobs.md, plan.md (Phase 16), progress.md.
* **13:47-14:40 10-15, status check (MacBook Pro; user: "check job status").** Read only: 37 jobs R, every log
  clean; 13 more ff_3.1_depth_data arms COMPLETE and healthy (14 of 32); its chains in block 2, ff_3.0_gdepth's in
  block 4; ff_2.1_nobp's arms all complete, ff_2.1_ctrl38's 23 of 32. TM4 at ff3.1's block-1 end
  (`tm4_ff31_20261015.txt`) and a new all-rung glycine scan (`~/watch_probe/gly_flip_rungs.py`): TM4 0.93-1.00 per
  group; GLY136/149 flip at some rungs under both force fields, 79HIS's GLY136 under ff3.1 at T 0.70. No watch
  started; WATCH_STATUS.md left at "Owner: none". Files modified: findings.md (1.35), plan.md (Phase 14 step 6),
  remote_jobs.md (header, owed checks, §1), progress.md.
* **14:00-14:30 10-15, scratchpad cleanup (user: "do a cleanup on the scratchpads, both in project root and
  examples").** By the findings 10.25 rule: 1,054 files deleted (224 GB; 219 GB of `.up` trajectories, the rest
  bulk npy/npz, build objects, logs and stale loose files), every script, input, README and small result kept,
  `polygly/` kept whole (resumable); `example/16.MARTINI/scratchpad/` (CLC-ec1 prep outputs) removed. Manifest
  `scratchpad/CLEANUP_20261015_manifest.txt`. Scratchpad 210 G -> 1.2 G. Files modified: remote_jobs.md
  ("Where things are"), progress.md.
* **14:05-14:20 10-15, ff2.1 deterministic-side-chain scoring (user: "submit the score now and keep record so
  another computer do not need to resubmit finished jobs").** New `setup/score/submit_unscored.sh` in
  `$P/training/ff21_nobp_dt009`: submits `score.sbatch` only for arms with `COMPLETE`, no json and no queued score
  job, logging to `submitted.log`. It submitted the 19 unscored complete arms (49406327-49406359, broadwl, the
  law); a rerun at once submitted nothing. Records on the cluster (the run's README "Scoring", `setup/notes/plan.md`
  P6b, `WATCH_STATUS.md` "Not watched here" and "Owed checks", each backed up `.bak_pre_score_20261015`) and in
  remote_jobs.md (nobp row, snapshot, owed checks, the HDX recipe's block-2 ids). Files modified: remote_jobs.md,
  progress.md.
* **14:30-14:45 10-15, ff3.1 glycine Ramachandran figures (user).** ff_3.1_depth_data's rama.dat copied from the
  cluster to `scratchpad/ff_3.1_depth_data/` (md5 4b6a9aaf, equal in both trees and to `rama_round_04.dat`). New
  `scratchpad/plot_gly_rama_ff31.py` -> `~/Downloads/gly_rama_ff31.png`, the ff3.0 figure's format: GLY|right|ALL
  (untrained) plus the mean change of the 37 GLY|X maps, each checked to be a sum of basin offsets (residual
  1e-6). `plot_gly_rama_nonhelix_zero.py` now draws four maps, ff3.1 lower left and BioEmu lower right, and writes
  `gly_rama_ff31_nonhelix_zero.png`; the other figures it and `plot_gly_rama_ff30.py` rewrote are unchanged in
  content. Files modified: findings.md (1.35), progress.md.
* **15:00-15:20 10-15, TM4 by DSSP (user: "I suspect you are not measuring it correctly").** Rebuilt
  `glpg_tm4_dssp_ff21_vs_ff30.png` from its own runs (runs_cov old-table seed 1, energies match) beside their
  fixed-table twins, and DSSP timelines of the three 79HIS REMD chains at T 0.80 (`~/Downloads/
  glpg_tm4_dssp_old_vs_fixed_table.png`, `glpg_tm4_dssp_remd_ff21_ff30_ff31.png`). ff2.1's unfolding in the
  user's figure is the old table (last quarter 0.644 -> 0.994 on the fixed table); in the REMD, DSSP puts ff3.1
  lowest (0.892 against 0.937 and 0.931), not the dihedral box's 0.97-0.98. Files modified: findings.md (1.35,
  10.27), remote_jobs.md (owed checks), progress.md; WATCH_STATUS.md on the cluster.
* **15:50-16:10 10-15, ff3.2 submitted (user: "same as ff3.1, except that it starts with ff3.0 instead of
  ff2.1"; Mac Studio, not the owner).** `$P/training/ff32_depth_data`: ff31_depth_data's trainer byte for
  byte, ff_3.0_gdepth's released files as `init_param` and as the depths' base library (plan.md Phase 17).
  Initialised and checked against ff31_depth_data's start and ff3.0's release (`checks/ff32_deploy_20261015/`),
  `--test-only`, then midway2 49406711, R from 16:02 (successor 49406726). `submit_new.sh` and `stock.sh` know its `f32_` tag (backups
  `.bak_pre_ff32_20261015`). Files modified: plan.md (Phase 17), remote_jobs.md, progress.md; on the cluster
  `submit_new.sh`, `~/watch_probe/stock.sh`, WATCH_STATUS.md.
* **00:13-01:45 10-17, the Mac Studio takes the jobs (user: "check remote job status, and start monitoring").**
  Watch cron `82dc7422`; caught up two unwatched days: ff3.2 steps 4-38 (steric catch-up job 49410086, no flag),
  panels f32_00 and f32_01, ff3.1 TM4 at block 2's end (findings 1.35), the ff3.1 TM4 HDX (midway3 60979543,
  60979545; findings 5.5), the nobp/ctrl38 scoring to 31 of 32 pairs and the training-set RMSD comparison (the
  run's `setup/notes/findings.md`), ff3.2's glycine maps (findings 1.21). Files modified: findings.md, plan.md
  (Phase 16), remote_jobs.md, progress.md; on the cluster WATCH_STATUS.md and the nobp run's notes.
* **01:50-10:10 10-17, the hybrid side-chain placement defect (findings 4.17; plan.md Phase 18).** Found by an
  agent's code survey for a deterministic-side-chain seed, then measured in the engine: every hybrid seed placed
  its side-chain beads at `rcf[li[li]]`. User: "fix the bug and rerun the affected simulations"; cancel all 8
  running chains, rerun all 10. Cancelled (STOP files first) and kept the old chains; fixed
  `py/martini_prepare_system_lib.py` (and the cluster copies); repaired 4 base and 8 per-force-field seeds and
  built ctrl38 and nobp 79HIS seeds (`checks/placefix_20261017/`); verified on the production engine (placement
  to 3.6e-6 A, 1-body exact, only the intended datasets changed; the fixed preparation reproduces the repaired
  seeds); submitted 49411315-49411324 after `--test-only`. single_state.py compares the libraries' geometry
  at the float32 the seed stores: the trained files differ from ff2.1's at float64 by 1 ulp (2.2e-16). Files modified:
  py/martini_prepare_system_lib.py, findings.md (4.17, 1.35, 5.5), plan.md (Phase 18, Known Errors),
  remote_jobs.md, progress.md; on the cluster the seeds, the two preparation copies, WATCH_STATUS.md.
* **10:30-23:00 10-17, held-out panels for ff_2.1_nobp against ff_2.1_ctrl38 (plan.md Phase 19; findings 1.21).**
  User: "Proceed" (the 44-domain panel), then "find something new into the validation set" (dt 0.009, nobp and
  ctrl38). The 44 domains (`nb_01`, `ctrl38`, which reproduces `c9_01` bitwise) resolve nothing. Built a second panel
  of 89 mdCATH domains (selection, range-request fetch, `panel.py prepare`, mdCATH 320 K references; `panel.py`
  gains an optional time step and `aa_mdcath`, backup `.bak_pre_mdcath_20261017`); 1gqoB00 excluded for its broken
  stored start. Result: no-BP folds 0.264 against 0.365 (-0.101 [-0.155, -0.051], p 7e-5; lower in 59 of 89), mean
  Q and CA-RMSD worse too, in every length bin and fold class. Files modified: plan.md (Phase 19), findings.md
  (1.21), remote_jobs.md, progress.md; on the cluster `ff3_selection/` (panel.py, `panel_mdcath.sbatch`,
  `panel_compare.py`, `mdcath_candidates/`, `{domains,aa,runs}_mdcath/`) and the nobp run's `setup/notes/`.
* **10-18, audit of the no-BP training protocol (user: check the protocol and that no-BP's training has ended;
  plan.md Phase 20).** Three read-only agents and the convergence tools, scripts and logs in
  `/project/trsosnic/yinhan/checks/nobp_audit_20261018/{A,B,C,main}/`. The comparison is clean: single-state library,
  configs, engine derivatives at n_rot 1, matching to the control, releases and the mdCATH numbers all check out; no
  no-BP training job can start. Two findings: `convergence_gate.py` fails no-BP's `rot` group at step 38 (control
  passes), so the verdict is at matched length, not at no-BP's fixed point (findings 1.21 "Audit"); and the angular
  gradient of the side-chain pair term is misattributed in master's engine, in every ConDiv run (findings 3.12). No
  job submitted or cancelled. Files modified: findings.md (1.21, 3.12), plan.md (Phase 20, Known Errors), progress.md.
* **08:00-08:50 10-18, the angular-gradient fix (user: "fix the issue and continue training if needed ... resolve
  them as well ... check ff3.1 as well"; plan.md Phase 21; findings 3.12).** Fixed `src/interaction_graph.h` (the
  [t2,t1] derivative of a symmetric interaction is the edge's derivative with the elements exchanged); built on
  midway2 from `$P/src` plus the fix and tested against the installed build (`checks/angfix_20261018/`): simulation,
  energies, forces and every other derivative bitwise equal, the angular entries now equal to central differences.
  Stopped ff3.2 in step 69 and the two step-37 glpG chains (STOP files); installed the build; continued no-BP and its
  control from step 38 with a gate-driven end (`ff21_nobp_dt009/setup/after_gate.sh`: gate, then the mdCATH panel or
  one more epoch to 114); reran ff3.1 and ff3.2 from step 0 as `_angfix` copies of their sources. Queued 49415222,
  49415223, 49415231, 49415233. The watch session (`upside2-md-f1`) told and took them. Files modified:
  src/interaction_graph.h, findings.md (3.12), plan.md (Phases 17, 21, Known Errors), remote_jobs.md, progress.md; on
  the cluster `$P/obj`, `$P/src/interaction_graph.h`, the two runs' `after_training.sbatch`, `submit_new.sh`,
  `~/watch_probe/stock.sh`, two new run directories, two STOP files.
* **11:00-12:10 10-18, `/project` disk reclaim (user: "cleanup the disk a little bit ... unused training
  trajectories, trajectories from wrong ff").** Surveyed `/project/trsosnic` (591 G free against ~0.75 T of growth due
  within a week) and proposed three groups; the user chose all three. Deleted 1.1 T on midway2 from manifests
  (`checks/disk_reclaim_20261018/`): half-written steps of stopped trainings (86 G), replica and HDX-copy trajectories
  of the defective glpG runs (821 G), 14 finished panels' `.up` (228 G); 1.7 T free after. Verified: every manifest
  file gone and nothing else; `panel.py select` reproduces the f32_02 table byte for byte from the kept npz; live logs
  clean; HDX readouts intact. The delete script stopped after group 1 on its first run (a `while` count under
  `pipefail` returned 1); fixed and rerun. Files modified: remote_jobs.md (§1 "Disk" and the affected records),
  findings.md (4.17), plan.md (Phase 18), progress.md.
- 21:40 10-18, the corrected glpG chains' block 1 TM4 read (Mac Studio, watch): all eight chains closed block 1
  18:52-20:39 with KE/1.5kT 0.989-1.008; `tm4_ff31.py`, `gly_flip_rungs.py` on midway2 login1 (40 s, 150 s) and DSSP
  through block-1 copies of the block-2 scripts (ff_2.1_sc dropped, its files deleted 10-18), mdtraj locally. TM4
  holds under both force fields; no TM4 glycine at phi > 0 at T 0.70. Files modified: findings.md (1.35), plan.md
  (Phase 18 item 7), remote_jobs.md (§1, resume checklist), progress.md.
- 21:55 10-18, handoff written for a cold start on another computer (user: "write everything on md files and expect
  me to run the monitoring on another computer tomorrow"). remote_jobs.md: header and §1 snapshot to 21:45 (21 jobs),
  "What happens next" (the no-BP gate ~00:00 10-19, the arms' chunks, the three PD trainings, glpG block 2), "Resume
  here" rewritten (reading list, owed checks with commands, the planned move and how to stop the Mac Studio's watch
  first, open questions), watch step 4 now covers the two ff2.1 gates and the mdCATH comparison, and the watch prompt
  rewritten for the running jobs. findings.md: 1.21 (the side-chain model moves the H-bond margin, opposite ways in
  no-BP and its control), 10.28 (name runs by what they are). WATCH_STATUS.md on the cluster carries the planned
  move. The Mac Studio's watch keeps running until the user closes its session.
