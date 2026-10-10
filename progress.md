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
