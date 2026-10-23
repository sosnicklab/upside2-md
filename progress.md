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
- 00:20 10-19 on, large-protein de novo investigation (plan.md Phase 22; user: "de novo folding on large proteins ...
  performs poorly"; hypothesis: co-translational folding, S. Fried). Read-only on the cluster. Training set 50-150
  residues (456), 44-domain panel 50-75. On the mdCATH native-start panel (ctrl38) the folded-fraction drop at 151-250
  is a metric artifact: Upside's mean Q sits 0.13 below all-atom in every bin and the all-atom Q spread narrows with N;
  with a fixed-offset cutoff there is no trend, Rg stays at 1.01 of native. Local test (Mac Studio, 00:38-05:00):
  T4 lysozyme, RNase HI, apoflavodoxin (in-vitro refoldable, 155-179 residues) and the homeodomain control under
  ff_2.1 at T 0.80, native x3, bulk x6, vectorial growth x6, 60,000 tu each, with quenches. Native starts hold; every
  de novo run, the control included, ends in a compact globule without native tertiary contacts, 150-180 E_up above
  native after quench (mostly missing H-bonds); vectorial equals bulk. A midway2 move was staged but not submitted
  (our own pending trainings put any new job at 10-20 09:45). Follow-up 05:20-06:37 at T 0.88 and 0.95 (2LZM, 2RN2,
  1ENH): no coil folds; RNase HI's native holds at 0.88, the control's melts. Proposal sent. Study directory moved to
  `scratchpad/large_protein_20261019/` (76 G, gitignored). Files modified: plan.md (Phase 22), findings.md (11e),
  up.md (§2.4, stale bb_env claim).
- 08:40 10-19 on, the `environment_coverage_hb` typing defect fixed (user: "Fix the bug, and write your findings in
  upside3.md"). The writer stores per-residue types; the engine sizes `aa_types` from its dataset and requires every
  id2 to index it. Pre-fix and fixed builds from one source (`envfix/`): old configs bitwise equal (8 configs x 3
  frames, 208 arrays); both builds match a numpy coverage reference to 5e-7; finite differences ratio 1.0001; energy
  change only in backbone desolvation, -0.1 to +1.3 E_up; the pre-fix engine rejects new configs. Installed in local
  `obj/` (backups `.bak_pre_envfix_20261019`), not on the clusters. A first dynamics check died: its Pool lacked a
  `__main__` guard under macOS spawn; the rerun (5,000 tu native starts) shows KE/(1.5 kT) 1.007-1.010 and Q, Rg
  within run noise of the pre-fix runs. upside3.md written (findings and the Upside 3 plan, gated on the
  side-chain verdict); findings 11e shortened to a pointer. Files modified: src/environment.cpp, py/upside_config.py,
  up.md, upside3.md (new), findings.md, plan.md, progress.md.
- 14:00 10-19 on, Upside 3 planning on the MacBook Pro (plan.md Phase 23; the user's five goals). Three read-only
  research agents. Their reports are in `scratchpad/upside3_planning_20261019/`; the agents could not write files, so
  the main session saved them.
  * **Audit.** Defect inventory, 18 master and 7 branch-only audit units, and a writer/reader registry. A first pass
    found:
    * the Mac/cluster finite-math flag split;
    * an ASan out-of-bounds read in `DihedralCoord`;
    * `Concat`'s uninitialised member;
    * unconfirmed candidates in `CurvatureChange`.
  * **Units.** 45 C++ and about 67 Python unit sites. E_up's kJ/mol value is FF2's 298 K / 0.85 calibration, and the
    master analysis scripts carry three conflicting kelvin conventions.
  * **GPU.** Per-node profile (M3 Max) with no single hot spot; the hybrid spends 86% in the dry-MARTINI pair loop.
    Recommended route: single-source kernels with many systems per launch, behind a prototype go/no-go. midway2's
    `gpu2` has K80s only.

  Also found: the MacBook Pro's `obj/` is the 10-02 build. upside3.md was rewritten for all five goals. Files
  modified: upside3.md, plan.md, progress.md, findings.md. Engine and Python sources were not changed.
- 14:10 10-19 on, Upside 3 audit started on martini-dev. The user made martini-dev Upside 3, with fixes there only
  (plan.md Phase 24).
  * Local `obj/` rebuilt from the branch.
  * Harness written in `tests/` (`fixtures.py`, `node_fd.py`, `regress.py`), plus `get_param_size` in the C API.
  * The first regression baseline did not reproduce: 2 of 6 dumps differed in the membrane fixture, one with NaN.
    ASan through Python localized it to ledger A2. `MembraneCBPotential`/`MembraneSurfCBPotential::get_param_deriv`
    wrote the right-hand spline derivative a whole table past `deriv` (master code).
  * Fixed A2 and A1 (`DihedralCoord` reading dummy atom -1).
  * Verification: ASan clean on all static records; the membrane parameter derivative's finite-difference median
    error is 1e-6 (83-97 of 220 directions failed before); regression changes only that array, and three dumps
    are bitwise equal.
  * The REMD reviewer and six unit reviewers are running (batch 1: interaction graph, environment, BP, H-bond,
    placement, writers).

  Files modified: src/coord_basic.cpp, src/membrane_potential.cpp, src/engine_c_library.{h,cpp},
  py/upside_engine.py, tests/ (new), plan.md, findings.md (3.14, 10.29), upside3.md.
- 15:20-18:00 10-19, audit fixes on martini-dev (plan.md Phase 24).
  * CLAUDE.md Development Rules reworded at the user's word: martini-dev is Upside 3, with deliberate departures
    only.
  * Reviewers reported: REMD, environment, interaction graph, BP, H-bond, writers, placement.
  * Fixed and verified A3-A28 (findings 3.14):
    * REMD: criterion 1 removed, swap scheduling, the momentum-swap reason (master heats the replicas), the guard
      removed, hashed per-system seeds, inner_step taken from the integrator, the barostat seed, NPT replica
      exchange refused;
    * pair-list cache invalidation, NonlinearCoupling, the lateral-pressure mapping;
    * writers: trans-cis maps, fix-rotamer, rot-scale, Rama and H-bond across chain breaks, the restraint group,
      integer index arrays;
    * engine: the doubled H-bond count, rotamer-count check, bad-solve count, extensive-logging names, pair-list
      capacity, eigensolver tolerance (+2% runtime), Mac build flags.
  * `tests/config_diff.py` added.
  * Regression baseline: `scratchpad/upside3_audit/w5.npz` on `fixtures_w2`.
  * Running:
    * the membrane and engine-core reviewers;
    * a Python option-fix agent;
    * the BP floor patch.
- 18:00-20:30 10-19, the user's three approvals applied (findings 3.14 A29-A34).
  * A29: exact backbone-burial derivatives, the backbone weights' gradient added to `envw` in both ConDiv copies.
  * A30: exact H-bond boundary derivatives, held frozen in the trainer's new `hbr` field.
  * A31: the coverage stubs raise "not implemented".
  * A32: dead nodes and options deleted.
  * A33: the angular cosine held to [-1,1].
  * A34: the BP floor patch.

  Also: the Python option batch (13 options fixed by an agent); the dead `rl_cavity_radial.py` and
  `write_sidechain_radial` deleted; make-unbound seeded; example 07's AFM argument fixed.

  Reviews reported: membrane (8 findings) and engine core (9). Their verified patches are in preparation by
  their reviewers. The usage limit stopped all agents at 20:00; they resumed at 20:20.

  Regression baseline: `w10.npz` on `fixtures_w3`.
- 21:00 10-19 to 01:40 10-20, audit fixes continued.
  * Applied the engine-core patch (A42-A49), the coordinate patch (A51-A55), the I/O fixes (A56-A60) and the
    analysis patch (A61-A65).
  * My own A26 regression (capacity check too strict, A50) fixed.
  * My CMake slicing error found and repaired, lesson findings 10.31.
  * Analysis: the hybrid HDX accessibility was off by one residue; the arrays to regenerate are listed under plan.md
    Known Errors.
  * Dry-MARTINI engine review: the choline charge is missing from the side-chain/lipid table, and the prep barostat
    never scaled lipids. Both are recorded under Known Errors; the patch is in preparation, with the martini.h5
    rebuild the user approved.
  * The usage limit stopped agents twice (20:00, 01:00).
- 01:40 to 02:20 10-20, hybrid coupling and stages: the review reported 4 defects; applied its patch (A67-A69). The
  barostat finding went to the dry-MARTINI patch.
  * Files: `src/martini_fix_rigid.cpp`, `martini_hybrid.cpp`, `martini_internal.h`, `martini_stage_params.cpp`,
    `martini.h`, `main.cpp`, `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, `up.md`, CLAUDE.md
    (group list).
  * I added: a singular rigid-group inertia is refused instead of silently zeroing the rotation.
  * Checks on the repo build: regression w18 against w17, 0 of 545 arrays differ; 7.0 MD with minimisation
    identical bit for bit to the reviewer's build; 6.6 differs from it by one float ulp after one minimisation
    iteration (fast-math reordering around my edit; each build reproduces itself exactly); "Production" refused.
  * Memory note on the hybrid `pos` layout corrected (4i+3 is O; BB beads follow the backbone block).
- 02:20 to 08:10 10-20, hybrid preparation: the review found 7 defects in every hybrid result (Known Errors);
  applied its patch and its minimiser fix (A70-A73). The usage limit stopped agents 03:40 to 06:20.
  * Files: `py/martini_prepare_system.py`, `py/martini_prepare_system_lib.py`, `example/16.MARTINI/` (readme, run.py,
    two shell drivers), `src/martini_hybrid.cpp`, `martini_internal.h`, `martini_potential.cpp` (minimiser),
    `main.cpp` (help text), up.md.
  * Mine: atom roles required for a rigid prep; the minimiser's step floor reported as a stop.
  * Checks: regression w19 against w18, 0 of 545 differ; a patched 1rkl 6.3 stage bitwise equal to the reviewer's
    build. The reviewer's end-to-end 1rkl prep: KE/1.5kT 0.97-1.00 from 6.4, production 1.011, net charge 0.
  * User 10-20: per-lipid diffusion targets, yes. Applied the follow-up (overlap_settling deleted, per-lipid table,
    carriers summed per lipid type; DOPC builds bitwise unchanged). POPE/POPG values await the literature search.
  * Cluster (midway3, read-only): August glpG preps did not overheat (KE/1.5kT 0.83-0.87) but carry the restraint
    offset; the live REMD seeds carry the fixed rotamer 1-body.
- 08:10 to 09:00 10-20, dry-MARTINI, examples, swap set, docs.
  * Applied the dry-MARTINI patch (A74-A76): barostat, box, tables, bond spline, dead code; the side-chain/lipid
    table's choline charge, mirrored angle, CB offset, cutoff tail and cubic spline. Installed the rebuilt
    `parameters/ff_2.1/martini.h5` (md5 `cc90a933...`); old kept as `martini.h5.bak_20261020`.
  * The examples agent's fixes were already in the repo (A78). Its leftover REMD item fixed by me (A77): the swap set
    from the round; split run 10 of 10 frames equal to the continuous run (old build 1 of 10).
  * POPE/POPG diffusion values from the literature search (findings 2.5): POPE 8.6, POPG 11.5 um^2/s.
  * up.md: stale CGL nodes, groups and outputs removed (no code reads them); CLAUDE.md file-format list updated.
  * Regression w21 against w20: 0 of 545 differ.
  * Running: a fresh 1rkl prep with everything (box test), the angle-table follow-up, the hybrid analysis review.
- 12:30 to 13:00 10-20 (after the 09:30-12:30 usage limit): the 1rkl box test finished (plan.md Phase 24: healthy
  except the expanding box); applied the angle-table patch (A79), regression w22 0 of 545; smoke run on the new
  binary KE/1.5kT 0.967. Agents resumed: box diagnosis, analysis patch.
- 13:00 to 13:20 10-20: applied the hybrid analysis patch (A80; 16 files in py/ and example/). Pass 1 of the audit
  has now reviewed every unit. Open: the box expansion (diagnosis running).
  * apply_fix_rigid_md's unread force zeroing deleted (poison test bitwise), index ranges validated at build; w23 0 of 545; stale src/main.cpp.orig from a patch removed.
- 13:40 10-20: box expansion diagnosed (no engine defect; CHARMM-GUI step5 bilayers start interdigitated, DOPC 93%, POPE/POPG 88%); leaflet-separation fix and area validation running. Live glpG seeds: about 71.5 A^2 box per lipid, P-P 42.5 A (constructed template).
- 13:55 10-20: user offered CHARMM-GUI packages and any resource; GROMACS reference runs started (DOPC and POPE/POPG dry-MARTINI, 1.1 and 1.2 nm); leaflet-translation prep patch on hold.
- 17:30 10-20: usage limit 14:00-17:30 stopped the GROMACS reference and the Upside bilayer validation; both resumed (no processes were left running).
- 23:15 10-20: the 17:30-18:50 session (usage limit until 22:30) reported nothing. Both agents spent it re-reading
  very large contexts while polling for simulations that had finished (findings 10.30). Every GROMACS reference run
  and every Upside bilayer validation run completed; the analysis now runs as a detached job
  (`scratchpad/detach.sh`), and the small outputs are copied to `scratchpad/bilayer_ref/`. Old agents not resumed.
  * Results: GROMACS shows the same early box jump from CHARMM-GUI's overlapped start and recovers over ~12 ns, so the
    Upside prep box expansion is the start, not the engine (plan.md Phase 24, findings 8 reference). Upside pure
    bilayers (0.55 ns) not converged. Decision: prep from GROMACS-equilibrated membranes; leaflet-translation patch dropped.
  * 23:29 detached campaign `scratchpad/bilayer_ref/upside_from_gmx_eq/run.sh` (10 threads, ~4-5 h): Upside pure
    bilayers from the GROMACS equilibrium, glpG and 1rkl preps from GROMACS-equilibrated against step5 membranes.
    DONE 02:34 10-21: DOPC holds the GROMACS equilibrium, POPE/POPG 1.3% smaller at 0.73 ns; 1rkl healthy from the
    equilibrated membrane, interdigitated in production from step5; glpG fails the solvation gate from both
    (plan.md Phase 24, findings 8). Files modified: findings.md (8), plan.md (Phase 24).
- 22:30-23:30 10-20, status and convergence review (user: "check remote job status"; then "check if the finished ffs
  needs no more trainings" and "check if the deterministic / probabilistic side chain have a conclusion"; plan.md
  Phase 25). Status read only: the control converged at step 76 and `mdc_ctrl_03` folds 0.353; the reruns at steps
  55/54 (ff3.1 rerun margin +0.0987 at step 52, a hold; a steric flag on each, beta' + pPII' at one step); glpG block 3
  running; all benchmark arms complete. One job submitted, `panel_compare` 49435922 (4 min, read only). New read-only
  tools in `checks/converge_20261020/` (`gate_window.py`, `gram_probe.py`, `ang_pull.py`, `depth_pull.py`,
  `paired_npz.py`); `gate_window.py` reproduces the recorded gates and `paired_npz.py` the step-37 `panel_compare`
  exactly. Results: ff_3.1_depth_data and ff31_depth_all pass the gate and their depths show no pull; the defect left
  no detectable pull at the control's or ff3.0's fixed point; ff_3.0_gdepth's glycine alpha_L depth is pulled one way
  every epoch (under the cut over the last); no-BP's ~110 angular coefficients are still pulled at step 76; at the
  fixed points no-BP folds -0.111 [-0.160, -0.065] (p 1e-5). Files modified: findings.md (1.21, 1.36), plan.md
  (Phases 21, 25, Known Errors), upside3.md (status), remote_jobs.md; WATCH_STATUS.md on the cluster.
- 10-21 00:30-03:10, the side-chain comparison on one metric (user: email drafts to Tobin and John; "I prefer to use a
  single metric (TM score?) for all tests"; plan.md Phase 26). Jobs: `tm_mdcath` 49436371 (TM of the converged mdCATH
  runs), the 44-domain panel rerun at dt 0.009 for both converged force fields (49436374, 49436375; all runs finite,
  KE/1.5kT 1.00-1.01), `tm_panel` 49436485, and the last Peng pair's score (49436488). Mean TM: mdCATH -0.034 (p 7e-8),
  44 domains -0.026 (p 0.007), Peng natives -0.029 (p 0.18), de novo +0.010; the TM > 0.5 shares do not differ.
  References checked on the web (Peng JCTC 2022, Charron Nat Chem 2025, CATH 4.2 NAR 2019, mdCATH Sci Data 2024).
  Lesson recorded: findings 10.32 (short results emails). Files modified: findings.md (1.21, 10.32), plan.md
  (Phase 26), remote_jobs.md, progress.md; new `tm_panel.py`, `tm_panel.sbatch`, `panel44_dt009.sbatch` in
  `checks/converge_20261020/`.
- 01:20 and 02:00 10-21: merged the user's pushes 9871ca81 (the 22:30-23:30 session above) and 138e6f12 (TM-score
  comparison) at the user's request: HEAD moved to each with the local edits kept unstaged. The other machine numbers
  from the committed plan.md, so its Phases 23 and 24 became 25 and 26 (plan.md, progress.md, remote_jobs.md) and its
  lesson 10.29 became 10.32; upside3.md's status takes the verdict line. Pre-merge backups in
  `scratchpad/pre_pull_20261021{,b}/`.
- 02:40-03:35 10-21, "continue debug" (plan.md Phase 24). Fixed and verified: A81 (HDX peptide numbering), A82
  (zero force where the angular cosine is clamped; the null-select build proves the 138 moved arrays are code
  placement), A83 (the hybrid dropped every ALA and GLY side chain from the protein's side-chain node: glpG 42 of
  210 residues), A84 (dead `Const1D`/`Const2D`, duplicate `WallSpring`). Harness: `tests/registry.py` (69 nodes),
  `tests/line_integral.py`, node_fd with measured noise and a third step; fixtures regenerated (`fixtures_w4`),
  baseline `w26.npz`. The rotamer FD candidates are refuted (float32 noise, kinks). glpG's bare belt: the gate
  counts side chains as gap and the insertion void heals slowly (prep NPT 4.6 ps against CHARMM-GUI's 10 ns); NPT
  convergence test running. Engine changes verified in `scratchpad/upside3_audit/clampfix_home`; obj/ not yet
  rebuilt. Files modified: src/bead_interaction.h, src/martini_potential.cpp, src/coord_const.cpp,
  src/sidechain_radial.cpp, py/martini_prepare_system_lib.py, py/upside_nodes.py,
  example/00.AnalysisScripts/4.calc_D_uptake.py, tests/node_fd.py, tests/registry.py (new), tests/line_integral.py
  (new), up.md, findings.md (3.14 A81-A84, 8, 10.12), plan.md.
- 10-21 03:00-03:35, the Boltz-2 HDX sampler, Phase 0 (user: "develop a modified Boltz to replace upside (for HDX
  sampling)"; "we are building it on Boltz, not BioEmu"; "run on midway3 if it is the only one available").
  `diffusion_dynamics.md` rewritten around Boltz-2 with HDX as the goal: what the equilibrium HDX estimator needs from
  a sampler (findings 5.1-5.2), Boltz-2's conditioning at `b1ebfc4` (an `md` tag, no temperature input, mdCATH in its
  MD training data). midway2 cannot host it (K80s only for the group; its one V100 node down since 2024-06). Set up
  on `/beagle3` (`/project` metadata 60x slower): Boltz at `b1ebfc4`, venv with torch 2.14.1+cu126 (PyPI's CUDA 13
  build lacks sm_70), weights, MSAs and references. First test queued on the midway3 `gpu` partition: 61216915-17,
  one each on A100, V100, RTX 6000 (ubiquitin and EHEE_rd2_0005; 5 default and 100 `md` samples). Lessons: findings
  10.34 (Boltz, never BioEmu) and 10.14 (venvs on `/beagle3`, the CUDA 13 torch wheel, `pkill -f` over ssh). Files
  modified: diffusion_dynamics.md, findings.md, remote_jobs.md, progress.md; cluster scripts in `$B` with copies in
  `scratchpad/boltz_hdx/`.
- 03:35-08:40 10-21, debug continued (plan.md Phase 24).
  * Fixed and verified:
    * A85, logged kinetic energy degrees of freedom;
    * A86, a refused config no longer segfaults Python;
    * A87, environment weights mode 1 deleted;
    * A88-A93, Python writer batch from the option fixtures (fix-rotamer on chainless `.chi`, dead options and
      writerless placement nodes deleted, input checks that fire, example 07's frame, integer group ids).
    * Refuted: the heavy-atom coverage node is dead (HDX reads it).
  * Harness:
    * fixtures for every writer option group (`fixtures_w6`, 28 configs);
    * node_fd ties NonlinearCoupling;
    * `regress.py compare` lists one-sided keys apart;
    * `invariance.py` takes its SE from pymbar's statistical inefficiency.
  * obj/ rebuilt with A82-A93; regression equals the o12 baseline bitwise (w30/w31, 0 of 4073).
  * Equipartition passes in both systems. The hybrid protein "deficit" was seed 2's noise draw sampled by two
    builds (same-seed r 0.92-0.96); the lipid FAILs were block SEs shorter than the 55-frame correlation.
  * glpG: prep NPT converges the box in 30k steps and the belt plateaus at 8 A by 60k. Healthy interfaces have at
    most 1% of belt sites beyond 8 A, glpG's prep 45-51%, so the void is real; production continuation running.
  * D9 documented (needs an engine-side exclusion id).
  * Pass 2 started on the six trained-node units. The session limit stopped all six reviewers at about 07:45, and
    two were resumed at 08:31 (lesson in findings 10.30).
  * Files modified: src/engine_c_library.cpp, src/environment.cpp, src/placement.cpp, src/main.cpp, py/upside_config.py,
    py/run_upside.py, py/advanced_config.py, example/07, example/09, tests/{fixtures,node_fd,regress,invariance}.py,
    up.md, findings.md (3.14 A85-A93, 10.30, 10.33), plan.md.
- 10-21 08:30-08:45, the Boltz-2 HDX sampler: first test read, BioEmu2 code allowed (user: "Yes it can use code
  from bioemu2"). The three jobs ran early (04:37-07:44), about 7 min each. All three GPU types run bf16 with the
  cuEquivariance kernels; about 1 s per ubiquitin sample. Default samples are near-native with clean geometry;
  `md` samples at `step_scale` 1.0 are wider (ubiquitin pairwise CA-RMSD 2.0 A) but have bond-length scatter 4x
  crystal spread and clashes in 86-97% of samples; tag and step scale were changed together, so the cause is open.
  The 2x2 separating test (user: "Yes.") is queued on a V100, `tag_stepscale_test.sbatch`. Phase 3's lag conditioning now ports BioEmu2 `f99e0d8` code,
  paths verified in a local read-only clone. Files modified: diffusion_dynamics.md (status, sections 1 and 3, Phase
  0 readout, Phases 1 and 3, risks, references), findings.md 10.34, remote_jobs.md, progress.md, memory
  `feedback_boltz_not_bioemu.md`.
- 10-21 11:10-12:00, remote status, the Peng numbers and released ff2.1 (plan.md Phase 27). Read-only pass: both
  `_angfix` reruns at step 76 (the ff3.2 rerun's margin -0.0024 at 75), steric at 75 clean beyond the known beta' +
  pPII' flag, 8 glpG chains in block 3. Checked against Peng's SI: both step-37 side-chain arms sit below the
  published FF2 (0.55 / 0.42) because alpha3d and homeodomain unfold at their cold rungs; the scorer and protocol
  reproduce FF1 and gave the retired ff3.0 FF2-level numbers (findings 1.21). The probabilistic arm is a retrained
  ff2.1 (no trained file equal to the release), so released ff2.1 was submitted on all three tests: mdCATH and
  44-domain panels at 0.009 with dependent TM jobs, Peng's 30 missing arms, and a self-requeueing collector that
  scores them and writes the three-way table. Then (user: does it affect the other trained force fields?) every
  epoch end of every run on the 44-domain panel against ff2.1: all fold 0.13-0.23 less, mostly from epoch 0
  (findings 1.37); the owed panels and the Peng scores of ff3.0 / ff3.1 were queued with jobs that write the final
  tables. Files modified: findings.md 1.21 and 1.37, plan.md Phase 27, remote_jobs.md; new
  `scratchpad/ff21_released_20261021/` (cluster copy `checks/ff21_released_20261021/`).
- 08:40-10:55 10-21, Upside 3 audit continued (plan.md Phase 24; user: "git pull it and merge", "proceed with your
  suggestions", "record progress"). Pulled 9a88a139 (other machine: Boltz HDX sampler, TM-score phase): HEAD moved,
  local edits unstaged, its lesson 10.30 renumbered 10.34 and its "Phase 24" kept as 26 (backup
  `scratchpad/pre_pull_20261021c/`).
  * Index oracle (`tests/index_oracle.py`) clean but for D9; found A94; up.md corrected in 16 rows.
  * Pass 2 reported for BP, environment, interaction graph and H-bond. Applied and verified:
    * A95, the BP stop test;
    * A96, environment input checks;
    * A97, dead radial, contact and linear-coupling nodes and the exclude-residue options;
    * A98-A100, interaction-graph sens, slot-30 derivative and input validation;
    * A101, the H-bond cutoff per pair.
  * Placement and writers were stopped by the session limit.
  * Equipartition passes with the SE from pymbar's g (lesson 10.33).
  * glpG: the bundle interior inflated the old gate; the lipid face was really bare at 18% and production heals it
    (findings 8b). Gate and 6.6 length changed (A102). The first prep with them passed only at 9.9%, so it is not yet
    validated.
  * obj/ holds A97; the A101 build is in `scratchpad/upside3_audit/bp_home`, baseline `w36_hb.npz`.
  * A pause of local compute is armed for 14:05-14:15 (user).
  * Files modified: src/{rotamer,environment,sigmoid_energy,sidechain_radial,hbond}.cpp,
    src/{interaction_graph,bead_interaction}.h, py/{upside_config,run_upside,martini_prepare_system}.py,
    example/16.MARTINI/{run_sim_hybrid.sh,run.py,readme.md}, tests/invariance.py, up.md, findings.md (3.14 A94-A102,
    8b, 10.30, 10.33), plan.md (Phase 24 Resume here), upside3.md, CLAUDE.md.
- 10-21 12:30-13:40, ff3.1.1 and ff3.2.1 (plan.md Phases 28-29). Single-file resets of the control's epoch 0: the
  loss is spread over side-chain table, backbone burial and H-bond (findings 1.37). Released ff2.1 at dt 0.009 folds
  as at 0.015, so the loss holds at the training time step. Peng and Jumper stopped on a fixed budget on purpose
  (findings 1.23). ff3.2 rerun converged; its glycine map plotted (`scratchpad/plot_gly_rama_ff32_nonhelix_zero.py`).
  The depth maps stay at their starts because per-map priors hold a pull that is consistent only pooled (findings
  1.34). Mode `pooled` added to `training/ff31/rama_basin.py` (one shared offset per central residue and basin plus
  per-map deviations), checked against `data` to 1e-17, width 0.5 by 10-fold cross-validation; ff3.1.1 and ff3.2.1
  built from the angular-fix reruns with only that file changed and submitted (user: everything else trains).
  Files modified: training/ff31/rama_basin.py, findings.md 1.23/1.34/1.36/1.37, plan.md Phases 27-29,
  remote_jobs.md; new `scratchpad/ff311_20261021/`, `scratchpad/ff_3.2_angfix/` (cluster `checks/ff311_20261021/`).
- 10-21 13:38-14:10, ff3.1.1 / ff3.2.1 revised (user: "resolve the accuracy lost bug as well in this run"): the
  every-term chains cancelled 3 min in; rebuilt as `ff311_depths_only` / `ff321_depths_only`, ff2.1's files with every
  Adam rate zero (`freeze_groups.py`) and only the pooled depths training, from NDRD's and ff3.0's map; submitted
  13:42 with their collector. The full 44-domain table and released ff2.1's TM there (0.681 against 0.650 and 0.623)
  recorded. remote_jobs.md regrouped by the user's four lines of work (header, §1, the watch prompt; finished rows to
  "Finished runs"; the obsolete HDX spec removed). Files modified: findings.md 1.37, plan.md Phases 27-29,
  remote_jobs.md; new `scratchpad/ff311_20261021/{freeze_groups.py,build_depths_only.sh}`.

- 10-21 14:23-14:50, Upside 3 audit continued (plan.md Phase 24 Resume here). obj/ rebuilt with A1-A101, regression
  equals w36_hb (w37_obj, 0 of 4073). Pass-2 reviewers for placement and config writers restarted from the updated
  brief. glpG gate run analysed by the gate's own fractions: the lipid face heals to 0.7% within 20,000 production
  steps; the first prep's 6.6 ended at 13% against the earlier prep's 6%, the two diverging from the first
  minimisation through A83's 42 ALA/GLY rows; seed replicates launched. D5 applied (the MARTINI pair spline clamped to
  the analytic end slopes; `solve_clamped_1d_spline` takes the slopes, the stage writer stores
  `combined_energy_end_slopes`): reproducer 2.7e-2 to 1.5e-7 E_up/A at the cutoff, fixtures_w6 unchanged; hybrid
  fixtures regenerating on the D5 build for its node_fd and hybrid regression. Lesson 10.35 (user: concise edits in
  master's style). Files modified: src/{spline.cpp,spline.h,martini_potential.cpp}, py/martini_prepare_system_lib.py,
  up.md, findings.md 10.35, plan.md Phase 24.
- 10-21 14:50-16:25, Upside 3 audit continued. Both pass-2 reviews in (placement, config writers; reports saved). The
  user's four decisions taken by measurement (plan Phase 24 item 5) and applied with every queued fix that was
  ready, one at a time, each reproduced before and verified after: D5 closed on regenerated 1rkl hybrid fixtures
  (A103), generators deleted (A104), `membranelateral_potential` deleted (A105, its forces sum to up to 0.1 E_up/A:
  no energy exists), chain breaks (A106), insertion codes (A107), placement D1-D3 (A108-A110), writers items 1, 3-12
  (A111-A119). `obj/` rebuilt 16:18 (A101 kept as `obj/*.bak_a101`), regression equal to its source build;
  `tests/regress.py` runs hybrids without pivot moves; fixture baseline `fixtures_w13`, hybrid baseline
  `h41_obj.npz`. glpG seed replicates: gate 9.2% pass, 16.1% fail; the gap sits at the belt ends and wanders
  7-16% under the rigid protein; s13's refused state continued in production to test whether it heals. Files
  modified: src/{eig,placement,membrane_potential,main}.cpp, src/CMakeLists_{M1,Other}.txt (generators deleted),
  py/{upside_config,advanced_config,upside_nodes,run_upside,PDB_to_initial_structure}.py, tests/{fixtures,
  index_oracle,registry,regress}.py, example/08 (lateral script, table, readme), up.md, upside3.md, findings.md
  (3.14 A103-A119, 10.35), plan.md Phase 24.
- 10-21 16:25-16:50, Upside 3 audit continued. git pull (user: "git pull and continue"): fast-forward to 3224a0bf;
  findings.md and progress.md, edited both locally and upstream, merged by hand (upstream's lesson and three Boltz-HDX
  entries kept), no stash. s13's refused glpG state heals in production (16.1% to 4.5%): the rigid-stage gate fails
  healthy preps; options to the user (plan Phase 24 item 3). Fixed and verified one at a time: writers items 2 and 15
  (A120-A121), `HBondEnergy` input checks (A122), BP D4 log split (A123), the dead-code batch (A124); `obj/` holds to
  A123. Reader preconditions (Distance ids, WallSpring wall type and ids, Surface included residues) under test.
  Files modified: src/{main,hbond,rotamer,bead_interaction.h,interaction_graph.h,spline.h,spline,sigmoid_energy,
  coord_basic,spring,surface}.cpp, py/{generate_restart_config,run_upside,upside_config,upside_nodes,
  PDB_to_initial_structure}.py, tests/regress.py, example/08 scripts 4-5, up.md, findings.md (3.14 A120-A124, 8b),
  plan.md Phase 24.
- 10-21 16:50-17:40, Upside 3 audit continued. Gate option (a) applied (user): the glpG solvation gate measures
  production's steps 15,000-20,000 and `run.py` starts REMD from that gated production (A126). Coordinate-node pass 2
  reported; fixed one at a time, each reproduced first: MovingConst3D's clock under g-JF and mv's slow level (A127),
  26 unchecked coordinate-reader inputs and the empty group (A128), reversed residue ranges (A129), Const3D and
  MovingConst3D carrying their values with the auto tip written from the config (A130), eight writerless nodes
  deleted (A131), dim1 against true widths (A132; corrects A124/A125's padded-width checks). `src/` to A132, `obj/`
  to A125. Files modified: src/{coord_basic,coord_const,coord_map,coord_operator,main,spring,sigmoid_energy}.cpp,
  py/{martini_prepare_system,advanced_config,upside_config,upside_nodes,generate_restart_config}.py,
  example/{06.PullingSimulation/0.run.py,16.MARTINI/run.py,16.MARTINI/readme.md}, tests/{fixtures,index_oracle}.py,
  up.md, findings.md (3.14 A126-A132, 8b), plan.md Phase 24.
- 10-21 17:40-18:10, the engine-graph pass 2 (P1-01) reported; D2-D5, D7 and D6's documentation fixed (A133-A137),
  D1 `mag` of a zero vector (A138; reproduced with builds without -ffast-math: nan before, finite after; the hybrid
  moves at rounding). `obj/` rebuilt with every fix to A138; new baselines w45_obj/h45_obj on fixtures_w16. Second
  reviews of P1-10 (driver, REMD, thermostat) and the dry-MARTINI engine started. Files modified: src/{main,
  deriv_engine.h,vector_math.h}, up.md 1.3, findings.md 3.14, plan.md Phase 24, pass2/brief.md.
- 10-21 21:45-22:35, pass 2 of the driver (P1-10) and the dry-MARTINI engine (P2-01, -03, -04) reported
  (`pass2/driver/report.md`, `pass2/martini/report.md`). Fixed and verified in two builds: MARTINI reader
  preconditions and the dead PO4 z reference (A139, A140), double sums of the MARTINI potentials (A141; hybrid 20
  of 196 arrays move, all potential values), rigid-group degrees of freedom in the logged KE (A142); the annealed
  continuation temperature, the synchronized stop, the frame taken before the round's MC move and re-centering,
  option checks, replica labels across restarts, exact rejected exchanges (A143-A148; fixtures_w16 218 of 4070
  arrays move, all exchange-run dynamics or frames on the MC grid). `obj/` holds every fix to A148 (`*.bak_a138`,
  `*.bak_a142`); baselines w47_obj/h47_obj. Second reviews of P2-02 (hybrid coupling) and P1-08 (local backbone)
  started. Files modified: src/{main,martini_potential,martini_brownian,martini_hybrid,martini_fix_rigid}.cpp,
  src/{martini.h,martini_internal.h}, py/generate_restart_config.py, examples 02, 03, 04, 06, 08 (5 scripts), 10,
  13 (2) continuations, up.md 1.1, 1.5, 1.7, findings.md 3.14, plan.md Phase 24, pass2/brief.md.
- 10-21 23:00-10-22 00:10, driver pass 2 item 1 (A149): the OU thermostat step moved into the leapfrog's
  `integration_cycle`, on the momentum at the coordinates; g-JF dispatched from `main.cpp`. First the four kick-drift
  loops merged into `leapfrog_step` (bitwise on fixtures_w16, the hybrid, v and mv runs). Harmonic tethers at tau
  0.05: v 0.920 to 1.0008 (exact 1.00008), mv slow level 0.792 to 1.0010 (exact 1.00073); REMD seeds 1-3,
  equipartition and pivot MC PASS; fixtures_w16 224 of 4070 arrays move (seeded-run dynamics), the hybrid bitwise.
  New `tests/invariance.py thermostat` (21 s, PASS). `obj/` to A149 (`*.bak_a148`); baselines w48/h48. Files
  modified: src/{deriv_engine.cpp,deriv_engine.h,main.cpp}, tests/invariance.py, findings.md 3.14, plan.md
  Phase 24, upside3.md, pass2/brief.md.
- 10-22 00:40-01:30, local backbone pass 2 (P1-08) reported (`pass2/localbb/report.md`); fixed A150 (the rama maps'
  `set_param`, refused since A57/A134, so a ConDiv step could not finish with the repo engine; the midway2 campaigns
  run their own older binary), A151 (13 unchecked local-term indices), A152-A153 (`--secstr-bias` residue range,
  non-proline trans-cis tables), A154 (dead spring block, dead `@sheet_eps`, up.md log columns). Reproducers match
  the reviewer's, node_fd 194 rows clean on 11 fixtures, regression bitwise on w48/h48; `obj/` to A154
  (`*.bak_a149`). D3 (`--secstr-bias` not periodic in phi) is with the user. P1-09 membrane review started. Files
  modified: src/{rama_map_pot,backbone_steric,tension,spring}.cpp, py/upside_config.py, up.md 1.3 and 1.7,
  findings.md 3.14, plan.md Phase 24, upside3.md, pass2/brief.md.
- 10-22 02:10-02:50, after the usage limit: the three stopped reviewers restarted fresh (P2-02 told to re-derive
  its first attempt's saved work, P1-11; P1-09 waits). P1-17 build read by the coordinator (leads in plan.md).
  `--secstr-bias` decided for the user (A155): one periodic phi window, phi < 0 continued to +120, each window the
  edge sigmoids summed over periodic images; the engine energy smooth through +-180 (forces off by 0.1 E_up/A at
  most, from 9.7); only `chig_rama` moves. The rebuilt fixtures are `fixtures_w17`. Files modified:
  py/upside_config.py, findings.md 3.14, plan.md Phase 24, upside3.md.
- 10-22 03:00-04:00, Monte Carlo pass 2 (P1-11) reported (`pass2/mc/report.md`); fixed A156 (dead `prob_mod` that
  broke jump detailed balance), A157 (10 unchecked MC indices), A158 (pivot bins drawn and scored on one double
  cdf; deliberate departure). A 2x drop of pivot acceptance in a short MD+pivot chig control was chain noise: the
  full-potential pivot-only chains agree (0.243 +- 0.028 against 0.249 +- 0.033). `obj/` to A158 (`*.bak_a155`);
  baseline w50 on fixtures_w17. up.md 1.7 lists `pivot_stats`, `jump_stats`. P1-09 membrane review started. Files
  modified: src/{monte_carlo_sampler.cpp,monte_carlo_sampler.h,random.h,main.cpp}, up.md 1.7, findings.md 3.14,
  plan.md Phase 24, upside3.md, pass2/brief.md.
- 10-22 02:30-03:10, hybrid coupling pass 2 (P2-02) reported (`pass2/hybrid/report.md`); fixed A159 (the rigid
  protein re-referenced every round, drifting 1e-2 A per stage), A160 (the PO4 hold on the round clock,
  `@po4_z_hold_end_round`; positions-only starts hold nothing; `generate_restart_config` carries the box), A161
  (barostat centre over real particles, bond spline a function of the box: NPT continuations bitwise), A162 (silent
  hybrid defaults refused, `preprod_protein_mode` removed), A163 (`AngleSpring` guard removed). Reviewer's
  reproducers reproduced; on converted fixtures (`fixtures_hybrid_w2c`) the A158 engine gives h48 and the A163 one
  moves only stage 7.0's seeded run; fixtures_w17 bitwise. `obj/` to A163 (`*.bak_a158`); baselines w51 (= w50),
  h49. The hybrid fixtures regenerating through the full prep (`fixtures_hybrid_w3/`, detached). Files modified:
  src/{martini_hybrid,martini_fix_rigid,martini_stage_params,martini_brownian,martini_potential,box,main}.cpp,
  src/{martini.h,martini_internal.h,box.h}, py/{martini_prepare_system,martini_prepare_system_lib,
  generate_restart_config}.py, tests/{invariance,registry}.py, up.md 1.4-1.5, findings.md 3.14, plan.md, upside3.md.
- 10-22 03:20-04:15, membrane and surface pass 2 (P1-09) reported (`pass2/membrane/report.md`); fixed A164 (the
  surface rotation scan: scores now rotation invariant; deliberate departure for the channel potential), A165 (the
  curvature changer needs a curved membrane and samples R above the half-thickness), A166 (11 unchecked readers),
  A167 (`--memb-scale` on the channel's inner tables), A168 (excluded residues' H/O rows), A169, A170 (dead surface
  columns), A171 (dead inputs), A172 (the channel potential flat-membrane only; `1rkl_channel` rebuilt flat), A173
  (up.md). M2 and M3 (ff2.1 table properties) left for the membrane retraining. The reviewer's changer patch held
  stale hunks that would have reverted A160; its curvature parts applied by hand. Verified: rotation invariance on
  four channel configs, FD 127 rows clean, the engine on the old set bitwise against w51 but for the surface
  width, the rebuilt set moving only 1rkl_surface (its H/O rows), hybrid bitwise. `obj/` to A173 (`*.bak_a163`);
  `fixtures_w18`, baseline w52. Files modified: src/{membrane_potential,surface,main}.cpp, py/upside_config.py,
  tests/fixtures.py, up.md 1.3 and 2.5, findings.md 3.14, plan.md, upside3.md, pass2/brief.md.
- 10-22 04:20-07:30, `git pull` (user: "git pull and merge with the recent change"): fast-forward to 6a320de5 (the
  Boltz-2 notes), my edits to findings.md, progress.md, remote_jobs.md restored and re-applied, nothing staged. The
  hybrid fixtures regenerated through the full prep on the A163 engine (`fixtures_hybrid_w3`): minimisers run their
  budgets, KE/1.5kT 6.6 1.003, production 0.99, C-N max 1.86 A, node_fd clean but one known rotamer flag; new hybrid
  baseline h50. Shipped-parameter pass 2 (P1-18) reported (`pass2/params/report.md`): documentation fixed (A174); D1,
  the SC/environment table's MARTINI side-chain beads placed on Upside's trained rotamer point and vector (ARG's
  charged bead 0.5-2 A from CB on the backbone side; table energies off by up to 21 E_up) and D7 dead parameter
  content held for the user, D1 after the P2-05 reviewer's independent check. The P1-12 reviewer stopped at the
  usage limit and was restarted fresh; P2-05 hybrid prep started. Files modified: up.md, training/{ConDiv.py,
  README.md}, py/martini_build_tables.py (comment), findings.md 3.14, plan.md, upside3.md, pass2/brief.md.
- 10-22 07:30-09:00, I/O and C API pass 2 (P1-12) reported and fixed (A175-A182: lost messages and a locked file
  after a refused in-process run, exact integer conversion on read, negative round / array string attribute /
  Brownian defaults refused, every run logs its starting round, dead dense logger, rama parameters as the model at
  rounding, up.md 1.7); reproducers equal the reviewer's after-logs, regression equal to its build bitwise; baselines
  w53, h51. Hybrid prep pass 2 (P2-05) reported and fixed (A183-A188: molecules wrapped whole, multi-chain molecule
  ids, the soft core from r = 0, a negative cache_buffer refused, run.py's POPE/POPG thickness, numbering gaps
  refused); regression bitwise. Its D1 confirms the side-chain table's bead placement independently; that, the table's
  averaging temperature, dead parameter content and the ion-count thickness asked of the user. `obj/` to A188
  (`*.bak_a182`). Files modified: src/{h5_support.h,h5_support.cpp,engine_c_library.cpp,main.cpp,martini_brownian.cpp,
  martini_hybrid.cpp,rama_map_pot.cpp,state_logger.h,martini_potential.cpp}, py/{run_upside,martini_prepare_system,
  martini_prepare_system_lib}.py, example/16.MARTINI/run.py, up.md 1.7, findings.md 2.4, 2.5, 3.14, plan.md,
  upside3.md, pass2/brief.md.
- 10-22 09:00-09:40, the user's decisions (all as recommended): rebuild the side-chain table from ideal chi-row
  conformers through martinize's mapping, averaged at the simulation temperature; delete the dead content (done,
  A189); the ion count from the P-P thickness (done, A190); then "Use DSSP" for the backbone types after I reported
  that it changes helical proteins too (glpG 84 of 210; done, A191, `dssp_bb/compare.log`). The table rebuild runs in
  an agent (`sc_rebuild/`); P2-07 ConDiv review started. Files modified: py/{martini_prepare_system,
  martini_prepare_system_lib}.py, example/16.MARTINI/readme.md, up.md 2.7, 6; deleted parameters/ff_2.0/,
  parameters/packing/, parameters/dryMARTINI/DDM.pdb.
- 10-22 10:40-12:30, the side-chain table rebuilt by an agent from ideal chi-row conformers through martinize's
  mapping, one free energy at 0.8647 (A192; beads against 1146 crystal side chains median 0.18 A, old 1.38; particle
  tables bitwise); installed with the old `martini.h5` backed up. The hybrid fixtures re-prepared end to end on A188
  + A189-A192 (`fixtures_hybrid_w4`): healthy (KE production 1.007, C-N max 1.82 A); two rotamer FD flags settled as
  noise by the line integral; baseline h52. ConDiv pass 2 (P2-07) reported: the trainer runs right on the repo engine
  (a two-protein step end to end, every group's derivative against FD); its D1 (the SARW kept the backbone desolvation
  term) checked against Kleinmann's `condiv2.py` on midway2 (a separate `bb_scale` never passed) and decided by the
  user for the SI (A193); A194 gate, A195 epoch length from `initialize`, A196 README. Files modified:
  py/martini_build_tables.py, parameters/ff_2.1/martini.h5, training/{ConDiv.py,ff31/ConDiv.py,README.md,
  train_chain.sbatch}, up.md 2.6, findings.md 3.14, plan.md, upside3.md, remote_jobs.md 0b, pass2/brief.md.
- 10-22 12:35-13:00, parameter glue pass 2 (P1-14) and hybrid analysis pass 2 (P2-06) reported; applied and verified
  together (bitwise on w53/h52; the wrapper runs a config without `/input/sequence`): A197 dead engine API, A198 the
  unread `/input/sequence`, A199 dead Adam hook and `_read_comp`, A200 the REMD FES tool (flat FES on hybrid data),
  A201 a REMD HDX analysis refused without `HDX_T_UP` (it ran at the bottom rung), A202 walker labels across
  position-only continuations, A203 dG errors NaN where none is computed. `obj/` to A203 (`*.bak_a196`). The user had
  me decide the glue's D1 and D6: torch out of the glue and trainers by closed-form gradients, the dead `rotscalar`
  block removed (an agent, `notorch/`); P1-17 build review started. Files modified: src/{engine_c_library.h,
  engine_c_library.cpp,deriv_engine.h,deriv_engine.cpp,main.cpp}, py/{upside_engine,rotamer_parameter_estimation}.py,
  tests/{regress,invariance,line_integral,node_fd}.py, training/{ConDiv.py,ff31/ConDiv.py}, example/00.AnalysisScripts/
  helpers/{analyze_remd,calc_hdx_ht}.py, example/16.MARTINI/{run_hdx_analysis.sh,run.py,readme.md}, up.md 1.1,
  findings.md 3.14, plan.md, upside3.md, remote_jobs.md 0b.
- 10-22 13:00-13:45, torch removal and build pass 2 (P1-17) reported; both checked and applied. A204: the glue and
  both trainers in numpy with closed-form gradients, `rotscalar` removed (31,420 to 31,334 coordinates); I reran the
  agent's torch comparison and import-order tests (logs identical: tables to 3.1e-16, gradients to 3.9e-16 of scale,
  one libomp in every order); its seeded two-protein ConDiv step matched torch's bitwise except `rot` at 1.4e-13.
  A205-A211 from the build review: cluster CMake OpenMP REQUIRED and dead lines, `engine_c_library.cpp` once, the
  Mac Eigen option honoured (eigen@3), the unused `omp.h` include, `source_arm`/`source_sh` without the empty
  PYTHONPATH entry (a stray `tables.py` in the run directory was imported) and the stale python@3.10 block,
  `training/env.sh` refusing outside a tree, `install_python_env.sh` listing what the code imports (jax kept: pymbar's
  backend). Built in `buildfix/home`: bitwise on w53 and h52; `obj/` to A211 (`*.bak_a204`); the local `source.sh`
  regenerated from the template (old copy `buildfix/old/source.sh`). The build review's U1 confirmed by `nm` on
  midway2: every cluster binary carries `set_fast_math` (flush-to-zero set on load), left for the first midway2
  build. Started P1-15 analysis and P1-16 examples reviewers. Files modified: py/rotamer_parameter_estimation.py,
  training/{ConDiv.py,ff31/ConDiv.py,README.md,env.sh}, src/{CMakeLists_M1.txt,CMakeLists_Other.txt,main.cpp},
  source_arm, source_sh, source.sh (untracked), install_python_env.sh, findings.md 3.14, plan.md, upside3.md,
  remote_jobs.md 0b, pass2/brief.md, pass2/build/report.md. Then A212, the dead `NONCOMPACT_SIGMOID` branch out of
  `src/vector_math.h` (bitwise on w53 and h52, `buildfix/results2.txt`); `obj/` to A212 (`*.bak_a211`).
- 10-22 13:47, the build review's B9 (user: "decide for me for that question"): `.devcontainer/` deleted (A213). It
  clones and runs GitHub master, so in Upside 3 a developer would edit the workspace and run Upside 2 without an
  error; nothing referred to it, and with no container runtime here a repaired one could not be tested. Copy in
  `scratchpad/upside3_audit/buildfix/old/devcontainer/`. Files: .devcontainer/ (deleted), findings.md 3.14,
  plan.md, upside3.md, progress.md.
- 10-22 13:50-17:30, analysis pass 2 (P1-15) reported after a usage-limit stop (both reviewers resumed at 17:05); all
  eight defects reproduced by me on an unpatched copy of example 00 against a regenerated chig REMD, then fixed in the
  repo (`scratchpad/upside3_audit/afix/run.sh`, `results.txt`): A214 `analyze_remd.py` (mixed ladders refused, the
  top bin closed, empty bins +inf; the same floor out of `calc_hdx_ht.py`'s H-bond plot), A215 the uptake curve grid,
  A216 PDB numbering and colours in 5.analyze, A217 the temperature fit's glpG tag exclusion made a setting with a
  refusal, A218 6.generate_hx_plots without glpG defaults, A219 `getkchem` pKD2O at T, A220 README and a dead MDTraj
  loader. A200's record corrected: its FES maximum 21.9 was 3 empty bins at the floor (real 4.16). Files modified:
  example/00.AnalysisScripts/{0.run_HXMS.py,4.calc_D_uptake.py,5.analyze_D_uptake.py,6.generate_hx_plots.py,
  analysis.sh,README.md,helpers/{analyze_remd,calc_hdx_ht,hxfunctions_clean}.py}, py/mdtraj_upside.py, findings.md
  3.14, plan.md, upside3.md, progress.md, pass2/brief.md, pass2/analysis/report.md.
- 10-22 17:30-17:55, examples pass 2 (P1-16) reported; pass 2 is complete for every unit. Reproduced and applied:
  A221 example 04's dG sentinel and unclipped log (a synthetic always-protected residue: 1000 and 20.087 kcal/mol, the
  latter set by round-off; NaN after), A222 the AFM tip clock resumed by an integer step count (on the A212 build a
  moving-tip continuation drifts 8.8 and 7.4 A from the continuous run at 50 + 50 and 300 + 300 time units; bitwise
  after, moving and still; old pulling configs refused; `efix/results.txt`), A223 the intensive-memory switch
  tri-state, A224 AFM and tension units in the help (E_up, not kT), A225 example readmes and comments (edited by me).
  Fixtures rebuilt by the writers: only the two AFM configs change (`step_initial`), so fixtures_w19 and its
  baseline w54 (equal to w53 bitwise); h52 bitwise. `obj/` to A225 (`*.bak_a220`). A130 and A172 records corrected.
  Files modified: src/coord_const.cpp, py/{upside_nodes,generate_restart_config,run_upside,advanced_config}.py,
  up.md, example/{01,02,05,06,07,08,09,12,13}/ readmes and scripts, example/04.HDX/4.calc_HDX.py, findings.md 3.14,
  plan.md, upside3.md, progress.md, remote_jobs.md 0b, pass2/brief.md, pass2/examples/report.md.
- 10-22 17:57-19:10, the audit's closing criteria (upside3.md 2.2 item 5), met. One driver (`closing/driver.sh`):
  invariance all PASS; node_fd on all 30 configs with every flagged rotamer direction line-integrated to noise (the
  line integral made to follow node_fd's tied direction and take a width, A230); sanitizers clean on a rerun that
  preloads ASan into Python.app (the first lane aborted everywhere: the python3.11 launcher drops DYLD_*). The
  registry, moved to the current fixtures, left 14 unmatched rows; a reviewer dispositioned all of them (no real
  disagreement) and found three latent defects. Applied: A226 the MARTINI box required (a box-less hybrid ran
  without periodic images at exit 0) and the non-periodic path removed (its release records move at float32
  rounding from fast-math code generation alone: built without -ffast-math the two sources are bitwise; node_fd on
  both hybrid stages clean), A227 dead writer outputs and the curvature changer's documented node, A228 RamaMap2's
  unused log_pot, A229 the writerless old membrane node (decided for the user, superseding A41), A230 the test tools.
  Baselines w55 on fixtures_w20 (bitwise w54) and h53 on fixtures_hybrid_w5; registry 44 of 50 matched, the rest
  read outside the engine. `obj/` to A230 (`*.bak_a225`). Also: `git pull` fast-forwarded to c42d570a (the other
  session's Boltz-2 note in remote_jobs.md), my uncommitted edits reapplied. Files modified:
  src/{martini_potential,rama_map_pot,membrane_potential,main}.cpp, py/{martini_prepare_system_lib,upside_config}.py,
  tests/{line_integral,registry,index_oracle}.py, up.md, findings.md 3.14, plan.md, upside3.md, progress.md,
  pass2/brief.md, pass2/registry/report.md.
- 10-22 19:20-22:30 (user: "keep ff2.1 useable for now"; then "git pull and merge. also fix D9 with (a). submit Peng 16
  validation protein with ff2.1 and glpG hybrid simulation with ff2.1, on midway2, as a validation for our entire
  debug"). `git pull`: fast-forward to 605901d2 (the other session's glpG ff3.0 HDX notes), my edits to the four shared
  files reapplied identically. D9 (an agent, then reproduced by me): A231, a required `sequence_id` exclusion key in every
  node that excludes by sequence distance; single-chain and hybrid records bitwise, only 1dfn's and 1dfn_cavity's seeded
  runs move, the oracle's 7 cross-chain mismatches gone, old configs refused; `obj/` to A231, baseline w56 on
  fixtures_w21, registry 44 of 50 as before. Validation: an Upside 3 tree `/beagle3/trsosnic/yinhan/upside3` (separate
  from the shared deployment, which week-long campaigns read), built on midway2 with GCC 10.1, its own venv with mdtraj
  1.9.9 from source (the prep's DSSP; glpG's DSSP equal to the Mac's); Peng's 32 arms and the glpG prep plus 4 REMD blocks
  submitted 22:26 after a smoke test (remote_jobs.md §1 line 3). Found: `example/16.MARTINI/run.py`'s Stage 2 cannot run
  (200,000 time units of 48 replicas in one 35:45 job), and the campaign's `run_remd.py` fails on Upside 3 and carries a
  NaN rollback guard; the validation uses a block driver without one. Files modified: src/{interaction_graph.h,
  bead_interaction.h,rotamer.cpp,environment.cpp,backbone_steric.cpp}, py/{upside_config,martini_prepare_system_lib}.py,
  tests/{index_oracle,registry}.py, up.md, findings.md 3.14, plan.md, upside3.md, progress.md, remote_jobs.md,
  pass2/brief.md; on the cluster the new tree and `upside3_validation/`.
- 10-22 22:30-23:05 (user: "as we changed the upside core, the job scripts should also be updated to support it").
  Inventory of the job scripts: examples 01-15 and ConDiv already run on Upside 3 (P1-16's runs; ConDiv writes fresh
  configs each step); `example/16.MARTINI/run.py` could not (A232) and `run_sim_1afo_full.sh` ran the pre-audit tree
  (A233). `run.py` rewritten to submit its stages as jobs, REMD in `--time-limit` blocks that queue the next; two
  1rkl blocks run locally (a 2-replica run found the swap-set index bug, and a failed first block that the counter must
  not count). The glpG validation relaunched through it (prep 49448466, block 1 49448467); the stand-in driver and its
  partial prep removed. The submission showed numpy's zero-subnormal warning: the cluster library's flush-to-zero
  reaches Python. Files modified: example/16.MARTINI/{run.py,readme.md,run_sim_1afo_full.sh}, findings.md 3.14, plan.md,
  upside3.md, progress.md, remote_jobs.md, pass2/brief.md.
- 10-22 23:10 (user: "yes. add it to claude.md"): CLAUDE.md now names the Upside 3 tree beside the shared deployment
  (`/beagle3/trsosnic/yinhan/upside3`, `env_upside3.sh`, why it exists, how it is synced and its venv built), in the
  Slurm environment list and wrapper example, and in the cluster compile section. Files modified: CLAUDE.md,
  progress.md.
- 10-21 16:10, the Boltz-2 HDX sampler: tag by step-scale test read (61231966, finished 09:03). `step_scale` 1.0
  causes the bond-length scatter (C-N sd 0.031 to 0.036 A with no tag, 0.054 to 0.055 with `md`; 0.008 to 0.009
  at 1.5 with either tag). The `md` tag widens the ensemble at both scales. Contacts under 2.5 A are mostly
  over-short salt bridges and H-bonds; one ubiquitin sample has Arg42 and Arg72 superimposed. The Phase 1 setting
  (sharpened 1.5 against the unscaled 1.0 with poor geometry) is open. Files modified: diffusion_dynamics.md
  (status, Phase 0, risks), remote_jobs.md, progress.md.
- 10-21 16:15-16:40, the Boltz-2 HDX sampler: protection flags on the tag by step-scale samples (user: "proceed
  with your plan"). Local, no GPU: `example/04.HDX` configs (ff_2.1) for ubiquitin and EHEE, Boltz backbones written
  as Upside trajectories, `py/get_protection_state.py` unchanged, from a scratch engine tree (the repo's `obj/` is
  older than its `py/`). Noise matched to step scale 1.0's bond spread moves per-residue protected fractions by at
  most 0.03; `md` 1.0 against 1.5 by at most 0.16, with fewer never-open amides at 1.0. Phase 1 will use `md` at
  `step_scale` 1.0. Files: diffusion_dynamics.md (status, Phase 0, Phase 1, risks), findings.md 10.14, progress.md;
  new `scratchpad/boltz_hdx/flags/`.
- 10-21 16:23, the Boltz-2 HDX sampler: Phase 1's first job submitted (user: "submit phase 1"), 61254153 on one
  midway3 V100, 1,000 `md` samples per protein at `step_scale` 1.0 (`phase1_hdx.sbatch`). Files: diffusion_dynamics.md,
  remote_jobs.md, progress.md.
- 10-22 03:15-03:30, the Boltz-2 HDX sampler: Phase 1's first job read (user: "check boltz job status", then
  "proceed next step"). 61254153 gave EHEE's 1,000 samples, matching Phase 0's 100; ubiquitin ran out of GPU memory
  because `diffusionv2.py` cut 1,000 samples into one chunk whatever `--max_parallel_samples` said, and Boltz exited 0.
  Fixed in our clone `$B/boltz` (`sample_ids.split(max_parallel_samples)`, uncommitted); wrappers now fail a skipped
  batch. Ubiquitin rerun with an EHEE chunking check submitted, 61279169. Files: remote_jobs.md, diffusion_dynamics.md
  (status, Phase 1), findings.md 10.14, progress.md; new `scratchpad/boltz_hdx/{phase1_ubq_rerun.sbatch,
  compare_chunks.py}`, `phase1_hdx.sbatch` changed.
- 10-22 12:05, the Boltz-2 HDX sampler: ubiquitin rerun read (61279169, finished 04:01). Chunking leaves the samples
  unchanged (EHEE, chunks of 10 against one batch: at most 0.001 A, PDB precision); ubiquitin's 1,000 `md` samples
  match Phase 0's 100, with 2 broken C-terminal C-N bonds (2.05, 3.71 A) at Phase 0's rate. Phase 1 now has both
  1,000-sample sets; flags next. Files: remote_jobs.md (both Phase 1 jobs to "Finished runs"), diffusion_dynamics.md
  (status, Phase 1), findings.md 10.14, progress.md.
- 10-22 12:10-12:50, the Boltz-2 HDX sampler: Phase 1 dG profiles (user: "proceed to next step"). Protection flags
  on both 1,000-sample sets (Phase 0's flags reproduced byte for byte), uniform-weight dG at the experiment's RT, and
  the comparison with Peng 2022 SI Tables S5/S7 (experiment) and S8/S10 (Upside FF2), transcribed and checked against
  the PDF. Every measured amide but one lies above Boltz's bound RT ln(1000) = 3.75-4.09 kcal/mol; Boltz opens 2 of
  39, Upside 11. Rank order and over-protection are not testable; the ubiquitin data are psWT (F45W/H68N), our
  samples wild type. Files: diffusion_dynamics.md (status, Phase 1), remote_jobs.md, findings.md 11b, progress.md;
  in `scratchpad/boltz_hdx/flags/`: `samples_to_traj.py`, `run_flags.sh` extended, new `dg_phase1.py`,
  `plot_dg_phase1.py`, `reference/`, `dg/`.
- 10-22 12:20-12:40, the Boltz-2 HDX sampler (user: "update global CLAUDE.md. also decide the two questions for
  me"). Global CLAUDE.md's plotting section now names the working venv (Python 3.10, matplotlib 3.10.0). Decided:
  a Ub psWT ensemble for a like-for-like ubiquitin comparison, and rank order tested on HEEH_rd4_0097 (measured dG
  2.2-4.3 kcal/mol, resolvable at 10^4 samples) rather than 10^6 samples of ubiquitin. `prep.py` extended (psWT and
  HEEH MSAs, references), `check_samples.py` subsamples pairwise RMSD above 1,000 (output unchanged below);
  61298446 submitted. Files: ~/.claude/CLAUDE.md, diffusion_dynamics.md (Phase 1), remote_jobs.md, findings.md 11b,
  progress.md; `scratchpad/boltz_hdx/{prep.py,check_samples.py,phase1_heeh_pswt.sbatch}`,
  `flags/reference/heeh_experiment.csv`, `flags/pdb/{ubqps,heeh}.pdb`.

- 10-22 18:20-19:20, glpG ff3.0 HDX (user: "check if we have a glpG simulation result on midway2 with ff3.0 after
  martini repaired"; then "generate delta G plot into ~/Downloads"). The four ff_3.0_gdepth chains on the 10-17
  repaired seeds exist, blocks 1-3 closed, block 4 running, logs clean. midway2's venvs lack mdtraj/pymbar/matplotlib,
  so HDX ran on midway3 caslake: block 2 and block 3 apart, 8 jobs, all 28/28 (one rerun after an interleaved
  log line in `extract_block.py`, fixed). Block 3 against block 2 disagrees by RMS 0.65-0.91 kcal/mol against
  jackknife sigma ~0.1, with helix-wide shifts and no ladder round trip (findings 5.6). Figures in
  `~/Downloads/glpG_HDX_ff3.0_blocks_20261022/`, redrawn in the project's `plot_ref_style.py` format after the user
  rejected a new layout, then as one figure per variant from blocks 2-3 pooled (4 more midway3 jobs, 28/28; lesson
  findings 10.36). Files modified: findings.md 5.6 and 10.36, plan.md Phase 18, remote_jobs.md;
  new `checks/hdx_ff30_placefix_20261022/` on the cluster.
- 10-23 01:00-01:10, the Boltz-2 HDX sampler: 61298446 read. Ub psWT's 1,000 samples match the wild type's
  statistics; HEEH's single 10^4-sample run ran out of GPU memory in Boltz's confidence step after 1 h 52 min, and
  the wrapper's new check reported FAIL. Cause traced to `confidencev2.py:121-143` (every sample's PAE/PDE logits
  concatenated on the GPU). Resubmitted as 10 seeds of 1,000 (61328353, `phase1_heeh_seeds.sbatch`). The Boltz dG
  figures will use `plot_ref_style.py` unchanged (findings 10.36). Files: remote_jobs.md, diffusion_dynamics.md,
  findings.md 10.14, progress.md; new `scratchpad/boltz_hdx/phase1_heeh_seeds.sbatch`.
