# Remote jobs on midway2/midway3: status and handbook

**Current state (2026-11-02 16:45 CST, the MacBook Pro; line 5 at 08:40 11-04 and ff3.2 at 17:50 11-03, the Mac Studio; user 10-31: "1. Deploy ff3.1, run validations (Peng, glpG, and
other tests we ran) 2. Fix the folding issue and continue working on ff3.2").** Owner: none since 11-04 (the Mac Studio's session, owner from 17:41 11-02, stopped after
its pass of 09:26 11-03; user 11-04: no watch until the user resumes on the Mac Studio); before it none since 21:55 10-27, when the MacBook Pro's session (owner from 12:46 10-26) handed off and quit
("Handoff"). Status passes in between, none the owner and none with a cron: a
Mac Studio session (08:06-13:45 10-29; also the user's "fix them all" and "run both"), the MacBook Pro's of 14:30
10-29, 17:10 10-30, 14:05 and 23:25 10-31, 14:07 11-01, 12:20 and 16:45 11-02, and the Mac Studio's of 01:30 and 08:50 11-01
and 00:15, 07:30, 08:30 and 17:20 11-02 (line 5 only) and 17:40 11-02 (ff3.2 only). **The MacBook Pro handed off 16:50 11-02** (user: "dump everything into md files. I
will start on another computer"): these records, its working-tree patch, its memory and its local-only scripts are in
the handoff bundle (addendum of 11-02 in its `README.txt`; scripts in `scratchpad/handoff_20261102/`). **midway2's fault, 21:19 11-01 to 09:43 11-02** (§1):
its Slurm answered srun and sbatch with "Unexpected message received" and its login node opened no session, while
jobs already running kept writing; a banner or a Duo login proved nothing until a command returned (08:31). At 09:43
the controller returned and stamped every job that ended meanwhile TIMEOUT at 09:43; their logs show their own ends.
It broke ff3.2's link 49737476 (insurance 49737477 took over 09:45) and 30 of ff3.1's 32 Peng chains, restarted 12:18
11-02 at the user's word (§1 line 2). The next computer's "continue monitoring" takes the jobs ("Resume here"). The Mac
Studio's uncommitted records are in the handoff bundle (`records/`, `repo_working_tree_mac_studio_20261029.patch`
against `3362e2a`) and its session tools in the bundle's `scratchpad/handoff_20261029/` (README). **ff3.1 means
`ff31_nativefix` at step 76 since 16:49 10-31** (released as `ff_3.1`, under validation; plan.md Phase 37); the 10-26
merge is retired (`.bak_pre_nativefix_20261031`), as were the older sets before it (`.bak_pre_merge_20261026`).
**ff3.2 (`ff32_nativefix`) is finished and stops: step 95 reached 03:10 10-31, its gate there CONVERGED (dhb p 0.027,
cut 0.005), and its epoch-4 panel `n32_04` folds 0.447 against released ff2.1's 0.615** (-0.169 [-0.215, -0.126]), so
by plan.md Phase 33's rule it gets no full validation; what the loss is: findings 1.37, plan.md Phase 33 item 9 (§1
line 2). Its fix (plan.md Phase 38): with the DSE target off epoch 0 folds 0.575 (with it 0.449, released ff2.1 0.615), so
`ff32_dse0` was ff3.2 from 11-01, every term trained (checked: `params_check.txt`), each epoch end queuing its own
panel, unfolded shift, Rg check and next epoch, stopping at step 76 (FF2's budget, no convergence gate; user 11-02,
plan.md "End of training and release"). Epoch 1 (`z32_01`, read 17:40 11-02) folds 0.578, -0.038 [-0.072, -0.005]
under released ff2.1, and the unfolded chain compacts (Rg at T 1.1 0.939 of ff2.1's files', epoch 2 0.923). The
remaining loss's fix, the native-basin restraint, ran as the one-epoch arm `ff32_basin` and was read 03:05 11-03 by
the rule fixed before its results (plan.md Phase 38 item 10): it meets both parts, `w32_00` folds 0.630 against
`z32_00`'s 0.575 with mean TM level. **Since 11-03 ff3.2 is `ff32_basin` trained on** (user: "Ok. Proceed"; Phase
38 item 11), tags `w32_EE`, from step 19 to 76: link 49757705 R from 11:28 11-03, step 31 at 17:33 (30 min a
step; H-bond margin +0.098 to +0.106 since step 19, under +0.10 at steps 24-26 and 31, user notified). `ff32_dse0`
trains on to 76 as the comparison: step 69 at 17:37; its `z32_02`, read 17:50 11-03, folds 0.599, level with
released ff2.1 (-0.016 [-0.059, +0.026]) (§1 line 2). Running on midway2
besides: ff3.1's validation (panels and TM read, level with released ff2.1; 32 Peng arms and their collector; its
glpG block 1 on Upside 3 done, TM4 DSSP 0.74 at T 0.70; §1 line 2). Upside 3's validation has no job: glpG ff3.0 done (3 of 3 blocks, TM4 held),
every Peng arm scored (de novo 16 of 16 below Upside 2's, -0.030, p 0.004), the 8 seed replicates scored and read 16:50 11-02 (top7's gap beyond its seed spread, alpha3d's within). Released ff2.1's Peng is DONE, 32 of 32 (§1 line 1).
Upside 3 on the GPU (line 5): soluble proteins on ff2.1 and ff3.1 are validated and closed (11-03, user: "close
what is validated"; plan.md Phase 35 Resume here). The hybrid and the mdCATH panel run as plan.md Phase 35 item 12
steps 3-5 (user 11-03: "proceed to next validations"; build g9a, installed 21:54 11-02): glpG V2 by seeds passed
(22:10 11-03, Bonferroni 0.94); glpG's TM4 block on a V100 agrees with the CPU's (23:40 11-03); the mdCATH panel's
CPU reference (49777036-43) is complete; the probe (A100, after the V100 could not hold batch 4's 20 GB) put the
host work at 19%, so it now runs in parallel (build g9b, installed 00:36 11-04, bitwise g9a's on a V100); the GPU
panel's 8 jobs completed by 08:29 11-04 and pass the CPU reference's checks, the folded fraction agreeing; the
TM-score comparison, step 5's deciding test, passes on both force fields (49803361-62, read 10:50 11-04; §1 line 5). Boltz-2's E1 particle-count test at low noise (user 11-02: "Particle scaling, full"; §1 line 4), after
61805263 did not remove the collapse at 256 particles: 62073874-78 COMPLETED by 15:33; the two 16,384 arms 62073879-80
ran to their end but lost their results to the quota below, and the readout 62073881 FAILED. Read on the Mac Studio
11-03 from the finished runs: removal fails at 1,024 and 4,096, and E1 is closed (user 11-03; diffusion_dynamics.md
Phase 2 E1, Phase 3 next). `/beagle3/trsosnic` passed
its group inode soft quota and the grace expired: no file could be created there 19:46-21:51 11-02, until the reclaim
of 21:51 (disk section). Each line of work has a collector, a
dependency chain or an item under "Owed by hand". `WATCH_STATUS.md` on the cluster is newer than this file whenever a
pass ran after it. Job ids, logs and data paths are in §1.

The user's four lines of work (10-21), and what each has on the clusters:
1. **Deterministic vs probabilistic side chain** (plan.md Phases 20-27; findings 1.21, 1.37). Both trainings are
   finished and converged; the deterministic side chain folds less on every held-out test, and both retrained force
   fields fold less than released ff2.1. Released ff2.1 on the same three tests decides whether ff2.1 itself stays: the
   44 domains and mdCATH are done (released ff2.1 ahead of both on each); Peng is DONE, 32 of 32 scored. The email to Tobin and John waits on these numbers.
2. **Ramachandran map training** (plan.md Phases 14-17, 21, 28-30, 32; findings 1.34, 1.37). ff3.1 (the merge) is
   level with released ff2.1 on both soluble panels and above the control on both; Peng 13 of 32 scored (19 complete). The ff3.1
   training debug (Phase 32: does the trainer recover ff3.0's glycine pair) runs two arms to step 76. Finished,
   awaiting review: the ff3.1 and ff3.2 reruns on the fixed engine, whose maps barely moved and which, as every run
   trained from ff2.1 with all terms free, fold 0.15-0.24 less than ff2.1. The 8 Upside 2 glpG REMD chains on ff3.0 and
   the old ff3.1 finished all five blocks 10-24; ff3.0's TM4 checks are owed.
3. **Upside 3 design, debug stage** (upside3.md; plan.md Phase 24 and its "Resume here"). The audit's closing criteria
   are met locally (ledger A1-A231). Its validation runs on midway2 from the Upside 3 tree
   `/beagle3/trsosnic/yinhan/upside3` since 22:26 10-22 (§1 line 3), all finished: Peng's 16 proteins with ff2.1 (32
   of 32 scored; native level with Upside 2, de novo -0.030 [-0.049, -0.014] over 16), the Peng seed replicates
   (unscored), and the glpG hybrid with ff2.1 (4 blocks), the merge's ff3.1 (1) and ff3.0 (3).
   The TM4 test is answered (10-27, findings 3.15): TM4's loss under ff2.1 is ff2.1's own. At T 0.70 ff3.0 holds TM4
   (DSSP 0.85, 0.82, 0.85 over blocks 1-3), the merge's ff3.1 does not (0.42), and ff2.1 loses it (0.44 to 0.05).
   The released ff3.1 reads 0.74 in its one block, falling within it (halves 0.78 / 0.70; line 2).
4. **Boltz-2 HDX sampler** (diffusion_dynamics.md). Phases 0 and 1, the BioEmu2 side test and the controlled 10^5
   tests are finished and read out (§1 "Finished runs"). Phase 2's E1 pilot passed 10-28; E1's tuning (user 10-28 "run
   E1") shows the Feynman-Kac estimate failing above kappa 1, and more particles (1,024, 4,096; read 11-03) did not
   fix it, so E1 is closed (user 11-03). Phase 3 replaces the steered sampler with replica exchange whose move is
   Boltz-2's own renoise-denoise; its toy test T0 passed on the Mac 11-03, and T1, the same move on Boltz-2 (Ub psWT),
   failed its parity check as 62199777 on a code defect (a 5-step schedule), passed it as 62262317 and failed after it
   on a second (a strided index); both fixed 11-04, it reruns as 62263853 (§1 line 4).

Open and not in the user's list:
- ff3.0's four Upside 2 glpG chains: TM4 over blocks 1-5 read 14:56 11-04 (findings 1.35), then their trajectories
  deleted (user 11-04: ff3.0 is transitional; "Disk"). ff3.0's HDX on blocks 2-3 is done (findings 5.6); ff3.1's was dropped by the user 10-23, and the user closed the TM4 HDX question 22:15 10-23
  (plan.md Phase 16). The fraying turn, Val145 and the weak protein protection are read on Upside 3's glpG runs.
- Large-protein de novo folding (plan.md Phase 22) is done; its proposal (compact misfolds as negative examples,
  training past 150 residues, replica exchange for the large-protein benchmarks) awaits the user.
- The BP validation deck for Tobin (§0c) is not sent.
- Awaiting the user: the large directories on `/project/trsosnic` (§1 "Disk") and the deferred HDX post-processing
  fixes (findings 5.5).

The ConDiv runs, all on midway2 broadwl at dt 0.009 (`$P/training/<run>`, `$P` =
`/project/trsosnic/yinhan/upside2-md-mdw2`; "trained" is the port's rates, hb 0.01, rot 0.0125):

| run (tags) | what it is | state |
|---|---|---|
| ff21_nobp_dt009 | ff2.1 retrained with one side-chain state per residue type (the deterministic side chain) | converged at step 76 10-19; step 37 released as ff_2.1_nobp |
| ff21_ctrl_dt009 (`c9_EE`) | ff2.1 retrained the same way with BP (the probabilistic control) | converged at step 76 10-20; step 37 released as ff_2.1_ctrl38 |
| ff21_bpmax_dt009 (`mdc_bpmax_EE`) | the same with BP's most probable state per residue (the most-probable-state side chain); in `/project/trsosnic/yinhan/upside2-md-bpmax`, not `$P` | running from 12:42 11-03: link 49766810 to 38 on the pre-fix engine, then the fixed engine to 76 (§1 line 1) |
| ff30_gdepth_dt009 (`d9_EE`) | ff3.0: NDRD plus one trained glycine depth pair, every term trained | released as ff_3.0_gdepth 10-09 |
| ff31_depth_data (`fd_EE`) | the old ff3.1: 1,442 per-map basin depths, prior 0.2, every term trained | released as ff_3.1_depth_data 10-13; retired 10-26 (`.bak_pre_merge_20261026`) |
| ff31_depth_all (`fa_EE`) | ff3.1 with all 5,040 depths, prior 0.3 | finished 10-13, not validated |
| ff32_depth_data (`f32_EE`) | ff3.1's depths started from ff3.0 | stopped by the user in step 69 10-18 |
| ff31_depth_data_angfix (`a31_EE`), ff32_depth_data_angfix (`a32_EE`) | the ff3.1 and ff3.2 runs on the fixed engine | finished step 76 10-21; `a31_EE`'s panel trajectories deleted 11-04 (folded and Q tables kept, "Disk"); `a32_EE`'s deleted 11-04 at the user's word (tables kept) |
| ff311_depth_pooled, ff321_depth_pooled | ff3.1.1 / ff3.2.1 with every term trained | cancelled 13:38 10-21, 3 min in (user's revision); kept as `<run>.cancelled_20261021` |
| ff311_depths_only (`p311_EE`, ff3.1.1), ff321_depths_only (`p321_EE`, ff3.2.1) | ff2.1 with only its depths trained, pooled prior, from NDRD's and ff3.0's map | finished step 114 10-25; merged into ff3.1 and released as `ff_3.1` 10-26 (line 2) |
| ff311_glypair, ff311_glypair_mb2 | the ff3.1 training debug: ff3.1.1 with glycine as ff3.0's shared pair, two minibatch orders | `ff311_glypair` running from 23:05 10-27, `_mb2` queued; target 76 (line 2) |

Dropped by the user 10-11 and 10-12: the BioEmu-map line, ff30_gdepth_si, the H-bond-frozen arms, polygly and the
local TM4 glycine-flip test (what was deleted and what kept: §1 "Finished runs"). midway2 has no GPU that can run
Boltz-2.

CPU work runs on midway2 broadwl unless the user names midway3 caslake; midway3's login node may be used to move or
read files on `/project` and `/beagle3`. Since 10-07 10:48 midway2-login2's `/project` hangs; route `/project` work
through login1 (§8). §1 is live job state; §0, §2 and §5-§8 are the handbook.

---

## Resume here, from any computer (2026-10-27)

### "continue jobs" or "continue monitoring": what a new session does, in order

Claude's memory does not travel between computers; everything needed is in this file, plan.md and findings.md (pulled
from git) and on the cluster. If this computer's repo lacks the last owner's edits (the user has not committed them
yet), read the handoff bundle's `records/` and working-tree patches ("Where things are") and tell the user.
1. **Read** this file's header and §1 (live state; each watch pass leaves its readouts there), then:
   - plan.md Phases 24 (Upside 3, with its "Resume here"), 27 (released ff2.1), 30 (the merged ff3.1), 32 (the ff3.1
     training debug) and Known Errors;
   - findings.md 1.34 (the depth maps and the merge), 1.37 (the retraining loss), 3.14 (the audit ledger, A234), 3.15
     (Upside 3 against Upside 2: glpG and Peng), and the rules in §10;
   - diffusion_dynamics.md Phase 1, section 7 (the controlled 10^5 tests) and Phase 2 E1.
2. **Connect** to both clusters (§0; CLAUDE.md "Cluster SSH"). Check each socket first. If
   `scratchpad/mdw2_master.exp`, `mdw3_master.exp` and `~/.bin/ssh_mdw3` exist on this computer, run each expect script
   once (one Duo push each); otherwise ask the user to open the masters. The expect scripts carry the password and are
   never copied to the cluster.
3. **Take stock, read only.** One pass of the watch prompt ("The watch"), with `sacct -S <WATCH_STATUS.md's last pass>`
   on midway2 (the last pass was the MacBook Pro's, 14:30 10-29; the Mac Studio's scripts `st_1229.sh` and
   `st_peng4.sh` are in the bundle's `scratchpad/handoff_20261029/`). Tell the user what started, finished or failed
   since.
4. **Owed by hand** (`WATCH_STATUS.md` is newer if the two differ; details in §1):
   - **Boltz Phase 3's T1, 62263853** (§1 line 4; the second rerun, after 62199777 and 62262317): when it ends, pull `t1.npz` and read it on the Mac with
     `scratchpad/boltz_hdx/addon/t1_read.py` by the rule in diffusion_dynamics.md Phase 3 T1; report to the user.
     Boltz E1's lost 16,384 arms are not rerun (E1 closed 11-03). The `/beagle3` quota was reclaimed 21:51 11-02 (disk section), and every job that started while
     it blocked was swept: ff3.1's collector 49748039 and its two scorings ran 21:45-21:48 and wrote to `/project`
     (scores present), ff3.2's 49748605 writes only to `/project`; the A100 probe reran as 62120973 (read).
   - Upside 3 glpG (§1 line 3): every block read; ff2.1 (4 of 4) and ff3.0 (3 of 3, block 3 read 00:09 11-02) are
     done. Do not run `glpg_ff30/run.py`, which would queue another block.
   - Upside 3's Peng arms: all 32 scored against Upside 2, `compare_u3_u2.py` rerun 00:12 11-02 over the 16 de novo
     pairs (§1 line 3).
   - The Peng seed replicates: scored (`sc_seed_*` 49743788-49743795 COMPLETED by 12:41 11-02) and read 16:50 11-02
     with `peng_seeds/score/seed_spread.py` (`373e6733`): Upside 3 - Upside 2 over seeds 1, 101, 202 is top7 -0.087 +-
     0.031, alpha3d -0.086 +- 0.064 (plan.md Phase 24 item 2b 3, findings 3.15). Nothing owed.
   - ff3.2 on the fixed native target is closed (all its panels read; stops by plan.md Phase 33's rule). The round
     readout for the other runs is the bundle's `scratchpad/handoff_20261030/round4_read.sh` with the epoch changed,
     and `handoff_20261029/round2/gly_LR_runs.py` on the round's `rama_round_NN.dat`, run on the Mac (its `sys.path`
     line names the Mac's `scratchpad/ff31_compare_r6/`, which is in the bundle).
   - ff3.1 (`ff31_nativefix` step 76, released 16:49 10-31; plan.md Phase 37): its validation (§1 line 2). Both
     panels and glpG's preparation read 23:40 10-31, the mdCATH TM files 01:30 11-01. Owed as each ends: each Peng arm's
     final frames at `[bench] COMPLETE` (the collector scores them). glpG block 1 is done and read (07:30 11-02).
     ff3.1 is transitional (user 11-04): the scored arms' replicas, both panels' trajectories and glpG's REMD replicas
     were deleted 11-04 after their readouts were archived ("Disk"); every scored arm's replicas may go as it is
     scored.
   - The all-term loss's fix (plan.md Phase 38): `ff32_dse0` and `ff32_sarwsi` queue their own panels (`z32_00`,
     `s32_00`) and unfolded-shift checks, read 14:10 11-01 (§1 line 2). Since 11-03 ff3.2 is `ff32_basin` trained on
     (user; plan.md Phase 38 item 11), and `ff32_dse0` trains on to 76 as its comparison. For each epoch end of both
     read the panel (`w32_EE` against `ff21_released`, `w32_00` and `z32_EE`; `z32_EE` against `ff21_released` and
     `z32_00`), `peng_full_<tag>.txt` (`checks/ff32_fix_20261031/`), the Rg check (`checks/ff32_rest_20261102/
     rg_check_ff32_basin.txt`, `checks/ff32_fix_20261031/rg_check_ff32_dse0.txt`), for `w32_EE` also
     `decomp2_w32_EE.txt` and `crystal_view_w32_EE.txt` (`checks/ff32_rest_20261102/`, against `w32_00`'s -0.40 and
     -1.76), and the round's glycine depths; at 76 the epoch choice of plan.md Phase 38 item 7. Read so far: `z32_00`
     to `z32_03`, `w32_00` to `w32_02` (`w32_01`, `w32_02` read 15:30 11-04, plan.md Phase 38 item 11); `w32_03` and
     the epoch choice are read by the detached job 49816524 into `checks/ff32_epoch3_20261104/summary.txt` (about
     20:30 11-04): read it, then validation of the chosen epoch or, if below released ff2.1, the fixes of item 11.
     ff3.2's restart points are archived read-only in `checks/ff32_checkpoints_20261104/` (README). The panel read,
     from `ff3_selection/` with `training/ff30_gly/env.sh`: `python3 -I
     /project/trsosnic/yinhan/checks/converge_20261020/paired_npz.py domains aa runs <tag> <reference tag>`.
     The glycine depths: copy the round file
     (`$P/training/<ff32_basin|ff32_dse0>/run_output/rama_round_NN.dat`) to a Mac and run the bundle's
     `scratchpad/handoff_20261102/gly/gly_depths.py "ff3.2 rNN=<file>"` and `gly_LR_runs.py` the same way, from the
     repo `.venv` (they read ff2.1's, ff3.0's and ff3.1's libraries from `/Users/yinhan/Documents/upside2-md/parameters`;
     round 1 reproduces -0.877 / -0.640).
   - The ff3.1 training debug: both arms at 76 (`_mb2` 11:44 10-30, `ff311_glypair` 03:35 10-31), round 4's glycine
     pairs read (§1 line 2), not yet judged by plan.md Phase 32's rule.
   - Upside 3 on the GPU (line 5): item 12 steps 3-5 passed (step 5 read 10:50 11-04); step 6, the report and
     `upside3.md` 10, is owed (plan.md Phase 35).
   - Read by the collectors and reported when they appear: ff3.1's `peng_ff31.txt` (new ff3.1's collector 49692241,
     `checks/ff31_nativefix_release_20261031/`; the merge's stopped at 25 of 32).
     Released ff2.1's `peng_ff21.txt` is DONE (32 of 32, §1 line 1); its summary lines are NaN, since the retired
     `ff_3.1_depth_data` column scores nothing and the summary takes only proteins scored in every column.
5. **Find the owner** ("Handoff"): `cat` the cluster's `checks/r4_epochs/tm4_local/WATCH_STATUS.md`. A move between
   computers (user, 10-26: "be prepared to be shut down and continue this session on another computer"; done Mac
   Studio to MacBook Pro 12:46 10-26, and the MacBook Pro's handoff 21:55 10-27, owner none): the user's "continue
   monitoring" on the next computer is the take-over. If `WATCH_STATUS.md` names an owner whose pass is under 2 h old,
   ask the user to confirm that session is closed before starting a watch.
6. **Take the jobs:** start the watch ("Starting the watch"; no submitter is needed), then rewrite `WATCH_STATUS.md`
   with this computer as owner.
7. **Open questions for the user.** Nothing is blocking. Awaiting the user's word, not to be done without it:
   - **Boltz Phase 2:** requesting the Ferrari et al. data (the user only) matters only if E1 passes; E1 proper is held
     on the tuning readout; the diagnosis 61755382 could not be read (no reproduction at the same seed), and the
     reproducibility job 61786565 showed exact repeats are impossible on these GPUs (spread of `log_z` at kappa 1:
     s.d. 1.6 over 8 runs). Read statistically (user 10-29), the diagnosis points to a defect in the weights at low
     power; the toy test with a known normalizing constant (12:55 10-29) then found `sampler.py`'s weights and
     `log_z` right, pointing to the estimator's variance with Boltz-2's network. Resampling at low noise only
     (user 10-29: "run both") passed the toy but not the GPU (61805263, 14:40 10-29: ln Z-hat -15.05 and -20.26
     against the -11.1 line), so the limit does not rescue the estimator at kappa 2 with 256 particles; whether
     more particles do is the test running since 11-02 (§1 line 4), whose outcome decides the next question.
   - The force field for the Upside 3 hybrid HDX (findings 3.15), and the DDM micelle route (plan.md Phase 31).
   - The large-protein de novo proposal (plan.md Phase 22 review).
   - The large directories still on `/project/trsosnic` (§1 "Disk").
   - The deferred HDX post-processing fixes (plan.md Known Errors; findings 5.5).
   - The BP validation deck for Tobin (§0c), unsent.
   - **Upside 3 training terms** (`upside3.md` 8.2; plan.md Phase 23 item 5): H1 is done (findings 5.7); H2-H5 and the
     baselines of 8.3-8.4 are planned and run later (user 10-28: "record the plannings. We will run it later"). H2 is ready once ff3.1's
     ubiquitin benchmark arms are COMPLETE, but it also waits for the user's word.

### Handoff: two computers, one owner (2026-10-26)

**Owner now: none, and no watch runs** (user 11-04: "Mac Studio stopped, but do not monitor. I will resume on Mac
studio later today"). The Mac Studio owned the jobs from 17:41 11-02 (session `upside2-md-a0`, watch cron `1b3ab805`;
user: "plan that, while still monitoring") to its last pass 09:26 11-03; the MacBook Pro's status pass and fixes of
10:40-11:45 11-04 were not a watch. A second Mac Studio session, `upside2-md-fe`, reads ff3.2's jobs by hand and owns their
record entries (line 2's ff3.2 bullets and rows, plan.md Phase 38); it runs no watch. Before, no owner from 21:55 10-27
(user: "Record this session's progress and I will quit and work on another computer"). The MacBook Pro owned the jobs 12:46 10-26 to its last pass 21:32 10-27 and deleted its watch cron at the
handoff; the Mac Studio owned them 23:55 10-24 to 10:29 10-26, and its non-owner session of 08:06-13:45 10-29 handed
off at the user's move to another computer (no cron to delete). The next computer's "continue monitoring" takes them
("Resume here" step 6).

**One owner at a time.** Another computer may open a session at any time; it must know what the owner is doing and
not do it a second time, since two submitter crons would race on the midway3 ledgers. After "Resume here" steps 1-4, a
session `cat`s `WATCH_STATUS.md` over the midway2 socket and then:
- **"Owner: none": take the jobs** ("Resume here" step 6). Its shell and paths may differ: `~/watch_probe/` and
  `checks/` on the cluster hold every tool the watch uses.
- **An owner on another computer whose last pass is less than 2 h old: stand by.** Do not start a watch or submit
  anything; reporting status read only and answering the user are fine. Take over only on the user's word ("take
  over"). If its last pass is more than 2 h old, tell the user and ask whether that session has stopped; do not assume
  it. Taking over: if the old owner's session is still open, the user has it run its stop steps first; then open the
  midway2 socket (§0), start the watch ("Starting the watch") and rewrite `WATCH_STATUS.md` with the new owner.

**Owner's session restart.** The crons are session-only and end with their session; the computer stays the owner, and
the midway2 socket master is its own process, which survives unless the computer restarts. The new session reads §1
and `WATCH_STATUS.md` and checks the socket (the expect script once if it is dead), runs one watch pass at once for the
time since the last pass, starts the watch again, and rewrites `WATCH_STATUS.md` (last pass; any restart note removed).

**The owner's stop steps** (when the user says "hand off", "stop" or closes the session): finish the pass in hand,
update §1 and `WATCH_STATUS.md` (owner: none), then CronDelete the watch and any submitter.

**Live status on the cluster:** `/project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md`. The owner
rewrites it at every watch pass (prompt step 4; "Round-4 watch" step 5): owner, time of the last pass, the newest step
of each run, and what is owed. This file reaches another computer only through the user's commits, so it can be hours
old there; `WATCH_STATUS.md` is the current word.

### Where things are

**What travels with git and what does not.** Tracked: `plan.md`, this file, `findings.md`, `progress.md`, `training/`,
`src/`, `py/`, `up.md` and `parameters/common/rama31.dat` (the BioEmu-fitted library). Commit here, pull there. Not in
git (`.gitignore` has `*scratchpad*`): the `scratchpad/` tools, including `mdw2_master.exp`, `mdw3_master.exp` and the
local TM4 test (the MacBook Pro's scratchpad lost its run data 10-15 at the user's "cleanup", by the findings 10.25
rule: scripts, inputs and small results kept, manifest `scratchpad/CLEANUP_20261015_manifest.txt`). Their
authoritative copies are on the cluster:
- `/project/trsosnic/yinhan/checks/handoff_20261026/` (`README.txt`), the handoff bundle: `records/` (this file,
  plan.md, findings.md, progress.md and diffusion_dynamics.md, refreshed by the watch after each record change); the
  Mac Studio's 10-26 copies of the same, its uncommitted working tree as `repo_working_tree.patch` (against `b1e70d0`)
  and the MacBook Pro's as `repo_working_tree_macbook_20261027.patch` (against `3005a3fd`); the cron prompts (the 10-27
  `cron_watch.txt`); the git-ignored readout tools (`scratchpad/`, scripts and small inputs; in `scratchpad/boltz_hdx/`
  the Phase 2 add-on `addon/` without the MacBook's CPU sandbox, and `flags/urea_classes.py` with its log; in
  `scratchpad/hdx_term_20261028/` H1's scripts and readout, `upside3.md` 8.2); Claude's memory directories (`memory/`,
  `memory_macbook_pro_20261027/`).
- `/project/trsosnic/yinhan/ff3_selection/`: the panel tools (`panel.py`, `panel.sbatch`, `submit_new.sh`,
  `slurm.args`).
- Under `/project/trsosnic/yinhan/checks/`: `r4_epochs/<tag>` (extracted checkpoints); `r4_epochs/tm4_local`
  (`WATCH_STATUS.md` and the retired local TM4 test: README, `scripts/`, `seed/`, `runs/` pre-fix, `runs_cov/`);
  `cov_fix_20261007` (the tested patch-script fix); `fz_init_20261007` and `dt009_init_20261006` (the initial force
  fields, shown byte-identical); `r4val_20261005` (the automatic-validation tests); `broken_replica_20261006` (the
  comparisons behind findings 1.21); `gly_bioemu_map` (the BioEmu library, its fit and the push probe, README);
  `hbg_deploy_20261002` (the worked binary deploy); `ff31_deploy_20261009` (the old ff3.1's trainer files as deployed,
  the patched `check_step.py` with a branch for per-map offsets, and the test; `training/ff31` and `ff31-gly` are
  September's retired FF1-form glycine retrains); `tm1_hbmem_20261009` (the 10-09 glpG TM1 work, findings 1.27 and
  4.6-4.9, plan.md Phase 13, with a README; its `mac_studio_20261010/` holds the follow-up of findings 4.10-4.14);
  `sc_table_20261010` (the converged side-chain/lipid table's deployment, findings 4.15: staged martini.h5, builder,
  `patch_glpg.py`, `deploy.sh`, `repatch.sh` and `check_repatch.py` for the 8 re-patched glpG seeds; every replaced file
  kept as `.bak_pre_sc600x12_20261010`). `$P/py/martini_build_tables.py` keeps the old sampling (campaign rule): never
  rebuild martini.h5 in `$P`.

Claude's memory is per computer; the rules that matter are in this file and findings.md §10. The rotamer-BP
validation figure and deck are in §0c.

**The cluster's `training/` is not the repo's layout, and must stay so** (rule in §1b). `$P/training` keeps the
pre-10-02 layout: `ConDiv.py`, `train_chain.sbatch`, `check_step.py`, `extract_ff.py`, `convergence_gate.py`,
`gate_or_continue.sh`, `validate_ff.sh` and `patch_glpg.py`, with each run's own trainer copy in `<run>/trainer/` where
it differs. The runs' `after_training.sbatch` calls `gate_or_continue.sh`, which runs `validate_ff.sh`. Do not sync the
repo's `training/` over it. The repo's `ConDiv.py` and `train_chain.sbatch` carry the same dt 0.009 and 54-step links.

### The watch

**One watch at a time.** A watch is a session cron: it dies with its Claude session and expires after 7 days. No
watch runs since the Mac Studio's session stopped (its cron `1b3ab805`, session `upside2-md-a0`, hourly at :23 from
17:41 11-02, last pass 09:26 11-03; user 11-04: "do not monitor"); its prompt is "Hourly watch prompt" below with
ff3.2's readouts left to `upside2-md-fe`. None ran 21:55 10-27 to 17:41 11-02 (the MacBook Pro's cron `4ff7a1bc`, 12:46 10-26 to its last pass 21:32 10-27, was deleted at
the handoff; user 10-24: "continue monitoring midway2 and midway3"). A watch is read only: it submits, resubmits and
changes nothing, never reopens a socket, and reports to the user; after a record change it copies the five records
into the handoff bundle's `records/`, so the bundle stays current if the owner shuts down first. Each pass records what
changed in §1, so the newest §1 tells a new watch where the last one stopped. The collectors do the routine checks of
the training and validation runs (§1). No submitter is needed: the last one (cron `0fd74c14`) submitted the controlled
10^5 tests on midway3 and deleted itself 00:39 10-27 when neither midway3 ledger had a `pending` or `hold` line.

**Starting the watch on a computer** (set up as in "Handoff"): CronCreate hourly at `23 * * * *` with the watch prompt
below, its owner in step 4 replaced with this computer. The prompt is also the handoff bundle's `cron_watch.txt`
(10-27, md5 `9969518a`; the Mac Studio's 10-26 prompt is `cron_watch.txt.bak_mac_studio_20261026`), except step 3,
which there still reads the finished E1 pilot 61663101; use the step 3 below. Its midway2 half
is one read-only probe, `~/watch_probe/watch_mdw2.sh` on midway2 (md5 `cf834d4b` since 10-31, when its Peng table and TM
files moved to the new ff3.1's `checks/ff31_nativefix_release_20261031/`; `0098c032` 10-26 to 10-31, kept as
`.bak_pre_nativefix_20261031`; older
`.bak_pre_ff31_20261026`): the queue by job name (Peng arms of every force field grouped), every job that left it in
the last 80 minutes (each Peng chunk with its `[bench] resubmitted` or `COMPLETE` line), two of the three Upside 3 glpG
runs (ff2.1 in `glpg/`, ff3.1 in `glpg_ff31/`: REMD error lines and `block_count`; ff3.0's `glpg_ff30/` is read by the
prompt), the Peng tables `peng_ff21.txt` and `peng_ff31.txt`, ff3.1's TM files, and `/project` space and inodes. Other
probes there: `~/watch_probe/stock.sh` and `~/watch_probe/steps.sh <run>:<step> ...` (a training run's "Round-4 watch"
steps 2-3), `~/watch_probe/logs.sh` and `final_frames.py` (steps 2 and 2c).
**The probe's 80 minutes cover one pass.** When `WATCH_STATUS.md`'s last pass is older (a cron recreated, the session
busy: the cron fires only when it is idle), the pass also runs `sacct -X -u yinhanw -S <last pass>` on midway2 and
reads every job that left the queue since. On 10-25 three skipped passes (10:51-12:51) hid ff3.2.1's NODE_FAIL of
09:56 until the 13:39 pass found the queue counts one arm short.

**Hourly watch prompt** (10-27, for the next owner):
"Hourly watch of midway2 and midway3 (user 10-24 23:45: "continue monitoring midway2 and midway3"; remote_jobs.md "The watch"). Job ids, logs and next actions are in remote_jobs.md §1 and its job table. READ ONLY: never submit, cancel, resubmit or change a job, never edit a cluster file other than WATCH_STATUS.md and the handoff bundle's records/, and never run mdw2_master.exp / mdw3_master.exp (each launch sends a Duo push). Always pass -o BatchMode=yes and check a socket before each call through it. Run cluster Python from our own directories with python3 -I, never from /tmp (findings 10.41).
1. Sockets: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu` and the same with ~/.ssh/cm-mdw3.sock and midway3.rcc.uchicago.edu. If one is dead, tell the user, skip that cluster, and do not reconnect (/beagle3 and /project stay readable through the other).
2. midway2: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu 'bash ~/watch_probe/watch_mdw2.sh' < /dev/null` (read-only probe, md5 cf834d4b). The probe lists only the last 80 minutes: if the "Last pass" line of /project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md is older than that, also run `sacct -X -n -P -u yinhanw -S <that last-pass time> -o JobID,JobName,State,End,ExitCode` on midway2 and read every job that left the queue since. Look for: any job that left the queue not COMPLETED (read its log tail: Peng arms /beagle3/trsosnic/yinhan/upside3_validation/peng/logs/<prot>_<kind>_ff_2.1_<id>.out for u3_*, /beagle3/trsosnic/yinhan/upside3_validation/peng_seeds/logs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>_<id>.out for the seed replicates u3s*/u2s*, /beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<kind>_<ff>_<id>.out for b_ff_2.1_* and b_ff_3.1_*; collectors /project/trsosnic/yinhan/checks/ff31_merged_20261026/logs/ and checks/ff21_released_20261021/; Upside 3 glpG REMD logs /beagle3/trsosnic/yinhan/upside3_validation/<glpg|glpg_ff31|glpg_ff30>/outputs/martini_glpG-RKRK-79HIS_hybrid/REMD/glpG-RKRK-79HIS.remd_block<k>.log); a Peng chunk COMPLETED without a "[bench] resubmitted" or "[bench] COMPLETE" line; on "[bench] COMPLETE", the arm's final frames (`source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh; export HDF5_USE_FILE_LOCKING=FALSE; nice python3 -I ~/watch_probe/final_frames.py <run dir>`: every coordinate finite, KE/1.5kT near 1.00); error lines in the newest REMD block log of each of the three glpG runs (the probe reads glpg and glpg_ff31; for glpg_ff30 `grep -ciE "error|nan|traceback|abort"` on its newest block log); a REMD block ending (block_count rising, or a block job leaving the queue: ff2.1 block 3 49584589, ff3.1 block 1 49573354, ff3.0 block 1 49589829 with block 2 49634605 queued behind it): report it; its readout is owed by hand (remote_jobs.md §1 line 3: TM4 DSSP, KE and C-N at T 0.70, against Upside 2's ff3.0 0.96); peng_ff21.txt or peng_ff31.txt existing (report its table); /project/trsosnic under 200G free.
3. midway3: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu 'df -i /beagle3/trsosnic | tail -1; squeue -h -u yinhanw -o "%i %j %T %S %r"' < /dev/null`. The live midway3 jobs are the E1 jobs of remote_jobs.md §1 line 4 (job table); when one leaves the queue, report its State and the `[e1] ... FAILED` lines of its log /beagle3/trsosnic/yinhan/boltz_hdx/logs/<job name>_<id>.out (the readout is owed by hand on the MacBook Pro). Report IFree under 30000.
4. Record: when a job finished, failed or changed state, update remote_jobs.md §1 (the line it belongs to and the job table; snapshot time) briefly, and set the "Last pass" line of /project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md on midway2 to this pass's time (owner <this computer>, watch cron in this session). Do not record unchanged state. After any record change, refresh the handoff bundle's copies from the repo root so a session on another computer reads them (remote_jobs.md "Resume here"), but first read WATCH_STATUS.md: if another computer passed since this one's last sync, back up `records/` as `records.bak_<computer>_<date>` before the rsync (findings 10.46): `rsync -a -e "ssh -o BatchMode=yes -S $HOME/.ssh/cm-mdw2.sock" remote_jobs.md plan.md findings.md progress.md diffusion_dynamics.md yinhanw@midway2.rcc.uchicago.edu:/project/trsosnic/yinhan/checks/handoff_20261026/records/`.
5. Reply: if nothing needs the user, one line (time; midway2 queue counts; midway3 states; disk and inodes). Otherwise a short report of what changed and what is owed. Name runs by what they are, with the tag in parentheses."

**Round-4 watch.** The by-hand checks of training runs, which the collectors do where one runs; a run without one
(the ff3.1 training debug) is read by hand. `<run>` is a run directory under `$P/training`; `submit_new.sh` maps runs
to their tags.
1. **Socket** (§0): run `scratchpad/mdw2_master.exp` once if it is dead, unless §1 records that the previous pass
   already tried; still down, record it in §1, notify once and stop. Never loop on connection attempts.
2. **Queue and logs.** `squeue -u yinhanw`; `bash /project/trsosnic/yinhan/ff3_selection/submit_new.sh` submits each
   new epoch end's panel for the runs it lists, without asking (user, 10-08 11:05). Check each run's newest link log for
   `WORKER_FAIL`, `Traceback`, `STOPPED`, `never started`. Benchmark and glpG logs: `ssh ... 'bash ~/watch_probe/logs.sh'
   | grep -v -E "done=[0-9]+ / target" | grep -v -E "^== .*errs 0$|^\[remd\] chunk [0-9]+ steps"`; no output (exit 1)
   is clean; a chain's `block k/5` and `calibration` lines are expected.
2b. **Disk.** `df -h /project/trsosnic`; record it in §1 "Disk" when it moves by more than 100 G, and tell the user
   before it falls under 200 G (10-11: it reached 88 G unnoticed).
2c. **A job that left the queue:** `sacct -n -X -j <id> -o JobID,JobName%45,State,End,Elapsed,ExitCode` and its log
   tail. A benchmark arm (log `/beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<kind>_<ff>_<id>.out`) stops each
   chunk at its 10.5 h wall with `[bench] resubmitted: <id>`: check the new job carries the law (`scontrol show job
   <id> | grep -o "ExcNodeList=[^ ]*"`). On `[bench] COMPLETE`, its final frames: `source
   /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh; export HDF5_USE_FILE_LOCKING=FALSE; nice python3
   ~/watch_probe/final_frames.py /beagle3/trsosnic/yinhan/ff3_benchmark/runs/<prot>_<kind>_<ff>` (finite, KE/1.5kT near
   1.00). An Upside 2 glpG chain (log `popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<id>.out`) ends a block with `[remd] block
   complete` after `Submitted batch job <next>`: read the last `avg_kinetic_energy/1.5kT` line and the next job's
   `block k/5` line. Record new ids in §1 and `WATCH_STATUS.md`.
3. **Every step finished since the last pass:** `cd $P/training && source <run>/env.sh && python3 check_step.py <run>
   <epoch_xx_minibatch_yy>` under `ulimit -v 4000000` (24 of 24 proteins); the KE scan `grep -m1
   ^avg_kinetic_energy/1.5kT run_output/<step>/*.run.*.up.output`, flagging a maximum above 1.2 (findings 1.21); record
   the shared margin E_other - E_alpha (start +0.192; ff3.2's +0.082), dhb, sheet mean, rot rms, step time, helical /
   left glycine alpha_L free/restrained, and for depth runs dL - dR (start +0.554) and `run_output/rama_rounds.txt`. A
   trained-H-bond run whose margin falls below +0.10: notify the user once and hold; nothing is cancelled without the
   user. **The steric check of depth runs** (user 10-10, after Jumper; plan.md Phase 14 step 5, findings 1.33): `nice
   python3 /project/trsosnic/yinhan/checks/ff31_steric/steric_check.py <run> [<step>] >
   /project/trsosnic/yinhan/checks/ff31_steric/<run>_<step>.txt` with the run's env under `ulimit -v 4000000` (3 min,
   360 MB; repo `training/ff31/steric_check.py`, md5 `e35e3d98...`): record the side-chain clash probability now /
   start, the largest map and residue ratios and the non-glycine free beta' + pPII'. A `STERIC FLAG` (a coil map or
   training residue whose clash probability more than doubled against NDRD, the training residues' mean side-chain
   clash probability up by more than 10%, or free non-glycine beta' + pPII' above the run's steps 0-18): notify the user
   once, with the table; the run keeps training.
4. **A panel or a gate.**
   - A panel with 176 npz in `ff3_selection/runs/<tag>/`: `panel.py select aa domains runs ff21_released ff21_awh
     gdepth_start <tags>` with the run's env (70 s, 300 MB); notify the table. Each rerun's panel goes beside its
     source's (`a31_EE` with `fd_EE`, `a32_EE` with `f32_EE`, `c9_EE` with `c9_00`/`c9_01`). Panels whose `.up` files
     were deleted (10-18) compare on folded fraction and Q only, by `checks/converge_20261020/paired_npz.py DOMAIN_DIR
     AA_DIR RUNS_ROOT TAG_A TAG_B [CSV]` (panel_compare's folded and Q from the npz, reproduced exactly; from
     `ff3_selection/` with `training/ff30_gly/env.sh`), never by `panel_compare.py`.
   - A gate: at the target `train_chain.sbatch` cancels its queued successor and submits the run's
     `after_training.sbatch`, which writes `<run>/gate_step<N>.txt` (`convergence_gate.py` over the last epoch); not
     converged, one more epoch (`train_chain.sbatch <run> <N+19>`) up to step 114, then "stopped for review". For the
     ff2.1 runs a converged gate put the epoch-end files in `ff3_selection/params/<tag>` and submitted the mdCATH panel
     (`ff3sel_mdc_<tag>`, 356 runs in `ff3_selection/runs_mdcath/<tag>/`); a gate that runs `validate_ff.sh` (§1b)
     queues 32 Peng arms and 4 glpG chains: read their first logs. Notify the gate either way, with the failing groups.
     A run with a target and no gate (the ff3.1 training debug) stops for the user.
   - Two mdCATH panels: `cd ff3_selection && sbatch $(cat slurm.args) panel_compare.sbatch domains_mdcath aa_mdcath
     runs_mdcath <tag_a> <tag_b> mdcath_candidates/candidates_mdcath.csv` (read-only; `panel_compare_<id>.out`: folded,
     mean Q, CA-RMSD, paired, by length bin and fold class); tags without `.up` files by `python3 panel.py select
     aa_mdcath domains_mdcath runs_mdcath <tags>` on login1.
5. **Update §1 briefly, and rewrite the cluster's `WATCH_STATUS.md`** ("Handoff"): owner and computer, time of this
   pass, the newest step of each run with its margin, and anything owed or awaiting the user.

## 0. Connect first (needs a Duo push on the user's phone)

Password + Duo is the only method. Key-based ssh is refused on both clusters: tested 2026-09-17 with a correctly
installed `authorized_keys` (600/700) and again with `PubkeyAcceptedAlgorithms=+ssh-rsa`; the server advertises
`publickey` but does not honour it for this account. Do not try again. The inert key entry left in midway2's
`~/.ssh/authorized_keys` is harmless.

**midway2** (training, selection panels, benchmark arms, all job submission):
```bash
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp    # USER MUST APPROVE DUO
```
By hand from a new computer: `ssh -M -S ~/.ssh/cm-mdw2.sock -o ControlPersist=8h -o ServerAliveInterval=30
yinhanw@midway2.rcc.uchicago.edu` (password and Duo). To let Claude open it unattended, copy
`scratchpad/mdw2_master.exp` over and put the password line `set password "..."` in `~/.bin/ssh_mdw3`. If midway2
blocks or throttles this IP, tunnel through midway3. It needs a live `~/.ssh/cm-mdw3.sock`, proxies the TCP leg
through it (`ProxyCommand=ssh -S cm-mdw3.sock -W`), and holds `~/.ssh/cm-mdw2.sock` only while the process runs, so
start it in the background:
```bash
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_via_mdw3.exp &   # USER MUST APPROVE DUO
```

**midway3** (file moves and reads on `/project` and `/beagle3`, GPU jobs):
```bash
ssh -o BatchMode=yes -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp    # if not: USER MUST APPROVE DUO
```

**Every routine call: check the socket first, and always pass `-o BatchMode=yes`.**
```bash
if ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu >/dev/null 2>&1; then
    ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<cmd>'
else
    echo "master down - do NOT retry blindly; see throttle note"
fi
```
Why both (learned 2026-09-17 and 09-19): a plain `ssh -S sock host 'cmd'` on a dead master falls back to a fresh
connection, offers keys, tries password auth twice non-interactively and hits `Too many authentication failures`.
`BatchMode=yes` alone protects only when the socket file is gone; if the master died but `~/.ssh/cm-mdw2.sock` is still
on disk, ssh tries the stale socket and then falls through to a real attempt anyway (`Permission denied
(publickey,...)`). A few of those and the host answers `Connection closed by 128.135.112.69 port 22`: the IP throttle,
self-inflicted by a status check, after which the next expect launch spends a Duo push before it discovers the block.

Rules that follow:
* **When throttled, stop.** The throttle clears by itself in tens of minutes. Wait at least 30 min before one single
  retry. The way around it meanwhile is the midway3 tunnel, which needs a live `~/.ssh/cm-mdw3.sock` and therefore its
  own Duo approval.
* **Run the expect script at most once per attempt**, and never let a monitoring loop retry it more than once per
  tick. Each launch is a Duo push on the user's phone (sent twice by mistake on 2026-09-18).
* **Any automated credential send is one-shot.** On 2026-09-18 `mdw2_hold.exp` used `exp_continue` on the password
  prompt, re-sent a rejected password on every re-prompt, and one launch spent the account's failed-attempt budget.
  Send the password at most once and treat a second prompt as a hard stop (`PASSWORD_REJECTED`); match `Permission
  denied` and `Too many authentication failures` and exit; pass `-o PubkeyAuthentication=no -o IdentitiesOnly=yes`
  (every key the client offers wastes an attempt against `MaxAuthTries`) and `-o NumberOfPasswordPrompts=1`. An expect
  script has no dry-run check: `expect -n -c 'source ...'` still connects.
* **Minimise connections rather than trying to remove Duo.** The training chain insures itself (§8), so a status check
  every few hours is enough. No cluster-side monitor runs; `/project/trsosnic/yinhan/STATUS.md` is stale from
  2026-09-27 23:54.
* **Recurring work on the cluster.** `scrontab` is disabled, `crontab` is denied, and `pi-trsosnic` has no
  association with the `cron` partition (`sbatch -p cron` is refused). Use a 1-core `broadwl` job of at most 36 h (7
  days hits `QOSMaxWallDurationPerJobLimit`) that queues its successor at the start
  (`--dependency=afterany:$SLURM_JOB_ID`); one that resubmitted at the end of its loop overran the wall and died
  silently.
* **Laptop sleep is not what drops the master.** With `caffeinate -is` holding `PreventSystemSleep` it still died
  after ~40 min; `mdw2_master.exp` now uses `ServerAliveInterval=30 ServerAliveCountMax=20 TCPKeepAlive=yes` (~10 min
  tolerance). When the socket keeps dropping, look first for our own load on the login node (`ps -u yinhanw
  --sort=-rss`, findings 10.6). **Never run anything heavy on a login node** (findings 10.6).
* zsh does not word-split unquoted variables: `M2="ssh -S sock host"; $M2 'cmd'` runs silently and prints nothing,
  which is not a connection failure. Write the `ssh` call out in full.
* `timeout` does not exist on this Mac; do not wrap `expect` in it.

**Paths.** `~/project` on midway2 is a symlink to `/project/trsosnic` (not `/project/trsosnic/yinhan/`); glpG data is
at `~/project/yinhan/popepopg_REMD_mdw2/`. `/project` is the same filesystem on midway2 and midway3 (both mount
`midway3_cap`), so checkpoints, chain logs and force-field directories under `$P` can be read from midway3 without a
midway2 login. Only `squeue`/`sacct` need midway2: the clusters have separate Slurm controllers and accounting
databases, and `sacct --clusters=all` on midway3 does not see midway2 jobs (nor midway2's `squeue` midway3's). Python
env on midway2: `source /software/modules/init/bash && module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1 && export
HDF5_USE_FILE_LOCKING=FALSE`. While login2's `/project` hangs, run `/project` work as `... 'ssh -o BatchMode=yes
midway2-login1 bash -s' < script.sh` (§8). `~/project` on midway3 is a symlink to `/project/trsosnic/yinhan/` (yinhan,
not yinhanw); load the python env with `source ~/project/NP-1AO6/env.sh` before any h5py work there.

---

## 0a. `broadwl-lc` is empty but unusable: its nodes cannot see `/project` (2026-09-18)

Do not send work there, however idle it looks. Its nodes advertise `AvailableFeatures=lc,e5-2680v4,64GB,noib`: no
InfiniBand, and `/project` and `/beagle3` are served over it (a probe job on `midway2-0213` reported
`PROJECT_MISSING`, `PROJECT_NOT_WRITABLE`, `BEAGLE3_MISSING`). Every input we use is on `/project`, so the partition
is dead for this project, the same trap as `/cds3`. The failure is silent: a job whose `--output` path is on
`/project` dies at once with `ExitCode 0:53`, `Elapsed 00:00:00` and no log file, because it cannot create the log; a
`--wrap="hostname"` probe of the same shape fails identically, so it is not the script.

Related traps: `sinfo -F`'s "0 idle" counts a partly filled `mix` node as allocated, so use `sinfo -p <part> -o
"%.8t %.6D %.10C"` for A/I/O/T cores; and `/tmp` on a compute node is node-local, so a probe writing there "succeeds"
and leaves nothing the login node can read (probe with `--output` under `$HOME`, which is shared). A job_submit plugin
silently rewrites `#SBATCH --partition=broadwl-lc` to `broadwl` while honouring every other directive; command-line
`-p` is honoured. Verify placement with `scontrol show job <id> | grep Partition`.

---

## 0b. Shared Upside deployment on beagle3 (2026-09-09)

`/beagle3/trsosnic/yinhan/upside2-md` is the one tree and binary both clusters use; how and why it works is in
CLAUDE.md, "Shared Upside Deployment". It was deployed by rsync from `$P` (`src py obj parameters cmake example` and
the install scripts), excluding `training/`, which is ConDiv campaign data rather than code. Benchmark output goes to
`/beagle3` so it does not compete with `/project`. A release writes `parameters/<ff>` in both trees. The old
`ff_3.0_trained` / `ff_3.0_trained_rf` directories are still on the cluster but no live script references them
(checked 2026-09-28).

**Before deploying an Upside 3 engine (ledger A139 onward) to either tree** (2026-10-22): the cluster trees run
pre-audit binaries, and the live hybrid files are written for them.
* From A160 a production `.up` needs `/input/hybrid_control@po4_z_hold_end_round` (the burn-in's round + 150) or the
  engine refuses it, and a positions-only start (each REMD segment) no longer holds PO4 z; from A162
  `preprod_protein_mode` is gone and every barostat attribute, the stage label and `/input/mass` are required. Live
  glpG chains must finish on their own binary or have their files re-prepared.
* A191-A192 change every hybrid's model (DSSP backbone types; `parameters/ff_2.1/martini.h5`'s side-chain table
  rebuilt from martinize's mapping, md5 7496dd90): deploying them means re-preparing the hybrids from stage 0, not
  reusing their stage files. From A222 a pulling config needs `MovingConst3D@step_initial` (rewrite it with the current
  writers).
* ConDiv needs A150 (rama `set_param`) in its engine. From the glue review's D6 (10-22) the side-chain latent drops its
  86 dead `rotscalar` coordinates (31420 to 31334), so a ConDiv checkpoint written before it does not resume under the
  new `py/rotamer_parameter_estimation.py`: deploy it only to a tree no running campaign reads (CLAUDE.md's rule), and
  start new campaigns fresh. The trainers no longer import torch (A204).
* Deploy `py/upside_engine.py` and `obj/libupside.so` together (A197 changed the C API's signatures).
* The cluster CMake file changed (A205, A206; verified with GCC 16 on the Mac, not yet with midway2's GCC 10.1), so an
  Upside 3 build on midway2 checks its link lines.
* **`set_fast_math`, open:** every cluster binary (`obj/upside` and `obj/libupside.so` in both trees, `nm`, 10-22; and
  the Upside 3 tree's build, §1 line 3) carries it. GCC before 13 links `crtfastmath.o` whenever `-ffast-math` is on
  the link line, `-shared` included, so loading the library sets flush-to-zero and denormals-are-zero in the whole
  Python process, and the cluster executable runs with both while the Mac's does not. The fix to measure is
  `-ffast-math` on the compile lines only, with a regress compare of the two midway2 builds.

### Other Upside copies (inventory 2026-09-09)

Owned by yinhanw:

| path | size | branch | last touched | disposition |
|---|---|---|---|---|
| `/project/trsosnic/yinhan/upside2-md-mdw2` | 8.3 G | martini-dev | live | ACTIVE, the ConDiv training tree (`$P`) |
| `/beagle3/trsosnic/yinhan/upside2-md` | 6.4 G | martini-dev | live | the shared deployment |
| `/scratch/midway2/yinhanw/upside2-md-water-diff` | 13 G | water-diff | 2025-10-30 | idle; see below |
| `/scratch/midway2/yinhanw/upside2-md` | 176 M | master | 2026-01-28 | keep, master |
| `/home/yinhanw/upside2-md` | 256 M | master | 2025-11-19 | keep, master |

The Upside 3 tree `/beagle3/trsosnic/yinhan/upside3` is described in CLAUDE.md and §1 line 3.

**Not ours, never touch:** `/beagle3/trsosnic/upside2-md` and `/project2/trsosnic/upside2-md` (bayhi),
`/project2/trsosnic/software/upside2-md` (nffaruk), `/project2/trsosnic/pengxd/*` (pengxd), plus copies under baxa,
avmolina, ruofan, tobin, schwartznw, yiheng, zonganw, simoneritchey, amz and bayhi.

`upside2-md-water-diff` is the only deletion candidate, and it is not purely a code copy. Its source is safe (HEAD
`c44a404` is on `origin/water-diff` with 0 uncommitted files), but 11 of its 13 GB is
`example/16.MARTINI/outputs/water_T*`, water-diffusion output from 2025-10-30 that exists nowhere else. Left in place
pending the user's decision.

---

## 0c. Rotamer belief-propagation stopping-test bug (2026-09-26)

**The bug** (reported by John Jumper via Tobin). `NodeHolder::max_deviation` in `src/rotamer.cpp` took
`max(cur_belief - old_belief, dev)` with no absolute value. Beliefs are rescaled each iteration so the dominant rotamer
stays at 1, so an update that only sharpens (others fall) gives dev = 0 and the solver stops as if converged. It has
been there since the original rotamer commit: ff2.1 and ff30_basin (ff_3.0) were trained with it.

**The fix** is `fabsf(...)` on that line (279): `master` in `288d1fec` (pushed by the user 2026-09-26); `martini-dev`
by merging master in (PR #48, `11044069`, 2026-09-26); both cluster trees (`$P` and `/beagle3/trsosnic/yinhan/upside2-md`)
installed by the ff3.0 gate 2026-09-30 02:24:59. Both trees carry `BPFIX_INSTALLED` and keep the old binaries as
`obj/*.bak_pre_bpfix_20260930-022459` (gate log `ff30_basin/ff30b-gate_49131949.out`, `[bpfix]` lines).

**Static frames** (local, rerun on fresh frames 2026-10-02): 8 ConDiv training proteins (2xf6, 3h36, 1afh, 2lpn,
3mwz, 1a62, 3dm8, 1w4s; 50-146 res), ff2.1, 91 frames each (native + 30 each from 2000-unit runs at T 0.80/0.95/1.10),
buggy and fixed library at the production tol 1e-3 against the fixed library at tol 1e-6 (max_iter 20000). The bug
cuts the solve short in 338 of 728 frames (46%), mostly by 2 iterations. Energy error |E - E_ref| median 7.5e-4 E_up
(fixed 4.3e-4), 99th percentile 0.0095, max 0.157 (1afh frame 21, a real early stop, 0.2 kT at T 0.80); the fixed
solver's own max is 0.044 (the 1e-3 tolerance alone); E_bug < E_ref in 95% of frames, mean -0.0007 E_up. RMS force
error / RMS force median 8.3e-5, max 2.9e-3 (fixed 4.6e-5, 5.3e-4). The ConDiv contrastive rotamer gradient (native
minus decoy mean of dE/d pair parameters over the 8 proteins): bug - fix 1.05e-3 of its norm, turned by 1e-3 rad,
against a sampling noise of 7.2%, so 68x below the noise; per coordinate below 0.08 of its own noise in 99%. Scope:
local, at ff2.1's parameters, rotamer pair group only; the direction further training would take, not the end point of
a long retraining. Conclusion given to the user: ff2.1 and ff3.0 do not need retraining.

**Simulation validation for Tobin** (16 arms complete 2026-09-29, analysed; not yet reported to Tobin, the user's
call): ff_2.1, native start, 14-replica REMD with the Peng benchmark protocol (Table S2 ladder and duration, dt 0.009,
frame 100), proteinG / homeodomain / WWdomain / NTL9 x {bug, fix} x seeds {1, 2}, read as bug-minus-fix against the
seed spread (z = delta / sigma_seed), since the two binaries' trajectories decorrelate within a few hundred time units.
Output `bp_validation/analysis_20260929.txt`. WWdomain is the only protein with a consistent difference: at T
0.764-0.879 the bug arms are less native (Q lower by 0.02-0.04, RMSD higher by 0.1-0.3 A), T_mid 0.859 bug vs 0.861
fix, both bug seeds below both fix seeds in Q from 0.700 to 0.879; suggestive, not established at two seeds, and its
implied ddG (~0.15 E_up) is ~200x the median static-frame |dE| and not explained. NTL9 goes the same way at low T but
rests on one arm; proteinG and homeodomain are not resolved (|z| < 1.5). The 224 cells match the no-effect null
(Student-t, 2 dof: |z| > 2 in 41 against 41.1 expected), so no single z is evidence by itself. The fix costs 0-5% in
speed.

Files in `/beagle3/trsosnic/yinhan/bp_validation`: `obj_bug/`, `obj_fix/` (deployment source built twice, differing
only in the fix; `obj_bug` reproduces the old deployed binary bitwise), `bp_run.py` (copy of `ff3_benchmark/bench_run.py`:
binary, seed, ff_2.1 native), `bp.sbatch`, `submit_all.sh`; logs `logs/<prot>_<build>_s<seed>_<jobid>.out`, data
`runs/<prot>_<build>_s<seed>/`. `analyse_bp.py [protein ...]` (midway2 login node, after `source
/beagle3/trsosnic/yinhan/upside2-md/env_shared.sh`) gives per-T Q, CA-RMSD, Rg and E per build, bug-fix vs the seed
spread, melting midpoints and speed, and writes `obs_<protein>.npz`. **Figure and deck (2026-10-02, unsent):**
`bp_validation/static/bp_validation.{png,pdf}` from `plot_bp.py` (11 panels: static-frame errors and cut-short
iterations, Q(T) per seed, Q_bug - Q_fix with the seed sigma; red = bug, blue = fix; SciencePlots style files in
`styles/`); the 2-slide deck `bp_validation_slides.pptx` was built by `scratchpad/bp_validation/make_slides.py` on the
computer that made it. `static/` holds everything the Mac had (325 files, md5 manifest identical): the figure inputs
(`remd/obs_*.npz`), the static-frame harness (`make_src.py`, `build_libs.sh`, `gen.py`, `eval.py` run as `bug 1e-3
1000 bug`, `fix 1e-3 1000 fix`, `fix 1e-6 20000 ref`), its sources and arm64 libraries, and the results
(`sys/<code>/`, `res_{bug,fix,ref}.npz`). Replot from any computer by copying `static/` and running `python
plot_bp.py` in it.

---

## 0d. midway2 node exclusions: one law (2026-10-08)

**The only midway2 node list is `/project/trsosnic/yinhan/slurm/midway2.args`** (user, 10-08):
`--partition=broadwl --exclude=midway2-[0003,0010-0011,0027,0033-0035,0037,0060,0080,0085,0088,0103,0116,0151,0316-0327,0342-0345]`.
Each node's record of failures is in `README.md` beside it; §8 has the history. The later additions, each at the
user's word: midway2-0088 10-09 04:12, midway2-0103 and 0116 10-09 10:30, midway2-0151 10-13 10:29 (two NODE_FAILs
10-10), and midway2-[0316-0327] 11-02 18:10 (every busy core held at 1.2 GHz, §1 line 2); backups
`*.bak_pre_0088_20261009`, `*.bak_pre_0103_20261009`, `*.bak_pre_0151_20261013`, `*.bak_pre_0316_20261102`. Each was checked
with `--test-only`: the wrapper refuses `-w` on the new node and accepts a neighbour, and `update_pending.sh` gave the
pending jobs the new list. Before 10-08 the list was copied by hand into 13 `slurm.args`, `submit_remd.sh`,
`bench.sbatch` and a dozen static `#SBATCH` lines in three versions, and a queued chain never saw a node added later.

Where it is read, always at submission time:
* Every midway2 `slurm.args` is a symlink to it (the round-4 runs, ff30_glyhb, ff30_gly, ff30_basin, `ff3_selection/`,
  `polygly/prod/`; backups `slurm.args.bak_pre_law_20261008`). `train_chain.sbatch` reads `slurm.args` at every link,
  so a running chain takes a change at its next link; its gate, `validate_ff.sh`, `submit_new.sh` and `prod.sbatch`
  read it too.
* `popepopg_REMD_mdw2/submit_remd.sh` takes `$(cat <law>)`; every REMD block resubmits through it.
* `/beagle3/.../ff3_benchmark/bench.sbatch`: on midway2 a chunk resubmits with the law; on caslake it carries its own
  partition and `--exclude`. Its static `#SBATCH --exclude` is gone.
* `~/bin/sbatch` on midway2 (a link to `slurm/sbatch_wrapper`): `~/bin` comes before Slurm on PATH in login shells,
  ssh commands and jobs, so every other submission from a midway2 host gets the law's `--exclude` first, and a
  script's own `#SBATCH --exclude` is overridden. sbatch keeps the last `--exclude` given, so a caller's own
  `--exclude`, or the real `/software/slurm-current-el7-x86_64/bin/sbatch` by path, is the deliberate way around it.
* Left alone: ff30_glyprobe's `slurm.args` (extra flags; cancelled, no gate), the caslake runs' files, and finished
  one-off scripts with static lists (`ff3_benchmark/scoring/`, `gly_peptides/`, `polygly/collapse.sbatch`, `checks/`),
  which the wrapper overrides on midway2.

**Adding a node.** A node is excluded for its own record (a second failure on it), not for being near another. Edit the
one line of `midway2.args`, add the node's row to `README.md`, then run `bash /project/trsosnic/yinhan/slurm/update_pending.sh`
on a midway2 login node: it sets every pending job's ExcNodeList to the new list and prints each. Running jobs keep
their nodes. A new run's `slurm.args` is `ln -s /project/trsosnic/yinhan/slurm/midway2.args slurm.args` (`cp -a` keeps
the link; plain `cp` copies the text, which then drifts).

**Checked 10-08** (no test jobs; `--test-only` and stubs): the law parses; through the wrapper `-w midway2-0027` and
`-w midway2-0343` are refused and `-w midway2-0026` accepted, while the real sbatch by path accepts 0027, so the
refusal is the wrapper's; a script with its own old `#SBATCH --exclude` is still refused 0027; all 13 `slurm.args`
resolve to the law; with a stub sbatch, `submit_remd.sh` passes the law, and `bench.sbatch`'s resubmit takes the law on
a midway2 host and its own flags on a midway3 one; `update_pending.sh` set all 13 pending jobs, their dependencies
intact although squeue briefly showed reason "None"; real submissions (a panel, and a chain successor submitted from a
running link) carry 0027 in ExcNodeList.

## 1. Current jobs

Snapshot **2026-11-02 16:45 CST** (the MacBook Pro; status pass, not the owner, no cron; ff3.2's entries in line
2 updated 17:40 by the Mac Studio's ff3.2 session, the rest as of 12:20). Owner none since the
MacBook Pro's handoff 21:55 10-27 ("Handoff"). **midway2 was faulty 21:19 11-01 to 09:43 11-02**: srun answered
"Unable to confirm allocation ... Unexpected message received" and sbatch "Batch job submission failed: Unexpected
message received"; from 23:57 its login node opened no session, and compute nodes lost the controller too ("Unable to
send job complete message"). Jobs already running kept writing. At 09:43 the controller returned and stamped every job
that ended meanwhile TIMEOUT at 09:43 (`sacct -S 2026-11-01T21:00`: the 32 ff3.1 Peng chunks and both Upside 3 glpG
blocks, whose logs end normally); ff3.2's link 49737476 FAILED, its insurance 49737477 R from 09:45. The fault broke
the resubmission of 30 ff3.1 Peng chunks (23:22-23:28 11-01 and 03:16 11-02), restarted 12:18 11-02 (line 2). Midway2's
login answered again 12:16 (one Duo push). Earlier, between 23:25 10-31 and 14:10 11-01, ff3.1's
alpha3d native 49712331 lost midway2-0156 (NODE_FAIL 06:11, requeued and resumed); with cspA de novo's midway2-0189
(20:11 10-31), each node's first failure on record, so neither is excluded. `/project/trsosnic` 312 G free, 626,074
inodes free; `/beagle3` 118,631 inodes free (12:20 11-02); at 17:00 207 G and 621,037 inodes free, `/beagle3` 118,152
inodes free (105 G used in under 5 h; the watch reports below 200 G). Paths: `$P` =
`/project/trsosnic/yinhan/upside2-md-mdw2`; `checks/`, `ff3_selection/` and `popepopg_REMD_mdw2/` are under
`/project/trsosnic/yinhan/`; `ff3_benchmark/` and `upside3_validation/` under `/beagle3/trsosnic/yinhan/`.

**Line 1: deterministic vs probabilistic side chain**
- **Released ff2.1 on the three side-chain tests** (user 10-21: "run validation with released ff2.1 so we know if we
  should keep ff2.1 or not"; plan.md Phase 27; findings 1.37), submitted 10-21 by
  `checks/ff21_released_20261021/submit.sh` (ids in `submitted.log` there); files `ff3_selection/params/ff21_released`
  (md5 equal to `/beagle3` `ff_2.1`), each test's own protocol, binary and domains at dt 0.009.
  - 44 domains (`p44_ff21_released`): released ff2.1 folds 0.603 against the control's fixed point 0.430 and no-BP's
    0.363; mean TM 0.681 / 0.650 / 0.623 (`p44_dt009_folded.txt`, `tm_p44_*.txt`).
  - mdCATH (`runs_mdcath/mdc_ff21_released/`): on 89 domains mean TM 0.727 against the control's 0.702 (+0.025
    [+0.018, +0.033], p 6e-09; higher on 70 domains) and no-BP's 0.667 (+0.060, p 3.5e-14; higher on 80); TM > 0.5
    0.972 / 0.946 / 0.920 (`tm_mdcath_mdc_ff21_released_vs_<mdc_ctrl_03|mdc_nobp_03>.txt`).
  - Peng: 30 arms `b_ff_2.1_<prot>_<kind>` plus lambda's two arms of 10-04, each resubmitting itself in 10.5 h chunks
    (logs `ff3_benchmark/logs/<prot>_<kind>_ff_2.1_<id>.out`, runs `ff3_benchmark/runs/<prot>_<kind>_ff_2.1/`). 31 of 32
    scored, the last proteinL de novo 12:55 10-27; the five de novo arms finished from 11:58 10-26 read final frames
    finite, KE/1.5kT 1.000-1.007 over 14 replicas. Ubiquitin de novo, the last, COMPLETE 16:33 10-29 (49639439):
    final frames finite, KE/1.5kT 1.002-1.007 over 14 replicas.
  - **32 of 32 scored**: collector `ff21_collect` wrote `checks/ff21_released_20261021/peng_ff21.txt` ending DONE (its
    scores in `ff21_nobp_dt009/setup/score/ff_2.1/`). Mean TM at the bottom rung over the 16 proteins, released ff2.1 /
    control ctrl38 / no-BP / ff3.0: native 0.552 / 0.486 / 0.457 / 0.480 (Peng's published FF2 0.55), de novo 0.427 /
    0.299 / 0.309 / 0.305 (0.42). These are means of the table's rows; its summary lines read NaN with no paired test,
    since they take only proteins scored in every column and the retired `ff_3.1_depth_data` column scores none. Then
    plan.md Phase 27 step 4.
- **Both side-chain trainings are finished** ("Finished runs"); no jobs.
- **A third arm, the most-probable-state side chain `ff21_bpmax_dt009`** (user 11-03, the PI's request; plan.md
  Phase 39): BP runs, each residue keeps only its state of largest marginal. Trained from released ff2.1 (init_param
  md5-equal to `ff_2.1`) with the control's inputs, trainer and minibatch order, and the earlier arms' engine history
  (user 11-03: "the program should be the same state as the other tests"): steps 0-37 on the engine before the
  angular-gradient fix, 38 onward on the fixed one, both with the mode.
  - **Running: link 49766810** (R from 12:42 11-03 on midway2-[0017-0018,...]; target 38, no gate; log
    `/project/trsosnic/yinhan/upside2-md-bpmax/training/condiv-train_49766810.out`), insurance successor 49766811
    (afterany). At 38 `after_training.sbatch` -> `setup/at_target.sh` -> `setup/switch_engine.sh` switches the tree's
    `obj` symlink from `obj_preangfix` to `obj_angfix` and submits the link to 76 (`setup/engine_history.log`); from 76
    the gate and, converged, the mdCATH panel `mdc_bpmax_EE`; not, one more epoch, up to 114. Its workers were seen
    running `obj_preangfix/upside` and its `libupside.so` at 12:45.
  - Setup 49766481 COMPLETED 12:29 (step-0 parameters, solver and minibatch order identical to the control's). The
    first link 49766483 (target 76, fixed engine from step 0) ran 12:33-12:37 and was cancelled in step 0 to match
    the engine history; its successor 49766605 cancelled before starting; partial step deleted (no checkpoint). Its
    336 worker configs all carried `most_probable_state = 1`.
  - **Its own tree, not `$P`:** run dir `/project/trsosnic/yinhan/upside2-md-bpmax/training/ff21_bpmax_dt009`
    (README there). The tree's `obj` (a symlink, see above) carries the patched engine and its `py/upside_config.py`
    writes `most_probable_state = 1`, so every config written there is in the mode; `$P` is untouched. **Any
    evaluation of this force field must run in that tree** (source its run's `env.sh`): in `$P` or `/beagle3` it would
    silently run as BP. The mdCATH panel does (`ff3_selection/panel_mdcath_bpmax.sbatch`); the 44-domain panel at dt
    0.009 and TM-scoring against `mdc_ctrl_03` / `mdc_nobp_03` and `p44_ctrl_03` / `p44_nobp_03` are owed by hand, the
    same way. Never switch `obj` by hand while a link of this run is in the queue.
  - The engine change is not in git (user 11-03): patch and tests in the repo's git-ignored
    `scratchpad/bpmax_20261103/`, the cluster's checks in `/project/trsosnic/yinhan/checks/bpmax_20261103/`
    (`cluster_parity.sh`, `cluster_parity2.sh`: each patched build in BP mode is bitwise the binary the earlier arms
    ran at the same steps).
  - **Next:** after step 0 (about 36 min from 12:42), the step's log shows no worker failure; at 38 (about 23 h), the
    engine switch in `gate_<id>.out` and `engine_history.log`, and the link to 76 queued. The mode runs 1-2% hot by
    design (findings 1.38); a protein-step's KE/1.5kT above 1.05 is not expected and is to be investigated.

**Line 2: Ramachandran map training**
- **ff3.1 is `ff31_nativefix` at step 76, released 16:49 10-31 and under validation** (user 10-31: "Deploy ff3.1, run
  validations (Peng, glpG, and other tests we ran)", the retrain chosen over the merge; plan.md Phase 37). Scripts,
  logs and ids in `/project/trsosnic/yinhan/checks/ff31_nativefix_release_20261031/` (`release.sh`, `release_apply.log`,
  `submit.sh`, `submitted.log`, `collect.sbatch`, `tm.sbatch`). Files: ff_2.1's five and the run's `rama_round_04.dat`
  (md5 `7693f10a`), in `$P/parameters`, the shared `/beagle3/.../upside2-md/parameters`, the Upside 3 tree's
  `parameters` and `ff3_selection/params`, md5-verified, and in the Mac repo. The GPU tree (`upside3_gpu`, line 5) was
  not touched. Validation, the merge's tests:
  - panels `p44_ff31nf` and `mdc_ff31nf` (dt 0.009) COMPLETED 18:02 and 23:09 10-31, read 23:40 (`paired_npz.py`;
    findings 1.37): 44 domains fold 0.631 against released ff2.1's 0.603 (+0.028 [-0.004, +0.059] p 0.047), mean TM
    +0.001 [-0.006, +0.008]; mdCATH's 89 fold 0.515 against 0.488 (+0.027 [+0.000, +0.053] p 0.079); level with the
    merge on both, +0.20 / +0.16 folded over the control. TM files `tm_p44_*_ff31nf_vs_*.txt` in the check directory
    (49692208) and `tm_mdcath_*` (49692207 COMPLETED 00:23 11-01, read 01:30): mdCATH mean TM 0.731 against
    released ff2.1's 0.727 (+0.004 [-0.003, +0.012] p 0.22, higher on 52 of 89), the merge's 0.733 (-0.001), the
    control's 0.702 (+0.030 [+0.022, +0.037]);
  - Peng: 32 arms `b_ff_3.1_<prot>_<kind>` **49692209-49692240**, each resubmitting itself (runs
    `ff3_benchmark/runs/<prot>_<kind>_ff_3.1/`, logs `ff3_benchmark/logs/<prot>_<kind>_ff_3.1_[0-9]*.out`); collector
    `ff31_collect` (every 6 h; scores into `score/ff_3.1/`, writes `peng_ff31.txt` in the check directory, ending
    DONE; 6 of 32 scored after pass 49748728, 03:49 11-03, next **49750446**). Two chunks lost their node and were
    requeued by Slurm, each resuming from its checkpoint: cspA de novo 49692216 (midway2-0189, 20:11 10-31) and alpha3d
    native 49712331 (midway2-0156, 06:11 11-01). Second chunks ended 12:52-12:56 11-01, mean 0.38 of the target done
    (0.12-0.89). BBA native and BBL native COMPLETE 15:41 and 15:51 11-01 and scored; final frames read 12:17 11-02:
    finite, KE/1.5kT 1.002-1.008 over 14 replicas each. **30 chains broke on midway2's fault**, 28 at 23:22-23:28 11-01
    and alpha3d native 49738144 and cspA de novo 49738143 (offset by their requeues) at 03:16 11-02: each chunk ended
    normally (`exit=0`) and its resubmission failed ("sbatch: error: Batch job submission failed: Unexpected message
    received"). **Restarted 12:18 11-02** (user: "Restart 30 Peng chains") by
    `checks/ff31_nativefix_release_20261031/restart_broken.sh submit` (md5 `c77df969`; its dry run listed the 30 and
    no other): **49743748-49743779**, ids in `restart_broken.log`, 27 R and 3 PD at 12:20; each arm resumes from its
    checkpoint and resubmits itself as before;
  - glpG on Upside 3: `upside3_validation/glpg_ff31/` (the merge's `run.py`, one REMD block). Prep COMPLETED 20:57
    10-31; `glpg_u3_prep_tm4.py` (`checks/ff31_nativefix_release_20261031/glpg_prep_tm4.txt`): TM4 DSSP alpha through
    stage 7.0 by tenth 1.00 0.98 1.00 0.99 1.00 0.90 0.91 0.91 0.93 0.87 (mean 0.95, last frame 0.82, box 1.00), where
    the merge fell to 0.36 and ff3.0 settled near 0.82 (findings 3.15). REMD block 1 49692243 ended 07:19 11-02 (log
    `REMD/glpG-RKRK-79HIS.remd_block1.log` ends "finished", final KE/1.5kT 0.996-1.003, no NaN, ROLLBACK or error
    line; `block_count` 1 of 1, so the run is done; Slurm stamped it TIMEOUT at 09:43, the controller's return). **Read 07:30 11-02** on the
    midway3 login (`glpg_u3_block1.py` on `/output` of run.0 and run.23; `glpg_ff31/glpg_ff31_block1.txt`): at T 0.70
    KE/1.5kT protein 1.004, lipid 1.005; C-N 1.312 +- 0.120 A; TM1 dihedral-box 0.98; **TM4 DSSP alpha 0.74** (halves
    0.78 / 0.70, last quarter 0.64; 135 at 0.34, 146-151 at 0.48-0.70), at T 0.795 0.71; ff3.0 0.85, the merge 0.42.
- **The 10-26 merge, retired 16:49 10-31** (plan.md Phase 30; `rama.dat` `5b66c78a`). Its results stand as recorded:
  level with released ff2.1 on both panels (p44 mean TM +0.006 [-0.001, +0.012], mdCATH +0.005 [-0.001, +0.012]),
  above the control on both; Peng 25 of 32 scored when its last 5 chunks and collector were cancelled 16:47 10-31; on
  Upside 3's glpG it lost TM4 in block 1 (DSSP 0.42; line 3). Kept: `<tree>/backup/ff_3.1.bak_pre_nativefix_20261031`,
  `ff3_selection/params_superseded/bak_pre_nativefix_20261031/ff_3.1`, its Peng runs and logs in
  `ff3_benchmark/{runs,logs}_superseded/ff_3.1_bak_pre_nativefix_20261031/` (replicas deleted 11-03), scores
  `score/ff_3.1.bak_pre_nativefix_20261031`, glpG `upside3_validation/glpg_ff31.bak_pre_nativefix_20261031` (REMD
  replicas deleted 11-03; "Disk"), its panels
  `p44_ff31` / `mdc_ff31` and its scripts in `checks/ff31_merged_20261026/`. A log glob `<prot>_<kind>_ff_3.1_*.out`
  also matches the retired `ff_3.1_depth_data` logs; use `ff_3.1_[0-9]*.out`.
- **The fix of the all-term loss, step 1** (user 10-31: "Fix the folding issue and continue working on ff3.2"; plan.md
  Phase 38; findings 1.37): `ff32_dse0` (dse_weight 0) **49692528** R from 16:57 and `ff32_sarwsi` (the SI's SARW,
  backbone burial off) **49692529** R from 17:02 10-31 (steps 9 and 6 of 19 at 23:25, no error line), one epoch each, built as `ff32_nativefix` with one trainer line
  changed (`checks/ff32_fix_20261031/`: `build.sh`, `build.log`, `compare_init.py`, `sarw_check/`, `submitted.log`).
  Both reached 19 (03:42 and 10:12 11-01; 19 minibatches each, no error line); `after_training.sbatch` ran each
  44-domain panel (`z32_00`, `s32_00`, dt 0.015) and `peng_full.sbatch` (`peng_full_<tag>.txt` in the check
  directory). **Read 14:10 11-01** (`paired_npz.py`; findings 1.37): `z32_00` (DSE off) folds 0.575, +0.127 [+0.091,
  +0.163] over `n32_00` and -0.040 [-0.076, -0.007] p 0.012 under released ff2.1, unfolded shift +3.80 E_up (2 of 16
  negative; `n32_00` -11.9, 16 of 16); `s32_00` (the SI's SARW) folds 0.446, level with `n32_00` (-0.003 [-0.036,
  +0.029]), shift -9.93 (16 of 16). Logs `$P/training/<run>/condiv-train_<id>.out`.
- **ff3.2 is `ff32_dse0` trained on** (user 11-01: "Proceed to continue work on ff3.2. Double check it is training
  all parameters, not holding anything fixed"; plan.md Phase 38 items 4-6). Every group the engine gives a gradient
  for moved in epoch 0 (`checks/ff32_fix_20261031/params_check.txt`; `param_moves.py`, `ckpt_rates.py`). From step 19
  at 14:43 11-01, one epoch a link: link **49737476** FAILED 21:19 11-01 after 13 of epoch 1's 19 minibatches (step
  32, minibatch 13), when every worker's srun was refused on midway2's fault ("Unable to confirm allocation for job
  49737476: Unexpected message received"; `run_output/epoch_01_minibatch_13/*.srun`); the 10-28 failure was the same
  kind. Insurance **49737477** (midway2-[0316-0327]) resumed from minibatch 12's checkpoint 09:45 11-02 and
  COMPLETED at the epoch end, step 38, 16:47 11-02 (minibatches 13-18 at 3,841-4,581 s, 24 of 24 proteins each, KE/1.5kT
  0.994-1.014, no error line; H-bond margin +0.131 at step 19, +0.107 at 29, +0.126 at 38, above the +0.10 hold); at
  its target it cancelled its own insurance 49738281. Epoch 2's link 49748334 ran from 16:49 11-02, again on
  midway2-[0316-0327], until it was cancelled 18:12 (below). At each epoch
  end the run's `after_training.sbatch` (md5 `11bd82e9` since 12:33 11-02, log
  `$P/training/ff32_dse0/after_training_<id>.out`) submits panel `z32_EE` (dt 0.015, `ff3_selection/runs/z32_EE/`)
  and `peng_full.sbatch` on it (`peng_full_z32_EE.txt`), rewrites `rg_check_ff32_dse0.txt` (free-replica Rg against
  `ff31_nativefix`), then the next epoch; at step 76 it submits `z32_03` and stops, with no convergence gate (user
  11-02; the gated version kept as `.bak_pre_budget_20261102`). Its submissions are logged in
  `checks/ff32_fix_20261031/submitted.log`. Since 09:45 11-02 a step takes about 70 min against 30-37 before the
  fault, each worker's simulation 2.4 times longer (median wall 1,192 s at epoch 0 minibatch 13, 2,885 s at epoch
  1's). **The cause is those nodes' clock** (read only, `ssh` to the job's nodes: `/proc/cpuinfo` MHz on the cores
  `taskset` gives each `upside` worker; probe `node_clock.sh`, in the bundle's `scratchpad/handoff_20261102/`). At
  18:08 11-02, both workers running on each, every busy core of midway2-0317 to 0327 read 1,199 MHz, the minimum
  P-state, and 0316's 991 and 1,562, against 2,500-2,700 on `ff32_basin`'s nodes (0248, 0262) in the same minutes,
  with the same governor (`powersave`), `max_perf_pct` 100 and turbo on; the ratio is the workers' 2.4. The link
  before the fault on 12 other nodes (49737476, midway2-[0170-0171,0365-0374]) ran 1,553-2,264 s a step, so the
  12-node layout is not the cause. **At the user's word (18:10 11-02) midway2-[0316-0327] are excluded** (§0d;
  `midway2.args` and its README row, backups `*.bak_pre_0316_20261102`; `--test-only`: the wrapper refuses 0320 and
  0327 and accepts 0328, the real sbatch accepts 0320; `update_pending.sh` set 49748335, 49747516 and 49748039),
  and epoch 2's link 49748334 was cancelled 18:12:12 after step 39's checkpoint (18:05, 4,594 s; it loads), 5 min
  into step 40. Its insurance **49748335** waited 59 min and runs from 19:11 on 22 nodes outside the list
  (midway2-[0026,0028,0030,0096,0098-0100,0113,0152,0156,0237-0239,0254,0270,0272,0295,0365,0376,0379,0384,0388]),
  resumed at step 39 for 18 steps; at 19:20 every worker's cores read 1,961-2,777 MHz at 11-12 of 14 cores busy.
  It ran step 40 (1,768 s, checkpoint 19:41) and was ended 19:49:22, 8 min into step 41, by a node failure
  (NODE_FAIL; its log: "CANCELLED ... DUE TO NODE FAILURE"; no node of its list drained at 20:35). Its insurance
  **49748605** started 19:49:31 on midway2-[0026,0028,0098-0100,0102,0112-0113,0160-0161,0170-0171], resumed at
  step 40 from `epoch_02_minibatch_01/checkpoint.pkl` (`.chain_starts` 38, 39, 40) and ran step 41 in 1,681 s; its
  successor **49748621** PD afterany carries the list. The ff3.2 jobs create no file on `/beagle3`, whose group
  inode quota is spent since about 19:30 ("Owed by hand"): training, panels and Rg checks write under `/project`,
  and `peng_full.py` only reads released ff2.1's frames from `/beagle3/trsosnic/yinhan/ff3_benchmark/runs`.
  `/project` is 480 G under its soft quota since the 11-03 reclaim ("Disk"). **Step 57 reached 05:50 11-03**
  (steps 48-57 31-43 min each); its `after_training.sbatch` queued panel **`z32_02` 49754402** and shift
  **49754413** (afterok) at 05:50 and epoch 3's first link **49754482** (target 76) at 05:52, and wrote the Rg check:
  at T 1.1 0.923 [0.868, 0.964] of `ff31_nativefix`'s, below 1 in 89% (epoch 1 0.939), 0.937 of its own epoch 0; at
  T 0.8 1.004. The panel and shift COMPLETED 12:12 and 12:16; the link R from 11:41 on 16 nodes, step 69 at 17:37
  (epoch 3 at 26-32 min a step, 24 of 24 proteins, KE/1.5kT at most 1.019, margin +0.131 at step 57 rising to
  +0.159 at 69), successor 49765255 PD afterany; step 76 about 21:05. **`z32_02` read 17:50 11-03**
  (`paired_npz.py`): folds 0.599, against released ff2.1 -0.016 [-0.059, +0.026] p 0.33 (25 lower, 18 higher), the
  first of this run's epochs level with it; against `z32_00` +0.024 [-0.016, +0.065] p 0.41, `z32_01` +0.021
  [-0.024, +0.066], `w32_00` -0.031 [-0.069, +0.007] p 0.088. Unfolded shift +6.38 E_up, 0 of 16 negative. Since 11-03 this run is the comparison for ff3.2 (`ff32_basin`,
  below); its epoch choice (plan.md Phase 38 item 7) is read beside ff3.2's, and the full validation, about 3 days
  (Peng's 32 arms in 10.5 h chunks over 2.5-3 days, glpG on Upside 3 about 4 h of preparation and 34 h of REMD,
  mdCATH about 6 h), goes to the epoch chosen for ff3.2 if it is level with released ff2.1.
  - **Epoch 0's glycine and folding, read 12:45 11-02** (plan.md Phase 38 item 7): the pooled glycine map's alpha_L -
    alpha_R -1.064 (start) to -0.877 (ff3.0 -0.836, released ff3.1 -0.781); `z32_00` folds 0.575 and keeps helical and
    left-handed glycines in their native basins level with released ff2.1.
  - **The remaining loss** (user 11-02: "find the source of the remaining loss and try to resolve it"; plan.md Phase 38
    items 8-10; `/project/trsosnic/yinhan/checks/ff32_rest_20261102/`, ids in `submitted.log`): the read-only jobs
    (`diagnose` 49744922, `decomp3` 49745312, `crystal_view` 49745333) COMPLETED 13:05-13:15 and are read (findings
    1.37 "The remaining loss"). Six one-file reset panels `z32_00_{rot21,bb21,hb21,env21,sheet21,rama0}`
    49744916-49744921 COMPLETED about 14:10 (176 npz each, `ff3_selection/runs/<tag>/`), read 14:20
    (`checks/ff32_rest_20261102/resets_z32_00.txt`; findings 1.37): side chain +0.054, burial +0.041, backbone burial
    -0.096 against `z32_00`.
  - **Epoch 1, read 17:40 11-02** (the Mac Studio; plan.md Phase 38 item 6; after-training job 49748329, panel
    49748331 and shift 49748332 COMPLETED 16:49-17:36): `z32_01` folds 0.578, against released ff2.1 -0.038 [-0.072,
    -0.005] p 0.038 (29 lower, 14 higher), against `z32_00` +0.002 [-0.033, +0.039] p 0.94, against `n32_01` +0.163.
    Unfolded shift +5.61 E_up (1 of 16 negative; `z32_00` +3.80). Rg at T 1.1 against `ff31_nativefix` 0.939
    [0.887, 0.981] (epoch 0 0.978), against its own epoch 0 0.957 [0.907, 0.995], below 1 in 77%; at T 0.8 1.004,
    level. Round 2's glycine map: alpha_L - alpha_R -0.791 (round 1 -0.877, released ff3.1 -0.781),
    -ln(P_aL/P_aR) -0.554 (round 1 -0.640).
  - **The fix's test arm `ff32_basin`** (user 11-02 chose the native-basin restraint; plan.md Phase 38 item 10): the
    restrained replica held by a one-sided wall at 1.2 d0 on every native CA contact instead of the 1.0 A elastic
    network, DSE off, otherwise `ff32_dse0`; one epoch from the start. First link **49746562** submitted 13:54, R
    from about 14:25 11-02 on 24 nodes (midway2-[0248-0249,0258,0261-0264,0266-0269,0283-0290,0368-0372]); minibatches
    0-4 done by 17:13 (1748, 1693, 1728, 2013, 1618 s, no error line; step 5: 24 of 24 proteins, KE/1.5kT max
    1.014, margin +0.188); steps 5-18 1,783-2,386 s; **COMPLETED 00:48:37 11-03** (19 of 19 steps with 24 of 24
    proteins, no error line; log `$P/training/ff32_basin/condiv-train_49746562.out`), insurance 49747516 cancelled by
    the chain. Its `after_training.sbatch` (md5 `a5ffaef1`) 49748993 COMPLETED 00:50: panel `w32_00` 49748994
    COMPLETED 02:36 (176 npz, `ff3_selection/runs/w32_00/`), shift 49748995 COMPLETED 02:51
    (`checks/ff32_fix_20261031/peng_full_w32_00.txt`), and `checks/ff32_rest_20261102/rg_check_ff32_basin.txt`: Rg at
    T 1.1 against `ff31_nativefix` 0.943 [0.869, 0.985], below 1 in 84% (`z32_00` 0.978, `z32_01` 0.939); at T 0.8
    0.996. Step 19 by `check_step.py`: 24 of 24 proteins, KE/1.5kT 0.998-1.013, margin +0.106 (`z32_00` +0.131); the
    restrained replica's median RMSD 1.23-1.53 A over the 19 steps (`ff32_dse0`'s epoch 0 0.85-0.91; free 2.1-3.0).
    Item 10's remaining readouts (`checks/ff32_rest_20261102/submitted.log`): `mech_w32_00.sbatch` (md5 `523447a3`)
    49749015 COMPLETED 02:48 (`decomp2_w32_00.txt`, `crystal_view_w32_00.txt`; its first submission 49749009 FAILED
    in 82 s on the half-written panel, which both scripts read) and `tm_panel.sbatch domains runs w32_00 z32_00`
    49749014 COMPLETED 02:56 (`checks/converge_20261020/tm_panel_49749014.out`; series cached in `tm_runs/w32_00/`).
    **Read 03:05 11-03 by item 10's rule: (a) and (b) both met** (plan.md Phase 38 item 10; findings 1.37):
    crystal-minus-native-like change -1.76 E_up and loose -0.40 on the 32 domains with 20 loose frames in `w32_00`
    (`z32_00` on the same 32: -3.69 and -1.13; smaller in 32 and 31 of 32); on 3g0m's fixed frames the 3-7 A dip of
    `z32_00` (-0.71 to -0.82 against the restrained mean) is gone (+0.75 to +2.17, rising with RMSD). Folded 0.630
    against `z32_00` 0.575, +0.055 [+0.022, +0.091] p 0.0039 (30 higher, 13 lower), against `ff21_released` +0.015
    [-0.015, +0.045], against `p44_ff31nf` -0.001; mean TM 0.684 against `z32_00` 0.679, +0.006 [-0.002, +0.013]
    p 0.23 (`p44_ff21_released` 0.681, `p44_ff31nf` 0.682). Unfolded shift +6.99 E_up, 0 of 16 negative (`z32_00`
    +3.80). Reported to the user, who decides.
    **Adopted (user 11-03, "Ok. Proceed"): ff3.2 is this run trained on to 76** (plan.md Phase 38 item 11), driver
    and minibatch order unchanged. `after_training.sbatch` is now `ff32_dse0`'s per-epoch form (md5 `451ae243`; the
    one-epoch `a5ffaef1` kept as `after_training.sbatch.bak_pre_epochs_20261103`, the same pair in
    `checks/ff32_rest_20261102/` as `after_training_ff32_basin.sbatch`): at each epoch end it queues panel `w32_EE`,
    `peng_full.sbatch` and `mech.sbatch` (md5 `46766551`; `decomp2.py`, `crystal_view.py` into `decomp2_w32_EE.txt`,
    `crystal_view_w32_EE.txt`, log `logs/mech_<id>.out`) afterok it, rewrites `rg_check_ff32_basin.txt`, then the
    next epoch until 76. Checked by a stub run of every branch (`stub_epochs.sh`, `stub_epochs.txt`: step 38 queues
    `w32_01` and trains on to 57, a rerun skips the panel, step 42 is refused, step 76 stops) and `--test-only`; no
    test job. First link **49757705** (`train_chain.sbatch` target 38, from `epoch_00_minibatch_18`) submitted
    07:39:28 11-03, PD on priority; at 08:05 `ff32_dse0`'s epoch-3 link 49754482, 74 points above it from age alone, was set to nice 100 (user: "Decide for me"), so this link takes the first 12-node slot (`submitted.log`). R from 11:28 11-03 on 23 nodes,
    successor 49764882 PD afterany. Read 17:45 by `check_step.py`: step 31 at 17:33, steps 20-31 at 1,603-1,989 s
    (mean 30 min), 24 of 24 proteins each, KE/1.5kT 0.994-1.015, restrained median RMSD 1.21-1.50 A; glycines over
    epoch 1's 12 steps, free, helical in alpha_R 0.855, left in alpha_L 0.894. **H-bond margin** +0.106 at step 19,
    +0.098 to +0.106 since, under +0.10 at steps 24-26 (lowest +0.0977) and 31 (+0.0997); dhb -0.384 to -0.375,
    sheet mean 0.079 to 0.061. The user was notified 17:50 11-03; the run trains on, as the earlier all-term runs did
    under the line (user 10-28), and the margin does not set the folded loss (findings 1.21). At ~30 min a step,
    step 38 about 21:05 11-03; step 76 at 16:00 11-04 with no queue wait, about 01:00 11-05 with two waits like
    this morning's 3.8 h (broadwl full, 79 nodes of one-node jobs pending above our priority at 18:35), 07:00 at 6 h.
  `ff32_nativefix` (`n32_EE`, with the DSE) stays stopped as the comparison.
- **ff3.1.1 and ff3.2.1, finished and merged** (user 10-21: "apply the fix and train ff3.1.1 and ff3.2.1", then "try
  to resolve the accuracy lost bug as well in this run"; plan.md Phases 29-30; findings 1.34, 1.37):
  `ff311_depths_only` (`p311_EE`, from NDRD's map) and `ff321_depths_only` (`p321_EE`, from ff3.0's map), ff2.1 with only
  its depths trained. Built by `checks/ff311_20261021/build_depths_only.sh` (logs `build_ff311_depths_only.log`,
  `build_ff321_depths_only.log`): ff2.1's files; `trainer/rama_basin.py` the repo's with MODE `pooled` (md5
  `54816673...`, shared width 0.5); every Adam rate zero (`freeze_groups.py`; the trained rates kept as
  `initial_checkpoint.pkl.trained_rates`); ff3.1's rerun's minibatch order and held-out set. Steps 0-113 done 04:56
  (ff3.1.1) and 17:19 10-25 (ff3.2.1); step 76's library is each run's `run_output/rama_round_04.dat`, step 114's
  `rama_round_06.dat`. `submit_new.sh` does not list them (backup `.bak_pre_phase30_20261023`).
  - Collector `pooled_collect` (`checks/ff311_20261021/collect_pooled.sbatch`, md5 `560cc098`; `status.txt`, `logs/`)
    ended DONE: 114 steps each, 24 of 24 proteins at every step, KE/1.5kT at most 1.019 and 1.018. Steric flags (free
    non-glycine beta' + pPII' above each run's steps 0-18 maximum): ff3.1.1 (0.0043) steps 33, 39, 71, 77, 90, 96, 109;
    ff3.2.1 (0.0050) step 90 (0.0051); each minibatch 1 or 14, which read as high in every epoch and in the restrained
    runs (findings 1.33). No clash flag. Its first version failed on `set -eo pipefail` and a `while` pipeline (findings
    10.14).
  - Glycine's alpha_L shared offset, common frame, by round 0-6: ff3.1.1 0, +0.157, +0.217, +0.293, +0.358, +0.367,
    +0.386; ff3.2.1 +0.517, +0.555, +0.556, +0.549, +0.537, +0.504, +0.458 (-0.059 + ff3.0's +0.517); gap 0.179 at step
    76, 0.072 at round 6. In `rama_rounds.txt` "|c| > 0.1: N of M" counts every (map, basin) offset over M trainable
    depths, so N can exceed M (ff3.1.1 round 6: 1726 of 1517).
  - Agreement check (10-26; plan.md Phase 30 item 2): `checks/ff311_20261021/dump_rama_state.py` (md5 `72306af1`) wrote
    rounds 0-6 to `rama_state_pooled.npz`, and `compare_pooled_runs.py` (`48567a34`, run on the Mac from
    `scratchpad/ff31_compare_r6/`) wrote `compare_pooled_runs.txt`: every term agrees to the runs' noise except
    glycine's alpha_L; at round 6 (a) met, (b) not, short only by a constant the data cannot move. The user chose the
    merge.
- **ff3.1 and ff3.2 on the fixed native target** (user 10-28: "try a version with everything trainable, start with
  current ff3.1, call it ff3.2. Fix ff3.1 according to your plan"; plan.md Phase 33; findings 1.34). The runs' own
  `ConDiv.py` (md5 `54383004`) adds a 20 deg phi/psi spring per residue to the native-restrained replica
  (`add_native_rama_restraint`); both start from the released ff3.1 (round-0 library md5 `5b66c78a`) in ff3.1.1's
  minibatch order. `ff31_nativefix` (ff3.1: every Adam group held, every depth trainable) **49635325** and
  `ff32_nativefix` (ff3.2: every group at the trainer's rates) **49635326**, submitted 02:26 10-28 to target 76 through
  `train_chain.sbatch` (15 nodes a link; Slurm's start estimate 10:23 10-30). Run dirs `$P/training/<run>`, logs
  `<run>/condiv-train_<id>.out`. Scripts, `build.log` and the static check `/project/trsosnic/yinhan/checks/ff31_fix_20261028/`
  (`ConDiv.py`, `check_restraint.py`, `set_start.py`, `compare_start.py`, `native_fidelity.py`, `submitted.log`).
  **Gate read 10:50 10-28: FAIL as written** (glycine beta' natives 0.656 in ff3.2, 0.660 in ff3.1, need 0.75; alpha_R
  0.998), but an ideal 20 deg spring gives those glycines 0.726 at most (`gate_ceiling.py`; findings 1.34), so the
  threshold cannot be met. Both first links FAILED 10:41 10-28 when every worker's srun lost the Slurm controller
  ("Unable to contact slurm controller"), ff3.2 in step 11 and ff3.1 in step 16, from checkpoints of steps 10 and 15;
  the successor links 49635669 (ff3.2) and 49635327 (ff3.1, nice 10) were held, then released 10:58 10-28 on the
  user's "proceed with your best decision" (plan.md Phase 33 item 4). ff3.2 has priority over ff3.1 and over the
  debug chains (user 10-28; those at nice 10): when a link starts, give the insurance link it submits the same nice.
  Then read each round end as plan.md Phase 33 lists.
  **At 22:10 10-28:** ff3.2 49635669 R from 11:12 10-28 (insurance 49636501, nice 0, PD afterany), in step 32, about
  28 min a step, 24 of 24 proteins every step; round 1 read (plan.md Phase 33 item 6): glycine's pair reached ff3.0's,
  and the **H-bond margin is +0.059 at step 31, under the +0.10 line** (notified the user 22:10 10-28; it keeps
  training; user 10-28 "proceed with your suggestion"). Its epoch-end panels `n32_EE` now come from `submit_new.sh`
  (backup `.bak_pre_n32_20261028`). `n32_00` (49638827, 23:06 10-28) folds 0.449, -0.167 [-0.206, -0.129] against
  released ff2.1 and level with `a32_00`: the fix does not cure epoch 0's loss (findings 1.37).
  **At 09:30 10-29:** ff3.2 in step 56 (epoch 2 minibatch 17), 24 of 24 proteins every finished step, no error line,
  KE/1.5kT at most 1.016 (steps 33-55). Round 2 read (plan.md Phase 33 item 6; `rama_round_02.dat` md5 `6c9ce4f1`).
  Its margin (from the log's `hb` line, E_other - E_alpha) came back from +0.059 (step 31) to +0.080 at epoch 1's end,
  fell to +0.036 (step 50) and is +0.051 at step 54. The epoch-1 panel `n32_01` **49640339** was submitted 09:22 10-29
  by `submit_new.sh` (PD, the law's exclusions carried). ff3.1 49635327 PD (Resources, nice 10, Slurm's estimate
  15:09 10-29), still at step 16.
  **At 12:30 10-29:** ff3.2 in step 62 (epoch 3 minibatch 5), 24 of 24 proteins every finished step, no error
  line, KE/1.5kT at most 1.015; round 3 read (plan.md Phase 33 item 6; `rama_round_03.dat` md5 `67dd261c`); margin
  +0.058 at epoch 2's end and +0.061 to +0.067 over steps 57-61. `n32_01` and ff3.1 still PD (Priority).
  **At 14:30 10-29:** ff3.2's link 49635669 ended COMPLETED 13:58 10-29 after its 54 steps (11-64): 65 of 76 steps done,
  24 of 24 proteins every step (two srun launch failures relaunched in step 22), KE/1.5kT at most 1.012 over steps
  60-64, no error line. Margin (hb[2] - hb[0]) +0.061 / +0.066 / +0.067 / +0.066 / +0.064 / +0.062 / +0.061 / +0.060
  over steps 57-64, dhb -0.556 -> -0.574; E_beta - E_alpha -0.042 at step 64, from -0.061 at step 58.
  **At 17:10 10-30: ff3.2 is at its target, step 76, and stops for the user** (no gate; no link of the run queued).
  Its last link 49636501 (R 04:36-11:58 10-30, after ff3.1's and the debug's links, which started 02:03 and 03:17)
  ran steps 65-75 and ended `target 76 reached`: 24 of 24 proteins every step, KE/1.5kT at most 1.013, no error line.
  Margin +0.058 (step 65), +0.042 (step 69), +0.065 (step 75); dhb -0.575 -> -0.605 and E_beta - E_alpha -0.036 ->
  -0.071 over those steps. Round 4 read (plan.md Phase 33 item 6; `rama_round_04.dat` md5 `96e265d2`). Its
  convergence gate, run by hand (`checks/ff31_fix_20261028/gate_step76_ff32_nativefix.txt`): NOT CONVERGED, only
  `dhb` pulled (p 0.0045, cut 0.0050). The epoch-1
  panel `n32_01` (49640339, ran 02:14-03:05 10-30) read (Phase 33 item 7): folded 0.414 against released ff2.1's
  0.615.
  **Extended to step 95 17:45 10-30 and finished there** (user 10-30: "Go ahead and continue training ff3.2"; plan.md
  Phase 33 item 8): link 49670425 (`train_chain.sbatch ff32_nativefix 95` from the run dir, 12 nodes, nice 0; logged in
  `checks/ff31_fix_20261028/submitted.log`) ran steps 76-94 18:09 10-30 to 03:10 10-31 and ended `target 95 reached`:
  24 of 24 proteins every step, KE/1.5kT at most 1.013, no error line. Margin +0.064 (step 76) rising to +0.095 (step
  94), still under the +0.10 line; dhb fell every step, -0.608 -> -0.647; E_beta - E_alpha -0.080 -> -0.076. Round 5
  read 14:05 10-31 (plan.md Phase 33 item 6; `rama_round_05.dat` md5 `ad998fe5`). **Gate at step 95, run by hand
  (`checks/ff31_fix_20261028/gate_step95_ff32_nativefix.txt`): CONVERGED**, dhb p 0.027 against the cut 0.005, every
  other group 0.71-1.00; by item 8's rule ff3.2 is not extended to 114. Panels `n32_02` (49670431) and `n32_03`
  (49670432) ran 18:09-18:53 10-30 and are read (Phase 33 item 7). `submit_new.sh` (dry-run first: only `n32_04`)
  submitted the epoch-4 panel `n32_04` 49686492 14:00 10-31; COMPLETED 14:46, no error line: folds 0.447, -0.169
  [-0.215, -0.126] against released ff2.1, so ff3.2 stops (Phase 33 item 8). The diagnosis of the loss (user 10-31:
  "Check what is wrong in ff3.2 that fails to fold protein"; findings 1.37, Phase 33 item 9) is read only, scripts and
  outputs in `/project/trsosnic/yinhan/checks/ff32_fold_20261031/`; its one job, the valgrind check 49687255
  (`repro_pd.sbatch`), COMPLETED 14:39 10-31, no engine error under memcheck.
  ff3.1 (`ff31_nativefix`): reached its target 76 at 11:20 10-31 (links 49635327 to 08:07 and insurance 49662306,
  24 of 24 proteins every step, no error line; its H-bond terms are held, so its margin stays +0.192); no gate, stops
  for the user. Round 4 library md5 `7693f10a`: glycine shared aR -0.171, aL +0.495; -ln(P_aL / P_aR) -0.545 (ff3.0
  -0.598). Rounds 1-4 read 10-31 (plan.md Phase 33 item 6).
- **The ff3.1 training debug** (user 10-27: "setup your debug for ff3.1 training"; plan.md Phase 32; findings 1.34):
  ff3.1.1 as built (`build_depths_only.sh`, every Adam group held) with glycine parameterised as ff3.0 is, its shared
  alpha_R and alpha_L only (`restrict_gly.py`); expected, if the training is sound, to reach ff3.0's pair (alpha_R
  -0.107, alpha_L +0.517). `ff311_glypair` (ff3.1.1's minibatch order) **49634971** and `ff311_glypair_mb2` (its own
  order, same 456 proteins and held-out set) **49634972**, submitted 23:03 10-27 to target 76 through
  `train_chain.sbatch` (15 nodes a link). `ff311_glypair` R from 23:05 10-27 (successor 49634973 afterany, the chain's
  insurance link); `_mb2` R from ~23:30 10-27 (successor 49635107). Both first links died in the 10:41 10-28 outage;
  49634973 and 49635107 R from 10:58 10-28 (insurance 49636476, 49636479 PD afterany, nice 10). At 22:10 10-28: 33 and
  32 steps done, 24 of 24 proteins every step, no error line; round 1 read (plan.md Phase 32 item 3). At 09:30 10-29
  `_mb2` is in step 56 (epoch 2 minibatch 17), 24 of 24 proteins every finished step, KE/1.5kT at most 1.013, no
  error line; round 2 read (Phase 32 item 3; `rama_round_02.dat` md5 `51e4dc72`). At 12:30 10-29 `_mb2` is in
  step 61 (epoch 3 minibatch 4), 24 of 24 proteins every step, KE/1.5kT at most 1.013, no error line; round 3 read
  (`rama_round_03.dat` md5 `f65a593e`).
  **`ff311_glypair` hung in step 34 from 22:41 10-28** (epoch 1 minibatch 14, started 21:47). 23 of 24 workers
  wrote their `divergence.pkl` by 22:41; 3g0m's (srun step 49634973.301 on midway2-0391) ran its simulation to
  22:14, and then one of `compute_divergence`'s seven Pool children printed glibc `double free or corruption (out)`
  in `RotamerSidechain<PosQuadSplineInteraction>`'s destructor, called by `free_deriv_engine` from an `Upside`
  object's finalizer (`epoch_01_minibatch_14/3g0m.output_worker`). Pool replaced the dead child (pid 12942, 51 s
  younger than the other six) and `pool.map` waited for the lost task, all eight processes in futex waits at 0.0%
  CPU (§8). 3g0m's epoch-0 step passed; no other worker of the four runs printed the error. Cause not identified (an
  earlier heap overrun in that child is the reading glibc's message supports; ConDiv's proteins are under 1,000
  residues, so A26's intensive-memory overflow does not apply). The same message, in the same destructor, came once
  on the midway2 login node 10-31 from a read-only script that called `get_param_deriv` on many engines and freed
  them (`checks/ff32_fold_20261031/grad_split.py`, pre-audit `$P/obj`); a memcheck of 4 configs with the same calls
  (49687255, `repro_pd.py`) found no invalid access in the engine, so it is still not localised. **Cancelled 09:20 10-29** (user 10-29 "fix them
  all"): 3g0m's 19 step files (base, `initial.npy`, 14 replica `.up`, `run.0.up.output`, `output_worker`, `srun`) and
  the link's log were copied first, md5-checked, to `checks/ff31_debug_20261027/hang_3g0m_step34/` (2.1 G); then
  `scancel 49634973.301`. The master logged `3g0m WORKER_FAIL` and `1 of 24 workers failed`, the link ended FAILED
  09:21, and its insurance **49636476** (nice 10) is PD (Priority) to rerun step 34 from step 33's checkpoint in the
  same directory. **The rerun did not repeat it:** 49636476 R from 03:17 10-30 reran step 34, and 3g0m's
  `output_worker` is empty with its `divergence.pkl` written 03:45. That link then lost its node (NODE_FAIL,
  midway2-0097, 21:56 10-30, in step 65); its insurance 49662466 (R 21:58 10-30) reran step 65 from step 64's
  checkpoint and reached the target 76 at 03:35 10-31, 24 of 24 proteins every step, no error line. Round 4:
  `rama_round_04.dat` md5 `15087283`, glycine pair alpha_R +0.069, alpha_L +0.456, -ln(P_aL / P_aR) -0.834.
  `_mb2`'s 49636479 (R 04:46-11:44 10-30) reached its target 76: round 4's glycine pair
  alpha_R +0.012, alpha_L +0.388 (ff3.0's -0.107 / +0.517), epoch-3 pull t -0.2 and +0.6; `rama_round_04.dat` md5
  `c4f2ac9d`, -ln(P_aL / P_aR) -0.844 (round 3 -0.863). Both arms stop for the user. Run dirs `$P/training/<run>`, logs
  `<run>/condiv-train_<id>.out`. Scripts, build log and the initial-checkpoint check
  `/project/trsosnic/yinhan/checks/ff31_debug_20261027/` (`build.sh`, `restrict_gly.py`, `compare_init.py`,
  `compare_init.txt`, `pooled_pull.py`, `submitted.log`). No collector: read each round end by hand with
  `pooled_pull.py <run> <epoch>` and the repo's `scratchpad/ff31_gly_pull_20261027/gly_LR.py` on `rama_round_NN.dat`,
  against plan.md Phase 32's rule; check each step's "24 of 24 proteins" and the log's error lines.
- **Finished, awaiting the user's review:** the fixed-engine reruns `ff31_depth_data_angfix` (`a31_EE`) and
  `ff32_depth_data_angfix` (`a32_EE`), step 76 10-21: margins +0.109 and -0.002 at step 75; a `STERIC FLAG` on both
  from the free beta' + pPII' alone (0.0052, 0.0053 against steps 0-18's 0.0045), clash ratios small; the ff3.2 rerun
  at a fixed point by every test (findings 1.36); their 44-domain panels at step 76 0.420 and 0.432 against ff2.1's
  0.615 (findings 1.37).
- **The 8 corrected glpG chains on Upside 2, finished 10-24** (user 10-17: "fix the bug and rerun the affected
  simulations"; plan.md Phase 18; findings 4.17): the side-chain placement fixed in the preparation and every live
  seed, 28 replicas T 0.70-0.90, five 36 h blocks (`submit_remd.sh`), ff_3.0_gdepth and the old ff3.1 x {79HIS,
  79HIS_S115T, 79ALA, 79ALA_S115T}, on the pre-Upside-3 engine (`$P/obj`, md5 `7136d99a...`). Logs
  `popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<jobid>.out`, run dirs `popepopg_REMD_mdw2/<V>.<ff>/`, seeds
  `seeds/<V>.<ff>.up` (built and verified in `checks/placefix_20261017/`); the old ff3.1's were retired 10-26
  (`*.bak_pre_merge_20261026`). Every block closed clean (no ROLLBACK, DESTROYED or NaN; block 5's final KE/1.5kT
  0.986-1.009); block 1's TM4 in findings 1.35. They run force fields that carry the retraining loss (findings 1.37).
  At the block ends the T 0.70 slot's protein potential still falls (ff3.0 79HIS -1357, -1453, -1531, -1681) at steady
  Rg 19-21 A, one slot's instantaneous value, read beside findings 5.6's block 2 / block 3 disagreement. Groups in
  `run.<r>.up` (each block opens with its 287-frame calibration group ending at t 18.0; `checks/glpg_groups_20261022.py`),
  block 2 / block 3: 79HIS 19-37 / 38-56, 79HIS_S115T 18-36 / 37-54, 79ALA 19-37 / 38-55, 79ALA_S115T 19-37 / 38-56.
  **Read 14:56 11-04** (`checks/tm4_ff30_blocks_20261104/tm4_ff30_blocks_20261104.txt`, job 49814556; block 1
  reproduces the 10-18 readout, the block ranges these): TM4 at T 0.70 over blocks 1-5 (phi/psi 134-151; DSSP alpha 135-151): 79HIS holds (0.98-1.00; 0.92-0.96), 79HIS_S115T
  holds (0.97-0.99; 0.91-0.94; GLY149 at phi > 0 rising to 0.16 in block 5), 79ALA frays from block 4 (0.91; DSSP 0.74,
  0.71; weakest 141 and 149, GLY136/143/149 at 0.10-0.19), and 79ALA_S115T's GLY149 sits at phi > 0 in 0.26-0.48 of
  frames from block 2 (DSSP 0.65-0.75). So under ff3.0 TM4 frays in the H145A (79ALA) backgrounds and holds in the
  79HIS ones. The replicas were then
  deleted ("Disk").
- **HDX dG on the four ff3.0 chains, blocks 2 and 3** (user 10-22: "generate delta G plot into ~/Downloads, and I will
  see if we need longer run"; on midway3 caslake, since midway2's venvs lack mdtraj, pymbar and matplotlib). Readout
  findings 5.6: the blocks disagree well beyond their error bars, so one block does not converge the dG. Scripts,
  prechecks, ladder travel and job ids in `C2=/project/trsosnic/yinhan/checks/hdx_ff30_placefix_20261022/`
  (`submitted.txt`; `hdx_block.sbatch` is the recipe's `hdx_block2.sbatch` with BLOCK a parameter, `extract_block.py`
  md5 `cede3aed` writes each log line in one call). Work dirs `<V>.ff_3.0_gdepth/hdx_sc600x12_ff_3.0_gdepth_block<k>/`
  for blocks 2, 3 and pooled 2-3 (frames per rung 11,021-11,554; resolution limit at T 0.85 6.25-6.29). The four
  figures `glpG-RKRK-<V>_POPEPOPG_dG_vs_residue.png` (user: "only 4 figures, one for each protein variant"; PS, pooled
  blocks 2-3, through the pipeline's `plot_ref_style.py` unchanged, md5 `fad2166e`, `--temperatures 0.75,0.80,0.85`),
  with `ref_style_stats.txt`, `summary.txt` (block 2 against block 3) and `remd_mixing.txt`, are in `$C2/figures/` and
  the MacBook Pro's `~/Downloads/glpG_HDX_ff3.0_blocks_20261022/`; each figure's source is its work dir's
  `results/glpG-RKRK-<V>_dG_profiles.npz`. **The HDX work dirs are kept** (user 10-22: "do not cleanup for now"); their
  frame copies (`block2/`, `block3/`, `block2-3/`, `joined/`, 388 G) were deleted 11-03 at the user's word (disk
  section), so the 12 hold 29 G; `hdx_block.sbatch` (`extract_block.py`) rebuilds them from the chains into a new
  work dir. Their membrane term carries the A61 accessibility shift (plan.md Known Errors; findings 5.6): 61383990
  `hdx_accfix` re-ran only the dG step with the shift undone (`$C2/accfix_ps.py` md5 `b5186e2d` into `<work
  dir>/acc_shift_corrected/results/`, then the pipeline's `4.calc_D_uptake.py`; `$C2/accfix_readout.py` `27be5c1c`
  writes `$C2/tm4_accfix_blocks2-3.txt`). Peptide 140-144's amides 142-144 stay censored at 298 K in 79ALA and
  79ALA_S115T, agreeing with HXMS on that peptide in 79ALA; Val145 comes down to 4.02 (79ALA). The corrected npz's
  seven amides before a proline or the chain end (labels 13 27 58 61 127 151 209) keep the shifted value, so no figure
  is drawn from it.
- **A61-fixed accessibility and censored display for those four HDX sets** (user 11-01: the figures' TM bands carry
  lone finite points the implicit figures lack; the depth-aware rules were scored and set aside, the user chose
  "A61 fix + censor"; findings 5.8). `C3=/project/trsosnic/yinhan/checks/hdx_acc_flicker_20261101/`; the 10-22 work
  dirs are only read.
  - 49726840 `hdx_accgeom` (midway2, done 09:51): per-frame amide geometry, `$C3/<variant>/geom_<r>.npz` (~2 G),
    scored by `$C3/score_rules.py` (`score_79HIS.txt`, figures `$C3/fig/`), MBAR checks `mbar_check_79HIS.txt`,
    `overlap_{implicit,hybrid}.txt`.
  - 49727925 `hdx_a61acc` (midway2, done 10:43): views by the repo's `martini_hdx_project.py`, accessibility by its
    `martini_hdx_membrane_accessibility.py` (Upside 3 tree, md5-equal to the repo), PS by `combine_hdx_protection.py`;
    log `$C3/a61_acc_49727925.out`, step logs `$C3/a61/logs/`, output `$C3/a61/<variant>/{views,acc,block2-3}`. Its
    block split failed on a too-strict Energy check (each work dir carries its own mean reference); `split_blocks.py`
    was fixed and rerun on the login node, writing `$C3/a61/<variant>/{block2,block3}` (336 of 336 PS).
  - The dG step ran on the Mac Studio (caslake would have started it 14:50): `$C3/diag/a61_local.sh`, log
    `$C3/diag/a61_local.log`, profiles, figures and `summary.txt` in `$C3/a61/mac_dG/`, copied to the 11-02 deck's
    `figs/hdx_ff30_a61_20261101/`. Slides 17-18 rebuilt from it.
  - 49737262 `hdx_coreacc` (midway2, done 15:05; user 11-01 afternoon: TM4 and TM5 still show finite points):
    accessibility with the core test (plan.md Known Errors) by `$C3/code/martini_hdx_membrane_accessibility.py` (md5
    `28b7f111`; the repo's copy differs only by one docstring line, and the Upside 3 tree still holds the 64e9c698
    version, sync pending) on the a61 views, radii 7.00 and 8.12 A in all four; then PS and blocks as above (336 of
    336). `$C3/core_acc.sbatch`, log `$C3/core_acc_49737262.out`, step logs `$C3/core/logs/`, output
    `$C3/core/<variant>/{acc,block2-3,block2,block3}`. `$C3/code/check_core.py` checked the new test against the
    geometry arrays (agreement 1.00000). The dG step ran on the Mac Studio; its script, log, profiles, figures and the
    band-hole and event checks are archived in `$C3/core/mac_dG/`, copied to the 11-02 deck's
    `figs/hdx_ff30_core_20261101/`; slides 17-18 rebuilt. Nothing of this runs.
  - `$C3` holds 18 G (views 9.4 G, geometry ~2 G); the views are needed only to recompute the accessibility. Deleting
    them waits on the user, as the other HDX work dirs do.

**Line 3: Upside 3 validation of the debug (midway2; user 10-22: "submit Peng 16 validation protein with ff2.1 and glpG
hybrid simulation with ff2.1, on midway2, as a validation for our entire debug"; plan.md Phase 24 2b).**
- The tree: `/beagle3/trsosnic/yinhan/upside3`, the Mac working tree at A231 (rsync of `src py parameters cmake tests
  example/16.MARTINI`, no `.git`), built on the midway2 login node into `obj/` (GCC 10.1; upside md5 `670c9c87`,
  libupside.so `0a4b83d3`; both still carry `set_fast_math`, 0b), environment `env_upside3.sh`, its own `.venv` (the
  shared venv's versions plus mdtraj 1.9.9 from source; glpG's DSSP string md5-equal to the Mac's mdtraj 1.11). Its
  ff_2.1 and common rama files are md5-equal to the shared tree's.
- **Peng:** 32 arms `u3_<prot>_<kind>` in 10.5 h chunks, the benchmark's own `bench_run.py` and `bench.sbatch` with
  only the base directory and env file changed, in `upside3_validation/peng/` (logs `logs/<prot>_<kind>_ff_2.1_<id>.out`,
  runs `runs/<prot>_<kind>_ff_2.1/`, ids `submitted.log`; `submit_peng.sh` is idempotent). The control is released
  ff2.1 on the pre-audit binary, line 1's arms. COMPLETE, 30 of 32 at 12:40 10-29: all 16 native, and BBA, WWdomain,
  homeodomain, top7, lambda, alpha3d, proteinB, NTL9, BBL, gpW, NuG2, proteinL, hyp (00:06 10-29) and cspA (08:24;
  final frames finite, KE/1.5kT 0.995-1.011 and C-N mean 1.32-1.34 A over the last 50 frames of 14 replicas) de novo.
  proteinG de novo COMPLETE 10-30 and ubiquitin de novo 12:12 11-01 (final frames finite, KE/1.5kT 1.003-1.007 over
  14 replicas, cold rung Rg 11.1 A), so all 32 are COMPLETE. No collector scores them.
- **Peng scored against Upside 2** (findings 3.15): `peng/score/score_compare.py` (md5 `e17d977f`; reproduces released
  ff2.1's scores exactly), `score.sbatch`, `submit_scores.sh` (idempotent; jobs `sc_cmp_<upside3|upside2>_<arm>`, ids
  in `peng/score/submitted.log`) into `peng/score/ff_2.1_upside3/` and `ff_2.1_upside2/`, read by `compare_u3_u2.py`
  (`compare_u3_u2.txt`). Native, 16 proteins: Upside 3 0.564 against Upside 2 0.552 (+0.012 [-0.001, +0.028], p 0.56).
  De novo, all 16 pairs (rerun 00:12 11-02; the 10-26 table kept as `compare_u3_u2.txt.bak_20261026`): 0.396 against
  0.427, -0.030 [-0.049, -0.014], p 0.0042, lower in 14 of 16, one run per arm per engine (on 9 pairs 10-27 it was
  -0.025, median -0.005, carried by top7 -0.108 and alpha3d -0.071); reweighting each engine's frames by the other's
  Hamiltonian predicted -0.002 on those 9 (scripts and outputs `/project/trsosnic/yinhan/checks/u3_vs_u2_20261027/`,
  whose `pull_cold.py` + `score_local.py` reproduce `score_compare.py`).
- **Peng seed replicates**, the test of that gap: top7 and alpha3d de novo with ff_2.1, seeds 101 and 202, on Upside 3
  and on Upside 2 (the original arms ran seed 1), 8 arms `u<3|2>s<seed>_<prot>_denovo` **49634596-49634603**, submitted
  and R from 20:32 10-27 from the Mac Studio (second chunks R since 04:13-06:03 10-29: Upside 3 alpha3d 2.12 M / 1.88 M
  of 3.04 M and top7 0.90 / 1.19 of 2.15 for seeds 101 / 202, Upside 2 alpha3d 1.42 / 1.66 and top7 1.00 / 0.97) (`submit_seeds.sh`, idempotent; targets 2.15 M and 3.04 M time units). In
  `upside3_validation/peng_seeds/` (logs `logs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>_<id>.out`, runs
  `runs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>/`, ids `submitted.log`). `bench_run.py` (`c79004af`) is
  `peng/bench_run.py` with the seed and the engine name (checked against `UPSIDE_HOME`) as arguments; `bench.sbatch`
  (`5ab60bd9`) sources `env_upside3.sh` or `env_shared.sh` by `ENGINE`; Upside 2 seeds each system seed + index, so 101
  and 202 share no stream with seed 1. `ff3_benchmark`'s scripts are untouched (ff3.1's arms read them every chunk).
  All 8 COMPLETE (the last 09:17 10-31); scored by `score/submit_seed_scores.sh` (`f23c8d41`), jobs `sc_seed_*`
  49743788-49743795 submitted 12:18 11-02 (user: "Score 8 seed replicates"; ids in `score/submitted.log`, logs
  `score/logs/`), all COMPLETED by 12:41. Read 16:50 11-02 by `score/seed_spread.py` (the spread over seeds 1, 101,
  202 within each engine against the Upside 3 - Upside 2 difference): top7 de novo -0.087 +- 0.031 SE, beyond its
  seed spread (sd 0.023 and 0.049); alpha3d de novo -0.086 +- 0.064, within it (sd 0.064 and 0.091).
- **glpG-RKRK-79HIS hybrid with ff2.1** (`upside3_validation/glpg/`): `run.py` is the tree's `example/16.MARTINI/run.py`
  (A232) with only its three input paths set to the cluster copies; 48 replicas T 0.70-0.90, 4 blocks each queueing
  the next (wall 35:45). Prep log `outputs/martini_glpG-RKRK-79HIS_hybrid/glpG-RKRK-79HIS.prep.log`; run dir and block
  logs `outputs/martini_glpG-RKRK-79HIS_hybrid/REMD/` (`glpG-RKRK-79HIS.remd_block<k>.log`). Re-running `python3 run.py`
  there submits the next block or reports the count; raise `n_blocks` to extend. Loading the cluster libupside.so sets
  flush-to-zero in the Python process (numpy's "smallest subnormal is zero" warning; 0b's `set_fast_math`).
  - Prep COMPLETED 02:54 10-23: KE/1.5kT 0.997 and 0.976 in the stage 6.0 / 6.1 minimisations, 0.586-0.771 over the
    500-step MD stages 6.2-6.5 (each 4.5 time units against the thermostat's 5.0), 1.002 in 6.6, 1.018 in the 7.0
    burn-in, 1.005 in 7.0; core checks' largest pair energy 1.97-5.52 kT (limit ~32 kT); stage-7.0 solvation gate
    passed (0.0% of 39 belt sites without an environment bead, 0.4% without an acyl-tail bead within 8 A; limit 10%);
    box 100.1 x 100.1 x 123.7 A, 260 lipids.
  - Block 1 ran to its time limit (step 239,976, t 2,159.8; KE/1.5kT 0.997-1.004), then `run.py` failed to write
    `REMD/block_count` on the `/beagle3` quota (14:11 10-24, "Disk"); `block_count` was set to 1 by hand and block 2
    submitted 23:53 10-24 (user: "Resubmit both"). **Running the `remd` stage without `block_count` starts over as block
    1 and copies the seed over every replica file.** Block 1 read 10-26 (findings 3.15; `glpg/glpg_u3_block1.py`,
    `glpg_u3_prep_tm4.py`, outputs `glpg_u3_vs_u2_block1.txt`, `glpg_u3_prep_tm4.txt`): KE/1.5kT 1.00, C-N and five TM
    helices as Upside 2's hybrid; TM4 lost (DSSP alpha 0.44 at T 0.70) from preparation stage 7.0, ff2.1's own (the
    ff3.0 test, below). HDX waits on the user's choice of force field.
  - Block 2 COMPLETED 12:59 10-26 at the time limit (step 235,532, t 2,120; KE/1.5kT 0.996-1.003; no NaN, ROLLBACK or
    error line). Block 3 COMPLETED 21:41 10-28; block 4 of 4 **49638780** R from 21:42 10-28 (step 69,993 at 08:10
    10-29, no error line).
  - **Blocks 2 and 3 read 09:30 10-29** (`glpg/glpg_u3_block1.py` on `output_previous_1` and `_2` of run.0 and
    run.23; output `glpg/glpg_u3_block23.txt`; findings 3.15). At T 0.70, block 2 / block 3: KE/1.5kT protein 0.999 /
    1.003, lipid 1.000 / 1.001; C-N 1.310 +- 0.119 A in both (block 1 1.312 +- 0.120); Rg 19.1 / 19.3 A; TM1
    dihedral-box 0.96 / 0.98; **TM4 135-151 DSSP alpha 0.26 / 0.14** (block 1 0.44; block 3's last quarter 0.02), at
    T 0.795 0.34 / 0.21. TM4 keeps unwinding under ff2.1; the other helices hold (DSSP 0.70-0.96).
  - **Block 4 49638780 COMPLETED 09:00 10-30, read 23:50 10-31** (`/output` of run.0 and run.23;
    `glpg/glpg_u3_block4.txt`): KE/1.5kT protein 0.998, lipid 1.002; C-N 1.311 +- 0.119 A; TM1 dihedral-box 0.98;
    **TM4 DSSP alpha 0.05** (every residue at most 0.08), at T 0.795 0.17. TM4 is gone under ff2.1; the run is done.
- **glpG with ff3.1, the TM4 test** (user 10-26: "diagnose upside3 ... first double check if it is because you run it
  with ff2.1 instead of ff3.0 (now you should run it with ff3.1)"; findings 3.15). A234 in the Upside 3 tree 10:10
  10-26 (only these files): `py/martini_prepare_system_lib.py` `200ed3ac`, `example/16.MARTINI/run.py` `48109b7e`,
  `run_sim_hybrid.sh` `dfa30180`, `readme.md` `fad05dca`, and `py/martini_prepare_system.py` (since A235's);
  `parameters/ff_3.1` from line 2. `upside3_validation/glpg_ff31/run.py` (`3faa5db4`) is `glpg/run.py` with
  `force_field = "ff_3.1"` and `n_blocks = 1`, the same prep seed 2026 and REMD seed, so stage 7.0 differs from the
  ff2.1 run only in the force field. Outputs `glpg_ff31/outputs/martini_glpG-RKRK-79HIS_hybrid/`. Prep 49573353
  COMPLETED 13:37 10-27: **TM4 lost** (`glpg/glpg_u3_prep_tm4.py`, md5 `457ddb59`: stage 7.0 TM4 DSSP alpha last frame
  0.59, last tenth 0.36; ff2.1's last tenth 0.41; the courses in findings 3.15), so the map does not rescue TM4 and
  the run stays at 1 block. REMD block 1
  **49573354** COMPLETED 05:20 10-29 at the time limit (log `REMD/glpG-RKRK-79HIS.remd_block1.log`; final KE/1.5kT
  0.998-1.002, no NaN, ROLLBACK or error line; `block_count` 1). With `n_blocks` 1 nothing follows. **Read 09:30
  10-29** (`glpg_u3_block1.py` on `/output` of run.0 and run.23; `glpg_ff31/glpg_ff31_block1.txt`): at T 0.70
  KE/1.5kT protein 1.004, lipid 1.006; C-N 1.312 +- 0.120 A (min 0.68, the lowest single value of the three runs;
  max 1.80); TM1 dihedral-box 0.996; **TM4 DSSP alpha 0.42** (last quarter 0.32), 135-145 near 0.7 and 146-151 at
  0.00-0.02, so TM4's C-terminal half is gone; at T 0.795 0.39. ff2.1's block 1 read 0.44 / 0.39.
- **glpG with ff3.0, the matched TM4 test** (user 10-26: "ff3.0 surely fix the TM4 issue even before the bug fix in
  upside2 hybrid, so we want to run it to completely be sure", then "do submit ff3.0"; findings 3.15). A235: the
  preparation compares only the three `sidechain.h5` datasets the side-chain/lipid table is built from, so ff_2.1,
  ff_3.0 and ff_3.1 pass; in the Upside 3 tree 14:10 10-26, `py/martini_prepare_system.py` `802e1a90` (old
  `.bak_pre_A235_20261026`, `5e6b4627`). The tree's `parameters/ff_3.0` is the released `ff_3.0_gdepth` (all six files
  md5-equal), which the placement-fixed Upside 2 chains ran. `upside3_validation/glpg_ff30/run.py` (`db1cf3b4`) is
  `glpg_ff31/run.py` with `force_field = "ff_3.0"`, same seeds; outputs and prep log in
  `glpg_ff30/outputs/martini_glpG-RKRK-79HIS_hybrid/`. Prep 49589828 COMPLETED 14:38 10-27: **TM4 held** (last frame
  0.88, last tenth 0.82; TM1 0.95): the ff2.1 run's loss is ff2.1's own with correct placement, not Upside 3's. REMD
  block 1 **49589829** COMPLETED 05:21 10-29 at the time limit (log `REMD/glpG-RKRK-79HIS.remd_block1.log`; final
  KE/1.5kT 0.997-1.004, no NaN, ROLLBACK or error line; `block_count` 1). **Read 09:30 10-29** (`glpg_u3_block1.py`
  on `/output` of run.0 and run.23 while block 2 is queued; `glpg_ff30/glpg_ff30_block1.txt`): at T 0.70 KE/1.5kT
  protein 1.006, lipid 1.004; C-N 1.312 +- 0.119 A; TM1 dihedral-box 0.985; **TM4 135-151 DSSP alpha 0.85** (halves
  0.84 / 0.85, last quarter 0.84) against Upside 2's placement-fixed ff3.0 chains' 0.96; the shortfall is 135 and
  136, never helical (137 0.80), while the helix DSSP finds at 137-151 reads 0.96; at T 0.795 0.81 (Upside 2 0.93).
  ff3.0 holds TM4 through block 1 on Upside 3, two residues shorter at its N-terminal end.
  - **Block 2 49634605 COMPLETED 13:21 10-31** (R from 02:03 10-30; log `REMD/glpG-RKRK-79HIS.remd_block2.log`, no
    error line). **Read 23:50 10-31** (`output_previous_1` of run.0 and run.23; `glpg_ff30/glpg_ff30_block2.txt`): at
    T 0.70 KE/1.5kT protein 0.997, lipid 1.002; C-N 1.311 +- 0.119 A; TM1 dihedral-box 0.995; **TM4 DSSP alpha 0.82**
    (halves 0.80 / 0.84), 137-151 at 0.93 and 135-136 never helical, as in block 1; at T 0.795 0.85. ff3.0 holds TM4
    through block 2. **Block 3 49686445 ended 23:38 11-01 at its time limit** (t 2,090; KE/1.5kT 0.996-1.004, no
    error line; `block_count` 3 of 3, so the run is done; Slurm reads TIMEOUT, stamped 09:43 11-02). **Read 00:09 11-02**
    (`/output` of run.0 and run.23 after the block ended; `glpg_ff30/glpg_ff30_block3.txt`): at T 0.70 KE/1.5kT
    protein 1.001, lipid 1.001; C-N 1.311 +- 0.119 A; TM1 dihedral-box 0.988; **TM4 DSSP alpha 0.85** (halves 0.84 /
    0.85), 138-151 at 0.98, 135-136 never helical; at T 0.795 0.82. ff3.0 holds TM4 through all three blocks.
    `n_blocks` was raised 1 -> 3 by hand
    from the Mac Studio 20:32 10-27 (old `run.py.bak_pre_nblocks3_20261027`). Block 1 loaded `n_blocks = 1` and does
    not queue block 2, so 49634605 was submitted with `--dependency=afterok:49589829` and the sbatch line of `run.py`'s
    `submit()`; it reads `block_count` 1 and queues block 3, so blocks 2-3 pool for HDX as Upside 2's ff3.0 chains did
    (findings 5.6). **Do not run `python3 run.py` there before block 1 ends:** `block_count` does not exist yet, so it
    would submit a second block 1 that copies the seed over every replica file. Running `run.py` now would queue
    another block.
- Next action: at each REMD block end read KE per class, C-N, TM helicity (DSSP) against findings' ff2.1 glpG numbers,
  and TM4's protein protection at 298 K by HDX through the fixed analysis code (plan.md Phase 24 2b 4); when the Peng
  arms finish, score them as line 1's collector does and compare with the control and Peng's FF2.

**Line 4: Boltz-2 HDX sampler (midway3).** Every earlier set is finished and read out: Phase 1's sets, the first
BioEmu2 job and the controlled 10^5 tests (user 10-25: "propose and apply fix"), Boltz-2 Ub psWT with the query-only
MSA and BioEmu2 EHEE and Ub psWT ("Finished runs"; results in `diffusion_dynamics.md` Phase 1 and section 7). Now
Phase 2 (user 10-27: no retraining of Boltz, small trained add-ons allowed), a learned tilt sampled through Boltz-2's
own steering. `$B` = `/beagle3/trsosnic/yinhan/boltz_hdx`.
- **E1 pilot, `boltz_addon_e1_pilot_ubqps_v100` 61663101** (midway3 gpu, V100; COMPLETED 03:11-03:15 10-28 on
  midway3-0279, read out 10:50 10-28): `$B/addon/e1.sbatch` running `e1_rare.py` (with `sampler.py`, `opening.py`)
  once per line of `$B/addon/e1_pilot_ubqps.txt` on Ub psWT's processed input
  (`runs/boltz_phase1_ubqps_1e5_v100_61428012/ubqps_md_s10/s42/boltz_results_ubqps/processed`). The sampler parity
  check on the GPU says "equal" (steering off against Boltz's own sampler, 8 samples, seed 7; 132 s), as it did in the
  MacBook Pro's CPU sandbox (`scratchpad/boltz_hdx/addon/sandbox/`, on that computer only: a Python 3.12 uv venv, a
  copy of `$B/boltz`, `$B/cache`'s weights and CCD molecules, and the processed input; another computer rebuilds it
  from those). The 64 unguided samples (106 s) agree with the pooled 10^5 flags: 852 openings over the 41 amides
  the pool opens at least 10 times, against 854.5 expected, and no amide's count is off by a two-sided binomial test
  at p < 0.05; longest C-N at most 1.56 A. Log `$B/logs/boltz_addon_e1_pilot_ubqps_v100_61663101.out`; output
  `$B/addon/runs/boltz_addon_e1_pilot_ubqps_v100_61663101/unguided64.npz` (md5 `753ac79f`), local copy and
  `readout.txt` in `scratchpad/boltz_hdx/addon/runs/<same>/`; ledger `$B/addon/submitted.log`. The pilot's code is
  kept in `$B/addon/pilot_61663101/`.
- **E1 tuning, `boltz_addon_e1_tune_ubqps_v100` 61701305** (midway3 gpu, V100, 2:30; submitted 11:11 10-28, user
  10-28 "run E1"): the same `e1.sbatch` on the same input with `$B/addon/e1_tune_ubqps.txt`, 12 steered runs of 256
  particles (V70 and residue 61, kappa 0.5-4, guidance 0.01 and 0.05), through the estimator-faithful potential
  (`protection.py`, `opening.py`; `diffusion_dynamics.md` Phase 2 E1). Code and `inputs/{ubqps,heeh}_protection.npz`
  md5-equal with the MacBook Pro's `addon/` at 11:10 10-28; fixes go into `$B/addon/` only between jobs. Log
  `$B/logs/boltz_addon_e1_tune_ubqps_v100_61701305.out`; outputs `$B/addon/runs/boltz_addon_e1_tune_ubqps_v100_61701305/`.
  Finished 18:06 10-28: the log ends `[done]`, no `[e1] ... FAILED` line. Read 10-28 on the Mac Studio
  (`diffusion_dynamics.md` Phase 2 E1): the rule picks kappa 2, guidance 0.05, but the Feynman-Kac estimate fails
  above kappa 1 (I61's `log_z` 22-76 nats under its provable floor), so **E1 proper is held**.
- **E1 diagnosis, `boltz_addon_e1_diag_ubqps_v100` 61755382** (midway3 gpu, V100, 3:00; submitted 23:41 10-28, user
  10-28 "go ahead on boltz job"; COMPLETED 23:43-01:04 10-29 on midway3-0277, six runs, no FAILED line): `e1.sbatch` with `$B/addon/e1_diag_ubqps.txt`, six
  runs on I61 at guidance 0.05 (kappa 1 and 2; 256 and 1,024 particles; seeds 1 and 2), the sampler now recording a
  per-resampling `trace` (`sampler.py` `83cdb782`, `e1_rare.py` `7f0e1ee5`; old `.bak_pre_trace_20261028`). Log
  `$B/logs/boltz_addon_e1_diag_ubqps_v100_61755382.out`; outputs `$B/addon/runs/boltz_addon_e1_diag_ubqps_v100_61755382/`.
  Reading fixed before the data: `diffusion_dynamics.md` Phase 2 E1, "The diagnosis job" (first, the two repeated
  runs must reproduce the tuning's `log_z` -7.93 and -30.49). **Read 08:20 10-29: rule 1 fails** (-6.33 and -28.35,
  samples different wholesale at the same seed and settings), so nothing else is read; whether Boltz-2 reproduces a
  seed across jobs on these GPUs is not known (the nodes differed: tuning 0280, diagnosis 0277). npz, log and
  `readout.txt` in the Mac Studio's `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_diag_ubqps_v100_61755382/`.
  E1 proper stays held.
- **E1 reproducibility, `boltz_addon_e1_repro_ubqps_v100` 61786565** (midway3 gpu, V100, 1:30; submitted 09:24 10-29,
  user 10-29 "fix them all"; COMPLETED 09:29-10:08 10-29 on midway3-0277, the diagnosis's node): `$B/addon/e1_repro.sbatch`
  runs the tuning's I61 kappa 1 guidance 0.05 run (256 particles, seed 1) six times, each in its own process:
  `new_a/b` with the current code, `old_a/b` with the tuning's (`$B/addon/repro_old/`: the `.bak_pre_trace_20261028`
  `e1_rare.py` `20bdd6b7` and `sampler.py` `9ae0eb7e`, with `opening.py` and `protection.py` md5-equal to `addon/`'s),
  `det_a/b` with the current code under PyTorch's deterministic mode (warn only). Log
  `$B/logs/boltz_addon_e1_repro_ubqps_v100_61786565.out`; outputs `$B/addon/runs/boltz_addon_e1_repro_ubqps_v100_61786565/`.
  Read with `python repro_compare.py <that dir> runs/boltz_addon_e1_tune_ubqps_v100_61701305/r61_k1_g0p05.npz
  runs/boltz_addon_e1_diag_ubqps_v100_61755382/r61_k1_g0p05_n256_s1.npz` (any numpy; `addon/`), by
  `diffusion_dynamics.md` Phase 2 E1's rules ("The reproducibility job"), and the `det` runs' deterministic-mode
  warnings in the log. **Read 10:15 10-29: rule 1.** Every pair differs wholesale (current and tuning code alike);
  the runs part at step 0, before any resampling; deterministic mode (no warning) shrinks the step-0 difference
  about a thousandfold and leaves it, so the source lies outside what it governs (not identified). Eight `log_z` at
  the setting: mean -8.26, s.d. 1.58. npz, log and `readout.txt` in the Mac Studio's
  `scratchpad/boltz_hdx/addon/runs/boltz_addon_e1_repro_ubqps_v100_61786565/` (readout also on `/beagle3`).
- **The diagnosis 61755382's rules 2-3, read statistically 10:40 10-29** (user 10-29: "go with your
  recommendation"; rule and readout in `diffusion_dynamics.md` Phase 2 E1; `addon/diag_rules23.py`, output
  `runs/boltz_addon_e1_diag_ubqps_v100_61755382/rules23.txt`, both also on `/beagle3`): the kappa 2 shortfall does
  not shrink from 256 to 1,024 particles (means -31.0 and -44.1), which by the rule points to a defect in the
  weights, at low power (the two 1,024 seeds 26.5 nats apart). Step 0 carries 26 nats in every run; the loss sits in
  the potential increments, the guidance log-ratio summing to about -0.1. E1 stays held; the next step is the user's.
- **The toy test, 12:40-12:55 10-29, on the Mac Studio's CPU, no job** (user 10-29: "yes"; rule and readout in
  `diffusion_dynamics.md` Phase 2 E1): `addon/toy.py` runs `sampler.py` unchanged on Boltz's own `AtomDiffusion`
  (the midway3 tree's Python files in `addon/sandbox/`) with a Gaussian-ubiquitin denoiser and a known normalizing
  constant from 10^6 plain samples; `toy_readout.py`, output `runs/toy/readout.txt`. All four settings (kappa 1 and
  2, guidance 0 and 0.05) pass: mean Z-hat / Z within 2 standard errors of 1 at 1,024 and 4,096 particles, and the
  shortfall shrinks with the particle count. `sampler.py`'s weights and `log_z` are right, which points I61's
  shortfall to the estimator's variance with Boltz-2's network (in I61 the selected particles' energies rise again
  between resamplings, in the toy they do not). E1 stays held; the next step is the user's.
- **Resampling at low noise only** (user 10-29: "run both"; rule and readouts in `diffusion_dynamics.md` Phase 2 E1):
  `sampler.py`'s `resampling_t_max`, set by `e1_rare.py --resample-below 16` (sigma_data), resamples only where
  t_hat <= 16 A; inert by default (toy runs bitwise equal). The toy passes with it (`runs/toy_low/readout.txt`).
  **61805263** (`boltz_addon_e1_lownoise_ubqps_v100`, V100, 1 h; submitted 13:19 10-29 from `$B/addon`, task file
  `e1_lownoise_ubqps.txt`, logged in `addon/submitted.log`; log
  `$B/logs/boltz_addon_e1_lownoise_ubqps_v100_61805263.out`, outputs
  `addon/runs/boltz_addon_e1_lownoise_ubqps_v100_61805263/`): I61 at kappa 2, guidance 0.05, 256 particles, seeds 1
  and 2. Read: the collapse is removed if both seeds' ln Z-hat >= -11.1 (old schedule -28.4, -33.7). **Ended
  COMPLETED 13:32 10-29** (log `[done]`, no FAILED line); readout 14:40 10-29 (`lownoise_read.py`, run from midway2):
  ln Z-hat -15.05 and -20.26, **not removed**; the first resampling at t_hat 13.8 charges -25.3 and -22.7 nats and
  is not repaid; calibration criterion not met (Feynman-Kac open fraction 2.7e-7 and 1.6e-9 against 3.8e-4). E1
  stays held; the next step is the user's (diffusion_dynamics.md Phase 2 E1).
  Code `sampler.py` `21df0b41`, `e1_rare.py` `e4314698` (previous as `.bak_pre_lownoise_20261029`).
- **E1 particle count at low noise** (user 11-02: "run the GPU jobs. We need to continue working on trying to make
  Boltz work for HDX", then "Particle scaling, full"; rule fixed before the data in `diffusion_dynamics.md` Phase 2
  E1, "Particle count at low noise"): the 61805263 setting (I61, kappa 2, guidance 0.05, `--resample-below 16`) at
  more particles, code unchanged. Submitted 13:09 11-02 from `$B/addon`, `e1.sbatch` with task files
  `addon/scale/<tag>.txt`, V100 (`--constraint=v100`), logged in `addon/submitted.log`: **62073874** (`n1024`, seeds
  1-4 in turn, 3 h limit), **62073875-62073878** (`n4096_s1`-`s4`, 3 h), **62073879-62073880** (`n16384_s1`-`s2`,
  10 h; about 5 h each by the 256 and 1,024 timings). Logs `$B/logs/boltz_addon_e1_scale_<tag>_v100_<id>.out`,
  outputs `addon/runs/boltz_addon_e1_scale_<tag>_v100_<id>/`. Readout **62073881** (`addon/scale/read.sbatch`,
  caslake, afterany on the seven) runs `scale_read.py` (`5f959e3d`) on them and 61805263's two 256 runs into
  `addon/runs/scale_readout.txt`, with any `[e1] ... FAILED` line appended. All PD at 13:09 (Slurm's estimate for
  the first 16,384 job 19:32). Read the file when it appears; the rule's verdict is its last line. Outcome:
  62073874-78 COMPLETED 14:54-15:33; 62073879-80 ran their samplers to the end (17,585 s and 18,259 s) and both failed
  at `np.savez` of the result, 19:48 and 20:00 (`OSError: [Errno 122] Disk quota exceeded`: the `/beagle3` inode soft
  quota, grace expired), so `[e1] ... FAILED` and no result file; the readout 62073881 FAILED in 1 s at 20:01 with no
  log. Read 11-03 on the Mac Studio from the finished runs (npz pulled, `scale_read.py` there; readout copied back to
  `addon/runs/scale_readout.txt`): removal fails at 1,024 and 4,096; the rule's other branch needs the 16,384 arms.
  E1 closed 11-03 (user: "proceed with your recommendation"); the 16,384 arms are not rerun.
- **Phase 3 T1, invariance of Boltz-2's own move** (user 11-03: "proceed", to T1 as proposed; rule fixed before the
  data in `diffusion_dynamics.md` Phase 3 T1): **62263853** `boltz_remc_t1_ubqps_v100`, submitted 11:27 11-04 from
  `$B/addon` with `sbatch --job-name=boltz_remc_t1_ubqps_v100 --partition=gpu --qos=gpu --account=pi-trsosnic
  --gres=gpu:1 --constraint=v100 --time=06:00:00 t1.sbatch` (`--test-only` passed first), logged in
  `addon/submitted.log`. Code `remc.py` `e5e9dfb2`, `remc_boltz.py` `7556f0ae`, `t1.sbatch`
  `8902a64b` (copies in the MacBook Pro's `scratchpad/boltz_hdx/addon/`; the Mac Studio's are the old ones). The first
  submission, 62199777 (18:13 11-03, code `6bf12acf` / `7b21e83b`), FAILED 19:31 at the parity check (`DIFFERENT, max
  |d| 36.9 A`; the sampler against itself equal): `remc.schedule` read the module's `num_sampling_steps`, Boltz2's
  default of 5, while the sampler ran 200, so every move would have started at step 3-4 of 5. Fixed in `remc_boltz.py`
  (backups `.bak_pre_nsteps_20261104`; diffusion_dynamics.md Phase 3 T1). Its rerun 62262317 (10:57-11:05 11-04)
  passed the check (`equal`) and FAILED at its first observation: `atom_idx`, a field of Boltz's 43-byte residue
  records, is a strided view torch cannot index with; copied to a contiguous int64 array (backup
  `.bak_pre_stride_20261104`), and the whole run path tested on the Mac with Boltz stubbed before 62263853. The job extracts the packed `.venv` (`checks/disk_reclaim_20261102/boltz_venv.tar`)
  and `cache/mols.tar` to the node's `/tmp` and deletes them at its end, so `/beagle3` gets only the log and the
  result; it reruns `boltz predict` on `inputs/ubqps.yaml` for the processed input (the 11-03 cleanup deleted the old
  one), then the parity check (`[check]` line; the job stops if not equal), then 256 chains x 100 sweeps with 1,024
  same-trunk plain samples (about 3.3 h; `[remc] sweep N` every 10 sweeps). Log
  `$B/logs/boltz_remc_t1_ubqps_v100_62263853.out`; data `$B/addon/runs/boltz_remc_t1_ubqps_v100_62263853/`
  (`predict.log`, `t1.npz`). Any later Boltz GPU job reuses `t1.sbatch`'s node-local extraction instead of unpacking
  on `/beagle3`.
- The readout tools are git-ignored, `scratchpad/boltz_hdx/flags/` and `scratchpad/bioemu2_hdx/readout/`, copied into
  the bundle's `scratchpad/` (scripts and small inputs; the large samples are on the cluster). The flag code's scratch
  engine is built on the MacBook Pro (`scratchpad/boltz_hdx/flags/upside_tree/`: `src`, `parameters/common`,
  `parameters/ff_2.1` and `obj` at repo 6d295b0, 10-27); another computer rebuilds it from `git archive 6d295b0` and
  checks the earlier flags reproduce.

**Shared:** the engine fix (findings 3.12; plan.md Phase 21): `$P/obj/upside` and `libupside.so` from
`checks/angfix_20261018/obj` (md5 `7136d99a...`, `3d9c6366...`), installed 08:35 10-18 with backups
`*.bak_pre_angfix_20261018`; `$P/src/interaction_graph.h` patched (backup alongside). `/beagle3`'s build is not
changed. Energies, forces, simulations and every derivative but the pair term's angular slots are bitwise those of the
old build.

Finished and cancelled rows are deleted; their lessons are in §8 and "Finished runs".

| JobID / where | what | state | next action |
|---|---|---|---|
| **62263853** (`boltz_remc_t1_ubqps_v100`, midway3 gpu, V100) | Boltz Phase 3 T1: Boltz-2's own renoise-denoise move at lambda 0 on Ub psWT, 256 chains x 100 sweeps, with the parity check and 1,024 same-trunk plain samples (line 4); the second rerun (62199777 parity failed, 62262317 strided index; both fixed) | PD from 11:27 11-04; log `$B/logs/boltz_remc_t1_ubqps_v100_62263853.out`; `[remc] sweep N` every 10 sweeps, about 3.3 h | pull `t1.npz`, read with `addon/t1_read.py` by Phase 3 T1's rule, report to the user |
| `$P/training/ff311_glypair*` (no job) | The ff3.1 training debug, `ff311_glypair` and `ff311_glypair_mb2` (plan.md Phase 32; line 2), on the old native target | both at their target 76 (`_mb2` 11:44 10-30, `ff311_glypair` 03:35 10-31 after a NODE_FAIL in step 65) | both arms' round 4 by Phase 32's rule; `ff311_glypair`'s round 3 |
| 30 restarted `b_ff_3.1_<prot>_<kind>`, now chunks **49760923-49760946**, collector next **49750446** (`ff31_collect`) | ff3.1's Peng arms on the pre-audit Upside 2 engine and their collector (plan.md Phase 37; line 2) | BBA and BBL native COMPLETE and scored (final frames finite, KE/1.5kT 1.002-1.008); the 30 broken chains restarted 12:18 11-02: BBA de novo COMPLETE 13:49, WWdomain de novo 18:05, WWdomain native 18:35 and NTL9 native 22:10 (final frames finite, KE/1.5kT 1.000-1.007, 1.003-1.008, 1.004-1.009 and 1.004-1.007, 14 replicas each), WWdomain native and de novo scored 21:45-21:48, NTL9 native 04:23 11-03 by pass 49748728 (`<TM>` 0.199 at T 0.680, released ff2.1 0.249); the other 26 chunks ended 22:51 and resubmitted themselves (49748835-49748860, logs written), none on the throttled midway2-[0316-0327]; of them gpW native 49748839 COMPLETE 05:55 11-03 and alpha3d native 49748860 07:57 (final frames finite, KE/1.5kT 1.002-1.006 and 1.002-1.006, 14 replicas each); the other 24 ended at their wall time 09:22 11-03 and resubmitted themselves (49760923-49760946), PD on Priority at 09:26 (broadwl 179 of 180 nodes allocated); collector next 49750446 (begin 09:49) | each arm's final frames at `[bench] COMPLETE`; `peng_ff31.txt` when DONE |
| **49801508** (`condiv-train`, `ff32_basin` epoch 3 to step 76), insurance **49801596** | ff3.2: `ff32_basin`, the native-basin restraint, DSE off (plan.md Phase 38 item 11; line 2) | R from 07:37 11-04, step 70 at 15:08, 30 min a step, 24 of 24 proteins; step 76 about 18:10 | at 76 its `after_training_ff32_basin.sbatch` submits panel `w32_03`, its unfolded shift and gaps, then stops |
| **49816524** (`ff32_e3_wait`, midway2 1 core, 16 h) | ff3.2's epoch-3 read, detached (user 11-04: "proceed with your suggestion, but make sure the mistakes can be corrected later"; `checks/ff32_epoch3_20261104/`, README) | R from 15:48 11-04, waiting for step 76's after-training line | submits `read.sbatch` afterok on `w32_03`'s panel, shift and gaps; it writes `summary.txt` (the epoch choice and verdict by `choose.py`'s rule, fixed before the data, then the item 11 readouts) and archives the step-76 checkpoint. Read `summary.txt`, then by its verdict: validation, or the fixes |

**Line 5: Upside 3 on the GPU (plan.md Phase 35; midway2 `gpu2` and midway3 `gpu`, never `beagle3`, user 10-29).**
- The tree: `/beagle3/trsosnic/yinhan/upside3_gpu`, a development tree beside the Upside 3 tree (which its validation,
  line 3, reads), filled from the Mac by `scratchpad/gpu_work/sync.sh` (`src py tests cmake`; `--params` adds
  `parameters example/16.MARTINI`). Built on the midway2 login node with `gcc/10.1.0 cuda/11.5`, one binary for sm_37
  (K80) to sm_86. GPU jobs pass the partition on the command line: midway2 `--partition=gpu2 --qos=gpu2`, midway3
  `--partition=gpu --qos=gpu` (midway3.md), account pi-trsosnic, 1 GPU, under 1 h. Environment
  `env_upside3_gpu.sh` (the Upside 3 tree's venv, this tree's code).
- Job scripts in `checks/g0/`: `v1.sbatch CONFIG...` (the G0 probe, then `tests/device_v1.py`), `v2run.sbatch WORK ARM
  SEED_BASE CONFIG N_REPLICA DURATION OBJ [upside option ...]` (one seed of one V2 arm per array task, writing
  `WORK/ARM.SEED.npz`; a GPU arm adds `--gpu` and its partition with `--gres=gpu:1`), `v2cmp.sbatch WORK ARM_A
  ARM_B...` (`device_dyn.py compare`), `prof.sbatch OBJ N_COPY DURATION CONFIG [option ...]` (nvprof), `pair.sbatch`,
  `bitwise.sbatch`, `mom.sbatch`, `scale.sbatch`, `diag.sbatch SCRIPT ARGS`; configs in `checks/configs/`; logs
  `<gpu>_<id>.out`, `v1_<gpu>_<id>.out`, `v2run_<arm>_<id>_<task>.out`, `diag_<id>.out`. Installed `obj/` is build g8a
  (g7b plus A240; 17:22 11-02); beside it `obj_g7b` (g6c plus A239, the g-JF position update as x + b dt/m (p + dt/2
  f + beta/2); installed 10-30 12:15 to 11-02), the earlier builds `obj_g1a` to `obj_g6c` and `obj_hb6b` (no CUDA:
  the host backend, before A239). Install a new one with
  `install_gpu.sh TAG` only when no job of the tree runs. The hybrid fixtures are `checks/configs/1rkl.stage_7.{0,1}.up`.
  The ff3.1 configs `checks/configs/{chig,1ubq,1dfn}_ff31.up` come from `scratchpad/gpu_work/ff31_configs.py`; since
  00:20 11-01 they are the released ff3.1's (only `rama_pot` and `proposal_pot` differ from the 10-30 ones, which are
  in `checks/configs/bak_pre_nativefix_20261031/`), and the tree's `parameters/ff_3.1/rama.dat` is the release
  (`7693f10a`; the old `5b66c78a` kept as `rama.dat.bak_pre_nativefix_20261031`; the other five files were equal).
  midway3's `gpu` QOS takes at most 12 submitted jobs per user, each array task counted.
  Resumed on the Mac Studio 00:05 11-01 ("resume GPU jobs"; repo HEAD `7f31f72` carries the code, equal to the tree's
  `src py tests cmake`; tools from the bundle `/project/trsosnic/yinhan/checks/handoff_gpu_20261031/`). No job of
  line 5 ran 00:31 11-01 to 13:40 11-02.
- **The soluble remainder of the ladder** (user 11-02: "continue the validation", then "Soluble remainder"; plan.md
  Phase 35 item 11, its readings fixed there first). Harness changes, host-backend tested on the MacBook Pro
  (`scratchpad/gpu_work/dev13_home`): `tests/device_v1.py --frames` (md5 `688aefa3`), `tests/device_dyn.py nve`
  (`53912dc7`), `g0/v2run.sbatch` now takes MODE first (`12ee22a2`), new `g0/v2seq.sbatch` (seeds in turn on one GPU,
  for midway3's 12-job limit), `g0/v1x_prep.py`, `g0/v3_dist.py`, `g0/v3dist.sbatch`, and `g0/v3_analyse.py`
  refactored to share its frames (`24722ed9`); every replaced file kept as `.bak_pre_soluble_20261102`. Submitted
  13:40 11-02, all from `checks/g0/`:
  - V1 at scale 62076233 (midway3 V100; log `v1x_v100_62076233.out`), ran 15:06-15:07 11-02, state FAILED because
    `device_v1.py` exits 1 on a failure: 16 Peng proteins (ff2.1, Upside 3's Peng runs) at BP tol 1e-3 and 1e-6, 10
    frames each from `checks/configs/v1x/` (`v1x_prep.log` lists each frame's CA-RMSD and Rg). 30 of 32 configs pass;
    BBA and NTL9 fail `hbond_energy`'s output at one frame, at both tolerances (read and diagnosed: plan.md Phase 35
    item 11, findings 11g); fixed as A240 (findings 3.14), rerun below.
  - V2 energy drift, `device_dyn.py nve`, 4 copies at T 0.85 for 5000 tu: CPU **49746182** (1ubq full, seeds
    5001-5016, `v2n_full/`) and **49746183** (1ubq ff3.1, 6001-6016, `v2n_ff31/`) on broadwl, logs
    `v2run_cpunve*_<id>_<task>.out` (all 32 COMPLETED, about 400 s per seed); V100 **62076234** (seeds 5201-5208) and **62076235** (6201-6208), logs
    `v2seq_v100nve*_<id>.out` (COMPLETED 16:48 11-02, 755 s per seed). Read 16:55 by `device_dyn.py compare v2n_full
    cpu v100` and the same for `v2n_ff31`: Bonferroni 0.65 and 1, so the engines agree; both drift about 90 E_up per
    1000 tu (plan.md Phase 35 item 11, Known Errors).
  - V3 distributions 49746180 (broadwl; log `v3dist_49746180.out`) COMPLETED 11-02 in 8 minutes: the refactor
    reproduces `v3/analysis.txt`; `v3/distributions.txt` passes (Bonferroni 1, plan.md Phase 35 item 11).
- **A240 on the GPU tree** (user 11-02: "fix the bug you find"): the four files `src/hbond.cpp`,
  `src/hbond_interaction.h`, `src/Float4.h`, `src/device_hbond.cu` rsynced alone from the Mac 16:46 11-02 (md5 equal;
  `sync.sh` was not run, since it would have sent the MacBook Pro's older ff3.1 configs over the cluster's), built as
  `obj_g8a` on the midway2 login node (`build_gpu.sh g8a`, no errors). Before its install `home_g8a` held copies of
  `py` and `tests` with `obj` linked to `obj_g8a`, and `g0/v1home.sbatch HOME [option ...] CONFIG...` (`514face5`) ran
  V1 through it. Submitted 16:49 11-02 from `checks/g0/` (ids also in `g0/submitted_a240.log`), both COMPLETED
  17:16-17:17 on midway3-0280 (V100) and read 17:20 by plan.md Phase 35 item 11's V1 rule: **62096441** (log
  `v1x_g8a_v100_62096441.out`; `--frames ../configs/v1x/frames.npz ../configs/v1x/*.up`, the 32 Peng configs) passes
  every quantity at both tolerances, and **62096444** (log `v1_g8a_v100_62096444.out`; the 14 standard configs of
  61891693: chig, 1ubq, 1dfn, 1ubq_rotamer, chig_rama, chig_hbond, chig_springs, 1dfn_backbone, 1ubq_norot,
  1rkl.stage_7.0, 1rkl.stage_7.1, and the three ff3.1 configs) passes every quantity. **g8a installed** 17:22 11-02
  (Mac Studio; no job of the tree on either cluster; `install_gpu.sh g8a`, `obj` = `obj_g8a` by `diff -r`, upside
  md5 `22d8dfc2`, libupside.so `fdc1721c`); g7b stays as `obj_g7b`.
- **Item 12 steps 1-2, build g9a** (user 11-02: "plan that, while still monitoring", then steps 1-2 approved;
  plan.md Phase 35 item 12). g9a is g8a plus batches of different proteins (every system's pair slab sized by the
  largest system's `max_n_edge`, each held to its own; `rama_map_pot`'s maps per system, equal maps uploaded once):
  six `src` files and `tests/device_v1.py --batch` rsynced alone from the Mac 18:14 11-02 (md5 equal; the tree's
  other `src py tests cmake` files are HEAD's but for the old non-GPU `py/martini_hdx_membrane_accessibility.py`),
  with the hybrid `example/16.MARTINI/run.py` and its `remd_engine` option (the tree's old copy kept as
  `run.py.bak_pre_gpuopt_20261102`), built 18:16 on the midway2 login node; `home_g9a` (copies of `py`, `tests`;
  `obj` linked to `obj_g9a`). On the Mac's host backend first: V1 of a one-protein batch identical to HEAD's,
  mixed batches pass every quantity, `--gpu` trajectories bitwise equal to HEAD's, CPU regress 0 of 4260 arrays.
  glpG's V1 configs `checks/configs/v1g/` (`g0/v1g_prep.py`: replica slot 0's `/input` of `glpg_ff31` block 1 copied
  whole by h5py, at BP tol 1e-3 and 1e-6; 10 frames of block 1 across the ladder). Submitted 18:20 11-02 from
  `checks/g0/` (ids in `g0/submitted_g9a.log`): **62105832** (`bitwise_g9a_1ubq_v100_<id>.out`, 14 copies of
  1ubq with rotamer, 100 tu, obj (g8a) against obj_g9a; COMPLETED 19:26 on a V100: 101 frames of 14 systems bitwise equal,
  76.1 and 76.6 us/system/step), **62105835** (`bitwise_g9a_1rkl_v100_<id>.out`, 4 copies of
  1rkl 7.0, 40 tu), **62105836** and **62105838** (`v1b_g9a_ff21_v100_<id>.out`, `v1b_g9a_ff31_v100_<id>.out`: chig,
  1ubq and 1dfn in one batch), **62105839** (`v1g_g9a_v100_<id>.out`, glpG's frames by item 11's V1 rule), and by
  `g0/glpg_probe.sh HOME GPU` (a sandbox `g0/glpg_probe_<gpu>/` holding glpg_ff31's stage 7.0 checkpoint, whose
  `run.py` with `remd_engine = "gpu"` submits one 15-minute REMD block of 48 replicas itself and logs the GPU's
  memory to `gpu_mem_<id>.csv`): **62105864** (V100) and **62105868** (A100), logs
  `glpg_probe_<gpu>/outputs/martini_glpG-RKRK-79HIS_hybrid/PROBE/glpG-RKRK-79HIS.remd_block1.log`. Read 20:30 11-02,
  all on one V100: 62105832 and 62105835 bitwise equal (1rkl 41 frames of 4 systems); 62105836 and 62105838 pass every
  quantity (192 for each of chig, 1ubq, 1dfn in one process, ff2.1 and ff3.1); 62105839 passes glpG's 1040 quantities
  at both BP tolerances. The V100 probe 62105864 ran its 15-minute block of 48 replicas (903 s, 794.7 us/system/step,
  850 tu an hour; GPU memory 4,382 of 16,384 MiB; KE/1.5kT 0.993-1.017), then FAILED writing `block_count` (Errno 122,
  the `/beagle3` quota). The A100 probe 62105868 started 20:37 and FAILED in 2 s with no
  log and no file written (Slurm could not create its log: the quota). After the reclaim: `install_gpu.sh g9a` 21:54
  (`obj` equal to `obj_g9a`, md5 `e928f0fa`; g8a stays as `obj_g8a`); the A100 sandbox (its inputs only) removed and
  `glpg_probe.sh home_g9a a100` resubmitted 21:56 as **62120973** (`g0/submitted_g9a.log`), COMPLETED 23:43: 48
  replicas, 906 s, 628.9 us/system/step, 1,100 tu an hour (1.27 times the V100), GPU memory 4,495 of 40,960 MiB,
  KE/1.5kT 0.994-1.014, `block_count` written. Both probes' replica files deleted after reading (logs, `block_count`
  and `gpu_mem_<id>.csv` kept).
- **Item 12 steps 3-5, submitted 18:25-18:45 11-03** (user 11-03: "close what is validated, and proceed to next
  validations"; soluble proteins closed in plan.md Phase 35; V100 chosen for step 4, `run.py`'s default). Step 3, glpG
  V2: `checks/configs/glpg.stage_7.0.up` (md5 `20ee3555`), `glpg_ff31`'s stage 7.0 checkpoint promoted by
  `martini_prepare_system.promote_md_output_state_to_input` as `run.py` does, its input `pos` and `mom` bitwise the
  CPU REMD's `remd_start.up`; 49776751 (CPU, about 4,470 s a seed) and 62200421 (V100, 19:01-19:52, about 376 s a
  seed) COMPLETED, all 24 npz written; compared on the midway2 login 22:10 (`g0/v2cmp_glpg_20261103.txt`):
  **passes**, Bonferroni 0.94 over 48 comparisons, largest 2.6 SE. Step 4: `g0/glpg_tm4.sh` (md5 `569f65bd`) wrote
  `g0/glpg_tm4_v100/run.py` (only force field, block count, engine, GPU and wall time differ from the tree's, all
  other lines equal, so the settings are `glpg_ff31`'s) and `run.py` submitted 62200634. Step 5: `checks/mdc/`
  (`mdc_run.py`, `mdc.sbatch`, `submitted.log`), tested on the midway2 login node (2 domains, 20 tu, both force
  fields, `panel._analyze`'s engine check passing); 49777036-43 and 62201376. Ids in `g0/submitted_g9a.log` and
  `checks/mdc/submitted.log`. Step 4's 62200634 COMPLETED 23:03 (2,790 tu, 845 us per system step) and agrees with the
  CPU block (read 23:40: `g0/glpg_tm4_read.py` and `g0/glpg_tm4_rungs.py`, outputs `g0/glpg_tm4_*.txt`). The probe
  62201376 FAILED 22:50 (lesson: a job's temporary directory takes upside's log with it, so `mdc_run.py` now prints
  the log's end on a failure; the cause, about 20 GB of device memory for batch 4, was measured on the Mac's host
  backend); resubmitted on an A100 as 62225254, 23:49, COMPLETED 00:17 11-04 (19,755 MiB, 68.5 us per system step,
  host work 19.1%). g9b (the host loop parallel, plan.md item 12 step 5) built on midway2 and installed 00:36 11-04;
  62227671-79 submitted 00:38. 62227671 COMPLETED 03:46 (`g0/g9b_bitwise_62227671.out`: 21 frames of 8 systems
  bitwise equal; 216 against 237 us per system step). The panel 62227672-79 COMPLETED 04:57-08:29 (1.2-1.9 h a
  batch; 39.4-39.6, 51.4-51.5, 65.2-65.5 us per system step on V100s for batches 1-3, 60.4-60.8 on the A100 for
  batch 4): 356 runs per force field to t 9990, no non-finite value, KE/1.5kT after t 2000 1.006 (1.001-1.011 per
  run), the last frames' longest C-N 1.85 A (ff2.1) and 1.84 A (ff3.1). Folded fraction against the CPU reference
  (`paired_npz.py`, `checks/mdc/paired_u3_<ff>_gpu_cpu.txt`): ff2.1 -0.002 [-0.023, +0.019] p 0.91, ff3.1 -0.003
  [-0.030, +0.023] p 0.92. TM-score 49803361-62 (submitted 08:33) COMPLETED 09:27-10:37 11-04 and **passes**
  (read 10:50 on the MacBook Pro): mean TM, GPU - CPU over the 89 domains, ff2.1 -0.000 [-0.005, +0.005] p 0.86, ff3.1
  +0.000 [-0.005, +0.006] p 0.91 (logs `checks/mdc/logs/tm_mdc_u3_<ff>_gpu_cpu_<id>.out`). No job of line 5 runs.
  The panel's trajectories (`upside3_gpu/checks/mdc/runs/u3_<ff>_{gpu,cpu}/*.up`, 197 GB) were deleted 11-04 after
  this read; npz, logs and the TM series kept ("Disk").
- V1 of the released ff3.1, 61986517 (V100, log `checks/g0/v1_v100_61986517.out`, COMPLETED 00:30 11-01): the three
  configs pass every quantity (384 each), the same worst error per config as 61891693 on the earlier tables; the host
  backend passes them too (Mac Studio).
- Read so far (findings 11g and A239, plan.md Phase 35 items 6-10): V3 on g7b (CPU 49668720-54, one ladder each,
  seeds 1-4; V100 61895691-98, 4 ladders each, seed 101; all 40 in one chunk; data `v3/runs/` 7.4 GB, logs
  `v3/logs/`, analysis 49703078 into `v3/analysis.txt`): Bonferroni 0.12 to 1 per protein and force field, 0.93
  overall, T_mid within 3 K. V1 on g7b (49666580 K80, 61891693 V100, 61891694
  RTX6000, 61891695 A100; logs `v1_<gpu>_<id>.out`): every quantity of the 11 fixtures and the 3 ff3.1 configs passes
  on all four, the same worst error per config on each. The protein V2 by seeds on g7b (49666588-93 cpu and k80,
  61891699, 61892576, 61893444 v100; seeds x001-x016 cpu, x101-x116 k80, x201-x208 v100, x = 2 without rotamer, 3
  full, 4 ff3.1; `v2p_norot/`, `v2p_full/`, `v2p_ff31/`): Bonferroni 0.087 to 1 per config and GPU. V2 of the
  hybrid by seeds (REMD of 1rkl 7.0, 4 replicas, 200 tu; `v2s_1rkl70/`, npz per seed) found the GPUs' lipids 0.35% hot: 40 GPU seeds (49662280, 61856214,
  49662978, 49662979 with nvcc `-fmad=false`) against 24 CPU and host-backend seeds (49662279, 49662631-2), 4.5 SE.
  GPU runs repeat bit for bit (49662440-1); the g-JF coefficients agree (49662465, `coef/`); the forces along a GPU
  run are right (49663123-4); damping and noise agree over 20 rounds (49663131); paired runs diverge as t^2 and
  cannot show heating (49662457). Without exchange (4 copies at T 0.80, `ns_1rkl70/`) the K80 stayed +0.42% hot
  (49663136-7). The cause was the force term rounding away (A239): the K80 with the step summed first (49663236,
  k80fix) +0.1 SE. V2 on g7b (49663908 cpu 1001-1016, 49663909 k80 1101-1116, 61867766 v100 1201-1208;
  `v2s7_1rkl70/`): Bonferroni 1 on both GPUs, lipid KE within 1.8 SE in every replica. The old host backend without
  exchange (49663219, 10 of 16 seeds; 6 hit the 3 h limit on shared nodes): with the CPU (+1.0 SE), not the old K80
  (-2.4 SE). Speed: g6a, g6b, g6c bitwise equal to g5b (49662431-3, 49662459-61, 49662644-6); per system step 245,
  526 us and 6.5 ms (g5b 360, 867 us, 14.5 ms; 49662436, 49662462-4, 49662647); against batch size 49662667-8 (K80)
  and 61859739-40 (V100), findings 11g. Profiles 49662308-10 (nvprof, K80): time in kernels whose one thread loops
  over hundreds of partners (plan.md item 9).
- Earlier: 61818142 (A100, g7b): the probe and the bonded V1 pass. The K80 and V100 probes and V1 runs, and the protein V2
  runs 49650220, 49650478 (K80) and 61822367, 61822883 (V100): ensembles within 2.9 block SE, the GPUs slower per
  system step than 4 CPU cores. 61818143 (RTX6000) passed the probe and the bonded V1. 49661582 (K80, g5a) and
  61849316 (V100, g5a): V1 passes every node of 11 fixtures, the hybrid's included. Hybrid V2 by one seed and block SE:
  49661583 (K80) up to 4.7 SE, 61849317 (V100) up to 5.3 SE at other observables and signs; GPU 15.3 and 4.5 against
  CPU 24.3 and 18.0 ms per system step. Its control 49662130 (broadwl, two CPU runs at seeds 7 and 8) differs by up to
  7.6 block SE (lipid D, replica 3; potential 5.3), so the block SE understates the hybrid's error and V2 now takes it
  from seeds. 49660377 (K80, g4b): every node passes but chig_rama's total potential, whose rounding-scale draw came
  out exactly 0 (findings 11g; the test now takes 4). 61822213 and 61846223 cancelled.

### TM4 HDX on ff_3.1_depth_data: the recipe (done 10-17, findings 5.5; user 10-14 02:10, "Do the analysis on ff3.1 as well")

The HDX readout of one glpG REMD block, as run on ff_3.0_gdepth's block 2 (10-14), ff_3.1_depth_data's (10-17, 79ALA
then 79HIS) and, with BLOCK a parameter, ff3.0's blocks 2-3 (line 2). Scripts in
`C=/project/trsosnic/yinhan/checks/hdx_ff30_20261014/` (md5: `hdx_block2.sbatch` `019af27d`, `prechecks.py`
`7ede2a7a`, `extract_block.py` `0f8781f6`, `readout.py` `5e78ea57`, `martini_hdx_project.py` `8317fbe1`, the repo's;
`plot_tm4.py` `1c91ddbd`; `*.ff30_run_20261014` are the ff3.0 run's own copies).
1. **When:** the chain's block log (`popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<id>.out`) ends with `[remd] block complete`.
2. **The block's groups.** n_k = `grep -c "^avg_kinetic_energy/1.5kT"` in block k's log (one calibration run plus its
   chunks). `output_previous_0` is the seed's own group, so block 2 is FIRST = n1 + 1 to LAST = n1 + n2. Check that
   group FIRST of `run.0.up` has 287 frames whose `time` ends at 18.0 (the 2000-step calibration); read live files only
   with `HDF5_USE_FILE_LOCKING=FALSE`. Measured on ff_3.0_gdepth 79HIS: n = 17, 18, 19, so block 2 = 18-35 and block 3
   = 36-54.
3. **Pre-checks on midway2 login1** (read only, ~2 min): `cd $C; source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh;
   export HDF5_USE_FILE_LOCKING=FALSE; ulimit -v 4000000; nice python3 prechecks.py <ff> glpG-RKRK-<V> FIRST LAST
   <block jobid> | tee prechecks_<ff>_<V>.txt`. It asserts the log's run count equals the group count (ff_3.0_gdepth:
   KE/1.5kT 0.998-1.002, swap acceptance 0.36-0.53, no round trip).
4. **Submit from midway3** (caslake: the only place with mdtraj, pymbar and matplotlib). Tell the user, open the midway3
   socket once (one Duo push), and per variant: `sbatch --job-name=hdx_<V>_<tag> --output=$C/hdx_glpG-RKRK-<V>_%j.out
   --partition=caslake --account=pi-trsosnic --exclude=midway3-0014 --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
   --time=08:00:00 --export=ALL,PDB_ID=glpG-RKRK-<V>,FF=<ff>,FIRST=<g>,LAST=<g> $C/hdx_block2.sbatch`. It refuses to
   run if its work dir exists (~5 min each).
5. **Check:** the log has `block 2 copies verified: 28/28`, `replicas done: 28/28` and `HDX DONE`, no `FAIL rep` or
   `Traceback`. Work dir `popepopg_REMD_mdw2/<V>.<ff>/hdx_sc600x12_<ff>/` (~30 G each; `results/` PS,
   `protein_only/results/` PS_protein, `results/<V>_tm4_readout.txt`). Quote only censored-aware dG (the readout's ">
   bound"), never the standard plots' unclipped values.
6. **Read and report:** TM4 130-151 at 298 K (T 0.85) against the reference run (findings 5.5). Figures: `plot_tm4.py
   NPZ_DIR OUT_DIR` on a computer with SciencePlots (fix its `SP` path), with `<V>_PS.npz` = `results/<V>_dG_profiles.npz`
   and `<V>_PS_protein.npz` = `protein_only/results/<V>_dG_profiles.npz`. The ff3.0 figures are on the Mac Studio in
   `~/Downloads/glpG_TM4_HDX_ff3.0_gdepth_block2/`. Then findings 5.5, plan.md Phase 16, this file and
   `WATCH_STATUS.md`.

### Finished runs whose files are still used

**ConDiv trainings** (`$P/training/<run>`, logs `<run>/condiv-train_<id>.out`; results in findings 1.21, 1.26, 1.30,
1.34, 1.36, 1.37; ff3.1.1 / ff3.2.1 and the reruns are in §1 line 2):
* ff30_gdepth_dt009 (released ff_3.0_gdepth 10-09, §1b; 32 Peng arms complete, scored 10-21); ff31_depth_all (not
  validated); ff31_depth_data (released ff_3.1_depth_data 10-13, retired 10-26; 32 Peng arms complete 10-20,
  their replicas deleted 11-03, "Disk").
* ff32_depth_data, stopped by the user 08:33 10-18 in step 69 (plan.md Phases 17 and 21): trained with the
  angular-gradient defect (findings 3.12); steps 0-68 and panels f32_00-02 are the defect side of the comparison with
  its rerun; readouts to step 68 in WATCH_STATUS.md of 10-18 07:42 and progress.md; the partial step 69
  (`run_output/epoch_03_minibatch_12`, 23 of 24 proteins) is left as it stopped.
* ff21_nobp_dt009 and ff21_ctrl_dt009: step 38 gave the 38-step comparison (findings 1.21), its audit (plan.md Phase
  20) and the releases `ff_2.1_nobp`, `ff_2.1_ctrl38` (step-37 panels and TM4 sets `c9_00`, `c9_01` in findings 1.22
  and 1.24); ff21_ctrl_dt009 is ff30_bio_dt009 with ff2.1's `rama.dat`, so it starts from ff2.1 exactly
  (`checks/ctrl_dt009_init_20261006`). Continued on the fixed engine, gates CONVERGED at step 76; mdCATH panels
  `mdc_nobp_03`, `mdc_ctrl_03`; their 32 Peng arms scored as ff_2.1_nobp / ff_2.1_ctrl38
  (`ff21_nobp_dt009/setup/score/compare_20261021.txt`).
* Retired, kept for comparison: **ff30_basin** (released ff_3.0 09-30, failed; plan.md Phase 4, findings 1.8-1.14):
  gates `ff30b-gate_<jobid>.out`, `gate_step{76,95,114}.txt`, round libraries and `run_output/rama_rounds.txt`, the
  release in `release_20260930-022510/`, the three 09-28 rewinds in `rewound_20260928_{0851,1225,1608}/` (pre-change
  files `$P/training/backup_pre_{basin,ggsheet,reduced}_20260928/`), pre-proline scripts and logs
  `/project/trsosnic/yinhan/checks/prepro_*`. **ff30** (full-map glycine row, cancelled 09-28 at step 223): last
  checkpoint `ff30/run_output/epoch_11_minibatch_12`; extract with `$P/training/backup_pre_basin_20260928/extract_ff.py`.
  **ff30_gly** (stopped in epoch 3): checkpoint `epoch_03_minibatch_07`, resumable; panels e00-e02; run_output synced
  from the local Mac run (marker `run_output/FROM_LOCAL`; the cluster's own as `run_output.superseded_20261001-111755`).
  Its `upside_input/` is a hardlink copy of ff30_basin's except `rama.dat` (a copy of `parameters/common/rama31.dat`,
  md5 `fc479d45...`, built by `training/build_gly_library.py`): never edit a hardlinked file in place. **ff30_glyhb**
  (stopped by the user 10-04 at step 63, release held, findings 1.17): checkpoint
  `$P/training/ff30_glyhb/run_output/epoch_03_minibatch_06`. **ff30_glyprobe** (cancelled 09-30, findings 1.15): no
  `after_training.sbatch`, so nothing can resume or release from it.

**Panels and analyses** (`/project/trsosnic/yinhan/ff3_selection` unless stated):
* Controls `runs/ff21_released/`, `runs/ff21_awh/`; all-atom `aa/` (44 domains); `runs/e00/`; diagnostics
  `runs/e00_hb21/`, `runs/fp_e00/` (the ff21-fixedpoint epoch-0 control) and the single-group resets
  `runs/e00_{bb21,env21,rot21,sheet21}/` (findings 1.17); round 4's baseline `runs/gdepth_start/`. `panel.py select aa
  domains runs ff21_released ff21_awh gdepth_start` (43 domains, 4hwiB01 dropped; folded, then helix / beta /
  gly_helix / gly_left): ff21_released 0.615, -0.018 / -0.017 / -0.081 / +0.029; ff21_awh 0.585, -0.021 / -0.019 /
  -0.081 / -0.065; gdepth_start 0.603, -0.019 / -0.015 / -0.041 / +0.008 (bootstrap SE 0.004-0.005 for helix and beta,
  0.043-0.049 for gly_helix, 0.019-0.028 for gly_left); gdepth_start is not resolved from ff21_released in any class.
  Finished panels' `.up` trajectories were deleted 10-11 and 10-18 ("Disk"); npz, `.up.output` and `params/` kept.
* Analyses of 10-20 and 10-21 (`checks/converge_20261020/`, `checks/ff21_released_20261021/`, `checks/ff311_20261021/`):
  TM-score of the side-chain fixed points (plan.md Phase 26); the epoch-end panels `c9_02`, `c9_03`, `a31_00`-`03`,
  `a32_00`-`03` and `panel44_trajectory_final.txt`; Peng scores of ff3.0 / ff3.1 (`peng_trained.txt`); the single-file
  resets `c9_00_*21` (`resets_c9_00.txt`); the ff3.2 rerun's convergence checks; the shared-width cross-validation
  (`cv_pool_*.txt`).
* The lambda ff_2.1 benchmark (findings 11d; `/beagle3/trsosnic/yinhan/ff3_benchmark/runs/lambda_*_ff_2.1/`, WT and
  G46A/G48A, native and de novo; line 1 reuses two of its arms).

**glpG (Upside 2), stopped or superseded:**
* 79HIS ff_2.1_ctrl38 and ff_2.1_nobp, stopped by the user 08:33 10-18 with STOP files (their releases superseded):
  `popepopg_REMD_mdw2/glpG-RKRK-79HIS.ff_2.1_{nobp,ctrl38}/`, seeds and logs kept.
* ff2.1_sc, cancelled by the user 10-14 in block 2 (block 1 whole): `popepopg_REMD_mdw2/glpG-RKRK-79HIS.ff_2.1_sc/`
  (counters kept), logs `logs/remd.glpG-RKRK-79HIS.ff_2.1_sc.{49223635,49299095}.out`, seed
  `seeds/glpG-RKRK-79HIS.ff_2.1_sc.up`; its HDX was not run. The kill left every replica file readable (Upside traps
  SIGTERM, `main.cpp:1252-1404`; KillWait 30 s). **To stop any chain at a block's end:** `run_remd.py` starts a chunk
  only while remaining - 2100 > 1.25 x the last chunk's wall and reads STOP at that loop's top and before the resubmit
  (lines 190-210), so the last chunk is the one whose remaining R < 2100 + 2.25 x its own wall; take the walls from the
  block in question (block 2's ran 6,745-7,051 s against block 1's ~7,200 s).
* TM4 HDX on ff_3.0_gdepth's block 2 (10-14, the recipe above; readout findings 5.5): launcher
  `checks/hdx_ff30_20261014/hdx_block2.sbatch` (`092d3e29...`), logs `hdx_<V>_<jobid>.out`, pre-checks
  `prechecks_20261014.txt`; results in `popepopg_REMD_mdw2/<V>.ff_3.0_gdepth.placedefect_20261017/hdx_sc600x12_ff_3.0_gdepth/`
  (`results/` with `<V>_tm4_readout.txt`, `protein_only/results/`); HDX topology from ff_3.0_gdepth (md5 `dae228ee...`
  79ALA, `f4c0ea21...` 79HIS).
* The ff_3.0 validation of 09-30 (ff30_basin's release), cancelled by the user after ~18 h because ff_3.0 failed
  (findings 1.14-1.15). Its chains ran ff_3.0's map, H-bond and pair tables with the FF1-form ff_3.0's coverage tables,
  so TM4's glycine flips there belong to that mixture (findings 3.11). Partial data in
  `ff3_benchmark/runs/<prot>_<mode>_ff_3.0/` (trajectories deleted 11-02, inputs and logs kept); superseded
  benchmark `runs_superseded/ff_3.0_20260930-022510` (the same), pre-release seeds `seeds/*.bak_pre_ff_3.0_20260930-022510`. **Each base
  glpG variant directory (`popepopg_REMD_mdw2/glpG-RKRK-*/`) has a `STOP` file**, so `run_remd.py` ends any chain there
  at once: remove it, and reset `block_count`, before validating a new force field there.
* The local glpG TM4 test (one checkpoint patched into the hybrid glpG 79HIS seed, T 0.80, 4000 tu, 12 seeds), retired
  10-11 (user: "We have a working ff and we are trying to expand it so that metric no longer matters and no longer
  needs to be tested locally"): tools in `scratchpad/ff3_local_test/` (method in findings 1.25 and 4.15), record and
  tables in `checks/r4_epochs/tm4_local` (`runs_cov/`; trajectories deleted 10-18); results in findings 1.21-1.26,
  1.31, 4.15. Every set predates the placement fix (findings 4.17), and every set before 10-07 11:22 is invalid (the
  FF1-form coverage tables, findings 3.11). Its total-potential jumps above 3000 (findings 1.24) are not explained.

**Dropped by the user and deleted** (manifests in `checks/disk_reclaim_20261011/`; results in findings 1.21-1.25,
1.31, 1.32, 4.15; the `tm4_compare_cov_*.txt` and `select_last_*.txt` tables kept):
* The three H-bond-frozen arms (ff30_gdepth_fz, ff30_bio_fz, ff21_ctrl_fz): stopped 10-09 ("meaningless at this point:
  the general training at dt 0.009 already made a force field that keeps TM4 stable"), deleted 10-11 ("They are wrong
  and no longer needed").
* ff30_gdepth_si, the SI-rate twin of ff30_gdepth_dt009 (findings 1.23), stopped at step 72 of 76 with nothing
  released, deleted 10-12 ("since ff3.0 is working, it is the only branch we need (and its expansion to all residues).
  delete ff_3.0_gdepth_si").
* The BioEmu-map line (ff30_bio, ff30_bio_dt009 released `ff_3.0_bio`, ff30_bio_si released `ff_3.0_bio_si`), dropped
  10-11 ("we decided to go with basin depth training path (as it is expandable to all residues)"): runs, parameter
  sets in both trees, glpG chains, benchmark runs, panels and TM4 sets deleted. Not touched: the BioEmu glycine library
  itself (repo `parameters/common/rama31.dat`, `checks/gly_bioemu_map/`).
* polygly (plan.md Phase 12), stopped 10-11 ("I don't think poly GLY matters any more as our fix repaired tm4 issue in
  GLY flip") and kept whole (213 MB): `polygly/prod/STOP` written, replicas checkpointed at ~97 ns, resumable from
  `prod.cpt` after deleting STOP; analysis in `polygly/prod_analysis/` (gmx 2024.4 needs `module load gcc/10.1.0`), the
  collapse stage in `polygly/collapse/repN/`.

**midway3: Boltz-2 HDX sampler** (`$B` = `/beagle3/trsosnic/yinhan/boltz_hdx`; readouts in `diffusion_dynamics.md`
Phases 0-1). Every job 1 GPU on `gpu`/QOS `gpu`, `pi-trsosnic`; midway2's `squeue` does not see midway3 jobs, so check
from the midway3 socket. `$B`: `boltz/` at `b1ebfc4` (with the chunk fix to Boltz's sample batching, uncommitted, `git
diff` there; findings 10.14), `.venv` (torch 2.14.1+cu126), `cache/` (weights, CCD), `inputs/` (YAMLs, MSAs,
references, from `prep.py`), scripts `install.sh`, `prep.py`, `check_samples.py` and one sbatch per set (local copies
in `scratchpad/boltz_hdx/`). Logs `$B/logs/<job name>_<id>.out`; data `$B/runs/<job name>_<id>/` (per run `<tag>.log`,
`<tag>.gpu_mib`, predictions under `<tag>/`, `<tag>.npz` with every sample's N/CA/C).
* Phase 0 (10-21): 61216915 A100, 61216916 V100, 61216917 RTX 6000 (`boltz_t0_<gpu>_<id>`, `first_test.sbatch`); the
  tag by step-scale test 61231966 `boltz_tag_stepscale_v100` (`tag_stepscale_test.sbatch`).
* Phase 1, V100, `md` at `step_scale` 1.0, 1,000 samples per Boltz run: EHEE's set `ehee_md_s10` in 61254153
  `boltz_phase1_md_s10_v100`; ubiquitin's set `ubq_md_s10` (chunks of 100) and the chunking check (`ehee_c1a`,
  `ehee_c1b`, `ehee_c10`) in 61279169 `boltz_phase1_ubq_rerun_v100`; Ub psWT's 1,000 in 61298446
  `boltz_phase1_heeh_pswt_v100` (references `inputs/1UBQ_psWT_ref.pdb`, `inputs/HEEH_rd4_0097_ref.pdb`); HEEH's 10^4 as
  10 seeds in 61328353 (TIMEOUT with seeds 3-7) and 61368734 (seeds 8-12), all in
  `$B/runs/boltz_phase1_heeh_seeds_v100_61328353/`.
* EHEE and Ub psWT at 10^5, and the MSA control (Ub psWT with EHEE's query-only MSA: `$B/inputs/ubqps_nomsa.yaml`
  `1a881a28`, `inputs/msa/ubqps_nomsa_0.csv` `53ebcb6c`): 10 jobs per set of 10 seeds (13-112),
  `$B/phase1_seeds.sbatch` through the ledger `$B/submitted_1e5.txt` (ids there) and `$B/submit_1e5.sh` (`8f31186c`);
  each seed packed into one tar by `$B/pack_predictions.sh`, read by `$B/check_samples.py`, every seed `ok`. The MSA
  control's gate passed 10-26 (CA-RMSD median 1.54 A, Rg 11.34 A); read out 10-27: MSA depth is ruled out. Data: each
  ledger line's OUT, else `$B/runs/boltz_phase1_<name>_1e5_v100_<id>/` with its `<name>_md_s10.npz` (10^4 samples) and
  one tar per seed.
* Local copies (md5-equal, `samples/jobs/`) and pooled sets in `scratchpad/boltz_hdx/flags/samples/`; `pool_seeds.py`
  pools a protein's 10 npz in ledger order.

**midway3: BioEmu2 HDX side test** (`diffusion_dynamics.md` section 7, kept apart from the Boltz work). `$E =
/project/trsosnic/yinhan/bioemu2_hdx`: BioEmu2 at `06be369`, weights `microsoft/bioemu2` `dd4b0c1b` in `$E/cache/hf`
with `$E/ckpt_full/` (`trimul_variant: full`), `install.sh`, `prep.py`, `sample.sbatch`, `submit.sh`,
`bioemu2_samples.py`. 61445028 (EHEE and Ub psWT 10^3, HEEH 10^4): log `$E/logs/bioemu2_hdx_v100_61445028.out`, data
`$E/runs/bioemu2_hdx_v100_61445028/<name>_n<N>/` and `<name>_n<N>.npz`. EHEE and Ub psWT at 10^5 (seeds 1000-1900,
through `$E/submitted_1e5.txt`), every frame physical but one of Ub psWT seed 1700's (frame 1622, "isolated_atoms", kept
in `samples.xtc`): logs in `$E/runs/bioemu2_hdx_<ehee|ubqps>_1e5_v100_<id>/`, whose samples were deleted 11-03
(disk section). The 10^5 readout files are at `$E/readout/`, uploaded by the MacBook Pro 16:25 11-03 from its
scratchpad and md5-checked: `flags/{ehee,ubqps}_bioemu2_n100000.npy` (md5 `567aa284`, `5d06181a`), `flags/{ehee,
ubqps}.resid` and `samples/ehee_bioemu2_n100000.npz` (`1fe494fe`). The Mac Studio pulled them md5-equal 11-03 into
`scratchpad/bioemu2_hdx/readout/`, where the readout is. Ub psWT's 10^5 samples npz was not uploaded and is on the
MacBook Pro only.

### Disk: from midway2, `df` on the subdirectory is the only number you get for `/project`

```
df -h /project/trsosnic     -> the FILESET (size, used, free), correct
df -ih /project/trsosnic    -> fileset inodes
df -h /project              -> 6.3P total, the whole device, useless
rcchelp quota               -> home, scratch, project2 group only; no /project or /beagle3 row
mmlsquota                   -> "File system project is not known", the GPFS client is not here
quota   (midway3 login)     -> /beagle3 and /project group blocks and files with soft quota, limit and grace: the
                               only place the soft quota and its grace show
```

`df` on the subdirectory reports the fileset because GPFS `--filesetdf` is on; `df` on the mount point reports the
6.3 PB device. That distinction is the whole trap.

`df` counts free space against the hard limit, but past the soft quota a grace period runs, and when it ends the soft
quota is the limit: creation fails while `df -i` still shows ~118,000 free (findings 10.39).

Now: **`/beagle3/trsosnic` files 656,016 of a 1,268,000 soft quota, 611,984 under** (`df -i` 14:42 11-04, after
`anaconda/`'s deletion below; limit 1,394,800; `quota` lags `df -i` by up to an hour). Ours are about 118,000 (`yinhan/` 121,705
by `du --inodes` 11:03 11-04, about 92,000 of them the venvs and interpreters live work uses); the rest are other
members': `anaconda/` 507,205 (owner uid 1105467662, a removed account), `schwartznw/` 339,897, `yiheng/` 180,347,
about 15,000 in smaller directories, before `anaconda/` was deleted 11-04 (below). **`anaconda/`** (checked 11-04, read only): an Anaconda3 2022.05 base (Python 3.9.18) made under
bayhi's account (gone on midway3), group-writable; no process, job or script outside bayhi's trees uses it now; last
activity a conda command 08-25; members installed into it until 2026-02 (jbolger tensorflow, MDAnalysis; simoneritchey
mrcfile; schwartznw jupyter; ours antiorm, db, sqlitedict). 449 packages, 161 beyond stock 2022.05; site-packages
154,251 files. Ours depends on it only through `/beagle3/trsosnic/yinhan/upside2-md/source.sh` (anaconda first on
`PATH`, used by `run_hdx_analysis.sh`), and for us it is already broken: `~/.local`'s numpy 2.0.2 shadows its 1.26.4,
so pandas, tables and mdtraj fail unless `PYTHONNOUSERSITE=1`. **Its replacement** (user 11-04: public, not under
one member's name; a new venv, since the shared deployment's `.venv` would have to change numpy, tables, scipy and
prody under a running campaign): `source /project/trsosnic/shared/anaconda_py39.sh` on either cluster, a venv on
`python-3.9.18` (a copy of the deployment's `pyrt`) holding anaconda's 426 Python packages at their versions, conda's
tooling left out, about 166,000 files (9.7 GB) on `/project`; differences, build (`build/build.sh`, `build.log`) and
checks in `/project/trsosnic/shared/README.md`. Every package imported in its own process: 420 of 429 on both
clusters (fail: pycurl, which has no Python 3.9 wheel; db, streamlit and a check artifact, as in anaconda), after
gcc 10.1.0's libstdc++ and conda's unixODBC were added to its interpreter's `lib`; on midway2 its tensorflow stack
imports where anaconda's does not. It reproduces the HDX tools' results of `.venv_el8_py311_bak` (md5-equal
`PDB_to_initial_structure.py` outputs, bitwise trajectory, MBAR to 10 decimals) on midway3. **Backup**:
`/project/trsosnic/shared/archive/beagle3_anaconda_20261104.tar` (job 49809806: 23 GB, 637,439 entries equal to the
tree's, no error; a manifest of every path with its access time written first, since reading the tree resets them;
`tar --compare` found no difference in the 65 min it ran before it was stopped for the deletion), restore by
`tar -xpf <tar> -C /beagle3/trsosnic`. **Deleted 14:03-14:42 11-04** (user: "delete it now, as we already have tar
file"; job 49812826, `shared/archive/rm_anaconda.sbatch`, log `rm_anaconda_49812826.out`), after a last check that
no process or job on either cluster used it and nothing in it changed since the tar: `/beagle3` files 1,163,216 to
656,016. Left: `etc/jupyter/labconfig/page_config.json` and its directories (schwartznw's, not group-writable).
`/beagle3/trsosnic/yinhan/upside2-md/source.sh` still puts the deleted path on `PATH`; remove that line when the
ff3.1 Peng campaign ends. Blocks 1.2 T free by `df` after the 11-04 reclaim. It was over the soft quota from 10-27 at the latest (IFree 125,183 then) until the grace expired 11-02
between 19:30 and 19:46, and no file could be created until 21:51. Blocks 4.39 T of 5.00 T. **`/project/trsosnic`
2.92 T of a 3.49 T soft quota** by `quota` 14:47 11-04, before ff3.0's deletion; 1.2 T free by `df` after it
(14:58; limit 3.84 T; files about 661,000 of 1,035,800, the shared environment's 166,000 included). `/project2/trsosnic` 1.45 T of 1.49 T; nothing in these campaigns writes there. Midway3 home: 25 G of a
30 G quota (§8).

**`/beagle3` inode exhaustion, 12:20-14:30 10-24** (findings 10.39): the first submission of Boltz's 10^5 tests filled
the inode quota. Jobs writing there failed: Upside 3's glpG REMD could not write `block_count` (`OSError: [Errno 122]
Disk quota exceeded`), and four Peng chunks died in 9 s with no log (exit 0:53: Slurm could not create the log); both
were resubmitted (§1 line 3).

**ff3.0 and ff3.1 are transitional** (user 11-04): "once ff3.2 is complete and validated, all ff3.0 and ff3.1
results can completely gone". Until then only what is analysed and read by nothing running or owed is removed. Left
for that day: the trainings in `$P/training` (`ff31_nativefix`, which ff3.2's Rg check reads; `ff311_glypair`, whose
step-34 replicas `crystal_view.py` reads at each `ff32_basin` epoch end; `_mb2`, `ff31_depth_data_angfix`,
`ff311_depths_only`, `ff321_depths_only`, `ff31_depth_all`, `ff31_depth_data`; about 29 GB and 91,000 files), the 9
running ff3.1 Peng arms, every remaining ff3.0/ff3.1 npz and score, and the released `ff_3.0_gdepth` / `ff_3.1` files.

**Reclaimed at the user's word** (manifests, scripts and logs kept):
* 11-04 15:13, ff3.1's analysed data, 129.9 GB and 1,606 files (`checks/disk_reclaim_20261104_ff31/`, `m13`-`m16`;
  user: "ff3.1 is also transitional ... check how far ff3.2 is, and if any existing ff3.1 can be analyzed and
  removed"; `/project` 1.2 T to 1.3 T free). Readouts first to `checks/archive_1102_20261104/ff31/` and the 1102
  report folder (the panels' TM tables, the paired folded/Q tables of `p44_ff31nf`, `mdc_ff31nf` and `a31_00`-`03`
  against released ff2.1, the 23 Peng scores, glpG's preparation and block 1 TM4). Deleted: the replicas of the 23
  scored ff3.1 Peng arms (18.3 GB; `input/`, logs, `COMPLETE`, scores kept; the 9 running arms untouched),
  `p44_ff31nf` / `mdc_ff31nf`'s `.up` (60.1 GB; npz and TM series kept), Upside 3 `glpg_ff31`'s REMD replicas (8.4
  GB; checkpoints kept) and `a31_00`-`03`'s `.up` (43.1 GB; npz kept; no TM-score, as `ff21_released` has no
  trajectories left to pair with). Then `a32_00`-`03`'s `.up` (704 files, 43.0 GB; user: "delete a32"; `m17`,
  tables in the same archive), and the 9 running ff3.1 arms kept (user: "keep the current running jobs").
* 11-04 14:58, ff3.0's trajectories, 260.3 GB and 1,002 files (`checks/disk_reclaim_20261104_ff30/`, `m8`-`m12`,
  `reclaim.log`; user: "ff3.0 is transitional ff, so once we have the results analyzed and on the slides, we can
  comfortably delete the trajectory files"; `/project` 972 G to 1.2 T free). First the one owed readout, TM4 over
  blocks 1-5 (above), then the readouts copied to `checks/archive_1102_20261104/ff30/` and the 1102 report folder
  (TM4 tables, Upside 3 `glpg_ff30` blocks 1-3, the 12 HDX `results/`). Deleted: the four Upside 2 chains' replicas
  (192.2 GB; `seeds/`, logs, `block_count` kept) and their HDX `outputs/popg_remd` (21.1 GB; `results/`, `inputs/`
  kept), the `.placedefect_20261017` HDX outputs (2.8 GB), Upside 3 `glpg_ff30`'s REMD replicas (9.4 GB;
  checkpoints and `remd_start.up` kept) and the `ff_3.0_gdepth` Peng replicas (34.6 GB; `input/`, scores kept).
  The released force field files `ff_3.0_gdepth` themselves are kept.
* 11-04 14:03-14:42, `/beagle3` files 1,163,216 to 656,016: `/beagle3/trsosnic/anaconda` deleted (507,200 files),
  backed up first as one tar on `/project` and replaced by `/project/trsosnic/shared/anaconda_py39` (above).
* 11-04 11:42-11:44, 559 GB and 7,420 files: `/beagle3` 294.5 GB (3,544 files), `/project` 264.8 GB (3,876)
  (`checks/disk_reclaim_20261104/`: `reclaim.sh` md5 `81723eb8`, `m1`-`m7` manifests, `reclaim.log`; user: "Make sure
  we have the results analysed and recorded in /Users/yinhan/Documents/2026/reports/GroupMeetings/1102, then proceed
  to delete results already done analysed"). Readouts first copied to `checks/archive_1102_20261104/` and to the
  report folder's `results_archive_20261104/` (README there: dataset, result, what was deleted). Deleted, trajectory
  `.up` files only, each with one link and unwritten for an hour: the GPU mdCATH panel and its CPU reference
  (`upside3_gpu/checks/mdc/runs/`, 197.4 GB), the GPU glpG TM4 block's replicas (8.6 GB; `remd_start.up` and the
  stage 7.0 checkpoint kept), V3's runs (7.6 GB), the `ff_2.1_nobp` and `ff_2.1_ctrl38` Peng replicas (69.1 GB;
  `input/`, `run.log`, `COMPLETE`, scores kept), `bp_validation/runs` (11.8 GB; `static/`, `obs_*.npz`, analyses
  kept), the panels `c9_00_{bb21,env21,hb21,rot21,sheet21}`, `c9_02`, `c9_03`, `n32_00`-`04`, `s32_00`,
  `z32_00_{rot21,bb21,hb21,env21,sheet21,rama0}` (204.7 GB) and the retired merge's `p44_ff31`, `mdc_ff31` (60.1
  GB); every panel's npz kept, so `paired_npz.py` and `panel.py select` still read them. The first run stopped after
  `m1` (its check loop returned 1 under `pipefail`); fixed and rerun, 0 of every manifest remain. Kept as not
  finished or still read: the Upside 3 Peng runs (top7's gap open), the released ff2.1 Peng arms, the reference
  panels (released ff2.1, released ff3.1, ctrl38 and no-BP epoch 3), `w32_*`, `z32_00`-`03`, `a31_*`, `a32_*`, the
  ff3.0 Upside 2 glpG chains and the Upside 3 glpG runs.
* 11-03 08:05-08:12, `/beagle3` files, 1,254,102 to 1,163,475 by `df -i`
  (`checks/disk_reclaim_20261102/reclaim_boltz_phase2.sh` md5 `7a15ee01`, `manifest_boltz_mols.txt`,
  `manifest_boltz_venv.txt`; user: "Delete some files to clear the file limit", "Maybe tar some of them then delete",
  then chose "Boltz CCD molecule cache" and "Boltz Phase 2 venv"). Deleted `boltz_hdx/cache/mols/` (45,228 files,
  1.9 G), the extracted copy of `cache/mols.tar` beside it (`tar --compare`: only the owner differs, same count;
  `mols.tar` kept). Packed `boltz_hdx/.venv` (45,422 files) into `boltz_venv.tar` there (7.3 G, compared equal and
  counted), then deleted it. **Before any Boltz E1 rerun, restore both** (`e1_rare.py` reads `cache/mols` directly;
  90,650 files back): `tar -xf /beagle3/trsosnic/yinhan/boltz_hdx/cache/mols.tar -C
  /beagle3/trsosnic/yinhan/boltz_hdx/cache` and `tar -xf
  /project/trsosnic/yinhan/checks/disk_reclaim_20261102/boltz_venv.tar -C /beagle3/trsosnic/yinhan`. Offered and
  kept: the HDX analysis venv `upside2-md/.venv_el8_py311_bak` (46,248 files).
* 11-03 00:50-00:57, 486 G and 2,300 files, `/project` 3,468 to 3,094 GiB and `/beagle3` 1,392 files
  (`checks/disk_reclaim_20261102/reclaim_ff3x_20261103.sh` md5 `7035382e`, `manifest_ff3x_20261103.txt`; user: "We
  can delete some ff3.0 and ff3.1 trajectories to free some space", then chose "Retired ff3.1 sets" and "ff3.0 HDX
  frame copies"; the manifest checked first: every name `*.run.<r>.up`, no hardlink, none written in 24 h). Deleted:
  the replica files of the three retired ff3.1 Peng benchmarks in `ff3_benchmark/runs_superseded/`
  (`ff_3.1_bak_pre_nativefix_20261031` 30.6 G, `ff_3.1_depth_data_bak_pre_merge_20261026` 32.2 G,
  `ff_3.1_trained_bak_pre_merge_20261026` 1.9 G; each arm's `input/`, `run.log` and `COMPLETE` kept); the REMD
  replicas of the merge's Upside 3 glpG `upside3_validation/glpg_ff31.bak_pre_nativefix_20261031` (7.8 G; its
  preparation checkpoints, `remd_start.up` and readouts kept); the 3-minute step 0 of the cancelled
  `ff311_depth_pooled.cancelled_20261021` (25.3 G; `base.up` and logs kept); the frame copies `block2/`, `block3/`,
  `block2-3/` and `joined/` of the 12 ff3.0 glpG HDX work dirs (388.5 G; `results/`, `outputs/popg_remd/`,
  `inputs/`, `pdb/`, `protein_only/`, `acc_shift_corrected/` kept, 29 G; checked after: 12 of 12 `dG_profiles.npz`,
  336 of 336 outputs). The four ff3.0 chains' replicas (112 files, 179 G) stay, so the recipe's
  `hdx_block.sbatch` (`extract_block.py`) rebuilds any copy into a new work dir. Offered and kept then: those chains
  (TM4 at blocks 2-5 owed) and the released `ff_3.0_gdepth` Peng arms; both deleted 11-04 (above).
* 11-03 00:24-00:38, old Boltz and BioEmu2 (`reclaim_boltz_bioemu2.sh`, `manifest_boltz_runs.txt`,
  `manifest_bioemu2.txt`; user: "clean up old boltz and bioemu2 files, if no longer in use", then chose BioEmu2's env
  and weights, its samples, and Boltz Phase 0-1 runs; both finished and read out, `diffusion_dynamics.md` Phase 1 and
  section 7): `/beagle3` `boltz_hdx/runs/` 23.8 G and 7,231 files (39 runs, each run's `.log` and `.gpu_mib` kept);
  `/project` `bioemu2_hdx/` 11 G and about 43,000 files (samples of its 21 runs, logs kept; `.venv`, `python/`,
  `tools/`, `cache/hf`, rebuilt by its `install.sh`; `cache/msa` kept). Boltz Phase 2 (`boltz/`, `cache/`, `.venv`,
  `addon/`, `inputs/`) untouched. The Mac Studio keeps the pooled Phase 1 sets (`scratchpad/boltz_hdx/flags/samples/`)
  and BioEmu2's 10^3/10^4 sets (`scratchpad/bioemu2_hdx/readout/samples/`); the query-only-MSA Ub psWT 10^5 set is
  gone. BioEmu2's 10^5 flags and EHEE's 10^5 samples survived on the MacBook Pro and are now at
  `$E/readout/` (uploaded 16:25 11-03; BioEmu2 side-test paragraph); Ub psWT's 10^5 samples npz is on the MacBook Pro
  only. The read A100 probe's replicas deleted 00:34.
* 11-02 22:56, `/beagle3`, 44 G and 896 files (`reclaim_ff30_peng.sh`, `manifest_ff30_peng.txt`; user: "some ff3.0
  trajectories can be cleaned", then "Failed + superseded Peng"): the replica `.up` files of the cancelled 09-30
  ff_3.0 validation's 32 Peng arms (`ff3_benchmark/runs/*_ff_3.0`) and of `runs_superseded/ff_3.0_20260930-022510`;
  each arm's `input/`, `run.log` and `COMPLETE` kept. The ff3.0 glpG chains, their HDX work dirs and the released
  `ff_3.0_gdepth` Peng arms were offered and kept.
* 11-02 22:40, `/project` blocks, 180 G (`checks/disk_reclaim_20261102/reclaim_project.sh`, manifest
  `manifest_ff31_depth_data_chains.txt`; user: delete the "Retired ff3.1 glpG chains"): the 112 replica `.up` files of
  `popepopg_REMD_mdw2/<V>.ff_3.1_depth_data.bak_pre_merge_20261026/` (the old ff3.1 on Upside 2, retired 10-26, every
  readout done); their `block_count`, `block_jobid`, logs, seeds and the `*.placedefect_20261017.bak_pre_merge_20261026`
  directories kept.
* 11-02 21:51-21:54, `/beagle3` files, from 1,276,989 to 1,263,519 (`checks/disk_reclaim_20261102/`:
  `reclaim_beagle3.sh`, `reclaim.log`, a manifest per target; "do a cleanup of the disk, check what is no longer
  used"): deleted `ff3_benchmark/sigmoid_discarded` and `armM_discarded` (§1b "Safe to delete", 1,409 files, 7.6 G)
  and the read V100 glpG probe's replicas (7.8 G); packed into tars there, each compared with its directory and
  counted before removal, `chi1_perplexity` (done 10-12, 5,148 files, `chi1_perplexity.tar`), `pylibs` (4,322,
  `pylibs.tar`; the FF1 `.xtc` read's mdtraj, findings: unpack and add to `PYTHONPATH` to use it) and the GPU tree's
  builds `obj_g1a`-`obj_g6c`, `obj_hb6b` (2,542, `gpu_builds_pre_g7b.tar`, 3.6 G). Restore any by `tar -xf <tar> -C
  /beagle3/trsosnic/yinhan`.
* 10-11, from 88 G free to 1,041 G (`checks/disk_reclaim_20261011/`; "delete the trajectory file but try to keep the
  setup scripts so I can regenerate them"): the 39 finished panels' `.up` trajectories (397 G; `select` reproduces
  b9_03's table byte for byte); NP-1AO6's run data (`prod_ff3`, `prod_ff3_armM_20260910-084115`, `prod_counterion`,
  `tscan`), keeping every top-level script (`build_np_ff3.py`/`.sbatch`, `build_tscan.py`, `np_hybrid.py`,
  `run_np_prod.py`, `np_prod.sbatch`, `submit_np.sh`, ...), though `prod_counterion/onset.py` and `onset2.py` had no
  copy and are gone; the four `glpG-*.ff_3.0_cancelled_20260930` directories' replica files and counters.
* 10-18, 1.1 T (`checks/disk_reclaim_20261018/`: `build_manifests.sh`, `delete_manifests.sh`, `delete.log`): half-written
  steps of stopped trainings (86 G; ff3.2's partial step 69 left); the replica `.up` files of the defective glpG runs
  and their HDX `joined/` and `block2/` copies (821 G: the 8 `.placedefect_20261017` chains, `popepopg_REMD/*_pre_cbfix`
  and `*_rigidbug`, the 79HIS ff_2.1_nobp, ctrl38 and ff_2.1_sc chains, `armtest_{M,R}`, `tm4_local/runs` and
  `runs_cov`; HDX `results/`, `outputs/`, `inputs/`, `pdb/`, `protein_only/`, STOP files, counters, logs, seeds and
  tables kept); the `.up` trajectories of 14 finished panels (228 G: `runs/{ctrl38,nb_01,fa_01-03,fd_00-03,f32_00-02}`,
  `runs_mdcath/{mdc_nobp,mdc_ctrl38}`; `panel.py select` reproduces `select_last_f32_02_login1_20261017.txt` byte for
  byte, but CA-RMSD (`panel_compare.py`) can no longer be recomputed for these tags).

**Large directories awaiting the user:** `popepopg_REMD`'s four variant directories (355 G, `hdx/` on the §1b list),
`gly_peptides` (§1b), `glpG_DDM_micelle_REMD` (73 G, §1b), `implicit_79HIS_*` (70 G), `popepopg_glpG` (9 G). What
grows: each panel ~11 G (mdCATH ~47 G), each
training step in flight ~38 G of replica files (deleted at the step's end), each block-2 HDX ~28 G.

---

## 1b. Data and decisions that outlive their jobs

* **ff_3.0 is ff30_gdepth_dt009's release** (user, 10-09 10:40): `$P/parameters/ff_3.0_gdepth` and
  `/beagle3/trsosnic/yinhan/upside2-md/parameters/ff_3.0_gdepth` (from `ff30_gdepth_dt009/release_20261009-085341`, =
  `checks/r4_epochs/d9_03`), copied to the repo's local `parameters/ff_3.0/` (md5 verified). The cluster names stay as
  they are.
* **Data that exists only on the cluster, do not delete:**
  * AWH glycine surfaces, `/project/trsosnic/yinhan/gly_peptides/` (rep1, `awh_amber99sb-ildn_rep2/`,
    `evidence_diffusion_bug/`).
  * `popepopg_REMD_mdw2/run_remd.py` and `NP-1AO6/run_np_prod.py` exist only there, with no git history.
  * `popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`, the pre-ff3 TM baseline.
  * HDX: keep `popepopg_REMD/<V>/hdx/` (pre-fix). `hdx_postfix/` and the earlier ff3 results at
    `popepopg_REMD_mdw2/<V>/hdx_10k/` are on neither `/project` nor `/beagle3` (searched to depth 5, 2026-09-30);
    `hdx_results/` holds only the DDM-campaign figures and MBAR files.
* **Deletion decisions waiting on the user:** midway3 DDM data `glpG_DDM_micelle_REMD/`, `glpG_DDM_REMD/`; rockfish
  post-fix glpG data `/scratch4/rherna21/ywang268/upside2-md-rf/popepopg_REMD/<variant>/` (~72 GB, needed only for an
  independent-seed cross-check), before scratch is reaped; `upside2-md-water-diff` (§0b).
* **Open, unverified, recorded so they are not lost:**
  * The HDX topology uses ff_2.1 tables on an ff_3.0 trajectory.
  * `NP-1AO6/np_footprint.py` has CB-offset and minimum-image bugs, with no record of a fix.
  * The one-shot upgrade line in `NP-1AO6/np_prod.sbatch` is still to delete.
  * An NP orientation at T = 0.70-0.75 was never run.
  * A latent stride-3/4 issue was noted at `upside_config.py:815` `_input_phi`.
  * The HDX resubmit block uses `HDX_WORK=.../hdx` without `HDX_N=28` and would overwrite the pre-fix baseline: fix it
    before any HDX rerun.
* **NDRD library files must never be copied to the cluster** (licence).
* **Do not overwrite the cluster's `py/` or `training/` from the local repo during these campaigns.** The repo's
  `training/` is now five files: `extract_ff.py`, `convergence_gate.py` and `rama_basin.py` became `ConDiv.py extract` /
  `gate` and functions, and `train_chain.sbatch` runs the gate itself instead of submitting `after_training.sbatch` ->
  `gate_or_continue.sh` (plan.md Phase 9). The runs in `$P/training` need the midway2 tree's old layout ("Where things
  are"), so the new files go to the cluster only for a new run. The midway2 tree also still runs on the campaign files
  the repo moved to `scratchpad/redistribution_cleanup_20261001/`: releases use `$P/training/{validate_ff.sh,patch_glpg.py}`,
  the glpG HDX jobs (`popepopg_REMD/hdx_*.sbatch`) call `$R/py/martini_remd_concat.py`, and `NP-1AO6/build_np_ff3.py`
  imports `martini_inject_coverage`. The repo's `training/env.sh` and `train_chain.sbatch` are site-neutral (no midway2
  module loads, no `--account`, since 10-02); the midway2 tree and every run directory keep their own copies with both.
* **The midway2 tree's training scripts as deployed:**
  * Glycine offsets since 2026-10-02 (`checks/hbg_deploy_20261002`, logs `deploy.log`, `finish.log`, `parity.log`):
    `src/hbond.cpp`, `py/upside_config.py`, `training/{ConDiv,extract_ff,check_step,patch_glpg}.py`,
    `training/README.md` and `obj/{upside,libupside.so}`, each with a `.bak_pre_hbg_20261002`; on 1ga3 the new build is
    bitwise the old one with ff2.1's 12-entry file. The `/beagle3` shared deployment was not updated with it on 10-02;
    it was to be synced from `$P` before ff3.0's benchmark or glpG runs.
  * dt 0.009 (user, 2026-10-06; backups `.bak_pre_dt009_20261006`): `$P/training/ConDiv.py` `time_step=0.009`,
    `train_chain.sbatch` STEPS_PER_LINK 54 (90 at dt 0.015); the repo's the same.
  * Round 4 (backups `.bak_pre_gdepth_20261005`, `.bak_pre_bio_20261005`): `extract_ff.py` writes a depth-training
    checkpoint's trained library and finds the run's own trainer for an initial checkpoint; `check_step.py` skips an
    empty glycine-offset field and prints the glycine depth.
  * `training/ConDiv.py` is the FF2 trainer (findings 9v); the FF1-form trainer and helpers are in
    `$P/backup_training_ff1form/`. The offset training was removed 10-01 (superseded files, `verify_rama_basin.py` and
    `build_rama_from_awh.py` among them, in `$P/backup/training_pre_glyfix_20261001`); old run directories keep their
    own ConDiv copies, so the analyses of ff30_basin and ff30_glyprobe still work.
  * **`validate_ff.sh` and automatic validation** (user, 2026-10-05; backups `.bak_pre_r4val_20261005`): a converged
    gate (`gate_or_continue.sh` case 0) calls `validate_ff.sh <run> <ff_name> <checkpoint>` (an epoch end only); a
    validation failure fails the gate job. It extracts the checkpoint through the run's own `expand_param`, backs up
    and overwrites the force field in `$P/parameters` and the `/beagle3` deployment (md5-verified), moves the superseded
    benchmark directories to `runs_superseded/`, submits the 32 Peng arms (12 h chunks, `runs/<tag>/input` made first)
    and the 4 glpG chains: it gates the patch on a pristine ff_2.1 seed, patches the campaign's `seeds/<V>.up` (left
    untouched) into `seeds/<V>.<ff>.up`, runs variant `<V>.<ff>`, and stops before anything is released if `<V>.<ff>`
    exists. Every submission takes the run's `slurm.args`, and glpG follows the partition (broadwl
    `popepopg_REMD_mdw2`; caslake `popepopg_REMD_mdw3`, its `run_remd.py` and `env.sh` md5-identical to mdw2's, its
    `remd.sbatch` sourcing `env_shared.sh` then `env.sh`, its `submit_remd.sh` on `--partition=caslake
    --exclude=midway3-0014 --mem=56G`). Per candidate ~95k core-hours, ~165 G on `/project` and ~50 G on `/beagle3`.
    Tested in `checks/r4val_20261005/`. **Submitted from midway3 it fails** ("not every benchmark arm was submitted":
    the release happens, validation does not): release from midway2. Choosing which candidate becomes a release stays
    the user's.
  * `ff3_benchmark/bench.sbatch`: `BENCH_RESUBMIT` carries the job's own `ExcNodeList`. `bench_run.py` lost its type-0
    burial override and `rama3.dat` fallback 09-24 (backup `bench_run.py.bak_ff1form_20260924`), so it cannot
    benchmark the old FF1-form ff_3.0.
  * Panels from `submit_new.sh` go to broadwl (`ff3_selection/slurm.args`); one sent to caslake needs `--mem`
    (`panel.sbatch` sets none). The gate (`after_training.sbatch`, 2 CPUs) needs none.

---

## 2. THE TWO CAMPAIGNS ARE DIFFERENT SIMULATIONS

Conflating them has caused several errors, including a threshold copied from NP that killed a
healthy 6 h glpG block. **Never transfer settings, thresholds, or analysis between them.**

| | **NP** (`np_1AO6_prod`) | **glpG** (`remd_glpG-*`) |
|---|---|---|
| method | regular MD, 6 independent trajectories, single T=0.8647, no exchange | **REMD**, 28 replicas, T ladder 0.70–0.90, configuration exchange |
| purpose | nanoparticle adsorption footprinting (K190 exposure) | **HDX** protection factors / ΔG |
| system | 1AO6 albumin 578 res + 5 nm MPA-AuNP, 8608 atoms, box 300 Å | glpG 210 res in a POPE/POPG bilayer. **Read the atom count and box from the seed**: two generations exist, 4949 atoms / 279 lipids / box 99.77² × 180 Å and an older 4709 / 261 / 99.869² × 123.697 Å |
| composition | PROTEIN 2890 + GOLD 887 + MPA 203 + ION 4628 (K+ 2423 / Cl- 2205, 0.15 M KCl) | PROTEIN 1050 + LIPID (13 beads each) + ions regenerated at 0.15 M; the counts follow the seed generation |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: ions, lipids and the 630 protein N/CA/C sites are on the single-stage g-JF **Brownian** path (4529 of 4949 atoms on the current seeds); the other 420 protein atoms are not in `/input/brownian` |
| timestep | **0.001**, freely settable at runtime | **0.009, HARD-LOCKED** by `/input/brownian/numerical_time_step`; `martini_brownian.cpp:100` throws on mismatch. Friction is tuned against it for lipid D=11.5 µm²/s; **do not change it** |
| detection | non-finite positions OR ≥5 stretched bonds | non-finite potential (whole chunk) OR ≥5 stretched bonds |

---

## 5. How to check TM helix health (TM1 and TM4)

Run this after the first block of any new seed generation or force field, and at every block after.
TM1's C-cap (GLY49) and TM4 (N-cap loop 131-133) are where the old glycine bias failed.

```bash
cd /project/trsosnic/yinhan/checks
source /software/modules/init/bash; module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1
export HDF5_USE_FILE_LOCKING=FALSE
python3 glpg_tm_windows.py      # TM1 30-48 / TM4 134-151 per rotated group, T 0.70 and 0.80, seed first
python3 glpg_ff30_health.py     # also GLY49/GLY133 phi and stretched C-N by residue and replica
python3 gly_tm4_flip.py         # fraction of frames with TM4's helical glycines at phi > 0
```

Both TM scripts read only rotated `output_previous_*` groups (never the live `/output`) and skip
`output_previous_0`, the seed's own block. Pass: TM1 and TM4 mean helix fraction > 0.8 at T 0.70
(phi in [-130,-20], psi in [-90,15]), with no downward trend across blocks.

**Three rules that earlier versions of this section got wrong** (findings.md, memory
`glpg-vtf-reading-traps`):
* **Dihedral sign.** Use the negated-atan2 form of `check_tm4_ss.py` (as in the two scripts above),
  validated against BioPython. The helper formerly printed here and `check_seeds_current.py` return
  -phi_std: a healthy helix reads as unfolded.
* **Windows.** The crystal seed has 131-133 and 152 non-helical, so 131-152 caps TM4 at 0.818; TM1's
  29-49 caps at 0.952. Use 134-151 and 30-48.
* **No glycine-phi criterion.** GLY49 and GLY133 are helix caps at phi_std +94 and +142 in the crystal
  seed; positive phi there is native. `check_seeds_current.py` (map mirror symmetry of those two
  glycines) is invalid twice over: sign-flipped, and the GLY|X maps are not mirror-symmetric
  (ff_3.0's by training, rama31.dat's by measurement), so it reports such a seed BROKEN. Do not gate
  a submission on it.

---

## 5b. How to check health CORRECTLY (general bond/energy check)

**`isfinite` is not a health check.** At a real glpG failure the environment coordinates were
±4.65e12 Å, numerically finite and physically destroyed. And in a forced NP tear the protein reached
431 broken bonds with the potential still finite at +3e5, so no energy-based test fires at all.

Use the broken-bond count (torn 279-431). The "healthy 0-2" figure holds for the cold rungs only.
Measured on glpG 2026-09-30 (`checks/glpg_cn_control_20260930.log`, every 5th frame): at T 0.70 no
frame has 3 or more C-N above 2.0 A, at T 0.80 ~0.1% do, and at T 0.88-0.90 1.5-2.7% under ff_3.0
(`inner_steps` 4, max 10 in one frame) against 1.8-5.1% in the pre-ff3 campaign (`inner_steps` 1,
max 13). So local transient tearing at the hot rungs predates ff_3.0; it is an open integrator
question for the hybrid, not a force-field regression, and it does not accumulate (sampled frames
return to 0, and the logs had no ROLLBACK). The Peng arms (pure Upside) never exceed 1 in any frame.

```python
import h5py, numpy as np
def n_broken(up, group="output", frame=-1):
    with h5py.File(up, "r") as h:
        pos = np.asarray(h[group]["pos"][frame, 0])
        if not np.isfinite(pos).all(): return -1          # -1 => non-finite
        nm  = np.array([x.decode() for x in h["/input/atom_names"][:]])
        pc  = np.array([x.decode() for x in h["/input/particle_class"][:]])
        rid = h["/input/residue_ids"][:]
    prot = np.where(pc == "PROTEIN")[0]
    C = {int(rid[a]): a for a in prot if nm[a] == "C"}
    N = {int(rid[a]): a for a in prot if nm[a] == "N"}
    res = [r for r in sorted(C) if (r+1) in N]
    cn = np.linalg.norm(pos[[C[r] for r in res]] - pos[[N[r+1] for r in res]], axis=1)
    return int((cn > 2.0).sum())
```

Cheap whole-chunk scan (glpG failures go fully NaN, so the scalar potential catches them):
```python
v = np.asarray(h[group]["potential"][:]).reshape(-1);  bad = int((~np.isfinite(v)).sum())
```
Also useful: protein Rg, and `avg_kinetic_energy/1.5kT` at the end of a log (healthy ≈ 1.0).

---

## 6. The detection gate and rollback (what "DESTROYED" / "ROLLBACK" in a log means)

- NP `health()`: non-finite positions OR `n_stretched >= NP_CN_COUNT` (5) → **ends the chain**
- glpG `destroyed()`: non-finite potential anywhere in the chunk OR `n_stretched >= REMD_CN_COUNT` (5) → **rolls back that replica** (as of 2026-08-12 driver)

**NP**: a gate trip looks like `[np] DESTROYED ...` followed by `no resubmit`; the job exits
rc=0/COMPLETED short of its wall limit. A COMPLETED job that did not resubmit means the gate fired;
check the log.

**glpG**: a rollback looks like `[remd] ROLLBACK #N filename: reason` followed by
`rolled back M/N replicas; continuing chain`. The chain does not terminate. The NaN output is rotated
to `output_previous_N` as normal history; the rolled-back replica restarts the next chunk from its
pre-chunk positions. A replica that repeatedly blows up is rolled back repeatedly, not dropped from
the ladder, so a high rollback count on the same file indicates a persistent physics problem.

**Rollback mechanism**: before each chunk `run_remd.py` snapshots `/input/pos` of every replica.
On NaN detection it overwrites the last `output/pos` frame and `output/potential[-1]` with the
pre-chunk values so that `reseed()` on the next iteration picks up the clean state.

**These driver scripts are not in git.** They live on the cluster at
`~/project/yinhan/popepopg_REMD_mdw2/run_remd.py` (midway2) and `~/project/NP-1AO6/run_np_prod.py`,
with no version history. Edit directly on the cluster. A running job keeps the version it loaded at
start; edits take effect at the next block.

---

## 7. If the chain terminates: manual rollback procedure

**glpG (if the chain stops rather than rolling back):** patch the last output frame of each NaN file
with the last finite frame from `output_previous_0` (end of block 1), then resubmit.

```python
import h5py, numpy as np
from pathlib import Path
run_dir = Path("~/project/yinhan/popepopg_REMD_mdw2/<variant>").expanduser()
for fn in sorted(run_dir.glob("*.run.*.up")):
    with h5py.File(str(fn), "r") as h5:
        n_bad = int((~np.isfinite(np.asarray(h5["/output/potential"][:]))).sum())
        has_prev = "output_previous_0" in h5
    if n_bad:
        with h5py.File(str(fn), "r+") as h5:
            prev_pos = np.asarray(h5["/output_previous_0/pos"][-1, 0, :, :])
            last = h5["/output/pos"].shape[0] - 1
            h5["/output/pos"][last, 0, :, :] = prev_pos
            h5["/output/potential"][last, 0] = 0.0
# then: bash ~/project/yinhan/popepopg_REMD_mdw2/submit_remd.sh <variant>
```

**NP**: `run_np_prod.py` `reseed()` is idempotent, so a config with no `/output` but a valid
`/input/mom` (`restart_valid=1`) restarts fine. To reset a chain: `echo 0 > prod/block_count`, `rm -f prod/STOP`.

---

## 8. Known lessons (environment)

Slurm and filesystem:
* **A Slurm controller fault breaks self-resubmitting chains silently** (2026-11-01, 21:19 onward on midway2):
  srun and sbatch answered "Unexpected message received". A Peng chunk then ends `exit=0` and still prints
  `[bench] resubmitted:`, but followed by `sbatch: error: ...` instead of a job id, and nothing follows it (30 of
  ff3.1's 32 arms); a training link's workers all fail at once, and only its queued insurance link can carry on, once
  the controller answers. Read
  the id on every `[bench] resubmitted` line, not only that the line exists.
* **A training chain insures itself** (the design that worked unattended on ff30_basin): every link queues its
  `afterany` successor before it trains; `--no-requeue`; the successor resumes from the newest checkpoint after
  deleting the half-written step; a worker whose `srun` never started is relaunched (twice at most); three links from
  the same step stop the chain. The gate job itself has no successor.
* **A dead Pool child hangs a training step until the link's wall** (2026-10-29, `ff311_glypair` step 34): the
  worker's `compute_divergence` runs in a `multiprocessing.Pool`; when a child dies (there, glibc aborting on a
  double free), Pool starts a replacement and `pool.map` waits forever for the lost task, so the srun step stays
  RUNNING at 0% CPU and the chain's insurance never fires. Seen as one protein with no `divergence.pkl` while the
  step's other workers finished hours earlier; the child's message is in `<protein>.output_worker`. The way out:
  copy the protein's step files aside (the rerun overwrites them), then `scancel <link>.<srun step>` for that
  worker alone; the master records `WORKER_FAIL`, the link exits FAILED within a minute, and the insurance link
  reruns the step from the previous checkpoint (done 09:20 10-29).
* **A 54-step link times out when steps run long** (2026-10-07): `STEPS_PER_LINK=54` assumes ~36 min a step;
  ff30_gdepth_dt009's took 40-115 min, so its first link hit the 36 h wall in step 35. The successor resumes from the
  newest `checkpoint.pkl`; the half-written step is lost, and the run waits in the queue as at any link end.
* **Size wall-time requests to the work; long requests do not backfill on broadwl** (2026-10-01). Two 1-node panel
  jobs asking 24 h were estimated to start 1.5-2 days out; cut to 4 h with `scontrol update JobId=<id>
  TimeLimit=04:00:00`, they started within the hour. 36 h benchmark arms were estimated 2.5 days out, so the benchmark
  runs 12 h jobs with ~10.5 h chunks (`bench.sbatch`, `bench_run.py` `CHUNK_SEC`; backups `*.bak_pre_12h_20261001`).
* **A script that submits another job must unset `SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU` and `SLURM_MEM_PER_GPU`**
  after its `#SBATCH` block (or match the child's memory-request type), or every `srun` of the child dies with `... are
  mutually exclusive`. Findings 10.1.
* **Slurm snapshots the batch script at submission.** A requeue reruns it with its original arguments, and an edit does
  not reach queued jobs; a running job builds its resubmit command from Slurm's spooled copy, so only a directive in
  the file at submission reaches the next chunk. Chain scripts carry `--no-requeue` (check `scontrol show job <id> |
  grep Requeue`) and resolve the newest checkpoint at run time; after an edit, cancel and resubmit the queued
  successors (`scontrol write batch_script <id>` shows what a queued job will run). Findings 10.2.
* **Cancel a chain's PENDING successor before its running link**, then check `squeue`; the other order starts the
  successor.
* **A link killed by SIGBUS on several nodes at once is a transient `/project` outage**: 0-byte `*.output_worker`
  files, no traceback, and a successor that dies in its first second. Resume; exclude nothing. Findings 10.3.
* **One login node's `/project` can hang while the other's works** (2026-10-07 10:48): on midway2-login2 every lookup
  under `/project/trsosnic/yinhan` sat in D state (`wchan` `lookup_slow`) while login1's answered and the jobs kept
  writing. `timeout` cannot kill a D-state process and `ps` on it hangs too, so each probe stacks another stuck command.
  Test with `timeout 15 ls <dir>` once, then route `/project` work through login1: `ssh ... midway2.rcc.uchicago.edu
  'ssh -o BatchMode=yes midway2-login1 bash -s' < script.sh` (host-based, no Duo), and rsync with `-e "ssh -o
  BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu ssh -o BatchMode=yes"` and host
  `midway2-login1`. login1's `/project` was slow too (a 70 s panel select took 4 m 44 s).
* **A failed worker fails the step** (`ConDiv.py` raises): step 81 of ff30 once ran on 5 of 24 proteins after `srun:
  Job credential expired`, and a partial step is a different objective. That credential error at launch is a race,
  not a node fault (24 steps issued together lost 4 on healthy nodes); the trainer relaunches a worker that never
  started, up to twice. When both relaunches land on one node the link ends; the step is discarded and the successor
  resumes from the step before. That happened twice on midway2-0088 ("Invalid job credential", 10-07 14:40 and 10-09
  03:42), its record for the law (§0d).
* **Node exclusions are per record of failures** (the list lives only in the law, §0d, each node's record in its
  `README.md`). The broadwl NODE_FAILs of 10-02 to 10-09 were one fault on a fixed set of nodes: jobs of several users
  died in the same second, nearly all holding one of the 9 broadwl nodes whose slurmd started at 2026-10-01 18:10:54
  (midway2-[0010-0011,0033-0035,0037,0060,0080,0085]). Those nodes compute; their link to the controller drops, so the
  controller declares NODE_FAIL and kills what runs there. An apparent 8 h 10 min period was disproved. The 9 nodes
  and 0003 were excluded 10-03, and 0103 and 0116 after the same fault on 10-09 (user's word each time). A group event
  (many nodes in the same second) is not one node's record.
* **`sbatch --test-only`'s start estimate is no guide to the real wait**: it predicted 13:36 on both clusters, then
  the midway3 submission sat PENDING (Resources) while midway2 started at once.
* **midway3 gives a new job no start estimate for its first 20 minutes**: `bf_min_age_reserve=1200` in its
  `SchedulerParameters`, so backfill neither reserves for it nor fills `squeue --start` (N/A, reason Priority) until it
  has pended 1200 s. Not a sign the job is unschedulable.
* Training steps ran 300-450 s slower than the engine accounts for (2026-09-29), a different slowest worker each time;
  most likely `/project` I/O.
* **`/beagle3` was badly degraded for small-file writes on 2026-09-28**: 200 one-line files took 51 s from midway3 and
  503 s from midway2, against 0.06-0.27 s on `/project` and home; a `pip install torch==2.6.0` into the shared venv
  took ~2 h for that reason. Its inode quota can fill (10-24, §1 "Disk").
* **Do not tell the clusters apart by `/software/modules/init/bash`**: it exists on midway3 too. The cluster
  `training/env.sh` tests whether the tree venv's interpreter exists (verified 2026-09-28: midway2 gets the tree venv,
  midway3 the `/beagle3` venv, both torch 2.6.0+cpu, numpy 1.23.5, scipy 1.13.1, tables 3.8.0). Never pipe `source
  env.sh`: a pipe runs it in a subshell and the environment is lost.
* **Delete a half-written minibatch directory before resuming**: a stale `divergence.pkl` is reused.
* **Hand-check success-only end-of-chain branches**, which never run until the end.
* **Pre-create `runs/<tag>/input` before a mass submission** on `/beagle3` (a makedirs race); a partial aggregate mean
  is not a result.
* **Midway3 home quota is 30 G** (21 G used on 2026-09-28, 28.6 G on 09-23). Jobs can fail oddly if home fills.
* **`sacct` on midway3 can hold zombie RUNNING rows** (53233848 and 53233852 from 2026-08-11 still show under `-S
  <today>`, while `sacct -j` says `COMPLETED`); `squeue` is the authority.
* **Exclude midway3-0014.**

Login nodes and shells:
* **`/tmp` on a login node is per node, and reconnects land on different ones.** A `nohup` job launched from
  midway2-login2 writes a `/tmp/<log>` that login1 cannot see, so a later check reports the job gone. Put anything you
  will read again on `/beagle3` or `/project`, and prefer `sbatch` over `nohup` for work that must outlive an ssh
  session. Do not run scripts from `/tmp` either.
* **A stray module in the working directory shadows the stdlib.** Another user's `/tmp/inspect.py`, or an `inspect.py`
  left in a scratch directory, broke `import numpy` with a circular-import error. Name scratch probes something that is
  not a stdlib module.
* **`pgrep -f <name>` matches your own shell command.** `ssh host 'pgrep -f bm_final && echo running'` reports
  "running" because the remote `bash -c` line contains the string, so a dead job looks alive (two false reports on
  2026-09-12). Match the interpreter (`ps -u $USER -o cmd | grep python3`), check the output file, or use a Slurm job
  id.
* **The shared venv is Python 3.9 / NumPy 1.23**: `str | Path` needs `from __future__ import annotations`.
* **midway2's site GROMACS is unusable**: use `/project/trsosnic/yinhan/gmx2024` (`gmx_build/build_gmx.sbatch`). The
  HDX pipeline runs only on midway3 (`.venv_el8_py311_bak`, with pymbar and matplotlib); a renamed venv's
  `bin/activate` still hardcodes its old `VIRTUAL_ENV`.

Simulations:
* **NP dt hard limit: 0.001.** dt=0.005 caused backbone blow-ups during unfolding (large-amplitude spring instability
  at t>250, proven by A/B). Never raise above 0.001.
* **Do not transfer thresholds between NP and glpG.** A CN_MAX borrowed from NP false-positived on a healthy glpG
  chunk (2.52 Å vs healthy max 2.659 Å) and cost a 6 h block.
* **glpG NaN propagation via REMD exchange.** A single blow-up in one replica spreads to every replica via exchange
  within ~60 steps (IEEE 754: `NaN < 0.f` = false). The NaN cascade fix (`!isfinite(lboltz_diff)`) and the per-chunk
  rollback driver address this.
* **A green exit code means nothing** for a self-submitting REMD job. Check the log for DESTROYED/ROLLBACK counts and
  verify physical observables.
* **glpG requeues do not advance `block_count`** (fixed 2026-09-30 13:20). `remd.sbatch` has no `--no-requeue`, so a
  NODE_FAIL requeue reran the job under the same id and `run_remd.py` counted it as a new block; it now keeps the
  block's job id in `<V>/block_jobid` and increments only for a new id (sandbox-tested: fresh, requeued and successor
  jobs; backups `*.bak_pre_requeue_20260930`). The aborted attempts' output is kept as short `output_previous_N` groups
  (100-287 frames), all finite (`checks/glpg_fragments_20260930.log`).
* **glpG: patch the LIVE seeds** (they carry `inner_steps = 4` on `/input/brownian`); prove the method by an ff_2.1
  round trip on a pristine seed. A release must remove each ~150 GB variant directory before resubmitting it.
* **glpG REMD: skip `output_previous_0`** (the frozen seed block). Do not gate a resubmission on
  `check_seeds_current.py` (§5).
* **Identify a `.up`'s force field by least-squares match of its baked tables to `sidechain.h5`** (scale 1.000000),
  not by hash. `martini.h5` comes from ff_2.1 by design.
* **Live dynamics are shown by internal Kabsch CA-RMSD growth** (frozen = 0.000 A), never by reading `current_stage`;
  read only `output_previous_*` of a live `.up`, with `HDF5_USE_FILE_LOCKING=FALSE`.
* **Runs branched from a shared checkpoint do not test reproducibility**; branch at step 0. Independently trained force
  fields are different Hamiltonians and are never pooled.
* **NP:** the box must exceed the molecule's maximum extent plus the 12 A cutoff; `NP_DT = 0.001`; the Rg in
  `np.<jobid>.out` is not minimum-imaged (measure adsorption by backbone contacts with the MPA shell, atoms 3750:3950);
  a restart frame must pass a stricter test than detection; `NP_TEMP` is Peng's FF2 calibration and must be redone for
  a new force field.
* **HDX figures:** keep continuous off-scale excursions and `Y_LIMITS=(-20,30)`, and never render censored amides as
  carets or bounds (user-directed). Cluster HDX scripts drift from the repo; re-upload before trusting results.
* The glpG detergent campaign (closed 2026-08-13, model retired 2026-09-09) used a constant `--seed`, so a rollback
  re-ran the identical failing chunk; the driver now takes a per-chunk seed.
