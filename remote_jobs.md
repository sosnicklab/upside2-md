# Remote jobs on midway2/midway3: status and handbook

**Current state (2026-10-27 23:10; midway3 read 10:50 10-28 from the MacBook Pro).** Owner: none since 21:55 10-27, when the MacBook Pro's session (owner from 12:46
10-26) handed off and quit ("Handoff"); no watch runs. A Mac Studio session is active; the user's "continue
monitoring" there takes the jobs ("Resume here"). ff3.1 means the merge of ff3.1.1 and ff3.2.1, released 10-26 as
`ff_3.1` and under validation; the old ff3.1 and ff3.2 sets are retired (`.bak_pre_merge_20261026`; §1 line 2).
Running on midway2: released ff2.1's last Peng arm, ff3.1's other 29 Peng arms, the ff3.1 training debug (from 23:05
10-27), Upside 3's validation (7 Peng de novo arms, the 8 Peng seed replicates, the glpG hybrid with ff2.1, ff3.1 and
ff3.0) and two collectors. Queued on midway3: the Boltz-2 E1 tuning job 61701305 (the pilot 61663101 passed). Each line of work has a collector, a
dependency chain or an item under "Owed by hand". `WATCH_STATUS.md` on the cluster is newer than this file whenever a
pass ran after it. Job ids, logs and data paths are in §1.

The user's four lines of work (10-21), and what each has on the clusters:
1. **Deterministic vs probabilistic side chain** (plan.md Phases 20-27; findings 1.21, 1.37). Both trainings are
   finished and converged; the deterministic side chain folds less on every held-out test, and both retrained force
   fields fold less than released ff2.1. Released ff2.1 on the same three tests decides whether ff2.1 itself stays: the
   44 domains and mdCATH are done (released ff2.1 ahead of both on each); Peng 31 of 32 scored, ubiquitin de novo
   running. The email to Tobin and John waits on these numbers.
2. **Ramachandran map training** (plan.md Phases 14-17, 21, 28-30, 32; findings 1.34, 1.37). ff3.1 (the merge) is
   level with released ff2.1 on both soluble panels and above the control on both; Peng 3 of 32 scored. The ff3.1
   training debug (Phase 32: does the trainer recover ff3.0's glycine pair) runs two arms to step 76. Finished,
   awaiting review: the ff3.1 and ff3.2 reruns on the fixed engine, whose maps barely moved and which, as every run
   trained from ff2.1 with all terms free, fold 0.15-0.24 less than ff2.1. The 8 Upside 2 glpG REMD chains on ff3.0 and
   the old ff3.1 finished all five blocks 10-24; ff3.0's TM4 checks are owed.
3. **Upside 3 design, debug stage** (upside3.md; plan.md Phase 24 and its "Resume here"). The audit's closing criteria
   are met locally (ledger A1-A231). Its validation runs on midway2 from the Upside 3 tree
   `/beagle3/trsosnic/yinhan/upside3` since 22:26 10-22 (§1 line 3): Peng's 16 proteins with ff2.1 (25 of 32 arms
   complete), the Peng seed replicates, and the glpG hybrid with ff2.1 (REMD block 3), ff3.1 and ff3.0 (block 1 each).
   The TM4 test is answered (10-27, findings 3.15): TM4's loss under ff2.1 is ff2.1's own.
4. **Boltz-2 HDX sampler** (diffusion_dynamics.md). Phases 0 and 1, the BioEmu2 side test and the controlled 10^5
   tests are finished and read out (§1 "Finished runs"). Phase 2's E1 pilot passed 10-28; E1 runs (user 10-28 "run
   E1"), its tuning job queued on midway3 (§1 line 4).

Open and not in the user's list:
- ff3.0's four Upside 2 glpG chains: TM4 at the ends of blocks 2-5 is owed (§1 line 2). ff3.0's HDX on blocks 2-3 is
  done (findings 5.6); ff3.1's was dropped by the user 10-23, and the user closed the TM4 HDX question 22:15 10-23
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
| ff30_gdepth_dt009 (`d9_EE`) | ff3.0: NDRD plus one trained glycine depth pair, every term trained | released as ff_3.0_gdepth 10-09 |
| ff31_depth_data (`fd_EE`) | the old ff3.1: 1,442 per-map basin depths, prior 0.2, every term trained | released as ff_3.1_depth_data 10-13; retired 10-26 (`.bak_pre_merge_20261026`) |
| ff31_depth_all (`fa_EE`) | ff3.1 with all 5,040 depths, prior 0.3 | finished 10-13, not validated |
| ff32_depth_data (`f32_EE`) | ff3.1's depths started from ff3.0 | stopped by the user in step 69 10-18 |
| ff31_depth_data_angfix (`a31_EE`), ff32_depth_data_angfix (`a32_EE`) | the ff3.1 and ff3.2 runs on the fixed engine | finished step 76 10-21, awaiting the user's review |
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
   on midway2 (the MacBook Pro's last pass was 21:32 10-27). Tell the user what started, finished or failed since.
4. **Owed by hand** (`WATCH_STATUS.md` is newer if the two differ; details in §1):
   - Upside 3 glpG with ff3.1 and ff3.0: TM4, KE and C-N at T 0.70 of REMD block 1 when each ends (49573354,
     49589829; against Upside 2's ff3.0 0.96). ff3.0's block 2 49634605 waits on its block 1 (afterok): never run
     `glpg_ff30/run.py` by hand before block 1 ends, since with no `block_count` it would submit a second block 1 over
     the replica files (§1 line 3).
   - Upside 3 glpG with ff2.1: REMD block 2's readout (C-N, DSSP; KE read); block 3 49584589 runs.
   - Upside 3's Peng de novo arms: rescore against Upside 2 as they finish (`peng/score/submit_scores.sh`, idempotent).
   - The Peng seed replicates (49634596-49634603): when all 8 are COMPLETE, score them and read the seed spread against
     the engine gap (§1 line 3).
   - ff3.2 and ff3.1 on the fixed native target (target 76; ff3.2 first): links 49635669 and 49635327 PD, released
     10:58 10-28; nice 10 on ff3.1's next insurance link when it starts; each round end (§1 line 2; plan.md Phase 33).
   - The ff3.1 training debug, 49634971 and 49634972 (target 76): at each round end, glycine's pair against ff3.0's
     by plan.md Phase 32's rule (§1 line 2).
   - ff3.0's four Upside 2 glpG chains: TM4 at the ends of blocks 2-5 (§1 line 2).
   - Read by the collectors and reported when they appear: ff3.1's `peng_ff31.txt` (3 of 32 scored 22:12 10-27) and
     released ff2.1's `peng_ff21.txt` (31 of 32).
5. **Find the owner** ("Handoff"): `cat` the cluster's `checks/r4_epochs/tm4_local/WATCH_STATUS.md`. A move between
   computers (user, 10-26: "be prepared to be shut down and continue this session on another computer"; done Mac
   Studio to MacBook Pro 12:46 10-26, and the MacBook Pro's handoff 21:55 10-27, owner none): the user's "continue
   monitoring" on the next computer is the take-over. If `WATCH_STATUS.md` names an owner whose pass is under 2 h old,
   ask the user to confirm that session is closed before starting a watch.
6. **Take the jobs:** start the watch ("Starting the watch"; no submitter is needed), then rewrite `WATCH_STATUS.md`
   with this computer as owner.
7. **Open questions for the user.** Nothing is blocking. Awaiting the user's word, not to be done without it:
   - **Boltz Phase 2:** requesting the Ferrari et al. data (the user only) matters only if E1 passes (E1 runs since
     10-28, §1 line 4).
   - The force field for the Upside 3 hybrid HDX (findings 3.15), and the DDM micelle route (plan.md Phase 31).
   - The large-protein de novo proposal (plan.md Phase 22 review).
   - The large directories still on `/project/trsosnic` (§1 "Disk").
   - The deferred HDX post-processing fixes (plan.md Known Errors; findings 5.5).
   - The BP validation deck for Tobin (§0c), unsent.
   - **Upside 3 training terms** (`upside3.md` 4.1; plan.md Phase 23 item 5): M1 is done (findings 5.7); M2-M6 and T1
     are planned and run later (user 10-28: "record the plannings. We will run it later"). M2 is ready once ff3.1's
     ubiquitin benchmark arms are COMPLETE, but it also waits for the user's word.

### Handoff: two computers, one owner (2026-10-26)

**Owner now: none** (since 21:55 10-27, user: "Record this session's progress and I will quit and work on another
computer"). The MacBook Pro owned the jobs 12:46 10-26 to its last pass 21:32 10-27 and deleted its watch cron at the
handoff; the Mac Studio owned them 23:55 10-24 to 10:29 10-26 and has an active session now. The next computer's
"continue monitoring" takes them ("Resume here" step 6).

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
  `scratchpad/hdx_term_20261028/` M1's scripts and readout, `upside3.md` 4.1); Claude's memory directories (`memory/`,
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

**One watch at a time.** A watch is a session cron: it dies with its Claude session and expires after 7 days. No watch
runs since 21:55 10-27 (the MacBook Pro's cron `4ff7a1bc`, 12:46 10-26 to its last pass 21:32 10-27, was deleted at
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
is one read-only probe, `~/watch_probe/watch_mdw2.sh` on midway2 (md5 `0098c032` since 10-26, old
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
2. midway2: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu 'bash ~/watch_probe/watch_mdw2.sh' < /dev/null` (read-only probe, md5 0098c032). The probe lists only the last 80 minutes: if the "Last pass" line of /project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md is older than that, also run `sacct -X -n -P -u yinhanw -S <that last-pass time> -o JobID,JobName,State,End,ExitCode` on midway2 and read every job that left the queue since. Look for: any job that left the queue not COMPLETED (read its log tail: Peng arms /beagle3/trsosnic/yinhan/upside3_validation/peng/logs/<prot>_<kind>_ff_2.1_<id>.out for u3_*, /beagle3/trsosnic/yinhan/upside3_validation/peng_seeds/logs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>_<id>.out for the seed replicates u3s*/u2s*, /beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<kind>_<ff>_<id>.out for b_ff_2.1_* and b_ff_3.1_*; collectors /project/trsosnic/yinhan/checks/ff31_merged_20261026/logs/ and checks/ff21_released_20261021/; Upside 3 glpG REMD logs /beagle3/trsosnic/yinhan/upside3_validation/<glpg|glpg_ff31|glpg_ff30>/outputs/martini_glpG-RKRK-79HIS_hybrid/REMD/glpG-RKRK-79HIS.remd_block<k>.log); a Peng chunk COMPLETED without a "[bench] resubmitted" or "[bench] COMPLETE" line; on "[bench] COMPLETE", the arm's final frames (`source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh; export HDF5_USE_FILE_LOCKING=FALSE; nice python3 -I ~/watch_probe/final_frames.py <run dir>`: every coordinate finite, KE/1.5kT near 1.00); error lines in the newest REMD block log of each of the three glpG runs (the probe reads glpg and glpg_ff31; for glpg_ff30 `grep -ciE "error|nan|traceback|abort"` on its newest block log); a REMD block ending (block_count rising, or a block job leaving the queue: ff2.1 block 3 49584589, ff3.1 block 1 49573354, ff3.0 block 1 49589829 with block 2 49634605 queued behind it): report it; its readout is owed by hand (remote_jobs.md §1 line 3: TM4 DSSP, KE and C-N at T 0.70, against Upside 2's ff3.0 0.96); peng_ff21.txt or peng_ff31.txt existing (report its table); /project/trsosnic under 200G free.
3. midway3: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu 'df -i /beagle3/trsosnic | tail -1; squeue -h -u yinhanw -o "%i %j %T %S %r"' < /dev/null`. The live midway3 jobs are the E1 jobs of remote_jobs.md §1 line 4 (job table); when one leaves the queue, report its State and the `[e1] ... FAILED` lines of its log /beagle3/trsosnic/yinhan/boltz_hdx/logs/<job name>_<id>.out (the readout is owed by hand on the MacBook Pro). Report IFree under 30000.
4. Record: when a job finished, failed or changed state, update remote_jobs.md §1 (the line it belongs to and the job table; snapshot time) briefly, and set the "Last pass" line of /project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md on midway2 to this pass's time (owner <this computer>, watch cron in this session). Do not record unchanged state. After any record change, refresh the handoff bundle's copies from the repo root so a session on another computer reads them (remote_jobs.md "Resume here"): `rsync -a -e "ssh -o BatchMode=yes -S $HOME/.ssh/cm-mdw2.sock" remote_jobs.md plan.md findings.md progress.md diffusion_dynamics.md yinhanw@midway2.rcc.uchicago.edu:/project/trsosnic/yinhan/checks/handoff_20261026/records/`.
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
`--partition=broadwl --exclude=midway2-[0003,0010-0011,0027,0033-0035,0037,0060,0080,0085,0088,0103,0116,0151,0342-0345]`.
Each node's record of failures is in `README.md` beside it; §8 has the history. The later additions, each at the
user's word: midway2-0088 10-09 04:12, midway2-0103 and 0116 10-09 10:30, and midway2-0151 10-13 10:29 (two NODE_FAILs
10-10); backups `*.bak_pre_0088_20261009`, `*.bak_pre_0103_20261009`, `*.bak_pre_0151_20261013`. Each was checked
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

Snapshot **2026-10-27 23:10 CDT** (midway2 from `squeue`; midway3 as of 21:32, after which its socket master ended and
was not reopened). Owner none since the MacBook Pro's handoff 21:55 ("Handoff"). Nothing failed since 00:35 10-27.
`/project/trsosnic` 527 G free and `/beagle3/trsosnic` 125,183 inodes free (22:25 10-27). Paths: `$P` =
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
    finite, KE/1.5kT 1.000-1.007 over 14 replicas. Ubiquitin de novo, the last, runs (49630881 at 23:10 10-27).
  - Collector `ff21_collect` (every 6 h; next pass **49633770**, 00:47 10-28): scores each finished arm into
    `ff21_nobp_dt009/setup/score/ff_2.1/` and, with all 32, writes `checks/ff21_released_20261021/peng_ff21.txt` ending
    DONE (released ff2.1 against the control, no-BP, ff3.0 and ff3.1, beside Peng's published FF2). Then plan.md Phase
    27 step 4.
- **Both side-chain trainings are finished** ("Finished runs"); no jobs.

**Line 2: Ramachandran map training**
- **ff3.1, the merge of ff3.1.1 and ff3.2.1, released 10-26 and under validation** (user 10-26: "merge into ff3.1 now
  and run validation. Remove any existing ff3.1, ff3.2, ff3.1.1 and ff3.1.2 (you can keep a backup, but from now on
  ff3.1 refers to the merged version)"; plan.md Phase 30; findings 1.34).
  - Merge on the Mac, `scratchpad/ff31_merge/merge_ff31.py` (output `merge_ff31.txt`, offsets `merge_ff31.npz`): the
    mean of the two runs' round-6 offsets in the NDRD frame, written by the runs' own `rama_basin.write_library`;
    `rama.dat` md5 `5b66c78a`, the other five files ff_2.1's byte for byte. Glycine figures
    `~/Downloads/gly_rama_ff31_nonhelix_zero.png`, `gly_rama_ff21_ff30_ff31_page.pdf`: pooled map alpha_L -2.574,
    alpha_R -1.510 (ff2.1 -3.022 / -1.562, ff3.0 -2.505 / -1.669).
  - Released by `checks/ff31_merged_20261026/release.sh` (log `release.log`) to `$P/parameters`, the shared
    `/beagle3/.../upside2-md/parameters`, the Upside 3 tree's `parameters` and `ff3_selection/params`, md5-verified. It
    moved aside with suffix `.bak_pre_merge_20261026` (104 moves): `ff_3.1_depth_data` and `ff_3.1_trained` in both
    Upside 2 trees (to `<tree>/backup/`), the parameter sets `a31_*`, `a32_*`, `f32_*`, `fa_*`, `fd_*` (to
    `ff3_selection/params_superseded/`), their 64 Peng runs (to `ff3_benchmark/runs_superseded/<ff>_bak_pre_merge_20261026/`),
    `score/ff_3.1_depth_data`, and the old ff3.1's 8 glpG chain directories and 8 seeds in `popepopg_REMD_mdw2`; also
    the MacBook Pro's scratch `rama.dat` copies. Training run directories stay as provenance.
  - Validation: the soluble tests released ff2.1 ran, on the same pre-audit binaries (`checks/ff31_merged_20261026/submit.sh`,
    ids in `submitted.log`; 10:03 10-26). The glpG validation runs on Upside 3 only (user's choice; line 3).
    - Panels done, no error line: `p44_ff31` (176 runs) and `mdc_ff31` (356); TM jobs `tm_ff31.sbatch` (outputs
      `tm_<panel>_<tag>_vs_<ref>.txt`). p44 against released ff2.1: mean TM 0.686 / 0.681, +0.006 [-0.001, +0.012] p
      0.036, TM > 0.5 0.896 / 0.886 p 0.25, 18 domains lower and 26 higher; against ctrl_03 +0.037 [+0.027, +0.047] p
      1.4e-10, 3 lower and 41 higher. mdCATH against released ff2.1: 0.733 / 0.727, +0.005 [-0.001, +0.012] p 0.032,
      TM > 0.5 0.974 / 0.972 p 0.39, 40 lower and 49 higher, only the mixed class moving (+0.012 [+0.002, +0.022]);
      against ctrl_03 +0.031 [+0.024, +0.037] p 1.1e-12, 13 lower and 76 higher. ff3.1 is level with released ff2.1
      on both panels and above the control on both.
    - Peng: 32 arms `b_ff_3.1_<prot>_<kind>` (runs `ff3_benchmark/runs/<prot>_<kind>_ff_3.1/`, logs
      `ff3_benchmark/logs/`). COMPLETE and scored: BBA native (10:21 10-27), BBL native (17:39), BBA de novo (19:39);
      final frames finite, KE/1.5kT 0.999-1.008 over 14 replicas (BBA de novo's longest C-N 1.68 A). The other 29 at
      23:10 10-27 (d de novo, n native): alpha3d d 49633005, n 49633194; BBL d 49633003; cspA d 49634607, n 49633753;
      gpW d 49634594, n 49634065; homeodomain d 49634608, n 49634611; hyp d 49634614, n 49634718; lambda d 49634854, n
      49634664; NTL9 d 49634855, n 49634858; NuG2 d 49634861, n 49634860; proteinB d 49627508, n 49627503; proteinG d
      49627505, n 49626771; proteinL d 49627504, n 49627496; top7 d 49627507, n 49627506; ubiquitin d 49629329, n
      49628060; WWdomain d 49629721, n 49629328.
    - Collector `ff31_collect` (every 6 h; next pass **49634812**, 04:07 10-28; logs `checks/ff31_merged_20261026/logs/`):
      scores into `score/ff_3.1/` and writes `peng_ff31.txt` (ff3.1 against released ff2.1 and ff3.0), ending DONE.
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
- **The ff3.1 training debug** (user 10-27: "setup your debug for ff3.1 training"; plan.md Phase 32; findings 1.34):
  ff3.1.1 as built (`build_depths_only.sh`, every Adam group held) with glycine parameterised as ff3.0 is, its shared
  alpha_R and alpha_L only (`restrict_gly.py`); expected, if the training is sound, to reach ff3.0's pair (alpha_R
  -0.107, alpha_L +0.517). `ff311_glypair` (ff3.1.1's minibatch order) **49634971** and `ff311_glypair_mb2` (its own
  order, same 456 proteins and held-out set) **49634972**, submitted 23:03 10-27 to target 76 through
  `train_chain.sbatch` (15 nodes a link). `ff311_glypair` R from 23:05 10-27 (successor 49634973 afterany, the chain's
  insurance link); `_mb2` R from ~23:30 10-27 (successor 49635107). At 01:32 10-28: 5 and 2 steps done, about 30 min
  a step, 24 of 24 proteins every step, no error line. Run dirs `$P/training/<run>`, logs
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
  **Owed:** TM4 at the ends of blocks 2, 3, 4 and 5 of ff3.0's four chains (`tm4_ff31.py`, `gly_flip_rungs.py`, DSSP).
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
  `results/glpG-RKRK-<V>_dG_profiles.npz`. **Keep every HDX work dir** (user 10-22: "do not cleanup for now"): the 12
  hold 417 G. Their membrane term carries the A61 accessibility shift (plan.md Known Errors; findings 5.6): 61383990
  `hdx_accfix` re-ran only the dG step with the shift undone (`$C2/accfix_ps.py` md5 `b5186e2d` into `<work
  dir>/acc_shift_corrected/results/`, then the pipeline's `4.calc_D_uptake.py`; `$C2/accfix_readout.py` `27be5c1c`
  writes `$C2/tm4_accfix_blocks2-3.txt`). Peptide 140-144's amides 142-144 stay censored at 298 K in 79ALA and
  79ALA_S115T, agreeing with HXMS on that peptide in 79ALA; Val145 comes down to 4.02 (79ALA). The corrected npz's
  seven amides before a proline or the chain end (labels 13 27 58 61 127 151 209) keep the shifted value, so no figure
  is drawn from it.

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
  ff2.1 on the pre-audit binary, line 1's arms. COMPLETE, 25 of 32: all 16 native, and BBA, WWdomain, homeodomain,
  top7, lambda, alpha3d, proteinB, NTL9 and BBL de novo. Running de novo at 23:10 10-27: cspA 49629882, gpW 49630880,
  hyp 49630574, NuG2 49630879, proteinG 49629732, proteinL 49630878, ubiquitin 49631645. No collector scores them.
- **Peng scored against Upside 2** (findings 3.15): `peng/score/score_compare.py` (md5 `e17d977f`; reproduces released
  ff2.1's scores exactly), `score.sbatch`, `submit_scores.sh` (idempotent; jobs `sc_cmp_<upside3|upside2>_<arm>`, ids
  in `peng/score/submitted.log`) into `peng/score/ff_2.1_upside3/` and `ff_2.1_upside2/`, read by `compare_u3_u2.py`
  (`compare_u3_u2.txt`). Native, 16 proteins: Upside 3 0.564 against Upside 2 0.552 (+0.012 [-0.001, +0.028], p 0.56).
  De novo over 9 pairs: -0.025 [-0.051, -0.003], median -0.005, carried by top7 (-0.108) and alpha3d (-0.071), one run
  per arm per engine; reweighting each engine's frames by the other's Hamiltonian predicts -0.002 (scripts and outputs
  `/project/trsosnic/yinhan/checks/u3_vs_u2_20261027/`, whose `pull_cold.py` + `score_local.py` reproduce
  `score_compare.py`). Rescore the de novo arms as the rest finish.
- **Peng seed replicates**, the test of that gap: top7 and alpha3d de novo with ff_2.1, seeds 101 and 202, on Upside 3
  and on Upside 2 (the original arms ran seed 1), 8 arms `u<3|2>s<seed>_<prot>_denovo` **49634596-49634603**, submitted
  and R from 20:32 10-27 from the Mac Studio (`submit_seeds.sh`, idempotent; targets 2.15 M and 3.04 M time units). In
  `upside3_validation/peng_seeds/` (logs `logs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>_<id>.out`, runs
  `runs/upside<3|2>_<prot>_denovo_ff_2.1_seed<S>/`, ids `submitted.log`). `bench_run.py` (`c79004af`) is
  `peng/bench_run.py` with the seed and the engine name (checked against `UPSIDE_HOME`) as arguments; `bench.sbatch`
  (`5ab60bd9`) sources `env_upside3.sh` or `env_shared.sh` by `ENGINE`; Upside 2 seeds each system seed + index, so 101
  and 202 share no stream with seed 1. `ff3_benchmark`'s scripts are untouched (ff3.1's arms read them every chunk).
  When all 8 are COMPLETE: score them and set the spread over seeds 1, 101, 202 within each engine against the Upside
  3 - Upside 2 difference.
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
    error line); its readout (C-N, DSSP) is owed. Block 3 of 4 **49584589** R from 10:24 10-27.
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
  **49573354** R from 18:03 10-27 (log `REMD/glpG-RKRK-79HIS.remd_block1.log`, no error line at step 2,222).
- **glpG with ff3.0, the matched TM4 test** (user 10-26: "ff3.0 surely fix the TM4 issue even before the bug fix in
  upside2 hybrid, so we want to run it to completely be sure", then "do submit ff3.0"; findings 3.15). A235: the
  preparation compares only the three `sidechain.h5` datasets the side-chain/lipid table is built from, so ff_2.1,
  ff_3.0 and ff_3.1 pass; in the Upside 3 tree 14:10 10-26, `py/martini_prepare_system.py` `802e1a90` (old
  `.bak_pre_A235_20261026`, `5e6b4627`). The tree's `parameters/ff_3.0` is the released `ff_3.0_gdepth` (all six files
  md5-equal), which the placement-fixed Upside 2 chains ran. `upside3_validation/glpg_ff30/run.py` (`db1cf3b4`) is
  `glpg_ff31/run.py` with `force_field = "ff_3.0"`, same seeds; outputs and prep log in
  `glpg_ff30/outputs/martini_glpG-RKRK-79HIS_hybrid/`. Prep 49589828 COMPLETED 14:38 10-27: **TM4 held** (last frame
  0.88, last tenth 0.82; TM1 0.95): the ff2.1 run's loss is ff2.1's own with correct placement, not Upside 3's. REMD
  block 1 **49589829** R from 18:04 10-27 (log `REMD/glpG-RKRK-79HIS.remd_block1.log`, no error line at step 2,222), to
  be read at T 0.70 against Upside 2's placement-fixed ff3.0 chains (DSSP 0.96, block 1).
  - **Block 2 49634605 PD afterok 49589829** (log `REMD/glpG-RKRK-79HIS.remd_block2.log`): `n_blocks` 1 -> 3 by hand
    from the Mac Studio 20:32 10-27 (old `run.py.bak_pre_nblocks3_20261027`). Block 1 loaded `n_blocks = 1` and does
    not queue block 2, so 49634605 was submitted with `--dependency=afterok:49589829` and the sbatch line of `run.py`'s
    `submit()`; it reads `block_count` 1 and queues block 3, so blocks 2-3 pool for HDX as Upside 2's ff3.0 chains did
    (findings 5.6). **Do not run `python3 run.py` there before block 1 ends:** `block_count` does not exist yet, so it
    would submit a second block 1 that copies the seed over every replica file. If block 1 fails, 49634605 stays
    pending (DependencyNeverSatisfied): cancel it.
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
  Next: the npz to the MacBook Pro, `e1_readout.sh` on them, the setting by the task file's rule, then E1's calibration
  and support runs at 1,024 particles and two seeds.
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
| **61701305** (midway3 gpu V100) | Boltz-2 Phase 2 E1 tuning: 12 steered runs of 256 particles on Ub psWT V70 and residue 61 (line 4) | PD from 11:11 10-28 | readout on the MacBook Pro, the setting by `e1_tune_ubqps.txt`'s rule, then E1 proper |
| **49584589** (midway2) | Upside 3 glpG hybrid with ff2.1, REMD block 3 of 4 (`upside3_validation/glpg`; line 3) | R from 10:24 10-27 | block 2's readout (C-N, DSSP) owed; HDX waits on the user's choice of force field |
| **49573354** | Upside 3 glpG with ff3.1, REMD block 1 (`glpg_ff31`); TM4 lost in prep, `n_blocks` stays 1 | R from 18:03 10-27 | readout when it ends: TM4, KE, C-N at T 0.70 |
| **49589829**, **49634605** | Upside 3 glpG with ff3.0, REMD blocks 1 and 2 (`glpg_ff30`); TM4 held in prep, `n_blocks` 3 | 49589829 R from 18:04 10-27; 49634605 PD afterok | block 1's TM4 at T 0.70 against Upside 2's 0.96. Never run its `run.py` before block 1 ends; if block 1 fails, cancel 49634605 |
| cspA 49629882, gpW 49630880, hyp 49630574, NuG2 49630879, proteinG 49629732, proteinL 49630878, ubiquitin 49631645 (`u3_<prot>_denovo`) | Upside 3 Peng arms with ff2.1, the 7 de novo left (`upside3_validation/peng`; line 3) | R, 10.5 h chunks resubmitting themselves | rescore against Upside 2 as each finishes |
| **49634596-49634603** (`u<3\|2>s<seed>_<prot>_denovo`) | Peng seed replicates: top7 and alpha3d de novo, seeds 101 and 202, both engines (`upside3_validation/peng_seeds`; line 3) | R from 20:32 10-27 | all 8 COMPLETE: score with the seed-1 arms; seed spread against the engine gap |
| **49635669 (ff3.2), 49635327 (ff3.1) PD** (midway2) | **ff3.2 and ff3.1 on the fixed native target** (plan.md Phase 33): `ff32_nativefix` (every term; priority over ff3.1, user 10-28) and `ff31_nativefix` (depths only), from the released ff3.1 | first links 49635326/49635325 FAILED 10:41 10-28 in a Slurm controller outage (ff3.2 after step 10, ff3.1 after step 15); successor links held 10:50 and released 10:58 10-28 (user: "proceed with your best decision"), ff3.1's at nice 10; PD (Resources), partition full | when each starts, give the insurance link it submits the same nice (ff3.1's 10, ff3.2's 0); then each round end (line 2) |
| **49634973**, **49635107** R (successors **49636476**, **49636479** PD afterany; all nice 10) (`condiv-train`) | The ff3.1 training debug, `ff311_glypair` and `ff311_glypair_mb2` (plan.md Phase 32; line 2), on the old native target | first links 49634971/49634972 FAILED 10:41 10-28 in the same outage; second links R from 10:58:56 10-28 from their last checkpoints; ranked behind ff3.2 | each round end: glycine's pair, its pull and -ln(P_aL/P_aR) against ff3.0's |
| **49630881** (`b_ff_2.1_ubiquitin_denovo`), collector **49633770** (`ff21_collect`) | Released ff2.1's last Peng arm (31 of 32 scored) and its collector (plan.md Phase 27; line 1) | arm R; collector PD, 00:47 10-28 | `peng_ff21.txt` (DONE); then Phase 27 step 4 |
| 29 `b_ff_3.1_<prot>_<kind>` chunks (ids in line 2), collector **49634812** (`ff31_collect`) | ff3.1's Peng arms (3 of 32 scored) and their collector (plan.md Phase 30; line 2) | arms R; collector PD, 04:07 10-28 | `peng_ff31.txt` when DONE |

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
  validated); ff31_depth_data (released ff_3.1_depth_data 10-13, retired 10-26; 32 Peng arms complete 10-20).
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
  `ff3_benchmark/runs/<prot>_<mode>_ff_3.0/` and `popepopg_REMD_mdw2/<V>/`; superseded benchmark
  `runs_superseded/ff_3.0_20260930-022510`, pre-release seeds `seeds/*.bak_pre_ff_3.0_20260930-022510`. **Each base
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
in `samples.xtc`): data `$E/runs/bioemu2_hdx_<ehee|ubqps>_1e5_v100_<id>/`. Local copies (md5-equal,
`readout/samples/jobs/`), the pooled sets `readout/samples/{ehee,ubqps}_bioemu2_n100000.npz` and the readout in
`scratchpad/bioemu2_hdx/readout/`.

### Disk: from midway2, `df` on the subdirectory is the only number you get for `/project`

```
df -h /project/trsosnic     -> the FILESET (size, used, free), correct
df -ih /project/trsosnic    -> fileset inodes
df -h /project              -> 6.3P total, the whole device, useless
rcchelp quota               -> home, scratch, project2 group only; no /project or /beagle3 row
mmlsquota                   -> "File system project is not known", the GPFS client is not here
```

`df` on the subdirectory reports the fileset because GPFS `--filesetdf` is on; `df` on the mount point reports the
6.3 PB device. That distinction is the whole trap.

Now: `/project/trsosnic` 527 G free (22:25 10-27; tell the user before it falls under 200 G); `/beagle3/trsosnic`
125,183 inodes free (the watch reports IFree under 30000). `/project2/trsosnic` (group) was at 1.45 T of a 1.49 T
soft quota (97%) on 09-23; nothing in these campaigns writes there. Midway3 home: 30 G quota (§8).

**`/beagle3` inode exhaustion, 12:20-14:30 10-24** (findings 10.39): the first submission of Boltz's 10^5 tests filled
the inode quota. Jobs writing there failed: Upside 3's glpG REMD could not write `block_count` (`OSError: [Errno 122]
Disk quota exceeded`), and four Peng chunks died in 9 s with no log (exit 0:53: Slurm could not create the log); both
were resubmitted (§1 line 3).

**Reclaimed at the user's word** (manifests, scripts and logs kept):
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
`gly_peptides` (§1b), `glpG_DDM_micelle_REMD` (73 G, §1b), `implicit_79HIS_*` (70 G), `popepopg_glpG` (9 G). The 12
ff3.0 HDX work dirs (417 G) are kept at the user's word (§1 line 2). What grows: each panel ~11 G (mdCATH ~47 G), each
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
* **Safe to delete:** `/beagle3/trsosnic/yinhan/ff3_benchmark/sigmoid_discarded/` (4 G) and `armM_discarded/`.
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
* **A training chain insures itself** (the design that worked unattended on ff30_basin): every link queues its
  `afterany` successor before it trains; `--no-requeue`; the successor resumes from the newest checkpoint after
  deleting the half-written step; a worker whose `srun` never started is relaunched (twice at most); three links from
  the same step stop the chain. The gate job itself has no successor.
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
