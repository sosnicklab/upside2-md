# Remote jobs on midway2/midway3: status and handbook

**Current state (2026-10-21 13:50).** No watch runs (the owner of record is the Mac Studio, idle since its 21:50 10-18
pass). The running jobs look after themselves: each line of work below has a collector job or a dependency chain that
writes its own result file, so nothing needs polling, and `submit_new.sh` is idempotent if a watch also runs it.
`WATCH_STATUS.md` on the cluster is newer than this file whenever a pass ran after it. Job ids, logs and data paths
are in §1.

The user's four lines of work (10-21), and what each has on the clusters:
1. **Deterministic vs probabilistic side chain** (plan.md Phases 20-27; findings 1.21, 1.37). Both trainings are
   finished and converged; the deterministic side chain folds less on every held-out test, and both retrained
   force fields fold less than released ff2.1. **Running:** released ff2.1 on the same three tests, to decide whether
   ff2.1 itself stays: the 44 domains are done, mdCATH reads out ~20:00 10-21, Peng's 30 arms in about a week. The
   email to Tobin and John waits on these numbers.
2. **Ramachandran map training** (plan.md Phases 14-17, 21, 28-29; findings 1.34, 1.37). **Running:** ff3.1.1 and
   ff3.2.1, ff2.1 with only its basin depths trained under the pooled prior (from NDRD's and ff3.0's map); the first
   test of the accuracy fix is their epoch-0 panels (~02:00 10-22); step 76 ~10-23. **Finished, awaiting review:**
   the ff3.1 and ff3.2 reruns on the fixed engine, whose maps barely moved and which, as every run trained from
   ff2.1 with all terms free, fold 0.15-0.24 less than ff2.1. **Also running in this line:** the 8 glpG REMD chains
   on the released ff3.0 and ff3.1 (block 3 of 5), the TM4 and HDX validation.
3. **Upside 3 design, debug stage** (upside3.md; plan.md Phase 24 and its "Resume here"). Local work, no cluster jobs.
4. **Boltz-2 HDX sampler** (diffusion_dynamics.md). The tag by step-scale test on midway3 finished 09:02 10-21; its
   readout into Phase 0 is owed. Nothing running.

Open and not in the user's list:
- The 8 glpG chains run force fields that carry the retraining loss; whether to stop them when block 3 closes
  (~18:00 10-21 to ~01:30 10-22) is the user's. Their TM4 check at block 2's end and the HDX on block 2 are owed.
- Large-protein de novo folding (plan.md Phase 22) is done; its proposal (compact misfolds as negative examples,
  training past 150 residues, replica exchange for the large-protein benchmarks) awaits the user.
- The BP validation deck for Tobin (§0c) is not sent.
- Awaiting the user: the large directories on `/project/trsosnic` (§1 "Disk") and the deferred HDX post-processing
  fixes (findings 5.5).

The ConDiv runs, all on midway2 broadwl at dt 0.009 (`$P/training/<run>`; "trained" is the port's rates, hb 0.01, rot
0.0125):

| run (tags) | what it is | state |
|---|---|---|
| ff21_nobp_dt009 | ff2.1 retrained with one side-chain state per residue type (the deterministic side chain) | converged at step 76 10-19; step 37 released as ff_2.1_nobp |
| ff21_ctrl_dt009 (`c9_EE`) | ff2.1 retrained the same way with BP (the probabilistic control) | converged at step 76 10-20; step 37 released as ff_2.1_ctrl38 |
| ff30_gdepth_dt009 (`d9_EE`) | ff3.0: NDRD plus one trained glycine depth pair, every term trained | released as ff_3.0_gdepth 10-09 |
| ff31_depth_data (`fd_EE`) | ff3.1: 1,442 per-map basin depths, prior 0.2, every term trained | released as ff_3.1_depth_data 10-13 |
| ff31_depth_all (`fa_EE`) | ff3.1 with all 5,040 depths, prior 0.3 | finished 10-13, not validated |
| ff32_depth_data (`f32_EE`) | ff3.1's depths started from ff3.0 | stopped by the user in step 69 10-18 |
| ff31_depth_data_angfix (`a31_EE`), ff32_depth_data_angfix (`a32_EE`) | the ff3.1 and ff3.2 runs on the fixed engine | finished step 76 10-21, awaiting the user's review |
| ff311_depth_pooled, ff321_depth_pooled | ff3.1.1 / ff3.2.1 with every term trained | cancelled 13:38 10-21, 3 min in (user's revision); kept as `<run>.cancelled_20261021` |
| **ff311_depths_only (`p311_EE`, ff3.1.1)**, **ff321_depths_only (`p321_EE`, ff3.2.1)** | ff2.1 with only its depths trained, pooled prior, from NDRD's and ff3.0's map | **running from 13:42 10-21** |

Dropped by the user and deleted 10-11 and 10-12: the BioEmu-map line, ff30_gdepth_si, the H-bond-frozen arms,
polygly and the local TM4 glycine-flip test (§1 "Finished runs"). midway2 has no GPU that can run Boltz-2.

CPU work runs on midway2 broadwl unless the user names midway3 caslake. midway3's login node may be
used to move or read files on `/project` and `/beagle3`. Since 10-07 10:48 midway2-login2's
`/project` hangs; route `/project` work through login1 (§8).

Written so a fresh session can pick up cold. §1 is live job state; §0, §2 and §5-§8 are the
handbook.

---

## Resume here, from any computer (2026-10-21)

### "continue jobs": what a new session does, in order

Claude's memory does not travel between computers; everything needed is in this file, plan.md and
findings.md (pulled from git) and on the cluster.
1. **Read** this file's header table and §1 (live state; each watch pass leaves its readouts
   there, and "What happens next" lists the coming events), then:
   - plan.md Phases 18 (the glpG placement fix), 21 (the angular-gradient fix), 27-29 (released ff2.1, the
     retraining loss, ff3.1.1 and ff3.2.1) and Known Errors;
   - findings.md 1.37 (the retraining loss), 1.34 (the depth maps), 3.12, 4.17 (the side-chain placement defect), 1.21 (the
     deterministic side chain against its control: panels, audit, H-bond margin), 1.35 (TM4), 5.5 (TM4 HDX),
     and the rules in §10: 10.16-10.19 (no test jobs, same data and workflow, re-read the queue, TM4 only
     in the hybrid), 10.20 and 10.27 (DSSP beside any helix readout), 10.28 (name runs by what they are).
2. **Connect** to midway2 ("Connect" below).
   - Check the socket first. If `scratchpad/mdw2_master.exp` and `~/.bin/ssh_mdw3` exist on this
     computer, run the expect script once (one Duo push); otherwise ask the user to open the master.
   - While login2's `/project` hangs, run `/project` work through midway2-login1 (§8). Test with one
     `ls` first.
3. **Take stock, read only.** Watch steps 1-2 and 2c ("The watch"): `squeue`, the log scan, and every
   job that left the queue since §1's snapshot. Tell the user what started, finished or failed since
   then.
4. **Owed checks** (as of 13:50 10-21; `WATCH_STATUS.md` lists them too, and is newer if the two differ). The
   collectors do the routine work for released ff2.1's tests and for ff3.1.1 / ff3.2.1: read their result files (§1).
   Owed by hand:
   - **TM4 at the block-2 end of the 8 corrected glpG chains** (findings 1.35 is the method and holds block 1's
     numbers, 10-18; the defective chains' numbers are no baseline): on midway2 login1 in
     `checks/tm4_ff31_20261014/`, `source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh`, `export
     HDF5_USE_FILE_LOCKING=FALSE`, `ulimit -v 6000000`, `nice python3 tm4_ff31.py > tm4_ff31_<date>.txt` (40 s)
     and `nice python3 ~/watch_probe/gly_flip_rungs.py > gly_flip_rungs_<date>.txt` (150 s); both list
     ff_3.1_depth_data against ff_3.0_gdepth in the four variants. DSSP alpha beside them (findings 10.20,
     10.27): on midway2 `nice python3 - FIRST LAST < extract_bb_dssp_block1.py > bb_dssp_<block>.npz` (79HIS's T
     0.80 rung, both force fields; FIRST LAST are the block's `output_previous_` groups, block 1 was 1-17), then
     `python tm4_dssp_block1.py bb_dssp_<block>.npz` on a computer with mdtraj (the repo `.venv` has it), and
     upload the output beside the others. Block 2 closed on all eight 07:55-14:56 10-20; block 3 ends ~36 h after.
   - **TM4 HDX on the corrected chains' block 2** (the user's 10-14 requests: ff_3.0_gdepth and
     ff_3.1_depth_data, 79ALA then 79HIS): §1 "TM4 HDX on ff_3.1_depth_data: the recipe" with each chain's
     own group counts. Block 2 has closed.
   - **The Boltz-2 test's readout** into `diffusion_dynamics.md` Phase 0 (§1 line 4).
5. **Find the owner** ("Handoff" below): `cat` the cluster's
   `checks/r4_epochs/tm4_local/WATCH_STATUS.md`. **Planned move (user, 21:55 10-18: "expect me to run the
   monitoring on another computer tomorrow"):** the user's opening the session on the next computer on 10-19 is
   the "take over" for it. Session `upside2-md-f1` (cron `82dc7422`) is closed, so the Mac Studio runs no watch;
   the passes since (10:24 10-19, 13:02 10-19, 22:30 10-20, 11:10-13:50 10-21) started none; the last, at the
   user's word, submitted released ff2.1's validation and ff3.1.1 / ff3.2.1, each with its own collector. If `WATCH_STATUS.md` names a
   pass under 2 h old, ask the user to confirm no watch runs on the Mac Studio before starting one. Two watches
   would both submit panels. With "Owner: none", or the user's confirmation, go to step 6.
6. **Take the jobs:** run one watch pass now (steps 1-5 of "Round-4 watch"), start the watch
   ("Starting the watch"), and rewrite `WATCH_STATUS.md` with this computer as owner.
7. **Open questions for the user.** Nothing is blocking. Awaiting the user's word, not to be done
   without it:
   - Whether to stop the 8 glpG chains when block 3 closes (they run ff3.0 and ff3.1, which carry the retraining
     loss; findings 1.37).
   - The large-protein de novo proposal (plan.md Phase 22 review).
   - The large directories still on `/project/trsosnic` (§1 "Disk": popepopg_REMD variant dirs 355 G with
     `hdx/` kept, glpG_DDM_micelle_REMD 73 G, implicit_79HIS_* 70 G, popepopg_glpG 9 G).
   - The deferred HDX post-processing fixes (plan.md Known Errors; findings 5.5).

### Handoff: two computers, one owner (2026-10-14)

**Owner now: the Mac Studio** (10-17 00:13, at the user's "check remote job status, and start
monitoring"), **until the user moves the watch to another computer on 10-19** (21:55 10-18; "Resume here"
step 5). Its watch cron `82dc7422` ended with session `upside2-md-f1` after the 21:50 10-18 pass. It held the jobs before, from 10-08 09:15 until the user closed its
session 10-14 ~16:00 (cron `c618e779` deleted then); none owned them between. The rules below hold
with whichever computer owns them: another computer may open a session
at any time; it must know what the owner is doing and not do it a second time, since two watches
would both submit panels. **One owner at a time.** (The local TM4 glycine-flip driver
the owner also ran was retired 10-11; its record is §1 "Local glpG TM4 test".)

**Owner's session restart.** The watch cron is session-only and ends with its session; the computer
stays the owner, and the midway2 socket master is its own process, which survives unless the
computer restarts. The new session on that computer:
1. Reads §1 and `WATCH_STATUS.md`, and checks the socket (the expect script once if it is dead).
2. Runs one watch pass at once (watch steps 1-5) for the steps since the last pass.
3. Starts the watch again: CronCreate hourly at `23 * * * *` with the prompt in "The watch" below.
4. Rewrites `WATCH_STATUS.md` (last pass; any restart note removed).

**Live status on the cluster:** `/project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md`.
The owner rewrites it at every watch pass (watch step 5): owner, time of the last pass, the newest
step of each run, and what is owed. This file (`remote_jobs.md`) reaches another computer only
through the user's commits, so it can be hours old there; `WATCH_STATUS.md` is the current word.

**What a session does on "continue jobs"** (after "Resume here" steps 1-4):
1. `cat` `WATCH_STATUS.md` over the midway2 socket.
2. **"Owner: none": take the jobs** ("Resume here" step 6): one pass, start the watch, write this
   computer into `WATCH_STATUS.md` as owner. Its shell and paths may differ: `~/watch_probe/` and
   `checks/` hold every tool the watch uses, on the cluster.
3. **An owner on another computer whose last pass is less than 2 h old: stand by.** Do not start a
   watch or submit anything; reporting status read only and answering the user are fine. Take over
   only on the user's word ("take over"). If its last pass is more than 2 h old, tell the user and
   ask whether that session has stopped; do not assume it. Taking over: if the old owner's session
   is still open, the user has it run its stop steps first; then open the midway2 socket ("Connect"
   below), start the watch ("Starting the watch") and rewrite `WATCH_STATUS.md` with the new owner.
   The probes `~/watch_probe/stock.sh` and `~/watch_probe/steps.sh <run>:<step> ...` on midway2 do the
   training runs' watch steps 2-3, `~/watch_probe/logs.sh` and `final_frames.py` steps 2 and 2c.

**The owner's stop steps** (when the user says "hand off", "stop" or closes the session): finish the
pass in hand, update §1 and `WATCH_STATUS.md` (owner: none), then CronDelete the watch cron.

### Where things are

**What travels with git and what does not.** Tracked: `plan.md`, this file, `findings.md`,
`progress.md`, `training/`, `src/`, `py/`, `up.md` and `parameters/common/rama31.dat` (the
BioEmu-fitted library). Commit here, pull there. Not in git (`.gitignore` has `*scratchpad*`): the
`scratchpad/` tools, including `mdw2_master.exp`, `mdw3_master.exp` and the local TM4 test (the
MacBook Pro's scratchpad lost its run data 10-15 at the user's "cleanup", by the findings 10.25 rule:
scripts, inputs and small results kept, manifest `scratchpad/CLEANUP_20261015_manifest.txt`). Their
authoritative copies are on the cluster:
- `ff3_selection/` holds the panel tools (`panel.py`, `panel.sbatch`, `submit_new.sh`,
  `slurm.args`).
- Under `/project/trsosnic/yinhan/checks/`:
  - `r4_epochs/<tag>`: extracted checkpoints.
  - `r4_epochs/tm4_local`: the TM4 test (README, `scripts/`, `seed/`, `runs/` pre-fix, `runs_cov/`).
  - `cov_fix_20261007`: the tested patch-script fix.
  - `fz_init_20261007` and `dt009_init_20261006`: the initial force fields, shown byte-identical.
  - `r4val_20261005`: the automatic-validation tests.
  - `broken_replica_20261006`: the comparisons behind findings 1.21.
  - `gly_bioemu_map`: the BioEmu library, its fit and the push probe, with a README.
  - `hbg_deploy_20261002`: the worked binary deploy.
  - `ff31_deploy_20261009`: ff3.1's trainer files as deployed, the patched `check_step.py` (a branch
    for per-map offsets; on a round-4 step its output is identical to the old one's) and the test.
    The older `training/ff31` and `ff31-gly` are September's retired FF1-form glycine retrains.
  - `tm1_hbmem_20261009`: the MacBook Pro's 10-09 glpG TM1 work (findings 1.27 and 4.6-4.9, plan.md
    Phase 13): scripts, tables, the membrane H-bond term runs' logs and positions, the replay, and
    a README of what stayed on the MacBook Pro (those `.up` files were deleted in the 10-15 scratchpad
    cleanup; scripts and tables kept). Its `mac_studio_20261010/` holds the Mac Studio's
    follow-up (findings 4.10-4.14): the fixed-table inputs, both fixed-table sets, the kick-probe
    analysis and `README_mac_studio.txt`.
  - `sc_table_20261010`: the converged side-chain/lipid table's deployment (findings 4.15): the
    staged martini.h5, builder and `patch_glpg.py`, `deploy.sh`/`.log`, and `repatch.sh`/`.log` with
    `check_repatch.py` for the 8 re-patched glpG seeds. Every replaced file is kept beside itself as
    `.bak_pre_sc600x12_20261010`. `$P/py/martini_build_tables.py` keeps the old sampling (campaign
    rule): never rebuild martini.h5 in `$P`.

Claude's memory is per computer; the rules that matter are in this file and findings.md §10.

**The rotamer-BP validation figure** (for Tobin, unsent) is in
`/beagle3/trsosnic/yinhan/bp_validation/static/`, with its data and `plot_bp.py`. The 2-slide deck
`bp_validation_slides.pptx` was built by `scratchpad/bp_validation/make_slides.py` on the computer
that made it (§0c).

**The cluster's `training/` is not the repo's layout, and must stay so while round 4 trains.**
`$P/training` keeps the pre-10-02 layout:
- `ConDiv.py`, `train_chain.sbatch`, `check_step.py`, `extract_ff.py`, `convergence_gate.py`,
  `gate_or_continue.sh`, `validate_ff.sh` and `patch_glpg.py`;
- each run's own trainer copy in `<run>/trainer/` where it differs.

The runs' `after_training.sbatch` calls `gate_or_continue.sh`, which runs `validate_ff.sh`. Do not
sync the repo's `training/` over it. The repo's `ConDiv.py` and `train_chain.sbatch` carry the same
dt 0.009 and 54-step links.

**Connect.** From a new computer, open the master yourself, with password and Duo:
```bash
ssh -M -S ~/.ssh/cm-mdw2.sock -o ControlPersist=8h -o ServerAliveInterval=30 yinhanw@midway2.rcc.uchicago.edu
```
Later commands reuse it: `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu
'<cmd>'`. To let Claude open it unattended, copy `scratchpad/mdw2_master.exp` over (§0) and put the
password line `set password "..."` in `~/.bin/ssh_mdw3`. While login2's `/project` hangs, run
`/project` work as `... 'ssh -o BatchMode=yes midway2-login1 bash -s' < script.sh` (§8). **Never run
anything heavy on a login node** (findings 10.6).

### The watch

**One watch at a time.** A watch is a session cron: it dies with its Claude session and expires
after 7 days. **None runs now** (13:50 10-21): the Mac Studio's, cron `82dc7422`, ended with its session after the
21:50 10-18 pass, and no computer has started one since. The jobs of 10-21 do not need it: their collectors run the
per-step checks and submit their panels (§1). A watch, if started, adds the glpG chains' block ends, failures and
disk; `submit_new.sh` is idempotent, so it and a collector cannot both submit one panel. Each pass records its
readouts in §1, so the newest §1 tells a new watch where the last one stopped.

**Starting the watch on a computer** (set up as in the Handoff): CronCreate hourly at an off-minute
(e.g. `23 * * * *`) with the prompt below.

**Round-4 watch prompt** (rewritten 10-21 for the jobs now running; the rules are unchanged):
"Round-4 remote job watch, the only monitor (remote_jobs.md "The watch"). Two trainings run, ff3.1.1 and ff3.2.1
(ff311_depths_only, p311_EE; ff321_depths_only, p321_EE): ff2.1 with only its basin depths trained, 0 to 76, no gate,
stopping for the user (plan.md Phase 29). Their collector `pooled_collect` does their per-step checks and panels and
writes checks/ff311_20261021/status.txt: read it, and report the epoch-0 panels against ff21_released. Released
ff2.1's validation (plan.md Phase 27) runs its own collector `ff21_collect`; read
checks/ff21_released_20261021/ for new files. The 8 corrected glpG REMD chains on ff3.0 (ff_3.0_gdepth) and ff3.1
(ff_3.1_depth_data) run five blocks (plan.md Phase 18): TM4 at each block end, HDX on block 2 (§1 recipe).
Job ids are in §1. Follow "Round-4 watch" steps 1-5 exactly; read §1 first. In reports, name runs by
what they are with the tag in parentheses (findings 10.28). Obey:
- CLAUDE.md: read-only git, no guards, no physics changes, nothing heavy on a login node;
- findings 10.16 (no test jobs), 10.17 (same data, same workflow) and 10.18 (re-read the queue
  before any decision on job state).
Never cancel, resubmit or change a job without the user's approval.
Notify only for: failures, a panel result, a gate, a new margin hold in a trained-H-bond run
(below +0.10), a steric flag, a disk warning, or a dead socket.
Otherwise stay quiet."

**Round-4 watch.** `<run>` is any of the running or queued run directories under `$P/training` (§1).
`submit_new.sh` maps them to their tags.
1. **Socket.** `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu`.
   If it fails, run `scratchpad/mdw2_master.exp` once, unless §1 records that the previous pass
   already tried. If it is still down, record that in §1, notify once and stop the pass. Never loop
   on connection attempts.
2. **Queue and logs.** Run `squeue -u yinhanw` and `bash
   /project/trsosnic/yinhan/ff3_selection/submit_new.sh`, which submits each new epoch end's panel
   without asking (user, 10-08 11:05). Check the newest link log of each run for
   `WORKER_FAIL`, `Traceback`, `STOPPED`, `never started`. The benchmark arms' and glpG chains' logs:
   `ssh ... 'bash ~/watch_probe/logs.sh' | grep -v -E "done=[0-9]+ / target" | grep -v -E "^== .*errs
   0$|^\[remd\] chunk [0-9]+ steps"`; no output (exit 1) is clean. A chain starting a block prints
   its `block k/5` and `calibration` lines, which are expected.
2c. **A job that left the queue** (compare the ids with §1): `sacct -n -X -j <id> -o
   JobID,JobName%45,State,End,Elapsed,ExitCode` and its log tail.
   - Benchmark arm (log `/beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<kind>_<ff>_<id>.out`):
     each chunk stops at its 10.5 h wall and prints `[bench] resubmitted: <id>`; check the new job
     carries the law (`scontrol show job <id> | grep -o "ExcNodeList=[^ ]*"`). On `[bench] COMPLETE`,
     read its final frames: `source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh; export
     HDF5_USE_FILE_LOCKING=FALSE; nice python3 ~/watch_probe/final_frames.py
     /beagle3/trsosnic/yinhan/ff3_benchmark/runs/<prot>_<kind>_<ff>` (KE/1.5kT near 1.00, finite).
   - glpG chain (log `popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<id>.out`): a whole block ends with
     `[remd] block complete` after `Submitted batch job <next>`; read the last `avg_kinetic_energy/1.5kT`
     line and the next job's `block k/5` line.
   - Record the new ids in §1 and `WATCH_STATUS.md`.
2b. **Disk.** `df -h /project/trsosnic`; record the free space in §1 "Disk" when it moves by more than
   100 G, and tell the user before it falls under 200 G (10-11: it reached 88 G unnoticed).
3. **Every step finished since the last pass.**
   - Run `cd $P/training && source <run>/env.sh && python3 check_step.py <run>
     <epoch_xx_minibatch_yy>` under `ulimit -v 4000000`.
   - KE scan of every protein: `grep -m1 ^avg_kinetic_energy/1.5kT run_output/<step>/*.run.*.up.output`.
     Flag a maximum above 1.2 (findings 1.21).
   - Record the shared margin E_other - E_alpha (start +0.192; ff3.2's +0.082), dhb, sheet mean, rot rms, step time,
     and helical / left glycine alpha_L free/restrained. For the gdepth runs also record dL - dR
     (start +0.554) and `run_output/rama_rounds.txt`.
   - **A trained-H-bond run whose margin falls below +0.10:** notify the user once and hold. The run
     keeps training; nothing is cancelled without the user.
   - **The ff3.1 and ff3.2 runs: the steric check** (user 10-10, after Jumper; plan.md Phase 14 step 5,
     findings 1.33). `cd $P/training && source <run>/env.sh && nice python3
     /project/trsosnic/yinhan/checks/ff31_steric/steric_check.py <run> <step> >
     /project/trsosnic/yinhan/checks/ff31_steric/<run>_<step>.txt` under `ulimit -v 4000000`
     (without `<step>` it reads the newest finished step)
     (3 min, 360 MB; repo `training/ff31/steric_check.py`, md5 `e35e3d98...`). Record the
     side-chain clash probability now / start, the largest map and residue ratios, and the
     non-glycine free beta' + pPII'. A `STERIC FLAG` line (a coil map or training residue whose
     clash probability more than doubled against NDRD, the training residues' mean side-chain clash
     probability up by more than 10%, or a step's free non-glycine beta' + pPII' above the run's
     steps 0-18): notify the user once, with the table. The run keeps training; nothing changes
     without the user.
**The local glpG TM4 glycine-flip (stability) test is retired** (user, 10-11 23:10: "We have a working ff and we are trying
to expand it so that metric no longer matters and no longer needs to be tested locally"). ff2.1
holds TM4 on the fixed side-chain/lipid table (findings 4.15), so new epoch ends are not extracted,
patched or queued for it. The tools stay in `scratchpad/ff3_local_test/` (`run_glpg.sh`,
`tm4_queue.sh`, `tm4_compare.py`; method in findings 1.25 and 4.15) and every set's logs and
tables are in `tm4_local/runs_cov/` (the trajectories were deleted 10-18). **TM4 HDX is not retired** (user 23:15: "Whether the delta g of tm4 should go to
infinity or not is still not resolved"): ff_3.0_gdepth's block 2 is done and ff_3.1_depth_data's is owed
(§1 "TM4 HDX on ff_3.1_depth_data: the recipe"; the ff2.1_sc chain was cancelled by the user 10-14 01:11).
4. **A panel or a gate.**
   - A panel whose `ff3_selection/runs/<tag>/` holds 176 npz: run `panel.py select aa domains runs
     ff21_released ff21_awh gdepth_start <tags>` with the run's env (70 s, 300 MB; bio_start was
     deleted with the BioEmu line 10-11 23:45) and
     notify the table. Put each rerun's panel beside its source's: `a31_EE` with `fd_EE`, `a32_EE` with
     `f32_EE`, `c9_EE` with `c9_00`/`c9_01`. Those old panels' `.up` files were deleted 10-18, so they compare
     on folded fraction and Q only, never by `panel_compare.py`: paired by
     `checks/converge_20261020/paired_npz.py DOMAIN_DIR AA_DIR RUNS_ROOT TAG_A TAG_B [CSV]` (panel_compare's folded
     and Q from the npz, reproduced exactly; from `ff3_selection/` with `training/ff30_gly/env.sh`, a minute).
   - **Gates of the two ff2.1 runs** (`ff21_nobp_dt009`, `ff21_ctrl_dt009`). At the target,
     `train_chain.sbatch` cancels its queued successor link and submits the run's `after_training.sbatch`
     (job `nobp-gate` / `ctrl-gate`, 1 h; logs `ff21_nobp_dt009/setup/gate_<id>.out`,
     `ff21_ctrl_dt009/gate_<id>.out`), which runs `ff21_nobp_dt009/setup/after_gate.sh`: the report
     `<run>/gate_step<N>.txt` (`convergence_gate.py` over the last epoch). Converged: the epoch-end files
     go to `ff3_selection/params/<tag>` and the mdCATH panel is submitted, tag `mdc_nobp_EE` / `mdc_ctrl_EE`
     (`EE` the epoch, 03 at step 76; job `ff3sel_mdc_<tag>`, log `ff3_selection/logs/ff3sel_mdc_<tag>_<id>.out`,
     runs `ff3_selection/runs_mdcath/<tag>/`, 356 runs, ~1.6 h for nobp, ~6 h for ctrl). Not converged: one
     more epoch (`train_chain.sbatch <run> <N+19>`), gate again, up to step 114; then "stopped for review".
     Notify the gate either way, with the failing groups.
   - **An mdCATH panel with 356 npz:** once both `mdc_nobp_EE` and `mdc_ctrl_EE` are done, `cd ff3_selection &&
     sbatch $(cat slurm.args) panel_compare.sbatch domains_mdcath aa_mdcath runs_mdcath mdc_nobp_EE mdc_ctrl_EE
     mdcath_candidates/candidates_mdcath.csv` (read-only, as 49412845; output `panel_compare_<id>.out`:
     folded, mean Q, CA-RMSD, paired, by length bin and fold class) and notify it beside the step-37 result
     (findings 1.21: 0.264 against 0.365). Until the control's is done, the first one alone: on login1 with
     the run's env, `python3 panel.py select aa_mdcath domains_mdcath runs_mdcath mdc_ctrl38 mdc_nobp <new tag>`
     (from the npz; the step-37 tags have no `.up` since 10-18, so `panel_compare.py` cannot read them).
   - A gate: read `gate_step<N>.txt` and the gate log. Converged means `validate_ff.sh` ran inside
     it, so confirm 32 Peng arms and 4 glpG chains for the candidate were queued, and read their
     first logs.
5. **Update §1 briefly, and rewrite the cluster's `WATCH_STATUS.md`** ("Handoff"): owner and
   computer, time of this pass, the newest step of each run with its margin, and anything owed or
   awaiting the user.

**Validation of the round-4 candidates is automatic** (§1 "Automatic validation"). Each converged
gate runs `validate_ff.sh`, whose glpG seeds now carry the candidate's coverage tables. Choosing
which candidate becomes ff_3.0 stays the user's.

## 0. Connect first (needs a Duo push on the user's phone)

Password + Duo is the only method. Key-based ssh is refused on both clusters: tested 2026-09-17 with
a correctly installed `authorized_keys` (600/700) and again with `PubkeyAcceptedAlgorithms=+ssh-rsa`;
the server advertises `publickey` but does not honour it for this account. Do not try again. An inert
key entry was left in midway2's `~/.ssh/authorized_keys`; it is harmless. The ControlMaster socket
expires roughly hourly, so expect to redo this most sessions.

**midway2** (training, selection panels, benchmark arms):
```bash
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp    # USER MUST APPROVE DUO
```
If midway2 blocks or throttles this IP, tunnel through midway3 as a fallback. It needs a live
`~/.ssh/cm-mdw3.sock`, proxies the TCP leg through it (`ProxyCommand=ssh -S cm-mdw3.sock -W`), and
holds `~/.ssh/cm-mdw2.sock` only while the process runs, so start it in the background:
```bash
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_via_mdw3.exp &   # USER MUST APPROVE DUO
```

**Every routine call: check the socket first, and always pass `-o BatchMode=yes`.**
```bash
if ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu >/dev/null 2>&1; then
    ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<cmd>'
else
    echo "master down - do NOT retry blindly; see throttle note"
fi
```
Why both (learned 2026-09-17 and 09-19): a plain `ssh -S sock host 'cmd'` on a dead master falls back
to a fresh connection, offers keys, tries password auth twice non-interactively and hits `Too many
authentication failures`. `BatchMode=yes` alone protects only when the socket file is gone; if the
master died but `~/.ssh/cm-mdw2.sock` is still on disk, ssh tries the stale socket and then falls
through to a real attempt anyway (`Permission denied (publickey,...)`). A few of those and the host
answers `Connection closed by 128.135.112.69 port 22`: the IP throttle, self-inflicted by a status
check, after which the next expect launch spends a Duo push before it discovers the block.

Rules that follow:
* **When throttled, stop.** The throttle clears by itself in tens of minutes. Wait at least 30 min
  before one single retry. The way around it meanwhile is the midway3 tunnel, which needs a live
  `~/.ssh/cm-mdw3.sock` and therefore its own Duo approval.
* **Run the expect script at most once per attempt**, and never let a monitoring loop retry it more
  than once per tick. Each launch is a Duo push on the user's phone (sent twice by mistake on
  2026-09-18).
* **Any automated credential send is one-shot.** On 2026-09-18 `mdw2_hold.exp` used `exp_continue`
  on the password prompt, re-sent a rejected password on every re-prompt, and one launch spent the
  account's failed-attempt budget. Send the password at most once and treat a second prompt as a hard
  stop (`PASSWORD_REJECTED`); match `Permission denied` and `Too many authentication failures` and
  exit; pass `-o PubkeyAuthentication=no -o IdentitiesOnly=yes` (every key the client offers wastes
  an attempt against `MaxAuthTries`) and `-o NumberOfPasswordPrompts=1`. An expect script has no
  dry-run check: `expect -n -c 'source ...'` still connects.
* **Minimise connections rather than trying to remove Duo.** The training chain insures itself (§1),
  so a status check every few hours is enough. No cluster-side monitor runs;
  `/project/trsosnic/yinhan/STATUS.md` is stale from 2026-09-27 23:54.
* **Recurring work on the cluster.** `scrontab` is disabled, `crontab` is denied, and `pi-trsosnic`
  has no association with the `cron` partition (`sbatch -p cron` is refused). Use a 1-core `broadwl`
  job of at most 36 h (7 days hits `QOSMaxWallDurationPerJobLimit`) that queues its successor at the
  start (`--dependency=afterany:$SLURM_JOB_ID`); one that resubmitted at the end of its loop overran
  the wall and died silently.
* **Laptop sleep is not what drops the master.** With `caffeinate -is` holding `PreventSystemSleep`
  it still died after ~40 min; `mdw2_master.exp` now uses `ServerAliveInterval=30
  ServerAliveCountMax=20 TCPKeepAlive=yes` (~10 min tolerance). When the socket keeps dropping, look
  first for our own load on the login node (`ps -u yinhanw --sort=-rss`, findings 10.6).
* zsh does not word-split unquoted variables: `M2="ssh -S sock host"; $M2 'cmd'` runs silently and
  prints nothing, which is not a connection failure. Write the `ssh` call out in full.
* `timeout` does not exist on this Mac; do not wrap `expect` in it.

**Paths.** `~/project` on midway2 is a symlink to `/project/trsosnic` (not `/project/trsosnic/yinhan/`);
glpG data is at `~/project/yinhan/popepopg_REMD_mdw2/`. `/project` is the same filesystem on midway2
and midway3 (both mount `midway3_cap`), so checkpoints, chain logs and force-field directories under
`/project/trsosnic/yinhan/upside2-md-mdw2/` can be read from midway3 without a midway2 login. Only
`squeue`/`sacct` need midway2: the clusters have separate Slurm controllers and accounting databases,
and `sacct --clusters=all` on midway3 does not see midway2 jobs. Python env on midway2:
`source /software/modules/init/bash && module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1 && export HDF5_USE_FILE_LOCKING=FALSE`

**midway3**:
```bash
ssh -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu '<command>'
```
`~/project` on midway3 is a symlink to `/project/trsosnic/yinhan/` (yinhan, not yinhanw). Load the
python env with `source ~/project/NP-1AO6/env.sh` before any h5py work.

---

## 0a. `broadwl-lc` is empty but unusable: its nodes cannot see `/project` (2026-09-18)

Do not send work there, however idle it looks. Its nodes advertise
`AvailableFeatures=lc,e5-2680v4,64GB,noib`: no InfiniBand, and `/project` and `/beagle3` are served
over it (a probe job on `midway2-0213` reported `PROJECT_MISSING`, `PROJECT_NOT_WRITABLE`,
`BEAGLE3_MISSING`). Every input we use is on `/project`, so the partition is dead for this project,
the same trap as `/cds3`.

The failure is silent: a job whose `--output` path is on `/project` dies at once with
`ExitCode 0:53`, `Elapsed 00:00:00` and no log file, because it cannot create the log. A
`--wrap="hostname"` probe of the same shape fails identically, so it is not the script.

Related traps: `sinfo -F`'s "0 idle" counts a partly filled `mix` node as allocated, so use
`sinfo -p <part> -o "%.8t %.6D %.10C"` for A/I/O/T cores; and `/tmp` on a compute node is node-local,
so a probe writing there "succeeds" and leaves nothing the login node can read (probe with `--output`
under `$HOME`, which is shared). A job_submit plugin silently rewrites `#SBATCH --partition=broadwl-lc`
to `broadwl` while honouring every other directive; command-line `-p` is honoured. Verify placement
with `scontrol show job <id> | grep Partition`.

---

## 0b. Shared Upside deployment on beagle3 (2026-09-09)

`/beagle3/trsosnic/yinhan/upside2-md` is the one tree and binary both clusters use; how and why it
works is in CLAUDE.md, "Shared Upside Deployment". It was deployed by rsync from
`/project/trsosnic/yinhan/upside2-md-mdw2` (`src py obj parameters cmake example` and the install
scripts), excluding `training/`, which is ConDiv campaign data rather than code. Benchmark output
goes to `/beagle3` so it does not compete with `/project`.

Naming: a release writes `parameters/ff_3.0` in both trees. The old `ff_3.0_trained` /
`ff_3.0_trained_rf` directories are still on the cluster but no live script references them
(checked 2026-09-28).

### Other Upside copies (inventory 2026-09-09)

Owned by yinhanw:

| path | size | branch | last touched | disposition |
|---|---|---|---|---|
| `/project/trsosnic/yinhan/upside2-md-mdw2` | 8.3 G | martini-dev | live | ACTIVE, `ff30_glyhb` training |
| `/beagle3/trsosnic/yinhan/upside2-md` | 6.4 G | martini-dev | live | the shared deployment |
| `/scratch/midway2/yinhanw/upside2-md-water-diff` | 13 G | water-diff | 2025-10-30 | idle; see below |
| `/scratch/midway2/yinhanw/upside2-md` | 176 M | master | 2026-01-28 | keep, master |
| `/home/yinhanw/upside2-md` | 256 M | master | 2025-11-19 | keep, master |

**Not ours, never touch:** `/beagle3/trsosnic/upside2-md` and `/project2/trsosnic/upside2-md`
(bayhi), `/project2/trsosnic/software/upside2-md` (nffaruk), `/project2/trsosnic/pengxd/*`
(pengxd), plus copies under baxa, avmolina, ruofan, tobin, schwartznw, yiheng, zonganw,
simoneritchey, amz and bayhi.

`upside2-md-water-diff` is the only deletion candidate, and it is not purely a code copy. Its source
is safe (HEAD `c44a404` is on `origin/water-diff` with 0 uncommitted files), but 11 of its 13 GB is
`example/16.MARTINI/outputs/water_T*`, water-diffusion output from 2025-10-30 that exists nowhere
else. Left in place pending the user's decision.

---

## 0c. Rotamer belief-propagation stopping-test bug (2026-09-26)

**The bug** (reported by John Jumper via Tobin). `NodeHolder::max_deviation` in `src/rotamer.cpp`
took `max(cur_belief - old_belief, dev)` with no absolute value. Beliefs are rescaled each iteration
so the dominant rotamer stays at 1, so an update that only sharpens (others fall) gives dev = 0 and
the solver stops as if converged. It has been there since the original rotamer commit: ff2.1 and
ff30_basin (ff_3.0) were trained with it.

**The fix** is `fabsf(...)` on that line (279). Status by tree:

| tree | status |
|---|---|
| `master` (GitHub) | fixed in `288d1fec`, pushed by the user 2026-09-26 |
| `martini-dev` | fixed; master merged in by PR #48 (`11044069`) 2026-09-26 |
| `/project/trsosnic/yinhan/upside2-md-mdw2` (midway2) | fixed; installed by the ff3.0 gate 2026-09-30 02:24:59 |
| `/beagle3/trsosnic/yinhan/upside2-md` (midway2 + midway3) | fixed; installed by the ff3.0 gate 2026-09-30 02:24:59 |

Both trees carry `BPFIX_INSTALLED` and keep the old binaries as
`obj/*.bak_pre_bpfix_20260930-022459` (gate log `ff30_basin/ff30b-gate_49131949.out`, `[bpfix]`
lines).

**Static-frame result (local, run 2026-09-26, rerun on fresh frames 2026-10-02).** 8 ConDiv training
proteins (2xf6, 3h36, 1afh, 2lpn, 3mwz, 1a62, 3dm8, 1w4s; 50-146 res), ff2.1, 91 frames each (native
+ 30 each from 2000-unit runs at T 0.80/0.95/1.10), buggy and fixed library at the production tol
1e-3 against the fixed library at tol 1e-6 (max_iter 20000). The 10-02 rerun: the bug cuts the solve
short in 338 of 728 frames (46%), by 2 iterations in 243, 4 in 75, 6 in 18, 8 and 12 in one each.
Energy error |E - E_ref| median 7.5e-4 E_up (fixed 4.3e-4), 99th percentile 0.0095, max 0.157 (1afh
frame 21: bug stops at 8 iterations, fix at 12, fix error 3e-4); the fixed solver's own max is 0.044
(1afh frame 10, both stop at 10, the 1e-3 tolerance alone). The 09-26 frames had max 0.029; the 0.157
frame is a real early stop that those frames did not sample, still 0.2 kT at T 0.80. RMS force error
/ RMS force median 8.3e-5, max 2.9e-3 (fixed 4.6e-5, 5.3e-4). ConDiv contrastive rotamer gradient
(native minus decoy mean of dE/d pair parameters, summed over the 8 proteins; recomputed 10-05 from
these frames by `make_slides.py`): bug - fix 1.05e-3 of its norm, turned by 1e-3 rad, against a
sampling noise of 7.2% (standard error, interleaved halves; first/second halves give 10.7%), so 68x
below the noise. Per coordinate (17327 live of 21600): below 0.08 of its own noise in 99%, above it
in 0.03% (CYS-ASN, TYR-VAL, CYS-TYR, one distance knot each, max 2.9x, under 2% of that gradient).
Scope: local, at ff2.1's parameters, rotamer pair group only; it gives the direction further
training would take, not the end point of a long retraining (a bias accumulates, the curvature was
not measured). The 1.4-2.5e-4 recorded here before (definition not recorded, likely the 09-26 frames)
does not reproduce on these frames. Energy bias one-signed: E_bug < E_ref in 95% of frames,
mean -0.0007 E_up. Conclusion given to the user: ff2.1 and ff3.0 do not need retraining.

**Simulation validation for Tobin: all 16 arms complete 2026-09-29 03:50, analysed; not yet
reported to Tobin (the user's call).** Output `bp_validation/analysis_20260929.txt`, 25-32k
frames/replica after equilibration. Result:
* WWdomain is the only protein with a consistent difference: at T 0.764-0.879 the bug arms are
  less native, Q lower by 0.02-0.04 (z -1.8 to -12), RMSD higher by 0.1-0.3 A, E higher by 1-3; T_mid
  0.859 bug vs 0.861 fix (~0.7 K). Per seed (checked 2026-09-29), both bug seeds lie below both fix
  seeds in Q at every T from 0.700 to 0.879; at 0.852 bug 0.511 / 0.505 against fix 0.547 / 0.550.
  At 0.700 the gap is one bug seed (0.842 against 0.89).
* NTL9 goes the same way at low T (Q 0.25 bug vs 0.30 fix, z -2 to -4; RMSD +1.0-1.3 A, z 1.1-1.8),
  but per seed that rests on one arm: fix s2 sits at RMSD 6.5-6.9 A while the other three sit at
  8.0-9.0 A. Its folded population is low under ff_2.1 here and T_mid is not resolved (0.821 vs
  0.815). Not resolved.
* proteinG and homeodomain: not resolved, |z| < 1.5 through the folded and transition range; T_mid
  0.861 vs 0.868 and 0.931 vs 0.929, inside the seed spread.
* Unfolded-T z of 2-5 sit on seed sd of 1e-4 in Q and 0.01-0.03 A in Rg; the differences are that
  small, and with 2-dof sigmas across 224 cells such z values are expected by chance.
* With no effect, z is Student-t with 2 dof, and the 224 cells match that null (checked 2026-10-01:
  |z| > 2 in 41 against 41.1 expected, |z| > 9.92 in 2 against 2.2). No single z, WWdomain's -12
  included, is evidence by itself. WWdomain rests on the sign holding through its transition and on
  the per-seed separation: suggestive, not established at two seeds. Its implied ddG (~0.15 E_up,
  van 't Hoff from dT_mid 0.002) is ~200x the median static-frame |dE| (7.5e-4) and matched only by
  the single worst frame (0.157); not explained.
* Speed: the fix costs 0-5% (proteinG 17.6 vs 18.6, WWdomain 29.3 vs 29.8, NTL9 28.0 vs 28.0
  time units/s); node-to-node variation is not controlled.

Protocol: ff_2.1, native start, 14-replica REMD with the Peng benchmark protocol (Table S2 ladder
and duration, dt 0.009, frame 100), proteinG / homeodomain / WWdomain / NTL9 x {bug, fix} x seeds
{1, 2}. Trajectories of the two binaries decorrelate within a few hundred time units, so the test is
statistical: bug-minus-fix against the seed-to-seed spread (z = delta / sigma_seed). All 16 checked
healthy at the first 100 frames (C-N 1.33-1.34 +- 0.13 A, KE/1.5kT 0.99-1.00).
* dir `/beagle3/trsosnic/yinhan/bp_validation`: `obj_bug/`, `obj_fix/` (deployment source built
  twice, differing only in the fix; `obj_bug` reproduces the old deployed binary bitwise), `bp_run.py`
  (copy of `ff3_benchmark/bench_run.py`: binary, seed, ff_2.1 native), `bp.sbatch`, `submit_all.sh`
* logs `logs/<prot>_<build>_s<seed>_<jobid>.out`, data `runs/<prot>_<build>_s<seed>/`, completion
  markers `runs/*/COMPLETE`
* analyse on the midway2 login node after `source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh`:
  `python3 analyse_bp.py [protein ...]` gives per-T Q, CA-RMSD, Rg and E per build, bug-fix vs the
  seed spread, melting midpoints and speed. Two seeds give only a rough sigma: |z| ~ 1 is "not
  resolved"; only |z| well above 2 across the ladder is an effect.

**Figure for Tobin (done 2026-10-02 15:48).** `~/Downloads/bp_validation.{png,pdf}` on the Mac
that made it, and `bp_validation/static/bp_validation.{png,pdf}`. `plot_bp.py`, 11 panels: a-c
static frames (ECDF of |E - E_ref|, ECDF of relative force error, histogram of iterations cut short),
d-g Q(T) per seed for the four REMD proteins, h-k Q_bug - Q_fix with the seed sigma. Red = bug, blue
= fix throughout; SciencePlots 2.1.1 style files bundled in `styles/` (identical rcParams to the
global CLAUDE.md path, so the PNG is bitwise the same). T_mid in d-g is from the two-seed mean
Q(T), so it differs by up to 0.001 from the per-seed values above.
* REMD arrays: `analyse_bp.py` now also writes `obs_<protein>.npz` (per arm and T: Q, RMSD, Rg, E;
  original kept as `analyse_bp.py.bak_pre_npz`). Rerun 10-02 on the midway2 login node (user's call:
  compute-node job 49143602 had a 20:11 start estimate and was cancelled; 4 niced processes under
  `ulimit -v 4 GB`, 2.5 min); `analysis_rerun.txt` is identical to `analysis_20260929.txt`.
  `dump.sbatch` is the compute-node version of the same run.
* **`bp_validation/static/` holds everything the Mac had** (325 files, md5 manifest identical to
  the Mac's `scratchpad/bp_validation/` on 10-02 16:10), so nothing of this work is on one computer:
  * figure: `plot_bp.py`, `styles/`, `remd/obs_*.npz`, `bp_validation.{png,pdf}`
  * static-frame harness, in run order: `make_src.py` (writes `src_bug/` and `src_fix/` from a
    source tree: the `read last_iter` BP-iteration readout in both, and the fabsf removed in
    src_bug), `build_libs.sh` (the Apple Silicon `install_M1.sh` recipe into `obj_bug/`,
    `obj_fix/`), `gen.py <code>` (frames), `eval.py <bug|fix> <tol> <max_iter> <tag>` (run as
    `bug 1e-3 1000 bug`, `fix 1e-3 1000 fix`, `fix 1e-6 20000 ref`)
  * as used: `src_bug/`, `src_fix/` (complete source trees), `obj_{bug,fix}/libupside.dylib`
    (arm64 macOS, loadable by `eval.py` on another Apple Silicon Mac without rebuilding)
  * results: `sys/<code>/{frames,init}.npy`, `base.up`, `res_{bug,fix,ref}.npz`,
    `T*.up.output`; logs `gen_<code>.log`, `{bug,fix,ref}.log`
  * not copied, and not needed: the frame-run trajectories `sys/*/T*.up` (their frames are in
    `frames.npy`), `eval_*.up` (rewritten by `eval.py`), the cmake build trees.
  Replot from any computer by copying `static/` and running `python plot_bp.py` in it.

---

## 0d. midway2 node exclusions: one law (2026-10-08)

**The only midway2 node list is `/project/trsosnic/yinhan/slurm/midway2.args`** (user, 10-08):
`--partition=broadwl --exclude=midway2-[0003,0010-0011,0027,0033-0035,0037,0060,0080,0085,0088,0103,0116,0151,0342-0345]`.
Each node's record of failures is in `README.md` beside it; §8 has the history. midway2-0088 was
added 10-09 04:12, midway2-0103 and 0116 10-09 10:30, and midway2-0151 10-13 10:29 (two NODE_FAILs
10-10: 49218464 at 09:52:43, link 49216034 at 20:52:43), each at the user's word (backups
`*.bak_pre_0088_20261009`, `*.bak_pre_0103_20261009`, `*.bak_pre_0151_20261013`). For 0151 the
wrapper refused `--test-only -w midway2-0151` and accepted 0152; no job was pending. For the second, `--test-only` accepted the
list, the wrapper refused `-w midway2-0103` and `-w midway2-0116` as it refuses 0088 and accepted
0110, and `update_pending.sh` gave all 19 pending jobs the new list. Before 10-08 the
list was copied by hand into 13 `slurm.args`, `submit_remd.sh`, `bench.sbatch` and a dozen static
`#SBATCH` lines in three versions, and a queued chain never saw a node added later.

Where it is read, always at submission time:
* Every midway2 `slurm.args` is a symlink to it (the round-4 runs, ff30_glyhb, ff30_gly,
  ff30_basin, `ff3_selection/`, `polygly/prod/`; backups `slurm.args.bak_pre_law_20261008`).
  `train_chain.sbatch` reads `slurm.args` at every link, so a running chain takes a change at its
  next link; its gate, `validate_ff.sh`, `submit_new.sh` and `prod.sbatch` read it too.
* `popepopg_REMD_mdw2/submit_remd.sh` takes `$(cat <law>)`; every REMD block resubmits through it.
* `/beagle3/.../ff3_benchmark/bench.sbatch`: on midway2 a chunk resubmits with the law; on caslake
  it carries its own partition and `--exclude`. Its static `#SBATCH --exclude` is gone.
* `~/bin/sbatch` on midway2 (a link to `slurm/sbatch_wrapper`): `~/bin` comes before Slurm on
  PATH in login shells, ssh commands and jobs, so every other submission from a midway2 host gets
  the law's `--exclude` first, and a script's own `#SBATCH --exclude` is overridden. sbatch keeps
  the last `--exclude` given, so a caller's own `--exclude`, or the real
  `/software/slurm-current-el7-x86_64/bin/sbatch` by path, is the deliberate way around it.
* Left alone: ff30_glyprobe's `slurm.args` (extra flags; cancelled, no gate), the caslake runs'
  files, and finished one-off scripts with static lists (`ff3_benchmark/scoring/`, `gly_peptides/`,
  `polygly/collapse.sbatch`, `checks/`), which the wrapper overrides on midway2.

**Adding a node.** A node is excluded for its own record (a second failure on it), not for being
near another. Edit the one line of `midway2.args`, add the node's row to `README.md`, then run
`bash /project/trsosnic/yinhan/slurm/update_pending.sh` on a midway2 login node: it sets every
pending job's ExcNodeList to the new list and prints each. Running jobs keep their nodes. A new
run's `slurm.args` is `ln -s /project/trsosnic/yinhan/slurm/midway2.args slurm.args` (`cp -a`
keeps the link; plain `cp` copies the text, which then drifts).

**Checked 10-08 09:35-09:40** (no test jobs; `--test-only` and stubs):
* The law parses.
* Through the wrapper, `-w midway2-0027` and `-w midway2-0343` are refused and `-w midway2-0026` is
  accepted. The real sbatch by path accepts 0027, so the refusal is the wrapper's. A caller's own
  `--exclude` replaces the law.
* A script with its own old `#SBATCH --exclude` is still refused 0027.
* All 13 `slurm.args` resolve to the law.
* With a stub sbatch, `submit_remd.sh` passes the law, and `bench.sbatch`'s resubmit takes the law
  on a midway2 host and its own flags on a midway3 one.
* `update_pending.sh` set all 13 pending jobs. Their dependencies were intact, although squeue
  briefly showed reason "None".
* End to end on real submissions: bs_01's panel 49208587 (11:43) and the first chain successor
  submitted from a running link, ff30_gdepth_fz's 49214121 (21:35), carry 0027 in ExcNodeList.

## 1. Current jobs

Snapshot **2026-10-21 13:50 CDT, verified live against `squeue` on midway2** (MacBook Pro; no watch runs). 46 jobs on
midway2: released ff2.1's 30 Peng arms, its mdCATH panel and TM job and its collector; ff3.1.1 and ff3.2.1's two links,
their successors and their collector; the 8 glpG chains. midway3: nothing known to run (socket closed this pass).

**Line 1: deterministic vs probabilistic side chain**
- **Released ff2.1 on the three side-chain tests** (user 10-21: "run validation with released ff2.1 so we know if
  we should keep ff2.1 or not"; plan.md Phase 27; findings 1.37). Submitted 11:43 10-21 by
  `checks/ff21_released_20261021/submit.sh` (ids in `submitted.log` there); files `ff3_selection/params/ff21_released`
  (md5 equal to `/beagle3` `ff_2.1`), each test's own protocol, binary and domains at dt 0.009.
  - 44 domains: `p44_ff21_released` 49438695 and its TM job 49438697 done 12:59. Released ff2.1 folds 0.603 against
    the control's fixed point 0.430 and no-BP's 0.363; mean TM 0.681 / 0.650 / 0.623 (`p44_dt009_folded.txt`,
    `tm_p44_*.txt`).
  - mdCATH: `mdc_ff21_released` **49438694** R (log `ff3_selection/logs/ff3sel_mdc_49438694.out`, runs
    `runs_mdcath/mdc_ff21_released/`, ~6 h), then **49438696** (afterok) writes
    `tm_mdcath_mdc_ff21_released_vs_<mdc_ctrl_03|mdc_nobp_03>.txt`.
  - Peng: 30 arms `b_ff_2.1_<prot>_<kind>` (first chunks 49438698-49438727, each resubmitting itself in 10.5 h chunks;
    logs `ff3_benchmark/logs/<prot>_<kind>_ff_2.1_<id>.out`, runs `ff3_benchmark/runs/<prot>_<kind>_ff_2.1/`; lambda's
    two arms of 10-04 reused), about a week. Collector `ff21_collect` (every 6 h, queued **49438731**): scores each
    finished arm into `ff21_nobp_dt009/setup/score/ff_2.1/` and, with all 32, writes
    `checks/ff21_released_20261021/peng_ff21.txt` ending DONE: released ff2.1 against the control, no-BP, ff3.0 and
    ff3.1 (all scored), beside Peng's published FF2.
- **Both side-chain trainings are finished** ("Finished runs"); no jobs.

**Line 2: Ramachandran map training**
- **ff3.1.1 and ff3.2.1: ff2.1 with only its depths trained** (user 10-21: "apply the fix and train ff3.1.1 and
  ff3.2.1", then "try to resolve the accuracy lost bug as well in this run"; plan.md Phase 29 and its revised
  decision; findings 1.34, 1.37). `ff311_depths_only` (`p311_EE`, from NDRD's map) **49439398** and
  `ff321_depths_only` (`p321_EE`, from ff3.0's map) **49439399**, R from 13:42 10-21 on midway2-[0232-0243] and
  [0246-0257], 54-step links with successors queued afterany. Built by `checks/ff311_20261021/build_depths_only.sh`
  (logs `build_ff311_depths_only.log`, `build_ff321_depths_only.log`): ff2.1's files, md5-equal to the release;
  `trainer/rama_basin.py` the repo's with MODE `pooled` (md5 `54816673...`, shared width 0.5); every Adam rate zero
  (`freeze_groups.py`; the trained rates kept as `initial_checkpoint.pkl.trained_rates`); ff3.1's rerun's minibatch
  order and held-out set; `compare_init_pooled.py`: equal to the source in every other field. Target 76, no gate,
  stop for the user; step 76 ~02:00-10:00 10-23. Logs `<run>/condiv-train_<id>.out`. Collector `pooled_collect`
  **49439400** (first pass 16:42, then every 4 h): `check_step.py` per step into `checks/ff311_20261021/<run>_steps.txt`,
  the steric check into `checks/ff31_steric/<run>_<step>.txt`, the epoch-end panels by `submit_new.sh`, and
  `checks/ff311_20261021/status.txt`, ending DONE at step 76. **First read: the epoch-0 panels `p311_00` and
  `p321_00` against `ff21_released`** (~02:00 10-22), the test of the accuracy fix.
- **Finished, awaiting the user's review:** the fixed-engine reruns `ff31_depth_data_angfix` (`a31_EE`) and
  `ff32_depth_data_angfix` (`a32_EE`), step 76 10:35 / 10:43 10-21: margins +0.109 and -0.002 at step 75; a `STERIC
  FLAG` on both from the free beta' + pPII' alone (0.0052, 0.0053 against steps 0-18's 0.0045), clash ratios small;
  the ff3.2 rerun at a fixed point by every test (findings 1.36); their 44-domain panels at step 76 0.420 and 0.432
  against ff2.1's 0.615 (findings 1.37).
- **The 8 corrected glpG chains** (plan.md Phase 18, findings 4.17), each on its own node; logs
  `popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<jobid>.out`, run dirs `popepopg_REMD_mdw2/<V>.<ff>/`, seeds
  `seeds/<V>.<ff>.up` (built and verified in `checks/placefix_20261017/`). Blocks 1 and 2 closed 18:52-20:39 10-18 and
  07:55-14:56 10-20 (KE/1.5kT 0.989-1.011; block 1's TM4 in findings 1.35). **Block 3/5, all R**, logs clean at 11:10
  10-21 (KE/1.5kT 0.988-1.014): ff3.0 79HIS_S115T 49428215 from 07:55 10-20, 79ALA 49428473 from 08:23; ff3.1 79HIS
  49432319 from 14:29; from 14:36 ff3.1 79ALA_S115T 49432546, 79ALA 49432586, ff3.0 79ALA_S115T 49432556, 79HIS
  49432585; ff3.1 79HIS_S115T 49432786 from 14:56. They run force fields that carry the retraining loss (findings
  1.37); stopping them at block 3's end is the user's. Owed: the TM4 check at block 2's end and the HDX on block 2.

**Line 3: Upside 3 design (debug stage).** No cluster jobs; local work (plan.md Phase 24 "Resume here", upside3.md).

**Line 4: Boltz-2 HDX sampler (midway3)**
- **Tag by step-scale test** (user 10-21, "Yes."; `diffusion_dynamics.md` Phase 0): **61231966**
  `boltz_tag_stepscale_v100` (1 V100, `--partition=gpu --qos=gpu`), submitted 08:41 10-21; its log ends with the eight
  cells' `check_samples.py` summaries and `[done]` at 09:02 (Slurm state not read). Script
  `$B/tag_stepscale_test.sbatch` (`$B = /beagle3/trsosnic/yinhan/boltz_hdx`; local copy `scratchpad/boltz_hdx/`): ubiquitin
  and EHEE_rd2_0005 at {no tag, `md`} x {`step_scale` 1.5, 1.0}, 100 samples each, seed 2. Log
  `$B/logs/boltz_tag_stepscale_v100_61231966.out`, data `$B/runs/boltz_tag_stepscale_v100_61231966/`. **Owed:** the
  readout (peptide bond scatter and clash counts across the four cells) in `diffusion_dynamics.md` Phase 0.

**Shared**
- **The engine fix** (findings 3.12; plan.md Phase 21): `$P/obj/upside` and `libupside.so` from
  `checks/angfix_20261018/obj` (md5 `7136d99a...`, `3d9c6366...`), installed 08:35 10-18 with backups
  `*.bak_pre_angfix_20261018`; `$P/src/interaction_graph.h` patched (backup alongside). `/beagle3`'s build is not
  changed. Energies, forces, simulations and every derivative but the pair term's angular slots are bitwise those of
  the old build.
- Disk: `/project/trsosnic` 1.4 T free at 11:10 10-21 (1.6 T at 10:24 10-19; "Disk" below); `/beagle3/trsosnic`
  1.2 T at 22:30 10-20.
- **What happens next** (estimates; Slurm's start times move):
  - glpG block 3 closes ~18:00-19:00 10-21 (ff3.0 79HIS_S115T, 79ALA) and ~00:00-01:30 10-22 (the other six).
  - Released ff2.1's mdCATH TM comparison ~20:00 10-21.
  - ff3.1.1 / ff3.2.1: epoch 0 ends ~23:00 10-21; its panels are submitted by the collector's next pass and read
    ~02:00 10-22; step 76 ~10-23.
  - Released ff2.1's Peng table (`peng_ff21.txt`) in about a week.

Finished and cancelled rows are deleted; their lessons are in §8 and "Finished runs".

| JobID / where | what | state | next action |
|---|---|---|---|
| **block 3: 49428215, 49428473, 49432319, 49432546, 49432556, 49432585, 49432586, 49432786 R** (midway2; block 1 49411315-49411322 R 09:58 10-17 to 20:39 10-18; block 2 49417749-49418152 to 14:56 10-20; 49411323-49411324 stopped 08:33 10-18) | **The corrected glpG REMD chains** (user 10-17: "fix the bug and rerun the affected simulations"; plan.md Phase 18; findings 4.17): the side-chain placement fixed in the preparation and every live seed, and the 10 chains rerun from block 1, 28 replicas T 0.70-0.90, five 36 h blocks (`submit_remd.sh`). ff_3.0_gdepth: 79HIS 49411315, 79HIS_S115T 49411316, 79ALA 49411317, 79ALA_S115T 49411318; ff_3.1_depth_data: 79HIS 49411319, 79HIS_S115T 49411320, 79ALA 49411321, 79ALA_S115T 49411322; 79HIS ff_2.1_ctrl38 49411323 and ff_2.1_nobp 49411324 (the deterministic side chain, `single_state.py`: one row per side-chain residue at the kept state, 1-body 0) were stopped by the user 08:33 10-18 (STOP files first, directories kept): their step-37 releases are superseded (plan.md Phase 21). Seeds `popepopg_REMD_mdw2/seeds/<V>.<ff>.up`: base seeds repaired (`fix_placement.py`), per-force-field seeds by `patch_glpg.py` on the repaired base as `validate_ff.sh` makes them (each equal to its old seed repaired); built and verified in `checks/placefix_20261017/` (placement at `rcf[li]` to 3.6e-6 A, 1-body exact, side-chain/lipid rows on, only the intended datasets changed). Logs `popepopg_REMD_mdw2/logs/remd.<V>.<ff>.<jobid>.out` | Block 1 closed in all eight 18:52-20:39 10-18 (KE/1.5kT 0.989-1.008; TM4 read 21:40, findings 1.35: TM4 holds, no glycine at phi > 0 at T 0.70); block 2 queued at once, all carrying the law. **First chunks healthy** (read 10:40): calibration 615-682 ms/step (the old chains 645-694; nobp fastest, 615); at chunk 1 step ~1,500-1,950 replicas 0, 13, 27 at Rg 20.0-21.8 A, backbone H-bonds 184-221, protein potential -731 to -855 (T 0.70), -589 to -750 (0.79), -471 to -610 (0.90), nobp least negative (-731, -589, -471); logs clean. Block 2 closed in all eight 07:55-14:56 10-20 (KE/1.5kT 0.989-1.011), block 3 queued at once, R, logs clean at 22:30 10-20 | the watch: KE/1.5kT at each chunk end, then first chunks (calibration, replicas 0, 13, 27: Rg, protein potential, KE against the old chains'), then each block end; TM4 under every force field at each block end (`tm4_ff31.py`, `gly_flip_rungs.py`, DSSP) and the HDX on block 2 (recipe below) |
| **49438694 R, 49438696 PD; Peng arms resubmitting (first chunks 49438698-49438727); collector 49438731** (midway2) | **Released ff2.1 on the three side-chain tests** (plan.md Phase 27): the mdCATH and 44-domain panels at dt 0.009 and Peng's benchmark, on the same protocols, binaries and domains as the retrained arms, to decide whether to keep ff2.1 itself | submitted 11:43 10-21; 44 domains done 12:59 | read `checks/ff21_released_20261021/tm_mdcath_*.txt` (~20:00 10-21) and `peng_ff21.txt` (DONE, ~a week); then plan.md Phase 27 step 4 |
| **49439398, 49439399 R; successors afterany; collector 49439400** (midway2) | **ff311_depths_only and ff321_depths_only (ff3.1.1, ff3.2.1)**: ff2.1 with only its basin depths trained under the pooled prior, from NDRD's and ff3.0's map (plan.md Phase 29) | first links R from 13:42 10-21 | the epoch-0 panels `p311_00` / `p321_00` against `ff21_released` (~02:00 10-22); `checks/ff311_20261021/status.txt`; at step 76 the maps and `depth_data_vs_prior.py`; stop for the user |

### TM4 HDX on ff_3.1_depth_data: the recipe (done 10-17, findings 5.5; user 10-14 02:10, "Do the analysis on ff3.1 as well")

The ff_3.0_gdepth run of 10-14 ("Finished runs") repeated on ff_3.1_depth_data's block 2, 79ALA
then 79HIS; not the S115T chains. Block 1 was 49350831 / 49350829; block 2 is 49402903 / 49402866
(R from 23:37 and 22:54 10-14), the `<block-2 jobid>` of steps 1 and 3. Scripts in
`C=/project/trsosnic/yinhan/checks/hdx_ff30_20261014/` (md5: `hdx_block2.sbatch` `019af27d`, `prechecks.py`
`7ede2a7a`, `extract_block.py` `0f8781f6`, `readout.py` `5e78ea57`, `martini_hdx_project.py` `8317fbe1`, the
repo's; `plot_tm4.py` `1c91ddbd`; `*.ff30_run_20261014` are the ff3.0 run's own copies).
1. **When:** the chain's block-2 log (`popepopg_REMD_mdw2/logs/remd.<V>.ff_3.1_depth_data.<id>.out`) ends
   with `[remd] block complete`. Block 1 ends ~23:30 10-14 to ~01:00 10-15; block 2 takes ~36 h after it starts.
2. **Block 2's groups.** n_k = `grep -c "^avg_kinetic_energy/1.5kT"` in block k's log (one calibration
   run plus its chunks). `output_previous_0` is the seed's own group, so FIRST = n1 + 1 and LAST = n1 + n2.
   Check that group FIRST of `run.0.up` has 287 frames whose `time` ends at 18.0 (the 2000-step
   calibration); read live files only with `HDF5_USE_FILE_LOCKING=FALSE`. Measured on ff_3.0_gdepth 79HIS:
   n = 17, 18, 19, so block 2 = 18-35 and block 3 = 36-54, each first group 287 frames to t 18.0.
3. **Pre-checks on midway2 login1** (read only, ~2 min): `cd $C; source
   /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh; export HDF5_USE_FILE_LOCKING=FALSE; ulimit -v 4000000;
   nice python3 prechecks.py ff_3.1_depth_data glpG-RKRK-<V> FIRST LAST <block-2 jobid> | tee
   prechecks_ff31_<V>.txt`. It asserts the log's run count equals the group count. ff_3.0_gdepth gave
   KE/1.5kT 0.998-1.002, swap acceptance 0.36-0.53, no round trip.
4. **Submit from midway3** (caslake: the only place with mdtraj, pymbar and matplotlib; the user's
   "Run it on midway2 if possible" could not be met). Tell the user, then open the midway3 socket once
   (`expect scratchpad/mdw3_master.exp`, one Duo push), and per variant:
   `sbatch --job-name=hdx_<V>_ff31 --output=$C/hdx_glpG-RKRK-<V>_%j.out --partition=caslake
   --account=pi-trsosnic --exclude=midway3-0014 --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
   --time=08:00:00 --export=ALL,PDB_ID=glpG-RKRK-<V>,FF=ff_3.1_depth_data,FIRST=<g>,LAST=<g>
   $C/hdx_block2.sbatch`. It refuses to run if its work dir exists. ff_3.0_gdepth's took 5 min each.
5. **Check:** the log has `block 2 copies verified: 28/28`, `replicas done: 28/28` and `HDX DONE`, no
   `FAIL rep` or `Traceback`. Work dir `popepopg_REMD_mdw2/<V>.ff_3.1_depth_data/hdx_sc600x12_ff_3.1_depth_data/`
   (~30 G each; `results/` PS, `protein_only/results/` PS_protein, `results/<V>_tm4_readout.txt`). Quote
   only censored-aware dG (the readout's "> bound"), never the standard plots' unclipped values.
6. **Read and report:** TM4 130-151 at 298 K (T 0.85) against ff_3.0_gdepth's (findings 5.5: core
   censored above 5.85, only 130-134 and 149-151 resolved). Figures: `plot_tm4.py NPZ_DIR OUT_DIR` on a
   computer with SciencePlots (fix its `SP` path), with `<V>_PS.npz` = `results/<V>_dG_profiles.npz` and
   `<V>_PS_protein.npz` = `protein_only/results/<V>_dG_profiles.npz`; adapt it to draw ff3.1 beside ff3.0.
   The ff3.0 figures are on the Mac Studio in `~/Downloads/glpG_TM4_HDX_ff3.0_gdepth_block2/`. Then
   findings 5.5, plan.md Phase 16, this file and `WATCH_STATUS.md`.

**The three H-bond-frozen arms were stopped 10-09 23:09 by the user** ("meaningless at this point:
the general training at dt 0.009 already made a force field that keeps TM4 stable"): successor
49214121 cancelled first, then ff30_gdepth_fz's running link 49206070 and the pending successors
49206071 (ff30_bio_fz) and 49214161 (ff21_ctrl_fz); re-checked 10-11 23:30, no job runs from an
`_fz` directory. They reached step 58, 45 and 10, every step frozen at margin +0.1919, 24 of 24,
KE/1.5kT at most 1.017. **The three run directories were deleted entirely 10-11 23:33-23:40 at
the user's word** ("They are wrong and no longer needed"), with their first-link logs
`training/condiv-train_{49201673,49204264,49206070}.out`; file lists in
`checks/disk_reclaim_20261011/manifest_<run>_{steps,all}.txt`. Their derived products (extracted
ends bz_00-02 and dz_00-01, their panels and TM4 sets) were deleted 10-11 23:45 with the BioEmu line
("Finished runs"), and `submit_new.sh` no longer names them. Their TM4 and panel results stand in
findings 1.21-1.25. The dz_01 panel 49218034 was cancelled 23:16 (user).

**Local glpG TM4 test** (the hybrid glpG 79HIS seed in POPE/POPG dry-MARTINI, one checkpoint
patched in, T 0.80, 4000 tu, 12 seeds; `scratchpad/ff3_local_test`, record in
`checks/r4_epochs/tm4_local`). **Retired 10-11 23:10 (user)**. The sets kept in `tm4_local/runs_cov/`
(and the Mac's `runs/`) are the depth path's and the controls': ff21_released, gdepth_start, d9_00-03,
c9_00-01, fp_e00, and on the fixed table ff21_released_sc, gdepth_start_sc, d9_01_sc, d9_02_sc,
d9_03_sc, fa_00_sc, fd_00_sc. The BioEmu line's and the frozen arms' sets (b9_EE, bs_EE, ff21_bioT1_6,
bz_EE, dz_EE) were deleted 10-11 23:45 with their runs, and ff30_gdepth_si's (ds_00, ds_01_sc,
ds_02_sc) on 10-12; their `tm4_compare_cov_*.txt` tables stay, and findings 1.21-1.26, 1.31, 4.15
hold their numbers. Every set predates the placement fix (findings 4.17); their cluster `.up` trajectories in
`runs/` and `runs_cov/` were deleted 10-18 (user; logs and tables kept).
* **Invalid, kept as a record:** every set before 10-07 11:22. That covers ff21_released,
  ff21_bioT1_6, gdepth_start, fp_e00 and b9_00 at 12 seeds, and b00, b01, b01m12, d00 and d01m09 at
  3 seeds (cluster `tm4_local/runs/`). Each paired its force field with the retired FF1-form
  ff_3.0's coverage tables (findings 3.11). Stopped sets (b9_01 at t 1200 on the MacBook Pro; b9_01
  and d9_00 at t ~1210 on the Mac Studio, `runs_stopped_20261007/`) are in no table.
* **Total-potential jumps above 3000** persist on fixed inputs: 1 to 4 seeds in every set (findings
  1.24), with no frame dropped and KE/1.5kT normal. Their cause is not identified.
* `tm4_local.py` also tables an unfinished set; `tm4_compare.py` refuses one.

**The lambda ff_2.1 benchmark is complete** (all four arms at target: WT de novo 10-04 14:17, WT native
16:17, G46A/G48A de novo 10:37, G46A/G48A native 10-05 01:24; `/beagle3/trsosnic/yinhan/ff3_benchmark/runs/lambda_*_ff_2.1/`):
WT native unfolds from ~1 M tu (10.2-10.3 A) and de novo never folds, helix 3 last-block aR 0.61 / 0.57;
G46A/G48A keeps helix 3 (0.97 native, 0.80-0.88 de novo), native on a ~6.9 A plateau, de novo never
folds (findings 11d). No cluster jobs run.

**ff30_glyhb was stopped 10-04 21:58 by the user's decision** (link 49145374 at step 63, successor
49161468 and panel h02 49174801 cancelled; the chain script submits its successor only at a link's
start, so the successor was cancelled first). Newest checkpoint
`$P/training/ff30_glyhb/run_output/epoch_03_minibatch_06`; at step 63 shared margin +0.031,
glycine margin -0.309, hbg [+0.083 +0.204 -0.257]. Panels h00 / h01 (37 domains): folded 0.460 /
0.469, helix -0.029 / -0.030, gly_helix -0.130 / -0.122, gly_left -0.075 / -0.054, against ff2.1's
0.615, -0.016, -0.114, +0.032; "RELEASE HELD: worse than ff2.1 in helix" (findings 1.17). The local
glpG test that showed its epoch-2 checkpoint flipping GLY143 (findings 1.19) used pre-fix inputs
(findings 3.11), so that result is invalid.

The gradient split is closed (12:52): threshold test 49139877 finished 24 of 24, both threshold
statements give the same DSE (findings 1.17); the held remainder (49138763_[88-95],
49139665_[25-47]) was cancelled.

**ff30_gly was stopped 11:51:19** (49135913 and its successor 49136311 cancelled, approved plan) in
epoch 3 minibatch 8; its newest checkpoint is `epoch_03_minibatch_07`, so it can be resumed if ever
needed. Its panels e00, e01, e02 stay as the comparison.

**The midway2 tree carries the glycine offsets since 2026-10-02 11:46** (deploy dir
`/project/trsosnic/yinhan/checks/hbg_deploy_20261002`, logs `deploy.log`, `finish.log`,
`parity.log`): `src/hbond.cpp`, `py/upside_config.py`, `training/{ConDiv,extract_ff,check_step,
patch_glpg}.py`, `training/README.md`, each with a `.bak_pre_hbg_20261002`; `obj/{upside,
libupside.so}` replaced by rename from `obj_hbg/` (backups `.bak_pre_hbg_20261002`). On 1ga3 the new
build is bitwise the old one with ff2.1's 12-entry file (engine and a 200-unit run), and zero
offsets are bitwise 12 entries. `py/upside_config.py` there had been stale (it still wrote the
`_ALL` sheet arrays removed 09-24); it is now the repo's. The `/beagle3` shared deployment is NOT
updated yet: the benchmark arms run from it, and it must be synced from `$P` before ff3.0's
benchmark or glpG runs.

**dt 0.009 restart (user, 2026-10-06 09:46; backups `.bak_pre_dt009_20261006`)**:
* `$P/training/ConDiv.py` (the ff30_bio trainer) has `time_step=0.009` and a docstring item for
  it; `ff30_gdepth_dt009/trainer/ConDiv.py` the same. The repo's `training/ConDiv.py` matches.
* `$P/training/train_chain.sbatch` has STEPS_PER_LINK 54, against 90 at dt 0.015; the repo's the
  same.
* `ff3_selection/submit_new.sh` maps the new runs to tags `b9_EE` and `d9_EE`, and since 10-06 23:58
  the control `ff21_ctrl_dt009` to `c9_EE` (backup `.bak_pre_ctrl_20261006`).

**Watch**: the dt 0.015 runs' cron `02ff3878` was deleted with them (10-06 09:41); see "One watch
at a time" for the current one. Epoch panels from `submit_new.sh` still go to broadwl
(`ff3_selection/slurm.args`); one sent to caslake needs `--mem` (panel.sbatch sets none). The gate
(`after_training.sbatch`, 2 CPUs, no `--mem`) projects the same start on caslake with or without
`--mem=4G` (`--test-only` 19:08: both 10-06 06:36; a fresh 8-node link 09:14), so it needs none.
Shared-script changes for round 4 (backups `.bak_pre_gdepth_20261005` and `.bak_pre_bio_20261005`):
`$P/training/extract_ff.py` (writes a depth-training checkpoint's trained library; finds the run's own
trainer for an initial checkpoint; byte-identical output for the other runs, tested),
`$P/training/check_step.py` (empty glycine-offset field skipped; glycine depth printed), and
`ff3_selection/submit_new.sh` (ff30_bio as bEE, ff30_gdepth as dEE). The
ff30_glyhb selection watch (`19174818`) was deleted with the run.

**Automatic validation (user, 2026-10-05; deployed 20:05, backups `.bak_pre_r4val_20261005`; since 10-06 the trainings run on broadwl, so validation follows their `slurm.args` there, with glpG in `popepopg_REMD_mdw2`).** A
converged gate runs `validate_ff.sh` on the epoch-end checkpoint it judged; a run at 13 epochs
unconverged still stops for review. What changed:
* `$P/training/gate_or_continue.sh`: case 0 calls `validate_ff.sh <run> <ff_name> <checkpoint>`
  (refusing a checkpoint that is not an epoch end); a failure of validation fails the gate job.
* `$P/training/validate_ff.sh`: every submission takes the run's `slurm.args`, and the glpG campaign
  follows the partition (broadwl `popepopg_REMD_mdw2`, caslake `popepopg_REMD_mdw3`). Peng arms
  get the script's 12 h (not 36 h) with `runs/<tag>/input` made first. glpG is per candidate:
  `popepopg_REMD_mdw2/seeds/<V>.up` (left untouched) is patched into `$GM/seeds/<V>.<ff>.up` and
  runs as variant `<V>.<ff>`; it stops before anything is released if `$GM/<V>.<ff>` exists. No
  shared seed or replica directory is overwritten or deleted any more.
* `ff30_bio/after_training.sbatch` and `ff30_gdepth/after_training.sbatch` name `ff_3.0_bio` and
  `ff_3.0_gdepth`.
* `/beagle3/.../ff3_benchmark/bench.sbatch`: `BENCH_RESUBMIT` carries the job's own `ExcNodeList`
  (its in-script exclusions name midway2 nodes, which caslake rejects).
* New `/project/trsosnic/yinhan/popepopg_REMD_mdw3/`: `run_remd.py` and `env.sh` md5-identical to
  mdw2's; `remd.sbatch` sources `env_shared.sh`, then `env.sh` (`$P`'s binary and py);
  `submit_remd.sh` has mdw2's block arguments on `--partition=caslake --exclude=midway3-0014
  --mem=56G`.
Per candidate that is 32 Peng arms (14 CPUs, 12 h chunks) and 4 glpG chains (28 CPUs, five 36 h
blocks), ~95k core-hours; ~165 G on `/project` and ~50 G on `/beagle3`. Tested 10-05 in
`/project/trsosnic/yinhan/checks/r4val_20261005/` (gate branches with stubs; `validate_ff.sh` end
to end with `--test-only`, 32 + 4 PASSED in 29 s; base-seed independence; engine parity). Two short
real caslake test jobs were cancelled at the user's correction (findings 10.16); every `ff_r4valtest`
artifact is removed. The first production arm and block are the runtime check: the watch reads them.

`$P` = `/project/trsosnic/yinhan/upside2-md-mdw2`. ff30_gly's cluster run_output holds the local
Mac run's synced steps (marker `run_output/FROM_LOCAL`); the cluster's own initialised run_output is
kept as `run_output.superseded_20261001-111755`. The Mac Studio's local run
(`training/ff30_gly_local`) was stopped by the user (confirmed 2026-10-02 00:25); its steps after
the cluster start (10-01 12:16:47) are discarded.

**Round-4 epoch panels (dt 0.015)** `b00` (49192284), `d00` (49192501) and `b01` (49194315, COMPLETED 10-06 ~10:20; folded 0.400, helix -0.034, findings 1.22) are done; `b00` and `d00` by 04:17 on 10-06: folded 0.465 / 0.471, each start dominates its epoch in helix (table in findings 1.22). `b00`'s local glpG TM4 (Mac `scratchpad/ff3_local_test/runs/79HIS_b00_*`): no glycine flip, TM4 0.98 / 0.92 / 0.91 / 0.92. `d00`'s (finished 04:47): GLY149 flips in one seed of three from block 2, GLY143 transiently in one; TM4 0.95 / 0.90 / 0.92 / 0.92 (findings 1.22). At the user's request (06:48) the newest checkpoints are under the same test: `b01m12` (ff30_bio step 32, margin +0.071) and `d01m09` (ff30_gdepth step 29, margin +0.081), extracted to `checks/r4_epochs/`, Mac runs `79HIS_{b01m12,d01m09}_T080_s{1,2,3}` finished 08:26: TM4 last block 0.87 / 0.87, b01m12 no flips, d01m09 GLY149 in one seed (findings 1.22).

**Selection panel data** (`/project/trsosnic/yinhan/ff3_selection`): controls complete,
`runs/ff21_released/` (49136636) and `runs/ff21_awh/` (49136637 NODE_FAIL on midway2-0085, finished
as 49136943, 17:47); all-atom `aa/` (44 domains); `runs/e00/` (49137240, COMPLETED 19:37);
diagnostics `runs/e00_hb21/`, `runs/fp_e00/` (the ff21-fixedpoint epoch-0 control) and the
single-group resets `runs/e00_{bb21,env21,rot21,sheet21}/` (all COMPLETED 10-02 by 10:55), findings
1.17. Candidates `e01`, ... are submitted by the watch. Round 4's baselines `runs/bio_start/`
(60113478) and `runs/gdepth_start/` (60113190, both caslake, COMPLETED 10-05 17:18) give, with
`panel.py select aa domains runs ff21_released ff21_awh bio_start gdepth_start` (43 domains, 4hwiB01
dropped; folded, then helix / beta / gly_helix / gly_left): ff21_released 0.615, -0.018 / -0.017 /
-0.081 / +0.029; ff21_awh 0.585, -0.021 / -0.019 / -0.081 / -0.065; bio_start 0.602, -0.022 /
-0.018 / -0.050 / +0.003; gdepth_start 0.603, -0.019 / -0.015 / -0.041 / +0.008 (bootstrap SE
0.004-0.005 for helix and beta, 0.043-0.049 for gly_helix, 0.019-0.028 for gly_left). Neither start
is resolved from ff21_released in any class. The select run takes 70 s and 300 MB on a login node.

**midway2-0085 failed twice on 2026-10-01** (NODE_FAIL of panel jobs 49136637 at 16:31 and
49137818 at 22:13), the record that excluded midway2-0037. Panel jobs use `ff3_selection/slurm.args`
and ff30_glyhb's `slurm.args` carries it too; ff30_gly's never did.

**ff30_gly, what must not be forgotten:**
* `upside_input/` is a hardlink copy of ff30_basin's, except `rama.dat`, which is a fresh copy
  of `parameters/common/rama31.dat` (md5 `fc479d45...`, built by `training/build_gly_library.py`,
  log `checks/build_gly_library_20261001.log`). Never edit a hardlinked file in place.
* **The gate no longer releases (changed 2026-10-01 ~13:45, user).** `after_training.sbatch` gates
  with max 13 epochs. Converged: it stops with nothing released and lists the epoch-end
  checkpoints; the one to release is chosen by the selection panel (plan.md Phase 8), then
  `bash ../validate_ff.sh . ff_3.0 run_output/epoch_EE_minibatch_18/checkpoint.pkl` from the run
  dir releases it, backing up the failed basin-offset release, and submits the Peng arms and glpG
  to broadwl. Not converged: trains one more epoch, as before. Backups
  `$P/training/{gate_or_continue.sh,validate_ff.sh}.bak_pre_selection_20261001` and
  `ff30_gly/after_training.sbatch.bak_pre_selection_20261001`; tested on midway2 in a sandbox
  (converged, not converged, and a missing checkpoint, which stops before any release).
  **From midway3 those submissions fail**, and `validate_ff.sh` stops at "not every benchmark arm
  was submitted": the release happens, validation does not. Release from midway2.
* The trainer code in the midway2 tree has had the offset training removed. The superseded files,
  including `verify_rama_basin.py` and `build_rama_from_awh.py`, are in
  `$P/backup/training_pre_glyfix_20261001`. Old run directories keep their own ConDiv copies,
  so the analyses of ff30_basin and ff30_glyprobe still work.
* The ff_3.0 glpG chains were moved aside to
  `popepopg_REMD_mdw2/<V>.ff_3.0_cancelled_20260930` (STOP files inside), so a release cannot
  delete them. `validate_ff.sh` would `rm -rf popepopg_REMD_mdw2/<V>`, and `run_remd.py`
  recreates the directory.

**ff_3.0 validation cancelled by the user at 20:34 on 09-30**, after ~18 h, because ff_3.0 failed
(findings 1.14-1.15): Peng arms 49131950-49131981, glpG chains 49131982-49131985. Checked afterwards
that no benchmark wrapper resubmitted. **Each glpG variant directory has a `STOP` file**
(`popepopg_REMD_mdw2/glpG-RKRK-*/STOP`), so `run_remd.py` will end any chain at once: remove it, and
reset `block_count`, before validating a new force field there. The partial ff_3.0 data stay in
`ff3_benchmark/runs/<prot>_<mode>_ff_3.0/` and `popepopg_REMD_mdw2/<V>/` for comparison. What they
showed by 19:00 (`checks/glpg_tm_windows_20260930_1900.log`, `gly_tm4_series_20260930_1900.log`):
TM4's helical glycines GLY136/143/149 flipped to phi > 0 (up to 0.23-0.64 of frames at T 0.80 in all
four variants, and at T 0.70 in 79ALA_S115T), and TM1 (30-48, no glycine) also declined at T 0.70
(79HIS 0.99 -> 0.89-0.90). In the pre-ff3 campaign GLY143 never flipped, and T 0.70 gave TM4b 0.991
and TM1 1.000 (findings 3.10c). **These chains did not run ff_3.0 as released** (findings 3.11).
The 09-30 02:25 deploy rewrote the live seeds with ff_3.0's map, H-bond and pair tables, but left
the FF1-form ff_3.0's coverage tables in them. The flips belong to that mixture.

The glycine probe (49133133, `training/ff30_glyprobe`) and its successor 49133135 were cancelled by
the user at 20:30 on 09-30, at 15 of 19 steps, once the direction was settled (findings 1.15); the
directory has no `after_training.sbatch`, so nothing can resume or release from it.

Logs and data: Peng `/beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<mode>_ff_3.0_<jobid>.out`,
`runs/<prot>_<mode>_ff_3.0/`; glpG `/project/trsosnic/yinhan/popepopg_REMD_mdw2/logs/remd.<V>.<jobid>.out`,
`popepopg_REMD_mdw2/<V>/`. The superseded ff_3.0 benchmark is in `runs_superseded/ff_3.0_20260930-022510`,
the pre-release glpG seeds are `seeds/*.bak_pre_ff_3.0_20260930-022510`.

**glpG requeues no longer advance `block_count`** (fixed 2026-09-30 13:20). `remd.sbatch` has no
`--no-requeue`, so a NODE_FAIL requeue reran the job under the same id and `run_remd.py` counted it
as a new block. `run_remd.py` now keeps the block's job id in `<V>/block_jobid` and increments only
for a new id (sandbox-tested: fresh, requeued and successor jobs). Backups
`run_remd.py.bak_pre_requeue_20260930`, `submit_remd.sh.bak_pre_requeue_20260930`,
`<V>/block_count.bak_pre_requeue_20260930`. `submit_remd.sh` carries the training chain's node
exclusions (`midway2-0003,[0010-0011],0037,[0342-0345]`). The aborted attempts' output is kept as
short `output_previous_N` groups (100-287 frames), all finite (`checks/glpg_fragments_20260930.log`).

### Finished runs whose files are still used

* **ConDiv runs finished 10-09 to 10-21** (files kept in `$P/training/<run>`; results in findings 1.21, 1.26, 1.30,
  1.34, 1.36, 1.37): ff30_gdepth_dt009 (released ff_3.0_gdepth 10-09; links 49194447/49194449; its 32 Peng arms
  complete, scored 10-21), ff31_depth_all (links 49218340, 49219752; not validated), ff31_depth_data (links 49218341,
  49221772; released ff_3.1_depth_data 10-13; 32 Peng arms complete 10-20), ff32_depth_data (stopped by the user
  08:33 10-18 in step 69, link 49406726 and successor 49410635; steps 0-68 kept), ff21_nobp_dt009 and ff21_ctrl_dt009
  (continued from step 38 on the fixed engine, links 49415222 / 49415223, gates 49418965 / 49430101 CONVERGED at step
  76; mdCATH panels `mdc_nobp_03` 49419028 and `mdc_ctrl_03` 49430102; their 32 Peng arms complete and scored as
  ff_2.1_nobp / ff_2.1_ctrl38, `ff21_nobp_dt009/setup/score/compare_20261021.txt`), and the fixed-engine reruns
  ff31_depth_data_angfix / ff32_depth_data_angfix (built by `checks/angfix_20261018/build_rerun.sh`; links 49415231,
  49415233, 49424034, 49424035; step 76 10:35 / 10:43 10-21).
* **glpG chains stopped by the user 08:33 10-18:** 79HIS ff_2.1_ctrl38 49411323 and ff_2.1_nobp 49411324 (STOP files
  first, directories kept): their step-37 releases are superseded.
* **Analyses of 10-20 and 10-21, all done** (files in `checks/converge_20261020/`, `checks/ff21_released_20261021/`,
  `checks/ff311_20261021/`): TM-score of the side-chain fixed points (`tm_mdcath` 49436371, `p44_dt009` 49436374 /
  49436375, `tm_panel` 49436485, plan.md Phase 26); the epoch-end panels `c9_02`, `c9_03`, `a31_00`-`03`,
  `a32_00`-`03` and the table `panel44_trajectory_final.txt`; Peng scores of ff3.0 / ff3.1 (`peng_trained.txt`); the
  single-file resets `c9_00_*21` (`resets_c9_00.txt`); the ff3.2 rerun's convergence checks; the shared-width
  cross-validation (`cv_pool_*.txt`).
* **midway3: Boltz-2 HDX sampler, Phase 0 first test, finished 10-21** (`diffusion_dynamics.md` Phase 0, where the
  readout is): 61216915 A100 (07:37-07:44, midway3-0294), 61216916 V100 (04:37-04:44, 0279), 61216917 RTX 6000
  (05:17-05:24, 0283), all COMPLETED in about 7 min, each 1 GPU on `gpu`/QOS `gpu`, `pi-trsosnic`. Directory
  `$B = /beagle3/trsosnic/yinhan/boltz_hdx`: `boltz/` at `b1ebfc4`, `.venv` (torch 2.14.1+cu126), `cache/` (weights,
  CCD), `inputs/` (YAMLs, MSAs, references), scripts `install.sh`, `prep.py`, `first_test.sbatch`,
  `check_samples.py` (local copies in `scratchpad/boltz_hdx/`). Logs `$B/logs/boltz_t0_<gpu>_<id>.out`; data
  `$B/runs/boltz_t0_<gpu>_<id>/` (per run `<tag>.log`, `<tag>.gpu_mib`, predictions under `<tag>/`, `<tag>.npz`
  with every sample's N/CA/C). midway2's `squeue` does not see midway3 jobs; check from the midway3 socket.

* **ff3.2 `ff32_depth_data` stopped by the user 08:33 10-18 in step 69** (link 49406726, successor 49410635; plan.md
  Phases 17 and 21): trained with the angular-gradient defect (findings 3.12). Steps 0-68 and panels f32_00-02 are the
  defect side of the comparison with its rerun `ff32_depth_data_angfix`; its readouts to step 68 are in WATCH_STATUS.md
  of 10-18 07:42 and progress.md. The partial step 69 (`run_output/epoch_03_minibatch_12`, 23 of 24 proteins) is left as it stopped.
* **ff21_nobp_dt009 and ff21_ctrl_dt009 at step 38** (10-12 and 10-07): the 38-step comparison (findings 1.21) and its
  audit (plan.md Phase 20); releases `ff_2.1_nobp`, `ff_2.1_ctrl38`; their 79HIS glpG chains stopped 08:33 10-18 with
  STOP files (`popepopg_REMD_mdw2/glpG-RKRK-79HIS.ff_2.1_{nobp,ctrl38}/`; replica files deleted 10-18, seeds and logs kept).

* **TM4 HDX on ff_3.0_gdepth's block 2 COMPLETED 02:04 10-14** (user 10-14 01:19, confirmed directly: "Yes, 79ALA then
  79HIS"; on midway3 caslake because no midway2 environment has mdtraj, pymbar and matplotlib): 79ALA 60668572 (groups 20-38,
  5,723 frames per rung), 79HIS 60668573 (18-35, 5,421), 5 min each. Launcher `checks/hdx_ff30_20261014/hdx_block2.sbatch`
  (`092d3e29...`) with `extract_block.py`, the repo's `martini_hdx_project.py` and `readout.py`; logs
  `checks/hdx_ff30_20261014/hdx_<V>_<jobid>.out`; pre-checks `prechecks_20261014.txt`. **Results**:
  `popepopg_REMD_mdw2/<V>.ff_3.0_gdepth.placedefect_20261017/hdx_sc600x12_ff_3.0_gdepth/` since 10-17 (findings 4.17;
  its `block2/` and `joined/` copies deleted 10-18 with the chain's replica files), `results/` for PS with `<V>_tm4_readout.txt`, `protein_only/results/` for PS_protein); HDX topology from ff_3.0_gdepth
  (rama `common/rama.dat`; md5 `dae228ee...` 79ALA, `f4c0ea21...` 79HIS). Readout in findings 5.5: TM4's core censored at
  every rung (above 5.85 at 298 K), only the helix ends resolve; no ladder round trip in block 2.
* **glpG REMD ff2.1_sc cancelled by the user 10-14 01:11** ("you can cancel ff2.1 job"): `scancel 49299095` in block 2's
  chunk at step 6,930 of 10,494 (block 1 49223635 whole). Run dir `popepopg_REMD_mdw2/glpG-RKRK-79HIS.ff_2.1_sc/` (replica
  files deleted 10-18, counters kept), logs `logs/remd.glpG-RKRK-79HIS.ff_2.1_sc.{49223635,49299095}.out`, seed `seeds/glpG-RKRK-79HIS.ff_2.1_sc.up`.
  All 28 replica files open read-only after the kill: `output_previous_0..29` with finite last frames and a
  100-frame `output` (the 00:34 flush; Upside traps SIGTERM, `main.cpp:1252-1404`, and KillWait is 30 s). Its HDX
  (a spec was written 10-11) was not run. Lesson for stopping any chain at a block's end: `run_remd.py` starts a chunk only
  while remaining - 2100 > 1.25 x the last chunk's wall and reads STOP at that loop's top and before the resubmit
  (lines 190-210), so the last chunk is the one whose remaining R < 2100 + 2.25 x its own wall; take the walls
  from the block in question (block 2's ran 6,745-7,051 s against block 1's ~7,200 s).
* **ff21_ctrl_dt009 (49200579) COMPLETED 10-07 21:25 at its target 38** (`target 38 reached`; the chain cancelled successor 49200622; no gate). The matched control: ff30_bio_dt009 with ff2.1's `rama.dat`, so it starts from ff2.1 exactly (`checks/ctrl_dt009_init_20261006`); `$P/training/ff21_ctrl_dt009`, log `condiv-train_49200579.out`. Steps 0-37, 24 of 24 each, KE/1.5kT at most 1.017. Margin +0.132 at step 14, +0.135 at step 18, +0.125-0.130 over steps 19-31, then +0.134 / +0.139 / +0.147 / +0.152 / +0.155 / **+0.157** at steps 32-37; dhb -0.468, sheet mean 0.174 at step 37. Panel **`c9_00` folded 0.390, helix -0.032** (findings 1.22); TM4 `c9_00` unwound 7, flipped 1 against ff21_released's 9 and 5, not resolved (findings 1.24). `c9_01` (step 37) extracted to `checks/r4_epochs/c9_01` (21:47); **panel `c9_01` (49206076, 22:50-23:38): folded 0.405, helix -0.029** on 34 domains, c9_00 0.390, -0.031 on the same set (findings 1.22); **TM4 `c9_01`** (Mac Studio 04:19-06:15): unwound 10, flipped 6 of 12, against ff21_released's 9 and 5 (p 1.00) and c9_00's 7 and 1 (flipped two-sided p 0.069), not resolved (findings 1.24).
* **ff30_gdepth_si deleted by the user 10-12 00:25** ("since ff3.0 is working, it is the only branch we need (and its expansion to all residues). delete ff_3.0_gdepth_si"; then "clean up to save disk space. no need to keep input files"): the SI-rate twin of ff30_gdepth_dt009 (hb 0.005, dhb 0.0025, sheet 0.0075, rot 0.0125; findings 1.23), stopped at step 72 of 76, before its target and gate, so nothing was released (link 49219651 and successor 49222183 cancelled 00:24). Margin on hold since step 48, +0.1004 at step 70; glycine dL - dR +0.575 after round 3. Deleted (manifest `checks/disk_reclaim_20261011/manifest_ff30_gdepth_si_20261012.txt`, 54.5 G, mostly a half-written step): the run directory and its first-link log, the ds_00-02 extracted ends, panels and TM4 sets (ds_00, ds_01_sc, ds_02_sc) on the cluster and the Mac. Its results stand in findings 1.21, 1.24 and 4.15; the comparison tables are kept. `submit_new.sh` drops it (backup `.bak_pre_gdepth_si_drop_20261012`).
* **The BioEmu-map line was dropped by the user 10-11 23:50** ("the training based on bio emu Rama map is no longer needed ... we decided to go with basin depth training path (as it is expandable to all residues)"): ff30_bio (dt 0.015), ff30_bio_dt009 (b9_EE, released `ff_3.0_bio`) and ff30_bio_si (bs_EE, released `ff_3.0_bio_si`), all finished. Their 57 running validation jobs were cancelled (64 benchmark arms' chunks and the 8 glpG chains, STOP first; `checks/disk_reclaim_20261011/cancelled_bio_jobs_20261011.txt`). **Deleted entirely** (user 10-11 23:40: "Delete them as well. Also no need to leave input file for this case"; manifests in `checks/disk_reclaim_20261011/`): the three run directories and their first-link logs, `parameters/ff_3.0_bio` and `ff_3.0_bio_si` in both trees, the 8 glpG chains with their seeds and logs, the 64 benchmark run directories and logs on `/beagle3`; and the derived products of this line and of the frozen arms: extracted ends `checks/r4_epochs/{b00,b01,b01m12,b9_*,bs_*,bio_start,bz_*,dz_*}`, their panels (`ff3_selection/runs`, `params`, logs, `.submitted`), their TM4 sets on the cluster and the Mac. Freed ~230 G on `/project`, ~40 G on `/beagle3`, ~45 G on the Mac. Kept as the record: the `tm4_compare_cov_*.txt` and `select_last_*.txt` tables; findings 1.21, 1.31, 1.32, 4.15. `submit_new.sh` lists only the live lines since 10-11 23:47 (backup `.bak_pre_bio_drop_20261012`); panel tables now use ff21_released, ff21_awh and gdepth_start as baselines. Not touched: the BioEmu glycine library itself (repo `parameters/common/rama31.dat`, `checks/gly_bioemu_map/`).
* **polygly production stopped by the user 10-11 23:05** ("I don't think poly GLY matters any more as our fix repaired tm4 issue in GLY flip"; plan.md Phase 12): `polygly/prod/STOP` written, successor 49223804 then link 3 49219705 cancelled; every replica checkpointed at 23:05, ~97 ns each (link 1 49204428 to 48.3 ns at 33.8-34.0 ns/day, link 2 49206266 to 94.4-94.7 ns at 32.4-32.6 ns/day); `gmx mindist -pi` minimum 2.31-3.57 nm over both links (`polygly/prod_analysis/`; gmx 2024.4 needs `module load gcc/10.1.0`). All of `polygly/` is kept (213 MB: protein-only trajectories, checkpoints, builds, scripts), so it can resume from `prod.cpt` after deleting STOP.
* **polygly collapse (49186415) COMPLETED 10-07 12:48** (exit 0, 34.4 h, stopped by `-maxh 34.5`): 4 replicas of Ac-(Gly)20-NHMe, amber99sb-ildn / TIP3P, 300 K, reached 19.79-20.00 ns of the 30 ns cap. Minimisation had stopped at Fmax 1380 kJ/mol/nm (atom 124), above its 100 target. Replicas are in `polygly/collapse/repN/`, log `polygly/logs/collapse_49186415.out`; their last frames start production 49204428.

**ff30_basin** (released ff_3.0, 2026-09-28 to 09-30; plan.md Phase 4, findings 1.8-1.14):
* Where: `$P/training/ff30_basin`; logs `condiv-train_<jobid>.out`, gates `ff30b-gate_<jobid>.out`
  and `gate_step{76,95,114}.txt`; checkpoints `run_output/epoch_EE_minibatch_MM`, round libraries
  `run_output/rama_round_0{0..6}.dat` and the round summary `run_output/rama_rounds.txt`; the release
  files in `release_20260930-022510/`. Converged at the step-114 gate (rama p 0.0158, prior-limited
  drift, findings 1.14). Pre-proline scripts and logs: `/project/trsosnic/yinhan/checks/prepro_*.py`,
  `prepro_{leverage,control,left}_20260930.log`.
* Release (gate log): hb [-1.878 -1.872 -1.798], dhb -0.617, bb scale -0.351, sheet mean 0.247;
  glpG round-trip gate 3.6e-15 on a pristine ff_2.1 seed.
* Rewound three times on 2026-09-28, each time recomputing round 1 from the epoch-0 simulations:
  `rewound_20260928_0851/` (MAP step), `rewound_20260928_1225/` (GLY|GLY sheet symmetry, findings
  1.11), `rewound_20260928_1608/` (reduced to 60 maps, findings 1.12-1.13; all-map step-18
  checkpoint `checkpoint_epoch00_mb18_all_maps.pkl`). Pre-change files in
  `$P/training/backup_pre_{basin,ggsheet,reduced}_20260928/`.
* Recovery design that worked unattended: every link queues its `afterany` successor before it
  trains; `--no-requeue`; the successor resumes from the newest checkpoint after deleting the
  half-written step; a worker whose `srun` never started is relaunched (twice at most); three links
  from the same step stop the chain. The gate job itself has no successor.

**ff30** (full-map glycine row) was cancelled 2026-09-28 00:25 at step 223; the run dir and its
last checkpoint `ff30/run_output/epoch_11_minibatch_12` are kept for comparison (the half-written
`epoch_11_minibatch_13`, 40 G, was deleted 10-18); extract with
`$P/training/backup_pre_basin_20260928/extract_ff.py`, since the current one does not read that
layout.

**The FF2 trainer** on midway2 is `training/ConDiv.py` (findings 9v); the old FF1-form trainer and
helpers are in `$P/backup_training_ff1form/`. A full-protocol step waits for its slowest worker,
~21-22 min (5vhg, 150 residues, 1242 s; job 49056799, `worker_test/`), so `STEPS_PER_LINK = 90`
fits a 36 h wall.

**`validate_ff.sh`** (in `$P/training`) extracts a checkpoint through the run's own `expand_param`,
backs up and overwrites `parameters/ff_3.0` in `$P` and in the /beagle3 deployment (md5-verified),
moves the superseded `runs/*_ff_3.0` benchmark directories to `runs_superseded/`, submits the 32 Peng
arms (~7 days), gates the glpG patch on a pristine ff_2.1 seed, patches the 4 live seeds and submits
the 4 REMD chains (~7.5 days). It submits to broadwl, so it runs from midway2.

**`bench_run.py`** lost its type-0 burial override and `rama3.dat` fallback on 2026-09-24 (backup
`bench_run.py.bak_ff1form_20260924`), so it cannot benchmark the old FF1-form ff_3.0.

### Disk: from midway2, `df` on the subdirectory is the only number you get for `/project`

```
df -h /project/trsosnic     -> the FILESET (size, used, free), correct
df -ih /project/trsosnic    -> fileset inodes
df -h /project              -> 6.3P total, the whole device, useless
rcchelp quota               -> home, scratch, project2 group only; no /project or /beagle3 row
mmlsquota                   -> "File system project is not known", the GPFS client is not here
```

`df` on the subdirectory reports the fileset because GPFS `--filesetdf` is on; `df` on the mount
point reports the 6.3 PB device. That distinction is the whole trap.

| fileset | 2026-09-09 | 2026-09-23 | 2026-09-28 | note |
|---|---|---|---|---|
| `/project/trsosnic` | 1514 G free | 445 G free | 953 G free; **965 G free 09-30 12:00** | `training/` is ~108 G, of which `ff30` 54 G; the four glpG variant directories were emptied at the release and are refilling (~150 G each at the end of the last campaign) |
| `/beagle3/trsosnic` | - | 1.4 T free (4.2 T of 5.5 T) | 1.4 T free | badly degraded for small-file writes on 09-28, see §8 |
| `/project2/trsosnic` (group) | - | 1.45 T of a 1.49 T soft quota, **97%** | - | nothing in these campaigns writes there |
| midway3 home | - | 28.6 of 30 G | 21 G | |

**2026-10-11 17:20: `/project/trsosnic` 88 G free of 3.9 T (98%), inodes 32%.** Ours: training
313 G, ff3_selection 422 G (84 panel tags of ~11 G each, of which the 176 `.up` trajectories are
~11 G and the npz `panel.py select` reads ~70 MB), popepopg_REMD 638 G (pre-fix baselines, kept),
NP-1AO6 491 G (deletion awaiting the user), popepopg_REMD_mdw2 199 G, checks 113 G, polygly 0.2 G;
other users the rest. Growth: each of the 12 running glpG chains adds ~0.22 GB/h (seed 5.9 G,
then ff_3.0_gdepth 13 G at 31 h, ff_3.0_bio_si 12 G at 29 h, ff_3.0_bio 9 G at 14 h), ~63 G a day,
so the 88 G lasts ~1.4 days, less each new panel (~11 G) and each training step in flight (~37 G
of replica files, deleted at the step's end). Their remaining blocks need ~410 G.
**Reclaimed 17:26-17:29 at the user's word** ("delete the trajectory file but try to keep the setup
scripts so I can regenerate them"), manifests and log in `checks/disk_reclaim_20261011/`: the 39
finished panels' `.up` trajectories, 6,864 files, 397 G (npz, `.up.output`, `params/`, `domains/` and
`panel.py` kept; `select` reproduces b9_03's table byte for byte); NP-1AO6's `prod_ff3` (471 G),
`prod_ff3_armM_20260910-084115` (14 G), `prod_counterion` (3 G) and `tscan` (1 G), keeping every top-level
script (`build_np_ff3.py`/`.sbatch`, `build_tscan.py`, `np_hybrid.py`, `run_np_prod.py`, `np_prod.sbatch`,
`submit_np.sh`, ...); the four `glpG-*.ff_3.0_cancelled_20260930` directories (40 G, replica files and
counters only). The NP directories also held each orientation's built `system.MARTINI.pdb` and
`hybrid_mapping.h5` (regenerated by the build scripts) and two small analysis scripts,
`prod_counterion/onset.py` and `onset2.py` (1.6 and 1.3 KB), which had no copy and are gone.
**17:42: 1,041 G free (74%).**

**2026-10-18 11:57-11:59: reclaimed 1.1 T at the user's word** (591 G free before, 1.7 T after, 59%; ours ~0.75 T of
2.3 T used). Manifests, scripts and log in `checks/disk_reclaim_20261018/` (`build_manifests.sh`,
`delete_manifests.sh`, `delete.log`):
1. Half-written steps of stopped trainings, no checkpoint (86 G, whole directories): `ff30` `epoch_11_minibatch_13`,
   `ff30_gly` `epoch_03_minibatch_08` and `run_output.superseded_20261001-111755/epoch_00_minibatch_00`, `ff30_gdepth`
   `epoch_01_minibatch_18`, `ff30_glyhb` `epoch_03_minibatch_07`. ff3.2's partial step 69 was left.
2. Replica `.up` files of the defective glpG runs, and their HDX `joined/` and `block2/` copies (821 G): the 8
   `.placedefect_20261017` chains (findings 4.17), `popepopg_REMD/*_pre_cbfix` and `*_rigidbug` (August), the 79HIS
   ff_2.1_nobp, ctrl38 and ff_2.1_sc chains, `armtest_{M,R}` and their failed-venv copies, `tm4_local/runs` and
   `runs_cov`. Kept: HDX `results/`, `outputs/`, `inputs/`, `pdb/`, `protein_only/`; STOP files, counters, logs, seeds,
   tables.
3. The `.up` trajectories of 14 finished panels (228 G): `runs/{ctrl38,nb_01,fa_01-03,fd_00-03,f32_00-02}`,
   `runs_mdcath/{mdc_nobp,mdc_ctrl38}`. npz, `.up.output`, `params/` kept; `panel.py select` reproduces
   `select_last_f32_02_login1_20261017.txt` byte for byte (`select_f32_02_after.txt`). CA-RMSD (`panel_compare.py`)
   can no longer be recomputed for these tags.

Not touched: `popepopg_REMD`'s four variant directories (355 G, `hdx/` on the §1b list), `gly_peptides` (§1b),
`glpG_DDM_micelle_REMD` (73 G, §1b, waiting on the user), `implicit_79HIS_*` (70 G), `popepopg_glpG` (9 G).
What grows from here: the 8 corrected glpG chains (~38 G each over blocks 2-5), four trainings' steps in
flight (~38 G each), each panel ~11 G (mdCATH ~47 G), each block-2 HDX ~28 G.

---

## 1b. Data and decisions that outlive their jobs

The campaign sections these came from were removed on 2026-09-28 (they are in git history). What
must not be forgotten:

* **ff_3.0 is ff30_gdepth_dt009's release** (user, 10-09 10:40): `$P/parameters/ff_3.0_gdepth`
  and `/beagle3/trsosnic/yinhan/upside2-md/parameters/ff_3.0_gdepth` (from
  `ff30_gdepth_dt009/release_20261009-085341`, = `checks/r4_epochs/d9_03`), copied to the repo's
  local `parameters/ff_3.0/` (md5 verified). The cluster names stay as they are.
* **Data that exists only on the cluster, do not delete:**
  * AWH glycine surfaces, `/project/trsosnic/yinhan/gly_peptides/` (rep1, `awh_amber99sb-ildn_rep2/`,
    `evidence_diffusion_bug/`).
  * `popepopg_REMD_mdw2/run_remd.py` and `NP-1AO6/run_np_prod.py` exist only there, with no git history.
  * `popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`, the pre-ff3 TM baseline.
  * HDX: keep `popepopg_REMD/<V>/hdx/` (pre-fix). `hdx_postfix/` and the earlier ff3 results at
    `popepopg_REMD_mdw2/<V>/hdx_10k/`, both recorded here before, are on neither `/project` nor
    `/beagle3` (searched to depth 5, 2026-09-30); `hdx_results/` holds only the DDM-campaign figures
    and MBAR files.
* **Deletion decisions waiting on the user:** (NP-1AO6's run data were deleted 10-11, scripts kept;
  §1 "Disk")
  * midway3 DDM data `glpG_DDM_micelle_REMD/`, `glpG_DDM_REMD/`.
  * rockfish post-fix glpG data `/scratch4/rherna21/ywang268/upside2-md-rf/popepopg_REMD/<variant>/`
    (~72 GB, needed only for an independent-seed cross-check), before scratch is reaped.
* **Safe to delete:** `/beagle3/trsosnic/yinhan/ff3_benchmark/sigmoid_discarded/` (4 G) and
  `armM_discarded/`.
* **Open, unverified, recorded so they are not lost:**
  * The HDX topology uses ff_2.1 tables on an ff_3.0 trajectory.
  * `NP-1AO6/np_footprint.py` has CB-offset and minimum-image bugs, with no record of a fix.
  * The one-shot upgrade line in `NP-1AO6/np_prod.sbatch` is still to delete.
  * An NP orientation at T = 0.70-0.75 was never run.
  * A latent stride-3/4 issue was noted at `upside_config.py:815` `_input_phi`.
  * The HDX resubmit block uses `HDX_WORK=.../hdx` without `HDX_N=28` and would overwrite the pre-fix
    baseline: fix it before any HDX rerun.
* **NDRD library files must never be copied to the cluster** (licence).
* **Do not overwrite the cluster's `py/` or `training/` from the local repo during these campaigns.**
  The repo's `training/` is now five files: `extract_ff.py`, `convergence_gate.py` and `rama_basin.py`
  became `ConDiv.py extract` / `gate` and functions, and `train_chain.sbatch` runs the gate itself
  instead of submitting `after_training.sbatch` -> `gate_or_continue.sh` (plan.md Phase 9). The
  running ff30_glyhb chain needs the midway2 tree's old layout (`after_training.sbatch`,
  `gate_or_continue.sh`, `extract_ff.py`, `check_step.py`, `rama_basin.py` beside its run copy), so
  the new files go to the cluster only for a new run. The midway2 tree also still runs on the
  campaign files the repo moved to `scratchpad/redistribution_cleanup_20261001/`: ff30_gly's release uses `$P/training/{validate_ff.sh,patch_glpg.py}`, the glpG HDX
  jobs (`popepopg_REMD/hdx_*.sbatch`) call `$R/py/martini_remd_concat.py`, and
  `NP-1AO6/build_np_ff3.py` imports `martini_inject_coverage`. On 2026-10-02 the repo's
  `training/env.sh` and `train_chain.sbatch` became site-neutral (no midway2 module loads, no
  `--account`); the midway2 tree and every run directory keep their own copies with both, which is
  what the running chains need.

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
* **A 54-step link times out when steps run long** (2026-10-07). `STEPS_PER_LINK=54` assumes ~36
  min a step; ff30_gdepth_dt009's took 40-115 min, so its first link (49194447) hit the 36 h wall
  in step 35. The successor resumes from the newest `checkpoint.pkl`; the half-written step is
  lost, and the run waits in the queue as at any link end. Expect the same of ff30_gdepth_si
  (32-63 min a step).
* **Size wall-time requests to the work; long requests do not backfill on broadwl** (2026-10-01).
  Two 1-node panel jobs asking 24 h were estimated to start 1.5-2 days out; cut to 4 h with
  `scontrol update JobId=<id> TimeLimit=04:00:00`, they started within the hour. 36 h benchmark arms
  were estimated 2.5 days out, so the benchmark now runs 12 h jobs with ~10.5 h chunks
  (`bench.sbatch`, `bench_run.py` `CHUNK_SEC`; backups `*.bak_pre_12h_20261001`). A queued job keeps
  the batch script it was submitted with, so a changed resubmit line needs a fresh submission.
* **A script that submits another job must unset `SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU` and
  `SLURM_MEM_PER_GPU`** after its `#SBATCH` block (or match the child's memory-request type), or
  every `srun` of the child dies with `... are mutually exclusive`. Findings 10.1.
* **Slurm snapshots the batch script at submission.** A requeue reruns it with its original
  arguments, and an edit does not reach queued jobs. Chain scripts carry `--no-requeue` (check
  `scontrol show job <id> | grep Requeue`) and resolve the newest checkpoint at run time; after an
  edit, cancel and resubmit the queued successors (`scontrol write batch_script <id>` shows what a
  queued job will run). Findings 10.2.
* **Cancel a chain's PENDING successor before its running link**, then check `squeue`; the other
  order starts the successor.
* **A link killed by SIGBUS on several nodes at once is a transient `/project` outage**: 0-byte
  `*.output_worker` files, no traceback, and a successor that dies in its first second. Resume;
  exclude nothing. Findings 10.3.
* **One login node's `/project` can hang while the other's works** (2026-10-07 10:48). On
  midway2-login2 every lookup under `/project/trsosnic/yinhan` sat in D state (`wchan`
  `lookup_slow`), while `/beagle3`, `/software`, home and login1's `/project` answered and the jobs
  kept writing. `timeout` cannot kill a D-state process, and `ps` on it hangs too, so each probe
  stacks another stuck command. Test with `timeout 15 ls <dir>` once, then route `/project` work
  through login1 over the socket: `ssh ... midway2.rcc.uchicago.edu 'ssh -o BatchMode=yes
  midway2-login1 bash -s' < script.sh` (host-based, no Duo), and rsync with `-e "ssh -o
  BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu ssh -o BatchMode=yes"` and
  host `midway2-login1`. login1's `/project` was slow too (a 70 s panel select took 4 m 44 s).
* **A failed worker fails the step** (`ConDiv.py` raises): step 81 of ff30 once ran on 5 of 24
  proteins after `srun: Job credential expired`, and a partial step is a different objective. That
  credential error at launch is a race, not a node fault (24 steps issued together lost 4 on healthy
  nodes); the trainer relaunches a worker that never started, up to twice. **When both relaunches
  land on one node the link ends:** ff30_bio_fz step 5 (10-07 14:40, link 49201673) lost 7 of 24
  launches; 2cwr and 4exo failed twice more on midway2-0088 ("Invalid job credential"). The step
  was discarded, steps 0-4 stand, and the successor resumes from step 4 but re-enters the queue.
  It is one event on 0088, so no exclusion yet; a second on the same node is the record to exclude
  it. **The second came 10-09 03:42:** ff21_ctrl_fz step 11 (link 49205852) lost 8 launches, and
  2r2y, 3jyz and 4qbo failed all three on 0088 ("Invalid job credential"), and the link FAILED at
  04:11. 0088 was added to the law at 04:12 (user, §0d).
* **10-09 09:07-09:17: the controller-link fault again (below), on excluded nodes and two new
  ones.** Other users' broadwl jobs NODE_FAILed at 09:07:42 (on midway2-[0033-0035,...,0342-0345]),
  09:12:42 (0080) and 09:17:42 (0085), all excluded nodes. ff30_bio_dt009's link 49206580 and
  ff30_gdepth_si's link 49205173 died at 09:17:42 on disjoint 16-node sets. Afterwards one node of
  each shows MIXED+NOT_RESPONDING: midway2-0103 (slurmd started 10-01 09:10:33) in 49206580 and
  midway2-0116 (09-23 12:11:13) in 49205173; the other 30 are back in use. Five-minute spacing to
  the second is the 10-02/03 signature. Neither node had an earlier record; both were
  added to the law at the user's word (10:30, §0d). Each link lost its running step; the successors resume from
  steps 64 and 21.
* **Node exclusions are per record of failures; since 10-08 the list lives only in the law
  (§0d).** midway2-0027 was added then for two NODE_FAILs as a batch host (00:12, 08:32); a third,
  ff30_gdepth_si's link at 11:07, had started on it before the law.
  midway2-[0010-0011] (three NODE_FAILs in one
  campaign) and midway2-0037 (two links in three hours on 2026-09-26) are excluded in `slurm.args`;
  two later link failures named no node and excluded nothing. On 2026-10-02 20:40:36-37 midway2-0003 and
  0033-0035 failed in the same second (link 49141995, lambda 49144354), with ~20 broadwl nodes not
  responding at 20:50: a group event rather than one node's record, so nothing was added. On 2026-10-03
  04:50:36 link 49142816 died NODE_FAIL on midway2-0323 or 0366 (not reallocated, both back in service at
  05:20); one failure each, nothing added. **Superseded 2026-10-03 13:20: the NODE_FAILs are one
  fault on a fixed set of nodes.** At 10-02 20:40:37, 10-03 04:50:36 and 13:00:38 broadwl jobs of
  several users died in the same second, and 32 of the 33 such jobs held at least one of the 9 of 212
  broadwl nodes whose slurmd started at 2026-10-01 18:10:54: midway2-[0010-0011,0033-0035,0037,0060,0080,0085]
  (the exception was midway2-0003, DOWN for its own reasons). No node rebooted or restarted slurmd at
  the events; the workers' MaxRSS was 1.1-1.5 GB. Our exclusions already carry 0010-0011, 0037 and 0085
  but not 0033-0035, 0060 and 0080, which each sat in a failed job at all three events. The first three
  events were 8 h 10 min apart to 3 s, but the period is **disproved**: at 17:25:37-38 jobs on 0003 and
  0010 died (ours and wdenault's), and on 0011 at 17:30:37, all three then shown `mixed*` (not
  responding), and none came at the 21:10 the period predicted. The fixed-node part holds, the timing does not. At 17:48 our job on 0011 was writing
  frames at 12.5 tu/s while Slurm showed 0011 `mixed*`: the nodes compute, and it is their link to the
  controller that drops, so the controller declares NODE_FAIL and kills what runs there. **Excluded 10-03 14:25 (user-approved):** both
  `ff30_glyhb/slurm.args` and `ff3_selection/slurm.args` now read `--partition=broadwl
  --exclude=midway2-[0003,0010-0011,0033-0035,0037,0060,0080,0085,0342-0345]` (backups
  `.bak_pre_nf_20261003`, format checked with `sbatch --test-only`), and the pending 49145374,
  49145390 and 49145792 carry the same `ExcNodeList`. The lambda benchmark got the 9 nodes and 0003 at
  10-03 17:45 (user-approved) as a `#SBATCH --exclude` line in `ff3_benchmark/bench.sbatch` (backup
  `.bak_pre_nf_20261003`) rather than in `BENCH_RESUBMIT`: a running job builds its resubmit command
  from Slurm's spooled copy of the script, so only a directive in the file reaches the next chunk.
  Checked with `sbatch --test-only -w`: an excluded node is refused, an allowed one accepted. The four
  running lambda jobs keep their nodes and requeue themselves on a NODE_FAIL (`Requeue=1`).
* **`sbatch --test-only`'s start estimate is no guide to the real wait**: it predicted 13:36 on both
  clusters, then the midway3 submission sat PENDING (Resources) while midway2 started at once.
* **midway3 gives a new job no start estimate for its first 20 minutes**: `bf_min_age_reserve=1200`
  in its `SchedulerParameters`, so backfill neither reserves for it nor fills `squeue --start`
  (N/A, reason Priority) until it has pended 1200 s. Not a sign the job is unschedulable.
* Training steps ran 300-450 s slower than the engine accounts for (2026-09-29), a different slowest
  worker each time; most likely `/project` I/O.
* **`/beagle3` was badly degraded for small-file writes on 2026-09-28**: 200 one-line files took 51 s
  from midway3 and 503 s from midway2, against 0.06-0.27 s on `/project` and home; a `pip install
  torch==2.6.0` into the shared venv took ~2 h for that reason.
* **Do not tell the clusters apart by `/software/modules/init/bash`**: it exists on midway3 too. The
  cluster `training/env.sh` tests whether the tree venv's interpreter exists (verified 2026-09-28:
  midway2 gets the tree venv, midway3 the `/beagle3` venv, both torch 2.6.0+cpu, numpy 1.23.5, scipy
  1.13.1, tables 3.8.0). Never pipe `source env.sh`: a pipe runs it in a subshell and the
  environment is lost.
* **Delete a half-written minibatch directory before resuming**: a stale `divergence.pkl` is reused.
* **Hand-check success-only end-of-chain branches**, which never run until the end.
* **Pre-create `runs/<tag>/input` before a mass submission** on `/beagle3` (a makedirs race); a
  partial aggregate mean is not a result.
* **Midway3 home quota is 30 G** (21 G used on 2026-09-28, 28.6 G on 09-23). Jobs can fail oddly if
  home fills.
* **`sacct` on midway3 can hold zombie RUNNING rows** (53233848 and 53233852 from 2026-08-11 still
  show under `-S <today>`, while `sacct -j` says `COMPLETED`); `squeue` is the authority.
* **Exclude midway3-0014.**

Login nodes and shells:
* **`/tmp` on a login node is per node, and reconnects land on different ones.** A `nohup` job
  launched from midway2-login2 writes a `/tmp/<log>` that login1 cannot see, so a later check reports
  the job gone. Put anything you will read again on `/beagle3` or `/project`, and prefer `sbatch` over
  `nohup` for work that must outlive an ssh session. Do not run scripts from `/tmp` either.
* **A stray module in the working directory shadows the stdlib.** Another user's `/tmp/inspect.py`,
  or an `inspect.py` left in a scratch directory, broke `import numpy` with a circular-import error.
  Name scratch probes something that is not a stdlib module.
* **`pgrep -f <name>` matches your own shell command.** `ssh host 'pgrep -f bm_final && echo running'`
  reports "running" because the remote `bash -c` line contains the string, so a dead job looks alive
  (two false reports on 2026-09-12). Match the interpreter (`ps -u $USER -o cmd | grep python3`),
  check the output file, or use a Slurm job id.
* **The shared venv is Python 3.9 / NumPy 1.23**: `str | Path` needs `from __future__ import
  annotations`.
* **midway2's site GROMACS is unusable**: use `/project/trsosnic/yinhan/gmx2024`
  (`gmx_build/build_gmx.sbatch`). The HDX pipeline runs only on midway3 (`.venv_el8_py311_bak`, with
  pymbar and matplotlib); a renamed venv's `bin/activate` still hardcodes its old `VIRTUAL_ENV`.

Simulations:
* **NP dt hard limit: 0.001.** dt=0.005 caused backbone blow-ups during unfolding (large-amplitude
  spring instability at t>250, proven by A/B). Never raise above 0.001.
* **Do not transfer thresholds between NP and glpG.** A CN_MAX borrowed from NP false-positived on a
  healthy glpG chunk (2.52 Å vs healthy max 2.659 Å) and cost a 6 h block.
* **glpG NaN propagation via REMD exchange.** A single blow-up in one replica spreads to every
  replica via exchange within ~60 steps (IEEE 754: `NaN < 0.f` = false). The NaN cascade fix
  (`!isfinite(lboltz_diff)`) and the per-chunk rollback driver address this.
* **A green exit code means nothing** for a self-submitting REMD job. Check the log for
  DESTROYED/ROLLBACK counts and verify physical observables.
* **glpG: patch the LIVE seeds** (they carry `inner_steps = 4` on `/input/brownian`); prove the method
  by an ff_2.1 round trip on a pristine seed. A release must remove each ~150 GB variant directory
  before resubmitting it.
* **glpG REMD: skip `output_previous_0`** (the frozen seed block). Do not gate a resubmission on
  `check_seeds_current.py` (§5).
* **Identify a `.up`'s force field by least-squares match of its baked tables to `sidechain.h5`**
  (scale 1.000000), not by hash. `martini.h5` comes from ff_2.1 by design.
* **Live dynamics are shown by internal Kabsch CA-RMSD growth** (frozen = 0.000 A), never by reading
  `current_stage`; read only `output_previous_*` of a live `.up`, with `HDF5_USE_FILE_LOCKING=FALSE`.
* **Runs branched from a shared checkpoint do not test reproducibility**; branch at step 0.
  Independently trained force fields are different Hamiltonians and are never pooled.
* **NP:** the box must exceed the molecule's maximum extent plus the 12 A cutoff; `NP_DT = 0.001`; the
  Rg in `np.<jobid>.out` is not minimum-imaged (measure adsorption by backbone contacts with the MPA
  shell, atoms 3750:3950); a restart frame must pass a stricter test than detection; `NP_TEMP` is
  Peng's FF2 calibration and must be redone for a new force field.
* **HDX figures:** keep continuous off-scale excursions and `Y_LIMITS=(-20,30)`, and never render
  censored amides as carets or bounds (user-directed). Cluster HDX scripts drift from the repo;
  re-upload before trusting results.
* The glpG detergent campaign (closed 2026-08-13, model retired 2026-09-09) used a constant `--seed`,
  so a rollback re-ran the identical failing chunk; the driver now takes a per-chunk seed.
