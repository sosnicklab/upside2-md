# Remote jobs on midway2/midway3: status and handbook

**Current state (2026-10-08 10:48; the Mac Studio keeps running the watch and the local TM4 queue
(user, 10-08); another computer stands by until the user moves the jobs: "Handoff" below, and the
live `WATCH_STATUS.md` on the cluster).** The glpG TM4
repair training, round 4 (plan.md Phase 11), runs on midway2 broadwl at dt 0.009, every run from
ff2.1. The user's question is which run moves glpG TM4 toward stable. It is judged only by the
hybrid glpG TM4 test of each checkpoint against its own start (findings 10.19):

| run (tags) | glycine library | `hb`, `dhb`, `sheet` | target | TM4 start |
|---|---|---|---|---|
| ff30_bio_dt009 (`b9_EE`) | BioEmu-fitted, frozen | trained | 76, gate | bio_start (`ff21_bioT1_6`) |
| ff30_gdepth_dt009 (`d9_EE`) | NDRD, glycine depth trained | trained | 76, gate | gdepth_start |
| ff21_ctrl_dt009 (`c9_EE`) | ff2.1's (control) | trained | 38, no gate | ff21_released |
| ff30_bio_fz (`bz_EE`) | BioEmu-fitted, frozen | frozen | 76, gate | bio_start |
| ff30_gdepth_fz (`dz_EE`) | NDRD, glycine depth trained | frozen | 76, gate | gdepth_start |
| ff21_ctrl_fz (`cz_EE`) | ff2.1's (control) | frozen | 38, no gate | ff21_released |
| ff30_bio_si (`bs_EE`) | BioEmu-fitted, frozen | trained at the SI's rate | 76, gate | bio_start |
| ff30_gdepth_si (`ds_EE`) | NDRD, glycine depth trained | trained at the SI's rate | 76, gate | gdepth_start |

"Trained" is the port's rate (hb 0.01). "The SI's rate" (user, 10-07; findings 1.23) is every group
but rot at half that (hb 0.005, as FF2 was trained); rot is 0.0125 in every run. The SI-rate runs
started 10-07 16:34 (bio) and 17:20 (gdepth); both runs' rates are confirmed. ff30_gdepth_fz and
ff21_ctrl_fz were cancelled before they started (10-07 13:12). ff21_ctrl_fz was resubmitted at 13:35
and ff30_gdepth_fz at 21:46 (user: run every frozen arm, though freezing H-bond and sheet departs
from the group's workflow). ff21_ctrl_fz against c9 is the cleanest test of "H-bond drift
destabilizes TM4", with no glycine change. ff21_ctrl_dt009 (c9) finished at its target at 21:25.

Job ids, states and readouts are in §1. ff30_bio_dt009 has been on margin hold since step 14,
ff30_gdepth_dt009 since step 33 (+0.097, 10-07 20:18) and ff30_bio_si since step 26 (+0.097,
10-08 06:06); the user was told of each, and all train on. Links that ended and resume by chain design: ff30_bio_dt009's first at step 53 (normal) and its
second on a node failure in step 64 (10-08 08:32; successor queued), ff30_gdepth_dt009's at its wall time in step 35 (successor queued),
ff30_bio_fz's on a node fault at step 5 (resumed 10-07 21:46), ff30_bio_si's first at step 53 (normal, 10-08 20:10; successor queued), ff21_ctrl_fz's on a node failure
in step 2 (10-08 00:12; successor queued), and ff30_gdepth_si's on a node failure in step 21 (10-08
11:07, midway2-0027 again; successor queued). ff30_gdepth_fz's first link runs from 21:34 10-08.

**TM4 on fixed inputs so far** (findings 1.24): bio_start 5 of 12 unwound, b9_00 7 of 24, b9_01 6,
bs_00 6, bs_01 5 (2 flipped, against b9_01's 5; p 0.37), bz_00 8 of 24, bz_01 6 (2 flipped; TM4 0.846
against b9_01's 0.833, so the freeze does not stop the epoch-1 lean), b9_02 8 (4 flipped; TM4 0.830, so
epoch 2 keeps b9_01's lean), ds_00 9 (TM4 0.775, the lowest mean; d9_00 5, p 0.21), ff21_released 9, c9_00 16 of 24 (flipped 1 in seeds 1-12, 8 in 13-24, so its 12-seed lean was chance), c9_01 10, fp_e00 5, gdepth_start 7, d9_00 5,
d9_01 6. No pair is resolved. bz_00's first half flipped none of 12 (one-sided p 0.047 against
bio_start's and b9_00's 4); its seeds 13-24 flipped 4 of 12, so that was chance. At 24 seeds each,
bz_00 against b9_00 is 8 against 7 unwound and 4 against 5 flipped (p 1.00); only TM1, in an added
test, is higher in bz_00 (0.891 against 0.813, p 0.003). Until 10-07 11:22, `patch_glpg.py` left the
retired FF1-form ff_3.0's coverage tables in every glpG input (findings 3.11), so every earlier TM4
count is invalid. The queue runs locally ("Handoff" below; §1 "Local glpG TM4 test"). Panels:
every epoch end so far loses folding against its start (findings 1.22), the control as much as
any (c9_00 0.390, c9_01 0.405 against ff2.1's 0.615). b9_02 0.448, d9_01 0.457, d9_02 0.439 and bs_00 0.457
(findings 1.21): each is dominated by gdepth_start in helix, and folding falls about 0.015 an epoch
after epoch 0's 0.125. The frozen and SI-rate runs fold less than their port-rate twins at the same
step (bz_00 0.450, bs_01 0.430, ds_00 0.412 against 0.478, 0.465, 0.477); d9_00 dominates ds_00 in
helix and b9_01 dominates bs_01 in helical glycine.

Also on midway2:
- polygly production 49204428 (`polygly/prod.sbatch`, self-chaining; plan.md Phase 12), running
  since 10-07 22:41 at 33.7 ns/day per replica, so 1 us takes ~30 days. Its collapse (49186415)
  finished 10-07 12:48.
- The BP validation deck for Tobin (§0c) is not yet sent.

CPU work runs on midway2 broadwl unless the user names midway3 caslake. midway3's login node may be
used to move or read files on `/project` and `/beagle3`. Since 10-07 10:48 midway2-login2's
`/project` hangs; route `/project` work through login1 (§8).

Written so a fresh session can pick up cold. §1 is live job state; §0, §2 and §5-§8 are the
handbook.

---

## Resume here, from any computer (2026-10-08)

### "continue jobs": what a new session does, in order

Claude's memory does not travel between computers; everything needed is in this file, plan.md and
findings.md (pulled from git) and on the cluster.
1. **Read** this file's header table and §1 (live state; each watch pass leaves its readouts
   there), then:
   - plan.md Phase 11 and Known Errors;
   - findings.md 1.24 (TM4 on fixed inputs), 3.11 (the coverage defect) and 10.16-10.19 (the
     rules: no test jobs, same data and workflow, re-read the queue, TM4 only in the hybrid).
2. **Connect** to midway2 ("Connect" below).
   - Check the socket first. If `scratchpad/mdw2_master.exp` and `~/.bin/ssh_mdw3` exist on this
     computer, run the expect script once (one Duo push); otherwise ask the user to open the master.
   - While login2's `/project` hangs, run `/project` work through midway2-login1 (§8). Test with one
     `ls` first.
3. **Take stock, read only.** Run `squeue -u yinhanw`. For every run in the header table, find its
   newest `run_output/epoch_*` and check its link log for `WORKER_FAIL`, `Traceback`, `STOPPED`.
   Tell the user what started, finished or failed since the last pass in §1, and the current
   margins.
4. **Owed checks, once each** (as of 10-09 21:27; `WATCH_STATUS.md` "Owed checks" is the current list):
   - Successors' first lines: ff30_bio_fz 49206071 `step 59 of 76, resuming from
     run_output/epoch_03_minibatch_01/checkpoint.pkl`; ff21_ctrl_fz 49214161 `step 11 of 38, resuming
     from run_output/epoch_00_minibatch_10/checkpoint.pkl`; ff30_bio_dt009 49215780 `step 65 of 76,
     resuming from run_output/epoch_03_minibatch_07/checkpoint.pkl`; ff30_gdepth_si 49216034 `step 22
     of 76, resuming from run_output/epoch_01_minibatch_02/checkpoint.pkl`; polygly 49206266's first
     lines. ff30_gdepth_fz's link 49206070 ends after step 53 (~04:00 10-10); 49214121 follows.
   - Panels at 176 npz (watch step 6): d9_03 49216387, bs_02 49216513, dz_01 49218034.
   - ff_3.0_gdepth validation: the 32 benchmark arms (R in their second chunk, logs clean at 21:27) and the 4 glpG chains
     49216140-49216143 (PD) once they start.
   - **ff30_bio_dt009's gate** (`ff30bio9-gate`) after step 75, ~7 h after 49215780 starts: read
     `gate_step<N>.txt` and the gate log (watch step 6). Converged means validation was queued:
     32 Peng arms and 4 glpG chains for `ff_3.0_bio`. Not converged means one more epoch.
5. **Find the owner** ("Handoff" below): read `checks/r4_epochs/tm4_local/WATCH_STATUS.md` on the
   cluster. While the Mac Studio owns the jobs and its last pass is under 2 h old, stand by: report
   status read only, start nothing. A new session on the Mac Studio itself after a restart follows
   "Mac Studio session restart" in "Handoff" instead.
6. **Only on the user's "take over":** set up the local TM4 test ("Setting up the next computer"),
   run the TM4 queue from `WATCH_STATUS.md` (a set left running is rerun unless `runs_cov/` has
   it), start the watch ("The watch"), and rewrite `WATCH_STATUS.md` with the new owner.
7. **Open questions for the user.** Nothing is blocking. Queue waits on broadwl run to days; midway3
   caslake is an option only on the user's word, after checking its current estimate. midway2-0103
   and 0116 joined the law 10-09 10:30 (user; §8, 09:17). **Asked 10-09 21:30:** 25 of the 32
   ff_3.0_gdepth benchmark chunks run on the pre-10:30 list, without 0103 and 0116, because
   `bench.sbatch:25` resubmits each chunk with its own job's `ExcNodeList`, so the gap carries to
   every later chunk; 49218194 (proteinB native) runs on midway2-0103. The fix would be
   `bench.sbatch` taking the law file on midway2; nothing changed without the user. Offered to
   the user (10-09 05:50): seeds 13-24 of bz_02 and b9_02, whose 12-seed primary test resolved at
   p 0.046 with an interval reaching -0.010 (findings 1.25). midway2-0088 joined the law 10-09 04:12 (user; §0d). The user made DSSP alpha-helix the
   primary TM4 test (10-09). The 10-07 offer of 12 more ff21_released seeds, to test c9_00's flip difference (1
   against 5 of 12), is withdrawn: c9_00's own seeds 13-24 flipped 8, which removes it (findings
   1.24). The 24-seed frozen comparisons were approved 10-08 09:10, and the midway2 exclusion law
   (§0d, midway2-0027 added) was built at the user's word 10-08 09:40.

### Handoff: two computers, one owner (2026-10-08)

**The Mac Studio keeps the jobs** (user, 10-08 09:15, replacing the planned 10:00 move). It runs
the watch (session cron `45cc9e4a`, hourly) and the local TM4 driver (`scratchpad/ff3_local_test`).
Another computer, such as the user's work computer, may open a session at any time. It must know
what the Mac Studio is doing and must not do it a second time: two watches would both submit
panels, and two TM4 drivers would run the same sets.

**One owner at a time.** The owner runs the watch and the TM4 driver; every other computer only
reads. **Owner now: the Mac Studio.**

**Mac Studio session restart.** The watch cron is session-only and ends with its session; the Mac
Studio stays the owner. The TM4 driver runs under `nohup` and the midway2 socket master is its own
process, so both survive unless the computer restarts. The new session on the Mac Studio:
1. Reads §1 and `WATCH_STATUS.md`, and checks the socket (the expect script once if it is dead).
2. Runs one watch pass at once (watch steps 1-7) for the steps since the last pass.
3. Starts the watch again: CronCreate hourly at `23 * * * *` with the prompt in "The watch" below.
4. Rewrites `WATCH_STATUS.md` (last pass; any restart note removed).
The last restart (10-09 21:20, user: an upgrade) ended cron `cad6ca59`; the new session ran the
pass at 21:27 and started cron `45cc9e4a` at 21:28.

**Live status on the cluster:** `/project/trsosnic/yinhan/checks/r4_epochs/tm4_local/WATCH_STATUS.md`.
The owner rewrites it at every watch pass (watch step 7): owner, time of the last pass, the TM4 set
running and the queue behind it, and what is owed. This file (`remote_jobs.md`) reaches another
computer only through the user's commits, so it can be hours old there; `WATCH_STATUS.md` is the
current word.

**What another computer does on "continue jobs"** (after "Resume here" steps 1-4):
1. `cat` `WATCH_STATUS.md` over the midway2 socket.
2. **If the owner is the Mac Studio and its last pass is less than 2 h old: stand by.** Do not
   start a watch or a TM4 driver; do not submit, extract or patch anything. Reporting status read
   only and answering the user's questions are fine. To be ready to take over, it may set up
   `scratchpad/ff3_local_test` and run the parity check ("Setting up the next computer", steps 1-4).
3. **Take over only on the user's word** ("take over"). If the last pass is more than 2 h old, tell
   the user and ask whether the Mac Studio has stopped; do not assume it. Taking over:
   1. If the Mac Studio's session is still open, the user has it run its stop steps (below) first.
   2. Set up as in "Setting up the next computer".
   3. Write the queue from `WATCH_STATUS.md` into `tm4_queue.txt`. A set listed as running is
      rerun unless `tm4_local/runs_cov/` has all its seeds and its table (a rerun reproduces it
      frame for frame), and its lines go first. Then start the driver (setup step 5).
   4. Start the watch ("Starting the watch") and rewrite `WATCH_STATUS.md` with the new owner.

**The Mac Studio's stop steps** (when the user says "hand off" or "stop"):
1. Move the remaining lines of `scratchpad/ff3_local_test/tm4_queue.txt` (and `tm4_queue.next`
   while it exists) into `WATCH_STATUS.md`'s queue and empty the file. The driver runs under
   `nohup`, outlives the session, and with the file empty exits after its current set.
2. Copy every finished set and its table to `tm4_local/runs_cov/` (watch step 5).
3. Record the set still running and when it ends. If the session is open when it ends, compare and
   copy it as usual; if not, the next owner reruns it.
4. Update §1, this section and `WATCH_STATUS.md` (owner: none), then CronDelete `45cc9e4a`.

**24 seeds for the frozen comparisons** (user, 10-08 09:10). bz_00 against b9_00 and cz_00 against
c9_00 run seeds 13-24 as well as 1-12. A 24-seed set is two queue lines, `<tag>` (seeds 1-12) and
`<tag> 13 24` (log `run_<tag>_s13-24.log`); `tm4_queue.sh` reads both forms (10-08 09:19, md5
`eac953d0...`). `tm4_compare.py` (10-09, with the DSSP primary test; md5 in "Setting up the next computer") reads a set's seeds 1-12 or
1-24 from `runs_cov/` and `runs/` together and gives each set its own n in the Fisher table; it
reproduces the earlier 12-seed tables exactly, and a mock 24-seed set (b9_00 plus c9_00 seeds) gave
the expected 11 unwound of 24. A set whose seeds 13-24 have started is read at 12 seeds with a third
argument, `tm4_compare.py <tag> <ref> 12` (it reproduces b9_00's table with and without it). bz_00
was read so at 12:05; the 24-seed comparison follows once both halves of bz_00 and b9_00 are done. Copy seeds 13-24 and their table to
`runs_cov/` like any set. The patched inputs of b9_00 and c9_00 are the files seeds 1-12 ran on
(every `/input` dataset equal; `ref_pos` holds NaN placeholders, so compare it NaN-aware).

**The TM4 queue** (one 12-seed half at a time, T 0.80, 4000 tu, `run_glpg.sh 10 <first> <last>
<tag>`; about 2 h each on the Mac Studio's cores). Each tag is patched from `checks/r4_epochs/<tag>`
(ff21_released from `parameters/ff_2.1`) with the fixed `patch_glpg.py`. Done on fixed inputs:
bio_start (`ff21_bioT1_6`), b9_00, bz_00 and c9_00 (24 seeds each; `tm4_compare_cov_bz_00_24_vs_b9_00.txt`,
`tm4_compare_cov_b9_00_24.txt`, `tm4_compare_cov_c9_00_24.txt`), b9_01, b9_02, bs_00, bs_01, bs_02, bz_01, ds_00, ff21_released,
c9_01, fp_e00, gdepth_start, d9_00, d9_01, d9_02, d9_03, bz_02, dz_00, dz_01. State at 21:27 10-09:
1. Nothing runs: the driver finished dz_01 at 20:45 and exited with `tm4_queue.txt` empty; the
   next epoch end restarts it (watch step 4).
2. New epoch ends as they come. Each is extracted with `$P/training/extract_ff.py`, its panel
   goes in through `submit_new.sh`, and it is patched and added to `tm4_queue.txt` (watch step 4).
   A half-trained (`_01`) end goes behind the 24-seed lines and ahead of the other epoch
   ends. `cz_00` runs as `cz_00` and `cz_00 13 24`. `submit_new.sh` skips a checkpoint written
   less than 5 min ago, so a pass that finds one that new runs it again later. Expected (09:47 10-09):
   `cz_00` (step 18, against ff21_released and c9_00 at 24 seeds) ~4 h after successor 49214161
   starts; `b9_03` (step 75, the end of training) ~7 h after successor 49215780 starts; `ds_01`
   (step 37, against gdepth_start and d9_01) ~10-13 h after successor 49216034 starts. Both of the
   last two links ended by NODE_FAIL at 09:17 (§8). `dz_02` (step 56) ~2 h after successor
   49214121 starts (the running link 49206070 ends after step 53, ~04:00 10-10).

Reading: each checkpoint against its own start (header table). **The primary test is TM4's
secondary structure** (user, 10-09; findings 1.25): each seed's DSSP alpha-helix fraction of 135-151
over the last block, set against reference by a two-sided Mann-Whitney, resolved at p < 0.05, with a
bootstrap 95% interval on the difference (`tm4_compare.py`, "## primary test"). bz, cz, dz and the
SI-rate runs are also compared with their port-rate twin at the same step (bz_00 and bs_00 with
b9_00, cz_00 with c9_00, dz_00 and ds_00 with d9_00). The direction is read across epoch ends.
- **What it can resolve.** Seeds spread with SD 0.28 in last-block alpha. A true difference of 0.10,
  0.15, 0.20 or 0.30 is detected with power 18%, 31%, 47% or 73% at 12 seeds per set, 34%, 58%, 80%
  or 96% at 24, and 56%, 86%, 97% or 100% at 48. The differences seen so far are 0.01-0.15.
- **Secondary tests**, kept so earlier results stay comparable: the dihedral-box counts
  pre-registered 10-06, unwound (last-block TM4 < 0.90) and flipped (a TM4 glycine at phi > 0 in
  > 0.25 of the last block), by Fisher's exact test; and the DSSP fraction in helix of any kind.
  The box accepts pi-bulges and frayed residues, and DSSP alpha runs 0.10-0.22 below it in every
  set. c9_00's two halves, on one input, flipped 1 and 8 of 12.
- **24 seeds** run for the frozen comparisons (above). Any other extension is the user's decision,
  given the power figures.

**Setting up the next computer**
1. Pull the repo (the user commits these files first; Claude never touches git); build `.venv` and
   `obj/upside`.
2. Open the midway2 socket ("Connect" below).
3. Make `scratchpad/ff3_local_test` from `checks/r4_epochs/tm4_local` (not in git):
   - `scripts/` (fixed `patch_glpg.py` md5 `705891194b57b2139f1e817b80b0135a`, `run_glpg.sh`,
     `tm4_local.py` md5 `2701e5d0...`, `tm4_compare.py` md5 `ce50f0c6...`, `tm4_queue.sh` md5 `eac953d0...`) and
     `events_vs_tm4.py` (md5 `ef9ce96c...`).
   - `seed/glpG-RKRK-79HIS.live.up` into `seeds/` (the watch prompt names `seeds/`).
   - `runs_cov/` (~14 GB): the finished sets `tm4_compare.py` reads, seeds 13-24 included.
   - `checks/r4_epochs/<tag>` into `ff/<tag>` for every tag still to run, and `bio_start` into
     `ff/ff21_bioT1_6` for the parity run. Patch each:
     `python3 scripts/patch_glpg.py --ff ff/<tag> --seed seeds/glpG-RKRK-79HIS.live.up --out
     patched/79HIS_<tag>.up`. Its output lists `hbond_coverage changed by ...` (~19.8 / ~14.2).
   - `run_glpg.sh` and `tm4_queue.sh` set
     `L=/Users/yinhan/Documents/upside2-md/scratchpad/ff3_local_test`; change it if the repo lives
     elsewhere. `tm4_queue.sh` uses macOS `sed -i ''` (GNU: `sed -i`). `run_glpg.sh` runs the 12
     seeds in parallel, one core each, so a set needs 12 free cores to take ~2 h.
4. **Engine parity**, once per new computer. Run the patched `ff21_bioT1_6` with `run_glpg.sh`'s
   flags for 200 tu (`--duration 200`, `--seed 1`, same integrator flags) and compare its frame
   lines, from `time` on, with the first 21 of `runs_cov/79HIS_ff21_bioT1_6_T080_s1.log` (t 0-200):
   `diff <(grep ^step NEW.log | sed 's/^step [0-9]* \/ [0-9]* //') <(grep ^step
   runs_cov/79HIS_ff21_bioT1_6_T080_s1.log | head -n 21 | sed 's/^step [0-9]* \/ [0-9]* //')`. The
   Mac Studio matched exactly (10-07 16:29). A mismatch means its sets are not comparable.
5. Write the queue's remaining lines (`<tag>` or `<tag> 13 24`) into `tm4_queue.txt` and start the driver from
   `scratchpad/ff3_local_test`: `nohup caffeinate -is bash scripts/tm4_queue.sh >> tm4_queue.log
   2>&1 &` (on Linux, without `caffeinate -is`).
6. Start the watch ("Starting the watch" below). The probes `~/watch_probe/stock.sh` (newest steps,
   fail patterns) and `~/watch_probe/steps.sh <run>:<step> ...` (check_step.py and the KE scan) on
   midway2 do watch steps 2-3.

### Where things are

**What travels with git and what does not.** Tracked: `plan.md`, this file, `findings.md`,
`progress.md`, `training/`, `src/`, `py/`, `up.md` and `parameters/common/rama31.dat` (the
BioEmu-fitted library). Commit here, pull there. Not in git (`.gitignore` has `*scratchpad*`): the
`scratchpad/` tools, including `mdw2_master.exp`, `mdw3_master.exp` and the local TM4 test. Their
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
  - `tm1_hbmem_20261009`: the MacBook Pro's 10-09 glpG TM1 work (findings 1.27 and 4.6-4.9, plan.md
    Phase 13): scripts, tables, the membrane H-bond term runs' logs and positions, the replay, and
    a README of what stayed on the MacBook Pro.

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
after 7 days. The current one runs on the Mac Studio, hourly at :23 (cron `45cc9e4a`, created
10-09 21:28 after the upgrade restart; it expires 10-16, so a session still open then starts
another). It runs until the user moves the jobs ("Handoff" above); another computer starts its own
only after the Mac Studio's has ended. Two watches would both submit panels. Each pass records its
readouts in §1, so the newest §1 tells a new watch where the last one stopped.

**Starting the watch on a computer** (set up as in the Handoff): CronCreate hourly at an off-minute
(e.g. `23 * * * *`) with the prompt below.

**Round-4 watch prompt:**
"Round-4 remote job watch, the only monitor (remote_jobs.md "The watch"). Runs on midway2 broadwl,
dt 0.009: ff30_bio_dt009 (b9_EE), ff30_gdepth_dt009 (d9_EE), ff30_bio_fz (bz_EE) and
ff30_gdepth_fz (dz_EE) (H-bond and sheet frozen), ff30_bio_si (bs_EE) and ff30_gdepth_si (ds_EE)
(SI learning rates), target 76 with the gate; the controls ff21_ctrl_dt009 (c9_EE, finished) and
ff21_ctrl_fz (cz_EE, frozen), target 38, no gate. Job ids are in §1. Also polygly production
(polygly_prod, self-chaining): at its first link's end read each replica's prod.log Performance
line, and check `gmx mindist -pi` stays above 1.0 nm. Follow "Round-4 watch" steps 1-7 exactly;
read §1 first. Obey:
- CLAUDE.md: read-only git, no guards, no physics changes, nothing heavy on a login node;
- findings 10.16 (no test jobs), 10.17 (same data, same workflow), 10.18 (re-read the queue before
  any decision on job state) and 10.19 (TM4 is judged only by the hybrid glpG test).
Never cancel, resubmit or change a job without the user's approval.
For an SI-rate run's first finished step, confirm from its link log that the solver rates are hb
0.005, dhb 0.0025, sheet 0.0075 and rot 0.0125, and report it once.
TM4: run the queue in remote_jobs.md "Handoff", then each new epoch end, half-trained first. On
this computer the queue runs through scripts/tm4_queue.sh and tm4_queue.txt (watch step 4); the
seed directory is seeds/.
- One 12-seed set at a time, on inputs patched with the fixed patch_glpg.py (md5 70589119...).
- Compare each checkpoint with its own start, and bz / cz / dz / bs / ds with their port-rate twin
  at the same step, by the primary DSSP test (findings 1.25; it replaced Fisher's exact test on
  10-09 at the user's word).
- Scan finished seeds for total-potential jumps above 3000 between frames (events_vs_tm4.py).
- Copy each finished set and its table to checks/r4_epochs/tm4_local/runs_cov/.
Notify only for: failures, a TM4 or panel result, a gate, polygly ending, a new margin hold in a
trained-H-bond run (below +0.10), an SI-rate run starting, or a dead socket.
Otherwise stay quiet."

**Round-4 watch.** `<run>` is any of the seven running or queued run directories under `$P/training`.
`submit_new.sh` maps them to their tags.
1. **Socket.** `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu`.
   If it fails, run `scratchpad/mdw2_master.exp` once, unless §1 records that the previous pass
   already tried. If it is still down, record that in §1, notify once and stop the pass. Never loop
   on connection attempts.
2. **Queue and logs.** Run `squeue -u yinhanw` and `bash
   /project/trsosnic/yinhan/ff3_selection/submit_new.sh`, which submits each new epoch end's panel
   without asking (user, 10-08 11:05). Check the newest link log of each run for
   `WORKER_FAIL`, `Traceback`, `STOPPED`, `never started`.
3. **Every step finished since the last pass.**
   - Run `cd $P/training && source <run>/env.sh && python3 check_step.py <run>
     <epoch_xx_minibatch_yy>` under `ulimit -v 4000000`.
   - KE scan of every protein: `grep -m1 ^avg_kinetic_energy/1.5kT run_output/<step>/*.run.*.up.output`.
     Flag a maximum above 1.2 (findings 1.21).
   - Record the shared margin E_other - E_alpha (start +0.192), dhb, sheet mean, rot rms, step time,
     and helical / left glycine alpha_L free/restrained. For the gdepth runs also record dL - dR
     (start +0.554) and `run_output/rama_rounds.txt`.
   - **A trained-H-bond run whose margin falls below +0.10:** notify the user once and hold. The run
     keeps training; nothing is cancelled without the user.
4. **TM4**, in `scratchpad/ff3_local_test` with the repo `.venv` and `source.sh`.
   - For a new `epoch_EE_minibatch_18`, extract it: `$P/training/extract_ff.py
     <run>/run_output/epoch_EE_minibatch_18/checkpoint.pkl /project/trsosnic/yinhan/checks/r4_epochs/<tag>`.
   - rsync it to `ff/<tag>` and run `patch_glpg.py --ff ff/<tag> --seed seeds/glpG-RKRK-79HIS.live.up
     --out patched/79HIS_<tag>.up`.
   - Add `<tag>` to `tm4_queue.txt` in the order "The TM4 queue" gives (a half-trained end behind
     the 24-seed lines, ahead of other epoch ends; `cz_00` as `cz_00` and `cz_00 13 24`). If
     `tm4_queue.log` ends in `queue empty`, restart the driver:
     `nohup caffeinate -is bash scripts/tm4_queue.sh >> tm4_queue.log 2>&1 &`. It runs one line at
     a time: `<tag>` as `run_glpg.sh 10 1 12 <tag>` logged to `run_<tag>.log`, `<tag> 13 24` as
     seeds 13-24 logged to `run_<tag>_s13-24.log`.
5. **A set whose log says `set finished`.**
   - From the repo root, `python3 scratchpad/ff3_local_test/scripts/tm4_compare.py <tag> <ref> >
     scratchpad/ff3_local_test/tm4_compare_cov_<tag>.txt`, with `<ref>` the tag's start (header
     table); a second table against the twin or the epoch before is
     `tm4_compare_cov_<tag>_vs_<ref>.txt`, and a 24-seed table `tm4_compare_cov_<tag>_24.txt`. It
     refuses a set with an unfinished log or seeds other than 1-12 or 1-24; a first half read while
     its seeds 13-24 run takes a third argument, 12. Its output holds:
     - per-seed `avg_kinetic_energy/1.5kT` (~1.0) and `tm4_local.py`'s tables;
     - TM4's DSSP alpha-helix and any-helix fractions (135-151) per seed and block, and alpha per
       residue;
     - the primary test, DSSP alpha-helix by a two-sided Mann-Whitney with a bootstrap interval
       (findings 1.25);
     - the secondary tests: Fisher's exact test on unwound and flipped seeds, with the counts that
       would resolve, and one-sided Mann-Whitney tests on the dihedral TM4, TM1 and any-helix;
     - the energy jumps (`events_vs_tm4.py`).
     It reproduces the MacBook Pro's b9_01 table exactly. Read the direction across epoch ends.
   - Record the result in findings.md, copy the runs (seeds 13-24 too) and the tables to
     `tm4_local/runs_cov/`, and notify.
6. **A panel or a gate.**
   - A panel whose `ff3_selection/runs/<tag>/` holds 176 npz: run `panel.py select aa domains runs
     ff21_released ff21_awh bio_start gdepth_start <tags>` with the run's env (70 s, 300 MB) and
     notify the table.
   - A gate: read `gate_step<N>.txt` and the gate log. Converged means `validate_ff.sh` ran inside
     it, so confirm 32 Peng arms and 4 glpG chains for the candidate were queued, and read their
     first logs.
7. **Update §1 briefly, and rewrite the cluster's `WATCH_STATUS.md`** ("Handoff"): owner and
   computer, time of this pass, the TM4 line running and its expected end, the queue in order, the
   newest step of each run with its margin, and anything owed or awaiting the user.

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
`--partition=broadwl --exclude=midway2-[0003,0010-0011,0027,0033-0035,0037,0060,0080,0085,0088,0103,0116,0342-0345]`.
Each node's record of failures is in `README.md` beside it; §8 has the history. midway2-0088 was
added 10-09 04:12, and midway2-0103 and 0116 10-09 10:30, each at the user's word (backups
`*.bak_pre_0088_20261009`, `*.bak_pre_0103_20261009`). For the second, `--test-only` accepted the
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

Snapshot **2026-10-09 21:27 CDT, verified live against `squeue` on midway2**, the Mac Studio's first watch pass after its 21:20 restart (the socket landed on midway2-login1, whose `/project` answers; login2's hung from 10:48, §8) (no jobs on midway3). Finished and cancelled
rows are deleted; their lessons are in §8.

| JobID / where | what | state | next action |
|---|---|---|---|
| **49215780** (midway2) | **ff30_bio_dt009**: ff30_bio restarted at dt 0.009, otherwise identical (trainer `$P/training/ConDiv.py`, which differs from ff30_bio's run copy only in dt; initial force field byte-identical, `checks/dt009_init_20261006`); `$P/training/ff30_bio_dt009`, target 76, 54 steps per link, gate `ff30bio9-gate` up to 13 epochs (converged: validated as `ff_3.0_bio` on broadwl); first-link log `$P/training/condiv-train_49194446.out` (submitted from `training/`), later links' logs in the run dir | First link 49194446 COMPLETED 10-07 16:02 after its 54 steps (0-53, `epoch_02_minibatch_15`); link 2, 49194448 (from 00:12 10-08; `step 54 of 76, resuming from run_output/epoch_02_minibatch_15/checkpoint.pkl`), **NODE_FAIL 10-08 08:32** in step 64 (`epoch_03_minibatch_07` half-written): its batch host midway2-0027 failed, as it did for ff21_ctrl_fz's link at 00:12 (0027 is MIXED again at 08:45, and is ff30_gdepth_si's batch host; excluded by the law from 09:40, §0d, and every pending job carries it). Steps 54-63 done. Link 3, 49206580 (from 08:23 10-09, 16 nodes; `step 64 of 76, resuming from run_output/epoch_03_minibatch_06/checkpoint.pkl`, `this link runs 12 steps`, owed check done), finished step 64 (08:57, margin +0.007, 24 of 24, KE/1.5kT at most 1.008) and **NODE_FAILed 10-09 09:17:42 in step 65** (`epoch_03_minibatch_08` half-written), in the same second as ff30_gdepth_si's link and another user's job on midway2-0085; its node midway2-0103 shows NOT_RESPONDING afterwards (§8). **Successor 49215780 PD (Priority, estimate 10-10 22:05 at 11:47)**, carrying the law, resumes from step 64 by chain design with 11 steps left (owed: its first line `step 65 of 76, resuming from run_output/epoch_03_minibatch_07/checkpoint.pkl`). Steps 29-63 min each, except 106 min for step 46, which relaunched five workers that never started, as step 45 did six; 24 of 24 returned every step; KE/1.5kT at most 1.017 over every replica of every step (all 1,104 protein-steps to step 45, findings 1.21; at most 1.015 at steps 46-53). **HOLD since step 14** (user told). Margin +0.192 start, +0.096 step 14, +0.072 step 24, **+0.050 at step 37** (`b9_01`), +0.036 / +0.034 / +0.032 at steps 43-45, +0.027 / +0.022 / +0.018 / +0.014 / +0.014 / +0.015 / +0.014 / +0.015 at steps 46-53, +0.017 / +0.018 / +0.020 / +0.016 / +0.011 / +0.005 / +0.002 / +0.003 / +0.001 / +0.002 at steps 54-63 (01:15-07:46, the last link-2 step; KE/1.5kT at most 1.014): E_alpha -1.961 to -1.828, E_other -1.769 to -1.816; dhb -0.406 to -0.492, sheet mean 0.161 to 0.199. Helical GLY aL free 0.02-0.11 per step, restrained <= 0.002. Panels: `b9_00` folded 0.478, helix -0.030; **`b9_01` 0.465, -0.031** (bio_start 0.602, -0.019; findings 1.21); **`b9_02` 0.448, -0.026** on 34 domains (49207052, COMPLETED 21:28), dominated by gdepth_start in helix (findings 1.21). **TM4 `b9_01` on fixed inputs** (finished 15:14): unwound 6, flipped 5 of 12, against bio_start's 5 and 4, Fisher p 1.00 for both; last-block TM4 mean 0.833 against 0.899 (findings 1.24; runs in `tm4_local/runs_cov/`). **TM4 `b9_02`** (Mac Studio 19:49-21:45): unwound 8, flipped 4 of 12, against bio_start's 5 and 4 (p 0.41, 1.00) and b9_01's 6 and 5 (p 0.68, 1.00); last-block TM4 0.830 against 0.899 and 0.833, so epoch 2 keeps the epoch-1 lean (findings 1.24) | the watch every step; `b9_03` (step 75): extract, panel, patch, TM4; the gate at step 76 (`ff30bio9-gate`) ~7 h after 49215780 starts |
| **49216108-49216143** (midway2) | **ff30_gdepth_dt009**: ff30_gdepth restarted at dt 0.009, otherwise identical (trainer `ff30_gdepth_dt009/trainer/`, dt the only change; rama_round_00 and initial force field byte-identical); `$P/training/ff30_gdepth_dt009`, target 76, 54 steps per link, gate `ff30gdep9-gate` (converged: `ff_3.0_gdepth`); first-link log `$P/training/condiv-train_49194447.out` | **Target 76 reached: link 49194449 COMPLETED 10-09 08:49** after step 75 (`epoch_03_minibatch_18`: margin +0.082, dhb -0.551, sheet mean 0.228; depth round 4 offsets aR -0.1068 aL +0.5169, dL - dR +0.624; steps 70-75 24 of 24, KE/1.5kT at most 1.013); the chain cancelled successor 49206658 at 08:49. **Gate 49216104 CONVERGED** (08:54; `ff30_gdepth_dt009/gate_step76.txt`, log `ff30gdep9-gate_49216104.out`): steps 58-76, every group p > 0.005 (lowest dhb 0.045, then bbenve 0.79). It released `epoch_03_minibatch_18` as **ff_3.0_gdepth** to `$P/parameters/ff_3.0_gdepth` and `/beagle3/trsosnic/yinhan/upside2-md/parameters/ff_3.0_gdepth` (md5 verified; `release_20261009-085341`, rama.dat = `rama_round_04.dat`), byte-identical to `checks/r4_epochs/d9_03`. `validate_ff.sh` submitted the 32 benchmark arms 49216108-49216139 (all R at 11:47; the first at 40k of 2.7 M tu with sane T, Rg and potential; logs `/beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<kind>_ff_3.0_gdepth_<jobid>.out`; at 21:27 all 32 R, 31 in a second chunk 49218083-49218237 resubmitted 19:34-21:14 and ubiquitin de novo 49216137 still in its first, 0.22-1.31 M tu done, no error in any log; 25 of the chunks carry the pre-10:30 exclusion list without 0103 and 0116, and 49218194 runs on 0103, see "Resume here" step 7) and the 4 glpG chains 49216140-49216143 (79HIS and 79ALA, each also S115T; PD, estimates 10-10 08:52-22:05 at 11:47; logs `popepopg_REMD_mdw2/logs/remd.<variant>.ff_3.0_gdepth.<jobid>.out`), all carrying the law. Before: **First link 49194447 TIMEOUT 10-07 21:46** at its 36 h wall, in step 35 (`epoch_01_minibatch_16` half-written): its steps took 40-115 min, not the ~36 min the 54-step link assumed (§8). Steps 0-34 done (`epoch_01_minibatch_15` 21:13); workers relaunched after a failed launch at step 16 (2nmu) and step 31 (3tjy), 24 of 24 returned every step; KE/1.5kT at most 1.016 over every replica of all 840 protein-steps. **Successor 49194449 R since 00:22 10-08** (18 nodes; log `ff30_gdepth_dt009/condiv-train_49194449.out`: `step 35 of 76, resuming from run_output/epoch_01_minibatch_15/checkpoint.pkl`, `this link runs 41 steps`; at the 36-70 min of steps 57-71 it reaches step 75 ~09:00-10:15 10-09, before its 12:22 wall; if it times out, successor 49206658 finishes the run); its successor 49206658 queued. Margin +0.131 at step 13, **+0.141 at step 18** (epoch-0 end), +0.140 at step 24, then falling: +0.138 / +0.130 / +0.121 / +0.114 / +0.110 / +0.108 / +0.105 / +0.101 / +0.097 / +0.096 at steps 25-34, +0.097 / +0.096 / +0.097 / +0.101 / +0.106 / +0.109 / +0.111 / +0.108 / +0.104 / +0.099 / +0.091 / +0.084 / +0.078 / +0.077 / +0.075 / +0.075 / +0.074 / +0.073 / +0.073 / +0.070 / +0.066 / +0.066 / +0.067 / +0.072 / +0.075 / +0.078 / +0.075 / +0.074 / +0.074 / +0.072 / +0.073 / +0.073 / +0.080 / +0.087 / +0.091 / +0.091 / +0.089 / +0.089 / +0.085 / **+0.082** at steps 35-74 (00:57-08:16 10-09, steps 38-56 on round 2's library and 57-74 on round 3's; steps 57-74 took 36-70 min, step 67 relaunched 1pz4 once; step 60 relaunched 1wjg once, 24 of 24 returned), 30-54 min each; 24 of 24, KE/1.5kT at most 1.015; E_alpha -1.854, E_other -1.775 at step 47; dhb -0.488, sheet mean 0.222, rot rms 0.107). **HOLD since step 33** (user told 20:30); it trains on. Depth round 1 (`rama_rounds.txt`): free-native gap aR -0.002 aL +0.004, dL - dR +0.567; round 2 (epoch-1 end): gap aR +0.000 aL +0.007, offsets aR -0.0965 aL +0.4808, dL - dR +0.577; round 3 (epoch-2 end): gap aR -0.001 aL +0.013, offsets aR -0.1002 aL +0.5017, dL - dR +0.602. Panel **`d9_00` 0.477, helix -0.024** (gdepth_start 0.603, -0.018, dominates it in helix); **`d9_01` 0.457, -0.025** on 34 domains (49207053, COMPLETED 21:35), dominated by gdepth_start in helix (findings 1.21); **TM4 `d9_00` on fixed inputs** (Mac Studio 22:25-00:26): unwound 5, flipped 2 of 12, against gdepth_start's 7 and 2, Fisher p 0.68 and 1.00, not resolved; last-block TM4 mean 0.900 against 0.849 (findings 1.24). TM4 **`d9_01`** (Mac Studio 06:15-08:12): unwound 6, flipped 4, against gdepth_start's 7 and 2 and d9_00's 5 and 2, not resolved (findings 1.24); DSSP alpha 0.801 at d9_00 and 0.615 at d9_01, against 0.657 (findings 1.25) | the watch: benchmark arms and glpG chains; **`d9_03`** (step 75, the released ff_3.0_gdepth): extracted 09:50, panel 49216387 (09:47, carries the law; estimate 10-10 22:05), patched (its glpG patch matches the gate's own to the printed digits); **TM4 `d9_03`** (Mac Studio 10:46-12:43), primary test: DSSP alpha 0.726 against gdepth_start's 0.657 (+0.069 [-0.165, +0.300], p 0.64) and d9_02's 0.765 (p 1.00), not resolved; 5 unwound, 4 flipped, no energy jump (findings 1.26). d9_02's panel and TM4 are in findings 1.21 and 1.25. Choosing ff_3.0 stays the user's |
| **49206071** (midway2) | **ff30_bio_fz**: ff30_bio_dt009 with hb, dhb, hbg and sheet at learning rate 0 (user, 10-07; plan.md Phase 11 revised decision), trainer copy `ff30_bio_fz/trainer/ConDiv.py` (md5 `a97e95a6...`; differs from `$P/training/ConDiv.py` only in that rate line, its comment and three docstring items); initial force field byte-identical to ff30_bio_dt009's (`checks/fz_init_20261007`); target 76, gate `ff30biofz-gate` up to 13 epochs (converged: validated as `ff_3.0_bio_fz`); first-link log `$P/training/condiv-train_49201673.out` | **First link FAILED 10-07 14:40** (exit 1, 4 h 35 min): at step 5, 2cwr and 4exo could not start on midway2-0088 (`srun: Invalid job credential`, a node-side Slurm fault; 22 of 24 workers had finished), and after two relaunches each the step raised `2 of 24 workers failed` (`condiv-train_49201673.out`). Steps 0-4 done (11:56-14:03), 24 of 24 each, KE/1.5kT at most 1.013; hb, dhb, sheet unchanged through step 4 (margin +0.1919). Link 2, 49202419 (from 21:46 10-07, `step 5 of 76, resuming from run_output/epoch_00_minibatch_04/checkpoint.pkl`), **COMPLETED 10-09 02:55** after its 54 steps (5-58, `epoch_03_minibatch_01`; 28-53 min each), 24 of 24 every step, KE/1.5kT at most 1.017, margin +0.1919 (frozen). **Successor 49206071 PD (Priority, estimate 10-10 07:57 at 11:47)**; it resumes from step 58 by chain design with 17 steps left (owed: its first line `step 59 of 76, resuming from run_output/epoch_03_minibatch_01/checkpoint.pkl`). **`bz_01`** (step 37, 13:58) extracted: hbond.h5, sheet and rama.dat md5-identical to bz_00's (so bio_start's); panel 49210456 COMPLETED 23:14: **`bz_01` 0.461, helix -0.028** (38 domains; b9_01 0.465, -0.031; bz_00 0.450), neither it nor b9_01 dominating, gdepth_start dominating it in helix (findings 1.21); patched, **TM4 `bz_01`** (Mac Studio 17:53-19:49): unwound 6, flipped 2 of 12, against bio_start's 5 and 4 (p 1.00, 0.64) and b9_01's 6 and 5 (p 1.00, 0.37); last-block TM4 0.846 against 0.899 and 0.833, so the freeze does not stop the epoch-1 lean; s12 ends at 0.22, the lowest on fixed inputs (findings 1.24). **`bz_00`** (step 18) extracted: its hbond.h5, sheet and rama.dat are md5-identical to bio_start's (`checks/r4_epochs`); panel 49207340 COMPLETED 21:47: **`bz_00` 0.450, helix -0.030** (36 domains; b9_00 0.478, -0.030), not dominated by b9_00, dominated by gdepth_start in helix (findings 1.21); **TM4 `bz_00`, 24 seeds** (Mac Studio 10:08-14:00): seeds 1-12 unwound 3, flipped 0 (one-sided p 0.047 against bio_start's and b9_00's 4 flips), seeds 13-24 unwound 5, flipped 4; at 24 seeds 8 and 4 against bio_start's 5 and 4 of 12, p 0.45 and 0.24, not resolved (findings 1.24). Against b9_00 at 24 seeds (16:05): 8 against 7 unwound, 4 against 5 flipped, p 1.00; TM1 higher in bz_00 (0.891 against 0.813, an added test, p 0.003). **Freeze check passed** (12:05; `checks/fz_init_20261007/bz_step00`): hbond.h5 and sheet md5-identical to ff2.1's, hb / dhb / sheet mean unchanged (margin +0.1919); sidechain.h5, environment.h5, bb_env.dat moved; rama.dat the frozen BioEmu library | the watch every step; **`bz_02`** (step 56, 01:42): extracted (hbond.h5, sheet and rama.dat md5-identical to bz_01's), panel 49214736 COMPLETED 09:01: **`bz_02` 0.424, helix -0.031** (42 domains), dominated by gdepth_start in helix; against its twin b9_02 0.424 / 0.448 on 40 domains, neither dominating (findings 1.21), patched; **TM4 `bz_02`** (Mac Studio 03:35-05:33), primary test: DSSP alpha 0.856 against its twin b9_02's 0.694 (+0.162 [-0.010, +0.330], two-sided p 0.046), **resolved, bz_02 higher**, the first resolved primary test and only just (the interval reaches below zero; one of 32 tables read); against bio_start's 0.769, p 0.30, not resolved; secondary counts 4 unwound, 3 flipped (findings 1.25) |
| **49205039** (midway2) | **ff30_bio_si**: ff30_bio_dt009 with every group but rot at the SI's learning rate (user, 10-07; plan.md Phase 11; findings 1.23): `alpha_scale` 0.25, rot base doubled so rot stays 0.0125, hence hb 0.005, dhb 0.0025, sheet 0.0075, burial groups halved. Trainer copy `ff30_bio_si/trainer/ConDiv.py` (md5 `a1943099...`; differs from `$P/training/ConDiv.py` only in those lines, their comment and docstring); initial force field byte-identical to ff30_bio_dt009's (`checks/si_init_20261007`, build log and trainer diff there); target 76, gate `ff30biosi-gate` up to 13 epochs (converged: validated as `ff_3.0_bio_si`); first-link log `$P/training/condiv-train_49204255.out` | First link 49204255 (from 16:34 10-07, 14 nodes) **COMPLETED 10-08 20:10** after its 54 steps (0-53, `epoch_02_minibatch_15`; 26-41 min each). **Successor 49205039 R since 08:49 10-09** (16 nodes; log `ff30_bio_si/condiv-train_49205039.out`: `step 54 of 76, resuming from run_output/epoch_02_minibatch_15/checkpoint.pkl`, `this link runs 22 steps`, owed check done), carrying the law; steps 54-67 done 09:32-20:24 (margin +0.060, +0.058, +0.055, +0.053, +0.050, +0.047, +0.043, +0.040, +0.038, +0.037, +0.038, +0.040, +0.040, +0.042; 24 of 24, KE/1.5kT at most 1.016); **`bs_02`** (step 56) extracted 10:58 (rama.dat md5-identical to bs_01's; hbond.h5 and sheet moved), panel 49216513 (10:57, carries the law; estimate 10-10 22:05), patched, **TM4 `bs_02`** (Mac Studio 12:43-14:38), primary test: DSSP alpha 0.684 against bio_start's 0.769 (p 0.58) and its twin b9_02's 0.694 (p 0.93), not resolved; strong at 135-138, weak at 149-151 (b9_02 the reverse); 5 unwound, 2 flipped, no energy jump (findings 1.25); its successor 49216105 queued with the law; **`bs_00`** (step 18, 02:09) extracted, panel 49207054 COMPLETED 21:39: **`bs_00` 0.457, helix -0.028** on 34 domains, dominated by gdepth_start in helix and not by its twin b9_00 (findings 1.21); **`bs_01`** (step 37, 11:18) extracted 11:44 (rama.dat md5-identical to bs_00's, hbond.h5 moved), panel 49208587 COMPLETED 22:41: **`bs_01` 0.430, helix -0.033** (36 domains; b9_01 0.465, -0.031), dominated by its twin b9_01 in helical glycine (findings 1.21); **TM4 `bs_01`** (Mac Studio 15:56-17:53): unwound 5, flipped 2 of 12, against bio_start's 5 and 4 (p 1.00, 0.64) and its twin b9_01's 6 and 5 (p 1.00, 0.37), not resolved, leaning better than b9_01 (last-block TM4 0.879 against 0.833; findings 1.24); **TM4 `bs_00`** (Mac Studio 08:12-10:08): unwound 6, flipped 5 of 12, against bio_start's 5 and 4 (p 1.00) and its twin b9_00's 4 and 4 (p 0.68, 1.00), not resolved, leaning worse than b9_00 (findings 1.24); 24 of 24 each, KE/1.5kT at most 1.015. **Rates confirmed** at step 0: checkpoint `initial_alpha` hb 0.005, dhb 0.0025, sheet 0.0075, rot 0.0125, and the first step moved each parameter by exactly that (ff30_bio_dt009's step 0: 0.01, 0.005, 0.015, 0.0125); the link log prints no rates. Margin +0.192, +0.190, +0.186, +0.181, +0.176, +0.169, +0.162, +0.155, +0.149, +0.142, +0.135, +0.131, +0.127, +0.125, +0.122, +0.121, +0.119, +0.117, +0.115, +0.115, +0.113, +0.111, +0.109, +0.107, +0.104, +0.100 (+0.1001), **+0.097**, +0.095, +0.092, +0.090, +0.088, +0.087, +0.088, +0.088, +0.088, +0.086, +0.086, +0.084, +0.083, +0.081, +0.079, +0.076, +0.072, +0.069, +0.066, +0.063, +0.062, +0.060, +0.060, +0.059, +0.059, +0.060, +0.060, +0.060 at steps 0-53, +0.060, +0.058, +0.055, +0.053 at steps 54-57. **HOLD since step 26** (06:06; user told 06:50); it trains on (ff30_bio_dt009 +0.144 at step 3, +0.126 at step 7); dhb -0.416, sheet mean 0.171 at step 31 | watch as the dt 0.009 runs; tags `bs_EE`; TM4 against bio_start and ff30_bio_dt009 at the same step |
| **49216034** (midway2) | **ff30_gdepth_si**: ff30_gdepth_dt009 at the same SI rates, trainer copy `ff30_gdepth_si/trainer/ConDiv.py` (md5 `90ae572c...`); initial force field and `rama_round_00.dat` byte-identical to ff30_gdepth_dt009's; depth start the same (alpha_R -0.0918, alpha_L +0.4621); target 76, gate `ff30gdepsi-gate` (converged: `ff_3.0_gdepth_si`); first-link log `$P/training/condiv-train_49204256.out` | Link 49204256 (R from 17:20 10-07, 24 nodes) **NODE_FAIL 10-08 11:07** in step 21 (`epoch_01_minibatch_02` half-written): its batch host midway2-0027 failed, the third such failure today; the link started before the law excluded 0027. Steps 0-20 done (18:00-10:05; 29-63 min). Link 2, 49205173 (from 08:40 10-09, 16 nodes; `step 21 of 76, resuming from run_output/epoch_01_minibatch_01/checkpoint.pkl`, owed check done), finished step 21 (09:09, margin +0.128, 24 of 24, KE/1.5kT at most 1.014) and **NODE_FAILed 10-09 09:17:42 in step 22** (`epoch_01_minibatch_03` half-written), with ff30_bio_dt009's link; its node midway2-0116 shows NOT_RESPONDING afterwards (§8). **Successor 49216034 PD (Priority, estimate 10-10 22:05 at 11:47)**, carrying the law, resumes from step 21 by chain design (owed: its first line `step 22 of 76, resuming from run_output/epoch_01_minibatch_02/checkpoint.pkl`); **`ds_00`** (step 18) extracted, **TM4 `ds_00`** (Mac Studio 21:45-23:41): unwound 9, flipped 3 of 12, against gdepth_start's 7 and 2 (p 0.67, 1.00) and its twin d9_00's 5 and 2 (p 0.21, 1.00), not resolved; last-block TM4 0.775, the lowest mean of any set (findings 1.24); panel 49207993 COMPLETED 22:39: **`ds_00` 0.412, helix -0.030** (36 domains; d9_00 0.477, -0.024), dominated by its twin d9_00 in helix (findings 1.21), TM4 queued. Depth round 1 (`rama_rounds.txt`): gap aR +0.001 aL +0.008, offsets aR -0.0882 aL +0.4754, dL - dR +0.564; 24 of 24, KE/1.5kT at most 1.015. **Rates confirmed** at step 0, as ff30_bio_si's (`initial_alpha` hb 0.005, dhb 0.0025, sheet 0.0075, rot 0.0125, and the first step moved each parameter by exactly that; ff30_gdepth_dt009's step 0: 0.01, 0.005, 0.015, 0.0125). Margin +0.192, +0.191, +0.185, +0.179, +0.175, +0.170, +0.168, +0.165, +0.162, +0.160, +0.157, +0.155, +0.152, +0.149, +0.147, +0.144, +0.142, +0.140, +0.136, +0.134, +0.132, +0.128 at steps 0-21 (ff30_gdepth_dt009 +0.183 at step 3); depth round 0, dL - dR +0.554 | watch as the dt 0.009 runs; tags `ds_EE`; TM4 against gdepth_start and ff30_gdepth_dt009 |
| **49214161** (midway2) | **ff21_ctrl_fz**: the frozen arm's matched control, resubmitted unchanged after its 13:12 cancellation (user, 10-07 13:35): ff21_ctrl_dt009 with hb, dhb, hbg and sheet at learning rate 0 (trainer copy md5 `a97e95a6...`, the same as ff30_bio_fz's), ff2.1's `rama.dat`; initial force field is ff2.1 exactly (`checks/fz_init_20261007`); **target 38, no gate**; first-link log `$P/training/condiv-train_49204264.out` | **First link 49204264 NODE_FAIL 10-08 00:12** in step 2: its batch node midway2-0027 failed (`CANCELLED ... DUE TO NODE FAILURE`; 0027 is allocated again, now as ff30_bio_dt009's batch host). Steps 0-1 done (22:10, 23:17; 54 and 67 min), 24 of 24, KE/1.5kT at most 1.012, margin +0.1919 (frozen). **Successor 49205852 R since 22:28 10-08** (12 nodes; log `ff21_ctrl_fz/condiv-train_49205852.out`: `step 2 of 38, resuming from run_output/epoch_00_minibatch_01/checkpoint.pkl`, `this link runs 36 steps`, owed check done); its successor 49214161 queued 22:28, carrying the law. Steps 2-10 done (23:01-03:42 10-09; 28-45 min each; step 10 relaunched 12 workers, 5 twice, and returned 24 of 24). **Link 49205852 FAILED 10-09 04:11 in step 11** (exit 1, `RuntimeError: 3 of 24 workers failed: 2r2y 3jyz 4qbo`): 8 workers' srun launches failed at 03:42, and those three failed all three on midway2-0088 (`Invalid job credential`). Steps 0-10 stand. **Successor 49214161 PD (Priority, estimate 10-10 22:05 at 11:47)**, which resumes from step 10 by chain design; midway2-0088, the node's second such event (§8), is excluded by the law since 04:12 (user), and 49214161 carries it (owed: its first line `step 11 of 38, resuming from run_output/epoch_00_minibatch_10/checkpoint.pkl`), 24 of 24, KE/1.5kT at most 1.015, margin +0.1919. **Freeze check passed** (22:12; `checks/fz_init_20261007/cz_step00`, log `extract_cz_step00.log`): hbond.h5 and sheet md5-identical to `checks/fz_init_20261007/ff21_ctrl_fz` and to `$P/parameters/ff_2.1`, rama.dat unchanged; sidechain.h5, environment.h5, bb_env.dat moved; check_step: hb / dhb / sheet mean unchanged, margin +0.1919 | the watch every step; tags `cz_EE`; `cz_00` (step 18) ~4 h after 49214161 starts: TM4 as `cz_00` and `cz_00 13 24`, against ff21_released and c9_00 at 24 seeds |
| **49206070** (midway2) | **ff30_gdepth_fz**: ff30_gdepth_dt009 with hb, dhb, hbg and sheet at learning rate 0, resubmitted unchanged after its 13:12 cancellation (user, 10-07 21:45: run every frozen arm, though it departs from the group's workflow). Trainer copy `ff30_gdepth_fz/trainer/ConDiv.py` (md5 `84704d98...`), which differs from `ff30_gdepth_dt009/trainer/ConDiv.py` only in that rate line, its comment and one docstring item; `init_param/` (all five files, so hbond.h5 and sheet are ff2.1's), `rama_round_00.dat`, `rama_basin.py`, `pdb_list`, `env.sh` and the `upside_input` listing identical to ff30_gdepth_dt009's (checked 21:45); target 76, gate `ff30gdepfz-gate` (converged: `ff_3.0_gdepth_fz`); first-link log `$P/training/condiv-train_49206070.out` | **First link 49206070 R since 21:34 10-08** (17 nodes; log `training/condiv-train_49206070.out`: `step 0 of 76, resuming from run_output/initial_checkpoint.pkl`, `this link runs 54 steps`); successor 49214121 queued 21:35, carrying the law with 0027. Steps 0-42 done (22:14-21:14 10-09; 24-54 min; step 2 relaunched 2krc and 4bou once, step 10 2kph and 2ejx, step 25 2j6b, step 27 4ous and 1tig; 24 of 24 every step), 24 of 24, KE/1.5kT at most 1.016; hb, dhb and sheet mean unchanged (margin +0.1919). Depth round 1 at the epoch-0 end (`rama_rounds.txt`): gap aR -0.002 aL +0.005, offsets aR -0.0965 aL +0.4704, dL - dR +0.567 (ff30_gdepth_dt009's round 1: +0.567); round 2 at the epoch-1 end: gap aR +0.000 aL +0.004, offsets aR -0.0960 aL +0.4775, dL - dR +0.574 (ff30_gdepth_dt009's: +0.577). **Freeze check passed** (22:50; `checks/fz_init_20261007/dz_step00`, log `extract_dz_step00.log`): hbond.h5 and sheet md5-identical to `checks/fz_init_20261007/ff30_gdepth_fz` and to `$P/parameters/ff_2.1`, rama.dat unchanged (depth round 0, offsets aR -0.0918 aL +0.4621); sidechain.h5, environment.h5, bb_env.dat moved | the watch every step; tags `dz_EE`; **`dz_01`** (step 37, 18:36): extracted 18:48 (hbond.h5 and sheet md5-identical to gdepth_start's; rama.dat is its own `rama_round_02.dat`), panel 49218034 (18:46, carries the law), patched (against dz_00's patched input only rama, side-chain pair and coverage tables differ), **TM4 `dz_01`** (Mac Studio 18:49-20:45), primary test: DSSP alpha 0.780 against gdepth_start's 0.657 (+0.123 [-0.100, +0.336], p 0.58) and its twin d9_01's 0.615 (+0.165 [-0.056, +0.385], p 0.069), not resolved; strong at 144-149, weakest at 150-151; 5 unwound, 3 flipped (all GLY143), two energy jumps (findings 1.25); **`dz_00`** (step 18, 08:23): extracted (hbond.h5 and sheet md5-identical to gdepth_start's; rama.dat is its own `rama_round_01.dat`), panel 49216103 COMPLETED 09:36: **`dz_00` 0.432, helix -0.033** (41 domains), dominated by gdepth_start in helix and by its twin d9_00 (0.477, -0.024 on 38 domains; findings 1.21); patched (against gdepth_start's patched input only rama, side-chain pair and coverage tables differ), **TM4 `dz_00`** (Mac Studio 08:50-10:46), primary test: DSSP alpha 0.791 against gdepth_start's 0.657 (+0.134 [-0.074, +0.346], p 0.47) and its twin d9_00's 0.801 (p 0.80), not resolved; strong at 136-141 and weakest at 143-151, the reverse of d9_00; secondary counts 4 unwound, 5 flipped; s11's jump of 23501 is the largest on fixed inputs (findings 1.24-1.25) |
| **49206266** (midway2) | **polygly production** (plan.md Phase 12, step 4; user 10-07): Ac-(Gly)20-NHMe, 4 replicas from the collapse's last frames in 7.5 nm dodecahedra (~28,950 atoms each; `polygly/scripts/build_prod.sh`, log `polygly/logs/build_prod.log`), 1 us cap, 7 threads each on one broadwl node. `polygly/prod.sbatch` minimises, runs 500 ps NPT (seeds 20261021-24), then production. It self-chains: successor queued at link start, cancelled when all four reach 1 us; three links from the same steps, or `polygly/prod/STOP`, end it. Logs `polygly/logs/prod_<jobid>.out`; replicas `polygly/prod/repN/` | **Link 1, 49204428, COMPLETED 10-09 09:12** at its `-maxh 34.5`: 24.08-24.19 M steps (48.2-48.4 ns) per replica, `prod.log` Performance 33.99 / 33.91 / 33.85 / 33.84 ns/day; `gmx mindist -pi` over the link (4816-4838 frames; `polygly/prod_analysis/link1_49204428/`): minimum image distance 2.35 / 2.97 / 3.29 / 3.57 nm (rep1-4), above 1.0 nm. **Successor 49206266 PD (Priority)**, carrying the law, continues from the checkpoints; ~20 more links to 1 us. Link 1: **R since ~22:41 10-07** on midway2-0612 (`polygly/logs/prod_49204428.out`: GROMACS 2024.4, AVX2_256, starting steps 0 0 0 0 of 500000000); successor 49206266 queued. Minimisation (`mdrun_min.log`) stopped at its 5000 steps without reaching Fmax < 100: Fmax 1326 / 1411 / 223 / 543 kJ/mol/nm on protein atoms 110 / 117 / 131 / 47, potential -4.758e5 to -4.764e5, as the collapse's did (1380, atom 124) before it ran cleanly. 500 ps NPT equilibration finished 23:03 in all four (no LINCS or SETTLE warnings; `equil.log` Performance 33.7-33.8 ns/day per replica); production at 1.1 ns by 23:47, so 1 us takes ~30 days (~20 links). `gmx mindist -pi` on each replica's first 1.1 ns (protein-only `prot.tpr`, 112-113 frames): minimum image distance 3.54 / 4.53 / 4.55 / 4.94 nm, above 1.0 nm | the watch: the next link's first lines; at each link end the Performance line and mindist -pi; the blank at each ~200 ns (~4 links) |

**Local glpG TM4 test** (the hybrid glpG 79HIS seed in POPE/POPG dry-MARTINI, one checkpoint
patched in, T 0.80, 4000 tu, 12 seeds; `scratchpad/ff3_local_test`, record in
`checks/r4_epochs/tm4_local`).
* **Valid sets (coverage-fixed inputs, from 10-07 11:22):** `ff21_bioT1_6` (bio_start) and `b9_01`
  on the MacBook Pro; `ff21_released`, `gdepth_start`, `d9_00`, `b9_00` (24 seeds), `c9_00` (24 seeds), `fp_e00`, `c9_01`, `d9_01`, `d9_02`, `d9_03`, `bs_00`, `bs_01`, `bs_02`, `bz_00` (24 seeds), `bz_01`, `bz_02`, `b9_02`, `ds_00`, `dz_00` and `dz_01` on the Mac Studio, whose engine
  and inputs reproduce the MacBook Pro's frame for frame. The rest is the queue in "Handoff". Each
  is in `tm4_local/runs_cov/`, with its table `tm4_compare_cov_<set>.txt`.
* **Mac Studio, from 10-07 16:20, the owner** ("Handoff"). Driver `scripts/tm4_queue.sh` under
  `nohup caffeinate -is`; it pops the first line of `tm4_queue.txt` (`<tag>` or `<tag> 13 24`),
  runs it and logs to `run_<tag>.log` or `run_<tag>_s13-24.log`; its own log is `tm4_queue.log`.
  21:27 10-09: idle; the driver finished `dz_01` at 20:45 and exited with `tm4_queue.txt` empty. Its pre-fix sets are in `runs_precov_20261007/`.
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

* **ff21_ctrl_dt009 (49200579) COMPLETED 10-07 21:25 at its target 38** (`target 38 reached`; the chain cancelled successor 49200622; no gate). The matched control: ff30_bio_dt009 with ff2.1's `rama.dat`, so it starts from ff2.1 exactly (`checks/ctrl_dt009_init_20261006`); `$P/training/ff21_ctrl_dt009`, log `condiv-train_49200579.out`. Steps 0-37, 24 of 24 each, KE/1.5kT at most 1.017. Margin +0.132 at step 14, +0.135 at step 18, +0.125-0.130 over steps 19-31, then +0.134 / +0.139 / +0.147 / +0.152 / +0.155 / **+0.157** at steps 32-37; dhb -0.468, sheet mean 0.174 at step 37. Panel **`c9_00` folded 0.390, helix -0.032** (findings 1.22); TM4 `c9_00` unwound 7, flipped 1 against ff21_released's 9 and 5, not resolved (findings 1.24). `c9_01` (step 37) extracted to `checks/r4_epochs/c9_01` (21:47); **panel `c9_01` (49206076, 22:50-23:38): folded 0.405, helix -0.029** on 34 domains, c9_00 0.390, -0.031 on the same set (findings 1.22); **TM4 `c9_01`** (Mac Studio 04:19-06:15): unwound 10, flipped 6 of 12, against ff21_released's 9 and 5 (p 1.00) and c9_00's 7 and 1 (flipped two-sided p 0.069), not resolved (findings 1.24).
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

**ff30** (full-map glycine row) was cancelled 2026-09-28 00:25 at step 223; the run dir (54 G) and
last checkpoint `ff30/run_output/epoch_11_minibatch_13` are kept for comparison; extract with
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

`/project` is the one to watch: ~965 G free (2026-09-30); the glpG REMD trees hold ~1.26 T
(`popepopg_REMD_mdw2`'s four variant directories ~150 GB each, `popepopg_REMD` another 638 G), and
`NP-1AO6` at ~0.5 T is the reclaim.

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
* **Deletion decisions waiting on the user:**
  * `NP-1AO6/` data, including `prod_ff3/` (~500 GB).
  * midway3 DDM data `glpG_DDM_micelle_REMD/`, `glpG_DDM_REMD/`.
  * rockfish post-fix glpG data `/scratch4/rherna21/ywang268/upside2-md-rf/popepopg_REMD/<variant>/`
    (~72 GB, needed only for an independent-seed cross-check), before scratch is reaped.
* **Safe to delete:** `/beagle3/trsosnic/yinhan/ff3_benchmark/sigmoid_discarded/` (4 G) and
  `armM_discarded/`.
* **Open, unverified, recorded so they are not lost:**
  * The HDX topology uses ff_2.1 tables on an ff_3.0 trajectory.
  * TM4 shields only ~5 residues in HDX, not checked against the crystal structure.
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
