# Remote jobs on midway2/midway3: status and handbook

**Current state (2026-10-02 15:00).** The ff3.0 retrain is **ff30_glyhb** (plan.md Phase 8,
revised 10-02), chain link 49141995 on midway2 broadwl since 13:25:41, from ff2.1 with glycine's own
H-bond offsets and a 10x smaller side-chain step; steps 0 and 1 done and healthy. ff30_gly (the AWH
library alone) was stopped at 11:51 in epoch 3. Run only on midway2 broadwl (user). midway3's login
node may be used to move or read files on the shared `/project` and `/beagle3` (user, 2026-10-02).
The BP validation is analysed (§0c) and not yet reported to Tobin; its summary figure is being
made (job 49143602 dumps the per-seed REMD arrays, §0c "Figure").

Written so a fresh session can pick up cold: how to connect, check health correctly and react to a
failure. Job state in §1 is live; finished jobs are kept only where their files are still used (§1)
or as one-line lessons (§8). rockfish is gone as a host, so glpG is midway2-only.

---

## Resume here, from any computer (2026-10-02)

**What travels with git and what does not.** `plan.md`, this file, `findings.md`, `progress.md`,
`training/`, `src/`, `py/` and `up.md` are tracked; commit them here and pull them there. Not in git
(`.gitignore` has `*scratchpad*`): `scratchpad/mdw2_master.exp`, `mdw3_master.exp` and the local
copies of the panel tools. The authoritative panel tools are on the cluster, in
`/project/trsosnic/yinhan/ff3_selection` (`panel.py`, `panel.sbatch`, `submit_new.sh`,
`slurm.args`), and the deploy and test scripts of 10-02 are in
`/project/trsosnic/yinhan/checks/hbg_deploy_20261002`. Claude's memory is per computer; the rules
that matter here are in this file and `findings.md` §10.

**The cluster's `training/` is not the repo's layout, and must stay so until ff30_glyhb is done.**
The repo was cut to five files on 10-02 (`extract` and `gate` became ConDiv commands,
`gate_or_continue.sh` folded into `train_chain.sbatch`, `check_step.py` moved to
`scratchpad/redistribution_cleanup_20261002/training_pre_merge/`). `$P/training` on midway2 keeps the
old layout, deployed 10-02 11:46 with the glycine term: `ConDiv.py` (the copy ff30_glyhb's
`run_output/ConDiv.py` was made from), `check_step.py`, `extract_ff.py`, `convergence_gate.py`,
`gate_or_continue.sh`, `validate_ff.sh`, `patch_glpg.py`. ff30_glyhb's `after_training.sbatch` calls
`$P/training/gate_or_continue.sh`, and the watch calls `$P/training/check_step.py`, so **do not sync
the repo's `training/` over `$P/training`** before ff30_glyhb is released; a synced tree would end
the run at its gate with a missing script.

**Connect.** The simplest way from a new computer is to open the master yourself, password and Duo:
```bash
ssh -M -S ~/.ssh/cm-mdw2.sock -o ControlPersist=8h -o ServerAliveInterval=30 yinhanw@midway2.rcc.uchicago.edu
```
and leave it open (or exit; ControlPersist keeps the socket). Every later command reuses it:
`ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<cmd>'`. To let Claude
open it unattended, copy `scratchpad/mdw2_master.exp` over (§0 shows what it does) and put the
password line `set password "..."` in `~/.bin/ssh_mdw3` there. **Never run anything heavy on a login
node** (findings 10.6): a 60 GB analysis on 10-02 got every session on midway2-login1 killed five
times in 40 min. If the socket keeps dropping, run `ps -u yinhanw --sort=-rss` there first.

**One watch at a time.** This computer's watch is session cron `6733412b` (:17 and :47); it dies
with that Claude session. If the other computer takes over, stop this one (or end the session), and
start the same watch there with the prompt below, every 30 min. `submit_new.sh` is idempotent, so
two watches would not double-submit, but they would both report.

**The watch, by hand or as the cron prompt:**
1. `ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu`; if down,
   reconnect (above). Never loop on connection attempts, and never run the expect script from an
   unattended watch.
2. On midway2: `bash /project/trsosnic/yinhan/ff3_selection/submit_new.sh` (submits the panel `hEE`
   for each finished epoch-end checkpoint of ff30_glyhb, lists jobs and finished panel runs).
3. Link log: newest `/project/trsosnic/yinhan/upside2-md-mdw2/training/ff30_glyhb/condiv-train_*.out`,
   for `WORKER_FAIL`, `Traceback`, `STOPPED`, `never started`.
4. Each new finished step: `cd /project/trsosnic/yinhan/upside2-md-mdw2/training && source
   ff30_glyhb/env.sh && python3 check_step.py ff30_glyhb` (all finite, KE/1.5kT ~1.0-1.05, restrained
   RMSD ~1 A, unfolded target 24 of 24; the glycine offsets and margins; the glycine readout).
5. When `ff3_selection/runs/<tag>/` holds 176 npz: `cd /project/trsosnic/yinhan/ff3_selection &&
   source /project/trsosnic/yinhan/upside2-md-mdw2/training/ff30_glyhb/env.sh && python3 panel.py
   select aa domains runs ff21_released ff21_awh e00 e01 e02 h00 [h01 ...]`. It prints the table,
   the chosen checkpoint and whether the release is held. Nothing is released by the watch: the
   decision goes to the user.
6. Update §1 below.

**What happens next, and when** (steps take 37-43 min, so an epoch of 19 steps takes ~12.5 h):

| when (CDT, estimates) | event | what to do |
|---|---|---|
| every ~40 min | a training step | watch items 3-4 |
| 10-02 ~19:00-20:00 | lambda ff2.1 chunk 2 ends (4 arms, self-resubmitting; ~450k of 2.53M time units per chunk) | `cd /project/trsosnic/yinhan/checks/lambda_ff21 && python3 lambda_check.py <arm> ...` and `lambda_packing.py <arm> ...`, arms `lambda_{native,denovo}_ff_2.1`, `lambda_G46A_G48A_{native,denovo}_ff_2.1`; on a compute node or with `ulimit -v` |
| 10-03 ~01:00-02:30 | epoch 0 ends (`epoch_00_minibatch_18`) | the watch submits panel `h00` (1 node, 1-3 h) |
| 10-03 morning | `h00` table | **first verdict on the glycine term**, see below |
| 10-03 ~14:00-15:00 | epoch 1 ends, panel `h01` | as for h00 |
| 10-04 ~01:30 | link 49141995 hits its 36 h wall at ~step 54; successor 49142816 continues | check it started and resumed from the newest checkpoint |
| 10-04 ~04:00 and ~16:00 or later | epochs 2 and 3 end (h02, h03); at step 76 `after_training.sbatch` (`ff30h-gate`) runs the convergence gate | converged: it stops and lists checkpoints, choose by panel; not converged: one more epoch, up to 13 |
| 10-04 to 10-05 | lambda ff2.1 arms finish | final lambda readout (findings 11d) |

**Reading the h tables.** Compare each hEE against ff21_released and ff21_awh in the same table (37
domains common to e00-e02 at 13:10; the set can change as tags are added). The aim: helical glycines
(gly_helix) at least as close to all-atom as ff2.1 (-0.114 at 13:10) while left-handed glycines
(gly_left) stay near ff2.1's (+0.032), and helix and folded no worse than ff2.1 (-0.016, 0.615).
ff30_gly reached gly_helix -0.112, gly_left -0.042, helix -0.025, folded 0.501 at e02, and was
held for helix. **Trend to watch on every step:** check_step's glycine margin
`E_other - E_alpha` for glycines (start +0.192, +0.139 after step 1). A margin that keeps falling
means glycines are pushed toward alpha_L, the direction that cost helical glycines before.

**After ff3.0 is chosen** (plan.md Phase 8): release from midway2 with
`bash ../validate_ff.sh . ff_3.0 run_output/epoch_EE_minibatch_18/checkpoint.pkl` in the run dir;
sync `/beagle3/trsosnic/yinhan/upside2-md` from `$P` first (it lacks the glycine term); remove the
`STOP` files in `popepopg_REMD_mdw2/glpG-RKRK-*/` and reset `block_count` before glpG; the glpG seed
takes the new term through `patch_glpg.py`; re-simulate lambda for its packing (findings 11d).

---

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
/ RMS force median 8.3e-5, max 2.9e-3 (fixed 4.6e-5, 5.3e-4); ConDiv contrastive rotamer-gradient
error 1.4-2.5e-4 relative (fixed 6-8e-5). Energy bias one-signed: E_bug < E_ref in 95% of frames,
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

**Figure for Tobin (in progress 2026-10-02).** `plot_bp.py` makes one 11-panel figure: a-c static
frames (ECDF of |E - E_ref|, ECDF of relative force error, histogram of iterations cut short), d-g
Q(T) per seed for the four REMD proteins, h-k Q_bug - Q_fix with the seed sigma. Red = bug, blue =
fix in every panel. SciencePlots by path (global CLAUDE.md plotting standard).
* REMD arrays: job **49143602** (`dump.sbatch`, broadwl, 4 cores, 40 min) reruns
  `analyse_bp.py` per protein. `analyse_bp.py` now also writes `obs_<protein>.npz` (per-arm, per-T
  Q, RMSD, Rg, E; original kept as `analyse_bp.py.bak_pre_npz`). Log `logs/dump_<jobid>.out`, text
  `analysis_rerun.txt`, which must reproduce `analysis_20260929.txt` (a reproducibility check).
* Static-frame harness and results: `bp_validation/static/` (`gen.py` frames, `eval.py` one solve per
  library and tolerance, `plot_bp.py`, `sys/<code>/{frames,init}.npy`, `base.up`, `res_{bug,fix,ref}.npz`).
  `eval.py` needs Mac builds of the bug and fix libraries with a `read last_iter` readout added
  to `rotamer.cpp` (iteration count); these and the working copy are in the gitignored
  `scratchpad/bp_validation/` of the Mac that made them, not in git.
* To finish from any computer: when 49143602 is done, copy `bp_validation/obs_*.npz` into
  `static/remd/`, copy `static/` locally, and run `python plot_bp.py` there (it writes
  `bp_validation.png` and `.pdf`; it loads SciencePlots from the path in the global CLAUDE.md).

---

## 1. Current jobs

Snapshot **2026-10-02 15:24 CDT, verified live against `squeue` on midway2.** Finished and cancelled
rows are deleted; their lessons are in §8.

| JobID / where | what | state | next action |
|---|---|---|---|
| **49141995** | **ff30_glyhb: ff3.0 from ff2.1 with glycine's own H-bond basin offsets (`hbg`, hbond.h5 entries 12-14, from zero), side-chain lr 10x smaller, AWH glycine library** (plan.md Phase 8), chain link on broadwl, `$P/training/ff30_glyhb`, target 76 steps, gate up to 13 epochs (`after_training.sbatch`, job name `ff30h-gate`); log `condiv-train_<jobid>.out` | R since 13:25:41 on 19 nodes; step 0 done 14:02 in 2196 s (worker median 1229 s), healthy (`check_step.py`: 24 of 24 finite, KE/1.5kT 1.006-1.041, restrained RMSD median 1.03 A, unfolded target 24 of 24); step 1 done 14:46 in 2586 s, healthy (24 of 24 finite, KE/1.5kT 1.008-1.046); hbg after step 1 [+0.019 +0.019 -0.019], glycine margin E_other - E_alpha +0.139 (start +0.192), shared margin +0.177; link log clean; at ~35 min a step, the 36 h wall holds ~60 of the 76 steps and the successor finishes them | `check_step.py ff30_glyhb` on every new step; panels h00, h01, ... by the watch |
| 49142816 | its insurance successor (`afterany:49141995`) | PD (Dependency) | none; resumes from the newest checkpoint if the link dies |
| 49139947, 49139944 | Peng benchmark lambda under **ff_2.1**, native and de novo arms, chunk 2 (Table S2 ladder 0.780-0.980, self-resubmitting 12 h jobs); `/beagle3/trsosnic/yinhan/ff3_benchmark/runs/lambda_{native,denovo}_ff_2.1/`, logs `logs/lambda_<kind>_ff_2.1_<jobid>.out` | R | after each chunk: `checks/lambda_ff21/lambda_check.py` and `lambda_packing.py` |
| 49139945, 49140670 | the same arms for **lambda G46A/G48A** (helix 3's glycines to Ala), chunk 2; `runs/lambda_G46A_G48A_{native,denovo}_ff_2.1/` | R | compare with wild type at matched time |
| 49143602 | **bp_dump**: per-seed REMD arrays for the BP-validation figure (§0c "Figure"); `/beagle3/trsosnic/yinhan/bp_validation/dump.sbatch`, reruns `analyse_bp.py` on the 16 finished arms, writes `obs_<protein>.npz` and `analysis_rerun.txt`; log `logs/dump_<jobid>.out` | PD (Priority) since 15:03 | check `analysis_rerun.txt` matches `analysis_20260929.txt`, then copy `obs_*.npz` to `static/remd/` and run `plot_bp.py` |

**Panel e02 finished 13:05** (49140630). `select` with ff21_released ff21_awh e00 e01 e02 (37
domains common to all; 7 dropped for too few folded frames in some candidate): e02 folded 0.501,
helix -0.025, beta -0.017, gly_helix -0.112, gly_left -0.042 against ff2.1's 0.615, -0.016, -0.017,
-0.114, +0.032; "RELEASE HELD: worse than ff2.1 in helix". ff30_gly is stopped; these are the
comparison for the h tags.

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

**Selection watch** is session cron `6733412b` (:17 and :47, replaces `39a1fa34`), watching
ff30_glyhb with tags h00, h01, ...; `submit_new.sh` backup `.bak_pre_glyhb_20261002`. It submits the
panel for each finished epoch-end checkpoint (marker `runs/<tag>.submitted`), checks the link log,
and runs `panel.py select` when data are complete. Session-only, expires 10-09; never runs a Duo
script. The release decision goes to the user, by the rule in `panel.py`.

`$P` = `/project/trsosnic/yinhan/upside2-md-mdw2`. ff30_gly's cluster run_output holds the local
Mac run's synced steps (marker `run_output/FROM_LOCAL`); the cluster's own initialised run_output is
kept as `run_output.superseded_20261001-111755`. The Mac Studio's local run
(`training/ff30_gly_local`) was stopped by the user (confirmed 2026-10-02 00:25); its steps after
the cluster start (10-01 12:16:47) are discarded.

**Selection panel data** (`/project/trsosnic/yinhan/ff3_selection`): controls complete,
`runs/ff21_released/` (49136636) and `runs/ff21_awh/` (49136637 NODE_FAIL on midway2-0085, finished
as 49136943, 17:47); all-atom `aa/` (44 domains); `runs/e00/` (49137240, COMPLETED 19:37);
diagnostics `runs/e00_hb21/`, `runs/fp_e00/` (the ff21-fixedpoint epoch-0 control) and the
single-group resets `runs/e00_{bb21,env21,rot21,sheet21}/` (all COMPLETED 10-02 by 10:55), findings
1.17. Candidates `e01`, ... are submitted by the watch.

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
and TM1 1.000 (findings 3.10c).

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
* **A failed worker fails the step** (`ConDiv.py` raises): step 81 of ff30 once ran on 5 of 24
  proteins after `srun: Job credential expired`, and a partial step is a different objective. That
  credential error at launch is a race, not a node fault (24 steps issued together lost 4 on healthy
  nodes); the trainer relaunches a worker that never started.
* **Node exclusions are per record of failures.** midway2-[0010-0011] (three NODE_FAILs in one
  campaign) and midway2-0037 (two links in three hours on 2026-09-26) are excluded in `slurm.args`;
  two later link failures named no node and excluded nothing.
* **`sbatch --test-only`'s start estimate is no guide to the real wait**: it predicted 13:36 on both
  clusters, then the midway3 submission sat PENDING (Resources) while midway2 started at once.
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
