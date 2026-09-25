# Remote jobs on midway2/midway3 — status and handbook

Snapshot: **2026-09-19 09:00 CDT.** Eight jobs (six long-running plus `scmiss` and `ff21-restart`), all healthy, none needing intervention.
Glycine AWH campaign relaunched after a settings defect made the first attempt unconverged
(49032988/49033468 cancelled; 49033509/49033947/49033985 are the corrected run; see Campaign 6).
The go/no-go is answered: per-neighbour structure is noise, only the neighbour-average is usable,
so the three AWH jobs run out their walls for that average and nothing further is spent on them.
ff3.0C training stays cancelled. The PI email and both figures in `~/Downloads`
are marked `UNCONVERGED_DO_NOT_USE` / `HOLD`. Earlier state, still current unless noted:
rockfish is GONE as a host, so glpG is **midway2-only** and the two-cluster replicate design is
retired. glpG post-fix healthy; NP running clean on the rollback driver. **midway2 IP block lifted**
as of 2026-09-12; key-based ssh is refused, every connection costs a Duo push. See findings.md
3.10-3.10d for the temperature-fix diagnosis.**
Written so a fresh session can pick up cold. Everything needed to connect, check health correctly,
and react to a failure is here. Job state below is live; superseded jobs are not listed, only
summarised in §8 where they carry a lesson.

---

## 0a. `broadwl-lc` is EMPTY BUT UNUSABLE: its nodes cannot see `/project` (2026-09-18)

**Do not send work there, however idle it looks.** Measured from a job on `midway2-0213`:
`PROJECT_MISSING`, `PROJECT_NOT_WRITABLE`, `BEAGLE3_MISSING`. The nodes advertise
`AvailableFeatures=lc,e5-2680v4,64GB,`**`noib`** - no InfiniBand, and `/project` and `/beagle3` are
served over it. Every input we use (training set, `.venv`, `libupside.so`, checkpoints) is on
`/project`, so the partition is dead for this project. This is the same trap as `/cds3`.

**The failure signature is silent and easy to misread.** A job whose `--output` path is on
`/project` dies **instantly** with `ExitCode 0:53`, `Elapsed 00:00:00`, and **no log file at all**,
because it cannot create the log. Nothing says "filesystem". Do not read that as a bad script: a
`--wrap="hostname"` probe of the same shape fails identically.

**Two things that wasted time here, worth not repeating.** `sinfo -F`'s "0 idle" counts a
partially-filled `mix` node as allocated, so it does **not** mean no cores are free; use
`sinfo -p <part> -o "%.8t %.6D %.10C"` for A/I/O/T cores. And **`/tmp` on a compute node is
node-local** - a probe writing there "succeeds" and leaves nothing the login node can read, which
looks like another silent failure. Probe with `--output` under `$HOME`, which is shared.

**A job_submit plugin silently rewrites `#SBATCH --partition=broadwl-lc` to `broadwl`**, with no
warning, while honouring every other directive. Command-line `-p` is honoured. Verify placement
with `scontrol show job <id> | grep Partition` rather than trusting the script. (Recorded because
the same rewrite may apply to other partitions.)

For context on why this was attractive: at 22:50 `broadwl` had **224 pending jobs, 220 outranking
ours**, 173 nodes `mix`, zero idle; `broadwl-lc` had 18 idle nodes and 504 free cores.

---

## 0. Connect first (needs a Duo push on the user's phone)

Key-based auth is NOT enabled; password + Duo is the only method. The ControlMaster socket expires
roughly hourly, so expect to redo this most sessions.

**midway2** (POPE/POPG REMD campaign): Direct connection confirmed working 2026-09-12 -- the IP block
that was present through 2026-09-10 has been lifted. Use `mdw2_master.exp` directly:
```bash
ssh -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<command>'
```
If the IP block returns and direct access fails, tunnel through midway3 as a fallback:
```bash
# 1. port-forward midway2:22 to localhost:2222 over the existing midway3 master (no Duo)
ssh -o BatchMode=yes -f -N -L 2222:128.135.112.69:22 \
    -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu
# 2. open the midway2 master over that forward           # USER MUST APPROVE DUO
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_via_tunnel.exp
```

**Key-based ssh is NOT accepted by midway2. Do not try again.** (Tested 2026-09-17, and
confirmed by the user.) The public key was installed in `~/.ssh/authorized_keys` on midway2 with
correct 600/700 permissions, and a key-only connection was still refused:

```
debug1: Offering public key: ~/.ssh/midway3 RSA SHA256:7gJxuwcx...
debug1: Authentications that can continue: publickey,gssapi-keyex,...,password,keyboard-interactive
debug2: we did not send a packet, disable method
```

Retested with `PubkeyAcceptedAlgorithms=+ssh-rsa` in case the legacy SHA-1 signature was the
blocker; still refused. The server advertises `publickey` but does not honour `authorized_keys`
for this account, so **every new connection costs a Duo push**. An inert key entry was left in
`~/.ssh/authorized_keys` on midway2; harmless, remove it if it ever bothers you.

Consequences, and the mitigations that actually work:

* **Minimise connections rather than trying to remove Duo.** The cluster-side monitor
  (`/project/trsosnic/yinhan/monitor.sh`, run by `monitor_loop.sbatch` on one `broadwl` core)
  refreshes `/project/trsosnic/yinhan/STATUS.md` every 30 min. It only reports; the training chain
  insures itself, and a dead chain shows as **CHAIN DOWN** in STATUS.md. Nothing depends on a live laptop connection, so a monitoring loop can run every few
  hours instead of hourly.
* **`scrontab` is disabled and `crontab` is denied** on this cluster, and `pi-trsosnic` has **no
  association with the `cron` partition** (`sbatch -p cron` is refused; the old monitor, believed
  to be on cron, actually ran on `broadwl`). A 1-core `broadwl` job is the way to schedule recurring
  work. A 7-day request hits `QOSMaxWallDurationPerJobLimit`; 36 h works, **provided the successor
  is queued at the start** (`--dependency=afterany:$SLURM_JOB_ID`). The first monitor resubmitted at
  the end of its loop, overran the wall and died silently on 2026-09-19.
* **Disabling laptop sleep is not sufficient.** With `caffeinate -is` holding
  `PreventSystemSleep`, the master still died after ~40 min, so network blips rather than sleep
  are tearing it down. `mdw2_master.exp` now uses `ServerAliveInterval=30
  ServerAliveCountMax=20 TCPKeepAlive=yes` (~10 min tolerance, was 3 min).
* **Run the expect script at most once per attempt.** Retrying it while diagnosing sends a second
  Duo push to the user's phone, which is intrusive and was done by mistake on 2026-09-18.

**`-o BatchMode=yes` is NOT sufficient on its own — it only protects when the socket FILE is gone.**
(Learned the hard way 2026-09-19.) If the master dies but `~/.ssh/cm-mdw2.sock` is still on disk,
ssh tries the stale socket, fails, and then **falls through to a real connection attempt**, spending
an authentication attempt and printing `Permission denied (publickey,...)`. A few of those and the
host starts answering `Connection closed by 128.135.112.69 port 22`, which is the IP throttle, and
the next `mdw2_hold.exp` launch dies with `MASTER_DIED_EARLY` **after** it has already sent a Duo
push. That is how a single unguarded status call costs a push and locks you out for tens of minutes.

**So: run `-O check` FIRST and only issue the real command if it succeeds.** Never send a bare
`ssh -S sock host 'cmd'` on the assumption BatchMode will catch a dead master:

```bash
if ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu >/dev/null 2>&1; then
    ssh -o BatchMode=yes -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<cmd>'
else
    echo "master down - do NOT retry blindly; see throttle note"
fi
```

**When throttled, STOP.** The throttle clears on its own in tens of minutes. Wait at least 30 min
before one single retry. Do not relaunch the hold script to "see if it works" — each launch spends
a Duo push on the user's phone before it discovers the throttle.

**ALWAYS put `-o BatchMode=yes` on routine ssh calls.** (Learned 2026-09-17.) A plain
`ssh -S ~/.ssh/cm-mdw2.sock host 'cmd'` does **not** fail when the socket is dead: it silently
falls back to a fresh connection, offers the key, then tries password auth twice non-interactively,
hits `Received disconnect ... Too many authentication failures`, and after a couple of those the
host starts closing connections immediately (`Connection closed by 128.135.112.69 port 22`). That
is the IP throttle described below, self-inflicted by a status check. With `BatchMode=yes` the same
call fails instantly and harmlessly with `Control socket connect: No such file or directory`, which
is the signal to run the expect script. Check the socket first, and never let a monitoring loop
retry the expect script more than once per tick.

The throttle appears to clear on its own in tens of minutes. The documented way around it while it
lasts is the midway3 tunnel (`scratchpad/mdw2_via_mdw3.exp`), which needs a live
`~/.ssh/cm-mdw3.sock` and therefore its own Duo approval.

Two zsh/tooling traps that cost time here:
* **zsh does not word-split unquoted variables.** `M2="ssh -S sock host"; $M2 'cmd'` runs silently
  and produces NOTHING, it is not a connection failure. Write the `ssh` call out in full.
* `timeout` does not exist on this Mac; do not wrap `expect` in it.

`~/project` on midway2 is a symlink to `/project/trsosnic` (not `/project/trsosnic/yinhan/`).
Data is at `~/project/yinhan/popepopg_REMD_mdw2/`.

**`/project` is the SAME filesystem on midway2 and midway3** (both mount `midway3_cap`), so
`/project/trsosnic/yinhan/upside2-md-mdw2/` is fully readable from midway3. Checkpoint progress,
chain logs, `.ff_installed`/`.production_relaunched` flags and force-field directories can all be
checked **without a midway2 login at all**, useful when the IP block or a login outage is in the
way. Only `squeue`/`sacct` need midway2 itself; the two clusters have separate Slurm controllers and
separate accounting databases, so `sacct --clusters=all` on midway3 does **not** see midway2 jobs.
Python env: `source /software/modules/init/bash && module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1 && export HDF5_USE_FILE_LOCKING=FALSE`

**midway3** (NP campaign):
```bash
ssh -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu '<command>'
```
`~/project` on midway3 is a symlink to `/project/trsosnic/yinhan/` (note: **yinhan**, not yinhanw).
Load the python env with `source ~/project/NP-1AO6/env.sh` before any h5py work.

**rockfish** (JHU ARCH — added 2026-09-07 as the outage-proof host for ff3.0 training):
```bash
ssh -o BatchMode=yes -o ConnectTimeout=25 rockfish '<command>'   # key auth, NO Duo
# ~/.ssh/rockfish, user ywang268. Strip the banner by emitting your own marker first:
ssh -o BatchMode=yes rockfish 'echo === ; <command>' | sed -n '/^===/,$p'
# A BINARY file needs base64, not cat: the banner is on stdout and prepends 1659 bytes.
ssh -o BatchMode=yes rockfish 'echo __B64__; base64 <file>' 2>/dev/null \
  | sed -n '/^__B64__$/,$p' | tail -n +2 | base64 -d > local_file   # then check md5 both sides
```
Key-based, so this host needs no interactive second factor and can be polled freely. The login
banner (survey notice + quota tables) is **not** a clean shell: it breaks `rsync` outright
(`protocol incompatibility`) and `scp` with `Received message too long`, so transfer with
`tar czf - … | ssh rockfish "tar xzf - -C …"` or `cat file | ssh rockfish "cat > dest"`. Host key was
verified over two independent network paths (this Mac and midway3) before being recorded, since ARCH
publishes no fingerprints:
`ED25519 SHA256:V58d1zhfocFT/JR90J3HqMw6uJTLhn+Nc58NfnoJqM8`.

Paths: repo `/scratch4/rherna21/ywang268/upside2-md-rf` (`$RF`), scratch4 group quota 15 TB with
4.4 TB free. `module` needs `source /etc/profile.d/modules.sh` first — `lmod.sh` does not exist and
without it `module load` silently does nothing.

---

## 0b. Shared Upside deployment on beagle3 (2026-09-09)

**One tree, one binary, both clusters: `/beagle3/trsosnic/yinhan/upside2-md`.**
`source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh` sets everything up identically on
midway2 and midway3; verified by running the binary on each.

Why this works, all measured rather than assumed:

* **`/beagle3` is visible and writable from the COMPUTE NODES of both clusters.** A probe job on
  midway2-0291 and midway3-0036 confirmed it, along with `/project` and `/project2`.
  **`/cds3` is login-node only** and cannot be used for jobs.
* **One binary serves both.** It was compiled on midway2's Broadwell with `-march=native`
  (`src/CMakeLists_Other.txt:8`) and contains **zero AVX-512** (`zmm` register uses = 0), so it runs
  on midway3's Cascade Lake. **Compiling on midway3 instead would be a trap**: `-march=native` there
  emits AVX-512 and the binary would die with SIGILL on broadwl.
* **Module requirements differ and the env file handles both.** `hdf5/1.14.3+oneapi-2023.1` exists
  on both and supplies the `libhdf5.so.310` this binary links against. midway2 *additionally* needs
  `gcc/10.1.0`, because its system libstdc++ lacks `GLIBCXX_3.4.20` and the binary will not load
  without it; midway3 does not need gcc.

Deployed by rsync from `/project/trsosnic/yinhan/upside2-md-mdw2` (430 MB: `src py obj parameters
cmake example` and the install scripts), **excluding `training/`**, which is 6.6 GB of ConDiv
campaign data that does not belong in a code deployment. `parameters/ff_3.0_trained/sidechain.h5`
verified `c67351ca...` on beagle3.

**Filesystem headroom, and a trap.** `df` is misleading on `/project2`: it shows 787 T free while
the trsosnic *group* quota there allows only **195 G** more. By group headroom the usable
filesystems are `/project` (1516 G), `/beagle3` (1434 G), `/cds3` (789 G, unusable from compute
nodes) and `/project2` (195 G). beagle3 was chosen so benchmark output does not compete with the
live glpG campaign on `/project`.

**Naming:** the cluster keeps `ff_3.0_trained` / `ff_3.0_trained_rf` because the running chain
scripts reference those paths. The repo uses the plain `ff_3.0` slot. Do not rename on the cluster
while the arm test is running.

### Other Upside copies (inventory 2026-09-09)

Owned by yinhanw:

| path | size | branch | last touched | disposition |
|---|---|---|---|---|
| `/project/trsosnic/yinhan/upside2-md-mdw2` | 8.3 G | martini-dev | today | **ACTIVE**, all 3 running jobs use it |
| `/beagle3/trsosnic/yinhan/upside2-md` | 6.4 G | martini-dev | today | **the shared deployment**, refreshed from stale 2026-08-27 |
| `/scratch/midway2/yinhanw/upside2-md-water-diff` | 13 G | water-diff | 2025-10-30 | idle; see below |
| `/scratch/midway2/yinhanw/upside2-md` | 176 M | **master** | 2026-01-28 | keep, master |
| `/home/yinhanw/upside2-md` | 256 M | **master** | 2025-11-19 | keep, master |

**Not ours, never touch:** `/beagle3/trsosnic/upside2-md` and `/project2/trsosnic/upside2-md`
(bayhi), `/project2/trsosnic/software/upside2-md` (nffaruk), `/project2/trsosnic/pengxd/*`
(pengxd), plus copies under baxa, avmolina, ruofan, tobin, schwartznw, yiheng, zonganw,
simoneritchey, amz and bayhi.

**`upside2-md-water-diff` is the only deletion candidate and it is NOT purely a code copy.** Its
source is safe: HEAD `c44a404` is reachable from `origin/water-diff` on GitHub with 0 uncommitted
files. But **11 of its 13 GB is `example/16.MARTINI/outputs/water_T*`**, water-diffusion simulation
output from 2025-10-30 that is not in git and exists nowhere else. Deleting the directory discards
that data. Left in place pending a decision.

## 1. Current jobs

Snapshot **2026-09-24 21:50 CDT, verified live against `squeue`.** Finished and cancelled rows are
deleted; only lessons worth reusing are kept, below the table. Phase 1 extension (49073314) completed
2026-09-24 in 2:06; its successor 49073315 was cancelled as planned.

| JobID | what | where / state | next action |
|---|---|---|---|
| **49074120** | **Phase 2: ff3.0 from ff2.1**, `training/ff30`, 76 steps, `TRAIN_GLY = True` | R 9:43 of 36 h, 13 nodes; 26/76 steps done (`epoch_01_minibatch_07` running), ~22 min/step, no errors in log | step 76 ~**2026-09-25 16:00 CDT**; then `after_training.sbatch` -> `gate_or_continue.sh`: converged -> `validate_ff.sh ff_3.0`; not -> +1 epoch (~7 h), up to 13 |
| 49074122 | its insurance successor (`train_chain.sbatch`) | PD, `afterany:49074120` | resumes the chain only if 49074120 dies before step 76 |
| 49082029 | `ff30_monitor`, `/project/trsosnic/yinhan/monitor_loop.sbatch` -> `monitor.sh` | R on midway2-0461, 36 h | writes `STATUS.md` every 30 min (queue, step, errors, GLY-row line from `training/ff30/analysis/gly_status.py`); log `monitor_loop.out` |
| 49082147 | its successor | PD, `afterany:49082029` | takes over at the wall; to stop the monitor cancel **both** |

Log `/project/trsosnic/yinhan/upside2-md-mdw2/training/ff30/condiv-train_49074120.out`; checkpoints
`.../training/ff30/run_output/epoch_EE_minibatch_MM`. Monitor rewritten for Phase 2 and restarted
2026-09-24 21:48; the ff3.0C/AWH version is kept as `monitor.sh.bak_ff30C_awh`.

**dhb diagnosis (steps 20-24), settled.** The native-state and 0.3 x unfolded-state gradients
nearly cancel at ff2.1 for every H-bond parameter (dhb: NSE mean +19.8, lambda*DSE mean -15.3, net
+4.5 +/- 9 changing sign; backbone scale -114 vs +127), which is what ff2.1 having been trained
with lambda = 0.3 predicts. The reweighting fix reduces rather than causes the pull (the port's
weights add ~+1.8 each step). dhb is a mildly unconverged ff2.1 parameter, not a port error.
**Glycine gate passed on midway2** (job 49073372, 6.6 min) before phase 2 was initialised.

Step 20 lost 4pqz to `srun: Invalid job credential` on the first launch at job start, a Slurm
race; every later step ran 24/24. If it recurs at every link start, add a settle delay before the
first launch.

**Phase 1 (19 steps, job 49056803) COMPLETED 0:0 at 09:15**; report in
`fixedpoint_report_19steps.txt`. Every trained file updates; 8 of 9 groups at a fixed point; **dhb
(second-H-bond term) is not**: -0.406 -> -0.448, t = +2.98, still drifting. The hand-off worked:
target reached, successor cancelled, report job 49073246 submitted and completed in 37 s.

**`run_output/ConDiv.py` in this run carries TEMPORARY DIAGNOSTIC code** (marked `DIAGNOSTIC`); the
clean copy is `run_output/ConDiv.py.clean` (md5 = repo `training/ConDiv.py`). Restore it, or
discard this run, before this directory is used for anything else. Phase 2 initialises a fresh run
from the repo trainer, so it does not inherit the diagnostic.

**Nothing else is queued.** The 32 Peng arms and 4 glpG chains of the FF1-form ff3.0 were cancelled
2026-09-24 ~01:45 at the user's request; that force field is superseded (findings 9t-9v).

### The FF2 trainer on midway2 (2026-09-24)

`training/ConDiv.py` is now the FF2 dual-target trainer (findings 9v, plan.md). Deployed md5-matched
to `$P/training/`, with `extract_ff.py`, `patch_glpg.py`, `validate_ff.sh`, `check_converged.py`,
`train_chain.sbatch`; the old trainer and helpers are in `$P/backup_training_ff1form/`.

**Measured step cost, full protocol** (`worker_test/`, job 49056799): 5vhg, 150 residues, the
largest in the set, 1242 s; 1ean, 114 residues, 897 s. A step waits for its slowest worker, so
**~21-22 min per step**; `STEPS_PER_LINK = 90` fits a 36 h wall. Both unfolded properly on the SI's
ladder (mean Rg 14 -> 46 A and 13 -> 36 A across T = 0.8-1.1), so the DSE target is real.

**Schedule (estimates):**

| stage | length | expected end |
|---|---|---|
| Phase 1, 19 steps from ff2.1 | ~7 h | 2026-09-24 ~09:30 CDT |
| read the report; glycine gate on midway2; initialise phase 2 | ~1-2 h | 09-24 ~12:00 |
| Phase 2, 76 steps, `TRAIN_GLY = True`, ff2.1 -> ff3.0 | ~28 h + queue | **2026-09-25 ~18:00-24:00** |
| validation, auto-submitted: 32 Peng arms | ~7 days (longest arms ~51 k time units/h) | ~10-02 |
| validation, auto-submitted: 4 glpG REMD chains, 5 x 36 h blocks | ~7.5 days | ~10-03 |

**Phase 2's hand-off.** Its run directory gets an `after_training.sbatch` that runs
`bash $P/training/validate_ff.sh "$SLURM_SUBMIT_DIR" ff_3.0`: extract through the run's own
`expand_param`, back up and overwrite `parameters/ff_3.0` in `$P` and in the /beagle3 deployment
(md5-verified), move the superseded `runs/*_ff_3.0` benchmark directories to `runs_superseded/`,
submit the 32 arms, gate the glpG patch on a pristine ff_2.1 seed, patch the 4 live seeds, clear
their replicas and submit the 4 chains. `train_chain.sbatch` submits it itself when the target is
reached and cancels its insurance successor, so no 336-CPU job has to queue just to call sbatch.

**`bench_run.py` changed 2026-09-24** (backup `bench_run.py.bak_ff1form_20260924`): the type-0
burial override for ff_3.0 and the `rama3.dat` fallback are gone, since the new ff3.0 is FF2-form.
Nothing may benchmark the old FF1-form ff_3.0 with it.

### The 2026-09-23 training failure: an infrastructure kill, not a defect

Link **49047139** ran 42 minibatches cleanly (steps 449 -> 491) and then **FAILED with exit 7**
8 h 36 m in, far inside its 36 h wall. What the evidence says:

* All 12 workers of minibatch 35 were killed by **signal 7 (SIGBUS)**, `sacct` steps `.492`-`.503`
  all `CANCELLED 0:7`, spread across **four different nodes** (midway2-[0276-0279]). A bad node
  kills 3 tasks, not 12, so this is not node-local and the nodes were **not** added to `--exclude`.
* The physics was healthy at the moment of death: every worker was near frame 1985/4000 with
  Rg 14.5 A, ~110 hbonds and potential around -200. Not a blow-up.
* **No traceback anywhere**, and all 12 `*.output_worker` files are 0 bytes. The job's own `.out`
  stopped being written at 11:22 while the workers kept writing until 11:26-11:31, and the job was
  not reaped until 12:15.
* **Disk was not the cause, checked directly.** `/project` trsosnic fileset: 3.5 T of 3.9 T, 445 G
  free, inodes 216 K of 1.1 M (20%). A 200 MB write+delete on `/project` succeeded at 1.7 GB/s.
  `/project2` group is the tight one at 1.45 T of a 1.49 T soft quota (97%), but nothing in this
  campaign writes there. `rcchelp quota` is the tool that reports the group numbers; plain `df` on
  the mount point shows the whole 6.3 P filesystem and tells you nothing, while `df` on the
  **subdirectory** does report the fileset.

* **The successor settles it.** The chain queues its replacement with `--dependency=afterany`, so
  this should have been survivable. 49053769 started 12:18:02 and was `CANCELLED` the same second
  with zero elapsed, and its `.out` file was never created: it died before it could open its own
  output file. Those nodes could neither read nor create files on `/project` between 11:26 and
  12:18.

SIGBUS on mmap'd HDF5 across four nodes, plus a batch step that cannot create its output file
52 minutes later, is a transient `/project` outage. **No GPFS log was available**, so this is
inferred from symptoms rather than confirmed at the source. The chain logic itself is sound and
needed no change; recovery was to resume.

### The validation handoff was broken and would have fired nothing (found and fixed 2026-09-23)

`train_gly.sbatch` ended its "target met" branch with `sbatch "$T/validate_ff31.sbatch"`, but `$T`
is the **run directory** `training/ff31-gly/` and the script has always lived one level up in
`training/`. The `||` fallback would have printed "WARNING: validation was NOT submitted" into a
log nobody was watching, and the 32 Peng arms and 4 glpG chains would simply never have been
queued. Fixed by defining `UP=` and calling `$UP/training/validate_ff31.sbatch`.

**Why it survived this long: that branch runs exactly once, at the very end of a four-day chain.**
Every other line of the script had been exercised twelve times. A code path that only executes on
success, at the end, is untested by construction; check it by hand before the run that will use it.

**And the fix alone was not enough, which is the reusable part.** Slurm **snapshots the batch
script at submit time**. The successor 49056522 had already been queued by the running job, so it
still carried the old path and editing the file on disk changed nothing for it. Dump what a queued
job will actually execute:

```
scontrol write batch_script <jobid> /tmp/js.sh && grep -n <the-thing-you-fixed> /tmp/js.sh
```

That showed 49056522 still on `$T/validate_ff31.sbatch`, so it was cancelled and resubmitted as
49056550 with `--dependency=afterany:49056521`. Cancelling a **pending successor** while the
running link continues is the safe direction; the dangerous one, cancelling the running link first
and leaving the successor to start on the same `run_output`, is what produced two concurrent
writers earlier in this campaign.

**Pre-flight check that was run after the fix**, all present: `bench.sbatch` and `logs/` under
`$B`, `submit_remd.sh` under `$GM`, `ff_2.1/bb_env.dat`, the pristine `glpG-RKRK-79HIS` handoff
seed, all four live glpG seeds, and `extract_ff31.py`/`patch_glpg_ff31.py`/`rama_gly_gradient.py`
md5-identical to the repo copies. `parameters/ff_3.1_trained` correctly does not exist yet.

### Two things fixed during the recovery

**`train_gly.sbatch` overshot its own target.** The `STEP >= TARGET` test only runs at link start,
and the link then always ran `STEPS_PER_LINK=150`. Resuming at 491 of a 500 target would have run
to 641, roughly 29 h of pointless training before validation could fire. The last link now runs
`TARGET - STEP` steps. Backup at `train_gly.sbatch.bak_overshoot`.

**The aborted minibatch directory had to go before resuming.** `run_output/epoch_12_minibatch_35`
held 13 GB of half-written `.h5` and no `checkpoint.pkl` and no `divergence.pkl`. Deleting it took
`run_output` from 19 G to 6.1 G. Leaving a partial minibatch directory in place is a real hazard:
`main_worker` reads `<name>.divergence.pkl` by path, so a stale one from an earlier attempt would
be picked up as a fresh result for a worker that failed.

**Pending cleanup, safe to defer.** `$P/py/rama_gly_gradient.py` on the cluster is a hardlink to
`$P/training/rama_gly_gradient.py`, left from the move out of `py/`. One inode, two names, so the
two cannot drift; `ff31-gly/env.sh` now carries both directories on `PYTHONPATH`, so the `py/`
name can be deleted whenever nothing is running.

**Validation is armed and will submit itself.** When `train_gly.sbatch` sees step >= 500 it
runs `training/validate_ff31.sbatch`, which extracts `parameters/ff_3.1_trained` and submits 32
benchmark arms (16 Peng proteins x native/denovo) to `broadwl`, each self-chaining via
`bench.sbatch`. Nothing needs doing by hand. If the log line
"WARNING: validation was NOT submitted" appears, run `validate_ff31.sbatch` manually.

**`parameters/ff_3.1_trained` does not exist yet, deliberately.** A dry-run copy was created to
test the pipeline and then deleted, so that nothing can benchmark a mid-training force field by
accident. `validate_ff31.sbatch` creates it from the final checkpoint.

**Disk headroom for the glpG half of the validation, measured 2026-09-23.** Each of the four
variant directories in `popepopg_REMD_mdw2` holds ~150 GB of ff3.0 replicas, 604 GB for the four,
against only 445 G free on `/project`. This is safe only because `validate_ff31.sbatch` does
`rm -rf "$GM/$V"` immediately before submitting that variant's replacement, so the space is
returned ahead of the write. Do not reorder that loop.

### glpG for ff3.1: patched into the midway2 tree, submits itself

`validate_ff31.sbatch` patches the trained force field into the four seeds in
**`popepopg_REMD_mdw2`** and submits all four REMD chains on **broadwl**, 28 replicas each
(`REMD_N` follows `--cpus-per-task=28`), `REMD_T_HI=0.90`, `REMD_MAX_BLOCKS=5`. No midway3 step.

**An earlier version of this section was wrong.** It said glpG needed `caslake` and therefore
midway3. That came from reading only `popepopg_REMD/` (48 replicas, caslake) and concluding from
`sinfo` that midway2 could not host it. **`popepopg_REMD_mdw2` already exists and has run all
four variants to COMPLETED on broadwl** (jobs 48890657-60, ~1d10h each). Check job history before
concluding a cluster cannot run something.

**Seeds are the LIVE ones, not the pristine ff2.1 backups**, because only they carry
`inner_steps = 4` on `/input/brownian`. That temperature fix is worth TM4 helix fraction 0.893
against 0.346 (findings 3.10a). It is an HDF5 *attribute*, so patching arrays preserves it;
patching the pristine backup instead would silently throw it away.

**Measured scope.** The live seed differs from the pristine ff_2.1 seed in exactly two arrays,
`rama_map_pot/rama_pot` and `rotamer/.../interaction_param`, which carry ff3.0. Those are
precisely what the patch overwrites, and ff3.0 is superseded. Patching touches exactly three
arrays (those two plus `hbond_energy/parameters`) and nothing else; verified array-by-array.

| check | result |
|---|---|
| method gate: round trip ff_2.1 into a **pristine** seed | **3.55e-15** |
| round trip into a **live** seed | fails at 3.26, **as it must** -- live seeds carry ff3.0 |
| arrays changed by the patch | exactly 3, the intended ones |
| `inner_steps` after patching | **4**, preserved |
| engine, live ff3.0 seed | -24957.682 finite, rama 56.338, hbond -408.530 |
| engine, ff3.1 patched | -24965.998 finite, rama 68.505, hbond -425.492 |

The gate runs on a pristine seed because that is the only place a round trip *can* succeed; a
failure on a live seed would carry no information. Method proven once, then applied.

Seeds are backed up to `*.bak_pre_ff31_<stamp>` before patching, and a failed patch restores
them. Stale ff3.0 replica directories are removed so `block_count` restarts.

**Criterion: helix stability over time, TM4 above all.**

**Watching Track A.** `grep '^gly' ff31gly_*.out` prints `dG(aR->aL)`, `|A| rms` and the
`GLY|GLY` asymmetry every step. The handedness starts at exactly 0 and the number to compare
against is the AWH's **-0.303 nats**. **`GLY|GLY asymmetry` must stay `0.00e+00`**; if it ever
moves, the symmetric re-projection in `backprop_deriv` has broken and the run is invalid.
Progress: `find run_output -name checkpoint.pkl -path '*epoch_*' | wc -l` against 500.

**Also watch `hb`.** It is trained, unconstrained, and starts at 1.0. Adam's first step is
scale-invariant, so early movement of exactly +/-0.00255 per step is just `alpha` and means
nothing. Sustained drift does: 500 steps in one direction reaches **-0.27** and inverts the sign
of every hydrogen bond. **Outside roughly 0.8-1.2, stop and diagnose.** Do not add a clamp; the
original trainer has none and a runaway is real information about the model.

The earlier chain on the measured map (49037907/08) was stopped: that map is now the Track B
*reference*, and Track A learns its own from a symmetric start so the comparison is not circular.

**`training/ff31-gly/` starts from the ff2.1 original `rama.dat`, NOT `rama31.dat`.** Verified by
md5. Seeding Track A from the measured map would make the whole comparison circular, and it is a
one-character mistake to make when cloning a run directory.

**`run_output`: 13 MB per completed step, so ~6.5 GB over 500, but a transient peak of ~8.4 GB.**
Measured. A completed minibatch directory holds only the parameter files; the in-progress one also
holds 12 proteins x 8 replica `.h5` trajectories (158 MB each for the largest protein), which each
worker deletes on success. **So `du` on `run_output` mid-step reads ~9 GB and looks alarming; it
is not cumulative.** Check `du -sm epoch_*` and compare a finished directory against the live one
before concluding anything. The 35 MB rama library each step writes is also deleted once its
workers exit.

**`training/ff31-gly/` starts from the ff2.1 original `rama.dat`, NOT `rama31.dat`.** Verified by
md5. Seeding Track A from the measured map would make the whole comparison circular, and it is a
one-character mistake to make when cloning a run directory.

### Disk: the AWH campaign was going to cost ~100 GB, and most of it was unread output

`gly_peptides` was 14 GB for 20 systems at 100 ns, i.e. **645 MB per system per 100 ns**. Scaled to
40 systems at 400 ns that is **~100 GB**. Broken down for one system:

| file | size / 100 ns | read by anything? |
|---|---|---|
| `awh.part0001.edr` | 96 MB | yes, but only its AWH frames |
| `awh.part0001_pullx.xvg` | 27 MB | **no** |
| `awh.part0001_pullf.xvg` | 27 MB | **no** |
| `awh.part0001.xtc` | 2 MB | no, but cheap |

Two changes to `awh_template.mdp`, both output-only, neither touching the dynamics:
* **`nstenergy = 5000 -> 50000`**, matching `awh-nstout`. AWH data is written as part of energy
  frames and `awh-nstout` must be a multiple of `nstenergy`, so equality keeps every AWH frame and
  drops the 10x redundant plain-energy frames.
* **`pull-nstxout = 0`, `pull-nstfout = 0`.** The pull code still runs and still drives the AWH
  bias; only its output files are suppressed. The PMF comes from `gmx awh` on the `.edr`.

**Measured after the change, not projected**: the new systems write **0.93 MB/ns** against the
extensions' 1.62 MB/ns, a 43% cut. Verified in the tpr (`gmx check -e` reports 16 frames at a
100 ps timestep, matching `nstenergy = awh-nstout = 50000`, and no pull files exist).

Campaign cost is therefore about **21 GB** on top of the current 12 GB: 30 new systems x 400 ns at
0.93 MB/ns is 11 GB, and 20 extending systems x 300 more ns at 1.62 MB/ns is 10 GB. That is small
against `/project`, where the retired NP campaign alone holds ~500 GB, so **no further action is
warranted**.

**What the remaining bytes are, in case it ever does matter.** At 90 KB per energy frame for a
14-atom peptide, the `.edr` is almost entirely the AWH 2D grid, not the plain energy terms. So the
`nstenergy` change mattered much less than removing the pull files did, and the only lever left is
`awh-nstout`. Raising it 50000 -> 500000 would cut the dominant term 10x and still give 400
convergence snapshots over 400 ns, which is far more than the ~10 ns spacing any analysis has
used. It was not done: it needs a third rebuild of all 30 systems to save ~10 GB.

The three new groups were cancelled and resubmitted so they build under the new template; the 20
already-running systems keep their large `part0001` files, and their closed `part0001_pull?.xvg`
were deleted for an immediate **1.2 GB** back. Their output settings **cannot** be changed
mid-flight: `grompp -t` would reset the AWH bias history and throw away 100 ns of learning, and
`convert-tpr` cannot alter output frequency. Backup: `awh_template.mdp.bak_verbose`.

**The RCC quota tools do not work from midway2**: `mmlsquota` reports "File system project is not
known to the GPFS cluster" for project, project2 and beagle3, and `rcchelp quota` dies with a
`TypeError` inside `/project2/rcc/rupat/bin/quota.py`. `df` shows `/project` 77% used and
`/beagle3` 64%, but those are whole-filesystem numbers, not the group quota. Until a working tool
turns up, track our own footprint with `du` rather than trusting either.

**ff3.1 training chain details.** ff_2.1 init params plus `rama31.dat`, whose coil GLY row is the
AWH-measured dipeptide surface and holds **no library data**, trained with the
strict-modernization ConDiv with `hb` and `sheet` unfrozen. The successor is queued BEFORE
training starts, so a wall kill cannot silently end the chain. Progress:
`find run_output -name checkpoint.pkl -path '*epoch_*' | wc -l` against 500. Rebuild the map with
`py/build_rama_from_awh.py` when the 400 ns extensions land, and restart if the row moves.
Superseded attempts, all cancelled when the map construction changed: 49037796/97/804 (uniform
lambda), 49037810/11 (artifact subtraction), 49037815/16 (`S_library + A_measured`, rejected
because `S_library` is ff3.0), 49037514. The 49037815 output is kept as `run_output.bak_libsym`.
See findings.md 9i and 9l.

**Cancelling a self-chaining job: cancel the PENDING successor FIRST, then the running link.**
Doing it the other way spawns a chain you do not know about. `scancel` on the running link
immediately satisfies the successor's `afterany` dependency, so the successor starts before the
second `scancel` in the same command lands, and it queues a successor of its own. That happened on
2026-09-19: two chains ended up writing the same `run_output/epoch_00_minibatch_00` concurrently.
Both were killed, the partial epoch directories deleted, and one chain resubmitted. Check
`squeue` after any chain cancellation rather than assuming it worked.

**`awh_extend.sbatch` does not self-chain; `awh_batch.sbatch` does.** The two extension jobs were
submitted with the former and would have stopped dead at their 36 h wall around 200 ns, with no
successor and no error. Continuations are now queued on `afterany` (49038531, 49038532) using
`awh_batch.sbatch`, which takes an optional third argument for the replica subdirectory so it can
drive `awh_amber99sb-ildn_rep2` as well as rep1. **Anything submitted with `awh_extend.sbatch`
needs a successor queued by hand.**

**midway2-0080 added to the exclude list** after a NODE_FAIL took 49037918 down at 5:07. The chain
self-healed: 49037921 resumed from `awh.cpt`, and the affected systems show both
`awh.part0001.log` and `awh.part0002.log`, so nothing restarted from zero. The list is now
`midway2-0003,midway2-[0010-0011],midway2-0080,midway2-[0342-0345]`.

**Reading AWH progress: glob `awh.part*.log`, never `awh*.log`.** `gmx awh` leaves an
`awhtool.log` in the system directory, it sorts last under `ls -v`, and it contains no
`Step  Time` records. A progress script that takes the last `awh*.log` therefore reports **0 ns**
for any system that has ever been analysed, which on 2026-09-19 made the `LG` blank look dead when
it was running normally at 119 ns.

**AWH extension details, and how to analyse it.** `awh_extend.sbatch REPLICA_DIR [TARGET_PS]` in
`/project/trsosnic/yinhan/gly_peptides/` runs `convert-tpr -until` then `mdrun -cpi -noappend` for
all 10 dipeptides. **`gmx awh` must read the LAST part file**, `ls -v awh.part*.edr | tail -1`; the
AWH state is cumulative, so the final part carries the whole PMF and reading `awh.edr` silently
gives the old 100 ns answer.

**Why extend at all.** The achiral `GLY|GLY` blank is the convergence criterion, not a wall time.
It must read 0, reads rms 0.032 at 100 ns, and decays as `1/sqrt(t)` (0.233 at 10 ns) while the
signal converges to 0.071. 400 ns should halve it. Stop on the blank.

**Node exclusions for every job in `gly_peptides`:**
`midway2-0003,midway2-[0010-0011],midway2-[0342-0345]`. Three NODE_FAILs came from the adjacent
0010-0011 pair with `ExitCode 0:0`. AWH builds survive a NODE_FAIL, so a resubmit skips to grompp
and resumes from `awh.cpt`.

**Completed campaigns, kept for their conclusions only:**
* **Campaign 6 (AWH, 49033509 / 49033985)** delivered the answer, **-0.26 E_up**, and the finding
  that per-neighbour structure does not reproduce (Spearman rho +0.048). Both were cancelled once
  their 10 dipeptides had passed 100 ns; 49037819/20 continue from their checkpoints.
* **Campaign 7 (lambda diagnosis, 49035287/94/300)** is complete: the failure is in helix crossing
  angles, and no glycine map change rescues it (-0.024 E_up).
* **Campaign 5 (ff3.0C training, `training/gly-ctx`)** is cancelled and must not be resurrected;
  the AWH measurement contradicted its founding premise, and the library's per-pair ordering is
  anti-correlated with the measurement (r = -0.540).
* **`49032235 np_1AO6_prod`** cancelled 2026-09-18, the last job still simulating with ff3.0. Six
  `np.run.N.up` intact (~500 GB in `NP-1AO6/prod_ff3/`), resumable, and the obvious disk reclaim.
* **`49032512` pentapeptides** ran to its wall with no successor. Not usable as a cross-check:
  Gly5 read -0.224 against an exact 0.

**The first AWH attempt measured nothing, and was cancelled 2026-09-18.** `awh.mdp` had
`awh1-dimN-diffusion = 5e-5 rad^2/ps` while AWH's own friction metric implies `D ~ 0.77` (10-90%
range 0.33-3.4), about **15000x too small**, and `awh1-error-init = 5` kJ/mol against a true
surface range of 20-40. AWH left the initial stage at t = 13.34 ns holding ~1.0 kJ/mol of PMF
range and then crept up at 0.07 kJ/mol/ns, reaching 2.0 kJ/mol at 25.8 ns. The surface was
therefore near-flat, which forces every chirality observable to its achiral value, so the
"handedness is zero" result was an artifact rather than a measurement. Withdrawn in `findings.md`.
Evidence kept in `gly_peptides/evidence_diffusion_bug/`.

Relaunch carries `diffusion = 0.5` and `error-init = 30`. Neither enters the free-energy
estimator, so only the convergence rate changed. Backups: `awh.sbatch.bak_pre_diffusion`,
`awh14.sbatch.bak_pre_diffusion`. Build products (`equil.gro`, topologies, index files) were kept,
so both jobs skip straight to grompp; production files were cleared so nothing resumes from the
old AWH state (1.2 G -> 26 M).

**Replica 2 launched 2026-09-18 (49033985), and it is the error bar.** `awh_rep2.sbatch`, dir
`awh_amber99sb-ildn_rep2/`. The only achiral control, LG, **cannot certify a chiral value**: an
achiral surface has errors that cancel by symmetry, which is exactly the error mode afflicting the
chiral systems. LG held a 0.032 E_up spread over 24-33 ns while LA drifted monotonically -0.459 ->
-0.678 over 15-33 ns, and LP and LL drifted the opposite way toward zero. So between-replica spread
is the only valid uncertainty available, and no individual `dG` should be quoted without it.
Independence: identical topology and box, fresh 200 ps equilibration from replica 1's `min.gro`
with `gen-seed = 20260918`, then `gen-seed = 20260918` and `awh-seed = 776611` in production
(replica 1 used an auto-generated `awh-seed = 1034933739`). Everything else, including
`diffusion = 0.5` and `error-init = 30`, is copied from replica 1's own `awh.mdp`.

**Convergence test.** Watch `range` in `awh_an.py` output, the largest PMF over sampled cells.
Below 15 kJ/mol the row prints `UNCONV` and its `dG` means nothing. **Coverage fraction is not a
convergence test** and looked healthy all through the failed run.

**The fix is confirmed, and the answer changed.** At 15 ns, 49033509 shows range 50-66 kJ/mol and
100% coverage, against 2.0 kJ/mol at 25.8 ns before. The handedness is no longer zero. Achiral
controls: LG +0.005, GGGGG -0.016 E_up, so the noise floor is ~0.02. Chiral neighbours: LA -0.348,
LL -0.398, LR -0.267, LP -0.162, LE -0.092, LT +0.097, LV +0.053, LD -0.017, LM -0.045. Several are
far outside the control band, so **XGX asymmetry is real but smaller than the library's** (-0.43 to
-1.40), which would make ff3.0's exact zero wrong as well, just less wrong than ff2.1. Magnitudes
are provisional at 15 ns; time-stability has not been checked yet. **SAGAS is not an achiral
control** (it has L residues), so its +0.289 is a real signal and earlier briefs mislabelled it;
only LG and GGGGG are achiral.

**Tier 2 is back on.** The 40-system neighbour table (`Ac-X-Gly-NHMe` and `Ac-Gly-X-NHMe` for 20 X,
already built in `gly_peptides/xg/` and passing `pdb2gmx` in both force fields) is the route to a
measured replacement for the GLY row of the coil group, since a converged AWH PMF is `-ln P(phi,psi)`
in kT, which is what `dimer_pot[GLY, dir, X]` stores. For the library's native 5 deg grid the umbrella
needs `force-constant ~ 330` rather than 128; 128 gives GROMACS' 46-point grid and is kept for the
cheap handedness re-run.

**NDRD provenance, settled with the licensed data.** The user obtained all four releases. Our
`coil` group is **`NDRD_TCB`**, identified exactly (correlation 1.00000, max deviation 0.0000
against `GLY|ALL` both directions). Central-GLY `dG(aR->aL)` by variant: Conly -1.876,
Tonly -0.839, TCB -0.965, TCBIG -0.410. **No variant is near the measured zero**, and the purest
coil set is the *most* biased, which refutes the turn-occupancy explanation and kills the idea of
switching to `Conly`. Full detail and the withdrawn claims are in `findings.md`.
Licence note: the NDRD files may not be redistributed outside the lab group, so they were analysed
locally and never copied to the cluster.

#### THE SITE GROMACS ON MIDWAY2 IS UNUSABLE; BUILD YOUR OWN (2026-09-17)

Do not spend time on the `gromacs/*` modules. Measured:

* **`gromacs/2024.1` dies with SIGILL**, even `--version`. It was built for a newer SIMD, and
  **midway2 has no AVX-512 on any partition** (all Broadwell E5-2680v4/E5-2690v4 or older), so no
  partition can run it.
* **`gromacs/2021.1` and `2019.3`** need five stacked modules just to load
  (`intel/19.1.1 intelmpi/2019.up7+intel-19.1.1 cuda/11.5 mkl gcc/10.1.0` for 2021.1; `cuda` and
  `mkl` satisfy the CUDA build's link deps, `gcc/10.1.0` supplies `GLIBCXX_3.4.21`), and then their
  MPI-built tools crash anyway: `pdb2gmx` gives SIGSEGV under `mpirun` and **SIGFPE** on a second
  run, and under `srun` Intel MPI fails in `PMPI_Init_thread`/`MPIDU_bc_table_create`.

**Working build recipe** (verified 2026-09-17, `/project/trsosnic/yinhan/gmx_build/build_gmx.sbatch`):

```
module load cmake/3.26 gcc/10.1.0          # 3.11 default is too old for GROMACS 2024
cmake ../gromacs-2024.4 -DCMAKE_INSTALL_PREFIX=/project/trsosnic/yinhan/gmx2024 \
  -DGMX_MPI=OFF -DGMX_THREAD_MPI=ON -DGMX_SIMD=AVX2_256 -DGMX_GPU=OFF -DGMX_DOUBLE=OFF \
  -DGMX_BUILD_OWN_FFTW=ON -DGMX_BUILD_OWN_FFTW_URL=file://$B/fftw-3.3.8.tar.gz -DBUILD_TESTING=OFF
```

Gives `GROMACS 2024.4, mixed precision, thread_mpi, AVX2_256, fftw-3.3.8`. Thread-MPI is the point:
no MPI initialisation at all, so 12 independent `gmx mdrun -nt 2` processes just work, no `mpirun`
and no `srun` wrapping.

Four traps, all hit:

1. **GROMACS pins FFTW by MD5.** `GMX_BUILD_OWN_FFTW_URL` does not relax the check, so the tarball
   must be **fftw-3.3.8** (`8aac833c943d8e90d51b697b27d4384d`). 3.3.10 fails verification at 12%.
2. **Download the tarballs on the login node**, which has outbound internet; compute nodes do not.
3. **Do not `source GMXRC.bash` under `set -u`** -- it references an unbound `GMXLDLIB` and aborts
   the script, which is why job 49032104 reported FAILED after installing successfully. Set
   `LD_LIBRARY_PATH=$PREFIX/lib64:$PREFIX/lib` by hand instead.
4. **Run `pdb2gmx` from a clean working directory.** It scans the cwd for force-field directories
   and throws `filesystem error: status: Permission denied` on any socket it cannot stat; a Cursor
   ssh socket in `/tmp` was enough to make it exit 1.

#### THE BLOCKER ANYONE RETRAINING WILL HIT FIRST (found 2026-09-16)

**`obj/libupside.so` was rebuilt 2026-09-10 20:55; gly-sym finished 2026-09-09 17:15. No ConDiv
run works against the current library without a patch, and there is no backup of the old `.so`**
(only of the `upside` executable, `obj/upside.pre_tempfix_20260910`).

Symptom: every worker dies with
`ERROR: Wrong number of parameters, expected 760 but got 360` from
`engine_c_library.cpp:110`, raised as `RuntimeError: Unable to get param deriv`.

Cause: `NonlinearCoupling` holds `coeff` (20x18 = 360) and `weights` (400). When the config's
`number_independent_weights > 1` (it is 20, from `upside_config`'s `--environment-weights-number`
default) `get_param_deriv` returns **both**, 760 values. `ConDiv.py` asked for `coeff.shape` alone.

Fix, applied in `gly-ctx/ConDiv.py` with the original kept as `ConDiv.py.bak_pre_envderiv`: ask for
the full 760 and slice the first 360. `environment.cpp` writes the coeff derivatives to
`deriv[ctype*n_coeff + starting_bin + i]`, i.e. indices below 360, and appends the weights
derivatives after, so the slice recovers exactly the gradient ConDiv used before and comparability
with ff3.0 is preserved. `weights` is not trained and is byte-identical between ff_2.1 and ff_3.0.

**Also required before submitting any new training run:** `gly-ctx/preflight.py`, which builds the
per-residue Rama maps for all 456 training proteins against a candidate library and checks they are
finite and in range. It catches a malformed library in two minutes instead of after a Slurm failure.
A first ff3.0C library was rejected by it: subtracting an antisymmetric component is **not**
range-preserving (averaging is), and doing it to the sheet group, whose `GLY|GLY` antisymmetry is
66 E_up out of a 3.8-73.3 range, drove the maps to 0.5-96.3 E_up.

### Campaign 1: glpG production, post-fix, midway2 only (updated 2026-09-12)

| JobIDs | variant | state | block | note |
|---|---|---|---|---|
| 49006599 | 79HIS | RUNNING 5h43m | 7/9 | temp 0.998-1.036, no rollbacks |
| 49006600 | 79HIS_S115T | RUNNING 2h47m | 7/9 | temp 0.997-1.016, clean |
| 49006601 | 79ALA | RUNNING ~2m | 8/9 | just started |
| 49006602 | 79ALA_S115T | RUNNING 11h3m | 6/9 | temp healthy |
| 49010274-49010277 | all 4 | PENDING (afterany) | next block | the LAST block; see below |

**The chain is confirmed on ff_3.0, measured not inferred (2026-09-13).** The tables baked into
`glpG-RKRK-79HIS.run.0.up` match `parameters/ff_3.0/sidechain.h5` at a least-squares scale of
**1.000000** with a residual of 1.6e-6, against 1.107 and 21 for `ff_2.1`, so all 33 output blocks in
the file are post-retraining. This had to be checked because `run_remd.py:63-65` only copies a replica
from the seed when the file is absent, and a naive verbatim hash sweep reports **ff_2.1** (the
dry-MARTINI `martini.h5` was never retrained and `ff_3.0/martini.h5` does not exist, so the seed pulls
it from `ff_2.1` by design). Procedure in findings.md 6.6.

**glpG VTF re-delivered 2026-09-14** to
`~/Documents/2026/reports/GroupMeetings/0914/glpG_RKRK_79HIS_run0_remd.vtf` (3146 frames, 321 MB,
blocks 1-54 at stride 5, `output_previous_0` skipped because it is the frozen seed block, internal RMSD
0.000). It replaces the 2026-09-13 file, which showed the protein a full box length out of the bilayer
in 159 of 1822 frames -- a periodic-image fault in `martini_extract_vtf.py`, not in the trajectory
(findings 3.8). Rebuild it with
`python3 ~/project/yinhan/extract_glpg_vtf.py <variant> <replica> <out_dir>` on midway2, under
`module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1` with `HDF5_USE_FILE_LOCKING=FALSE`; it imports
`/project/trsosnic/yinhan/upside2-md-mdw2/py/martini_extract_vtf.py`, which now carries the fix
(pre-fix copy kept as `.bak_pre_pbcfix`).

Verified before handover: protein COM exactly 0 in all 3146 frames, **0 displaced frames**,
protein-lipid xy centroid separation mean 0.70 A / max 2.20 A, every declared bond under 10 A except
the known residue-210 C-O. The earlier physical checks still stand: TM4 alpha fraction mean 0.791 with
no frame below 0.50, TM4 centroid within 7.7 A of the bilayer midplane, peptide C-N mean 1.324 A.

**THE CHAIN ENDS AFTER THE PENDING DEPENDENTS, verified against the files on 2026-09-13.**
`submit_remd.sh:33` exports `REMD_MAX_BLOCKS=5`, and `run_remd.py:204` resubmits only while
`blk < MAX_BLOCKS`. Block counts are **79HIS 7, 79HIS_S115T 7, 79ALA 10, 79ALA_S115T 7** -- every one
already past 5. So the four running jobs print `no resubmit`, the four `afterany` dependents each run
one further block, and then glpG stops.

**Nothing expires, and an earlier note here claiming otherwise was wrong.** It said `MAX_BLOCKS` had to
be raised *before the dependents finish* or the campaign would need restarting. That is not how this
works: `49010274`-`49010277` were submitted with `--export=ALL,...,REMD_MAX_BLOCKS=5` already baked into
their job records, so editing `submit_remd.sh` now cannot reach them, and they will not resubmit either
way. Continuing the campaign means editing `MAX_BLOCKS` and running `submit_remd.sh <V>` again, which
picks up the existing replicas and can be done at any point after the dependents finish. There is no
deadline to beat.

This supersedes the notes below claiming `MAX_BLOCKS=4` pinned with `=9` dependents; the file was edited
after those were written, and 79ALA at block 10 is past 9 in any case.

Timing as of 2026-09-13 01:00 CDT: `REMD_WALL_SEC=129600` (36 h) per block, so `49006599` has ~13 h left
on its block and its dependent would then run into Monday. There is enough data in hand for the Monday
meeting without extending anything.

**rockfish is retired as a host.** Access came from a former PI's allocation and was stopped on
2026-09-11; `30791125`-`30791128` were cancelled after 14 h 33 m. That leaves ~14.5 h of post-fix
data on an identical protocol at
`/scratch4/rherna21/ywang268/upside2-md-rf/popepopg_REMD/<variant>/`, currently **not retrieved**.
Copying it is read-only and does not consume the allocation, so it is still recoverable; decide
before that scratch is reaped. Do not submit anything else there.

**The fix is confirmed in production** (79HIS, all 28 rungs, completed chunk):
protein temperature **1.016-1.031** across the ladder against 1.10 (cold) to 1.84 (hot) before,
lipids 1.007, zero bad frames in that chunk, TM4 healthy everywhere (cold rungs TM4a 0.998,
TM4b 0.944; ladder means TM4a 0.930, TM4b 0.894) against a pre-fix TM4_full of 0.589-0.727 where
**0 of 60** replica measurements passed the > 0.8 criterion. Blow-ups: **434x fewer** bad frames
(findings 3.10d). All of this at a ladder ceiling that was *raised* 0.82 -> 0.90.

**Four NODE_FAIL requeues, and they are not this job's fault.** Slurm requeued all four variants at
identical seconds (21:55:47, 01:05:47, 01:35:47, ~11:35) on *different* nodes, which is an
infrastructure event, not hardware and not the change: the pre-fix job ran 12 h 52 m with none and
the 32 benchmark jobs were untouched. Each requeue re-runs `run_remd.py`, which increments
`block_count` and **truncates the log** (same `%j` output file), so only the current block's log
survives.

**Two traps this exposed, both now guarded:**
* `block_count` is read once at job start, so every requeue burns a block. It reached 4 of 4, which
  would have ended the chain silently. It is now pinned at **4** so the running jobs (which carry
  `MAX_BLOCKS=4` from their submit-time `--export`) can never resubmit, and continuation is handled
  by the `afterany` dependents with `MAX_BLOCKS=9`.
* **Do not "fix" that by resetting `block_count` to 0.** A later requeue would then read 0, get
  `blk = 1 < 4`, and resubmit itself *alongside* the dependent -- two jobs writing the same replica
  files. That mistake was made and reverted on 2026-09-11.
* Archive from the cluster that will run the jobs. Archiving `block_count` from midway3 while
  relaunching on midway2 left the compute nodes reading a stale copy on the shared `/project`, so all
  four variants started at block 2. The replica `.up` files were correctly absent and rebuilt from the
  seeds, so only the counter was affected.

**The two campaigns use DIFFERENT glycine maps, verified by hash 2026-09-12.** Testing each map against
both mirrors: NP `prod_ff3` is symmetric under the **correct** mirror `i -> (-i) mod n` (max err 0.0000,
and 3.77 under a plain reversal), while glpG post-fix AND glpG pre_ff3 are symmetric only under the
**off-by-one** mirror (0.0000 under reversal, 3.35 under the correct one). ALA/SER/HIS controls sit at
~11 in both, so the test is sound. So NP carries the glycine treatment ff_3.0 was actually trained
against and glpG does not. **Do not compare the two campaigns on anything glycine-sensitive**, and do not
"fix" glpG to match without re-reading the note below.

Force field otherwise confirmed identical and genuinely ff3.0 in both: NP's
`nonlinear_coupling_environment/coeff` is byte-equal to `ff_3.0/environment.h5:energies` (`2c5619ee9b12`;
ff_2.1 is `ce82e8b9`), and NP's `rotamer/pair_interaction/interaction_param` is byte-equal to the glpG
post-fix config (`8636da4601c8`). The `ff_2.1` paths still in `build_np_ff3.py` are the documented split
-- bead geometry, rotamer counts, `bb_env.dat` and the rama reference state stay at ff_2.1 because
ff_3.0 retrained only `sidechain.h5`'s interaction tables and one array in `environment.h5`. Nothing
there is stale, and the NP run does NOT need restarting on force-field grounds.

**Still open, deliberately not changed:** all 23 GLY Ramachandran maps carry an off-by-one mirror.
Correcting it makes TM4 **worse** (TM4a 0.786 -> 0.531), because only the buggy map favours the
right-handed helix; see findings 3.10c. TM4's residual fraying at the hot rung is a force-field
property, not a bug.

### Campaign 2: ff3.0 re-benchmark of Peng et al. JCTC 2022 (launched 2026-09-09 ~23:55 CDT)

**Scoring job 49010900 COMPLETED 0:0 after 2 h 54 m (2026-09-13).** 21 arms scored: 15 of 16 from
native (alpha3D still simulating) and 6 of 16 de novo. Table at
`ff3_benchmark/scoring/score_arms.json`, log `bm_score.49010900.out`, copy of the text table in the
deck at `0914/figs/score_arms.txt`.

| | ff3.0 | FF2 published |
|---|---|---|
| native, mean TM | **0.583** | 0.55 |
| native, mean Ca-RMSD | **3.78 A** | 4.0 A |
| de novo, mean TM | 0.417 (7 of 16) | 0.42 |
| de novo, mean Ca-RMSD | 6.18 A (7 of 16) | 6.1 A |

These are like-for-like: Fig. S5 reports mean TM and mean Ca-RMSD, which is what the script computes.
The lowest-Ca-RMSD comparison against Fig. S4 is the per-protein one and mean lowest goes 2.35 -> 1.49 A
over the 15 natives. **ff3.0 is lower on 9 of 15** (ubiquitin, NuG2, NTL9, hyp, top7, proteinL, cspA,
proteinG, lambda) and the only regression is gpW, 1.4 -> 2.0. Sorted by FF2 difficulty the pattern is
clean: level on the six FF2 already folded to ~1 A, better on every harder one.

**The de novo aggregate moved when one arm landed, which is the proof that a partial mean is not a
result.** At 6 of 16 it read TM 0.450 / 6.01 A and looked like a win; NTL9 de novo then finished at
TM 0.22 and pulled it to **0.417 / 6.18 A**, level with or marginally behind FF2's 0.42 / 6.1 A. One
arm shifted the mean by 0.033 and flipped its sign relative to the baseline. Quote the native
aggregate (15 of 16); give de novo as "indistinguishable so far" and expect it to keep moving.

**Rescored 2026-09-13 22:xx as job 49012360** (`score_arms_dist.py`, COMPLETED 0:0 in 2 h 44 m, 22/22
arms, zero errors). Reproduced the previous run's values exactly where they overlap, e.g. proteinG
native 0.792 / 2.37 / 0.68. It additionally writes per-frame cold-rung TM and Ca-RMSD to
`scoring/dist/<arm>.npz`, which is what the paper-format figures are built from.

The per-protein residue mapping is calibrated by requiring the native arm's frame 0 to score 0 RMSD
against the reference; `hyp` and `cspA` need non-zero offsets and that check is what caught the mapping
bug that made hyp read TM 0.267.

**Live progress 2026-09-12 ~08:55 CDT.** **7/32 COMPLETE** (BBA_native, BBA_denovo, BBL_native,
gpW_native, NTL9_native, WWdomain_native, WWdomain_denovo); 25 RUNNING at 12-18.5 h of a 36 h wall.
Still running: hyp, NuG2, alpha3d, BBL_denovo, gpW_denovo, homeodomain, NTL9_denovo, proteinG, lambda,
proteinL, ubiquitin, top7, proteinB, cspA. De novo arms set the campaign length; ubiquitin_denovo is
the longest and is still days out.

The long de novo arms set the campaign length. `ubiquitin_denovo` has done 561 k of 8.03 M time
units in 11 h, a rate of ~51 k/h, so it needs roughly **6-7 more days** and about 4 more wall blocks.
Expect the native arms to finish within 1-2 days and the de novo arms about a week out. Percentages
above are against each job's *remaining* target, not the full Table S2 duration, so they reset at
each resubmission; measure absolute progress from `run.log`'s first column plus what is archived.

32 jobs on midway2 broadwl, **relaunched a second time 2026-09-10 ~09:15 CDT** once ff_3.0 became
the arm-R tables; the 03:50 launch had been built against arm M and was cancelled 38 min in (parked
in `armM_discarded/`). All 32 RUNNING, and all 32 configs verified against the deployed libraries:
`nonlinear_coupling_environment` present, `sigmoid_coupling_environment` absent, `coeff` equal to
`ff_3.0/environment.h5:energies` and the rotamer pair table equal to `ff_3.0/sidechain.h5`.
16 proteins x {native, de novo}, 14-replica REMD, per-protein temperature ladders and durations
verbatim from Table S2, dt 0.009, frames every 100 time units, exchange every 10.

The first launch (`49002902`-`49002933`, cancelled after ~3 h) built the wrong environment coupling;
its `runs/` and `logs/` are parked in `sigmoid_discarded/` (4.0 G) and can be deleted.

* **Working dir** `/beagle3/trsosnic/yinhan/ff3_benchmark`, runs under `runs/<prot>_<kind>_ff_3.0/`,
  logs in `logs/`. Structures in `pdb/` came from Xiangda Peng's
  `/project2/trsosnic/pengxd/Data/test-set-15`, each verified against the residue count in Fig S4.
* **They resubmit themselves.** Several arms need far more than the 36 h wall (ubiquitin de novo is
  8.03 M time units, days of it), so each job runs under a `--time-limit`, archives its output to
  `output_previous_N`, reseeds from the last frame and re-queues for the remainder. Progress is
  measured from frames actually written. A run is finished when `runs/<tag>/COMPLETE` appears.
* **Baseline to beat** (Fig S5, 16-protein means): FF2 from-native TM 0.55 / Ca-RMSD 4.0 A, de novo
  TM 0.42 / 6.1 A. Per-protein FF2 values and the terminal-residue exclusions are in `plan.md`.
  **Peng's own `ff_2.1` is byte-identical to ours** (`799dc192...`), so the published numbers are a
  valid baseline rather than an approximate one.
* **Weakest FF2 cases, where ff3.0 has room:** lambda (5.4 A), hyp (11.1 A centroid), protein G
  (4.2/7.9), cspA (4.2/5.5), NTL9 (8.2 centroid).
* **Not yet built: TM-score.** Ca-RMSD can use `example/01.GettingStarted/calc_rmsd.py`; TM-score
  needs the standard binary or an implementation before the comparison is complete.
* **The environment coupling was wrong on the first launch, and is fixed.** ff_3.0 retrained
  exactly one array in `environment.h5`, `energies`, and only `nonlinear_coupling_environment`
  reads it; `scale`, `center` and `sharpness` are byte-identical to ff_2.1. `ConDiv.py` trained
  with `environment_potential_type = 0` (nonlinear), but `upside_config`'s CLI default is `1`
  (sigmoid) and `bench_run.py` did not override it, so the first 32 configs never read the
  retrained table and measured ff_3.0's sidechain against ff_2.1's environment. `bench_run.py` now
  passes `environment_potential_type = 0` for ff_3.0.
* **Pre-create `runs/<tag>/input` before submitting.** Eleven of the first 32 resubmissions died
  instantly in `os.makedirs(inp_dir, exist_ok=True)` with `FileNotFoundError`: 32 jobs racing to
  create nested directories on beagle3 at once. The directories are now pre-created, which removes
  the race; the eleven were resubmitted with `--exclude=midway2-0003`.
* **No ff_2.1 control was run.** The comparison is against the paper's published numbers, so it
  inherits any difference between their analysis and ours. An internal ff_2.1 arm would remove that
  and costs another 32 jobs.

### Campaign 3: NP 1AO6 + MPA-AuNP rebuild on ff3.0 (midway2, started 2026-09-10 ~02:51 CDT)

| JobID | what | log |
|---|---|---|
| 49003158 | `np_build_ff3`, builds the six cardinal-orientation configs | `NP-1AO6/build_ff3_49003158.out` |

Rebuilt from scratch rather than patched: the old replicas carry 98 accumulated output groups
(246 GB) and their coordinates are old-force-field unfolded.

* **Working dir** `/project/trsosnic/yinhan/NP-1AO6`, building into `prod_ff3/`. Driver
  `build_np_ff3.py`, checker `verify_np_ff3.py`, run under `build_np_ff3.sbatch` against the shared
  beagle3 Upside.
* **Built with `--environment-type nonlinear`**, which is the coupling ff_3.0 was actually trained
  on, so unlike the 32 benchmark jobs this one does read the retrained `energies` table.
  `--coverage on` matches both glpG arms and the replicas these configs replace, and
  `--rama-library rama3.dat` carries the GLY-symmetric maps.
* **Production is not submitted yet.** `submit_np.sh` still points at `prod/` and the old
  `upside2-md-mdw2`; repoint both at `prod_ff3/` and the beagle3 Upside before launching, and keep
  `NP_DT=0.001`.
* **Four bugs were fixed to get the build to run**, all real and all previously masked because the
  driver had never been executed since its 09-08 edit: `str | Path` annotations that need Python
  3.10 while the shared venv is 3.9 (`py/martini_itp_reader.py`,
  `py/martini_prepare_system_lib.py`, both now `from __future__ import annotations`);
  `np.exceptions.VisibleDeprecationWarning`, which needs NumPy >= 1.25 against the cluster's 1.23
  (`py/martini_prepare_system_lib.py`, now version-tolerant); `SYSTEM[...]` where the dict is
  `SPEC` (`build_np_ff3.py`, four references); and an exact-equality rotamer-geometry check that
  fired on 5e-13 A of float64 round-trip noise, now compared against a 1e-9 A tolerance that also
  prints the measured deviation every run.

### Campaign 4: NP production on ff3.0 -- rollback-on-tear driver, block 2 running (2026-09-12)

| JobID | what | status |
|---|---|---|
| 49003839 | block 1 of 8 | COMPLETED 2026-09-11 11:06, orientation 1 tore at frame 2948; chain rolled back |
| **49009692** | block 2 of 8 | RUNNING 9h05m elapsed, **all 6 orientations: 0 stretched bonds**, ~28h wall remaining |

### NP re-measured 2026-09-13 on 69 blocks: the campaign still cannot test Carlson, and is now worse

Measured with a correct minimum image on all six orientations (`NP-1AO6/np_state.py`, completed blocks
only). Native 1AO6 Rg(CA) is ~27 A and every run starts at 26.0-26.7, so the metric is anchored.

| run | Rg first | Rg last | Rg max | bound % | COM-NP last | contacting residues |
|---|---|---|---|---|---|---|
| 0 | 26.3 | 74.3 | 75.3 | 10.7% | 91.4 A | 503-522 |
| 1 | 26.5 | **130.4** | 130.4 | 97.1% | 143.7 A | 35-293 |
| 2 | 26.0 | 38.6 | 39.4 | 97.7% | 42.8 A | 38-577 |
| **3** | 26.5 | **108.7** | 108.7 | **0.1%** | **185.6 A** | never (1 frame) |
| 4 | 26.7 | 114.4 | 114.4 | 100.0% | 164.3 A | 498-579 |
| 5 | 26.1 | 41.4 | 42.4 | 97.6% | 32.8 A | 12-438 |

**Run 3 is the control that settles it.** It never binds -- 0.1% of frames in contact, centre of mass
185 A from the particle -- and it still expands from 26.5 to **108.7 A**. At T = 0.8647 this force field
denatures albumin in bulk, with no nanoparticle involved. Every "NP-induced deformation" statement is
therefore unsupported, and Carlson's "opens and exposes its centre" cannot be tested here because the
opening happens without the particle.

**This is a regression against block 2.** On 2026-09-12 run 3 was the compact unbound case at Rg 33 A and
run 0 was 64 A. Both have since blown past 100 A and 74 A. The free protein is not drifting to a swollen
equilibrium, it is still unfolding.

**Run 1 shows what "bound" now means.** 97.1% of frames in contact, Rg 130 A, COM 143 A away: an
essentially fully extended chain anchored to the particle at one end. That is not adsorption of a folded
protein and its contact set is not a footprint of albumin.

**Verdict unchanged and now unimprovable by more sampling:** blocks 4-8 add frames to denatured chains.
The blocker is protein stability, and the one number never revalidated is `run_np_prod.py:28`
`NP_TEMP = 0.8647`, inherited from Peng's FF2 calibration (T = 0.86 assigned 298 K) and never redone for
ff3.0. The cheap decisive test is still one orientation at T = 0.70-0.75 with run 3's unbound geometry as
the control. **Continuing the current production is spending compute on an uninterpretable result.**

### NP block-2 analysis, 2026-09-12: two findings that outrank the Carlson comparison

**0. It is NOT the ion concentration; measured 2026-09-12.** The hypothesis was that the NP box carried
too much salt and denatured albumin. Measured from the config rather than the build script: box 300^3 A,
**K+ 2423 = 0.149 M**, Cl- 2205 = 0.136 M, against 2439 pairs for exactly 0.15 M. The Cl- deficit is
counterion balance for the anionic MPA coating and albumin. For contrast the glpG system, which holds its
fold, runs **Na+ 201 = 0.186 M** nominal, i.e. *more* concentrated. Local condensation does not explain it
either: K+ enrichment within 8 A of protein beads, first frame -> last, is run 0 `0.33 -> 1.74` (Rg 69.6),
run 3 `0.27 -> 1.17` (Rg 53.1), run 2 `0.18 -> 3.00` (Rg 34.2). **The correlation runs backwards** -- the
most ion-enriched run is the most compact -- and 1.2-3.0x counterion enrichment around a charged protein
at 0.15 M is ordinary polyelectrolyte behaviour.

**The temperature is the number to question instead.** `run_np_prod.py:28` hardcodes
`NP_TEMP = 0.8647`, measured in the trajectory as the single fixed T = 0.8647 T_up = **303 K**. That value
comes from Peng's calibration, where **T = 0.86 was assigned 298 K for FF2**, and that calibration has
**not been redone for ff3.0**. If the scale moved at all, the NP is simply running hot; glpG's ladder by
comparison starts at 0.70 and its TM4 still frays at the 0.90 end. Cheap test: one orientation at
T = 0.70-0.75, and see whether albumin stays compact. Also note run 3 kept expanding through the day,
Rg 33 -> 53 A, so by 2026-09-12 evening all six trajectories are expanding, bound and unbound alike.

**1. The protein unfolds WITHOUT the nanoparticle, so deformation cannot be attributed to it.**
Independently verified: run 0 reaches **Rg(CA) 64.1 A while 165 A from the NP centre**, having never
bound (0.3% of frames in contact). Run 3 is also unbound (206 A) but stays compact at 33 A. So at
T = 0.8647 (303 K) this force field unfolds isolated albumin in bulk, which removes the premise of the
campaign: the bound runs' deformation is not demonstrably NP-induced. **Treat this as a protein-stability
question in ff3.0 before spending more NP compute.**

**2. Only 4 of 6 orientations ever adsorb, and they are frozen footprints, not an ensemble.**
Bound-frame fractions: run 0 **0.3%**, run 1 94.3%, run 2 95.4%, run 3 **0.2%**, run 4 99.9%, run 5 95.1%.
Pairwise Spearman between per-run lysine profiles falls to **-0.05** (run 4 vs run 5); **no lysine exceeds
0.2 in all four** adsorbed runs. Each run binds one contiguous subdomain and never reorients in 2.7 M
steps: run 1 -> IA, run 2 -> IIIB+IB, run 4 -> IIIB, run 5 -> IIA/IIB. **A contact fraction of 0.000
therefore means "that face was never presented", not "never contacts"** - which is exactly K190's
situation. Blocks 3-8 add frames, not faces. Testing a specific lysine needs a run seeded with that face
toward the particle, or a method that lets the protein reorient.

**Not the old PBC failure.** Intra-protein MARTINI pairs are skipped at runtime
(`skip_pair_if_intra_protein`) and all BB-NP pairs within 12 A use image (0,0,0), so these Rg values are
real chain extension rather than wrap-around into a second image. The 230 A artifact is gone.

**`np_footprint.py` has two real bugs; any earlier footprint conclusion used the broken version.**
* It sets `CB_PLACEMENT = [0, 0.94375626, 1.2068012]`, which is CB relative to **CA**, and adds it to the
  **N/CA/C centroid**. The config's own centroid-frame value is `[-0.0198, 1.5117, 1.2068]`, so it is off
  by **0.568 A on every residue** - the same sidechain-anchor error already recorded as fixed elsewhere,
  which never propagated here. Rebuilding CB from the config's own `affine_alignment/ref_geom` agrees with
  the Upside engine to 2e-4 A.
* It applies **no minimum image** in the 300 A cubic box while coordinates are stored unwrapped (ions have
  drifted to +/-3000 A).

**Carlson et al. 2025 claims, measured over 21,941 healthy frames (CB within 8.0 A of GOLD/MPA,
minimum-imaged; 67 of 22,008 frames carrying a torn peptide bond excluded):**

| claim | measured | verdict |
|---|---|---|
| K190 most protected | **0.0000**, rank 43/58, closest approach ever **14.0 A** against 3.9-4.4 A for genuine contacts | contradicted by this data; **INSUFFICIENT DATA** for the ensemble claim, since only 4 faces docked |
| K525 protected | 0.0000 at 8 A, 0.0011 at 10 A, closest 9.1 A | **INSUFFICIENT** - a near-miss, not a null |
| K541 protected | **0.12** bound-normalised, rank 18/58, 0.46 in run 4 | **SUPPORTS**, weakly |
| K12, K73 protected | K73 **0.15** (rank 16), K12 **0.09** (rank 22) | **SUPPORTS**, weakly, but in different single orientations |
| "opens and exposes its center" | opens: Rg 26.7 -> 32/39/50/74 A, kappa2 0.08 -> 0.19-0.77. Centre: contact-weighted native r_COM **31.6 A vs 25.3 A** unweighted, Pearson **+0.598**; f>0.30 set mean r_COM 39.7 A, never-contacted set 20.6 A | opening **SUPPORTS** but is not NP-attributable (see finding 1); "exposes its center" **CONTRADICTS** - the footprint is measurably peripheral |

Top sites are K573 0.44, K560 0.32, K574 0.28, K545 0.25, K93 0.25. 18 lysines sit at exactly 0.0000.
**Do not quote the pooled six-run ranking as an ensemble result**; it is the union of four frozen faces,
and dropping the two unusable runs (0 and 1) collapses it to the IIIB patch plus the IIA patch.

**Unfolding is the objective and it is working.** Over block 1 run.1's protein Rg rose smoothly
52.1 -> 55.0 A with covalent bonds entirely normal (mean C-N 1.327-1.347 A, max 1.68-1.83 A, zero
bonds over 2.0 A) for the first **2947 of 3011** frames, and the other five orientations are clean
throughout. Do not mistake a rising Rg for damage here.

**The tear is a separate event.** Unfolding does not stretch a covalent bond. Onset frame **2948**
(t = 82.54); all **63** frames to the end stayed torn with no recovery, the potential stepping from a
healthy -11611 to -10306 E_up; the final frame carried 7 C-N over 2.0 A (worst 3.22 A) plus a CA-C at
3.20 A, in one contiguous stretch, residues 439-450. At `k = 48` a 3.22 A peptide bond is ~100 kT, so
it is not thermal. It is the same mass-1-backbone-versus-MARTINI-core tear as glpG, rare here because
dt = 0.001 puts the one-step-kick radius at ~2.33 A against glpG's 3.23 A. NP logged
`avg_kinetic_energy/1.5kT` = 0.998-1.005 on all six and runs at the 0.8647 design temperature, so it
never had the glpG temperature defect and neither of the 09-10 fixes applies to it (no
`/input/brownian`).

**`run_np_prod.py` now rolls back instead of ending the chain** (backup:
`run_np_prod.py.bak_pre_rollback_20260911`). `reseed` restarts each orientation from its last healthy
output frame; torn frames stay in the file as history and are simply never used as a restart point,
so detection is unchanged and nothing is masked. A system that cannot produce a healthy frame for
`NP_MAX_STRIKES` (3) consecutive chunks still ends the chain. Verified by dry-run against the real
data before deploying: run.1 restarts at frame 2947 discarding exactly the 63 torn frames, and the
five clean orientations restart at the last frame exactly as before.

**The trap that dry-run caught, and it matters.** The restart criterion must be STRICTER than the
detection criterion. `CN_COUNT = 5` separates healthy (never more than 2 stretched bonds) from a
catastrophic tear (279-431), so it is the right *trigger*, but selecting a restart frame with it
picked frame 3004 -- which still carried 4 stretched bonds -- and would have propagated the very
damage the gate exists to stop. Restart eligibility is therefore `RESTART_CN_MAX = 0`: a restart
point must have no stretched bond at all. **Do not conflate the two thresholds.** dt stays pinned at
0.001 in `np_prod.sbatch` and the 2.0 A gate is untouched.

### Campaign 6: HDX dG on the grown trajectory (submitted 2026-09-13, midway3)

| JobID | variant | work dir |
|---|---|---|
| 59041160 | `glpG-RKRK-79HIS` | `popepopg_REMD_mdw2/<V>/hdx_10k/` |
| 59041161 | `glpG-RKRK-79HIS_S115T` | same pattern |
| 59041162 | `glpG-RKRK-79ALA` | same pattern |
| 59041163 | `glpG-RKRK-79ALA_S115T` | same pattern |

Same launcher and same settings as Campaign 5 (`HDX_N=28`, `HDX_DISCARD=500`, `HDX_LIVE=1`), so the
**only** variable changed is how much trajectory is available: usable frames per replica went
**6,121 -> ~9,600-10,000** (+60%) as the chain ran on. Writes to a new `hdx_10k/` so the Campaign 5
results that are currently in the 09/14 deck are not overwritten.

**This pipeline runs on midway3 only, and that is not a preference.** The venv it needs,
`.venv_el8_py311_bak`, has `bin/python3 -> /software/python-3.11.9-el8-x86_64/bin/python3.11`, and
`/software` is per-cluster: that target exists on midway3 (el8) and not on midway2 (el7), so on midway2
`python3` falls through to `/usr/bin/python3`, which has no h5py, and the job dies in step 1. The shared
`env_shared.sh` venv is not a substitute because it carries no pymbar or matplotlib. This is the
documented exception to the midway2 default.

**All four COMPLETED 0:0 on 2026-09-13, 7-8 min each**, verified rather than taken from the exit code:
28/28 replicas every variant, 9,465 / 9,782 / 9,956 / 9,723 frames per replica, the off-temperature seed
block correctly dropped in all of them, and zero FAIL/Traceback lines. Four 79HIS replicas each dropped
one rolled-back frame, which is the existing detector working.

**Result: the conclusions are unchanged, so the profiles are converged.** Off-scale amides at T = 0.85
go 69->64, 66->62, 69->65 and 81->68 of 203; the ESS-based resolution limit deepens uniformly from ~5.9
to ~6.2 kcal/mol, which is why the counts fall; resolved medians wander a few tenths with no trend; TM4
censoring is unchanged at 0-1 of 21. **The earlier "79ALA_S115T is mildly tighter" reading was sampling
noise and is retracted** -- at the larger frame count it is back with the others. Detail in findings 5.3d.

**Open question, deliberately not changed for this run.** The HDX analysis topology is built with
**ff_2.1** tables (`1.config.py` passes `parameters/ff_2.1/{sidechain,environment,hbond}.h5`) while the
trajectory it analyses is ff_3.0. That topology is used for projection and for the geometric protection
state, and the MBAR energies come from the joined trajectory rather than from it, so it should not affect
the result -- but it has not been verified, and it was left alone here so that Campaign 6 differs from
Campaign 5 in frame count alone. Changing two things at once would make the comparison uninterpretable.

### Campaign 5: post-fix HDX dG on midway3 (submitted 2026-09-12 ~07:40 CDT)

| JobID | variant | work dir |
|---|---|---|
| 58910978 | `glpG-RKRK-79HIS` | `popepopg_REMD_mdw2/<V>/hdx_postfix/` |
| 58910979 | `glpG-RKRK-79HIS_S115T` | same pattern |
| 58910980 | `glpG-RKRK-79ALA` | same pattern |
| 58910981 | `glpG-RKRK-79ALA_S115T` | same pattern |

**All four COMPLETED 0:0 on 2026-09-12 07:43**, ~5 min each. Input verified purely post-fix:
6921 frames/replica - 300 (seed block `output_previous_0`, skipped) - 500 (discard) = **6121 frames**,
which is exactly what the join produced, 0 non-finite potentials, 28/28 replicas, x4 variants.

Reads the **live** midway2 replicas over the shared `/project` (23-24 REMD blocks) with `HDX_LIVE=1`, so
this is a snapshot, not the final dataset; re-run after the chain reaches block 9.
Logs `popepopg_REMD/logs/hdxpf.<jobid>.out`. Launcher `popepopg_REMD/hdx_cluster.sbatch`.

**Result: global protection rose, but TM4 did NOT become non-exchanging.** Non-exchanging residues at
T=0.85 went **43/203 (pre-fix) -> 66-81/203 (post-fix)**, so the temperature/GLY fixes did stabilise the
protein. TM4 (131-152) gained new off-scale spikes near 120-123 and 143-147 and its resolved median rose
only 2.2-2.8 -> 2.6-3.0 kcal/mol, staying **0-2 of 21 censored** against TM1's **16-18 of 21**.

**That is not TM4 being unfolded.** Measured on the same post-fix trajectory, TM4 is 0.86 helical by
phi/psi and in the protected state in **97.9%** of frames (`PS_combined` 0.979, `PS_protein` 0.974)
against TM1's 0.991/0.948. Censoring requires `mean_pf >= 1 - 1/ESS` (~0.9998); TM4 sits at 0.97-0.99
while 16/21 TM1 amides sit at ~1.000, so TM4 resolves to a finite ~3 kcal/mol instead of going off scale.

**Resolved 2026-09-12: TM4 IS well protected, and the regional median was the misleading statistic.**
Per residue at T=0.85 (pooled 28 replicas x 6120 frames, `dg_limit` 5.90):

| TM4 residues | dG (kcal/mol) | 1 - protection | reading |
|---|---|---|---|
| 129-134 (N-cap/interface) | 0.15-1.23 | 8e-2 to 4e-1 | genuinely weak, drags the median down |
| 135-139 | 2.58-3.14 | 4e-3 to 1.8e-2 | intermediate |
| **140-146 (helix core)** | **3.46-4.85** | **5e-4 to 2e-3** | **strongly protected** |
| 147 | 1000 (sentinel) | 0.00e+00 | the one censored TM4 amide |
| 149-151 (C-cap) | 1.15-1.44 | 3e-2 to 5e-2 | interfacial |

So "TM4 median 3.0" averaged a protected core (~4.3) with weak caps (~1.0). **TM4's core is stable and
protected; it is not unfolded and not anomalous.**

**Why it does not go off scale, quantitatively.** Censoring needs `1 - pf < 1/ESS`, which at this ESS is
**4.7e-05**. TM4's best core amide (144) sits at 4.8e-04 -- real protection, but 10x above the resolution
limit, so it resolves at 4.73 rather than censoring. TM1 by contrast is `1.000000` exactly (**zero**
exchange events in 171,360 frames) for residues 36-42, which is why it saturates the sentinel. The
contrast is 99.95%-protected versus literally-never-exchanging, not folded versus unfolded.

**Mechanism, decomposed per frame 2026-09-12 and NOT what was first written here.** Protection combines as
`protection_t = 1 - (1 - pp_t) * acc_t`, so an amide counts as exchanged only when protein protection
fails **and** it is water-accessible in the same frame. Counting both over 171,382 samples:

| | H-bond/burial flicker `pp_fail` | `acc` | exchanged |
|---|---|---|---|
| TM1 30-48 | **0.0369** | **0.0021** | 3.13e-05 -> spikes |
| TM4 135-151 | **0.0322** | **0.4230** | 1.11e-02 -> does not spike |

**TM4's backbone protection is slightly BETTER than TM1's.** The two differ by ~200x in `acc` alone, and
that is the whole explanation. Per residue it is unambiguous: **res 36 flickers 7.0%**, the worst of any
amide, but `acc=0.0000` so it logs **0** exchange events and spikes; **res 140 flickers only 1.2%**, six
times less, but `acc=0.195` so it logs **371** events and resolves at ~4 kcal/mol.

**Therefore the off-scale plateau is largely a lipid-burial map, not a protection map.** Whenever
`acc_t = 0`, protection is identically 1 regardless of H-bond state, so helicity barely enters. TM1 does
**not** "never exchange" -- its H-bonds break 3.7% of the time; the lipid simply hides every break. Do not
read a spike as evidence of secondary structure.

**Correction to an earlier entry: within the TM4 core it IS a helical face.** Across the full 131-152
window the 21-residue end effect dominates the variance, which is what was measured first. In the core,
`acc` for 140-147 runs 0.20, 0.39, 0.42, 0.016, 0.039, 0.654, 0.066, 0.001 -- 141/145 (i, i+4) both
exposed, 143/147 (i, i+4) both buried, i.e. ~3.6-4 periodicity. One TM4 face is lipid-facing and the other
points into the protein interior and the catalytic cavity. Eight residues, so suggestive not settled.
**The one structural oddity left:** TM4's lipid-shielded stretch is only ~5 residues (143-147) against
TM1's 19 (30-48). A helix crossing a ~30 A hydrophobic core should shield ~20, so in this model TM4 sits
shallow/tilted. That is plausible for glpG, whose catalytic cavity is water-filled near the midplane
(S115 is the catalytic serine in this numbering, hence the S115T variants), but it has **not** been
checked against the crystal structure. That comparison is the remaining open item -- a structural
question about the model, no longer a suspected analysis bug.

**Figures regenerated 2026-09-12 with T=0.90 added, and NOTHING else changed.** The only edit to
`plot_ref_style.py` is the `--temperatures` default, now `0.75,0.80,0.85,0.90`; T=0.90 leaves the fewest
amides off scale (73/203 against 105 at T=0.75). The `.npz` was not recomputed. Local repo copy and both
cluster copies kept byte-identical.

**A "censored amides as bounds" rendering was tried the same day and REVERTED at the user's direction.**
It replaced the off-scale excursions with hollow carets sitting on each temperature's resolution limit,
broke the profile line across them, and retightened the axis to `(-4, 8.6)`. Both changes were wrong and
the reasons are worth keeping:
* **Breaking the line fragments the profile.** The continuous excursion is what makes each temperature
  read as one curve; gapping it at every censored amide turns the figure into disconnected islands.
* **Collapsing every censored amide onto `dg_limit` asserts they are all equal to ~6 kcal/mol**, when the
  actual statement is "unmeasurably large". It also caps the visible range at the limit, which is a worse
  distortion than running the excursions off the top.
The off-scale-excursion rendering in the docstring is a deliberate, documented choice ("reads as one
continuous excursion rather than a capped plateau ... that is how these profiles are conventionally
read"). **Do not replace it.** Keep `Y_LIMITS = (-20, 30)`.

**Still true, and the real caveat on the numbers:** the markers above
`0.001987 * temp_scale * ln(ESS)` (5.16-5.96 kcal/mol here) rest on a reweighted `1-p_f` below one
effective frame. They order amides correctly but are lower bounds, not quotable dG. That is the
docstring's own warning and it applies to every value above ~6 on the figure, including the 18-19
kcal/mol excursion feet.

### Rockfish pooling evaluated and REJECTED (2026-09-12)

It is genuinely poolable, checked rather than assumed: all **112 force-field datasets byte-identical** to
midway2 (so rockfish production did install the reference `c67351ca...` FF, as the plan required), same
friction 0.252482 / 4529 Brownian atoms / box 99.768^2 x 180 / 28-rung ladder with run.0 = T=0.70, and
healthy `avg_kinetic_energy/1.5kT` 1.001-1.034. Rockfish is still SSH-reachable despite the allocation
being stopped; data at `/scratch4/rherna21/ywang268/upside2-md-rf/popepopg_REMD/<variant>/`, 28 replicas,
17-18 GB each.

**But the payoff is 4 residues.** Rockfish holds only 2715 frames/replica against midway2's 6921, so
+31% samples raises ESS 21,418 -> 28,057, the limit 5.90 -> 6.06, and resolves 127 -> **131 of 203**. Not
worth ~72 GB off a dying allocation.
**The plot was never sampling-limited.** Of 84 censored amides at T=0.85, **59 have exactly zero exchange
events** in 171,382 samples, so they are censored at any achievable ESS; only 25 are data-limited, and
resolving those needs ESS > 171k, i.e. ~8x more simulation. The ceiling grows as `ln(ESS)`, so pooling
cannot fix a sparse plot -- rendering and rung choice did.
Rockfish remains useful for exactly one thing: an **independent-seed convergence cross-check** under an
identical Hamiltonian. Its `/scratch4` may be reaped, so retrieve it only if that check is wanted.

**Writes to `hdx_postfix/`, NOT `hdx/`.** The Sep-4 `hdx/` results are the pre-fix dG baseline and are
deliberately preserved; do not point a rerun at `hdx/`.

**`HDX_N=28` must be passed explicitly.** `hdx_cluster.sbatch` defaults to `N=48` (the retired midway3
48-replica ladder). The midway2 ladder is 28, and a wrong N makes step 2's `replicas done: n/N` check
fail the job.

**Two rounds failed in 2 s each (`58910584`-`58910592`, then `58910849`-`58910852`). Two independent
causes, both now fixed; `env.sh` backup is `env.sh.bak_pre_venvfix_20260912`.**

*Cause 1 -- the shared `.venv` was swapped out from under this pipeline.* `env.sh` activated
`$UPSIDE_HOME/.venv`, which was **rebuilt 2026-09-10 from the portable `pyrt` 3.9 interpreter**. That
venv carries numpy/scipy/h5py/tables/prody/Bio but **NOT pymbar and NOT matplotlib**, so the HDX
pipeline cannot run in it at all, and its interpreter additionally needs `pyrt/lib` on
`LD_LIBRARY_PATH` (which `env.sh` never set) or it dies on `libpython3.9.so.1.0`. The right environment
is `.venv_el8_py311_bak` (el8 python 3.11.9: pymbar 4.0.3, matplotlib 3.11.1, scipy 1.17.1, prody 2.6.1,
Bio 1.87) -- the venv `env.sh`'s own `module load python/3.11.9` was written for, and the one that
produced the Sep-4 results.
**Do not "re-unify" this on the pyrt 3.9 venv.** That venv exists for binary portability across both
clusters; HDX is pure Python, runs only on midway3, and needs the fuller 3.11 stack.

*Cause 2 -- and this is the trap: **that backup venv's `activate` is poisoned by the rename.*** It was
created as `.venv` and later renamed, so `bin/activate:38` still hardcodes
`VIRTUAL_ENV="/beagle3/.../upside2-md/.venv"` and line 42 does `PATH="$VIRTUAL_ENV/bin:$PATH"`.
**Sourcing `.venv_el8_py311_bak/bin/activate` therefore puts the CURRENT pyrt 3.9 `.venv` on PATH** --
it silently activates the very environment you were trying to avoid, which is why round 2 failed with
the identical `libpython3.9.so.1.0` error as round 1. `env.sh` now sets `VIRTUAL_ENV`/`PATH` by hand and
never sources that activate. Verified on a compute node before resubmitting: python 3.11.9 from the
right venv, all eight HDX imports OK, step-1 script runnable.
**A renamed venv's `activate` is not relocatable.** Check `grep VIRTUAL_ENV= <venv>/bin/activate`
against the directory it actually lives in before trusting any `.bak` venv.

**`midway3-0014` is excluded, and it is genuinely broken, not flaky.** A node probe shows it carries
only **2** beagle3 mount entries against 4 on healthy nodes, and
`/beagle3/trsosnic/yinhan/upside2-md` is simply **absent** there, so any job of ours that lands on it
cannot see the deployment. It failed this way on 2026-09-01 (job 57033313) and again 2026-09-12;
0019/0026/0050/0053/0061 are all fine. Treat it like `midway2-0003`.
**midway3-0201 is NOT bad** -- it sees beagle3 correctly; its failures were entirely Cause 1/2 and it
needs no exclusion.

**A known real bug still lurks in step 1 for some variants.** Job 57033314 (2026-09-01) died in
`upside_config.py:815 _input_phi` with `IndexError: index 642 out of bounds for axis 0 with size 630`:
the HDX topology's `input/pos` is stride-3 (N/CA/C, 210x3 = 630) while `_input_phi` indexes stride-4.
The Sep-4 run got past it, so it is variant/-path-dependent rather than universal; if a variant fails
there again, fix the stride in `_input_phi`, do not skip the variant.

### midway3, otherwise idle apart from the HDX jobs above

`squeue -u yinhanw` on midway3 is **empty** (checked 2026-09-09 08:14). Its last activity of any
kind was 2026-09-04 (`hdx_glpG-*` COMPLETED 08:56, then two `upside-gly-sym` attempts that FAILED
and were cancelled). The NP campaign is stopped, and the glpG detergent campaign that used to live
here is retired along with the DDM model (2026-09-09), so nothing is expected to run on midway3 for
glpG at all.

**Trap: midway3's accounting database holds stale `RUNNING` rows.** `sacct --clusters=all` reports
`53233848 remd_glpG-RKRK-79HIS RUNNING 29-02:54:01` and `53233852 ... RUNNING 29-02:52:52`. These
jobs are **not running**, they are August records that never received a final state, most likely
because the controller lost them. `squeue` is the authority for what is live; treat any multi-week
`Elapsed` in `sacct` as a zombie row, not a long job.

### glpG production is IDLE — all four variants finished 2026-09-05

| JobID | variant | elapsed |
|---|---|---|
| 48974448 | `mdw2_glpG-RKRK-79HIS` | 21:04:56 |
| 48974449 | `mdw2_glpG-RKRK-79HIS_S115T` | 1-04:34:07 |
| 48974450 | `mdw2_glpG-RKRK-79ALA` | 1-05:23:26 |
| 48974451 | `mdw2_glpG-RKRK-79ALA_S115T` | 1-05:03:45 |

Two consequences for `decide_and_launch.sbatch`, which was written assuming production was live:
its `scancel mdw2_glpG*` now has nothing to cancel, and its replica-deletion step would delete
**completed pre-ff3 baseline trajectories** (132 GB under `popepopg_REMD_mdw2/`) rather than a
half-finished run. Re-read that script before letting it fire.

The baseline *measurement* is already saved as text and is not at risk:
`popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt` (2026-09-05 13:09), replica run.0 = T=0.70 rung. It
records the failure ff_3.0 has to beat — 79HIS TM1 `mean_helix=0.906` / **TM4 `0.645`**, GLY49
`median_phi=+75.9` and GLY133 `+67.7` (alphaL) — against the pass criterion of TM `mean_helix > 0.8`
and negative (alphaR) GLY phi.

### The two training runs share history, so the arm test is a weaker check than it looks

rockfish resumed from midway2's half-trained checkpoint rather than training from scratch, so the
two runs are identical up to the branch at **step ~274** and have diverged only over the 222 steps
since. Their separation grows at exactly the same sqrt(t) rate as one run's own drift
(rel/sqrt(steps-since-branch) settles to 0.0176, against 0.0169-0.0200 within a run), which is the
signature of two independent random walks from a common point. Full detail and the table are in
`findings.md`.

Consequence for reading the arm-test verdict: the two arms are **more alike than two independent
trainings would be** (extrapolating the same rate to a full 500 independent steps gives ~40%
separation, against the 26% actually measured). If arms M and R agree on TM1/TM4, that shows this
pair agrees; it does NOT show that a retraining is reproducible. A real reproducibility test needs a
run branched at step 0.

### The two training runs produce DIFFERENT force fields, and only one may be used

Rockfish resumed from the Mac's step-274 checkpoint and midway2 from its own step-338; both then
advanced independently, so their step-500 parameters differ by the accumulated stochastic difference
of ~150 contrastive-divergence steps. They are two samples of the same training procedure, not two
copies of one result.

**Production on the two clusters must therefore install the SAME force field, or the two glpG
campaigns are not poolable** — they would be different Hamiltonians, and combining their
trajectories or their HDX ΔG estimates would be meaningless.

Decision recorded: **midway2's force field is the reference.** Two reasons — its lineage is the
unbroken one (its Adam state was never exported and re-imported), and `check_continue.sbatch`
installs it there automatically whatever else happens. So after the chain fires, copy
`midway2:parameters/ff_3.0_trained/sidechain.h5` to rockfish and patch the rockfish seeds from
**that** file, not from rockfish's own training output. Verify with the SHA-256 digest check that
already caught nothing wrong on the checkpoint transfer.

Rockfish's own step-500 force field is still worth keeping as an independent replicate of the
training procedure — it is the only evidence available about run-to-run spread in the trained
parameters — but it must not be mixed into production.

### NP rebuild at a 500 A box — staged, waiting on the force field

`scratchpad/ff3_retraining/np/build_and_launch_np.sh <sidechain.h5> <environment.h5> [--launch]`
does the whole thing: build six runs on rockfish, verify, transfer, launch on midway2.

* **Build on rockfish, run on midway2.** `build_np_ff3.py` needs ONE interpreter with both h5py and
  pytables. Rockfish's venv has both; midway2 has them split (venv -> pytables, module python ->
  h5py), and pip-installing into midway2's venv would disturb an environment a queued training job
  depends on.
* **Pass midway2's force field**, not rockfish's — see the divergence note above.
* **Box 500 A** via the new `--box-len`. A molecule sees its own image once its extent exceeds
  L minus the 12 A cutoff, so the old 200 A box was honest only to ~188 A while albumin reached
  Rg 230.9 A: that structural readout was PBC-contaminated. 500 A is honest to ~488 A, at ~11292 ion
  pairs against 723 at 200 A, so the systems are ~3.5x larger and correspondingly slower per step.
* **The replicas go in `prod/`.** `np_prod.sbatch` hardcodes `NP_RUN_DIR="$BASE/prod"` and globs
  `np.run.*.up`; `submit_np.sh` logs to `$B/prod`. Using `prod/` leaves both untouched, and it holds
  no replicas since the deletion. The script **refuses to proceed if any `np.run.*.up` is still
  there**, so two box sizes cannot end up in one run directory, and it writes `prod/BOX_LEN_A`
  because the box is not visible in the directory name.
* Every transfer is md5-verified per replica.
* NP resources, unchanged: 1 node, **6 cpus-per-task**, 36 h wall, self-resubmitting to MAX_BLOCKS,
  and `NP_DT=0.001` — a hard limit for this system, since unfolding drives backbone bonds to
  large-amplitude oscillation where accuracy, not MARTINI LJ stability, sets the step.

### The local reference set is on a DIFFERENT Mac, not this one

`planned_job.md` lists `scratchpad/ff3_retraining/` as the local artifact set (four-arm experiment
record, `verify_ff.py`, `converge.py`, `drift.py`, `handoff_to_mdw2.sh`, copies of every deployed
chain script, `ConDiv_original.py`). **None of it exists in this working copy**, and neither does
`training/gly-sym/` or `scratchpad/rf.sh`: this checkout last committed 2026-09-03 and received the
Sep 7-8 work as two `pull: Fast-forward`s, which carry the `.md` files but not a gitignored
directory. That work, and the Mac trainer, ran on another machine.

Consequence: from this machine the **cluster copies are the only readable copies** of the chain
scripts, and the four-arm experiment record is unreadable here. Do not plan a step that reads
`scratchpad/ff3_retraining/*` without first checking it is present.


### Disk: from midway2, `df` on the subdirectory is the ONLY number you get for `/project`

**Rewritten 2026-09-23 after re-measuring.** The earlier version of this section told you to read
the `Midway3 GPFS mounted at /project` row out of `rcchelp quota`, and to use a `check_quota.py`
helper. Both instructions are dead: on midway2 `rcchelp quota` now emits 14 lines covering only
**home, scratch and project2**, with no `/project` and no `/beagle3` row, and no `check_quota.py`
exists anywhere in the repo or on either cluster.

What actually works, all three verified 2026-09-23:

```
df -h /project/trsosnic     -> 3.9T total  3.5T used   445G free  89%   (the FILESET, correct)
df -ih /project/trsosnic    -> 1.1M inodes 216K used   898K free  20%
df -h /project              -> 6.3P total                                (the whole device, useless)
rcchelp quota               -> home, scratch, project2 group only
mmlsquota                   -> "File system project is not known", the GPFS client is not here
```

`df` on the **subdirectory** reports the fileset because GPFS `--filesetdf` is on; `df` on the
**mount point** reports the 6.3 PB device. That distinction is the whole trap.

**Current headroom, and the trend, which is the part that matters:**

| fileset | 2026-09-09 | 2026-09-23 | note |
|---|---|---|---|
| `/project/trsosnic` | 1514 G free | **445 G free** | ~1.1 T consumed in two weeks |
| `/beagle3/trsosnic` | - | 1.4 T free (4.2 T of 5.5 T) | the benchmark tree is only 42 G |
| `/project2/trsosnic` (group) | - | 1.45 T of a 1.49 T soft quota, **97%** | nothing in this campaign writes there |

`/project` is the one to watch. The four glpG variant directories in `popepopg_REMD_mdw2` hold
~150 GB each and `popepopg_REMD` holds another 638 G; those two trees plus `NP-1AO6` at 491 G are
1.75 T of the 3.5 T used.

### sbatch propagates the submitter's environment — this broke the whole chain once

**2026-09-06, found the hard way.** `check_continue.sbatch` requests `#SBATCH --mem=8G`, which puts
`SLURM_MEM_PER_NODE` in its environment. `sbatch` hands that environment to the job it submits, and
`srun_mdw2.sh` requests `--mem-per-cpu`, which sets `SLURM_MEM_PER_CPU`. With both present every
worker launch died instantly:

```
srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive.
```

All 12 workers failed, `run_minibatch` raised `All jobs failed`, and the job exited in under a
minute. `check_continue` resubmitted, the new job died the same way, three times — then the
forward-progress guard correctly aborted the chain at `stall 3/3` with training frozen at 236/600.

**This path had never executed before.** Every training job until then was submitted by hand from an
interactive shell, which has no `SLURM_MEM_*` set. The one job `check_continue` did submit was
cancelled before it ran. So the chain's only resubmit path was broken from the day it was written and
nothing revealed it. Had it stayed hidden, the unattended run would have aborted and produced nothing.

Fixed by unsetting the three mutually-exclusive variables after the `#SBATCH` block in every script
that submits another job — `check_continue.sbatch`, `run_arm_test.sbatch`, `decide_and_launch.sbatch`.
Unsetting them does not change the running job's own allocation; Slurm has already granted it.

Verified by execution, not inspection: `check_continue` was rerun with no dependency so it performed a
real resubmission, and the resulting job showed 0 `mutually exclusive` errors, 12/12 worker outputs,
96 replica `.h5` files (12 workers x 8 replicas), 0 `WORKER_FAIL`, and `MinMemoryCPU=2000M` as the
only memory variable.

**The general rule: a nested `sbatch` inherits `SLURM_*` from the job that calls it.** Before adding
any new submitting script, sanitize that environment or match the memory-request *type* of the child.
`remd.sbatch`, `np_prod.sbatch` and `armtest_remd.sbatch` set no `--mem` at all, so they were never
exposed — but they are equally reliant on the caller not leaking a conflicting pair.

### Slurm snapshots the batch script at submission — two bugs came from this

**A requeue restarts training from a STALE checkpoint and destroys newer ones.** Job 48977118 was
requeued by Slurm after a node failure (0010 → 0011). A requeue reruns the batch script *verbatim
with its original arguments*, so it resumed from the `epoch_02_minibatch_07` it had been submitted
with, and `main_loop`'s `rmtree` deleted checkpoints 08–21 on the way back up. Net progress was zero
from 15:32 to 17:37 and the wall clock resets on every requeue, so it would never have terminated on
its own. Fixed two ways in `srun_mdw2.sh`: `#SBATCH --no-requeue`, and the script now resolves the
newest checkpoint on disk at run time and ignores a staler argument (it logs when they differ).
Verify with `scontrol show job <id> | grep Requeue` — must be `Requeue=0`.

**An edit to a `.sbatch` does not reach jobs already queued.** `check_continue` 48977123 was submitted
at 12:47 and ran at 17:37, using its 12:47 snapshot — so it submitted 200 steps even though the file
on disk had said `STEPS_PER_JOB=180` for hours. 200 × 637 s = 35 h 23 min against a 36 h wall, and a
wall-limit kill is the trigger for the latent NaN path. **After editing any script in the chain,
cancel and resubmit the already-queued jobs that use it, or the edit silently does nothing.**

### midway2 — POPE/POPG REMD (correct seeds, all 6 GLY fixed)

**No glpG or NP job is running.** All of them ran to `COMPLETED` on 2026-09-05 and nothing resubmitted
them, because the only thing that will is the chain's `decide_and_launch.sbatch`:

| JobID | Name | Ended | Outcome |
|---|---|---|---|
| 48974448 | `mdw2_glpG-RKRK-79HIS` | 2026-09-05T16:24 | COMPLETED after 21:04:56. The other three variants finished alongside it. |
| 48974470 | `np_1AO6_prod` | 2026-09-05T15:59 | COMPLETED after 1-03:19:47. |

This is the intended resting state, not a fault: production is deliberately idle on the old force
field until the trained tables exist. The consequence is that `decide_and_launch.sbatch` will find
nothing for its `scancel mdw2_glpG*` to cancel, which is harmless — the resubmit step is what matters.
The pre-install TM baseline is already captured in `popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`
(79HIS TM4 mean_helix 0.645, 79HIS_S115T 0.474, 79ALA 0.571), and the four `seeds/*.up` files are
in place for the coverage injection.

**NP is out of the chain** and must be started by hand from `NP-1AO6/build_np_ff3.py`; it is being
rebuilt from scratch rather than patched, because its replicas carry 98 accumulated output groups
(246 GB) of old-FF unfolded coordinates.

**Health measured 2026-09-03 13:45** (not inferred from exit codes):
* Protein is live, not frozen: `potential[:,0]` std = 58 to 621 across recent groups (frozen signature is 0.000).
* glpG Rg = 17.9 to 22.0 Å across all 28 replicas of all 4 variants; hbonds 85 to 202. Physically sane.
* No `nan`/`inf`/error lines in any of the 5 logs.
* **Top of the ladder is straining.** Positive `protein_potential` excursions (vs ~-1400 typical) are
  confined to replicas 26 (T=0.89) and 27 (T=0.90): 68 to 156 such frames per variant, peaking at
  +2905 with hbonds down to 88. These are the same replicas that trip the peptide-bond ROLLBACK
  (`final-frame N of 209 peptide bonds > 2.0 A`). Consistent with the known dt-too-large / LJ-core
  instability at the hot end, not with a healthy ladder. Do not widen the rollback gate; the gate is
  reporting a real event.
* **NP Rg is far above native albumin (about 27 to 30 Å)**: run0 64.9, run1 95.3, run2 72.4, run3 118.4,
  run4 48.9, run5 44.5 Å, and all six grew over the block. Same unphysical expansion flagged before
  (previous peak 230.9 Å). NP footprint conclusions remain unsupported until this is explained.

**79HIS requeue loop, fixed 2026-09-03 14:00.** Job 48971711 was requeued 5 times, every time
`NODE_FAIL` on `midway2-0003`, which `scontrol show node` reports as `ALLOCATED+NOT_RESPONDING`
(a hung node that slurmctld kept re-allocating). It completed zero chunks across 4.5 h.

What was wrong and what was done:
* **The dead node kept being reselected.** `--exclude=midway2-0003` added to `submit_remd.sh`, so it
  propagates to every self-resubmission of all four variants (backup: `submit_remd.sh.bak_pre_exclude`).
* **Node failures were burning the block budget.** `block_count` is a plain file in the run dir that
  `run_remd.py` increments at *every process start*, not per completed chunk, so 6 dead starts had
  advanced 79HIS to 6/12 of `MAX_BLOCKS` with nothing to show. Reset to 0, putting it level with its
  siblings. **Any future NODE_FAIL requeue needs the same reset, or the chain silently ends early.**
* **Verified safe before resubmitting**: `midway2-0003` was `NOT_RESPONDING`, which can also mean a
  network partition with the process still writing. Confirmed no writer was active (log and h5 mtimes
  and log size unchanged over 75 s) and that all 28 replica files open with finite `input/pos`. The
  killed chunk's unflushed HDF5 buffers were discarded, so each file sits at its last consistent
  post-calibration state and `reseed()` continues cleanly from there.
* Left alone: `sc_env_transition_step_start` is 18000 for 79HIS against 39096-40618 for its siblings.
  That counter tracks cumulative simulated steps and correctly reflects that 79HIS has run far fewer,
  so it is consistent, not corrupt. It does undercount slightly, because a requeue does not set
  `REMD_LAST_CHUNK_STEPS` and `reseed()` then credits the default 2000 steps rather than the steps the
  killed chunk actually ran. Not hand-patched: editing that attr is a physics-schedule change.

**Seeds corrected 2026-09-03**: Restored `bak_broken_gly` seeds (properly equilibrated, `activation_stage=production`, `current_stage=production`, has `output` group), applied GLY symmetrization to all 6 helical residues (GLY49, GLY104, GLY128, GLY133, GLY156, GLY180), and expanded box Z from 123.7 Å to 180 Å (updated `martini_potential.z_len` in all 4 seeds). Old replica files deleted; block_count reset.

**VTF extraction note**: `output_previous_0` in each replica is the old seed's `output` group (single-T stage-7 data), moved there by `reseed()`. Skip it; REMD data starts at `output_previous_1`. Script: `/home/yinhanw/project/yinhan/extract_glpg_vtf.py` (skips block 0, wraps atoms into box, unwraps protein backbone chain).

28 replicas per variant (midway2 REMD uses 28-core nodes). Data: `/project/trsosnic/yinhan/popepopg_REMD_mdw2/<variant>/`.
Logs: `.../logs/remd.<V>.<jobid>.out`. Self-submitting via `bash submit_remd.sh $V`.

**CRITICAL — rigid-stage fix 2026-09-02**: Seeds were built with `current_stage = minimization` and
`preprod_protein_mode = rigid_body`. The `set_stage_label(seed, "production")` step was never applied,
so all prior trajectory data (output_previous_0 through output_previous_13 for 79HIS; similar blocks for
other variants) was produced with a completely frozen protein: `rama_map_potential std = 0.000` across all
frames. Fixed by `fix_stage.py`: patched `input/stage_parameters.current_stage` to `production` in all
116 files (4 seeds + 28 replicas × 4 variants), reset `block_count` to 0.
Seed backups: `seeds/<V>.up.bak_rigid_stage`.

**DATA EXCLUSION**: All trajectory data produced before 2026-09-02 (output_previous_* groups) in
`popepopg_REMD_mdw2/` has a rigid protein and must NOT be used for any conformational analysis or HDX.

**Verify fix after first chunk** (compare `rama_map_potential std` to NP control = 0.773):
```python
import h5py, numpy as np, os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
with h5py.File("/project/trsosnic/yinhan/popepopg_REMD_mdw2/glpG-RKRK-79HIS/glpG-RKRK-79HIS.run.0.up", "r") as h:
    for grp in sorted([k for k in h if k.startswith("output_previous_")], key=lambda k: int(k.split("_")[-1])):
        pot = h[grp]["potential"][:, 0]
        print(f"{grp}: rama_std={np.std(pot):.3f}")
```

Seeds: rebuilt from `hybrid_prep/` MARTINI structures using the corrected `upside_config.py`. Installed to
`/project/trsosnic/yinhan/popepopg_REMD_mdw2/seeds/<V>.up` (old seeds backed up as `.bak_alphaL_restart_2026-09-01`).
Protein is in alphaR (phi_std ≈ -60°, confirmed by BioPython). All 6 helical GLY sym_err=0.000000.
Cluster upside_config.py also patched at `/project/trsosnic/yinhan/upside2-md-mdw2/py/upside_config.py` (backup: `.bak_phi_criterion`).

**GLY map fix 2026-09-01 (complete)**: `inject_backbone_nodes` had a phi convention bug — it uses
a reversed-b0 `_input_phi` formula (phi_up = phi_std ± 180°) but compared against `(-150,-20)` which
captures alphaL GLY, not alphaR. All 6 helical GLY residues (GLY49, GLY104, GLY128, GLY133, GLY156,
GLY180) were never symmetrized during setup. Their maps were biased against alphaR by +0.65–1.65 E_up,
progressively destabilizing alphaR TM helices (manifests as H-bond loss at bilayer midplane — Bug 2).
Fix: `py/upside_config.py` criterion `(-150 <= phi <= -20)` → `(30 <= phi <= 160)` in both
`write_rama_map_pot` and `write_rama_map_pot2`.

**Before resubmitting again, ALWAYS verify seeds:**
```bash
source /software/modules/init/bash
module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1
export HDF5_USE_FILE_LOCKING=FALSE
cd /project/trsosnic/yinhan/popepopg_REMD_mdw2
python3 check_seeds_current.py   # must show sym_err < 0.001 for all 6 helical GLY
```

### midway3 (caslake) — context notes

**NP GLY fix 2026-09-01**: 1AO6 albumin had 4 broken helical GLY maps: GLY11, GLY81, GLY203, GLY244
(sym_err 3.0–5.3 E_up, alphaR bias +1.1 to +1.8 E_up). Prior job 56987370 cancelled.
Fixed all 6 seed files with the same periodic-mirror symmetrization as the glpG fix.
Backups at `prod/np.run.*.up.bak_broken_gly`. Check script: `~/project/NP-1AO6/check_np_gly.py`.

**HDX jobs (57041866-69, 57033313-16) — CANCELLED / DONE.** Data they analysed came from the broken
runs with biased GLY49/GLY133 maps. Do not use those HDX results. Rerun HDX after the new REMD
(57068431-34) accumulates sufficient data.

**sbatch fixed 2026-08-31.** RCC resolved the slurmctld spool issue. NP job self-submitted successfully as 56565841. REMD jobs will chain normally once they start. Tmux session `remd_popg` on midway3-login1 may still hold srun processes for REMD jobs from the workaround period — verify before relaunching if jobs start and immediately fail. REMD logs: `~/project/popepopg_REMD/logs/remd.srun.<V>.<jobid>.out`.

**REMD local bypass completed 2026-08-28/29.** Jobs 55758021-24 had been PENDING (Priority queue) since 2026-08-27. Cancelled all 4. Ran block 1 locally on Mac M1 Pro (finished 2026-08-29 ~02:12 EDT); extracted last-frame positions (~10 MB); uploaded and reconstructed 48 × 4 minimal `.up` files on cluster with `cluster_setup.py`. All 192 replica files verified 48/48 per variant with `block_count=1`.

**REMD clean restart 2026-08-27.** Cancelled pending jobs 55118057/55119260/55122401/55122658. Deleted wrong-map trajectory data (57-76 chunks generated with biased G136 Ramachandran map, ~228 GB freed). Reset `block_count` to 0 in all 4 variant dirs. Fresh REMD jobs 55758021-24 will reinitialize 48 replicas from `seeds/` at startup; seeds were patched correctly on 2026-08-21 with the original stride-4 code.


**CRITICAL — post-fix trajectories do NOT reproduce the paper's claims.** Block-2 measurement (2026-08-19, 18 246 frames, 3.2% adsorbed-and-compact):

| Paper claim (Carlson et al.) | Post-fix block 2 | Status |
|---|---|---|
| K190 is the most protected site (water + TES common site) | 0.000 contact — lowest of all 58 Lys | **CONTRADICTS paper** |
| K525 protected (water) | 0.000 contact | Contradicts paper |
| K541 protected (TES) | 0.033 contact | Contradicts paper |
| K12, K73 also protected | 0.000, 0.005 | Contradicts paper |
| "Opens and exposes its center" | Dominant patch: Lys313/Glu311/Asp314/Asp562/Lys560 (max 0.33) | Inconsistent |

Rg reaches 230.9 Å (run.3, block 3) — well past the paper's DPD-simulation range. No `np_footprint.npz` exists yet for the current block-3 trajectories (~231k frames/run). None of the paper's five target lysines are contacted; the footprint must be rerun after block 3 completes before drawing any updated conclusions.

**NaN cascade fix** is deployed in the running binary (`src/main.cpp:~517`, added `!isfinite(lboltz_diff)` to exchange rejection). **GLY Ramachandran stride bug** fixed 2026-08-27 (stride-4 restored, all 192 active replica maps patched). Both are in the current deployed binary; no action needed.

**Verify real dynamics after any stage-label or seed change:** internal CA RMSD (Kabsch-superposed) should grow across frames. A frozen protein gives 0.000–0.001 Å; a live one reaches 1.4 Å within 100 frames. Reading `current_stage` is not sufficient — measure the RMSD.

**`HDX_LIVE=1` is mandatory while the REMD jobs run.** `run_remd.py` renames `output` to the lowest free
`output_previous_<n>` at every chunk, so the live group can vanish between being listed and being read; the
flag makes `martini_remd_concat.py` skip it, at a cost of ~200 of 5043 frames.

### The rigid-protein bug — fixed 2026-08-17, keep the lesson (findings 116)

`/input/stage_parameters.current_stage` was **`production_handoff`** in the four seeds, and
`martini_hybrid.cpp:637-641` holds the protein as a rigid group whenever `preprod_protein_mode == rigid_body`
**and** the stage is not exactly `production`. Both were true, so 28–30 h × 4 ran with a frozen protein:
internal CA RMSD **0.000 Å** across whole chunks, H-bond count constant to ±0.01, identical at every rung.

**Fixed.** `set_stage_label(seed, "production")` on all four (backups at `seeds/<V>.up.bak_production_handoff`),
verified `production` / `production` / `rigid_body`, relaunched as 53441123–26.

**The lesson, which is why this stays here.** The 2026-08-15 note declaring these jobs healthy checked
`hybrid_interface_active_stage` (lines 644-648), which *does* accept `production_handoff`. That is a different
predicate on the same string with the opposite accept set. One string, two gates. Reading a gate is not
verification — measure the protein's internal CA RMSD instead.


**Confirmed loading on the new binary 2026-08-15 00:53**: all four reached `[remd] block 1/12` and
`calibration 2000 steps x 48` with `n_atom 4949`, no `expected 1 arguments but got 2`. The seed migration
held.

**These four were verified rather than resubmitted (2026-08-15).** The BB proxy now reads its backbone O
from `infer_H_O`, so `martini_hybrid_position` takes two arguments and every config written by the older
prep fails to load. Checked, in order: `env.sh` points `UPSIDE_HOME` at the rebuilt
`~/beagle3/yinhan/upside2-md`; the four seeds were upgraded in place with
`py/martini_upgrade_hybrid_args.py` and one was proved to load and run on the new binary; the variant
directories are empty, so `run_remd.py` materialises replicas from those upgraded seeds and there are no
stale copies and no `STOP` files; the seeds carry an `/output`, which the first-block `reseed()` requires;
their stage is `production_handoff`, which `martini_hybrid.cpp:646-647` treats as active, so the SC-env
interface is on — **this line is the error findings 116 corrected; that gate is not the one that decides
whether the protein moves.** The 36 h `REMD_WALL_SEC` equals the Slurm limit but the chunk loop guards with
`remaining() - MARGIN > last_wall*1.25` (35 min plus a chunk), so the chain still has room to reseed and
resubmit. **No resubmission needed.**

The NP job cannot be patched the same way — it holds its six replicas open for the whole block — so
`np_prod.sbatch` now runs the same upgrade in the gap before `upside` reopens them. **Confirmed working
2026-08-15 00:25**: the block rolled over to 53372760, all six reported `upgraded`, and the block then
loaded every config under the new binary. That line is a one-shot migration — delete it once a block
reports "already upgraded".

**Reading a `.up` that a running job owns.** `run_np_prod.py` reseeds by renaming `output` to the lowest
free `output_previous_<n>`, so on a live file that group can vanish between being listed and being read;
a first footprint attempt died on exactly that (`KeyError: object 'output' doesn't exist`). Rotated blocks
are never renamed again, so an analysis should read only `output_previous_*` and skip the live group —
that costs ~1.4% of frames and removes the race. Set `HDF5_USE_FILE_LOCKING=FALSE` to open at all, and
never let `upside_engine` touch a live file: it opens read-write.

Launcher `popepopg_REMD/{submit_remd.sh,remd.sbatch,run_remd.py,env.sh}`, seeds in `popepopg_REMD/seeds/`,
logs `popepopg_REMD/logs/remd.<jobid>.out`, data `popepopg_REMD/<variant>/`. Seeds are the stage-7.0
checkpoints from `popepopg_glpG/<variant>/checkpoints/`: 4709 atoms (1050 protein / 3393 lipid = 261×13 /
266 ions), box 99.869² × 123.697 Å, 7.48 M pairs, ions regenerated at 0.15 M with
`--membrane-thickness-angstrom 42.75`. This copy of `run_remd.py` uses a **per-chunk RNG seed**; the DDM
copy passes a constant one, which is what made its rollback re-run the identical failing chunk.

Seeds: stage-7.0 checkpoints from `popepopg_glpG/<variant>/checkpoints/`: 4709 atoms (1050 protein / 3393 lipid / 266 ions), box 99.869² × 123.697 Å. `run_remd.py` uses a **per-chunk RNG seed** (DDM copy used constant seed → rollback re-ran the identical failing chunk; now fixed).

## 2. THE TWO CAMPAIGNS ARE DIFFERENT SIMULATIONS

Conflating them caused several errors this session, including a threshold copied from NP that killed a
healthy 6 h glpG block. **Never transfer settings, thresholds, or analysis between them.**

| | **NP** (`np_1AO6_prod`) | **glpG** (`remd_glpG-*`) |
|---|---|---|
| method | regular MD, 6 independent trajectories, single T=0.8647, no exchange | **REMD**, 28 replicas, T ladder 0.70–0.90, configuration exchange |
| purpose | nanoparticle adsorption footprinting (K190 exposure) | **HDX** protection factors / ΔG |
| system | 1AO6 albumin 578 res + 5 nm MPA-AuNP, 8608 atoms, box 300 Å | glpG 210 res in a POPE/POPG bilayer. **Read the atom count and box from the seed**: two generations exist, 4949 atoms / 279 lipids / box 99.77² × 180 Å and an older 4709 / 261 / 99.869² × 123.697 Å |
| composition | PROTEIN 2890 + GOLD 887 + MPA 203 + ION 4628 (K+ 2423 / Cl- 2205, 0.15 M KCl) | PROTEIN 1050 + LIPID (13 beads each) + ions regenerated at 0.15 M; the counts follow the seed generation |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: ions, lipids and the 630 protein N/CA/C sites are on the single-stage g-JF **Brownian** path (4529 of 4949 atoms on the current seeds); the other 420 protein atoms are not in `/input/brownian` |
| timestep | **0.001**, freely settable at runtime | **0.009, HARD-LOCKED** by `/input/brownian/numerical_time_step`; `martini_brownian.cpp:100` throws on mismatch. Friction is tuned against it for lipid D=11.5 µm²/s — **do not change it** |
| detection | non-finite positions OR ≥5 stretched bonds | non-finite potential (whole chunk) OR ≥5 stretched bonds |

---

## 3. NP campaign — `np_1AO6_prod`

**Unfolding is the expected result, not a failure.** 1AO6 albumin spreads on the MPA-AuNP surface,
so a large Rg must **not** be reported as a blow-up. Judge health on non-finite frames, peptide C-N
bonds, and `avg_kinetic_energy/1.5kT`. (Contrast glpG, where Rg ~19 Å is the health signal.)

**The Rg printed in `np.<jobid>.out` is NOT periodic-image corrected, so do not read the campaign's
result off it** (measured 2026-09-18 on block 6). It is computed on the stored coordinates directly,
so on any face where the adsorbed chain straddles a box boundary it is inflated by roughly the number
of box lengths spanned. Logged vs minimum-image Rg, referencing every backbone atom to the MPA shell
centre: run0 123.1/128.8, run1 **184.0/76.4**, run2 96.9/102.9, run3 **332.5/118.7**, run4
**170.6/73.4**, run5 80.8/80.3. Three of six faces are inflated 2.3-2.8x, and run3's headline
"Rg 332 Å" in a 300 Å box is a boundary crossing, not more spreading. The minimum-image values are
themselves only reliable where the chain stays inside half a box; run0 (max 172 Å) and run5 (164 Å)
exceed that, so treat those two as unresolved rather than agreeing.

**The protein is adsorbed on all six faces**, which is what the footprint analysis needs and what Rg
failed to show. Backbone atoms within 8 Å of the MPA shell: run0 257, run1 1128, run2 913, run3 399,
run4 780, run5 1243, of 2312. Nothing has escaped the nanoparticle. Spreading is still real
(minimum-image Rg 73-129 Å against native albumin's ~27 Å), just smaller than the log implies.
**Use the contact count, not Rg, as the adsorption observable.** Atom layout in `/input/pos`:
backbone `0:2312` (578 res, stride 4), protein sidechain beads `2312:~3200`, Au core `~3200:3750`,
MPA carboxylate shell `3750:3950` (use this to locate the NP, it is unambiguous), ions `3950:8608`.
Box is 300 Å cubic (`martini_potential` attrs `x_len/y_len/z_len`).

**Dir** `~/project/NP-1AO6/` — `prod/` holds `np.run.{0..5}.up` + `np.<jobid>.out`, `block_count`. **Current configs are the envfull+300Å rebuild** (protein-protein terms injected, 4628 ions at 0.15 M KCl, 8608 atoms); block_count reset to 0.
**Driver** `run_np_prod.py` · **sbatch** `np_prod.sbatch` (sets `NP_DT=0.001`) · **submit** `submit_np.sh`
**Footprint analysis** `np_footprint.py` → `np_footprint.npz` (must be run; no current npz for block-3 data).
**Orientation map** (cardinal Euler faces; see `scratchpad/NP-footprinting/orientation_map.txt`):

```
np.run.0 = 0-0-0   (yaw=0.0,  pitch=0.0,  roll=0)
np.run.1 = 90-0-0  (yaw=90.0, pitch=0.0,  roll=0)
np.run.2 = 180-0-0 (yaw=180.0,pitch=0.0,  roll=0)
np.run.3 = 270-0-0 (yaw=270.0,pitch=0.0,  roll=0)
np.run.4 = 0-90-0  (yaw=0.0,  pitch=90.0, roll=0)
np.run.5 = 0-270-0 (yaw=0.0,  pitch=270.0,roll=0)
```

Self-resubmits up to `NP_MAX_BLOCKS=8`; chunks ~104 time units, ~1765 t.u. per 36 h block.
Local build source: `scratchpad/NP-footprinting/` (`build_all.py`, `np_hybrid.py`, six face dirs).

**Status check:**
```bash
f=$(ls -t ~/project/NP-1AO6/prod/np.*.out | head -1)   # newest block log
grep "^\[np\]" $f | tail -20          # driver trace + per-chunk bond counts
grep -ic nan $f                        # expect 0
```

---

## 4. glpG campaign — `remd_glpG-*` (HDX)

**Dir** `~/project/yinhan/popepopg_REMD_mdw2/` on **midway2** — one subdir per variant, each with 28
`*.run.N.up`, `remd.<jobid>.out`, `block_count`.
**Driver** `run_remd.py` · **sbatch** `remd.sbatch` · **submit** `submit_remd.sh <variant>`
**Variants:** `glpG-RKRK-79HIS`, `glpG-RKRK-79HIS_S115T`, `glpG-RKRK-79ALA`, `glpG-RKRK-79ALA_S115T`

Config: 28 replicas, T 0.70–0.90, `REMD_DT=0.009` (hard-locked), `REMD_MAX_BLOCKS=5` as of 2026-09-12
(read from `submit_remd.sh:33`; all four block counts already exceed it, see §1) for the ff3.0
campaign. The per-run flags live in `submit_remd.sh` on the cluster; that file is the authority.

**The midway3 detergent campaign is gone.** `glpG_DDM_micelle_REMD/` (and the older lamellar
`glpG_DDM_REMD/`) belonged to the DDM model, retired 2026-09-09. Its data is not deleted, but nothing
in this handbook points at it any more and no number from it is used as evidence.

### HDX analysis for the POPE/POPG campaign — RUNNING (submitted 2026-09-01)

Jobs 57033313–57033316, reading from midway2 block-1 data (`popepopg_REMD_mdw2/`, ~7000–7400 frames/replica, `HDX_LIVE=1`). Outputs land in `popepopg_REMD_mdw2/<variant>/hdx/`. Plot: `results/<V>_POPEPOPG_dG_vs_residue.png`.

To resubmit after more data:
```bash
B=/home/yinhanw/project/popepopg_REMD; MDW2=/project/trsosnic/yinhan/popepopg_REMD_mdw2
for V in glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T; do
  sbatch --job-name=hdx_${V} --output=$B/logs/hdx.%j.out --partition=caslake --account=pi-trsosnic \
    --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --time=04:00:00 \
    --export=ALL,PDB_ID=${V},HDX_SRC=${MDW2}/${V},HDX_LIVE=1,HDX_WORK=${MDW2}/${V}/hdx \
    $B/hdx_cluster.sbatch
done
```

It runs `py/martini_remd_concat.py` first — **required**, because the chained driver rotates `output` to
`output_previous_<n>` every chunk, so `/output` alone is just the last ~300 frames. The concat joins the
chunks in restart order, renumbers `time`, and drops whole any chunk with a non-finite potential (a
rolled-back chunk is not a sample). Verified against a 3-chunk file: boundary frames in place, all 16
datasets carried, strided output identical to the naive slice. Then the same path as the local run —
`example/00.AnalysisScripts` + `write_hybrid_energy.py` + `plot_ref_style.py`.

**The cluster copy of the analysis pipeline drifts from the repo.** On 2026-08-15 the cluster's
`calc_hdx_ht.py` and `4.calc_D_uptake.py` still lacked the reference subtraction from findings 91, so an
HDX run there would have silently returned uniform MBAR weights. Re-upload the analysis scripts before
trusting any cluster-side result; only the C++ build is kept current by `install.sh`.

**CRITICAL — exchange recirculates a destroyed configuration.** `run_remd.py` reseeds each replica from
`output/pos[-1]`. If that frame is destroyed the next chunk starts destroyed, then an exchange swaps a
healthy configuration in and the bad one moves to another slot. Consequences:

- A single bad configuration appears in many slot files over time (observed walking 3 → 13 → 23 → 33).
- **Per-replica "this file is clean" is meaningless.** Count destroyed configurations *per frame across
  all 48 slots*; the count was 1, then 2, and grew.
- It never self-clears. That is why a detection gate exists.

**Status check:**
```bash
for v in glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T; do
  d=~/project/yinhan/popepopg_REMD_mdw2/$v; f=$(ls -t $d/remd.*.out | head -1)
  echo "$v block=$(cat $d/block_count) nan=$(grep -ic nan $f) DESTROYED=$(grep -c DESTROYED $f)"
done
```

---

## 5. How to check TM helix health (TM1 and TM4)

This check is run frequently after any seed change or after the first trajectory chunk completes.
TM1 (GLY49 at C-cap) and TM4 (GLY133 at N-cap) are the two helices most sensitive to the GLY
Ramachandran bias bug and must be verified independently from the global health check.

### Step 1 — verify seeds before submitting (run once per seed generation)

```bash
cd /project/trsosnic/yinhan/popepopg_REMD_mdw2
module load python/3.11.9
export HDF5_USE_FILE_LOCKING=FALSE
python3 check_seeds_current.py
```

Expected output for a healthy seed (all 6 helical GLY, sym_err=0):
```
glpG-RKRK-79HIS:
  GLY49:  phi=-94.1  aR=0.965 aL=0.965 sym_err=0.000000  OK
  GLY104: phi=-88.6  aR=1.012 aL=1.012 sym_err=0.000000  OK
  GLY128: phi=-80.0  aR=1.581 aL=1.581 sym_err=0.000000  OK
  GLY133: phi=-141.6 aR=1.121 aL=1.121 sym_err=0.000000  OK
  GLY156: phi=-69.6  aR=1.300 aL=1.300 sym_err=0.000000  OK
  GLY180: phi=-68.5  aR=1.194 aL=1.194 sym_err=0.000000  OK
All seeds OK.
```

A broken seed shows `sym_err ~ 3–6` and `aL < aR`. **Do not submit if any seed shows BROKEN.**

**Current state (2026-09-01)**: Seeds rebuilt from `hybrid_prep/` MARTINI structures using corrected
upside_config.py. Protein is in alphaR (BioPython phi_std ≈ -60° for TM1 body). All 6 helical GLY
maps symmetric. helix_fraction using phi_std ∈ [-130,-20] should be nonzero from the start.

**WARNING — dihedral sign error in the VTF analysis script above**: the homemade `dihedral()` function
returns `-phi_std`. For alphaR residues (true phi_std ≈ -60°), it reports ≈ +60°. Use BioPython's
`calc_dihedral` for any phi verification instead of the function defined above.

### Step 2 — check helix health from a VTF trajectory (after first chunk)

Extract a VTF for the T=0.70 replica (slot 0) and analyse phi/psi:

```python
import numpy as np, re

def parse_vtf(vtf_path):
    """Return atoms list and positions array from a VTF trajectory."""
    atoms = []
    with open(vtf_path) as f:
        for line in f:
            m = re.match(r"atom\s+(\d+)\s+name\s+(\S+)\s+resid\s+(\d+).*chain\s+(\S+)", line)
            if m:
                atoms.append({"aid": int(m.group(1)), "name": m.group(2),
                               "resid": int(m.group(3)), "chain": m.group(4)})
            elif line.startswith("timestep"):
                break
    n_atoms = max(a["aid"] for a in atoms) + 1
    frames = []
    pos = np.zeros((n_atoms, 3)); count = 0; in_frame = False
    with open(vtf_path) as f:
        for line in f:
            if line.startswith("timestep"):
                if in_frame: frames.append(pos.copy())
                pos[:] = 0; count = 0; in_frame = True; continue
            if in_frame:
                if line.startswith(("pbc","bond","atom","#")) or not line.strip(): continue
                parts = line.split()
                if len(parts) >= 3:
                    try: pos[count] = [float(x) for x in parts[:3]]; count += 1
                    except ValueError: pass
    if in_frame and count > 0: frames.append(pos.copy())
    return atoms, np.array(frames)

def dihedral(a, b, c, d):
    b1=b-a; b2=c-b; b3=d-c
    n1=np.cross(b1,b2); n2=np.cross(b2,b3)
    l1=np.linalg.norm(n1); l2=np.linalg.norm(n2)
    if l1<1e-10 or l2<1e-10: return np.nan
    n1/=l1; n2/=l2
    m1=np.cross(n1,b2/np.linalg.norm(b2))
    return np.degrees(np.arctan2(np.dot(m1,n2),np.dot(n1,n2)))

def helix_fraction(atoms, frames, chain="A", res_range=(131, 152)):
    """Fraction of frames where res_range is helical (phi in [-130,-20] AND psi in [-90,15])."""
    n_by_r  = {a["resid"]: a["aid"] for a in atoms if a["chain"]==chain and a["name"]=="N"}
    ca_by_r = {a["resid"]: a["aid"] for a in atoms if a["chain"]==chain and a["name"]=="CA"}
    c_by_r  = {a["resid"]: a["aid"] for a in atoms if a["chain"]==chain and a["name"]=="C"}
    res_list = [r for r in range(res_range[0], res_range[1]+1)
                if r in n_by_r and r in ca_by_r and r in c_by_r]
    hel_frac = {}
    for r in res_list:
        phis = []; psis = []
        for pos in frames:
            if r-1 not in c_by_r: phis.append(np.nan); psis.append(np.nan); continue
            phi = dihedral(pos[c_by_r[r-1]], pos[n_by_r[r]], pos[ca_by_r[r]], pos[c_by_r[r]])
            psi = dihedral(pos[n_by_r[r]], pos[ca_by_r[r]], pos[c_by_r[r]],
                           pos[n_by_r[r+1]] if r+1 in n_by_r else pos[c_by_r[r]]) if r+1 in n_by_r else np.nan
            phis.append(phi); psis.append(psi)
        phis = np.array(phis); psis = np.array(psis)
        hel_frac[r] = float(np.mean(
            (-130<=phis) & (phis<=-20) & (-90<=psis) & (psis<=15)))
    return hel_frac

# Usage:
atoms, frames = parse_vtf("/path/to/glpG_79HIS_T0.70_slot0.vtf")
# TM4 health (residues 131-152, GLY133 at N-cap)
tm4 = helix_fraction(atoms, frames, res_range=(131, 152))
print("TM4 helix fraction per residue:", {r: f"{v:.2f}" for r, v in tm4.items()})
print("TM4 mean:", np.mean(list(tm4.values())))
# TM1 health (residues 29-49, GLY49 at C-cap)
tm1 = helix_fraction(atoms, frames, res_range=(29, 49))
print("TM1 mean:", np.mean(list(tm1.values())))
# GLY133 phi distribution
# expect: mostly in [-130, -20] for a stable TM4 N-cap
```

**Pass criteria (TM helix healthy):**
- TM4 mean helix fraction > 0.8 across residues 131–152
- TM1 mean helix fraction > 0.8 across residues 29–49
- GLY49 phi stays in [-130°, -20°] for >80% of frames
- GLY133 phi stays in [-150°, -20°] for >80% of frames

**Fail signal (biased maps still active):**
- TM helix fraction near 0 — the helix collapsed
- GLY49 or GLY133 phi drifting to +60° (alphaL) — the map is pushing it left-handed

---

## 5b. How to check health CORRECTLY (general bond/energy check)

**`isfinite` is not a health check.** At a real glpG failure the environment coordinates were
**±4.65e12 Å** — numerically finite, physically destroyed. And in a forced NP tear the protein reached
**431 broken bonds with the potential still finite at +3e5**, so no energy-based test fires at all.

Use the broken-bond **count** (healthy 0–2; torn 279–431 — a two-order-of-magnitude gap):

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

**NP**: a gate trip looks like `[np] DESTROYED ...` followed by `no resubmit`; job exits rc=0/COMPLETED
short of its wall limit. A COMPLETED job that did not resubmit means the gate fired — check the log.

**glpG (updated driver)**: a rollback looks like `[remd] ROLLBACK #N filename: reason` followed by
`rolled back M/48 replicas; continuing chain`. The chain does NOT terminate. The NaN output is rotated
to `output_previous_N` as normal history; the rolled-back replica restarts the next chunk from its
pre-chunk positions. A replica that repeatedly blows up gets rolled back repeatedly; it does not get
dropped from the ladder. Watch for high rollback counts on the same file — that indicates a replica with
a persistent physics problem that won't self-correct.

**Rollback mechanism**: before each chunk `run_remd.py` snapshots `/input/pos` (3.6 MB total for 48
replicas). On NaN detection it overwrites the last `output/pos` frame and `output/potential[-1]` with
the pre-chunk values so that `reseed()` on the next iteration picks up the clean state.

**These driver scripts are NOT in git.** They live on the cluster at
`~/project/yinhan/popepopg_REMD_mdw2/run_remd.py` (midway2) and `~/project/NP-1AO6/run_np_prod.py`.
No version history exists for them. Edit directly on the cluster.
A running job keeps the version it loaded at start; edits take effect at the **next block**.

---

## 7. If the chain terminates — manual rollback procedure

**glpG (old driver, or gate fired before rollback logic):** patch last output frame of each NaN file
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

- **NP dt hard limit: 0.001.** dt=0.005 caused backbone blow-ups during unfolding (large-amplitude spring instability at t>250, proven by A/B). Never raise above 0.001.
- **Do not transfer thresholds between NP and glpG.** A CN_MAX borrowed from NP false-positived on a healthy glpG chunk (2.52 Å vs healthy max 2.659 Å) and cost a 6 h block. The two jobs have different physics.
- **glpG NaN propagation via REMD exchange.** A single blow-up in one replica spreads to all 48 via exchange within ~60 steps (IEEE 754: `NaN < 0.f` = false). The NaN cascade fix (`!isfinite(lboltz_diff)`) and the per-chunk rollback driver address this.
- **A green exit code means nothing** for a self-submitting REMD job. Check the log for DESTROYED/ROLLBACK counts and verify physical observables.
- **Midway3 home quota**: 28.6 G of 30 G. Jobs can fail oddly if home fills.
- **Do not run scripts from `/tmp`** on the login node (another user's `/tmp/inspect.py` shadows stdlib).
- **`/tmp` on a login node is PER-NODE, and reconnects land on different ones.** A `nohup` job launched
  from midway2-login2 writes `/tmp/<log>` that is invisible from login1, so a later check reports the log
  missing and the job gone even if it is alive elsewhere. Put anything you intend to read again on
  `/beagle3` or `/project`, and prefer `sbatch` over `nohup` for work that must outlive an ssh session:
  a compute-node job survives disconnects and does not compete with other users on a login node.
- **`pgrep -f <name>` matches your own shell command.** `ssh host 'pgrep -f bm_final && echo running'`
  reports "running" because the remote `bash -c` command line contains the string, so a dead job looks
  alive indefinitely. Match the interpreter instead (`ps -u $USER -o cmd | grep python3`), or check for the
  output file, or use a Slurm job id and `squeue`/`sacct`. This produced two false "still running" reports
  on 2026-09-12.
- **A stray module in the working directory shadows the stdlib.** An `inspect.py` left in a scratch
  directory broke `import numpy` with a circular-import error, because numpy imports `inspect`. The same
  hazard as the `/tmp` note above but it applies to any working directory; name scratch probes something
  that is not a stdlib module.
- **glpG detergent campaign (closed 2026-08-13, model retired 2026-09-09).** Kept for one lesson only: it used a constant `--seed`, so a rollback re-ran the identical failing chunk deterministically, which is why the driver now takes a per-chunk seed. All four variants failed at block 2–3. Its HDX ΔG output is in `~/Downloads/glpG_DDM_micelle_HDX_dG/` and is no longer cited.
