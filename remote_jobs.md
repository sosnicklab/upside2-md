# Remote jobs on midway2/midway3 — status and handbook

**2026-10-01 11:20: the ff3.0 retrain (plan.md Phase 8) is QUEUED ON MIDWAY2 BROADWL (chain
49135913, estimated start 10-02 22:16) and TRAINS LOCALLY on the Mac Studio until that job starts.**
The 1.2M SU allocation reached the midway2 scheduler at ~11:15 on 10-01. Every half hour,
`sync_to_midway2.sh` copies each new local step into the waiting cluster run, staged, converted
and moved in atomically. When the chain is RUNNING, the script stops the Mac and writes
`HANDED_OVER`. Do not use midway3 or any partition but broadwl for this (user). BP validation
analysed (§0c), not yet reported to Tobin.

Written so a fresh session can pick up cold. Everything needed to connect, check health correctly,
and react to a failure is here. Job state below is live; superseded jobs are not listed, only
summarised in §8 where they carry a lesson. Earlier state still current: rockfish is gone as a host,
so glpG is midway2-only; key-based ssh is refused on both clusters, every connection costs a Duo push.

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

---

## 0. Connect first (needs a Duo push on the user's phone)

Key-based auth is NOT enabled; password + Duo is the only method. The ControlMaster socket expires
roughly hourly, so expect to redo this most sessions.

**midway2** (training, BP validation, and the glpG chains at the ff3.0 release):
```bash
ssh -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<command>'
```
If midway2 blocks or throttles this IP, tunnel through midway3 as a fallback:
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

* **Minimise connections rather than trying to remove Duo.** The training chain and the BP arms
  insure themselves (§1), so nothing depends on a live laptop connection and a status check every
  few hours is enough. **No cluster-side monitor runs now**: `monitor.sh` watched only `ff30` and was
  stopped with it, so `/project/trsosnic/yinhan/STATUS.md` is stale from 2026-09-27 23:54.
* **`scrontab` is disabled and `crontab` is denied** on this cluster, and `pi-trsosnic` has **no
  association with the `cron` partition** (`sbatch -p cron` is refused). A 1-core `broadwl` job is
  the way to schedule recurring work. A 7-day request hits `QOSMaxWallDurationPerJobLimit`; 36 h
  works, **provided the successor is queued at the start** (`--dependency=afterany:$SLURM_JOB_ID`):
  a monitor that resubmitted at the end of its loop overran the wall and died silently.
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
the next expect launch dies **after** it has already sent a Duo push. That is how a single
unguarded status call costs a push and locks you out for tens of minutes.

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
before one single retry. Do not relaunch the expect script to "see if it works" — each launch
spends a Duo push on the user's phone before it discovers the throttle. The way around it while it
lasts is the midway3 tunnel above, which needs a live `~/.ssh/cm-mdw3.sock` and therefore its own
Duo approval.

**ALWAYS put `-o BatchMode=yes` on routine ssh calls.** (Learned 2026-09-17.) A plain
`ssh -S ~/.ssh/cm-mdw2.sock host 'cmd'` does **not** fail when the socket is dead: it silently
falls back to a fresh connection, offers the key, then tries password auth twice non-interactively,
hits `Received disconnect ... Too many authentication failures`, and after a couple of those the
host starts closing connections immediately (`Connection closed by 128.135.112.69 port 22`). That
is the IP throttle described above, self-inflicted by a status check. With `BatchMode=yes` the same
call fails instantly and harmlessly with `Control socket connect: No such file or directory`, which
is the signal to run the expect script. Check the socket first, and never let a monitoring loop
retry the expect script more than once per tick.

Two zsh/tooling traps that cost time here:
* **zsh does not word-split unquoted variables.** `M2="ssh -S sock host"; $M2 'cmd'` runs silently
  and produces NOTHING, it is not a connection failure. Write the `ssh` call out in full.
* `timeout` does not exist on this Mac; do not wrap `expect` in it.

`~/project` on midway2 is a symlink to `/project/trsosnic` (not `/project/trsosnic/yinhan/`).
Data is at `~/project/yinhan/popepopg_REMD_mdw2/`.

**`/project` is the SAME filesystem on midway2 and midway3** (both mount `midway3_cap`), so
`/project/trsosnic/yinhan/upside2-md-mdw2/` is fully readable from midway3. Checkpoint progress,
chain logs and force-field directories can all be checked **without a midway2 login at all**,
useful when an IP block or a login outage is in the way. Only `squeue`/`sacct` need midway2
itself; the two clusters have separate Slurm controllers and separate accounting databases, so
`sacct --clusters=all` on midway3 does **not** see midway2 jobs.
Python env: `source /software/modules/init/bash && module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1 && export HDF5_USE_FILE_LOCKING=FALSE`

**midway3**:
```bash
ssh -S ~/.ssh/cm-mdw3.sock -O check yinhanw@midway3.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw3_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu '<command>'
```
`~/project` on midway3 is a symlink to `/project/trsosnic/yinhan/` (note: **yinhan**, not yinhanw).
Load the python env with `source ~/project/NP-1AO6/env.sh` before any h5py work.

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

Deployed by rsync from `/project/trsosnic/yinhan/upside2-md-mdw2` (`src py obj parameters cmake
example` and the install scripts), **excluding `training/`**, which is ConDiv campaign data that
does not belong in a code deployment. Benchmark output goes to `/beagle3` so it does not compete
with `/project`.

**Naming:** the release writes `parameters/ff_3.0` in both trees. The old `ff_3.0_trained` /
`ff_3.0_trained_rf` directories are still on the cluster but no live script references them
(checked 2026-09-28).

### Other Upside copies (inventory 2026-09-09)

Owned by yinhanw:

| path | size | branch | last touched | disposition |
|---|---|---|---|---|
| `/project/trsosnic/yinhan/upside2-md-mdw2` | 8.3 G | martini-dev | live | **ACTIVE**, `ff30_basin` training and the ff3.0 release |
| `/beagle3/trsosnic/yinhan/upside2-md` | 6.4 G | martini-dev | live | **the shared deployment** |
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

## 0c. Rotamer belief-propagation stopping-test bug: fix, deployment, validation (2026-09-26)

**The bug** (reported by John Jumper via Tobin). `NodeHolder::max_deviation` in `src/rotamer.cpp`
took `max(cur_belief - old_belief, dev)` with no absolute value. Beliefs are rescaled each iteration
so the dominant rotamer stays at 1, so an update that only *sharpens* (others fall) gives dev = 0 and
the solver stops as if converged. It has been there since the original rotamer commit: ff2.1 was
trained with it, and the running ff3.0 training uses it.

**The fix** is `fabsf(...)` on that line (279). Status by tree:

| tree | status |
|---|---|
| `master` (GitHub) | fixed in `288d1fec`, pushed by the user 2026-09-26 |
| `martini-dev` | fixed; master merged in by PR #48 (`11044069`) 2026-09-26 |
| `/project/trsosnic/yinhan/upside2-md-mdw2` (midway2) | fixed; installed by the ff3.0 gate 2026-09-30 02:24:59 |
| `/beagle3/trsosnic/yinhan/upside2-md` (midway2 + midway3) | fixed; installed by the ff3.0 gate 2026-09-30 02:24:59 |

**Static-frame result (local, 2026-09-26).** 8 training proteins (50-146 res), ff2.1, 91 frames each
(native + T 0.80/0.95/1.10), old vs fixed library vs a tol=1e-6 reference. The early stop fires in
~45% of frames but the error is small: |dE| <= 0.029 E_up, force error 7e-5 (median) / 3e-3 (max) of
the RMS force, ConDiv rotamer-gradient error ~1e-4 relative; the energy bias is one-signed (converged
is higher, mean 0.0009 E_up). Conclusion given to the user: ff2.1 and ff3.0 do not need retraining.

**Simulation validation for Tobin: all 16 arms COMPLETE 2026-09-29 03:50, analysed; not yet
reported to Tobin (the user's call).** Output `bp_validation/analysis_20260929.txt`, 25-32k
frames/replica after equilibration. Result:
* WWdomain is the only protein resolved against the seed spread: at T 0.764-0.879 the bug arms are
  less native, Q lower by 0.02-0.04 (z -1.8 to -12), RMSD higher by 0.1-0.3 A, E higher by 1-3; T_mid
  0.859 bug vs 0.861 fix (~0.7 K).
  Per seed (checked 2026-09-29), both bug seeds lie below both fix seeds in Q at every T from 0.700
  to 0.879; at 0.852 bug 0.511 / 0.505 against fix 0.547 / 0.550. At 0.700 the gap is one bug seed
  (0.842 against 0.89).
* NTL9 goes the same way at low T (Q 0.25 bug vs 0.30 fix, z -2 to -4; RMSD +1.0-1.3 A, z 1.1-1.8),
  but per seed that rests on one arm: fix s2 sits at RMSD 6.5-6.9 A while the other three sit at
  8.0-9.0 A. Its folded population is low under ff_2.1 here and T_mid is not resolved (0.821 vs
  0.815). Not resolved.
* proteinG and homeodomain: not resolved, |z| < 1.5 through the folded and transition range; T_mid
  0.861 vs 0.868 and 0.931 vs 0.929, inside the seed spread.
* Unfolded-T z of 2-5 sit on seed sd of 1e-4 in Q and 0.01-0.03 A in Rg; the differences are that
  small, and with 2-dof sigmas across 224 cells such z values are expected by chance.
* Speed: the fix costs 0-5% (proteinG 17.6 vs 18.6, WWdomain 29.3 vs 29.8, NTL9 28.0 vs 28.0
  time units/s); node-to-node variation is not controlled.

ff_2.1, native-start
14-replica REMD with the Peng benchmark protocol (Table S2 ladder and duration, dt 0.009, frame 100),
proteinG / homeodomain / WWdomain / NTL9 x {bug, fix} x seeds {1, 2}. Trajectories of the two
binaries decorrelate within a few hundred time units, so the test is statistical: bug-minus-fix
against the seed-to-seed spread (z = delta / sigma_seed). All 16 checked healthy at the first 100
frames (C-N 1.33-1.34 +- 0.13 A as in the ff_3.0 benchmark, KE/1.5kT 0.99-1.00).
* dir `/beagle3/trsosnic/yinhan/bp_validation`: `obj_bug/`, `obj_fix/` (deployment source built
  twice, differing only in the fix; `obj_bug` reproduces the deployed binary bitwise), `bp_run.py`
  (copy of `ff3_benchmark/bench_run.py`: binary, seed, ff_2.1 native), `bp.sbatch`, `submit_all.sh`
* logs `logs/<prot>_<build>_s<seed>_<jobid>.out`, data `runs/<prot>_<build>_s<seed>/`; each arm
  resubmits itself every ~29 h under the same job name, so track with `squeue -u yinhanw | grep bp_`
  and completion with `ls runs/*/COMPLETE`
* analyse on the midway2 login node after `source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh`:
  `python3 analyse_bp.py [protein ...]` gives per-T Q, CA-RMSD, Rg and E per build, bug-fix vs the
  seed spread, melting midpoints and speed. Two seeds give only a rough sigma: |z| ~ 1 is "not
  resolved"; only |z| well above 2 across the ladder is an effect
* local static-frame harness (gen.py, eval.py, cmp.py) was in the session scratchpad, not kept

**Deployment: fired by the ff3.0 gate 2026-09-30 02:24:59** (`ff30_basin/ff30b-gate_49131949.out`,
`[bpfix]` lines). `ff30_basin` trained on the old binary; its release validation (Peng arms, glpG
chains) runs on the fixed one. Both trees carry `BPFIX_INSTALLED` and the old binaries as
`obj/*.bak_pre_bpfix_20260930-022459`; `src/rotamer.cpp` has the `fabsf`. The gate script was
restored to the repo's version (md5 0a653ade). The BP test arms were unaffected (`obj_bug`/`obj_fix`).

---

## 1. Current jobs

Snapshot **2026-10-01 11:20 CDT, verified live against `squeue` on midway2.** Finished and cancelled
rows are deleted; only lessons worth reusing are kept below.

| JobID / where | what | state | next action |
|---|---|---|---|
| **49135913** | **ff3.0 retrain, AWH glycine library fixed** (plan.md Phase 8), chain link on broadwl, `$P/training/ff30_gly` | PD (Resources), est. start 10-02 22:16, 12 nodes; resumes from the newest step in its `run_output` when it starts | at start: confirm which step it resumed from; first cluster step without WORKER_FAIL; `check_step.py` |
| local Mac, `training/ff30_gly_local` | the same run, until the cluster job starts | driver PID in `driver.pid`, running since 09:16 under `caffeinate -i` with `CONDIV_LOCAL_WORKERS=2`; log `train_local.log`; ~94 min per step | `sync_to_midway2.sh` every half hour; stopped by it when 49135913 runs |

`$P` = `/project/trsosnic/yinhan/upside2-md-mdw2`. The cluster run_output now holds the local run
(marker `run_output/FROM_LOCAL`). The cluster's own initialised run_output is kept as
`run_output.superseded_20261001-111755`.

**Watch and hand-over.** Session cron job `3869da4e` runs at :07 and :37. It runs the sync,
`check_step.py` on every new step, and re-arms a live monitor. It is session-only: it dies with the
Claude session and expires on 10-08, so a cold session must re-create it.

**`bash training/ff30_gly_local/sync_to_midway2.sh`** (on the Mac) is idempotent:
* chain RUNNING: it stops the local driver and workers and writes `HANDED_OVER`;
* otherwise it stages each new complete local step in `$P/training/ff30_gly/sync_stage`, runs
  `training/move_run.py` there (paths, NumPy 2 -> 1.23.5; findings 10.0f), and renames the step
  into `run_output`, so a starting job never reads a half-copied or unconverted checkpoint;
* no chain job and midway2 accepts one: it submits the chain.

A local step that finishes after the cluster has started is discarded. Both run dirs hold
identical `init_param/` (ff2.1, md5-checked) and `upside_input/`, including `rama.dat` =
`parameters/common/rama31.dat` (md5 `fc479d45...`).

**ff30_gly, what must not be forgotten:**
* `upside_input/` is a hardlink copy of ff30_basin's, except `rama.dat`, which is a fresh copy
  of `parameters/common/rama31.dat` (md5 `fc479d45...`, built by `training/build_gly_library.py`,
  log `checks/build_gly_library_20261001.log`). Never edit a hardlinked file in place.
* `after_training.sbatch` gates with max 13 epochs. On convergence it releases `ff_3.0`, backing
  up the failed basin-offset release, and then tries to submit the Peng arms and glpG to broadwl.
  **From midway3 those submissions fail**, and `validate_ff.sh` stops at "not every benchmark arm
  was submitted": the release happens, validation does not. Submit validation by hand once midway2
  has an allocation, or adapt `validate_ff.sh` to midway3.
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
showed by 19:00 is kept here:

**TM4 under ff_3.0 at 19:00 (`checks/glpg_tm_windows_20260930_1900.log`,
`gly_tm4_series_20260930_1900.log`; 9-13 completed chunks of ~95 time units per variant).**
* **Handbook pass (T 0.70, TM4 134-151 > 0.8):** still met in all four, latest 0.93-1.00; 79ALA_S115T
  dipped to 0.82-0.86 over groups 5-9 and recovered to 0.93.
* **T 0.80, TM4 falling with time:** 79ALA 0.98 -> 0.59-0.82 over its last five groups (fails),
  79HIS 0.97 -> 0.82-0.86, 79ALA_S115T latest 0.84 (residue 141 at 0.45), 79HIS_S115T 0.84-0.96.
* **The mechanism is TM4's helical glycines at phi > 0.**
  * At T 0.80 all four variants flip: GLY136/143/149 up to 0.23-0.64 of frames in recent groups.
  * At T 0.70 79ALA_S115T flips GLY143, then GLY149 (0.26-0.74), and 79HIS has just started
    (GLY136 and GLY143 0.07 in its latest group).
  * The other two variants are clean at T 0.70.
* **TM1 (30-48, no glycine) also declines at T 0.70:** 79HIS 0.99 -> 0.89-0.90, 79HIS_S115T 1.00 ->
  0.91-0.93.
* **Baseline:** in the pre-ff3 campaign GLY143 never flipped, and T 0.70 gave TM4b 0.991 and TM1 1.000
  (findings 3.10c, slightly different windows).

**The glycine probe (49133133, `training/ff30_glyprobe`) and its insurance successor 49133135 were
cancelled by the user at 20:30 on 09-30, at 15 of 19 steps**, once the direction was settled (findings
1.15, `checks/glyprobe_final_partial_20260930.log`). The directory has no `after_training.sbatch`, so
nothing can resume or release from it; keep it for the record.

Logs and data: Peng `/beagle3/trsosnic/yinhan/ff3_benchmark/logs/<prot>_<mode>_ff_3.0_<jobid>.out`,
`runs/<prot>_<mode>_ff_3.0/`; glpG `/project/trsosnic/yinhan/popepopg_REMD_mdw2/logs/remd.<V>.<jobid>.out`,
`popepopg_REMD_mdw2/<V>/`. The superseded ff_3.0 benchmark is in `runs_superseded/ff_3.0_20260930-022510`,
the pre-release glpG seeds are `seeds/*.bak_pre_ff_3.0_20260930-022510`.

**A Slurm requeue used to count as a glpG block (fixed 2026-09-30 13:20).** `remd.sbatch` has no
`--no-requeue`, so after a NODE_FAIL Slurm restarts the job under the same id, and `run_remd.py`
incremented `block_count` at every start. Three chains lost 2-3 of their 5 blocks that way. Now
`run_remd.py` keeps the job id of the block in `<V>/block_jobid` and increments only for a new id
(sandbox-tested: fresh, requeued and successor jobs); the three counters were reset to 1 and each
variant's `block_jobid` set to its running job, so every chain gets the current block plus four.
Backups `run_remd.py.bak_pre_requeue_20260930`, `submit_remd.sh.bak_pre_requeue_20260930`,
`<V>/block_count.bak_pre_requeue_20260930`. The aborted attempts' output is kept as short
`output_previous_N` groups (100-287 frames), all finite (`checks/glpg_fragments_20260930.log`).
`submit_remd.sh` now carries the training chain's node exclusions
(`midway2-0003,[0010-0011],0037,[0342-0345]`) from the next block; the failures were simultaneous on
several nodes, so they look like cluster events, but 0010 and 0037 hosted two of them.

### ff30_basin: finished, released ff_3.0 (2026-09-28 to 09-30)

**What it was** (plan.md, findings 1.8-1.14): ff2.1's own FF2 workflow, plus 158 trained basin
offsets on 60 Ramachandran coil maps (GLY|X alpha_R, alpha_L, beta; GLY|GLY helix and beta, tied to
their mirrors; X|right|PRO alpha_R, beta), updated once per epoch by a damped MAP Newton step; the
other 780 maps NDRD; 46 proteins held out of the offset update. Steps 0-113 (epochs 0-5), links
49125955 -> 49131476, no failed step after the third rewind.

* **Where:** `$P=/project/trsosnic/yinhan/upside2-md-mdw2`, run dir `$P/training/ff30_basin`; logs
  `condiv-train_<jobid>.out`, gates `ff30b-gate_<jobid>.out` and `gate_step{76,95,114}.txt`;
  checkpoints `run_output/epoch_EE_minibatch_MM`, round libraries `run_output/rama_round_0{0..6}.dat`,
  per-protein `*.divergence.pkl` kept in every step directory; the release files in
  `release_20260930-022510/`.
* **Gates:** step 76 and 95 not converged (rama p 0.0011, 0.0006), step 114 converged (every group
  p > 0.005, rama 0.0158). The rama pass is prior-limited drift of the X|right|PRO offsets, not
  closure (findings 1.14).
* **`run_output/rama_rounds.txt`** (mismatch over the trained maps, training / held out; largest
  step; largest offset; GLY|X mean alpha_L - alpha_R offset): round 1 0.0403 / 0.0814, 0.44, 0.44,
  +0.060; 2 0.0408 / 0.0861, 0.40, 0.75, +0.132; 3 0.0365 / 0.0777, 0.34, 1.09, +0.177; 4 0.0348 /
  0.0712, 0.27, 1.36, +0.221; 5 0.0364 / 0.0826, 0.27, 1.63, +0.260; 6 0.0338 / 0.0827, 0.26, 1.89,
  +0.306. The held-out figure sits at its sample-size floor (findings 1.12).
* **Release (gate log):** hb [-1.878 -1.872 -1.798], dhb -0.617, bb scale -0.351, sheet mean 0.247;
  glpG round-trip gate 3.6e-15 on a pristine ff_2.1 seed; the live seeds changed by rama_pot 3.42,
  hbond 0.19, pair 17.3 (max).
* **Pre-proline analysis on this run:** `/project/trsosnic/yinhan/checks/prepro_*.py`, logs
  `prepro_{leverage,control,left}_20260930.log` (findings 1.14).
* **Rewound three times on 2026-09-28**, each time recomputing round 1 from the epoch-0
  simulations: `rewound_20260928_0851/` (Dirichlet log-ratio step moved empty basins by up to 1.76
  nats; replaced by the MAP step), `rewound_20260928_1225/` (GLY|GLY sheet entry not mirror-symmetric,
  findings 1.11), `rewound_20260928_1608/` (all 840 maps noise-limited; reduced to the 60, findings
  1.12-1.13; all-map step-18 checkpoint `checkpoint_epoch00_mb18_all_maps.pkl`). Pre-change files in
  `$P/training/backup_pre_{basin,ggsheet,reduced}_20260928/`.
* **Cluster flags live in `<run>/slurm.args`** (the midway2 broadwl line with node exclusions);
  `env.sh` picks midway2's tree venv or the `/beagle3` venv on midway3.
* **Recovery design that worked unattended:** every link queues its `afterany` successor before it
  trains; `--no-requeue`; the successor resumes from the newest checkpoint after deleting the
  half-written step; a worker whose `srun` never started is relaunched (twice at most); three links
  from the same step stop the chain. The gate job itself has no successor.
* **`sbatch --test-only`'s start estimate is no guide to the real wait**: it predicted 13:36 on both
  clusters, then the midway3 submission sat PENDING (Resources) while midway2 started at once.
* **`srun: Job credential expired` at launch is a race, not a node fault**: the first link lost 4 of
  24 workers 2-3 min after the 24 steps were issued together, on nodes whose other steps ran fine.
* **Steps slower than the engine accounts for** (2026-09-29): 300-450 s extra per step, from a
  different slowest worker each time, engine time unchanged; most likely `/project` I/O.

**`/beagle3` is badly degraded for small-file writes (measured 2026-09-28 ~01:10).** 200 one-line
files: `/beagle3` 51 s from midway3 and 503 s from midway2, against 0.06-0.27 s on `/project` and
home. `pip install torch==2.6.0` into the shared `/beagle3` venv took from 00:24 to ~02:15 for that
reason; it finished, and midway3 can now run the trainer.

**`env.sh` must not tell the clusters apart by `/software/modules/init/bash`: it exists on midway3
too.** The first version did, so on midway3 it activated the tree's venv, whose interpreter lives in
midway2's `/software`, and silently fell through to `~/.bin/python3` (torch 2.8, numpy 2.2). It now
tests whether the tree venv's interpreter exists. Verified 2026-09-28: midway2 gets the tree venv and
midway3 the `/beagle3` venv, both torch 2.6.0+cpu, numpy 1.23.5, scipy 1.13.1, tables 3.8.0 and this
tree's engine. `env_shared.sh` uses the same test only to choose a module init, which is harmless.
Never pipe `source env.sh`: a pipe runs it in a subshell and the environment is lost.

**ff30 (full-map glycine row) CANCELLED 2026-09-28 00:25 at step 223.** Its glycine group failed
every gate from 76 to 209 (p = 0) and its GLY|X handedness had
drifted back to -0.82. The run dir (54 G) and the last checkpoint `ff30/run_output/epoch_11_minibatch_13`
are kept for comparison; extract with the old extractor,
`$P/training/backup_pre_basin_20260928/extract_ff.py`, since the current one reads basin offsets.
Two lessons from its chain carry over to every run:
* **A failed worker fails the step** (`ConDiv.py` raises): step 81 once ran on 5 of 24 proteins after
  an `srun: Job credential expired`, and a partial step is a different objective.
* **Node exclusions** in the midway2 `slurm.args` carry midway2-[0010-0011] (three NODE_FAILs in one
  campaign) and midway2-0037 (two links in three hours on 2026-09-26). Two later link failures
  (2026-09-27 00:45, 03:30) named no node; nothing was excluded for them.

### The FF2 trainer on midway2 (2026-09-24)

`training/ConDiv.py` is the FF2 dual-target trainer (findings 9v, plan.md); the old FF1-form trainer
and helpers are in `$P/backup_training_ff1form/`.

**Measured step cost, full protocol** (`worker_test/`, job 49056799): 5vhg, 150 residues, the
largest in the set, 1242 s; 1ean, 114 residues, 897 s. A step waits for its slowest worker, so
**~21-22 min per step**; `STEPS_PER_LINK = 90` fits a 36 h wall. Both unfolded properly on the SI's
ladder (mean Rg 14 -> 46 A and 13 -> 36 A across T = 0.8-1.1), so the DSE target is real.

**What `validate_ff.sh` does at release** (run by the gate at convergence, see ff30_basin above):
extract through the run's own `expand_param`, back up and overwrite `parameters/ff_3.0` in `$P`
and in the /beagle3 deployment (md5-verified), move the superseded `runs/*_ff_3.0` benchmark
directories to `runs_superseded/`, submit the 32 Peng arms (~7 days), gate the glpG patch on a
pristine ff_2.1 seed, patch the 4 live seeds, clear their replicas and submit the 4 REMD chains
(~7.5 days). It submits to broadwl, so it runs from midway2.

**`bench_run.py` changed 2026-09-24** (backup `bench_run.py.bak_ff1form_20260924`): the type-0
burial override for ff_3.0 and the `rama3.dat` fallback are gone, since the new ff3.0 is FF2-form.
Nothing may benchmark the old FF1-form ff_3.0 with it.

### A link killed by SIGBUS on several nodes at once is a transient `/project` outage

On 2026-09-23 all 12 workers of one step died with **signal 7 (SIGBUS)** across four nodes, with
healthy physics, no traceback and 0-byte `*.output_worker` files, and the `afterany` successor died
in the same second it started, before it could create its `.out`. Disk space was checked and was not
the cause. A bad node kills its own 3 tasks, not 12 across four nodes, so nothing was excluded; the
chain logic needed no change and recovery was to resume. No GPFS log was available, so this is
inferred from symptoms.

### Disk: from midway2, `df` on the subdirectory is the ONLY number you get for `/project`

```
df -h /project/trsosnic     -> the FILESET (size, used, free), correct
df -ih /project/trsosnic    -> fileset inodes
df -h /project              -> 6.3P total, the whole device, useless
rcchelp quota               -> home, scratch, project2 group only; no /project or /beagle3 row
mmlsquota                   -> "File system project is not known", the GPFS client is not here
```

`df` on the **subdirectory** reports the fileset because GPFS `--filesetdf` is on; `df` on the
**mount point** reports the 6.3 PB device. That distinction is the whole trap.

**Current headroom, and the trend, which is the part that matters:**

| fileset | 2026-09-09 | 2026-09-23 | 2026-09-28 | note |
|---|---|---|---|---|
| `/project/trsosnic` | 1514 G free | 445 G free | 953 G free; **965 G free 09-30 12:00** | `training/` is ~108 G, of which `ff30` 54 G; the four glpG variant directories were emptied at the release and are refilling (~150 G each at the end of the last campaign) |
| `/beagle3/trsosnic` | - | 1.4 T free (4.2 T of 5.5 T) | 1.4 T free | badly degraded for small-file writes on 09-28, see §1 |
| `/project2/trsosnic` (group) | - | 1.45 T of a 1.49 T soft quota, **97%** | - | nothing in these campaigns writes there |
| midway3 home | - | 28.6 of 30 G | 21 G | |

`/project` is the one to watch. The four glpG variant directories in `popepopg_REMD_mdw2` hold
~150 GB each and `popepopg_REMD` holds another 638 G; those two trees plus `NP-1AO6` at 491 G are
1.75 T of the 3.5 T used.

### A nested sbatch inherits `SLURM_*` from the job that calls it

A submitting script that requests `--mem` passes `SLURM_MEM_PER_NODE` to the job it submits; if that
job's `srun` requests `--mem-per-cpu`, every worker launch dies instantly with
`srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive.`
This broke a training chain's only resubmit path on 2026-09-06, unseen because every earlier job had
been submitted from an interactive shell. **Any script that submits another job must unset those
three variables after its `#SBATCH` block**, or match the child's memory-request type. Unsetting them
does not change the running job's own allocation.

### Slurm snapshots the batch script at submission

* **A requeue reruns the script verbatim with its original arguments**, so a requeued training link
  resumed from the stale checkpoint it was submitted with and `main_loop`'s `rmtree` deleted the
  newer ones. Chain scripts carry `#SBATCH --no-requeue` and resolve the newest checkpoint on disk at
  run time. Verify with `scontrol show job <id> | grep Requeue` (must be `Requeue=0`).
* **An edit to a `.sbatch` does not reach jobs already queued.** After editing any script in a chain,
  cancel and resubmit the queued jobs that use it; check what a queued job will run with
  `scontrol write batch_script <id>`.

## 1b. Data and decisions that outlive their jobs

The campaign sections these came from were removed on 2026-09-28 (they described jobs long finished
as running; they are in git history). What must not be forgotten:

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
| timestep | **0.001**, freely settable at runtime | **0.009, HARD-LOCKED** by `/input/brownian/numerical_time_step`; `martini_brownian.cpp:100` throws on mismatch. Friction is tuned against it for lipid D=11.5 µm²/s — **do not change it** |
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
```

Both read only rotated `output_previous_*` groups (never the live `/output`) and skip
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
  glycines) is invalid twice over: sign-flipped, and ff_3.0's GLY|X maps are asymmetric by design,
  so it reports every ff_3.0 seed BROKEN. Do not gate a submission on it.

**ff_3.0, first ~10 h (2026-09-30, `checks/glpg_tm_windows_20260930.log`):** seed TM1 = TM4 = 1.00; at
T 0.70 TM1 0.90-1.00 in all four variants, TM4 0.97-1.00 in three; soft spots are 79ALA_S115T TM4
0.82-0.86 over its last three groups at T 0.70 (residue 148 at 0.55) and 79HIS_S115T TM4 0.95 -> 0.84
at T 0.80 (residue 136 at 0.29). TM4's helical glycines already leave the helix for phi > 0 in
some replicas (GLY143 36% of frames in 79ALA_S115T at T 0.70; `checks/gly_tm4_flip.py`): ff_3.0's
glycine term favours alpha_L at every glycine, most at GLY136 (findings 1.14). Too early to judge;
run `gly_tm4_flip.py` with the TM check after block 1.


## 5b. How to check health CORRECTLY (general bond/energy check)

**`isfinite` is not a health check.** At a real glpG failure the environment coordinates were
**±4.65e12 Å** — numerically finite, physically destroyed. And in a forced NP tear the protein reached
**431 broken bonds with the potential still finite at +3e5**, so no energy-based test fires at all.

Use the broken-bond **count** (torn 279–431). The "healthy 0–2" figure holds for the cold rungs
only. Measured on glpG 2026-09-30 (`checks/glpg_cn_control_20260930.log`, every 5th frame): at
T 0.70 no frame has 3 or more C-N above 2.0 A, at T 0.80 ~0.1% do, and at T 0.88-0.90 1.5-2.7% under
ff_3.0 (`inner_steps` 4, max 10 in one frame) against 1.8-5.1% in the pre-ff3 campaign (`inner_steps`
1, max 13). So local transient tearing at the hot rungs predates ff_3.0; it is a known open
integrator question for the hybrid, not a force-field regression; it does not accumulate (sampled
frames return to 0, and the current logs have no ROLLBACK). The Peng arms (pure Upside) never exceed 1 in any frame.

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
`rolled back M/N replicas; continuing chain`. The chain does NOT terminate. The NaN output is rotated
to `output_previous_N` as normal history; the rolled-back replica restarts the next chunk from its
pre-chunk positions. A replica that repeatedly blows up gets rolled back repeatedly; it does not get
dropped from the ladder. Watch for high rollback counts on the same file — that indicates a replica with
a persistent physics problem that won't self-correct.

**Rollback mechanism**: before each chunk `run_remd.py` snapshots `/input/pos` of every replica.
On NaN detection it overwrites the last `output/pos` frame and `output/potential[-1]` with
the pre-chunk values so that `reseed()` on the next iteration picks up the clean state.

**These driver scripts are NOT in git.** They live on the cluster at
`~/project/yinhan/popepopg_REMD_mdw2/run_remd.py` (midway2) and `~/project/NP-1AO6/run_np_prod.py`.
No version history exists for them. Edit directly on the cluster.
A running job keeps the version it loaded at start; edits take effect at the **next block**.

---

## 7. If the chain terminates — manual rollback procedure

**glpG (if the chain stops rather than rolling back):** patch last output frame of each NaN file
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
- **glpG NaN propagation via REMD exchange.** A single blow-up in one replica spreads to every replica via exchange within ~60 steps (IEEE 754: `NaN < 0.f` = false). The NaN cascade fix (`!isfinite(lboltz_diff)`) and the per-chunk rollback driver address this.
- **A green exit code means nothing** for a self-submitting REMD job. Check the log for DESTROYED/ROLLBACK counts and verify physical observables.
- **Midway3 home quota is 30 G** (21 G used on 2026-09-28, 28.6 G on 09-23). Jobs can fail oddly if home fills.
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

Lessons kept from the campaigns removed on 2026-09-28, one line each:
- **Cancel a chain's PENDING successor before its running link**, then check `squeue`; the other
  order starts the successor.
- **Delete a half-written minibatch directory before resuming**: a stale `divergence.pkl` is reused.
- **Hand-check success-only end-of-chain branches**, which never run until the end.
- **glpG: patch the LIVE seeds** (they carry `inner_steps = 4` on `/input/brownian`); prove the method
  by an ff_2.1 round trip on a pristine seed. A release must remove each ~150 GB variant directory
  before resubmitting it.
- **Identify a `.up`'s force field by least-squares match of its baked tables to `sidechain.h5`**
  (scale 1.000000), not by hash. `martini.h5` comes from ff_2.1 by design.
- **Pre-create `runs/<tag>/input` before a mass submission** on `/beagle3` (a makedirs race); a
  partial aggregate mean is not a result.
- **Runs branched from a shared checkpoint do not test reproducibility**; branch at step 0.
  Independently trained force fields are different Hamiltonians and are never pooled.
- **midway2's site GROMACS is unusable**: use `/project/trsosnic/yinhan/gmx2024`
  (`gmx_build/build_gmx.sbatch`). The HDX pipeline runs only on midway3 (`.venv_el8_py311_bak`, with
  pymbar and matplotlib); a renamed venv's `bin/activate` still hardcodes its old `VIRTUAL_ENV`.
- **The shared venv is Python 3.9 / NumPy 1.23**: `str | Path` needs `from __future__ import
  annotations`.
- **Live dynamics are shown by internal Kabsch CA-RMSD growth** (frozen = 0.000 A), never by reading
  `current_stage`; read only `output_previous_*` of a live `.up`, with `HDF5_USE_FILE_LOCKING=FALSE`.
- **glpG REMD: skip `output_previous_0`** (the frozen seed block), and run `check_seeds_current.py` (six
  helical GLY, sym_err < 0.001) before any resubmission.
- **NP:** the box must exceed the molecule's maximum extent plus the 12 A cutoff; `NP_DT = 0.001`; the Rg
  in `np.<jobid>.out` is not minimum-imaged (measure adsorption by backbone contacts with the MPA shell,
  atoms 3750:3950); a restart frame must pass a stricter test than detection; `NP_TEMP` is Peng's FF2
  calibration and must be redone for a new force field.
- **HDX figures:** keep continuous off-scale excursions and `Y_LIMITS=(-20,30)`, and never render
  censored amides as carets or bounds (user-directed). Cluster HDX scripts drift from the repo;
  re-upload before trusting results.
- **`sacct` on midway3 can hold zombie RUNNING rows** (53233848 and 53233852 from 2026-08-11 still
  show under `-S <today>`, while `sacct -j` says `COMPLETED`); `squeue` is the authority.
- **Exclude midway3-0014.**
