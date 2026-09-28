# Remote jobs on midway2/midway3 — status and handbook

**2026-09-28 16:35: ff3.0 retrain RUNNING on midway2 (`training/ff30_basin`, §1), with basin
offsets on 60 Ramachandran maps only (GLY|X, GLY|GLY, X|right|PRO; 158 offsets; the other 780 maps
NDRD), resumed from step 19 at 16:10 after the third rewind of 09-28. The full-map glycine chain
`ff30` was cancelled at step 223. At convergence the gate installs the rotamer-BP fix (§0c),
releases ff_3.0 and submits the Peng benchmark and the glpG chains unattended. BP validation: 8 of
16 arms complete. `/beagle3` was badly degraded for small-file writes on 09-28 (§1).**

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
| `/project/trsosnic/yinhan/upside2-md-mdw2` (midway2) | **still old binary**; installs at the ff3.0 gate, below |
| `/beagle3/trsosnic/yinhan/upside2-md` (midway2 + midway3) | **still old binary**; installs at the ff3.0 gate, below |

**Static-frame result (local, 2026-09-26).** 8 training proteins (50-146 res), ff2.1, 91 frames each
(native + T 0.80/0.95/1.10), old vs fixed library vs a tol=1e-6 reference. The early stop fires in
~45% of frames but the error is small: |dE| <= 0.029 E_up, force error 7e-5 (median) / 3e-3 (max) of
the RMS force, ConDiv rotamer-gradient error ~1e-4 relative; the energy bias is one-signed (converged
is higher, mean 0.0009 E_up). Conclusion given to the user: ff2.1 and ff3.0 do not need retraining.

**Simulation validation for Tobin: jobs 49120982-49120997 (16 arms, running).** ff_2.1, native-start
14-replica REMD with the Peng benchmark protocol (Table S2 ladder and duration, dt 0.009, frame 100),
proteinG / homeodomain / WWdomain / NTL9 x {bug, fix} x seeds {1, 2}. Trajectories of the two
binaries decorrelate within a few hundred time units, so the test is statistical: bug-minus-fix
against the seed-to-seed spread (z = delta / sigma_seed). All 16 checked healthy at the first 100
frames (C-N 1.33-1.34 +- 0.13 A as in the ff_3.0 benchmark, KE/1.5kT 0.99-1.00). Expected done
~1.5 days (WW, NTL9) and 2.5-3 days (proteinG, homeodomain) from 2026-09-26 18:55.
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

**Deployment: armed in the ff3.0 gate, NOT yet fired (re-armed 2026-09-28 09:00).** The cluster copy
of `upside2-md-mdw2/training/gate_or_continue.sh` carries ONE TEMPORARY line in its converged branch,
before `validate_ff.sh`: `bash /project/trsosnic/yinhan/bpfix_stage/install_bpfix.sh`, so
`ff30_basin` trains on the old binary and all its validation runs on the fixed one.
`bpfix_stage/gate_or_continue.sh.orig` is the repo's current gate (md5 0a653ade), which the script
restores afterwards; the 09-26 `.orig` is in `$P/training/backup_pre_basin_20260928/`. **It must not
be run by hand while `ff30_basin` trains**, since it swaps `$P/obj` under the run.
* staged builds `/project/trsosnic/yinhan/bpfix_stage/mdw2/obj` and
  `/beagle3/trsosnic/yinhan/bpfix_stage/beagle3/obj`: each tree's own source plus the one line, on the
  tree's filesystem so install is a rename; the same recipe unpatched (`obj_ref`) reproduced each
  tree's binary bitwise
* the script checks both trees first (source against `src_manifest.md5`, staged binaries against
  `obj_manifest.md5`) and installs nothing on any mismatch; keeps old binaries as
  `obj/*.bak_pre_bpfix_<stamp>`; writes `<tree>/BPFIX_INSTALLED`; always restores the gate script
  from `bpfix_stage/gate_or_continue.sh.orig`. Sandbox-tested: install, repeat no-op,
  refuse-on-change
* **do not edit `src/` in either tree until it has run**, or it will refuse
* afterwards: `cat <tree>/BPFIX_INSTALLED` and the `[bpfix]` lines of its output; then update the
  status table above
* the BP test is unaffected: its arms use `obj_bug`/`obj_fix`, and the deployed library its Python
  imports is replaced by rename, so a loaded copy keeps its inode

---

## 1. Current jobs

Snapshot **2026-09-28 17:25 CDT, verified live against `squeue` on midway2** (midway3 last checked
10:40, nothing of ours queued or running there). Finished and cancelled rows are deleted; only
lessons worth reusing are kept, below the table.

| JobID | what | where / state | next action |
|---|---|---|---|
| **49127867** | **ff3.0 retrain, basin offsets on 60 maps (158 offsets)**, `training/ff30_basin`, steps 19 -> 76 (epochs 1-3), resumed from step 19 after the reduction to the trained set (third rewind, below) | R since 16:10, midway2-[0247-0258], broadwl; steps 19-21 done in 1253-1263 s (24/24 workers each, no failures, DSE 24/24, median RMSD 0.91-1.00 / 2.69-2.99 A), step 22 running since ~17:13; live step-19 worker maps equal the new library exactly, residues reading no trained map exactly NDRD | 54 steps at ~1258 s ~19 h, inside the wall; step 76 ~2026-09-29 12:00, then `after_training.sbatch` gates it |
| 49127868 | its insurance successor | PD, `afterany:49127867` | resumes only if 49127867 dies |
| **49125873, 49125874, 49125875, 49125877, 49125878, 49125879, 49125880, 49126021** (8 running) | **BP stopping-test validation**, 16 arms `bp_<prot>_<bug\|fix>_s<1\|2>` (§0c); 8 of 16 `COMPLETE` (NTL9 and WWdomain, all 4 each, NTL9 at ~28 time units/s); still 8 at 17:25, all 8 remaining arms running. Remaining at the start of this link: proteinG 1.38-1.55 M of 3.37 M, homeodomain 1.50-1.52 M of 3.12 M | R, 1 node each, self-resubmitting every ~29 h under the same names | measured 13:42 on this link's own rate (proteinG 17.6-18.9, homeodomain 15.2-15.4 time units/s), proteinG finishes 2026-09-28 20:45 to 09-29 00:40 and homeodomain 09-29 03:15-04:05, all inside this link (homeodomain needs 27.1-27.9 h of the ~29.2 h per link); then run `analyse_bp.py` and report to Tobin. They write to `/beagle3`; see the degradation note below |

### ff30_basin: ff3.0 from ff2.1 with per-pair basin offsets (started 2026-09-28)

**What it is** (plan.md, findings 1.8-1.13): ff2.1's own FF2 workflow, plus 158 trained basin
offsets on 60 Ramachandran coil maps: GLY|X (alpha_R, alpha_L, beta), GLY|GLY (helix and beta, tied
to their mirrors) and X|right|PRO (alpha_R, beta); the other 780 maps are NDRD unchanged. Each map is
its own parameter set, updated once per epoch by a damped Newton step matching free to
native-restrained basin populations, no DSE term on them. 46 proteins are held out of that update.
Epoch 0 ran at ff2.1 with every offset zero; round 1 was recomputed for the 60-map set from those
simulations. `verify_rama_basin.py` PASSES on the run (2026-09-28 16:05).

* **Where:** `$P=/project/trsosnic/yinhan/upside2-md-mdw2`, run dir `$P/training/ff30_basin`
  (`init_param` = ff_2.1 md5-verified; `upside_input` hardlinked from `ff30`, its `rama.dat` =
  `parameters/common/rama.dat`). Log `ff30_basin/condiv-train_<jobid>.out`; per-worker errors only in
  `run_output/epoch_*/<code>.output_worker`; checkpoints `run_output/epoch_EE_minibatch_MM`.
* **What to watch:** `run_output/rama_rounds.txt`, one line per epoch: sites of the trained maps,
  free-native mismatch over those maps (the fraction of residue time in a different basin) for
  training and held-out proteins, largest step and offset, mean alpha_L - alpha_R offset over the
  GLY|X maps. Round 1: 7,528 / 882 sites, mismatch 0.0403 training / 0.0814 held out, largest step
  0.44, GLY|X mean(dL-dR) +0.060. Round 2 comes after step 37 (~2026-09-28 23:00). **Read the
  held-out mismatch against its sample-size floor**, not as falling or stalling: with only 46
  proteins it sits near its floor whatever the offsets do (findings 1.12).
* **Cluster flags live in `ff30_basin/slurm.args`** (now the midway2 broadwl line with the node
  exclusions) and every submission reads them; to move the run to midway3, write
  `--partition=caslake` there. `env.sh` picks the Python by cluster: midway2 the tree's `.venv`,
  midway3 the `/beagle3` shared env.
* **At step 76** `after_training.sbatch` runs `gate_or_continue.sh ff30_basin ff_3.0 13`: not
  converged -> one more epoch through `slurm.args`, up to 13 epochs (step 247), then it stops for
  review; **converged -> the armed line installs the rotamer-BP fix, then `validate_ff.sh` releases
  ff_3.0 to both trees and submits the 32 Peng arms and the 4 glpG REMD chains**, all unattended.
  Dry-run 2026-09-28 on the round-1 checkpoint: extraction (six files, GLY|GLY asymmetry 0), glpG
  round-trip gate 3.6e-15 on the pristine ff_2.1 seed, a live seed patched into scratch, and
  `sbatch --test-only` accepting a Peng arm, a glpG chain, the gate and a training link. The four
  variant directories it deletes hold only the one-block ff3.1 test replicas (28 each, 5.7-5.9 G).
* **Recovery without a login:** every link queues its insurance successor (`afterany`) before it
  trains, so a node failure, a wall kill or a failed step is resumed by the successor from the newest
  checkpoint, whose half-written step directory `main_loop` deletes first; `--no-requeue` stops Slurm
  replaying a link with stale arguments; a worker `srun` never started is relaunched; three links in
  a row starting from the same step stop the chain rather than loop. The one gap: the gate job itself
  has no successor, so if it dies the chain waits for a login (`sbatch $(cat slurm.args)
  after_training.sbatch` from the run dir).
* **Round 1 was rewound (2026-09-28 08:51).** Its log-ratio step with Dirichlet pseudo-counts moved
  offsets by up to 1.76 nats in basins with no residues in either ensemble. The chain was stopped
  during step 20; steps 19-20 and the old-rule checkpoint, code, library and round log are in
  `ff30_basin/rewound_20260928_0851/` (outside `run_output/`, so the chain cannot resume from them);
  round 1 was recomputed from the same epoch-0 simulations with the MAP step (largest 0.44) using the
  run's own updated `rama_basin.py` and `ConDiv.py`, and training resumed from step 19.
* **Rewound a third time, 2026-09-28 16:08, to train only 60 maps (findings 1.12-1.13, plan.md).**
  Training all 840 maps was noise-limited on 456 proteins; literature and ff2.1's own error point to
  GLY|X, GLY|GLY and X|right|PRO, 158 offsets; the other 780 maps stay NDRD. 49127867's predecessor
  49126777 was stopped during step 29; steps 19-29, the old code and library, `rama_rounds.txt`,
  `.chain_starts` and the previous step-18 checkpoint (`checkpoint_epoch00_mb18_all_maps.pkl`) are in
  `ff30_basin/rewound_20260928_1608/`; the unfinished step 29's 315 replica files (42 G, no
  checkpoint) were deleted. Round 1 was recomputed from the epoch-0 simulations with the run's own
  `rama_basin.update` and `ConDiv.finish_round` (`/project/trsosnic/yinhan/checks/recompute_round1.py`;
  same held-out proteins; trained offsets equal the previous ones to 2e-16; every non-rama part of
  the checkpoint identical); `extract_ff.py` releases from it, `verify_rama_basin.py` PASSES. Replaced
  files in `$P/training/backup_pre_reduced_20260928/`; deployed md5s: `ConDiv.py` 83e309da,
  `rama_basin.py` 0dc41a35 (both also in `run_output/`), `verify_rama_basin.py` 1657ef0e, `README.md`
  94f507a4. The round log's mismatch is now over the trained maps only (0.0403 / 0.0814 held out),
  not comparable with the earlier all-map figure.
* **Rewound again, 2026-09-28 12:25, for the GLY|GLY fix (findings 1.11).** The engine's map for a
  glycine between glycines was not mirror-symmetric, through the raw NDRD sheet GLY|GLY entry.
  `rama_basin.py` now symmetrises that entry too, by the probability mean for coil and sheet.
  49126332 was stopped during step 29; steps 19-29, the old code, the old `rama_round_01.dat`,
  `rama_rounds.txt` and `.chain_starts` are in `ff30_basin/rewound_20260928_1225/`.
  `rama_round_01.dat` was rewritten from the step-19 checkpoint (same offsets, bitwise identical
  outside GLY|GLY, GLY|GLY asymmetry 66.5 -> 0), the round log's "X|GLY" label corrected to GLY|X,
  `verify_rama_basin.py` PASSES on the run, and three configs with GGG or terminal-GG glycines built
  through `upside_config` are exactly symmetric. Replaced files are in
  `$P/training/backup_pre_ggsheet_20260928/`. The two unfinished steps' replica trajectories and
  configs (step 29 here, step 20 in `rewound_20260928_0851/`, 87 G, no worker of either had
  finished, no checkpoint, referenced by nothing) were deleted at 12:45; logs and parameter files
  are kept and each directory has a `README.txt`. Deployed md5s: `ConDiv.py` fd9115d5, `rama_basin.py`
  0b1bec68 (both also in `run_output/`), `extract_ff.py` 9f9d3d3b, `verify_rama_basin.py` 12cbb4d6,
  `README.md` f8445629.
* **Deployed 2026-09-28**, md5-verified against the repo: `ConDiv.py`, `rama_basin.py`,
  `verify_rama_basin.py`, `extract_ff.py`, `convergence_gate.py`, `gate_or_continue.sh`,
  `train_chain.sbatch`, `env.sh`, `README.md` in `$P/training/`. The replaced files, the removed
  `rama_gly_gradient.py` / `verify_gly_gradient.py` and the stale `$P/py/rama_gly_gradient.py` are in
  `$P/training/backup_pre_basin_20260928/`.
* **First try on midway3 (59635685) was cancelled after 10 min PENDING (Resources)**, as the user
  asked, and resubmitted on midway2, where it started at once. `sbatch --test-only` had predicted a
  13:36 start on both clusters, so its estimate is no guide to the real wait.
* **The first midway2 link (49125951) lost 4 of 24 workers at launch** to `srun: error: Task launch
  for StepId=... failed: Job credential expired`, 2-3 min after the 24 steps were issued together,
  on nodes whose other steps ran fine: a launch race, not a node fault. It was cancelled before the
  step completed, with its successor, and the trainer now relaunches a worker whose `srun` never
  started it (`<code>.srun` says `Task launch`), as soon as that is seen, up to twice; a worker that
  ran and failed still fails the step. The run was re-initialised with that code (no step had
  finished; its log is in `ff30_basin/aborted_20260928_0126/`). A relaunch shows in the link log as
  `<code> never started (srun launch failed), relaunching`.

**`/beagle3` is badly degraded for small-file writes (measured 2026-09-28 ~01:10).** 200 one-line
files: `/beagle3` 51 s from midway3 and 503 s from midway2, against 0.06-0.27 s on `/project` and
home. `pip install torch==2.6.0` into the shared `/beagle3` venv took from 00:24 to ~02:15 for that
reason; it finished, and midway3 can now run the trainer. The BP arms write their runs to
`/beagle3` and may be slowed.

**`env.sh` must not tell the clusters apart by `/software/modules/init/bash`: it exists on midway3
too.** The first version did, so on midway3 it activated the tree's venv, whose interpreter lives in
midway2's `/software`, and silently fell through to `~/.bin/python3` (torch 2.8, numpy 2.2). It now
tests whether the tree venv's interpreter exists. Verified 2026-09-28: midway2 gets the tree venv and
midway3 the `/beagle3` venv, both torch 2.6.0+cpu, numpy 1.23.5, scipy 1.13.1, tables 3.8.0 and this
tree's engine. `env_shared.sh` uses the same test only to choose a module init, which is harmless.
Never pipe `source env.sh`: a pipe runs it in a subshell and the environment is lost.

**ff30 (full-map glycine row) CANCELLED 2026-09-28 00:25 at step 223** (49125428 and successor
49125429, and its monitor 49120529 / 49125766, which watched only ff30; `STATUS.md` is stale from
then). Its glycine group failed every gate from 76 to 209 (p = 0) and its GLY|X handedness had
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

`training/ConDiv.py` is now the FF2 dual-target trainer (findings 9v, plan.md). Deployed md5-matched
to `$P/training/`, with `extract_ff.py`, `patch_glpg.py`, `validate_ff.sh`, `check_converged.py`,
`train_chain.sbatch`; the old trainer and helpers are in `$P/backup_training_ff1form/`.

**Measured step cost, full protocol** (`worker_test/`, job 49056799): 5vhg, 150 residues, the
largest in the set, 1242 s; 1ean, 114 residues, 897 s. A step waits for its slowest worker, so
**~21-22 min per step**; `STEPS_PER_LINK = 90` fits a 36 h wall. Both unfolded properly on the SI's
ladder (mean Rg 14 -> 46 A and 13 -> 36 A across T = 0.8-1.1), so the DSE target is real.

**What `validate_ff.sh` does at release** (now run by hand, see ff30_basin above): extract through
the run's own `expand_param`, back up and overwrite `parameters/ff_3.0` in `$P` and in the /beagle3
deployment (md5-verified), move the superseded `runs/*_ff_3.0` benchmark directories to
`runs_superseded/`, submit the 32 Peng arms (~7 days), gate the glpG patch on a pristine ff_2.1
seed, patch the 4 live seeds, clear their replicas and submit the 4 REMD chains (~7.5 days). It
submits to broadwl, so it runs from midway2.

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

| fileset | 2026-09-09 | 2026-09-23 | 2026-09-28 | note |
|---|---|---|---|---|
| `/project/trsosnic` | 1514 G free | 445 G free | **953 G free** at 12:45 (inodes 21%) | `training/` is ~108 G, of which `ff30` 54 G; a finished `ff30_basin` step keeps 35 MB, a running one holds ~40 G of replicas until its workers finish, and a step killed mid-run leaves them behind |
| `/beagle3/trsosnic` | - | 1.4 T free (4.2 T of 5.5 T) | 1.4 T free | badly degraded for small-file writes on 09-28, see §1 |
| `/project2/trsosnic` (group) | - | 1.45 T of a 1.49 T soft quota, **97%** | - | nothing in these campaigns writes there |
| midway3 home | - | 28.6 of 30 G | 21 G | |

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

## 1b. Data and decisions that outlive their jobs

The campaign sections these came from were removed on 2026-09-28 (they described jobs long finished
as running; they are in git history). What must not be forgotten:

* **Data that exists only on the cluster, do not delete:**
  * AWH glycine surfaces, `/project/trsosnic/yinhan/gly_peptides/` (rep1, `awh_amber99sb-ildn_rep2/`,
    `evidence_diffusion_bug/`).
  * `popepopg_REMD_mdw2/run_remd.py` and `NP-1AO6/run_np_prod.py` exist only there, with no git history.
  * `popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`, the pre-ff3 TM baseline.
  * HDX: keep `hdx/` (pre-fix) and `hdx_postfix/`; ff3 results are in `popepopg_REMD_mdw2/<V>/hdx_10k/`.
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
- **Slurm snapshots the script at submission**: check a queued job with `scontrol write batch_script
  <id>`, and hand-check success-only end-of-chain branches, which never run until the end.
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
- **`sacct` on midway3 can hold zombie RUNNING rows**; `squeue` is the authority.
- **Exclude midway3-0014.**
