# Remote jobs on midway2/midway3 — status and handbook

Snapshot: **2026-09-08 ~15:30 CDT. The RCC storage outage is OVER. ff3.0 training runs on TWO hosts — midway2 `48988330` and rockfish `30725720`, both at step ~355/500, finishing 2026-09-09 afternoon. The install chain is RE-ARMED as `48994406` (`ff-check-cont`, afterany:48988330) and now hands off to a retargeted arm test that compares the two trained force fields instead of two coverage recipes. Measured 2026-09-08: the two runs differ by 15-17% rms on every trained table, a third of the whole training signal, so that test is not a formality. glpG production is idle (all four variants COMPLETED 2026-09-05), so the launch has nothing to cancel. OWED: extract rockfish's step-500 force field into `parameters/ff_3.0_trained_rf/` with a matching `STEP` file, or the test runs one-armed.**
Written so a fresh session can pick up cold. Everything needed to connect, check health correctly,
and react to a failure is here. Job state below is live; superseded jobs are not listed, only
summarised in §8 where they carry a lesson.

---

## 0. Connect first (needs a Duo push on the user's phone)

Key-based auth is NOT enabled; password + Duo is the only method. The ControlMaster socket expires
roughly hourly, so expect to redo this most sessions.

**midway2** (POPE/POPG REMD campaign):
```bash
ssh -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu   # alive?
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_master.exp    # if not: USER MUST APPROVE DUO
ssh -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<command>'
```
`~/project` on midway2 is a symlink to `/project/trsosnic` (not `/project/trsosnic/yinhan/`).
Data is at `~/project/yinhan/popepopg_REMD_mdw2/`.
Python env: `source /software/modules/init/bash && module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1 && export HDF5_USE_FILE_LOCKING=FALSE`

**midway3** (NP campaign and glpG-DDM micelle):
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

## 1. Current jobs

Snapshot **2026-09-08 ~14:00 CDT (verified live against `squeue` on both hosts)**.

### midway2 — gly-sym core FF training

**RUNNING again.** The GPFS outage cleared and training resumed 2026-09-08 10:41 CDT.

| JobID | Name | State | Notes |
|---|---|---|---|
| **48988330** | `upside-gly-sym` | **RUNNING** on `midway2-[0027,0150,0260,0269]` since 2026-09-08 10:41 CDT, 1-12:00:00 wall (ends 2026-09-09 22:41 CDT) | Resumed from `run_output/epoch_08_minibatch_33/checkpoint.pkl` = step 338, "for 162 steps". Newest written checkpoint `epoch_09_minibatch_14` = **step 357** as of 15:20 CDT. 0 failures. Log `upside-gly-sym_48988330.out`, run dir `run_output/`. |
| **48994406** | `ff-check-cont` | PENDING (`afterany:48988330`) | The armed install chain. See below. |

It resumed into the **existing** `run_output`, not the fresh directory the recovery plan called for,
and the previously wedged inodes gave no trouble. `fs-probe 48988264` returned **COMPLETED 0:0**
after the 16 probes before it all read TIMEOUT; that is what unblocked the resume.

`.ff_installed` and `.production_relaunched` are absent, so no destructive step fired during the
outage and production was never at risk. `parameters/ff_3.0_trained/` does not exist yet.

**`continue_mdw2.sbatch` deliberately does not install the force field** — its own header says
"Deliberately NOT check_continue.sbatch: that script installs the trained force field". The install
scripts are all still present in `training/gly-sym/` (`check_continue.sbatch`,
`decide_and_launch.sbatch`, `extract_ff.py`, `patch_seeds.py`, `check_hybrid_up.py`,
`check_quota.py`, `push_progress.sh`, `tm_health.py`), but **nothing is armed to call them**.

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

### The install chain is re-armed, and what it now does

Armed 2026-09-08 as **`48994406`** (`ff-check-cont`, `--dependency=afterany:48988330`).
`48988387` (`continue_mdw2.sbatch`) was cancelled once the submission was confirmed: it only ever
continued training, which `check_continue.sbatch` also does, so keeping both was two chain drivers
racing for one `run_output`.

**`48994331` was the first attempt and had to be replaced** — a live instance of the snapshot trap
below. It was submitted before `check_continue.sbatch` gained the `FF_DIR/STEP` stamp, so it carried
the older snapshot; `run_arm_test.sbatch` is submitted at *runtime* and would therefore have picked
up the newer version that refuses to run without that file, and the chain would have died at the
handoff with nothing launched. The rule that catches this: **after editing any script in the chain,
resubmit every queued job that runs it.** The re-arm now greps the on-disk script for the feature
before resubmitting against it, and always submits the replacement before cancelling the original.

What fires when 48988330 ends:

```
48988330 ends --afterany--> 48994406 check_continue.sbatch
    step < 500 ?  resubmit srun_mdw2.sh for the remainder, re-arm, exit
    step = 500 ?  extract_ff.py    -> parameters/ff_3.0_trained/{sidechain,environment}.h5
                  write FF_DIR/STEP = 500
                  compare_ff.py    -> logs midway2 vs rockfish, gates nothing
                  back up 4 glpG seeds as .up.bak_pre_ff3_<stamp>
                  submit run_arm_test.sbatch
                      arm M = coverage nodes + midway2's ff        (always)
                      arm R = coverage nodes + rockfish's ff       (only if STEP matches)
                      12 h each, 28 replicas, then afterany -->
                  decide_and_launch.sbatch
                      decide_arm.py: health gate, then higher mean(TM1,TM4);
                                     inside 0.05 scatter midway2 wins as primary
                      check_quota.py 120 -> abort if short, nothing touched
                      verify the seed backups exist -> abort
                      inject coverage + patch the 4 real seeds -> abort, restore from backup
                      check_hybrid_up.py --require -> abort
                      archive $V/run.0.up -> $V/pre_ff3/, delete rungs 1-27, keep $V/hdx/
                      block_count=0, submit_remd.sh x4 (REMD_MAX_BLOCKS=4)
                      remove armtest_M / armtest_R, log quota
                      push_progress.sh  (EXIT trap: runs on every path)
```

**The safety property is unchanged and is now stronger:** every failure before the archive step
leaves the pre-ff3 data, the seed backups and the queue untouched. There is no running production to
cancel, so the old script's most destructive stretch no longer exists; what replaced it moves
`run.0` aside before deleting anything.

Scripts changed on midway2 (originals kept as `.bak_pre_ffcompare`, `submit_remd.sh.bak_pre_maxblocks`),
each uploaded and md5-verified end to end:

| file | change |
|---|---|
| `check_continue.sbatch` | step derived from the newest written checkpoint's directory name, not a file count; stamps `FF_DIR/STEP`; runs `compare_ff.py`; hands off to `run_arm_test.sbatch` |
| `run_arm_test.sbatch` | arms M/R differ by force field, not by node recipe; arm R optional and refused unless its `STEP` matches |
| `decide_arm.py` | M/R labels, midway2 wins ties as declared primary, honest one-armed n=1 verdict |
| `decide_and_launch.sbatch` | reads the verdict again instead of a hardcoded `ARM_B`; coverage always injected; archives `run.0` instead of deleting the ladder; dead `NP_BASE`/`NP_NODES` removed |
| `compare_ff.py` | new; array-by-array diff of two trained `sidechain.h5` |
| `submit_remd.sh` | `REMD_MAX_BLOCKS=4`, because 12 needs ~690 GB and the quota allows ~284 |

**Verified before arming, not assumed:** `bash -n` on all four sbatch scripts, `py_compile` on both
python files, `compare_ff.py` against `ff_2.1` versus itself (all `rel_rms` 0.0000) and against a
missing file (exit 1, the one-armed path), `decide_arm.py` on unbuilt arms (exit 2, `NO_WINNER`),
and `extract_ff.py` dry-run on rockfish writing a valid `sidechain.h5` from its live checkpoint.

### OWED: deliver rockfish's step-500 force field

midway2 cannot reach rockfish from inside a Slurm job, so this is a manual step, and without it the
arm test runs one-armed at n=1. When rockfish `30725720` reaches step 500 (`epoch_13_minibatch_05`):

```bash
RF=/scratch4/rherna21/ywang268/upside2-md-rf
D=$RF/training/gly-sym
ssh rockfish "cd $RF && source /etc/profile.d/modules.sh && module load gcc/9.3.0 hdf5/1.10.7 && \
  source $RF/.venv/bin/activate && export HDF5_USE_FILE_LOCKING=FALSE && \
  python3 $D/extract_ff.py $D/run_output_ff3/epoch_13_minibatch_05/checkpoint.pkl /tmp/ff3_rf_500"
# then transfer via base64 -- see the banner warning below -- and on midway2:
#   mkdir -p parameters/ff_3.0_trained_rf && cp sidechain.h5 environment.h5 there
#   echo 500 > parameters/ff_3.0_trained_rf/STEP
```

`extract_ff.py` and `compare_ff.py` are already uploaded to rockfish's `training/gly-sym/` and the
extraction was dry-run successfully there on 2026-09-08.

**Do not create `parameters/ff_3.0_trained_rf/` early with a mid-training force field.** The `STEP`
file exists to stop exactly that: the two runs are as far apart from each other as ~85 steps of
training moves them, so a step-354 rockfish file against a step-500 midway2 file would measure the
step gap and the trajectory gap together and separate neither.

**A binary file cannot be piped out of rockfish with `cat`.** The login banner lands on stdout, so
`ssh rockfish 'cat file' > local` prepends 1659 bytes and the HDF5 superblock is gone — the failure
reads as `file signature not found`, not as a truncated transfer. This is the same hazard already
recorded for `rsync` and `scp`, and it applies to `cat` too. Use a marker and base64, and verify:

```bash
ssh -o BatchMode=yes rockfish 'echo __B64__; base64 <file>' 2>/dev/null \
  | sed -n '/^__B64__$/,$p' | tail -n +2 | base64 -d > local_file
# then compare md5 against `ssh rockfish 'echo __M__; md5sum <file>'`
```

### Which trainer is the deliverable — decided 2026-09-08: both run, the arm test picks

midway2 is the **declared primary**: the seeds, the install chain, the production data and the
baseline all live there, so it wins any tie inside the 0.05 TM scatter. Both trainers run to 500 and
both force fields are measured, because the pair is worth more than the saved node-hours — it is the
only convergence evidence available, and rockfish's cost is charged to `rherna21`, 1,150,000
core-hours essentially unused.

The measurement that settled it (full detail in `findings.md`): at step ~355 the two runs had drifted
the **same distance** from ff_2.1 to within 0.1-1.7%, but in directions 23 degrees apart, leaving
them 33-40% of that drift away from each other. Training reproduces how far it moves and not where.
A single run therefore does not pin the force field to better than a third of what it changed, which
is why the pre-production check was restored rather than skipped.

### Outage lessons worth keeping (the outage itself is over)

Direct SSH from this Mac works again (`nc -z midway2.rcc.uchicago.edu 22` succeeds), so the midway3
`-L 2222` tunnel and `mdw2_via_tunnel.exp` are no longer needed; one push to
`scratchpad/mdw2_master.exp` is enough.

* **A returning login node is not a healthy filesystem.** Port 22 reopened while `/project` still
  stalled. Gate on a timed `ls` of the run directory *and* `parameters/`, or on a Slurm probe.
* **Probe a wedged filesystem with a Slurm job, never interactively.** A GPFS read in `D` state is
  uninterruptible, so `timeout 20 ls <wedged path>` neither returns nor dies: it hangs the shell and
  holds a ControlMaster channel until the mux refuses new sessions. Read the verdict from the probe
  job's exit code instead — `11` upside_input, `12` srun_mdw2.sh, `13` run_output, `14` the
  checkpoint, `COMPLETED 0:0` healthy, `TIMEOUT` still wedged:

```bash
sbatch --parsable -p broadwl -A pi-trsosnic -N1 -n1 -t 5 --job-name=fs-probe --output=/dev/null \
  --wrap="timeout 30 ls $D/upside_input >/dev/null || exit 11; ... exit 0"
sacct -j <jid> -o State,ExitCode -n
```

* **Do not let ssh fall through to password auth.** Wedged sessions hold mux channels, the session
  cap is hit, and ssh then tries a password — the RCC IP-ban trigger, and the ban took ~3 h to clear.
  Every cluster call must carry `-o BatchMode=yes -o PasswordAuthentication=no
  -o NumberOfPasswordPrompts=0`.
* **The damage was not uniform.** A probe from 00:25 CDT returned `run dir OK` and `parameters OK`
  while `ls ~` and the wedged minibatch directory hung, so specific inodes stall rather than the
  whole fileset going offline. Group quota just before the hang was 1.45 T / 1.49 T soft with no
  grace, so it was not a quota condition.
* **Connecting to a specific login node.** `scratchpad/rcc_master.exp <host> <socket> [host-key-alias]`
  opens a socket to any RCC node, reading the password from `~/.bin/ssh_mdw3`. All four
  `midway3-login[1-4]` present the same ED25519 key as `midway3.rcc.uchicago.edu`
  (`SHA256:DFRZlrKTqj6XjN78r5j/rEFqEKlz2yQpZGSuW9MPPr0`), so pass that as the alias rather than
  editing `known_hosts`. Without it the connection stops at a host-key prompt, the script answers
  *that* prompt, and **no Duo push is ever sent** — which looks exactly like a dead push.
### midway2: the original chain is re-armed (user, 2026-09-08 14:18)

The user cancelled the training-only link 48988387 and submitted **48994406 `ff-check-cont`**
(`check_continue.sbatch`, `--dependency=afterany:48988330`). There is no crontab on the account; this
was a deliberate manual action.

Consequence, which differs from a training-only link: when 48988330 ends, 48994406 continues to
`MAX_STEPS=500`, then extracts the force field to `parameters/ff_3.0_trained`, backs up the seeds,
writes `.ff_installed`, and chains into `decide_and_launch.sbatch`, which patches the four real seeds
with ARM_B and **relaunches glpG production on midway2 automatically**. It also runs
`push_progress.sh`, which commits and pushes to git. Both sentinels were absent when checked, so it
will fire cleanly.

* It counts checkpoint **files** (367) while the true step is **369** — an offset of 2 from
  checkpoints an old requeue deleted — so it trains ~2 steps past 500. Harmless.
* It does **not** touch NP, by design. The 500 A NP rebuild stays manual.
* It knows nothing about rockfish. Rockfish's glpG launch stays manual, via
  `scratchpad/ff3_retraining/install_and_launch_rf.sh`.
* `install_and_launch_mdw2.sh` is now a fallback for the case where the chain aborts.

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

### Disk: 246 GB reclaimed 2026-09-08 for the ff3.0 production campaign

Deleted with the user's approval: the six `NP-1AO6/prod/np.run.[0-5].up` replicas, ~41 GB each.
They were superseded on three independent counts — built on ff_2.1, run in a **200 A box** whose
PBC-honest protein extent is only ~188 A, and structurally unfolded to Rg 230.9 A, i.e. already
self-interacting through the boundary. `NP-1AO6` went from ~251 GB to **4.9 GB**; `prod/` keeps
`block_count`, `k190_results`, `logs_predate_fixes` and the job logs (953 MB).

Headroom was the reason: 1.45 T of a 1.49 T soft / 1.64 T hard group quota left only ~195 GB, and
4 relaunched glpG variants plus 6 NP runs at a larger box do not fit in that. Exceeding the hard
limit gives ENOSPC, which can corrupt an HDF5 file mid-write. **The quota number lags a large delete
by minutes — verify with `du`, not `rcchelp quota`.**

Small `.bak_*` files were deliberately left: `build_np_ff3.py.bak_pre_rama3`,
`np_hybrid.py.bak_pre_ffparam`, `footprint.sbatch.bak_hardcoded` and
`build_np_ff3.py.bak_pre_boxlen_*`. They are 32 KB each and are the record of which defect each fix
addressed.

### rockfish (JHU ARCH) — ff3.0 training, the outage-proof host

Set up 2026-09-07/08 because midway2 is unusable (GPFS down) and the Mac alone would not finish
until Thursday evening. **Continues the half-trained force field; it does not restart training.**

| item | value |
|---|---|
| account | `rherna21` — 1,150,000 core-hours allocated, **0.0 used** this quarter (ends 2026-09-30); this job needs ~2,400 |
| partition | `parallel`, **3-day wall**, so ONE job covers every remaining step — no chain, no stall guard, no install link |
| nodes | 2 x 48-core Xeon Gold 6248R @3.00 GHz (Cascade Lake) = the 96 cores a minibatch needs |
| repo | `/scratch4/rherna21/ywang268/upside2-md-rf` |
| build | `module load gcc/9.3.0 cmake/3.27.7 eigen/3.4.0 hdf5/1.10.7`, then `cmake ../src/ -DEIGEN3_INCLUDE_DIR=/data/apps/extern/spack_on/gcc/9.3.0/eigen/3.4.0-jof5yd7ppm3sdsc2dfd5t2idrzcwsz7b/include/eigen3 && make -j 12` |
| python | `module load python/3.11.9`, repo `.venv` with numpy, scipy, tables, h5py, torch (CPU wheel) |
| data | `training/gly-sym/upside_input` — 1827 files, 259 MB, transferred and counted |

**Why 2 nodes and not more:** `minibatch_size = 12` and `UPSIDE_TRAIN_NTHREADS = 8`, and
`n_threads` sets **both** the OpenMP thread count and the REMD replica count (`ConDiv.py:54,397`),
so it cannot be lowered without truncating the reference 8-replica ladder. 12 x 8 = 96 cores is the
algorithm's ceiling; extra nodes would idle.

**Rama library verified before anything ran.** `upside_input/rama.dat` md5 `932649af145b366d12791abd64915fd5`,
identical to `parameters/common/rama3.dat` (the GLY-symmetric library) and NOT the old asymmetric
`rama.dat` (`996a607b…`). Training on the wrong one reproduces the alphaL-biased GLY failure.

**Queue estimate at setup time:** `sbatch --test-only` for 2 nodes x 30 h answered
`Job 30725644 to start at 2026-09-08T12:23:38 ... on nodes c[132-133]`, i.e. ~11 h out
(cluster clock is EDT, so 11:23 CDT Tue). 1038 jobs were pending in `parallel`.

| JobID | Name | State | Notes |
|---|---|---|---|
| **30725720** | `upside-ff3` | **RUNNING** on `c[625,655]` since 2026-09-07 23:49 CDT, 48 h wall (ends 2026-09-09 ~23:00 CDT) | Resumed at **step 275**; newest written checkpoint `epoch_09_minibatch_11` = **step 354** as of 2026-09-08 15:00 CDT. 0 failures. `extract_ff.py` and `compare_ff.py` are uploaded here and the extraction is dry-run verified. Log `upside-ff3_30725720.out`, run dir `run_output_ff3/`. |
| 30726466 | `upside-ff3-cont` | PENDING (`afterany:30725720`) | Chain link: resumes the newest checkpoint if the wall kills the trainer, exits immediately if step >= 500. Self-perpetuating (queues its own successor before training starts) with a 3-stall abort guard. |

Retired 2026-09-08 09:12 CDT at the user's request: trainer B `30725855` + link `30726467`,
CANCELLED at step 326, run dir `run_output_ff3b/` left on disk. It existed to cover "one run goes
bad"; that role is now filled by midway2 being back.

**`continue_rf.sbatch` installs nothing.** rockfish holds only `train_rf.sbatch` and
`continue_rf.sbatch` — no `extract_ff.py`, no `decide_and_launch.sbatch`, and no `.ff_installed` or
`RESUME_STEP` flag files. At step 500 this host stops and waits.

**Measured rate: 516 s/minibatch** (first step, 0 failures, Median RMSD 0.86 / 2.32). That is 1.1x
midway2's 565 s and 2.2x the Mac's 1141 s — the 390 s projected from clock speed was optimistic.
225 steps x 516 s = **32.3 h, finishing Wed 2026-09-09 ~08:00 CDT**, inside the 48 h wall with 16 h
spare and a full day before Thursday ends.

**The 30 h wall was not enough, which is why 30725696 was replaced.** 226 x 516 s = 32.4 h, so that
job would have been wall-killed around step 483 — 17 steps short, after 30 h of compute. A running
job's `TimeLimit` cannot be raised by a user, so the fix was cancel and resubmit at 48 h (partition
max is 3 days). Only the in-flight minibatch was lost, because step 275's checkpoint was already on
disk. **Check the measured `seconds elapsed` against the requested wall after the first step of any
new run** — the estimate from clock speed was off by 32%.

**The queue estimate was wrong in our favour.** `--test-only` predicted a 2026-09-08T12:23 start;
the job actually started **within seconds** of submission. Do not plan around `--test-only`.

**Superseded: 30725696** (30 h wall, see above) and **30725686**, submitted 12 minutes earlier and cancelled. It ran on the *old* snapshot of
`train_rf.sbatch`, which took argument 2 as a literal step count, so its log read `for 500 steps` —
it would have trained 274 → 774, overshot the MAX_STEPS=500 target by 274 steps and been wall-killed
at 30 h. Slurm snapshots the batch script at submission, so editing the file could not fix the
running job; it had to be cancelled and resubmitted. The corrected script derives the count at run
time from `RESUME_STEP` and logs `resuming at step 274 for 226 steps -> 500`. **Check that line in
the log after every submission** — it is the only place the arithmetic is visible.

Also cancelled: 30725687, a 1-step validation job on `express`. It was pending with
reason `PartitionConfig` (express refuses 2-node requests) and became pointless once the real job
started immediately.

**The transfer preserved the half-trained force field exactly, and this was verified, not assumed.**
`extract_ff.py` run on both the Mac's source checkpoint and the adopted Rockfish checkpoint gives
identical SHA-256 digests on all four trained arrays:

| array | sha256 (first 24) | sum |
|---|---|---|
| `pair_interaction` | `ef24b0b911287425efd289ca` | +1.6294386169e+04 |
| `coverage_interaction` | `191d6fbddeba8837d4f86c81` | +7.4490769816e+03 |
| `hydrophobe_interaction` | `58d4010f74b611279ca4bc1e` | +8.9675347942e+03 |
| `environment/energies` | `5ed10d16729d92c698707cf7` | +6.2365761830e+01 |

Adoption used `scratchpad/ff3_retraining/rf_adopt_checkpoint.py`, a straight path rewrite of the
Mac's `checkpoint.pkl` (md5 verified end to end) — **no `pack_param` refit**, so the latent vector is
the Mac's exactly and `solver step_num 5` carried the Adam moments over. Both sides run numpy 2, so
the pickle transfers directly; the `.npz` export path is only needed for a numpy 1.x destination.

**The submit script derives the step count at run time**, from a `RESUME_STEP` file next to a
hand-built checkpoint, or from the directory name (`epoch_XX_minibatch_YY` -> `XX*38 + YY + 1`) when
resuming from inside a run. So a resubmission always runs exactly the steps still missing, and the
log states its own arithmetic (`resuming at step 275 for 225 steps -> 500`). Read that line.

**Transfers to Rockfish must not use `rsync` or `scp`.** The login banner makes the shell unclean:
`rsync` dies with `protocol incompatibility` and `scp` with `Received message too long`. Use
`tar czf - … | ssh rockfish "tar xzf - -C …"` or `cat file | ssh rockfish "cat > dest"`, and strip
the banner from command output with a self-emitted marker (see §0).

**Do not use the system `python3` on the login node** for checkpoint surgery — it cannot read pickle
protocol 5 (`unsupported pickle protocol: 5`). Activate `$RF/.venv` first.

### Checkpoint bookkeeping and cross-host transfer (lessons; no local run is active)

The Mac trainer that carried training through the outage was stopped 2026-09-08 09:12 CDT at step
304, after 35 minibatches and 0 failures. Its step-275 checkpoint is what seeded rockfish. What is
worth keeping from it:

* **Count written `checkpoint.pkl` files, never minibatch directories.** `main_loop` creates
  `epoch_XX_minibatch_YY/` when a step *starts* and writes `checkpoint.pkl` when it *ends*, so
  globbing directories reports a step that has not happened. The step number of a written checkpoint
  is `XX*38 + YY + 1`.
* **Worker processes do not die with the driver.** `kill <driver>` leaves the 12 `ConDiv.py worker`
  subprocesses and their `upside` children running on every core. Stop them in order: `kill
  <driver>`, then `pkill -f 'ConDiv.py worker'`, then `pkill -x upside`, and confirm `pgrep -cf
  ConDiv.py` and `pgrep -cx upside` both read 0.
* **A checkpoint moves between hosts by path rewrite, not by refit.**
  `rf_adopt_checkpoint.py` rewrote the paths inside the Mac's `checkpoint.pkl` with **no `pack_param`
  refit**, so the latent vector was preserved exactly and `solver step_num` carried the Adam moments
  over. Verified rather than assumed: `extract_ff.py` on the source and the adopted checkpoint gave
  identical SHA-256 digests on all four trained arrays (`pair_interaction`
  `ef24b0b911287425efd289ca`, `coverage_interaction` `191d6fbddeba8837d4f86c81`,
  `hydrophobe_interaction` `58d4010f74b611279ca4bc1e`, `environment/energies`
  `5ed10d16729d92c698707cf7`). Both sides run numpy 2, so the pickle transfers directly; the `.npz`
  export detour is only needed for a numpy 1.x destination, whose array pickles cannot reference
  `numpy._core`.

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


### Disk: the GPFS group quota is the binding limit, NOT `df`

`df -h /project` reports over a terabyte free and is **misleading**. The real constraint is the
per-group GPFS quota, visible only via `rcchelp quota`:

```
trsosnic   blocks (group)   used 1.45T   soft 1.49T   hard 1.64T
           files  (group)   used 596023  soft 728600  hard 801460
```

Exceeding the hard limit fails writes with ENOSPC, which can corrupt an HDF5 file mid-write — the
worst possible failure for an unattended run. Check `rcchelp quota` before anything that writes tens
of GB. Note the quota accounting updates on a timer, so it lags a large delete by minutes; verify
with `du` instead of waiting for the number to move.

**Footprints measured 2026-09-07:** training `run_output` ~17 MB per minibatch (~10 GB for a full
600-step run); a glpG seed is 191 MB, so a 28-replica arm costs 5.2 GB before trajectory growth and
the two-arm test ~14 GB all-in; each glpG variant's live replicas run 22.5 GB and the whole variant
directory reaches 43-53 GB once trajectories accumulate. `decide_and_launch` deletes the live
replicas before resubmitting, so the relaunch is net **-69 GB**, not a cost.

**2026-09-07 cleanup, 69 GB reclaimed.** Deleted three sets of replica trajectory backups whose data
this file already rules out — `.bak_alphaL_run` (112 files, 27.5 G, alphaL era), `.bak_broken_gly_2`
(112 files, 25.8 G, pre-GLY-fix) and NP `.bak_broken_gly` (6 files, 14.9 G). All seed-level backups
were kept (`.bak_rigid_stage`, `.bak_broken_gly` seeds, `.bak_alphaL_restart`, 2.25 GB total) since
they are the cheap record of seed state. Take care with the patterns: glpG has a *seed* set named
`.bak_broken_gly` that must survive while NP's *replica* set of the same name is the 14.9 GB target.

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

### The unattended chain (built 2026-09-05) — DISARMED as of 2026-09-08

**Nothing below is queued any more.** The GPFS outage broke the chain: `ff-check-cont` (48981236)
was cancelled with its trainer, and the links that replaced it on both hosts
(`continue_mdw2.sbatch`, `continue_rf.sbatch`) only resume training and explicitly refuse to install.
Every script named here still exists in `training/gly-sym/` and is still the intended design; it has
to be re-armed by hand, and `decide_and_launch.sbatch` needs re-reading first because glpG
production is now COMPLETED rather than running (see §1). The description is kept as the design of
record.

The original rationale, still valid: the decision logic belongs in Slurm scripts rather than in
monitoring, because a Claude session exists only while the Mac is on.

```
upside-gly-sym --afterany--> ff-check-cont  (loops to 600 minibatches)
                                  | install branch: extract_ff.py -> parameters/ff_3.0_trained/,
                                  |                 back up 4 seeds, write .ff_installed
                                  v
                            armtest-build (run_arm_test.sbatch)
                                  |-> armtest_A_<V>  (28 cpu, 12 h)  --+
                                  |-> armtest_B_<V>  (28 cpu, 12 h)  --+
                                  +-> arm-decide  --dependency=afterany:A:B  <-+
                            arm-decide (decide_and_launch.sbatch)
                               ARM_A / ARM_B -> patch the 4 real seeds, verify, scancel mdw2_glpG*,
                                                wipe replicas, block_count=0, submit_remd.sh x4
                               NO_WINNER     -> exit 3, cancel nothing, launch nothing
                               always        -> push_progress.sh (EXIT trap)
```

`decide_and_launch` is submitted by `run_arm_test`, not by `ff-check-cont`, because the two arm job
ids do not exist until `run_arm_test` runs. If `run_arm_test` dies before queueing it, its `die()`
writes the reason to `arm_verdict.txt` and runs `push_progress.sh`, so the failure still reaches git.

**Arms.** A = plain hybrid + trained `pair_interaction` only. B = `martini_inject_coverage.py` adds
`hbond_coverage` + `hbond_coverage_hydrophobe` (+ prerequisite `placement_fixed_point_vector_scalar`),
then the full trained FF. Measured on the 79HIS seed: arm A energy −25107.28 (25 nodes), arm B
−24995.55 (28 nodes) — the two coverage terms contribute **+111.73 E_up**. Full `/input` diff vs the
pristine seed: arm A changes exactly one dataset; arm B changes that one and adds exactly 17.

**Assertions before anything launches** (and again on the production seeds before any `scancel`):
`exclude_intra_protein_martini == 1`, `current_stage == production`, and the three hybrid interface
nodes present. Proven to fire: a tampered copy with the flag cleared exits 1 — and its energy jumps
to **+51070** vs −24996, which is what dry-MARTINI acting inside the protein looks like.

**Two toolchain traps, already handled.** The venv has pytables but *no h5py*; the module python has
h5py but *no pytables*. Both sbatch scripts capture `PY_H5` before activating the venv and `PY`
after, and assert both imports up front. Do not pip-install into the venv while training is live.

**`decide_arm.py` frame selection.** It skips `output_previous_0` — that group is the seed's own
pre-test `/output`, rotated there by the first `reseed()`, and is old-force-field data (~17% of
frames). It also orders groups numerically, since a string sort puts `output_previous_10` before
`output_previous_9`. Verified on real data: exit 2 / `NO_WINNER` with no trajectories, exit 0 with a
winner otherwise; ~14 s for 12 replicas, so minutes for 56.

**NP is now IN the chain** (enabled 2026-09-05 at the user's request). It runs after the glpG
relaunch, so a `NO_WINNER` verdict (exit 3 at line 93) leaves both systems untouched. NP already
carries both coverage nodes, so all three trained tables install directly and
`martini_inject_coverage.py` is never run on it — the arm decision does not apply to NP. Its
environment is `sigmoid_coupling_environment` while the trained environment is nonlinear, so **NP
keeps its old environment table**; that is a known limitation, not an oversight. Patch and verify
happen before the `scancel`, so a failure cannot leave NP cancelled with nothing resubmitted.
Verified on a live replica: `check_hybrid_up.py --require "$NP_NODES"` passes (8608 atoms, 36 nodes,
energy −11159.36, interface intact, intra-protein MARTINI excluded, production stage).

**Legacy-ff jobs stopped 2026-09-05 15:03.** All four glpG variants and NP were running on the old
force field, so their output would have been discarded at install time. `STOP` files were written
into each run dir (`popepopg_REMD_mdw2/<V>/STOP` and `NP-1AO6/prod/STOP`), which `run_remd.py:185`
and `run_np_prod.py:131` poll between chunks — each job finishes its current chunk and exits without
resubmitting, rather than being killed mid-write. `scancel` was not used because the midway2 login
node was refusing connections; the STOP path is the cleaner stop anyway. `decide_and_launch.sbatch`
removes every `STOP` before resubmitting (glpG at line 150, NP in its own block), so the relaunch is
not blocked by them. **If you ever restart these by hand, delete the STOP file first or the driver
will exit immediately.**

**Monitoring without midway2.** `/project` is shared between midway2 and midway3, so training
progress, logs and run dirs are all readable from the midway3 socket when the midway2 login node is
down. Only `squeue`/`sbatch` need midway2. The chain is unaffected by a login-node outage: every
`sbatch` in it is issued from inside a running job on a compute node.

**Backups before anything destructive:** seeds and NP replicas copied to `*.bak_pre_ff3_<stamp>`
with `cp -n`. Reverting = copy them back.

**GLY symmetry: there are TWO rama libraries, and picking the wrong one is the whole hazard.**

| file | GLY asymmetry | belongs to |
|---|---|---|
| `parameters/common/rama.dat` | coil 6.31, sheet 67.47 E_up | the **old** force field — asymmetric *by design*, leave it alone |
| `parameters/common/rama3.dat` | **0.0000 / 0.0000** | the **GLY-symmetric** force field — use this for anything trained |
| `training/gly-sym/upside_input/rama.dat` | 0.0000 / 0.0000 | training's copy, byte-identical to `rama3.dat` |

`rama3.dat` and the training copy have the same md5 (`932649af…`); `rama.dat` is `996a607b…`.

**Every existing config is already correct.** glpG seeds, glpG live replicas and NP live replicas all
measure worst GLY `sym_err = 0.000000` (chiral controls 11.2 E_up, so the metric does detect real
chirality). The arm test copies existing seeds, so it is unaffected.

**The real defect was in the NP builder, not the library.** `np_hybrid.upside_ff_paths()` hardcodes
`common/rama.dat`, so a trained-ff build would have silently produced αL-biased GLY maps — the
documented TM4 failure mode. Fixed 2026-09-07: `build_np_ff3.py` takes `--rama-library`, defaulting
to `rama3.dat`, and `np_hybrid.build_system(spec, work_dir, ff=None)` now accepts an ff override
(default unchanged, so `build_k190*.py` still work).

Note the trap that override exposed: overriding the dict in the *caller* was not enough, because
`build_system` re-derived the paths internally. The symmetry check would have printed
"GLY rama asymmetry 0" while the build used the asymmetric file. **Verify that a force-field
override reaches the writer, not just the pre-flight check.**

**Do NOT "fix" `parameters/common/rama.dat`.** It is deliberately the old force field's library.
Overwriting it with the symmetric maps (attempted and reverted 2026-09-07) silently changes the old
force field for every consumer — `py/martini_prepare_system.py` and all of
`example/08.MembraneSimulation/` read it.

**Two grid conventions, and they differ — this cost a false alarm.** The GLY mirror
`(phi,psi) -> (-phi,-psi)` is:
* `m[::-1, ::-1]` for the `.up` `rama_map_pot/rama_pot` maps (bin centres at `-180+(i+0.5)*5`);
* `roll(m[::-1, ::-1], 1)` for the library's `dimer_pot` (bin centres at `-180+i*5`).

Using the library convention on a `.up` map fabricates ~3.3 E_up of asymmetry on perfectly good
seeds. Always validate the metric against a chiral control (ALA/SER/HIS must show ~10-11 E_up) and
against `dG(aR->aL)`, which is the physically meaningful quantity and reads exactly 0.000 on a
correctly symmetrized map.

**CWD fix:** `srun_mdw2.sh` now `cd "$PROJECT_ROOT"` before running ConDiv; the checkpoint stores relative paths from the project root.

**Init (2026-09-04 11:13):** pack_param converged to loss=54.1821 (palindrome floor). 38 minibatches × 12 proteins = 456 proteins. GLY dp1 forced palindromic by `rotamer_parameter_estimation.py`.

**Seven bugs fixed 2026-09-04:**
1. `--slurmd-debug=0` in the srun worker launch is restricted to root/SlurmUser on midway2 Slurm. Removed from ConDiv.py srun args.
2. `shutil.copy` does not set the execute bit on the copied `run_output/ConDiv.py`. Workers launched by srun with just `worker_path` got `execve() Permission denied`. Fixed: pass `sys.executable, worker_path`.
3. `compute_divergence` referenced `sigmoid_coupling_environment` (wrong intermediate fix) → reverted to `nonlinear_coupling_environment`. The env spline type must be 0 (nonlinear) for training, not 1 (sigmoid). Added `environment_potential_type=0` to kwargs in ConDiv.py worker call.
4. `py/run_upside.py` checked `if environment_potential_type:` which is `False` for integer `0`. `--environment-potential-type=0` was never passed to upside_config.py, so all base.h5 files got the sigmoid node instead of nonlinear. Fixed: `if environment_potential_type is not None:`. Applied both locally and on the cluster.
5. `compute_divergence` passed `coeff.shape = (20, 18)` → product 360 to `get_param_deriv`, but the C++ `NonlinearCoupling::get_param_deriv()` returns `coeff_grad (360) + weights_grad (400) = 760` elements concatenated. Fixed: read both `coeff.shape` and `weights.shape[0]`, call `get_param_deriv((760,), ...)`, then slice `full_grad[:360].reshape(20, 18)` to recover just the coeff gradient.
6. `d_obj` in `backprop_deriv`: `rot_expect`, `hb_expect`, and `hyd_expect` were 1D tensors of shape `(n_knot_sc-2,)` / `(n_knot_hb-2,)`. PyTorch broadcasts 1D → last dim, but `rot_e`/`cov_e`/`hyd_e` have shape `(..., n_dist-2, n_ang-2, n_ang-2)`, so the subtraction `rot_e - rot_expect` tried to match dim 4 (13, angular) against the expectation (10, distance). Fixed: `.view(-1, 1, 1)` on all three expectations so they align with the distance dimension (dim 2). Applied locally and to remote `run_output/ConDiv.py`.
7. Same fix as bug 6 but was NOT applied to the remote *main* ConDiv.py at `training/gly-sym/ConDiv.py`. The worker copy (`run_output/ConDiv.py`) only runs `compute_divergence`; `backprop_deriv`/`d_obj` run in the main process from `training/gly-sym/ConDiv.py`. Patched the correct file.

**Submit rule:** always submit gly-sym (and all other training/simulation jobs) through the **midway2 SSH socket**. Submitting via the midway3 socket routes to caslake even when the sbatch script requests broadwl.

**Resubmit:** `cd ~/project/yinhan/upside2-md-mdw2/training/gly-sym && sbatch srun_mdw2.sh training/gly-sym/run_output/<latest_checkpoint.pkl> 200`

Note: checkpoint path must be relative to PROJECT_ROOT (the repo root), not to the gly-sym directory, because srun_mdw2.sh does `cd "$PROJECT_ROOT"` before running ConDiv.

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
| method | regular MD, 6 independent trajectories, single T=0.8647, no exchange | **REMD**, 48 replicas, T ladder 0.70–0.90, configuration exchange |
| purpose | nanoparticle adsorption footprinting (K190 exposure) | **HDX** protection factors / ΔG |
| system | 1AO6 albumin 578 res + 5 nm MPA-AuNP, 8608 atoms, box 300 Å | glpG 210 res in a DDM micelle, 3156 atoms, box 137 Å |
| composition | PROTEIN 2890 + GOLD 887 + MPA 203 + ION 4628 (K+ 2423 / Cl- 2205, 0.15 M KCl) | PROTEIN 1050 + LIPID 1674 + ION 432 |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: 2736 atoms (ION+LIPID+630 protein) overdamped **Brownian**; 420 protein atoms velocity-Verlet |
| timestep | **0.001**, freely settable at runtime | **0.009, HARD-LOCKED** by `/input/brownian/numerical_time_step`; `martini_brownian.cpp:100` throws on mismatch. Friction is tuned against it for lipid D=11.5 µm²/s — **do not change it** |
| detection | non-finite positions OR ≥5 stretched bonds | non-finite potential (whole chunk) OR ≥5 stretched bonds |

---

## 3. NP campaign — `np_1AO6_prod`

**Unfolding is the expected result, not a failure.** 1AO6 albumin spreads on the MPA-AuNP surface. Rising Rg (currently up to 230.9 Å on run.3, block 3) is the intended observable and must **not** be reported as a blow-up. Judge health on non-finite frames, peptide C–N bonds, and `avg_kinetic_energy/1.5kT`. (Contrast glpG, where Rg ~19 Å is the health signal.)

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

**Dir** `~/project/glpG_DDM_micelle_REMD/` — one subdir per variant, each with 48 `*.run.N.up`,
`remd.<jobid>.out`, `block_count`.
**NOTE the directory name**: `glpG_DDM_micelle_REMD`. The older `glpG_DDM_REMD` (lamellar, 72 G) is gone.
**Driver** `run_remd.py` · **sbatch** `remd.sbatch` · **submit** `submit_remd.sh <variant>`
**Variants:** `glpG-RKRK-79HIS`, `glpG-RKRK-79HIS_S115T`, `glpG-RKRK-79ALA`, `glpG-RKRK-79ALA_S115T`

Config: 48 replicas, T 0.70–0.90, `REMD_DT=0.009`, `--replica-interval 0.09`, `--exchange-criterion 0`,
swap sets A=(0-1,2-3,…) B=(1-2,3-4,…), 300 frames/chunk, `REMD_MAX_BLOCKS=12`.

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
  d=~/project/glpG_DDM_micelle_REMD/$v; f=$(ls -t $d/remd.*.out | head -1)
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
`~/project/glpG_DDM_micelle_REMD/run_remd.py` and `~/project/NP-1AO6/run_np_prod.py`.
No version history exists for them. Edit directly on the cluster.
A running job keeps the version it loaded at start; edits take effect at the **next block**.

---

## 7. If the chain terminates — manual rollback procedure

**glpG (old driver, or gate fired before rollback logic):** patch last output frame of each NaN file
with the last finite frame from `output_previous_0` (end of block 1), then resubmit.

```python
import h5py, numpy as np
from pathlib import Path
run_dir = Path("~/project/glpG_DDM_micelle_REMD/<variant>").expanduser()
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
# then: bash ~/project/glpG_DDM_micelle_REMD/submit_remd.sh <variant>
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
- **glpG-DDM micelle campaign (closed 2026-08-13).** All four variants failed at block 2–3; used a constant `--seed` causing rollback to re-run the identical failing chunk deterministically. MBAR fix (findings 91) also required before that data was usable. HDX ΔG results in `~/Downloads/glpG_DDM_micelle_HDX_dG/`.
