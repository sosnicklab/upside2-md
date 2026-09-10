# Remote jobs on midway2/midway3 — status and handbook

Snapshot: **2026-09-10 ~00:00 CDT. Two campaigns are live. (1) The glpG two-arm force-field test runs on BOTH clusters as independent replicates: midway2 `49002097`/`49002098` with decider `49002099`, rockfish `30768485`/`30768486` with verdict-only decider `30768512`. Verdicts ~05:20 (rockfish) and ~05:50 (midway2) CDT 2026-09-10. (2) The **ff3.0 re-benchmark of Peng et al. JCTC 2022** launched on midway2 broadwl: 32 jobs, `49002902`-`49002933`, 16 proteins x native/de novo, 14-replica REMD each, self-resubmitting until their Table S2 durations complete. Upside now has ONE shared deployment at `/beagle3/trsosnic/yinhan/upside2-md` usable from both clusters via `env_shared.sh`; see §0b. ff3.0 is deployed in the repo as `parameters/ff_3.0`. midway3 has no other work.**
Written so a fresh session can pick up cold. Everything needed to connect, check health correctly,
and react to a failure is here. Job state below is live; superseded jobs are not listed, only
summarised in §8 where they carry a lesson.

---

## 0. Connect first (needs a Duo push on the user's phone)

Key-based auth is NOT enabled; password + Duo is the only method. The ControlMaster socket expires
roughly hourly, so expect to redo this most sessions.

**midway2** (POPE/POPG REMD campaign): **this Mac's IP is BLOCKED on midway2's login nodes.**
`mdw2_master.exp` fails with `Connection refused` on port 22 while midway3 answers normally, so the
refusal is an address block, not an outage and not a Duo problem. Reach midway2 through midway3 in
two steps (confirmed working 2026-09-09; connect midway3 first, then):
```bash
ssh -S ~/.ssh/cm-mdw2.sock -O check yinhanw@midway2.rcc.uchicago.edu   # alive?
# 1. port-forward midway2:22 to localhost:2222 over the existing midway3 master (no Duo)
ssh -o BatchMode=yes -f -N -L 2222:128.135.112.69:22 \
    -S ~/.ssh/cm-mdw3.sock yinhanw@midway3.rcc.uchicago.edu
# 2. open the midway2 master over that forward           # USER MUST APPROVE DUO
expect /Users/yinhan/Documents/upside2-md/scratchpad/mdw2_via_tunnel.exp
ssh -S ~/.ssh/cm-mdw2.sock yinhanw@midway2.rcc.uchicago.edu '<command>'
```
**Use `mdw2_via_tunnel.exp`, not `mdw2_via_mdw3.exp`.** The ProxyCommand variant authenticates and
prints "ssh backgrounded", then the master dies without ever binding its socket, it looks like a
success and leaves you with no socket. The script's own header records this; I ignored it and burned
a Duo push. The forward must be opened *before* the expect script runs.

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

Snapshot **2026-09-10 ~09:20 CDT (verified live against `squeue`)**.

### Campaign 1: glpG production on ff3.0, lowered ladder ceiling (both clusters, 2026-09-10 ~09:10)

| host | JobIDs | variants | wall |
|---|---|---|---|
| midway2 | 49003548-49003551 | all four | 36 h x 4 blocks |
| rockfish | 30775667-30775670 | all four | 48 h |

**The arm test is decided: R wins**, so `parameters/ff_3.0` now carries the **rockfish** training
run's tables (`sidechain.h5` `cfe4ba5e...`, `environment.h5` `301f8418...`, both step 500),
byte-identical in the repo, on beagle3 and on rockfish. The arm-M tables are kept beside them as
`*.bak_armM_20260910-*`. Scored on the segments the build actually treats as helix rather than
`decide_arm.py`'s capped window, over replicas 0-2:

| measure | arm M | arm R | R - M |
|---|---|---|---|
| TM4a 134-139 | 0.979 | 0.987 | +0.008 |
| TM4b 141-155 | 0.884 | 0.946 | **+0.062** |
| TM1 30-43 | 0.873 | 0.944 | **+0.071** |
| TM3 99-115 | 0.931 | 0.931 | 0.000 |

**Ladder ceiling lowered 0.90 -> 0.82, from measurement.** Pooling all 8428 frames of the arm-R run
by their logged instantaneous temperature, TM4's 134-139 segment holds 0.981 below T = 0.72 and
collapses to **0.201** in the 0.88-0.90 bin, while TM3 and TM6 barely move. Whole-fold helix content
sits on a plateau at 87% of the coldest bin from 0.78 to 0.84, then falls to 72% at 0.90. 0.82 is the
top of that plateau. `REMD_N` stays 28 on both clusters, so the narrower span also tightens the
spacing (0.0043 at the cold end, was 0.0071) and raises exchange acceptance.

* Set through `REMD_T_HI=0.82` in each cluster's submit script; `run_remd.py` builds
  `linspace(sqrt(T_LO), sqrt(T_HI), N)**2` and needed no change. Verified in the running job:
  ladder starts 0.7, 0.704275 and ends 0.815386, 0.82 over 28 replicas.
* **Seeds were pristine ff_2.1 with no coverage nodes**; the cancelled decider never patched them.
  Both clusters were patched by hand with `martini_inject_coverage.py` then `patch_seeds.py` against
  the arm-R `sidechain.h5`, backed up first as `*.bak_pre_armR_20260910-*`. All four seed energies
  agree exactly across clusters (-24957.68 / -24843.17 / -24903.80 / -24806.57), so the two runs are
  true replicates. `rockfish` had no `patch_seeds.py`; it was copied from midway2.
* **No environment nodes were added.** The arm test that chose R ran without them, and adding them
  locally traded TM4a (-0.116) against TM4b (+0.063) with no net gain, so the ladder ceiling is the
  only variable that changed.
* midway2's 130 GB of pre-ff3 replicas were moved to `<variant>/pre_ff3/` (same filesystem, nothing
  deleted); `hdx/` was left alone. rockfish had no prior production data.

### Campaign 2: ff3.0 re-benchmark of Peng et al. JCTC 2022 (launched 2026-09-09 ~23:55 CDT)

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

### Campaign 4: NP production on ff3.0 (midway2, launched 2026-09-10 04:31 CDT)

| JobID | what | log |
|---|---|---|
| 49003318 | `np_1AO6_prod`, block 1 of 8, self-resubmitting | `NP-1AO6/prod_ff3/np.49003318.out` |

Build `49003158` finished rc=0 and all six configs passed `verify_np_ff3.py` (composition, ion
counts, 300 A box, `exclude_intra_protein_martini` 1, environment nodes present, no
`sigmoid_coupling_environment`, coverage nodes wired into the rotamer, finite energy and gradient).

* **dt is pinned to 0.001 in `np_prod.sbatch`.** `run_np_prod.py` defaults `NP_DT` to 0.009, which
  is wrong for NP: the short step is required for accuracy while the backbone unfolds.
* `np_prod.sbatch` and `submit_np.sh` now point at `prod_ff3/` and source the shared beagle3
  `env_shared.sh` instead of the old per-cluster `upside2-md-mdw2`. The `martini_upgrade_hybrid_args`
  call was dropped: it exists to migrate old replicas, and these configs are freshly built.
* The old `prod/` replicas (246 GB, 98 output groups, old-force-field coordinates) are untouched and
  can be deleted once `prod_ff3` has produced a block.

### midway3, IDLE, nothing queued or running

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


### Disk: check BOTH the group quota and the fileset, either can bind

**Corrected 2026-09-09.** This section previously said the group quota is the binding limit and
`df` is misleading, and quoted `used 1.45T / hard 1.64T`. **Those are `/project2`'s numbers, not
`/project`'s**, a different filesystem, and not where any of this data lives. `rcchelp quota`
reports four separate `trsosnic` group rows and the one that governs `/project/trsosnic` is the
`Midway3 GPFS mounted at /project` section (see §1 for the full table).

The two limits for `/project/trsosnic`, measured 2026-09-09:

```
group quota (rcchelp, "mounted at /project" section):  used 2.36T  soft 3.49T  hard 3.84T -> 1516 G free
fileset (statvfs / df on /project/trsosnic):           3929 G total, 2.4T used            -> 1514 G free
```

They agree to within 2 G, so neither is misleading here and **`df` on the fileset is a sound
number**, what was misleading was reading the wrong quota row. The safe check is `min()` of the
two, since a fileset smaller than the group limit inverts which one binds. `check_quota.py` now
does exactly this; use it rather than reading `rcchelp` by eye.

Do **not** run `df` or `statvfs` on `/project` itself: that is the whole 6.3 PB `midway3_cap`
device and reports 1.7 PB free. Use `/project/trsosnic`.

Exceeding the hard limit fails writes with ENOSPC, which can corrupt an HDF5 file mid-write, the
worst possible failure for an unattended run. Note the quota accounting updates on a timer, so it
lags a large delete by minutes; verify with `du` instead of waiting for the number to move.

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

**Dir** `~/project/yinhan/popepopg_REMD_mdw2/` on **midway2** — one subdir per variant, each with 28
`*.run.N.up`, `remd.<jobid>.out`, `block_count`.
**Driver** `run_remd.py` · **sbatch** `remd.sbatch` · **submit** `submit_remd.sh <variant>`
**Variants:** `glpG-RKRK-79HIS`, `glpG-RKRK-79HIS_S115T`, `glpG-RKRK-79ALA`, `glpG-RKRK-79ALA_S115T`

Config: 28 replicas, T 0.70–0.90, `REMD_DT=0.009` (hard-locked), `REMD_MAX_BLOCKS=4` for the ff3.0
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
- **glpG detergent campaign (closed 2026-08-13, model retired 2026-09-09).** Kept for one lesson only: it used a constant `--seed`, so a rollback re-ran the identical failing chunk deterministically, which is why the driver now takes a per-chunk seed. All four variants failed at block 2–3. Its HDX ΔG output is in `~/Downloads/glpG_DDM_micelle_HDX_dG/` and is no longer cited.
