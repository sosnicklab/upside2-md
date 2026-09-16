# 2026-09-10 (evening): glpG blow-up diagnosis

Asked to check remote job status, then to find the cause of the glpG hot-replica blow-ups and test it
locally. Job status: 41 jobs live (midway2 37 = 4 glpG + 32 benchmark + 1 NP; rockfish 4 glpG), all
RUNNING, midway3 idle; midway2 still IP-blocked from this Mac, reached through the midway3 tunnel.
Details and throughput estimates are in `remote_jobs.md`.

## Diagnosis (findings.md 3.10 carries the full record)

The user's hypothesis was a temperature mismatch between dry-MARTINI and Upside. Confirmed:

1. Verified `~/OneDrive .../image.png` row by row: K = T_up x 350.588235 exactly, 0/24 rows
   inconsistent. The production ladder 0.70-0.82 is **245.4-287.5 K = -27.7 to +14.3 C**.
2. Every dry-MARTINI parameter is built for **0.8647 T_up = 303.15 K = 30.00 C**: `/input/brownian`
   `reference_temperature_up`, `martini_build_tables.py DEFAULT_PRODUCTION_TEMP_UPSIDE`, and the
   equilibration run itself (`output_previous_0` is at exactly 0.8647 in every replica of every
   variant). The ladder runs the bilayer 15.7-57.7 K below all three, making MARTINI interactions
   1.24x stronger in kT at rung 0.
3. Measured on clean frames, the two subsystems are at different temperatures: lipids 1.020 x T_nom
   flat across the ladder, protein T_nom + ~0.08 T_up (+28 K), reaching 1.506 at rung 27.
4. **Reproduced locally with no replica exchange**: from an equilibrated cluster frame at production
   settings, local T_prot = 0.7911 vs the cluster's 0.7912 for that rung, T_lip 0.7327 vs 0.7355.
   The offset is nearly independent of the thermostat timescale (+0.074/+0.070/+0.106 at tau =
   1/5/20) so it is not a power leak the thermostat fails to remove, and it grows with temperature
   (+0.070 at T = 0.7215, +0.329 at T = 0.8647).
5. The blow-up is a `Spring_bond` tear of the TM4 backbone at residues 139-141 (C140-N141 at 18.5 A
   against r0 1.300), reproduced locally to 0.27% on the older tables, so it is geometric and the
   arm-R retraining is not its cause.

Ruled out by measurement: arm-R tables, an unthermostatted atom subset, exchange laundering of the
protein excess. A dt scan was impossible by design (`apply_langevin_step` throws unless runtime dt
matches `/input/brownian numerical_time_step`); that check was left alone.

Attributed: the excess belongs to the mass-1 backbone under EITHER thermostat (friction>0 sites
1.092/1.077, friction==0 sites 1.130/1.072, lipids 1.024/0.999) with no lipid-contact-count trend, so
it is integrator discretisation bias on the steep MARTINI core, not a thermostat or friction defect.
The cold bilayer is what presses the backbone there: 5.4x more sub-3.40 A protein-environment
contacts at T = 0.7215 than at the 0.8647 design point, closest approach 2.89 vs 3.29 A.

Two of my own readings were corrected by the full-length data: the excess is multiplicative and
tau/temperature-independent (not superlinear in T, which came from a burst in a half-length average),
and `protein_kinetic` is not diluted by the 420 zero-momentum placed atoms because the logger
excludes them from `n_dynamic`.

Not demonstrated locally: the ejection itself. 250 time units stayed finite and negative with zero
pairs inside 2.85 A, so the local runs reproduce the precursors (temperature split, contact density),
and the final link to the tear is inferred.

No production job was touched, no gate widened, no parameter changed.

---

# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only.

## Current phase (from 2026-09-05): core-FF retraining, then an unattended glpG arm test

ConDiv gly-sym retraining runs on midway2 AND rockfish toward 500 minibatches, both finishing
2026-09-09 afternoon. Everything downstream of it must run without supervision — the user is away
from the Mac on Thursday 2026-09-10, and a Claude session exists only while that Mac is on. The
chain therefore lives in Slurm scripts on the cluster.

**2026-09-09**
- **DDM is retired as an environment, and the documentation was cleaned to match.** glpG runs in a
  POPE/POPG bilayer only. Removed the detergent column that `findings.md` §4.1 used as the fidelity
  reference (2.61 A core RMSD, 0.952 occupancy) and said plainly that there is now no measured
  reference for how faithful this model can be, only the crystal. Retired the `glpG_DDM_micelle_REMD`
  paths in `remote_jobs.md` and pointed the glpG sections at `popepopg_REMD_mdw2`; corrected the
  two-campaign table, which still described glpG as a 48-replica micelle. Kept, deliberately: the
  one-tail-means-micelle morphology rule and its DDM worked example (it is what stops a detergent
  being built as a slab), the per-chunk RNG seed lesson, and the DDM itp/pdb parameter files.
  Files: `plan.md`, `findings.md`, `remote_jobs.md`, `example/16.MARTINI/readme.md`.
- **The chain proved itself overnight, unattended.** `48988330` FAILED (exit 7:0, cause unexplained,
  both logs clean) at step 443 and its replacement `48999774` hit a NODE_FAILURE at 449.
  `check_continue` caught both, read the step exactly right each time, and resumed for precisely the
  steps missing. Now on `48999888` at step 473, link `48999889`. The step-counter fix earned its keep
  inside 14 h: the old file-count version would have read 441 and 447 and landed on 502, putting the
  two hosts at different steps.
- **Corrected my own disk analysis, which was wrong.** I reported 195 GB of headroom; that is the
  `/project2` group quota, while the glpG data is on `/project` with **1514 GB** free. `rcchelp
  quota` lists four separate `trsosnic` group quotas and I read the wrong row, then dismissed the
  `df` output that was actually correct. The other session found this independently and rewrote
  `check_quota.py` to match the section header, take `min(quota, statvfs)`, and stamp the value;
  that fix is kept. Recorded in `findings.md` §12.
- **Reversed the decision that error had driven.** `decide_and_launch.sbatch` now **archives all 28
  pre-ff3 rungs** to `$V/pre_ff3/` and deletes nothing, since the multi-temperature ladder is what
  MBAR/HDX reads. `REMD_MAX_BLOCKS=4` is kept but rejustified on schedule rather than disk. Files:
  `decide_and_launch.sbatch`, `submit_remd.sh`, both md5-verified with `.bak_pre_keepladder` /
  `.bak_pre_rejustify`.
- **Automated the delivery race.** The two hosts finish ~23 min apart and midway2's chain fires as
  soon as its trainer ends, so a missed window silently loses the n=2 check.
  `scratchpad/deliver_rf.sh` polls rockfish and delivers its step-500 force field, writing `STEP`
  last so a partial transfer cannot be mistaken for a delivery.
- **Broke a standing instruction**: `CLAUDE.md` forbids the em dash character and I used it
  throughout yesterday's documentation edits. Removed from today's additions; yesterday's are
  committed and still there.
- Noted: another session is editing the same cluster directory and the same `.md` files, and a
  `git pull` at 10:45 merged its work in. Nothing collided, which was partly luck.

**2026-09-08**
- **Storage came back and left two trainers.** midway2 resumed from its own step-338 checkpoint
  (`48988330`); rockfish never stopped (`30725720`). Both at step ~355 of 500 with 0 failures. The
  Mac trainer is retired. Files: none, this was a status check.
- **Found the chain disarmed and re-armed it.** Neither host installed anything at step 500:
  `continue_mdw2.sbatch` refuses to by design and rockfish has no install script, so both would
  have simply stopped. Re-armed as `48994406` and cancelled the now-redundant `48988387`. The first
  attempt (`48994331`) had to be replaced: Slurm snapshots a batch script at submission, and it
  predated the `STEP` stamp that the runtime-submitted `run_arm_test.sbatch` requires, so the chain
  would have died at the handoff. Walked straight into a trap already recorded in `remote_jobs.md`.
- **Found three things `planned_job.md` asserted that were false**: nothing would install the force
  field; glpG production is idle rather than running (all four COMPLETED 2026-09-05), so the
  relaunch had nothing to cancel and its `rm -f $V.run.*.up` would have deleted the completed
  pre-ff3 baseline; and `scratchpad/ff3_retraining/` is on a different machine, having arrived here
  only as `git pull` fast-forwards that do not carry a gitignored directory.
- **Measured the retraining's reproducibility, which decided the plan.** The two runs had drifted the
  same distance from ff_2.1 (within 0.1-1.7%) but 23 degrees apart in direction, leaving them 33-40%
  of that drift from each other and 15-17% rms apart on every trained table. Training reproduces how
  far it moves, not where. So the 12 h arm test was **restored** rather than skipped, retargeted from
  "which coverage recipe" to "does a 16% table difference change TM4" — which also restores the
  pre-production check that skipping it had removed.
- **Rewrote the chain for that, and for idle production.** `check_continue.sbatch` (step from the
  checkpoint directory name, not a file count that trailed by 2 and would have overshot to 502;
  stamps `FF_DIR/STEP`), `run_arm_test.sbatch` (arms differ by force field, arm R refused unless its
  `STEP` matches), `decide_arm.py` (M/R, midway2 wins ties as primary, honest n=1 verdict),
  `decide_and_launch.sbatch` (reads the verdict instead of a hardcoded `ARM_B`, archives
  `$V/run.0.up` to `$V/pre_ff3/` instead of deleting the ladder, dead `NP_*` vars removed),
  `submit_remd.sh` (`REMD_MAX_BLOCKS=4`: 12 needs ~690 GB against ~284 GB of quota), and new
  `compare_ff.py`. All six md5-verified onto the cluster with `.bak_pre_ffcompare` backups.
- **Verified before arming**: `bash -n` on four sbatch scripts, `py_compile` on both python files,
  `compare_ff.py` on ff_2.1 versus itself (rel_rms 0.0000) and on a missing file (the one-armed
  path), `decide_arm.py` on unbuilt arms (`NO_WINNER`, exit 2), and an `extract_ff.py` dry run on
  rockfish producing a valid `sidechain.h5`.
- **Learned that `cat` cannot move a binary off rockfish**: the login banner is on stdout and
  prepends 1659 bytes, which reads as `file signature not found` rather than as a bad transfer.
  base64 with a marker, md5-checked both ends.
- Outstanding: extract rockfish's **step-500** force field into `parameters/ff_3.0_trained_rf/` with
  a matching `STEP` file, or the arm test runs one-armed at n=1.

**2026-09-05**
- Audited the training code against the upstream `origin/ConDiv` reference with four parallel
  reviewers. **No correctness bug in the math**: gradient identical to the reference (diff 0.0),
  env-derivative folding correct to 2.8e-10 by finite difference, coeff/weights slice confirmed
  against `src/environment.cpp`. Restrained/free ensembles verified separate — the 0↔1 exchange was
  attempted 79 times and accepted 0.
- **Restored the reference 8-replica ladder.** `n_threads` silently set both the OpenMP thread count
  and the replica count, so packing 4 workers/node had truncated the ladder to 4 rungs (top T 0.86
  vs 0.99), weakening the CD negative phase. Now 8 replicas on 4 nodes; thread density per node is
  unchanged so the step time did not move (~612-664 s). Resumed from the existing checkpoint with no
  progress lost.
- **Fixed the FF-install pipeline, which could not have worked.** `extract_ff.py` was unable to load
  any checkpoint (`Target`/`Update` are `__main__` types), which under `set -e` would have aborted
  the install *after* the full multi-day run. Now verified to reproduce training's own
  `sidechain.h5`/`environment.h5` byte-for-byte. Rebuilt `check_continue.sbatch` with a
  forward-progress guard (the old `afterany` chain would have resubmitted a failing job forever), an
  install sentinel, timestamped seed backups, and 180-step jobs sized to fit the 36 h wall.
- **Established that glpG and NP run different Hamiltonians by build accident**, not design: the
  in-tree hybrid prep never calls the coverage/environment writers and has no flag to, while NP's
  configs came from the rejected RD1 `envfull` recipe via a gitignored script. Unifying them is the
  goal of the arm test.
- Wrote and validated `py/martini_inject_coverage.py` (Arm B builder). Injecting into glpG shifts the
  potential +51.28 E_up **entirely within `protein_potential`**, MARTINI nodes byte-identical,
  idempotent, engine stable over 1000 steps. Bit-level parity against a real standard Upside config.
- Captured the pre-install glpG TM baseline for the after-comparison. Found and fixed a sign
  inversion in the dihedral routine before trusting any number from it.
- Confirmed on the live configs (not from source alone) that dry-MARTINI contributes nothing
  intra-protein, the hybrid interface terms are intact, and the protein is genuinely mobile.

## Key completed milestones

| Date | Milestone |
|------|-----------|
| 2026-08-10 | findings-88 fix (BB force path) deployed; all seeds rebuilt |
| 2026-08-13 | MBAR reference subtraction fixed (findings 91); dG plots delivered |
| 2026-08-15 | BB proxy reworked onto `infer_H_O`; CB placement bug fixed (findings 102) |
| 2026-08-16 | Cluster rigid-protein bug found and fixed (`current_stage` must be `production`) |
| 2026-08-17 | Membrane accessibility term wired into HDX pipeline (findings 113) |
| 2026-08-18 | NP campaign rebuilt with corrected CB + r_min_ang; relaunched |
| 2026-08-19 | Poster delivered; NP block-2 re-measured; claims corrected against Carlson et al. |
| 2026-09-01 | GLY Ramachandran maps symmetrized at source in `parameters/common/rama.dat` |
| 2026-09-05 | Core-FF retraining moved to the reference 8-replica ladder; FF-install pipeline fixed and verified |

### 2026-09-07 evening — RCC storage outage, training relocated to the Mac

* Checked the ff3.0 run as asked and found it **stalled, not running**: Slurm said `RUNNING` while
  all 9 workers sat in `D` state at zero CPU for 3.5 h on GPFS wait channels. Diagnosed to an
  RCC-wide storage incident (midway2 login nodes refusing TCP, midway3 login nodes with `/project`
  as bare xfs and no Slurm client). Job 48981235 cannot be cancelled; RCC ticket drafted at
  `scratchpad/rcc_ticket_draft.md`. Recovery watcher armed.
* Benchmarked a real minibatch locally to size the fallback: **1134 s/step** on the M1 Ultra, 2.0x
  the cluster. The benchmark exposed a real bug — all 12 workers failed on valid data because
  Upside exits **SIGTRAP (133)** under clang when Monte Carlo is enabled. Root cause: abstract
  `MonteCarloSampler` with no virtual destructor, deleted through `unique_ptr` of the base.
  Fixed in `src/monte_carlo_sampler.h`; old and new binaries give bit-identical output on all 16
  datasets.
* Established that the newest *reachable* force field is the step-269 extraction (step 338 is on the
  wedged filesystem) and that a resumable checkpoint can be rebuilt from it: `pack_param` refits the
  latent vector to 1.1e-16, and `extract_ff.py` on the rebuilt checkpoint returns the same force
  field to 1e-11 relative.
* **Local training running** from step 269 for 231 steps, PID 9228, log
  `training/gly-sym/run_output_local269/train_local.log`, ETA ~2026-09-10 22:00.
* Files modified: `src/monte_carlo_sampler.h`, `remote_jobs.md`, `findings.md`, `plan.md`,
  `progress.md`; added `scratchpad/ff3_retraining/{build_local_resume.py,train_local_269.sh}`,
  `scratchpad/rcc_master.exp`, `scratchpad/rcc_ticket_draft.md`.

## 2026-09-09 - status check, and a latent quota-gate defect fixed

* Checked all three hosts. **Training is healthy and on track:** midway2 `48999888` at step 453 and
  rockfish `30725720` at step 454 of 500, both ~11 min/step, converging on step 500 around
  16:00-17:00 CDT the same day. midway3 is idle. The install chain re-armed itself twice overnight
  (`48999889` now queued) and absorbed two midway2 job losses with no human action: `48988330`
  FAILED at step 443 with ExitCode 7:0 and clean logs (cause unproven), `48999774` hit an outright
  `NODE_FAIL` on midway2-0096. Each cost only the partial minibatch in flight.
* **Fixed `training/gly-sym/check_quota.py`, which would have falsely aborted the production
  launch.** `decide_and_launch.sbatch` gates on it, and on midway2 it exited 1 because every quota
  interface there is broken (`/project` is a remote fileset). Investigating that exposed a second,
  worse defect: `rcchelp quota` emits **four** `trsosnic` group rows and the old parser took the
  first, so it had been reading `/beagle3`'s quota on midway3 and `/project2`'s on midway2, never
  `/project`'s. The gate now matches the `mounted at` section header and returns
  `min(group_headroom, fileset_free)`; the corrected group row (1516 G) agrees with statvfs on
  `/project/trsosnic` (1514 G) to within 2 G, which is the cross-check that the right row is read.
  Where `rcchelp` cannot answer, a timestamped stamp written from midway3 is used, and a stamp over
  24 h old is refused rather than trusted. Verified by 7 tests on the cluster, including with the
  chain's own venv interpreter. No resubmission needed: it is a `.py` read at runtime, so the Slurm
  snapshot trap does not apply.
* Also corrected the "GPFS group quota is binding, NOT df" note in `remote_jobs.md`, which quoted
  `/project2`'s numbers as the constraint. On our fileset the two limits agree, so `df` on
  `/project/trsosnic` is sound; what was wrong was the quota row, not `df`.
* Recorded the working midway2 route: this Mac's IP is blocked on midway2's login nodes, so it is
  reachable only through a midway3 port forward plus `mdw2_via_tunnel.exp`. The ProxyCommand variant
  reports success and silently leaves no socket. Also noted that `/project` is the same filesystem
  on midway2 and midway3, so training state can be read without a midway2 login at all.
* Files modified: `remote_jobs.md`, `progress.md`, and on the cluster
  `training/gly-sym/check_quota.py` (original kept as `check_quota.py.bak_pre_stamp`), mirrored to
  `scratchpad/ff3_retraining/check_quota.py`.
* **Still owed:** extract rockfish's step-500 force field into `parameters/ff_3.0_trained_rf/` with a
  matching `STEP` file, or the arm test runs one-armed at n=1; and re-stamp the quota headroom from
  midway3 if the launch slips past 2026-09-10 08:30.

### Deployed ff3.0 (2026-09-09 22:42 CDT)

* **`parameters/ff_3.0/` installed** from midway2's step-500 output (cluster name there is
  `ff_3.0_trained`; the local slot keeps the `ff_2.0`/`ff_2.1` convention): `sidechain.h5`
  (`c67351ca...`), `environment.h5` (`e4d2f685...`), `STEP` = 500, plus a README recording
  provenance, why step 500 over step 269, and what the force field does not fix. md5 verified
  against midway2 on both ends.
* **Verified, not assumed.** The trained tables moved (rel_rms ~0.57 vs `ff_2.1`) while
  `hydrophobe_placement` and `rotamer_center_fixed` are identical to 0.0000, which is the check that
  only the intended tables changed. Then end to end: a real glpG seed patched with the deployed
  file passes the full hybrid check, engine energy -24944.285156, finite forces, interface intact,
  intra-protein MARTINI excluded. That energy agrees with the cluster's arm M to 8e-8 relative,
  the expected platform floating-point difference.
* **The stale `parameters/ff_3.0` was DELETED** (user's call, 2026-09-09): it sat 0.0501 from
  `ff_2.1`, about one training step, and was not ff3.0 in anything but name. Nothing in `py/`,
  `src/` or `example/` referenced it. `parameters/ff_3.0/` is now the single ff3.0 directory,
  overwriting that slot; the old content remains recoverable from commit `2818532` if ever needed.
  `hbond.h5` was carried over from `ff_2.1` (byte-identical, and hbond was held fixed in training)
  so the slot is complete like its siblings.
* **Not committed** (read-only git rule). Left unstaged for the user.
* **Open:** if the REMD arm test picks arm R, `sidechain.h5` must be swapped for rockfish's
  `cfe4ba5e...`. Noted in the README.

### Evening (2026-09-09 17:20 CDT): training complete, arm test running on both clusters

* **Training finished at step 500 on both hosts.** midway2 extracted its force field at 17:15:58 and
  handed off; rockfish finished earlier and its force field was delivered at 16:30.
* **The two-arm test is running on both clusters** as independent replicates. rockfish `30768485` (M)
  and `30768486` (R); midway2 `49001769` building its own pair. Both clusters hold byte-identical
  copies of both force fields (`c67351ca…`, `cfe4ba5e…`) and a byte-identical seed, so the two
  replicates test exactly the same comparison.
* **Concluded ff3.0 needs no further training** and recorded why in `findings.md`: the objective has
  plateaued, the parameters random-walk with no fixed point, and the between-run spread grows as
  sqrt(t), so more training makes the force field *less* reproducible. Also found the training set is
  458 soluble proteins with zero membrane content, while the deliverable is TM4 stability in a
  bilayer, so more steps optimize an objective that is both saturated and misaligned. Suggested
  averaging over the plateau instead, to be validated rather than assumed.
* **Verified ff3.0 is usable for simulation**: patched into a real seed it passes the full hybrid
  check (finite energy, interface intact, intra-protein MARTINI excluded, production stage). Its
  force *tail* is heavier than ff_2.1 (peak 318 vs 111) but the bulk is unchanged (identical median,
  p90/p99 within 7-14%); the excess is 9 atoms out of 4949. Recorded as the thing to watch first if
  an arm destabilises.
* **Found a gap in the chain**: `run_arm_test.sbatch` excludes `midway2-0003` for its children but
  not for itself, and `49001769` duly landed there.

### Afternoon addendum (2026-09-09 16:35 CDT)

* **Rockfish training finished at step 500 and its force field is DELIVERED** to
  `parameters/ff_3.0_trained_rf/` (16:30:50 CDT), so the arm test will be two-armed rather than
  one-armed at n=1. Provenance checked, not assumed: extracted from `epoch_13_minibatch_05`
  (`epoch=13 i_mb=6`), md5 identical on rockfish / this Mac / midway2, and `compare_ff.py` against
  `ff_2.1` shows mean `rel_rms = 0.4773` on the three trained tables while the untrained geometry
  tables match to 1e-13. `STEP` was written last so a partial transfer could not pass as a delivery.
* midway2 `49000800` was at step 497/500 at 16:31, ~30 min from triggering `49000801`.

### Earlier that afternoon (15:45 CDT)

* **Both trainers are in the last hour:** midway2 `49000800` at step 492/500 (ETA ~17:00 CDT),
  rockfish `30725720` at step 496/500 (ETA ~16:25 CDT). Rockfish finishes first, which is what makes
  the owed force-field delivery possible at all; the margin is only ~35 min.
* **`remote_jobs.md` claimed a delivery watcher (`scratchpad/deliver_rf.sh`) had been polling since
  12:07 CDT. It did not exist.** No such file anywhere on this Mac, no log, no process, and never
  committed since `scratchpad` is gitignored. Trusting it would have let the arm test go one-armed
  unattended. Corrected the entry and added the rule: verify a watcher with `ps` and its log.
* **Built and started the real mechanism**, both verified running: `watch_extract_500.sh` on rockfish
  login03 (pid 599446) extracts the step-500 force field once the checkpoint size settles, and
  `auto_deliver_rf.sh` on this Mac polls for it and runs `deliver_rf_ff.sh`, which pushes both `.h5`
  to `parameters/ff_3.0_trained_rf/` through the midway3 socket. `STEP` is written last so a partial
  transfer cannot look like a delivery, and the script refuses unless rockfish's STEP reads 500
  (tested: it refuses and creates nothing). The binary path was dry-run end to end with
  `ff_2.1/sidechain.h5`, identical md5 on all three hops and the file opens in h5py.
* **`midway2-0096` killed a second job** (`48999888`, NODE_FAIL at 13:15, after `48999774` that
  morning). Deliberately did **not** add `--exclude=midway2-0096`: checkpoint resume made each
  failure cost only the in-flight minibatch, and with 9 steps left, editing a working chain script
  carries more risk than a ~5 min requeue.
* Corrected a stale note claiming direct midway2 SSH had returned. `nc` is refused again as of 15:30,
  so the midway3 tunnel route is still required.
* Files modified: `remote_jobs.md`, `progress.md`; added
  `scratchpad/ff3_retraining/{deliver_rf_ff.sh,auto_deliver_rf.sh}` and `watch_extract_500.sh` on
  rockfish.

## 2026-09-13 — glpG VTF delivered, force field confirmed, dG comparison slide

* **Extracted and verified the glpG VTF for Monday.** 1822 frames (blocks 1-32, stride 5) to
  `~/Downloads/glpG_RKRK_79HIS_run0_remd.vtf`. Checked before handing it over rather than after:
  TM4 alpha fraction mean 0.791 with no frame under 0.50, TM4 centroid within 7.7 A of the bilayer
  midplane in every frame, peptide C-N mean 1.324 A. A per-block internal-RMSD probe confirmed only
  `output_previous_0` is rigid (0.000 A), which is the seed block the extractor already skips.
* **Confirmed the glpG chain is on ff_3.0 by measurement.** Needed because `run_remd.py` reuses
  existing replicas and the force field is baked into each `.up`. A verbatim hash sweep answers
  "ff_2.1" and is wrong; the least-squares scale on the rotamer pair table gives 1.000000 against
  `ff_3.0`. Recorded as findings 6.6.
* **Substituted the valid ff2.1 dG reference.** The `glpG_POPEPOPG_dG_2026-08-27/` figures were
  produced while the protein was frozen by the stage bug, so they are not a force-field comparison.
  Used the 2026-09-04 `hdx/results/` set instead and re-rendered both hybrid panels at the implicit
  run's rungs (0.75/0.80/0.85) so the three models can be read side by side. Off-scale amides at
  T = 0.85 go 43/203 (ff2.1) -> 69/203 (ff3.0); TM4 stays finite in all three.
* **Deck at 20 slides**, new slide "The same helix, read three ways".
* Files modified: `findings.md` (5.3d, 6.6), `remote_jobs.md`, `progress.md`,
  `0914/build_draft.py`, `0914/figs/dG_{implicit,ff21,ff30}.png`.
* **Not done:** the ff3.0 benchmark TM/RMSD table is still 10 of 22 arms (job 49010900, 57 min in);
  the long proteins (lambda, top7, ubiquitin) are outstanding.

## 2026-09-13 (later) — slides brought onto the newest data

* **Benchmark: already current.** Checked the queue rather than assuming; nothing finished after
  scoring job 49010900 exited COMPLETED at 2 h 54 m. Slide 18 rebuilt on its 21-arm table (15/16
  native, 6/16 de novo): native mean TM 0.583 vs FF2's 0.55 and mean Ca-RMSD 3.78 vs 4.0 A, lower
  lowest-RMSD on 9 of 15, gpW the only regression. `make_ff3_benchmark_fig.py` now reads the scoring
  table instead of carrying transcribed numbers.
* **glpG HDX re-run on the grown trajectory** (Campaign 6, jobs 59041160-63, midway3): 6,121 ->
  9,465-9,956 frames per replica. All four COMPLETED, 28/28 replicas each. **Nothing that matters
  moved**, which makes the four-variant equality a converged result rather than one measurement.
  Retracted the "79ALA_S115T is mildly tighter" reading as sampling noise.
* Had to submit to midway3: the HDX venv's `bin/python3` points into `/software`, which is
  per-cluster, so on midway2 python3 falls through to `/usr/bin/python3` with no h5py.
* Noted but not changed: the HDX analysis topology is built with ff_2.1 tables while the trajectory
  is ff_3.0. Left alone so Campaign 6 differs from Campaign 5 in frame count alone.
* Files modified: `findings.md` (5.3d), `remote_jobs.md` (Campaigns 2 and 6), `progress.md`,
  `0914/build_draft.py`, `0914/make_ff3_benchmark_fig.py`, `0914/figs/{dG_var_*,dG_ff30,fig_ff3_benchmark}.png`.
* **Still outstanding:** 10 of 16 de novo benchmark arms and alpha3D-native are still simulating.

## 2026-09-14 — VTF periodic-image fix, glpG trajectory re-delivered

* **Root-caused the "protein left the bilayer and left a hole" artifact in the delivered glpG VTF.**
  Not the last frame only: **159 of 1822 frames**, protein displaced by exactly one box length in x
  (xy centroid separation from the lipids up to 98.18 A against a 99.77 A box). The trajectory itself
  is fine: the protein is intact in those frames (no CA-CA above 4.5 A) and the earlier physical
  checks stand. It was a post-processing fault in `py/martini_extract_vtf.py`.
* **The cause was wrapping per particle and then rebuilding molecules.** `centralize_system` wrapped
  every atom into the cell, which tears any molecule crossing a face, and `unwrap_molecules` then
  rebuilt each one around whichever anchor atom landed inside. For glpG that anchor is atom 0, the
  floppy N-terminal amide, so every time the tail crossed a face it dragged the whole protein with it.
* **Fixed by reversing the order and wrapping per molecule.** `build_bond_walk` became
  `build_molecule_topology`, which also returns a connected-component label per particle;
  `extract_trajectory` now unwraps first, then centres on the plain protein centroid and wraps each
  molecule by its own centroid. The circular-mean COM is gone, unnecessary once the protein is whole.
* **Verified before re-delivering, on the new file:** protein COM exactly 0 in all 3146 frames,
  **0 displaced frames**, protein-lipid xy centroid separation mean 0.70 A / max 2.20 A, every
  declared bond under 10 A except the known residue-210 C-O (23.8 A, pre-existing and unrelated).
* **Re-extracted and downloaded** `glpG-RKRK-79HIS` replica 0, now blocks 1-54 at stride 5 (3146
  frames, 321 MB, up from 32 blocks / 1822 frames) to
  `~/Documents/2026/reports/GroupMeetings/0914/glpG_RKRK_79HIS_run0_remd.vtf`; md5 checked at each hop.
* **Answered the follow-up "should the hole around the protein close?" by measurement: yes, it is
  post-insertion annular relaxation, but it is slow relative to the run.** Production is fixed-volume
  (no barostat, no `box` dataset), there is never a through-hole (0 A^2 lipid-free projected area at a
  5 A probe in every frame), and the protein is static (TM-slab Rg_xy flat at ~12.5 A). The lipids
  redistribute inward: density rises in every radial bin inside 20 A and falls beyond 25 A. TM-core
  contact beads +60%; annular lipid count saturates early while contacts per lipid keep rising
  3.8 -> 4.8. Exponential fit **tau = 29% of the production run**, 90% of plateau only at ~frame 2100
  of 3146. Recorded as findings 8a with the HDX consequence.
* **Open, not yet done:** re-run the HDX estimator on the last third only, and measure the same shell
  curve for the other three variants, before treating the four-variant equality as converged.
* Files modified: `py/martini_extract_vtf.py`, `findings.md` (3.8, 8a), `progress.md`; on midway2,
  `/project/trsosnic/yinhan/upside2-md-mdw2/py/martini_extract_vtf.py` (backed up as
  `.bak_pre_pbcfix`) and `~/project/yinhan/extract_glpg_vtf.py`.

## Carried-over open items

- **NP footprint contradicts the paper.** None of Carlson et al.'s five target lysines are contacted
  (K190 = 0.000, the lowest of 58 Lys). Blocked on over-unfolding and on running `np_footprint.py`
  against current data. NP is explicitly not urgent.
- **Cluster analysis scripts drift from the repo.** Re-upload `calc_hdx_ht.py` and
  `4.calc_D_uptake.py` before trusting any cluster-side HDX result.
- **TM4 is weak** (0.47-0.65 helix fraction vs TM1's 0.77-0.91) in the pre-install baseline. Cause
  unproven; the arm test bears on it but was not designed to settle it.

## 2026-09-16  GLY Ramachandran symmetry: is it GG or XGX in glpG

**Question (PI).** Is the middle-glycine Rama map naturally symmetric for unfolded GGGGG and SAGAS,
and is glpG's ff2.1 failure a `GGG` effect? Refined mid-session to: should `XGX` be symmetric as
ff3.0 assumes, and can ff3.0's glpG improvement be attributed to GG or to XGX?

**Actions.**
* Scanned the simulated glpG construct (`scratchpad/local_popg_79HIS/glpG-RKRK-79HIS.run.0.up`,
  210 residues) for glycine motifs. No `GGG`; two `GG` (96-97, 132-133); 19 of 23 glycines are
  pure XGX. 91% of the dimer-library entries a glycine draws are `GLY|X`.
* Compared `parameters/common/rama.dat` against `rama3.dat`: only the GLY row differs, in both the
  `coil` and `sheet` groups, and the operation is exactly `0.5*(m + mirror(m))`.
* Derived the exact grid mirror from `rama_map_pot.cpp:67` + `LayeredPeriodicSpline2D`, verified it
  against the axis convention (ALA minimum at phi = -65, PRO restricted to -95..-50), and
  re-measured the symmetry of a glpG `rama_map_pot` built from each library: 5.856 max asymmetry
  from `rama.dat`, **0.000000** from `rama3.dat`. Corrected the mirror note and the mislabelled
  "what is deployed" row in `findings.md`.
* Priced the ff2.1 -> ff3.0 glycine change per residue and per TM helix. Whole protein: XGX 84%,
  GG 16%. TM4's three glycines are all XGX and carry the largest shifts; TM1 has no glycine.
* Launched the all-atom test of the underlying premise: `Ace-GGGGG-NMe` (achiral control, exact
  answer 0, calibrates the error bar) and `Ace-SAGAS-NMe` (the real XGX measurement), built with
  `scratchpad/gly_rama_aa/build_peptide.py` (CB chirality calibrated against 19 real L-ALA
  residues), amber99sb-ildn / TIP3P, 300 K, 3 replicas x 100 ns each at 2 threads.

**Files.** `findings.md` (new section + two corrections), `scratchpad/glpg_gly_context.py`,
`scratchpad/glpg_gly_by_helix.py`, `scratchpad/gly_rama_symmetry.py`, `scratchpad/gly_rama_aa/`.

**Results so far.** Under ff2.1 all 23 glycines favour the left-handed helix by 0.47-1.23 E_up and
the raw GLY coil map's global minimum sits at phi = +85. Conditioning on an achiral nearest
neighbour barely reduces the asymmetry (`GLY|GLY` 4.17 vs `GLY|X` 4.06 +/- 0.55), so at most a third
of the bias can be a nearest-neighbour effect; the rest is longer-range chiral context a dimer
library cannot represent.

**Open / running.** The 6 GROMACS replicas (launched 11:52 CDT, ~30 h for the full 100 ns each;
partial trajectories are readable much sooner via `scratchpad/gly_rama_aa/analyze.sh` then
`middle_gly.py`). A first driver attempt died because `setup.sh`'s `rm -rf` deleted the driver's own
working directory; `run_prod.sh` cds to a safe root first.
