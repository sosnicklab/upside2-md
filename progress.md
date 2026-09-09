# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only.

## Current phase (from 2026-09-05): core-FF retraining, then an unattended glpG arm test

ConDiv gly-sym retraining runs on midway2 AND rockfish toward 500 minibatches, both finishing
2026-09-09 afternoon. Everything downstream of it must run without supervision — the user is away
from the Mac on Thursday 2026-09-10, and a Claude session exists only while that Mac is on. The
chain therefore lives in Slurm scripts on the cluster.

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

## Carried-over open items

- **NP footprint contradicts the paper.** None of Carlson et al.'s five target lysines are contacted
  (K190 = 0.000, the lowest of 58 Lys). Blocked on over-unfolding and on running `np_footprint.py`
  against current data. NP is explicitly not urgent.
- **Cluster analysis scripts drift from the repo.** Re-upload `calc_hdx_ht.py` and
  `4.calc_D_uptake.py` before trusting any cluster-side HDX result.
- **TM4 is weak** (0.47-0.65 helix fraction vs TM1's 0.77-0.91) in the pre-install baseline. Cause
  unproven; the arm test bears on it but was not designed to settle it.
