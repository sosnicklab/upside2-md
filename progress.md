# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only.

## Current phase (from 2026-09-05): core-FF retraining, then an unattended glpG arm test

ConDiv gly-sym retraining is running on midway2 toward 600 minibatches. Everything downstream of it
must run without supervision — the user is away from the Mac on Thursday 2026-09-10, and a Claude
session exists only while that Mac is on. The chain therefore lives in Slurm scripts on the cluster.

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

## Carried-over open items

- **NP footprint contradicts the paper.** None of Carlson et al.'s five target lysines are contacted
  (K190 = 0.000, the lowest of 58 Lys). Blocked on over-unfolding and on running `np_footprint.py`
  against current data. NP is explicitly not urgent.
- **Cluster analysis scripts drift from the repo.** Re-upload `calc_hdx_ht.py` and
  `4.calc_D_uptake.py` before trusting any cluster-side HDX result.
- **TM4 is weak** (0.47-0.65 helix fraction vs TM1's 0.77-0.91) in the pre-install baseline. Cause
  unproven; the arm test bears on it but was not designed to settle it.
