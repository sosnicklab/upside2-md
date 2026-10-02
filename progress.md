# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`.

---

## 2026-09-18 to 09-30: basin-offset ff_3.0 (superseded)

* Track B: 2D AWH on 40 capped glycine dipeptides (GROMACS, 400 ns) finished: left neighbours
  -0.20 nats, right -0.11; ff14SB agrees with ff99SB-ILDN to 0.045 nats (findings 9r). Track A, a
  learned glycine row ("ff3.1"), was invalidated with the trainer it ran on.
* 09-24: the ConDiv port in use trained FF1's functional form (spline burial, no backbone
  desolvation, no unfolded-state objective; findings 9t, 9u). Replaced by Kleinmann's port of Peng's
  FF2 trainer restored to the SI (findings 9v); the Phase 1 fixed-point check from ff2.1 left 8 of 9
  groups at a fixed point.
* 09-27/28: Ramachandran maps redesigned as per-pair basin offsets on NDRD (findings 1.8-1.13);
  `training/rama_basin.py` implemented and verified; `ff30_basin` started on midway2 and was rewound
  to step 19 three times (MAP Newton step, GLY|GLY sheet symmetry, the reduced 158-offset set).
* 09-30 02:25: the gate released ff_3.0 to both cluster trees with the rotamer-BP fix installed.
  The pre-proline free-native gap turned out to be composition, not a map defect (findings 1.14).
  glpG `run_remd.py` now counts a block per job id, so a Slurm requeue no longer burns one.
* 09-30: the equal-depth glycine probe pulled back toward alpha_L (findings 1.15). At 20:34 the user
  cancelled all ff_3.0 validation (Peng arms, glpG chains; STOP files placed): its glycine maps favour
  alpha_L everywhere and TM4's helical glycines flipped.

## 2026-10-01: glycine measurements for the replacement design

* Coverage of a glycine H-bond term (`checks/gly_mismatch_by_hbond.py`): of the glycine-specific
  helical loss, 62-66% is in glycines with their own H-bond, 25-32% in glycines spanned by a
  short-range bond, 2-9% in neither; natively left-handed glycines show no glycine-specific deficit;
  every helical glycine in glpG's seed, TM4's included, has its own bond (findings 1.16).
* Per-type misses in the same native basin and H-bond state follow helix and beta propensity, so a
  trained per-type correction would flatten real propensity. Residue counts per basin cannot correct
  the map: the implied glycine handedness depends on the reference residue (-2.61 to -1.16) and on
  burial (-1.95 to -0.17) (findings 1.16).
* Proline: central prolines hold their basins better than any residue; residues before a proline
  lose +0.024 (extended) and +0.059 (helical). User correction: glycine and pre-proline are separate,
  and pre-proline is not an established problem; Phase 6 parked (findings 10.8).
* Experiment: our AWH Gly-Gly surface against the GGG populations (Andrews 2020 SI): pPII 0.32 vs
  0.46, alpha 0.09 vs 0.06; in the matched system ff14SB is within 0.06, so no change. No handedness
  data exists. The helix benchmark is dropped (Pace & Scholtz is an 11-system average).

## 2026-10-01: Phase 8 implemented and training

* Library: `training/build_gly_library.py` (new); `parameters/common/rama31.dat` rebuilt (GLY|GLY
  symmetrised, GLY|right|PRO its own surface, other GLY|X pooled, coil = sheet, reference correction
  subtracted). Builder checks and an end-to-end `.up` check pass; the midway2 build is
  byte-identical.
* Trainer: offset training removed from `ConDiv.py`, `extract_ff.py` and `convergence_gate.py`;
  `rama_basin.py` reduced to basin populations; `verify_rama_basin.py` deleted. README, up.md,
  architecture.md, GLY_sym.md updated. Deployed to midway2 with backups.
* `ff30_gly` initialised from ff2.1. midway2 refused the chain (`AssocMaxCpuPerJobLimit`) until the
  allocation reached the scheduler at ~11:15. A start on midway3 `amd` was cancelled by the user
  after 14 min: no CPU jobs on the GPU allocation (findings 10.5).
* Local run: `CONDIV_LOCAL_WORKERS` added to ConDiv's launcher; new `training/check_step.py`
  (per-step physical health) and `move_run.py` (checkpoint path remap). The Mac Studio ran
  09:16-12:16 and the cluster chain resumed at step 1 from the local step 0. Checkpoints written
  under NumPy 2 did not load under the cluster's 1.23.5; `move_run.py` now converts them
  (findings 10.10).
* Repo cleanup for redistribution (plan.md Phase 9): campaign files moved from `py/` and
  `training/` to `scratchpad/redistribution_cleanup_20261001/`; `env.sh` lost its midway3 branch,
  `gate_or_continue.sh` stops at convergence. Compile, env and gate-branch checks pass.
* Release by checkpoint selection (user): the midway2 gate stops at convergence with nothing
  released, and `validate_ff.sh` takes the chosen checkpoint. Sandbox-tested on midway2.
* Selection panel built and running: `ff3_selection/panel.py` (prepare, run, aa, select) and
  `panel.sbatch`; 44 domains prepared (2,823 residues, 198 glycines). Tested locally (phi/psi equal
  the engine's `rama_coord`) and end to end on midway2.
* Bottom-up glycine fallback proposed (plan.md; findings 1.17). Charron et al. 2025 read: their 50
  CATH domains are ff99SB-ILDN/TIP3P/300 K, the map's own force field.

## 2026-10-02: glycine H-bond offsets and damped side-chain step (Phase 8 revised, approved)

* Diagnosis behind the revision (findings 1.17): single-group resets of ff30_gly's epoch 0 put the
  helix and fold weakening on the side-chain pair update, a noise-driven random walk under Adam, and
  the helical-glycine loss on the shared H-bond drift. Both readings of the DSE threshold give the
  same DSE.
* Engine: `HBondEnergy` takes 12 + 3 per residue class parameters with a per-residue
  `residue_class` (`src/hbond.cpp`). A 12-entry config is bitwise equal to the old build (energy,
  forces, hbond, rotamer and coverage derivatives); zero offsets bitwise equal to 12 entries; offset
  derivatives and forces match finite differences (1ubq, local). The first layout changed a shared
  sum's last bit under -ffast-math (findings 10.12).
* Config writer and tools: `py/upside_config.py` (class_restype, residue_class, hb_scale on the
  offsets); `training/ConDiv.py` (field `hbg`, SARW zeroes it, rot lr 0.025); `extract_ff.py`,
  `check_step.py`, `training/README.md`, `up.md`; `patch_glpg.py` carries the offsets into a seed.
* One real ConDiv step, local, 1ga3 (8 glycines), 15-entry init: hbg gradient [-0.91 -0.52 +0.81],
  first Adam step +-0.01, rot step rms 0.0103 at lr 0.0125, hbond.h5 written back with `GLY`.
* Deployed to midway2 by two detached, gated scripts (`checks/hbg_deploy_20261002`), since the
  master socket kept dropping (our own login-node load, findings 10.6): files with backups, build
  in `obj_hbg`, parity on 1ga3 (engine and a 200-unit run bitwise against the installed build; zero
  offsets bitwise 12 entries), install by rename, `ff30_glyhb` initialised (19 x 24, hbg 0, lr rot
  0.0125, hbg 0.01). ff30_gly stopped 11:51, the new chain submitted 11:52. A staging slip left
  `upside_config.py` non-executable for ~2 min (findings 10.11); nothing ran in the window. The
  selection watch now follows ff30_glyhb (tags h00, h01, ...).

## 2026-10-02: documentation and training/ cleanup (user request)

* `findings.md` 6,819 -> 3,886 lines: superseded history (FF1-form trainer, the symmetrised ff3.0,
  ff3.0C, Track A, the AWH convergence saga, the pre-temperature-fix TM4 hunt) condensed into
  1.18, 4.4, 4.5, 6.7, 9e, 9s, 9w and 11d; duplicates kept once; §3 and §9 reordered; §10 renumbered
  10.1-10.13; every label cited elsewhere still resolves. `plan.md` 364 -> 224, `progress.md`
  250 -> 107, `remote_jobs.md` 822 -> 688, `GLY_sym.md` 537 -> 364, `architecture.md` restated
  against Phase 8. Memory: 4 notes deleted, 11 rewritten, index rebuilt.
* `training/`: `env.sh`, `train_chain.sbatch`, README and the ConDiv docstring made site-neutral
  (plan.md Phase 9); `py/__pycache__` and `training/__pycache__` deleted. Checks: `bash -n` on the
  three shell scripts, `py_compile` of ConDiv.py, `source training/env.sh` then importing
  `upside_engine`, `run_upside` and `rama_basin`. The midway2 tree is untouched.
* `CLAUDE.md` (user request): the compile section now builds on midway2 only, with the modules,
  cmake and Eigen paths of the 10-02 deploy, in a fresh `obj_<tag>` installed by rename; the Slurm
  section uses `env_shared.sh` or the training tree's `env.sh` instead of the old `python/3.11.9`
  recipe; the SSH section covers both clusters with the check-first rule; the hard-coded-identity
  rule lost its dated incident note (findings 1.7). `remote_jobs.md` §0's tunnel fallback now names
  the script that exists, `scratchpad/mdw2_via_mdw3.exp`.
* Scripts cut to what the workflows need (user request; plan.md Phase 9): `training/` 12 -> 5 files
  (`rama_basin.py`, `extract_ff.py`, `convergence_gate.py` into `ConDiv.py` as functions and the
  `extract` / `gate` commands; `gate_or_continue.sh` into `train_chain.sbatch`; `check_converged.py`
  deleted; `build_gly_library.py`, `check_step.py` to scratchpad); `py/martini_gen_params.py` into
  `martini_build_tables.py`'s `__main__`; root `analyze_sc_pairs.py` to scratchpad.
  `py/martini_hdx_project.py` now copies the coupled hybrid potential instead of re-scoring a
  protein-only one (the root cause `write_hybrid_energy.py` patched by hand). `example/00.AnalysisScripts`
  is left exactly at HEAD (user: those files follow master), with `combine_hdx_protection.py` and the
  accessibility script unchanged. Tests: midway2 jobs 49143118 and 49143223
  (`/project/trsosnic/yinhan/checks/merge_test_20261002`) gave byte-identical `extract` output and an
  identical gate table and verdict; basin functions bitwise equal; a local 1rkl projection keeps the
  source potential exactly.
* MPI check (user question): Upside needs no MPI and cannot use it (findings 10.13).
  `example/01.GettingStarted/0.run_mpi.py` was broken by design (two `mpirun` ranks collide on the
  same files, 0 frames, exit 0) and `run_1UBQ_local.sh` was `0.run.sh` with another `pdb_id`; both,
  and their 1UBQ PDB copy (master keeps one in `example/07.MoreRestraints/pdb/`), moved to scratchpad,
  so `example/01.GettingStarted` matches master. `module load openmpi` dropped from
  `example/16.MARTINI/run_sim_1afo_full.sh`.
* 12:50-15:00: threshold test closed (both SI statements give the same DSE; findings 1.17), held
  split tasks cancelled; e02 panel held for helix; ff30_glyhb started 13:25:41, steps 0-1 healthy,
  hbg [+0.019 +0.019 -0.019] after step 1. `remote_jobs.md` gained "Resume here, from any computer":
  what is and is not in git, how to connect, the watch by hand, the timeline to release, how to read
  the h tables, and the cluster/repo `training/` layout warning.

## 2026-10-02: rotamer-BP validation figure (for Tobin)

* All 16 REMD arms confirmed complete. Static-frame harness rebuilt from the 09-26 transcript (its
  scratchpad was gone) and rerun on fresh frames: early stop in 46% of 728 frames, |dE| median
  7.5e-4 E_up, max 0.157 (one real early stop in 1afh, 0.2 kT), force error median 8e-5 of RMS force.
  REMD arrays rerun on the midway2 login node (user's call; the compute-node job would have started
  ~20:11), output identical to `analysis_20260929.txt`. 11-panel figure in
  `~/Downloads/bp_validation.{png,pdf}`; data and scripts in `bp_validation/static/` (remote_jobs.md §0c).
