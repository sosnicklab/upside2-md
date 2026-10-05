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
* MPI check (user question): Upside needs no MPI and cannot use it (findings 10.14).
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
  `~/Downloads/bp_validation.{png,pdf}`. Everything the Mac had (scripts, the source patch as
  `make_src.py`, the build recipe as `build_libs.sh`, sources, libraries, SciencePlots styles,
  results; 325 files) mirrored to `bp_validation/static/`, md5 manifest identical (remote_jobs.md §0c).

## 2026-10-04: ff3.0 round-3 check, local (plan.md Phase 10)

* Panels h00 and h01 (ff30_glyhb epochs 0 and 1) both worse than the run's start in helix; release
  held. Lambda ff2.1 arms finished (three at target, G46A/G48A native in its last chunk): wild-type
  helix 3 settles near 60% helical from both starts, G46A/G48A holds it (findings 11d).
* User: the TM4 cause (placement read as energy) is settled; continue training; test the current
  stage locally on glpG TM4 and lambda; propose a design for both. New rule in plan.md Key Decisions:
  probe every new trainable term from a balanced start before committing epochs.
* Phase 10 A done: local engine rebuilt (`obj/*.bak_sep10_20261004` kept), parity with both cluster
  builds by initial potential (lambda ff2.1 -196.48, lambda e02 -193.92, glpG ff2.1 -24914.84); e02
  extracted on midway2; patch gate passed 2.6e-6; seeds patched for ff21_released, ff21_awh, e02,
  e02_gly0 (`scratchpad/ff3_local_test/`).
* Phase 10 B started 20:09: 12 glpG 79HIS runs, T 0.80, 4000 time units, 3 seeds per force field.
  The first launch died silently (`source.sh` fails under `set -u`); relaunched with `set -eo pipefail`.
* Charron et al.'s archive finished downloading 21:16 (md5 verified). The CATH and octapeptide h5
  were extracted to `~/Downloads/charron2025/`.
* User correction: the problem is the Rama map; an H-bond-term design (achiral glycine offsets, a
  bottom-up E_gly fit) was rejected, and its runs and analysis were stopped (findings 10.13; memory
  fix-the-rama-map). Phase 11 was rewritten around the map.
* The octapeptides were run with adaptive sampling (about 100 segments of about 10 ns per peptide).
  The L residues' alpha_L decays inside the segments; glycine's populations do not change. Interior
  glycines give ln(aR/aL) = -0.70 +- 0.04, against -0.11 for rama31 and -1.1 to -1.7 for NDRD
  (findings 1.19). Four literature surveys were added to the same section.
* Phase 10 B interim (2,000 tu): ff21_awh holds TM4; e02 flips GLY143.
* Started (`scratchpad/ff3_local_test/`):
  * `ff21_oct`, rama31 with the octapeptide surface for GLY|X. Patched from the live seed; only the
    23 glycine maps differ from ff21_awh. glpG TM4 runs, 3 seeds, from 21:22.
  * The map-fitting step M1 on 376 octapeptides in Upside (`ibi_run.py`, `ibi_compare.py`), from
    rama31, from 21:30.
  * lambda REMD under ff21_awh, queued to start when Phase 10 B ends.

## 2026-10-04/05 night: BioEmu glycine map and the push probe (plan.md Phase 11)

* Literature (six sub-agent searches, five paywalled papers read in full from the user's
  downloads): no experiment resolves glycine's alpha_R/alpha_L next to L residues; gas-phase QM of
  Ac-Ala-Gly-Ala-NHMe (RHF/3-21G) finds no helical glycine minima at all; Childers 2016 (Levitt
  force field) finds glycine between L-alanines indistinguishable from GGGGG (achiral). One agent put
  the user's email into an Unpaywall query, once; reported.
* Upside never designed the Rama transition region for any residue (Jumper thesis 4.3.1; NDRD tails
  6-10 E_up for all types; contrastive divergence cannot train it); decision: leave the top as the
  library has it (findings 1.19).
* BioEmu plain MD (Zenodo 15641199) replaces Charron's adaptive frames as the target (glycine biased
  by 0.16-0.28 there); all residue types compared with NDRD, only glycine's error flips helices in
  the panel. training_a_cg_model.zip deleted (extracts kept, user).
* Fit of every central-glycine entry in Upside: unit convention corrected to T_up = 1 (the first fit
  at 0.8557 was 1.169x too small; it and the probe built on it were discarded); GLY|GLY and
  GLY|right|PRO now fitted on their own BioEmu contexts (user); a sheet-entry bug found and fixed
  after passes 1-3 (findings 1.19); corrected passes 4-6 running.
* Push probe set up locally and verified against midway2 (inputs, ff2.1 files, engine energy on
  1ga3 identical); reference at the rama31 start d -0.105 [-0.150, -0.057]; midway2 would start ~39 h
  later. Gate (user): CI entirely below 0 -> train with the frozen map; TM4 checked per epoch.
* midway2: ff30_glyhb chain and panel h02 cancelled 21:58 (user); lambda G46A/G48A native complete
  01:24 (helix 3 aR 0.97, ~6.9 A plateau; findings 11d); no cluster jobs, watch cron deleted.
* Fit converged after the sheet fix (map 6: X-G-Y within 0.01 per basin of BioEmu in Upside). Push
  probe on map 6, 72 proteins: d -0.014 [-0.041, +0.012] (rama31 start -0.105); TM4 on the untrained
  start: no glycine flips, TM4 0.95. By the pre-agreed rule "do not train"; the user decided to train
  since TM4 holds (findings 1.19).
* ff30_bio set up on midway2 (ff2.1 init, BioEmu library fresh in upside_input, round 3's trainer),
  update path and extract_ff.py tested with the empty glycine-offset field, submitted 07:44 as
  49179449 (estimate 10-07 03:40); panel bio_start 49179451; submit_new.sh repointed to ff30_bio
  (tags bEE; backup .bak_pre_bio_20261005).
* Library installed as `parameters/common/rama31.dat` (AWH version backed up); up.md 2.8, the GLY
  note, training/README.md, ConDiv.py docstring, GLY_sym.md and architecture.md updated. Session
  scripts, data, logs and probe files copied to `/project/trsosnic/yinhan/checks/gly_bioemu_map/`
  for use from another computer (remote_jobs.md "Round 4 ... from another computer").
* Run 2 (user): ff30_gdepth, glycine's alpha_R / alpha_L depths trainable on ff2.1's own library,
  one pooled pair on the 37 GLY|X maps (user's choice), started at BioEmu's weights (c_aR -0.0918,
  c_aL +0.4621), so no outside data enter the training. Trainer: round 3's ConDiv with round 2's
  offset hooks ported; module rewritten for the pooled pair. Tested on midway2 (round-0 library,
  one end-of-round update, extraction). Fixed on the way: extract_ff.py imported the wrong trainer
  for initial checkpoints (findings 10.14); check_step.py would have crashed on an empty glycine-
  offset field (both round-4 runs). Submitted 49179457 with start panel 49179456; submit_new.sh
  serves both runs. Trainer files and gdepth_start_offsets.py in `checks/gly_bioemu_map/`.
