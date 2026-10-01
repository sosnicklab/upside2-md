# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`.

Condensed 2026-09-24 when the ff2.1-workflow retrain began; the superseded campaign is summarised
below and its detail is in `findings.md` 9i-9s and git history.

---

## 2026-09-18 to 09-23 — the glycine campaign (superseded)

* **Measurement (Track B):** 2D AWH on 40 capped glycine dipeptides, GROMACS, 400 ns. Finished;
  left neighbours -0.20 nats, right -0.11, context-averaged -0.154 (`rama31.dat`, rebuilt 09-24).
* **Learned map (Track A, "ff3.1"):** the central-glycine row trained as one shared `S + A` map
  with the ConDiv port of the time; 500 steps, converged to -0.885, shipped briefly as ff3.0.
  Invalidated 09-24 when that port turned out to be FF1's trainer (below). Removed from the tree.
* **Benchmark and glpG** for it were launched, then cancelled 09-24.

## 2026-09-24 — the trainer was FF1's; replaced with ff2.1's own

* Found that the port trained FF1's functional form: spline burial where ff2.1 uses a sigmoid
  (-12.9 vs -47.5 on native ubiquitin), no backbone desolvation term, no unfolded-state objective,
  burial weights never trained (findings 9t, 9u).
* Located the only FF2 dual-target trainer, Kleinmann's Python 3 port of Peng's code, and adapted
  it into `training/ConDiv.py`: torch for Theano, numpy for mdtraj, protocol restored to the SI
  (lambda 0.3, 12 replicas 0.8-1.1, 8000 time units, 24 per minibatch), replica reweighting and a
  DSE-dropping guard fixed (findings 9v). Glycine row rewritten as 42 separate maps, GLY|GLY
  mirror-symmetric, behind `TRAIN_GLY`; gate `training/verify_gly_gradient.py` passes.
* New `extract_ff.py`, `patch_glpg.py`, `validate_ff.sh`; `check_converged.py` and
  `train_chain.sbatch` rewritten (the chain now submits validation itself). `bench_run.py` type-0
  override removed. The local Mac binary traps at exit with MC moves on; worker tests run on
  midway2.
* **Phase 1, ff2.1 fixed-point check** (19 + 5 steps from ff2.1): every trained file updates; 8 of
  9 groups at a fixed point; dhb's residual pull is the small leftover of the native and 0.3 x
  unfolded gradients nearly cancelling, the balance ff2.1's own training leaves.
* **Phase 2 queued:** ff3.0 from ff2.1, 76 steps, glycine row trained, auto-release and validation.
* Cleanup for the push: the invalid `parameters/ff_3.0` and `common/rama3.dat` removed (backups in
  `backup/`), the old trainer's scripts and `init.sh` removed, the unused `_ALL` sheet pair removed
  from `upside_config.py` (`write_rama_map_pot` identical to master again), `verify_gly_gradient.py`
  tracked, stale `CLAUDE.md`, `up.md`, `GLY_sym.md`, `architecture.md` passages corrected,
  `rama31.dat` rebuilt from the finished AWH data with both Gly-Gly blanks excluded.

## 2026-09-27 — Ramachandran map redesign (design only, no code changed)

* Diagnosed why the map fails: NDRD maps are loop-site statistics of folded proteins, so they carry
  evolutionary placement; Upside adds them as a pure energy to every residue. Measured: GLY
  ln(aR/aL) -1.19 in NDRD against -0.58 over all residues of the 456 natives (findings 1.9).
* Rejected, each for a measured or verified reason: a shared correction transferred across maps
  (violates per-pair independence, findings 1.8); a per-map reference ratio from Upside's own runs
  (that is ConDiv's native term); experimental dipeptide J-couplings (blind to aR vs pPII); peptide
  all-atom maps (cost); removing the map (no remaining term has local sterics or proline's phi);
  an AlphaFold-built map (same placement, more data).
* Chosen, with the user: per-pair basin offsets on the NDRD maps, updated once per training round
  from native-restrained against free basin populations. Written into `plan.md` with three
  decisions left to confirm.
* Docs corrected: NDRD's 44,112 is the whole loop set, not glycines (`GLY_sym.md`, `findings.md`);
  the memory note's mirror index and unreproduced reference-state claim.

## 2026-09-28 — basin offsets implemented, deployed, training started

* Implemented `training/rama_basin.py` (basins partitioning the torus, 13 deg edges, per-map
  renormalisation, so each offset is a weight on its basin; 840 maps, 4,234 offsets; GLY|GLY exact)
  and `training/verify_rama_basin.py`; rewired `ConDiv.py`, `extract_ff.py`, `convergence_gate.py`,
  `gate_or_continue.sh` (stops for review), `train_chain.sbatch` (`slurm.args`), `env.sh` (per-cluster
  Python); removed `rama_gly_gradient.py`, `verify_gly_gradient.py`. Docs: README, up.md,
  architecture.md, .gitignore.
* Tests: reach gate PASS locally and on midway3; synthetic round (glycine offsets move the right
  way, untouched maps stay 0, gate flags the pull); within-basin spread of the energy change 0.03
  nats (median) at offsets of +-1. User corrections folded in: `other` basin, weight-factor design.
* Cluster: ff30 glycine chain and its monitor cancelled at step 223; files deployed md5-verified with
  backups; torch added to the shared /beagle3 venv (slow: /beagle3 degraded); midway3 attempt
  cancelled after 10 min pending; running on midway2. First link lost 4 workers to a Slurm launch
  race; the trainer now relaunches never-started workers; restarted cleanly.
* 09-28 morning: status check found round 1's log-ratio step moving empty basins by up to 1.76 nats;
  replaced by the MAP Newton step with a Gaussian prior (largest 0.44 on the same data), rewound to
  step 19 and resumed (49126332). Gate restored to release and validate automatically, BP fix re-armed
  before release; release path dry-run passed; stray macOS `._*` files removed from the cluster tree.
* 09-28 midday, user's GLY check (findings 1.11): coil GLY|GLY exact, but the engine map of a glycine
  between glycines was not, through the raw NDRD sheet entry; GLY|X deepened alpha_R vs NDRD in 30 of
  38 maps (mean +0.06), all still alpha_L-favoured. Fixed: `rama_basin.py` symmetrises the sheet
  GLY|GLY entry and uses the probability mean for coil and sheet; `extract_ff.py` and
  `verify_rama_basin.py` gate on both; "X|GLY" label corrected to GLY|X. Tests: local (all other
  entries bitwise unchanged, engine map exact), verifier PASS on the run, three GGG/terminal-GG
  configs built through `upside_config` exact. Stopped the chain during step 29, set steps 19-29
  aside, rewrote the round-1 library, resumed from step 19 (49126777). Docs: README, up.md, plan.
* 09-28 afternoon: audit of what is trained and what 456 proteins resolve (findings 1.12): per-pair
  offsets noise-limited, ~70% of the noise from the protein set. Literature survey (sub-agent) plus
  ff2.1's own misses by residue type and neighbour (findings 1.13): pre-proline alpha_R is the largest
  miss, made by the left/right mixture; glycine next; others small. Deployed the 158-offset set
  (GLY|X, GLY|GLY, X|right|PRO; the rest NDRD): code and verifier rewritten, local tests, verifier
  PASS on the run, round 1 recomputed from epoch 0, resumed from step 19 (49127867). Literature
  recorded with read status in findings 1.13.

## 2026-09-30: ff_3.0 released; pre-proline offsets diagnosed

* Status: `ff30_basin` converged at the step-114 gate (02:25), BP fix installed in both trees, ff_3.0
  released, 32 Peng arms and 4 glpG chains running. Three glpG chains had `block_count` advanced by
  Slurm requeues after NODE_FAILs; reset proposed, not done (remote_jobs.md §1).
* Pre-proline (findings 1.14): the mixture caps the X|right|PRO offsets (map aR 0.174 -> 0.144 at
  best); six rounds left the free-native gap at +0.023 while the offsets grew ~0.08 per round on
  average (max 1.89); the gap is composition (89% extended sites), not a pre-proline defect.
  Scripts `checks/prepro_leverage.py`, `prepro_control.py` (new), `prepro_residual.py`,
  `prepro_rules.py`, `prepro_left.py` rerun; logs `*_20260930.log`.
* Files: `findings.md` (1.13 corrected, 1.14 new), `plan.md` (Phase 4 done, Phase 5 running, Phase 6
  pending), `remote_jobs.md` (status, job table, ff30_basin record, BP-fix deployment).
* Validation fixes (13:20): glpG `run_remd.py` counts a block per job id (a Slurm requeue no longer
  advances `block_count`), three counters reset 3/4/3 -> 1, `submit_remd.sh` given the training node
  exclusions; backups `*.bak_pre_requeue_20260930`. Checks, all logged in `checks/*_20260930.log`:
  aborted-attempt fragments finite; TM1/TM4 on the valid windows 0.90-1.00 at T 0.70 except
  79ALA_S115T TM4 0.82-0.86; stretched C-N at T 0.88-0.90 predates ff_3.0 (pre-ff3 campaign 1.8-5.1%
  of frames, ff_3.0 1.5-2.7%); all 32 Peng arms clean. Handbook §5 rewritten (sign, windows, glycine
  criterion; `check_seeds_current.py` invalid for ff_3.0). `hdx_postfix/` and `hdx_10k/` not on disk.

## 2026-09-30 evening: glycine probe settled, ff_3.0 validation cancelled, replacement proposed

* Glycine probe (Phase 7) stopped by the user at 15 of 19 steps: from equal depth the data pull back
  toward alpha_L, -0.035 [-0.048, -0.022], 35 of 38 maps (findings 1.15). Then, also at the user's
  request, all ff_3.0 validation was cancelled (32 Peng arms, 4 glpG chains; glpG `STOP` files placed,
  no resubmission). Nothing of ours is running on midway2.
* Measurements for the replacement design (findings 1.16; scripts and logs in
  `/project/trsosnic/yinhan/checks/`, `gly_native_hbond.py`, `prepro_rightonly_natives.py`):
  * native glycines are H-bonded in both basins (alpha_R 0.83, alpha_L 0.74), in different
    `hbond_energy` branches;
  * the AWH Ac-Gly-Pro-NHMe surface has alpha_R 0.006, which `rama31.dat`'s pooled surface erases;
  * `rama_map_pot_ref` reshapes a glycine map (helical 0.105 -> 0.145);
  * right-only for pre-proline residues can be written into the library's weights (local test
    exact to 1e-5, other residues bitwise unchanged), at +0.84 E_up median for the 6.7% of
    pre-proline residues that are natively helical.
* Literature survey (sub-agent) recorded in findings 1.16.
* Files: `findings.md` (1.15 result, 1.16 new), `plan.md` (Phases 5 and 7 closed), `remote_jobs.md`
  (status, empty job table, STOP files).

## 2026-10-01: coverage of the glycine H-bond term, measured

* User concern: a glycine H-bond term reaches only H-bonded glycines. Joined each native glycine's
  H-bond state with its free/restrained basin populations, against non-glycine controls
  (`checks/gly_mismatch_by_hbond.py`, log `gly_mismatch_by_hbond_20261001.log`): of the
  glycine-specific helical loss, 62-66% is in glycines with their own H-bond, 25-32% in glycines
  spanned by a short-range bond, 2-9% in neither; natively left-handed glycines show no
  glycine-specific deficit. Every helical glycine in glpG's seed, TM4's included, has its own bond.
  Glycine already has a trained side-chain bead (`GLY_0`). Recorded in findings 1.16.
* Per-type check for a trained selection correction (`checks/type_mismatch_by_context.py`, log
  `type_mismatch_by_context_20261001.log`): per-type misses in the same native basin and H-bond state
  follow helix and beta propensity (GLY, SER, ASN fray most in helices; GLU, ALA, LEU least; VAL, ILE
  hold strands), so a trained per-type correction would flatten real propensity (findings 1.16).
* Residue counts per basin as a selection correction (`checks/aa_basin_counts.py`, log
  `aa_basin_counts_20261001.log`): glycine fills 54% of alpha_L sites, but the glycine handedness the
  counts imply depends on the reference residue (-2.61 to -1.16) and on burial (-1.95 exposed to
  -0.17 buried), so counts cannot correct the map (findings 1.16).
* User corrections: glycine and pre-proline are separate, and pre-proline is not an established
  problem. Phase 6 parked, lesson findings 10.0d, two memory notes. Experimental data checked
  (findings 1.16): the GGG populations and J-couplings exist (Andrews 2020 SI). Our AWH Gly-Gly
  surface vs experiment: pPII 0.32 vs 0.46, alpha 0.09 vs 0.06; in the matched system ff14SB is
  within 0.06, so no change. No handedness data exists. The Pace & Scholtz scale is only an
  11-system average, so the helix benchmark is dropped. Revised glycine plan written as Phase 8
  (proposed).
* Proline question (`checks/pro_mismatch_by_context.py`, log `pro_mismatch_by_context_20261001.log`):
  central prolines hold their basins better than any residue; residues before a proline, against the
  same type elsewhere, lose +0.024 (extended) and +0.059 (helical, which right-only would worsen). No
  proline fix proposed; Phase 6 stays parked (findings 1.16).

## 2026-10-01: Phase 8 implemented up to training; training blocked by the midway2 allocation

* Library: `training/build_gly_library.py` (new, tracked). `parameters/common/rama31.dat`
  rebuilt from the AWH surfaces (old build in `backup/`). Glycine row:
  * GLY|GLY: LG+RG symmetrised;
  * GLY|right|PRO: RP;
  * other GLY|X: 37 pooled;
  * coil = sheet; reference correction subtracted.

  The builder's checks pass: non-glycine identical, glycine = measured to 1e-6, G-G-G symmetric.
  End-to-end `.up` check passes. The midway2 build is byte-identical.
* Trainer: offset training removed from `ConDiv.py`, `extract_ff.py` and `convergence_gate.py`;
  `rama_basin.py` reduced to basin populations; `verify_rama_basin.py` deleted. README, up.md,
  architecture.md, GLY_sym.md and the .gitignore note updated. Deployed to midway2 (md5 verified,
  backups in `$P/backup/training_pre_glyfix_20261001`).
* Run `training/ff30_gly` set up and initialised from ff2.1 with the new library. The ff_3.0 glpG
  chains were moved aside so a release cannot delete them.
* Blocked: `pi-trsosnic` has no midway2 allocation since 10-01, so sbatch fails with
  `AssocMaxCpuPerJobLimit`. A midway3 connection attempt timed out on Duo.
* Training was started on midway3 `amd` (59834233) after both clusters refused for lack of
  allocation. The user cancelled it after 14 min (~76 core-hours): no CPU jobs on the GPU
  allocation (lesson findings 10.0e, memory note). The run directory has been restored to the
  midway2 flags and waits for the allocation. A Slack request to Tobin was drafted, with usage
  (Jul-Sep about 450k core-hours, September about 280k).
* Local run: a full worker ran cleanly on the Mac (1ga3, 345 s, 1.6 GB; the old MC SIGTRAP is gone).
  Added `CONDIV_LOCAL_WORKERS` to ConDiv's launcher (no Slurm: at most N at once) and two scripts:
  `training/move_run.py` (checkpoint path remap for the transfer, tested both ways) and
  `training/check_step.py` (per-step physical health; on ff_3.0's epoch 5 it reproduces findings
  1.15). Deployed to midway2. `training/ff30_gly_local` started 09:16 with two workers under
  caffeinate. Session cron `4b09fd3d` checks the allocation and transfers the run when it returns.
