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
