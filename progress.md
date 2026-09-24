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
