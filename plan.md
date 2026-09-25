# ff3.0: ff2.1's own training workflow, then the glycine row

## Project Goal

Train ff3.0 from ff2.1 with **exactly ff2.1's training workflow**, modernised only, plus two project
additions: the GLY|GLY maps held mirror-symmetric and the central glycine row of the Ramachandran
library trained. The result is released as `parameters/ff_3.0`, which is absent from the tree
until then: both earlier versions were trained with FF1's workflow and were removed.

## Architecture and Key Decisions

**Why the trainer was replaced (findings 9t, 9u).** The previous port reproduced Peng's FF1 Theano
trainer, not FF2's. It trained the spline burial table (type 0) where ff2.1 uses the sigmoid (type
1; -12.9 vs -47.5 on native ubiquitin), never built the backbone desolvation node, had no
unfolded-state objective, and never trained the 400 burial weights.

**Source: the only FF2 dual-target trainer that exists.** O. Kleinmann's Python 3 port of Peng's
code, `/project2/trsosnic/okleinmann/condiv/condiv2.py`, adapted into `training/ConDiv.py`. Every
difference from it is listed and justified in that file's docstring: Theano -> torch and mdtraj ->
numpy (modernisation); lambda, replicas, length and minibatch restored to the SI (the port had
drifted and ran with lambda = 0); the replica reweighting exponent fixed (it was T0 times the
correct one); a guard that silently dropped the DSE term removed.

**Protocol (Peng et al. 2022 SI).** Per protein and step: one native-restrained replica, 12 free
replicas at T = 0.8 to 1.1, one SARW replica; 8000 time units, second half analysed; contrast =
NSE + 0.3 * DSE; 456 proteins in 19 minibatches of 24; 4 epochs = 76 steps; 24 workers x 14 cores.

**Exactly ff2.1's parameter set, by the user's choice (2026-09-24).** Trained: rot (pair,
coverage, hydrophobe), the 20 x 3 sigmoid burial parameters and 400 weights, the backbone term's
scale, the three H-bond energies and the second-H-bond term, the 20 sheet values. Not trained,
exactly as in ff2.1's training, because the engine returns no derivative: the backbone term's
center, sharpness and hbond weight (commented out in master's `get_param_deriv`), and `hbond.h5`
entries 4-11. Non-glycine Ramachandran rows stay the library.

**Glycine row, phase 2 only (`TRAIN_GLY`).** All 42 finite maps of the central-GLY coil row are
separate parameters starting from ff2.1's; the two GLY|GLY maps are projected mirror-symmetric at
the start and after every update; updates are Fourier band-limited per map. The analytic gradient
is gated by `training/verify_gly_gradient.py` (passed 2026-09-24 on 1bgf: spline 1e-6, map
reconstruction 3e-6, directional finite differences 6e-5 to 1.3e-3).

**Revised 2026-09-25: glycine step size doubled from step 77, 0.005 -> 0.01 (after the global
factor).** The glycine alpha is the project's own choice, not part of ff2.1's workflow, and it was
the limit: over steps 58-76 each cell moved in a consistent direction (sign consistency 0.70
against 0.23 for noise) while Adam's per-step utilisation sat near the noise floor (0.37 against
0.33), the signature of a weak steady pull whose drift scales with alpha. Only `solver.alpha.gly`
in the step-77 checkpoint was changed (backup `checkpoint.pkl.bak_gly_alpha_0.005`); the nine ff2.1
groups keep their rates and had passed the gate. The destination is unchanged; only the rate.

## Execution Phases

### Phase 1 - validate the trainer on ff2.1 (DONE 2026-09-24)
- [x] Adapt the port; switchable glycine row; `check_converged.py`, `train_chain.sbatch`,
      `extract_ff.py`, `patch_glpg.py`, `validate_ff.sh` written for it
- [x] Glycine gradient gate passes, locally and on midway2
- [x] Smoke worker and full-length timing on midway2: ~21 min/step (5vhg, 150 res, 1242 s)
- [x] One epoch (19 steps) from ff2.1: every trained file updates; 8 of 9 groups at a fixed point
- [x] The ninth, dhb (second-H-bond term), diagnosed over 5 more steps: the native and 0.3 x
      unfolded gradients nearly cancel at ff2.1 (the balance lambda = 0.3 training leaves), and
      the reweighting fix reduces rather than causes the residual. A mildly unconverged ff2.1
      parameter, not a port error (findings 9v)

### Phase 2 - train ff3.0 from ff2.1 (QUEUED 2026-09-24 ~11:55)
- [x] `TRAIN_GLY = True`; `training/ff30` initialised from ff2.1; gate re-run on midway2
- [x] `bench_run.py`'s type-0 override and `rama3.dat` fallback removed
- [ ] 76 steps done 2026-09-25 (49074120, 26:49); gate at 76: all groups ok except `gly` (p = 0).
      Epoch 5 was cancelled after step 77 to double the glycine alpha (above) and resumes from 77
      as 49119234, target 95. The convergence gate
      (`convergence_gate.py`, exact sign-flip test per group over the last epoch, family-wise 5%):
      converged -> `validate_ff.sh ff_3.0` releases to both trees and starts the 32 Peng arms and 4
      glpG chains; not converged -> one more epoch (~7 h) and judge again, up to 13 epochs, then
      stop for review. Calibrated on ff2.1: all groups p 0.74-1.0 except dhb, p = 0.002
- [x] Unattended path audited and dry-run 2026-09-25 (remote_jobs.md); `train_chain.sbatch` stops
      after three links in a row fail at the same step
- [ ] Copy the released `parameters/ff_3.0` into the local repo

### Phase 3 - validation (NOT STARTED)
- [ ] Peng benchmark, 16 proteins x native/de novo, scored on the last third, paired against ff2.1
- [ ] glpG, four variants: helix stability over time, TM4 above all

## Known Errors / Blockers

* **The local Mac `obj/upside` traps (SIGTRAP, exit 133) at exit whenever Monte Carlo pivot moves
  are on**, even for one system, after every frame completes. The midway2 binary ran
  `mc_interval = 5` for 500 steps cleanly, so trainer tests run on midway2. Not yet diagnosed; the
  Mac binary dates from 2026-08-24.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
* **`/project` has ~445 G free.** The glpG REMD trees hold ~1.26 T; `NP-1AO6` ~0.5 T is the
  obvious reclaim.
* **Quote ~20% uncertainty on the AWH glycine handedness** (findings 9s); the like-for-like target
  is -0.15 to -0.31 depending on whether neighbour effects add.
