# ff3.1: glycine Ramachandran handedness set by measurement

## Project Goal

ff2.1's rama library gives glycine a strong left-handed bias; ff3.0 removed it entirely by
mirror-symmetrising the whole central-GLY row. **Both are wrong, in opposite directions.** An AWH
measurement on capped `Ac-X-Gly-NHMe` dipeptides puts the true neighbour-averaged
`dG(aR->aL)` at **-0.26 E_up**, against the library's -1.32 and ff3.0's exact 0.

Build **ff3.1** = ff2.1 with the glycine handedness scaled to the measurement, retrain the core
force field against it, and test it on the 32-arm benchmark. ff3.0 is retired and is not a
reference for anything.

## Architecture and Key Decisions

**The correction is one scalar on an existing axis.** `M_new = S + lambda*A`, where `S` and `A` are
the mirror-symmetric and antisymmetric parts of the library's GLY row under
`(phi,psi) -> (-phi,-psi)`. `lambda = 1` is ff2.1, `lambda = 0` is ff3.0. `dG` is linear in lambda
(`-1.318*lambda`), so the calibration is a division, not a fit. **`lambda = 0.20`**; the two AWH
replicas imply 0.197 and 0.206, and the systematic uncertainty makes lambda good to only +/-0.07,
so more digits would be false precision.

**`GLY|GLY` is forced fully symmetric (`lambda = 0`).** A glycine flanked by glycine is achiral by
construction, so its true value is known exactly without measuring. This is physics, not fitting,
and the AWH control LG confirms it.

**Only the coil group is corrected, and that is evidence-based.** The sheet library **passes** its
own achiral control (`GLY|GLY` = +0.0011) where the coil library fails it badly (-0.7095). Sheet's
apparent handedness is a ratio of near-zeros because both helical basins are empty there. Scaling
it would be scaling noise.

**The library's own internals corroborate the measurement.** Coil `GLY|GLY` reads -0.71/-0.97 where
physics demands 0. That is an estimate of the library's systematic error requiring no simulation,
and it agrees with the AWH's independent estimate of ~1.06.

**The trainer was repaired first, and verified at ff2.1 before use.** The Theano -> PyTorch port had
silently dropped `hb` and `sheet` training; both are restored (`hb` as a multiplicative scale on
`hbond.parameters[:4]`, `sheet` as a common offset on the 20 per-type values). Restarting ConDiv
from ff2.1 gives gradients indistinguishable from minibatch noise on all four parameters, so ff2.1
is a stationary point of the objective and the port reproduces the original's optimum.

**Keep the per-neighbour structure, scaled.** The AWH could not resolve neighbour dependence
(replica-to-replica noise 0.178 against a scaled spread of ~0.11), but the library's PDB statistics
for it are solid. Scaling `A` shrinks the spread 5x along with the amplitude rather than erasing it.

## Execution Phases

### Phase 1 - measure the true glycine handedness (DONE)
- [x] 2D AWH on (phi,psi) of the central glycine, 10 capped dipeptides, GROMACS 2024.4, 300 K
- [x] Two independent replicas to their own plateaus: **-0.261** (62-100 ns) and **-0.273**
      (62-75 ns), agreeing to 0.012
- [x] Per-neighbour structure shown to be noise (Spearman rho +0.048 between replicas)

### Phase 2 - repair and verify the trainer (DONE)
- [x] Strict modernization of the Theano original; GLY palindrome and epsilon guards removed
- [x] `pack_param`'s `< 1.6e-4` residual gate restored (ff2.1 packs to 3e-30)
- [x] `hb` and `sheet` unfrozen; env param-deriv request fixed (coeff+weights = 760, slice 360)
- [x] ff2.1 confirmed a stationary point: sign-flip test p >= 0.22 on all four parameters

### Phase 3 - build ff3.1's library (DONE)
- [x] `parameters/common/rama31.dat`; only the coil GLY row differs from `rama.dat`
- [x] `GLY|GLY` = -0.000000 both directions, exactly self-mirror-symmetric
- [x] 8-neighbour dG **-0.2648**; end-to-end in a built config GLY goes -1.234 -> -0.248 while
      non-GLY is bit-identical

### Phase 4 - train ff3.1 (RUNNING)
- [x] `training/ff31/`, ff2.1 init params, `rama31.dat`, corrected trainer
- [ ] 500 minibatches, jobs 49037796 -> 49037797 -> ..., **~9 days over ~8 chain links**
- [ ] `extract_ff.py` -> `parameters/ff_3.1_trained`

### Phase 5 - benchmark (NOT STARTED)
- [ ] Re-run the 32 arms and score on the **last third**, not the whole run
- [ ] **The falsifiable prediction:** ff3.0 loses on de novo arms (-0.030) while gaining on native
      (+0.039), paired p = 0.021. If that is because zeroing the glycine alpha_L bias removed turn
      nucleation, restoring 20% should recover de novo while keeping native. If the de novo arms do
      not move, that explanation is wrong.

## Known Errors / Blockers

* **Clone training directories from `gly-ctx` or `ff31`, never from `gly-sym`.** `gly-sym` has 498
  checkpoints and looks like the best template, but its `ConDiv.py` predates the 2026-09-10
  `libupside.so` rebuild and still has the env param-deriv bug that kills every worker.
* **Sheet training costs 2.8x.** A step is ~26 min instead of ~9.4, entirely from the two extra
  passes over all 250 frames. Dropping sheet returns to the old throughput.
* **rep2's achiral control sits at +0.094 where rep1's decayed to -0.015.** Unexplained. It does
  not propagate into the neighbour-average (controls differ by 0.109, chiral averages agree to
  0.012), but it is the honest floor on quoting an uncertainty below ~0.1.
* **The site GROMACS on midway2 is unusable** (2024.1 SIGILLs, 2021.1 MPI tools SIGFPE). Use
  `/project/trsosnic/yinhan/gmx2024`.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log.
* **`/project` is ~77% full.** The retired NP campaign holds ~500 GB in `NP-1AO6/prod_ff3/` and is
  the obvious reclaim if space is needed.
