# ff3.1: the glycine Ramachandran map, learned and measured

## Project Goal

ff2.1's rama library gives glycine a strong left-handed bias (`dG(aR->aL) = -1.24 nats`); ff3.0
removed it entirely by mirror-symmetrising the row. Both are wrong. **Two independent routes to the
right answer are now running in parallel, and they check each other:**

| | route | what it produces |
|---|---|---|
| **Track A** | ConDiv learns the glycine map from the 456-protein training set | a map consistent with Upside's own force field and the native ensemble |
| **Track B** | 2D AWH on all 40 capped `X-Gly` / `Gly-X` dipeptides | an independent physical map from explicit-solvent MD |

Neither is imported into the other. When Track A converges, its map is compared against Track B's
to decide whether further training is needed. Agreement is evidence; disagreement localises which
term is wrong.

Both routes replaced an earlier plan to set the map from Track B alone. That map was built and is
kept as `parameters/common/rama31.dat`, but training on it was stopped: a map imported from
dipeptide MD is still an assumption about protein interiors, and it cannot be validated against
the thing it was imported into.

## Architecture and Key Decisions

### Track A: the glycine map as a trained parameter

**Parameterisation. Two 72x72 maps, `S` symmetric and `A` antisymmetric under
`(phi,psi) -> (-phi,-psi)`:**

* `X|GLY`, every neighbour and both directions: `S + A`, unconstrained in total
* `GLY|GLY`: `S` alone

**One map serves every neighbour, and that is a real loss taken deliberately.** The library's
glycine row resolves neighbour dependence (spread 0.264 nats) and left/right asymmetry (0.679
nats); this parameterisation discards both at step 0. Keeping them would mean 40 x 5,184
parameters against ~30,000 glycine samples per minibatch, and Track B cannot resolve per-pair
structure either (S/N 1.48 against 3.89 for the average). Whether one map is enough is precisely
what the 40-context AWH campaign is for.

That *is* the constraint the physics demands and nothing more. A glycine flanked by L-amino acids
sits in a chiral environment and may be biased; a glycine flanked by glycines has no chirality
source, so its map must be mirror-symmetric. Writing the pair as `S` and `A` makes `GLY|GLY = S`
symmetric by construction, with no projection step during training.

**Starting point: `A = 0`, `S` = the library's symmetrised coil GLY row.** So training begins
assuming no handedness at all and has to discover it from protein data. **This is what makes the
Track B comparison a real test**: if the learned `A` lands near the measured one, two methods
sharing no input agree. Seeding from the measurement would make the comparison circular. The sheet
group is untouched, as before.

**Why this fixes what the library gets wrong.** The library is `-lnP` over glycines in folded
structures, so it carries `E_local + E_fold`, and Upside then adds its own `hbond + env +
sidechain` model of the fold: the fold is counted twice. Contrastive divergence fits the local term
that reproduces the native ensemble **given the rest of the force field**, which subtracts exactly
the double-counted part. That is the defect, addressed at its cause rather than corrected for.

**The gradient is exact and needs no C++ change.** Established by measurement:

* `rama_map_pot` evaluates a **periodic interpolating bicubic spline**, and
  `solve_periodic_2d_spline` is a tensor product of 1D periodic solves
  (`src/spline.cpp:262`). So the map-to-energy operator is **linear, separable and
  translation-invariant**: `E = sum_r sum_ij map_r[i,j] b(x_r - i) b(y_r - j)` with one 1D
  cardinal function `b`.
* Therefore `dE/d(map_r[i,j]) = sum_frames b(x_r - i) b(y_r - j)`, a smoothed 2D histogram of the
  glycine `(phi,psi)` samples. No finite differences, so no per-parameter cost: finite differencing
  5,184 map values would cost 10,369x a divergence, which is why `sheet` is a single scalar.
* **Verified against the engine**: reconstructing `rama_map_pot` in Python from `rama_coord` plus
  the config's `rama_pot` arrays gives -5.332561 against the engine's -5.332566, a 5e-6 float32
  round-off difference.

**Chain rule from the per-residue config map back to `S` and `A`.** `read_weighted_maps` mixes the
left and right neighbour maps and then the coil and sheet maps, both as weighted log-sum-exp, so
`dC_r/dS` and `dC_r/dA` are softmax factors. Two simplifications make this cheap: glycine's
`dimer_weight` entries are all exactly 1.0, and all `X|GLY` maps are identical, so for a glycine
with no glycine neighbour the inner mixture returns `S + A` exactly and only the coil-versus-sheet
factor survives.

**Smoothing.** The exact gradient is already spread over a 4x4 node neighbourhood by the spline
basis, but 5,184 free values against ~30,000 glycine samples per minibatch will still be noisy.
The gradient is projected onto a truncated 2D Fourier basis on the torus before the Adam step, so
the learned correction is smooth by construction and the truncation order is a knob rather than a
rewrite. The starting map keeps its sharp forbidden-region structure, since only the correction is
band-limited.

### Track B: complete the dipeptide set

All 40 contexts, not the 10 that were run. `awh_batch.sbatch` builds and runs any group of systems
and self-chains to a target time; the 40 PDBs already existed in `gly_peptides/xg/`.

**Why the other 30 are needed.** The existing map asserts no neighbour dependence and no left/right
asymmetry, because only 10 left-neighbour systems were ever run and per-pair S/N was 1.48. The
library resolves both effects (0.264 and 0.679 nats). Running all 40 either resolves them or rules
them out.

**The achiral `GLY|GLY` blank is the convergence criterion, not a wall time.** It must read 0,
reads rms 0.032 at 100 ns, and falls as `1/sqrt(t)` while the signal converges: 400 ns should
halve it.

**Units, settled.** A library map holds `-lnP` normalised so `sum(exp(-E)) = 1`, exact on every map
in `rama.dat`. A PMF in kJ/mol therefore enters as `PMF / kT(300 K) = PMF / 2.494339`, **not**
`PMF / 2.914952774272`. The two differ by 1.169x.

## Execution Phases

### Phase 1 - measure the glycine handedness (Track B, in progress)
- [x] 2D AWH on (phi,psi) of the central glycine, 10 `Ac-X-Gly-NHMe` dipeptides, GROMACS 2024.4
- [x] Two independent replicas agreeing to 0.012; per-neighbour structure not resolved
      (Spearman rho +0.048)
- [x] Achiral blank validated as sampling noise: uncorrelated between replicas (r = +0.157) and
      decaying as `1/sqrt(t)` (0.233 -> 0.032) while the signal converges (0.080 -> 0.071)
- [ ] **Extend the original 10 to 400 ns** (49037819, 49037820)
- [ ] **Run the missing 30 contexts to 400 ns**: 10 left neighbours (49037909) and all 20 right
      neighbours (49037910, 49037911), self-chaining
- [ ] Rebuild the reference map with `py/build_rama_from_awh.py` and report per-neighbour and
      left/right structure with error bars

### Phase 2 - repair and verify the trainer (DONE)
- [x] Strict modernization of the Theano original; GLY palindrome and epsilon guards removed
- [x] `pack_param`'s `< 1.6e-4` residual gate restored (ff2.1 packs to 3e-30)
- [x] `hb` and `sheet` unfrozen; env param-deriv request fixed (coeff+weights = 760, slice 360)
- [x] ff2.1 confirmed a stationary point: sign-flip test p >= 0.22 on all four parameters

### Phase 3 - make the glycine map trainable (Track A, STARTING)
- [x] Establish the gradient is analytic: linear separable spline, verified to 5e-6 against the
      engine
- [x] `training/rama_gly_gradient.py`: cardinal-function weights, per-residue chain rule to `S` and `A`
      by torch autograd, Fourier projection, library writer, symmetric start
- [x] **Gradient verified against finite differences through the whole pipeline**
      (`training/verify_gly_gradient.py`). It failed first at 37%, which found the missing
      Boltzmann shift in `write_rama_map_pot`; it now agrees to **3.8e-5** at the sweep optimum,
      with the rise at small eps explained by float32 library storage
- [x] **Terminal glycines added to the gradient (2026-09-19).** They were skipped, but a terminal
      glycine still gets a glycine map, so 3.0% of the glycine gradient was missing and biased
      toward flexible chain ends. The gate missed it because 1a62 has no terminal glycine; the
      verifier now reports which branches a test protein exercises and is run on proteins chosen
      to cover terminal and glycine-adjacent cases
- [x] Symmetric starting map: `rama_gly_gradient.symmetric_start`, coil GLY row symmetrised and
      renormalised, `A = 0`, sheet untouched
- [x] `training/ConDiv.py`: `gly` added to `Update`, accumulated per frame in
      `compute_divergence`, library written in `expand_param` and deleted once the workers exit,
      update band-limited and re-projected in `backprop_deriv`, `initial_alpha.gly = 0.02`
- [x] `training/ff31-gly/` set up from the **ff2.1 original** `rama.dat`, not the measured map.
      `initialize` gives `pack_param residual = 3.38e-30` unchanged, `dG(aR->aL) = +0.0000`,
      `|A| = 0`, `GLY|GLY asymmetry = 0`
- [x] Two-minibatch smoke test (49037934) passed: gradient finite, map moves ~`alpha` per step,
      `GLY|GLY` asymmetry `0.00e+00` throughout, and **686-698 s/minibatch against 650-712 with
      the map fixed, so training it is free**
- [x] `initial_alpha.gly = 0.02` kept: it puts the measured -0.303 about 60 steps away and the
      library's -1.24 about 240, so 500 steps can reach either and then sit

### Phase 4 - train (RUNNING)
- [ ] 500 minibatches from ff_2.1 with the symmetric starting map, ~700 s/step so ~4 days over
      ~4 chain links. **Restarted from step 0 after the terminal-glycine fix**; the first attempt
      (49037939) reached step 24 under the incomplete gradient and was discarded. Watch `gly dG(aR->aL)` against the
      AWH's -0.303, and `GLY|GLY asymmetry` which must stay 0
- [ ] `extract_ff.py` -> `parameters/ff_3.1_trained`

### Phase 5 - compare the two tracks, then benchmark (NOT STARTED)
- [ ] Learned `A` against measured `A`: basin dG, full-surface correlation, per-neighbour structure
- [ ] Re-run the 32 benchmark arms and score on the **last third**, not the whole run
- [ ] **The falsifiable prediction:** ff3.0 loses on de novo arms (-0.030) while gaining on native
      (+0.039), paired p = 0.021. If that is because zeroing the glycine alpha_L bias removed turn
      nucleation, restoring the measured part of it should recover de novo while keeping native. If
      the de novo arms do not move, that explanation is wrong.

## Known Errors / Blockers

* **Clone training directories from `gly-ctx` or `ff31`, never from `gly-sym`.** `gly-sym` has 498
  checkpoints and looks like the best template, but its `ConDiv.py` predates the 2026-09-10
  `libupside.so` rebuild and still has the env param-deriv bug that kills every worker.
* **Sheet training is effectively free**, contrary to an earlier note here. Measured
  625-705 s/minibatch with sheet trained against 650-712 historically with it frozen. A step is
  ~11 min, so 500 steps is ~3.8 days.
* **`parameters/common/rama31.dat` is NOT the training map any more.** It is the Track B
  reference, built from rep1 at 100 ns and rep2 at ~76 ns, and it will be rebuilt when the
  extensions land. Track A trains from a symmetric starting map instead, so that the comparison
  between the two is not circular.
* **Track A's gradient must be verified against finite differences before any long run.** The
  analytic route is only valid because the map enters linearly; a mistake in the chain rule
  through the coil/sheet mixture would train silently in a wrong direction for days. The
  engine-reconstruction check (5e-6) validates the spline half, not the mixture half.
* **5,184 free map values is a lot against ~30,000 glycine samples per minibatch.** Fourier
  truncation is the first line of defence. If the learned map still comes out noisy, lower the
  truncation order rather than adding a penalty term on top.
* **Quote ~20% uncertainty on the handedness, not +/-0.01.** The replica agreement of 0.012 is not
  the error bar. The achiral blank still carries a +0.039 basin asymmetry against a signal of
  -0.194, and a single surface is only S/N 2.2 at 100 ns. The 400 ns extension should roughly
  halve this.
* **The site GROMACS on midway2 is unusable** (2024.1 SIGILLs, 2021.1 MPI tools SIGFPE). Use
  `/project/trsosnic/yinhan/gmx2024`.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log.
* **`/project` is ~77% full.** The retired NP campaign holds ~500 GB in `NP-1AO6/prod_ff3/` and is
  the obvious reclaim if space is needed.
