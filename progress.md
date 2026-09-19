# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`.

Condensed 2026-09-19 when the ff3.1 phase began. Earlier day-by-day detail was removed where
`findings.md` already carries the record.

---

## 2026-09-19 (later) — two independent routes to the glycine map

The user stopped training on the measured map and set a two-track plan: **learn** the glycine map
from the protein training set, **measure** it from GROMACS in parallel, then compare. Training on
an imported dipeptide map was the thing being abandoned, because a map imported from small-molecule
MD is still an assumption about protein interiors and cannot be validated against what it was
imported into.

**Track B, the measurement, is running at full scope.** All 40 neighbour contexts, not the 10 ever
run: the missing 10 left neighbours (49037918) and all 20 right neighbours (49037919, 49037920),
self-chaining to 400 ns, plus the original 10 extending (49037819/20). The 40 dipeptide PDBs
already existed in `gly_peptides/xg/`. This resolves or rules out the two structures ff3.1 asserts
away, per-neighbour dependence and left/right asymmetry, which the library puts at 0.264 and
0.679 nats.

**Track A is training** (49037939, ~4 days). The two-minibatch smoke test passed every invariant:
the gradient is finite, the map moves by about `alpha` per step, `GLY|GLY` asymmetry holds at
`0.00e+00`, and at **686-698 s/minibatch against 650-712 with the map fixed, training it is
free** - the payoff for the gradient being analytic rather than finite-differenced.

**And the first two steps already say something.** Starting from *exactly zero* handedness, the
protein data pushes the map left-handed (-0.0080 then -0.0128 nats). That is the same sign as the
library (-1.24) and the AWH measurement (-0.303): **three sources sharing no input agree that
glycine is left-handed**, and ff3.0's assertion that the answer is zero is the one value none of
them support. The magnitude means nothing at n=2, since Adam's first step is scale-invariant.

**Track A, the trained map, is implemented up to its gate.** The parameterisation is two 72x72
maps, `S` symmetric and `A` antisymmetric, with `X|GLY = S + A` and `GLY|GLY = S`, so glycine's
molecular symmetry holds exactly where it applies and nowhere else. Training starts from `A = 0`
and a symmetrised library row, so the handedness is entirely learned and the comparison against
Track B is not circular.

**The gradient is analytic, which is the only reason a 5,184-value trainable map is affordable.**
`rama_map_pot` is a periodic interpolating bicubic spline and `solve_periodic_2d_spline` is a
tensor product of 1D solves, so the map-to-energy operator is linear, separable and translation
invariant: `dE/d(map[i,j])` is a spline-smoothed 2D histogram of the glycine (phi,psi) samples.
Finite differencing instead would cost 10,369x a divergence. Reconstructing `rama_map_pot` in
Python agrees with the engine to 5e-6.

**The verification gate immediately earned its keep.** `training/verify_gly_gradient.py` checks
the analytic gradient against finite differences taken through the whole pipeline, and it failed
on the first run at 37% error. Cause: `write_rama_map_pot` ends with a per-map Boltzmann-weighted
shift, `rama_pot -= (rama_pot*exp(-rama_pot)).sum()`, which is constant per map and therefore
invisible in forces and in every basin difference, but map-dependent and so carrying gradient.
With it included the reconstructed per-residue map matches the config to 1.6e-6, the library's
float32 resolution. Documented in `up.md` 2.8a. A second smaller trap in the same area: glycine's
`dimer_weight` is **not** 1.0 (0.908 / 0.948), and printing it with `precision=0` rounds it to 1.

**Disk was going to cost ~100 GB and most of it was never read.** 645 MB per system per 100 ns, of
which 54 MB is `pullx`/`pullf` that nothing opens and most of the 96 MB `.edr` is plain energy
frames at 10x the AWH output rate. Set `nstenergy` equal to `awh-nstout` and disabled the pull
output files, both output-only with no effect on the dynamics. Measured afterwards at **0.93 MB/ns
against 1.62**, so the campaign costs ~21 GB rather than ~100. Most of what is left is the AWH 2D
grid inside the `.edr` at 90 KB per frame, so `awh-nstout` is the only lever left and it is not
worth a third rebuild. Deleted closed `part0001_pull?.xvg` for 1.2 GB back. The RCC quota tools
are both broken on midway2, so footprint is tracked with `du`.

---

## 2026-09-18/19 — the glycine question answered, ff3.1 built and training

**Measured the true glycine handedness.** 2D AWH on (phi,psi) of the central glycine in ten capped
`Ac-X-Gly-NHMe` dipeptides. Two independent replicas reached their own plateaus and agree:
**-0.305** (rep1, 62-100 ns) and **-0.319** (rep2, 62-75 ns), i.e. **-0.30 nats** against the
library's -1.24 and ff3.0's exact 0. Per-neighbour structure is not resolved: the dG orderings
correlate at Spearman rho +0.048, and even the per-neighbour 2D surfaces only reach signal-to-noise
1.48 against the averaged surface's 3.89. So the deliverable is one surface, not a 40-entry row.

**Validated the measurement against its own blank, after the challenge that a non-symmetric
`Gly-Gly` means the simulation is wrong.** It is sampling noise, shown two ways: the residual is
uncorrelated between replicas (r = +0.157, where a real bug would reproduce), and it decays as
`1/sqrt(t)` from 0.233 to 0.032 E_up while the signal converges to 0.071. An artifact would decay
too. Blank subtraction was rejected, since subtracting a noisy estimate of zero adds noise. The
honest error bar is ~20% on the basin dG, not the +/-0.01 the replica agreement suggests, so both
replicas were extended from 100 to 400 ns (49037819/20) with the blank as the stopping criterion.

**Found that the library fails its own achiral control.** Coil `GLY|GLY` reads -0.71/-0.97 where
physics demands exactly 0. An estimate of the library's systematic error needing no simulation, and
it agrees with the AWH's independent estimate of ~1.06. The sheet group, by contrast, passes
(+0.0011), which is why ff3.1 corrects coil only.

**Completed the 32-arm benchmark.** ff3.0 gains on native (+0.039) and loses on de novo (-0.030);
paired across 16 proteins, 13 of 16 positive, p = 0.021. lambda's failure is lambda-specific, not a
topology class: alpha3d is the same fold class at L=73 and folds fine from native (3.57 A).

**Repaired the trainer.** The Theano -> PyTorch port was faithful in its maths but had silently
dropped `hb` and `sheet` training and added a GLY palindrome, two epsilon guards, and a weakened
`pack_param` convergence gate. All reverted or restored. Verified ff2.1 is a stationary point of
the repaired objective (sign-flip test p >= 0.22 on all four parameters), so the port reproduces
the original's optimum.

**Built and launched ff3.1**, after four rejected constructions. Final: **the library's
central-glycine coil row is discarded and rebuilt from the AWH surfaces**, containing no library
data. Two maps: `X|GLY` = the measured surface, `GLY|GLY` = its exactly symmetric part. Nothing
fitted. `dG(aR->aL)` is -0.303 nats against the library's -1.24, and a glycine between two
glycines gets no handedness at all, which falls out of the left/right map mixture rather than
being imposed. Builder saved as `py/build_rama_from_awh.py`. Training chain 49037907 -> 49037908,
500 steps, ~3.8 days.

The route there was driven entirely by user challenges, and the last one overturned my own
argument. A uniform `S + 0.20*A` was built first, then `S + 0.4245*(A - A_GG)` when the manual
`GLY|GLY` exception showed the multiplicative model was the wrong shape, then `S_library +
A_measured` when "can't you get the map from the GROMACS runs?" proved right about the
antisymmetric part. That third one was rejected too, on the observation that **`S_library` is
exactly ff3.0**, so the construction was ff3.0 plus a patch, on a map already established as
wrong. Replacing the row outright removed the last library dependence.

**Two of my own claims had to be withdrawn to get there**, both recorded in `findings.md` and
`GLY_sym.md`:
* "The library's symmetric part must be kept; it differs from the measured one by rms 1.752 E_up,
  aR basins 0.193 against 0.089." Those were per-pair *antisymmetric* statistics, not symmetric
  parts. Measured properly the two surfaces correlate at **r = +0.867**, their aR basins agree to
  **0.207 nats**, and the map mean is unchanged to three decimals by the swap. The only large
  disagreement is the aL basin, 1.087 nats, which is the handedness error itself.
* The units choice was inverted. A library map holds `-lnP` normalised to `sum(exp(-E)) = 1`
  (verified exactly on every map), so an AWH PMF must be divided by `kT(300 K) = 2.494339`, not by
  2.914952774272 kJ/mol per E_up. Same measurement, stated in the convention the file uses.

**Checked that replacing the row is safe rather than assuming it.** Everything outside the coil
GLY row is bit-identical including the sheet group and both weight arrays; all 42 GLY maps
normalise to 1.000000; the map mean holds at 11.577 so glycine's weight against the other 19
types does not shift; the aR->aL saddle *falls* 6.32 -> 4.10 nats so glycine samples more freely;
and the engine gives a finite energy on 1a62. The measured surface's higher global maximum (25.9
against 18.3) sits in a forbidden corner no path crosses, and the library cannot represent that
region anyway, since an empty bin among 44,112 glycines is censored at about ln N = 10.7.

**Retired ff3.0** at the user's direction, and cancelled `np_1AO6_prod`, the only job still
simulating with it.

**Corrections made this session, all recorded in `findings.md`:**
* "`hb` has never been trained" — wrong. The Theano original trains it at lr 0.02; the port dropped
  it. The proposed one-line learning-rate edit would have been a no-op, since the gradient is
  hardcoded to zero in two other places and there is no writer.
* "All ConDiv training has been broken since the 2026-09-10 rebuild" — wrong. The env param-deriv
  bug was already diagnosed and fixed in `gly-ctx`, which trained 141 steps afterwards. It bit
  again only because `ff21-restart` was cloned from the un-patched `gly-sym`.
* "The sheet gradient is systematically nonzero" (t = +3.91 at n=3) — collapsed to t = -0.01 at
  n=6. Third instance this session of a partial-n statistic misleading.
* ff3.0 cannot be cited as precedent for anything, since it is the thing being replaced.

---

## 2026-09-16/17 — is ff3.0's glycine symmetrisation right?

The PI asked whether forcing every glycine map mirror-symmetric is justified, given only `GLY|GLY`
is achiral. Established that glpG has no GGG and that all three TM4 glycines are XGX, so ff3.0's
improvement acts entirely through the context where symmetry is *not* forced by any argument.
ff3.0B (GG-only) and ff3.0C (context-subtracted) were both built and both cancelled — ff3.0B's
perturbation was below the retrain reproducibility floor, ff3.0C's founding premise was
contradicted by the measurement. This is what motivated measuring the true value directly.

## 2026-09-18 — why the lambda benchmark arm fails

Every helix holds its own shape (local Ca-RMSD 0.62-2.84 A) while the assembly is 8.6 A out; the
error is in crossing angles (H0-H3 28 -> 69 deg, H1-H4 22 -> 86 deg) with centroid distances right
to within 4 A. H2 (`QSGVGALFN`) unravels first. Discrimination between folded and misfolded sits
outside the local term (+0.205 non-rama against +0.075 rama). The candidate glycine map does
nothing for lambda: -0.024 E_up on the equilibrated last third.

**Two method results that outlived the lambda question.** The native benchmark arms are not at
equilibrium — lambda's decays 6.41 -> 10.39 A and `BURN = 2000` discards only one of twelve blocks,
so every quoted native number mixes decay with equilibrium; re-score on the last third before
ranking any force field. And a reweighting is only as meaningful as the stationarity of the
ensemble it reweights.

## 2026-09-05 to 2026-09-14 — ff3.0 training, deployment and glpG delivery

ConDiv `gly-sym` retraining to 500 minibatches on midway2 and rockfish; the self-resubmitting chain
proved itself unattended through a FAILED link and a NODE_FAILURE, resuming at exactly the right
step each time. ff3.0 deployed 2026-09-09. glpG VTF trajectories delivered and re-delivered after a
periodic-image fix (unwrap first, then wrap by molecule centroid — a per-atom wrap threw the protein
a box length out of the bilayer). DDM retired as an environment; glpG runs in POPE/POPG only, which
leaves the crystal as the only fidelity reference.

## Carried-over open items

* **AWH extensions to 400 ns are running** (49037819/20). When they land, re-derive `A_measured`
  and rebuild `rama31.dat` if it shifts; read the last `awh.part*.edr`, the AWH state is cumulative.
* **The PI email and both figures in `~/Downloads`** are marked `UNCONVERGED_DO_NOT_USE` / `HOLD`.
  They argued the handedness is zero, which is now known false; they need rewriting from the
  -0.26 result.
* **glpG production is paused** on ff3.0 configs. Anything resumed there should wait for ff3.1.
