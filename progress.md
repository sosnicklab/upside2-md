# Progress log

High-level execution diary. Job ids, states and log paths live in `remote_jobs.md` only;
technical findings live in `findings.md`; technical direction lives in `plan.md`.

Condensed 2026-09-19 when the ff3.1 phase began. Earlier day-by-day detail was removed where
`findings.md` already carries the record.

---

## 2026-09-18/19 — the glycine question answered, ff3.1 built and training

**Measured the true glycine handedness.** 2D AWH on (phi,psi) of the central glycine in ten capped
`Ac-X-Gly-NHMe` dipeptides. Two independent replicas reached their own plateaus and agree:
**-0.261** (rep1, 62-100 ns) and **-0.273** (rep2, 62-75 ns), i.e. **-0.26 E_up** against the
library's -1.32 and ff3.0's exact 0. Per-neighbour structure is noise (Spearman rho +0.048 between
replicas), so the deliverable is one number, not a 40-entry row.

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

**Built and launched ff3.1.** `rama31.dat` = ff2.1 with the coil GLY antisymmetry scaled to 0.20
and `GLY|GLY` forced symmetric. Verified at library, config and trainer level. Training chain
49037796 -> 49037797 -> ..., 500 steps, ~9 days.

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

* **rep2's achiral control sits at +0.094 where rep1's decayed to -0.015.** Unexplained. It does not
  propagate into the neighbour-average, but it is the floor on quoting an uncertainty below ~0.1.
* **The PI email and both figures in `~/Downloads`** are marked `UNCONVERGED_DO_NOT_USE` / `HOLD`.
  They argued the handedness is zero, which is now known false; they need rewriting from the
  -0.26 result.
* **glpG production is paused** on ff3.0 configs. Anything resumed there should wait for ff3.1.
