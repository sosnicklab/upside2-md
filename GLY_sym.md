# Glycine Ramachandran handedness: the problem, the measurement, and ff3.1

Rewritten 2026-09-19. The previous version of this file made three claims that our own
measurements and its own cited literature contradict; they are listed at the end under
[Corrections](#corrections-to-the-previous-version) so the error is on the record rather than
silently deleted.

---

## 1. Why glycine handedness matters, and the exact limit of the symmetry argument

Glycine has no Cβ. Swapping its two Hα substituents maps the molecule onto itself, so **an isolated
glycine residue is achiral** and its (φ,ψ) free-energy surface must be symmetric under
`(φ,ψ) → (−φ,−ψ)`. Every other residue has a Cβ that breaks this.

**But that argument only reaches as far as the molecule's own symmetry.** A glycine inside a
protein is flanked by L-amino acids. The local environment is chiral, so the *conditional*
distribution `P(φ,ψ | glycine, L-neighbours)` has no symmetry requirement at all. The symmetry is
exact only when the neighbours are themselves achiral:

| context | is symmetry required? |
|---|---|
| glycine in an achiral environment (`Gly-Gly-Gly`, poly-Gly) | **yes, exactly** |
| glycine with any L-amino-acid neighbour (`X-Gly-X`) | **no** |

This distinction is the whole substance of the problem. It is also confirmed from two independent
directions. Simulation studies of glycine oligomers find that for Gly₃/Gly₁₀ the αR and αL
populations "would approach equality" in the limit of sufficient sampling
([Force field dependent solution properties of glycine oligomers](https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/)).
And the continuous-chirality analysis of 4,366 glycine residues from 160 high-resolution
structures finds that **glycine is "practically always conformationally chiral", with chirality
magnitudes similar to genuinely chiral amino acids**
([Baruch-Shpigler et al. 2017](https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525)).

So: achiral molecule, chiral conformation, chiral neighbourhood. Forcing the map symmetric
everywhere throws away something real.

---

## 2. What is actually wrong with ff2.1's glycine row

Upside's local term is the neighbour-dependent Ramachandran library **NDRD_TCB** of
[Ting, Wang, Shapovalov, Mitra, Jordan & Dunbrack, *PLoS Comput Biol* **6**(4): e1000763 (2010)](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000763).

Measured on our own basins (`αR: φ∈[−100,−40], ψ∈[−60,10]`; `αL: φ∈[40,100], ψ∈[−10,60]`; exact
mirror images, 195 grid cells each), the library's central-glycine row gives a neighbour-averaged
`dG(αR→αL)` of **−1.24 nats** — a large left-handed preference. (Quoted as −1.32 in earlier
text, which used the E_up energy conversion; see the units discussion in §5.)

**It is not sampling noise.** This is the previous version's central error. The library carries
44,112 glycine residues; binomial noise on `P(φ>0)` is 0.0023 while the library sits 0.155 away
from 0.5, about **68σ**. Per-neighbour entries at ~1,500 residues have σ = 0.012 against an
observed neighbour spread of ±0.10. The asymmetry is a real, well-resolved feature of the PDB.

**It is not wrong about the PDB either.** Measured independently from the 16 ff3.0 benchmark
native structures, **42 of 67 interior glycines sit at φ > 0 (63%)**, against the library's 66%.

**The defect is an ensemble mismatch that double-counts.** The library is `−ln P` over residues in
*folded* structures, so it carries `E_local + E_fold,effective`. Upside adds that to `hbond + env +
sidechain`, which are its own model of the fold. The fold's influence on (φ,ψ) is therefore counted
twice. Worse, the term is context-free: it applies a 63% average αL preference to *every* glycine,
which is right for glycines in left-handed turns and actively wrong for glycines inside helices.

**The library's handedness decomposes cleanly by how much chiral context is removed**, and that is
what identifies the double-counted part:

| what is measured | chiral influences present | dG(αR→αL) |
|---|---|---|
| NDRD coil, averaged over neighbours | neighbour side chain **+** fold **+** longer-range sequence | **−1.24** |
| NDRD coil, `GLY\|GLY` entry | fold **+** longer-range sequence | **−0.71** (left), −0.97 (right) |
| AWH capped `Ac-Gly-Gly-NHMe` | none | **0** (measured 0.03 ± noise) |
| AWH capped `Ac-X-Gly-NHMe` | neighbour side chain only | **−0.303** |

Each row removes a source of chirality and the handedness falls. **That monotonic ordering is the
signature of real physics at every level, not of a corrupted dataset**, and it separates the
library's handedness into a part Upside should carry locally and a part it already models
elsewhere. Roughly −0.30 is intrinsic to a glycine next to an L residue, and the remaining ~−0.9
is fold and longer-range context: precisely what `hbond`, `env` and `sidechain` compute for
themselves.

**Correction to an earlier version of this file.** It claimed `GLY|GLY` "must read exactly 0 by
the symmetry argument in §1" and called −0.71 a failed control. That over-applied §1. The symmetry
argument needs the **whole environment** achiral; NDRD's `GLY|GLY` is conditioned only on the
immediate neighbour being glycine, while the rest of the chain is still L-amino acids in a chiral
fold. Only the capped dipeptide is required to read 0, and that is our AWH `LG` control, not a
library entry. Nothing downstream depended on the discarded claim: ff3.1 corrects the coil group
because the sheet group's helical basins are essentially empty (α_R 1.6e-10, α_L 8.7e-16), so its
apparent handedness is a ratio of near-zeros.

---

## 3. Why this applies to glycine and to no other residue

Because glycine is the only residue whose library handedness cannot be local physics. Every other
residue's αL/αR imbalance is explained by Cβ sterics, and the library's own numbers order exactly
that way:

| residue | αL fraction | dG (kT) |
|---|---|---|
| **GLY** | **30.95%** | **−1.238** (the only negative) |
| ASN | 13.4% | +0.35 |
| HIS | 8.0% | ... |
| THR | 0.71% | ... |
| VAL | 0.58% | ... |
| ILE | 0.26% | ... |
| PRO | 0.00% | +12.6 |

All nineteen non-glycine residues are positive (+0.35 to +12.6) and sort by side-chain branching,
which is what local sterics predict. Glycine is the sole outlier, 7× larger than the next, and the
only one where a Cβ-based explanation is unavailable. The ensemble mismatch applies to all 20 by
construction, but only for glycine is it large *and* separable from real physics.

**ASN is the only other residue worth ever checking.**

---

## 4. Why ff3.0's fix over-corrected

ff3.0 replaced the whole central-GLY row with `0.5*(m + mirror(m))` in both the coil and sheet
groups, forcing every glycine map symmetric. That sets `dG(αR→αL)` to **exactly 0** in all contexts.

This has a real precedent: Rosetta's `-score:symmetric_gly_tables` does the same thing, for the same
stated reason — its PDB-derived glycine tables are asymmetric and the flag makes glycine "equally
likely to be in the D- or L-regions". But note two things about that precedent. It is **off by
default**, and its documented purpose is **mixed D/L peptide design**, where D/L equivalence is
genuinely required
([Rosetta D-amino acids](https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids),
[simple_cycpep_predict](https://docs.rosettacommons.org/docs/latest/structure_prediction/simple_cycpep_predict)).
For an all-L protein it is the wrong default, for exactly the reason in §1.

Atomistic force fields sidestep the question entirely by deriving backbone torsions from quantum
chemistry on the glycine dipeptide rather than PDB statistics, which yields a symmetric surface
without any explicit symmetrisation. **The problem is specific to statistical potentials**, and
that is the argument for ff3.1's route: take the glycine surface from a dipeptide calculation, as
the atomistic force fields do, rather than repairing a PDB-derived one.

---

## 5. How we are resolving it: two independent routes

Setting the map from the dipeptide measurement alone was built, and then stopped. A map imported
from small-molecule MD is still an assumption about protein interiors, and it cannot be validated
against the thing it was imported into. So there are now two routes, sharing no input:

| | route | produces |
|---|---|---|
| **Track A** | ConDiv learns the glycine map from the 456-protein training set | a map consistent with Upside's own force field and the native ensemble |
| **Track B** | 2D AWH on all 40 capped `X-Gly` / `Gly-X` dipeptides | an independent physical map from explicit-solvent MD |

**Track A starts from zero handedness**, a symmetrised library row with `A = 0`, so what it learns
is not seeded from Track B and the comparison between them is a real test. Agreement is evidence;
disagreement localises which term is wrong.

### The measurement

2D AWH (GROMACS 2024.4, TIP3P, 300 K) on (φ,ψ) of the central glycine in capped
`Ac-X-Gly-NHMe` dipeptides, ten X, with `Ac-Gly-Gly-NHMe` as the achiral control that must read 0.
Two fully independent replicas with different seeds, each run to its own converged plateau:

| | window | neighbour-averaged dG(αR→αL) | achiral control |
|---|---|---|---|
| replica 1 | 62–100 ns (n=381) | **−0.305 nats** | −0.015 |
| replica 2 | 62–75 ns (n=133) | **−0.319 nats** | +0.094 |

**The answer is −0.303 nats (≈ −0.76 kJ/mol, −0.18 kcal/mol, −0.26 E_up): small, left-handed,
and definitely not zero.** Neither ff2.1 (−1.24) nor ff3.0 (0) is right; the library overstates by
~4×.

Two things the campaign established that are as important as the number:

* **Per-neighbour structure is noise.** The two replicas' neighbour orderings correlate at
  Spearman ρ = +0.048 (p = 0.91), and the disagreement (0.178) did not shrink with sampling. The
  deliverable is one number, not a 40-entry row.
* **The achiral control is necessary but nowhere near sufficient.** Its errors cancel by symmetry —
  exactly the error mode afflicting the chiral systems — so only independent replicas at matched
  sampling bound the uncertainty. This matters because sampling asymmetry is a *known* artifact in
  glycine simulations ([glycine oligomers study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/)),
  and our own first attempt produced a spurious "handedness is zero" from an unconverged flat
  surface.

### Track B's map: the measured surface, held as the reference

**The library's central-glycine coil row is discarded and rebuilt from the AWH surfaces, with no
library data in it at all.** This is `parameters/common/rama31.dat`. It was briefly the training
map and is now the **reference** Track A will be compared against. Built by `py/build_rama_from_awh.py`; the row holds exactly two distinct maps:

| entry | content |
|---|---|
| `X\|GLY`, every neighbour, both directions | the measured surface, `S_meas + A_meas` |
| `GLY\|GLY` | `S_meas` alone, antisymmetric part exactly 0 |

`S_meas` averages the symmetric part of all 19 surfaces (10 dipeptides x 2 replicas, less the one
replica that had not reached LR); `A_meas` averages the antisymmetric part of the 17 chiral ones,
since the achiral blank measures zero rather than a neighbour effect. Both are interpolated from
the AWH 46x46 grid to the library's 72x72 with periodic cubic splines, and parity is re-imposed
exactly on the target grid afterwards. **Nothing is fitted.**

The two-map split is not a special case bolted on: it is §1's symmetry argument applied where it
holds. A glycine flanked by two L-amino acids may be biased; a glycine flanked by glycines may not.
The measurement agrees, the blank's antisymmetric part being 0.029 rms against 0.083 for the chiral
average. The consequence falls out on its own, without a rule: `read_weighted_maps` mixes the left
and right neighbour maps, so `Leu-Gly-Gly` gets **half** the handedness (−0.146) and `Gly-Gly-Gly`
gets **none**.

Result: `dG(αR→αL)` = **−0.303 nats** for `X|Gly` against the library's −1.24, and **exactly 0**
for `Gly|Gly`.

### Units: the map holds −lnP, and that fixes the conversion

**A library map stores `−ln P`, normalised so that `sum(exp(−E)) = 1` over the 72×72 grid.** That
is not an interpretation: every map in `rama.dat` sums to `1.000000`, and `mixture_potential` in
`upside_config.py` states the requirement in its docstring. An AWH PMF is a free energy in kJ/mol
at the run temperature, so the conversion into that slot is

```
E[nats] = PMF[kJ/mol] / kT(300 K),      kT(300 K) = 2.494339 kJ/mol
```

**not** `PMF / 2.914952774272` (kJ/mol per E_up). The two differ by 1.169×. The rule is that a
replacement must be in the same units as the thing it replaces, and the arbitrary additive offset
of the PMF is removed by the normalisation.

This correctly supersedes the earlier choice of the E_up factor, which was made when only an
antisymmetric *energy* correction was being injected and the target was a dG quoted in E_up. Once
the whole surface occupies a `−lnP` slot, the log-probability conversion is the only consistent
one. The handedness is therefore quoted as **−0.303 nats**, the same measurement as the earlier
−0.26 E_up, differing only by which convention it is expressed in.

Upside's own definition is what makes these two so close: `kT = 1 E_up` at `T_up = 1` (350.588 K),
so `kT(300 K) = 0.8557 E_up` and the whole fork is that one factor. Note that ConDiv trains at
`T_up ≈ 0.80`; the library has exactly the same nominal-temperature mismatch, so matching its
convention keeps the two comparable.

### What changed, and what did not

Against `rama.dat`, measured on the built library:

| | ff2.1 | ff3.1 |
|---|---|---|
| `dG(αR→αL)`, mean over 20 left neighbours | −1.238 | **−0.288** |
| αR basin free energy | 2.388 | 2.595 |
| αL basin free energy | 1.206 | 2.293 |
| β basin | 2.044 | 1.346 |
| αR→αL saddle (bottleneck path) | 6.32 | **4.10** |
| map mean over the grid | 11.580 | 11.577 |
| everything outside the coil GLY row | — | bit-identical |

Three things worth reading off that table. **Glycine's αR basin is essentially untouched** (0.207
nats), so the change is not a wholesale re-scaling of glycine: it is the αL basin, 1.087 nats too
deep in the library, coming back to where it was measured. **The map's mean is unchanged to three
decimals**, so glycine's weight relative to the other 19 residue types does not shift. And the
**barrier between the basins falls**, so glycine samples more freely rather than less, despite the
measured surface having a higher global maximum (25.9 against 18.3). That maximum sits in a
genuinely forbidden corner no path crosses, and the library cannot represent it anyway: with
44,112 glycines over 5,184 bins, an empty bin is censored at roughly `ln N ≈ 10.7` above the mean
no matter how forbidden it really is.

### What is given up by using one surface for every neighbour

Honest accounting, because this is the real cost of the change:

* **Per-neighbour structure.** The library resolves it (spread 0.264 nats across its 42 neighbour
  maps in the populated region, on ~1,500 glycines each). ff3.1 asserts none. Our own measurement
  cannot resolve it: per-pair replica agreement is r = 0.41–0.92 at S/N **1.48**, against r = 0.968
  and S/N **3.89** for the average. Using noisy per-pair surfaces would inject roughly 40% noise
  into every map, which is worse than asserting the measured mean.
* **Left/right asymmetry.** The library's left- and right-neighbour glycine maps differ by rms
  0.679 nats. Every dipeptide run has X on the N-terminal side, so the right direction has never
  been measured and receives the same surface.
* Both losses are smaller than the 1.087 nats of αL error being removed, which is why the trade is
  worth making. Neither is permanent: measuring them needs the other 11 left neighbours and all 21
  right neighbours, 32 more AWH systems.

**The library's per-pair ordering could not be kept even if we wanted it**, because it is
*anti-correlated* with the measurement (r = −0.540, rep1 −0.792, rep2 −0.239) while the two
replicas agree with each other at +0.515. It points the wrong way rather than merely being
unresolved, which is also why ff3.0C, built to preserve exactly that structure, scored worse than
plain ff3.0.

* **Coil group only.** The sheet group's helical basins are essentially empty (α_R 1.6e-10,
  α_L 8.7e-16), so its apparent handedness is a ratio of near-zeros, and a capped dipeptide in water measures nothing about a residue in a
  β-sheet. Leaving it alone is not a patch; it is declining to substitute a measurement of a
  different thing.

### Track A: the same two maps, learned instead of measured

Identical parameterisation, `X|Gly = S + A` and `Gly|Gly = S`, so the symmetry constraint is the
same. The difference is where the numbers come from: contrastive divergence against the
456-protein training set, starting at `A = 0`.

**This addresses the §2 defect at its cause rather than correcting for it.** The library
double-counts the fold because it is `−lnP` over folded structures and Upside then adds its own
`hbond + env + sidechain` model of the same fold. Contrastive divergence fits the local term that
reproduces the native ensemble *given the rest of the force field*, which is exactly the
subtraction that double-counting needs.

**The gradient is analytic**, which is the only reason a 5,184-value trainable map is affordable.
`rama_map_pot` is a periodic interpolating bicubic spline built from a tensor product of 1D
solves, so the map enters the energy linearly and separably and `dE/d(map[i,j])` is a
spline-smoothed 2D histogram of the glycine (φ,ψ) samples. Finite differencing would cost 10,369×
a divergence. See `findings.md` §9m; the check that it is the real gradient is
`training/verify_gly_gradient.py`, which failed at 37% on its first run and found a genuine
omission.

**Both tracks assert one map for all 40 neighbour contexts**, for the same reason: neither can
resolve per-pair structure. That is the one assumption they share, and the 40-context AWH campaign
is what tests it.

### Constructions that were built and rejected

Four were built before this one, and the progression is worth keeping because each was killed by a
specific objection rather than by taste:

| | measured mean | `GLY\|GLY` = 0 | fitted params | library data in the GLY row |
|---|---|---|---|---|
| uniform `S + 0.196·A` | yes | no (−0.140) | 1 | yes |
| `S + 0.4245·(A − A_GG)` | yes | yes | 1 | yes |
| shift only `S + A − 1.49·A_GG` | yes | **no (+0.352)** | 1 | yes |
| `S_library + A_measured` | yes | yes | **0** | yes, the symmetric part |
| **measured surface, whole row** | **yes** | **yes** | **0** | **none** |

The first three each need a fitted parameter, and no one-parameter family can satisfy the measured
mean, `GLY|GLY = 0` and the library's per-pair spread at once. That trilemma is itself evidence
that the library's glycine row is not a simple additive artifact of an otherwise sound map.

The fourth removed the fitted parameter but kept the library's symmetric part, and it has a defect
that is obvious once stated: **`S_library` is exactly ff3.0**, the fully mirror-symmetrised map. So
`S_library + A_measured` is ff3.0 plus a correction term, which is a patch on a construction
already established as wrong, no matter that the correction itself is measured. Replacing the row
outright is both cleaner and better supported, since the objection that motivated keeping the
library's symmetric part does not survive measurement (see §5, *What changed, and what did not*).

---

## 6. Additional simulations and checks performed

| what | why | outcome |
|---|---|---|
| 2D AWH, 10 dipeptides × 2 replicas, 100 ns | measure the surface | **−0.303 nats**, replicas agree to 0.012 |
| `Ac-Gly-Gly-NHMe` achiral control | must read 0 | rms 0.032–0.035; **sampling noise, not a defect** (see below) |
| blank vs time, and rep1 vs rep2 | is the residual an artifact? | decays 0.233→0.032 as `1/√t`, uncorrelated across replicas (r = +0.157) |
| AWH extension to 400 ns, both replicas | halve the blank | running (49037819/20) |
| library vs measured, full surface | can the library be discarded? | **yes**: they correlate at r = +0.867 and their αR basins agree to 0.207 nats; the disagreement is the αL basin, 1.087 nats |
| saddle barriers and map mean after replacement | does glycine change weight or freeze? | no: mean 11.580 → 11.577, αR→αL barrier 6.32 → 4.10 nats |
| AWH under ff14SB as well as ff99SB-ILDN | force-field dependence | **agree to 0.045 nats**: -0.246 against -0.291, both blanks on zero |
| unbiased Gly₅ / SAGAS pentapeptides | independent cross-check | **not usable**: Gly₅ reads −0.224 against an exact 0 |
| 32-arm folding benchmark, ff3.0 vs FF2 | what ff3.0 actually did | native **+0.039**, de novo **−0.030**, paired p = 0.021 |
| ConDiv restart from ff2.1 | is the trainer trustworthy? | ff2.1 is a stationary point (sign-flip p ≥ 0.22) |

**The benchmark result gives ff3.1 a falsifiable prediction.** ff3.0 gains on native arms and loses
on de novo ones. If that is because zeroing the glycine αL bias removed turn nucleation — glycine
being the classic turn residue, and de novo folding having to build turns from an extended chain —
then restoring the measured part of it should recover de novo performance while keeping the native
gains. **If the de novo arms do not move, that explanation is wrong.**

### The achiral blank: resolved, and it is the convergence criterion

The `Ac-Gly-Gly-NHMe` control must read exactly 0. It does not: rep1's antisymmetric surface reads
rms 0.032 and rep2's 0.035 at 100 ns, and rep2's basin asymmetry sits at +0.094 where rep1's
decayed to −0.015. If that residual were a defect in the simulation, none of the chiral surfaces
could be used either. It was tested rather than assumed, and it is sampling noise:

* **It does not reproduce between replicas: r = +0.157.** A chirality error in the topology, a
  biased start, or a broken AWH grid would give the *same* residual from any seed. Uncorrelated
  residuals are what finite sampling looks like.
* **It decays while the signal converges.** Blank rms: 0.233 (10 ns) → 0.051 (30) → 0.042 (70) →
  **0.032 (100)**, roughly `1/√t`. The 8-neighbour signal over the same window: 0.080 (40) → 0.075
  (70) → **0.071 (100)**. An artifact would decay too. This contrast is the strongest single piece
  of evidence that the measured handedness is real.

**Blank subtraction does not help** and was rejected. With the error random rather than systematic,
the blank is a noisy estimate of zero; subtracting it adds noise (blank rms 0.0255) instead of
removing bias.

**The honest error bar is wider than first quoted.** A single surface is only S/N 2.2 at 100 ns.
`A_measured` averages 16 surfaces, where the replica-to-replica difference is 0.018 against a
signal of 0.069 (S/N 3.9), but the blank's residual basin asymmetry of +0.039 against the signal's
−0.194 puts roughly **20% uncertainty on the basin dG** — not the ±0.01 that the replica agreement
alone suggests.

**So the blank, not a fixed wall time, is the stopping criterion.** All 10 dipeptides in both
replicas are extended from 100 ns to **400 ns** (jobs 49037819, 49037820), which should halve the
blank and take a single surface to S/N ~4.4. `A_measured` will be re-derived from the extended data
and `rama31.dat` rebuilt if it moves materially. Training on the current map proceeds in parallel,
since the construction is unchanged and only the surface it is built from can shift.

---

## 7. Implementation notes

The chirality operation on the 72×72 grid (φ,ψ from −180° in 5° steps) is `(φ,ψ) → (−φ,−ψ)`, which
is a reversal **with a roll**, because index `i` maps to `(−i) mod 72`:

```python
mirror = lambda m: np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)
```

A plain `m[::-1, ::-1]` is off by one bin and wrong by ~2.3 in practice. This reproduces ff3.0's
`rama3.dat` from `rama.dat` bit-exactly, which is how the convention was confirmed.

Exactly 4.5% of the coil array is NaN, and that is not scattered unpopulated bins: it is **one
whole neighbour column, `CPR`**, which is never read, because `read_rama_maps_and_weights` maps a
cis-proline *neighbour* onto `PRO`. `build_rama_from_awh.py` leaves that column untouched so the
file's structure is unchanged. A naive `abs(a-b) > tol` diff still reports "no difference" across
it, because NaN comparisons are False.

**Verification run on the built library.** Everything outside the coil GLY row is bit-identical,
including the sheet group and both weight arrays; the NaN mask is preserved; all 42 GLY maps
normalise to `1.000000`; the row holds exactly two distinct maps; `GLY|GLY` is antisymmetric-part
zero to machine precision. End to end through `read_weighted_maps`, non-glycine residues are
bit-identical and glycines move from ≈ −1.3 to −0.303, with `Leu-Gly-Gly` at −0.146 and a glycine
between two glycines at 0. In the engine, `1a62` gives a finite total energy of −226.72 against
ff2.1's −228.88.

A **second, separate** glycine symmetrisation also existed: ff3.0's trainer forced the GLY row of
the *rotamer pair-interaction angular profile* palindromic, which is not the rama map and is not
documented anywhere in the ff3.0 write-up. It has been removed; see `findings.md` §9f.1.

---

## Corrections to the previous version

1. **"PDB-derived GLY maps are asymmetric due to finite sampling noise."** False. The asymmetry is
   ~68σ against binomial noise on 44,112 residues, and it matches an independent count of the
   benchmark natives (63% vs 66%). It is a real property of folded proteins; the defect is
   double-counting the fold, not noise.
2. **"Glycine's Ramachandran distribution is intrinsically symmetric about φ=0."** True only for an
   isolated glycine or an achiral neighbourhood. With L-neighbours the conditional distribution has
   no symmetry requirement, and we measure it to be −0.303 nats.
3. **Baruch-Shpigler 2017 was cited as confirming "GLY occupies both φ-quadrants essentially
   equally".** The paper finds the opposite: glycine is *"practically always conformationally
   chiral"*. The citation supported the reverse of its own conclusion.
4. **"The spurious asymmetry applied a net torque on TM4 and drove the helix out; symmetrizing
   fixed it."** Not supported. TM4 melts in **all four** glpG variants including wild type, the
   behaviour **predates ff3.0**, glycine density does **not** predict which helix fails (r = +0.23,
   wrong sign; TM3 is 17.9% glycine at 0.928 occupancy), and ff3.0 only *slows* the decay rather
   than preventing it. See `findings.md`.
5. **"Symmetrize ALL GLY maps unconditionally."** This is ff3.0, now retired: it replaces a −1.24
   error with a −0.30 one rather than removing it.

## Corrections made while building ff3.1

6. **"The library's symmetric part must be kept, because it and the measured one differ by rms
   1.752 E_up with αR basins of 0.193 against 0.089, 25× the handedness correction."** Wrong, and
   it was the argument for a construction that has now been replaced. Those numbers were a
   *per-pair antisymmetric* statistic, not the symmetric parts. Measured properly, the library and
   the dipeptide surfaces correlate at **r = +0.867**, their αR basins agree to **0.207 nats**, and
   the map mean is unchanged to three decimals by the swap. The genuine disagreement is the αL
   basin, 1.087 nats, which is the handedness error itself. The objection that scoped ff3.1 to the
   antisymmetric part did not survive being measured.
7. **The units choice was inverted.** The earlier text gave `PMF / 2.914952774272` (kJ/mol per
   E_up) as correct and `PMF / kT(300 K)` as the 17% error. It is the other way round for this
   purpose: a library map holds `−lnP` normalised to `sum(exp(−E)) = 1`, verified exactly on every
   map in `rama.dat`, so an AWH PMF entering that slot must be divided by `kT`. The measurement did
   not change; −0.26 E_up and −0.303 nats are the same surface in two conventions, and the second
   is the one the file is written in.

## References

- Ting D, Wang G, Shapovalov M, Mitra R, Jordan MI, Dunbrack RL Jr (2010). *Neighbor-Dependent
  Ramachandran Probability Distributions of Amino Acids Developed from a Hierarchical Dirichlet
  Process Model.* PLoS Comput Biol 6(4): e1000763. — the library Upside uses (NDRD_TCB).
  https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000763
- Baruch-Shpigler Y, Wang H, Tuvi-Arad I, Avnir D (2017). *Chiral Ramachandran Plots I: Glycine.*
  Biochemistry 56(42): 5635–5643. — glycine is practically always conformationally chiral.
  https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525
- Shelar A, Chattopadhyay A (2005). *The Ramachandran plots of glycine and pre-proline.*
  BMC Struct Biol 5:14. https://pmc.ncbi.nlm.nih.gov/articles/PMC1201153/
- *Force field dependent solution properties of glycine oligomers* (2015), PMC4450816. — αR/αL
  populations approach equality for poly-Gly only in the limit of sufficient sampling; sampling
  asymmetry is a known artifact. https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/
- Rosetta `-score:symmetric_gly_tables`. — same fix, **off by default**, intended for mixed D/L
  peptides. https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids
- Lovell SC et al. (2003). *Structure validation by Cα geometry: φ,ψ and Cβ deviation.* Proteins
  50:437–450. — glycine's distinct Ramachandran as the validation baseline.
- Hovmöller S, Zhou T, Ohlson T (2002). *Conformations of amino acids in proteins.*
  Acta Crystallogr D 58:768–776.
