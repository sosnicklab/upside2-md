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
[Ting, Wang, Shapovalov, Mitra, Jordan & Dunbrack, *PLoS Comput Biol* **6**(4): e1000763 (2010)](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000763),
identified in our tree by direct comparison (correlation 1.00000, max deviation 0.0000). `TCB` =
Turn, Coil and Bridge, i.e. it is fitted to **loop** residues.

Measured on our own basins (`αR: φ∈[−100,−40], ψ∈[−60,10]`; `αL: φ∈[40,100], ψ∈[−10,60]`; exact
mirror images, 195 grid cells each), the library's central-glycine row gives a neighbour-averaged
`dG(αR→αL)` of **−1.32 E_up** — a large left-handed preference.

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

**The library fails its own achiral control, which proves the point without any simulation.**
`GLY|GLY` must read exactly 0 by the symmetry argument in §1. Measured:

| entry | dG(αR→αL) | required |
|---|---|---|
| coil `GLY\|GLY`, left neighbour | **−0.7095** | 0 |
| coil `GLY\|GLY`, right neighbour | **−0.9717** | 0 |
| sheet `GLY\|GLY`, left neighbour | +0.0011 | 0 |

**More than half of the library's glycine handedness is present in the one context where it must be
zero.** That is an internal estimate of the systematic error, and it agrees with the independent
simulation estimate of ~1.06 (below). Note the sheet group *passes*, which is why ff3.1 corrects
the coil group only.

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
without any explicit symmetrisation. **The problem is specific to statistical potentials.**

---

## 5. How we resolved it: measure the number, then set it

### The measurement

2D AWH (GROMACS 2024.4, TIP3P, 300 K) on (φ,ψ) of the central glycine in capped
`Ac-X-Gly-NHMe` dipeptides, ten X, with `Ac-Gly-Gly-NHMe` as the achiral control that must read 0.
Two fully independent replicas with different seeds, each run to its own converged plateau:

| | window | neighbour-averaged dG(αR→αL) | achiral control |
|---|---|---|---|
| replica 1 | 62–100 ns (n=381) | **−0.261 E_up** | −0.015 |
| replica 2 | 62–75 ns (n=133) | **−0.273 E_up** | +0.094 |

**The answer is −0.26 E_up (≈ −0.76 kJ/mol, −0.18 kcal/mol): small, left-handed, and definitely
not zero.** Neither ff2.1 (−1.32) nor ff3.0 (0) is right; the library overstates by ~4×.

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

### The correction: ff3.1

`M_new = S + λ·A`, where `S` and `A` are the mirror-symmetric and antisymmetric parts of the GLY
row. `λ = 1` is ff2.1, `λ = 0` is ff3.0, and `dG` is linear in λ (`−1.318·λ`), so calibration is a
division rather than a fit.

* **λ = 0.20.** The two replicas imply 0.197 and 0.206. The systematic uncertainty makes λ good to
  only ±0.07, so more digits would be false precision.
* **`GLY|GLY` forced to λ = 0**, fully symmetric — physics, not fitting, and the one entry whose
  true value is known exactly without measuring.
* **Coil group only.** The sheet group passes its own achiral control (+0.0011) and its helical
  basins are essentially empty, so scaling it would be scaling noise.
* The library's per-neighbour structure is **kept but scaled** (spread shrinks 5×, −1.54..−0.87 →
  −0.31..−0.17). Our measurement could not resolve neighbour dependence, which is not a reason to
  erase what the PDB statistics do resolve.

Built as `parameters/common/rama31.dat`. Verified at three levels: only the coil GLY row differs
from `rama.dat`; `GLY|GLY` reads −0.000000 and is exactly self-mirror-symmetric; and end-to-end in
a built config, GLY residues move from −1.234 to −0.248 while non-GLY residues are bit-identical.

---

## 6. Additional simulations and checks performed

| what | why | outcome |
|---|---|---|
| 2D AWH, 10 dipeptides × 2 replicas, 100 ns | measure the true value | **−0.26 E_up**, replicas agree to 0.012 |
| `Ac-Gly-Gly-NHMe` achiral control | must read 0 | rep1 −0.015 ✓, rep2 +0.094 (**open anomaly**) |
| AWH under ff14SB as well as ff99SB-ILDN | force-field dependence | in progress |
| unbiased Gly₅ / SAGAS pentapeptides | independent cross-check | **not usable**: Gly₅ reads −0.224 against an exact 0 |
| 32-arm folding benchmark, ff3.0 vs FF2 | what ff3.0 actually did | native **+0.039**, de novo **−0.030**, paired p = 0.021 |
| ConDiv restart from ff2.1 | is the trainer trustworthy? | ff2.1 is a stationary point (sign-flip p ≥ 0.22) |

**The benchmark result gives ff3.1 a falsifiable prediction.** ff3.0 gains on native arms and loses
on de novo ones. If that is because zeroing the glycine αL bias removed turn nucleation — glycine
being the classic turn residue, and de novo folding having to build turns from an extended chain —
then restoring 20% of it should recover de novo performance while keeping the native gains. **If
the de novo arms do not move, that explanation is wrong.**

### Open

`rep2`'s achiral control sits at +0.094 where `rep1`'s decayed to −0.015. Unexplained. It does not
propagate into the answer (the controls differ by 0.109 while the chiral averages agree to 0.012),
but it is the honest floor on quoting any uncertainty below ~0.1 E_up.

---

## 7. Implementation notes

The chirality operation on the 72×72 grid (φ,ψ from −180° in 5° steps) is `(φ,ψ) → (−φ,−ψ)`, which
is a reversal **with a roll**, because index `i` maps to `(−i) mod 72`:

```python
mirror = lambda m: np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)
```

A plain `m[::-1, ::-1]` is off by one bin and wrong by ~2.3 in practice. This reproduces ff3.0's
`rama3.dat` from `rama.dat` bit-exactly, which is how the convention was confirmed. About 4.5% of
the coil array is NaN (unpopulated bins) and the NaN pattern is itself mirror-symmetric, so no
special handling is needed — but a naive `abs(a-b) > tol` diff reports "no difference" because NaN
comparisons are False.

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
   no symmetry requirement, and we measure it to be −0.26 E_up.
3. **Baruch-Shpigler 2017 was cited as confirming "GLY occupies both φ-quadrants essentially
   equally".** The paper finds the opposite: glycine is *"practically always conformationally
   chiral"*. The citation supported the reverse of its own conclusion.
4. **"The spurious asymmetry applied a net torque on TM4 and drove the helix out; symmetrizing
   fixed it."** Not supported. TM4 melts in **all four** glpG variants including wild type, the
   behaviour **predates ff3.0**, glycine density does **not** predict which helix fails (r = +0.23,
   wrong sign; TM3 is 17.9% glycine at 0.928 occupancy), and ff3.0 only *slows* the decay rather
   than preventing it. See `findings.md`.
5. **"Symmetrize ALL GLY maps unconditionally."** This is ff3.0, now retired: it replaces a −1.32
   error with a −0.26 one rather than removing it.

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
