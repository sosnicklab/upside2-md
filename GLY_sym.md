# Glycine Ramachandran handedness: the problem, the measurement, and ff3.0

## 0. Where this stands (2026-10-02)

* **Training a context-free map cannot fix glycine.** One map serves helical glycines, which need
  less alpha_L, and natively left-handed loop glycines, which need more. A map trained against
  native structures therefore relearns where evolution placed glycine. Three attempts all show
  this: the full-map row, ff3.0's basin offsets (released 09-30, failed in glpG, validation
  cancelled) and a probe started from equal depth (findings 1.14-1.15).
* **Folded-protein data cannot separate energy from selection**, either through residue counts
  per basin or through per-type context corrections (findings 1.16).
* **ff3.0 takes glycine's row from the measurement and does not train it.**
  `parameters/common/rama31.dat` was built on 2026-10-01 by `build_gly_library.py` (§5).
* **With that map the trainer kept loop glycines left-handed by bending the shared H-bond, sheet
  and side-chain terms, and helical glycines paid.** So plan.md Phase 8 gives glycine its own
  offsets on the three H-bond basin energies, trained by ConDiv with everything else
  (findings 1.17).
* **The experimental check:** our Gly-Gly surface agrees with the GGG spectroscopic distribution
  as closely as ff14SB does on the matched system. No experiment resolves handedness (findings
  1.16).

---

## 1. Why glycine handedness matters, and the exact limit of the symmetry argument

Glycine has no Cβ. Swapping its two Hα substituents maps the molecule onto itself, so an isolated
glycine residue is achiral and its (φ,ψ) free-energy surface must be symmetric under
`(φ,ψ) → (−φ,−ψ)`. Every other residue has a Cβ that breaks this.

But that argument only reaches as far as the molecule's own symmetry. A glycine inside a protein is
flanked by L-amino acids. The local environment is chiral, so the conditional distribution
`P(φ,ψ | glycine, L-neighbours)` has no symmetry requirement at all. The symmetry is exact only
when the neighbours are themselves achiral:

| context | is symmetry required? |
|---|---|
| glycine in an achiral environment (`Gly-Gly-Gly`, poly-Gly) | yes, exactly |
| glycine with any L-amino-acid neighbour (`X-Gly-X`) | no |

This distinction is the substance of the problem, and it is confirmed from two independent
directions. Simulation studies of glycine oligomers find that for Gly₃/Gly₁₀ the αR and αL
populations "would approach equality" in the limit of sufficient sampling
([Force field dependent solution properties of glycine oligomers](https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/)).
And the continuous-chirality analysis of 4,366 glycine residues from 160 high-resolution
structures finds that glycine is "practically always conformationally chiral", with chirality
magnitudes similar to genuinely chiral amino acids
([Baruch-Shpigler et al. 2017](https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525)).

So: achiral molecule, chiral conformation, chiral neighbourhood. Forcing the map symmetric
everywhere throws away something real.

---

## 2. The physical cause: alpha_L glycine is placed, not preferred

Before any measurement, the mechanism.

**Glycine is achiral at C-alpha.** Two hydrogens; swap them and the molecule is unchanged. A
left-handed backbone conformation is intrinsically no more favourable to glycine than a
right-handed one. So glycine's strong occupancy of alpha_L in the PDB (the left-handed
alpha-helical basin at phi ~ +60 deg, psi ~ +45 deg, the exact mirror of alpha_R) cannot be a
property of the residue. It has to come from the environment, through this chain:

1. **L-amino acids build right-handed elements.** Right-handed alpha-helices, right-twisted
   beta-sheets. This is the only fundamental chirality in the system.
2. **Connecting those elements compactly requires occasional left-handed positions.** A chain that
   is uniformly right-handed cannot turn back on itself efficiently. Type I' and II' beta-turns,
   and the left-handed bridge, each need a residue sitting in alpha_L.
3. **alpha_L is sterically forbidden to anything with a C-beta.** Rotating to positive phi drives
   C-beta into the preceding carbonyl. That is the classic Ramachandran exclusion.
4. **Glycine has no C-beta, so alpha_L costs it nothing.** Asn and Asp manage it too, helped by
   side-chain-to-backbone hydrogen bonds, but glycine dominates.
5. **So evolution places glycine wherever a fold needs alpha_L.**

Glycine's alpha_L population is therefore a placement effect: a statement about where glycines are
put, not about what glycine prefers.

### Why that breaks when it is used as a local energy

The statistic is context-specific; the term is context-free. NDRD averages over every environment
and hands Upside one number saying "glycine prefers alpha_L", which is then applied to every
glycine everywhere. Two failure modes result, pointing opposite ways:

* **At turn glycines the force field pushes twice.** The hydrogen bonding and packing that create
  the turn are computed explicitly by `hbond` and `env`; the map adds the statistical shadow of
  those same interactions on top.
* **At helix glycines it pushes the wrong way.** A glycine inside a helix belongs in alpha_R like
  everything around it. None of the turn statistics apply, but the alpha_L bias is applied anyway,
  working against the helical hydrogen bonding that should hold it.

### Why glycine and essentially no other residue

For the other nineteen, C-beta forbids alpha_L regardless of context. Their map value is the same
in a turn as in a helix, it is set by local sterics, and it is genuine local physics that belongs
in a local term. Glycine is the one residue whose alpha_L occupancy is decided by architecture
rather than by its own atoms, so it is the one residue whose map entry is largely a record of where
the fold put it.

### The part that is real

A capped dipeptide with a single L-neighbour still shows a modest left-handed bias, and that one is
genuine local physics: the neighbour's C-beta and carbonyl make the immediate environment weakly
chiral with no fold involved. That belongs in the map. The disagreement is over how much of the
library's value is this, and how much is architecture.

---

## 2a. What the numbers say

The library's handedness decomposes by how much chiral context is removed, and that is what
identifies the double-counted part (dG(αR→αL) in nats):

| what is measured | chiral influences present | dG(αR→αL) |
|---|---|---|
| NDRD coil, averaged over neighbours | neighbour side chain + fold + longer-range sequence | −1.24 |
| NDRD coil, `GLY\|GLY` entry | fold + longer-range sequence | −0.71 (left), −0.97 (right) |
| AWH capped `Ac-X-Gly-NHMe` / `Ac-Gly-X-NHMe`, 400 ns | neighbour side chain only | −0.20 (left), −0.11 (right) |
| AWH capped `Ac-Gly-Gly-NHMe` | none | 0 within noise |

Each row removes a source of chirality and the handedness falls. That monotonic ordering is the
signature of real physics at every level, not of a corrupted dataset, and it separates the
library's handedness into a part Upside should carry locally and a part it already models
elsewhere: about −0.2 is intrinsic to a glycine next to an L residue, and the remaining ~−1.0 is
fold and longer-range context, which `hbond`, `env` and `sidechain` compute for themselves.

Only the capped dipeptide is required to read 0. The symmetry argument in §1 needs the whole
environment achiral; NDRD's `GLY|GLY` is conditioned only on the immediate neighbour being glycine,
while the rest of the chain is still L-amino acids in a chiral fold, so its −0.71 is fold context
and not a failed control. The table uses the coil group because the sheet group's helical basins
are essentially empty (α_R 1.6e-10, α_L 8.7e-16), so its apparent handedness is a ratio of
near-zeros.

---

## 3. Why this applies to glycine and to no other residue

Every other residue's αL/αR imbalance is explained by Cβ sterics, and the library's own numbers
order exactly that way. Glycine holds 30.95% of its coil-map weight in αL, dG(αR→αL) = −1.238 kT,
the only negative value; the other nineteen are all positive, from +0.35 (Asn, 13.4% αL) to
+12.6 (Pro, 0%), sorted by side-chain branching, which is what local sterics predict (table in
findings 1.18). The ensemble mismatch applies to all 20 by construction, but only for glycine is
it large and separable from real physics. Asn is the only other residue worth ever checking.

---

## 4. Why the first ff3.0's symmetrisation over-corrected

The first ff3.0 (retired) replaced the whole central-GLY row with `0.5*(m + mirror(m))` in both the
coil and sheet groups, which sets `dG(αR→αL)` to exactly 0 in every context. Rosetta's
`-score:symmetric_gly_tables` does the same, but it is off by default and documented for mixed D/L
peptide design, where D/L equivalence is genuinely required
([Rosetta D-amino acids](https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids),
[simple_cycpep_predict](https://docs.rosettacommons.org/docs/latest/structure_prediction/simple_cycpep_predict)).
For an all-L protein it is the wrong default, for the reason in §1.

Atomistic force fields avoid the question by deriving backbone torsions from quantum chemistry on
the glycine dipeptide rather than from PDB statistics. The problem is specific to statistical
potentials, and that is the argument for the measured-map route: take the glycine surface from a
dipeptide calculation, as the atomistic force fields do, rather than repairing a PDB-derived one.

---

## 5. How it is resolved

### The measurement

2D AWH (GROMACS 2024.4, amber99sb-ildn, TIP3P, 300 K) on (φ,ψ) of the central glycine in capped
dipeptides, all 40 contexts: `Ac-X-Gly-NHMe` for a left neighbour and `Ac-Gly-X-NHMe` for a right
one, with `Ac-Gly-Gly-NHMe` as the achiral control that must read 0. Replica 1 ran every context to
400 ns; an independent replica 2 with different seeds is the reproducibility check.

* **Result:** left neighbours −0.20 nats, right −0.11, context-averaged −0.154, the like-for-like
  number for Upside's mixture of a residue's left and right maps (findings 9s). Neither ff2.1
  (−1.24) nor the symmetrised ff3.0 (0) is right. An early 100 ns estimate from ten left
  neighbours, −0.30, is superseded.
* **It does not depend on the force field:** ff14SB agrees with ff99SB-ILDN to 0.045 nats, both
  blanks on zero (findings 9r). Both are AMBER; CHARMM36m has not been run.
* **Per-neighbour structure is noise.** The two replicas' neighbour orderings correlate at
  Spearman ρ = +0.048 (p = 0.91), and the disagreement did not shrink with sampling. Per-pair
  replica agreement is r = 0.41–0.92 at S/N 1.48, against r = 0.968 and S/N 3.89 for the
  neighbour average. The one context that is clearly different is a proline on the right: its ring
  empties both helical basins.
* **The achiral control is necessary but nowhere near sufficient.** Its errors cancel by symmetry,
  which is exactly the error mode of the chiral systems, so only independent replicas at matched
  sampling bound the uncertainty. Sampling asymmetry is a known artifact in glycine simulations
  ([glycine oligomers study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/)), and our own first
  attempt produced a spurious "handedness is zero" from an unconverged flat surface (findings 12a).

### The map: the measured surface is ff3.0's glycine row

The library's central-glycine row, coil and sheet, is discarded and rebuilt from the AWH surfaces,
with no library data in it. `parameters/common/rama31.dat` is built by
`build_gly_library.py` (a one-time builder, kept in scratchpad because its AWH input exists only on
the cluster) from replica 1 at 400 ns. Each neighbour class gets its own
measurement:

| entry | content |
|---|---|
| `GLY\|GLY` | the two Gly-Gly surfaces (`LG`, `RG`) pooled and made exactly mirror-symmetric |
| `GLY\|right\|PRO` | the Ac-Gly-Pro-NHMe surface (`RP`) alone: the ring empties both helical basins (alpha_R 0.006, alpha_L 0.004) |
| every other `GLY\|X`, both directions | the remaining 37 surfaces pooled |

* Pooling and symmetrising average probabilities, not energies. Each surface is interpolated from
  the AWH 46x46 grid to the library's 72x72 with periodic cubic splines, in energy.
* Each entry is stored with Upside's reference correction (`rama_map_pot_ref`) subtracted, so the
  engine applies the measured surface itself.
* The same entry is written to the sheet group, because the sheet maps are strand statistics,
  which is selection again. Nothing is fitted.

**Result, as applied by the engine (checked through `upside_config`):**
* X-G-Y glycines: ln(aR/aL) = −0.120 at T = 1, against −0.95 to −1.37 for NDRD.
* A glycine between glycines: exactly 0.
* `Leu-Gly-Gly`: −0.064, from the unchanged left/right mixture.
* Every non-glycine residue is identical to ff2.1's.

**The 2026-09-24 build had two defects, both fixed in the 10-01 build:**
1. It gave `GLY|right|PRO` the pooled surface (alpha_R 0.105), which erased the proline clash.
2. It omitted the reference correction, which reshapes a glycine map: helical basins 0.105 ->
   0.145.

Against the library, the measured surface of the earlier 100 ns build correlated at r = +0.867 over
the whole grid, its αR basin agreed to 0.207 nats and the map mean was unchanged to three
decimals; the disagreement sat in the αL basin, 1.087 nats, which is the handedness error itself.
The library cannot represent the forbidden corners anyway: its TCB set holds 44,112 residues of all
types (Ting et al. 2010), glycine only a share of them, so over 5,184 bins an empty glycine bin is
censored below `ln 44112 ≈ 10.7` above the mean however forbidden it really is.

### Units: the map holds −lnP, and that fixes the conversion

A library map stores `−ln P`, normalised so that `sum(exp(−E)) = 1` over the 72×72 grid. Every map
in `rama.dat` sums to `1.000000`, and `mixture_potential` in `upside_config.py` states the
requirement in its docstring. An AWH PMF is a free energy in kJ/mol at the run temperature, so the
conversion into that slot is

```
E[nats] = PMF[kJ/mol] / kT(300 K),      kT(300 K) = 2.494339 kJ/mol
```

and not `PMF / 2.914952774272` (kJ/mol per E_up); the two differ by 1.169×. A replacement must be
in the units of the thing it replaces, and the normalisation removes the PMF's arbitrary additive
offset.

Upside's own definition is what makes the two so close: `kT = 1 E_up` at `T_up = 1` (350.588 K), so
`kT(300 K) = 0.8557 E_up` and the whole difference is that one factor. ConDiv trains at
`T_up ≈ 0.80`; the library has the same nominal-temperature mismatch, so matching its convention
keeps the two comparable.

### What is given up by one surface per neighbour class

* **Per-neighbour structure.** The library resolves it (spread 0.264 nats across its 42 neighbour
  maps in the populated region, on ~1,500 glycines each). `rama31.dat` asserts none, because the
  measurement cannot resolve it (S/N 1.48 per pair): per-pair surfaces would inject roughly 40%
  noise into every map, which is worse than asserting the measured mean.
* **Left/right asymmetry.** The pool merges left and right contexts, which the measurement puts at
  −0.20 and −0.11; the library's left and right glycine maps differ by rms 0.679 nats.
* Both losses are smaller than the αL error being removed.
* The library's per-pair ordering could not be kept even if wanted: on the 100 ns data it is
  anti-correlated with the measurement (r = −0.540; rep1 −0.792, rep2 −0.239) while the two
  replicas agree with each other at +0.515. It points the wrong way rather than merely being
  unresolved.

A glycine map learned by ConDiv (Track A) was tried with the FF1-form trainer and settled at
−0.885, re-deriving a PDB-like handedness; it is retired with that trainer (findings 9e, 9s).

---

## 6. Simulations and checks performed

| what | why | outcome |
|---|---|---|
| 2D AWH, 10 dipeptides × 2 replicas, 100 ns | does the neighbour average reproduce? | yes: −0.305 against −0.319 nats; the per-neighbour ordering does not (Spearman +0.048) |
| `Ac-Gly-Gly-NHMe` achiral control | must read 0 | rms 0.032–0.035 at 100 ns; sampling noise, not a defect (see below) |
| blank vs time, and rep1 vs rep2 | is the residual an artifact? | decays 0.233→0.032 as `1/√t`, uncorrelated across replicas (r = +0.157) |
| AWH to 400 ns, all 40 contexts | halve the blank, add right neighbours | done (replica 1): left −0.20, right −0.11, context-averaged −0.154 |
| AWH under ff14SB as well as ff99SB-ILDN | force-field dependence | agree to 0.045 nats: −0.246 against −0.291, both blanks on zero |
| unbiased Gly₅ / SAGAS pentapeptides | independent cross-check | not usable: Gly₅ reads −0.224 against an exact 0 |

The retired symmetrised ff3.0 gained on native benchmark arms and lost on de novo ones (findings
9d); if that came from removing glycine's turn-nucleating αL bias, the measured map should not
show it, which the Peng validation in plan.md Phase 8 tests.

### The achiral blank: resolved, and it is the convergence criterion

The `Ac-Gly-Gly-NHMe` control must read exactly 0. At 100 ns it did not: rep1's antisymmetric
surface read rms 0.032 and rep2's 0.035, and rep2's basin asymmetry sat at +0.094 where rep1's had
decayed to −0.015. If that residual were a defect in the simulation, none of the chiral surfaces
could be used either. It was tested rather than assumed, and it is sampling noise:

* **It does not reproduce between replicas: r = +0.157.** A chirality error in the topology, a
  biased start, or a broken AWH grid would give the same residual from any seed. Uncorrelated
  residuals are what finite sampling looks like.
* **It decays while the signal converges.** Blank rms: 0.233 (10 ns) → 0.051 (30) → 0.042 (70) →
  0.032 (100), roughly `1/√t`. The 8-neighbour signal over the same window: 0.080 (40) → 0.075
  (70) → 0.071 (100). An artifact would decay too. This contrast is the strongest single piece of
  evidence that the measured handedness is real.

Blank subtraction does not help and was rejected. With the error random rather than systematic,
the blank is a noisy estimate of zero; subtracting it adds noise (blank rms 0.0255) instead of
removing bias.

The error bar at 100 ns was wider than the replica agreement suggests. A single surface was only
S/N 2.2; the 16-surface average had a replica-to-replica difference of 0.018 against a signal of
0.069 (S/N 3.9), but the blank's residual basin asymmetry of +0.039 against the signal's −0.194 put
roughly 20% uncertainty on the basin dG, not ±0.01.

So the blank, not a fixed wall time, is the stopping criterion. The dipeptides were extended from
100 ns to 400 ns on it, which should halve the blank and take a single surface to S/N ~4.4, and
`rama31.dat` is built from the extended data.

---

## 7. Implementation notes

The chirality operation on the 72×72 grid (φ,ψ from −180° in 5° steps) is `(φ,ψ) → (−φ,−ψ)`, which
is a reversal with a roll, because index `i` maps to `(−i) mod 72`:

```python
mirror = lambda m: np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)
```

A plain `m[::-1, ::-1]` is off by one bin and wrong by ~2.3 in practice. This reproduced the
first ff3.0's `rama3.dat` (since removed; in git history) from `rama.dat` bit-exactly, which is how
the convention was confirmed.

Exactly 4.5% of the coil array is NaN, and that is not scattered unpopulated bins: it is one whole
neighbour column, `CPR`, which is never read, because `read_rama_maps_and_weights` maps a
cis-proline neighbour onto `PRO`. `build_gly_library.py` leaves that column untouched so the file's
structure is unchanged. A naive `abs(a-b) > tol` diff still reports "no difference" across it,
because NaN comparisons are False.

**Verification, run by the builder itself on every build (2026-10-01):** through
`read_weighted_maps` with ff2.1's sheet energies, every non-glycine residue's map is identical to
the source library's; a glycine whose maps are all one entry equals that measured surface, plus
the reference correction, to 1e-6 (X-G-Y and terminal glycines the pool, G-G-G and terminal G-G the
symmetric map); and those reading `GLY|GLY` are mirror-symmetric to 2e-6. An end-to-end
`upside_config` build of `AGATGVGGGLKGPS` with and without the new library leaves every
non-glycine residue's `rama_pot` identical. The library built on midway2 from that tree's files is
byte-identical to the local build.

A second, separate glycine symmetrisation also existed: the first ff3.0's trainer forced the GLY
row of the rotamer pair-interaction angular profile palindromic, which is not the rama map and was
not documented in the ff3.0 write-up. It has been removed; see `findings.md` 9e.

---

## References

- Ting D, Wang G, Shapovalov M, Mitra R, Jordan MI, Dunbrack RL Jr (2010). *Neighbor-Dependent
  Ramachandran Probability Distributions of Amino Acids Developed from a Hierarchical Dirichlet
  Process Model.* PLoS Comput Biol 6(4): e1000763. The library Upside uses (NDRD_TCB).
  https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000763
- Baruch-Shpigler Y, Wang H, Tuvi-Arad I, Avnir D (2017). *Chiral Ramachandran Plots I: Glycine.*
  Biochemistry 56(42): 5635–5643. Glycine is practically always conformationally chiral.
  https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525
- Shelar A, Chattopadhyay A (2005). *The Ramachandran plots of glycine and pre-proline.*
  BMC Struct Biol 5:14. https://pmc.ncbi.nlm.nih.gov/articles/PMC1201153/
- *Force field dependent solution properties of glycine oligomers* (2015), PMC4450816. αR/αL
  populations approach equality for poly-Gly only in the limit of sufficient sampling; sampling
  asymmetry is a known artifact. https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/
- Rosetta `-score:symmetric_gly_tables`. The same fix, off by default, intended for mixed D/L
  peptides. https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids
- Lovell SC et al. (2003). *Structure validation by Cα geometry: φ,ψ and Cβ deviation.* Proteins
  50:437–450. Glycine's distinct Ramachandran as the validation baseline.
- Hovmöller S, Zhou T, Ohlson T (2002). *Conformations of amino acids in proteins.*
  Acta Crystallogr D 58:768–776.
