# Glycine Ramachandran handedness: the problem, the measurement, and ff3.0

## 0. Where this stands

This record is largely historical. It keeps the symmetry argument (§1-§3), why the first ff3.0's symmetrisation was
wrong (§4), and the AWH dipeptide measurement and the library built from it on 2026-10-01 (§5-§7). The current
direction is elsewhere:
* **ff3.0's glycine map** is NDRD plus one trained depth pair on the 37 central-glycine coil maps (alpha_R -0.1068,
  alpha_L +0.5169; findings 1.27, up.md 2.8). Round 4 (plan.md Phase 11) ran both of its options, a glycine map
  fitted to selection-free data and frozen (the `ff30_bio` training, dropped 10-11 with the frozen arms, user) and
  only glycine's depth trained on ff2.1's map, which is ff3.0's. The depth training for ff3.1 is reviewed in findings
  1.34.
* **`parameters/common/rama31.dat`** holds, since 2026-10-05, the glycine row fitted to BioEmu's plain-MD
  octapeptides (findings 1.19, up.md 2.8), the fixed library of `ff30_bio`, not ff3.0's. The AWH build of §5 is in
  git history and `backup/rama31.dat.bak_pre_bioemu_20261005`.

Lessons that still hold (2026-10-07):
* **Training a context-free map cannot fix glycine.** One map serves helical glycines, which need less alpha_L, and
  natively left-handed loop glycines, which need more, so a map trained against native structures relearns where
  evolution placed glycine. Three attempts show this: the full-map row, ff3.0's basin offsets of 09-30 (withdrawn;
  its glpG run had mixed coverage tables, findings 3.11) and a probe started from equal depth (findings 1.14-1.15).
* **Folded-protein data cannot separate energy from selection,** through residue counts per basin or per-type
  context corrections (findings 1.16).
* With the AWH map fixed, the trainer kept loop glycines left-handed by bending the shared H-bond, sheet and
  side-chain terms, and helical glycines paid; glycine's own offsets on the three H-bond basin energies (round 3)
  went left-handed (findings 1.17).
* **The experimental check:** our Gly-Gly surface agrees with the GGG spectroscopic distribution as closely as
  ff14SB does on the matched system. No experiment resolves handedness (findings 1.16).

---

## 1. Why glycine handedness matters, and the exact limit of the symmetry argument

Glycine has no Cβ, so an isolated glycine is achiral and its (φ,ψ) free-energy surface must be symmetric under
`(φ,ψ) → (−φ,−ψ)`. A glycine inside a protein is flanked by L-amino acids, so the conditional distribution
`P(φ,ψ | glycine, L-neighbours)` has no symmetry requirement. The symmetry is exact only in an achiral environment
(`Gly-Gly-Gly`, poly-Gly); with any L neighbour (`X-Gly-X`) it is not required.

Two independent directions confirm it. Simulations of Gly₃/Gly₁₀ find the αR and αL populations "would approach
equality" in the limit of sufficient sampling
([glycine oligomers](https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/)), and the continuous-chirality analysis of
4,366 glycine residues from 160 high-resolution structures finds glycine "practically always conformationally chiral"
([Baruch-Shpigler et al. 2017](https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525)). Forcing the map symmetric
everywhere throws away something real.

---

## 2. The physical cause: alpha_L glycine is placed, not preferred

Glycine's strong alpha_L occupancy in the PDB (phi ~ +60 deg, psi ~ +45 deg, the mirror of alpha_R) cannot be a
property of the achiral residue. L-amino acids build right-handed helices and right-twisted sheets; connecting them
compactly needs occasional left-handed positions (type I' and II' turns, the left-handed bridge); alpha_L is
sterically forbidden to anything with a C-beta, which positive phi drives into the preceding carbonyl; glycine has
none (Asn and Asp manage it with side-chain-to-backbone H-bonds), so evolution places glycine wherever a fold needs
alpha_L. Its alpha_L population is a placement effect, a statement about where glycines are put.

Used as a context-free local energy, the statistic fails both ways. At turn glycines the map adds the statistical
shadow of the H-bonding and packing that `hbond` and `env` already compute; at helix glycines it pushes against the
helical H-bonds. For the other nineteen residues C-beta forbids alpha_L in every context, so their map value is local
physics. The real part of glycine's bias: a capped dipeptide with a single L neighbour still shows a modest
left-handed bias from the neighbour's C-beta and carbonyl, and that belongs in the map.

---

## 2a. What the numbers say

The library's handedness by how much chiral context is removed (dG(αR→αL) in nats):

| what is measured | chiral influences present | dG(αR→αL) |
|---|---|---|
| NDRD coil, averaged over neighbours | neighbour side chain + fold + longer-range sequence | −1.24 |
| NDRD coil, `GLY\|GLY` entry | fold + longer-range sequence | −0.71 (left), −0.97 (right) |
| AWH capped `Ac-X-Gly-NHMe` / `Ac-Gly-X-NHMe`, 400 ns | neighbour side chain only | −0.20 (left), −0.11 (right) |
| AWH capped `Ac-Gly-Gly-NHMe` | none | 0 within noise |

Each row removes a source of chirality and the handedness falls, the signature of real physics at every level: about
−0.2 is intrinsic to a glycine next to an L residue, and the remaining ~−1.0 is fold and longer-range context, which
`hbond`, `env` and `sidechain` compute for themselves. Only the capped dipeptide must read 0: NDRD's `GLY|GLY` is
conditioned only on the immediate neighbour, with the rest of the chain L-amino acids in a chiral fold, so its −0.71
is fold context, not a failed control. The coil group is used because the sheet group's helical basins are
essentially empty (α_R 1.6e-10, α_L 8.7e-16).

---

## 3. Why this applies to glycine and to no other residue

Glycine holds 30.95% of its coil-map weight in αL, dG(αR→αL) = −1.238 kT, the only negative value; the other nineteen
are positive, from +0.35 (Asn, 13.4% αL) to +12.6 (Pro, 0%), sorted by side-chain branching as local sterics predict
(table in findings 1.18). Only for glycine is the ensemble mismatch large and separable from real physics; Asn is the
only other residue worth ever checking.

---

## 4. Why the first ff3.0's symmetrisation over-corrected

The first ff3.0 (retired) replaced the whole central-GLY row with `0.5*(m + mirror(m))` in coil and sheet, setting
`dG(αR→αL)` to exactly 0 in every context. Rosetta's `-score:symmetric_gly_tables` does the same, but is off by
default and meant for mixed D/L peptide design
([Rosetta D-amino acids](https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids));
for an all-L protein it is wrong (§1). Atomistic force fields avoid the question by deriving backbone torsions from
quantum chemistry on the glycine dipeptide, so the problem is specific to statistical potentials. The retired
symmetrised ff3.0 gained on native benchmark arms and lost on de novo ones (findings 9d).

---

## 5. How it is resolved (the AWH build of 2026-10-01)

**The measurement.** 2D AWH (GROMACS 2024.4, amber99sb-ildn, TIP3P, 300 K) on (φ,ψ) of the central glycine in capped
dipeptides, all 40 contexts: `Ac-X-Gly-NHMe` for a left neighbour, `Ac-Gly-X-NHMe` for a right one, and
`Ac-Gly-Gly-NHMe` as the achiral control that must read 0. Replica 1 ran every context to 400 ns; replica 2, with
other seeds, is the reproducibility check.
* Left neighbours −0.20 nats, right −0.11, context-averaged −0.154, the like-for-like number for Upside's mixture of
  a residue's left and right maps (findings 9s). Neither ff2.1 (−1.24) nor the symmetrised ff3.0 (0) is right.
* It does not depend on the force field within AMBER: ff14SB agrees with ff99SB-ILDN to 0.045 nats (−0.246 against
  −0.291), both blanks on zero (findings 9r). CHARMM36m has not been run.
* Per-neighbour structure is noise. The two replicas' neighbour orderings correlate at Spearman ρ = +0.048
  (p = 0.91), and the disagreement did not shrink with sampling; per-pair replica agreement is r = 0.41–0.92 at
  S/N 1.48, against r = 0.968 and S/N 3.89 for the neighbour average. The one clearly different context is a proline
  on the right, whose ring empties both helical basins.

**The map.** `build_gly_library.py` (a one-time builder, kept in scratchpad because its AWH input exists only on the
cluster) rebuilt the central-glycine row, coil and sheet, from replica 1 at 400 ns with no library data in it:
`GLY|GLY` the two Gly-Gly surfaces pooled and made exactly mirror-symmetric, `GLY|right|PRO` the Ac-Gly-Pro-NHMe
surface alone (alpha_R 0.006, alpha_L 0.004), every other `GLY|X` the remaining 37 pooled. Pooling averages
probabilities, not energies. Each entry is stored with `rama_map_pot_ref` subtracted, so the engine applies the
measured surface (omitting it, as the 09-24 build did, reshapes a glycine map: helical basins 0.105 -> 0.145), and
the sheet group gets the same entry, because sheet maps are strand statistics, selection again. As applied by the
engine: X-G-Y glycines ln(aR/aL) = −0.120 at T = 1, against −0.95 to −1.37 for NDRD; a glycine between glycines
exactly 0; every non-glycine residue identical to ff2.1's.

**Units: the map holds −lnP.** A library map stores `−ln P`, normalised so that `sum(exp(−E)) = 1` over the 72×72
grid (every map in `rama.dat` sums to `1.000000`; `mixture_potential` in `upside_config.py` states it). An AWH PMF in
kJ/mol at the run temperature converts as `E[nats] = PMF[kJ/mol] / kT(300 K)`, `kT(300 K) = 2.494339 kJ/mol`, not
`PMF / 2.914952774272` (kJ/mol per E_up); the two differ by 1.169×, and the normalisation removes the PMF's
additive offset. Since `kT = 1 E_up` at `T_up = 1` (350.588 K), `kT(300 K) = 0.8557 E_up`. ConDiv trains at
`T_up ≈ 0.80`; the library has the same nominal-temperature mismatch, so matching its convention keeps the two
comparable.

**Given up by one surface per neighbour class:** the library's per-neighbour structure (spread 0.264 nats across its
42 neighbour maps) and its left/right difference (rms 0.679 nats). The measurement cannot resolve per-pair surfaces
(S/N 1.48, roughly 40% noise per map), and on the 100 ns data the library's per-pair ordering is anti-correlated with
it (r = −0.540) while the replicas agree at +0.515; both losses are smaller than the αL error removed. The library
cannot represent the forbidden corners either: its TCB set holds 44,112 residues of all types (Ting et al. 2010), so
an empty glycine bin is censored below `ln 44112 ≈ 10.7` above the mean.

A glycine map learned by ConDiv with the FF1-form trainer (Track A) settled at −0.885, re-deriving a PDB-like
handedness; it is retired with that trainer (findings 9e, 9s).

---

## 6. Simulations and checks performed

Besides §5: 10 dipeptides × 2 replicas at 100 ns reproduced the neighbour average (−0.305 against −0.319 nats) but
not the per-neighbour ordering (Spearman +0.048); unbiased Gly₅ / SAGAS pentapeptides were not usable as a cross-check
(Gly₅ reads −0.224 against an exact 0).

**The achiral blank is the convergence criterion.** At 100 ns the control read rms 0.032 (rep1) and 0.035 (rep2),
not 0. It is sampling noise, tested rather than assumed: the residuals do not reproduce between replicas (r = +0.157),
where a chirality error in the topology, a biased start or a broken AWH grid would; and they decay roughly as `1/√t`
(0.233 at 10 ns, 0.051 at 30, 0.042 at 70, 0.032 at 100) while the 8-neighbour signal converges (0.080 at 40, 0.075
at 70, 0.071 at 100). Blank subtraction was rejected: subtracting a noisy estimate of zero adds noise (blank rms
0.0255). The blank's residual basin asymmetry (+0.039 against the signal's −0.194) put roughly 20% uncertainty on the
basin dG at 100 ns, so the dipeptides were extended to 400 ns on it. The control is necessary but not sufficient: its
errors cancel by symmetry, which is the error mode of the chiral systems, so only independent replicas at matched
sampling bound the uncertainty. Our own first attempt produced a spurious "handedness is zero" from an unconverged
flat surface (findings 12a).

---

## 7. Implementation notes

The chirality operation on the 72×72 grid (φ,ψ from −180° in 5° steps) is `(φ,ψ) → (−φ,−ψ)`, a reversal with a roll,
because index `i` maps to `(−i) mod 72`:

```python
mirror = lambda m: np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)
```

A plain `m[::-1, ::-1]` is off by one bin and wrong by ~2.3 in practice. This reproduced the first ff3.0's
`rama3.dat` (in git history) from `rama.dat` bit-exactly, which confirmed the convention.

Exactly 4.5% of the coil array is NaN: one whole neighbour column, `CPR`, which is never read, because
`read_rama_maps_and_weights` maps a cis-proline neighbour onto `PRO`. Builders leave it untouched. A naive
`abs(a-b) > tol` diff reports "no difference" across it, because NaN comparisons are False.

The first ff3.0's trainer also forced the GLY row of the rotamer pair-interaction angular profile palindromic, a
second, undocumented symmetrisation; it has been removed (findings 9e).

---

## References

- Ting D, Wang G, Shapovalov M, Mitra R, Jordan MI, Dunbrack RL Jr (2010). *Neighbor-Dependent Ramachandran
  Probability Distributions of Amino Acids Developed from a Hierarchical Dirichlet Process Model.* PLoS Comput Biol
  6(4): e1000763. The library Upside uses (NDRD_TCB).
  https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000763
- Baruch-Shpigler Y, Wang H, Tuvi-Arad I, Avnir D (2017). *Chiral Ramachandran Plots I: Glycine.* Biochemistry
  56(42): 5635–5643. https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525
- Shelar A, Chattopadhyay A (2005). *The Ramachandran plots of glycine and pre-proline.* BMC Struct Biol 5:14.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC1201153/
- *Force field dependent solution properties of glycine oligomers* (2015), PMC4450816: αR/αL approach equality for
  poly-Gly only with sufficient sampling; sampling asymmetry is a known artifact.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4450816/
- Rosetta `-score:symmetric_gly_tables`, off by default, for mixed D/L peptides.
  https://docs.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/D-Amino-Acids ;
  https://docs.rosettacommons.org/docs/latest/structure_prediction/simple_cycpep_predict
- Lovell SC et al. (2003). *Structure validation by Cα geometry: φ,ψ and Cβ deviation.* Proteins 50:437–450.
- Hovmöller S, Zhou T, Ohlson T (2002). *Conformations of amino acids in proteins.* Acta Crystallogr D 58:768–776.
