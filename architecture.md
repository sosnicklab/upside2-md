# Proposed architecture change: a context-resolved glycine Ramachandran term

**Status: PROPOSED, NOT APPROVED, NOT STARTED, and superseded in part.** plan.md Phase 8 adds a
smaller context term: glycine-specific offsets on the three H-bond basin energies, trained by
ConDiv, which read the glycine's own H-bond state and basin rather than its neighbours' (phi, psi).
This document is the larger fallback if that fails. Section 6 is the decision gate; read it before
writing any code.

Written 2026-09-20, updated 2026-10-02. The physics is in `GLY_sym.md`, the measurements in
`findings.md` 1.15-1.17, and the retired learned-map work in findings 9e.

---

## 1. The defect, physically

Glycine's occupancy of alpha_L in the PDB is a placement effect of protein architecture, not a
property of the residue. Glycine is achiral at C-alpha. L-amino acids build right-handed elements;
connecting them compactly requires occasional left-handed backbone positions; C-beta sterically
forbids positive phi for the other nineteen; so glycine is the residue that gets placed wherever a
fold needs alpha_L.

The statistic NDRD records is therefore context-specific, while the rama term that consumes it is
context-free. Two failure modes result, pointing opposite ways:

* **At turn glycines the force field pushes twice.** The hydrogen bonding and packing that create
  the turn are computed explicitly by `hbond` and `env`; the map adds the statistical shadow of
  those same interactions on top.
* **At helix glycines it pushes the wrong way.** A glycine inside a helix belongs in alpha_R like
  its neighbours. None of the turn statistics apply, but the alpha_L bias is applied anyway,
  against the helical hydrogen bonding that should hold it.

This is specific to glycine because for the other nineteen, C-beta forbids alpha_L regardless of
context: their map value is the same in a turn as in a helix, it is set by local sterics, and it is
genuine local physics that belongs in a local term.

## 2. Why the current architecture cannot express the fix

Sort every ff2.1 term by two properties:

| term | residue-type specific? | (phi,psi) resolved? |
|---|---|---|
| `env` | yes, 8 coefficients per type | no, a function of burial |
| `sidechain` | yes, per-type pair interactions | no, and absent entirely for GLY |
| `hbond` | no, 3 branch energies shared by all 20 | yes, selected by (phi,psi) |
| rama map | yes | yes |

**The rama map is the only term that is both.** So when the fold's influence on glycine backbone
conformation has to be represented somewhere, the map is the only place it can go. NDRD did not
sloppily include fold effects; the architecture leaves them nowhere else to live.

Two further limits compound it:

* **`hbond` classifies alpha_L as a turn.** Its turn branch is `phi in (0, 165) deg` with no psi
  condition, so a genuine alpha_L helix, a single alpha_L residue in a type II turn, and a glycine
  at phi>0 inside an otherwise right-handed helix all receive the same `E_other` (ff2.1: -1.769,
  0.192 E_up weaker than `E_alpha`). A left-handed helix's i->i+4 hydrogen bonds are the exact
  mirror of a right-handed helix's and geometrically identical, so this is a conformational prior,
  not a geometry measurement. Phase 8's offsets give glycine its own three branch energies, which
  makes `hbond` residue-type specific for glycine, but still blind to what the residue is part of.
* **The coil/sheet blend is static.** Upside does carry two structural hypotheses per residue,
  blended as `coil_w*exp(-coil) + sheet_w*exp(-(sheet+offset))`, but that mixture is computed at
  config-build time and baked into `rama_pot`. It never updates as the simulation runs.

So Upside distinguishes where a residue's dihedrals are, never what the residue is part of.
Nothing reads the neighbours' conformations.

## 3. The proposed change

Make the antisymmetric part of the glycine map depend on the conformations of its neighbours:

```
E_rama(i) = S(phi_i, psi_i) + w(neighbours of i) * A(phi_i, psi_i)
```

* `S` symmetric and `A` antisymmetric under `(phi,psi) -> (-phi,-psi)`: the mirror-symmetric
  and antisymmetric projections of the measured `GLY|X` map (`parameters/common/rama31.dat`). The
  BioEmu-fitted surface (since 2026-10-05) is moderately asymmetric (ln(aR/aL) about -0.55 in
  Upside on the octapeptides); the AWH dipeptide surface it replaced was mildly so (-0.12).
* `w` is a smooth scalar in roughly `[0,1]`, near 0 when both neighbours sit in alpha_R
  (helical context) and near 1 otherwise.

Three to five new parameters: the alpha_R window centre and width, and the suppression depth. Not
a new library, not per-neighbour-type maps.

**Why neighbour (phi,psi) and not hydrogen bonding.** An H-bond-driven detector re-introduces the
double counting in explicit form, since `hbond` already computes exactly that. Neighbour (phi,psi)
is data no other term reads, so it is new information. Phase 8 tries the H-bond route first
because it needs no new coupling between residues.

## 4. Why NOT a naive context detector

The obvious version of this idea makes de novo folding worse, and the reason is easy to get
backwards.

A turn must form before the hydrogen bonds that stabilise it exist. In the unfolded state nothing
pushes a glycine toward alpha_L: no turn geometry, no H-bond partner, no packing. Only the local
term acts. ff2.1's context-free bias applies to every glycine all the time, including in a fully
extended chain, so turns nucleate. The context-free map is wrong as physics but useful as a
nucleation prior.

A symmetric "detect the context and use the right value" term switches that off exactly where it
is needed: an extended chain is not a turn, so the detector reports coil and the drive weakens.
The more precisely the map is made context-aware, the more thoroughly it removes the prior that
gets folding started. The native panel already shows the non-local terms' limit in a folded
setting: with the AWH map and ff2.1 untrained, natively left-handed glycines lose about 9 points of
alpha_L (findings 1.17).

The proposal above is therefore deliberately one-sided: suppress `A` only in helical context, keep
it everywhere else including extended chain. That targets the case where the bias is demonstrably
wrong while preserving nucleation.

## 5. Implementation sketch

Not a specification. Enough to scope the work.

**C++, `src/rama_map_pot.cpp`.** `RamaMapPot::compute_value` currently reads one residue:
`load_vec<2>(ramac, p.residue)`, evaluates one spline layer, and accumulates
`rama_sens(0/1, p.residue)`. It would need to

* read `p.residue - 1` and `p.residue + 1` from the same `rama` CoordNode (already available),
* evaluate two layers per residue, `S` and `A`, and return `S + w*A`,
* accumulate sensitivity on the neighbours too, from `dw/dphi_{i±1} * A`.

The node already supports layers via `rama_map_id`, so storing `S` and `A` as separate layers is
natural. Handle chain termini explicitly: a terminal glycine has one neighbour, and `w` must be
defined for that case rather than reading past the array.

**Python, `py/upside_config.py`.** `write_rama_map_pot` would emit `S` and `A` layers for glycine
residues instead of a single pre-mixed map. Everything else keeps one layer with `A = 0`, so
non-glycine behaviour is bit-identical.

**Training.** The glycine map itself stays fixed; only `w`'s parameters would be trained, through
an analytic derivative in ConDiv. A gate would have to show that `w` reaches exactly the glycines
whose neighbours it reads, on proteins that exercise helical, turn and terminal glycines. An
indexing slip fails silently; the gate is the only thing that catches it.

**Constraints that still bind.** Master-branch parity for every existing configuration, no guards,
spline tables remain exact representations of their potential, and `GLY|GLY` must stay
mirror-symmetric by construction.

## 6. Decision gate

Do not build this unless the smaller term fails. Order of operations:

1. Train ff3.0 as plan.md Phase 8 specifies (`ff30_glyhb`): ff2.1's workflow from ff2.1, the AWH
   glycine map fixed, glycine's H-bond offsets trained, the side-chain step damped.
2. Validate it: the selection panel at every epoch end (helical and natively left-handed glycines
   against all-atom), lambda's helix 3 (wild type against G46A/G48A), glpG TM4's helical glycines,
   and the Peng benchmark paired against ff2.1.
3. **Trigger: helical glycines still lose their helix, or natively left-handed glycines still lose
   alpha_L, with the H-bond offsets trained.** That says the glycine's own H-bond state is not
   enough context, which is what neighbour (phi,psi) adds. If both classes pass, delete this
   document.

The 32-arm benchmark statistics that first motivated this concerned the retired symmetrised ff3.0
(findings 9d).

**What is measured already, and what is not.** Native glycines are H-bonded in both basins, in
different `hbond` branches, and 62-66% of the glycine-specific helical loss sits in glycines with
their own H-bond (findings 1.16), which is why Phase 8 tries the H-bond term first. If the trigger
fires, split the remaining glycine loss by whether the neighbours are in alpha_R before writing any
code: that number says how much `w` would have to do.

## 7. Risks

* **A new coordinate coupling.** The map would exert forces on neighbouring residues, which the
  optimiser can exploit in ways a purely local term cannot. Watch for `w` drifting to a degenerate
  value that buys energy without physical meaning.
* **One-sidedness is a modelling assumption**, not a derived result. Test it against a prediction
  made before seeing the trained result, not fitted to it afterwards.
* **It only treats glycine.** The same conflation applies to all 20 residue types; for the other
  nineteen the contaminated fraction is small because C-beta sterics dominate, but the argument
  does not stop at glycine and the machinery would generalise. Resist scope creep until glycine is
  demonstrably fixed.
* **Training cost.** The gradient stays analytic, so the cost is the gate plus a fresh run of the
  FF2 workflow.
