# Proposed architecture change: a context-resolved glycine Ramachandran term

**Status: PROPOSED, NOT APPROVED, NOT STARTED.** This is a conditional design. The trigger that
would justify building it has not fired, and may not. Section 6 is the decision gate; read it
before writing any code.

Written 2026-09-20 from the ff3.1 glycine investigation. Background and evidence live in
`GLY_sym.md` (the physics) and `findings.md` 9i, 9m-9r (the measurements).

---

## 1. The defect, physically

Glycine's occupancy of **alpha_L** in the PDB is a **placement effect of protein architecture**,
not a property of the residue. Glycine is achiral at C-alpha. L-amino acids build right-handed
elements; connecting them compactly requires occasional left-handed backbone positions; C-beta
sterically forbids positive phi for the other nineteen; so glycine is the residue that gets placed
wherever a fold needs alpha_L.

The statistic NDRD records is therefore **context-specific**, while the rama term that consumes it
is **context-free**. Two failure modes result, pointing opposite ways:

* **At turn glycines the force field pushes twice.** The hydrogen bonding and packing that create
  the turn are computed explicitly by `hbond` and `env`; the map adds the statistical shadow of
  those same interactions on top.
* **At helix glycines it pushes the wrong way.** A glycine inside a helix belongs in alpha_R like
  its neighbours. None of the turn statistics apply, but the alpha_L bias is applied anyway,
  against the helical hydrogen bonding that should hold it.

This is specific to glycine because for the other nineteen, C-beta forbids alpha_L **regardless of
context**: their map value is the same in a turn as in a helix, it is set by local sterics, and it
is genuine local physics that belongs in a local term.

## 2. Why the current architecture cannot express the fix

Sort every ff2.1 term by two properties:

| term | residue-type specific? | (phi,psi) resolved? |
|---|---|---|
| `env` | yes, 8 coefficients per type | no, a function of burial |
| `sidechain` | yes, per-type pair interactions | no, and **absent entirely for GLY** |
| `hbond` | **no**, 3 branch energies shared by all 20 | yes, selected by (phi,psi) |
| **rama map** | **yes** | **yes** |

**The rama map is the only term that is both.** So when the fold's influence on glycine backbone
conformation has to be represented somewhere, the map is the only place it can go. NDRD did not
sloppily include fold effects; the architecture leaves them nowhere else to live.

Two further limits compound it:

* **`hbond` classifies alpha_L as a turn.** Its turn branch is `phi in (0, 165) deg` with no psi
  condition, so a genuine alpha_L helix, a single alpha_L residue in a type II turn, and a glycine
  at phi>0 inside an otherwise right-handed helix all receive the same `E_other = -1.769`, which is
  0.192 E_up weaker than `E_alpha`. A left-handed helix's i->i+4 hydrogen bonds are the exact
  mirror of a right-handed helix's and geometrically identical, so this is a conformational prior,
  not a geometry measurement.
* **The coil/sheet blend is static.** Upside does carry two structural hypotheses per residue,
  blended as `coil_w*exp(-coil) + sheet_w*exp(-(sheet+offset))`, but that mixture is computed at
  config-build time and baked into `rama_pot`. It never updates as the simulation runs.

So Upside distinguishes **where a residue's dihedrals are**, never **what the residue is part of**.
Nothing reads the neighbours' conformations.

## 3. The proposed change

Make the **antisymmetric part** of the glycine map depend on the conformations of its neighbours:

```
E_rama(i) = S(phi_i, psi_i) + w(neighbours of i) * A(phi_i, psi_i)
```

* `S` symmetric and `A` antisymmetric under `(phi,psi) -> (-phi,-psi)`, the parameterisation
  already used by the trainable map (`training/rama_gly_gradient.py`).
* `w` is a smooth scalar in roughly `[0,1]`, near **0** when both neighbours sit in alpha_R
  (helical context) and near **1** otherwise.

Three to five new parameters: the alpha_R window centre and width, and the suppression depth. Not
a new library, not per-neighbour-type maps.

**Why neighbour (phi,psi) and not hydrogen bonding.** An H-bond-driven detector would re-introduce
the double counting in explicit form, since `hbond` already computes exactly that. Neighbour
(phi,psi) is data **no other term reads**, so it is genuinely new information.

## 4. Why NOT a naive context detector

The obvious version of this idea makes de novo folding **worse**, and the reason is worth stating
because it is easy to get backwards.

A turn must form **before** the hydrogen bonds that stabilise it exist. In the unfolded state
nothing pushes a glycine toward alpha_L: no turn geometry, no H-bond partner, no packing. Only the
local term acts. ff2.1's context-free bias applies to every glycine all the time, including in a
fully extended chain, so turns nucleate. **The context-free map is wrong as physics but useful as
a nucleation prior.**

A symmetric "detect the context and use the right value" term switches that off exactly where it
is needed: an extended chain is not a turn, so the detector reports coil and the drive weakens.
**The more precisely you make the map context-aware, the more thoroughly you remove the prior that
gets folding started.**

The proposal above is therefore deliberately **one-sided**: suppress `A` only in helical context,
keep it everywhere else including extended chain. That targets the case where the bias is
demonstrably wrong while preserving nucleation.

## 5. Implementation sketch

Not a specification. Enough to scope the work.

**C++, `src/rama_map_pot.cpp`.** `RamaMapPot::compute_value` currently reads one residue:
`load_vec<2>(ramac, p.residue)`, evaluates one spline layer, and accumulates
`rama_sens(0/1, p.residue)`. It would need to

* read `p.residue - 1` and `p.residue + 1` from the same `rama` CoordNode (already available),
* evaluate **two** layers per residue, `S` and `A`, and return `S + w*A`,
* accumulate sensitivity on the neighbours too, from `dw/dphi_{i±1} * A`.

The node already supports layers via `rama_map_id`, so storing `S` and `A` as separate layers is
natural. Handle chain termini explicitly: a terminal glycine has one neighbour, and `w` must be
defined for that case rather than reading past the array.

**Python, `py/upside_config.py`.** `write_rama_map_pot` would emit `S` and `A` layers for glycine
residues instead of a single pre-mixed map. Everything else keeps one layer with `A = 0`, so
non-glycine behaviour is bit-identical.

**Training, `training/rama_gly_gradient.py`.** The analytic gradient extends directly: `dE/dA`
picks up the factor `w`, and new derivatives appear with respect to `w`'s own parameters. The
existing finite-difference gate (`training/verify_gly_gradient.py`) must be extended to cover
them, and must be run on proteins that exercise helical, turn and terminal glycines. **An analytic
gradient fails silently; the gate is the only thing that catches it.**

**Constraints that still bind.** Master-branch parity for every existing configuration, no guards,
spline tables remain exact representations of their potential, and `GLY|GLY` must stay
mirror-symmetric by construction.

## 6. Decision gate

**Do not build this until the trigger fires.** It is not established that ff2.1 fails on glycine:
in the 32-arm ff3.0 benchmark only the native/de-novo **asymmetry** is significant
(13/16, p = 0.021); neither arm alone reaches p = 0.21, and `hyp_denovo` at +0.166 is an outright
counterexample. Adding machinery to fix an undemonstrated problem is how force fields accumulate
cruft.

Order of operations:

1. Finish Track A training, extract `parameters/ff_3.1_trained`.
2. Run the 32-arm benchmark, scoring on the **last third** of each arm.
3. **Trigger: de novo arms regress again while native arms hold.** That is evidence the fold terms
   cannot nucleate turns on their own, which is exactly what this change addresses. If a corrected
   uniform map holds both arms, the architecture was always sufficient and this document should be
   deleted.

**Cheap measurement that should precede any code.** Split the 456 training natives' glycines by
whether their neighbours are in alpha_R, and compare alpha_L occupancy between the groups. If
helix-context glycines are essentially never in alpha_L while turn-context ones frequently are,
that turns "the architecture conflates two things" into a number and quantifies how much `w` would
have to do. It needs no force-field change and a few minutes of reading structures.

## 7. Risks

* **A new coordinate coupling.** The map would exert forces on neighbouring residues, which the
  optimiser can exploit in ways a purely local term cannot. Watch for `w` drifting to a degenerate
  value that buys energy without physical meaning.
* **One-sidedness is a modelling assumption**, not a derived result. It is chosen to match the
  benchmark's native/de-novo asymmetry, so it must be tested against a prediction made *before*
  seeing the ff3.1 result, not fitted to it afterwards.
* **It only treats glycine.** The same conflation applies to all 20 residue types; for the other
  nineteen the contaminated fraction is small because C-beta sterics dominate, but the argument
  does not stop at glycine and the machinery would generalise. Resist scope creep until glycine is
  demonstrably fixed.
* **Training cost.** The gradient stays analytic, so the cost is bounded, but the gate must be
  re-run and a fresh 500-step run is ~4 days.
