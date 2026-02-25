# Infrastructure Plan: Hybrid-Cycle to AA-Backbone-in-Explicit-Lipid Runtime

## Summary
Implement Upside infrastructure so:
1. `MARTINI protein` is used only to generate/relax initial training systems.
2. Runtime physics for target simulations is driven by **Upside all-atom backbone + explicit lipid bilayer** interactions from a trained **RBM radial FF**.
3. All trained FF parameters are stored in a single file: `explicit_bilayer.h5`.

## Core runtime architecture change
Add two explicit operating modes in hybrid control:

1. `bootstrap_martini_protein` (training prep mode)
- MARTINI protein beads present.
- protein rigid; lipids relax.
- used only for generating stable initial states.

2. `aa_backbone_explicit_lipid` (target mode)
- MARTINI protein beads are retained only as optional inert carriers for compatibility, but **excluded from force/energy/integration**.
- active protein degrees of freedom are Upside AA carriers/backbone atoms.
- all protein-lipid interactions come from trained RBM cross potential.

This removes dependence on MARTINI protein interaction physics in the final model.

## New potential infrastructure in Upside
Implement a new node in `src/martini.cpp` (or split file), e.g. `martini_rbm_cross_potential`:

1. Inputs
- protein interaction sites:
  - BB classes: `(AA, N/CA/C/O)` (80 classes)
  - SC classes (existing Upside SC representation)
- explicit lipid/environment MARTINI bead classes.
- radial basis features (bins/spline basis) per pair.

2. Energy
- RBM-style radial coupling (Upside-style learnable node, not LJ/Coulomb table).

3. Parameter API
- fully implement `get_param`, `set_param`, `get_param_deriv`.
- deterministic flatten/unflatten order.
- gradient-compatible with existing `engine.get_param_deriv`.

4. Gating/filtering in integration cycle
- when mode = `aa_backbone_explicit_lipid`, skip MARTINI protein-protein and MARTINI protein-lipid pair contributions entirely.
- keep lipid-lipid MARTINI interactions.
- keep Upside internal protein terms.
- apply RBM cross node for protein-lipid interactions.

## `explicit_bilayer.h5` as single artifact
Store only FF/state needed for reuse:

1. `/meta`
2. `/classes` (`martini_types`, `sc_types`, `bb_types`, radial basis metadata)
3. `/rbm/sc` (W, hidden biases, optional visible biases)
4. `/rbm/bb` (W, hidden biases, optional visible biases)
5. `/views`
- explicit effective tables:
  - MARTINI type vs SC type vs radial basis
  - MARTINI type vs BB `(AA,N/CA/C/O)` vs radial basis
6. `/priors` and bounds
7. optional `/system_maps` for reproducible class indexing

## Training workflow (unchanged objective, updated infrastructure use)
1. Read IDs from `/Users/yinhan/Documents/Train/upside_input/list`.
2. Download OPM-oriented proteins.
3. Build dry-MARTINI protein and pack in dry-MARTINI bilayer.
4. Run rigid-protein lipid relaxation (`bootstrap_martini_protein`).
5. Switch to training configs where contrasts/derivatives are taken from `martini_rbm_cross_potential`.
6. CD updates write directly into `explicit_bilayer.h5`.

## Integration-cycle implementation tasks
1. Add mode switch handling in `hybrid_control` parsing and stage update logic.
2. Add atom activity masks:
- `active_protein_atoms` (Upside AA backbone/SC carriers)
- `inactive_martini_protein_atoms` in target mode.
3. In MARTINI pair loop:
- hard-skip inactive MARTINI protein atoms.
4. In force projection / BB refresh paths:
- disable MARTINI-proxy refresh/projection when target mode is active.
5. Ensure logging reports mode and active atom counts each stage.

## Validation/test plan
1. Unit
- RBM node param roundtrip and FD gradient checks.
- mode gating tests (pairs correctly skipped/kept by mode).
2. Integration
- bootstrap run relaxes lipids with rigid MARTINI protein.
- target mode run shows no MARTINI-protein force contribution.
- RBM cross term contributes finite energy/forces and updates in training.
3. Regression
- old MARTINI workflows unchanged when new mode/node disabled.

## Assumptions locked
1. Final production physics should not depend on MARTINI protein interaction terms.
2. MARTINI protein may remain in files as inert compatibility carriers unless explicitly stripped.
3. BB granularity remains AA-specific `N/CA/C/O`.
4. Output remains single-file `explicit_bilayer.h5`.
