# ConDiv_symlay Training Potential Design

This note records the design principles for the `ConDiv_symlay` membrane-training workflow as currently implemented.

## Objective

The goal is to train an implicit membrane potential that emulates a fixed DOPC bilayer while keeping the training surface compatible with the existing Upside membrane machinery.

The design therefore uses:

1. the existing implicit membrane spline potentials as the trainable runtime model
2. a hard projection into a symmetric DOPC topology-slot subspace
3. a soft dry-MARTINI-based teacher that is distance-resolved and exportable as explicit pair tables

## Runtime Potential

The trainable runtime membrane terms remain the standard Upside terms:

- `cb_membrane_potential`
- `hb_membrane_potential`

The optimizer updates the outer spline coefficient blocks:

- `cb` with shape `(20, 2, 18)`
- `hb` with shape `(2, 2, 18)`

The inner blocks:

- `icb`
- `ihb`

are still carried in the parameter object and written to `membrane.h5`, but the direct MD differentiation path in this workflow is through the standard outer membrane nodes.

Conceptually:

```text
E_CB(r) = sum_b g_b(burial_r) * f_aa,b(z_r)

E_HB(v) = (1 - p_v)^2 * f_type,0(z_v)
        + (1 - (1 - p_v)^2) * f_type,1(z_v)
```

where `z` is the residue or hydrogen-bond depth in the membrane frame.

## Total Training Update

The training loop keeps the original ConDiv divergence gradient and adds the membrane regularizer gradient:

```text
g_total = g_ConDiv + g_pair
theta_tmp = theta + Adam(g_total)
theta_next = P(theta_tmp)
```

where:

- `g_ConDiv` is the original restrained-vs-free trajectory contrast gradient
- `g_pair` is the dry-MARTINI teacher regularization gradient
- `P` is the hard symmetric-layer projector

So the soft constraint influences the step, but the hard membrane subspace is enforced after every update.

## Hard Constraint: Symmetric DOPC Layer Projection

The hard constraint is a projection onto an even, mirrored, layer-anchored spline family derived from the DOPC slot manifest.

For each membrane curve `f(z)`:

```text
1. choose symmetric support
   Z = ceil(max(|z_min|, |z_max|, max_slot + margin))
   support = [-Z, Z]

2. even symmetrize
   f_even(z) = 0.5 * (f(z) + f(-z))

3. define DOPC anchors
   a = [0, d_1, d_2, ..., d_K, Z]

4. assign anchor values by dense-grid Voronoi averaging

5. reconstruct an even piecewise-linear target over |z|

6. solve back to clamped-spline coefficients
```

This enforces:

- even symmetry in `z`
- symmetric support
- a shape tied to the DOPC topology-slot depths instead of arbitrary free-form spline variation

The full constrained slot sequence is:

```text
Q0-Qa-Na-C1-C3-C1-C1-C1-C1-C3-C1-Na-Qa-Q0
```

Repeated dry types are preserved as distinct positional layers in the hard projection.

## Soft Constraint: Distance-Resolved Dry-MARTINI Teacher

The soft constraint is not a direct replacement for the implicit membrane objective. It is a regularizer that forces the trained implicit membrane to stay consistent with a dry-MARTINI-inspired explicit-bilayer picture.

### Covered DOPC Anchor Types

The regularizer first pools the repeated DOPC slots back to the covered dry types:

```text
Q0, Qa, Na, C1, C3
```

Pooling uses the topology multiplicity carried by the manifest, so repeated slots remain properly weighted before being collapsed to these covered anchor types.

### Reference Anchor Values

At initialization, the seeded implicit membrane is sampled at the DOPC slot depths. Those sampled values are pooled to the covered dry types and treated as the reference anchor values for three channels:

- `cb_mean`
- `hb_donor`
- `hb_acceptor`

These are fixed for one training run.

### All-Type Radial Pair Tables

For each channel `a` and dry-MARTINI type `t`, the workflow constructs a radial pair table:

```text
V_a,t(r)
```

The spanning rule is role-specific:

- `cb_mean` uses the mean of dry proxy rows `C1` and `C2`
- `hb_donor` uses dry proxy row `Qd`
- `hb_acceptor` uses dry proxy row `Qa`

The table shape is built from:

1. a dry-MARTINI descriptor from the role-specific pair `epsilon`
2. a distance scale from the role-specific pair `sigma`
3. interpolation from the covered DOPC anchor types to all dry types
4. a smooth radial control-point curve on a fixed radial grid

So the teacher spans all dry types while retaining distance dependence.

### Explicit Bilayer Projection

Those radial tables are projected through an explicit DOPC shell kernel measured from the template bilayer:

```text
s_a(c) = sum_t sum_m K[c,t,m] * V_a,t(r_m)
```

where:

- `c` is one covered DOPC type
- `t` runs over all dry-MARTINI types
- `r_m` is the radial grid
- `K[c,t,m]` is the explicit bilayer shell kernel

This gives a fixed teacher score `s_a(c)` for each covered DOPC type.

### Loss Form

During training, the current implicit membrane is sampled at the slot depths, pooled back to the covered types, and compared to the teacher trend with the same affine residual form as the older regularizer:

```text
L_pair^a = 0.5 * sum_c w_c * ( y_a(c) - (beta0_a + beta1_a * s_a(c)) )^2
```

where:

- `y_a(c)` is the current pooled implicit membrane value
- `s_a(c)` is the fixed projected teacher score
- `w_c` is the covered-type weight
- `beta0_a, beta1_a` are fitted affine trend parameters for that channel

Total pair loss:

```text
L_pair = lambda_cb * L_pair^cb
       + lambda_donor * L_pair^donor
       + lambda_acceptor * L_pair^acceptor
```

with the workflow defaults:

```text
lambda_cb = 0.05
lambda_donor = 0.10
lambda_acceptor = 0.10
```

This keeps the regularizer on the membrane mean channels while allowing a direct export path to explicit pair tables through `V_a,t(r)`.

## Bilayer Geometry Assumptions

### Fixed Bilayer Species

The workflow is designed to emulate one fixed DOPC bilayer, not a target-specific slab.

Accordingly:

- membrane thickness is one global workflow value
- the current default is `30.2 A`
- per-target `.thickness` files are not used by the active `ConDiv_symlay` training workflow

### Protein Orientation Relative to Bilayer

The protein-vs-bilayer frame is assumed to already be encoded in the stored training coordinates.

That means:

- the workflow does not infer a new membrane alignment per target
- the native `*.initial.pkl` coordinates are used directly as the initial structure
- z recentering is disabled during simulation so the membrane frame is preserved

Conceptually:

```text
bilayer midplane = z = 0
bilayer half-thickness = H / 2, with one fixed DOPC thickness H
protein depth relative to bilayer = z already present in the training coordinates
```

So the design intentionally separates:

- fixed DOPC slab geometry
- target-specific protein placement within that slab frame

## Why This Design

This combination was chosen to satisfy four constraints at once:

1. keep the runtime training surface stable by staying on the existing implicit membrane nodes
2. enforce DOPC symmetry and layering exactly through a hard projection
3. make the soft teacher physically span all dry-MARTINI types with distance dependence
4. keep a conversion path from the learned implicit membrane back to explicit pair tables

In short:

```text
train implicit membrane splines directly
constrain them hard to a DOPC layered subspace
regularize them softly against a dry-MARTINI explicit-bilayer teacher
preserve the protein-bilayer frame already present in the data
```
