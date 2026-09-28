"""Per-pair basin offsets on the Ramachandran library, trained by matching basin populations.

WHAT IS TRAINED. Every directional map of the library's coil group, k = (central residue,
direction, neighbour), keeps its NDRD values as a fixed base and gets a few smooth basin offsets:

    E_k(phi,psi) = E_k,base(phi,psi) + sum_b c_k,b * w_b(phi,psi)

Each map is renormalised afterwards, as the NDRD maps are (sum exp(-E) = 1), so an offset is a
weight factor on its basin's probability: the shape inside every basin is the base's, and only the
depth of each basin, its frequency, is trained. That is where the fold placement the NDRD statistics
carry shows up (glycine's alpha_L excess, the extra helix). The basins partition the torus, so no
probability can move into an untrained region: alpha_R, alpha_L, beta, pPII and `other` (phi > 0
outside alpha_L) per map, and for a central glycine, which populates `other`, that region split into
its two mirror halves beta' and pPII'.

EACH PAIR IS ITS OWN PARAMETER SET. Offsets are indexed by (central, direction, neighbour) and are
never tied, pooled or shared between maps: the offsets of ALA|GLY act on no map but ALA|GLY. They
are added to the coil entry of their own pair only, never to the sheet group, because a central
cis-proline reads PRO's sheet entry and an offset there would act on both. The one constraint lives
inside single maps: GLY|GLY has an achiral pair, so its base is symmetrised and each of its offsets
is held equal to its mirror basin's. Its sheet entry is symmetrised too, since upside_config mixes
every residue's coil map with its sheet map, and NDRD's GLY|GLY sheet map holds 94% of its weight at
phi < 0. Symmetrising averages a map's probabilities with its mirror's, i.e. pools every site with
its mirror image.

HOW THEY ARE UPDATED. Once per training round (one epoch), per map and basin, the basin population
of the native-restrained replica is compared with that of the free replicas over every residue that
reads the map. Each offset takes a damped Newton step on the maximum a posteriori objective: the
native basin counts under the model's basin populations, with a Gaussian prior of width `SIGMA`
on the offset, centred on zero, i.e. on the map's own NDRD values:

    c += eta * [T0 N (p_free - p_native) - T0^2 c / sigma^2] / [N p (1 - p) + T0^2 / sigma^2]

with N the residues reading the map and p the mean of the two populations. A basin the free
simulation over-populates is raised, one it under-populates is lowered. Where a basin holds many
residues this is the log-ratio step T0 ln(p_free / p_native); where it holds almost none, the
prior bounds the step and an offset with no evidence decays back to zero. A plain log-ratio of
two near-zero populations is counting noise and gave steps of 1.8 nats in forbidden basins in the
first round of the first run. Nothing is borrowed from another map. There is no unfolded-state
(DSE) term on the offsets: in ConDiv the SARW replica keeps the rama term, so that comparison would
be against the map itself, not against data.
"""

import os

import numpy as np
import tables as tb

BASINS = ('alpha_R', 'alpha_L', 'beta', 'pPII', "beta'", "pPII'", 'other')
N_BASIN = len(BASINS)
# the basin each basin maps onto under (phi,psi) -> (-phi,-psi); `other` is active only on maps
# without a central glycine, where no mirror tie applies
MIRROR = (1, 0, 4, 5, 2, 3, 6)
# Logistic edge scale. At 3 deg the 10-90% transition is 13 deg, two to three grid cells, the
# sharpest the engine's spline represents cleanly; an offset is then realised at a median 98% of
# its value where a basin's probability lies, and 8% of a map's probability sits in transitions.
# The 10 deg scale of secstr_bias gives 44 deg edges, which leave no flat interior in an 80 deg
# basin and tilt its shape instead of scaling it.
EDGE = np.deg2rad(3.)
SIGMA = 1.                  # prior width on every offset, nats
ETA = 0.5


# ---------------------------------------------------------------------------
# Basins
# ---------------------------------------------------------------------------

def _wrap(x):
    """Angle difference into [-pi, pi)."""
    return (np.asarray(x) + np.pi) % (2. * np.pi) - np.pi


def _sigmoid(x):
    return 1. / (1. + np.exp(-x))


def _arc(x, a, b):
    """Smooth periodic indicator of the arc from a to b (radians, anticlockwise).

    A sigmoid of the circular distance from the arc's midpoint, so it is continuous everywhere,
    including across phi = +-180, and an arc mirrors exactly onto the arc (-b, -a).
    """
    half = 0.5 * ((b - a) % (2. * np.pi))
    return _sigmoid((half - np.abs(_wrap(x - (a + half)))) / EDGE)


def basin_weights(phi, psi):
    """w_b(phi,psi) for every basin, shape (..., N_BASIN). Angles in radians."""
    d = np.deg2rad
    phi, psi = np.asarray(phi, dtype=float), np.asarray(psi, dtype=float)
    helix_R = _arc(psi, d(-100.), d(50.))       # alpha_R's psi band and its complement
    helix_L = _arc(psi, d(-50.), d(100.))       # the mirror band
    beta_m = _arc(phi, d(100.), np.pi) * (1. - helix_L)
    ppii_m = _arc(phi, 0., d(100.)) * (1. - helix_L)
    return np.stack([
        _arc(phi, -np.pi, 0.) * helix_R,                  # alpha_R
        _arc(phi, 0., np.pi) * helix_L,                   # alpha_L
        _arc(phi, -np.pi, d(-100.)) * (1. - helix_R),     # beta
        _arc(phi, d(-100.), 0.) * (1. - helix_R),         # pPII
        beta_m,                                           # beta'
        ppii_m,                                           # pPII'
        beta_m + ppii_m,                                  # other
    ], axis=-1)


def grid_basins(n_grid):
    """Basin weights on the library grid: node i sits at -180 + i * 360/n_grid degrees."""
    t = -np.pi + 2. * np.pi * np.arange(n_grid) / n_grid
    phi, psi = np.meshgrid(t, t, indexing='ij')
    return basin_weights(phi, psi)                         # (n_grid, n_grid, N_BASIN)


def mirror(m):
    """(phi,psi) -> (-phi,-psi) on the grid: a reversal WITH a roll, because node i maps to (-i) % n."""
    return np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)


def symmetrise(m):
    """The mirror-symmetric map whose probabilities are the mean of m's and its mirror's.

    Averaging the energies instead would take the geometric mean of the probabilities, which
    empties any basin whose mirror is empty (the sheet map's pPII). Exact: every point and its
    mirror are given the same pair of numbers, and logaddexp is symmetric in them.
    """
    return np.log(2.) - np.logaddexp(-m, -mirror(m))


# ---------------------------------------------------------------------------
# The parameter set: one offset vector per directional map
# ---------------------------------------------------------------------------

def _decode(attr):
    return [x.decode() if isinstance(x, bytes) else x for x in attr]


_base_cache = {}


def base_maps(state):
    """The fixed base of every key: the source library's map, GLY|GLY symmetrised.

    Rebuilt from the source library rather than kept in the state, which goes into every
    checkpoint; cached per library.
    """
    src = state['source']
    if src not in _base_cache:
        with tb.open_file(src) as t:
            pot = t.root.coil.dimer_pot[:].astype(float)
        base = np.stack([pot[c, d, n] for c, d, n in state['keys']])
        g = state['is_gg']
        base[g] = symmetrise(base[g])
        _base_cache[src] = base
    return _base_cache[src]


def init_state(library):
    """Offsets at zero on the library's own maps, GLY|GLY symmetrised.

    Keys are every (central, direction, neighbour) that `read_rama_maps_and_weights` can read in
    mixture mode: every central residue type of the coil group (cis-proline included), both
    directions, and the 20 standard neighbours (a cis-proline neighbour is read as PRO, and ALL is
    used only by the product rule).
    """
    with tb.open_file(library) as t:
        coil = t.root.coil
        restype = _decode(coil._v_attrs.restype)
        dirs = _decode(coil._v_attrs.dir)
        pot = coil.dimer_pot[:]
    neighbours = [n for n, r in enumerate(restype) if r not in ('ALL', 'CPR')]
    keys = np.array([(c, d, n) for c in range(pot.shape[0]) for d in range(len(dirs))
                     for n in neighbours if np.isfinite(pot[c, d, n]).all()], dtype=int)
    gly = restype.index('GLY')
    is_gg = (keys[:, 0] == gly) & (keys[:, 2] == gly)
    state = dict(source=os.path.abspath(library), keys=keys, restype=restype, dirs=dirs,
                 is_gg=is_gg)
    base = base_maps(state)

    active = np.zeros((len(keys), N_BASIN), dtype=bool)
    active[:, :4] = True
    central_gly = keys[:, 0] == gly
    active[central_gly, 4:6] = True
    active[~central_gly, 6] = True
    state.update(active=active, offset=np.zeros((len(keys), N_BASIN)), round=0, history=[])
    return state


def independent(state):
    """Offsets that are free parameters: the active ones, counting each GLY|GLY mirror pair once."""
    sel = state['active'].copy()
    sel[np.ix_(state['is_gg'], [1, 4, 5])] = False
    return sel


def n_param(state):
    return int(independent(state).sum())


def _normalise(m):
    """The library's convention: sum(exp(-E)) == 1 per map."""
    lo = m.min(axis=(-2, -1), keepdims=True)
    return m + np.log(np.exp(-(m - lo)).sum(axis=(-2, -1), keepdims=True)) - lo


def maps(state):
    """Every key's current map, shape (n_key, n_grid, n_grid), normalised like the NDRD maps.

    Renormalising makes each offset a weight factor on its basin: probability moves between the
    basins of one map, never between a map and its partner in upside_config's left/right mixture,
    which mixes the maps without renormalising them first. GLY|GLY is symmetric by construction
    (symmetric base, tied offsets on mirror basins); the projection only removes summation-order
    rounding, so the written maps are exactly mirror-symmetric.
    """
    base = base_maps(state)
    m = base + np.einsum('kb,ijb->kij', state['offset'], grid_basins(base.shape[-1]))
    g = state['is_gg']
    m[g] = 0.5 * (m[g] + mirror(m[g]))
    return _normalise(m)


def write_library(state, source, out):
    """Copy `source`, replace each key's coil entry with its current, normalised map and
    symmetrise the GLY|GLY sheet entries.

    Nothing else changes: the rest of the sheet group, the ALL and cis-proline neighbour columns and
    both weight arrays come through as they were. Only a central glycine next to a glycine reads the
    GLY|GLY sheet entry (a cis-proline reads PRO's), so symmetrising it acts on no other residue.
    """
    import shutil
    shutil.copy(source, out)
    m = maps(state)
    with tb.open_file(out, 'a') as t:
        pot = t.root.coil.dimer_pot[:]
        for k, (c, d, n) in enumerate(state['keys']):
            pot[c, d, n] = m[k]
        t.root.coil.dimer_pot[:] = pot
        sheet = t.root.sheet.dimer_pot[:]
        g = _decode(t.root.sheet._v_attrs.restype).index('GLY')
        sheet[g, :, g] = _normalise(symmetrise(sheet[g, :, g].astype(float)))
        t.root.sheet.dimer_pot[:] = sheet
    return out


def gly_gly_asymmetry(library):
    """Largest |E - mirror(E)| over the GLY|GLY coil and sheet entries of a library file."""
    with tb.open_file(library) as t:
        asym = []
        for grp in (t.root.coil, t.root.sheet):
            g = _decode(grp._v_attrs.restype).index('GLY')
            m = grp.dimer_pot[:][g, :, g].astype(float)
            asym.append(np.abs(m - mirror(m)).max())
    return max(asym)


def residue_keys(state, seq):
    """For each residue, the key index of every coil map it reads, as `read_rama_maps_and_weights`.

    A terminal residue reads one direction, an interior one both. A cis-proline is CPR only as the
    central residue; as a neighbour it is read as PRO.
    """
    ridx = {r: i for i, r in enumerate(state['restype'])}
    didx = {d: i for i, d in enumerate(state['dirs'])}
    kidx = {tuple(k): i for i, k in enumerate(state['keys'])}

    def key(c, d, n):
        return kidx[(ridx[c], didx[d], ridx['PRO' if n == 'CPR' else n])]

    seq = list(seq)
    out = [[key(seq[0], 'right', seq[1])]]
    for i in range(1, len(seq) - 1):
        out.append([key(seq[i], 'left', seq[i - 1]), key(seq[i], 'right', seq[i + 1])])
    out.append([key(seq[-1], 'left', seq[-2])])
    return out


# ---------------------------------------------------------------------------
# Statistics and the update
# ---------------------------------------------------------------------------

def residue_populations(rama_coord, weights=None):
    """Per-residue basin populations of one ensemble.

    rama_coord: (n_frame, n_res, 2), (phi, psi) in radians as the engine outputs them.
    weights: per-frame weights summing to one, or None for a plain mean.
    Returns (n_res, N_BASIN).
    """
    w = basin_weights(rama_coord[..., 0], rama_coord[..., 1])
    if weights is None:
        return w.mean(axis=0)
    return np.tensordot(np.asarray(weights, dtype=float), w, axes=1)


def empty_accumulator(state):
    n = len(state['keys'])
    return dict(n=np.zeros(n), free=np.zeros((n, N_BASIN)), native=np.zeros((n, N_BASIN)))


def accumulate(acc, state, seq, free, native):
    """Add one protein's per-residue populations to every map its residues read."""
    for i, ks in enumerate(residue_keys(state, seq)):
        for k in ks:
            acc['n'][k] += 1.
            acc['free'][k] += free[i]
            acc['native'][k] += native[i]


def _tie(state, x):
    """Within each GLY|GLY map, give a basin and its mirror the sum over both, i.e. the union."""
    x = np.array(x, dtype=float)
    g = state['is_gg']
    x[g] = x[g] + x[g][:, MIRROR]
    return x


def step_gradient(state, acc):
    """One step's share of the objective's gradient, in residue counts, for the convergence gate.

    The data part, free minus native counts, plus this step's share of the prior's pull, so that
    at the fixed point the steps of an epoch are noise about zero.
    """
    pull = state['T0'] * state['offset'] / SIGMA ** 2 / state['steps_per_round']
    return (_tie(state, acc['free'] - acc['native']) - pull)[independent(state)]


def update(state, acc):
    """One round's damped Newton step on every offset. Returns the step taken, (n_key, N_BASIN)."""
    T0 = state['T0']
    n = acc['n'][:, None]
    free, native = _tie(state, acc['free']), _tie(state, acc['native'])
    p = 0.5 * (free + native) / np.maximum(n, 1.)
    grad = T0 * (free - native) - T0 ** 2 * state['offset'] / SIGMA ** 2
    hess = n * p * (1. - p) + T0 ** 2 / SIGMA ** 2
    step = np.where(state['active'], ETA * grad / hess, 0.)
    state['offset'] = state['offset'] + step
    state['round'] += 1
    return step


def summary(state, acc_train, acc_held):
    """A round's mismatch at the offsets it ran with: the fraction of residue time in a different
    basin between the free and native ensembles (total variation), averaged over residues."""
    def mismatch(acc):
        d = 0.5 * np.abs(np.where(state['active'], acc['free'] - acc['native'], 0.)).sum(axis=1)
        return float(d.sum() / max(acc['n'].sum(), 1.))
    return dict(sites_train=float(acc_train['n'].sum()), sites_held=float(acc_held['n'].sum()),
                mismatch_train=mismatch(acc_train), mismatch_held=mismatch(acc_held))


def offset_summary(state):
    """Largest offset, and the mean alpha_L minus alpha_R offset over the GLY|X maps (central
    glycine, neighbour not glycine)."""
    gly = state['restype'].index('GLY')
    xg = (state['keys'][:, 0] == gly) & ~state['is_gg']
    return dict(max_offset=float(np.abs(state['offset']).max()),
                gly_dL_minus_dR=float((state['offset'][xg, 1] - state['offset'][xg, 0]).mean()))
