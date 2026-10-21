"""Every coil map's basin depths on the ff2.1 Ramachandran library, trained by matching basin
populations under a Gaussian prior centred on NDRD (ff3.1; repo plan.md Phases 14 and 29, findings
1.28-1.29, 1.34).

THE PARAMETERS. Every coil map the engine reads (central residue, cis-proline included; direction;
neighbour, a cis-proline neighbour read as PRO, as `read_rama_maps_and_weights` does) keeps its NDRD
values as a fixed base and gets one offset per basin b of the six below:

    E_m(phi,psi) = E_m,base(phi,psi) + sum_b c_mb w_b(phi,psi),    c_mb = q_X(m),b + d_mb

renormalised afterwards (sum exp(-E) = 1, the library's convention). An offset is a weight factor
on its basin's probability; the shape inside every basin stays NDRD's. d_mb is the map's own
deviation; q_Xb is shared by the maps of central residue X, and only in mode 'pooled'. The sheet
group is written unchanged.

THREE RUNS, one file (MODE below):
  'all'     every deviation, prior width SIGMA_ALL (the cross-validated width, findings 1.29);
  'data'    only the deviations with at least MIN_NATIVE native sites in the training proteins'
            native structures, at the narrower width SIGMA_DATA; the rest stay at NDRD;
  'pooled'  'data''s deviations, plus one shared offset per (central residue, basin) where that
            residue has MIN_NATIVE native sites in the basin, at width SIGMA_POOL (cross-validated,
            plan.md Phase 29). The training proteins pull a map's basin through each residue that
            reads it, a pull that is consistent only summed over a residue's maps (findings 1.34),
            so in 'data' each map's own prior holds it at the start. A map stays its own conditional
            (findings rule 1.8): the shared offset is an estimate from the residue's own maps, as
            NDRD's hierarchical prior pulls a rare pair toward its pooled map, never one residue's
            correction applied to another; a map with enough data departs from it through d_mb.
            Pre-proline maps (X|right|PRO) keep only their deviations, as ff3.0's glycine pair
            left them out (findings 1.14).
Modes 'all' and 'data' have no shared offsets (q = 0), so each map has only its own prior.

NO OUTSIDE DATA IN THE TRAINING. The base is ff2.1's own library and the offsets move only by the
training set; they start at zero, i.e. at NDRD.

HOW THEY ARE UPDATED. Once per training round (one epoch), from every residue of the round's
training proteins: the maps it reads (an interior residue its left and right map at share 1/2 each,
the library's mixture; a terminal residue its one map at share 1) and its basin populations in the
native-restrained replica and in the free ensemble. One damped Newton step on the maximum a
posteriori objective, solved jointly over all shared offsets and deviations by conjugate gradient,
with z = (q, d) and c = P z:

    [P^T (sum_r J_r^T Cov_r J_r) P + Lambda] dz = P^T T0 sum_r J_r^T (p_free - p_native)_r - Lambda z
    z += ETA dz,    Lambda = T0^2/sigma^2 for d, T0^2/SIGMA_POOL^2 for q

with J_r the residue's shares on its maps' offsets and Cov_r = diag(p) - p p^T at
p = (p_free + p_native)/2, the first-order response of its basin populations. A basin the free
simulation over-populates is raised, one it under-populates is lowered. ETA < 1 also absorbs the
coil/sheet mixture, which dilutes a coil map's share. There is no unfolded-state (DSE) term: the
SARW replica keeps the rama term, so it would compare the map with itself. A tenth of the proteins
is held out of the update (ConDiv.py); each round reports how much of their map-level gap the step
would remove, the check that the offsets transfer between proteins.
"""

import os
import shutil

import numpy as np
import tables as tb

MODE = 'all'                                       # 'all', 'data' or 'pooled' (see the top)
SIGMA_ALL = 0.3                                    # prior width, nats (findings 1.29)
SIGMA_DATA = 0.2
SIGMA_POOL = 0.5                                   # the shared offsets' width, cross-validated (plan.md Phase 29)
MIN_NATIVE = 20.                                   # native sites a 'data' depth or a shared offset needs
ETA = 0.5
BASINS = ('alpha_R', 'alpha_L', 'beta', 'pPII', "beta'", "pPII'", 'other')
N_BASIN = len(BASINS)
N_OFF = 6                                          # the six basins; 'other' is beta' + pPII'
# Logistic edge scale: at 3 deg the 10-90% transition is 13 deg, two to three grid cells, the
# sharpest the engine's spline represents cleanly (round 2).
EDGE = np.deg2rad(3.)
CG_TOL, CG_MAX_ITER = 1e-10, 2000


# ---------------------------------------------------------------------------
# Basins
# ---------------------------------------------------------------------------

def _wrap(x):
    """Angle difference into [-pi, pi)."""
    return (np.asarray(x) + np.pi) % (2. * np.pi) - np.pi


def _sigmoid(x):
    return 1. / (1. + np.exp(-x))


def _arc(x, a, b):
    """Smooth periodic indicator of the arc from a to b (radians, anticlockwise): a sigmoid of the
    circular distance from the arc's midpoint, continuous across phi = +-180, and an arc mirrors
    exactly onto the arc (-b, -a)."""
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


def residue_populations(rama_coord, weights=None):
    """Per-residue basin populations, (n_res, N_BASIN), of an ensemble of (n_frame, n_res, 2)
    (phi, psi) in radians, with per-frame weights summing to one or None for a plain mean."""
    w = basin_weights(rama_coord[..., 0], rama_coord[..., 1])
    if weights is None:
        return w.mean(axis=0)
    return np.tensordot(np.asarray(weights, dtype=float), w, axes=1)


def native_dihedrals(pos):
    """(phi, psi) of every residue of an N, CA, C backbone (3 n_res, 3), radians; the first phi
    and the last psi are NaN."""
    n, ca, c = pos[0::3], pos[1::3], pos[2::3]

    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
        v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
        w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
        return np.arctan2((np.cross(b1, v) * w).sum(-1), (v * w).sum(-1))
    phi = np.full(len(n), np.nan)
    psi = np.full(len(n), np.nan)
    phi[1:] = dih(c[:-1], n[1:], ca[1:], c[1:])
    psi[:-1] = dih(n[:-1], ca[:-1], c[:-1], n[1:])
    return phi, psi


# ---------------------------------------------------------------------------
# The parameters: one deviation per (coil map, basin), one shared offset per (central residue, basin)
# ---------------------------------------------------------------------------

def _decode(attr):
    return [x.decode() if isinstance(x, bytes) else x for x in attr]


def sigma(state):
    """The deviations' prior width."""
    return SIGMA_ALL if state['mode'] == 'all' else SIGMA_DATA


def residue_reads(state, seq):
    """For each residue, its two map slots (map index, share), index -1 where there is none: an
    interior residue reads its left and right map at 1/2 each, a terminal residue one map at 1."""
    seq = list(seq)
    n = len(seq)
    m = np.full((n, 2), -1, dtype=int)
    s = np.zeros((n, 2))
    for i, c in enumerate(seq):
        for slot, (d, j) in enumerate(((state['left'], i - 1), (state['right'], i + 1))):
            if 0 <= j < n:
                m[i, slot] = state['index'][(state['ridx'][c], d, state['ridx']['PRO' if seq[j] == 'CPR' else seq[j]])]
        s[i] = (0.5, 0.5) if 0 < i < n - 1 else (m[i] >= 0).astype(float)
    return m, s


def init_state(library, natives, mode=MODE):
    """Offsets zero (NDRD) on every coil map the engine reads; natives, the training proteins'
    (sequence, N CA C positions), give each (map, basin)'s native sites, which set the trainable
    deviations in modes 'data' and 'pooled', and each (central residue, basin)'s, which set the
    trainable shared offsets in mode 'pooled'."""
    with tb.open_file(library) as t:
        restype = _decode(t.root.coil._v_attrs.restype)
        dirs = _decode(t.root.coil._v_attrs.dir)
        n_central = t.root.coil.dimer_pot.shape[0]
    keys = [(c, d, n) for c in range(n_central) for d in range(len(dirs))
            for n, r in enumerate(restype) if r not in ('CPR', 'ALL')]
    right = dirs.index('right')
    state = dict(source=os.path.abspath(library), mode=mode, keys=np.array(keys, dtype=int),
                 restype=restype, dirs=dirs, ridx={r: i for i, r in enumerate(restype)},
                 index={k: i for i, k in enumerate(keys)}, left=dirs.index('left'), right=right,
                 n_central=n_central, offset=np.zeros((len(keys), N_OFF)),
                 deviation=np.zeros((len(keys), N_OFF)), pooled=np.zeros((n_central, N_OFF)),
                 pool_member=np.array([not (d == right and restype[n] == 'PRO') for c, d, n in keys]),
                 round=0, history=[])
    native = np.zeros((len(keys), N_OFF))
    native_central = np.zeros((n_central, N_OFF))
    for seq, pos in natives:
        m, _ = residue_reads(state, seq)
        phi, psi = native_dihedrals(pos)
        ok = np.isfinite(phi) & np.isfinite(psi)
        w = basin_weights(phi[ok], psi[ok])[:, :N_OFF]
        for slot in range(2):
            mm = m[ok, slot]
            np.add.at(native, mm[mm >= 0], w[mm >= 0])
        np.add.at(native_central, [state['ridx'][c] for c in np.asarray(list(seq))[ok]], w)
    state['native_sites'] = native
    state['native_sites_central'] = native_central
    state['trainable'] = np.ones_like(native, dtype=bool) if mode == 'all' else native >= MIN_NATIVE
    state['pool_trainable'] = (native_central >= MIN_NATIVE) if mode == 'pooled' else np.zeros_like(native_central, dtype=bool)
    return state


def n_param(state):
    return int(state['trainable'].sum() + state['pool_trainable'].sum())


def _spread(state, q):
    """P's shared part: each map's central residue's shared offsets, zero on pre-proline maps."""
    return state['pool_member'][:, None] * q[state['keys'][:, 0]]


def _collect(state, g):
    """The transpose of _spread: a per-map vector summed onto the maps' central residues."""
    out = np.zeros((state['n_central'], N_OFF))
    np.add.at(out, state['keys'][:, 0], state['pool_member'][:, None] * g)
    return out


def _normalise(m):
    """The library's convention: sum(exp(-E)) == 1 per map."""
    lo = m.min(axis=(-2, -1), keepdims=True)
    return m + np.log(np.exp(-(m - lo)).sum(axis=(-2, -1), keepdims=True)) - lo


def write_library(state, source, out):
    """Copy `source` and replace every coil map the engine reads by its base plus its offsets,
    renormalised. Nothing else changes."""
    shutil.copy(source, out)
    with tb.open_file(state['source']) as t:
        base = t.root.coil.dimer_pot[:].astype(float)
    W = grid_basins(base.shape[-1])[..., :N_OFF]
    with tb.open_file(out, 'a') as t:
        pot = t.root.coil.dimer_pot[:]
        for k, (c, d, n) in enumerate(state['keys']):
            pot[c, d, n] = _normalise(base[c, d, n] + np.einsum('b,ijb->ij', state['offset'][k], W))
        t.root.coil.dimer_pot[:] = pot
    return out


# ---------------------------------------------------------------------------
# Statistics and the update
# ---------------------------------------------------------------------------

def empty_accumulator(state):
    return dict(n=0., m=np.zeros((0, 2), dtype=int), s=np.zeros((0, 2)), p=np.zeros((0, N_OFF)),
                gap=np.zeros((0, N_OFF)))


def accumulate(acc, state, seq, free, native):
    """Add one protein's residues: map slots, shares, mean and free-minus-native populations."""
    m, s = residue_reads(state, seq)
    free, native = np.asarray(free)[:, :N_OFF], np.asarray(native)[:, :N_OFF]
    merge(acc, dict(n=float((m >= 0).sum()), m=m, s=s, p=0.5 * (free + native), gap=free - native))


def merge(acc, add):
    acc['n'] += add['n']
    for f in ('m', 's', 'p', 'gap'):
        acc[f] = np.concatenate([acc[f], add[f]])


def _scatter(acc, v, n_map):
    """sum_r J_r^T v_r: each residue's per-basin vector onto its maps, by share."""
    out = np.zeros((n_map, N_OFF))
    for slot in range(2):
        ok = acc['m'][:, slot] >= 0
        np.add.at(out, acc['m'][ok, slot], acc['s'][ok, slot, None] * v[ok])
    return out


def _gather(acc, c):
    """J_r c for every residue: its maps' offsets, by share."""
    out = np.zeros((len(acc['m']), N_OFF))
    for slot in range(2):
        ok = acc['m'][:, slot] >= 0
        out[ok] += acc['s'][ok, slot, None] * c[acc['m'][ok, slot]]
    return out


def _cov(p, v):
    return p * v - p * (p * v).sum(1, keepdims=True)


def data_gradient(state, acc):
    """T0 sum_r J_r^T (p_free - p_native)_r, the data part of the objective's gradient."""
    return state['T0'] * _scatter(acc, acc['gap'], len(state['keys']))


def _lambdas(state):
    """The priors' curvatures on the deviations and on the shared offsets."""
    return state['T0'] ** 2 / sigma(state) ** 2, state['T0'] ** 2 / SIGMA_POOL ** 2


def newton_step(state, acc):
    """The MAP Newton step (dq, dd) on the trainable shared offsets and deviations (the rest stay zero)."""
    mq, md = state['pool_trainable'].astype(float), state['trainable'].astype(float)
    lam_d, lam_q = _lambdas(state)
    n_q = mq.size

    def split(x):
        return x[:n_q].reshape(mq.shape), x[n_q:].reshape(md.shape)

    def A(x):
        q, d = split(x)
        h = _scatter(acc, _cov(acc['p'], _gather(acc, _spread(state, q) + d)), len(state['keys']))
        return np.concatenate([(mq * (_collect(state, h) + lam_q * q)).ravel(), (md * (h + lam_d * d)).ravel()])

    g = data_gradient(state, acc)
    b = np.concatenate([(mq * (_collect(state, g) - lam_q * state['pooled'])).ravel(),
                        (md * (g - lam_d * state['deviation'])).ravel()])
    x = np.zeros_like(b)
    r = b.copy()
    d = r.copy()
    rr, b2 = (r * r).sum(), (b * b).sum()
    for _ in range(CG_MAX_ITER):
        if rr <= CG_TOL * b2:
            break
        Ad = A(d)
        a = rr / (d * Ad).sum()
        x += a * d
        r -= a * Ad
        rr_new = (r * r).sum()
        d = r + (rr_new / rr) * d
        rr = rr_new
    return split(x)


def map_level_gap(state, acc, change=None):
    """sum over maps of (the map's summed gap)^2 / its reads, after the predicted first-order
    response to an offset change, if given: the part of the gap an offset can act on."""
    gap = acc['gap'] if change is None else acc['gap'] - _cov(acc['p'], _gather(acc, change)) / state['T0']
    g = _scatter(acc, gap, len(state['keys']))
    n = _scatter(acc, np.ones_like(gap), len(state['keys']))[:, 0]
    return float((g ** 2 / np.maximum(n, 1e-9)[:, None]).sum())


def step_gradient(state, acc):
    """One step's share of the objective's gradient, on the deviations (n_map, N_OFF) and on the
    shared offsets (n_central, N_OFF): the data part, free minus native, plus this step's share of
    each prior's pull, zero where not trainable."""
    lam_d, lam_q = _lambdas(state)
    g = data_gradient(state, acc)
    n = state['steps_per_round']
    return (state['trainable'] * (g - lam_d * state['deviation'] / n),
            state['pool_trainable'] * (_collect(state, g) - lam_q * state['pooled'] / n))


def step_record(state, step_train, step_held):
    """What a step keeps for the record (rama_step.npz)."""
    grad, grad_pooled = step_gradient(state, step_train)
    return dict(grad=grad, grad_pooled=grad_pooled, train_n=step_train['n'], held_n=step_held['n'],
                train_gap=_scatter(step_train, step_train['gap'], len(state['keys'])),
                held_gap=_scatter(step_held, step_held['gap'], len(state['keys'])))


def movable(state):
    """The (map, basin) depths that can move: a trainable deviation, or a shared offset reaching the map."""
    return state['trainable'] | (state['pool_member'][:, None] & state['pool_trainable'][state['keys'][:, 0]])


def update(state, acc_train, acc_held):
    """One round's damped Newton step. Returns the step on the maps' offsets and the round's summary,
    which reports the held-out proteins' map-level gap before and after the step's predicted response."""
    dq, dd = newton_step(state, acc_train)
    dq, dd = ETA * dq, ETA * dd
    step = _spread(state, dq) + dd
    held_before = map_level_gap(state, acc_held)
    held_after = map_level_gap(state, acc_held, step)
    train_before = map_level_gap(state, acc_train)
    state['pooled'] = state['pooled'] + dq
    state['deviation'] = state['deviation'] + dd
    state['offset'] = _spread(state, state['pooled']) + state['deviation']
    state['round'] += 1
    return step, dict(reads_train=acc_train['n'], reads_held=acc_held['n'],
                      train_gap=train_before / max(acc_train['n'], 1.),
                      held_gap=held_before / max(acc_held['n'], 1.),
                      held_left=held_after / max(held_before, 1e-30),
                      step_rms=float(np.sqrt((step[movable(state)] ** 2).mean())),
                      step_max=float(np.abs(step).max()))


def offset_summary(state):
    """Offsets by basin (rms over the depths that can move), counts, glycine's GLY|X alpha_R and
    alpha_L means for comparison with ff3.0's pooled pair, and the shared offsets by basin (rms over
    the trainable ones) with glycine's own."""
    c, t, q, tq = state['offset'], movable(state), state['pooled'], state['pool_trainable']
    gly = [k for k, (cc, d, n) in enumerate(state['keys']) if state['restype'][cc] == 'GLY'
           and state['restype'][n] != 'GLY' and not (d == state['right'] and state['restype'][n] == 'PRO')]
    g = state['ridx']['GLY']
    return dict(n_trainable=n_param(state), n_above_01=int((np.abs(c) > 0.1).sum()),
                rms=[float(np.sqrt((c[t[:, b], b] ** 2).mean())) if t[:, b].any() else 0. for b in range(N_OFF)],
                gly_aR=float(c[gly, 0].mean()), gly_aL=float(c[gly, 1].mean()), max=float(np.abs(c).max()),
                pool_rms=[float(np.sqrt((q[tq[:, b], b] ** 2).mean())) if tq[:, b].any() else 0. for b in range(N_OFF)],
                pool_gly_aR=float(q[g, 0]), pool_gly_aL=float(q[g, 1]))


def round_line(h):
    """A round's line in rama_rounds.txt."""
    return ('round %2i  reads %6.0f / %5.0f held out  map-level gap %.5f / held out %.5f  held-out '
            'gap left after the step %.4f  step rms %.4f max %.4f  offsets rms %s  |c| > 0.1: %i of %i  '
            'GLY|X aR %+.4f aL %+.4f  shared rms %s  GLY shared aR %+.4f aL %+.4f\n'
            % (h['round'], h['reads_train'], h['reads_held'], h['train_gap'], h['held_gap'],
               h['held_left'], h['step_rms'], h['step_max'],
               ' '.join('%s %.4f' % (b, v) for b, v in zip(BASINS, h['rms'])), h['n_above_01'],
               h['n_trainable'], h['gly_aR'], h['gly_aL'],
               ' '.join('%s %.4f' % (b, v) for b, v in zip(BASINS, h['pool_rms'])),
               h['pool_gly_aR'], h['pool_gly_aL']))


def describe(state):
    """One line for the trainer's parameter print and check_step.py."""
    o = offset_summary(state)
    widths = ('sigma %.2f, shared %.2f' % (sigma(state), SIGMA_POOL) if state['mode'] == 'pooled'
              else 'sigma %.2f' % sigma(state))
    return ('round %i, mode %s, %s: %i trainable depths on %i coil maps, %i with |c| > 0.1, '
            'max |c| %.4f; rms by basin %s; GLY|X mean alpha_R %+.4f alpha_L %+.4f; shared GLY alpha_R %+.4f '
            'alpha_L %+.4f'
            % (state['round'], state['mode'], widths, o['n_trainable'], len(state['keys']),
               o['n_above_01'], o['max'], ' '.join('%s %.4f' % (b, v) for b, v in zip(BASINS, o['rms'])),
               o['gly_aR'], o['gly_aL'], o['pool_gly_aR'], o['pool_gly_aL']))
