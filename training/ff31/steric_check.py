"""Steric check of a per-map depth run (ff3.1): does the trained Ramachandran library put more
weight where a residue's side chain and backbone clash? (repo plan.md Phase 14 step 5)

    python3 steric_check.py <run_dir> [<epoch_xx_minibatch_yy>]

Upside has no excluded volume between a residue's CB and its neighbours' backbone: the map keeps
them apart. The reference is the hard-sphere Ramachandran map. The fragment
CA-C(=O)-N(H)-CA(CB)-C(=O)-N(H)-CA is built at trans omega from ideal geometry (Engh & Huber
1991); a (phi, psi) node of the library grid clashes when an atom pair whose distance depends on
phi or psi is closer than Ramachandran's outer limit (Ramachandran & Sasisekharan 1968). Glycine
has no CB; before a proline the next residue's CD takes the amide H's place. A central proline is
read with the general fragment (its ring is not built). Two clash classes:
  side chain  a pair with the CB, or with the next proline's CD: the risk the watch guards;
  backbone    the other pairs. The rigid geometry puts a quarter of NDRD's mass here (alpha_R's
              upper phi edge by C(i-1)..C(i), the bridge by N(i)..H(i+1)), which proteins occupy
              by widening N-CA-C, so only its change is read.

The step's checkpoint holds the force field as it stands: the library its next step reads and its
sheet energies. Three readouts against the start:
  1. every coil map: its probability on clashing nodes, the step's library against NDRD;
  2. every training residue's engine map (left/right mixture, then the sheet group at the sheet
     energies, as upside_config builds rama_map_pot): now, at the start (NDRD and the initial
     sheet energies), and NDRD at the step's sheet energies, which separates the depths' part;
  3. every finished step's free and native-restrained phi > 0 populations of non-glycine residues
     (divergence.pkl), what the simulations did.
"""

import glob
import os
import pickle as cp
import sys
from types import SimpleNamespace

import numpy as np
import tables as tb

RUN = os.path.abspath(sys.argv[1])
OUT = os.path.join(RUN, 'run_output')
sys.path.insert(0, OUT)
import ConDiv  # noqa: E402  the run's own copy
import __main__  # noqa: E402
for _n in dir(ConDiv):
    if _n[0].isupper():
        setattr(__main__, _n, getattr(ConDiv, _n))
import upside_config as uc  # noqa: E402

deg = np.deg2rad

# Outer-limit contact distances, A, by element pair.
LIMIT = {('C', 'C'): 3.0, ('C', 'O'): 2.7, ('C', 'N'): 2.8, ('N', 'N'): 2.6, ('N', 'O'): 2.6,
         ('O', 'O'): 2.6, ('C', 'H'): 2.2, ('H', 'N'): 2.2, ('H', 'O'): 2.2, ('H', 'H'): 1.9}
SIDE = {'CB1', 'CD2'}
CLASSES = ('side chain', 'backbone')
RATIO_FLAG = 2.     # a map or residue whose clash probability more than doubles (either class)
MEAN_FLAG = 1.1     # the training residues' mean side-chain clash probability up by more than 10%


# ---------------------------------------------------------------------------
# The hard-sphere reference
# ---------------------------------------------------------------------------

def _place(a, b, c, bond, angle, torsion):
    """The atom bonded to c at `bond`, with angle b-c-d and dihedral a-b-c-d (radians)."""
    bc = c - b
    bc /= np.linalg.norm(bc, axis=-1, keepdims=True)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    m = np.cross(n, bc)
    torsion = np.asarray(torsion)[..., None]
    return c + bond * (-np.cos(angle) * bc + np.sin(angle) * (np.cos(torsion) * m + np.sin(torsion) * n))


def fragment(phi, psi, has_cb, pre_pro):
    """Atom positions (name -> (n, 3)) of the fragment at each (phi, psi); a name's first letter
    is its element."""
    g = np.ones(phi.shape + (1,))
    ca0 = g * np.zeros(3)
    c0 = g * np.array([1.525, 0., 0.])
    n1 = c0 + g * 1.329 * np.array([-np.cos(deg(116.2)), np.sin(deg(116.2)), 0.])
    ca1 = _place(ca0, c0, n1, 1.458, deg(121.7), np.pi)
    at = dict(CA0=ca0, C0=c0, N1=n1, CA1=ca1)
    at['O0'] = _place(ca1, n1, c0, 1.231, deg(123.0), 0.)
    at['H1'] = _place(at['O0'], c0, n1, 1.0, deg(119.15), np.pi)
    at['C1'] = _place(c0, n1, ca1, 1.525, deg(111.2), phi)
    at['N2'] = _place(n1, ca1, at['C1'], 1.329, deg(116.2), psi)
    at['CA2'] = _place(ca1, at['C1'], at['N2'], 1.458, deg(121.7), np.pi)
    at['O1'] = _place(at['CA2'], at['N2'], at['C1'], 1.231, deg(123.0), 0.)
    if pre_pro:
        at['CD2'] = _place(at['O1'], at['C1'], at['N2'], 1.473, deg(125.0), np.pi)
    else:
        at['H2'] = _place(at['O1'], at['C1'], at['N2'], 1.0, deg(119.15), np.pi)
    if has_cb:
        # The ideal CB of an L residue from N, CA, C.
        b, c = ca1 - n1, at['C1'] - ca1
        at['CB1'] = -0.58273431 * np.cross(b, c) + 0.56802827 * b - 0.54067466 * c + ca1
    return at


def clash_maps(n_grid, has_cb, pre_pro):
    """(2, n_grid, n_grid) [class, phi, psi]: the library grid's side-chain and backbone clashes."""
    t = -np.pi + 2. * np.pi * np.arange(n_grid) / n_grid
    phi, psi = np.meshgrid(t, t, indexing='ij')
    at = fragment(phi.ravel(), psi.ravel(), has_cb, pre_pro)
    names = sorted(at)
    clash = np.zeros((2, phi.size), dtype=bool)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = np.linalg.norm(at[a] - at[b], axis=-1)
            if d.std() > 1e-3:                      # the pair moves with phi or psi
                clash[0 if SIDE & {a, b} else 1] |= d < LIMIT[tuple(sorted((a[0], b[0])))]
    return clash.reshape((2,) + phi.shape)


# ---------------------------------------------------------------------------
# Libraries and engine maps
# ---------------------------------------------------------------------------

def _decode(attr):
    return [x.decode() if isinstance(x, bytes) else x for x in attr]


def read_library(path):
    """The coil and sheet groups in memory, in the form read_rama_maps_and_weights reads."""
    with tb.open_file(path) as t:
        groups = {}
        for name in ('coil', 'sheet'):
            g = getattr(t.root, name)
            groups[name] = SimpleNamespace(
                _v_attrs=SimpleNamespace(restype=_decode(g._v_attrs.restype), dir=_decode(g._v_attrs.dir)),
                dimer_pot=g.dimer_pot[:], dimer_weight=g.dimer_weight[:])
    return groups


def engine_maps(seq, lib, sheet_energy):
    """Each residue's normalised probability map, as write_rama_map_pot builds it."""
    coil, coil_w = uc.read_rama_maps_and_weights(seq, lib['coil'])
    sheet, sheet_w = uc.read_rama_maps_and_weights(seq, lib['sheet'], allow_CPR=False)
    srid = {r: i for i, r in enumerate(lib['sheet']._v_attrs.restype)}
    mix = np.array([sheet_energy[srid['PRO' if r == 'CPR' else r]] for r in seq])
    pot = uc.mixture_potential([coil_w, sheet_w * np.exp(-mix)], [coil, sheet])
    p = np.exp(-(pot - pot.min(axis=(-2, -1), keepdims=True)))
    return p / p.sum(axis=(-2, -1), keepdims=True)


def clash_probability(p, refs):
    """(..., 2) probability of each clash class, for maps p (..., n, n) and references refs."""
    return np.einsum('...ij,...kij->...k', p, refs)


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

done = sorted(d for d in glob.glob(os.path.join(OUT, 'epoch_*_minibatch_*'))
              if os.path.exists(os.path.join(d, 'checkpoint.pkl')))
step = os.path.join(OUT, sys.argv[2]) if len(sys.argv) > 2 else done[-1]
state = cp.load(open(os.path.join(step, 'checkpoint.pkl'), 'rb'))
init = cp.load(open(os.path.join(OUT, 'initial_checkpoint.pkl'), 'rb'))
rama = state['rama']
restype = rama['restype']
seqs = {nm: list(uc.read_fasta(open(t.fasta))) for mb in state['minibatches'] for nm, t in mb}

now_lib, ndrd_lib = read_library(rama['library']), read_library(rama['source'])
n_grid = now_lib['coil'].dimer_pot.shape[-1]
REF = {(cb, pp): clash_maps(n_grid, cb, pp) for cb in (True, False) for pp in (True, False)}
flags = []

print(f'{RUN}: {os.path.basename(step)}, depth round {rama["round"]}, library {rama["library"]}\n')
print('hard-sphere reference, share of the plane clashing (side chain / backbone): '
      + ', '.join(f'{"CB" if cb else "no CB"}{" pre-Pro" if pp else ""} {r[0].mean():.3f} / {r[1].mean():.3f}'
                  for (cb, pp), r in REF.items()))

# 1. Every coil map.
pm = []
for c, d, n in rama['keys']:
    ref = REF[restype[c] != 'GLY', d == rama['right'] and restype[n] == 'PRO']
    for lib in (now_lib, ndrd_lib):
        p = np.exp(-lib['coil'].dimer_pot[c, d, n].astype(float))
        pm.append(clash_probability(p / p.sum(), ref))
pm = np.array(pm).reshape(len(rama['keys']), 2, 2)          # [map, now / NDRD, class]
ratio = pm[:, 0] / np.maximum(pm[:, 1], 1e-300)
print(f'\n1. coil maps ({len(pm)}), clash probability, step against NDRD:')
for k, cls in enumerate(CLASSES):
    r = ratio[:, k]
    print(f'   {cls:10s} NDRD median {np.median(pm[:, 1, k]):.4f} (max {pm[:, 1, k].max():.4f}); step / NDRD '
          f'median {np.median(r):.4f}, max {r.max():.4f}, maps above {RATIO_FLAG:g}x: {(r > RATIO_FLAG).sum()}')
    for j in np.argsort(-r)[:3]:
        c, d, n = rama['keys'][j]
        print(f'      {restype[c]} {rama["dirs"][d]} {restype[n]}: {pm[j, 1, k]:.4f} -> {pm[j, 0, k]:.4f} '
              f'(x{r[j]:.3f}); offsets {" ".join(f"{x:+.3f}" for x in rama["offset"][j])}')
    if (r > RATIO_FLAG).any():
        flags.append(f'{(r > RATIO_FLAG).sum()} coil maps above {RATIO_FLAG:g}x ({cls})')

# 2. Every training residue's engine map.
sheet_now, sheet_init = np.asarray(state['param'].sheet), np.asarray(init['param'].sheet)
names, pr = [], []
for code, seq in sorted(seqs.items()):
    refs = np.array([REF[r != 'GLY', i + 1 < len(seq) and seq[i + 1] in ('PRO', 'CPR')]
                     for i, r in enumerate(seq)])
    pr.append(np.stack([clash_probability(engine_maps(seq, lib, sh), refs)
                        for lib, sh in ((now_lib, sheet_now), (ndrd_lib, sheet_init), (ndrd_lib, sheet_now))], 1))
    names += [(code, i, '-'.join(seq[max(i - 1, 0):i + 2])) for i in range(len(seq))]
pr = np.concatenate(pr)                                       # [residue, now / start / NDRD now, class]
m = pr.mean(axis=0)
rr = pr[:, 0] / np.maximum(pr[:, 1], 1e-300)
print(f'\n2. engine maps of {len(seqs)} training proteins ({len(pr)} residues), mean clash probability '
      f'now / start / NDRD at the step\'s sheet energies:')
for k, cls in enumerate(CLASSES):
    print(f'   {cls:10s} {m[0, k]:.5f} / {m[1, k]:.5f} / {m[2, k]:.5f}; now / start {m[0, k] / m[1, k]:.4f} '
          f'(depths {m[0, k] / m[2, k]:.4f}, sheet {m[2, k] / m[1, k]:.4f}); per residue median '
          f'{np.median(rr[:, k]):.4f}, max {rr[:, k].max():.4f}, above {RATIO_FLAG:g}x: {(rr[:, k] > RATIO_FLAG).sum()}')
    for j in np.argsort(-rr[:, k])[:3]:
        print(f'      {names[j][0]} {names[j][1]:3d} {names[j][2]}: {pr[j, 1, k]:.4f} -> {pr[j, 0, k]:.4f} '
              f'(x{rr[j, k]:.3f})')
    if (rr[:, k] > RATIO_FLAG).any():
        flags.append(f'{(rr[:, k] > RATIO_FLAG).sum()} residues above {RATIO_FLAG:g}x ({cls})')
if m[0, 0] / m[1, 0] > MEAN_FLAG:
    flags.append(f'training residues\' mean side-chain clash probability x{m[0, 0] / m[1, 0]:.3f}')

# 3. What the simulations did.
print('\n3. non-glycine residues, free / restrained basin population (divergence.pkl); steps 0-18 '
      'ran on NDRD\n   step  alpha_L          beta\' + pPII\'')
series = []
for d in done:
    fr, na = [], []
    for f in glob.glob(os.path.join(d, '*.divergence.pkl')):
        s, dv = seqs[os.path.basename(f).split('.')[0]], cp.load(open(f, 'rb'))
        keep = np.array([r != 'GLY' for r in s])
        fr.append(np.asarray(dv['rama_free'])[keep])
        na.append(np.asarray(dv['rama_native'])[keep])
    fr, na = np.concatenate(fr), np.concatenate(na)
    series.append((fr[:, 1].mean(), na[:, 1].mean(), fr[:, 6].mean(), na[:, 6].mean()))
    print(f'   {len(series) - 1:4d}  {series[-1][0]:.4f} / {series[-1][1]:.4f}  {series[-1][2]:.4f} / '
          f'{series[-1][3]:.4f}')
series = np.array(series)
if len(series) > 19 and series[19:, 2].max() > series[:19, 2].max():
    flags.append(f"free beta' + pPII' {series[19:, 2].max():.4f} above steps 0-18's maximum "
                 f"{series[:19, 2].max():.4f}")

print('\nSTERIC FLAG: ' + '; '.join(flags) if flags else '\nno steric flag')
