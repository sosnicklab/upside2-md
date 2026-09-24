"""Exact gradient of the Upside energy with respect to a trainable glycine Ramachandran map.

WHY THIS IS POSSIBLE WITHOUT FINITE DIFFERENCES. `rama_map_pot` evaluates a periodic
interpolating bicubic spline, and `solve_periodic_2d_spline` (src/spline.cpp:262) is a tensor
product of 1D periodic solves. The map-to-energy operator is therefore linear, separable and
translation invariant on the grid:

    E = sum_r sum_ij  map_r[i,j] * b(x_r - i) * b(y_r - j)

with a single 1D cardinal function `b`, the periodic cubic spline interpolating a unit impulse.
So dE/d(map_r[i,j]) is just a spline-smoothed 2D histogram of the glycine (phi,psi) samples.
Finite differencing instead would cost 2 extra passes per parameter, i.e. 10,369x a divergence
for a 72x72 map, which is why `sheet` is trained as a single scalar.

PARAMETERISATION. Every finite map in the library's central-glycine coil row is its own
parameter: one map per (direction, neighbour), 42 in the ff2.1 layout (2 directions x 21 neighbour
types; the CPR column is NaN and never read, because a cis-proline neighbour maps onto PRO). The
two GLY|GLY maps are held exactly mirror-symmetric under (phi,psi) -> (-phi,-psi), each on its
own, because a glycine flanked by glycines has no chirality source; every X|GLY map is left free,
because an L-amino-acid neighbour makes the environment chiral. The constraint is kept by
projecting both the starting maps and every update, so it holds to machine precision throughout.

THE CHAIN RULE IS DONE BY AUTOGRAD, DELIBERATELY. Between the library maps and the per-residue
map the engine actually sees, `read_weighted_maps` applies a normalisation, a weighted log-sum-exp
over the left and right neighbour maps, another normalisation, and a second log-sum-exp against
the sheet map. Hand-differentiating that is four chances to get a softmax factor or a sign wrong,
and a wrong sign trains in the wrong direction for days without ever failing. The forward pass is
reimplemented in torch and autograd supplies the adjoint.
"""

import numpy as np
import tables as tb
import torch
from scipy.interpolate import CubicSpline


# ---------------------------------------------------------------------------
# The spline half: dE/d(per-residue map)
# ---------------------------------------------------------------------------

def cardinal_function(n_grid):
    """Periodic cubic spline interpolating a unit impulse at node 0, as a callable of grid units."""
    d = np.zeros(n_grid)
    d[0] = 1.
    return CubicSpline(np.arange(n_grid + 1), np.r_[d, d[0]], bc_type='periodic')


def grid_coords(rama_coord, n_grid):
    """Radians to grid units, matching RamaMapPot's scale and shift exactly."""
    scale = n_grid * (0.5 / np.pi - 1e-7)
    return (np.asarray(rama_coord) + np.pi) * scale


def sample_weights(u, n_grid, cardinal):
    """Weight vectors w[s,i] = b(u_s - i) for samples u (grid units). Shape (n_sample, n_grid)."""
    off = (np.asarray(u)[:, None] - np.arange(n_grid)[None, :]) % n_grid
    return cardinal(off)


def map_gradient(coords, n_grid, cardinal=None, weights=None):
    """dE/d(map[i,j]) accumulated over samples.

    coords: (n_sample, 2) array of (phi, psi) in radians for ONE residue over frames.
    weights: optional per-sample weights; an ensemble average is `weights` summing to one.
    Returns (n_grid, n_grid). Summed with the given weights (1 each by default).
    """
    if cardinal is None:
        cardinal = cardinal_function(n_grid)
    g = grid_coords(coords, n_grid)
    wx = sample_weights(g[:, 0], n_grid, cardinal)
    wy = sample_weights(g[:, 1], n_grid, cardinal)
    if weights is not None:
        wy = wy * np.asarray(weights, dtype=float)[:, None]
    return wx.T @ wy


def energy_from_map(coords, rama_pot, cardinal=None):
    """Reimplementation of RamaMapPot, for validating against the engine."""
    n_grid = rama_pot.shape[-1]
    return float((map_gradient(coords, n_grid, cardinal) * rama_pot).sum())


# ---------------------------------------------------------------------------
# The mixture half: d(per-residue map)/d(row maps), by autograd
# ---------------------------------------------------------------------------

def _normalize(p):
    """The library's convention, sum(exp(-E)) == 1, as upside_config applies it."""
    return p + torch.logsumexp(-p.reshape(-1), 0)


def _mixture(weights, potentials):
    """upside_config.mixture_potential, in torch. weights are scalars, potentials (n,g,g)."""
    w = torch.as_tensor(weights, dtype=potentials.dtype)
    w = w / w.sum()
    shifted = potentials - torch.log(w)[:, None, None]
    return -torch.logsumexp(-shifted, 0)


class GlycineMapChain:
    """Forward map from the glycine row's maps to the per-residue maps the engine sees.

    Holds everything about one protein that does not change during training: which residues are
    glycine, which row map each of their neighbour branches reads, the coil and sheet mixing
    weights, and the sheet maps. `keys` is the (direction, neighbour) index of every parameter
    map, as returned by `row_keys`.
    """

    def __init__(self, seq, rama_library, sheet_mixing_energy, keys, central='GLY'):
        seq = list(seq)
        self.seq = seq
        self.central = central
        with tb.open_file(rama_library) as t:
            coil, sheet = t.root.coil, t.root.sheet
            self.coil_restype = _decode(coil._v_attrs.restype)
            self.coil_dir = _decode(coil._v_attrs.dir)
            sheet_restype = _decode(sheet._v_attrs.restype)
            sheet_dir = _decode(sheet._v_attrs.dir)
            cw, sw = coil.dimer_weight[:], sheet.dimer_weight[:]
            spot = sheet.dimer_pot[:]
        sheet_mix = np.loadtxt(sheet_mixing_energy)
        self.n_grid = spot.shape[-1]

        gi_c = self.coil_restype.index(central)
        gi_s = sheet_restype.index(central)
        sidx = {r: i for i, r in enumerate(sheet_restype)}
        cidx = {r: i for i, r in enumerate(self.coil_restype)}
        kidx = {(int(d), int(n)): k for k, (d, n) in enumerate(keys)}

        # read_rama_maps_and_weights reads a cis-proline NEIGHBOUR as PRO in both groups
        def c_nb(name):
            return cidx['PRO'] if name == 'CPR' else cidx[name]

        def s_nb(name):
            return sidx['PRO'] if name == 'CPR' else sidx[name]

        self.residues = []
        for i, c in enumerate(seq):
            if c != central:
                continue
            # A terminal glycine reads one direction and an interior one mixes two, as in
            # read_rama_maps_and_weights. Skipping termini once left 3.0% of the glycine gradient
            # missing across the training set, and a biased 3%: chain ends are more flexible.
            if i == 0:
                pairs = [('right', seq[1])]
            elif i == len(seq) - 1:
                pairs = [('left', seq[-2])]
            else:
                pairs = [('left', seq[i - 1]), ('right', seq[i + 1])]

            dirs = [self.coil_dir.index(d) for d, _ in pairs]
            cw_i = [float(cw[gi_c, d, c_nb(nb)]) for d, (_, nb) in zip(dirs, pairs)]
            sw_i = [float(sw[gi_s, sheet_dir.index(d), s_nb(nb)]) for d, nb in pairs]
            sp_i = np.stack([spot[gi_s, sheet_dir.index(d), s_nb(nb)] for d, nb in pairs])
            self.residues.append(dict(
                index=i,
                map_k=[kidx[(d, c_nb(nb))] for d, (_, nb) in zip(dirs, pairs)],
                coil_w=cw_i,
                sheet_w=sw_i,
                sheet_pot=torch.as_tensor(sp_i, dtype=torch.float64),
                sheet_mix=float(sheet_mix[sidx[central]]),
                # the per-residue weight is the mean over the directions used
                coil_mix_w=float(np.mean(cw_i)),
                sheet_mix_w=float(np.mean(sw_i)),
            ))

    def residue_map(self, G, r):
        """The map upside_config writes for glycine residue r, as a differentiable tensor.

        Every stage of `read_weighted_maps` plus the Boltzmann-weighted shift that
        `write_rama_map_pot` applies last. That shift is a constant per map, so it changes no
        force, but it enters the total energy, so it carries gradient. Leaving it out made the
        reconstructed map wrong by a constant 7.23 and the gradient wrong by 37%.
        """
        branches = torch.stack([G[k] for k in r['map_k']])
        coil = _normalize(_mixture(r['coil_w'], branches))
        sheet = _normalize(_mixture(r['sheet_w'], r['sheet_pot']))
        m = _mixture([r['coil_mix_w'], r['sheet_mix_w'] * np.exp(-r['sheet_mix'])],
                     torch.stack([coil, sheet]))
        return m - (m * torch.exp(-m)).sum()

    def backprop(self, G, residue_grads):
        """Given dE/d(per-residue map) for each glycine, return dE/dG, shape of G, as numpy."""
        Gt = torch.as_tensor(np.asarray(G), dtype=torch.float64).clone().requires_grad_(True)
        total = Gt.new_zeros(())
        for r in self.residues:
            g = residue_grads.get(r['index'])
            if g is None:
                continue
            total = total + (self.residue_map(Gt, r)
                             * torch.as_tensor(g, dtype=torch.float64)).sum()
        if not total.requires_grad:
            return np.zeros_like(np.asarray(G, dtype=float))
        total.backward()
        return Gt.grad.numpy()


def _decode(attr):
    return [x.decode() if isinstance(x, bytes) else x for x in attr]


# ---------------------------------------------------------------------------
# Symmetry
# ---------------------------------------------------------------------------

def mirror(m):
    """(phi,psi) -> (-phi,-psi): a reversal WITH a roll, because index i maps to (-i) % n."""
    return np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)


def project_symmetric(m):
    return 0.5 * (m + mirror(m))


def project_antisymmetric(m):
    return 0.5 * (m - mirror(m))


# ---------------------------------------------------------------------------
# The glycine row as a parameter: keys, reading, writing, the constraint
# ---------------------------------------------------------------------------

def row_keys(library, central='GLY'):
    """(direction, neighbour) indices of every finite map in the central row, and which are GLY|GLY."""
    with tb.open_file(library) as t:
        coil = t.root.coil
        restype = _decode(coil._v_attrs.restype)
        g = restype.index(central)
        row = coil.dimer_pot[g]
    keys = np.array([(d, n) for d in range(row.shape[0]) for n in range(row.shape[1])
                     if np.isfinite(row[d, n]).all()], dtype=int)
    return keys, keys[:, 1] == g


def read_row(library, keys, central='GLY'):
    """The row's maps in `keys` order, shape (n_map, n_grid, n_grid)."""
    with tb.open_file(library) as t:
        coil = t.root.coil
        g = _decode(coil._v_attrs.restype).index(central)
        row = coil.dimer_pot[g]
    return np.stack([row[d, n] for d, n in keys]).astype(float)


def write_row(source, out, keys, G, central='GLY'):
    """Copy `source` and replace its central row's maps with G.

    Nothing else in the file changes: the sheet group, every other central residue, both weight
    arrays and the all-NaN CPR neighbour column come through as they were. No renormalisation is
    applied, because `read_rama_maps_and_weights` normalises the mixture itself and a shift here
    would move the inner left/right mixing in a way the gradient does not model.
    """
    import shutil
    shutil.copy(source, out)
    with tb.open_file(out, 'a') as t:
        coil = t.root.coil
        g = _decode(coil._v_attrs.restype).index(central)
        pot = coil.dimer_pot[:]
        for k, (d, n) in enumerate(keys):
            pot[g, d, n] = G[k]
        coil.dimer_pot[:] = pot
    return out


def constrain(G, is_gg):
    """Project the GLY|GLY maps onto their mirror-symmetric part; X|GLY maps pass through."""
    G = np.array(G, dtype=float)
    G[is_gg] = project_symmetric(G[is_gg])
    return G


def start_row(library, central='GLY'):
    """The starting row: the library's own maps with each GLY|GLY map symmetrised.

    ff2.1's GLY|GLY maps are not mirror-symmetric (ff2.1 dG(aR->aL) -0.72 on an achiral pair), so
    they are projected once; every X|GLY map starts exactly at ff2.1's value.
    """
    keys, is_gg = row_keys(library, central)
    return keys, is_gg, constrain(read_row(library, keys, central), is_gg)


def handedness(m):
    """dG(aR->aL) in nats on the basins every glycine number in this project uses."""
    n = m.shape[-1]
    t = np.arange(-180., 180., 360. / n)
    phi, psi = np.meshgrid(t, t, indexing='ij')

    def basin(a, b, c, d):
        k = (phi >= a) & (phi <= b) & (psi >= c) & (psi <= d)
        return -np.log(np.exp(-m[..., k]).sum(-1))
    return basin(40., 100., -10., 60.) - basin(-100., -40., -60., 10.)


# ---------------------------------------------------------------------------
# Smoothing
# ---------------------------------------------------------------------------

def fourier_lowpass(m, order):
    """Keep 2D Fourier modes with |kx| <= order and |ky| <= order on the torus.

    The exact gradient is already spread over a 4x4 node neighbourhood by the spline basis, but
    5,184 free values against ~30,000 glycine samples per minibatch is still noisy. Band-limiting
    the *update* keeps the learned correction smooth while the starting map keeps its sharp
    forbidden-region structure.
    """
    if order is None:
        return m
    n = m.shape[-1]
    f = np.fft.fft2(m)
    keep = np.zeros((n, n), dtype=bool)
    k = np.fft.fftfreq(n, 1. / n).astype(int)
    kx = np.abs(k)[:, None] <= order
    ky = np.abs(k)[None, :] <= order
    keep[kx & ky] = True
    return np.real(np.fft.ifft2(f * keep))
