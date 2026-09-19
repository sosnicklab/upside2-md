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

PARAMETERISATION. Two maps, `S` symmetric and `A` antisymmetric under (phi,psi) -> (-phi,-psi):

    X|GLY  (every neighbour, both directions) = S + A
    GLY|GLY                                   = S

That is exactly the constraint the physics demands. A glycine flanked by L-amino acids sits in a
chiral environment and may be biased; a glycine flanked by glycines has no chirality source, so
its map must be mirror symmetric. Writing the pair as S and A makes GLY|GLY symmetric by
construction, with no projection step during training, while S + A stays unconstrained.

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


def map_gradient(coords, n_grid, cardinal=None):
    """dE/d(map[i,j]) accumulated over samples.

    coords: (n_sample, 2) array of (phi, psi) in radians for ONE residue over frames.
    Returns (n_grid, n_grid). Summed, not averaged; the caller decides the normalisation.
    """
    if cardinal is None:
        cardinal = cardinal_function(n_grid)
    g = grid_coords(coords, n_grid)
    wx = sample_weights(g[:, 0], n_grid, cardinal)
    wy = sample_weights(g[:, 1], n_grid, cardinal)
    return wx.T @ wy


def energy_from_map(coords, rama_pot, cardinal=None):
    """Reimplementation of RamaMapPot, for validating against the engine."""
    n_grid = rama_pot.shape[-1]
    return float((map_gradient(coords, n_grid, cardinal) * rama_pot).sum())


# ---------------------------------------------------------------------------
# The mixture half: d(per-residue map)/d(S, A), by autograd
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
    """Forward map from (S, A) to the per-residue maps the engine sees, for glycine residues.

    Holds everything about one protein that does not change during training: which residues are
    glycine, their neighbours, the coil and sheet mixing weights, and the sheet maps.
    """

    def __init__(self, seq, rama_library, sheet_mixing_energy, central='GLY'):
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

        def s_nb(name):                      # sheet group has no CPR
            return sidx['PRO'] if name == 'CPR' else sidx[name]

        self.residues = []
        for i, c in enumerate(seq):
            if c != central:
                continue
            # read_rama_maps_and_weights mixes left and right for an interior residue but uses a
            # SINGLE direction at each terminus, and a terminal glycine still gets a glycine map,
            # so its energy still depends on (S, A). Skipping termini left 3.0% of the glycine
            # gradient missing across the 456-protein training set, and a biased 3%: chain ends
            # are more flexible than the interior. Both cases are handled here by carrying a list
            # of (direction, neighbour) pairs that has one entry at a terminus and two inside.
            if i == 0:
                pairs = [('right', seq[1])]
            elif i == len(seq) - 1:
                pairs = [('left', seq[-2])]
            else:
                pairs = [('left', seq[i - 1]), ('right', seq[i + 1])]

            cw_i = [float(cw[gi_c, self.coil_dir.index(d), cidx['PRO'] if nb == 'CPR' else cidx[nb]])
                    for d, nb in pairs]
            sw_i = [float(sw[gi_s, sheet_dir.index(d), s_nb(nb)]) for d, nb in pairs]
            sp_i = np.stack([spot[gi_s, sheet_dir.index(d), s_nb(nb)] for d, nb in pairs])
            # the per-residue weight is the mean over directions used, which for one direction
            # is just that weight
            self.residues.append(dict(
                index=i,
                is_central=[nb == central for _, nb in pairs],
                coil_w=cw_i,
                sheet_w=sw_i,
                sheet_pot=torch.as_tensor(sp_i, dtype=torch.float64),
                sheet_mix=float(sheet_mix[sidx[central]]),
                coil_mix_w=float(np.mean(cw_i)),
                sheet_mix_w=float(np.mean(sw_i)),
            ))

    def residue_map(self, S, A, r):
        """The map upside_config writes for glycine residue r, as a differentiable tensor.

        Every stage of `read_weighted_maps` plus the Boltzmann-weighted shift that
        `write_rama_map_pot` applies last. That shift is a constant per map, so it changes no
        force and is invisible in any basin-to-basin difference, but it is map-dependent and it
        enters the total energy, so it carries gradient. Leaving it out made the reconstructed
        map wrong by a constant 7.23 and the dS gradient wrong by 37%.
        """
        v_xg = S + A
        v_gg = S
        branches = torch.stack([v_gg if is_c else v_xg for is_c in r['is_central']])
        coil = _normalize(_mixture(r['coil_w'], branches))
        sheet = _normalize(_mixture(r['sheet_w'], r['sheet_pot']))
        m = _mixture([r['coil_mix_w'], r['sheet_mix_w'] * np.exp(-r['sheet_mix'])],
                     torch.stack([coil, sheet]))
        return m - (m * torch.exp(-m)).sum()

    def backprop(self, S, A, residue_grads):
        """Given dE/d(per-residue map) for each glycine, return (dE/dS, dE/dA) as numpy."""
        St = torch.as_tensor(S, dtype=torch.float64).clone().requires_grad_(True)
        At = torch.as_tensor(A, dtype=torch.float64).clone().requires_grad_(True)
        total = St.new_zeros(())
        for r in self.residues:
            g = residue_grads.get(r['index'])
            if g is None:
                continue
            total = total + (self.residue_map(St, At, r)
                             * torch.as_tensor(g, dtype=torch.float64)).sum()
        if not total.requires_grad:
            z = np.zeros_like(np.asarray(S, dtype=float))
            return z, z.copy()
        total.backward()
        return St.grad.numpy(), At.grad.numpy()


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
# Writing the trained maps back into a library
# ---------------------------------------------------------------------------

def write_gly_library(source, out, S, A, central='GLY'):
    """Copy `source` and set its central-`central` coil row to X|GLY = S+A, GLY|GLY = S.

    Nothing else is touched: the sheet group, every other central residue, both weight arrays and
    the all-NaN CPR neighbour column (never read, since a cis-proline *neighbour* maps to PRO)
    come through unchanged. No renormalisation is applied, because `read_rama_maps_and_weights`
    normalises the mixture itself and an extra shift here would change the inner left/right
    mixing weights in a way the gradient does not model.
    """
    import shutil
    shutil.copy(source, out)
    xg = np.asarray(S) + np.asarray(A)
    gg = np.asarray(S)
    with tb.open_file(out, 'a') as t:
        coil = t.root.coil
        restype = _decode(coil._v_attrs.restype)
        g = restype.index(central)
        pot = coil.dimer_pot[:]
        for d in range(pot.shape[1]):
            for n in range(pot.shape[2]):
                if np.isnan(pot[g, d, n]).all():
                    continue
                pot[g, d, n] = gg if n == g else xg
        coil.dimer_pot[:] = pot
    return out


def read_gly_maps(library, central='GLY'):
    """Recover (S, A) from any library whose central-`central` coil row holds two maps.

    `S` is the symmetric part of the `GLY|GLY` map and `A` is the **antisymmetric part of an
    `X|GLY` map**, rather than the difference of the two. Differencing looks equivalent and is
    not: `build_rama_from_awh.py` normalises each map separately, so the two carry different
    additive constants and `X|GLY - GLY|GLY` returns `A` plus an offset. That offset is symmetric,
    so projecting removes it, and for a library written by `write_gly_library` (which does not
    renormalise) the two routes agree exactly anyway.

    This is also how a worker gets its current parameters: they are already in the library file it
    was handed, so no 5,184-value array has to be threaded through the command line.
    """
    with tb.open_file(library) as t:
        coil = t.root.coil
        restype = _decode(coil._v_attrs.restype)
        g = restype.index(central)
        S = coil.dimer_pot[g, 0, g][:]
        xg = None
        for n in range(coil.dimer_pot.shape[2]):
            if n == g:
                continue
            m = coil.dimer_pot[g, 0, n][:]
            if np.isfinite(m).all():
                xg = m
                break
    if xg is None:
        raise ValueError(f'{library} has no finite X|{central} coil map')
    return (project_symmetric(np.asarray(S, dtype=float)),
            project_antisymmetric(np.asarray(xg, dtype=float)))


def symmetric_start(source, central='GLY'):
    """The starting S for training: the library's coil row symmetrised, with A = 0.

    Training then begins assuming no glycine handedness and has to find it in the protein data,
    which is what makes the comparison against the AWH measurement a real test rather than a
    circular one.
    """
    with tb.open_file(source) as t:
        coil = t.root.coil
        restype = _decode(coil._v_attrs.restype)
        g = restype.index(central)
        m = coil.dimer_pot[g, 0, restype.index('ALL') if 'ALL' in restype else -1][:]
    if not np.isfinite(m).all():
        raise ValueError('starting map has non-finite entries')
    S = project_symmetric(m)
    S = S + np.log(np.exp(-S).sum())      # cosmetic: match the library's sum(exp(-E)) == 1.
    return S, np.zeros_like(S)            # A constant shift cancels in every downstream mixture.


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
