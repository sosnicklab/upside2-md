"""Ramachandran basins, and per-residue basin populations of an ensemble.

A diagnostic, not a parameter: every training step records, per residue, the basin populations of
the native-restrained replica and of the free ensemble (`rama_native`, `rama_free` in each
`<code>.divergence.pkl`). Split by each residue's own native basin, they show where the free
simulation leaves the native conformation, for instance helical glycines visiting alpha_L
(findings 1.11-1.16).

THE BASINS partition the torus: alpha_R (phi < 0, -100 < psi < 50), beta (phi < -100, psi outside
that band), pPII (-100 < phi < 0, psi outside it), alpha_L (the mirror of alpha_R), and phi > 0
outside alpha_L split into beta' and pPII', the mirrors of beta and pPII; `other` is their union.
Edges are logistic with a 3 deg scale (13 deg from 10% to 90%), continuous across phi = +-180, and
every basin mirrors exactly onto its partner under (phi, psi) -> (-phi, -psi).
"""

import numpy as np

BASINS = ('alpha_R', 'alpha_L', 'beta', 'pPII', "beta'", "pPII'", 'other')
N_BASIN = len(BASINS)
EDGE = np.deg2rad(3.)


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
