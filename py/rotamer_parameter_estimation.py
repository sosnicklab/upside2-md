import os
import pickle

import numpy as np
import scipy.optimize as opt

# ---------------------------------------------------------------------------
# Dimensionality constants
# ---------------------------------------------------------------------------
n_fix        = 3
n_rotpos     = 86
n_restype    = 20
n_cov        = 8    # coverage interaction types (backbone-group × donor/acceptor)
n_hyd        = 12   # hydrophobe interaction types (n_fix × 4)

n_knot_angular = 15
n_knot_sc      = 12
n_knot_hb      = 10
hb_dr          = 0.625
sc_dr          = 0.7


# ---------------------------------------------------------------------------
# The latent vector ConDiv trains, block by block. Angular knots are logits
# (the tables hold their sigmoids), radial knots the interior knots of a
# clamped spline (the pair table's symmetrised over the two residue types),
# directions unnormalised.
# ---------------------------------------------------------------------------

latent_shapes = [
    (n_restype, n_restype, n_knot_angular),   # pair dp1 (dp2 is its transpose)
    (n_restype, n_restype, n_knot_sc - 3),    # pair uni
    (n_restype, n_restype, n_knot_sc - 3),    # pair direc
    (n_cov, n_restype, n_knot_angular),       # coverage dp1
    (n_cov, n_restype, n_knot_angular),       # coverage dp2
    (n_cov, n_restype, n_knot_hb - 3),        # coverage uni
    (n_cov, n_restype, n_knot_hb - 3),        # coverage direc
    (n_hyd, n_restype, n_knot_angular),       # hydrophobe dp1
    (n_hyd, n_restype, n_knot_angular),       # hydrophobe dp2
    (n_hyd, n_restype, n_knot_hb - 3),        # hydrophobe uni
    (n_hyd, n_restype, n_knot_hb - 3),        # hydrophobe direc
    (n_fix, 3), (n_fix, 3),                   # hydrophobe placement point, direction
    (n_rotpos, 3), (n_rotpos, 3),             # rotamer centre point, direction
]
n_param = sum(int(np.prod(s)) for s in latent_shapes)


def _blocks(lparam):
    lparam = np.asarray(lparam, dtype=np.float64)
    if lparam.shape != (n_param,):
        raise ValueError(f'a side-chain latent vector has {n_param} coordinates, not shape {lparam.shape} '
                         f'(a ConDiv checkpoint written before 2026-10-22 has 31420: training/README.md)')
    edges = np.cumsum([0] + [int(np.prod(s)) for s in latent_shapes])
    return [lparam[a:b].reshape(s) for a, b, s in zip(edges[:-1], edges[1:], latent_shapes)]


def _sigmoid(x):
    return 1. / (1. + np.exp(-x))


def _symm(x):
    """0.5 (x + x^T) over the two residue-type axes, so E(A,B) = E(B,A); its own adjoint."""
    return 0.5 * (x + x.transpose(1, 0, 2))


def _clamp_spline(middle):
    """The knots of a clamped spline from its interior ones: the left end clamped to
    middle[1], the right end to zero."""
    return np.concatenate([middle[..., 1:2], middle, -0.5 * middle[..., -1:], middle[..., -1:]], axis=-1)


def _clamp_spline_vjp(d):
    d_middle = d[..., 1:-2].copy()
    d_middle[..., 1] += d[..., 0]
    d_middle[..., -1] += d[..., -1] - 0.5 * d[..., -2]
    return d_middle


def _norm(u):
    return np.sqrt((u ** 2).sum(axis=-1, keepdims=True))


def _normalize_vjp(u, d):
    """The adjoint of u -> u/|u|: the part of d perpendicular to u, over |u|."""
    y = u / _norm(u)
    return (d - y * (y * d).sum(axis=-1, keepdims=True)) / _norm(u)


# ---------------------------------------------------------------------------
# Unpacking: latent vector -> side-chain tables, and its gradient
# ---------------------------------------------------------------------------

def unpack_params(lparam):
    """The side-chain tables of a latent vector: pair_interaction (20, 20, 54),
    coverage_interaction (8, 20, 50), hydrophobe_interaction (12, 20, 50),
    hydrophobe_placement (3, 7) and rotamer_center_fixed (86, 6)."""
    (rot_dp, rot_uni, rot_dir, cov_dp1, cov_dp2, cov_uni, cov_dir, hyd_dp1, hyd_dp2, hyd_uni, hyd_dir,
     hydpl_com, hydpl_dir, rotpos_com, rotpos_dir) = _blocks(lparam)

    rot_dp = _sigmoid(rot_dp)
    rot = np.concatenate([rot_dp, rot_dp.transpose(1, 0, 2),
                          _clamp_spline(_symm(rot_uni)), _clamp_spline(_symm(rot_dir))], axis=-1)
    cov = np.concatenate([_sigmoid(cov_dp1), _sigmoid(cov_dp2),
                          _clamp_spline(cov_uni), _clamp_spline(cov_dir)], axis=-1)
    hyd = np.concatenate([_sigmoid(hyd_dp1), _sigmoid(hyd_dp2),
                          _clamp_spline(hyd_uni), _clamp_spline(hyd_dir)], axis=-1)
    hydpl = np.concatenate([hydpl_com, hydpl_dir / _norm(hydpl_dir), np.zeros((n_fix, 1))], axis=-1)
    rotpos = np.concatenate([rotpos_com, rotpos_dir / _norm(rotpos_dir)], axis=-1)
    return rot, cov, hyd, hydpl, rotpos


def unpack_params_vjp(lparam, d_rot, d_cov, d_hyd, d_hydpl, d_rotpos):
    """The gradient in lparam of sum(d_x * x) over the tables x of unpack_params(lparam)."""
    (rot_dp, _, _, cov_dp1, cov_dp2, _, _, hyd_dp1, hyd_dp2, _, _,
     _, hydpl_dir, _, rotpos_dir) = _blocks(lparam)
    na, nsc, nhb = n_knot_angular, n_knot_sc, n_knot_hb

    def sigmoid_vjp(x, d):
        s = _sigmoid(x)
        return d * s * (1. - s)

    def table_vjp(dp1, dp2, d):
        return [sigmoid_vjp(dp1, d[..., :na]), sigmoid_vjp(dp2, d[..., na:2*na]),
                _clamp_spline_vjp(d[..., 2*na:2*na + nhb]), _clamp_spline_vjp(d[..., 2*na + nhb:])]

    grads = ([sigmoid_vjp(rot_dp, d_rot[..., :na] + d_rot[..., na:2*na].transpose(1, 0, 2)),
              _symm(_clamp_spline_vjp(d_rot[..., 2*na:2*na + nsc])),
              _symm(_clamp_spline_vjp(d_rot[..., 2*na + nsc:]))] +
             table_vjp(cov_dp1, cov_dp2, d_cov) +
             table_vjp(hyd_dp1, hyd_dp2, d_hyd) +
             [d_hydpl[:, :3], _normalize_vjp(hydpl_dir, d_hydpl[:, 3:6]),
              d_rotpos[:, :3], _normalize_vjp(rotpos_dir, d_rotpos[:, 3:6])])
    return np.concatenate([g.ravel() for g in grads])


# ---------------------------------------------------------------------------
# pack_param: find latent vector that reproduces given parameter arrays
# ---------------------------------------------------------------------------

def _logit(p, eps=1e-6):
    """Inverse sigmoid, clamped to avoid inf."""
    p = np.clip(p, eps, 1. - eps)
    return np.log(p / (1. - p))


def _clamp_middle(full_spline):
    """The interior knots of a _clamp_spline output, full_spline[1:n+1]."""
    return full_spline[..., 1:-2]


def _init_x0(rotp, covp, hydp, hydplp, rotposp):
    """Construct a warm-start latent vector from structured parameter arrays.

    This is an approximate inversion of unpack_params.  It only sets the
    starting point for pack_param's solve; the residual check there is what
    establishes that the result is right.  The Theano original started from a
    flat 0.5 vector instead.
    """
    na, nsc, nhb = n_knot_angular, n_knot_sc, n_knot_hb
    parts = [
        _logit(rotp[:, :, :na]),                       # pair dp1
        _clamp_middle(rotp[:, :, 2*na:2*na + nsc]),    # pair uni
        _clamp_middle(rotp[:, :, 2*na + nsc:]),        # pair direc
        _logit(covp[:, :, :na]),
        _logit(covp[:, :, na:2*na]),
        _clamp_middle(covp[:, :, 2*na:2*na + nhb]),
        _clamp_middle(covp[:, :, 2*na + nhb:]),
        _logit(hydp[:, :, :na]),
        _logit(hydp[:, :, na:2*na]),
        _clamp_middle(hydp[:, :, 2*na:2*na + nhb]),
        _clamp_middle(hydp[:, :, 2*na + nhb:]),
        np.asarray(hydplp)[:, :3],
        np.asarray(hydplp)[:, 3:6],    # target directions are unit vectors; set equal
        np.asarray(rotposp)[:, :3],
        np.asarray(rotposp)[:, 3:6],
    ]
    x0 = np.concatenate([p.ravel() for p in parts])
    assert len(x0) == n_param, f'x0 length {len(x0)} != n_param {n_param}'
    return x0


def pack_param(rotp, covp, hydp, hydplp, rotposp):
    """Inverse of unpack_params: find latent vector x s.t. unpack_params(x) ≈ inputs."""
    targets = [np.asarray(t, dtype=np.float64) for t in [rotp, covp, hydp, hydplp, rotposp]]

    def obj_and_grad(x):
        resid = [a - b for a, b in zip(unpack_params(x), targets)]
        return sum((r ** 2).sum() for r in resid), unpack_params_vjp(x, *[2. * r for r in resid])

    x0     = _init_x0(rotp, covp, hydp, hydplp, rotposp)
    result = opt.minimize(
        obj_and_grad, x0, jac=True, method='L-BFGS-B',
        options={'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-10},
    )
    # Gate on the residual, as the Theano original does: the unpacking must
    # reproduce the input arrays, and a stationary point that does not is a
    # failure, not a convergence.
    if not (result.fun < 1.6e-4):
        raise ValueError(
            f'Failed to converge: residual {result.fun:.4g} >= 1.6e-4 '
            f'({result.message})'
        )
    print(f'pack_param residual = {result.fun:.6g}')
    return result.x


# ---------------------------------------------------------------------------
# Quadratic spline energy (used for regularization in ConDiv d_obj)
# ---------------------------------------------------------------------------

def _grid_values(sp):
    """A quadratic spline's values at its grid points 1..n-2."""
    return (1. / 6.) * sp[..., :-2] + (2. / 3.) * sp[..., 1:-1] + (1. / 6.) * sp[..., 2:]


def _grid_values_vjp(d):
    sp = np.zeros(d.shape[:-1] + (d.shape[-1] + 2,))
    sp[..., :-2] += (1. / 6.) * d
    sp[..., 1:-1] += (2. / 3.) * d
    sp[..., 2:] += (1. / 6.) * d
    return sp


def _split_knots(params, n_knots):
    edges = np.cumsum([0] + list(n_knots))
    return [params[..., a:b] for a, b in zip(edges[:-1], edges[1:])]


def quadspline_energy(params, n_knots):
    """Evaluate the quadratic spline potential on its grid.

    Args:
        params:  array of shape (..., sum(n_knots))
        n_knots: tuple (n1, n2, n3, n4) giving knot counts for (dp1, dp2, uni, direc)

    Returns:
        array of shape (..., n3-2, n1-2, n2-2): uni(r) + dp1(a) dp2(b) direc(r)
    """
    dp1, dp2, uni, direc = (_grid_values(x) for x in _split_knots(params, n_knots))
    return (
        uni  [..., :, None, None] +
        dp1  [..., None, :, None] *
        dp2  [..., None, None, :] *
        direc[..., :, None, None]
    )


def quadspline_energy_vjp(params, n_knots, d_energy):
    """The gradient in params of sum(d_energy * quadspline_energy(params, n_knots))."""
    dp1, dp2, _, direc = (_grid_values(x) for x in _split_knots(params, n_knots))
    return np.concatenate([_grid_values_vjp(d) for d in (
        np.einsum('...rab,...b,...r->...a', d_energy, dp2, direc),
        np.einsum('...rab,...a,...r->...b', d_energy, dp1, direc),
        d_energy.sum(axis=(-2, -1)),
        np.einsum('...rab,...a,...b->...r', d_energy, dp1, dp2),
    )], axis=-1)


# ---------------------------------------------------------------------------
# Adam optimizer
# ---------------------------------------------------------------------------

def _read_comp(x, i):
    """A setting shared by every component, or component i of a per-component one."""
    return x if np.isscalar(x) else x[i]


class AdamSolver:
    """Adam optimizer for a fixed number of independent parameter components.

    Each component may be a scalar or an array; arithmetic is element-wise.
    alpha is roughly the maximum step size.
    """

    def __init__(self, n_comp, alpha=1e-2, beta1=0.8, beta2=0.96, epsilon=1e-6):
        self.n_comp   = n_comp
        self.alpha    = alpha
        self.beta1    = beta1
        self.beta2    = beta2
        self.epsilon  = epsilon
        self.step_num = 0
        self.grad1    = [0. for _ in range(n_comp)]
        self.grad2    = [0. for _ in range(n_comp)]

    def update_step(self, grad):
        """Consume gradient list, return list of parameter updates."""
        r = _read_comp
        self.step_num += 1
        t = self.step_num
        updates = [None] * self.n_comp

        for i, g in enumerate(grad):
            b1 = r(self.beta1, i)
            self.grad1[i] = b1 * self.grad1[i] + (1. - b1) * g
            m_hat = self.grad1[i] / (1. - b1 ** t)

            b2 = r(self.beta2, i)
            self.grad2[i] = b2 * self.grad2[i] + (1. - b2) * g ** 2
            v_hat = self.grad2[i] / (1. - b2 ** t)

            updates[i] = -r(self.alpha, i) * m_hat / (np.sqrt(v_hat) + r(self.epsilon, i))

        return updates

    def log_state(self, direc):
        with open(os.path.join(direc, 'solver_state.pkl'), 'wb') as f:
            pickle.dump(
                dict(step_num=self.step_num, grad1=self.grad1,
                     grad2=self.grad2, solver=str(self)),
                f, -1,
            )

    def __repr__(self):
        return (f'AdamSolver({self.n_comp}, alpha={self.alpha!r}, '
                f'beta1={self.beta1!r}, beta2={self.beta2!r}, '
                f'epsilon={self.epsilon!r})')

    def __str__(self):
        return (f'AdamSolver({self.n_comp}, alpha={self.alpha}, '
                f'beta1={self.beta1}, beta2={self.beta2}, '
                f'epsilon={self.epsilon})')
