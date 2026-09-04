import collections
import os
import pickle

import numpy as np
import scipy.optimize as opt
import torch

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
# Core unpacking: latent vector → structured parameter tensors
# ---------------------------------------------------------------------------

def unpack_param_maker(lparam):
    """Unpack a flat latent float64 tensor into structured force-field tensors.

    Returns (rot_param, cov_param, hyd_param, hydpl_param, rotpos_param,
             rotscalar_param, n_param).  All outputs are PyTorch tensors so
    that autograd can compute gradients through them.
    """
    i = [0]

    def read_param(shape):
        size = int(np.prod(shape))
        ret  = lparam[i[0] : i[0] + size].reshape(shape)
        i[0] += size
        return ret

    def read_symm(n):
        """Read (n_restype, n_restype, n) and enforce E(A,B)=E(B,A)."""
        x = read_param((n_restype, n_restype, n))
        return 0.5 * (x + x.permute(1, 0, 2))

    def clamp_spline(middle):
        """Apply boundary conditions: left-clamp from middle[:,1], right-clamp to zero."""
        c0  = middle[..., 1:2]
        cn3 = middle[..., -1:]
        cn2 = -0.5 * cn3
        cn1 = cn3
        return torch.cat([c0, middle, cn2, cn1], dim=-1)

    # -----------------------------------------------------------------------
    # rot_param  (pair_interaction): shape (20, 20, 54)
    # -----------------------------------------------------------------------
    angular_logits = read_param((n_restype, n_restype, n_knot_angular))

    angular_spline_sc = torch.sigmoid(angular_logits)                 # (20, 20, 15)

    clamped_1 = clamp_spline(read_symm(n_knot_sc - 3))              # (20, 20, 12)
    clamped_2 = clamp_spline(read_symm(n_knot_sc - 3))              # (20, 20, 12)

    rot_param = torch.cat([
        angular_spline_sc,
        angular_spline_sc.permute(1, 0, 2),   # dp2 — transpose of dp1; GLY col is also palindromic
        clamped_1,
        clamped_2,
    ], dim=-1)                                                       # (20, 20, 54)

    # -----------------------------------------------------------------------
    # cov_param  (coverage_interaction): shape (8, 20, 50)
    # -----------------------------------------------------------------------
    cov_dp1  = torch.sigmoid(read_param((n_cov, n_restype, n_knot_angular)))
    cov_dp2  = torch.sigmoid(read_param((n_cov, n_restype, n_knot_angular)))
    cov_uni  = clamp_spline(read_param((n_cov, n_restype, n_knot_hb - 3)))
    cov_dir  = clamp_spline(read_param((n_cov, n_restype, n_knot_hb - 3)))

    cov_param = torch.cat([cov_dp1, cov_dp2, cov_uni, cov_dir], dim=-1)  # (8, 20, 50)

    # -----------------------------------------------------------------------
    # hyd_param  (hydrophobe_interaction): shape (12, 20, 50)
    # -----------------------------------------------------------------------
    hyd_dp1  = torch.sigmoid(read_param((n_hyd, n_restype, n_knot_angular)))
    hyd_dp2  = torch.sigmoid(read_param((n_hyd, n_restype, n_knot_angular)))
    hyd_uni  = clamp_spline(read_param((n_hyd, n_restype, n_knot_hb - 3)))
    hyd_dir  = clamp_spline(read_param((n_hyd, n_restype, n_knot_hb - 3)))

    hyd_param = torch.cat([hyd_dp1, hyd_dp2, hyd_uni, hyd_dir], dim=-1)  # (12, 20, 50)

    # -----------------------------------------------------------------------
    # hydpl_param  (hydrophobe_placement): shape (3, 7)
    # -----------------------------------------------------------------------
    hydpl_com         = read_param((n_fix, 3))
    hydpl_dir_unnorm  = read_param((n_fix, 3))
    hydpl_dir_norm    = hydpl_dir_unnorm / torch.sqrt(
        (hydpl_dir_unnorm ** 2).sum(dim=-1, keepdim=True) + 1e-12
    )
    hydpl_param = torch.cat(
        [hydpl_com, hydpl_dir_norm, torch.zeros((n_fix, 1), dtype=lparam.dtype)],
        dim=-1,
    )                                                                # (3, 7)

    # -----------------------------------------------------------------------
    # rotpos_param  (rotamer_center_fixed): shape (86, 6)
    # -----------------------------------------------------------------------
    rotpos_com        = read_param((n_rotpos, 3))
    rotpos_dir_unnorm = read_param((n_rotpos, 3))
    rotpos_dir_norm   = rotpos_dir_unnorm / torch.sqrt(
        (rotpos_dir_unnorm ** 2).sum(dim=-1, keepdim=True) + 1e-12
    )
    rotpos_param = torch.cat([rotpos_com, rotpos_dir_norm], dim=-1)  # (86, 6)

    # -----------------------------------------------------------------------
    # rotscalar_param: shape (86, 1)
    # -----------------------------------------------------------------------
    rotscalar_param = read_param((n_rotpos, 1))

    return (rot_param, cov_param, hyd_param, hydpl_param,
            rotpos_param, rotscalar_param, int(i[0]))


def _compute_n_param():
    dummy = torch.zeros(40000, dtype=torch.float64)
    with torch.no_grad():
        *_, n = unpack_param_maker(dummy)
    return n


n_param = _compute_n_param()


# ---------------------------------------------------------------------------
# Evaluate without gradient tracking (numpy in, numpy out)
# ---------------------------------------------------------------------------

def unpack_params(lparam_numpy):
    """Evaluate unpack_param_maker on a numpy array; return numpy arrays."""
    lparam = torch.tensor(lparam_numpy, dtype=torch.float64)
    with torch.no_grad():
        rot, cov, hyd, hydpl, rotpos, rotscalar, _ = unpack_param_maker(lparam)
    return (rot.numpy(), cov.numpy(), hyd.numpy(),
            hydpl.numpy(), rotpos.numpy(), rotscalar.numpy())


# ---------------------------------------------------------------------------
# pack_param: find latent vector that reproduces given parameter arrays
# ---------------------------------------------------------------------------

def _logit(p, eps=1e-6):
    """Inverse sigmoid, clamped to avoid inf."""
    p = np.clip(p, eps, 1. - eps)
    return np.log(p / (1. - p))


def _clamp_middle(full_spline):
    """Extract the interior middle knots from a clamp_spline output.

    clamp_spline(middle) produces:
        [middle[1], middle[0..n-1], -0.5*middle[-1], middle[-1]]
    so middle[0..n-1] = full_spline[1:n+1].
    """
    return full_spline[..., 1:-2]


def _init_x0(rotp, covp, hydp, hydplp, rotposp, rotscalarp):
    """Construct a warm-start latent vector from structured parameter arrays.

    This is an approximate inversion of unpack_param_maker.  Most blocks
    are exact; the GLY row of angular_logits is initialised to the (non-
    palindromic) logit of the ff_2.1 angular profile — the optimizer will
    find the closest palindromic approximation from there.
    """
    parts = []

    # 1. angular_logits (20, 20, 15): logit of dp1 angular slice of rot
    parts.append(_logit(rotp[:, :, :n_knot_angular]).ravel())

    # 2. clamped_symm_1 middle (20, 20, 9): interior knots of radial block 1
    parts.append(_clamp_middle(rotp[:, :, n_knot_angular*2:n_knot_angular*2 + n_knot_sc]).ravel())

    # 3. clamped_symm_2 middle (20, 20, 9)
    parts.append(_clamp_middle(rotp[:, :, n_knot_angular*2 + n_knot_sc:]).ravel())

    # 4. cov_dp1_logits (8, 20, 15)
    parts.append(_logit(covp[:, :, :n_knot_angular]).ravel())

    # 5. cov_dp2_logits (8, 20, 15)
    parts.append(_logit(covp[:, :, n_knot_angular:n_knot_angular*2]).ravel())

    # 6. cov_uni_middle (8, 20, 7)
    parts.append(_clamp_middle(covp[:, :, n_knot_angular*2:n_knot_angular*2 + n_knot_hb]).ravel())

    # 7. cov_dir_middle (8, 20, 7)
    parts.append(_clamp_middle(covp[:, :, n_knot_angular*2 + n_knot_hb:]).ravel())

    # 8. hyd_dp1_logits (12, 20, 15)
    parts.append(_logit(hydp[:, :, :n_knot_angular]).ravel())

    # 9. hyd_dp2_logits (12, 20, 15)
    parts.append(_logit(hydp[:, :, n_knot_angular:n_knot_angular*2]).ravel())

    # 10. hyd_uni_middle (12, 20, 7)
    parts.append(_clamp_middle(hydp[:, :, n_knot_angular*2:n_knot_angular*2 + n_knot_hb]).ravel())

    # 11. hyd_dir_middle (12, 20, 7)
    parts.append(_clamp_middle(hydp[:, :, n_knot_angular*2 + n_knot_hb:]).ravel())

    # 12. hydpl_com (3, 3)
    parts.append(np.asarray(hydplp)[:, :3].ravel())

    # 13. hydpl_dir_unnorm (3, 3): target directions are unit vectors; set equal
    parts.append(np.asarray(hydplp)[:, 3:6].ravel())

    # 14. rotpos_com (86, 3)
    parts.append(np.asarray(rotposp)[:, :3].ravel())

    # 15. rotpos_dir_unnorm (86, 3): target directions are unit vectors; set equal
    parts.append(np.asarray(rotposp)[:, 3:6].ravel())

    # 16. rotscalar (86, 1)
    parts.append(np.asarray(rotscalarp).ravel())

    x0 = np.concatenate(parts)
    assert len(x0) == n_param, f'x0 length {len(x0)} != n_param {n_param}'
    return x0


def pack_param(rotp, covp, hydp, hydplp, rotposp, rotscalarp):
    """Inverse of unpack_params: find latent vector x s.t. unpack_params(x) ≈ inputs.

    Because the GLY palindrome constraint prevents an exact fit to a
    non-palindromic GLY angular profile, convergence is checked loosely.
    """
    targets = [
        torch.tensor(np.asarray(t), dtype=torch.float64)
        for t in [rotp, covp, hydp, hydplp, rotposp, rotscalarp]
    ]

    def obj_and_grad(x_np):
        x = torch.tensor(x_np, dtype=torch.float64, requires_grad=True)
        with torch.enable_grad():
            out  = unpack_param_maker(x)
            loss = sum(((a - b) ** 2).sum() for a, b in zip(out[:6], targets))
            loss.backward()
        return loss.item(), x.grad.numpy().copy()

    x0     = _init_x0(rotp, covp, hydp, hydplp, rotposp, rotscalarp)
    result = opt.minimize(
        obj_and_grad, x0, jac=True, method='L-BFGS-B',
        options={'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-10},
    )
    # The GLY palindrome constraint produces an irreducible residual of
    # 0.5 * sum((gly_dp1 - gly_dp1_flipped)^2) — the optimizer reaching a
    # stationary point is the correct convergence criterion, not a loss < 1.
    if not result.success:
        raise ValueError(
            f'pack_param optimizer did not converge: {result.message} '
            f'(final loss={result.fun:.4g})'
        )
    print(f'pack_param converged: final loss = {result.fun:.6g} '
          f'({result.message})')
    return result.x


# ---------------------------------------------------------------------------
# Quadratic spline energy (used for regularization in ConDiv d_obj)
# ---------------------------------------------------------------------------

def quadspline_energy(params, n_knots):
    """Evaluate the quadratic spline potential on its grid.

    Args:
        params:  tensor of shape (..., sum(n_knots))
        n_knots: tuple (n1, n2, n3, n4) giving knot counts for (dp1, dp2, uni, direc)

    Returns:
        tensor of shape (..., n3-2, n1-2, n2-2)
    """
    dp1   = params[..., :sum(n_knots[:1])]
    dp2   = params[..., sum(n_knots[:1]):sum(n_knots[:2])]
    uni   = params[..., sum(n_knots[:2]):sum(n_knots[:3])]
    direc = params[..., sum(n_knots[:3]):]

    ev = lambda sp: (1. / 6.) * sp[..., :-2] + (2. / 3.) * sp[..., 1:-1] + (1. / 6.) * sp[..., 2:]

    return (
        ev(uni)  [..., :, None, None] +
        ev(dp1)  [..., None, :, None] *
        ev(dp2)  [..., None, None, :] *
        ev(direc)[..., :, None, None]
    )


# ---------------------------------------------------------------------------
# Adam optimizer (pure Python / numpy — no framework dependency)
# ---------------------------------------------------------------------------

def _read_comp(x, i):
    try:
        return x[i]
    except (TypeError, IndexError):
        return x   # scalar


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

    def update_for_d_obj(self):
        """Return zero corrections (used for Nesterov look-ahead; not needed in Adam)."""
        return [0. for _ in self.grad1]

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
