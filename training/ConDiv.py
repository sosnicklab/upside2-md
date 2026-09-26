#!/usr/bin/env python3
"""FF2 dual-target contrastive-divergence training for the Upside core force field.

This is ff2.1's own training workflow: Peng et al., JCTC 2022, SI "Parameterization by Contrastive
Divergence". The code is O. Kleinmann's Python 3 port of Peng's trainer
(/project2/trsosnic/okleinmann/condiv/condiv2.py, 2025-12-16, github nnamnielk/condiv4upside2),
modernised and restored to the published protocol. Differences from that file, each for a reason:

  * modernised: the Theano regulariser is the same objective in torch; mdtraj is replaced by the
    equivalent numpy (C-alpha RMSD and unweighted Rg) on Upside's own output; site paths (conda,
    srun, LD_LIBRARY_PATH, rama files under /home/okleinmann) come from UPSIDE_HOME and env.sh.
  * protocol restored to the SI where the port had drifted: DSE weight lambda = 0.3 (port 0.0, so
    the port never used its second objective), 12 free replicas from T = 0.8 to 1.1 (port 6, up
    to ~0.97, where many proteins never unfold), 8000 time units (port 1000), minibatch 24 (port
    21), 4 epochs.
  * what the engine gives no gradient for stays at ff2.1, exactly as in ff2.1's training: the
    backbone term's center, sharpness and hbond weight (their derivatives are commented out in
    BackboneSigmoidCoupling::get_param_deriv, in master too) and hbond.h5 entries 4-11, the rama
    boundaries and sharpnesses. bb_env.dat still trains through its scale.
  * two numerical defects fixed: the replica reweighting exponent was E*(T0-Ti)/Ti, which is
    T0 times the correct E*(1/Ti - 1/T0); and a dE < -200 clamp stood in for a normalisation,
    now an exact max-subtracted one. The port's guard that silently dropped the DSE term when a
    replica's final energy exceeded 1000 is removed: a blown-up replica must fail, not vanish.
    For the same reason a failed worker fails the step; the port summed whatever returned.
  * added, switched by TRAIN_GLY: the central-glycine row of the Ramachandran library, 42 maps
    each trained on its own, GLY|GLY held mirror-symmetric (rama_gly_gradient.py). With it off
    the library is read unchanged, as in ff2.1's training.

One training step, per protein (main_worker):
  - 14 systems in one replica-exchange run: a native-restrained replica, 12 free replicas and one
    self-avoiding-random-walk (SARW) replica with H-bond, side-chain burial and rotamer pair
    energies scaled to zero. All start from the native structure.
  - NSE target: <dV/da>_native - <dV/da>_free, the free ensemble being the three coldest free
    replicas reweighted to T0 and mixed 0.6/0.3/0.1.
  - DSE target: <dV/da>_SARW - <dV/da>_unfolded, the unfolded ensemble being the frames of the two
    replicas bracketing the Rg midpoint (the Tm estimate) whose Rg exceeds
    0.67*Rg(coldest) + 0.33*Rg(hottest).
  - contrast = NSE + lambda * DSE, summed over the minibatch, fed to Adam.
"""

import sys
import os

# The script runs from its run directory's copy, next to the rama_gly_gradient.py it was
# initialised with, which therefore takes precedence over the one in training/.
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(os.environ['UPSIDE_HOME'], 'py'), _here):
    if _p not in sys.path:
        sys.path.insert(0, _p)

is_worker = __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'worker'

import collections
import pickle as cp
import re
import shutil
import socket
import subprocess as sp
import time
from base64 import b64decode, b64encode
from multiprocessing import Pool

import numpy as np
import tables as tb
import torch

import run_upside as ru
import upside_engine as ue
import upside_config as uc
import rama_gly_gradient as rgg

if not is_worker:
    import rotamer_parameter_estimation as rp

np.set_printoptions(precision=2, suppress=True)

# ---------------------------------------------------------------------------
# Protocol (SI values unless noted)
# ---------------------------------------------------------------------------
n_free          = 12                                     # SI: 12 replicas ...
T_free          = 0.8 * (1.1 / 0.8) ** (np.arange(n_free) / (n_free - 1.))   # ... from 0.8 to 1.1
n_system        = n_free + 2                             # + native-restrained + SARW
n_threads       = n_system                               # one core per system
native_restraint_strength = 1. / 3. ** 2                 # holds the native near 1.0 A RMSD
rmsd_k          = 10                                     # residues trimmed from each end for RMSD
minibatch_size  = 24                                     # SI: 456 proteins in 19 batches
n_epoch         = 4                                      # SI: 76 iterations, 4 cycles
sim_time        = 8000.                                  # SI: 8000 Upside time per step
frame_interval  = 10.                                    # the port's: 100 frames per 1000
equil_fraction  = 0.5                                    # SI: second half of each trajectory
dse_weight      = 0.3                                    # SI: lambda
sarw_scale      = 0.0                                    # SARW: H-bond, burial, rotamer off
nse_mix         = (0.60, 0.30, 0.10)                     # the port's mix of the 3 coldest replicas
alpha_scale     = 0.5                                    # the port's global learning-rate factor
GLY_FOURIER_ORDER = 8                                    # smoothing of each glycine map update

# False reproduces ff2.1's own workflow, which is run first from ff2.1 to check that every file
# updates and every group plateaus. True adds the GLY|GLY symmetry and trains the glycine row.
# It is read once, at initialisation, and kept in the run's state.
TRAIN_GLY = True

resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
hydrophobicity_order = ['ASP', 'GLU', 'LYS', 'HIS', 'ARG', 'GLY', 'ASN', 'GLN', 'ALA', 'SER',
                        'THR', 'PRO', 'CYS', 'VAL', 'MET', 'TYR', 'ILE', 'LEU', 'PHE', 'TRP']

Target = collections.namedtuple('Target', 'fasta native native_path init_path n_res chi')
FIELDS = 'enve envc envs envw bbenve bbenvc bbenvs bbenvw cov rot hyd hb dhb sheet gly'.split()
UpdateBase = collections.namedtuple('UpdateBase', FIELDS)


class Update(UpdateBase):
    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False
        if is_seq:
            assert len(self) == len(other)
            ret = [None if a is None or b is None else op(a, b) for a, b in zip(self, other)]
        else:
            ret = [None if a is None else op(a, other) for a in self]
        return Update(*ret)

    def __add__(self, other):     return self._do_binary(other, lambda a, b: a + b)
    def __sub__(self, other):     return self._do_binary(other, lambda a, b: a - b)
    def __mul__(self, other):     return self._do_binary(other, lambda a, b: a * b)
    def __truediv__(self, other): return self._do_binary(other, lambda a, b: a / b)


def hb_split(parameter):
    """hbond.h5 `parameter` (12) -> (hb: entries 0-2 and 4-11, dhb: entry 3), as the port."""
    parameter = np.asarray(parameter)
    return parameter[[0, 1, 2] + list(range(4, 12))], parameter[3:4]


def hb_join(hb, dhb):
    out = np.zeros(12)
    out[:3], out[3], out[4:] = hb[:3], dhb[0], hb[3:]
    return out


# ---------------------------------------------------------------------------
# Driver-side: parameters, files, regulariser
# ---------------------------------------------------------------------------

if not is_worker:

    def print_param(param, row):
        print('hb    %s   dhb %.4f' % (np.array2string(param.hb[:3], precision=4), param.dhb[0]))
        print('sheet mean %.4f' % np.mean(param.sheet))
        print('bb env scale %.4f  center %.4f  sharpness %.4f  hbond_weight %.4f'
              % (param.bbenve, param.bbenvc, param.bbenvs, param.bbenvw))
        if param.gly is not None:
            dg = rgg.handedness(param.gly)
            gg = param.gly[row['is_gg']]
            print('gly   X|GLY dG(aR->aL) mean %+.4f [%+.3f, %+.3f]   GLY|GLY dG %s   '
                  'asymmetry %.2e'
                  % (dg[~row['is_gg']].mean(), dg[~row['is_gg']].min(), dg[~row['is_gg']].max(),
                     np.array2string(dg[row['is_gg']], precision=4),
                     np.abs(gg - rgg.mirror(gg)).max()))
        print('env   scale center sharpness weight')
        env = dict(zip(resnames, np.vstack([param.enve, param.envc, param.envs,
                                            param.envw[:20]]).T))
        for r in hydrophobicity_order:
            print('   ', r, env[r])

    def _d_obj_fn():
        """The port's Theano regulariser and coupling objective, in torch."""
        def d_obj(lparam_np, rot_coupling_np, cov_coupling_np, hyd_coupling_np, reg_scale):
            lparam = torch.tensor(lparam_np, requires_grad=True, dtype=torch.float64)
            rot_p, cov_p, hyd_p, _, _, _, _ = rp.unpack_param_maker(lparam)

            hb_n = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_hb, rp.n_knot_hb)
            sc_n = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_sc, rp.n_knot_sc)
            rot_e = rp.quadspline_energy(rot_p, sc_n)
            cov_e = rp.quadspline_energy(cov_p, hb_n)
            hyd_e = rp.quadspline_energy(hyd_p, hb_n)

            def cutoff(x, scale):
                return 1. / (1. + torch.exp(x / scale))

            def lower_bound(x, lb):
                return torch.where(x < lb, (x - lb) ** 2, torch.zeros_like(x))

            def student_t(t, scale, nu):
                return 0.5 * (nu + 1.) * torch.log(1. + (1. / (nu * scale ** 2)) * t ** 2)

            r_rot = torch.tensor(np.arange(1, rp.n_knot_sc - 1) * rp.sc_dr)
            r_hb = torch.tensor(np.arange(1, rp.n_knot_hb - 1) * rp.hb_dr)
            rot_expect = (5. * cutoff(r_rot - 2., 0.2)).view(-1, 1, 1)
            cov_expect = (5. * cutoff(r_hb - 2., 0.2)).view(-1, 1, 1)
            hyd_expect = (5. * cutoff(r_hb - 1., 0.2)).view(-1, 1, 1)

            reg = (student_t(rot_e - rot_expect, 200., 3.).sum() / reg_scale +
                   student_t(cov_e - cov_expect, 200., 3.).sum() / reg_scale +
                   student_t(hyd_e - hyd_expect, 200., 3.).sum() / reg_scale +
                   lower_bound(rot_e, -6.).sum() +
                   lower_bound(cov_e, -6.).sum() +
                   lower_bound(hyd_e, -6.).sum())
            obj = ((torch.tensor(rot_coupling_np) * rot_p).sum() +
                   (torch.tensor(cov_coupling_np) * cov_p).sum() +
                   (torch.tensor(hyd_coupling_np) * hyd_p).sum() + reg)
            obj.backward()
            return lparam.grad.numpy()
        return d_obj

    d_obj = _d_obj_fn()

    def get_init_param(init_dir, rama_library, train_gly):
        files = dict(
            env   = os.path.join(init_dir, 'environment.h5'),
            bbenv = os.path.join(init_dir, 'bb_env.dat'),
            rot   = os.path.join(init_dir, 'sidechain.h5'),
            hb    = os.path.join(init_dir, 'hbond.h5'),
            sheet = os.path.join(init_dir, 'sheet'),
            rama  = rama_library,
        )
        with tb.open_file(files['rot']) as t:
            rotp, covp, hydp = (t.root.pair_interaction[:], t.root.coverage_interaction[:],
                                t.root.hydrophobe_interaction[:])
            hydplp, rotposp = t.root.hydrophobe_placement[:], t.root.rotamer_center_fixed[:]
        with tb.open_file(files['env']) as t:
            enve, envc, envs, envw = (t.root.scale[:], t.root.center[:], t.root.sharpness[:],
                                      t.root.weights[:])
        bbe, bbc, bbs, bbw = np.loadtxt(files['bbenv'])
        with tb.open_file(files['hb']) as t:
            hb, dhb = hb_split(t.root.parameter[:])
        if train_gly:
            keys, is_gg, gly = rgg.start_row(rama_library)
            row = dict(keys=keys, is_gg=is_gg)
        else:
            gly, row = None, None

        param = Update(*([None] * len(FIELDS)))._replace(
            enve=enve, envc=envc, envs=envs, envw=envw,
            bbenve=bbe, bbenvc=bbc, bbenvs=bbs, bbenvw=bbw,
            rot=rp.pack_param(rotp, covp, hydp, hydplp, rotposp, np.zeros((rotposp.shape[0], 1))),
            hb=hb, dhb=dhb, sheet=np.loadtxt(files['sheet']), gly=gly)
        return param, files, row

    def expand_param(params, orig, new, row):
        """Write params to the files a worker (or a deployment) reads."""
        rotp, covp, hydp, hydplp, rotposp, _ = rp.unpack_params(params.rot)
        shutil.copyfile(orig['rot'], new['rot'])
        with tb.open_file(new['rot'], 'a') as t:
            t.root.pair_interaction[:] = rotp
            t.root.coverage_interaction[:] = covp
            t.root.hydrophobe_interaction[:] = hydp
            t.root.hydrophobe_placement[:] = hydplp
            t.root.rotamer_center_fixed[:] = rotposp

        shutil.copyfile(orig['env'], new['env'])
        with tb.open_file(new['env'], 'a') as t:
            t.root.scale[:] = params.enve
            t.root.center[:] = params.envc
            t.root.sharpness[:] = params.envs
            t.root.weights[:] = params.envw

        np.savetxt(new['bbenv'], [params.bbenve, params.bbenvc, params.bbenvs, params.bbenvw])

        shutil.copyfile(orig['hb'], new['hb'])
        with tb.open_file(new['hb'], 'a') as t:
            t.root.parameter[:] = hb_join(params.hb, params.dhb)

        np.savetxt(new['sheet'], params.sheet)

        # The library is 35 MB, so it is written only where asked: for the workers when the
        # glycine row is trained, never for the per-step record, since param.gly is in every
        # checkpoint. Untrained, the workers read the original library directly.
        if 'rama' in new and params.gly is not None:
            rgg.write_row(orig['rama'], new['rama'], row['keys'], params.gly)

    def backprop_deriv(param, deriv, reg_scale, row):
        # Each glycine map's update is band-limited so the learned correction is smooth, and the
        # GLY|GLY maps are re-projected so they stay exactly mirror-symmetric. Untrained, the
        # component is a zero gradient on a None parameter, which Update carries as None.
        if deriv.gly is None:
            gly = 0.
        else:
            g = np.stack([rgg.fourier_lowpass(m, GLY_FOURIER_ORDER) for m in deriv.gly])
            gly = rgg.constrain(g, row['is_gg'])
        return deriv._replace(
            rot=d_obj(param.rot, deriv.rot, deriv.cov, deriv.hyd, reg_scale),
            cov=0., hyd=0., gly=gly)


# ---------------------------------------------------------------------------
# Minibatch
# ---------------------------------------------------------------------------

def run_minibatch(worker_path, param, init_files, direc, minibatch, solver, reg_scale, row):
    os.makedirs(direc, exist_ok=True)
    print(direc)

    train_gly = param.gly is not None
    d_obj_files = dict((k, os.path.join(direc, 'nesterov_temp__' + os.path.basename(v)))
                       for k, v in init_files.items())
    if not train_gly:
        d_obj_files['rama'] = init_files['rama']
    expand_param(param + solver.update_for_d_obj(), init_files, d_obj_files, row)

    has_slurm = shutil.which('srun') is not None
    files_arg = b64encode(cp.dumps(dict(files=d_obj_files, train_gly=train_gly))).decode('ascii')
    jobs = collections.OrderedDict()
    for nm, t in minibatch[::-1]:
        args = [sys.executable, worker_path, 'worker', nm, direc, t.fasta, t.init_path,
                str(t.n_res), t.chi, files_arg]
        if has_slurm:
            args = ['srun', '--nodes=1', '--ntasks=1', '--cpus-per-task=%i' % n_threads,
                    '--output=%s/%s.output_worker' % (direc, nm)] + args
            jobs[nm] = sp.Popen(args, close_fds=True)
        else:
            jobs[nm] = sp.Popen(args, close_fds=True,
                                stdout=open('%s/%s.output_worker' % (direc, nm), 'w'),
                                stderr=sp.STDOUT)

    rmsd, change, no_dse, failed = dict(), [], [], []
    for nm, j in jobs.items():
        if j.wait() != 0:
            print(nm, 'WORKER_FAIL')
            failed.append(nm)
            continue
        with open('%s/%s.divergence.pkl' % (direc, nm), 'rb') as f:
            div = cp.load(f)
        rmsd[nm] = (div['rmsd_restrain'], div['rmsd'])
        change.append(div['contrast'])
        if not div['has_dse']:
            no_dse.append(nm)
    if train_gly:
        os.remove(d_obj_files['rama'])     # every worker has exited
    # A step from part of the minibatch is a different objective, so it is never taken: the step
    # fails and the chain's successor repeats it from the last checkpoint.
    if failed:
        raise RuntimeError('%i of %i workers failed: %s' % (len(failed), len(jobs), ' '.join(failed)))

    d_param = backprop_deriv(
        param, Update(*[None if x[0] is None else np.sum(x, axis=0) for x in zip(*change)]),
        reg_scale, row)

    with open('%s/rmsd.pkl' % direc, 'wb') as f:
        cp.dump(rmsd, f, -1)
    print('\nMedian RMSD %.2f %.2f' % tuple(np.median(np.array(list(rmsd.values())), axis=0)))
    print('DSE target from %i of %i proteins%s' % (len(change) - len(no_dse), len(change),
          ('; no unfolded ensemble: ' + ' '.join(no_dse)) if no_dse else ''))

    new_files = dict((k, os.path.join(direc, os.path.basename(v)))
                     for k, v in init_files.items() if k != 'rama')
    new_param = param + solver.update_step(d_param)
    expand_param(new_param, init_files, new_files, row)
    solver.log_state(direc)
    print()
    print_param(new_param, row)
    return new_param


# ---------------------------------------------------------------------------
# Worker: one protein, one step
# ---------------------------------------------------------------------------

def zero_for_sarw(config, scale):
    """The port's apply_param_scale(hb, env, rot): H-bond energies, side-chain burial and rotamer
    pair energies scaled by `scale`, leaving backbone geometry, rama and the backbone term."""
    with tb.open_file(config, 'a') as t:
        p = t.root.input.potential
        p.hbond_energy.parameters[:4] *= scale
        p.sigmoid_coupling_environment.scale[:] *= scale
        p.hbond_coverage.interaction_param[:] *= scale
        p.hbond_coverage_hydrophobe.interaction_param[:] *= scale
        p.rotamer.pair_interaction.interaction_param[:] *= scale


def read_output(config, start):
    with tb.open_file(config) as t:
        return t.root.output.pos[start:, 0]


def ca(pos):
    return pos[..., 1::3, :]


def radius_of_gyration(pos):
    x = ca(pos)
    return np.sqrt(((x - x.mean(axis=-2, keepdims=True)) ** 2).sum(-1).mean(-1))


def compute_divergence(args):
    """Per-frame dV/da for every trained parameter, under the free (unrestrained) Hamiltonian.

    Returns an Update of per-frame arrays; `gly` holds the glycines' (phi,psi) per frame, since
    the glycine-map derivative is linear in the samples and is formed after frame weighting.
    """
    config_base, traj, sheet_fd, start, gly_idx = args
    pos = read_output(traj, start)
    with tb.open_file(config_base) as t:
        p = t.root.input.potential
        seq = [s.decode() if isinstance(s, bytes) else s for s in t.root.input.sequence[:]]
        sheet_restype = [r.decode() if isinstance(r, bytes) else r
                         for r in p.rama_map_pot._v_attrs.restype]
        eps = p.rama_map_pot._v_attrs.sheet_eps
        present = sorted(set('PRO' if s == 'CPR' else s for s in seq))
        more = {r: p.rama_map_pot['more_sheet_rama_pot_' + r][:] for r in present}
        less = {r: p.rama_map_pot['less_sheet_rama_pot_' + r][:] for r in present}
        hb_n = p.hbond_energy.parameters.shape
        rot_s = p.rotamer.pair_interaction.interaction_param.shape
        cov_s = p.hbond_coverage.interaction_param.shape
        hyd_s = p.hbond_coverage_hydrophobe.interaction_param.shape
        n_env = sum(p.sigmoid_coupling_environment[k].shape[0]
                    for k in ('scale', 'center', 'sharpness', 'weights'))
        n_type = p.sigmoid_coupling_environment.scale.shape[0]

    engine = ue.Upside(config_base)
    c = Update(*[[] for _ in FIELDS])
    for x in pos:
        engine.energy(x)
        c.rot.append(engine.get_param_deriv(rot_s, 'rotamer'))
        c.cov.append(engine.get_param_deriv(cov_s, 'hbond_coverage'))
        c.hyd.append(engine.get_param_deriv(hyd_s, 'hbond_coverage_hydrophobe'))
        env = engine.get_param_deriv((n_env,), 'sigmoid_coupling_environment').ravel()
        c.enve.append(env[:n_type])
        c.envc.append(env[n_type:2 * n_type])
        c.envs.append(env[2 * n_type:3 * n_type])
        c.envw.append(env[3 * n_type:])
        bb = engine.get_param_deriv((4,), 'bb_sigmoid_coupling_environment').ravel()
        for k, f in enumerate((c.bbenve, c.bbenvc, c.bbenvs, c.bbenvw)):
            f.append(bb[k])
        hb, dhb = hb_split(engine.get_param_deriv(hb_n, 'hbond_energy').ravel())
        c.hb.append(hb)
        c.dhb.append(dhb)
        c.sheet.append(np.zeros(len(sheet_restype)))
        c.gly.append(engine.get_output('rama_coord')[gly_idx])

    if sheet_fd:
        # Central difference on each present residue type's sheet mixing energy, as the port.
        for r in present:
            rid = sheet_restype.index(r)
            for arr, sign in ((more[r], +1.), (less[r], -1.)):
                engine.set_param(arr, 'rama_map_pot')
                for i, x in enumerate(pos):
                    engine.energy(x)
                    c.sheet[i][rid] += sign * engine.get_output('rama_map_pot')[0, 0] / (2. * eps)
    return Update(*[np.array(x) for x in c])


def glycine_contrast(seq, files, div, w, sel, gly_idx):
    """The glycine row's contrast, with the same ensembles, signs and lambda as every other field.

    The map derivative is linear in the samples, so each ensemble's per-residue derivative is a
    weighted spline histogram of its (phi,psi); they are combined first and pushed through the
    left/right and coil/sheet mixtures once.
    """
    gi = FIELDS.index('gly')
    keys, _ = rgg.row_keys(files['rama'])
    G = rgg.read_row(files['rama'], keys)
    chain = rgg.GlycineMapChain(seq, files['rama'], files['sheet'], keys)
    n_grid = G.shape[-1]
    card = rgg.cardinal_function(n_grid)
    n_nat = len(div[0][gi])
    terms = [(div[0][gi], np.full(n_nat, 1. / n_nat), 1.),
             (np.concatenate([div[i][gi] for i in (1, 2, 3)]), w, -1.)]
    if sel is not None:
        n_s = len(div[6][gi])
        u = np.concatenate([div[4][gi][sel[0]], div[5][gi][sel[1]]])
        terms += [(div[6][gi], np.full(n_s, 1. / n_s), dse_weight),
                  (u, np.full(len(u), 1. / len(u)), -dse_weight)]
    res_grad = {}
    for j_, r in enumerate(gly_idx):
        res_grad[r] = sum(sign * rgg.map_gradient(coords[:, j_], n_grid, card, weights)
                          for coords, weights, sign in terms)
    return chain.backprop(G, res_grad)


def main_worker():
    tstart = time.time()
    code, direc, fasta, init_path, n_res, chi = sys.argv[2:8]
    n_res = int(float(n_res))
    payload = cp.loads(b64decode(sys.argv[8].encode('ascii')))
    files, train_gly = payload['files'], payload['train_gly']
    n_frame = int(sim_time / frame_interval)
    start = int(equil_fraction * n_frame)

    init = cp.load(open(init_path, 'rb'), encoding='latin1') if init_path.endswith('.pkl') \
        else np.load(init_path)
    init = init[:, :, 0] if init.ndim == 3 else init
    init_npy = os.path.join(direc, code + '.initial.npy')
    np.save(init_npy, init)

    kwargs = dict(
        environment_potential=files['env'], environment_potential_type=1,
        bb_environment_potential=files['bbenv'],
        rotamer_interaction=files['rot'], rotamer_placement=files['rot'],
        dynamic_rotamer_1body=True,
        hbond_energy=files['hb'],
        rama_sheet_mix_energy=files['sheet'], rama_param_deriv=True,
        rama_library=files['rama'],
        reference_state_rama=os.path.join(os.environ['UPSIDE_HOME'], 'parameters', 'common',
                                          'rama_reference.pkl'),
        initial_structure=init_npy,
    )
    # restrained native at T0, 12 free replicas, SARW at the hottest temperature
    T = np.concatenate([T_free[:1], T_free, T_free[-1:]])
    config_base = os.path.abspath('%s/%s.base.up' % (direc, code))
    configs = [config_base.replace('.base.up', '.run.%i.up' % i) for i in range(n_system)]
    ru.upside_config(fasta, config_base, **kwargs)
    for c_ in configs:
        shutil.copyfile(config_base, c_)
    ru.advanced_config(configs[0], restraint_groups=['0-%i' % (n_res - 1)],
                       restraint_spring_constant=native_restraint_strength)
    zero_for_sarw(configs[-1], sarw_scale)

    # The port's schedule: start at 5% of T, relax, anneal up to T between t = 96 and 400, then
    # sample at T. Swaps run among the restrained and free replicas; the SARW replica is apart.
    j = ru.run_upside('', configs, sim_time, frame_interval, n_threads=n_threads,
                      temperature=T * 0.05, swap_sets=ru.swap_table2d(n_system - 1, 1),
                      mc_interval=5., replica_interval=5., time_step=0.015,
                      anneal_factor=20., anneal_start=96., anneal_end=400.)
    if j.job.wait() != 0:
        raise RuntimeError('RUN_FAIL')

    seq = uc.read_fasta(open(fasta))
    gly_idx = [i for i, s in enumerate(seq) if s == 'GLY'] if train_gly else []
    trajs = [read_output(c_, start) for c_ in configs]
    rg = [radius_of_gyration(x) for x in trajs]
    rmsd = [ru.traj_rmsd(ca(x)[:, rmsd_k:-rmsd_k], ca(init)[rmsd_k:-rmsd_k]) for x in trajs[:2]]

    # Tm estimate: the pair of free replicas bracketing the midpoint of the mean-Rg ladder.
    mRg = np.array([r.mean() for r in rg[1:n_free + 1]])
    has_dse = mRg[-1] > mRg[0]
    if has_dse:
        mid_Rg = 0.5 * (mRg[0] + mRg[-1])
        left_Rg = 0.67 * mRg[0] + 0.33 * mRg[-1]
        rid_left = np.where(mRg < mid_Rg)[0][-1] + 1          # config index
        rid_right = rid_left + 1
        sel = [np.where(rg[r] > left_Rg)[0] for r in (rid_left, rid_right)]
        has_dse = sum(len(s) for s in sel) > 0

    jobs = [(configs[1], configs[i], True, start, gly_idx) for i in range(4)]
    if has_dse:
        jobs += [(configs[1], configs[r], False, start, gly_idx) for r in (rid_left, rid_right)]
        jobs += [(configs[1], configs[-1], False, start, gly_idx)]
    with Pool(processes=len(jobs)) as pool:
        div = pool.map(compute_divergence, jobs)

    # NSE free ensemble: replicas 1-3 each reweighted to T0 exactly, then mixed 0.6/0.3/0.1.
    engine = ue.Upside(configs[1])
    w = []
    for k, i in enumerate((1, 2, 3)):
        E = np.array([engine.energy(x) for x in trajs[i]])
        logw = E * (1. / T[i] - 1. / T[0])
        wi = np.exp(logw - logw.max())
        w.append(nse_mix[k] * wi / wi.sum())
    w = np.concatenate(w)

    # Ensembles as (per-frame arrays, frame weights summing to one), per field.
    def mean(d):
        return [x.mean(axis=0) for x in d]

    native = mean(div[0])
    free = [np.tensordot(w, np.concatenate([div[i][f] for i in (1, 2, 3)]), axes=1)
            for f in range(len(FIELDS))]
    contrast = [a - b for a, b in zip(native, free)]
    if has_dse:
        unf = [np.concatenate([div[4][f][sel[0]], div[5][f][sel[1]]]).mean(axis=0)
               for f in range(len(FIELDS))]
        sarw = mean(div[6])
        contrast = [c_ + dse_weight * (s - u) for c_, s, u in zip(contrast, sarw, unf)]

    # Glycine rows: the map derivative is linear in the frame samples, so each ensemble's
    # per-residue derivative is a weighted 2D histogram, combined with the same signs and
    # lambda as above and pushed through the mixture once.
    gi = FIELDS.index('gly')
    contrast[gi] = None
    if train_gly:
        contrast[gi] = glycine_contrast(seq, files, div, w, sel if has_dse else None, gly_idx)

    divergence = dict(
        contrast=Update(*contrast),
        rmsd_restrain=rmsd[0].mean(), rmsd=rmsd[1].mean(),
        has_dse=bool(has_dse), mean_rg=mRg, walltime=time.time() - tstart)

    for fn in [config_base] + configs + [init_npy]:
        os.remove(fn)
    with open('%s/%s.divergence.pkl' % (direc, code), 'wb') as f:
        cp.dump(divergence, f, -1)


# ---------------------------------------------------------------------------
# Main loop and initialisation
# ---------------------------------------------------------------------------

def main_loop(state_bytes, max_iter):
    for _ in range(max_iter):
        state = cp.loads(state_bytes)
        print('#########################################')
        print('####      EPOCH %2i MINIBATCH %2i      ####' % (state['epoch'], state['i_mb']))
        print('#########################################\n')
        sys.stdout.flush()
        tstart = time.time()
        state['mb_direc'] = os.path.join(state['base_dir'], 'epoch_%02i_minibatch_%02i'
                                         % (state['epoch'], state['i_mb']))
        if os.path.exists(state['mb_direc']):
            shutil.rmtree(state['mb_direc'])
        state['param'] = run_minibatch(state['worker_path'], state['param'],
                                       state['init_param_files'], state['mb_direc'],
                                       state['minibatches'][state['i_mb']], state['solver'],
                                       float(state['n_prot']), state['row'])
        print('\n%.0f seconds elapsed this minibatch' % (time.time() - tstart))
        sys.stdout.flush()
        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1
        state_bytes = cp.dumps(state, -1)
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'), 'wb') as f:
            f.write(state_bytes)


def main_initialize(init_dir, protein_dir, protein_list, base_dir):
    # absolute, so a checkpoint's file paths resolve from any working directory
    init_dir, protein_dir, base_dir = (os.path.abspath(x) for x in (init_dir, protein_dir, base_dir))
    os.makedirs(base_dir, exist_ok=True)
    state = dict(init_dir=init_dir, base_dir=base_dir)
    # the checkpoint carries the exact code that produced it
    state['worker_path'] = os.path.join(base_dir, 'ConDiv.py')
    shutil.copy(__file__, state['worker_path'])
    shutil.copy(os.path.join(_here, 'rama_gly_gradient.py'), base_dir)

    names = [x.split()[0] for x in open(protein_list)]
    assert names[0] == 'prot'
    training_set, excluded = dict(), []
    for code in sorted(names[1:]):
        base = os.path.join(protein_dir, code)
        native = cp.load(open(base + '.initial.pkl', 'rb'), encoding='latin1')[:, :, 0]
        if np.sqrt((np.diff(native, axis=0) ** 2).sum(-1)).max() < 2.:
            training_set[code] = Target(base + '.fasta', native, base + '.initial.pkl',
                                        base + '.initial.pkl', len(native) // 3, base + '.chi')
        else:
            excluded.append(code)
    print('Excluded %i proteins due to chain breaks' % len(excluded))

    tl = sorted(training_set.items(), key=lambda x: (x[1].n_res, x[0]))
    np.random.shuffle(tl)
    tl = tl[:len(tl) - len(tl) % minibatch_size]
    n_mb = len(tl) // minibatch_size
    state['minibatches'] = [tl[i::n_mb] for i in range(n_mb)]
    state['n_prot'] = n_mb * minibatch_size
    print('Constructed %i minibatches of size %i (%i proteins)'
          % (n_mb, minibatch_size, state['n_prot']))

    state['train_gly'] = TRAIN_GLY
    state['param'], state['init_param_files'], state['row'] = get_init_param(
        init_dir, os.path.join(protein_dir, 'rama.dat'), TRAIN_GLY)

    # The port's learning rates, times its global factor. bbenvc/bbenvs/bbenvw are 0 because the
    # engine returns no derivative for them, as in ff2.1's training. gly is 0.01 so that after
    # the factor each map moves ~0.005 nats per step.
    state['initial_alpha'] = Update(
        enve=0.10, envc=0.05, envs=0.02, envw=0.10,
        bbenve=0.05, bbenvc=0.00, bbenvs=0.00, bbenvw=0.00,
        cov=0., rot=0.25, hyd=0., hb=0.02, dhb=0.01, sheet=0.03, gly=0.01) * alpha_scale
    state['solver'] = rp.AdamSolver(len(FIELDS), alpha=state['initial_alpha'])
    state['epoch'], state['i_mb'] = 0, 0
    print('\nOptimizing with solver', state['solver'], '\n')
    print_param(state['param'], state['row'])
    return state


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'worker':
        main_worker()
    elif len(sys.argv) == 4 and sys.argv[1] == 'restart':
        print('Running as PID %i on host %s' % (os.getpid(), socket.gethostname()))
        main_loop(open(sys.argv[2], 'rb').read(), int(sys.argv[3]))
    elif len(sys.argv) == 6 and sys.argv[1] == 'initialize':
        state = main_initialize(*sys.argv[2:])
        with open(os.path.join(state['base_dir'], 'initial_checkpoint.pkl'), 'wb') as f:
            cp.dump(state, f, -1)
    else:
        print('Usage:')
        print('  ConDiv.py initialize <init_param_dir> <upside_input_dir> <pdb_list> <output_dir>')
        print('  ConDiv.py restart <checkpoint.pkl> <n_steps>')
        raise SystemExit(1)
