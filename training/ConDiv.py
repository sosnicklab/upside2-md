#!/usr/bin/env python3
"""FF2 dual-target contrastive-divergence training for the Upside core force field.

This is ff2.1's own training workflow: Peng et al., JCTC 2022, SI "Parameterization by Contrastive
Divergence". The code is O. Kleinmann's Python 3 port of Peng's trainer
(condiv2.py, 2025-12-16, github nnamnielk/condiv4upside2),
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
  * the Ramachandran library is a fixed input, `upside_input/rama.dat`, never trained. For ff3.0
    it is parameters/common/rama31.dat, whose central-glycine row is fitted in Upside to BioEmu's
    plain-MD octapeptides (up.md 2.8, findings 1.19): the PDB row is part local energy and part
    evolutionary placement, and a map trained against native structures relearns the placement
    (findings 1.15).
  * two changes for ff3.0 (findings 1.17). Glycine has its own offsets on the three H-bond basin
    energies (hbond.h5 entries 12-14, field `hbg`), trained from zero with everything else: with
    the shared energies alone, keeping loop glycines left-handed costs every helical glycine its
    helix. The side-chain (`rot`) learning rate is 10x the port's smaller: at the port's rate
    Adam's normalised step random-walks the 31,420 pair coefficients, of which ~6% carry signal,
    and the walk weakens helices and folds.
  * recorded, not trained: per-residue basin populations of the native-restrained replica and of
    the free ensemble (basin_weights below), the diagnostic that shows where the free simulation
    leaves the native conformation.

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
  - per-residue basin populations of the native-restrained replica and of the free ensemble (the
    NSE mix), written to the step's divergence file.

Commands (`initialize` copies this file into the run's output directory, and every later command
runs that copy, so a run is never continued, judged or extracted by later code):
  ConDiv.py initialize <init_param_dir> <upside_input_dir> <pdb_list> <output_dir>
  ConDiv.py restart <checkpoint.pkl> <n_steps>
  ConDiv.py gate <run_dir>                      exit 0 converged, 3 not, anything else a failure
  ConDiv.py extract <checkpoint.pkl> <out_dir>  the force-field files of one checkpoint
"""

import sys
import os

sys.path.insert(0, os.path.join(os.environ['UPSIDE_HOME'], 'py'))

is_worker = __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'worker'

import collections
import glob
import itertools
import pickle as cp
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
# Workers run at once on a machine without Slurm. One worker keeps ~10 cores busy (measured on an
# M1 Ultra, 1ga3: 345 s wall, 3411 s CPU), so a 20-core machine runs two.
local_workers   = int(os.environ.get('CONDIV_LOCAL_WORKERS', '1'))

resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
hydrophobicity_order = ['ASP', 'GLU', 'LYS', 'HIS', 'ARG', 'GLY', 'ASN', 'GLN', 'ALA', 'SER',
                        'THR', 'PRO', 'CYS', 'VAL', 'MET', 'TYR', 'ILE', 'LEU', 'PHE', 'TRP']

Target = collections.namedtuple('Target', 'fasta native native_path init_path n_res chi')
FIELDS = 'enve envc envs envw bbenve bbenvc bbenvs bbenvw cov rot hyd hb dhb hbg sheet'.split()
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
    """hbond.h5 `parameter` -> (hb: entries 0-2 and 4-11, dhb: entry 3, as the port; hbg: the
    residue-class offsets, entries 12 on)."""
    parameter = np.asarray(parameter)
    return parameter[[0, 1, 2] + list(range(4, 12))], parameter[3:4], parameter[12:]


def hb_join(hb, dhb, hbg):
    out = np.zeros(12 + len(hbg))
    out[:3], out[3], out[4:12], out[12:] = hb[:3], dhb[0], hb[3:], hbg
    return out


# ---------------------------------------------------------------------------
# Ramachandran basins: a diagnostic, recorded per residue every step
# ---------------------------------------------------------------------------
# The basins partition the torus: alpha_R (phi < 0, -100 < psi < 50), beta (phi < -100, psi outside
# that band), pPII (-100 < phi < 0, psi outside it), alpha_L (the mirror of alpha_R), and phi > 0
# outside alpha_L split into beta' and pPII', the mirrors of beta and pPII; `other` is their union.
# Edges are logistic with a 3 deg scale (13 deg from 10% to 90%), continuous across phi = +-180,
# and every basin mirrors exactly onto its partner under (phi, psi) -> (-phi, -psi).

BASINS = ('alpha_R', 'alpha_L', 'beta', 'pPII', "beta'", "pPII'", 'other')
BASIN_EDGE = np.deg2rad(3.)


def _arc(x, a, b):
    """Smooth periodic indicator of the arc from a to b (radians, anticlockwise): a sigmoid of the
    circular distance from the arc's midpoint, so an arc mirrors exactly onto the arc (-b, -a)."""
    half = 0.5 * ((b - a) % (2. * np.pi))
    dist = (np.asarray(x) - (a + half) + np.pi) % (2. * np.pi) - np.pi
    return 1. / (1. + np.exp(-(half - np.abs(dist)) / BASIN_EDGE))


def basin_weights(phi, psi):
    """w_b(phi,psi) for every basin, shape (..., len(BASINS)). Angles in radians."""
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
    """Per-residue basin populations, (n_res, len(BASINS)), of an ensemble of (n_frame, n_res, 2)
    (phi, psi) in radians, with per-frame weights summing to one or None for a plain mean."""
    w = basin_weights(rama_coord[..., 0], rama_coord[..., 1])
    if weights is None:
        return w.mean(axis=0)
    return np.tensordot(np.asarray(weights, dtype=float), w, axes=1)


# ---------------------------------------------------------------------------
# Driver-side: parameters, files, regulariser
# ---------------------------------------------------------------------------

if not is_worker:

    def print_param(param):
        print('hb    %s   dhb %.4f   hbg %s' % (np.array2string(param.hb[:3], precision=4),
                                              param.dhb[0], np.array2string(param.hbg, precision=4)))
        print('sheet mean %.4f' % np.mean(param.sheet))
        print('bb env scale %.4f  center %.4f  sharpness %.4f  hbond_weight %.4f'
              % (param.bbenve, param.bbenvc, param.bbenvs, param.bbenvw))
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

    def get_init_param(init_dir, rama_library):
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
            hb, dhb, hbg = hb_split(t.root.parameter[:])

        param = Update(*([None] * len(FIELDS)))._replace(
            enve=enve, envc=envc, envs=envs, envw=envw,
            bbenve=bbe, bbenvc=bbc, bbenvs=bbs, bbenvw=bbw,
            rot=rp.pack_param(rotp, covp, hydp, hydplp, rotposp, np.zeros((rotposp.shape[0], 1))),
            hb=hb, dhb=dhb, hbg=hbg, sheet=np.loadtxt(files['sheet']))
        return param, files

    def expand_param(params, orig, new):
        """Write the trained parameters to the files a worker (or a deployment) reads. The rama
        library is not trained and is never rewritten."""
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
            t.root.parameter[:] = hb_join(params.hb, params.dhb, params.hbg)

        np.savetxt(new['sheet'], params.sheet)

    def backprop_deriv(param, deriv, reg_scale):
        return deriv._replace(
            rot=d_obj(param.rot, deriv.rot, deriv.cov, deriv.hyd, reg_scale),
            cov=0., hyd=0.)


# ---------------------------------------------------------------------------
# Minibatch
# ---------------------------------------------------------------------------

def run_minibatch(worker_path, param, init_files, direc, minibatch, solver, reg_scale):
    os.makedirs(direc, exist_ok=True)
    print(direc)

    # the workers read the fixed library; every other file is written for this step
    d_obj_files = dict((k, os.path.join(direc, 'nesterov_temp__' + os.path.basename(v)))
                       for k, v in init_files.items() if k != 'rama')
    d_obj_files['rama'] = init_files['rama']
    expand_param(param + solver.update_for_d_obj(), init_files, d_obj_files)

    has_slurm = shutil.which('srun') is not None
    files_arg = b64encode(cp.dumps(dict(files=d_obj_files))).decode('ascii')
    targets = dict(minibatch)

    def launch(nm):
        t = targets[nm]
        args = [sys.executable, worker_path, 'worker', nm, direc, t.fasta, t.init_path,
                str(t.n_res), t.chi, files_arg]
        if has_slurm:
            # srun's own messages go to <code>.srun, the worker's output to <code>.output_worker
            args = ['srun', '--nodes=1', '--ntasks=1', '--cpus-per-task=%i' % n_threads,
                    '--output=%s/%s.output_worker' % (direc, nm)] + args
            return sp.Popen(args, close_fds=True, stderr=open('%s/%s.srun' % (direc, nm), 'w'))
        return sp.Popen(args, close_fds=True, stdout=open('%s/%s.output_worker' % (direc, nm), 'w'),
                        stderr=sp.STDOUT)

    def never_launched(nm):
        path = '%s/%s.srun' % (direc, nm)
        return has_slurm and os.path.exists(path) and 'Task launch' in open(path).read()

    # Under Slurm every worker of the minibatch starts at once, srun placing each on its own cores.
    # Without it the workers share one machine, so at most `local_workers` run at a time.
    # A worker whose srun never started it ('Task launch ... failed', e.g. an expired job
    # credential when 24 steps start at once) has not run, so it is launched again as soon as that
    # is seen, up to twice; a worker that ran and failed is never relaunched.
    queue = [nm for nm, _ in minibatch[::-1]]
    limit = len(queue) if has_slurm else local_workers
    jobs, tries, rc = collections.OrderedDict(), dict(), dict()
    while len(rc) < len(minibatch):
        while queue and len(jobs) - len(rc) < limit:
            nm = queue.pop(0)
            jobs[nm], tries[nm] = launch(nm), 1
        for nm, j in list(jobs.items()):
            if nm in rc or j.poll() is None:
                continue
            if j.returncode != 0 and never_launched(nm) and tries[nm] < 3:
                print('%s never started (srun launch failed), relaunching' % nm)
                sys.stdout.flush()
                tries[nm] += 1
                jobs[nm] = launch(nm)
            else:
                rc[nm] = j.returncode
        time.sleep(5.)

    rmsd, change, no_dse, failed = dict(), [], [], []
    for nm in jobs:
        if rc[nm] != 0:
            print(nm, 'WORKER_FAIL')
            failed.append(nm)
            continue
        with open('%s/%s.divergence.pkl' % (direc, nm), 'rb') as f:
            div = cp.load(f)
        rmsd[nm] = (div['rmsd_restrain'], div['rmsd'])
        change.append(div['contrast'])
        if not div['has_dse']:
            no_dse.append(nm)
    # A step from part of the minibatch is a different objective, so it is never taken: the step
    # fails and the chain's successor repeats it from the last checkpoint.
    if failed:
        raise RuntimeError('%i of %i workers failed: %s' % (len(failed), len(jobs), ' '.join(failed)))

    d_param = backprop_deriv(
        param, Update(*[None if x[0] is None else np.sum(x, axis=0) for x in zip(*change)]),
        reg_scale)

    with open('%s/rmsd.pkl' % direc, 'wb') as f:
        cp.dump(rmsd, f, -1)
    print('\nMedian RMSD %.2f %.2f' % tuple(np.median(np.array(list(rmsd.values())), axis=0)))
    print('DSE target from %i of %i proteins%s' % (len(change) - len(no_dse), len(change),
          ('; no unfolded ensemble: ' + ' '.join(no_dse)) if no_dse else ''))

    new_files = dict((k, os.path.join(direc, os.path.basename(v)))
                     for k, v in init_files.items() if k != 'rama')
    new_param = param + solver.update_step(d_param)
    expand_param(new_param, init_files, new_files)
    solver.log_state(direc)
    print()
    print_param(new_param)
    return new_param


# ---------------------------------------------------------------------------
# Worker: one protein, one step
# ---------------------------------------------------------------------------

def zero_for_sarw(config, scale):
    """The port's apply_param_scale(hb, env, rot): H-bond energies (with any residue-class
    offsets), side-chain burial and rotamer pair energies scaled by `scale`, leaving backbone
    geometry, rama and the backbone term."""
    with tb.open_file(config, 'a') as t:
        p = t.root.input.potential
        p.hbond_energy.parameters[:4] *= scale
        p.hbond_energy.parameters[12:] *= scale
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

    Returns an Update of per-frame arrays and every residue's (phi,psi) per frame, from which the
    per-residue basin populations are formed after frame weighting.
    """
    config_base, traj, sheet_fd, start = args
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
    rama_coord = []
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
        hb, dhb, hbg = hb_split(engine.get_param_deriv(hb_n, 'hbond_energy').ravel())
        c.hb.append(hb)
        c.dhb.append(dhb)
        c.hbg.append(hbg)
        c.sheet.append(np.zeros(len(sheet_restype)))
        rama_coord.append(engine.get_output('rama_coord'))

    if sheet_fd:
        # Central difference on each present residue type's sheet mixing energy, as the port.
        for r in present:
            rid = sheet_restype.index(r)
            for arr, sign in ((more[r], +1.), (less[r], -1.)):
                engine.set_param(arr, 'rama_map_pot')
                for i, x in enumerate(pos):
                    engine.energy(x)
                    c.sheet[i][rid] += sign * engine.get_output('rama_map_pot')[0, 0] / (2. * eps)
    return Update(*[np.array(x) for x in c]), np.array(rama_coord)


def main_worker():
    tstart = time.time()
    code, direc, fasta, init_path, n_res, chi = sys.argv[2:8]
    n_res = int(float(n_res))
    files = cp.loads(b64decode(sys.argv[8].encode('ascii')))['files']
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

    jobs = [(configs[1], configs[i], True, start) for i in range(4)]
    if has_dse:
        jobs += [(configs[1], configs[r], False, start) for r in (rid_left, rid_right)]
        jobs += [(configs[1], configs[-1], False, start)]
    with Pool(processes=len(jobs)) as pool:
        div, coords = zip(*pool.map(compute_divergence, jobs))

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

    # The NSE's own two ensembles as per-residue basin populations, recorded for diagnosis.
    divergence = dict(
        contrast=Update(*contrast),
        rama_native=residue_populations(coords[0]),
        rama_free=residue_populations(np.concatenate(coords[1:4]), w),
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
                                       float(state['n_prot']))
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

    state['param'], state['init_param_files'] = get_init_param(
        init_dir, os.path.join(protein_dir, 'rama.dat'))
    with tb.open_file(state['init_param_files']['rama']) as t:
        origin = t.root._v_attrs.glycine_row if 'glycine_row' in t.root._v_attrs else b'as source'
    print('rama library (fixed): %s\n  glycine row: %s'
          % (state['init_param_files']['rama'],
             origin.decode() if isinstance(origin, bytes) else origin))

    # The port's learning rates, times its global factor, except rot, 10x smaller (see the top);
    # the glycine offsets take hb's. bbenvc/bbenvs/bbenvw are 0 because the engine returns no
    # derivative for them, as in ff2.1's training.
    state['initial_alpha'] = Update(
        enve=0.10, envc=0.05, envs=0.02, envw=0.10,
        bbenve=0.05, bbenvc=0.00, bbenvs=0.00, bbenvw=0.00,
        cov=0., rot=0.025, hyd=0., hb=0.02, dhb=0.01, hbg=0.02, sheet=0.03) * alpha_scale
    state['solver'] = rp.AdamSolver(len(FIELDS), alpha=state['initial_alpha'])
    state['epoch'], state['i_mb'] = 0, 0
    print('\nOptimizing with solver', state['solver'], '\n')
    print_param(state['param'])
    return state


# ---------------------------------------------------------------------------
# Convergence gate and extraction, run from the run's own copy
# ---------------------------------------------------------------------------

def main_gate(run_dir, alpha=0.05):
    """Is every trained group at a fixed point over the last full epoch? Exit 0 if so, 3 if not
    (including when there is less than one epoch to judge).

    The last full epoch is the window because every protein has then contributed exactly once, so
    the window-mean gradient is the full-training-set gradient. For each Adam group the raw
    per-step gradient is recovered from the Adam state, g_t = (grad1_t - b1*grad1_{t-1}) / (1-b1).
    At a fixed point each step's gradient is noise symmetric about zero, so ||sum_t g_t|| should be
    unremarkable among all 2^n sign flips of the steps; with K the Gram matrix of the steps,
    ||sum_t s_t g_t||^2 = s^T K s, so the permutation p-value is exact. Every group must have
    p > alpha / (number of groups). Parameter movement is not used: Adam's steps are
    scale-invariant.
    """
    out = os.path.join(os.path.abspath(run_dir), 'run_output')
    steps = sorted(d for d in glob.glob(os.path.join(out, 'epoch_*_minibatch_*'))
                   if os.path.exists(os.path.join(d, 'checkpoint.pkl')))
    state = cp.load(open(os.path.join(steps[-1], 'checkpoint.pkl'), 'rb'))
    n = len(state['minibatches'])
    if len(steps) < n:
        print('only %d steps, fewer than one epoch (%d); nothing to judge yet' % (len(steps), n))
        return 3
    b1 = state['solver'].beta1
    trained = [f for i, f in enumerate(FIELDS) if np.any(np.asarray(state['initial_alpha'][i]) != 0)
               and getattr(state['param'], f) is not None]

    def grad1(d):
        st = cp.load(open(os.path.join(d, 'solver_state.pkl'), 'rb'), encoding='latin1')
        return [np.asarray(st['grad1'][FIELDS.index(f)], dtype=float) for f in trained]

    prev = grad1(steps[-n - 1]) if len(steps) > n else [0. for _ in trained]
    G = {f: [] for f in trained}
    for d in steps[-n:]:
        g1 = grad1(d)
        for k, f in enumerate(trained):
            G[f].append(np.atleast_1d((g1[k] - b1 * prev[k]) / (1. - b1)).ravel())
        prev = g1

    signs = np.array([(1,) + s for s in itertools.product((1, -1), repeat=n - 1)], dtype=float)
    cut = alpha / len(trained)
    print('%s: steps %d-%d (last epoch), %d groups, pass if every p > %.4f\n'
          % (run_dir, len(steps) - n + 1, len(steps), len(trained), cut))
    print('   %-7s %7s  %10s  %8s' % ('group', 'size', 'p', 'verdict'))
    ok = True
    for f in trained:
        X = np.array(G[f])
        K = X @ X.T
        p = (np.einsum('pi,ij,pj->p', signs, K, signs) >= K.sum() * (1 - 1e-12)).mean()
        ok &= p > cut
        print('   %-7s %7d  %10.4f  %8s' % (f, X.shape[1], p, 'ok' if p > cut else 'PULLED'))
    print('\n' + ('CONVERGED: every group is consistent with a fixed point.' if ok else
                  'NOT CONVERGED: at least one group still has a systematic pull.'))
    return 0 if ok else 3


def main_extract(checkpoint, out):
    """Write the force field of a checkpoint: sidechain.h5, environment.h5, bb_env.dat, hbond.h5
    and sheet through expand_param, as every step writes them, plus the run's fixed rama.dat."""
    os.makedirs(out, exist_ok=True)
    state = cp.load(open(checkpoint, 'rb'))
    param, init = state['param'], state['init_param_files']
    new = dict(rot='sidechain.h5', env='environment.h5', bbenv='bb_env.dat', hb='hbond.h5',
               sheet='sheet', rama='rama.dat')
    new = dict((k, os.path.join(out, v)) for k, v in new.items())
    expand_param(param, init, new)
    shutil.copyfile(init['rama'], new['rama'])

    print('checkpoint %s\n  step %d, next epoch %d minibatch %d'
          % (checkpoint, state['solver'].step_num, state['epoch'], state['i_mb']))
    print_param(param)
    for k, v in sorted(new.items()):
        print('  %-6s %-16s %8.2f MB' % (k, os.path.basename(v), os.path.getsize(v) / 1e6))


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
    elif len(sys.argv) == 3 and sys.argv[1] == 'gate':
        raise SystemExit(main_gate(sys.argv[2]))
    elif len(sys.argv) == 4 and sys.argv[1] == 'extract':
        main_extract(os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]))
    else:
        raise SystemExit(__doc__[__doc__.index('Commands'):])
