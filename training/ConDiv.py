#!/usr/bin/env python3
"""Contrastive divergence training for the Upside core force field.

Trains rotation (pair_interaction / coverage / hydrophobe / placement),
environment, hbond and sheet, which is the same set the Theano original trained.

The two scalars needed remapping because their nodes changed shape since the
original.  hbond went from a single energy to a 12-entry table, so hb is now a
multiplicative scale on parameters[:4]; the potential is exactly linear in those
four (verified to 7e-7 against the engine), so dE/ds = E/s as before.  Sheet
mixing went from one value to one per residue type, so sheet is a common offset
added to all of them, differenced against the more_/less_sheet_rama_pot_ALL pair.

No GLY symmetry correction is applied: rotamer_parameter_estimation.py is a
strict modernization of the Theano original, so the glycine angular profile is
learned free.  This run restarts from ff_2.1 with ff_2.1's own (unsymmetrized)
rama.dat to test whether the port reproduces the original's fixed point.
"""

import sys
import os

# ---------------------------------------------------------------------------
# Path setup — must come before any project-module imports
# ---------------------------------------------------------------------------
_project_py = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'py')
if _project_py not in sys.path:
    sys.path.insert(0, _project_py)

import torch  # noqa: F401 — verifies torch is available before heavy imports

is_worker = (
    __name__ == '__main__'
    and len(sys.argv) > 1
    and sys.argv[1] == 'worker'
)

import numpy as np
import subprocess as sp
import tables as tb
import pickle as cp
import shutil
import collections
import time
import socket

import run_upside as ru
import upside_engine as ue

if not is_worker:
    import rotamer_parameter_estimation as rp
import rama_gly_gradient as rgg
import upside_config as uc

np.set_printoptions(precision=2, suppress=True)

# ---------------------------------------------------------------------------
# Training hyper-parameters
# ---------------------------------------------------------------------------
n_threads           = int(os.environ.get('UPSIDE_TRAIN_NTHREADS', '8'))
native_restraint_k  = 1. / 3.**2   # hold native within ~1.0 Å RMSD
rmsd_k              = 15            # residues to trim from each end for RMSD
minibatch_size      = 12

resnames = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL',
]

hydrophobicity_order = [
    'ASP', 'GLU', 'LYS', 'HIS', 'ARG', 'GLY', 'ASN', 'GLN',
    'ALA', 'SER', 'THR', 'PRO', 'CYS', 'VAL', 'MET', 'TYR',
    'ILE', 'LEU', 'PHE', 'TRP',
]

Target     = collections.namedtuple('Target',     'fasta native native_path init_path n_res chi')
UpdateBase = collections.namedtuple('UpdateBase', 'env cov rot hyd hb sheet gly')


class Update(UpdateBase):
    def __init__(self, *args, **kwargs):
        super(Update, self).__init__()

    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False

        if is_seq:
            assert len(self) == len(other)
            ret = []
            for a, b in zip(self, other):
                if a is None or b is None:
                    ret.append(None)
                else:
                    ret.append(op(a, b))
        else:
            ret = [None if a is None else op(a, other) for a in self]
        return Update(*ret)

    def __add__(self, other):  return self._do_binary(other, lambda a, b: a + b)
    def __sub__(self, other):  return self._do_binary(other, lambda a, b: a - b)
    def __mul__(self, other):  return self._do_binary(other, lambda a, b: a * b)
    def __truediv__(self, other): return self._do_binary(other, lambda a, b: a / b)


# ---------------------------------------------------------------------------
# Non-worker utilities
# ---------------------------------------------------------------------------

if not is_worker:

    def print_param(param):
        print('hb    %.6f  (scale on hbond energies)' % param.hb)
        print('sheet %.6f  (common offset on sheet mixing)' % param.sheet)
        S, A = param.gly
        print('gly   dG(aR->aL) %+.4f nats   |A| rms %.5f   GLY|GLY asymmetry %.2e'
              % (_gly_handedness(S + A), np.sqrt((A ** 2).mean()),
                 np.abs(S - rgg.mirror(S)).max()))
        print('env')
        env_dict = dict(zip(resnames, param.env[:, 1::2]))
        for r in hydrophobicity_order:
            print('   ', r, env_dict[r])

    def _d_obj_fn():
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
                return torch.where(x < lb, 1e0 * (x - lb) ** 2, torch.zeros_like(x))

            def student_t(t, scale, nu):
                return 0.5 * (nu + 1.) * torch.log(1. + (1. / (nu * scale**2)) * t**2)

            r_rot = torch.tensor(np.arange(1, rp.n_knot_sc - 1) * rp.sc_dr)
            r_hb  = torch.tensor(np.arange(1, rp.n_knot_hb - 1) * rp.hb_dr)

            # view(-1, 1, 1): distance axis is dim 2 in (*_e); 1D tensor broadcasts to last dim by default
            rot_expect = (5. * cutoff(r_rot - 2., 0.2)).view(-1, 1, 1)
            hb_expect  = (5. * cutoff(r_hb  - 2., 0.2)).view(-1, 1, 1)
            hyd_expect = (5. * cutoff(r_hb  - 1., 0.2)).view(-1, 1, 1)

            reg = (
                student_t(rot_e - rot_expect, 200., 3.).sum() / reg_scale +
                student_t(cov_e - hb_expect,  200., 3.).sum() / reg_scale +
                student_t(hyd_e - hyd_expect, 200., 3.).sum() / reg_scale +
                lower_bound(rot_e, -6.).sum() +
                lower_bound(cov_e, -6.).sum() +
                lower_bound(hyd_e, -6.).sum()
            )

            obj = (
                (torch.tensor(rot_coupling_np) * rot_p).sum() +
                (torch.tensor(cov_coupling_np) * cov_p).sum() +
                (torch.tensor(hyd_coupling_np) * hyd_p).sum() +
                reg
            )
            obj.backward()
            return lparam.grad.numpy()

        return d_obj

    d_obj = _d_obj_fn()

    # Fourier truncation order for the glycine map update.  Modes with |kx|,|ky| <= this survive;
    # on a 72-node grid that is features down to 360/(2*8) = 22 degrees, which is finer than any
    # Ramachandran basin and far smoother than per-node noise.
    GLY_FOURIER_ORDER = 8

    def _gly_handedness(m):
        """dG(aR->aL) on the same basins every other glycine number in this project uses."""
        n = m.shape[0]
        t = np.arange(-180., 180., 360. / n)
        phi, psi = np.meshgrid(t, t, indexing='ij')
        def basin(a, b, c, d):
            k = (phi >= a) & (phi <= b) & (psi >= c) & (psi <= d)
            return -np.log(np.exp(-m[k]).sum())
        return basin(40., 100., -10., 60.) - basin(-100., -40., -60., 10.)

    def get_init_param(init_dir, rama_library):
        init_param_files = dict(
            env   = os.path.join(init_dir, 'environment.h5'),
            rot   = os.path.join(init_dir, 'sidechain.h5'),
            hb    = os.path.join(init_dir, 'hbond'),
            sheet = os.path.join(init_dir, 'sheet'),
            rama  = rama_library,
        )

        with tb.open_file(init_param_files['rot']) as t:
            rotp      = t.root.pair_interaction[:]
            covp      = t.root.coverage_interaction[:]
            hydp      = t.root.hydrophobe_interaction[:]
            hydplp    = t.root.hydrophobe_placement[:]
            rotposp   = t.root.rotamer_center_fixed[:]
            rotscalarp = np.zeros(rotposp.shape[0])

        with tb.open_file(init_param_files['env']) as t:
            env = t.root.energies[:, :-1]

        # hb is a multiplicative scale on hbond_energy.parameters[:4] and sheet is a common
        # additive offset on every entry of the sheet mixing file, so both start at their
        # identity values and the init files below are exactly ff_2.1.
        # gly is the pair (S, A): the central-glycine coil map split into its mirror-symmetric
        # and antisymmetric parts, with X|GLY = S + A and GLY|GLY = S.  Training starts at A = 0
        # on a symmetrised library row, so the handedness is entirely learned and the comparison
        # against the AWH measurement is not circular.
        gly = np.stack(rgg.symmetric_start(rama_library))

        param = Update(*([None] * 7))._replace(
            env   = env,
            rot   = rp.pack_param(rotp, covp, hydp, hydplp, rotposp, rotscalarp[:, None]),
            hb    = 1.0,
            sheet = 0.0,
            gly   = gly,
        )
        return param, init_param_files

    def expand_param(params, orig_param_files, new_param_files):
        """Write current params to the files the worker will read."""
        rotp, covp, hydp, hydplp, rotposp, _ = rp.unpack_params(params.rot)

        shutil.copyfile(orig_param_files['rot'], new_param_files['rot'])
        with tb.open_file(new_param_files['rot'], 'a') as t:
            t.root.pair_interaction[:]       = rotp
            t.root.coverage_interaction[:]   = covp
            t.root.hydrophobe_interaction[:] = hydp
            t.root.hydrophobe_placement[:]   = hydplp
            t.root.rotamer_center_fixed[:]   = rotposp

        shutil.copyfile(orig_param_files['env'], new_param_files['env'])
        with tb.open_file(new_param_files['env'], 'a') as t:
            tmp = np.zeros(t.root.energies.shape)
            tmp[:, :-1] = params.env
            tmp[:, -1]  = tmp[:, -3]
            t.root.energies[:] = tmp

        # hbond: scale the four energy entries.  The potential is exactly linear in
        # parameters[:4] (E_alpha, E_beta, E_other, E_bias), verified to 7e-7 against the engine,
        # which is what makes dE/ds = E/s in compute_divergence correct.  Entries 4..11 are the
        # rama boundary and sharpness values and are not energies, so they are left alone.
        shutil.copyfile(orig_param_files['hb'], new_param_files['hb'])
        with tb.open_file(new_param_files['hb'], 'a') as t:
            hbp = t.root.parameter[:]
            hbp[:4] = hbp[:4] * params.hb
            t.root.parameter[:] = hbp

        # sheet: one common offset added to all 20 per-residue-type mixing energies, which is the
        # single scalar the Theano original trained.
        np.savetxt(new_param_files['sheet'],
                   np.loadtxt(orig_param_files['sheet']) + params.sheet)

        # gly: rewrite the rama library's central-glycine coil row from the current (S, A).
        # Nothing else in the 35 MB file changes, and it is not renormalised, because
        # read_rama_maps_and_weights normalises the mixture itself and a shift here would move
        # the inner left/right weights in a way the gradient does not model.
        # Only written when asked.  A library is 35 MB, and keeping one per minibatch would add
        # 35 GB to run_output over 500 steps for no benefit: param.gly is in every checkpoint, so
        # the library can be regenerated with write_gly_library whenever it is actually wanted.
        if 'rama' in new_param_files:
            rgg.write_gly_library(orig_param_files['rama'], new_param_files['rama'],
                                  params.gly[0], params.gly[1])

    def backprop_deriv(param, deriv_update, reg_scale):
        envd = deriv_update.env[:, :-1].copy()
        envd[:, -2] += deriv_update.env[:, -1]
        # The raw map gradient is already spread over a 4x4 node neighbourhood by the spline
        # basis, but 5,184 free values against ~30,000 glycine samples per minibatch is still
        # noisy.  Band-limit the update so the learned correction is smooth by construction; the
        # starting map keeps its sharp forbidden-region structure because only the update is
        # filtered.  Re-project the symmetry afterwards, since S must stay symmetric for
        # GLY|GLY to stay achiral.
        gS, gA = deriv_update.gly
        gS = rgg.project_symmetric(rgg.fourier_lowpass(gS, GLY_FOURIER_ORDER))
        gA = rgg.project_antisymmetric(rgg.fourier_lowpass(gA, GLY_FOURIER_ORDER))
        return deriv_update._replace(
            rot   = d_obj(param.rot, deriv_update.rot, deriv_update.cov,
                          deriv_update.hyd, reg_scale),
            cov   = 0.,
            hyd   = 0.,
            env   = envd,
            gly   = np.stack([gS, gA]),
        )


# ---------------------------------------------------------------------------
# Divergence computation (run inside worker process)
# ---------------------------------------------------------------------------

def compute_divergence(config_base, pos, hb_scale, gly_chain, gly_param):
    try:
        with tb.open_file(config_base) as t:
            rot_shape  = t.root.input.potential.rotamer.pair_interaction.interaction_param.shape
            cov_shape  = t.root.input.potential.hbond_coverage.interaction_param.shape
            hyd_shape  = t.root.input.potential.hbond_coverage_hydrophobe.interaction_param.shape
            env_shape  = t.root.input.potential.nonlinear_coupling_environment.coeff.shape
            # The engine's parameter vector for this node is coeff THEN weights (verified
            # 2026-09-19 against get_param: 360 + 400 = 760).  Asking for coeff.shape alone gets
            # "Wrong number of parameters, expected 760 but got 360" and kills every worker --
            # which is how the 2026-09-10 libupside.so rebuild silently broke all training.
            # Request the full vector and keep the coeff slice; weights were never trained.
            env_nw     = t.root.input.potential.nonlinear_coupling_environment.weights.shape[0]
            more_sheet = t.root.input.potential.rama_map_pot.more_sheet_rama_pot_ALL[:]
            less_sheet = t.root.input.potential.rama_map_pot.less_sheet_rama_pot_ALL[:]
            sheet_eps  = t.root.input.potential.rama_map_pot._v_attrs.sheet_eps
    except Exception as e:
        print(os.path.basename(config_base)[:5], 'ANALYSIS_FAIL', e)
        return None

    sheet_scale = 1. / (2. * sheet_eps)

    engine   = ue.Upside(config_base)
    contrast = Update([], [], [], [], [], [], [])

    # dE/d(map) is analytic: rama_map_pot is a periodic interpolating bicubic spline and
    # solve_periodic_2d_spline is a tensor product of 1D solves, so the map enters the energy
    # linearly and separably and the derivative is a spline-smoothed 2D histogram of the glycine
    # (phi,psi) samples.  Finite differencing 5,184 map values would cost 10,369x a divergence.
    n_grid   = gly_chain.n_grid
    cardinal = rgg.cardinal_function(n_grid)
    gly_idx  = [r['index'] for r in gly_chain.residues]

    # Central difference on the common sheet mixing offset, as in the Theano original: evaluate
    # every frame with the map shifted +eps, then again with -eps.  The other derivatives are
    # collected in the +eps pass, which is also what the original did; eps is 5e-4 so the bias
    # that introduces is negligible, and matching the original's structure is the point here.
    engine.set_param(more_sheet, 'rama_map_pot')
    for i in range(pos.shape[0]):
        engine.energy(pos[i])
        contrast.rot.append(engine.get_param_deriv(rot_shape, 'rotamer'))
        contrast.cov.append(engine.get_param_deriv(cov_shape, 'hbond_coverage'))
        contrast.hyd.append(engine.get_param_deriv(hyd_shape, 'hbond_coverage_hydrophobe'))
        env_full = engine.get_param_deriv((int(np.prod(env_shape)) + env_nw,),
                                          'nonlinear_coupling_environment')
        contrast.env.append(env_full.ravel()[:int(np.prod(env_shape))].reshape(env_shape))
        # Logarithmic derivative of a scale factor: the hbond potential is exactly linear in
        # parameters[:4], so dE/ds = E/s.
        contrast.hb.append(engine.get_output('hbond_energy')[0, 0] / hb_scale)
        contrast.sheet.append(engine.get_output('rama_map_pot')[0, 0] * sheet_scale)
        if gly_idx:
            rc = engine.get_output('rama_coord')
            res_grad = {j: rgg.map_gradient(rc[j][None, :], n_grid, cardinal) for j in gly_idx}
            contrast.gly.append(np.stack(gly_chain.backprop(gly_param[0], gly_param[1], res_grad)))
        else:
            contrast.gly.append(np.zeros((2, n_grid, n_grid)))

    engine.set_param(less_sheet, 'rama_map_pot')
    for i in range(pos.shape[0]):
        engine.energy(pos[i])
        contrast.sheet[i] -= engine.get_output('rama_map_pot')[0, 0] * sheet_scale

    return Update(*[np.array(x) for x in contrast])


# ---------------------------------------------------------------------------
# Minibatch runner
# ---------------------------------------------------------------------------

def run_minibatch(worker_path, param, init_param_files, direc, minibatch,
                  solver, reg_scale, sim_time):
    os.makedirs(direc, exist_ok=True)
    print(direc)

    # Build nesterov look-ahead parameter files
    d_obj_files = dict(
        rot   = os.path.join(direc, 'nesterov_temp__sidechain.h5'),
        env   = os.path.join(direc, 'nesterov_temp__environment.h5'),
        hb    = os.path.join(direc, 'nesterov_temp__hbond'),
        sheet = os.path.join(direc, 'nesterov_temp__sheet'),
        rama  = os.path.join(direc, 'nesterov_temp__rama.dat'),
    )
    d_obj_param = param + solver.update_for_d_obj()
    expand_param(d_obj_param, init_param_files, d_obj_files)

    with open(os.path.join(direc, 'sim_time'), 'w') as f:
        print(sim_time, file=f)

    param_files_arg = cp.dumps(d_obj_files, protocol=0).decode('latin1')
    has_slurm = shutil.which('srun') is not None
    jobs = collections.OrderedDict()

    for nm, t in minibatch[::-1]:
        if has_slurm:
            args = [
                'srun', '--nodes=1', '--ntasks=1',
                '--cpus-per-task=%i' % n_threads,
                '--output=%s/%s.output_worker' % (direc, nm),
                sys.executable, worker_path,
                'worker', nm, direc, t.fasta, t.init_path,
                str(t.n_res), t.chi, param_files_arg, str(sim_time),
                str(d_obj_param.hb),
            ]
            output_dest = None
        else:
            args = [
                sys.executable, worker_path,
                'worker', nm, direc, t.fasta, t.init_path,
                str(t.n_res), t.chi, param_files_arg, str(sim_time),
                str(d_obj_param.hb),
            ]
            output_dest = open('%s/%s.output_worker' % (direc, nm), 'w')

        jobs[nm] = sp.Popen(args, close_fds=True,
                            stdout=output_dest, stderr=sp.STDOUT)

    rmsd, change = dict(), []
    for nm, j in jobs.items():
        if j.wait() != 0:
            print(nm, 'WORKER_FAIL')
            continue
        if not has_slurm and j.stdout:
            j.stdout.close()
        try:
            with open('%s/%s.divergence.pkl' % (direc, nm), 'rb') as f:
                div = cp.load(f)
                rmsd[nm]  = (div['rmsd_restrain'], div['rmsd'])
                change.append(div['contrast'])
        except Exception as e:
            print(nm, 'RESULT_READ_FAIL', e)

    # every worker has exited, so the 35 MB library they shared can go
    if os.path.exists(d_obj_files['rama']):
        os.remove(d_obj_files['rama'])

    if not change:
        raise RuntimeError('All jobs failed')

    d_param = backprop_deriv(
        param,
        Update(*[np.sum(x, axis=0) for x in zip(*change)]),
        reg_scale,
    )

    with open('%s/rmsd.pkl' % direc, 'wb') as f:
        cp.dump(rmsd, f, -1)
    print()
    print('Median RMSD %.2f %.2f' % tuple(np.median(np.array(list(rmsd.values())), axis=0)))

    new_files = dict(
        rot   = os.path.join(direc, 'sidechain.h5'),
        env   = os.path.join(direc, 'environment.h5'),
        hb    = os.path.join(direc, 'hbond'),
        sheet = os.path.join(direc, 'sheet'),
    )
    new_param = param + solver.update_step(d_param)
    expand_param(new_param, init_param_files, new_files)
    solver.log_state(direc)

    print()
    print_param(new_param)

    return new_param


# ---------------------------------------------------------------------------
# Worker entry point
# ---------------------------------------------------------------------------

def main_worker():
    assert is_worker
    tstart   = time.time()
    code     = sys.argv[2]
    direc    = sys.argv[3]
    fasta    = sys.argv[4]
    init_path = sys.argv[5]
    n_res    = int(sys.argv[6])
    chi      = sys.argv[7]

    param_files = cp.loads(sys.argv[8].encode('latin1'))
    sim_time    = float(sys.argv[9])
    # Current hbond scale.  The config is built from the already-scaled file, so the derivative
    # of the hbond potential with respect to this scale is E/hb_scale.
    hb_scale    = float(sys.argv[10])

    n_frame        = 250.
    frame_interval = int(sim_time / n_frame)
    input_dir      = os.path.dirname(init_path)

    # upside_config uses np.load() which cannot read pickle files in modern numpy.
    # Convert .pkl initial structures to a temporary .npy file in the work directory.
    if init_path.endswith('.pkl'):
        npy_init_path = '%s/%s.init.npy' % (direc, code)
        with open(init_path, 'rb') as f:
            init_coords = cp.load(f, encoding='latin1')
        np.save(npy_init_path, init_coords[:, :, 0])
    else:
        npy_init_path = init_path

    kwargs = dict(
        environment_potential      = param_files["env"],
        environment_potential_type = 0,
        rotamer_interaction   = param_files['rot'],
        rotamer_placement     = param_files['rot'],
        initial_structure     = npy_init_path,
        hbond_energy          = param_files['hb'],
        rama_sheet_mix_energy = param_files['sheet'],
        dynamic_rotamer_1body = True,
        rama_library          = param_files['rama'],
        rama_param_deriv      = True,   # writes more_/less_sheet_rama_pot_* and sheet_eps
        reference_state_rama  = os.path.join(input_dir, 'rama_reference.pkl'),
    )

    T = 0.80 * (1. + np.sqrt(100. / n_res) * 0.020 * np.arange(n_threads - 1)) ** 2
    T = np.concatenate((T[0:1], T))

    try:
        config_base = '%s/%s.base.h5' % (direc, code)
        ru.upside_config(fasta, config_base, **kwargs)
        configs = ['%s/%s.run.%i.h5' % (direc, code, i) for i in range(len(T))]
        for i in range(len(T)):
            shutil.copyfile(config_base, configs[i])
        ru.advanced_config(configs[0],
                           restraint_groups=['0-%i' % (n_res - 1)],
                           restraint_spring_constant=native_restraint_k)
    except tb.NoSuchNodeError:
        raise RuntimeError('CONFIG_FAIL')
    except Exception as e:
        print('Configuration failed:', e)
        raise RuntimeError('CONFIG_FAIL')

    j = ru.run_upside(
        '', configs, sim_time, frame_interval,
        n_threads=n_threads, temperature=T,
        swap_sets=ru.swap_table2d(1, len(T)),
        mc_interval=5., replica_interval=10.,
    )
    if j.job.wait() != 0:
        raise RuntimeError('RUN_FAIL')

    with tb.open_file(configs[0]) as t:
        target       = t.root.input.pos[:, :, 0]
        pos_restrain = t.root.output.pos[int(n_frame / 2):, 0]

    with tb.open_file(configs[1]) as t:
        pos_free = t.root.output.pos[int(n_frame / 2):, 0]

    # The current (S, A) are already in the library this worker was handed, so they are read
    # back from it rather than threaded through the command line as a 5,184-value array.
    gly_chain = rgg.GlycineMapChain(uc.read_fasta(open(fasta)), param_files['rama'],
                                    param_files['sheet'])
    gly_param = rgg.read_gly_maps(param_files['rama'])

    alldiv = compute_divergence(
        config_base,
        np.concatenate([pos_restrain, pos_free], axis=0),
        hb_scale, gly_chain, gly_param,
    )
    if alldiv is None:
        raise RuntimeError('DIVERGENCE_FAIL')

    n = len(pos_restrain)
    contrast = Update(*[x[:n].mean(axis=0) - x[n:].mean(axis=0) for x in alldiv])

    divergence = dict(
        contrast      = contrast,
        rmsd_restrain = ru.traj_rmsd(pos_restrain[:, rmsd_k:-rmsd_k],
                                     target[rmsd_k:-rmsd_k]).mean(),
        rmsd          = ru.traj_rmsd(pos_free[:, rmsd_k:-rmsd_k],
                                     target[rmsd_k:-rmsd_k]).mean(),
        walltime      = time.time() - tstart,
    )

    for fn in [config_base] + configs:
        os.remove(fn)
    if init_path.endswith('.pkl') and os.path.exists(npy_init_path):
        os.remove(npy_init_path)

    with open('%s/%s.divergence.pkl' % (direc, code), 'wb') as f:
        cp.dump(divergence, f, -1)


# ---------------------------------------------------------------------------
# Main loop (restart from checkpoint)
# ---------------------------------------------------------------------------

def main_loop(state_bytes, max_iter):
    for _ in range(max_iter):
        state = cp.loads(state_bytes)

        epoch = state['epoch']
        i_mb  = state['i_mb']
        print('#########################################')
        print('####      EPOCH %2i MINIBATCH %2i      ####' % (epoch, i_mb))
        print('#########################################\n')
        sys.stdout.flush()

        tstart = time.time()
        state['mb_direc'] = os.path.join(
            state['base_dir'],
            'epoch_%02i_minibatch_%02i' % (epoch, i_mb),
        )
        if os.path.exists(state['mb_direc']):
            shutil.rmtree(state['mb_direc'], ignore_errors=True)

        state['param'] = run_minibatch(
            state['worker_path'], state['param'],
            state['init_param_files'], state['mb_direc'],
            state['minibatches'][i_mb],
            state['solver'], float(state['n_prot']), state['sim_time'],
        )
        print('\n%.0f seconds elapsed this minibatch' % (time.time() - tstart))
        sys.stdout.flush()

        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1

        state_bytes = cp.dumps(state, -1)
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'), 'wb') as f:
            f.write(state_bytes)


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

def main_initialize(args):
    init_dir, protein_dir, protein_list, base_dir = args
    os.makedirs(base_dir, exist_ok=True)

    state = dict()
    state['init_dir']  = init_dir
    state['base_dir']  = base_dir

    # Copy the worker script into the output directory so each checkpoint
    # carries the exact code that produced it.
    state['worker_path'] = os.path.join(base_dir, 'ConDiv.py')
    shutil.copy(__file__, state['worker_path'])

    # Build training set
    if protein_list != 'cached':
        print('Reading training set')
        with open(protein_list) as f:
            protein_names = [x.split()[0] for x in f]
        assert protein_names[0] == 'prot'
        protein_names = protein_names[1:]

        training_set, excluded = dict(), []
        for code in sorted(protein_names):
            base = os.path.join(protein_dir, code)
            with open(base + '.initial.pkl', 'rb') as f:
                native_pos = cp.load(f, encoding='latin1')[:, :, 0]
            n_res = len(native_pos) // 3
            max_sep = np.sqrt(np.sum(np.diff(native_pos, axis=0) ** 2, axis=-1)).max()
            if max_sep < 2.:
                training_set[code] = Target(
                    base + '.fasta', native_pos, base + '.initial.pkl',
                    base + '.initial.pkl', n_res, base + '.chi',
                )
            else:
                excluded.append(code)
                print(code)

        print('Excluded %i proteins due to chain breaks' % len(excluded))
        with open(os.path.join(base_dir, 'cd_training.pkl'), 'wb') as f:
            cp.dump(training_set, f, -1)
    else:
        with open(os.path.join(base_dir, 'cd_training.pkl'), 'rb') as f:
            training_set = cp.load(f, encoding='latin1')

    training_list = sorted(training_set.items(), key=lambda x: (x[1].n_res, x[0]))
    np.random.shuffle(training_list)

    excess = len(training_list) % minibatch_size
    if excess:
        training_list = training_list[:-excess]
    n_mb = len(training_list) // minibatch_size
    state['minibatches'] = [training_list[i::n_mb] for i in range(n_mb)]
    state['n_prot']      = n_mb * minibatch_size
    print('Constructed %i minibatches of size %i (%i proteins)' %
          (n_mb, minibatch_size, state['n_prot']))

    # Load initial parameters
    if init_dir != 'cached':
        print('Loading initial parameters...')
        state['param'], state['init_param_files'] = get_init_param(
            init_dir, os.path.join(protein_dir, 'rama.dat'))
        print('Initial parameters loaded.')
        with open(os.path.join(base_dir, 'condiv_init.pkl'), 'wb') as f:
            cp.dump((init_dir, state['param'], state['init_param_files']), f, -1)
    else:
        with open(os.path.join(base_dir, 'condiv_init.pkl'), 'rb') as f:
            init_dir, state['param'], state['init_param_files'] = cp.load(f, encoding='latin1')

    # Optimizer: rot, env, hb and sheet, the same four the Theano original trained
    state['initial_alpha'] = Update(
        env   = 0.1,
        cov   = 0.,    # included in rot latent vector, not a separate variable
        rot   = 0.5,
        hyd   = 0.,    # included in rot latent vector
        # The original trained hb at 0.02 on an absolute per-hbond energy of magnitude ~2.1.
        # Here hb is a multiplicative scale starting at 1.0, so dividing by that reference energy
        # makes one step move the hbond energies by the same absolute amount as the original.
        hb    = 0.02 / 1.96,
        sheet = 0.03,  # the original's value, on the same common-offset parameter
        # The map is in nats with basins a few units deep, so a step of ~0.005 nats lets 500
        # minibatches move it by up to ~2.5 where the gradient is consistent: the right order for
        # a handedness that runs from 0 to -1.24, and small enough that a wrong direction shows up
        # long before it does damage.
        gly   = 0.02,
    ) * 0.25

    state['solver']   = rp.AdamSolver(len(state['initial_alpha']),
                                      alpha=state['initial_alpha'])
    state['sim_time'] = 1000. * 4
    state['epoch']    = 0
    state['i_mb']     = 0

    print('\nOptimizing with solver', state['solver'])
    print()
    print_param(state['param'])

    return state


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'worker':
        main_worker()

    elif len(sys.argv) > 1 and sys.argv[1] == 'restart':
        assert len(sys.argv) == 4, 'Usage: ConDiv.py restart <checkpoint.pkl> <n_steps>'
        print('Running as PID %i on host %s' % (os.getpid(), socket.gethostname()))
        with open(sys.argv[2], 'rb') as f:
            checkpoint = f.read()
        main_loop(checkpoint, int(sys.argv[3]))

    elif len(sys.argv) > 1 and sys.argv[1] == 'initialize':
        assert len(sys.argv) == 6, \
            'Usage: ConDiv.py initialize <init_param_dir> <upside_input_dir> <pdb_list> <output_dir>'
        initial_state = main_initialize(sys.argv[2:])
        with open(os.path.join(initial_state['base_dir'], 'initial_checkpoint.pkl'), 'wb') as f:
            cp.dump(initial_state, f, -1)

    else:
        print('Usage:')
        print('  ConDiv.py initialize <init_param_dir> <protein_dir> <pdb_list> <output_dir>')
        print('  ConDiv.py restart <checkpoint.pkl> <n_steps>')
        raise SystemExit(1)
