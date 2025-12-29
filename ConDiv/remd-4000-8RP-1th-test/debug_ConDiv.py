#!/usr/bin/env python3
import sys, os
import numpy as np
import subprocess as sp
import tables as tb
import pickle as cp
from glob import glob
import re
import shutil
import collections
import time
import socket
import torch

# Ensure local modules are found
rp_path = os.path.abspath(os.path.join(os.getcwd(), "../..", "src")) 
if rp_path not in sys.path:
    sys.path.append(rp_path)

is_worker = __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'worker'

if not is_worker:
    import rotamer_parameter_estimation as rp

import run_upside as ru
import upside_engine as ue

np.set_printoptions(precision=2, suppress=True)

## --- DEBUG SETTINGS ---
n_threads = 1        # Use 1 thread to avoid Mac OpenMP crashes
minibatch_size = 1   # Process 1 protein at a time
debug_sim_time = 500.0
## ----------------------

native_restraint_strength = 1./3.**2
rmsd_k = 15

resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL']

hydrophobicity_order = ['ASP', 'GLU', 'LYS', 'HIS', 'ARG', 'GLY', 'ASN', 'GLN',
                        'ALA', 'SER', 'THR', 'PRO', 'CYS', 'VAL', 'MET', 'TYR',
                        'ILE', 'LEU', 'PHE', 'TRP']

Target = collections.namedtuple('Target', 'fasta native native_path init_path n_res chi'.split())
UpdateBase = collections.namedtuple('UpdateBase', 'env cov rot hyd hb sheet'.split())

class Update(UpdateBase):
    def __init__(self, *args, **kwargs):
        pass

    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False

        if is_seq:
            ret = []
            assert len(self) == len(other)
            for a, b in zip(self, other):
                if a is None or b is None:
                    ret.append(None)
                else:
                    ret.append(op(a, b))
        else:
            ret = [None if a is None else op(a, other) for a in self]
        return Update(*ret)

    def __add__(self, other): return self._do_binary(other, lambda a, b: a + b)
    def __sub__(self, other): return self._do_binary(other, lambda a, b: a - b)
    def __mul__(self, other): return self._do_binary(other, lambda a, b: a * b)
    def __truediv__(self, other): return self._do_binary(other, lambda a, b: a / b)

def print_param(param):
    print(f'hb    {float(param.hb):.3f}')
    if np.ndim(param.sheet) == 0:
        print(f'sheet {float(param.sheet):.3f}')
    else:
        print(f'sheet (mean) {np.mean(param.sheet):.3f}')

    print('env')
    env_data = param.env
    if hasattr(env_data, 'detach'): 
        env_data = env_data.detach().numpy()
        
    try:
        env_dict = dict(zip(resnames, env_data[:, 1::2]))
        for r in hydrophobicity_order:
            print('   ', r, env_dict[r])
    except IndexError:
        print('    (Environment parameters not fully initialized yet)')

def get_d_obj_torch():
    def student_t_neglog(t, scale, nu):
        return 0.5 * (nu + 1.) * torch.log(1. + (1. / (nu * scale**2)) * t**2)

    def cutoff_func(x, scale):
        return 1. / (1. + torch.exp(x / scale))

    def lower_bound(x, lb):
        return torch.where(x < lb, 1e0 * (x - lb)**2, torch.zeros_like(x))

    r_for_cov_energy = torch.arange(1, rp.n_knot_hb - 1) * rp.hb_dr
    r_for_rot_energy = torch.arange(1, rp.n_knot_sc - 1) * rp.sc_dr
    r_for_hyd_energy = torch.arange(1, rp.n_knot_hb - 1) * rp.hb_dr
    
    cov_expect = 5. * cutoff_func(r_for_cov_energy - 2., 0.2)
    rot_expect = 5. * cutoff_func(r_for_rot_energy - 2., 0.2)
    hyd_expect = 5. * cutoff_func(r_for_hyd_energy - 1., 0.2)

    def compute_grad(lparam_np, rot_coupling, cov_coupling, hyd_coupling, reg_scale):
        lparam = torch.tensor(lparam_np, dtype=torch.float64, requires_grad=True)
        rot_c = torch.tensor(rot_coupling, dtype=torch.float64)
        cov_c = torch.tensor(cov_coupling, dtype=torch.float64)
        hyd_c = torch.tensor(hyd_coupling, dtype=torch.float64)
        
        rot_expr, cov_expr, hyd_expr, _, _, _ = rp.unpack_params_pytorch(lparam)
        
        hb_n_knot = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_hb, rp.n_knot_hb)
        sc_n_knot = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_sc, rp.n_knot_sc)
        
        cov_energy = rp.quadspline_energy_torch(cov_expr, hb_n_knot)
        rot_energy = rp.quadspline_energy_torch(rot_expr, sc_n_knot)
        hyd_energy = rp.quadspline_energy_torch(hyd_expr, hb_n_knot)

        cov_reg = student_t_neglog(cov_energy - cov_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale
        rot_reg = student_t_neglog(rot_energy - rot_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale
        hyd_reg = student_t_neglog(hyd_energy - hyd_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale

        reg_expr = (cov_reg + hyd_reg + rot_reg +
                    lower_bound(cov_energy, -6.).sum() +
                    lower_bound(rot_energy, -6.).sum() +
                    lower_bound(hyd_energy, -6.).sum())

        obj_expr = ((rot_c * rot_expr).sum() +
                    (cov_c * cov_expr).sum() +
                    (hyd_c * hyd_expr).sum() +
                    reg_expr)

        obj_expr.backward()
        return lparam.grad.detach().numpy()

    return compute_grad

if not is_worker:
    d_obj = get_d_obj_torch()

    def get_init_param(init_dir):
        init_param_files = dict(
                env = os.path.join(init_dir, 'environment.h5'),
                rot = os.path.join(init_dir, 'sidechain.h5'),
                hb  = os.path.join(init_dir, 'hbond.h5'),
                sheet  = os.path.join(init_dir, 'sheet'))

        with tb.open_file(init_param_files['rot']) as t:
            rotp = t.root.pair_interaction[:]
            covp = t.root.coverage_interaction[:]
            hydp = t.root.hydrophobe_interaction[:]
            hydplp = t.root.hydrophobe_placement[:]
            rotposp = t.root.rotamer_center_fixed[:]
            rotscalarp = np.zeros((rotposp.shape[0],))

        with tb.open_file(init_param_files['env']) as t:
            env = t.root.energies[:,:-1]

        with tb.open_file(init_param_files['hb']) as t:
            hb_val = t.root.parameter[0] 

        sheet = np.loadtxt(init_param_files['sheet'])

        param = Update(*([None]*6))._replace(
                env = env,
                rot = rp.pack_param(rotp, covp, hydp, hydplp, rotposp, rotscalarp[:,None]),
                hb  = hb_val,
                sheet = sheet)
        return param, init_param_files

    def expand_param(params, orig_param_files, new_param_files):
        rotp, covp, hydp, hydplp, rotposp, rotscalarp = rp.unpack_params(params.rot)

        shutil.copyfile(orig_param_files['rot'], new_param_files['rot'])
        with tb.open_file(new_param_files['rot'], 'a') as t:
            t.root.pair_interaction[:]       = rotp
            t.root.coverage_interaction[:]   = covp
            t.root.hydrophobe_interaction[:] = hydp
            t.root.hydrophobe_placement[:]   = hydplp
            t.root.rotamer_center_fixed[:]   = rotposp

        shutil.copyfile(orig_param_files['hb'], new_param_files['hb'])
        with tb.open_file(new_param_files['hb'], 'a') as t:
            p = t.root.parameter[:]
            p[0] = params.hb
            t.root.parameter[:] = p
        
        if np.ndim(params.sheet) == 0:
            with open(new_param_files['sheet'], 'w') as f: print(params.sheet, file=f)
        else:
            np.savetxt(new_param_files['sheet'], params.sheet)

        shutil.copyfile(orig_param_files['env'], new_param_files['env'])
        with tb.open_file(new_param_files['env'], 'a') as t:
            tmp = np.zeros(t.root.energies.shape)
            tmp[:,:-1] = params.env
            tmp[:,-1]  = tmp[:,-3]
            t.root.energies[:] = tmp

    def backprop_deriv(param, deriv_update, reg_scale):
        envd = deriv_update.env[:,:-1].copy()
        envd[:,-2] += deriv_update.env[:,-1]
        
        rot_grad = d_obj(param.rot, deriv_update.rot, deriv_update.cov, deriv_update.hyd, reg_scale)
        
        return deriv_update._replace(
                rot = rot_grad,
                cov = 0.,
                hyd = 0.,
                env = envd)

def run_minibatch(worker_path, param, initial_param_files, direc, minibatch, solver, reg_scale, sim_time):
    if not os.path.exists(direc): os.mkdir(direc)
    print(f"Running minibatch in {direc}")

    d_obj_param_files = dict([(k,os.path.join(direc, 'nesterov_temp__'+os.path.basename(x))) 
        for k,x in initial_param_files.items()])
    expand_param(param+solver.update_for_d_obj(), initial_param_files, d_obj_param_files)

    with open(os.path.join(direc,'sim_time'),'w') as f:
            print(sim_time, file=f)

    jobs = collections.OrderedDict()
    for nm,t in minibatch[::-1]:
        args = [sys.executable, worker_path,
                'worker', nm, direc, t.fasta, t.init_path, str(t.n_res), t.chi,
                cp.dumps(d_obj_param_files, protocol=0).decode('latin1'),
                str(sim_time)]
        
        pickled_files_hex = cp.dumps(d_obj_param_files).hex()
        args[9] = pickled_files_hex

        outfile = open(f'{direc}/{nm}.output_worker', 'w')
        print(f"Launching worker for {nm}...")
        jobs[nm] = sp.Popen(args, close_fds=True, stdout=outfile, stderr=sp.STDOUT)

    rmsd = dict()
    change = []
    
    for nm,j in jobs.items():
        retcode = j.wait()
        if retcode != 0:
            print(f"WARNING: Worker for {nm} exited with code {retcode} (ignoring segfault). Checking for output...")
        
        try:
            with open(f'{direc}/{nm}.divergence.pkl', 'rb') as f:
                divergence = cp.load(f)
                rmsd[nm] = (divergence['rmsd_restrain'], divergence['rmsd'])
                change.append(divergence['contrast'])
        except FileNotFoundError:
            print(nm, 'NO_OUTPUT (Real Failure)')
            continue
            
    if not change:
        raise RuntimeError('All jobs failed (No output generated)')

    print("Workers finished. Computing gradients...")
    d_param = backprop_deriv(param, Update(*[np.sum(x,axis=0) for x in zip(*change)]), reg_scale)

    with open(f'{direc}/rmsd.pkl', 'wb') as f:
        cp.dump(rmsd, f, -1)
    print()
    print('Median RMSD %.2f %.2f' % tuple(np.median(np.array(list(rmsd.values())), axis=0)))

    new_param_files = dict([(k,os.path.join(direc, os.path.basename(x))) for k,x in initial_param_files.items()])
    new_param = param+solver.update_step(d_param)
    expand_param(new_param, initial_param_files, new_param_files)
    solver.log_state(direc)

    print()
    print_param(new_param)

    return new_param


def compute_divergence(config_base, pos):
    try:
        with tb.open_file(config_base) as t:
            seq        = t.root.input.sequence[:]
            
            rama_node = t.root.input.potential.rama_map_pot
            sheet_restypes = rama_node._v_attrs.restype
            if isinstance(sheet_restypes[0], bytes):
                sheet_restypes = [x.decode('utf-8') for x in sheet_restypes]
                
            eps = rama_node._v_attrs.sheet_eps
            sheet_scale = 1./(2.*eps)
            
            more_sheets = {}
            less_sheets = {}
            for res in sheet_restypes:
                try:
                    m = getattr(rama_node, f'more_sheet_rama_pot_{res}')[:]
                    l = getattr(rama_node, f'less_sheet_rama_pot_{res}')[:]
                    more_sheets[res] = m
                    less_sheets[res] = l
                except tb.NoSuchNodeError:
                    more_sheets[res] = None
                    less_sheets[res] = None

            hb_strength = t.root.input.potential.hbond_energy.parameters[0]

            rot_param_shape = t.root.input.potential.rotamer.pair_interaction.interaction_param.shape
            cov_param_shape = t.root.input.potential.hbond_coverage.interaction_param.shape
            hyd_param_shape = t.root.input.potential.hbond_coverage_hydrophobe.interaction_param.shape
            
            # Check for env group
            if not hasattr(t.root.input.potential, 'nonlinear_coupling_environment'):
                print(f"ANALYSIS_FAIL: Missing nonlinear_coupling_environment in {config_base}")
                return None
            
            env_param_shape = t.root.input.potential.nonlinear_coupling_environment.coeff.shape
            env_weights_shape = t.root.input.potential.nonlinear_coupling_environment.weights.shape
    except Exception as e:
        print(os.path.basename(config_base)[:5],'ANALYSIS_FAIL',e)
        return None

    engine = ue.Upside(config_base)
    contrast = Update([],[],[],[],[],[])

    for i in range(pos.shape[0]):
        engine.energy(pos[i]) 
        
        contrast.cov.append(engine.get_param_deriv(cov_param_shape, 'hbond_coverage'))
        contrast.rot.append(engine.get_param_deriv(rot_param_shape, 'rotamer'))
        contrast.hyd.append(engine.get_param_deriv(hyd_param_shape, 'hbond_coverage_hydrophobe'))
        contrast.hb .append(engine.get_output('hbond_energy')[0,0]/hb_strength)
        
        # Handle env derivatives including weights (which we ignore)
        total_env_params = np.prod(env_param_shape) + np.prod(env_weights_shape)
        env_full = engine.get_param_deriv((int(total_env_params),), 'nonlinear_coupling_environment')
        contrast.env.append(env_full[:np.prod(env_param_shape)].reshape(env_param_shape))
        
        current_sheet_grad = np.zeros(len(sheet_restypes))
        
        for idx, res in enumerate(sheet_restypes):
            m = more_sheets[res]
            l = less_sheets[res]
            if m is None: continue 
                
            engine.set_param(m, 'rama_map_pot')
            engine.energy(pos[i])
            e_plus = engine.get_output('rama_map_pot')[0,0]
            
            engine.set_param(l, 'rama_map_pot')
            engine.energy(pos[i])
            e_minus = engine.get_output('rama_map_pot')[0,0]
            
            current_sheet_grad[idx] = (e_plus - e_minus) * sheet_scale
            
        contrast.sheet.append(current_sheet_grad)

    contrast = Update(*[np.array(x) for x in contrast])
    return contrast


def main_worker():
    tstart = time.time()
    assert is_worker
    code        = sys.argv[2]
    direc       = sys.argv[3]
    fasta       = sys.argv[4]
    init_path   = sys.argv[5]
    n_res       = int(sys.argv[6])
    chi         = sys.argv[7]
    param_files = cp.loads(bytes.fromhex(sys.argv[8]))
    sim_time    = float(sys.argv[9])

    n_frame = 250.
    frame_interval = int(sim_time / n_frame)
    if frame_interval < 1: frame_interval = 1

    hb = param_files['hb']
    sheet_mix = param_files['sheet']
    project_root = os.path.abspath(os.path.join(os.getcwd(), "../.."))

    kwargs = dict(
            environment_potential = param_files['env'],
            rotamer_interaction   = param_files['rot'],
            rotamer_placement     = param_files['rot'],
            initial_structure     = init_path,
            hbond_energy          = hb,
            dynamic_rotamer_1body = True,
            rama_sheet_mix_energy = sheet_mix,
            rama_param_deriv      = True,
            environment_potential_type = "0",  
            rama_library          = os.path.join(project_root, "parameters/common/rama.dat"),
            reference_state_rama  = os.path.join(project_root, "parameters/common/rama_reference.pkl"),
    )
    
    if n_threads < 2:
        T = np.array([0.80, 0.80])
    else:
        T = 0.80 * (1. + np.sqrt(100./n_res)*0.020*np.arange(n_threads-1))**2
        T = np.concatenate((T[0:1],T))

    try:
        config_base = '%s/%s.base.h5' % (direc,code)
        ru.upside_config(fasta, config_base, **kwargs)

        # --- FIX: Patch metadata by reading from source file ---
        try:
            # 1. Read inv_dx from source
            src_env_path = param_files['env']
            with tb.open_file(src_env_path, 'r') as t_src:
                # Assuming source structure is typically root.energies
                if hasattr(t_src.root.energies._v_attrs, 'inv_dx'):
                    src_inv_dx = t_src.root.energies._v_attrs.inv_dx
                else:
                    # Fallback calculation: 18 knots (17 intervals) over 12A range
                    src_inv_dx = 17.0 / 12.0 # ~1.4167
            
            # 2. Write inv_dx to config_base
            with tb.open_file(config_base, 'r+') as t_dst:
                env_node = t_dst.root.input.potential.nonlinear_coupling_environment.coeff
                env_node._v_attrs.inv_dx = src_inv_dx
                print(f"Patched {code} env metadata: inv_dx set to {src_inv_dx}")
                
        except Exception as e:
            print(f"Warning: Failed to patch env metadata for {code}: {e}")
        # -----------------------------------------------------

        configs = [re.sub(r'\.base\.h5','.run.%i.h5'%i_rs, config_base) for i_rs in range(len(T))]
        for i in range(1,len(T)):
            shutil.copyfile(config_base, configs[i])

        upside_config_script = os.path.join(os.path.dirname(ru.__file__), 'upside_config.py')
        
        cmd = [sys.executable, upside_config_script]
        key_map = {
            'environment_potential': '--environment-potential',
            'rotamer_interaction': '--rotamer-interaction',
            'rotamer_placement': '--rotamer-placement',
            'initial_structure': '--initial-structure',
            'hbond_energy': '--hbond-energy',
            'rama_sheet_mix_energy': '--rama-sheet-mixing-energy',
            'rama_library': '--rama-library',
            'reference_state_rama': '--reference-state-rama',
            'dynamic_rotamer_1body': '--dynamic-rotamer-1body',
            'rama_param_deriv': '--rama-param-deriv',
            'environment_potential_type': '--environment-potential-type'
        }

        cmd.append(f'--fasta={fasta}')
        cmd.append(f'--output={configs[0]}')
        
        for k, v in kwargs.items():
            if k in key_map:
                if v is True:
                    cmd.append(key_map[k])
                elif v is not None:
                    cmd.append(f'{key_map[k]}={v}')

        cmd.append(f'--restraint-group=0-{n_res-1}')
        cmd.append(f'--restraint-spring={native_restraint_strength}')

        try:
            sp.check_output(cmd, stderr=sp.STDOUT)
        except sp.CalledProcessError as e:
            print("Error generating restrained config:")
            print(e.output.decode('ascii'))
            raise RuntimeError('CONFIG_FAIL_RESTRAINT')
    except tb.NoSuchNodeError:
        raise RuntimeError('CONFIG_FAIL')
    
    j = ru.run_upside('', configs, sim_time, frame_interval, n_threads=n_threads, temperature=T,
            swap_sets=ru.swap_table2d(1,len(T)), mc_interval=5., replica_interval=10.)
    
    j.job.wait() # Ignore segfault

    swap_stats = []

    with tb.open_file(configs[0]) as t:
        target = t.root.input.pos[:,:,0]
        o = t.root.output
        pos_restrain = o.pos[int(n_frame/2):,0]

    with tb.open_file(configs[1]) as t:
        o = t.root.output
        pos_free     = o.pos[int(n_frame/2):,0]
        swap_stats.extend(o.replica_cumulative_swaps[-1] - o.replica_cumulative_swaps[int(n_frame/2)])

    for nrep in range(3,len(T),2):
        with tb.open_file(configs[nrep]) as t:
            o = t.root.output
            swap_stats.extend(o.replica_cumulative_swaps[-1] - o.replica_cumulative_swaps[int(n_frame/2)])

    divergence = dict()
    alldiv = compute_divergence(config_base, np.concatenate([pos_restrain,pos_free],axis=0))
    divergence['contrast'] = Update(*[x[:len(pos_restrain)].mean(axis=0) -
                                      x[len(pos_restrain):].mean(axis=0) for x in alldiv])

    divergence['rmsd_restrain']= ru.traj_rmsd(pos_restrain[:,rmsd_k:-rmsd_k], target[rmsd_k:-rmsd_k]).mean()
    divergence['rmsd']         = ru.traj_rmsd(pos_free    [:,rmsd_k:-rmsd_k], target[rmsd_k:-rmsd_k]).mean()
    divergence['walltime'] = time.time()-tstart
    divergence['swap_stats'] = np.array(swap_stats)

    for fn in [config_base] + configs:
        if os.path.exists(fn): os.remove(fn)

    with open('%s/%s.divergence.pkl' % (direc,code),'wb') as f:
        cp.dump(divergence, f, -1)


def main_loop(state_str, max_iter):
    def main_loop_iteration(state):
        print('#########################################')
        print('####      EPOCH %2i MINIBATCH %2i      ####' % (state['epoch'],state['i_mb']))
        print('#########################################')
        sys.stdout.flush()
    
        tstart = time.time()
        state['mb_direc'] = os.path.join(state['base_dir'],'epoch_%02i_minibatch_%02i'%(state['epoch'],state['i_mb']))
        if os.path.exists(state['mb_direc']):
            shutil.rmtree(state['mb_direc'], ignore_errors=True)
        
        state['param'] = run_minibatch(state['worker_path'], state['param'], state['init_param_files'],
                state['mb_direc'], state['minibatches'][state['i_mb']],
                state['solver'], state['n_prot']*1., state['sim_time'])
        print()
        print('%.0f seconds elapsed this minibatch'%(time.time()-tstart))
        sys.stdout.flush()

        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1

    for i in range(max_iter):
        if hasattr(state_str, 'encode'): 
             if os.path.exists(state_str):
                 with open(state_str, 'rb') as f:
                     state = cp.load(f)
             else:
                 raise FileNotFoundError(f"Checkpoint file not found: {state_str}\n")
        else:
             state = cp.loads(state_str)

        state['sim_time'] = debug_sim_time 

        main_loop_iteration(state)
        
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'),'wb') as f:
            cp.dump(state, f, -1)
        
        state_str = cp.dumps(state, -1)


def main_initialize(args):
    state = dict()
    state['init_dir'], state['protein_dir'], protein_list, state['base_dir'], = args
    if not os.path.exists(state['base_dir']): os.mkdir(state['base_dir'])

    state['worker_path'] = os.path.join(state['base_dir'],'ConDiv.py')
    shutil.copy(__file__, state['worker_path'])

    if protein_list != 'cached':
        print('Reading training set')
        with open(protein_list) as f:
            protein_names = [x.split()[0] for x in f]
            assert protein_names[0]=='prot'
            protein_names = protein_names[1:]
                
        training_set = dict()
        excluded_prot = []
        for code in sorted(protein_names):
            base = os.path.join(state['protein_dir'], code)

            with open(base+'.initial.pkl', 'rb') as f:
                native_pos = cp.load(f, encoding='latin1')[:,:,0] 
                n_res = len(native_pos)//3

            max_sep = np.sqrt(np.sum(np.diff(native_pos,axis=0)**2,axis=-1)).max()
            if max_sep < 2.:
                training_set[code] = Target(base+'.fasta', native_pos, base+'.initial.pkl',
                    base+'.initial.pkl', n_res, base+'.chi')
            else:
                excluded_prot.append(code)
                print(code)

        print('Excluded %i proteins due to chain breaks'%len(excluded_prot))
        with open(os.path.join(state['base_dir'], 'cd_training.pkl'),'wb') as f: cp.dump(training_set, f, -1)
    else:
        with open(os.path.join(state['base_dir'], 'cd_training.pkl'), 'rb') as f:
            training_set = cp.load(f)

    training_list = sorted(list(training_set.items()), key=lambda x: (x[1].n_res,x[0]))
    np.random.shuffle(training_list)
    
    # DEBUG: Reduce to 1 protein
    print("DEBUG MODE: Reducing training set to 1 protein.")
    training_list = training_list[:1]

    state['minibatches'] = [training_list] 
    state['n_prot'] = len(training_list)
    
    print('Constructed %i minibatches of size %i (%i proteins)' % (1, 1, state['n_prot']))


    if state['init_dir'] != 'cached':
        state['param'], state['init_param_files'] = get_init_param(state['init_dir'])
        with open(os.path.join(state['base_dir'], 'condiv_init.pkl'),'wb') as f:
            cp.dump((state['init_dir'],state['param'],state['init_param_files']),f,-1)
    else:
        with open(os.path.join(state['base_dir'], 'condiv_init.pkl'), 'rb') as f:
            state['init_dir'],state['param'],state['init_param_files'] = cp.load(f)

    state['initial_alpha'] = Update(*[0.1, 0., 0.5, 0., 0.02, 0.03]) * 0.25
    state['solver'] = rp.AdamSolver(len(state['initial_alpha']), alpha=state['initial_alpha']) 
    
    state['sim_time'] = debug_sim_time

    print('Optimizing with solver', state['solver'])
    print_param(state['param'])

    state['epoch'] = 0
    state['i_mb'] = 0
    return state

if __name__ == '__main__':
    if sys.argv[1] == 'worker':
        main_worker()

    elif sys.argv[1] == 'restart':
        print('Running as PID %i on host %s' % (os.getpid(), socket.gethostname()))
        main_loop(sys.argv[2], int(sys.argv[3]))

    elif sys.argv[1] == 'initialize':
        initial_state = main_initialize(sys.argv[2:])
        with open(os.path.join(initial_state['base_dir'],'initial_checkpoint.pkl'),'wb') as f:
            cp.dump(initial_state, f, -1)

    else:
        raise RuntimeError('Illegal mode %s'%sys.argv[1])