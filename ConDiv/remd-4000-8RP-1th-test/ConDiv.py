#!/usr/bin/env python3
import sys, os
os.environ['OMP_NUM_THREADS'] = '1'  # Prevent threading conflicts

import numpy as np
import subprocess as sp
import tables as tb
import pickle as cp  # Replaces cPickle
from glob import glob
import re
import shutil
import collections
import time
import socket
import torch  # Replaces Theano

class ExplosionError(Exception):
    pass

# Ensure local modules are found
# Ensure local modules are found
rp_path = os.path.abspath(os.path.join(os.getcwd(), "../..", "src")) 
if rp_path not in sys.path:
    sys.path.append(rp_path)

# Mock check for worker to avoid imports if running strictly as worker
# (Adjust based on how your modules are structured in Py3)
is_worker = __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'worker'

if not is_worker:
    import rotamer_parameter_estimation as rp

import run_upside as ru
import upside_engine as ue

np.set_printoptions(precision=2, suppress=True)

## Important parameters
n_threads = 12
native_restraint_strength = 1./3.**2
rmsd_k = 15
minibatch_size = 24
max_parallel_jobs = 12

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
        # Handle namedtuple init in Py3 cleanly
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
    def __truediv__(self, other): return self._do_binary(other, lambda a, b: a / b) # Py3 div

def print_param(param):
    # Handle hb (might be numpy scalar or float)
    print(f'hb    {float(param.hb):.3f}')
    
    # Handle sheet (Scalar vs Vector)
    if np.ndim(param.sheet) == 0:
        print(f'sheet {float(param.sheet):.3f}')
    else:
        print(f'sheet (mean) {np.mean(param.sheet):.3f}')

    print('env')
    
    # Check if param.env is a tensor or numpy array
    env_data = param.env
    if hasattr(env_data, 'detach'): 
        env_data = env_data.detach().numpy()
        
    # Ensure env_data has expected dimensions before zipping
    # (Just a safeguard against shape mismatches during init)
    try:
        env_dict = dict(zip(resnames, env_data[:, 1::2]))
        for r in hydrophobicity_order:
            print('   ', r, env_dict[r])
    except IndexError:
        print('    (Environment parameters not fully initialized yet)')

def get_d_obj_torch():
    """
    Replaces get_d_obj (Theano) with PyTorch logic.
    Assumes rp modules now provide PyTorch Tensors or compatible numpy arrays.
    """
    
    # NOTE: You must update 'rp' to expose these as torch.Tensors or conversion will be needed here.
    # For now, we wrap them assuming they are numpy arrays.
    def to_t(x): return torch.tensor(x, dtype=torch.float64, requires_grad=False)
    
    # We assume rp.unpack_*_expr were Theano variables. 
    # In PyTorch, we need the logic that generated those expressions.
    # Since we don't have rp source, we assume rp.unpack_params can be used 
    # to reconstruct the inputs for the energy function.

    def student_t_neglog(t, scale, nu):
        return 0.5 * (nu + 1.) * torch.log(1. + (1. / (nu * scale**2)) * t**2)

    def cutoff_func(x, scale):
        return 1. / (1. + torch.exp(x / scale))

    def lower_bound(x, lb):
        return torch.where(x < lb, 1e0 * (x - lb)**2, torch.zeros_like(x))

    # Grid setup (Fixed constants)
    r_for_cov_energy = torch.arange(1, rp.n_knot_hb - 1) * rp.hb_dr
    r_for_rot_energy = torch.arange(1, rp.n_knot_sc - 1) * rp.sc_dr
    r_for_hyd_energy = torch.arange(1, rp.n_knot_hb - 1) * rp.hb_dr
    
    cov_expect = 5. * cutoff_func(r_for_cov_energy - 2., 0.2)
    rot_expect = 5. * cutoff_func(r_for_rot_energy - 2., 0.2)
    hyd_expect = 5. * cutoff_func(r_for_hyd_energy - 1., 0.2)

    def compute_grad(lparam_np, rot_coupling, cov_coupling, hyd_coupling, reg_scale):
        # Convert inputs to torch tensors
        lparam = torch.tensor(lparam_np, dtype=torch.float64, requires_grad=True)
        rot_c = torch.tensor(rot_coupling, dtype=torch.float64)
        cov_c = torch.tensor(cov_coupling, dtype=torch.float64)
        hyd_c = torch.tensor(hyd_coupling, dtype=torch.float64)
        
        # Unpack parameters using rp logic (needs rp to support torch or we manually slice)
        # Note: This relies on rp.unpack_params returning Tensors if we pass it a Tensor
        # If rp is purely numpy, this section requires the splines to be rewritten in Torch.
        # Placeholder: We assume rp.quadspline_energy works with Torch tensors.
        
        # 1. Unpack params from flat lparam
        # This part depends heavily on rp.unpack_params implementation.
        # Assuming it splits lparam into the specific interaction tensors.
        #unpacked = rp.unpack_params_pytorch(lparam) # You must implement this in rp
        #cov_expr, rot_expr, hyd_expr = unpacked['cov'], unpacked['rot'], unpacked['hyd']

        rot_expr, cov_expr, hyd_expr, _, _, _ = rp.unpack_params_pytorch(lparam)
        
        # 2. Calculate Energies
        # You need a torch implementation of quadspline_energy in rp
        hb_n_knot = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_hb, rp.n_knot_hb)
        sc_n_knot = (rp.n_knot_angular, rp.n_knot_angular, rp.n_knot_sc, rp.n_knot_sc)
        
        cov_energy = rp.quadspline_energy_torch(cov_expr, hb_n_knot)
        rot_energy = rp.quadspline_energy_torch(rot_expr, sc_n_knot)
        hyd_energy = rp.quadspline_energy_torch(hyd_expr, hb_n_knot)

        # 3. Regularization
        cov_reg = student_t_neglog(cov_energy - cov_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale
        rot_reg = student_t_neglog(rot_energy - rot_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale
        hyd_reg = student_t_neglog(hyd_energy - hyd_expect[None, None, :, None, None], 200., 3.).sum() / reg_scale

        reg_expr = (cov_reg + hyd_reg + rot_reg +
                    lower_bound(cov_energy, -6.).sum() +
                    lower_bound(rot_energy, -6.).sum() +
                    lower_bound(hyd_energy, -6.).sum())

        # 4. Objective
        obj_expr = ((rot_c * rot_expr).sum() +
                    (cov_c * cov_expr).sum() +
                    (hyd_c * hyd_expr).sum() +
                    reg_expr)

        # 5. Backward
        obj_expr.backward()
        return lparam.grad.detach().numpy()

    return compute_grad


if not is_worker:
    # d_obj = get_d_obj() # OLD THEANO
    # d_obj = get_d_obj_torch() # NEW TORCH (Requires rp updates)
    # For now, we mock d_obj so the script compiles, but it will fail at runtime if rp isn't fixed.
    print("WARNING: Using PyTorch placeholder. Ensure rotamer_parameter_estimation is updated.")
    try:
        d_obj = get_d_obj_torch()
    except Exception as e:
        print(f"Could not load Torch objective: {e}")
        d_obj = None

    def get_init_param(init_dir):
        init_param_files = dict(
                env = os.path.join(init_dir, 'environment.h5'),
                rot = os.path.join(init_dir, 'sidechain.h5'),
                hb  = os.path.join(init_dir, 'hbond.h5'),
                sheet  = os.path.join(init_dir, 'sheet'))

        # Read Rotamers
        with tb.open_file(init_param_files['rot']) as t:
            rotp = t.root.pair_interaction[:]
            covp = t.root.coverage_interaction[:]
            hydp = t.root.hydrophobe_interaction[:]
            hydplp = t.root.hydrophobe_placement[:]
            rotposp = t.root.rotamer_center_fixed[:]
            rotscalarp = np.zeros((rotposp.shape[0],))

        # Read Environment
        with tb.open_file(init_param_files['env']) as t:
            env = t.root.energies[:,:-1]

        # Read HBond (H5 format)
        with tb.open_file(init_param_files['hb']) as t:
            hb_val = t.root.parameter[0] 

        # --- CHANGED: Read Sheet using numpy (handles vectors) ---
        sheet = np.loadtxt(init_param_files['sheet'])
        # -------------------------------------------------------

        param = Update(*([None]*6))._replace(
                env = env,
                rot = rp.pack_param(rotp, covp, hydp, hydplp, rotposp, rotscalarp[:,None]),
                hb  = hb_val,
                sheet = sheet)
        return param, init_param_files


    def expand_param(params, orig_param_files, new_param_files):
        rotp, covp, hydp, hydplp, rotposp, rotscalarp = rp.unpack_params(params.rot)

        # Write Rotamer H5
        shutil.copyfile(orig_param_files['rot'], new_param_files['rot'])
        with tb.open_file(new_param_files['rot'], 'a') as t:
            t.root.pair_interaction[:]       = rotp
            t.root.coverage_interaction[:]   = covp
            t.root.hydrophobe_interaction[:] = hydp
            t.root.hydrophobe_placement[:]   = hydplp
            t.root.rotamer_center_fixed[:]   = rotposp

        # Write HBond H5
        shutil.copyfile(orig_param_files['hb'], new_param_files['hb'])
        with tb.open_file(new_param_files['hb'], 'a') as t:
            p = t.root.parameter[:]
            p[0] = params.hb
            t.root.parameter[:] = p
        
        # --- CHANGED: Write Sheet using numpy (handles vectors) ---
        # If it's a scalar, make it iterable for savetxt or just save it
        if np.ndim(params.sheet) == 0:
            with open(new_param_files['sheet'], 'w') as f: print(params.sheet, file=f)
        else:
            np.savetxt(new_param_files['sheet'], params.sheet)
        # --------------------------------------------------------

        # Write Environment H5
        shutil.copyfile(orig_param_files['env'], new_param_files['env'])
        with tb.open_file(new_param_files['env'], 'a') as t:
            tmp = np.zeros(t.root.energies.shape)
            tmp[:,:-1] = params.env
            tmp[:,-1]  = tmp[:,-3]
            t.root.energies[:] = tmp


    def backprop_deriv(param, deriv_update, reg_scale):
        # Robust handling for env array dimensions
        if hasattr(deriv_update.env, 'ndim') and deriv_update.env.ndim > 1:
            envd = deriv_update.env[:, :-1].copy()
            envd[:, -2] += deriv_update.env[:, -1]
        else:
            envd = np.copy(deriv_update.env)
        
        # d_obj is now the torch function
        rot_grad = d_obj(param.rot, deriv_update.rot, deriv_update.cov, deriv_update.hyd, reg_scale)

        # Cap gradients at +/- 1.0 (or 5.0). 
        # This prevents infinite repulsion forces from updating the parameters.
        print(f"DEBUG: Max rot_grad before clip: {np.max(np.abs(rot_grad)):.4e}")
        rot_grad = np.clip(rot_grad, -1.0, 1.0)
        
        # Debug: check gradient stats
        print(f"DEBUG: rot_grad mean={np.mean(rot_grad):.4e} min={np.min(rot_grad):.4e} max={np.max(rot_grad):.4e}")
        if np.isnan(rot_grad).any():
             print("ERROR: NaNs detected in rot_grad!")
             raise ExplosionError("NaNs in gradient")

        return deriv_update._replace(
                rot = rot_grad,
                cov = 0.,
                hyd = 0.,
                env = envd)


def run_minibatch(worker_path, param, initial_param_files, direc, minibatch, solver, reg_scale, sim_time):
    if not os.path.exists(direc): os.mkdir(direc)
    print(direc)
    print()

    d_obj_param_files = dict([(k,os.path.join(direc, 'nesterov_temp__'+os.path.basename(x))) 
        for k,x in initial_param_files.items()])
    expand_param(param+solver.update_for_d_obj(), initial_param_files, d_obj_param_files)

    with open(os.path.join(direc,'sim_time'),'w') as f:
            print(sim_time, file=f)

    # REPLACEMENT: Use subprocess directly (Local execution) instead of srun
    jobs = collections.OrderedDict()

    rmsd = dict()
    change = []
    
    # Split the minibatch into manageable chunks
    # This prevents running out of memory/cores
    for i in range(0, len(minibatch), max_parallel_jobs):
        chunk = minibatch[i : i + max_parallel_jobs]
        jobs = collections.OrderedDict()
        
        print(f"--- Running chunk {i//max_parallel_jobs + 1} ({len(chunk)} proteins) ---")

        # 1. LAUNCH CHUNK
        for nm,t in chunk:
            args = [sys.executable, worker_path,
                    'worker', nm, direc, t.fasta, t.init_path, str(t.n_res), t.chi,
                    cp.dumps(d_obj_param_files, protocol=0).decode('latin1'),
                    str(sim_time)]
            
            pickled_files_hex = cp.dumps(d_obj_param_files).hex()
            args[9] = pickled_files_hex

            outfile = open(f'{direc}/{nm}.output_worker', 'w')
            j = sp.Popen(args, close_fds=True, stdout=outfile, stderr=sp.STDOUT)
            jobs[nm] = j

        # 2. COLLECT CHUNK (Wait for these to finish before starting next chunk)
        for nm, j in jobs.items():
            retcode = j.wait()
            if retcode != 0:
                print(f"WARNING: Worker for {nm} exited with code {retcode}.")
            try:
                # print(f"Reading output for {nm}...") # Optional: Comment out to reduce clutter
                with open(f'{direc}/{nm}.divergence.pkl', 'rb') as f:
                    divergence = cp.load(f)

                    # Explosion Check
                    grad_norm = np.linalg.norm(divergence['contrast'].rot) 
                    if grad_norm > 1000.0 or np.isnan(grad_norm):
                        raise ExplosionError(f"Gradient explosion detected in {nm}")

                    rmsd[nm] = (divergence['rmsd_restrain'], divergence['rmsd'])
                    change.append(divergence['contrast'])
            except FileNotFoundError:
                # raise ExplosionError(f"Worker crashed for {nm}")
                # For safety, let's treat missing files as explosions to trigger retry
                raise ExplosionError(f"Worker crashed (no output) for {nm}")
            except Exception as e:
                raise ExplosionError(f"Error reading output for {nm}: {e}")
        
        # End of Chunk Loop - The code loops back to launch the next 12
            
    if not change:
        raise RuntimeError('All jobs failed')

    d_param = backprop_deriv(param, Update(*[np.sum(x,axis=0) for x in zip(*change)]), reg_scale)

    # --- NEW: Global Norm Rescaling (Preserves Physics Direction) ---
    # 1. Calculate the total norm of the gradient vector
    total_norm_sq = 0.0
    for field in ['cov', 'rot', 'hyd', 'env', 'hb', 'sheet']:
        val = getattr(d_param, field)
        if val is not None:
            total_norm_sq += np.sum(val**2)
    total_norm = np.sqrt(total_norm_sq)

    # 2. Rescale if the norm exceeds a safe threshold (e.g., 1.0)
    #    This allows you to use a high Learning Rate safely.
    #    If the gradient is huge (1000.0), it gets scaled down to 1.0,
    #    but the RATIO between parameters remains perfect.
    max_grad_norm = 1.0
    if total_norm > max_grad_norm:
        scale_factor = max_grad_norm / total_norm
        print(f"WARNING: Rescaling massive gradient (Norm: {total_norm:.2f} -> {max_grad_norm})")
        
        for field in ['cov', 'rot', 'hyd', 'env', 'hb', 'sheet']:
            val = getattr(d_param, field)
            if val is not None:
                scaled_val = val * scale_factor
                d_param = d_param._replace(**{field: scaled_val})
    # -------------------------------------------------------------

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

    print(f"Starting compute_divergence for {os.path.basename(config_base)} with {pos.shape[0]} frames...")

    for i in range(pos.shape[0]):
        if i % 50 == 0: print(f"Processing frame {i}/{pos.shape[0]}...")
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
    # Decode hex pickle from argument
    param_files = cp.loads(bytes.fromhex(sys.argv[8]))
    sim_time    = float(sys.argv[9])

    n_frame = 250.
    frame_interval = int(sim_time / n_frame)
    if frame_interval < 1: frame_interval = 1

    #with open(param_files['hb'])    as f: hb        = float(f.read())
    #with open(param_files['sheet']) as f: sheet_mix = float(f.read())
    hb = param_files['hb']
    sheet_mix = param_files['sheet']
    
    # Calculate project root for correct pathing
    project_root = os.path.abspath(os.path.join(os.getcwd(), "../.."))

    # Map ConDiv keys to run_upside.py keys
    kwargs = dict(
            environment_potential = param_files['env'],    # Changed from 'environment'
            rotamer_interaction   = param_files['rot'],    # Changed from 'rotamer_interaction_param'
            rotamer_placement     = param_files['rot'],    # Changed from 'placement'
            initial_structure     = init_path,             # Changed from 'init'
            hbond_energy          = hb,                    # Changed from 'hbond'
            dynamic_rotamer_1body = True,                  # Changed from 'dynamic_1body'
            rama_sheet_mix_energy = sheet_mix,             # Changed from 'sheet_mix_energy'
            rama_param_deriv      = True,
            environment_potential_type = "0",  
            
            # Paths to common files
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

        # Build command manually to bypass ru.upside_config wrapper limitations
        # 1. Base arguments from kwargs
        upside_config_script = os.path.join(os.path.dirname(ru.__file__), 'upside_config.py')
        
        cmd = [sys.executable, upside_config_script]
        
        # Add mapped kwargs as flags
        # Note: We need to convert the keys back to CLI flags. 
        # This mapping must match what we fixed in Step 13.
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

        # Add standard flags
        cmd.append(f'--fasta={fasta}')
        cmd.append(f'--output={configs[0]}')
        
        for k, v in kwargs.items():
            if k in key_map:
                if v is True: # Bool flags like dynamic-rotamer-1body
                    cmd.append(key_map[k])
                elif v is not None:
                    cmd.append(f'{key_map[k]}={v}')

        # Add Restraint flags
        # Restraining all residues 0 to N-1
        cmd.append(f'--restraint-group=0-{n_res-1}')
        cmd.append(f'--restraint-spring={native_restraint_strength}')

        # Run it
        try:
            sp.check_output(cmd, stderr=sp.STDOUT)
        except sp.CalledProcessError as e:
            print("Error generating restrained config:")
            print(e.output.decode('ascii'))
            raise RuntimeError('CONFIG_FAIL_RESTRAINT')
    except tb.NoSuchNodeError:
        raise RuntimeError('CONFIG_FAIL')
    
    # Run Upside
    # n_threads set globally at top.
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
        pos_free = o.pos[int(n_frame/2):,0]
        # Safely compute swap stats if enough frames are available
        if len(o.replica_cumulative_swaps) > int(n_frame/2):
            swap_stats.extend(o.replica_cumulative_swaps[-1] - o.replica_cumulative_swaps[int(n_frame/2)])
        else:
            # Fallback: use zeros of appropriate shape
            swap_stats.extend(np.zeros_like(o.replica_cumulative_swaps[0]))

    for nrep in range(3, len(T), 2):
        with tb.open_file(configs[nrep]) as t:
            o = t.root.output
            if len(o.replica_cumulative_swaps) > int(n_frame/2):
                swap_stats.extend(o.replica_cumulative_swaps[-1] - o.replica_cumulative_swaps[int(n_frame/2)])
            else:
                swap_stats.extend(np.zeros_like(o.replica_cumulative_swaps[0]))

    divergence = dict()
    alldiv = compute_divergence(config_base, np.concatenate([pos_restrain,pos_free],axis=0))
    divergence['contrast'] = Update(*[x[:len(pos_restrain)].mean(axis=0) -
                                      x[len(pos_restrain):].mean(axis=0) for x in alldiv])

    divergence['rmsd_restrain']= ru.traj_rmsd(pos_restrain[:,rmsd_k:-rmsd_k], target[rmsd_k:-rmsd_k]).mean()
    divergence['rmsd']         = ru.traj_rmsd(pos_free    [:,rmsd_k:-rmsd_k], target[rmsd_k:-rmsd_k]).mean()
    divergence['walltime'] = time.time()-tstart
    divergence['swap_stats'] = np.array(swap_stats)

    for fn in [config_base] + configs:
        os.remove(fn)

    with open('%s/%s.divergence.pkl' % (direc,code),'wb') as f:
        cp.dump(divergence, f, -1)

"""
def main_loop(state_str, max_iter):
    def main_loop_iteration(state):
        print('#########################################')
        print('####      EPOCH %2i MINIBATCH %2i      ####' % (state['epoch'],state['i_mb']))
        print('#########################################')
        print()
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
        print()
        sys.stdout.flush()

        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1

    for i in range(max_iter):
        if hasattr(state_str, 'encode'): # Handle string path input
             if os.path.exists(state_str):
                 with open(state_str, 'rb') as f:
                     state = cp.load(f)
             else:
                 # --- NEW ERROR HANDLING ---
                 raise FileNotFoundError(f"Checkpoint file not found: {state_str}\n"
                                         "Did you run the 'initialize' mode first?")
                 # --------------------------
        else:
             state = cp.loads(state_str)

        main_loop_iteration(state)
        
        # Save state
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'),'wb') as f:
            cp.dump(state, f, -1)
        
        # Reload state from bytes for next iter to ensure consistency
        state_str = cp.dumps(state, -1)
"""
def main_loop(state_str, max_iter):
    # Load state initially
    if hasattr(state_str, 'encode'): 
         if os.path.exists(state_str):
             with open(state_str, 'rb') as f: state = cp.load(f)
         else: raise FileNotFoundError(f"Checkpoint {state_str} not found")
    else: state = cp.loads(state_str)

    # We need a copy of the "safe" state to revert to if things blow up
    safe_state_param = state['param'] 

    for i in range(max_iter):
        
        # --- RETRY LOOP ---
        while True:
            try:
                print('#########################################')
                print('####      EPOCH %2i MINIBATCH %2i      ####' % (state['epoch'],state['i_mb']))
                print('#########################################')
                
                tstart = time.time()
                state['mb_direc'] = os.path.join(state['base_dir'],'epoch_%02i_minibatch_%02i'%(state['epoch'],state['i_mb']))
                if os.path.exists(state['mb_direc']): shutil.rmtree(state['mb_direc'], ignore_errors=True)
                                
                # Curriculum: Ramp up time to prevent early explosions
                if state['epoch'] < 5:
                    current_sim_time = 500.0
                elif state['epoch'] < 10:
                    current_sim_time = 1000.0
                else:
                    current_sim_time = state['sim_time']
                print(f"Current Sim Time: {current_sim_time}")

                # Try to run the minibatch
                # NOTE: We pass 'safe_state_param', not 'state['param']' to ensure we start from a clean slate
                new_param = run_minibatch(state['worker_path'], safe_state_param, state['init_param_files'],
                        state['mb_direc'], state['minibatches'][state['i_mb']],
                        state['solver'], state['n_prot']*1., current_sim_time)
                
                # If we get here, it succeeded! 
                state['param'] = new_param       # Update the official param
                safe_state_param = new_param     # Update our "safe" backup
                break                            # Exit the Retry Loop
                
            except ExplosionError as e:
                print(f"\n!!! EXPLOSION DETECTED: {e} !!!")
                print("!!! ACTION: Reducing Learning Rate by 50% and Retrying Batch !!!\n")
                
                # 1. Adjust Parameters (Reduce Learning Rate)
                # Assuming state['solver'] is an AdamSolver instance
                # We multiply its alpha (learning rate) by 0.5
                state['solver'].alpha = state['solver'].alpha * 0.5
                print(f"NEW LEARNING RATE: {state['solver'].alpha}")
                
                # 2. Restart (Cleanup)
                # We do NOT update state['i_mb']. We simply loop back.
                # 'safe_state_param' has NOT been updated, so we automatically revert to the start of the batch.
                if os.path.exists(state['mb_direc']):
                    shutil.rmtree(state['mb_direc'], ignore_errors=True)
                
                # Optional: If Learning rate gets too small, abort completely
                if np.max(state['solver'].alpha) < 1e-6:
                    raise RuntimeError("Learning rate dropped too low. Training is impossible.")
                
                time.sleep(2) # Breathe before retrying
        # ------------------

        print('%.0f seconds elapsed this minibatch'%(time.time()-tstart))
        sys.stdout.flush()

        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1
        
        # Save checkpoint
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'),'wb') as f:
            cp.dump(state, f, -1)


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
                native_pos = cp.load(f, encoding='latin1')[:,:,0] # Handle Py2 numpy pickles
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

    minibatch_excess = len(training_list)%minibatch_size
    if minibatch_excess: training_list = training_list[:-minibatch_excess]
    n_minibatch = len(training_list)//minibatch_size
    state['minibatches'] = [training_list[i::n_minibatch] for i in range(n_minibatch)]
    state['n_prot'] = n_minibatch*minibatch_size
    print('Constructed %i minibatches of size %i (%i proteins)' % (n_minibatch, minibatch_size, state['n_prot']))


    if state['init_dir'] != 'cached':
        print('about to get init')
        state['param'], state['init_param_files'] = get_init_param(state['init_dir'])
        print('found init')
        with open(os.path.join(state['base_dir'], 'condiv_init.pkl'),'wb') as f:
            cp.dump((state['init_dir'],state['param'],state['init_param_files']),f,-1)
    else:
        with open(os.path.join(state['base_dir'], 'condiv_init.pkl'), 'rb') as f:
            state['init_dir'],state['param'],state['init_param_files'] = cp.load(f)

    state['initial_alpha'] = Update(*[
            0.1, 0., 0.5, 0., 0.02, 0.03])
    #state['initial_alpha'] = state['initial_alpha'] * 0.025
    state['initial_alpha'] = state['initial_alpha'] * 0.05
    state['solver'] = rp.AdamSolver(len(state['initial_alpha']), alpha=state['initial_alpha']) 
    state['sim_time'] = 4000

    print()
    print('Optimizing with solver', state['solver'])
    print()
    print_param(state['param'])
    print()

    state['epoch'] = 0
    state['i_mb'] = 0
    return state


def main_debug():
    """Debug mode: Run with 1 protein for 1000.0 simulation time (about 100 steps)"""
    
    # Set hardcoded paths as requested in the prompt
    init_dir = "/Users/yinhan/Documents/upside2-md/parameters/ff_2.1"
    protein_dir = "/Users/yinhan/Documents/upside2-md/training_set"
    protein_list = os.path.join(os.path.dirname(__file__), "pdb_list")  # Use the pdb_list in current directory
    base_dir = os.path.join(os.path.dirname(__file__), "debug_mode_run")
    
    # Call main_initialize with these parameters
    args = [init_dir, protein_dir, protein_list, base_dir]
    state = main_initialize(args)
    
    # Modify state for debug mode
    state['minibatch_size'] = 1
    state['sim_time'] = 1000.0  # 100 steps
    
    # Load only one protein from the training set
    print("\n=== DEBUG MODE: Loading only 1 protein from training set ===")
    # Take the first protein from the first minibatch
    first_protein = state['minibatches'][0][0]
    state['minibatches'] = [[first_protein]]  # Create a single minibatch with just one protein
    
    # Save the modified state
    debug_checkpoint = os.path.join(base_dir, "debug_initial_checkpoint.pkl")
    with open(debug_checkpoint, 'wb') as f:
        cp.dump(state, f, -1)
    
    print("\n=== DEBUG MODE: Running simulation with 100 steps ===")
    print(f"Debug checkpoint saved at: {debug_checkpoint}")
    print(f"Simulation time set to: {state['sim_time']}")
    print(f"Number of proteins: 1")
    
    # Run one iteration of the main loop
    main_loop(cp.dumps(state, -1), 1)


if __name__ == '__main__':
    # Run debug mode by default if no mode is specified
    if len(sys.argv) < 2 or sys.argv[1] not in ['worker', 'restart', 'initialize']:
        main_debug()
    elif sys.argv[1] == 'worker':
        main_worker()

    elif sys.argv[1] == 'restart':
        assert len(sys.argv[1:]) == 3
        print('Running as PID %i on host %s' % (os.getpid(), socket.gethostname()))
        # Restart expects a checkpoint FILE path now
        main_loop(sys.argv[2], int(sys.argv[3]))

    elif sys.argv[1] == 'initialize':
        initial_state = main_initialize(sys.argv[2:])
        with open(os.path.join(initial_state['base_dir'],'initial_checkpoint.pkl'),'wb') as f:
            cp.dump(initial_state, f, -1)
    if sys.argv[1] == 'worker':
        main_worker()

    elif sys.argv[1] == 'restart':
        assert len(sys.argv[1:]) == 3
        print('Running as PID %i on host %s' % (os.getpid(), socket.gethostname()))
        # Restart expects a checkpoint FILE path now
        main_loop(sys.argv[2], int(sys.argv[3]))

    elif sys.argv[1] == 'initialize':
        initial_state = main_initialize(sys.argv[2:])
        with open(os.path.join(initial_state['base_dir'],'initial_checkpoint.pkl'),'wb') as f:
            cp.dump(initial_state, f, -1)

    else:
        raise RuntimeError('Illegal mode %s'%sys.argv[1])