#!/usr/bin/env python3
import sys, os
import numpy as np
import pickle as cp
import shutil
import collections
import time

# --- Setup Paths (Same as ConDiv.py) ---
rp_path = os.path.abspath(os.path.join(os.getcwd(), "../..", "src")) 
if rp_path not in sys.path:
    sys.path.append(rp_path)

import rotamer_parameter_estimation as rp

np.set_printoptions(precision=2, suppress=True)

# --- Definitions required for Unpickling ---
resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL']

hydrophobicity_order = ['ASP', 'GLU', 'LYS', 'HIS', 'ARG', 'GLY', 'ASN', 'GLN',
                        'ALA', 'SER', 'THR', 'PRO', 'CYS', 'VAL', 'MET', 'TYR',
                        'ILE', 'LEU', 'PHE', 'TRP']

Target = collections.namedtuple('Target', 'fasta native native_path init_path n_res chi'.split())
UpdateBase = collections.namedtuple('UpdateBase', 'env cov rot hyd hb sheet'.split())

class Update(UpdateBase):
    def __init__(self, *args, **kwargs): pass
    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False
        if is_seq:
            ret = []
            for a, b in zip(self, other):
                if a is None or b is None: ret.append(None)
                else: ret.append(op(a, b))
        else:
            ret = [None if a is None else op(a, other) for a in self]
        return Update(*ret)
    def __add__(self, other): return self._do_binary(other, lambda a, b: a + b)
    def __sub__(self, other): return self._do_binary(other, lambda a, b: a - b)
    def __mul__(self, other): return self._do_binary(other, lambda a, b: a * b)
    def __truediv__(self, other): return self._do_binary(other, lambda a, b: a / b)

# --- Helper Functions ---

def print_param(param):
    print(f'hb    {float(param.hb):.3f}')
    if np.ndim(param.sheet) == 0:
        print(f'sheet {float(param.sheet):.3f}')
    else:
        print(f'sheet (mean) {np.mean(param.sheet):.3f}')

    print('env')
    env_data = param.env
    if hasattr(env_data, 'detach'): env_data = env_data.detach().numpy()
        
    try:
        env_dict = dict(zip(resnames, env_data[:, 1::2]))
        for r in hydrophobicity_order:
            print('   ', r, env_dict[r])
    except IndexError:
        print('    (Environment parameters not fully initialized yet)')

def expand_param(params, orig_param_files, new_param_files):
    import tables as tb
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

# --- Main Logic ---

def main_zero_descent(state_str, n_steps=100):
    """
    Quick Descent Algorithm to zero out all force fields.
    """
    print("\n=== STARTING QUICK DESCENT TO ZERO FIELDS ===")
    
    # Load State
    if hasattr(state_str, 'encode') and os.path.exists(state_str):
        print(f"Loading checkpoint: {state_str}")
        with open(state_str, 'rb') as f: state = cp.load(f)
    else:
        state = cp.loads(state_str)

    # Increase learning rate for rapid descent
    quick_alpha = state['initial_alpha'] * 5.0 
    state['solver'].alpha = quick_alpha
    print(f"Solver Alpha boosted for descent: {state['solver'].alpha}")

    for i in range(n_steps):
        param = state['param']
        
        # Gradient = Current Parameters (L2 minimization)
        # FIX: The solver crashes if it receives 'None' for cov/hyd.
        # We replace any None fields with 0.0 float.
        d_param_clean = [p if p is not None else 0.0 for p in param]
        d_param = Update(*d_param_clean)
        
        # Update via Adam
        update_vec = state['solver'].update_step(d_param)
        
        # Apply (Negative update for minimization)
        new_param = param + update_vec
        state['param'] = new_param
        
        # Check progress
        if hasattr(new_param.rot, 'detach'): rot_val = new_param.rot.detach().numpy()
        else: rot_val = new_param.rot
             
        norm_rot = np.linalg.norm(rot_val)
        print(f"Step {i+1}/{n_steps}: Rotamer Norm = {norm_rot:.4f}")
        
        if norm_rot < 1e-4:
            print("Parameters effectively zeroed. Stopping early.")
            break

    # Save Results
    print("Saving zeroed parameters...")
    new_param_files = dict([(k,os.path.join(state['base_dir'], 'zeroed_'+os.path.basename(x))) for k,x in state['init_param_files'].items()])
    expand_param(new_param, state['init_param_files'], new_param_files)
    
    output_checkpoint = os.path.join(state['base_dir'], 'checkpoint_zeroed.pkl')
    with open(output_checkpoint, 'wb') as f:
        cp.dump(state, f, -1)
    
    print()
    print_param(new_param)
    print(f"=== COMPLETE. Saved to {output_checkpoint} ===")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ConDiv_zero.py <checkpoint_file> [steps]")
        sys.exit(1)
    
    checkpoint_file = sys.argv[1]
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    
    main_zero_descent(checkpoint_file, steps)