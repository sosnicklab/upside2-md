#!/usr/bin/env python
import sys,os

rp_path = '/home/pengxd/software/source/upside-devel/src'
if rp_path not in sys.path:
    sys.path.append(rp_path)

os.environ['THEANO_FLAGS'] = '' # 'cxx=,optimizer=fast_compile'
is_worker = __name__ == '__main__' and sys.argv[1]=='worker'
if not is_worker:
    import rotamer_parameter_estimation as rp

import numpy as np
import subprocess as sp
import tables as tb
import cPickle as cp
from glob import glob
import re
import shutil
import collections
import time

import run_upside as ru
import upside_engine as ue
import socket

np.set_printoptions(precision=2, suppress=True)

## Important parameters
n_threads = 8
native_restraint_strength = 1./3.**2  # weak restraints aiming to hold at about 1.0A RMSD
rmsd_k = 15   # number of atoms to cut
minibatch_size = 12

resnames = ['ALA', 'ARG', 'ASN', 'ASP',
            'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS',
            'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL']

hydrophobicity_order = ['ASP', 'GLU', 'LYS', 'HIS',
                        'ARG', 'GLY', 'ASN', 'GLN',
                        'ALA', 'SER', 'THR', 'PRO',
                        'CYS', 'VAL', 'MET', 'TYR',
                        'ILE', 'LEU', 'PHE', 'TRP']

Target = collections.namedtuple('Target', 'fasta native native_path init_path n_res chi'.split())
UpdateBase = collections.namedtuple('UpdateBase', 'env cov rot hyd hb sheet'.split())
class Update(UpdateBase):
    def __init__(self, *args):
        super(Update, self).__init__(*args)

    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False

        if is_seq:
            ret = []
            assert len(self) == len(other)
            for a,b in zip(self,other):
                if a is None or b is None:   # treat None as an non-existent value and propagate
                    ret.append(None)
                else:
                    ret.append(op(a,b))
        else:
            ret = [None if a is None else op(a,other) for a in self]
        return Update(*ret)

    def __add__(self, other):
        return self._do_binary(other, lambda a,b: a+b)
    def __sub__(self, other):
        return self._do_binary(other, lambda a,b: a-b)
    def __mul__(self, other):
        return self._do_binary(other, lambda a,b: a*b)
    def __div__(self, other):
        return self._do_binary(other, lambda a,b: a/b)


def print_param(param):
    print 'hb    %.3f'% param.hb
    print 'sheet %.3f'% param.sheet
    
    print 'env'
    env_dict = dict(zip(resnames, param.env[:,1::2]))
    for r in hydrophobicity_order:   # easier to understand in hydrophobicity order
        print '   ',r,env_dict[r]

def get_d_obj():
    import theano.tensor as T
    import theano

    cov_wide   = rp.unpack_cov_expr[:,:,2*rp.n_knot_angular:2*rp.n_knot_angular+rp.n_knot_hb]
    cov_narrow = rp.unpack_cov_expr[:,:,2*rp.n_knot_angular+rp.n_knot_hb:]
    cov_rHN    = rp.unpack_cov_expr[:,:,:rp.n_knot_angular]
    cov_rSC    = rp.unpack_cov_expr[:,:,rp.n_knot_angular:2*rp.n_knot_angular]
    
    hyd_wide   = rp.unpack_hyd_expr[:,:,2*rp.n_knot_angular:2*rp.n_knot_angular+rp.n_knot_hb]
    hyd_narrow = rp.unpack_hyd_expr[:,:,2*rp.n_knot_angular+rp.n_knot_hb:]
    hyd_rHN    = rp.unpack_hyd_expr[:,:,:rp.n_knot_angular]
    hyd_rSC    = rp.unpack_hyd_expr[:,:,rp.n_knot_angular:2*rp.n_knot_angular]
    
    rot_wide   = rp.unpack_rot_expr[:,:,2*rp.n_knot_angular:2*rp.n_knot_angular+rp.n_knot_sc]
    rot_narrow = rp.unpack_rot_expr[:,:,2*rp.n_knot_angular+rp.n_knot_sc:]
    rot_rSC1   = rp.unpack_rot_expr[:,:,:rp.n_knot_angular]
    rot_rSC2   = rp.unpack_rot_expr[:,:,rp.n_knot_angular:2*rp.n_knot_angular]

    def student_t_neglog(t, scale, nu):
        return 0.5*(nu+1.)*T.log(1.+(1./(nu*scale**2)) * t**2)

    hb_n_knot   = (rp.n_knot_angular,rp.n_knot_angular, rp.n_knot_hb, rp.n_knot_hb)
    sc_n_knot   = (rp.n_knot_angular,rp.n_knot_angular, rp.n_knot_sc, rp.n_knot_sc)
    cov_energy  = rp.quadspline_energy(rp.unpack_cov_expr, hb_n_knot)
    rot_energy  = rp.quadspline_energy(rp.unpack_rot_expr, sc_n_knot)
    hyd_energy  = rp.quadspline_energy(rp.unpack_hyd_expr, hb_n_knot)

    r_for_cov_energy = np.arange(1,rp.n_knot_hb-1)*rp.hb_dr
    r_for_rot_energy = np.arange(1,rp.n_knot_sc-1)*rp.sc_dr
    r_for_hyd_energy = np.arange(1,rp.n_knot_hb-1)*rp.hb_dr
    
    def cutoff_func(x, scale, module=np):
        return 1./(1.+module.exp(x/scale))
    
    def lower_bound(x, lb):
        return T.where(x<lb, 1e0*(x-lb)**2, 0.*x)
    
    cov_expect = 5.*cutoff_func(r_for_cov_energy-2., 0.2)
    rot_expect = 5.*cutoff_func(r_for_rot_energy-2., 0.2)
    hyd_expect = 5.*cutoff_func(r_for_hyd_energy-1., 0.2)
    
    reg_scale = T.dscalar('reg_scale')
    cov_reg = student_t_neglog(cov_energy - cov_expect[None,None,:,None,None], 200., 3.).sum()/reg_scale
    rot_reg = student_t_neglog(rot_energy - rot_expect[None,None,:,None,None], 200., 3.).sum()/reg_scale
    hyd_reg = student_t_neglog(hyd_energy - hyd_expect[None,None,:,None,None], 200., 3.).sum()/reg_scale
    
    reg_expr = (
            cov_reg + 
            hyd_reg +
            rot_reg + 
    
            # allow no energies below -5. to avoid weird pathologies
            lower_bound(cov_energy, -6.).sum() + 
            lower_bound(rot_energy, -6.).sum() + 
            lower_bound(hyd_energy, -6.).sum() + 
    
            0. #8
            )  #6

    # For ease of use (and avoiding double evaluation of simulations for derivs), we don't want 
    # to make Upside run within the theano evaulation.  We still want theano to handle the deriv
    # propagation for us.  The way to do this is to define an auxiliary objective f2 instead of the 
    # real objective f using the identity (d/dx)(f(g(x))) = (d/dx)(c*g(x)) if c == f'(g(x)).  I call
    # the c the "coupling".  
    rot_coupling = T.dtensor3()
    cov_coupling = T.dtensor3()
    hyd_coupling = T.dtensor3()

    obj_expr = ((rot_coupling*rp.unpack_rot_expr).sum() +
                (cov_coupling*rp.unpack_cov_expr).sum() +
                (hyd_coupling*rp.unpack_hyd_expr).sum() +
                reg_expr)

    d_obj = theano.function([rp.lparam, rot_coupling, cov_coupling, hyd_coupling, reg_scale],
            T.grad(obj_expr, rp.lparam))
    return d_obj


if not is_worker:
    d_obj = get_d_obj()

    def get_init_param(init_dir):
        init_param_files = dict(
                env = os.path.join(init_dir, 'environment.h5'),
                rot = os.path.join(init_dir, 'sidechain.h5'),
                hb  = os.path.join(init_dir, 'hbond'),
                sheet  = os.path.join(init_dir, 'sheet'))

        with tb.open_file(init_param_files['rot']) as t:
            rotp = t.root.pair_interaction[:]
            covp = t.root.coverage_interaction[:]
            hydp = t.root.hydrophobe_interaction[:]
            hydplp = t.root.hydrophobe_placement[:]
            rotposp = t.root.rotamer_center_fixed[:]
            rotscalarp = np.zeros((rotposp.shape[0],)) #t.root.rotamer_prob_fixed[:]

        with tb.open_file(init_param_files['env']) as t:
            env = t.root.energies[:,:-1]  # don't need last elem due to zero-deriv clamping

        with open(init_param_files['hb']) as f: hb = float(f.read())
        with open(init_param_files['sheet']) as f: sheet = float(f.read())

        param = Update(*([None]*6))._replace(
                env = env,
                #rot = rp.pack_param(rotp, covp, hydp, hydplp, rotposp),
                rot = rp.pack_param(rotp, covp, hydp, hydplp, rotposp, rotscalarp[:,None]),
                hb  = hb,
                sheet = sheet)
        return param, init_param_files


    def expand_param(params, orig_param_files, new_param_files):
        rotp, covp, hydp, hydplp, rotposp, rotscalarp = rp.unpack_params(params.rot)

        shutil.copyfile(orig_param_files['rot'], new_param_files['rot'])  # both rot and cov are in the same file
        with tb.open_file(new_param_files['rot'],'a') as t:
            t.root.pair_interaction[:]       = rotp
            t.root.coverage_interaction[:]   = covp
            t.root.hydrophobe_interaction[:] = hydp
            t.root.hydrophobe_placement[:]   = hydplp
            t.root.rotamer_center_fixed[:]   = rotposp
            # t.root.rotamer_prob_fixed[:]     = rotscalarp

        # read sheet and hb from file here and update
        with open(new_param_files['hb'],'w') as f:    print >>f, params.hb
        with open(new_param_files['sheet'],'w') as f: print >>f, params.sheet

        shutil.copyfile(orig_param_files['env'], new_param_files['env'])
        with tb.open_file(new_param_files['env'],'a') as t:
            tmp = np.zeros(t.root.energies.shape)
            tmp[:,:-1] = params.env
            tmp[:,-1]  = tmp[:,-3]   # impose zero clamp on right side only
            t.root.energies[:] = tmp


    def backprop_deriv(param, deriv_update, reg_scale):
        # place all of the derivatives in rot for good measure
        envd = deriv_update.env[:,:-1].copy()
        envd[:,-2] += deriv_update.env[:,-1]   # from zero clamp condition
        return deriv_update._replace(
                rot = d_obj(param.rot, deriv_update.rot, deriv_update.cov, deriv_update.hyd, reg_scale),
                cov = 0.,
                hyd = 0.,
                env = envd)


def run_minibatch(worker_path, param, initial_param_files, direc, minibatch, solver, reg_scale, sim_time):
    if not os.path.exists(direc): os.mkdir(direc)
    print direc
    print

    d_obj_param_files = dict([(k,os.path.join(direc, 'nesterov_temp__'+os.path.basename(x))) 
        for k,x in initial_param_files.items()])
    expand_param(param+solver.update_for_d_obj(), initial_param_files, d_obj_param_files)

    with open(os.path.join(direc,'sim_time'),'w') as f:
            print >>f, sim_time

    jobs = collections.OrderedDict()
    for nm,t in minibatch[::-1]:
        args = ['srun', '--nodes=1', '--ntasks=1', '--cpus-per-task=%i'%n_threads, '--slurmd-debug=0', 
                '--output=%s/%s.output_worker'%(direc,nm), 
                worker_path,
                'worker', nm, direc, t.fasta, t.init_path, str(t.n_res), t.chi,
                cp.dumps(d_obj_param_files), str(sim_time)]

        jobs[nm] = sp.Popen(args, close_fds=True)

    rmsd = dict()
    change = []
    for nm,j in jobs.items():
        if j.wait() != 0: 
            print nm, 'WORKER_FAIL'
            continue
        with open('%s/%s.divergence.pkl'%(direc,nm)) as f:
            divergence = cp.load(f)
            rmsd[nm] = (divergence['rmsd_restrain'], divergence['rmsd'])
            change.append(divergence['contrast'])
    if not change:
        raise RuntimeError('All jobs failed')

    d_param = backprop_deriv(param, Update(*[np.sum(x,axis=0) for x in zip(*change)]), reg_scale)

    with open('%s/rmsd.pkl' % (direc,),'w') as f:
        cp.dump(rmsd, f, -1)
    print
    print 'Median RMSD %.2f %.2f' % tuple(np.median(np.array(rmsd.values()), axis=0))

    ## Update the parameters
    new_param_files = dict([(k,os.path.join(direc, os.path.basename(x))) for k,x in initial_param_files.items()])
    new_param = param+solver.update_step(d_param)
    expand_param(new_param, initial_param_files, new_param_files)
    solver.log_state(direc)

    print
    print_param(new_param)

    return new_param


def compute_divergence(config_base, pos):
    try:
        with tb.open_file(config_base) as t:
            seq        = t.root.input.sequence[:]
            more_sheet = t.root.input.potential.rama_map_pot.more_sheet_rama_pot[:]
            less_sheet = t.root.input.potential.rama_map_pot.less_sheet_rama_pot[:]
            eps        = t.root.input.potential.rama_map_pot._v_attrs.sheet_eps
            sheet_scale = 1./(2.*eps)
            hb_strength = t.root.input.potential.hbond_energy._v_attrs.protein_hbond_energy
            rot_param_shape = t.root.input.potential.rotamer.pair_interaction.interaction_param.shape
            cov_param_shape = t.root.input.potential.hbond_coverage.interaction_param.shape
            hyd_param_shape = t.root.input.potential.hbond_coverage_hydrophobe.interaction_param.shape
            env_param_shape = t.root.input.potential.nonlinear_coupling_environment.coeff.shape
    except Exception as e:
        print(os.path.basename(config_base)[:5],'ANALYSIS_FAIL',e)
        return None

    #engine = ue.Upside(pos.shape[1], config_base)
    engine = ue.Upside(config_base)
    contrast = Update([],[],[],[],[],[])

    engine.set_param(more_sheet, 'rama_map_pot')
    for i in range(pos.shape[0]):
        engine.energy(pos[i])

        # logarithmic derivative of scale factor
        contrast.cov.append(engine.get_param_deriv(cov_param_shape, 'hbond_coverage'))
        contrast.rot.append(engine.get_param_deriv(rot_param_shape, 'rotamer'))
        contrast.hyd.append(engine.get_param_deriv(hyd_param_shape, 'hbond_coverage_hydrophobe'))

        contrast.hb   .append(engine.get_output('hbond_energy')[0,0]/hb_strength)
        contrast.sheet.append(engine.get_output('rama_map_pot')[0,0]*sheet_scale)

        contrast.env.append(engine.get_param_deriv(env_param_shape, 'nonlinear_coupling_environment'))

    engine.set_param(less_sheet, 'rama_map_pot')
    for i in range(pos.shape[0]):
        engine.energy(pos[i])
        # reverse sign for finite difference
        contrast.sheet[i] -= engine.get_output('rama_map_pot')[0,0]*sheet_scale

    # convert to numpy arrays
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
    param_files = cp.loads(sys.argv[8])
    sim_time    = float(sys.argv[9])

    n_frame = 250.
    frame_interval = int(sim_time / n_frame)

    with open(param_files['hb'])    as f: hb        = float(f.read())
    with open(param_files['sheet']) as f: sheet_mix = float(f.read())
    
    ## Configure the files
    kwargs = dict(
            environment               = param_files['env'],
            rotamer_interaction_param = param_files['rot'],
            placement                 = param_files['rot'],  # placement file is currently the same as the interaction file
            init = init_path,
            hbond = hb,
            dynamic_1body=True,
            sheet_mix_energy=sheet_mix,
            #reference_rama=ru.params_dir+'reference_state_rama.pkl',)
            rama_pot = "/home/pengxd/software/upside-18.10.08/parameters/common/rama.dat",
            reference_rama='/home/pengxd/software/upside-18.10.08/parameters/common/rama_reference.pkl',)
    
    # T = [0.80,0.80, 0.86, 0.91]
    
    # T = 0.80 * (1. + np.sqrt(100./n_res)*0.025*np.arange(n_threads-1))**2
    T = 0.80 * (1. + np.sqrt(100./n_res)*0.020*np.arange(n_threads-1))**2
    T = np.concatenate((T[0:1],T))  # double the first temperature for the restrained replica

    try:
        config_base = '%s/%s.base.h5' % (direc,code)
        ru.upside_config(fasta, config_base, **kwargs)

        configs = [re.sub(r'\.base\.h5','.run.%i.h5'%i_rs, config_base) for i_rs in range(len(T))]
        for i in range(1,len(T)):  # just copy to get different temperatures
            shutil.copyfile(config_base, configs[i])

        # now create the near-native simulation config
        # I will also restrain the chi's here
        kwargs['restraint_groups'] = ['0-%i'%(n_res-1)]
        kwargs['restraint_spring'] = native_restraint_strength
        # kwargs['fix_rotamer']      = chi
        ru.upside_config(fasta, configs[0], **kwargs)
    except tb.NoSuchNodeError:
        raise RuntimeError('CONFIG_FAIL')
    
    ## Launch the jobs
    j = ru.run_upside('', configs, sim_time, frame_interval, n_threads=n_threads, temperature=T,
            swap_sets=ru.swap_table2d(1,len(T)), mc_interval=5., replica_interval=10.)
    if j.job.wait() != 0: raise RuntimeError('RUN_FAIL')

    # record swap_stats only for second half of trajectory and only between free replicas
    swap_stats = []

    ## Compute divergence of the structures
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

    ## Cleanup files on success (let's just leave for triage on failure)
    for fn in [config_base] + configs:
        os.remove(fn)

    ## Write output
    with open('%s/%s.divergence.pkl' % (direc,code),'w') as f:
        cp.dump(divergence, f, -1)


def main_loop(state_str, max_iter):
    def main_loop_iteration(state):
        print '#########################################'
        print '####      EPOCH %2i MINIBATCH %2i      ####' % (state['epoch'],state['i_mb'])
        print '#########################################'
        print
        # decay = lambda x: 1./(1.+0.5*x)
        # epoch_real = state['epoch']*1. + state['i_mb']*(1./len(state['minibatches']))
        # print 'alpha_scale %.2f'%decay(epoch_real)
        # state['solver'].alpha = state['initial_alpha'] * decay(epoch_real)
        sys.stdout.flush()
    
        tstart = time.time()
        state['mb_direc'] = os.path.join(state['base_dir'],'epoch_%02i_minibatch_%02i'%(state['epoch'],state['i_mb']))
        if os.path.exists(state['mb_direc']):  # possibly cleanup an earlier failed execution
            shutil.rmtree(state['mb_direc'], ignore_errors=True)
        state['param'] = run_minibatch(state['worker_path'], state['param'], state['init_param_files'],
                state['mb_direc'], state['minibatches'][state['i_mb']],
                state['solver'], state['n_prot']*1., state['sim_time'])
        print
        print '%.0f seconds elapsed this minibatch'%(time.time()-tstart)
        print
        sys.stdout.flush()

        # increment counters for next minibatch
        state['i_mb'] += 1
        if state['i_mb'] >= len(state['minibatches']):
            state['i_mb'] = 0
            state['epoch'] += 1
            # state['sim_time'] *= 2 # np.sqrt(3.)  # progressively lengthen simulations to explore new areas of phase space

    for i in range(max_iter):
        # Load the state from a string on every iteration to ensure that we are pickling
        #  everything correctly and thus checkpointing is equivalent to not stopping
        state = cp.loads(state_str)
        main_loop_iteration(state)   # mutates state
        state_str = cp.dumps(state,-1)
        with open(os.path.join(state['mb_direc'], 'checkpoint.pkl'),'w') as f:
            f.write(state_str)


def main_initialize(args):

    state = dict()
    state['init_dir'], state['protein_dir'], protein_list, state['base_dir'], = args
    if not os.path.exists(state['base_dir']): os.mkdir(state['base_dir'])

    # Copy this file to the working directory to ensure that it is unchanged during worker invocations
    if not __file__: raise RuntimeError('No file name available')
    state['worker_path'] = os.path.join(state['base_dir'],'ConDiv.py')
    shutil.copy(__file__, state['worker_path'])  # copy preserve execute permission

    ## Read training set
    if protein_list != 'cached':
        print 'Reading training set'
        with open(protein_list) as f:
            protein_names = [x.split()[0] for x in f]
            assert protein_names[0]=='prot'
            protein_names = protein_names[1:]
                
        training_set = dict()
        excluded_prot = []
        for code in sorted(protein_names):
            base = os.path.join(state['protein_dir'], code)

            with open(base+'.initial.pkl') as f:
                native_pos = cp.load(f)[:,:,0]
                n_res = len(native_pos)/3

            max_sep = np.sqrt(np.sum(np.diff(native_pos,axis=0)**2,axis=-1)).max()
            if max_sep < 2.:  # no chain breaks
                training_set[code] = Target(base+'.fasta', native_pos, base+'.initial.pkl',
                    base+'.initial.pkl', n_res, base+'.chi')
            else:
                excluded_prot.append(code)
                print code

        print 'Excluded %i proteins due to chain breaks'%len(excluded_prot)
        with open(os.path.join(state['base_dir'], 'cd_training.pkl'),'w') as f: cp.dump(training_set, f, -1)
    else:
        training_set = cp.load(open(os.path.join(state['base_dir'], 'cd_training.pkl')))

    ## Construct minibatches
    # ensure each minibatch has roughly the same mix of protein sizes for variance reduction
    training_list = sorted(training_set.items(), key=lambda x: (x[1].n_res,x[0]))
    np.random.shuffle(training_list)

    minibatch_excess = len(training_list)%minibatch_size
    if minibatch_excess: training_list = training_list[:-minibatch_excess]
    n_minibatch = len(training_list)//minibatch_size
    state['minibatches'] = [training_list[i::n_minibatch] for i in range(n_minibatch)]
    state['n_prot'] = n_minibatch*minibatch_size
    print 'Constructed %i minibatches of size %i (%i proteins)' % (n_minibatch, minibatch_size, state['n_prot'])


    if state['init_dir'] != 'cached':
        print 'about to get init'
        state['param'], state['init_param_files'] = get_init_param(state['init_dir'])
        print 'found init'
        with open(os.path.join(state['base_dir'], 'condiv_init.pkl'),'w') as f:
            cp.dump((state['init_dir'],state['param'],state['init_param_files']),f,-1)
    else:
        state['init_dir'],state['param'],state['init_param_files'] = cp.load(open(os.path.join(state['base_dir'], 'condiv_init.pkl')))

    state['initial_alpha'] = Update(*[
            0.1,    # env
            0.,     # cov (not needed during backprop because it is in rot)
            0.5,    # rot 
            0.,     # hyd (part of hyd) 
            0.02,   # hb
            0.03])  # sheet
    state['initial_alpha'] = state['initial_alpha'] * 0.25
    state['solver'] = rp.AdamSolver(len(state['initial_alpha']), alpha=state['initial_alpha']) 
    state['sim_time'] = 1000.*4  # about 0.5 hour

    print
    print 'Optimizing with solver', state['solver']
    print

    print
    print_param(state['param'])
    print

    state['epoch'] = 0
    state['i_mb'] = 0
    return state

if __name__ == '__main__':
    if sys.argv[1] == 'worker':
        main_worker()

    elif sys.argv[1] == 'restart':
        assert len(sys.argv[1:]) == 3
        print 'Running as PID %i on host %s' % (os.getpid(), socket.gethostname())
        main_loop(open(sys.argv[2]).read(), int(sys.argv[3]))  # first is checkpoint, 2nd is max_iter

    elif sys.argv[1] == 'initialize':
        initial_state = main_initialize(sys.argv[2:])  # this will dump an initial state file in the working directory
        with open(os.path.join(initial_state['base_dir'],'initial_checkpoint.pkl'),'w') as f:
            cp.dump(initial_state, f, -1)

    else:
        raise RuntimeError('Illegal mode %s.  Please see contrastive_divergence.py for details'%sys.argv[1])
