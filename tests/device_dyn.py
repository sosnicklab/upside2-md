#!/usr/bin/env python
# V2 of the device engine (upside3.md 10.2): one config's short replica-exchange run on the CPU engine and with --gpu,
# compared as ensembles, since trajectories part within a few hundred rounds (findings 11g).  An arm is an engine with
# its options; each runs from the same start at several seeds.  Per replica over the second half of the frames:
# kinetic energy / (1.5 kT), the potential, the backbone bond lengths N-CA, CA-C and C-N and the backbone angles
# N-CA-C, CA-C-N and C-N-CA (degrees); for a hybrid config also
# the protein's and the lipids' kinetic energies, the lipid beads' lateral diffusion coefficient (the xy mean square
# displacement over half the window, per frame pair, over 4 times the lag) and the bilayer thickness (twice the PO4
# beads' mean |z| from their mean).  Each mean is averaged over an arm's seeds with its standard error from their
# spread, and two arms' difference is given in standard errors with its two-sided p (Welch); they agree when the
# smallest p times the number of comparisons (Bonferroni) is above 0.05.  The error comes from independent seeds
# because block averages of one hybrid run understate it: two CPU runs at seeds 7 and 8 differed by up to 7.6 block
# standard errors (findings 11g).  A run keeps its trajectories in a temporary directory and writes its means to
# WORK_DIR/ARM.SEED.npz, which compare reads.  OBJ is the build, a directory of UPSIDE_HOME (obj, or obj_<tag>).
# usage: device_dyn.py run WORK_DIR ARM SEED CONFIG N_REPLICA DURATION OBJ [upside option ...]   (--gpu for the device)
#        device_dyn.py compare WORK_DIR ARM_A ARM_B
import sys, os, re, shutil, tempfile
import subprocess as sp
import numpy as np
import tables as tb
from scipy import stats

upside_path = os.environ['UPSIDE_HOME']

def run(work, arm, seed, config, n_rep, duration, obj, extra):
    temps = ','.join('%.4f' % t for t in np.linspace(0.80, 0.95, n_rep))
    swaps = ['--swap-set', ','.join('%d-%d' % (i, i+1) for i in range(0, n_rep-1, 2))]
    if n_rep > 2: swaps += ['--swap-set', ','.join('%d-%d' % (i, i+1) for i in range(1, n_rep-1, 2))]
    cmd = [os.path.join(upside_path, obj, 'upside'), '--duration', str(duration), '--frame-interval', str(duration/200.),
           '--seed', str(seed), '--temperature', temps, '--replica-interval', '1', '--disable-recentering'] + swaps + extra
    with tempfile.TemporaryDirectory() as d:
        copies = [os.path.join(d, '%d.up' % i) for i in range(n_rep)]
        for c in copies: shutil.copyfile(config, c)
        log = sp.run(cmd + copies, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, check=True).stdout
        obs = [observables(c) for c in copies]
    names = list(obs[0])
    os.makedirs(work, exist_ok=True)
    np.savez(os.path.join(work, '%s.%d.npz' % (arm, seed)), names=names,
             means=np.array([[o[k].mean() for k in names] for o in obs]))
    speed = re.search(r'finished in ([0-9.]+) seconds \(([0-9.]+) us/systems/step', log)
    print('%s seed %d: %s s, %s us per system per step' % (arm, seed, speed.group(1), speed.group(2)))

def observables(path):
    with tb.open_file(path) as t:
        o, i = t.root.output, t.root.input
        n = len(o.pos)
        x = o.pos[n//2:, 0]
        T = o.temperature[n//2:, 0]
        obs = {'KE/1.5kT': o.kinetic[n//2:, 0]/(1.5*T), 'potential': o.potential[n//2:, 0]}
        hybrid = 'hybrid_control' in i
        if hybrid:
            obs['protein KE/1.5kT'] = o.protein_kinetic[n//2:, 0]/(1.5*T)
            obs['lipid KE/1.5kT'] = o.lipid_kinetic[n//2:, 0]/(1.5*T)
            n_res = len(i.sequence)
            lipid = i.particle_class[:] == b'LIPID'
            po4 = lipid & (i.atom_names[:] == b'PO4')
            t_frame = o.time[n//2+1] - o.time[n//2]
    if hybrid:  # the backbone N, CA, C, O of each residue first (stride 4)
        N, CA, C = x[:, 0:4*n_res:4], x[:, 1:4*n_res:4], x[:, 2:4*n_res:4]
        lag = len(x)//2
        dxy = x[lag:, lipid, :2] - x[:-lag, lipid, :2]
        obs['lipid D'] = (dxy**2).sum(-1).mean(1)/(4.*lag*t_frame)
        z = x[:, po4, 2]
        obs['thickness'] = 2.*np.abs(z - z.mean(1, keepdims=True)).mean(1)
    else:
        N, CA, C = x[:, 0::3], x[:, 1::3], x[:, 2::3]
    obs['N-CA'] = np.linalg.norm(CA-N, axis=-1).mean(1)
    obs['CA-C'] = np.linalg.norm(C-CA, axis=-1).mean(1)
    obs['C-N']  = np.linalg.norm(N[:,1:]-C[:,:-1], axis=-1).mean(1)
    obs['N-CA-C'] = angle(N, CA, C)
    obs['CA-C-N'] = angle(CA[:,:-1], C[:,:-1], N[:,1:])
    obs['C-N-CA'] = angle(C[:,:-1], N[:,1:], CA[:,1:])
    return obs

def angle(a, b, c):
    # mean over the chain of the angle at b, per frame
    u, v = a-b, c-b
    cos = (u*v).sum(-1)/(np.linalg.norm(u, axis=-1)*np.linalg.norm(v, axis=-1))
    return np.degrees(np.arccos(cos)).mean(1)

def seed_means(work, arm):
    # the arm's seeds and observables, and each seed's means by replica and observable
    files = [f for f in os.listdir(work) if f.split('.')[0] == arm and f.endswith('.npz')]
    if len(files) < 2:
        raise SystemExit('arm %s has %d seeds in %s; its standard error needs 2 or more' % (arm, len(files), work))
    files.sort(key=lambda f: int(f.split('.')[1]))
    runs = [np.load(os.path.join(work, f)) for f in files]
    return [int(f.split('.')[1]) for f in files], list(runs[0]['names']), np.array([r['means'] for r in runs])

def compare(work, arm_a, arm_b):
    seeds_a, names, a = seed_means(work, arm_a)
    seeds_b, _, b = seed_means(work, arm_b)
    print('%s: seeds %s; %s: seeds %s' % (arm_a, seeds_a, arm_b, seeds_b))
    t, p = stats.ttest_ind(b, a, axis=0, equal_var=False)
    for r in range(a.shape[1]):
        print('replica %d: ' % r + ', '.join('%s %.4g %.4g (%+.1f SE, p %.2g)' %
              (k, a[:, r, j].mean(), b[:, r, j].mean(), t[r, j], p[r, j]) for j, k in enumerate(names)))
    print('largest difference %.1f standard errors; smallest p %.2g of %d comparisons, Bonferroni %.2g' %
          (abs(t).max(), p.min(), t.size, min(1., p.min()*t.size)))

def main():
    if sys.argv[1] == 'run':
        work, arm, seed, config, n_rep, duration, obj = sys.argv[2:9]
        run(work, arm, int(seed), config, int(n_rep), float(duration), obj, sys.argv[9:])
    elif sys.argv[1] == 'compare':
        compare(*sys.argv[2:5])
    else:
        raise SystemExit('usage: device_dyn.py run|compare ...')

if __name__ == '__main__':
    main()
