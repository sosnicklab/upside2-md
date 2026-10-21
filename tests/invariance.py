#!/usr/bin/env python
# Statistical invariance tests of the engine's sampling.  Each prints its measurements and PASS or FAIL with the
# criterion, which was fixed before the data were seen.
#   equipartition  mean kinetic energy per degree of freedom over kT/2 for every atom class, under the engine's own
#                  thermostat: N, CA, C of a protein-only config; PROTEIN, LIPID, ION of a hybrid, which restarts
#                  from /input/mom with the production command.  Fixed atoms and z-held atoms are not counted.
#                  The hybrid must be an equilibrated production state: the first production frame after the
#                  rigid handoff carries no internal protein momentum and heats the bilayer as it relaxes.  The SE
#                  comes from the statistical inefficiency g of each class's series, pooled over the copies: under
#                  the diffusion-calibrated friction (m/gamma 35 time units) lipid kinetic energy stays correlated
#                  for about 55 frames, longer than a block of a few hundred frames can absorb.
#   kinetic        the same analysis of a trajectory already recorded with --record-momentum (read only)
#   mc             pivot Monte Carlo on a copy of a protein-only config holding only rama_map_pot, so each
#                  residue's (phi,psi) must follow exp(-E_rama/kT); the integrator is made inert by a 1e-7 time step
#   remd           temperature replica exchange of harmonically tethered atoms, whose potential energy at T is
#                  Gamma(3N/2, T): per-replica energy statistics and the exchange acceptance
# usage: invariance.py equipartition WORK_DIR CONFIG
#        invariance.py kinetic RECORD.up
#        invariance.py mc WORK_DIR CONFIG
#        invariance.py remd WORK_DIR
import sys, os, shutil, time
import subprocess as sp
import numpy as np
import tables as tb
from scipy import stats, special
from pymbar import timeseries

upside_path = os.environ['UPSIDE_HOME']
n_thread = 4

def run(configs, args):
    env = dict(os.environ, OMP_NUM_THREADS=str(n_thread))
    t0 = time.time()
    sp.check_call([upside_path + '/obj/upside'] + configs + args, stdout=sp.DEVNULL, env=env)
    print('  engine: %d systems, %.0f s' % (len(configs), time.time()-t0))

def copies(config, work, name, n):
    paths = [os.path.join(work, '%s.%d.up' % (name, i)) for i in range(n)]
    for p in paths: shutil.copyfile(config, p)
    return paths

def block_means(x, n_block):
    return np.array([y.mean(axis=0) for y in np.array_split(x, n_block)])

def mean_se(blocks):
    return blocks.mean(axis=0), blocks.std(axis=0, ddof=1)/np.sqrt(len(blocks))

def correlated_mean_se(series):
    # mean and SE of equally distributed correlated series, one per trajectory, and their pooled statistical
    # inefficiency g, the number of frames per independent sample
    g = timeseries.statistical_inefficiency_multiple(series)
    x = np.concatenate(series)
    return x.mean(), x.std()*np.sqrt(g/len(x)), g

def verdict(ok):
    return 'PASS' if ok else 'FAIL'


# ---- equipartition ----

def kinetic_classes(t):
    # atom classes and the dynamic mask; for a hybrid also the z-held PO4 atoms and the rounds they are held for
    inp = t.root.input
    n_atom = inp.pos.shape[0]
    if 'particle_class' not in inp:
        return np.array(['N', 'CA', 'C'])[np.arange(n_atom) % 3], np.ones(n_atom, bool), None
    cls = inp.particle_class[:].astype(str)
    roles = inp.atom_roles[:].astype(str)
    protein = inp.hybrid_env_topology.protein_membership[:] >= 0
    # the O and BB slots are placed by martini_hybrid_position and held fixed (up.md 5)
    dynamic = ~protein | np.isin(roles, ['N', 'CA', 'C'])
    ctrl = inp.hybrid_control._v_attrs
    hold = None
    if inp.stage_parameters._v_attrs.current_stage.decode() == 'production' and ctrl.sc_env_po4_z_hold_steps > 0:
        # the logged momentum of round start+r is held when transition_start + r <= hold_steps
        start = int(inp.mom._v_attrs.round)
        transition = int(ctrl.sc_env_transition_step_start) if 'sc_env_transition_step_start' in ctrl else 0
        hold = (~protein & (roles == 'PO4'), start + int(ctrl.sc_env_po4_z_hold_steps) - transition)
    return cls, dynamic, hold

def halfstep_excess(paths, cls, mass, names, time_step, n_frame=200):
    # The leapfrog logs the momentum half a step behind the coordinates, p = p(t) + (h/2) F + O(h^2), and
    # <p.F> = 0 in equilibrium, so its kinetic energy per degree of freedom exceeds kT/2 by (h^2/4) <F^2>/(m kT).
    sys.path.insert(0, upside_path + '/py')
    import upside_engine as ue
    f2 = []
    for path in paths:
        engine = ue.Upside(path, quiet=True)
        with tb.open_file(path) as t:
            pos = t.root.output.pos[:, 0]
            temp = t.root.output.temperature[0, 0]
        for i in np.linspace(len(pos)//10, len(pos)-1, n_frame//len(paths)).astype(int):
            f2.append(np.sum(engine.deriv(pos[i]).astype('f8')**2, axis=1)/mass/temp)
        del engine
    f2 = np.mean(f2, axis=0)
    return dict((c, 0.25*time_step**2*f2[cls == c].mean()/3.) for c in names)

def kinetic_report(paths, burn_fraction, time_step=None):
    # time_step: the leapfrog's, whose half-step kinetic excess is the expectation; None expects exactly kT/2
    with tb.open_file(paths[0]) as t:
        cls, dynamic, hold = kinetic_classes(t)
        mass = t.root.input.mass[:] if 'mass' in t.root.input else np.ones(len(cls))
    names = list(dict.fromkeys(cls))
    expect = halfstep_excess(paths, cls, mass, names, time_step) if time_step else dict((c, 0.) for c in names)
    rows, n_frame = [], 0
    for path in paths:
        with tb.open_file(path) as t:
            p = t.root.output.mom[:, 0]
            temp = t.root.output.temperature[:, 0]
            rounds = t.root.output.round[:]
        first = int(burn_fraction*len(p))
        p, temp, rounds = p[first:], temp[first:], rounds[first:]
        n_frame += len(p)
        counted = np.repeat(dynamic[None, :, None], len(p), axis=0).repeat(3, axis=2)
        if hold is not None:
            counted[np.ix_(rounds <= hold[1], hold[0], [2])] = False
        stray = np.sum(counted & (p == 0.))
        if stray: print('  %s: %d counted momentum components are exactly zero' % (path, stray))
        e2 = p.astype('f8')**2/mass[None, :, None]*counted
        row = []
        for c in names:
            s = cls == c
            row.append(e2[:, s].sum(axis=(1, 2))/(counted[:, s].sum(axis=(1, 2))*temp))
            row.extend(e2[:, s, d].sum(axis=1)/(counted[:, s, d].sum(axis=1)*temp) for d in range(3))
        rows.append(row)
    m = np.concatenate(rows, axis=1).mean(axis=1)
    print('  %d frames from %d trajectories, SE from g; criterion: FAIL if |KE/(0.5 kT) - expected| > 3 SE or > 2%%'
          % (n_frame, len(paths)))
    print('  %-8s %6s %6s %9s %9s %8s %8s %6s    %s' % ('class', 'atoms', 'dof', 'KE/0.5kT', 'expected', 'SE',
                                                       'dev/SE', 'g', 'x y z'))
    ok = True
    for i, c in enumerate(names):
        r, e, g = correlated_mean_se([row[4*i] for row in rows])
        x = 1. + expect[c]
        bad = abs(r-x) > 3.*e or abs(r-x) > 0.02
        ok &= not bad
        n_dyn = np.sum((cls == c) & dynamic)
        print('  %-8s %6d %6d %9.4f %9.4f %8.4f %8.1f %6.1f    %.4f %.4f %.4f  %s' % (c, n_dyn, 3*n_dyn, r, x, e,
              (r-x)/e, g, m[4*i+1], m[4*i+2], m[4*i+3], verdict(not bad)))
    return ok

def equipartition(work, config):
    with tb.open_file(config) as t:
        hybrid = 'brownian' in t.root.input
    paths = copies(config, work, 'equipartition', 4)
    if hybrid:
        # the production command of the workflow (Stage 7.0: MD), restarted from the fixture's momentum
        # 45000 rounds give each copy about 16 g of the lipid class
        run(paths, ['--duration-steps', '45000', '--frame-interval', '0.45', '--temperature', '0.8647',
                    '--time-step', '0.009', '--integrator', 'v', '--thermostat-timescale', '5.0',
                    '--thermostat-interval', '-1', '--seed', '1', '--disable-recentering', '--record-momentum',
                    '--restart-using-momentum'])
        return kinetic_report(paths, 0.)
    run(paths, ['--duration-steps', '2000000', '--frame-interval', '1', '--temperature', '0.85', '--seed', '1',
                '--disable-recentering', '--record-momentum'])
    return kinetic_report(paths, 0.1, time_step=0.009)


# ---- pivot Monte Carlo ----

def dihedral(a, b, c, d):
    b0, b1, b2 = a-b, c-b, d-c
    b1 = b1/np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0*b1, axis=-1, keepdims=True)*b1
    w = b2 - np.sum(b2*b1, axis=-1, keepdims=True)*b1
    return np.arctan2(np.sum(np.cross(b1, v)*w, axis=-1), np.sum(v*w, axis=-1))

def rama_energy(grid, phi, psi):
    # the engine's LayeredPeriodicSpline2D: the periodic cubic interpolating spline through the 72x72 grid, as
    # cubic B-spline coefficients from the grid's Fourier transform
    n = grid.shape[0]
    h = (4. + 2.*np.cos(2.*np.pi*np.arange(n)/n))/6.
    c = np.real(np.fft.ifft2(np.fft.fft2(grid)/(h[:, None]*h[None, :])))
    scale = n*(0.5/np.pi - 1e-7)
    x, y = (phi+np.pi)*scale, (psi+np.pi)*scale
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    def basis(f): return [(1.-f)**3/6., (3.*f**3-6.*f**2+4.)/6., (-3.*f**3+3.*f**2+3.*f+1.)/6., f**3/6.]
    bx, by = basis(x-ix), basis(y-iy)
    return sum(bx[a]*by[b]*c[(ix+a-1) % n, (iy+b-1) % n] for a in range(4) for b in range(4))

def bin_probabilities(grid, temperature, n_bin, n_sub):
    # exp(-E/T) integrated over n_bin x n_bin cells by the midpoint rule on n_sub x n_sub points per cell
    m = n_bin*n_sub
    a = -np.pi + (np.arange(m)+0.5)*2.*np.pi/m
    e = rama_energy(grid, *np.meshgrid(a, a, indexing='ij'))
    w = np.exp(-(e-e.min())/temperature).reshape(n_bin, n_sub, n_bin, n_sub).sum(axis=(1, 3))
    return w/w.sum()

def chi_square(idx, prob, n_batch):
    # goodness of fit of the bin indices idx to prob.  Bins expecting fewer than 5 counts are pooled into one cell,
    # itself merged into the smallest other cell if it still expects fewer than 5.  The statistical inefficiency g
    # is the between-batch homogeneity statistic of the same cells over its degrees of freedom.
    n = len(idx)
    keep = prob*n >= 5.
    cell = np.where(keep, np.cumsum(keep)-1, keep.sum())
    if prob[~keep].sum()*n < 5.:
        cell[~keep] = np.argmin(np.where(keep, prob, np.inf))
    k = cell.max()+1
    obs = np.bincount(cell[idx], minlength=k)
    exp = np.bincount(cell, weights=prob, minlength=k)*n
    x2 = np.sum((obs-exp)**2/exp)
    kl = np.sum(obs[obs > 0]*np.log(obs[obs > 0]/exp[obs > 0]))/n
    batches = np.array([np.bincount(cell[b], minlength=k) for b in np.array_split(idx, n_batch)])
    expect_b = batches.sum(axis=1)[:, None]*obs[None, :]/n
    seen = obs > 0
    g = np.sum((batches-expect_b)[:, seen]**2/expect_b[:, seen])/((n_batch-1)*(seen.sum()-1))
    return x2, k-1, g, kl

def mc(work, config):
    temps = [0.7, 1.3]
    n_bin, n_sub = 24, 24
    paths = copies(config, work, 'mc', len(temps))
    for path in paths:
        with tb.open_file(path, 'a') as t:
            pot = t.root.input.potential
            for node in list(pot._v_children):
                if node not in ('Dihedral_phi', 'Dihedral_psi', 'rama_coord', 'rama_map_pot'):
                    t.remove_node(pot, node, recursive=True)
    # one pivot attempt per round, a frame every 100 rounds; the 1e-7 time step and a thermostat that renews the
    # momentum each round leave the atoms in place, so only the pivot moves sample
    run(paths, ['--duration-steps', '20000000', '--time-step', '1e-7', '--monte-carlo-interval', '3e-7',
                '--frame-interval', '3e-5', '--thermostat-timescale', '3e-7', '--temperature',
                ','.join(map(str, temps)), '--seed', '1', '--disable-recentering'])
    with tb.open_file(paths[0]) as t:
        rama_id = t.root.input.potential.rama_coord.id[:]
        node = t.root.input.potential.rama_map_pot
        row = dict((int(r), i) for i, r in enumerate(node.residue_id[:]))
        maps = node.rama_pot[:][node.rama_map_id[:]]
        pivot = t.root.input.pivot_moves.pivot_atom[:]
    residues = [int(np.where((rama_id == q).all(axis=1))[0][0]) for q in pivot]
    print('  pivot residues %s; %d x %d bins of %.0f deg; criterion: FAIL if p(X2/g, dof) < 1e-4'
          % (residues, n_bin, n_bin, 360./n_bin))
    print('  %5s %4s %8s %5s %10s %6s %9s %9s %9s %12s' % ('T', 'res', 'samples', 'dof', 'X2', 'g', 'p',
          'KL', 'KL null', 'X2 at T=1.0'))
    ok = True
    for path, temp in zip(paths, temps):
        with tb.open_file(path) as t:
            pos = t.root.output.pos[:, 0].astype('f8')
            e_log = t.root.output.rama_map_potential[:]
            acc = t.root.output.pivot_stats[:].sum(axis=0)
        pos, e_log = pos[len(pos)//10:], e_log[len(e_log)//10:]
        bond = lambda x: np.linalg.norm(x[1:]-x[:-1], axis=-1)
        print('  T %.2f: pivot acceptance %.4f, largest bond length change over the run %.1e A'
              % (temp, acc[0]/acc[1], np.abs(bond(pos[-1])-bond(pos[0])).max()))
        for r in residues:
            q = rama_id[r]
            phi, psi = dihedral(*[pos[:, q[i]] for i in range(4)]), dihedral(*[pos[:, q[i+1]] for i in range(4)])
            # the reference energy must be the engine's: compare with its logged per-residue energy
            grid = maps[row[r]]
            de = np.abs(rama_energy(grid, phi, psi) - e_log[:, row[r]]).max()
            if de > 1e-3: raise RuntimeError('reference rama energy differs from the engine by %.2e' % de)
            w = 2.*np.pi/n_bin
            idx = (np.floor((phi+np.pi)/w).astype(int) % n_bin)*n_bin + np.floor((psi+np.pi)/w).astype(int) % n_bin
            prob = bin_probabilities(grid, temp, n_bin, n_sub).ravel()
            x2, dof, g, kl = chi_square(idx, prob, 20)
            p = stats.chi2.sf(x2/max(g, 1.), dof)
            x2_wrong = chi_square(idx, bin_probabilities(grid, 1.0, n_bin, n_sub).ravel(), 20)[0]
            ok &= p >= 1e-4
            print('  %5.2f %4d %8d %5d %10.1f %6.3f %9.2e %9.2e %9.2e %12.0f  %s' % (temp, r, len(idx), dof, x2, g,
                  p, kl, max(g, 1.)*dof/(2.*len(idx)), x2_wrong, verdict(p >= 1e-4)))
    return ok


# ---- replica exchange ----

def harmonic_config(path, n_atom, k):
    # every atom tethered to the origin by k/2 |x|^2, written as one Spring node per Cartesian component of pos
    with tb.open_file(path, 'w') as t:
        inp = t.create_group(t.root, 'input')
        t.create_array(inp, 'pos', obj=np.random.default_rng(1).normal(0., 1., (n_atom, 3, 1)).astype('f4'))
        pot = t.create_group(inp, 'potential')
        for d in range(3):
            g = t.create_group(pot, 'Spring_tether_%s' % 'xyz'[d])
            g._v_attrs.arguments = np.array([b'pos'])
            g._v_attrs.integrator_level = 0
            g._v_attrs.dim1 = d
            g._v_attrs.pbc = 0
            g._v_attrs.box_len = 0.
            t.create_array(g, 'id', obj=np.arange(n_atom))
            t.create_array(g, 'equil_dist', obj=np.zeros(n_atom))
            t.create_array(g, 'spring_const', obj=k*np.ones(n_atom))

def remd(work):
    n_atom, k = 10, 4.
    temps = [0.6, 0.75, 1.0, 1.4]
    a = 1.5*n_atom   # U/T is Gamma(a, 1) and so is the kinetic energy
    n_block = 20
    harmonic_config(os.path.join(work, 'harmonic.up'), n_atom, k)
    paths = copies(os.path.join(work, 'harmonic.up'), work, 'remd', len(temps))
    # exchange every 20 rounds alternating the two swap sets, a frame every 40, so each pair is attempted once
    # per frame; the thermostat (tau 5) is slow against the exchanges, so the exchanges shape the distributions
    run(paths, ['--duration-steps', '12000000', '--time-step', '0.009', '--frame-interval', '1.08',
                '--replica-interval', '0.54', '--swap-set', '0-1,2-3', '--swap-set', '1-2', '--thermostat-timescale',
                '5', '--temperature', ','.join(map(str, temps)), '--seed', '1', '--disable-recentering'])
    pot, kin, swaps = [], [], []
    for path in paths:
        with tb.open_file(path) as t:
            o = t.root.output
            pot.append(o.potential[:, 0])
            kin.append(o.kinetic[:, 0])
            swaps.append(dict(zip(o.replica_swap_partner[:], np.moveaxis(o.replica_cumulative_swaps[:], 1, 0))))
    first = len(pot[0])//10
    pot, kin = np.array(pot)[:, first:], np.array(kin)[:, first:]
    print('  %d atoms, k %.1f, %d frames per replica after 10%% burn-in, %d blocks; criterion: FAIL if |dev| > 4 SE'
          % (n_atom, k, pot.shape[1], n_block))
    print('  %5s %18s %18s %18s' % ('T', '<U>/(aT)', 'var U/(aT^2)', '<KE>/(1.5T)'))
    ok = True
    for i, temp in enumerate(temps):
        u = pot[i]/(a*temp)
        cols = [mean_se(block_means(u, n_block)),
                mean_se(np.array([a*np.var(b) for b in np.array_split(u, n_block)])),
                mean_se(block_means(kin[i]/(1.5*temp), n_block))]
        line = '  %5.2f' % temp
        for m, se in cols:
            ok &= abs(m-1.) <= 4.*se
            line += ' %8.4f +- %.4f%s' % (m, se, '' if abs(m-1.) <= 4.*se else '*')
        print(line)
    print('  %9s %22s %22s %9s  %s' % ('pair', 'measured acceptance', 'from sampled U', 'analytic', 'attempts'))
    for i in range(len(temps)-1):
        c = swaps[i][i+1]
        hits = np.diff(c[:, 0])[first:]
        tries = np.diff(c[:, 1])[first:]
        if not np.all(tries == 1): raise RuntimeError('expected one attempt per pair per frame')
        meas, meas_se = mean_se(block_means(hits.astype('f8'), n_block))
        db = 1./temps[i] - 1./temps[i+1]
        shifts = range(1, 400, 8)
        pred = np.array([np.mean([np.exp(np.minimum(0., db*(ui-np.roll(uj, s)))).mean() for s in shifts])
                         for ui, uj in zip(np.array_split(pot[i], n_block), np.array_split(pot[i+1], n_block))])
        pred, pred_se = mean_se(pred)
        exact = 2.*special.betainc(a, a, temps[i]/(temps[i]+temps[i+1]))
        good = abs(meas-pred) <= 4.*np.hypot(meas_se, pred_se) and abs(meas-exact) <= 4.*meas_se
        ok &= good
        print('  %4.2f-%4.2f %12.4f +- %.4f %12.4f +- %.4f %9.4f  %8d  %s' % (temps[i], temps[i+1], meas, meas_se,
              pred, pred_se, exact, len(hits), verdict(good)))
    return ok


def main():
    t0 = time.time()
    test = sys.argv[1]
    if test == 'kinetic':
        ok = kinetic_report([sys.argv[2]], 0.)
    else:
        work = sys.argv[2]
        os.makedirs(work, exist_ok=True)
        if   test == 'equipartition': ok = equipartition(work, sys.argv[3])
        elif test == 'mc':            ok = mc(work, sys.argv[3])
        elif test == 'remd':          ok = remd(work)
        else: raise ValueError('test must be equipartition, kinetic, mc or remd')
    print('%s %s (%.0f s)' % (test, verdict(ok), time.time()-t0))

if __name__ == '__main__':
    main()
