#!/usr/bin/env python
# V1 of the device engine (upside3.md 10.2): every node's output (a potential node's is its potential) and
# sensitivity, and each system's potential and derivative, from the device engine against the CPU engine on fixed
# frames.  Each config runs as two systems of one device engine at different frames, so the per-system offsets are
# exercised.  The frames are regress.py's: the native, and Gaussian noise of 0.2 and 1.0 A on it; with --frames
# FRAMES.npz, the npz's array (frame, atom, xyz) under the config's file name without .up.  Each frame's line gives the
# potential's difference in E_up and the derivative's error, which a BP tolerance comparison reads.
#
# Rounding: the device differs from the CPU in the rounding of every operation, and a quantity that is a small
# difference of large ones (a spring's energy near its rest length) amplifies that.  So each quantity's error, its
# RMS difference relative to the CPU's RMS, is set beside its rounding scale: the RMS, over N_MOVE moves, of the same
# measure between the CPU at the frame and at the frame with every coordinate moved by ULPS (default 4) units in the
# last place, in a random direction.  One move is a noisy estimate where a quantity is a small sum of large terms:
# chig_rama's native potential, 0.316 E_up from terms of ~5, moved by 1e-6 to 5e-5 E_up over six moves on the
# broadwl build, and by exactly 0 on the one move the test drew (10-29).  One coordinate ulp is too small a scale:
# the device's omega dihedrals differ from the CPU's by one ulp of pi, about 2.5 coordinate ulps, which moves chig's
# native omega-spring energy (0.0044 E_up, each omega 0.01 rad from rest) by 2.7e-5 relative, 1e-7 E_up (K80, 10-29;
# findings 11g).  A quantity passes if its error is at most TOL (default 1e-5) or at most ULP_FACTOR (default 10)
# times its rounding scale.
#
# A rigid frame's quaternion (affine_alignment's last 4 columns) is compared up to its sign: q and -q are one rotation,
# every consumer of the frame is even in q, and the CPU's eigensolver, which iterates four residues together, can end
# at the other sign than the device's, which iterates each residue alone.
# A frame whose CPU potential or derivative reaches 2^26, the device's fixed-point range (device.h), is past every
# physical state (dry-MARTINI beads overlapped by the 1.0 A noise: 2.7e9 E_up for 1rkl's hybrid); there a term reads
# NaN on the device by design, so that a blow-up shows.  Such a frame passes if the device's potential or derivative
# is not finite, and its quantities are not compared.
#
# With --batch the configs run together as one device engine, a system each (different proteins of one force field
# share a process), each system at the configs' frame f in turn, and each is tested against its own CPU engine; the
# configs need the same node names and the same number of frames.
#
# Prints each config's worst quantity and every failure, and exits 1 if any fails.  A config with a node that has no
# device implementation is reported as skipped (the engine names the nodes on stderr).
# usage: device_v1.py [--tol TOL] [--ulps N] [--ulp-factor F] [--frames FRAMES.npz] [--batch] CONFIG [CONFIG ...]
import sys, os
import numpy as np
import tables as tb

upside_path = os.environ['UPSIDE_HOME']
ACCUM_RANGE = 2.**26  # device.h's accum_limit
N_MOVE = 4
sys.path.insert(0, upside_path + '/py')
import upside_engine as ue

def rel_rms(ref, x):
    ref = np.asarray(ref, dtype='f8'); x = np.asarray(x, dtype='f8')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(ref)): return np.inf
    scale = np.sqrt(np.mean(ref**2))
    diff = np.sqrt(np.mean((x-ref)**2))
    return diff/scale if scale > 0. else diff

def frames(config, given):
    if given is not None:
        return list(given[os.path.basename(config)[:-len('.up')]].astype('f4'))
    with tb.open_file(config) as t:
        x0 = t.root.input.pos[:,:,0].astype('f4')
    rng = np.random.default_rng(1)
    return [x0, (x0 + rng.normal(0., 0.2, x0.shape)).astype('f4'), (x0 + rng.normal(0., 1.0, x0.shape)).astype('f4')]

def ulp_moved(x, n_ulp, rng):
    return (x + np.where(rng.random(x.shape) < 0.5, 1., -1.)*n_ulp*np.spacing(x)).astype('f4')

def same_sign_quaternion(ref, x):
    """x with each row's quaternion (columns 3-6) of ref's sign"""
    x = x.copy()
    flip = np.sum(ref[:,3:7]*x[:,3:7], axis=1) < 0.
    x[flip,3:7] *= -1.
    return x

def cpu_record(engine, nodes, x):
    deriv = engine.deriv(x)  # leaves the engine's node outputs and sensitivities at x
    rec = {'potential': [engine.energy(x)], 'deriv': deriv}
    for n in nodes:  # a potential node's output is its potential, its sensitivity 1
        rec[n+'/output'] = engine.get_output(n)
        rec[n+'/sens'] = engine.get_sens(n)
    return rec

def node_names(config):
    with tb.open_file(config) as t:
        return [g._v_name for g in t.root.input.potential]

def compare(dev, cpu, nodes, x, f, n_ulp, rng, results):
    """Each system's quantities at x (one frame per system, frame number f) from the device engine against the
    system's CPU engine, appended to results[ns] as (name, error, rounding scale)"""
    pot, deriv = dev.evaluate(x)
    for ns in range(len(cpu)):
        d = {'potential': [pot[ns]], 'deriv': deriv[ns]}
        for n in nodes:
            d[n+'/output'] = dev.get_output(ns, n)
            d[n+'/sens'] = dev.get_sens(ns, n)
        ref = cpu_record(cpu[ns], nodes, x[ns])
        moved = [cpu_record(cpu[ns], nodes, ulp_moved(x[ns], n_ulp, rng)) for i in range(N_MOVE)]
        if max(abs(ref['potential'][0]), np.abs(ref['deriv']).max()) >= ACCUM_RANGE:
            shown = not (np.isfinite(d['potential'][0]) and np.all(np.isfinite(d['deriv'])))
            print('    frame%d/system%d: CPU potential %.3g, largest derivative %.3g, past the device range; '
                  'device not finite: %s' % (f, ns, ref['potential'][0], np.abs(ref['deriv']).max(), shown))
            results[ns].append(('frame%d/system%d/past range' % (f, ns), 0. if shown else np.inf, 0.))
            continue
        print('    frame%d/system%d: potential %.6g, device - CPU %+.2e E_up; derivative error %.2e of its RMS' % (
            f, ns, ref['potential'][0], d['potential'][0] - ref['potential'][0], rel_rms(ref['deriv'], d['deriv'])))
        for n in nodes:
            if n.startswith('affine_alignment'):
                d[n+'/output'] = same_sign_quaternion(ref[n+'/output'], d[n+'/output'])
                for m in moved: m[n+'/output'] = same_sign_quaternion(ref[n+'/output'], m[n+'/output'])
        for k in ref:
            scale = np.sqrt(np.mean([rel_rms(ref[k], m[k])**2 for m in moved]))
            results[ns].append(('frame%d/system%d/%s' % (f, ns, k), rel_rms(ref[k], d[k]), scale))

def report(config, results, tol, ulp_factor):
    passed = lambda err, scale: err <= tol or err <= ulp_factor*scale
    failed = [r for r in results if not passed(r[1], r[2])]
    worst = max(results, key=lambda r: r[1]/max(tol, ulp_factor*r[2]) if np.isfinite(r[1]) else 1e300)
    print('%s: %d quantities, worst %s error %.2e (rounding scale %.2e), %d failed' % (
        os.path.basename(config), len(results), worst[0], worst[1], worst[2], len(failed)))
    for name, err, scale in failed: print('    FAIL %s error %.3e rounding scale %.3e' % (name, err, scale))
    return not failed

def check(config, tol, n_ulp, ulp_factor, given_frames):
    """A config as two systems of one device engine, at frames f and f+1"""
    nodes = node_names(config)
    cpu = [ue.Upside(config), ue.Upside(config)]
    try:
        dev = ue.DeviceUpside(cpu)
    except RuntimeError:
        print('%s: skipped, refused by the device engine' % os.path.basename(config))
        return True
    fr = frames(config, given_frames)
    rng = np.random.default_rng(2)
    results = []  # both systems' quantities in one list, in the order they are computed
    for f in range(len(fr)):
        compare(dev, cpu, nodes, [fr[f], fr[(f+1)%len(fr)]], f, n_ulp, rng, [results, results])
    del dev
    return report(config, results, tol, ulp_factor)

def check_batch(configs, tol, n_ulp, ulp_factor, given_frames):
    """The configs as one device engine, a system each, every system at its config's frame f"""
    nodes = node_names(configs[0])
    if any(node_names(c) != nodes for c in configs): raise ValueError('--batch needs configs of the same nodes')
    cpu = [ue.Upside(c) for c in configs]
    try:
        dev = ue.DeviceUpside(cpu)
    except RuntimeError:
        print('batch of %d configs: refused by the device engine' % len(configs))
        return False
    fr = [frames(c, given_frames) for c in configs]
    if len(set(map(len, fr))) != 1: raise ValueError('--batch needs the same number of frames per config')
    rng = np.random.default_rng(2)
    results = [[] for c in configs]
    for f in range(len(fr[0])):
        compare(dev, cpu, nodes, [x[f] for x in fr], f, n_ulp, rng, results)
    del dev
    return all([report(c, r, tol, ulp_factor) for c, r in zip(configs, results)])

def main():
    args = sys.argv[1:]
    tol, n_ulp, ulp_factor, given_frames, batch = 1e-5, 4, 10., None, False
    while args[0].startswith('--'):
        if args[0] == '--batch':
            batch = True
            args = args[1:]
            continue
        if args[0] == '--tol': tol = float(args[1])
        elif args[0] == '--ulps': n_ulp = int(args[1])
        elif args[0] == '--ulp-factor': ulp_factor = float(args[1])
        elif args[0] == '--frames': given_frames = np.load(args[1])
        else: raise ValueError('unknown option ' + args[0])
        args = args[2:]
    print(ue.device_description())
    if batch: ok = [check_batch(args, tol, n_ulp, ulp_factor, given_frames)]
    else:     ok = [check(c, tol, n_ulp, ulp_factor, given_frames) for c in args]
    sys.exit(0 if all(ok) else 1)

if __name__ == '__main__':
    main()
