#!/usr/bin/env python
# Bitwise regression record of the engine.  "dump" evaluates every node's output, the force and every parameter
# derivative on fixed frames of each config, and runs a short seeded replica-exchange simulation with pivot moves;
# "compare" lists every array that differs between two records, and apart from them the keys only one record holds.
# usage: regress.py dump WORK_DIR RECORD.npz CONFIG [CONFIG ...]
#        regress.py compare RECORD_A.npz RECORD_B.npz
import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb

upside_path = os.environ['UPSIDE_HOME']
sys.path.insert(0, upside_path + '/py')
import upside_engine as ue

def static_record(config, record):
    name = os.path.basename(config)[:-3]
    with tb.open_file(config) as t:
        x0 = t.root.input.pos[:,:,0].astype('f4')
        nodes = [g._v_name for g in t.root.input.potential]
    rng = np.random.default_rng(1)
    frames = [x0, (x0 + rng.normal(0., 0.2, x0.shape)).astype('f4'), (x0 + rng.normal(0., 1.0, x0.shape)).astype('f4')]
    engine = ue.Upside(config, quiet=True)
    for i, x in enumerate(frames):
        key = '%s/frame%d/' % (name, i)
        record[key+'deriv'] = engine.deriv(x)
        for node in nodes:
            record[key+node] = engine.get_output(node)
            n_param = engine.get_param_size(node)
            if not n_param: continue
            try:
                record[key+node+'/param_deriv'] = engine.get_param_deriv((n_param,), node)
            except RuntimeError:  # not implemented for this node
                pass
    del engine

def dynamics_record(work, config, record):
    name = os.path.basename(config)[:-3]
    copies = [os.path.join(work, '%s.run.%d.up' % (name, i)) for i in range(2)]
    for c in copies: shutil.copyfile(config, c)
    sp.check_call([upside_path + '/obj/upside', '--duration', '20', '--frame-interval', '1', '--seed', '1',
                   '--temperature', '0.80,0.85', '--replica-interval', '1', '--swap-set', '0-1',
                   '--monte-carlo-interval', '2', '--disable-recentering'] + copies, stdout=sp.DEVNULL)
    for i, c in enumerate(copies):
        with tb.open_file(c) as t:
            for node in ['pos', 'potential', 'kinetic', 'replica_index']:
                record['%s/run%d/%s' % (name, i, node)] = t.get_node('/output/' + node)[:]

def main():
    if sys.argv[1] == 'dump':
        work, out = sys.argv[2], sys.argv[3]
        os.makedirs(work, exist_ok=True)
        record = {}
        for config in sys.argv[4:]:
            static_record(config, record)
            dynamics_record(work, config, record)
        np.savez(out, **record)
        print('%d arrays in %s' % (len(record), out))
    elif sys.argv[1] == 'compare':
        a, b = np.load(sys.argv[2]), np.load(sys.argv[3])
        only_a, only_b = sorted(set(a.files) - set(b.files)), sorted(set(b.files) - set(a.files))
        for k in only_a: print('%-70s only in A' % k)
        for k in only_b: print('%-70s only in B' % k)
        both = sorted(set(a.files) & set(b.files))
        n_diff = 0
        for k in both:
            if a[k].shape != b[k].shape:
                print('%-70s shape %s vs %s' % (k, a[k].shape, b[k].shape)); n_diff += 1
            elif a[k].tobytes() != b[k].tobytes():
                d = np.abs(a[k].astype('f8') - b[k].astype('f8'))
                print('%-70s max |diff| %.3e  (%d of %d differ)' % (k, np.nanmax(d), np.sum(d != 0), d.size)); n_diff += 1
        print('%d of %d arrays in both records differ; %d only in A, %d only in B' % (n_diff, len(both), len(only_a), len(only_b)))
    else:
        raise ValueError('mode must be dump or compare')

if __name__ == '__main__':
    main()
