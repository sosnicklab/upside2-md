"""Line integral of a parameter derivative against the energy, at the exact frames node_fd.py drew: the test
that settles a direction node_fd flags.

Usage: line_integral.py WORK_DIR CONFIG [CONFIG ...] -- CASE [CASE ...], with node_fd's own arguments before --
CASE is config_name:term:frame:node:param[:width], e.g. 1qhj_memb:rotamer:1:placement_fixed_point_vector_only:459;
width replaces node_fd's step scale min(|p|, 1) + 0.01 where the energy change over it is below float32 resolution.
node_fd runs unchanged in-process so its random stream, hence its frames, are reproduced; its check_term is
wrapped to keep the frames. Then r(s) = E(s) - E(0) - integral of the analytic derivative over the realised
float32 path, on 801 points over node_fd's widest stencil. A wrong derivative leaves a drift; a correct one noise.
"""
import os
import sys

import numpy as np

sys.argv, cases = sys.argv[:sys.argv.index('--')], sys.argv[sys.argv.index('--') + 1:]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import node_fd

frames = {}
orig = node_fd.check_term


def keep(config, term, path, fr):
    frames[(os.path.basename(config)[:-3], term)] = (path, [x.copy() for x in fr])
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        orig(config, term, path, fr)


node_fd.check_term = keep
node_fd.main()
ue = node_fd.ue
for case in cases:
    name, term, frame, node, k, *width = case.split(':')
    path, fr = frames[(name, term)]
    nodes = node_fd.single_term_config([a for a in sys.argv[2:] if os.path.basename(a)[:-3] == name][0], term, path)
    eng = ue.Upside(path)
    x = fr[int(frame)].astype('f4')
    n = eng.get_param_size(node)
    p0 = eng.get_param((n,), node).astype('f8')
    k = int(k)
    # the direction node_fd's param k checks: tied as the node ties it, so the engine accepts every step
    d, scale = node_fd.tie(node, nodes[node], np.eye(1, n, k)[0])[0], node_fd.tie(node, nodes[node], p0)[1]
    w = float(width[0]) if width else min(abs(p0[k]), 1.) + 1e-2
    s = np.linspace(-0.4 * w, 0.4 * w, 801)
    E, A, S = [], [], []
    for si in s:
        pf = (p0 + si * d).astype('f4')
        eng.set_param(pf, node)
        E.append(float(eng.energy(x)))
        A.append(scale * eng.get_param_deriv((n,), node).astype('f8'))
        S.append(pf.astype('f8') - p0)
    eng.set_param(p0.astype('f4'), node)
    E, A, S = map(np.array, (E, A, S))
    integral = np.concatenate([[0.], np.cumsum(np.sum(0.5 * (A[1:] + A[:-1]) * np.diff(S, axis=0), axis=1))])
    r = (E - E[400]) - (integral - integral[400])
    t = S @ d / (d @ d)   # the realised step along the direction
    a0 = A[400] @ d / np.sqrt(d @ d)
    print(case)
    for win in (0.05, 0.2, 1.0):
        m = np.abs(s) <= win * 0.4 * w + 1e-12
        slope = np.polyfit(t[m], r[m], 1)[0] / np.sqrt(d @ d)
        print('   |s| <= %.3f: drift slope % .2e (analytic at 0 % .3e, relative %.1e), max|r| %.1e' % (
            win * 0.4 * w, slope, a0, abs(slope) / (abs(a0) + 1e-30), np.abs(r[m]).max()))
