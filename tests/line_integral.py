"""Line integral of a parameter derivative against the energy, at the exact frames node_fd.py drew: the test
that settles a direction node_fd flags.

Usage: line_integral.py WORK_DIR CONFIG [CONFIG ...] -- CASE [CASE ...], with node_fd's own arguments before --
CASE is config_name:term:frame:node:param, e.g. 1qhj_memb:rotamer:1:placement_fixed_point_vector_only:459.
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
    name, term, frame, node, k = case.split(':')
    path, fr = frames[(name, term)]
    node_fd.single_term_config([a for a in sys.argv[2:] if os.path.basename(a)[:-3] == name][0], term, path)
    eng = ue.Upside(path, quiet=True)
    x = fr[int(frame)].astype('f4')
    n = eng.get_param_size(node)
    p0 = eng.get_param((n,), node).astype('f8')
    k = int(k)
    w = min(abs(p0[k]), 1.) + 1e-2
    s = np.linspace(-0.4 * w, 0.4 * w, 801)
    E, A, S = [], [], []
    for si in s:
        p = p0.copy()
        p[k] += si
        pf = p.astype('f4')
        eng.set_param(pf, node)
        E.append(float(eng.energy(x)))
        A.append(float(eng.get_param_deriv((n,), node)[k]))
        S.append(float(pf[k]) - p0[k])
    eng.set_param(p0.astype('f4'), node)
    E, A, S = map(np.array, (E, A, S))
    integral = np.concatenate([[0.], np.cumsum(0.5 * (A[1:] + A[:-1]) * np.diff(S))])
    r = (E - E[400]) - (integral - integral[400])
    print(case)
    for win in (0.05, 0.2, 1.0):
        m = np.abs(s) <= win * 0.4 * w + 1e-12
        slope = np.polyfit(S[m], r[m], 1)[0]
        print('   |s| <= %.3f: drift slope % .2e (analytic at 0 % .3e, relative %.1e), max|r| %.1e' % (
            win * 0.4 * w, slope, A[400], abs(slope) / (abs(A[400]) + 1e-30), np.abs(r[m]).max()))
