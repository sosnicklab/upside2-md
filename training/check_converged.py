"""Has a ConDiv run updated every file, and has every parameter group plateaued?

    python3 check_converged.py <run_dir>

Three reports, for the FF2 dual-target trainer (training/ConDiv.py), read with the run's own copy
of ConDiv.py so the field order is the one that wrote the files:

1. FILES. The newest step's written parameter files against init_param, dataset by dataset. Every
   trained file must differ; a byte-identical one means that group is not being trained.

2. FIXED POINT. The raw per-step gradient of each trained group, recovered from the Adam state:
   grad1_t = b1*grad1_{t-1} + (1-b1)*g_t, so g_t = (grad1_t - b1*grad1_{t-1}) / (1-b1). Adam's
   step is scale-invariant, so "the parameters barely moved" is NOT the signal; the gradient is.
   At a fixed point the expected gradient is zero and steps are pure noise:
     vectors: ||mean g|| / mean|g| near 1/sqrt(n), pairwise cosine near 0;
     scalars: mean g indistinguishable from 0, |t| under about 2.

3. PLATEAU. Parameter drift over the second half of the steps against the first half. A group
   that has settled drifts less late than early, or no more than its step-to-step jitter.

Judge nothing on fewer than one epoch of steps: the glycine map has stalled for ~20 steps and
then moved again, and short windows have been misread as plateaus before.
"""
import glob
import os
import pickle
import sys

import numpy as np
import tables as tb

RUN = os.path.abspath(sys.argv[1])
OUT = os.path.join(RUN, 'run_output')
sys.path.insert(0, OUT)
import ConDiv  # noqa: E402  the run's own copy
import __main__  # noqa: E402
for _n in dir(ConDiv):
    if _n[0].isupper():
        setattr(__main__, _n, getattr(ConDiv, _n))

FIELDS = ConDiv.FIELDS
steps = sorted(d for d in glob.glob(os.path.join(OUT, 'epoch_*_minibatch_*'))
               if os.path.exists(os.path.join(d, 'checkpoint.pkl')))
if not steps:
    sys.exit('no completed step under %s yet' % OUT)
state = pickle.load(open(os.path.join(steps[-1], 'checkpoint.pkl'), 'rb'))
B1 = state['solver'].beta1
alpha = state['initial_alpha']
trained = [f for i, f in enumerate(FIELDS)
           if np.any(np.asarray(alpha[i]) != 0) and getattr(state['param'], f) is not None]
print('%s: %d steps, %d minibatches per epoch; trained groups: %s\n'
      % (RUN, len(steps), len(state['minibatches']), ' '.join(trained)))

# --- 1. files --------------------------------------------------------------------------------
init = state['init_param_files']
last = steps[-1]
print('1. FILES (%s vs init_param)' % os.path.basename(last))
for key, src in sorted(init.items()):
    if key == 'rama':
        continue
    new = os.path.join(last, os.path.basename(src))
    if src.endswith('.h5'):
        with tb.open_file(src) as a, tb.open_file(new) as b:
            for node in a.walk_nodes('/', 'Array'):
                x, y = np.asarray(node[:]), np.asarray(b.get_node(node._v_pathname)[:])
                if x.dtype.kind in 'fi':
                    d = np.abs(np.asarray(y, float) - np.asarray(x, float)).max()
                    print('   %-15s %-28s %s' % (os.path.basename(src), node._v_pathname,
                                                 'max|change| %.3g' % d if d else 'UNCHANGED'))
    else:
        d = np.abs(np.loadtxt(new) - np.loadtxt(src)).max()
        print('   %-15s %-28s %s' % (os.path.basename(src), '', 'max|change| %.3g' % d
                                     if d else 'UNCHANGED'))
print()

# --- 2. fixed point ---------------------------------------------------------------------------
prev = {f: 0. for f in trained}
grads = {f: [] for f in trained}
params = {f: [] for f in trained}
for d in steps:
    st = pickle.load(open(os.path.join(d, 'solver_state.pkl'), 'rb'), encoding='latin1')
    ck = pickle.load(open(os.path.join(d, 'checkpoint.pkl'), 'rb'))
    for f in trained:
        g1 = np.asarray(st['grad1'][FIELDS.index(f)], dtype=float)
        grads[f].append(np.atleast_1d((g1 - B1 * np.asarray(prev[f])) / (1. - B1)).ravel())
        prev[f] = g1
        params[f].append(np.atleast_1d(np.asarray(getattr(ck['param'], f), float)).ravel())

n = len(steps)
print('2. FIXED POINT (n = %d; pure noise gives ratio %.3f, cosine 0)' % (n, 1. / np.sqrt(n)))
print('   %-7s %6s  %-44s' % ('group', 'size', 'statistic'))
for f in trained:
    G = np.array(grads[f])
    if G.shape[1] == 1:
        v = G[:, 0]
        t = v.mean() / (v.std(ddof=1) / np.sqrt(n)) if n > 2 and v.std() > 0 else np.nan
        print('   %-7s %6d  mean %+.4g  sd %.4g  t %+.2f' % (f, 1, v.mean(), v.std(ddof=1), t))
        continue
    norms = np.linalg.norm(G, axis=1)
    ratio = np.linalg.norm(G.mean(0)) / norms.mean()
    iu = np.triu_indices(n, 1)
    cos = (G @ G.T / np.outer(norms, norms))[iu] if n > 1 else np.array([np.nan])
    t = cos.mean() / (cos.std(ddof=1) / np.sqrt(len(cos))) if len(cos) > 2 else np.nan
    print('   %-7s %6d  ratio %.3f   cosine %+.3f (t %+.2f)' % (f, G.shape[1], ratio,
                                                                 np.nanmean(cos), t))
print()

# --- 3. plateau ------------------------------------------------------------------------------
print('3. PLATEAU (rms parameter drift: first half, second half, step-to-step jitter)')
h = n // 2
for f in trained:
    P = np.array(params[f])
    if n < 4:
        print('   %-7s needs at least 4 steps' % f)
        continue
    early = np.sqrt(((P[h] - P[0]) ** 2).mean())
    late = np.sqrt(((P[-1] - P[h]) ** 2).mean())
    jitter = np.sqrt((np.diff(P, axis=0) ** 2).mean())
    # a random walk of step size `jitter` drifts about jitter*sqrt(steps) with no systematic pull
    walk = jitter * np.sqrt(n - 1 - h)
    print('   %-7s %.3g   %.3g   %.3g%s' % (f, early, late, jitter,
          '   still drifting' if late > early and late > 2 * walk else ''))
