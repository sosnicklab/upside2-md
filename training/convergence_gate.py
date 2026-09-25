"""Is a ConDiv run at a fixed point? The automatic decision between validating and training on.

    python3 convergence_gate.py <run_dir>

Exit 0 if every trained group is consistent with a fixed point over the last full epoch, 3 if not
(including when there is less than one epoch to judge). Any other code is a failure of the gate
itself, and must stop whatever acts on its answer rather than read as "not converged".

WINDOW. The last full epoch: every protein has contributed exactly once, so the window-mean
gradient is the full-training-set gradient. Nothing shorter is judged; short windows have been
misread as plateaus before.

TEST. For each trained group, the raw per-step gradient g_t is recovered from the Adam state
(grad1_t = b1*grad1_{t-1} + (1-b1)*g_t). At a fixed point each step's gradient is noise symmetric
about zero, so the statistic ||sum_t g_t|| should be unremarkable among all 2^n sign flips of the
steps. The permutation p-value is exact: ||sum_t s_t g_t||^2 = s^T K s with K the n x n Gram
matrix of the steps, so all 2^(n-1) patterns (s and -s are the same) are enumerated without ever
forming a sum over the parameters, which keeps the 217,728-value glycine row cheap.

DECISION. Every group's p above 0.05 / (number of groups): a family-wise 5% chance of failing a
force field that is truly at a fixed point. Adam's steps are scale-invariant, so parameter movement
is not used at all.
"""
import glob
import itertools
import os
import pickle
import sys

import numpy as np

ALPHA = 0.05

RUN = os.path.abspath(sys.argv[1])
OUT = os.path.join(RUN, 'run_output')
sys.path.insert(0, OUT)
import ConDiv  # noqa: E402  the run's own copy
import __main__  # noqa: E402
for _n in dir(ConDiv):
    if _n[0].isupper():
        setattr(__main__, _n, getattr(ConDiv, _n))

steps = sorted(d for d in glob.glob(os.path.join(OUT, 'epoch_*_minibatch_*'))
               if os.path.exists(os.path.join(d, 'checkpoint.pkl')))
state = pickle.load(open(os.path.join(steps[-1], 'checkpoint.pkl'), 'rb'))
n = len(state['minibatches'])
if len(steps) < n:
    print('only %d steps, fewer than one epoch (%d); nothing to judge yet' % (len(steps), n))
    sys.exit(3)

B1 = state['solver'].beta1
FIELDS = ConDiv.FIELDS
trained = [f for i, f in enumerate(FIELDS) if np.any(np.asarray(state['initial_alpha'][i]) != 0)
           and getattr(state['param'], f) is not None]


def grad1(d):
    st = pickle.load(open(os.path.join(d, 'solver_state.pkl'), 'rb'), encoding='latin1')
    return [np.asarray(st['grad1'][FIELDS.index(f)], dtype=float) for f in trained]


window = steps[-n:]
prev = grad1(steps[-n - 1]) if len(steps) > n else [0. for _ in trained]
G = {f: [] for f in trained}
for d in window:
    g1 = grad1(d)
    for k, f in enumerate(trained):
        G[f].append(np.atleast_1d((g1[k] - B1 * prev[k]) / (1. - B1)).ravel())
    prev = g1

signs = np.array([(1,) + s for s in itertools.product((1, -1), repeat=n - 1)], dtype=float)
cut = ALPHA / len(trained)
print('%s: steps %d-%d (last epoch), %d groups, pass if every p > %.4f\n'
      % (RUN, len(steps) - n + 1, len(steps), len(trained), cut))
print('   %-7s %7s  %10s  %8s' % ('group', 'size', 'p', 'verdict'))
ok = True
for f in trained:
    X = np.array(G[f])
    K = X @ X.T
    stat = np.einsum('pi,ij,pj->p', signs, K, signs)
    obs = K.sum()
    p = (stat >= obs * (1 - 1e-12)).mean()
    passed = p > cut
    ok &= passed
    print('   %-7s %7d  %10.4f  %8s' % (f, X.shape[1], p, 'ok' if passed else 'PULLED'))
print('\n' + ('CONVERGED: every group is consistent with a fixed point.' if ok else
              'NOT CONVERGED: at least one group still has a systematic pull.'))
sys.exit(0 if ok else 3)
