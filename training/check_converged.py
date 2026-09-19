"""Is ff_2.1 a fixed point of the strict-modernization trainer?

Adam is scale-invariant on its first step (update = -alpha*g/(|g|+eps)), so "the parameters
barely moved" is NOT the signal even at a perfect optimum.  The signal is the RAW gradient,
recovered from the Adam accumulators in solver_state.pkl:

    grad1_t = b1*grad1_{t-1} + (1-b1)*g_t   ->   g_t = (grad1_t - b1*grad1_{t-1})/(1-b1)

Update field order is (env, cov, rot, hyd, hb, sheet).  cov and hyd are folded into rot.

VECTOR parameters (env, rot): at a true optimum the expected gradient is zero, so gradients from
different minibatches are pure noise -> pairwise cosine ~ 0 and ||mean g||/mean|g| ~ 1/sqrt(n).

SCALAR parameters (hb, sheet): the test is sharper -- is the mean gradient distinguishable from
zero?  Report t = mean/(sd/sqrt(n)).  |t| under about 2 is consistent with a fixed point.

CAVEAT when reading hb and sheet: ff_2.1s 12-entry hbond.h5 and 20-value sheet file were NOT
produced by the original trainer (whose outputs were the single scalars -2.112 and -0.268); they
came from the node rewrite.  A nonzero hb/sheet gradient therefore means those hand-set values are
not at a ConDiv optimum, which is NOT the same as the port being wrong.  rot and env are the clean
port-fidelity test.
"""
import glob, os, pickle, sys
import numpy as np

RUN = sys.argv[1] if len(sys.argv) > 1 else "run_output"
B1 = 0.8
IDX = {0: "env", 2: "rot", 4: "hb", 5: "sheet"}

dirs = sorted(glob.glob(os.path.join(RUN, "epoch_*_minibatch_*")))
dirs = [d for d in dirs if os.path.exists(os.path.join(d, "solver_state.pkl"))]
if not dirs:
    sys.exit("no solver_state.pkl under %s yet" % RUN)

prev = {i: 0.0 for i in IDX}
grads = {i: [] for i in IDX}
print("steps found: %d\n" % len(dirs))
for d in dirs:
    with open(os.path.join(d, "solver_state.pkl"), "rb") as f:
        st = pickle.load(f, encoding="latin1")
    row = []
    for i in IDX:
        g1 = np.asarray(st["grad1"][i], dtype=float)
        g = (g1 - B1 * np.asarray(prev[i], dtype=float)) / (1.0 - B1)
        grads[i].append(np.atleast_1d(g).ravel())
        prev[i] = g1
        row.append("%s=%.4g" % (IDX[i], np.linalg.norm(grads[i][-1])))
    print("  %-26s %s" % (os.path.basename(d), "  ".join(row)))
print()

for i, name in IDX.items():
    G = np.array(grads[i]); n = len(G)
    if G.shape[1] == 1:                       # scalar parameter
        v = G[:, 0]
        sd = v.std(ddof=1) if n > 1 else float("nan")
        t = v.mean() / (sd / np.sqrt(n)) if n > 1 and sd > 0 else float("nan")
        print("%s (scalar):  n=%d  mean %+0.5g  sd %0.5g" % (name, n, v.mean(), sd))
        print("   t vs zero                  %+0.2f   (|t| < ~2 is consistent with a fixed point)\n" % t)
        continue
    norms = np.linalg.norm(G, axis=1)
    ratio = np.linalg.norm(G.mean(axis=0)) / norms.mean()
    cos = [float(G[a] @ G[b] / (np.linalg.norm(G[a]) * np.linalg.norm(G[b])))
           for a in range(n) for b in range(a + 1, n)
           if np.linalg.norm(G[a]) > 0 and np.linalg.norm(G[b]) > 0]
    cos = np.array(cos) if cos else np.array([np.nan])
    print("%s (vector):  n=%d" % (name, n))
    print("   mean |g|                    %.4g   (spread %.4g)" % (norms.mean(), norms.std()))
    print("   ||mean g|| / mean|g|        %.3f   (pure noise ~ %.3f = 1/sqrt(n))" % (ratio, 1/np.sqrt(n)))
    print("   pairwise cosine             mean %+.3f  sd %.3f  n_pairs %d" % (np.nanmean(cos), np.nanstd(cos), len(cos)))
    if len(cos) > 1:
        print("   cosine t-stat vs 0          %+.2f" % (np.nanmean(cos) / (np.nanstd(cos, ddof=1) / np.sqrt(len(cos)))))
    print()
