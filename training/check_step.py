"""Health of a ConDiv run's newest finished step, from physical observables, not exit codes.

    python3 check_step.py <run_dir> [<epoch_xx_minibatch_yy>]

Per protein of the step: the replicas' avg_kinetic_energy/1.5kT (1.0 when the thermostat holds the
temperature the replica was given), their final radius of gyration, the restrained and free
C-alpha RMSD, whether the unfolded-state target was found, and whether every recorded gradient and
basin population is finite. Then the trained parameters against the initial ones, and the glycine
readout of plan.md Phase 8: free against restrained basin populations of glycines that are helical
or left-handed in the restrained replica, for this step and for the epoch so far.
"""

import glob
import os
import pickle as cp
import re
import sys

import numpy as np

RUN = os.path.abspath(sys.argv[1])
OUT = os.path.join(RUN, 'run_output')
sys.path.insert(0, OUT)
import ConDiv  # noqa: E402  the run's own copy
import __main__  # noqa: E402
for _n in dir(ConDiv):
    if _n[0].isupper():
        setattr(__main__, _n, getattr(ConDiv, _n))
import upside_config as uc  # noqa: E402

done = sorted(d for d in glob.glob(os.path.join(OUT, 'epoch_*_minibatch_*'))
              if os.path.exists(os.path.join(d, 'checkpoint.pkl')))
if len(sys.argv) > 2:
    step = os.path.join(OUT, sys.argv[2])
elif done:
    step = done[-1]
else:
    sys.exit('no finished step yet')
state = cp.load(open(os.path.join(step, 'checkpoint.pkl'), 'rb'))
init = cp.load(open(os.path.join(OUT, 'initial_checkpoint.pkl'), 'rb'))
seqs = {nm: list(uc.read_fasta(open(t.fasta))) for mb in state['minibatches'] for nm, t in mb}


def kinetic(path):
    for line in open(path):
        if line.startswith('avg_kinetic_energy/1.5kT'):
            return np.array(line.split()[1:], dtype=float)
    return None


def final_rg(path):
    rg = {}
    pat = re.compile(r'^\s*\d+ / \d+ elapsed\s+(\d+) .*Rg\s+([-\d.naN]+) A')
    for line in open(path):
        m = pat.match(line)
        if m:
            rg[int(m.group(1))] = float(m.group(2))
    return np.array([rg[k] for k in sorted(rg)]) if rg else None


def glycine_pops(steps):
    """(class, free, native) basin populations of glycines with non-GLY flanks, by restrained basin."""
    rows = []
    for d in steps:
        for f in glob.glob(os.path.join(d, '*.divergence.pkl')):
            code = os.path.basename(f).split('.')[0]
            s, dv = seqs[code], cp.load(open(f, 'rb'))
            fr, na = np.asarray(dv['rama_free']), np.asarray(dv['rama_native'])
            for i in range(1, len(s) - 1):
                if s[i] != 'GLY' or 'GLY' in (s[i - 1], s[i + 1]):
                    continue
                cls = 'helical' if na[i, 0] > .5 else 'left' if na[i, 1] > .5 else None
                if cls:
                    rows.append((cls, fr[i, :2], na[i, :2]))
    return rows


print(f'{RUN}: step {os.path.basename(step)}, {len(done)} steps finished\n')
bad, ke, rgs, rmsd = [], [], [], []
for f in sorted(glob.glob(os.path.join(step, '*.divergence.pkl'))):
    code = os.path.basename(f).split('.')[0]
    dv = cp.load(open(f, 'rb'))
    finite = all(np.isfinite(np.asarray(x, float)).all() for x in dv['contrast'] if x is not None) \
        and np.isfinite(dv['rama_free']).all() and np.isfinite(dv['rama_native']).all()
    log = os.path.join(step, code + '.run.0.up.output')
    k = kinetic(log) if os.path.exists(log) else None
    r = final_rg(log) if os.path.exists(log) else None
    if not finite or k is None or not np.isfinite(k).all():
        bad.append(code)
    if k is not None:
        ke.append(k)
    if r is not None:
        rgs.append((code, r))
    rmsd.append((dv['rmsd_restrain'], dv['rmsd'], dv['has_dse'], dv['walltime']))
n_div = len(rmsd)
r = np.array([x[:2] for x in rmsd])
print(f'proteins with results {n_div} of {len(state["minibatches"][0])}; non-finite or missing: '
      f'{", ".join(bad) or "none"}')
print(f'RMSD restrained median {np.median(r[:, 0]):.2f} max {r[:, 0].max():.2f} A; free median '
      f'{np.median(r[:, 1]):.2f} max {r[:, 1].max():.2f} A; unfolded-state target in '
      f'{sum(x[2] for x in rmsd)} of {n_div}; worker wall median {np.median([x[3] for x in rmsd]):.0f} s')
if ke:
    k = np.array(ke)
    print(f'avg_kinetic_energy/1.5kT over all replicas: min {k.min():.3f} median {np.median(k):.3f} '
          f'max {k.max():.3f}  (restrained replica {np.median(k[:, 0]):.3f}, SARW {np.median(k[:, -1]):.3f})')
if rgs:
    cold = np.array([x[1][1] for x in rgs])
    hot = np.array([x[1][-2] for x in rgs])
    print(f'final Rg, coldest free replica median {np.median(cold):.1f} A, hottest free {np.median(hot):.1f} A')

p0, p = init['param'], state['param']
print(f'\nhb E_alpha E_beta E_other {np.array2string(p.hb[:3], precision=4)} (start '
      f'{np.array2string(p0.hb[:3], precision=4)}); margin E_other - E_alpha '
      f'{p.hb[2] - p.hb[0]:+.4f} (start {p0.hb[2] - p0.hb[0]:+.4f})')
print(f'dhb {p.dhb[0]:.4f} (start {p0.dhb[0]:.4f}); sheet mean {np.mean(p.sheet):.4f} (start '
      f'{np.mean(p0.sheet):.4f}); bb env scale {p.bbenve:.4f} (start {p0.bbenve:.4f})')
print(f'env scale change max {np.abs(p.enve - p0.enve).max():.4f}; rot change rms '
      f'{np.sqrt(np.mean((np.asarray(p.rot) - np.asarray(p0.rot)) ** 2)):.4f}')

epoch = os.path.basename(step).split('_minibatch_')[0]
this_epoch = [d for d in done if os.path.basename(d).startswith(epoch)]
for label, steps in (('this step', [step]), (f'{epoch} so far ({len(this_epoch)} steps)', this_epoch)):
    rows = glycine_pops(steps)
    print(f'\nglycines (non-GLY flanks), {label}:  free / restrained')
    for cls in ('helical', 'left'):
        a = [x for x in rows if x[0] == cls]
        if a:
            fr, na = np.array([x[1] for x in a]), np.array([x[2] for x in a])
            print(f'  {cls:8s} n {len(a):4d}  alpha_R {fr[:, 0].mean():.3f} / {na[:, 0].mean():.3f}'
                  f'   alpha_L {fr[:, 1].mean():.3f} / {na[:, 1].mean():.3f}')
