#!/usr/bin/env python
"""Temperature-REMD analysis for Upside replica trajectories.

Upside runs coordinate-swap REMD: one .up file per slot, each at a FIXED temperature
(output/temperature[0,0], in Upside units so kT == that value in E_up). This tool reads all
replica files, reports per-temperature observables, and reconstructs the free energy F(CV) at a
target temperature by MBAR-reweighting the pooled samples (falls back to a per-temperature
histogram if pymbar is unavailable).

Usage:
  python analyze_remd.py 'outputs/sim/prot.run.*.up' out_base --cv-atoms 12 87 [--target-kT 0.8647]
  python analyze_remd.py 'outputs/sim/prot.run.*.up' out_base --cv rg      [--target-kT 0.8647]

CV options: --cv-atoms i j (distance between two atoms, Angstrom) or --cv rg (CA radius of gyration).
Outputs out_base_remd_obs.csv (T, <E>, <CV>), out_base_remd_fes.npy (cv, F), out_base_remd.png.
"""
import argparse
import glob
import numpy as np
import tables as tb


def _output_groups(t):
    i = 0
    while 'output_previous_%i' % i in t.root:
        yield t.get_node('/output_previous_%i' % i)
        i += 1
    if 'output' in t.root:
        yield t.get_node('/output')


def _read_replica(path, cv_mode, cv_atoms):
    """Return (kT, energy[n], cv[n]) for one replica file."""
    ener, pos_parts = [], []
    with tb.open_file(path) as t:
        kT = float(t.root.output.temperature[0, 0])
        for k, g in enumerate(_output_groups(t)):
            e = g.potential[:, 0]
            p = g.pos[:, 0]                       # (n, n_atom, 3)
            ener.append(e if k == 0 else e[1:])
            pos_parts.append(p if k == 0 else p[1:])
    energy = np.concatenate(ener)
    pos = np.concatenate(pos_parts)
    if cv_mode == 'dist':
        i, j = cv_atoms
        cv = np.linalg.norm(pos[:, i] - pos[:, j], axis=1)
    elif cv_mode == 'rg':
        c = pos - pos.mean(axis=1, keepdims=True)
        cv = np.sqrt((c ** 2).sum(axis=2).mean(axis=1))
    else:
        raise ValueError('unknown cv mode ' + cv_mode)
    return kT, energy, cv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('pattern', help='glob for per-replica .up files (quote it)')
    p.add_argument('output_base')
    p.add_argument('--cv', choices=['rg', 'dist'], default='rg')
    p.add_argument('--cv-atoms', type=int, nargs=2, default=None, help='two atom indices for --cv dist')
    p.add_argument('--target-kT', type=float, default=None, help='target kT (E_up) for the FES; default = lowest replica')
    p.add_argument('--bins', type=int, default=40)
    args = p.parse_args()

    cv_mode = 'dist' if args.cv == 'dist' else 'rg'
    if cv_mode == 'dist' and not args.cv_atoms:
        p.error('--cv dist requires --cv-atoms i j')

    files = sorted(glob.glob(args.pattern))
    if not files:
        p.error('no files match %s' % args.pattern)

    kTs, energies, cvs = [], [], []
    for f in files:
        kT, e, cv = _read_replica(f, cv_mode, args.cv_atoms)
        kTs.append(kT); energies.append(e); cvs.append(cv)
    kTs = np.array(kTs)
    order = np.argsort(kTs)
    kTs, energies, cvs = kTs[order], [energies[i] for i in order], [cvs[i] for i in order]

    with open('%s_remd_obs.csv' % args.output_base, 'w') as fh:
        fh.write('kT_Eup,mean_energy,mean_cv,frames\n')
        for kT, e, cv in zip(kTs, energies, cvs):
            fh.write('%.5f,%.4f,%.4f,%d\n' % (kT, e.mean(), cv.mean(), len(e)))
    print('REMD: %d replicas, kT %.4f .. %.4f' % (len(files), kTs.min(), kTs.max()))
    for kT, e, cv in zip(kTs, energies, cvs):
        print('  kT=%.4f  <E>=%.2f  <CV>=%.3f  (n=%d)' % (kT, e.mean(), cv.mean(), len(e)))

    target_kT = args.target_kT if args.target_kT else kTs.min()
    cv_all = np.concatenate(cvs)
    edges = np.linspace(cv_all.min(), cv_all.max(), args.bins + 1)
    ctr = 0.5 * (edges[:-1] + edges[1:])

    fes = None
    try:
        from pymbar import MBAR
        N_k = np.array([len(e) for e in energies])
        E = np.concatenate(energies)
        u_kn = np.array([E / kT for kT in kTs])          # reduced potentials, all frames at each T
        mbar = MBAR(u_kn, N_k)
        try:
            w = mbar.weights()                           # pymbar >= 4: (N, K) sample weights per state
        except AttributeError:
            w = mbar.getWeights()                        # pymbar 3
        k_target = int(np.argmin(np.abs(kTs - target_kT)))
        wt = w[:, k_target]
        hist = np.array([wt[(cv_all >= edges[b]) & (cv_all < edges[b + 1])].sum() for b in range(args.bins)])
        hist = np.maximum(hist, 1e-12)
        fes = -target_kT * np.log(hist); fes -= fes.min()
        print('  MBAR FES at kT=%.4f reconstructed' % target_kT)
    except Exception as e:
        print('  (pymbar unavailable -> per-temperature histogram FES: %s)' % e)
        k_target = int(np.argmin(np.abs(kTs - target_kT)))
        h, _ = np.histogram(cvs[k_target], bins=edges, density=True)
        h = np.maximum(h, 1e-12)
        fes = -kTs[k_target] * np.log(h); fes -= fes.min()

    np.save('%s_remd_fes.npy' % args.output_base, np.column_stack([ctr, fes]))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for kT, cv in zip(kTs, cvs):
            ax[0].hist(cv, bins=edges, histtype='step', density=True, label='kT=%.3f' % kT)
        ax[0].set_xlabel('CV'); ax[0].set_ylabel('P(CV)'); ax[0].set_title('per-replica CV'); ax[0].legend(fontsize=6)
        ax[1].plot(ctr, fes, '-o', ms=3)
        ax[1].set_xlabel('CV'); ax[1].set_ylabel('F (E_up)'); ax[1].set_title('FES at kT=%.4f' % target_kT)
        fig.tight_layout(); fig.savefig('%s_remd.png' % args.output_base, dpi=130)
        print('  wrote %s_remd.png' % args.output_base)
    except Exception as e:
        print('  (plot skipped: %s)' % e)


if __name__ == '__main__':
    main()
