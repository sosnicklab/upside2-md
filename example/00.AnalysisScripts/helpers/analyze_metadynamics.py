#!/usr/bin/env python
"""Reconstruct the free-energy surface from an Upside metadynamics trajectory.

The metadynamics node logs two quantities per frame:
  output/metadynamics_cv    (n_frame, 2)          [CV value, instantaneous bias]
  output/metadynamics_bias  (n_frame, grid_bins)  accumulated bias on a fixed CV grid

The last frame of metadynamics_bias is the converged bias V(s); the free energy is
  well-tempered (bias_factor gamma > 1):  F(s) = -(gamma/(gamma-1)) V(s)
  standard      (gamma <= 1):             F(s) = -V(s)

Usage:
  python analyze_metadynamics.py traj.up out_base [--target-kT 0.8647]
Outputs out_base_fes.npy (columns: cv, F in E_up), out_base_cv.npy (CV time series),
and out_base_metad.png.
"""
import argparse
import numpy as np
import tables as tb


def _output_groups(t):
    i = 0
    while 'output_previous_%i' % i in t.root:
        yield t.get_node('/output_previous_%i' % i)
        i += 1
    if 'output' in t.root:
        yield t.get_node('/output')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_h5', help='metadynamics .up trajectory')
    p.add_argument('output_base', help='output name prefix')
    p.add_argument('--target-kT', type=float, default=None,
                   help='kT (E_up) for the E_up->kT axis in the plot; default reads none')
    args = p.parse_args()

    with tb.open_file(args.input_h5) as t:
        g = t.root.input.potential.metadynamics._v_attrs
        grid_min, grid_max = float(g.grid_min), float(g.grid_max)
        grid_bins = int(g.grid_bins)
        gamma = float(g.bias_factor)

        cv_parts, bias_last = [], None
        for grp in _output_groups(t):
            if 'metadynamics_cv' not in grp:
                continue
            cv = grp.metadynamics_cv[:]
            cv_parts.append(cv if not cv_parts else cv[1:])   # drop shared restart frame
            bias_last = grp.metadynamics_bias[-1]
        if bias_last is None:
            raise RuntimeError('no metadynamics output found in %s' % args.input_h5)
        cv_ts = np.concatenate(cv_parts)                      # (n_frame, 2)

    s = grid_min + (grid_max - grid_min) * (np.arange(grid_bins) + 0.5) / grid_bins
    fes = -(gamma / (gamma - 1.0)) * bias_last if gamma > 1.0 else -bias_last
    fes -= fes.min()

    np.save('%s_fes.npy' % args.output_base, np.column_stack([s, fes]))
    np.save('%s_cv.npy' % args.output_base, cv_ts)

    span = cv_ts[:, 0].max() - cv_ts[:, 0].min()
    print('metadynamics analysis: %s' % args.input_h5)
    print('  bias_factor gamma = %.3g (%s)' % (gamma, 'well-tempered' if gamma > 1 else 'standard'))
    print('  CV explored %.3f .. %.3f (span %.3f)' % (cv_ts[:, 0].min(), cv_ts[:, 0].max(), span))
    print('  final bias max %.3f E_up ; FES barrier %.3f E_up' % (bias_last.max(), fes.max()))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(cv_ts[:, 0], lw=0.6)
        ax[0].set_xlabel('frame'); ax[0].set_ylabel('CV'); ax[0].set_title('CV vs time')
        ax[1].plot(s, fes, '-')
        ax[1].set_xlabel('CV'); ax[1].set_ylabel('F (E_up)'); ax[1].set_title('free-energy surface')
        if args.target_kT:
            ax2 = ax[1].twinx(); ax2.set_ylim(np.array(ax[1].get_ylim()) / args.target_kT)
            ax2.set_ylabel('F (kT)')
        fig.tight_layout(); fig.savefig('%s_metad.png' % args.output_base, dpi=130)
        print('  wrote %s_metad.png' % args.output_base)
    except Exception as e:
        print('  (plot skipped: %s)' % e)


if __name__ == '__main__':
    main()
