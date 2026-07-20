"""Membrane-aware HDX protection state for the Upside/dry-MARTINI hybrid.

An amide NH is exchange-competent (unprotected) only when it is BOTH
  (i)  open           -- its backbone hydrogen bond is broken, and
  (ii) water-accessible -- it lies outside the bilayer hydrophobic core.
Otherwise it is protected. In a dry (waterless) membrane the accessibility
term is essential: a hydrogen-bond-broken amide buried between the lipid
phosphate planes cannot exchange and must be scored protected, or every
transmembrane segment is wrongly reported as exchanging.

Accessibility is approximated from bilayer depth. Per frame the phosphate
(PO4) beads define the two leaflet planes; an amide N between them is in the
hydrophobic core (protected), beyond them it is in the aqueous region
(accessible). The full hybrid trajectory is read from /output/pos and the
Upside engine is run over all beads to obtain the backbone hydrogen-bond score.

The saved array follows the HX-MS pipeline convention: shape (n_frames,
n_donor), 1 = protected, so it is a drop-in for steps 4-6 of
example/00.AnalysisScripts. With no lipids present it reduces to the
hydrogen-bond-only criterion.
"""

import argparse
import os
import sys

import numpy as np
import tables as tb

upside_path = os.environ['UPSIDE_HOME']
sys.path.insert(0, os.path.expanduser(upside_path + '/py'))
import upside_engine as ue


def leaflet_planes(z_po4):
    """Bilayer midplane and the two phosphate-plane z-positions.

    PO4 z is bimodal (one cluster per leaflet), so the mean is the midplane;
    the leaflet planes are the mean of the beads above and below it.
    """
    z0 = z_po4.mean()
    upper = z_po4[z_po4 > z0]
    lower = z_po4[z_po4 < z0]
    return z0, (lower.mean() if lower.size else z0), (upper.mean() if upper.size else z0)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('top_h5', help='hybrid run file holding /input and /output (.up)')
    p.add_argument('output_npy', help='output protection-state .npy (1 = protected)')
    p.add_argument('--traj-h5', default=None,
                   help='trajectory file with /output/pos (default: top_h5)')
    p.add_argument('--start', type=int, default=0, help='(default 0) first frame')
    p.add_argument('--stride', type=int, default=1, help='(default 1) frame stride')
    p.add_argument('--hbond-cut', type=float, default=0.01,
                   help='(default 0.01) backbone H-bond score above which the NH is bonded')
    p.add_argument('--core-margin', type=float, default=0.0,
                   help='(default 0.0) Angstrom inset from each phosphate plane; '
                        'positive shrinks the protected core so more amides can exchange')
    p.add_argument('--residue', default=None, help='write donor residue ids to this file')
    p.add_argument('--report-raw-data', action='store_true',
                   help='also save H-bond score, signed depth, and buried flag')
    args = p.parse_args()

    traj_h5 = args.traj_h5 or args.top_h5
    for path, desc in [(args.top_h5, 'topology'), (traj_h5, 'trajectory')]:
        if not os.path.isfile(path):
            raise FileNotFoundError('cannot read %s file: %s' % (desc, path))

    with tb.open_file(args.top_h5, 'r') as t:
        donors = t.root.input.potential.infer_H_O.donors
        donor_res = donors.residue[:]
        donor_N = donors.id[:, 1].astype(int)  # middle id column is the amide N atom
        particle_class = t.root.input.particle_class[:]
        atom_names = np.array([x.decode() for x in t.root.input.atom_names[:]])
    n_donor = donor_res.size
    po4 = np.where((particle_class == b'LIPID') & (atom_names == 'PO4'))[0]

    with tb.open_file(traj_h5, 'r') as t:
        pos = np.asarray(t.root.output.pos[:]).squeeze()  # (n_frame, n_atom, 3)
    if pos.ndim == 2:
        pos = pos[None]
    frames = range(args.start, pos.shape[0], args.stride)
    print('{} frames are used ({} donors, {} PO4 beads)'.format(
        len(frames), n_donor, po4.size))

    engine = ue.Upside(args.top_h5)

    PS = np.zeros((len(frames), n_donor), dtype=np.float32)
    Hbond = np.zeros_like(PS)
    Depth = np.zeros_like(PS)
    Buried = np.zeros_like(PS)
    for row, i in enumerate(frames):
        xyz = pos[i].astype(np.float32)
        engine.energy(xyz)
        hb = engine.get_output('protein_hbond')[:n_donor, 6]
        bonded = hb > args.hbond_cut

        z = xyz[:, 2]
        z_N = z[donor_N]
        if po4.size:
            z0, lo, up = leaflet_planes(z[po4])
            buried = (z_N > lo + args.core_margin) & (z_N < up - args.core_margin)
        else:
            z0 = 0.0
            buried = np.zeros(n_donor, dtype=bool)

        PS[row] = bonded | buried  # protected = H-bonded OR buried in the lipid core
        Hbond[row] = hb
        Depth[row] = z_N - z0
        Buried[row] = buried

    base = args.output_npy[:-4] if args.output_npy.endswith('.npy') else args.output_npy
    np.save(base, PS)

    exch = 1.0 - PS.mean(axis=0)                 # membrane criterion open fraction
    exch_hb_only = (Hbond <= args.hbond_cut).mean(axis=0)  # naive H-bond-only open fraction
    print('per-donor exchange (open) probability:')
    print('  H-bond-only mean {:.3f}  |  open-AND-accessible mean {:.3f}'.format(
        exch_hb_only.mean(), exch.mean()))
    print('  amides rescued by membrane burial (open but not accessible): '
          '{:.3f} mean fraction of frames'.format((exch_hb_only - exch).mean()))

    if args.report_raw_data:
        np.save('{}_Hbond.npy'.format(base), Hbond)
        np.save('{}_Depth.npy'.format(base), Depth)
        np.save('{}_Buried.npy'.format(base), Buried)
    if args.residue:
        np.savetxt(args.residue, donor_res, fmt='%i')


if __name__ == '__main__':
    main()
