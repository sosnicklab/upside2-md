"""Build a Ramachandran library whose central-glycine row is the AWH-measured free energy.

    python3 build_gly_library.py <source_library> <awh_dir> <rama_reference.pkl> <sheet> <out_library>

WHY. A library map is -ln P over residues of folded PDB structures, so for glycine it is part local
energy and part selection: evolution places glycine where a fold needs a left-handed residue.
Upside applies the map as pure energy, which pushes helical glycines toward alpha_L. Data from
folded proteins cannot separate the two parts (findings 1.15-1.16), so glycine's row is replaced by
a measurement in which no fold selects anything: 2D AWH on (phi, psi) of the central glycine in the
capped dipeptides Ac-X-Gly-NHMe (contexts L<X>) and Ac-Gly-X-NHMe (R<X>), 300 K, explicit water.
Every other row of the library is copied unchanged.

THE ROW, one map per neighbour class, each from its own measurement:
  * GLY|GLY: the two Gly-Gly surfaces (LG, RG), pooled and made mirror-symmetric. An achiral
    residue between achiral neighbours has no handedness.
  * GLY|right|PRO: the Ac-Gly-Pro-NHMe surface (RP) alone. The proline ring empties both helical
    basins (alpha_R 0.006, alpha_L 0.004), twenty-fold below every other context.
  * every other GLY|X, both directions: the remaining 37 surfaces pooled. Per-neighbour structure
    does not reproduce between replicas (findings, 2026-09-18 verdict); the pool carries the
    measured left-handed bias of an L-neighbour.
Pooling and symmetrising average probabilities, never energies: averaging energies is a geometric
mean of probabilities, which shrinks any basin that varies between the surfaces.

THE SAME ENTRY IN COIL AND SHEET. upside_config mixes every residue's coil map with its sheet map by
the sheet mixing energy. The sheet maps are statistics of strand residues, selection again, so a
central glycine reads the measured entry in both groups and the mixture leaves it unchanged.

THE REFERENCE CORRECTION. Every config path adds `rama_map_pot_ref` (log rama_reference.pkl, mean
removed) to every residue as an energy (ConDiv, the benchmark, the MARTINI hybrid, the examples).
A measured surface needs no reference state, so each glycine entry is stored with that correction
subtracted, and the engine applies the measured surface exactly. A library entry is therefore not
itself mirror-symmetric for GLY|GLY; the engine's total is.

UNITS. A map holds -ln P normalised so that sum(exp(-E)) == 1 on the 72 x 72 grid (node i at
-180 + 5i degrees), so a PMF enters as PMF / kT(300 K), not divided by the kJ/mol-per-E_up factor.
Each AWH surface (46 x 46) is interpolated with periodic cubic splines in energy before pooling.

Each entry is one directional map. upside_config still combines a residue's left and right maps by
its mixture rule, unchanged, so a glycine before a proline reads the RP entry mixed with its left
neighbour's entry.

CHECKS, through upside_config's own `read_weighted_maps`, before anything is written as final:
every non-glycine residue's map identical to the source library's; every glycine whose maps are
all one entry equal to that measured surface up to a constant (X-G-Y and a terminal glycine read
the pooled surface, G-G-G and a terminal G-G the symmetric one); those that read GLY|GLY exactly
mirror-symmetric. Any failure exits nonzero and removes the output.
"""

import glob
import os
import pickle
import shutil
import sys

import numpy as np
import tables as tb
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, os.path.join(os.environ['UPSIDE_HOME'], 'py'))
import upside_config as uc  # noqa: E402

KT = 8.31446261815324e-3 * 300.0                       # kJ/mol at the AWH temperature
ONE = 'ACDEFGHIKLMNPQRSTVWY'                          # the AWH contexts are L<X> and R<X>
BLANKS = ('LG', 'RG')
PRE_PRO = 'RP'
TOL = 1e-3                                             # E_up; the library is float32


def read_awh(path):
    """A `gmx awh` 2D surface as a (phi, psi) array of PMF in kJ/mol on its periodic grid."""
    d = np.loadtxt(path, comments=('#', '@'))
    phi, psi = np.unique(d[:, 0]), np.unique(d[:, 1])
    n = len(phi)
    if len(psi) != n or not np.allclose(phi, -180. + 360. * np.arange(n) / n, atol=1e-3):
        raise ValueError(f'{path}: not a periodic square grid starting at -180')
    m = np.full((n, n), np.nan)
    i, j = np.searchsorted(phi, d[:, 0]), np.searchsorted(psi, d[:, 1])
    m[i, j] = d[:, 2]
    if not np.isfinite(m).all():
        raise ValueError(f'{path}: incomplete surface')
    return m


def latest(system_dir):
    files = glob.glob(os.path.join(system_dir, 'fe_t*.xvg'))
    if not files:
        raise FileNotFoundError(f'no AWH surface in {system_dir}')
    return max(files, key=lambda p: int(p.split('_t')[-1].split('.')[0]))


def to_grid(m, n_out):
    """Periodic cubic interpolation of a square angular grid onto n_out x n_out."""
    n = m.shape[0]
    x = -180. + 360. * np.arange(n) / n
    xx = np.concatenate([x - 360., x, x + 360.])
    t = -180. + 360. * np.arange(n_out) / n_out
    return RectBivariateSpline(xx, xx, np.tile(m, (3, 3)), kx=3, ky=3, s=0)(t, t)


def normalise(e):
    lo = e.min()
    return e - lo + np.log(np.exp(-(e - lo)).sum())


def mirror(m):
    """(phi,psi) -> (-phi,-psi) on the grid: a reversal with a roll, because node i maps to (-i) % n."""
    return np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)


def pool(energies):
    """-ln of the mean probability of several normalised surfaces."""
    return normalise(-np.log(np.mean([np.exp(-e) for e in energies], axis=0)))


def symmetrise(e):
    """-ln of the mean of a surface's probabilities and its mirror's; exactly symmetric."""
    return normalise(np.log(2.) - np.logaddexp(-e, -mirror(e)))


def handedness(e):
    """ln(P(alpha_R) / P(alpha_L)) of a map at T = 1, in the boxes used throughout findings."""
    t = -180. + 360. * np.arange(e.shape[0]) / e.shape[0]
    phi, psi = np.meshgrid(t, t, indexing='ij')
    p = np.exp(-(e - e.min()))
    ar = (phi < 0) & (phi > -160) & (psi > -100) & (psi < 50)
    al = (phi > 0) & (phi < 160) & (psi < 100) & (psi > -50)
    return np.log(p[ar].sum() / p[al].sum())


def decode(a):
    return [x.decode() if isinstance(x, bytes) else x for x in a]


def measured_surfaces(awh_dir, n_grid):
    surfaces = {}
    for one in ONE:
        for side in 'LR':
            path = latest(os.path.join(awh_dir, side + one))
            surfaces[side + one] = normalise(to_grid(read_awh(path), n_grid) / KT)
            ns = int(path.split('_t')[-1].split('.')[0]) / 1000.
            print(f'  {side + one}  {ns:6.1f} ns  ln(aR/aL) {handedness(surfaces[side + one]):+.3f}')
    return surfaces


def write_row(group, entries):
    """Replace every finite central-glycine entry of one library group."""
    restype, dirs = decode(group._v_attrs.restype), decode(group._v_attrs.dir)
    g = restype.index('GLY')
    pot = group.dimer_pot[:]
    for d, direction in enumerate(dirs):
        for n, neighbour in enumerate(restype):
            if np.isnan(pot[g, d, n]).all():
                continue                                # a neighbour column no residue reads
            if neighbour == 'GLY':
                pot[g, d, n] = entries['gg']
            elif neighbour == 'PRO' and direction == 'right':
                pot[g, d, n] = entries['pre_pro']
            else:
                pot[g, d, n] = entries['pool']
    group.dimer_pot[:] = pot


def check(source, out, ref, sheet_file, measured):
    """The library through upside_config, against the source library and the measurements."""
    with tb.open_file(source) as t:
        sheet_restype = decode(t.root.sheet._v_attrs.restype)
    sheet_energy = np.loadtxt(sheet_file)
    tests = ['ALA GLY ALA THR GLY VAL ARG GLY GLU ASN GLY ILE',
             'GLY GLY GLY ALA GLY GLY LEU',
             'ALA GLY PRO SER GLY PRO LYS',
             'VAL ALA LEU ASP PRO GLU PHE TRP']
    failures = []
    for text in tests:
        seq = np.array(text.split())
        sheet = sheet_energy[[sheet_restype.index(s) for s in seq]]
        old = uc.read_weighted_maps(seq, source, sheet)
        new = uc.read_weighted_maps(seq, out, sheet)
        print(f'\n  {text}')
        for i, aa in enumerate(seq):
            if aa != 'GLY':
                if not np.array_equal(old[i], new[i]):
                    failures.append(f'{text}: non-glycine residue {i} {aa} changed')
                continue
            total = new[i] + ref
            left = seq[i - 1] if i > 0 else None
            right = seq[i + 1] if i + 1 < len(seq) else None
            ctx = f'{left or "-"}-GLY-{right or "-"}'
            read = set()                                # the entry of each map this glycine reads
            if left is not None:
                read.add('gg' if left == 'GLY' else 'pool')
            if right is not None:
                read.add('gg' if right == 'GLY' else 'pre_pro' if right == 'PRO' else 'pool')
            note = 'a mixture of two measured entries'
            if len(read) == 1:
                dev = total - measured[read.pop()]
                dev = np.abs(dev - dev.mean()).max()
                if dev > TOL:
                    failures.append(f'{text}: {ctx} at {i} differs from its measured surface ({dev:.2e})')
                note = f'max|engine - measured| {dev:.1e}'
            if 'GLY' in (left, right) and set((left, right)) <= {'GLY', None}:
                asym = np.abs(total - mirror(total)).max()
                if asym > TOL:
                    failures.append(f'{text}: {ctx} at {i} not mirror-symmetric ({asym:.2e})')
                note += f', max|E - mirror E| {asym:.1e}'
            print(f'    {i:2d} {ctx:16s} ln(aR/aL) {handedness(total):+.3f} '
                  f'(NDRD {handedness(old[i] + ref):+.3f})  {note}')
    return failures


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    source, awh_dir, reference, sheet_file, out = sys.argv[1:]

    with tb.open_file(source) as t:
        n_grid = t.root.coil.dimer_pot.shape[-1]
    ref = np.log(pickle.load(open(reference, 'rb'), encoding='latin1'))
    ref -= ref.mean()
    if ref.shape != (n_grid, n_grid):
        sys.exit(f'reference correction {ref.shape} does not match the {n_grid}-node library grid')

    print(f'AWH surfaces in {awh_dir}:')
    s = measured_surfaces(awh_dir, n_grid)
    pooled = [k for k in s if k not in BLANKS + (PRE_PRO,)]
    measured = dict(gg=symmetrise(pool([s[k] for k in BLANKS])),
                    pre_pro=s[PRE_PRO],
                    pool=pool([s[k] for k in pooled]))
    print(f'\nmeasured maps:  GLY|X pooled over {len(pooled)} contexts ln(aR/aL) '
          f'{handedness(measured["pool"]):+.3f};  GLY|GLY {handedness(measured["gg"]):+.3f};  '
          f'GLY|right|PRO {handedness(measured["pre_pro"]):+.3f}')

    entries = {k: normalise(v - ref) for k, v in measured.items()}
    shutil.copy(source, out)
    with tb.open_file(out, 'a') as t:
        write_row(t.root.coil, entries)
        write_row(t.root.sheet, entries)
        t.root._v_attrs.glycine_row = np.bytes_(
            'central-glycine coil and sheet entries: AWH capped dipeptides (ff99SB-ILDN, 300 K), '
            'GLY|GLY from LG+RG symmetrised, GLY|right|PRO from RP, other GLY|X pooled; stored '
            'with rama_reference subtracted; built by training/build_gly_library.py from '
            + os.path.basename(source))

    failures = check(source, out, ref, sheet_file, measured)
    if failures:
        os.remove(out)
        sys.exit('\nFAILED, output removed:\n  ' + '\n  '.join(failures))
    print(f'\nPASS: written {out}')


if __name__ == '__main__':
    main()
