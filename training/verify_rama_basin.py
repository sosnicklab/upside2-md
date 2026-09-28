"""Gate for the rama basin offsets: does each trained map's offset reach exactly its own residues?

The offsets are the only Ramachandran parameters trained, one set per trained (central, direction,
neighbour) map, and each must act on no residue but those that read that map. An indexing slip
(a direction swapped, a cis-proline neighbour not read as PRO, a terminal residue given two maps)
would train one pair's offsets on another pair's residues and never raise anything. So the whole
path the trainer uses is checked, library file -> upside_config -> the per-residue maps in the
.up file:

  1. basins and trained set: continuous across phi = +-180, exactly mirror-symmetric on the grid,
     each map's basin set partitions the torus, and the trained set is exactly GLY|X (alpha_R,
     alpha_L, beta), GLY|GLY (helix and beta, tied) and X|right|PRO (alpha_R, beta);
  2. writer: one map's offset changes that map's coil entry and nothing else in the library; under
     random offsets every other coil entry and the sheet group but GLY|GLY are bitwise the source,
     every written map is normalised like the NDRD maps (an offset is a weight factor on its
     basin), and GLY|GLY's coil and sheet entries stay exactly mirror-symmetric;
  3. reach: for a map of each trained class, and for both termini (the basin under test is switched
     on there), perturbing that map's offsets changes the per-residue maps of exactly the residues
     `residue_keys` says read it and raises each inside the basin it was raised in; and with random
     offsets on every trained map, a residue that reads only untrained maps keeps its NDRD map
     exactly;
  4. engine map: a glycine that reads only GLY|GLY maps, terminal or between two glycines, gets an
     exactly mirror-symmetric map from upside_config's coil/sheet mixture.

    python3 verify_rama_basin.py <training_dir> [protein_code]

Without a code, the first training protein with a glycine pair and a pre-proline residue is used.
"""

import os
import shutil
import sys
import tempfile

import numpy as np
import pickle as cp
import tables as tb

sys.path.insert(0, os.path.join(os.environ.get('UPSIDE_HOME', '..'), 'py'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rama_basin as rb                 # noqa: E402
import run_upside as ru                 # noqa: E402
from upside_config import read_fasta, read_weighted_maps    # noqa: E402


def pick_protein(D):
    codes = [l.split()[0] for l in open(os.path.join(D, 'pdb_list'))][1:]
    for c in codes:
        f = os.path.join(D, 'upside_input', c + '.fasta')
        seq = list(read_fasta(open(f))) if os.path.exists(f) else []
        gg = any(a == b == 'GLY' for a, b in zip(seq[1:-1], seq[2:-1]))
        pre_pro = any(a != 'GLY' and b == 'PRO' for a, b in zip(seq[1:-1], seq[2:]))
        if gg and pre_pro:
            return c
    sys.exit('no training protein has both a glycine pair and a pre-proline residue')


def rama_pot(D, code, library, work, init_npy):
    P = os.path.join(D, 'init_param')
    out = os.path.join(work, 'c.up')
    ru.upside_config(
        os.path.join(D, 'upside_input', code + '.fasta'), out,
        environment_potential=os.path.join(P, 'environment.h5'), environment_potential_type=1,
        bb_environment_potential=os.path.join(P, 'bb_env.dat'),
        rotamer_interaction=os.path.join(P, 'sidechain.h5'),
        rotamer_placement=os.path.join(P, 'sidechain.h5'), initial_structure=init_npy,
        hbond_energy=os.path.join(P, 'hbond.h5'), rama_sheet_mix_energy=os.path.join(P, 'sheet'),
        dynamic_rotamer_1body=True, rama_library=library, rama_param_deriv=True,
        reference_state_rama=os.path.join(D, 'upside_input', 'rama_reference.pkl'))
    with tb.open_file(out) as t:
        g = t.root.input.potential.rama_map_pot
        return g.rama_pot[:][g.rama_map_id[:]]


def main():
    D = os.path.abspath(sys.argv[1])
    code = sys.argv[2] if len(sys.argv) > 2 else pick_protein(D)
    work = tempfile.mkdtemp()
    rs = np.random.RandomState(0)
    src = os.path.join(D, 'upside_input', 'rama.dat')
    ok = True

    # --- 1: basins and the trained set --------------------------------------------------------
    w = rb.grid_basins(72)
    mirror_err = max(np.abs(rb.mirror(w[..., b]) - w[..., m]).max() for b, m in enumerate(rb.MIRROR[:6]))
    seam = max(np.abs(rb.basin_weights(np.pi - 1e-9, p) - rb.basin_weights(-np.pi, p)).max()
               for p in np.linspace(-np.pi, np.pi, 13))
    cover = [w[..., list(p)].sum(-1) for p in (rb.PARTITION, rb.PARTITION_GLY)]
    st = rb.init_state(src)
    keys, restype, dirs = st['keys'], st['restype'], st['dirs']
    gly, pro, right = restype.index('GLY'), restype.index('PRO'), dirs.index('right')
    spec = np.zeros_like(st['active'])
    for k, (c, d, n) in enumerate(keys):
        if c == gly:
            spec[k, [0, 1, 2, 4] if n == gly else [0, 1, 2]] = True
        elif d == right and n == pro:
            spec[k, [0, 2]] = True
    trained = st['active'].any(axis=1)
    as_spec = np.array_equal(st['active'], spec)
    print(f'1. basins      mirror error {mirror_err:.1e}, jump across phi = +-180 {seam:.1e}; basin '
          f'sets cover the torus {cover[0].min():.3f}-{cover[0].max():.3f} (five), '
          f'{cover[1].min():.3f}-{cover[1].max():.3f} (central glycine); trained {trained.sum()} of '
          f'{len(keys)} maps, {rb.n_param(st)} offsets, exactly GLY|X, GLY|GLY, X|right|PRO: {as_spec}')
    ok &= mirror_err < 1e-12 and seam < 1e-6 and min(c.min() for c in cover) > 0.99 and as_spec

    # --- 2: writer -----------------------------------------------------------------------------
    with tb.open_file(src) as t:
        src_coil, src_sheet = t.root.coil.dimer_pot[:], t.root.sheet.dimer_pot[:]
        sheet_restype = rb._decode(t.root.sheet._v_attrs.restype)
    gs = sheet_restype.index('GLY')
    sheet_other = np.ones(src_sheet.shape[:3], dtype=bool)
    sheet_other[gs, :, gs] = False
    coil_other = np.ones(src_coil.shape[:3], dtype=bool)
    for c, d, n in keys[trained]:
        coil_other[c, d, n] = False

    def coil_of(state, name):
        with tb.open_file(rb.write_library(state, src, os.path.join(work, name))) as t:
            return t.root.coil.dimer_pot[:], t.root.sheet.dimer_pot[:]
    coil0, _ = coil_of(st, 'zero.dat')
    k = int(rs.choice(np.where(trained & ~st['is_gg'])[0]))
    st['offset'][k, 2] = 0.8
    coil, sheet = coil_of(st, 'one.dat')
    moved = {tuple(x) for x in np.argwhere(np.any(np.nan_to_num(coil - coil0) != 0, axis=(-1, -2)))}
    expect = {tuple(keys[k])}
    zero_dev = max(np.abs(coil0[c, d, n].astype(float) - src_coil[c, d, n]).max()
                   for c, d, n in keys[trained & ~st['is_gg']])
    norm = max(abs(np.log(np.exp(-m.astype(float)).sum()))
               for m in [coil[c, d, n] for c, d, n in keys[trained]] + list(sheet[gs, :, gs]))
    st['offset'][:] = 0.
    acc = rb.empty_accumulator(st)
    acc['n'][:] = 50.
    acc['free'][:] = rs.rand(len(keys), rb.N_BASIN) * 50.
    acc['native'][:] = rs.rand(len(keys), rb.N_BASIN) * 50.
    st['T0'], st['steps_per_round'] = 0.8, 19
    rb.update(st, acc)
    lib = rb.write_library(st, src, os.path.join(work, 'rand.dat'))
    with tb.open_file(lib) as t:
        rand_coil, rand_sheet = t.root.coil.dimer_pot[:], t.root.sheet.dimer_pot[:]
    untouched = (np.array_equal(rand_coil[coil_other], src_coil[coil_other], equal_nan=True) and
                 np.array_equal(rand_sheet[sheet_other], src_sheet[sheet_other], equal_nan=True))
    asym = rb.gly_gly_asymmetry(lib)
    print(f'2. writer      one offset moved {sorted(tuple(int(i) for i in x) for x in moved)} '
          f'(expected {tuple(int(i) for i in keys[k])}); under random offsets everything but the '
          f'trained entries and GLY|GLY sheet is bitwise the source: {untouched}; zero offsets '
          f'reproduce the source to {zero_dev:.1e}; every written map normalised to {norm:.1e}; '
          f'GLY|GLY coil and sheet asymmetry {asym:.1e}')
    ok &= moved == expect and untouched and zero_dev < 1e-4 and norm < 1e-4 and asym == 0.

    # --- 3: reach, through upside_config -------------------------------------------------------
    st = rb.init_state(src)
    seq = list(read_fasta(open(os.path.join(D, 'upside_input', code + '.fasta'))))
    rk = rb.residue_keys(st, seq)
    init = cp.load(open(os.path.join(D, 'upside_input', code + '.initial.pkl'), 'rb'),
                   encoding='latin1')
    init_npy = os.path.join(work, 'init.npy')
    np.save(init_npy, init[:, :, 0])
    base_pot = rama_pot(D, code, rb.write_library(st, src, os.path.join(work, 'b.dat')), work, init_npy)

    use = {}
    for i, ks in enumerate(rk):
        for kk in ks:
            use.setdefault(kk, set()).add(i)
    first = lambda cond: next(kk for kk in sorted(use, key=lambda kk: -len(use[kk])) if cond(kk))
    cases = [('GLY|X', first(lambda kk: keys[kk][0] == gly and not st['is_gg'][kk]), 1),
             ('GLY|GLY', first(lambda kk: st['is_gg'][kk]), 0),
             ('pre-proline', first(lambda kk: trained[kk] and keys[kk][0] != gly), 0),
             ('first residue', rk[0][0], 3),
             ('last residue', rk[-1][0], 2)]
    print(f'3. reach       {code}, {len(seq)} residues')
    for name, kk, b in cases:
        sc = rb.init_state(src)
        sc['active'][kk, b] = True
        sc['offset'][kk, b] = 0.5
        if sc['is_gg'][kk]:
            sc['active'][kk, rb.MIRROR[b]] = True
            sc['offset'][kk, rb.MIRROR[b]] = 0.5
        pot = rama_pot(D, code, rb.write_library(sc, src, os.path.join(work, 'p.dat')), work, init_npy)
        d = pot - base_pot
        changed = set(np.where(np.abs(d).max(axis=(1, 2)) > 1e-6)[0])
        inside = w[..., b] > 0.9
        contrast = min((d[i][inside].mean() - d[i][w[..., b] < 0.1].mean()) for i in use[kk])
        ok &= changed == use[kk] and contrast > 0.
        c, dd, n = keys[kk]
        print(f'   {name:>16}  {restype[c]}|{dirs[dd]}|{restype[n]} {rb.BASINS[b]:>8}'
              f'{"" if trained[kk] else " (switched on)"}: {len(changed)} residues changed, '
              f'{len(use[kk])} read it, same set {changed == use[kk]}; weakest in-basin rise '
              f'{contrast:+.3f}')
    sr = rb.init_state(src)
    sr['offset'] = rs.randn(*sr['offset'].shape) * 0.5 * sr['active']
    g = sr['is_gg']
    sr['offset'][g] = 0.5 * (sr['offset'][g] + sr['offset'][g][:, rb.MIRROR])
    d = np.abs(rama_pot(D, code, rb.write_library(sr, src, os.path.join(work, 'r.dat')), work,
                        init_npy) - base_pot).max(axis=(1, 2))
    reads = np.array([any(trained[kk] for kk in ks) for ks in rk])
    print(f'   random offsets on every trained map: {reads.sum()} residues read one, all changed: '
          f'{bool((d[reads] > 1e-6).all())}; {(~reads).sum()} read none, largest change '
          f'{d[~reads].max():.1e}')
    ok &= bool((d[reads] > 1e-6).all()) and d[~reads].max() == 0.

    # --- 4: the map the engine gets, through upside_config's coil/sheet mixture ----------------
    seq = ['GLY', 'GLY', 'GLY', 'ALA']
    sheet_E = np.loadtxt(os.path.join(D, 'init_param', 'sheet'))
    pots = read_weighted_maps(seq, lib, sheet_E[[sheet_restype.index(s) for s in seq]])
    res_asym = [float(np.abs(p - rb.mirror(p)).max()) for p in pots[:2]]
    print(f'4. engine map  G-G-G-A under random offsets, |E - mirror(E)| of the terminal glycine '
          f'{res_asym[0]:.1e}, of the glycine between glycines {res_asym[1]:.1e}')
    ok &= max(res_asym) == 0.

    shutil.rmtree(work)
    print('\n' + ('PASS: every offset reaches exactly its own residues.' if ok else
                  'FAIL: see the parts above.'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
