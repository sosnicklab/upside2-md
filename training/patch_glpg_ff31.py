"""Patch a trained force field into an existing glpG hybrid seed config.

WHY PATCH RATHER THAN REBUILD. Re-running `popepopg_prep.sbatch` would redo bilayer packing and
multi-stage equilibration, producing a DIFFERENT starting structure. Since the question is
whether TM4 stays helical, a changed starting configuration would confound the force-field
comparison with an initial-condition change on exactly the observable of interest. Patching keeps
the membrane, the packing, the ions and the starting coordinates bit-identical, so the force
field is the only thing that differs. `patch_gly_rama_v2.py` set the precedent in this tree.

WHAT ACTUALLY NEEDS PATCHING. The hybrid config has no `nonlinear_coupling_environment`,
`hbond_coverage` or `hbond_coverage_hydrophobe` node -- it uses the MARTINI environment instead.
So of the five parameters ConDiv trains, `env` does not apply at all, and of the sidechain arrays
only `pair_interaction` is used. `rotamer_center_fixed` and `hydrophobe_placement` move by 2e-16,
machine precision, so the `placement_*` nodes are left alone.

    rama_map_pot/rama_pot     <- recomputed from the trained rama.dat and sheet
    hbond_energy/parameters   <- [:4] scaled by the trained hb
    rotamer/interaction_param <- the trained pair_interaction

    python3 patch_glpg_ff31.py --ff DIR --seed IN.up --out OUT.up [--verify-roundtrip FF21_DIR]

`--verify-roundtrip` patches the ORIGINAL force field back in and checks the result reproduces
the untouched seed. That is the gate: it proves the rama_pot recomputation here matches what the
preparation pipeline actually wrote, rather than being a plausible reimplementation.
"""

import argparse
import os
import shutil
import sys

import numpy as np
import tables as tb

sys.path.insert(0, os.path.join(os.environ['UPSIDE_HOME'], 'py'))
import upside_config as uc          # noqa: E402


def rebuild_rama_pot(seq, rama_library, sheet_file):
    """Reproduce write_rama_map_pot's `rama_pot`, including the final Boltzmann shift."""
    with tb.open_file(rama_library) as tr:
        sheet_restype = [s.decode() if isinstance(s, bytes) else s
                         for s in tr.root.sheet._v_attrs.restype]
    rid = {x: i for i, x in enumerate(sheet_restype)}
    values = np.loadtxt(sheet_file)
    if values.size != len(sheet_restype):
        raise ValueError(f'sheet file has {values.size} values, expected {len(sheet_restype)}')
    sheet = np.array([values[rid['PRO' if s == 'CPR' else s]] for s in seq])

    pot = uc.read_weighted_maps(np.asarray(seq), rama_library, sheet, 'mixture')
    # write_rama_map_pot's last line; a constant per map, invisible in forces but part of the
    # stored values, so omitting it would make the patched file differ from a built one.
    pot = pot - (pot * np.exp(-pot)).sum(axis=(-2, -1), keepdims=True)
    return pot


def rama_library_for(ff_dir):
    """A trained force field ships its own rama.dat; ff_2.1 and earlier use common/rama.dat.

    Same rule as bench_run.py. Getting this wrong silently pairs trained sidechain and hbond
    parameters with ff_2.1's glycine row, which is the one thing ff3.1 changes.
    """
    own = os.path.join(ff_dir, 'rama.dat')
    if os.path.exists(own):
        return own
    return os.path.join(os.environ['UPSIDE_HOME'], 'parameters', 'common', 'rama.dat')


def patch(ff_dir, seed, out):
    shutil.copy(seed, out)
    with tb.open_file(out, 'a') as t:
        seq = [s.decode() if isinstance(s, bytes) else s for s in t.root.input.sequence[:]]
        p = t.root.input.potential

        new_pot = rebuild_rama_pot(seq, rama_library_for(ff_dir),
                                   os.path.join(ff_dir, 'sheet'))
        if new_pot.shape != p.rama_map_pot.rama_pot.shape:
            raise ValueError(f'rama_pot shape {new_pot.shape} != '
                             f'{p.rama_map_pot.rama_pot.shape}')
        d_rama = float(np.abs(np.asarray(p.rama_map_pot.rama_pot[:]) - new_pot).max())
        p.rama_map_pot.rama_pot[:] = new_pot

        with tb.open_file(os.path.join(ff_dir, 'hbond.h5')) as h:
            hbp = h.root.parameter[:]
        old_hb = np.asarray(p.hbond_energy.parameters[:])
        if hbp.shape != old_hb.shape:
            raise ValueError('hbond parameter shape mismatch')
        d_hb = float(np.abs(old_hb - hbp).max())
        p.hbond_energy.parameters[:] = hbp

        with tb.open_file(os.path.join(ff_dir, 'sidechain.h5')) as s:
            pair = s.root.pair_interaction[:]
        # interaction_param sits in a subgroup (rotamer/pair_interaction/), and the hybrid and
        # soluble configs nest it differently, so locate it rather than assuming a path.
        cand = [nd for nd in p.rotamer._f_walknodes('Array')
                if nd._v_name == 'interaction_param']
        if len(cand) != 1:
            raise ValueError(f'expected one interaction_param under rotamer, found {len(cand)}')
        node = cand[0]
        old_pair = np.asarray(node[:])
        if pair.shape != old_pair.shape:
            raise ValueError(f'pair_interaction {pair.shape} != {old_pair.shape}')
        d_pair = float(np.abs(old_pair - pair).max())
        node[:] = pair
    return d_rama, d_hb, d_pair


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ff', required=True)
    ap.add_argument('--seed', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--verify-roundtrip', default=None,
                    help='force-field directory the seed was ORIGINALLY built with')
    a = ap.parse_args()

    if a.verify_roundtrip:
        tmp = a.out + '.roundtrip'
        patch(a.verify_roundtrip, a.seed, tmp)
        with tb.open_file(a.seed) as f0, tb.open_file(tmp) as f1:
            worst = 0.0
            for nd in ('rama_map_pot/rama_pot', 'hbond_energy/parameters',
                       'rotamer/interaction_param'):
                g, n = nd.split('/')
                pick = lambda f: [q for q in getattr(f.root.input.potential, g)
                                  ._f_walknodes('Array') if q._v_name == n][0]
                x = np.asarray(pick(f0)[:]); y = np.asarray(pick(f1)[:])
                d = float(np.abs(x - y).max())
                worst = max(worst, d)
                print(f'  roundtrip {nd:28s} max|diff| {d:.3e}')
        os.remove(tmp)
        if worst > 1e-4:
            sys.exit(f'ROUNDTRIP FAILED: {worst:.3e}. The recomputation does not reproduce what '
                     f'the preparation pipeline wrote; do not patch.')
        print(f'  roundtrip PASSED (worst {worst:.3e})\n')

    d_rama, d_hb, d_pair = patch(a.ff, a.seed, a.out)
    print(f'patched {os.path.basename(a.out)}')
    print(f'  rama_pot          changed by max {d_rama:.4f}')
    print(f'  hbond parameters  changed by max {d_hb:.4f}')
    print(f'  pair_interaction  changed by max {d_pair:.4f}')


if __name__ == '__main__':
    main()
