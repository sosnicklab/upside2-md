"""Gate for the trainable glycine row: is the analytic gradient the real one?

The gradient in `rama_gly_gradient` is analytic rather than finite-differenced, which is the only
reason 42 trainable 72x72 maps are affordable at all. Analytic gradients fail silently: a wrong
softmax factor, a flipped sign, or a branch reading the wrong neighbour's map trains steadily in
the wrong direction and never raises anything. So it is checked end to end against the pipeline
the trainer uses, library file -> upside_config -> engine, in three parts:

  1. the Python spline reproduces the engine's rama energy;
  2. every glycine's reconstructed per-residue map equals the rama_pot upside_config wrote, with
     all 42 row maps perturbed so that they differ, which pins the direction/neighbour indexing;
  3. directional finite differences match the analytic gradient, for all maps at once, for the
     GLY|GLY maps alone, and for the single X|GLY map most used by the protein.

    python3 verify_gly_gradient.py <training_dir> [protein_code]

Without a code, the first training protein that has a terminal glycine and a Gly-Gly pair is used,
so that every branch of the mixture is exercised.
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
import rama_gly_gradient as rg          # noqa: E402
import run_upside as ru                 # noqa: E402
import upside_engine as ue              # noqa: E402
from upside_config import read_fasta       # noqa: E402


def read_seq(fasta):
    """The sequence exactly as the trainer and upside_config read it."""
    return list(read_fasta(open(fasta)))


def branches(seq):
    n_int = sum(1 for i in range(1, len(seq) - 1) if seq[i] == 'GLY')
    n_term = int(seq[0] == 'GLY') + int(seq[-1] == 'GLY')
    n_gg = sum(1 for i in range(1, len(seq) - 1)
               if seq[i] == 'GLY' and (seq[i - 1] == 'GLY' or seq[i + 1] == 'GLY'))
    return n_int, n_term, n_gg


def pick_protein(D):
    codes = [l.split()[0] for l in open(os.path.join(D, 'pdb_list'))][1:]
    for c in codes:
        f = os.path.join(D, 'upside_input', c + '.fasta')
        if os.path.exists(f) and all(branches(read_seq(f))):
            return c
    sys.exit('no training protein has interior, terminal and Gly-Gly glycines')


def build(D, code, library, work, tag, init_npy):
    P = os.path.join(D, 'init_param')
    out = os.path.join(work, tag + '.up')
    ru.upside_config(
        os.path.join(D, 'upside_input', code + '.fasta'), out,
        environment_potential=os.path.join(P, 'environment.h5'),
        environment_potential_type=1,
        bb_environment_potential=os.path.join(P, 'bb_env.dat'),
        rotamer_interaction=os.path.join(P, 'sidechain.h5'),
        rotamer_placement=os.path.join(P, 'sidechain.h5'),
        initial_structure=init_npy,
        hbond_energy=os.path.join(P, 'hbond.h5'),
        rama_sheet_mix_energy=os.path.join(P, 'sheet'),
        dynamic_rotamer_1body=True, rama_library=library, rama_param_deriv=True,
        reference_state_rama=os.path.join(D, 'upside_input', 'rama_reference.pkl'))
    return out


def rama_energy(config, pos):
    e = ue.Upside(config)
    e.energy(pos)
    return float(e.get_output('rama_map_pot')[0, 0]), e.get_output('rama_coord')


def main():
    D = os.path.abspath(sys.argv[1])
    code = sys.argv[2] if len(sys.argv) > 2 else pick_protein(D)
    work = tempfile.mkdtemp()
    rs = np.random.RandomState(0)

    seq = read_seq(os.path.join(D, 'upside_input', code + '.fasta'))
    n_int, n_term, n_gg = branches(seq)
    print(f'{code}: {len(seq)} residues, {n_int} interior glycines, '
          f'{n_term} terminal, {n_gg} with a glycine neighbour')
    for what, n in (('terminal', n_term), ('glycine-neighbour', n_gg)):
        if n == 0:
            print(f'  NOTE: no {what} glycine here, so that branch is NOT exercised by this run')
    if n_int == 0:
        sys.exit(f'{code} has no interior glycine; pick another protein')

    init = cp.load(open(os.path.join(D, 'upside_input', code + '.initial.pkl'), 'rb'),
                   encoding='latin1')
    init_npy = os.path.join(work, 'init.npy')
    np.save(init_npy, init[:, :, 0])
    pos = np.load(init_npy)

    src = os.path.join(D, 'upside_input', 'rama.dat')
    sheet_file = os.path.join(D, 'init_param', 'sheet')

    # A non-trivial point: the starting row plus a smooth random perturbation of every map, so no
    # two branches coincide and a branch reading the wrong neighbour cannot pass by accident.
    keys, is_gg, G0 = rg.start_row(src)
    dG = np.stack([rg.fourier_lowpass(rs.randn(*G0.shape[1:]), 8) for _ in keys])
    G = rg.constrain(G0 + 0.3 * dG / np.abs(dG).max(), is_gg)

    lib0 = rg.write_row(src, os.path.join(work, 'lib0.dat'), keys, G)
    cfg0 = build(D, code, lib0, work, 'base', init_npy)
    E0, coord = rama_energy(cfg0, pos)

    ok = True
    # --- 1: does the Python spline reproduce the engine? ---------------------------------
    with tb.open_file(cfg0) as t:
        g = t.root.input.potential.rama_map_pot
        P, rid, mid = g.rama_pot[:], g.residue_id[:], g.rama_map_id[:]
    card = rg.cardinal_function(P.shape[-1])
    recon = sum(rg.energy_from_map(coord[rid[k]][None, :], P[mid[k]], card) for k in range(len(rid)))
    print(f'\n1. spline      engine {E0:.6f}   python {recon:.6f}   diff {abs(recon - E0):.2e}')
    ok &= abs(recon - E0) < 1e-3

    # --- 2: is every glycine's reconstructed map the one upside_config wrote? ------------
    import torch
    chain = rg.GlycineMapChain(seq, lib0, sheet_file, keys)
    Gt = torch.as_tensor(G, dtype=torch.float64)
    worst = 0.
    res_map = {int(r): P[m] for r, m in zip(rid, mid)}
    for r in chain.residues:
        diff = np.abs(chain.residue_map(Gt, r).numpy() - res_map[r['index']]).max()
        worst = max(worst, diff)
    print(f'2. maps        {len(chain.residues)} glycines, worst |python - rama_pot| {worst:.2e}')
    ok &= worst < 1e-3

    # --- 3: does the chain rule through the mixture match finite differences? -------------
    res_grad = {r['index']: rg.map_gradient(coord[r['index']][None, :], P.shape[-1], card)
                for r in chain.residues}
    grad = chain.backprop(G, res_grad)
    used = np.bincount([k for r in chain.residues for k in r['map_k']], minlength=len(keys))

    # The library stores dimer_pot as float32, so a perturbation of eps has a rounding error of
    # about 5e-7 per cell against a map value of order 10. Too small an eps and the finite
    # difference is rounding noise; too large and it picks up the real curvature of the
    # log-sum-exp mixtures. The analytic value has to sit inside the bowl, not match every eps.
    eps_list = (3e-2, 1e-2, 3e-3, 1e-3)
    print(f'\n3. {"direction":>14} {"analytic":>12} | '
          + ' '.join(f'{"eps=" + f"{e:g}":>22}' for e in eps_list))
    masks = [('all maps', np.ones(len(keys), bool)),
             ('GLY|GLY only', is_gg),
             ('one X|GLY', np.arange(len(keys)) == np.argmax(np.where(is_gg, -1, used)))]
    for name, sel in masks:
        if not (used[sel] > 0).any():
            print(f'   {name:>14}   not used by {code}, skipped')
            continue
        d = np.zeros_like(G)
        d[sel] = np.stack([rg.fourier_lowpass(rs.randn(*G.shape[1:]), 8) for _ in range(sel.sum())])
        d = rg.constrain(d, is_gg)
        d /= np.abs(d).max()
        analytic = float((grad * d).sum())
        cells, best = [], np.inf
        for eps in eps_list:
            vals = []
            for s_ in (+1, -1):
                lib = rg.write_row(src, os.path.join(work, f'l{s_}.dat'), keys, G + s_ * eps * d)
                vals.append(rama_energy(build(D, code, lib, work, f'c{s_}', init_npy), pos)[0])
            fd = (vals[0] - vals[1]) / (2 * eps)
            rel = abs(analytic - fd) / max(1e-12, abs(fd))
            best = min(best, rel)
            cells.append(f'{fd:12.6f} ({rel:7.1e})')
        ok &= best < 3e-3
        print(f'   {name:>14} {analytic:12.6f} | ' + ' '.join(f'{c:>22}' for c in cells))

    shutil.rmtree(work)
    print('\n' + ('PASS: the analytic gradient is the real one.' if ok else
                  'FAIL: see the parts above.'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
