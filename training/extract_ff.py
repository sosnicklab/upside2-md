"""Write a complete force field from a ConDiv checkpoint.

    python3 extract_ff.py <checkpoint.pkl> <output_dir>

Calls the producing run's own `expand_param`, so the files are byte-identical to what training
writes each step and there is no second implementation to drift. Produces, in output_dir:

    sidechain.h5    rot (pair, coverage, hydrophobe interactions and placements)
    environment.h5  sigmoid burial scale, center, sharpness and the 400 weights
    bb_env.dat      backbone desolvation term
    hbond.h5        the twelve H-bond parameters
    sheet           the 20 sheet mixing energies
    rama.dat        the trained library if the run trained the glycine row, otherwise a copy of
                    the library it read
"""

import os
import pickle as cp
import shutil
import sys

import numpy as np


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    ckpt, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    os.makedirs(out, exist_ok=True)

    # Import the ConDiv (and rama_gly_gradient) that PRODUCED this checkpoint: `initialize` puts
    # both in run_output, and the pickles reference the trainer's classes.
    run_output = os.path.dirname(os.path.dirname(ckpt))
    sys.path.insert(0, run_output)
    import ConDiv
    import __main__
    for n in dir(ConDiv):
        if n[0].isupper():
            setattr(__main__, n, getattr(ConDiv, n))
    rgg = ConDiv.rgg

    state = cp.load(open(ckpt, 'rb'))
    param, init, row = state['param'], state['init_param_files'], state['row']
    for k, v in init.items():
        if not os.path.exists(v):
            sys.exit(f'initial parameter file missing: {k} -> {v}')

    new = dict(rot=os.path.join(out, 'sidechain.h5'), env=os.path.join(out, 'environment.h5'),
               bbenv=os.path.join(out, 'bb_env.dat'), hb=os.path.join(out, 'hbond.h5'),
               sheet=os.path.join(out, 'sheet'), rama=os.path.join(out, 'rama.dat'))
    ConDiv.expand_param(param, init, new, row)
    if param.gly is None:
        shutil.copyfile(init['rama'], new['rama'])

    print(f'checkpoint {ckpt}')
    print(f'  step {state["solver"].step_num}, next epoch {state["epoch"]} minibatch {state["i_mb"]}')
    print(f'  hb {np.array2string(param.hb[:3], precision=4)}  dhb {param.dhb[0]:.4f}'
          f'  bb scale {param.bbenve:.4f}  sheet mean {np.mean(param.sheet):.4f}')
    if param.gly is not None:
        gg = param.gly[row['is_gg']]
        asym = np.abs(gg - rgg.mirror(gg)).max()
        print(f'  gly X|GLY dG(aR->aL) mean {rgg.handedness(param.gly)[~row["is_gg"]].mean():+.4f}'
              f'  GLY|GLY asymmetry {asym:.2e}  (must be 0)')
        if asym > 1e-9:
            sys.exit('FAILED: GLY|GLY is not mirror-symmetric')
    else:
        print('  glycine row not trained; rama.dat is the library the run read')

    missing = [k for k, v in new.items() if not os.path.getsize(v)]
    if missing:
        sys.exit(f'FAILED to write: {missing}')
    print(f'\nwritten to {out}:')
    for k, v in sorted(new.items()):
        print(f'  {k:6s} {os.path.basename(v):16s} {os.path.getsize(v) / 1e6:8.2f} MB')


if __name__ == '__main__':
    main()
