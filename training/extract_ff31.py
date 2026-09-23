"""Write a complete force field from an ff3.1 ConDiv checkpoint.

`training/gly-sym/extract_ff.py` writes only `sidechain.h5` and `environment.h5`, which was
enough when `hb`, `sheet` and the glycine map were frozen. ff3.1 trains all five, so extracting
only two would silently validate a force field that is part trained and part ff_2.1.

This calls ConDiv's own `expand_param`, so the output is byte-identical to what training writes
each minibatch -- there is no second implementation to drift.

    python3 extract_ff31.py <checkpoint.pkl> <output_dir>

Produces, in output_dir:
    sidechain.h5      rot
    environment.h5    env
    hbond.h5          hb, as a scale on parameter[:4]
    sheet             sheet, as a common offset
    rama.dat          the trained glycine coil row (S, A); everything else is ff_2.1's
"""

import os
import pickle as cp
import shutil
import sys

import numpy as np


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    ckpt, out = os.path.abspath(sys.argv[1]), sys.argv[2]
    os.makedirs(out, exist_ok=True)

    # Import the ConDiv that PRODUCED this checkpoint, not whatever copy is nearest. `initialize`
    # writes one into run_output alongside the checkpoints precisely so a result can be re-read
    # with the code that made it, and the pickles reference its classes.
    run_output = os.path.dirname(os.path.dirname(ckpt))
    sys.path.insert(0, os.path.join(os.environ['UPSIDE_HOME'], 'py'))
    sys.path.insert(0, os.path.join(os.environ['UPSIDE_HOME'], 'training'))
    sys.path.insert(0, run_output)
    global ConDiv, rgg
    import ConDiv
    import rama_gly_gradient as rgg
    # ConDiv checkpoints were pickled from __main__, so its classes must resolve there.
    import __main__
    for _n in dir(ConDiv):
        if _n[0].isupper():
            setattr(__main__, _n, getattr(ConDiv, _n))

    with open(ckpt, 'rb') as f:
        state = cp.load(f)
    param = state['param']
    orig = dict(state['init_param_files'])

    # init_param_files may hold paths relative to the training directory
    base = os.path.dirname(run_output)
    for k, v in orig.items():
        if not os.path.isabs(v):
            orig[k] = os.path.join(base, v)
        if not os.path.exists(orig[k]):
            sys.exit(f'initial parameter file missing: {k} -> {orig[k]}')

    new = dict(
        rot=os.path.join(out, 'sidechain.h5'),
        env=os.path.join(out, 'environment.h5'),
        hb=os.path.join(out, 'hbond.h5'),
        sheet=os.path.join(out, 'sheet'),
        rama=os.path.join(out, 'rama.dat'),
    )
    ConDiv.expand_param(param, orig, new)

    S, A = param.gly
    print(f'checkpoint : {ckpt}')
    print(f'  step_num {state["solver"].step_num}, epoch {state["epoch"]}, i_mb {state["i_mb"]}')
    print(f'  hb    {param.hb:.6f}')
    print(f'  sheet {param.sheet:+.6f}')
    print(f'  gly   GLY|GLY asymmetry {np.abs(S - rgg.mirror(S)).max():.2e}  (must be 0)')
    print(f'\nwritten to {out}:')
    for k, v in sorted(new.items()):
        print(f'  {k:6s} {os.path.basename(v):16s} {os.path.getsize(v) / 1e6:8.2f} MB')

    # A benchmark run must not silently fall back to ff_2.1 for anything.
    missing = [k for k, v in new.items() if not os.path.exists(v)]
    if missing:
        sys.exit(f'FAILED to write: {missing}')
    if np.abs(S - rgg.mirror(S)).max() > 1e-9:
        sys.exit('FAILED: GLY|GLY is not symmetric')


if __name__ == '__main__':
    main()
