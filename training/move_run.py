"""Point a ConDiv checkpoint at a run directory on another machine.

    python3 move_run.py <checkpoint.pkl> <old_run_dir> <new_run_dir> <out_checkpoint.pkl>

A checkpoint records absolute paths: the run's base and initial-parameter directories, the
trainer copy the workers execute, every initial parameter file including the rama library, and
each protein's fasta, native and chi files. Training continued elsewhere (local Mac to midway2, or
back) needs those paths under the new run directory, which must hold the same `init_param/`,
`upside_input/` and the copied `run_output/`. Every string in the state that starts with the old
run directory is rewritten; the script fails if any path still names the old one, or if a rewritten
path does not exist where it is run, so run it on the machine that continues the training.
"""

import collections
import os
import pickle as cp
import sys


def rewrite(x, old, new):
    """Every str under x that starts with old, with that prefix replaced; containers rebuilt."""
    if isinstance(x, str):
        return new + x[len(old):] if x.startswith(old) else x
    if isinstance(x, tuple) and hasattr(x, '_fields'):
        return type(x)(*[rewrite(v, old, new) for v in x])
    if isinstance(x, (list, tuple)):
        return type(x)(rewrite(v, old, new) for v in x)
    if isinstance(x, set):
        return {rewrite(v, old, new) for v in x}
    if isinstance(x, dict):
        return type(x)((k, rewrite(v, old, new)) for k, v in x.items())
    return x


def strings(x):
    if isinstance(x, str):
        yield x
    elif isinstance(x, (list, tuple, set)):
        for v in x:
            yield from strings(v)
    elif isinstance(x, dict):
        for v in x.values():
            yield from strings(v)


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    ckpt, old, new, out = sys.argv[1:]
    old, new = old.rstrip('/'), os.path.abspath(new)

    # the producing trainer's classes, from the copy that travels with run_output
    run_output = os.path.dirname(os.path.abspath(ckpt))
    while not os.path.exists(os.path.join(run_output, 'ConDiv.py')):
        if run_output == os.path.dirname(run_output):
            sys.exit(f'no ConDiv.py above {ckpt}')
        run_output = os.path.dirname(run_output)
    sys.path.insert(0, run_output)
    import ConDiv
    import __main__
    for n in dir(ConDiv):
        if n[0].isupper():
            setattr(__main__, n, getattr(ConDiv, n))

    state = cp.load(open(ckpt, 'rb'))
    moved = collections.OrderedDict((k, rewrite(v, old, new)) for k, v in state.items())

    # every absolute path a checkpoint holds lies under its run directory
    stray = [s for s in strings(dict(moved)) if s.startswith('/') and not s.startswith(new + '/')]
    if stray:
        sys.exit(f'FAILED: {len(stray)} paths not under {new} (is {old} the old run dir?), '
                 f'e.g. {stray[0]}')
    paths = [s for s in strings(dict(moved)) if s.startswith(new + '/')]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit(f'FAILED: {len(missing)} rewritten paths do not exist here, e.g. {missing[0]}')
    with open(out, 'wb') as f:
        cp.dump(dict(moved), f, -1)
    print(f'{ckpt}\n  {len(paths)} paths moved from {old} to {new}\n  written {out}')


if __name__ == '__main__':
    main()
