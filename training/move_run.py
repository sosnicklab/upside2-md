"""Make a copied ConDiv run_output native to the run directory and machine that continue it.

    python3 move_run.py <run_output> <old_run_dir> <new_run_dir>

Run it on the machine that continues the training. <run_output> is either <new_run_dir>/run_output
itself or a staging copy whose contents are moved into <new_run_dir>/run_output afterwards; a path
under <new_run_dir>/run_output is checked in the copy being converted, every other path where it
is (the new run dir's `init_param/` and `upside_input/`). Every pickle under <run_output>
(checkpoints, solver states, divergence and rmsd files) is read and written back in place:

  * PATHS. A checkpoint records absolute paths: the base and initial-parameter directories, the
    trainer copy the workers execute, every initial parameter file including the rama library, and
    each protein's fasta, native and chi files. Every string starting with the old run directory is
    rewritten to the new one.
  * NUMPY. Pickles written under NumPy 2 name `numpy._core`, which NumPy 1 (midway2's 1.23.5) lacks;
    there the same functions live in `numpy.core`. They are read with that module name and written
    back by the NumPy that runs this, so the cluster's gate and analyses can read every step,
    including those trained elsewhere. NumPy 2 reads NumPy 1 pickles as they are.

It fails, writing nothing, if any absolute path is left outside the new run directory or a
rewritten path does not exist here.
"""

import os
import pickle
import sys

import numpy as np

NUMPY_1 = int(np.__version__.split('.')[0]) < 2


class Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if NUMPY_1 and (module == 'numpy._core' or module.startswith('numpy._core.')):
            module = 'numpy.core' + module[len('numpy._core'):]
        return super().find_class(module, name)


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
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    run_output, old, new = os.path.abspath(sys.argv[1]), sys.argv[2].rstrip('/'), os.path.abspath(sys.argv[3])
    target = os.path.join(new, 'run_output') + '/'

    def exists(p):
        """A path under the new run_output is looked for in the copy being converted."""
        return os.path.exists(os.path.join(run_output, p[len(target):]) if p.startswith(target) else p)

    # the producing trainer's classes, from the copy that travels with run_output
    sys.path.insert(0, run_output)
    import ConDiv
    import __main__
    for n in dir(ConDiv):
        if n[0].isupper():
            setattr(__main__, n, getattr(ConDiv, n))

    files = sorted(os.path.join(d, f) for d, _, fs in os.walk(run_output) for f in fs
                   if f.endswith('.pkl'))
    moved, n_paths = {}, 0
    for f in files:
        with open(f, 'rb') as fh:
            obj = Unpickler(fh).load()
        obj = rewrite(obj, old, new)
        # every absolute path a ConDiv pickle holds lies under its run directory
        stray = [s for s in strings(obj) if s.startswith('/') and not s.startswith(new + '/')]
        if stray:
            sys.exit(f'FAILED, nothing written: {f} holds {len(stray)} paths not under {new} '
                     f'(is {old} the old run dir?), e.g. {stray[0]}')
        paths = [s for s in strings(obj) if s.startswith(new + '/')]
        missing = [p for p in paths if not exists(p)]
        if missing:
            sys.exit(f'FAILED, nothing written: {f} names {len(missing)} paths that do not exist '
                     f'here, e.g. {missing[0]}')
        moved[f], n_paths = obj, n_paths + len(paths)

    for f, obj in moved.items():
        with open(f + '.moved', 'wb') as fh:
            pickle.dump(obj, fh, -1)
        os.replace(f + '.moved', f)
    print(f'{run_output}: {len(files)} pickles rewritten for numpy {np.__version__}, '
          f'{n_paths} paths moved from {old} to {new}')


if __name__ == '__main__':
    main()
