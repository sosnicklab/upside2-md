# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field with **ff2.1's own workflow**: Peng
et al., JCTC 2022, SI "Parameterization by Contrastive Divergence". This directory holds the
**infrastructure only**, no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the FF2 dual-target trainer, adapted from O. Kleinmann's Python 3 port of Peng's code (`/project2/trsosnic/okleinmann/condiv/condiv2.py`); its docstring lists every difference and why |
| `rama_gly_gradient.py` | analytic gradient of the energy w.r.t. the glycine row of the Ramachandran library, plus its reader/writer and the GLY\|GLY symmetry constraint. Used only with `TRAIN_GLY = True` |
| `check_converged.py` | has a run updated every file, is every group at a fixed point, has it plateaued? |
| `train_chain.sbatch` | self-chaining Slurm job; submits `<run>/after_training.sbatch` when the target is reached |
| `extract_ff.py` | a checkpoint -> the six parameter files, through the run's own `expand_param` |
| `patch_glpg.py` | patch a force field into a glpG hybrid seed without rebuilding it |
| `validate_ff.sh` | release a trained force field and submit the Peng benchmark and glpG validation |
| `convergence_gate.py` | exact sign-flip test of every trained group over the last epoch: exit 0 converged, 3 not |
| `gate_or_continue.sh` | run by a run's `after_training.sbatch`: gate, then release and validate, or train one more epoch |
| `env.sh` | module/venv/PYTHONPATH setup, derives `PROJECT_ROOT` from its own location |
| `pdb_list` | the 456-protein training-set manifest (a list, not data) |

## What a run directory needs

`training/<name>/` with:

```
init_param/     environment.h5, bb_env.dat, sidechain.h5, hbond.h5, sheet   (parameters/ff_2.1)
upside_input/   per protein: <code>.fasta, <code>.initial.pkl, <code>.chi
                plus rama.dat (the library the run reads) and rama_reference.pkl
pdb_list        copy from here
env.sh          copy of this directory's env.sh, adjusted if the tree differs
```

`upside_input/` is ~265 MB and is **not** in the repo. Hardlink it from an existing run
(`cp -al`) rather than copying. For training from ff2.1, `upside_input/rama.dat` is
`parameters/common/rama.dat`.

## Running

```bash
cd training/myrun && source env.sh
python3 ../ConDiv.py initialize init_param upside_input pdb_list run_output
sbatch ../train_chain.sbatch . 76                 # 4 epochs of 19 minibatches, self-chaining
python3 ../check_converged.py .
```

**What happens at the target.** `train_chain.sbatch` submits `<run>/after_training.sbatch`. For
a release run that calls `gate_or_continue.sh <run> <ff_name> <max_epochs>`: `convergence_gate.py`
judges the last full epoch, and a converged run is released and validated by `validate_ff.sh`,
while an unconverged one is trained one more epoch and judged again, up to `<max_epochs>`, after
which it stops for review. A failure of the gate itself stops everything.

`initialize` copies `ConDiv.py` and `rama_gly_gradient.py` into `run_output/`, and every later
step, the driver included, runs that copy: a run is never continued by later code. `TRAIN_GLY` is
read at initialisation and kept in the run's state.

## One training step

Per protein, one replica-exchange run of 14 systems for 8000 time units, all starting from the
native: a native-restrained replica, 12 free replicas at T = 0.8 to 1.1, and a self-avoiding
random walk (SARW) with H-bond, side-chain burial and rotamer pair energies scaled to zero. From
the second half:

* **NSE** = `<dV/da>_native - <dV/da>_free`, the free ensemble being the three coldest free
  replicas each reweighted exactly to T0 and mixed 0.6/0.3/0.1;
* **DSE** = `<dV/da>_SARW - <dV/da>_unfolded`, the unfolded ensemble being the frames of the two
  replicas bracketing the Rg midpoint (the Tm estimate) with Rg above 0.67 Rg(coldest) + 0.33
  Rg(hottest). A protein with no such frames contributes NSE only, and the step log says which;
* contrast = NSE + 0.3 DSE, summed over the 24 proteins of the minibatch, into Adam.

## What it trains

Exactly ff2.1's set: `rot` (pair, coverage and hydrophobe interactions), the sigmoid burial
`scale`, `center`, `sharpness` for 20 types and the 400 weights, the backbone term's `scale`, the
three secondary-structure H-bond energies and the second-H-bond term, and the 20 sheet mixing
energies (by central differences per residue type present).

**Not trained, as in ff2.1's own training**, because the engine returns no derivative: the
backbone term's `center`, `sharpness` and `hbond_weight` (commented out in
`BackboneSigmoidCoupling::get_param_deriv`, in master too) and `hbond.h5` entries 4-11, the rama
boundaries and sharpnesses. Their learning rates are 0 so `check_converged.py` does not list them.

**With `TRAIN_GLY = True`**, additionally the central-glycine coil row: all 42 finite maps, each
its own parameter, starting from ff2.1's; the two GLY|GLY maps are projected mirror-symmetric at
the start and after every update, and each update is Fourier band-limited. The gradient is
analytic: `rama_map_pot` is a periodic interpolating bicubic spline, so the map enters the energy
linearly and `dE/d(map)` is a spline-smoothed histogram of the glycine `(phi,psi)` samples. The
chain rule back through the left/right and coil/sheet mixtures is torch autograd. **Run
`verify_gly_gradient.py` after touching any of it.**

## Traps

* **`--ntasks` must equal the minibatch size (24) and `--cpus-per-task` the 14 systems.**
* **The local Mac binary traps (SIGTRAP) at exit whenever Monte Carlo pivot moves are on**, after
  every frame completes, so every worker reports `RUN_FAIL` locally. Test workers on midway2, whose
  binary runs them cleanly.
* **Opening a `.up` with PyTables before constructing `ue.Upside` makes the engine fail to
  initialise.** Construct the engine first, then read arrays.
* **`run_output/ConDiv.py` must exist.** Only `initialize` makes it, so a hand-made `run_output/`
  leaves every worker dying with `can't open file '.../run_output/ConDiv.py'`. The per-worker
  reason is in `run_output/epoch_*/<code>.output_worker`, never in the Slurm log.
* **A rama library is 35 MB, so it is never kept per minibatch.** It is written for the workers
  and deleted once they exit; the glycine row lives in `param.gly` in every checkpoint.
* **Do not judge convergence by parameter movement.** Adam's steps are scale-invariant. Read the
  raw gradients, which `check_converged.py` recovers from the Adam accumulators.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log.

## `verify_gly_gradient.py`, and why it is not optional

```bash
python3 verify_gly_gradient.py <training_dir> [protein_code]
```

It checks the glycine gradient end to end, library file -> `upside_config` -> engine: the spline,
every glycine's reconstructed per-residue map against the `rama_pot` upside_config wrote (with all
42 maps perturbed, which pins the direction/neighbour indexing), and directional finite
differences for all maps, the GLY|GLY maps and one X|GLY map. **An analytic gradient fails
silently**: on its first run an earlier version failed at 37% and found the per-map Boltzmann
shift that `write_rama_map_pot` applies last (`up.md` 2.8a).

**Read the eps sweep, not a single column.** The library stores `dimer_pot` as float32, so small
eps is rounding noise and large eps picks up curvature from the log-sum-exp mixtures; the
agreement is the minimum of the bowl.

## Reading `check_converged.py`

The sound fixed-point statistic is `||mean g|| / mean|g|` against `1/sqrt(n)`. **The
pairwise-cosine t-statistic is anti-conservative**, because it treats the `n(n-1)/2` pairs as
independent when they share vectors, so do not let it carry a conclusion. A systematic drift shows
as *positive* cosine; negative means oscillation about a minimum. Scalar gradients are
heavy-tailed, and nothing should be judged on less than one epoch.
