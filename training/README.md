# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field. This directory holds the
**infrastructure only** — no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the trainer |
| `rama_gly_gradient.py` | analytic gradient of the energy w.r.t. the glycine Ramachandran map, plus the library reader/writer. Lives here rather than in `py/`, which is shared Upside infrastructure, and is imported by `ConDiv.py` via `$PROJECT_ROOT/training` on `PYTHONPATH` (set by `env.sh`) |
| `check_converged.py` | is a given force field a stationary point of the objective? |
| `env.sh` | module/venv/PYTHONPATH setup, derives `PROJECT_ROOT` from its own location |
| `init.sh` | pack a force field into the latent vector, build minibatches |
| `train_chain.sbatch` | self-chaining Slurm job |
| `pdb_list` | the 456-protein training-set manifest (a list, not data) |

## What a run directory needs

`training/<name>/` with:

```
init_param/     environment.h5, sidechain.h5, hbond, sheet   (e.g. from parameters/ff_2.1)
upside_input/   per protein: <code>.fasta, <code>.initial.pkl, <code>.chi
                plus rama.dat and rama_reference.pkl
pdb_list        copy from here
```

`upside_input/` is ~265 MB and is **not** in the repo. Hardlink it from an existing run
(`cp -al`) rather than copying, then replace `rama.dat` if the run uses a different library.
`init_param/hbond` is `parameters/ff_X/hbond.h5` renamed; the others keep their names.

`upside_input/rama.dat` is whichever library the run trains against, copied under that name
(ff3.1 uses `parameters/common/rama31.dat`). `rama_reference.pkl` is a single 72x72 density
applied to **every** residue, so it does not depend on the library's per-residue maps and is not
regenerated when one changes.

## Running

```bash
bash init.sh training/myrun                       # check the printed pack_param residual
sbatch train_chain.sbatch training/myrun 500      # 500 minibatches, self-chaining
python3 check_converged.py training/myrun/run_output
```

Progress is `find run_output -name checkpoint.pkl -path '*epoch_*' | wc -l` against the target.
38 minibatches per epoch, 12 proteins each.

## What this trainer trains

`rot`, `env`, `hb`, `sheet` and `gly`. The first four are what the Theano original trained; `gly`
is new. Two of the original four needed remapping because their nodes changed shape since:

* **`hb`** is a multiplicative scale on `hbond_energy.parameters[:4]`. The potential is exactly
  linear in those four (measured to 7e-7 against the engine), so `dE/ds = E/s`, which is the
  original's formula unchanged.
* **`sheet`** is a single common offset on all 20 per-residue-type mixing energies, differenced
  against the `more_/less_sheet_rama_pot_ALL` pair that `write_rama_map_pot` emits under
  `--rama-param-deriv`. Training all 20 types separately would cost ~41x the divergence instead
  of 3x, because each finite-difference direction is two extra passes over all 250 frames.

* **`gly`** is the central-glycine coil Ramachandran map, as a pair `(S, A)` of 72x72 arrays that
  are mirror-symmetric and mirror-antisymmetric under `(phi,psi) -> (-phi,-psi)`. The library's
  glycine coil row is written from them as `X|GLY = S + A` and `GLY|GLY = S`, so a glycine
  flanked by glycines stays achiral by construction while every other context is free. Training
  starts at `A = 0` on a symmetrised row, so the handedness is learned rather than assumed.

  **One map serves every neighbour.** Per-neighbour maps would be 40 x 5,184 parameters against
  ~30,000 glycine samples per minibatch, and the AWH measurement cannot resolve per-pair structure
  either (S/N 1.48). So the row loses the library's neighbour dependence at step 0. Whether that
  is acceptable is exactly what the 40-context AWH campaign is measuring.

  Its gradient is **analytic, not finite-differenced**: `rama_map_pot` is a periodic interpolating
  bicubic spline built from a tensor product of 1D solves, so the map enters the energy linearly
  and separably and `dE/d(map[i,j])` is a spline-smoothed 2D histogram of the glycine `(phi,psi)`
  samples. Finite differencing 5,184 values would cost 10,369x a divergence. The chain rule from
  the per-residue map back to `(S, A)` runs through four log-sum-exp stages and is done by torch
  autograd rather than by hand. **Run `verify_gly_gradient.py` after touching any of it** (see
  below).

Sheet training is **effectively free**: measured 625-705 s/minibatch with it on, against 650-712
historically with it off. A step is ~11 min on 96 CPUs, so 500 minibatches is ~3.8 days.

## Traps

* **`--ntasks` must equal the minibatch size (12).** Each task is one protein worker.
* **The env parameter vector is `coeff` then `weights` (360 + 400 = 760).** Requesting
  `coeff.shape` alone fails with `Wrong number of parameters, expected 760 but got 360` and kills
  every worker. A 2026-09-10 `libupside.so` rebuild introduced this and broke all training for
  nine days before anyone ran a job; `compute_divergence` now requests the full vector and slices.
* **Opening a `.up` with PyTables before constructing `ue.Upside` makes the engine fail to
  initialize.** Construct the engine first, then read arrays.
* **`run_output/ConDiv.py` must exist.** Only `initialize` makes that copy
  (`state['worker_path']`), so restarting from an existing `initial_checkpoint.pkl` in a
  hand-made directory leaves every worker dying with `can't open file '.../run_output/ConDiv.py'`
  and the link failing on `All jobs failed`. The real reason is in
  `run_output/epoch_*/<code>.output_worker`, never in the Slurm log, which shows only `exit code 2`.
* **A rama library is 35 MB, so it is never kept per minibatch.** `expand_param` writes one for
  the workers to share and `run_minibatch` deletes it once they exit; keeping one per step would
  add 35 GB over a 500-step run. The parameters themselves live in `param.gly` in every
  checkpoint, so the library is regenerated with `rgg.write_gly_library` when it is wanted.
* **Do not judge convergence by parameter movement.** Adam's first step is scale-invariant: the
  update is `-alpha*g/(|g|+eps)`, so parameters move by ~`alpha` regardless of gradient size.
  Read the raw gradients, which `check_converged.py` recovers from the Adam accumulators.
* **A hand-made `run_output/` needs `ConDiv.py` copied into it.** Only `initialize` makes that
  copy (`state['worker_path']`), so restarting from an existing `initial_checkpoint.pkl` in a
  fresh directory leaves every worker exiting with
  `can't open file '.../run_output/ConDiv.py'` and the link dying on `All jobs failed`. The
  per-worker reason is in `run_output/epoch_*/<code>.output_worker`, not in the Slurm log, which
  shows only `exit code 2`.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log, because they cannot create the output file.

## `verify_gly_gradient.py`, and why it is not optional

```bash
python3 verify_gly_gradient.py <training_dir> [protein_code]
```

It checks the analytic `gly` gradient against finite differences taken through the entire
pipeline: library file -> `upside_config` -> engine. **An analytic gradient fails silently** — a
missing stage or a flipped softmax factor trains steadily in the wrong direction for days without
raising anything. On its first run it failed at 37% and found a real omission, the per-map
Boltzmann shift that `write_rama_map_pot` applies last (`up.md` 2.8a).

**Read the eps sweep, not a single column.** The library stores `dimer_pot` as float32, so small
eps is rounding noise and large eps picks up real curvature from the log-sum-exp mixtures; the
agreement is the minimum of the bowl. At the optimum it is 3.8e-5. Judging by `eps=1e-3` alone
would have called a correct gradient a 1-4% failure.

## Reading `check_converged.py`

The sound statistic is `||mean g|| / mean|g|` against `1/sqrt(n)`, tested exactly by sign-flipping
the per-step gradients. **The pairwise-cosine t-statistic it also prints is anti-conservative** —
it treats the `n(n-1)/2` pairs as independent when they share vectors — so do not let it carry a
conclusion. And a systematic drift shows as *positive* cosine; negative means oscillation about a
minimum.

Expect scalar gradients (`hb`, `sheet`) to be heavy-tailed. A signal at n=3 can vanish by n=6; do
not read a partial run.
