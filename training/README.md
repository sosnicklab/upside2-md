# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field. This directory holds the
**infrastructure only** — no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the trainer |
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

## Running

```bash
bash init.sh training/myrun                       # check the printed pack_param residual
sbatch train_chain.sbatch training/myrun 500      # 500 minibatches, self-chaining
python3 check_converged.py training/myrun/run_output
```

Progress is `find run_output -name checkpoint.pkl -path '*epoch_*' | wc -l` against the target.
38 minibatches per epoch, 12 proteins each.

## What this trainer trains

`rot`, `env`, `hb` and `sheet` — the same four the Theano original trained. Two of them needed
remapping because their nodes changed shape since:

* **`hb`** is a multiplicative scale on `hbond_energy.parameters[:4]`. The potential is exactly
  linear in those four (measured to 7e-7 against the engine), so `dE/ds = E/s`, which is the
  original's formula unchanged.
* **`sheet`** is a single common offset on all 20 per-residue-type mixing energies, differenced
  against the `more_/less_sheet_rama_pot_ALL` pair that `write_rama_map_pot` emits under
  `--rama-param-deriv`. Training all 20 types separately would cost ~41x the divergence instead
  of 3x, because each finite-difference direction is two extra passes over all 250 frames.

Sheet training is what makes a step ~26 min instead of ~9. Dropping it roughly triples throughput.

## Traps

* **`--ntasks` must equal the minibatch size (12).** Each task is one protein worker.
* **The env parameter vector is `coeff` then `weights` (360 + 400 = 760).** Requesting
  `coeff.shape` alone fails with `Wrong number of parameters, expected 760 but got 360` and kills
  every worker. A 2026-09-10 `libupside.so` rebuild introduced this and broke all training for
  nine days before anyone ran a job; `compute_divergence` now requests the full vector and slices.
* **Opening a `.up` with PyTables before constructing `ue.Upside` makes the engine fail to
  initialize.** Construct the engine first, then read arrays.
* **Do not judge convergence by parameter movement.** Adam's first step is scale-invariant: the
  update is `-alpha*g/(|g|+eps)`, so parameters move by ~`alpha` regardless of gradient size.
  Read the raw gradients, which `check_converged.py` recovers from the Adam accumulators.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log, because they cannot create the output file.

## Reading `check_converged.py`

The sound statistic is `||mean g|| / mean|g|` against `1/sqrt(n)`, tested exactly by sign-flipping
the per-step gradients. **The pairwise-cosine t-statistic it also prints is anti-conservative** —
it treats the `n(n-1)/2` pairs as independent when they share vectors — so do not let it carry a
conclusion. And a systematic drift shows as *positive* cosine; negative means oscillation about a
minimum.

Expect scalar gradients (`hb`, `sheet`) to be heavy-tailed. A signal at n=3 can vanish by n=6; do
not read a partial run.
