# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field with **ff2.1's own workflow**: Peng
et al., JCTC 2022, SI "Parameterization by Contrastive Divergence". This directory holds the
**infrastructure only**, no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the FF2 dual-target trainer, adapted from O. Kleinmann's Python 3 port of Peng's code (`/project2/trsosnic/okleinmann/condiv/condiv2.py`); its docstring lists every difference and why |
| `rama_basin.py` | the Ramachandran basins and per-residue basin populations, recorded every step as a diagnostic (not a parameter) |
| `build_gly_library.py` | builds the library ff3.0 trains with: the central-glycine row replaced by the AWH-measured free energy of capped glycine dipeptides, every other row unchanged; checks itself through `upside_config` |
| `check_converged.py` | has a run updated every file, is every group at a fixed point, has it plateaued? |
| `train_chain.sbatch` | self-chaining Slurm job; submits `<run>/after_training.sbatch` when the target is reached |
| `extract_ff.py` | a checkpoint -> the six parameter files, through the run's own `expand_param` |
| `convergence_gate.py` | exact sign-flip test of every trained group over the last epoch: exit 0 converged, 3 not |
| `gate_or_continue.sh` | run by a run's `after_training.sbatch`: gate, then stop for review, or train one more epoch |
| `env.sh` | the Python (midway2: this tree's `.venv` with its modules; locally the repo `.venv`), always this tree's `py/` and `obj/`; finds `PROJECT_ROOT` from its own location |
| `pdb_list` | the 456-protein training-set manifest (a list, not data) |

## What a run directory needs

`training/<name>/` with:

```
init_param/     environment.h5, bb_env.dat, sidechain.h5, hbond.h5, sheet   (parameters/ff_2.1)
upside_input/   per protein: <code>.fasta, <code>.initial.pkl, <code>.chi
                plus rama.dat (the fixed library the run reads) and rama_reference.pkl
pdb_list        copy from here
env.sh          copy of this directory's env.sh, adjusted if the tree differs
slurm.args      the cluster's sbatch flags, given on the command line of every submission:
                midway2  --partition=broadwl --exclude=<the nodes listed in train_chain.sbatch>
```

`upside_input/` is ~265 MB and is **not** in the repo. Hardlink it from an existing run
(`cp -al`) rather than copying, then replace `rama.dat` by a fresh copy (never edit a hardlinked
file in place). For ff2.1's own library it is `parameters/common/rama.dat`; for ff3.0 it is
`parameters/common/rama31.dat`, the output of `build_gly_library.py`.

## Running

```bash
cd training/myrun && source env.sh
python3 ../ConDiv.py initialize init_param upside_input pdb_list run_output
sbatch $(cat slurm.args) ../train_chain.sbatch . 76   # 4 epochs of 19 minibatches, self-chaining
python3 ../check_converged.py .
```

**What happens at the target.** `train_chain.sbatch` submits `<run>/after_training.sbatch`, which
calls `gate_or_continue.sh <run> <ff_name> <max_epochs>`: `convergence_gate.py` judges the last
full epoch. A converged run stops and prints the `extract_ff.py` command that writes its newest
checkpoint to `parameters/<ff_name>`; an unconverged one is trained one more epoch and judged
again, up to `<max_epochs>`, after which it stops for review. A failure of the gate itself stops
everything.

`initialize` copies `ConDiv.py` and `rama_basin.py` into `run_output/`, and every later step, the
driver included, runs that copy: a run is never continued by later code.

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
* contrast = NSE + 0.3 DSE, summed over the 24 proteins of the minibatch, into Adam;
* per residue, the basin populations of the native-restrained replica and of the NSE's free
  ensemble, written to `<code>.divergence.pkl` (`rama_native`, `rama_free`) for diagnosis.

## What it trains

Exactly ff2.1's set: `rot` (pair, coverage and hydrophobe interactions), the sigmoid burial
`scale`, `center`, `sharpness` for 20 types and the 400 weights, the backbone term's `scale`, the
three secondary-structure H-bond energies and the second-H-bond term, and the 20 sheet mixing
energies (by central differences per residue type present).

**Not trained, as in ff2.1's own training**, because the engine returns no derivative: the
backbone term's `center`, `sharpness` and `hbond_weight` (commented out in
`BackboneSigmoidCoupling::get_param_deriv`, in master too) and `hbond.h5` entries 4-11, the rama
boundaries and sharpnesses. Their learning rates are 0 so `check_converged.py` does not list them.

**The Ramachandran library is fixed, not trained.** ff3.0's differs from ff2.1's only in the
central-glycine row, which `build_gly_library.py` takes from AWH on capped glycine dipeptides: the
PDB row is part local energy and part evolutionary placement (glycine is put where a fold needs a
left-handed residue), Upside applies a map as pure energy, and a map trained against native
structures relearns the placement (findings 1.15-1.16). The non-local terms above are trained around
it, so they must place glycine where a fold needs it.

## Traps

* **`--ntasks` must equal the minibatch size (24) and `--cpus-per-task` the 14 systems.**
* **Running without Slurm** (a local machine), the driver starts at most `CONDIV_LOCAL_WORKERS`
  workers at a time (default 1). One worker keeps ~10 cores busy, so an M1 Ultra runs two:
  `CONDIV_LOCAL_WORKERS=2 nohup caffeinate -i python run_output/ConDiv.py restart <checkpoint> <n_steps> > train_local.log 2>&1 &`.
  The Mac binary no longer traps at exit with Monte Carlo moves on (findings 3.9): a full worker
  ran there cleanly on 2026-10-01 (1ga3, 345 s).
* **Opening a `.up` with PyTables before constructing `ue.Upside` makes the engine fail to
  initialise.** Construct the engine first, then read arrays.
* **A worker that `srun` never starts is relaunched; one that ran and failed fails the step.** 24
  steps issued at once sometimes lose a few to `Task launch ... failed: Job credential expired` on
  healthy nodes. `srun`'s own messages go to `<code>.srun`; a `Task launch` failure there relaunches
  that worker at once, up to twice, and the link log says `never started ..., relaunching`.
* **`run_output/ConDiv.py` must exist.** Only `initialize` makes it, so a hand-made `run_output/`
  leaves every worker dying with `can't open file '.../run_output/ConDiv.py'`. The per-worker
  reason is in `run_output/epoch_*/<code>.output_worker`, never in the Slurm log.
* **Do not judge convergence by parameter movement.** Adam's steps are scale-invariant. Read the
  raw gradients, which `check_converged.py` recovers from the Adam accumulators.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log.

## Reading `check_converged.py`

The sound fixed-point statistic is `||mean g|| / mean|g|` against `1/sqrt(n)`. **The
pairwise-cosine t-statistic is anti-conservative**, because it treats the `n(n-1)/2` pairs as
independent when they share vectors, so do not let it carry a conclusion. A systematic drift shows
as *positive* cosine; negative means oscillation about a minimum. Scalar gradients are
heavy-tailed, and nothing should be judged on less than one epoch.
