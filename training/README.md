# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field with **ff2.1's own workflow**: Peng
et al., JCTC 2022, SI "Parameterization by Contrastive Divergence". This directory holds the
**infrastructure only**, no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the FF2 dual-target trainer, adapted from O. Kleinmann's Python 3 port of Peng's code (github `nnamnielk/condiv4upside2`); its docstring lists every difference and why. Commands: `initialize`, `restart`, `gate` (convergence test), `extract` (a checkpoint's force-field files) |
| `train_chain.sbatch` | self-chaining Slurm job; at the target it runs the gate and either stops or trains one more epoch |
| `env.sh` | this tree's `.venv`, `py/` and `obj/` (load any site modules first); finds `PROJECT_ROOT` from its own location |
| `pdb_list` | the 456-protein training-set manifest (a list, not data) |

## What a run directory needs

`training/<name>/` with:

```
init_param/     environment.h5, bb_env.dat, sidechain.h5, hbond.h5, sheet   (parameters/ff_2.1;
                for ff3.0 hbond.h5 carries glycine's three offsets, see What it trains)
upside_input/   per protein: <code>.fasta, <code>.initial.pkl, <code>.chi
                plus rama.dat (the fixed library the run reads) and rama_reference.pkl
pdb_list        copy from here
env.sh          copy of this directory's env.sh, adjusted if the tree differs
slurm.args      the cluster's sbatch flags, given on the command line of every submission:
                --account=<account> --partition=<partition> --exclude=<nodes>
```

`upside_input/` is ~265 MB and is **not** in the repo. Hardlink it from an existing run
(`cp -al`) rather than copying, then replace `rama.dat` by a fresh copy (never edit a hardlinked
file in place). For ff2.1's own library it is `parameters/common/rama.dat`; for ff3.0 it is
`parameters/common/rama31.dat`.

## Running

```bash
cd training/myrun && source env.sh
python3 ../ConDiv.py initialize init_param upside_input pdb_list run_output
sbatch $(cat slurm.args) ../train_chain.sbatch . 76 13   # 4 epochs of 19 minibatches, at most 13
```

`initialize` copies `ConDiv.py` into `run_output/`, and every later command runs that copy
(`python3 run_output/ConDiv.py gate .`, `... extract <checkpoint> <out_dir>`), so a run is never
continued, judged or extracted by later code.

**What happens at the target.** The link that reaches it runs `ConDiv.py gate` on the last full
epoch (report in `gate_step<N>.txt`). The gate recovers each group's raw per-step gradients from the
Adam state and asks whether their sum is unremarkable among all sign flips of the steps, an exact
permutation test; every group must pass at a family-wise 5%. Converged: the chain stops and lists
its epoch-end checkpoints. The one to release is chosen by simulating them, because the fixed-point
test does not say which iterate simulates best (the last one is only where Adam stopped), and
`ConDiv.py extract` writes it to `parameters/<ff_name>`. Not converged: the chain trains one more
epoch and judges again, up to the third argument (`<max_epochs>`; default no extension), then stops
for review. A failure of the gate itself stops the chain.

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

ff2.1's set: `rot` (pair, coverage and hydrophobe interactions), the sigmoid burial `scale`,
`center`, `sharpness` for 20 types and the 400 weights, the backbone term's `scale`, the three
secondary-structure H-bond energies and the second-H-bond term, and the 20 sheet mixing energies
(by central differences per residue type present).

ff3.0 adds one group and damps one (findings 1.17):

* **`hbg`, glycine's own offsets on the three H-bond basin energies**, `hbond.h5` entries 12-14.
  The file names the class in `class_restype` (`[b'GLY']`), and `upside_config` then writes a
  per-residue `residue_class`; a 12-entry file has neither and is unchanged. The offsets start at
  zero (an init file is ff2.1's twelve values, then `0, 0, 0`), train at `hb`'s rate, and are zeroed
  with the shared energies in the SARW replica. Without them, the shared energies can keep loop
  glycines left-handed only at the cost of every helical glycine's helix.
* **`rot` at 0.025, 10x the port's smaller.** At the port's rate Adam's normalised step
  random-walks the 31,420 pair coefficients, of which ~6% carry any signal, and the walk weakens
  helices and folds; putting the epoch-0 table back restored them.

**Not trained, as in ff2.1's own training**, because the engine returns no derivative: the
backbone term's `center`, `sharpness` and `hbond_weight` (commented out in
`BackboneSigmoidCoupling::get_param_deriv`, in master too) and `hbond.h5` entries 4-11, the rama
boundaries and sharpnesses. Their learning rates are 0, so the gate does not test them.

**The Ramachandran library is fixed, not trained.** ff3.0's differs from ff2.1's only in the
central-glycine row, taken from AWH on capped glycine dipeptides (GLY_sym.md §5): the
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
* **Do not judge convergence by parameter movement, or on less than one epoch.** Adam's steps are
  scale-invariant, so the gate reads raw gradients; scalar gradients are heavy-tailed, and windows
  shorter than an epoch have been misread as plateaus.
