# ConDiv core force-field training

Contrastive-divergence training of the Upside core force field with **ff2.1's own workflow**: Peng
et al., JCTC 2022, SI "Parameterization by Contrastive Divergence". This directory holds the
**infrastructure only**, no training data. `.gitignore` keeps `training/*` out of the repo and
re-includes just these files, so a run directory created here stays untracked.

## Files

| file | what it is |
|---|---|
| `ConDiv.py` | the FF2 dual-target trainer, adapted from O. Kleinmann's Python 3 port of Peng's code (`/project2/trsosnic/okleinmann/condiv/condiv2.py`); its docstring lists every difference and why |
| `rama_basin.py` | the Ramachandran basin offsets: one offset set per trained (central, direction, neighbour) map on its fixed NDRD base (the glycine-centred and pre-proline maps), the library writer, per-residue basin populations and the once-per-epoch update |
| `verify_rama_basin.py` | gate: the trained set is as specified, every untrained map stays NDRD, and every offset reaches exactly the residues that read its map, through `upside_config` |
| `check_converged.py` | has a run updated every file, is every group at a fixed point, has it plateaued? |
| `train_chain.sbatch` | self-chaining Slurm job; submits `<run>/after_training.sbatch` when the target is reached |
| `extract_ff.py` | a checkpoint -> the six parameter files, through the run's own `expand_param` |
| `patch_glpg.py` | patch a force field into a glpG hybrid seed without rebuilding it |
| `validate_ff.sh` | release a trained force field and submit the Peng benchmark and glpG validation |
| `convergence_gate.py` | exact sign-flip test of every trained group over the last epoch: exit 0 converged, 3 not |
| `gate_or_continue.sh` | run by a run's `after_training.sbatch`: gate, then stop for review, or train one more epoch |
| `env.sh` | per-cluster Python with identical package versions (midway2: this tree's `.venv`; midway3: the shared /beagle3 venv; locally the repo `.venv`), always this tree's `py/` and `obj/`; finds `PROJECT_ROOT` from its own location |
| `pdb_list` | the 456-protein training-set manifest (a list, not data) |

## What a run directory needs

`training/<name>/` with:

```
init_param/     environment.h5, bb_env.dat, sidechain.h5, hbond.h5, sheet   (parameters/ff_2.1)
upside_input/   per protein: <code>.fasta, <code>.initial.pkl, <code>.chi
                plus rama.dat (the library the run reads) and rama_reference.pkl
pdb_list        copy from here
env.sh          copy of this directory's env.sh, adjusted if the tree differs
slurm.args      the cluster's sbatch flags, given on the command line of every submission:
                midway3  --partition=caslake
                midway2  --partition=broadwl --exclude=<the nodes listed in train_chain.sbatch>
```

`upside_input/` is ~265 MB and is **not** in the repo. Hardlink it from an existing run
(`cp -al`) rather than copying. For training from ff2.1, `upside_input/rama.dat` is
`parameters/common/rama.dat`.

## Running

```bash
cd training/myrun && source env.sh
python3 ../verify_rama_basin.py .                 # must PASS before anything is trained
python3 ../ConDiv.py initialize init_param upside_input pdb_list run_output
sbatch $(cat slurm.args) ../train_chain.sbatch . 76   # 4 epochs of 19 minibatches, self-chaining
python3 ../check_converged.py .
```

**What happens at the target.** `train_chain.sbatch` submits `<run>/after_training.sbatch`, which
calls `gate_or_continue.sh <run> <ff_name> <max_epochs>`: `convergence_gate.py` judges the last
full epoch. A converged run is released and validated by `validate_ff.sh` (the Peng benchmark and
the glpG chains, submitted to midway2's broadwl, so the run must be on midway2); an unconverged one
is trained one more epoch and judged again, up to `<max_epochs>`, after which it stops for review.
A failure of the gate itself stops everything. The rama offsets' training and held-out mismatch is
in `run_output/rama_rounds.txt`, one line per epoch.

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
  ensemble, added to every map the residue reads, for the rama offsets (no DSE term).

## What it trains

Exactly ff2.1's set: `rot` (pair, coverage and hydrophobe interactions), the sigmoid burial
`scale`, `center`, `sharpness` for 20 types and the 400 weights, the backbone term's `scale`, the
three secondary-structure H-bond energies and the second-H-bond term, and the 20 sheet mixing
energies (by central differences per residue type present).

**Not trained, as in ff2.1's own training**, because the engine returns no derivative: the
backbone term's `center`, `sharpness` and `hbond_weight` (commented out in
`BackboneSigmoidCoupling::get_param_deriv`, in master too) and `hbond.h5` entries 4-11, the rama
boundaries and sharpnesses. Their learning rates are 0 so `check_converged.py` does not list them.

**Additionally, the Ramachandran basin offsets.** A trained coil map of the library, k = (central,
direction, neighbour), keeps its NDRD values as a fixed base and gets offsets on some of its
basins, after which the map is renormalised as the NDRD maps are. So an offset is a weight factor
on its basin: the shape inside each basin is the base's, and only the basins' depths, their
frequencies, are trained. The basins partition the torus, with 13 deg edges: alpha_R, alpha_L,
beta, pPII and `other` (phi > 0 outside alpha_L), and for a central glycine `other` split into its
mirror halves beta' and pPII'; a trained map's untrained basins together are its reference. Only
the maps where the literature and ff2.1's own error both show a local defect are trained, 158
offsets on 60 maps: GLY|X (alpha_R, alpha_L, beta), GLY|GLY (helix and beta) and X|right|PRO
(alpha_R, beta; the ring's C-delta clashes with residue i in alpha_R). The other 780 maps stay at
NDRD, because their neighbour effects are small and 456 proteins cannot resolve them per pair
(findings 1.12-1.13). Each map is its own parameter set: nothing is tied or pooled across maps, and
a map's offsets act only on the residues that read it. GLY|GLY's base is
symmetrised, as the mean of its probabilities and its mirror's, and each of its offsets is held
equal to its mirror basin's; its sheet entry, which `upside_config` mixes into the same residues,
is symmetrised the same way, so a glycine between glycines gets an exactly mirror-symmetric map.
At the end of every epoch each offset takes a damped Newton step on the MAP objective: the native basin counts
over every residue reading the map under the model's populations, with a Gaussian prior of width 1
nat on the offset. Where a basin holds many residues this is `0.5 * T0 * ln(p_free / p_native)`;
where it holds almost none the prior bounds the step, and an offset with no evidence decays to 0. There is no DSE term on the offsets: the SARW
replica keeps the rama term, so it would compare the unfolded ensemble with the map itself. 46
proteins (10%) are held out of the update and reported beside it each round in
`run_output/rama_rounds.txt`. **Run `verify_rama_basin.py` after touching any of it.**

## Traps

* **`--ntasks` must equal the minibatch size (24) and `--cpus-per-task` the 14 systems.**
* **The local Mac binary traps (SIGTRAP) at exit whenever Monte Carlo pivot moves are on**, after
  every frame completes, so every worker reports `RUN_FAIL` locally. Test workers on midway2, whose
  binary runs them cleanly.
* **Opening a `.up` with PyTables before constructing `ue.Upside` makes the engine fail to
  initialise.** Construct the engine first, then read arrays.
* **A worker that `srun` never starts is relaunched; one that ran and failed fails the step.** 24
  steps issued at once sometimes lose a few to `Task launch ... failed: Job credential expired` on
  healthy nodes. `srun`'s own messages go to `<code>.srun`; a `Task launch` failure there relaunches
  that worker at once, up to twice, and the link log says `never started ..., relaunching`.
* **`run_output/ConDiv.py` must exist.** Only `initialize` makes it, so a hand-made `run_output/`
  leaves every worker dying with `can't open file '.../run_output/ConDiv.py'`. The per-worker
  reason is in `run_output/epoch_*/<code>.output_worker`, never in the Slurm log.
* **A rama library is 35 MB, so it is never kept per minibatch.** One is written per epoch,
  `run_output/rama_round_NN.dat`, for that epoch's workers; the offsets themselves are in every
  checkpoint.
* **Do not judge convergence by parameter movement.** Adam's steps are scale-invariant. Read the
  raw gradients, which `check_converged.py` recovers from the Adam accumulators.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`; jobs die instantly
  with `ExitCode 0:53` and no log.

## `verify_rama_basin.py`, and why it is not optional

```bash
python3 verify_rama_basin.py <training_dir> [protein_code]
```

It checks the whole path, library file -> `upside_config` -> the per-residue maps in the `.up`
file: the basins are continuous across phi = +-180, exactly mirror-symmetric and partition the
torus, and the trained set is exactly the 60 maps above; one map's offset changes that map's coil
entry and nothing else, under random offsets every untrained map and the sheet group but GLY|GLY
are bitwise the source, every written map is normalised, and GLY|GLY's coil and sheet entries stay
exactly symmetric; for a GLY|X, a GLY|GLY and a pre-proline map and for the first and last
residues' maps, perturbing the offsets changes exactly the residues `residue_keys` says read that
map, raising each inside the basin it was raised in; a residue that reads only untrained maps keeps
its NDRD map exactly; and a glycine that reads only GLY|GLY maps gets an exactly symmetric map from
the coil/sheet mixture. **An indexing slip fails silently**: it trains one
pair's offsets on another pair's residues and never raises anything.

## Reading `check_converged.py`

The sound fixed-point statistic is `||mean g|| / mean|g|` against `1/sqrt(n)`. **The
pairwise-cosine t-statistic is anti-conservative**, because it treats the `n(n-1)/2` pairs as
independent when they share vectors, so do not let it carry a conclusion. A systematic drift shows
as *positive* cosine; negative means oscillation about a minimum. Scalar gradients are
heavy-tailed, and nothing should be judged on less than one epoch.
