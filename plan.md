# ff3.0: ff2.1's own training workflow, a physics glycine map and glycine H-bond offsets

## Project Goal

Train ff3.0 from ff2.1 with exactly ff2.1's training workflow (Peng et al. 2022 SI), modernised
only, plus the smallest changes that fix glycine:

* the central-glycine Ramachandran row is fixed to the AWH free energy of capped glycine
  dipeptides (`parameters/common/rama31.dat`, built by `training/build_gly_library.py`); every
  other map stays NDRD, bitwise;
* glycine gets its own offsets on the three H-bond basin energies (`hbg`), trained by ConDiv from
  zero;
* the side-chain (`rot`) learning rate is 10x smaller.

The release is a chosen epoch-end checkpoint, written to `parameters/ff_3.0`.

Why: the NDRD maps are statistics over loop sites of folded proteins, so they carry evolutionary
placement (glycine is put where a fold needs alpha_L) as well as local physics, and Upside adds them
as a pure energy to every residue. Glycine is where placement reverses the handedness: ln(aR/aL) is
-1.19 in NDRD and -0.58 over all residues of the 456 training natives (findings 1.9). Data from
folded proteins cannot separate the energy from the selection (findings 1.15-1.16), so glycine's
map comes from a source where no fold selects anything, and its placement in folds must come from
the non-local terms.

## Architecture and Key Decisions

**Trainer (unchanged).** `training/ConDiv.py`, O. Kleinmann's Python 3 port of Peng's FF2 dual-target
trainer, restored to the Peng et al. 2022 SI: per protein and step one native-restrained replica, 12
free replicas at T = 0.8 to 1.1, one SARW replica; 8000 time units, second half analysed; contrast =
NSE + 0.3 * DSE, with the DSE threshold as coded (both SI statements give the same DSE, findings
1.17); 456 proteins in 19 minibatches of 24. Every difference from the port is justified in the
file's docstring.

**ff2.1's parameter set (user's choice 2026-09-24).** rot (pair, coverage, hydrophobe), the 20 x 3
sigmoid burial parameters and 400 weights, the backbone term's scale, the three H-bond energies and
the second-H-bond term, the 20 sheet values. Not trained, as in ff2.1, because the engine returns no
derivative: the backbone term's center, sharpness and hbond weight, and `hbond.h5` entries 4-11.

**The glycine row is a fixed physics map (user, 2026-10-01).** In `rama31.dat`:
* each GLY|X entry is the pooled AWH surface;
* GLY|GLY is its symmetric part only, so a glycine between glycines is exactly mirror-symmetric;
* GLY|right|PRO is its own measured surface, since it differs ~20-fold from the pool;
* each entry is stored with `rama_map_pot_ref` subtracted, so the engine applies the measured
  surface exactly, and the same entry goes in the coil and sheet groups;
* every other entry stays bitwise as in `rama.dat`. Nothing in the map is trained.

**Each map is its own conditional distribution (user, findings 1.8).** A correction, symmetry or
pattern measured on one (central, direction, neighbour) map is never transferred to another, and
pooled maps are not used to derive a per-map correction.

**Glycine's own H-bond basin offsets (user, 2026-10-02; findings 1.17).** With the physics map, the
trainer keeps loop glycines left-handed by bending the shared H-bond, sheet and side-chain terms,
and helical glycines pay (lambda helix 3, glpG TM4). FF2's H-bond energy already scores each
residue's H-bonds by that residue's own basin, so a glycine-only set of its three energies is the
smallest context-aware glycine term. `hbond_energy` takes an optional 15-entry `parameters` (12 +
glycine's dE_alpha, dE_beta, dE_other) with a per-residue `residue_class`; a 12-entry config is
unchanged bitwise. The offsets start at zero, train at `hb`'s rate, and the SARW replica zeroes them
with the shared energies.

**Side-chain step 10x smaller (user, 2026-10-02; findings 1.17).** Every retrain from ff2.1 weakened
helices and folds because the side-chain pair table (31,420 coefficients, ~6% with any signal)
random-walks under Adam's normalised step; putting the epoch-0 table back restored the panel. At a
tenth of the rate, per-step noise averages out while a consistent pull still accumulates; every file
stays trained.

**Release by checkpoint selection (user, 2026-10-01; findings 1.17).** The gate stops at convergence
with nothing released (`gate_or_continue.sh`, as in the repo), and `validate_ff.sh <run> ff_3.0
<checkpoint>` releases the chosen one. Training is unchanged; not-converged still trains one more
epoch. The selection panel and the final validation stay disjoint: the panel chooses, Peng and glpG
judge.

**Pre-proline is separate and parked (user, 2026-10-01).** It came from analysis, not from a
simulation failure (Phase 6).

**Retired decisions.**
* Full-map glycine row (289 Fourier modes per map): failed the convergence gate at every
  checkpoint; cancelled 2026-09-28 (Phase 2).
* Per-pair basin offsets on the NDRD maps, the ff_3.0 released 2026-09-30: noise-limited per pair
  (findings 1.12); the released maps favoured alpha_L at every glycine and TM4's helical glycines
  flipped in glpG; a context-free glycine map trained against natives relearns their placement
  (findings 1.14-1.15).
* Ting et al.'s product combining rule: it divides by glycine's pooled map and breaks GLY|GLY
  symmetry (findings 1.13).
* No DSE term on the maps: the SARW replica keeps the rama term, so a DSE term on a map compares the
  unfolded ensemble with the map's own distribution, not with data. Moot now that the map is fixed.

**Bottom-up fallback (PROPOSED 2026-10-01, not approved; findings 1.17).** Used only if Phase 8
validation fails. The user's constraints: bottom-up only for glycine, every other residue and term as
current Upside, training allowed, design change kept to the minimum.
* Measure in-context glycine dG(alpha_L - alpha_R) in all-atom, in the map's own force field and
  temperature (ff99SB-ILDN, 300 K): Charron et al.'s CATH domains for loop, turn and beta glycines
  (excluding the six D-amino-acid domains and the decoy frames), new AWH only for helical glycines
  whose flips are too rare there.
* Compare with Upside's dG for the same glycines at T0, counting only frames whose surroundings are
  native-like. If Upside matches within error in every class, there is nothing to fix.
* Otherwise fit the glycine H-bond offsets to the all-atom dG (the gradient is exact from one Upside
  run), then one ConDiv epoch around them and refit; never fit protein data into glycine's map.
* Not in scope: refitting any non-glycine term to all-atom data, force matching (ff2.1 was never
  fitted to all-atom forces, so glycine would absorb their mismatch), and training the offsets
  against the native-restrained replica (its ~0.99-in-basin target would flatten glycine's
  helix-breaking propensity, findings 1.16).

## Execution Phases

### Phase 1 - validate the trainer on ff2.1 (DONE 2026-09-24)
Port adapted and gated; one epoch from ff2.1 leaves 8 of 9 groups at a fixed point, and the ninth
(dhb) is a mildly unconverged ff2.1 parameter, not a port error (findings 9v).

### Phase 2 - full-map glycine row (CANCELLED 2026-09-28 at step 223)
The glycine group failed the gate at every checkpoint (p = 0) and GLY|X handedness drifted back to
-0.82. Last checkpoint `training/ff30/run_output/epoch_11_minibatch_13`; extract it with the
pre-basin `extract_ff.py` in `training/backup_pre_basin_20260928/`.

### Phase 3 - basin offsets, implementation and tests (DONE 2026-09-28; offset code removed 2026-10-01)
`training/rama_basin.py` (basins partitioning the torus, per-map renormalisation), its verifier and
the ConDiv wiring. Since Phase 8, `rama_basin.py` only records per-residue basin populations.

### Phase 4 - train ff_3.0 with basin offsets (DONE 2026-09-30, midway2)
`training/ff30_basin` from ff2.1; the gate declared convergence at step 114 and released ff_3.0 to
both cluster trees on 09-30 02:25 with the rotamer-BP fix. The rama group's pass was prior-limited
drift, not closure (findings 1.14). Superseded by Phase 8.

### Phase 5 - validation of the basin-offset ff_3.0 (CANCELLED 2026-09-30 20:34, user)
Its glycine maps favour alpha_L everywhere and TM4's helical glycines flipped to phi > 0 in glpG.
The Peng arms and glpG chains were cancelled after ~18 h; partial data kept (remote_jobs.md §1).

### Phase 6 - pre-proline (PARKED 2026-10-01, user)
No action unless evidence appears. Against type-matched controls, residues before a proline show only
small excess losses (+0.024 extended, +0.059 helical, which the right-only rule would worsen), and
central prolines hold their basins better than any residue (findings 1.16).

### Phase 7 - glycine handedness probe (DONE 2026-09-30; findings 1.15)
From equal alpha_R/alpha_L depth the data pull glycine back toward alpha_L, -0.035 [-0.048, -0.022]
per residue read, 35 of 38 maps: a context-free glycine map cannot be trained right-handed.

### Phase 8 - ff3.0 with the physics glycine map and glycine H-bond offsets (APPROVED 2026-10-01, revised 2026-10-02; TRAINING)
Run `training/ff30_glyhb` on midway2, from ff2.1, submitted 2026-10-02 11:52; job state in
remote_jobs.md §1. It replaces `ff30_gly` (the AWH library alone), stopped 2026-10-02 11:51 in
epoch 3; its panels e00-e02 are the comparison.

- [x] **Library.** `parameters/common/rama31.dat`, built by `training/build_gly_library.py` (old
  build in `backup/`; midway2 build byte-identical). Checks run by the builder through
  `upside_config` and confirmed end to end in a `.up` file: the engine's glycine map equals the
  measured surface to 1e-6, the middle glycine of G-G-G is symmetric to 2e-6, every non-glycine
  residue is identical to ff2.1.
- [x] **Trainer.** Offset training removed from `ConDiv.py`, `extract_ff.py` and
  `convergence_gate.py`; `verify_rama_basin.py` deleted. Deployed to midway2 with backups
  (`backup/training_pre_glyfix_20261001`).
- [x] **Experimental check of the symmetric part (findings 1.16).** In the matched system the force
  field is within 0.06 of the GGG experiment in pPII and 0.01 in alpha, so no change is needed. No
  experiment resolves handedness; it rests on two force fields agreeing to 0.045 nats.
- [x] **Release by checkpoint selection.** The midway2 gate stops at convergence with nothing
  released; `validate_ff.sh` releases the chosen checkpoint.
- [x] **Selection panel built.** The 44 L-only CATH domains of Charron et al. (2,823 residues, 198
  glycines), all-atom ff99SB-ILDN at 300 K, against Upside native-start runs of each epoch-end
  checkpoint and of ff2.1 as control, all at T0 = 0.8. Code
  `/project/trsosnic/yinhan/ff3_selection/panel.py` (copy in `scratchpad/ff3_selection`).
- [x] **Engine.** Optional glycine offsets in `HBondEnergy` (value, angle forces, parameter
  derivative). The class accumulation sits after the shared loops: placed beside them, -ffast-math
  reordered the shared sums and changed their last bit (findings 10.12).
- [x] **Config writer.** A 15-entry `hbond.h5` writes the offsets and `residue_class`;
  `patch_glpg.py` carries them into a hybrid seed.
- [x] **Trainer for ff30_glyhb.** Field `hbg`, `zero_for_sarw`, rot learning rate / 10;
  `extract_ff.py`, `check_step.py`, README, up.md.
- [x] **Tests.** 12-entry config bitwise equal to the old build (energy, forces, all derivatives);
  zero offsets bitwise equal to 12 entries; offset derivatives and forces match finite differences;
  one real trainer step (1ga3, local): hbg gets a gradient and moves, rot steps 10x smaller.
- [x] **Deployed** to midway2 (bitwise parity on 1ga3, engine and binary); ff30_glyhb initialised
  from ff2.1 and submitted.
- [ ] **Every finished step** checked with `training/check_step.py`: kinetic-energy ratio, RMSD,
  unfolded-state target, finiteness, parameter drift, glycine readout. Baseline from ff_3.0's epoch
  5: helical glycines alpha_L 0.102 free / 0.005 restrained, H-bond margin 0.080.
- [ ] **Selection panels h00, h01, ...** at each epoch end.
  - Per residue, basin populations (`rama_basin.py` basins), counted on both sides only in frames
    whose residues i-4..i+4 are native-like.
  - Reported by class with a bootstrap over domains: natively helical and beta non-glycine
    residues; helical glycines (alpha_L); natively left-handed glycines (alpha_L). A beta residue is
    scored on the extended region (beta + pPII), since the -100 deg line cuts native strands and
    ff2.1's beta-basin "miss" is a phi shift (findings 1.17).
  - Selection panel and final validation stay disjoint: the panel chooses, Peng and glpG judge.
  - **Rule.** The newest epoch-end checkpoint is the default. An earlier one replaces it only if
    it dominates: no class worse and at least one of helix, helical glycine or left-handed glycine
    better, on |error| by a paired bootstrap over domains (95% interval excluding zero); if several
    dominate, the latest. Why: the training objective is the design, and a 198-glycine panel picking
    the top scorer would mostly pick noise; the panel overrides only on clear evidence. The release
    is held for the user if the choice is worse than released ff2.1 in helix or helical-glycine
    retention, the ff_3.0 failure.
- [ ] **Then lambda, Peng, glpG.** Sync `/beagle3` from `$P` first, and build the glpG seed through
  `patch_glpg.py` so it carries the new term.
  - lambda: wild type toward G46A/G48A. Its helix 2 / helix 4 misorientation has no single
    side-chain cause under ff2.1 (findings 11d), so it is not trained against; re-test it under
    ff3.0: crossing angle, the Q33-F51 dock, helix 1-helix 2 contacts.
- [ ] **Validation.**
  - Free ensembles by native basin: helical glycines' alpha_L clearly below ff_3.0's 0.101, and
    their helical loss not below Ser/Asn's (~0.08), so propensity is not flattened. Natively
    left-handed glycines lose no more alpha_L than the Asn/Asp that sit there (~0.17).
  - glpG TM4 glycine flips against ff_3.0.
  - Peng benchmark paired against ff2.1, de novo arms included.
- Dropped: a host-guest helix benchmark and any correction calibrated to experiment. The Pace &
  Scholtz scale is an 11-system average, and the per-host data are not available to us.
- If validation fails, decide a context term then, with the evidence (bottom-up fallback above). Do
  not add one in advance.

### Phase 9 - repo cleanup for redistribution (DONE 2026-10-01 and 10-02, user)
`py/` and `training/` keep only what another Upside user could run. Personal, cluster-path and
campaign-specific files moved to `scratchpad/redistribution_cleanup_20261001/` (gitignored); the
cluster trees keep their own copies, so running jobs are unaffected (remote_jobs.md §1b). Verified:
every kept file compiles or parses, `env.sh` resolves the repo `.venv`, and all four gate branches
(converged, retrain, max epochs, gate failure) behave in a sandbox. On 10-02 the remaining midway2
specifics went too: `training/env.sh` activates the repo `.venv` only (site modules are loaded
before sourcing it), `train_chain.sbatch` carries no account (account, partition and exclusions go
in `slurm.args`), and README and the ConDiv docstring cite the trainer's GitHub source instead of a
cluster path; originals in `scratchpad/redistribution_cleanup_20261002/training/`.

## Known Errors / Blockers

* **Run only on midway2 broadwl** (user, 10-01). No `amd`, `beagle3` or GPU partitions: CPU jobs
  must not run on the group's GPU allocation. No midway3 jobs; its login node may be used to move or
  read files on `/project` and `/beagle3` (user, 10-02).
* **Beta for the PI's sheet modelling.** No residue type has a significant beta miss at ff2.1, and
  per-pair beta has no signal in these data; per-type beta is ff2.1's sheet mixing energy, trained.
* **lambda's helix 2 / helix 4 packing misorientation** is a known ff2.1 limitation with no single
  side-chain cause (findings 11d); re-tested after ff3.0, not trained against.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
