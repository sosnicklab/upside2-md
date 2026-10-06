# ff3.0: ff2.1's own training workflow, a physics glycine map and glycine H-bond offsets

## Project Goal

Train ff3.0 from ff2.1 with exactly ff2.1's training workflow (Peng et al. 2022 SI), modernised
only, plus the smallest changes that fix glycine:

* the central-glycine Ramachandran row is fixed to the AWH free energy of capped glycine
  dipeptides (`parameters/common/rama31.dat`, built once by `build_gly_library.py` from the
  cluster-only AWH data; the builder is kept in scratchpad, not in the repo); every
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

**Probe a balanced start before any training design (user, 2026-10-04).** Rounds 2 and 3 both
drifted back to ff2.1's alpha_L preference only after days of epochs. Every new trainable term starts
at an estimate of its balanced position (for glycine: alpha_L as deep as alpha_R) and is reported
first with the data's push direction and size from one step (or the gradient split), and whether that
push is placement or energy; a push back toward the known-bad side means the design changes before
any epochs run.

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
with nothing released (in the repo `train_chain.sbatch` runs `ConDiv.py gate` at the target), and
`validate_ff.sh <run> ff_3.0
<checkpoint>` releases the chosen one. Training is unchanged; not-converged still trains one more
epoch. The selection panel and the final validation stay disjoint: the panel chooses, Peng and glpG
judge. **Revised for round 4 (user, 2026-10-05):** a converged gate validates its own final
epoch-end checkpoint at once, on the cluster the run trains on (its `slurm.args`, caslake for both),
as the candidates `ff_3.0_bio` and `ff_3.0_gdepth`; a run that reaches 13 epochs unconverged still
stops for review. The epoch panels still run, and which candidate becomes ff_3.0 stays the user's.

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
the ConDiv wiring. Since Phase 8 only the per-residue basin populations remain, now functions in
`ConDiv.py`.

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

- [x] **Library.** `parameters/common/rama31.dat`, built by `build_gly_library.py` (old
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
- [ ] **Every finished step** checked with the midway2 tree's `check_step.py`: kinetic-energy ratio, RMSD,
  unfolded-state target, finiteness, parameter drift, glycine readout. Baseline from ff_3.0's epoch
  5: helical glycines alpha_L 0.102 free / 0.005 restrained, H-bond margin 0.080.
- [ ] **Selection panels h00, h01, ...** at each epoch end.
  - Per residue, basin populations (the `ConDiv.py` basins), counted on both sides only in frames
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
cluster path; originals in `scratchpad/redistribution_cleanup_20261002/training/`. Then the
directory was cut from twelve files to five: `rama_basin.py`, `extract_ff.py` and
`convergence_gate.py` became functions and the `extract` / `gate` commands of `ConDiv.py`;
`gate_or_continue.sh` folded into `train_chain.sbatch`; `check_converged.py` (duplicated the gate) was
deleted; `build_gly_library.py` (one-time) and `check_step.py` (campaign monitor) moved to
`scratchpad/redistribution_cleanup_20261002/training_pre_merge/`. In `py/`, `martini_gen_params.py`
became the `__main__` of `martini_build_tables.py`. Checked on midway2 data: `extract` writes files
byte-identical to `extract_ff.py`, and `gate` gives the same p-values and verdict as
`convergence_gate.py` on ff30_gly's three epochs; the basin functions are bitwise equal.

### Phase 10 - local test of ff3.0 at its current stage on glpG TM4 and lambda helix 3 (STARTED 2026-10-04, user request)
The cause of the helical-glycine failures is settled (user; memory gly-map-energy-plus-selection):
glycine's deeper alpha_L basin in the PDB map is placement, which Upside applies as energy. Round 1
(symmetrised map) fixed TM4 but is physically wrong; round 2 (trained basin depth) was pushed back to
alpha_L by the training data; round 3 is ff30_glyhb (physics map fixed, glycine H-bond offsets
trained). By epoch 2 its offsets favour glycine's left-handed H-bond by 0.35 E_up, so the question is
whether training is re-learning placement through the offsets, and what design gives Upside only the
energy part for glycine while fixing both glpG TM4 (GLY136, 143, 149) and lambda helix 3 (G46, G48).
ff30_glyhb keeps training on midway2 meanwhile (user). All runs on the Mac Studio (M1 Ultra, 16
performance cores); midway2 is only read from.

Force fields, one protocol, so differences are the force field's:
* `ff21`: released ff2.1 (PDB glycine map), the known-bad control.
* `ff21_awh`: ff2.1 with the physics glycine map, untrained: round 3's starting point.
* `e02`: ff30_glyhb `epoch_02_minibatch_18`, the newest epoch end (also panel h02).
* `e02_gly0`: e02 with the three glycine offsets zeroed, everything else trained.
ff21_awh against e02 says whether two epochs of training broke what the physics map fixed; e02_gly0
says whether the offsets carry it or the other trained changes (the shared H-bond drift, sheet,
side chains) do.

- [x] **A. Setup (done 10-04 20:10).** Rebuilt locally (old binaries kept as `obj/*.bak_sep10_20261004`); initial potentials match the cluster builds: lambda ff2.1 -196.48 and e02 -193.92 (midway2 tree), glpG 79HIS ff2.1 -24914.84 (`/beagle3`); patch gate on the pristine seed passed (2.6e-6). Inputs in `scratchpad/ff3_local_test/`. Back up `obj/` binaries, rebuild with `make` in `obj/` (the Sep 10 binary
  predates the offsets and silently ignores them); check `residue_class` is in the binary, and that
  1ga3 energies under ff2.1 and under e02 match the midway2 build. From midway2: extract e02 with
  `$P/training/extract_ff.py`; copy it, ff30_glyhb's 15-entry init `hbond.h5`,
  `$P/training/patch_glpg.py`, the pristine seed
  `popepopg_REMD/seeds/glpG-RKRK-79HIS.up.bak_production_handoff`, `checks/gly_tm4_flip.py`,
  `checks/glpg_tm_windows.py` and `ff3_benchmark/bench_run.py` to `scratchpad/ff3_local_test/`.
- [ ] **B. glpG TM4 (the focus; 12 runs started 10-04 20:09, `scripts/run_glpg_md.sh`, analysis `scripts/tm4_local.py`).** Patch the seed per force field with `patch_glpg.py` (its gate
  must reproduce ff2.1 within 1e-4); it changes only the protein's rama, H-bond and rotamer tables,
  so the membrane and the SC-env and BB-env interactions stay as built, dt 0.009 and inner_steps 4
  as the seed sets. Single-temperature MD of 79HIS at T 0.80 (where ff_3.0's flips were strongest),
  3 seeds for each of the four force fields, ~4000 time units each, 12 runs at once (~6 h); then T
  0.70 if needed.
  Read: phi of GLY136, 143, 149 (phi > 0 is a flip), TM4 (134-151) and TM1 (30-48) helix fraction by
  the cluster scripts' rules, as time series. Single-temperature runs are harsher than REMD rungs, so
  force fields are compared with each other, not with REMD numbers.
- [ ] **C. lambda helix 3.** Wild type from native, the Table S2 14-replica ladder, dt 0.009,
  `bench_run.py`'s recipe, ~400k time units for e02 (and ff21_awh or e02_gly0 as B suggests),
  against the ff2.1 cluster arm at matched time (helix 3 alpha_R 0.97 -> 0.60 by ~300k, 11d). After
  B, for cores (~8-10 h per arm).
- [ ] **D. Design.** Within the user's rule (Upside gets only the energy part for glycine; no
  glycine-specific basin term fitted to natives, which relearns placement): attribute any failure;
  then test whether glycine offsets fitted bottom-up to all-atom in-context dG (the plan's fallback)
  can keep helical glycines helical without breaking left-handed ones, by reweighting the B and C
  frames over the offsets (per-glycine H-bond counts by basin). Propose a design with that
  evidence; no training change without user approval.
  - Data for the bottom-up fit: Charron et al.'s `training_a_cg_model.zip` (Zenodo 15465782, 65.4 GB,
    md5 619f1b18b04db5502e4cab5b75f977a5), downloading to `~/Downloads` since 10-04 19:41 (user).
    Its frames are mapped to N, CA, CB, C, O, which is all Upside's backbone H-bond term needs (it
    infers H and O from N, CA, C), so Upside's per-glycine H-bond counts by basin can be evaluated
    on the all-atom ensemble itself. Use only the real frames (not the 0.5 A decoys) of the 44 all-L
    domains (not the six D-residue ones, findings 1.17).

### Phase 11 - round 4: a selection-free glycine map, fitted bottom-up and frozen (dt 0.015 runs STOPPED 2026-10-06 09:40, user; RESTARTED at dt 0.009 on midway2 as ff30_bio_dt009 and ff30_gdepth_dt009)
The defect is the map (user, 2026-10-04; findings 10.13): Upside uses glycine's PDB coil statistic as
energy, and its alpha_L excess is fold selection. A design through H-bond terms was rejected: a map
defect is fixed on the map. Rounds 2 and 3 showed that any glycine parameter trained against natives
takes the alpha_L preference back, so the map is fitted to selection-free data and stays frozen.

Design, as built (findings 1.19; up.md 2.8):
1. Target: glycine's in-chain distribution in BioEmu's plain-MD octapeptides (Zenodo 15641199;
   amber ff99sb-ildn, 300 K; sequences simulated alone, so no fold selects a conformation),
   glycines at residues 2-5. Charron's adaptive frames of the same peptides are biased (glycine by
   0.16-0.28 toward alpha_L) and are not the target.
2. Every central-glycine entry fitted by iterative Boltzmann inversion in Upside on the same
   octapeptides with ff2.1's other terms, so Upside's own sterics, side chains and H-bonds in local
   context count once: pooled GLY|X on X-G-Y, GLY|left|GLY on G-G-Y, GLY|right|GLY on X-G-G,
   GLY|right|PRO on X-G-P (user: use the data where they exist). **Units: the library convention**
   (up.md 2.8): a map reproduces its source at T_up = 1, so Upside runs at T_up = 1 against the
   300 K target; a fit at T_up 0.8557 scales glycine 1.169x too small and was discarded.
3. The transition region is not designed in Upside for any residue and contrastive divergence
   cannot train it, so outside the data (< 10 all-atom frames per cell) each entry keeps NDRD's own
   coil top with its barrier height above alpha_R (user). Glycine's coil and sheet entries are the
   same map.
4. Frozen in ConDiv (the library is a fixed input), ff2.1's 12-entry hbond.h5, no glycine-specific
   trainable parameter anywhere; side-chain learning rate 10x smaller as in round 3.

Done:
- [x] Library: map 6 of the fit, now `parameters/common/rama31.dat` (the AWH version is in
  `backup/rama31.dat.bak_pre_bioemu_20261005` and git history). Its own pass matches BioEmu's X-G-Y
  glycine within 0.01 per basin (ln aR/aL -0.55 vs -0.535). Scripts, data, logs:
  `/project/trsosnic/yinhan/checks/gly_bioemu_map/` (README).
- [x] Push probe P1 (local, ff30_glyhb's exact worker, 72 proteins, no update): d -0.014 [-0.041,
  +0.012], against -0.105 [-0.150, -0.057] from rama31's start; the training data come to rest at
  about BioEmu's L-R difference. By the pre-agreed rule (CI including zero = do not train) this
  was "do not train".
- [x] glpG TM4 on the untrained start (ff2.1 terms + library): no glycine flips, TM4 0.95.
- [x] Training started on the user's decision (2026-10-05): TM4 is stable on this start and no push
  is left for other terms to absorb. ff30_bio initialized and submitted 07:44 (job 49179449, 76
  steps, gate up to 13 epochs); the empty glycine-offset field was tested through one parameter
  update and through extract_ff.py first. Panel of the untrained start (`bio_start`, job 49179451).
- [x] Non-glycine maps against BioEmu, measured in Upside on all 1,100 octapeptides plus the panel
  by native class (2026-10-05, findings 1.19): kept NDRD (recommendation to the user). Their gaps
  to ff99sb-ildn are as large as glycine's was, but the reference is weaker for L residues and the
  fit's direction (more alpha_R, more alpha_L) would worsen the folded-protein errors that exist.

**Run 2, ff30_gdepth (user, 2026-10-05): the claim is a force field trained on the original training
set alone, with no outside data, that keeps glpG TM4 stable; BioEmu is only the external check.**
* Library: ff2.1's own (NDRD), with glycine's alpha_R and alpha_L basin depths trainable: one pooled
  offset pair on all 37 GLY|X maps (user's choice, the quantity the push probe measured); GLY|GLY,
  GLY|right|PRO, every other map and the sheet group stay NDRD (user).
* Start: the offsets at which the GLY|X maps hold BioEmu's alpha_R / alpha_L weight on average
  (c_aR -0.0918, c_aL +0.4621; `gdepth_start_offsets.py`), the only place BioEmu enters.
* Update: round 2's rule (damped Newton step per epoch matching free to native-restrained
  populations over the GLY|X reads, Gaussian prior sigma 1 nat on NDRD, no DSE, 10% of proteins held
  out); log `run_output/rama_rounds.txt`. Everything else as ff30_bio (ff2.1 init, 12-entry hbond,
  rot lr 10x smaller). The convergence gate judges only the Adam groups; the depth's convergence is
  read from the rounds log.
* Trainer: round 3's ConDiv.py with the offset hooks of round 2 ported (`$P/training/ff30_gdepth/
  trainer/`, copied into run_output at initialisation; also in `gly_bioemu_map/gdepth_trainer/`).
  Tested before submission: round 0's library (untrained entries and sheet identical to NDRD; GLY|X
  at the target weights), one end-of-round update from 24 real proteins, extraction.
* Reading: per epoch the depth (dL - dR, start +0.554) and TM4 on the extracted checkpoint; at the
  end, where dL - dR settles against BioEmu's start (external check), and the panel against
  `gdepth_start` and ff30_bio.
* Shared scripts changed for it (backups `.bak_pre_gdepth_20261005`): `extract_ff.py` writes a
  depth-training checkpoint's trained library and finds the run's own trainer for the initial
  checkpoint too; `check_step.py` skips the empty glycine-offset field and prints the depth.

**Revised 2026-10-06 (user): both runs stopped and restarted at dt 0.009 on midway2.** At dt 0.015
the first runs lost what the library had gained:
- Both shared H-bond margins fell below +0.10 (ff30_bio +0.050 at step 37).
- Both epoch-0 panels lost folding (0.60 to 0.47) and helix (findings 1.22).
- glpG TM4 fell with every checkpoint (bio_start 0.95, b00 0.92, b01m12 0.87).
- Free replicas were destroyed in 6 protein-steps.

The time step is the protocol's one integration setting that departs from Upside's standard:
- 0.015 came with the FF2 trainer on 09-24 as the port's schedule.
- Master, the FF1 trainer and the Peng benchmark use 0.009.
- FF1-form runs at 0.009 had 0 destroyed replicas in 13,714 protein-steps, against 20 in ~13,800
  for FF2 runs at 0.015 (findings 1.21).

So the restart puts the model's own integration step back; it does not tune physics to avoid a
crash. Everything else is unchanged:
- the ff2.1 start, the libraries and START offsets;
- the side-chain lr / 10 and the 76-step target;
- the gate and the automatic validation (candidates `ff_3.0_bio`, `ff_3.0_gdepth`; on broadwl, glpG
  from `popepopg_REMD_mdw2`).

A step takes ~1.67x longer (8000 time units at 0.009), so `train_chain.sbatch` links run 54 steps
instead of 90.
- [x] Trainers at dt 0.009: `$P/training/ConDiv.py` (ff30_bio's) and a copy of ff30_gdepth's
  `trainer/`; repo `training/ConDiv.py` the same. Backups `.bak_pre_dt009_20261006`.
- [x] New run dirs `ff30_bio_dt009` and `ff30_gdepth_dt009`, copied from the stopped runs:
  - `upside_input` as hardlinks;
  - midway2 `slurm.args`;
  - `after_training.sbatch` with the candidate names.

  Then initialize, and check that the initial force fields are byte-identical to the stopped
  runs' and that the trainers differ only in dt.
- [x] `train_chain.sbatch` STEPS_PER_LINK 54 (cluster and repo); `submit_new.sh` tags `b9_EE`
  and `d9_EE`.
- [x] Submitted 10-06 09:46 on broadwl, 49194446 / 49194447; both started at once from
  `initial_checkpoint.pkl`, 54 steps per link. Watch cron `250534b9`. The trainers differ from the
  stopped runs' only in dt, and the initial force fields are byte-identical
  (`checks/dt009_init_20261006`).

Next (remote_jobs.md has the watch):
- [ ] Every step: check_step.py (all finite, KE/1.5kT, restrained RMSD, unfolded target); watch the
  shared H-bond margin E_other - E_alpha (start +0.192) and the helical / left-handed glycine
  populations. **Hold for the user if the shared margin falls below +0.10** (round 3's failure mode:
  the alpha_L pull moving into shared terms).
- [ ] Each epoch end: selection panel `b9_EE` / `d9_EE` (submit_new.sh), and locally the glpG TM4
  test on the extracted checkpoint (**12 seeds**, T 0.80, 4000 tu, `run_glpg.sh` after
  `patch_glpg.py`, on the watch's computer). Revised 10-06 (user) from 3 seeds: TM4 is lost seed by
  seed, and 3 seeds resolve only a start at one in three going to three in three, while 12 resolve
  one in four against three in four (findings 1.22). Each checkpoint's unwound and flipped seed
  counts are compared with its own start at 12 seeds (`bio_start` = `ff21_bioT1_6`, `gdepth_start`)
  by Fisher's exact test, with ff21_released and the same-step dt 0.015 checkpoint (b00, d00, ...)
  as context; direction is read across the epoch ends. Nothing is released without the user.
  - [x] gdepth_start seeds 1-3 (no local test of it existed): one seed of three unwinds, as at d00.
  - [ ] Seeds 4-12 of both starts, before the first epoch-end set.
- [x] Automatic validation at convergence (user, 2026-10-05; caslake; Peng's 32 arms, with lambda
  helix 3, and glpG's 4 REMD variants per candidate). Deployed 10-05 20:05 (remote_jobs.md §1
  "Automatic validation"):
  - [x] `gate_or_continue.sh`: a converged gate runs `validate_ff.sh` on its epoch-end checkpoint.
  - [x] `validate_ff.sh`: submissions from the run's `slurm.args`; Peng at the script's 12 h with
    `runs/<tag>/input` made first; glpG per candidate (`<V>.<FF>` seeds and directories, nothing
    shared overwritten or deleted).
  - [x] `popepopg_REMD_mdw3/` for glpG on caslake (`env_shared.sh` + `$P`'s binary; mdw2's
    `run_remd.py` unchanged); `bench.sbatch` resubmits with the job's own node exclusions.
  - [x] `after_training.sbatch` of each run names its candidate.
  - [x] Tests: sandbox gate branches; `validate_ff.sh` end to end with `--test-only`; base-seed
    independence; engine parity of the two trees on a Peng config (bitwise).
  - Short real test jobs dropped (user: the pipeline already ran in production; findings 10.16); the
    first production Peng arm and glpG block on caslake are read by the watch.
  - Not changed, flagged: `run_remd.py` rolls back NaN replicas and continues, which the NO GUARDS
    rule forbids; it is copied unchanged so glpG stays comparable with the ff_2.1 and ff_3.0
    campaigns. Whether to remove it is the user's decision.
- Known limitations: X-G-P's helical basins keep NDRD values (too few BioEmu frames); a poly-Gly run
  gets the average of the G-G-Y and X-G-G entries, which carry their L-neighbour contexts (Upside's
  maps see only nearest neighbours); ff99sb-ildn may overstate in-chain alpha_L (findings 1.19);
  lambda's helix 2/helix 4 packing error stays a known limitation (11d).

### Phase 12 - poly-Gly: all-atom Ac-(Gly)20-NHMe as Upside's reference (STARTED 2026-10-05, PI's suggestion; collapse 49186415 queued)
Where does a chain with no chiral residue go, and does Upside send it to the same place? With
neutral caps the all-atom chain is achiral, so at equilibrium every glycine has alpha_R = alpha_L,
beta = beta' and pPII = pPII' exactly; what the run measures in those pairs is its sampling error
(the blank, findings 6.7). Upside's GLY|left|GLY and GLY|right|GLY entries carry L-neighbour
contexts (NDRD; in the BioEmu library, fitted on G-G-Y and X-G-G, findings 1.19), so Upside's
poly-Gly can be chiral: this is the one system where the glycine map's own handedness shows with
nothing else present, alongside what the chain does globally (collapse, pPII, helices, hairpins).

Decisions (user, 2026-10-05): Ac-(Gly)20-NHMe; amber99sb-ildn / TIP3P at 300 K (BioEmu's force
field, so round 4's target, and our AWH dipeptides'); midway2 broadwl; compared with Upside under
ff2.1 and bio_start, later each round-4 epoch.

All-atom, the `gly_peptides` protocol and its mdp files (PME 1.0 nm, h-bond constraints, 2 fs,
v-rescale 0.1 ps, C-rescale 1 bar; `gmx2024` with `gcc/10.1.0`):
1. Build with `build_peptide.py` at phi = psi = 180: the fully extended chain is its own mirror
   image, so the start carries no handedness. pdb2gmx `-ignh`, dodecahedron `-d 1.0` around the
   ~8 nm chain (~70k atoms).
2. Collapse: min, 500 ps NPT, 4 replicas with distinct seeds for one 36 h link (4 x 7 threads,
   ~20 ns each, estimated). Rg(t) shows whether the chain has relaxed; this phase is discarded.
3. Production box: each replica's last frame in one common dodecahedron whose image distance is the
   largest chain diameter of step 2's second half plus 2.4 nm; solvate, min, 500 ps NPT. The
   protein's minimum periodic-image distance (`gmx mindist -pi`) must stay above the 1.0 nm
   cutoff through production; if it does not, the box is too small and production is redone.
4. Production: 4 replicas, 1 us target, protein frames every 10 ps, chained 36 h links (successor
   queued at link start, as `train_chain.sbatch`); 4 x 7 threads on one node; the first link
   measures ns/day. Estimate ~55 ns/day per replica, ~18 days, ~12k core-hours.
5. Converged when the blank is zero within its standard error over replicas, holds between the
   halves of each replica, and the replicas agree on the Rg distribution; extended otherwise.

Upside, `ibi_run.py`'s recipe (the force field's rama, sheet, hbond, side chains, environment and
bb_env; common rama_reference; dt 0.009; T_up = 1, the library's convention against a 300 K
ensemble, up.md 2.8; a frame every 10 tu):
6. G20 from the 20 glycines' N, CA, C of the same extended build. ff21_released; bio_start
   (`extract_ff.py` on ff30_bio's initial checkpoint); each round-4 epoch as the watch extracts it.
   8 seeds each, length set from a measured short run; on the Mac, cost stated first (findings 10.7).

Analysis, one script for both (all-atom xtc via mdtraj, Upside output via tables; interior glycines
2-19): per-residue basin populations with `rama_basin.py`'s mirror-exact basins (copied, so the
definition lives in one place); ln(aR/aL) per residue and pooled; Rg over N, CA, C; CA1-CA20
distance; helix runs by hand (>= 4 consecutive residues in alpha_R or in alpha_L); CA contact map.
Figure: all-atom gray, ff2.1 red, bio_start blue.

Files: `/project/trsosnic/yinhan/polygly/` on the cluster holds build, mdp, sbatch, scripts and a
README (authoritative); the Mac's `scratchpad/polygly/` holds the Upside runs and the analysis.
Jobs are recorded in remote_jobs.md.

- [x] Build and collapse job (`collapse.sbatch`, 64,976 atoms, grompp clean, physics diffed
  identical to `gly_peptides/prod.mdp`), submitted 10-05 15:20 as 49186415
- [ ] Production box from the measured diameter; production chain submitted
- [x] Upside: an all-glycine config builds (2.8e6 tu/h on one core); ff21_released and bio_start,
  8 seeds x 200,000 tu each: both chiral, ln(aR/aL) -0.72 and -0.45 (findings 1.20)
- [ ] Analysis script; Upside results first, all-atom at each ~200 ns
- [ ] Round-4 epochs as they come

## Known Errors / Blockers

* **The midway2 tree's `training/` keeps the pre-merge layout until ff30_glyhb is released.**
  ff30_glyhb's gate calls `$P/training/gate_or_continue.sh` and its monitor is
  `$P/training/check_step.py`, both gone from the repo since the 10-02 merge (Phase 9). Do not
  sync the repo's `training/` to the cluster before then (remote_jobs.md, "Resume here").
* **CPU work runs on midway2 broadwl unless the user names midway3 caslake** (user, 10-01; on
  10-05 the user moved both round-4 trainings and their start panels to caslake for earlier
  starts, and chose caslake for their automatic validation). No `amd`, `beagle3` or GPU partitions: CPU jobs must not run on the group's GPU
  allocation. midway3's login node may be used to move or read files on `/project` and `/beagle3`
  (user, 10-02).
* **Beta for the PI's sheet modelling.** No residue type has a significant beta miss at ff2.1, and
  per-pair beta has no signal in these data; per-type beta is ff2.1's sheet mixing energy, trained.
* **lambda's helix 2 / helix 4 packing misorientation** is a known ff2.1 limitation with no single
  side-chain cause (findings 11d); re-tested after ff3.0, not trained against.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
