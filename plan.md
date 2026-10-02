# ff3.0: ff2.1's own training workflow, with the Ramachandran maps trained as per-pair basin offsets

## Project Goal

Train ff3.0 from ff2.1 with **exactly ff2.1's training workflow**, modernised only, plus one project
addition: the Ramachandran maps where the literature and ff2.1's own error both point to a local
defect, the glycine-centred and pre-proline maps, get a few trained basin offsets each; every other
map stays at NDRD. The result is released as `parameters/ff_3.0`, which is absent from the tree
until then.

Why: the NDRD maps are statistics over loop sites of folded proteins, so they carry evolutionary
placement (glycine is put where a fold needs alpha_L) as well as local physics, and Upside adds them
as a pure energy to every residue. Glycine's alpha_L excess is mostly a property of that loop-only
site selection: ln(aR/aL) is -1.19 in NDRD and -0.58 over all residues of the 456 training natives
(findings 1.9).

## Architecture and Key Decisions

**Trainer (unchanged).** `training/ConDiv.py`, O. Kleinmann's Python 3 port of Peng's FF2 dual-target
trainer, restored to the Peng et al. 2022 SI: per protein and step one native-restrained replica, 12
free replicas at T = 0.8 to 1.1, one SARW replica; 8000 time units, second half analysed; contrast =
NSE + 0.3 * DSE; 456 proteins in 19 minibatches of 24. Every difference from the port is justified in
the file's docstring.

**ff2.1's parameter set (unchanged, user's choice 2026-09-24).** rot (pair, coverage, hydrophobe), the
20 x 3 sigmoid burial parameters and 400 weights, the backbone term's scale, the three H-bond
energies and the second-H-bond term, the 20 sheet values. Not trained, as in ff2.1, because the
engine returns no derivative: the backbone term's center, sharpness and hbond weight, and `hbond.h5`
entries 4-11. About 42,000 values, mostly spline coefficients of smooth 1D functions held by priors.

**Revised 2026-09-27: the Ramachandran maps are trained as per-pair basin offsets, replacing the
full-map glycine row.** The glycine row trained 289 Fourier modes per map on a median of 169 native
sites per map (12,138 values on 3,204 glycines) and has failed the convergence gate at every
checkpoint through step 209 while every other group passes. Training all 800 maps that way would be
231,200 values on a median of 107 sites per map, fewer than one site per value. So:

* **Parameter.** For each directional map k = (central X, direction, neighbour Y):
  `E_k(phi,psi) = E_k,NDRD(phi,psi) + sum_b c_k,b * w_b(phi,psi)`. The NDRD map is the fixed base
  (steric boundaries and within-basin shape); only the basin balance is trained.
* **Basins, and each offset is a weight on one basin (user, 2026-09-28).** The map is renormalised
  after the offsets are added, as the NDRD maps are, so the shape inside every basin is the base's
  and only its depth, its frequency, is trained. The basins partition the torus, so no probability
  can move into an untrained region: alpha_R (phi < 0, -100 < psi < 50), beta (phi < -100, psi
  outside that band), pPII (-100 < phi < 0, psi outside it), alpha_L (mirror of alpha_R) and
  `other` (phi > 0 outside alpha_L). For a central glycine, which populates `other`, it is split
  into its mirror halves beta' and pPII' so that GLY|GLY can pair each basin with its mirror.
  Edges are logistic with a 3 deg scale (13 deg from 10% to 90%, the sharpest the 5 deg grid
  represents): measured on the NDRD maps, an offset is realised at a median 98% of its value where
  a basin's probability lies, and with offsets of +-1 the energy change inside a basin spreads by
  a median 0.03 nats around its single depth shift. The 10 deg edges of `secstr_bias` would tilt
  the shape instead (44 deg transitions, no flat interior in beta).
* **Which maps are trained. Revised 2026-09-28 evening (user; findings 1.12, 1.13).** Training all
  840 maps (3,394 identifiable offsets) was noise-limited: split-half reliability 0.09, the true
  per-pair corrections (SD ~0.04 nats) below one epoch's noise (0.07-0.08), ~70% of that noise from
  the protein set itself, so the fixed 456 proteins, kept for the comparison with ff2.1, cannot
  resolve it. The literature ranks pre-proline first (the CD(i+1) clash removes alpha_R) and glycine
  next (no CB; the PDB map carries placement), and finds other neighbour effects small; ff2.1's
  free simulations miss most on exactly these (pre-PRO alpha_R +3.2 points, z 9.7; GLY alpha_L +2.1,
  z 4.8). Trained, 158 offsets on 60 maps: GLY|X (38 maps) alpha_R, alpha_L, beta; GLY|GLY (2)
  helix and beta, each tied to its mirror; X|right|PRO (20, every central type but GLY) alpha_R and
  beta. In each trained map the untrained basins together are the reference, whose weight follows
  from the normalisation, so the basins still partition the torus. Every other map is written
  unchanged. The mixture only half-transmits the pre-proline map, so its offsets fix that miss in
  part; the product rule of Ting et al. would fix it, but divides by glycine's alpha_L-rich pooled
  map and so breaks GLY|GLY symmetry (findings 1.13), and is not adopted.
* **Each pair is its own parameter set (user's rule, findings 1.8).** Offsets are indexed by
  (central, direction, neighbour) and are never tied, pooled or shared: the offsets of ALA|GLY are
  never used for GLY|GLY, or for any other map. The prior of each map is centred on that map's own
  NDRD values. An offset of map k is applied to the coil entry of pair k only, so it acts only on
  residues whose left or right pair is k; never to the sheet group, since a central cis-proline
  reads PRO's sheet entry and an offset there would act on both.
* **Target and update (user's design).** After each full round over the training proteins, for every
  map and basin: basin population of the native-restrained replica against the free replicas (the
  three coldest, reweighted to T0 and mixed 0.6/0.3/0.1, as ConDiv's NSE), accumulated over all
  residues that use the map. **Revised 2026-09-28:** each offset takes a damped (0.5) Newton step on
  the MAP objective, the native basin counts under the model's populations with a Gaussian prior of
  width 1 nat on the offset (the penalty first proposed, which the first implementation had replaced
  with Dirichlet pseudo-counts). With many residues in a basin this is `0.5 T0 ln(p_free/p_native)`;
  with almost none the prior bounds the step and an unsupported offset decays to zero. The
  pseudo-count version gave steps of 1.76 nats in round 1 in basins with no residues in either
  ensemble (a ratio of two near-zero numbers); on the same data the MAP step's largest is 0.44 and
  none exceeds 0.5. The restrained replica is the target because it has the same thermal breadth
  as the free replicas; a static crystal histogram would sharpen the wells to cancel it.
* **Data.** Over all residues of the 456 natives: median 107 sites per map, 24 maps under 20, 12
  under 10; cis-proline maps about 5 each (115 cis-prolines). About 21 sites per offset for the
  median map, which resolves basin populations to ~0.2-0.3 nats.
* **No new engine or config code.** Trained offsets are written into a library file of the same
  format as `rama.dat`, read by the unchanged `upside_config.py`.
* **Held-out check.** 10% of the 456 proteins are kept out of the map update and simulated each
  round. Its mismatch must be read against its own sample-size floor (Known Errors).

**GLY|GLY stays mirror-symmetric (user, 2026-09-27).** Its NDRD base is symmetrized, and within that
one map each offset is held equal to its mirror basin's (alpha_R with alpha_L, beta with beta',
pPII with pPII'). The constraint lives inside one map, so it does not break per-pair independence.
**Revised 2026-09-28 (findings 1.11):** the sheet GLY|GLY entry is symmetrized too, because
`upside_config` mixes every coil map with its sheet map and NDRD's GLY|GLY sheet map holds 94% of
its weight at phi < 0 (a glycine between glycines got ln(beta/beta') = +0.14 to +0.32). And
symmetrizing now averages probabilities, pooling every site with its mirror image, instead of
energies: the energy mean is a geometric mean of probabilities, which would empty the sheet map's
pPII on both sides, and gave the coil maps' helical basins 0.205 and 0.229 of the probability
each, against 0.213 and 0.246 for the probability mean.

**No DSE term on the offsets (user, 2026-09-27).** In ConDiv the SARW
replica keeps the rama term and runs at the hottest temperature, so the DSE term on a map compares
the unfolded ensemble with the map's own distribution, not with data. The offsets are then shaped by
the native-state term and the prior only; the DSE term still trains every nonlocal group as in ff2.1.

## Execution Phases

### Phase 1 - validate the trainer on ff2.1 (DONE 2026-09-24)
Port adapted and gated; one epoch from ff2.1 leaves 8 of 9 groups at a fixed point, and the ninth
(dhb) is a mildly unconverged ff2.1 parameter, not a port error (findings 9v).

### Phase 2 - full-map glycine row (CANCELLED 2026-09-28 at step 223)
Cancelled at the user's discretion: its glycine group failed the gate at every checkpoint (p = 0),
its GLY|X handedness had drifted back to -0.82 (the library is -0.97), and it put the DSE term on
the maps. Last checkpoint `training/ff30/run_output/epoch_11_minibatch_13`, kept for comparison;
extract it with the pre-basin `extract_ff.py` in `training/backup_pre_basin_20260928/`.

### Phase 3 - basin offsets, implementation and tests (DONE 2026-09-28)
- [x] `training/rama_basin.py`: basins partitioning the torus (continuous across phi = +-180,
      exactly mirror-symmetric on the grid), 840 maps and 4,234 offsets keyed by (central,
      direction, neighbour), library writer (coil entries only, every map renormalised), per-residue
      populations, Dirichlet-prior log-ratio update
- [x] `training/verify_rama_basin.py` PASSES locally and on midway3: one offset moves only its own
      map's entry, reaches exactly the residues that read it (most-read map, both termini,
      GLY|GLY, GLY|X) and raises them inside its basin; GLY|GLY exactly symmetric
- [x] `ConDiv.py`: full-map glycine path removed; per-step accumulation, `rama_step.npz` for the
      gate, round update at each epoch end, 46 proteins held out; `extract_ff.py`,
      `convergence_gate.py` (offsets as a sign-flip-tested group) and a synthetic round checked
      locally: glycine offsets move the right way, untouched maps stay exactly 0, the gate flags a
      systematic pull
- [x] Cluster-portable: `env.sh` picks midway2's tree venv where its interpreter exists, else the
      shared /beagle3 venv (torch 2.6.0+cpu added, the same versions); both verified to import the
      tree's engine; partition and node exclusions per run in `slurm.args`; a converged gate stops
      for review instead of releasing

### Phase 4 - train ff3.0 from ff2.1 with basin offsets (DONE 2026-09-30, midway2)
Trained from ff2.1 with all offsets zero (`training/ff30_basin`), three rewinds to step 19 on
2026-09-28 (MAP step, GLY|GLY sheet symmetry, the 158-offset set; remote_jobs.md). The gate declared
convergence at step 114 (end of epoch 5, six offset rounds) and released ff_3.0 to both cluster
trees on 2026-09-30 02:25, after installing the rotamer-BP fix. Training mismatch over the trained
maps 0.040 -> 0.034; held-out 0.081-0.086 -> 0.083, noise about its floor (Known Errors). The rama
group's pass is prior-limited drift, not closure, for the X|right|PRO offsets (findings 1.14).
- [ ] Copy `parameters/ff_3.0` into the local repo

### Phase 5 - validation of ff_3.0 (CANCELLED 2026-09-30 20:34 by the user)
ff_3.0 failed: its glycine maps favour alpha_L everywhere, TM4's helical glycines flipped to phi > 0
in glpG, and training a context-free glycine map cannot fix it (Phase 7). The 32 Peng arms and 4 glpG
chains were cancelled after ~18 h; partial data kept (remote_jobs.md §1). A replacement is proposed
in findings 1.16, decision pending.

### Phase 6 - pre-proline (PARKED 2026-10-01, user)
Not an established problem. It came from AI analysis (findings 1.13-1.14), not from a simulation
failure: the class gap ConDiv saw was composition, and the mixture rule's extra pre-proline alpha_R
has not been shown to change any observable. Kept separate from the glycine work; no action unless
evidence appears. Retraining from ff2.1 with no offsets returns pre-proline residues to ff2.1.
Checked 2026-10-01 against type-matched controls (findings 1.16):
- Central prolines hold their basins better than any residue.
- Residues before a proline show only small excess losses: +0.024 extended, +0.059 helical.
- The helical excess would get worse under the right-only rule.
- No physics map for proline.

### Phase 7 - glycine handedness probe (DONE 2026-09-30; findings 1.15)
One epoch from the ff_3.0 checkpoint with every GLY|X map's alpha_R and alpha_L at equal depth.
Stopped by the user at 15 of 19 steps once the direction was settled: the data pull glycine back
toward alpha_L, -0.035 [-0.048, -0.022] per residue read, 35 of 38 maps. Training a context-free
glycine map relearns the natives' placement and cannot make glycine right-handed.

### Phase 8 - glycine map from physics (APPROVED 2026-10-01; TRAINING, local then midway2)
The problem, in the user's words: Upside applies the PDB Ramachandran map as pure energy, when it
is part local energy and part selection bias. Glycine is where selection reverses the handedness.

Data from folded proteins cannot separate the two (findings 1.15-1.16). So Upside gets only the
energy part, from a source where no fold selects anything: the AWH capped-peptide surfaces.
Glycine's placement in folds must then come from the non-local terms.

Glycine only. Pre-proline is separate and parked (Phase 6).

- [x] **Library.** `parameters/common/rama31.dat`, built by `training/build_gly_library.py`
  (old build in `backup/`; midway2 build byte-identical). Glycine row from the AWH surfaces:
  - each GLY|X entry: the pooled surface;
  - GLY|GLY: the symmetric part only;
  - GLY|right|PRO: its own measured surface, since it differs ~20-fold from the pool;
  - store each entry with `rama_map_pot_ref` subtracted, so the engine applies the measured
    surface exactly;
  - put the same entry in the coil and sheet groups;
  - every other entry stays bitwise as in `rama.dat`.
- [x] **Map checks.** All pass, run by the builder through `upside_config` and confirmed end to
  end in a `.up` file:
  - the engine's glycine map equals the measured surface to 1e-6;
  - the middle glycine of G-G-G is symmetric to 2e-6;
  - every non-glycine residue is identical to ff2.1.
- [x] **Trainer.** The offset training was removed from `ConDiv.py`, `extract_ff.py` and
  `convergence_gate.py`, and `verify_rama_basin.py` was deleted. `rama_basin.py` now only records
  per-residue basin populations. The library is a fixed input. Deployed to midway2 with backups
  (`backup/training_pre_glyfix_20261001`).
- [x] **Experimental check of the symmetric part (findings 1.16).**
  - Result: in the matched system the force field is within 0.06 of the GGG experiment in pPII and
    0.01 in alpha, so no change is needed.
  - No experiment resolves handedness. It rests on two force fields agreeing to 0.045 nats.
- [ ] **Training.** Retrain from ff2.1 with ff2.1's workflow: glycine map fixed, no map offsets,
  no new parameters.
  - Run dir `training/ff30_gly` on midway2: initialised (starts at ff2.1's values), with its gate
    job written.
  - Not submitted: midway2 refused the chain (`AssocMaxCpuPerJobLimit`). A start on midway3 `amd` was
    cancelled by the user (no CPU jobs on the GPU allocation).
  - **Running locally since 2026-10-01 09:16**: `training/ff30_gly_local` on the Mac Studio, two
    workers at a time (`CONDIV_LOCAL_WORKERS=2`), 94 min per step.
  - **Running on midway2 since 12:16:47** (chain 49135913), resumed at step 1 from the local step 0.
    The Mac Studio run must be stopped by `sync_to_midway2.sh` there (remote_jobs.md §1).
  - Hand-over procedure and the half-hourly watch: remote_jobs.md §1.
  - Every finished step is checked with `training/check_step.py`: kinetic-energy ratio, RMSD,
    unfolded-state target, finiteness, parameter drift, glycine readout. Baseline from ff_3.0's
    epoch 5: helical glycines alpha_L 0.102 free / 0.005 restrained, H-bond margin 0.080.
- [x] **Release by checkpoint selection, not the last iterate (user, 2026-10-01; findings 1.17).**
  The midway2 gate now stops at convergence with nothing released (`gate_or_continue.sh`, as in
  the repo), and `validate_ff.sh <run> ff_3.0 <checkpoint>` releases the chosen one. Training is
  unchanged; not-converged still trains one more epoch.
- [ ] **Selection panel (rule decided 2026-10-01, user delegated; running).** The 44 L-only CATH
  domains of Charron et al. (2,823 residues, 198 glycines), all-atom ff99SB-ILDN at 300 K, against
  Upside native-start runs of each epoch-end checkpoint and of ff2.1 as control, all at T0 = 0.8.
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
  - Code `/project/trsosnic/yinhan/ff3_selection/panel.py` (copy in `scratchpad/ff3_selection`),
    jobs and watch in remote_jobs.md §1.
- [ ] **Validation.**
  - Free ensembles by native basin: helical glycines' alpha_L clearly below ff_3.0's 0.101, and
    their helical loss not below Ser/Asn's (~0.08), so propensity is not flattened. Natively
    left-handed glycines lose no more alpha_L than the Asn/Asp that sit there (~0.17).
  - glpG TM4 glycine flips against ff_3.0.
  - Peng benchmark paired against ff2.1, de novo arms included.
- Dropped: a host-guest helix benchmark and any correction calibrated to experiment. The Pace &
  Scholtz scale is an 11-system average, and the per-host data are not available to us.
- If validation fails, decide a context term then, with the evidence. Do not add one in advance.

### Phase 8 revised - glycine H-bond term and a damped side-chain step (APPROVED 2026-10-02, user)
Why (findings 1.17): with the physics map, the trainer keeps loop glycines left-handed by bending
the shared H-bond, sheet and side-chain terms, and helical glycines pay (lambda helix 3, glpG TM4).
Separately, every retrain from ff2.1 weakens helices and folds because the side-chain pair table
(31,420 coefficients, ~6% with any signal) random-walks under Adam's normalised step; putting it
back restores the panel. **Revised decisions:**
- **Glycine gets its own offsets on the three H-bond basin energies**, trained by ConDiv with
  everything else, starting at zero. `hbond_energy` takes an optional 15-entry `parameters` (12 +
  glycine's dE_alpha, dE_beta, dE_other) with a per-residue `residue_class`; a 12-entry config is
  unchanged bitwise. The SARW replica zeroes the offsets with the shared energies.
- **The side-chain (`rot`) learning rate is 10x smaller**, so per-step noise averages out while a
  consistent pull still accumulates; every file stays trained.
- Unchanged: the AWH glycine library, Peng's objective (lambda 0.3, threshold as coded), the
  protein set, the protocol.
- New run from ff2.1 (`training/ff30_glyhb`); ff30_gly is stopped when it is ready to start.
- Validation: the panel at every epoch end (helical and left-handed glycines toward all-atom, helix
  and folded fraction no worse than ff2.1), then lambda (WT toward G46A/G48A), Peng, glpG; the glpG
  seed patch must carry the new term.
- [x] Engine: optional glycine offsets in `HBondEnergy` (value, angle forces, parameter derivative).
  The class accumulation sits after the shared loops: placed beside them, -ffast-math reordered the
  shared sums and changed their last bit.
- [x] Config writer: a 15-entry `hbond.h5` writes the offsets and `residue_class`; `patch_glpg.py`
  carries them into a hybrid seed.
- [x] Trainer: field `hbg`, `zero_for_sarw`, rot learning rate / 10; `extract_ff.py`,
  `check_step.py`, README, up.md.
- [x] Tests: 12-entry config bitwise equal to the old build (energy, forces, all derivatives);
  zero offsets bitwise equal to 12 entries; offset derivatives and forces match finite differences.
- [x] Test: one real trainer step (1ga3, local): hbg gets a gradient and moves, rot steps 10x smaller.
- [x] Deployed to midway2 (bitwise parity on 1ga3, engine and binary), ff30_gly stopped, ff30_glyhb
  initialised from ff2.1 and submitted (2026-10-02 11:52).
- [ ] Panels h00, h01, ... at each epoch end; then lambda, Peng, glpG (sync `/beagle3` from `$P`
  first, and the glpG seed through `patch_glpg.py`). Lambda's helix 2 / helix 4 misorientation has
  no single side-chain cause under ff2.1 (findings, lambda update 2026-10-02), so it is not trained
  against; re-test it under ff3.0: crossing angle, the Q33-F51 dock, helix 1-helix 2 contacts.

### Phase 8 fallback - glycine context from all-atom physics (PROPOSED 2026-10-01, not approved)
Used only if Phase 8 validation fails. The user's constraints: bottom-up only for glycine, every
other residue and term as current Upside, training allowed, design change kept to the minimum.

Why a context term: with ff2.1 unchanged and the AWH map (ff30_gly step 0, findings 1.17), helical
glycines still sit in alpha_L 0.076 of the time and natively left-handed ones fall to 0.744. The
local map is physics now; what is missing is a glycine term that knows its context, which no
context-free map can supply (1.15). The smallest such term is a glycine-only offset on FF2's three
H-bond basin energies (1.15-1.16): Upside already scores each residue's H-bonds by its own basin.

- [ ] **Measure the target in all-atom (bottom-up).** In-context glycine free energies
  dG(alpha_L - alpha_R), from AWH on the glycine's (phi, psi) inside native proteins, in the map's
  own force field and temperature (ff99SB-ILDN, 300 K). Glycines: internal helical, helix C-cap,
  natively left-handed loop, beta; ~10 training proteins, ~30 glycines.
  - Our own all-atom data are peptides only (`gly_peptides/`: capped dipeptides, GGGGG, SAGAS), so
    they cannot give context. They are the no-fold baseline: the engine applies them exactly, so
    in-protein dG minus dipeptide dG is the part Upside's context terms must supply. Their AWH
    setup is reused for the biased runs.
  - **Unbiased part from public data in the same force field (findings 1.17):** Charron et al.'s
    50 CATH domains, amber99sb-ildn + TIP3P, 300 K, 4 x 0.5 us from native. Gives loop, turn and
    beta glycines directly, and a flip count for helical ones. Exclude the six D-amino-acid domains
    and the decoy frames.
  - New AWH runs only for helical glycines whose flips are too rare in those 2 us.
- [ ] **Compare with Upside first.** Upside's same dG for the same glycines, from free simulation
  at the trainer's T0. If Upside matches within error in every class, there is nothing to fix here.
  - Count only Upside frames whose surroundings are native-like, as the all-atom runs are. If not,
    Upside's extra fraying gets absorbed into the glycine offsets.
  - Fit the protein data only into the context-indexed offsets, never into glycine's map. Folded
    proteins carry placement physically: Charron's glycine phi prior, inverted from data that
    include the 50 native domains, has P(phi > 0) 0.603 against 0.500 for our dipeptides
    (findings 1.17).
- [ ] **Add the term only if the gap is there and follows the own-H-bond class.**
  - Engine: an optional per-residue class and per-class offsets on E_alpha, E_beta, E_other in
    `hbond_energy`. They go in their own dataset, outside `parameters`, so ConDiv's 12-entry layout
    and every existing config are unchanged. The config writer writes them only if `hbond.h5`
    carries them.
  - Values: fitted to the all-atom dG, not trained by ConDiv. The gradient is exact:
    d dG_j / d theta = <dV/dtheta>_alpha_L,j - <dV/dtheta>_alpha_R,j from one Upside run. Then
    one ConDiv epoch around them with everything else, and refit; stop when they no longer move.
- [ ] **Validate** as Phase 8, plus the all-atom dG of held-out glycines.

Not in this plan, on purpose: refitting any non-glycine term to all-atom data, force matching
(ff2.1 was never fitted to all-atom forces, so glycine would absorb their mismatch), and training
the offsets against the native-restrained replica (its ~0.99-in-basin target would flatten
glycine's helix-breaking propensity, 1.16).

### Phase 9 - repo cleanup for redistribution (DONE 2026-10-01, user)

`py/` and `training/` keep only what another Upside user could run. Personal, cluster-path and
campaign-specific files moved to `scratchpad/redistribution_cleanup_20261001/` (gitignored); the
cluster trees keep their own copies, so running jobs are unaffected (remote_jobs.md §1b).
- [x] Deleted `py/__pycache__`, `training/__pycache__`.
- [x] Moved from `py/`: `martini_upgrade_hybrid_args.py` (one-shot migration),
  `martini_protection_state.py` (superseded HDX criterion), `martini_inject_coverage.py` (NP-only
  retrofit), `martini_remd_concat.py` (glue for the cluster-only `run_remd.py`), `tm_score.py`
  (unused).
- [x] Moved from `training/`: `validate_ff.sh`, `patch_glpg.py`, `move_run.py`.
- [x] `training/env.sh`: midway3 `/beagle3` branch removed. `gate_or_continue.sh`: at convergence
  it stops and prints the `extract_ff.py` command for `parameters/<ff_name>`.
- [x] `training/README.md`, `train_chain.sbatch` comments, `architecture.md`, `up.md` updated.
  Verified: every kept file compiles or parses, `env.sh` resolves the repo `.venv`, all four gate
  branches (converged, retrain, max epochs, gate failure) checked in a sandbox, and no tracked
  file names a moved one.

## Known Errors / Blockers

* **Run only on midway2 broadwl** (user, 10-01). The 1.2M SU allocation reached the midway2
  scheduler at ~11:15 on 10-01. No `amd`, `beagle3` or GPU partitions: CPU jobs must not run on the
  group's GPU allocation. No midway3.
* **The held-out check in Phase 4 is ill-posed as written (findings 1.12).** The held-out mismatch
  (0.0597) is its own sample-size floor (random 46-protein training subsets: 0.0585 +- 0.0021), so it
  cannot fall with training. It must be compared against that floor, or replaced by a split-sample
  statistic.
* **Beta for the PI's sheet modelling.** No residue type has a significant beta miss at ff2.1, and
  per-pair beta has no signal in these data; per-type beta is ff2.1's sheet mixing energy, trained.
* **Left and right offsets of one central residue are nearly degenerate.** Every interior residue
  reads one left and one right map, so raising all left maps of a residue type and lowering its
  right maps changes little. The per-map prior fixes the split; watch that pair of directions in
  the gate.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
* **`/project` has ~965 G free (2026-09-30).** The glpG REMD trees hold ~1.26 T; `NP-1AO6` ~0.5 T is the reclaim.
