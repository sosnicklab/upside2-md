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

### Phase 5 - validation (RUNNING since 2026-09-30 02:25, midway2)
- [ ] Peng benchmark, 16 proteins x native/de novo, scored on the last third, paired against ff2.1
- [ ] glpG, four variants: helix stability over time, TM4 above all

### Phase 6 - pre-proline (DECISION PENDING, findings 1.14)
- [ ] Decide on the X|right|PRO offsets: they have no leverage through the mixture and their target
      signal is composition; proposed: drop the 40 from the trained set at the next revision
- [ ] Before any rule change: reweight existing ff_2.1 trajectories to the right-only pre-proline
      rule and measure helical and extended pre-proline alpha_R (no new simulation)

### Phase 7 - glycine handedness probe (RUNNING, user request 2026-09-30; findings 1.15)
Which way do the data pull glycine when alpha_R and alpha_L start at equal depth? One epoch
branched from the ff_3.0 checkpoint (step 114) into `training/ff30_glyprobe`, every GLY|X map's
alpha_R and alpha_L offsets set so the two basins hold equal probability, everything else as
released. The answer is the data term of the offset update (free minus native basin counts over
the epoch), read apart from the prior, which pulls every offset back toward NDRD.
- [x] Branch checkpoint, verify (only GLY|X coil maps differ from round 6; each has aR = aL), submit
      one epoch with no gate (job 49133133, 2026-09-30 13:28; analysis validated on the ff_3.0 run)
- [ ] Read the data pull per map and in aggregate, extended and helical glycines apart; record.
      Partial (12 of 19 steps): toward alpha_L, -0.031 [-0.046, -0.016], 32 of 38 maps (findings
      1.15); confirm on the full epoch

## Known Errors / Blockers

* **The held-out check in Phase 4 is ill-posed as written (findings 1.12).** The held-out mismatch
  (0.0597) is its own sample-size floor (random 46-protein training subsets: 0.0585 +- 0.0021), so it
  cannot fall with training. It must be compared against that floor, or replaced by a split-sample
  statistic.
* **Pre-proline is not fixed by the offsets, and ConDiv cannot measure it (findings 1.14).** The
  mixture caps the right map's effect (pre-PRO alpha_R of the rama term 0.174 at best to 0.144), and
  six rounds left the free-native gap at +0.023. That gap is the class's composition (89% extended
  residues, which visit alpha_R in every class); within matched native conformations pre-proline
  residues already visit alpha_R less than others. The physical defect, three times NDRD's own
  pre-proline alpha_R through the mixture, lives in loops and unfolded chains, which the
  native-restrained target does not sample. The candidate fix is a combining rule, right-map-only
  for residues followed by PRO (GLY|GLY stays exact), justified by physics and NDRD data, not by
  the training target; its risk is the 7% of pre-proline residues that are helical.
* **Beta for the PI's sheet modelling.** No residue type has a significant beta miss at ff2.1, and
  per-pair beta has no signal in these data; per-type beta is ff2.1's sheet mixing energy, trained.
* **Left and right offsets of one central residue are nearly degenerate.** Every interior residue
  reads one left and one right map, so raising all left maps of a residue type and lowering its
  right maps changes little. The per-map prior fixes the split; watch that pair of directions in
  the gate.
* **The local Mac `obj/upside` traps (SIGTRAP, exit 133) at exit whenever Monte Carlo pivot moves
  are on**, even for one system. The midway2 binary runs them cleanly, so trainer tests run there.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
* **`/project` has ~965 G free (2026-09-30).** The glpG REMD trees hold ~1.26 T; `NP-1AO6` ~0.5 T is the reclaim.
