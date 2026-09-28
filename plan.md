# ff3.0: ff2.1's own training workflow, with the Ramachandran maps trained as per-pair basin offsets

## Project Goal

Train ff3.0 from ff2.1 with **exactly ff2.1's training workflow**, modernised only, plus one project
addition: every Ramachandran library map gets a small set of trained basin offsets, so the map
supplies only what Upside's other terms do not already explain at the training natives. The result
is released as `parameters/ff_3.0`, which is absent from the tree until then.

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
  the shape instead (44 deg transitions, no flat interior in beta). 5 offsets per map, 6 for a
  central glycine: 4,234 over the 840 maps.
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
  round; overfitting shows as held-out agreement stalling while training agreement improves.

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

### Phase 4 - train ff3.0 from ff2.1 with basin offsets (RUNNING since 2026-09-28 01:33, midway2)
- [x] Initialise from ff2.1 with all offsets zero (`training/ff30_basin`); first target 76 steps,
      then one epoch at a time up to 13 epochs. Step 0 done: 24/24 workers, 4,688 map-sites, free
      and native basin populations within ~0.02 of each other at ff2.1
- [x] First offset update at the end of epoch 0 (step 19). Its noise-driven steps were found at the
      status check; the chain was stopped during step 20, steps 19-20 set aside
      (`ff30_basin/rewound_20260928_0851/`), round 1 recomputed from the same epoch-0 data with the
      MAP step (free-native mismatch 0.025 training, 0.060 held out), and training resumed from
      step 19 at 08:56
- [x] Release path dry-run on the round-1 checkpoint: extraction, glpG round-trip gate (3.6e-15),
      live-seed patch, and `sbatch --test-only` for a Peng arm, a glpG chain, the gate and a
      training link. At convergence the gate installs the rotamer-BP fix, releases ff_3.0 and
      submits the Peng benchmark and the glpG chains unattended
- [x] GLY|GLY made symmetric in the map the engine gets (findings 1.11): sheet entry symmetrized,
      probability mean for coil and sheet, `extract_ff.py` and `verify_rama_basin.py` (new step 4,
      the per-residue map) check both; verifier PASSES on the run. Chain stopped during step 29 on
      2026-09-28 12:25, steps 19-29 (29 unfinished) set aside (`ff30_basin/rewound_20260928_1225/`), round-1
      library rewritten from the same offsets (identical outside GLY|GLY), resumed from step 19.
      Epoch 0 and round 1's offsets came from the old GLY|GLY maps; round 2 onward uses the new ones
- [ ] Held-out agreement per round; stop and review if it stalls while training agreement improves
- [ ] Release `parameters/ff_3.0` and copy it into the local repo

### Phase 5 - validation (NOT STARTED)
- [ ] Peng benchmark, 16 proteins x native/de novo, scored on the last third, paired against ff2.1
- [ ] glpG, four variants: helix stability over time, TM4 above all

## Known Errors / Blockers

* **Left and right offsets of one central residue are nearly degenerate.** Every interior residue
  reads one left and one right map, so raising all left maps of a residue type and lowering its
  right maps changes little. The per-map prior fixes the split; watch that pair of directions in
  the gate.
* **The local Mac `obj/upside` traps (SIGTRAP, exit 133) at exit whenever Monte Carlo pivot moves
  are on**, even for one system. The midway2 binary runs them cleanly, so trainer tests run there.
* **Do not use `broadwl-lc`.** Its nodes are `noib` and cannot see `/project`.
* **`/project` has ~445 G free.** The glpG REMD trees hold ~1.26 T; `NP-1AO6` ~0.5 T is the reclaim.
