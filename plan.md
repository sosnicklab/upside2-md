# ff3.0: repair glpG TM4 by fixing glycine's Ramachandran map within ff2.1's own training workflow

## Project Goal

Train ff3.0 from ff2.1 with ff2.1's own training workflow (Peng et al. 2022 SI), modernised only.
The one design change is glycine's. The test is whether glpG's TM4 stays helical in the hybrid
dry-MARTINI glpG system, under the in-training force field (findings 10.19).

**Why glycine.** TM4 (residues 134-151) holds three helical glycines: GLY136, 143 and 149. Upside
applies NDRD's glycine Ramachandran map as a local energy at every glycine. That map is a
statistic over loop sites in folded proteins, so its alpha_L excess is fold selection (glycine is
placed wherever a fold needs alpha_L), not local physics: ln(aR/aL) is -1.19 in NDRD against -0.58
over all residues of the 456 training natives (findings 1.9). Data from folded proteins cannot
separate the energy from the selection (findings 1.15-1.16). So glycine's map comes from a source
where no fold selects anything, and its placement in folds has to come from the non-local terms.

**Fixed by the user.** ff2.1's training data and workflow are kept, so the old and the new force
field compare directly (findings 10.17). A drift is attributed against a control trained the same
way without the glycine change, never fixed by retuning the workflow. TM4 is judged only by the
hybrid glpG test; the soluble selection panel, the H-bond margin and the training set's composition
are not TM4 evidence (findings 10.19).

## Architecture and Key Decisions

**Trainer.** `training/ConDiv.py` is O. Kleinmann's Python 3 port of Peng's FF2 dual-target trainer,
restored to the SI (findings 9v). Per protein and step it runs:
- one native-restrained replica, 12 free replicas at T 0.8 to 1.1 and one SARW replica;
- 8000 time units, with the second half analysed;
- contrast NSE + 0.3 * DSE;
- 456 proteins in 19 minibatches of 24.

Every difference from the port is justified in the file's docstring.
- **dt 0.009** (user, 10-06): Upside's standard step, as in master, the FF1 trainer and the Peng
  benchmark. The port's 0.015 destroyed free replicas, and dt 0.009 has destroyed none in
  2,136 protein-steps (findings 1.21).
- **Side-chain learning rate 10x smaller** (user, 10-02): the 31,420-coefficient pair table
  random-walks under Adam's normalised step (findings 1.17).
- **Learning rates otherwise the port's** (`alpha_scale` 0.5). This makes `hb`, `dhb`, `sheet` and
  the burial groups twice the rate the SI implies (fine-tuning x 0.25; findings 1.23). The SI-rate
  arm below tests it.

**ff2.1's parameter set** (user, 09-24):
- rot: pair, coverage, hydrophobe;
- burial: 20 x 3 sigmoid parameters and 400 weights;
- the backbone-desolvation scale;
- the three H-bond energies and the second-H-bond term;
- the 20 sheet values.

Left untrained, as in ff2.1, because the engine returns no derivative for them: the backbone term's
center, sharpness and hbond weight, and `hbond.h5` entries 4-11.

**Round 4's glycine designs** (user, 10-04/05; findings 1.19):
* **bio:** glycine's maps are fitted bottom-up to BioEmu's plain-MD octapeptides (ff99sb-ildn,
  300 K; no fold present) by iterative Boltzmann inversion in Upside, with ff2.1's other terms, at
  T_up = 1 (up.md 2.8).
  - Fitted entries: pooled GLY|X, GLY|left|GLY, GLY|right|GLY, GLY|right|PRO.
  - Cells with fewer than 10 all-atom frames keep NDRD's coil top and barrier.
  - Stored as `parameters/common/rama31.dat` and frozen in ConDiv. Every other map stays NDRD.
* **gdepth:** ff2.1's own NDRD library with one pooled glycine alpha_R / alpha_L depth pair on the
  37 GLY|X maps.
  - It starts where the maps hold BioEmu's alpha_R / alpha_L weights; that start is the only place
    BioEmu enters.
  - It is trained per epoch by a damped Newton step on free against restrained populations, with a
    Gaussian prior of 1 nat on NDRD.
  - The claim to test: training on the original data alone keeps TM4 stable.

**H-bond and sheet frozen, as a second arm** (user, 10-07). `ff30_bio_fz` runs with `hb`, `dhb`,
`hbg` and `sheet` at learning rate 0, so `hbond.h5` and `sheet` stay ff2.1's. So does its control
`ff21_ctrl_fz` (ff2.1's library; against `ff21_ctrl_dt009` it tests the H-bond drift with no glycine
change). The gdepth twin was cancelled before it started (10-07 13:12, user left it to judgement);
its queue start was a day away, and the SI-rate arm takes the gdepth rate question. Why the freeze:
* Round 1 kept TM4 stable while its FF1-form port trained neither term (findings 9e). That round
  also had a symmetric glycine map, no DSE term and FF1's burial, so it cannot credit the freeze.
* The DSE term pulls every H-bond energy weaker: SARW has no H-bonds (findings 1.17, 1.23).
* In the hybrid, glpG has no burial or backbone-desolvation node, so a weaker E_alpha reaches TM4
  with nothing to offset it.

Neither Jumper's nor Peng's soluble training froze `hb` (findings 1.23). The frozen values are
Peng's ConDiv-trained ones, and ff2.0 to ff2.1 kept `hbond.h5` and `sheet` byte-identical.
Frozen means frozen for the whole run. Releasing them later would return to the joint fixed point,
because ConDiv's step depends only on the current parameters.

**SI-rate arm** (user, 10-07). `ff30_bio_si` and `ff30_gdepth_si` are the two designs with every
group but rot at the SI's learning rate (hb 0.005, as FF2 was trained; findings 1.23), otherwise
their port-rate twins. With the full-rate and frozen runs, bio has three H-bond rates on one
design: full, half and zero.

**Control** (user, 10-06): `ff21_ctrl_dt009` has ff2.1's own `rama.dat`, the same trainer, and
target 38 (the half-trained point), with no gate. It separates the workflow's own drift from the
glycine change.

**The TM4 test** (user, 10-04/06; findings 1.22, 3.11):
- Setup: glpG 79HIS in the POPE/POPG dry-MARTINI hybrid seed. `patch_glpg.py` writes the
  checkpoint's Rama map, H-bond energies, pair table and both coverage tables into it. The membrane,
  ions and coordinates stay bit-identical.
- Run: T 0.80, 4000 tu, Verlet dt 0.009, 12 seeds.
- Counts: a seed is unwound if its last-block TM4 helix fraction is below 0.90. It is flipped if a
  TM4 glycine has phi > 0 in more than 0.25 of the last block.
- Each checkpoint is compared with its own start by Fisher's exact test. A frozen run is also
  compared with its unfrozen twin.
- Recipe: remote_jobs.md "Round-4 watch" step 4.

**Release by checkpoint selection** (user, 10-01/05).
- At convergence the gate validates its final epoch-end checkpoint automatically:
  - `validate_ff.sh` runs Peng's 32 arms and glpG's 4 REMD variants;
  - the candidates are `ff_3.0_bio`, `ff_3.0_gdepth`, `ff_3.0_bio_fz`, `ff_3.0_bio_si` and
    `ff_3.0_gdepth_si`.
- A run that reaches 13 epochs unconverged stops for review.
- The selection panel (44 all-L CATH domains against all-atom) chooses. Peng and glpG judge.
- Which candidate becomes ff_3.0 is the user's.

**Retired designs** (one line each; detail in findings):
* Full-map glycine row (Phase 2): failed the convergence gate at every checkpoint.
* Per-pair basin offsets on the NDRD maps, the ff_3.0 released 09-30: favoured alpha_L at every
  glycine; noise-limited per pair (1.12, 1.14-1.15).
* Physics glycine map (AWH dipeptides) plus glycine H-bond offsets, round 3 (`ff30_glyhb`): panels
  lost helix; the offsets went left-handed (1.17).
* Ting et al.'s product combining rule: breaks GLY|GLY symmetry (1.13).
* A design through H-bond terms for glycine (user, 10-04): a map defect is fixed on the map.
* Pre-proline correction: parked; no simulation failure behind it (1.16).

## Execution Phases

### Done or closed (detail in findings and progress.md)
- Phase 1: the trainer validated on ff2.1 (09-24; 9v).
- Phase 2: the full-map glycine row, cancelled 09-28 (9w).
- Phases 3-5: the basin-offset ff_3.0, released 09-30 and withdrawn the same day (1.14-1.15).
  - Its glpG validation ran on live seeds that kept the FF1-form ff_3.0's coverage tables (3.11),
    so its "TM4 glycines flipped" result is of a mixed force field.
- Phase 6: pre-proline, parked.
- Phase 7: the glycine handedness probe. From equal depth the data pull glycine back toward alpha_L
  (1.15).
- Phase 8: round 3, `ff30_glyhb`, stopped 10-04 (1.17).
- Phase 9: repo cleanup for redistribution, 10-01/02. `training/` has five files; the cluster keeps
  the old layout (remote_jobs.md "Where things are").
- Phase 10: the local glpG test of round 3. Its TM4 counts used pre-fix inputs and are invalid
  (3.11).

### Phase 11 - round 4 (TRAINING on midway2 broadwl since 2026-10-06; frozen arm since 10-07)

Done:
- [x] BioEmu-fitted library in `rama31.dat`. Its own pass matches BioEmu's X-G-Y glycine within
  0.01 per basin (scripts and data in `checks/gly_bioemu_map/`).
- [x] Push probe from the library (72 proteins): d -0.014 [-0.041, +0.012]. The training data come
  to rest near BioEmu's L-R difference.
- [x] Non-glycine maps kept NDRD (findings 1.19).
- [x] gdepth trainer and its start, with tests (round 2's offset hooks ported).
- [x] dt 0.015 runs stopped 10-06 09:40 and restarted at dt 0.009, otherwise identical (initial
  force fields byte-identical).
- [x] Control `ff21_ctrl_dt009` submitted.
- [x] Frozen arm submitted. ff30_bio_fz's step 0 keeps `hbond.h5` and `sheet` md5-identical to
  ff2.1's. ff30_gdepth_fz was cancelled unstarted (10-07 13:12) and resubmitted unchanged at 21:50
  (user: run every frozen arm, though it departs from the group's workflow); ff21_ctrl_fz was
  cancelled with it and resubmitted (user).
- [x] SI-rate arm `ff30_bio_si` and `ff30_gdepth_si` built, initialised and submitted (10-07
  13:15). Initial force fields and gdepth's round-0 library are byte-identical to their twins'.
  The trainers differ only in the rate lines (`checks/si_init_20261007`).
- [x] Automatic validation deployed (gate, `validate_ff.sh`, `after_training.sbatch` with
  candidate names).
- [x] `patch_glpg.py` coverage fix (findings 3.11), local and on midway2.
- [x] Pre-fix TM4 runs moved aside; local inputs re-patched.

- [x] midway2 node exclusions as one law (user, 10-08 09:40; built and tested by 09:40): `/project/trsosnic/yinhan/slurm/
  midway2.args`, read by every submission path (slurm.args symlinks, `submit_remd.sh`,
  `bench.sbatch`), a `~/bin/sbatch` wrapper for the rest, `update_pending.sh` for queued jobs;
  midway2-0027 added (remote_jobs.md §0d).
- [x] **ff30_gdepth_dt009 converged and is the local ff_3.0** (gate 10-09 08:54, findings 1.26; user,
  10-09 10:40: "pull ff30_gdepth_dt009 down locally as ff3.0", to show). `parameters/ff_3.0/` holds
  the six released files, md5-identical to `$P/parameters/ff_3.0_gdepth` and to `d9_03`. Its glycine
  library is its own `rama.dat`: the example scripts read `parameters/common/rama.dat`, and
  `py/martini_prepare_system.py:1922-1924` takes rama, sheet and hbond from ff_2.1, so a run on
  ff_3.0 passes its files explicitly or is patched with `patch_glpg.py --ff parameters/ff_3.0`.
  `martini.h5` and `membrane.h5` are not trained and stay ff_2.1's. Its validation runs on.

Next:
- [ ] **TM4 on fixed inputs, 12 seeds; the frozen comparisons at 24** (user, 10-08 09:10).
  bz_00 against b9_00 and cz_00 against c9_00 run seeds 13-24 as well. A 24-seed set runs as two
  12-seed halves (queue line `<tag> 13 24`), and `tm4_compare.py` reads seeds 1-24 from
  `runs_cov/` and `runs/` together, with each set's own n in the Fisher table.
  - Done: bio_start 5 of 12 unwound, b9_00 7 of 24, b9_01 6, b9_02 8, bs_00 6, bs_01 5, bz_00 8 of 24, bz_01 6, bz_02 4, ds_00 9, d9_02 4, dz_00 4, ff21_released 9,
    c9_00 16 of 24, c9_01 10, fp_e00 5, gdepth_start 7, d9_00 5, d9_01 6 (findings 1.24). No pair is
    resolved. bz_00's first half flipped none of 12 (one-sided p 0.047); its seeds 13-24 flipped 4,
    bio_start's count, and against b9_00 at 24 seeds each it is 8 against 7 unwound and 4 against 5
    flipped, so freezing H-bond and sheet shows no TM4 change at epoch 0. All three trainings tested
    at two epoch ends lean better or even at epoch 0 and worse at epoch 1, the control included.
    bs_00 leans worse than its twin b9_00; bs_01 leans better than b9_01 (2 flipped against 5,
    p 0.37), the first half-trained end not leaning worse than its epoch 0. bz_01 leans as b9_01
    does (6 unwound each, TM4 0.846 against 0.833): freezing H-bond and sheet does not stop the
    epoch-1 lean, and in the glpG test bz differs from bio_start only in the side-chain tables.
    b9_02 keeps b9_01's lean at epoch 2 (TM4 0.830 against 0.833, 8 unwound against 6, p 0.68).
    ds_00 leans worse than its twin d9_00 (9 unwound against 5, p 0.21; TM4 0.775, the lowest mean).
    c9_00's seeds 13-24 flipped 8 of 12 against seeds 1-12's 1 (same input, p 0.009 between
    halves), so its 12-seed lean was chance; at 24 seeds it is even with ff2.1 (16 and 9 against 9
    and 5 of 12, TM4 0.843 against 0.827). A lean in one 12-seed set is a reason for 24 seeds.
  - [x] **TM4 secondary structure by DSSP** (user, 10-09 02:15: "make sure you are testing the
    stability of the secondary structure of TM4"; findings 1.25). The count scores helix from
    phi/psi boxes alone; `tm4_local.py` and `tm4_compare.py` now also give DSSP alpha-helix and
    any-helix fractions of 135-151 per block, alpha per residue, and a Mann-Whitney on each, with
    the model's own carbonyl O and mdtraj's DSSP. All 28 tables regenerated, earlier lines
    unchanged. DSSP alpha is 0.10-0.22 below the box in every set, the order of the sets is nearly
    the same, and nothing resolves.
  - [x] **DSSP alpha-helix is the primary TM4 test** (user, 10-09 03:20): each seed's last-block
    alpha of 135-151, two-sided Mann-Whitney, resolved at p < 0.05, with a bootstrap interval;
    the dihedral counts are secondary. All tables carry it. No in-training checkpoint is resolved
    from its start; d9_02 leans highest (0.765 against gdepth_start's 0.657, p 0.58). bz_02 against
    its twin b9_02 is the first resolved primary test, only just: 0.856 against 0.694, p 0.046, with
    the bootstrap interval reaching -0.010, one of 32 tables read; against bio_start p 0.30. At 12
    seeds a true 0.10-0.15 difference is detected 18-31% of the time (findings 1.25).
  - Locally, on the Mac Studio, which keeps the jobs (user, 10-08 09:15; another computer stands
    by until told to take over, remote_jobs.md "Handoff"), in this order: d9_03 (from 10:46, the released ff_3.0_gdepth), bs_02, then every new epoch end, half-trained (step 37) first;
    cz_00 at 24 seeds. Its engine and inputs reproduce the MacBook Pro's runs.
- [ ] **Answer the user's question: which run moves TM4 toward stable.** Compare each checkpoint
  with its own start, each frozen run with its unfrozen twin at the same step, and read the direction
  across epoch ends. With 12 seeds only large changes resolve (findings 1.22); run 24 seeds where a
  difference falls between.
- [ ] Every step: `check_step.py`, the KE scan, the margin and glycine readouts. A trained-H-bond run
  whose margin falls below +0.10 holds for the user; ff30_bio_dt009 has been on hold since step 14
  and trains on.
- [ ] Each epoch end: selection panel (`submit_new.sh`).
- [ ] At a converged gate: automatic validation. The user chooses ff_3.0.

### Phase 12 - poly-Gly: all-atom Ac-(Gly)20-NHMe as Upside's reference (STARTED 2026-10-05, PI's suggestion)
With neutral caps the chain is achiral, so alpha_R = alpha_L exactly at equilibrium. It shows the
glycine map's own handedness with nothing else present.
- Setup: amber99sb-ildn / TIP3P at 300 K on midway2. Files and README in
  `/project/trsosnic/yinhan/polygly/`; Upside runs and analysis in the Mac's `scratchpad/polygly/`.
- Steps: collapse (4 replicas, discarded), then a production box sized from the measured diameter,
  then 4 x 1 us production (~12k core-hours). It is converged when the blank is zero within its
  error and the replicas agree on Rg.

Status:
- [x] Collapse: job 49186415 COMPLETED 10-07 12:48 at 19.79-20.00 ns of the 30 ns cap, at 14 ns/day
  per replica. Every replica was compact (Rg < 0.75 nm) within 0.4-3.7 ns. Second-half Rg was
  0.67-0.73 nm, with reopenings up to 19% of frames. The minimum periodic-image distance was
  5.2 nm (`polygly/collapse_analysis/`).
- [x] Production boxes (user, 10-07; `scripts/build_prod.sh`): each replica's last frame in a 7.5 nm
  dodecahedron. 7.5 nm is the largest diameter after the first ns, 5.14 nm, plus 2.4 nm; the
  plan's second-half rule gives 6.9 nm. About 28,950 atoms and 9,600 waters each, against 64,976
  before. grompp is clean.
- [x] Production chain `prod.sbatch` (minimise, 500 ps NPT with seeds 20261021-24, then 1 us,
  chained 36 h links, frames every 10 ps). Submitted 10-07 14:10. Expected ~30 ns/day per replica,
  about a month to 1 us; the first link measures it.
- [x] Upside G20 under ff21_released and bio_start, 8 x 200,000 tu: both chiral, ln(aR/aL) -0.72 and
  -0.45 (findings 1.20).
- [ ] First link: equilibration logs, ns/day, `gmx mindist -pi` above 1.0 nm. Then the blank at each
  ~200 ns (step 5).
- [ ] Analysis script.
- [ ] Round-4 epochs.

## Known Errors / Blockers

* **Every glpG TM4 count before 2026-10-07 11:22 is invalid** (findings 3.11). That covers:
  - every local 3- and 12-seed set (Phase 10, Phase 11);
  - the 09-30 ff_3.0 glpG validation (Phases 3-5);
  - any other glpG run from the live seeds under a force field other than the FF1-form ff_3.0.

  Each of these used the FF1-form ff_3.0's coverage tables. `patch_glpg.py` is fixed (local, `$P/training/`,
  `tm4_local/scripts/`, md5 `70589119...`); the automatic validation uses the fixed copy. Another
  computer must update its local copy before use (remote_jobs.md "Setting up the next computer").
* **Total-potential jumps in local glpG runs** (findings 1.22): 13 of 105 pre-fix seeds, one blow-up.
  Cause not identified. Check the fixed-input runs for them; do not drop frames.
* **midway2-login2's `/project` hangs** (since 10-07 10:48). Route `/project` work through login1
  (remote_jobs.md §8).
* **The cluster's `training/` keeps the pre-10-02 layout** while round 4 trains. Do not sync the
  repo's `training/` over it.
* **`run_remd.py` rolls back NaN replicas and continues**, which the NO GUARDS rule forbids. It is
  kept unchanged so glpG stays comparable across campaigns. Whether to remove it is the user's call.
* **CPU work runs on midway2 broadwl** unless the user names midway3 caslake. No `amd`, `beagle3`
  or GPU partitions, and not `broadwl-lc` (its nodes cannot see `/project`).
* **lambda's helix 2 / helix 4 misorientation** is a known ff2.1 limitation with no single
  side-chain cause (findings 11d). It is re-tested after ff3.0, not trained against.
