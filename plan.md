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

**TM1** (user, 10-09). glpG's TM1 loses structure under every force field, ff2.1 included
(findings 4.6). It is addressed in the hybrid itself, by correct physics and without retraining
(Phase 13).

**ff3.1** (user, 10-09, after talking with Jumper). Train the relative depth of every basin of every
Ramachandran map, not only glycine's pair, as offsets from NDRD under an explicit L2 prior, with no
new input simulations (route 1 of findings 1.28; Phase 14).

**ff3.2** (user, 10-15). ff3.1's `data` training started from ff3.0 instead of ff2.1 (Phase 17).

**The basin-depth path, not the BioEmu map** (user, 10-11 23:50: "we decided to go with basin depth
training path (as it is expandable to all residues)"). The BioEmu-map line (ff30_bio, ff30_bio_dt009,
ff30_bio_si; released ff_3.0_bio, ff_3.0_bio_si) is dropped and its validation cancelled. **ff3.0 is
ff_3.0_gdepth** (user, 10-12 00:25: "since ff3.0 is working, it is the only branch we need (and its
expansion to all residues)"): its SI-rate twin ff30_gdepth_si was stopped at step 72 before its gate
and deleted, and ff3.1 (Phase 14) extends the depth training to every residue. The
H-bond-frozen arms, polygly (stopped) and the local TM4 glycine-flip test (retired) went the same
day, and the BioEmu line and the frozen arms were deleted with their inputs and derived products
(10-11 23:45, user; remote_jobs.md §1 "Finished runs"): on the fixed side-chain/lipid table ff2.1 holds TM4 (findings
4.15). TM4 HDX has the local readout only; its REMD test was cancelled 10-14 (Phase 16).

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

### Phase 11 - round 4 (TRAINING on midway2 broadwl since 2026-10-06; BioEmu and frozen arms dropped 10-11, user)

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
  cancelled with it and resubmitted (user). **All three frozen arms stopped 10-09 23:09 (user):**
  the dt 0.009 training already gives a force field that keeps TM4 stable, so the frozen
  comparison no longer decides anything (ff30_bio_fz at step 58, ff30_gdepth_fz at 45, ff21_ctrl_fz
  at 10; remote_jobs.md §1). Their TM4 and panel results stand (findings 1.21-1.25).
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
- [x] **ff30_bio_si converged** (gate 10-10 02:25, steps 58-76, lowest dhb p 0.0125) and was
  released as `ff_3.0_bio_si` (= bs_03; the frozen BioEmu library, SI rates); its 32 benchmark arms
  and 4 glpG chains are queued (remote_jobs.md §1). Its TM4 is level with ff_3.0_gdepth's (0.739
  against 0.726, p 1.00; findings 1.31). Choosing between them stays the user's.

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
    from its start; d9_02 leans highest (0.765 against gdepth_start's 0.657, p 0.58), and the
    released d9_03 is 0.726 (p 0.64; against d9_02 p 1.00; findings 1.26). bz_02 against
    its twin b9_02 is the first resolved primary test, only just: 0.856 against 0.694, p 0.046, with
    the bootstrap interval reaching -0.010, one of 32 tables read; against bio_start p 0.30. At 12
    seeds a true 0.10-0.15 difference is detected 18-31% of the time (findings 1.25).
  - Locally, on the Mac Studio, which keeps the jobs (user, 10-08 09:15; another computer stands
    by until told to take over, remote_jobs.md "Handoff"), every new epoch end as it comes (dz_01 done 20:45; the queue is empty), half-trained (step 37) first;
    cz_00 at 24 seeds. Its engine and inputs reproduce the MacBook Pro's runs.
- [ ] **Answer the user's question: which run moves TM4 toward stable.** Compare each checkpoint
  with its own start, each frozen run with its unfrozen twin at the same step, and read the direction
  across epoch ends. With 12 seeds only large changes resolve (findings 1.22); run 24 seeds where a
  difference falls between. On the fixed side-chain/lipid table (10-10) ff2.1 itself holds TM4
  (0.944 against 0.610 on the old table, p 0.002), every start reads alike and no released
  candidate improves on it (findings 4.15); what follows from that is the user's.
- [ ] Every step: `check_step.py`, the KE scan, the margin and glycine readouts. A trained-H-bond run
  whose margin falls below +0.10 holds for the user; ff30_bio_dt009 was on hold from step 14 to its
  end.
- [ ] Each epoch end: selection panel (`submit_new.sh`).
- [ ] At a converged gate: automatic validation. The user chooses ff_3.0. Converged: ff_3.0_gdepth
  (10-09 08:54), ff_3.0_bio_si (10-10 02:25), ff_3.0_bio (10-10 16:46; findings 1.26, 1.31, 1.32);
  ff30_gdepth_si still trains.

### Phase 12 - poly-Gly: all-atom Ac-(Gly)20-NHMe as Upside's reference (STOPPED 2026-10-11, user; started 10-05, PI's suggestion)
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
- [x] Links 1 and 2: 33.8-34.0, then 32.4-32.6 ns/day; `gmx mindist -pi` 2.31-3.57 nm, above 1.0 nm.
- **Stopped by the user at ~97 ns per replica (10-11 23:05):** the question it served, the glycine map
  behind TM4's glycine flips, was settled by the side-chain/lipid table fix (findings 4.15: ff2.1 holds
  TM4 on the fixed table, 0 of 12 flipped). The blank, the analysis script and the round-4 epochs are
  not pursued. The data are kept (`polygly/`, 213 MB), and the chain resumes from `prod.cpt` if the
  glycine map's own handedness (Upside G20 ln(aR/aL) -0.72, findings 1.20) is ever taken up.

### Phase 13 - glpG TM1: hold it in the hybrid by correct physics (STARTED 2026-10-09, user)
Hold glpG's TM1 in the hybrid without retraining the force field (user, 10-09): no scaling or
disabling of BB-env or SC-env, no tuned parameters, no guards.

**Diagnosis** (done; findings 4.6-4.7). TM1's middle (32-38) sits at the bilayer centre and loses
structure in every force field, ff2.1 included. Lipid gain is not the driver. Opening it costs only
Upside's soluble-calibrated H-bond energy: the hybrid has no cost for an unpaired backbone NH or CO
in the acyl core, and its MARTINI backbone typing is fixed at the starting secondary structure.
Upside's own implicit-membrane model has that cost (`hb_membrane_potential`), and on the stored
frames it would add +7 to +10 E_up to the opening.

**Approach** (user chose A, 2026-10-09 13:45):
* **A. Add Upside's membrane backbone H-bond term, its H-bond-state-dependent part only.**
  - E = sum over backbone donors and acceptors of (1-p)^2 [f_unpaired(z) - f_paired(z)], p from the
    hybrid's existing `protein_hbond`, f from `parameters/ff_2.1/membrane.h5` hb_energy unchanged.
    The paired baseline is left out because MARTINI's helix-typed BB already describes a paired
    backbone against lipid.
  - First test is config-only: the engine's existing `hb_membrane_potential` node, written into a
    copy of the patched glpG input with coeff = [hb_energy[t,0] - hb_energy[t,1], 0], use_curvature
    0, the bilayer centred at z = 0 (it stays within 1-2 A over 4000 tu). No C++ change, master
    parity untouched.
  - Half-thickness from our own bilayer, not borrowed: the glycerol-ester (GL1/GL2) mid-plane, the
    hydrocarbon-core boundary Upside's membrane thickness denotes; the C1 plane (12.4 A) would give a
    smaller penalty (`tm1_hbmem.txt` has the sensitivity).
  - Production would need the reference to follow the bilayer's centre (a node reading the lipid
    mid-plane) and a builder option in `py/martini_prepare_system*.py`; decided after the test.
* **B. Make each BB bead's MARTINI type follow its H-bond state** (N0 paired, P5 unpaired). MARTINI's
  own endpoints, but a new C++ coupling in `martini_potential` with derivatives through
  `protein_hbond`; on the stored frames it would charge about +38 E_up for the opening, far more
  than A.

Status (A):
- [x] Half-thickness 15.9 A: mid-plane of GL1 (16.42 A) and GL2 (15.40 A) over the second half of the
  local ff2.1 run; the first acyl beads sit at 12.1-12.7 A (seed: GL1 17.6, GL2 16.3).
- [x] Term built into copies of both inputs (`scratchpad/glpg_tm1_hbmem/add_hbmem.py`). Verified on
  stored frames (`verify_hbmem.txt`): E_new - E_old equals the node value and an independent
  evaluation to 1e-3 E_up; node forces match finite differences to ~1e-3. +86.5 E_up at t 0.
- [x] Test, ff2.1 + term, seeds 1-12 (finished 15:39; findings 4.8): every TM helix more helical
  (TM1 +0.11, TM4 primary +0.18), but TM4's backbone tears in 2 seeds and the protein potential
  jumps > 3000 E_up in 5. Not acceptable as it stands.
Second lead, transient backbone excursions (findings 4.9):
- [x] Dense-frame replay of ff2.1 seed 1 to t 420 (`kick/replay.sh`): bitwise reproduction; onset
  between t 407.70 and 407.97 at MET34's backbone; no MARTINI pair spike before it; no term's force
  on 33-35 abnormal at the stored frames (`kick/force_by_node.txt`). The impulse lies inside 10
  steps that no stored frame shows.

**Local plan on the Mac Studio** (user, 10-09 21:40: resolve TM1 in the local ff3.0, d9_03 =
ff_3.0_gdepth, by local simulation). The excursions come first: TM1's middle opens at one (findings
4.9), the unexplained TM4-test jumps may share their cause (Known Errors), and so may option A's
tears. Option A on ff3.0 is the second line, and runs the same night on the idle cores (user, 10-09
21:55): step 5's set with `obj/upside` while step 2 runs; any tear in it is replayed exactly with
the probe later, and the set is rerun if step 3 changes the engine.
* **Readout, fixed before any new run.** TM1 primary: each seed's last-block DSSP alpha of 29-48,
  two-sided Mann-Whitney against the 12 local d9_03 seeds (`ff3_local_test/runs/`, not rerun),
  resolved at p < 0.05, with a bootstrap 95% interval. Also TM1 32-38 i->i+4 H-bonds, TM4's primary,
  and health: frames with a C-N above 2 A, total-potential jumps above 3000, KE/1.5kT, Rg and the
  protein's depth. TM4-test flags unchanged (`run_glpg.sh`: T 0.80, 4000 tu, dt 0.009).
* **Cores.** The TM4 queue comes first (12 cores a set; bs_03 expected ~02:00 10-10). TM1 work uses
  at most 16 cores minus the running TM4 seeds, and never replaces `obj/`, the TM4 driver's binary.
1. [x] Set up `scratchpad/glpg_tm1/` from midway2's `checks/tm1_hbmem_20261009/` (script paths
   fixed); TM1 baseline table for d9_03 and ff21_released from the local runs.
2. [x] Localise the impulse in ff3.0 (findings 4.10): the side-chain/lipid table
   (`martini.h5` `sc_table`) is sampled on 24 directions and one bead frame, and its one-node ridges
   give single-evaluation spikes up to 33,260 E_up/A on CA34/CA35. Method: a temporary probe build: a copy of `src/` in
   `scratchpad/glpg_tm1/src_probe`, built into its own `obj_probe/` (the tracked `src/` is not
   touched). An environment-gated ring buffer keeps the last ~40 steps' positions and momenta and
   writes them out when a protein backbone atom moves more than a set distance in one step; it only
   reads. First check that it reproduces a stored d9_03 log frame for frame. Then run d9_03 s10 (TM1
   excursion by t 140), s7 (t 250), s8 (t 260) and s11 (t 440) to just past their events, and
   evaluate every potential node's force on 33-36 at each written step, with the rotamer's
   side-chain/lipid 1-body table split out and the integrator stage that carries the spike named.
3. [x] Fix at the root, tested locally: the table built with converged azimuthal sampling (600 directions x 12
   frames, `kick_ff30/sc_table_vectorized.py`, the builder's own computation vectorised), tested in
   a copy of `patched/79HIS_d9_03.up` with only `rotamer_full_energy_eup` replaced.
   `parameters/ff_2.1/martini.h5`, `py/martini_build_tables.py` and the cluster stay unchanged until
   the user decides: the round-4 TM4 sets all ran on the old table.
4. [x] Fixed ff3.0, 12 seeds (findings 4.13): TM1 0.931 against 0.685, p 0.002, resolved; TM4
   0.924 against 0.726, p 0.10; no evaluation above 1000 E_up/A; one C-N excursion frame against 33.
5. [x] Option A on ff3.0, 12 seeds each. On the old table (findings 4.11): TM1 0.818 against
   0.685, p 0.073, not resolved. On the fixed table (findings 4.14, `runs_sc600x12_hbmem/`, probe
   on): TM1 0.964 (p < 0.001) and TM4 0.959 against 0.726 (p 0.009), both resolved; against the
   table alone TM1 +0.033 (p 0.049) and TM4 p 0.21; no C-N excursion, no jump; three evaluations
   above 1000 E_up/A, protein side-chain contacts outside TM1 and TM4. ff2.1 + term's tears are the
   table's kicks on MET142 (findings 4.8). The literature (findings 4.12) supports the term's
   magnitude and forbids combining it with option B (double counting).
6. [ ] Delete the probe build and its source. The report is done (10-10): the table fix is in
   production (Phase 15); the term is not adopted (user, 10-10 11:20; findings 4.14).

### Phase 14 - ff3.1: every basin depth of every map, under an L2 prior (STARTED 2026-10-09, user)
Route 1 of findings 1.28 (user, 10-09: no new simulations for input data). Each coil map (central,
direction, neighbour) gets one offset per basin of `rama_basin.py`'s six (alpha_R, alpha_L, beta,
pPII and the two phi > 0 extended regions; five free per map after renormalisation), trained on that
map's own reads with its own Gaussian prior centred on NDRD. No pooling across maps (rule 1.8);
the sheet group stays as in ff3.0. The update is the depth rounds' damped Newton step with the
prior, which is the MAP (ridge) estimate, not Adam (findings 1.28: Adam moves noise-only
coordinates at its full rate).
1. [x] **Probe on recorded data, no job** (findings 1.29): cross-validation finds a finite width
   (0.3 nat), but per-map depths transfer less well than pooled ones and what transfers is a global
   alpha_R deepening. Go or no-go is with the user. Method: ff21_ctrl_dt009 (c9: ff2.1's NDRD maps unchanged, 38
   steps, each protein twice) saved every residue's native-restrained and free basin populations in
   its `divergence.pkl`. For every (map, basin): the free-native gap over the reads the map serves,
   with each residue's left/right share; its protein-bootstrap error; split-half reliability (protein
   halves, and epoch 0 against epoch 1); the MAP step for a range of prior widths, with the width
   chosen by cross-validation over protein folds of the linear response. Report how many depths carry
   signal at that width, which, and whether the push points to the native placement of the map's
   class (1.15). Go or no-go to the user.
2. [x] Trainer (`training/ff31/`, findings 1.30): `rama_basin.py` generalised to per-(map, basin)
   offsets with a trainable mask; its reads reproduce round 4's round-1 glycine statistics exactly
   on recorded data, and its first step was simulated on c9's epoch 0 in both modes.
3. [x] Both variants (user, 10-09 23:05) from ff2.1 with every offset at 0 (user): `ff31_depth_all`
   (5,040 depths, prior 0.3) and `ff31_depth_data` (1,442, prior 0.2), midway2 49218340 and
   49218341, submitted 23:40, target 76, same 456 proteins and workflow, no gate. The push at the
   first epoch end (`rama_rounds.txt`, held-out gap left) is read before the user lets more epochs
   count.
4. [x] Judged by the selection panels and the steric check; **reviewed 10-13 (user; findings 1.34):
   ff31_depth_all is worse than ff_3.0_gdepth (folded -0.058, paired) and keeps losing; ff31_depth_data
   ties it (folded -0.012, not resolved; helix and beta lean better by about one SE). Neither earns
   replacing ff_3.0_gdepth on the panel; whether fd_03 is released for the Peng benchmark is the
   user's.** Both runs reached target 76 on 10-13 (ff31_depth_data 03:37, ff31_depth_all 03:53). Four depth
   rounds closed 1.6% (`all`) and 1.3% (`data`) of the held-out gap, `all`'s last round none
   (findings 1.30); folding fell over the first three epochs and held in the last, `all` 0.471, 0.444,
   0.404, 0.406 (fa_00-03) and `data` 0.479, 0.472, 0.449, 0.452 (fd_00-03), against ff2.1's 0.615 and
   gdepth_start's 0.603, which dominates every end in helix; at the target `data` leads `all`, 0.452
   against 0.406, helix -0.024 against -0.030 (findings 1.21); no steric flag except
   `data`'s step 39, shown to be scatter (findings 1.33). The choice is the user's. The local glpG TM4 glycine-flip test was
   retired by the user 10-11 23:10 (ff2.1 holds TM4 on the fixed table, findings 4.15; "we have a
   working ff and we are trying to expand it"), after fa_00_sc and fd_00_sc. Whether TM4's HDX dG
   should be effectively infinite is unresolved (user, 23:15); its REMD test was cancelled 10-14
   (Phase 16).
5. [x] **Steric watch** (user, 10-10 21:10, after Jumper: a moving map may open conformations where
   the side chain and backbone clash; findings 1.33). Upside has no excluded volume between a
   residue's CB and its neighbours' backbone, so the map is what keeps them apart.
   `training/ff31/steric_check.py` (deployed `checks/ff31_steric/`) measures the step's library
   against a hard-sphere Ramachandran map (outer-limit contacts, ideal geometry, CB, amide H, a
   following proline's CD; glycine without CB), as side-chain and backbone clash probability of
   every coil map and of every training residue's engine map against NDRD, and the free ensembles'
   non-glycine phi > 0 populations per step. The watch runs it at every ff3.1 step (cron
   `c618e779`) and notifies on a `STERIC FLAG`; the jobs change only with the user. Step 8: no flag.
6. [ ] **Validate ff31_depth_data's step 75 (fd_03) as `ff_3.1_depth_data`** (user, 10-13 12:00: "Use
   Peng validation set and glpG TM4 to validate ff3.1"; ff31_depth_all is not validated, findings
   1.34). Launched 12:10 10-13: the established `validate_ff.sh` on `epoch_03_minibatch_18` from a
   2-CPU launcher job, as the gates run it. Released to both trees, every file equal to the panel's
   `params/fd_03`; 32 Peng arms with ff_3.0_gdepth's targets; 4 glpG REMD chains whose seeds differ from
   ff_3.0_gdepth's only in the five force-field datasets, so they carry the fixed side-chain/lipid
   table (remote_jobs.md §1). Judged against ff_3.0_gdepth's arms and chains, TM4 above all. The name
   avoids the retired `ff_3.1_trained`. **TM4 stays helical** through block 1 and block 2's first
   chunks of all four chains (0.93-1.00 per group; user, 10-14 01:00; findings 1.35), as under
   ff_3.0_gdepth on the same chunks; GLY136 and GLY149 flip at some rungs under both, and only 79HIS's
   GLY136 under ff3.1 does so at T 0.70 (0.66 of frames). Peng arms: 14 of 32 complete by 10-15 13:50,
   all healthy; 18 run.

Risks recorded before the probe (findings 1.28): the depths can learn native placement and absorb
H-bond and side-chain physics, which a prior does not prevent; 46 held-out proteins sit at the
sample-size floor (1.12), so the width needs cross-validation over many folds.

### Phase 15 - the side-chain/lipid table at converged sampling, in production (DONE 2026-10-10, user)
The user's word (10-10 08:30): apply Phase 13's table fix to the force-field files so future
simulations run with it. The fix is a sampling density, not a parameter: the same potential, Boltzmann
average and grids, at 600 target directions x 12 side-chain bead frames instead of 24 x 1
(findings 4.10, 4.13). Not part of it: the membrane H-bond term (option A, not adopted; findings 4.14).
1. [x] Builder `py/martini_build_tables.py`: the per-sample loop computed as arrays, sampling
   600 x 12 by default (explicit arguments; the env overrides removed). At 24 x 1 it reproduces the
   old table (9 of 5.1 M cells one float32 step off); at 600 x 12 it equals the tested table to
   float32 rounding (findings 4.15).
2. [x] `parameters/ff_2.1/martini.h5` rebuilt (md5 `2e612261...`; old kept as
   `.bak_pre_sc600x12_20261010`): `particles` identical, `sc_table` converged. Every new hybrid
   preparation reads it.
3. [x] `up.md` 2.6 documents the sampling; findings 4.15 the build and its checks.
4. [x] Seed-patched runs (user 10-10 08:45: everywhere). `patch_glpg.py` (md5 `97ca3294...`) copies
   martini.h5's table into the seed. The 8 queued glpG chains (ff_3.0_gdepth, ff_3.0_bio_si) were
   re-patched and checked. TM4 sets run as `<tag>_sc` from now on.
5. [x] Cluster: martini.h5 in `$P` and `/beagle3/.../upside2-md`, the builder in `/beagle3` and
   `patch_glpg.py` in `$P/training/` and `tm4_local/scripts/`, md5 checked
   (`checks/sc_table_20261010/`). `$P/py/martini_build_tables.py` stays old (campaign rule).
6. [x] TM4 references on the fixed table, run locally 08:50-22:18 10-10: d9_03_sc, bs_03_sc,
   ff21_bioT1_6_sc, gdepth_start_sc, ff21_released_sc, d9_01_sc and b9_03_sc. d9_03_sc reproduces
   Phase 13's table-only set within noise; ff2.1 holds TM4 on the fixed table, every start reads
   alike, and none of the three released candidates is resolved from its start (findings 4.15).
   Each later epoch end is patched with the new tool and compared with its `_sc` start and twin.

### Phase 16 - glpG TM4 HDX on the fixed table (local readout done; ff2.1 REMD cancelled; ff3.0 REMD readout done 2026-10-14; ff3.1 REMD readout done 2026-10-17)

The only hybrid HDX set (ff2.1, 09-04) gives TM4 3.5-4.9 kcal/mol against the H145A peptide
140-144's 6.1-7.4, through protein protection failing (`pp_fail` 0.032 on 135-151 at the T 0.85
rung, findings 5.3). It predates the 09-10 temperature fix and the 10-10 table (findings 5.5). The
test repeats its force field and temperature with both fixes: if TM4's `pp_fail` falls, the gap
was the stale inputs; if not, the residual candidates of findings 5.5 stand. User, 10-11: "go
ahead and run this session's job first", ahead of the TM4 queue.

* Inputs: `ff3_local_test/patched/79HIS_ff21_released_sc.up` (ff2.1, fixed table, the 10-10
  `ff21_released_sc` set's input), 12 seeds, T 0.85 (298 K), 4000 tu, frame interval 2, the
  `run_glpg.sh` integrator. Runs and analysis in `scratchpad/tm4_hdx_local/`.
* Analysis: the repo's HDX steps unchanged (projection, `get_protection_state.py`,
  `martini_hdx_membrane_accessibility.py`) on the ff2.1 HDX topology; no MBAR, no uptake curves.
  Readout per TM4 amide 135-151: `pp_fail`, `acc`, exchanged fraction, dG at 298 K with the
  censored bound; 142-144 against the peptide. The 10-10 T 0.80 `ff21_released_sc` and
  `d9_03_sc` runs are read the same way at no cost.
* Limits: 4000 tu from the crystal is not the REMD's equilibrium, and 79HIS is not the
  experiment's 79ALA; the midway3 REMD analysis (remote_jobs.md §1) stays the matched readout.

- [x] 12 seeds run 07:09-09:04 10-11: KE/1.5kT 0.999-1.004; TM4 DSSP alpha 0.94-0.95 per block,
  0 unwound, 0 flipped (`results/tm4_local_T085.txt`).
- [x] HDX topology (the stale set's, md5 166c9d12) and per-seed arrays; tail-contact radius 7.00 A
  on every T 0.85 seed.
- [x] Readout (findings 5.5): ff2.1's TM4 137-148 lose protein protection in at most 0.06% of frames
  at 298 K, against 10-25% in the surviving pre-fix set; 142-144 never exchange (dG > 5.8, the
  frame-count bound). The pre-fix lability does not reproduce on current inputs. No REMD readout is
  pending: the ff2.1_sc chain below was cancelled by the user 10-14 01:11. HDX on the ff3.0 chains
  stays declined (user 10-13 10:20, "Do not do hdx on ff3.0"; a 20:25 "ff3.0 only" was withdrawn ~21:00, when the user
  told upside2-md-98 directly to keep this plan), so no 79ALA readout is planned. If it is ever
  approved: ff_3.0_gdepth block 2, 79ALA then 79HIS, skip S115T.
- [-] ff2.1 REMD on the fixed table, cancelled (user, 10-11 17:20: "just run REMD on remote to see if the
  problem still exist", in place of a local control while the Mac Studio is busy). The ff3.0 chains
  cannot make this test: the pre-fix lability was measured under ff2.1, 79HIS, and ff3.0 changes
  TM4's three glycines. One 79HIS chain from the pristine ff_2.1 seed with the fixed table, the
  ff3.0 chains' protocol, HDX after block 2 as in remote_jobs.md §1; readout against the pre-fix
  set's 0.851 rung (138-148 `pp_fail` 0.10-0.25, 142-144 dG 4.1-4.4) and the HXMS 6.1-7.4.
  Requested of upside2-md-98, which submits and records it in remote_jobs.md. **Cancelled by the user
  10-14 01:11** ("you can cancel ff2.1 job"): `glpG-RKRK-79HIS.ff_2.1_sc` ran from 17:41 10-11, block 1
  whole and block 2 to chunk step 6,930 of 10,494; its files are kept and read back, and its HDX
  (with the pre-HDX checks: KE/1.5kT, mixing, TM4 by half of block 2) was not run (remote_jobs.md §1
  "Finished runs"). Its seed differs from the pristine ff2.1 seed only in the four table arrays,
  and from the local T 0.85 input only in rama_pot (at most 2.6e-6), so the REMD and local runs
  share inputs. 28 replicas over T 0.70-0.90, not the pre-fix set's 48, so rungs are compared by
  temperature. One block resolves the pre-fix 4.1-4.4 (censored bound ~5.0 kcal/mol at one 298 K
  rung) but not the HXMS 6.1-7.4.
- [ ] ff3.0 REMD readout (user 10-14 01:19: "check if you have enough data from ff3.0 or ff3.1 ...
  Check if the TM4 still stretched down from infinity now"). ff_3.0_gdepth block 2 is whole in all
  four chains (17-18 chunks, ~5,100-5,400 frames per replica, bound ~5.1 kcal/mol at one 298 K
  rung); ff_3.1_depth_data has block 1 only until ~10-16. 79ALA then 79HIS, the spec of the ff2.1
  item above. It reverses the 10-13 decline; the user confirmed it to upside2-md-98 ("Yes, 79ALA
  then 79HIS"), which runs it (pre-checks: `checks/hdx_ff30_20261014/`). Resolved if 137-147 stay censored above ~5.1 at 298 K and 148-151 fray no more
  than locally (ff3.0 T 0.80: 148-151 at 1.6-2.4, mostly one seed); HXMS stays a sampling limit.
  **This phase, the request and the deferred fixes of findings 5.5 moved to upside2-md-98 10-14
  01:35 (user: "append jobs of this session to another session"); upside2-md-87 closed.**
  - [x] Pre-checks on block 2 (`prechecks.py`, `prechecks_20261014.txt`): KE/1.5kT 0.998-1.002 per
    slot; swap acceptance 0.36-0.53 per pair; no 0 -> 27 -> 0 round trip in either variant (walkers
    span 8-23 slots; net ~1.7 slots per chunk against ~21 accepted swaps per frame interval: each
    configuration stays near the temperature its energy fits). Each rung is still canonical, but its
    frames come from the few walkers near it, so the readout says so.
  - [x] One caslake job per variant (60668572, 60668573, 01:59-02:04 10-14; the user asked for midway2, where no
    environment has mdtraj, pymbar and matplotlib), `checks/hdx_ff30_20261014/hdx_block2.sbatch`: hdx_cluster.sbatch's
    pipeline with block 2 copied out of the live file first (`extract_block.py`: rotated groups
    only, every dataset compared; /input from the seed, pos = the block's first frame), concat
    discard 0, the HDX topology from ff_3.0_gdepth (rama from `common/rama.dat`, as 1.config.py
    always does), the repo's `martini_hdx_project.py` (coupled potential; the cluster copy re-scores
    a protein-only one), the stability dG on PS and on PS_protein, the reference-style plot, and
    `readout.py` (dG or "> bound" at 0.80/0.85/0.90 for TM4 130-151, pp_fail/acc/exchanged at the
    nearest replica, 135-151 by half at 0.85). Output `<chain>/hdx_sc600x12_ff_3.0_gdepth/`
    (~40 GB per variant). 79ALA groups 20-38, 79HIS 18-35.
  - [x] Readout (findings 5.5): TM4's core is censored at every rung, above 5.85 at 298 K (79ALA's 137 and
    139 at the bound), 142-144 against the pre-fix 4.1-4.4; only 130-134 and 149-151 resolve. 79HIS
    passes the closing criterion; 79ALA passes on 137-147, its 150 slightly more frayed than locally.
    Whether TM4's dG is finite, as HXMS says, needs ~3e5 independent frames per rung, and block 2's
    walkers do not cross the ladder.
  - [x] **The same readout on ff_3.1_depth_data** (user 10-14 02:10: "Do the analysis on ff3.1 as well"), 79ALA
    then 79HIS, block 2 = groups 19-37, caslake 60979543 and 60979545 00:19-00:25 10-17 (findings 5.5): the
    core 140-147 censored at 298 K in both, as under ff3.0; the 135-139 turn resolves at 4.0-5.6 and 79HIS's 148
    at 2.72, where ff3.0's are censored; 79ALA's 137-148 lose protein protection in up to 4.6% of frames at the
    0.853 rung (ff3.0 0.37%), unevenly between halves. 79HIS passes the closing criterion on 137-147; 79ALA fails
    it at 138-139. Figures `~/Downloads/glpG_TM4_HDX_ff3.1_depth_data_block2/` (Mac Studio). The recipe as run:
    remote_jobs.md §1 "TM4 HDX on ff_3.1_depth_data: the recipe": block 2's
    groups from the logs' run counts (output_previous_0 is the seed's; each block's first group is its
    calibration, 287 frames to t 18.0; `prechecks.py` asserts the log's run count matches); `prechecks.py
    ff_3.1_depth_data <V> FIRST LAST JOBID` on midway2; then `hdx_block2.sbatch` with `FF=ff_3.1_depth_data`
    on caslake (midway3 socket: a Duo push), the HDX topology from ff_3.1_depth_data; then TM4 at 298 K
    against ff_3.0_gdepth's, one figure per variant. Launcher and pre-check take the force field and
    block as arguments (the ff3.0 run's copies kept as `*.ff30_run_20261014`; the refactored pre-check
    reproduces its output exactly). The 10-14 session moved this phase's tools off its scratchpad onto
    the cluster (`plot_tm4.py` beside the launcher).

### Phase 17 - ff3.2: ff3.1's depth training started from ff3.0 (STOPPED 2026-10-18 in step 69; rerun in Phase 21)
The user's word (10-15): "same as ff3.1, except that it starts with ff3.0 instead of ff2.1". It is
ff31_depth_data's run with one input changed, the start: ff_3.0_gdepth's six released files in place
of ff2.1's. Decisions, each the plainest reading of that:
- **Mode `data`** (1,442 depths with >= 20 native sites, prior 0.2), the ff3.1 variant released as
  ff_3.1_depth_data; `all` lost to ff3.0 on the panel (findings 1.34).
- **The trainer is ff31_depth_data's, byte-identical** (ConDiv.py `df3e2497...`, rama_basin.py
  `235a1bf7...`). It takes its base library from `upside_input/rama.dat` and its Adam groups from
  `init_param/`, so the start is changed by those files alone. The offsets start at zero on ff3.0's
  library and each depth's prior is centred there, as ff3.1's were on its start (NDRD), so ff3.0's
  glycine depths are the prior's centre and are not pulled back to NDRD. The trainer's printed
  "NDRD" and the steric check's "against NDRD" mean the base library, here ff3.0's.
- Same 456 proteins, held-out set (seeded), dt 0.009, target 76, no gate; its own minibatch shuffle,
  as every run had. Epoch ends tagged `f32_EE` for the panel.
- ff_3.0_gdepth ended training at margin +0.082, below the watch's +0.10 hold line, so ff3.2 starts
  under it. Whether that line applies to a run that starts below it is the user's.
1. [x] Run directory `$P/training/ff32_depth_data` (10-15 15:57): ff31_depth_data's trainer, env.sh and
   pdb_list; `init_param` = ff_3.0_gdepth's five Adam-group files; `upside_input` = ff31's protein files
   (hard links) with ff_3.0_gdepth's `rama.dat`; every file cmp-checked (`checks/ff32_deploy_20261015/`).
2. [x] Initialised: 456 proteins, 1,442 trainable, 46 held out; held-out set and trainable mask equal
   to ff31_depth_data's; `rama_round_00.dat` array-identical to ff3.0's library; printed start equal
   to ff30_gdepth_dt009's step 75 (dhb -0.5513, sheet mean 0.2283, margin +0.082).
3. [x] Submitted 16:01 after `--test-only`: midway2 49406711 (remote_jobs.md §1); `ff32_depth_data:f32_`
   registered in `submit_new.sh` and `~/watch_probe/stock.sh`.
4. [x] Stopped by the user 08:33 10-18 in step 69 (steps 0-68 complete) for the angular-gradient defect; rerun on the fixed engine as
   `ff32_depth_data_angfix` (Phase 21). Its readouts to step 68 and panels f32_00-02 are the defect side.

### Phase 18 - the hybrid side-chain placement defect: fix and rerun glpG (STARTED 2026-10-17, user)
Findings 4.17: the hybrid preparation stores `placement_data` gathered per row but `layer_index` as the
absolute library row, so the engine places every side-chain bead at `rcf[li[li]]`. User 10-17: "fix the
bug and rerun the affected simulations"; cancel the 8 running chains; rerun all 10 (ff3.0 and ff3.1 in
the 4 variants, and 79HIS under ff_2.1_nobp and its matched control ff_2.1_ctrl38). Decisions:
- **The layout is upside_config's (and master's):** `placement_data` is the full library table (86 rows;
  the 1-body scalar likewise) and `layer_index` the absolute row. The preparation writes that; the
  seeds are repaired to it. The engine is correct and does not change.
- **Repair, do not rebuild:** only the two `placement_data` datasets change, so membrane, packing, ions,
  coordinates and every other dataset stay bit-identical (patch_glpg.py's reason). The repair recomputes
  the full tables from the build-time `ff_2.1/sidechain.h5` by the preparation's own rule and refuses
  unless the seed's gathered rows equal them exactly.
- **The base seeds are repaired too** (`popepopg_REMD_mdw2/seeds/glpG-RKRK-<V>.up`, backups
  `.bak_pre_placefix_20261017`), so `validate_ff.sh` patches correct seeds from now on. The per-force-field
  seeds are made the pipeline's way, `patch_glpg.py` on the repaired base, and must equal the repaired old
  per-force-field seeds.
- **The deterministic seed** (79HIS, ff_2.1_nobp) is `patch_glpg.py` on the repaired base, then one row per
  side-chain residue at the kept state: placement, 1-body scalar, pair `index/type/id`, both coverage
  nodes' side-chain rows, and `martini_sc_table_1body`'s rows with `row_rotamer_index` = the kept state's
  index in its type. The side-chain/lipid and backbone/lipid couplings stay on, unchanged in form.
- **The old chains are kept**, renamed `<V>.<ff>.placedefect_20261017` with their STOP files, so the reruns
  take the canonical names the HDX scripts expect. Their trajectories were deleted 10-18 (user; HDX results kept).
1. [x] Cancel the 8 chains (09:40 10-17, STOP files first; none resubmitted).
2. [x] Fix `py/martini_prepare_system_lib.py` (`load_sidechain_rotamer_payload`); the same hunk in the `$P`
   and `/beagle3` copies (`checks/placefix_20261017/patch_prep_lib.py`). The fixed function reproduces the
   repaired 79HIS and 79ALA seeds' side-chain datasets exactly (`test_prep_payload.py`, run locally: the
   preparation needs Python 3.10+).
3. [x] Seed tools in `checks/placefix_20261017/`: `fix_placement.py`, `single_state.py`, `verify_seed.py`,
   `build.sh` (`build.log`), `verify_all.sh` (`verify.log`). single_state compares the libraries' geometry at the
   float32 the seed stores: the trained files carry ff2.1's float64 geometry to 1 ulp (2.2e-16).
4. [x] Built: 4 base seeds, 8 ff3.0/ff3.1 seeds (each equal, dataset for dataset, to its old seed repaired),
   ctrl38 79HIS, nobp 79HIS (747 rows -> 168, 1-body 0 at every kept state).
5. [x] Verified (`verify.log`): only the intended datasets differ; placement at `rcf[li]` to 3.6e-6 A (old: up
   to 7.9 A); 1-body exact (old: up to 22 E_up off); side-chain/lipid rows active; nobp's 168 rows equal the
   multi-state seed's at the same (residue, state). Every multi-state force field keeps ff2.1's
   `rotamer_prob` exactly, so `patch_glpg.py` leaving the 1-body table is right for them.
6. [x] Old chains renamed `.placedefect_20261017`; 12 live seeds kept as `.bak_pre_placefix_20261017` and the
   14 new ones installed (md5 checked); 10 chains `--test-only`, then submitted 49411315-49411324 (all R
   at once).
7. [ ] Records: findings 4.17 (done), remote_jobs.md, WATCH_STATUS.md; the watch reads the new chains' first
   blocks (calibration, replicas 0, 13, 27), and the TM4 checks and HDX are redone on them. Block 1 closed in all
   eight 18:52-20:39 10-18 and its TM4 read is done (findings 1.35: TM4 holds, no glycine at phi > 0 at T 0.70);
   owed: TM4 at each later block end, HDX on block 2.

### Phase 19 - a second native-state panel from mdCATH (DONE 2026-10-17, user)
The Peng benchmark's 16 natives and the 44-domain panel (60 proteins) do not resolve ff_2.1_nobp from its control;
~75-100 native-start proteins are needed (nobp run's `setup/notes/findings.md`). User 10-17: "find something new
into the validation set"; time step 0.009, as trained; ff_2.1_nobp and ff_2.1_ctrl38. Decisions:
- **The 90 domains** of `ff3_selection/mdcath_candidates/` (README there): mdCATH whole-chain CATH 4.2.0 S20 domains,
  30 per length bin 50-100 / 101-150 / 151-250, helical / beta / mixed balanced; not in the training or panel PDB ids
  and under 30% identity to them; X-ray at 2.5 A or better with one protein entity; fewer than 4 cysteines; folded at
  320 K in all 5 all-atom replicas.
- **The all-atom reference is mdCATH's 320 K** (CHARMM22*, TIP3P; 5 replicas x ~500 frames at 1 ns), not Charron's
  amber99sb-ildn at 300 K: paired comparisons between Upside force fields hold; absolute folded fractions are not
  compared with the 44-domain panel's.
- **A separate panel** in `ff3_selection/{domains,aa,runs}_mdcath/`, so the 44-domain panel and its history are
  untouched. `panel.py` gains an optional time step on `run` (default 0.015, so every existing and automatic panel runs
  as before) and an `aa_mdcath` builder beside `aa`; a `panel_mdcath.sbatch` beside `panel.sbatch`.
- **Fetch only what is used:** each domain's starting structure and its 320 K coordinates read from the Hub by range
  requests, backbone kept; nothing else of the 0.2-2.5 GB files.
1. [x] Fetched (Mac, 10-17): `mdcath_fetch.py` (in `mdcath_candidates/`, log `fetch.log`), 90 domains, 980 MB kept
   (`pdb/`, `backbone/`); 320 K frames 420-500 per replica.
2. [x] `panel.py prepare` into `domains_mdcath/` (`prepare.log`): 89 domains. 1gqoB00 excluded (`excluded/README.txt`):
   mdCATH's stored starting structure has C-N 3.1 and 4.8 A at its C-terminus, its trajectory at most 1.49 A.
3. [x] `panel.py aa_mdcath` into `aa_mdcath/` (`aa.log`): every sequence equal; Q median per domain 0.891 (10th
   percentile 0.796, lowest 0.636) against Charron's 0.910 (0.824, 0.675); 17 of 445 replicas under Q 0.7.
4. [x] `panel.py` (backup `.bak_pre_mdcath_20261017`): `run` takes an optional time step (default 0.015, so every
   existing and automatic panel runs as before), `aa_mdcath` beside `aa`; `panel_mdcath.sbatch` beside `panel.sbatch`.
5. [x] `mdc_nobp` 49412017 and `mdc_ctrl38` 49412018 (released files copied, md5 equal) at 0.009, 356 runs each, into
   `runs_mdcath/`, COMPLETED 17:15 and 22:28 10-17; read by `panel_compare.sbatch` (49412845,
   `panel_compare_49412845.out`) and `panel.py select` (`select_mdcath_ctrl38_nobp_login1_20261017.txt`).

**Result (findings 1.21):** the deterministic side chain loses folding on held-out proteins. Over the 89 domains
no-BP folds 0.264 against 0.365 (-0.101 [-0.155, -0.051], Wilcoxon p 7e-5; lower in 59 domains, higher in 19), mean Q
-0.035 (p 2e-7), CA-RMSD +0.21 A (p 0.004), in the same direction in every length bin and fold class. ff_2.1 keeps the
probabilistic side chain. Open: whether the 44-domain panel's -0.016 differs from this because of its proteins or its
time step 0.015; not needed for the decision.

### Phase 20 - audit the no-BP training protocol, then close the branch (STARTED 2026-10-18, user)
User 10-18: no-BP performs below ff2.1; check that the training protocol was correct and that its training
has ended. Read only: nothing is resubmitted or cancelled without the user's word.
1. [x] The single-state library, configs and releases: correct (findings 1.21 "Audit").
2. [x] The engine and trainer at n_rot 1: correct, derivatives to float precision. Found the shared angular-gradient
   defect (findings 3.12), the same in both runs and in master; not fixed.
3. [x] Matching to ff21_ctrl_dt009: identical but for the side chain.
4. [x] The mdCATH verdict recomputes exactly.
5. [x] No no-BP training job is queued or can start (target 38 met; a resubmission exits). Still running from the
   comparison: glpG 79HIS ff_2.1_nobp 49411324 and ff_2.1_ctrl38 49411323, ctrl38's ubiquitin de novo 49413536.
6. [ ] **The gate fails no-BP's `rot` group at step 38** (findings 1.21 "Audit"): the verdict holds at matched
   length, not at no-BP's fixed point. Close as is, or train on to the gate: the user's call, with the three
   running jobs and whether 3.12 is fixed. Answered by Phase 21.

### Phase 21 - fix the angular-gradient defect and retrain what it touches (STARTED 2026-10-18, user)
User 10-18: "fix the issue and continue training if needed. also check if other on-going training also have the
defect, resolve them as well. re-run the validations if necessary." Then, on ff3.2: "In addition to this, check ff3.1
as well"; stop the step-37 glpG chains. Decisions:
- **The fix is in the engine** (`src/interaction_graph.h`). For a symmetric interaction, p[t2,t1] is p[t1,t2] with
  the two elements' roles exchanged (`is_compatible`), so the derivative added to [t2,t1] is the edge's parameter
  derivative evaluated at p[t2,t1] with the elements exchanged, in place of the [t1,t2] derivative copied unchanged.
  Fully symmetric parameter sets (radial splines; the pair term's radial entries) get bit-identical derivatives.
  The pair term's angular entries change: off-diagonal to the right values, diagonal redistributed between the two
  slots with the sum the trainer reads unchanged to rounding. Energies and forces do not use this path.
- **The trainer's uniform 2x on pair coefficients stays.** It sets the regulariser's relative weight and has since
  ff2.1's own training (master's `rotamer_parameter_estimation.py` builds the same layout).
- **Install only when no training is mid-step on the old library:** after ff3.2's last step, ~14:30 10-18.
- **No-BP and its control continue from step 38 on the fixed engine, each to its own fixed point:** at step 76
  the gate; not converged, one more epoch, up to 6 epochs (step 114), then stop for review. At the converged epoch
  end, the mdCATH panel (dt 0.009). No new release names in `parameters/`.
- **ff3.2 was the only training running: stopped in step 69 and rerun from step 0; ff3.1 (`ff31_depth_data`, the
  released one) rerun too.** Each rerun keeps its source's trainer, start, inputs, minibatch order and held-out set, so
  it differs from the source only by the fix; target 76, no gate, as the sources; epoch-end panels (`a31_EE`, `a32_EE`)
  and steric checks as theirs. ff3.0 and ff2.1 were trained with the defect and are not retrained unless the user asks;
  ff3.2's rerun starts from ff3.0 as ff3.2 did.
1. [x] Fix in the repo (`src/interaction_graph.h`, 32 lines); compiles locally.
2. [x] Built on midway2 from `$P/src` plus the fix (`checks/angfix_20261018/obj`, md5 libupside.so `3d9c6366...`, upside
   `7136d99a...`). Against the installed build (`parity.log`, `parity_glpg.log`): energy, forces and every other
   derivative ConDiv reads bitwise equal on single- and multi-state 1aba at two frames; radial pair slots bitwise;
   200 tu binary runs bitwise; hybrid glpG seeds' energy and forces bitwise. With the fix (`fd.log`): each angular
   entry against tied central differences ratio 1.000 (5.8e-6 median error single-state, 9.6e-4 multi-state); the
   trainer's angular-logit gradient 2.000x, the same as every radial coefficient (before: median 5.1 and 1.11).
3. [x] Checked in (10-18): ff3.2 stopped and rerun, ff3.1 rerun as well, the step-37 glpG chains stopped.
4. [x] Stopped 08:33 10-18: ff3.2 49406726 and successor 49410635 (steps 0-68 kept, step 69 partial), glpG 79HIS ff_2.1_nobp 49411324 and
   ff_2.1_ctrl38 49411323 (STOP files first). Installed 08:35 (backups `.bak_pre_angfix_20261018`); the installed
   files cmp-equal to the tested build.
5. [x] No-BP 49415222 and control 49415223, 38 -> 76; `after_training.sbatch` (no-BP's old one backed up) runs
   `ff21_nobp_dt009/setup/after_gate.sh`: the gate, then the mdCATH panel `mdc_<nobp|ctrl>_EE` or one more epoch.
6. [x] Reruns `ff31_depth_data_angfix` 49415231 and `ff32_depth_data_angfix` 49415233 (`checks/angfix_20261018/
   build_rerun.sh`; initial checkpoints equal to their sources' in every field but `native_sites`, summation order
   2e-15, which sets only the identical trainable mask); `a31_`/`a32_` in `submit_new.sh`, all four runs in `stock.sh`.
7. [ ] Validations as each run ends: the gates and mdCATH panels (no-BP, control); `a31_EE`/`a32_EE` panels against
   `fd_EE`/`f32_EE` on folded fraction and Q (the old panels' `.up` files were deleted 10-18, so CA-RMSD cannot be
   compared against them; the same holds for the step-38 `mdc_nobp`/`mdc_ctrl38`). Records done (findings 3.12, remote_jobs.md); the watch session told (WATCH_STATUS.md is its).

### Phase 22 - why de novo folding fails on large proteins (STARTED 2026-10-19, user)
User 10-19: Upside folds large proteins de novo poorly; investigate, read the literature and other force fields'
designs, propose, run any test needed. The user's hypothesis: folding is co-translational (S. Fried, JHU: about a
third of E. coli proteins do not refold in vitro, more often the multi-domain ones), so a thermodynamic force field
cannot fold them from a denatured chain. Facts so far: the 456 training proteins are 50-150 residues (max 150), the
44-domain panel 50-75, the Peng de novo set at most 92; on the mdCATH native-start panel ff_2.1_ctrl38 folds 0.42,
0.44 and 0.24 of frames in the 50-100, 101-150 and 151-250 bins.
**Design.** De novo failure is thermodynamic (the native basin is not Upside's free-energy minimum, so no route reaches
it) or kinetic (it is the minimum but bulk search, or a missing co-translational route, does not find it). Each test
separates the two; the co-translational test is run only on proteins whose native basin Upside holds.
1. [x] Size artifacts in the panel metrics: the 151-250 drop is the cutoff, not drift (findings 11e).
2. [x] Native basin stability against size: answered by item 1 and item 4's native starts; the 89-domain rerun dropped.
3. [x] Native energy gap: quenched native against the de novo globules, by term (findings 11e).
4. [x] 2LZM, 2RN2, 1YOB and the 1ENH control at T 0.80: native x3, bulk x6, vectorial x6. Compact non-native globules
   150-180 E_up above native; vectorial equals bulk (findings 11e).
5. [x] Literature and other force fields; code audit (findings 11e).
6. [x] Coil and native starts at T 0.88 and 0.95 (2LZM, 2RN2, 1ENH): no coil folds; RNase HI's native holds at 0.88,
   the control's melts (findings 11e).
7. [x] Proposal to the user (10-19 06:50).
8. [x] Fix `environment_coverage_hb`'s covering-residue types (user 10-19: "Fix the bug"). The writer stores one type
   per residue, `aa_types[r]` the type of residue r, as `environment_coverage_sc` already does; the engine sizes
   `aa_types` from its dataset and requires every `id2` to index it. A file written before the fix holds per-site
   types and runs bitwise as master; a newly written file gets the right types, so new configs depart from master
   by intent. Tests: old files bitwise before and after; the new output against a numpy coverage reference; finite
   differences; ff_2.1's energy change on natives and globules. Repo and local `obj/` only: the cluster trees are not
   touched while trainings run, and ff_2.1 and every training so far ran with the defect.
   Done 10-19 09:20: old files bitwise, numpy reference 5e-7, finite differences 1.0001, energy change at most
   1.3 E_up, dynamics normal (findings 3.13).
9. [x] `upside3.md`: the findings and the proposal; Upside 3 starts after the side-chain verdict (Phase 21 item 7).

**Review.** For these in-vitro refolders Upside's de novo failure is search into compact non-native globules that
it ranks well above the native; N-to-C growth does not help. Proposed, for the user's decision: compact misfolds as
negative examples in ConDiv, training past 150 residues, and replica exchange for large-protein benchmarks. Not
measured: the globules' free energy against the native's, and a size contrast under one protocol (the 1ENH control
does not fold at constant T in 60,000 tu); both need replica exchange.

## Known Errors / Blockers

* **Every ConDiv run before 2026-10-18 08:35 trained the side-chain pair term's off-diagonal angular coefficients on
  a mixed gradient** (findings 3.12; fixed in the engine then): ff2.1's retrainings to step 37, ff3.0, ff3.1, ff3.2.
  The fix's effect is measured by the reruns (Phase 21); ff3.0 and ff2.1 are not retrained.

* **Every hybrid glpG result before 2026-10-17 ran with misplaced side-chain beads** (findings 4.17; Phase 18):
  TM1/TM4 helicity, glycine flips and HDX, under every force field, on the cluster and locally.

* **Every glpG TM4 count before 2026-10-07 11:22 is invalid** (findings 3.11). That covers:
  - every local 3- and 12-seed set (Phase 10, Phase 11);
  - the 09-30 ff_3.0 glpG validation (Phases 3-5);
  - any other glpG run from the live seeds under a force field other than the FF1-form ff_3.0.

  Each of these used the FF1-form ff_3.0's coverage tables. `patch_glpg.py` is fixed (local, `$P/training/`,
  `tm4_local/scripts/`; md5 `97ca3294...` since 10-10, which also carries the side-chain/lipid table);
  the automatic validation uses it. Another computer must update its local copy before use
  (remote_jobs.md "Setting up the next computer").
* **Every hybrid run before 2026-10-10 08:39 used the under-sampled side-chain/lipid table**
  (findings 4.10, 4.15): every round-4 TM4 set and both candidates' first glpG patches. Their
  comparisons with each other stand; their absolute TM1 and TM4 values do not, and they are not
  compared with `_sc` sets.
* **Total-potential jumps in local glpG runs** (findings 1.22): 13 of 105 pre-fix seeds, one blow-up.
  The likely cause is the old table's kicks (findings 4.10; none in the 24 fixed-table seeds of 4.13
  and 4.14), unproven for the jumps themselves. Check the `_sc` sets for them; do not drop frames. One-frame backbone
  excursions below the 3000 cut sit where the helices fail (findings 4.9, Phase 13's second lead).
* **The TM1 term is not adopted** (user, 10-10 11:20; findings 4.14). Its f was trained in Upside's
  implicit membrane and it is anchored at z = 0; any later test of it watches the protein's depth
  against the explicit bilayer (`analyze_hbmem.py` health block).
* **midway2-login2's `/project` hangs** (since 10-07 10:48). Route `/project` work through login1
  (remote_jobs.md §8).
* **The cluster's `training/` keeps the pre-10-02 layout** while round 4 trains. Do not sync the
  repo's `training/` over it.
* **`run_remd.py` rolls back NaN replicas and continues**, which the NO GUARDS rule forbids. It is
  kept unchanged so glpG stays comparable across campaigns. Whether to remove it is the user's call.
* **CPU work runs on midway2 broadwl** unless the user names midway3 caslake. No `amd`, `beagle3`
  or GPU partitions, and not `broadwl-lc` (its nodes cannot see `/project`).
* **HDX post-processing items await the user's decision** (user, 10-11; findings 5.5): the
  `4.calc_D_uptake.py` PF = 1000 cap on never-exchanging amides (a master line), three docstring/README
  statements that call TM helices non-exchanging, and the unmeasured interface boundary of the 7 A
  tail-contact criterion. Nothing is changed until then; do not compare uptake with the HXMS data before.
* **lambda's helix 2 / helix 4 misorientation** is a known ff2.1 limitation with no single
  side-chain cause (findings 11d). It is re-tested after ff3.0, not trained against.
