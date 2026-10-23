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
   eight 18:52-20:39 10-18 and its TM4 read is done (findings 1.35: TM4 holds, no glycine at phi > 0 at T 0.70).
   HDX on ff3.0's blocks 2 and 3, all four variants, done 10-22 (findings 5.6: the blocks disagree beyond their
   error bars, so a longer run is the user's decision); owed: TM4 at each later block end, ff3.1's HDX.

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
   No-BP and the control both converged at step 76 (10-19, 10-20); `mdc_nobp_03` folds 0.241, `mdc_ctrl_03` 0.353,
   paired -0.111, p 1e-5 (findings 1.21; Phase 25). Owed: the `a31_EE`/`a32_EE` panels (none submitted by 10-20).

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

### Phase 23 - plan Upside 3 (STARTED 2026-10-19, user)
User 10-19: Upside 3 starts after the side-chain verdict, with five goals:
* a debug of the whole repo, all bugs found by us as Jumper did on master (this branch, master's code first);
* unit alignment;
* new parameters chosen after the verdict;
* large proteins (Phase 22);
* GPU support with one force field for both backends.

The plan, the order and the decisions waiting on the user are in `upside3.md`. The research reports are in
`scratchpad/upside3_planning_20261019/`.
1. [x] Defect inventory, audit units and method; first-pass warnings, sanitizers and probes (`upside3.md` 2).
2. [x] Unit map and options (`upside3.md` 3). E_up's kJ/mol value is FF2's temperature calibration (298 K / 0.85),
   so A+ is recommended: reduced units inside the engine, standard units at every input and output.
3. [x] GPU profile, routes, the CPU-GPU equivalence ladder and hardware (`upside3.md` 6).
4. [ ] The user's decisions (`upside3.md` 7), then the audit.
   * Decided 10-19: Upside 3 is a new release, so no list of fixes and their effects on ff2.1 is published. Each fix
     is proven correct, but none is measured against ff2.1 or master. ff2.1 is not re-validated on the audited engine.

### Phase 24 - Upside 3 audit: find and fix every bug (STARTED 2026-10-19, user)
User 10-19: "the fix only goes into upside 3. You can modify on current branch (martini-dev) as this will be the
'upside 3'. Later in the stage we should make it refuse to run ff2.1 (not yet ...). You can start the debug progress
now, on this branch." Code concise, simple and logical, in master's style (never example 00's).

**Resume here ("continue upside3 progress"; state at 2026-10-22 19:10, Mac with the scratchpad).** Working files
(fixtures, baselines, build trees, reviewer outputs) are under the gitignored `scratchpad/upside3_audit/`.
Running on midway2 (item 2b; ids, logs and next actions in remote_jobs.md §1 line 3): the Peng arms and the glpG chain.
1. [x] `obj/` holds every fix to A231 (earlier binaries `obj/*.bak_a230`, `*.bak_a225`, `*.bak_a220`, `*.bak_a211`, `*.bak_a204`,
   `*.bak_a196`, `*.bak_a182`, `*.bak_a174`, `*.bak_a163`, `*.bak_a158`, `*.bak_a155`, `*.bak_a149`, `*.bak_a148`,
   `*.bak_a142`, `*.bak_a138`, `*.bak_a125`, `*.bak_a123`, `*.bak_a121`, `*.bak_a110`, `*.bak_a101`, `*.bak_20261021`).
   Baselines `w56_obj.npz` on `fixtures_w21` (the writers' current set: fixtures_w20 plus A231's `sequence_id`; w55 on
   fixtures_w20 was fixtures_w19 less A227's dead items, bitwise w54; older sets are refused in part: pre-A130 anchors, pre-A172 channel configs, pre-A222 pulling configs)
   and `h53_obj.npz` on `fixtures_hybrid_w6` (fixtures_hybrid_w4 less A227's 13 items plus A231's 4 datasets; h53 is h52 moved at float32
   rounding by A226's fast-math code generation. The w4 hybrids were regenerated 10-22 through the full prep with
   A189-A192: minimisers run their budgets; KE/1.5kT 6.6 1.008, burn-in 1.000, production 1.007; peptide C-N mean
   1.319, max 1.82 A; CA Rg 12.3-13.2 A; PO4 hold end 102150).
2. [x] Pass 2: every unit has both reviews (`pass2/<unit>/report.md`, all 18 P1 and 7 P2 units). Fixed from pass 2:
   placement and writers A108-A121, coordinate nodes A127-A132, engine graph A133-A138, dry-MARTINI engine (P2-01, -03,
   -04) A139-A142, driver/REMD/thermostat (P1-10) A143-A149, local backbone (P1-08) A150-A155 (D3 decided for the user,
   item 5), Monte Carlo (P1-11) A156-A158, hybrid coupling (P2-02) A159-A163, membrane and surface (P1-09) A164-A173
   (its M2, M3 left for the membrane retraining: findings A173), shipped parameters (P1-18) A174 (its D1 and D7 for
   the user, item 5), I/O and C API (P1-12) A175-A182, hybrid preparation (P2-05) A183-A188 (its D1, D6, D8 and the
   membrane-thickness definition decided, item 6), ConDiv trainer (P2-07) A193-A196 (D1 decided by the user: the SI's
   SARW), parameter glue (P1-14) A197-A199 and A204 (D1 and D6 decided for the user 10-22: torch replaced by
   closed-form numpy gradients, so a process holds one OpenMP runtime; the 86 dead `rotscalar` coordinates removed,
   which stops older checkpoints resuming in a newer tree, remote_jobs.md 0b), hybrid analysis (P2-06) A200-A203,
   build (P1-17) A205-A213, analysis (P1-15) A214-A220, examples (P1-16) A221-A225.
2a. [x] The closing criteria (upside3.md 2.2 item 5) on the final engine, met 10-22 (`scratchpad/upside3_audit/closing/`):
   * invariance: every test passes (equipartition on chig and on the hybrid production state, pivot MC, REMD,
     thermostat v / mv fast / mv slow; `closing/inv.summary`).
   * node_fd: all 30 configs; every parameter derivative exact or "not implemented" (A29-A31); every flagged
     direction, all in the rotamer node, line-integrates to noise (`closing/{fd,li,li2}.summary`, `gfix/li_75.log`),
     the pair-table ones once `tests/line_integral.py` follows node_fd's tied direction (A230).
   * sanitizers: an ASan + UBSan build, `tests/regress.py dump` on all 30 configs through the library and the binary:
     every config exit 0, no report, 4260 arrays (`closing/san.summary`; the preload goes to Python.app, whose
     `python3.11` launcher drops DYLD_*).
   * registry (`registry_20261022.md`, on fixtures_w20 and fixtures_hybrid_w5): 44 of 50 node types matched; the 6
     written-but-never-read rows hold only items read outside the engine's construction, each named in
     `pass2/registry/report.md`; its fixes A226-A229.
   ff2.1 stays usable (user 10-22: "keep ff2.1 useable for now"); refusing it waits.
2b. [ ] Validation of the debug on midway2 (user 10-22: "fix D9 with (a). submit Peng 16 validation protein with ff2.1
   and glpG hybrid simulation with ff2.1, on midway2, as a validation for our entire debug"):
   1. [x] D9 by design (a), A231: a required per-element `sequence_id`, the exclusion key of every node that
      excludes by sequence distance, written by every writer. `obj/` to A231 (`*.bak_a230`); baselines w56 on
      fixtures_w21, h53 on fixtures_hybrid_w6 (bitwise); registry 44 of 50 matched as before.
   2. [x] An Upside 3 tree on the cluster, `/beagle3/trsosnic/yinhan/upside3` (compute nodes of both clusters see
      /beagle3). It is separate from the shared deployment because running campaigns read that tree for a week (the
      released-ff2.1 Peng arms are this validation's control) and CLAUDE.md forbids deploying under them. Code
      rsynced from the Mac working tree; built on the midway2 login node (GCC 10.1, the cluster CMake file of
      A205-A206 checked there); its own `.venv` from the shared `pyrt/` with mdtraj, which the hybrid prep's DSSP
      needs and the shared venv lacks; `env_upside3.sh`.
   3. [ ] (submitted 22:26 10-22, ~a week) Peng: the 16 proteins x {native, denovo}, the benchmark's own protocol (`ff3_benchmark/bench_run.py`:
      14-replica REMD on Table S2's ladders, dt 0.009, frame 100), ff_2.1, in
      `/beagle3/trsosnic/yinhan/upside3_validation/peng/`; compared with the released-ff2.1 arms on the pre-audit
      binary (`ff3_benchmark/runs/<prot>_<kind>_ff_2.1/`) and Peng's FF2 (Table S2, Fig S4).
   4. [ ] (submitted 22:58 10-22 through `example/16.MARTINI/run.py`) glpG-RKRK-79HIS hybrid with ff2.1: the repo's job
      script, rewritten for the Upside 3 core (A232; user 10-22: "the job scripts should also be updated to support
      it"): the preparation as a job, then 4 REMD blocks of 48 replicas on T 0.70-0.90, each stopped at
      `--time-limit` and queueing the next. `run_sim_1afo_full.sh` sources the Upside 3 tree too (A233). The campaign's
      `run_remd.py` stays as it is under the running campaigns (Upside 2; it also carries a NaN-rollback guard).
   Next after it: refusing ff2.1 (the step `upside3.md` names after the audit), then the open items below.
   Open from the build review (`pass2/build/report.md`): the cluster binaries set flush-to-zero and
   denormals-are-zero in every process that loads the library (U1, `set_fast_math` found by `nm` in both trees),
   measured at the first midway2 build (remote_jobs.md 0b); `sse2neon.h` is an unpinned Homebrew header whose
   defaults set the Mac's `_mm_rcp_ps` accuracy (U2). Its B9 decided for the user 10-22: `.devcontainer/` deleted
   (A213); a container for Upside 3 is written when there is a runtime to test it on.
   Open, performance only: the g-JF round evaluates the forces once more than it needs (MARTINI review U4, about
   20% of a hybrid's cost); skipping it needs to know that nothing moved the coordinates since the last evaluation
   (barostat, rigid fit, MC, exchange), so it waits until correctness work is done. U6, a hot rigid protein in 6.6,
   is a transient.
3. [x] The glpG protocol validation (findings 8b): the rigid-stage gate failed one healthy prep in three, and every
   prep's flexible production closes the belt ends, so the gate measures production's steps 15,000-20,000 (A126,
   user 10-21: option (a)); `run.py` starts REMD from the gated production's last frame. Not yet rerun end to end
   through the workflow (the function is tested on recorded productions, `glpg_gate_20261021/gate_a/test.log`).
4. [ ] Open fixes. Done 10-21 and in the ledger: D5 as A103 (installed in `obj/` 16:18 with A104-A110; A101's
   binaries kept as `obj/*.bak_a101`), the four decisions below as A104-A107, placement D1-D3 as A108-A110, writers
   items 1-5, 7-12 and 15 as A111-A121 (6 and 13 are A106-A107), the `HBondEnergy` input checks as A122, BP D4 as A123, the dead-code batch as A124, the reader preconditions as
   A125 (all in `obj/`). The 1rkl hybrid fixtures are regenerated (`scratchpad/upside3_audit/
   fixtures_hybrid_w2/1rkl.stage_7.{0,1}.up`); `tests/regress.py` now runs a hybrid's seeded REMD without pivot moves
   (the engine refuses pivot moves with the PO4 hold's fixed atoms). Fixture baseline `fixtures_w14`. Left:
   * D9 (design in `pass2/igraph/report.md`): the reviewer's sketch reads an optional `sequence_id`, falling back to
     the ids, which the clean-slate rule discourages; a required dataset means regenerating every config. Only
     chain ends in contact are affected. Kept for a session with the user.
5. Decided 10-21 (user: "make decisions for me, do whatever test you need to"):
   1. `membranelateral_potential` is deleted (node, `--membrane-lateral-potential`, the `run_upside` passthrough,
      `main.cpp`'s z-centering list, index oracle entry, up.md, example 08's `2.lateral_pressure.run.py`,
      `lateral.dat`, its readme lines). It has no energy to report: its forces do not sum to zero, which any
      energy of the xy-uniform membrane must (`scratchpad/upside3_audit/lateral/net_force.py`, `.log`: 1qhj net
      force up to 0.096 E_up/A of 1.4-2.8 total, z torque up to 1.27 E_up, closed rigid loop work up to 0.11 E_up;
      1rkl 0.017, 0.15, 0.03). Each side's push is spread over its CBs at their own heights, so per z-slice the
      sides do not cancel, where a pressure P(z) on a closed outline gives zero net force. So it is not physical,
      and logged energy, pivot MC and replica exchange cannot see it. Only example 08 uses it. A lateral-pressure
      term, if wanted, is a new energy, the integral of P(z) A(z) dz over the protein's cross-section (upside3.md
      proposal).
   2. The rotamer-library generators (`compute_rotamer_centers.cpp`, `compute_rotamer_pos.cpp`,
      `generate_from_rotamer.{h,cpp}`, their CMake targets) are deleted: nothing builds or calls them, no Upside 3
      plan rebuilds libraries from Dunbrack, and master keeps them as the libraries' provenance (noted in up.md).
   3. Chain breaks: writers patch 6 applied as written. Every shipped workflow (tests/fixtures.py, examples 03-16)
      converts with `--record-chain-breaks` just before configuring, so only a reused single-chain input directory
      is affected, and it now fails loudly; re-run the converter.
   4. Insertion codes are refused by `PDB_to_initial_structure.py` with a request to renumber: every consumer of
      the PDB numbering (HDX peptide ranges, hybrid `bb_resseq`, plot labels, `martini_hdx_project`) is
      integer-only, and none of the project's 46 PDBs has one.
   Still with the user: `upside3.md` 7.
   6. Decided by the user 10-22 09:05 (all four as recommended):
      * the side-chain/environment table is rebuilt: each rotamer state's MARTINI beads from ideal-geometry heavy
        atoms at its library (chi1, chi2) rows, chi3/chi4 sampled, mapped by martinize's martini22 centres, the table
        averaged over those conformers (parameter review D1, prep review D1, confirmed against 1146 crystal side
        chains: the old table put the beads on Upside's trained rotamer point and vector);
      * its azimuthal average at the simulation temperature (0.8647 T_up), not 25 T_up (prep review U1);
      * dead content deleted: `parameters/ff_2.0/`, `parameters/packing/`, `dryMARTINI/DDM.pdb`, and the unread
        martini.h5 datasets (`rotamer_probability_fixed`, `type_charge`, `type_order`, `@coulomb_k_eup`) with the
        prep's requirement of the first;
      * the ion count uses the P-P thickness, the prep's fallback measuring it from the PO4 beads.
      With them (user 10-22 09:20: "Use DSSP"), the backbone bead types by DSSP through martinize's
      ssClassification (prep D8); it changes helical proteins too (glpG 84 of 210 residues: the old detector merged
      helices across kinks). Done: A189 (dead content), A190 (P-P ion count), A191 (DSSP types), A192 (the table rebuilt
      and installed 10:45, `martini.h5` md5 7496dd90; old one `martini.h5.bak_pre_martinize_beads_20261022`).
      The hybrid fixtures re-prepared end to end on them (`fixtures_hybrid_w4/`, 10:50-11:55): healthy, the
      hybrid baseline h52 (item 1). Then `martini.h5` is rebuilt (old one
      backed up, overwritten), the hybrid fixtures re-prepared end to end, and the cluster's live hybrid files
      re-prepared before Upside 3 is deployed there (remote_jobs.md 0b).
   5. `--secstr-bias` (local backbone pass 2 D3; user 10-22: "make decision for me"): its windows were non-periodic
      in phi, stepping by the full bias across +-180. Decided: one phi window for both basins, phi < 0 continued past
      +-180 to +120 (the empty gap between beta and alpha_L, so the beta region, which crosses +-180, is biased
      whole), and the psi windows on the same form; each window is the original pair of 10 degree edge sigmoids
      summed over its periodic images, so it is periodic and smooth. A cosine-form window was rejected: on a 120
      degree gap with 10 degree edges it leaves 3.5% of the bias on alpha_L (A155).
6. Git is the user's. The commit should include `tests/` and `example/15.SpatialTransformation/calc_rmsd.py`.

**Revised decisions** (CLAUDE.md Development Rules reworded to match, 10-19):
* martini-dev is Upside 3. Fixes go only here; master is not changed.
* ff2.1 keeps running on this branch for development. Refusing it comes later.
* The harness lives in `tests/` (new). Its generated fixtures and outputs go to the gitignored
  `scratchpad/upside3_audit/`.

Method in `upside3.md` 2.2.
1. [x] Local `obj/` rebuilt from the branch (old binaries `obj/*.bak_20261002`).
2. [x] `tests/`:
   * `fixtures.py`: ff_2.1 configs for chig, 1dfn, 1ubq and membrane 1qhj, plus one config per writer option group
     (10-21): every node writer of `upside_config`, `advanced_config` and `upside_nodes` that a shipped library can
     build, 44 of 44, in 28 configs (`fixtures_w6`; report `scratchpad/upside3_audit/option_fixtures_20261021.md`);
   * `node_fd.py`: each potential leaf alone; forces and every upstream parameter derivative against
     fourth-order differences on the realized float32 steps; tie maps for the pair table and sigmoid weights;
   * `regress.py`: node outputs, forces and parameter derivatives on 3 frames, plus a seeded 2-replica
     REMD-and-pivot run; bitwise compare.
   * `line_integral.py` (10-21): the line integral of a flagged parameter derivative at node_fd's own frames.
   * `invariance.py` (10-21): pivot MC detailed balance (8 residues, T 0.7 and 1.3, smallest p 0.06; against the
     wrong temperature X^2 15,700-56,000) and REMD on a harmonic ladder (per-replica moments within 2.5 SE, acceptance
     0.5466/0.4374/0.3641 against 0.5449/0.4354/0.3619) pass. Equipartition: chig reads +0.4%, the leapfrog
     half-step momentum's (dt^2/4)<F^2>/(m kT), which the test expects from the frames' forces; the SE now comes
     from the pooled statistical inefficiency g (pymbar), since lipid kinetic energy stays correlated for about
     55 frames under the diffusion-calibrated friction (findings 10.33). With the current build: chig N, CA, C
     within 2.8 SE (g 2.5); the 1rkl production record and the 7.1 continuation (4 x 45,000 rounds) PROTEIN
     0.9970 +- 0.0017, LIPID 0.9999 +- 0.0020 (g 68), ION 1.0097 +- 0.0090, all PASS
     (`invariance/equi_g.log`). The earlier PROTEIN 0.995 was seed 2's noise draw sampled twice (same-seed runs
     of two builds correlate at r 0.92-0.96 frame by frame; A82 is not involved, `ab82.py`, `pair_seed.log`);
     8 x 90,000 rounds on seed 3 give 1.0005 +- 0.0006. Found on the way: A85.
     Hybrid fixtures `scratchpad/upside3_audit/fixtures_hybrid/1rkl.stage_7.{0,1}.up`.
   * Sanitizers (10-21): an ASan + UBSan build (`scratchpad/upside3_audit/build_asan`) driven by
     `scratchpad/upside3_audit/sanitize_20261021/run.sh`: every node's static evaluation on `fixtures_w4` and both
     A83 hybrids, the seeded REMD-and-pivot run on each fixture, 200 production steps of each hybrid. Clean.
   * `index_oracle.py` (10-21): an independent numpy reconstruction of every writer's index arrays, compared slot
     by slot with the meaning the C++ reader gives them; 69 of 69 injected single-index errors caught. Clean on
     every fixture except D9's exclusions; found A94; up.md corrected in 16 rows (report
     `scratchpad/upside3_audit/index_oracle_20261021.md`). The hybrid fixtures predate A83 (rows for residues 7,
     11, 23 missing): regenerate them through the current prep before they serve as a baseline [ ].
   * `registry.py` (10-21): each registered node's C++ reads against its Python writers, statically and from HDF5
     traces of the engine; 69 node types, 40 matched, the rest benign or fixed (A83, A84; report
     `scratchpad/upside3_audit/registry_20261021.md`).

   `get_param_size` was added to the C API and `upside_engine.py`.
3. [x] Regression baseline: the first baseline was not reproducible (2 of 6 dumps differed). That was ledger A2's
   heap overflow. After A1 and A2 the record reproduces; the baseline is `scratchpad/upside3_audit/fix1_1.npz`.
4. [ ] Review pass 1: one reviewer per P1 unit. Each claim comes with a reproducer.
   * REMD reviewer running (user 10-19: Jumper's Claude found a REMD issue).
   * Harness candidates, resolved 10-21 (`scratchpad/upside3_audit/fd_candidates_20261021.md`): the rotamer force
     mismatches are float32 energy noise (about 60 ulp, node_fd assumes 30); the placement and coverage parameter
     mismatches are stencils crossing a kink (a clamped cosine, a compact-support edge). One real defect behind
     them, A82. node_fd now measures the energy noise along each direction (fourth differences at h/8, at least one
     ulp) and adds a third step h/4; on the four `fixtures_w4` configs it flags 3 placement directions, each a kink
     closer than h/4, and `tests/line_integral.py` settles them (residual at noise, 5-12 times below the claimed
     error; `fd_w5.txt`).
     obj/ rebuilt 10-21 07:23 from src with A82-A93 (old binaries `obj/*.bak_20261020`); its regression equals
     the `o12_home` baseline w30 bitwise (w31_obj, 0 of 4073).
   Reported: REMD, environment, interaction graph, BP, H-bond, writers, placement, membrane, hybrid coupling,
   hybrid preparation, dry-MARTINI engine, examples, hybrid analysis: every pass-1 unit.
   Fixed and verified: ledger A1-A101 (findings 3.14), plus the Python option batch (13 options repaired). The user's 10-19 approvals
   are applied: exact derivatives or "not implemented" for every parameter, dead nodes and options deleted, and the
   angular cosine held to [-1,1].

   Remaining queue:
   * applied 10-20: dry-MARTINI engine side and the `martini.h5` rebuild (A74-A76), examples (A78), the swap-set
     restart (A77).
   * the box test (fresh 1rkl prep, every patch, new table, seed 11; `scratchpad/.../fullstack_1rkl/`): exit 0;
     KE/1.5kT 0.55-0.82 over 6.2-6.6 (recovering from minimisation), 0.982 burn-in, 0.991 production; peptide C-N
     mean 1.334 A (sd 0.134 against 0.116-0.120 in pure-protein runs); lipid bond sd per bond 0.440 A against the
     harmonic 0.449 (the old 8.5 A stretches gone with A75). The protein holds 2.2 H-bonds while rigid, as pure
     ff2.1 scores the 1rkl NMR backbone, and relaxes to helix within 50 steps of release. Box expansion (100.2 to 110.5 A over 6.2-6.6):
     explained, no engine defect. GROMACS running CHARMM-GUI's own protocol on the same step5 bilayers does the same
     thing: DOPC area 69.7 to 88 A^2 within 50 ps, POPE/POPG 60.4 to 74, then recovery over about 12 ns (2 ns without
     the PO4 restraint) as the overlapped leaflets separate. Upside's prep covers about 5 ps (2500 steps of 1.8 fs),
     so it only sees the jump. GROMACS equilibria at the published 1.2 nm (second half of 50-100 ns; full table in
     `scratchpad/bilayer_ref/results_20261020.txt`): DOPC 64.8 A^2 / P-P 43.8 A with 0.15 M salt, 62.8 / 44.9
     without; POPE/POPG 61.8 / 41.8 with or without salt; CHARMM-GUI's 1.1 nm gives +1.9 to +2.4 A^2 and the
     Arnarez force switch +3 A^2. This replaces findings 8's 65.6 A^2 / 43.3 A (old barostat).
     Upside pure bilayers from leaflet-shifted starts (0.55 ns each) have not converged: DOPC 59.6-62.5, POPE/POPG
     57.6-60.4 A^2 depending on the starting area, thicker and less overlapped than the GROMACS equilibrium.
     Decision: the prep starts from a GROMACS-equilibrated membrane of the same potential (no geometric correction;
     the leaflet-translation patch is dropped).
     Validated 10-21 (`scratchpad/bilayer_ref/upside_from_gmx_eq/results.txt`, findings 8): started at the GROMACS
     equilibrium, Upside holds it for DOPC (63.3 against 62.9 A^2, P-P 44.7 against 44.9); POPE/POPG ends 1.3%
     smaller and 0.4 A thicker, still drifting at 0.73 ns. 1rkl from the equilibrated 2 x 2 DOPC template keeps its
     membrane through prep and production (tail ends past the midplane 0.26-0.28, P-P 41.8-44.7, C-N max 1.89 A,
     KE/1.5kT 1.001); from step5 it enters production interdigitated (0.77 past, P-P 31.7, box 110.5 A). glpG fails
     the stage-7.0 solvation gate from either membrane (111 of 136 and 119 of 138 belt sites bare, nearest
     environment bead 8.2 A on average), while the equilibrated membrane itself stays at P-P 41, 0.30 past: the gap
     is at the protein-lipid interface, cause open.
     glpG's bare belt, 10-21 (`scratchpad/glpg_belt_20261021/`): the gate measures from backbone sites, so the side
     chains between them and the lipids count as gap (on 1rkl's lone helix in production 40% of sites read bare at
     6 A); and the void the insertion leaves heals slowly: at fixed box, production closes the belt's nearest
     bead from 8.23 to 7.14 A over 50,000 steps, still falling, protein intact (C-N max 1.95 A). The prep's NPT
     totals 4.6 ps (5 x 500 steps) against CHARMM-GUI's 10 ns (5 x 2 ns, tau_p 4 ps), so the box cannot relax.
     Measured 10-21 (findings 8b): the 6 A rule fails every healthy state and counts glpG's bundle interior as
     vacuum; on lipid-facing sites (burial 40-60) healthy states have 1-4% beyond 8 A, the rejected 7.0 state 18%,
     the 6.6 NPT plateau 5% (box converged by 30,000 steps, lipid face by 60,000), production at 200,000 steps 4%.
     Applied 10-21 (user: "proceed with your suggestions"; ledger A102): the gate takes the lipid-facing belt sites
     over the rigid burn-in's last quarter and fails above 10% beyond 8 A (single frames are too noisy, about 30
     sites); 6.6 runs 60,000 steps. On the recorded states it fails the rejected one (17.9%) and passes the relaxed
     glpG states (3.6-5.5%) and 1rkl (0.4%). The first full prep with it passed at 9.9%: not yet validated (Resume
     here, item 3).
     [ ] a longer POPE/POPG run to settle the 1.3%; [ ] equilibrated membranes as the prep input.
   * applied 10-20: the MARTINI angle term as a spline table, and the side-chain table's CB direction normalised in
     the energy as in its gradient (node_fd now passes); regression w22: 0 of 545.
   * done 10-20: the unread fixed-atom force zeroing in `apply_fix_rigid_md` deleted (a 1e30-poisoned variant left
     6.6 and 7.0 bitwise unchanged), fixed-atom indices validated where the lists are built; regression w23 0 of 545.
   * hybrid preparation: applied (A70-A73) with its follow-up: the overdamped overlap_settling friction deleted
     (every prep MD stage at dry-MARTINI's own Langevin friction, mass / 4 ps; lipids now move 1.1-1.3 A per stage
     in 6.4-6.6, not 0.12); per-lipid diffusion targets (user 10-20: yes) in `LIPID_LATERAL_DIFFUSION_UM2_S`: DOPC 11.5, POPE 8.6,
     POPG 11.5 um^2/s, from ratios to DOPC (findings 2.5; no same-method set exists). Ions keep the DOPC clock for
     now. Before the dry-MARTINI patch the prep 1rkl box expanded 100.2 to 104.5 A (about 82 A^2 per lipid against
     65.5); the 1rkl test above re-measures it, and a glpG prep through production is still to run.
   * hybrid analysis: applied (A80); peptide numbering (A81, 10-21). The local drivers
     `scratchpad/local_popg_79HIS/{hdx,acc}.sh` are retired, not repointed: they call the deleted
     `write_hybrid_energy.py`, their trajectories are gone, and `run_hdx_analysis.sh` runs a ladder
     (`HYBRID_SOURCE_PATTERN` with `{replica}`).
   * D9: the sequence exclusions (|i-j| <= 1 or 2 by residue id) of `environment_coverage`, `protein_hbond`,
     `hbond_coverage`, `hbond_coverage_hydrophobe` and the `hb_environment_coverage` nodes also exclude the last
     residues of one chain against the first of the next (1dfn: 169 element pairs over 7 nodes, index oracle). On 1dfn's native dimer those pairs sit 14-16 A apart (CA), so the fixture is
     unaffected; the effect needs chain ends in contact. A53's writer-only `sequence_id` does not carry over:
     `EnvironmentCoverage` also indexes `aa_types` by `id2` (`environment.cpp:223`). Design open: the
     interaction graph keeps the index id for edges and takes a separate sequence id for its exclusion test.
   * From pass 2 (reports `scratchpad/upside3_audit/pass2/<unit>/report.md`), open [ ]:
     * the dry-MARTINI pair spline (`LayeredClampedSpline1D`) fits zero end slopes, so the force is 0 at the 1.2 nm
       cutoff where the potential-shifted LJ has up to 2.8e-2 E_up/A: against CLAUDE.md's spline-table rule; fit the
       grids with the analytic end slope (igraph D5, P2 MARTINI unit);
     * per-residue rotamer logs put a 1-rotamer/multi-rotamer pair's energy on one residue (BP D4; sums right,
       nothing reads them);
     * `HBondEnergy` input ranges unchecked (no writer produces bad values);
     * dead code: `run_upside` `read_hb`, `read_constant_hb`, `compute_topology`; `ProteinHBondInteraction::param_deriv`
       (returns -1); `write_infer_H_O`'s multichain check; `interaction_graph.h` `message()`,
       `PosDistSplineInteraction`, the `PARAM_7A/10A_CUTOFF` branches, `spline.cpp` `evaluate_*`,
       `deBoor2d_coeff_deriv`; `rotamer.cpp` `energy_cap*`, commented-out blocks, unread freshness flags, a doubled
       `edges11.reset()`; `SigmoidEnergy`'s `dim1` assert and its writerless `use_cutoff` 0 branch.
   Not changed on purpose: the H-bond 3.5 A cutoff (findings 3.14).
   Regression baseline: `scratchpad/upside3_audit/w36_hb.npz` on `fixtures_w6` (after A101, built in
   `scratchpad/upside3_audit/bp_home`; obj/ still holds A97 while the 10-21 glpG validation prep runs, install by
   rename after it [ ]; `regress.py compare`
   lists keys found in only one record apart from the arrays that differ).
   For the user: `membranelateral_potential` (master too) applies an xy force on surface CBs and reports potential 0,
   so the logged energy, pivot MC and replica exchange do not see it; make it the energy whose gradient it is, or
   keep it as a documented external force.
5. [ ] Confirm each claim by measurement, and record refuted claims.
6. [ ] Fix confirmed bugs one at a time. Each fix needs its reproducer passing, its node passing finite differences,
   and every other regression output bitwise unchanged.
7. [ ] Review pass 2, independent, on the trained-node units; then the P2 units.
   * Started 10-21 07:40: interaction graph, environment, BP, H-bond, placement, config writers, one reviewer each
     from `scratchpad/upside3_audit/pass2/brief.md`; reports in `pass2/<unit>/report.md`, no fixes applied by them.
8. [ ] Later, not now: refuse ff2.1.

### Phase 25 - are the finished force fields done training; the side-chain verdict (DONE 2026-10-20, user)
User 10-20: "check if the finished ffs needs no more trainings (converge, etc)"; "check if the deterministic /
probabilistic side chain have a conclusion". Read only, but for one read-only `panel_compare` job (1 CPU, 4 min).
Working files in `/project/trsosnic/yinhan/checks/converge_20261020/`.
1. [x] The gate on `ff31_depth_data` (ff_3.1_depth_data) and `ff31_depth_all` over steps 57-75: both CONVERGED
   (lowest dhb 0.033, 0.032).
2. [x] Their depths: no pull, whole, by basin or on the glycine maps (`depth_pull.py`). ff_3.0_gdepth's glycine
   alpha_L depth is pulled one way in all four epochs but passes the cut over the last (p 0.10).
3. [x] The defect left no detectable pull: the control's first fixed-engine epoch and the ff3.2 rerun's epoch 0 (the
   corrected gradient at ff3.0's parameters) show no angular pull; `mdc_ctrl_03` equals `mdc_ctrl38` (-0.012, p 0.67).
   Found on the way: no-BP's ~110 off-diagonal angular coefficients are still pulled at step 76, which its
   whole-group `rot` test passes.
4. [x] `panel_compare` 49435922 (4 min): no-BP -0.111 [-0.160, -0.065] folded, p 1e-5, at the fixed points.

**Review (findings 1.21, 1.36).** Released and finished force fields need no further training by the gate's rule:
ff_3.0_gdepth, ff_3.1_depth_data, ff31_depth_all and the control; the defect's direct test for ff3.1 is its rerun
(step 76 on 10-21). Two exceptions, both reported to the user: ff_3.0_gdepth's glycine alpha_L depth drifts one way
(+0.015 nats a round), and no-BP's ~110 angular coefficients are not at a fixed point. The side-chain verdict is in:
at both fixed points the deterministic side chain folds 0.241 against 0.353 on 89 held-out domains. Choosing the
side-chain model for Upside 3 (upside3.md gate) is the user's.

### Phase 26 - the side-chain comparison on one metric, TM-score (DONE 2026-10-21, user)
User 10-21: "I prefer to use a single metric (TM score?) for all tests". TM-score as the Peng benchmark is scored
(`ff3_benchmark/scoring/tmscore.py`: fixed correspondence, superposition optimised, step 2, normalised by the full
chain): per domain, 4 seeds pooled, frames after T_EQUIL, mean TM and the fraction of frames at TM > 0.5.
1. [x] mdCATH at the fixed points (`mdc_nobp_03`, `mdc_ctrl_03`, the only panel tags whose trajectories are kept):
   `checks/converge_20261020/tm_panel.py`, one broadwl node, every 10th frame (80 a seed; 2.7 s a frame at 200
   residues).
2. [x] The 44-domain panel: `nb_01`/`ctrl38` kept no trajectories. User 10-21: run both fixed points. Tags
   `p44_nobp_03`, `p44_ctrl_03` in `ff3_selection/runs/`, params `params/mdc_<nobp|ctrl>_03` (the gate's copies),
   at the training time step 0.009 as the mdCATH panel (the earlier 44-domain runs used 0.015); one node each
   (`checks/converge_20261020/panel44_dt009.sbatch`), then `tm_panel.py`.
3. [x] Peng's benchmark: TM for the 38-step releases, the last pair (ctrl38 ubiquitin de novo) scored 10-21; the
   fixed points not run (user: "44-domain panel only").

**Review (findings 1.21).** Mean TM to the native: mdCATH -0.034 (p 7e-8), the 44 domains -0.026 (p 0.007), Peng
natives -0.029 (p 0.18), de novo +0.010 (p 0.49); TM > 0.5 shares equal within error.

### Phase 27 - released ff2.1 on the same three tests (STARTED 2026-10-21, user)
User 10-21: "run validation with released ff2.1 so we know if we should keep ff2.1 or not". Both side-chain arms are
ConDiv retrainings of ff2.1 (`ff_2.1_ctrl38` shares none of ff2.1's trained files, only its Ramachandran library), and
both sit below Peng's published FF2 on his benchmark (findings 1.21), so neither answers whether ff2.1 itself should
stay. Decisions:
- **The released files** are `ff3_selection/params/ff21_released` (md5 equal to `/beagle3` `parameters/ff_2.1` and
  `common/rama.dat`); every run uses each test's own protocol, tool and binary, unchanged, so ff2.1 pairs with the
  retrained arms on the same domains, seeds and time step.
- **Time step 0.009 everywhere**, as the retrained arms' panels and Peng's protocol ran; the 44-domain `ff21_released`
  panel of 10-01 ran at 0.015 and kept no trajectories, so it is rerun.
- **Peng benchmark:** the 30 arms not yet run with `bench.sbatch` as `validate_peng.sh` submits them; the two lambda
  arms of 10-01 to 10-04 (11d) are the same protocol and binary (`/beagle3` `obj/` unchanged since 09-30) and complete.
- **Analysis runs by itself** (no agent polls): the panels' TM jobs wait on their panels by dependency; a 1-CPU
  collector job wakes every 6 h, scores each finished arm with `score_arms_ff.py`, and when all 32 are scored writes
  the three-way table and stops. Files in `/project/trsosnic/yinhan/checks/ff21_released_20261021/`.
1. [ ] mdCATH panel `mdc_ff21_released` (89 domains x 4 seeds) and the 44-domain panel `p44_ff21_released`.
2. [ ] TM against the fixed points `mdc_ctrl_03` / `mdc_nobp_03` and `p44_ctrl_03` / `p44_nobp_03` (`tm_panel.py`).
   44 domains done: released ff2.1 0.681, control 0.650, no-BP 0.623 (findings 1.37); mdCATH running.
3. [ ] Peng's 30 remaining arms; scored and compared against the step-37 `ff_2.1_ctrl38` and `ff_2.1_nobp` and
   against Peng's published FF2 (0.55 native, 0.42 de novo).
4. [x] User 10-21: does the loss reach the other trained force fields? 44-domain panel at 0.015, every epoch end
   against `ff21_released` (`panel44_trajectory.py`): every run of this workflow folds 0.13-0.23 less, most of it in
   epoch 0 (findings 1.37).
5. [ ] The same table with the owed panels (`c9_02`, `c9_03`, `a31_00`-`03`, `a32_00`-`03`; job `panel44_traj`
   writes `panel44_trajectory_final.txt`), and Peng's arms of ff_3.0_gdepth and ff_3.1_depth_data scored
   (`peng_trained.txt`; the collector's final table includes them).
6. [ ] Review: keep ff2.1, or a retrained force field, by the three tests.

### Phase 28 - is retraining's folding loss a workflow defect (STARTED 2026-10-21, user)
User 10-21: "is this a defect in our training workflow? can you fix it?" Known: every run loses 0.13-0.23 folded in
epoch 0 (findings 1.37); in round 3 ff2.1's side-chain table put back restored it (1.17 A), but round 4 trains that
table at a tenth of the rate and loses more. Peng's FF2 ran 76 iterations with no convergence test and his trainer is
lost (9v), so a deviation of our port cannot be read off his code; the SI-rate twin (`ff30_gdepth_si`, half our
hb/env/sheet/dhb rates) lost as much (`ds_02` 0.406), so the rates set speed, not the loss. Decisions:
- **Attribute before fixing**, by the round-3 method on the round-4 control's epoch 0 (`c9_00`, ff2.1's library, so
  every difference from ff2.1 is training): one file at a time put back to ff2.1's, 44-domain panel at 0.015.
- **No workflow change without the user.** Every candidate fix (rates, a prior toward ff2.1, freezing groups, early
  stopping on a validation panel) departs from the FF2 workflow the project fixes (findings 10.17, memory
  `condiv-must-be-ff21-workflow`); a port defect against the SI is fixed directly.
1. [x] Panels `c9_00_hb21`, `c9_00_rot21`, `c9_00_env21`, `c9_00_bb21`, `c9_00_sheet21` (`resets_c9_00.txt` in
   `checks/ff21_released_20261021/`): no single file restores ff2.1's folding; the side-chain table (+0.089), the
   backbone burial (+0.081) and the H-bond (+0.068) each give back part and together about all of the 0.225 lost;
   environment (+0.028) and sheet (+0.035) are not resolved.
2. [x] Each group's course over the control's epoch 0: hb and env pulled one way from step 0, rot spreads as noise;
   training proteins' RMSD flat over four epochs (findings 1.37).
3. [x] The ensemble definitions against the SI: ours are Kleinmann's (Peng's code), so no port defect found
   (findings 1.37). The loss holds at the training time step 0.009: released ff2.1 0.603 against the control's
   fixed point 0.430 (all 44 domains, findings 1.37).
4. [x] How Peng and Jumper chose where to stop (user 10-21): a fixed budget short of convergence on purpose, Jumper
   having seen accuracy fall when training continued (findings 1.23). Our gate is our own addition, and every run
   started from ff2.1, so each is a second budget on top of FF2's.
5. [ ] Fix options to the user, with the file that carries the loss (user 10-21: "check if we should correct current
   trained / training ffs"). Every force field trained from ff2.1 by this workflow is affected; nothing is training
   now; the 8 glpG chains run ff3.0 and ff3.1.
6. [x] Why ff3.1 and ff3.2 end at their starting maps (user 10-21): each depth sits where its 0.2-nat prior,
   centred on the start, balances the data; the glycine pull is consistent only pooled over the 37 maps. From ff2.1
   the data pull glycine alpha_L shallower (t +5.6), from ff3.0's map not at all (findings 1.34). A retraining
   that should move the maps needs a pooled prior (one offset per central residue and basin, maps as deviations),
   which 1.29 found to transfer best; the user's decision.

### Phase 29 - ff3.1.1 and ff3.2.1: the depth training with a pooled prior (STARTED 2026-10-21, user)
User 10-21: "apply the fix and train ff3.1.1 and ff3.2.1", after Phase 28 step 6 found each depth held at its start by
its own 0.2-nat prior while the data's pull is clear only pooled. The user chose (10-21, asked): everything else trains
as in ff3.1, so the pooled prior is the only change. Decisions:
- **The model** (`training/ff31/rama_basin.py`, mode `pooled`): a map's offset on a basin is its central residue's
  shared offset plus the map's own deviation, c_mb = q_X(m),b + d_mb. The deviations are `data` mode's, at its width
  0.2 and on its 1,442 depths. The shared offsets pool the maps of one central residue only (rule 1.8: central GLY
  says nothing about central ALA; the coupling is estimation, as NDRD's hierarchical prior), each under its own
  Gaussian prior, trainable where the residue's pooled native sites reach `data`'s threshold of 20. Pre-proline maps
  (X|right|PRO) keep only their own deviations, as ff3.0's glycine pair left them out (findings 1.14). Same damped
  Newton step per epoch, solved jointly over shared offsets and deviations.
- **The shared offsets' width is cross-validated**, as the per-map widths were (findings 1.29): 10 folds over the
  proteins of each source run's epoch 0, the trainer's own step and held-out map-level gap, widths 0 (= `data`) to 2.
- **Each run is its source with one change**: ff3.1.1 = `ff31_depth_data_angfix` (start ff2.1), ff3.2.1 =
  `ff32_depth_data_angfix` (start ff3.0), built as `build_rerun.sh` builds (same `ConDiv.py` md5 `df3e2497...`,
  inputs, minibatch order, held-out set, fixed engine), only `rama_basin.py` replaced. Target 76, no gate, stop for
  the user; epoch-end panels and per-step steric checks as for the reruns.
- **Verification before launch, no test jobs** (findings 10.16): with the shared offsets off, mode `pooled` gives
  `data`'s step exactly, on recorded epoch-0 data and against the source run's own round 1; P and its transpose are
  adjoint; the initial checkpoints equal the sources' but for the new depth fields.
1. [x] `rama_basin.py` mode `pooled` (findings 1.34).
2. [x] Equivalence and adjointness checks: `data` reproduced to 1e-17, adjoint to 4e-16.
3. [x] Cross-validated width of the shared offsets: 0.5 (findings 1.34).
4. [x] Built and initialised `ff311_depth_pooled` and `ff321_depth_pooled`: equal to their sources in every field
   but the pooled depth state; submitted 13:34 10-21 (remote_jobs.md §1).
5. [ ] At step 76: maps against their sources (`plot_gly_rama_*_nonhelix_zero.py`), `depth_data_vs_prior.py`, panels.
   Both runs reached step 76 10-23. Glycine maps done 10-23 for the user's review (findings 1.34: ff3.1.1's alpha_L
   +0.41 from ff2.1, ff3.2.1 at ff3.0's map); panels and the per-step checks wait on it, since their collector died
   on its first pass (remote_jobs.md §1 line 2). Two-start comparison done 10-23 (findings 1.34): every map but
   glycine's alpha_L agrees to the runs' noise; that gap went 0.517 to 0.179 and still closes ~0.07 an epoch.

**Revised decision (10-21 13:40, user: "besides rama map not moving bug, try to resolve the accuracy lost bug as well
in this run").** The first chains (`ff311_depth_pooled`, `ff321_depth_pooled`, every group trained) were cancelled 3
minutes into step 0 and kept as `<run>.cancelled_20261021`. Reasons, from Phase 28: the loss is spread over the
co-trained side-chain table, backbone burial and H-bond (findings 1.37); holding some groups while the rest train
lost more folding before (the 10-07 frozen arms, findings 1.21); and Peng and Jumper stopped those terms on a fixed
budget that ff2.1 has already had (findings 1.23). So:
- **Every Adam group is held at ff2.1** (rate 0, `checks/ff311_20261021/freeze_groups.py` on the initial checkpoint;
  the trainer is unchanged). Only the depths train, under the pooled prior: ff3.1.1 = ff2.1 plus depths from NDRD's
  map, ff3.2.1 = ff2.1 plus depths from ff3.0's map. ff3.0's other files carry the loss (folded 0.464), so ff3.2.1
  takes ff2.1's; the two runs differ only in the starting map, and share ff3.1's rerun's minibatch order and
  held-out set.
- **Runs** `ff311_depths_only` and `ff321_depths_only` (`build_depths_only.sh`), tags `p311_EE`, `p321_EE`.
- **The test of the fix** was to be the epoch-0 panel (`p311_00`, `p321_00`) against `ff21_released`.
6. [x] Superseded 10-23 by Phase 30 (user): the per-run panels are not run; the two runs train on until they agree,
   and the merged force field is what gets validated.

### Phase 30 - merge ff3.1.1 and ff3.2.1 into ff3.1 and validate it (STARTED 2026-10-23, user)
User 10-23: "train 3.1.1 and 3.2.1 further until they agree, then we can merge them into one (overwrite 3.1) and run
validation on it"; chosen: glpG hybrid included in the validation, merge by averaging. At step 76 every term agrees
between the runs to their own noise except glycine's alpha_L, whose shared-offset gap closes ~0.07 an epoch
(findings 1.34). Decisions:
- **Agreement rule, fixed before the data:** at one round, (a) glycine's alpha_L shared-offset gap is no larger than
  the rms of the other trainable shared-offset gaps (the runs' noise, 0.092 at round 4), and (b) the glycine maps'
  gap rms over the data-reached depths is no larger than the non-glycine maps'. Train in epoch pairs: 114, then 133;
  at 152 stop and report.
- **Merge:** the mean of the two runs' offsets per (map, basin) in the common frame (NDRD plus offsets, ff3.2.1's
  including ff3.0's pair), written by `rama_basin.write_library`; the other five files are ff2.1's.
- **Release as ff3.1:** overwrite `parameters/ff_3.1_depth_data` (backups) through `validate_ff.sh`, its third
  argument extended to take a ready force-field directory; the old ff3.1's glpG chains, seeds and Peng scores move
  aside as `.bak_pre_merge_<date>`, after the running chains' block 5 ends.
- **Validation:** released ff2.1's three soluble tests on the same binaries and protocols (44-domain and mdCATH
  panels at dt 0.009 with TM comparisons against `ff21_released`, Peng's 32 arms) and validate_ff.sh's four glpG
  REMD chains (TM4 at each block end, HDX on block 2).
1. [x] Dropped p311/p321 from `submit_new.sh` (backup `.bak_pre_phase30_20261023`); both runs continued to step 114,
   49451890 / 49451891, 15:33 10-23 (remote_jobs.md §1 line 2).
2. [ ] At 114: KE scan, `check_step.py`, the agreement rule (`compare_pooled_runs.py` over every round); extend if
   not met.
3. [ ] Merge, check, glycine figure to ~/Downloads; move the old ff3.1 results aside; release by `validate_ff.sh`.
4. [ ] Panels, TM jobs and the Peng collector (`checks/ff31_merged_<date>/submit.sh`); read each as it finishes.

## Known Errors / Blockers

* **Every hybrid dry-MARTINI result so far carries these defects** (audit 10-20; fixes in preparation, plan.md Phase 24).
  * The side-chain/lipid interface table treats the choline bead NC3 (+1) as neutral: the builder takes charge from
    the bead type (Q0 gives 0), not the topology. ASP-NC3 misses -2.15, -1.46, -0.91 E_up at 5, 6, 7 A; LYS-NC3 at
    6 A has the wrong sign (-0.38 for +0.84).
  * The prep stages' NPT barostat never scaled a lipid: z-fixed phosphates were treated as fully fixed. So the box
    was never relaxed from CHARMM-GUI's.
  * The hybrid rotamer 1-body is not ff2.1's: no library has `rotamer_prob_fixed`, so the prep silently writes
    -log of the unweighted phi/psi-grid average instead of the phi/psi-dependent term. Relative rotamer energies in
    a helix move by up to 9.1 (LEU), 7.8 (ILE), 7.4 (PHE) E_up. The live glpG REMD seeds
    (`popepopg_REMD_mdw2/seeds/`, checked 10-20) carry `placement_fixed_scalar` and no `placement_scalar`, while
    master's ff2.1 workflows use `--dynamic-rotamer-1body`: the running chains have the fixed 1-body.
  * The prep never removes CHARMM-GUI's bead overlaps (its hard minimisation is off by default), and the PO4
    restraint's reference sits 53 A off in z, so only the hard z freeze holds the heads. Local preps of 1rkl (4)
    and glpG (1) heat to KE/1.5kT 26-33 in stage 6.3 and the box jumps by 6-7 A; in 1rkl's production the area per
    lipid is 91-100 A^2 against this force field's 65.5. P-P stays at the template's 36.4 A (force field: 43.3).
    The cluster glpG preps (`popepopg_glpG/<V>/logs`, 08-14, read 10-20) did not overheat: KE/1.5kT 0.83-0.87 in
    every stage. They carry the restraint offset (logged potential +326k E_up at 6.2, +153k at 6.3, scaling with
    the force constant) with the heads frozen, so their thickness and box defects apply.
    `scratchpad/relax_popepopg.sh` calls the deleted `prepare --mode bilayer`.
  * The ion count neutralises residue charges the protein does not carry: net +3 e (1rkl), -6 e (glpG).
  * The same table read its angle mirrored in the engine and placed the CB 0.57 A off (A76): a bead 4.5 A from a
    LEU CB on the side-chain side got the backbone side's energy.
  * Every ALA and GLY side chain was missing from the protein's side-chain node, its placement and the H-bond
    coverage (A83, found 10-21): glpG 42 of 210 residues, 1rkl 3 of 31. Their pair energies with every other side
    chain were absent, under every force field, on the cluster and locally.
  * The prep's NPT stages total 4.6 ps against CHARMM-GUI's 10 ns, so production starts at an unrelaxed box: 1rkl's
    lipids 10% over their equilibrium area, glpG with an unhealed gap at its belt (Phase 24).
  * Fixed in the repo 10-20 (A70-A76): `parameters/ff_2.1/martini.h5` rebuilt (md5 `cc90a933...`; the old
    `2e612261...` kept as `martini.h5.bak_20261020`; ff_3.0 configs use it too). Every prepared hybrid `.up` must be
    re-prepared; the patched engine refuses the old table, and since A103 (10-21) also any stage file without
    `combined_energy_end_slopes`. The cluster trees running the campaign are unchanged.
* **Hybrid HDX lipid accessibility made before 2026-10-20 is shifted by one residue** (findings 3.14 A61, fixed 10-20).
  glpG, numbered from 1, gave each amide its N-terminal neighbour's accessibility. To regenerate (on the cluster, by
  the campaign's drivers):
  * `/project/trsosnic/yinhan/popepopg_REMD/<V>/hdx*/results/<V>_popg_remd_<i>_ACC.npy`;
  * the combined `_PS.npy` (`_PS_protein.npy` stays valid);
  * every calc_hdx_ht output in those directories (`_dG_profiles.npz`, `.csv`, `_ranked.csv`, `_<T>.csv`, PNGs);
  * anything `scratchpad/local_popg_79HIS/acc.sh` produced; findings 5.3's `acc` values.
  Local drivers that call the deleted `example/00` helpers need repointing: `scratchpad/local_popg_79HIS/hdx.sh`,
  `hdx_cluster.sbatch`, three `scratchpad/np_slurm/*.sh` (`get_info` is now `py/get_info_from_upside_traj.py`; the
  `write_hybrid_energy` step is dropped, since the projection copies the coupled potential and the MBAR consumers
  subtract the mean).

* **Every ConDiv run before 2026-10-18 08:35 trained the side-chain pair term's off-diagonal angular coefficients on
  a mixed gradient** (findings 3.12; fixed in the engine then): ff2.1's retrainings to step 37, ff3.0, ff3.1, ff3.2.
  The corrected gradient shows no pull at ff3.0's or the control's fixed point (Phase 25); ff3.1's rerun is the
  direct test for ff3.1 (Phase 21); ff3.0 and ff2.1 are not retrained.

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
* **HDX post-processing items await the user's decision** (user, 10-11; findings 5.5): three docstring/README
  statements that call TM helices non-exchanging, and the unmeasured interface boundary of the 7 A tail-contact
  criterion. The PF = 1000 uptake cap was replaced in Upside 3 by the per-temperature ESS bound with censored
  amides (A64, 10-20, under the user's 10-19 direction to fix every bug in Upside 3; master unchanged), told to
  the user 10-20. Do not compare uptake with the HXMS data before the remaining items are settled.
* **lambda's helix 2 / helix 4 misorientation** is a known ff2.1 limitation with no single
  side-chain cause (findings 11d). It is re-tested after ff3.0, not trained against.
