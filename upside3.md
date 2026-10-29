# Upside 3

**Status.** The engine audit (goal 1) met its closing criteria on 2026-10-22 (2.2 item 5; plan.md Phase 24 item 2a);
the defect ledger is `findings.md` 3.14. Refusing ff2.1 is next. The other goals are proposals.

The user set five goals on 2026-10-19 (section 1). The protein-model work is built on the side-chain model that wins
the verdict, deterministic (no-BP) against probabilistic (BP) side chains at their converged epochs (plan.md Phase 21
item 7). Measured 10-20 (findings 1.21): at both fixed points BP folds 0.353 against no-BP's 0.241 on 89 held-out
mdCATH domains (p 1e-5); the choice is the user's. The audit and the units work do not depend on it. Everything
measured below is ff_2.1, the probabilistic side chain, unless named.

Records: plan.md Phases 22-24. The planning reports for goals 1, 2 and 5 (`audit_report.md`, `units_report.md`,
`gpu_report.md`, with scripts and logs) are in `scratchpad/upside3_planning_20261019/` (gitignored); the large-protein
investigation is in `scratchpad/large_protein_20261019/` on the Mac Studio (76 G with trajectories). Section 8 lists
both.

## 1. Goals and order

1. **Debug the whole repo.** Find every bug ourselves, as Jumper did with Claude on master; master's code first
   (section 2).
2. **Unit alignment.** Run Upside in standard units without conversion, if that can be done (section 3).
3. **New parameters,** chosen after the side-chain verdict (section 4).
4. **Fold large proteins de novo** (section 5).
5. **GPU support,** with one force field for CPU and GPU, never trained twice (section 6).

| step | work | waits on |
|---|---|---|
| 0 | side-chain verdict | the user's choice (measured 10-20, findings 1.21) |
| 1 | audit and fixes (goal 1); units (goal 2) | nothing |
| 2 | large-protein baseline (goal 4, Phase A) | steps 0 and 1 |
| 3 | training changes (goal 4, Phase B) and new parameters (goal 3), one change per training | step 2 |
| 4 | validation, including Upside 3's own temperature calibration (section 3) | step 3 |
| G | GPU (goal 5), in parallel from the end of step 1 | step 1 for the refactor; steps 0 and 3 for the BP kernels and any new term |

Why this order: every Upside 3 training has to run on the audited engine, since each defect found was carried by the
trainings that ran before it (2.1) and a training is weeks of cluster time; and the GPU port is written once, against
the audited engine in its final units and with its final set of terms.

**Upside 3 is a new release, not a patch series on ff2.1** (user, 10-19). No list of fixed bugs and their effects on
ff2.1 is published, so no fix is measured against ff2.1 or master for attribution. Each fix is still proven correct
(2.2); the debug as a whole is validated on ff2.1 runs (plan.md Phase 24 item 2b), and Upside 3 as a whole in 5.6
Phase C. Master and ff2.1 stay as released.

## 2. Engine audit (goal 1)

### 2.1 Defects known before the audit; status

The closing criteria (2.2 item 5) were met on 10-22 (`scratchpad/upside3_audit/closing/`). Every confirmed defect,
with its effect and verification, is in the ledger (`findings.md` 3.14); open items are in plan.md Phase 24.

Defects that trained force fields carry, found before the audit:
* the BP stop test (`rotamer.cpp:279`, missing `fabs`; Jumper's report, fixed 2026-09-26): ff2.1 and ff_3.0;
* findings 3.12, the pair term's angular parameter gradient (`interaction_graph.h`): every ConDiv run before 10-18;
* findings 3.13, `environment_coverage_hb`'s covering-residue types (section 5.4): ff2.1 and every training so far.

The first pass (10-19, `audit_report.md`) seeded the audit: build flags that compiled `isfinite` out on the Mac only,
an AddressSanitizer read of atom -1 in `DihedralCoord`, three bugs in `Concat`, candidates in `CurvatureChange` and
the thermostat, 708 compiler warnings with no real bug beyond `Concat`'s, and no test suite (the harnesses of 2.2
item 1 are now `tests/`). Its items went into the audit's units; the confirmed ones are ledger entries, for example
A1 (`DihedralCoord`), A25 (`make_restraint_group`), A28 (the Mac build flags), A29-A30 (the incomplete backbone
desolvation and H-bond parameter derivatives), A37 (`CurvatureChange`), A55 (`Concat`), A62 and A64 (the HDX
analysis factors), A68 (the hybrid stage-string trap). Still open: `run_remd.py` rolls NaN replicas back and continues, a guard under NO GUARDS, kept in
the running Upside 2 campaign (plan.md Known Errors; section 7 item 4).

### 2.2 Method

1. **Harnesses first:**
   * fixtures covering every `advanced_config` option, plus a regenerated hybrid seed;
   * per-node finite differences of forces and of parameter derivatives, with tie maps, against a float64 oracle;
   * a numpy oracle of each writer's index meaning;
   * a bitwise regression harness, so a fix is shown to change only what it means to;
   * a sanitizer and build-flag matrix;
   * invariance and statistical tests: equipartition, MC detailed balance, replica exchange on a harmonic ladder.
2. **Units of review.**
   * 18 units cover every master file exactly once (P1): engine graph, coordinate nodes, placement, BP, interaction
     graph and splines, H-bond, environment, local backbone terms, membrane, driver and replica exchange, MC, I/O,
     config writers, parameter glue, analysis, examples, build, and shipped parameters.
   * 7 units cover the branch-only code (P2): MARTINI potential, hybrid coupling, g-JF integrator, box and barostat,
     hybrid preparation, hybrid analysis, and ConDiv.
   * A writer/reader registry lists 80 node prefixes, so each node is reviewed from its C++ reader and its Python
     writer together. Findings 3.13 and 4.17 were writer/reader mismatches; findings 3.12 sat in a
     parameter-derivative path that energies and forces never use.
3. **Two independent reviewers per unit.** The trained nodes go first: interaction graph, environment, BP, H-bond,
   placement, config writer.
4. **A claimed bug counts only after someone other than its finder reproduces it,** with a failing measurement or a
   minimal reproducer. Refuted claims are recorded too.
5. **A unit closes** when both reviews are in; every node passes force and parameter finite differences and the
   invariance tests; the sanitizers are clean; and every registry pair is checked.
6. **Fixes go into the Upside 3 tree one at a time.** Each fix needs its reproducer to pass, its node to pass the
   finite-difference tests, and the regression harness to show every other path bitwise unchanged. Its effect on
   ff2.1 is not measured (section 1).

## 3. Units (goal 2)

### 3.1 What was found (`units_report.md`)

**Upside's units are reduced.** kT at T_up = 1 is 1 E_up; lengths are in A; every backbone atom has mass 1 m_up
(12 g/mol). The implied time unit is 1 A x sqrt(12 g/mol / 2.915 kJ/mol) = 0.203 ps, so dt 0.009 is 1.83 fs.
Master's README says Upside times "should not be interpreted as picoseconds". The hybrid runs its own calibrated clock
of 40 ps per step, whose origin is not recorded.

**The kJ/mol value of E_up is a calibration, not a unit.** 350.588 K is 298 K / 0.85: Peng 2022 SI eq. S14 matches
T_up 0.85 to 298 K for FF2 (FF1's match was 1.0, which gives 2.479 kJ/mol per E_up). Writing the parameter files in
kJ/mol would bake one force field's fit into every number, and a retrained Upside 3 needs its own calibration. The
MARTINI conversion factor is this same calibration: the hybrid's protein-lipid balance depends on it.

**Some terms are evaluated at a fixed reference kT, at any simulation temperature:** the Ramachandran `-ln P` and its
reference state, the sheet weights, the rotamer `-ln p`, the BP free energy (`rotamer.cpp:252, 300, 466, 896`) and
`weighted_pos` (`environment.cpp:518`). Any change of units keeps them at kT_ref; replacing kT_ref with kB T_sim would
change the model.

**Where units sit in the code:** in C++, 45 implicit sites in 20 files (its one explicit constant, `FORCE_UNIT`, was
corrected in A38 and deleted with its only user in A170); in Python, about 67 groups of sites in about 30 files.

**Defects found.** The conflicting kelvin conventions (0.593 kcal/mol and 41.4 pN per E_up/A, which imply 298 K and
300 K) were fixed in the audit (A62), and up.md's `/output` errors too (A182). Left for the units work:
`martini_prepare_system_lib.py` hard-codes `*10.0`, bypassing the explicit length conversion, and several conversion
literals are duplicated.

### 3.2 Options

| | engine | parameter files | dry-MARTINI conversion | ConDiv and past results | old `.up` files |
|---|---|---|---|---|---|
| A+: reduced units inside, standard units at every input and output | unchanged, bitwise | unchanged | stays, with the factor taken from the force field's recorded calibration | unchanged | run as before |
| B: kJ/mol, K, A, amu | about 28 sites; not a coherent set (the time unit is 0.1 ps) | every energy x 2.915 | length factor only | learning rates rescaled; every E_up number relabelled | converter and format marker |
| C: GROMACS units (nm, ps, kJ/mol, amu, K) | about 42 sites, every length included | energies and lengths | none, except the calibration | as B | as B |
| D: kcal/mol, A | as B | as B | still both factors | as B | as B |

**Recommended: A+.** Conversion moves out of the engine into one units `Config` module; each force field records its
calibration (`T_ref_K`) in its parameter directory; every `.up` file carries `/input/units`; input and analysis take
and report K, kJ/mol or kcal/mol, bar and pN through that module. The recommendation rests on the calibration (kJ/mol
in the files would carry one force field's fit), not on A+ leaving the engine bitwise unchanged, which a new release
no longer requires. C is the coherent target only if Upside 3 is retrained from scratch with no comparison to ff_2.1,
or the hybrid becomes the main product; even then the kT_ref terms keep an explicit reference temperature.

### 3.3 Plan

1. [ ] The units module and `/input/units`.
2. [ ] Every Python site moved onto the module.
3. [ ] The hard-coded `*10.0` and the duplicated literals (3.1).
4. [ ] Verification:
   * engine bitwise on the 8-config x 3-frame test;
   * each derived constant equal to the literal it replaces to 1e-12, choosing explicitly where repo literals
     disagree (14.521180763676 is off by 1.3e-7);
   * configs regenerated by the old and new writers byte-identical in `/input`.
5. [ ] Upside 3's temperature calibration at validation, by Peng's eq. S14 procedure. Rebuild the MARTINI tables if it
   moves.

## 4. New parameters (goal 3)

Chosen after the side-chain verdict and the Phase A baseline (section 5.6). Each one:

* answers a measured defect, with its predicted effect written down before training;
* goes in the term that carries the defect (findings 10.13);
* is trained as the only change against a control (findings 10.17);
* passes the audit's parameter-derivative harness before any training;
* is a spline table read from `.h5`, so the GPU kernels can carry it too.

Candidates on record:

* **Existing parameters held frozen:** H-bond entries 4-11 and backbone desolvation's centre and sharpness. They have
  exact derivatives since A29-A30 and are held frozen by the trainer (`hbr` at learning rate 0; `bbenvc`, `bbenvs`),
  so training them needs no new term.
* **From goal 4 (5.7, Phase B3):** one table of burial and pairing for the backbone polar groups, contacts that depend
  on their surroundings, burial weights by who covers a residue, the secondary-structure balance in compact states,
  and temperature-dependent burial.
* **The glycine context term** (`architecture.md`), only if its trigger fires.
* **The side-chain model's own parameters,** which follow the verdict.
* **Training terms for HDX, large proteins and IDRs** (4.1); what disordered regions need (4.2).

### 4.1 Training terms for HDX, then large proteins, then IDRs (2026-10-28)

User 10-28: change what the training fits, with HDX first, large-protein performance second and IDR sequences third.
The energy's functional form stays spline tables. The change is in the objective: terms added to ConDiv's contrast.

**Why ConDiv does not fit HDX.** An amide's dG_HX is the free energy of the state in which it is open,
-RT ln K_op. ConDiv's objective holds no free energy. NSE asks that free replicas started at the native stay where the
native-restrained replica is, and DSE asks that unfolded frames match a chain with no attractions. Neither sets how
stable the fold is or which part of it opens first, so HDX agreement is whatever that fixed point leaves (Upside FF2
on Peng's three proteins: censored Spearman 0.61, 0.43 and 0.32).

**The terms.** Each is built from the per-frame dV/da the engine already returns (`ConDiv.py` `compute_divergence`):
* HDX, per measured amide i at the experiment's temperature, with MBAR weights over the replica ladder:
  d dG_i / da = <dV/da>_open,i - <dV/da>_closed,i, a contrast of the same kind as NSE. The loss is
  sum_i (dG_i - dG_i,exp)^2; an amide the model censors enters one-sided. Global stability is the same contrast between
  unfolded and folded frames.
* IDR size: d<Rg>/da = -(<Rg dV/da> - <Rg><dV/da>) / T (4.2, I1).
* Large proteins: B1 and B2 of 5.6 already change the training (decoys through a replay buffer; crops past 150
  residues).

Each is trained as the only change against a control (section 4).

**Conditions, from M1** (ubiquitin, findings 5.7):
1. **Its own sampling.** One ConDiv step (12 x 8000 tu) is 1/319 of the REMD that resolves ubiquitin's amides to
   5.6 kcal/mol at 273 K. Each HDX protein therefore gets its own REMD, reused across updates by reweighting to the
   new parameters while its ESS holds (DiffTRe, Thaler and Zavadlav 2021; CALVADOS, Tesei et al. 2021), and is rerun
   when the ESS falls.
2. **Signal only below the bound.** Resolving experiment's 8.5 kcal/mol at 273 K would take about 200 such runs. The
   term acts on local and subglobal openings, and on small, marginally stable domains (HEEH_rd4_0097 at 2.2-4.3
   kcal/mol; the 28-64-residue domains of Ferrari et al.).
3. **A fixed temperature mapping.** Peng reads FF2 at 300 K against an experiment at 273 K, at the temperature where
   the model's global stability matches. A term trained at a chosen temperature therefore learns nothing about
   stability. The mapping is fixed with the units (section 3) before any HDX training.
4. **A frozen definition of open.** The estimator reads the trained burial weights. ff2.1 to ff3.0 moves dG by
   0.12 kcal/mol RMS through the definition alone. The estimator keeps ff2.1's tables, so training cannot redefine
   "open" instead of changing the ensemble.
5. **A ladder that mixes.** Each protein's ladder spans its own measured transition, with walkers making round trips
   across it. ff3.0's ubiquitin ladder makes 3 round trips and its dG is not converged.

**What the term has to move, by M1.**
* ff2.1 opens the seven amides experiment assigns to global unfolding by fluctuations of the folded state
  (3.5-7.0 kcal/mol, against 7.7-8.5).
* ff3.0's ubiquitin melts about 0.06 T_up lower than ff2.1's, which the benchmark's bottom-rung TM-score hides.

So the term must raise the folded state's opening free energies at those amides while holding global stability.

Steps (user 10-28: "record the plannings. We will run it later"; none starts without the user's word):
- [x] M1 (10-28): ubiquitin on the Peng-benchmark REMD; findings 5.7, `scratchpad/hdx_term_20261028/` (section 8).
- [ ] M2: M1's readout on ff3.1's ubiquitin arms (`/beagle3/trsosnic/yinhan/ff3_benchmark/runs/ubiquitin_{native,
  denovo}_ff_3.1`, running 10-28) once both hold a COMPLETE file. No new run.
  * Replace the arm and force-field names in `fetch.sh`, `build_inputs.sh`, `flags.sh`, `run.sh` and `analyze.py`
    (`ARMS`, `DENOVO`). The HDX topology takes `parameters/ff_3.1`. Then run `fetch.sh` and `run.sh`.
  * Keep M1's rules unchanged, and read against M1's ff2.1.
- [ ] M3: Upside 3 REMD of EHEE_rd2_0005, HEEH_rd4_0097 and Ub psWT, each under ff2.1 and the base force field
  (after plan.md Phase 33), from the native, on the Upside 3 tree (`env_upside3.sh`) through midway2. Read with M1's
  rules at each experiment's temperature (EHEE 298 K, HEEH 278 K, Ub 273 K) and at Peng's comparison temperatures.
  * **Protocol: Peng's for FF2** (SI pp. 12-13; Table S3).
    * EHEE and HEEH: 14 replicas on 0.80-0.995, 4,000,000 tu, the first 500 frames dropped.
    * Ub: three independent runs, 14 replicas on 0.80-1.02, at least 2,400,000 tu each, the first half dropped.
    * Frames every 100 tu in all of them.
    * The ladder changes only if a force field's transition leaves the range, and the change is recorded.
  * **Inputs.**
    * Peng's design models: `/project2/trsosnic/share/HEEH_rd4_0097.pdb` and
      `/project/trsosnic/yinhan/bioemu2_hdx/inputs/EHEE_rd2_0005_ref.pdb`.
    * Ub psWT: `/project/trsosnic/yinhan/bioemu2_hdx/inputs/1UBQ_psWT_ref.pdb`, 1UBQ with residues 45 and 68 renamed,
      as Peng started from 1UBQ.
  * **Cost, estimated from M1's arm.** That arm (one node, 14 CPUs, 55 us per system-step) took six 10.5 h chunks for
    2.43M tu. So each Ub run is about 60 node-hours, and EHEE and HEEH about 50 each at half the size and 1.6 times
    the length. That makes about 280 node-hours per force field, about 570 for both.
  * Run through the benchmark's `bench_run.py` chunking, or an equivalent driver, with the jobs recorded in
    `remote_jobs.md`.
  * This gives the baseline a trained term must beat, and the three proteins stay held out of training.
  * No Upside run of EHEE or HEEH exists on either cluster.
- [ ] M4: training data. Ferrari et al. 2025: opening free energies of 5,778 domains, requested by the user through
  the authors' form, as for `diffusion_dynamics.md` Phase 2 E2. Tsuboyama et al. 2023 for global stability. Families
  are checked against M3's test set.
- [ ] M5 (large proteins): A3 of 5.6 as written.
  * The benchmark's homeodomain arms cannot stand in. Their ladder (0.87-1.00) is above the transition: native-start
    TM-score 0.19-0.31 at the bottom rung.
  * Ubiquitin under ff3.0 is the small-protein case of the same search failure: de novo never below 4.4 A in 7.2M tu,
    while the native start refolds.
- [ ] M6 (IDRs): Rg of IDPs with published SAXS Rg under ff2.1 and the base force field, at each experiment's
  temperature on the calibrated scale; the scaling exponent from intrachain distances; the Rg statistical
  inefficiency measured, so the mean is read over many decorrelation times. Running on the Mac since 10-28 with ff2.1
  and ff3.1 (4.2, B0); the base force field follows plan.md Phase 33.
- [ ] T1: the HDX term in the trainer, after M3 and the temperature mapping.
  * It must reach censored Spearman at least FF2's on M3's three proteins, with a block bootstrap
    (`diffusion_dynamics.md` Phase 2 E4's rule), and the native panels must not lose.
  * The Boltz add-on of that Phase 2 uses the same data and test, and both routes are read on one rule.

### 4.2 Disordered regions (2026-10-28)

User 10-28: "for upside3, investigate the possibility of enhancing intrinsic disordered region simulation accuracy";
then "record everything in upside3.md. We will implement this later." IDRs are third in 4.1's order. Literature:
`lit/review.md` (section 8), with the basis of every claim marked.

**What a disordered chain gets from the model.**
* **Training never sees one.** The 456 proteins are folded. NSE contrasts native-started free replicas with the
  restrained native. DSE (lambda 0.3) contrasts a near-Tm replica's unfolded frames (Rg above 0.67 Rg_coldest +
  0.33 Rg_hottest) with a SARW replica: Upside with the H-bond, burial, H-bond coverage and side-chain pair energies
  scaled to 0 (`training/ConDiv.py` `zero_for_sarw`), at the ladder's hottest temperature. The SARW's only excluded
  volume is the backbone's, and its only sequence dependence is the Ramachandran maps. Peng 2022 calls the DSE a
  regularizer giving "the minimum interaction energy that can still generate both a stable NSE and a disordered DSE".
* **The SARW is the high-denaturant end.** PNt is 51.3 A in water and 62.0 A in 4 M Gdn, which equals
  1.9 A x 334^0.6 (Riback 2017); FhuA's plug and reduced RNase A sit 12% and 7% below their 2 M Gdn sizes. The DSE acts
  only on near-Tm frames, so how far it shapes a disordered chain at 298 K is measured in B0, not assumed.
* **Electrostatics are implicit in the trained terms; there is no explicit electrostatic term and no ionic-strength
  input** (user 10-28). Charge effects are carried by the trained terms, among them the side-chain pair table to about
  7 A. **No explicit electrostatic term is added** (user 10-28): it would count the trained terms' electrostatics
  twice. So salt dependence is outside the model, and a charge-dependent error is fixed through the trained terms and
  the objective.
* **Local structure** comes from the NDRD coil-library maps (loop segments of folded proteins; per-map depths trained
  in ff3.1) and the H-bond energies, both fitted on folded proteins.

**Prior work with Upside.** No paper compares Upside with an IDP's SAXS or NMR data in water.
* Jumper 2018: thermally unfolded training proteins against 1.9 A N^0.6, 15% smaller at the heat-capacity peak and
  10% larger at high temperature.
* Riback 2019: the folding model by REMD at 280-320 K on PNt and five other unfolded proteins gives nu 0.4-0.6 with
  temperature and Rg proportional to Ree; PNt's measured Rg is not compared.
* Jeon et al. 2026 (PNAS 123:e2518125123): Upside2 REMD frames pooled over T 0.7-1.1 were the best single library for
  chemical-shift reweighting of 413 BMRB IDRs, and carried transient helices that CALVADOS2's lacked. Pooled over
  temperatures, this shows that Upside reaches the conformations, not that its weights at 298 K are right.

**What made other models right for IDPs** (`lit/review.md` section 2).
1. **The protein-solvent balance,** the first-order error in all-atom, Martini and residue-level models alike, and
   almost always over-compaction. In the a99SB-disp benchmark, force fields changed on the water side had mean IDP Rg
   errors of 9-13% against 22-36% for the rest, and a99SB-disp 6%. a03ws scales protein-water LJ by 1.10, Martini 3
   protein-water by 1.10. In Upside this balance sits in the burial, backbone desolvation and side-chain pair
   attraction.
2. **Local backbone statistics, which couple to collapse.** Martini3-IDP cut its Rg error from 1.06 to 0.39 nm through
   bonded terms alone. In a99SB-disp and ff19SB, helicity and compaction compensate each other. C36m removed spurious
   alpha_L.
3. **Residue specificity,** such as Arg against Lys and Tyr against Phe (HPS-Urry, Mpipi, CALVADOS's per-residue
   lambda).
4. Debye-Huckel at the experiment's ionic strength, in every residue-level model with good Rg: not open to Upside
   (above).

How folded proteins were protected:
* Residue-level models hold folded domains by restraints or rigid bodies. That is not open to Upside, which folds de
  novo with one energy.
* a99SB-disp's water change melted helices and folded proteins until the backbone carbonyl-amide LJ and the torsions
  were refit.
* SOP-IDP needs folded-protein energy scales 1.7-4.6 times its IDP ones.

So an IDR change is trained in the same objective as the native terms, and judged on the native panels as well.

Fitting to ensemble data has two precedents:
* FB-HPS (Dannenhoffer-Lafage and Best 2021) fitted 20 lambda by gradient descent on d<Rg>/dlambda =
  -(<Rg dU/dlambda> - <Rg><dU/dlambda>)/kT, with an L2 prior toward its start. That is the estimator of 4.1's IDR
  term.
* CALVADOS fits by Bayesian reweighting. It reuses frames while their effective fraction stays above 60%, under a
  prior over 70 hydrophobicity scales.

**B0, the baseline (4.1 M6).** Mac, `scratchpad/idp_20261028/`, launched 10-28 11:56.
* **Panel** (`panel.fasta`; table, sources and flags in `lit/review.md` 3.1-3.2): 16 SAXS entries of 24-441 residues,
  14 distinct sequences and conditions. They span net charge -43 to +15, kappa 0.08-0.78, proline 0-15% and aromatics
  0-15%:
  * Hst5, ACTR, Ash1;
  * p27 Cv14 and Cv78: same composition, 29.4 and 22.1 A;
  * PTMA;
  * A1-LCD* at 50 and 500 mM;
  * A1-LCD WT-NLS, -12F+12Y, WT+NLS and +7K+12D;
  * alpha-synuclein, CoRNID, PNt and tau 2N4R.

  Seven constructs are CALVADOS 2's curation, not checked against their primary papers.
* **Runs.** Constant T = T_K / 350.588 (the Peng calibration; mapping 278-298 K linearly is unproven), from a
  random-dihedral chain, pivot MC every 5 tu, dt 0.009, 30,000 tu per seed, t < 5,000 dropped. Stage A is ff2.1 x 6
  seeds and the SARW x 3; stage B is ff3.1 x 6. The two A1-LCD* entries share one set of runs.
* **Readout.**
  * `results.txt`: CA sqrt(<Rg^2>), the quantity a Guinier or MFF fit estimates, with the seed SE and both halves; nu
    from R(|i-j|) for |i-j| > 5; per-residue helicity (inside four or more consecutive alpha residues) and PPII.
  * `stat_ineff.txt`: the Rg autocorrelation time and n_eff.
  * CA Rg sits a few percent below an all-atom Rg, unmeasured here.
* **Pilot** (alpha-synuclein, 35.5 +- 1.0 A measured):
  * ff2.1: 35.6 and 31.0 A in two seeds, tau 340-1,300 tu, 2,500 tu block means 24-40 A.
  * SARW: 40.1 A, tau 110 tu.
* **Reading rules, fixed before the data.**
  * **Summary set.** The summary statistics use the 14 entries other than A1-LCD*. A1-LCD* is out of scope for salt
    and reported apart.
  * **Global offset:** the geometric-mean ratio Rg/Rg_exp, with a bootstrap over entries. It is called an offset if
    its 95% interval excludes 1.
  * **Sequence dependence beyond length:** Pearson r of log(Rg/N^0.6), model against experiment, and the rms log
    error after the offset is removed.
  * **Pairs,** against the measured ratios:
    * Cv78/Cv14: 0.753.
    * +7K+12D / WT+NLS: 1.131.
    * -12F+12Y / WT-NLS: 0.943.
    * WT+NLS / WT-NLS: 0.936.

    A pair is reproduced if its ratio lies on the measured side of 1 and within 2 combined SE.
  * **SARW:** its ratio places each force field between the training's DSE target and experiment.
  * **ACTR helicity:** read per region against the NMR helices. A number is compared only once the wild type's
    helicity is taken from Iesmantavicius 2014 Fig. 1.
  * **Convergence** follows findings 10.37: the seeds of one entry must agree. An entry whose seed spread exceeds its
    distance from experiment is reported as unresolved.
* **Routing.**
  * An offset with good sequence dependence points to the protein-solvent balance: I1, and I2 if the error leans
    toward the SARW.
  * Wrong pairs or a poor r point to residue specificity in the trained terms: I1 with variants in its training set.
  * Wrong helicity points to I3.

**B0 result: ff2.1, ff3.1 and the SARW** (stage A finished 10-28 16:10, stage B 19:37; `read_b0.py` applies the
rules above, output `read_b0.txt`).
* **Mechanism.** Under both force fields some disordered chains spend long stretches as collapsed globules, with Rg
  within 1.3 times 2.2 A N^0.38 (fraction of frames, ff2.1 / ff3.1):
  * ACTR, at 278 K: 77% / 68%.
  * The A1-LCD variants: 32-61% / 35-50%, except +7K+12D at 5% / 6%.
  * PNt: 16% / 23%. CoRNID: 2% / 18%. Alpha-synuclein: 8% / 5%.
  * p27, PTMA, Ash1 and tau: at most 1%.

  The switch between globule and coil is slow. For these entries the Rg autocorrelation time is 730-2,250 tu,
  against 50-180 tu for the SARW. So their seeds disagree: under ff2.1, PNt seeds span 23.9-51.6 A and
  alpha-synuclein seeds 20.0-33.7 A.
* **Offset.**
  * ff2.1: 0.886, 95% bootstrap 0.79-0.99.
  * ff3.1: 0.850, bootstrap 0.77-0.95.
  * SARW: 1.20, bootstrap 1.11-1.30.

  Both force fields have an offset by the rule, a compaction of 11-15%.
* **Sequence dependence beyond length.**
  * ff2.1: r 0.46, rms log error 0.21 after the offset.
  * ff3.1: r 0.49, rms log error 0.20.
  * SARW: r 0.34, rms log error 0.16.

  The misses whose seeds agree go opposite ways:
  * ACTR: 14.1 A (ff2.1) and 14.5 A (ff3.1) against 26.3 A.
  * p27 Cv78: 32.0 A and 30.9 A against 22.1 A.
* **ff3.1 against ff2.1.** Entry by entry, ff3.1's seed range overlaps ff2.1's, with two exceptions:
  * CoRNID: 31.3 A (seeds 25.0-40.0) against ff2.1's 47.7 A (43.6-50.2) and 47.0 A measured.
  * Hst5: 12.3 A against 13.0 A, measured 13.8 A.

  Apart from these two, ff3.1's trained map depths leave disordered-chain dimensions where ff2.1 had them.
* **Pairs.**
  * Cv78/Cv14: 1.056 +- 0.023 (ff2.1) and 1.076 +- 0.045 (ff3.1) against 0.753. Both are on the wrong side of 1:
    the charge-segregated variant stays as expanded as the mixed one.
  * The A1-LCD pairs carry model SE of 0.07-0.12 against measured differences of 6-13%, so this sampling cannot
    resolve them. -12F+12Y/WT-NLS falls within 2 SE under both force fields. WT+NLS/WT-NLS is on the wrong side
    under ff2.1 (1.081 +- 0.094) and within 2 SE under ff3.1 (0.900 +- 0.088). +7K+12D/WT+NLS is within 2 SE under
    ff2.1 (1.284 +- 0.120) and 2.2 SE high under ff3.1 (1.368 +- 0.107).
  * The SARW reproduces no pair, which is expected from a model without side-chain pair energies.
* **SARW placement.** For A1-LCD, PNt, alpha-synuclein and tau, both force fields sit below experiment and the SARW
  above it. For Cv78, both sit between them, and Cv14 is near experiment in both. Where the SARW matches experiment (Hst5, ACTR, Ash1, PTMA), both are
  smaller. So the force-field error does not lean toward the SARW.
* **ACTR helix 1** (1044-1054 in Iesmantavicius 2014's text numbering; their Fig. 1 labels the same residues one
  higher). The wild type's secondary structure propensity, read from their Fig. 1b, is 0.39 at 278 K.
  * ff2.1: 0.33 +- 0.10, with seeds 0.00-0.73.
  * ff3.1: 0.25 +- 0.08, with seeds 0.00-0.50.
  * SARW: 0.05.

  Each collapsed chain holds its own helix content, so helix 1 is unresolved under both force fields. Stage A's
  reading of "close to measured" came from the seed mean alone and is withdrawn. The other ACTR helices are not
  compared.
* **Convergence.** Unresolved by the rule (seed range above the miss): 11 of 16 ff2.1 entries (n_eff 53-580) and 10
  of 16 ff3.1 entries (n_eff 39-413). Every SARW entry is resolved except ACTR and Ash1, which sit on experiment.
* **Routing.** Both force fields make the same sequence-dependent error: ACTR and the A1-LCD variants collapse, and
  the charge-segregated Cv78 does not compact. By the rules this points to residue specificity in the trained terms,
  that is, I1 with variants such as A1-LCD and p27 in its training set. I2 is not indicated. I3 is undecided, since
  helix 1 is unresolved.
* **Open.**
  * ACTR is the only entry at 278 K (T 0.793), so its collapse may come from the temperature mapping rather than the
    sequence (step B1).
  * The A1-LCD pairs, ACTR's helicity and the unresolved entries need longer sampling or REMD before they can be
    read (step B2).
  * ACTR's 26.3 A is CALVADOS 2's value. The local copy of Kjaergaard 2010 holds only the abstract, so its SAXS
    source is unchecked (`lit/review.md` 3.2).

**Candidates.** Each changes what training fits or an existing trained term, never adds a term, and follows
section 4's rules.
* **I1. An IDR size term in ConDiv's contrast** (4.1).
  * The loss is sum_k (<Rg>_k - Rg_exp,k)^2 / sigma_k^2, with the gradient above from the per-frame dV/da that
    `compute_divergence` already returns.
  * Each training IDP needs its own sampling, since its Rg decorrelates over 10^2-10^3 tu and one step's 8,000 tu is
    too short. The sampling is reused by reweighting while the effective fraction holds (CALVADOS's 60%) and rerun
    when it falls.
  * An L2 prior holds the parameters toward the start (FB-HPS).
  * The training IDPs span composition, patterning and length, since FB-HPS found that a missing class biases the
    fit, and they are kept apart from B0's test entries.
  * The native terms stay in the same contrast.
  * Intrachain distances (PRE) can be added to the loss later.
* **I2. The DSE target.** The SARW pulls near-Tm unfolded frames toward high-denaturant dimensions. It is also the
  Sosnick lab's picture of an expanded denatured state with minimal residual structure, and it serves HDX (Peng 2022).
  B0 decides whether it is a suspect:
  * If ff2.1 sits at SAXS in water, the SARW is not hurting disordered chains.
  * If disordered chains come out expanded toward the SARW, the target is weighed against I1. A target from
    measurements in water is the alternative, but the DSE is an ensemble contrast and SAXS gives one number per
    protein, so that form is still to be designed.

  B0's answer: neither case holds. Both force fields sit on the compact side of experiment, away from the SARW, so
  the SARW target is not what puts disordered chains off.
* **I3. Local backbone balance.** Helicity of IDRs comes from the map depths, which ff3.1 trains, and the H-bond
  energies. B0's ff2.1 against ff3.1 is its first test. Its design is shared with 5.7's D3.
  * B0's answer: ff3.1's map depths leave dimensions where ff2.1 had them, except CoRNID and Hst5. ACTR helix 1 is
    unresolved, because the seeds of the collapsed chain range 0.00-0.73.
  * Targets: ACTR's helix-1 series, eight mutants at 23-69% helicity and the wild type at 0.39, at 278 K (5 C)
    (Iesmantavicius 2014). Reading it needs ACTR expanded and sampled first (B1, B2).
  * Per-residue helix and strand populations from BMRB chemical shifts, reweighted as in Jeon 2026.
* **Not taken:** an explicit electrostatic term (user 10-28, double counting); restraints on folded parts.

**Steps** (recorded 10-28; the user will take them up in a later session).
- [x] B0: the baseline above, read and recorded (10-28 11:56-19:37, MacBook Pro).
- [ ] B1. ACTR's temperature (Mac, about 30 min including the re-analysis). Run ACTR at 298 K under ff2.1 and ff3.1
  with 6 seeds each, and read it against the 278 K runs.
  * If ACTR opens at 298 K, the collapse comes from the T_K / 350.588 mapping, which reaches every IDP measured
    below 298 K.
  * If it stays collapsed, the cause is the sequence.

  The test compares the model with itself, so it needs no measured Rg at 298 K. The new entry is kept out of
  `read_b0.py`'s summary set.
- [ ] B2. Sampling, before I1, whose targets need the same sampling. The seed spread falls as the square root of
  run length, so 4x longer runs (120,000 tu) halve it.
  * Running the whole panel that long takes about 13 h per force field on the Mac. PNt would still be marginal, with
    a seed range of 23-28 A against a 12-14 A miss.
  * REMD on midway2 is the alternative for the long chains.
  * The choice is the user's.
- [ ] B0 on the base force field (plan.md Phase 33), once released.
- [ ] I1-I3: designs and the training IDP set, after B0 and in 4.1's order (after HDX and large proteins). B0
  points to I1.

**Resume here** (a cold session doing B1, B2 or a rerun):
* **Where the data are.** Everything is in `scratchpad/idp_20261028/` on this MacBook Pro only. It is gitignored,
  with 2.6 GB of trajectories. A session on another computer copies it first, or reruns B0, which takes about 7.7 h
  on 11 workers.
* **Environment.** Run `source .venv/bin/activate && source source.sh` from the repo root. `panel_run.py`,
  `stat_ineff.py` and `read_b0.py` need `UPSIDE_HOME`.
* **Adding an entry.**
  * Add it to `panel.fasta`; the header is `name|Rg_A|err_A|T_K|ionic_strength_M|reference`.
  * `panel_run.py` skips finished runs, so the rerun simulates only the new entry's seeds, then re-analyses every
    entry and rewrites `results.txt`.
  * Launch it detached with `scratchpad/detach.sh LOG python panel_run.py`. It takes 11 workers, and the Mac has
    10 performance cores.
  * Wait with one background waiter for the `DONE` line.
* **Longer runs.** Changing `DURATION` marks every run unfinished, and `run_idp.py` has no continuation, so every
  run restarts from a random coil. Longer runs therefore go in their own directory.
* **After any stage.** Run `python stat_ineff.py` and then `python read_b0.py`, and record the result in this
  section. `stat_ineff.py` covers only the models on `results.txt`'s `# seeds:` line, because files of a running
  stage are locked.

## 5. Large proteins (goal 4)

Investigation: plan.md Phase 22, 2026-10-19.

### 5.1 The problem

Upside folds large proteins de novo poorly (user, 10-19). Jumper 2018 (PLoS CB e1006578): all 20 test proteins under
100 residues come within 5 A of the native from an unfolded start, the CASP11 targets of 115-178 residues do not, and
T0855-D1 (115) never found its fold; ConDiv's short runs "optimiz[e] only over fluctuations in or near the native
well". Peng 2022 (FF2, JCTC 18:550): reversible folding only for the 40- and 43-residue mini-proteins; ubiquitin (76)
did not refold within a CPU-week because of long-lived misfolded states. Nothing above 150 residues has been trained
or benchmarked de novo: the 456 training proteins are 50-150 residues (median 114), the 44-domain panel 50-75, the
Peng de novo set at most 92.

### 5.2 What was measured

**The native-start size trend is the metric.** On the 89-domain mdCATH native-start panel, ff_2.1_ctrl38 folds 0.42,
0.44 and 0.24 of frames in the 50-100, 101-150 and 151-250 bins, "folded" being Q at or above the all-atom 5th
percentile, a cutoff that tightens with N as the all-atom Q spread narrows. Upside's mean Q sits about 0.13 below
all-atom in every bin and its Rg at 1.01 of native; with a fixed-offset cutoff the trend disappears. So Upside holds
large natives; the problem is reaching them.

**De novo test** (ff_2.1 released, 60,000 tu per run, Mac Studio). T4 lysozyme 2LZM (164 residues), RNase HI 2RN2
(155) and A. vinelandii flavodoxin II 1YOB (179, apo chain) all refold in vitro without a ribosome; the engrailed
homeodomain 1ENH (54) is the control. Starts: native, with the panel's anneal (T 0.05 to T over t 96-400); bulk, a
random-dihedral chain with the same anneal; vectorial, 30 residues annealed for 2,000 tu, then 5 residues added every
1,000 tu at T (each new segment PPII and turned away from the chain), then the full chain to 60,000 tu. Every last
frame was quenched for 1,000 tu at T 0.3. Q is the CA native-contact fraction (|i-j| >= 4, d < 8 A, formed below
1.2 d0; nonlocal Q only |i-j| > 12). H-bonds are counted from H and O placed on the backbone (H...O < 2.5 A), which
reproduces the crystal counts.

T 0.80, native x3, bulk x6, vectorial x6:

| | Q_end native / bulk / vectorial | nonlocal Q de novo, max | helical H-bonds de novo / native start | E_quench native / best of 12 de novo |
|---|---|---|---|---|
| 1ENH (54) | 0.74 / 0.43 / 0.40 | 0.20 | 0.64 | -271 / -264 |
| 1YOB (179) | 0.60 / 0.10 / 0.09 | 0.07 | 0.62 | -805 / -703 |
| 2LZM (164) | 0.62 / 0.15 / 0.15 | 0.11 | 0.36 | -793 / -684 |
| 2RN2 (155) | 0.68 / 0.11 / 0.11 | 0.09 | 0.56 | -736 / -645 |

T 0.88 and 0.95, native x1, bulk x3:

| T | 1ENH | 2LZM | 2RN2 |
|---|---|---|---|
| 0.88 | native melts (Q 0.64 to 0.17, Rg 1.9x native); coils <= 0.11 | native partly holds (0.48); coils 0.14-0.16 | native holds (0.66, nonlocal 0.62); coils 0.05-0.08 |
| 0.95 | all unfolded | all unfolded | native melting (0.26); coils <= 0.05 |

No coil forms native tertiary structure at any of the three temperatures. At T 0.80 every de novo run ends in a
compact, helix-poor globule (Rg 0.9-1.3 of native for the large proteins, 1.1-1.7 for the control); T4 lysozyme's
carry 1.5-5x the crystal's nonlocal, sheet-type H-bonds. Native starts keep crystal helix and H-bond content; their Q
drifts by up to about 0.2 over 60,000 tu and mostly plateaus.

**Energy by term, quenched globule minus quenched native,** bulk / vectorial, each the mean over the last 10 frames of
every quench (6, 6 and 3 runs), in E_up. The terms sum to the engine total.

| term | 1ENH | 1YOB | 2LZM | 2RN2 |
|---|---|---|---|---|
| H-bond | +34 / +50 | +115 / +118 | +140 / +142 | +160 / +144 |
| rotamer (side-chain free energy) | +23 / +32 | +35 / +39 | +42 / +23 | +24 / +37 |
| side-chain burial | +2 / +3 | +20 / +20 | +13 / +15 | +13 / +15 |
| Rama | -5 / -5 | -7 / -19 | -20 / -27 | -19 / -32 |
| backbone desolvation | -4 / -7 | -17 / -17 | -18 / -19 | -16 / -15 |
| total (bulk / vectorial) | +49 / +73 | +155 / +150 | +179 / +159 | +163 / +149 |

Most of the gap is backbone H-bonds the globule never forms. The globule's side chains recover 0.87 (1ENH), 0.91, 0.90
and 0.93 of the native's side-chain energy (rotamer plus burial, bulk), so collapse costs little. Rama and backbone
desolvation favour the globule, by 24-47 E_up in the large proteins.

**The co-translational hypothesis** (the user's: folding is co-translational, S. Fried, JHU: about a third of E. coli
proteins do not refold in vitro, so a thermodynamic force field cannot fold large proteins from a denatured chain)
does not explain these failures. N-to-C growth equals bulk collapse in every measure (slightly worse on the
homeodomain), and all three large proteins refold in vitro without a ribosome. It may hold for obligately
co-translational folders. The growth here has no ribosome tunnel or surface, and its rate (200 tu per residue) is not
mapped to real time.

**Not shown by these runs:** the globules' free energy against the native's (only energies were compared); a size
contrast under one protocol (the control does not fold at constant T in 60,000 tu at any of the three temperatures,
so the size dependence rests on Jumper 2018 and the user's observation); generality (four proteins, one force field,
constant temperature).

### 5.3 Reading

For these in-vitro refolders, Upside's failure is search. RNase HI's native is stable at T 0.88, and its coils do not
find it at 0.80 or 0.88. Upside ranks the native 150-180 E_up below the globules on average, and 90-110 below the best
of them. Yet the chain collapses on side-chain burial, which already pays 87-93% of the native's side-chain energy;
inside the globule it cannot form native H-bonds, and Rama and backbone desolvation reward the non-H-bonded compact
backbone.

Mechanisms, ranked by the evidence (literature in `lit/review.md`, code in `code_audit/audit.md`):

1. **No negative design against compact misfolds.** ConDiv contrasts a native-restrained replica with 12 free
   replicas started at the native (8,000 tu, T 0.8-1.1) and an expanded self-avoiding-walk ensemble. Compact
   non-native globules never enter the contrast, and their number grows with N (Bryngelson and Wolynes 1987).
   AWSEM (T_f/T_g, Davtyan 2012) and OPEP (energy gap, Sterpone 2014) handle this with decoys.
2. **The training range and weighting.** Training stops at 150 residues. Each protein's contrast is summed, not
   normalised per residue (`ConDiv.py:403-405`), so 110-150-residue proteins dominate. The fixed 12-replica ladder
   loses acceptance with N (0.63 / 0.55 / 0.49 across the training bins, 0.40 predicted at 200). Airas and Zhang 2025
   found that adding larger structures to training is what fixed size transfer.
3. **Backbone desolvation acts as a surface term.** It rewards exposed NH and CO (scale -0.30). The native's exposed
   fraction falls from 0.27 (50-100 residues) to 0.19 (251-400), while almost every group of an unfolded chain is
   exposed, so its pull toward unfolded grows with N. Here it also favours the globules, by 15-19 E_up.
4. **Weak cooperativity.** Jumper 2018 reports a transition "not fully cooperative"; additive pair terms fail the
   calorimetric two-state test (Chan 2000).

Ruled out, or not size-dependent: explicit size terms or neighbour-list capacities (none exist); BP convergence
(14-28 iterations from 76 to 271 residues); the burial window (the same coverage range at every size, residues above
it rising from 5.3% to 10%); the pivot sampler (it loses efficiency with N but cannot bias).

### 5.4 Engine defect fixed on 2026-10-19 (findings 3.13)

`environment_coverage_hb` gave each covering side chain the type of the H/O site whose index equals its residue
number, so only 7-9% got their own type and each NH/CO's burial summed scrambled weights (0.56 ASP to 2.12 TRP). The
same code is in master; it does not depend on size. On ff_2.1 the effect is small (backbone desolvation -0.1 to
+0.2 E_up at the natives, +0.05 to +1.3 at the globules, no other term moving; native-start dynamics within
single-run noise), so it is not a cause of the large-protein failure. ff_2.1 and every training so far ran with it;
Upside 3 trains with it fixed. Fix and verification: findings 3.13, `envfix/` (section 8).

### 5.5 Not planned: a co-translational protocol

It did not help here, and the large test proteins fold in vitro without one.

### 5.6 Plan

**Phase A: baseline on the chosen side chain and the audited engine.**

1. [ ] Deploy the audited engine to the Upside 3 training tree, a tree of its own: a fresh midway2 build that passes
   the audit harnesses and is bitwise equal to the locally tested build on their fixtures, installed by rename,
   carrying file modes (findings 10.11).
2. [ ] Fix the de novo benchmark: 8-12 single-domain proteins of 50-250 residues with published in-vitro refolding,
   over the alpha, beta and alpha/beta classes, the four here included; scored by nonlocal Q and a size-normalised
   RMSD (RMSD_100) or TM-score, from coil, with replica exchange.
3. [ ] Replica exchange on RNase HI and the homeodomain: 16 replicas over T 0.80-1.00, from coil and native starts
   (about 5 h and 15 h on the Mac, or about 27 h on midway2). It answers 5.2's open questions: if the globules beat
   the native in free energy, the energy ranking must change; if the native wins but coils never reach it, the funnel
   shape is the target.
4. [ ] **Repeat 5.2 and 5.7 with ff3.1, once its validation is in** (plan.md Phase 30 item 4; user 10-27). 5.2 and 5.7
   are ff_2.1's, so no 5.7 change is trained before this. The same four proteins and protocol (`sim/landscape.py`:
   native x3, bulk x6, vectorial x6 at T 0.80, 60,000 tu, then quench), on the audited engine with configs from the
   current writer; the 5.2 energy-by-term table, then the six 5.7 scripts, the engine-reading two with the current
   `py/upside_engine.py` in place of `eng_1019`. Each 5.7 row is read against ff_2.1's; a design candidate keeps its
   place in 5.7's order only if its target still holds.

**Phase B: training and force-field changes, one at a time, each scored on A2 and on the native panels (44-domain,
mdCATH).** User 10-26: the fix must include a change in the force field's design, not only in its training. B1 and B2
combine with any change of B3.

1. [ ] **Compact misfolds as negative examples, through a replay buffer.** Each training protein keeps a buffer of
   its own states: globules collapsed from coil under the current force field (the bulk protocol at T 0.80, then
   quench) and its free replicas' end states from earlier visits. Some of the 12 free replicas start from the buffer
   instead of the native, and the buffer refreshes from their end states, so decoys are made once, not every epoch.
   This keeps ConDiv's objective and reduces the short-run bias Jumper 2018 names (free replicas started at the native
   explore only its well); it is the standard remedy in energy-based model training (persistent CD, Tieleman 2008;
   replay buffer, Du and Mordatch 2019; Nijkamp et al. 2020). The share of buffer starts is the knob, set on held-out
   proteins. The alternative form is an explicit native-against-decoy gap term next to the SARW term. Targets: the
   native-globule gap per residue; the globule's side-chain energy fraction (0.87-0.93); Rama plus desolvation
   favouring globules (24-47 E_up). Hypothesis, not tested: findings 1.37's loss (training buys nothing on its own
   proteins while held-out folding falls) may come from the same bias, since a native-started 8,000 tu run carries
   little folding information.
2. [ ] **Train past 150 residues, on crops in a frozen context,** the physical version of AlphaFold2's 256-residue
   crops fine-tuned at 384 (AlphaFold3 adds spatial crops), so an example's cost does not grow with the chain.
   * Add mdCATH domains of 150-400 residues, kept apart from the panel's held-out domains. In each, a crop of 60-90
     residues moves, contiguous or a compact subdomain (a foldon); the rest is frozen at the native through
     `/input/fix_rigid/atom_indices` (`src/martini_fix_rigid.cpp:316`), which any config can carry. The
     native-restrained replica restrains the crop as now.
   * The crop sees large-protein cores and subdomain interfaces (the share of residues above the burial window rises
     from 5.3% to 10% with size, `code_audit/audit.md` 1a), and on the hotter rungs it can partly unfold and repack,
     so the trainer also sees local misfolds in a native context.
   * To be measured: how far a crop moves in 8,000 tu (the 54-residue control does not refold in 60,000 tu), and
     whether replica acceptance follows the crop's size. Each step still costs a whole-protein evaluation, but runs
     need not get longer. The SARW replica needs a crop form, or is skipped for crops.
   * Normalise each protein's contrast per residue (`training/ConDiv.py:396`), and space each ladder by 1/sqrt(N) of
     its moving residues.
3. [ ] **Force-field changes (5.7),** each trained as the only change against a control, in the order A3 sets (5.7,
   Order).

**Phase C: validation.** The native panels must not lose, and the de novo benchmark must gain. Re-check the hybrid
(glpG), which runs on the protein force field. Calibrate Upside 3's temperature scale (3.3).

### 5.7 The force field term by term (2026-10-26)

User 10-26: review every term and every possible design change, not one term. Each term is read on 5.2's T 0.80
quenched end states; "globule" means the bulk and vectorial ones, and ranges are over the three large proteins.
Scripts and their `*_T080.txt` outputs are in `scratchpad/large_protein_20261019/analysis/`.

**Every number in this section is ff_2.1's,** run on the 10-19 engine, before the audit's fixes, from configs written
before findings 3.13's fix. They are repeated with ff3.1 on the audited engine at Phase A4 (user 10-27), and the
candidates are re-ranked on those numbers before any is trained. ff3.1 differs from ff2.1 only in its Ramachandran
library (the other five files are ff2.1's), so the Ramachandran row and D3 are the most likely to move; the engine
fixes can move any row.

**What the trap is.** The globule is 150-180 E_up above the native, so it wins on entropy if at all; it is stable
against the unfolded chain, and the chain does not leave it. A design change can **widen the native's lead** (needed
if A3 finds the globules ahead in free energy) or **make non-native compact states less stable relative to the
unfolded chain while the native keeps its stability** (needed if the native wins but coils do not reach it). The
trap's depth is the side-chain energy a globule collects without native structure.

**Every term, globule against native.** The term gaps are 5.2's; everything else was measured on 10-26.

| term (node) | form | globule minus native | what the end states show |
|---|---|---|---|
| H-bond (`hbond_energy`) | one energy per bond by its residue's Rama class (-1.96, -1.95, -1.77), plus a bias on the per-peptide count; no dependence on burial or on neighbouring bonds | +115 to +160 | the globules form 0.36-0.62 of the native start's helical bonds; 2LZM's form 1.5-5x the crystal's nonlocal, sheet-type bonds |
| side-chain free energy (`rotamer`, with its pair table, rotamer 1-body, `hbond_coverage` and `hbond_coverage_hydrophobe`) | BP over one-bead rotamers; pair splines in r and two angles, to about 7 A | +23 to +42 | with burial, the bulk globules keep 0.90-0.93 of the native's side-chain energy, at 0.54-0.72 of its hydrophobic/polar segregation (CB counts within 8 A) |
| side-chain burial (`sigmoid_coupling_environment`) | per residue, a sigmoid of the type-weighted side-chain count within 9 A; 20 weights shared by every centre type | +13 to +20 | all from hydrophobic residues, whose coverage is 3.4-3.6 against 4.1-4.7; polar residues' burial energy is equal within 3.1 E_up; the hydrophobic share of a hydrophobic residue's coverage is 0.57-0.62 against 0.61-0.70 |
| backbone desolvation (`bb_sigmoid_coupling_environment`) | per NH/CO, -0.30 x sigmoid(side-chain coverage + 2 x pairing - 2) | -15 to -19 | the globules have 0.33-0.48 more unpaired groups per residue, mostly exposed; buried and paired groups are 0.41-0.47 per residue against 0.74-0.80; buried and unpaired 0.17-0.21 against 0.09-0.16 |
| Ramachandran (`rama_map_pot`, `rama_map_pot_ref`, sheet mixing) | per residue, a neighbour-dependent coil-library map at kT_ref, with trained basin offsets (ff3.1) | -7 to -32 | native-helical residues spend 0.44-0.49 of globule frames in the strand/PPII basin; by type the order follows Chou-Fasman (Spearman +0.64 over 15 types, p 0.01; ILE, VAL 0.62-0.64), but even GLU, ALA and LYS spend 0.34-0.49 |
| backbone sterics, covalent springs | walls between N, CA, C and CB; bond, angle and omega springs | small | not a design question |
| absent: electrostatics beyond the pair table | | | Debye-Huckel between charged CBs beyond 7 A (eps_r 80, Debye length 10 A): native and globule within 1 E_up, inside the globules' own spread of 0.8-2.7 |
| absent: burial beyond 9 A | | | hydrophobic/polar segregation at 14 A: the bulk globules keep 0.56-0.72 of the native's, as at 8 A (0.59-0.72); vectorial ones 0.48-0.55 against 0.54-0.61 |

Kinetics, from the bulk runs: the chain does not collapse before forming helices. For the first 2,000-5,000 tu it
stays at 1.5-2.0x the native Rg with 0.12-0.40 of its helical bonds; it then collapses, to 1.1-1.3x by 20,000 tu,
with 0.32-0.54.

**What the forms cannot express.**
* Burial and pairing are interchangeable in backbone desolvation: a buried H-bond gains 0.3 E_up over an exposed one,
  and a buried unpaired group pays nothing beyond the lost exposure reward.
* No side-chain energy depends on how dense a contact's surroundings are; the burial term counts side chains, not
  which ones.
* No term couples neighbouring residues' conformations. The maps condition only on the neighbours' identity, so helix
  nucleation and propagation come only from the H-bond geometry.
* No energy depends on temperature except through the kT_ref terms (3.1); the hydrophobic effect's heat capacity,
  which grows with N, is absent.

**Design candidates.** Each is a spline table or the parameters of one, with exact parameter derivatives through the
audit harness before any training (section 4).

1. **D1. Backbone polar groups: one table of burial and pairing.** A 2D spline E(side-chain coverage, pairing) per
   NH/CO replaces backbone desolvation's sigmoid of a sum, giving the H-bond a dependence on burial and setting the
   four corners (exposed or buried, paired or unpaired) separately. Physics: backbone H-bonds in nonpolar
   surroundings are worth more than exposed ones (amide-to-ester mutations: Deechongkit et al. 2004, Nature 430:101;
   Gao et al. 2009, Nat Struct Mol Biol 16:684). Target: buried paired groups (0.74-0.80 per residue in the native
   against 0.41-0.47), the largest structural contrast measured; the buried-unpaired corner alone differs by at most 0.11 per residue (0.01-0.02 in 1YOB). The same
   table answers 5.3's third mechanism and is the soluble side of the hybrid's TM1 opening (findings 4.7). It mainly
   widens the native's lead, but also rewards the globules' non-native ladders, so it is judged on 2LZM.
2. **D2. Contacts that depend on their surroundings,** the force-field form of AlphaFold2's triangle updates (a
   pair's energy depends on the third bodies around it), established as AWSEM's density-switched contacts (Papoian et
   al. 2004, PNAS 101:3352; Davtyan et al. 2012). In Upside, each side-chain pair energy is scaled by a spline of the
   two residues' `environment_coverage_sc`, which uses the rotamer prior, not the BP solution (`code_audit/audit.md`
   1a), so it is fixed before BP runs and adds no circular dependence. Target: the 0.90-0.93 side-chain energy at
   0.54-0.72 segregation, the globules' hydrophobic contacts at coverage 3.4-3.6 against 4.1-4.7. It makes the trap shallower, because a
   contact made at core density counts more than one made in a loose globule.
3. **D2b. Who covers a residue.** The burial node already reads 400 weights (`number_independent_weights` 400);
   ff_2.1 trains 20 shared ones, and untying them by centre class needs no engine change. Backbone desolvation keeps
   its own 20 weights, since its writer refuses unequal rows. Target: modest, the hydrophobic share of 0.57-0.62
   against 0.61-0.70.
4. **D3. Secondary-structure balance in compact states.** The order by type is right, but the offset toward
   strand/PPII is not. Either a nearest-neighbour term on consecutive residues' helix scores (the `hbond_energy` basin
   sigmoids), a Lifson-Roig-like nucleation and propagation no current term can express, or per-type alpha/beta
   energies through `hbond.h5`'s residue-class offsets, which the engine already reads (2.2 of `up.md`) and no file
   uses. Before training it needs a reference for the right offset: the helicity of these proteins' isolated helix
   peptides, or all-atom compact states (mdCATH high-temperature frames).
5. **D4. Temperature-dependent burial.** A heat-capacity coefficient per type on the burial and pair energies. It
   bears on the size dependence of T_m and on replica ladders (lit/review.md 1), less on the trap at fixed T, and
   belongs with the temperature calibration (3.3 item 5); there is no measurement here.

**Not taken, by measurement:** a penalty for buried unpaired NH/CO alone (0.01-0.11 groups per residue); burial beyond
9 A (no better at 14 A than at 8 A); long-range electrostatics (within 1 E_up beyond 7 A); a desolvation barrier in
the pair table, meant to slow collapse, which does not outrun helix formation here.

**Not testable on end states, and large changes:** more beads per side chain, for packing finer than one bead; a
backbone frame-frame orientation term (AlphaFold2's residue frames), where the side-chain beads' angular pair term
already carries residue orientation and 2LZM's non-native ladders are the only hint.

**Not taken: structure priors.** MSA or language-model contacts keep AlphaFold2 accurate on large proteins, but would
make Upside a structure predictor and bias its pathways. Diffusion samplers (Boltz) as Monte Carlo moves have no
proposal density for the acceptance test, and import a structure prior too.

**Order.**
* First, A4 re-measures this section with ff3.1; the rules below apply to its numbers, not to ff_2.1's.
* If A3 finds the globules ahead in free energy: D1 first (the largest contrast), then D2b.
* If the native wins but coils do not reach it: D2 first (the trap's depth), then D3.
* Either way, D4 goes with the temperature calibration.
* The changes are made one at a time, each against the control, and scored as in Phase B. Every training so far has
  lost 0.13-0.23 folded fraction on leaving ff2.1 (findings 1.37); each change is judged against that loss.

## 6. GPU (goal 5)

### 6.1 Measured (`gpu_report.md`)

Per force evaluation, ff_2.1 native starts, one Apple M3 Max core, the 10-19 source; 3 evaluations per step.

| term | 1ubq 76 | 2qke 108 | glpG 210 | 1AO6 578 |
|---|---|---|---|---|
| BP solve | 23.9% | 21.8% | 20.3% | 21.6% |
| side-chain pair edges and edge table | 26.8% | 26.2% | 25.7% | 25.3% |
| burial and environment | 18.9% | 19.5% | 21.5% | 20.7% |
| H-bond | 13.2% | 14.5% | 17.6% | 19.0% |
| side-chain placement | 6.4% | 6.8% | 5.9% | 4.7% |
| Rama, sterics, springs | 4.5% | 4.8% | 4.4% | 3.4% |
| integrator and thermostat | 1.5% | 1.5% | 1.4% | 1.1% |
| untimed (rotamer derivative scatter) | 4.8% | 4.7% | 3.1% | 4.1% |
| ms per evaluation | 0.262 | 0.338 | 0.750 | 2.63 |

* The protein model has no single hot spot, so a speed-up needs nearly every node on the GPU. The hybrid is different:
  in 1rkl in DOPC (3,812 particles) the dry-MARTINI pair loop is 86% of 17.5 ms per step, and every Upside protein
  node together 1.5%.
* One system cannot fill a GPU: an evaluation is about 150 kernel launches at 3-5 us each, about the CPU's whole
  0.75 ms for glpG. So the GPU must run many systems per launch, a replica-exchange ladder or a ConDiv minibatch
  (24 proteins x 14 replicas = 336 systems).
* Already suited: the graph is static, all arithmetic is float32, one OpenMP thread runs per system, and the RNG is
  counter-based Threefry, which runs under CUDA. Not suited: BP's data-dependent iteration count, its per-node
  products of edge messages, and per-step allocations and map lookups in several nodes.
* The trainer's parameter derivatives cost about 1/1000 of the simulation, so they can stay on the CPU.

### 6.2 One force field on two backends

The force field already runs on numerically different backends without retraining: the Mac and cluster builds differ
in their reciprocal square-root refinement (and differed in the finite-math flag until A28), and between el7 and el8
forces are bitwise equal but trajectories differ from step 6.

A GPU is a third backend under the same standard. It uses the same spline coefficient arrays, fitted on the host and
uploaded, never refitted; BP's same initialisation, damping, tolerance, stop test and synchronous schedule, each system
stopping at its own test (another schedule can reach another fixed point, a model change); float32 kernels with
deterministic fixed-point force accumulation; and the same Threefry streams. The tolerance has a precedent: the user
ruled that the BP stop-test fix needs no retraining, and it changed static-frame energies by a median of 7.5e-4 E_up,
forces by 8e-5 of their RMS, the ConDiv gradient by 0.1% and REMD T_mid by at most 0.002. The GPU should sit well
inside that.

Verification ladder:

* **L0.** Every parameter array, table, index and RNG key on the device is bitwise equal to the host's.
* **L1. Per node, on static frames** (about 20 proteins x 10 frames, including misfolded globules and hybrid frames):
  energy within 1e-5 relative; force RMS within 1e-5; parameter derivatives within 1e-4; BP compared at a tight
  tolerance and at the production one.
* **L2. Short dynamics:** KE/(1.5 kT), drift, C-N bond and angle distributions, per-term means; lipid diffusion and
  area per lipid in the hybrid.
* **L3. Ensembles:** REMD melting on the Peng set (the `bp_validation` harness), Q and Rg distributions, the mdCATH
  panel, glpG TM4.
* **L4. Training:** a ConDiv gradient from GPU trajectories within that step's own sampling noise.

### 6.3 Route

| route | effort (estimate) | risk |
|---|---|---|
| a. a CUDA copy beside the CPU code | 4-8 person-months | high: two implementations drift apart |
| b. one implementation of each potential, compiled for CPU and GPU | 6-10 person-months | medium: CPU parity has to be re-proved |
| c. many systems per kernel launch | 2-3 more person-months | medium; needed by any route for a speed-up |
| d. OpenMM plugin, or a JAX/PyTorch rewrite | 4-12 person-months | high: OpenMM has no BP, no derived coordinates, no batching, and refits its own splines; a rewrite is a second implementation |
| d2. GPU for the dry-MARTINI pair loop only | 1-2 person-months | low; the hybrid only, at most about 7x |

**Recommended: b with c.** `InteractionGraph<IType>` is already an evaluator pattern; templating the evaluators on the
scalar type gives one implementation for both backends and removes today's split between the SIMD `compute_edge` and
the scalar parameter derivative.

1. [ ] G1, after the audit closes: refactor the CPU engine into the single-source form, bitwise equal to the audited
   CPU engine, with the audit harnesses as its test.
2. [ ] G2, prototype: side-chain pair edges, environment coverage and H-bond coverage on the GPU, batched. Go/no-go:
   one GPU evaluates a batch of 336 ConDiv-sized systems in about 2 ms, which matches 336 broadwl cores.
3. [ ] G3, the remaining nodes, after any new terms (section 4) are fixed.
4. [ ] G4, the BP kernels, after the verdict.
5. [ ] G5, the ladder L0-L4.

The CPU engine stays the reference, and development happens on the Mac against the host backend (the Mac's GPU cannot
run CUDA). d2 is a separate option, worth taking only if the hybrid is in Upside 3's scope.

**The side-chain verdict changes the size of the work.** With BP, the rotamer node is about half the protein time and
holds the three hardest kernels; it is also the only node whose CPU-GPU difference can exceed rounding, because
rounding can move the iteration at which the stop test fires. A deterministic side chain is an ordinary pair term,
saving about a quarter to a third of the GPU effort. That is no reason to choose it: it folds worse (0.264 against
0.365 at matched length, 0.241 converged; findings 1.21).

### 6.4 Hardware

* **midway2, measured.** Our account can use `gpu2`: 5 nodes of 4 Tesla K80; QOS of 4 GPUs, 36 h, 10 jobs; billing
  counts CPUs only. The K80 is Kepler, which CUDA 12 dropped, so it needs CUDA 11.x (`cuda/11.5` is the default
  module). The V100 node has been down since 2024-06. No Kokkos or OpenMM modules. Fit for correctness testing, not
  for speed.
* **midway3, from RCC documentation, not measured.** The `gpu` partition has 1 node of 4 A100, 5 of 4 V100 and 5 of
  4 RTX6000.
* **beagle3, from RCC documentation.** 22 nodes of 4 A100 and 22 of 4 A40. RCC's documentation says beagle3 uses no
  SUs, while CLAUDE.md says it bills the group's GPU allocation. Ask before any use.

## 7. Decisions for the user

Answered 10-19: Upside 3 is a new release, with no fix measured against ff2.1 (section 1); martini-dev is Upside 3,
fixes go only there, and master is not changed; ff2.1 keeps running for development and Upside 3 later refuses it;
the audit starts at once (plan.md Phase 24). So the restraint-pairing defect needed no ruling of its own: it was fixed
in Upside 3 like every other confirmed defect (A25).

Still open:

1. Units: A+ (recommended) or C. Also, should `martini.h5` move to native units, with conversion at config time? That
   would change the CLAUDE.md unit contract.
2. GPU: the route (b with c, prototype first); which GPUs (midway3 `gpu` under which account; beagle3 once its
   billing is clear); whether the hybrid's d2 is in scope.
3. Goal 4: B1's form (5.6: the replay buffer proposed, the gap term the alternative); which 5.7 change goes first is
   set by A3.
4. Items in the running Upside 2 campaign that await a NO GUARDS ruling: the `run_remd.py` rollback; the HDX analysis
   items of 2.1; three checks in the hybrid code.
5. Training terms (4.1): the base force field for M3 and T1; the HDX temperature mapping, decided with item 1; the
   request for Ferrari et al.'s data, which only the user can make.
6. Disordered regions (4.2). B0 points to I1, with variants in its training set. It does not point to I2, and I3
   is undecided. Still open: whether to design I1 and when, and whether to first spend Mac time on ACTR at 298 K
   and on longer sampling of the unresolved entries.

## 8. Files

`scratchpad/upside3_planning_20261019/` (this MacBook Pro):

| path | what |
|---|---|
| `audit_report.md`, `audit_work/` | defect inventory, audit units, writer/reader registry (`registry.txt`), warning, cppcheck and sanitizer logs, the parameter-derivative prototype (`harness/param_fd2.py`), the finite-math probe (`probe/fm.cpp`) |
| `units_report.md`, `units_work/` | every unit site with file:line, the parameter-file table, options and verification |
| `gpu_report.md`, `gpu_work/` | profile scripts and logs, replica scaling, the hybrid profile, `web_notes.md` with the sources of every documentation claim |

`scratchpad/hdx_term_20261028/` (this MacBook Pro; 4.1 M1):

| path | what |
|---|---|
| `fetch.sh`, `run.sh`, `build_inputs.sh`, `meta.py`, `flags.sh`, `analyze.py`, `fig.py` | copy of the arms from midway2; driver; HDX topologies from `1.config.py`; per-frame energy, walker and CA-RMSD; the stock estimator per replica; MBAR readout with its rules in the docstring; figure |
| `results/results.txt`, `results/ubq_dG_measured.png`, `results/summary.npz` | the readout of findings 5.7 |
| `data/` (4.4 G), `results/*__top_*` (2.0 G) | the three arms (`/beagle3/trsosnic/yinhan/ff3_benchmark/runs/`), copied by `fetch.sh`; flags and raw estimator arrays, regenerated by `run.sh` |

Copies of the scripts and the readout: `/project/trsosnic/yinhan/checks/handoff_20261026/scratchpad/hdx_term_20261028/`.

`scratchpad/idp_20261028/` (this MacBook Pro; 4.2):

| path | what |
|---|---|
| `lit/review.md`, `lit/panel.fasta`, `lit/make_panel.py`, `lit/dl*/` | literature review with the basis of each claim; the panel's sequences and SAXS values; its builder (`python make_panel.py lit`); the full texts, the CALVADOS 2 tables and Bremer 2022's SI it read |
| `run_idp.py`, `panel_run.py`, `stat_ineff.py`, `read_b0.py` | one run from a sequence (ff directory or the ConDiv SARW); the panel driver with its analysis (two stages); the Rg statistical inefficiency of the models `results.txt` has analysed; 4.2's reading rules |
| `panel.fasta`, `results.txt`, `results.npz`, `stat_ineff.txt`, `read_b0.txt`, `panel.log` | the panel as run; B0's readout and per-residue arrays; the reading by the rules; the driver's log |
| `runs/`, `pilot/` | trajectories (one directory per sequence); the alpha-synuclein pilot (`conv_pilot.py`, `pilot.txt`) |

`scratchpad/large_protein_20261019/` (the Mac Studio):

| path | what |
|---|---|
| `sim/landscape.py`, `sim/common.py` | native, bulk and vectorial runs with quench |
| `sim/analyze_landscape.py`; `analysis/{ss_check,hb_check,terms}.py` | Q, Rg, quench; helix, local and nonlocal Q; H-bonds; energy by term |
| `analysis_T0{80,88,95}.txt`, `analysis/*_T0{80,88,95}.txt` | outputs |
| `runs_T080`, `runs_T088`, `runs_T095` | trajectories |
| `envfix/` | the 3.13 fix's builds (pre and fixed), `coverage_test.py`, `compare.py`, `fd_test.py`, `dyn_check.py`, `input_diff.py` |
| `lit/review.md`, `code_audit/audit.md` | literature review, code audit |
| `analysis/{buried_unpaired,burial_composition,longrange_burial,electrostatics,helix_to_strand,collapse_vs_helix}.py`, `analysis/*_T080.txt` | the term-by-term review (5.7); the first two need `analysis/eng_1019/`, the engine binding that matches the 10-19 `obj/` build |
