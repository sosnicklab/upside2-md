# Upside 3

**Status: the engine audit (goal 1) is under way, plan.md Phase 24, ledger A1-A102 fixed and verified (`findings.md` 3.14); the other goals are proposals.** The user set five goals on 2026-10-19 (section 1). The protein-model work starts
after the side-chain verdict, deterministic (no-BP) against probabilistic (BP) side chains at their converged epochs
(plan.md Phase 21 item 7, findings 1.21), and is built on the model that wins. The measurement is in (10-20, findings
1.21): at both fixed points BP folds 0.353 against no-BP's 0.241 on 89 held-out mdCATH domains (p 1e-5); the choice
is the user's. The engine audit and the units work do not depend on that verdict. Everything measured below is
ff_2.1, the probabilistic side chain, unless named.

Records: plan.md Phases 22-24. The planning reports for goals 1, 2 and 5, with their scripts and logs, are in
`scratchpad/upside3_planning_20261019/` (`audit_report.md`, `units_report.md`, `gpu_report.md`; gitignored). The
large-protein investigation is in `scratchpad/large_protein_20261019/` on the Mac Studio (76 G with trajectories).

## 1. Goals and order

The user's goals:

1. **Debug the whole repo.** Find every bug ourselves, as Jumper did with Claude on master. Audit this branch, master's
   code first (section 2).
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

Why this order:

* Every Upside 3 training has to run on the audited engine. Each defect found so far was carried by trainings that ran
  before it was found (section 2.1), and a training is weeks of cluster time.
* The GPU port is written once, against the audited engine in its final units and with its final set of terms.
* The audit is read-only and local, and it is the longest item (section 2.2). Starting it before the verdict is the
  user's call.

**Upside 3 is a new release, not a patch series on ff2.1** (user, 10-19).
* No list of fixed bugs and their effects on ff2.1 will be published, so no fix is measured against ff2.1 or master
  for attribution.
* Each fix is still proven correct (section 2.2).
* Upside 3 as a whole is validated in section 5.6, Phase C.
* Master and ff2.1 stay as released.

## 2. Engine audit (goal 1)

### 2.1 What is known

Fixed defects that trained force fields carry:

* the BP stop test (`rotamer.cpp:279`, missing `fabs`; Jumper's report, fixed 2026-09-26): ff2.1 and ff_3.0;
* findings 3.12, the pair term's angular parameter gradient (`interaction_graph.h`): every ConDiv run before 10-18;
* findings 3.13, `environment_coverage_hb`'s covering-residue types (`environment.cpp`, `upside_config.py`): ff2.1
  and every training so far (section 5.4).

Open. Each is fixed in Upside 3 once confirmed. Master and the running Upside 2 campaign keep them until the user
rules otherwise (section 7):

* `make_restraint_group` (`py/advanced_config.py:560`): `x/3 != y/3` is float division in Python 3. It skips only
  self-pairs, where it means to skip same-residue pairs. It is the same in master and affects only ConDiv's
  native-restrained replica.
* `run_remd.py` rolls NaN replicas back and continues, a guard under the NO GUARDS rule (plan.md Known Errors).
* HDX analysis, all master lines:
  * the PF = 1000 cap in `4.calc_D_uptake.py`;
  * the `legacy_T_range` default of 1.14 (about 400 K);
  * the 0.593 kcal/mol factor (section 3.1);
  * MBAR without the shared reference-energy subtraction in two example scripts.
* Two parameter derivatives are incomplete, so these parameters cannot be trained:
  * `HBondEnergy::get_param_deriv` returns entries 0-3 only;
  * `BackboneSigmoidCoupling::get_param_deriv` returns `scale` only.
* The stage-string trap (`martini_hybrid.cpp:637-647`): the seeds were fixed, the code was not.

First pass, 10-19, MacBook Pro (`audit_report.md`):

* **Build flags differ by platform.** Master builds with `-ffast-math` on both platforms, which compiles out
  `isfinite` tests. This branch added `-fno-finite-math-only` to the cluster's `CMakeLists_Other.txt` (08-31) and not
  to `CMakeLists_M1.txt`.
  * A probe under the Mac flags returns `isfinite(NaN)` true.
  * So the replica-swap rejection at `main.cpp:519` is probably compiled out on the Mac.
  * Not yet shown in an engine run.
* **AddressSanitizer** (chig REMD with pivot moves, 1dfn with jump moves): `DihedralCoord` (`coord_basic.cpp:343`)
  reads atom -1 before testing for the dummy. The value is not used. UBSan is clean.
* **`Concat`** (`coord_map.cpp:174, 205`):
  * the constructor reads its own uninitialised member;
  * the size check compares the wrong quantity;
  * `propagate_deriv` swaps two indices.

  No config uses it today.
* **Candidates, unconfirmed:**
  * `CurvatureChange` (`main.cpp:553-610`):
    * its multiplicative proposal has no Hastings factor;
    * it draws from an unseeded `rand()`;
    * it shares the replica-exchange random stream.
  * the thermostat treats mass-0 atoms as mass 1;
  * an HDF5 handle is closed twice on a load error;
  * node types are matched by longest prefix;
  * `MonteCarloSampler` has no virtual destructor.
* **Warnings:** 708 unique warnings under `-Wall -Wextra -Wshadow -Wconversion`, nearly all size and float
  conversions. The only real bug among them is `Concat`'s.
* **No test suite.** The only built-in check is `--potential-deriv-agreement`, one whole-system RMS that round-off
  dominates at this size (findings 12b).
  * There are no per-node force or parameter-derivative tests.
  * There is no parity harness, and no CTest, pytest or CI.
  * The finite-difference and parity scripts written for past fixes are single-use, spread over three machines.

### 2.2 Method

1. **Harnesses first** (3-4 days):
   * fixtures covering every `advanced_config` option, plus a regenerated hybrid seed (the August checkpoints no
     longer load);
   * per-node finite differences of forces and of parameter derivatives, with tie maps, against a float64 oracle;
   * a numpy oracle of each writer's index meaning;
   * a bitwise regression harness, so a fix is shown to change only what it means to;
   * a sanitizer and build-flag matrix;
   * invariance and statistical tests: equipartition, MC detailed balance, replica exchange on a harmonic ladder.
2. **Units of review.**
   * 18 units cover every master file exactly once (P1). They are engine graph, coordinate nodes, placement, BP,
     interaction graph and splines, H-bond, environment, local backbone terms, membrane, driver and replica exchange,
     MC, I/O, config writers, parameter glue, analysis, examples, build, and shipped parameters.
   * 7 units cover the branch-only code (P2). They are MARTINI potential, hybrid coupling, g-JF integrator, box and
     barostat, hybrid preparation, hybrid analysis, and ConDiv.
   * A writer/reader registry lists 80 node prefixes, so each node is reviewed from its C++ reader and its Python
     writer together. Findings 3.13 and 4.17 were writer/reader mismatches. Findings 3.12 sat in a
     parameter-derivative path that energies and forces never use.
3. **Two independent reviewers per unit.** The trained nodes go first: interaction graph, environment, BP, H-bond,
   placement, config writer.
4. **A claimed bug counts only after someone other than its finder reproduces it.** It needs a failing measurement
   or a minimal reproducer. Refuted claims are recorded too.
5. **A unit closes** when:
   * both reviews are in;
   * every node passes force and parameter finite differences and the invariance tests;
   * the sanitizers are clean;
   * every registry pair is checked.
6. **Fixes go into the Upside 3 tree one at a time.** Each fix needs three things:
   * its reproducer passes;
   * its node passes the finite-difference tests;
   * the regression harness shows every other path bitwise unchanged.

   Its effect on ff2.1 is not measured (section 1).
7. **Effort:** about 110 reviewer-hours for P1 and 45 for P2 (the audit agent's estimate).

## 3. Units (goal 2)

### 3.1 What was found (`units_report.md`)

**Upside's units are reduced.**
* kT at T_up = 1 is 1 E_up.
* Lengths are in A.
* Every backbone atom has mass 1 m_up (12 g/mol).
* The implied time unit is 1 A x sqrt(12 g/mol / 2.915 kJ/mol) = 0.203 ps, so dt 0.009 is 1.83 fs.
* Master's README says Upside times "should not be interpreted as picoseconds".
* The hybrid runs its own calibrated clock of 40 ps per step, whose origin is not recorded.

**The kJ/mol value of E_up is a calibration, not a unit.**
* 350.588 K is 298 K / 0.85. Peng 2022 SI eq. S14 matches T_up 0.85 to 298 K for FF2. FF1's match was 1.0, which
  gives 2.479 kJ/mol per E_up.
* Writing the parameter files in kJ/mol would bake one force field's fit into every number.
* A retrained Upside 3 needs its own calibration.
* The MARTINI conversion factor is this same calibration: the hybrid's protein-lipid balance depends on it.

**Some terms are evaluated at a fixed reference kT, at any simulation temperature:**
* the Ramachandran `-ln P` and its reference state;
* the sheet weights;
* the rotamer `-ln p`;
* the BP free energy (`rotamer.cpp:252, 300, 466, 896`);
* `weighted_pos` (`environment.cpp:518`).

Any change of units keeps them at kT_ref. Replacing kT_ref with kB T_sim would change the model.

**Where units sit in the code:**
* C++: one conversion constant (`FORCE_UNIT`, `surface.cpp:15`), plus 45 implicit sites in 20 files.
* Python: about 67 groups of sites in about 30 files.

**Defects found, all for the audit:**
* The kelvin conventions conflict, in master:
  * `2.mbar_meltingCurve_freeEnergy.py:109` and `calc_hdx_ht.py:337` use 0.593 kcal/mol per T_up, which implies
    298 K;
  * `1.get_force.py:15` uses 41.4 pN per E_up/A, which is kT at 300 K; under the 350.6 K convention the factor is
    48.4;
  * `FORCE_UNIT` is 1 bar only if E_up = kT(298 K).
* `martini_prepare_system_lib.py:1968, 2017, 2026` hard-code `*10.0`, bypassing the explicit length conversion.
* Several conversion literals are duplicated.
* `up.md` errors: `/output/kinetic` is the per-atom mean, not the total; it cites a script that does not exist; it
  lists retired nodes.

### 3.2 Options

| | engine | parameter files | dry-MARTINI conversion | ConDiv and past results | old `.up` files |
|---|---|---|---|---|---|
| A+: reduced units inside, standard units at every input and output | unchanged, bitwise | unchanged | stays, with the factor taken from the force field's recorded calibration | unchanged | run as before |
| B: kJ/mol, K, A, amu | about 28 sites; not a coherent set (the time unit is 0.1 ps) | every energy x 2.915 | length factor only | learning rates rescaled; every E_up number relabelled | converter and format marker |
| C: GROMACS units (nm, ps, kJ/mol, amu, K) | about 42 sites, every length included | energies and lengths | none, except the calibration | as B | as B |
| D: kcal/mol, A | as B | as B | still both factors | as B | as B |

**Recommended: A+.**
* Conversion moves out of the engine into one units `Config` module.
* Each force field records its calibration (`T_ref_K`) in its parameter directory.
* Every `.up` file carries `/input/units`.
* Input and analysis take and report K, kJ/mol or kcal/mol, bar and pN through that module.

A+ also leaves the engine bitwise unchanged. Since Upside 3 is a new release that is no longer required, so the
recommendation rests on the calibration: kJ/mol in the files would carry one force field's fit. C is the coherent
target only if Upside 3 is retrained from scratch with no comparison to ff_2.1, or the hybrid becomes the main
product. Even then the kT_ref terms keep an explicit reference temperature.

### 3.3 Plan

1. [ ] The units module and `/input/units`.
2. [ ] Every Python site moved onto the module.
3. [ ] Fix the conflicting factors and the hard-coded `*10.0` above.
4. [ ] Verification:
   * engine bitwise on the 8-config x 3-frame test;
   * each derived constant equal to the literal it replaces to 1e-12, choosing explicitly where repo literals
     disagree (14.521180763676 is off by 1.3e-7);
   * configs regenerated by the old and new writers byte-identical in `/input`.
5. [ ] Upside 3's temperature calibration at validation, by Peng's eq. S14 procedure. Rebuild the MARTINI tables if it
   moves.

## 4. New parameters (goal 3)

These are chosen after the side-chain verdict and the Phase A baseline (section 5.6). Each one:

* answers a measured defect, with its predicted effect written down before training;
* goes in the term that carries the defect (findings 10.13);
* is trained as the only change against a control (findings 10.17);
* passes the audit's parameter-derivative harness before any training;
* is a spline table read from `.h5`, so the GPU kernels can carry it too.

Candidates on record:

* **Parameters that already exist but cannot be trained:** H-bond entries 4-11 and backbone desolvation's centre and
  sharpness (section 2.1). These need derivatives, not new terms.
* **From goal 4, only if training changes leave a gap (section 5.6, B3):** the surface balance of backbone
  desolvation, and a cooperativity term.
* **The glycine context term** (`architecture.md`), only if its trigger fires.
* **The side-chain model's own parameters,** which follow the verdict.

## 5. Large proteins (goal 4)

Investigation: plan.md Phase 22, 2026-10-19.

### 5.1 The problem

* Upside folds large proteins de novo poorly (user, 10-19).
* Jumper 2018 (PLoS CB e1006578): all 20 test proteins under 100 residues come within 5 A of the native from an
  unfolded start. The CASP11 targets of 115-178 residues do not; T0855-D1 (115) never found its fold. The authors
  write that ConDiv's short runs "optimiz[e] only over fluctuations in or near the native well".
* Peng 2022 (FF2, JCTC 18:550): reversible folding only for the 40- and 43-residue mini-proteins. Ubiquitin (76) did
  not refold within a CPU-week because of long-lived misfolded states.
* Nothing above 150 residues has been trained or benchmarked de novo:
  * the 456 training proteins are 50-150 residues (median 114);
  * the 44-domain panel is 50-75;
  * the Peng de novo set is at most 92.

### 5.2 What was measured

**The native-start size trend is the metric.**
* On the 89-domain mdCATH native-start panel, ff_2.1_ctrl38 folds 0.42, 0.44 and 0.24 of frames in the 50-100,
  101-150 and 151-250 bins.
* "Folded" is Q at or above the all-atom 5th percentile. That cutoff tightens with N, because the all-atom Q spread
  narrows.
* Upside's mean Q sits about 0.13 below all-atom in every bin, and its Rg at 1.01 of native. With a fixed-offset
  cutoff the trend disappears.
* So Upside holds large natives; the problem is reaching them.

**De novo test.**
* Proteins: T4 lysozyme 2LZM (164 residues), RNase HI 2RN2 (155) and A. vinelandii flavodoxin II 1YOB (179, apo
  chain), all three refolded in vitro without a ribosome. The engrailed homeodomain 1ENH (54) is the control.
* ff_2.1 (released), 60,000 tu per run, on the Mac Studio.
* Starts:
  * native, with the panel's anneal (T 0.05 to T over t 96-400);
  * bulk: a random-dihedral chain, same anneal;
  * vectorial: 30 residues annealed for 2,000 tu, then 5 residues added every 1,000 tu at T, each new segment PPII
    and turned away from the chain, then the full chain to 60,000 tu.
* Every last frame was quenched for 1,000 tu at T 0.3.
* Q is the CA native-contact fraction: |i-j| >= 4, d < 8 A, formed below 1.2 d0. Nonlocal Q counts only |i-j| > 12.
* H-bonds were counted from H and O placed on the backbone (H...O < 2.5 A), which reproduces the crystal counts.

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

No coil forms native tertiary structure at any of the three temperatures.

**What the chains end in.**
* At T 0.80 every de novo run ends in a compact globule with almost no native tertiary contacts. Rg is 0.9-1.3 of
  native for the large proteins, 1.1-1.7 for the control.
* The globules are helix-poor: they form 0.36 (T4 lysozyme) to 0.64 (homeodomain) of the native start's helical
  H-bonds.
* T4 lysozyme's globules carry 1.5-5x the crystal's nonlocal, sheet-type H-bonds.
* Native starts keep crystal helix and H-bond content. Their Q drifts by up to about 0.2 over 60,000 tu and mostly
  plateaus.

**Energy by term, quenched globule minus quenched native.** Bulk / vectorial against native, each the mean over the
last 10 frames of every quench (6, 6 and 3 runs), in E_up. The terms sum to the engine total.

| term | 1ENH | 1YOB | 2LZM | 2RN2 |
|---|---|---|---|---|
| H-bond | +34 / +50 | +115 / +118 | +140 / +142 | +160 / +144 |
| rotamer (side-chain free energy) | +23 / +32 | +35 / +39 | +42 / +23 | +24 / +37 |
| side-chain burial | +2 / +3 | +20 / +20 | +13 / +15 | +13 / +15 |
| Rama | -5 / -5 | -7 / -19 | -20 / -27 | -19 / -32 |
| backbone desolvation | -4 / -7 | -17 / -17 | -18 / -19 | -16 / -15 |
| total (bulk / vectorial) | +49 / +73 | +155 / +150 | +179 / +159 | +163 / +149 |

* Most of the gap is backbone H-bonds the globule never forms.
* The globule's side chains recover 0.87 (1ENH), 0.91, 0.90 and 0.93 of the native's side-chain energy (rotamer plus
  burial, bulk), so collapse costs little.
* Rama and backbone desolvation favour the globule, by 24-47 E_up in the large proteins.

**The co-translational hypothesis.** The user's hypothesis was that folding is co-translational (S. Fried, JHU:
about a third of E. coli proteins do not refold in vitro), so a thermodynamic force field cannot fold large proteins
from a denatured chain. On these proteins:
* N-to-C growth equals bulk collapse in every measure, and is slightly worse on the homeodomain.
* All three large proteins refold in vitro without a ribosome, so Upside fails at something the real chains do on
  their own.

The hypothesis may hold for obligately co-translational folders, but it does not explain these failures. The growth
here has no ribosome tunnel or surface, and its rate (200 tu per residue) is not mapped to real time.

**What these runs do not show:**
* the globules' free energy against the native's: only energies were compared, not entropies;
* a size contrast under one protocol: the control does not fold at constant T in 60,000 tu at any of the three
  temperatures, so the size dependence rests on Jumper 2018 and the user's observation;
* generality: four proteins, one force field, constant temperature.

### 5.3 Reading

For these in-vitro refolders, Upside's failure is search.
* RNase HI's native is stable at T 0.88, and its coils do not find it at 0.80 or 0.88.
* Upside ranks the native 150-180 E_up below the globules on average, and 90-110 below the best of them.
* Yet the chain collapses on side-chain burial before it forms secondary structure. Burial already pays 87-93% of
  the native's side-chain energy.
* Inside the globule the chain cannot form native H-bonds, and Rama and backbone desolvation reward the
  non-H-bonded compact backbone.

Mechanisms, ranked by the evidence (literature in `lit/review.md`, code in `code_audit/audit.md`):

1. **No negative design against compact misfolds.**
   * ConDiv contrasts a native-restrained replica with 12 free replicas started at the native (8,000 tu, T 0.8-1.1)
     and an expanded self-avoiding-walk ensemble.
   * Compact non-native globules never enter the contrast, and their number grows with N (Bryngelson and Wolynes
     1987).
   * Decoy-based training is how AWSEM (T_f/T_g, Davtyan 2012) and OPEP (energy gap, Sterpone 2014) handle this.
2. **The training range and weighting.**
   * Training stops at 150 residues.
   * Each protein's contrast is summed, not normalised per residue (`ConDiv.py:403-405`), so 110-150-residue proteins
     dominate.
   * The fixed 12-replica ladder loses acceptance with N: 0.63 / 0.55 / 0.49 across the training bins, and 0.40
     predicted at 200.
   * Airas and Zhang 2025 found that adding larger structures to training is what fixed size transfer.
3. **Backbone desolvation acts as a surface term.**
   * It rewards exposed NH and CO (scale -0.30).
   * The native's fraction of exposed groups falls from 0.27 (50-100 residues) to 0.19 (251-400), while almost every
     group of an unfolded chain is exposed. So the term's pull toward unfolded grows with N.
   * Here it also favours the globules, by 15-19 E_up.
4. **Weak cooperativity.** Jumper 2018 reports a transition "not fully cooperative". Additive pair terms fail the
   calorimetric two-state test (Chan 2000).

Ruled out, or not size-dependent:
* explicit size terms or neighbour-list capacities: none exist;
* BP convergence: 14-28 iterations from 76 to 271 residues;
* the burial window: the same coverage range at every size, with residues above it rising from 5.3% to 10%;
* the pivot sampler: it loses efficiency with N but cannot bias.

### 5.4 Engine defect fixed on 2026-10-19 (findings 3.13)

`environment_coverage_hb` gives the side-chain coverage of backbone H and O that feeds backbone desolvation. It
labelled each covering side chain with the type of the H/O site whose index equals that side chain's residue number.
The writer stored one type per site (`py/upside_config.py`), and the engine read `aa_types[id2]` with `id2` the
residue (`src/environment.cpp`).

* Only 7-9% of covering side chains got their own residue's type. The coverage weights run from 0.56 (ASP) to 2.12
  (TRP), so each NH/CO's burial was a sum over scrambled weights.
* The same code is in master, and the defect does not depend on size.
* **Fix.**
  * The writer stores one type per residue.
  * The engine sizes `aa_types` from its dataset and requires every `id2` to index it.
  * The side-chain burial node already held per-residue types and is unchanged for every allowed config (the parser
    blocks residue exclusions).
* **Old files run as master.** Files written before the fix keep their per-site arrays and run exactly as master.
  * All energies, terms, coverages and forces are bitwise equal on 8 configs at 3 frames each.
  * An engine without the fix refuses a new file with a size error, so a new config cannot run silently on an old
    binary.
* **Verified.**
  * Both builds match a numpy reference of the coverage to 5e-7, each with the types it should use.
  * The finite-difference derivative of the energy change matches the change in force (ratio 1.0001 at h = 0.005, 40
    coordinates on two globules).
* **Effect on ff_2.1 is small.**
  * Backbone desolvation changes by -0.1 to +0.2 E_up at the natives and +0.05 to +1.3 at the globules, with no other
    term moving.
  * So the defect is not a cause of the large-protein failure.
* **Native-start dynamics.** With the fixed engine on new configs at T 0.80 (t 2,500-5,000 tu):
  * KE/(1.5 kT) is 1.007-1.010, and the peptide C-N bond 1.33-1.34 A, as before;
  * Q, Rg and potential fall within single-run noise of the overnight runs: Q 0.83 / 0.69 / 0.81 / 0.74 against
    0.84 / 0.67 / 0.82 / 0.74 for 1ENH, 2LZM, 2RN2, 1YOB.
* **Where it is installed.**
  * In the repo and in the Mac Studio's `obj/` (backups `*.bak_pre_envfix_20261019`).
  * The MacBook Pro's `obj/` is the 10-02 build, which carries neither this fix nor 3.12's. It refuses configs from
    the current writer.
  * It is not deployed to `$P` or `/beagle3`, because running trainings write their configs with the old writer.
  * ff_2.1 and every training so far ran with the defect; Upside 3 trains with it fixed.

### 5.5 Not planned: a co-translational protocol

It did not help here, and the large test proteins fold in vitro without one.

### 5.6 Plan

**Phase A: baseline on the chosen side chain and the audited engine.**

1. [ ] Deploy the audited engine to the Upside 3 training tree, a tree of its own:
   * a fresh midway2 build;
   * check that it passes the audit harnesses and is bitwise equal to the locally tested build on their fixtures;
   * install by rename, carrying file modes (findings 10.11).
2. [ ] Fix the de novo benchmark.
   * 8-12 single-domain proteins of 50-250 residues with published in-vitro refolding, spread over the alpha, beta and
     alpha/beta classes, the four here included.
   * Score by nonlocal Q and a size-normalised RMSD (RMSD_100) or TM-score, from coil, with replica exchange.
3. [ ] Replica exchange on RNase HI and the homeodomain.
   * 16 replicas over T 0.80-1.00, from coil and native starts.
   * About 5 h and 15 h on the Mac, or about 27 h on midway2.
   * This answers 5.2's two open questions:
     * if the globules beat the native in free energy, the energy ranking must change;
     * if the native wins but coils never reach it, the funnel shape is the target.
4. [ ] Repeat 5.2 under the chosen side-chain model. The globule's side-chain energy fraction, and the Rama and
   desolvation preference for globules, are Phase B's targets.

**Phase B: training changes, one at a time, each scored on A2 and on the native panels (44-domain, mdCATH).**

1. [ ] **Compact misfolds as negative examples.** Make decoys by collapsing each training protein from coil at T 0.80
   (the bulk protocol, then quench), regenerated from the current force field each epoch. Two forms:
   * start some free replicas from decoys, which keeps ConDiv's objective and fixes the sampling bias Jumper 2018
     names;
   * add an explicit native-against-decoy gap term next to the SARW term.

   The first is the smaller change. Targets:
   * the native-globule gap per residue;
   * the globule's side-chain energy fraction (now 0.87-0.93);
   * Rama plus desolvation favouring globules (now 24-47 E_up in the large proteins).
2. [ ] **Train past 150 residues.**
   * Add mdCATH domains of 150-250 residues, kept apart from the panel's held-out domains.
   * Normalise each protein's contrast per residue.
   * Space large proteins' replica ladders by 1/sqrt(N).
3. [ ] Only if B1-B2 leave a gap: the new parameters of section 4.

**Phase C: validation.**
* The native panels must not lose, and the de novo benchmark must gain.
* Re-check the hybrid (glpG), which runs on the protein force field.
* Calibrate Upside 3's temperature scale (section 3.3).

## 6. GPU (goal 5)

### 6.1 Measured (`gpu_report.md`)

Per force evaluation, ff_2.1 native starts, one Apple M3 Max core, today's source. There are 3 evaluations per step.

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

* **The protein model has no single hot spot.** A speed-up needs nearly every node on the GPU.
* **The hybrid is different.** In the hybrid (1rkl in DOPC, 3,812 particles), the dry-MARTINI pair loop is 86% of
  17.5 ms per step, and every Upside protein node together is 1.5%.
* **One system at a time cannot fill a GPU.** An evaluation is about 150 kernel launches at 3-5 us each, which is
  about the CPU's whole 0.75 ms for glpG. So the GPU must run many systems per launch:
  * a replica-exchange ladder;
  * a ConDiv minibatch, which is 24 proteins x 14 replicas = 336 systems.
* **What already suits a GPU:**
  * the graph is static;
  * all arithmetic is float32;
  * one OpenMP thread runs per system;
  * the RNG is counter-based Threefry, which runs under CUDA.
* **What does not:**
  * BP's data-dependent iteration count;
  * its per-node products of edge messages;
  * per-step allocations and map lookups in several nodes.
* **The trainer's parameter derivatives cost about 1/1000 of the simulation,** so they can stay on the CPU.

### 6.2 One force field on two backends

The force field already runs on two numerically different backends without retraining:
* The Mac and cluster builds differ in their reciprocal square-root refinement and in the finite-math flag (2.1).
* Between el7 and el8, forces are bitwise equal but trajectories differ from step 6.

A GPU is a third backend under the same standard. It must:
* use the same spline coefficient arrays, fitted on the host and uploaded, never refitted;
* use BP's same initialisation, damping, tolerance, stop test and synchronous schedule, with each system stopping at
  its own test (another schedule can reach another fixed point, which would be a model change);
* run float32 kernels with deterministic fixed-point force accumulation;
* use the same Threefry streams.

The tolerance has a precedent. The user ruled that the BP stop-test fix needs no retraining. That fix changed:
* static-frame energies by a median of 7.5e-4 E_up;
* forces by 8e-5 of their RMS;
* the ConDiv gradient by 0.1%;
* REMD T_mid by at most 0.002.

The GPU should sit well inside that.

Verification ladder:

* **L0.** Every parameter array, table, index and RNG key on the device is bitwise equal to the host's.
* **L1. Per node, on static frames** (about 20 proteins x 10 frames, including misfolded globules and hybrid frames):
  * energy within 1e-5 relative;
  * force RMS within 1e-5;
  * parameter derivatives within 1e-4;
  * BP compared at a tight tolerance and at the production one.
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

**Recommended: b with c.** `InteractionGraph<IType>` is already an evaluator pattern. Templating the evaluators on the
scalar type gives one implementation for both backends, and it also removes today's split between the SIMD
`compute_edge` and the scalar parameter derivative.

1. [ ] G1, after the audit closes:
   * refactor the CPU engine into the single-source form, bitwise equal to the audited CPU engine;
   * use the audit harnesses as its test.
2. [ ] G2, prototype:
   * side-chain pair edges, environment coverage and H-bond coverage on the GPU, batched;
   * go/no-go: one GPU evaluates a batch of 336 ConDiv-sized systems in about 2 ms, which matches 336 broadwl cores.
3. [ ] G3, the remaining nodes, after any new terms (section 4) are fixed.
4. [ ] G4, the BP kernels, after the verdict.
5. [ ] G5, the ladder L0-L4.

Throughout:
* The CPU engine stays the reference.
* Development happens on the Mac against the host backend; the Mac's GPU cannot run CUDA.
* d2 is a separate option, worth taking only if the hybrid is in Upside 3's scope.

**The side-chain verdict changes the size of the work.**
* With BP, the rotamer node is about half the protein time and holds the three hardest kernels. It is also the only
  node whose CPU-GPU difference can exceed rounding, because rounding can move the iteration at which the stop test
  fires.
* A deterministic side chain is an ordinary pair term, which saves about a quarter to a third of the GPU effort.
* That saving is no reason to choose it: it folds worse (0.264 against 0.365 at matched length, 0.241 converged;
  findings 1.21).

### 6.4 Hardware

* **midway2, measured.** Our account can use `gpu2`:
  * 5 nodes of 4 Tesla K80; QOS of 4 GPUs, 36 h, 10 jobs; billing counts CPUs only.
  * The K80 is Kepler, which CUDA 12 dropped, so it needs CUDA 11.x (`cuda/11.5` is the default module).
  * The V100 node has been down since 2024-06.
  * There are no Kokkos or OpenMM modules.
  * Fit for correctness testing, not for speed.
* **midway3, from RCC documentation, not measured.** The `gpu` partition has 1 node of 4 A100, 5 of 4 V100 and 5 of
  4 RTX6000.
* **beagle3, from RCC documentation.** 22 nodes of 4 A100 and 22 of 4 A40. RCC's documentation says beagle3 uses no
  SUs, while CLAUDE.md says it bills the group's GPU allocation. Ask before any use.

## 7. Decisions for the user

Answered by the user on 10-19:
* Upside 3 is a new release. No fix is measured against ff2.1, and ff2.1 is not re-validated on the audited engine
  (section 1).
* martini-dev is Upside 3. Fixes go only there, and master is not changed.
* ff2.1 keeps running for development. Later, Upside 3 refuses it.
* The audit starts now (plan.md Phase 24).

What follows from that: the restraint-pairing defect needs no ruling of its own. It is fixed in Upside 3 like every
other confirmed defect.

Still open:

1. Units: A+ (recommended) or C. Also, should `martini.h5` move to native units, with conversion at config time? That
   would change the CLAUDE.md unit contract.
2. GPU:
   * the route (b with c, prototype first);
   * which GPUs: midway3 `gpu` under which account; beagle3 once its billing is clear;
   * whether the hybrid's d2 is in scope.
3. B1's form (goal 4).
4. Items in the running Upside 2 campaign that await a NO GUARDS ruling:
   * the `run_remd.py` rollback;
   * the HDX items (2.1);
   * three checks in the hybrid code.

## 8. Files

`scratchpad/upside3_planning_20261019/` (this MacBook Pro):

| path | what |
|---|---|
| `audit_report.md`, `audit_work/` | defect inventory, audit units, writer/reader registry (`registry.txt`), warning, cppcheck and sanitizer logs, the parameter-derivative prototype (`harness/param_fd2.py`), the finite-math probe (`probe/fm.cpp`) |
| `units_report.md`, `units_work/` | every unit site with file:line, the parameter-file table, options and verification |
| `gpu_report.md`, `gpu_work/` | profile scripts and logs, replica scaling, the hybrid profile, `web_notes.md` with the sources of every documentation claim |

`scratchpad/large_protein_20261019/` (the Mac Studio):

| path | what |
|---|---|
| `sim/landscape.py`, `sim/common.py` | native, bulk and vectorial runs with quench |
| `sim/analyze_landscape.py`; `analysis/{ss_check,hb_check,terms}.py` | Q, Rg, quench; helix, local and nonlocal Q; H-bonds; energy by term |
| `analysis_T0{80,88,95}.txt`, `analysis/*_T0{80,88,95}.txt` | outputs |
| `runs_T080`, `runs_T088`, `runs_T095` | trajectories |
| `envfix/` | the fix's builds (pre and fixed), `coverage_test.py`, `compare.py`, `fd_test.py`, `dyn_check.py`, `input_diff.py` |
| `lit/review.md`, `code_audit/audit.md` | literature review, code audit |
