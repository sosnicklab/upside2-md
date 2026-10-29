# Upside 3

**Status (2026-10-29).**
* The engine audit (goal 1) closed on 2026-10-22 (4.3); the defect ledger is `findings.md` 3.14. The harness reran
  clean on HEAD on 10-29 (plan.md Phase 24 item 2c).
* Its validation on ff2.1, the gate of stage S1, is running on midway2 (4.4): 29 of 32 Peng arms complete at 08:25
  10-29, glpG REMD block 4 of 4.
* No Upside 3 training has run. Every ConDiv chain so far, Phase 33's included, runs on the pre-audit engine (5.4).
* Waiting on the user: the side-chain model and the starting parameters of the base force field (section 7), and the
  other decisions of section 11.
* Not blocked by any other step: S2's cluster sync and training tree (5.4), and the baselines that do not depend on the
  base (H2, B1), each started on the user's word as section 8 records. The units module and G1 wait on the units
  option.

Records: plan.md Phases 21-34 (history and status), `findings.md` 3.14 (ledger), `remote_jobs.md` (cluster jobs).
Files are listed in section 12.

## 1. Scope

### 1.1 What Upside 3 is

martini-dev is Upside 3, a new version that branches off master (user, 10-19). Master is the record of Upside 2 and
the coding-style reference; it is not a parity target and is not changed.

* **A new release, not a patch series on ff2.1** (user, 10-19). No list of fixed defects and their effects on ff2.1 is
  published, so no fix is measured against ff2.1 or master for attribution. Each fix is still proven correct (2.1);
  the debug as a whole is validated on ff2.1 runs (4.4), and Upside 3 as a whole at stage S6 (9.1). Master and ff2.1
  stay as released.
* **Clean slate.** No legacy layers, deprecated wrappers or options kept only for compatibility (CLAUDE.md).
* **ff2.1 during development.** ff2.1 configs keep running while Upside 3 is developed (user 10-22: "keep ff2.1
  useable for now"). Upside 3 refuses ff2.1 at release (S7, 9.3).

### 1.2 Goals and priorities

The user set five goals on 10-19:

1. **Debug the whole repo.** Find every bug ourselves, as Jumper did with Claude on master; master's code first
   (section 4).
2. **Unit alignment.** Run Upside in standard units without conversion, if that can be done (section 5).
3. **New parameters,** chosen after the side-chain verdict (section 8).
4. **Fold large proteins de novo** (8.3).
5. **GPU support,** with one force field for CPU and GPU, never trained twice (section 10). User 10-29: the GPU
   "should produce the same results as the current ffs, without retraining or training another version of ff".

Two later rulings shape goals 3 and 4:

* User 10-26: the large-protein fix must include a change in the force field's design, not only in its training.
* User 10-28: change what the training fits, "the highest priority is HDX. The one following it is large protein
  performance, then by IDR sequences".

Goals 3 and 4 are therefore one program of training changes (section 8), taken in that order.

Why the stages of section 3 run in their order:
* Every Upside 3 training runs on the audited engine, since each defect found was carried by the trainings that ran
  before it (4.1), and a training is weeks of cluster time.
* The GPU port is written against the audited engine in its final units, so that it is written once. Its target is
  the current force fields (10.2), so it does not wait for the training changes.

### 1.3 The hybrid

The dry-MARTINI hybrid (the Upside protein in an explicit dry-MARTINI membrane) is in the Upside 3 tree but is not
one of the five goals.

* It was audited with the engine: the seven branch-only review units (4.2).
* The user uses it as a validation case (10-22: "glpG hybrid simulation with ff2.1 ... as a validation for our
  entire debug").
* The user asked for hybrid work on Upside 3: the job scripts (10-22, A232) and glpG in a DDM micelle with ff3.1
  (10-26, plan.md Phase 31).

Whether it is a product goal is open (11, item 14). The answer decides whether the GPU's d2 route is in scope (10.3)
and whether `martini.h5` moves to native units (11, item 2).

### 1.4 Definitions

**Force fields.**

| name | what | where |
|---|---|---|
| ff2.1 | the released force field of the FF2 lineage (Peng 2022), probabilistic (BP) side chains | `parameters/ff_2.1` |
| c9, ff_2.1_ctrl38 | `ff21_ctrl_dt009`, ConDiv's control retraining of ff2.1 (BP); ctrl38 is its step-38 force field | `$P/training/` |
| `ff21_nobp_dt009` | the same with deterministic (no-BP) side chains | `$P/training/` |
| ff3.0 | ff2.1 with glycine's Ramachandran map repaired | `parameters/ff_3.0` |
| ff3.1 | released 10-26: the merge of ff3.1.1 and ff3.2.1; its `rama.dat` (`5b66c78a`) carries trained basin depths for every map; the other five files are ff2.1's | `parameters/ff_3.1` |
| `ff31_nativefix`, `ff32_nativefix` | Phase 33's trainings from ff3.1 on the fixed native target: depths only, and every term (user 10-28: "call it ff3.2") | `$P/training/` |
| the base | the starting parameters of Upside 3's training program; not chosen yet (section 7) | |

`$P` is the ConDiv tree `/project/trsosnic/yinhan/upside2-md-mdw2`.

**Labels.**

| label | meaning | section |
|---|---|---|
| S1-S7, G | the stages of the roadmap and the GPU track | 3 |
| H1-H5 | HDX steps | 8.2 |
| L1-L6 | large-protein steps | 8.3 |
| D1-D4 | the large-protein design changes, inside L6 | 8.3 |
| B0-B2 | IDR baseline steps | 8.4 |
| I1-I3 | IDR candidates | 8.4 |
| G1-G5 | GPU steps | 10.3 |
| V0-V3 | the GPU verification ladder | 10.2 |
| A+, B, C, D | units options | 5.2 |
| A1-A237 | the defect ledger (`findings.md` 3.14) | 4 |
| P1, P2 | review units: the 18 that cover master's code, the 7 that cover branch-only code | 4.2 |

**Terms.** An *energy term* is a node of the potential, evaluated from spline tables. An *objective term* is a part of
ConDiv's training objective: the native-state contrast (NSE), the denatured-state contrast (DSE), and the ones section
8 proposes.

## 2. Standards for every step

### 2.1 Engine changes

Every change to the engine or to a config writer follows the audit's method (4.2):

* **A claimed defect counts only after someone other than its finder reproduces it,** with a failing measurement or a
  minimal reproducer. Refuted claims are recorded too.
* **Fixes go in one at a time.** Each needs:
  * its reproducer to pass;
  * its node to pass `tests/node_fd.py` (forces and parameter derivatives against fourth-order differences; a
    flagged direction is settled by `tests/line_integral.py`);
  * `tests/regress.py compare` against the build before the change to show every other path bitwise unchanged.
* **Every parameter derivative is exact or raises "not implemented"** (A29-A31).
* **NO GUARDS, spline tables only, and each table equal to its published form** (CLAUDE.md).

### 2.2 Training changes

A training change is a new energy form, a parameter set released from freezing, or a new objective term. Each one:

1. **Answers a measured defect,** with its predicted effect written down before training. A candidate whose defect is
   not yet measured stays a candidate (D4 waits for its measurement).
2. **Goes in the term that carries the defect** (findings 10.13).
3. **Is a spline table read from `.h5`, if it is an energy,** so the GPU kernels can carry it too. New forms are
   allowed: D1-D3 are new forms.
4. **Passes the audit's parameter-derivative harness** before any training.
5. **Runs on the audited engine,** in the Upside 3 training tree (5.4).
6. **Names the parameter groups it trains.** Every other Adam group is held at its start. This is the user's Phase 29
   decision ("every Adam group held at ff2.1"), under which ff3.1's depth training left the native panels level with
   released ff2.1 (7.3). A change that must co-train the trained tables says so; 7.3 sets what that needs.
7. **Is trained as the only change against a control** (findings 10.17), and is judged on the native panels as well as
   on its own target.

**Accumulating accepted changes (proposal).**
* Changes are accepted one at a time, and the control of each change is the last accepted force field. The accepted
  changes accumulate into the release.
* Candidates of different tracks may train in parallel against the same control. When more than one passes, they are
  combined, and the combination is trained and tested again before it is accepted.
* Precedents (literature read 10-29):
  * a99SB-disp changed its parameters in sequence, each step refit and the final set tested on held-out proteins
    (Robustelli et al. 2018, PNAS, full text). The sequence was a torsion refit, then water dispersion, then the
    backbone O-H Lennard-Jones together with the torsions.
  * CALVADOS 2 tested its cutoff change alone before refitting its lambda set (Tesei and Lindorff-Larsen 2022, Open
    Res Eur, full text in part).
* Not taken: joint fitting of every target at once, as ForceBalance does with a weighted sum and a regularization (Wang
  et al. 2014, J Phys Chem Lett, abstract and manual). A ConDiv training takes weeks, and a joint fit cannot say which
  change caused a loss.
* What the held-out panels catch: Greener 2024 (Chem Sci, full text) fitted 108 parameters jointly on 4 folded and 4
  disordered proteins, and held-out folded proteins lost secondary structure (HEWL, BPTI).

**Acceptance rule (proposal; today only H5 has a threshold).** A change is accepted when:
* its target moves the predicted way, with a 95% interval that excludes zero; and
* on each native panel (2.3), the lower 95% bound of the paired difference against its control is at or above -δ. δ is
  fixed before the training, from the paired difference between two seeds of the control.

### 2.3 Test sets and reading rules

| set | size | observable | used for | where defined |
|---|---|---|---|---|
| 44-domain panel | 44 Charron CATH domains of 50-75 residues, 4 seeds; at each epoch end 176 runs at dt 0.015 | native start: folded fraction (Q at or above the all-atom 5th percentile) and TM-score, against all-atom ff99SB-ILDN at 300 K | native panel; every training's epoch ends | findings 1.37; `ff3_selection/` |
| mdCATH panel | 89 domains, 30 per bin of 50-100, 101-150 and 151-250 residues, 4 seeds | the same, against CHARMM22* at 320 K | native panel; the size trend | plan.md Phase 19 |
| Peng benchmark | 16 proteins, native and de novo, 14-replica REMD on Peng's Table S2 ladders | TM-score; Peng's FF2 0.55 native, 0.42 de novo | the debug validation; folding | findings 11b, 1.21 |
| de novo benchmark | 8-12 single-domain proteins of 50-250 residues with published in-vitro refolding (L1, to be built) | nonlocal Q and RMSD_100 or TM-score from coil, with REMD | large proteins | 8.3 |
| HDX test set | EHEE_rd2_0005, HEEH_rd4_0097, Ub psWT (H3) | censored Spearman of dG_HX; Upside FF2 0.61, 0.43 and 0.32 | HDX | 8.2 |
| IDR panel | 16 SAXS entries, 14 sequences (B0) | Rg against experiment, sequence dependence, pairs | IDRs | 8.4 |
| glpG hybrid | glpG-RKRK-79HIS in POPE/POPG, 48-replica REMD on T 0.70-0.90 | KE/1.5kT, C-N bonds, DSSP per TM helix, TM4 protection | the hybrid | plan.md Phase 24 item 2b |
| ConDiv held-out | 46 of the 456 training proteins (seed 20260928), held out of the depth update only; at its sample-size floor | gradient pull | training diagnostics | `training/ff31/ConDiv.py` |

* **Training sets stay apart from every test set,** checked by family (H4, I1).
* **Convergence follows findings 10.37:** the seeds of one entry must agree, and an entry whose seed spread exceeds its
  distance from experiment is reported as unresolved.
* **Reading rules are fixed before the data.** B0's rules (8.4) are the model.

## 3. Roadmap

| stage | work | waits on | output | gate | detail |
|---|---|---|---|---|---|
| S1 | engine audit; its validation on ff2.1 | nothing | the audited engine, A1-A237 | closed 10-22; validation gate below | 4 |
| S2 | units module; Upside 3 training tree; cluster tree synced | S1; the units module also on the units option (11, item 2) | a training tree running the audited engine and trainer | the tree's build bitwise equal to the locally tested build on the harness fixtures | 5 |
| S3 | base force field: side-chain model and starting parameters; the base's temperature calibration; the control training of 7.3 if a change needs it | S2; the user's choices (11, items 3-4) | the base, the training temperature mapping | 7.3 | 6, 7 |
| S4 | baselines (measurements, no training): H2-H4, L1-L3, B0-B2 | S3; H2, L3 and B1 can run earlier on ff2.1 and ff3.1 | each track's target numbers | each baseline's reading rules | 8 |
| S5 | training changes, one at a time: H5; L4-L6; I1-I3; the rest | S4, per track | accepted changes | 2.2 | 8 |
| S6 | validation and Upside 3's own calibration | S5 | the Upside 3 force field, its calibration, rebuilt MARTINI tables | 9.1 | 9 |
| S7 | release: Upside 3 refuses ff2.1 | S6 | Upside 3 | 9.3 | 9 |
| G | GPU, G1-G5 | G1: S1 and the units option; G2-G5: G1 | the current force fields on the GPU, with the CPU's results | V0-V3 | 10 |

**S1's validation gate (proposal; the records name the comparators but set no threshold).**
* Peng, native starts: the paired TM-score difference against the released-ff2.1 arms has an interval that includes
  zero. Measured 10-27: +0.012 [-0.001, +0.028].
* Peng, de novo: the difference lies within the spread of the seed replicates (running), or is diagnosed.
* glpG: KE/1.5kT, C-N bonds and TM helicity as Upside 2's, except TM4, whose loss is ff2.1's own (findings 3.15).

**What runs in parallel.**
* S2 needs no decision but the units option.
* S4's baselines run alongside each other: H3 and L2 on midway2, B1 and B2 on the Mac.
* S5's trainings follow 2.2's accumulation rule.
* G runs beside S2-S7. A term that S5 accepts is ported when it is accepted (10.3).

## 4. Engine audit (S1)

### 4.1 Defects known before the audit

Defects that trained force fields carry, found before the audit:
* the BP stop test (`rotamer.cpp`, the belief deviation without `fabs`; Jumper's report, fixed 2026-09-26): ff2.1 and
  ff3.0;
* findings 3.12, the pair term's angular parameter gradient (`interaction_graph.h`): every ConDiv run before 10-18;
* findings 3.13, `environment_coverage_hb`'s covering-residue types: ff2.1 and every training so far (below).

**The 3.13 defect** (fixed 2026-10-19).
* `environment_coverage_hb` gave each covering side chain the type of the H/O site whose index equals its residue
  number. So only 7-9% got their own type, and each NH/CO's burial summed scrambled weights (0.56 ASP to 2.12 TRP).
* The same code is in master; it does not depend on size.
* On ff2.1 the effect is small: backbone desolvation -0.1 to +0.2 E_up at the natives and +0.05 to +1.3 at the
  globules, no other term moving, native-start dynamics within single-run noise. So it is not a cause of the
  large-protein failure (8.3).
* ff2.1 and every training so far ran with it; Upside 3 trains with it fixed. Fix and verification: findings 3.13,
  `envfix/` (section 12).

**The first pass** (10-19, `audit_report.md`) seeded the audit. It found:
* build flags that compiled `isfinite` out on the Mac only;
* an AddressSanitizer read of atom -1 in `DihedralCoord`;
* three bugs in `Concat`;
* candidates in `CurvatureChange` and the thermostat;
* 708 compiler warnings with no real bug beyond `Concat`'s;
* no test suite. The harnesses of 4.2 are now `tests/`.

Its items went into the audit's units. The confirmed ones are ledger entries, for example A1 (`DihedralCoord`), A25
(`make_restraint_group`), A28 (the Mac build flags), A29-A30 (the incomplete backbone desolvation and H-bond parameter
derivatives), A37 (`CurvatureChange`), A55 (`Concat`), A62 and A64 (the HDX analysis factors), and A68 (the hybrid
stage-string trap).

### 4.2 Method

1. **Harnesses first:**
   * fixtures covering every `advanced_config` option, plus a regenerated hybrid seed;
   * per-node finite differences of forces and of parameter derivatives, with tie maps, against a float64 oracle;
   * a numpy oracle of each writer's index meaning;
   * a bitwise regression harness, so a fix is shown to change only what it means to;
   * a sanitizer and build-flag matrix;
   * invariance and statistical tests: equipartition, MC detailed balance, replica exchange on a harmonic ladder.
2. **Units of review.**
   * P1, 18 units covering every master file exactly once: engine graph, coordinate nodes, placement, BP, interaction
     graph and splines, H-bond, environment, local backbone terms, membrane, driver and replica exchange, MC, I/O,
     config writers, parameter glue, analysis, examples, build, and shipped parameters.
   * P2, 7 units covering the branch-only code: MARTINI potential, hybrid coupling, g-JF integrator, box and
     barostat, hybrid preparation, hybrid analysis, and ConDiv.
   * A writer/reader registry lists 80 node prefixes, so each node is reviewed from its C++ reader and its Python
     writer together. Findings 3.13 and 4.17 were writer/reader mismatches; findings 3.12 sat in a parameter-derivative
     path that energies and forces never use.
3. **Two independent reviewers per unit.** The trained nodes go first: interaction graph, environment, BP, H-bond,
   placement, config writer.
4. **Claims, fixes and verification** as 2.1.

### 4.3 Closure

A unit closes when:
* both reviews are in;
* every node passes force and parameter finite differences and the invariance tests;
* the sanitizers are clean; and
* every registry pair is checked.

Every unit closed on 10-22 (`scratchpad/upside3_audit/closing/`; plan.md Phase 24 item 2a). On 10-29 the harness
reran on HEAD on the Mac Studio (plan.md Phase 24 item 2c):
* regress was bitwise against the 10-27 records;
* the invariance tests passed;
* the index oracle was clean after its own fix (A237);
* node_fd's flagged directions line-integrated to noise, once the channel potential's steps were removed (4.5).

### 4.4 Validation of the debug on ff2.1

User 10-22: "submit Peng 16 validation protein with ff2.1 and glpG hybrid simulation with ff2.1, on midway2, as a
validation for our entire debug".
* Both run on the Upside 3 tree (`/beagle3/trsosnic/yinhan/upside3`).
* The comparators are the released-ff2.1 arms on the pre-audit binary and Peng's FF2 (Table S2, Fig S4).

Results so far (findings 3.14 "The core against Upside 2", 3.15; 10-27):
* **Hamiltonian at identical coordinates.** Only `bb_sigmoid_coupling_environment` (the 3.13 fix) moves at the E_up
  scale (median 0.20, max 1.05 E_up).
* **Peng, native starts.** Mean TM 0.564 against 0.552, +0.012 [-0.001, +0.028], Wilcoxon p 0.52. The engines agree.
* **Peng, de novo** (9 of 16 then). 0.437 against 0.461, -0.025 [-0.051, -0.003], p 0.051, carried by top7 (-0.108)
  and alpha3d (-0.071), both still rising between run halves.
  * Reweighting each engine's frames by the other's Hamiltonian predicts -0.002, so the energy change does not explain
    the gap.
  * The seed replicates (seeds 1, 101, 202) measure the run-to-run spread.
* **Open:** on BBA's random start `backbone_pairs` differs by 6.3 E_up and `rotamer` by 4.0, not diagnosed.
* **glpG.** It matches Upside 2 except TM4, which ff2.1 loses on Upside 3 already in preparation stage 7.0 (DSSP alpha
  0.44 at T 0.70 against 0.93-0.96).
  * ff3.0 holds TM4 on Upside 3 (0.88 at the end of stage 7.0), so the loss is ff2.1's own.
  * Which force field the hybrid's HDX uses is the user's call.

State at 08:25 10-29 (`remote_jobs.md`): 29 of 32 Peng arms complete; glpG blocks 1-3 done and block 4 running,
blocks 2-3 not read.
S1's gate is in section 3.

### 4.5 Open engine items

| item | what | blocks | where |
|---|---|---|---|
| U1 | the cluster binaries set flush-to-zero and denormals-are-zero in every process that loads the library (`set_fast_math`); the proposed fix is `-ffast-math` on compile lines only, checked by a regress compare of two midway2 builds | to be decided before S2's training build, so every Upside 3 training shares one build | plan.md Phase 24 item 2b |
| U2 | `sse2neon.h` is an unpinned Homebrew header whose defaults set the Mac's `_mm_rcp_ps` accuracy | Mac builds' reproducibility | same |
| U4 | the g-JF round evaluates the forces once more than it needs, about 20% of a hybrid's cost; skipping it needs to know that nothing moved the coordinates since the last evaluation | performance only; after correctness work | same |
| cluster tree | the Upside 3 tree carries A231, A234 and A235; A236 (`dihedral_spring`) and A237 are local only | S2 | `remote_jobs.md` |
| BBA | 6.3 E_up in `backbone_pairs` on one random start | S1's gate if it bears on the de novo gap | 4.4 |
| channel potential | energy steps with no force (below) | the decision of 11, item 15 | findings 3.14, after the ledger |
| container | A213 deleted `.devcontainer/`; a container for Upside 3 is written when there is a runtime to test it on | nothing | plan.md Phase 24 |

**Items of the running Upside 2 campaign that await a NO GUARDS ruling** (11, item 13):
* `run_remd.py` rolls NaN replicas back and continues, a guard under NO GUARDS. It is kept unchanged so that glpG
  stays comparable across campaigns, and it is not in the repo.
* The HDX analysis items of findings 5.5.
* Three checks in the hybrid code.

**The channel membrane potential** (`--channel-membrane-potential`).
* **Cause.** Its `surface` score counts the orientations, of 20 about the membrane normal, in which a residue's CB
  faces the lipid. The score has no derivative (`Surface::propagate_deriv` in `surface.cpp` is empty), so
  `cb_surf_membrane_potential` and `hb_surf_membrane_potential` change their energy with no force behind the change.
  MD feels none of it, while the logged energy, pivot MC and replica exchange read it, the mirror of A105.
* **Along each CB placement parameter** on 1rkl_channel: 88-149 steps over node_fd's stencil range, single steps up to
  1.4 E_up. With the steps removed the energy line-integrates to noise, so the derivatives are right.
* **In dynamics** (10-29; 1rkl_channel at T 0.85, dt 0.009, 180 time units, positions every step):
  * The score changes on 19,988 of 20,000 steps, a median of 8 of 31 residues per step.
  * The part of the channel energy no force carries changes by 0.22 E_up rms per step, 0.88 over 0.9 time units and
    1.43 (1.7 kT) over 9 time units. The channel energy's own standard deviation over the run is 1.66 E_up.
  * A smooth term on the same frames (`rama_map_pot`, Simpson's rule per step) leaves 0.012 and 0.33 E_up at the
    same lags.
  * The energy is exactly invariant under xy translation, as a flat membrane's must be.
  * One fixture at one temperature (`channel_steps.py`, `channel_simpson.py`, section 12).
* **Options for the user:** keep it as master's model; give the score a smooth form with a derivative; or retire it
  with the membrane retraining (8.5).

## 5. Foundations (S2): units and the training tree

### 5.1 Units: what was found (`units_report.md`)

**Upside's units are reduced.**
* kT at T_up = 1 is 1 E_up; lengths are in A; every backbone atom has mass 1 m_up (12 g/mol).
* The implied time unit is 1 A x sqrt(12 g/mol / 2.915 kJ/mol) = 0.203 ps, so dt 0.009 is 1.83 fs. Master's README
  says Upside times "should not be interpreted as picoseconds".
* The hybrid runs its own calibrated clock of 40 ps per step, whose origin is not recorded.

**The kJ/mol value of E_up is a calibration, not a unit.**
* 350.588 K is 298 K / 0.85, Peng's temperature calibration of FF2 (section 6). kB x 350.588 K = 2.915 kJ/mol is the
  E_up of CLAUDE.md's table.
* Writing the parameter files in kJ/mol would bake one force field's fit into every number, and a retrained Upside 3
  needs its own calibration.
* The MARTINI conversion factor is this same calibration: the hybrid's protein-lipid balance depends on it.
* Some terms are evaluated at a fixed reference kT at any simulation temperature (6.1). Any change of units keeps them
  there.

**Where units sit in the code.**
* In C++: 45 implicit sites in 20 files. Its one explicit constant, `FORCE_UNIT`, was corrected in A38 and deleted
  with its only user in A170.
* In Python: about 67 groups of sites in about 30 files.

**Defects found.**
* Fixed in the audit: the conflicting kelvin conventions (0.593 kcal/mol and 41.4 pN per E_up/A, which imply 298 K and
  300 K; A62), and up.md's `/output` errors (A182).
* Left for the units work: `martini_prepare_system_lib.py` hard-codes `*10.0`, bypassing the explicit length
  conversion, and several conversion literals are duplicated.

### 5.2 Options

| | engine | parameter files | dry-MARTINI conversion | ConDiv and past results | old `.up` files |
|---|---|---|---|---|---|
| A+: reduced units inside, standard units at every input and output | unchanged, bitwise | unchanged | stays, with the factor taken from the force field's recorded calibration | unchanged | run as before |
| B: kJ/mol, K, A, amu | about 28 sites; not a coherent set (the time unit is 0.1 ps) | every energy x 2.915 | length factor only | learning rates rescaled; every E_up number relabelled | converter and format marker |
| C: GROMACS units (nm, ps, kJ/mol, amu, K) | about 42 sites, every length included | energies and lengths | none, except the calibration | as B | as B |
| D: kcal/mol, A | as B | as B | still both factors | as B | as B |

**Recommended: A+.**
* Conversion moves out of the engine into one units `Config` module. Each force field records its calibration
  (`T_ref_K`) in its parameter directory, every `.up` file carries `/input/units`, and input and analysis take and
  report K, kJ/mol or kcal/mol, bar and pN through that module.
* The recommendation rests on the calibration (kJ/mol in the files would carry one force field's fit). It does not
  rest on A+ leaving the engine bitwise unchanged, which a new release no longer requires.
* C is the coherent target only if Upside 3 is retrained from scratch with no comparison to ff2.1, or if the hybrid
  becomes the main product. Even then, the kT_ref terms keep an explicit reference temperature.
* The choice is the user's (11, item 2). Under C the engine's sites change, and G1 waits for them.

### 5.3 Plan (A+)

1. [ ] The units module and `/input/units`. Each force-field directory carries its calibration record (`T_ref_K`,
   from its measured c, section 6). The same record is what S7 checks (9.3).
2. [ ] Every Python site moved onto the module.
3. [ ] The hard-coded `*10.0` and the duplicated literals (5.1).
4. [ ] Verification:
   * engine bitwise on the 8-config x 3-frame test;
   * each derived constant equal to the literal it replaces to 1e-12, choosing explicitly where repo literals
     disagree (14.521180763676 is off by 1.3e-7 relative);
   * configs regenerated by the old and new writers byte-identical in `/input`.

### 5.4 The Upside 3 training tree

**State.**
* No Upside 3 training tree exists.
* Every ConDiv run is in `$P/training/`. `$P/obj` is the pre-audit build plus only the 3.12 fix (md5
  `7136d99a`/`3d9c6366`; `remote_jobs.md` 0b), so every training so far, Phase 33's included, ran on the pre-audit
  engine.
* The audited trainer (A193-A199, A204; torch replaced by closed-form numpy, the 86 dead `rotscalar` coordinates
  removed) refuses older checkpoints. So the audited engine goes only to a tree no running campaign reads, and Upside 3
  campaigns start fresh (`remote_jobs.md` 0b).
* The Upside 3 tree on `/beagle3` is rsynced without `training/` and serves validation (CLAUDE.md "Shared Upside
  Deployment").

**Steps.**
1. [ ] Decide U1 (4.5).
2. [ ] Sync A236, A237 and later to the Upside 3 tree, and build on midway2 (CLAUDE.md "Compiling Upside on the
   cluster").
3. [ ] The training tree, a tree of its own:
   * the repo's `training/` with the audited trainer;
   * a fresh midway2 build that passes the audit harnesses and is bitwise equal to the locally tested build on their
     fixtures;
   * installed by rename and carrying file modes (findings 10.11).
4. [ ] If S3 keeps Phase 33's native target, move its `add_native_rama_restraint` from the runs' own `ConDiv.py`
   (md5 `54383004`) into the repo's trainer, through 2.1.

## 6. Temperature scale (S3, S6)

### 6.1 What sets temperature in Upside

* kT at T_up = 1 is 1 E_up.
* **Some terms are evaluated at a fixed reference kT, at any simulation temperature:**
  * the Ramachandran `-ln P` and its reference state;
  * the sheet weights;
  * the rotamer `-ln p` and the BP free energy (`rotamer.cpp`);
  * `weighted_pos` (`environment.cpp`).

  Any change of units keeps them at kT_ref; replacing kT_ref with kB T_sim would change the model.
* **No other energy depends on temperature.** The hydrophobic effect's heat capacity, which grows with N, is absent
  (D4).

### 6.2 Peng's calibration of FF2

Peng et al. 2022, J Chem Theory Comput 18:550, "Prediction and Validation of a Protein's Free Energy Surface Using
Hydrogen Exchange and (Importantly) Its Denaturant Dependence". SI section "Temperature calibration and Stability
prediction", read 10-29 from the local PDF (findings 11b gives its path).

* **Eq. S14:** c x R_up x T_up = 0.593 kcal/mol at 298 K, with R_up = 1. c is kcal/mol per E_up.
* **Eq. S16:** c = (0.593 dS_up + dG_exp) / dH_up.
  * dH and dS come from two-state fits of T-REMD melting curves (H-bond count or Rg) of 13 Rosetta-designed
    mini-proteins (Table S1), averaged.
  * dG_exp is their mean measured stability at 298 K, 2.71 kcal/mol.
* **For FF2,** c is 0.700 (Rg) and 0.702 (H-bonds), so 298 K is T_up 0.593 / 0.700 = 0.847. The text says 0.85; Fig.
  S9's caption says 0.86.
* **FF1 and FF1+.** The text maps 298 K to T_up 1.0 for FF1 (which makes E_up 2.479 kJ/mol) and 0.94 for FF1+, but
  Table S1's c gives 0.95 and 1.0, the reverse.
* **What is matched.** The observable is global stability, not Tm and not HDX. The main text calls it the single
  Upside temperature that best reproduces the experimental stabilities, with a mean error of 1.3 kcal/mol (<dTm> 7 K;
  literature agent, full text).
* **HDX readings.** Peng then read each HDX protein at its own temperature (293-305 K), chosen so that the simulated
  stability matched the HDX-derived dG_global (literature agent, full text).

### 6.3 The three mappings

1. **The training mapping (S3; proposal).**
   * Why: HDX and IDR targets are measured in kelvin, so the mapping from kelvin to T_up is fixed before H5 trains
     (8.2, condition 3). A term trained at a temperature chosen to match stability learns nothing about stability.
   * How: measure the base's c by eqs. S14-S16, then read every target at T_up = 0.593 (T_K / 298 K) / c_base.
   * The mapping stays fixed through S5, so a term that changes stability changes the model's stability instead of
     the temperature scale.
2. **Upside 3's own calibration (S6).** The released force field's c is measured the same way.
   * Its `T_ref_K` goes into its parameter directory (5.3 item 1).
   * The MARTINI tables are rebuilt if it moves, since the dry-MARTINI conversion uses the same calibration.
3. **Per-experiment readings.** Peng's practice, a temperature per protein matched to dG_global, removes global
   stability from the comparison and answers a different question. It may be reported beside the fixed mapping, never
   in its place. H1 read ubiquitin at 273 and 300 K.

Two points bear on the mappings:
* **HEEH_rd4_0097** is one of the 13 calibration proteins and one of the three HDX test proteins (H3). Proposal:
  report c with and without it.
* **B0's IDR runs** used T_K / 350.588, the FF2 mapping. Whether a linear mapping holds from 278 to 298 K is unproven,
  and B1 tests it on ACTR (8.4).

### 6.4 Steps

1. [ ] The base's c (S3): 14-replica T-REMD of Peng's 13 mini-proteins under the base, two-state fits, eq. S16. Peng's
   ladders and lengths are in his SI; the cost follows H3's estimate per protein.
2. [ ] Upside 3's c (S6), the same procedure on the release candidate; then the MARTINI tables if c moved.

## 7. The base force field (S3)

### 7.1 The side-chain model

Measured 10-20 and 10-21 (findings 1.21):
* **At the two fixed points** (both continued from step 38, at step 76): mdCATH folded fraction 0.241 for
  deterministic (no-BP) side chains against 0.353 for probabilistic (BP).
  * The difference is -0.111 [-0.160, -0.065], Wilcoxon p 1e-5, lower in 54 domains and higher in 18.
  * Folded falls in every length bin.
* **On TM-score at the fixed points:** mdCATH 0.667 against 0.702, and the 44 domains 0.623 against 0.650.
* **At matched length,** earlier: 0.264 against 0.365.
* **Verdict:** the deterministic side chain loses folding on held-out proteins.

Every current force field (ff3.0, ff3.1, the Phase 33 runs) uses BP. The choice is the user's (plan.md Phase 25:
"Choosing the side-chain model for Upside 3 ... is the user's"), and it has not been made (11, item 3). The GPU needs
the BP kernels whichever model is chosen, because the current force fields use BP (10.3).

### 7.2 The starting parameters

The candidates are released ff2.1, ff3.1, and a Phase 33 result once it ends (11, item 4).
* **ff3.1** is level with released ff2.1 on the native panels: mean TM +0.006 [-0.001, +0.012] on the 44 domains and
  +0.005 [-0.001, +0.012] on mdCATH, and above the control on both (+0.037, +0.031).
  * On Upside 3, glpG loses TM4 under ff3.1 as under ff2.1 (A234; findings 3.15).
  * Its Peng arms are running: 15 of 32 complete at 08:25 10-29 (plan.md Phase 30 item 4).
  * H2, its ubiquitin HDX, feeds this choice.
* **ff3.0** holds glpG's TM4 on Upside 3.
* **Phase 33's ff3.2** read after round 1 (10-28):
  * glycine's pair is back at ff3.0's balance;
  * its H-bond margin fell to +0.059, under the +0.10 line;
  * its epoch-0 panel folds 0.449, -0.167 [-0.206, -0.129] against released ff2.1.
* **ff31_nativefix** is at step 16.

All of these trained on the pre-audit engine. Their parameters can start Upside 3, but none is an Upside 3 training.

### 7.3 The trainer's fixed point

* **findings 1.37:** every force field this workflow trains from ff2.1 with its co-trained tables folds worse than
  released ff2.1, by 0.13 to 0.23 folded fraction on the held-out panels.
* **Phase 28** found no port defect. The loss is spread over three co-trained files: the side-chain table gives back
  +0.089, the backbone burial +0.081 and the H-bond +0.068. The learning rates set its speed, not its size, and the
  training proteins do not fold better either.
* **Why.** Peng and Jumper stopped on a fixed budget short of convergence on purpose, Jumper having seen accuracy fall
  when training continued (findings 1.23). So each retraining from ff2.1 is a second budget on top of FF2's.
* **Phase 33's native target does not cure it:** ff3.2's epoch 0, above.
* **The user's Phase 29 decision** holds every Adam group at ff2.1 and trains only the new parameters. Under it,
  ff3.1's depths left the panels level with ff2.1 (7.2).

Consequences for Upside 3:
* **A change that trains only its own parameters** (2.2 rule 6) is compared against its start. No control training is
  needed for the trainer's loss.
* **A change that co-trains the trained tables inherits the loss** unless the trainer is changed. Examples: an
  objective term that moves the burial or pair tables, D1 (replacing backbone desolvation's sigmoid), or ff3.2's
  all-term design. Phase 33's ff3.2 is the measurement of that case.
* **Proposal (S3's gate):** before the first change that co-trains the trained tables, a control training of those
  tables on the audited engine is read on the held-out panels. If it loses as findings 1.37's runs did, those changes
  wait for the trainer.
  * L4 is the leading candidate for the trainer's fix: its hypothesis, the short-run bias of native-started 8,000 time
    unit runs, is the same.
  * Every trainer change needs the user's word (plan.md Phase 28: every candidate fix departs from the FF2 workflow
    the project fixes; findings 10.17).

### 7.4 Steps

1. [ ] The user chooses the side-chain model and the starting parameters (11, items 3-4).
2. [ ] The base's c (6.4 item 1), which fixes the training mapping.
3. [ ] The baselines of S4 on the base (section 8).
4. [ ] If a change co-trains the trained tables: the control training and its gate (7.3).

## 8. Training changes (S4, S5)

### 8.1 Candidates and order

| track | step | kind | answers | stage | waits on | status |
|---|---|---|---|---|---|---|
| HDX | H1 | baseline, ubiquitin under ff2.1 and ff3.0 | | S4 | | done 10-28 |
| HDX | H2 | baseline, ubiquitin under ff3.1 | | S4 | ff3.1's arms complete | open |
| HDX | H3 | baseline REMD of the three HDX test proteins | | S4 | the base | open |
| HDX | H4 | training data | | S4 | the user (data use, 11 item 9) | open |
| HDX | H5 | objective term | the folded state's openings (H1) | S5 | H3, H4, the training mapping | open |
| large | L1 | the de novo benchmark and the base's run on it | | S4 | the base | open |
| large | L2 | REMD free energies, globule against native | | S4 | the base | open |
| large | L3 | 8.3's measurements repeated with ff3.1 | | S4 | ff3.1's validation | open |
| large | L4 | decoys through a replay buffer | no negative design (8.3) | S5 | L1, L2 | open |
| large | L5 | training past 150 residues on crops | the training range | S5 | L1, its pilot | open |
| large | L6 | design changes D1-D4 | the trap | S5 | L2, L3 | open |
| IDR | B0 | baseline under ff2.1, ff3.1 and the SARW | | S4 | | done 10-28 |
| IDR | B1 | ACTR at 298 K | the temperature mapping | S4 | | open (Mac) |
| IDR | B2 | longer sampling or REMD | the unresolved entries | S4 | the user's choice | open |
| IDR | B0 on the base | | | S4 | the base | open |
| IDR | I1 | objective term, IDR size | residue specificity (B0) | S5 | B2; after HDX and large proteins | open |
| IDR | I2 | the DSE target | not indicated by B0 | | | |
| IDR | I3 | local backbone balance, with D3 | undecided (B0) | S5 | B1, B2 | open |
| other | | 8.5 | | S5 | a measured defect each | |

Order within S5 is the user's priority: HDX, then large proteins, then IDRs. Inside large proteins, L2's branches set
D1-D4's order (8.3). The trainer's gate (7.3) comes before any change that co-trains the trained tables.

### 8.2 HDX

User 10-28: change what the training fits, with HDX first, large-protein performance second and IDR sequences third.
The energy's functional form stays spline tables; the change is in the objective, terms added to ConDiv's contrast.

**Why ConDiv does not fit HDX.**
* An amide's dG_HX is the free energy of the state in which it is open, -RT ln K_op. ConDiv's objective holds no free
  energy.
* NSE asks that free replicas started at the native stay where the native-restrained replica is. DSE asks that
  unfolded frames match a chain with no attractions.
* Neither sets how stable the fold is or which part of it opens first, so HDX agreement is whatever that fixed point
  leaves. Upside FF2 on Peng's three proteins: censored Spearman 0.61, 0.43 and 0.32.

**The terms.** Each is built from the per-frame dV/da the engine already returns (`ConDiv.py` `compute_divergence`).
* **HDX, per measured amide i** at the experiment's temperature, with MBAR weights over the replica ladder:
  d dG_i / da = <dV/da>_open,i - <dV/da>_closed,i, a contrast of the same kind as NSE.
  * The loss is sum_i (dG_i - dG_i,exp)^2.
  * An amide the model censors enters one-sided.
* **Global stability** is the same contrast between unfolded and folded frames. Folded and unfolded are read with the
  native panels' Q cutoff (2.3).
* Each term names the parameter groups it trains (2.2 rule 6).

**The data, and what a per-amide term can use** (literature read 10-29).
* **Ferrari et al.** Ramos Ferrari AJ, Dixit S, ... Guttman M, Rocklin GJ, "Large-scale discovery, analysis and
  design of protein energy landscapes", Nature 2026, doi:10.1038/s41586-026-10465-z; bioRxiv 2025,
  doi:10.1101/2025.03.20.644235.
  * Intact-protein mHDX-MS gives a distribution of opening free energies for each of 5,778 domains of 28-64 residues.
  * The authors write that the measurements "do not resolve which rates stem from which residues". Site-resolved data
    come from HDX-NMR on 13 domains.
  * 4,464 of the domains overlap Tsuboyama's set.
  * The data are public on Zenodo (doi:10.5281/zenodo.14983481, CC BY 4.0, 81 GB, registration requested).
  * Read: the preprint in part and the Zenodo record; the Nature page was not opened.
* **Tsuboyama et al. 2023** (Nature 620:434, full text in part): about 776,000 unfolding free energies by cDNA display
  proteolysis.
  * Coverage: 331 natural and 148 designed domains of 40-72 residues and their mutants, at room temperature, with a
    dynamic range of about 5 kcal/mol.
  * Data: Zenodo record 7992926.
* **So the per-amide term fits residue-resolved data only:** Peng's three proteins, which are the test set, and
  Ferrari's 13 HDX-NMR domains.
* **For the 5,778 domains, a domain-level term (proposal):** the model's distribution of opening free energies per
  domain against the measured one, with global stability from Tsuboyama. Its form is still to be designed.
* **Precedent.**
  * HDX has been used to reweight simulated ensembles: Best and Vendruscolo 2006, Structure (secondary source); HDXer,
    Bradshaw et al. 2020, Biophys J (abstract); Wan et al. 2020, JCTC (abstract).
  * For one protein, HDX fitted a protein-specific coarse-grained model's interaction strengths (Zhu, Li and Wang
    2021, Chin Phys B, AICG2+ on IkBa, abstract).
  * No transferable force field trained on HDX was found; the search was not exhaustive.

**Conditions, from H1** (ubiquitin, findings 5.7):
1. **Its own sampling.**
   * One ConDiv step (12 x 8000 tu) is 1/319 of the REMD that resolves ubiquitin's amides to 5.6 kcal/mol at 273 K.
   * So each HDX protein gets its own REMD, reused across updates by reweighting to the new parameters while its ESS
     holds (DiffTRe, Thaler and Zavadlav 2021; CALVADOS, Tesei et al. 2021), and rerun when the ESS falls.
2. **Signal only below the bound.**
   * Resolving experiment's 8.5 kcal/mol at 273 K would take about 200 such runs.
   * So the term acts on local and subglobal openings, and on small, marginally stable domains: HEEH_rd4_0097 at
     2.2-4.3 kcal/mol, and the 28-64-residue domains of Ferrari et al.
3. **A fixed temperature mapping** (6.3 item 1).
   * Peng reads FF2 at 300 K against an experiment at 273 K, at the temperature where the model's global stability
     matches. A term trained at a chosen temperature therefore learns nothing about stability.
4. **A frozen definition of open.**
   * The estimator reads the trained burial weights. ff2.1 to ff3.0 moves dG by 0.12 kcal/mol RMS through the
     definition alone.
   * The estimator keeps ff2.1's tables, also after D2b retrains the burial weights, so training cannot redefine
     "open" instead of changing the ensemble.
5. **A ladder that mixes.** Each protein's ladder spans its own measured transition, with walkers making round trips
   across it. ff3.0's ubiquitin ladder makes 3 round trips and its dG is not converged.

**What the term has to move, by H1.**
* ff2.1 opens the seven amides experiment assigns to global unfolding by fluctuations of the folded state: 3.5-7.0
  kcal/mol, against 7.7-8.5.
* ff3.0's ubiquitin melts about 0.06 T_up lower than ff2.1's, which the benchmark's bottom-rung TM-score hides.
* So the term must raise the folded state's opening free energies at those amides while holding global stability.

**H1's readout:**
* the pipeline reproduces Peng's FF2 (Spearman 0.90);
* a full benchmark REMD resolves 14-15 of the 24 measured amides at 273 K (bound 5.6 kcal/mol);
* the "open" definition alone shifts dG by 0.12 kcal/mol (condition 4).

**Steps** (user 10-28: "record the plannings. We will run it later"; none starts without the user's word).
- [x] H1 (10-28): ubiquitin on the Peng-benchmark REMD; findings 5.7, `scratchpad/hdx_term_20261028/` (section 12).
- [ ] H2: H1's readout on ff3.1's ubiquitin arms, once both hold a COMPLETE file. No new run.
  * The arms are `/beagle3/trsosnic/yinhan/ff3_benchmark/runs/ubiquitin_{native,denovo}_ff_3.1`.
  * H1's rules are kept unchanged, and H2 is read against H1's ff2.1.
  * It feeds the choice of the base (7.2). The script edits it needs are in section 12.
- [ ] H3: Upside 3 REMD of EHEE_rd2_0005, HEEH_rd4_0097 and Ub psWT, each under ff2.1 and under the base, from the
  native, on the Upside 3 tree (`env_upside3.sh`) through midway2.
  * **Reading.** With H1's rules, at each experiment's temperature on the training mapping (EHEE 298 K, HEEH 278 K,
    Ub 273 K), and at Peng's comparison temperatures.
  * **Protocol: Peng's for FF2** (SI pp. 12-13; Table S3):
    * EHEE and HEEH: 14 replicas on 0.80-0.995, 4,000,000 tu, the first 500 frames dropped.
    * Ub: three independent runs, 14 replicas on 0.80-1.02, at least 2,400,000 tu each, the first half dropped.
    * Frames every 100 tu in all of them.
    * The ladder changes only if a force field's transition leaves the range, and the change is recorded.
  * **Inputs.**
    * Peng's design models: `/project2/trsosnic/share/HEEH_rd4_0097.pdb` and
      `/project/trsosnic/yinhan/bioemu2_hdx/inputs/EHEE_rd2_0005_ref.pdb`.
    * Ub psWT: `/project/trsosnic/yinhan/bioemu2_hdx/inputs/1UBQ_psWT_ref.pdb`, 1UBQ with residues 45 and 68
      renamed, as Peng started from 1UBQ.
  * **Cost, estimated from H1's arm.** That arm (one node, 14 CPUs, 55 us per system-step) took six 10.5 h chunks for
    2.43M tu. So each Ub run is about 60 node-hours, and EHEE and HEEH about 50 each at half the size and 1.6 times
    the length: about 280 node-hours per force field, about 570 for both.
  * **Running.** Through the benchmark's `bench_run.py` chunking, or an equivalent driver, with the jobs recorded in
    `remote_jobs.md`.
  * **Purpose.** It gives the baseline a trained term must beat, and the three proteins stay held out of training.
    No Upside run of EHEE or HEEH exists on either cluster.
- [ ] H4: the training data.
  * Ferrari et al.'s 5,778 domains (public on Zenodo; the records say the user requested them through the authors'
    form, as for `diffusion_dynamics.md` Phase 2 E2). How they enter the objective is 11, item 9.
  * Tsuboyama et al. 2023 for global stability, processed into per-domain targets with the folded and unfolded
    definitions above.
  * Families are checked against H3's test set.
- [ ] H5: the HDX term in the trainer, after H3, H4 and the training mapping (6.3).
  * It must reach censored Spearman at least FF2's on H3's three proteins, with a block bootstrap
    (`diffusion_dynamics.md` Phase 2 E4's rule), and pass 2.2's acceptance rule on the native panels.
  * The Boltz add-on of `diffusion_dynamics.md` Phase 2 uses the same data and test, and both routes are read on one
    rule.

### 8.3 Large proteins

Investigation: plan.md Phase 22, 2026-10-19.

#### 8.3.1 The problem

Upside folds large proteins de novo poorly (user, 10-19).
* **Jumper 2018** (PLoS CB e1006578): all 20 test proteins under 100 residues come within 5 A of the native from an
  unfolded start. The CASP11 targets of 115-178 residues do not, and T0855-D1 (115) never found its fold. ConDiv's
  short runs "optimiz[e] only over fluctuations in or near the native well".
* **Peng 2022** (FF2, JCTC 18:550): reversible folding only for the 40- and 43-residue mini-proteins. Ubiquitin (76)
  did not refold within a CPU-week because of long-lived misfolded states.
* **Nothing above 150 residues has been trained or benchmarked de novo.** The 456 training proteins are 50-150 residues
  (median 114), the 44-domain panel 50-75, the Peng de novo set at most 92.

#### 8.3.2 What was measured (ff2.1)

**The native-start size trend is the metric.**
* On the 89-domain mdCATH panel, ff_2.1_ctrl38 folds 0.42, 0.44 and 0.24 of frames in the 50-100, 101-150 and 151-250
  bins. "Folded" is Q at or above the all-atom 5th percentile, a cutoff that tightens with N as the all-atom Q spread
  narrows.
* Upside's mean Q sits about 0.13 below all-atom in every bin and its Rg at 1.01 of native. With a fixed-offset cutoff
  the trend disappears.
* So Upside holds large natives; the problem is reaching them.

**De novo test** (ff2.1 released, 60,000 tu per run, Mac Studio).
* **Proteins.** T4 lysozyme 2LZM (164 residues), RNase HI 2RN2 (155) and A. vinelandii flavodoxin II 1YOB (179, apo
  chain) all refold in vitro without a ribosome. The engrailed homeodomain 1ENH (54) is the control.
* **Starts.**
  * native, with the panel's anneal (T 0.05 to T over t 96-400);
  * bulk, a random-dihedral chain with the same anneal;
  * vectorial: 30 residues annealed for 2,000 tu, then 5 residues added every 1,000 tu at T (each new segment PPII
    and turned away from the chain), then the full chain to 60,000 tu.
* **Quench.** Every last frame was quenched for 1,000 tu at T 0.3.
* **Measures.** Q is the CA native-contact fraction (|i-j| >= 4, d < 8 A, formed below 1.2 d0; nonlocal Q only
  |i-j| > 12). H-bonds are counted from H and O placed on the backbone (H...O < 2.5 A), which reproduces the crystal
  counts.

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

* **No coil forms native tertiary structure at any of the three temperatures.**
* At T 0.80 every de novo run ends in a compact, helix-poor globule: Rg 0.9-1.3 of native for the large proteins and
  1.1-1.7 for the control. T4 lysozyme's globules carry 1.5-5x the crystal's nonlocal, sheet-type H-bonds.
* Native starts keep crystal helix and H-bond content. Their Q drifts by up to about 0.2 over 60,000 tu and mostly
  plateaus.

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

* Most of the gap is backbone H-bonds the globule never forms.
* The globule's side chains recover 0.87 (1ENH), 0.91, 0.90 and 0.93 of the native's side-chain energy (rotamer plus
  burial, bulk), so collapse costs little.
* Rama and backbone desolvation favour the globule, by 24-47 E_up in the large proteins.

**The co-translational hypothesis does not explain these failures.**
* The hypothesis is the user's: folding is co-translational (S. Fried, JHU: about a third of E. coli proteins do not
  refold in vitro), so a thermodynamic force field cannot fold large proteins from a denatured chain.
* Here N-to-C growth equals bulk collapse in every measure (slightly worse on the homeodomain), and all three large
  proteins refold in vitro without a ribosome.
* It may hold for obligately co-translational folders. The growth here has no ribosome tunnel or surface, and its rate
  (200 tu per residue) is not mapped to real time.

**Not shown by these runs:**
* the globules' free energy against the native's (only energies were compared; L2 measures it);
* a size contrast under one protocol: the control does not fold at constant T in 60,000 tu at any of the three
  temperatures, so the size dependence rests on Jumper 2018 and the user's observation;
* generality: four proteins, one force field, constant temperature.

#### 8.3.3 Reading

For these in-vitro refolders, Upside's failure is search.
* RNase HI's native is stable at T 0.88, and its coils do not find it at 0.80 or 0.88.
* Upside ranks the native 150-180 E_up below the globules on average, and 90-110 below the best of them.
* Yet the chain collapses on side-chain burial, which already pays 87-93% of the native's side-chain energy. Inside the
  globule it cannot form native H-bonds, and Rama and backbone desolvation reward the non-H-bonded compact backbone.

Mechanisms, ranked by the evidence (literature in `lit/review.md`, code in `code_audit/audit.md`):

1. **No negative design against compact misfolds.**
   * ConDiv contrasts a native-restrained replica with 12 free replicas started at the native (8,000 tu, T 0.8-1.1)
     and an expanded self-avoiding-walk ensemble.
   * Compact non-native globules never enter the contrast, and their number grows with N (Bryngelson and Wolynes
     1987).
   * AWSEM (T_f/T_g, Davtyan 2012) and OPEP (energy gap, Sterpone 2014) handle this with decoys.
2. **The training range and weighting.**
   * Training stops at 150 residues.
   * Each protein's contrast is summed, not normalised per residue (`ConDiv.py` `run_minibatch`), so 110-150-residue
     proteins dominate.
   * The fixed 12-replica ladder loses acceptance with N: 0.63 / 0.55 / 0.49 across the training bins, 0.40 predicted
     at 200.
   * Airas and Zhang 2025 found that adding larger structures to training is what fixed size transfer.
3. **Backbone desolvation acts as a surface term.**
   * It rewards exposed NH and CO (scale -0.30).
   * The native's exposed fraction falls from 0.27 (50-100 residues) to 0.19 (251-400), while almost every group of an
     unfolded chain is exposed, so its pull toward unfolded grows with N.
   * Here it also favours the globules, by 15-19 E_up.
4. **Weak cooperativity.** Jumper 2018 reports a transition "not fully cooperative"; additive pair terms fail the
   calorimetric two-state test (Chan 2000).

Ruled out, or not size-dependent:
* explicit size terms or neighbour-list capacities (none exist);
* BP convergence (14-28 iterations from 76 to 271 residues);
* the burial window (the same coverage range at every size, residues above it rising from 5.3% to 10%);
* the pivot sampler (it loses efficiency with N but cannot bias);
* the 3.13 defect (4.1).

#### 8.3.4 The force field term by term (2026-10-26)

User 10-26: review every term and every possible design change, not one term.
* Each term is read on 8.3.2's T 0.80 quenched end states. "Globule" means the bulk and vectorial ones, and ranges
  are over the three large proteins.
* Scripts and their `*_T080.txt` outputs are in `scratchpad/large_protein_20261019/analysis/`.

**Every number in this subsection is ff2.1's,** run on the 10-19 engine, before the audit's fixes, from configs
written before findings 3.13's fix.
* L3 repeats them with ff3.1 on the audited engine (user 10-27), and the candidates are re-ranked on those numbers
  before any is trained.
* ff3.1 differs from ff2.1 only in its Ramachandran library, so the Ramachandran row and D3 are the most likely to
  move; the engine fixes can move any row.

**What the trap is.**
* The globule is 150-180 E_up above the native, so it wins on entropy if at all. It is stable against the unfolded
  chain, and the chain does not leave it.
* A design change can **widen the native's lead** (needed if L2 finds the globules ahead in free energy), or **make
  non-native compact states less stable relative to the unfolded chain while the native keeps its stability** (needed
  if the native wins but coils do not reach it).
* The trap's depth is the side-chain energy a globule collects without native structure.

**Every term, globule against native.** The term gaps are 8.3.2's; everything else was measured on 10-26.

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

**Kinetics, from the bulk runs: the chain does not collapse before forming helices.**
* For the first 2,000-5,000 tu it stays at 1.5-2.0x the native Rg with 0.12-0.40 of its helical bonds.
* It then collapses, to 1.1-1.3x by 20,000 tu, with 0.32-0.54.

**What the forms cannot express.**
* Burial and pairing are interchangeable in backbone desolvation: a buried H-bond gains 0.3 E_up over an exposed one,
  and a buried unpaired group pays nothing beyond the lost exposure reward.
* No side-chain energy depends on how dense a contact's surroundings are; the burial term counts side chains, not
  which ones.
* No term couples neighbouring residues' conformations. The maps condition only on the neighbours' identity, so helix
  nucleation and propagation come only from the H-bond geometry.
* No energy depends on temperature except through the kT_ref terms (6.1). The hydrophobic effect's heat capacity,
  which grows with N, is absent.

**Not taken, by measurement:**
* a penalty for buried unpaired NH/CO alone (0.01-0.11 groups per residue);
* burial beyond 9 A (no better at 14 A than at 8 A);
* long-range electrostatics (within 1 E_up beyond 7 A);
* a desolvation barrier in the pair table, meant to slow collapse, which does not outrun helix formation here.

**Not testable on end states, and large changes:**
* more beads per side chain, for packing finer than one bead;
* a backbone frame-frame orientation term (AlphaFold2's residue frames). The side-chain beads' angular pair term
  already carries residue orientation, and 2LZM's non-native ladders are the only hint.

**Not taken: structure priors.**
* MSA or language-model contacts keep AlphaFold2 accurate on large proteins, but would make Upside a structure
  predictor and bias its pathways.
* Diffusion samplers (Boltz) as Monte Carlo moves have no proposal density for the acceptance test, and import a
  structure prior too.

**Not planned: a co-translational protocol.** It did not help here, and the large test proteins fold in vitro without
one.

#### 8.3.5 Steps

**Baselines (S4).**

1. [ ] **L1. The de novo benchmark.**
   * 8-12 single-domain proteins of 50-250 residues with published in-vitro refolding, over the alpha, beta and
     alpha/beta classes, the four of 8.3.2 included.
   * Scored by nonlocal Q and a size-normalised RMSD (RMSD_100) or TM-score, from coil, with replica exchange; ladders
     as H3's condition 5.
   * The base's run on it is the baseline every later change is scored on.
2. [ ] **L2. Replica exchange on RNase HI and the homeodomain.**
   * Protocol: 16 replicas over T 0.80-1.00 (adjusted so the ladder spans each protein's transition, H3's condition 5),
     from coil and native starts. About 5 h and 15 h on the Mac, or about 27 h on midway2.
   * Force field: the base. If L2 runs before S3, it runs on ff3.1 and is read again on the base.
   * Estimator: MBAR over the ladder gives the free energy of the native basin (Q at or above the panel's folded
     cutoff) against the compact non-native basin (Rg within 1.3x native, nonlocal Q below 0.2), at T 0.80, with a
     block bootstrap.
   * Three branches:
     * native favoured, its interval excluding zero: the funnel shape is the target, and D2 goes first;
     * globules favoured: the energy ranking must change, and D1 goes first;
     * unresolved: sampling is extended, and no design change is chosen on it.
   * Existing runs cannot stand in.
     * The Peng benchmark's homeodomain arms have a ladder (0.87-1.00) above the transition: native-start TM-score
       0.19-0.31 at the bottom rung.
     * Ubiquitin under ff3.0 is the small-protein case of the same search failure: de novo never below 4.4 A in 7.2M
       tu, while the native start refolds.
3. [ ] **L3. Repeat 8.3.2 and 8.3.4 with ff3.1,** once its validation is in (plan.md Phase 30 item 4; user 10-27).
   * The same four proteins and protocol (`sim/landscape.py`: native x3, bulk x6, vectorial x6 at T 0.80, 60,000 tu,
     then quench), on the audited engine with configs from the current writer.
   * The 8.3.2 energy-by-term table first, then the six 8.3.4 scripts. The two that read the engine run with the
     current `py/upside_engine.py` in place of `eng_1019`.
   * Each 8.3.4 row is read against ff2.1's. A design candidate keeps its place in the order only if its target still
     holds.
   * If the base is not ff3.1, the same is read on the base.

**Changes (S5).** Each is trained as the only change, scored on L1 and on the native panels (44-domain, mdCATH), and
accepted by 2.2's rule. L4 and L5 combine with any change of L6.

4. [ ] **L4. Compact misfolds as negative examples, through a replay buffer.**
   * **Mechanism.** Each training protein keeps a buffer of its own states: globules collapsed from coil under the
     current force field (the bulk protocol at T 0.80, then quench), and its free replicas' end states from earlier
     visits. Some of the 12 free replicas start from the buffer instead of the native, and the buffer refreshes from
     their end states, so decoys are made once, not every epoch.
   * **Basis.** This keeps ConDiv's objective and reduces the short-run bias Jumper 2018 names (free replicas started
     at the native explore only its well). It is the standard remedy in energy-based model training: persistent CD,
     Tieleman 2008; replay buffer, Du and Mordatch 2019; Nijkamp et al. 2020.
   * **Knob.** The share of buffer starts, set on held-out proteins.
   * **Alternative form.** An explicit native-against-decoy gap term next to the SARW term (11, item 8).
   * **Targets.** The native-globule gap per residue; the globule's side-chain energy fraction (0.87-0.93); Rama plus
     desolvation favouring globules (24-47 E_up).
   * **Hypothesis, not tested.** findings 1.37's loss (training buys nothing on its own proteins while held-out folding
     falls) may come from the same bias, since a native-started 8,000 tu run carries little folding information. So
     L4 is also the leading candidate for the trainer's loss (7.3).
5. [ ] **L5. Train past 150 residues, on crops in a frozen context.** This is the physical version of AlphaFold2's
   256-residue crops fine-tuned at 384 (AlphaFold3 adds spatial crops), so an example's cost does not grow with the
   chain.
   * **Data.** mdCATH domains of 150-400 residues, kept apart from the panel's held-out domains.
   * **Crops.** In each domain, a crop of 60-90 residues moves, contiguous or a compact subdomain (a foldon). The rest
     is frozen at the native through `/input/fix_rigid/atom_indices` (`martini_fix_rigid.cpp`
     `register_fix_rigid_for_engine`), which any config can carry. The native-restrained replica restrains the crop
     as now.
   * **What the crop sees.** Large-protein cores and subdomain interfaces: the share of residues above the burial
     window rises from 5.3% to 10% with size (`code_audit/audit.md` 1a). On the hotter rungs it can partly unfold and
     repack, so the trainer also sees local misfolds in a native context.
   * **L5a, the pilot,** before the training:
     * how far a crop moves in 8,000 tu (the 54-residue control does not refold in 60,000 tu);
     * whether replica acceptance follows the crop's size;
     * the SARW replica's crop form, or the decision to skip it for crops.

     Each step still costs a whole-protein evaluation, but runs need not get longer.
   * **Weighting.** Each protein's contrast is normalised per residue (`ConDiv.py` `run_minibatch`), and each ladder
     is spaced by 1/sqrt(N) of its moving residues.
6. [ ] **L6. Force-field design changes D1-D4,** in the order below. Each is a spline table or the parameters of one,
   with exact parameter derivatives through the audit harness before any training (2.2).
   1. **D1. Backbone polar groups: one table of burial and pairing.**
      * **Form.** A 2D spline E(side-chain coverage, pairing) per NH/CO replaces backbone desolvation's sigmoid of a
        sum. It gives the H-bond a dependence on burial and sets the four corners (exposed or buried, paired or
        unpaired) separately.
      * **Physics.** Backbone H-bonds in nonpolar surroundings are worth more than exposed ones (amide-to-ester
        mutations: Deechongkit et al. 2004, Nature 430:101; Gao et al. 2009, Nat Struct Mol Biol 16:684).
      * **Target.** Buried paired groups, 0.74-0.80 per residue in the native against 0.41-0.47, the largest
        structural contrast measured. The buried-unpaired corner alone differs by at most 0.11 per residue (0.01-0.02
        in 1YOB).
      * **Reach.** The same table answers 8.3.3's third mechanism and is the soluble side of the hybrid's TM1 opening
        (findings 4.7).
      * **Risk.** It mainly widens the native's lead, but it also rewards the globules' non-native ladders, so it is
        judged on 2LZM.
      * **Training.** It replaces a trained table, so it co-trains (7.3).
   2. **D2. Contacts that depend on their surroundings.** This is the force-field form of AlphaFold2's triangle
      updates (a pair's energy depends on the third bodies around it), established as AWSEM's density-switched
      contacts (Papoian et al. 2004, PNAS 101:3352; Davtyan et al. 2012).
      * **Form.** Each side-chain pair energy is scaled by a spline of the two residues' `environment_coverage_sc`.
        That coverage uses the rotamer prior, not the BP solution (`code_audit/audit.md` 1a), so it is fixed before
        BP runs and adds no circular dependence.
      * **Target.** The 0.90-0.93 side-chain energy at 0.54-0.72 segregation: the globules' hydrophobic contacts sit
        at coverage 3.4-3.6 against 4.1-4.7.
      * **Effect.** It makes the trap shallower, because a contact made at core density counts more than one made in
        a loose globule.
   3. **D2b. Who covers a residue.**
      * The burial node already reads 400 weights (`number_independent_weights` 400). ff2.1 trains 20 shared ones, and
        untying them by centre class needs no engine change.
      * Backbone desolvation keeps its own 20 weights, since its writer refuses unequal rows.
      * Target, modest: the hydrophobic share of 0.57-0.62 against 0.61-0.70.
   4. **D3. Secondary-structure balance in compact states.** The order by type is right, but the offset toward
      strand/PPII is not.
      * **Forms.** Either a nearest-neighbour term on consecutive residues' helix scores (the `hbond_energy` basin
        sigmoids), a Lifson-Roig-like nucleation and propagation no current term can express; or per-type alpha/beta
        energies through `hbond.h5`'s residue-class offsets, which the engine already reads (up.md 2.2) and no file
        uses.
      * **Prerequisite step, D3a:** a reference for the right offset, taken before training. Either the helicity of
        these proteins' isolated helix peptides, or all-atom compact states (mdCATH high-temperature frames).
   5. **D4. Temperature-dependent burial.**
      * **Form.** A heat-capacity coefficient per type on the burial and pair energies.
      * **Reach.** It bears on the size dependence of T_m and on replica ladders (`lit/review.md` 1), less on the trap
        at fixed T.
      * **Timing.** It belongs with the temperature calibration (section 6). There is no measurement yet, so by 2.2
        rule 1 it stays a candidate until one is made.

**Order.**
* First, L3 re-measures 8.3.4 with ff3.1. The rules below apply to its numbers, not to ff2.1's.
* By L2's branch:
  * globules ahead in free energy: D1 first (the largest contrast), then D2b;
  * native wins but coils do not reach it: D2 first (the trap's depth), then D3;
  * unresolved: no design change is chosen until L2 resolves.
* Either way, D4 goes with the temperature calibration.
* The changes are made one at a time, each against the control, and scored as above. Every training so far has lost
  0.13-0.23 folded fraction on leaving ff2.1 (findings 1.37); 7.3 sets what a co-training change needs first.

### 8.4 Disordered regions

User 10-28: "for upside3, investigate the possibility of enhancing intrinsic disordered region simulation accuracy";
then "record everything in upside3.md. We will implement this later." IDRs are third in the user's order. Literature:
`lit/review.md` (section 12), with the basis of every claim marked.

#### 8.4.1 What a disordered chain gets from the model

* **Training never sees one.**
  * The 456 proteins are folded. NSE contrasts native-started free replicas with the restrained native.
  * DSE (lambda 0.3) contrasts a near-Tm replica's unfolded frames (Rg above 0.67 Rg_coldest + 0.33 Rg_hottest) with a
    SARW replica at the ladder's hottest temperature. The SARW is Upside with the H-bond, burial, H-bond coverage and
    side-chain pair energies scaled to 0 (`training/ConDiv.py` `zero_for_sarw`).
  * The SARW's only excluded volume is the backbone's, and its only sequence dependence is the Ramachandran maps.
  * Peng 2022 calls the DSE a regularizer giving "the minimum interaction energy that can still generate both a stable
    NSE and a disordered DSE".
* **The SARW is the high-denaturant end.**
  * PNt is 51.3 A in water and 62.0 A in 4 M Gdn, which equals 1.9 A x 334^0.6 (Riback 2017).
  * FhuA's plug and reduced RNase A sit 12% and 7% below their 2 M Gdn sizes.
  * The DSE acts only on near-Tm frames, so how far it shapes a disordered chain at 298 K is measured in B0, not
    assumed.
* **Electrostatics are implicit in the trained terms** (user 10-28).
  * There is no explicit electrostatic term and no ionic-strength input. Charge effects are carried by the trained
    terms, among them the side-chain pair table to about 7 A.
  * No explicit electrostatic term is added (user 10-28): it would count the trained terms' electrostatics twice. So
    salt dependence is outside the model, and a charge-dependent error is fixed through the trained terms and the
    objective.
* **Local structure** comes from the NDRD coil-library maps (loop segments of folded proteins; per-map depths trained
  in ff3.1) and the H-bond energies, both fitted on folded proteins.

#### 8.4.2 Prior work with Upside

No paper compares Upside with an IDP's SAXS or NMR data in water.
* **Jumper 2018:** thermally unfolded training proteins against 1.9 A N^0.6, 15% smaller at the heat-capacity peak
  and 10% larger at high temperature.
* **Riback 2019:** the folding model by REMD at 280-320 K on PNt and five other unfolded proteins gives nu 0.4-0.6 with
  temperature, and Rg proportional to Ree. PNt's measured Rg is not compared.
* **Jeon et al. 2026** (PNAS 123:e2518125123): Upside2 REMD frames pooled over T 0.7-1.1 were the best single library
  for chemical-shift reweighting of 413 BMRB IDRs, and carried transient helices that CALVADOS2's lacked. Pooled over
  temperatures, this shows that Upside reaches the conformations, not that its weights at 298 K are right.

#### 8.4.3 What made other models right for IDPs (`lit/review.md` section 2)

1. **The protein-solvent balance,** the first-order error in all-atom, Martini and residue-level models alike, and
   almost always over-compaction.
   * In the a99SB-disp benchmark, force fields changed on the water side had mean IDP Rg errors of 9-13%, against
     22-36% for the rest and 6% for a99SB-disp.
   * a03ws scales protein-water LJ by 1.10, and Martini 3 scales protein-water by 1.10.
   * In Upside this balance sits in the burial, backbone desolvation and side-chain pair attraction.
2. **Local backbone statistics, which couple to collapse.**
   * Martini3-IDP cut its Rg error from 1.06 to 0.39 nm through bonded terms alone.
   * In a99SB-disp and ff19SB, helicity and compaction compensate each other.
   * C36m removed spurious alpha_L.
3. **Residue specificity,** such as Arg against Lys and Tyr against Phe (HPS-Urry, Mpipi, CALVADOS's per-residue
   lambda).
4. **Debye-Huckel at the experiment's ionic strength,** in every residue-level model with good Rg. This is not open
   to Upside (8.4.1).

How folded proteins were protected:
* Residue-level models hold folded domains by restraints or rigid bodies. That is not open to Upside, which folds de
  novo with one energy.
* a99SB-disp's water change melted helices and folded proteins until the backbone carbonyl-amide LJ and the torsions
  were refit.
* SOP-IDP needs folded-protein energy scales 1.7-4.6 times its IDP ones.

So an IDR change is trained in the same objective as the native terms, and judged on the native panels as well.

Fitting to ensemble data has two precedents:
* **FB-HPS** (Dannenhoffer-Lafage and Best 2021) fitted 20 lambda by gradient descent on d<Rg>/dlambda =
  -(<Rg dU/dlambda> - <Rg><dU/dlambda>)/kT, with an L2 prior toward its start. That is I1's estimator.
* **CALVADOS** fits by Bayesian reweighting. It reuses frames while their effective fraction stays above 60%, under a
  prior over 70 hydrophobicity scales.

#### 8.4.4 B0, the baseline

Mac, `scratchpad/idp_20261028/`, launched 10-28 11:56.

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
* **Runs.**
  * Constant T = T_K / 350.588, the FF2 mapping (6.3), from a random-dihedral chain.
  * Pivot MC every 5 tu, dt 0.009, 30,000 tu per seed, t < 5,000 dropped.
  * Batch A is ff2.1 x 6 seeds and the SARW x 3; batch B is ff3.1 x 6.
  * The two A1-LCD* entries share one set of runs.
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
  * **Convergence** follows findings 10.37 (2.3).
* **Routing.**
  * An offset with good sequence dependence points to the protein-solvent balance: I1, and I2 if the error leans
    toward the SARW.
  * Wrong pairs or a poor r point to residue specificity in the trained terms: I1 with variants in its training set.
  * Wrong helicity points to I3.

**B0 result: ff2.1, ff3.1 and the SARW** (batch A finished 10-28 16:10, batch B 19:37; `read_b0.py` applies the rules
above, output `read_b0.txt`).
* **Mechanism.** Under both force fields some disordered chains spend long stretches as collapsed globules, with Rg
  within 1.3 times 2.2 A N^0.38. Fraction of frames, ff2.1 / ff3.1:
  * ACTR, at 278 K: 77% / 68%.
  * The A1-LCD variants: 32-61% / 35-50%, except +7K+12D at 5% / 6%.
  * PNt: 16% / 23%. CoRNID: 2% / 18%. Alpha-synuclein: 8% / 5%.
  * p27, PTMA, Ash1 and tau: at most 1%.

  The switch between globule and coil is slow. For these entries the Rg autocorrelation time is 730-2,250 tu, against
  50-180 tu for the SARW. So their seeds disagree: under ff2.1, PNt seeds span 23.9-51.6 A and alpha-synuclein seeds
  20.0-33.7 A.
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
  * Cv78/Cv14: 1.056 +- 0.023 (ff2.1) and 1.076 +- 0.045 (ff3.1) against 0.753. Both are on the wrong side of 1: the
    charge-segregated variant stays as expanded as the mixed one.
  * The A1-LCD pairs carry model SE of 0.07-0.12 against measured differences of 6-13%, so this sampling cannot
    resolve them.
    * -12F+12Y/WT-NLS falls within 2 SE under both force fields.
    * WT+NLS/WT-NLS is on the wrong side under ff2.1 (1.081 +- 0.094) and within 2 SE under ff3.1 (0.900 +- 0.088).
    * +7K+12D/WT+NLS is within 2 SE under ff2.1 (1.284 +- 0.120) and 2.2 SE high under ff3.1 (1.368 +- 0.107).
  * The SARW reproduces no pair, which is expected from a model without side-chain pair energies.
* **SARW placement.**
  * For A1-LCD, PNt, alpha-synuclein and tau, both force fields sit below experiment and the SARW above it.
  * For Cv78, both sit between them, and Cv14 is near experiment in both.
  * Where the SARW matches experiment (Hst5, ACTR, Ash1, PTMA), both are smaller.

  So the force-field error does not lean toward the SARW.
* **ACTR helix 1** (1044-1054 in Iesmantavicius 2014's text numbering; their Fig. 1 labels the same residues one
  higher). The wild type's secondary structure propensity, read from their Fig. 1b, is 0.39 at 278 K.
  * ff2.1: 0.33 +- 0.10, with seeds 0.00-0.73.
  * ff3.1: 0.25 +- 0.08, with seeds 0.00-0.50.
  * SARW: 0.05.

  Each collapsed chain holds its own helix content, so helix 1 is unresolved under both force fields. Batch A's
  reading of "close to measured" came from the seed mean alone and is withdrawn. The other ACTR helices are not
  compared.
* **Convergence.** Unresolved by the rule (seed range above the miss): 11 of 16 ff2.1 entries (n_eff 53-580) and 10 of
  16 ff3.1 entries (n_eff 39-413). Every SARW entry is resolved except ACTR and Ash1, which sit on experiment.
* **Routing.** Both force fields make the same sequence-dependent error: ACTR and the A1-LCD variants collapse, and the
  charge-segregated Cv78 does not compact. By the rules this points to residue specificity in the trained terms, that
  is, I1 with variants such as A1-LCD and p27 in its training set. I2 is not indicated. I3 is undecided, since helix 1
  is unresolved.
* **Open.**
  * ACTR is the only entry at 278 K (T 0.793), so its collapse may come from the temperature mapping rather than the
    sequence (B1).
  * The A1-LCD pairs, ACTR's helicity and the unresolved entries need longer sampling or REMD before they can be read
    (B2).
  * ACTR's 26.3 A is CALVADOS 2's value. The local copy of Kjaergaard 2010 holds only the abstract, so its SAXS source
    is unchecked (`lit/review.md` 3.2).

#### 8.4.5 Candidates

I1 and I2 change what training fits; I3 changes the local backbone terms, sharing D3's forms. Each follows 2.2's rules.

* **I1. An IDR size term in ConDiv's contrast.**
  * The loss is sum_k (<Rg>_k - Rg_exp,k)^2 / sigma_k^2, with gradient d<Rg>/da = -(<Rg dV/da> - <Rg><dV/da>) / T,
    from the per-frame dV/da that `compute_divergence` already returns.
  * Each training IDP needs its own sampling, since its Rg decorrelates over 10^2-10^3 tu and one step's 8,000 tu is
    too short. The sampling is reused by reweighting while the effective fraction holds (CALVADOS's 60%), and rerun
    when it falls.
  * An L2 prior holds the parameters toward the start (FB-HPS).
  * The training IDPs span composition, patterning and length, since FB-HPS found that a missing class biases the fit,
    and they are kept apart from B0's test entries.
  * The native terms stay in the same contrast.
  * Intrachain distances (PRE) can be added to the loss later.
* **I2. The DSE target.**
  * The SARW pulls near-Tm unfolded frames toward high-denaturant dimensions. It is also the Sosnick lab's picture of
    an expanded denatured state with minimal residual structure, and it serves HDX (Peng 2022).
  * B0 decides whether it is a suspect. If ff2.1 sits at SAXS in water, the SARW is not hurting disordered chains. If
    disordered chains come out expanded toward the SARW, the target is weighed against I1.
  * A target from measurements in water is the alternative. The DSE is an ensemble contrast and SAXS gives one number
    per protein, so that form is still to be designed.
  * B0's answer: neither case holds. Both force fields sit on the compact side of experiment, away from the SARW, so
    the SARW target is not what puts disordered chains off.
* **I3. Local backbone balance.**
  * Helicity of IDRs comes from the map depths, which ff3.1 trains, and from the H-bond energies. Its design is shared
    with D3 (8.3.5).
  * B0's answer: ff3.1's map depths leave dimensions where ff2.1 had them, except CoRNID and Hst5. ACTR helix 1 is
    unresolved, because the seeds of the collapsed chain range 0.00-0.73.
  * Targets: ACTR's helix-1 series, eight mutants at 23-69% helicity and the wild type at 0.39, at 278 K (5 C)
    (Iesmantavicius 2014). Reading it needs ACTR expanded and sampled first (B1, B2).
  * Also: per-residue helix and strand populations from BMRB chemical shifts, reweighted as in Jeon 2026.
* **Not taken:** an explicit electrostatic term (user 10-28, double counting); restraints on folded parts.

#### 8.4.6 Steps

Recorded 10-28; the user will take them up in a later session.

- [x] B0: the baseline above, read and recorded (10-28 11:56-19:37, MacBook Pro).
- [ ] B1. ACTR's temperature (Mac, about 30 min including the re-analysis). Run ACTR at 298 K under ff2.1 and ff3.1
  with 6 seeds each, and read it against the 278 K runs.
  * If ACTR opens at 298 K, the collapse comes from the T_K / 350.588 mapping, which reaches every IDP measured below
    298 K.
  * If it stays collapsed, the cause is the sequence.

  The test compares the model with itself, so it needs no measured Rg at 298 K. The new entry is kept out of
  `read_b0.py`'s summary set.
- [ ] B2. Sampling, before I1, whose targets need the same sampling. The seed spread falls as the square root of run
  length, so 4x longer runs (120,000 tu) halve it.
  * Running the whole panel that long takes about 13 h per force field on the Mac. PNt would still be marginal, with a
    seed range of 23-28 A against a 12-14 A miss.
  * REMD on midway2 is the alternative for the long chains.
  * The choice is the user's (11, item 10).
- [ ] B0 on the base (S4), once it is chosen.
- [ ] I1-I3: designs and the training IDP set, after B0 and in the user's order (after HDX and large proteins). B0
  points to I1.

How to rerun or extend B0 is in section 12.

### 8.5 Other candidates

* **Parameters held frozen today:** H-bond entries 4-11, and backbone desolvation's centre and sharpness.
  * They have exact derivatives since A29-A30 and are held frozen by the trainer (`hbr` at learning rate 0; `bbenvc`,
    `bbenvs`), so training them needs no new term.
  * No measured defect is assigned to them yet (2.2 rule 1).
* **The glycine context term** (`architecture.md`; proposed, not approved, not started).
  * Glycine's antisymmetric map part would depend on the neighbours' phi/psi.
  * Its trigger has not fired: TM4's helical glycines still flipping, or natively left-handed glycines losing
    alpha_L, under every round-4 candidate. ff3.0 holds TM4.
  * `architecture.md` predates Upside 3 (10-07; it still requires master parity).
* **The side-chain model's own parameters,** which follow the choice of 7.1.
* **The implicit membrane model.** A173 left two properties of the trained membrane tables for an Upside 3 membrane
  retraining that no plan yet contains:
  * the hard channel switch at +-half_thickness;
  * the end slopes of `cb_energy` and `hb_energy`.

  The channel potential's decision (4.5) belongs with it. So does a lateral-pressure term, if one is wanted after
  A105 deleted `membranelateral_potential` (a force with no energy): it would be a new energy, the integral of
  P(z) A(z) dz over the protein's cross-section (plan.md Phase 24 item 5.1).

## 9. Validation, calibration and release (S6, S7)

### 9.1 Validation

The native panels must not lose and the de novo benchmark must gain. Read against released ff2.1, the release
candidate is tested on every set of 2.3:

* native panels (44-domain, mdCATH), by 2.2's paired rule against released ff2.1, not only against the last
  control;
* L1's de novo benchmark and the Peng benchmark;
* H3's HDX test set by H5's rule, and the IDR panel by B0's rules;
* glpG, which runs on the protein force field: KE, C-N, DSSP per TM helix, TM4;
* the GPU ladder's V3 on the release candidate, if the GPU exists by then (10.2).

### 9.2 Calibration

Upside 3's own c (6.3 item 2, 6.4 item 2). Its `T_ref_K` goes into the force field's parameter directory. The MARTINI
tables are rebuilt if c moved.

### 9.3 Refusing ff2.1

* **The user's words.** 10-19: "Later in the stage we should make it refuse to run ff2.1 (not yet ...)". 10-22: "keep
  ff2.1 useable for now". The records say nothing of how.
* **When.** At release, after S6. Until then H3, B0's comparisons and S1's validation run ff2.1 on Upside 3.
* **Constraint.** ff3.1 and ff3.2 share five files with ff2.1, and every force field's hybrid reads
  `ff_2.1/martini.h5` (A234). So a refusal by file content would refuse them too.
* **How (proposal, 11 item 11).**
  * Each Upside 3 force-field directory carries the units record of 5.3 item 1 (`T_ref_K`), which ff2.1's directory
    does not.
  * The config writer requires it and copies it into `/input/units`, and the engine refuses a config without
    `/input/units`.
  * The check is by directory, so the shared files are not caught.

## 10. GPU (G)

### 10.1 Measured (`gpu_report.md`)

Per force evaluation, ff2.1 native starts, one Apple M3 Max core, the 10-19 source; 3 evaluations per step.

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

* **No single hot spot.** The protein model has none, so a speed-up needs nearly every node on the GPU. The hybrid is
  different: in 1rkl in DOPC (3,812 particles) the dry-MARTINI pair loop is 86% of 17.5 ms per step, and every Upside
  protein node together 1.5%.
* **One system cannot fill a GPU.** An evaluation is about 150 kernel launches at 3-5 us each, about the CPU's whole
  0.75 ms for glpG. So the GPU must run many systems per launch: replica-exchange ladders, or a test panel's runs
  together (the 44-domain panel's 176).
* **Already suited:** the graph is static, all arithmetic is float32, one OpenMP thread runs per system, and the RNG
  is counter-based Threefry, which runs under CUDA.
* **Not suited:** BP's data-dependent iteration count, its per-node products of edge messages, and per-step
  allocations and map lookups in several nodes.

### 10.2 The standard: the current force fields, with the CPU's results

User 10-29: the GPU "should produce the same results as the current ffs, without retraining or training another
version of ff". So:
* It runs the force fields the CPU engine runs (ff2.1 and ff3.1 now, 1.4) from the same parameter files. No force
  field is trained for it, retrained on it, or adjusted to it.
* Its reference is the audited CPU engine. The pre-audit engine's results differ from it by the audit's fixes, which
  the GPU does not reproduce.
* It computes energies and forces. Training stays on the CPU, so the trainer's parameter derivatives (about 1/1000
  of the simulation's cost) are not ported.

"The same results" means what it means between today's CPU builds, which already run the force field on numerically
different backends without retraining:
* the Mac and cluster builds differ in their reciprocal square-root refinement, and differed in the finite-math flag
  until A28;
* between el7 and el8, forces are bitwise equal but trajectories differ from step 6.

So trajectories cannot match step for step. Each evaluation matches to rounding, and the ensembles match within
sampling error. A GPU is a third backend under that standard. It uses:
* the same spline coefficient arrays, fitted on the host and uploaded, never refitted;
* BP's same initialisation, damping, tolerance, stop test and synchronous schedule, each system stopping at its own
  test (another schedule can reach another fixed point, a model change);
* float32 kernels with deterministic fixed-point force accumulation;
* the same Threefry streams.

The tolerance has a precedent. The user ruled that the BP stop-test fix needs no retraining. That fix changed
static-frame energies by a median of 7.5e-4 E_up, forces by 8e-5 of their RMS, the ConDiv gradient by 0.1%, and REMD
T_mid by at most 0.002. The GPU should sit well inside that. The open build items U1 and U2 (4.5) are of this kind.

Verification ladder:
* **V0.** Every parameter array, table, index and RNG key on the device is bitwise equal to the host's.
* **V1. Per node, on static frames** (about 20 proteins x 10 frames, including misfolded globules and hybrid frames):
  energy within 1e-5 relative; force RMS within 1e-5; BP compared at a tight tolerance and at the production one.
* **V2. Short dynamics:** KE/(1.5 kT), drift, C-N bond and angle distributions, per-term means; lipid diffusion and
  area per lipid in the hybrid.
* **V3. Ensembles,** on ff2.1 and ff3.1: REMD melting on the Peng set (the `bp_validation` harness), Q and Rg
  distributions, the mdCATH panel, glpG TM4.

### 10.3 Route

| route | effort (estimate) | risk |
|---|---|---|
| a. a CUDA copy beside the CPU code | 4-8 person-months | high: two implementations drift apart |
| b. one implementation of each potential, compiled for CPU and GPU | 6-10 person-months | medium: CPU parity has to be re-proved |
| c. many systems per kernel launch | 2-3 more person-months | medium; needed by any route for a speed-up |
| d. OpenMM plugin, or a JAX/PyTorch rewrite | 4-12 person-months | high: OpenMM has no BP, no derived coordinates and no batching, and its refitted splines would be another force field (10.2); a rewrite is a second implementation |
| d2. GPU for the dry-MARTINI pair loop only | 1-2 person-months | low; the hybrid only, at most about 7x |

**Recommended: b with c.** `InteractionGraph<IType>` is already an evaluator pattern. Templating the evaluators on the
scalar type gives one implementation for both backends, and removes today's split between the SIMD `compute_edge`
and the scalar parameter derivative.

1. [ ] G1: refactor the CPU engine into the single-source form, bitwise equal to the audited CPU engine, with the
   audit harnesses as its test. It waits on S1 (done) and the units option: A+ leaves the engine unchanged, while C
   changes about 42 sites first.
2. [ ] G2, prototype: side-chain pair edges, environment coverage and H-bond coverage on the GPU, batched. Go/no-go:
   one GPU evaluates a batch of 336 systems of ConDiv's protein size (24 proteins x 14 replicas) in about 2 ms, which
   matches 336 broadwl cores. Production supplies batches of that size: a panel's runs together (the 44-domain
   panel's 176, mdCATH's 356) or several REMD ladders.
3. [ ] G3, the remaining nodes of the current force fields.
4. [ ] G4, the BP kernels, which every current force field needs.
5. [ ] G5, the ladder V0-V3.

* **Where development happens.** The CPU engine stays the reference, and development happens on the Mac against the
  host backend (the Mac's GPU cannot run CUDA).
* **d2** is a separate option, worth taking only if the hybrid is in Upside 3's scope (1.3).
* **Terms added later.** A term that S5 accepts is ported when it is accepted, and passes the same ladder on the force
  field that carries it. The GPU does not wait for S5.

**BP is the hardest part of the work.**
* With BP, the rotamer node is about half the protein time and holds the three hardest kernels. It is also the only
  node whose CPU-GPU difference can exceed rounding, because rounding can move the iteration at which the stop test
  fires; V1 compares it at two tolerances for that reason.
* Every current force field uses BP, so G4 is needed whichever side-chain model S3 chooses. A deterministic side chain
  is an ordinary pair term; it would save about a quarter to a third of the GPU effort only for a force field built
  on it, and that is no reason to choose it: it folds worse (7.1).

### 10.4 Hardware

* **midway2, measured.**
  * Our account can use `gpu2`: 5 nodes of 4 Tesla K80; QOS of 4 GPUs, 36 h, 10 jobs; billing counts CPUs only.
  * The K80 is Kepler, which CUDA 12 dropped, so it needs CUDA 11.x (`cuda/11.5` is the default module).
  * The V100 node has been down since 2024-06. There are no Kokkos or OpenMM modules.
  * Fit for correctness testing, not for speed.
* **midway3, from RCC documentation, not measured.** The `gpu` partition has 1 node of 4 A100, 5 of 4 V100 and 5 of 4
  RTX6000.
* **beagle3, from RCC documentation.** 22 nodes of 4 A100 and 22 of 4 A40. RCC's documentation says beagle3 uses no
  SUs, while CLAUDE.md says it bills the group's GPU allocation. Ask before any use.

## 11. Decisions for the user

**Answered.**
* 10-19: Upside 3 is a new release, with no fix measured against ff2.1 (1.1). martini-dev is Upside 3, fixes go only
  there, and master is not changed. ff2.1 keeps running for development, and Upside 3 later refuses it. The audit
  starts at once (plan.md Phase 24).
  * So the restraint-pairing defect needed no ruling of its own: it was fixed in Upside 3 like every other confirmed
    defect (A25).
* 10-21 and 10-22: the audit's own decisions are in plan.md Phase 24 item 5.
* 10-28: HDX first, then large proteins, then IDRs (1.2). No explicit electrostatic term (8.4.1).
* 10-29: the GPU produces the current force fields' results, with no retraining and no force field of its own
  (10.2).

**Open, by stage.**

| # | stage | decision | section |
|---|---|---|---|
| 1 | S1 | the debug validation's gate (proposal in section 3) | 3, 4.4 |
| 2 | S2 | units: A+ (recommended) or C; and whether `martini.h5` moves to native units, with conversion at config time, which would change the CLAUDE.md unit contract | 5.2 |
| 3 | S3 | the side-chain model: BP (folds better, 7.1) or deterministic | 7.1 |
| 4 | S3 | the starting parameters (ff2.1, ff3.1, a Phase 33 result), and which native target the Upside 3 trainer keeps | 7.2, 5.4 |
| 5 | S3 | the gate for changes that co-train the trained tables, and what follows if the control loses | 7.3 |
| 6 | S5 | the accumulation and acceptance rules (proposals) | 2.2 |
| 7 | S5 | the large-protein program as a whole, proposed to the user 10-19 and not yet approved (`remote_jobs.md`) | 8.3 |
| 8 | S5 | L4's form: the replay buffer (proposed) or the gap term | 8.3.5 |
| 9 | S5 | HDX: the training mapping (6.3); how Ferrari et al.'s domain-level data enter (8.2); the base for H3 follows item 4 | 6.3, 8.2 |
| 10 | S4, S5 | IDRs: whether to design I1 and when; Mac time on B1, and B2's form (longer runs or REMD) | 8.4 |
| 11 | S7 | the refusal design (proposal) | 9.3 |
| 12 | G | the route (b with c, prototype first); which GPUs (midway3 `gpu` under which account; beagle3 once its billing is clear); whether d2 is in scope | 10 |
| 13 | | the items of the running Upside 2 campaign that await a NO GUARDS ruling | 4.5 |
| 14 | | the hybrid's scope: a product goal, or a validation case | 1.3 |
| 15 | | the channel membrane potential: keep, smooth, or retire | 4.5 |

## 12. Records and files

`scratchpad/upside3_planning_20261019/` (MacBook Pro; gitignored): the planning reports for goals 1, 2 and 5.

| path | what |
|---|---|
| `audit_report.md`, `audit_work/` | defect inventory, audit units, writer/reader registry (`registry.txt`), warning, cppcheck and sanitizer logs, the parameter-derivative prototype (`harness/param_fd2.py`), the finite-math probe (`probe/fm.cpp`) |
| `units_report.md`, `units_work/` | every unit site with file:line, the parameter-file table, options and verification |
| `gpu_report.md`, `gpu_work/` | profile scripts and logs, replica scaling, the hybrid profile, `web_notes.md` with the sources of every documentation claim |

`scratchpad/upside3_audit/` (MacBook Pro; gitignored): fixtures, baselines, build trees, reviewer outputs (`pass2/`),
and the closing checks (`closing/`); plan.md Phase 24 "Resume here" lists them.

`scratchpad/hdx_term_20261028/` (MacBook Pro; H1):

| path | what |
|---|---|
| `fetch.sh`, `run.sh`, `build_inputs.sh`, `meta.py`, `flags.sh`, `analyze.py`, `fig.py` | copy of the arms from midway2; driver; HDX topologies from `1.config.py`; per-frame energy, walker and CA-RMSD; the stock estimator per replica; MBAR readout with its rules in the docstring; figure |
| `results/results.txt`, `results/ubq_dG_measured.png`, `results/summary.npz` | the readout of findings 5.7 |
| `data/` (4.4 G), `results/*__top_*` (2.0 G) | the three arms (`/beagle3/trsosnic/yinhan/ff3_benchmark/runs/`), copied by `fetch.sh`; flags and raw estimator arrays, regenerated by `run.sh` |

Copies of the scripts and the readout: `/project/trsosnic/yinhan/checks/handoff_20261026/scratchpad/hdx_term_20261028/`.

**Running H2.** Replace the arm and force-field names in `fetch.sh`, `build_inputs.sh`, `flags.sh`, `run.sh` and
`analyze.py` (`ARMS`, `DENOVO`). The HDX topology takes `parameters/ff_3.1`. Then run `fetch.sh` and `run.sh`.

`scratchpad/idp_20261028/` (MacBook Pro; B0):

| path | what |
|---|---|
| `lit/review.md`, `lit/panel.fasta`, `lit/make_panel.py`, `lit/dl*/` | literature review with the basis of each claim; the panel's sequences and SAXS values; its builder (`python make_panel.py lit`); the full texts, the CALVADOS 2 tables and Bremer 2022's SI it read |
| `run_idp.py`, `panel_run.py`, `stat_ineff.py`, `read_b0.py` | one run from a sequence (ff directory or the ConDiv SARW); the panel driver with its analysis (two batches); the Rg statistical inefficiency of the models `results.txt` has analysed; B0's reading rules |
| `panel.fasta`, `results.txt`, `results.npz`, `stat_ineff.txt`, `read_b0.txt`, `panel.log` | the panel as run; B0's readout and per-residue arrays; the reading by the rules; the driver's log |
| `runs/`, `pilot/` | trajectories (one directory per sequence); the alpha-synuclein pilot (`conv_pilot.py`, `pilot.txt`) |

**Resume here (B1, B2 or a rerun of B0, from a cold session).**
* **Where the data are.** Everything is in `scratchpad/idp_20261028/` on the MacBook Pro only. It is gitignored, with
  2.6 GB of trajectories. A session on another computer copies it first, or reruns B0, which takes about 7.7 h on 11
  workers.
* **Environment.** Run `source .venv/bin/activate && source source.sh` from the repo root. `panel_run.py`,
  `stat_ineff.py` and `read_b0.py` need `UPSIDE_HOME`.
* **Adding an entry.**
  * Add it to `panel.fasta`; the header is `name|Rg_A|err_A|T_K|ionic_strength_M|reference`.
  * `panel_run.py` skips finished runs, so the rerun simulates only the new entry's seeds, then re-analyses every
    entry and rewrites `results.txt`.
  * Launch it detached with `scratchpad/detach.sh LOG python panel_run.py`. It takes 11 workers, and the Mac has 10
    performance cores.
  * Wait with one background waiter for the `DONE` line.
* **Longer runs.** Changing `DURATION` marks every run unfinished, and `run_idp.py` has no continuation, so every run
  restarts from a random coil. Longer runs therefore go in their own directory.
* **After any batch.** Run `python stat_ineff.py` and then `python read_b0.py`, and record the result in 8.4.
  `stat_ineff.py` covers only the models on `results.txt`'s `# seeds:` line, because files of a running batch are
  locked.

`scratchpad/large_protein_20261019/` (the Mac Studio; 76 G with trajectories):

| path | what |
|---|---|
| `sim/landscape.py`, `sim/common.py` | native, bulk and vectorial runs with quench |
| `sim/analyze_landscape.py`; `analysis/{ss_check,hb_check,terms}.py` | Q, Rg, quench; helix, local and nonlocal Q; H-bonds; energy by term |
| `analysis_T0{80,88,95}.txt`, `analysis/*_T0{80,88,95}.txt` | outputs |
| `runs_T080`, `runs_T088`, `runs_T095` | trajectories |
| `envfix/` | the 3.13 fix's builds (pre and fixed), `coverage_test.py`, `compare.py`, `fd_test.py`, `dyn_check.py`, `input_diff.py` |
| `lit/review.md`, `code_audit/audit.md` | literature review, code audit |
| `analysis/{buried_unpaired,burial_composition,longrange_burial,electrostatics,helix_to_strand,collapse_vs_helix}.py`, `analysis/*_T080.txt` | the term-by-term review (8.3.4); the first two need `analysis/eng_1019/`, the engine binding that matches the 10-19 `obj/` build |

The channel-potential measurement of 4.5 (10-29, Mac Studio session scratchpad, not kept in the repo):
`channel_steps.py` runs the fixture and evaluates the channel terms per frame; `channel_simpson.py` gives the
force-less part by Simpson's rule against `rama_map_pot`.
