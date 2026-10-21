# Upside 3: folding large proteins de novo

**Status: proposal, not started.** Upside 3 starts after the side-chain verdict, deterministic (no-BP) against
probabilistic (BP) side chains at their converged epochs (plan.md Phase 21 item 7, findings 1.21), and is built on
the model that wins. The measurement is in (10-20, findings 1.21): at both fixed points BP folds 0.353 against
no-BP's 0.241 on 89 held-out mdCATH domains (p 1e-5); the choice is the user's. Everything measured below is ff_2.1, the probabilistic side chain.

Investigation: plan.md Phase 22, 2026-10-19. Data, scripts, the literature review (`lit/review.md`) and the code
audit (`code_audit/audit.md`) are in `scratchpad/large_protein_20261019/` (gitignored, 76 G with trajectories).

## 1. The problem

* Upside folds large proteins de novo poorly (user, 10-19).
* Jumper 2018 (PLoS CB e1006578): all 20 test proteins under 100 residues come within 5 A of the native from an
  unfolded start. The CASP11 targets of 115-178 residues do not: T0855-D1 (115) never found its fold. The authors
  write that ConDiv's short runs "optimiz[e] only over fluctuations in or near the native well".
* Peng 2022 (FF2, JCTC 18:550): reversible folding only for the 40- and 43-residue mini-proteins; ubiquitin (76)
  did not refold within a CPU-week because of long-lived misfolded states.
* Nothing above 150 residues has been trained or benchmarked de novo: the 456 training proteins are 50-150 residues
  (median 114), the 44-domain panel 50-75, the Peng de novo set at most 92.

## 2. What was measured

### 2.1 The native-start size trend is the metric

On the 89-domain mdCATH native-start panel, ff_2.1_ctrl38 folds 0.42, 0.44 and 0.24 of frames in the 50-100,
101-150 and 151-250 bins. "Folded" is Q at or above the all-atom 5th percentile, a cutoff that tightens with N
because the all-atom Q spread narrows. Upside's mean Q sits about 0.13 below all-atom in every bin and its Rg at
1.01 of native, and with a fixed-offset cutoff the trend disappears. Upside holds large natives; the problem is
reaching them.

### 2.2 De novo test

T4 lysozyme 2LZM (164 residues), RNase HI 2RN2 (155) and A. vinelandii flavodoxin II 1YOB (179, apo chain), all
three refolded in vitro without a ribosome, and the engrailed homeodomain 1ENH (54) as the control. ff_2.1
(released), 60,000 tu per run, Mac Studio. Starts:

* native, with the panel's anneal (T 0.05 to T over t 96-400);
* bulk: a random-dihedral chain, same anneal;
* vectorial: 30 residues annealed for 2,000 tu, then 5 residues added every 1,000 tu at T, each new segment PPII and
  turned away from the chain, then the full chain to 60,000 tu.

Every last frame was quenched for 1,000 tu at T 0.3. Q is the CA native-contact fraction (|i-j| >= 4, d < 8 A,
formed below 1.2 d0); nonlocal Q counts only |i-j| > 12. H-bonds were counted from H and O placed on the backbone
(H...O < 2.5 A), which reproduces the crystal counts.

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

### 2.3 What the chains end in

At T 0.80 every de novo run ends in a compact globule with almost no native tertiary contacts: Rg 0.9-1.3 of native
for the large proteins, 1.1-1.7 for the control. The globules are helix-poor: they form 0.36 (T4 lysozyme) to 0.64 (homeodomain) of the
native start's helical H-bonds. T4 lysozyme's carry 1.5-5x the crystal's nonlocal, sheet-type H-bonds. Native starts
keep crystal helix and H-bond content; their Q drifts by up to about 0.2 over 60,000 tu and mostly plateaus.

### 2.4 Energy by term, quenched globule minus quenched native

Bulk / vectorial against native, each the mean over the last 10 frames of every quench (6, 6 and 3 runs), E_up;
the terms sum to the engine total.

| term | 1ENH | 1YOB | 2LZM | 2RN2 |
|---|---|---|---|---|
| H-bond | +34 / +50 | +115 / +118 | +140 / +142 | +160 / +144 |
| rotamer (side-chain free energy) | +23 / +32 | +35 / +39 | +42 / +23 | +24 / +37 |
| side-chain burial | +2 / +3 | +20 / +20 | +13 / +15 | +13 / +15 |
| Rama | -5 / -5 | -7 / -19 | -20 / -27 | -19 / -32 |
| backbone desolvation | -4 / -7 | -17 / -17 | -18 / -19 | -16 / -15 |
| total (bulk / vectorial) | +49 / +73 | +155 / +150 | +179 / +159 | +163 / +149 |

* Most of the gap is backbone H-bonds the globule never forms.
* The globule's side chains recover 0.87 (1ENH), 0.91, 0.90 and 0.93 of the native's side-chain energy (rotamer
  plus burial, bulk), so collapse costs little.
* Rama and backbone desolvation favour the globule, by 24-47 E_up in the large proteins.

### 2.5 The co-translational hypothesis

The user's hypothesis was that folding is co-translational (S. Fried, JHU: about a third of E. coli proteins do not
refold in vitro), so a thermodynamic force field cannot fold large proteins from a denatured chain. On these
proteins:

* N-to-C growth equals bulk collapse in every measure, and is slightly worse on the homeodomain.
* All three large proteins refold in vitro without a ribosome, so Upside fails at something the real chains do on
  their own.

The hypothesis may hold for obligately co-translational folders, but it does not explain these failures. The growth
here has no ribosome tunnel or surface, and its rate (200 tu per residue) is not mapped to real time.

### 2.6 What these runs do not show

* The globules' free energy against the native's: only energies were compared, not entropies.
* A size contrast under one protocol. The control does not fold at constant T in 60,000 tu at any of the three
  temperatures, so the size dependence rests on Jumper 2018 and the user's observation.
* Generality: four proteins, one force field, constant temperature.

## 3. Reading

For these in-vitro refolders Upside's failure is search. RNase HI's native is stable at T 0.88, and its coils do
not find it at 0.80 or 0.88. Upside ranks the native 150-180 E_up below the globules on average and 90-110 below the
best of them, but the chain collapses
on side-chain burial, which already pays 87-93% of the native's side-chain energy, before it forms secondary
structure. It then cannot form native H-bonds inside the globule, and Rama and backbone desolvation reward the
non-H-bonded compact backbone.

Mechanisms, ranked by the evidence (literature in `lit/review.md`, code in `code_audit/audit.md`):

1. **No negative design against compact misfolds.** ConDiv contrasts a native-restrained replica with 12 free
   replicas started at the native (8,000 tu, T 0.8-1.1) and an expanded self-avoiding-walk ensemble. Compact
   non-native globules never enter the contrast, and their number grows with N (Bryngelson and Wolynes 1987).
   Decoy-based training is how AWSEM (T_f/T_g, Davtyan 2012) and OPEP (energy gap, Sterpone 2014) handle this.
2. **The training range and weighting.** Training stops at 150 residues. Each protein's contrast is summed, not
   normalised per residue (`ConDiv.py:403-405`), so 110-150-residue proteins dominate. The fixed 12-replica ladder
   loses acceptance with N (0.63 / 0.55 / 0.49 across the training bins; 0.40 predicted at 200). Airas and Zhang
   2025 found that adding larger structures to training is what fixed size transfer.
3. **Backbone desolvation acts as a surface term.** It rewards exposed NH and CO (scale -0.30). The native's
   fraction of exposed groups falls from 0.27 (50-100 residues) to 0.19 (251-400), while almost every group of an
   unfolded chain is exposed, so the term's pull toward unfolded grows with N. Here it also favours the globules by 15-19 E_up.
4. **Weak cooperativity.** Jumper 2018 reports a transition "not fully cooperative"; additive pair terms fail the
   calorimetric two-state test (Chan 2000).

Ruled out or not size-dependent:

* explicit size terms or neighbour-list capacities (none exist);
* BP convergence (14-28 iterations from 76 to 271 residues);
* the burial window (the same coverage range at every size; residues above it rise from 5.3% to 10%);
* the pivot sampler (loses efficiency with N but cannot bias).

## 4. Engine defect fixed on 2026-10-19

`environment_coverage_hb` (the side-chain coverage of backbone H and O that feeds backbone desolvation) labelled
each covering side chain with the type of the H/O site whose index equals that side chain's residue number. The
writer stored one type per site (`py/upside_config.py`), and the engine read `aa_types[id2]` with `id2` the residue
(`src/environment.cpp`).

* Only 7-9% of covering side chains got their own residue's type. The coverage weights run from 0.56 (ASP) to 2.12
  (TRP), so each NH/CO's burial was a sum over scrambled weights.
* The same code is in master, and the defect does not depend on size.
* **Fix.** The writer stores one type per residue. The engine sizes `aa_types` from its dataset and requires every
  `id2` to index it. The side-chain burial node already held per-residue types and is unchanged for every allowed
  config (residue exclusions are blocked by the parser).
* **Old files run as master.** Files written before the fix keep their per-site arrays and run exactly as master: all
  energies, terms, coverages and forces are bitwise equal on 8 configs at 3 frames each. An engine without the fix
  refuses a new file with a size error, so a new config cannot run silently on an old cluster binary.
* **Verified.** Both builds match a numpy reference of the coverage to 5e-7, each with the types it should use. The
  finite-difference derivative of the energy change matches the change in force (ratio 1.0001 at h = 0.005, 40
  coordinates on two globules).
* **Effect on ff_2.1 is small.** Backbone desolvation changes by -0.1 to +0.2 E_up at the natives and +0.05 to +1.3 at
  the globules, with no other term moving. So the defect is not a cause of the large-protein failure.
* **Native-start dynamics.** With the fixed engine on new configs at T 0.80 (t 2,500-5,000 tu), KE/(1.5 kT) is
  1.007-1.010 and the peptide C-N bond 1.33-1.34 A as before. Q, Rg and potential fall within single-run noise of the
  overnight runs (Q 0.83 / 0.69 / 0.81 / 0.74 against 0.84 / 0.67 / 0.82 / 0.74 for 1ENH, 2LZM, 2RN2, 1YOB).
* **Where it is installed.** In the repo and the local `obj/` (backups `*.bak_pre_envfix_20261019`). It is not
  deployed to `$P` or `/beagle3`, because running trainings write their configs with the old writer. ff_2.1 and
  every training so far ran with the defect, and Upside 3 trains with it fixed.

Also found by the audit and not fixed (the user's call): `make_restraint_group` (`py/advanced_config.py:560`)
means to skip same-residue restraint pairs, but `x/3 != y/3` is float division in Python 3 and skips only
self-pairs. It is the same in master and affects only ConDiv's native-restrained replica.

## 5. Plan

Gate: the side-chain verdict (above).

**Phase A: baseline on the chosen side chain and the fixed engine.**

1. [ ] Deploy the fixed engine to the Upside 3 training tree by Phase 21's procedure: a fresh midway2 build, bitwise
   parity on pre-fix configs, then install.
2. [ ] Fix the de novo benchmark: 8-12 single-domain proteins of 50-250 residues with published in-vitro refolding,
   spread over the alpha, beta and alpha/beta classes (the four here included). Score by nonlocal Q and a
   size-normalised RMSD (RMSD_100) or TM-score, from coil, with replica exchange.
3. [ ] Replica exchange on RNase HI and the homeodomain: 16 replicas over T 0.80-1.00, coil and native starts.
   About 5 h and 15 h on the Mac, or about 27 h on midway2. This answers 2.6's two open questions:
   * if the globules beat the native in free energy, the energy ranking must change;
   * if the native wins but coils never reach it, the funnel shape is the target.
4. [ ] Repeat 2.2-2.4 under the chosen side-chain model. The globule's side-chain energy fraction and the Rama and
   desolvation preference for globules are Phase B's targets.

**Phase B: training changes, one at a time, each scored on A2 and on the native panels (44-domain, mdCATH).**

1. [ ] **Compact misfolds as negative examples.** Make decoys by collapsing each training protein from coil at
   T 0.80 (the bulk protocol, then quench), regenerated from the current force field each epoch. Two forms:
   * start some free replicas from decoys, which keeps ConDiv's objective and fixes the sampling bias Jumper 2018
     names;
   * add an explicit native-against-decoy gap term next to the SARW term.

   The first is the smaller change. Targets: the native-globule gap per residue, the globule's side-chain energy
   fraction (now 0.87-0.93), and Rama plus desolvation favouring globules (now 24-47 E_up in the large proteins).
2. [ ] **Train past 150 residues.** Add mdCATH domains of 150-250 residues, kept apart from the panel's held-out
   domains. Normalise each protein's contrast per residue, and space large proteins' replica ladders by
   1/sqrt(N) (the audit's 6a and 6c).
3. [ ] Only if B1-B2 leave a gap: the desolvation term's surface balance, and a cooperativity term.

**Phase C: validation.** The native panels must not lose, and the de novo benchmark must gain. Re-check the
hybrid (glpG), which runs on the protein force field.

**Not planned: a co-translational protocol.** It did not help here, and the large test proteins fold in vitro
without one.

Open decisions for the user:

* the side-chain model (the gate);
* B1's form;
* the restraint-pairing defect;
* whether ff_2.1 is re-validated on the fixed engine.

## 6. Files (`scratchpad/large_protein_20261019/`)

| path | what |
|---|---|
| `sim/landscape.py`, `sim/common.py` | native, bulk and vectorial runs with quench |
| `sim/analyze_landscape.py`; `analysis/{ss_check,hb_check,terms}.py` | Q, Rg, quench; helix, local and nonlocal Q; H-bonds; energy by term |
| `analysis_T0{80,88,95}.txt`, `analysis/*_T080.txt`, `analysis/*_T088.txt`, `analysis/*_T095.txt` | outputs |
| `runs_T080`, `runs_T088`, `runs_T095` | trajectories |
| `envfix/` | the fix's builds (pre and fixed), `coverage_test.py`, `compare.py`, `fd_test.py`, `dyn_check.py`, `input_diff.py` |
| `lit/review.md`, `code_audit/audit.md` | literature review, code audit |
