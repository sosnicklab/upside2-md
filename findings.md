# Findings

Consolidated and reorganised by subject on 2026-09-05. This file had grown into an append-only log of 84
numbered updates spanning 2026-07-17 to 2026-09-03, much of it development narrative for code that has
since been rewritten or replaced. That narrative has been removed. What remains is organised by subject
rather than by date and is meant to read as one document: the standing rules and the measurements that
justify them, how the hybrid is put together, the defects whose root causes are established, the questions
still open, the diagnostic procedures worth reusing, and a closing list of claims that turned out to be
wrong. Where a finding is cited elsewhere in the repo by its old update number, that number is kept inline
(for example "findings 103") so a grep still lands on the right passage.

---

## ff3.0 DOES fix the HDX-relevant observable for glpG: backbone H-bond retention (2026-09-09)

Backbone H-bond count from the run logs, 3 seeds each, same local protocol. This is the quantity
HDX actually reports on, and the one the pre-ff3 failure was recorded against (occupancy 0.844).

| arm | start | mid | end | retained |
|---|---|---|---|---|
| control, ff_2.1 no coverage | 194.6 | 145.6 | 124.5 | **0.640** |
| trained 269 + coverage | 194.6 | 183.5 | 137.0 | 0.704 |
| trained 500 + coverage | 194.6 | 192.0 | 170.1 | **0.874** |
| recorded pre-ff3 REMD failure | | | | 0.844 |

**Step 500 retains 0.874 of the backbone H-bonds against 0.640 for `ff_2.1`, and is the only arm
that clears the 0.844 failure baseline.** It also beats step 269 by a wide margin (0.874 vs 0.704),
which is the opposite of the TM4-helix-fraction ranking and agrees with the core-RMSD ranking.

Caveat, as everywhere in this local test: single-temperature MD, not the REMD the 0.844 came from,
so the comparison to that number is indicative rather than strict. The between-arm ordering under
identical conditions is the solid part.

This matters more than the fold drift for the deliverable. HDX measures H-bond opening, not
absolute tertiary packing, so a bundle that loosens while its H-bonds stay closed can still give
correct protection factors. Of the four observables measured, ff3.0 step 500 is now better than
ff_2.1 on all of them and better than step 269 on three.

## Why the glpG bundle splay is NOT a bug, and what is left to try (2026-09-09)

Two candidate defects were tested and both came back clean, so the splay is thermodynamics, not a
broken table or exclusion.

**The intercalating lipid is acyl tail, not headgroup.** Bead types wedged between two or more TM
helices, enrichment against their abundance in the system:

| role | share in bundle | share in system | enrichment |
|---|---|---|---|
| acyl tail | 83.4% | 69.2% | **1.20x** |
| glycerol | 10.7% | 15.4% | 0.69x |
| headgroup | 5.9% | 15.4% | **0.39x** |

Per bead the most enriched is C4A at 2.41x, a deep tail bead; the most excluded are GL0 at 0.11x
and NH3 at 0.21x. Charged and polar headgroups are being kept out of the hydrophobic interior
exactly as they should be. This is physically correct hydrophobic solvation of the helix surfaces,
not a spline table putting headgroups in the core.

**The coverage terms do not count lipid as burial.** `hbond_coverage` takes
`id2` = 747 sidechain beads and `hbond_coverage_hydrophobe` takes `id1` = 630 backbone atoms with
the same 747, i.e. protein only; no MARTINI bead enters either. So a helix solvated by lipid
correctly reads as *exposed* rather than buried, and the coverage term is not rewarding the
splayed state.

**What is left.** With both mechanisms ruled out, the imbalance stands: protein-lipid attraction is
full-strength dry-MARTINI while protein-protein packing comes from a core force field trained on
456 soluble proteins with zero membrane content, and lipid wins the competition for the same
hydrophobic surfaces. The principled fix is to include membrane proteins in the ConDiv training
targets so packing is calibrated against lipid competition. Scaling SC-env or BB-env down, or
adding an orientational term to hold the bundle, is forbidden and would break the model.

## THE ACTUAL DEFECT: glpG's hybrid config has NO backbone-environment term at all (2026-09-10)

While implementing the MARTINI-type coverage fix I read the coverage code properly and found my
mechanism was misframed, and then found a much simpler and better-supported defect.

### My coverage story was wrong

`HBondCoverage` is built as `CoordNode(get_dset_size(1,grp,"index2")[0], 1)` and accumulates
`output(0, edge_indices2[ne])` (`hbond.cpp:438-460`). Its output is therefore **per sidechain
bead**, and it is a 1-body cost handed to the rotamer solver for *placing a sidechain where it
covers a backbone H-bond*. It is **not** a measure of how shielded a backbone H-bond is from
water. So "make hbond_coverage lipid-aware" does not implement "lipid shields the H-bond", and the
48% collapse I measured at TM4 is a sidechain-placement quantity, not a desolvation one.

### What is actually missing

Comparing the node lists directly:

| | soluble benchmark config | glpG hybrid config |
|---|---|---|
| `bb_sigmoid_coupling_environment` | present | **ABSENT** |
| `environment_coverage_hb` / `_sc` | present | **ABSENT** |
| `hb_environment_coverage_hn` / `_oc` | present | **ABSENT** |
| `cat_pos_bb_coverage`, `weighted_pos`, `sigmoid_coupling_environment` | present | **ABSENT** |

**glpG has no environment-dependent backbone term whatsoever.** Its `hbond_energy` is purely
geometric: a backbone H-bond is worth the same whether buried in the protein core or dangling in
bulk. In the soluble force field `bb_sigmoid_coupling_environment` and the two
`hb_environment_coverage` nodes supply exactly that burial dependence.

The cause is one line: `upside_config.py:2444`, `if args.environment_potential:` gates every one of
those nodes, and the hybrid glpG builder never passes `--environment-potential`. My soluble
benchmark configs DO pass it, which is why they have the nodes and glpG does not.

### Why this explains everything measured

* TM4 has the lowest intrinsic helix propensity of the six helices (23% glycine against TM1's 5%).
  With no burial bonus for its backbone H-bonds it has nothing but intrinsic propensity to hold it,
  so it is the first to melt while TM1 survives.
* Temperature independent across the whole ladder: a **missing energy term**, not a thermal effect.
* Present in all four variants and predating ff3.0: the hybrid builder has always omitted it.
* Retraining the rotamer tables (ff3.0) helps a little but cannot substitute for an absent
  backbone term, which is exactly the partial improvement observed.
* It also explains the long-standing note that "the trained environment table has no home in glpG":
  the environment nodes are not there to receive it.

### The fix, and the caveat on it

Rebuild the glpG seeds passing `--environment-potential parameters/ff_3.0/environment.h5` and
`--bb-environment-potential .../bb_env.dat`. Both are already trained and already shipped. **No new
parameters, no C++, no retraining.**

The one thing to check first: the omission may have been deliberate. `planned_job.md` records
"the environment node would double-count explicit dry-MARTINI lipid". That worry is real for the
implicit *membrane* potential, but the environment coverage counts **protein neighbours only**, so
it does not see lipid and cannot double-count it. That should be verified rather than assumed
before deploying anywhere.

This also puts the user's MARTINI-type idea in the right place: once the environment term exists,
whether lipid should contribute to *environment* coverage is the meaningful question, and the
MARTINI Qa/Qd/C1 classification is the parameter-free way to answer it.

## Carried over from planned_job.md before deleting it (2026-09-10)

That file tracked the ff3.0 training chain, which is finished, deployed and superseded. Everything
in it was either done, recorded elsewhere, or wrong. These two were recorded nowhere else.

**The two Ramachandran libraries, and the mirror that differs between them.**
`parameters/common/rama.dat` is the OLD force field and is asymmetric in GLY *by design*; do not
"fix" it. The GLY-symmetric library ff_3.0 must run against is `parameters/common/rama3.dat`. The
correct GLY mirror is **`m[::-1,::-1]` for `.up` `rama_map_pot` maps but `roll(m[::-1,::-1], 1)` for
the library's own `dimer_pot`**. Using the wrong one fabricates about 3.3 E_up of asymmetry on
perfectly good seeds. Always verify against a chiral control (ALA/SER/HIS must still show 10-11 E_up
of asymmetry, or the mirror hit everything) and against `dG(aR->aL)`, which reads 0.000 on a correct
GLY map.

**One GLY experiment is still untested:** whether GLY asymmetry ever contributed *independently* to
TM4's instability. The control is a de-symmetrized run, a few hours locally. Worth doing only if the
glycine thread is picked up again; note that the 2026-09-10 measurements argue against glycine being
the TM4 cause at all, and that the asymmetric reference state is the live GLY defect instead.

## The MARTINI-typed H-bond correction makes things worse; not deployed (2026-09-10, settled)

Built, tested, refuted. The idea was sound and maps cleanly onto the code: dry-MARTINI's own
particle typing says where water is not, so let the beads that can neither donate nor accept an
H-bond (C1 and C3, the 2511-bead acyl core) count toward the backbone environment coordinate, while
the donor/acceptor headgroups (Qa PO4, Qd NH3, Na GL1/GL2, and P4 GL0, which is MARTINI's own water
bead) do not. It needs no C++: `hbbb_coverage` takes bare positions, `Add` sums coverages
elementwise, and node types resolve by prefix, so the apolar channel enters `bl` through the already
trained `hbond_weight`. Nothing was refit.

Three seeds, same protocol and seeds as the other two arms, scored on the build's real helical
segments:

| arm | TM4a 134-139 | TM4b 141-155 | TM1 30-43 |
|---|---|---|---|
| armB500, no environment | 0.807 | 0.805 | 0.911 |
| + backbone environment | 0.691 | 0.868 | 0.900 |
| + MARTINI apolar channel | **0.464** | **0.717** | **0.810** |

Every measure got worse, by -0.089 to -0.228. The likely reason is a sign trap: the coverage is
non-negative and enters `bl` with a positive weight, and raising `bl` raises the energy, so the
channel acts as an effective **repulsion between backbone and acyl chains**, which is the opposite
of the intended "the core is dry, stop charging for desolvation". Making it act in the intended
direction needs a signed term, which is the same wall the lipid-coverage scan hit.

Not deployed anywhere. Also worth recording: scored on proper segments rather than the old capped
window, the plain backbone-environment term is not neutral either, it trades TM4a (-0.116) against
TM4b (+0.063). Neither route helps.

## Where TM4 actually loses helix, and what it is not (2026-09-10)

Four hypotheses died on measurement this round, and the survivor is a much narrower target.

**TM4 is the least lipid-exposed helix in the protein, not the most.** Backbone beads with a lipid
neighbour inside 12 A: TM1 10.7, TM3 7.0, TM5 6.5, TM6 4.1, TM2 4.4, **TM4 2.0**, and 7 of TM4's 22
backbone beads have no lipid within 12 A at all. The 134-139 stretch sits at the bilayer midplane
(|z| = 0.5-2.1 A) and still contacts nothing: its total dry-MARTINI energy against all 3627 lipid
beads is 0.000. TM4 is packed in the middle of the six-helix bundle. Every protein-lipid explanation
is therefore excluded, which is why the backbone-environment term, the lipid-coverage scan and the
MARTINI-typed channel all returned nothing.

**The MARTINI backbone typing is inert here even though it looks wrong.** dry-MARTINI types TM4's
134-137 as `Nd` (helix N-cap, eps 2.7 kJ/mol against C1) and 131-133 as `P5` (coil, 0.5 kJ/mol)
against `N0`'s 3.5, so those beads are nominally under-stabilised in the acyl core. Re-typing them
all `N0` changes the total by **-0.3 kJ/mol**, because there is no lipid near them to interact with.

**Much of the headline number was a window artifact.** TM4 was being scored over 131-152, but the
build never treats 130-133 (`P5`) or 140 (`C5`) as helix. Scored on the build's own helical
segments, TM4's two pieces sit at 0.807 (134-139) and 0.805 (141-155), against a protein mean near
0.82. TM4 is not an outlier by that measure; the worst segment in the protein is the proline-rich
60-77 at 0.536.

**The real loss is sharp and local.** Per-residue helix occupancy, first 10% of frames vs last 10%,
three seeds: 139 TYR 1.000 -> 0.350, 140 ALA 1.000 -> 0.333, 141 LEU 1.000 -> 0.333, 142 MET
1.000 -> 0.533, and 134/135/137 lose 0.32-0.43, while **143-146 and 149 hold at 1.000** and 153-155
hold. So the unwinding is confined to 134-142, the deeply buried midplane half, and the half nearer
the headgroups is untouched.

**It is not input strain.** In the seed structure TM4's i->i+4 O...N distances run 2.81-3.12 A
continuously from 134 to 149, a pristine helix. Residue 140 in particular is perfectly helical
(O140...N144 = 2.85 A) despite being typed `C5`, so that SS call is a mis-assignment rather than a
real break, and it is energetically inert anyway.

**A separate real defect, found on the way and worth fixing on its own merits.** `rama_map_pot`'s
GLY maps in the built config are exactly symmetric (max asymmetry 0.000000), so that fix holds. But
the shared reference state `rama_map_pot_ref` is a single map applied to all 210 residues and is
chirally asymmetric by 1.07, worth **+0.32 E_up** between the two helical basins. It is added
directly as energy (`*pot += value`, `rama_map_pot.cpp`). For every non-GLY residue the residue map
swamps it (net alphaR favoured by 2.07). For GLY, whose own map is now exactly symmetric, it is the
only chirality term left, so **every glycine is pushed toward alphaL by +0.37 E_up**. That is
measurable: GLY residues sit below their own helix's mean occupancy in every helix. It is the same
class of bug as the GLY symmetrization fix, in the one place that fix did not reach. It is **not**
the TM4 cause: glycine density does not predict which helix fails (r = +0.23, wrong sign, and TM3
carries 17.9% GLY at 0.928 occupancy).

**What is left.** TM4 is the most protein-buried helix (13.55 backbone coverage against a 9.66
protein mean), and burial correlates *positively* with helix stability across segments (r = +0.52),
so it is not a buried-backbone penalty either. That leaves the sidechain and rotamer terms as the
dominant environment for this helix, which is exactly what the M-vs-R arm test varies.

## The backbone-environment term does NOT fix TM4 (2026-09-10, settled)

Tested and refuted. The hybrid glpG config genuinely lacks every environment node while all 32
soluble benchmark configs carry them, so restoring them was a real modelling gap to close. It does
not close the TM4 gap.

**Measurement.** `armB500_s{1234,2345,3456}` (coverage, no environment) against the same seed with
`write_weighted_pos` + `write_bb_environment` injected, identical protocol, identical three seeds,
300000 steps at T=0.70, none diverged, 401 frames each:

| arm | TM4 helix fraction | TM1 |
|---|---|---|
| armB500, no environment | 0.673 [0.543-0.744] | 0.879 |
| + backbone environment  | 0.693 [0.623-0.813] | 0.876 |

**+0.020 with fully overlapping ranges on n=3.** Seed scatter is +-0.10, five times the effect. The
protocol is not insensitive: the same three seeds resolved control 0.441 vs armB269 0.782 on this
exact system, so it detects +0.34 comfortably. It detects nothing here.

**Why, measured rather than argued.** My first explanation was that the burial coordinate counts
protein neighbours only, so a TM helix reads as solvent-exposed. That is wrong, and the measurement
says the opposite. Counting neighbours within the channel's own 6 A radius on the folded seed:

| region | protein | apolar lipid | polar lipid |
|---|---|---|---|
| TM1 29-49 | 10.77 | 0.50 | 0.01 |
| TM4 131-152 | **13.69** | **0.02** | 0.04 |
| whole protein | 9.66 | 0.15 | 0.02 |

`bb_sigmoid_coupling_environment` is `scale * compact_sigmoid(bl - center, sharpness)` with
scale -0.301, center 2.0, sharpness 0.5, and `compact_sigmoid` reversed, so the solvation credit is
fully on below bl = 0 and fully off above bl = 4. TM4 sits at 13.7, the **most buried region of the
protein**, three times past saturation. Its credit is already zero, which is why restoring the term
moved TM4 by 0.020 and why the whole config's energy moved by only 32 E_up. TM4's backbone is
shielded by its own sidechains and by the rest of the six-helix bundle; it barely touches acyl
chains at all (0.02 apolar neighbours, against TM1's 0.50).

So the environment term is not mispricing TM4, it is saturated and inert there. Any fix that works
by feeding the existing burial coordinate is therefore dead on arrival for TM4, including making
lipid count toward burial: the MARTINI-typed channel adds +0.0 to TM4 and +1.0 to TM1 on the folded
seed. Whatever destabilises TM4 is not the backbone environment term.

**Not deployed.** The gate was that it had to fix TM4 locally before going to rockfish or midway2.
It did not, so glpG production was left alone on both clusters. The term is still correct physics
that a hybrid config should carry, and the NP ff3.0 rebuild includes it for that reason, but it is
not the TM4 answer and must not be sold as one.

## Lipid-coverage scan: what was built, and the sign caveat (2026-09-10, running)

Implemented and launched. Three findings about feasibility, each correcting an earlier estimate of
mine, and one caveat that limits what the result can mean.

**It needs no C++.** Three properties of the engine make a lipid coverage channel a config-only
change: node types resolve by **prefix** (`deriv_engine.cpp:598`), so a group named
`hbbb_coverage_lipid` is instantiated as the registered `hbbb_coverage`; `rotamer` takes a
**variable-length** `prob_nodes` list (`rotamer.cpp:1174-1178`) and uses each node's output
**directly as a 1-body energy** (`rotamer.cpp:868-870`); and `hbbb_coverage`'s interaction is
`HbondEnvironmentCoverageInteraction2` with **n_dim2 = 3**, so its second group can be bare
positions. The trained `hbond_coverage` type cannot be used: it needs n_dim2 = 6, a direction
vector, and lipid beads have none. Each prob node must have n_elem = 747, the sidechain bead count
(`rotamer.cpp:702`), which group1 satisfies.

Built as `scratchpad/lipid_coverage_test/inject_lipid_coverage.py` (gitignored): group1 = the 747
sidechain beads, group2 = the 3627 LIPID beads with ids offset by 100000 so the built-in
|id1-id2|<=2 exclusion cannot fire, `interaction_param` (1,1,4) = r0, r_sharpness, dot0=-2 (angular
sigmoid ~1, making it a pure radial burial count), dot_sharpness.

**Sign caveat, and it matters.** This interaction returns a **non-negative** count and it is used
directly as energy, so the channel can only *penalise* lipid proximity, never reward it. Measured
at the seed: baseline E = -24944.285, and with the channel +15.0 (r0=4), +72.2 (r0=6), +183.0
(r0=8), all finite, |deriv|max 310-351. So this is **not** the hypothesised fix, which would
*remove* a spurious desolvation penalty from lipid-solvated H-bonds. What it actually tests is
whether biasing sidechains toward protein burial and away from lipid rescues TM4. That is a
related but distinct proposition, and the parameters are chosen rather than trained.

**Readout is the trend, not any value.** r0 = 4.0, 6.0, 8.0 against the existing `armB500_s1234`
baseline, same protocol and seed. Monotonic TM4 improvement with r0 would support the mechanism and
justify doing it properly with a trained, signed term; flat or non-monotonic refutes this route.

## Coverage test (2026-09-10): the hypothesis survives, narrowed to the HYDROPHOBE term

Ran the cheap precursor to the lipid-aware-coverage experiment. Two things came out of it, one a
correction to my own cost estimate.

### Correction: the fix is NOT a Python edit

I estimated "2 hours, extend the neighbour list in `martini_inject_coverage.py`". Reading the code,
that is wrong. Both coverage nodes take **exactly two argument nodes**,
`["protein_hbond", <sidechain placement node>]`, and `index2` indexes directly into the sidechain
placement node's bead list. Group2 *is* that node's output. Lipid beads cannot be appended from
Python without either giving the node a third input group or putting MARTINI beads into the
sidechain placement node, which would corrupt the rotamer solver. **Making coverage lipid-aware is
an engine change in `src/environment.cpp`, not a script edit.**

### What was measured instead, and what it shows

`upside_engine.get_output()` exposes both coverage nodes directly, so the *premise* is testable
without touching C++: does the coverage associated with TM4 actually collapse, and does it collapse
more than for a helix that stays folded? Mean of 3 seeds, first 6 frames against last 6:

| term / helix | early | late | change |
|---|---|---|---|
| `hbond_coverage` TM4 | 0.383 | 0.536 | +0.153 (rises) |
| `hbond_coverage` TM1 | 0.580 | 0.792 | +0.212 (rises) |
| **`hbond_coverage_hydrophobe` TM4** | 3.948 | 2.047 | **-1.901, a 48.2% collapse** |
| `hbond_coverage_hydrophobe` TM1 | 3.359 | 2.848 | -0.511, 15.2% |

against helix outcomes of TM4 0.97 -> 0.67 and TM1 0.97 -> 0.92 over the same window.

**The plain `hbond_coverage` term rises for both helices and does not discriminate.** The
**hydrophobe** term collapses **3x more for TM4 than TM1**, and that is the one term whose
discrimination matches the helix outcome. It is also exactly the term lipid should contribute to:
acyl tails are hydrophobic and do bury the backbone, yet contribute zero here.

### Caveats, stated because they bound the claim

* `get_output` returns shape **(747, 1)** for both nodes, i.e. one value per *sidechain bead*, and
  both coverage nodes are **arguments to the `rotamer` node**, not direct multipliers on the
  backbone H-bond energy. So the route from this term to backbone helix stability is indirect,
  through the coupled sidechain solver. The correlation is real; the causal chain is not proven.
* An earlier probe of mine reported `hbond_coverage` TM4 dropping 24.6%. That used only column 0 of
  the output and is superseded by the table above, which averages the full output with a consistent
  `id2` mapping for both terms.

### Verdict

The hypothesis is **not refuted and is now sharper**: it is the hydrophobic backbone-burial
coverage, not the sidechain-competition coverage, that fails at TM4. That is worth the engine work.
The next step is a third input group on `HbondEnvironmentCoverageInteraction` in
`src/environment.cpp` so MARTINI beads can contribute hydrophobic coverage, then rerun the
three-arm test. Note the trained tables were fitted with protein-only coverage, so a positive
result is evidence about the mechanism, not a deployable force field.

## RE-ANALYSIS: the real cause of TM4 unfolding (2026-09-10). Supersedes the splay-causes-melt story.

Redone from scratch on the WT production trajectory (13117 frames) and the REMD ladder. Two of my
earlier claims do not survive.

### It is NOT thermal melting

WT glpG, ff_2.1 production, TM cores in the last block of each replica across the ladder:

| replica | T | TM1 30-48 | TM4 134-151 |
|---|---|---|---|
| 0 | 0.700 | 0.943 | **0.432** |
| 6 | 0.742 | 0.910 | 0.521 |
| 12 | 0.786 | 0.840 | 0.534 |
| 18 | 0.831 | 0.712 | 0.363 |
| 24 | 0.877 | 0.818 | 0.494 |
| 27 | 0.900 | 0.918 | **0.531** |

**TM4 is flat across the whole ladder, 0.36-0.57, with no monotonic temperature dependence, and is
no better at the coldest rung than the hottest.** A helix melting thermally would be high at
T=0.700 and low at T=0.900. TM1 meanwhile behaves normally. So the partially melted TM4 is what
this force field prefers at *every* temperature sampled: the native helix is not the free-energy
minimum, and no amount of cooling or sampling fixes that.

### It IS reversible, and it has NOT equilibrated

Same trajectory, 13117 frames at T=0.70: starts at 1.000, and by tenths runs
0.99, 1.00, 0.97, 0.88, 0.82, 0.72, 0.63, 0.47, 0.31, 0.39. Overall mean 0.716, sd 0.278,
min 0.056, max 1.000. After the initial drop, **47.8% of frames are still above 0.8** and it
recovers above 0.8 repeatedly.

So TM4 is not destroyed; it interconverts between helical and partly melted, with the population
drifting toward melted. And it is **still drifting after 13117 frames**, so even the long
production run is not equilibrated.

### What I got wrong: splay does not demonstrably cause the melt

I previously wrote that the bundle opens, TM4 loses packing, and the helix then melts. Testing the
ordering properly does not support it. Cross-correlating the raw time series gave
"lipid gain leads helix loss by 99 time units, r = 0.64", but both series are monotone trends and
that number is an artifact of the shared trend. **On first differences, which remove the trend,
every coupling collapses to r = 0.21-0.26 and the lipid-helix lag falls to +/-9 time units, less
than one sampling interval.** No causal direction is resolvable. The earlier per-residue
correlation was already weak (+0.21), and the melt is not a discrete event either: the largest
single-step drops are -0.17 to -0.39 at unrelated times in different seeds.

### What the cause actually is, and what remains hypothesis

**Established:** the force field gives TM4's native helix insufficient free-energy preference over
a partially melted alternative, at all temperatures in the ladder. Everything else has been
excluded: not GLY asymmetry (verified 0.000000 in every seed, production and live), not the
mutations (TM4 is byte-identical in all four variants and all four decay), not ff3.0 (it predates
it, and ff3.0 only slows it), not the capped measurement window (the decay is real on 134-151),
not temperature, not irreversible damage.

**Hypothesis for why TM4 and not TM1**, still untested: TM4 has the lowest intrinsic helix
propensity of the six helices (5 glycines in 22 residues, 23%, against TM1's 1 in 21) *and* is the
most protein-buried at the start (0.22 lipid beads/residue against TM1's 0.55), so it depends most
on the burial-dependent coverage terms, which count only protein and are blind to lipid. A helix
that is both intrinsically floppy and maximally dependent on a burial term that mis-reads its
environment is the one that gives way first. Testable by making coverage lipid-aware and repeating
the three-arm run.

## The GLY fix is intact; the earlier "TM4 is fixed" conclusion was a convergence error (2026-09-10)

Two separate questions, separate answers.

### The GLY symmetrization is correctly applied and is not the problem

Checked with the correct `.up` mirror `m[::-1,::-1]` on every relevant config: the four production
seeds, the production replica that actually decayed (`glpG-RKRK-79HIS.run.0.up`), the two live
ff3.0 arm-test seeds, and the local three-arm configs. **All PASS: 23 GLY maps, max asymmetry
0.000000**, with the chiral control (SER/HIS/ALA) at 11.10-11.15 E_up, exactly the expected 10-11.

So GLY symmetry is real and present, and TM4 still melts. This does not contradict the earlier
record, which already said GLY symmetry was **necessary but not sufficient** (ff_2.1 with symmetric
GLY still gave TM4 0.441). It was never the claimed fix.

### What was actually wrong: means compared across arms that had not equilibrated

The claimed fix was trained pair + coverage, ARM B at 0.782 against a 0.441 control. Both are
**trajectory means**. Broken into quintiles, TM4 core 134-151, 3 seeds:

| arm | q1 | q2 | q3 | q4 | q5 | mean | 2nd-half slope /1000 t.u. |
|---|---|---|---|---|---|---|---|
| control ff_2.1 | 0.87 | 0.51 | 0.43 | 0.43 | 0.41 | 0.528 | **-0.022** |
| trained 269 + coverage | 0.99 | 0.97 | 0.87 | 0.84 | 0.81 | 0.898 | -0.055 |
| trained 500 + coverage | 0.97 | 0.88 | 0.85 | 0.72 | 0.67 | 0.817 | **-0.135** |

**The control had already bottomed out by q3 and is flat thereafter. The trained arms were still
falling, the step-500 arm six times faster than the control.** Comparing means, or any fixed-time
value, between a converged arm and two non-converged ones systematically flatters the
non-converged ones: part of their higher score is simply "has not finished decaying yet".

The independent check that this is the right diagnosis: my local ff_2.1 control ends at **0.41** and
the ff_2.1 **production** run, 25-40x longer, ends at **0.43**. ff_2.1 had genuinely converged
inside the short run. No production-length ff3.0 run exists yet, so its endpoint is unmeasured, and
extrapolating the -0.135 slope from 0.67 reaches the control's 0.41 within roughly another 1900
time units.

**What survives:** the trained tables plus coverage are better than `ff_2.1` at every time point,
and that improvement is real. **What does not survive:** "TM4 is fixed", or any absolute claim that
it clears a threshold. On this evidence ff3.0 *delays* TM4 loss rather than preventing it.

### The tooling bakes the error in

`final_analysis.py` and `decide_arm.py` both compute helix fraction as `.mean()` over every frame,
with no convergence or trend check. **The arm test running tonight inherits this**: it will rank
arm M against arm R on trajectory means of runs that are almost certainly still decaying. For an
M-vs-R *relative* comparison that is tolerable if both decay alike, but no absolute "TM4 passed"
should be read out of it. Report the quintile trend and the second-half slope alongside any mean.

## TM4 melts in ALL FOUR glpG variants, including the wild type, and predates ff3.0 (2026-09-10)

Asked whether the TM4 unfolding is variant-specific, on the expectation that wild-type 79HIS at
least should hold. It is not variant-specific, and the wild type is not spared.

**First: TM4 is byte-identical in all four variants.** The mutations are at position 79 (HIS/ALA)
and 115 (SER/THR); TM4 is residues 134-151, `LeuThrGlyValValTyrAlaLeuMetGlyTyrValTrpLeuArgGlyGluArg`
in every one. There is no sequence basis for one variant's TM4 behaving differently.

**Second, from the ff_2.1 PRODUCTION runs that produced the pre-ff3 baseline** (28-replica REMD,
replica run.0 = T=0.70 rung, first output block against last), core windows TM1 30-48 and TM4
134-151:

| variant | blocks | TM1 first | TM4 first | **TM4 last** | TM1 last |
|---|---|---|---|---|---|
| glpG-RKRK-79HIS (wild type) | 51 | 1.000 | 1.000 | **0.432** | 0.943 |
| glpG-RKRK-79HIS_S115T | 32 | 1.000 | 1.000 | **0.350** | 0.771 |
| glpG-RKRK-79ALA | 32 | 1.000 | 1.000 | **0.601** | 0.556 |
| glpG-RKRK-79ALA_S115T | 32 | 1.000 | 1.000 | **0.368** | 0.828 |

**Every variant starts at TM4 = 1.000 and every one decays to 0.35-0.60.** The wild type ends at
0.432, no better than the mutants; 79ALA is in fact the *best* at 0.601. TM1 largely holds
(0.77-0.94) except in 79ALA.

**This predates ff3.0 entirely.** These are the ff_2.1 runs from before the retraining, so the TM4
melt is a property of the glpG/dry-MARTINI hybrid model rather than of the retrained force field or
of any mutation. ff3.0 slows it (WT local test: 0.97 -> 0.67 against ff_2.1's 0.87 -> 0.41) but the
phenomenon is the same one, and the expectation that wild-type TM4 should be stable has never been
met by this model.

Caveat on comparing the numbers directly: the production runs are 28-replica REMD and the local
three-arm tests are single-temperature MD, so absolute values are not interchangeable. The
within-table comparison across the four variants is like for like.

A matching ff3.0 single-temperature test of the other three variants was launched 2026-09-10 to
confirm the ordering carries over; the WT arm of that test is the one already reported.

## Why TM4 still partially unfolds under ff3.0 (diagnosed 2026-09-10, from the VTF)

User inspected the ff3.0 trajectory and saw TM4 better than ff_2.1 but still partly unfolded. That
is correct, and the trajectory mean (0.817) hides it: TM4 starts essentially perfect and decays.

**It is a progressive melt in the MIDDLE of TM4, not end-fraying.** Per residue, first fifth of the
run against the last fifth, mean of 3 seeds, with the i->i+4 backbone H-bond distance at the end:

| res | early | late | i->i+4 O...N late |
|---|---|---|---|
| 139 TYR | 0.992 | **0.432** | 6.44 A |
| 140 ALA | 0.992 | **0.333** | 6.36 A |
| 141 LEU | 0.996 | **0.342** | 5.22 A |
| 142 MET | 0.996 | 0.543 | 4.32 A |

A formed helix H-bond is ~3.0 A, and 146->150 is still 3.57 A, so the break is local to 139-142
with a second weak patch at 134-137. **The seed has 134-151 fully helical**, so none of this is
inherited from the starting structure.

**The cause is the bundle opening, not anything specific to TM4's sequence.** TM4 begins
**completely protein-buried, 0.0 lipid beads per residue**, and ends with lipid at every residue
while losing **40% of its packing contacts** to the other five helices (10.2 -> 6.1 per residue).

Same process in every arm; ff3.0 slows it without stopping it:

| arm | TM4 helix early -> late | lipid beads early -> late | packing early -> late |
|---|---|---|---|
| control ff_2.1 | 0.87 -> 0.41 | 2 -> 31 | 41 -> 31 |
| trained 269 + coverage | 0.99 -> 0.81 | 2 -> 22 | 42 -> 34 |
| trained 500 + coverage | 0.97 -> 0.67 | 1 -> 20 | 44 -> 36 |

ff3.0 keeps more packing and admits ~35% less lipid than ff_2.1, which is the visible improvement,
but the mechanism is untouched.

**Honest limit on the causal claim.** Across the 18 core residues the per-residue correlation
between helix loss and packing loss is only **+0.21**, and with lipid gain **+0.09**. So the link is
global, TM4 as a whole loses support and melts somewhere in the middle, and is NOT a residue-level
"this contact broke so this turn opened".

### The specific candidate defect: the coverage term is blind to lipid

`hbond_coverage` and `hbond_coverage_hydrophobe` take `id2` = 747 **protein sidechain beads** and
`id1` = 630 protein backbone atoms. No MARTINI bead enters either (verified in the config). Those
terms modulate backbone H-bond strength by burial.

In a soluble protein, trained on 456 soluble proteins, "not covered by protein" always means
"exposed to water", where the backbone H-bond competes with water and should be weakened. In a
bilayer it can instead mean "surrounded by acyl tails", which shield an H-bond at least as well as
protein does. The model cannot tell those apart, so as TM4 goes from protein-buried to
lipid-solvated its coverage collapses and its H-bonds are penalised as if they had been dunked in
water. That is exactly when and where the helix melts.

It also predicts the observed ordering: TM4 is the most buried helix at the start (0.22 lipid
beads/residue against TM1's 0.55) so it suffers the largest coverage change, and it is TM4 rather
than TM1 that degrades (TM4 0.817 vs TM1 0.919).

**This is a hypothesis consistent with the data, not a proven cause.** It is directly testable:
include MARTINI beads in the coverage calculation and rerun the same three-arm test. If TM4 holds,
the mechanism is confirmed. Note that would be a model-structure change requiring retraining and
validation, not a parameter tweak, and must not be done by scaling SC-env or BB-env.

## The remaining glpG problem is NOT TM4: the helix bundle splays and lipid wedges into it (2026-09-09)

TM4's helix is fine. Decomposing the residual core CA-RMSD on the 9 local trajectories shows the
defect is tertiary, not secondary.

**Individual helices are intact; the bundle is not.** CA-RMSD of each helical segment fitted to
itself, against the whole-core value:

| arm | whole core | mean per-helix | worst helix |
|---|---|---|---|
| control | 7.47 | 2.96 | 5.40 |
| trained 269 | 5.33 | 2.07 | 2.98 |
| trained 500 | **4.89** | **1.71** | 2.38 (TM4) |

Every helix holds its own shape to 1.2-2.4 A while the assembly is 4.9 A out. Dropping the two
short peripheral helices (19-26, 50-57) changes almost nothing (4.96 -> 4.77), so this is the TM
core repacking, not floppy termini.

**The bundle splays laterally and flattens.** Against the seed (mean inter-helix centroid distance
16.19 A, Rg_xy 11.16, Rg_z 9.39):

| arm | inter-helix dist | Rg in-plane | Rg along normal |
|---|---|---|---|
| control | +4.36 | +4.34 | -1.62 |
| trained 269 | +1.72 | +2.50 | -1.26 |
| trained 500 | +2.32 | +2.42 | -1.16 |

The helices spread apart in the membrane plane and the bundle pancakes along the normal.

**Lipid wedges into the bundle.** Dry-MARTINI beads within 7 A of two or more different TM helices,
i.e. sitting in inter-helical space:

| | seed | control | trained 269 | trained 500 |
|---|---|---|---|---|
| intercalated lipid beads | **5** | 46.2 | 41.0 | **36.1** |

A 7-fold invasion. ff3.0 cuts it by ~25% against the control but does not come close to restoring
crystal packing.

**Nothing has equilibrated except the lateral splay.** An earlier version of this note said the
step-500 core RMSD had "plateaued". That was wrong: it came from an arbitrary quintile threshold.
Fitting the actual slope over the second half of each run, per 1000 time units, mean +/- sd over
the 3 seeds:

| metric | control | trained 269 | trained 500 |
|---|---|---|---|
| whole-core RMSD | +0.440 +/- 0.545 | +0.309 +/- 0.077 | **+0.123 +/- 0.068** |
| per-helix RMSD | +0.384 +/- 0.164 | +0.526 +/- 0.208 | +0.279 +/- 0.205 |
| Rg in-plane | +0.278 +/- 0.390 | **-0.033 +/- 0.258** | **-0.034 +/- 0.156** |
| intercalated lipid | +6.150 +/- 3.578 | +5.500 +/- 4.604 | +5.751 +/- 2.909 |

Only **Rg in-plane** is genuinely flat, and only for the two trained arms: the lateral splaying
does stop, at about +2.4 A of spread. Everything else is still moving at 2700 time units.

The step-500 core RMSD rises at +0.123 A per 1000 time units, roughly 2.5x slower than step 269 but
**not zero**. Its per-helix RMSD is also still climbing. So step 500 slows the fold degradation, it
does not halt it, and 4.9 A is a snapshot of an ongoing process rather than an endpoint.

**The most important number here: ff3.0 does not slow lipid invasion at all.** The intercalation
rate is statistically identical across all three arms, +5.5 to +6.2 beads per 1000 time units,
including the untrained control. The trained force field starts the invasion later and so ends with
a lower count (36 against 46), but the *rate* is untouched. Lipid is still entering the bundle at
the end of every run, with no sign of saturating.

This also corrects the ranking claim made from TM4 helix fraction: step 500 is better than step 269
on fold drift rate, but neither is stable, and longer sampling should be expected to make both
worse. The 12 h REMD arms sample far longer than 2700 time units.

**Ruled out: that TM4 is preferentially attacked.** TM4 is the *least* lipid-exposed helix, in the
seed (0.22 beads/residue against TM1's 0.55) and after the run (+1.00, the smallest increase of the
six). Correlation between helix glycine content and lipid-contact increase is **negative**,
r = -0.41. TM4 is a buried, central, glycine-rich helix with small side chains: it barely touches
lipid itself, and depends almost entirely on packing against its neighbours. That is why it reads
as fragile when the bundle loosens, even though nothing is attacking it directly.

### Hypothesis for the cause, NOT yet tested

An asymmetry in the hybrid model. `exclude_intra_protein_martini = 1`, so dry-MARTINI supplies
**zero** intra-protein attraction and all helix-helix packing comes from the Upside core force
field, which is trained on **458 soluble proteins with zero membrane content**. Protein-lipid
attraction, meanwhile, is full-strength dry-MARTINI. Lipid therefore competes for the same
hydrophobic surfaces with an interaction that was never balanced against the protein-protein term,
and wins, prying the bundle open. Consistent with training reporting median RMSD ~0.93 A on its
soluble set while this membrane protein settles 4.9 A out.

This is a hypothesis. It cannot be tested by scaling SC-env or BB-env down, which the project
forbids and which would break the physical model. A legitimate test would compare the same protein
against a different surrounding phase, or measure whether the splay tracks protein-lipid contact
energy frame by frame.

## Is TM4 resolved? NOT YET. Helix fraction says yes, fold integrity does not (2026-09-09)

Measured on the same 9 local trajectories, against the metrics the failure was actually recorded
on rather than helix fraction alone.

| arm | TM4 core 134-151 | TM1 core 30-48 | helical-core CA-RMSD | GLY49 phi | GLY133 phi |
|---|---|---|---|---|---|
| control | 0.528 | 0.832 | 7.62 A | +41.3 | +18.2 |
| trained 269 + coverage | 0.898 | 0.867 | 5.66 A | +93.9 | -18.6 |
| trained 500 + coverage | 0.817 | 0.919 | **5.07 A** | +42.5 | +31.8 |
| recorded pre-ff3 REMD failure | - | - | **4.15 A** | +75.9 | +67.7 |

**The two proxies disagree, and they disagree about which force field is better.** Step 500 beats
step 269 on core RMSD (5.07 vs 5.66), TM1 (0.919 vs 0.867) and Rg (20.35 vs 19.48 against a 20.4 A
crystal), and loses only on TM4 helix fraction (0.817 vs 0.898). Picking a winner on TM4 alone
would invert the ranking given by the other three.

**Core CA-RMSD does not clear the bar.** Every arm, including the best, sits above the 4.15 A that
was recorded as the *failure*. Two caveats that stop this being conclusive: this is
single-temperature MD at T=0.70 for 2700 time units, not the 28-replica REMD the 4.15 A came from,
and the core here is the 162 residues helical in the seed, which may be a larger and looser
selection than whatever the 4.15 A was measured over. **The two numbers may not be commensurable;
re-measure with a matched core definition before treating the comparison as real.** What is solid
is the internal contrast: ff3.0 cuts core drift from 7.62 to 5.07 A.

**The GLY phi pass criterion is not met by any arm, and is itself suspect.**
`BASELINE_TM_pre_ff3.txt` requires GLY49 and GLY133 median phi negative (alphaR). No arm achieves
both. But **the crystal-derived seed itself has GLY49 at +94.1 and GLY133 at +141.6**, i.e. it
fails the criterion at t=0, before any dynamics. Both are helix caps (section above), and a cap
glycine in a left-handed or extended conformation is normal. This is the same class of error as the
TM4 window: a criterion that the native starting structure does not satisfy. **Do not certify or
condemn a force field on this criterion until it is re-derived from the reference structure.**

### Verdict

Not resolved, and specifically:

* **Resolved:** the old force field's genuine TM4 weakness (core helix 0.528) is fixed, to 0.817-0.898
  against TM1's 0.867-0.919, i.e. near parity. And most of the apparent "TM4 is half of TM1" gap is
  the capped measurement window, now diagnosed.
* **Not established:** that the *fold* is fixed. Core CA-RMSD is the metric the deliverable depends
  on and no arm demonstrably beats the failure baseline on it.
* **Not usable as evidence:** the GLY phi criterion, until re-derived.
* **Thin:** n=3, one temperature, one variant, local MD. The 28-replica REMD arm tests are the real
  measurement and had not reported when this was written.

## What actually makes TM4 look unstable (diagnosed 2026-09-09)

Asked what destabilises TM4 besides GLY symmetry. Diagnosed on the 9 local trajectories. The
dominant cause is **the definition of the metric**, not the physics, and two plausible physical
causes were tested and ruled out.

### 1. The TM4 window includes a 3-residue non-helical cap. This is most of the effect.

`TM4 = (131, 152)` in `tm_health.py`, `decide_arm.py` and `final_analysis.py`. Dihedrals of the
**crystal-derived seed itself**, before any dynamics:

| res | aa | phi | psi | helical? |
|---|---|---|---|---|
| 131 | PHE | -134.9 | +161.5 | no |
| 132 | GLY | -160.4 | +147.6 | no |
| 133 | GLY | +141.6 | +146.7 | no |
| 134 | LEU | -76.2 | -0.6 | **yes, the helix starts here** |
| 151 | ARG | -89.3 | -25.2 | yes |
| 152 | ASP | -140.8 | +73.6 | no |

So 4 of the 22 residues in the window are loop/cap in the starting structure and cannot be
helical. **The metric is capped at 18/22 = 0.818**, against a stated pass criterion of >0.8. A
perfect TM4 helix scores 0.818. `tm_health.py`'s own comment says "GLY133 at N-cap", so the cap was
known and included anyway.

TM1 has the same problem, but far milder: only GLY49 is a cap, ceiling 20/21 = 0.952. **That
asymmetry alone manufactures much of the "TM4 is half of TM1" impression.**

Measured over the true helix instead:

| arm | TM4 131-152 | TM4 **134-151** | TM1 29-49 | TM1 **30-48** |
|---|---|---|---|---|
| control, ff_2.1 no coverage | 0.441 | 0.528 | 0.800 | 0.832 |
| trained 269 + coverage | 0.782 | **0.898** | 0.832 | 0.867 |
| trained 500 + coverage | 0.673 | **0.817** | 0.879 | 0.919 |

As a fraction of the achievable ceiling, armB269 reaches **96%** and armB500 **82%**. armB269 was
essentially at the maximum the metric permits.

### 2. Ruled out: buried charge

TM4 ends ARG148, GLY149, GLU150, ARG151, ASP152, and ARG148/GLU150 sit 11.1 A from the bilayer
centre, inside the acyl region (headgroup planes at +19.0/-21.1 A, thickness 40.1 A). But this is
**not special to TM4**: 10 charged residues are buried below 12 A across the protein, only 2 of them
in TM4. And they are paired, not naked: ARG148-GLU150 4.1 A, ARG151-ASP152 3.4 A. Not a strain.

### 3. Ruled out: glycine content as an internal helix breaker

TM4 is glycine-rich, 5 of 22 (23%) at 132, 133, 136, 143, 149, against 1 of 21 (5%) for TM1, which
is a real sequence difference and the obvious suspect. But within the true helix core the glycines
are **not** the weak positions:

| arm | helix fraction at GLY | at non-GLY | deficit |
|---|---|---|---|
| control | 0.486 | 0.536 | -0.050 |
| armB269 | 0.854 | 0.907 | -0.053 |
| armB500 | 0.911 | 0.799 | **+0.113** |

Under the final force field the glycines are the *better* positions. The two glycines that never go
helical, 132 and 133, are the cap from section 1. So glycine content makes TM4 *susceptible*, and
GLY rama symmetry was needed to stop the alphaL bias, but internal glycine is not what limits TM4
now.

### What is left

Under the old force field TM4 was genuinely weak: core 0.528 against TM1's 0.832. The trained
tables plus coverage nodes fix that, to 0.898 (step 269) and 0.817 (step 500) against TM1's 0.867
and 0.919, i.e. near parity. What remains is the step-269 -> step-500 regression measured above,
which is force-field drift, not a TM4-specific defect.

### Consequences

* **Report TM4 on 134-151, not 131-152**, or report the ceiling alongside it. As written, the
  >0.8 criterion demands ~100% helicity of every core residue.
* **`decide_arm.py` scores the arms on the same capped window**, so both cluster arms are being
  judged against a 0.818 ceiling. That does not bias the M-vs-R *comparison*, since both arms share
  it, but any absolute "TM4 passed/failed" read from the arm test is measured against the wrong
  scale.
* Do not "fix" this by widening the helix criterion. The window is what is wrong, and the seed's own
  dihedrals say where the helix starts.

## Local TM4 check of the FINAL step-500 force field (2026-09-09): NOT resolved, and worse than step 269

Three arms, three matched RNG seeds, the identical protocol used for the 2026-09-07 four-arm run
(same seed file `6a8285d1...`, dt 0.009, T=0.70, `--disable-recentering`, 300000 steps, 400 frames).
All 9 runs finished clean, 0 non-finite frames.

| arm | diverged | TM4 helix | TM1 helix | Rg mean | max C-N |
|---|---|---|---|---|---|
| CONTROL, ff_2.1, no coverage | 0/3 | 0.441 [0.298-0.633] | 0.800 | 20.15 | 4.16 |
| ARM B269, trained step-269 + coverage | 0/3 | **0.782** [0.657-0.863] | 0.832 | 19.48 | 11.16 |
| ARM B500, trained step-500 + coverage | 0/3 | **0.673** [0.543-0.744] | 0.879 | 20.35 | 4.06 |

**The two anchors reproduced the 2026-09-07 reference exactly**, to three decimals including the
replicate ranges, so the binary rebuild of 2026-09-07 21:34 (the `MonteCarloSampler` virtual
destructor) changed nothing and ARM B500 is directly comparable.

**Training from step 269 to 500 made TM4 worse, consistently.** Paired on matched seeds:

| observable | per-seed change 269 -> 500 | mean | t(2) |
|---|---|---|---|
| TM4 | -0.114, -0.131, -0.082 | **-0.109** | -7.6 |
| TM1 | +0.040, +0.063, +0.039 | +0.047 | +6.0 |
| Rg | +0.17, +1.11, +1.31 | +0.86 A | +2.5 |

Every seed moves the same way in every observable, so this is a real effect at n=3, not scatter.

**It is a trade, not a pure regression.** Step 500 is *better* on everything except TM4: TM1 rises
0.832 -> 0.879, the over-compaction flagged in `planned_job.md` is gone (Rg 19.48 -> 20.35 against a
20.4 A crystal, essentially exact), and the worst peptide C-N bond drops 11.16 -> 4.06 A, i.e. the
one stretched bond in ARM B269 is gone.

**This is the random walk showing up in an observable.** The training objective is flat between
step 269 and 500, so the parameters diffuse; TM4 and TM1 then diffuse in *opposite* directions.
There is no monotone "more training is better", which is the same conclusion the parameter-drift
analysis reached, now confirmed downstream. It is also the strongest argument yet for averaging over
the plateau rather than picking whatever step training happened to stop at.

**Do NOT read 0.673 against the >0.8 REMD criterion directly.** That criterion comes from
`BASELINE_TM_pre_ff3.txt`, which is 28-replica REMD; this is single-temperature MD at T=0.70 and is
harsher, as the control shows (0.441 here against 0.645 for 79HIS at the T=0.70 rung of the REMD
baseline). The trustworthy result is the *paired, relative* one: step 500 is 0.109 below step 269 on
TM4 under identical conditions. Whether either clears 0.8 is what the 28-replica arm tests measure.

## ff3.0 convergence: the objective converged, the parameters never will (measured 2026-09-09)

**The parameter vector performs a random walk, not a descent.** Measured on midway2's own run by
comparing `pair_interaction` at step 496 against earlier steps of the same run:

| lag (steps) | rel_rms drift | rel_rms / sqrt(lag) |
|---|---|---|
| 12 | 0.0486 | 0.0140 |
| 25 | 0.0682 | 0.0136 |
| 50 | 0.0943 | 0.0133 |
| 100 | 0.1444 | 0.0144 |
| 200 | 0.2241 | 0.0158 |

`rel_rms / sqrt(lag)` is constant to within 15% over a 16-fold range of lag. Drift therefore grows
as sqrt(steps) with no fixed point: **training longer does not converge the parameters, it just
diffuses further.** The mild rise at lag 200 is a weak systematic component on top of the diffusion.

**The training objective, by contrast, has plateaued.** Median RMSD across rockfish's 272 logged
minibatches, by quarter: 0.9337, 0.9296, 0.9190, 0.9322 (col 1, noise +/- 0.035) and 2.188, 2.125,
2.085, 2.152 (col 2, noise +/- 0.23). The slope over the last half is +0.0084 and +0.1266 per 100
minibatches, i.e. zero within noise and if anything slightly *worse*. Fit quality saturated well
before step 500.

**Consequence: "ff3.0" is a sample from a distribution, not a unique answer.** Compared at the
*same* step 496, midway2's and rockfish's force fields differ by mean `rel_rms = 0.2115` across the
three trained tables (pair 0.2118, coverage 0.1963, hydrophobe 0.2264), while the untrained
`hydrophobe_placement` and `rotamer_center_fixed` agree to 1e-16.

**The two runs are NOT independent replicates, and that spread is a LOWER bound (corrected
2026-09-09 by the user).** rockfish did not train from scratch: it resumed from midway2's
half-trained checkpoint, so the runs share all history up to the branch at **step ~274** and have
diverged only over the 222 steps since. Measured from that branch, on `pair_interaction`:

| step | steps since branch | rel_rms | rel/sqrt(steps since branch) |
|---|---|---|---|
| 275 | 1 | 0.0757 | 0.0757 |
| 300 | 26 | 0.1410 | 0.0276 |
| 338 | 64 | 0.1729 | 0.0216 |
| 400 | 126 | 0.2136 | 0.0190 |
| 496 | 222 | 0.2627 | 0.0176 |

Against the within-run drift of one run over the same lags: rel/sqrt(lag) = 0.0169, 0.0170, 0.0183,
0.0200. **The two runs separate at the same diffusive rate as a single run wanders**: after the
branch they are two independent walks, and nothing more.

Extrapolating that rate to a full independent training (500 steps rather than 222) gives
`0.0176 * sqrt(500) = 0.39`, so two force fields trained independently from `ff_2.1` should differ
by roughly **40%**, not 26%. The measured spread is small only because the runs share 274 steps.

(The one anomaly: at step 275, a single update after the branch, they already differ by 0.0757,
far above `sqrt(1) * 0.0176`. The first update after a resume is outsized, most likely because the
Nesterov/momentum state is not carried across the resume. Worth knowing before reading any
single-step comparison near a restart.)

**ff3.0 loads and runs, but its force tail is heavier (measured 2026-09-09 on rockfish).** The
step-500 rockfish force field was patched into a real glpG seed (coverage nodes injected, then the
trained tables) and put through `check_hybrid_up.py --require`. It passes: hybrid interface intact,
intra-protein MARTINI still excluded, production stage, finite energy. Same seed, three configs:

| config | nodes | energy | median \|f\| | p90 | p99 | p99.9 | max | atoms \|f\|>50 |
|---|---|---|---|---|---|---|---|---|
| unpatched seed | 25 | -25076.96 | 0.045 | 8.29 | 23.1 | 41.5 | 65.1 | 2 |
| coverage + ff_2.1 | 28 | -25025.68 | 0.045 | 8.52 | 27.8 | 48.8 | 111.2 | 5 |
| coverage + ff3.0_rf | 28 | -24957.68 | 0.045 | 9.14 | 31.6 | 97.4 | 318.3 | 14 |

Energies agree within 0.5%. The **bulk is unchanged**: the median is identical and p90/p99 are only
7% and 14% above the `ff_2.1` control, so ff3.0 is not globally stiffer. The difference is entirely
in the extreme tail, 9 atoms out of 4949 that are hot under ff3.0 and not under `ff_2.1`, with the
peak going 111 -> 318. This tracks the tables themselves, whose repulsive maxima widened
(pair 27.8 -> 34.0, hydrophobe 21.2 -> 28.6).

Note the coverage nodes alone already double the peak force (65 -> 111) before any retraining, so
part of this is the coverage recipe rather than ff3.0.

This is at the *seed* configuration, which is unrelaxed, and a handful of stiff contacts normally
relax out in the first steps. It is not a failure and not a reason to change anything. It is the
specific thing to watch in the arm test's first block: if an arm destabilises, these ~9 atoms are
where to look first.

**What follows from this:**
* Cutting MAX_STEPS from 600 to 500 for the deadline cost nothing measurable. The objective was
  already flat, and 100 more steps would only have moved the parameters another ~14% at random.
* **Neither trained force field can be preferred on training grounds**, they are statistically
  equivalent fits. Only a downstream physical observable can choose, which is exactly what the
  arm test (TM1/TM4 helix fraction) does. Do not skip it and do not decide by inspecting the tables.
* **Read the arm test as the weaker statement it is.** Its two arms share 274 of 500 training steps,
  so they are more alike than two independent trainings would be. If the arms agree on TM1/TM4, that
  shows this pair agrees, NOT that any two retrainings would. A genuine reproducibility test needs a
  run branched at step 0.
* Any future claim that a retrain "improved" the force field must clear the 21% reproducibility
  floor before it means anything.


## 1. Standing rules and the measurements behind them

### 1.1 A spline table must BE the published potential

The MARTINI table shipped before 2026-08-05 was **not dry-MARTINI**. It stored bare LJ plus bare `1/r`
Coulomb, hard-truncated at 1.2 nm. Published dry Martini is run with `coulombtype = reaction-field`,
`epsilon_r = 15`, `epsilon_rf = 0` (infinite, i.e. conducting) and `vdw-modifier = Potential-shift`, so both
terms reach the cutoff smoothly. The stored table therefore had a step at `r_c` of `k*qq/r_c` = **2.65 E_up
for a charged pair, about 3.4 kT**, verified analytically (LJ -0.057 + Coulomb -2.648 = -2.705, matching the
stored -2.7049 exactly).

This is the correction flagged as outstanding in findings 75 and pinned down while hunting the glpG micelle
runaway (findings 79-80); `py/martini_build_tables.py` cites those numbers in its comments. Fixed there, in
BOTH nonbonded builders (the particle-particle grids and the SC-env/BB-env `_pair_energy_and_grads`, which
carried the same bare form including its analytic gradient).
`scratchpad/verify_table_matches_drymartini.py` asserts the equivalence by rebuilding the reference in
native kJ/mol + nm and converting:

| table | max deviation from the reference form | rows non-zero at the cutoff |
|---|---|---|
| old (`martini.h5.bak.pre-reactionfield`) | **3.95 E_up** | 81 / 81 |
| regenerated | 2.2e-11 E_up (round-off) | **0 / 81** |

Neutral pairs change by a constant only (forces bit-identical); charged pairs change by 3.1-3.9 E_up across
the sampled range, because the reaction field is a genuine r-dependent term. So this alters results for
anything with charges: ions, PC/PG headgroups, charged residues. Robertson et al.'s
`step8_production.mdp` independently confirms the same contract (`reaction-field`, `epsilon_r = 15`,
`epsilon_rf = 0`, `Potential-shift-verlet`, `rvdw = rcoulomb = 1.2`, `ref_t = 303.15 K` = our 0.8647 T_up).

Two things to remember. GROMACS `epsilon_rf = 0` means **infinity** (conducting boundary), so
`(eps_rf - eps_r)/(2 eps_rf + eps_r) = 1/2`, NOT -1; taking the 0 literally makes the charged triples look
~100% wrong at large r. And the current table is verified correct, equalling the analytic form to max
relative error **3.6e-11 across all 23 (eps, sigma, qq) triples** over r = 2.5-11.9 A, with the builder
asserting it at build time rather than assuming it. Do not re-litigate the table when hunting an
instability.

### 1.2 NO GUARDS

User directive, and now a rule in CLAUDE.md: no guard, anywhere, for any numerical problem. Removed from the
engine on 2026-08-05:

* `src/main.cpp` -- the blow-up guard that aborted the run on a non-finite potential or kinetic energy, plus
  its `blew_up` / `blew_up_ns` / `blew_up_round` state, the loop break and the FATAL block.
  `compute_logged_kinetic_energy` stays; it is still used for ordinary logging.
* `src/martini_potential.cpp` -- three silent-skip masks, which were the harmful ones. The pair loop's
  `if(isfinite(pot) && isfinite(force_mag))` dropped non-finite pair contributions outright; the main.cpp
  comment even admitted this by noting its kinetic check existed to catch "diverging momenta even when
  non-finite pair forces are masked out of the potential". The SC-env path had the same pattern twice.
* Kept in `update_martini_node_boxes`: `!(scale_xy > 0.f) || !(scale_z > 0.f)`. A box length being positive
  is a domain precondition, not a masked numerical error. The `isfinite` clause beside it was removed.

**Operational consequence, stated plainly:** a divergence now propagates instead of being stopped. NaN will
enter the log, be exchanged between replicas, and be written into restarts. That is the intended trade, the
evidence survives instead of the run, but it means a bad run must be caught by monitoring.

Three existing checks are in scope and were left alone pending a ruling, because none masks a numerical
error: `assert_environment_solvation` (prep-time, fails a build whose belt faces vacuum), the `np_hybrid.py`
ion assertions (neutrality and salt composition of a built system), and `run_np_prod.py`'s `health()`
peptide C-N check, which ends a chain rather than propagate a torn system. The last is closest to a guard
and is the one most worth a decision.

### 1.3 Timesteps are calibrated quantities, not free parameters

**NP-1AO6: dt = 0.001, and it must not be raised.** Two independent constraints, both measured:

1. *MARTINI LJ-core stability.* Velocity-Verlet is stable only for `dt < 2/omega`. Taking
   `omega = sqrt(U''/mu)` from the **tabulated** grid curvature, the limit collapses as pairs approach:
   `dt_max = 0.0209 @ 4.0 A, 0.00808 @ 3.5 A, 0.00430 @ 3.2 A, 0.00278 @ 3.0 A, 0.00186 @ 2.84 A`. The
   closest protein-environment approach the system samples is **2.84 A**, so dt = 0.005 was 2.7x over the
   limit, which is why failure was stochastic (16x spread in onset across six faces). Proven by an A/B
   restart from the identical pre-tear frame: dt = 0.005 gave 42 broken bonds and Epot -7988 -> +9130 with
   avg KE/1.5kT = 5.34; dt = 0.001 gave 0 broken bonds, Epot -7988 -> -8146, avg KE/1.5kT = 1.006.
2. *Backbone spring accuracy at large amplitude during unfolding.* After the findings-88 interface fix,
   dt = 0.005 still destroyed runs 1 and 4 at t ~ 250, at residues 119-122 which are 70+ A from the NP
   surface with no MARTINI pair anywhere near (nearest ion 8.6+ A). The spring linear stability criterion is
   satisfied throughout (`omega*dt = sqrt(48/0.5)*0.005 = 0.049 << 2`), but unfolding drives backbone bonds
   to 0.3-0.5 A amplitude, 2-3x the thermal `sqrt(kT/k) = 0.14 A`, and at that amplitude the bond force
   reverses sign over a half-period of ~0.32 t_u. With 54 force evaluations per frame interval the reversal
   is tracked too coarsely; a coincident alignment of CA-C and C-N spring forces put 14.4 E_up/A on C(121),
   injecting 12 E_up over 54 steps and stretching CA-C to 2.58 A.

   So "the springs are NOT the constraint, dt/dt_max = 0.035" holds only for small-amplitude thermal
   oscillation. **Do not raise NP dt above 0.001** even though the MARTINI contact limit allows much more,
   and note that a 50 t_u validation run is far too short to sample the unfolding regime (failure at
   t > 250).

**glpG hybrid: dt is hard-locked at 0.009** (`martini_brownian.cpp:100` throws on mismatch) because the
friction is tuned against it for a target lipid diffusion of 11.5 um^2/s. Changing one silently invalidates
the other. The stability problem there is handled by sub-stepping the integrator, not by changing dt (see
2.4).

Mass repartitioning would buy stability (protein sites are 1 m_up against 6 for a MARTINI bead, and
`dt_max ~ sqrt(mu)`) and is exact for an equilibrium observable such as an HDX free energy, since masses do
not enter configurational averages. But `/input/brownian` ties `numerical_time_step 0.009`,
`target_lipid_diffusion_um2_s 11.5` and `bare_particle_friction_up 0.169` together, so it voids the
friction/diffusion calibration. It is a physics decision, not a fix to apply unilaterally.

### 1.4 The dry-MARTINI unit contract

Native-to-Upside unit conversion happens ONCE, in the Python h5 build, and the engine does no unit math.

* The main `martini_potential` node's `coulomb_k` is READ but NEVER used in compute: the LJ+Coulomb energy
  is fully baked into `combined_energy_grids` by `convert_stage` (the particles group is already in eup:
  `unique_eps_eup`, `unique_sig_ang`, `combined_energy_grids`). The old conversion attrs on the config node
  were consumed only by the Python softening builder, not the engine, which is why moving the conversion was
  bit-identical on the bilayer (-11039.306641, |diff| = 0.000000). The single baked attr is
  `coulomb_constant`.
* `sc_table` conversion: keep the whole build/read logic in native units and convert + rename ONLY at the
  final `create_dataset` step (`grid_ang = grid_nm*10`, `*_energy_eup = *_kj_mol/2.914952774272`). The C++
  tail-subtraction shift is retained because it is a physics zeroing at the cutoff applied identically
  before and after; `(native[ig]-tail)/E == eup[ig]-eup_tail`, so the float32 result is unchanged (bake was
  exact, maxabs 0.0).
* Dimensionless datasets (`angular_profile`, `rotamer_angular_profile`, `cos_theta_grid`) are NOT converted.
* The engine reads `grid_ang` / `*_energy_eup` / `coulomb_constant` and does zero arithmetic on them. It is
  also analytic-LJ/Coulomb-free: `martini_potential` eval is purely spline (`combined_spline`), the
  node-level epsilon/sigma/coulomb_k reads and the analytic PairParam coefficients are gone, and
  lj_cutoff/coul_cutoff are unified into one cutoff attr.

### 1.5 dry-MARTINI membranes are NVT

A tensionless barostat is the wrong ensemble for these systems. Three reasons, and the second settles it
generally:

1. **A defect is a stress sink.** Lateral tension relaxes through a pore or an under-filled patch instead of
   through the lipid area, so zero measured tension stops marking the intact bilayer's equilibrium and the
   barostat compresses the intact regions past it.
2. **Implicit solvent has no solvent virial.** Dry MARTINI carries no water, so the lateral pressure the
   barostat reads is missing the solvent contribution entirely and is not the physical one. Dry Martini
   (Arnarez et al., JCTC 2015) is specified NVT for exactly this reason.
3. Papers that use semiisotropic Parrinello-Rahman on these lipids do so legitimately because **their
   systems are wet**; the water is merely stripped from the frames they deposit, which invites the wrong
   inference.

Measured consequence on our own POPE/POPG tile run both ways: under the tensionless xy barostat it condensed
to **APL 56.0 A^2 against the reference 61.7 (-9.4%)**, with tail core +1.6% and head-head +0.9%, i.e.
over-compressed and correspondingly thicker. The area is therefore an **input matched to an equilibrated
reference of the same lipids at the same temperature**, and validation moves to the area-sensitive
structural observables measured AT that area.

Before adopting an ensemble from a paper you are reproducing, check whether their system has the degrees of
freedom that make it valid.

### 1.6 Do not transfer settings, thresholds or analysis between simulations

NP and glpG are different simulations and must be analysed separately. Conflating them cost a healthy 6 h
glpG block (a threshold borrowed from NP false-positived) and produced a wrong root-cause writeup.

| | NP (`np_1AO6_prod`) | glpG (`remd_glpG-*`) |
|---|---|---|
| method | regular MD, 6 independent trajectories, single T=0.8647 | **REMD**, 48 replicas T=0.70-0.90 with configuration exchange |
| purpose | NP adsorption footprinting | **HDX** protection factors |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: ions + lipids + protein backbone overdamped Brownian; only the remaining protein atoms velocity-Verlet |
| timestep | free; fixed at 0.001 | **hard-locked 0.009**, friction tuned against it |

An energy test can never guard NP: in a forced tear the protein reached 431 broken bonds with the potential
still finite at +3e5. glpG's blow-up by contrast goes fully NaN within one 46-step interval, so there
`isfinite` suffices. The two jobs need different detectors.

### 1.7 No hard-coded protein or system identity

Scripts under `py/` are shared infrastructure. A hard-coded id does not fail loudly; it silently attaches
one protein's metadata to another system's trajectory. Instances found and fixed:

* `martini_extract_vtf.infer_pdb_id` returned `"1rkl"` for any non-bilayer input, which mislabelled the
  1AO6 nanoparticle VTFs; `martini_prepare_system` defaulted `--pdb-id` to `1rkl` and `--run-dir` to
  `outputs/martini_test_1rkl_hybrid` (fixed 2026-08-09).
* `martini_hdx_membrane_accessibility.py` defaulted `--tail-bead-names` to `C1,C2,C3`, the DDM tails. On
  POPE/POPG it found zero tail beads and raised, which was luck: a lipid containing a bead called `C1` would
  have been scored on the wrong subset silently. Tails are now detected from the trajectory by the MARTINI
  apolar naming pattern `^[CD]\d[A-Z]?$`, which covers a single-tailed detergent and a two-tailed lipid
  alike.
* The same class of error can arrive through a wildcard rather than a default: `4.calc_D_uptake.py` and
  `5.analyze_D_uptake.py` located the experimental HXMS arrays with `glob.glob(..._d_norm_peps_*.npy)` and
  took `matches[0]`, so the requested `protein_state` was ignored whenever more than one state matched. This
  dataset has four, and it was silently comparing against `pd9 SUB` while `pd9` was requested. **A glob that
  can match more than one identity must be resolved by the identity, not by `[0]`.** Both now call
  `helpers/function.py:select_state_file`, which raises naming the state and the available files.

---

## 2. The hybrid model: what each side supplies

### 2.1 Replaced by design

The `rotamer` node's 1-body input is deliberately swapped from Upside's implicit-solvent coupling to the
explicit MARTINI SC-env table:

| | rotamer arguments |
|---|---|
| standard Upside | `placement_fixed_point_vector_only`, `placement_scalar`, **`hbond_coverage`**, **`hbond_coverage_hydrophobe`** |
| hybrid | `placement_fixed_point_vector_only`, `placement_fixed_scalar`, **`martini_sc_table_1body`** |

Upside's implicit bilayer (`membrane.h5` via `--membrane-potential` / `write_membrane_potential{,3,4}` /
`write_membrane_lateral_potential`, plus `--membrane-thickness`) is correctly absent: an explicit-lipid
model should not also carry an implicit slab.

### 2.2 Absent and NOT replaced (findings 101)

MARTINI supplies only protein-environment interactions, so three protein-protein terms of stock Upside are
simply gone:

* `sigmoid_coupling_environment` <- `environment_coverage_sc`: the many-body **protein self-burial** term.
  MARTINI's 1-body rewards lipid contact, not helix-helix packing.
* `bb_sigmoid_coupling_environment` <- `environment_coverage_hb` + `cat_pos_bb_coverage`, and
  `hb_environment_coverage_hn/oc`: backbone burial coupling.
* `hbond_coverage` / `hbond_coverage_hydrophobe`: sidechain-to-backbone-H-bond competition, which in
  standard Upside is solved **inside the rotamer solver**. MARTINI has no H-bond concept at all.

These are core force field, not niche: **all 24** master example configurations pass `environment.h5` +
`bb_env.dat`, including **all six** `08.MembraneSimulation` scripts, which use them *alongside*
`membrane.h5`. (An earlier grep missed the membrane example because it searched only the
`--environment-potential` CLI form and not the `environment_potential=` kwargs form.)

Restoring them is not a flag flip: it changes the `rotamer` node's arity and requires deciding how Upside's
protein-burial 1-body composes with MARTINI's lipid 1-body inside the rotamer solver, a C++ interface
question. Scale check on real structures: the self-burial term disfavours the drifted state by only
**8.0 E_up = 5.6 kcal/mol**. The experiment that restored them is findings 103 in section 4.1, and it did
not recover the fold.

### 2.3 The backbone interface: sites, force routing, CB placement

**Only BB is on the protein side of the pair list.** Dry-MARTINI represents a residue backbone as one BB
bead, and N/CA/C/O are the atoms that bead stands for. Earlier builds put all five sites in the pair list at
full epsilon, counting the backbone-environment interaction five times over; in energy the over-count was
only 1.15x, because O's placement lands it inside env repulsive cores and the four atom sites largely cancel
(+2117.5 against -2295.7 E_up at one measured frame). Current builds and cluster seeds carry one BB site per
residue. A corollary worth remembering: the backbone `O` being the closest protein atom to the environment
(2.97-3.32 A, against 4.0-4.4 A for BB) is *expected*, not a defect, because O carries no MARTINI
interaction and nothing repels it.

**BB is a derived site built from nodes Upside already differentiates.** `HybridPositionNode` takes `pos`
and `infer_H_O`, and BB is the mass-weighted centre of N/CA/C/O with weights [14,12,12,16]/54, where the O is
*Upside's own derived carbonyl O* rather than one rebuilt from a stored local frame. Every term is then a
linear combination of node outputs, so `propagate_deriv` is a constant-weight split: N/CA/C shares go to
`pos.sens`, the O share goes into `infer_H_O`'s sensitivity, and its chain rule carries it back to CA, C and
the next residue's N. No hand-written placement Jacobian is needed, and the 184 lines of frame/Jacobian
helpers were deleted rather than kept as a fallback.

Newton's third law holds by construction: `martini_potential` runs on the hybrid node's output, a
pass-through copy of `pos` with only the BB and O slots overwritten, so an environment particle's
sensitivity returns to `pos.sens` unchanged while the protein side is redistributed by weights summing to 1.
Measured on 1rkl: force remaining on the BB and O slots is exactly **0**, and the total force sums to
**1.69e-09** of the largest single force.

The C-terminal residue has no acceptor: `infer_H_O` builds a carbonyl O from the *next* residue's N, so it
emits `n_res - 1` acceptors. That one BB is the N/CA/C mass centre with weights renormalised, and its O slot
is left as it came in. This is not a guard, there is genuinely no such site, and no MARTINI term reads that
slot.

**CB placement omitted the frame-origin subtraction (findings 102).** `affine_alignment` builds each residue
frame with its **origin at the centroid of N/CA/C** (`src/eig.cpp`: `center = (atom1+atom2+atom3)/3`), so a
point expressed in that frame must be given relative to the centroid. `upside_config.write_environment` does
that; the hybrid prep stored the raw CB coordinate:

| | placement_data (frame coords) |
|---|---|
| standard Upside | `[-0.019807,  1.511741, 1.206801]`  = CB - centroid |
| hybrid | `[ 0.000000,  0.943756, 1.206801]`  = CB |
| difference | `[ 0.019807, -0.567985, 0.000000]`  = exactly the centroid |

The frame is orthonormal, so the displacement is **0.568 A in Cartesian space for every residue**, moving
the site that anchors the entire sidechain-environment term (`martini_sc_table_1body` takes
`placement_fixed_point_vector_only_CB` as its sidechain input). At a MARTINI bead sigma of 4.7 A that is
~12% of a contact radius. Fixed 2026-08-15: `CB_PLACEMENT` derives from an explicit reference geometry with
the frame origin subtracted, matching `upside_config.write_environment` to 5.8e-8. `CB_VECTOR` is unchanged
because CB-CA is a difference, and `martini_build_tables.py` needs no change because it uses CB only as a
relative origin. The total potential of an identical configuration moves by **~+680 E_up**, confirming the
displacement was materially biasing SC-env energetics.

Note how it was found: by comparing the hybrid against a *standard* Upside config of the same protein,
array by array. The bonded and hydrogen-bonding core came back bit-identical (`rama_map_pot`,
`hbond_energy`, `protein_hbond`, `backbone_pairs` all max|delta| = 0), which is what made the two arrays
that did differ worth reading rather than dismissing as index remapping. Configs predating the fix: the four
cluster REMD chains and their seeds, and the local production seed behind the delivered HDX figure; the NP
campaign was rebuilt (section 9).

### 2.4 The integrator the hybrid actually runs (findings 123)

`DerivEngine::integration_cycle` (`src/deriv_engine.cpp:396`) begins with

```
if(martini_brownian::has_brownian(this)) {
    compute(DerivMode);
    martini_brownian::apply_langevin_step(this, mom, dt);
    return;
}
```

so as soon as `/input/brownian` exists, which it does for every hybrid config, the function returns before
reaching the **three-stage Predescu et al. (2012) integrator** the same function uses otherwise
(`mom_update = {1.5-3a, 1.5-3a, 6a}`, `pos_update = {3b, 3-6b, 3b}`, one force evaluation per stage). Stock
Upside examples run that three-stage scheme at dt = 0.009; the hybrid ran **one** g-JF stage at dt = 0.009:

`x <- x + (b dt/m) p + (b dt^2/2m) f + (b dt/2m) beta`,  `b = 1/(1 + alpha dt / 2m)`

One stage means a force spike is committed to displacement with no intermediate force re-evaluation. Same
dt, different stability.

`/input/brownian` covers **4529 of 4949 atoms**: all 272 ions, all 3627 lipids, and **630 of 1050 protein
atoms** (the N/CA/C backbone), with friction 0, 0.1692, 0.3384 and 0.5075, interface-dependent and higher
near lipid within a 12 A cutoff. So the protein backbone is inside the single-stage Langevin path and
carries the interface friction, and cannot simply be removed from the list. Measured single-step
displacements `b dt^2 |F| / 2m`: **protein max 0.010 A** (99.9th 0.006, median 0.0005) against **lipid max
0.0003 A**, i.e. ~30x, from one sixth the mass and stiffer forces.

`--integrator mv` does not help despite its name: `build_integrator_levels` makes `integrator_level == 1`
the **slow** set integrated at `dt * inner_step` and level 0 the fast set at `dt`. It is a cost optimisation
for expensive *smooth* terms, and every potential node in the hybrid config takes the default level, so
`mv`'s slow set is empty.

**Fix implemented (findings 124): RESPA-style g-JF inner sub-stepping.** `n_inner_steps = N` wraps the
position update, inner force evaluation and momentum update in a loop of N inner steps at `dt_i = dt/N`. The
outer `dt` seen by the engine, the `numerical_time_step` check and the friction/diffusion calibration all
remain at 0.009; only the g-JF integrator is sub-stepped. In `src/martini_brownian.cpp`:
`BrownianRuntime::n_inner_steps` (default 1, backward compatible), read from the `/input/brownian` attribute
`inner_steps` and throwing if < 1; `apply_langevin_step` loops with a full position update, an
`engine->compute(DerivMode)` and a full momentum update per iteration, indexing the invocation counter by
`outer_invoc x N + inner` to keep random streams distinct. No other file changed.

| r [A] | F [E_up/A] | kick N_inner=1 [A] | kick N_inner=9 [A] |
|---|---|---|---|
| 2.853 | 1.27e5 | 5.1 | 0.064 |
| 2.584 | 4.64e5 | 18.8 | 0.232 |
| 2.432 | 1.02e6 | 41.3 | 0.510 |
| 1.783 | 5.60e7 | 2268 | 28.0 |

At every approach distance that occurs thermally (>= 2.43 A), N_inner=9 keeps the kick below 1 A; the
1.783 A row is for completeness, since the LJ potential there is ~10^7 kT. Thermostat and diffusion are
preserved, because dissipation per outer step is `prod_{i=1..N} (1 - alpha dt_i/(2m)) ~= 1 - alpha dt/(2m)`,
identical to N=1, leaving `D = kT/alpha` unchanged. Local timing (1 replica, 4000 steps): 12.9 ms/step at
N=1 against 55.4 at N=9, a **4.3x overhead** rather than 9x, because `engine->compute()` is only ~41% of
step time. Thermodynamics validated: `avg_kinetic_energy/1.5kT` 1.011 baseline and 1.002 fixed, potentials
~-22 000 E_up, Rg ~20.5 A in both arms. No blow-up was captured locally, because the 79HIS configs are in a
conformational state that does not sample 2.43 A protein-MARTINI contacts in that window; the 79ALA variants
have the susceptible conformation.

An engineered local A/B did not produce a clean contrast, because the test script used the wrong BB formula
(`martini_hybrid_position` uses `infer_H_O`-derived O positions, and only N/CA/C renormalised for the
C-terminal residue) and because moving one atom by >3 A in a dense bilayer overlaps several neighbours at
once. A clean local A/B needs a toy two-body system or a pre-failure frame from the cluster.

Cluster deployment: rebuild the binary, then patch each running config with
`f['/input/brownian'].attrs['inner_steps'] = np.int32(9)`; the binary reads `inner_steps = 1` when the
attribute is absent, so it is backward compatible with every existing config.

### 2.5 The friction clock

The live calibration is the sub-molecular fallback: `D_raw = 4*11.5 = 46 um^2/s`, `dt_raw = 40/4 = 10 ps`,
`D_bead,up = D_raw*1e-4*dt_raw/.009 = 5.1111 A^2/U`, and `alpha_bead = kT/D_bead,up = 0.1691804`. Each
environment bead receives this friction; a real protein N/CA/C carrier receives `n_contact*alpha_bead`,
where `n_contact` is the number of lipid beads inside the existing 12 A spline cutoff. Counts are refreshed
after stage handoff, minimization promotion and production continuation, then held fixed during each segment
so the SDE does not silently acquire position-dependent multiplicative noise.

Measured molecular DOPC diffusion under it is only **0.013-0.015 um^2/s against the 11.5 target**. This is
reported as a failed molecular target, not hidden. Name the calibrated observable in H5 and in the paper,
and never report a friction-calibrated trajectory as having the target lateral diffusion.

Two related facts. `1 T_up = 350.588235 K`, so `--temperature 0.8647` is 303.15 K; the MARTINI factor four
changes time, not thermodynamic temperature. And one temperature controls the whole system: the workflow
assigns the DOPC friction reference directly from the single authoritative `TEMPERATURE` and overwrites any
independently supplied value, because calibrating friction at one kT while driving its noise at another
changes the nominal diffusion.

**Never accept a kinetic calibration from temperature and structural stability alone.** An earlier mapping
(`tau_up = .0036`, giving `alpha*dt/(2m) = 1.25` for a mass-6 bead) produced trajectories that were
effectively frozen (0.0081 A of drift-removed COM motion per saved frame) while still showing a
thermal-looking momentum distribution and retained secondary structure. Gate every friction change on
protein displacement and whole-molecule lipid MSD in the saved trajectory, and inspect the VTF rather than
only the H5 statistics.

---

## 3. Known defects: root causes and fixes

### 3.1 The LJ core table was force-free, and particles reached it (findings 90, 92, 93)

`py/martini_build_tables.py` evaluated the grid at `r = max(r, 0.1*sig)` on a domain starting at r = 0, so
below 0.1*sig (0.47-0.60 A) the tabulated potential was a **constant** and therefore exerted **no force at
all**. Per-step instrumentation (`UPSIDE_MARTINI_PAIR_DIAG`, two 48-replica jobs with exchange disabled so a
blow-up stays in the slot that made it, 340 000 steps, 9.2 h) measured what that allows:

| | job A | job B |
|---|---|---|
| reported approaches < 1 A | 4742 | 5032 |
| approaches < 0.6 A (inside the floored plateau) | 1446 | 1432 |
| closest approach | **0.0355 A** | **0.0466 A** |
| largest force delivered anywhere | 3.4e11 E_up/A | 1.0e11 E_up/A |

At 0.0355 A the true dry-MARTINI LJ force is **5.8e29 E_up/A**, so the table's largest delivered force
anywhere in the run was about **18 orders of magnitude too weak**. The clearest single event: a pair at
0.0804 A while the whole box's maximum force was 7.28 E_up/A, i.e. that pair felt nothing. Offending pairs
are **environment-environment** (LIPID-LIPID and LIPID-ION), which is why `lipid_kinetic` is what explodes
while the protein KE is merely NaN.

An earlier reading (findings 90) argued the catastrophic region was ~500 kT out of reach, computing one-step
displacement for an *inertial* mass-72 bead (0.0058 A at r = 3 A). That was wrong: ION and LIPID are
integrated as **overdamped Brownian**, whose step is proportional to the force, so a large force gives a
large displacement that can overshoot *through* a partner, and once inside the force-free plateau nothing
ejects it. Entry by overshoot is inferred; the missing exit is measured. **Check which integrator governs
the particles before computing a stability margin for them.**

Findings 90 also ruled out, each by measurement: a stale pair list (`cache_buffer` = 2.0 A with
`pairlist_needs_rebuild` before every force evaluation, so the list is at most one step old); minimum image
(`simulation_box::minimum_image` uses `roundf(dr/box)`, correct for arbitrarily large separations, which
matters because unwrapped ion positions reach 489 A in a 137.4 A box where a single-shift implementation
would be wrong); and the table build formula (the grid reproduces its own analytic expression to 2.3e-12).
The event itself is a **single-frame catastrophe from a fully healthy state**: potential -7070 with KE
1.51/1.37 in one frame, non-finite with `lipid_kinetic` = 6.246e18 in the next ~60 steps later, then
random-walking the ladder by exchange and destroying every slot it lands in, which is why 48/48 replicas end
up destroyed from a single event.

Fix, two coupled changes because the domain was declared in one place and assumed in another:
* `py/martini_build_tables.py`: floor removed, `PARTICLES_R_MIN_A` 0.0 -> 0.3, grid built vectorised over
  the true potential, plus an assertion that every grid point equals the analytic form to 1e-12 relative.
  0.3 A is far inside anything reachable; the core there is ~5e17 E_up/A.
* `src/martini_potential.cpp`: `r_min`/`r_max` were **hardcoded to [0, 12]** and ignored the
  `r_min_ang`/`r_max_ang` attributes the builder already wrote, so changing the builder's domain alone would
  have mapped every distance onto the wrong knot. It now reads and validates the domain from the table.
  Old `.up` files still read correctly, since their own attrs say `r_min = 0`. Relatedly,
  `inject_particles_table` restated the domain instead of copying it, and the C++ hardcoded the 1000-point
  grid size in four places; the grid geometry is now declared once by the builder and carried through.

Verified on a real glpG-DDM system with the corrected table: initial potential -7811.76, min pair distance
4.0349 A, max force 33.3441 E_up/A, **identical to the old table**, so the change is a no-op everywhere the
system actually samples, while the core now delivers 8.5e15 E_up/A at 0.35 A where the old table delivered
~0.

**The fix removed the consequence, not the entry (findings 93).** On the corrected table the same system
still reached **0.2078 A**, with 32 approaches under 1.0 A and 6 under 0.3 A, and forces up to
1.27e15 E_up/A: under the old table those pairs coasted through force-free, now they receive an enormous
impulse and the run dies. The residual dead zone below the new `r_min = 0.3 A` was entered, which is exactly
the risk recorded, because the clamped spline still returns a constant with zero derivative below its
domain.

**And the core is where the cascade ends, not where it starts.** One capture with diagnostics running gives
the causal order directly: two environment beads reach 1.83 A with max force 3.4e6 (inside the valid table
domain, where the tabulated force is correct and simply huge), then 1.59 A on a BB-proxy/lipid pair, then
1.05 A at 5.5e10, and only 1072 diagnostic steps after that first enormous force does anything reach
0.245 A. Removing the floor was right and was never going to prevent this.

One more defect the first corrected-table run exposed (findings 93): **the MD loop destroyed its own error
messages.** `src/main.cpp:1252` integrates systems under `#pragma omp parallel for` with no exception
trapping, and an exception cannot leave an OpenMP structured block, so any `throw string(...)` from the
engine mid-run called `std::terminate`: a local POPE/POPG run died at step 155 640 reporting only
`libc++abi: terminating due to uncaught exception`. The setup loop already traps for exactly this reason;
the integration loop did not. It now traps per system, keeps the first message, and rethrows once serial.
A diagnostic that cannot report is worse than none.

### 3.2 The blow-up mechanism: a 1 m_up protein site ejected by the MARTINI wall (findings 122)

Localised by pulling a 10-frame window around the onset of a `glpG-RKRK-79ALA` event off the cluster and
re-running it against the local engine.

At frame 105 the total potential is +2.536e5 while Rg is 19.6 A and |pos|max 130 A, so global observables
see nothing. The excess is entirely in `Spring_bond` (280 -> 274 347 E_up), and per-bond it is three bonds
of one residue: **CA of residue 170 sits 80.7 A from its own C (r0 1.526) and 69.8 A from its own N
(r0 1.453)**, worth 150 522 + 112 250 + 8 502 E_up. One atom has been ejected; the rest of the protein is
intact. Re-evaluating the recorded coordinates locally reproduces the recorded potential to 0.03 E_up.

Every protein site carries **mass 1 m_up** while every lipid and ion bead carries **6**, so the same force
throws a backbone site six times as far. The one-step kick `F dt^2 / m` at dt = 0.009 (note the g-JF update
carries a factor b/2, so the true prediction is about half these numbers; the order of magnitude is what
matters):

| separation | steepest pair force | dt for a 1 A kick |
|---|---|---|
| 2.853 A | 1.27e5 E_up/A | 0.0028 |
| 2.584 A | 4.64e5 | 0.0015 |
| **2.432 A** | **1.02e6** | **0.00099** |
| 1.783 A | 5.60e7 | 0.00013 |

The run lives there: over the ten frames the closest interaction-list pair is 2.49-2.74 A in every clean
frame (1.78 A in the bad one), with ~169 pairs per frame inside 3.40 A, ~15 inside 2.85 A and **0.2 per
frame inside 2.43 A**. An ejection is not an accident, it is the expected outcome of continuous sampling at
that separation whenever the bead that happens to be there is a protein site rather than a lipid one. The
steepest tabulated force is 5.2e17 E_up/A at the 0.3 A inner edge, which is where a replica's 8e18 potential
and |pos| 6.9e10 come from: once a pair is driven to the inner edge the kick is unbounded.

Consistency checks that make this an explanation rather than a story: the ejected atom is a protein site,
the lightest species; `martini_hybrid_position` rises with it (1700 -> 12 988) as the ejected site drags its
proxy; and the propagation matches the exchange arithmetic exactly (37 steps/frame with
`--replica-interval 0.09` = 10 steps gives 3.7 exchanges per frame, and the wreck moves ~4 rungs per frame).
The bond period is ~100 steps and the thermostat timescale 555 steps, so an 80 A excursion cannot relax in
37 steps: the clean frame that follows is a *different configuration swapped in*, not recovery.

**Global observables are the wrong instrument for a local failure.** Rg, |pos|max and the peptide C-N scan
all passed on a frame carrying +2.5e5 E_up, because one atom in 4949 was 80 A out of place. Term
decomposition found it in one step where three rounds of Rg-and-C-N checking had not.

### 3.3 REMD launders a wreck around the ladder, and the driver rolls it back

All four POPE/POPG variants carry non-finite potentials in about 0.15% of frames, in nearly every replica
file, arriving in pairs two frames apart on the exchange period with clean frames on either side.

**I first concluded this was an output defect, and that was wrong.** `run_remd.py` says what actually
happens, in its own docstring and in `destroyed()`: any non-finite potential in a chunk is treated as a
blow-up, that replica is rolled back to its pre-chunk positions, and the NaN chunk is rotated to
`output_previous_N` as normal historical data. So a replica genuinely blows up; exchange carries the wrecked
configuration around the ladder, which is why NaN appears in nearly every replica file and why the
neighbouring frames look clean (those slots held *other*, healthy configurations at the time); and at chunk
end the driver rolls the affected replicas back and keeps the chunk as history. `grep ROLLBACK` on the logs
confirms it and shows the events cluster by chunk rather than by replica, one chunk rolling back all 48
replicas and one replica needing five consecutive rollbacks.

The analysis is nevertheless sound, and not because the NaN are rare: `martini_remd_concat.py` keeps a frame
only if its potential is **finite and negative**, a physical test rather than a finiteness test. A condensed
bilayer plus protein sits near -2.2e4 E_up, so a positive total means overlapping cores, and its docstring
records a replica that stayed finite for 96 frames at +1.9e6 before reaching NaN, exactly the ramp a
NaN-only filter would have kept.

Two lessons: **a clean neighbourhood is not evidence of a clean trajectory** (read the producer before
explaining its output), and **check the log the tool already writes before inferring a mechanism from the
data**.

## 3.10 The glpG blow-ups: the two subsystems are not at the same temperature (2026-09-10)

Diagnosed after 15 rollbacks appeared in block 1 of the four production chains. The user's hypothesis
(a temperature mismatch between dry-MARTINI and Upside) is confirmed and quantified below. The dt
lock prevented one planned test: `apply_langevin_step` throws when the runtime dt differs from
`/input/brownian numerical_time_step`, so a dt scan is impossible without rebuilding the node, and
that check was left alone.

**The Upside temperature conversion is exactly K = T_up x 350.588235.** The reference table in
`~/OneDrive - The University of Chicago/image.png` was checked row by row and is internally
consistent on all 24 rows (max |dK| < 1e-5 K, |dC| <= 0.05 from its own rounding). Using it:

| | T_up | K | C |
|---|---|---|---|
| ladder rung 0 (coldest) | 0.700 | 245.4 | -27.7 |
| ladder rung 27 (hottest) | 0.820 | 287.5 | +14.3 |
| dry-MARTINI reference | 0.8647 | 303.15 | +30.00 |

**Every dry-MARTINI parameter is built for 0.8647, and the ladder runs entirely below it.** Three
independent places carry that same number, so it is the model's design temperature and not a stray
constant: `/input/brownian` `reference_temperature_up` = 0.8647, where the friction is fixed for a
target lipid diffusion of 11.5 um^2/s (and `bare_particle_diffusion_up` = 5.1111 = 0.8647/0.16918
confirms D = kT/gamma is evaluated there); `py/martini_build_tables.py`
`DEFAULT_PRODUCTION_TEMP_UPSIDE` = 0.8647; and the equilibration itself, since `output_previous_0`
of every replica of every variant is a single-temperature run at exactly 0.8647, after which
production dropped to the ladder. The bilayer is therefore run **15.7 to 57.7 K below** the
temperature its friction, its tables and its equilibration all assume. Because MARTINI energies are
fixed in absolute units, at rung 0 every MARTINI interaction is 0.8647/0.70 = **1.24x stronger in
units of kT** than at the reference, which over-condenses the bilayer rather than merely slowing it.
The two scales are not interchangeable: 0.7-0.9 is a sensible reduced-temperature folding range for
Upside's trained statistical potential, where it carries no Kelvin meaning, but the same number is
handed to the MARTINI subsystem as a literal kT in the Brownian noise.

**Measured, the protein and the lipids sit at different temperatures.** On clean frames only (finite
and negative potential), 79ALA, all 28 rungs: **lipids track their set point** at T_lip/T_nom = 1.020
flat across the ladder, while **the protein does not**, sitting at T_nom + ~0.08 T_up (about +28 K)
and reaching 1.506 at rung 27 against a nominal 0.820. The rungs that rolled back (27, 18, 17, 14,
26) are the ones with the largest excess. Two controls make the offset real rather than an artefact
of the wreck: it is present at the same ~1.10 ratio in cold rungs 0-13 that never rolled back, and
the lipids occupy the *same slots* under the *same* exchange yet stay on target.

`protein_kinetic` is a trustworthy instrument here, which was worth checking: the 420 O/sidechain
slots placed by the `placement_*` nodes carry **exactly zero** momentum, but they are also excluded
from the logger's `n_dynamic` count, so the logged value is the mean over the 630 dynamic backbone
sites and is not diluted. (A naive per-atom average over all 1050 `PROTEIN` slots gives 0.661 x T_nom
instead of 1.136 and is simply wrong.)

**Reproduced locally on one system with no replica exchange**, started from an equilibrated cluster
frame. At production settings (T = 0.7215, tau = 5) the local run gives T_prot = 0.7911 against the
cluster's 0.7912 for that same rung, and T_lip 0.7327 against 0.7355. Exchange is therefore not
involved at all; the excess is intrinsic to the hybrid integration. Two scans characterise it:

| | T_nom | T_prot | T_lip | T_prot/T_nom |
|---|---|---|---|---|
| tau = 1 | 0.7215 | 0.8487 | 0.7349 | 1.176 |
| tau = 5 | 0.7215 | 0.8197 | 0.7489 | 1.136 |
| tau = 20 | 0.7215 | 0.8600 | 0.7412 | 1.192 |
| tau = 5 | 0.8647 | 0.9821 | 0.8796 | 1.136 |

The excess is **multiplicative and independent of both the thermostat timescale (over a 20x range)
and the temperature** (identical 1.136 ratio at 0.7215 and 0.8647). It is therefore not a power leak
the thermostat fails to remove, which would scale as P*tau. A first reading of partial data suggested
it grew superlinearly with temperature; that was a transient burst contaminating a half-length
average and is wrong.

**It is the mass-1 backbone, under either thermostat.** With momentum logging, split by thermostat
mechanism:

| group | T/T_nom at 0.7215 | T/T_nom at 0.8647 |
|---|---|---|
| backbone, friction > 0 (g-JF Brownian) | 1.092 | 1.077 |
| backbone, friction == 0 (OU thermostat) | 1.130 | 1.072 |
| lipid/ion (g-JF Brownian) | 1.024 | 0.999 |

Both protein subsets are hot by the same amount and there is no trend with lipid-contact count
(1.04-1.18, scattered), so the interface friction and the OU thermostat are both exonerated: the bias
belongs to the mass-1 protein sites themselves. That is integrator discretisation bias, which is
multiplicative and tau-independent exactly as measured. It is not the explicit springs, whose
stiffest mode (`Spring_angle`, k = 175) gives only (omega*dt)^2/4 = 0.7% at dt = 0.009; the steep
MARTINI pair core is the only curvature in the system large enough, and it reaches the protein
through `martini_hybrid_position`, since `martini_potential` takes the proxy positions as its
argument.

**And the cold bilayer is what presses the backbone into that core.** Measured on the two local runs,
minimum-image protein-backbone-to-environment distances (raw coordinates, so a comparative proxy for
the true proxy-mediated distance rather than the exact interacting one):

| | mean of per-frame minimum | closest seen | pairs/frame < 3.40 A |
|---|---|---|---|
| T = 0.7215 (production) | 3.512 A | **2.887 A** | **0.38** |
| T = 0.8647 (design point) | 3.623 A | 3.289 A | 0.07 |

**5.4x more sub-3.40 A contacts at the production temperature**, with the closest approach falling
from 3.29 to 2.89 A. Per 3.2 the pair force at 2.853 A is 1.27e5 E_up/A, where dt for a 1 A one-step
kick on a mass-1 site is 0.0028, so dt = 0.009 is already 3.2x too large there. The causal chain is
therefore: the bilayer is run 16-58 K below its design point, which strengthens every MARTINI
interaction by up to 1.24x in kT and over-condenses the environment; that drives protein backbone
sites measurably closer into the steep core; and mass-1 sites at dt = 0.009 integrate that core
inaccurately, giving a standing 8-13% kinetic excess with intermittent bursts to 1.5-1.8, until one
site is ejected and tears the TM4 backbone.

**What the local runs did NOT do: produce a blow-up.** Over 250 time units they stayed finite and
negative, with zero pairs inside 2.85 A. They reproduce the standing temperature split and the
contact-density shift, which are the precursors; the ejection itself is a rare event (3.2 measured
0.2 pairs/frame inside 2.43 A only in long cluster runs). So the last link, that the increased
contact density is what produces the ejections, is consistent with everything measured but is
inferred rather than demonstrated here.

**The blow-up itself is a backbone tear in TM4, and not a force-field-table defect.** Decomposing a
finite-but-positive onset frame (79ALA r27, `output_previous_2`, frame 208, -19765 -> +1683 ->
+43240 E_up) by node puts essentially all of the excess in `Spring_bond` (1272 -> 21724 -> 62558),
with `Spring_angle` +655 and `Spring_omega` +332 and every MARTINI term flat. Per bond it is
consecutive backbone bonds of residues 139-141: at frame 208 `C140-N141` is at 18.5 A against
r0 = 1.300, `CA139-C139` at 16.7 A and `N139-CA139` at 16.3 A. That is TM4. Re-evaluating the
recorded coordinates on the local engine reproduces the recorded total to **0.27%** (+43356 vs
+43240) while using the *older* local tables rather than the deployed arm-R ones, so the catastrophe
is geometric and the arm-R retraining is not its cause. Note also that the frames immediately before
onset are already far out of equilibrium: peptide bonds sit at 2.0-2.9 A against r0 = 1.3, roughly
60 kT of bond strain, where equipartition at T = 0.82 with k = 48 allows 0.13 A rms.

Ruled out by measurement, each: the arm-R tables (above); an unthermostatted subset, since
`stochastic_mask` is set only where `friction > 0` (`martini_brownian.cpp:78`) and the OU thermostat
therefore still reaches every friction-zero atom (`thermostat.cpp:31`), so each atom is thermostatted
by exactly one mechanism and both target the same kT; and exchange laundering of the protein excess.

**A measurement trap worth keeping: lipid diffusion cannot be measured from production output.**
Replica exchange puts a different configuration in a slot every exchange interval, so consecutive
frames of a production chunk are not a trajectory. Measured on PO4 beads, the apparent lateral D
*falls* with lag in production (9.07, 5.14, 2.56, 1.43, 0.82 A^2/time_up at lags 1-16), the signature
of frame-to-frame discontinuity, while the exchange-free 0.8647 equilibration behaves like a real
trajectory and *rises* with lag (0.025 -> 0.103). Use a continuous single-temperature run for any
transport observable.

### 3.10a The fix, what it verifiably does, and what it does not (2026-09-10)

The ladder was made authoritative and dry-MARTINI brought to it (plan.md, Revised decision
2026-09-10). Two changes, both verified; one deliberate non-change; and one claim that could **not**
be tested.

**Change 1: friction follows the replica temperature** (`src/martini_brownian.cpp`). Friction is
built as `gamma = kT_ref/D_target` with T_ref = 0.8647, so at rung 0 the realised lipid diffusion was
`D_target * 0.70/0.8647`, 19% below the 11.5 um^2/s the node exists to deliver. The runtime now
scales gamma by `T/T_ref`, giving `D = kT/gamma(T) = D_target` at every rung and through exchange.
Keyed on `reference_temperature_up`, so a config that never declared one is untouched. This changes
**only transport, not thermodynamics** -- friction does not enter the Boltzmann distribution, so
potential statistics and exchange acceptance are unaffected. Verified: at T = 0.8647 the new binary
is **bitwise identical** to the old over 200 steps (scale is exactly 1 there), and at T = 0.7215 it
differs (max |dpot| = 40.5 E_up), so the scaling engages where it should and nowhere else.

**Change 2: `inner_steps` = 4 by default** (`py/martini_prepare_system_lib.py`, env
`UPSIDE_MARTINI_INNER_STEPS`). N substeps of `dt/N` inside each outer step, noise using `dt_i` so FDT
holds. The **outer dt stays 0.009**, so `numerical_time_step`, the friction/dt lock and the 40
ps-per-step clock are all untouched; it changes no force field and no parameter, it integrates the
same equations more accurately. The capability already existed in C++ and no Python code had ever
written the attribute. Measured on glpG-RKRK-79ALA from an equilibrated frame at T = 0.7215:

| inner_steps | backbone (Brownian) | backbone (OU) | lipid | logged T_prot | excess | wall |
|---|---|---|---|---|---|---|
| 1 | 1.092 | 1.130 | 1.024 | 1.101 | +10.1% | 262 s |
| 2 | 1.034 | 1.039 | 1.018 | 1.035 | +3.5% | 333 s (1.27x) |
| 4 | 1.010 | 1.011 | 1.009 | 1.011 | **+1.1%** | 538 s (2.05x) |

**The temperature mismatch is resolved**: at N = 4 both protein thermostat groups and the lipids sit
within ~1% of the set point, so the two subsystems are finally at the same temperature. Cost is far
below the naive Nx because the baseline already does two force evaluations per step and a substep
adds one, so the ratio is `(N+1)/2`.

**Non-change: MARTINI epsilons are not rescaled.** Solute tempering (`eps * T/T_ref`, preserving
`eps/kT`) is the textbook fix for the remaining defect and is forbidden here, because a spline table
must equal the published dry-MARTINI form exactly. So the over-condensation stands: at rung 0 every
MARTINI interaction is still 1.24x too strong in kT, and the 5.4x excess of sub-3.40 A
protein-environment contacts is unchanged.

**Not demonstrated: that any of this stops the blow-ups.** Two separate reasons, both quantitative.

*The unsafe window only shrinks as 1/sqrt(N).* Differentiating the deployed
`combined_energy_grids`, the separation below which one step throws a mass-1 site more than 1 A is:

| inner_steps | dt_i | r(kick > 1.0 A) | r(kick > 0.3 A) |
|---|---|---|---|
| 1 | 0.00900 | 3.228 A | 3.544 A |
| 2 | 0.00450 | 2.900 | 3.181 |
| 4 | 0.00225 | 2.607 | 2.865 |
| 8 | 0.00112 | 2.350 | 2.572 |
| 64 | 0.00014 | 1.705 | 1.869 |

At N = 1 the run sits continuously inside its own unsafe window (0.38 pairs/frame below 3.40 A),
which is the mechanism. N = 4 shrinks it to 2.607 A but 3.2 measured ~0.2 pairs/frame inside 2.43 A
on long runs, still inside. **Covering the measured close-approach population needs N = 8**; the
1.78 A approaches of 3.2 would need N ~ 64.

*And the direct test was underpowered by 50x.* Starting from the recorded last clean frame of the
real 79ALA r27 event (frame 207, potential -19765, one backbone bond already at 6.27 A, one frame
before the +43240 catastrophe), four seeds at T = 0.82 with N = 1 and four with N = 8: **all eight
held**, relaxing to -21100..-21550 with no non-finite or positive frame. The N = 1 arm did not
reproduce the tear, so the comparison carries no information. The rate explains why: the cluster
shows 15 rollbacks over 11.4M replica-steps, one per ~762k, and this test sampled 13.3k steps, 1.8%
of one expected waiting time. A properly powered local test is ~4 h at N = 1 and ~14.5 h at N = 8.
**Do not read the eight held runs as evidence the fix works.**

### 3.10b The 0.90 ceiling is well below the unfolding transition, and TM4's loss is local (2026-09-10)

`reports/GroupMeetings/0323/group_meeting_03_23.pptx` slide 8 ("Phase transition (defolding)") measured
glpG's thermal transition under the **implicit membrane** model. The transition is sharp between
T = 1.05 and 1.10 -- mean CA-RMSD 17.6 -> 28.4 A and mean Rg 23.1 -> 35.7 A -- and plateaus by 1.15
(34.7 A / 42.5 A). At T = 0.90 the protein sits on the smooth pre-transition baseline at RMSD 10.0 A,
Rg 18.4 A. **T = 0.90 is therefore a legitimate ladder ceiling, far below unfolding.**

The dry-MARTINI hybrid agrees, which is a useful cross-model check. Placing the local wild-type runs
on the same axes (CA-RMSD to seed over t = 500-1000, protein-only Rg):

| run | T | CA-RMSD | Rg |
|---|---|---|---|
| fixed, T = 0.70 | 0.70 | 8.67 A | 20.17 A |
| fixed, T = 0.90 seed 1 | 0.90 | 8.38 A | 19.02 A |
| fixed, T = 0.90 seed 2 | 0.90 | 8.64 A | 19.97 A |
| unfixed, T = 0.90 | 0.90 | 9.91 A | 20.51 A |

All four land essentially on the implicit-membrane value at 0.90 and nowhere near the post-transition
28-35 A / 35-43 A. The fix also lowers RMSD slightly (8.4-8.6 against 9.9 unfixed).

**Consequence for the TM4 diagnosis.** The helix-fraction loss measured at T = 0.90 is **not** thermal
unfolding: RMSD and Rg are native-like and, in the run that lost the most helix (seed 2: TM4a
0.945 -> 0.556, TM1 1.000 -> 0.646), Rg *fell* 21.2 -> 19.5 A while the potential *fell* -20776 ->
-21192 E_up. That is the protein settling into a more compact, lower-energy, less-helical state, not
melting. So the residual TM4 problem belongs to the hybrid environment coupling or to helix propensity
in the bilayer (the GLY maps are a live suspect, 3.10a), not to temperature. An earlier suggestion in
this session to consider reverting the ceiling to 0.82/0.86 on the strength of the T = 0.90 helix
numbers was wrong and is withdrawn.

Two limits on how far the slide-8 result transfers: it measured RMSD and Rg only, so it cannot certify
TM4's *helicity* at 0.90, and it used the implicit membrane, so its transition temperature does not
carry over quantitatively to the dry-MARTINI hybrid.

### 3.10c What actually caused TM4 to be unstable, and what the GLY maps do (2026-09-10)

Six hypotheses were eliminated by measurement, in this order:

| hypothesis | test | verdict |
|---|---|---|
| global thermal unfolding | CA-RMSD 8.4-8.6 A, Rg 19-20 A; transition is at T = 1.05-1.10 (3.10b) | ruled out |
| GLY143 mid-helix alphaL flip | phi stays -67 to -88, h = 0.91-1.00 for the whole run | ruled out, it never flips |
| GLY133 / GLY49 cap sign | phi = +71 at GLY133 in the *healthy* T = 0.70 run (TM4b 0.991) | ruled out, no correlation |
| TM4a at the bilayer interface | \|z\| = 1.47 A, the most *central* segment (TM3 = 8.58 A) | ruled out |
| arm-R force-field tables | recorded blow-up reproduced to 0.27% on the older tables | ruled out |
| GLY Ramachandran mis-symmetrisation | controlled run, correct mirror vs buggy (below) | ruled out, correcting it is **worse** |

**Primary cause: the temperature mismatch (3.10, 3.10a).** The protein ran ~10% above its set point,
so the nominal 0.70-0.90 ladder drove it at ~0.77-0.99 T_up, about +28 K, which is across TM4's
fraying range while TM3 stays well below its own. With `inner_steps = 4` the cold end is now
*healthier than the seed*: at T = 0.70, TM4a 0.999, TM4b **0.991**, TM1 1.000 against seed values of
1.000 / 0.933 / 0.952. The pre-fix cluster gave TM4_full 0.589-0.727 at the same rungs, and 0 of 60
replica measurements passed the > 0.8 criterion.

**The GLY symmetrisation is a real defect but not this defect, and correcting it makes TM4 worse.**
Controlled test, wild type, T = 0.90, `inner_steps = 4`, two seeds per arm, t = 500-1000, everything
identical but the 23 GLY maps:

| GLY maps | TM4a | TM4b | TM1 | TM3 |
|---|---|---|---|---|
| buggy off-by-one mirror (what is deployed) | 0.786 | 0.840 | 0.810 | 0.819 |
| correct periodic mirror | **0.531** | **0.717** | 0.849 | 0.845 |

This is what the energetics predicted. For GLY143 the mirror penalty `E(alphaL) - E(alphaR)` is
**+0.571 E_up** under the buggy mirror, **0.000** under the correct one, and **-0.644** in the raw
library map. Only the buggy map favours the right-handed helix; the correct symmetrisation is exactly
neutral and the raw library actively prefers left-handed. So "fixing" the symmetrisation removes a
spurious ~0.57 E_up (0.6 kT) per-glycine helix bias that TM4 had been leaning on.

**Therefore the residual TM4 fraying at T = 0.90 is a force-field property, not a bug.** Glycine in
this Rama library carries no right-handed helix preference, and TM4 is the glycine-dense helix: 3 in
17 residues (1.76 per 10) against TM1's **zero** and TM3's 2 near its ends. TM4 is simply the marginal
helix, and it frays from its N-terminal turn (loss begins at 135/136, then 134) at the hot end of the
ladder, which is what a REMD hot end is for. Glycine genuinely is a helix breaker, so this may be
correct physics rather than something to repair.

**Statistical caveat, stated because the effect is not large relative to the scatter:** n = 2 per arm,
and within-arm spread is comparable to the between-arm difference (buggy TM4a = 0.990 and 0.583;
correct TM4a = 0.663 and 0.400). What the test establishes firmly is only the negative: **no run with
corrected maps beat the best run with the buggy maps**, so correcting the symmetrisation is not a TM4
fix. Deciding what the GLY map *should* be is a force-field question -- validate the library's GLY
dimer maps against PDB glycine statistics -- not a patch to apply to seeds.

### 3.10d The fix eliminates the blow-ups: 434x fewer bad frames (2026-09-11)

This is the measurement that was missing when the fix was deployed. Earlier attempts to demonstrate
blow-up prevention locally were underpowered by ~50x (3.10a); the production runs settle it. Counting
every production frame with a non-finite **or** positive potential across all 112 replica files
(4 variants x 28 rungs), skipping the rigid-protein equilibration group:

| | production frames | non-finite or positive | rate |
|---|---|---|---|
| pre-fix (archived `pre_tempfix_20260910/`) | 246 120 | 1 604 | **0.6517%** |
| post-fix (`inner_steps = 4`, ceiling 0.90) | 481 040 | **7** | **0.0015%** |

A **434-fold reduction on nearly twice the data**. At the pre-fix rate the post-fix runs would have
carried ~3 135 bad frames; 7 were observed. The logged rollback count agrees: 15 rollbacks in the
pre-fix block 1 against 0 in the 12 visible post-fix chunks. Note the ceiling was simultaneously
*raised* 0.82 -> 0.90, so this is not a temperature-lowering artefact.

**Not zero, and that is expected.** 7 frames remain, consistent with the arithmetic in 3.10a: at
`inner_steps = 4` the one-step-kick radius is 2.607 A, which still does not cover the ~2.43 A
approaches long runs reach. `inner_steps = 8` would cover that population at ~3.6x cost. The residual
rate is low enough that the rollback machinery absorbs it.

### 3.4 The four cluster POPE/POPG jobs were simulating a RIGID protein (findings 116)

The cluster HDX came out empty (188 of 203 amides off scale, resolved values to -53.9 kcal/mol). Not the
estimator, not the membrane term, not equilibration:

* **MBAR is healthy.** ESS 6529 of 48 624 (13.4%), top single-frame weight 0.0002, f_k spread 97.3, and
  ladder overlap *better* than the local run's (adjacent-rung mean gap / std 0.25 against 1.54).
* **Equilibration is not the cause.** The potential drifts -2.1 to -2.5 sigma over the run, but the last
  quarter drifts only -0.13 sigma and re-running on that tail alone gave the same degenerate profile.
* **The membrane term is not the cause.** Protein-only p_f == 1 exactly for 169 of 203 amides; adding the
  lipid term takes it to 171.
* **The protein has no internal dynamics.** In the *raw* hybrid trajectory, bypassing every analysis step,
  the internal CA RMSD between frames separated by whole chunks is **0.000-0.001 A**. Projected Rg is
  17.60 +- 0.000 and H-bond count 194.6 +- **0.01**, identical at T = 0.70 and T = 0.90. The coordinates do
  move (per-atom std 1.0 A), as a rigid body. Local control: Rg 17.67 +- 0.228, RMSD 2.46 +- 0.563,
  H-bond 169.7 +- 12.07.

**Root cause: `/input/stage_parameters.current_stage`.**

| | local seed | cluster seeds |
|---|---|---|
| `current_stage` | `production` | **`production_handoff`** |
| `activation_stage` | `production` | `production` |
| `preprod_protein_mode` | `rigid_body` | `rigid_body` |

`martini_hybrid.cpp:637-641`, `enforce_preprod_rigid_stage` returns `preprod_rigid && (stage != "production")`
and then calls `martini_fix_rigid::set_dynamic_rigid_groups`. `production_handoff` is not `production`, so
the protein is held as a rigid group. Both configs set `preprod_protein_mode = rigid_body`, so the only
thing separating a live protein from a frozen one is that string.

**Why the earlier verification missed it.** The record said "their stage is `production_handoff`, which
`martini_hybrid.cpp:646-647` treats as active, so the SC-env interface is on." That is true and it is the
wrong gate: `hybrid_interface_active_stage` accepts `production_handoff`, but `enforce_preprod_rigid_stage`
is a *different* predicate on the same string with the opposite accept set. When a stage string is checked
anywhere, enumerate every site that reads it, and verify the conclusion dynamically (does the protein's
internal RMSD change?) rather than by reading one gate.

The bilayer hole seen in the cluster trajectory and the -2.5 sigma potential drift are both downstream of
this: lipids relaxing around a protein that cannot respond. Fix is one attribute,
`set_stage_label(seed, "production")`, which `martini_prepare_system.py:1303` already does.

### 3.5 MBAR silently returns uniform weights for a hybrid coupled potential (findings 91)

`helpers/calc_hdx_ht.py` and `4.calc_D_uptake.py` built `beta[l] * cE0[k]` from raw energies with no
reference subtraction. For the protein-only potentials they were written against, O(1e2-1e3), that is fine.
A hybrid trajectory's `Energy.npy` is the full coupled-system potential, ~-7.6e3 E_up for a protein plus a
DDM micelle, so `beta*U` reaches -1.2e4, `exp(-u)` overflows, and the solver never leaves `f_k = 0`.

Measured on 48 states x 423 frames: raw gave f_k spread **0.000**, neighbour overlap **0.0000** and **0/423**
columns carrying weight at every target temperature except the bottom rung; mean-subtracted gave f_k spread
**71.07** rising monotonically, neighbour overlap 0.115-0.128, all 423 columns weighted, ESS 671-3378 of
20304.

The failure mode is the dangerous part: `f_k = 0` makes every weight equal, so the estimator returns an
unweighted average over the whole pooled ladder while reporting it as a reweighted ensemble at one
temperature. It raises nothing. Tell-tales are a gradient norm of exactly `sqrt(n_rep-1) * n_frames`,
`max_delta` of exactly 0, and ESS of exactly `n_rep * n_frames`. Two variants had already produced
plausible-looking dG plots this way.

Fixed in both files by referencing `cE0` to its pooled mean, which is exact (f_k shifts by `-beta_k*C` and
`exp(beta_target*C)` cancels in the normalisation) and leaves the protein-only path numerically identical.
`03.TrajectoryAnalysis/2.mbar_meltingCurve_freeEnergy.py` and `04.HDX/4.calc_HDX.py` have the same
construction but are byte-identical to master and only ever fed protein-only energies, so they were left
alone under master parity. Also note: passing the pymbar-3 style 3D `u_kln` under the installed pymbar 4.0.3
is *not* a bug; 3D and 2D give identical `f_k`.

**When a solver reports a gradient that is an exact function of the array shape rather than of the data, it
has not solved anything.** Check that before reading any number downstream of it.

### 3.6 Trajectory assembly: segments are not safe to concatenate

* **Production seeds carry an equilibration `/output`.** `run_remd.py` materialises replicas by copying the
  production seed, and those seeds already hold 300-400 frames from the handoff stage at a single
  temperature; on the first reseed that output rotates to `output_previous_0`, so the oldest chunk of every
  replica is equilibration. The concat joined oldest-first and `get_info_from_upside_traj.py` takes the
  temperature from the first frame, so all 48 replicas were labelled T = 0.8647: one state 48 times, exactly
  the degenerate uniform-weight condition of findings 91, and pymbar said so ("States 45 and 47 have the
  same energies") without failing. Fixed: the production temperature is read from the newest chunk and any
  chunk disagreeing is dropped and named.
* **Dropping destroyed chunks whole throws away most of a good trajectory.** One replica blew up 2333 frames
  into 2448, so a chunk-level rule would have cost 2332 good frames to remove ~100 bad ones, silently, while
  reporting a plausible frame count. Now filtered per frame on finite AND negative total potential. (A chunk
  also carries whole-run records that are not per-frame, such as `replica_swap_partner`; those are detected
  by comparing the first axis against the frame count and left out with a note.)
* The local run was clean while the cluster was broken purely because of how replicas were materialised
  (`warm_start.py` builds from a `seed.up` with no `/output`). **Testing one does not test the other.**

### 3.7 Analysis-code defects found in the shipped workflow

* **`k_chem` defaults to ~400 K.** `4.calc_D_uptake.py` defaults `legacy_T_range` to `[1.14]`, so exchange
  rates are evaluated at T_up = 1.14 ~ 400 K regardless of the trajectory's ladder. Base catalysis then
  dominates, `k_chem` reaches 4.9e5 s^-1, every amide is fully exchanged in milliseconds long before the
  first experimental time point at 60 s, and every normalised curve is the same step (a constant COF of
  27877.86 for all 63 peptides). At `legacy_T_range=0.85` (298 K, the experiment's rung) `k_chem` is
  ~10 s^-1 for a mid-chain serine and the 63 curves are distinct. Check `k_chem` in
  `<pdb>_percentD_feats.csv`, not the exit code. The failure surfaced 200 lines downstream as matplotlib's
  `TwoSlopeNorm ... must be in ascending order`, the only alarm the workflow raises for a degenerate COF
  set, which is why that exception was left unfixed.
* **`_DG_Hbond.png` free-energy scale is 15% low.** `calc_hdx_ht.py:337` forms `g = -0.593*np.log(hist)*t`
  with `t` in Upside reduced temperature. 0.593 kcal/mol is kT at 298 K, i.e. at T_up = 0.85, so the correct
  factor is `kB*t*350.588 = 0.6966*t`, and the shipped expression is low by 0.851 at every temperature. Left
  alone for master parity; the poster version is computed correctly in `make_hbond_landscape_figure.py`,
  which also drops bins carrying less than one effective frame (without that, the 245 K curve reads
  82 kcal/mol where the reweighting has no support at all).
* **COF is not the readable form of the uptake comparison.** It is the integral of the squared derivative
  of a curve normalised by `(max - first)`, so a nearly flat experimental curve is divided by a small span
  and inflates by orders of magnitude (experimental 21 to 9.1e4, simulated 129 to 2.5e4), and a linear R^2
  is then dominated by that normalisation. The rank correlation is the usable statistic (Spearman 0.30,
  p = 0.023, n = 57); rank within each dataset rather than putting both on one scale.
* Two latent incompatibilities in the uptake path: `5.analyze_D_uptake.py` sliced `<pdb>_<sim>_<i>_T.npy` by
  frame although it is written as a 0-d array (now `np.atleast_1d`), and `helpers/write_hybrid_energy.py`
  emitted a flat `Energy.npy` where the path indexes `[:, 0]` (now `reshape(-1, 1)`).
* **The Python engine computed a different model than the binary.** `engine_c_library.cpp` never called
  `load_masses_for_engine` / `register_fix_rigid_for_engine` / `register_stage_params_for_engine` /
  `register_hybrid_for_engine`, which `main.cpp` does, so any analysis through `upside_engine` evaluated a
  system with the hybrid interface inactive: **-12827 vs -18172 E_up** on the same coordinates. Every
  `get_output`/`energy` result taken through the Python engine before 2026-08-15 is suspect.
* When a node's arity changes, the migration is part of the change: a two-argument
  `martini_hybrid_position` cannot load a config declaring one, so every existing `.up` became unloadable
  the moment the binary was replaced. A config held open by a running job can only be migrated in the gap
  before the next block, so the migration has to live in that job's submit script.
* A silent fallback is worse than a missing input. Two from one session: a `0.0` belt half-thickness that
  made a solvation gate accuse a correctly inserted protein, and a metadata-PDB lookup hardcoded to
  `example/16.MARTINI/pdb/<id>.MARTINI.pdb` while the workflow writes it to `<run_dir>/hybrid_prep/`, so
  **every VTF the workflow had written labelled its lipids `UNK`** with positions intact, i.e. a trajectory
  that looked complete and was unselectable by lipid. `find_martini_metadata_pdb` now searches the run
  directories implied by the `.up` paths passed in, each derived from an explicit argument.

### 3.8 Two VTF-generation bugs, both fixed at the root (findings 126)

**Bug 1: the library never unwrapped molecules, so every VTF had torn lipids.** `extract_trajectory` wrapped
every particle into the box via `centralize_system` and then wrote the frame; nothing unwrapped, so any
molecule straddling a periodic face was left split across the cell, which VMD renders as bonds shooting
across the box. Measured on the same seed file, 400 frames, 4187 bonds:

| declared bond length | before | after |
| --- | --- | --- |
| mean | 6.89 A | 3.85 A |
| max | **141.00 A** | **6.91 A** |
| instances > 50 A | 54502 of 1674800 (3.254%) | **0** |

The 141 A worst case was a PO4-GL1 bond *inside one lipid*, so a protein-only integrity check passes while
the file is unusable. Fixed by adding `build_bond_walk` and `unwrap_molecules` next to `centralize_system`
and calling them in `extract_trajectory` after centring: walk the connected components of the declared bond
graph once, then apply a minimum-image displacement child-relative-to-parent per frame. No per-species
knowledge, so it covers protein, lipids and ions. **Check the declared bonds across all frames, not just the
protein backbone, when validating a VTF.**

**Bug 2: mode 1 on a hybrid system emits the protein twice, and VMD then rejects both copies.**
`build_mode1_mapping` emits the MARTINI-side protein (1050 atoms: N/CA/C/O x 210 plus 210 BB beads, all
named `PRO` by `infer_residue_names_from_class`, with zero bonds) *and* the appended all-atom backbone (840
atoms, real residue names, 839 bonds). VMD saw 210 unbonded pseudo-proline residues sharing resids 1..210
with the real chain and `atomselect protein` returned nothing usable; dropping only the duplicated N/CA/C/O
was not enough either, because the 210 BB beads still traced a ghost backbone under `not protein`.

This was never a bug to patch downstream: mode 1 is the wrong mode for a hybrid system.
`build_mode2_mapping` keeps `protein_membership < 0` (every environment particle, no protein particles) plus
the sequence-named all-atom backbone, and the library CLI already auto-detects it when
`input/hybrid_env_topology/protein_membership` exists; the analysis scripts were calling
`build_mode1_mapping` directly and bypassing that, and both were switched. Mode 2 reproduces byte-identical
atom and bond records to a hand-built "mode 1 minus the protein particles" and verifies in VMD (`protein` ->
840 atoms with names {C, CA, N, O} only, `name BB` -> 0). Coordinates can differ by exactly one box length,
since the unwrap anchor per molecule depends on atom ordering; each molecule is intact either way.

Selection notes: `lipid` returns 0 because `infer_residue_names_from_class` leaves lipids `UNK` (it only
assigns a lipid name for `cls == "OTHER"`), so select them with `chain X`, and do **not** relabel them
`DOPC` the way that function does, since this is a POPE/POPG system (the split is recoverable from the
head-group bead name, `NH3` = POPE, `GL0` = POPG). Residue numbering is trustworthy: `input/sequence`
position 79 is HIS or ALA and position 115 is SER or THR exactly as the variant names imply.

**The C-terminal carbonyl O is an unconstrained particle, a data property and not an extraction bug.**
Measured C-O distance over 180 frames: residues 1..209 mean 1.24 A and **max 1.24 A**, a rigid constraint,
while residue 210 starts at 1.17 A and escapes monotonically to 22.5 A. Exactly 1 of 210 residues is
affected and the N/CA/C backbone is unaffected. Harmless for the backbone physics and for HDX, but hide it
when rendering: `protein and not (resid 210 and name O)`.

---

### 3.9 Upside traps on exit under clang whenever Monte Carlo is enabled

`MonteCarloSampler` (`src/monte_carlo_sampler.h:12`) is abstract — it declares
`propose_random_move` pure virtual — but has **no virtual destructor**, while
`MultipleMonteCarloSampler` holds `std::vector<std::unique_ptr<MonteCarloSampler>>` and therefore
deletes `PivotSampler`/`JumpSampler` through the abstract base. That is guaranteed UB, and clang on
arm64 compiles the delete to a trap: the process dies with **SIGTRAP (exit 133)** in
`~MonteCarloSampler`, *after* the run has finished and flushed. GCC on midway2 does not trap, which
is why the same code trains fine on the cluster and fails on the Mac.

Bisected to `--monte-carlo-interval` alone: a single config with no REMD traps, and the same run
without that flag exits 0. Every local Upside run using MC has been exiting nonzero all along, and
`parameters/ff_2.1`-era code in `upside2-md-master` has the identical defect, so this is longstanding
rather than something this branch introduced.

The consequence was severe and silent in the wrong direction: ConDiv's worker checks
`j.job.wait() != 0` and raises `RUN_FAIL`, so **all 12 workers of a local training minibatch reported
`WORKER_FAIL` on completely valid data** — 250/250 frames written, Rg 13.2-13.8 A, potentials
negative, the temperature ladder correct. `run_minibatch` then raised `All jobs failed`. Read that
way round, an exit code was condemning good physics.

Fixed by giving the base class a virtual destructor. Verified by execution, not inspection: the
three previously-trapping invocations exit 0, and on an identical 160-time-unit run the old and new
binaries produce **bit-identical output across all 16 datasets** (`pos`, `potential`, `kinetic`,
`hbond`, `pivot_stats`, `rama_map_potential`, ...), so results are unchanged and master parity in
results holds. The trap sat purely in the teardown path.


## 4. What the hybrid gets wrong, and what is still open

### 4.1 Fold fidelity: a ~4.5 A helical-core deviation, cause unidentified (findings 100, 103)

The protein does not hold its tertiary helix packing. Helical-core CA-RMSD from the crystal **plateaus at
4.1-4.4 A** in POPE/POPG (3.48 -> 4.67 within the first segment, then flat across two more, so equilibrium
rather than drift).

| at T = 0.70, crystal-bonded amides | POPE/POPG |
|---|---|
| H-bond occupancy (median) | 0.844 |
| burial fails | 3.4% |
| both fail -> exposure | 3.34% |
| implied raw dG_open | 1.99 |
| helical-core CA-RMSD | 4.15 A |
| CA-Rg (crystal 20.43 A) | 20.77 |

**The detergent column is retired (2026-09-09).** These numbers were originally quoted against a DDM
micelle, which scored 2.61 A core RMSD and 0.952 occupancy. DDM is no longer an environment of this model,
so that column is not evidence for anything and is not carried here. It was never a clean comparison in any
case: that campaign died at block 2-3 so it had less time to drift, and its Rg was 1.8 A *below* the crystal
while POPE/POPG matches it, so it was compacted rather than more faithful. The consequence to keep in mind
is that **there is now no measured reference for how faithful this model can be**, only the crystal, so the
size of the deficit is stated against the crystal and nothing calibrates how much of it is recoverable.

**Ruled out by measurement (findings 100):**
* *Integrator.* The `avg_kinetic_energy/1.5kT` excess is +2-3% and dt-independent, far too small to produce
  a 15% unbonded population; covalent geometry is intact (worst C-N 1.78 A, 0 broken); the temperature
  dependence is weak (H-bond loss 27% -> 16% from 315 K to 245 K, ~1.7x, what a ~2 kcal/mol opening free
  energy predicts). The cross-check that used to close this bullet, the same integrator scoring 2.6 A core
  RMSD in a detergent micelle, is retired with DDM; the dt scan in section 4.2 is what carries the argument
  now.
* *H-bond assignment.* On the crystal geometry Upside's H-bond score agrees with the DSSP electrostatic
  criterion to within 8% inside helices (DSSP 86.5%, Upside 78.4%), and the 12 disagreements are marginal.
  ~16% of DSSP-helical amides are helix N-termini with no i-4 partner, so 86.5% is near the ceiling.
* *Lipid voids / bilayer prep.* Every backbone site has environment beads within 8 A (0 exceptions);
  nearest-environment median 5.2 A, max 7.3 A; coordination 8.86 within 8 A (9.59 in the TM belt).
* *Hydrophobic mismatch.* PO4-PO4 thickness 38.0 +- 0.1 A, acyl core 25.4 A against glpG's 28.2 A belt, a
  mismatch of only -2.8 A.
* *The burial threshold.* Burial failure is almost perfectly nested inside H-bond failure (3.4% vs 3.34% of
  frames), so the cut value is not what is binding.

**RD1: restoring the absent environment terms does not recover the fold (findings 103).** Three arms rerun
on the CB-corrected placement so the two effects were separable, at comparable step counts (750-780 k steps
each, single system, T = 0.70, all from one identical starting configuration):

| arm | restored | helical-core CA-RMSD (A) | Rg (A) |
|---|---|---|---|
| `base` | nothing (CB fix only) | 4.61 +- 0.09 | 20.36 |
| `env` | protein self-burial | 4.71 +- 0.12 | 20.43 |
| `envfull` | all non-membrane Upside terms | 4.53 +- 0.08 | 20.16 |

`env - base` is +0.10 A and `envfull - base` is -0.08 A, both inside the run-to-run scatter, so no arm
repaired anything. The -1.5 A improvement this was originally scored against came from the retired
detergent comparison; there is no calibrated target now. Rg stays at 20.2-20.4 against a crystal value of 20.43 in
every arm, so nothing over-compacted either: the predicted failure mode did not occur, but neither did the
intended repair. The CB correction also did not improve fold fidelity (`base` at 4.61 A is no better than
the 4.15 A previously measured, though the two are not directly comparable, 4.15 A being a 16-replica REMD
segment average and 4.61 A a single longer trajectory).

So the cause remains **unidentified**. Ruled out so far: the integrator, the H-bond assignment, lipid
packing and voids, hydrophobic mismatch, the burial threshold, the CB placement, and the absent
protein-protein environment terms. **Not** tested: the rotamer 1-body representation
(`placement_fixed_scalar`, fixed rotamer probabilities, versus standard Upside's rama-dependent
`placement_scalar`), which is the one remaining node-level difference from a standard configuration.
Consequence for the deliverable: do not add the environment terms to production.

### 4.2 Helical H-bond occupancy and non-cooperative opening

Measured directly on the raw per-donor scores (`get_protection_state.py --report-raw-data`):

| quantity | helix interior | loop |
|---|---|---|
| backbone H-bond occupancy, starting structure | **0.954** | 0.221 |
| backbone H-bond occupancy, trajectory | **0.681** | 0.226 |

The loops are unchanged; the loss is entirely inside the helices, ~27 percentage points of it, with 24 of
108 helix-interior donors below 0.5 occupancy. It is **not a thresholding artefact** (the raw score is
sharply bimodal: 60% of helix-interior frames above 0.5, 29% below 0.001, only 4.6% within [0.001, 0.05] of
the 0.01 criterion) and **not a starting-structure problem** (0.954 at t = 0).

The openings are also close to independent rather than cooperative. Cooperativity here is
`P(neighbour open | this open) / P(open)`, where 1.0 is independent opening; this metric **must be
normalised**, because the raw isolated-open fraction is higher in stock Upside (0.730) than in the hybrid
(0.528) purely because stock opens 5% of the time and the hybrid 29%. Never compare a count of isolated
events between two systems with different event rates.

| arm | protein KE excess | helix occupancy | P(open) | cooperativity |
|---|---|---|---|---|
| **stock Upside, no environment** | **+1.19%** | **0.950** | **0.050** | **4.18x** |
| hybrid baseline (dt 0.009) | +5.32% | 0.711 | 0.289 | 1.80x |
| hybrid, uniform max gamma | +2.01% | 0.763 | 0.237 | 2.12x |

So the hybrid opens the helical H-bond network **6x more often and roughly half as cooperatively** as stock
Upside on the identical construct. Real local unfolding is cooperative, a turn or segment opening together;
59% of the hybrid's openings are one amide alone with both neighbours still bonded, and that is what
produces log-amplified residue-to-residue scatter in a dG profile.

**The thermostat defect is real, is finite-dt error, and is not the cause.** `thermostat.cpp:31` makes the
global OU thermostat skip every atom with gamma > 0, so each protein backbone atom's only thermostat is its
own Langevin friction, set proportional to lipid-contact count (0 to 8.80, i.e. 0 to 50x the bare value)
while lipids are uniform at 0.169: buried helix cores drain slowest and run hottest, exactly where
protection should be highest. The g-JF propagator itself is correct. A dt scan at **matched physical time**
(270 t_up, steps scaled as 1/dt) settles the cause:

| dt | steps | protein excess | lipid excess | helix occupancy | cooperativity |
|---|---|---|---|---|---|
| 0.00900 | 30 000 | +5.32% | +1.02% | 0.711 | 1.80x |
| 0.00450 | 60 000 | +2.39% | -0.88% | 0.765 | 1.71x |
| 0.00225 | 120 000 | **+1.07%** | **+0.02%** | 0.758 | **1.53x** |

The excess falls by a consistent factor of 2.23 per halving and the lipids land exactly on target, so the
heat source is the discretisation and gamma's heterogeneity only matters because there is a source for it to
drain unevenly. **But dt does not recover the fold**: occupancy plateaus at ~0.76 against stock's 0.950 and
cooperativity does not improve. Therefore the excess opening and the non-cooperativity come from the hybrid
Hamiltonian itself, and no timestep, friction or thermostat change will remove them. This also retires the
"correction factor in analysis" idea: the defect is not a mislabelled temperature, it is missing cooperative
structure in the sampled ensemble.

Caveats to keep with these numbers: 270 t_up is short and occupancy is still decaying from the seed's 0.954,
so ~0.76 is an upper bound on how bad it is; and the stock arm has no membrane, so glpG's TM helices sit in
vacuum (Rg 16.6 A) which over-stabilises intramolecular H-bonds, making 0.950 an upper bound and the gap an
upper bound with it. The clean comparison needs stock Upside with its implicit-membrane terms, which the
hybrid topology does not carry (section 2.2).

Note also that PS hides most of this: 79% of broken-H-bond frames are still called protected because burial
exceeds `criterion3 = 5.0`, so helix-interior PS reads 0.936 while H-bond occupancy is 0.681. **When a
composite indicator (`A OR B`) is used on a system where `B` is nearly always true, the indicator stops
measuring `A`.** For a compact membrane protein nearly every backbone amide is buried by the protein's own
atoms, so quote the H-bond occupancy profile next to the dG profile.

### 4.3 Lipid packing against the protein is a seed property

The lipid-shielded fraction of amide N (>= 1 tail bead within 7.00 A, the criterion of section 5.3) differs
between two POPE/POPG datasets by almost a factor of two, and it traces to the seeds:

| dataset | seed, last frame | run |
|---|---|---|
| local 16-replica | **0.857** | flat 0.87 from the first block |
| cluster 48-replica | **0.433** | climbing 0.32 -> 0.50 over 14 chunks, plateauing near 0.50 |

Same protein, same 4949 atoms, same 279 lipids, same box to 0.5 A, same criterion, same code. It is **not
insertion and not a broken protein**: 87.1% of CA sit within 20 A of the bilayer midplane in the cluster run
and that is constant across all 16 chunks, Rg is 20.4 A (the crystal value), and consecutive CA-CA are
3.79-3.90 A with zero exceptions. The 12 CA reading 30-60 A from the midplane are the N-terminal tail
extended into vacuum, which is outside the bilayer in the local run too. The protein is in the membrane at
the right depth; what is missing is lipid *packed against its surface*. An HDX profile from the cluster data
will therefore show far fewer +inf amides than the local one for reasons that are neither the force field
nor the analysis, and the two datasets must not be compared as if they differed only in ladder size.

**Check the environment's contact with the protein, not just the protein's position in the environment.**
Insertion depth and Rg were both correct and constant here while half the protein-lipid contact was missing.

### TM4 needs the retrained tables AND the coverage nodes — neither alone is enough (2026-09-07)

Measured locally, four arms from one pristine glpG seed, three paired replicates each (seeds
1234/2345/3456), T=0.70, dt=0.009, 300 k steps (2700 time units). The trained force field was a
mid-training snapshot at step 269/500 of the ConDiv retraining of ff_2.1.

| arm | diverged | TM4 helix fraction | TM1 helix | Rg mean |
|---|---|---|---|---|
| CONTROL — ff_2.1, no coverage nodes | 0/3 | 0.441 [0.298-0.633] | 0.800 | 20.15 |
| ARM A — trained `pair_interaction` only | **1/3** | 0.588 [0.495-0.668] | 0.783 | 20.23 |
| ARM C — ff_2.1 tables + coverage nodes | 0/3 | 0.562 [0.318-0.766] | 0.814 | 20.30 |
| ARM B — trained pair + coverage nodes | 0/3 | **0.782 [0.657-0.863]** | 0.832 | 19.48 |

**The two causes are synergistic, not additive.** Restoring the coverage nodes with *old* tables
(ARM C, 0.562) and installing the *new* tables without the nodes (ARM A, 0.588) both land inside the
control's replicate spread — neither fixes TM4. Only both together (ARM B, 0.782) clears it, with
every replicate above 0.65 and the worst beating the control's best. TM4 goes from roughly half of
TM1 to parity with it.

This is what "the pair term was co-trained with coverage" predicts: the trained pair table is only
correct in the presence of the partners it was optimized against, and glpG sets those to zero.

**Corollary for the hybrid builder:** restoring `hbond_coverage` + `hbond_coverage_hydrophobe` is
worth doing only *together* with the retrained tables. That is why RD1 (findings 103) measured no
benefit — it restored the terms with old parameters, which is exactly ARM C.

**ARM A diverges deterministically on one seed.** Seed 1234 blew up at t=2490.8 in two independent
rounds with an identical signature: potential jumps -23462.7 -> -22142.0 (+1320 E_up in one frame
interval), then seven consecutive peptide C-N bonds across residues 136-145 (inside TM4) stretch to
2.10-2.73 A against a 1.33 A equilibrium, and the box goes NaN one frame later. Seeds 2345 and 3456
ran the full 2700 clean. 1/3 against 0/9 for the other arms is suggestive, not significant, but a
blow-up is a hard failure rather than a graded observable, and it starts in the TM4 backbone.

**Caveats.** n=3 with wide spreads (control TM4 ranged 0.298-0.633 across seeds); single-temperature,
single-replica-per-seed, no REMD; mid-training force field. `rama_map_potential` std varied 360-1955
across runs, far more than expected and **unexplained** — do not read that column as a health metric
until it is understood.

### A stopping criterion for the retraining, measured rather than guessed (2026-09-07)

`MAX_STEPS` was originally a guessed heuristic. Cumulative rms drift of `pair_interaction` from
ff_2.1, sampled across the run:

| step (approx) | cumulative drift | added since previous |
|---|---|---|
| 29 | 0.487 | +0.394 |
| 59 | 0.637 | +0.150 |
| 97 | 0.745 | +0.108 |
| 158 | 0.834 | +0.089 |
| 187 | 0.908 | +0.074 |
| 217 | 0.983 | +0.075 |
| 246 | 1.057 | +0.075 |
| 269 | 1.188 | +0.063 |

Against an initial rms of 2.965, 269 steps moved `pair_interaction` **40%**, `coverage_interaction`
45% and `hydrophobe_interaction` 41%. The drift decelerates sharply over the first ~60 steps and then
settles into a slow near-linear crawl of ~0.07 per 30 steps — it does **not** asymptote to zero, so
there is no natural convergence point. The first ~60 steps carry the bulk of the refinement; the tail
is the least productive part. Use this table, not a round number, to justify where to stop.

### Retraining reproduces how far it moves, not where — the between-run scatter is ~a third of the training signal (2026-09-08)

The GPFS outage split training across two hosts, which produced an accident worth more than the
inconvenience cost: two runs from a common ancestor at step 269 that took different routes to the
same step. midway2 continued from its own step-338 checkpoint; the rockfish lineage lost steps
270-338 and retrained them. Comparing their extracted `sidechain.h5` at step 355 and 354:

| array | drift_M | drift_R | between | between/drift |
|---|---|---|---|---|
| `pair_interaction` | 1.3902 | 1.3883 | 0.5449 | **39%** |
| `coverage_interaction` | 1.6664 | 1.6951 | 0.5587 | **33%** |
| `hydrophobe_interaction` | 1.3487 | 1.3415 | 0.5434 | **40%** |

`drift_*` is `rms(trained - ff_2.1)`; `between` is `rms(rockfish - midway2)`.

Two things are true at once, and only reading both gets it right:

* **The magnitude of training is highly reproducible.** The two runs drifted the same distance from
  ff_2.1 to within 0.1-1.7%. The drift also matches the table above (1.188 at step 269 to 1.39 at
  355 is the recorded ~0.07 per 30 steps), so both runs are on the same slow crawl.
* **The direction is not.** Separated by 0.39 of the distance each travelled, the two drift vectors
  differ by about 23 degrees. Relative to the tables themselves the disagreement is 15-17% rms.

So a single ConDiv run does not determine the force field to better than ~a third of what training
changed. **`MAX_STEPS` is not the only thing that needed a measured justification: the run itself
has a reproducibility scale, and it is large.** Consequences:

* A force field quoted from one run should carry this scatter. Two runs agreeing on an observable is
  evidence; one run is a sample.
* It is the reason the arm test was restored rather than skipped (`plan.md`). The question it answers
  is no longer "which coverage recipe" but "does a 16% table difference change TM4". If it does not,
  the TM number is real; if it does, 500 steps is not convergence for the deliverable.
* Never pair force fields from different steps when comparing runs. The step difference and the
  trajectory difference are the same size here, so a mismatched pair measures neither.
  `run_arm_test.sbatch` now refuses arm R unless `ff_3.0_trained_rf/STEP` matches
  `ff_3.0_trained/STEP`.

Not yet known: whether the 23-degree spread shrinks with more steps, or is the stationary noise of
the contrastive-divergence estimator. The drift table says the magnitude crawls without asymptote,
which argues for stationary noise, but that is an inference and has not been measured.

---

## 5. HDX: what the estimator measures and how to read it

### 5.1 The estimator is equilibrium, not kinetic

The `example/00.AnalysisScripts` uptake path does not read simulated elapsed time as exchange time. For each
amide donor, `get_protection_state.py` assigns a binary protected flag from backbone H-bond score, an
Asp/Glu side-chain-contact proxy and backbone/side-chain burial; `4.calc_D_uptake.py` MBAR-reweights those
flags to `p_protected`, computes sequence-, pD- and temperature-dependent intrinsic `k_chem`, then applies
the EX2-like `k_obs = k_chem * (1-p_protected)` and `D(t) = 1-exp(-k_obs*t)` in experimental seconds. So a
wrong friction or time mapping does not rescale the HDX time axis; it matters because it controls
decorrelation and the ability to sample opening/closing equilibria.

`dG` is `RT log(p/(1-p))`, not helix occupancy. At `T_up = 0.70`, 2, 3 and 5 kcal/mol mean roughly 98.4%,
99.79% and 99.9965% protection. Exact `p = 1` from a finite trajectory is censored, not a measured
1000-kcal/mol value, so a defensible plot must separate censored markers from finite dG points.

Two representation traps: donor IDs stored in `.resid` are zero-based while the VTF/PDB is one-based, and
`T.npy` is in Upside kT (values such as 0.85), *not* Kelvin, contrary to
`example/00.AnalysisScripts/README.md`; feeding 303 instead of 0.8647 invalidates MBAR.

The hybrid feeds this machinery through a projection rather than a fork: the adapter builds a protein-only
HDX-view H5 whose `/input` comes from the ordinary `-HDX.up` config and whose positions are N/CA/C mapped
through `hybrid_bb_map/atom_indices[:, :3]`, with the full coupled potential, temperature and H-bond logs
copied from the hybrid group, so the stock tools see their native `3*n_res` contract while MBAR still uses
the correct protein+bilayer Hamiltonian. The stock protein-only `PS.npy` is kept beside any combined PS so a
membrane correction stays observable and reversible.

### 5.2 Resolution limits, clips and sentinels

Master keeps two conventions on purpose: the **T-slice** uses the clipped step-6 path
(`mean_pf >= 0.99999 -> sentinel`, a ceiling of `0.001987 * 298 * ln(99999)` = 6.82 kcal/mol at T = 0.85) and
the **full-temperature** profile uses the unclipped jscripts path. They are not interchangeable.

The clip coincides with the estimator's statistical limit. With a binary protection state the smallest
resolvable `(1-p_f)` is ~1/ESS, so the largest supportable dG is `0.001987 * T_scale * ln(ESS)`. On the
CB-corrected ladder (85 760 pooled frames, ESS 10 697 = 12.5% at T = 0.85) that is 5.5 kcal/mol against a
clip at 6.8. Removing the clip let dG reach 19.6, but the values above the limit are noise:

| dG band | n | median effective frames carrying (1-p_f) |
|---|---|---|
| < 2 | 111 | 1526 |
| 2-5 | 58 | 45 |
| 5-8 | 12 | 1.4 |
| 8-12 | 3 | **1.00** |
| > 12 | 5 | **1.01** |

**100% of residues above 8 kcal/mol had their value set by fewer than two effective frames.** The jackknife
concealed it: `jk_pf` is clipped to `1-1e-6`, so the reported error is exactly 0.0 for anything above
~8 kcal/mol, i.e. maximum uncertainty displayed as maximum confidence. Measured per-temperature limits on
this ladder:

| T | resolution limit (kcal/mol) | residues at the bound |
|---|---|---|
| 0.70 | 4.33 | 84 of 203 |
| 0.75 | 4.81 | 69 |
| 0.80 | 5.13 | 58 |
| 0.85 | **5.49** | 27 |
| 0.90 | 5.57 | 28 |

An amide reaching the bound is **right-censored**: the data say "at least this protected". The bound is now
carried as a dotted line per temperature on the figure rather than by truncating the data, and
`plot_ref_style.py` renders the unclipped profile (`calc_hdx_ht.py`'s `_DG_res_T_slice.png` keeps master's
clipped convention untouched, for master parity; the npz stores both arrays so both renderings come from one
MBAR solve). Smooth 10-20 kcal/mol bands are not reachable from a binary protection state at any feasible
sampling, since dG = 20 needs `(1-p_f) ~ 2e-15`.

One plotting trap worth keeping: the out-of-range sentinels are `+1000.0` and `-100.0`, which **are
finite**, so `finite_mask = np.isfinite(...)` let every unresolved residue into the connected `errorbar`
series and the join to its in-range neighbours painted a full-height vertical line. Excluding them entirely
goes too far, since an excursion off the top of the axis is how a non-exchanging amide is conventionally
read; the settled behaviour draws the line over all values and draws error bars only where dG is resolved.
A sentinel encoded as a large finite number silently passes an `isfinite` filter, and here it had to be
excluded in two places because the same array feeds both the line and the markers. (Separately,
`plot_ref_style` now reports "the reweighting resolved nothing at this temperature" instead of crashing on
an empty array when every residue is off-scale.)

### 5.3 The membrane term, and how its criterion was calibrated (findings 113)

`get_protection_state.py` scores an amide protected only if H-bonded or buried by *protein*. A TM amide
facing lipid is buried by neither, so it drops out of the protected state whenever its backbone H-bond
flickers, even though there is no water there to exchange with. Master handles this with
`--use-TM-region`, reading the `surface` node; the hybrid HDX topology is a protein-only Upside config with
no membrane potential, so that node does not exist. The designed replacement is
`combine_hdx_protection.py --water-accessibility` fed by `py/martini_hdx_membrane_accessibility.py`, and the
driver was never calling it, so `PS.npy` was bit-identical to `PS_protein.npy`.

Measured on one replica (T = 0.844, 1350 frames, **unweighted**, so nothing here is estimator behaviour):

| bilayer-embedded helical amides (lipid-shielded in >=90% of frames) | at +inf | finite |
|---|---|---|
| protein-only PS | 35 of 92 | 57 |
| + lipid shielding | **79 of 92** | 13 |

Contiguous `+inf` runs go from 8, 6, 5, 4, 4, 3 ... to 17, 16, 15, 14, 9, 7: needles become blocks, which is
the reference figure's pattern arrived at from the same trajectory. Helix 104-126 is the clean example, with
the four residues that sit *outside* the bilayer keeping their finite values while 108-126 all go to `+inf`.

**Both numbers in the criterion are measured from the bilayer, not chosen:**
* **Radius = 7.00 A**, the flat first minimum of the intermolecular tail-tail g(r) (first peak 5.12 A). The
  minimum spans two 0.25 A bins whose ordering flips with noise, so the radius is the mean of the bins
  within 5% of the minimum, which gives 7.000 A on **all 16 replicas** where a bare argmin flips between
  6.88 and 7.12.
* **Threshold = 1 contact**, because at that radius a phosphate bead has a **median of 0** intermolecular
  tail neighbours (ester 2, first tail bead 6, terminal tail 11). The first tail contact is therefore the
  first step inside the head groups, i.e. the bilayer boundary itself.

`--cutoff` and `--min-contacts` are gone from `martini_hdx_membrane_accessibility.py`, so the criterion has
no free parameter.

A bare slab is the wrong shape and was rejected on measurement, not taste:

| criterion | helical amides at +inf / 148 | loop+turn at +inf / 55 | longest contiguous +inf runs |
|---|---|---|---|
| none (protein-only) | 41 | 0 | 8, 6, 5, 4, 4, 3 |
| slab between PO4 leaflet planes | 136 | 34 | 51, 46, 42, 22, 9 |
| slab between ester (GL1/GL2) planes | 117 | 22 | 42, 35, 22, 20, 19 |
| >=1 tail bead within 7.0 A | 115 | 10 | 22, 20, 18, 16, 13, 11 |
| >=3 tail beads within 7.0 A | 83 | 4 | 16, 15, 14, 12, 9, 7 |

A slab merges the helices into 40-50-residue blocks and erases the peak/valley structure, because glpG's
interfacial loops sit inside it too. Master does not use a bare slab either: `src/surface.cpp` selects TM
residues by `tm_min < z < tm_max` and ANDs that with a lipid-facing surface-exposure calculation, so a
residue lining the protein's own polar interior is not protected. The local tail-contact test is the
coarse-grained analogue of that conjunction.

With the membrane term in place the unclipped profile rises *continuously* into its excursions instead of
jumping, so the clip is not needed and only truncates. On the 16-replica ladder (177 k pooled frames,
unclipped, T = 0.85): 85 of 203 off scale with **81 of them in helices**, contiguous off-scale runs of 15,
15, 14, 12, 9, 8, a resolved range of -0.53 to 21.3 kcal/mol, helix-interior / N-cap / loop medians of
4.44 / 4.26 / 1.63, and exactly one helix-interior amide below zero.

What is left is genuinely the fold, and it is small: 13 of 92 deep helical amides stay finite at
2.3-4.2 kcal/mol. The measurements of section 4.2 stand, but they are no longer the explanation for the
figure: with the lipid term in place, a flickering H-bond on a bilayer-embedded amide is invisible to HDX,
which is physically correct.

**Update 2026-09-12: the off-scale plateau is a lipid-burial map, not a protection map.** Decomposing the
conjunction per frame over 171,382 samples (`protection_t = 1 - (1 - pp_t) * acc_t`, so an amide counts as
exchanged only when protein protection fails *and* it is water-accessible in the same frame):

| | `pp_fail` (H-bond/burial flicker) | `acc` | exchanged |
|---|---|---|---|
| TM1 30-48 | **0.0369** | **0.0021** | 3.13e-05 -> saturates |
| TM4 135-151 | **0.0322** | **0.4230** | 1.11e-02 -> stays finite |

The two helices flicker at the **same rate** -- TM4 slightly less -- and differ by ~200x in `acc` alone.
Per residue: **res 36 flickers 7.0%**, the worst of any amide, yet `acc = 0.0000` so it logs **0** exchange
events and saturates; **res 140 flickers 1.2%**, six times less, yet `acc = 0.195` so it logs **371** and
resolves at ~4 kcal/mol. Because `acc_t = 0` makes protection identically 1 whatever the H-bond is doing,
saturation is decided by lipid contact and helicity barely enters. **A `+inf` run is therefore not evidence
of secondary structure**, and the claim that TM1 "never exchanges" is wrong as a structural statement: its
H-bonds break 3.7% of the time and the lipid hides every break. This is the designed behaviour (a flickering
H-bond with no water present is correctly invisible), but it means the figure's protection map is only as
good as the 7 A tail-contact test, which asks solely "is a lipid tail nearby" and leans on the protein-burial
term to separate a protein-interior amide from a solvent-exposed one.

Within the TM4 core the exposure is a **helical face**, not the 21-residue depth gradient first recorded from
the full 131-152 window: `acc` for 140-147 runs 0.20, 0.39, 0.42, 0.016, 0.039, 0.654, 0.066, 0.001, so
141/145 (i, i+4) are both exposed and 143/147 both buried. One face is lipid-packed, the other points into
the protein interior and the catalytic cavity. **Open:** whether 141/142/145 line a genuinely water-filled
cavity (hybrid TM4 is then right, and better than implicit, which cannot represent a face at all) or are
packed against neighbouring helices (hybrid is then under-protecting them). Eight residues, so the face is
suggestive rather than settled.

**The same argument applies to any implicit-versus-hybrid comparison.** In the implicit model membrane
burial is part of the force field, so a burial-based protection state sees it for free; applying the
protein-only criterion to both strips the hybrid of protection the implicit model gets. (`--use-TM-region`
cannot even the two up: neither the legacy HDX topology nor the implicit configs carry a `surface` node,
which `upside_config.py --surface` alone creates.) Using each model's own membrane term inverts the result:

| T | n | Spearman | median implicit | median hybrid | hybrid never open |
|---|---|---|---|---|---|
| 263 K | 128 | 0.668 | 1.40 | 1.95 | 65 / 203 |
| 280 K | 130 | 0.691 | 1.28 | 1.91 | 61 / 203 |
| 298 K | 138 | 0.745 | 1.17 | 1.94 | 56 / 203 |

The hybrid is the **more** protective model, its helices saturate as they should, and the 298 K agreement is
better than the protein-only version (Spearman 0.745 against 0.676). Note the hybrid is also the more
censored of the two (56 of 203 unresolved against 32, ceilings 6.2 against 8.0 kcal/mol), so a sampling
contribution cannot be excluded. **"Hold the analysis fixed" is not "hold the physics fixed":** two models
can require different analysis in order to measure the same quantity.

### 5.3b Implicit vs hybrid: what each one actually calls "protected" (verified 2026-09-12)

Established by reading the scripts and the saved arrays, not the records, after a session in which the
`.md` notes were twice misleading on this point.

**The two protection rules are the same logical form.** Master, in `get_protection_state.py`:
`PS = HB1 + HB2 + BL`, then `if use_TM_region: PS += Su`, then `PS[PS>1] = 1` -- an **OR** over
H-bond, sidechain H-bond, protein burial and lipid-surface exposure. Ours, in
`combine_hdx_protection.py:31-32`: `exchange = (1 - pp) * acc; protection = 1 - exchange`, which is
`pp OR (not acc)` = `pp OR lipid-shielded`. So the hybrid's combination is master's, and it also
validates finiteness, the [0,1] range and shape equality. **The combination rule is not the problem.**

**But the implicit run never applied any membrane term at all.** Verified three ways:
`hdx_implicit.sbatch` calls `get_protection_state.py` bare; **nothing anywhere in the repo passes
`--use-TM-region`** (the only hits are its own `add_argument` and a docstring mention in
`martini_hdx_membrane_accessibility.py`); and the saved implicit results
(`/project/trsosnic/yinhan/implicit_79HIS_run3/results`, 48 replicas) contain only `_PS_protein.npy`,
`_Hbond.npy`, `_Energy.npy`, `_T.npy` -- **no `_ACC.npy` and no combined `_PS.npy`**. The implicit path
also never calls `calc_hdx_ht.py` or `plot_ref_style.py`; it runs its own inline pymbar block and saves
`_implicit_plain_pf.npy`. **So the lipid credit belongs to the hybrid, not the implicit run** -- the
reverse of the natural guess.

**Consequence, and it is the important one.** Because `acc_t = 0` makes protection identically 1
whatever the H-bond is doing, saturation is decided by lipid contact and helicity barely enters:

| | `pp_fail` | `acc` | exchanged |
|---|---|---|---|
| TM1 30-48 | 0.0369 | 0.0021 | 3.13e-05 -> saturates |
| TM4 135-151 | 0.0322 | 0.4230 | 1.11e-02 -> finite |

The two helices flicker at the **same rate, TM4 slightly less**, and differ ~200x in `acc` alone.
Res 36 flickers **7.0%**, the worst of any amide, but `acc = 0.0000` so it logs **0** events and
saturates; res 140 flickers **1.2%**, six times less, but `acc = 0.195` so it logs **371** and resolves
near 4 kcal/mol. **A `+inf` run is therefore a lipid-contact map, not a fold map**, and "TM1 never
exchanges" is false as a structural statement.

**The TM4 signal is real, and must not be suppressed.** It was proposed to mark TM4 water-inaccessible
because it is protein-buried. Rejected on measurement: the state-B frames (`pp=0 AND acc=1`) are
**temperature-activated** -- **zero** events across the 8 coldest rungs, 71-89% in the hottest 7, res 140
climbing monotonically 0 -> 123. A cutoff-margin artifact would appear at all temperatures. Also, `BL`
already protects these amides in 98-99% of frames, so an override would change only the 1-2% that *is*
the opening. And it would hard-code one protein's identity into shared `py/`.
**What the proposal did correct:** TM4 141/142/145 are **not** cavity-lining. They carry 23-24 protein
heavy atoms within 6 A -- indistinguishable from the deeply buried 143/147 (22.6, 22.3) -- and their
nearest lipid bead is a tail. They sit at **6.87-7.34 A** from the nearest tail against 4.90-5.16 A for
143/147, i.e. straddling the 7.00 A shell, which is why `acc` lands near 0.5. So they are protein-packed
in the hydrophobic region at the edge of tail contact, and the signal is local thermal opening.
Open: a cutoff **sensitivity check** at 7.5 and 8.0 A, reported as robustness, not a recalibration --
the 7.00 A is measured from this system's own g(r) and has no free parameter. Also open: res 143's 333
events are **non-monotonic** in T (76/93/83 at T=0.823-0.838, then 0-3 at the hottest rungs), which a
genuine exchange signal should not be. Do not quote res 143.

**Neither model has ever been compared to experiment.** `calc_hdx_ht.py:52-54` loads
`<pdb_dir>/<pdb_id>_{HXMS,NMR,NMR_MS}.csv` through `load_optional_numeric_csv`, and the `r_square` at
line 575 is computed only if one is present. Both `hdx/pdb/` and `hdx_postfix/pdb/` hold **only** the
`.pdb`; a search of `/project/trsosnic/yinhan` and `/beagle3/trsosnic/yinhan` for `*_NMR.csv`,
`*_HXMS.csv`, `*_NMR_MS.csv` and `*NMR_compare*` returns nothing, and no HDX log contains `r_square`.
**Therefore the Spearman column in 5.3's table cannot be agreement with experiment** -- there is no
experimental array in the project for it to correlate against -- and it must be the implicit-to-hybrid
correlation, i.e. how well the two models agree with *each other*. It was cited once in this session as
evidence of accuracy; that was wrong. Dropping a `<V>_NMR.csv` into `pdb/` makes the r-squared and the
scatter appear with no code change, which is the cheapest route to an actual accuracy statement.
Two cautions for that comparison: the implicit arrays are dated **2026-08-19**, predating the GLY,
rigid-stage and temperature fixes, so the implicit side must be re-run or the comparison confounds model
with three bug fixes; and with 56-84 of 203 amides censored the regression can only use the resolved
subset, which is biased toward the least protected amides.

**Which model makes sense, and it depends on the claim.** For lipid-dependent protection, cavity access,
PE-vs-PG or the RKRK variants, the hybrid is the only option -- an implicit potential is a function of z
and cannot distinguish two lipids that differ only in headgroup, nor represent the lipid-facing /
inward-facing asymmetry measured on TM4. For *intrinsic fold stability*, `pp_fail` is the right
observable and the implicit model is cleaner, cheaper, better sampled and internally self-consistent,
while the hybrid's lipid term actively obscures it. The counter-argument to keep in view: the upside
environment and burial terms were **trained against implicit solvent**, so the hybrid is a chimera
precisely in the terms that decide protection. The hybrid is a superset in practice -- it saves both
`PS_protein` and `PS` from one trajectory -- so state which of the two any figure quotes. Conflating
them is what made TM4 look paradoxical.

**Added later the same day, and it retracts a claim made above and on the 09/14 slide.** Between writing
the paragraphs above and the end of the session I argued, from the implicit model's many identical
near-zero rates, that implicit showed Englander's *cooperative-unfolding* signature while the hybrid
showed the native local-fluctuation one. **That was wrong**, for a plain reason: Englander's convergence
is onto the **finite unfolding rate of a cooperative unit**, whereas implicit's rates converge at
**zero**, which means uniformly frozen -- the opposite of unfolding. "Same value" is not "converged onto
an unfolding rate".

Measured directly instead, one replica per model at matched temperature (implicit T=0.851 of 48 rungs,
hybrid T=0.853 of 28; identical ladder span 0.700-0.900, mean 0.798, so this is not a ladder artifact),
counting how many interior amides of a helix are open in the **same frame**:

| open amides in one helix, same frame | implicit | hybrid |
|---|---|---|
| 0 | 77.2% | 69.8% |
| 1 | 16.1% | 12.8% |
| 2 | 4.9% | 9.7% |
| 3+ | 1.8% | 7.2% |
| **mean given >=1 open** | **1.39** | **1.99** |

Both distributions fall off monotonically with **no second peak at large counts**, so **neither model
opens a helix as a unit** and the all-or-nothing picture is wrong for both. Correlation
`P(i,j open)/(P(i)P(j))` by sequence separation 1..5: implicit `3.0, 2.7, 2.8, 1.6, 0.69`; hybrid
`4.3, 3.1, 2.6, 3.0, 3.0`. So **the hybrid is the more cooperative and longer-ranged of the two**, and on
the one-H-bond-at-a-time criterion **implicit is the closer match**. Figure: `fig_helix_opening.png` in
the 09/14 deck, generator `make_helix_opening_fig.py`, data `figs/coop_out.npz`.

Across all **108 helix-interior amides** of the ten DSSP helices, what actually differs is *freezing*,
not cooperativity: implicit has **65%** below one event per thousand frames and **35%** at exactly zero,
against **24%** and **9%** for the hybrid. But it is not uniform -- in **3 of 10** helices the hybrid is
the more rigid one, overall means are close (0.026 vs 0.034), and implicit fails outright where the
hybrid does not (res 126 **0.577** vs 0.051; res 150 0.212 vs 0.062; res 25 0.134 vs 0.014). So
"implicit keeps all helical regions completely rigid" is also **false** and should not be said; the
defensible statement is that implicit is close to all-or-nothing per amide while the hybrid breathes a
few percent throughout.

**Net position: do not rank the two models.** Implicit wins the textbook-mechanism comparisons
(one-at-a-time, short-ranged, clean frayed-terminus/rigid-core profile) and internal self-consistency.
The hybrid wins the only *quantitative* comparison to a measured number -- its TM4 core at 3.5-4.9
kcal/mol sits inside the 3-4 (poly-Ala) to 5-6 (Leu) range Langosch reports for TM-helix cores, whereas
implicit's frozen amides imply >5.4 kcal/mol at best and >7.8 if frames were independent, at or above the
top of that range. That bound depends on the **effective** number of independent frames, which has not
been measured; **an autocorrelation-time estimate on the implicit protection state is the one live
quantitative discriminator** and is the next thing to run. Both errors I made today ran in the same
direction, toward flattering the hybrid, which is worth remembering when reading any model-ranking claim
in this file.

### 5.4 Interpreting a per-residue profile

* **Do not use a bundled secondary-structure annotation as ground truth.**
  `hybrid_bb_map/bb_secondary_structure` is an idealised `CCCC1111HHHH...` pattern with helix segments of 41
  and 50 residues, which is not glpG's topology; DSSP on the same structure gives ten helices of 5-23
  residues. Re-scored against DSSP, 11 apparently broken helical donors resolve into 5 that are not in a
  helix at all, 8 that lie within four residues of a helix N-terminus (an amide in the first four positions
  has no i-4 carbonyl to bond to, so fast exchange there is what experimental HDX measures), and only **3
  genuine mid-helix breaks**. Those same three are the only ones of the eleven that were H-bonded in the
  prepared starting structure and lost it during the run, which is what makes the partition credible. The
  bundled annotation turned 2% of donors into an apparent 8% failure and sent me looking for a force-field
  cause for three days.
* **To decide whether a feature in a reweighted observable is physics or estimator, recompute it with
  uniform weights.** Anything that survives is the trajectory's. Of the 23 up-excursions at T = 0.85, 14
  amides never open in any of 86 333 raw frames and register zero transitions, so those are the trajectory's
  own statement; at the cold rungs most are estimator artefacts (78 at T = 0.70, of which 64 spurious),
  which is why the deliverable is plotted at the production temperature where reweighting is near-identity.
  Cold rungs earn their keep as ladder rungs for exchange, not as reporting temperatures.
* Helical donors read median dG 2.04 with 9% negative against loop donors' -0.05 with 51% negative, so the
  profile does track secondary structure.

---

## 6. Reusable diagnostic procedures

### 6.1 GLY Ramachandran maps: the rule, the check, and the failed fixes

**Rule: symmetrize ALL GLY maps unconditionally.** GLY has no beta-carbon, so its intrinsic Ramachandran
potential is symmetric under `(phi, psi) -> (-phi, -psi)`. Context-dependent maps trained on PDB data break
this as a statistical artifact: GLY in helical positions sees predominantly helical neighbours in the
database, so the raw map over-populates alphaL and places the global minimum at `phi_std ~ +85 deg` instead
of -85. GLY residues in TM helices then preferentially sample alphaL during simulation, breaking helix
H-bonds. This is the root cause of the repeated TM4 instability in glpG REMD. TM4 is the most vulnerable
helix (GLY132, GLY133, GLY136, GLY143, GLY149) and TM1's C-cap GLY49 is next; measured on a biased map,
GLY49 had alphaR 1.72 against alphaL 0.21 E_up and GLY133 alphaR 1.81 against alphaL 0.43 E_up, so both
helices were strongly destabilized. The initial hybrid backbone positions do come from the PDB and are
helical; the maps then penalise that.

Correct fix, applied 2026-09-04 to `py/upside_config.py`:

```python
for i, aa in enumerate(seq):
    if aa == 'GLY':
        m = rama_pot[i]
        rama_pot[i] = 0.5 * (m + m[np.ix_(idx_phi, idx_psi)])
```

Unconditional, no phi criterion, no alphaR/alphaL guard, applied in all three places: `write_rama_map_pot`,
and both GLY loops in `write_rama_map_pot2` (trans and cis variants). The `input_pos_override` parameter
added to support the discarded phi criterion was removed with it. Any future change to this block must
preserve the unconditional form.

**Four failed fixes, each of which failed silently:**
1. Using 1-indexed residue numbers as 0-indexed h5 array lookups: fixed non-GLY residues while leaving every
   GLY untouched.
2. A phi filter. First `[-130, -20]` deg, which missed GLY133 at -141.6; then `[-150, -20]`; then in Upside
   convention `[30, 160]`. Every range boundary misses some GLY residue. There must be no range.
3. An `alphaR energy > alphaL energy` guard, which likewise causes silent misses.
4. The wrong mirror formula. `m[::-1, ::-1]` maps index i to `n-1-i`, not `(-i) % n`, so it is off by one
   for all i > 0 on a periodic grid and creates a NEW asymmetry (alphaR preferred by ~0.58 E_up) instead of
   a neutral map. Worse, the `[::-1, ::-1]` symmetry metric reads 0.0000 on that broken map.

Two further checking traps: reading phi from initial hybrid positions requires the stride-4 backbone layout
(N=4i, CA=4i+1, C=4i+2, proxy=4i+3), and reading it as stride-3 gives garbage angles, so
`inject_backbone_nodes` reads actual N/CA/C from `hybrid_bb_map/atom_indices`; and grid indices must use the
`int()`/floor convention that `inject_backbone_nodes` uses, not `round`, since the two differ (25 vs 24 on a
72-point grid).

**How to check a seed file:**

```python
import h5py, numpy as np

fn = "path/to/seed.up"
with h5py.File(fn, "r") as h:
    rama_pot = np.array(h["/input/potential/rama_map_pot/rama_pot"])
    rama_resid = np.array(h["/input/potential/rama_map_pot/residue_id"])

n = rama_pot.shape[1]
mirror = (-np.arange(n)) % n  # CORRECT periodic mirror

for r1 in [49, 133]:   # 1-indexed GLY residues at TM1 C-cap and TM4 N-cap
    r0 = r1 - 1
    ridx = np.where(rama_resid == r0)[0]
    m = rama_pot[ridx[0]]
    sym_err = float(np.max(np.abs(m - m[np.ix_(mirror, mirror)])))
    # int() convention matches inject_backbone_nodes
    i_aR = int((-60. + 180.) / (360. / n)) % n   # = 24 for n=72
    j_aR = int((-45. + 180.) / (360. / n)) % n   # = 27 for n=72
    i_aL = int((+60. + 180.) / (360. / n)) % n   # = 48 for n=72
    j_aL = int((+45. + 180.) / (360. / n)) % n   # = 45 for n=72
    print(f"GLY{r1}: sym_err={sym_err:.4f} aR={m[i_aR,j_aR]:.3f} aL={m[i_aL,j_aL]:.3f}")
    # CORRECT state: sym_err < 0.001, aR == aL (within float precision)
    # BROKEN state: sym_err >> 0 (typically 3-4 E_up), aL << aR
```

Run it on every GLY row, not just 49 and 133. Before the 2026-09-04 fix all 23 GLY residues had symmetry
error 3.3-4.3 E_up and GLY132 had its global minimum at phi_std = +85. Remote h5 files already on the
cluster are fixed separately by `fix_gly_maps.py`, which symmetrizes by sequence lookup.

### 6.2 Verifying TM4 stability in a VTF trajectory

After any seed re-preparation, extract a short VTF and check:

1. **GLY133 phi** must remain in the helical range [-130, -20] deg over the trajectory. A drift to +60
   (alphaL) means the map is still biased.
2. **TM4 helix fraction**, residues 131-152 (1-indexed), with the criterion phi in [-130, -20] AND psi in
   [-90, +15]. Expect > 0.8 for a stable TM helix.
3. **The Ramachandran map check above**, run on the seed BEFORE submitting; every `sym_err` must be < 0.001.

A VTF has backbone atoms N/CA/C/O per residue for protein chain A, so atom indices start at 0 with
N, CA, C, O of residue 1, then residue 2. The phi angle for residue i uses C(i-1), N(i), CA(i), C(i). Local
verification after the 2026-09-04 fix (79HIS seed, 300 frames, T = 0.70) gave TM4 residues 134-151 at 1.000
helix fraction; GLY132/133 at the N-cap read 0.0, which is expected for helix-cap positions and not a
defect.

### 6.3 Detecting a frozen or rigid protein

Kinetic temperature, finite coordinates and retained secondary structure can all pass while the protein does
not move. Measure motion directly:

* **Internal CA RMSD between frames separated by whole chunks.** 0.000-0.001 A means a rigid body (section
  3.4); a live control gives 2.46 +- 0.563 A.
* **Standard deviations of projected observables.** Rg 17.60 +- 0.000 and H-bond count 194.6 +- 0.01 are the
  signature; a live run gives +- 0.228 and +- 12.07.
* **Drift-removed molecular COM motion per saved frame** for the lipids: 0.0081 A/frame was the frozen case,
  against 3.4-3.7 A net RMS when working.
* Read `/input/stage_parameters.current_stage` and confirm it is exactly `production`.

### 6.4 Health gates: count broken peptide bonds, do not use a worst-bond threshold

Measured on a forced NP tear (dt = 0.01), recording when each candidate criterion first fires:

```
FIRST FIRING TIME PER CRITERION
   maxCN (>3.5 A)                   t=168.0
   count (>=5 bonds >2 A)           t=168.0
   potential (non-finite or >1e6)   never
   |coord| (>1e4 A)                 never
```

* **`maxCN` is redundant**: the count fires at the identical frame, and the worst-bond threshold is the
  fragile one (healthy max 2.659 A against a torn 3.93 A is a knife-edge; at 2.5 A it false-fired on a
  healthy chunk and cost a 6 h glpG block).
* **An energy test can never guard NP** (431 broken bonds with the potential still finite at +3e5), while
  glpG's blow-up goes fully NaN within one 46-step interval, so there `isfinite` suffices. Different jobs
  need different detectors.
* **The count is a robust discriminant, not a tuned knob**: healthy frames have 0-2 stretched bonds
  (verified on all six NP systems and on healthy glpG), a torn one has 279-431, and any cut between 3 and
  200 behaves identically.
* Both drivers carry no invented magnitudes: `CN_MAX`, `POT_MAX` and `COORD_MAX` are deleted. NP fires on
  non-finite positions OR >= 5 stretched bonds; glpG on a non-finite potential anywhere in the chunk OR >= 5
  stretched bonds in the final frame.

### 6.5 Detecting a destroyed run when `isfinite` passes

* **`isfinite` is not a health check.** A finite 52 511-frame NP file was physically destroyed: protein Rg
  26.5 -> 76.7 A, max peptide C-N 56 A with ~300 bonds over 2 A, potential -12 000 -> +1e5. In glpG the
  environment coordinates at a failed frame reach +-4.65e12 A, numerically finite and physically destroyed.
* **The sign of the total potential is the cheap physical test.** A condensed bilayer plus protein sits near
  -2.2e4 E_up, so a positive total means core overlap. That is physics, not a tuned cut, and it is what
  `martini_remd_concat.py` and the reseed script filter on. In one capture the total went -2.17e4 -> +5.4e3
  -> recovered -> +1.98e6 -> NaN, oscillating for ~96 frames, which is why a single spot check can miss it.
  The peptide C-N count passes on those frames, because the protein is intact and the damage is in the
  lipids: the two tests catch different failures and both are worth running.
* **A collapsing adaptive chunk size is a free tell**: blown-up coordinates wreck the neighbour lists, so
  steps/second craters (394.9 -> 76.0 time units per chunk once). Gate a self-resubmit on a physical health
  check so a dead run stops instead of chaining.
* **Output-group order in a restarted `.up` is ascending**: `output_previous_0` is the oldest and `output`
  the newest, so a scanner starting at `previous_1` skips the oldest chunk.
* **Do not read a diagnostic's post-NaN lines as evidence.** The pair diagnostic then prints `max_force 0`
  with indices `-1`, because `NaN > max` is false and nothing updates the accumulator, and a repeated
  `min_dist 0.0010` is not a physical contact. Only the pre-NaN escalation is a measurement.
* **Term decomposition beats global observables for a local failure** (section 3.2).

---

## 7. System preparation

### 7.1 Environment morphology is derived from topology, not chosen

**Status (2026-09-09): DDM is retired as an environment of this model, and glpG runs in a POPE/POPG
bilayer only.** Everything in this section is kept, because it is the rule that stops any single-tail
detergent from being built as a slab, and DDM is the worked example the rule was derived on. Nothing here
is a live production path.

`derive_environment_morphology` counts acyl chains as connected components of the apolar (`C1`-`C5`) bond
subgraph in the lipid ITP: one tail means micelle, two or more means bilayer. DDM resolves to micelle,
DOPC/POPC/POPG to bilayer, and detergent+lipid mixtures are rejected. Deriving it from topology rather than
a name table or a CLI flag makes the unphysical combination unreachable. Note `parse_lipid_from_itp` already
returns **0-based** bead indices; subtracting 1 again silently split DDM's single tail into two and reported
it as a bilayer.

This matters because **a DDM lamellar slab cannot solvate GlpG** (findings 76, cited as "Update 76" from
`example/16.MARTINI/readme.md`). A CHARMM-GUI DDM slab has a tail core of 12.7-13.9 A against a 28-30 A
protein TM belt, and 50% of TM backbone CA had a polar maltose bead as their nearest detergent bead. TM4,
the most buried helix, was the only one that failed (helicity 0.90 -> 0.17 across the run, internal RMSD
4.45 A, against 0.84-0.96 and 1.6-2.5 A for the other five), and every residue that unwound sat at or beyond
the edge of that core while no residue inside it did. This is not tunable: lamellar thickness is
`2*V_tail/APL`, and `V(C12) ~ 324 A^3` means a 28 A core needs APL ~ 23 A^2, which a maltose head
(>= 40 A^2) cannot reach, which is why DDM forms micelles. Experiment settles which result is right: the
HXMS TM4 peptide 140-144 reaches only 50% deuteration at 24 h (dG_op ~ 10-11 kcal/mol), so a trajectory that
unwinds TM4 in a sub-microsecond segment is inconsistent with the data as well as with the implicit model.
Note also that `--membrane-thickness-angstrom` does not set slab geometry; it is only read for ion counting.

Building the micelle taught three things, each found by measurement: seed from the shell VOLUME rather than
from convex-hull support points (32 molecules against 186); fill innermost outward, since random-order
seeding lets molecules seeded far out block the contact layer; and do not take a packing distance from a
CHARMM-GUI step5 template, which is pre-minimization and contains bead pairs 0.24 A apart (the force field's
own `2^(1/6) sigma_max` = 5.276 A is the correct spacing).

**A packed-state thickness span cannot be gated on.** A gate comparing the environment's 5-95 percentile
tail-bead z span against OPM's hydrophobic thickness fails DOPC (20.4 A), POPE/POPG (20.6 A) and DDM
(11.0 A) alike against a 22.9 A limit, because CHARMM-GUI templates are laterally compressed and a clipped
percentile is not a relaxed hydrophobic thickness. What ships instead is `assert_environment_solvation` at
the production handoff, on equilibrated coordinates and per belt residue: hard-fail on vacuum (any belt site
with no environment bead within 2x the contact distance), and REPORT acyl-tail reach and local tail-core
thickness without gating on them, since on a post-damage snapshot both recover. **When a build-time gate and
the thing it guards are measured differently the gate is worthless**: measure both on the same state or
demote it to a report, and check a new geometric criterion against the paths it must NOT break.

### 7.2 A periodic tile's box is the tile

`prepare_bilayer_structure` parsed the template's CRYST1 into `bilayer_box`, used it only for tiling, and
then sized the actual box from the **lipid coordinate extent** with `force_square_xy=True`. For a
rectangular tile that squares the box up to the longer edge: an 84.41 x 73.10 A tile became 84.26 x 84.26, a
15% stretch along y that opens a vacuum stripe and moves the area per lipid from 61.70 to 71.0 A^2. Fixed by
using CRYST1 as the lateral box (`force_xy_box=target_xy`, `force_square_xy=False`) and letting
`force_xy_box` be an (x, y) pair; a template without CRYST1 is now a hard error. This is the same bug class
as an earlier "box inflated to 83.3 A for a 75 A tile", which had been patched by wrapping molecules into
the tile so extent ~= box. **When a fix works by making a wrong input look right, the cause is still there
and resurfaces on the first input that breaks the coincidence. Fix the derivation, not the input.**

### 7.3 Lipid bead models must match the ITP before anything else

POPE/POPG had never actually run through this pipeline, because the bead models disagree:

| source | POPE beads | tails |
|---|---|---|
| Robertson `last_frame.pdb` | 12: NH3 PO4 GL1 GL2 C1A **D2A** C3A C4A / C1B C2B C3B C4B | 4 + 4, unsaturation in tail A |
| CHARMM-GUI Martini Maker | 12: identical to the above | 4 + 4 |
| our `dry_martini_v2.1_lipids.itp` | 13: NH3 PO4 GL1 GL2 C1A C2A C3A C4A / C1B C2B **D3B** C4B C5B | 4 + 5, unsaturation in tail B |

Different bead count, different tail lengths, different unsaturation position: not a renaming. DOPC matches
exactly between CHARMM-GUI and our itp, which is why DOPC and DDM worked while POPE/POPG silently never did.
Adding the 12-bead model to the itp is rejected (it means running a lipid the force field was not
parameterised for) and substituting DOPC is rejected (the result is about composition), so the ITP is the
authority and the coordinates must be generated for it. **A lipid being present in the ITP does not mean the
pipeline can build it**, and one `diff` of the two bead-name lists would have found this in the first minute
rather than after the orientation, cropping and composition work.

Two facts from writing that builder. **Rigid-rod conformers cannot start at the target area**: nearest
intermolecular bead distance was min 0.20 / median 2.05 A at APL 69.4 against the reference bilayer's
4.21 / 4.60 A at APL 61.7, because at 69.4 A^2 the lipid axes are 8.33 A apart while a rigid conformer is
~6 A wide at every height. Four attempts inside the rigid-rod picture failed; the route that works is to
start deliberately LOOSE and condense with a pure-bilayer run under the xy barostat, then use the
equilibrated tile as the template (valid only for a pure bilayer, since with a protein the box is pinned by
its footprint). And **compute the comparison metric on THEIR data before touching your own builder**: four
iterations went into tuning against intuition with the reference bilayer already in hand. A deposited
composition may also not be the nominal one, 1814 POPE : 959 POPG being 1.892:1 rather than 2:1.

### 7.4 Ions and box size (NP-1AO6)

Final: **neutralizing counterions only, no bulk salt** (user decision 2026-08-04). 218 K+, 0 Cl-, cancelling
MPA (-203) + protein (-15), for 4198 particles; `salt_molar`, `estimate_salt_pairs` and the free-salt
assertion are gone from the NP path.

Salt molarity was never the defect (every build measured 0.148-0.150 M). The defect was **box volume**:
`box_len` was derived from the rotated protein's reach *about the gold COM*, and because the frozen NP is
pinned at the box centre while the protein adsorbs to one side, that doubles the protein's lever arm, so
boxes came out 232-284 A for a complex only 122-148 A across and correct ion density times inflated volume
gave too many ions. Now complex-centred at a fixed 200 A. It survived repeated rebuilds because nothing
asserted the *built* composition, so every rebuild re-derived the count correctly from an unexamined
premise; `build_system` now asserts exact charge neutrality and zero salt pairs from the placed ion counts
and rejects a box too small for the complex, all negative-tested. The 200 A box is itself too small for the
states the model visits (two orientations exceeded it at Rg 152 and 206 A) and the box is near-vacuum, so
enlarging it costs almost nothing.
---

## 8. Bilayer physics measured with this force field

Static structure validates the force field and is the strongest positive physical claim available
(g-JF, 128-DOPC, T = 0.8647, 1500 frames): APL 65.6 A^2 (experiment ~67.4, -3%), P-P thickness 43.3 A
(experiment 37-39, +10-15% thick), chain order P2 tail-average 0.30 (MARTINI DOPC 0.2-0.3), tilt
22.5 +- 12.3 deg, director S = 0.74. The bilayer is fluid but somewhat over-ordered and thick. Report that,
do not twist it. The P-P excess is **headgroup projection**: the hydrophobic core d_c (C1A/C1B
leaflet-leaflet) is 27.8 A = 2.78 nm, matching DOPC's 2.7-2.9, so the functionally important
TM-hydrophobic-matching thickness is correct.

**Transport is not claimable.** Diffusion, reorientation and flip-flop are cage-escape processes needing
physical friction, of which the g-JF has almost none. MSD(lag) on a long bilayer run has a local exponent
falling from 0.51 to 0.27, i.e. toward caging rather than Fickian, so no diffusive plateau exists on this
window and any apparent match of D to 11.5 um^2/s at one lag is a coincidental crossing on a falling curve.
An overdamped control on the same box is also sub-diffusive (alpha 0.35), so this is a property of a small
128-lipid patch plus the real short-time regime, and a curve-match of the CG director rotational ACF to
all-atom CHARMM36 overlays poorly (rms 0.14). **No single effective-time factor maps CG lipid time to
physical time.** What is claimable is correct thermodynamics with a nominal sampling time (which is what
REMD needs) and correct static mechanics; not any single-scalar transport timescale.

**Elastic moduli need a fluctuation-correct barostat.** `box.cpp`'s "Parrinello-Rahman" path is a damped
relaxation scheme (0.95/step box-velocity damping plus tight scale clamps), not an extended-Lagrangian
barostat, so the area barely fluctuates (sd 0.95 A^2 against ~81 expected at K_A ~ 265) and K_A from
fluctuations is nonsense. The Monte Carlo barostat added for this (`BarostatType::MonteCarlo`, type 2) is
exact-NPT Metropolis scaling molecule COMs via `/input/molecule_ids`, so dU is intermolecular; it restores
fluctuations (sd 0.95 -> 29) but K_A is still undersampled, because area relaxation is slow and coupled to
the sub-diffusive lipids. Mean APL under a barostat equals the NVT APL, so the FF's zero-tension APL is
~65.5.

Three engine facts established alongside. REMD was verified correct for these lipids (per-slot T reaches the
lipid integrator at the thermostat cadence, the lipid potential is in the swap criterion, no momentum
rescaling is needed under a coordinate swap, and Arrhenius gamma(T) only sets timescale), and a sweep from
T_up 0.5 to 1.5 is stable at every rung. One latent hazard, flagged and not fixed: the `mv` RESPA integrator
overload does not call `apply_brownian_step` nor skip `brownian_mask`, so `mv` plus brownian lipids would
silently mis-integrate them; MARTINI uses `v`, so it is not triggered. And `effective_time_factor` is NEVER
read by the engine, it is analysis-only metadata by design.
---

## 9. The nanoparticle campaign (1AO6 + MPA-AuNP)

The pre-fix campaign is not interpretable and its site conclusions are retracted. It predated **two**
simulation fixes, not one:

| | old NP config | corrected | where |
|---|---|---|---|
| LJ core table inner knot | `r_min_ang = 0.00` | **0.30** | section 3.1 (findings 92/93) |
| CB placement (centroid-relative) | `[0.0000, 0.9438, 1.2068]` | **`[-0.0198, 1.5117, 1.2068]`** | section 2.3 (findings 102) |
| `martini_hybrid_position` arity | 2 args | 2 args | migrated, OK |
| `current_stage` | `production` | `production` | OK, not rigid |

Both bear directly on what the campaign measures: the old table is force-free at short range and particles
were shown to reach it, which on an adsorbing surface is exactly where they go; and the CB offset displaces
every sidechain site by 0.568 A, and `martini_sc_table_1body` anchored there is both the term that drives
adsorption and the site at which the footprint is scored. A re-run needed a rebuild, not a patch, since the
tables have to be regenerated. **A migration script fixes what it says it fixes**: the arity migration
passing had been read as evidence the NP configs were current, and they were two generations behind.

The rebuilt run is a different simulation in every respect that matters:

| | pre-fix campaign | rebuilt |
|---|---|---|
| Rg, median / max | 85.9 / **209.0 A** (exceeds the 200 A box) | **48.3 / 78.2 A** |
| adsorbed **and** compact | 236 of 12 312 = 1.9% | **1736 of 25 924 = 6.7%** |
| dominant orientation's share of that window | **71%** | **26%** |
| residues in contact per compact frame | **105.5 of 578** | **16.8** |
| residues with contact frequency > 0.3 | 136 | **2** |

The pre-fix "footprint" was the protein smeared over the particle. **A contact-frequency footprint needs a
sanity check on how much is in contact, not only where**: 105 of 578 residues touching a 5 nm particle
should have been read as a smeared protein, and that one number would have invalidated the site ranking a
week earlier.

In the rebuilt window none of the paper's five lysines is contacted (K12 0.000, K73 0.004, K190 0.000,
K525 0.000, K541 0.029), where the pre-fix run gave K525 0.542 and K541 0.678. Only the K190 result
survives, and unchanged: 0.000 in both. What the rebuilt run contacts instead is centred elsewhere (Lys313
0.34, Glu311 0.31, Asp314 0.28, Asp562 0.25, Lys560 0.25). Re-measured at block 2: unchanged, with Rg median
risen 48.3 -> 62.7 A, so albumin is still unravelling. **This is provisional and must not be quoted**: the
highest per-residue contact frequency is only 0.341, so no pose is yet preferred.

Method notes that remain valid:
* **Choose the analysis window before looking.** Once albumin spreads, contact discriminates nothing (core
  and surface mean contact 0.272 vs 0.261); restricted to adsorbed and still-compact frames
  (Rg < 1.25x native) it separates properly (0.128 vs 0.245). The experiment labels albumin whose CD still
  shows its secondary structure, so that window is the only comparable state, and orientations reaching
  Rg 152-172 A in a 200 A box self-interact through the periodic image and are unusable outright.
* **Score contact at CB**, because that is where `martini_sc_table_1body` anchors. The reconstruction
  (Kabsch fit of the stored `affine_alignment` reference onto N/CA/C, then the fixed CB offset) matches the
  engine to 1.5e-4 A. A CB cutoff under-counts lysine by ~6 A of sidechain, so rank lysine against lysine.
* **Replicates beat length.** The compact fraction fell from 3.2% to 1.9% as trajectories grew, so running
  longer moves away from the informative window; many shorter independent runs sample it better at the same
  cost. Pick a test orientation by its historical onset too: cumulative time to first tearing differed by
  16x across the six faces (90-0-0 t~115 up to 0-0-0 t~1873), so a passing run on 0-0-0 proves almost
  nothing while the same run on 90-0-0 is the cheap decisive test.
* Selecting `residue_ids == r` without the `particle_class == "PROTEIN"` mask picks up GOLD/MPA/ION beads
  that share residue numbering, which once reported a protein-gold distance of exactly 0.00 A. And "distance
  to gold" taken at the argmax-C-N residue per frame is meaningless while the argmax wanders: pin the
  residue first, then track it.

---

## 10. Cluster and operational lessons

* **A wedged GPFS makes a dead job look healthy, and `squeue` will not tell you (2026-09-07).** Job
  48981235 was reported `RUNNING` for 3.5 h while all nine of its workers sat in `D` state at
  `00:00:00` CPU, wchan `cxiWaitEventWait` / `lookup_slow`, having never started their compute
  binary. The honest probes are per-process, not per-job: `sstat -a -j <id>` (a step whose `AveCPU`
  does not climb is not computing), then `ps -o pid,stat,time,etime,wchan` on the allocated nodes.
  `D` state is uninterruptible, so `timeout 10 ls <wedged dir>` does **not** return — it leaks a
  process and, over an SSH ControlMaster, burns a session channel until the mux refuses new
  sessions and ssh falls through to password auth. That fall-through is the RCC-ban trigger, so pin
  `-o BatchMode=yes -o PasswordAuthentication=no -o NumberOfPasswordPrompts=0` on every cluster
  call before probing anything that might hang.
* **RCC's GPFS serves midway2, midway3 and beagle3, so "try the other cluster" is not a fallback for
  a storage incident.** During the 2026-09-07 outage midway2's login nodes refused TCP while
  midway3's login nodes had lost `/home`, `/project`, `/project2`, `/scratch` and `/software`
  outright — `stat -f /project` reported **xfs**, and with `/software` gone there was no `squeue` or
  `sbatch` in PATH at all. Distinguish "this node's mount is stalled" from "the filesystem is gone"
  with `stat -f`, and check a second login node before concluding either.
* **All four `midway3-login[1-4]` share `midway3.rcc.uchicago.edu`'s host key, and skipping that
  detail costs a Duo push.** Connecting by the per-node hostname stops at an unknown-host-key
  prompt; an expect script then answers *that* prompt with the password and no push is ever sent,
  which is indistinguishable from a Duo failure. Pass
  `-o HostKeyAlias=midway3.rcc.uchicago.edu` (`scratchpad/rcc_master.exp` does) instead of editing
  `known_hosts`.

* **A Slurm NODE_FAIL requeue silently eats the REMD block budget (findings 125).** `run_remd.py` derives
  its block number from a plain `block_count` file incremented at **every process start**, and nothing
  decrements it when a start produces no data; the chain stops once `blk >= MAX_BLOCKS` (12). One job was
  requeued 5 times, every time `NODE_FAIL` on a node `scontrol` reported as `ALLOCATED+NOT_RESPONDING`, each
  incarnation burning ~15 min of calibration and dying, so after 4.5 h the variant was at block 6/12 having
  completed **zero** chunks while its three siblings were at block 1/12 with 2 chunks each. It is hard to
  see because `squeue` shows the job `RUNNING` and `sacct -X` shows only the newest incarnation (use
  `sacct -j <id> --duplicates` or `scontrol show job <id> | grep Restarts`), because the log is truncated on
  each requeue, and because `grep -c "chunk done"` is the only honest progress metric. Rules: reset
  `block_count` to the genuinely completed blocks after any requeue; exclude the failed node (`--exclude=`
  in the submit script, so it propagates to self-resubmissions) or Slurm re-allocates it indefinitely; and
  treat `NOT_RESPONDING` as unproven death, confirming static log size and h5 mtimes over a dwell and finite
  `input/pos` in every replica before resubmitting. An abrupt kill loses unflushed HDF5 buffers, which is
  safe here: the file reverts to its last consistent state and `reseed()` resumes from it.
* **`~/cds3` is `/cds3/trsosnic/yinhan`, which compute nodes cannot read.** A job reading it fails with
  `FileNotFoundError` on every file while the login node lists them happily. Stage to `/project` first;
  `~/project` is a symlink to `/project/trsosnic/yinhan` and works from compute nodes.
* **Never pipe a script that performs writes into `head`.** `python3 reseed.py <16 files> | head -3` made
  `head` exit after three lines, the pipe close, and the producer take SIGPIPE, so only the first three
  replicas were reseeded and thirteen re-simulated ~41 000 steps. Use `tail`, which drains its input.
* **A single-replica smoke test bounds loading and setup, not stability.** A rare non-recovering force spike
  needs ~1e4 steps across 48 replicas to appear; a 400-step seed test could never have caught it.
* **Compute the expected event count BEFORE running an A/B on a rare stochastic failure.** Two tests were
  spent on the reaction-field question without discriminating power, one confounded by ongoing structural
  relaxation and one 10x under-exposed (the observed failure rate was ~2 events per 2e6 replica-steps while
  each arm was 1.9e5 replica-steps, i.e. ~10% chance of a single event even in the defective arm).
* **`--thermostat-interval -1` does NOT mean NVE.** `main.cpp` computes
  `thermostat_interval = max(1., round(arg / (inner_step*dt)))`, so -1 clamps to **1** and the thermostat
  fires every step. Notes describing it as "effectively NVE" were wrong.
* **Recovery is not a guard.** A supervisor that runs the ladder in rounds against a wall-clock deadline and
  reseeds a destroyed replica between rounds clamps nothing, widens no threshold, and drops the destroyed
  frames rather than repairing them. That is the same rollback-and-continue design `run_remd.py` uses.
* Removing whole Python functions programmatically: a naive "def line -> next column-0 line" scan mis-cuts
  multi-line signatures whose closing `)` sits at column 0. Use `ast` (`node.lineno..node.end_lineno`), or
  verify with `py_compile` after each bulk removal.
* Check the units and axis directions of any workflow figure before showing it. Three shipped-analysis
  presentation defects turned up while building one poster: the `_DG_Hbond.png` scale (section 3.7), the
  ESS-censoring confusion, and `_Tm_curve.png`'s inverted hydrogen-bond axis.

---

## 10a. TM4 is flat between training steps 269 and 404 (measured 2026-09-09)

ARM_B repeated locally on the step-404 force field, paired with the step-269 four-arm test: same
seed file (`glpG-RKRK-79HIS.up`, md5 `6a8285d1...`), same protocol (300k steps, dt 0.009
hard-locked, T=0.70, `--disable-recentering`, frame-interval 6.75), same three RNG seeds.

| force field | TM4 mean helix fraction | Rg mean | diverged |
|---|---|---|---|
| ff_2.1, no coverage | 0.441 [0.298-0.633] | ~20.4 A | 0/3 |
| trained step 269 + coverage | 0.782 [0.657-0.863] | 19.48 A | 0/3 |
| **trained step 404 + coverage** | **0.709 [0.641, 0.812]** | **19.62 A** | **0/3** |

Per replicate 0.641 / 0.812 / 0.673; TM1 0.823 mean. **135 further training steps bought no
measurable TM4.** The ranges overlap heavily and n=3 at one temperature cannot resolve 0.07, so the
claim is "flat", not "worse". It is worth knowing because the force field itself moved a lot over
those steps — |dpair| 5.65, |dcoverage| 9.34 against init-to-269 magnitudes of 12.40 and 20.25 — so
the parameters were still changing while this observable was not.

**Rg stays ~0.8 A compact** (19.62 vs the 20.4 A crystal), unchanged from 19.48 at step 269. The
over-burying risk from training the environment against implicit solvent has neither resolved nor
worsened.

**One replicate had a full recovery from a large excursion.** s1234 reached potential **+11269 E_up**
from a -24948 start, stretched one peptide C-N to **9.83 A** (17 bonds over 2 A at some point) and
logged `avg_KE/1.5kT` **1.117** against a 1.000 target, then returned to -23739 with mean C-N
1.320 A and zero stretched bonds in the final frame. The other two stayed clean (max C-N 4.97 and
3.17 A, KE 1.030 and 1.042). This is the signature class of the documented blow-up mechanism
surviving rather than propagating, so treat a single such excursion in production as a warning, not
proof of failure.

**Two traps in the local test harness.** `analyze.py` ignores argv: it hardcodes four arm names and
resolves configs from `run2/` and logs from `logs2/` relative to its own file, so it must be driven
by staging that layout (symlinks are enough) and overriding `ARMS`; its final per-arm block then
crashes on a format string that assumes four arms, while the per-run table is complete. Separately,
**Upside overwrites `/output` rather than appending**: a fresh run on a production seed replaces the
seed's frames, and the way to tell them apart is the time spacing, not the frame count — the seed
carried 300 frames at 0.45 spacing and the run wrote 401 at 6.75, so an unwary "skip the first 300
frames" would have discarded three quarters of the new data.

## 10b. A ConDiv checkpoint can be rebuilt from an extracted force field

Measured 2026-09-07, when the outage left the step-269 `sidechain.h5`/`environment.h5` as the newest
reachable force field and no checkpoint at all.

`expand_param` writes five of the six `unpack_params` blocks to `sidechain.h5` and **discards the
sixth (`rotscalar`)**, so the h5 is not a complete parameter record. It is still enough:

* `pack_param` (`py/rotamer_parameter_estimation.py:257`) is an L-BFGS-B refit rather than an
  analytic inverse, but on trained tables it is effectively exact — final loss **1.95e-18**,
  reproducing pair/coverage/placement/centre to **1.1e-16** and hydrophobe to 1.4e-9 (5.6e-11
  relative, against a trained signal of 10-35 in float64 tables). The palindrome floor of 54.18 seen
  at init does **not** reappear, because the GLY row of a trained table is already palindromic.
* `rotscalar` is identically zero at init (`|x|max = 0.000000`) and is not part of the deployed force
  field, so borrowing it from the init checkpoint is exact for what the simulation reads.
* `params.env = energies[:, :-1]` recovers exactly, because `expand_param` sets the last column to a
  copy of the third-from-last — assert `energies[:, -1] == energies[:, -3]` to confirm.
* Adam state is **not** recoverable. `alpha` is constant (rot 0.125, env 0.025) and the bias
  correction applies from step 1, so a fresh solver costs a transient rather than a mis-scaled step;
  with `beta2 = 0.96` the second-moment memory is only ~25 steps, so the transient is short.

Acceptance test that matters: run `extract_ff.py` on the rebuilt checkpoint and diff its output
against the force field you started from. Builder: `scratchpad/ff3_retraining/build_local_resume.py`.

**The stale worker copy is the trap here.** `state['worker_path']` in a local checkpoint may point at
`run_output/ConDiv.py`, which on this Mac predates the env-derivative fix (it calls
`get_param_deriv(env_shape, ...)` with 360 elements against a node returning 760) and the
`environment_potential_type = 0` fix, so it would build a sigmoid environment node and then fail in
`compute_divergence`. Point `worker_path` at the all-fixes `training/gly-sym/ConDiv.py` and assert
the fix markers are present in the copy.

## 11. Reference: the two glpG PDBs disagree about TM4

Two sources are in use and they disagree, which matters when shading helical regions on a figure:
`glpg_oriented.pdb` (crystal-derived, oriented post hoc) and the representative structure from the membrane
REMD simulation. The sequence is identical in the TM4 region. In the crystal PDB residues 135-140 sit at
z = +2 to +9 A, splayed toward the extracellular surface, so their backbone H-bond geometry fails DSSP; in
the membrane-equilibrated structure they sit 6-10 A deeper (z = -0.6 to +2.5 A), properly threaded through
the bilayer, and DSSP assigns them as helix. The simulation PDB agrees with the 2IC8 literature boundaries
(construct offset -66) to within 1-2 residues at every helix, while the crystal PDB gives a broken TM2
(83-96 against 82-103) and a truncated TM4 (141-149 against 135-151). **Use the simulation PDB's DSSP for
figure annotation.** Simulation-PDB helices: TM1 29-48, TM2 82-103, TM3 105-127, TM4 135-151, TM5 161-175,
TM6 185-207.

Construct mapping, verified against UniProt P09391 and the experimental HDX file: construct index + 66 =
E. coli GlpG numbering; the base construct is already the catalytically dead S201T (construct 135),
`79HIS`/`79ALA` is WT H145 vs H145A, `S115T` is S181T, and RKRK is the C-terminal tag at 207-210. No proline
is in the donor list (Upside excludes all six) and there are no chain breaks or resseq gaps.
---

## 12. Claims that turned out to be wrong

One line each: what was believed, what is true, and why it is worth keeping.

* **findings 87 (withdrawn by findings 88, cited elsewhere as findings-88):** the glpG blow-ups were a
  timestep failure at protein-lipid contacts. Wrong: `omega*dt` had been computed for contacts that were all O sites, whose force the engine
  discarded entirely (`propagate_deriv` marked O derived and redistributed only BB's gradient), so the
  stiffness measured belonged to pairs exerting no force at all. More force was thrown away than delivered
  (ratio 2.12), leaving 3590 E_up/A of net one-sided force per step acting ON THE ENVIRONMENT, matching the
  blow-up's first symptom. When the stiffest interaction in a diagnosis is also the most suspicious one,
  check that it is connected to the dynamics before building a theory on its magnitude.
* **findings 88's own prediction:** that fixing the discarded O force would collapse the +2-3%
  `avg_kinetic_energy/1.5kT` excess to 1.000. Refuted by measurement (the excess stayed); the real cause was
  finite-dt error (section 4.2).
* **The engine's `--potential-deriv-agreement` as a correctness gate:** not evidence of anything at this
  system size. The metric is `sqrt(sum (fd-analytic)^2 / sum fd^2)` with an FD step of 1e-3 A against a
  float-precision total potential of order 1e4 E_up, so it is round-off dominated; pre-fix and post-fix
  binaries gave 0.41960 and 0.41959 on the same file. It is a developer probe and its own help says so.
* **findings 90 (corrected by findings 92):** the table's force-free core was ~500 kT out of reach. Wrong,
  because the margin was computed for an inertial particle while ION and LIPID are overdamped Brownian.
  Check which integrator governs the particles before computing a stability margin for them.
* **findings 92/93 as the trigger:** the zero-force core was named as what starts a blow-up. It is where the
  cascade ends up; the trigger is two environment beads reaching 1.83 A, inside the valid table domain.
* **findings 95 (retracted by findings 118):** K525 and K541 supported as reduced-labelling sites. That
  support was the defect, produced by one orientation pressing a C-terminal run onto the surface while the
  protein unravelled on an uncorrected table and CB placement.
* **findings 96 (retracted by findings 104):** the `mean_pf >= 0.99999` clip called a bug and removed from
  the T-slice. It is master's deliberate convention and it encodes the estimator's statistical limit; the
  values it hid rest on fewer than two effective frames. Before calling a threshold in inherited analysis
  code a bug, diff it against the reference implementation and ask what statistical limit it might be
  encoding.
* **findings 104's follow-on (corrected by findings 113/114):** that the unclipped rendering caused the
  discrete spikes. Causality backwards: the spikes were the missing lipid-shielding term, and with it in
  place the unclipped profile rises continuously. What stands from 104/109 is the resolution caveat, not the
  clip.
* **findings 106 item 4 (corrected by findings 113, wrong by 7x):** the missing membrane term measured as
  "6 amides in the hydrophobic core", concluding that fixing it would not improve the figure. Both the
  >50%-exchanging cutoff and the |z - midplane| < 10 A criterion were wrong for the purpose, since dG is
  logarithmic and tail contact rather than midplane distance decides where water is; the real size is 44
  extra `+inf` donors and it *is* the explanation for the figure.
* **findings 107 (narrowed by findings 108):** "the fold defect visible in HDX is 3 amides out of 148" is
  true of the *protection state* and badly understates the fold problem, because PS is `H-bond OR burial`
  and burial does almost all the work in a membrane protein.
* **findings 121's first inference:** the cluster NaN were an output artifact, since they recur on the
  exchange period with clean neighbours. Every observation was real and the inference was not; the driver's
  own log says the replicas genuinely blew up and were rolled back.
* **The implicit-vs-hybrid comparison (originally recorded as a second Update 124):** both axes used the
  protein-only protection state, on the reasoning that holding the analysis fixed makes the axes compare
  models. Wrong for this pair of models, and it inverted the sign of the result (section 5.3). Retracted
  with it: the claim that the hybrid's most protected amides "read low", its attribution to the loose-helix
  defect, and the argument against a sampling explanation.
* **The zero-variance cluster HDX read as "needs more blocks":** attributed to 48 short descendants of one
  seed. The real cause was that the protein was rigid (findings 116).
* **The stage-7 freeze:** accepted on canonical kinetic energy and retained secondary structure, both of
  which a high-friction g-JF process satisfies while its coordinates barely move.
* **"The group quota leaves 195 GB, so the pre-ff3 ladder has to be deleted" (2026-09-08, corrected
  2026-09-09):** the glpG data is on `/project`, which has **1514 GB** free; 195 GB is the
  `/project2` group quota, a different filesystem. `rcchelp quota` reports **four** separate
  `trsosnic` group block quotas (`/beagle3`, `/project`, `/project2`, `/cds3`), so the section
  header (`mounted at <path>`) has to be matched against the mount the data is on. Taking the first
  `trsosnic blocks` row reads `/beagle3` on midway3 and `/project2` on midway2, neither of which is
  the right filesystem. `df -h` on the data path gave the correct 1.5 T and was dismissed as
  "the whole filesystem, not the group quota". On `/project` the fileset *is* the group's 3.9 T
  allocation, so `df` and the quota agree there (1514 G vs 1515 G), and that agreement is the
  cross-check to run. Consequence of the error: an argument for deleting 89 GB of completed
  baseline trajectories that did not need deleting. `check_quota.py` now matches the section
  header, takes `min(group quota, statvfs)`, and stamps the value to `QUOTA_HEADROOM_GB` because
  `rcchelp quota` only answers fully on midway3 (on midway2 `/project` is a remote fileset and every
  quota interface fails partway). A stamp older than 24 h is refused rather than used.
* **The BB-env PMF:** built to fix a protein "kick" that was a setup artifact (a non-standard timestep
  inherited from the abandoned CGL plus under-resolved lipids driving a displacement cap), and the PMF then
  caused the drift it was meant to prevent. Rule out setup artifacts, timestep and sub-step resolution above
  all, before building a corrective force-field term.
* **Rewriting `plot_ref_style.py` to draw censored amides as bounds (2026-09-12, reverted same day):**
  asked to make a sparse-looking dG figure more informative, I replaced the off-scale excursions with
  hollow carets on each temperature's resolution limit, broke the profile line across every censored
  amide, and retightened the axis from `(-20,30)` to `(-4,8.6)`. Both ideas were wrong. Breaking the line
  fragments a profile that is read as one continuous curve per temperature, and putting every censored
  amide at `dg_limit` asserts that unmeasurably-different values are all equal to ~6 kcal/mol while
  capping the visible range -- a worse distortion than the excursion it replaced. The excursion rendering
  was a **deliberate choice already argued in the file's own docstring** ("reads as one continuous
  excursion rather than a capped plateau ... that is how these profiles are conventionally read"), and I
  overrode it and presented the result as an improvement. Only the `--temperatures` default (adding
  T=0.90) survived.
  **Rules taken from it.** When a file documents *why* it does something, that rationale outranks my
  judgement about how the output should look; change it only if the user asks or the rationale is
  demonstrably false, and say which it is. Never collapse right-censored values onto one ceiling value,
  and never introduce gaps into a curve read as continuous. And when the complaint is "this looks like it
  lacks data", fix what is plotted (here: which rungs are drawn) before restyling how it is drawn.
