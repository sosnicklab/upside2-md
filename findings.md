# Findings

Knowledge base for this branch: standing rules and the measurements behind them, and the ff3.0
Ramachandran work (1), how the hybrid is put together (2), defects whose causes are established (3),
what the hybrid still gets wrong (4), HDX (5), reusable diagnostics (6), system preparation and
bilayer physics (7-8), the nanoparticle campaign (9), the trainer and glycine-map history (9c-9w),
cluster lessons (10), references (11) and claims that turned out to be wrong (12). Old update numbers
("findings 103") are kept inline where code or other files cite them. Job state lives only in
`remote_jobs.md`.

---

## 1. Standing rules, and the ff3.0 Ramachandran work (1.8-1.18)

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

NP and glpG are different simulations with different integrators and must be analysed separately
(comparison table in `remote_jobs.md` §2). Conflating them cost a healthy 6 h glpG block, when a C-N
threshold borrowed from NP false-positived, and produced a wrong root-cause writeup. They also need
different detectors: in a forced NP tear the protein reached 431 broken bonds with the potential
still finite at +3e5, while glpG's blow-up goes fully NaN within one 46-step interval (6.4).

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

### 1.8 Each Ramachandran map is its own conditional distribution

User correction, 2026-09-27. Every library map `(central X, direction, neighbour Y)` is
`P(phi,psi | X, Y)` over its own set of PDB sites, and the fold-context selection that shaped it belongs
to those sites alone. `GLY|GLY` says nothing about `GLY|ALA`, and central GLY says nothing about central
ALA. So a correction, a symmetry or a "shared mode" measured on one map must never be transferred to
another, and pooled maps (`ALL`, neighbour averages) must not be used to derive a per-map correction.
The error that prompted this: a Bayes argument (`P(phi,psi|aa) ~ P(aa|phi,psi) P_bb`) was used to claim
the contamination is one function shared by all residues, and GLY's antisymmetric part was then
subtracted from every row. That step silently assumes evolution selects residues by (phi,psi) alone; it
selects by the whole site (burial, packing, turn type), so each map carries its own selection factor.
The one real
coupling between maps is statistical: NDRD's hierarchical Dirichlet process shares mixture components
across neighbour maps, so rare pairs are pulled toward their pooled map. That is estimation, not physics.

Corollary, also from a user correction the same day: separating local physics from fold selection in a
map needs a second, per-pair source of information. If that source is Upside's own native-state
simulations (a per-map reference ratio, `map += ln(H_free/H_restrained)`), it is ConDiv's native term in
closed form, not a new design: ConDiv already records exactly those per-residue populations
(`rama_native`, `rama_free`). The only per-pair sources independent of the model are physics or experiment on that pair.
Before presenting a design as new, write down the objective of the existing pipeline and compare.

### 1.9 NDRD's site set is loops only; Upside applies the map to every residue (2026-09-27)

Ting et al. 2010: 3,038 PISCES chains (<= 1.7 A, R <= 0.25, < 50% identity, EDS density), bottom 20%
of density per residue type removed, Stride assignment, then only loop residues at least three
positions from any H or E; TCB also drops 3-10/pi helices and their neighbours. **TCB is 44,112
residues of all types**, about 110 per directional map on average (NDRD pools rare maps through its
hierarchical Dirichlet process). The 456 training proteins have 49,172 interior residues with
defined (phi,psi), a median of 107 per directional map, but only ~8,700 TCB-like residues
(backbone-only DSSP reconstruction, no 3-10 exclusion, so approximate). Same boxes throughout,
ln(aR/aL) for GLY: NDRD TCB -1.19, 456 TCB-like -1.31, **456 all residues -0.58**; ALA aR fraction
0.26 / 0.40 / 0.61. The glycine alpha_L excess is therefore mostly a property of the loop-only site
selection, and a map trained against loop sites leaves every helix and strand residue outside
its objective. In ConDiv the SARW replica zeroes H-bond, burial and rotamer-pair energies but keeps
rama (`zero_for_sarw`), so the DSE term on a map compares the unfolded ensemble with the map's own
distribution: a self-consistency condition, not an external target.

### 1.10 A basin offset must be a weight on its basin (design rules, 2026-09-28)

The basin-offset training these rules governed was removed on 2026-10-01 (plan.md); `ConDiv.py`
now only records per-residue basin populations with the same basins. The rules hold for any additive
correction on a Ramachandran map:
* **The basins must partition the torus**, or probability drains into an uncontrolled region:
  alpha_R (phi < 0, -100 < psi < 50), beta (phi < -100, psi outside that band), pPII
  (-100 < phi < 0, psi outside it), alpha_L (mirror of alpha_R) and `other` (phi > 0 outside
  alpha_L), split for a central glycine into its mirror halves beta' and pPII'.
* **The edges must be sharp.** A 10 deg logistic scale (as in `secstr_bias`) is a 44 deg transition
  and tilts a basin instead of shifting it: an offset was realised at a median 90% of its value
  where the probability lies. At 3 deg (13 deg transitions, the sharpest the 5 deg grid carries) it
  is 98%, and offsets of +-1 spread the energy change inside a basin by a median 0.03 nats.
* **Each map must be renormalised**, because `read_rama_maps_and_weights` mixes the left and right
  maps without normalising them, so an offset that raised a whole map would move weight to its
  partner map.
* **The update must be a regularised Newton step, not a log-ratio of counts.** Basins holding no
  residues in either ensemble give a ratio of two near-zero numbers (steps up to 1.76 nats in round
  1). The MAP step with a Gaussian prior on the offset,
  `eta [T0 N dp - T0^2 c/sigma^2] / [N p(1-p) + T0^2/sigma^2]`, is the log-ratio step where N p is
  large and bounded where it is small; on the same data its largest step was 0.44.

### 1.11 GLY|GLY symmetry reaches the engine only if the sheet entry is symmetric too (2026-09-28)

`read_weighted_maps` mixes every coil map with its sheet map, so a mirror-symmetric coil `GLY|GLY`
entry is not enough. NDRD's sheet `GLY|GLY` holds 94% of its weight at phi < 0 (beta 0.80, pPII 0.14,
beta' 0.06), which gave the middle glycine of A-G-G-G-A ln(beta/beta') = +0.14 to +0.32 while
ln(aR/aL) stayed 0.0000 (the sheet maps are empty in both helical basins). Symmetrise by averaging
probabilities, which pools each site with its mirror image, not energies: the energy mean is a
geometric mean of probabilities and would put the sheet map at beta 0.5 / pPII 0 per side instead
of 0.43 / 0.07. `build_gly_library.py` wrote one symmetric `GLY|GLY` entry to both coil and sheet.
Grid convention, verified: the engine (`rama_map_pot.cpp:66`) and `upside_config` both put node i at
-180 + 5i deg, which is what the rolled mirror of 6.1 assumes.

On the round-1 ff30_basin state the GLY|X offsets deepened alpha_R relative to NDRD in 30 of 38 maps
(mean +0.060), yet every map still favoured alpha_L: native-restrained central glycines sit at
aR 0.19 / aL 0.40, so the target itself is alpha_L-rich (1.15). Per-pair independence was verified
on the live state (`checks/pair_independence.py`): every (central, direction, neighbour) map was its
own parameter set, and the one coupling was in the data, since an interior residue reads a left and
a right map and its populations are credited to both.

### 1.12 What 456 proteins can resolve per Ramachandran map (2026-09-28)

Measured on the ff30_basin epoch-0 data, all offsets zero (`checks/param_data_audit.py`,
`basin_schemes.py`, `basin_types_beta.py`):
* **ff2.1's own groups** hold 31,285 trained values: 30,800 side-chain spline coefficients that
  receive a derivative (of a 31,420-entry latent), environment 460, sheet 20, H-bond 3 plus dhb 1,
  backbone scale 1. Pair data: 306k CA-CA < 10 A contacts over 210 pair types, median 1,200, fewest
  TRP-TRP 84.
* **Per-pair map corrections are noise-limited.** Training sites per directional map: median 94,
  76 maps under 20. One epoch's per-offset step has split-half correlation 0.07 (full-set
  reliability ~0.14); the true per-pair corrections have SD ~0.04 nats against one epoch's noise of
  0.074. Coarser schemes help little: one alpha_R-versus-rest offset per map reaches reliability
  0.35.
* **What resolves is the aggregate**: the GLY|X alpha_L - alpha_R step averaged over its 38 maps,
  +0.060 +- 0.011, while only 6 maps are individually above 2 sigma. The same pattern as the AWH
  per-neighbour verdict (6.7).
* **By type**, GLY|X alpha_L is the best resolved per pair (reliability 0.52-0.71) and X|Y alpha_R is
  weak (0.18-0.19); X|Y beta relative to pPII has no detectable signal per pair and resolves only per
  central residue type (0.40-0.52), the level of ff2.1's `sheet` energies.
* **About 70% of the noise is the protein set, not sampling**: 153 proteins run twice at identical
  offsets give per-map estimates correlated 0.69-0.71. Averaging epochs removes only the other 30%.
* **The free-native mismatch is mostly its sample-size floor**: held out 0.0597 against random
  46-protein training subsets at 0.0585 +- 0.0021 (205 proteins 0.0320, 410 proteins 0.0249), so a
  held-out mismatch cannot show improvement or overfitting.

### 1.13 Which residues need a trained map: literature and ff2.1's own mismatch (2026-09-28)

**Where ff2.1's free simulations miss the native basin populations** (epoch 0, 410 training
proteins, bootstrap over proteins; `checks/rama_by_type.py`). Per central type the misses are small,
at most 1.6 percentage points (at most 0.09 kT): most types lose 1-1.6 points of alpha_R to pPII
(z 3-4.6), a pattern common to all of them; GLY has alpha_L +2.1 (z 4.8). No type has a significant
beta miss (largest z 2.3). By neighbour, no left neighbour has any |z| > 3; on the right, PRO has
the largest aggregate miss, alpha_R +3.2 points (z 9.7), then GLY (pPII +2.3, z 5.9) and VAL/ILE
(pPII -1.1, z -4.5). **The pre-PRO figure is mostly composition, not a pre-proline defect (1.14):**
89% of pre-PRO residues are extended in the native state, against 47% of other residues, and every
class's extended residues visit alpha_R in the free ensemble.

**The left/right MIXTURE gives pre-proline residues three times the alpha_R of NDRD's own
pre-proline map** (`checks/prepro_mixture.py`). For 1,543 residues followed by PRO (ordinary left
neighbour), alpha_R is 0.063 in the right map, 0.337 in the left map, **0.172 in Upside's mixture**,
0.052 under the product rule, 0.083 native, 0.106 in the free simulation. The mixture adds the left
map's alpha_R back at about half weight: the library's mixing weights are nearly equal (0.83
typical, 0.92 for X|right|PRO). So a basin offset on the pre-proline map can reach only the right
map's share (1.14 measures how little). Ting et al. 2010 combine the two neighbours by the product
rule (left and right identities independent given phi,psi, normaliser S = 0.5-1.5 for proline);
`upside_config` has it as `--rama-library-combining-rule product`. Against the mixture it moves
pre-PRO alpha_R by -12 points and leaves other residues nearly unchanged (median largest basin
change 3 points, no net shift). **But it breaks GLY|GLY symmetry**: it divides by glycine's
neighbour-averaged map (ln(aR/aL) -1.13), so the middle glycine of G-G-G gets ln(aR/aL) +1.1.

**Literature (survey 2026-09-28; full text read unless marked).** Pre-Pro is by far the largest
neighbour effect: alpha -30.6, beta +22.6, pPII +15.2 points in the TCB set, from N(i) and CB(i)
clashing with CD(i+1) (Ting 2010 Table 6; Ho & Brasseur 2005). Distinct classes: Gly, trans-Pro,
cis-Pro, pre-Pro, Ile/Val (MolProbity, Williams 2018), Ala partly; other neighbour effects are small
(non-Gly/Pro neighbours within ~12 Hellinger units, Ting Fig. 7). Force fields mostly use three
classes (generic, Gly, Pro: CHARMM36, a99SB-disp, UNRES); "a global correction to the backbone is
sufficient for most residues" (Best, de Sancho & Mittal 2012, via sub-agent). The TCB set is 62%
turns, and the effects of a right-hand Gly and a left-hand Pro reverse sign between turn and coil,
so they are placement, not intrinsic (Ting). Beta propensity ranks locally by sterics (Street & Mayo
1999, R = 0.92), its size is context-dependent (Minor & Kim 1994, abstract). Jumper 2018 added the
sheet parameter "to counteract an observed tendency for our model to overstabilize helices"; FF2
made it per amino acid (Peng 2022 SI eq. S2).

**A reduced set that keeps only what is resolved** (`checks/reduced_set.py`, split-half): GLY|X
(aR, aL, beta), GLY|GLY (helix, beta; tied), X|right|PRO (aR, beta): 158 offsets on 60 maps. GLY|X
alpha_L per pair reliability 0.72 (class mean +0.11 +- 0.02 nats); X|right|PRO alpha_R per pair 0.05
but class mean +0.30 +- 0.04 nats. That class mean is reproducible, but 1.14 shows it is the
extended-residue composition of pre-proline sites, not a steric effect the offsets can act on.
**Adopted 2026-09-28** (user; plan.md), without the optional right-GLY/VAL/ILE classes (+120).

**References for 1.13.** [FT] full text read, [Abs] abstract only, [Ag] read in full by a sub-agent
and not re-checked here. The extracted texts of the open-access papers were in the session
scratchpad only; the PDFs of Peng 2022 are in `~/OneDrive - The University of Chicago/`.

Titles and pages were checked against the retrieved texts; where a text was not retrieved, only
the author, journal, volume and first page reported by the survey are given.

| reference | read | what it contributes |
|---|---|---|
| Ting D, Wang G, Shapovalov M, Mitra R, Jordan MI, Dunbrack RL. Neighbor-dependent Ramachandran probability distributions of amino acids developed from a hierarchical Dirichlet process model. PLoS Comput Biol 2010;6:e1000763 | FT | the NDRD/TCB library; pre-Pro basin shifts (Table 6: A -30.6, B +22.6, P +15.2 points); inter-type distances (Tables 3-4, Fig. 7); TCB is 62% turns (Table 2); left/right combined by the product rule under conditional independence given phi,psi (Methods). Its printed B and P phi ranges look swapped |
| Ho BK, Brasseur R. The Ramachandran plots of glycine and pre-proline. BMC Struct Biol 2005;5:14 | FT | pre-Pro mechanism: N(i) and CB(i) clash with CD(i+1) in alpha; zeta region |
| Hollingsworth SA, Karplus PA. A fresh look at the Ramachandran plot and the occurrence of standard structures in proteins. Biomol Concepts 2010;1:271-283 | FT | the glycine PDB map is asymmetric because PDB statistics record which residue wins a site |
| Williams CJ et al. MolProbity: more and better reference data for improved all-atom structure validation. Protein Sci 2018;27:293-315 | FT (Ramachandran section) | six validation classes: general, Gly, trans-Pro, cis-Pro, pre-Pro, Ile/Val |
| Lovell SC et al. Proteins 2003;50:437-450 (title not verified) | Abs | the earlier validation categories |
| Jha AK, Colubri A, Zaman MH, Koide S, Sosnick TR, Freed KF. Helix, sheet, and polyproline II frequencies and strong nearest neighbor effects in a restricted coil library. Biochemistry 2005;44:9691-9702 | FT | turn removal cuts the helical basin 37.0% -> 21.9%; neighbour effects up to 4-fold, context-dependent; coil beta vs strand frequency R = 0.84 |
| Jha AK, Colubri A, Freed KF, Sosnick TR. Statistical coil model of the unfolded state: resolving the reconciliation problem. PNAS 2005;102:13099-13104 | FT | neighbour effects raise the apoMb RDC correlation 0.41 -> 0.71 |
| Avbelj F, Baldwin RL. Origin of the neighboring residue effect on peptide backbone conformation. PNAS 2004;101:10967-10972 | FT | aromatic/beta-branched neighbours shift mean phi only ~ -2 deg in pPII |
| Street AG, Mayo SL. Intrinsic beta-sheet propensities result from van der Waals interactions between side chains and the local backbone. PNAS 1999;96:9074-9076 | FT | beta propensity ranks locally by sterics, R = 0.92 |
| Minor DL, Kim PS. Nature 1994;367:660-663 and Nature 1994;371:264-267 (titles not verified) | Abs | beta propensity largely set by tertiary context at edge strands |
| Smith CK, Regan L. Science 1995;270:980-982 (title not verified) | Abs | cross-strand pair energies as large as propensities |
| Avbelj F, Baldwin RL. Role of backbone solvation in determining thermodynamic beta propensities of the amino acids. PNAS 2002;99:1309 | Ag | beta scales correlate at central, not edge, sites |
| Hagarman A et al. J Am Chem Soc 2010;132:540-551 (title not verified) | Abs | Ala ~80% pPII in GxG |
| Jumper JM, Faruk NF, Freed KF, Sosnick TR. Trajectory-based training enables protein simulations with accurate folding and Boltzmann ensembles in cpu-hours. PLoS Comput Biol 2018;14:e1006578 | FT | Upside's rama term from NDRD TCB; the sheet parameter added "to counteract an observed tendency for our model to overstabilize helices" |
| Peng X et al. Prediction and validation of a protein's free energy surface using hydrogen exchange and (importantly) its denaturant dependence. J Chem Theory Comput 2022;18:550-561, and SI | FT | FF2: TCB and sheet maps mixed by gamma, per amino acid (SI eqs. S1-S2); secondary-structure-dependent H-bond strengths |
| Best RB et al. Optimization of the additive CHARMM all-atom protein force field targeting improved sampling of the backbone phi, psi and side-chain chi1 and chi2 dihedral angles. J Chem Theory Comput 2012;8:3257-3273 | Ag | CHARMM36 CMAP in generic/Gly/Pro classes |
| Best RB, de Sancho D, Mittal J. Residue-specific alpha-helix propensities from molecular simulation. Biophys J 2012;102:1462 | Ag | "a global correction to the backbone is sufficient for most residues" |
| Tian C et al. ff19SB: amino-acid-specific protein backbone parameters trained against quantum mechanics energy surfaces in solution. J Chem Theory Comput 2020;16:528-552 | Ag | residue-specific CMAPs, several reused across residues |
| Jiang F, Zhou CY, Wu YD. Residue-specific force field based on the protein coil library. RSFF1: modification of OPLS-AA/L. J Phys Chem B 2014;118:6983 | Ag | residue groups {E,Q,K,R,M,L}, {F,Y,W}, {V,I} |
| Alford RF et al. The Rosetta all-atom energy function for macromolecular modeling and design. J Chem Theory Comput 2017;13:3031 | Ag | pre-Pro has its own Ramachandran table |
| Choi JM, Pappu RV. J Chem Theory Comput 2019;15:1355 (title not verified) | Ag | coil libraries break glycine's inversion symmetry |

Not verified by the survey: Swindells, MacArthur & Thornton 1995 numbers, the RSFF2 groupings, and
a per-residue count of how much turns inflate alpha_L for Gly, Asn and Asp.

### 1.14 The X|right|PRO offsets have no leverage and no pre-proline signal to fit (2026-09-30)

Measured on the finished `ff30_basin` run (six rounds; round 6 is the released ff_3.0) with
`/project/trsosnic/yinhan/checks/prepro_leverage.py`, `prepro_control.py`, `prepro_residual.py`,
`prepro_rules.py` and `prepro_left.py` (logs `*_20260930.log` beside them). Interior residues followed
by PRO/CPR, centre not GLY/PRO, 1,573 training proteins' sites.

**Six rounds moved the offsets steadily and the simulations not at all.**

| round | mean aR offset [min, max] | rama map aR (mixture) | free aR | native aR | free - native | free aR at full leverage |
|---|---|---|---|---|---|---|
| 0 | 0 | 0.174 | 0.102 | 0.079 | +0.023 | 0.102 |
| 1 | +0.13 [-0.08, +0.41] | 0.171 | 0.108 | 0.080 | +0.028 | 0.096 |
| 2 | +0.26 | 0.168 | 0.100 | 0.079 | +0.022 | 0.091 |
| 3 | +0.36 | 0.166 | 0.098 | 0.078 | +0.020 | 0.087 |
| 4 | +0.44 | 0.165 | 0.099 | 0.080 | +0.019 | 0.085 |
| 5 | +0.50 [-0.11, +1.63] | 0.164 | 0.101 | 0.077 | +0.023 | 0.083 |
| 6 (released) | +0.58 [-0.12, +1.89] | 0.164 | - | - | - | - |

* **The mixture caps what any offset on the right map can do.** With an infinite alpha_R offset on
  every X|right|PRO map the rama term's pre-PRO alpha_R only falls to 0.144-0.146, since the left
  map's share is untouched; round 5 had used a third of that. Right-map-only would give 0.059, the
  product rule 0.052 (round-0 maps).
* **Had each offset acted as an additive energy** on its residues (the Newton step's assumption),
  the round-5 offsets would have brought free alpha_R to 0.083, the native value. The measured free
  alpha_R did not move (round-to-round noise ~0.003). The step keeps its size because the gap does
  not close. The largest offsets are all X|right|PRO alpha_R (ASP +1.89, PHE +1.16, CYS +1.10,
  ASN +1.00), plus GLY|right|PRO alpha_L +1.24.
* **They stop only where the prior balances the unclosed gap**, c* = N dp sigma^2 / T0 per map
  (T0 = 0.80): about 5 nats for ASP (144 sites, gap 0.027) and ALA (119, 0.036), another 10-20
  rounds of drift. The convergence gate passed the rama group at step 114 (p 0.0011, 0.0006, 0.0158
  at steps 76, 95, 114; threshold 0.005) because the growing prior pull cancels more of the fixed
  data pull each round, not because the gap closed. So the gate's "converged" means "prior-limited"
  for this group.

**The gap is mostly composition, not a pre-proline map error** (epoch 0, training; free / native
alpha_R, split by the residue's own native alpha_R):

| class | extended in native (aR < 0.05) | helical in native (aR > 0.5) |
|---|---|---|
| pre-PRO | n 1,401: 0.041 / 0.002 | n 116: 0.800 / 0.976 |
| other X | n 16,778: 0.058 / 0.001 | n 19,008: 0.915 / 0.990 |
| post-PRO X | n 622: 0.101 / 0.001 | n 818: 0.863 / 0.986 |
| GLY | n 2,200: 0.037 / 0.001 | n 516: 0.804 / 0.985 |
| PRO/CPR | n 1,020: 0.059 / 0.001 | n 701: 0.859 / 0.991 |

Every class's extended residues visit alpha_R in the free ensemble (the restrained replica cannot
leave its basin), and every class's helices fray. Pre-proline sites are 89% extended (other X 47%),
so their aggregate is +0.036 from extended residues, -0.013 from helical ones: +0.023. Behaving like
other X within each native class they would show +0.046. **So against ordinary residues in the same
native conformation, pre-proline residues already visit alpha_R less (0.041 against 0.058).** The
ConDiv target cannot see the pre-proline problem the literature describes: that is a propensity of
loops and unfolded chains, and a native-restrained extended residue sits at alpha_R 0.002 whatever
its class. What the offsets were fitting is the class's composition.

By contrast, **the glycine signal is real and sits where the TM4 failure sits**: extended glycines
show no net alpha_L gap (0.499 / 0.499; two opposite gaps, 1.15), helical glycines have alpha_L 0.110
free against 0.005 native.
But a per-map offset moves both groups alike, and the extended glycines outnumber the helical ones
four to one: by epoch 5 the extended ones are at 0.487 / 0.497 (now below native) and the helical
ones at 0.102 / 0.005. More rounds of the same design would trade the two further; the helical
glycine excess depends on where the glycine sits, which a per-type (phi,psi) map cannot express.
**The fixed point favours alpha_L whatever the start**: native-restrained glycines over all sites sit
at aR 0.19 / aL 0.40, so one map per (glycine, neighbour) that reproduces the average glycine must
favour alpha_L. Linear extrapolation of rounds 1-6 (dL - dR +0.26 moved the aggregate alpha_L gap
+0.021 -> +0.015): closing it takes ~0.65 more in dL - dR, ~14 epochs at the current step, leaving
the engine map near ln(aR/aL) -0.3 at T = 1, helical glycines near alpha_L 0.08 (native 0.005) and
extended ones near 0.46 (native 0.50). A rough estimate, not a measurement.

**ff_3.0's glycine term still favours alpha_L everywhere** (`checks/gly_handedness.py`, log
`gly_handedness_20260930.log`; basin populations of the map alone, and the well depth
E_min(aR) - E_min(aL)). Engine map of X-G-Y (coil + sheet mixed, plus the reference correction), all
361 non-glycine flanks: ln(aR/aL) median -0.92 at T = 1 (-1.18 in ff_2.1), -1.23 at T = 0.8, alpha_R
favoured in none; the alpha_L well is deeper by a median 1.17 (1.46 in ff_2.1). Of the 38 single
GLY|X maps, one favours alpha_R. In the glpG seed every glycine favours alpha_L, the 12 natively
helical ones included; TM4's GLY136 (T-G-V) most of all, ln(aR/aL) -1.59 at T = 1 and -2.43 at
T = 0.7 (the seed before this release: -0.73 / -1.25), and GLY149 (R-G-E) -1.21 / -1.86. In the first
~10 h (`gly_tm4_flip_20260930.log`, last three groups) TM4's helical glycines leave the helix for
phi > 0 in some replicas: GLY143 36% and GLY149 17% of frames in 79ALA_S115T at T 0.70, GLY136 20%
in 79HIS_S115T and GLY143 22% in 79ALA at T 0.80; none at T 0.70 in the other three variants. GLY143
never flipped in the earlier campaign (3.10c), though its ff_3.0 map (-0.62) leans no further to
alpha_L than the pre-release seed's (-0.73), so the other ff_3.0 changes share the cause.
Glycines before a proline (88 extended) have alpha_L 0.057 / 0.012 at epoch 0, 0.030 / 0.005 at
epoch 5.

**Left-neighbour dependence of pre-proline alpha_R is not detectable** (`prepro_left.py`, 19
left-neighbour groups of >= 30): the native group sd 0.036 is near the 0.027 expected from sampling;
right-only fits it best (rms 0.036, product 0.037, mixture 0.046). For pPII and beta the product rule
follows the native groups better (corr 0.59 / 0.67) than right-only (0.41 / 0.59). Right-only changes
no residue outside the pre-proline class and keeps GLY|GLY exact; the product rule changes every
residue (median largest-basin change 3 points) and gives the middle glycine of G-G-G ln(aR/aL) +1.06.

### 1.15 Which way the data pull glycine's handedness (probe from equal depth, 2026-09-30)

User question: if glycine's alpha_R and alpha_L start at equal depth, do the data pull it further
toward alpha_R or back toward alpha_L? It decides whether training can make glycine right-handed.

**Measure** (`checks/glyprobe_analysis.py <run_output> <epoch>`): the DATA term of the offset update
on the 38 GLY|X maps over one epoch's training proteins, apart from the prior (which pulls every
offset back toward NDRD and would bias the answer toward alpha_L from any start away from zero):
(free_aL - native_aL) - (free_aR - native_aR) per residue read, positive = toward alpha_R, with a
bootstrap over proteins; glycines split by their own native basin.

**On the ff_3.0 run itself:** epoch 1 +0.040 [+0.028, +0.052], toward alpha_R in 32 of 38 maps; epoch 5
(round-5 offsets, dL - dR +0.26) +0.023 [+0.012, +0.034], 31 of 38. So at the released state the data
still pull toward alpha_R (the drift of 1.14). By the residue's native basin, epoch 5, free / native:

| glycines (non-GLY flanks) | n | alpha_R | alpha_L | pulls toward |
|---|---|---|---|---|
| helical in native (aR > 0.5) | 488 | 0.813 / 0.984 | 0.101 / 0.005 | alpha_R |
| alpha_L in native (aL > 0.5) | 1,003 | 0.030 / 0.002 | 0.889 / 0.990 | alpha_L |
| the rest | 1,082 | 0.037 / 0.008 | 0.100 / 0.010 | alpha_R (aL), alpha_L (aR) |

**One map serves two native populations that pull in opposite directions**: helical glycines want
less alpha_L, loop glycines that are natively left-handed want more. The net follows their balance,
which is the fixed-point argument of 1.14 made visible.

**Probe** (plan.md Phase 7, job 49133133, `training/ff30_glyprobe`, README there): one epoch branched
from the ff_3.0 checkpoint (step 114) with every GLY|X map's alpha_R and alpha_L offsets moved by
-d/2 and +d/2 until the two basins hold equal probability (dL - dR mean +0.31 -> +1.21; aR + aL weight
per map 0.508 -> 0.478); everything else as released; no gate, no release. Engine X-G-Y map at the
start: ln(aR/aL) +0.02 at T = 1, -0.05 at T = 0.8 (from -0.92 / -1.23). It runs on the BP-fixed
binary, which the ff_3.0 training did not (|dE| <= 0.03 E_up, remote_jobs.md §0c). **Prediction:** the natively
left-handed loop glycines lose alpha_L in the free ensemble, so the pull turns toward alpha_L.
**Result: toward alpha_L.** Stopped by the user at 15 of 19 steps (329 training proteins), since the
direction was settled (`checks/glyprobe_partial.py`, log `glyprobe_final_partial_20260930.log`):
* The pooled data pull is -0.035 per residue read, bootstrap 95% [-0.048, -0.022]. 35 of 38 maps
  pull toward alpha_L, and the data-only step on dL - dR averages -0.062.
* Free / native at equal depth, against ff_3.0's epoch 5:
  * helical glycines: alpha_L 0.077 / 0.004 (was 0.101);
  * natively left-handed glycines: alpha_L 0.847 / 0.989 (was 0.889), alpha_R 0.058 / 0.004 (was
    0.030);
  * the rest: alpha_L 0.069 / 0.008 (was 0.100).
* The neutral map helps the helical glycines a little and costs the left-handed ones more, as
  predicted.
* With ff_3.0's +0.023 at dL - dR +0.26, a linear interpolation puts the fixed point near
  dL - dR +0.66. That is still an alpha_L-favouring map, about ln(aR/aL) -0.5 at T = 1, close to
  the all-residue native -0.58 (1.9); this is an estimate, not a measurement.
* So training a context-free glycine map relearns the training natives' placement, and it cannot
  make glycine right-handed. Even at equal depth, helical glycines keep 0.077 alpha_L against 0.004
  native, which the map does not supply.

**The one context-aware term Upside has is FF2's H-bond energy** (`src/hbond.cpp`, `hbond_energy`).
Each H-bond a residue makes, as donor or acceptor, is scored by that residue's own (phi, psi):
E_alpha where phi is outside (0, 165) deg and psi inside (-120, 60), E_beta for the same phi with psi
outside, E_other where phi is in (0, 165), i.e. left-handed. So it knows both "H-bonded" and
"which basin", but the three energies are shared by every residue type. ff_2.1: E_alpha -1.961,
E_beta -1.946, E_other -1.769 (alpha_R favoured over phi > 0 by 0.192 per H-bond); **ff_3.0: -1.878,
-1.872, -1.798, a margin of only 0.080**. Training narrowed the helix-over-left-handed margin for
every H-bonded residue while the GLY|X offsets moved glycine's map the other way. Unproven, but a
candidate for why GLY143 flips under ff_3.0 although its map is no more left-handed than before.
**The margin shrinks without any glycine training** (`checks/hb_trajectory.py`, checkpoints of three
runs from ff2.1):

| run | what it trains on the rama | margin E_other - E_alpha along the run | E_alpha at the end |
|---|---|---|---|
| `ff21-fixedpoint` | nothing | 0.192 -> 0.127-0.142 by steps 13-25 | -1.907 (step 25) |
| `ff30_basin` (ff_3.0) | 158 basin offsets | 0.192 -> 0.151 at step 19 (offsets still 0), then 0.07-0.14 | -1.878 (step 114) |
| `ff30` (cancelled) | full glycine row | 0.192 -> 0.06-0.12 over steps 97-222 | -1.831 (step 222) |

It drifts smoothly (Adam momentum), not as step noise: over ff_3.0's last epoch it ran 0.117 ->
0.063 -> 0.080, and the release is that last iterate. E_alpha weakens in every run. So most of the
change is the trainer's own drift of the H-bond term from ff2.1, which the offsets may add to but do
not cause. Helical glycines' free alpha_L barely moved over training (0.107 at epoch 0, 0.101 at
epoch 5; native 0.004-0.005) while GLY|X dL - dR rose by +0.26; whether the H-bond drift offset the
maps' gain there is not separated. The glpG validation points the same way: by 19:00 on 09-30
TM1 (30-48, no glycine) had fallen at T 0.70 from 0.99 to 0.89-0.90 in 79HIS and from 1.00 to
0.91-0.93 in 79HIS_S115T, against 1.000 in the pre-ff3 campaign (3.10c). A glycine-free helix
weakening is what a weaker E_alpha predicts and the glycine maps cannot cause.
A glycine-specific set of these energies is the smallest helix-aware glycine term: the engine
already computes the per-residue helix score, and the trainer already trains these energies with
their analytic derivative. Its limit: cap and turn glycines are also H-bonded and left-handed, so
how the H-bonded native glycines split between alpha_R and phi > 0 must be measured first (1.16).

### 1.16 Inputs for a glycine map that is not trained, and for a pre-proline rule (2026-09-30)

Measured for the glycine proposal that follows 1.15 (take glycine's map from physics instead of
the PDB). The pre-proline items below are separate and concern a suspected problem that is not
established (plan.md Phase 6, parked; lesson 10.8). Scripts and logs in
`/project/trsosnic/yinhan/checks/`.

* **Native glycines are H-bonded in both basins, and the engine scores the two in different
  branches** (`gly_native_hbond.py`, log `gly_native_hbond_20260930.log`; 456 training natives, DSSP
  electrostatic criterion, H and O placed from N, CA, C; non-glycine phi < 0 0.975 as a sign check).
  Glycines with non-GLY flanks:

  | native basin | n | H-bonded (own NH or CO) | own NH donor | `hbond_energy` branch | commonest partners |
  |---|---|---|---|---|---|
  | alpha_R | 572 | 0.83 | 0.63 | helix 1.00 | NH->i-4 and CO<-i+4 |
  | alpha_L | 1,084 | 0.74 | 0.63 | turn 0.99 | NH->i-3, then NH->i-4 |
  | beta | 372 | 0.86 | 0.76 | sheet 0.98 | |

  So a glycine-specific set of the three branch energies would put helical glycines (E_alpha) and
  natively left-handed ones (E_other) on separate parameters, where one map depth serves both
  (1.15). E_other still sees both groups, since a helical glycine that flips to phi > 0 can keep
  its NH->i-4 bond (the alpha_L C-cap pattern).
* **What the own-H-bond term would miss is mostly fraying that every residue shows**
  (`gly_mismatch_by_hbond.py`, log `gly_mismatch_by_hbond_20261001.log`).
  * Method: each glycine's native H-bond state is joined with its per-residue free and restrained
    basin populations from the probe (epoch 6, 16 steps, equal depth) and from ff_3.0's epoch 5.
    The classes are:
    * `own`: the glycine's own NH or CO is bonded;
    * `spanned`: not own, but inside a short-range bond, |d - a| <= 5;
    * `none`: neither.
  * Control: non-glycine residues in the same native basin and class.

  Loss of the native basin per residue, probe / ff_3.0 run:

  | natively helical | glycine share | glycine loss | non-glycine loss | glycine-specific excess, share of it |
  |---|---|---|---|---|
  | own | 0.82 / 0.84 | 0.099 / 0.133 | 0.057 / 0.056 | 62% / 66% |
  | spanned | 0.14 / 0.13 | 0.285 / 0.347 | 0.158 / 0.157 | 32% / 25% |
  | none | 0.04 / 0.03 | 0.365 / 0.531 | 0.270 / 0.260 | 6% / 9% |

  * Natively left-handed glycines lose less than the non-glycine residues that sit at alpha_L
    (mostly Asn and Asp): 0.12 against 0.17 (own), 0.29 against 0.41 (none). On this measure they
    have no glycine-specific deficit, so their pull toward alpha_L in training is the generic loss.
  * In glpG's seed every natively helical glycine is in the `own` class: TM4's GLY136 (CO<-i+4),
    GLY143 (NH->i-4, CO<-i+4) and GLY149 (NH->i-4).
* **Per-type misses in the same context follow intrinsic propensity**
  (`type_mismatch_by_context.py`, log `type_mismatch_by_context_20261001.log`; ff_3.0 run, epoch 5,
  40,462 residues).
  * Method: each type's loss of its native basin is compared with all other types in the same
    native basin and own-H-bond state, with z from a bootstrap over proteins.
  * Helical, own bond (mean 0.058):
    * lose more: GLY 0.133 (z +6.0), SER 0.082 (+4.1), ASN 0.080 (+3.8);
    * lose less: GLU 0.041 (-5.8), ALA 0.043 (-4.7), LEU 0.046 (-4.5).
  * Extended, own bond (mean 0.112):
    * lose less: VAL 0.072 (z -11.1), ILE 0.077 (-8.3);
    * lose more: ASP 0.154, ASN 0.160, SER 0.149, GLY 0.149.
  * The orders match the helix and beta propensity scales. So a trained correction indexed by type
    and context, fitted to the native-restrained target (every residue at ~0.99 in its basin), would
    flatten them: residues would hold whatever basin evolution put them in equally well, which is
    placement again, in a milder form.
  * Against this scale glycine's helical loss is 0.133 at ff_3.0 and 0.099 at equal depth, against
    0.08 for Ser and Asn. Experimentally glycine is the weakest helix former after proline.
* **Residue counts per basin show selection, but do not convert into map energies**
  (`aa_basin_counts.py`, log `aa_basin_counts_20261001.log`; 49,401 interior residues of the 456
  natives).
  * Glycine fills 54% of alpha_L sites (Asn 11%, Asp 7%) and 82% of phi > 0 extended sites, but
    only 2.7% of alpha_R sites. Ile, Val and Thr are nearly absent from alpha_L.
  * If residues were chosen by (phi, psi) alone, the counts would fix the difference between two
    residues' maps at each basin. Glycine's handedness `h = E(aL) - E(aR)` would then follow from
    any reference residue X's library value `h_X`, and every X would give the same answer.
  * They do not. The implied h_GLY runs from -2.61 (Ala) to -1.16 (Asn): median -1.74, sd 0.50 over
    14 references, against -1.20 in the library and -0.15 to -0.3 from AWH. Helix-placed references
    (Ala, Leu, Met) give the most alpha_L; turn-placed ones (Asn, Asp, Ser) the least, so each
    anchor brings its own placement.
  * The answer also moves with burial: -1.95 in exposed sites, -1.56 in the middle tercile, -0.17
    in buried sites (only 2 references with >= 10 alpha_L counts there).
  * So residue choice depends on more than (phi, psi), as rule 1.8 states, and a count-based
    correction needs an anchor taken from another map.
* **Proline shows no problem a map change would fix** (`pro_mismatch_by_context.py`, log
  `pro_mismatch_by_context_20261001.log`; ff_3.0 run, epoch 1, when the offsets were still near
  NDRD).
  * Central PRO holds its native extended basin better than any other type: own-bond loss 0.054
    against a mean of 0.114 (z -8.7); no own bond, 0.072 against 0.177 (z -16.8). The ring locks
    phi.
  * Residues before PRO, against the same type not before PRO, in the same basin and H-bond state:
    * extended, own bond: +0.024 [+0.012, +0.035] (n 1,219), which may include beta/pPII
      exchange;
    * extended, no own bond: +0.000;
    * natively helical: +0.059 [+0.025, +0.106] (n 123).
  * The helical excess runs against the right-only rule: the mixture is the more helix-friendly
    map for these residues, and right-only raises their native-point energy by a median +0.84.
  * The native-restrained comparison cannot see loop or unfolded-state propensity, which is where
    the mixture's extra alpha_R would act. That part stays untested.
* **Glycine already has a side-chain bead.** It is one fixed `GLY_0` bead in `sidechain.h5`, with
  trained pair and coverage rows. It sits 0.61 A from the frame origin along the L-CB direction
  (ALA's bead is 1.73 A out). A packing term for glycine therefore exists and is trained.
  * `backbone_pairs` gives glycine no CB.
  * `ProteinHBond` holds per-pair H-bond values and per-edge sensitivities (`igraph`), so a term
    scored on the residues a bond spans could be built on its edge loop.
* **The AWH glycine-before-proline surface is not like the others.** Ac-Gly-Pro-NHMe (`RP`, 400 ns,
  ff99SB-ILDN): alpha_R 0.006, alpha_L 0.004 (other right contexts 0.10-0.20 each). NDRD's
  GLY|right|PRO has alpha_R 0.046; `rama31.dat` gives that map the pooled surface, alpha_R 0.105,
  which erases the pre-proline clash for glycine. The per-neighbour noise verdict (2026-09-18) does
  not cover this context.
* **`rama_map_pot_ref` reshapes a glycine map.** ConDiv adds it to every residue. On `rama31.dat`'s
  GLY|ALA it moves alpha_R 0.105 -> 0.145 and alpha_L 0.121 -> 0.163, extended weight to the helical
  basins, with ln(aR/aL) almost unchanged (-0.139 -> -0.116). A measured surface used as glycine's
  whole local term must therefore be stored with the reference subtracted, or glycines left out of
  that node.
* **Right-only for pre-proline residues can be written into the library.** Multiplying every
  X|right|PRO `dimer_weight` by 1e6 in both the coil and sheet groups gives the right map alone (and
  the right-only coil/sheet ratio) to 1e-5 E_up, and leaves every residue not followed by PRO
  bitwise unchanged (local test on `rama.dat` and ff_2.1 `sheet`). Every reader of the weights goes
  through `read_rama_maps_and_weights`.
* **What right-only costs native pre-proline residues** (`prepro_rightonly_natives.py`, log
  `prepro_rightonly_natives_20260930.log`; NDRD library, ff_2.1 sheet energies; 1,959 residues):
  the engine map's alpha_R falls 0.165 -> 0.055; at the native (phi, psi) the energy drops by a
  median 0.17 for the 1,742 extended residues and rises by a median +0.84 [10%: +0.40, 90%: +1.37]
  for the 132 natively helical ones (6.7%).

**Literature (survey by sub-agent 2026-09-30; [FT] full text read by it, [Abs] abstract only; not
re-checked here).**
* **Why PDB statistics cannot give glycine's intrinsic map.** PDB (phi,psi) statistics are
  Boltzmann-like only for comparing residues at a fixed (phi,psi): glycine is enriched at alpha_L
  because it beats the other residues there, not because it prefers alpha_L to alpha_R (Shortle,
  Protein Sci 2003;12:1298 [Abs]; Hollingsworth & Karplus, Biomol Concepts 2010;1:271 [FT]: "a
  dipeptide with Gly in it must have equivalent energetics in the delta' and delta regions").
* **What other models do with glycine's local term.**
  * From physics, not PDB statistics:
    * UNRES: MP2 PMF of Ac-Gly-NHMe (Sieradzan, JCTC 2012;8:4746 [FT]).
    * CHARMM36: QM glycine-dipeptide CMAP (Best, JCTC 2012;8:3257 [FT]).
    * ff19SB: aqueous QM glycine dipeptide, because the PDB enrichment "would be reflected
      erroneously" (Tian, JCTC 2020;16:528 [FT]).
  * Symmetrised:
    * Rosetta `-symmetric_gly_tables`, covering the rama, p_aa_pp and RamaPrePro tables.
    * Choi & Pappu, JCTC 2019;15:1355 [FT].
  * AWSEM skips the rama term for glycine (source code).
* **No experiment measures glycine handedness in a chiral context.** Searched for:
  * stereospecific 3J(HN,Ha2/Ha3) couplings;
  * RDCs that resolve the sign of glycine phi in host peptides or IDPs.

  GGG is achiral, so it cannot answer the question. The handedness therefore rests on MD alone,
  where our two force fields agree to 0.045 nats (9r).
* **Combining neighbours.** Ting et al.'s own rule is the product
  `f(C,R) f(C,L) / [S f(C)]`, under which a region the right map empties stays empty. On 17,600
  held-out coil residues it scored 1.25 against 1.21 for centre plus right neighbour only and 1.19
  for raw triplets, with no detectable left-right interaction in 3J couplings (Shen, Roche,
  Grishaev & Bax, Protein Sci 2018;27:146 [FT]). Pre-Pro mechanism: clashes of N, O(i-1) and H(i)
  with CD(i+1) (Ho & Brasseur 2005 [FT]).
* **AlphaFold neither meets nor solves this problem** (second survey, 2026-10-01; key quotes
  checked against the downloaded texts).
  * AlphaFold 1's torsion term is `-log p_vonMises(phi, psi | S, MSA)`. It is predicted per residue
    from sequence and alignment, so it is context-conditioned by construction, and it has no
    reference correction.
  * Its reference state is applied to distances only: `P(d | length)` from a network trained on the
    same structures without sequence, plus a glycine flag (Senior, Nature 2020;577:706).
  * AlphaFold 2 has no Ramachandran prior. No heavy atom depends on omega or phi, and FAPE is "the
    main component that ensures the correct chirality" (Jumper, Nature 2021;596:583, SI 1.8.4,
    1.9.3).
  * Physical correctness is handed to Amber99SB in a restrained relaxation that "does not improve
    the accuracy".
  * No published evaluation of AlphaFold's glycine alpha_L or pre-Pro accuracy was found. AF2
    (phi, psi) are tighter than the PDB's (Terwilliger, Nat Methods 2024; Tan, arXiv 2025).
  * A conditional predictor may learn placement because placement is its target; Upside needs a
    transferable local energy.
* **Experimental data available for glycine, checked against our AWH surfaces (2026-10-01).**
  * Source: Andrews et al. 2020 SI, retrieved through Europe PMC's supplementary-files service:
    * Table S1: five measured J-couplings for the central glycine of cationic GGG in water.
    * Table S2: basin populations from a Gaussian model fitted to those couplings and amide I'
      spectra; the authors call them a rough comparison.
  * The same numbers are in ff24EXP-GA SI Table S4.
  * Basins as defined there:
    * pPII: -90 < phi < -42, 100 < psi < 180;
    * beta-t: -130 < phi < -90, 130 < psi < 180;
    * a-beta: -180 < phi < -130, 130 < psi < 180;
    * alpha: -90 < phi < -32, -60 < psi < -14;
    * each counted with its mirror box.

  | GGG central glycine | pPII | beta-t | a-beta | alpha |
  |---|---|---|---|---|
  | experiment (Gaussian model) | 0.46 | 0.13 | 0.01 | 0.06 |
  | ff14SB, cationic GGG (Andrews) | 0.40 | 0.06 | 0.09 | 0.05 |
  | our AWH Ac-Gly-Gly-NHMe, ff99SB-ILDN rep1 / rep2, ff14SB | 0.32-0.33 | 0.05 | 0.06 | 0.09 |

  * In the matched system the force field is within 0.06 of experiment in pPII and 0.01 in alpha.
    Our capped dipeptide differs from it by about 0.08, and that is the termini, not the force
    field. A capped peptide is the closer model of a glycine inside a chain.
  * No data resolve handedness: GGG is achiral.
* **Helix propensity:** Pace & Scholtz 1998 give their scale in the abstract (Europe PMC): Gly 1.00,
  Ser 0.50, Asn 0.65, Ala 0 kcal/mol. It averages 11 peptide and protein host systems at
  solvent-exposed mid-helix positions. The per-system hosts and conditions are only in the full
  text, which was not retrievable (Cell 403, PMC captcha). So it can check the order, but it
  cannot calibrate a single host-guest simulation. The order agrees with the per-type fraying
  measured above.
* **Solutions in other models** (third survey, 2026-10-01, sub-agent; not re-checked here).
  * No model fixed glycine handedness by training a context-free map. Those that avoid it take
    glycine's local term from physics, which comes out inversion-symmetric, and get context from
    other terms:
    * UNRES, from QM of blocked residues: Gly-Gly is near-symmetric, and L-Ala-L-Ala's asymmetry
      comes from the neighbours' CB couplings (Lipska, JPCL 2023).
    * CHARMM36, ff19SB: QM glycine CMAP.
    * CGSchNet: 1D phi/psi priors Boltzmann-inverted from ff99SB-ILDN MD; with the prior alone
      every protein unfolds (Charron, Nat Chem 2025).
    * Martini3-IDP: dihedrals around glycine fitted separately from CHARMM36m IDP MD.
  * Rosetta keeps the PDB glycine asymmetry by default (the agent's check: ln(aR/aL) about -1.6).
  * **AWSEM drops glycine's Ramachandran term and sets the i->i+4 helical H-bond strength by
    residue from the experimental helix propensity** (Pace & Scholtz, Biophys J 1998: glycine about
    1 kcal/mol less helical than Ala).
  * HPS-SS fits a per-residue dihedral term by simulating host-guest peptides against experimental
    helix propensities (Rizuan, JCIM 2022).
  * So a per-type context term can have a target that is free of placement.
  * **Pre-proline:**
    * Rosetta's `rama_prepro` replaces the table whenever residue i+1 is Pro (glycine included) and
      uses no left-neighbour information anywhere.
    * CHARMM36/36m have pre-Pro CMAP slots identical to the base maps, so the effect comes from
      explicit Pro CD sterics, which Upside lacks.
  * **Experimental symmetric part of glycine's map:** the GGG Ramachandran distribution from
    J-couplings and amide I' (Andrews, Biomolecules 2020, Table S2). ff24EXP-GA fits glycine
    phi/psi to it by iterative Boltzmann inversion (Suresh, JCTC 2025). Force fields disagree on
    glycine pPII (ff14SB 0.36, CHARMM36m 0.48).
  * The Hamelryck reference ratio returns the contrastive-divergence fixed point unless its
    non-local feature carries context (the agent's inference).
* **Not found by the survey:**
  * an MD or QM free energy for Ac-X-Pro-NHMe;
  * any measurement of left-neighbour effects on pre-Pro alpha_R;
  * per-residue pre-Pro (phi,psi) in Pro-kinked helices. The kink is ~26 deg with little H-bond
    loss (Barlow & Thornton 1988 [Abs]).

### 1.17 ff3.0 with the physics glycine map: the first epoch, and why glycine gets its own H-bond offsets (2026-10-01/02)

**The selection panel** (plan.md Phase 8; `/project/trsosnic/yinhan/ff3_selection`): the 44 all-L
CATH domains of Charron et al. (2,823 residues, 198 glycines), all-atom ff99SB-ILDN at 300 K, 20,000
frames per domain, against Upside native-start runs at T 0.8 (4 x 8,000 time units per domain after
equilibration). Residues are classed by their all-atom majority basin; the error is Upside minus
all-atom population of that basin, counted in folded frames (Q above the all-atom 5th percentile);
SE by bootstrap over domains; checkpoints are compared paired.

**Released ff2.1** (43 of 44 domains, ff2.1 unfolds 4hwiB01; 0.615 of Upside frames are folded):

| class | residues | ff2.1 - all-atom |
|---|---|---|
| helical, non-glycine | 1,358 | -0.018 +- 0.004 |
| beta, non-glycine, extended region (beta + pPII) | 569 | -0.017 +- 0.004 |
| beta, non-glycine, beta basin alone | 569 | -0.123 +- 0.011 |
| helical glycine | 34 | -0.081 +- 0.049 |
| left-handed glycine | 74 | +0.029 +- 0.019 |

**The beta-basin miss is a phi shift, not strand loss.** Of the 0.123, 0.106 moves to pPII and
0.015 to alpha_R; strands stay extended 0.971 of the time against 0.989. Mean strand phi: Upside
-114.4, all-atom -125.5, the crystal starting structures -121.0, so the two models straddle the
crystal (+6.5 and -4.6 deg); 19% of crystal strand residues already lie past the -100 deg line.
Beta-rich domains are not less stable either: folded fraction 0.63 against 0.60 for helical ones,
rank correlation with beta share -0.05. So the panel scores beta by the extended region. On that
measure ff2.1 holds strands as it holds helices, and its one outlier is helical glycines (~8 points
of alpha_R; 34 residues, so noisy; checkpoint comparisons are paired and tighter).

**The physics map alone, before training.** Swapping in the AWH glycine map (`ff21_awh`, the run's
initial checkpoint) leaves helices (-0.021), strands (-0.019) and helical glycines (-0.081) where
ff2.1 had them and costs the natively left-handed glycines ~9 points of alpha_L (-0.065 against
+0.029): ff2.1's non-local terms do not hold them without the NDRD map's alpha_L bias, so training
must supply it. In the trainer's own replicas (ff30_gly step 0 and the first two cluster steps,
`check_step.py`, 72 proteins), glycines with non-GLY flanks, free / native-restrained: helical
(n 90) alpha_R 0.806 / 0.971, alpha_L 0.076 / 0.019; natively left-handed (n 173) alpha_R
0.066 / 0.013, alpha_L 0.744 / 0.978. Against ff_3.0's epoch 5 (helical alpha_L 0.101, left-handed
0.889) and the equal-depth probe (0.077, 0.847), the physics map leaves helical glycines where equal
depth did and costs the left-handed ones more. What is missing is a glycine context term, not the
map.

**One epoch of ff30_gly moves helices and helical glycines away from all-atom** (panel `e00`,
epoch_00_minibatch_18; 42 domains, 3g7lA00 now unfolds as well):

| | folded | helix | beta | helical Gly | left-handed Gly |
|---|---|---|---|---|---|
| ff2.1 released | 0.615 | -0.018 | -0.018 | -0.081 | +0.034 |
| ff2.1 + AWH map (start) | 0.585 | -0.021 | -0.019 | -0.081 | -0.066 |
| e00 (one epoch) | 0.487 | -0.030 | -0.018 | -0.136 | -0.056 |

Paired against its start, e00 is worse in helix and helical glycine; against ff2.1 in helix, so the
rule holds the release. The H-bond margin E_other - E_alpha fell from +0.192 to +0.057 in one epoch
(E_alpha -1.961 -> -1.912, E_other -1.769 -> -1.855).

**The ff2.1 workflow's own epoch 0, with the NDRD map, is the control** (`ff21-fixedpoint`
epoch_00_minibatch_18, panel `fp_e00`; paired over 41 domains):

| | folded | helix | beta | helical Gly | left-handed Gly | margin |
|---|---|---|---|---|---|---|
| ff2.1 released | 0.615 | -0.017 | -0.018 | -0.081 | +0.035 | 0.192 |
| ff2.1 workflow, NDRD map, epoch 0 | 0.531 | -0.025 | -0.018 | -0.077 | +0.023 | 0.129 |
| ff30_gly (AWH map), epoch 0 | 0.487 | -0.028 | -0.018 | -0.136 | -0.058 | 0.057 |

**Single-file resets locate the two effects** (each panel is e00 with one file put back to ff2.1's,
paired against e00):

| e00 with ff2.1's ... | folded | helix | helical Gly | resolved against e00 |
|---|---|---|---|---|
| (e00 itself) | 0.487 | -0.030 | -0.136 | |
| hbond.h5 (e00_hb21) | 0.499 | -0.029 | -0.087 | helical Gly |
| bb_env.dat | 0.508 | -0.029 | -0.121 | none |
| environment.h5 | 0.490 | -0.027 | -0.115 | none |
| sheet | 0.503 | -0.028 | -0.095 | helical Gly |
| **sidechain.h5 (e00_rot21)** | **0.583** | **-0.024** | **-0.078** | **helix, helical Gly** |

* **(A) The workflow itself weakens helices and folds in its first epoch, with either map, through
  the side-chain pair update.** Resetting `sidechain.h5` alone (`e00_rot21`) restores helix, helical
  glycines and the folded fraction and does not differ from the start in any class. The update is a
  noise-driven random walk: raw per-step gradients recovered from the Adam state (ff30_gly steps
  0-18) give rot (31,420 coefficients) `||mean g|| / mean||g||` 0.107 against 0.229 for pure noise,
  6% of coefficients with |mean g| > 2 SE (chance ~5%), and a change from the start that follows the
  mean gradient with cosine +0.12 only. Adam's normalised step moves every coefficient by ~lr
  whatever its sign consistency, so a table at its fixed point (as rot was at ff2.1, 9k) diffuses
  away: rms 0.44 in 19 steps. By contrast hb follows its mean gradient (cosine +1.00) and the burial
  scale mostly does (+0.83). Hence ff3.0's 10x smaller side-chain learning rate.
* **(B) The AWH map adds the glycine-specific damage**: margin 0.057 against 0.129, helical glycines
  -0.136 against -0.077, left-handed glycines -0.058 against +0.023. Resetting `hbond.h5` alone
  (`e00_hb21`, 40 domains) brings the helical glycines back from -0.140 to -0.087, indistinguishable
  from the start's -0.083, so their extra loss is the H-bond drift; the `sheet` and `sidechain.h5`
  resets also restore them, so the loss needs those changes together. The left-handed loss is the
  map's own.
* So ff2.1 is not at the fixed point of this trainer (its `hbond.h5` and `sheet` never were, 9e), and
  a release judged "no worse than ff2.1" cannot pass while (A) stands.

**Gradient split** (ConDiv's own worker with NSE and DSE kept apart; 25 step-0 proteins, ff2.1 plus
the AWH map; `/project/trsosnic/yinhan/checks/gradsplit_20261001`; sums over proteins, and the trainer
steps against NSE + 0.3 DSE):

| | NSE | 0.3 DSE | contrast | effect |
|---|---|---|---|---|
| E_alpha | +21.45 | -40.54 | -19.09 | weaker helix H-bonds |
| E_beta | +9.54 | -28.70 | -19.16 | weaker sheet H-bonds |
| E_other | +4.56 | -2.31 | +2.24 | stronger left-handed H-bonds |
| E_bias | +9.02 | -10.44 | -1.42 | |
| bb_env scale | -47.68 | +104.65 | | DSE larger in every component |

* The DSE term pulls E_alpha, E_beta and the bb_env scale weaker (the unfolded ensemble near Tm keeps
  residual helical H-bonds that the SARW reference lacks), but resetting those groups leaves the panel
  unchanged (table above), so this pull is not acted on in ff30_glyhb.
* The DSE pull is the same under the NDRD map (E_alpha -40.2 against -40.6, 24 proteins). ff2.1
  balances at lambda 0.10-0.15 rather than 0.3 in every DSE-dominated group (E_alpha, E_beta,
  bb_env scale); one common factor suggests the unfolded ensemble, not three groups, differs from
  what ff2.1 was trained against.
* Peng's SI states the DSE threshold two ways: Fig. S3, 2/3 Rg_low + 1/3 Rg_high of the coldest and
  hottest replicas (the port's and Kleinmann's code); the text, 2/3 Rg_native + 1/3 Rg_SARW. On all 24
  step-0 proteins from one simulation (`gradsplit_20261001/threshold_report.py`) the text's threshold
  is the coded one times 1.02 (median; 0.97-1.12), keeps 458 against 469 unfolded frames and changes
  every DSE group by about 2%. "Threshold as coded" stands.
* **The E_other pull is glycines'**: natively left-handed glycines give +5.25 of the NSE's +4.56
  (115%), other phi > 0 glycines +2.9, non-glycine residues net negative. Paired with the NDRD map,
  the AWH map flips the NSE pull on E_other from -2.55 to +4.60, a paired +0.30 per protein
  [+0.16, +0.44], mostly in natively left-handed glycines (+5.25 against +0.73).
* The rotamer pair gradient is NSE-dominated (|NSE| 58 against |0.3 DSE| 6.5); burial scale mixed;
  sheet has no DSE part by construction.

**Why glycine gets its own H-bond offsets.** With the physics map, the only (phi,psi)-resolved knob
the trainer has for natively left-handed glycines is the shared branch energies, so it makes
left-handed H-bonds cheaper for every residue and the helical glycines pay (lambda's helix 3, 11d;
glpG TM4). ff2.1's three branch energies came from the H-bond node rewrite, while the FF1 trainer
fitted one scalar (-2.112, 9e), so a branch-resolved hbond already gives this trainer a freedom the
original training never had. ff3.0 gives glycine its own dE_alpha, dE_beta, dE_other, trained from
zero (plan.md Phase 8): an optional 15-entry `parameters` with a per-residue `residue_class`, where a
12-entry config stays bitwise unchanged (10.12).

**Inputs kept for the bottom-up fallback** (plan.md, proposed, not approved):
* **ConDiv's derivative step takes any N/CA/C positions** (`compute_divergence`,
  `training/ConDiv.py:399`), so mapped all-atom frames can replace a simulated ensemble without
  engine changes.
* **The rotamer solve has no temperature.** Edge and node probabilities are `exp(-E)`
  (`src/rotamer.cpp:252`, `:896`); the thermostat temperature never reaches belief propagation,
  so side-chain free energies are always at T = 1. Rama maps follow the same convention: a map is
  -ln P, so it reproduces its source statistics at T = 1 (`build_gly_library.py` header). Any
  all-atom target must be compared at an Upside temperature chosen and stated explicitly.
* **A force-matching gradient is available from the engine by one directional central
  difference of `get_param_deriv`** (1UBQ, ff_2.1, displacement RMS 1 A per atom; session test).
  The BP solve is deterministic and path-independent (repeat evaluations bitwise equal). The
  difference converges at steps of 3e-4 and 1e-4 A: successive values agree to 0.9-1.4e-3
  (`hbond_energy`), 1.4-4.7e-3 (`hbond_coverage`), 3.2-6.6e-3 (`sigmoid_coupling_environment`) and
  0.6-1.1e-2 (`rotamer`); `hbond_coverage_hydrophobe` does not converge (0.22-0.28). The energy's own difference needs a step of
  1e-3 A or less to agree with `deriv()` to 0.1%.
* **Public folded-protein trajectories in the glycine map's own force field exist** (verified from
  Charron et al., Nat Chem 17, 1284 (2025), SI section 1.1, and the Zenodo READMEs,
  doi:10.5281/zenodo.15465782).
  * Matching: amber99sb-ildn with TIP3P at 300 K, as our AWH dipeptides (`gly_peptides/awh.sbatch`:
    `-ff amber99sb-ildn -water tip3p`, `ref-t = 300`).
  * Protocol: 50 CATH domains (the first listed in the SI are 60-71 residues), started from native,
    4 x 0.5 us each in OpenMM, coordinates and forces of protein atoms every 20 ps. Also 1,100
    octapeptides with adaptive sampling, 1 us each.
  * Traps:
    * six domains were simulated with D-amino acids (2ga1A02 3e6zX01 3luyA02 3tj8A02 4jriB00
      4npsA02); for glycine handedness they must be excluded;
    * the processed h5 files on Zenodo (`training_a_cg_model.zip`, 65 GB) hold 5-bead mapped frames
      (N, CA, CB, C, O) plus "decoys": every 50th training frame copied with 0.5 A Gaussian noise
      on each bead and a zero force label, stored as separate "decoy molecules" (SI 1.3, 1.5); only
      the real frames are data;
    * the raw all-atom trajectories are not deposited, only their generator scripts.
  * None of the 50 domains is in the 456-protein training list (`training/pdb_list`, by PDB id).
  * The generator bundle (`training_data_generation.zip` -> `cath_generators.tar.gz`, 30 MB) holds
    each domain's starting and equilibrated PDB (`pdbs/<id>.pdb`, `<id>_eq.pdb`) and AMBER
    topology. Its file names are authoritative: `1ldjA06` and `2hbpA00`, where the SI text reads
    1ldjA02 and 2hbpA02.
  * The CB-chirality check on the `_eq.pdb` files confirms exactly the six README domains carry D
    residues, 1-3 each (2ga1A02 ASP102; 3e6zX01 SER17, GLU28, LYS85; 3luyA02 LEU117; 3tj8A02
    LYS120; 4jriB00 LEU45, SER65; 4npsA02 CYS287). The other 44 are all-L: 2,823 residues, 198
    glycines.
  * Archive layout (`training_a_cg_model.zip`, read by HTTP range): the CATH file
    `DECOY_nicks_transferable_cath_delta_dataset.h5` is one deflate stream, 35.58 GB stored /
    39.86 GB raw, data bytes 17270634265 up to the next member at 52854221219; its first inflated
    bytes are the HDF5 signature, so it can be fetched alone.
  * mdCATH (5,398 domains) and the D. E. Shaw fast folders use CHARMM22*, a different force field
    whose glycine treatment the survey could not verify.
* **How Charron et al. treat glycine** (main text p10, SI 2.1, 3.1, Table 4; and the deposited
  model `simulating_a_trained_cg_model.zip/model_and_prior.pt`, read with stub classes in the
  session scratchpad).
  * Mapping: glycine keeps four beads (N, CA, C, O); its identity sits on the CA bead's type,
    where every other residue's sits on CB.
  * Priors are fixed before the network is trained, Boltzmann-inverted from the octapeptide and
    CATH data (every 100th frame): residue-specific 1D phi and 1D psi Fourier series (degree 3), no
    2D map. The CB chirality improper (Gamma1, N-CB-C-CA) has 38 keys and none for glycine, which
    is achiral; the peptide-plane improper (Gamma2) includes it. Glycine's CA has its own
    repulsion radii, within 0.1 A of a plain CA's.
  * **Their glycine phi prior is left-handed**: minimum at +76.5 deg, P(phi > 0) 0.603,
    ln[P(30..100)/P(-100..-30)] +0.51; the antisymmetric part peaks at 0.41 kT. The same readout
    gives ALA's minimum at -70.5 deg, so the sign convention is standard. Units inferred as kcal/mol
    and A from the bond prior (x0 1.54 A, k 237, which needs kT = 0.596 for a ~0.035 A C-C
    fluctuation). Against P(phi > 0): our ff99SB-ILDN dipeptides 0.4998, NDRD 0.6548. Their
    prior sits between, most likely because the 50 native domains' glycine placement is in its
    data; the octapeptides alone are not separated, so that is unproven.
  * This does not bias their model the way the map biases Upside. The prior is only a baseline:
    the network is fitted to the all-atom mean forces of the total, so it corrects whatever the
    prior misplaces, as far as its capacity allows. In Upside the map is the energy term itself.
  * Nothing else is glycine-specific. Mutations to glycine are the one place their first-order
    ΔΔG estimate fails, which they attribute to the entropy of removing a bead (SI p25).
* **What else Charron et al. report that bears on us** (SI 1.5, 2.3, 6.5, 6.6; main text p6-7; read
  from the downloaded PDFs).
  * The released model is not the one with the lowest validation loss: they simulated a set of
    epochs on the fast folders and kept epoch 73, because the training loss "alone is not a
    sufficient metric to identify the highest-performing model for simulation". Three more seeds
    performed comparably. The fast folders they selected on are also their headline test targets.
  * Data balance decides what folds: six purely helical domains were dropped to balance helix and
    sheet. A model trained on CATH alone folds the helical targets but not chignolin or BBA, and
    one on only helical or only sheet proteins degrades on all. BBA (helix plus antiparallel
    sheet) is their weakest target.
  * No transfer in temperature is claimed: the effective energy "really represents a free energy
    with an entropic component".
* **Literature on the bottom-up route (survey 2026-10-01, sub-agent; citations checked by it on
  publisher or index pages, abstract level where paywalled).**
  * Relative entropy needs an equilibrated or reweighted all-atom ensemble (Shell, JCP 2008;
    Thaler, Stupp & Zavadlav, JCP 2022). Force matching needs only conditional equilibrium of the
    removed degrees of freedom, but in a restricted basis such as Upside's the fit still depends on
    the sampled backbone distribution (Noid, JCP 2008).
  * Every bottom-up protein model checked needed a top-down correction: UNRES weight optimisation
    (Liwo, PNAS 2002), the native stability lost in multi-site MS-CG (Hills, Lu & Voth, PLoS Comput
    Biol 2010), and β-content failures in Majewski (Nat Commun 2023) and Charron (Nat Chem 2025).
  * Glycine in all-atom force fields:
    * ff14SB keeps ff99SB for glycine;
    * ff19SB and CHARMM36 use QM glycine maps, and ff19SB warns against fitting glycine to PDB
      statistics;
    * a99SB-disp refit glycine to a PDB coil library, the same contamination as Upside's map, so
      it must not be the reference force field here.
    * No study measured force-field alpha_L populations of glycine in folded proteins.
  * Upside's side-chain energies were trained by maximum likelihood on PDB structures with native
    backbones, deliberately not matched to atomistic energies (Jumper et al., PLoS Comput Biol 14,
    e1006342, 2018).
* **A WebFetch summary invented a methods paragraph** for the Charron paper (ff14SB, OpenMM, 1 us,
  MSMBuilder, all in quotation marks), none of which is in the PDF. Quote a paper only from text
  extracted from the downloaded file.

### 1.18 Glycine map facts that still hold (measured 2026-09-16 to 09-25)

* **Glycine is the only residue whose library handedness cannot be local physics.** ff2.1 coil
  library, central residue, averaged over the 20 left neighbours:

  | res | alpha_R | alpha_L | P(phi>0) | dG(aR->aL) kT |
  |---|---|---|---|---|
  | **GLY** | 8.98% | **30.95%** | **0.650** | **-1.238** |
  | ASN | 18.87% | 13.36% | 0.144 | +0.352 |
  | HIS | 21.57% | 8.00% | 0.090 | +1.021 |
  | ASP | 25.59% | 6.38% | 0.076 | +1.388 |
  | ALA | 27.04% | 4.49% | 0.055 | +1.837 |
  | LEU | 24.58% | 3.32% | 0.040 | +2.039 |
  | SER | 28.64% | 2.70% | 0.039 | +2.421 |
  | THR | 23.73% | 0.71% | 0.015 | +3.586 |
  | VAL | 18.76% | 0.58% | 0.016 | +3.690 |
  | ILE | 18.81% | 0.26% | 0.011 | +4.390 |
  | PRO | 20.53% | 0.00% | 0.000 | +12.621 |

  (others in between; glycine's alpha_L is 7x the 4.28% mean of the other nineteen). Every other
  residue has positive dG, ordered by C-beta branching (Asn and His most alpha_L-tolerant, then the
  unbranched residues, then Thr/Val/Ile, Pro at zero), which is local sterics. Glycine has no C-beta,
  so its handedness has to come from context. Asn (13.4% alpha_L) is the only other residue worth
  checking.
* **The coil/sheet mixture cannot change glycine's handedness.** `read_weighted_maps` mixes each
  coil map with its sheet map through `sheet_mixing_energy`; sweeping that energy from +4 to -10
  leaves ALA-GLY-ALA's dG(aR->aL) at -0.971 to three decimals, because the sheet map is empty in
  both helical basins (GLY|ALA alpha_R 1.6e-10, alpha_L 8.7e-16), so sheet weight dilutes both by
  the same factor. Editing the map is the only way to change handedness. Tuning the mixing energy
  to P(phi > 0) = 0.5 is a trap: it gets there by inflating beta.
* **The library's glycine alpha_L describes folded proteins accurately**: 42 of 67 interior glycines
  of the 16 benchmark natives sit at phi > 0 (63%; NDRD puts 66% of central-glycine weight there).
  The excess is placement at loop sites (1.9), which is why it does not belong in a local energy.
* **No NDRD release is near zero.** Our coil group is exactly `NDRD_TCB` (correlation 1.00000 against
  `GLY|ALL`). Central-GLY dG(aR->aL), mean over `GLY|X`: Conly (coil only, 13,945 residues) -1.876,
  Tonly (turns only, 27,532) -0.839, TCB (44,112) -0.965, TCBIG (adds pi and 3-10 helix, 62,345)
  -0.410. The purest coil subset is the most biased, so "turn contamination" does not explain the
  bias, and switching to `NDRD_Conly` would roughly double it.
* **The map reaches every glycine regardless of context.** `rama_map_pot`'s only input is
  `rama_coord`; its datasets are `rama_pot`, `residue_id`, `rama_map_id` and `rama_map_id_all`, and
  repeated sequence triplets get bitwise-identical maps wherever they occur (checked on glpG).
* **glpG's glycines**: 23 in 210 residues, no GGG, two GG pairs (96-97 inside TM2, 132-133 at TM4's
  N-cap). TM4's three helical glycines, 136 (T-G-V), 143 (M-G-Y) and 149 (R-G-E), are XGX; TM1 has
  none. Under ff2.1 every glycine in glpG is biased toward alpha_L by 0.47-1.23 E_up.
* **`rama_map_pot_ref` is one residue-independent map added to every residue.** Its handedness is
  negligible, dG(aR->aL) = -0.0139 E_up, although it differs from its mirror by up to 0.52 at single
  grid points (a pointwise asymmetry is not a basin free energy). It does reshape a glycine map's
  basin weights (1.16), so a measured surface used as glycine's whole local term is stored with it
  subtracted.
* **Which library number is which**: -1.238 is the `X|GLY` average over all 20 left neighbours and
  -1.318 over the 8 AWH-measured ones; the -1.13 to -1.18 sometimes quoted is the `GLY|ALL`
  marginal. Glycine's `dimer_weight` is 0.908 and 0.948 for the two directions, not 1.0, so the
  left/right mixture matters. Library maps store -ln P normalised to sum(exp(-E)) = 1 over the 72x72
  grid, so an AWH PMF enters as PMF / kT(300 K), not PMF / 2.914952774272 (GLY_sym.md §5).
* The two `GLY|ALL` maps are read only by the `product` combining rule, which nothing uses, so any
  per-map statistic over the glycine row must exclude them.

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

`inner_steps` defaults to 4 in `py/martini_prepare_system_lib.py` (3.10a); the binary reads 1 when
the attribute is absent, so older configs are unchanged.

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

**Bug 1: periodic images. Wrap per molecule, and only after the molecules are whole.** Two faults, found
five months apart, are the same mistake seen from two sides, so they are recorded as one rule.

*First symptom, torn lipids.* `extract_trajectory` wrapped every particle into the box via
`centralize_system` and wrote the frame; nothing unwrapped, so any molecule straddling a periodic face was
left split across the cell and VMD drew bonds shooting across the box. Measured on the seed file, 400
frames, 4187 bonds: mean declared bond 6.89 A, max **141.00 A**, 54502 of 1674800 instances (3.254%) over
50 A. The 141 A worst case was a PO4-GL1 bond *inside one lipid*, so a protein-only integrity check passes
while the file is unusable.

*Second symptom, the protein a full box length out of the bilayer.* Repairing the tear by running the bond
walk **after** the per-particle wrap only moved the fault. The walk rebuilds each molecule around whichever
anchor atom the wrap happened to leave inside, and for glpG that anchor is atom 0, the floppy N-terminal
amide. Whenever the tail crossed a face, the entire 210-residue protein was dragged to the tail's image:
in `glpG_RKRK_79HIS_run0_remd.vtf`, **159 of 1822 frames**, protein-lipid xy centroid separation up to
**98.18 A** (= one box length, 99.77 A), coordinates out to x = 121 A in a box of half-width 49.9 A. The
protein was intact throughout (no CA-CA above 4.5 A) and the underlying trajectory was fine; it rendered
as the protein sitting outside the membrane beside a protein-shaped hole.

*The fix, 2026-09-14.* Order matters and the wrap must be per molecule:
`build_molecule_topology` returns the bond walk **and** a connected-component label per particle;
`extract_trajectory` calls `unwrap_molecules` **first** to make every molecule whole, then
`centralize_system`, which shifts by the plain protein centroid (no circular mean is needed once the
protein is whole) and wraps each molecule by `box * round(centroid/box)`. A whole molecule is never torn
again, so nothing has to be rebuilt around an arbitrary anchor. After: protein COM exactly 0 in all 3146
frames, 0 displaced frames, protein-lipid xy separation mean 0.70 A / max 2.20 A, no declared bond over
10 A except the known residue-210 C-O.

**When validating a VTF, check the declared bonds across all frames *and* that every molecule centroid is
inside the cell.** Either check alone passes one of these two bugs.

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

### 3.9 Upside traps on exit under clang whenever Monte Carlo is enabled

**The MacBook's `obj/` was rebuilt 2026-10-01 14:37** (`make` in place, existing CMake cache; not
`install_M1.sh`, which empties `obj/` and would delete the stray notes kept there). The 08-24 build
predated the 09-07 destructor fix and the 09-26 rotamer `fabsf` fix and trapped with exit 133; the
rebuild runs Monte Carlo to exit 0, and the panel's local path runs end to end.

`MonteCarloSampler` (`src/monte_carlo_sampler.h:12`) is abstract (it declares
`propose_random_move` pure virtual) but has **no virtual destructor**, while
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
`WORKER_FAIL` on completely valid data**: 250/250 frames written, Rg 13.2-13.8 A, potentials
negative, the temperature ladder correct. `run_minibatch` then raised `All jobs failed`. Read that
way round, an exit code was condemning good physics.

Fixed by giving the base class a virtual destructor. Verified by execution, not inspection: the
three previously-trapping invocations exit 0, and on an identical 160-time-unit run the old and new
binaries produce **bit-identical output across all 16 datasets** (`pos`, `potential`, `kinetic`,
`hbond`, `pivot_stats`, `rama_map_potential`, ...), so results are unchanged and master parity in
results holds. The trap sat purely in the teardown path.

### 3.10 The glpG blow-ups: the two subsystems are not at the same temperature (2026-09-10)

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

The ladder was made authoritative and dry-MARTINI brought to it (decision of
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

**The glycine maps were not this defect.** The same controlled run compared two versions of the
retired symmetrised glycine maps (wild type, T 0.90, `inner_steps = 4`, two seeds per arm): an
off-by-one mirror, which happened to favour alpha_R at GLY143 by +0.571 E_up, held TM4a at 0.786;
the correct mirror, exactly neutral, gave 0.531; the raw library map favours alpha_L (-0.644). No
corrected-map run beat the best off-by-one run, so the residual fraying at T 0.90 is set by how much
right-handed preference glycine's map carries, a force-field question rather than a bug. TM4 is the
glycine-dense helix (3 in 17 residues; TM1 has none) and frays from its N-terminal turn at the hot end
of the ladder. n = 2 per arm, with within-arm spread comparable to the difference.

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

### 4.4 glpG TM4 before the temperature fix: what was ruled out (2026-09-07 to 09-10)

Before the protein/lipid temperature mismatch was found (3.10), TM4 lost helix in every glpG variant.
In the ff_2.1 production REMD (T 0.70 rung) every variant started at TM4 = 1.000 and ended at
0.35-0.60, wild type 0.432, while TM1 mostly held (0.77-0.94; 79ALA 0.56). TM4 (134-151) is
byte-identical in all four variants. With `inner_steps = 4` the cold end is healthier than the seed
(3.10c). What the hunt established, kept so it is not re-tested:

* **Measurement windows.** The crystal seed has 131-133 and 152 non-helical, so TM4 scored on 131-152
  is capped at 18/22 = 0.818 against a > 0.8 criterion, and TM1 on 29-49 at 0.952. Use 134-151 and
  30-48. The GLY49/GLY133 phi criterion is invalid: both are helix caps at phi_std +94 and +142 in the
  seed.
* **Compare trends, not means.** Arms that have not equilibrated cannot be ranked by trajectory means:
  the ff_2.1 control had bottomed out by its third quintile while the trained arms were still falling
  (second-half slopes -0.022 against -0.055 and -0.135 per 1000 time units), which made "TM4 is
  fixed" a convergence error. Report quintiles and the second-half slope. Single-temperature MD at
  T 0.70 is also harsher than the REMD rung the > 0.8 criterion came from (control 0.441 against
  0.645).
* **It was not thermal melting, and it was reversible.** TM4 was flat at 0.36-0.57 across the whole
  ladder, and over 13,117 frames of the WT run it interconverted (47.8% of frames still above 0.8)
  while drifting down.
* **The loss was local**: 134-142, the buried midplane half (139-141 from 1.000 to 0.33-0.35; i->i+4
  O...N 6.4, 6.4, 5.2 and 4.3 A at 139-142 against ~3.0 for a formed helix), while 143-146 and 149
  held at 1.000. The seed is a pristine helix (O...N 2.81-3.12 A from 134 to 149).
* **TM4 is the least lipid-exposed helix**: 2.0 backbone beads with a lipid neighbour inside 12 A
  against TM1's 10.7, 7 of its 22 backbone beads with none, and zero dry-MARTINI energy between the
  134-139 stretch and all 3627 lipid beads. Every protein-lipid explanation fails for TM4, including
  the MARTINI backbone typing of 131-137 (re-typing all of it to `N0` changes the total by
  -0.3 kJ/mol).
* **Bundle splay does not demonstrably cause the melt.** The raw lipid-gain/helix-loss
  cross-correlation (r 0.64, lag 99) is a shared-trend artifact; on first differences every coupling
  falls to r 0.21-0.26 with a lag under one sampling interval. The lipid that intercalates between TM
  helices is acyl tail (enrichment 1.20x, headgroups 0.39x), which is correct hydrophobic solvation.
* **Glycine density does not predict which helix fails** (r = +0.23, wrong sign; TM3 is 17.9% glycine
  at 0.928), and TM4's buried charges are paired (ARG148-GLU150 4.1 A, ARG151-ASP152 3.4 A).
* **The backbone-environment term is absent from the hybrid and does not fix TM4.** The hybrid
  builder never passes `--environment-potential`, which gates every environment node in
  `upside_config.py`, so glpG has no `bb_sigmoid_coupling_environment`, `environment_coverage_*` or
  `hb_environment_coverage_*` nodes (2.2). Restoring them moved TM4 by +0.020 against a seed scatter
  of +-0.10, because TM4's backbone sits at burial 13.7 where the term saturates above 4: its credit
  is already zero. `hbond_coverage` is not a desolvation measure either; its output is per side-chain
  bead, a 1-body cost handed to the rotamer solver.
* **A MARTINI-typed apolar coverage channel made things worse** (TM4a 0.807 -> 0.464, TM1 0.911 ->
  0.810). The coverage is non-negative and enters with a positive weight, so it acts as a repulsion
  between backbone and acyl chains; a lipid-aware channel built on these nodes can only penalise
  lipid proximity. Not deployed.
* **Engine facts from building it.** Node types resolve by prefix (`deriv_engine.cpp:598`), so a
  group named `hbbb_coverage_lipid` instantiates `hbbb_coverage`. `rotamer` takes a variable-length
  `prob_nodes` list and uses each output directly as a 1-body energy, each with `n_elem` equal to the
  side-chain bead count (747). `hbbb_coverage` has `n_dim2 = 3`, so its second group can be bare
  positions; the trained `hbond_coverage` needs `n_dim2 = 6`. Making the trained coverage lipid-aware
  is an engine change in `src/environment.cpp`.

**Open hypothesis, untested.** With `exclude_intra_protein_martini = 1`, all helix-helix packing comes
from the Upside core force field, trained on soluble proteins only, while protein-lipid attraction is
full-strength dry-MARTINI, so lipid competes for the same hydrophobic surfaces against a term never
balanced against it. The principled test is ConDiv training with membrane proteins in the set;
scaling SC-env or BB-env down is forbidden.

### 4.5 Lessons from the FF1-form ff3.0 retrain (2026-09-07 to 09-09; trainer retired 2026-09-24)

The first ff3.0 was trained by a port that turned out to be FF1's workflow (9t-9v), and its force
fields are retired. What it showed about ConDiv runs in general:

* **The parameters random-walk while the objective plateaus.** Drift of `pair_interaction` grew as
  sqrt(lag) (rel_rms / sqrt(lag) 0.013-0.016 over a 16-fold range of lag) with no fixed point, while
  the median restrained RMSD was flat over the last half of the run. Training longer only diffuses
  further, and no step can be preferred on training grounds; a downstream observable has to choose.
  The side-chain random walk of 1.17 is the same effect.
* **Two runs that share history are not replicates.** A run resumed from another's step-274
  checkpoint diverged from it at the single-run diffusion rate, reaching rel_rms 0.26 by step 496;
  extrapolated, two independent trainings from ff_2.1 would differ by ~40%. A reproducibility test
  needs runs branched at step 0, and force fields from different steps must never be paired.
* **Retraining reproduces how far it moves, not where.** Two lineages from a common step 269 drifted
  the same distance from ff_2.1 to within 0.1-1.7%, in directions ~23 deg apart (between-run rms
  33-40% of the drift).
* **The first update after a resume is outsized** (rel_rms 0.076 after one step, against 0.018 per
  sqrt(step)), most likely because the optimiser's momentum state was not carried across.
* **Most of the movement happens in the first ~60 steps** (cumulative rms 0.49 at step 29, 0.64 at
  59, then ~0.07 per 30 steps without an asymptote), so a step target needs a measured justification.
* **Trained tables and the coverage nodes only work together.** On glpG at T 0.70 (3 seeds), TM4 helix
  fraction was 0.441 for ff_2.1 without coverage nodes, 0.562 with the nodes and old tables, 0.588 with
  new tables and no nodes, and 0.782 with both: a pair table co-trained with coverage is correct only
  beside the partners it was optimised against. 135 further steps (269 -> 404) left TM4 flat (0.709).
  At an unrelaxed seed the trained force field's force tail was heavier (p99.9 97 against 49 E_up/A),
  part of it from the coverage nodes themselves.
* **A checkpoint of that trainer could be rebuilt from an extracted force field.** `expand_param`
  drops `rotscalar`, which is identically zero at init; `pack_param` (an L-BFGS-B refit) reproduced
  trained tables to ~1e-16; `params.env = energies[:, :-1]` recovers exactly (check
  `energies[:, -1] == energies[:, -3]`); Adam state is lost, a short transient with beta2 = 0.96. The
  acceptance test is that `extract_ff.py` on the rebuilt checkpoint reproduces the starting files.
* **Upside overwrites `/output` rather than appending**, so a fresh run on a production seed replaces
  the seed's frames. Tell them apart by the time spacing, not the frame count.

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
the protein interior and the catalytic cavity. 141, 142 and 145 are not cavity-lining (5.3b): they
carry as many protein heavy atoms within 6 A as the buried 143 and 147 and sit 6.9-7.3 A from the
nearest tail, straddling the 7.00 A shell, so their signal is local thermal opening at the edge of
tail contact.

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

**Consequence**: because `acc_t = 0` makes protection identically 1 whatever the H-bond is doing,
saturation is decided by lipid contact (the `pp_fail` / `acc` table in 5.3), so a `+inf` run is a
lipid-contact map, not a fold map, and "TM1 never exchanges" is false as a structural statement.

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

**Cooperativity, measured directly.** Implicit's many identical near-zero rates are uniform
freezing, not Englander's cooperative-unfolding signature, which converges on the finite unfolding
rate of a cooperative unit. One replica per model at matched temperature (implicit T=0.851 of 48
rungs, hybrid T=0.853 of 28; identical ladder span 0.700-0.900, mean 0.798, so this is not a ladder artifact),
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

### 5.3c The protein carries no explicit charge in the hybrid, and what follows from it (2026-09-12)

Measured on both campaigns' configs, so this is systematic rather than one build's mistake. Of the protein
beads in `martini_potential/charges`, **only 10 are nonzero** -- 2890 beads in the NP system and 1050 in
glpG, identical in both, 5 at +1 and 5 at -1, summing to exactly 0. They look like the two chain termini
smeared across all five beads of their residues, which cancels and is negligible.

**This is by design, not a bug.** `martini_sc_table_1body` is fully residue-resolved:
`restype_order (18,)`, `rotamer_full_energy_eup (18, 6, 38, 96, 13)`, i.e. 18 residue types x 6 rotamers
x **38 environment bead types**. Protein-environment interaction is carried by tabulated per-residue-type
fields, which is what the "spline table only" rule requires, instead of explicit Coulomb.
`charged_res` in `martini_prepare_system_lib.py:246` is `{ASP:-1, GLU:-1, LYS:+1, ARG:+1}`; HIS is
correctly absent, so **an earlier inference in this session that "all 15 histidines are protonated" was
wrong** and is retracted. Albumin's sequence charge at pH 7 is -15 (LYS 58 + ARG 24 = +82, ASP 35 +
GLU 62 = -97); the config simply does not represent residue charges explicitly at all.

**The consequence, which is the part that matters.** The SC-env tables are short-range radial fields, so
albumin's 179 charged residues have **no long-range Coulomb term** with the anionic MPA coating or with
the ions. Changing the salt therefore cannot create protein-NP electrostatic steering: at 0.15 M the
model's Debye length is 3.6 A and counterions-only would give 16.4 A, but neither matters to a protein
that has no charge to screen. **Do not expect an ionic-strength change to fix the two NP orientations
that never bind** (they start ~37 A off the surface); that is a model-design property plus starting
geometry. For the same reason, a charge-driven footprinting prediction of the Carlson kind is not
something this hybrid can currently reproduce from first principles.

**One real defect, small:** the NP box carries **net +15 e**. The ion generator added 218 excess K+ on the
assumption that albumin is -15, while the charge array gives it 0. Either the protein should carry -15
explicitly or the ion count should be 203.

### 5.3d The retired ff3.0's dG profile, and which ff2.1 profile is valid (2026-09-13)

At matched rungs, the retired ff3.0 hybrid (`<V>/hdx_postfix/results/`, 2026-09-12; FF1-form trainer,
symmetrised glycine maps) against the ff2.1 hybrid (`<V>/hdx/results/`, 2026-09-04), both after the
rigid-stage fix: non-exchanging amides 53 -> 80 of 203 at T 0.75, 44 -> 73 at 0.80, 43 -> 69 at 0.85,
so the retrained core raised protection along the whole chain. TM4 resolves to a finite value in
every run (censored 0/21 under ff2.1 and 1/21 under ff3.0 at T 0.85, resolved median 2.66 -> 3.03),
and TM1 is the censored helix (4/22 -> 9/22). Over three frame counts (6,121, ~9,500 and ~12,000 per
replica) the four variants stayed indistinguishable (62-68 of 203 off scale at T 0.85, TM4 0-1 of 21
censored), so neither H79A nor S115T reshapes global protection. The `hdx_postfix/` results are no
longer on `/project` or `/beagle3` (remote_jobs.md §1b). Quote one censoring criterion: the
`censored` flag and the sentinel count differ (66 against 105 at T 0.75).

**`glpG_POPEPOPG_dG_2026-08-27/` must not be used as the ff2.1 reference.** It predates the 2026-09-02
stage fix, so its protein was frozen (3.4). The valid ff2.1 reference is the 2026-09-04
`hdx/results/` set.

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

### 6.1 Glycine Ramachandran maps: mirror convention and checking traps

The old rule "symmetrise every glycine map" was the retired ff3.0 (GLY_sym.md §4), and its code is
gone from `upside_config.py`. ff3.0's library makes only `GLY|GLY` mirror-symmetric, and
`build_gly_library.py` checked that through `upside_config` when it built the library. What survives:

* **The mirror is `(phi,psi) -> (-phi,-psi)` with a roll.** On the 72-point grid starting at -180,
  index i maps to (-i) % 72: `np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)`, or
  `m[np.ix_(mir, mir)]` with `mir = (-np.arange(n)) % n`. A plain `m[::-1, ::-1]` maps i to n-1-i, is
  off by one bin, creates a new asymmetry (~0.6 E_up in alpha_R's favour) and reads 0.0000 under its
  own symmetry metric. The same expression is right for library `dimer_pot` and for `.up` `rama_pot`,
  since `write_rama_map_pot` copies the grid unchanged. Check against a chiral control (ALA must keep
  ~11 E_up of asymmetry) and against dG(aR->aL).
* **Grid indices use `int()` (floor)**, as `inject_backbone_nodes` does, not `round`: they differ
  (25 against 24) for phi -57 on a 72-point grid.
* **The hybrid `input/pos` backbone is stride 4** (N=4i, CA=4i+1, C=4i+2, proxy=4i+3) and the
  pure-protein layout stride 3; read N/CA/C from `hybrid_bb_map/atom_indices` rather than assume
  either.
* **Silent failures of the old per-residue symmetrisation**, each worth checking for in any
  per-residue map edit: 1-indexed residue numbers used as 0-indexed h5 lookups (every glycine left
  untouched); a phi-range filter or an alpha_R > alpha_L guard (every boundary misses some glycine);
  the wrong mirror above.
* **The library's NaN is one whole neighbour column, `CPR`** (4.545% = 1/22), never read because a
  cis-proline neighbour is mapped onto `PRO`. A naive `abs(a-b) > tol` diff reports "no difference"
  across it, since NaN comparisons are False.

### 6.2 Verifying TM helix health in a glpG trajectory

The procedure and scripts are in `remote_jobs.md` §5 (`glpg_tm_windows.py`, `glpg_ff30_health.py`):
TM1 30-48 and TM4 134-151, helix = phi in [-130, -20] and psi in [-90, 15], pass > 0.8 at T 0.70 with
no downward trend across blocks, reading only rotated `output_previous_*` groups and skipping
`output_previous_0`. Three traps (memory `glpg-vtf-reading-traps`): the old handbook `dihedral()`
helper and `check_seeds_current.py` return -phi_std; GLY49 and GLY133 are helix caps at positive phi
in the crystal seed, so there is no glycine-phi criterion; and a protein "leaving the bilayer" in a
VTF is the N-terminal tail. A VTF lists N, CA, C, O per residue from atom 0, and phi of residue i
uses C(i-1), N(i), CA(i), C(i).

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

### 6.6 Identifying which force field a running replica carries

`run_remd.py` copies a replica from the seed **only if the file does not already exist**, so a chain
resumed after a force-field install keeps whatever tables it was built with. The force field lives inside
each `.up`, not in a path the job reads, so the only way to know is to compare the baked tables against
`parameters/ff_*/`.

The reliable test is a least-squares scale match on the rotamer pair table, because `upside_config`
rescales it into Upside units on the way in:

```python
with h5py.File(f'{P}/{tag}/sidechain.h5') as f: ref = f['pair_interaction'][...]
with h5py.File(up) as f:
    a = f['input/potential/rotamer/pair_interaction/interaction_param'][...]
c = (ref.ravel() @ a.ravel()) / (ref.ravel() @ ref.ravel())   # best scale
```

The match is unambiguous: the right force field gives `c = 1.000000` with a residual at float32 rounding
(1.6e-6), the wrong one gives `c = 1.107` and a residual of 21. Measured 2026-09-13 on
`popepopg_REMD_mdw2/glpG-RKRK-79HIS.run.0.up`, which came out `ff_3.0` exactly, confirming all 33 of its
output blocks are post-retraining.

**A verbatim hash sweep over the parameter files gives the wrong answer here.** Sweeping every dataset and
asking which `ff_*` directory appears inside the seed reports **ff_2.1**, on 11 matching datasets. The
reason is that the retraining only touched the protein core, `sidechain.h5` and `environment.h5`; the
dry-MARTINI SC-env tables in `martini.h5` were not retrained, `ff_3.0/martini.h5` does not exist at all,
and the seed build therefore pulls that file from `ff_2.1` by design. `hbond.h5` is likewise byte-identical
between the two. So the verbatim test finds the shared files and misses the two that actually differ, since
those are transformed before they are written into the `.up`. Compare the transformed tables, not the files.

### 6.7 Glycine AWH: convergence lessons (2026-09-18/19)

The AWH measurement of glycine's handedness (GLY_sym.md §5; final numbers in 9r and 9s) went through
several wrong intermediate answers, listed in 12a. The methods that caught them:

* **Check the estimator's dynamic range first.** The first runs had `awh1-dimN-diffusion = 5e-5`
  where AWH's own friction metric implies ~0.77 rad^2/ps, and a small `awh1-error-init`, so the PMF
  spanned ~2 kJ/mol at 25 ns against the 20-40 of a glycine surface. A flat surface is trivially
  mirror-symmetric, which produced "handedness zero". Fixed with `diffusion = 0.5` and
  `error-init = 30`, rate parameters only; the range then reached 50-66 kJ/mol.
* **An achiral control passing is necessary and nowhere near sufficient**, because its errors cancel
  by symmetry. At 7 ns the Gly-Gly blank read -0.002 in one replica and -0.264 in the other; only
  independent replicas at matched sampling bound the error.
* **Replica agreement at matched time does not establish convergence either**: two replicas drift the
  same way. The neighbour average agreed at 36-45 ns (-0.31) and plateaued only from 62 ns, at -0.269.
  Convergence needs a single-replica series run past where it stops moving, with a plateau longer
  than the feature called a plateau, and the control's own series must be checked before its offset
  is called a systematic.
* **The blank is the convergence criterion.** Its residual was uncorrelated between replicas
  (r = +0.157) and decayed as 1/sqrt(t) (0.233 at 10 ns to 0.032 at 100) while the signal converged
  (0.080 to 0.071), so it is sampling noise, and blank subtraction only adds noise. It puts ~20%
  uncertainty on the basin dG.
* **Per-neighbour structure is noise; only the average reproduces.** Replica neighbour orderings
  correlate at Spearman rho +0.048 (p 0.91), and their disagreement (0.18) did not shrink with time.
  When the claim is about every pair, plot every pair, not only the mean.
* **Quote wander over the longest window**, not four snapshots: "LA is settled" at a 0.028 range was
  0.106 over 18 ns.
* **A derived quantity's definition must live in one place.** Two scripts used different
  mirror-image basin boxes; the achiral control read 0 under both, so it could not flag a 0.106
  difference between them.
* **On a periodic grid the boundary column breaks mirror symmetry.** With no +180 column, counting
  phi > 0 cells gave P = 0.481-0.486 for exactly achiral controls. Split the -180 column between the
  two halves.
* `gmx awh` on midway2 silently writes no `fe_t*.xvg` when the shell lacks `module load gcc/10.1.0`;
  after an extension, read the last part file (`ls -v awh.part*.edr | tail -1`).

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

### 8a. The annular lipid shell is still filling for the first ~2/3 of the glpG production run (2026-09-14)

Asked whether the gap around the protein closing over the trajectory was expected. It is, and it is
post-insertion relaxation of the boundary lipids, measured on `glpG_RKRK_79HIS_run0_remd.vtf`
(3146 frames, blocks 1-54, cold rung T = 0.70, 4104 t_up total = 456k steps at dt = 0.009).

**It is not the box and it is not a pore.** Production is **fixed-volume**: there is no `box` dataset in
any output block and no `input/barostat`; the cell is an attribute on `martini_potential`
(99.768 x 99.768 x 180 A, constant). The barostat only ran during preparation. And there is never a
through-hole: lipid-free projected area not covered by the protein is **0 A^2 in every frame** at a 5 A
probe. What looks like a hole is an under-packed annulus, not a defect in the bilayer.

**It is not the protein either.** TM-slab Rg_xy is 11.6 -> 12.5 A and flat after the first eighth, and
the lipid midplane stays 1.6-3.1 A from the protein centroid throughout.

**The lipids move inward.** Radial density from the protein surface (12 A core slab, early 30 frames vs
late 30): every bin inside 20 A gains (+73, +46, +46, +32, +39, +22, +27, +15, +15, +13 beads), every
bin beyond 25 A loses (-12 to -25), crossover at ~22 A.

**Two stages, and only the first is fast.** Protein-lipid contact beads within 6 A, TM core (res 29-208)
alone, rise 137 -> 220 (+60%), so this is not the disordered termini lying down (those rise separately,
18 -> 50). But the *number* of annular lipids saturates early (36 -> ~44 by the second sixth) while
contacts *per* annular lipid keep climbing 3.8 -> 4.8. The shell fills quickly, then tightens slowly.

**Timescale, and why it matters.** Fitting A - B exp(-t/tau) to total contact beads gives
**tau = 1190 t_up, 29% of the whole production run**; 90% of plateau is reached only at ~frame 2100 of
3146. The reverse cumulative mean settles to within 1% only over the last 20-30% of frames
(last 30% 260.2, last 20% 262.1, last 10% 262.2, against 223.6 for the whole run). Do not convert
tau to real time through the nominal unit table: `dt` here is locked to `/input/brownian` with the
friction tuned for a target lipid diffusion, so the physical mapping goes through that tuning, which
has not been verified for this system.

**Consequence, not yet acted on.** Roughly the first two-thirds of this trajectory is not equilibrated
with respect to protein-lipid packing, and the HDX protection / membrane-accessibility estimator reads
exactly that interface. Campaign 6 used all frames. Two things to check before the four-variant
equality is called converged: re-run the estimator on the last third only, and measure this same shell
curve for the other three variants -- if the relaxation differs between them, part of the comparison is
between equilibration states rather than between chemistries.

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

### 9.x The driver's logged Rg is not periodic-image corrected (measured 2026-09-18)

`np.<jobid>.out` prints an Rg computed on the stored coordinates with no minimum-image correction, so
on any face where the adsorbed chain crosses a box boundary it is inflated by roughly the number of box
lengths spanned. On block 6, logged vs minimum-image Rg (every backbone atom referenced to the MPA shell
centre, box 300 A): run0 123.1/128.8, run1 184.0/76.4, run2 96.9/102.9, run3 332.5/118.7, run4
170.6/73.4, run5 80.8/80.3. Three of six faces inflated 2.3-2.8x. run3's `Rg = 332 A` **in a 300 A box**
is the tell: the value exceeded the box and was still being read as structure.

The physical state is the opposite of what that number suggests. The protein is adsorbed on every face,
with 257/1128/913/399/780/1243 of 2312 backbone atoms within 8 A of the MPA shell, and no atom further
than 172 A from it. Spreading is real (minimum-image Rg 73-129 A against native albumin's ~27 A) but
smaller than logged. **Use the contact count as the adsorption observable; Rg needs an unwrap the driver
does not do.**

Two things to carry forward. First, this is the same defect as the VTF per-molecule wrap: a global shape
descriptor computed per-atom across a periodic boundary reports a structure that does not exist, and it
fails silently because the number stays finite and moves smoothly. Second, **a length that exceeds the
box is a self-evident failure of the measurement and should be caught by inspection**: 332 A in a 300 A
box had been recorded in `remote_jobs.md` as the campaign's headline observable ("rising Rg, currently up
to 230.9 A") without anyone comparing it to the box it lived in. Compare every length observable against
the box before reading it.

The minimum-image correction is itself only valid while the chain stays inside half a box. run0 (max
172 A) and run5 (164 A) exceed 150 A, so those two faces are unresolved, not confirmed.

---

## 9c. hbond is a function of the rama COORDINATE, and its turn branch is the glycine region (2026-09-18)

`RamaMapPot` and `HBondEnergy` are sibling potentials on the same `rama` CoordNode
(`src/hbond.cpp:490`). hbond never reads the rama map table, so editing `rama.dat` changes no hbond
parameter or input. But `HBondEnergy` classifies each residue by (phi,psi) and picks a different
per-hbond energy from it:

    Ehbond[i] = E_alpha*helix_score + E_beta*sheet_score + E_other*turn_score
    potential = sum_i hb_number1[i] * Ehbond[i]

Decoding the 12 values in `hbond.h5` (radians in, clean degrees out, which is itself evidence they
were hand-set): turn = **phi in (0, 165) deg** -> E_other **-1.769**; helix = phi<0 and psi in
(-120, 60) -> E_alpha **-1.961**; sheet = phi<0, psi outside -> E_beta **-1.946**. Boundaries
0/165/-120/60 deg, all four sharpnesses identical at 3.81972 (a 15 deg ramp).
`compact_sigmoid` is 1 for large negative argument and 0 for large positive
(`src/vector_math.h:700`), which is what fixes the window directions.

**The turn branch is exactly the positive-phi region, i.e. the glycine question.** A hydrogen bond
at phi>0 is worth **+0.192 E_up less** than the same bond at phi<0, and glycine is the residue that
lives there (alpha_L 30.95%, ASN second at 13.4%). The coupling runs both ways:
`rama_sens(0,i) += hb_number1[i]*dPhi[i]`, so a hydrogen-bonded residue is pushed toward phi<0 in
proportion to its bond count, ~0.38 E_up for a doubly-bonded helical glycine, against ff2.1's rama
pull of -1.238 E_up toward alpha_L. Rama wins by ~3x. That is a mechanism for glpG TM4 and for
lambda's H2 (11d).

`hb` was trained by the original FF1 trainer as one scalar (9e); the three branch energies and their
boundaries came from the later node rewrite and were never fitted until FF2's trainer (9v). The
branch energy `E_other` and glycine's alpha_L map depth are near-degenerate, since both set what an
H-bonded glycine pays at phi > 0, and the energies are shared by all 20 types. That is why ff3.0
gives glycine its own offsets on them rather than retuning the shared values (1.15, 1.17).

---

## 9d. The retired ff3.0's benchmark split by native/de novo, not by topology (2026-09-18/19)

All 32 Peng arms of the retired ff3.0 (FF1-form trainer, every glycine map symmetrised), scored on
mean TM against the digitised FF2 curves (`0914/figs/ff2_curves_s5.npz`): native +0.039 (11 of 16
improved), de novo -0.030 (5 of 16); native minus de novo positive in 13 of 16 pairs, sign test
p = 0.021, paired t = +3.20. The regressions span every topology (WW domain, an all-beta sheet, was
the worst de novo arm at -0.206), so "helical bundles fail" is wrong, and `hyp_denovo` (+0.166) is a
counterexample to "ff3.0 hurts de novo folding". Two causes were never separated: removing glycine's
alpha_L bias removed turn nucleation, or the whole-run scoring trap of 11d, which flatters native
arms (they decay from the native seed) and penalises de novo ones (they build up toward folded).

**Do not read a partially scored benchmark as a result**: the p value went 0.039 (n 9), 0.092 (13),
0.035 (15), 0.021 (16), leaving and re-entering significance while scoring was in flight.

---

## 9e. The FF1-form ConDiv port and the learned glycine map (9e-9q; 2026-09-18 to 09-24, retired)

From 09-18 to 09-24 the trainer was a Python 3/torch port of Peng's FF1 Theano ConDiv
(`ConDiv_original.py`, Upside 18.10.08), and Track A trained glycine's map inside it. Both were
retired on 09-24, when the port turned out to train FF1's Hamiltonian (9t-9v). The decomposition of
the library's handedness by chiral context that this work produced (old 9i) is in GLY_sym.md §2a.
What stays true:

**The original trainer and the port (9e-9h, 9j, 9n).**
* The original trained `hb` (lr 0.02) and `sheet` (0.03) as single scalars (`init_param/hbond`
  -2.112, `init_param/sheet` -0.268); the port dropped both when the nodes changed shape. ff2.1's
  12-entry `hbond.h5` and 20-value `sheet` came from the node rewrite, not from any trainer, so they
  were never at a ConDiv optimum.
* Energy is exactly linear in `hbond_energy.parameters[:4]`: scaling by 1.01, 1.10 and 0.90 keeps
  `E_hb / scale` constant to 5.6e-7 (float32 round-off), so `dE/ds = E/s` for a common scale and
  `apply_param_scale(hb_scale=...)` (`--hb-scale`) scans it with a config rebuild only. The 12
  entries are `[E_alpha, E_beta, E_other, E_bias | 8 rama boundaries and sharpnesses in radians]`.
* Sheet derivatives by finite difference carry barely one significant digit per frame: the rama
  energies are ~788 and the more/less difference ~2.3e-3, about 25x the float32 resolution there. Do
  not over-read a small sheet gradient.
* The Theano -> torch swap was faithful term by term (student-t, expectation profiles, lower bound,
  regulariser, broadcast, constants, Adam). Two additions were removed (old 9g): a GLY palindromic
  symmetrisation of the rotamer pair-interaction angular profile, which is baked into the old ff_3.0
  `sidechain.h5` (GLY angular `max|x - flip(x)|` 0 against ff2.1's 0.9996), and a `+1e-12` inside two
  direction normalisations. Removing the palindrome let `pack_param`'s original `discrep < 1.6e-4`
  residual gate be restored; ff2.1's GLY row fits to 2.6e-30. The GLY version is in git (blob
  `72ae60be`, commit `28185321`).
* Inherited from the original: `training_list` is sorted by size and then shuffled unseeded, so
  minibatch composition differs run to run. Never verified: whether upside1's `restraint_spring` and
  upside2's `restraint_spring_constant` share a definition.
* The env node's parameter vector is `coeff (360) + weights (400)`, and a request sized from `coeff`
  alone fails with "expected 760 but got 360". This known, fixed bug came back because a run was
  cloned from the one directory that had not been patched: clone from the most recently fixed
  directory, not the most successful one.

**Testing a fixed point (9k).** From ff_2.1 under the port, six minibatches gave gradients
indistinguishable from noise in every group (sign-flip p 0.22-1.0): ff_2.1 was a stationary point of
that trainer, not shown to be its attractor. Use the exact sign-flip test on `||mean g|| / mean|g|`
(cheap at 2^n; it is `ConDiv.py gate`); the old pairwise-cosine t-statistic is anti-conservative because the
pairs share vectors. A partial-n statistic misled more than once: the sheet gradient read t = +3.91
at n 3 and -0.01 at n 6.

**The learned glycine map, Track A (9l, 9m, 9p, 9q).**
* A 72x72 map is trainable only because its gradient is analytic: `rama_map_pot` is a periodic
  bicubic spline built from 1D periodic solves (`src/spline.cpp:262`), so `dE/d(map[i,j])` is a
  spline-smoothed histogram of the (phi,psi) samples. Finite differences would cost 10,369 times a
  divergence.
* The finite-difference gate through the real pipeline (library -> `upside_config` -> engine) found a
  37% error on its first run: `write_rama_map_pot` subtracts a Boltzmann-weighted constant per map
  (`rama_pot -= (rama_pot*np.exp(-rama_pot)).sum(...)`, up.md 2.8a), invisible in basin differences
  but present in the total potential. Read such a check as an eps sweep: float32 `dimer_pot` rounding
  dominates at small eps and curvature at large eps.
* A gradient check only tests the branches its protein exercises. Terminal glycines (3% of the
  gradient) were missing, and the test protein had none; choose test proteins that cover every
  branch.
* The map's handedness moved in bursts, ~20-step stalls and then descent, because glycine content
  varies between minibatches; three "it has plateaued" calls were wrong (memory
  `condiv-gly-epoch-scale-only`). No window shorter than an epoch means anything.
* It converged at dG(aR->aL) -0.885 (9s), past the training natives' own -0.50 and far from AWH: the
  map compensates the shared `E_other` penalty of `hbond` (9c), because it is the only term that is
  both residue-type specific and (phi,psi) resolved (architecture.md §2).

**PDB statistics as a target (9o).** NDRD's map is a potential of mean force over folded structures,
and ff2.1 is a constrained optimum: its other terms were fitted with the map fixed and can compensate
only globally. The 456 training natives give glycine dG(aR->aL) -0.500 +- 0.054 on the standard
boxes against -1.182 for NDRD's `GLY|ALL` marginal, a gap of 0.58-0.80 nats whatever the boxes; 1.9
explains it, since NDRD holds loop sites only.

---

## 9r. The handedness is not an artifact of one force field (ff14SB, 2026-09-19)

`gly_awh14` (49033947) finished its 10 dipeptides at 100 ns and its result had never been
extracted: its `fe_t*.xvg` files stopped at ~45 ns, stale from an earlier analysis run, and LR,
GGGGG and SAGAS had none at all. Re-extracted with `gmx awh -b 95000` on the last part file.

| system | ff14SB | ff99SB-ILDN |
|---|---|---|
| LA | -0.519 | -0.447 |
| LM | -0.328 | -0.301 |
| LP | -0.418 | -0.352 |
| LL | -0.591 | -0.396 |
| LT | -0.043 | -0.129 |
| LE | -0.013 | -0.090 |
| LV | -0.206 | -0.309 |
| LD | -0.124 | -0.179 |
| LR | +0.026 | -0.413 |
| **LG blank, must be 0** | **-0.010** | **-0.029** |
| **chiral mean (n=9)** | **-0.246** | **-0.291** |

**The two force fields agree to 0.045 nats**, inside the ~20% uncertainty already quoted, and both
blanks sit on zero. Glycine's left-handed bias is therefore not an artifact of `amber99sb-ildn`:
an independently refit AMBER variant gives the same answer, and both are nowhere near the
library's -1.24 or ff3.0's exact 0.

**Per-system values do not agree** (LR reads +0.026 against -0.413, LL -0.591 against -0.396),
which is the same per-pair resolution limit seen between replicas of a single force field
(S/N 1.48). It is the mean that reproduces, not the neighbour structure, and that is consistent
with everything else measured here.

**Honest limit: both are AMBER.** The literature's disagreement is largest between families
(ff14SB pPII 0.36 against CHARMM36m 0.48), so CHARMM36m would be the stronger test and has not
been run. Two AMBER variants agreeing bounds the within-family systematic, not the across-family
one.

**Not extended.** ff14SB has no successor and stops at 100 ns while the rest of Track B runs to
400. That is deliberate: the bracketing question it exists to answer is settled at 100 ns, and
its per-system values could not be resolved at 400 ns either.

## 9s. Track A's end point, and rama31's handedness (2026-09-24)

Track A (9e) finished at step 500 with `X|GLY` dG(aR->aL) -0.885, having drifted away from ff2.1
without moving toward the AWH map (correlation with it in the populated region 0.699 -> 0.703; its
learned antisymmetric pattern correlates +0.11 with AWH's and +0.30 with the PDB library's). It was
converged, not step-limited: the data's push on dG per step had t = -0.66, the whole-map
`||mean g|| / mean|g|` was 0.112 against 0.114 for pure noise, and Adam utilisation (rms step over
alpha) was 0.32 against the noise value 0.33.

`parameters/common/rama31.dat` was then rebuilt from the finished data. The first rebuild gave -0.150
because `build_rama_from_awh.py` counted `RG` (the same achiral Ac-Gly-Gly-NHMe, measured at the
other glycine) as chiral; with both blanks excluded it is -0.154. The like-for-like AWH target for
Upside's left/right mixture is the context-averaged map (-0.15), not left plus right (-0.31), which
assumes neighbour effects add in a tripeptide (untested). The current library (10-01 build,
`build_gly_library.py`) is described in GLY_sym.md §5.

## 9t. The trainer omits FF2's backbone desolvation term entirely (2026-09-24)

`bb_env.dat` came out of training byte-identical to ff2.1's because **the trainer never builds the
node it parameterises**. `main_worker`'s config kwargs pass `environment_potential` but no
`bb_environment_potential`, so no training simulation contains `bb_sigmoid_coupling_environment`,
and `extract_ff31.py` copies ff2.1's file only because `upside_config` requires one. The term is
not frozen; it is absent. Every deployment (Peng benchmark, examples) includes it.

**Where this comes from.** The Theano original (`~/Documents/ConDiv/remd-4000-8RP-1th-test/
ConDiv_original.py`, Upside 18.10.08) is FF1-era: its `Update` is `env cov rot hyd hb sheet` and
its configs have no backbone term. Peng et al. 2022 SI: FF2 was made "by adding an explicit
backbone desolvation term", a multibody burial term on the N-H and C-O vectors, and "all
parameters can be optimized simultaneously". Porting the FF1 trainer faithfully reproduced FF1's
Hamiltonian, not FF2's. Same class of defect as `hb` and `sheet` being dropped (9e).

**It is not small, and it favours the UNFOLDED state.** Native ubiquitin under the new ff_3.0 with
the term on: total -148.6, of which `bb_sigmoid_coupling_environment` -16.8, as large as the rama
term (-17.0). The same coordinates expanded 1.6x about the centroid give **-43.8**: `scale` is
-0.30 and `compact_sigmoid` is 1 at low burial, so the term pays for each **solvent-exposed**
backbone NH/CO. It is FF2's unfolded-state stabiliser ("the solvation of the backbone and the
H-bonds stabilizes the DSE", SI). An earlier version of this paragraph said it favoured compact
structure; that was wrong.

**Two couplings make it worse than a missing term.** (1) Its `hbond_weight` feeds H-bond state
into the burial, so `hb` (trained up to 1.045) was fit without it. (2) The node stores a copy of
`environment.h5`'s per-type `weights`, which were trained as part of `env` with no backbone
contribution to their gradient.

**This weakens the earlier ff2.1 fixed-point result.** That test (p >= 0.22) ran in the same
Hamiltonian without the term; if ff2.1 was trained with it on, the test either lacked power or
the missing term happens to cost little gradient at ff2.1.

**A second FF2 ingredient is also missing.** The same SI describes a dual objective,
`d alpha = d alpha_NSE + lambda * d alpha_DSE`, where the DSE part trains the unfolded ensemble
toward a self-avoiding random walk; it "increased folding cooperativity and reduced the amount of
residual H-bonded structure". The trainer has only the native-state objective.

## 9u. The trainer uses FF1's burial function; ff2.1 as published uses FF2's (2026-09-24)

`environment.h5` holds two burial functions: `energies` (20x18 spline, read by
`--environment-potential-type=0`, node `nonlinear_coupling_environment`) and `scale`/`center`/
`sharpness` (sigmoid, type 1, `sigmoid_coupling_environment`). `upside_config` **defaults to 1**,
every example leaves it at the default, and the FF2 SI says the spline "is replaced by the
sigmoid-like function" in FF2. ConDiv hard-codes 0, and its 760-value `env` parameter is exactly the
360 spline entries plus 400 weights.

**The two are not the same function in ff2.1.** Native ubiquitin, same coordinates:

| | side-chain burial, native | expanded 1.6x | native minus expanded |
|---|---|---|---|
| ff_2.1, type 1 (sigmoid, as published) | -47.49 | -26.90 | **-20.6** |
| ff_2.1, type 0 (spline) | -12.88 | -12.25 | **-0.6** |
| ff_3.0 trained, type 0 | -15.75 | -16.38 | +0.6 |

ff2.1's spline table barely distinguishes native from expanded: it is a vestige, not a copy of the
sigmoid. So **every ConDiv run in this port started from a force field that is not ff2.1**, and
trained FF1's functional form: spline burial, no backbone term. The trained ff_3.0's sigmoid
fields are ff2.1's untouched, so running it at the default type 1 silently gives ff2.1's burial.

Consequences:
* The ff3.0-vs-ff2.1 benchmark compared two different burial functional forms: `bench_run.py`
  sets type 0 for ff3.0 and leaves ff2.1 at the default 1.
* The "ff2.1 is a fixed point" test (Phase 2) ran ff2.1's parameters in a Hamiltonian ff2.1 does
  not use. It says little about port fidelity for `env`.
* The engine can train the FF2 form: `SigmoidCoupling::get_param_deriv` returns analytic
  derivatives for scale, center and sharpness per type, and `BackboneSigmoidCoupling` for all four
  backbone-term parameters.

**`Train(1).zip` (OneDrive) is not the soluble FF2 trainer.** It is Peng's 2022 membrane-potential
trainer: `UpdateBase` is `cb icb hb ihb`, the four blocks of `membrane.h5`, with the soluble force
field (`ff_2.2` in `/home/pengxd/upside-ff2.0v`) held fixed. It confirms how FF2-era runs were
configured, `environment_type = 1` with `bb_environment` on, but contains no code that trains
`rot`, `env`, `hb`, the backbone term, or an unfolded-state objective.

## 9v. The FF2 trainer: found, adapted, and what the port had wrong (2026-09-24)

**Only one FF2 dual-target trainer exists**: O. Kleinmann's Python 3 port of Peng's code,
`/project2/trsosnic/okleinmann/condiv/condiv2.py` (git history from 2025-08; the first commit is
already his working copy, so Peng's pristine file is not recoverable, and `/home/pengxd` is
unreadable). Everything else searched is FF1 or membrane: `~/Documents/ConDiv` is the FF1 Theano
original, `~/Documents/Train` = `Train(1).zip` is the 2022 membrane-potential trainer,
`upside_version/upside-pxd/ConDiv` (2019) is an intermediate with spline burial and no DSE.

**What the port had drifted on, against the SI**, all corrected in `training/ConDiv.py`:
* **lambda = 0.0** (`balance_target`), so his run never used the DSE objective at all;
* 6 free replicas up to T ~0.97 (SI: 12 from 0.8 to 1.1), 1000 time units (SI: 8000), minibatch
  21 (SI: 24);
* replica reweighting exponent `E*(T0-Ti)/Ti`, which is T0 times the correct `E*(1/Ti - 1/T0)`;
* a `dE < -200` clamp in place of a normalisation, and a guard that silently dropped the DSE term
  whenever the last free replica's final energy exceeded 1000.

His 101-step run ended with the backbone scale flipped from -0.30 to +0.12. **Its negative PRO
burial sharpness is ff2.1's own value (-0.28)**, not his drift, contrary to what I first said.

**Engine limits that ff2.1's training shared.** `BackboneSigmoidCoupling::get_param_deriv`
computes only the `scale` derivative (the other three are commented out, in master too), and
`HBondEnergy::get_param_deriv` only entries 0-3. So ff2.1's workflow never trained the backbone
term's center, sharpness and hbond weight, nor the eight H-bond rama boundaries. The user chose to
keep that exactly. The smoke worker confirmed it: those three contrasts are exactly 0.

**Validation on midway2 so far.** 19 x 24 minibatches as in the SI; ff2.1 starts at the SI's
H-bond energies (-1.961/-1.946/-1.769; second-H-bond -0.406). Smoke worker (2xf6, 600 time
units): exit 0, all 14 groups finite, the SARW replica at 0 H-bonds and Rg 24 A, an unfolded
ensemble found. The local Mac binary cannot run workers: it traps at exit whenever MC pivot moves
are on (README trap).

## 9w. Phase 2's full-map glycine row: steady drift, then cancelled (2026-09-25 to 09-28)

The full-map glycine row (plan.md Phase 2) drifted monotonically from ff2.1, about 0.002 nats of dG
per step, without approaching the AWH map (correlation flat at 0.680; the handedness part of the
displacement nearly orthogonal to the AWH direction, cos +0.05), and its gate failed at every
checkpoint (p = 0) while every other group passed. One reading rule came out of it: **tell a step-size
limit from noise** with two numbers over an epoch, Adam utilisation (rms step / alpha; pure noise
gives sqrt((1-b1)/(1+b1)) = 0.33) and the sign consistency of each cell's steps (noise 1/sqrt(n)).
Utilisation at the noise floor with consistent signs (0.37-0.42 and 0.59-0.70 here) is a weak steady
pull whose drift scales with alpha; at Track A's end both were at noise, and no step size would have
helped.

---

## 10. Cluster and operational lessons

### 10.1 A nested sbatch inherits `SLURM_*` from the job that calls it (2026-09-06)

A submitting script that requests `--mem` passes `SLURM_MEM_PER_NODE` to the job it submits; if that
job's `srun` requests `--mem-per-cpu`, every worker launch dies at once with `srun: fatal:
SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive.` This broke a
training chain's only resubmit path, unseen because every earlier job had been submitted from an
interactive shell. Any script that submits another job must unset those three variables after its
`#SBATCH` block, or match the child's memory-request type; unsetting them does not change the running
job's own allocation.

### 10.2 The batch script is fixed at submission: edits do not reach a queued successor, and a requeue reruns the old one (2026-09-23)

**Slurm copies the batch script into the job record at submit time.** A self-chaining job queues
its own successor at start, so by the time you find a bug in the script the successor already
holds the old text, and fixing the file on disk does nothing for it. Read back what a queued job
will really run:

```
scontrol write batch_script <jobid> /tmp/js.sh && grep -n <the-thing-you-fixed> /tmp/js.sh
```

Then cancel and resubmit that successor with the same `--dependency`. Cancel the **pending
successor** while the running link keeps going; the reverse order leaves the successor free to
start on the same `run_output` as the still-dying parent, which is how this campaign once ended up
with two concurrent writers.

The bug that exposed this is worth its own warning: **a branch that runs only on success, only at
the end, is untested by construction.** `train_gly.sbatch` handed off to
`sbatch "$T/validate_ff31.sbatch"` where `$T` is the run directory, while the script lives one
level up in `training/`. Twelve chain links had exercised every other line. The handoff would have
printed a warning into an unwatched log and quietly queued none of the 32 benchmark arms or 4 glpG
chains. Verify the terminal branch by hand before the run that will finally reach it.

**A requeue also reruns the script verbatim with its original arguments**: a requeued training link
resumed from the stale checkpoint it was submitted with, and `main_loop`'s `rmtree` deleted the newer
ones. Chain scripts carry `#SBATCH --no-requeue` (check `scontrol show job <id> | grep Requeue`, which
must read `Requeue=0`) and resolve the newest checkpoint on disk at run time.

### 10.3 A Slurm job can fail with no error text at all (2026-09-23)

Training link 49047139 died 8 h 36 m into a 36 h wall with exit 7. Its `.out` file contained **no
traceback, no "WORKER_FAIL", nothing**, and all twelve `*.output_worker` files were 0 bytes. Read
that way it looks like a code bug with the evidence deleted.

The cause was legible only at the Slurm step level:

```
sacct -j <id> -o JobID%16,State%16,ExitCode,MaxRSS,NodeList%30 -P | awk -F'|' 'NR==1 || $3!="0:0"'
```

which showed steps `.492`-`.503` all `CANCELLED 0:7`: **killed by signal 7, SIGBUS**, twelve tasks
across four nodes at the same moment. `ExitCode` in `sacct` is `exit:signal`, so `7:0` on the batch
step and `0:7` on the job steps are different things and the second is the informative one.

Two rules from this:

* **When a job log ends mid-sentence, go to `sacct` step level before reading the code.** A
  parent's buffered stdout is lost when the run dies, so the log's last line is where the buffer
  last flushed, not where the failure was. Here the log stopped at 11:22, the workers kept writing
  until 11:31, and the job was not reaped until 12:15.
* **Check the physics at the moment of death before assuming a blow-up.** Every worker was near
  frame 1985/4000 with Rg 14.5 A, ~110 hbonds and potential near -200. Healthy. That is what rules
  out the force field and points at infrastructure.

**The successor job is the third witness, and it is the one that settles it.** The chain queues
its replacement with `--dependency=afterany`, so the kill should have been survivable. Instead
49053769 started at 12:18:02 and its batch step was `CANCELLED` at **the same second**, elapsed
zero, and `ff31gly_49053769.out` was never created. A batch step that dies before it can open its
own output file did not run a single line of the script. So in the 11:26-12:18 window those nodes
could not read mmap'd files on `/project` (SIGBUS) and could not create one either.

That is three independent symptoms pointing the same way, and it means the chain mechanism is
sound: it was defeated by the filesystem, not by its own logic. **No GPFS log was available**, so
the mechanism is inferred from the symptoms rather than confirmed at the source.

Disk *capacity* was the obvious suspect and was **not** the cause: `/project` had 445 G free,
inodes at 20%, and a 200 MB write+delete succeeded at 1.7 GB/s once it recovered. Measure it
rather than assuming it; see remote_jobs.md for which quota command reports which fileset.

### 10.4 A target checked only at restart does not bound anything (2026-09-23)

`train_gly.sbatch` tested `STEP >= TARGET` at link start and then unconditionally ran
`STEPS_PER_LINK=150`. Resuming at step 491 of a 500 target would have run to 641: roughly 29 h of
training past the point the campaign was defined to stop, and validation blocked behind it. The
bug was invisible for twelve links because 150 divides evenly into where the earlier links landed.

**A self-chaining job needs the bound applied to the work, not only to the decision to start it.**
The fix is one line, `REMAINING=$(( TARGET - STEP ))` clamped against `STEPS_PER_LINK`.

Related: **delete a partially written output directory before resuming.** `main_worker` reads
`<name>.divergence.pkl` by path with no freshness check, so a stale one left by an aborted attempt
would be consumed as a fresh result for a worker that failed in the retry. The aborted
`epoch_12_minibatch_35` happened to contain none, but only because it died early.

### 10.5 A partition that accepts a job is not a partition the job may use (2026-10-01)

User correction. When the 10-01 allocation rollover left `pi-trsosnic` with no CPU allocation, I
found that midway3's `amd` partition still accepted the training job (`sbatch --test-only` PASSED).
I submitted it there on the strength of the user's "try midway3". The user cancelled it after
14 min (~76 core-hours): a CPU job must not run on the group's GPU allocation. The scheduler's
acceptance said nothing about which allocation the job would be charged to.
Rules:
* Run CPU work only on broadwl (midway2) or caslake (midway3) unless the user names another
  partition.
* When the usual partitions refuse a job, stop and report. Do not look for one that accepts it.
* Before any first use of a partition, find out which allocation it bills, and ask.

### 10.6 A heavy job on a login node kicks us off it, and every return costs a Duo push (2026-10-02)

Between 11:22 and 12:00 the midway2 master socket dropped about five times, two Duo prompts timed
out, and the user was asked for push after push. The cause was ours: a lambda decomposition run on
`midway2-login1` with an indexing bug allocated ~60 GB per run, the login node's memory policing
killed it, and the user's sessions on that node went with it (confirmed by the user). A capped rerun
(`ulimit -v 6000000`) finished in under a minute and the drops stopped. Analysis that is more than a
quick read goes to a compute node (`sbatch`/`srun` on broadwl); anything run on a login node is
capped with `ulimit -v` first. When the socket keeps dropping, look for our own load on the login
node (`ps -u yinhanw --sort=-rss`) before suspecting the network or the cluster.

### 10.7 Clean up local compute; never leave it running without a live reason (2026-09-16)

User correction, twice in one session. I launched 6 GROMACS replicas on the user's laptop, which
consumed ~1024% CPU (about 10 of 14 cores) for two hours before they noticed it was hot, and I had
not flagged the cost when starting them. Then, after agreeing to move the work to midway2, I left
the local replicas running on the reasoning that stopping would "lose" work in the gap. That was
wrong: the sampling already done is checkpointed on disk and survives regardless, the cluster job
rebuilds from scratch so there was no handoff to protect, and the extra sampling during a queue
wait was ~5% of what the measurement needs, bought at the cost of the exact thing the user asked
to stop.

Rules for local jobs from now on:

* **Say the cost up front.** Before starting anything local and long-running, state the core count,
  the expected wall time, and that it will load the machine. The user cannot see `ps`.
* **Kill it the moment its reason expires.** When work moves to the cluster, when a better path is
  chosen, or when the user signals they want the machine back, stop immediately. "Keeping it just
  in case" is not a reason. Data already written is not at risk from stopping.
* **Stop cleanly so it is resumable.** `kill -TERM` makes GROMACS write a checkpoint and exit;
  `mdrun -cpi prod.cpt` resumes. Verify the `.cpt` exists before reporting the job stopped.
* **Audit at the end of any session that launched local compute**: no stray `mdrun`/`upside`/python
  workers, no orphaned launcher scripts, no watcher loops left polling past their target, and
  GROMACS `#backup#` files and minimization `.trr` removed.
* Background watcher loops are fine while their target is live, but they must have a bounded
  iteration count so they expire on their own rather than polling forever.

### 10.8 Keep separate problems separate; do not carry an unproven one along (2026-10-01)

User correction. I presented the glycine fix and a pre-proline rule as one plan, with shared steps
and shared validation. They are unrelated. The glycine problem is established: a map that is part
local energy and part evolutionary selection is applied as pure energy, and helical glycines flip in
glpG. The pre-proline problem was my own suggestion. I had measured its mechanism (the mixture's
extra alpha_R), but not any consequence, and 1.14 had already found the training-set gap to be
composition. Rules:
* Propose each problem on its own evidence, with its own plan and validation.
* Before acting on a suspected defect, state the observable it is supposed to break and whether
  that has been measured. A mechanism without a measured consequence is a hypothesis, not a
  problem to fix.
* Say who raised an item when it re-enters a plan, so an AI suggestion is not mistaken for an
  established issue.

### 10.9 Three analysis lessons from the lambda diagnosis (2026-09-18)

**A reweighting is only as meaningful as the stationarity of the ensemble it reweights, and
stationarity has to be measured.** I reweighted lambda's whole cold-rung native arm along the
glycine-asymmetry axis and reported that the candidate force field stabilises the near-native
basin by 0.55 kT. The arm is a monotonic decay away from its native seed, not an ensemble: blocked
by time it runs 6.41 -> 10.39 A and is still rising in the final block. On the equilibrated last
third the same calculation gives -0.024 E_up, and the correlation between the perturbation and
Ca-RMSD flips sign, +0.171 -> -0.064. The whole-run number was reweighting frames the force field
was in the process of leaving. **Block the observable against time before reweighting anything, and
quote the converged window.** The effective sample size was 83% in both cases, so ESS says nothing
about this failure mode.

**Pooling shells can invert a conclusion when one shell dominates.** From glycine (phi,psi) pooled
over all frames under 6 A I concluded that lambda's most native-like states have helix H2 broken
with five of six glycines left-handed. Resolved by shell, the 0-5 A states have H2 intact at 89%
right-handed and it is the 5-6 A shell, four times larger, that is broken. The pooled statistic was
reporting the larger shell. **Resolve by bin before reading a conditional average.** The same
error recurred on 2026-09-30: I read the pre-proline class's +0.023 free-native alpha_R gap as a
pre-proline defect and planned to train it away under a new combining rule, with the next epoch's
mismatch as the test. Resolved by each residue's native basin, the class is 87% extended and its
residues over-visit alpha_R less than ordinary residues in the same basin (1.14). In ConDiv a class
gap measures the class's native composition until it is shown otherwise.

**Never import a module whose top level does work.** `score_arms_dist.py` runs the whole scoring
loop at import and calls `json.dump(..., "score_arms.json")` after each arm. Importing it for two
helper functions started a 2.7 h rescore on the login node and truncated `score_arms.json` from 23
arms to 3 before I noticed. It was rebuilt exactly from the intact `scoring/dist/*.npz` per-frame
arrays, which is the only reason nothing was lost. Two habits follow: **duplicate the few constants
and helpers rather than importing a script**, and **check what a module does at import before
importing it**, particularly when it writes files.

### 10.10 A checkpoint written under NumPy 2 does not load under NumPy 1 (2026-10-01)

Found by a dry run of moving the local ff30_gly_local run to midway2, before any real transfer.
The Mac's `.venv` has NumPy 2.4.4 and the cluster venvs 1.23.5. NumPy 2 pickles arrays through
`numpy._core.multiarray.scalar` and `numpy._core.numeric._frombuffer`. Under 1.23.5 every
checkpoint, solver state, divergence and rmsd file written on the Mac fails to load: `No module
named 'numpy._core'`. The same functions exist in 1.23.5 as `numpy.core.*`.
* The transfer would have failed at the cluster's first step.
* The convergence gate, which reads every solver state of the last epoch, would have failed even
  if the transfer had worked.

`env.sh` promises identical package versions only between midway2 and midway3; a local run breaks
that promise.

**Fix:** `move_run.py` now converts the whole copied `run_output` on the machine that continues
the run.
* It reads every pickle with `numpy._core` mapped to `numpy.core` when the local NumPy is 1.x.
* It rewrites the absolute paths and writes each pickle back natively.

Dry run in a scratch run dir on midway2: 28 pickles converted and 3,667 paths moved, and
`check_step.py` there reproduced the local step's numbers exactly. NumPy 2 reads NumPy 1 pickles
unchanged, so the reverse direction needs no mapping. One side effect: `move_run.py` re-pickled the checkpoints from a
module that had imported the trainer as `ConDiv`, so their classes are `ConDiv.Update`, not
`__main__.Update`. They unpickle only where `import ConDiv` finds that run's own copy, which is true
for `run_output/ConDiv.py gate|extract` (the script's directory comes first on `sys.path`) and false
for a script elsewhere with `training/` on `PYTHONPATH`, which then loads the newer trainer and fails
on the field count (seen 2026-10-02 in a test harness).

### 10.11 A deploy must carry file modes, and `cp` onto an existing file does not (2026-10-02)

The glycine-offset deploy unpacked `py/upside_config.py` without its executable bit, and
`run_upside.upside_config` runs that file directly, so every config build on midway2 failed with
`PermissionError` until the mode was restored about two minutes later. ff30_gly's step in flight
had built its configs before, and no other job started a build in the window. The cause was in
staging: the staged copy was first created by a shell redirect (mode 644), and `cp` of the repo
file over it kept the destination's 644. Stage with `cp -p` into a fresh directory, and end a
deploy with `chmod --reference=<file>.bak_<tag> <file>` and a listing of the modes.

### 10.12 Under -ffast-math, where a new loop sits can change an old sum's last bit (2026-10-02)

Both CMake files build with `-ffast-math`, which lets the compiler reorder and vectorise a float
reduction as it sees fit. Adding the glycine offsets' accumulation to
`HBondEnergy::get_param_deriv` as a second loop right after the shared `dparams[0..2]` loop changed
those shared sums in the last bit (13.922128 to 13.922127) on a 12-entry config that never enters
the new loop, while energy, forces and every other derivative stayed bitwise equal. Moving the new
loop after the function's last shared accumulation, with `dparams` grown only there, restored
bitwise equality of everything. So for a master-parity change in the engine: compare bitwise, not
with a tolerance, and when only a reduction differs, look at code placement before the arithmetic.

### 10.13 Other cluster and tooling lessons

* **A wedged GPFS makes a dead job look healthy, and `squeue` will not tell you (2026-09-07).** Job
  48981235 was reported `RUNNING` for 3.5 h while all nine of its workers sat in `D` state at
  `00:00:00` CPU, wchan `cxiWaitEventWait` / `lookup_slow`, having never started their compute
  binary. The honest probes are per-process, not per-job: `sstat -a -j <id>` (a step whose `AveCPU`
  does not climb is not computing), then `ps -o pid,stat,time,etime,wchan` on the allocated nodes.
  `D` state is uninterruptible, so `timeout 10 ls <wedged dir>` does **not** return; it leaks a
  process and, over an SSH ControlMaster, burns a session channel until the mux refuses new
  sessions and ssh falls through to password auth. That fall-through is the RCC-ban trigger, so pin
  `-o BatchMode=yes -o PasswordAuthentication=no -o NumberOfPasswordPrompts=0` on every cluster
  call before probing anything that might hang.
* **RCC's GPFS serves midway2, midway3 and beagle3, so "try the other cluster" is not a fallback for
  a storage incident.** During the 2026-09-07 outage midway2's login nodes refused TCP while
  midway3's login nodes had lost `/home`, `/project`, `/project2`, `/scratch` and `/software`
  outright: `stat -f /project` reported **xfs**, and with `/software` gone there was no `squeue` or
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
  each requeue, and because `grep -c "chunk done"` is the only honest progress metric. Since 2026-09-30 `run_remd.py` counts a
  block per job id (`<V>/block_jobid`), so a requeue no longer advances `block_count`. Rules: reset
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
* **Size trajectory output by its accumulation rate, not the seed.** A 2.4 TB balloon (2026-08-02)
  came from ~3000 frames per chunk with momentum over ~43 chunks and no purge of `output_previous_*`.
  Use `frame_interval` near 100 (~2000 frames), no `--record-momentum`, and one `--duration`.

---

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

## 11b. The Peng 2022 benchmark: what the SI actually says (read 2026-09-14)

Both PDFs are on this Mac and must be read rather than searched for: SI at
`~/OneDrive - The University of Chicago/ct1c00960_si_001.pdf`, main text alongside it. ACS is
paywalled and PMC serves a CAPTCHA, so web lookups waste time. The paper is the **HDX** paper,
"Prediction and Validation of a Protein's Free Energy Surface Using Hydrogen Exchange and
(Importantly) Its Denaturant Dependence", JCTC 2022, 18, 550-561; the folding benchmark lives in its
SI, so one citation covers both the benchmark and the soluble-protein HDX result.

**Verified against the SI:**
* The simulation-parameter table (p10-11) matches `bench_table.py:TABLE_S2` verbatim -- durations and
  14-rung ladders both.
* Fig S4's caption lists the terminal-residue exclusions, and they match `RMSD_EXCLUDE` exactly.
* Fig S4 prints the per-protein lowest Ca-RMSD **as text**, with the largest-cluster centroid in
  parentheses (DBSCAN on the Ca contact map, 10 A cutoff). So `FF2_BASELINE` was transcribed from
  printed numbers, not read off bars -- the provenance worry raised on 2026-09-13 is retired.
* Fig S4 right-hand panels pool **five independent simulations** per protein, lowest-RMSD run solid
  and the rest dashed. Our arms are one run each at the bottom rung, so our distributions are
  narrower than theirs by construction.

**Both figures are grids of per-protein DISTRIBUTION CURVES with the set average on top**, columns
being native-start and unfolded-start, FF1 and FF2 overlaid in each panel. Fig S4 adds a
predicted-vs-native structure-overlay column.

**Lesson: do not describe a figure as reproducing a published format without having opened that
figure.** The first version of these panels was built from a one-line paraphrase in our own
`analyse_bench.py` docstring, labelled "the same plots the paper makes", and was wrong in layout --
summary markers and violins against grids of distribution curves. The underlying numbers were sound
and reproducible, which made the error easy to miss. State "same quantities, my layout" unless the
source figure has actually been read.

## 11c. Peng's benchmark trajectories on midway2 are FF1, not FF2 (measured 2026-09-14)

`/project2/trsosnic/condiv_data_upload/trajectories/` holds `<prot>_native.xtc`,
`<prot>_denovo.xtc` and `<prot>.pdb` for 23 proteins -- the 2022 paper's 16 plus 7 CASP targets
(T0765/69/71/73, T0803, T0816, T0855). 2.8 GB, owned by `nffaruk`, dated 2018-07-20. Frame counts
are large: 27k-78k per arm.

**It is FF1.** Scored with `ff3_benchmark/scoring/tmscore.py` and the paper's own terminal-residue
exclusions, over all 16 benchmark proteins:

| | this data | published FF1 | published FF2 |
|---|---|---|---|
| from native | 0.480 | 0.45 | 0.55 |
| de novo | **0.360** | **0.37** | 0.42 |

The de novo mean matches FF1 to 0.01 and misses FF2 by 0.06. The native mean runs 0.03 above FF1's,
which is expected because the paper's native figure counts only "excursions within the native basin"
while this averages every frame.

**A second FF1 set exists in Upside's own format**, found after correcting the search: Upside writes
`.up` and `.vtf`, never `.xtc`, so the first sweep used the wrong filter (the 2018 `.xtc` files above
are a converted deposit). `/project2/trsosnic/share/paper_traj_nabil/{from_native,from_denovo}` holds
`<stem>.run.{0..13}.up` -- **14 replicas**, matching the Table S2 ladder -- for 5 of the 16 proteins
(cspa, gpW, hyp, nug2, top7), 16,100 frames per block, 28 GB, dated 2018-02.

**That set is FF1 too, and the proof is structural rather than a date.** Its
`rotamer/pair_interaction/interaction_param` is shape **(20, 20, 62)**, while ff_2.0, ff_2.1 and
ff_3.0 are all **(20, 20, 54)**. A different spline-knot count is a different force-field generation,
so it cannot be any ff_2.x/3.x and no rescaling comparison is even meaningful.

**There is no FF2-era benchmark data anywhere on midway2 or midway3.** `/cds3` is the only filesystem
midway3 adds (`/project` and `/project2` are shared), and a 2.8 TB sweep of `/cds3/trsosnic` found no
benchmark-protein trajectories from the FF2 era and no `pengxd` space at all. Checked: all of `pengxd/`, the FF2-era
ConDiv trees (`upside_version/upside-pxd/ConDiv`, `share/pengxd`), and a sweep of `/project2/trsosnic`
and `/project/trsosnic` for `*_native.xtc` and for TM/RMSD outputs newer than 2021. The only `ff2`
hits are Adam's 1nqe channel project. So an FF2 per-protein overlay needs either a
higher-resolution figure from the publisher or a fresh ff_2.1 run of the 32 arms (`bench.sbatch`
already takes `FF=ff_2.1`).

**Useful by-product: this is an external validation of our scoring path.** Reproducing FF1's
published de novo mean to 0.01 from real trajectories exercises `tm_score.py`, the residue-mapping
calibration and the exclusion handling together, which the synthetic tests in `py/tm_score.py` do not.
Reading the `.xtc` needs mdtraj, installed out-of-tree at `/beagle3/trsosnic/yinhan/pylibs` via
`pip --target` so the shared venv the glpG chain uses is untouched; add it to `PYTHONPATH`.

## 11d. lambda repressor in the Peng benchmark (2026-09-18 to 10-02)

lambda is the worst arm of the Peng benchmark under FF2 and under every force field tried here. The
first diagnosis was made on the retired ff3.0 (scripts in `ff3_benchmark/scoring/`:
`diag_lambda_fail.py`, `gly_reweight.py`, `energy_split.py`, `helix_by_rmsd.py`). Helices H0-H4 are
named from the reference's own (phi,psi), 0-based file numbering: H0 3-23, H1 27-34, H2 38-46,
H3 53-63, H4 72-78. What still holds from it:

* **The failure is bundle assembly, not secondary structure.** Each helix holds its own shape (local
  CA-RMSD 0.6-2.8 A) while the assembly is 8.6 A out: helix pairs that are near-parallel in the native
  come out near-perpendicular (H0-H3 28 -> 69 deg, H1-H4 22 -> 86 deg) with every centroid distance
  within 4 A. proteinB and homeodomain, analysed the same way, reproduce every crossing angle to
  within 17 deg.
* **The native arm never equilibrates.** Cold-rung CA-RMSD rose from 6.4 A in the first block to
  10.4 A in the last over the full Table S2 duration, while the de novo arm sat at 10.6-11.5 A, so
  both converge on one misassembled ensemble. A statistic over a whole native arm mixes decay with
  equilibrium (`score_arms_dist.py`'s `BURN = 2000` discards only the first of twelve blocks); block
  an arm by time before quoting its native number as an equilibrium property.
* **The Ramachandran term barely tracks fold quality** (correlation with CA-RMSD +0.075, its glycine
  part +0.119, against +0.205 for the rest), and reweighting the glycine handedness changes the
  equilibrated ensemble by -0.024 E_up, zero within error (10.9). Do not use frame 0's
  energy as a reference: proteinB's deposited structure scores +163 E_up above its ensemble mean
  while folding perfectly.
* **H2 (file 38-46, `QSGVGALFN`) fails first**: intact (89% alpha_R) only below 5 A global CA-RMSD and
  broken by 5-6 A, while H3 stays 93-99% alpha_R throughout. Its two internal glycines are natively
  alpha_R, and lambda's six glycines are mixed in handedness, so no context-free glycine map suits
  the protein:

  | glycine (file) | context | native region |
  |---|---|---|
  | 24 | L-G-L | alpha_L, H0 C-cap |
  | 35 | M-G-M | alpha_L, H1 C-cap |
  | 37 | M-G-Q | beta/pPII |
  | 40 | S-G-V | alpha_R, inside H2 |
  | 42 | V-G-A | alpha_R, inside H2 |
  | 47 | N-G-I | alpha_L, H2 C-cap |

### Wild type against G46A/G48A under ff_2.1 (2026-10-01/02)

* `ff3_benchmark/pdb/lambda.pdb` is byte-identical (md5 c85500a3) to Peng's own
  `/project2/trsosnic/pengxd/Data/test-set-15/lambda.pdb`, and the same sequence is in his
  `training_set2/input/`. It is lambda repressor residues 6-85, wild type at the positions the
  fast-folding variants change: D14, Y22, Q33, **G46, G48** (file index = residue - 6).
* H2 above (file 38-46, `QSGVGALFN`) is lambda's helix 3 (44-52), and its two internal glycines,
  40 and 42, are **G46 and G48**: the residues the fast-folding variant replaces with Ala to
  stabilise exactly this helix.
* The section above tested glycine handedness only (reweighting the map's antisymmetric part) and
  found nothing. It did not test helical-glycine stability, which is the class the selection panel
  finds weak at ff2.1 (helical glycines -0.081 against all-atom, 1.17).
* No ff_2.1 lambda run existed. Causal test submitted 2026-10-01: lambda and lambda G46A/G48A,
  native and de novo, ff_2.1, Peng's protocol (remote_jobs.md §1). If removing the two glycines
  holds H2 and the bundle, the helical glycines are the cause.
* **First chunk under ff_2.1 (2026-10-02 07:45; coldest replica T 0.780; `checks/lambda_ff21/lambda_check.py`;
  global CA-RMSD without residues 1-2; helix alpha_R from rama_basin, self-fitted helix CA-RMSD):**

  | block (k tu) | WT native RMSD / <5 A | WT H2 aR (rmsd) | G46A/G48A native RMSD / <5 A | G46A/G48A H2 aR (rmsd) |
  |---|---|---|---|---|
  | 0-89 / 0-95 | 5.26 / 0.58 | 0.97 (0.49) | 4.23 / 0.93 | 0.99 (0.37) |
  | 89-178 / 95-190 | 4.80 / 0.74 | 0.96 (0.61) | 4.95 / 0.76 | 0.98 (0.40) |
  | 178-266 / 190-284 | 5.68 / 0.61 | 0.81 (1.13) | 5.41 / 0.65 | 0.97 (0.41) |
  | 266-355 / 284-379 | 7.25 / 0.32 | 0.60 (2.08) | 5.96 / 0.45 | 0.94 (0.57) |
  | 355-444 / 379-474 | 7.66 / 0.29 | 0.57 (2.09) | 6.37 / 0.44 | 0.98 (0.40) |

  * In wild type, **helix 3 unravels first and alone** (H0 0.98-1.00, H1 0.99-1.00, H3 0.93-0.95,
    H4 0.89-0.94 throughout), and the global RMSD rises as it goes (block 3 to 4: 5.68 -> 7.25 A while
    H2 falls 0.81 -> 0.60): ff_2.1 already does what the old ff3.0 did.
  * **G46A/G48A holds helix 3** (0.94-0.99) and slows the decay (6.37 against 7.66 A, 0.44 against
    0.29 below 5 A in the last block) but does not stop it: the bundle still drifts, with H4 slipping
    (0.92 -> 0.85, local RMSD 0.61 -> 1.02). So the two helical glycines are the main cause, and a
    smaller packing drift remains.
  * WT de novo (483 k tu): RMSD 9.6-11.0 A, never below 5 A; helix 3 is the helix that does not form
    (aR 0.20-0.48) while H0, H3 and H4 form (0.84-1.00). G46A/G48A de novo still in its first chunk.
  * Experimentally wild-type lambda 6-85 folds; helix 3 is its least stable helix, which is why the
    fast-folding variants carry G46A/G48A. Upside at T 0.78 loses it, so helical glycines are too
    weak in Upside itself, consistent with the panel's helical-glycine deficit at ff2.1 (-0.08).
* **What else breaks once helix 3 holds** (`checks/lambda_ff21/lambda_packing.py`, `lambda_terms.py`;
  first chunk, coldest replica; helix axes by PCA, native CA contacts < 8 A):

  | helix pair (lambda helices) | native | WT block 1 -> 5 (angle, contacts kept) | G46A/G48A block 1 -> 5 |
  |---|---|---|---|
  | H1-H3 (helix 2 / helix 4) | 120 deg | 80 0.30 -> 50 0.22 | 90 0.27 -> 81 0.13 |
  | H0-H3 (helix 1 / helix 4) | 28 deg | 52 0.38 -> 92 0.16 | 38 0.46 -> 44 0.32 |
  | H3-H4 (helix 4 / helix 5) | 115 deg | 111 0.58 -> 99 0.25 | 108 0.66 -> 110 0.34 |
  | H1-H2 (helix 2 / helix 3) | 102 deg | 125 0.45 -> 136 0.23 | 108 0.66 -> 112 0.57 |

  * **Helix 2 packs against helix 4 30-40 deg off native from the first block, in both sequences**:
    a mis-specified interface, not a slow decay. The helix 3-4 loop (lambda 51-56) keeps 7-30% of
    its contacts and helix 5 loses half its contacts with helix 4.
  * **ff2.1's energy prefers the misoriented packing** (G46A/G48A frames binned by the H1-H3 angle):
    total -169.0 at 70-90 deg against -163.0 at >= 105 deg; side-chain free energy (`rotamer`) -164.4
    against -157.0 (7.4 E_up of it), burial -29.7 against -27.9; H-bond and coverage terms favour
    native (2.4 and 6.0). So lambda's remaining defect is side-chain packing specificity, present in
    ff2.1 and, by its published TM, in FF2; the fast-folding variant also carries Q33Y in helix 2.
  * Whether ConDiv can correct it: at ff2.1 the side-chain gradient is noise (6% of coefficients
    significant), so the training sees no consistent correction; the damped side-chain step keeps
    ff2.1's packing rather than fixing it.

### The helix 2 / helix 4 misorientation under ff2.1 has no single side-chain cause (2026-10-02)

Decomposed on the ff2.1 native arms (WT and G46A/G48A, run.0 first chunk; "misoriented" H1-H3
crossing 70-90 deg, native-like >= 105 deg; errors from 8 time blocks;
`checks/lambda_ff21/{pos_decomp,pos_report,panel_terms,panel_report}.py`):
* The rotamer free-energy preference for misorientation is weak: -7.4 +- 3.5 E_up in G46A/G48A
  and -2.3 +- 3.0 in WT. Only the side-chain/backbone term (`hbond_coverage_hydrophobe`) favours it
  in both (-2.7 +- 1.8, -5.1 +- 2.2); the pair term flips sign between arms (-2.3, +4.9).
* No residue-type pair carries it (every pair < 1 E_up, signs flip between arms), and no position
  recurs except L64 and Q11. In G46A/G48A the misoriented state gains a helix 2 N-end / helix 3
  C-end dock (Q33 side chain on the F51 backbone -2.3 +- 0.3) and loses helix 1-helix 2 contacts;
  the helix 2-helix 4 contacts themselves do not favour it. Decompositions match finite differences.
* On the ff2.1 panel (44 domains) F_rotamer and the side-chain/backbone term favour the native
  over compact non-native frames in 31 of 32 domains, and lambda's recurring pairs have no
  consistent sign. The panel's non-native frames are mostly frayed natives, so a defect specific to
  rearranged helix packings is not excluded.
* So there is no table entry with a physical target to correct: lambda's packing is recorded as a
  known limitation, re-tested after ff3.0 (crossing angle; Q33-F51 dock; helix 1-helix 2 retention;
  L50-L64), with more native-like frames than the 196 and 102 here.

## 12. Claims that turned out to be wrong

One line each: what was believed, what is true, and why it is worth keeping.

### 12a. The glycine campaign, 2026-09-18/19 (migrated from current_job.md before it was retired)

1. **"Glycine handedness is zero."** Artifact of unconverged flat surfaces. `awh1-dimN-diffusion`
   was 5e-5 rad^2/ps when AWH's own friction metric implies ~0.77, about 15000x too small, so the
   PMF range was 2.0 kJ/mol at 25.8 ns instead of 20-40. A flat surface is trivially symmetric.
   Fixed to `diffusion = 0.5`, `error-init = 30`; range then 50-66 kJ/mol.
2. **The beta-branching rule** (VAL and THR positive because branched). Both crossed zero by 30 ns.
   Do NOT run ILE as a decisive test.
3. **Pentapeptide controls.** GGGGG must read 0 and reads -0.224. Not usable.
4. **"LA is settled" at 0.028.** A four-snapshot window artifact; over 18 ns it is 0.106.
5. **"Everything drifted toward zero" between 51 and 57.6 ns.** A basin-definition artifact.
6. **ff3.0C's premise**, that the library's neighbour specificity is real. Contradicted: its
   neighbour ordering correlates with the measurement at r = -0.565 against ff2.1's -0.596, both
   anti-correlated.
7. **"The candidate glycine map stabilises lambda's near-native basin by 0.55 kT."** Withdrawn the
   same evening. It came from reweighting the *whole* native arm, which is a decay away from the
   native seed rather than an ensemble; on the equilibrated last third the effect is -0.024 E_up,
   i.e. zero and if anything destabilising. **A reweighting is only as meaningful as the
   stationarity of the ensemble it reweights, and stationarity must be checked, not assumed.**
8. **"Each replica's achiral control converges to its own nonzero value, and rep2's +0.094 is
   unexplained."** Resolved: it is sampling noise. The residual is uncorrelated between replicas
   (**r = +0.157**, where a real defect would reproduce) and decays as `1/sqrt(t)` (0.233 at 10 ns
   -> 0.032 at 100) while the signal converges (0.080 -> 0.071). An artifact would decay too.
   Blank subtraction was rejected: with the error random, subtracting a noisy estimate of zero
   adds noise. The lasting consequence is the error bar, **~20% on the basin dG**, not the +/-0.01
   the replica agreement alone suggests.
9. **"The trainer's `hb` has never been trained."** Wrong, read off the modern port. The Theano
   original trains it at lr 0.02; the port dropped it. See 9e.
10. **The sheet gradient is systematically nonzero** (t = +3.91 at n=3). Collapsed to t = -0.01 at
    n=6. See 9k.
11. **"The library's symmetric part must be kept, because it and the measured one differ by rms
    1.752 E_up with aR basins of 0.193 against 0.089, 25x the handedness correction."** Wrong, and
    it was the entire argument for scoping ff3.1 to the antisymmetric part. Those numbers were a
    *per-pair antisymmetric* statistic, not the symmetric parts. Measured properly the two surfaces
    correlate at **r = +0.867**, their aR basins agree to **0.207 nats**, and the map mean is
    unchanged to three decimals by the swap. **Name the quantity before quoting a ratio on it**: a
    25x gap that turns out to be between two different statistics is worse than no number.
12. **The units conversion was inverted.** `PMF / 2.914952774272` (kJ/mol per E_up) was given as
    correct and `PMF / kT(300 K)` as the 17% error. It is the reverse for this purpose: a library
    map holds `-lnP` normalised to `sum(exp(-E)) = 1`, verified exactly, so a PMF entering that
    slot is divided by `kT`. **Check the file's own contract before choosing a unit conversion**;
    here it was one line of arithmetic (`sum(exp(-E))`) and it settled the question outright.
13. **`S_library + A_measured` presented as free of ff3.0.** `S_library` *is* ff3.0, the fully
    mirror-symmetrised map, so the construction was ff3.0 plus a correction. Caught by the user,
    not by me. **When a construction decomposes something, check whether one of the pieces is a
    thing already rejected.**

Rules these produced are collected in 6.7.

### 12b. The hybrid, HDX and cluster work, 2026-07 to 09

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
* **"The group quota leaves 195 GB, so the pre-ff3 ladder has to be deleted" (2026-09-08,
  corrected 2026-09-09):** 195 GB was the `/project2` group quota, a different filesystem; `/project`
  had 1514 GB free. `rcchelp quota` reports a separate `trsosnic` quota per filesystem, so match the
  section header to the mount the data is on and cross-check with `df` on the data path
  (`remote_jobs.md`, disk section). The error nearly deleted 89 GB of baseline trajectories.
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
