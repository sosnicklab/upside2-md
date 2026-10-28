# Findings

Knowledge base for this branch: standing rules and the measurements behind them, and the
Ramachandran work (1), how the hybrid is put together (2), defects whose causes are established (3),
what the hybrid still gets wrong (4), HDX (5), reusable diagnostics (6), system preparation and
bilayer physics (7-8), the nanoparticle campaign (9), the trainer and glycine-map history (9c-9w),
cluster lessons (10), references (11) and claims that turned out to be wrong (12). Old update numbers
("findings 103") are kept inline where code or other files cite them. Job state lives only in
`remote_jobs.md`.

---

## 1. Standing rules, and the Ramachandran work

### 1.1 A spline table must be the published potential

The MARTINI table shipped before 2026-08-05 was not dry-MARTINI: it stored bare LJ plus bare `1/r`
Coulomb, hard-truncated at 1.2 nm. Published dry Martini runs `coulombtype = reaction-field`,
`epsilon_r = 15`, `epsilon_rf = 0` (conducting) and `vdw-modifier = Potential-shift`, so both terms
reach the cutoff smoothly. The old table stepped at `r_c` by `k*qq/r_c` = **2.65 E_up for a charged
pair, about 3.4 kT** (analytic LJ -0.057 + Coulomb -2.648 = -2.705, the stored -2.7049).

Flagged as outstanding in findings 75 and pinned down in the glpG micelle runaway (findings 79-80);
`py/martini_build_tables.py` cites those numbers. Fixed there in both nonbonded builders (the
particle-particle grids and the SC-env/BB-env `_pair_energy_and_grads`, which carried the same bare
form and analytic gradient). `scratchpad/verify_table_matches_drymartini.py` rebuilds the reference
in native kJ/mol + nm and converts:

| table | max deviation from the reference form | rows non-zero at the cutoff |
|---|---|---|
| old (`martini.h5.bak.pre-reactionfield`) | 3.95 E_up | 81 / 81 |
| regenerated | 2.2e-11 E_up (round-off) | 0 / 81 |

Neutral pairs change by a constant only (forces bit-identical); charged pairs by 3.1-3.9 E_up over
the sampled range, since the reaction field is r-dependent, so results change for anything charged:
ions, PC/PG headgroups, charged residues. Robertson et al.'s `step8_production.mdp` uses the same
contract (`reaction-field`, `epsilon_r = 15`, `epsilon_rf = 0`, `Potential-shift-verlet`,
`rvdw = rcoulomb = 1.2`, `ref_t = 303.15 K` = our 0.8647 T_up).

GROMACS `epsilon_rf = 0` means infinity, so `(eps_rf - eps_r)/(2 eps_rf + eps_r) = 1/2`, not -1;
reading the 0 literally makes the charged triples look ~100% wrong at large r. The current table
equals the analytic form to max relative error **3.6e-11 over all 23 (eps, sigma, qq) triples** at
r = 2.5-11.9 A, and the builder asserts it at build time. Do not re-litigate the table when hunting
an instability.

### 1.2 NO GUARDS

User directive, now a CLAUDE.md rule: no guard anywhere for any numerical problem. Removed from the
engine on 2026-08-05:

* `src/main.cpp`: the blow-up guard that aborted on a non-finite potential or kinetic energy, with
  its `blew_up` / `blew_up_ns` / `blew_up_round` state, loop break and FATAL block.
  `compute_logged_kinetic_energy` stays for ordinary logging.
* `src/martini_potential.cpp`: three silent-skip masks, the harmful ones. The pair loop's
  `if(isfinite(pot) && isfinite(force_mag))` dropped non-finite pair contributions (the main.cpp
  comment admitted its kinetic check existed to catch "diverging momenta even when non-finite pair
  forces are masked out of the potential"); the SC-env path had the same pattern twice.
* Kept in `update_martini_node_boxes`: `!(scale_xy > 0.f) || !(scale_z > 0.f)`, a domain
  precondition (a box length is positive). The `isfinite` clause beside it was removed.

Consequence: a divergence now propagates. NaN enters the log, is exchanged between replicas and is
written into restarts; the evidence survives instead of the run, so a bad run must be caught by
monitoring.

Three checks were left alone pending a ruling, since none masks a numerical error:
`assert_environment_solvation` (prep-time, fails a build whose belt faces vacuum), the `np_hybrid.py`
ion assertions (neutrality and salt composition of a built system), and `run_np_prod.py`'s
`health()` peptide C-N check, which ends a chain rather than propagate a torn system. The last is
closest to a guard and most needs a decision.

### 1.3 Timesteps are calibrated quantities, not free parameters

**NP-1AO6: dt = 0.001, not to be raised.** Two measured constraints:

1. *MARTINI LJ-core stability.* Velocity-Verlet needs `dt < 2/omega`. With `omega = sqrt(U''/mu)`
   from the tabulated grid curvature, `dt_max = 0.0209 @ 4.0 A, 0.00808 @ 3.5 A, 0.00430 @ 3.2 A,
   0.00278 @ 3.0 A, 0.00186 @ 2.84 A`. The closest protein-environment approach sampled is 2.84 A, so
   dt = 0.005 was 2.7x over the limit, hence the stochastic failure (16x spread in onset across six
   faces). A/B restart from the identical pre-tear frame: dt = 0.005 gave 42 broken bonds, Epot
   -7988 -> +9130, avg KE/1.5kT = 5.34; dt = 0.001 gave 0 broken bonds, Epot -7988 -> -8146, avg
   KE/1.5kT = 1.006.
2. *Backbone spring accuracy at large amplitude during unfolding.* After the findings-88 interface
   fix, dt = 0.005 still destroyed runs 1 and 4 at t ~ 250, at residues 119-122, 70+ A from the NP
   surface with no MARTINI pair near (nearest ion 8.6+ A). The linear criterion holds
   (`omega*dt = sqrt(48/0.5)*0.005 = 0.049 << 2`), but unfolding drives backbone bonds to 0.3-0.5 A
   amplitude, 2-3x the thermal `sqrt(kT/k) = 0.14 A`, where the bond force reverses sign over a
   half-period of ~0.32 t_u, and 54 force evaluations per frame interval track that too coarsely.
   Coincident CA-C and C-N spring forces put 14.4 E_up/A on C(121), injecting 12 E_up over 54 steps
   and stretching CA-C to 2.58 A. So "the springs are not the constraint, dt/dt_max = 0.035" holds
   only for small thermal oscillation, and a 50 t_u validation run cannot sample the unfolding regime
   (failure at t > 250).

**glpG hybrid: dt is locked at 0.009** (`martini_brownian.cpp:100` throws on mismatch): the friction
is tuned against it for a target lipid diffusion of 11.5 um^2/s, so changing one silently
invalidates the other. Its stability problem is handled by sub-stepping the integrator (2.4).

Mass repartitioning would buy stability (protein sites are 1 m_up against 6 for a MARTINI bead;
`dt_max ~ sqrt(mu)`) and is exact for an equilibrium observable such as an HDX free energy, but
`/input/brownian` ties `numerical_time_step 0.009`, `target_lipid_diffusion_um2_s 11.5` and
`bare_particle_friction_up 0.169` together, so it voids the friction/diffusion calibration: a physics
decision, not a fix to apply unilaterally.

### 1.4 The dry-MARTINI unit contract

Native-to-Upside conversion happens once, in the Python h5 build; the engine does no unit math.

* The main `martini_potential` node's `coulomb_k` is read but never used in compute: LJ+Coulomb is
  baked into `combined_energy_grids` by `convert_stage` (the particles group is already in E_up:
  `unique_eps_eup`, `unique_sig_ang`, `combined_energy_grids`). The old conversion attrs on the
  config node fed only the Python softening builder, which is why moving the conversion was
  bit-identical on the bilayer (-11039.306641, |diff| = 0.000000). The single baked attr is
  `coulomb_constant`.
* `sc_table`: build and read in native units; convert and rename only at the final `create_dataset`
  (`grid_ang = grid_nm*10`, `*_energy_eup = *_kj_mol/2.914952774272`). The C++ tail-subtraction
  shift stays, a physics zeroing at the cutoff applied identically before and after:
  `(native[ig]-tail)/E == eup[ig]-eup_tail`, so the float32 result is unchanged (maxabs 0.0).
* Dimensionless datasets (`angular_profile`, `rotamer_angular_profile`, `cos_theta_grid`) are not
  converted.
* The engine reads `grid_ang` / `*_energy_eup` / `coulomb_constant` with no arithmetic on them and is
  analytic-LJ/Coulomb-free: `martini_potential` evaluates only splines (`combined_spline`), the
  node-level epsilon/sigma/coulomb_k reads and the analytic PairParam coefficients are gone, and
  lj_cutoff/coul_cutoff are one cutoff attr.

### 1.5 dry-MARTINI membranes are NVT

A tensionless barostat is the wrong ensemble for these systems; the second reason settles it
generally:

1. A defect is a stress sink. Lateral tension relaxes through a pore or an under-filled patch
   instead of the lipid area, so zero tension no longer marks the intact bilayer's equilibrium and
   the barostat compresses the intact regions past it.
2. Implicit solvent has no solvent virial. Dry MARTINI carries no water, so the lateral pressure the
   barostat reads lacks the solvent contribution and is not the physical one. Dry Martini (Arnarez et
   al., JCTC 2015) is specified NVT for this reason.
3. Papers that use semiisotropic Parrinello-Rahman on these lipids do so validly because their
   systems are wet; the water is only stripped from the frames they deposit.

On our POPE/POPG tile the tensionless xy barostat condensed the bilayer to **APL 56.0 A^2 against the
reference 61.7 (-9.4%)**, tail core +1.6%, head-head +0.9%: over-compressed and thicker. The area is
therefore an input matched to an equilibrated reference of the same lipids at the same temperature,
and validation uses the area-sensitive structural observables measured at that area. Before adopting
an ensemble from a paper, check that its system has the degrees of freedom that make it valid.

### 1.6 Do not transfer settings, thresholds or analysis between simulations

NP and glpG are different simulations with different integrators and are analysed separately
(comparison table in `remote_jobs.md` §2). Conflating them cost a healthy 6 h glpG block, when a C-N
threshold borrowed from NP false-positived, and produced a wrong root-cause writeup. They also need
different detectors: in a forced NP tear the protein reached 431 broken bonds with the potential
still finite at +3e5, while glpG's blow-up goes fully NaN within one 46-step interval (6.4).

### 1.7 No hard-coded protein or system identity

Scripts under `py/` are shared infrastructure. A hard-coded id does not fail loudly; it silently
attaches one protein's metadata to another system's trajectory. Instances found and fixed:

* `martini_extract_vtf.infer_pdb_id` returned `"1rkl"` for any non-bilayer input, which mislabelled
  the 1AO6 nanoparticle VTFs; `martini_prepare_system` defaulted `--pdb-id` to `1rkl` and `--run-dir`
  to `outputs/martini_test_1rkl_hybrid` (fixed 2026-08-09).
* `martini_hdx_membrane_accessibility.py` defaulted `--tail-bead-names` to `C1,C2,C3`, the DDM tails.
  On POPE/POPG it found no tail beads and raised, which was luck: a lipid with a bead called `C1`
  would have been scored on the wrong subset silently. Tails are now detected from the trajectory by
  the MARTINI apolar naming pattern `^[CD]\d[A-Z]?$`, which covers a single-tailed detergent and a
  two-tailed lipid alike.
* The same error can arrive through a wildcard: `4.calc_D_uptake.py` and `5.analyze_D_uptake.py`
  located the experimental HXMS arrays with `glob.glob(..._d_norm_peps_*.npy)` and took
  `matches[0]`, so the requested `protein_state` was ignored whenever more than one state matched.
  This dataset has four, and it was comparing against `pd9 SUB` while `pd9` was requested. **A glob
  that can match more than one identity must be resolved by the identity, not by `[0]`.** Both now
  call `helpers/function.py:select_state_file`, which raises naming the state and the available
  files.

### 1.8 Each Ramachandran map is its own conditional distribution

User correction, 2026-09-27. Every library map `(central X, direction, neighbour Y)` is
`P(phi,psi | X, Y)` over its own set of PDB sites, and the fold-context selection that shaped it
belongs to those sites alone. `GLY|GLY` says nothing about `GLY|ALA`, and central GLY says nothing
about central ALA. A correction, a symmetry or a "shared mode" measured on one map must never be
transferred to another, and pooled maps (`ALL`, neighbour averages) must not be used to derive a
per-map correction. The error that prompted this: a Bayes argument (`P(phi,psi|aa) ~ P(aa|phi,psi)
P_bb`) was used to claim the contamination is one function shared by all residues, and GLY's
antisymmetric part was then subtracted from every row. That assumes evolution selects residues by
(phi,psi) alone; it selects by the whole site (burial, packing, turn type), so each map carries its
own selection factor. The one real coupling between maps is statistical: NDRD's hierarchical
Dirichlet process shares mixture components across neighbour maps, so rare pairs are pulled toward
their pooled map. That is estimation, not physics.

Corollary, from a user correction the same day: separating local physics from fold selection in a
map needs a second, per-pair source of information. If that source is Upside's own native-state
simulations (a per-map reference ratio, `map += ln(H_free/H_restrained)`), it is ConDiv's native
term in closed form, not a new design: ConDiv already records those per-residue populations
(`rama_native`, `rama_free`). The only per-pair sources independent of the model are physics or
experiment on that pair. Before presenting a design as new, write down the objective of the existing
pipeline and compare.

### 1.9 NDRD's site set is loops only; Upside applies the map to every residue (2026-09-27)

Ting et al. 2010: 3,038 PISCES chains (<= 1.7 A, R <= 0.25, < 50% identity, EDS density), bottom 20%
of density per residue type removed, Stride assignment, then only loop residues at least three
positions from any H or E; TCB also drops 3-10/pi helices and their neighbours. TCB is 44,112
residues of all types, about 110 per directional map (NDRD pools rare maps through its hierarchical
Dirichlet process). The 456 training proteins have 49,172 interior residues with defined (phi,psi),
a median of 107 per directional map, but only ~8,700 TCB-like residues (backbone-only DSSP
reconstruction, no 3-10 exclusion, so approximate). Same boxes throughout, glycine ln(aR/aL): NDRD
TCB -1.19, 456 TCB-like -1.31, **456 all residues -0.58**; ALA aR fraction 0.26 / 0.40 / 0.61.
Glycine's alpha_L excess is therefore mostly a property of the loop-only site selection, and a map
trained against loop sites leaves every helix and strand residue outside its objective. In ConDiv
the SARW replica zeroes H-bond, burial and rotamer-pair energies but keeps rama (`zero_for_sarw`), so
the DSE term on a map compares the unfolded ensemble with the map's own distribution: a
self-consistency condition, not an external target.

### 1.10 A basin offset must be a weight on its basin (design rules, 2026-09-28)

The basin-offset training these rules governed was removed on 2026-10-01 (plan.md);
`ConDiv.py` only records per-residue basin populations with the same basins. The rules hold for any
additive correction on a Ramachandran map:

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
entry is not enough. NDRD's sheet `GLY|GLY` holds 94% of its weight at phi < 0 (beta 0.80, pPII
0.14, beta' 0.06), which gave the middle glycine of A-G-G-G-A ln(beta/beta') = +0.14 to +0.32 while
ln(aR/aL) stayed 0.0000 (the sheet maps are empty in both helical basins). Symmetrise by averaging
probabilities, which pools each site with its mirror image, not energies: the energy mean is a
geometric mean of probabilities and would put the sheet map at beta 0.5 / pPII 0 per side instead
of 0.43 / 0.07. `build_gly_library.py` wrote one symmetric `GLY|GLY` entry to both coil and sheet.
Grid convention, verified: the engine (`rama_map_pot.cpp:66`) and `upside_config` both put node i at
-180 + 5i deg, which the rolled mirror of 6.1 assumes.

Per-pair independence was verified on the round-1 ff30_basin state (`checks/pair_independence.py`):
every (central, direction, neighbour) map was its own parameter set, and the one coupling was in the
data, since an interior residue reads a left and a right map and its populations are credited to
both. On that state the GLY|X offsets deepened alpha_R relative to NDRD in 30 of 38 maps (mean
+0.060), yet every map still favoured alpha_L, because the target itself is alpha_L-rich (1.14).

### 1.12 What 456 proteins can resolve per Ramachandran map (2026-09-28)

Measured on the ff30_basin epoch-0 data, all offsets zero (`checks/param_data_audit.py`,
`basin_schemes.py`, `basin_types_beta.py`):

* **ff2.1's own groups** hold 31,285 trained values: 30,800 side-chain spline coefficients that
  receive a derivative (of a 31,420-entry latent), environment 460, sheet 20, H-bond 3 plus dhb 1,
  backbone scale 1. Pair data: 306k CA-CA < 10 A contacts over 210 pair types, median 1,200, fewest
  TRP-TRP 84.
* **Per-pair map corrections are noise-limited.** Training sites per directional map: median 94, 76
  maps under 20. One epoch's per-offset step has split-half correlation 0.07 (full-set reliability
  ~0.14); the true per-pair corrections have SD ~0.04 nats against one epoch's noise of 0.074.
  Coarser schemes help little: one alpha_R-versus-rest offset per map reaches reliability 0.35.
* **What resolves is the aggregate**: the GLY|X alpha_L - alpha_R step averaged over its 38 maps,
  +0.060 +- 0.011, with only 6 maps individually above 2 sigma, the pattern of the AWH per-neighbour
  verdict (6.7).
* **By type**, GLY|X alpha_L is the best resolved per pair (reliability 0.52-0.71) and X|Y alpha_R
  is weak (0.18-0.19); X|Y beta against pPII has no detectable per-pair signal and resolves only per
  central residue type (0.40-0.52), the level of ff2.1's `sheet` energies.
* **About 70% of the noise is the protein set, not sampling**: 153 proteins run twice at identical
  offsets give per-map estimates correlated 0.69-0.71, so averaging epochs removes only the other
  30%.
* **The free-native mismatch is mostly its sample-size floor**: held out 0.0597 against random
  46-protein training subsets at 0.0585 +- 0.0021 (205 proteins 0.0320, 410 proteins 0.0249), so a
  held-out mismatch cannot show improvement or overfitting.

### 1.13 Which residues need a trained map: ff2.1's own mismatch and the literature (2026-09-28)

**Where ff2.1's free simulations miss the native basin populations** (epoch 0, 410 training
proteins, bootstrap over proteins; `checks/rama_by_type.py`): per central type at most 1.6 points
(0.09 kT), mostly alpha_R lost to pPII (1-1.6 points, z 3-4.6) and GLY alpha_L +2.1 (z 4.8); no beta
miss (largest z 2.3); no left neighbour has |z| > 3; right neighbour PRO alpha_R +3.2 (z 9.7), GLY
pPII +2.3 (z 5.9), VAL/ILE pPII -1.1 (z -4.5). The pre-PRO figure is mostly composition, not a
pre-proline defect (1.14).

**The left/right mixture gives pre-proline residues three times the alpha_R of NDRD's own
pre-proline map** (`checks/prepro_mixture.py`; 1,543 residues before PRO, ordinary left neighbour):
alpha_R right map 0.063, left 0.337, Upside's mixture 0.172, product rule 0.052, native 0.083, free
0.106. The mixing weights are nearly equal (0.83 typical, 0.92 for X|right|PRO), so an offset on the
pre-proline map reaches only the right map's share (1.14). Ting's product rule
(`upside_config --rama-library-combining-rule product`; S = 0.5-1.5 for proline) moves pre-PRO
alpha_R -12 points and others little (median largest basin change 3 points), but breaks GLY|GLY
symmetry: dividing by glycine's neighbour-averaged map (ln(aR/aL) -1.13) gives the middle glycine of
G-G-G +1.1.

**The reduced set** that the basin-offset ff_3.0 trained (`checks/reduced_set.py`, split-half;
adopted by the user 2026-09-28): GLY|X (aR, aL, beta), GLY|GLY (helix, beta; tied), X|right|PRO
(aR, beta), 158 offsets on 60 maps. GLY|X alpha_L per-pair reliability 0.72 (class mean +0.11 +- 0.02
nats); X|right|PRO alpha_R per pair 0.05, class mean +0.30 +- 0.04 nats, reproducible but the
composition of pre-proline sites (1.14).

**Literature** (survey 2026-09-28; [Ag] read by a sub-agent, [Abs] abstract only, otherwise full
text; the Peng 2022 PDFs are in `~/OneDrive - The University of Chicago/`). Pre-Pro is by far the
largest neighbour effect: N(i) and CB(i) clash with CD(i+1) in alpha (Ho & Brasseur, BMC Struct Biol
2005;5:14), and Ting's Table 6 shifts A -30.6, B +22.6, P +15.2 points (Ting et al., PLoS Comput Biol
2010;6:e1000763; its printed B and P phi ranges look swapped). Distinct classes are MolProbity's six,
general, Gly, trans-Pro, cis-Pro, pre-Pro, Ile/Val (Williams, Protein Sci 2018;27:293), plus partly
Ala. Other neighbour effects are small (non-Gly/Pro within ~12 Hellinger units, Ting Fig. 7;
aromatic/beta-branched neighbours shift mean pPII phi ~ -2 deg, Avbelj & Baldwin, PNAS
2004;101:10967), though context-dependent up to 4-fold in a restricted coil library (Jha et al.,
Biochemistry 2005;44:9691, where turn removal cuts the helical basin 37.0% -> 21.9%; Jha et al., PNAS
2005;102:13099: apoMb RDC correlation 0.41 -> 0.71). TCB is 62% turns (Ting Table 2), and right-hand
Gly and left-hand Pro effects reverse sign between turn and coil: placement, not intrinsic. The
glycine PDB map is asymmetric because PDB statistics record which residue wins a site
(Hollingsworth & Karplus, Biomol Concepts 2010;1:271). Force fields mostly use three classes,
generic, Gly, Pro (CHARMM36 [Ag], a99SB-disp, UNRES; "a global correction to the backbone is
sufficient for most residues", Best, de Sancho & Mittal, Biophys J 2012;102:1462 [Ag]); ff19SB,
RSFF1 ({E,Q,K,R,M,L}, {F,Y,W}, {V,I}) and Rosetta's own pre-Pro table go further [Ag]; coil libraries
break glycine's inversion symmetry (Choi & Pappu, JCTC 2019;15:1355 [Ag]). Beta propensity ranks
locally by sterics (Street & Mayo, PNAS 1999;96:9074, R = 0.92) but at edge strands is set largely by
tertiary context (Minor & Kim, Nature 1994 [Abs]; Avbelj & Baldwin, PNAS 2002;99:1309 [Ag]). Upside's
rama term comes from NDRD TCB, the sheet parameter added "to counteract an observed tendency for our
model to overstabilize helices" (Jumper et al., PLoS Comput Biol 2018;14:e1006578); FF2 made it per
amino acid (Peng et al., JCTC 2022;18:550, SI eqs. S1-S2). Not verified by the survey: Swindells,
MacArthur & Thornton 1995 numbers, the RSFF2 groupings, and how much turns inflate alpha_L for Gly,
Asn and Asp.

### 1.14 The X|right|PRO offsets were capped by the mixture and had no pre-proline signal to fit (2026-09-30)

On the finished `ff30_basin` run (six rounds; round 6 is the basin-offset ff_3.0 released and
withdrawn on 09-30, which `ff_3.0` means in 1.14-1.17): interior residues before PRO/CPR, centre not
GLY/PRO, 1,573 sites. Scripts `prepro_{leverage,control,residual,rules,left}.py`, logs
`*_20260930.log`, in `/project/trsosnic/yinhan/checks/`.

* **Six rounds moved the offsets steadily while the simulations stayed put.** The mean X|right|PRO
  alpha_R offset rose from 0 to +0.50 [-0.11, +1.63] (round 5) and +0.58 [-0.12, +1.89] (round 6);
  the mixed map's alpha_R fell only 0.174 -> 0.164, and free alpha_R stayed 0.098-0.108 (noise
  ~0.003) against native 0.077-0.080. Acting as an additive energy (the Newton step's assumption),
  the round-5 offsets would have reached the native 0.083.
* **The mixture caps what any offset on the right map can do**: an infinite alpha_R offset on every
  X|right|PRO map takes pre-PRO alpha_R only to 0.144-0.146; right-map-only gives 0.059, the product
  rule 0.052 (round-0 maps).
* **"Converged" meant "prior-limited".** The offsets stop only where the prior balances the unclosed
  gap, c* = N dp sigma^2 / T0 per map (T0 = 0.80): ~5 nats for ASP (144 sites, gap 0.027) and ALA
  (119, 0.036). The gate passed the rama group at step 114 (p 0.0011, 0.0006, 0.0158 at steps 76,
  95, 114; threshold 0.005) because the growing prior pull cancels more of the fixed data pull, not
  because the gap closed.

**The gap is mostly composition, not a pre-proline map error** (epoch 0; free / native alpha_R by
the residue's own native alpha_R):

| class | extended in native (aR < 0.05) | helical in native (aR > 0.5) |
|---|---|---|
| pre-PRO | n 1,401: 0.041 / 0.002 | n 116: 0.800 / 0.976 |
| other X | n 16,778: 0.058 / 0.001 | n 19,008: 0.915 / 0.990 |
| post-PRO X | n 622: 0.101 / 0.001 | n 818: 0.863 / 0.986 |
| GLY | n 2,200: 0.037 / 0.001 | n 516: 0.804 / 0.985 |
| PRO/CPR | n 1,020: 0.059 / 0.001 | n 701: 0.859 / 0.991 |

In every class extended residues visit alpha_R when free and helices fray. Pre-proline sites are 89%
extended (other X 47%), so their +0.023 gap is +0.036 extended and -0.013 helical; against ordinary
residues in the same native conformation they already visit alpha_R less (0.041 against 0.058). The
literature's pre-proline effect is a loop and unfolded-chain propensity, which the ConDiv target
cannot see (a native-restrained extended residue sits at alpha_R 0.002 whatever its class); the
offsets were fitting the class's composition.

**The glycine signal is real, and one map cannot serve it.** Extended glycines have no net alpha_L
gap (0.499 / 0.499, two opposite gaps, 1.15), helical ones 0.110 / 0.005; a per-map offset moves both
alike, and extended glycines outnumber helical four to one. The helical excess depends on where the
glycine sits, which a per-type (phi,psi) map cannot express. **The fixed point favours alpha_L
whatever the start**: native-restrained glycines over all sites sit at aR 0.19 / aL 0.40, so one map
per (glycine, neighbour) that reproduces the average glycine must favour alpha_L. The released map
did (`checks/gly_handedness.py`): over all 361 non-glycine flanks ln(aR/aL) median -0.92 at T = 1
(ff_2.1 -1.18), -1.23 at T = 0.8, alpha_R favoured in none; in glpG GLY136 (T-G-V) -1.59 at T = 1
(pre-release seed -0.73), GLY149 (R-G-E) -1.21. Its 09-30 glpG validation, where TM4's helical
glycines left for phi > 0, ran on seeds that paired ff_3.0's map, H-bond and pair tables with the
retired FF1-form coverage tables, so those flips belong to a mixed force field (3.11).

**Left-neighbour dependence of pre-proline alpha_R is not detectable** (`prepro_left.py`, 19 groups
of >= 30): native group sd 0.036 against 0.027 from sampling; right-only fits best (rms 0.036,
product 0.037, mixture 0.046); for pPII and beta the product rule follows the groups better (corr
0.59 / 0.67, right-only 0.41 / 0.59). Right-only changes no residue outside the pre-proline class and
keeps GLY|GLY exact; the product rule changes every residue and gives G-G-G +1.06 (1.13).

### 1.15 Which way the data pull glycine's handedness (probe from equal depth, 2026-09-30)

User question: from equal alpha_R and alpha_L depth, do the data pull glycine toward alpha_R or back
to alpha_L? **Measure** (`checks/glyprobe_analysis.py <run_output> <epoch>`): the data term of the
update on the 38 GLY|X maps over one epoch, prior excluded (it pulls toward NDRD):
(free_aL - native_aL) - (free_aR - native_aR) per residue read, positive = toward alpha_R, bootstrap
over proteins. On ff30_basin: epoch 1 +0.040 [+0.028, +0.052], epoch 5 (dL - dR +0.26) +0.023
[+0.012, +0.034], toward alpha_R in 31-32 of 38 maps. Free / native by native basin, epoch 5 and the
probe:

| glycines (non-GLY flanks) | n | alpha_R, ep. 5 | alpha_L, ep. 5 | pulls toward, ep. 5 | probe |
|---|---|---|---|---|---|
| helical in native (aR > 0.5) | 488 | 0.813 / 0.984 | 0.101 / 0.005 | alpha_R | aL 0.077 / 0.004 |
| alpha_L in native (aL > 0.5) | 1,003 | 0.030 / 0.002 | 0.889 / 0.990 | alpha_L | aL 0.847 / 0.989, aR 0.058 / 0.004 |
| the rest | 1,082 | 0.037 / 0.008 | 0.100 / 0.010 | alpha_R (aL), alpha_L (aR) | aL 0.069 / 0.008 |

One map serves two native populations that pull in opposite directions: helical glycines want less
alpha_L, natively left-handed loop glycines more, and the net follows their balance.

**Probe** (plan.md Phase 7, `training/ff30_glyprobe`, README there): one epoch from the ff30_basin
checkpoint (step 114) with each GLY|X map's alpha_R and alpha_L offsets moved by -d/2, +d/2 to equal
basin probability (dL - dR mean +0.31 -> +1.21; engine X-G-Y ln(aR/aL) +0.02 at T = 1, -0.05 at
T = 0.8, from -0.92 / -1.23), on the BP-fixed binary (|dE| <= 0.03 E_up, remote_jobs.md §0c). As
predicted, the pull turned toward alpha_L (stopped by the user at 15 of 19 steps, 329 proteins;
`checks/glyprobe_partial.py`, log `glyprobe_final_partial_20260930.log`): pooled -0.035 per residue
read [-0.048, -0.022], 35 of 38 maps. Interpolating with ff30_basin's +0.023 at +0.26 puts the fixed point near dL - dR +0.66,
about ln(aR/aL) -0.5 at T = 1, close to the all-residue native -0.58 (1.9; an estimate). **So a
context-free glycine map trained against native structures relearns the natives' placement and
cannot make glycine right-handed**: even at equal depth helical glycines keep 0.077 alpha_L against
0.004 native.

**The one context-aware term Upside has is FF2's H-bond energy** (`src/hbond.cpp`, `hbond_energy`):
each H-bond is scored by its donor's or acceptor's own (phi, psi), E_alpha for phi outside (0, 165)
deg and psi in (-120, 60), E_beta for that phi and other psi, E_other for phi in (0, 165), the three
shared by every type. ff_2.1 -1.961 / -1.946 / -1.769, a margin E_other - E_alpha of 0.192 per
H-bond; ff_3.0 -1.878 / -1.872 / -1.798, margin 0.080. **The margin shrinks without any glycine
training** (`checks/hb_trajectory.py`, checkpoints of runs from ff2.1):

| run | trains on the rama | margin along the run | E_alpha at the end |
|---|---|---|---|
| `ff21-fixedpoint` | nothing | 0.192 -> 0.127-0.142 by steps 13-25 | -1.907 (step 25) |
| `ff30_basin` | 158 basin offsets | 0.192 -> 0.151 at step 19 (offsets still 0), then 0.07-0.14 | -1.878 (step 114) |
| `ff30` (cancelled) | full glycine row | 0.192 -> 0.06-0.12 over steps 97-222 | -1.831 (step 222) |
| `ff30_glyhb` (round 3) | glycine H-bond offsets, side-chain lr / 10 | 0.192 -> +0.011 at step 18, negative from step 21, -0.036 at step 29 | -1.863 (step 29) |

At ff30_glyhb's step 29 E_other (-1.899) is the cheapest H-bond for every residue and glycine's own
margin is -0.315; why it drifts fastest is unproven. The populations did not follow
(`checks/hb_drift_20261003/nongly_basins.py`): helical non-GLY non-PRO free alpha_R 0.917 and 0.918
over epochs 0 and 1, alpha_L 0.0045 and 0.0042; a few hundredths of an E_up per H-bond is small
against kT at T 0.8-1.0. The margin drifts smoothly (Adam momentum), and E_alpha weakens in every
run, so most of the change is the trainer's own drift of the H-bond term from ff2.1, which the
offsets may add to but do not cause. Helical glycines' free alpha_L barely moved (0.107 at epoch 0,
0.101 at epoch 5) while dL - dR rose +0.26; whether the H-bond drift offset the maps' gain is not
separated. A glycine-specific set of these energies was therefore the smallest helix-aware glycine
term (round 3, 1.17); its limit is that cap and turn glycines are H-bonded and left-handed too
(1.16).

### 1.16 Inputs for an untrained glycine map and for a pre-proline rule (2026-09-30)

Measured for the proposal after 1.15 (glycine's map from physics, not the PDB). The pre-proline
items concern a problem that is not established (plan.md Phase 6, parked; lesson 10.8). Scripts are
in `/project/trsosnic/yinhan/checks/`, each with its log `<script>_20260930.log` or `_20261001.log`.

* **Native glycines are H-bonded in both basins, and the engine scores the two in different
  branches** (`gly_native_hbond.py`; 456 training natives, DSSP electrostatic criterion, H and O
  placed from N, CA, C; non-GLY flanks):

  | native basin | n | H-bonded (own NH or CO) | own NH donor | `hbond_energy` branch | commonest partners |
  |---|---|---|---|---|---|
  | alpha_R | 572 | 0.83 | 0.63 | helix 1.00 | NH->i-4 and CO<-i+4 |
  | alpha_L | 1,084 | 0.74 | 0.63 | turn 0.99 | NH->i-3, then NH->i-4 |
  | beta | 372 | 0.86 | 0.76 | sheet 0.98 | |

  Glycine-specific branch energies would put helical (E_alpha) and natively left-handed glycines
  (E_other) on separate parameters, where one map depth serves both (1.15); a helical glycine flipped
  to phi > 0 can still keep its NH->i-4 bond (the alpha_L C-cap pattern). TM4's helical glycines all
  have their own bond in glpG's seed: GLY136 (CO<-i+4), GLY143 (NH->i-4, CO<-i+4), GLY149 (NH->i-4).
* **What an own-H-bond term would miss is mostly fraying that every residue shows**
  (`gly_mismatch_by_hbond.py`; probe epoch 6 and ff_3.0 epoch 5; `own`: own NH or CO bonded,
  `spanned`: inside a short-range bond |d - a| <= 5, `none`; control: non-glycines in the same native
  basin and class). Native-basin loss, probe / ff_3.0:

  | natively helical | glycine share | glycine loss | non-glycine loss | glycine-specific excess, share of it |
  |---|---|---|---|---|
  | own | 0.82 / 0.84 | 0.099 / 0.133 | 0.057 / 0.056 | 62% / 66% |
  | spanned | 0.14 / 0.13 | 0.285 / 0.347 | 0.158 / 0.157 | 32% / 25% |
  | none | 0.04 / 0.03 | 0.365 / 0.531 | 0.270 / 0.260 | 6% / 9% |

  Natively left-handed glycines lose less than the non-glycines at alpha_L (mostly Asn, Asp): 0.12
  against 0.17 (own), 0.29 against 0.41 (none), so their pull toward alpha_L in training is the
  generic loss.
* **Per-type misses in the same context follow intrinsic propensity** (`type_mismatch_by_context.py`;
  ff_3.0 epoch 5, 40,462 residues): helical own-bond loss GLY 0.133 (z +6.0), SER 0.082, ASN 0.080
  high, GLU 0.041, ALA 0.043 low (mean 0.058); extended VAL 0.072 (z -11.1), ILE 0.077 low (mean
  0.112). The orders match the helix and beta propensity scales (glycine is the weakest helix former
  after proline), so a type-and-context correction fitted to the native-restrained target (~0.99 in
  every basin) would flatten them: placement again, in a milder form.
* **Residue counts per basin show selection but do not convert into map energies**
  (`aa_basin_counts.py`; 49,401 interior residues): glycine fills 54% of alpha_L sites, 82% of
  phi > 0 extended sites, 2.7% of alpha_R. Were residues chosen by (phi, psi) alone, every reference
  X's library would give glycine the same `h = E(aL) - E(aR)`; instead h_GLY runs from -2.61 (Ala) to
  -1.16 (Asn), median -1.74, sd 0.50 over 14 references (library -1.20, AWH -0.15 to -0.3), largest
  from helix-placed anchors (Ala, Leu, Met), smallest from turn-placed (Asn, Asp, Ser), and -1.95
  exposed to -0.17 buried (2 references with >= 10 alpha_L counts). Selection depends on more than
  (phi, psi) (rule 1.8), so a count-based correction needs an anchor from another map.
* **Proline shows no problem a map change would fix** (`pro_mismatch_by_context.py`): central PRO
  holds its extended basin best of all types (own bond 0.054 against mean 0.114, z -8.7); pre-PRO
  against the same type elsewhere differs by +0.024 [+0.012, +0.035] extended (n 1,219) and +0.059
  [+0.025, +0.106] helical (n 123). Loop and unfolded-state propensity, where the mixture's extra
  alpha_R would act, is invisible to this comparison and stays untested.
* **Right-only for pre-proline residues can be written into the library**: every X|right|PRO
  `dimer_weight` x 1e6 in coil and sheet gives the right map to 1e-5 E_up, every residue not before
  PRO bitwise unchanged; all weight readers go through `read_rama_maps_and_weights`. On 1,959 native
  pre-proline residues (`prepro_rightonly_natives.py`) it takes engine alpha_R 0.165 -> 0.055 and
  costs the native point -0.17 (median) for 1,742 extended, +0.84 [10%: +0.40, 90%: +1.37] for 132
  helical.
* **Glycine already has a side-chain bead**: a fixed `GLY_0` in `sidechain.h5` with trained pair and
  coverage rows, 0.61 A from the frame origin along L-CB (ALA's bead 1.73 A out). `backbone_pairs`
  gives glycine no CB. `ProteinHBond` holds per-pair H-bond values and per-edge sensitivities
  (`igraph`), so a term on the residues a bond spans could be built on its edge loop.
* **The AWH glycine-before-proline surface is unlike the others**: Ac-Gly-Pro-NHMe (`RP`, 400 ns,
  ff99SB-ILDN) alpha_R 0.006, alpha_L 0.004 (other right contexts 0.10-0.20); NDRD's GLY|right|PRO
  has alpha_R 0.046, and the AWH build of `rama31.dat` gave it the pooled surface (0.105), erasing
  glycine's pre-proline clash.
* **`rama_map_pot_ref` reshapes a glycine map.** ConDiv adds it to every residue; on the AWH GLY|ALA
  it moves alpha_R 0.105 -> 0.145 and alpha_L 0.121 -> 0.163 (extended weight to the helical basins),
  ln(aR/aL) almost unchanged (-0.139 -> -0.116). A measured surface used as glycine's whole local
  term must be stored with the reference subtracted, or glycines left out of that node.

**Literature** (sub-agent surveys 2026-09-30 and 10-01; [Abs] abstract only, the rest full text read
by the agent; not re-checked unless stated).

* **PDB statistics cannot give glycine's intrinsic map**: they are Boltzmann-like only between
  residues at a fixed (phi,psi), so glycine's alpha_L enrichment says it beats the others there, not
  that it prefers alpha_L to alpha_R (Shortle, Protein Sci 2003;12:1298 [Abs]; Hollingsworth &
  Karplus 2010: "a dipeptide with Gly in it must have equivalent energetics in the delta' and delta
  regions").
* **No model fixed handedness by training a context-free map.** Models that avoid the problem take an
  inversion-symmetric glycine term from physics and get context from other terms: UNRES (MP2 PMF of
  Ac-Gly-NHMe, Sieradzan JCTC 2012;8:4746; Gly-Gly near-symmetric, Lipska JPCL 2023), CHARMM36 (QM
  glycine CMAP, Best 2012), ff19SB (aqueous QM dipeptide, as PDB enrichment "would be reflected
  erroneously", Tian 2020), CGSchNet priors from ff99SB-ILDN MD (with which alone every protein
  unfolds; Charron 2025), Martini3-IDP, ff24EXP-GA (iterative Boltzmann inversion to GGG; Suresh JCTC
  2025), Rosetta `-symmetric_gly_tables` (its default keeps the PDB asymmetry, ln(aR/aL) about -1.6)
  and Choi & Pappu 2019. AWSEM drops glycine's Ramachandran term and sets the i->i+4 helical H-bond
  strength per residue from experimental helix propensity, and HPS-SS fits a per-residue dihedral
  term to host-guest propensities (Rizuan, JCIM 2022), so a per-type context term can have a
  placement-free target. Force fields disagree on glycine pPII (ff14SB 0.36, CHARMM36m 0.48). The
  Hamelryck reference ratio returns the contrastive-divergence fixed point unless its non-local
  feature carries context (the agent's inference). For pre-proline, Rosetta's `rama_prepro` replaces
  the table whenever residue i+1 is Pro, with no left-neighbour information; CHARMM36/36m get the
  effect from explicit Pro CD sterics, which Upside lacks.
* **No experiment measures glycine handedness in a chiral context** (searched: stereospecific
  3J(HN,Ha2/Ha3), RDCs resolving the sign of glycine phi). GGG is achiral, so handedness rests on MD
  alone, where our two force fields agree to 0.045 nats (9r).
* **GGG experiment against our AWH surfaces (2026-10-01)**: Andrews et al., Biomolecules 2020, SI
  Tables S1-S2 (five J-couplings for the central glycine of cationic GGG in water; basin populations
  from a Gaussian model fitted to them and amide I' spectra, a rough comparison per the authors).
  Mirror boxes: pPII -90 < phi < -42, 100 < psi < 180; beta-t -130 < phi < -90, 130 < psi < 180;
  a-beta -180 < phi < -130, 130 < psi < 180; alpha -90 < phi < -32, -60 < psi < -14:

  | GGG central glycine | pPII | beta-t | a-beta | alpha |
  |---|---|---|---|---|
  | experiment (Gaussian model) | 0.46 | 0.13 | 0.01 | 0.06 |
  | ff14SB, cationic GGG (Andrews) | 0.40 | 0.06 | 0.09 | 0.05 |
  | our AWH Ac-Gly-Gly-NHMe, ff99SB-ILDN rep1 / rep2, ff14SB | 0.32-0.33 | 0.05 | 0.06 | 0.09 |

  The matched force field is within 0.06 of experiment in pPII and 0.01 in alpha; our capped
  dipeptide differs by about 0.08, from the termini, not the force field.
* **Combining neighbours.** Ting's product rule `f(C,R) f(C,L) / [S f(C)]` keeps empty any region the
  right map empties. On 17,600 held-out coil residues it scored 1.25, centre plus right only 1.21, raw
  triplets 1.19, with no detectable left-right interaction in 3J couplings (Shen, Roche, Grishaev &
  Bax, Protein Sci 2018;27:146). Pre-Pro mechanism: N, O(i-1) and H(i) clash with CD(i+1) (Ho &
  Brasseur 2005).
* **AlphaFold neither meets nor solves this problem**: AF1 predicts torsions per residue from
  sequence and MSA (Senior, Nature 2020;577:706); AF2 has no Ramachandran prior and leaves chirality
  to FAPE (Jumper, Nature 2021;596:583). A conditional predictor may learn placement because placement
  is its target; Upside needs a transferable local energy.
* **Helix propensity:** Pace & Scholtz, Biophys J 1998 [Abs]: Gly 1.00, Ser 0.50, Asn 0.65, Ala 0
  kcal/mol over 11 hosts at solvent-exposed mid-helix sites; it agrees in order with the per-type
  fraying above but cannot calibrate a single host-guest simulation.
* **Not found:** an MD or QM free energy for Ac-X-Pro-NHMe; a measurement of left-neighbour effects on
  pre-Pro alpha_R; per-residue pre-Pro (phi,psi) in Pro-kinked helices.

### 1.17 Round 3, the physics glycine map: the selection panel, the first epoch, and glycine's own H-bond offsets (2026-10-01/02)

**The selection panel** (plan.md Phase 8; `/project/trsosnic/yinhan/ff3_selection`): the 44 all-L
CATH domains of Charron et al. (2,823 residues, 198 glycines; ff99SB-ILDN, 300 K, 20,000 frames
each) against Upside native-start runs at T 0.8 (4 x 8,000 tu per domain after equilibration). The
error is Upside minus all-atom population of each residue's all-atom majority basin in folded frames
(Q above the all-atom 5th percentile); SE by bootstrap over domains; checkpoints compared paired.
Released ff2.1 (43 of 44 domains, ff2.1 unfolds 4hwiB01; 0.615 of Upside frames are folded):

| class | residues | ff2.1 - all-atom |
|---|---|---|
| helical, non-glycine | 1,358 | -0.018 +- 0.004 |
| beta, non-glycine, extended region (beta + pPII) | 569 | -0.017 +- 0.004 |
| beta, non-glycine, beta basin alone | 569 | -0.123 +- 0.011 |
| helical glycine | 34 | -0.081 +- 0.049 |
| left-handed glycine | 74 | +0.029 +- 0.019 |

**The beta-basin miss is a phi shift, not strand loss.** Of the 0.123, 0.106 moves to pPII and 0.015
to alpha_R; strands stay extended 0.971 of the time against 0.989; mean strand phi (Upside -114.4,
all-atom -125.5) straddles the crystal's -121.0; beta-rich domains are not less stable (folded 0.63
against 0.60). So the panel scores beta by the extended region, and ff2.1's one outlier is helical
glycines (~8 points of alpha_R, n 34, so noisy).

**The physics (AWH) map alone** (`ff21_awh`) keeps ff2.1's helix, strand and helical-glycine errors
but costs natively left-handed glycines ~9 points of alpha_L (-0.065 against +0.029): without the
NDRD map's bias ff2.1's non-local terms do not hold them. What is missing is a glycine context term,
not the map.

**One epoch of training (ff30_gly, panel `e00`) moves helices and helical glycines away from
all-atom.** The control is the ff2.1 workflow's own epoch 0 with the NDRD map (`ff21-fixedpoint`,
panel `fp_e00`). Margin = E_other - E_alpha:

| | domains | folded | helix | beta | helical Gly | left-handed Gly | margin |
|---|---|---|---|---|---|---|---|
| ff2.1 released | 42 | 0.615 | -0.018 | -0.018 | -0.081 | +0.034 | |
| ff2.1 + AWH map (start) | 42 | 0.585 | -0.021 | -0.019 | -0.081 | -0.066 | |
| e00 (one epoch) | 42 | 0.487 | -0.030 | -0.018 | -0.136 | -0.056 | |
| ff2.1 released | 41 | 0.615 | -0.017 | -0.018 | -0.081 | +0.035 | 0.192 |
| ff2.1 workflow, NDRD map, epoch 0 | 41 | 0.531 | -0.025 | -0.018 | -0.077 | +0.023 | 0.129 |
| ff30_gly (AWH map), epoch 0 | 41 | 0.487 | -0.028 | -0.018 | -0.136 | -0.058 | 0.057 |

Single-file resets (e00 with one file put back to ff2.1's, paired against e00):

| e00 with ff2.1's ... | folded | helix | helical Gly | resolved against e00 |
|---|---|---|---|---|
| hbond.h5 (e00_hb21) | 0.499 | -0.029 | -0.087 | helical Gly |
| bb_env.dat | 0.508 | -0.029 | -0.121 | none |
| environment.h5 | 0.490 | -0.027 | -0.115 | none |
| sheet | 0.503 | -0.028 | -0.095 | helical Gly |
| sidechain.h5 (e00_rot21) | 0.583 | -0.024 | -0.078 | helix, helical Gly |

* **(A) The workflow itself weakens helices and folds in its first epoch, with either map, through
  the side-chain pair update.** `e00_rot21` restores helix, helical glycines and folding and differs
  from the start in no class. The update is a noise-driven random walk (ff30_gly steps 0-18, 31,420
  rot coefficients): `||mean g|| / mean||g||` 0.107 against 0.229 for noise, 6% of coefficients
  beyond 2 SE (chance ~5%), cosine with the mean gradient +0.12 (hb +1.00, burial +0.83). Adam steps
  every coefficient by ~lr whatever its sign consistency, so a table at its fixed point diffuses: rms
  0.44 in 19 steps. Hence ff3.0's 10x smaller side-chain learning rate.
* **(B) The AWH map adds the glycine-specific damage.** `e00_hb21` returns helical glycines to the
  start's value, so their extra loss is the H-bond drift (the `sheet` and `sidechain.h5` resets also
  restore them, so it needs those changes together); the left-handed loss is the map's own.
* So ff2.1 is not at this trainer's fixed point (FF2's training stopped at 76 iterations with no
  convergence test, 1.23), and a release judged "no worse than ff2.1" cannot pass while (A) stands.

**Gradient split** (ConDiv's worker with NSE and DSE apart; 25 step-0 proteins, ff2.1 + AWH map;
`/project/trsosnic/yinhan/checks/gradsplit_20261001`; sums over proteins; the step is NSE + 0.3
DSE):

| | NSE | 0.3 DSE | contrast | effect |
|---|---|---|---|---|
| E_alpha | +21.45 | -40.54 | -19.09 | weaker helix H-bonds |
| E_beta | +9.54 | -28.70 | -19.16 | weaker sheet H-bonds |
| E_other | +4.56 | -2.31 | +2.24 | stronger left-handed H-bonds |
| E_bias | +9.02 | -10.44 | -1.42 | |
| bb_env scale | -47.68 | +104.65 | | DSE larger in every component |

* DSE pulls E_alpha, E_beta and bb_env weaker (the unfolded ensemble near Tm keeps residual helical
  H-bonds the SARW reference lacks), the same under the NDRD map (E_alpha -40.2 against -40.6), but
  resetting those groups leaves the panel unchanged, so ff30_glyhb does not act on it.
  ff2.1 balances at lambda 0.10-0.15, not 0.3, in all three; one common factor suggests the unfolded
  ensemble differs from what ff2.1 was trained against.
* Peng's DSE threshold as coded (SI Fig. S3, the port's and Kleinmann's code: 2/3 Rg_low + 1/3
  Rg_high of the coldest and hottest replicas) and the text's (2/3 Rg_native + 1/3 Rg_SARW) differ by
  ~2% in every DSE group (`gradsplit_20261001/threshold_report.py`). "Threshold as coded" stands.
* **The E_other pull is glycines'**: natively left-handed glycines give +5.25 of the NSE's +4.56,
  non-glycine residues net negative; the AWH map flips the NSE pull from -2.55 (NDRD) to +4.60
  (paired +0.30 per protein [+0.16, +0.44]). Rotamer pairs are NSE-dominated (|NSE| 58 against
  |0.3 DSE| 6.5); sheet has no DSE part by construction.

**Why glycine got its own H-bond offsets.** With the physics map, the trainer's only
(phi,psi)-resolved knob for natively left-handed glycines is the shared branch energies, so it makes
left-handed H-bonds cheaper for every residue and the helical glycines pay (lambda's helix 3, 11d).
The three branch energies are FF2's ConDiv freedom (1.23; the FF1 trainer fitted one scalar, -2.112,
9e), shared by all 20 types. Round 3 gave glycine its own dE_alpha, dE_beta, dE_other, trained from
zero (plan.md Phase 8): an optional 15-entry `parameters` with a per-residue `residue_class`; a
12-entry config stays bitwise unchanged (10.12). The released ff3.0 has none (its `hbond.h5` keeps 12
entries; `training/ConDiv.py`).

**They did not deliver.** On 37 domains (helix 1135, beta 522, helical glycine 26, left-handed
glycine 63 residues), folded / helix / gly_helix / gly_left: ff21_released 0.615 / -0.016 / -0.114 /
+0.032; ff21_awh, the run's start, 0.585 / -0.017 / -0.087 / -0.077; e02 (ff30_gly) 0.501 / -0.025 /
-0.112 / -0.042; ff30_glyhb epochs 0 and 1 (`h00`, `h01`) 0.460 / -0.029 / -0.130 / -0.075 and
0.469 / -0.030 / -0.122 / -0.054 (helix SE 0.005-0.007, glycine SE 0.02-0.07; beta -0.017 to -0.018
throughout). Both epochs are resolved worse than the start in helix, and neither glycine class moves
beyond the noise (26 helical glycines cannot resolve the change sought). The release rule never
picked an ff30_glyhb checkpoint ("RELEASE HELD for the user: worse than ff2.1 in helix"). The
training-ensemble readout (natively helical non-glycine free alpha_R 0.917 in both epochs) disagreed
because it averages over the epoch's changing parameters, while the panel tests the end checkpoint.
The shared H-bond margin fell +0.192 -> +0.011 -> -0.034, glycine's to -0.203 and -0.350: H-bond
drift is the leading, untested candidate; the side-chain cause is unlikely (rot change rms 0.058
after 19 steps against 0.44 in ff30_gly's 19, yet the same -0.030 helix loss).

**Engine and data facts gathered for a bottom-up fallback** (proposed 10-01, not approved):

* ConDiv's derivative step takes any N/CA/C positions (`compute_divergence`,
  `training/ConDiv.py:399`), so mapped all-atom frames can replace a simulated ensemble without engine
  changes.
* **The rotamer solve has no temperature.** Probabilities are `exp(-E)` (`src/rotamer.cpp:252`,
  `:896`), so side-chain free energies are at T = 1, as are Rama maps (-ln P). Compare any all-atom
  target at an explicitly stated Upside temperature.
* A force-matching gradient is available by one directional central difference of `get_param_deriv`
  (1UBQ, ff_2.1, 1 A RMS displacement): at steps of 3e-4 and 1e-4 A values agree to 0.9-1.4e-3
  (`hbond_energy`), 1.4-4.7e-3 (`hbond_coverage`), 3.2-6.6e-3 (`sigmoid_coupling_environment`),
  0.6-1.1e-2 (`rotamer`), not for `hbond_coverage_hydrophobe` (0.22-0.28); the energy itself needs
  <= 1e-3 A to agree with `deriv()` to 0.1%.
* **Charron et al.'s deposit** (Nat Chem 17, 1284 (2025), SI 1.1; Zenodo doi:10.5281/zenodo.15465782;
  amber99sb-ildn/TIP3P, 300 K, as our AWH dipeptides): 50 CATH domains, 4 x 0.5 us from native,
  coordinates and forces every 20 ps, none in the 456-protein `training/pdb_list`; 1,100
  adaptive-sampled octapeptides, 1 us each. Six domains carry D residues and are excluded (the
  README's six, confirmed by CB chirality on the `_eq.pdb` files: 2ga1A02 ASP102; 3e6zX01 SER17,
  GLU28, LYS85; 3luyA02 LEU117; 3tj8A02 LYS120; 4jriB00 LEU45, SER65; 4npsA02 CYS287); the other 44
  are the panel. `training_a_cg_model.zip` (65 GB) holds 5-bead frames plus "decoy molecules" (every
  50th frame, 0.5 A noise, zero force) that are not data; its CATH member
  `DECOY_nicks_transferable_cath_delta_dataset.h5` (bytes 17270634265 to 52854221219) is fetchable
  alone by HTTP range. Raw all-atom trajectories are not deposited; `training_data_generation.zip` ->
  `cath_generators.tar.gz` has `pdbs/<id>.pdb`, `<id>_eq.pdb` and AMBER topologies (named `1ldjA06`,
  `2hbpA00` where the SI reads 1ldjA02, 2hbpA02). mdCATH and the D. E. Shaw fast folders use
  CHARMM22*.
* Charron et al.'s glycine (`model_and_prior.pt`) has fixed 1D phi and psi priors, no 2D map and no
  chirality improper; the phi prior is left-handed (minimum +76.5 deg, P(phi > 0) 0.603, between our
  ff99SB-ILDN dipeptides' 0.4998 and NDRD's 0.6548), but their network, fitted to all-atom mean
  forces, corrects the prior, while in Upside the map is the energy term. Their release (epoch 73) was
  chosen by simulating epochs on the fast folders (training loss "alone is not a sufficient metric").
* Every bottom-up protein model the 10-01 survey checked needed a top-down correction (Liwo PNAS 2002;
  Hills, Lu & Voth PLoS Comput Biol 2010; Majewski Nat Commun 2023; Charron 2025); a99SB-disp's
  glycine is PDB-refit and must not be the reference; Upside's side-chain energies are PDB maximum
  likelihood, deliberately not matched to atomistic energies (Jumper et al., PLoS Comput Biol 14,
  e1006342, 2018).
* **A WebFetch summary invented a methods paragraph** for Charron (ff14SB, OpenMM, 1 us, MSMBuilder,
  in quotation marks): quote a paper only from the downloaded file.

### 1.18 Glycine map facts that still hold (measured 2026-09-16 to 09-25)

* **Glycine is the only residue whose library handedness cannot be local physics.** ff2.1 coil
  library, central residue, averaged over the 20 left neighbours (others in between):

  | res | alpha_R | alpha_L | P(phi>0) | dG(aR->aL) kT |
  |---|---|---|---|---|
  | GLY | 8.98% | 30.95% | 0.650 | -1.238 |
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

  Glycine's alpha_L is 7x the 4.28% mean of the other nineteen; every other residue has positive dG,
  ordered by C-beta branching (local sterics); glycine has no C-beta, so its handedness has to come
  from context. Asn (13.4% alpha_L) is the only other residue worth checking.
* **The coil/sheet mixture cannot change glycine's handedness.** Sweeping `sheet_mixing_energy`
  (`read_weighted_maps`) from +4 to -10 leaves ALA-GLY-ALA's dG(aR->aL) at -0.971 to three decimals,
  the sheet map being empty in both helical basins (GLY|ALA alpha_R 1.6e-10, alpha_L 8.7e-16).
  Tuning the mixing to P(phi > 0) = 0.5 is a trap that inflates beta.
* **The library's glycine alpha_L describes folded proteins accurately**: 42 of 67 interior glycines
  of the 16 benchmark natives sit at phi > 0 (63%; NDRD 66%), placement at loop sites (1.9), which
  does not belong in a local energy.
* **No NDRD release is near zero.** Our coil group is exactly `NDRD_TCB` (correlation 1.00000 against
  `GLY|ALL`); central-GLY dG(aR->aL) over `GLY|X` is -1.876 for Conly (coil only, 13,945 residues),
  -0.839 Tonly (turns only, 27,532), -0.965 TCB (44,112), -0.410 TCBIG (adds pi and 3-10 helix,
  62,345). The purest coil subset is the most biased, so "turn contamination" does not explain the
  bias.
* **The map reaches every glycine regardless of context.** `rama_map_pot` reads only `rama_coord`
  (datasets `rama_pot`, `residue_id`, `rama_map_id`, `rama_map_id_all`); repeated triplets get
  bitwise-identical maps (checked on glpG).
* **glpG's glycines**: 23 in 210, no GGG, GG pairs 96-97 (inside TM2) and 132-133 (TM4's N-cap);
  TM4's helical 136 (T-G-V), 143 (M-G-Y), 149 (R-G-E) are XGX; TM1 has none. ff2.1 biases every glpG
  glycine toward alpha_L by 0.47-1.23 E_up.
* **`rama_map_pot_ref` is one residue-independent map added to every residue**, with negligible
  handedness (dG(aR->aL) -0.0139 E_up; mirror differences up to 0.52 at single grid points are not a
  basin free energy), but it reshapes a glycine map's basin weights (1.16).
* **Which library number is which**: -1.238 is `X|GLY` over all 20 left neighbours, -1.318 over the
  8 AWH-measured ones, -1.13 to -1.18 the `GLY|ALL` marginal; glycine's `dimer_weight` is 0.908 and
  0.948 for the two directions, not 1.0, so the left/right mixture matters. Maps store -ln P with
  sum(exp(-E)) = 1 over 72x72, so an AWH PMF enters as PMF / kT(300 K), not PMF / 2.914952774272
  (GLY_sym.md §5). The two `GLY|ALL` maps feed only the unused `product` rule; exclude them from
  per-map statistics.

### 1.19 Glycine's in-chain handedness in all-atom peptides, and the BioEmu-fitted library (2026-10-04)

The BioEmu-map line built here (the current central-glycine row of `rama31.dat`, the ff30_bio runs)
was dropped by the user on 10-11 for the basin-depth path (plan.md); ff3.0 is NDRD plus one trained
depth pair (1.26-1.27). The measurements still hold.

**Octapeptides** (1,100 L-only, ff99sb-ildn, 300 K; `scratchpad/ff3_local_test/scripts/opep_*.py`,
`scripts/bioemu_phipsi.py`):

* **Charron's adaptive frames are seed-biased; BioEmu's plain MD is the target.** Charron's peptides
  (`training_a_cg_model.zip`, md5 verified) are ~100 consecutive ~10 ns segments, and L residues'
  alpha_L decays inside them (alanine 0.283 -> 0.163). BioEmu's plain MD of the same peptides (Lewis
  et al., Science 2025; Zenodo 15641199, `ONE_octapeptides`; run001-005 x 1 us, a frame per 10 ns,
  first 100 ns per file dropped) gives glycine -0.54 where the adaptive frames give -0.70, although
  glycine did not drift inside segments: a flat within-segment profile does not prove an unbiased
  ensemble.
* **Glycine inside a chain favours alpha_L.** ln(aR/aL) (`training/rama_basin.py` basins, SE over
  glycines): interior X-G-Y **-0.54 +- 0.04** (365), positions 2-5 -0.53 +- 0.05 (249); next to the
  N-terminal residue -0.09, the C-terminal -0.88. At residues 2-5 by context: G-G-Y -0.19 (27), X-G-G
  -0.87 (25), X-G-P +1.06 (10; aR 0.017, aL 0.006), G-G-G -0.38 (5), so GLY|left|GLY, GLY|right|GLY
  and GLY|right|PRO are fitted on their own context.
* **Fitting rules.** A fitted glycine entry is fitted once in the coil group with NDRD's coil top and
  copied to sheet; a sheet top from NDRD's strand map left extra beta no coil update could remove
  (0.135 against 0.078). A map reproduces its source at T_up = 1 (up.md 2.8), i.e. -ln P(300 K), so a
  fit runs Upside at T_up = 1 with update factor 1; at T_up = 0.8557 it scales the map by 0.8557, the
  1.169x error up.md warns about.
* **The reference force field puts L residues at phi > 0 more than the PDB does** (plain-MD
  P(phi > 0), positions 2-5: Ala 0.095, Ser 0.187, Lys 0.113, Leu 0.099, Asn 0.086; NDRD alanine
  alpha_L 0.045), and ff14SB and others overpopulate glycine's helical basins against GGG couplings
  (Andrews 2020), so ff99sb-ildn may overstate glycine's in-chain alpha_L preference; the target is
  only as good as this force field, which is also the panel's.
* **Against the library**: NDRD central glycine -1.15 pooled, -0.62 to -1.75 per neighbour map; of
  that ~-1.15, ~-0.6 is selection and ~-0.55 in-chain physics in this force field. Part of glycine's
  alpha_L excess is therefore real local physics (turns and caps with L neighbours), and the
  selection-free target is the in-chain distribution, not the isolated dipeptide's.

**The non-glycine maps stay NDRD** (user, 2026-10-04/05: all maps are PDB statistics;
`/project/trsosnic/yinhan/checks/gly_bioemu_map/allres/`):

* Measured in Upside as glycine's map was fitted (`allres_run.py`, `allres_compare.py`; ff2.1 terms,
  T_up = 1, residues 2-5, ln(P_upside / P_BioEmu) per basin), every L type is less helical (median
  -0.31), more extended (+0.14 to +0.62) and mostly lower in alpha_L (Ile -2.11 to Ala -0.97; His
  +1.09, Asn +0.43) than BioEmu: gaps that match or exceed glycine's under NDRD (alpha_R -0.35,
  alpha_L +0.46, extended -0.10), so size does not single glycine out.
* The reference is not better than NDRD for L residues (`allres/lit/`): ff99sb-ildn changed only
  side-chain torsions (Lindorff-Larsen 2010); ff99SB oversamples beta against pPII (Wickstrom 2009)
  and "uniformly over-emphasize[s] alpha R" (Beauchamp 2012); alanine is ~80% pPII in GxG (Hagarman
  2010), ~90% in trialanine (Graf 2007); on 256 dipeptides' 3J(HN,Ha) ff99sb-ildn scores r^2 0.56
  against 0.82 for coil-fitted RSFF2 (Li & Elcock 2015). Alanine alpha_L: ff99SB Ala5 ~4% (Best &
  Hummer 2009), coil libraries ~5%, Upside 0.035, BioEmu 0.093.
* In folded proteins (`allres/panel_all_types.py`) natively alpha_R non-glycine residues lose
  <= 0.03 of alpha_R, while natively extended ones leak to alpha_R (Phe -0.127 +- 0.040); a BioEmu fit
  deepens L alpha_R by 0.1-0.7 and alpha_L by up to 2 (Ser 0.024 -> 0.184), which would worsen the
  errors that exist.
* Glycine differs because selection reverses the sign of its local preference, the reference is
  credible for it (physics force fields build it near-symmetric), and its defect shows in folded
  proteins as helical-glycine loss (lambda's helix 3). The L residues' gap is mostly a shared
  helix-against-extended offset whose sign experiment does not support, and replacing it would change
  what ff2.1's trained terms were balanced against.

**Literature beyond 1.16** (sub-agent surveys, checked on publisher, PubMed or PMC pages): a PISCES
count (<= 1.2 A, 884 chains, by a sub-agent) finds alpha_L glycines 4,944 vs alpha_R 3,161, 25% at
helix C' and 30% at type II turn i+2, and ln(aL/aR) = -0.34 without those motifs; both are one-sided
because their mirrors need an L residue at phi > 0, the chirality transfer the octapeptides show.
Jumper thesis 4.3: the Rama term is a "naive Boltzmann inversion" of NDRD TCB (62.4% turns), and the
reference, a free alanine chain's density, is applied to glycine too; native-trained torsions give
glycine low energy at phi > 0 (Greener & Jones, PLoS ONE 2021). Experiments to reproduce: Gly -> D-Ala
at phi > 0 sites +0.6 to +1.9 kcal/mol (Anil JACS 2004); Ala over Gly in helix interiors 0.4-2
kcal/mol (Scott PNAS 2007); lambda G46A/G48A Tm +6.5 C (Liu, Gao & Gruebele JMB 2010; Y22W). 1LMB:
G46, G48 alpha_R, G41, G53 (helix 3 C') alpha_L; 2XOV: TM4's G202/G209/G215 (our 136/143/149)
alpha_R.

**The alpha_R/alpha_L barrier: the engine reads the map faithfully, but NDRD's top is a floor, not a
peak** (`scripts/rama_barrier.py`, `rama_spline_check.py`, `gly_flip_route.py`).

* `rama_map_pot`'s periodic spline has no cap or default: on-grid engine values equal map + reference
  to 4e-4; half-cell dips are at most 0.056 (NDRD) and 0.019 E_up (AWH).
* A flip crosses the peak on its own route, which in a helix is not the bare map's minimax path
  (user). The one crossing caught in a 10-tu frame (glpG GLY143, pre-fix inputs, 3.11) sat at
  (phi, psi) = (-21.5, -86.5), +8.65 on its map; flipped glycines dwell at phi 70-100, psi -50 to
  -110 before alpha_L.
* **NDRD has no data where glycine does not go and sits on a floor of ~9-11 E_up there.**
  GLY|left|ALA, relative to alpha_R (engine frame):

  | map | phi = 0 column min | E(0, -40), helical psi | E(0, 0) | via phi = 180 | alpha_L - alpha_R |
  |---|---|---|---|---|---|
  | NDRD (ff2.1) | 8.1 (psi +90) | 10.1 | 9.6 | 5.7 | -1.19 |
  | AWH (rama31) | 9.4 (psi -90) | 14.4 | 17.2 | 3.1 | -0.16 |
  | octapeptide surface (ff21_oct) | 6.6 | 6.9 | 7.7 | 2.8 | -0.86 |

  Crossing phi = 0 at helical psi costs 4.3 E_up less on NDRD than on the physics surface (~200x in
  rate at T 0.8) and need not unwind the helix.
* **Upside never designed the transition region, for any residue** (Jumper thesis 4.3, 4.3.1): the
  reference correction targets populations; unvisited regions are NDRD's smooth density tail (no cap;
  ~500 distinct values in each map's top 10%) at similar heights (phi = 0 column minimum above the map
  minimum from 5.8, Asp, to 9.1, Pro; Gly 7.8, Ala 7.1). Contrastive divergence cannot set barriers
  (its gradient on a cell is native-minus-free population, ~0 where nothing goes). Rigid scans with
  sterics hit chain-specific clashes and are not free energies; minimax barriers barely move (glycine
  via phi = 0: 8.96 with sterics, 9.20 without). Decision (user): leave the top as the library has
  it; a fitted glycine map keeps NDRD's glycine top with its barrier height above alpha_R unchanged.
* In all-atom (Charron's segments, 0.1 ns frames, residues 2-5) alpha_R <-> alpha_L crossings run
  ~134 per glycine per us; of those resolved 84% go through the extended region (21,954), 16% across
  phi = 0 (4,240); 42,813 fall inside one frame.

**The library and the round-4 decision (2026-10-05).** The library (`parameters/common/rama31.dat`,
described in up.md 2.8; data in `/project/trsosnic/yinhan/checks/gly_bioemu_map/`) reproduces BioEmu
in Upside (residues 2-5): X-G-Y aR / aL / beta 0.166 / 0.286 / 0.082 vs 0.168 / 0.286 / 0.078,
ln(aR/aL) -0.55 vs -0.535; G-G-Y -0.37 vs -0.19, X-G-G -0.83 vs -0.87 (within BioEmu's SE of
0.09-0.10); its engine X-G-Y map has ln(aR/aL) -0.58 (ff2.1 -1.10). The push probe from it (1.15's
data term, ff30_glyhb's worker, no update, the 72 proteins of its minibatches 0-2) gave d =
**-0.014 [-0.041, +0.012]** (441 X-G-Y glycines; natively helical +0.169, natively left-handed
-0.172): the data come to rest near BioEmu's L-R difference, where the two native populations cancel.
The user trained despite the CI including zero, since no push was left for other terms to absorb
(round 3's failure mode; the local glpG TM4 test that also supported it used pre-fix inputs, 3.11):
ff30_bio with the library frozen, and **ff30_gdepth** with glycine's alpha_R / alpha_L depths
trainable on ff2.1's own library (one pooled pair on the GLY|X maps, started at BioEmu's weights:
c_aR -0.0918, c_aL +0.4621, dL - dR +0.554), so no data outside the training set enter and BioEmu is
the external check. The gdepth line became ff3.0 (1.26-1.27).

### 1.20 Upside's poly-glycine is chiral; the all-atom chain cannot be (2026-10-05)

A chain of glycines with neutral caps has no chiral residue, so in any physical force field each
glycine's alpha_R and alpha_L populations are equal, and so are beta / beta' and pPII / pPII'. The
all-atom Ac-(Gly)20-NHMe reference (plan.md Phase 12) was stopped by the user at ~97 ns per replica,
its data kept. Upside G20 (`scratchpad/polygly/`: `upside_run.py`, ibi_run.py's recipe at T_up = 1,
8 seeds x 200,000 tu from the fully extended start, the first 10% dropped; `polygly_analysis.py`, the
training's mirror-exact basins, residues 2-19), SE over seeds:

| | alpha_R | alpha_L | ln(aR/aL) | ln(beta/beta') | ln(pPII/pPII') | >= 4 alpha_L in a row | Rg (N, CA, C) |
|---|---|---|---|---|---|---|---|
| ff2.1 (NDRD) | 0.138 | 0.283 | -0.720 +- 0.010 | +0.157 +- 0.004 | +0.139 +- 0.007 | 0.058 of frames | 9.47 A |
| bio_start (BioEmu library) | 0.135 | 0.211 | -0.446 +- 0.014 | -0.020 +- 0.005 | +0.069 +- 0.003 | 0.030 | 9.81 A |

Runs of >= 4 alpha_R occur in 0.008 of frames under both. The handedness is uniform along the chain
(interior and end residues within 0.04) and each seed's halves agree, so it is not an end effect or
drift. Health: all coordinates finite, KE/1.5kT 1.002-1.005, backbone bond spread 0.141 A against
equipartition's sqrt(T/k) = 0.144 A for Upside's k = 48 springs. The BioEmu library removes most of
NDRD's beta and pPII asymmetry but keeps an alpha_L excess near the mean of its G-G-Y and X-G-G fits
(-0.19, -0.87; 1.19). Which terms carry the chirality is not measured: the GLY|GLY map entries, the
reference-state correction (built from a free alanine chain, applied to glycine too; Jumper thesis
4.3) and any CB-dependent placement acting on glycine are candidates. A per-term energy difference
between frames and their mirror images would separate them.

---

### 1.21 ConDiv workers destroy free replicas at dt 0.015; the round-4 panels; the side-chain verdict (2026-10-06/21)

**The event.** At dt 0.015 about 0.14% of worker protein-steps destroy a free replica. ff30_gdepth step 12
(epoch_00_minibatch_12): protein 3jtz ended with replicas 11 and 12 at `avg_kinetic_energy/1.5kT` 59.3 and 176.9,
its other replicas at 1.007-1.029. Replica 8 was at potential 2164 already at t = 20, in the anneal at T 0.05
(start -75), and the highest replica potential rose 2e3, 3e3, 1e4, 9e4 over t 20-70 with T at most 0.055. That is
not thermal: exchange only permutes and pivots reject uphill moves at that T, so the energy enters through
integration. Breakdowns begin in the anneal (T 0.04-0.05, potentials 2e3-2e5) or mid-run (T 0.80-1.04), once in r0,
the native-restrained replica (3f5r step 33), which fits integration better than conformational wandering. Replica
exchange carries the configuration to r12, the hottest free replica (T 1.10), as in glpG (§3, "exchange carries the
wrecked replica"). It stays finite (no H-bonds, Rg 25-34 A, potential ~1e4), so `check_step.py` passes it
("non-finite or missing: none") and the step's divergence includes its frames; only the KE line's maximum shows it,
and the watch scans that maximum for every protein in every new step. The worker deletes its `.up` files; the
step's inputs (`nesterov_temp__*`, `rama_round_00.dat`) survive in the minibatch directory.
- **Rate.** Every dt 0.015 FF2 run combined: 20 in about 13,800 protein-steps (ff30 9/5311, basin 5/2736, gly
  2/1584, gdepth 4/792, glyhb 0/1560, bio 0/863, fixedpoint 0/599, glyprobe 0/383), consistent with one shared rate
  (gdepth the high tail, p 0.03). FF1-form runs at dt 0.009 (ff31-gly, gly-sym, gly-ctx, ff21-restart): 0 in 13,714
  (at 0.14%, p ~1e-9; ~1e-4 allowing for twice the duration). The protocols differ in several things, so this
  places the cause in the FF2 protocol or terms without isolating dt.
- **Not caused by our code or library modifications** (three read-only comparisons,
  `/project/trsosnic/yinhan/checks/broken_replica_20261006/` `agentA`, `agentB`, `agentC`). Every engine path active
  in a ConDiv worker is byte-identical to master or force-neutral there; the swap-momentum rescaling by
  sqrt(T_dest/T_src) (`main.cpp:477-488`, since 08-05) decides where a broken configuration's heat goes, not the
  onset (7/4320 against 4/3096 on the pre- and post-10-02 builds, p 0.77). `main_worker` is identical in all deployed
  FF2 trainers; against the FF1 original (`~/Documents/ConDiv`, a Dec 2025 FF1-form copy, not a git repository) it
  differs in dt (0.015 against 0.009), anneal (from 0.05 T, ramped to T over t 96-400), replica interval (5 against
  10), systems (14 against 8), duration (8000 against 4000 tu) and the FF2 terms; the 0.015 and the anneal came with
  the FF2 trainer on 09-24 as "the port's schedule", which Kleinmann's port, not available locally, cannot confirm.
  The libraries change any residue's rama gradient by at most 2.5 E_up/rad against NDRD, and BioEmu's steeper glycine
  alpha_L edge (grad 77, curvature 1014) had no events in ff30_bio; a map spans about 20 E_up and cannot hold the
  1e3-1e5 seen, and at dt 0.015 its stiffest wall gives omega*dt 1.1, below Verlet's limit of 2.
- **The time step was the integration problem.** On 10-06 (user) both round-4 runs restarted from ff2.1 at dt 0.009
  as `ff30_bio_dt009` and `ff30_gdepth_dt009`, with the control `ff21_ctrl_dt009` on the same trainer. At dt 0.015
  every replica ran 1-2% hot through integration error (KE/1.5kT median 1.017, maximum 1.03-1.06 in every healthy
  step); at 0.009, step 0 read min 0.996, median 1.005, max 1.014. The test, set beforehand, was no destroyed free
  replica in about 1,600 protein-steps: by 10-07 the three runs had 2,136 protein-steps (29,904 replicas) at dt
  0.009 with none, KE/1.5kT at most 1.017, where 0.14% predicts 3.0 (Poisson probability of none 0.05).
- **dt does not set the helix drift of 1.22.** The dt twins change order in the H-bond margin (ff30_bio_dt009
  +0.116 at step 9 against ff30_bio +0.157 and ff30_gdepth +0.128; ff30_gdepth_dt009 +0.165 at step 5 against its
  twin's +0.136), so run-to-run variation (about 0.04) is at least as large as any dt effect; on the epoch-0 panel
  b9_00 folds 0.478, helix -0.031, against b00's 0.465, -0.033.

**The round-4 panels** (the 44-domain all-atom panel, `ff3_selection/select_*_login1_<date>.txt`; folded fraction,
then each class's error; every run's per-epoch folded fraction on all 44 domains is in 1.37's table):
- **Epoch 0 carries the loss.** Both dt-0.009 designs lose about 0.125 of folded frames in epoch 0 (bio_start 0.602
  to b9_00 0.478, gdepth_start 0.603 to d9_00 0.477) and 0.013-0.020 in each of epochs 1 and 2, while helix error
  stays at -0.023 to -0.030 and the start dominates each end in helix. gdepth's margin held at +0.14 through epoch 0
  while bio's fell below +0.10 by step 14, with the same folded loss, so the margin does not set it. No glycine class is
  resolved from the starts (gly_helix SE 0.05-0.07, gly_left 0.02-0.03). The panel is solvated, so the
  side-chain/lipid table (4.15) does not enter it.
- **ff_3.0_gdepth (`d9_03`, step 75):** folded 0.464, helix -0.030, beta -0.021, gly_helix -0.121, gly_left -0.008
  on 41 domains, against gdepth_start's 0.603, -0.018, -0.013, -0.059, +0.002, which dominates it in helix, beta and
  helical glycine (`select_last_d9_03_login1_20261010.txt`). The three round-4 releases on their 36 shared domains
  (`select_last_b9_03_bs_03_d9_03_login1_20261011.txt`): folded b9_03 (ff_3.0_bio) 0.452, bs_03 (ff_3.0_bio_si)
  0.467, d9_03 0.464; helix -0.030, -0.025, -0.030; beta -0.020, -0.014, -0.023; bs_03 dominates d9_03 in helix and
  beta and is worse than ff2.1 in helix.
- **Slower or frozen H-bond and sheet training loses more folding than the port rates at the same step.** These
  runs train the side chains at the port rate (rot 0.0125) and the H-bond and sheet terms at half the rate (`_si`)
  or not at all (`_fz`). bz_00 (0.450 against b9_00's 0.478), bs_00, bs_01 and ds_00 fold 0.028, 0.021, 0.035 and
  0.065 less than their port-rate twins, dz_00 0.432 against d9_00's 0.477; ds_00 lost more than d9_00 at nearly the
  same margin (+0.136 against +0.141). The frozen run's gap closes at epoch 1 (bz_01 0.461 against b9_01's 0.465);
  the SI-rate gap stays (bs against b9 0.021, 0.035, 0.019 at steps 18-56). Folded fraction carries no error bar in
  the panel, so this is a consistent direction, not a test.

**The side-chain verdict: BP against no-BP** (user 10-17: ff_2.1_nobp, the deterministic side chain, against its
control ff_2.1_ctrl38 on held-out proteins; detail in `$P/training/ff21_nobp_dt009/setup/notes/findings.md`).
- **The 44-domain panel does not resolve it** (`nb_01`, `ctrl38`, step 37; the panel runs every tag at time step
  0.015): folded 0.389 against 0.405 (-0.016 [-0.098, +0.065], p 0.62); mean Q, CA-RMSD (+0.04 [-0.19, +0.28] A) and
  the local classes resolve nothing. `ctrl38` reproduces `c9_01` bitwise in all 176 runs, so the panel is
  deterministic.
- **mdCATH, 89 held-out domains at dt 0.009** (`mdc_nobp`, `mdc_ctrl38`; plan.md Phase 19;
  `ff3_selection/panel_compare_49412845.out`): folded 0.264 against 0.365, -0.101 [-0.155, -0.051], Wilcoxon p
  7e-5, lower in 59 domains and higher in 19 (sign test p 6e-6); mean Q -0.035 [-0.046, -0.024] (p 2e-7); CA-RMSD
  +0.21 [+0.03, +0.38] A (p 0.004); lower in every length bin and fold class, the largest losses whole folds (1krrA00
  0.05 against 0.93). This agrees with the training set's +0.67 A. Whether the 44-domain panel missed it through its
  proteins or its time step is not separated.
- **Audit** (10-18, plan.md Phase 20; `checks/nobp_audit_20261018/`): the runs differ only in the intended side
  chain. Every no-BP config (training step 0, 32 benchmark arms, 356 panel runs) is single-state with 1-body 0 and
  burial weight 1, every control config multi-state; the engine at n_rot 1 is a plain sum whose parameter
  derivatives match central differences to float precision; solver, initial vector, minibatches, invocation (dt
  0.009, ladder 0.8-1.1, 8000 tu), code and binary are identical; no protein-step above KE/1.5kT 1.05; releases
  equal their step-37 files; the mdCATH numbers recompute exactly from the per-run files. But no-BP was not at a
  fixed point at step 38: `convergence_gate.py` fails its `rot` group (pair table, coverage and hydrophobe; p <
  4e-6, the test's floor) where the control's passes (p 1.000), its gradient pointing the same way in both epochs
  (cosine +0.81). The bead geometry is not trained; at the 1aba native the kept states cost +91.9 (MET0 +23.6,
  LEU23-PHE30 +8.3) against BP's -75.4, which fits a pull to soften clashes the fixed geometry makes; that the pull
  is this is not shown. The audit also found the shared angular-gradient defect (3.12).
- **At the two fixed points** (both continued from step 38 on the fixed engine; at step 76 every group passes,
  `ff21_nobp_dt009/gate_step76.txt`, `ff21_ctrl_dt009/gate_step76.txt`; `ff3_selection/panel_compare_49435922.out`):
  mdCATH folded 0.241 against 0.353, -0.111 [-0.160, -0.065] (Wilcoxon p 1e-5), lower in 54 domains, higher in 18,
  equal in 17; mean Q -0.030 [-0.040, -0.019] (p 1e-6); CA-RMSD +0.26 [+0.12, +0.40] A (p 6e-4). Folded falls in
  every length bin (50-100 -0.087, 101-150 -0.109, 151-250 -0.137) and in helix (-0.151) and mixed (-0.135) folds;
  beta is not resolved (-0.034 [-0.097, +0.038], p 0.13). Neither run's 38 fixed-engine steps changed its own
  folding, and no-BP is still not at a fixed point in ~110 angular pair coefficients (1.36).
- **On TM-score** (user 10-21: one metric for every test; plan.md Phase 26; `checks/converge_20261020/tm_panel.py`
  with the Peng scorer's `tmscore.py`, every 10th frame after 2,000 tu, 4 seeds pooled): mdCATH at the fixed points
  0.667 against 0.702, -0.034 [-0.045, -0.024], p 7e-8, lower on 65 of 89, in every length bin and class (beta
  -0.013, p 0.05); the 44 Charron domains, both fixed points rerun at dt 0.009 (`p44_nobp_03`, `p44_ctrl_03`), 0.623
  against 0.650, -0.026 [-0.043, -0.010], p 0.007, lower on 29 of 44 (Q folded fraction 0.363 against 0.430, p
  0.06); Peng's benchmark at step 37, all 16 pairs, natives 0.457 against 0.486 (-0.029 [-0.070, +0.015], p 0.18),
  de novo 0.309 against 0.299 (p 0.49). The share of frames at TM > 0.5 does not differ (0.920 against 0.946, p
  0.24; 0.824 against 0.832, p 0.74): without BP the native fold is kept as often but sits further from the native,
  which the Q-based folded fraction, a stricter cutoff, reads as a larger loss.
- **Verdict: the deterministic side chain loses folding on held-out proteins,** at step 37 and at both fixed points.
- **Both side-chain models sit below Peng's published FF2 on his benchmark** (Peng SI Fig S5: FF2 <TM> 0.55 native,
  0.42 de novo; FF1 0.46, 0.37): control 0.486 / 0.299, no-BP 0.457 / 0.309 at step 37. The scoring is not the
  cause: the same scorer gives Peng's own FF1 trajectories 0.48 / 0.36 (11c). Against Fig S5's FF2 curves read by
  eye (means 0.53 / 0.41), the native gap is two proteins that unfold completely at their cold rungs: alpha3d (T
  0.90; control 0.17, median CA-RMSD 23 A, no-BP 0.35; FF2 ~0.62, retired ff3.0 0.61) and homeodomain (T 0.87; 0.20
  and 0.19; FF2 ~0.51, retired ff3.0 0.54); the other 14 agree with FF2 to a mean difference of ~0.00. De novo,
  alpha3d, top7, homeodomain, NuG2, protein G and protein L fail to fold where FF2 reaches 0.4-0.6. Both side-chain
  models share the deficit, so it is not the side chain; it matches the folding the retraining lost (1.37). The
  paired side-chain comparison stands; the absolute values are not FF2's.
- **The side-chain model moves the shared H-bond margin, in opposite directions** (`check_step.py`, read 10-18):
  E_other - E_alpha at steps 18 / 37 is no-BP +0.296 / +0.354 and the control +0.135 / +0.157 (start +0.192), a gap
  of ~0.2 against run-to-run variation of ~0.04. No-BP's rise is E_other weakening (-1.769 -> -1.620), so an H-bond
  made at phi > 0 costs ~0.15 E_up more, E_alpha nearly fixed; the control's fall is E_alpha weakening (-1.961 ->
  -1.854). On the fixed engine no-BP read +0.289 to +0.357 over steps 38-75. Whether the high margin (left-handed
  natives, glycine alpha_L and some Asn/Asp, pay more per H-bond) is part of no-BP's mdCATH loss is not tested; the
  margin has not predicted helix populations or panel folding before (1.15, the panels above, 1.24).

### 1.22 Round 4 after one epoch: the panels lose folding as round 3's did (2026-10-06)

The dt 0.015 runs ff30_bio (`b00`, `b01`) and ff30_gdepth (`d00`) were replaced by their dt 0.009 restarts (1.21).
Epoch-0 panel, `select ... ff21_released ff21_awh bio_start gdepth_start b00 d00`, 42 domains (3g7lA00 and 4hwiB01
dropped for too few folded frames):

| tag | folded | helix | beta | gly_helix | gly_left |
|---|---|---|---|---|---|
| ff21_released | 0.615 | -0.018 +- 0.004 | -0.018 +- 0.005 | -0.081 +- 0.048 | +0.034 +- 0.019 |
| bio_start | 0.602 | -0.022 +- 0.005 | -0.017 +- 0.004 | -0.050 +- 0.041 | +0.005 +- 0.020 |
| gdepth_start | 0.603 | -0.019 +- 0.005 | -0.015 +- 0.004 | -0.041 +- 0.044 | +0.013 +- 0.019 |
| b00 | 0.465 | -0.032 +- 0.006 | -0.018 +- 0.004 | -0.096 +- 0.045 | -0.018 +- 0.023 |
| d00 | 0.471 | -0.026 +- 0.005 | -0.017 +- 0.004 | -0.046 +- 0.038 | +0.009 +- 0.021 |

* **Each start dominates its own epoch 0 in helix** (paired bootstrap); neither glycine class is resolved. The
  folded fraction falls from 0.60 to 0.47 in both runs, as round 3's h00 and h01 fell to 0.460 and 0.469 (1.17):
  frozen BioEmu (b00) and trained depth (d00) end alike, and epoch 1 continues the loss (`b01` 0.400 on 38 domains).
* **The matched control loses at least as much** (`c9_00`: ff21_ctrl_dt009, ff2.1's own library, otherwise
  ff30_bio_dt009): folded 0.390, helix -0.032 on 35 domains, against b9_00 0.478 and d9_00 0.477; d9_00 dominates it
  in helix and beta; its second epoch `c9_01` folds 0.405 (`select_c9_01_login1_20261007.txt`). So the soluble
  panel's folding loss comes with the training itself, not with the glycine change. This is the panel, not TM4
  (10.19).
* **ff30_gdepth's first depth round barely moved:** dL - dR +0.554 to +0.549, free-native gap +0.005. The training
  data put almost no push on the BioEmu-weighted start.
* **The shared H-bond margin E_other - E_alpha** crossed +0.10 in both dt 0.015 runs within the first three steps of
  epoch 1 (ff30_bio +0.088 at step 23, ff30_gdepth +0.099 at step 21), mostly through E_alpha weakening (ff30_bio
  -1.961 to -1.862). With no glycine change it stayed at +0.127 to +0.142 (ff21-fixedpoint, steps 13-25, 1.17) and
  at +0.125 to +0.157 through epoch 1 (ff21_ctrl_dt009), while ff30_gdepth_dt009 held +0.14 through step 22 and
  ff30_bio_dt009 fell to +0.050 at step 37; the BioEmu run's extra drop at step 14 (0.036) is at the edge of
  run-to-run variation. Round 1's symmetrised map held TM4 under the FF1-form port, which trained no H-bond energy
  (9e: the port dropped `hb` and `sheet`). From round 2 the FF2 trainer fits E_other, shared by all 20 types and
  near-degenerate with glycine's alpha_L depth (9e), so a map without the selection leaves the training data a
  second route to it.
* **glpG TM4, local test: every count of 10-04 to 10-07 is invalid** (3.11): each run paired its force field with
  the FF1-form ff_3.0's coverage tables. The record is kept only as such (Mac `runs_precov_20261007/`; cluster
  `checks/r4_epochs/tm4_local/runs/`), and every reading drawn from it is withdrawn. The counting design stays
  valid. TM4 is lost seed by seed: a seed either holds 0.95-1.00 or unwinds to 0.55-0.89, so a set reports how many
  seeds unwound (last-block TM4 below 0.90) and how many flipped (a TM4 glycine at phi > 0 in more than 0.25 of the
  last block), criteria fixed before the 12-seed runs. By Fisher's exact test, 3 seeds resolve only one in three
  against three in three; 12 seeds resolve one in four against three in four (p 0.04); one in three against two in
  three needs 24 seeds (p 0.04; 12 give 0.22).
* **Energy jumps in the local glpG runs** (every run log; `events_vs_tm4.py`): between frames 10 tu apart the total
  potential normally changes by 80-150, but 13 of 105 pre-fix seeds show single jumps of 3,000-12,000, under every
  force field (ff21_released, bio_start and ff21_awh included). The worst, `fp_e00` seed 11 at t = 3630 (-22,510 to
  +17,828 in 30 tu, H-bonds 189 to 61, Rg 20.4 to 23.9 A), leaves a last block that is not a TM4 measurement. Six
  of the 45 unwound seeds carry a jump, so jumps do not drive the unwound counts; no frames were dropped. These runs
  are dt 0.009 Verlet in the dry-MARTINI bilayer, a different simulation from ConDiv's workers (1.21). The
  fixed-input sets of 1.24 show such jumps too; all ran on the old side-chain/lipid table, whose single-evaluation
  kicks are in 4.10, and the fixed-table reruns of 4.15 report none above 3000.

### 1.23 How Jumper and Peng trained the H-bond energy (read 2026-10-07)

Sources: Jumper's thesis (`~/OneDrive.../Jumper-thesis-final-Jan2017.pdf`, printed pages), the ConDiv and
side-chain papers (PLoS Comput Biol 2018), Peng's 2022 SI, and every trainer on disk or on midway2 (paths in memory
`condiv-ff2-trainer-source`).

**No soluble trainer froze the H-bond energy, or trained it in a separate stage.** In each, `hb` (with `sheet`, and
`dhb` from FF2 on) is a group in the same Adam solver as every other group, at its own fixed rate. The only stage
that applies to all groups is Jumper's fine-tuning: every rate x 0.25 after two epochs, run as a restart
(`test37_fixbead_finetune`).
* **Before ConDiv the H-bond was hand-set:** scanned from -2.4 to -1.5 and set to -1.8 kT for the side-chain model,
  "the only parameter in the model directly optimized for simulation accuracy" (thesis p. 30-31; side-chain paper
  p. 19). ConDiv started from there.
* **Jumper's ConDiv trained it jointly** ("a single parameter that is chosen by contrastive divergence", p. 47), on
  purpose: fixing known interactions to experiment "is inadvisable" (p. 66), and training all terms together
  prevents "hydrogen bond terms overwhelming the side chain interactions to make very long helices" (p. 43).
* **Peng's FF2 trained them jointly too.** SI p. 2: "ConDiv training return an H-bond energy for helix, strand and
  turn of -1.96, -1.95, and -1.77", and the second H-bond -0.41 "obtained using ConDiv training"; these are ff2.1's
  `hbond.h5` entries 0-3. The SI does not say how the 20 `sheet` values were obtained.
* **The only frozen soluble H-bond is in Peng's membrane trainer** (`~/Documents/Train/ConDiv.py`), where the whole
  soluble force field (`ff_2.2`) is fixed and only `membrane.h5` trains.
* **ff_2.0 to ff_2.1 changed only burial.** `hbond.h5`, `sheet` and `bb_env.dat` are byte-identical, `sidechain.h5`
  differs by at most 3e-4 (a ConDiv run at rot's rate would move it far more), and `environment.h5` by up to 2.8 in
  center, 1.3 in scale, 0.95 in weights. How ff2.1 was made is not documented.

**Where they stopped: a fixed budget, short of convergence on purpose** (read 10-21, user: "check Peng's paper on
how he determined where to stop").
* Peng states the ideal and runs a fixed budget. Main text p. 2: the iteration "continues until no energy parameter
  can be updated to produce a better NSE", and "to avoid overtraining, this procedure is run for 456 proteins in
  batches of 19". SI p. 4: "a total of 76 ConDiv iterations (4 cycles) were executed", at the initial rate of Ref
  [2]'s fine-tuning stage. No convergence test, no validation-based choice of checkpoint and no starting point for
  FF2's training are given; the 16 proteins are scored after training.
* Jumper (Ref [2], PLoS Comput Biol 2018, p. 7; thesis pp. 51-52, 66) ran two epochs at the large step, cut it 4x,
  and stopped early knowingly: "there are indications that the parameters have not yet converged. Earlier tests,
  however, showed that continuing the contrastive divergence until convergence does not necessarily produce better
  results". His mechanism: when large barriers surround the native states little relaxation occurs, and further
  fine-tuning "overtrains the model to distinguish between native and near-native structure at the expense of
  distinguishing against a more diverse ensemble"; early stopping acts as a regulariser (Duvenaud et al. 2016).
  Thesis p. 66: contrastive divergence "is noted to obtain peak accuracy before fully converging then declining
  accuracy afterwards"; his proposed remedy is persistent contrastive divergence.
* So our convergence gate (1.21, 1.26) is our own addition, against both authors' practice, and every run started
  from ff2.1, which had already had FF2's full budget: each of our 76-step runs is a second budget on top of it. The
  first-epoch loss (1.37) is the decline Jumper describes, and it fits his mechanism: the training proteins' free
  replicas barely relax (RMSD flat over four epochs, 1.37).

**What kept the H-bond in check was rate and early stopping, not freezing.** Effective rates (Adam moves a scalar by
about its rate per step when its gradient is consistent; Kleinmann's log shows hb[0] -1.9609 to -1.9709 on step 1):

| trainer | rot | env | hb | dhb | sheet | DSE lambda | length |
|---|---|---|---|---|---|---|---|
| Jumper FF1, large step (thesis p. 67) | 0.5 | 0.1 | 0.02 | - | 0.03 | none | 2 epochs |
| Jumper FF1, fine-tuning (x 0.25) | 0.125 | 0.025 | 0.005 | - | 0.0075 | none | to ~200 steps, stopped early on purpose (p. 51) |
| Peng FF2, per SI p. 4 ("the initial learning rate of the fine-tuning stage in Ref [2]") | 0.125 | 0.025 | 0.005 | ? | 0.0075 | 0.3 | 76 iterations (4 cycles), no convergence test |
| Peng 2019 intermediate (`upside-pxd/ConDiv`) | 0.125 | 0.05 | 0.025 | - | 0.001 (helix/strand/turn) | none | ? |
| Kleinmann port (`condiv2.py`, x 0.5) | 0.125 | 0.05 | 0.01 | 0.005 | 0.015 | 0.0 | ~101 steps |
| this trainer (`training/ConDiv.py:695`, x 0.5, rot / 10) | 0.0125 | 0.05 | 0.01 | 0.005 | 0.015 | 0.3 | to the gate, up to 13 epochs |

The FF2 row assumes Peng's per-group bases were Jumper's: his own file is lost (9v). Both ports share Kleinmann's
x 0.5, so the hb, dhb, sheet and env rates here are twice what the SI implies; rot alone agrees with the SI
(Kleinmann's base 0.25 x 0.5 = 0.125), before our own / 10. A rate sets how fast a term drifts, not where it settles.

**The direction of the H-bond drift follows the DSE term.**
* With no DSE, H-bonds got stronger. Jumper: -1.8 to about -2.24 by step ~75 at the large step, then back to about
  -2.11 in fine-tuning (thesis Fig. 3.2, p. 50; his "curious behavior", p. 52). Kleinmann at lambda 0: E_alpha
  -1.961 to -2.052, E_beta -1.946 to -2.046, E_other -1.769 to -1.836, dhb -0.406 to -0.320, over ~101 steps
  (`condiv_training_results.output`).
* With DSE at 0.3 they get weaker. Peng introduced the term as "a regularizer that simultaneously makes the residual
  interactions as weak as possible" (SI p. 4); the SARW reference has no H-bonds, so its gradient on each branch
  energy can only weaken it. This matches the gradient split (1.17: E_alpha NSE +21.45 against 0.3 DSE -40.54) and
  the margin fall in every round-4 run.
* Neither author reports helices being too weak. Their failure mode was over-stable helices (the `sheet` mixing
  exists "to counteract an observed tendency for our model to overstabilize helices", ConDiv p. 5). The training set
  had no membrane proteins.

**Consequence for the frozen runs.** `ff30_*_fz` freezing hb, dhb, hbg and sheet has no precedent in Jumper's or
Peng's soluble training. Its precedents are the FF1-form port, which dropped `hb` and `sheet` by accident (9e), and
the membrane trainer.

### 1.24 glpG TM4 on coverage-fixed inputs (2026-10-07)

The first TM4 counts on inputs patched by the fixed `patch_glpg.py` (md5 `70589119...`, 3.11). Every set ran on the
old side-chain/lipid table, and the fixed-table reruns of 4.15 replace their readings: on the fixed table ff2.1
itself holds TM4 (4.15, 1.35). The method and what the sets show about resolution stand.
- **Setup.** Each set is 12 seeds (24 where marked) of hybrid dry-MARTINI glpG (79HIS): Verlet dt 0.009, T 0.80,
  4000 tu, on the MacBook Pro or the Mac Studio, whose patched inputs agree in all 150 `/input` datasets and whose
  runs agree frame for frame. Runs in `checks/r4_epochs/tm4_local/runs_cov/`, tables in
  `tm4_local/tm4_compare_cov_<set>.txt` (`tm4_local.py`, `tm4_compare.py`).
- **Counts.** A seed is unwound if its last-block TM4 helix fraction, which scores each residue by its phi/psi
  alone, is below 0.90; flipped if GLY136, 143 or 149 has phi > 0 in more than 0.25 of the last block. The jump
  scan flags a total-potential change above 3000 between frames.
- **Test.** Each checkpoint is compared with its own start (10.19): until 10-09 by Fisher's exact test on the
  counts; from 10-09 by DSSP alpha-helix (1.25), the counts kept as secondary tests.
- **What enters the hybrid input.** A frozen-H-bond (`_fz`) checkpoint's patched input has its start's H-bond
  energy and Rama map exactly and differs from it only in three datasets from the trained `sidechain.h5` (the
  rotamer pair interactions and the two coverage tables); the trained `environment.h5` and `bb_env.dat` do not enter
  the hybrid input.

| set | checkpoint | seeds | unwound | flipped | dihedral box | DSSP alpha (1.25) |
|---|---|---|---|---|---|---|
| `ff21_bioT1_6` (bio_start) | start of ff30_bio_dt009, _fz, _si | 12 | 5 | 4 | 0.899 | 0.769 |
| `b9_00` | ff30_bio_dt009 step 18 | 24 | 7 | 5 | 0.914 | 0.815 |
| `b9_01` | ff30_bio_dt009 step 37 | 12 | 6 | 5 | 0.833 | 0.711 |
| `b9_02` | ff30_bio_dt009 step 56 | 12 | 8 | 4 | 0.830 | 0.694 |
| `bz_00` | ff30_bio_fz step 18 | 24 | 8 | 4 | 0.897 | 0.794 |
| `bz_01` | ff30_bio_fz step 37 | 12 | 6 | 2 | 0.846 | 0.739 |
| `bz_02` | ff30_bio_fz step 56 | 12 | 4 | 3 | 0.910 | 0.856 |
| `bs_00` | ff30_bio_si step 18 | 12 | 6 | 5 | 0.846 | 0.684 |
| `bs_01` | ff30_bio_si step 37 | 12 | 5 | 2 | 0.879 | 0.760 |
| `bs_02` | ff30_bio_si step 56 | 12 | 5 | 2 | 0.878 | 0.684 |
| `bs_03` | ff30_bio_si step 75 (ff_3.0_bio_si, 1.31) | 12 | 4 | 2 | | 0.739 |
| `gdepth_start` | start of ff30_gdepth_dt009, _si, _fz | 12 | 7 | 2 | 0.849 | 0.657 |
| `d9_00` | ff30_gdepth_dt009 step 18 | 12 | 5 | 2 | 0.900 | 0.801 |
| `d9_01` | ff30_gdepth_dt009 step 37 | 12 | 6 | 4 | 0.827 | 0.615 |
| `d9_02` | ff30_gdepth_dt009 step 56 | 12 | 4 | 3 | 0.899 | 0.765 |
| `d9_03` | ff30_gdepth_dt009 step 75 (ff_3.0_gdepth, 1.26) | 12 | 5 | 4 | 0.882 | 0.726 |
| `dz_00` | ff30_gdepth_fz step 18 | 12 | 4 | 5 | 0.912 | 0.791 |
| `dz_01` | ff30_gdepth_fz step 37 | 12 | 5 | 3 | 0.911 | 0.780 |
| `ds_00` | ff30_gdepth_si step 18 | 12 | 9 | 3 | 0.775 | 0.572 |
| `ff21_released` (ff2.1) | start of ff21_ctrl_dt009 | 12 | 9 | 5 | 0.827 | 0.610 |
| `c9_00` | ff21_ctrl_dt009 step 18 | 24 | 16 | 9 | 0.843 | 0.660 |
| `c9_01` | ff21_ctrl_dt009 step 37 (final) | 12 | 10 | 6 | 0.783 | 0.595 |
| `fp_e00` | ff21-fixedpoint step 18 (dt 0.015) | 12 | 5 | 4 | 0.902 | 0.730 |

Every set but `d9_02`, `d9_03` and `bs_02` carries 1-5 seeds with a jump above 3000, up to 201861 (`bs_03` s1,
1.31), in seeds that keep TM4 at 1.00 and in seeds that unwind alike; KE/1.5kT is 1.001-1.023 in every seed except
ds_00 s7 (1.080, after a jump of 9732 at t 3776).
- **No checkpoint is resolved from its start by the counts.** The one nominal result, bz_00's seeds 1-12 flipping 0
  against 4 (one-sided p 0.047), came from one of about 30 uncorrected count tests, and its seeds 13-24 flipped 4 of
  12, bio_start's own count. The continuous measures (last-block TM4 mean, TM1) were added tests, uncorrected; TM1
  has no glycine.
- **Two 12-seed halves of one force field differed by seven flips.** c9_00's seeds 1-12 and 13-24 ran the same
  patched input on the same binary with the same flags, only `--seed` differing, and flipped 1 and 8 of 12
  (two-sided Fisher p 0.009, the most extreme of three comparisons that were not planned; b9_00's halves flipped 4
  and 1, bz_00's 0 and 4). One set's flip count is a noisy reading of its force field; a lean seen in one 12-seed
  set is a reason to run 24 seeds, not a result.
- **The count depends on the cut** where seeds end near 0.90: seven of ff2.1's nine unwound seeds end at 0.85-0.89,
  and at a cut of 0.85 ff2.1 and bio_start both count 3.
- **The TM4 lean does not follow the H-bond margin.** bs_00 leaned worse than b9_00 with the higher margin (+0.115
  against below +0.10 since step 14), and c9_01 leaned below c9_00 with c9's margin at +0.125 to +0.157 through
  epoch 1, as b9_01 did with b9's at +0.050.
- **The size of difference a 12-seed set can resolve.** With bio_start at 5 of 12 unwound, a one-sided Fisher test
  reaches p < 0.05 only if a checkpoint unwinds 10 or more of 12 (worse) or none (better); for flips (4 of 12) 9 or
  more, or none. Against ff21_released's 9 of 12 the controls' checkpoints can resolve only as better, at 4 or fewer.
  A one- or two-seed difference cannot be resolved at 12 seeds per set.

### 1.25 TM4 secondary structure by DSSP (2026-10-09)

The user asked (10-09) that the TM4 test measure the stability of TM4's secondary structure. 1.24's helix fraction
and unwound count score each residue's phi/psi in a box (phi -130 to -20, psi -90 to 15), which a residue can
satisfy while its backbone H-bonds are broken.
- **What is added.** `tm4_local.py` runs DSSP (mdtraj 1.11) on N, CA, C and the frame's O, which is Upside's own
  carbonyl from `infer_H_O`, written by the hybrid into the O slot (`src/martini_hybrid.cpp:86`); only residue
  210's slot is unused: it keeps its input coordinate, up to 18.9 A from its C (3.8). The window is 135-151, TM4's
  DSSP helix at t 0 and in the simulation PDB (§11); residue 134 is coil at t 0. Two readouts per time block:
  alpha-helix (DSSP H) and helix of any kind (H, G or I). `tm4_compare.py` prints both per seed, alpha per residue
  over the last block, and a one-sided Mann-Whitney on each. The pre-DSSP tables are in
  `tm4_local/tables_pre_dssp_20261009/`.
- **Where the two disagree, DSSP is right about the H-bonds.** 42 of 228 finished seeds hold the box at 0.90 or more
  in the last block but less than 0.90 alpha-helix by DSSP, and 19 less than 0.80. Of four examined, three
  (ff21_released s4, b9_00 s19, bz_00 s22) have 0.28-0.38 of TM4's last-block residues in pi-helix (DSSP I), where
  the i->i+4 O...N distance is 3.9-5.3 A against 2.6-3.1 A along an intact seed; in the fourth (bz_01 s2) 147-151
  unwind into turn and coil. In all four the phi/psi stay inside the box.
- **TM4 keeps less secondary structure than the box showed, in every set** (1.24's table): DSSP alpha is 0.10-0.22
  below the box, and helix of any kind is below it too (0.77 against 0.86 over all seeds); 36 seeds keep less than
  half of TM4 in any helix over the last block, against 10 below 0.50 by the box. The order of the sets is nearly the
  same (Spearman 0.88 for alpha and 0.90 for any helix against the box, over 228 seeds), so no direction read in
  1.24 reverses.
- **TM4 is still losing helix at the end of the run.** Alpha falls block by block in nearly every set (gdepth_start
  0.95, 0.86, 0.75, 0.66; ff2.1 0.91, 0.84, 0.70, 0.61). 4000 tu does not reach a plateau, so the readout is how
  fast TM4 loses its helix, not where it settles (4.4, "Compare trends, not means").
- **A seed loses a stretch, not a residue** (last-block DSSP codes per residue, `scripts/tm4_residue_ss.py`; TM4
  135-151 is T G V V Y A L M G Y V W L R G E R). In d9_02, s6 loses 140-150 whole, s7 143-151, s8 135-139 in part
  and 148-151, s5 139-145, while the other eight keep TM4 at 0.80-0.99, so the per-residue means are set by 2-4 seeds
  per set. On the old table ff2.1 loses the midplane (alpha 0.40-0.52 at 141-143, pi-helix 0.17-0.19 at 140-143,
  turn 0.22-0.25 at 142-147).
- **The primary test from 10-09 (user).** Each seed's DSSP alpha-helix fraction of 135-151 over the last block is
  compared, checkpoint against reference, by a two-sided Mann-Whitney test, resolved at p < 0.05, with a bootstrap
  95% interval on the difference of means (`tm4_compare.py`, "## primary test"). The dihedral counts of 1.24 remain
  as secondary tests.
- **No in-training checkpoint keeps TM4 more helical than its start, resolved;** every interval spans zero. Along
  each training (start, then steps 18, 37, 56, 75): ff30_bio_dt009 0.769, 0.815, 0.711, 0.694; ff30_bio_fz 0.769,
  0.794, 0.739, 0.856; ff30_bio_si 0.769, 0.684, 0.760, 0.684, 0.739; ff30_gdepth_dt009 0.657, 0.801, 0.615, 0.765,
  0.726; ff30_gdepth_fz 0.657, 0.791, 0.780; ff30_gdepth_si 0.657, 0.572; ff21_ctrl_dt009 0.610, 0.660, 0.595.
- **bz_02 against its twin b9_02 is the only resolved primary test, and only just**
  (`tm4_compare_cov_bz_02_vs_b9_02.txt`): 0.856 against 0.694, +0.162 [-0.010, +0.330], two-sided p 0.046. The
  bootstrap interval reaches below zero, and of the 32 primary-test tables read one or two would pass p < 0.05 with
  no true difference. The twin gaps open as the port-rate twin loses helix and its margin falls while the frozen
  twin's stays at +0.192: bz minus b9 -0.021, +0.028, +0.162 at steps 18, 37, 56; dz minus d9 -0.010, +0.165 (p
  0.069) at steps 18, 37.
- **What the test can resolve.** Seeds' last-block alpha spreads with SD 0.20-0.35 per set (pooled 0.28), and
  bimodally: a seed keeps TM4 near 1.0 or loses much of it. Resampling the pooled residuals, a two-sided
  Mann-Whitney at p < 0.05 detects a true difference of 0.10, 0.15, 0.20 and 0.30 with power 18%, 31%, 47% and 73%
  at 12 seeds per set, 34%, 58%, 80% and 96% at 24, and 56%, 86%, 97% and 100% at 48. The differences seen
  (0.01-0.15) are below what 12 seeds resolve; 0.10-0.15 needs about 48 seeds per set.

### 1.26 The first converged round-4 gate: ff30_gdepth_dt009 at step 76, released as ff3.0 (2026-10-09)

- **The gate** (`ff30_gdepth_dt009/gate_step76.txt`) read steps 58-76 and passed every group at p > 0.005: hb 0.95,
  sheet 1.00, rot 1.00, the environment groups 0.99-1.00, bbenve 0.79, and dhb lowest at 0.045. It released
  `epoch_03_minibatch_18` (d9_03) as `ff_3.0_gdepth`, with depth round 4's glycine library (offsets aR -0.1068 aL
  +0.5169, dL - dR +0.624 against the start's +0.554).
- **What converged is the held run.** The released H-bond margin is +0.082, under the +0.10 hold since step 33
  (start +0.192); dhb -0.551 (start -0.406), sheet mean 0.228 (start 0.161). The gate tests a fixed point, not the
  margin.
- **The user took it as ff_3.0** (10-09): `parameters/ff_3.0/`, md5-identical to the release. Its trained glycine
  library is the `rama.dat` in that directory.
- Its panel is in 1.21 (folded 0.464 against gdepth_start's 0.603). Its old-table TM4 (1.24-1.25) is not resolved
  from gdepth_start (DSSP alpha 0.726 against 0.657, p 0.64), d9_02 (0.765, p 1.00) or ff2.1 (0.610, p 0.25).
  ff2.1 and ff3.0 patched into the live 79HIS seed and run with `run_glpg.sh`'s flags reproduce the TM4-test sets'
  seed-1 logs frame for frame (all 401 frames).

### 1.27 ff3.0's glycine map is NDRD plus one depth pair (checked against the file, 2026-10-09)

Measured on the MacBook Pro; files are in its `scratchpad/glpg_ff21_vs_ff30/` unless named otherwise (copied to
midway2's `checks/tm1_hbmem_20261009/`, 4.6).
- `parameters/ff_3.0/rama.dat` (md5-identical to `$P/parameters/ff_3.0_gdepth`) differs from `common/rama.dat` only
  in the 37 GLY|X coil maps; the sheet group, every other map and the pooled GLY|ALL entry are NDRD's.
- All 37 carry one identical shift (to 1.2e-6). Fitted to the trainer's own basin weights
  (`ff30_gdepth_dt009/trainer/rama_basin.py`) it is alpha_R -0.1068, alpha_L +0.5169, the round-4 offsets of 1.26
  (`rama_rounds.txt`); the trainer's write rule reproduces every stored map to 1e-6.
- **A label value is not a depth.** The library renormalises each map and the config writer subtracts the
  Boltzmann-weighted mean (`upside_config.py:841`), so raising alpha_L moves every other cell by a constant: -0.155
  from renormalisation and -0.167 from the mean, -0.322 in all. On the 2,492 cells the offsets do not touch, ff3.0
  minus ff2.1 is -0.321 to -0.322. PII - alpha_L moves by exactly -0.517 and PII - alpha_R by +0.107.
- For figures, zero each map at the Boltzmann-weighted mean outside the trainer's helical basins (helix weight <
  1e-4, 2,072 cells): ff3.0 then equals ff2.1 there to 3e-5. On that zero, alpha_L / alpha_R / PII are -3.0 / -1.6 /
  -0.6 (ff2.1), -2.5 / -1.7 / -0.6 (ff3.0), -1.2 / -0.9 / -0.9 (BioEmu-fitted). Scripts
  `scratchpad/plot_gly_rama_ff30.py` and `plot_gly_rama_nonhelix_zero.py`; figures `gly_rama_ff30.png`,
  `*_nonhelix_zero.png` and the US Letter page `gly_rama_nonhelix_zero_page.pdf`, in the MacBook Pro's `~/Downloads`.

### 1.28 Training every basin depth (an ff3.1 question): the parameter count is not the limit (2026-10-09)

The user's question (10-09, after talking with Jumper): ff3.0 trains one glycine depth pair; could ff3.1 train the
relative depth of every basin of every map, given 456 proteins? Jumper said the parameter count is no obstacle,
citing deep-learning practice (regularization, early stopping).
- **The count is no obstacle, and Upside is the precedent.** ff2.1's ConDiv already trains 31,285 parameters with a
  gradient (rot 30,800 of them) on 456 proteins; Jumper et al. 2018 (PLoS Comput Biol 14:e1006578, Discussion):
  "very large numbers of parameters (even ~20000 in our case) are no obstacle". That paper trained one Ramachandran
  parameter (the sheet term) on fixed NDRD maps. The literature agrees with conditions: gradient descent stopped
  early is close to ridge toward its start (Yao et al. 2007; Ali, Kolter, Tibshirani 2019); double descent needs the
  stopping time or penalty set on held-out data (Belkin et al. 2019; Nakkiran et al. 2020, 2021); Potts/DCA models
  fit 1e5-1e7 couplings to far fewer sequences only with L2 (Ekeberg et al. 2013). The data count is sites, not
  proteins: 49,401 interior residues, 98,802 map reads.
- **A per-map depth model.** 800 coil maps (840 with cis-proline) of 72 x 72; with the trainer's six basins
  (alpha_R, alpha_L, beta, pPII and the two phi > 0 extended regions) that is 4,000 depths for coil alone, 2,400 with
  four basins. The engine's `rama_map_pot` has no parameter derivative, so a depth gradient comes from basin
  populations, as `rama_basin.py`'s rounds do.
- **The data per depth are thin.** Native sites per (map, basin) in the 410 training proteins have median 40
  (alpha_R), 23 (beta), 16 (pPII) and 1 (alpha_L; 649 of 800 maps have fewer than 5). Per-pair resolution (1.12):
  one epoch's per-offset step has split-half correlation 0.07, true per-pair corrections SD ~0.04 nat against 0.074
  noise, 70% of it the protein set, so more epochs remove little. Under a penalty set on held-out data most depths
  stay at NDRD; the effective parameters are the aggregates the data resolve (the GLY|X alpha_L class, per-central
  beta against pPII, which the `sheet` energies already carry).
- **Regularization does not touch the bias.** A depth trained so the free ensemble matches the native-restrained
  one converges to where evolution placed the residue (1.15-1.16, rule 1.8; round 2's alpha_L and the pre-proline
  offsets fitting composition, 1.14). Jumper's thesis (2017, section 4.3, PDF page 94): "it is imperative to limit
  the capacity of the Ramachandran potential to avoid capturing protein physics that should be described in the
  hydrogen bonding and side chain interactions"; free depths give ConDiv a route to stabilise native secondary
  structure past the H-bond and side-chain terms.
- **Two practical limits.** Adam normalises each coordinate's step, so a depth that receives only noise still moves
  at the learning rate: early stopping is ridge only for plain gradient descent (Wilson et al. 2017; an inference
  for ConDiv). The 46-protein held-out set sits at its sample-size floor (1.12), so it cannot show overfitting; a
  penalty or stopping time needs a validation signal that resolves (split-half reproducibility of the depths, a
  larger held-out set, the panels).
- Hierarchical pooling (a per-central mean plus penalised per-neighbour deviations, as NDRD and DCA do) is the
  standard shrinkage, but it is what rule 1.8 forbids for corrections; that is the user's call. Mode `pooled` (1.34)
  is this shrinkage.
- Sources checked 10-09 (literature agent; the thesis and paper quotes read directly): Zhang et al. ICLR 2017;
  Belkin et al. PNAS 2019; Nakkiran et al. ICLR 2020 and 2021; Bartlett et al. PNAS 2020; Yao et al. Constr Approx
  2007; Ali et al. AISTATS 2019; Li et al. AISTATS 2020; Prechelt 1998; Wilson et al. NeurIPS 2017; Jumper et al.
  PLoS Comput Biol 2018 (e1006342, e1006578); Jumper thesis 2017; Peng et al. JCTC 2022; Hinton UTML TR 2010-003;
  Ekeberg et al. PRE 2013; Morcos et al. PNAS 2011; Ting et al. PLoS Comput Biol 2010.

### 1.29 The ff3.1 probe on recorded data: per-map depths do not transfer between proteins (2026-10-09)

Route 1 (user, 10-09; plan.md Phase 14 step 1), measured on ff21_ctrl_dt009's recorded basin populations (c9:
ff2.1's NDRD maps, 38 steps, each of the 456 proteins once per epoch; 98,802 interior residue-steps) with no new
simulation. `scratchpad/ff31_probe/probe_depths.py`, output `probe_depths.txt`, `probe_global_by_epoch.txt`.
- **The model.** One offset per (coil map, basin), 840 maps x 6 basins of `rama_basin.py`; an interior residue reads
  its left and right map at share 1/2 each; MAP under a Gaussian prior of width sigma centred on NDRD, no pooling,
  solved jointly to first order in the populations' response. The width is chosen by 10-fold cross-validation over
  proteins (both epochs of a protein in one fold).
- **There is signal beyond noise.** Per-map gaps with |z| > 2 against protein bootstrap: alpha_R 99, pPII 97,
  alpha_L 76, beta 56 of 777 maps with >= 40 reads, against 35 by chance.
- **Most of it is the proteins, not the maps.** Split-half correlation of the per-map gaps between random protein
  halves: alpha_R 0.11, beta 0.05, pPII 0.11, alpha_L 0.29; between epoch 0 and 1 of the same proteins 0.63-0.70. A
  map's gap is mostly a property of the proteins that read it (as 1.12 found for the glycine offsets).
- **Cross-validation.** Held-out squared gap left: per residue 0.988 at the best width (0.3 nat), 1.008 at 1, 1.049
  at 5 (overfitting above ~0.5). On the map-level gap (the part an offset can act on): per-map 0.970 at 0.3 and
  1.03-1.13 for sigma >= 1; one offset set per central residue (a pooled model, diagnostic only) 0.956; one global
  set 0.965. The pooled models transfer better at every width, so the per-neighbour detail adds variance at 456
  proteins.
- **What transfers is a global helix push.** At sigma 0.3 the per-map fit moves 1,012 of 5,040 depths by more than
  0.1 nat, and its mean over each central residue deepens alpha_R for all 21 (-0.01 to -0.13) and raises pPII. The
  global fit is alpha_R -0.18, alpha_L +0.18, pPII +0.07 to +0.10 nat, the same in epoch 0 as in epoch 1: the free
  ensemble under-populates alpha_R by 0.9 points per residue. Freed depths would stabilise helices through the map,
  the work hb and dhb do in the same training (Jumper's capacity warning, 1.28).
- **Glycine agrees with ff3.0.** The per-central glycine offsets raise alpha_L (+0.35 to +0.40) and deepen alpha_R,
  the direction ff3.0's own glycine pair took (1.27: aR -0.107, aL +0.517).
- Limits: first-order response; c9's free ensembles under drifting non-Ramachandran terms; the share of 1/2 for each
  map; the held-out gap carries the held-out proteins' own sampling noise, so the removable fraction is small for
  every model.

### 1.30 The ff3.1 per-map trainer, and its two runs (2026-10-09)

The user chose to try both route-1 variants (10-09). Trainer `training/ff31/` (ConDiv.py, rama_basin.py): round 4's
gdepth trainer with only its Ramachandran calls generalised; deployed as `ff31_depth_all` and `ff31_depth_data`
(`checks/ff31_deploy_20261009`).
- **The model.** One offset per (coil map, basin), 840 maps (21 centrals with cis-proline) x 6 basins, started at
  zero (NDRD), each under its own Gaussian prior; per epoch one damped (ETA 0.5) Newton step on the MAP objective,
  solved jointly over the maps by conjugate gradient, with shares 1/2 for an interior residue's two maps and 1 for a
  terminal residue's one. Mode `all` trains all 5,040 depths at width 0.3 (1.29's cross-validated width); mode
  `data` the 1,442 with at least 20 native sites in the 410 training proteins (alpha_R 604, beta 458, pPII 317,
  alpha_L 39 of which 34 are GLY|X, beta' 13, pPII' 11) at width 0.2.
- **Tests before deployment** (`scratchpad/ff31_trainer_test/`). Its read assignment and accumulator, summed the
  round-4 way over the 37 GLY|X maps, reproduce ff30_gdepth_dt009's round 1 from its own epoch-0 data exactly (reads
  5,421 / 611; gap alpha_R -0.0023, alpha_L +0.0042; held out -0.0186, +0.0193). Zero offsets reproduce every map to
  1e-6 with the sheet group and the weights unchanged; an offset of +0.5 scales its basin by exp(-0.5) before
  renormalising. The CG solve reaches relative residual 4e-6, and mode `data` leaves every untrainable depth at
  zero.
- **The rounds** (`rama_rounds.txt`; "gap left" is the held-out map-level gap remaining). The first came out as
  predicted from c9's epoch 0 (`all` step rms 0.030, 68 depths beyond 0.1, 0.989; `data` 0.026, 3, 0.994):

| run, round (step) | step rms | max | depths beyond 0.1 | gap left | GLY\|X alpha_R / alpha_L |
|---|---|---|---|---|---|
| `all` 1 (18) | 0.0305 | 0.266 | 66 | 0.990 | -0.012 / +0.066 |
| `all` 2 (37) | 0.0222 | 0.146 | 235 | 0.9937 | -0.026 / +0.097 |
| `all` 3 (56) | 0.0192 | 0.108 | 345 | 0.9992 | -0.031 / +0.107 |
| `all` 4 (75) | 0.0191 | 0.120 | 418 | 1.0001 | -0.029 / +0.118 |
| `data` 1 (18) | 0.0262 | 0.136 | 3 | 0.994 | -0.006 / +0.034 |
| `data` 2 (37) | 0.0185 | 0.092 | 26 | 0.9959 | -0.011 / +0.050 |
| `data` 3 (56) | 0.0150 | 0.059 | 39 | 0.9983 | -0.014 / +0.054 |
| `data` 4 (75) | 0.0146 | 0.077 | 58 | 0.9994 | -0.017 / +0.062 |

- Offsets rms at round 4: `all` alpha_R 0.085, pPII 0.071, beta 0.056, alpha_L 0.043, beta' 0.021, pPII' 0.022
  (`rama_round_04.dat` md5 `3ffc75bf...`); `data` alpha_R 0.049, alpha_L 0.091, beta 0.035, pPII 0.049, beta' 0.080,
  pPII' 0.055 (md5 `4b6a9aaf...`). Margins at the target: `all` +0.070 (under +0.10 from step 62 on, and at steps
  53-57), `data` +0.141 (start +0.192).
- **The depths fit the training proteins and carry almost nothing to held-out reads.** `all`'s training-read gap
  went 0.00059, 0.00050, 0.00043 over rounds 1-3 (0.00046 at round 4) while the held-out transfer per round fell
  0.97%, 0.63%, 0.08%, -0.01%, 1.6% in all; `data`'s fell 0.56%, 0.41%, 0.17%, 0.06%, 1.3% in all, with under a third
  of the depths moved.
- The runs had no `after_training.sbatch`: the convergence gate judges only the Adam groups, so they stopped at step
  76 for the user (review in 1.34).
- `check_step.py` prints `basin depth offsets ...` (the run's `rama_basin.describe`) for a per-map run; on a round-4
  step its output is identical to the old version's (backup `check_step.py.bak_pre_ff31_20261009`).

### 1.31 The second converged round-4 gate: ff30_bio_si at step 76 (2026-10-10; retired)

ff30_bio_si's gate (`ff30_bio_si/gate_step76.txt`, steps 58-76, every group p > 0.005, dhb lowest at 0.0125)
released `epoch_03_minibatch_18` (bs_03) as `ff_3.0_bio_si`, with the frozen BioEmu-fitted glycine library and an
H-bond margin of +0.058 (under the +0.10 hold since step 26; start +0.192). It was retired with the BioEmu line
(1.32) and its validation cancelled unfinished. Its panel is in 1.21; its old-table TM4 (1.24) is not resolved from
bio_start (DSSP alpha 0.739 against 0.769, p 0.89), bs_02 (0.684, p 0.71) or d9_03 (0.726, p 1.00). Its seed s1
carries the largest jump on fixed inputs: the protein potential goes from -1578 to +200,071 at t 1490, is +2,368 and
+530 in the next frames and back near -1,500 by t 1700; no NaN, and TM4 stays at 0.99. The run used the old
side-chain/lipid table, whose single-evaluation kicks are in 4.10; this event has not been probed.

### 1.32 The third converged round-4 gate: ff30_bio_dt009 at step 76; ff3.0 chosen (2026-10-10; retired)

ff30_bio_dt009's gate (`ff30_bio_dt009/gate_step76.txt`, steps 58-76, every group p > 0.005, dhb lowest at 0.13)
released `epoch_03_minibatch_18` (b9_03) as `ff_3.0_bio` to both trees (`release_20261010-164651`, file-identical to
`checks/r4_epochs/b9_03`). Its rama.dat is ff_3.0_bio_si's (md5 `69051b78...`): both carry the frozen BioEmu-fitted
glycine library untrained. The released H-bond margin is +0.013, under the +0.10 hold since step 14 (ff_3.0_bio_si
+0.058, ff_3.0_gdepth +0.082). It was the first gate whose own glpG patch carried the converged side-chain/lipid
table (4.15).
- **10-11 23:50: the user dropped the BioEmu-map line for the basin-depth path** (expandable to every residue,
  ff3.1). ff_3.0_bio's and ff_3.0_bio_si's validation was cancelled unfinished; their results stand as the record
  (plan.md, "The basin-depth path").
- **10-12 00:25: ff3.0 is ff_3.0_gdepth, the only branch kept** (user); ff30_gdepth_si was stopped at step 72,
  before its gate, and deleted with its derived sets.

### 1.33 The steric check on ff3.1's moving maps (2026-10-10)

Jumper's point, through the user (10-10): a Ramachandran map free to move can open conformations where the side
chain and the backbone clash. In Upside the map is the only term that keeps a residue's CB off its neighbours'
backbone (the rotamer table is side chain to side chain, 20 x 20, and the coverage and environment nodes carry no
excluded volume), so the watch measures it (`training/ff31/steric_check.py`, plan.md Phase 14 step 5; outputs
`checks/ff31_steric/<run>_<step>.txt`).
- **The reference** is the hard-sphere Ramachandran map at Ramachandran's outer limits, from ideal geometry with CB,
  amide H, and a following proline's CD; glycine without CB. Its CB is the standard ideal one: on 1AO6 it lands 0.08
  A (median) from the crystal CB, so the chirality is L. Clashing share of the plane: CB 0.48 side chain / 0.44
  backbone (21% allowed, the textbook ~22% for alanine), no CB 0.00 / 0.44.
- **NDRD's own mass sits partly in it:** per coil map median 0.030 on side-chain clashes (alanine 0.020, asparagine
  0.089, valine 0.005; most of it alpha_L, by C(i-1)..CB) and 0.227 on backbone ones (alpha_R's upper phi edge by
  C(i-1)..C(i), the bridge by N(i)..H(i+1)). Proteins occupy those by widening N-CA-C, so the backbone class is read
  by its change only, and the alarm is on the side chain.
- **Why the trainer cannot open a new region:** an offset scales its basin's probability as a whole and keeps NDRD's
  shape inside it, so a map's clash probability moves only by the weights of the basins that hold it.
- **The flags** (watch step 3): a coil map or training residue whose clash probability more than doubles against
  NDRD (an offset of about -0.7 nat, beyond 2 sigma of the 0.3 prior), the mean side-chain clash probability up by
  more than 10%, or a step's free non-glycine beta' + pPII' above the run's own steps 0-18 (NDRD). The 2x and 10% are
  set by policy; the third is measured on the run it guards. The check sees the CB, not side-chain atoms beyond it.
- **ff31_depth_all and ff31_depth_data raised no map or residue flag in any round.** `all`'s largest coil-map ratio
  grew x1.10, 1.16, 1.22, 1.24 over its four rounds (side chain, GLU right GLY; backbone at most x1.223, VAL left
  LEU); `data` moved its clash classes about half as far (at most x1.052 side chain, x1.104 backbone). At round 4
  the training residues' mean side-chain clash probability is x1.0043 (`all`; depths x0.9984) and x1.0057 (`data`;
  depths x0.9994), most of the rise the sheet energies'; backbone x1.0223 and x1.0146 (the depths' part x1.0166 and
  x1.0081). The prediction from c9's epoch 0 was x0.999 side chain, x1.007 backbone for one round.
- **The beta' + pPII' flag is minibatch composition, not the maps.** ff31_depth_data's step 39 crossed it (free
  0.0050, restrained 0.0045, against steps 0-18's maximum 0.0047) with the clash classes unmoved. A minibatch is the
  same 24 proteins in every epoch, and that minibatch read 0.0047, 0.0042, 0.0050, 0.0047 free over the four epochs,
  NDRD's included. ff3.1.1 and ff3.2.1 (steps 0-113 and 0-101, read 10-25), which share one minibatch order, are flat
  by minibatch across epochs: minibatches 1 and 14 read 0.0041-0.0051 in every epoch, epoch 0 on NDRD included,
  minibatches 0 and 2 about 0.0022, and the restrained simulations of the same minibatches read as high
  (0.0038-0.0058), so the excess is in those minibatches' proteins. Every flag of either run (ff3.1.1 against its
  0.0043 at steps 33, 39, 71, 77, 90, 96 and 109; ff3.2.1 against its 0.0050 at step 90, 0.0051 free, 0.0041
  restrained) is minibatch 1 or 14; both means are 0.0032 before and after step 18. A maximum over a growing number
  of draws set against one draw per minibatch will cross by chance; a per-minibatch comparison against epoch 0 would
  separate a real shift from this. The threshold is not changed (the user's call). Until 10-24 the collector's
  status line counted the files carrying the flag (each file judges every step so far), not flagged steps.

### 1.34 ff3.1: the per-map review, the pooled prior and the merge (2026-10-13/27)

**The per-map review at the target** (user 10-13; panel.py's own functions in `checks/ff31_review/ff31_review_pairs.py`,
output `pairs_20261013.txt`; 44 domains: folded fraction over every domain and each class's |error| over the domains
both keep (36-40), each with a 95% paired bootstrap over domains, which holds the force fields fixed and so measures
panel sampling, not how two trainings of one design differ).
- **Two trainings of one design can differ more than the designs do.** After epoch 0 the ff2.1 control's two runs,
  fp_e00 (dt 0.015) and c9_00 (dt 0.009), differ by +0.141 folded [+0.091, +0.192]; gdepth's two, d00 and d9_00, by
  -0.006 [-0.040, +0.030]. Each pair also differs in dt. d9, fa and fd all fold 0.04-0.09 more than c9 at steps 18
  and 37, less than the control's own twin gap, so c9 cannot credit folding to the depths.
- **At the target `all` is worse and `data` ties ff3.0.** Against d9_03: fa_03 folded -0.058 [-0.101, -0.015], helix
  +0.001 (not resolved), beta -0.008 [-0.016, -0.002]; fd_03 folded -0.012 [-0.058, +0.035], helix -0.006 [-0.012,
  -0.000], beta -0.007 [-0.013, -0.002], gly_helix -0.039 [-0.088, +0.007]; fd_03 against fa_03 folded +0.046
  [+0.015, +0.082], helix -0.006 [-0.011, -0.001]. Over 16 pairs and 5 measures, an interval that touches zero, as
  fd_03's helix does, is within what multiplicity produces. fa loses folding from epoch 0 to the target (-0.065
  [-0.117, -0.021]) while its depths grow and its margin falls to +0.070; fd (-0.027 [-0.058, +0.003]) and d9
  (-0.013 [-0.046, +0.021]) hold, d9 losing beta (+0.009 [+0.003, +0.016]).
- **The depths carry no transferable rule.** Four rounds closed 1.6% (`all`) and 1.3% (`data`) of the held-out gap,
  falling to zero by round 4 while the training gap kept closing (1.30), as the probe predicted (1.29). `data`
  reached nearly `all`'s transfer with 58 depths beyond 0.1 against 418, so `all`'s extra depths fit the training
  proteins. `all`'s falling margin fits the capacity warning of 1.28 (depths taking helix work from the H-bond
  term); that is not shown.
- **Read.** On the panel neither earns replacing ff_3.0_gdepth. The user ordered the Peng benchmark with glpG TM4
  (10-13), and fd_03 was validated as `ff_3.1_depth_data` (plan.md Phase 14 step 6; glpG in 1.35). It and the
  per-map ff3.2 sets are retired as `*.bak_pre_merge_20261026` since the merge below.

**What holds the maps where they start** (10-21, user: "ff3.2 and ff3.1 barely moved their rama map ... they are
suppose to move rama map"; `checks/ff21_released_20261021/depth_data_vs_prior.py`, output
`depth_data_vs_prior_20261021.txt`, on the two angular-fix reruns). The data part of each step's depth gradient
(`T0 x train_gap`; the recorded `grad` adds the prior's share), apart from the prior: every epoch end sits at its MAP
point, offset = sigma^2/T0^2 x the epoch's data pull, so the prior's pull back to the start balances the data.
Single depths are noise (|t| > 4: 0-4 of 1,442, 1.2 by chance); the data's consistency is in the pooled subsets.
Glycine's GLY|X alpha_L: from ff2.1's map (ff3.1) the data pull it shallower in every epoch (t +5.1 at epoch 0, +5.6
at epoch 3), and the 0.2-nat prior centred on ff2.1 stops it at +0.07; from ff3.0's map (ff3.2) the data have no
consistent pull (t +1.1, +0.6; alpha_R -0.04, +1.67). So the data point to ff3.0's glycine map, ff3.2 stayed because
it started there, and ff3.1 stayed because its per-map priors hold 37 maps apart against a pull that is clear only
pooled. 1.29's probe found pooled offsets (one set per central residue) transferring better than per-map ones at
every width.

**The pooled prior, mode `pooled` of `training/ff31/rama_basin.py`** (10-21, plan.md Phase 29; checks in
`checks/ff311_20261021/`). Each map's offset is its central residue's shared offset plus its own deviation (`data`'s
1,442 at 0.2); 75 of 126 shared offsets reach 20 pooled native sites (alpha_R 21, alpha_L 12, beta 19, pPII 21, beta'
1, pPII' 1); pre-proline maps keep only their deviations.
- **Checks** (`check_pooled.py` on `ff31_depth_data_angfix`'s epoch 0): mode `data` of the new file gives the old
  file's step to 1.4e-17 and the recorded round 1 to 7e-17; `pooled` with no shared offset trainable gives `data`'s
  step exactly; spread and collect are adjoint to 4e-16.
- **Width** (`cv_pool.py`, 10 folds over the 456 proteins, the trainer's own round-1 step): held-out map-level gap
  left, ff2.1 start 0.9942 (width 0, = `data`), 0.9802, 0.9758, 0.9747, 0.9740, 0.9738, 0.9737 at 0.1, 0.2, 0.3, 0.5,
  1, 2; ff3.0 start 0.9938, then 0.9822, 0.9791, 0.9786, 0.9783, 0.9783, 0.9784. Width 0.5: the ff3.0 start's minimum
  and within 3e-4 of the ff2.1 start's, the curve flat beyond it. One round then moves glycine's shared alpha_L
  +0.185 from ff2.1 and +0.035 from ff3.0's map, the second consistent with the data having little left to pull
  there (above).
- **The two runs** are `ff311_depths_only` (ff3.1.1, from ff2.1's map) and `ff321_depths_only` (ff3.2.1, from
  ff3.0's), every other group held at ff2.1 (rates 0; every other group bitwise ff2.1's at step 76,
  `compare_params_step76.txt`), so the library is the only test.
- **Gauge.** A basin fit fixes a map's six offsets only up to one shared constant ('other' is nearly empty on the
  grid, condition 1.5e7): compare trained runs by their recorded offsets, and fitted libraries only by contrasts or
  maps.
- **At step 76 the pooled prior moves glycine's map; the per-map prior did not** (10-23; libraries
  `rama_round_04.dat`, md5 `c0244f91`, `5286b7e0`; `scratchpad/gly_rama_pooled_maps.py`, figures
  `~/Downloads/gly_rama_ff311_nonhelix_zero.png`, `..._ff321_...`). ff3.1.1's shared glycine alpha_L went +0.157,
  +0.217, +0.293, +0.358 over rounds 1-4, still rising ~0.065 an epoch: 1.29's +0.35 to +0.40. Its pooled glycine map
  on the non-helical zero: alpha_L -3.022 (ff2.1) to -2.611, alpha_R -1.562 to -1.551, PII -0.616 to -0.581,
  P(phi > 0) 0.647 to 0.606; trained alpha_L - alpha_R per GLY|X map +0.400 (+0.171 to +0.582). ff3.1, the same data
  under per-map priors, moved alpha_L to -2.952 only. ff3.2.1 stays at ff3.0's map: alpha_L -2.505 to -2.490,
  alpha_R -1.669 to -1.624, largest pooled-map change 0.16, shared alpha_L +0.039 to +0.021 over the rounds.
- **The two starts meet, except glycine's alpha_L, which is still closing** (10-23, user: "check other rama map
  terms and other parameters"; `checks/ff311_20261021/compare_pooled_runs.py` and `.txt`, the depth state of rounds
  0-4 from the checkpoints in `rama_state_pooled.npz`). In one frame (NDRD plus offsets; ff3.0 is NDRD plus one pair
  on the 37 GLY|X maps, alpha_R -0.1068, alpha_L +0.5169; rebuilt step-76 maps equal the libraries to 2e-6), over
  the depths the data reach: the 801 non-glycine maps start equal and end apart by rms 0.067, 0.078, 0.087, 0.092
  over rounds 1-4 (by basin mean within +/-0.02), the two runs' own noise; the glycine maps start apart by rms 0.209
  and end at 0.093, that floor. Per basin on glycine: alpha_R -0.101 to -0.015; the four minor basins open from 0 to
  +0.03-0.09, noise size; alpha_L +0.490 to +0.152. The shared glycine alpha_L gap goes 0.517, 0.398, 0.339, 0.256,
  0.179, closing 0.06-0.08 a round, mostly by ff3.1.1 rising (ff3.2.1 +0.517 to +0.537): within the other residues'
  shared-offset spread (rms 0.092, max 0.279 for SER alpha_L) by size, but systematic by trend, so the meeting point
  is not reached at step 76. Basin populations per map, largest |ln ratio| between the runs: glycine maps median
  0.356 at the start, 0.104 at step 76; non-glycine maps 0.094 (95th percentile 0.263).
- **Movement from ff2.1** (user: "check how much other AA's rama map moved from ff2.1"): non-glycine maps,
  read-weighted mean ln(P/P_ff2.1) alpha_R +0.068 / +0.061 (ff3.1.1 / ff3.2.1), alpha_L +0.074 / +0.059, beta
  -0.033 / -0.033, pPII -0.052 / -0.047; per-map largest basin change median 0.128 / 0.097, 95th percentile 0.355 /
  0.365; alpha_R up most on ALA, GLN, PRO, TYR, VAL (+0.10 to +0.15), alpha_L (sparse) most on GLU and LYS (+0.24 to
  +0.45), the runs apart on SER and LEU alpha_L; glycine maps median 0.30 / 0.36. Their rms movement grew 0.08, 0.11,
  0.13, 0.14 over the rounds, slowing, against the runs' noise of 0.09. Both runs continued to agreement (plan.md
  Phase 30).
- **At step 114 (round 6) the glycine gap is down to the noise except for a constant the data cannot see** (10-26,
  user: "is ff3.1.1 and ff3.2.1 currently equivalent?"; the same script over rounds 0-6, `dump_rama_state.py` taking
  every epoch end both runs reached; rounds 0-4 of the new dump equal the old one key for key, rebuilt round-6 maps
  equal both libraries to 2e-6). Phase 30's rule: (a) glycine's alpha_L shared gap 0.517, 0.398, 0.339, 0.256, 0.179,
  0.136, 0.072 against the other residues' shared-gap rms 0.000, 0.067, 0.077, 0.086, 0.092, 0.116, 0.084, met first
  at round 6; (b) the glycine maps' gap rms 0.209, 0.160, 0.139, 0.110, 0.093, 0.088, 0.092 against the other maps'
  0.000, 0.067, 0.078, 0.087, 0.092, 0.116, 0.084, met at round 5 only. Every round the glycine maps' gap carries a
  shared shift of +0.063 to +0.064 across their basins, unchanged from round 0: it is the shared component of ff3.0's
  pair as the basin fit placed it (alpha_R -0.107, alpha_L +0.517, the four minor basins set to 0), a direction the
  fit cannot fix and the data cannot move. Without it the glycine maps' gap is 0.198, 0.146, 0.122, 0.088, 0.065,
  0.059, 0.064 against the other maps' 0.000, 0.046, 0.042, 0.051, 0.070, 0.084, 0.063, and in basin populations the
  glycine maps' largest |ln ratio| per map is median 0.094, 95th percentile 0.162 at round 6 against the other maps'
  0.068 and 0.237. So the runs agree to their noise in everything the data constrain; rule (b) as written measures
  that constant too and cannot be met by more training.
- **The merge** (user, 10-26): ff3.1 is the mean of the two runs' round-6 offsets, `parameters/ff_3.1` (rama.dat
  `5b66c78a`; `scratchpad/ff31_merge/merge_ff31.txt`): the 84 untrained maps NDRD's byte for byte, the 840 trained
  ones NDRD plus the merged offsets to 9.5e-7, 38 of 5,880 basin populations outside the two runs' range by at most
  1e-4 (the mean is taken in energy, not probability); the other five files are ff2.1's, which both step-114
  checkpoints extract to within 2.2e-16. Pooled glycine map on the non-helical zero: alpha_L -2.574, alpha_R -1.510,
  PII -0.579 (ff2.1 -3.022, -1.562, -0.616; ff3.0 -2.505, -1.669, -0.616).
- **Why ff3.1's glycine is not ff3.0's** (10-27, user: "I suspect there is a defect in our training"; read only from
  both runs' step records, `scratchpad/ff31_gly_pull_20261027/`). Glycine's -ln(P_alpha_L / P_alpha_R) over the 37
  GLY|X maps (`gly_LR.txt`, gauge free): ff3.0 -0.598; ff3.2.1 -0.600 to -0.628 through round 4, then -0.702,
  -0.778; ff3.1.1 -1.219 closing to -0.821 at round 4, then -0.865, -0.872; the merge -0.825, their round-6 mean. So
  both runs turned away from ff3.0 in rounds 5-6, after step 76's restart, which changed nothing (engine 10-18,
  `ConDiv.py`, `rama_basin.py` the same files; the step-54 restart had no effect either).
  * The data pull on glycine's shared offsets per epoch (`pooled_pull.txt`, T0 x free - native collected, t over the
    19 steps): at ff3.0's map (ff3.2.1, epochs 0-3) none on alpha_L (t +1.1, +0.1, -0.2, -0.2) or alpha_R, so the
    trainer keeps ff3.0's pair; glycine's beta' pulled deeper in every epoch of both runs (t -5.1 to -1.7, -4.4 to
    -1.2), not converged at step 114 (-0.31, -0.20), the only shared offset of 75 with a persistent pull (the others'
    |t| median 0.5-1.1). ff3.0 never trained beta'.
  * The beta' deficit is broad: glycines native in beta' (phi 100-180, 494 of 2,972) that free simulations lose to
    beta, pPII and alpha_L; the 20 largest carry 16-17% (`gly_residues.txt`).
  * Rounds 5-6: glycine's alpha_R gap moved from about -0.003 to +0.004 per residue (t ~1 per epoch) while the other
    residues' helix balance stayed matched (|gap| < 0.0015; `helix_balance.txt`). The two runs share the minibatch
    order, so once their maps met the same noise moved both: their agreement at step 114 overstates precision. The
    pull fell below |t| 2 both near -0.60 (ff3.2.1) and near -0.85 (ff3.1.1): the data do not separate ff3.0's
    alpha_L - alpha_R from ff3.1's.
  * No code defect in the trainer; the test with a known answer is the debug run of plan.md Phase 32.
  * **The defect is the training target** (10-28; `native_fidelity.txt`, `beta_prime_response.txt`). ConDiv's native
    side is a replica under the force field itself with a weak distance restraint (spring 1/3^2, ~1 A RMSD, the
    original trainer's), which holds the fold but not the backbone dihedrals. At ff2.1 (epoch 0) it keeps each residue
    in its native basin: non-glycine alpha_R 0.979, glycine alpha_L 0.971, glycine alpha_R 0.897, glycine pPII' 0.613,
    glycine beta' 0.489, non-glycine pPII' 0.401, non-glycine beta' 0.004 (free ensemble 0.931 ... 0.001). For the 494
    glycines native at phi 100-180, as ff3.1.1's beta' offset went 0 to -0.271 the restrained population went 0.342 to
    0.371 and the free 0.272 to 0.304: the gap stayed -0.070 to -0.067, so the depth cannot close it and drifts at the
    damped step until the prior stops it. The 55 such glycines in the held-out proteins show the same gap (-0.06 to
    -0.12), so it is neither noise nor too few data. The fix (plan.md Phase 33) holds phi/psi in the native replica.
  * **The fix's step 0** (10-28; `ff32_nativefix`, `ff31_nativefix`, the same 24 proteins as ff3.1.1's step 0; the
    ceiling from `gate_ceiling.py` in `checks/ff31_fix_20261028/` and the repo's `scratchpad/`). The native replica now
    holds each basin close to what a 20 deg spring can: glycine beta' 0.656 (old target 0.356) against an ideal 20 deg
    Gaussian about each native's 0.726; over ff3.2's steps 0-10 (n 211) 0.754 against 0.769; glycine pPII' 0.829 against
    0.816, non-glycine beta 0.885 against 0.818, pPII' 0.879 against 0.876, alpha_R 0.997 against 0.987 (the force
    field's own well adds to the spring). The free-native gap on glycine beta' is now -0.343 (was -0.070), a pull the
    depth can act on. The one class the springs do not hold is non-glycine beta' (n 18): 0.071 against 0.582, the map
    at phi 100-180 outweighing the spring. The gate's 0.75 for glycine beta' was set before this ceiling was computed
    and is above it at step 0 (half these glycines sit within a spring width of a basin edge), so its FAIL measures
    the basin's edge, not the springs.

### 1.35 ff_3.1_depth_data keeps glpG's TM4 helical; the cold-rung glycine flip is absent from the corrected chains (2026-10-14/18)

ff3.1 in this section is the per-map `ff_3.1_depth_data` (retired, 1.34), not the merged ff3.1. Asked by the user
10-14: `checks/tm4_ff31_20261014/tm4_ff31.py` (outputs `tm4_ff31_<date>.txt`), remote_jobs.md §5's criterion (phi in
[-130,-20], psi in [-90,15]) on rotated groups only, T 0.70 (run.0) and 0.80 (run.14), against ff_3.0_gdepth on the
same group indices, its seeds differing only in the force-field datasets. Glycine phi > 0 on every rung:
`~/watch_probe/gly_flip_rungs.py` (midway2).
- **The corrected chains' block 1** (4.17; 10-17 to 10-18, groups 1-17, ~5,100 frames per rung;
  `tm4_ff31_20261018.txt`, `gly_flip_rungs_20261018.txt`, `tm4_dssp_block1_20261018.txt`, the DSSP from
  `extract_bb_dssp_block1.py` and `tm4_dssp_block1.py`). Pooled T 0.70 TM4 134-151: ff3.1 0.994 / 0.961 / 0.991 /
  0.977 (79ALA, 79HIS, 79ALA_S115T, 79HIS_S115T), ff3.0 0.990 / 0.997 / 0.993 / 0.982. The weakest residue is 134,
  TM4's first, in all eight (ff3.1 79HIS 0.52 with 135 at 0.81; ff3.1 79HIS_S115T 0.67; ff3.0 79HIS_S115T 0.70). No
  TM4 glycine is at phi > 0 at T 0.70 in any chain; over groups 14-17 the flips sit on the middle and hot rungs under
  both force fields (GLY136 or GLY149, 0.05-0.37 on rungs 15-27), and GLY143 never flips. DSSP alpha 135-151 on
  79HIS's T 0.80 rung: ff3.1 0.909 (first / second half 0.929 / 0.889; 135 0.72, 150-151 0.65 / 0.62), ff3.0 0.941
  (0.944 / 0.938; 135 0.73, 151 0.79); the core 137-149 holds under both.
- **The defective chains** (misplaced side-chain beads, 4.17; replaced by the corrected ones) held TM4 under all
  three force fields (79HIS DSSP alpha on the T 0.80 rung over block 2: ff3.1 0.863, ff3.0 0.945, ff2.1 0.887, ff3.1
  fraying most at both ends; figure `~/Downloads/glpg_tm4_dssp_remd_ff21_ff30_ff31.png` from
  `scratchpad/glpg_ff21_vs_ff30/plot_tm4_dssp_fixed_table.py` on `~/watch_probe/extract_bb_dssp.py`'s backbones).
  There 79HIS's GLY136 under ff3.1 sat at phi > 0 in 43-91% of T 0.70 frames from group 14 (pooled rung 0 0.66), and
  by block 2 79ALA's GLY149 flipped at the cold rungs too (0.05-0.13 on rungs 0-18); ff3.0's 79HIS had no glycine at
  phi > 0 at T 0.70. The walkers do not cross the ladder (5.5), so a flip near one rung is one walker's and these
  fractions are not equilibrium populations. With one chain each, whether the placement fix removed the cold-rung
  flip or the defective chain's walker history made it is not established.
- **Warm-rung glycine flips belong to all three force fields,** ff3.0's penalty included: ff2.1 on the same protocol
  (the cancelled `glpG-RKRK-79HIS.ff_2.1_sc` chain, 29 groups, 8,653 frames per rung; `~/watch_probe/gly_flip_ff21sc.py`)
  never puts GLY136 at phi > 0 on rungs 0-1, and flips GLY149 (up to 0.20) and GLY136 (up to 0.07) at the warm
  rungs 17-26.
- **The table, not the force field, removed ff2.1's unfolding.** The seed-1 runs behind
  `glpg_tm4_dssp_ff21_vs_ff30.png` (10-09) are runs_cov's old-table ff21_released and d9_03 (initial energy
  -25017.30 and -24992.61, as the MacBook Pro's logs). Redrawn from those files the figure reproduces; their
  fixed-table twins, which differ only in the table, read DSSP alpha 135-151 last quarter ff2.1 0.644 -> 0.994 (old
  table: 142-147 at 0.00-0.24) and ff3.0 0.989 -> 0.984 (`~/Downloads/glpg_tm4_dssp_old_vs_fixed_table.png`). Over
  12 seeds that is 4.15's ff2.1 0.610 -> 0.944, unwound 9 -> 0 (p 0.002). No local ff3.1 set exists; its TM4
  evidence is the REMD above.
- **ff3.1's glycine map alone does not make 79HIS's GLY136 flip** (10-15, read only).
  - ff3.1 trains from ff2.1 and does not carry ff3.0's glycine alpha_L penalty: the mean of its 37 GLY|X maps puts
    alpha_L - alpha_R at +0.079 E_up against ff3.0's +0.624 (per map -0.12 to +0.32; figures
    `~/Downloads/gly_rama_ff31.png`, `gly_rama_nonhelix_zero_page.pdf`, MacBook Pro).
  - The potential each glycine gets, from the 79HIS seeds' `rama_pot` (`~/watch_probe/gly_rama_seeds.py`;
    F(alpha_L) - F(alpha_R) at T 1 with the trainer's basin weights): GLY136 -1.32 (ff2.1), -0.70 (ff3.0), -1.36
    (ff3.1); GLY143 -1.15, -0.53, -1.10; GLY149 -1.37, -0.76, -1.31. At 136 ff3.1 is ff2.1's map to 0.04 E_up.
  - ff3.1's other parameters are ff3.0's, not ff2.1's: the same 76-step retraining changed every other file
    (hbond.h5, sidechain.h5, environment.h5, bb_env.dat, sheet) about equally, e.g. the H-bond energies -1.961 /
    -1.946 / -1.769 (ff2.1), -1.815 / -1.829 / -1.734 (ff3.0), -1.822 / -1.840 / -1.681 (ff3.1). So ff2.1 has neither
    the retrained terms nor the penalty, ff3.0 both, ff3.1 the retrained terms without the penalty, and only ff3.1
    flipped GLY136 at T 0.70 (defective chains): consistent with the penalty offsetting what the retraining loosened
    at TM4's N-terminal end, in one chain each.
  - ff3.1's whole coil library is ff2.1's to within 0.23 E_up in any populated cell (659 of 840 maps moved, median
    largest change 0.033); basin populations per map within 0.92-1.23 (alpha_R), 0.84-1.09 (alpha_L), 0.91-1.16
    (beta), 0.79-1.10 (pPII); the sheet group unchanged. The per-map priors centred on ff2.1 hold glycine there
    (1.34), while the same data pooled per central residue push alpha_L up +0.35 to +0.40 (1.29); ff3.0's pair was
    set at its start from BioEmu's weights (dL - dR +0.554, 1.19) and only kept by training.

### 1.36 Which finished force fields are at a fixed point (2026-10-20)

User 10-20: "check if the finished ffs needs no more trainings". Tools in `checks/converge_20261020/`:
`gate_window.py` (the gate's exact test over any one-epoch window; also `rot` split into its off-diagonal angular
logits (5,700), diagonal ones (300), radial (7,200) and the rest (18,220), and the depth gradient from each step's
`rama_step.npz`, data part plus the prior's share), `gram_probe.py`, `ang_pull.py` (each coefficient's t over an
epoch), `depth_pull.py`. `gate_window.py` reproduces `ff21_ctrl_dt009/gate_step76.txt` exactly and the step-38
audit's dhb 0.003.
- **Reading the gate's p.** At a fixed point the large groups give p 1.0000: over an epoch the step gradients cancel
  (mean pairwise cosine about -0.03, epoch sum 0.47-0.75 of the orthogonal-noise size), the training set's summed
  pull being near zero while each minibatch's is not. A pull shows as a sum above the noise size: no-BP's dhb 3.4x
  at step 38, ff3.0's glycine depths 2.3x in epoch 2.
- **The whole-group gate can pass a pulled subset.** No-BP's `rot` passes at step 76 (p 0.44), but its off-diagonal
  angular logits alone fail (p < 4e-6) in all three epochs from step 19: 118 and 106 of 5,700 with |t| > 4 in the
  two fixed-engine epochs, where noise gives 4.8; epoch sum 1.59, 1.45, 1.32 of noise; epoch-mean gradient 1.77,
  2.16, 1.90, not falling. Most sit near spline value 0.97 and move inward (96 of 118, 89 of 106), 0.19 logit an
  epoch, near Adam's ceiling for the group (0.0125 a step), spline value up to 0.06 an epoch. The control has 3 and 2
  such coefficients. So no-BP is not at a fixed point in ~110 angular coefficients, and reaching it would take
  several epochs; its mdCATH folding did not change over steps 38-75 (`mdc_nobp_03` against `mdc_nobp` -0.022
  [-0.057, +0.009], p 0.67; `paired_npz.py`, which reproduces `panel_compare`'s folded and Q from the npz).
- **The angular-gradient defect (3.12) left no detectable pull at a defect-engine fixed point.** The control sat at
  one at step 37; its first fixed-engine epoch (38-56) passes every group, off-diagonal angular p 1.000, 3
  coefficients with |t| > 4 (noise 4.8); `mdc_ctrl_03` against `mdc_ctrl38` -0.012 [-0.034, +0.010], p 0.67, Q
  -0.001. The ff3.2 rerun starts from ff3.0's step-75 parameters on the fixed engine, so its epoch 0 is the corrected
  gradient at ff3.0's fixed point: off-diagonal angular p 1.000, none with |t| > 4.
- **ff_3.1_depth_data and `ff31_depth_all` pass the gate** over steps 57-75 (every group p > 0.005, lowest dhb 0.033
  and 0.032; also over 19-37 and 38-56), and their depths show no pull, whole (p 1.000, sum 0.72 of noise), by basin,
  or on the glycine maps (0-1 of 1,442 with |t| > 4, noise 1.2).
- **ff_3.0_gdepth's glycine alpha_L depth is pulled the same way in every epoch.** The epoch-mean gradient on alpha_L
  has one sign in all four epochs (t +0.9, +2.2, +3.5, +2.1; alpha_R within +/-0.8), and the round steps did not
  shrink (alpha_L offset +0.007, +0.012, +0.021, +0.015; 1.26). Over the last epoch alone it is under the gate's cut
  (p 0.10), so by the gate's rule it is converged; the free-native alpha_L gap at round 4 is +0.009 on the training
  reads and +0.008 held out. Its Adam groups pass (1.26).
- **The ff3.2 rerun on the fixed engine is at a fixed point at step 76** (10-21; `ff32_depth_data_angfix`, files in
  `checks/ff21_released_20261021/`): `convergence_gate.py` passes every group over steps 58-76; `gate_window.py`
  passes every group and every `rot` subset and the depths over both steps 38-56 and 57-75; `depth_pull.py` finds no
  depth above |t| 4 (1.2 expected) and every subset's sign p >= 0.81, glycine's maps included. Its pooled glycine
  map equals ff3.0's within 0.02 E_up at alpha_L, alpha_R and PII on the non-helical zero
  (`scratchpad/plot_gly_rama_ff32_nonhelix_zero.py`; alpha_L - alpha_R per GLY|X map +0.626).

### 1.37 Every force field this workflow trains from ff2.1 folds worse than released ff2.1 (2026-10-21)

User 10-21: does the retraining's loss (1.21) reach the other trained force fields? 44-domain panel at time step
0.015, as `ff21_released` ran; all 44 domains; folded fraction (Q at least the all-atom 5th percentile), 4 seeds
pooled, paired against released ff2.1's 0.615 (`checks/ff21_released_20261021/panel44_trajectory.py`, output
`panel44_trajectory_20261021.txt`). ff3.1 and ff3.2 here are the per-map runs (retired, 1.34), not the merged ff3.1:

| run (tag) | start | epoch 0 | 1 | 2 | 3 (step 76) |
|---|---|---|---|---|---|
| control, probabilistic side chain (`c9`) | 0.615 | 0.390 | 0.405 | 0.375 | 0.430 |
| no-BP (`nb`) | 0.615 | | 0.389 | | |
| ff3.0, glycine depth (`d9`) | 0.603 | 0.477 | 0.457 | 0.439 | 0.464 (released) |
| ff3.1, data depths (`fd`) | | 0.479 | 0.472 | 0.449 | 0.452 (released) |
| ff3.1, all depths (`fa`) | | 0.471 | 0.444 | 0.404 | 0.406 |
| ff3.2 (`f32`, stopped in step 69) | | 0.486 | 0.465 | 0.432 | |
| ff3.1 rerun, fixed engine (`a31`) | | 0.436 | 0.449 | 0.462 | 0.420 |
| ff3.2 rerun, fixed engine (`a32`) | | 0.452 | 0.481 | 0.435 | 0.432 |

* **Every epoch of every run folds less than ff2.1** (p < 1e-6 each): ff_3.0_gdepth -0.152 [-0.196, -0.110], lower on
  37 of 44; ff_3.1_depth_data -0.164 [-0.202, -0.126], lower on 43. The starts are at ff2.1's level (`gdepth_start`
  0.603, p 0.5; `ff21_awh` 0.585), so the loss is the training, not the library change.
* **Most of it is the first epoch** (-0.13 to -0.23). Later epochs add -0.01 to -0.07 in ff3.1 and ff3.2, ff3.0 gives
  part back at epoch 3, and the control is flat (0.390, 0.405; on mdCATH its step 76 equals step 37, 1.36). So the
  cost is leaving ff2.1, not reaching the fixed point.
* **The round-3 mechanism** (1.17 A): with ff2.1's side-chain pair table put back, epoch 0 folds as ff2.1's start
  (`e00` 0.487, `e00_rot21` 0.583, ff2.1 -0.033); the table random-walks under Adam on a noise-dominated gradient. Not
  tested on the round-4 runs.
* **It holds at the training time step.** Released ff2.1 at 0.009 (`p44_ff21_released`, all 44 domains,
  `paired_npz.py`) folds 0.603, as at 0.015 (0.615, +0.013 [-0.025, +0.048], p 0.35); the control at its fixed point
  (`p44_ctrl_03`) 0.430, -0.173 [-0.227, -0.124], p 2e-7, lower on 40 of 44; no-BP (`p44_nobp_03`) 0.363, -0.240
  [-0.314, -0.171], p 3e-7. On TM-score (`tm_p44_*.txt`) released ff2.1 0.681 against the control's 0.650 (+0.031
  [+0.020, +0.044], p 4e-7, higher on 37 of 44) and no-BP's 0.623 (+0.057); TM > 0.5 0.886 against 0.832 and 0.824.
  The angular-fix reruns lose as their sources did (table; `panel44_trajectory_final.txt`).
* **Peng's benchmark agrees** (`checks/ff21_released_20261021/peng_trained.txt`): released ff3.0 0.480 native / 0.305
  de novo, ff_3.1_depth_data 0.464 / 0.326, the control 0.486 / 0.299 and no-BP 0.457 / 0.309 (1.21), against FF2's
  0.55 / 0.42; alpha3d native 0.15-0.35 and homeodomain 0.19-0.31 in all four trained force fields. The retired
  ff3.0 (FF1-form trainer, which trained no `hb` or `sheet`, 9e) kept them (0.61, 0.54) and scored 0.585 / 0.401 on
  the same Table S2 protocol (`ff3_benchmark/scoring/score_arms.json`). A lead on the H-bond, not a proof: that
  trainer also differed in burial (9t, 9u). Released ff2.1 itself now runs this benchmark on the same protocol
  (plan.md Phase 27; 31 of 32 arms scored by 10-27, the table `checks/ff21_released_20261021/peng_ff21.txt` when the
  last arm ends): native 0.552 over the 16 proteins against FF2's 0.55 (3.15), so the pipeline reproduces FF2 there.
* **What moves in the control's epoch 0** (`check_step.py`, steps 0-18): E_alpha -1.961 to -1.877 at ~0.0045 a step,
  one sign throughout, and the environment and backbone-burial terms likewise (env largest change 0.05 a step, bb
  scale -0.301 to -0.404): consistent pulls, not Adam's first-step transient. The side-chain table spreads as noise
  (rms 0.058 at step 18, lr x sqrt(19) = 0.054).
* **The training proteins do not fold better either**: free-replica CA-RMSD over all 456, per epoch, 2.85 2.86 2.86
  2.89 A (control) and 2.86 2.81 2.84 2.83 (ff3.1 rerun), epoch 3 against 0 p 0.37 in both. Not overfitting: the
  objective buys nothing on its own proteins while held-out folding falls.
* **Which files carry it** (`c9_00` with one file put back to ff2.1's, 44-domain panel at 0.015, all 44 domains,
  `checks/ff21_released_20261021/resets_c9_00.txt`): the side-chain table +0.089 [+0.049, +0.135], the backbone
  burial +0.081 [+0.044, +0.123], the H-bond +0.068 [+0.027, +0.111]; environment +0.028 and sheet +0.035 not
  resolved. Each reset stays 0.14-0.20 below ff2.1, and the three resolved ones add to about the whole loss (0.238
  against 0.225), so it is spread over the co-trained files, not one. Holding some files while the rest train does
  not remove it: the 10-07 runs with H-bond and sheet frozen folded less than their twins (1.21, `bz_00` 0.450
  against 0.478, `dz_00` 0.432 against 0.477), the trained terms compensating.
* **No port defect found.** Where `ConDiv.py` departs from the SI's wording (the NSE mixes the three coldest replicas
  0.6/0.3/0.1; the DSE takes the two replicas bracketing the Rg midpoint with Rg above 0.67 Rg(coldest) + 0.33
  Rg(hottest)), Kleinmann's port of Peng's code does the same (`condiv2.py` 690-735); its other drifts were corrected
  to the SI (9v). Peng's own file is lost, and ff2.1 was never this objective's fixed point (76 iterations without a
  convergence test, burial refitted after training, 1.23), so restarting the procedure moves it.

## 2. The hybrid model: what each side supplies

### 2.1 Replaced by design

The `rotamer` node's 1-body input is deliberately swapped from Upside's implicit-solvent coupling to the
explicit MARTINI SC-env table:

| | rotamer arguments |
|---|---|
| standard Upside | `placement_fixed_point_vector_only`, `placement_scalar`, **`hbond_coverage`**, **`hbond_coverage_hydrophobe`** |
| hybrid | `placement_fixed_point_vector_only`, `placement_fixed_scalar`, **`martini_sc_table_1body`** |

This is the original design; the live hybrid's rotamer also takes `hbond_coverage` and
`hbond_coverage_hydrophobe` (3.11, A83), and its 1-body uses `placement_scalar` on `rama_coord` (A72).

Upside's implicit bilayer (`membrane.h5` via `--membrane-potential` / `write_membrane_potential{,3,4}` /
`write_membrane_lateral_potential`, plus `--membrane-thickness`) is correctly absent: an explicit-lipid
model should not also carry an implicit slab. One part of it has no explicit counterpart: its backbone
H-bond term's cost for an unpaired NH or CO in the acyl core (4.7). Restoring that part alone is tested
in 4.8.

### 2.2 Absent and NOT replaced (findings 101)

MARTINI supplies only protein-environment interactions, so these protein-protein terms of stock Upside are
gone:

* `sigmoid_coupling_environment` <- `environment_coverage_sc`: the many-body protein self-burial term.
  MARTINI's 1-body rewards lipid contact, not helix-helix packing.
* `bb_sigmoid_coupling_environment` <- `environment_coverage_hb` + `cat_pos_bb_coverage`, and
  `hb_environment_coverage_hn/oc`: backbone burial coupling.
* `hbond_coverage` / `hbond_coverage_hydrophobe`: side-chain-to-backbone-H-bond competition, solved inside
  the rotamer solver in standard Upside (since restored in the hybrid: 3.11, A83).

These are core force field: all 24 master example configurations pass `environment.h5` + `bb_env.dat`,
including all six `08.MembraneSimulation` scripts, which use them alongside `membrane.h5`.

Restoring them is not a flag flip: it changes the `rotamer` node's arity and requires deciding how Upside's
protein-burial 1-body composes with MARTINI's lipid 1-body inside the rotamer solver. Scale check on real
structures: the self-burial term disfavours the drifted state by only **8.0 E_up = 5.6 kcal/mol**. The
experiment that restored them is findings 103 in section 4.1, and it did not recover the fold.

### 2.3 The backbone interface: sites, force routing, CB placement

**Only BB is on the protein side of the pair list.** Dry-MARTINI represents a residue backbone as one BB
bead standing for N/CA/C/O, and current builds and seeds carry one BB site per residue. (Earlier builds put
all five sites in the pair list at full epsilon; the energy over-count was only 1.15x, because O lands inside
env repulsive cores: +2117.5 against -2295.7 E_up at one frame.) The backbone `O` being the closest protein
atom to the environment (2.97-3.32 A, against 4.0-4.4 A for BB) is expected: O carries no MARTINI
interaction and nothing repels it.

**BB is a derived site built from nodes Upside already differentiates.** `HybridPositionNode` takes `pos`
and `infer_H_O`; BB is the mass-weighted centre of N/CA/C/O with weights [14,12,12,16]/54, where O is
Upside's own derived carbonyl O. Every term is a linear combination of node outputs, so `propagate_deriv`
is a constant-weight split: N/CA/C shares go to `pos.sens`, the O share into `infer_H_O`'s sensitivity,
whose chain rule carries it to CA, C and the next residue's N. No hand-written placement Jacobian is needed.

Newton's third law holds by construction: `martini_potential` runs on the hybrid node's output, a
pass-through copy of `pos` with only the BB and O slots overwritten, so an environment particle's
sensitivity returns to `pos.sens` unchanged while the protein side is redistributed by weights summing to 1.
Measured on 1rkl: force remaining on the BB and O slots exactly **0**, total force **1.69e-09** of the
largest single force.

The C-terminal residue has no acceptor: `infer_H_O` builds a carbonyl O from the next residue's N, so it
emits `n_res - 1` acceptors. That BB is the N/CA/C mass centre with weights renormalised, and its O slot is
left as it came in. This is not a guard: there is no such site, and no MARTINI term reads that slot.

**CB placement omitted the frame-origin subtraction (findings 102).** `affine_alignment` builds each residue
frame with its origin at the centroid of N/CA/C (`src/eig.cpp`: `center = (atom1+atom2+atom3)/3`), so a
point in that frame must be given relative to the centroid. `upside_config.write_environment` does that;
the hybrid prep stored the raw CB coordinate:

| | placement_data (frame coords) |
|---|---|
| standard Upside | `[-0.019807,  1.511741, 1.206801]`  = CB - centroid |
| hybrid | `[ 0.000000,  0.943756, 1.206801]`  = CB |
| difference | `[ 0.019807, -0.567985, 0.000000]`  = exactly the centroid |

The displacement is **0.568 A in Cartesian space for every residue**, moving the site that anchors the
whole sidechain-environment term (`martini_sc_table_1body` reads `placement_fixed_point_vector_only_CB`),
~12% of a 4.7 A bead contact radius. Fixed 2026-08-15: `CB_PLACEMENT` derives from an explicit reference
geometry with the frame origin subtracted, matching `upside_config.write_environment` to 5.8e-8 (`CB_VECTOR`
and `martini_build_tables.py` need no change: both use CB only as a difference). The total potential of an
identical configuration moves by **~+680 E_up**. Found by comparing the hybrid against a standard Upside
config of the same protein array by array: the bonded and H-bond core was bit-identical (`rama_map_pot`,
`hbond_energy`, `protein_hbond`, `backbone_pairs` max|delta| = 0), which made the two differing arrays worth
reading. Configs predating the fix: the four cluster REMD chains and their seeds, and the local production
seed behind the delivered HDX figure; the NP campaign was rebuilt (section 9).

### 2.4 The integrator the hybrid actually runs (findings 123)

In the pre-audit engine `DerivEngine::integration_cycle` (`src/deriv_engine.cpp:396`) began with

```
if(martini_brownian::has_brownian(this)) {
    compute(DerivMode);
    martini_brownian::apply_langevin_step(this, mom, dt);
    return;
}
```

so with `/input/brownian` present (every hybrid config) it returned before the three-stage Predescu et al.
(2012) integrator (`mom_update = {1.5-3a, 1.5-3a, 6a}`, `pos_update = {3b, 3-6b, 3b}`) and ran **one** g-JF
stage at dt = 0.009:

`x <- x + (b dt/m) p + (b dt^2/2m) f + (b dt/2m) beta`,  `b = 1/(1 + alpha dt / 2m)`

One stage commits a force spike to displacement with no intermediate force re-evaluation. Later changes to
this path: A5, A48, A149, A177 (3.14).

`/input/brownian` covers **4529 of 4949 atoms**: all 272 ions, all 3627 lipids and **630 of 1050 protein
atoms** (N/CA/C), with friction 0, 0.1692, 0.3384 and 0.5075, higher near lipid within a 12 A cutoff, so the
protein backbone is on the single-stage Langevin path and carries the interface friction. Measured
single-step displacements `b dt^2 |F| / 2m`: **protein max 0.010 A** (99.9th 0.006, median 0.0005) against
**lipid max 0.0003 A**, ~30x, from one sixth the mass and stiffer forces.

`--integrator mv` does not help: `build_integrator_levels` makes `integrator_level == 1` the slow set at
`dt * inner_step`, a cost optimisation for expensive smooth terms, and every hybrid node takes the default
level, so `mv`'s slow set is empty.

**Fix (findings 124): RESPA-style g-JF inner sub-stepping.** `n_inner_steps = N` loops position update,
`engine->compute(DerivMode)` and momentum update N times at `dt_i = dt/N`. The outer dt, the
`numerical_time_step` check and the friction/diffusion calibration stay at 0.009. In
`src/martini_brownian.cpp`: `BrownianRuntime::n_inner_steps`, read from the `/input/brownian` attribute
`inner_steps` (refused if < 1), with the invocation counter `outer_invoc x N + inner` keeping random streams
distinct.

| r [A] | F [E_up/A] | kick N_inner=1 [A] | kick N_inner=9 [A] |
|---|---|---|---|
| 2.853 | 1.27e5 | 5.1 | 0.064 |
| 2.584 | 4.64e5 | 18.8 | 0.232 |
| 2.432 | 1.02e6 | 41.3 | 0.510 |
| 1.783 | 5.60e7 | 2268 | 28.0 |

At every thermally reached distance (>= 2.43 A) N_inner=9 keeps the kick below 1 A (at 1.783 A the LJ is
~10^7 kT). Dissipation per outer step is `prod_{i=1..N} (1 - alpha dt_i/(2m)) ~= 1 - alpha dt/(2m)`, as at
N=1, so `D = kT/alpha` is unchanged. Local timing: 12.9 ms/step at N=1 against 55.4 at N=9 (4.3x, because
`engine->compute()` is ~41% of step time). Thermodynamics: `avg_kinetic_energy/1.5kT` 1.011 baseline and
1.002 fixed, potentials ~-22 000 E_up, Rg ~20.5 A in both arms. No blow-up was captured locally (the 79HIS
configs do not sample 2.43 A contacts in that window); a clean local A/B needs a toy two-body system or a
pre-failure cluster frame.

`inner_steps` defaults to 4 in `py/martini_prepare_system_lib.py` (3.10a); since A177 the binary requires the
attribute (and `@reference_temperature_up`), so a config written without it is refused.

### 2.5 The friction clock

The live calibration is the sub-molecular fallback: `D_raw = 4*11.5 = 46 um^2/s`, `dt_raw = 40/4 = 10 ps`,
`D_bead,up = D_raw*1e-4*dt_raw/.009 = 5.1111 A^2/U`, and `alpha_bead = kT/D_bead,up = 0.1691804`. Each
environment bead receives this friction; a real protein N/CA/C carrier receives `n_contact*alpha_bead`,
where `n_contact` is the number of lipid beads inside the existing 12 A spline cutoff. Counts are refreshed
at each stage handoff and minimization promotion, then held fixed for the stage and every production
continuation (so a continuation is the continuous run), so the SDE does not silently acquire
position-dependent multiplicative noise.

**Per-lipid targets** (user 10-20). The 11.5 was copied from Lipid14's Table 9 (Dickson et al. 2014), which
quotes pfg-NMR at 303 K (Filippov, Oradd and Lindblom, Langmuir 2003); the same group reads 9.8 at 30 wt%
water. No study measures DOPC, POPE and POPG by one method. The best same-method evidence puts POPE at
0.75 +- 0.15 of DOPC (a fluid PE-rich bilayer against POPC, MAS pfg-NMR, Gawrisch lab) and POPG at
1.0 +- 0.3 (faster than PC at low salt, slowed by salt, unmeasured at physiological salt).
`LIPID_LATERAL_DIFFUSION_UM2_S` holds DOPC 11.5, POPE 8.6 and POPG 11.5 with these citations and refuses any
other lipid; ions keep the DOPC clock as an interim rule.

Measured molecular DOPC diffusion under it is only **0.013-0.015 um^2/s against the 11.5 target**, reported
as a failed molecular target. Name the calibrated observable in H5 and in the paper, and never report a
friction-calibrated trajectory as having the target lateral diffusion.

`1 T_up = 350.588235 K`, so `--temperature 0.8647` is 303.15 K; the MARTINI factor four changes time, not
thermodynamic temperature. One temperature controls the whole system: the workflow assigns the DOPC
friction reference from the single authoritative `TEMPERATURE` and overwrites any independently supplied
value, because calibrating friction at one kT while driving its noise at another changes the nominal
diffusion.

**Never accept a kinetic calibration from temperature and structural stability alone.** An earlier mapping
(`tau_up = .0036`, `alpha*dt/(2m) = 1.25` for a mass-6 bead) gave effectively frozen trajectories (0.0081 A
of drift-removed COM motion per saved frame) with a thermal-looking momentum distribution and retained
secondary structure. Gate every friction change on protein displacement and whole-molecule lipid MSD in the
saved trajectory, and inspect the VTF rather than only the H5 statistics.

---

## 3. Known defects: root causes and fixes

### 3.1 The LJ core table was force-free, and particles reached it (findings 90, 92, 93)

`py/martini_build_tables.py` evaluated the grid at `r = max(r, 0.1*sig)` on a domain starting at r = 0, so
below 0.1*sig (0.47-0.60 A) the tabulated potential was a constant and exerted no force. Per-step
instrumentation (`UPSIDE_MARTINI_PAIR_DIAG`, two 48-replica jobs with exchange disabled, 340 000 steps):

| | job A | job B |
|---|---|---|
| reported approaches < 1 A | 4742 | 5032 |
| approaches < 0.6 A (inside the floored plateau) | 1446 | 1432 |
| closest approach | **0.0355 A** | **0.0466 A** |
| largest force delivered anywhere | 3.4e11 E_up/A | 1.0e11 E_up/A |

At 0.0355 A the true dry-MARTINI LJ force is 5.8e29 E_up/A, about 18 orders of magnitude above anything the
table delivered. Offending pairs are environment-environment (LIPID-LIPID, LIPID-ION), which is why
`lipid_kinetic` explodes while the protein KE is merely NaN. The event is a single-frame catastrophe from a
healthy state (potential -7070, KE 1.51/1.37; non-finite with `lipid_kinetic` 6.246e18 ~60 steps later), and
exchange then carries it into every slot (48/48 replicas destroyed from one event).

Findings 90 had argued the core was ~500 kT out of reach from the one-step displacement of an inertial
mass-72 bead. That was wrong: ION and LIPID are integrated as overdamped Brownian, whose step is
proportional to the force, so a large force can overshoot through a partner, and inside the plateau nothing
ejects it. **Check which integrator governs the particles before computing a stability margin for them.**
Ruled out by measurement: a stale pair list (`cache_buffer` 2.0 A, rebuild check before every evaluation),
minimum image (`roundf(dr/box)`, correct for unwrapped ions at 489 A in a 137.4 A box) and the table
formula (2.3e-12).

Fix, two coupled changes:
* `py/martini_build_tables.py`: floor removed, `PARTICLES_R_MIN_A` 0.0 -> 0.3, grid built over the true
  potential and asserted equal to the analytic form to 1e-12 relative (the core at 0.3 A is ~5e17 E_up/A).
* `src/martini_potential.cpp`: `r_min`/`r_max` were hardcoded to [0, 12], ignoring the table's
  `r_min_ang`/`r_max_ang`; the domain and the grid size are now declared once by the builder and read from
  the table.

On glpG-DDM the corrected table gives initial potential -7811.76, min pair 4.0349 A and max force
33.3441 E_up/A, identical to the old table, while the core delivers 8.5e15 E_up/A at 0.35 A where the old
table delivered ~0.

**The fix removed the consequence, not the entry (findings 93).** The same system still reached 0.2078 A
(6 approaches under 0.3 A, forces to 1.27e15 E_up/A), now receiving an enormous impulse rather than
coasting; the clamped spline below `r_min` is still a constant. The core is where the cascade ends: one
capture shows two environment beads at 1.83 A (max force 3.4e6, inside the valid domain), then 1.59 A on a
BB-proxy/lipid pair, 1.05 A at 5.5e10, and only 1072 steps later 0.245 A.

**The MD loop destroyed its own error messages (findings 93).** `src/main.cpp:1252` integrated systems under
`#pragma omp parallel for` with no exception trapping, so any engine `throw` mid-run called `std::terminate`
(`libc++abi: terminating due to uncaught exception`). It now traps per system, keeps the first message and
rethrows once serial. A diagnostic that cannot report is worse than none.

### 3.2 The blow-up mechanism: a 1 m_up protein site ejected by the MARTINI wall (findings 122)

Localised by re-running a 10-frame window around the onset of a `glpG-RKRK-79ALA` event on the local engine.
At frame 105 the total potential is +2.536e5 while Rg (19.6 A) and |pos|max (130 A) look normal. The excess
is all `Spring_bond` (280 -> 274 347 E_up), in three bonds of one residue: **CA of residue 170 sits 80.7 A
from its own C (r0 1.526) and 69.8 A from its own N (r0 1.453)**, worth 150 522 + 112 250 + 8 502 E_up. The
recorded coordinates reproduce the recorded potential to 0.03 E_up.

Every protein site carries **mass 1 m_up** and every lipid and ion bead **6**, so the same force throws a
backbone site six times as far. dt for a 1 A one-step kick `F dt^2 / m` at the pair forces of 2.4's table
(the g-JF factor b/2 halves it; the order of magnitude matters): 0.0028 at 2.853 A, 0.0015 at 2.584 A,
**0.00099 at 2.432 A**, 0.00013 at 1.783 A, against dt = 0.009. The run lives there: the closest pair is
2.49-2.74 A in every clean frame (1.78 A in the bad one), with ~169 pairs per frame inside 3.40 A, ~15 inside
2.85 A and **0.2 per frame inside 2.43 A**, so an ejection is the expected outcome whenever the bead at that
separation is a protein site. The steepest tabulated force, 5.2e17 E_up/A at the 0.3 A inner edge, gives the
8e18 potentials and |pos| 6.9e10 seen once a pair is driven there.

Consistent with this: the ejected atom is the lightest species; `martini_hybrid_position` rises with it
(1700 -> 12 988); the wreck moves ~4 rungs per frame, matching 3.7 exchanges per frame (37 steps/frame,
`--replica-interval 0.09` = 10 steps); an 80 A excursion cannot relax in 37 steps (bond period ~100 steps,
thermostat 555), so the clean frame after it is a different configuration swapped in.

**Global observables are the wrong instrument for a local failure.** Rg, |pos|max and the peptide C-N scan
all passed a frame carrying +2.5e5 E_up; term decomposition found it in one step.

### 3.3 REMD launders a wreck around the ladder, and the driver rolls it back

All four POPE/POPG variants carry non-finite potentials in about 0.15% of frames, in nearly every replica
file, in pairs two frames apart on the exchange period with clean frames around them. The first reading, an
output defect, was wrong. `run_remd.py` (`destroyed()`) treats any non-finite potential in a chunk as a
blow-up, rolls that replica back to its pre-chunk positions and rotates the NaN chunk to
`output_previous_N`. A replica genuinely blows up and exchange carries the wreck around the ladder, so the
neighbouring frames held other, healthy configurations. `grep ROLLBACK` confirms it, with events clustered by
chunk (one chunk rolled back all 48 replicas; one replica needed five consecutive rollbacks).

The analysis is sound because `martini_remd_concat.py` keeps a frame only if its potential is **finite and
negative**: a condensed bilayer plus protein sits near -2.2e4 E_up, and its docstring records a replica
finite for 96 frames at +1.9e6 before NaN, which a NaN-only filter would have kept.

Lessons: **a clean neighbourhood is not evidence of a clean trajectory** (read the producer before
explaining its output), and **check the log the tool already writes before inferring a mechanism from the
data**.

### 3.4 The four cluster POPE/POPG jobs were simulating a RIGID protein (findings 116)

The cluster HDX came out empty (188 of 203 amides off scale, values to -53.9 kcal/mol). MBAR was healthy
(ESS 6529 of 48 624, top single-frame weight 0.0002, ladder overlap better than the local run's);
equilibration was not the cause (the last quarter alone gave the same profile), nor the membrane term
(protein-only p_f == 1 for 169 of 203 amides, 171 with the lipid term). The protein had no internal
dynamics: in the raw trajectory the internal CA RMSD between chunks is **0.000-0.001 A**, Rg 17.60 +- 0.000
and H-bond count 194.6 +- 0.01 at both T = 0.70 and 0.90, while it moves as a rigid body (per-atom std
1.0 A). Local control: Rg 17.67 +- 0.228, RMSD 2.46 +- 0.563, H-bond 169.7 +- 12.07.

**Root cause: `/input/stage_parameters.current_stage`.**

| | local seed | cluster seeds |
|---|---|---|
| `current_stage` | `production` | **`production_handoff`** |
| `activation_stage` | `production` | `production` |
| `preprod_protein_mode` | `rigid_body` | `rigid_body` |

`martini_hybrid.cpp:637-641` `enforce_preprod_rigid_stage` returns `preprod_rigid && (stage != "production")`
and calls `martini_fix_rigid::set_dynamic_rigid_groups`, so `production_handoff` held the protein rigid. The
earlier verification read only `hybrid_interface_active_stage` (`martini_hybrid.cpp:646-647`), which accepts
`production_handoff`: a different predicate on the same string with the opposite accept set. **When a stage
string is checked anywhere, enumerate every site that reads it, and verify the conclusion dynamically (does
the protein's internal RMSD change?).** The bilayer hole and the -2.5 sigma potential drift were downstream:
lipids relaxing around a protein that cannot respond. Fix: `set_stage_label(seed, "production")`, as
`martini_prepare_system.py:1303` already does. Upside 3 replaced these stage rules (A68, A162).

### 3.5 MBAR silently returns uniform weights for a hybrid coupled potential (findings 91)

`helpers/calc_hdx_ht.py` and `4.calc_D_uptake.py` built `beta[l] * cE0[k]` from raw energies with no
reference subtraction. That is fine for protein-only potentials, O(1e2-1e3). A hybrid `Energy.npy` is the
full coupled-system potential (~-7.6e3 E_up for a protein plus a DDM micelle), so `beta*U` reaches -1.2e4,
`exp(-u)` overflows and the solver never leaves `f_k = 0`.

Measured on 48 states x 423 frames: raw gave f_k spread **0.000**, neighbour overlap **0.0000** and **0/423**
columns weighted at every target temperature but the bottom rung; mean-subtracted gave f_k spread **71.07**
rising monotonically, overlap 0.115-0.128, all 423 columns weighted, ESS 671-3378 of 20304.

`f_k = 0` makes every weight equal, so the estimator returns an unweighted average over the pooled ladder
while reporting it as a reweighted ensemble at one temperature, and raises nothing. Tell-tales: a gradient
norm of exactly `sqrt(n_rep-1) * n_frames`, `max_delta` exactly 0, ESS exactly `n_rep * n_frames`. Two
variants had already produced plausible-looking dG plots this way.

Fixed in both files by referencing `cE0` to its pooled mean, which is exact (f_k shifts by `-beta_k*C` and
`exp(beta_target*C)` cancels in the normalisation) and leaves the protein-only path numerically identical.
`03.TrajectoryAnalysis/2.mbar_meltingCurve_freeEnergy.py` and `04.HDX/4.calc_HDX.py` have the same
construction but are only fed protein-only energies, so they were left unchanged. Passing the pymbar-3
style 3D `u_kln` under the installed pymbar 4.0.3 is not a bug (3D and 2D give identical `f_k`).

**When a solver reports a gradient that is an exact function of the array shape rather than of the data, it
has not solved anything.** Check that before reading any number downstream of it.

### 3.6 Trajectory assembly: segments are not safe to concatenate

* **Production seeds carry an equilibration `/output`.** `run_remd.py` materialises replicas by copying the
  production seed, which already holds 300-400 single-temperature handoff frames; on the first reseed they
  rotate to `output_previous_0`. The concat joined oldest-first and `get_info_from_upside_traj.py` takes the
  temperature from the first frame, so all 48 replicas were labelled T = 0.8647, the degenerate
  uniform-weight condition of 3.5, and pymbar only warned ("States 45 and 47 have the same energies").
  Fixed: the production temperature is read from the newest chunk, and any chunk disagreeing is dropped and
  named.
* **Dropping destroyed chunks whole throws away most of a good trajectory.** One replica blew up 2333 frames
  into 2448, so a chunk-level rule would have cost 2332 good frames to remove ~100 bad ones. Now filtered per
  frame on finite AND negative total potential; whole-run records such as `replica_swap_partner` are
  detected by their first axis and left out with a note.
* The local run was clean while the cluster was broken purely because of how replicas were materialised
  (`warm_start.py` builds from a `seed.up` with no `/output`). **Testing one does not test the other.**

### 3.7 Analysis-code defects found in the shipped workflow

* **`k_chem` defaults to ~400 K.** `4.calc_D_uptake.py` defaulted `legacy_T_range` to `[1.14]`, so exchange
  rates were evaluated at T_up = 1.14 ~ 400 K whatever the ladder: `k_chem` reached 4.9e5 s^-1, every amide
  exchanged in milliseconds, and every normalised curve was the same step (COF 27877.86 for all 63
  peptides). At `legacy_T_range=0.85` (298 K) `k_chem` is ~10 s^-1 for a mid-chain serine and the 63 curves
  are distinct. Check `k_chem` in `<pdb>_percentD_feats.csv`, not the exit code; the only alarm was
  matplotlib's `TwoSlopeNorm ... must be in ascending order` 200 lines downstream. Upside 3: A64.
* **`_DG_Hbond.png` free-energy scale is 15% low.** `calc_hdx_ht.py:337` used `g = -0.593*np.log(hist)*t`
  with `t` in reduced temperature; 0.593 kcal/mol is kT at T_up = 0.85, so the correct factor is
  `kB*t*350.588 = 0.6966*t` (low by 0.851). Fixed in Upside 3 by A62. `make_hbond_landscape_figure.py`
  computes it correctly and drops bins carrying less than one effective frame (without that, the 245 K curve
  reads 82 kcal/mol where the reweighting has no support).
* **COF is not the readable form of the uptake comparison.** It integrates the squared derivative of a curve
  normalised by `(max - first)`, so a nearly flat curve inflates by orders of magnitude (experimental 21 to
  9.1e4, simulated 129 to 2.5e4) and a linear R^2 is dominated by that normalisation. Use the rank
  correlation (Spearman 0.30, p = 0.023, n = 57), ranking within each dataset.
* `5.analyze_D_uptake.py` sliced the 0-d `<pdb>_<sim>_<i>_T.npy` by frame (now `np.atleast_1d`), and
  `helpers/write_hybrid_energy.py` emitted a flat `Energy.npy` where the path indexes `[:, 0]` (now
  `reshape(-1, 1)`).
* **The Python engine computed a different model than the binary.** `engine_c_library.cpp` never called
  `load_masses_for_engine` / `register_fix_rigid_for_engine` / `register_stage_params_for_engine` /
  `register_hybrid_for_engine`, as `main.cpp` does, so `upside_engine` evaluated the hybrid with its
  interface inactive: **-12827 vs -18172 E_up** on the same coordinates. Every `get_output`/`energy` result
  taken through the Python engine before 2026-08-15 is suspect.
* When a node's arity changes, the migration is part of the change: a two-argument
  `martini_hybrid_position` cannot load a config declaring one, so every existing `.up` became unloadable
  when the binary was replaced. A config held open by a running job can be migrated only between blocks, so
  the migration belongs in that job's submit script.
* A silent fallback is worse than a missing input: a `0.0` belt half-thickness made a solvation gate accuse a
  correctly inserted protein, and a metadata-PDB lookup hardcoded to `example/16.MARTINI/pdb/<id>.MARTINI.pdb`
  (the workflow writes it to `<run_dir>/hybrid_prep/`) labelled every VTF's lipids `UNK`.
  `find_martini_metadata_pdb` now searches the run directories implied by the `.up` paths passed in.

### 3.8 Three VTF-generation bugs, all fixed at the root (findings 126; the third 2026-10-09)

**Bug 1: periodic images. Wrap per molecule, and only after the molecules are whole.** Two faults, found
five months apart, are one mistake seen from two sides.
* *Torn lipids.* `extract_trajectory` wrapped every particle into the box via `centralize_system` and never
  unwrapped, so molecules straddling a face were split: on the seed file (400 frames, 4187 bonds) max bond
  **141.00 A**, 3.254% of bond instances over 50 A. The worst was a PO4-GL1 bond inside one lipid, so a
  protein-only integrity check passes while the file is unusable.
* *The protein a box length out of the bilayer.* Running the bond walk after the per-particle wrap only moved
  the fault: the walk rebuilds each molecule around whichever anchor the wrap left inside, for glpG atom 0,
  the floppy N-terminal amide. Whenever the tail crossed a face the whole protein followed its image
  (`glpG_RKRK_79HIS_run0_remd.vtf`: 159 of 1822 frames, protein-lipid xy separation up to 98.18 A = one box
  length), rendering the protein beside a protein-shaped hole while the trajectory itself was fine.
* *Fix, 2026-09-14.* `build_molecule_topology` returns the bond walk and a connected-component label per
  particle; `extract_trajectory` calls `unwrap_molecules` first, then `centralize_system` shifts by the
  protein centroid and wraps each molecule by `box * round(centroid/box)`. After: protein COM exactly 0 in
  all 3146 frames, protein-lipid xy separation mean 0.70 A / max 2.20 A, no declared bond over 10 A except
  the residue-210 C-O (bug 3).

**When validating a VTF, check the declared bonds across all frames and that every molecule centroid is
inside the cell.** Either check alone passes one of these two bugs.

**Bug 2: mode 1 on a hybrid system emits the protein twice, and VMD rejects both copies.**
`build_mode1_mapping` emits the MARTINI-side protein (N/CA/C/O plus BB beads, all named `PRO`, no bonds) and
the appended all-atom backbone with real residue names, so `atomselect protein` returned nothing usable.
Mode 1 is the wrong mode for a hybrid: `build_mode2_mapping` keeps `protein_membership < 0` (every
environment particle) plus the sequence-named all-atom backbone. The library CLI already auto-detects it from
`input/hybrid_env_topology/protein_membership`; the analysis scripts called `build_mode1_mapping` directly
and were switched. Verified in VMD (`protein` -> names {C, CA, N, O} only, `name BB` -> 0). Coordinates can
differ by one box length between modes, since the unwrap anchor depends on atom ordering.

Selection notes: `infer_residue_names_from_class` left lipids `UNK` (it names a lipid only for `cls == "OTHER"`,
and then `DOPC`), so select them with `chain X` (A80 lists the `UNK` labels among its fixes). Do not relabel them
`DOPC`: this is POPE/POPG (head-group bead `NH3` = POPE, `GL0` = POPG). Residue numbering
is trustworthy: `input/sequence` position 79 is HIS or ALA and 115 is SER or THR as the variant names imply.

**Bug 3: the writer drew a phantom atom, the C-terminal O slot (fixed 2026-10-09).** First read as an
unconstrained O escaping (12b); in fact the O slot of glpG residue 210 never moves (0.0000 A over 401
frames) while its C moves away. The engine writes a residue's O slot only where infer_H_O has an acceptor on
that residue's carbonyl C (`src/martini_hybrid.cpp`, `resolve_bb_o_hbond_elements`); the C-terminal residue
has none. Nothing reads the slot (no MARTINI pair, potential node or `/input/brownian` entry), but
`py/martini_extract_vtf.py` wrote it as an atom. Fixed there by the engine's rule:
`build_backbone_projection_map` records `carrier_present` (N, CA, C always; O where an acceptor sits on the
residue's C), followed by the atom list (`backbone_output_atoms`, shared by modes 1 and 2), the frame
assembly and the bonds (`backbone_bonds`, replacing `mode2_backbone_bonds`). glpG now has 839 backbone atoms
(209 O), bonds N-CA 210, CA-C 210, C-O 209, C-N 209; O positions equal the trajectory's to 0.001 A. In the
user's commit 14a27cc3. The regenerated `glpG_79HIS_ff3.0_T080_s1.vtf` (ff3.0, seed 1, no TM1 term) is in
the MacBook Pro's `~/Downloads`.

### 3.9 Upside traps on exit under clang whenever Monte Carlo is enabled

`MonteCarloSampler` (`src/monte_carlo_sampler.h:12`) is abstract but had no virtual destructor, while
`MultipleMonteCarloSampler` deletes `PivotSampler`/`JumpSampler` through
`std::vector<std::unique_ptr<MonteCarloSampler>>`. That is guaranteed UB, which clang on arm64 compiles to a
trap: **SIGTRAP (exit 133)** in `~MonteCarloSampler`, after the run has finished and flushed. GCC on midway2
does not trap. Bisected to `--monte-carlo-interval` alone; master has the identical defect.

The consequence was silent in the wrong direction: ConDiv's worker raises `RUN_FAIL` on
`j.job.wait() != 0`, so all 12 workers of a local minibatch reported `WORKER_FAIL` on valid data (250/250
frames, Rg 13.2-13.8 A, negative potentials, correct ladder). An exit code was condemning good physics.

Fixed (09-07) by a virtual destructor. The three trapping invocations exit 0, and on an identical
160-time-unit run old and new binaries give **bit-identical output across all 16 datasets**. A Mac `obj/`
built before the fix traps; rebuild with `make` in place, not `install_M1.sh`, which empties `obj/` and
deletes the notes kept there.

### 3.10 The glpG blow-ups: the two subsystems are not at the same temperature (2026-09-10)

Diagnosed after 15 rollbacks in block 1 of the four production chains; the user's hypothesis (a
dry-MARTINI/Upside temperature mismatch) is confirmed. A dt scan was impossible: `apply_langevin_step`
throws when the runtime dt differs from `/input/brownian numerical_time_step`, and that check stays.

**K = T_up x 350.588235** (the table `~/OneDrive - The University of Chicago/image.png` checked on all 24
rows). The ladder ran from rung 0, 0.700 (245.4 K), to rung 27, 0.820 (287.5 K); the dry-MARTINI reference
is 0.8647 (303.15 K). Every dry-MARTINI parameter is built for 0.8647: `/input/brownian`
`reference_temperature_up` (friction fixed for 11.5 um^2/s; `bare_particle_diffusion_up` 5.1111 =
0.8647/0.16918), `py/martini_build_tables.py` `DEFAULT_PRODUCTION_TEMP_UPSIDE`, and the single-temperature
equilibration (`output_previous_0`). The bilayer ran **15.7 to 57.7 K below** that, so at rung 0 every
MARTINI interaction is 0.8647/0.70 = **1.24x stronger in kT**, over-condensing the bilayer. The scales are
not interchangeable: 0.7-0.9 is a reduced folding range for Upside's statistical potential, but the same
number is a literal kT in the Brownian noise.

**Measured, protein and lipids sit at different temperatures.** On clean frames of 79ALA, all 28 rungs:
lipids at T_lip/T_nom = 1.020 flat; the protein at T_nom + ~0.08 T_up (about +28 K), 1.506 at rung 27. The
rolled-back rungs (27, 18, 17, 14, 26) have the largest excess; the same ~1.10 ratio appears in cold rungs
0-13 that never rolled back, and lipids in the same slots stay on target. `protein_kinetic` is the mean over
the 630 dynamic backbone sites (the 420 placed O/sidechain slots carry zero momentum and are outside
`n_dynamic`); a naive average over all 1050 `PROTEIN` slots (0.661 x T_nom) is wrong.

Reproduced locally on one system without exchange (production T = 0.7215, tau = 5: T_prot 0.7911 against
the cluster's 0.7912, T_lip 0.7327 against 0.7355), so the excess is intrinsic to the hybrid integration:

| | T_nom | T_prot | T_lip | T_prot/T_nom |
|---|---|---|---|---|
| tau = 1 | 0.7215 | 0.8487 | 0.7349 | 1.176 |
| tau = 5 | 0.7215 | 0.8197 | 0.7489 | 1.136 |
| tau = 20 | 0.7215 | 0.8600 | 0.7412 | 1.192 |
| tau = 5 | 0.8647 | 0.9821 | 0.8796 | 1.136 |

The excess is multiplicative and independent of the thermostat timescale (20x range) and of temperature
(1.136 at both), so not a power leak, which would scale as P*tau.

**It is the mass-1 backbone, under either thermostat:**

| group | T/T_nom at 0.7215 | T/T_nom at 0.8647 |
|---|---|---|
| backbone, friction > 0 (g-JF Brownian) | 1.092 | 1.077 |
| backbone, friction == 0 (OU thermostat) | 1.130 | 1.072 |
| lipid/ion (g-JF Brownian) | 1.024 | 0.999 |

No trend with lipid-contact count (1.04-1.18), so the interface friction and the OU thermostat are
exonerated: integrator discretisation bias on the mass-1 sites. Not the springs (the stiffest,
`Spring_angle` k = 175, gives (omega*dt)^2/4 = 0.7% at dt 0.009); the steep MARTINI pair core is the only
curvature large enough, and it reaches the protein through `martini_hybrid_position`.

**The cold bilayer presses the backbone into that core** (minimum-image backbone-environment distances on
raw coordinates, a comparative proxy):

| | mean of per-frame minimum | closest seen | pairs/frame < 3.40 A |
|---|---|---|---|
| T = 0.7215 (production) | 3.512 A | **2.887 A** | **0.38** |
| T = 0.8647 (design point) | 3.623 A | 3.289 A | 0.07 |

5.4x more sub-3.40 A contacts at the production temperature. At 2.853 A dt for a 1 A kick on a mass-1 site
is 0.0028 (3.2), so dt = 0.009 is 3.2x too large. The causal chain: the bilayer runs 16-58 K below its
design point; MARTINI is up to 1.24x stronger in kT and over-condenses the environment; backbone sites are
pressed into the steep core; mass-1 sites at dt 0.009 integrate it inaccurately (standing 8-13% kinetic
excess, bursts to 1.5-1.8) until one site is ejected and tears TM4's backbone. The last link is inferred:
local runs (250 time units) reproduced the temperature split and the contact shift but no blow-up.

**The blow-up is a TM4 backbone tear, not a table defect.** A finite onset frame (79ALA r27,
`output_previous_2`, frame 208, -19765 -> +1683 -> +43240 E_up) puts the excess in `Spring_bond` (1272 ->
21724 -> 62558; `Spring_angle` +655, `Spring_omega` +332, MARTINI flat), in consecutive backbone bonds of
residues 139-141 (`C140-N141` at 18.5 A against r0 = 1.300). The local engine on the older tables reproduces
the recorded total to **0.27%** (+43356 vs +43240), so the arm-R retraining is not the cause. Frames before
onset already have peptide bonds at 2.0-2.9 A (about 60 kT of strain; equipartition at T 0.82 with k = 48
allows 0.13 A rms). Also ruled out: an unthermostatted subset (`stochastic_mask` is set only where
`friction > 0`, `martini_brownian.cpp:78`, and the OU thermostat reaches every friction-zero atom,
`thermostat.cpp:31`), and exchange laundering of the protein excess.

**Lipid diffusion cannot be measured from production output.** Exchange swaps configurations into a slot
every interval, so the apparent PO4 lateral D falls with lag (9.07, 5.14, 2.56, 1.43, 0.82 A^2/time_up at
lags 1-16), while the exchange-free equilibration rises (0.025 -> 0.103). Use a continuous
single-temperature run for any transport observable.

### 3.10a The fix, what it verifiably does, and what it does not (2026-09-10)

The ladder was made authoritative and dry-MARTINI brought to it (decision of 2026-09-10).

**Change 1: friction follows the replica temperature** (`src/martini_brownian.cpp`). With
`gamma = kT_ref/D_target` at T_ref = 0.8647, rung 0's lipid diffusion was 19% below target; the runtime
scales gamma by `T/T_ref`, giving `D = kT/gamma(T) = D_target` at every rung and through exchange. Keyed on
`reference_temperature_up` (required since A177). Transport only: friction does not enter the Boltzmann
distribution. Verified: bitwise at T = 0.8647 over 200 steps, different at 0.7215 (max |dpot| 40.5 E_up).

**Change 2: `inner_steps` = 4 by default** (`py/martini_prepare_system_lib.py`, env
`UPSIDE_MARTINI_INNER_STEPS`): N substeps of dt/N with noise on dt_i, so FDT holds. The outer dt stays
0.009, so `numerical_time_step`, the friction/dt lock and the 40 ps-per-step clock are untouched; it
integrates the same equations more accurately. glpG-RKRK-79ALA from an equilibrated frame at T = 0.7215:

| inner_steps | backbone (Brownian) | backbone (OU) | lipid | logged T_prot | excess | wall |
|---|---|---|---|---|---|---|
| 1 | 1.092 | 1.130 | 1.024 | 1.101 | +10.1% | 262 s |
| 2 | 1.034 | 1.039 | 1.018 | 1.035 | +3.5% | 333 s (1.27x) |
| 4 | 1.010 | 1.011 | 1.009 | 1.011 | **+1.1%** | 538 s (2.05x) |

At N = 4 both protein groups and the lipids sit within ~1% of the set point. The cost ratio is `(N+1)/2`,
because the baseline already does two force evaluations per step.

**Non-change: MARTINI epsilons are not rescaled.** Solute tempering (`eps * T/T_ref`) is forbidden, because a
spline table must equal the published dry-MARTINI form. At rung 0 MARTINI is still 1.24x too strong in kT
and the 5.4x contact excess stands.

**Not shown locally: that this stops the blow-ups.** The unsafe window shrinks only as 1/sqrt(N) (from the
deployed `combined_energy_grids`, the separation below which one step throws a mass-1 site more than 1 A or
0.3 A):

| inner_steps | dt_i | r(kick > 1.0 A) | r(kick > 0.3 A) |
|---|---|---|---|
| 1 | 0.00900 | 3.228 A | 3.544 A |
| 2 | 0.00450 | 2.900 | 3.181 |
| 4 | 0.00225 | 2.607 | 2.865 |
| 8 | 0.00112 | 2.350 | 2.572 |
| 64 | 0.00014 | 1.705 | 1.869 |

N = 4 shrinks the window to 2.607 A, but long runs reach ~0.2 pairs/frame inside 2.43 A (3.2): covering
that population needs N = 8, the 1.78 A approaches N ~ 64. A direct test from the last clean frame of the
real 79ALA r27 event (frame 207), four seeds at N = 1 and four at N = 8, held in all eight and so carries no
information: at the cluster rate (15 rollbacks over 11.4M replica-steps) it sampled 1.8% of one expected
waiting time. Production data settle it (3.10d).

### 3.10b The 0.90 ceiling is well below the unfolding transition, and TM4's loss is local (2026-09-10)

`reports/GroupMeetings/0323/group_meeting_03_23.pptx` slide 8 (implicit membrane): glpG's thermal
transition is sharp between T = 1.05 and 1.10 (mean CA-RMSD 17.6 -> 28.4 A, Rg 23.1 -> 35.7 A), plateauing by
1.15; at T = 0.90 RMSD is 10.0 A and Rg 18.4 A. **T = 0.90 is a legitimate ladder ceiling.** Local wild-type
hybrid runs agree (CA-RMSD to seed over t = 500-1000 / protein Rg): fixed T 0.70 8.67 / 20.17 A; fixed
T 0.90 seed 1 8.38 / 19.02 A, seed 2 8.64 / 19.97 A; unfixed T 0.90 9.91 / 20.51 A.

So TM4's helix loss at T 0.90 is not thermal unfolding: in the run that lost most (seed 2: TM4a 0.945 ->
0.556, TM1 1.000 -> 0.646) Rg fell 21.2 -> 19.5 A and the potential -20776 -> -21192 E_up, a more compact,
lower-energy, less-helical state. A suggestion to lower the ceiling to 0.82/0.86 on the T 0.90 helix numbers
was wrong and is withdrawn. Limits: slide 8 measured RMSD and Rg only, not TM4's helicity, and used the
implicit membrane.

### 3.10c What actually caused TM4 to be unstable, and what the GLY maps do (2026-09-10)

| hypothesis | test | verdict |
|---|---|---|
| global thermal unfolding | CA-RMSD 8.4-8.6 A, Rg 19-20 A; transition at T = 1.05-1.10 (3.10b) | ruled out |
| GLY143 mid-helix alphaL flip | phi stays -67 to -88, h = 0.91-1.00 for the whole run | ruled out, it never flips |
| GLY133 / GLY49 cap sign | phi = +71 at GLY133 in the healthy T = 0.70 run (TM4b 0.991) | ruled out, no correlation |
| TM4a at the bilayer interface | \|z\| = 1.47 A, the most central segment (TM3 = 8.58 A) | ruled out |
| arm-R force-field tables | recorded blow-up reproduced to 0.27% on the older tables | ruled out |
| GLY Ramachandran mis-symmetrisation | controlled run, correct mirror vs buggy (below) | ruled out, correcting it is **worse** |

**Primary cause: the temperature mismatch (3.10, 3.10a).** The protein ran ~10% hot, driving the nominal
0.70-0.90 ladder at ~0.77-0.99 T_up, across TM4's fraying range while TM3 stays below its own. With
`inner_steps = 4` the cold end is healthier than the seed: at T 0.70 TM4a 0.999, TM4b **0.991**, TM1 1.000
against seed 1.000 / 0.933 / 0.952; pre-fix the cluster gave TM4_full 0.589-0.727 at the same rungs, with 0
of 60 replica measurements above 0.8.

**The glycine maps were not this defect.** With the retired symmetrised maps (wild type, T 0.90,
`inner_steps = 4`, two seeds per arm): an off-by-one mirror, which happened to favour alpha_R at GLY143 by
+0.571 E_up, held TM4a at 0.786; the correct, exactly neutral mirror gave 0.531; the raw library map favours
alpha_L (-0.644). The residual fraying at T 0.90 is set by how much right-handed preference glycine's map
carries, a force-field question. TM4 is the glycine-dense helix (3 in 17 residues; TM1 has none) and frays
from its N-terminal turn at the hot end. n = 2 per arm, within-arm spread comparable to the difference.

### 3.10d The fix eliminates the blow-ups: 434x fewer bad frames (2026-09-11)

Every production frame with a non-finite or positive potential, across all 112 replica files (4 variants x
28 rungs), skipping the rigid-protein equilibration group:

| | production frames | non-finite or positive | rate |
|---|---|---|---|
| pre-fix (archived `pre_tempfix_20260910/`) | 246 120 | 1 604 | **0.6517%** |
| post-fix (`inner_steps = 4`, ceiling 0.90) | 481 040 | **7** | **0.0015%** |

A 434-fold reduction on nearly twice the data (the pre-fix rate predicts ~3 135); rollbacks 15 in the pre-fix
block 1 against 0 in the 12 visible post-fix chunks. The ceiling was simultaneously raised 0.82 -> 0.90, so
this is not a temperature-lowering artefact. The 7 remaining frames match 3.10a: at `inner_steps = 4` the
kick radius 2.607 A does not cover the ~2.43 A approaches of long runs; `inner_steps = 8` would, at ~3.6x
cost. The rollback machinery absorbs the residual rate.

### 3.11 Every glpG TM4 test since the coverage nodes ran on the retired ff3.0's coverage tables (2026-10-07; FIXED)

**Defect.** `patch_glpg.py` (md5 `e640487c...`; `$P/training/`, `checks/r4_epochs/tm4_local/scripts/`, the
Mac copies) writes a force field's `rama_pot`, `hbond_energy/parameters` and rotamer `pair_interaction` into
a glpG seed, on its docstring's premise that the hybrid has no `hbond_coverage` or
`hbond_coverage_hydrophobe` node. The pristine seed has none; the live seeds have both as `rotamer`
arguments, and the patch left their `interaction_param`. The local seed `glpG-RKRK-79HIS.live.up` and the
four cluster live seeds (`popepopg_REMD_mdw2/seeds/glpG-RKRK-{79ALA,79ALA_S115T,79HIS,79HIS_S115T}.up`)
carry coverage and hydrophobe tables equal to the retired FF1-form ff_3.0 `sidechain.h5` (git `9723c12f`,
09-10) to 1.5e-6, beside the 09-30 ff_3.0's pair table and H-bond energies (-1.8777 / -1.8719 / -1.7981 /
-0.6171), against `$P/parameters/ff_3.0`. The 09-30 02:25 ff_3.0 deploy (backups `*.bak_pre_ff_3.0_20260930-022510`) replaced the map, H-bond
energies and pair table and left the coverage; ff_3.0's own coverage tables sum to 4992.7 against the seeds'
7821.6, so the 09-30 glpG validation (1.14-1.15) was of that mixture too. Those coverage tables differ from
ff2.1's by rms 1.94 against ff2.1's own rms 3.29, so every local TM4 force field, "ff21_released" included,
ran beside v1's coverage; 4.5 shows a pair table is right only beside the coverage it was trained with.
`validate_ff.sh` patches the same cluster seeds with the same script.

**Size** (rotamer node output at the seed's starting frame, `upside_engine`):

| force field | as tested | own coverage |
|---|---|---|
| bio_start (ff2.1 side chains) | -19.3 | -121.7 |
| b00 | -14.3 | -109.9 |
| b9_01 | -13.9 | -108.1 |

Total potential moves by the same ~100 E_up; which way this moves TM4 is not measured. Round 1 (FF1-form
ff_3.0) is the only TM4 result whose coverage matched its pair table; every 3- and 12-seed TM4 count in
1.21-1.22, ff2.1's 8 of 12 included, has the mismatch.

**Fix (user, 10-07 11:22).** The patch also writes `hbond_coverage/interaction_param <- coverage_interaction`
and `hbond_coverage_hydrophobe/interaction_param <- hydrophobe_interaction` wherever the seed has those
nodes, and its round trip checks them: md5 `70589119...`, in the Mac test, `$P/training/` and
`tm4_local/scripts/` (backups `.bak_pre_cov_20261007`). A computer with its own `scratchpad/ff3_local_test`
must copy the cluster's script and re-patch before using it (remote_jobs.md "Resume here"; the cluster
`tm4_local/README.txt`). A re-patched input differs from its pre-fix input only in the two tables (rotamer
energy -121.74 at the seed frame, as predicted); `validate_ff.sh`'s method gate passes as before (worst
3.6e-15). The pre-fix runs are kept as a record (`runs_precov_20261007/` on the Mac, `tm4_local/runs/` on the
cluster); each rerun differs from its pre-fix twin only in the coverage tables.

### 3.12 The side-chain pair term's angular parameter gradient mixes the two residues (2026-10-18; FIXED)

**Defect (master's engine).** `PosQuadSplineInteraction` (`src/bead_interaction.h`) stores, for types (a,b),
bead 1's angular spline in slots 0-14 and bead 2's in 15-29 (`quadspline_param_deriv`, :125-130), and
`is_compatible` (:211-214) requires p[b,a] to hold the same two splines swapped. `propagate_derivatives`
(`src/interaction_graph.h:546-548`) added each edge's derivative unchanged to both [t1,t2] and [t2,t1], without
the swap, so wherever a-b pairs occur in both sequence orders, [a,b] collects b's angular derivative in a's
slots. With `rotamer_parameter_estimation.py:60-72` building p[a,b] = [sig(L[a,b]), sig(L[b,a]), radial], the
trainer's gradient on L[a,b] is sig' (D_a + D_b) where it should be proportional to D_a.

**Measured** (1aba, central differences; `checks/nobp_audit_20261018/B/stage2.log`): the slot sum
c[a,b,s] + c[a,b,15+s] equals FD(a's spline) + FD(b's) to 1.9e-5 (single-state) and 4.6e-4 (multi-state),
while c[a,b,s] against FD(a's spline) alone ranges from -66 to +65. Trainer over true gradient on
off-diagonal angular logits: median 5.1 (single-state), 1.11 (multi-state, range -4046 to +894, some signs
flip). Diagonal angular and every radial coefficient come out at exactly 2x (2.000 +/- 5e-4); coverage and
hydrophobe 1.000.

**Scope.** 5,700 of the 31,420 `rot` coefficients (rate 0.0125) in every ConDiv run of this campaign (ff3.0,
ff3.1, ff3.2, the control and no-BP), probably ff2.1's original training too (same engine code; not checked).
Simulation energies and forces do not use this path. The misdirected logits moved as far as the radial
coefficients (`checks/angfix_20261018/ang_drift.txt`, `ang_drift_sigmoid.txt`; rms change, off-diagonal
logit / radial: control 0.083 / 0.092, no-BP 0.150 / 0.180, ff3.0 0.137 / 0.142, ff3.1 0.134 / 0.143, ff3.2
0.113 / 0.107), rms 0.007-0.012 in spline value (max 0.06-0.11) against a median value of 0.016. Its effect on
any force field is not measured. (The first reading, "barely moved", compared unlike units; under Adam a
coordinate's step size says nothing about whether its gradient is right.)

**Fix (user 10-18: "fix the issue"; plan.md Phase 21).** `interaction_graph.h` adds to [t2,t1] the parameter
derivative evaluated at p[t2,t1] with the two elements exchanged (`ExchangedParamDeriv`, symmetric
interactions only); for fully symmetric parameter sets it is bitwise the old value. On midway2
(`checks/angfix_20261018/`: `parity.log`, `parity_glpg.log`, `fd.log`): energy, forces and every derivative
ConDiv reads bitwise on single- and multi-state 1aba except the angular slots; 200 tu binary runs and the
hybrid glpG seeds bitwise. Every angular entry matches its tied central difference (ratio 1.000; median error
5.8e-6 single-state, 9.6e-4 multi-state), and the trainer's angular-logit gradient is 2.000x the true one, as
every radial coefficient is. That uniform 2x stays: it sets the regulariser's relative weight, as in ff2.1's
own training. Installed in `$P/obj` 08:35 10-18; ff3.0, ff3.1, ff3.2 and both ff2.1 retrainings before that
carry the defect.

### 3.13 `environment_coverage_hb` typed each covering side chain with another residue's type (2026-10-19; FIXED)

**Defect (master's code).** `write_bb_environment` (`py/upside_config.py`) stored `aa_types` as one type per
H/O site (`n_donor + n_acceptor` entries), while `EnvironmentCoverage` (`src/environment.cpp:325, 336`) reads
`aa_types[id2]` with `id2` the covering side chain's residue. Side chain j took the type of site j (donor
residue j+1, shifted one more per preceding proline, or for the last 1 + n_PRO residues an acceptor near the
N-terminus): only 7-9% of covering side chains got their own type (2LZM, 2RN2, 1YOB, 1ENH). The
backbone-desolvation burial `bl = sum_a w[a] cov[a] + 2 hbcov` therefore summed scrambled weights (ff_2.1's
w runs 0.56-2.12), independent of size. The side-chain burial node already stored per-residue types.

**Fix (user 10-19: "Fix the bug"; plan.md Phase 22 item 8).** Both coverage writers store one type per
residue; the engine sizes `aa_types` from its dataset and requires every `id2` to index it. A file written
before the fix runs as master; an engine without the fix refuses a new file ("dimensions of 'aa_types',
expected (105) but got (54)"). Measured (`scratchpad/large_protein_20261019/envfix/`):
* pre-fix files on the pre-fix and fixed builds: energy, every term, both coverages and the derivative
  bitwise (8 configs x 3 frames, 208 arrays); new and old configs differ in `environment_coverage_hb/aa_types`
  only;
* both builds match a numpy rebuild of the coverage to 5e-7, each with the types it should use;
* finite differences: ratio 1.0001 at h 0.005 (40 coordinates, 2LZM and 1YOB globules);
* ff_2.1's energy changes only in backbone desolvation: -0.1 to +0.2 E_up at natives, +0.05 to +1.3 at
  globules;
* 5,000 tu native starts: KE/(1.5 kT) 1.007-1.010, C-N 1.33-1.34 A, Q and Rg within single-run noise of the
  pre-fix runs.

Installed in the repo's `obj/` only (backups `*.bak_pre_envfix_20261019`). `$P` and `/beagle3` keep the old
writer and engine while trainings run; ff_2.1 and every training so far carry the defect. Upside 3 trains
with the fix (`upside3.md` section 4).

### 3.14 Upside 3 audit ledger (plan.md Phase 24; from 2026-10-19)

Confirmed defects, fixed on martini-dev (Upside 3) only. Each fix is verified three ways:
* its reproducer passes;
* `tests/node_fd.py` passes on the node;
* `tests/regress.py compare` shows only the intended arrays changed.

Effects on ff2.1 are not measured (upside3.md section 1).

| id | file:line | in master | defect | effect | verification |
|---|---|---|---|---|---|
| A1 | `coord_basic.cpp:343` `DihedralCoord` | yes | loads the dummy atom -1 of terminal phi/psi before testing for the dummy | out-of-bounds read; the value was unused | ASan clean; regression bitwise |
| A2 | `membrane_potential.cpp:477, 978` `MembraneCBPotential`, `MembraneSurfCBPotential` `get_param_deriv` | yes | the right-hand spline's derivative was written to `shift + np`, one table past the end of `deriv` (the energy uses the same coefficients for both sides) | half the membrane coefficient gradient lost; the write corrupted heap memory (other nodes' outputs, NaN in 2 of 6 runs). Only membrane-potential training called it | 83 and 97 of 220 directions failed before; FD median error 1e-6 after; ASan clean; only `cb_membrane_potential/param_deriv` changes |
| A3 | `main.cpp` `ReplicaExchange::attempt_swaps`, `--exchange-criterion` | yes | criterion 1 compared the two replicas' Boltzmann factors before the swap, which is not a Metropolis exchange rule | wrong ensemble whenever it was selected: z up to 12 in mean potential; nothing in the repo selected it | removed with its option and the `run_upside.py` passthrough; regression bitwise |
| A4 | `main.cpp` main loop `do_break` | yes | the `nr>last_start` test (older than the `(nr+1)` fix) made a 1-round replica interval swap every 2 rounds | swap frequency halved at interval 1 | attempts now exactly rounds/interval at 1, 2 and 10 rounds |
| A5 | `main.cpp` `inner_step` | branch | `inner_step` was set from `--duration-steps` rather than the integrator | protein runs under `--duration-steps` logged 1/3 of their time, with every interval 3x short. g-JF hybrid runs under `--duration` logged 3x their time (MARTINI workflows use `--duration-steps`, so production was right) | `inner_step` comes from `/input/brownian`; every system must share it (validated precondition); protein 300 rounds log 8.1 = 300 x 3 x 0.009; hybrid 10 rounds log 0.09; regression bitwise |
| A6 | `main.cpp` swap acceptance | branch | `isfinite(lboltz_diff)` guard, compiled out under the Mac flags and active on the cluster | a guard under NO GUARDS, and platform-dependent behaviour | removed; regression bitwise |
| A7 | `main.cpp` `coord_swap` | yes (master swaps coordinates only) | the integrator leaves momenta half a step behind the coordinates, so a coordinate-only swap breaks the momentum-force correlation | master's REMD heats the replicas: KE/1.5kT 1.011-1.013 at a swap every 37 rounds and 1.027-1.039 at every 10, against 1.002-1.006 for independent runs | this branch's momentum swap with sqrt(T_dest/T_src) rescaling (2026-08-05) measures 1.003-1.007; its comment, which gave the wrong reason, rewritten |
| A8 | `interaction_graph.h` `change_cache_buffer` | yes | the cache check uses `sqr(0.5*(cache_cutoff-cutoff))`, so a `set_param` that raises the cutoff kept a stale pair list | wrong energies after an in-process parameter change (1dfn coverage r0 8 to 12: -14.61 against -12.84 fresh) | a new buffer invalidates the cache (reviewer: 0 mismatches over every case); simulations bitwise |
| A9 | `environment.cpp:769` `NonlinearCoupling` | yes | assigned the coverage sensitivity instead of adding | another consumer's force on the same coverage lost; no writer combines them | `+=`; regression bitwise |
| A10 | `membrane_potential.cpp:1346`, `upside_config.py:2029` `membranelateral_potential` | yes | the writer mapped the profile's z range onto spline coordinate [0,1]; the spline's samples sit at 1..n | the lateral-pressure term read one clamped value at every z: zero force | the spline reproduces the input profile; forces 0.06-0.08 where they were 1e-10 |
| A11 | `main.cpp` per-system seed | yes | `seed + ns` | system 1 of seed s ran the noise streams of system 0 of seed s+1 | `system_seed()` hashes (seed, ns) with splitmix64; seed 5 system 1 and seed 6 system 0 now differ by 1.6 A; only the dynamics arrays change |
| A12 | `box.cpp` barostat RNG; `main.cpp` replica exchange | branch | the barostat RNG was seeded from the config's `mc_seed`, the same in every replica and independent of `--seed`; replica exchange swapped coordinates without the box | identical barostat noise in every replica; mixed boxes in NPT replica exchange (latent: production replica exchange runs NVT) | seeded from the system seed (`mc_seed` and `UPSIDE_NPT_MC_SEED` removed); replica exchange refuses the barostat as a validated precondition, shown on 1rkl stage 6.0 |
| A13 | `upside_config.py` trans-cis writers | yes | `seq_new2 = seq_new[:]` is a numpy view, so the trans map was built from the cis sequence; the proline arrays were partly in file order and partly residue-sorted | identical trans and cis maps (true difference up to 4.75 E_up on 2qke); energies depended on the input's line order (3.2-3.9 E_up) | one sorted reader `read_pro_state_file`, used by both writers, and a copy; written maps equal the expected trans and cis maps exactly; sorted and unsorted inputs give equal energy |
| A14 | `upside_config.py` `--fix-rotamer` | yes | `argmin` of a signed periodic chi2 distance | picked the farthest-below state (22 of 51 chi2 choices wrong on 1UBQ) | `argmin(abs(...))` |
| A15 | `upside_config.py` `--rot-scale` | yes | `[2*n_knot_ang:]` sliced axis 0 of a (20,20,54) table | the option did nothing | radial knots now scale exactly by the factor, angular knots unchanged, in all three tables |
| A16 | `upside_config.py` `read_rama_maps_and_weights` | yes | neighbour-dependent maps ignored chain breaks | a chain's first and last residues took the other chain's residue as neighbour (1dfn: 1.2-1.9 E_up per map) | the maps take `chain_first_residue`; dimer maps equal each chain's own; single chains bitwise unchanged |
| A17 | `upside_config.py` `write_infer_H_O`, `main`; `martini_prepare_system_lib.py` | yes | both residues at each break were excluded as donor and acceptor | each break lost a valid NH and C=O (1dfn: +4.08 E_up H-bond energy) | `hbond_residues()` excludes only the sites spanning the break; separated dimer's H-bond energy equals A+B (-68.8280 against -68.8283) |
| A18 | `hbond.cpp:644` | yes | each bond was counted at both its donor's and its acceptor's residue | the printed "hbonds" was twice the real count (stdout only) | halved |
| A19 | `rotamer.cpp` `calculate_n_elem` | yes | accepted 2, 4 or 5 rotamers, which have no node holder | segfault at the first energy call | throws unless 1, 3 or 6 |
| A20 | `rotamer.cpp` `compute_value` | yes | `n_bad_solve` counted solves that ended near the iteration budget | a solve that converged at iteration 18 was counted bad for max_iter 19-21, and a NaN deviation never was | counts a final deviation that is not `<= tol` |
| A21 | `upside_config.py` `write_rotamer_placement` | yes | the path without `--dynamic-rotamer-1body` reads `rotamer_prob_fixed`, which no shipped library has | `NoSuchNodeError` | a clear error naming the option |
| A22 | `upside_engine.py` `get_rotamer_graph` | yes | asked for values the rotamer node does not log | always raised; nothing called it | deleted |
| A23 | `placement.cpp`, `environment.cpp` loggers | yes | every instance of a node type logged under one fixed name | `--log-level extensive` failed on every standard config | named after the node's group (`group_name()` in `h5_support`); extensive runs on chig and 1qhj_memb, one dataset per node |
| A24 | `upside_config.py` (`n_atom/3`, `bb_index/3`, `hbbb_coverage` id2) | yes | Python 3 true division wrote index datasets as floats | right only because HDF5 truncates on read | integer division; engine record bitwise unchanged |
| A25 | `advanced_config.py:560` `make_restraint_group` (O1) | yes | `x/3 != y/3` is float division in Python 3 | same-residue pairs were restrained, including in ConDiv's native-restrained replica | `//`; 1ubq group: 445 pairs, none within a residue |
| A26 | `interaction_graph.h` `PairlistComputation` | yes | 4-wide stores into a buffer of `max_n_edge` with no count check | heap overflow when an `intensive_memory` budget is exceeded (compressed 1ubq) | throws "pair list needs more than its max_n_edge" first (validated precondition); regression bitwise |
| A27 | `eig.cpp:357` `AffineAlignment` | yes | eigensolver tolerance 1e-5 left frame errors near 5e-5 that switched on and off as residues crossed it | rotamer energy jumps of 1-4e-3 E_up the force does not see | tolerance 2e-7: largest path step 5.8e-4, path deviation at most 2.1e-4; quaternions change by at most 4.4e-5 after sign, forces by 5e-5 RMS; runtime +2% |
| A28 | `src/CMakeLists_M1.txt` | branch (master lacks it on both) | the Mac build used `-ffast-math` without `-fno-finite-math-only`, which the cluster build had | NaN and infinity tests compiled out on the Mac only, so the two platforms behaved differently | flags now match; static arrays change by at most 1.3e-5 relative |
| A29 | `environment.cpp` `BackboneSigmoidCoupling`; `upside_config.py` `write_bb_environment`; both `ConDiv.py` | yes | the backbone desolvation term reads the side-chain burial's trained weights but exposed no weight parameters; its center, sharpness and hbond weight derivatives were commented out (D1) | ConDiv trained `envw` on a partial gradient: the missing backbone part was 40-130% of the side-chain part in a native-vs-perturbed contrast, in an unrelated direction (ff2.1's training included) | scale, center, sharpness, hbond weight and 20 weights, all exact (FD 0 bad of 44 per frame); the writer stores one weight per residue type and refuses a library whose rows differ (old 400-entry files read their first row); ConDiv adds the backbone weight derivative to `envw` (FD of a shared weight -0.582 against -0.586; side-chain part alone -1.340) and refuses `number_independent_weights` other than 20; bbenvc/bbenvs/bbenvw stay frozen |
| A30 | `hbond.cpp` `HBondEnergy::get_param_deriv`; both `ConDiv.py` | yes | entries 4-11 (rama boundaries and sharpnesses) returned 0 while sitting in the trained `hb` field | frozen by accident | exact derivatives (0 bad of 32 per frame, 35 with class offsets); the trainer's new field `hbr` holds them at learning rate 0 |
| A31 | `environment.cpp` coverage interactions | yes | `param_deriv` returned zeros, marked "not implemented" | a silent zero gradient for any trainer that asked | throw "the coverage parameter derivative is not implemented"; the harness records it |
| A32 | `environment.cpp`, `sigmoid_energy.cpp`, `hbond.cpp`, `sidechain_radial.cpp`; `advanced_config.py`, `run_upside.py` | yes | nodes and options that no writer creates or no reader accepts: SphereCoverage, WeightedSphereCoverage (heap overflow, one-sided counts), uniform_transform, linear_coupling_with_inactivation, SigmoidLEnergy (wrong derivatives), CTerHBond (no forces), cooperation_contacts (writer and reader disagree), `--plumed` and the external pair tables (no reader), and a commented-out legacy block | dead or broken code | deleted (user 10-19); regression bitwise |
| A33 | `bead_interaction.h` angular splines | yes | the cosine of two unit vectors rounds past -1 (norms exceed 1 by up to 2.5e-7), which read `p[-1]` and wrote `d_param[-1]` | rare out-of-range read and write, the write clobbering a neighbouring edge's slot in training | `angular_coord()` holds the cosine to [-1,1] (user 10-19); in-range values are unchanged; forces move by at most 6e-7 relative where it applies |
| A34 | `rotamer.cpp` BP (R2) | yes | 1e-10 floors applied to probabilities that were not normalised: single-rotamer edges were folded in after conversion and edge probabilities never rescaled | below-Bethe free energy and wrong marginals in severe clashes (-74 E_up at a pair shift of 30, NaN at 40); never reached at ff2.1 natives | edges hold energies, single-rotamer edges fold in before conversion, each edge carries its log scale, floors 1e-30. Clash cases equal float64 Bethe to at most 5e-5 (one residual, +0.60 at shift 40 on multi-rotamer pairs, needs a product-form cavity); forests exact to 7e-7; normal frames within the BP tolerance |
| A35 | `membrane_potential.cpp` `MembraneSurfHBPotential::get_param_deriv` | yes | the derivative multiplied by `surfv` alone where the energy uses `surfv*compact_sigmoid(bl-2.5,2)` | wrong surface-membrane H-bond gradient (63 of 92 slots failed, ratios 1.1-3.7) | uses the energy's factor; 0 of 92 bad |
| A36 | `membrane_potential.cpp` `MembraneLateralPotential` | yes | with `--membrane-exclude-residues`, arrays indexed by position were read by residue id | heap overflow (SIGBUS) and other residues' burial | consistent indexing; ASan clean; forces on kept residues bitwise equal to the run without exclusions |
| A37 | `main.cpp` `CurvatureChange`; `random.h` | yes | multiplicative proposal without its Jacobian, with irreversible support; unseeded `rand()` shared across threads; acceptance drawn on the replica-exchange stream with the base seed; uninitialised centre pointer; stale energies | sampled a flat prior in curvature by accident (chi2 1835/10 against the energy's measure), ignored `--seed`, raced between replicas, segfaulted without the centre node | step symmetric in log\|z\| with the Jacobian for a measure flat in 1/R (a flat-in-R measure runs away to R = 1e8 A); own stream keyed by the system seed; centre found at setup or refused; fresh energies; exact restore: chi2 8.7/10 against exp(-E/T); bitwise per seed |
| A38 | `surface.cpp` `FORCE_UNIT` | yes | 2.428602e-5 is 1 bar in units of kT at 298 K, not E_up | lateral pressures 17.6% too strong | 2.065948e-5 E_up/A^3, derived in the comment |
| A39 | `membrane_potential.cpp` curvature gradients | yes | assumed the curvature centre on the z axis | wrong forces for an off-axis centre (latent: the writer puts it on the axis) | `xyz - ccenter`; off-axis force FD 0 bad; on-axis bitwise unchanged |
| A40 | `main.cpp` re-centering checks | yes | `is_prefix` arguments reversed and node names out of date | z-recentering was never refused: the membrane fixture's protein z centre was reset every frame under default flags | `is_prefix(type, name)` with the writers' current names (every membrane variant, surface, lateral, walls, tension, pulling, tethers, `restraint_position`), the cavity for any re-centering, and fixed atoms; the membrane fixture now refuses default flags, chig runs bitwise as before |
| A41 | `upside_config.py` `write_membrane_potential`, `--membrane-exposed-criterion`; `run_upside.py` | yes | the old membrane writer raised KeyError under Python 3, and the option only fed it | dead | deleted; the C++ node stays registered so old configs load (the node itself deleted by A229) |
| A42 | `main.cpp`, `thermostat`, `martini_brownian`, `box.cpp`; MARTINI workflow writers; example 13 | yes | every random counter (round, thermostat, Brownian, MC, barostat) restarted at 0 on a continuation | a same-seed continuation replayed the previous segment's noise (chig centre-of-mass walk correlated 0.999; hybrid protein displacements 0.76); the MARTINI workflow restarts stages with one seed | `/output/round` logged, `/input/mom@round` resumes every counter: chig continuations bitwise the continuous run (plain, off the frame grid, with pivots, under `mv`), the hybrid Brownian path bitwise; momentum without the attribute (pre-10-19) restarts at round 0 |
| A43 | `monte_carlo_sampler.cpp` | yes | an accepted pivot or jump moved the coordinates but left the half-step momentum correlated with the old force | MD heated: KE/1.5kT 1.040 with a pivot every round | half-kick correction on acceptance (`inner_step*dt/2` for the slow level under `mv`): 1.007 against 1.004 without pivots; pivot detailed balance unchanged (chi2/dof 288/291) |
| A44 | `main.cpp` annealing | yes | without `--anneal-end` the schedule ended at its start | a step from T to T*factor | `anneal_end = anneal_start + anneal_duration`: 0.8, 0.754, ..., 0.4 |
| A45 | `main.cpp` initial thermalization | yes | applied after the real thermostat interval was set, as a partial refresh | KE/1.5kT 0.0097 at frame 0, 0.63 over the first 5 time units | full refresh: 0.903 at frame 0 |
| A46 | `martini_masses.cpp`, `thermostat.cpp`, `main.cpp` | branch | masses were not validated, with silent m=1 fallbacks | a zero mass gave NaN positions and exit 0; a short mass array ran with m=1 for the rest | one positive mass per atom required; fallbacks removed |
| A47 | `main.cpp` logged kinetic energy | branch | fixed atoms were counted | KE/1.5kT read 0.902 with 3 of 30 atoms fixed | averaged over dynamic atoms |
| A48 | `deriv_engine.cpp` Predescu `integration_cycle`, `integration_stage`; `src/pivot_sampler.h` | yes | the kick order was not palindromic (not reversible); the header did not compile; neither had a user | dead | deleted |
| A49 | `main.cpp` small items | yes | MC, replica and curvature intervals truncated with `int()`; logged time in float; `--temperature` help claimed a default; final KE report used the final temperature under annealing; `--minimize` wrote a duplicate frame; MC moves moved fixed atoms | an MC interval of 0.75 at dt 0.025 ran every 9 rounds; time off by 0.013 at 5.4e5; etc. | `round()`; double; help text; per-frame temperature; no duplicate; MC refused with fixed atoms |
| A50 | `interaction_graph.h` capacity check (follow-up to A26) | branch | the check compared against `max_n_edge`, stricter than the buffers' `round_up(max_n_edge,4)` | every energy call on a tripeptide threw (steric pair list capacity 3) | compares against the buffer size: tripeptide 25.657316 as master; compressed intensive-memory 1ubq still refused |
| A51 | `coord_operator.cpp` `Multiply` | yes | `propagate_deriv` scaled each input's sensitivity by its own factor instead of the other's | the cis-proline Ramachandran force in trans-cis runs was off by the factor lambda_cis (work 0.837 against an energy change of 3.346 across the transition) | correct factor per input and dimension; fed-twice elements refused; FD and path work agree (3.34631) |
| A52 | `upside_config.py` `read_pro_state_file` | yes | a switching proline at residue 0 or a chain start has no omega dihedral, but the switch ids assume one per proline | that proline read another proline's switch from allocation padding | refused with a clear message |
| A53 | `backbone_steric.cpp`, `upside_config.py` `write_backbone_pair` | yes | the steric wall excluded neighbours by residue number across chain breaks | no wall between a chain's last and the next chain's first residue (2.5 A apart: 0.000 against 32.0) | optional `sequence_id`, offset by 2 per break, replaces `id` for the exclusion only; single chains bitwise unchanged |
| A54 | `main.cpp` moving constants | yes | `MovingConst3D` advances by its own `time_step`, not the integrator's dt | the AFM tip moved 1.8x too fast at `--time-step 0.005`; the help gave the velocity per time step | refused when they differ; help gives A per time unit |
| A55 | `coord_const.cpp`, `coord_map.cpp`, `coord_operator.cpp`, `upside_nodes.py`, `advanced_config.py`, `upside_config.py`, `up.md` | yes | broken nodes nothing writes (MovingConst1D/2D and WhirlingConst1D write out of bounds; concat segfaults; weight_pos truncates weights to int; sum_of_values miscounts repeated indices; Sum indexes the weights wrongly) and writers their readers reject (`write_angle_spring2`, `write_rama_coord`, `ConcatCoord`, the always-off `append_const3D_to_pos`) | dead or broken code; `Const1D/2D` writers raised on numpy arrays | deleted; `if values is not None` in the kept writers; up.md rows for `rama_coord`, `backbone_pairs` and `rama_map_pot_ref` corrected; regression bitwise |
| A56 | MARTINI per-engine registries (`martini_masses`, `martini_fix_rigid`, `martini_stage_params`, `martini_hybrid`, `martini_brownian`, `box.cpp`); `engine_c_library.cpp`; `main.cpp` | branch | state in global maps keyed by the engine's address was never removed when the engine was destroyed; stage registration wrote only when the config had a stage | an engine built at a freed engine's address inherited its stage: -13016 E_up where the fresh value is -18316 (4 of 4 reuses) | `martini_unregister_engine()` erases every module's entries, called by `free_deriv_engine` and at the end of `upside_main`; registration resets first; reuses give the fresh energy |
| A57 | `engine_c_library.cpp` `set_param` | yes | no size check | a half-length vector overflowed the heap in the membrane nodes; parameter-less nodes accepted anything | refuses a wrong size, as `get_param` does |
| A58 | `main.cpp` `upside_main` | yes | `received_signal` was never reset | after one interrupted in-process run, the next did 0 rounds and raised KeyboardInterrupt | reset at entry |
| A59 | `engine_c_library.cpp` `get_sens`; `h5_support.cpp` `read_attribute`; `timing.h` | yes | `get_sens` on a potential returned the potential; an array attribute was read into a scalar; the timer counted the warm-up call it divided out | wrong value; stack overflow on malformed input (ASan); profile averages wrong | 1 (dE/dE); refuse a non-scalar attribute; skip warm-up as the divisor assumes |
| A60 | `src/hmm.cpp`, `py/tensorflow_upside.py`, `py/predict_chi1.py`, `upside_engine.freeze_nodes`, `h5::open_file`, `h5::ensure_not_exist`, the Mac `upside_engine` MODULE target and its Python lookup | yes / branch | not compiled; Python 2 and unimportable; broken (bytes against str, nonexistent `constant_*` nodes); unused; a module without `PyInit` | dead | deleted; regression bitwise |
| A61 | `py/martini_hdx_membrane_accessibility.py` | branch | 0-based donor indices were matched to PDB residue numbers | glpG (numbered from 1) gave each amide its neighbour's lipid accessibility with no error; 1rkl (from 4) raised | donor to amide through `/input/hybrid_bb_map`, as `martini_hdx_project` does. Arrays made before 10-20 and the HDX outputs built on them must be regenerated (plan.md Known Errors) |
| A62 | `example/03` MBAR, `helpers/calc_hdx_ht.py`, `example/06` force | yes | kT taken as 0.593 kcal/mol per T_up (298 K) and 41.4 pN per E_up/A (300 K) | free energies 14.9% low, pulling forces 14.5% low | one convention, T_up 1 = 350.588 K: kT 0.6966 kcal/mol per T_up, 48.40 pN per E_up/A; the extension labelled in A |
| A63 | `py/get_info_from_upside_traj.py` | yes | summed acceptor and donor columns of `/output/hbond` | the H-bond count doubled | donors only |
| A64 | `example/00` HDX: `4.calc_D_uptake.py`, `helpers/calc_hdx_ht.py` | yes | PF 1000 for a fully protected residue; pKD2O at 298 K beside a 400 K Arrhenius factor; no check that the target T lies on the ladder; the dG 1000 sentinel passed `isfinite`; cis-proline `*` left in the sequence | the most protected residues exchanged fastest; base rate off by 219x at T_up 1.14; R2 and the best-T choice corrupted by sentinels; IndexError | the ESS bound; pKD2O at the Arrhenius T; refused outside the ladder; NaN for censored dG; `*` stripped. The default `legacy_T_range` 1.14 now stops ladders that end below it |
| A65 | `py/extract_vtf.py`, `py/mdtraj_upside.py`, `helpers/select_h5_frames.py`, `py/martini_extract_vtf.py` | yes / branch | frame selection off by one across restarts; restart time offset added twice since A42; duplicate restart frames; datasets paired by position; H written on cis-proline and O/H placed across chain breaks; NPT frames unwrapped with the last box; `load_upside_ref` returned None | misaligned frames (up to 15.6 A), wrong times (2997), wrong VTF atoms, torn bonds, crash | each fixed as in `mdtraj_upside`; the dead helpers duplicate `get_info`, `write_hybrid_energy`, `estimate_rama_distributions.py`, `print_traj_vtf`, `load_upside_rep` and the stubs deleted |
| A66 | `main.cpp` re-centering checks (follow-up to A40) | yes | xy-only re-centering was allowed with terms on absolute x and y | a 5 A x-shift changed an xyz nail's energy by 1158 and a curved membrane's by 0.63 while the run proceeded | any re-centering is refused with tethers, walls, pulling, restraints, the cavity, or a membrane with `use_curvature`; a flat membrane still allows xy-only re-centering |
| A67 | `martini_fix_rigid.cpp` rigid-body fit | branch | Horn's correlation matrix built transposed (the inverse rotation); 16 power iterations, which find the largest-magnitude eigenvalue, not the largest | in prep 6.x and the 7.0 burn-in the protein turned at 0.29, 0.79, 0.92 of its angular velocity about its axes; NVE drift -19 to +37 E_up per 3000 steps; 6.6 minimisation rotated it 0.0005 rad | S(a,b) = sum m r_a c_b as in `AffineAlignment`, exact `SelfAdjointEigenSolver`, one projection routine for MD, minimisation and placement, masses from `martini_masses`, singular inertia refused: rotation capture 1.00 per axis, drift -0.74 E_up, 6.6 rotates 0.085 rad; production bitwise |
| A68 | `martini_hybrid.cpp`, `martini_stage_params.cpp`, `main.cpp` hybrid stages | branch | the interface was active only when the stage label matched; `--minimize` switched every config to the "minimization" stage; an unknown label was accepted | plain `--minimize` on a production config, or `current_stage` "Production", turned the SC-env and BB-env interface off and froze the protein at exit 0 | each `activation_stage` runs only in its own stages (minimization: minimization; production: production_handoff, production) with the interface on, anything else refused; `--minimize` keeps the config's stage; `--minimize-preserve-stage` and the stage switch deleted; up.md stage table |
| A69 | `martini_fix_rigid.cpp` reader; hybrid dead code | branch | `catch(...)` around the fix_rigid read and silent skips of out-of-range indices; unused flags and fallbacks (`exclude_intra_protein_martini`, `virtual_backbone_com_mode`, `n_env`, `bb_residue_index`, the BB inference fallback, an overwritten force projection) | a wrongly shaped or out-of-range `fix_rigid/atom_indices` fixed nothing, exit 0 | refused at registration; `bb_atom_index` required; dead code deleted (poisoning the projection's derivatives with 1e30 leaves the trajectory bitwise); up.md `hybrid_bb_map`, `hybrid_control`, hybrid `pos` layout corrected, `hybrid_remap` removed |
| A70 | hybrid prep overlaps: `martini_prepare_system.py`, `martini_prepare_system_lib.py`; `martini_potential.cpp` `martini_run_minimization` | branch | CHARMM-GUI's hard minimisation off by default and the soft core (alpha 0.2) a flat plateau, so template overlaps survived into MD; the minimiser stepped in proportion to the gradient | 6.3 at KE/1.5kT 26-33 (1rkl x4, glpG), box +6-7 A; a 6.4 start at 1.1e13 E_up (two beads at 1.5e8 A) carried into the 7.0 seed; turning 6.1 on threw beads 1e14 A (NaN placements, heap overflow, SIGBUS) | CHARMM-GUI's order: 6.0 soft-core minimisation (lambda 0.01), 6.1 hard; a pair above kT ln(1e9 N_pairs) refused before each MD stage; GROMACS steep minimiser (bounded step, 1.2x / 0.2x): max pair energy at MD starts 1.1-7.9 kT (1rkl), 1.4-4.8 kT (glpG); 6.3 KE/1.5kT 0.64 / 0.58, 0.97-1.00 by 6.4 |
| A71 | PO4 restraint, `martini_prepare_system_lib.py`; prep z freeze, `martini_hybrid.cpp` | branch | the head restraint's reference was in the uncentred PDB frame while the engine also froze PO4 z | every PO4 53.23 A off its reference, adding 250,773 E_up to the logged potential; the freeze made the 200-to-10 ramp a no-op and pinned P-P at the template's 36.4 A | one mechanism: the harmonic z restraint referenced in the stage frame, ramped 200, 100, 50, 20, 10 kJ/mol/nm^2 over 6.2-6.6, none and no freeze at 7.0; prep freeze and `preprod_lipid_headgroup_roles` deleted; 6.2 logs real potential + 270.21 E_up (analytic 270.206) |
| A72 | hybrid prep composition and topology: ions, rotamer 1-body, ITP mapping | branch | ions neutralised residue charges the runtime protein does not carry (findings 5.3c); without `rotamer_prob_fixed` (no library has it) the rotamer 1-body fell back to -log of the unweighted phi/psi-grid average; ITP bonds mapped by bead position, out-of-range indices skipped | net charge +3 e (1rkl), -6 e (glpG); helix rotamer energies off by up to 9.1 E_up (LEU); a reordered lipid bonded NC3-GL1 at 7.04 A | ions from the charges the topology carries, net 0 asserted; `placement_scalar` on `rama_coord` as `--dynamic-rotamer-1body` (equal to the pure ff2.1 config on 117 rows); ITP index to bead name to PDB index, mismatches refused |
| A73 | hybrid prep restarts and tables | branch | continuations position-only at round 0 with the same seed, the round dropped through minimisation and handoffs, one seed for prep and REMD; soft prep tables with bare 1/r Coulomb; continuations re-injected tables from disk; the protein restraint held the virtual O/BB slots | segment displacement correlation 0.220; prep tables -2.6 E_up at 12 A where the potential is 0; a continuation would switch force field mid-trajectory | momentum and round carried, per-stage splitmix64 seeds (a split production bitwise the continuous run); soft tables the published RF + shifted LJ with softened r (lambda 0 bitwise, 0 at 12 A); continuations copy `/input`; N/CA/C restrained only; examples set POPE/POPG, midway2 args, T 0.8647, `env_shared.sh`; dead micelle, water and bilayer-mode paths deleted |
| A74 | `box.cpp` MC barostat | branch | a molecule with any z-held atom, and the protein (derived BB/O slots counted as fixed), was never scaled; molecule centres from atoms in different images; volume term N ln V; nodes got an incrementally rescaled box | prep NPT never moved a lipid or the protein (5e-8 A per 400 calls); a split dimer stretched to 8.43 A; ideal-gas <A> 49.8 for 100; node box 552 float steps off after 1e6 calls | z-held atoms block only z; derived slots excluded; whole-molecule centres with consistent image moves; (N+1) ln V; nodes get the absolute box: 400 calls move lipids 2.1 A, protein 0.19 A; <A> 99.8 |
| A75 | `martini_potential.cpp`, `martini_stage_params.cpp`, `main.cpp`, `deriv_engine.h` | branch | table n_grid read wrongly for non-1000 grids; bond spline domain shorter than minimum-image lengths, with an analytic fallback; dead `DihedralSpring`, `MartiniScTablePotential`, an unused angle spline, no-op stage bond/angle calls, a pressure logger, a diagnostic block and a box hook in the core header | wrong energies at 500/1500 points; bond force 0.005 for 25.30 at 10.70 A | n_grid from the table; the bond spline spans every minimum-image length and refits when the box grows; dead code deleted (`deriv_engine.h` back to master's form); regression bitwise |
| A76 | `martini_build_tables.py` side-chain/environment table; `martini_potential.cpp` reader | branch | targets keyed by bead type, so choline NC3 (+1) was neutral; the engine read the table's angle mirrored; the builder put the CB 0.57 A off; the radial grid ended 12 A from the CB though side-chain beads reach 5.2 A further; linear interpolation | ASP-NC3 missed -2.15 E_up at 5 A; a bead 4.5 A from a LEU CB on the side-chain side got the backbone side's value (-1.04 for +151); 1702 of 3192 layers ended with slope > 0.01 (worst -0.414 E_up, 0.341 E_up/A); force jumps up to 78 E_up/A at nodes | targets keyed by (type, topology charge); angle and CB fixed; grid to bead reach + 1.2 nm, ending at value and slope exactly 0; clamped cubic spline in r and theta (37 theta nodes); `parameters/ff_2.1/martini.h5` rebuilt (md5 `cc90a933...`, old kept as `.bak_20261020`); prepared hybrid `.up` files must be re-prepared |
| A77 | `main.cpp` `ReplicaExchange` | yes | the alternating swap set came from a counter that restarts at 0 | a continuation after an odd number of attempts proposed the other set: 1 of 10 frames equal to the continuous run | set chosen by round / replica_interval - 1: 10 of 10 equal; fresh runs unchanged |
| A78 | examples 01-15, `generate_restart_config.py`, `mdtraj_upside.py` | yes | `[ $is_native ]` true for "false"; continuations position-only with the same seed; example 12's slow step at dt 0.009; mts integrator on 08; stale artefacts and docs | frame-0 RMSD 5.8-6.4 A instead of 0; segment step correlation 0.97-1.00; KE/1.5kT 1.04-1.18 | each fixed; continuations carry momentum and round (bitwise the continuous run); up.md parameter-file rows corrected |
| A79 | `martini_potential.cpp` angle term and side-chain table | branch | the MARTINI angle was evaluated analytically (Spline Table Only rule); the side-chain table's energy used the unnormalised CB vector while its gradient assumed a unit vector | node_fd flagged 22 of 26 CB-placement parameter derivatives | angle as a fixed spline of delta^2/2 over delta = cos - cos0 in [-2, 2] padded (5.1e-7 relative to the analytic form, no fallback); CB direction used for theta (0 flagged). Regression bitwise |
| A80 | hybrid HDX analysis: `run_hdx_analysis.sh`, `00` chain, `martini_hdx_project.py`, `martini_hdx_membrane_accessibility.py`, `martini_extract_vtf.py`, `get_info_from_upside_traj.py`, `PDB_to_initial_structure.py` | yes / branch | the shipped driver never applied the lipid term; `stage_7.N` continuations read as replicas; no ladder check in `calc_hdx_ht.py`; contact radius re-measured per replica; one temperature on every group; no terminal corrections inside a multichain sequence; segments dropped from VTFs; sentinels in `percentD_feats.csv`; 0-based residue labels; projection checked only by residue count; static box under NPT; the OPM DUM chain counted; unknown ids plotted as LysM; lipids labelled UNK; a 50 A box fallback | `PS.npy` identical to `PS_protein.npy`; two segments of one run pooled as two states; profiles at T with ESS 1.1; radius 6.6-7.3 A across rungs; 1AFO chain-B rates 0.024x and 63x; 3 of 5 frames | each fixed; one pooled radius recorded beside the ACC arrays; PDB numbering carried; the lipid term changes 7 amides on 1rkl (residue 31: 0.000 to 0.337), 6 on glpG; protein-only outputs identical except the added numbering and censored columns |
| A81 | `00.AnalysisScripts/4.calc_D_uptake.py` `sanitize_peptides`, peptide uptake | branch | experimental peptide Start/End used as 1-based positions in the simulated sequence; a sequence mismatch printed "using bounds" and continued | right only for a structure numbered from 1: numbered from 87, every peptide was dropped, and an offset numbering would sum the wrong residues | Start/End mapped through the PDB numbering (`load_pdb_numbering`); peptides with residues the structure lacks skipped, a sequence mismatch stops; `scratchpad/upside3_audit/hdx_numbering/test_peptides.py` fails before, passes after; identical for structures numbered from 1 |
| A82 | `bead_interaction.h` `quadspline` | branch (follows A33's clamp) | where `angular_coord` holds a cosine at +/-1 the energy no longer depends on it, but the force kept the spline slope there | only off the unit sphere (rounding, or a trainer moving placement vectors): the force's line integral drifted from the energy by 1.3e-2 and 4.8e-2 over 0.8 (1ubq placement params 40, 501; `scratchpad/upside3_audit/fd_candidates_20261021.md`) | angular derivatives zero where the cosine is held: drift 7e-5 and 1.3e-4 (`fd_probe/probe_clamp.py`); 138 of 545 arrays move at rounding, all -ffast-math code placement (a build whose select never fires equals the fix bitwise); baseline `w24_clampB.npz` |
| A83 | `martini_prepare_system_lib.py` `load_sidechain_rotamer_payload`; `martini_potential.cpp` `martini_sc_table_1body` | branch | the hybrid dropped every residue without a side-chain/environment table (ALA, GLY: no MARTINI side chain) from the rotamer rows, which the side-chain node, its placement and the H-bond coverage share | every hybrid lacked those residues' Upside side chains, their pair energies and H-bond coverage (1rkl 3 of 31 residues, glpG 42 of 210; found by the registry, `scratchpad/upside3_audit/registry_20261021.md`) | every residue keeps its rows; table index -1 marks ALA and GLY (no environment term); any other residue without a table raises. `scratchpad/upside3_audit/alagly/repro.py`: rows 117 and 747 before, 120 and 789 after, identical to upside_config; table-less rows exactly 0; node_fd clean; ff2.1 regression bitwise |
| A84 | `coord_const.cpp` `Const1D`/`Const2D`, `upside_nodes.py` `Const1DCoord`/`Const2DCoord`; `sidechain_radial.cpp` `WallSpring` | yes | dead: no caller writes `Const1D` or `Const2D`, and `Const1D` without `initialized_by_coord` could not be built (`value` sized rank 1, checked rank 2); an unregistered second `struct WallSpring` beside `spring.cpp`'s broke the one-definition rule | none at run time | deleted; regression bitwise; fixtures regenerated from the current writers (`fixtures_w4`), baseline `w26.npz` equal to w25 |
| A85 | `main.cpp` `compute_logged_subset_kinetic_energy` | branch | the logged `kinetic`, `protein_kinetic` and `lipid_kinetic` gave every non-fixed atom 3 degrees of freedom, also PO4 atoms whose z is held for the first 150 production rounds | `avg_kinetic_energy/1.5kT`, the run health check, read 2.3% low during the hold; dynamics unaffected | degrees of freedom counted (2 for a z-held atom): logged over the DOF-correct value 0.9773 during the hold and 1.0000 after (old), 1.000000 in both (new); ff2.1 regression bitwise |
| A86 | `engine_c_library.cpp` `free_deriv_engine` | branch (A56) | dereferenced the NULL engine that `construct_deriv_engine` returns for a refused config, when `Upside.__del__` freed the half-built object | a config the engine refuses killed the Python process with SIGSEGV after its RuntimeError was caught (found by the option fixtures) | NULL is a no-op, as `delete` is: the reproducer (`option_fixtures/repro_refused_config_segfault.py`) exits 0 after the except block (-11 before); regression w28 to w29 bitwise |
| A87 | `environment.cpp` `NonlinearCoupling`, `SigmoidCoupling`; `upside_config.py --environment-weights-number` | yes | mode 1 returned 0 for all 400 weight derivatives while the energy uses the stored library weights, which the help claimed were all 1.0 | a silent zero for parameters the energy depends on (a weight slot's central difference 0.083 against 0); no trainer uses mode 1 (ConDiv refuses anything but 20) | mode 1 deleted: the engine accepts 20 or 400 and refuses 1, the writer offers 20 and 400; the 28 other option fixtures bitwise |
| A88 | `upside_config.py` `write_rotamer_placement` (`--fix-rotamer`) | yes | unpacked six fields per `.chi` row, but `PDB_to_initial_structure` leaves the chain column empty for a PDB without chain ids | `--fix-rotamer` raised ValueError on any chainless structure (chig) | the reader takes residue, restype, chi1, chi2 as fields 0, 1, -2, -1 (the writer is master's format); chig fixes all 10 residues; a chained fixture (1ubq_rotamer) bitwise equal |
| A89 | `upside_config.py --dynamic-rotamer-placement`; `placement.cpp` `placement_point_only`, `placement_point_vector_only` | yes | the option reads `rotamer_center`, which none of the 13 side-chain libraries holds; the two Rama-dependent point placements had no other writer | the option always ended in a NoSuchNodeError; two registered nodes nothing could build | option, its writer path and both registrations deleted; regression bitwise (w6 to w30, 0 of 4073) |
| A90 | `upside_config.py --remove-pivot`, `--rama-sheet-library`, `--trans-cis` table; `run_upside.py` `cavity_radius` | yes | parsed and never read (the trans-cis `init_state` column dropped by its reader) | dead inputs that looked like controls | deleted; the trans-cis table is `residue energy k width` (example 09's table and plotter follow); every fixture's datasets bitwise equal |
| A91 | `upside_config.py` input checks; `apply_spatial_transformation` | yes | checks that could not fire (`is None` on a store_true flag and on a '' default); writers reading nodes other options must supply; `--surface-method 2` written though the engine lacks it; `parser.error` called out of scope | invalid configs written (a channel potential with no `surface` node), or failures deep in a writer | validated preconditions: channel needs `--surface`; bb environment needs `--environment-potential` and `--hbond-energy`; both membrane potentials need `--environment-potential` type 1 and `--hbond-energy`; `--surface-method` 0 or 1; the spatial table raises ValueError; fixtures bitwise |
| A92 | `example/07.MoreRestraints/0.run.py` | yes | the structure was converted with `--disable-recentering`, leaving 1UBQ 37-47 A from its tables' reference points | the default wall pulled CA0 and CA5 with 127-148 E_up/A and stretched a backbone bond to 3.36 A within 20 time units | the conversion recentres (the Upside run keeps `--disable-recentering`, A66): start forces 35 and 11 E_up/A, max backbone bond 1.99 A |
| A93 | `advanced_config.py` `write_group_const_spring`, `write_group_nail_spring` | yes | atom indices parsed with `float()` (the A24 class) | `id` stored as float64, right only while HDF5 converts integral values | `int()`, as `write_group_dist_spring`; ids int64 with equal values, node_fd and regression unchanged |
| A94 | `upside_config.py` `write_rotamer_backbone` | yes | the pre-proline backbone class went to any residue followed by PRO, also across a chain break | a chain's last residue took the pre-proline coverage rows when the next chain starts with proline (1dfn with residue 30 PRO: 24 coverage outputs moved, up to 1.77; potential -8.130 against -7.955; found by `tests/index_oracle.py`) | the class from one per-residue list respecting `chain_first_residue`, the dead `pre_PRO` toggle gone; the oracle leaves only D9's exclusions; 28 rebuilt fixtures equal fixtures_w6 |
| A95 | `rotamer.cpp` stop test, `max_deviation`, constructor | yes (stop test); A20 (NaN) | the stop test compared the damped step, (1-damping) times the residual, with `tol`, enforcing tol/(1-damping) (at 1.0 the solve stopped at once); `std::max` dropped a NaN deviation, so A20's count never saw one; any damping and an `iteration_chunk_size` of 0 accepted | at the shipped 0.4 beliefs stopped up to 1.6x tol from converged (free energy up to 2.4e-3 E_up, forces 1.1e-2 of RMS); at 0.99 55-97x; at 1.0 1ubq's free energy 19.4 E_up high with `n_bad_solve` 0; a NaN solve counted 0 | residual = damped step / (1-damping), NaN kept, damping in [0,1), chunk size >= 1: errors at most 1.2x tol for damping 0-0.97, NaN solves counted, converged answer unchanged; node_fd chig 0 bad of 42; only the 315 rotamer-fed arrays move |
| A96 | `environment.cpp` `EnvironmentCoverage`, `BackboneSigmoidCoupling`; `sigmoid_energy.cpp`; `upside_config.py` `read_pro_state_file` | yes | `aa_types` not checked against `num_aa_types`; `bb_sigmoid_coupling` coverage sizes (n_bb x 20, n_bb) unchecked; `SigmoidEnergy` windows closed at both ends; trans-cis widths above 45 degrees accepted (the four 90-degree windows no longer meet) | a bad type moved coverage into the next CB's channel (ASan BUS error); mismatched coverages overflowed (energy -21.7 or -22.5 for -16.5); E doubled at omega exactly +-90; width 60 a discontinuous wall (dE 0.049 against force work 0.010) | validated preconditions, windows half-open [min, max): the malformed configs and width 60 refused; regression bitwise |
| A97 | `sidechain_radial.cpp` `radial`, `hbond_sc_radial`, `symm_fast_contact`, `nosymm_fast_contact`; `environment.cpp` `linear_coupling_uniform`; `upside_config.py --rotamer-exclude-residues`, `--env-exclude-residues` | yes | five registered nodes with no writer (one with an `is_compatible` that accepted asymmetric tables), and two options the writer refused at parse time while `write_rotamer_placement` and `write_environment` kept their branches and two unused parameters | dead code | deleted with the `run_upside` keywords; regression bitwise; 28 rebuilt fixtures equal fixtures_w6 |
| A98 | `interaction_graph.h` `get_param_deriv`, `propagate_derivatives` | yes | the parameter derivative re-ran the whole backward pass, adding the pair forces into the input nodes' sens again | `protein_hbond` sens 2.000x after one `get_param_deriv`, 3.000x after two (latent: nothing reads sens after it) | `get_param_deriv` accumulates from the stored edge sensitivities, `propagate_derivatives` does positions only: ratios 1.000000; only the `param_deriv` arrays of `rotamer`, `hbond_coverage`, `hbond_coverage_hydrophobe` move (median 1 ulp, max 1.9e-6), bitwise with FMA contraction off |
| A99 | `bead_interaction.h` `quadspline_param_deriv` | yes | the four coefficient blocks were written with `=`, so at a cosine of +1 the angular block's last bin (weight 0) overwrote the wide radial spline's coefficient 0 | the training gradient of radial coefficient 0 of the pair, coverage and hydrophobe tables read 0 for pairs closer than 1.4 A with a bead pointing at the other (exact 3.1e-2 at 1.0 A, 1/6 at 0.5 A); energies and forces unaffected | every block adds into the zeroed derivative: the harness gives 0 wrong of 54 slots in the targeted cases and of 20,000 random pairs (reproduced with the fix reversed) |
| A100 | `interaction_graph.h` constructor, `set_param` | yes | a refused `set_param` had already copied its parameters in; index and type values were not range-checked; a symmetric graph was not required to have a square table | after a refused call later energies used the refused table (8.669 for 5.942); `type1` 11 of 8, -1, or `index1` 22 of 17 were accepted and read outside the arrays | compatibility is checked on the candidate before it is copied; indices, types and the square table are validated preconditions: the refused call changes nothing, the four bad configs are refused; regression unchanged by this part |
| A101 | `hbond.cpp` `ProteinHBondInteraction::compute_edge` | yes | the 90-degree angular cutoff was tested per SIMD group of four edges, so an out-of-cone pair was evaluated whenever another pair of its group passed | the H-bond energy depended on edge order (the same config with rows reversed: 6.2e-3 E_up on chig_hbond, 4.6e-2 on 1ubq with `--loose-hbond-criteria`; below float32 at the default geometry) | every edge evaluated (the cutoff was only a shortcut): row order gives 0 difference, the engine equals the float64 no-cutoff oracle; node_fd chig_hbond 0 bad of 42; regression: forces up to 5.3e-5 E_up/A and the seeded runs move |
| A102 | `martini_prepare_system.py` `assert_environment_solvation`, `run_stage70_burnin`, `--eq-66-nsteps`; `example/16.MARTINI/run_sim_hybrid.sh`, `run.py`, readme | branch | the stage-7.0 gate failed any belt BB site without an environment bead within 6 A on one frame, counting the inside of a helix bundle; stage 6.6 ran 500 NPT steps | it failed every healthy state (38% of 1rkl's sites beyond 6 A) and could not tell glpG's real lipid-face gap from its interior; box and lipid face entered production unrelaxed (findings 8b) | the gate averages the lipid-facing sites (fewer than 60 backbone atoms within 10 A), failing above 10% beyond 8 A and warning at that fraction for acyl tails; 6.6 runs 60,000 steps (user 10-21): rejected glpG 17.9% fails, relaxed glpG 3.6-5.5% and 1rkl 0.4% pass. Its window, the burn-in's last quarter, was replaced by A126 |
| A103 | `spline.cpp` `solve_clamped_1d_spline`, `spline.h` `LayeredClampedSpline1D::fit_spline`, `martini_potential.cpp` pair table; `martini_prepare_system_lib.py` `inject_particles_table` | branch | the dry-MARTINI pair spline was clamped to zero slope at both ends, though its grid ends at the 1.2 nm cutoff, where the potential-shifted LJ keeps a force | the pair force fell to zero over the last spacings before 12 A, up to 2.7e-2 E_up/A off the published form (the stage-6.0 soft core likewise) | the clamped spline takes the analytic dV/dr at r_min and, from inside, at r_max (`combined_energy_end_slopes`; the soft core's from its last smooth piece). 13,041 pairs over all 81 types (`d5_pair_ends/`): last 10 spacings 1.5e-7 E_up/A (soft core 3.6e-5), the rest unchanged; ff2.1 regression 0 of 4073; hybrid (`fixtures_hybrid_w2/`): only forces (up to 3e-2 E_up/A), one pair energy (2e-3 E_up) and seeded runs move; node_fd flags are float32 noise. Hybrid `.up` files prepared before it are refused |
| A104 | `compute_rotamer_centers.cpp`, `compute_rotamer_pos.cpp`, `generate_from_rotamer.{h,cpp}` and their `CMakeLists_{M1,Other}.txt` targets | yes | rotamer-library generators that no longer build (clang rejects a `constexpr` `atanf`; Eigen's OpenMP symbols unlinked), with no caller, inputs at `/home/jumper/...`, 275 of 623 lines commented out (placement pass 2 D4) | dead code | deleted (plan Phase 24 decision 2); master keeps them as the libraries' provenance (up.md 2.1); fresh configure and build, regression 0 of 4073 |
| A105 | `membrane_potential.cpp` `MembraneLateralPotential`; `upside_config.py --membrane-lateral-potential`; `run_upside.py`; `main.cpp`; example 08 `2.lateral_pressure.run.py` | yes | an xy force with no energy on surface CBs, P(z) times each side's area share at each CB's own height, so the sides do not cancel (a pressure on a closed outline gives zero net force per slice) | net force on the protein up to 0.096 E_up/A, z torque up to 1.27 E_up on 1qhj, a closed rigid loop doing up to 0.11 E_up (`lateral/net_force.log`); the logged potential, pivot MC and replica exchange never saw it | deleted (decision 1); fixture `1rkl_lateral` became `1rkl_surface`; regression 0 of 3912; the old config refused |
| A106 | `PDB_to_initial_structure.py` chain breaks; `upside_config.py --chain-break-from-file` | yes | breaks were recorded only with `--record-chain-breaks`, and upside_config fell back to one chain when the file was missing | 1dfn converted without the flag: exit 0, n_chains 1, a 17.74 A C-N bond at equilibrium 1.300 A (6486 E_up); a mistyped path passed silently (writers pass 2 item 6) | breaks always found and refused without the flag; with it the file is written for one chain too; a missing file refused (decision 3): `chain_breaks/repro_{before,after}.log`; fixtures equal. A reused single-chain input directory must be re-converted |
| A107 | `PDB_to_initial_structure.py` `read_residues` | yes | insertion codes were dropped from the `.chi` residue numbers | residues 5 and 5A both written 5, where every reader of the numbering (HDX projection, labels) is integer-only (writers pass 2 item 13) | a residue with an insertion code is refused with a request to renumber (decision 4; none of the project's 46 PDBs has one); `insertion_codes/repro_{before,after}.log`; 28 fixtures rebuilt equal |
| A108 | `eig.cpp` `house()` | yes | sigma2 floored at 1e-20, so an exactly zero column got beta = 1 with v = (1, 0, ...): a projection, not a reflection, deleting one axis of the 4x4 problem (placement pass 2 D1) | a residue whose N/CA/C are its reference triangle turned exactly 180 degrees about y got a frame wrong by pi, its CB 3.87 A off; needs exact symmetry (1e-6 A of noise leaves 8e-8 rad) | no floor; a column with no tail is already reduced, beta = 0 (Golub and van Loan 5.1.1); dead `house_accurate`/`givens_accurate` deleted. `eig_fix/`: frame error 3.2e-18 rad; with FMA contraction off only 11 of 3912 arrays move (seeded runs through the changed lanes) |
| A109 | `placement.cpp` `layer_index`, `rama_residue`, `affine_residue`; `eig.cpp` `affine_alignment` `atoms` | yes | the placement and alignment index datasets were not range-checked (placement pass 2 D2) | a bad index was accepted and read past its table or into padding: chig energy 5.9670 for 5.9421, a `placement_scalar` 1-body 0.461 for 9.983, an alignment energy 18.87, one SIGSEGV | validated preconditions where each is read: the six bad configs refused (`bounds_fix/results.txt`); regression 0 of 3912 against A108's build |
| A110 | `eig.cpp` `AffineAlignment` reader | yes | the derivative takes dR/dx_a as ref_geom[a], true only when each residue's ref_geom is centred, which was neither checked nor documented (placement pass 2 D3) | a shifted ref_geom gives the right energy (1.4e-5 off) and wrong forces (up to 8.1 against a largest 11.0; analytic/FD 6.96/0.42); latent, every writer centres it | a ref_geom whose centroid is beyond 1e-4 A is refused (writers below 1e-7), up.md documents it: `ref_fix/repro_after.log` (shifted refused; as written unchanged); regression 0 of 3912 against A109's build |
| A111 | `advanced_config.py` nail and group nail (`add_to_const3D2`), `write_group_node`, the group springs; `upside_nodes.py` `Const3DCoord2` | yes (group nail branch-only) | nail anchors written with `initialized_by_coord 1`, taken from `/input/pos` at first evaluation, which every continuation rewrites; group numbers and the `--group` table's "residues" header unchecked (writers pass 2 items 1, 4) | a continuation re-anchored the nail to its restart frame (anchor 1.15 A, nail energy 1.81 to 0.00, positions 1.63 A and potential 7.6 off the continuous run; group nail 2.23 A); a group number past the table read unrelated memory (ASan); a headerless table lost its first group | anchors stored from the config's positions (`Const3DCoord2`; `add_to_const3D2` deleted); groups checked against `n_group`, the header required. `nails_groups/`: continuations bitwise the continuous run, four bad tables refused; only the anchor datasets of chig_groups, chig_restraints_1d2d, chig_restraints_xyz change |
| A112 | `upside_config.py` `write_surface_coord` | yes | `--surface-included-residues` was written unchecked (writers pass 2 item 5) | 4-40 on 31-residue 1rkl accepted at exit 0; ASan heap-buffer-overflow in `Surface::compute_value` (`surface.cpp:193`) | residues past the chain refused; `writer_fixes/surface_{before,after}.log` (4-26 clean under ASan, 4-40 refused) |
| A113 | `upside_config.py` `--trans-cis` | yes | `write_rama_map_pot2` takes neither `--secstr-bias` nor `--rama-param-deriv`, which were dropped at exit 0 (item 7) | the options changed 1 and added 14 datasets without `--trans-cis`, none with it | the combination refused (`writer_fixes/transcis_*`) |
| A114 | `advanced_config.py` wall readers, `write_cavity_radial` | yes | wall_type parsed as any float and stored float, where the engine reads an int and tests only 0 and 1 (item 8) | 0.5 became a lower wall; 2 or -1 an undocumented two-sided spring (18 and 8 E_up where 0 and 1 give 18/0 and 0/8) | 0 or 1 only, stored int; `writer_fixes/wall_after.log`; five fixtures change only wall_type's dtype, regression 0 of 778 |
| A115 | `advanced_config.py` main | yes | `--tension` with `--ask-before-using-AFM`, `--make-unbound` without a cavity radius, and `--heuristic-cavity-radius`, `--cavity-radius-from-config`, `--make-unbound` on one chain printed a message and wrote nothing at exit 0; an unreachable overwrite branch (item 9) | requested restraints silently absent | each an error; the dead branch deleted (`writer_fixes/adv_dropped_options_*`) |
| A116 | `upside_config.py` `random_initial_config` | yes | the bond-length table was shifted by one bond: atom i sits at b[i] from atom i-1, and the table held the bond starting at each atom (item 10) | the random start chain had N-CA 1.526, CA-C 1.300, C-N 1.453: 18.6 E_up of bond strain at chig's start (example 11) | lengths of the bonds ending at N, CA, C (1.300, 1.453, 1.526): strain 0.000 (`writer_fixes/random_chain_bonds_*`) |
| A117 | `PDB_to_initial_structure.py` `read_residues` | yes | selenomethionine's SE does not match the CD pattern (item 11) | chi2 nan for every MSE, so `--fix-rotamer` left them free silently | SE read as SD for MSE: chi2 -172.151 as MET (`writer_fixes/pdb_cases_*`, other cases unchanged) |
| A118 | `upside_config.py` `rotate_to_axis` (`--spatial-transform-from-table`) | yes | exactly opposite vectors returned the identity; a dot product rounded past +-1 gave NaN (item 12) | the protein pointed the wrong way; NaN coordinates in 11% of pairs within 1e-13 of (anti)parallel | a half turn for opposite vectors, the dot clipped: NaN 0 of 40,000, opposite maps a to -a (`writer_fixes/rotate_*`); fixtures unchanged. Pairs within 1e-13 of anti-parallel keep up to 3.5e-2 direction error from the cross product's axis, which real coordinates cannot reach |
| A119 | `upside_config.py` `write_rotamer_placement` `--fix-rotamer` | yes | the chi2 distance had period 360, where the library holds ASP, PHE and TYR on a 180-degree window (spans 137, 159, 157 against about 340 for the rest) because a half turn leaves those side chains unchanged (item 3) | the fixed state followed the PDB's naming of OD1/OD2, CD1/CD2, CE1/CE2: relabelling changed 61 of 77 states, 18 of 76 not the nearest modulo 180 | period pi for the three: 1 of 77 changes (1qhj TYR 21 at the 63.9 boundary, ring not planar), 0 of 76 differ (`writer_fixes/{relabel,symmetry}_*`); fixtures unchanged |
| A120 | `generate_restart_config.py` `main` | yes | three or more arguments did nothing at exit 0 (writers pass 2 item 15) | a mistyped call wrote no restart config and reported success | other argument counts refused with the usage; `writer_fixes/restart_args_{before,after}.log`; a plain restart still continues bitwise (`restart_check_after.log`) |
| A121 | `main.cpp` `CurvatureChange`; `generate_restart_config.py`; example 08 `4.curvature_dynamics1.run.py`, `5.curvature_dynamics2.run.py` | yes | the changer moved `Const3D_curvature_center` in memory only, logged it at LOG_DETAILED, and every continuation restarted from the config's centre (item 2) | a `--curvature-changer-interval` continuation left the continuous run (1rkl_curv: positions 2.5 A, potential 16 E_up off); the sampled curvature went unrecorded at the default log level | `/output/curvature_center` logged every frame and copied into the centre by each continuation: bitwise the continuous run with and without the changer (`curv_fix/restart_{after,plain_after}.log`); regression 0 of 3912 |
| A122 | `hbond.cpp` `HBondEnergy` reader | yes | `donor_resid`, `acceptor_resid` and `rama_resid` unchecked, the donor/acceptor count only an `assert` (H-bond pass 2) | a donor at residue 0 wrote `hb_number2[-1]`, an acceptor at the last residue `hb_number2[n_res-1]` (ASan, chig); a one-row `rama_resid` read padding at exit 0; no writer produces them | validated preconditions: donors 1..n_res-1, acceptors 0..n_res-2, `rama_resid` the rama_coord rows in order, donors plus acceptors the `protein_hbond` rows (`hbond_inputs/repro_asan_before.log`, `results.txt`: the three refused); regression and hybrid bitwise |
| A123 | `rotamer.cpp` `residue_free_energies`, `residue_entropy`, the 1xN edge conversion | yes | a pair between a 1-rotamer and a multi-rotamer residue is folded into the partner's node energy, and the per-residue logs (`rotamer_free_energy`, `rotamer_potential`, `get_value_by_name`) gave the whole pair to the partner, where other pairs split half and half (BP pass 2 D4); folded edges were still converted to unread probabilities | chig's 1-rotamer residues 0, 6, 9 missed -0.027, -0.091, +0.039 E_up; sums right | the logs give half of each moved edge's expected pair energy (under the partner's marginal) back to the 1-rotamer residue (`bp_d4/split_{before,after}.log`: 0.0872 to 0.0601, -0.0858 to -0.1772, 0.5788 to 0.6173, the sum still F); regression and hybrid bitwise |
| A124 | dead code: `rotamer.cpp` (an unused energy cap and its parameters, an uncalled `ensure_fresh_energy` and its flag, a doubled `edges11.reset()`, commented-out blocks), `bead_interaction.h` (`PosDistSplineInteraction`, the never-defined `PARAM_7A/10A_CUTOFF` variants), `interaction_graph.h` `message()`, `spline.{h,cpp}` (`deBoor2d_coeff_deriv`, `evaluate_1d_spline`, `evaluate_1d_spline_first_derivative`, `evaluate_2d_spline`), `sigmoid_energy.cpp` (`use_cutoff`, written 1 by every writer, and its unreachable branches), `hbond.cpp` `ProteinHBondInteraction::param_deriv`; `run_upside.py` (`compile`, `continue_sim`, `read_output`, `compute_com_dist`, `status`, `read_hb`, `read_constant_hb`, `compute_topology`, `display_structure`, their helpers and imports), `upside_nodes.check_node_exist`, `upside_config` (`import pdb`, `write_infer_H_O`'s multichain check, unreachable since A17, commented-out code), `PDB_to_initial_structure` commented-out code | yes | code nothing reaches; `ProteinHBondInteraction::param_deriv` returned -1 for every parameter; `SigmoidEnergy`'s `dim1` held only by an `assert` | none at run time (`run_upside.continue_sim` was a position-only continuation, the A42/A78 class) | deleted, each shown unused by grep (`py`, `example`, `tests`, `training`); the H-bond `param_deriv` raises "not implemented", `dim1` a validated precondition; regression and hybrid bitwise (`dead_code/results.txt`) |
| A125 | `coord_basic.cpp` Distance3D/2D/1D `id`; `spring.cpp` WallSpring `id`, `wall_type`, `dim1`; `surface.cpp` `included_list` | yes | these indices and the wall type were stored unchecked (`dim1` only an `assert`), behind the writer checks of A111, A112 and A114 (writers pass 2 suggestions) | broken configs ran at exit 0 reading padding or other rows: a Distance id at n_atom, a Const3D row past its values, a WallSpring id past its distances, wall_type 2 (a two-sided spring), surface residue 40 of 31 | validated preconditions where each is read: all six refused (`reader_checks/repro_before.log`, `results.txt`); regression and hybrid bitwise |
| A126 | `martini_prepare_system.py` `assert_environment_solvation`, `run_stage70_burnin`, `run_stage70_handoff`, `--prod-70-nsteps`; `example/16.MARTINI/run.py`, readme | branch | A102's gate read the last quarter of the rigid burn-in, where the belt ends wander 7-16% bare under the frozen protein whatever production then does (findings 8b) | it failed one healthy glpG prep in three (s13 16.1%, which production then closed to 4.5%) | the gate measures production's steps 15,000-20,000 (user 10-21: option (a)); stage 7.0's production runs at least 20,000 steps when the belt is defined (a production-free prep says it is ungated); `run.py` starts its replicas from the gated run's last frame. `glpg_gate_20261021/gate_a/test.log`: glpG 1.3% and 2.9%, 1rkl 0.6% pass, a void cut into the lipid face 93.3% fails; a fresh 10,000-step OPM workflow refused, 20,000 and a 5,000-step continuation accepted |
| A127 | `main.cpp` MovingConst3D check; `coord_const.cpp` MovingConst3D clock | yes (g-JF branch-only) | the tip's clock advances by time_step per force evaluation and A54 required only time_step == dt: right for Verlet and mv's fast level, but g-JF evaluates 1 + inner_steps times a step and mv's slow level once a round plus per MC move (coordinate pass 2 C1) | time_estimate / time: hybrid at inner_steps 4 5.0000 (the tip 5 times too fast), a slow-level tip 0.3333, with pivot moves 0.89 drifting; no writer combines them | MovingConst3D refused under `/input/brownian` and, under mv, off integrator level 0 (`coord_fix/moving_clock*.py`; Verlet and mv fast at ratio 1.0000); regression and hybrid bitwise |
| A128 | coordinate readers: `coord_basic.cpp` (Distance2D/1D dims, Angle, Dihedral, GroupCenter), `coord_const.cpp` (Const3D, MovingConst3D), `coord_map.cpp` (rama_coord, cat_pos, SigmoidCoord), `coord_operator.cpp` (Add, Mean, Multiply) | yes | indices and dimensions held only by compiled-out `assert`s against padded sizes, so a row just past the end read a zero padding row ASan does not see; Const3D's value check compared the dataset with itself; an empty group (a reversed range "9-5") accepted (coordinate pass 2 C3, C4) | 22 of 26 broken configs ran at exit 0 on wrong values (distances 0.0 from the padding lane, phantom atoms at the origin, lambda 1.0 from padding), the rest overflowed or crashed; the empty group sat at the origin, its force lost | validated preconditions where each value is read (widths by `check_elem_width`, Const3D value rows equal to id rows, GroupCenter borders in range, no empty group): 26 of 26 and the empty group refused (`coord_fix/`); jac_fd 86 of 86; regression and hybrid bitwise |
| A129 | `advanced_config.py` `parse_residue_list`, `parse_segments`; `upside_config.py` `parse_segments` | yes | a reversed range ("9-5") parsed to nothing (coordinate pass 2 C4, the writer side) | a group, wall, spring or nail table, `--hbond-exclude-residues`, `--membrane-exclude-residues` or `--surface-included-residues` with one silently lost those residues; a `--group` row of it gave an empty group at the origin whose force reached no atom | a reversed range refused by all three parsers; valid ranges parse as before; fixtures_w15 all 28 equal; the engine refuses an empty group as well (A128) |
| A130 | `coord_const.cpp` Const3D and MovingConst3D; `upside_nodes.py` `Const3DCoord1`, `MovingConst3DCoord1/2`; `advanced_config.py` `write_pulling`, `--AFM-time-initial`; `generate_restart_config.py`; example 06; `tests/fixtures.py`, `tests/index_oracle.py`, up.md | yes | an AFM tip without given positions was written with `initialized_by_coord 1`, its start taken from `/input/pos` at first evaluation, which a continuation rewrites (A111's mechanism for nails); Const3D kept the branch with no writer left (coordinate pass 2 C2, C5) | a continuation built as the help described put the tip at CA(T) + v*T, 2.968 A from the continuous run (the stationary tip 1.470 A); a fresh auto-tip run with a time offset started v*t0 off its CA | the nodes carry their values (Const3D `value`, MovingConst3D `velocities` and `start_pos`); a tip without a given start is written on its CA at time 0 (an offset refused; continuations by `generate_restart_config.py`); `initialized_by_coord` and `MovingConst3DCoord2` deleted. `c2_fix/`: the auto-tip continuation bitwise over 3 + 3 time units (longer: A222), the old auto-tip config refused; only `chig_afm_auto` moves (its offset 5 to 0); hybrid bitwise; jac_fd 86 of 86 |
| A131 | `coord_basic.cpp` VectorAngle; `coord_map.cpp` paste_pos, select_pos, slice, scaled_sum, backbone_featurizer, conv1d (and its `relu`); `coord_operator.cpp` Mean; the unused Eigen includes | yes | registered coordinate nodes with no writer in the repo, master, examples, tests or training (VectorAngle's writer went in A55), several with unchecked inputs of their own (coordinate pass 2 C5) | dead code | deleted, each shown writerless by grep for its registered name; regression 0 of 4070 on fixtures_w16, hybrid 0 of 196 (`c5_fix/results.txt`) |
| A132 | `spring.cpp` Spring and WallSpring, `sigmoid_energy.cpp` `dim1`; Spring `id` and `box_len` | yes (A124/A125's checks branch-only) | Spring's `dim1` and `pbc` box only `assert`s, its ids unchecked; A124's SigmoidEnergy and A125's WallSpring `dim1` checks compared against `row_width`, the padded width (4 for 3-wide positions) | Spring_bond with dim1 1 on the 1-wide Distance3D acted on the next bond's length; WallSpring and SigmoidEnergy with dim1 3 on pos read the padding lane; all at exit 0 | `dim1` against the input's true width (`elem_width`), Spring's ids and a positive pbc box validated preconditions: all three refused (`c5_fix/repro.py`); regression and hybrid bitwise |
| A133 | `main.cpp` `--restart-using-momentum` | yes | `/input/mom` was read with the dataset's own shape into a momentum array of exactly n_atom rows (engine pass 2 D2) | 40 rows on a 30-atom config: ASan heap-buffer-overflow write; 20 rows: exit 0 with atoms 20-29 on allocator contents | its shape (n_atom, 3, 1) checked before each read: both refused, the exact shape runs (`eng_fix/results.txt`) |
| A134 | `main.cpp` `--set-param` | yes (A57 incomplete) | A57 checked the library's `set_param` only; the engine's option passed the file's vector unchecked (D3) | 10 values for `cb_membrane_potential`'s 720: ASan container-overflow read | the parameter count checked: "expected 720 parameters but got 10" |
| A135 | `main.cpp` `--integrator`, `--inner-step` | yes | neither validated; anything but `mv` ran Verlet (D4) | `--integrator MV` ran Verlet at exit 0; `--inner-step 0` moved nothing and a `--duration` run never ended (n_round inf); -2 ran negative time | `v` or `mv` only, an inner step of at least 1: the four cases refused |
| A136 | `deriv_engine.h` `Node::integrator_level`; `main.cpp` final KE report | yes | the pos node's level was never initialised, then branched on (D5, no effect: its compute and sens are no-ops); the report averaged frames nf > n/2, one short of the second half (D7) | undefined behaviour without a visible effect; a 1- or 2-frame run printed nan on the health-check line | the level initialised to the reader's default 2; the report takes nf >= n/2 (chig, 2 frames: 0.895). A133-A136 together: regression 0 of 4070 on fixtures_w16, hybrid 0 of 196 |
| A137 | up.md 1.3; `deriv_engine.h` `integration_cycle`; `vector_math.h` `sigmoid` comments | yes | `@integrator_level` was documented on every node with "1 = evaluated at backbone level"; `integration_cycle` as "3 time steps"; `sigmoid` as 1/(1+exp(x)) (D6) | misleading documentation | corrected: optional on potentials, read only under `--integrator mv`, 1 the slow level; the cycle per integrator; 1/(1+exp(-x)) |
| A138 | `vector_math.h` `mag` (`a_sqrt`) | yes | `mag(v)` was a*rsqrt(a) with a = \|v\|^2: NaN (0 times inf) for a zero vector, 0 only through fast-math's rewrite to sqrt(a), on which Distance3D/2D rely at dist 0 (engine pass 2 D1) | since A111 every nail starts at distance 0 from its anchor: without that rewrite (strict IEEE, debug, a port) NaN potential and force at frame 0 (`d1_fix/`) | `mag` takes `sqrtf(mag2(a))`, `a_sqrt` deleted: without -ffast-math 14.7042, 15.4340; the normal build bitwise on fixtures_w16; hybrid 44 of 196 arrays move at rounding (about 1e-7 of scale), all downstream of `martini_sc_table_1body`'s `mag` |
| A139 | `martini_potential.cpp` `martini_potential` `pairs`, `dist_spring` and `angle_spring` `id`, `restraint_position` `restraint_indices`; `martini_brownian.cpp` `atom_index` | no | indices unchecked: the cell list finds a pair only as (lower, higher) by binary search, so a pair written j > i or out of order was never evaluated, and a first index >= n_atom wrote past `adj_row_start`; an id of n_atom read the zero padding row; Brownian atoms unchecked for repeats or coverage (dry-MARTINI pass 2 D1) | exit 0 on wrong values: a reversed pair -1.4174 against -2.0677, `dist_spring` id n_atom 713.55 (0.103), `angle_spring` 6.652 (2.144), `restraint_position` 28.126 (0.686); atoms listed twice diffused 4.4 times too fast, atoms left out never moved | pairs in range, i < j, sorted, distinct; ids in range; Brownian atoms distinct and, with the fixed atoms, covering every atom: the nine cases refused (`mfix/repro_*.log`); every workflow stage file 6.0-7.1 meets them (3681 listed and 62 fixed of 3743); fixtures and hybrid bitwise |
| A140 | `martini_hybrid.cpp` `initialize_sc_env_po4_z_reference`, `martini_internal.h` | no | a PO4 z reference filled on every pair evaluation through a `const_cast` and read by nothing; the hold works by z-fixed atoms (D2) | none | deleted; bitwise |
| A141 | `martini_potential.cpp` potentials of `martini_potential`, `dist_spring`, `angle_spring`, `restraint_position` | no | summed in a float running sum (about 6.5e5 pair energies, total about -2.3e4 E_up, in a membrane box) (D3) | off a float64 sum by 0.07-0.40 E_up, nearby energy differences by sd 0.095 (max 0.38 E_up); the prep minimiser rejected descent steps on that noise and stopped at max\|grad\| 1.3-2.7 E_up/A (6.0: 8 iterations to -24495.4, against 5000 to -24711.7 with double sums; 6.1: 260 to -24798.0 against 5000 to -24983.3); the MC barostat compares the same totals | each node sums in double: within 0.0012 E_up of float64, dU error sd 0.0006 (`pass2/martini/float_sum/`); hybrid: only potential values move (node energies up to 0.40 E_up, logged totals up to 0.54), forces and positions bitwise (`mfix/results.txt`). Deliberate departure: the prep minimiser now runs to its iteration budget |
| A142 | `main.cpp` `compute_logged_subset_kinetic_energy`; `martini_fix_rigid::get_rigid_groups` | no | the logged kinetic energy counted 3 degrees of freedom per atom of a rigid group, which has 6 (D4) | in 6.0-6.6 and the 7.0 burn-in, where the protein's N/CA/C are one rigid body (279 counted for 1rkl's 6), `protein_kinetic`/1.5kT read 0.04 and the health line read low by the protein's share (0.981 at 6.6) | each group counted with 6: from 1rkl's 6.6 input `protein_kinetic`/1.5kT 1.9913, equal to the rigid body's own from the momenta (a transient, pass 2 U6); logging only, fixtures and hybrid bitwise |
| A143 | `main.cpp` restart branch | no (A42's continuation) | an annealed continuation started at `--temperature` and reached the schedule only at its first thermostat round (driver pass 2 item 2) | until then the logged temperature, MC, exchange and curvature acceptance and the g-JF bath ran at T0: chig logged 0.85 against 0.6409; a hybrid's bath ran at 0.8647 instead of 0.4686 for 60 rounds, 0.42 A and 37 E_up off after 100 | the restart sets the temperature the last thermostat round before it set: every annealing case bitwise with the continuous run (`dfix/logs/restart_anneal*.log`) |
| A144 | `main.cpp` stop test | yes (the refusal to continue is A42's) | each system tested SIGTERM and `--time-limit` on its own rounds (nr%8 == ns%8) inside its block, and `passed_time_lim` was written by several threads (item 4) | a multi-system run stopped at different rounds per system (50680, 53017, 61634, 52267) and then could not be continued ("all systems must restart from the same round") | blocks end every 8 rounds or at an exchange; signal and time limit tested between blocks, after the swap: all systems stop together and continue (`dfix/logs/*_stagger.log`); 4 systems without exchange cost 6% more |
| A145 | `main.cpp` round order, final capture | yes (the exact continuation is A42/A78's) | a frame was logged after its round's MC move and re-centering, so a continuation from it applied both again; re-centering a centred frame is not the identity in float32 (items 5, 6) | a continuation repeated an accepted move (pivot 1.5e-5 A, jump 0.04 A); under default re-centering continuations left the continuous run by up to 1.1 A (chig) and 1.2 A (1ubq) within 10 frames | a frame is the round's state before its operations (logged, then re-centred, not at the final capture, then the MC move): continuations bitwise under default re-centering and from accepted and rejected rounds (`dfix/logs/restart_recenter.log`, `restart_mc.log`); only frames with an accepted move change; fresh frames no longer exactly centred (up to 0.96 A) |
| A146 | `main.cpp` options | yes | `--time-step`, `--frame-interval`, `--thermostat-timescale`, `--temperature`, `--anneal-factor` and the event intervals not validated (item 8) | `--time-step 0` never ended (an infinite round count); negative or NaN temperatures, timescales and anneal factors gave NaN energies at exit 0; negative intervals gave an event every round or none | refused with the option's name (`dfix/logs/options.log`) |
| A147 | `main.cpp` `ReplicaExchange`, restart; `generate_restart_config.py`; examples' `continue_sim` | yes | replica labels restarted at the identity (item 7; A77 carried the swap set only) | walker labels broke at every continuation segment (from [2,0,3,1], system 0 logged [0,0,2,2,3] for [2,2,3,3,1]); dynamics unaffected | `/input/mom@replica_index` written by `generate_restart_config` and the examples, read by the engine and checked to be a permutation: labels bitwise (`dfix/logs/restart_replica_index.log`) |
| A148 | `main.cpp` `ReplicaExchange::attempt_swaps` | no (A7's momentum swap) | a rejected exchange swapped and scaled the momenta back: sqrtf(T1/T2)*sqrtf(T2/T1) is not 1 in float32 (item 9) | each rejected attempt changed about half the momentum components by an ulp; at T 0.3/3.0 with 0 of 333 accepted the trajectories left the no-exchange run (2.8e-6 A at frame 1, 2.1 A by frame 10) | positions alone are swapped for the energy test; momenta and the label move only on acceptance: bitwise with the no-exchange run (`dfix/logs/rejected_swap.log`); the dynamics of every exchange run move (218 of 4070, hybrid 12 of 196), static outputs bitwise |
| A149 | `deriv_engine.cpp` `integration_cycle`, `main.cpp` round; `thermostat.cpp` | yes | the OU step acted on the stored momentum, which the leapfrog keeps half a kick behind the coordinates, so it also damped that half kick: the coordinates sample exp(-U/T') with T' about T/(1+dt/(2 tau)); under mv the slow half kick is inner_step times larger (driver pass 2 item 1) | configurational temperature 0.09% low at the defaults (tau 5), 0.9% at tau 0.5, 8% at tau 0.05 (tethers 0.92002 +- 0.00099, exact for the scheme 0.91937); mv slow-level tethers 0.79221 (exact 0.79210); the invariance REMD ladder's <U>/(aT) 0.99858 over 12 replica-seeds, seed 3 failing; kinetic energy and g-JF configurations unaffected | the step taken in `integration_cycle` after the first force evaluation, on the momentum at the coordinates, and restored after; g-JF dispatched from `main.cpp` with its thermostat; the duplicated kick-drift loops first merged into `leapfrog_step` (bitwise). After: tau 0.05 1.00079 +- 0.00111 (exact 1.00008), mv fast 1.00035 +- 0.00192, slow 1.00099 +- 0.00192 (exact 1.00073); REMD seeds 1-3 PASS (mean 0.99974); equipartition and pivot MC PASS (`tfix/thermo/`); new `tests/invariance.py thermostat`. Deliberate departure: every protein-only trajectory changes; static outputs and the hybrid bitwise |
| A150 | `rama_map_pot.cpp` RamaMapPot, RamaMapPot2 | no (caused by A57 and A134) | the maps had `set_param` inside `#ifdef PARAM_DERIV` and no `get_param`, while A57 and A134 check every `set_param` against `get_param().size()`, here 0 (local backbone pass 2 D1) | ConDiv's sheet derivative (`training/ConDiv.py:517`, `training/ff31/ConDiv.py:495`) raised "expected 0 but got 51840", so a ConDiv step could not finish on the repo engine; `--set-param rama_map_pot` refused (master ignored it silently). The midway2 campaigns run `$P`'s 10-18 binary, which predates A57 | the maps are the nodes' parameters (`get_param`, `set_param` outside PARAM_DERIV), `get_param_deriv` raises "not implemented": ConDiv difference -0.045776 as with master's library; setting the base map back restores the energy (`lbfix/repro_rama_set_param.log`); regression bitwise |
| A151 | `rama_map_pot.cpp`, `backbone_steric.cpp`, `tension.cpp` readers | yes | `residue_id`, `rama_map_id`, RamaMap2's lambda rows, `backbone_pairs` `id` and `n_atom`, `tension` `atom` and input widths unchecked (D2; the A125/A128/A132 class) | exit 0 on wrong energies (rama `residue_id` 10 of 10 rows 15.71 against 11.03; a phantom 4.0 E_up for `n_atom` 5; tension 0.688 against 0.3435), SIGSEGV for `rama_map_id` 10 or 4000; no writer produces them | validated preconditions: the 13 cases refused, controls unchanged (`lbfix/repro_unchecked_indices.log`); regression bitwise |
| A152 | `upside_config.py` `--secstr-bias` | yes | residue indices unchecked (D4) | residue -1 biased the last residue; one past the end a bare IndexError | refused with the residue count; `chig_rama` unchanged |
| A153 | `upside_config.py` `read_pro_state_file` | yes | a `--trans-cis` table naming a non-proline was accepted, or held only by a bare `assert` (D5) | without `--rama-library` THR got the trans-cis omega double wall in place of its omega spring at exit 0 | refused with the residue's name, the assert removed; `chig_transcis` unchanged |
| A154 | `spring.cpp`; `upside_config.py` `write_rama_map_pot`; up.md 1.3, 1.7 | yes | a commented-out MovingSpring block; `@sheet_eps` 0.01 written with every `--rama-sheet-mixing-energy` though only ConDiv reads it, with the maps `--rama-param-deriv` writes and its own 5e-4; the `rama_map_potential` log documented as (n_frame, n_res) (D6, D7) | none (dead code, a dead attribute); under `--trans-cis` the log has one column per `rama_map_pot` row, without the switching prolines | the block and the attribute deleted; up.md corrected; with A152-A154's writers 25 of 28 fixtures differ only by `rama_map_pot@sheet_eps`, 3 equal (`lbfix/fixtures_cmp.txt`) |
| A155 | `upside_config.py` `write_rama_map_pot` `--secstr-bias` | yes | the helix and sheet windows were built from linear sigmoids in phi (`sigmoid_lessthan(phi, 0)`: 1 at -180, 0 at +175), the helix psi window likewise in psi (local backbone pass 2 D3) | a bias stepped by its full size across phi +-180, inside the beta region: with -1.0 sheet bias on chig residue 1 at psi 130 the energy rose 1.0 E_up from phi -175 to +175, with spline ringing (+1.10 at 172.5) and forces off by up to 9.7 E_up/A; the helix psi step 3.3e-4 of the bias | (decided for the user 10-22) one phi window for both basins, phi < 0 continued past +-180 to +120 (the empty gap between beta and alpha_L), each window the original pair of 10-degree edge sigmoids summed over its periodic images: equal on both sides of +-180 (0.9975274), within 2.5e-3 of the old basins for phi in (-150, 60), 0.005 at alpha_L (0.002 before); force differences through +-180 at most 0.1 E_up/A (`secstr_fix/secstr_periodicity.log`); only `chig_rama` moves, its `rama_pot` and `pivot_moves/proposal_pot` on the three biased residues (`secstr_fix/fixtures_cmp.txt`; `fixtures_w17`) |
| A156 | `monte_carlo_sampler.cpp` `JumpSampler` | yes | an optional `jump_moves@prob_mod`, written by nothing and undocumented, was added to every jump's log acceptance, though both jump moves are symmetric proposals (MC pass 2 D1) | any nonzero value broke detailed balance: at 0.5 every midpoint and orientation observable of 1dfn failed (mid x at T 0.7 X2 2557 on 21 dof; the T 0.7 chain sampled close to T 1.0), acceptance 0.630 against 0.528 | the member and its read deleted, the proposal ratio 0 with the reason: the same config with the attribute passes all 21 observables at acceptance 0.5285 (`mcfix/jump_probmod_0.5.log`) |
| A157 | `monte_carlo_sampler.cpp` Pivot and Jump readers; `main.cpp` | yes | `pivot_atom`, `pivot_range`, `jump_moves/atom_range` unchecked against the system, empty groups accepted (D2; the A151 class) | exit 0 on wrong data (writes past the end; a different trajectory; a read past the end that made two identical runs differ; no-op moves counted as accepted), SIGSEGV for a range 400 past the end and for empty groups; no writer produces them | validated preconditions with the atom count passed to the samplers; the commented-out construction line deleted: the 10 cases refused naming the dataset, controls unchanged (`mcfix/repro_mc_indices.log`) |
| A158 | `monte_carlo_sampler.cpp` PivotSampler; `random.h` | yes | a bin was drawn by a float variate on a float32 cdf but scored with its exact probability p, so the chance s of drawing it is quantized at about 6e-8: per layer up to 2908 bins never drawn but scored, s/p up to 150 (D3) | a pivot-only chain is stationary on pi s/p: worst-residue total variation 6e-6 at T 0.7, 4.5e-5 at 1.0, 4e-4 at 1.3, 4.7e-3 at 2.0; undetectable at T <= 1.3 | the draw on a double cdf with a 53-bit variate (`uniform_open_closed_double`), scored by -log of the same cdf step: agreement 5e-11 relative, TV <= 3e-16; pivot detailed balance PASS (acceptance 0.7438 against 0.7439); full-potential pivot-only chig at T 0.8 acceptance 0.243 +- 0.028 against 0.249 +- 0.033 (`mcfix/stationary/results.txt`). Deliberate departure: pivot draws consume one more random number, so seeded pivot runs change |
| A159 | `martini_hybrid.cpp` `apply_stage_fixing`, `refresh_transition_holds_for_engine`; `martini_fix_rigid.cpp` | no | every round re-registered the rigid group with `initialized=false`, so each projection took the current, already rounded geometry as its reference; errors accumulated, given a sign by float32 quaternions whose w rounds to 1 for small rotations (hybrid pass 2 D1) | the rigid protein of 6.0-6.6 and the 7.0 burn-in drifted linearly, about 1.6e-7 A per round: internal RMSD 6.56e-4 A and Rg +6.5e-4 A after 4000 rounds | the stage constraints set once at registration, each round re-applying only the PO4 z hold; `clear_dynamic_rigid_groups` deleted: internal RMSD 8.6e-7 A after 4000 rounds, no growth (`hfix/logs/rigid_drift.log`); hybrid record bitwise |
| A160 | `martini_hybrid.cpp` PO4 hold; `main.cpp`; `martini_prepare_system{,_lib}.py` `start_production_po4_hold`, `set_initial_position`; `generate_restart_config.py`; `tests/invariance.py`; up.md 1.4, 1.5 | no (A42/A73 incomplete) | the hold counted rounds per run from `@sc_env_transition_step_start`, restart state outside `/input/mom@round` that only the workflow's continuation advanced; `generate_restart_config` also dropped the barostat box (D2, D3c) | a `generate_restart_config` continuation from round 102100 held PO4 z to 102250 instead of 102150 (0.93 A off after 300 rounds), one from 102200 froze it again (1.47 A off after 200 rounds); every positions-only start held 150 rounds (example 16's REMD start and every segment); a continuation started in the config's box (0.277 A off after 200 rounds) | the prep writes `/input/hybrid_control@po4_z_hold_end_round` once (the burn-in's round + 150); PO4 z is held only in a momentum continuation before that round, a positions-only start holds nothing (deliberate departure); the per-run counter, `sc_env_po4_z_hold_steps` and the `continuation` flag removed; `generate_restart_config` carries `/output/box`: both tools bitwise the continuous run inside and after the hold (`hfix/logs/hold_clock.log`). A production file without the attribute is refused: re-prepare stage 7.0 |
| A161 | `box.cpp` MC barostat; `martini_potential.cpp` `DistSpring` | no | the scaled molecule centre included the derived BB/O slots, whose stored positions the integrator never moves; a trial box past the bond spline's fitted reach refitted it, and the refit stayed when the trial was rejected (D3a, b) | moving the stale slots by 5 A left the energy identical but moved real particles (0.0039 A by round 100); the protein's Jacobian per axis 1.006 at s 1.01 against the s the acceptance assumes (about 0.02 bar); NPT continuations 1.1e-5 to 2.3e-5 A off with no move accepted | the centre over the molecule's own particles; the bond spline refit at every box change, from the current box alone: NPT continuations by both tools bitwise (`hfix/logs/npt_box_check_*.log`), the slot test bitwise (`barostat_slots.log`); hybrid record bitwise |
| A162 | `martini_stage_params.cpp`, `martini_hybrid.cpp` readers, `box.cpp` barostat settings, `martini_brownian.cpp`; `martini_prepare_system{,_lib}.py`; up.md 1.4 | no | silent defaults: no stage label meant "production"; `preprod_protein_mode` other than rigid/rigid_body left the protein flexible in the rigid stages (an inactive flag); `activation_stage` defaulted; the hybrid nodes ran without `/input/hybrid_control` (SC-env 0); a BB proxy under another role lost its BB-env pairs; barostat `target_p` 1.0 E_up/A^3 (about 4.8e4 bar), `interval` 0 (registered, never run), `mc_dmax` 0.004 against the writer's 0.01, a dead `mc_init` and `periodic_boundary_potential` fallback; `/input/brownian` without `/input/mass` gave mass 1 (D4) | each ran at exit 0 with another model: a handoff file without the label ran flexible (internal RMSD 0.437 A); without `hybrid_control` total -14687.0 against -19950.1 with the interface off; the box shrank 92.559 to 92.026 A in 200 rounds; lipid <p^2>/(3kT) 3.46 at mass 1 | the stage label, the string attributes and every barostat attribute required; `preprod_protein_mode` removed; every protein atom mapped under its role; hybrid nodes refused without `hybrid_control`; `/input/brownian` needs `/input/mass`: each case refused (`hfix/logs/silent_inputs2.log`, `stage_refusal.log`); hybrid record bitwise |
| A163 | `martini_potential.cpp` `AngleSpring` | no | a term with a zero-length arm was skipped and cos theta clamped to [-1, 1], a guard against CLAUDE.md (D5) | two coincident atoms dropped the term at a finite energy (1040.1593 against 1040.7942) | skip and clamp removed: the coincident case gives NaN as `dist_spring` does, value and gradient on the same cosine; node_fd on `angle_spring`, `dist_spring` 0 bad; hybrid record bitwise |
| A164 | `surface.cpp` rotation scan | yes | `surf_cb` was declared outside the 20-angle loop, so a residue marked at one rotation counted at every later one, and `RotateXY` turned only the first n_tm_res of the 3 n_tm_res backbone atoms (membrane pass 2 M1) | the surface score, read by the channel potential, depended on orientation: a 4.5 deg rotation about the normal changed the energy by 1.86 E_up on 1qhj_channel (76 of 144 scores moved, up to 1.0), 1.63 on 1qhj_channel_m0, 0.15 on 1rkl_channel_m0 | marks reset per rotation, the whole backbone rotated: rotations by 4.5 and 90 deg leave every score and the energy unchanged (noise 5e-4), both methods (`mbfix2/logs/surface_rotation.log`); 1qhj_channel's mean score 0.87 to 0.49, energy -588.28 to -628.03. Deliberate departure: the channel potential's results change; its `icb`/`ihb` tables were fit with the old scores |
| A165 | `main.cpp` `CurvatureChange`, `curved_membrane_node`; `upside_config.py` | no (A37's changer) | the changer ran on a membrane without `use_curvature`, walking a centre no term reads; its measure, flat in curvature 1/R, was improper as R -> 0 with the membrane energy bounded (M4, M5) | 1qhj_memb (flat) logged R from 1000 to 110-1168 at exit 0; 1rkl_curv from R 6 A stayed at 5-17 A, below the 12.4 A half-thickness (no inner leaflet), for 1000 time units | the changer needs a curved membrane and samples radii above its `@half_thickness` (lower proposals rejected, a lower start refused); the writer refuses `--curvature-radius` at or below half `--membrane-thickness` and `--curvature-sign` other than +-1: from 13.5 A R never below 12.4 over 400 frames; the curved run otherwise bitwise |
| A166 | `membrane_potential.cpp` every reader; `surface.cpp` | yes | `cb_index`, `res_type`, `hb_index`, `hb_type`, `hb_res`, surface `index_cb`, `@method_type`, `@exclusion` and the input widths unchecked (M6; the A151 class) | exit 0 on wrong values; `hb_type` 2 read heap memory (energies 5.6e7, -1.4e7, 6e36 between runs) | validated preconditions: the 11 cases refused; bitwise on the membrane fixtures |
| A167 | `upside_config.py` `apply_param_scale` | yes | `--memb-scale` scaled `coeff` but not `coeff_inner` of the channel nodes (M7) | half the channel potential unscaled | both scaled (ratio 2 at scale 2 on both) |
| A168 | `upside_config.py` membrane writers | yes | `--membrane-exclude-residues` removed the residues' CB rows but left their backbone H/O rows; the `hbond_exclude_residues` filter beside them never fired (`infer_H_O` already drops those rows) (M8, M11) | the excluded residues still contributed -5.54 E_up on 1rkl_surface, -3.14 on 1rkl_channel_m0 | their H/O rows leave with their CB; the dead filter and its parameter deleted: they contribute 0 (`mbfix2/results.txt`); 1rkl_surface's `hb_index` 59 to 50 rows |
| A169 | `upside_config.py` membrane options | yes | an out-of-range `--membrane-exclude-residues` was written and excluded nothing (M9) | silently no exclusion | refused with the residue count |
| A170 | `surface.cpp` output columns 1-2 | no (A105 incomplete) | `sfx`, `sfy`, `fx`, `fy`, the area counts and `FORCE_UNIT` fed only the deleted `membranelateral_potential`; reset only for current TM residues and averaged over stale `tm_list` entries (M10) | columns 1-2 history-dependent (65 of 151 residues kept a previous call's values); nothing read them | deleted, one output column: column 0 bitwise with the A164 build (`w_patched.npz` against `mbfix2/w_new_on_w17.npz`) |
| A171 | `membrane_potential.cpp`, `upside_config.py` | yes | `@curvature_radius` was read and overwritten from the centre's z at every evaluation; unused `fs`, `dfs`; a stray `print` (M11) | none | the attribute no longer written or read (`curvature_radius` a local from the centre), the rest deleted; bitwise |
| A172 | `membrane_potential.cpp` `cb_`/`hb_surf_membrane_potential`; `main.cpp`; `upside_config.py`; `tests/fixtures.py` | yes (latent) | a curved channel configuration placed the membrane on a sphere while `surface` found the lipid-facing residues in the flat slab (M12) | 11-15 of 1qhj's residues differed in TM membership at R 150-200; no example used it, the `1rkl_channel` fixture did | the channel potential is flat-membrane only: its curved branch, centre argument and `@use_curvature`/`@curvature_sign` removed, a curved membrane with `surface` refused by engine and writer; `1rkl_channel` rebuilt flat; every channel config written before it, flat ones included (example 08/1's), refused by arity |
| A173 | up.md 1.3, 2.5 | yes | `membrane.h5` documented wrong (`cb_energy` axis 1 a burial bin, not a leaflet; `hb_energy` donor/acceptor by unpaired/paired; `burial_nodes` a burial level; the z attributes missing) and no schema rows for the membrane nodes, `surface`, `cat_pos_bb_cb` (M13) | misleading documentation | corrected and added, with ff2.1's two left defects (M2 hard channel switch at +-half_thickness, M3 `cb_energy`/`hb_energy` end slopes): properties of the trained membrane tables, left for the Upside 3 membrane retraining |
| A174 | up.md 1.7, 2.2, 2.3, 2.6, 2.7, 2.8, 4; `training/ConDiv.py` docstring; `training/README.md`; `martini_build_tables.py` comment | yes (in part) | documentation wrong against the files (parameter pass 2 D2-D6): a default `--bilayer-pdb`/mode=bilayer that no longer exists; ff3.0 said to use `common/rama31.dat` and glycine H-bond offsets 12-14 (the released ff_3.0 has `common/rama.dat` with a trained glycine depth pair and 12 entries; rama31 and the offsets were the unreleased ff30_bio/hbg lines); the depth pair credited to `ConDiv.py`, which copies `rama.dat` unchanged (it came from the run's `trainer/rama_basin.py`); ff3.0 hybrids said possible by passing files (the prep hard-coded ff_2.1); the absent `py/build_rama_from_awh.py` and `run_upside.py:read_output` cited; the pair radial sections labelled (attractive, repulsive); `scale` "positive = buried-favored" (PRO favours burial with both negative); `type_charge` as the bead charges (unread, NC3 0); a comment on the particle tables' end slopes stale since A103 | misleading provenance and documentation | corrected against the files as the parameter review measured them (`pass2/params/report.md`) |
| A175 | `engine_c_library.cpp` `evaluate_deriv`; `main.cpp` setup refusal, `-o`; `h5_support` messages | yes | `evaluate_deriv` caught only `...`, losing the engine's message; a run refused during setup returned before releasing the logger, leaving its file open and locked in the calling process; one `-o` for several systems truncated one file; "attribute X not present" named no path (I/O pass 2 D1, D2, D10, D14) | Python saw "Unable to evaluate derivative"; after one refused in-process run every later library engine failed ("while adding 'protein_hbond', error -1"); "ERROR: error -1" | the message passed on; the logger reset before the return; one output per system required; the path named (`pass2/io/r2/*_mine.log` equal the reviewer's after-logs); bitwise |
| A176 | `h5_support.h` `traverse_dset`, `h5_support.cpp` `read_attribute` | yes | HDF5 converted any stored number to the type read silently (D3; the reader side of A24/A93) | at exit 0: a float index i+0.6 truncated, -0.6 to index 0, `Spring_bond@dim1` 0.7 to 0, `/input/mom@round` 40.7 to 40, int64 3e9 clipped to 2147483647 and the logged rounds wrapped negative | integer reads exact under HDF5's conversion-exception callback (`H5Pset_type_conv_cb`; integral floats pass), attributes read as stored then `H5Tconvert`ed under the same list: each case refused naming the dataset or attribute; no fixture or prep-stage file refused; bitwise |
| A177 | `main.cpp` restart round; `martini_hybrid.cpp` `read_string_attribute`; `martini_brownian.cpp` | no | a negative `/input/mom@round` accepted (`round_num` uint64); no scalar check on string attributes; `@inner_steps` defaulted to 1 and `@reference_temperature_up` to 0 (friction not scaled with T) (D4-D6; A162 class) | round -5 logged -5, -4, -3; an 8-element `@current_stage` overflowed a one-element buffer (ASan WRITE of 80 bytes); a 1rkl 7.1 file without either Brownian attribute ran another integrator or friction (0.48 A, 0.033 A off after 20 rounds) | refused: a negative round, a non-scalar string attribute, a missing Brownian attribute or a non-positive T_ref; bitwise |
| A178 | `main.cpp` frames | no (A49, A145 incomplete) | a run logged its starting round only on the frame grid, so a continuation from an off-grid final capture (`--duration-steps`) did not repeat the junction frame, which every segment reader drops; `--minimize --duration-steps 0` on a continuation logged the minimised state twice (D7, D8) | chig 25 + 25 rounds, frames every 10: joined rounds 0, 10, 20, 25, 40, 50 (30 lost), `mdtraj_upside` times 0.945 and 1.215 for rounds 40, 50; rounds [40, 40]; on-grid junctions (every shipped workflow) unaffected | every run logs its starting round, re-centring only on the grid; the minimised frame only when no round follows: joined 0, 10, 20, 25, 30, 40, 50 with the right times, bitwise the continuous run under re-centering; fresh runs and on-grid continuations bitwise |
| A179 | `run_upside.py` | yes | `OMP_NUM_THREADS` built into `env` and not passed to `Popen` (D9) | a caller at 3 threads asking for 1 got 3 (performance) | passed |
| A180 | `main.cpp`, `state_logger.h` dense logger | yes | `--dense-frame-interval` and the dense-logger machinery: nothing calls `add_dense_logger`, an empty loop ran every round (D11) | none | deleted; bitwise |
| A181 | `rama_map_pot.cpp` RamaMapPot, RamaMapPot2 | no (A150 incomplete) | the constructor fitted the spline to the file's float64 maps while `get_param` returned them as float32, so `set_param(get_param())` was not the identity (D12) | forces up to 9.5e-6 apart after a round trip; ConDiv sets perturbed maps through this path | the spline fitted from the parameters: the round trip bitwise on every rama node. Deliberate departure at rounding: fixtures_w18 357 of 4064 arrays move (statically only `rama_map_pot*`, `RamaMap2*`, `deriv`, max 3.4e-5; the rest seeded runs), hybrid 23 of 196; equal to the reviewer's build bitwise (`iofix/results.txt`) |
| A182 | up.md 1.7 | yes | the `/output` table wrong against the files (D13): `kinetic` is KE per three degrees of freedom, `temperature` the thermostat's, `hbond` (n_frame, n_donor+n_acceptor) counts per donor then acceptor, shapes given with `n_replica` where 1; `replica_cumulative_swaps`, `replica_swap_partner`, `tip_pos`, `time_estimate`, `protein_kinetic`, `lipid_kinetic`, `sigmoid_coupling`, `Const3D_<k>` and the extensive-level datasets missing; the engine deletes a previous `/output` (the continuation scripts rename it) | misleading documentation | corrected from the real outputs (`pass2/io/r2/output_shapes_repo.log`) |
| A183 | `martini_prepare_system_lib.py` `convert_stage` | no | each atom wrapped into the box alone, while `set_box_from_lipid_xy` puts the protein's lowest atom on the box edge whenever the lipid extent sets the box (the CLI's `--box-padding-z 0`, used by `run.py` and the glpG preps) (prep pass 2 D2) | 14 of 40 ion seeds wrapped glpG's N-terminal N: a 122.32 A N-CA bond at exit 0 (Upside's protein terms are not periodic: Spring_bond 1.862 to 2.615e5 E_up on the 1rkl 6.0 checkpoint); the shipped preps escaped by 0.13 A | each molecule wrapped whole by its centre: that seed's largest bond 1.549 A; 1rkl molecule ids equal, no protein position changes (`pass2/prep/patched/verify_python.log`) |
| A184 | `martini_prepare_system_lib.py` `_group_atoms_into_molecules` | no | consecutive atoms grouped, so with the BB proxies following every chain's N/CA/C/O each chain's proxies became a molecule of derived slots only (D3) | a multi-chain protein could not run NPT: 1AFO stopped at 6.2 round 10 ("a molecule of /input/molecule_ids holds only derived hybrid slots") | protein molecules keyed by chain: 1AFO two 180-atom molecules, 6.2 runs; 1rkl's molecule ids bitwise |
| A185 | `martini_prepare_system_lib.py` `inject_particles_table`, `soft_core_radius` | no | the soft core was tabulated on [0.3, 12] A, the spline returning the edge value with zero slope below, though the soft core is finite at r = 0 (D4; dry-MARTINI pass 2 U1) | pairs under 0.3 A felt no force in 6.0 where the soft core gives 118 to 2.6e5 E_up/A; reached: `run.py`'s CHARMM-GUI template has 2 pairs under 0.3 A and its 6.0 start 1 (no measured effect on the outcome: every pair 4.4 A or more after 6.0) | the soft core tabulated from r = 0 as (alpha sigma^6 lambda^2 + r^6)^(1/6): the engine's force equals the analytic one at 0.10, 0.20, 0.25 A (V3); the hard tables keep 0.3 A (unreachable after 6.0, and the core check refuses such a pair). Deliberate departure for 6.0's tables |
| A186 | `martini_potential.cpp` `martini_potential`, `martini_sc_table_1body` | no | a negative `@cache_buffer` clamped to 0 (D5; its Brownian part is A177) | silently rebuilt the pair list every call | refused; bitwise |
| A187 | `example/16.MARTINI/run.py` | no | `membrane_thickness = 48.8`, the retired DDM slab's (commit 0d4a1c57), used for POPE/POPG (D7) | 68 salt pairs on `run.py`'s route instead of 74 at the POPE/POPG P-P thickness | 41.8 A, findings 8's GROMACS P-P, as the 10-21 glpG prep passed. Open: the CLI fallback takes the bead z extent (62 A for 1rkl) and the readme says "z-extent"; one definition is for the user (plan.md Phase 24 item 5) |
| A188 | `martini_prepare_system_lib.py` `derive_chain_break_metadata`, terminal bead types; `run_stage70_handoff` docstring; the OPM error text; findings 2.4, 2.5 | no | a numbering gap inside a chain was counted as a second segment of that chain, which the mapping validator refused with a misleading message, and its ends typed as charged termini (D9); the docstring called production's lipids free (PO4 is held 150 rounds), findings 2.5 refreshed counts at a continuation (none are), the OPM message said "A micelle build" (D10) | a gap refused at stage 0 as "chain_counts length must equal number of chains" | a gap refused by name ("chain A skips from residue 14 to 16; model the missing residues before preparing a hybrid"), terminal types only at chain ends; the texts corrected; regression bitwise |
| A189 | `parameters/ff_2.0/`, `parameters/packing/`, `parameters/dryMARTINI/DDM.pdb`; up.md 2.7, 6 | yes (in part) | content nothing reads: ff_2.0 (no reader in py, examples, tests or training since its one reader went), packing (no reader; the engine refuses its tables), DDM.pdb (the prep refuses DDM) (parameter pass 2 D7) | none | deleted (user 10-22: "Delete all"); in git history |
| A190 | `martini_prepare_system.py` ion count; `example/16.MARTINI/readme.md` | no | without `--membrane-thickness-angstrom` the ion count took the packed lipids' bead z extent (62 A for 1rkl) while `run.py` and the glpG preps pass the P-P thickness (prep pass 2 D7) | about 37 salt pairs on 1rkl against about 51 at P-P | the ion count uses the P-P thickness (user 10-22): given, or measured as the distance between the two leaflets' mean PO4 heights of the packed lipids; the readme states it with findings 8's values (POPE/POPG 41.8 A, DOPC 43.8 A with salt) |
| A191 | `martini_prepare_system_lib.py` `map_backbone_types_from_structure` | no | the MARTINI backbone bead types came from a homemade helix-only detector (a phi/psi window or CA(i)-CA(i+-4) 4.8-6.8 A), so strands were typed coil and helices merged across kinks (prep pass 2 D8) | off martinize's types: 1UBQ 54 of 76 residues, 1dfn 48 of 60 (strands P5 for Nda), glpG 84 of 210 (one 41-residue "helix"), 1AFO 14 of 72, 1rkl 6 of 31 (`dssp_bb/compare.log`) | DSSP (mdtraj's implementation) on each chain's backbone, typed by martinize's `ssClassification`, as `martinize -dssp` does (user 10-22: "Use DSSP"); the old detector and its helpers deleted. Changes every hybrid's BB-env interface; verified with the end-to-end re-prep (A192) |
| A192 | `martini_build_tables.py` side-chain/environment table; `martini_prepare_system{,_lib}.py`; `parameters/ff_2.1/martini.h5` (old one `martini.h5.bak_pre_martinize_beads_20261022`); up.md 2.6 | no | the table placed an idealised bead chain on Upside's trained rotamer point, aligned with its trained angular vector (`sidechain.h5 rotamer_center_fixed[:, :6]`, free trained interaction sites), not on the side chain's atoms; its azimuthal average used Boltzmann weights at 25 T_up, no recorded reason; four unread datasets (parameter pass 2 D1, D7; prep pass 2 D1, D6, U1) | ARG's charged bead 0.5-2 A from CB on the backbone side, 2.5-5 A from martinize's centre; LYS, ASP, GLU, PHE, TYR, TRP off by 0.5-5 A; median energy errors ARG 2-6, ASP/GLU 1-4 E_up against beads from crystal side chains; maximum 1.09e18 E_up from the misplaced bead | (user 10-22) each state's beads from ideal-geometry heavy atoms (master's `generate_from_rotamer.cpp` geometry) at its library (chi1, chi2), chi3/chi4 at their staggered wells (ARG chi4 -85, 85, 180; GLU/GLN chi3 every 30 deg), syn-pentane pairs left out, PHE/TYR chi2 wrapped into +-90 as PDB naming does, 546 conformers, martinize's mass-weighted martini22 centres; one free energy at 0.8647 T_up over conformers and the 2400 x 24 azimuthal directions, `exp(-F/kT) = sum_c w_c <exp(-E_c/kT)>`. Against 1146 crystal side chains: best-rotation RMSD median 0.18 A (old 1.38); energies against the crystal ensembles' free energy median 0.07 E_up per row (old 2.31); 289 identical-ITP layers and the particle tables bitwise; 1rkl 7.0: only the SC table term moves (-34.3 E_up over residues, PHE10 -9.6), forces finite (`sc_rebuild/`). Every hybrid must be re-prepared (fixtures_hybrid_w4; the cluster's live files before Upside 3 is deployed there) |
| A193 | `training/ConDiv.py`, `training/ff31/ConDiv.py` `zero_for_sarw`; README | no (Kleinmann's port) | the SARW (unfolded reference) replica kept the backbone desolvation term (`bb_sigmoid_coupling_environment`), while the SI defines it "without H-bonds, side chain and solvation terms" and master's `apply_param_scale` scales the backbone term with `env_scale`; Kleinmann's `condiv2.py` has a separate `bb_scale` it never passes, which the repo copied (ConDiv pass 2 D1) | the DSE reference carried a solvation term the SI's does not: 1rkl SARW Rg 17.1 A with it against 16.0 without (EHEE 18.6 against 18.4); at the native -5.9 and -7.3 E_up | (user 10-22: "Follow the SI") the backbone term scaled with the others in both trainers: the SARW config's total equals `apply_param_scale(0)`'s (224.0758, 28.2651), the backbone term 0 (`pass2/condiv/logs/verify_sarw_patch.log`). New campaigns only; running ones keep their own copy |
| A194 | `training/ConDiv.py` `main_gate` | no | an empty group (`hbg` with a 12-entry `hbond.h5`) counted in the Bonferroni cut (D2) | the cut 0.0050 for 0.0056 (family-wise 4.5%); no recorded verdict changes | groups with no parameters skipped: 9 groups, cut 0.0056, every p unchanged |
| A195 | `training/train_chain.sbatch`; `training/ConDiv.py`, `training/ff31/ConDiv.py` initialize | no | the chain script hard-coded 19 minibatches per epoch, while `initialize` builds `len(training_set) // 24` after its chain-break exclusion (D5) | an 18-minibatch epoch would put the step count one ahead per epoch (19 finished steps read as 20) | `initialize` writes `run_output/minibatches_per_epoch`, the sbatch reads it, and the step count is the number of finished steps' checkpoints: 19 steps of an 18-minibatch epoch read 19 |
| A196 | `training/README.md` | no | the frozen backbone and H-bond parameters said to have no engine derivative (exact since A29/A30); `upside_input/rama_reference.pkl` listed though the workers read the tree's (D3, D4) | misleading documentation | corrected |
| A197 | `engine_c_library.{h,cpp}`, `upside_engine.py`, `deriv_engine.{h,cpp}`, `tests/{regress,invariance,line_integral,node_fd}.py` | yes | dead API: four clamped-spline helpers called by nothing (one read `a[-1]` on an empty array: ASan heap overflow at `spline.cpp:148`), `get_bounded_value_and_deriv_2d` declared nowhere and exported mangled, a `quiet` option never read, unused imports, an empty `try/except OSError: raise` (glue pass 2 D2) | none in use | removed (the wrapper and the library deploy together); bitwise on fixtures_w18 and the hybrid |
| A198 | `upside_engine.py` `Upside` | yes | it required `/input/sequence`, which nothing reads (D3) | a config the binary runs (exit 0) refused by the wrapper | the read removed: chig without it 5.9421 as with it (`glfix/results.txt`) |
| A199 | `rotamer_parameter_estimation.py` `AdamSolver`; both ConDivs | D4 no, D5 yes | `update_for_d_obj` always returned zeros, master's Nesterov hook for a solver the branch deleted, added to the parameters before the look-ahead files; `_read_comp` treated an IndexError as a scalar (D4, D5) | none (the files identical with and without the addition); a per-component array one short gave the last component an array-valued update | the hook and both additions removed; a short sequence raises; 20 Adam steps with ConDiv-shaped settings bitwise |
| A200 | `example/00.AnalysisScripts/helpers/analyze_remd.py` | no | MBAR on the raw coupled energies (-20,500 E_up) without a reference subtraction never left uniform weights; `--cv rg` was over every particle, lipids and ions included; a target between rungs took the nearest rung's weights; any MBAR exception fell back to a histogram (hybrid analysis pass 2 D1) | the saved FES identically 0 while the tool printed "MBAR FES at kT=0.8700 reconstructed"; Rg 42.6 A for a protein CA Rg of 12.9 | a pooled-mean reference, reweighting to the target, the protein's CA Rg, an off-ladder target refused, fallback only without pymbar: CV 12.90-12.94 A, ESS 96.9 of 117, FES 0 to 4.16 with 3 of 40 bins empty (`pass2/hdx/analyze_remd/patched.log`; empty bins: A214) |
| A201 | `example/16.MARTINI/run_hdx_analysis.sh`, readme 4 | no | without `HDX_T_UP` a REMD analysis took replica 0's temperature, a ladder's bottom rung, and the readme's 48-replica command set none (D2) | every uptake, stability, plot and CSV temperature at 0.70 T_up (245 K) instead of the 0.8647 target on `run.py`'s ladder | a `{replica}` pattern without `HDX_T_UP` refused; the readme command sets it |
| A202 | `main.cpp` replica label; `example/16.MARTINI/run.py`; up.md 1.1 | no (A147 incomplete) | the label was read only in the momentum branch, and `run.py`'s position-only REMD continuations wrote none (D3) | walker labels reset at every hybrid REMD segment (2 of 3 slots relabelled at a bitwise-equal junction frame); dynamics unaffected | read whenever `/input/mom` carries it, and `run.py` writes it: labels continue; bitwise on fixtures_w18 and the hybrid |
| A203 | `example/00.AnalysisScripts/helpers/calc_hdx_ht.py` `residue_dg_profile_with_error` | yes | `dG_err` started at zeros, so where no jackknife is computed (one replica, the hybrid default) the error read 0 (D4) | every resolved residue written and plotted +- 0.000 | NaN where no error is computed; dG and the censoring bitwise; computed errors bitwise |
| A204 | `rotamer_parameter_estimation.py`; both ConDivs; `training/README.md` | D1 no, D6 yes | the regulariser and the latent's gradients ran in torch, whose libomp joins the engine's in one process; the latent carried an 86-coordinate `rotscalar` block that nothing reads and d_obj gives no gradient (glue pass 2 D1, D6) | importing the engine before torch exits 139 on the Mac (`pass2/build/omp_import_order.sh`); every checkpoint carried the dead block | numpy with closed-form gradients (`unpack_params_vjp`, `quadspline_energy_vjp`, d_obj), the block removed (31,420 to 31,334 coordinates; an older checkpoint is refused by `restart` and `extract`): against recorded torch outputs the tables agree to 3.1e-16, the gradients to 3.9e-16 of their scale, the packed latents bitwise; a seeded two-protein ConDiv step bitwise but `rot` (1.4e-13, float32 tables bitwise); one libomp in every import order (`notorch/`) |
| A205 | `src/CMakeLists_Other.txt` | yes | `find_package(OpenMP QUIET)`; a `CMAKE_LD_FLAGS` that CMake does not have; commented-out header lines (build pass 2 B1, B2) | without OpenMP, configure passed and every compile failed on a flag named `NOTFOUND` (`pass2/build/repro_cmake/run.log`) | REQUIRED (configure stops with "Could NOT find OpenMP"), the dead lines removed; GCC 16 link lines unchanged |
| A206 | `src/CMakeLists_Other.txt` | no | `engine_c_library.cpp` in `ENGINE_SRC` and named again in `add_library` (B3) | none: CMake drops the duplicate, moving the object to the front of the library's link | listed once, as in the Mac file; GCC 16 builds before and after bitwise (`pass2/build/reg/compare_*_gccA_gccB.txt`) |
| A207 | `src/CMakeLists_M1.txt`, `source_arm` | yes | the Mac file overwrote `-DEIGEN3_INCLUDE_DIR` with a hard-coded Cellar path (eigen@3 3.4.1) and added `opt/eigen@3` again, so the option was ignored and the template's `brew --prefix eigen` (Eigen 5, which needs C++14) never showed (B4) | the build breaks once Homebrew moves eigen@3 past 3.4.1 | the option used, as SYSTEM, and the template passes eigen@3 |
| A208 | `src/main.cpp` | yes | `omp.h` included though no file calls an `omp_*` function; the Mac build found the header only through `source.sh`'s `CPLUS_INCLUDE_PATH` (B5) | `make` in a shell without `source.sh` stopped at `'omp.h' file not found` | the include removed; that compile passes in a clean environment |
| A209 | `source_arm`, `source_sh` (the generated `source.sh` with them) | yes | `py:$PYTHONPATH` with `PYTHONPATH` empty left an empty entry, the working directory, ahead of site-packages (also `CPLUS_INCLUDE_PATH`, `LIBRARY_PATH`, `LD_LIBRARY_PATH`); a stale block for an absent python@3.10 keg put `/bin` first on `PATH` (B7) | a `tables.py` in a run directory imported in place of PyTables | `${VAR:+:$VAR}`, as `training/env.sh` and `env_shared.sh` already do; the block removed; the local `source.sh` regenerated from the template (its old copy `buildfix/old/source.sh`) |
| A210 | `training/env.sh` | no | outside any tree it returned 0 with `UPSIDE_HOME=/` (B8) | a run failing later, far from the cause | refuses with a message (return 1); inside a tree the variables are unchanged |
| A211 | `install_python_env.sh` | partly | biopython (`martini_prepare_system.py`), numpy and scipy imported directly but installed only as other packages' dependencies; pydssp and colorcet installed though nothing imports them (build U5) | none today | the direct imports listed and tested, the two removed; jax kept and named as pymbar's solver backend, which pymbar loads whenever jax is installed |
| A212 | `src/vector_math.h` `compact_sigmoid`, `compact_sigmoid_cutoff` | yes | a `NONCOMPACT_SIGMOID` branch that no build defines (build pass 2, section 4) | none | removed; bitwise on fixtures_w18 and the hybrid (`buildfix/results2.txt`) |
| A213 | `.devcontainer/` | yes | the container is Upside 2's: its Dockerfile clones GitHub master into `/upside2-md` and builds it, `entrypoint.sh` points `UPSIDE_HOME`, `PATH` and `PYTHONPATH` at that clone, and `post-create.sh` builds the workspace as root; its conda environment pins numpy below 1.26 and lacks scikit-learn, pyhdx and jax (build pass 2 B9) | a developer in the container edits the workspace and runs master's code and binary, with no error | deleted (decided for the user 10-22: Upside 3 runs through `source.sh` and `env_shared.sh`; nothing referred to the directory; master keeps it for Upside 2; a copy is in `buildfix/old/devcontainer/`) |
| A214 | `example/00.AnalysisScripts/helpers/analyze_remd.py`; `helpers/calc_hdx_ht.py` H-bond FES plot | D1-D3 no, U3 yes | kT read from `/output` alone while every `output_previous_*` segment was pooled under it; half-open bins, so the frame at the CV maximum fell in no bin; empty bins floored at 1e-12, a finite -kT ln 1e-12 (analysis pass 2 D1-D3, U3) | a REMD continued on another ladder analysed at the last ladder's labels (2408 of 3208 frames sampled at 0.80-0.94 labelled 0.90-1.04, exit 0); the top bin 5.004 against 4.415 E_up; 15 empty bins of 200 read 20.81 E_up beside a largest estimate of 4.26 | a file sampled at more than one temperature refused; `np.histogram`; an empty bin is +inf: the top bin equals a float64 MBAR oracle, the empty bins inf, the hybrid file still runs (`afix/results.txt`) |
| A215 | `example/00.AnalysisScripts/0.run_HXMS.py` | yes | each fitted uptake curve sampled on a log-time grid of as many points as the CSV has peptides; 5.analyze's cooperativity factor and `normalize_se` work on that grid (D4) | the same peptide's COF 2203, 2635, 4446, 4455 with 4, 12, 40, 150 peptides in the CSV (converged 4455.3) | a fixed 200-point grid: 4455.0 in every case |
| A216 | `example/00.AnalysisScripts/5.analyze_D_uptake.py` | yes | peptide Start/End, PDB numbers since A81, used as FASTA positions; residue colours saved in peptide-row order beside colours in sorted order (D5) | for numbering from 87, `COF_values.csv` names residues the sequence does not have (0 of 32 letters match); with peptides out of order 30 of 32 colours belong to other residues | positions through `load_pdb_numbering` (PDB numbers exported as `pdb_resnum`), colours saved with their own residues: 32 of 32 both ways |
| A217 | `example/00.AnalysisScripts/4.calc_D_uptake.py`, `analysis.sh` | yes | the temperature fit skipped every peptide starting at or before residue 11, glpG's tag, for any protein (D6; P2-06's U6) | on chig every peptide skipped: error 0 at every T, `T_optimal` the search bound 0.93963, the EXP uptake written there, exit 0 | the setting `T_fit_exclude_through` (0; glpG's 11 in `analysis.sh`'s checkme block), and a fit with no peptide left refused |
| A218 | `example/00.AnalysisScripts/6.generate_hx_plots.py`, `analysis.sh`, README | yes | glpG's construct sequence, residue offset, helix ranges, temperatures, 296 K and pD 9 as defaults that `analysis.sh` never overrode (D7) | chig's plots computed from glpG's 210-residue sequence, exit 0 | the sequence from `inputs/<pdb>.fasta`, the experiment settings required, glpG's values in `analysis.sh`'s checkme block: chig's 10 residues reach `getkchem` |
| A219 | `example/00.AnalysisScripts/helpers/hxfunctions_clean.py` `getkchem` | yes | pKD2O held at 15.05, its 293 K value, at every temperature (D8) | intrinsic rates low by 0.883 at 296 K and 0.135 at 323 K against `4.calc_D_uptake.py`'s (A64) | pKD2O at T from the same table: ratio 1.000 |
| A220 | `example/00.AnalysisScripts/README.md`; `py/mdtraj_upside.py` | no | a metadynamics helper and node documented that do not exist; a `.up` loader registered with MDTraj that `md.load` cannot use and nothing calls (analysis pass 2, minor) | misleading documentation; dead code | removed; `load_upside_traj` unchanged |
| A221 | `example/04.HDX/4.calc_HDX.py` | yes | a residue protected in every frame took dG 1000 in the temperature loop, and the denaturant loop's log(p/(1-p)) was unclipped, the defect A64 fixed in example 00 only (examples pass 2 D1; analysis U5) | on EHEE_rd2_0005 18 of 39 residues plotted at 1000 kcal/mol over T and at 20.09 kcal/mol over [den], a value set by float64 round-off of 1 - sum(w); an exposed residue -inf | one `residue_dG`, censored to NaN beyond one effective frame as in `calc_hdx_ht.py`: resolved values unchanged to 1.6e-13 |
| A222 | `coord_const.cpp` MovingConst3D; `upside_nodes.py`; `generate_restart_config.py`; example 06; up.md | yes (A130 incomplete) | the tip clock time_initial + time_step*n restarted n at 0 in a continuation, moving time_initial to the logged float32 clock and shifting start_pos: two roundings where the continuous run has one (examples pass 2 D2) | a pulling continuation one float32 ulp off on the clock and tip, then 7.4-8.8 A from the continuous run within 50 + 50 time units (example 06's own path: 2.9-29 A over 1000 + 1000); the tip schedule right to 1e-6 A | a required `@step_initial` (0 in a fresh config) and `/output/tip_step`; a continuation keeps start_pos and time_initial and sets the step count: bitwise at 3 + 3, 50 + 50, 300 + 300, moving and still; an old pulling config refused; fixtures_w19 (`efix/results.txt`) |
| A223 | `run_upside.py` `upside_config` | yes | the wrapper passed `--intensive-memory` for True and nothing for False, and the writer turns the mode on above 1000 residues (D3) | example 11's documented False/True comparison wrote one config (7 edge budgets either way, 7.0 GB both) | None keeps the writer's choice; True and False pass their flag (0 budgets with False) |
| A224 | `advanced_config.py` help of `--ask-before-using-AFM`, `--tension` | yes | spring constants said to be in kT/A^2 (41.14 pN/A at 298 K) and tensions in kT/A; the engine's terms are in E_up (D4) | a spring constant read 15% low as pN/A | E_up/A^2 and E_up/A, 48.40 pN/A per E_up/A^2 (A62); the 5-field table described |
| A225 | readmes and comments of examples 01, 02, 05-09, 12, 13; example 06 `1.get_force.py` | yes | file names that do not exist, a done TODO, steps for time units, three membrane scripts undocumented, "Example 8" for the restraints example, "switch to 1dfn" where 1dfn is absent or refused, a column of tip extension labelled tip position (D5) | misleading documentation | corrected |
| A226 | `martini_potential.cpp` `martini_potential`, `martini_sc_table_1body`, `dist_spring`, `angle_spring` | no | without `@x_len/y_len/z_len`, `martini_potential` fell back to `@wall_*` attributes no writer has written since 2025-07 and then to a zero box, running the MARTINI pairs without periodic images; `dist_spring` and `angle_spring` kept the dead wall branch (registry closing D1) | a hybrid without its box ran at exit 0 with `martini_potential` -19580.78 against -23341.35 E_up; no repo writer produces one | the box required and positive in all four nodes, the fallbacks and the non-periodic pair-list path removed: a box-less hybrid refused; ff2.1 bitwise; the hybrid moves at float32 rounding (1e-7 to 2.5e-6 of each array's scale) from fast-math code generation alone (bitwise built without `-ffast-math`); node_fd clean (`gfix/`); hybrid baseline h53 |
| A227 | `martini_prepare_system_lib.py`; `upside_config.py`; up.md 1.4 | partly | datasets and attributes no reader uses: `martini_potential` `@initialized`, `@potential_type`, `@coefficient_indices_reference`, `@n_combined_triples`; `dist_spring`, `angle_spring`, `restraint_position` `@initialized`; `martini_sc_table_1body` `restype_order`, `target_order`; the hybrid's placement `rama_residue`, `beadtype_seq`, `id_seq`, `fix_rotamer`; ff2.1's placement `fix_rotamer` and `rama_residue` (left from master's deleted `placement_point_vector_only`), `placement_fixed_scalar` `rama_residue`, `rama_map_pot@restype` outside `--rama-param-deriv`, `cb_membrane_potential@half_thickness` (up.md said it bounds the curvature changer, which reads `hb_membrane_potential`'s) (registry closing D2 and section 4) | none (dead records; misleading documentation) | not written; up.md names the readers (the changer reads the hb node; `rama_coord/id` and `rama_map_id_all` feed the pivot writer); fixtures_w20 and fixtures_hybrid_w5, records bitwise; the registry matches 44 of 50 node types, the other 6 holding only items read outside the engine (`pass2/registry/report.md`) |
| A228 | `rama_map_pot.cpp` RamaMapPot2 | yes | `@log_pot` and a per-residue potential array read and filled, never used | none | removed; bitwise |
| A229 | `membrane_potential.cpp` `membrane_potential`; `main.cpp` | yes | the old membrane node that A41 kept so old configs would load has no writer and no possible input: master's writer raised KeyError under Python 3 and it reads a `membrane.h5` layout neither tree ships; no `.up` in either tree holds it (99 HDF5 files searched; registry closing D3) | dead reader | deleted (decided for the user 10-22, superseding A41's choice under the clean-slate rule); bitwise |
| A230 | `tests/line_integral.py`, `tests/registry.py`, `tests/index_oracle.py` | no | the line integral stepped one slot of a tied parameter, which the engine refuses as a broken angular tie, and could not resolve a derivative whose energy change is below float32 at the default width; the registry traced the 10-19 `fixtures_w3` set; the index oracle kept a branch for master's `placement_point_vector_only` | closing checks that could not run or ran on old configs | the line integral follows node_fd's tied direction (untied cases unchanged) and takes an optional width; the registry reads fixtures_w20 and fixtures_hybrid_w5; the branch removed |
| A231 | `interaction_graph.h`, `bead_interaction.h`, `rotamer.cpp`, `environment.cpp`, `backbone_steric.cpp`; `upside_config.py`, `martini_prepare_system_lib.py`; `tests/index_oracle.py`; up.md 6 | yes | the pair lists excluded sequence neighbours by the elements' residue ids, so residues 1-2 apart across a chain break were excluded as if bonded (D9; decided by the user 10-22, design (a)) | missing from every multi-chain run; on 1dfn `protein_hbond` 1 pair, `hbond_coverage` 33, `hbond_coverage_hydrophobe` 72, `environment_coverage_sc` 24, `environment_coverage_hb` 33, `hb_environment_coverage_hn`/`_oc` 3 each, all among residues 28-29 against 30-31 | every node excluding by sequence distance reads a required `sequence_id` (`sequence_id1/2` when asymmetric): the residue number plus 2 per preceding chain break (A53's convention, one `upside_config.sequence_id` shared with the hybrid preparation; widest window 2); the pair list carries only that key; a config without it refused by name. fixtures_w21, fixtures_hybrid_w6: single-chain records and the hybrid bitwise, only 1dfn's and 1dfn_cavity's seeded runs move; index oracle on the three 1dfn configs 7 mismatches to 0 (`scratchpad/upside3_audit/d9/`, `d9fix/`) |
| A232 | `example/16.MARTINI/run.py`, readme 2-3 | no | Stage 2 was one unchunked Slurm job asking 48 replicas for 200,000 time units under a 35:45 wall (about 170 days at the measured ~28 time units per hour); Stage 1 ran the whole preparation (hours of hybrid dynamics) in the calling process, a login node on the cluster; `swap_table2d` returns one swap set for two replicas and the script indexed a second; the readme described caslake (user 10-22: "as we changed the upside core, the job scripts should also be updated to support it") | the documented REMD could never finish, and was killed at the wall without a clean last frame | `run.py` submits the preparation, then REMD blocks that stop at `--time-limit` before the wall and queue the next up to `n_blocks`, continuing between blocks; a block counts only once its engine run succeeds; every swap set passed; each job sources `job_env_file`; `continue_sim` and `duration` gone (raise `n_blocks` to extend); `frame_interval` 10 (~100 frames per block). Two local 1rkl blocks: block 2 starts every slot at block 1's last frame with its walker and temperature (`scratchpad/upside3_audit/runpy_test2/`) |
| A233 | `example/16.MARTINI/run_sim_1afo_full.sh` | no | the batch job sourced the shared deployment's `env_shared.sh`, the pre-audit engine | an Upside 3 job script ran Upside 2 | sources the Upside 3 tree's `env_upside3.sh` |
| A234 | `martini_prepare_system.py` (`--force-field`, `normalize_hybrid_workflow_args`, `validate_hybrid_workflow_args`), `martini_prepare_system_lib.py` `inject_stage7_sc_table_nodes`; `example/16.MARTINI/run.py`, `run_sim_hybrid.sh`, readme | no | the hybrid preparation could run only ff2.1: it hard-coded `common/rama.dat` and ff_2.1's sheet, hbond.h5 and sidechain.h5 (user 10-26: run glpG with ff3.1) | no other force field could be tested on the corrected hybrid, TM4's loss among them (3.15) | a required `--force-field` (`FORCE_FIELD`) under `parameters/`: its rama.dat if it ships one (the benchmark's rule), else common's; its sheet, hbond.h5 and sidechain.h5, the last passed into the injection; martini.h5 stays ff_2.1's, built from ff_2.1's side chains, and a force field with other side chains is refused (refined by A235). Old and new injection on one stripped glpG stage file: ff_2.1 identical in all 614 `/input` entries, ff_3.1 differs only in `rama_map_pot/rama_pot` (`scratchpad/ff31_merge/test_force_field.py`) |
| A235 | `martini_prepare_system.py` `validate_hybrid_workflow_args`, `rotamer_library_differences` | no | A234's check compared the whole `sidechain.h5` with ff_2.1's, but `martini_build_tables.py` builds the side-chain/lipid table from three of its datasets, `restype_order`, `rotamer_start_stop_bead` and `restype_and_chi_and_state` (user 10-26: "do submit ff3.0") | ff3.0 was refused though its side chains build exactly ff_2.1's table (its retraining changed `pair_interaction`, `coverage_interaction`, `hydrophobe_interaction`, `hydrophobe_placement` and `rotamer_center_fixed`, none of which the builder reads), so the TM4 test against Upside 2's placement-fixed ff3.0 chains (3.15) could not run | the check compares those three datasets (`SIDECHAIN_TABLE_DATASETS`) and names any that differ: ff_2.1, ff_3.0 and ff_3.1 pass, locally and in the Upside 3 tree; a copy of ff_2.1's with one chi row or one bead count changed is refused, naming that dataset |
| A236 | `martini_potential.cpp` `DihedralSpring` (`dihedral_spring`) and the box updater; `main.cpp` dry-MARTINI terms; `martini_itp_reader.py` `[ dihedrals ]`, `[ constraints ]`; `martini_prepare_system_lib.py` `_load_martini_topology`, `_build_bonded_terms`, the writer; `dryMARTINI_itp/dry_martini_v2.1_DDM.itp`; up.md; `tests/index_oracle.py` | no (a new path) | the engine had no dihedral term, so the repo's DDM ITP dropped the published topology's three head dihedrals (GROMACS type 1); the reader's dihedral entry wrote its force constant over its third atom (key `k` twice); `[ constraints ]` were skipped without notice (user 10-26: "you should implement the missing dihedral that is needed") | DDM's maltose head would have run without its torsions, off the published model | `dihedral_spring`: V = k (1 + cos(n phi - phi_s)) from one canonical spline of 1 + cos(u), phi as GROMACS defines it; ITPs read with `FLEXIBLE` defined (constraints enter as their FLEXIBLE bonds; an active constraint or a dihedral of another type refused); the DDM ITP is cgmartini's `martini_v2.0_DDM_01.itp` verbatim (md5 `0d587b39`). Against GROMACS 2026.4 on 5 frames of 200 DDM (600 dihedrals over phi -180 to 180): dihedral, bond and angle energies to 1.2e-6, 3.1e-6, 7.8e-7 relative, dihedral forces to 1.3e-5, all bonded forces to 1.2e-4; `node_fd` 0 bad of 320 directions per frame; `regress` 0 of 4260 arrays differ (`fixtures_w21`, `fixtures_hybrid_w6`); `index_oracle` on 1rkl stage 7.0 0 of 120 mismatch (`scratchpad/ddm_micelle/dihedral_check/`) |

Deliberately not changed: the H-bond radial cutoff of 3.5 A truncates the factor at 3.4e-4, which is below the
1e-4-scale energy noise BP already carries. Moving it to 4.0 A would cut the jump to 6e-6, at about 1.5x the
donor-acceptor pairs in a node that takes 13-19% of runtime.

Every parameter derivative in Upside 3 is exact or raises "not implemented" (A29-A31). Refuted 10-21: `--use-heavy-atom-coverage`
writes a node no potential reads, but `get_protection_state.py` reads its output for HDX, so it is not dead. Which parameters ConDiv
trains is decided in its learning-rate table. Open items are in plan.md Phase 24.

### 3.15 Upside 3's glpG hybrid with ff2.1 loses TM4 when the protein is first released (2026-10-26; answered 10-27: ff2.1's own, not Upside 3's)

The validation run of plan.md Phase 24 2b.4 (`upside3_validation/glpg`, remote_jobs.md §1 line 3), read against
Upside 2's hybrid with the same code (`upside3_validation/glpg/glpg_u3_block1.py`, `glpg_u3_prep_tm4.py`; outputs
`glpg_u3_vs_u2_block1.txt`, `glpg_u3_prep_tm4.txt`). Upside 2's reference is the pre-audit ff3.0 79HIS chain's block 1
(groups 1-17, every other one) and the recorded ff2.1 numbers on the fixed table (1.35: the cancelled `ff2.1_sc`
REMD chain, whose files are deleted, and 4.15's local seeds). **No ff2.1 reference shares Upside 3's side-chain
placement** (user, 10-26): both ff2.1 sets ran with the placement defect (4.17, fixed 10-17), which misplaced the
beads read by the rotamer pair term and the side-chain/backbone H-bond coverage, the protein-internal packing that
holds a buried helix such as TM4 (on the old side-chain/lipid table ff2.1 lost TM4 on Upside 2: 0.610, 9 of 12 seeds
unwound; 4.15). Only ff3.0 and the old ff3.1 have placement-fixed Upside 2 runs, and they differ from ff2.1 in every
retrained term.
* **Same in both engines:** KE/1.5kT 1.00 for protein and lipid on every rung; peptide C-N 1.312 +- 0.120 A (Upside 3)
  against 1.333 +- 0.118 A (Upside 2), 0.78-1.94 A in both, so the backbone's bond scatter is the hybrid's on either
  engine; TM1 30-48 dihedral-box helix fraction 0.986 / 0.956 at T 0.70; the other helices by DSSP alpha 0.76-0.93
  against 0.70-0.98.
* **TM4 is not:** DSSP alpha of 135-151 over block 1 is 0.44 at T 0.70 and 0.39 at T 0.795 (Upside 2's ff3.0 chain
  0.96 and 0.93; ff2.1's REMD 0.93 at T 0.80 and its 12 local seeds 0.94 at T 0.85, both placement-defective); its
  dihedral-box fraction 0.76 and 0.69 (Upside 2 1.00 and 0.99; ff2.1 0.98). Every TM4 residue is below 0.8 at T 0.70.
* **It goes in preparation stage 7.0**, the first in which the protein moves: TM4 DSSP alpha is 1.00 at the input and
  in every frame of stages 6.0-6.6 (the protein rigid: H-bond count 97.3, Rg 20.4 A throughout), then by tenth of
  stage 7.0's 401 frames (60,000 steps, about 540 time units at T 0.8647) 1.00 0.99 0.98 0.98 0.93 0.82 0.91 0.85
  0.59 0.41, while TM1 stays at 1.00.
* The log's "hbonds" column reads 97 for the rigid protein against Upside 2's 195, but the summed per-frame `hbond`
  score is comparable (174 against 195 at block 1's first frame): the column's definition differs, not the protein.
* **Two readings were open.** TM4 is the glycine-dense helix (3 in 17) and was the first casualty of the subsystem
  temperature mismatch (3.10) and of a mismatched coverage table (3.11). The audit changed the hybrid where it reaches
  TM4, each a fix: the side-chain/lipid table now charges NC3, ALA and GLY keep their rotamer rows and H-bond coverage
  (A83), and the rotamer 1-body is ff2.1's phi/psi-dependent term (A72). The loss was either a fix revealing ff2.1's
  own TM4 in a correct hybrid or a new defect.
* **ff3.1 loses TM4 too** (user 10-26; A234; `upside3_validation/glpg_ff31`, same preparation and seeds). ff3.1 is
  ff2.1 with only its rama.dat trained, so it tests the Ramachandran map alone. TM4 DSSP alpha 1.00 through stages
  6.0-6.6, then by tenth of stage 7.0 1.00 0.97 0.94 0.82 0.83 0.70 0.60 0.50 0.43 0.36 (last frame 0.59, mean 0.72,
  dihedral-box fraction 0.83 at the last frame), TM1 1.00. It falls below 0.8 from the sixth tenth (ff2.1 from the
  ninth) to the same end level, so the trained map alone does not hold TM4 and ff3.1's run stays at one block.
* **ff3.0 holds TM4** (user 10-26: "do submit ff3.0"; A235; `upside3_validation/glpg_ff30`, same preparation and
  seeds). ff3.0 is the one current force field with a fully fixed Upside 2 run (TM4 DSSP alpha 0.96 at T 0.70), so
  between it and Upside 3 only the engine and the preparation change; the rule set before the run was that a loss
  means an Upside 3 cause and a hold means ff2.1's loss is its own. By tenth of stage 7.0: 0.99 0.99 0.91 0.96 0.84
  0.84 0.79 0.81 0.82 0.82 (last frame 0.88, mean 0.88, dihedral-box fraction 0.94 at the last frame; TM1 0.95),
  settling near 0.82 by the fifth tenth where ff2.1 and ff3.1 keep falling. On the same engine, preparation and seeds
  the loss under ff2.1 is therefore ff2.1's own once side-chain placement is correct, not an Upside 3 defect; since
  ff3.1 differs from ff2.1 only in its map, what holds TM4 lies in ff3.0's other retrained terms. One seed and 540 time
  units at T 0.8647; REMD block 1 (remote_jobs.md §1) gives TM4 at T 0.70, against Upside 2's ff3.0 chains' 0.96.
* **The core (pure-protein ff2.1) against Upside 2** (10-27, Mac Studio; scripts and outputs
  `/project/trsosnic/yinhan/checks/u3_vs_u2_20261027/`):
  * Hamiltonian at identical coordinates: each tree's own writer builds the 16 Peng configs as `bench_run.py` does,
    each engine evaluates every node on 11 frames per protein (native, two noisy copies, 8 Upside 2 de novo frames);
    engines master, the control (the shared deployment's source, the binary the released-ff2.1 arms ran) and Upside 3.
    Against the control: the springs bitwise; hbond, rama and backbone pairs within 6e-3 E_up; only
    `bb_sigmoid_coupling_environment` (3.13) at E_up scale, median 0.20, max 1.05; `rotamer` at most 0.007 but on one
    frame of 176, a different BP fixed point (unchanged at tol 1e-6), where Upside 3 equals master (the control
    differs from master the same way, the 09-26 BP fix); dE(native) - dE(misfold) +0.09 +- 0.26 E_up over 16
    proteins. Open: on BBA's random start (one CA-CA pair at 2.3 A) `backbone_pairs` differs by 6.3 E_up and
    `rotamer` by 4.0, not diagnosed.
  * Peng, `upside3_validation/peng/score/compare_u3_u2.txt` plus BBL, NTL9 and proteinB de novo scored on the Mac
    (`pull_cold.py` + `score_local.py` reproduce `score_compare.py`, WWdomain de novo 0.43264): native, 16 proteins,
    mean TM 0.564 against 0.552 (+0.012 [-0.001, +0.028], Wilcoxon p 0.52), the engines agree; de novo, 9 of 16, 0.437
    against 0.461 (-0.025 [-0.051, -0.003], p 0.051, median -0.005), carried by top7 (-0.108) and alpha3d (-0.071), both
    still rising between run halves, one run per arm per engine.
  * Reweighting each engine's 1,500 scored de novo frames by the other's Hamiltonian at the same coordinates (ESS
    1,250-1,475) predicts Upside 3 - Upside 2 of -0.002 (top7 +0.000, alpha3d -0.009): the energy change does not
    explain the de novo gap. The seed replicates (remote_jobs.md §1 line 3) measure the run-to-run spread.

## 4. What the hybrid gets wrong, and what is still open

Every glpG hybrid result in 4.1-4.16 ran with the misplaced side-chain beads of 4.17 (fixed 10-17) and on the
pre-audit engine and side-chain/lipid table (3.14 A76, A83, A192); before 10-10 also on the unconverged table of 4.10.

### 4.1 Fold fidelity: a ~4.5 A helical-core deviation, cause unidentified (findings 100, 103)

The protein does not hold its tertiary helix packing. Helical-core CA-RMSD from the crystal **plateaus at
4.1-4.4 A** in POPE/POPG (3.48 -> 4.67 within the first segment, then flat across two more: equilibrium, not
drift).

| at T = 0.70, crystal-bonded amides | POPE/POPG |
|---|---|
| H-bond occupancy (median) | 0.844 |
| burial fails | 3.4% |
| both fail -> exposure | 3.34% |
| implied raw dG_open | 1.99 |
| helical-core CA-RMSD | 4.15 A |
| CA-Rg (crystal 20.43 A) | 20.77 |

There is no measured reference for how faithful this model can be, only the crystal. The DDM micelle column
(2.61 A core RMSD, 0.952 occupancy) is retired with DDM (2026-09-09); it was never a clean comparison (that
campaign died at block 2-3, and its Rg sat 1.8 A below the crystal, so it was compacted rather than more faithful).

**Ruled out by measurement (findings 100):**
* *Integrator.* The `avg_kinetic_energy/1.5kT` excess is +2-3% and dt-independent, far too small for a 15%
  unbonded population; covalent geometry is intact (worst C-N 1.78 A, 0 broken); H-bond loss goes 27% -> 16% from
  315 K to 245 K (~1.7x), what a ~2 kcal/mol opening free energy predicts. The dt scan of 4.2 carries the argument.
* *H-bond assignment.* On the crystal geometry Upside's H-bond score agrees with the DSSP electrostatic criterion
  to within 8% inside helices (DSSP 86.5%, Upside 78.4%; the 12 disagreements are marginal). ~16% of DSSP-helical
  amides are helix N-termini with no i-4 partner, so 86.5% is near the ceiling.
* *Lipid voids / bilayer prep.* Every backbone site has environment beads within 8 A; nearest-environment median
  5.2 A, max 7.3 A; coordination 8.86 within 8 A (9.59 in the TM belt).
* *Hydrophobic mismatch.* PO4-PO4 thickness 38.0 +- 0.1 A, acyl core 25.4 A against glpG's 28.2 A belt: -2.8 A.
* *The burial threshold.* Burial failure is nested inside H-bond failure (3.4% vs 3.34% of frames), so the cut
  value is not what binds.

**RD1: restoring the absent environment terms does not recover the fold (findings 103).** Three arms on the
CB-corrected placement, 750-780 k steps each, single system, T = 0.70, one identical starting configuration:

| arm | restored | helical-core CA-RMSD (A) | Rg (A) |
|---|---|---|---|
| `base` | nothing (CB fix only) | 4.61 +- 0.09 | 20.36 |
| `env` | protein self-burial | 4.71 +- 0.12 | 20.43 |
| `envfull` | all non-membrane Upside terms | 4.53 +- 0.08 | 20.16 |

`env - base` is +0.10 A and `envfull - base` -0.08 A, both inside the run-to-run scatter: no arm repaired anything,
and none over-compacted (Rg 20.2-20.4 against the crystal's 20.43). The CB correction did not improve fidelity
either (`base` 4.61 A against the earlier 4.15 A; not directly comparable, a single longer trajectory against a
16-replica REMD segment average).

Stock Upside in its implicit membrane is no more faithful (4.16): four ff2.1 REMD runs from the simulation PDB
(1.27 A from the crystal) reach helical-core CA-RMSD 3.5 A in the first quarter of 600,000 tu and 4.8-6.1 A in the
last at T 0.70 (three of the four still rising), and the fraction of the start's CA contacts (8 A) kept falls from
0.45-0.49 to 0.23-0.29. Those runs use standard Upside's `placement_scalar` and no dry-MARTINI coupling, so neither
is needed for a deviation of this size.

So the cause is **unidentified**. Ruled out: the integrator, the H-bond assignment, lipid packing and voids,
hydrophobic mismatch, the burial threshold, the CB placement, and the absent protein-protein environment terms.
The one remaining node-level difference from a standard configuration, the rotamer 1-body representation
(`placement_fixed_scalar`, fixed rotamer probabilities, against standard Upside's rama-dependent
`placement_scalar`), is not needed for the deviation (the implicit runs); whether it changes its size is not
tested. Do not add the environment terms to production.

### 4.2 Helical H-bond occupancy and non-cooperative opening

Measured on the raw per-donor scores (`get_protection_state.py --report-raw-data`):

| quantity | helix interior | loop |
|---|---|---|
| backbone H-bond occupancy, starting structure | **0.954** | 0.221 |
| backbone H-bond occupancy, trajectory | **0.681** | 0.226 |

The loss is inside the helices (~27 points; 24 of 108 helix-interior donors below 0.5 occupancy); the loops are
unchanged. It is not a thresholding artefact (the raw score is bimodal: 60% of helix-interior frames above 0.5, 29%
below 0.001, 4.6% within [0.001, 0.05] of the 0.01 criterion) and not the starting structure (0.954 at t = 0).

Cooperativity is `P(neighbour open | this open) / P(open)`, 1.0 for independent opening. It **must be
normalised**: the raw isolated-open fraction is higher in stock Upside (0.730) than in the hybrid (0.528) only
because stock opens 5% of the time and the hybrid 29%. Never compare a count of isolated events between two
systems with different event rates.

| arm | protein KE excess | helix occupancy | P(open) | cooperativity |
|---|---|---|---|---|
| **stock Upside, no environment** | **+1.19%** | **0.950** | **0.050** | **4.18x** |
| hybrid baseline (dt 0.009) | +5.32% | 0.711 | 0.289 | 1.80x |
| hybrid, uniform max gamma | +2.01% | 0.763 | 0.237 | 2.12x |

The hybrid opens the helical H-bond network **6x more often and about half as cooperatively** as stock Upside on
the identical construct. Real local unfolding is cooperative; 59% of the hybrid's openings are one amide alone with
both neighbours still bonded, which produces log-amplified residue-to-residue scatter in a dG profile.

**The thermostat defect is real, is finite-dt error, and is not the cause.** `thermostat.cpp:31` makes the global
OU thermostat skip every atom with gamma > 0, so each protein backbone atom's only thermostat is its own Langevin
friction, proportional to lipid-contact count (0 to 8.80, i.e. 0 to 50x the bare value), while lipids are uniform at
0.169: buried helix cores drain slowest and run hottest. The g-JF propagator itself is correct. A dt scan at
**matched physical time** (270 t_up, steps scaled as 1/dt):

| dt | steps | protein excess | lipid excess | helix occupancy | cooperativity |
|---|---|---|---|---|---|
| 0.00900 | 30 000 | +5.32% | +1.02% | 0.711 | 1.80x |
| 0.00450 | 60 000 | +2.39% | -0.88% | 0.765 | 1.71x |
| 0.00225 | 120 000 | **+1.07%** | **+0.02%** | 0.758 | **1.53x** |

The excess falls by a factor of 2.23 per halving and the lipids land on target, so the heat source is the
discretisation. **dt does not recover the fold**: occupancy plateaus at ~0.76 against stock's 0.950 and
cooperativity does not improve. The excess opening and the non-cooperativity come from the hybrid Hamiltonian; no
timestep, friction or thermostat change removes them, and no analysis-side correction factor can either: the defect
is missing cooperative structure in the sampled ensemble, not a mislabelled temperature.

Limits of these numbers: 270 t_up is short and occupancy is still decaying from the seed's 0.954, so ~0.76 is an
upper bound; the stock arm has no membrane (glpG's TM helices in vacuum, Rg 16.6 A), which over-stabilises
intramolecular H-bonds, so 0.950 and the gap are upper bounds too. The clean comparison needs stock Upside with its
implicit-membrane terms, which the hybrid topology does not carry (2.2).

PS hides most of this: 79% of broken-H-bond frames are still called protected because burial exceeds
`criterion3 = 5.0`, so helix-interior PS reads 0.936 while H-bond occupancy is 0.681. **A composite indicator
(`A OR B`) on a system where `B` is nearly always true stops measuring `A`.** For a compact membrane protein nearly
every backbone amide is buried by the protein's own atoms, so quote the H-bond occupancy profile next to the dG
profile.

### 4.3 Lipid packing against the protein is a seed property

The lipid-shielded fraction of amide N (>= 1 tail bead within 7.00 A, 5.3) differs between two POPE/POPG datasets
by almost 2x, and it traces to the seeds: the local 16-replica seed's last frame reads **0.857** (run flat at 0.87
from the first block), the cluster 48-replica seed's **0.433** (run climbing 0.32 -> 0.50 over 14 chunks, plateauing
near 0.50). Same protein, 4949 atoms, 279 lipids, box to 0.5 A, criterion and code. It is not insertion and not a
broken protein: 87.1% of CA sit within 20 A of the bilayer midplane in the cluster run, constant across all 16
chunks; Rg 20.4 A (the crystal value); consecutive CA-CA 3.79-3.90 A with no exceptions; the 12 CA at 30-60 A from
the midplane are the N-terminal tail in vacuum, as in the local run. What is missing is lipid packed against the
surface. An HDX profile from the cluster data therefore shows far fewer +inf amides for reasons that are neither the
force field nor the analysis, and the two datasets must not be compared as if they differed only in ladder size.

**Check the environment's contact with the protein, not just the protein's position in the environment.**
Insertion depth and Rg were correct and constant here while half the protein-lipid contact was missing.

### 4.4 glpG TM4 before the temperature fix: what was ruled out (2026-09-07 to 09-10)

Before the protein/lipid temperature mismatch was found (3.10), TM4 lost helix in every glpG variant: in the ff_2.1
production REMD (T 0.70 rung) every variant went from TM4 = 1.000 to 0.35-0.60 (wild type 0.432) while TM1 mostly
held (0.77-0.94; 79ALA 0.56). TM4 (134-151) is byte-identical in all four variants. With `inner_steps = 4` the cold
end is healthier than the seed (3.10c). Kept so it is not re-tested:

* **Measurement windows.** The crystal seed has 131-133 and 152 non-helical, so TM4 scored on 131-152 is capped at
  18/22 = 0.818 against a > 0.8 criterion, and TM1 on 29-49 at 0.952. Use 134-151 and 30-48. The GLY49/GLY133 phi
  criterion is invalid: both are helix caps at phi_std +94 and +142 in the seed.
* **Compare trends, not means.** Arms that have not equilibrated cannot be ranked by trajectory means: the ff_2.1
  control had bottomed out by its third quintile while the trained arms were still falling (second-half slopes
  -0.022 against -0.055 and -0.135 per 1000 time units), which made "TM4 is fixed" a convergence error. Report
  quintiles and the second-half slope. Single-temperature MD at T 0.70 is also harsher than the REMD rung the > 0.8
  criterion came from (control 0.441 against 0.645).
* **What the loss was.** Not thermal melting (flat at 0.36-0.57 across the ladder) and reversible (over 13,117
  frames of the WT run 47.8% above 0.8 while drifting down). Local: 134-142, the buried midplane half (139-141 from
  1.000 to 0.33-0.35; i->i+4 O...N 4.3-6.4 A at 139-142 against ~3.0 formed), while 143-146 and 149 held at 1.000.
* **Not lipid.** TM4 is the least lipid-exposed helix (2.0 backbone beads with a lipid neighbour inside 12 A against
  TM1's 10.7; zero dry-MARTINI energy between 134-139 and all 3627 lipid beads), so every protein-lipid explanation
  fails for it, including the MARTINI BB typing of 131-137 (re-typing all of it to `N0` changes the total by -0.3
  kJ/mol). The raw lipid-gain/helix-loss cross-correlation (r 0.64, lag 99) is a shared-trend artefact (first
  differences r 0.21-0.26); lipid between TM helices is acyl tail (enrichment 1.20x, headgroups 0.39x).
* **Not glycine density** (r = +0.23, wrong sign; TM3 is 17.9% glycine at 0.928), nor unpaired buried charges
  (ARG148-GLU150 4.1 A, ARG151-ASP152 3.4 A).
* **The backbone-environment term is absent from the hybrid and does not fix TM4.** The hybrid builder never passes
  `--environment-potential`, which gates every environment node in `upside_config.py`, so glpG has no
  `bb_sigmoid_coupling_environment`, `environment_coverage_*` or `hb_environment_coverage_*` nodes (2.2). Restoring
  them moved TM4 by +0.020 against a seed scatter of +-0.10, because TM4's backbone sits at burial 13.7 where the
  term saturates above 4. `hbond_coverage` is not a desolvation measure; its output is per side-chain bead, a 1-body
  cost handed to the rotamer solver.
* **A MARTINI-typed apolar coverage channel made things worse** (TM4a 0.807 -> 0.464, TM1 0.911 -> 0.810). The
  coverage is non-negative with a positive weight, so it repels backbone from acyl chains; a lipid-aware channel on
  these nodes can only penalise lipid proximity. Not deployed.
* **Engine facts.** Node types resolve by prefix (`deriv_engine.cpp:598`), so a group named `hbbb_coverage_lipid`
  instantiates `hbbb_coverage`. `rotamer` takes a variable-length `prob_nodes` list and uses each output directly as
  a 1-body energy, each with `n_elem` equal to the side-chain bead count (747). `hbbb_coverage` has `n_dim2 = 3`, so
  its second group can be bare positions; the trained `hbond_coverage` needs `n_dim2 = 6`. A lipid-aware trained
  coverage is an engine change in `src/environment.cpp`.

**Open hypothesis, untested.** With `exclude_intra_protein_martini = 1`, all helix-helix packing comes from the
Upside core force field, trained on soluble proteins only, while protein-lipid attraction is full-strength
dry-MARTINI, so lipid competes for the same hydrophobic surfaces against a term never balanced against it. The
principled test is ConDiv training with membrane proteins in the set; scaling SC-env or BB-env down is forbidden.

### 4.5 Lessons from the FF1-form ff3.0 retrain (2026-09-07 to 09-09; trainer retired 2026-09-24)

The first ff3.0 was trained by a port that turned out to be FF1's workflow (9t-9v); its force fields are retired.
What it showed about ConDiv runs in general:

* **The parameters random-walk while the objective plateaus.** Drift of `pair_interaction` grew as sqrt(lag)
  (rel_rms / sqrt(lag) 0.013-0.016 over a 16-fold range of lag) with no fixed point, while the median restrained
  RMSD was flat over the last half. Training longer only diffuses further; a downstream observable has to choose the
  step. The side-chain random walk of 1.17 is the same effect.
* **Two runs that share history are not replicates.** A run resumed from another's step-274 checkpoint diverged
  from it at the single-run diffusion rate (rel_rms 0.26 by step 496); extrapolated, two independent trainings from
  ff_2.1 would differ by ~40%. A reproducibility test needs runs branched at step 0, and force fields from different
  steps must never be paired.
* **Retraining reproduces how far it moves, not where.** Two lineages from a common step 269 drifted the same
  distance from ff_2.1 to within 0.1-1.7%, in directions ~23 deg apart (between-run rms 33-40% of the drift).
* **The first update after a resume is outsized** (rel_rms 0.076 after one step, against 0.018 per sqrt(step)),
  most likely because the optimiser's momentum state was not carried across.
* **Most of the movement happens in the first ~60 steps** (cumulative rms 0.49 at step 29, 0.64 at 59, then ~0.07
  per 30 steps without an asymptote), so a step target needs a measured justification.
* **Trained tables and the coverage nodes only work together.** On glpG at T 0.70 (3 seeds), TM4 helix fraction was
  0.441 for ff_2.1 without coverage nodes, 0.562 with the nodes and old tables, 0.588 with new tables and no nodes,
  and 0.782 with both: a pair table co-trained with coverage is correct only beside the partners it was optimised
  against. 135 further steps (269 -> 404) left TM4 flat (0.709). At an unrelaxed seed the trained force field's
  force tail was heavier (p99.9 97 against 49 E_up/A), part of it from the coverage nodes themselves.
* A checkpoint of that trainer could be rebuilt from an extracted force field (`pack_param`, an L-BFGS-B refit,
  reproduced trained tables to ~1e-16; Adam state is lost, a short transient with beta2 = 0.96); the acceptance test
  is that `extract_ff.py` on the rebuilt checkpoint reproduces the starting files.
* **Upside overwrites `/output` rather than appending**, so a fresh run on a production seed replaces the seed's
  frames. Tell them apart by the time spacing, not the frame count.

### 4.6 Which glpG TM helices lose structure, and from what (2026-10-09)

On the old side-chain/lipid table, whose defect 4.13 and 4.15 show caused the TM1 and TM4 losses. TM4-test sets of
1.24 (12 seeds per set, last block, DSSP; any-helix = H, G or I; ff3.0 is d9_03, 1.26). Data and tables
(`sets_bb*.npz`, `sets_full.npz`, `tm_attribution.txt`, `tm1_tm2_hbonds.txt`) in the MacBook Pro's
`scratchpad/glpg_ff21_vs_ff30/`; the scripts and outputs of 4.6-4.9 copied to midway2's
`/project/trsosnic/yinhan/checks/tm1_hbmem_20261009/` (its README lists what was copied).
- **TM1 lost structure in every set** (any-helix 0.75-0.85), at its middle 32-38, which sits at the bilayer centre
  (CA z -7 to +3 A): ff2.1 unwound it (about 39% of i->i+4 H-bonds broken), ff3.0 made more of it a pi-bulge (20%
  i->i+5). TM1 has no glycine.
- **TM2 was lower in ff3.0 than ff2.1** (any-helix 0.79 against 0.91, p 0.04), at 82-85, a loss that comes with
  training of the non-glycine terms (gdepth_start 0.86 against the control c9_01's 0.79, p 0.01), and at 91-96.
  **TM2's start has no gap**: at 91-96 every i->i+4 O...N is 2.8-3.5 A and phi/psi are helical; DSSP calls 92-96 pi
  from one bifurcated O92...N97 contact at 3.44 A.
- TM5 was within noise of ff2.1 (0.88 against 0.92, p 0.34).

### 4.7 Nothing in the hybrid charges for an unpaired backbone in the acyl core (2026-10-09)

`tm1_diag.txt`, `tm1_sets.txt`, `tm1_hbmem.txt`, in the MacBook Pro's `scratchpad/glpg_ff21_vs_ff30/`. Measured on the
old table, whose kicks turned out to open TM1 (4.10, 4.13); the missing cost below still holds (4.12).
- **Lipid gain is not the driver.** Time-matched (last block, 12 seeds), TM1's middle has the same environment intact
  or open: ff2.1 about 18 tail beads within 7 A either way, BB-env energy -43.7 against -44.3 E_up (hb4 against
  BB-env r +0.06); ff3.0 -40.6 against -44.0. No headgroup or ion comes near.
- **What opening costs in the hybrid:** only Upside's H-bond energy, E_alpha -1.96 E_up per bond (about +8.4 E_up for
  the 4.3 bonds the middle loses). MARTINI's BB typing is fixed at the starting secondary structure (TM1: 17 N0, 3
  C5), so an opened backbone keeps its helix-typed lipid attraction (N0-C1 well -1.18 E_up; MARTINI's coil type P5
  has -0.17).
- **Upside's own membrane model charges for exactly this.** `hb_membrane_potential`
  (`src/membrane_potential.cpp:488`, `parameters/ff_2.1/membrane.h5` hb_energy) gives each backbone donor and acceptor
  (1-p)^2 f_unpaired(z) + (1-(1-p)^2) f_paired(z). f_unpaired - f_paired is about +2.2 (donor) and +1.8 (acceptor)
  E_up within ~5 A of the midplane and slightly negative in water. Its state-dependent part on the stored frames,
  with protein_hbond's p (0.74 intact, 0.35 open in ff2.1), would add +7.0 / +8.4 / +10.3 E_up to opening 32-38 at
  half-thickness 12.7 / 14.0 / 15.9 A (ff3.0 +4.4 / +5.3 / +6.4). The hybrid has no such term (2.1): it charges about
  half of what Upside's membrane model says opening TM1's middle costs.
- Making each BB bead's MARTINI type follow its H-bond state (N0 paired, P5 unpaired) would charge about +38 E_up
  for the opening on the stored frames, far more than the membrane term (plan.md Phase 13, option B).
- The bilayer centre sits within 1 A of z = 0 and drifts 1-2 A over 4000 tu (PO4 mid-plane).

### 4.8 Upside's membrane H-bond term on ff2.1 (old table, 12 seeds, 2026-10-09)

- **The term** (option A, chosen by the user 10-09; plan.md Phase 13). E = sum over backbone donors and acceptors of
  (1-p)^2 [f_unpaired(z) - f_paired(z)], p from the hybrid's `protein_hbond`, f from `parameters/ff_2.1/membrane.h5`
  hb_energy unchanged; the paired baseline is left out because MARTINI's helix-typed BB already describes a paired
  backbone against lipid. It is the engine's own `hb_membrane_potential` node, written into copies of the patched
  inputs by `scratchpad/glpg_tm1_hbmem/add_hbmem.py` (coeff = [hb_energy[t,0] - hb_energy[t,1], 0], use_curvature 0,
  bilayer centre at z = 0). Half-thickness 15.9 A: the mid-plane of GL1 (16.42 A) and GL2 (15.40 A) over the second
  half of the local ff2.1 run, the hydrocarbon-core boundary Upside's membrane thickness denotes. Energy and forces
  verified on stored frames to ~1e-3 (`verify_hbmem.txt`).
- **The test** (`scratchpad/glpg_tm1_hbmem/analyze_hbmem.txt`, MacBook Pro; positions `runs_hbmem_ff21_pos.npz`,
  copied to the cluster with the logs): against ff21_released, last-block DSSP any-helix (alpha), two-sided
  Mann-Whitney, 12 seeds each. TM1 0.79 -> 0.90 (p 0.03), TM4 0.71 -> 0.82 (p 0.02), TM5 0.92 -> 0.97 (p < 0.01),
  TM2/TM3/TM6 +0.01 to +0.04 (n.s.); TM4 primary 0.609 -> 0.788 (p 0.03).
- It failed the health check: protein potential jumps of 3,800-17,900 E_up in 5 of 12 seeds, two of them TM4
  backbone tears over 137-148 (C-N up to 9.5 A). **The tears are the side-chain/lipid table's kicks, not the term**
  (4.10): probe replays reproduce the runs bit for bit, each tear starts with one evaluation of 77,407 (s1) and
  153,463 (s2) E_up/A on MET142's CA, all of it the rotamer node's and gone with the table zeroed, and
  `hb_membrane_potential` is never the largest term (`kick_ff30/attribute_events_ff21_hbmem.txt`).
- The term on ff3.0 is 4.11 (old table) and 4.14 (converged table).

### 4.9 Transient backbone excursions sit where the helices fail (2026-10-09)

`scratchpad/glpg_tm1_hbmem/cn_events_vs_tm1.txt` and `kick/` (MacBook Pro); ff21_released and d9_03 on the old table,
24 seeds. Their cause is the table's single-evaluation kicks (4.10).
- 61 of 9,624 stored frames hold a peptide C-N above 2 A (up to 6.0 A; normal frames' maximum 1.6-1.8), each for one
  stored frame, clustered in TM1 33-36, TM4 137-146 and the 74-82 loop; TM1's middle i->i+4 H-bonds average 5.10 over
  the 5 frames before its 16 excursions and 3.83 after. Their potential jumps (e.g. +1477 E_up) are under the 3000
  threshold of the TM4 test's jump scan (1.24), so these events are not in its counts.
- Not a MARTINI core collision (closest BB-environment approach median 5.09 A against blow-up 2.43 A, 2.4). A replay
  of ff2.1 seed 1 with frames every 10 steps puts the onset inside one 10-step interval (MET34's backbone moves 3.6 A
  against 0.5-0.9 A for every other atom), with no MARTINI pair spike (`UPSIDE_MARTINI_PAIR_DIAG`) and no precursor in
  any term at the stored frames (`kick/force_by_node.py`).

### 4.10 The TM1-middle kicks come from the side-chain/lipid table's sampling (local ff3.0, 2026-10-09)

Mac Studio, `scratchpad/glpg_tm1/kick_ff30/` (scripts and outputs; plan.md Phase 13 step 2). The 600 x 12 table
built here was itself replaced in Upside 3 (3.14 A76, A192).
- **The probe.** A temporary read-only build (`scratchpad/glpg_tm1/src_probe`, `obj_probe/`; never in `src/`)
  records, after every derivative evaluation, the largest |dE/dx| over protein atoms 0-839 and keeps the last 120
  evaluations' positions, written out above 1000 E_up/A. It reproduces the stored d9_03 logs of s7, s8, s10 and s11
  frame for frame. Over s8 to t 300 the median is 85 E_up/A and the 99.9th percentile 550.
- **The events.** 16 above 1000 E_up/A by t 200-500 (four per seed, the probe's cap), up to 33,260 E_up/A, each one
  evaluation between normal ones (s8: 9.6, 8,305, 5.2 on CA34), on CA34 and CA35 (TM1's middle; both MET) except two
  on CA83 and one on CA169. One evaluation at 8,305 E_up/A gives a CA about 18 A/tu (thermal 1.6), the size of 4.9's
  excursions. ff2.1 + term's TM4 tears (4.8) start the same way on MET142, GLY143's neighbour; whether its kicks feed
  GLY143's flips in the TM4 test is not measured.
- **The term.** The engine evaluated offline at the written positions reproduces the recorded dE/dx to 5e-4 E_up/A
  (no hidden state). In 15 of 16 events the whole spike is the rotamer node's; with `martini_sc_table_1body`'s tables
  zeroed it falls to 0.7-36 E_up/A, and with BP converged to tol 1e-6 it stays (8,300). The sixteenth is Spring_omega
  pulling back a backbone an earlier kick tore.
- **The mechanism** (s8, evaluation 4688). MET34's rotamer 2 (weight 0.59) against a lipid tail bead (C1) at 4.90 A
  and angular coordinate 0.167, a grid node: value -0.61 E_up, dV/dcos -12,593, 1,504 E_up/A on the CB point; one
  evaluation later 0.83. The table slice has a ridge one node wide at cos theta = 0, 2,091 E_up at 4.9 A and 795,448
  at 3.5 A against 4.9 and -0.6 at the neighbouring nodes, and bilinear interpolation turns it into a wall across two
  cells.
- **The ridge is sampling, not the potential.** `parameters/ff_2.1/martini.h5` `sc_table` (copied unchanged into
  every hybrid input) was built by `py/martini_build_tables.py` `_build_sc_table_group` with 24 target directions
  over the whole sphere and one side-chain bead frame (`sidechain_bead_frame_count` 1), so each of the 13 cos nodes
  rests on about two directions, not the azimuthal average the table states (tempered Boltzmann at T 25). The
  builder's own `_run_sc_task` reproduces the stored slice at those settings (to 0.03 E_up); with 12 frames it gives
  8.67 E_up at 4.9 A (stored 2,091) and 11.87 at 4.0 A (stored 85,257), and 600 directions x 12 frames and 2400 x 24
  agree to 0.01 (`sc_table_convergence_MET_C1_rot2.txt`). Over the whole table the largest value is 1.1e9 E_up and
  145,841 steps between neighbouring cos nodes exceed 100 E_up.
- **A vectorised rebuild** (`sc_table_vectorized.py`) uses the builder's pair energy (equal to
  `_compute_pair_energy_and_gradient` to 7e-15 relative), anchor, binning and average; at 24 directions x 1 frame it
  reproduces all 5.1 M stored cells to 1.2e-7 relative, so only the sampling density changes.
- **The converged table removes the kicks** (`old_vs_fixed_traces.txt`). The whole table rebuilt at 600 directions x
  12 frames (finite; largest value 1.4e7 against 1.1e9; angular steps above 100 E_up at r >= 4 A 3,230 against
  25,912, those left near cos theta = +-1, mostly ARG, LYS, ASP and PRO), written into a copy of the d9_03 input with
  only `rotamer_full_energy_eup` changed. Over the replays' durations (s10 to t 200, s7 and s8 to 300, s11 to 500) the
  largest |dE/dx| on a protein atom is 187-236 E_up/A against 32,572-65,012, no evaluation exceeds 500 against
  152-1,343, the median falls from 134-212 to 66-69, and the largest displacement of a protein atom between
  evaluations from 0.16-0.33 A to 0.012 A.
- Every hybrid run before 10-10 read the old table, the TM4-test sets and the glpG validation chains included. It was
  replaced in production on 10-10 at the user's word (4.15).

### 4.11 The membrane H-bond term on ff3.0: TM1 leans more helical, no tears (2026-10-10)

Option A (4.8) on d9_03 with the old table, 12 seeds (`scratchpad/glpg_tm1/tm1_compare_ff30_hbmem_vs_d9_03.txt`).
TM1, the readout fixed before the run (last-block DSSP alpha of 29-48, two-sided Mann-Whitney): 0.818 against
d9_03's 0.685, +0.133 [-0.000, +0.259], p 0.073, not resolved; TM4 0.829 against 0.726 (p 0.49). No total-potential
jump above 3000. The table's kicks are still in these runs; 4.14 is the term on the fixed table.

### 4.12 Is TM1 fixable without training? What the literature says (2026-10-10)

The user's question (10-10). Literature agent, citations checked against abstracts or full text; numbers converted at
1 E_up = 0.697 kcal/mol.
- **Experiment does not show TM1 opening.** GlpG's crystal structures (Wang, Zhang, Ha 2006 Nature; Wu et al. 2006
  NSMB; Ben-Shem et al. 2007 PNAS) and solid-state NMR in liposomes (Shi et al. 2019 JACS) put the mobility in TM5,
  loop L4 and the gating loops. TM1-TM2 carry the folding nucleus (Paslawski et al. 2015 PNAS), folding runs N to C
  in hairpins (Choi et al. 2019 Science), the N-subdomain TM1-TM3 is the more stable (Guo et al. 2016 Nat Chem Biol;
  Gaffney et al. 2022 PNAS), and TM1 and TM2 are GlpG's most hydrophobic helices (Lu, Schafer, Wolynes 2018 Nat
  Commun). The one contrary hint: HDX-MS in native nanodiscs (Reading et al. 2017 Angew Chem) found a TM1 peptide
  (94-108, full-length numbering, not mapped to ours) relatively unprotected; it cannot say where along TM1.
- **The physical cost the hybrid lacks.** Water to alkane, a helical peptide unit costs 2.1 kcal/mol H-bonded and 6.4
  not (Ben-Tal et al. 1996 Biophys J); an amide H-bond is worth about -0.5 in water and -4.4 kcal/mol in alkane in
  helix geometry (Ben-Tal et al. 1997 J Phys Chem B); White & Wimley 1999 (Annu Rev Biophys) give 4-5 kcal/mol per
  broken H-bond in the core, against ~0.4 at the interface (Ladokhin & White 1999; Almeida et al. 2012). That is
  6.3-9 E_up per broken H-bond at the midplane; the hybrid charges 1.96 E_up (Upside's soluble-calibrated E_alpha),
  the water value, and MARTINI's helix-typed BB adds nothing (4.7). The undercharge favours each opened residue by
  ~e^4-e^6, and it sits in the coupling, which is why every force field shows it (4.6).
- **The two options' magnitude.** Option A (4.8, Upside's membrane term, a knowledge-based depth-dependent term from
  native TM statistics, Wang et al. 2018 Biophys J, about 2 RT per unpaired group) brings a broken H-bond at the
  midplane to about 6 E_up with the H-bond itself, the low end of the physical range. Option B (MARTINI's own BB
  typing by H-bond state, +38 E_up for TM1's middle, 4.7) is above it. Using both charges the same cost twice; only
  one may be used.
- **Coarse-grained precedent.** MARTINI keeps secondary structure as an input and does its H-bond bookkeeping by BB
  polarity (Monticelli et al. 2008 JCTC; dry MARTINI, Arnarez et al. 2015 JCTC); PLUM, which lets TM helices fold,
  keeps WALP helical but underestimates insertion by ~15 kcal/mol, the size of a missing backbone desolvation (Bereau
  et al. 2015 J Chem Phys; our reading).
- What would prove the mechanism in our runs: opened TM1 residues within ~10 A of the midplane with unpartnered NH
  and C=O (4.6-4.7 put 32-38 at CA z -7 to +3 A), the same residues closed near the interface, and TM1 holding once
  the cost is added on a defect-free table (plan.md Phase 13).

### 4.13 The converged side-chain/lipid table holds TM1, with no training and no new term (2026-10-10)

ff3.0 (d9_03) with only `rotamer_full_energy_eup` replaced by the 600 x 12 table (4.10), 12 seeds with
`run_glpg.sh`'s flags on the probe build, against the 12 local d9_03 seeds on the old table
(`scratchpad/glpg_tm1/tm1_compare_sc600x12_vs_d9_03.txt`, `runs_sc600x12/`).
- **TM1, the readout fixed before the run:** last-block DSSP alpha 0.931 against 0.685, +0.246 [+0.146, +0.341],
  two-sided Mann-Whitney p 0.002, **resolved**. Every seed keeps TM1 at 0.81 or more (old table: four seeds at
  0.43-0.60). Helix of any kind 0.960 against 0.848 (p 0.043); the middle's i->i+4 H-bonds 6.01 of 7 against 4.41 (p
  0.017).
- **TM4** 0.924 against 0.726 (+0.197 [+0.043, +0.376], p 0.10), not resolved.
- **Health.** No evaluation above 1000 E_up/A in 12 x 741,543 (largest 769); one frame in 4,812 with a C-N above 2 A
  against 33; no total-potential jump above 3000; KE/1.5kT 0.999-1.005.
- So TM1's loss under ff3.0 was the table defect, not missing physics and not the force field: the table is the
  potential the builder defines, computed with enough azimuthal samples, so this changes no parameter. The cost of
  unpaired backbone groups in the core (4.12) is still absent from the hybrid; 4.14 shows what it adds on top.
- Runs on the old table (the round-4 TM4 sets, the glpG validation chains before 10-10): their comparisons with each
  other stand, their absolute TM1 and TM4 values do not.

### 4.14 The membrane H-bond term on the converged table: TM1 and TM4 both held (2026-10-10)

Option A (4.8) on top of 4.13: d9_03 with the 600 x 12 table and the term
(`fixed/79HIS_d9_03_sc600x12_hbmem.up`), 12 seeds with `run_glpg.sh`'s flags on the probe build. Tables
`scratchpad/glpg_tm1/tm1_compare_sc600x12_hbmem_vs_d9_03.txt` and `_vs_sc600x12.txt`.
- **Against d9_03 on the old table:** TM1 last-block DSSP alpha 0.964 against 0.685 (+0.279 [+0.183, +0.373], p <
  0.001), helix of any kind 0.979 against 0.848 (p 0.003), middle i->i+4 H-bonds 6.89 of 7 against 4.41 (p 0.001).
  **TM4 0.959 against 0.726 (+0.232 [+0.085, +0.405], p 0.009), resolved**, the first resolved TM4 difference of any
  local set. Every seed keeps TM1 at 0.85 or more and TM4 at 0.84 or more.
- **Against the table fix alone (4.13):** TM1 0.964 against 0.931 (+0.033 [-0.003, +0.069], p 0.049), just resolved,
  with the interval reaching below zero; middle H-bonds 6.89 against 6.01 (p 0.001); TM4 0.959 against 0.924 (p
  0.21), not resolved. The table fix carries most of the gain; the term adds a small one to TM1's middle.
- **Health.** No total-potential jump above 3000, no frame with a C-N above 2 A, KE/1.5kT 0.998-1.005. Three
  evaluations above 1000 E_up/A in 12 x 741,543 (1,567 on ASN185's CA in s7, 1,283 on SER119's CA in s9), outside
  TM1 and TM4, are protein side-chain interactions: they stay with the side-chain/lipid table zeroed
  (`kick_ff30/attribute_events_sc600x12_hbmem.txt`).
- **Not adopted** (user, 10-10 11:20). With the table fixed its gain (TM1 +0.033) is below the spread between two
  table-only sets (0.931 and 0.882, 4.15); it is knowledge-based and imposes the bilayer's centre and thickness on the
  protein from outside, against an explicit bilayer that drifts; it may double-count MARTINI's backbone-lipid cost;
  and it would void every TM4 reference again. Reopen it only if a fixed-table run shows a TM helix opening in the
  core with unpaired backbone groups (4.12's test); a fix then belongs in the backbone-environment interaction, not a
  z profile. Never with option B (4.12).

### 4.15 The converged side-chain/lipid table in production (2026-10-10)

The user (10-10): apply 4.13's fix to the force-field files so future simulations run with it, everywhere: the
queued glpG chains re-patched and the TM4 references rerun on it (plan.md Phase 15). Records
`scratchpad/glpg_tm1/sc_table_prod/` and midway2's `checks/sc_table_20261010/`. In the repo the table was rebuilt
again for Upside 3 (3.14 A76, A192).
- **The builder** (`py/martini_build_tables.py`) took the sampling as arguments (600 directions x 12 bead frames,
  stored as `@target_direction_count`, `@sidechain_bead_frame_count`); at 24 x 1 it reproduced the old table but for
  9 of 5.1 M cells one float32 step off (`gate_old24.txt`), at 600 x 12 the TM1-tested table to float32 rounding,
  md5-identical between builds (`gate_prod.txt`). `$P/py/martini_build_tables.py` keeps the old sampling (the
  training tree is not overwritten during a campaign); never rebuild martini.h5 there.
- **`parameters/ff_2.1/martini.h5`** became md5 `2e612261...` (was `87e8cfbf...`, kept as
  `.bak_pre_sc600x12_20261010`), locally and in both cluster trees (`deploy.sh`, `deploy.log`); in the repo it is
  `martini.h5.bak_20261020` since A76.
- **Patched seeds.** `patch_glpg.py` (md5 `97ca3294...`, was `70589119...`; in `$P/training/` and
  `tm4_local/scripts/`) also copies the table into `martini_sc_table_1body` (`--martini`, default the file the
  preparation reads); with the old martini.h5 it reproduces the existing d9_03 TM4 input in all 150 datasets. The
  queued ff_3.0_gdepth and ff_3.0_bio_si chains were re-patched before they started: in each of the 8 seeds only the
  four table arrays changed (`check_repatch.py`, `repatch.log`; old seeds kept as `.bak_pre_sc600x12_20261010`).

**TM4 references on the fixed table.** Local `<tag>_sc` sets, 12 seeds each, each input differing from its old-table
input only in the table; old-table sets keep their tags and compare only among themselves. TM4 DSSP alpha, difference
with bootstrap interval, two-sided Mann-Whitney p; tables `tm4_compare_cov_<set>.txt` (against its start) and
`tm4_compare_cov_<set>_vs_<ref>.txt`. Every `_sc` set: no jump, KE/1.5kT 0.995-1.008, TM1 0.965-0.988.

| set (force field, step) | TM4 | against | difference [interval], p | unwound / flipped of 12 |
|---|---|---|---|---|
| ff21_released_sc (ff2.1) | 0.944 | ff21_released, old table, 0.610 | +0.334 [+0.160, +0.527], 0.002 | 0 / 0 against 9 / 5 |
| gdepth_start_sc | 0.947 | ff21_released_sc | p 0.82 | 1 / 0 |
| d9_01_sc (ff30_gdepth_dt009, step 37) | 0.925 | gdepth_start_sc | -0.023 [-0.074, +0.031], 0.45 | 0 / 0 |
| d9_02_sc (step 56) | 0.916 | gdepth_start_sc | -0.032 [-0.082, +0.021], 0.17 | 2 / 1 |
| d9_03_sc (ff_3.0_gdepth) | 0.884 | gdepth_start_sc | -0.063 [-0.124, -0.001], 0.13 | 2 / 1 |
| | | ff21_released_sc | -0.059 [-0.121, +0.001], 0.15 | |

- **On the fixed table ff2.1 itself holds TM4**, against its old-table set: dihedral TM4 0.985 against 0.827, TM1
  0.973 against 0.848; the old set's 135-151 profile (0.40-0.80) is 0.99-1.00 to 148, then 0.89, 0.63, 0.58. So the
  TM4 loss the round-4 training was asked to fix was, in this test, the table defect.
- **Every start reads alike and no candidate improves on ff2.1**; none is resolved from its start. The other lines,
  since retired, read: ff_3.0_bio_si (bs_03_sc) 0.951 and ff_3.0_bio (b9_03_sc) 0.950 against their start
  ff21_bioT1_6_sc's 0.939 (p 0.77, 0.36) and ff2.1 (p 0.91, 0.67); ff30_gdepth_si steps 37 and 56 (ds_01_sc,
  ds_02_sc) 0.894 and 0.889 (p 0.19, 0.47 against gdepth_start_sc); ff31_depth_all and ff31_depth_data step 18
  (fa_00_sc, fd_00_sc) 0.926 and 0.914 (p 0.60, 0.23 against ff2.1). b9_03_sc against d9_03_sc, +0.065 [-0.003,
  +0.131], p just under 0.05, is the only candidate pair at the edge of resolution.
- d9_03_sc reproduces 4.13's table-only set within noise (TM1 0.882 against 0.931, p 0.19; TM4 0.884 against 0.924, p
  0.24) but is not frame-identical: the production table differs from the tested one by float32 rounding and the
  test ran on the probe build (same initial energy, -24938.38).
- **Across ff30_gdepth_dt009's epochs TM4's C-terminal end frays step by step** (147-151, last block): start 1.00,
  1.00, 0.95, 0.73, 0.65; d9_01_sc 0.99, 0.99, 0.88, 0.68, 0.62; d9_02_sc 0.99, 0.92, 0.79, 0.47, 0.52; d9_03_sc
  0.92, 0.83, 0.72, 0.46, 0.42. No single step is resolved; the direction is the same at every epoch end. d9_03_sc's
  dihedral TM4 is 0.952 against its start's 0.987 (one-sided lower p 0.04). On the old table d9_03 led its start
  (0.726 against 0.657, p 0.64); both pairs are unresolved, so the sign change is not evidence.

### 4.16 TM5 stays helical but leaves TM6 for TM2, at T 0.70 and 0.80 alike, in both membrane models (2026-10-13)

The PI's question (for an H145A intermediate model): does TM5 stay folded and docked against TM6, and how different
are the T 0.70 and T 0.80 contact maps. Definitions: native CA-CA < 8 A, |i-j| >= 5 (`diag_contacts.py`) and CA-CA <
10 A (`tm1_diag.py` bundle), the PCA-axis packing flag > 5 A / > 35 deg (`diag_packing.py`), DSSP as
`tm_attribution.py`; reference the seed at t = 0 (6.2), also the implicit runs' first frame. Data:
* hybrid: the ff_3.0_gdepth chains (before 4.17's fix), rungs 0 (T 0.70) and 14 (T 0.8006), blocks 2-3 (chunks 20-42,
  ~90 tu each), plus the ff2.1_sc chain;
* implicit membrane: `implicit_79HIS_src{,2,3,4}` (stock Upside ff2.1 with `cb_`/`hb_membrane_potential`, WT, the old
  common `rama.dat`, 08-19), four independent 48-replica REMD runs, rungs 0 and 24 (T 0.799), 12,000 frames over
  600,000 tu, the first 500 dropped as `hdx_implicit.sbatch` does; O placed as `py/mdtraj_upside.py`, which reproduces
  the engine's O exactly. There is no ff3.0 or H145A implicit run.

Scripts, `contacts.json`, maps and figures: `scratchpad/tm5_tm6/` and midway2 `checks/tm5_tm6_contacts_20261013/`
(backbone npz there).
- **The start has a thin TM5-TM6 interface**: 5 CA pairs < 8 A in the seed (2 in the crystal, 1.27 A core CA-RMSD
  from the seed) and 20 < 10 A; TM6 packs mainly on TM4 (49) and TM3 (24).
- **TM5 stays helical.** Hybrid H145A any-helix 0.82 at both rungs (163-172 at 0.91-1.00, the ends 161-162 and
  173-175 fray; alpha 0.77 against 0.55, the T 0.80 loss at 167-173 going to other helix types). Implicit 0.99 and
  0.99 (alpha 0.99, 0.98), helical from 158 to 175. This matches the H145A HXMS TM5 peptide 162-169 (no uptake in 20
  d, 5.5).
- **TM5 leaves TM6 and packs on TM2, at both rungs and in both models.** CA pairs within 10 A of TM5, start / T 0.70
  / T 0.80: hybrid H145A TM6 20 / 3.6 / 5.9, TM4 15 / 4.5 / 6.5, TM2 0 / 17.0 / 16.8, L5 19 / 23.7 / 24.6; implicit
  TM6 20 / 6.3 / 5.4, TM4 15 / 18.1 / 21.5, TM2 0 / 28.8 / 28.9. TM5-TM6 native retention (10 A): hybrid H145A 0.117
  +- 0.015 and 0.188 +- 0.011; implicit 0.20 +- 0.02 and 0.16 +- 0.01 over the four runs. In the hybrid the loss
  happens within the first ~300 tu at both rungs. Other hybrid chains: WT ff3.0 0.17 / 0.13, H145A+S181T 0.21 / 0.23,
  WT+S181T 0.23 / 0.17, WT ff2.1_sc 0.34 / 0.26 (5 chunks). So the move is neither H145A's, nor ff3.0's, nor the
  dry-MARTINI coupling's.
- **In the implicit model TM5-TM6 is the weakest interface**: Q 0.20, the lowest of the nine TM pairs, against
  0.62-0.65 at T 0.70 for TM4-TM6, TM2-TM4, TM3-TM4 and TM2-TM3; in the hybrid every pair touching TM5 or TM6 is low
  (0.07-0.26) and only TM2-TM4, TM1-TM3, TM2-TM3 keep 0.41-0.68.
- **T 0.70 and T 0.80 maps differ no more than repeat runs at one temperature.** Hybrid: mean |dP| over present pairs
  0.15 between rungs against 0.18 between T 0.70's blocks 2 and 3; the largest differences are TM6's register on TM4,
  not TM5-TM6. Implicit, pooled: 0.05 between rungs against 0.07-0.09 between independent runs at T 0.70.
- **Every helix pair loses start packing in both models** (hybrid whole-protein Q 0.29 / 0.26 at 8 A, core CA-RMSD 1
  -> 4-5 A within ~300 tu; implicit 4.1). TM5's move onto TM2 comes with it, so neither model shows that H145A's TM5
  undocks from TM6 in the real protein.

### 4.17 Every hybrid seed places its side-chain beads at the wrong rotamer geometry (2026-10-17)

Found 10-17 by a read-only code survey, then measured. `load_sidechain_rotamer_payload`
(`py/martini_prepare_system_lib.py:3299`) wrote `placement_data` already gathered per row,
`rotamer_center_fixed[layer_index]` (one row per (residue, state), 747 for 79HIS), but stored `layer_index` as the
absolute library row (`lid = start + rel`, :3275; written :3717, :3727). The fixed placement node reads
`data[layer_idx]` (`src/placement.cpp:139-141`), so row r was placed at
`rotamer_center_fixed[layer_index[layer_index[r]]]`: a geometry taken from another row, in 740 of 747 rows from
another residue. upside_config writes the full 86-row library with absolute indices (`py/upside_config.py:1185`,
:1195, as master's), which is correct; the hybrid preparation since 7d332c5 (2026-04-03) was not.
`placement_fixed_scalar` (the 1-body placement energy) was indexed the same way.
- Measured on every seed the chains used (`seeds/glpG-RKRK-<V>.<ff>.up`, ff_2.1_sc, ff_3.0_gdepth, ff_3.1_depth_data,
  79HIS and 79ALA): bead offsets in the residue frame off by up to 5.59 A; in the engine (local live 79HIS seed) the
  placement misses the intended geometry by up to 7.7 A, 2.4 A on average.
- What it reaches: the side-chain bead positions read by the rotamer pair term and the side-chain/backbone H-bond
  coverage terms. The side-chain/lipid table (`martini_sc_table_1body`) evaluates at the separate CB node, which is
  per residue and correct. Protein-only configs (ConDiv training, the Peng benchmark) come from upside_config and are
  not affected.
- Every hybrid glpG result before 10-17 (TM1/TM4 helicity, glycine flips, HDX, under every force field) ran with it.
  Comparisons between force fields shared the defect, but it acts through each force field's own pair and coverage
  tables, so it does not cancel.
- **Fixed 10-17** (user: "fix the bug and rerun the affected simulations"; plan.md Phase 18). The preparation writes
  the whole library table (repo, and the hunk applied to the `$P` and `/beagle3` copies with
  `.bak_pre_placefix_20261017`). Seeds repaired by `checks/placefix_20261017/fix_placement.py`, which recomputes the
  library tables and refuses unless the seed's gathered rows equal them exactly: the four base seeds and the eight
  ff3.0/ff3.1 seeds; only the two `placement_data` datasets change, and `patch_glpg.py` on a repaired base equals the
  old seed repaired. Every repaired seed places its beads at the intended rotamer to 3.6e-6 A and its 1-body output
  equals the library table (the old: off by up to 22 E_up); frame-0 potential drops by 18-71 E_up.
- The old chains are kept as `<V>.<ff>.placedefect_20261017`, their HDX outputs (5.5) with them. Their replica files
  and HDX trajectory copies were deleted 10-18 (user); logs, STOP files and the HDX `results/`, `outputs/` and
  `inputs/` remain.

---

## 5. HDX: what the estimator measures and how to read it

### 5.1 The estimator is equilibrium, not kinetic

The `example/00.AnalysisScripts` uptake path does not read simulated elapsed time as exchange time. For each
amide donor, `get_protection_state.py` assigns a binary protected flag from backbone H-bond score, an Asp/Glu
side-chain-contact proxy and backbone/side-chain burial; `4.calc_D_uptake.py` MBAR-reweights those flags to
`p_protected`, computes sequence-, pD- and temperature-dependent intrinsic `k_chem`, then applies the EX2-like
`k_obs = k_chem * (1-p_protected)` and `D(t) = 1-exp(-k_obs*t)` in experimental seconds. So a wrong friction or time
mapping does not rescale the HDX time axis; it matters because it controls decorrelation and the ability to sample
opening/closing equilibria.

**The Rama map's undesigned transition region does not enter this estimator as a barrier (user question,
2026-10-05).** Under EX2, `PF = k_cl/k_op = 1/K_op`: both rates cross the same transition state, so its free energy
cancels and only the open-state population counts; REMD samples the Boltzmann distribution whatever the friction or
barrier heights, and the estimator pools frames without order (Persson & Halle, PNAS 2015, BPTI from O-state
populations; Peng 2022 the same estimator in Upside). The map's top (NDRD's untrained density tail, 1.19) can still
reach dG_HX in three ways:
1. As equilibrium weight, through open frames in cells with no library data. The top sits about 6-9 E_up (4-6
   kcal/mol) above each map's minimum, the range of the rare openings the TM4 comparison needs. Unproven; measurable
   by binning open frames by cell data count.
2. Through sampling: barrier heights set decorrelation time and kinetic traps. The check is the protection-state
   autocorrelation time, still unmeasured (5.3b).
3. Through the experiment's regime. Upside cannot supply physical k_op or k_cl, so EX2 has to be shown
   experimentally (Peng 2022: k_f against k_chem, EHEE_rd2_0005 1700 against 9-26 s^-1). No such check exists for
   GlpG, and the Sosnick-lab GlpG HDX-MS abstracts (Biophys J 2023, 2024; search snippets only) describe regions that
   "cooperatively unfold only at long time scales". If the TM4 peptides show EX1 (bimodal envelopes), their uptake
   reports k_op and the dG ~6.1-7.4 kcal/mol of 5.5 is not a free energy.

Order of checks: EX1/EX2 in the GlpG HXMS spectra before converting any TM4 peptide to dG; then the autocorrelation
time; then the open-frame Rama binning, which only matters once dG above ~5.5 kcal/mol is resolvable.
`4.calc_D_uptake.py:1046-1094` fits the target temperature to experimental uptake with `minimize_scalar`, which can
absorb thermodynamic error, so agreement after that fit is weaker evidence than agreement at a fixed temperature.
Full texts not yet read (requested from the user): Persson & Halle 2015 with SI, Peng 2022 SI (S11, S20-S22),
McAllister & Konermann 2015, Skinner 2012 (both), the GlpG Biophys J abstracts, Lin et al. JASMS 2025.

`dG` is `RT log(p/(1-p))`, not helix occupancy. At `T_up = 0.70`, 2, 3 and 5 kcal/mol mean roughly 98.4%, 99.79%
and 99.9965% protection. Exact `p = 1` from a finite trajectory is censored, not a measured 1000-kcal/mol value, so a
defensible plot must separate censored markers from finite dG points.

Two representation traps: donor IDs stored in `.resid` are zero-based while the VTF/PDB is one-based, and `T.npy`
is in Upside kT (values such as 0.85), *not* Kelvin, contrary to `example/00.AnalysisScripts/README.md`; feeding 303
instead of 0.8647 invalidates MBAR.

The hybrid feeds this machinery through a projection rather than a fork: the adapter builds a protein-only HDX-view
H5 whose `/input` comes from the ordinary `-HDX.up` config and whose positions are N/CA/C mapped through
`hybrid_bb_map/atom_indices[:, :3]`, with the full coupled potential, temperature and H-bond logs copied from the
hybrid group, so the stock tools see their native `3*n_res` contract while MBAR still uses the correct
protein+bilayer Hamiltonian. The stock protein-only `PS.npy` is kept beside any combined PS so a membrane correction
stays observable and reversible.

### 5.2 Resolution limits, clips and sentinels

Master keeps two conventions on purpose: the **T-slice** uses the clipped step-6 path
(`mean_pf >= 0.99999 -> sentinel`, a ceiling of `0.001987 * 298 * ln(99999)` = 6.82 kcal/mol at T = 0.85)
and the **full-temperature** profile uses the unclipped jscripts path. They are not interchangeable.

The clip coincides with the estimator's statistical limit. With a binary protection state the smallest resolvable
`(1-p_f)` is ~1/ESS, so the largest supportable dG is `0.001987 * T_scale * ln(ESS)`. On the CB-corrected ladder
(85 760 pooled frames, ESS 10 697 = 12.5% at T = 0.85) that is 5.5 kcal/mol against a clip at 6.8. Without the clip
dG reached 19.6, but the values above the limit are noise:

| dG band | n | median effective frames carrying (1-p_f) |
|---|---|---|
| < 2 | 111 | 1526 |
| 2-5 | 58 | 45 |
| 5-8 | 12 | 1.4 |
| 8-12 | 3 | **1.00** |
| > 12 | 5 | **1.01** |

Every residue above 8 kcal/mol had its value set by fewer than two effective frames, and the jackknife concealed it:
`jk_pf` is clipped to `1-1e-6`, so the reported error is exactly 0.0 above ~8 kcal/mol, maximum uncertainty
displayed as maximum confidence. Measured per-temperature limits on this ladder:

| T | resolution limit (kcal/mol) | residues at the bound |
|---|---|---|
| 0.70 | 4.33 | 84 of 203 |
| 0.75 | 4.81 | 69 |
| 0.80 | 5.13 | 58 |
| 0.85 | **5.49** | 27 |
| 0.90 | 5.57 | 28 |

An amide reaching the bound is **right-censored**: the data say "at least this protected". The bound is drawn as a
dotted line per temperature on the figure rather than by truncating the data. `plot_ref_style.py` renders the
unclipped profile; `calc_hdx_ht.py`'s `_DG_res_T_slice.png` keeps master's clipped convention, and the npz stores both
arrays so both renderings come from one MBAR solve. Smooth 10-20 kcal/mol bands are not reachable from a binary
protection state at any feasible sampling, since dG = 20 needs `(1-p_f) ~ 2e-15`.

Plotting trap: the out-of-range sentinels `+1000.0` and `-100.0` **are finite**, so `np.isfinite(...)` let every
unresolved residue into the connected `errorbar` series and the join painted full-height vertical lines. **A sentinel
encoded as a large finite number silently passes an `isfinite` filter**; here it had to be excluded in two places
because one array feeds both the line and the markers. The settled behaviour draws the line over all values (an
excursion off the top is how a non-exchanging amide is conventionally read) and error bars only where dG is
resolved; `plot_ref_style` reports "the reweighting resolved nothing at this temperature" when every residue is
off-scale.

### 5.3 The membrane term, and how its criterion was calibrated (findings 113)

The per-residue `acc` values in this section come from `martini_hdx_membrane_accessibility.py` before its off-by-one
fix (3.14 A61, 2026-10-20): glpG is numbered from 1, so each value belongs to the preceding residue. Re-measure them
before use, and recompute the HDX outputs built from those arrays (plan.md Known Errors). The criterion and its
calibration are unaffected.

`get_protection_state.py` scores an amide protected only if H-bonded or buried by *protein*. A TM amide facing lipid
is buried by neither, so it drops out of the protected state whenever its backbone H-bond flickers, though there is
no water there to exchange with. Master handles this with `--use-TM-region`, reading the `surface` node; the hybrid
HDX topology is a protein-only Upside config with no membrane potential, so that node does not exist. The
replacement is `combine_hdx_protection.py --water-accessibility` fed by `py/martini_hdx_membrane_accessibility.py`
(the driver once did not call it, which left `PS.npy` bit-identical to `PS_protein.npy`).

On one replica (T = 0.844, 1350 frames, unweighted, so nothing here is estimator behaviour) lipid shielding takes the
bilayer-embedded helical amides (lipid-shielded in >=90% of frames) at `+inf` from 35 to **79 of 92**, and the longest
contiguous `+inf` runs from 8, 6, 5, 4, 4, 3 to 17, 16, 15, 14, 9, 7: needles become blocks, the reference figure's
pattern. In helix 104-126 the four residues outside the bilayer keep their finite values while 108-126 go to `+inf`.

**Both numbers in the criterion are measured from the bilayer, not chosen:**
* **Radius = 7.00 A**, the flat first minimum of the intermolecular tail-tail g(r) (first peak 5.12 A). The minimum
  spans two 0.25 A bins whose ordering flips with noise, so the radius is the mean of the bins within 5% of the
  minimum, which gives 7.000 A on all 16 replicas where a bare argmin flips between 6.88 and 7.12.
* **Threshold = 1 contact**, because at that radius a phosphate bead has a median of 0 intermolecular tail neighbours
  (ester 2, first tail bead 6, terminal tail 11). The first tail contact is therefore the first step inside the head
  groups, i.e. the bilayer boundary itself.

`--cutoff` and `--min-contacts` are gone from `martini_hdx_membrane_accessibility.py`, so the criterion has no free
parameter. A bare slab is the wrong shape, rejected on measurement:

| criterion | helical amides at +inf / 148 | loop+turn at +inf / 55 | longest contiguous +inf runs |
|---|---|---|---|
| none (protein-only) | 41 | 0 | 8, 6, 5, 4, 4, 3 |
| slab between PO4 leaflet planes | 136 | 34 | 51, 46, 42, 22, 9 |
| slab between ester (GL1/GL2) planes | 117 | 22 | 42, 35, 22, 20, 19 |
| >=1 tail bead within 7.0 A | 115 | 10 | 22, 20, 18, 16, 13, 11 |
| >=3 tail beads within 7.0 A | 83 | 4 | 16, 15, 14, 12, 9, 7 |

A slab merges the helices into 40-50-residue blocks and erases the peak/valley structure, because glpG's interfacial
loops sit inside it too. Master does not use a bare slab either: `src/surface.cpp` selects TM residues by
`tm_min < z < tm_max` and ANDs that with a lipid-facing surface-exposure calculation, so a residue lining the
protein's own polar interior is not protected. The local tail-contact test is the coarse-grained analogue of that
conjunction.

With the membrane term in place the unclipped profile rises *continuously* into its excursions instead of jumping, so
the clip is not needed and only truncates. On the 16-replica ladder (177 k pooled frames, unclipped, T = 0.85): 85 of
203 off scale with 81 of them in helices, contiguous off-scale runs of 15, 15, 14, 12, 9, 8, a resolved range of
-0.53 to 21.3 kcal/mol, helix-interior / N-cap / loop medians of 4.44 / 4.26 / 1.63, and one helix-interior amide
below zero. What is left is the fold, and it is small: 13 of 92 deep helical amides stay finite at 2.3-4.2 kcal/mol.
The measurements of 4.2 stand but no longer explain the figure: with the lipid term in place, a flickering H-bond on
a bilayer-embedded amide is invisible to HDX, which is physically correct.

**The off-scale plateau is a lipid-burial map, not a protection map (2026-09-12).** Decomposing the conjunction per
frame over 171,382 samples (`protection_t = 1 - (1 - pp_t) * acc_t`: an amide exchanges only when protein protection
fails *and* it is water-accessible in the same frame):

| | `pp_fail` (H-bond/burial flicker) | `acc` | exchanged |
|---|---|---|---|
| TM1 30-48 | **0.0369** | **0.0021** | 3.13e-05 -> saturates |
| TM4 135-151 | **0.0322** | **0.4230** | 1.11e-02 -> stays finite |

The two helices flicker at the same rate (TM4 slightly less) and differ by ~200x in `acc` alone. Res 36 flickers
7.0%, the worst of any amide, yet `acc = 0.0000`, so it logs 0 exchange events and saturates; res 140 flickers 1.2%
yet `acc = 0.195`, so it logs 371 and resolves at ~4 kcal/mol. Because `acc_t = 0` makes protection identically 1
whatever the H-bond is doing, saturation is decided by lipid contact. **A `+inf` run is therefore not evidence of
secondary structure**, and "TM1 never exchanges" is wrong as a structural statement: its H-bonds break 3.7% of the
time and the lipid hides every break. This is the designed behaviour, but the protection map is only as good as the
7 A tail-contact test, which asks solely "is a lipid tail nearby" and leans on the protein-burial term to separate a
protein-interior amide from a solvent-exposed one.

Within the TM4 core the exposure is a helical face, not a depth gradient: `acc` for 140-147 runs 0.20, 0.39, 0.42,
0.016, 0.039, 0.654, 0.066, 0.001, so 141/145 (i, i+4) are both exposed and 143/147 both buried. One face is
lipid-packed, the other points into the protein interior and the catalytic cavity; 141, 142 and 145 are not
cavity-lining (5.3b).

**The same argument applies to any implicit-versus-hybrid comparison.** In the implicit model membrane burial is part
of the force field, so a burial-based protection state sees it for free; applying the protein-only criterion to both
strips the hybrid of protection the implicit model gets. (`--use-TM-region` cannot even the two up: neither the
legacy HDX topology nor the implicit configs carry a `surface` node, which `upside_config.py --surface` alone
creates.) Using each model's own membrane term inverts the result (Spearman is implicit against hybrid, 5.3b):

| T | n | Spearman | median implicit | median hybrid | hybrid never open |
|---|---|---|---|---|---|
| 263 K | 128 | 0.668 | 1.40 | 1.95 | 65 / 203 |
| 280 K | 130 | 0.691 | 1.28 | 1.91 | 61 / 203 |
| 298 K | 138 | 0.745 | 1.17 | 1.94 | 56 / 203 |

The hybrid is the **more** protective model and its helices saturate as they should; at 298 K the two models agree
better than in the protein-only version (Spearman 0.745 against 0.676). The hybrid is also the more censored (56 of
203 unresolved against 32, ceilings 6.2 against 8.0 kcal/mol), so a sampling contribution cannot be excluded. **"Hold
the analysis fixed" is not "hold the physics fixed":** two models can require different analysis in order to measure
the same quantity.

### 5.3b Implicit vs hybrid: what each one actually calls "protected" (verified 2026-09-12)

Established from the scripts and the saved arrays, not the records, after the `.md` notes twice misled on this point.

**The two protection rules have the same logical form.** Master's `get_protection_state.py`: `PS = HB1 + HB2 + BL`,
then `if use_TM_region: PS += Su`, then `PS[PS>1] = 1`, an OR over H-bond, side-chain H-bond, protein burial and
lipid-surface exposure. Ours, `combine_hdx_protection.py:31-32`:
`exchange = (1 - pp) * acc; protection = 1 - exchange`, i.e. `pp OR lipid-shielded`, with checks on finiteness, the
[0,1] range and shape equality. The combination rule is not the problem.

**The implicit run applied no membrane term.** `hdx_implicit.sbatch` calls `get_protection_state.py` bare; nothing in
the repo passes `--use-TM-region`; the saved implicit results (`/project/trsosnic/yinhan/implicit_79HIS_run3/results`,
48 replicas) contain only `_PS_protein.npy`, `_Hbond.npy`, `_Energy.npy`, `_T.npy`; and the implicit path runs its own
inline pymbar block (`_implicit_plain_pf.npy`), not `calc_hdx_ht.py` or `plot_ref_style.py`. So the lipid credit
belongs to the hybrid, the reverse of the natural guess, and a hybrid `+inf` run is a lipid-contact map (5.3).

**The TM4 signal is real and must not be suppressed.** Marking TM4 water-inaccessible because it is protein-buried was
rejected on measurement: the state-B frames (`pp=0 AND acc=1`) are temperature-activated (none across the 8 coldest
rungs, 71-89% in the hottest 7, res 140 climbing monotonically 0 -> 123), where a cutoff-margin artefact would appear
at all temperatures; `BL` already protects these amides in 98-99% of frames, so an override would change only the
1-2% that is the opening; and it would hard-code one protein's identity into shared `py/`. TM4 141/142/145 are **not**
cavity-lining: they carry 23-24 protein heavy atoms within 6 A, as many as the deeply buried 143/147 (22.6, 22.3), and
sit 6.87-7.34 A from the nearest lipid bead, a tail, against 4.90-5.16 A for 143/147, straddling the 7.00 A shell,
which is why `acc` lands near 0.5. The signal is local thermal opening at the edge of tail contact.
- Open: a cutoff sensitivity check at 7.5 and 8.0 A, reported as such and not as a recalibration (the 7.00 A is
  measured from this system's own g(r) and has no free parameter).
- Do not quote res 143: its 333 events are non-monotonic in T (76/93/83 at T=0.823-0.838, then 0-3 at the hottest
  rungs), which a genuine exchange signal should not be.

**The Spearman column of 5.3 is model against model, not accuracy.** `calc_hdx_ht.py:52-54` loads
`<pdb_dir>/<pdb_id>_{HXMS,NMR,NMR_MS}.csv` through `load_optional_numeric_csv`, and the `r_square` at line 575 is
computed only if one is present; on 09-12 `hdx/pdb/` and `hdx_postfix/pdb/` held only the `.pdb`, no such file existed
under `/project/trsosnic/yinhan` or `/beagle3/trsosnic/yinhan`, and no HDX log contained `r_square`. It was once cited
as evidence of accuracy, which was wrong. Dropping a `<V>_NMR.csv` into `pdb/` makes the r-squared and the scatter
appear with no code change. For that comparison the implicit arrays (2026-08-19, before the GLY, rigid-stage and
temperature fixes) must be re-run, and with 56-84 of 203 amides censored the regression can only use the resolved
subset, which is biased toward the least protected amides.

**Which model makes sense depends on the claim.** For lipid-dependent protection, cavity access, PE-vs-PG or the RKRK
variants, the hybrid is the only option: an implicit potential is a function of z and cannot distinguish two lipids
that differ only in headgroup, nor represent the lipid-facing / inward-facing asymmetry measured on TM4. For
intrinsic fold stability, `pp_fail` is the right observable and the implicit model is cleaner, cheaper, better
sampled and self-consistent, while the hybrid's lipid term obscures it. Against that, Upside's environment and burial
terms were trained against implicit solvent, so the hybrid is a chimera precisely in the terms that decide
protection. The hybrid saves both `PS_protein` and `PS` from one trajectory: state which one a figure quotes;
conflating them is what made TM4 look paradoxical.

**Cooperativity, measured directly.** One replica per model at matched temperature (implicit T=0.851 of 48 rungs,
hybrid T=0.853 of 28; same ladder span 0.700-0.900, mean 0.798), counting interior amides of a helix open in the same
frame:

| open amides in one helix, same frame | implicit | hybrid |
|---|---|---|
| 0 | 77.2% | 69.8% |
| 1 | 16.1% | 12.8% |
| 2 | 4.9% | 9.7% |
| 3+ | 1.8% | 7.2% |
| **mean given >=1 open** | **1.39** | **1.99** |

Both fall off monotonically with no second peak, so neither model opens a helix as a unit. Correlation
`P(i,j open)/(P(i)P(j))` by sequence separation 1..5: implicit `3.0, 2.7, 2.8, 1.6, 0.69`; hybrid
`4.3, 3.1, 2.6, 3.0, 3.0`. The hybrid is the more cooperative and longer-ranged; on the one-H-bond-at-a-time
criterion implicit is the closer match. Figure `fig_helix_opening.png` (09/14 deck), generator
`make_helix_opening_fig.py`, data `figs/coop_out.npz`. What differs across the 108 helix-interior amides is
freezing: implicit has 65% below one event per thousand frames and 35% at zero, against 24% and 9% for the
hybrid; but in 3 of 10 helices the hybrid is the more rigid and implicit fails outright where the hybrid does
not (res 126 0.577 vs 0.051), so implicit is close to all-or-nothing per amide while the hybrid breathes a few
percent throughout. "Implicit keeps all helical regions completely rigid" is false.

**Net position (09-12, pre-fix hybrid).** Implicit wins the textbook-mechanism comparisons (one-at-a-time,
short-ranged, frayed-terminus / rigid-core profile) and self-consistency. On TM4 the one GlpG measurement also favoured
implicit (5.5): Langosch's 3-4 (poly-Ala) to 5-6 (Leu) kcal/mol is measured in 80% TFE, with no lipid and no protein
packing, so it is not a TM4 reference; the lab's HXMS TM4 peptide 140-144 gives 6.1-7.4 kcal/mol, above the hybrid's
3.5-4.9 then and consistent with implicit's >5.4 (>7.8 if frames were independent). That bound depends on the
effective number of independent frames; **an autocorrelation-time estimate on the implicit protection state is still
the live quantitative discriminator.** Every error made in this comparison, the TFE reference included, ran toward
flattering the hybrid; remember it when reading any model-ranking claim in this file.

### 5.3c The protein carries no explicit charge in the hybrid, and what follows from it (2026-09-12)

Measured on both campaigns' configs, so systematic rather than one build's mistake. Of the protein beads in
`martini_potential/charges` only 10 are nonzero (2890 beads in the NP system and 1050 in glpG, identical in both), 5
at +1 and 5 at -1, summing to exactly 0: the two chain termini smeared across all five beads of their residues, which
cancels and is negligible.

**This is by design, not a bug.** `martini_sc_table_1body` is fully residue-resolved (`restype_order (18,)`,
`rotamer_full_energy_eup (18, 6, 38, 96, 13)` in the pre-audit table: 18 residue types x 6 rotamers x 38 environment
bead types). Protein-environment interaction is carried by tabulated per-residue-type fields, which is what the
"spline table only" rule requires, instead of explicit Coulomb. `charged_res` in `martini_prepare_system_lib.py:246`
is `{ASP:-1, GLU:-1, LYS:+1, ARG:+1}`; HIS is correctly absent. Albumin's sequence charge at pH 7 is -15 (LYS 58 + ARG
24 = +82, ASP 35 + GLU 62 = -97); the config does not represent residue charges explicitly.

**The consequence.** The SC-env tables are short-range radial fields, so albumin's 179 charged residues have no
long-range Coulomb term with the anionic MPA coating or with the ions, and changing the salt cannot create protein-NP
electrostatic steering: at 0.15 M the model's Debye length is 3.6 A and counterions-only would give 16.4 A, but
neither matters to a protein that has no charge to screen. Do not expect an ionic-strength change to fix the two NP
orientations that never bind (they start ~37 A off the surface); that is a model-design property plus starting
geometry. For the same reason a charge-driven footprinting prediction of the Carlson kind is not something this
hybrid can reproduce from first principles.

**One real defect, small:** the NP box carries net +15 e. The ion generator added 218 excess K+ on the assumption that
albumin is -15, while the charge array gives it 0. Either the protein should carry -15 explicitly or the ion count
should be 203.

### 5.3d The retired ff3.0's dG profile, and which ff2.1 profile is valid (2026-09-13)

The retired FF1-form ff3.0 (`<V>/hdx_postfix/results/`, 2026-09-12, symmetrised glycine maps; no longer on `/project`
or `/beagle3`, remote_jobs.md §1b) raised protection along the whole chain against the ff2.1 hybrid at matched rungs,
both after the rigid-stage fix: non-exchanging amides 53 -> 80 of 203 at T 0.75, 44 -> 73 at 0.80, 43 -> 69 at 0.85.
TM4 resolved finite in every run (censored 0/21 under ff2.1 and 1/21 under ff3.0 at T 0.85, resolved median 2.66 ->
3.03) and TM1 was the censored helix (4/22 -> 9/22); the four variants stayed indistinguishable over three frame
counts (62-68 of 203 off scale at T 0.85), so neither H79A nor S115T reshapes global protection. Quote one censoring
criterion: the `censored` flag and the sentinel count differ (66 against 105 at T 0.75).

**`glpG_POPEPOPG_dG_2026-08-27/` must not be used as the ff2.1 reference.** It predates the 2026-09-02 stage fix, so
its protein was frozen (3.4). The valid ff2.1 reference is the `<V>/hdx/results/` set (48 rungs, T 0.70-0.90). Its
79HIS files are dated 08-27 on `/project`, not 09-04 as first recorded; its protection arrays fluctuate, so it is not
the rigid-protein set of 3.4; which campaign it belongs to is not established.

### 5.4 Interpreting a per-residue profile

* **Do not use a bundled secondary-structure annotation as ground truth.** `hybrid_bb_map/bb_secondary_structure` is
  an idealised `CCCC1111HHHH...` pattern with helix segments of 41 and 50 residues, which is not glpG's topology;
  DSSP on the same structure gives ten helices of 5-23 residues. Re-scored against DSSP, 11 apparently broken helical
  donors resolve into 5 that are not in a helix at all, 8 within four residues of a helix N-terminus (an amide in the
  first four positions has no i-4 carbonyl to bond to, so fast exchange there is what experimental HDX measures), and
  only 3 genuine mid-helix breaks, the only ones of the eleven that were H-bonded in the prepared starting structure
  and lost it during the run. The bundled annotation turned 2% of donors into an apparent 8% failure and cost three
  days looking for a force-field cause.
* **To decide whether a feature in a reweighted observable is physics or estimator, recompute it with uniform
  weights.** Anything that survives is the trajectory's. Of the 23 up-excursions at T = 0.85, 14 amides never open
  in any of 86 333 raw frames and register zero transitions, so those are the trajectory's own statement; at the cold
  rungs most are estimator artefacts (78 at T = 0.70, of which 64 spurious), which is why the deliverable is plotted
  at the production temperature where reweighting is near-identity. Cold rungs earn their keep as ladder rungs for
  exchange, not as reporting temperatures.
* Helical donors read median dG 2.04 with 9% negative against loop donors' -0.05 with 51% negative, so the profile
  does track secondary structure.

### 5.5 TM-helix exchange in a bilayer: what the literature says (2026-10-10)

The user's question: do the TM spikes go to infinity, are the helices frozen, or do they open with water blocked.
Three literature agents; citations checked against abstract or full text unless marked.

**Experiment measures slow, finite exchange, and reports "none" only as the end of its window.** WALP16 in DMPC:
half-lives up to ~2x10^4 min, slowest in the bilayer core (Demmers et al. 2000 PNAS). M13 coat protein (SDS): 15-20
TM amides >=10^5 and V29-V33 ~10^6-fold slower than poly-DL-Ala (Henry & Sykes 1990 Biochemistry). APP C99 (LMPG, 45
C): TM half-lives hours to weeks, I45-I47 none in several weeks (Cao et al. 2017 JACS). bR helices: little exchange
after 48 h (Earnest et al. 1990 Biophys J); 71% very slow over 1.5-25 h (Downer et al. 1986 JBC). That is PF
~10^5-10^7, dG_HX ~7-9.5 kcal/mol at 298 K, a lower bound where nothing exchanged. GlpG: in SMALPs over 9 s-2.5 h,
TM2-6 protected and TM1 (94-108, E. coli numbering) relatively unprotected (Reading et al. 2017 Angew Chem); in
liposomes, many core amides still deuterated after 10 days in protonated buffer at pH 7.4, with TM5/TM6 exchanged
(Chevelkov et al. 2023 ACS Phys Chem Au). Hebling et al. 2010 Anal Chem is GGCX, not GlpG.

The lab's H145A HXMS (`~/OneDrive - The University of Chicago/glpG/Final H145A w reps ... all peps.csv`, pD 6, 10 s
to 20 d, 3 replicates; construct numbering = ours, verified: 79 A, 135 T, 188 H = H254):

| peptide | helix | %D 18 h | %D 29 h | %D 20 d |
|---|---|---|---|---|
| 140-144 ALMGY | TM4 | 26.5 | 41.6 | 97.6 |
| 109-113 / 110-118 | TM3 | 39.5 / 12.5 | 40.9 / 15.4 | 85.7 / 84.6 |
| 193-200 AVGLAMAF | TM6 | 52.8 | 62.8 | 71.8 |
| 35-64 (TM1 + loop) | TM1 | 13.2 | 14.3 | 19.9 |
| 162-169 GLIIFALI | TM5 | 1.6 | 1.3 | 1.7 (0.10 of 6 D, flat from 10 s) |

Not back-exchange corrected (the His-tag peptide 1-20 reads 13-20% at every time). TM4 140-144 half-exchanges at
~1.6 d (k_obs 4.8e-6 s^-1, steady from 18 h to 3 d); k_int for amides 142-144 from `4.calc_D_uptake.py`'s own tables
is 0.14-0.27 s^-1 at pD 6.0 and 20 C, up to 0.69-1.36 at pD 6.4 and 25 C, so PF ~3x10^4-3x10^5 and dG ~6.1-7.4
kcal/mol. The TM5 peptide is the only one with no measurable uptake; the overlapping 157-186 is too noisy to confirm
it (s.d. 11 D at 20 d) and it conflicts with the ssNMR's exchanged TM5, so it is unconfirmed.

**Exchange needs an opening and water/OH- at the amide, and the bilayer suppresses both.**
* Opening still matters inside the bilayer: Pro speeds exchange along the whole WALP TM segment (Demmers et al. 2001
  JBC); rates span three orders of magnitude along one APP TM helix (Cao 2017). Opening also costs more there (4.12).
* Access matters separately: protection follows helical periodicity at constant H-bond type (alamethicin in DOPC,
  Dempsey & Handcock 1996 Biophys J; melittin, Dempsey & Butler 1992 Biochemistry; M2, Tian et al. 2003 Protein Sci);
  a Gly packing defect raises local hydration (Hogel et al. 2018 Biochemistry). Langosch 2026 Chem Rev:
  membrane-embedded parts of a helix "are not accessible to the catalytic hydroxide ions", which is why that lab
  measures TM cores in 80% TFE (3-4 kcal/mol poly-Ala, 5-6 Leu), values with no lipid shielding.
* Core water is low, not zero: the Wiener, King & White 1991 Gaussian (centre 22.5 A, 1/e half-width 4.6 A)
  extrapolates to 7e-4 of peak at 10 A and 6e-7 at 5 A from the midplane (fit tail, not measured); Marrink &
  Berendsen 1994 put the water barrier at 26 kJ/mol (~1e-4 of bulk). No one has measured k_ch in the acyl core.

**What this means for our estimator.** `combine_hdx_protection.py` uses `exchange = (1-pp)*acc` with `acc` in {0,1},
so a lipid-shielded amide has k_ch = 0 exactly and its PF is infinite by construction. The logic has precedent (Peng
2022's burial term as the implicit-solvent stand-in for Halle's "no water"; master's `--use-TM-region`; the HDXer
tutorial, Lee et al. 2021 LiveCoMS, which allows lipid in N_C; Jia et al. 2023 JACS counting POPE atoms for XylE),
but none of these papers treats TM saturation, and the physics says k_ch in the core is reduced, not zero. So:
1. A TM `+inf` is a censored value, as experiment's "no exchange" is, and the measured TM-core range (6-9.5 kcal/mol)
   mostly lies above the estimator's ~5.5 kcal/mol limit (5.2), so most TM cores would be censored even with a
   nonzero `acc`.
2. A continuous `acc` needs a water profile for this bilayer. Dry MARTINI has no water, so the profile would come from
   outside this simulation, which 1.6 rules out unless decided deliberately. Not done.

**What the post-processing gets wrong in the light of this (recorded 2026-10-11; the user decides).** The protection
rule and the dG profile (`calc_hdx_ht.py`, `plot_ref_style.py`) stand.
* **Uptake cap.** `4.calc_D_uptake.py:499-500` set `k_obs = k_chem/1000` when `mean_pf == 1`, i.e. PF = 1000 (4.1
  kcal/mol at 298 K), so the most protected amides were predicted to exchange fastest among the protected ones (on
  glpG at T 0.85 ~85 of 203, nearly all TM, half-exchanged in minutes to hours against TM4's ~1.6 d). Upside 3 uses the
  per-temperature ESS bound instead, with censored amides carried as lower bounds (3.14 A64, 10-20); master keeps the
  cap. Even so, a censored amide's uptake is determined only to t ~ 0.1*ESS/k_int, about 10^3-10^4 s, so the multi-day
  HXMS points cannot test the TM cores at the present sampling.
* **Wording.** `py/martini_hdx_membrane_accessibility.py` ("the continuous +inf that a non-exchanging TM helix should
  show") and `helpers/plot_ref_style.py` ("a non-exchanging transmembrane helix") say non-exchanging where the
  literature says slower than the estimator resolves. `00.AnalysisScripts/README.md` (~113-117) rejects a
  phosphate-plane Boolean because "it saturates stable transmembrane amides", which is what the literature predicts,
  and does not mention the tail-contact criterion.
* **Interface boundary (a measurement to make, not a change).** The 7 A tail contact sits between the phosphate and
  ester beads, so an amide at the glycerol/carbonyl level counts as dry, where the Wiener-White profile still has a
  few percent of bulk water and non-H-bonded interfacial amides exchange near intrinsic rates (O'Neil & Sykes 1989;
  Huster 2003). It may over-protect helix ends. Check on existing trajectories: how many censored amides have an
  ester bead nearer than any deeper tail bead.

**TM4: protein-buried, finite in experiment, and more protected than the pre-fix hybrid said (user question,
2026-10-11).** The user's concern: implicit TM4 runs off-scale, hybrid TM4 does not, and TM4 is buried in the protein
with no bilayer contact. Numbering (10.38): residue numbers taken here from HDX outputs (dG, `pp_fail`, `acc`) are the
outputs' labels, 0-based Upside indices, one below the PDB numbering the HXMS peptides and the crystal use. Peptide
140-144 ALMGY is labels 139-143 and reports labels 141-143, so the pre-fix "142-144 at 4.1-4.4" below is PDB 143-145,
and every comparison of labels 142-144 with the peptide is one residue off.
* **TM4 is mostly protein-buried, but not sealed.** In the OPM-oriented E. coli structure in the lab folder
  (`glpg.pdb`, 91-271), relative SASA of construct 131-145 is <=0.05 (146-150: 0.09-0.56). Probe centres a 2.35 A
  tail-sized probe can reach in |z| < 13 A sit 4.1-6.3 A from the backbone N of 142, 143, 144, 146, 147 and 7.2-10.8 A
  from the rest of TM4, so the crystal predicts that ~5 TM4 amides register a tail contact under the 7 A rule, the
  ones the hybrid shields. 5.3's `acc` ~0 at 143, 144, 146, 147 came from the pre-A61 script, which stored PDB residue
  d's value in the column labelled d, so read as PDB numbers they name the right residues; on the placement-fixed
  ff3.0 chains at the 298 K rung (5.6), PDB 143, 144, 146 and 147 have a tail contact in 87-100% of frames, 142 in
  72-97%, and 145 and 148 in 28-77%. The hybrid's lipid term on TM4 follows the crystal; it is not lipid intrusion.
* **The off-scale implicit TM4 is not lipid.** The implicit analysis applied no membrane term (5.3b), so its TM4 is
  off-scale because H-bond or protein burial held in every frame: TM4 never opened, dG >5.4.
* **The finite pre-fix hybrid TM4 was protein protection failing, not lipid**: 5.3's 16-replica ladder (its `pp_fail` /
  `acc` table; no longer on `/project`) and the 48-rung ff2.1 set
  (`popepopg_REMD/glpG-RKRK-79HIS/hdx/results/`, 5.3d; at the 0.851 rung core 138-148 `pp_fail` 0.10-0.25, `acc`
  0.001-0.08, 142-144 at 4.1-4.4 kcal/mol). Its resolved 3.5-4.9 (2.7 median under ff2.1, 5.3d) means half-exchange
  in 8-170 min at 298 K and pD 6, against ~1.6 d measured, i.e. 1-4 kcal/mol too labile. Experiment (142-144 at
  ~6.1-7.4, above the ~5.5 limit) says a correct, well-sampled model shows TM4 off-scale or at the limit.
* **On the fixed table the lability does not reproduce (plan.md Phase 16, 2026-10-11; before 4.17's fix).** Local
  single-temperature runs, 79HIS, 12 seeds of 4000 tu, first 1000 tu dropped, ff2.1 HDX topology, the repo's HDX code
  (`scratchpad/tm4_hdx_local/`, readouts in `results/`). ff2.1 at T 0.85, 18,024 frames: 137-148 lose protein
  protection in at most 0.06% of frames, 135, 138 and 140-148 never exchange (dG > 5.8, counting every frame as
  independent); GLY136 resolves at 3.4; losses otherwise at 130-134 and 149-151 (1.0-2.9). ff3.0 (`d9_03_sc`, T 0.80)
  at most 0.012 in 135-147, its 148-150 fraying. So the pre-fix lability was the stale table and temperature, not the
  Upside core, as far as this sampling can tell. The matched ff2.1 REMD chain, `glpG-RKRK-79HIS.ff_2.1_sc`, was
  cancelled by the user (10-14) before its HDX; its files are since deleted (3.15).
* **REMD readouts on the misplaced-bead chains (4.17)**, block 2, 79ALA and 79HIS, 28 rungs T 0.70-0.90, the repo's HDX
  steps (`checks/hdx_ff30_20261014/`; outputs in
  `popepopg_REMD_mdw2/<V>.ff_3.0_gdepth.placedefect_20261017/hdx_sc600x12_ff_3.0_gdepth/` and
  `<V>.ff_3.1_depth_data.placedefect_20261017/hdx_sc600x12_ff_3.1_depth_data/`, `results/<V>_tm4_readout.txt`;
  figures on the Mac Studio in `~/Downloads/glpG_TM4_HDX_ff3.0_gdepth_block2/` and
  `~/Downloads/glpG_TM4_HDX_ff3.1_depth_data_block2/`). For ff_3.0_gdepth (10-14) 5.6 replaces it.
  - Mixing: KE/1.5kT 0.998-1.002 and swap acceptance 0.30-0.53, but no walker crosses the ladder (each spans 8-23 of 28
    slots): each rung is canonical, but its frames come from a few nearby walkers, so the ESS bounds count correlated
    frames and overstate what is resolved.
  - ff_3.0_gdepth, PS at T 0.85 (bound 5.85): the core censored (79HIS 137-148; 79ALA 135-147 except 137 and 139
    at the bound, 5.69 and 5.83, with 148 at 4.99); only the ends resolve. Unweighted, 137-148 lose protein
    protection in at most 0.37% of frames against the pre-fix set's 10-25%. TM4's dG finite, as experiment says, is
    beyond this sampling (kT ln ESS = 7.4 needs ~3e5 independent frames per rung). On upside2-md-87's criterion
    (137-147 censored above ~5.1 at 298 K, 148-151 no more frayed than locally) both variants pass, 79ALA with 150
    slightly more frayed.
  - ff_3.1_depth_data (10-17, the old ff3.1): the core 140-147 censored as under ff3.0, but the N-terminal turn
    resolves (79ALA 136 4.36, 138 4.47, 139 4.89; 79HIS 135 4.96, 136 4.03, 139 5.58) and 79HIS's 148 drops to 2.72;
    79ALA's protein-protection losses reach 4.6% and are not stationary over block 2. It frays TM4 more at both ends,
    by HDX as by DSSP (1.35); 79ALA fails upside2-md-87's criterion at 138-139.
* **A residual cannot come from the core's trained parameters alone:** the implicit run uses the same core and keeps
  TM4 closed. What differs is how each model embeds the core: the implicit membrane's backbone H-bond cost at the
  midplane, where TM4a sits (|z| 1.47 A, 3.10c; the term raised hybrid TM4 0.71 -> 0.82 on the old table, 4.8, and
  0.924 -> 0.959 n.s. on the converged one, 4.14); helix-helix packing balanced only against soluble proteins while
  lipid attraction is full dry-MARTINI (4.4's open hypothesis); and the untested rotamer 1-body node (4.1). The
  missing environment terms are ruled out for TM4 (4.4). The one core parameter that acts on TM4 specifically is
  glycine's map (3 Gly; 3.10c), the same in both models.

### 5.6 One REMD block does not converge glpG's dG: block 2 against block 3 on the placement-fixed ff3.0 chains (2026-10-22)

User 10-22: the ff3.0 glpG dG plot, "and I will see if we need longer run". The four ff_3.0_gdepth chains on the
repaired seeds (4.17), block 2 and block 3 each through the 10-14 pipeline on its own (5,370-5,777 frames per rung;
remote_jobs.md §1 for jobs, groups and scripts). The figures are blocks 2-3 pooled through the same pipeline, one per
variant in the project's format (`plot_ref_style.py`, T 0.75 / 0.80 / 0.85, PS), with `summary.txt` (the block numbers
below) in `~/Downloads/glpG_HDX_ff3.0_blocks_20261022/` on the MacBook Pro.
- Mixing (`checks/hdx_ff30_placefix_20261022/`): KE/1.5kT per slot 0.997-1.004, swap acceptance 0.33-0.52, no 0 -> 27
  -> 0 round trip in any block. Over blocks 2-3 together 1-2 of 28 walkers per variant cross the ladder once; the
  median walker spans 15 of 27 slots.
- PS at 298 K (bound 5.81-5.88): 63-74 amides, the TM cores, are censored in both blocks and 14-19 in one only. Of the
  113-124 resolved in both, block 3 - block 2 has RMS 0.68-0.91 kcal/mol (max 2.3-3.8, mean -0.07 to +0.15), while the
  combined jackknife sigma has median 0.11-0.14, so 57-76% of them differ by more than 2 sigma. PS_protein: 191-196
  resolved in both, RMS 0.65-0.88, 59-74% beyond 2 sigma.
- Outside the TM helices the difference averages to zero (-0.12 to +0.18; RMS 0.44-0.66). Under PS_protein whole
  helices move together: 79ALA TM1, TM2, TM3, TM6 -0.43 to -0.65, 79ALA_S115T TM3 -0.60, 79HIS TM2 -0.53.
- Read: the jackknife over the 28 slots understates the per-block uncertainty about 6 times; one block's per-amide
  uncertainty is about RMS / sqrt(2), 0.46-0.64 kcal/mol, and the helix-wide shifts with no ladder round trip say the
  blocks are not independent samples of one equilibrium. One block cannot resolve a variant difference under ~1
  kcal/mol at an amide. The censoring bound rises only as RT ln(ESS), 0.41 kcal/mol per doubling: if ESS grows with
  frames, two blocks reach ~6.3 (the pooled run: 6.25-6.29) and four ~6.7, while TM4 140-144's HXMS 6.1-7.4 (5.5) needs
  ~14 blocks' worth for its upper end.
- **TM4** (user 10-23: "it looks like the TM4 issue is resolved according to plots glpG_HDX_ff3.0_blocks_20261022 ...
  Double check"). Residues here are PDB numbers, the HDX outputs' label + 1 (10.38). The figures draw `dG_unclipped`:
  their 11-19 kcal/mol TM4 points lie above RT ln(ESS) and are lower bounds (`plot_ref_style.py`), so they read as the
  npz's censored values do. At 298 K, pooled, PS:
  - 141-144 are censored above 6.25-6.29 in all four variants, as at T 0.75 and 0.80 and in each block alone, and 145,
    147 and 148 too, except 79ALA's 145 (5.68) and 79ALA_S115T's 148 (5.62). The HXMS peptide 140-144 reports 142-144
    at 6.1-7.4 (5.5), so the model agrees with it.
  - The N-terminal turn resolves: 138 at 3.71 / 5.09 / 4.43 / 4.81 and 140 at 4.80 / censored / 3.82 / 5.13 (79HIS,
    79HIS_S115T, 79ALA, 79ALA_S115T), and 79ALA_S115T's 146 at 4.76. On the misplaced-bead chains (5.5) they were
    censored or at the bound, so plan.md Phase 16's rule (labels 137-147, PDB 138-148, censored above ~5.1) now fails
    at 138 and 140 in 79HIS and 79ALA. 149-152 read 1.3-4.6 where resolved.
  - The core is off-scale only through the membrane term. Under PS_protein 79ALA's 141-145 resolve at 1.77-4.56 (the
    peptide's 142-144 at 2.65, 1.77, 4.56); at the 0.853 rung their protein protection breaks in 0.05-3.7% of frames
    (79ALA_S115T's 141-142: 8-11%), and in 141-144 every such frame is lipid-shielded.
  - That term's accessibility is the pre-A61 one (the shared deployment's script, md5 `ca115325`): each amide carried
    its N-terminal neighbour's value. Undone from the stored arrays at the 0.853 rung, 141-144 still exchange in 0-4 of
    ~11,000 frames, but Val145 goes from 0-2 to 0 / 3 / 18 / 11 (79HIS, 79HIS_S115T, 79ALA, 79ALA_S115T), and
    79ALA_S115T's 139 from 0 to 44 and its 148 from 1 to 9. The undoing is exact wherever the next residue is a donor,
    all of TM4 135-151: the run is NVT, so the box is the input one, and the radius stays each replica's own (7.00 A;
    6.88 or 7.12 on 1-4 of 28).
  - Reweighted on those arrays (`hdx_accfix.sbatch`; `$C2/tm4_accfix_blocks2-3.txt`, in labels), at 298 K: the
    peptide's 142-144 stay censored in 79ALA, the experiment's variant (E. coli H145A), and in 79ALA_S115T; in the
    79HIS variants 143-144 do and 142 resolves (5.43, 4.00). 146-148 stay censored but for 79ALA_S115T's 148 (4.49).
    **Val145 comes down**: 4.02 in 79ALA, 4.14 in 79ALA_S115T, 5.00 in 79HIS_S115T, censored in 79HIS; its protein
    protection breaks as before (PS_protein 3.1-3.3), and with the shift undone it has no tail contact in 23-72% of
    frames. No HXMS peptide isolates 145: 143-154 reports 145-152 and 154 at once (10.2% D at 10 min, 15.4% at 1 h,
    19.7% at 6 h, 31% at 48 h, 46% at 20 d, not back-exchange corrected), and the model's 150-152 at 1.5-2.9 would
    exchange in minutes. Phase 16's rule fails in every variant: 79HIS at 138 (4.13), 79HIS_S115T at 142 and 145, 79ALA
    at 140 (4.06) and 145, 79ALA_S115T at 138, 139, 145 and 148 (4.1-4.7). The dG code is not touched by A62 (its
    residue kT is 0.001987 x 298 K at T 0.85; the 0.593 sits in the free-energy surface plot only).
  - Read: in the experiment's variant the hybrid under ff3.0 agrees with HXMS on TM4's one isolated peptide, with or
    without the accessibility fix. What it frays (the 138-140 turn, Val145, the C-terminal end) no peptide resolves. The
    core holds by lipid shielding while its protein protection breaks far more than in the implicit model, on an engine
    without A83's ALA and GLY side chains (plan.md Known Errors; TM4 135-151 has four). The user closed the TM4 HDX
    question on this reading 10-23; the protein protection is read next on Upside 3's glpG run.
- `extract_block.py`'s `print` wrote its line and newline separately under `PYTHONUNBUFFERED=1`, so two of 16 parallel
  replicas shared one line of `extract.log` and the `verified$` count failed a correct copy. It now writes each line in
  one call; the check is unchanged.

### 5.7 What Upside REMD on ubiquitin gives an HDX training term (M1, 2026-10-28)

User 10-28: HDX first, then large proteins, then IDRs, as terms in the training (`upside3.md` 4.1); "check if you can do
analysis with" the midway2 results. Data: the Peng-benchmark arms `ubiquitin_native_ff_2.1` and
`ubiquitin_native_ff_3.0_gdepth` (14 rungs 0.800-0.940, 2,429,200 tu, frames every 100 tu, 24,287 per rung; WT 1UBQ,
where experiment is Ub psWT F45W/H68N) and `ubiquitin_denovo_ff_3.0_gdepth` (8,028,200 tu), run by the Upside 2 engine;
the stock estimator (`get_protection_state.py`, topology from `1.config.py` with each arm's own force-field files,
md5-equal to the shared deployment's) on Upside 3's writer and engine; MBAR on the Mac. Rules fixed before the data
(`analyze.py` docstring): first 10% dropped; targets 273 K (experiment) and 300 K (Peng's FF2 Table S10) at
T_up = 0.85 T / 298 K; censored at RT ln(ESS), Kish ESS; folded below 5 A CA-RMSD over 2-73. Everything in
`scratchpad/hdx_term_20261028/` (`fetch.sh`, then `run.sh`; `results/results.txt`, `results/ubq_dG_measured.png`).
* **The pipeline reproduces Peng's FF2.** ff2.1 at 300 K against Table S10: Spearman 0.90 over the 24 measured amides;
  against experiment, censored Spearman 0.515 and Pearson 0.460 over 22 resolved, where the published values give
  0.492 and 0.455 (Table S14: 0.455).
* **Resolution.** 3.06e7 replica-tu give a Kish ESS of 3.1-3.3e4 at 273 K, bound 5.6 kcal/mol, with 14 (ff2.1) and 15
  (ff3.0) of the 24 measured amides resolved; experiment spans 4.9-8.5. At 300 K the bound is 6.9-7.3 and 22-24
  resolve. The rung-0 open count decorrelates in 2.3 (ff2.1) and 13.7 (ff3.0) frames, the rung-0 RMSD in 1,057 and 639,
  so the Kish ESS is an upper bound for anything slower than a local opening. Resolving 8.5 kcal/mol at 273 K needs ESS
  6.4e6, about 200 times this run; one ConDiv step (12 x 8000 tu) is 1/319 of it.
* **ff2.1 opens the globally exchanging amides locally.** Of every amide's open weight at 300 K, 0.81-1.00 comes from
  frames under 5 A, and 0.96-1.00 for the seven experiment assigns to global unfolding by their urea slope (I3, F4, V5,
  L15, D21, K27, I30 at 0.70-0.92 of m_global; `scratchpad/boltz_hdx/flags/dg/urea_classes.log`). The model opens them
  by fluctuations of the folded state at 3.5-7.0 kcal/mol (K27 censored at 7.3) where experiment needs global
  unfolding at 7.7-8.5, the failure Boltz-2 showed on the same protein (`diffusion_dynamics.md` Phase 1 urea item). The global unfolding itself is censored at 273 K (P(RMSD >= 5 A)
  2.0e-6, dG >= 5.64; experiment's Table S4 dG_global 8.556).
* **ff3.0's ubiquitin melts about 0.06 T_up lower than ff2.1's.** Folded fraction per rung: ff2.1 1.00 to 0.899, 0.99
  at 0.910, then 0.91, 0.34, 0.08; ff3.0 0.90 at 0.857 and 0.23 at 0.867, 0.00 from 0.910. A native start biases toward folded, so
  ff3.0's transition lies at or below about 0.865. The bottom-rung TM-score the benchmark reports (0.697 against
  ff2.1_ctrl38's 0.742, `checks/ff21_released_20261021/peng_trained.txt`) does not show it. ff3.0's ladder makes 3 round
  trips against ff2.1's 806, and its dG halves drift by 0.52 kcal/mol (ff2.1 0.00), so its per-amide values are not
  converged; at 300 K all 24 read 1.0-2.1 kcal/mol, opened by global unfolding (folded share at most 0.1 for 18 of 24).
  The cause is not measured; findings 1.37's folding loss under retraining is the candidate.
* **Search fails at 76 residues.** ff3.0's de novo arm makes 4,102 round trips in 7.2M tu and never comes below 4.4 A,
  while its native start refolds below 3 A eight times after passing 8 A.
* **The estimator's definition moves with the force field, by little here.** It reads the trained burial weights
  (`sigmoid_coupling_environment.weights`), which differ by up to 0.20 between ff2.1 and ff3.0; on the same ff3.0
  frames 0.39% of flags change and dG moves 0.12 kcal/mol RMS at 273 K (0.03 at 300 K).
* **The comparison temperature carries the stability error.** Peng reads FF2 at 300 K against an experiment at 273 K,
  the temperature chosen where the model's global stability matches (SI Table S4 and p. 15). At 273 K ff2.1's censored
  Spearman is 0.394, with 10 amides censored.

---

## 6. Reusable diagnostic procedures

### 6.1 Glycine Ramachandran maps: mirror convention and checking traps

"Symmetrise every glycine map" was the retired ff3.0 (GLY_sym.md §4), and its code is gone from `upside_config.py`;
ff3.0's library makes only `GLY|GLY` mirror-symmetric, checked by `build_gly_library.py`. What survives:
* **The mirror is `(phi,psi) -> (-phi,-psi)` with a roll.** On the 72-point grid starting at -180, index i maps to
  (-i) % 72: `np.roll(np.roll(m[..., ::-1, ::-1], 1, -2), 1, -1)`, or `m[np.ix_(mir, mir)]` with
  `mir = (-np.arange(n)) % n`. A plain `m[::-1, ::-1]` maps i to n-1-i: one bin off, it creates ~0.6 E_up of
  asymmetry in alpha_R's favour and reads 0.0000 under its own metric. The same expression serves library `dimer_pot`
  and `.up` `rama_pot`. Check against a chiral control (ALA keeps ~11 E_up of asymmetry) and against dG(aR->aL).
* **Grid indices use `int()` (floor)**, as `inject_backbone_nodes` does, not `round` (25 against 24 for phi -57).
* **The hybrid `input/pos` backbone is stride 4** (N=4i, CA=4i+1, C=4i+2, proxy=4i+3), pure protein stride 3; read
  N/CA/C from `hybrid_bb_map/atom_indices`.
* **Silent failures of per-residue map edits**: 1-indexed residue numbers used as 0-indexed h5 lookups; a phi-range
  filter or an alpha_R > alpha_L guard (every boundary misses some glycine); the wrong mirror above.
* **The library's NaN is one whole neighbour column, `CPR`** (1/22), never read because a cis-proline neighbour maps
  onto `PRO`; a naive `abs(a-b) > tol` diff reports no difference across it (NaN compares False).

### 6.2 Verifying TM helix health in a glpG trajectory

Procedure and scripts: `remote_jobs.md` §5 (`glpg_tm_windows.py`, `glpg_ff30_health.py`). TM1 30-48 and TM4 134-151,
helix = phi in [-130, -20] and psi in [-90, 15], pass > 0.8 at T 0.70 with no downward trend across blocks, reading
only rotated `output_previous_*` groups and skipping `output_previous_0`. Traps (memory `glpg-vtf-reading-traps`):
* the old handbook `dihedral()` helper and `check_seeds_current.py` return -phi_std;
* GLY49 and GLY133 are helix caps at positive phi in the crystal seed, so there is no glycine-phi criterion;
* a protein "leaving the bilayer" in a VTF is the N-terminal tail; a VTF lists N, CA, C, O per residue from atom 0,
  and phi of residue i uses C(i-1), N(i), CA(i), C(i);
* a REMD replica's `/input/pos` is not t = 0: `run_remd.py` rewrites it at every chunk restart, so contacts read from it
  differ by rung (TM2-TM4 21 against 31) and drift. Take the reference from the chain's seed
  `popepopg_REMD_mdw2/seeds/glpG-RKRK-<V>.<ff>.up` (= the simulation PDB, 11).

### 6.3 Detecting a frozen or rigid protein

Kinetic temperature, finite coordinates and retained secondary structure can all pass while the protein does not
move. Measure motion: internal CA RMSD between frames whole chunks apart (0.000-0.001 A is a rigid body, 3.4; live
2.46 +- 0.563 A); SDs of projected observables (Rg 17.60 +- 0.000 and H-bond count 194.6 +- 0.01 frozen, +- 0.228 and
+- 12.07 live); drift-removed lipid COM motion (0.0081 A/frame frozen, 3.4-3.7 A net RMS working); and that
`/input/stage_parameters.current_stage` reads exactly `production`.

### 6.4 Health gates: count broken peptide bonds, do not use a worst-bond threshold

On a forced NP tear (dt = 0.01) the count (>= 5 bonds > 2 A) and the worst bond (> 3.5 A) both first fired at t = 168.0;
a potential test (non-finite or > 1e6) and |coord| > 1e4 A never fired.
* **The worst-bond threshold is redundant and fragile**: healthy max 2.659 A against torn 3.93 A; at 2.5 A it
  false-fired on a healthy chunk and cost a 6 h glpG block.
* **An energy test cannot guard NP** (431 broken bonds at a finite +3e5), while glpG goes fully NaN within one 46-step
  interval, so `isfinite` suffices there. Different jobs need different detectors.
* **The count is not a tuned knob**: healthy frames have 0-2 stretched bonds (all six NP systems and healthy glpG), a
  torn one 279-431, and any cut between 3 and 200 behaves identically.
* No invented magnitudes (`CN_MAX`, `POT_MAX`, `COORD_MAX` deleted). NP fires on non-finite positions OR >= 5 stretched
  bonds; glpG on a non-finite potential anywhere in the chunk OR >= 5 stretched bonds in the final frame.

### 6.5 Detecting a destroyed run when `isfinite` passes

* **`isfinite` is not a health check**: a finite NP file was destroyed (Rg 26.5 -> 76.7 A, max peptide C-N 56 A,
  potential -12 000 -> +1e5), and glpG environment coordinates at a failed frame reach +-4.65e12 A.
* **The sign of the total potential is the cheap physical test**: bilayer plus protein sits near -2.2e4 E_up, so a
  positive total means core overlap (`martini_remd_concat.py` and the reseed script filter on it). One capture
  oscillated -2.17e4 -> +5.4e3 -> recovered -> +1.98e6 -> NaN over ~96 frames, so a spot check can miss it. The C-N
  count passes on such frames (the damage is in the lipids): run both tests.
* **A collapsing adaptive chunk size is a free tell** (394.9 -> 76.0 time units per chunk). Gate a self-resubmit on a
  physical health check so a dead run stops.
* Output groups in a restarted `.up` ascend: `output_previous_0` is the oldest, `output` the newest.
* **Post-NaN diagnostic lines are not evidence**: the pair diagnostic then prints `max_force 0` with indices `-1`
  (`NaN > max` is false) and a repeated `min_dist 0.0010`. Only the pre-NaN escalation is a measurement.
* Term decomposition beats global observables for a local failure (3.2).

### 6.6 Identifying which force field a running replica carries

`run_remd.py` copies a replica from the seed only if the file does not exist, so a chain resumed after a force-field
install keeps the tables it was built with. Test by a least-squares scale match on the rotamer pair table, which
`upside_config` rescales into Upside units:

```python
with h5py.File(f'{P}/{tag}/sidechain.h5') as f: ref = f['pair_interaction'][...]
with h5py.File(up) as f:
    a = f['input/potential/rotamer/pair_interaction/interaction_param'][...]
c = (ref.ravel() @ a.ravel()) / (ref.ravel() @ ref.ravel())   # best scale
```

The right force field gives `c = 1.000000` with a residual at float32 rounding (1.6e-6), the wrong one `c = 1.107` and a
residual of 21. **A verbatim hash sweep of the parameter files gives the wrong answer** (it reported ff_2.1 for an
ff_3.0 seed): the files that match are shared (`martini.h5` comes from `ff_2.1` by design, `hbond.h5` is identical),
while the two that differ, `sidechain.h5` and `environment.h5`, are transformed before they enter the `.up`. Compare the
transformed tables, not the files.

### 6.7 Glycine AWH: convergence lessons (2026-09-18/19)

The AWH measurement of glycine's handedness (GLY_sym.md §5; final numbers in 9r and 9s) went through wrong answers
(12a). The methods that caught them:
* **Check the estimator's dynamic range first.** `awh1-dimN-diffusion = 5e-5` (AWH's friction metric implies ~0.77
  rad^2/ps) and a small `awh1-error-init` gave a PMF spanning ~2 kJ/mol at 25 ns against the 20-40 of a glycine
  surface, and a flat surface is mirror-symmetric ("handedness zero"). `diffusion = 0.5`, `error-init = 30` gave 50-66.
* **An achiral control passing is necessary, not sufficient**, because its errors cancel by symmetry. Only independent
  replicas at matched sampling bound the error.
* **Replica agreement at matched time is not convergence**: replicas drift together (agreement at -0.31 at 36-45 ns,
  plateau only from 62 ns at -0.269). Run a single-replica series past where it stops moving, with a plateau longer than
  the feature, and check the control's own series before calling its offset systematic.
* **The blank is the convergence criterion.** Its residual was uncorrelated between replicas (r = +0.157) and decayed as
  1/sqrt(t) (0.233 at 10 ns to 0.032 at 100) while the signal converged (0.080 to 0.071): sampling noise, so blank
  subtraction only adds noise. It puts ~20% uncertainty on the basin dG.
* **Per-neighbour structure is noise; only the average reproduces** (replica orderings Spearman rho +0.048, p 0.91).
  When the claim is about every pair, plot every pair.
* **Quote wander over the longest window** ("LA is settled" at a 0.028 range was 0.106 over 18 ns).
* **A derived quantity's definition lives in one place**: two scripts' different basin boxes differed by 0.106, which
  the achiral control could not flag.
* **On a periodic grid the boundary column breaks mirror symmetry** (counting phi > 0 cells gave P = 0.481-0.486 for
  achiral controls). Split the -180 column between the two halves.
* `gmx awh` on midway2 silently writes no `fe_t*.xvg` without `module load gcc/10.1.0`; after an extension, read the
  last part file (`ls -v awh.part*.edr | tail -1`).

## 7. System preparation

### 7.1 Environment morphology is derived from topology, not chosen

DDM is retired as an environment of this model (2026-09-09) and glpG runs in a POPE/POPG bilayer only; this section
keeps the rule that stops a single-tail detergent from being built as a slab. The prep counts acyl chains as connected
components of the apolar (`C1`-`C5`) bond subgraph in the lipid ITP: one tail means micelle and is refused
(`require_bilayer_lipids`), two or more mean bilayer. Deciding from topology rather than a name table or a CLI flag
makes the unphysical combination unreachable. `parse_lipid_from_itp` already returns **0-based** bead indices;
subtracting 1 again split DDM's single tail in two and reported a bilayer.

**A DDM lamellar slab cannot solvate GlpG** (findings 76, cited as "Update 76" in `example/16.MARTINI/readme.md`): a
CHARMM-GUI DDM slab has a 12.7-13.9 A tail core against a 28-30 A TM belt, and TM4, the most buried helix, alone failed
(helicity 0.90 -> 0.17, internal RMSD 4.45 A), every unwound residue at or beyond the core's edge. It is not tunable:
lamellar thickness is `2*V_tail/APL`, and with `V(C12) ~ 324 A^3` a 28 A core needs APL ~ 23 A^2, which a maltose head
(>= 40 A^2) cannot reach. Experiment agrees (HXMS TM4 peptide 140-144 half-exchanges at ~1.6 d, 5.5).
`--membrane-thickness-angstrom` does not set slab geometry; it is read only for ion counting.

Building a micelle: seed from the shell volume, not convex-hull support points (32 molecules against 186); fill
innermost outward; take the packing distance from the force field's `2^(1/6) sigma_max` = 5.276 A, not a CHARMM-GUI
step5 template (pre-minimization, bead pairs 0.24 A apart).

**A packed-state thickness span cannot be gated on**: the 5-95 percentile tail-bead z span against OPM's hydrophobic
thickness fails DOPC (20.4 A), POPE/POPG (20.6 A) and DDM (11.0 A) alike against a 22.9 A limit, because CHARMM-GUI
templates are laterally compressed. The shipped check is `assert_environment_solvation` on equilibrated coordinates per
belt site (current form: 8b). **When a build-time gate and the thing it guards are measured differently, the gate is
worthless**: measure both on the same state or demote it to a report, and check a new geometric criterion against the
paths it must not break.

### 7.2 A periodic tile's box is the tile

The bilayer prep read the template's CRYST1 only for tiling and sized the box from the **lipid coordinate extent**,
squared to the longer edge: an 84.41 x 73.10 A tile became 84.26 x 84.26, a 15% stretch along y that opens a vacuum
stripe and moves the area per lipid from 61.70 to 71.0 A^2. Fixed by using CRYST1 as the lateral box; a template
without CRYST1 is a hard error. An earlier "box inflated to 83.3 A for a 75 A tile" was the same bug, patched by
wrapping molecules into the tile. **When a fix works by making a wrong input look right, the cause resurfaces on the
first input that breaks the coincidence. Fix the derivation, not the input.**

### 7.3 Lipid bead models must match the ITP before anything else

POPE/POPG had never run through this pipeline, because the bead models disagree:

| source | POPE beads | tails |
|---|---|---|
| Robertson `last_frame.pdb` | 12: NH3 PO4 GL1 GL2 C1A **D2A** C3A C4A / C1B C2B C3B C4B | 4 + 4, unsaturation in tail A |
| CHARMM-GUI Martini Maker | 12: identical to the above | 4 + 4 |
| our `dry_martini_v2.1_lipids.itp` | 13: NH3 PO4 GL1 GL2 C1A C2A C3A C4A / C1B C2B **D3B** C4B C5B | 4 + 5, unsaturation in tail B |

DOPC matches exactly, which is why DOPC and DDM worked. The ITP is the authority and coordinates are generated for it
(adding the 12-bead model would run a lipid the force field was not parameterised for; substituting DOPC loses the
composition). **A lipid being present in the ITP does not mean the pipeline can build it**; one `diff` of the two
bead-name lists would have found this in the first minute.
* **Rigid-rod conformers cannot start at the target area** (at APL 69.4 A^2 the lipid axes are 8.33 A apart while a
  rigid conformer is ~6 A wide). Start deliberately loose and condense with a pure-bilayer run under the xy barostat,
  then use the equilibrated tile as the template (pure bilayer only; with a protein the box is pinned by its footprint).
* **Compute the comparison metric on the reference data before tuning your own builder.**
* A deposited composition may not be the nominal one: 1814 POPE : 959 POPG is 1.892:1, not 2:1.

### 7.4 Ions and box size (NP-1AO6)

Final: **neutralizing counterions only, no bulk salt** (user decision 2026-08-04): 218 K+, 0 Cl-, cancelling MPA (-203)
+ protein (-15). Salt molarity was never the defect; box volume was, since `box_len` came from the protein's reach about
the gold COM, doubling the lever arm of a protein adsorbed to one side (boxes 232-284 A for a complex 122-148 A across).
Now complex-centred at a fixed 200 A, and `build_system` asserts exact neutrality and zero salt pairs from the placed ion
counts and rejects a box too small for the complex. It survived rebuilds because nothing asserted the built
composition. 200 A is itself too small for states the model visits (Rg 152 and 206 A), and enlarging a near-vacuum box
costs almost nothing.

### 7.5 Dry Martini forms the SAXS DDM micelle; around a frozen Martini glpG its maltose heads coat the belt (2026-10-26)

plan.md Phase 31 step 2: GROMACS 2026.4 on the Mac with the potential Upside reproduces (8: reaction field eps_r 15,
potential-shifted LJ to 1.2 nm, `sd` at 303.15 K), cgmartini's DDM topology, NVT, 200 ns each (`scratchpad/ddm_micelle/`).
* **DDM alone is right**: 140 DDM from random placement form a micelle of 124 (plus one of 16) within 20 ns, steady to
  200 ns; tail-core semi-axes 14.5, 25.3 and 29.3 A against SAXS's oblate 15.7-16.1 and 27.9-28.3 A (Oliver et al.
  2013, Table 2), 99.5% of maltose beads outside the core.
* **Around glpG the heads do not stay out**: 200 DDM around the prep's oriented CG glpG, held frozen, all bind (128 in
  the belt, 64 on the +z face: L1 loop 47-63, W130, F176-G180), tails within 0.94 nm of 95% of exposed belt beads, but a
  maltose bead is within 0.6 nm of 56% of them at |z| < 0.5 nm, 67% at 0.5-1.0 nm and 84% at 1.0-1.41 nm; A37, V41,
  V122, A190, V194, A197, M198 and V201 are touched by heads only. POPE/POPG puts no head at the belt core. TM4 is
  tail-covered (W146).
* Whether the coat comes from this protein surface in GROMACS's Martini 2.2 glpG or from 200 molecules exceeding what the
  belt holds is not measured. Production does not use that representation: the hybrid's protein meets the environment
  through its BB bead and the side-chain/lipid table, which covers every DDM bead type (C1, Na, P1, P2, P4 among its 13
  targets), so the arrangement that matters is the one Upside relaxes to.

## 8. Bilayer physics measured with this force field

**Reference (GROMACS 2026.4, CHARMM-GUI's own protocol, 2026-10-20).** The published potential (RF eps_r 15,
potential-shifted LJ, both 1.2 nm), SD tau-t 4 ps, C-rescale semi-isotropic 0 bar, dt 40 fs, second half of 50-100 ns
unrestrained: DOPC (1148) 64.8 A^2 / P-P 43.8 A with 0.15 M NaCl, 62.8 / 44.9 without; POPE/POPG (332) 61.8 / 41.8
with or without salt. CHARMM-GUI's mdp cutoff (1.1 nm) adds 1.9-2.4 A^2; the Arnarez 2015 force switch (0.9-1.2 nm)
adds 3 A^2. About 27-30% of tail ends lie past the midplane at equilibrium (the bilayer_ts.py definition); CHARMM-GUI's
step5 bilayers start at 83-93%, so GROMACS's own run jumps to 88 A^2 (DOPC) within 50 ps and needs about 12 ns (2 ns
without the PO4 restraint) to recover. Frames and scripts: `scratchpad/bilayer_ref/gmx/`.

**Upside against GROMACS from the same start (2026-10-21, `scratchpad/bilayer_ref/upside_from_gmx_eq/`).** Each
GROMACS-equilibrated frame (1.2 nm, no added salt) run in Upside on the same potential as tables, tensionless xy MC
barostat, 400,000 steps (0.73 ns), second half against GROMACS's second half: DOPC (288) 63.26 +- 0.11 A^2, P-P 44.66 A,
27.0% of tail ends past the midplane, against 62.88 +- 0.22, 44.88, 26.7%; POPE/POPG (332, counterions) 60.93 +- 0.12,
42.24, 27.5%, against 61.74 +- 0.12, 41.84, 27.0%. DOPC matches. POPE/POPG is 1.3% smaller and 0.4 A thicker, its area
still falling at the end (61.1 at 300,000 steps, 60.4 at 400,000); the block error ignores correlation beyond 0.1 ns,
so 0.73 ns cannot separate a difference in the potentials from a slow fluctuation. The minimisation before MD thins
both bilayers by 3.4-3.6 A in P-P (POPE/POPG 41.8 to 38.3), recovered within 8,000 steps; cause not examined.

**Chain order** (g-JF, 128-DOPC, old barostat, whose area and thickness the reference supersedes): P2 tail-average 0.30
(MARTINI DOPC 0.2-0.3), tilt 22.5 +- 12.3 deg, director S = 0.74; fluid but somewhat over-ordered. Report that, do not
twist it. Its P-P excess over experiment (37-39 A) was headgroup projection: the hydrophobic core d_c (C1A/C1B
leaflet-leaflet) was 27.8 A, matching DOPC's 2.7-2.9 nm.

**Transport is not claimable.** Diffusion, reorientation and flip-flop are cage-escape processes needing physical
friction, of which the g-JF has almost none: MSD(lag)'s local exponent falls from 0.51 to 0.27, so any match of D to
11.5 um^2/s at one lag is a coincidental crossing (an overdamped control is also sub-diffusive, alpha 0.35). **No single
effective-time factor maps CG lipid time to physical time.** Claimable: correct thermodynamics with a nominal sampling
time (what REMD needs) and correct static mechanics.

**Elastic moduli need a fluctuation-correct barostat.** The damped "Parrinello-Rahman" relaxation in `box.cpp` left the
area almost still (sd 0.95 A^2 against ~81 expected at K_A ~ 265). The Monte Carlo barostat (exact-NPT Metropolis
scaling molecule COMs via `/input/molecule_ids`, so dU is intermolecular) restores fluctuations (sd 29), but K_A is
still undersampled because area relaxation is slow. Mean APL under a barostat equals the NVT APL.

Engine facts: REMD is correct for these lipids (per-slot T reaches the lipid integrator at the thermostat cadence, the
lipid potential is in the swap criterion, no momentum rescaling is needed under a coordinate swap, Arrhenius gamma(T)
only sets timescale), stable from T_up 0.5 to 1.5. A latent hazard, flagged and not fixed: the `mv` RESPA overload does
not call `apply_brownian_step` nor skip `brownian_mask`, so `mv` plus brownian lipids would mis-integrate them; MARTINI
uses `v`. `effective_time_factor` is never read by the engine (analysis-only metadata).

### 8a. The annular lipid shell is still filling for the first ~2/3 of the glpG production run (2026-09-14)

The gap around the protein closing over the trajectory is post-insertion relaxation of the boundary lipids
(`glpG_RKRK_79HIS_run0_remd.vtf`, 3146 frames, cold rung T = 0.70, 456k steps at dt = 0.009). It is not the box
(production is fixed-volume, 99.768 x 99.768 x 180 A) nor a pore (lipid-free area outside the protein is 0 A^2 in every
frame at a 5 A probe), nor the protein (TM-slab Rg_xy flat after the first eighth). The lipids move inward (every bin
inside 20 A of the protein surface gains, beyond 25 A loses) in two stages: the number of annular lipids saturates
early (36 -> ~44) while contacts per lipid climb 3.8 -> 4.8 (TM-core contact beads 137 -> 220). A - B exp(-t/tau) gives
**tau = 1190 t_up, 29% of the run**; the reverse cumulative mean settles within 1% only over the last 20-30% (last 10%
262.2 against 223.6 for the whole run). Do not convert tau to real time through the unit table: `dt` is locked to
`/input/brownian`, whose friction is tuned for a target lipid diffusion not verified here.

**Open, not acted on**: the first two-thirds are not equilibrated in protein-lipid packing, the interface the HDX
protection / membrane-accessibility estimator reads, and Campaign 6 used all frames. Before calling the four-variant
equality converged, re-run the estimator on the last third and measure this shell curve for the other three variants.

### 8b. The stage-7.0 gate counts glpG's bundle interior as vacuum; the lipid-facing gap is real and heals (2026-10-21)

The old gate (`assert_environment_solvation`) failed a build if any belt BB site (|z - midplane| <= 14.1 A for glpG) had
no environment bead within 2 x `protein_lipid_cutoff` = 6 A. Measured on its own quantity, nearest environment bead per
belt site, split by burial (protein N, CA, C, O within 10 A; `scratchpad/glpg_belt_20261021/gate_calibration.py`,
`burial_split.py`):

| state | all belt: median / >8 A | lipid-facing (burial 40-60): median / >8 A | buried (>= 80): >8 A |
|---|---|---|---|
| 1rkl production, healthy (one helix) | 5.62 / 1% | 5.64 / 1% | none |
| glpG findings-100 production, healthy | 5.2 / 0% (all BB sites, max 7.3) | | |
| glpG 7.0 state the gate rejected | 8.03 / 51% | 6.86 / 18% | 69% |
| glpG 6.6 continued in NPT, plateau 60k-300k steps | 7.78 / 45% | 6.27 / 5% | 68% |
| glpG production, 50k steps from the 7.0 state | 6.92 / 30% | 5.89 / 5% | 39% |
| glpG production, 200k steps | 6.02 / 15% | 5.56 / 4% | 25% |

* The 6 A rule fails every healthy state (38% of 1rkl's healthy sites exceed it), and most of what it counts in glpG is
  the inside of the six-helix bundle, surrounded by protein.
* The lipid face had a real gap when production started (18% beyond 8 A against 1%). Prep NPT cannot finish it: the box
  converges in 30,000 steps (99.8 to 98.5 A) and the lipid face reaches 5% by 60,000, then stays, under the rigid
  protein and the PO4 restraint (the prep's NPT had been 5 x 500 = 2,500 steps). Flexible production closes it to 4% by
  200,000 steps (the slow filling of 8a), protein intact (peptide C-N max 1.84 A, KE/1.5kT 1.003).
* With 6.6 at 60,000 steps (A102; `scratchpad/glpg_gate_20261021/`) the gate passed at 9.9% against the 10% limit and
  production started with 36-42% of all belt sites beyond 8 A, cause not measured. Replicates (`replicates/`, tables in
  `compare/`, independent draws) read 9.2% (pass) and 16.1% (fail): under the rigid protein the lipid face wanders 7-16%
  over tens of thousands of steps, so one 10,000-step window against 10% fails about one prep in three. The bare sites
  are the belt ends, 10-14 A from the midplane on both sides; the belt is centred on the bilayer (+0.3 to +0.9 A); the
  bilayer within 15 A of the protein is 1.5 A thinner (P-P) than beyond 30 A.
* Flexible production closes the face from every gate state seen, protein intact: 9.9% to 0.7%, s12 9.2% to 0.4-2.9%,
  and the refused s13 state (`s13_prod/`) 16.1% to 4.5% over 20,000 steps, still falling (C-N max 1.93 A, KE/1.5kT
  0.999). The rigid-stage gate does not predict a damaged production; it reads the belt ends' contact under a frozen
  protein, for a reason not measured.
* **So the gate measures production's steps 15,000-20,000 (A126, user 10-21).**

## 9. The nanoparticle campaign (1AO6 + MPA-AuNP)

Retired; its run directories were deleted 2026-10-11 (10.25). The pre-fix campaign's site conclusions are retracted: it
predated the LJ core table's inner knot fix (`r_min_ang` 0.00 -> 0.30, 3.1, findings 92/93) and the CB placement fix
(`[0.0000, 0.9438, 1.2068]` -> `[-0.0198, 1.5117, 1.2068]`, 2.3, findings 102), the first leaving the short range where
adsorbing particles go force-free, the second moving every sidechain site, where `martini_sc_table_1body` both drives
adsorption and scores the footprint, by 0.568 A. The rebuilt run's contacts (highest frequency 0.341, none of the
paper's five lysines) are provisional and not to be quoted. Lessons:
* **A migration script fixes what it says it fixes**: the arity migration passing had been read as the configs being
  current, and they were two generations behind.
* **A contact-frequency footprint needs a check on how much is in contact, not only where**: 105.5 of 578 residues per
  compact frame touching a 5 nm particle (16.8 after the fixes) was the protein smeared over it.
* **Choose the analysis window before looking**: once albumin spreads, contact discriminates nothing (core and surface
  0.272 vs 0.261); in adsorbed, still-compact frames (Rg < 1.25x native, the state the experiment's CD labels) it does
  (0.128 vs 0.245).
* Score contact at CB, where the table anchors (Kabsch fit of the stored `affine_alignment` reference onto N/CA/C plus
  the fixed CB offset matches the engine to 1.5e-4 A), and rank lysine against lysine. Replicates beat length (the
  compact fraction fell 3.2% -> 1.9% as runs grew); pick a test orientation by its historical onset (time to first
  tearing differed 16x across faces).
* Selecting `residue_ids == r` without the `particle_class == "PROTEIN"` mask picks up GOLD/MPA/ION beads sharing
  residue numbering (a protein-gold distance of exactly 0.00 A).

### 9.x The driver's logged Rg is not periodic-image corrected (measured 2026-09-18)

`np.<jobid>.out` printed Rg with no minimum-image correction, inflated 2.3-2.8x on three of six faces; run3's 332.5 A
**in a 300 A box** (118.7 minimum-image) had been recorded in `remote_jobs.md` as the headline observable. **Compare
every length observable against the box before reading it**: a shape descriptor computed per atom across a periodic
boundary fails silently, finite and smooth (the same defect as the VTF per-molecule wrap). Minimum-image Rg is itself
valid only within half a box. Use a contact count as the adsorption observable.

---

## 9c. hbond is a function of the rama COORDINATE, and its turn branch is the glycine region (2026-09-18)

`RamaMapPot` and `HBondEnergy` are sibling potentials on the same `rama` CoordNode (`src/hbond.cpp:490`): hbond never
reads the map table, so editing `rama.dat` changes no hbond parameter, but each residue's per-hbond energy depends on
its (phi,psi), `Ehbond[i] = E_alpha*helix_score + E_beta*sheet_score + E_other*turn_score`, summed as
`sum_i hb_number1[i] * Ehbond[i]`. Decoded from `hbond.h5` (directions per `compact_sigmoid`, 1 for large negative
argument, `src/vector_math.h:700`): turn = phi in (0, 165) deg -> E_other -1.769; helix = phi<0, psi in (-120, 60) ->
E_alpha -1.961; sheet = phi<0, psi outside -> E_beta -1.946; all four sharpnesses 3.81972 (a 15 deg ramp).

The turn branch is the positive-phi region, the glycine question: an H-bond at phi>0 is worth 0.192 E_up less than at
phi<0, and glycine lives there (alpha_L 30.95%, ASN second at 13.4%). The coupling runs both ways:
`rama_sens(0,i) += hb_number1[i]*dPhi[i]` pushes an H-bonded residue toward phi<0 by ~0.38 E_up for a doubly-bonded
helical glycine, against ff2.1's -1.238 E_up rama pull toward alpha_L, so rama wins by ~3x: a mechanism for glpG TM4 and
lambda's H2 (11d). FF1's trainer fitted `hb` as one scalar; of the node's 12 entries FF2's ConDiv fitted 0-3 (branch
energies, second-H-bond term), while the boundaries 4-11 are hand-set and have no engine derivative (9e, 9v). `E_other`
and glycine's alpha_L map depth are near-degenerate (both set what an H-bonded glycine pays at phi > 0) and the energies
are shared by all 20 types, so ff3.0 gives glycine its own offsets on them (1.15, 1.17).

---

## 9d. The retired ff3.0's benchmark split by native/de novo, not by topology (2026-09-18/19)

The retired ff3.0 (FF1-form trainer, every glycine map symmetrised) against FF2 on the 32 Peng arms
(`0914/figs/ff2_curves_s5.npz`): native +0.039 mean TM, de novo -0.030, native minus de novo positive in 13 of 16 pairs
(sign test p = 0.021); regressions span every topology. **Do not read a partially scored benchmark**: p went 0.039
(n 9), 0.092 (13), 0.035 (15), 0.021 (16).

---

## 9e. The FF1-form ConDiv port and the learned glycine map (9e-9q; 2026-09-18 to 09-24, retired)

From 09-18 to 09-24 the trainer was a Python 3/torch port of Peng's FF1 Theano ConDiv (9t), with Track A training
glycine's map inside it; both were retired when the port turned out to train FF1's Hamiltonian (9t-9v). Old 9f-9q are
merged here; the library's handedness by chiral context (old 9i) is in GLY_sym.md §2a. Facts still used:

* **H-bond parameters.** The FF1 original trained `hb` and `sheet` as single scalars (`init_param/hbond` -2.112,
  `init_param/sheet` -0.268); the port dropped both. ff2.1's branch energies and second-H-bond term (`hbond.h5` entries
  0-3) were returned by FF2's ConDiv (SI p. 2), which stopped at 76 iterations without a convergence test (1.23);
  entries 4-11 are hand-set, and the SI does not say how the 20 `sheet` values were obtained. `hbond.h5` is
  `[E_alpha, E_beta, E_other, E_bias | 8 rama boundaries and sharpnesses in radians]`, and energy is exactly linear in
  `hbond_energy.parameters[:4]` (to 5.6e-7), so `apply_param_scale(hb_scale=...)` (`--hb-scale`) scans it with a config
  rebuild only.
* **The port's removed additions**: a `+1e-12` in two direction normalisations, and a GLY palindrome on the rotamer pair
  angular profile baked into the old ff_3.0 `sidechain.h5` (blob `72ae60be`, commit `28185321`), whose removal restored
  `pack_param`'s `discrep < 1.6e-4` gate. Inherited: `training_list` is shuffled unseeded, so minibatches differ run to
  run. A fixed bug returned through a clone of the one unpatched directory: **clone from the most recently fixed
  directory, not the most successful one.**
* **Testing a fixed point (old 9k).** From ff_2.1 under the port, six minibatches gave gradients indistinguishable from
  noise in every group (sign-flip p 0.22-1.0): a stationary point, not shown to be an attractor. Use the exact sign-flip
  test on `||mean g|| / mean|g|` (`ConDiv.py gate`); the pairwise-cosine t-statistic is anti-conservative, and partial n
  misleads (the sheet gradient read t = +3.91 at n 3, -0.01 at n 6).
* **Track A, the learned map (old 9l-9q).** The 72x72 map's gradient is analytic (`rama_map_pot` is a periodic bicubic
  spline, `src/spline.cpp:262`). Its finite-difference gate first found a 37% error, the Boltzmann-weighted constant
  `write_rama_map_pot` subtracts per map (up.md 2.8a), invisible in basin differences: read such a check as an eps sweep
  and cover every branch (terminal glycines were missed by a test protein without them). No window shorter than an
  epoch means anything (memory `condiv-gly-epoch-scale-only`). It converged at dG(aR->aL) -0.885 (9s), past the training
  natives' -0.50, compensating `hbond`'s shared `E_other` penalty (9c) as the only term both residue-type specific and
  (phi,psi) resolved (architecture.md §2).
* **PDB statistics as a target (old 9o).** The 456 training natives give glycine dG(aR->aL) -0.500 +- 0.054 against
  -1.182 for NDRD's `GLY|ALL` marginal (NDRD holds loop sites only, 1.9). NDRD's map is a PMF over folded structures,
  and ff2.1 is a constrained optimum (its other terms were fitted with the map fixed).

---

## 9r. The handedness is not an artifact of one force field (ff14SB, 2026-09-19)

`gly_awh14` (10 dipeptides, 100 ns), extracted with `gmx awh -b 95000` on the last part file:

| force field | LA | LM | LP | LL | LT | LE | LV | LD | LR | LG blank, must be 0 | chiral mean (n=9) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ff14SB | -0.519 | -0.328 | -0.418 | -0.591 | -0.043 | -0.013 | -0.206 | -0.124 | +0.026 | -0.010 | -0.246 |
| ff99SB-ILDN | -0.447 | -0.301 | -0.352 | -0.396 | -0.129 | -0.090 | -0.309 | -0.179 | -0.413 | -0.029 | -0.291 |

The means agree to 0.045 nats, inside the ~20% uncertainty, and both blanks sit on zero: glycine's left-handed bias is
not an artifact of `amber99sb-ildn`, and both are far from the library's -1.24 and ff3.0's exact 0. Per-system values do
not agree (LR +0.026 against -0.413): the mean reproduces, the neighbour structure does not. Both are AMBER, so this
bounds the within-family systematic only; CHARMM36m, the stronger test, has not been run.

## 9s. Track A's end point, and rama31's handedness (2026-09-24)

Track A (9e) ended at step 500 converged, not step-limited, with `X|GLY` dG(aR->aL) -0.885, drifted from ff2.1 but not
toward the AWH map (populated-region correlation 0.699 -> 0.703). `parameters/common/rama31.dat`, built from the AWH
data, gives -0.154 with both blanks excluded (-0.150 when `build_rama_from_awh.py` counted `RG`, the achiral
Ac-Gly-Gly-NHMe, as chiral). **The like-for-like AWH target for Upside's left/right mixture is the context-averaged map
(-0.15)**, not left plus right (-0.31), which assumes neighbour effects add in a tripeptide (untested). The current
library (10-01 build, `build_gly_library.py`): GLY_sym.md §5.

## 9t. The trainer omits FF2's backbone desolvation term entirely (2026-09-24)

The retired FF1-form port never built `bb_sigmoid_coupling_environment` (`main_worker` passes `environment_potential`
but no `bb_environment_potential`), so `bb_env.dat` came out byte-identical to ff2.1's, while every deployment includes
the term. The Theano original (`~/Documents/ConDiv/remd-4000-8RP-1th-test/ConDiv_original.py`, Upside 18.10.08) is
FF1-era (`Update` is `env cov rot hyd hb sheet`); the Peng 2022 SI made FF2 "by adding an explicit backbone desolvation
term", with "all parameters ... optimized simultaneously". The term is large and favours the unfolded state: native
ubiquitin, total -148.6, the term -16.8, as large as rama (-17.0), and -43.8 for the same coordinates expanded 1.6x; with
`scale` -0.30 it pays for each solvent-exposed backbone NH/CO. Its `hbond_weight` couples it to the H-bonds and it copies
`environment.h5`'s `weights`, so the port fitted `hb` and `env` without it. FF2's dual objective
`d alpha = d alpha_NSE + lambda * d alpha_DSE` was missing too. The replacement trainer is 9v.

## 9u. The trainer uses FF1's burial function; ff2.1 as published uses FF2's (2026-09-24)

`environment.h5` holds two burial functions: `energies` (20x18 spline, `--environment-potential-type=0`, node
`nonlinear_coupling_environment`) and `scale`/`center`/`sharpness` (sigmoid, type 1, `sigmoid_coupling_environment`).
`upside_config` defaults to 1, every example leaves it there, and the FF2 SI says the spline "is replaced by the
sigmoid-like function"; the FF1-form ConDiv hard-coded 0. Native ubiquitin, same coordinates:

| side-chain burial | native | expanded 1.6x | native minus expanded |
|---|---|---|---|
| ff_2.1, type 1 (sigmoid, as published) | -47.49 | -26.90 | -20.6 |
| ff_2.1, type 0 (spline) | -12.88 | -12.25 | -0.6 |
| ff_3.0 trained, type 0 | -15.75 | -16.38 | +0.6 |

ff2.1's spline is a vestige. So the port trained FF1's form from a start that was not ff2.1, the trained ff_3.0 at the
default type 1 silently gives ff2.1's burial, and the ff3.0-vs-ff2.1 benchmark compared two burial forms (`bench_run.py`
set type 0 for ff3.0). `SigmoidCoupling::get_param_deriv` covers scale, center and sharpness per type; the backbone
term's covers only `scale` (9v). `Train(1).zip` (OneDrive) is Peng's 2022 membrane-potential trainer (`UpdateBase` is
`cb icb hb ihb`, soluble `ff_2.2` fixed); it confirms FF2-era runs used `environment_type = 1` with `bb_environment` on.

## 9v. The FF2 trainer: found, adapted, and what the port had wrong (2026-09-24)

The only FF2 dual-target trainer is O. Kleinmann's Python 3 port of Peng's code,
`/project2/trsosnic/okleinmann/condiv/condiv2.py` (git history from 2025-08, first commit already his working copy, so
Peng's pristine file is not recoverable; `/home/pengxd` is unreadable). Everything else searched is FF1 or membrane:
`~/Documents/ConDiv` (9t), `~/Documents/Train` = `Train(1).zip` (9u), `upside_version/upside-pxd/ConDiv` (2019, spline
burial, no DSE). His port had drifted from the SI, all corrected in `training/ConDiv.py`: lambda = 0.0
(`balance_target`), so his run never used the DSE objective; 6 free replicas up to T ~0.97 (SI: 12 from 0.8 to 1.1),
1000 time units (SI: 8000), minibatch 21 (SI: 24); replica reweighting exponent `E*(T0-Ti)/Ti`, T0 times the correct
`E*(1/Ti - 1/T0)`; a `dE < -200` clamp in place of a normalisation; and a guard that silently dropped the DSE term
whenever the last free replica's final energy exceeded 1000. His 101-step run ended with the backbone scale flipped
from -0.30 to +0.12; its negative PRO burial sharpness is ff2.1's own value (-0.28), not his drift.

Engine limits shared by ff2.1's training: `BackboneSigmoidCoupling::get_param_deriv` computes only the `scale`
derivative (the other three are commented out, in master too) and `HBondEnergy::get_param_deriv` only entries 0-3, so
ff2.1's workflow never trained the backbone term's center, sharpness and hbond weight, nor the eight H-bond rama
boundaries. The user chose to keep that exactly; the smoke worker confirmed those three contrasts are exactly 0.

Setup as validated on midway2: 19 x 24 minibatches as in the SI; ff2.1 starts at the SI's H-bond energies
(-1.961/-1.946/-1.769; second-H-bond -0.406). The smoke worker (2xf6, 600 time units) had all 14 groups finite and the
SARW replica at 0 H-bonds and Rg 24 A. The local Mac binary cannot run workers: it traps at exit whenever MC pivot moves
are on (README trap).

## 9w. Phase 2's full-map glycine row: steady drift, then cancelled (2026-09-25 to 09-28)

The full-map glycine row (plan.md Phase 2) drifted monotonically from ff2.1 (~0.002 nats of dG per step) without
approaching the AWH map, its gate failing at every checkpoint; cancelled 09-28. **To tell a step-size limit from noise,
use two numbers over an epoch**: Adam utilisation (rms step / alpha; pure noise gives sqrt((1-b1)/(1+b1)) = 0.33) and
the sign consistency of each cell's steps (noise 1/sqrt(n)). Utilisation at the noise floor with consistent signs
(0.37-0.42 and 0.59-0.70 here) is a weak steady pull whose drift scales with alpha; at Track A's end both were at noise,
and no step size would have helped.

---

## 10. Cluster and operational lessons

### 10.1 A nested sbatch inherits `SLURM_*` from the job that calls it (2026-09-06)

A submitting script that requests `--mem` passes `SLURM_MEM_PER_NODE` to its child; if the child's `srun` requests
`--mem-per-cpu`, every launch dies (`SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually
exclusive`). It broke a training chain's only resubmit path. **Any script that submits another job unsets those three
variables after its `#SBATCH` block, or matches the child's memory-request type**; this does not change its own
allocation.

### 10.2 The batch script is fixed at submission: edits do not reach a queued successor, and a requeue reruns the old one (2026-09-23)

Slurm copies the batch script into the job record at submit time, so a self-chaining job's queued successor runs the
old text. Read it back with `scontrol write batch_script <jobid> /tmp/js.sh && grep -n <the-thing-you-fixed> /tmp/js.sh`,
then cancel and resubmit the successor with the same `--dependency`. **Cancel the pending successor while the running
link keeps going**; the reverse order lets it start on the dying parent's `run_output` (two concurrent writers, once).
* **A branch that runs only on success, only at the end, is untested by construction**: `train_gly.sbatch`'s handoff
  `sbatch "$T/validate_ff31.sbatch"` pointed at the run directory, not `training/`, after twelve links had exercised
  every other line. Verify the terminal branch by hand before the run that reaches it.
* **A requeue reruns the script verbatim with its original arguments**: a requeued link resumed from its stale
  checkpoint and `main_loop`'s `rmtree` deleted the newer ones. Chain scripts carry `#SBATCH --no-requeue`
  (`scontrol show job <id> | grep Requeue` reads `Requeue=0`) and resolve the newest checkpoint on disk at run time.

### 10.3 A Slurm job can fail with no error text at all (2026-09-23)

A training link died at 8 h 36 m with exit 7, no traceback and every `*.output_worker` at 0 bytes. Step level showed it:
`sacct -j <id> -o JobID%16,State%16,ExitCode,MaxRSS,NodeList%30 -P | awk -F'|' 'NR==1 || $3!="0:0"'` gave twelve steps
`CANCELLED 0:7`, SIGBUS on four nodes at once (`ExitCode` is `exit:signal`, so `0:7` on job steps is the informative
one). The successor's batch step was cancelled the second it started, before opening its `.out`: those nodes could
neither read mmap'd files on `/project` nor create one (inferred, no GPFS log; disk capacity measured and not the cause).
* **When a job log ends mid-sentence, go to `sacct` step level before reading the code**: a dying parent's buffered
  stdout is lost, so the log's last line is where the buffer last flushed.
* **Check the physics at the moment of death before assuming a blow-up** (here every worker was healthy, Rg 14.5 A,
  potential near -200), which points at infrastructure. A successor that dies before opening its output ran no line.

### 10.4 A target checked only at restart does not bound anything (2026-09-23)

`train_gly.sbatch` tested `STEP >= TARGET` at link start, then ran `STEPS_PER_LINK=150` regardless: resuming at 491 of
500 would have run to 641. **A self-chaining job applies the bound to the work, not only to the decision to start it**
(`REMAINING=$(( TARGET - STEP ))` clamped against `STEPS_PER_LINK`). **Delete a partially written output directory
before resuming**: `main_worker` reads `<name>.divergence.pkl` by path with no freshness check.

### 10.5 A partition that accepts a job is not a partition the job may use (2026-10-01)

User correction. With no CPU allocation after the 10-01 rollover, midway3's `amd` partition still accepted the training
job (`--test-only` PASSED); the user cancelled it after 14 min (~76 core-hours): a CPU job must not run on the group's
GPU allocation. Rules:
* Run CPU work only on broadwl (midway2) or caslake (midway3) unless the user names another partition.
* When the usual partitions refuse a job, stop and report. Do not look for one that accepts it.
* Before any first use of a partition, find out which allocation it bills, and ask.

### 10.6 A heavy job on a login node kicks us off it, and every return costs a Duo push (2026-10-02)

A lambda decomposition on `midway2-login1` with an indexing bug allocated ~60 GB; the node's memory policing killed it
and the user's sessions with it, and the socket dropped about five times, each return a Duo push. **Analysis more than
a quick read goes to a compute node; anything on a login node is capped with `ulimit -v` first (`ulimit -v 6000000`).
When the socket keeps dropping, look for our own load (`ps -u yinhanw --sort=-rss`) before suspecting the network.**

### 10.7 Clean up local compute; never leave it running without a live reason (2026-09-16)

User correction, twice: 6 GROMACS replicas ran ~1024% CPU on the laptop for two hours unflagged, and were then left
running after the work moved to midway2 (the sampling done was checkpointed and the cluster job rebuilt from scratch).
* **Say the cost up front**: core count, expected wall time, that it will load the machine. The user cannot see `ps`.
* **Kill it the moment its reason expires** (work moved, better path chosen, the user wants the machine). Data already
  written is not at risk from stopping.
* **Stop cleanly so it is resumable**: `kill -TERM` makes GROMACS checkpoint and exit; `mdrun -cpi prod.cpt` resumes.
  Verify the `.cpt` exists before reporting it stopped.
* **Audit at the end of any session that launched local compute**: no stray `mdrun`/`upside`/python workers, orphaned
  launchers or watcher loops past their target; GROMACS `#backup#` files and minimization `.trr` removed.
* Background watchers carry a bounded iteration count so they expire on their own.

### 10.8 Keep separate problems separate; do not carry an unproven one along (2026-10-01)

User correction. I presented the glycine fix and a pre-proline rule as one plan; the glycine problem is established,
the pre-proline one was my suggestion with a mechanism but no measured consequence (and 1.14 had found the gap to be
composition). Rules:
* Propose each problem on its own evidence, with its own plan and validation.
* Before acting on a suspected defect, state the observable it is supposed to break and whether that has been measured.
  A mechanism without a measured consequence is a hypothesis.
* Say who raised an item when it re-enters a plan, so an AI suggestion is not mistaken for an established issue.

### 10.9 Three analysis lessons from the lambda diagnosis (2026-09-18)

* **A reweighting is only as meaningful as the stationarity of the ensemble it reweights.** Reweighting lambda's whole
  native arm gave 0.55 kT of near-native stabilisation, but the arm is a decay from its seed (6.41 -> 10.39 A by block);
  on the equilibrated last third the effect is -0.024 E_up and the correlation with CA-RMSD flips sign (+0.171 ->
  -0.064), at an ESS of 83% in both. **Block the observable against time before reweighting, and quote the converged
  window.**
* **Pooling shells can invert a conclusion when one shell dominates**: pooled under 6 A, lambda's near-native states
  had H2 broken; by shell, 0-5 A has H2 intact (89% right-handed) and the four-times-larger 5-6 A shell is broken. It
  recurred 09-30, when the pre-proline class's +0.023 alpha_R gap resolved by native basin into composition (87%
  extended; 1.14). **Resolve by bin before reading a conditional average; in ConDiv a class gap measures the class's
  native composition until shown otherwise.**
* **Never import a module whose top level does work**: importing `score_arms_dist.py` for two helpers started a 2.7 h
  rescore on the login node and truncated `score_arms.json` from 23 arms to 3. Duplicate the few constants and helpers
  rather than importing a script, and check what a module does at import.

### 10.10 A checkpoint written under NumPy 2 does not load under NumPy 1 (2026-10-01)

The Mac's NumPy 2.4.4 pickles through `numpy._core.*`, so the cluster's 1.23.5 fails on every checkpoint, solver
state, divergence and rmsd file (`No module named 'numpy._core'`); `env.sh` promises identical versions only between
midway2 and midway3. `move_run.py` converts a copied `run_output` on the machine that continues the run (maps
`numpy._core` to `numpy.core` under 1.x, rewrites absolute paths, re-pickles natively; NumPy 2 reads NumPy 1 pickles
unchanged). The re-pickled classes are `ConDiv.Update`, so they unpickle only where `import ConDiv` finds that run's
own copy (`run_output/ConDiv.py gate|extract`), not from a script with `training/` on `PYTHONPATH`.

### 10.11 A deploy must carry file modes, and `cp` onto an existing file does not (2026-10-02)

A deploy left `py/upside_config.py` without its executable bit (the staged copy was first made by a shell redirect, mode
644, and `cp` over it kept 644), and every midway2 config build failed with `PermissionError` until restored. **Stage
with `cp -p` into a fresh directory, and end a deploy with `chmod --reference=<file>.bak_<tag> <file>` and a listing of
the modes.**

### 10.12 Under -ffast-math, where a new loop sits can change an old sum's last bit (2026-10-02)

Both CMake files build with `-ffast-math`, which lets the compiler reorder a float reduction. A new loop in
`HBondEnergy::get_param_deriv` right after the shared `dparams[0..2]` loop changed those sums in the last bit
(13.922128 to 13.922127) on a config that never enters it; moved after the last shared accumulation, everything was
bitwise equal. **For a parity change in the engine, compare bitwise, and when only a reduction differs, look at code
placement before the arithmetic.** Where no placement keeps the bits (A82), the proof is a null build, the same code
with a condition that never fires: if the fix equals it bitwise, every difference from the old record is code
generation.

### 10.13 Test whether a fix holds, and put the fix where the cause is (2026-10-04)

The helical-glycine cause is settled: glycine's deeper alpha_L basin in the PDB map is selection, which Upside applies
as energy (round 1, symmetrised map, fixed TM4 but is wrong physics; round 2, trained basin depth, was pushed back by
the data; round 3, physics map plus glycine H-bond offsets, is ff30_glyhb). User corrections the same day:
* **Before designing a diagnostic, read the memory and findings for the established mechanism, and aim the test at
  whether the current fix holds** (I had proposed runs to find out why the glycines flip).
* **A fix belongs to the term that carries the defect. Do not route a map problem through another term's parameters**:
  round 3 did, and the trainer used those parameters to put the alpha_L preference back.

### 10.14 Other cluster and tooling lessons

Shell and ssh:
* **`pgrep -f` / `pkill -f` match every process whose arguments hold the pattern**, including the loop or ssh command
  running them (10-05, 10-07, 10-24): `while pgrep -f "<pattern>"; do sleep; done` never ends when its own shell, or
  another waiter naming the same files, carries the pattern, and `pkill -f "while pgrep"` killed the ssh session issuing
  it (exit 255). Wait on a PID (`wait`, `kill -0 $pid`) or on the output files, kill by PID, and match with a bracket
  pattern (`ps -u $USER -o pid,args | grep "[p]ack_predictions"`).
* **`timeout` around `ssh <host> <cmd>` stops only the local client (10-23)**; the remote command runs on. Check `ps` on
  that host and stop the leftover before retrying.
* **A collector under `set -eo pipefail` dies on a pipeline whose last stage reports nothing (10-23)**:
  `x=$(ls | while read d; do [ -f $d/f ] && basename $d; done | tail -1)` returns 1 when the final directory lacks the
  file, as does `grep -l ... | sed` on no match (`pooled_collect`). Count in the loop itself, and give `grep` an
  explicit `|| true` where no match is a valid answer.
* **The Mac's zsh (10-25)**: never name a variable `path` (tied to `PATH`; every `scp` became "command not found"); zsh
  does not split an unquoted `$C`, so write an ssh command out; macOS `paste -s` needs a file argument
  (`paste -sd' ' -`); `ssh host "md5sum $files"` with `$files` empty waits on stdin, so give remote commands `< /dev/null`.
* **Check a file for a credential by count, never by printing lines (10-26)**: a filtered `grep -n "password"` printed
  the inline password of `scratchpad/mdw2_master.exp` into the transcript. Use `grep -c`, `grep -l`, or
  `grep -F -f <(extract the secret) | wc -l`, and never copy the `.exp` scripts anywhere.
* **Never pipe a script that performs writes into `head`**: `reseed.py <16 files> | head -3` SIGPIPEd the producer after
  three replicas. Use `tail`, which drains its input.
* **A poll that looks back a fixed window misses what ended while polls were skipped (10-25)**: an 80-minute look-back
  hid ff3.2.1's NODE_FAIL for four hours. Sweep from the last recorded pass (`sacct -S <last pass>`).
* Removing Python functions programmatically: a "def line -> next column-0 line" scan mis-cuts multi-line signatures.
  Use `ast` (`node.lineno..node.end_lineno`), or `py_compile` after each bulk removal.

Slurm and the two clusters:
* **`sbatch --test-only` counts against the gpu QOS's 12-job limit (10-24)**; test a new job's shape when a slot is free.
* **caslake rejects an `--exclude` naming midway2 nodes (10-05)**, so a self-resubmitting script that runs there passes
  its job's own `ExcNodeList` (`ff3_benchmark/bench.sbatch`). Check a shape on the other cluster with `--test-only`.
* **midway2's modules are not on midway3 (10-05)**: `python/3.9.18` and `openmpi/4.1.1+gcc-10.1.0` are absent
  (`hdf5/1.14.3+oneapi-2023.1` is on both); caslake jobs take python and hdf5 from `env_shared.sh`.
* **Compare clusters with `--test-only` of the same shape, never a held job's StartTime (10-05)**: a held job is not
  scheduled. caslake has no per-CPU memory default (`DefMemPerNode=UNLIMITED`; broadwl `DefMemPerCPU=2048`), so a job
  that sets no memory (`panel.sbatch`) is given `--mem` there.
* **A NODE_FAIL requeue silently ate the REMD block budget (findings 125)**: one job requeued 5 times on an
  `ALLOCATED+NOT_RESPONDING` node reached block 6/12 with zero chunks. `squeue` shows `RUNNING`, `sacct -X` only the
  newest incarnation (use `sacct -j <id> --duplicates` or `scontrol show job <id> | grep Restarts`), and the log is
  truncated per requeue; `grep -c "chunk done"` is the honest progress metric. Since 2026-09-30 `run_remd.py` counts a
  block per job id (`<V>/block_jobid`). Rules: reset `block_count` to the completed blocks after a requeue; exclude a
  failed node through the one exclusion law (CLAUDE.md, `remote_jobs.md` §0d), or Slurm re-allocates it indefinitely;
  treat `NOT_RESPONDING` as unproven death, confirming static log size and h5 mtimes over a dwell and finite
  `input/pos` in every replica before resubmitting. An abrupt kill loses only unflushed HDF5 buffers, and `reseed()`
  resumes from the last consistent state.
* **`~/cds3` (`/cds3/trsosnic/yinhan`) is unreadable from compute nodes**; stage to `/project` (`~/project` works).

Storage incidents:
* **A wedged GPFS makes a dead job look healthy (2026-09-07)**: `RUNNING` for 3.5 h with every worker in `D` state at
  `00:00:00` CPU (wchan `cxiWaitEventWait` / `lookup_slow`). Probe per process: `sstat -a -j <id>` (an `AveCPU` that
  does not climb is not computing), then `ps -o pid,stat,time,etime,wchan` on the nodes. `timeout 10 ls <wedged dir>`
  does not return from `D` state; over a ControlMaster the leaked processes burn channels until ssh falls through to
  password auth, the RCC-ban trigger. Pin `-o BatchMode=yes -o PasswordAuthentication=no -o NumberOfPasswordPrompts=0`
  on every cluster call before probing anything that might hang.
* **RCC's GPFS serves midway2, midway3 and beagle3, so the other cluster is no fallback in a storage incident** (09-07:
  midway3's logins had lost `/home`, `/project`, `/project2`, `/scratch` and `/software`; `stat -f /project` read xfs).
  Tell a stalled mount from a lost filesystem with `stat -f`, and check a second login node.
* **All four `midway3-login[1-4]` share `midway3.rcc.uchicago.edu`'s host key**: by per-node hostname an expect script
  answers the host-key prompt with the password and no push is sent, which looks like a Duo failure. Pass
  `-o HostKeyAlias=midway3.rcc.uchicago.edu` (`scratchpad/rcc_master.exp` does).

Environments and builds:
* **Upside has no MPI (2026-10-02)**: `mpirun -np 2 upside ... a.up b.up` starts two copies of the whole job (HDF5 lock
  errors, 0 frames, exit 0); parallel systems run as OpenMP threads in one process (`OMP_NUM_THREADS=2`). The binary
  loads `libmpi.so.12` only because `hdf5/1.14.3+oneapi-2023.1` is a Parallel HDF5 build.
* **The `/beagle3` deployment and `$P` give bitwise the same energy and forces for a 12-entry force field (10-05)**
  (lambda native E -199.4555816650, max |dF| 0; `checks/r4val_20261005/engine_parity/`); a config with glycine H-bond
  offsets needs the deployment synced.
* **The one midway2 binary gives the same forces on midway3 but not the same trajectory (10-05)**: energy, forces and
  derivatives are bitwise equal (`checks/mdw3_gdepth_env_20261005`), but a run first differs at step 6 by 3.7e-9 A and
  decorrelates by ~180; within a cluster runs are bitwise. The difference follows the OS image (el7 against el8). So a
  binary-run parity test is valid only on one cluster, and a run that changes cluster is a statistically equivalent
  continuation.
* **`extract_ff.py` took the trainer from two directories above the checkpoint (10-05)**, wrong for
  `run_output/initial_checkpoint.pkl`, where it silently imported `training/ConDiv.py`; it now uses the checkpoint's own
  directory when that holds ConDiv.py. Test an extraction on the initial checkpoint of every new kind of run.
* **`/input/potential/backbone_pairs/ref_pos` holds NaN by design** (glycine's CB, `write_backbone_pair`), so compare
  float arrays with `equal_nan=True`. With that, `patch_glpg.py` output does not depend on which force field the seed
  carried (all 191 `/input` nodes identical).
* **Build Python environments on `/beagle3`, not `/project` (10-21)**: 200 `touch` plus `rm` took 26.7 s on
  `/project/trsosnic` against 0.45 s on `/beagle3/trsosnic`, and a torch venv is tens of thousands of files.
* **midway3's `python/3.13.11+gcc-12.2.0` module has no `_ctypes` (10-24)**, so torch fails at import. Above 3.11 use
  uv's managed CPython (`UV_PYTHON_INSTALL_DIR` on shared storage), as the BioEmu2 venv does; `python/3.11.9` is sound.
* **`pip install torch` from PyPI gives a CUDA 13 build (10-21)** (`torch 2.14.1+cu130`, sm_75 to sm_120), which cannot
  run on midway3's V100s (sm_70). Install from `https://download.pytorch.org/whl/cu126`
  (`/beagle3/trsosnic/yinhan/boltz_hdx/install.sh`).
* **Boltz-2 at `b1ebfc4` ignores `--max_parallel_samples` and exits 0 when it runs out of memory (10-22)**
  (`diffusionv2.py:383` made one chunk whenever the sample count is a multiple of the setting; Boltz then logs "Number
  of failed examples: 1"). Our clone `$B/boltz` chunks by `max_parallel_samples` (no backbone atom moves by more than
  0.001 A at a fixed seed), and the wrappers fail a run whose log lacks "Number of failed examples: 0".
* **Boltz-2's confidence step bounds the samples one run can hold (10-23)**: `confidencev2.py:121-143` keeps every
  sample's L x L x 64 float32 PAE and PDE logits on the GPU and concatenates them. Budget about 4 x N x L^2 x 256 bytes
  per run and split larger sets across seeds.
* **An isolated Upside engine without touching the repo's `obj/` (10-21)**: `py/upside_engine.py` loads
  `<its own directory>/../obj/libupside.*`, so a copied `py/` beside a fresh `obj/` (cmake on `src/` with
  `install_M1.sh`'s flags) and a symlinked `parameters/`, with `UPSIDE_HOME` pointed there, is a complete engine; a
  symlinked `py/` does not work. Example `scratchpad/boltz_hdx/flags/upside_tree/`.
* **SciencePlots' `science` style sets `savefig.bbox: tight`, which crops a fixed page size (10-09)**; set
  `'savefig.bbox': 'standard'` for a print page.

Testing and run hygiene:
* **A single-replica smoke test bounds loading and setup, not stability**: a rare force spike needs ~1e4 steps across
  48 replicas to appear.
* **Compute the expected event count before an A/B on a rare stochastic failure**: one reaction-field test was 10x
  under-exposed (~2 events per 2e6 replica-steps against 1.9e5 per arm), another confounded by structural relaxation.
* **`--thermostat-interval -1` does not mean NVE**: `max(1., round(arg / (inner_step*dt)))` clamps it to 1, so the
  thermostat fires every step.
* **Recovery is not a guard**: reseeding a destroyed replica between rounds clamps nothing, widens no threshold and
  drops the destroyed frames, as `run_remd.py`'s rollback does.
* **Size trajectory output by its accumulation rate, not the seed**: a 2.4 TB balloon (2026-08-02) came from ~3000
  frames per chunk with momentum over ~43 chunks. Use `frame_interval` near 100 (~2000 frames), no `--record-momentum`,
  and one `--duration`.
* Check the units and axis directions of any workflow figure before showing it (the `_DG_Hbond.png` scale, 3.7, the
  ESS-censoring confusion and `_Tm_curve.png`'s inverted hydrogen-bond axis were found while building one poster).

### 10.15 A results slide states the comparison, not the bookkeeping (2026-10-05)

User corrections on the BP-validation slides (an 11-panel figure, a named worst frame and per-seed orderings, and a
"does the fixed engine reproduce" framing when the decision was whether to retrain ff2.1). **Frame the slides on the
decision they serve and lead with the evidence that bears on it most directly (for retraining, the training gradient).
A slide carries one comparison at the level of the whole set, the effect against its reference scale (energy error
against kT, T_mid shift against the seed spread); single frames and per-seed orderings stay off.** When the result is
not what the user expects (WW domain's 0.8 K shift against its 0.3 K seed spread), state it as measured and say so.

### 10.16 No smoke jobs for a pipeline that has already run in production (2026-10-05)

User correction: I submitted two short caslake test jobs for the round-4 validation, though the same pipelines had run
in production. **Verify a change to a proven pipeline with checks that submit nothing: `bash -n`, stubbed sandbox runs
of the changed branches, `--test-only` of every submission and comparisons of the files it writes. Do not spend
allocation on test jobs. If a real run on a new cluster looks needed, say why and ask first.**

### 10.17 ff3.0 is a controlled comparison: same data, same workflow (2026-10-06)

User correction: asked about TM4's drift, I suggested a denatured-state weight and noted that no membrane protein is in
the training set. **Never propose new training data, a different contrast weight or any other workflow change as a fix;
the only design change is glycine's (plan.md Project Goal). Attribute a drift by comparing against the same workflow
run without the glycine change (`ff21-fixedpoint`), not by retuning the workflow.** The control answers whether training puts
the selection back into other terms.

### 10.18 Re-read the queue before a decision that depends on job state (2026-10-07)

I called a control "pending, never run" from a `--test-only` projection 23 h out; it had started five minutes after
submission. **A `--test-only` start time is a projection, often a day too pessimistic on broadwl. Run `squeue`
immediately before presenting any choice that depends on whether a job is running, and again before acting on it.**

### 10.19 TM4 is judged only in the hybrid, with the in-training force field (2026-10-07)

User correction, repeating 10.17's mistake: I answered with the soluble panel and the H-bond margin. **For any TM4
question, the measurement is the 12-seed hybrid glpG test of the checkpoint against its own start (`patch_glpg.py` +
`run_glpg.sh` + `tm4_local.py`). Do not offer panel helix, margins or the training set's composition as TM4 evidence or
as a reason to discount a force field. When no hybrid result exists yet, say so and say when it will.**

### 10.20 A helix test scores secondary structure, not only dihedrals (2026-10-09)

User instruction (10-09): the TM4 test measures TM4's secondary structure. 1.24's phi/psi box accepts pi-bulges and
frayed residues; by DSSP, TM4 holds 0.10-0.22 less alpha-helix in every set (1.25). On 10-15 I still reported "TM4 holds,
ff3.1 0.97 against ff2.1 0.98" from phi/psi-only scripts (`tm4_ff31.py` and my own); DSSP alpha on the same frames read
0.89-0.94 and put ff3.1 last (1.35).
* **When a test asks whether a helix holds, measure its H-bond pattern (DSSP on the backbone, with the model's own
  carbonyl O), report it beside any dihedral readout, and look at a structure where the two disagree before trusting
  either.**
* **Before quoting a helix result, read which criterion the script scores**, inherited scripts included; a TM4 readout
  carries DSSP alpha, and a check that lacks it gets it added before it is reported.
* **When the user points at a figure that disagrees, first rebuild that figure from its own files** (on 10-15 that
  showed the disagreement was the old table).
* **A DSSP pi label is not a gap**: check the i->i+4 and i->i+5 O...N distances before reading a white band as a helix
  break (TM2's 92-96, 4.6).
* **Compare a figure with a reference only under the reference's criterion, and name the criterion in the companion's
  title** (a DSSP companion to a phi/psi slide differed on 31 residues and was read as a different structure).

### 10.21 Three analysis lessons from the glpG TM1 diagnosis (2026-10-09)

* **An energy-structure correlation along one trajectory is confounded with time** ("opening is rewarded by lipid" in
  one seed was tails accumulating as the run aged; time-matched frames over 12 seeds showed no difference, 4.7). Compare
  states at matched times across seeds before naming a driver.
* **Annotate from the structure we simulate, never from online data** (user, 10-09: our construct has no online
  counterpart). Derive helix ranges from DSSP on the seed's own PDB and check them in the script (11).
* **An interim health count is not a trend** (4.8's term set had fewer C-N excursions mid-run and the largest tears when
  finished). Report health only from finished runs.

### 10.22 When another computer owns the shared records, write this session's to separate files (2026-10-09)

User instruction. **When another computer owns the shared md files, write this session's records to separate
`*_<computer>.md` files; git must show `plan.md`, `findings.md`, `progress.md` and `remote_jobs.md` unchanged. The owner
merges them and deletes the separate files.**

### 10.23 Write a remote file only from a local file that is redirected in and checked non-empty (2026-10-09)

`ssh ... 'cat > $W.new && mv $W.new $W'` without `< file` replaced the live `WATCH_STATUS.md` with nothing. **Every
remote write goes `ssh ... '... cat > $W.new && [ -s $W.new ] && mv $W.new $W' < LOCAL_FILE`, with the remote md5 checked
against the local one afterwards.**

### 10.24 Recompute a comparator from its files before citing it (2026-10-11)

Findings 5.5 cited TM4 `pp_fail` 3.2% and `acc` 0.42 from "the one surviving set"; its own arrays give 0.20 and 0.11, the
3.2% came from another ladder (5.3), and the files are dated 08-27, not 09-04. **Before quoting a stale set as the
comparator, open its files, recompute the number at the cited condition, and record the path, the file dates and the
value read.** **A script that changes directory takes absolute paths only** (a scratch script that `cd`s into the repo
wrote its output into the repo root).

### 10.25 Deleting run data: remove the data, keep what regenerates it (2026-10-11)

NP-1AO6's four run directories were deleted whole; the user then asked to keep the setup scripts, and two analysis
scripts (`prod_counterion/onset.py`, `onset2.py`) had no other copy. **Before deleting a directory, list its non-data
files (`*.py`, `*.sh`, `*.sbatch`, configs, inputs) and remove only the trajectories and per-step outputs; write the
manifest first.**

### 10.26 Show the evidence before asking the user to decide (2026-10-13)

The user proposed switching glpG's TM4 HDX to ff3.0 because of "a known defect" in ff2.1 with the hybrid; the fixed-table
records say otherwise (4.15, 5.5), and I asked which runs to keep without showing that. **When a request rests on a
premise the records contradict, state the contradiction with its numbers before offering options, and put it in the
question itself.** A decision about another session's jobs is confirmed by the user in that session: relay it with its
evidence and record it only once that session confirms.

### 10.27 An inherited check is not exempt from 10.20 (2026-10-15)

See 10.20.

### 10.28 Name runs by what they are, the tag in parentheses (2026-10-18)

A table labelled only "no-BP" read as a missing run. The user's names: `ff21_nobp_dt009` (release `ff_2.1_nobp`) is "the
deterministic side-chain ff2.1", `ff21_ctrl_dt009` (`ff_2.1_ctrl38`) "the probabilistic side-chain ff2.1". **In every
report or table lead with the plain description and put the tag in parentheses**: "deterministic side-chain ff2.1
(`ff21_nobp_dt009`)", "ff3.1 rerun on the fixed engine (`ff31_depth_data_angfix`)". Internal tags (no-BP, ctrl38, a31_,
a32_) alone read as unknown jobs.

### 10.29 A debug goal means we find the bugs ourselves (2026-10-19)

The user's "overall debug of the current code" cited Jumper's Claude audit of master; I proposed asking for his list.
**When the user names someone else's result as the model for a task, the task is to reproduce that kind of result
independently, not to obtain theirs**; an outside list is at most a cross-check afterwards. Upside 3 is a new release so
that no list of fixes and their effects on ff2.1 has to be published: **for Upside 3 work, prove each fix correct
(reproducer, finite differences, nothing else changed); do not plan effect measurements against ff2.1 or master.**
Upside 2 work keeps the master-parity rule.

### 10.30 A long job finishes without an agent, and no agent waits on it (2026-10-20)

User correction: a session used up the allowance and reported nothing. Two long-lived agents (one re-reading ~650k
tokens per turn) spent 113 turns waiting in `until <done>; do sleep; done` loops cut by the 10-minute Bash timeout, each
wake-up re-reading the whole context (11.3M cache-creation tokens in 42 turns for one); the simulations finished, but
the analysis existed only in the agents' plans. Rules:
* A long job is one driver script that runs the whole chain, simulation through analysis, and writes a small results
  file ending in `DONE`. Launch it with `scratchpad/detach.sh LOG CMD...` (own session via `setsid`, `caffeinate -i`
  keeps the Mac awake), so it survives this session, its agents and a usage limit. Keep durable outputs in the repo
  `scratchpad/`, not `/private/tmp`.
* Wait with one `run_in_background` waiter that exits when `DONE` appears. Never wait in agent turns.
* Agents are short-lived: one task, write results to a file, report, end. An agent that starts a long job hands back
  the results path instead of waiting. Do not resume an agent with a large context; start a fresh one from a brief file.
* Run at most two or three heavy agents at once and start the next when one reports (10-21: six reviewers launched
  together, beside two other agents, hit the session limit within minutes).

### 10.31 Cut a delimited block by searching for its end from its start (2026-10-19)

Removing a `# --- START CHANGE --- ... # --- END CHANGE ---` block from `CMakeLists_M1.txt`, a search for the end marker
from the top found an earlier block's, and the slice duplicated half the file (two `upside` targets, which still built).
**When cutting text between markers, search for the end marker from the start marker's position, assert the counts, and
read the diff against HEAD before building.**

### 10.32 Emails to Tobin and John carry the experiment, not the session (2026-10-21)

A page-long draft for the PI and Jumper (library mechanics, a trainer bug, convergence probes) was cut by the user.
**For a results email: what was compared and why (one or two sentences), each test set with its size and what it is, the
numbers, and what we conclude.** Leave out bugs found along the way, gradient checks, job mechanics and file paths
unless asked.

### 10.33 A block SE needs blocks longer than the correlation, and a shared seed is one sample (2026-10-21)

The equipartition test's fixed block count used 60-120-frame blocks, but in the 1rkl hybrid lipid kinetic energy stays
correlated for about 55 frames (25 time units; the diffusion-calibrated friction gives m/gamma 35 time units), so eight
copies scattered 3-4 times beyond their own SEs and read FAIL. A 0.995 protein value seen with seed 2 in two builds was
one draw seen twice: protein carriers are damped in 0.25 time units, so their momenta follow the bath noise, a function
of (seed, atom, step) only (per-frame KE r 0.92-0.96 across builds); independent seeds show no deficit. **Take a
statistical test's SE from the measured statistical inefficiency (pymbar `statistical_inefficiency_multiple`), with
each series many g long. Compare builds or settings on independent seeds, or say that the runs are paired. Keep a
deciding run's trajectories until its conclusion is written.**

### 10.34 The diffusion sampler is built on Boltz, never BioEmu (2026-10-21)

The user: we build on Boltz; Jumper does not accept BioEmu because its data are selected. **No BioEmu weights, samples,
training data or benchmarks enter the diffusion work, as baseline or comparator.** BioEmu2's code may be used (user,
2026-10-21; MIT-licensed and itself adapted from Boltz, so Phase 3 ports its lag conditioning): the objection is to the
data. When a plan document names a base model, check it against the user's stated goal before running it. **One
exception (user, 2026-10-24, "Lift it, side test only"):** a side test of whether BioEmu2's weights give HDX results,
kept as its own record (`diffusion_dynamics.md` section 7) and never a baseline or comparator in the Boltz work.

### 10.35 Read the fix's diff as master's author would before verifying it (2026-10-21)

User, during the D5 fix: keep code edits concise, logical and in master's style. My draft rewrote Jumper's whole
comment in `solve_clamped_1d_spline` where two lines had changed, added a pointer, two ternaries and a second derivative
function, and stripped an unrelated trailing space. **After every fix, read `git diff` before any verification run and
cut it to the smallest change that reads as written in one pass: keep the original author's words where they still
hold, one function per quantity, no incidental whitespace.** Verify the trimmed code, not the draft.

### 10.36 Draw a quantity the project already plots in the project's own format (2026-10-22)

User, on the glpG ff3.0 dG figures: "it should be in same format as glpG-RKRK-79ALA_S115T_POPEPOPG_dG_vs_residue.png.
your format is completely wrong". The project draws residue dG with `example/00.AnalysisScripts/helpers/plot_ref_style.py`,
and the user keeps the 09-04 set in `~/Downloads` (T 0.75 / 0.80 / 0.85, 3000 x 2000). **Before drawing a figure of a
quantity the project already plots, find the user's existing figures of it (`~/Downloads`, the work dirs' `results/`)
and the script that made them, and reproduce that format with that script unchanged.** "There should be only 4
figures": **one figure per system, from all its post-equilibration data pooled**; a comparison the format lacks goes in
text unless the user asks. The SciencePlots rule governs new figures, not a format the user names. The repo's copy
crashed on every input from 10-21 (its legend reused the name `labels`, which the tick formatter reads) and was fixed
10-24 by renaming the legend's list; the cluster's older copy (md5 `fad2166e`) never had the clash.

### 10.37 For runs from different starts, "converged" means they agree with each other (2026-10-23)

User, on ff3.1.1 and ff3.2.1: "by converge I mean ff3.1.1 and 3.2.1 have similar rama map, while have different starting
point for GLY". **When the user asks whether runs that differ only in their start "converged", the test is between the
runs:** the gap at the start against the gap at the end, term by term, read against the gap of the terms that started
equal (the runs' own noise). The fixed-point gate (1.36) is a different question; say which one is being answered
before running either.

### 10.38 Name an HDX residue from the sequence before comparing it with experiment (2026-10-23)

The HDX outputs (`resid_*.txt`, the npz `res`, the cluster `plot_ref_style.py` axis, `readout.py`) label each amide by
its 0-based Upside index, the donor's residue, while glpG's PDB, construct, variants and HXMS peptides count from 1:
label = PDB - 1. I told the user label 144, Val145, was in peptide 140-144. **Before comparing an HDX residue with a
peptide, a crystal contact or a variant, print its residue name from the HDX topology's sequence and check it against
the peptide's sequence** (ALMGY is labels 139-143). The repo's `plot_ref_style.py` labels its axis with `pdb_resnum`;
the cluster copy does not.

### 10.39 Count a fleet's files against the inode quota before submitting it (2026-10-24)

20 Boltz-2 jobs were checked against `/beagle3`'s free space (1.2 T), not its inodes; at 5 files per sample they needed
1,000,000 inodes of a group fileset with ~60,000 free, filled them within an hour (`df -i` 214 of 1,394,800 free), and
every trsosnic job creating a file on `/beagle3` failed. **Before a fleet that writes many small files, multiply files
per unit by units and check it against `df -i` on the fileset as well as `df -h`**; a sampler that writes a file per
sample packs them inside the job, one tar per run (`phase1_seeds.sbatch` through `pack_predictions.sh`). **Judge such a
tar by comparing it with its directory (`tar --compare` and a file count) before removing anything, not by tar's exit
status**: GNU tar 1.30 exits 1 ("file changed as we read it") when a directory's timestamps move during the read, as on
`/beagle3` for a directory finished 3 s earlier, every file intact.
**After a quota incident, sweep every job of every line over the whole window, up to the moment the quota was freed,
and look at the jobs that start in it**: the check at 13:10 found no midway2 failure, but Upside 3's REMD block 1 and
four of its Peng chunks failed 13:12-14:11, found only at 23:50. A job whose log cannot be created fails in seconds with
no log (exit 0:53), so only `sacct -X -S <start>` sees it. A chained job that dies writing its own progress (`run.py`'s
`block_count`) does not queue its successor, and rerunning its stage without that file would redo block 1 over block 1's
data: restore the progress file to what the finished step would have written, then resubmit by the script's own path.

### 10.40 A reference run is valid only if it carries every fix the run it judges carries (2026-10-26)

User correction. I judged Upside 3's glpG TM4 loss against "Upside 2 with ff2.1 held TM4" (0.944 in 12 seeds, 0.93 in
REMD), but both runs predate the 10-17 side-chain placement fix (4.17), which reaches the packing a buried helix depends
on. **Before citing a run as the reference for another, list the hybrid defects fixed between them (3.10, 3.11, 4.15,
4.17 and the audit ledger 3.14) and state which the reference carried; if it carries one the judged run does not, it is
no reference, and the matched control is owed instead.**

### 10.41 Never run Python from a shared directory on a cluster (2026-10-27)

Own mistake. A stdin script run with `cd /tmp` on the midway3 login node imported another user's `/tmp/inspect.py` in
place of the standard library's `inspect` (the current directory is first on `sys.path` for `python -`) and ran under
our account until its fourth line failed on a permission error. **Run cluster Python from our own directory (`$B`, `$E`,
the run directory) and with `python -I`, which drops the current directory and the user site from `sys.path`; never
from `/tmp`, `/scratch` or any directory others can write.**

### 10.42 A trained force field carries physics it has no explicit term for (2026-10-28)

User correction: I wrote that Upside has no electrostatics and that charge-dependent IDP entries (charge patterning,
net-charge expansion) would fail by construction. Upside's core has electrostatics baked into its trained terms; what
it lacks is an explicit electrostatic term, and an ionic-strength input. **Say "no explicit X term", never "no X", for
a trained model; and predict failure by construction only where no term can respond at all (here salt, since nothing
reads an ionic strength). Whether the implicit physics is enough is a measurement, not a deduction from the term
list.**

---

## 11. Reference: the two glpG PDBs disagree about TM4

`glpg_oriented.pdb` (crystal-derived, oriented post hoc) puts residues 135-140 at z = +2 to +9 A, splayed toward the
extracellular surface, so DSSP fails them; the membrane-equilibrated REMD structure has them 6-10 A deeper, helical. The
simulation PDB agrees with the 2IC8 literature boundaries (construct offset -66) to within 1-2 residues at every helix,
while the crystal PDB gives a broken TM2 (83-96 against 82-103) and a truncated TM4 (141-149 against 135-151). **Use the
simulation PDB's DSSP for figure annotation**: TM1 29-48, TM2 82-103, TM3 105-127, TM4 135-151, TM5 161-175, TM6 185-207.

The structure we simulate (checked 2026-10-09, 10.21): `scratchpad/local_popg_79HIS/glpG-RKRK-79HIS.pdb`
(= `pdb_staging/glpG-RKRK-79HIS.pdb`, md5 `e8c3e1a5...`) is the live seed's t = 0 protein (CA RMSD 0.00 A), with the six
helices above and DSSP I at 92-96 in TM2 (not a gap, 4.6). CA depth spans 21-30 A against PO4 planes at +-21 A; 19-26
and 50-57 lie flat in the interfaces and are not TM. Its B-factor column (88-96 in TM2, 39-58 at the termini) reads like
a prediction confidence (unconfirmed).

Construct mapping, verified against UniProt P09391 and the experimental HDX file: construct index + 66 = E. coli GlpG
numbering; the base construct is already the catalytically dead S201T (construct 135), `79HIS`/`79ALA` is WT H145 vs
H145A, `S115T` is S181T, and RKRK is the C-terminal tag at 207-210. No proline is in the donor list (Upside excludes all
six) and there are no chain breaks or resseq gaps.

---

## 11b. The Peng 2022 benchmark: what the SI actually says (read 2026-09-14)

Both PDFs are on this Mac and must be read rather than searched for (ACS is paywalled, PMC serves a CAPTCHA): SI at
`~/OneDrive - The University of Chicago/ct1c00960_si_001.pdf`, main text alongside. The paper is the HDX paper,
"Prediction and Validation of a Protein's Free Energy Surface Using Hydrogen Exchange and (Importantly) Its Denaturant
Dependence", JCTC 2022, 18, 550-561; the folding benchmark is in its SI. Verified against it:
* The simulation-parameter table (p10-11) matches `bench_table.py:TABLE_S2` verbatim (durations and 14-rung ladders);
  Fig S4's terminal-residue exclusions match `RMSD_EXCLUDE`; `FF2_BASELINE` was transcribed from Fig S4's printed
  per-protein lowest Ca-RMSD (largest-cluster centroid in parentheses, DBSCAN on the Ca contact map, 10 A cutoff).
* Fig S4's right-hand panels pool five independent simulations per protein; our arms are one run each at the bottom
  rung, so our distributions are narrower by construction. Both figures are grids of per-protein distribution curves
  with the set average on top (native and unfolded starts, FF1 and FF2 overlaid).

**Do not describe a figure as reproducing a published format without having opened that figure**: the first panels were
built from a one-line paraphrase in `analyse_bench.py` and were wrong in layout. Say "same quantities, my layout" unless
the source figure has been read.

**Per-residue HDX dG are printed in the SI (read 2026-10-22).** Experiment: Table S5 (EHEE_rd2_0005, 15 amides, NMR HX
at 298 K, pD_read 7.1), Table S6 (HEEH_rd4_0097, 20 amides, 278 K, pD_read 4.6, 2.2 to 4.3 kcal/mol at 0 M) and Table S7
(ubiquitin, 24 amides, 273 K, pD_read 7.6), each over urea. The designs are numbered by their design sequences (40 and
43 residues, His tag removed; SI p. 24); Peng's HEEH design model is `/project2/trsosnic/share/HEEH_rd4_0097.pdb`. Upside
FF2: Tables S8 and S10, at 293 and 300 K on Upside's calibrated scale (Figs 5 and 9), 0 M correlations with experiment
0.631 and 0.455 (Table S14). Traps:
* The ubiquitin is **Ub psWT (F45W/H68N)**, in experiment and simulation, not 1UBQ's wild type.
* The paper never says why amides are absent from Tables S5 and S7; for EHEE most have assigned NMR peaks (Table S12,
  His-tag numbering, subtract 21). Do not read an absent amide as fast-exchanging.
* Table S9 (HEEH, simulation) repeats Table S8, an SI copy error; Upside FF2's HEEH agreement survives only as Table
  S14's summary (0 M: Pearson 0.185, Spearman 0.317, RMSE 1.06 kcal/mol).
Transcribed 0 M values: `scratchpad/boltz_hdx/flags/reference/`.

## 11c. Peng's benchmark trajectories on midway2 are FF1, not FF2 (measured 2026-09-14)

`/project2/trsosnic/condiv_data_upload/trajectories/` holds `<prot>_native.xtc`, `<prot>_denovo.xtc` and `<prot>.pdb`
for the paper's 16 proteins plus 7 CASP targets (2.8 GB, 2018-07-20). **It is FF1**: scored with
`ff3_benchmark/scoring/tmscore.py` and the paper's exclusions it gives 0.480 from native and 0.360 de novo, against
published FF1 0.45 / 0.37 and FF2 0.55 / 0.42 (the paper's native figure counts only excursions within the native
basin). That also validates our scoring path externally. Reading the `.xtc` needs mdtraj, installed out-of-tree at
`/beagle3/trsosnic/yinhan/pylibs` (`pip --target`; add to `PYTHONPATH`).

A second FF1 set in Upside's format (search for `.up`, not `.xtc`): `/project2/trsosnic/share/paper_traj_nabil/
{from_native,from_denovo}`, 14 replicas for cspa, gpW, hyp, nug2, top7 (28 GB, 2018-02); its
`rotamer/pair_interaction/interaction_param` is (20, 20, 62) against (20, 20, 54) for ff_2.0, ff_2.1 and ff_3.0.

**There is no FF2-era benchmark data on midway2 or midway3** (a 2.8 TB sweep of `/cds3/trsosnic`, `pengxd/`, the
FF2-era ConDiv trees, and `/project2/trsosnic` and `/project/trsosnic` for `*_native.xtc` and TM/RMSD outputs newer than
2021). An FF2 per-protein overlay needs a higher-resolution figure from the publisher or a fresh ff_2.1 run of the 32
arms (`bench.sbatch` takes `FF=ff_2.1`).

## 11d. lambda repressor in the Peng benchmark (2026-09-18 to 10-02)

lambda is the worst arm of the Peng benchmark under FF2 and every force field tried here. Helices are named from the
reference's own (phi,psi), 0-based file numbering (file index = residue - 6): H0 3-23, H1 27-34, H2 38-46, H3 53-63,
H4 72-78. From the first diagnosis, on the retired ff3.0 (`ff3_benchmark/scoring/`: `diag_lambda_fail.py`,
`gly_reweight.py`, `energy_split.py`, `helix_by_rmsd.py`):
* **The failure is bundle assembly, not secondary structure**: each helix keeps its shape (local CA-RMSD 0.6-2.8 A)
  while the assembly is 8.6 A out, near-parallel native pairs coming out near-perpendicular (H0-H3 28 -> 69 deg, H1-H4
  22 -> 86 deg). proteinB and homeodomain reproduce every crossing angle within 17 deg.
* **The native arm never equilibrates** (cold-rung CA-RMSD 6.4 A in the first block to 10.4 A in the last; de novo
  10.6-11.5 A), so `score_arms_dist.py`'s `BURN = 2000` (one of twelve blocks) mixes decay with equilibrium (10.9).
* **The Ramachandran term barely tracks fold quality** (correlation with CA-RMSD +0.075), and glycine-handedness
  reweighting changes the equilibrated ensemble by -0.024 E_up (10.9). Do not use frame 0's energy as a reference:
  proteinB's deposited structure scores +163 E_up above its ensemble mean while folding perfectly.
* **H2 (file 38-46, `QSGVGALFN`) fails first**: intact (89% alpha_R) only below 5 A global CA-RMSD, while H3 stays
  93-99% alpha_R. The six glycines are mixed in handedness, so no context-free glycine map suits lambda: file 24 L-G-L
  alpha_L (H0 C-cap), 35 M-G-M alpha_L (H1 C-cap), 37 M-G-Q beta/pPII, 40 S-G-V and 42 V-G-A alpha_R inside H2, 47 N-G-I
  alpha_L (H2 C-cap).

### Wild type against G46A/G48A under ff_2.1 (2026-10-01/02)

`ff3_benchmark/pdb/lambda.pdb` is byte-identical (md5 c85500a3) to Peng's
`/project2/trsosnic/pengxd/Data/test-set-15/lambda.pdb`: lambda repressor 6-85, wild type at the positions the
fast-folding variants change (D14, Y22, Q33, G46, G48). H2 above is lambda's helix 3 (44-52), and its glycines 40 and 42
are **G46 and G48**, which the fast-folding variant replaces with Ala to stabilise this helix. The causal test: both
sequences, native and de novo, ff_2.1, Peng's protocol, coldest replica T 0.780, one trajectory per arm
(`checks/lambda_ff21/lambda_check.py`).
* **Wild type: helix 3 unravels first and alone** (aR 0.97 -> 0.57 over the first chunk, the other helices 0.89-1.00),
  and the global RMSD rises with it: ff_2.1 already does what the old ff3.0 did. WT native is unfolded from ~1 M tu
  (10.2-10.4 A); WT de novo never goes below 5 A. **The two WT arms meet on helix 3 from opposite sides** (last block aR
  0.61 native, 0.57 de novo): about 60% formed at equilibrium.
* **G46A/G48A holds helix 3** (aR 0.97-0.99 in every native block), and its native bundle drifts to a plateau near 7 A
  (block RMSD 5.47, 5.72, 6.95, 6.88, 6.93 A; 0.31-0.44 of frames below 5 A over the last 1.5 M tu) instead of
  unfolding; H0-H2 native contacts stay lost (0.03-0.15). De novo it does not fold (at most 1% below 6 A, helix 3 aR
  0.80-0.88). So the two helical glycines are the main cause, and a smaller packing drift remains. Experimentally WT
  lambda 6-85 folds with helix 3 its least stable helix, so helical glycines are too weak in Upside itself (the
  selection panel's ff2.1 deficit, -0.081 against all-atom, 1.17).
* **What breaks once helix 3 holds** (`lambda_packing.py`, `lambda_terms.py`): helix 2 packs against helix 4 30-40 deg
  off native from the first block in both sequences (H1-H3 native 120 deg; WT 80 -> 50, G46A/G48A 90 -> 81), a
  mis-specified interface rather than a slow decay. ff2.1's energy prefers the misoriented packing (G46A/G48A: total
  -169.0 at 70-90 deg against -163.0 at >= 105 deg; `rotamer` -164.4 against -157.0; H-bond and coverage favour native).
  At ff2.1 the side-chain gradient is noise (6% of coefficients significant), so ConDiv sees no consistent correction.

### The helix 2 / helix 4 misorientation under ff2.1 has no single side-chain cause (2026-10-02)

Decomposed on both native arms (misoriented = H1-H3 70-90 deg, native-like >= 105 deg; errors from 8 time blocks;
`checks/lambda_ff21/{pos_decomp,pos_report,panel_terms,panel_report}.py`): the rotamer preference is weak (-7.4 +- 3.5
E_up in G46A/G48A, -2.3 +- 3.0 in WT), only `hbond_coverage_hydrophobe` favours misorientation in both (-2.7 +- 1.8,
-5.1 +- 2.2), and no residue-type pair carries it (each < 1 E_up, signs flip between arms); in G46A/G48A the misoriented
state gains a Q33-F51 dock (-2.3 +- 0.3). On the ff2.1 panel F_rotamer and the side-chain/backbone term favour the
native over compact non-native frames in 31 of 32 domains, though those frames are mostly frayed natives. No table
entry has a physical target to correct: lambda's packing is a known limitation, re-tested after ff3.0 (crossing angle;
Q33-F51 dock; helix 1-helix 2 retention; L50-L64), with more native-like frames than the 196 and 102 here.

## 11e. Large proteins de novo: compact non-native globules; N-to-C growth does not help (2026-10-19)

The full record, its tables and the Upside 3 plan are in `upside3.md` (sections 2-3; plan.md Phase 22; data in
`scratchpad/large_protein_20261019/`).
* The mdCATH native-start panel's drop at 151-250 residues is the folded-fraction cutoff, which tightens with N;
  Upside's Q offset from all-atom and its Rg do not change with size.
* T4 lysozyme, RNase HI and flavodoxin II (in-vitro refolders, 155-179 residues) and the homeodomain control under
  ff_2.1: native starts hold. At T 0.80 every coil, collapsed in bulk or grown N to C, ends in a compact, helix-poor
  globule without native tertiary contacts, 150-180 E_up above the quenched native on average, mostly in missing
  H-bonds. Its side chains recover 87-93% of the native's side-chain energy; Rama and backbone desolvation favour it.
  At 0.88 RNase HI's native holds and its coils still do not fold; at 0.95 all melt. N-to-C growth equals bulk.
* Reading: a search failure into states ConDiv never contrasts with the native. Not measured: free energies, and a
  size contrast under one protocol (the control does not fold at constant T either).
* The coverage-typing defect found on the way: 3.13.

## 11f. Upside 3 planning: units, GPU and the first audit pass (2026-10-19)

The full record is in `upside3.md` sections 2, 3 and 6 (reports in `scratchpad/upside3_planning_20261019/`).
* **E_up's kJ/mol value is a temperature calibration.** 350.588 K = 298 K / 0.85, and Peng 2022 SI eq. S14 matches
  T_up 0.85 to 298 K for FF2; FF1's match of 1.0 gives 2.479 kJ/mol. Converting parameter files to kJ/mol fixes one
  force field's fit into them, and the MARTINI factor is this same calibration.
* The Rama, rotamer and BP free-energy terms are evaluated at a fixed kT_ref = 1, at any simulation temperature.
* Master's analysis scripts disagree on the kelvin scale: 0.593 kcal/mol per T_up implies 298 K, and 41.4 pN is kT at
  300 K, where the 350.6 K convention gives 48.4 pN.
* **The cluster and Mac builds treat NaN differently.** Master has `-ffast-math` on both. This branch added
  `-fno-finite-math-only` only to `CMakeLists_Other.txt`. Under the Mac flags `isfinite(NaN)` is true (probe; not yet
  shown in an engine run).
* **Profile, ff_2.1 on one M3 Max core.** BP is 20-24% of the time, side-chain pair edges 25-27%, burial 19-22% and
  H-bond 13-19%, from 76 to 578 residues. In the hybrid, the dry-MARTINI pair loop is 86%.

## 12. Claims that turned out to be wrong

One or two lines each: what was believed, what is true, and the lesson.

### 12a. The glycine campaign, 2026-09-18/19

1. **"Glycine handedness is zero."** An unconverged flat surface (`awh1-dimN-diffusion` ~15000x too small); 6.7.
2. **The beta-branching rule** (VAL and THR positive because branched): both crossed zero by 30 ns. Do not run ILE as a
   decisive test.
3. **Pentapeptide controls**: GGGGG must read 0 and reads -0.224. Not usable.
4. **"LA is settled" at 0.028**: a four-snapshot window; 0.106 over 18 ns (6.7).
5. **"Everything drifted toward zero" between 51 and 57.6 ns**: a basin-definition artifact.
6. **ff3.0C's premise**, that the library's neighbour specificity is real: its neighbour ordering correlates with the
   measurement at r = -0.565 (ff2.1 -0.596), both anti-correlated.
7. **"The candidate glycine map stabilises lambda's near-native basin by 0.55 kT"**: reweighting of a non-stationary
   arm; -0.024 E_up on the equilibrated third (10.9).
8. **"Each replica's achiral control converges to its own nonzero value"**: sampling noise, uncorrelated between
   replicas and decaying as 1/sqrt(t); the lasting consequence is ~20% on the basin dG (6.7).
9. **"The trainer's `hb` has never been trained"**: read off the modern port; the Theano original trains it (9e).
10. **The sheet gradient is systematically nonzero** (t = +3.91 at n=3): t = -0.01 at n=6 (9e).
11. **"The library's symmetric part must be kept" (rms 1.752 E_up, 25x the handedness correction)**: those numbers were a
    per-pair antisymmetric statistic; the symmetric parts correlate at r = +0.867 with aR basins within 0.207 nats.
    **Name the quantity before quoting a ratio on it.**
12. **The units conversion was inverted**: a library map holds `-lnP` normalised to `sum(exp(-E)) = 1`, so a PMF entering
    it is divided by `kT`, not by 2.914952774272. **Check the file's own contract before choosing a unit conversion.**
13. **`S_library + A_measured` presented as free of ff3.0**: `S_library` is ff3.0's symmetrised map (caught by the user).
    **When a construction decomposes something, check whether one piece is a thing already rejected.**

### 12b. The hybrid, HDX and cluster work, 2026-07 to 09

* **findings 87 (withdrawn by findings 88, cited elsewhere as findings-88)**: glpG blow-ups as a timestep failure at
  protein-lipid contacts. The `omega*dt` contacts were O sites whose force the engine discarded, leaving 3590 E_up/A of
  net one-sided force on the environment. Check that the stiffest interaction is connected to the dynamics.
* **findings 88's prediction** that fixing the O force would bring `avg_kinetic_energy/1.5kT` to 1.000: refuted; the
  excess was finite-dt error (4.2).
* **`--potential-deriv-agreement` as a correctness gate**: round-off dominated at this size (pre- and post-fix 0.41960
  and 0.41959). A developer probe only.
* **findings 90 (corrected by 92)**: the force-free core "~500 kT out of reach", computed for an inertial particle while
  ION and LIPID are overdamped Brownian. Check which integrator governs the particles before a stability margin.
* **findings 92/93 as the trigger**: the zero-force core is where the cascade ends; the trigger is two environment beads
  at 1.83 A, inside the valid table domain.
* **findings 95 (retracted by 118)**: K525 and K541 as reduced-labelling sites, an artifact of the uncorrected table and
  CB placement (9).
* **findings 96 (retracted by 104)**: the `mean_pf >= 0.99999` clip called a bug; it is master's convention encoding the
  estimator's statistical limit. Diff an inherited threshold against the reference and ask what limit it encodes first.
* **findings 104's follow-on (corrected by 113/114)**: the spikes were the missing lipid-shielding term, not the
  unclipped rendering; what stands from 104/109 is the resolution limit.
* **findings 106 item 4 (corrected by 113, wrong by 7x)**: the missing membrane term as "6 amides"; with tail contact
  rather than midplane distance deciding where water is, it is 44 extra `+inf` donors and explains the figure.
* **findings 107 (narrowed by 108)**: "3 amides of 148" holds for the protection state only (PS is `H-bond OR burial`,
  and burial does most of the work in a membrane protein).
* **findings 121's first inference** (cluster NaN as an output artifact): the replicas blew up and were rolled back.
* **findings 126's C-terminal O** as an unconstrained particle: the engine never writes that slot (3.8 Bug 3). Check
  whether the engine writes a particle before calling it unconstrained.
* **The implicit-vs-hybrid comparison on the protein-only protection state** (a second Update 124) inverted the sign of
  the result (5.3); its follow-on claims are retracted with it.
* **The zero-variance cluster HDX read as "needs more blocks"**: the protein was rigid (findings 116).
* **The stage-7 freeze accepted on kinetic energy and secondary structure**, which a high-friction g-JF satisfies while
  barely moving (6.3).
* **"195 GB left, so the pre-ff3 ladder must go" (2026-09-08)**: that was the `/project2` quota; `/project` had 1514 GB
  free. Match the quota section to the mount the data is on and cross-check with `df` on the data path (it nearly deleted
  89 GB of baselines).
* **The BB-env PMF**, built to fix a "kick" that was a setup artifact (a timestep inherited from the abandoned CGL plus
  under-resolved lipids driving a displacement cap), then caused the drift it was meant to prevent. Rule out setup
  artifacts, timestep and sub-step resolution first, before building a corrective term.
* **Rewriting `plot_ref_style.py` to draw censored amides as bounds (2026-09-12, reverted)**: it overrode a rendering
  argued in the file's docstring. A file's documented rationale outranks my judgement of how the output should look
  (change it only if the user asks or the rationale is false, and say which); never collapse right-censored values onto
  one ceiling or break a curve read as continuous; when output "looks like it lacks data", fix what is plotted first.

### 12c. The ff3.0 TM4 repair, 2026-09-30 to 10-07

* **"ff_3.0 (09-30) flips TM4's helical glycines in glpG."** That run had ff_3.0's map, H-bond and pair tables beside
  the FF1-form ff_3.0's coverage tables (3.11); ff_3.0 as released was never measured. This observation started rounds
  3 and 4.
* **Every local TM4 count, 10-04 to 10-07**, ran on the same foreign coverage tables (3.11). **When a test patches a
  force field into a prebuilt input, list every node the force field feeds, and check that each one was replaced.**
* **"ff2.1's H-bond energies came from the H-bond node rewrite, not from a trainer"** (1.17): Peng's SI says ConDiv
  returned them (1.23). Read the SI before saying where a parameter came from.
