# Findings

## NP-1AO6 ion build: box was oversized, now counterions only (2026-08-04)

- FINAL: **neutralizing counterions only, no bulk salt** (user decision 2026-08-04). 218 K+, 0 Cl-,
  cancelling MPA (-203) + protein (-15). System is 4198 particles (was 7868). `salt_molar`,
  `estimate_salt_pairs`, and the free-salt assertion are removed from the NP path entirely.
- Salt molarity was never the defect. Every build measured 0.148--0.150 M, matching both references.
  The defect was **box volume**: `box_len` was derived from the rotated protein's reach *about the gold
  COM*, and because the frozen NP was pinned at the box centre while the protein adsorbs to one side,
  that doubles the protein's lever arm. Boxes came out 232--284 A for a complex only 122--148 A across.
  Correct ion density x inflated volume = too many ions. Now: **complex-centred, fixed 200 A** for all
  six faces (tightest is 180-0-0 at 198.1 needed). `write_up` had to change too -- it previously pinned
  the gold at the origin and wrapped to +/-box/2, which would wrap-split the protein in a 200 A box.
- The two references disagree on counterion bookkeeping, which matters at small box size:
  `hernandez.onedrive/scripts/lib.py:482` (the paper, LAMMPS) spends counterions **from** the salt
  budget (`N_ion_pairs -= 1` each); CHARMM-GUI dry-MARTINI (`step4.3_ion.out:604`,
  `calc npos = @npos + int(@conc*6.021e-4*@volumn)`) adds them **on top**. In a 400^3 box the choice is
  invisible (0.144 vs 0.150 M); in a 200^3 box it is 0.104 vs 0.150 M. Moot now that salt is dropped,
  but record it before reintroducing salt.
- The paper's own simulation is **wet and fully salted**: `0-0-0/input.data` has 519,306 W beads, 5781
  K+, 5548 Cl- (net charge -233; note its MARTINI 1AO6 carries -30, not the -15 our residue-based
  `infer_protein_charge_from_residues` returns). Ions are 2.1% of its particles vs ~49% of ours at the
  same molarity -- that ratio, not the molarity, is what made the ported system look wrong.
- Why it survived repeated rebuilds: nothing asserted the *built* composition, so every rebuild
  re-derived the count correctly from an unexamined premise and no check ever fired. `build_system` now
  asserts exact charge neutrality and zero salt pairs from the placed ion counts, and rejects a box too
  small for the complex. All guards negative-tested.

## NP-1AO6 dt=0.005 re-run blew up on all 6 faces (2026-08-03)

- The dt 0.009 → 0.005 reduction plus the graded soft-core ramp did **not** fix the backbone tearing. All
  six orientations are unusable; the failure is the same surface-driven one, only slower.
- `isfinite` is not a health check. `np.run.1.up` is finite for all 52,511 frames yet is physically
  destroyed: protein Rg drifts 26.5 → 76.7 A, max peptide C--N bond reaches 56 A with ~300 bonds over
  2 A, and the potential goes from -12,000 to +10^5. Screen on Rg, peptide C--N, and the sign of the
  potential instead.
- Output-group order in a restarted `.up` is **ascending**, not descending. `run_np_prod.py` renames
  `output` into the lowest free `output_previous_N`, so `output_previous_0` is oldest and `output` is
  newest. `prod/onset.py` scans `previous_1..17` and therefore skips `output_previous_0`.
- Every face had already died inside block 1 (job 52852032, which exited COMPLETED). Job 52881141 then
  ran 25 h on destroyed systems. The adaptive chunk size collapsing from 394.9 to 76.0 time units is the
  cheap tell -- blown-up coordinates wreck the neighbour lists and steps/second craters. Gate the
  self-resubmit on a physical health check so a dead run stops instead of chaining.
- NOT an unthermostatted-NVE drift. `main.cpp:795` computes
  `thermostat_interval = max(1., round(arg / (inner_step*dt)))`, so `--thermostat-interval -1` clamps to
  **1 -- the thermostat fires every step**. Production and pre-production use identical thermostat
  settings, and the stability run reports `avg_kinetic_energy/1.5kT = 0.970`. The energy injection is
  therefore not escaping for lack of a thermostat.
- Pre-production passes only because it is far too short. `run_stability` does 2000 steps (t=10) while
  the tearing starts at t~79--105, i.e. 15,800--21,000 steps. Any local gate has to run past t~150 to
  have a chance of reproducing the failure.
- **Pick the face by its historical onset, not arbitrarily.** Cumulative time to first tearing differs
  by 16x across the six: 90-0-0 t~115, 0-270-0 t~269, 270-0-0 t~948, 0-90-0 t~1059, 180-0-0 t~1134,
  0-0-0 t~1873. A t=400 run on 0-0-0 passing proves almost nothing (it needed ~1873); the same run on
  **90-0-0** is 3.5x past its onset and is the cheap decisive test. Cost this mistake once already.

## glpG 79HIS REMD HDX figure audit (2026-07-28)

- The downloaded July 23 ΔG plot belongs to the first REMD/HDX generation. Its block-2 continuation
  (`52519397`) starts with NaN potential energies in all 48 slots and remains NaN. The HDX extractor and
  plotting path did not reject non-finite source frames, and the MBAR attempt failed to converge.
- The downloaded 478-frame T=0.70 VTF is structurally finite. Exact stock-HDX re-scoring gives zero observed
  openings for 16 of construct residues 30--50; finite examples are residue 47 at `p=0.9393`,
  residue 48 at `p=0.9665`, residue 49 at `p=0.7929`, and residue 50 at `p=0.1862`.
- ΔG is `RT log(p/(1-p))`, not helix occupancy. At `T_up=0.70`, even 2, 3, and 5 kcal/mol mean roughly
  98.4%, 99.79%, and 99.9965% protection.
- Donor IDs stored in `.resid` are zero-based but the VTF/PDB is one-based. Plot x=29 is construct residue
  30. The construct has a 19-residue N-terminal tag, so native GlpG numbering has another offset.
- Exact `p=1` from a finite trajectory is censored, not a measured 1000-kcal/mol value. A defensible plot must
  separate censored markers from finite ΔG points and must not connect them with off-scale line segments.
- Stock PS is backbone H-bond OR acidic side-chain contact OR protein burial. It is not DSSP helicity, and the
  dry-DDM workflow has no calibrated water-accessibility correction. Do not enable a global membrane flag as
  a substitute.
- Temperature MBAR for the coupled hybrid requires the full sampled Hamiltonian. The current protein-only
  re-scoring in `martini_hdx_project.py` is unsuitable for hybrid MBAR; direct per-temperature estimates at
  actually sampled ladder temperatures avoid this extra weighting problem.

## Root Cause of Protein Tilt-to-Vertical (2026-07-16)

The tilt bug is caused by DOPC lipids adjacent to the protein collapsing their
tail-compression coordinate q toward the compressed endpoint, which removes the
orientational anchoring of the protein. This happens because the single-CGL
tail-compression physics is derived from an ISOLATED single DOPC molecule in
vacuum and omits the lateral packing cost of shortening a tail among neighbors.

Proven by three engine/trajectory probes (scripts in scratchpad):
- `tilt_scan2.py`: with the membrane frozen at its native (35 deg) state, the
  conservative SC-CGL+BB-CGL coupling energy is minimized near the native tilt
  (~40 deg). The force field does NOT intrinsically prefer vertical; the well is
  asymmetric (soft toward vertical, stiff toward lying flat).
- `ratchet_test.py`: the coupling-energy minimum tracks the protein's current
  tilt (frame 0 min ~38 deg, frame 500 min ~18 deg, frame 1000 min ~7 deg);
  returning to 35 deg from a relaxed low-tilt membrane costs +1000s E_up. The
  membrane relaxation ratchets the protein down.
- `qwell.py`: restoring energy (native -> vertical) vs uniform lipid q:
  q=0 -> +292 E_up, q=0.25 -> +182, q=0.5 -> +112, q=1.0 -> +36, q=1.2 -> +37.
  Extended (long) tails anchor the tilt; compressed (short) tails let the
  protein go vertical for free.
- `qvsdist.py`: measured q by distance from protein over the trajectory:
  near (<15A) q~0.93-0.99, mid (15-25A) q~0.48, bulk (>30A) q~0.40. Lipids
  touching the protein collapse to q~0.95 by frame 100.

## Code Locations (relevant)

- Runtime q mixing (matches design doc, compression-only, dw/dq=0 outside
  [0,1.209], no extension force):
  `src/martini_cg_lipid.cpp` `compute_single_compaction_mix` (1739-1822);
  SC overlay 4448-4611, target overlay 4857-5069; q dynamics/self-PMF
  `CGLipidCompactionSelfPotential` (2072-2096), `integrate_dynamic_compaction`
  (5170-5187, hard walls at self_coord_min/max = -0.362 / 1.209).
- Builder (isolated-single-lipid, no packing):
  SC relief `_sc_conditioned_tail_relaxation_effective_energy`
  (`py/martini_build_tables.py` 2178-2349); target relief
  `_target_conditioned_tail_relaxation_effective_energy` (2456-2626); self-PMF
  `_fit_compaction_self_pmf` (4961-5059); deformation is compression-only
  (`_apply_tail_axial_compaction` skips delta<=0 at 1734, radial floored 0 at
  1782); pair (two-lipid, packing-aware) reference
  `_pair_conditioned_tail_relaxation_effective_energy` (1794-2040).

## Engine probe note

`upside_engine.Upside.energy()` on a cold-loaded hybrid config returns ~5.2e11
because the `martini_potential` node (generic backbone/ion particle term)
mis-evaluates without the periodic box; irrelevant to tilt. The CGL coupling
terms evaluate correctly and are read per-node via `get_output`:
`cg_lipid_target`, `cg_lipid_rotamer_sc` (sum over placement beads),
`cg_lipid_pair`, `cg_lipid_compaction_self`.

## Update 2 (2026-07-16, after PoC runs): compaction is NOT the physical fix

Direct 1RKL continuation runs (25k/20k steps, same seed) refute the
compaction-stiffening idea:
- Stiffening the self-PMF toward EXTENDED (near-protein q pinned ~ -0.1) made the
  drift WORSE (35 -> 2.5 deg). Allowing compression (baseline) drifts less
  (35 -> 18.8 deg).
- Frozen-q sweep (all lipids at fixed q): tilt final = 5.4 deg (q=0), 8.6 (q=0.5),
  53.8 (q=1.0), 34.4 (q=1.2). Tilt is governed by effective protein-facing tail
  length: long tails -> vertical, short tails -> tilted. Native ~35 deg only near
  q~1.1-1.2.
- Free-q self-PMF well centered at q0=1.05 (near-protein -> 1.12) preserves tilt
  (35 -> 34.9 deg); q0=0.85 partial (35 -> 28 deg).

BUT the tilt-preserving q~1.0-1.2 is UNPHYSICAL:
- Measured compaction q of the reference conformers: isolated-DOPC MC samples
  mean 0.43; the full-resolution near-protein bilayer conformers
  (`interface_ref_bead_positions_nm`, 12 confs) mean 0.475 — SAME as isolated.
  So the physically correct near-protein tail state is q~0.45
  (head-to-tail-end extension ~19.5 A, sensible for a DOPC leaflet), NOT ~1.0
  (~12 A, too short).
- Source packed-neighbor tail-relaxation estimate also gives only modest
  compaction (q~0.2, crude/rigid).

Conclusion: at the physically correct lipid geometry (q~0.45) the model's
SC-CGL/BB-CGL coupling places the protein's fully-relaxed equilibrium tilt near
VERTICAL, not native. The compaction coordinate can force the tilt but only at
unphysical tail lengths (= parameter twisting, forbidden). So the real defect is
in the SC-CGL/BB-CGL BASE interaction content at the physical lipid state, not in
the tail-compaction q.

Static vs dynamic: with the membrane frozen in its native-adapted state the
coupling correctly prefers ~40 deg (tilt_scan2). The drift is the fluid membrane
relaxing (positional ratchet) to make the current tilt optimal, and the tilt
well is soft toward vertical / stiff toward lying flat, so it creeps vertical.

Strongest remaining lead: the `cg_lipid_target` angular (orientation) coupling
was manually reversed on 2026-07-11
(`active_h5_update=reversed_cg_lipid_target_angular_axis_20260711`,
`target_angular_sign_source=validated_against_full_lipid_1rkl_tilt`) — a prior
tilt-chasing manipulation of the base table's orientation dependence. This is the
next thing to investigate for a genuine (non-twisting) fix.

## Update 3 (2026-07-16): angular-axis flip ruled out; refined mechanism

- Flipping the ENTIRE cg_lipid_target angular axis (interaction_param + all 3
  delta layers, the mutually-consistent set) made 1RKL WORSE: 35 -> 8.2 deg
  (vs baseline 35 -> 18.8). So the current (2026-07-11 "reversed") orientation is
  the better one; un-reversing is not the fix. The target angular coupling is NOT
  the bug.
- Data check: in the shipped target table the base repulsive lobe and the
  delta-relief lobe sit on the SAME angular side (mutually consistent). Tensor
  layout is [radial][angular] (radial outer, angular inner), 91 x 321.

Refined mechanism (footprint vs mismatch balance):
- The static coupling at the native-adapted membrane prefers ~40 deg (restoring),
  but the fully-relaxed membrane makes any tilt optimal (positional ratchet:
  lipids translate to fill space). A vertical protein has the smallest in-plane
  footprint, so lipid-lipid packing favors it. That footprint/packing cost of
  tilting (from CGL-CGL) is real, and in reality it is balanced by the protein's
  hydrophobic-matching + interfacial-anchoring drive to tilt.
- Geometry: protein spans ~31 A along the normal (tilted) vs a ~26 A hydrophobic
  core; going vertical raises the along-normal span to ~40 A (MORE hydrophobic
  exposure), which hydrophobic matching should penalize. The protein drifts
  vertical anyway -> the model's hydrophobic-mismatch / interfacial-anchoring
  drive (SC-CGL/BB-CGL) is too weak relative to the CGL-CGL footprint cost.
- Suspected universal cause: the tempered projection at tau=25 T_up over-smooths
  the SC-CGL/target orientational specificity (the sharp tail-attractive vs
  head-repulsive distinction that anchors tilt). tau=25 was chosen for CGL-CGL to
  prevent aggregation; it may be inappropriate for the protein-facing tables.

Ruled out as physical fixes: compaction-q resting shift (needs unphysical q~1.1
= per-protein tuning), and the target angular-axis flip (worsens tilt).

Note: cg_lipid_potentials.tex claims 1RKL stays ~35.9 deg (validation), but the
current CG hybrid drifts to ~17 deg -> possible regression from a prior working
state (many dopc.h5 .bak files exist) OR aspirational/full-resolution claim.

## Update 4 (2026-07-16): packing-anchor hypothesis REFUTED; hydrophobic matching

Per user guidance (anchor = CGL-CGL lipid packing / footprint, "too fluid"),
attributed the drift energy and tested membrane density:
- Energy attribution along the drift: CGL-CGL (pair) energy drops as the protein
  goes vertical (rewards smaller footprint); SC-CGL is noise; BB-CGL mildly
  favors vertical. So lipid packing does reward vertical.
- Membrane density: the DOPC bilayer is a UNIFORM fluid (6x6 occupancy, no empty
  cells) at ~88 A^2/lipid at ~zero lateral tension (box frozen 110.45, pressure
  ~0), i.e. ~30% MORE EXPANDED than experimental DOPC (67 A^2). (Earlier "49
  A^2 over-packed" was a bad close-packed formula applied to a fluid; the
  density-based 88 is correct.)
- DECISIVE TEST - compress to experimental 67 A^2 (fixed box, rigid protein):
  made 1RKL WORSE, 35 -> 1.9 deg (vs baseline 35 -> 18.8). So DENSER/correct
  packing drives MORE vertical, not less. The packing-anchor hypothesis is
  refuted; the current too-loose membrane actually partially MASKS the bug.

Unifying mechanism across ALL tests: the protein goes vertical when the effective
hydrophobic membrane is THICK, and holds tilt only when it is THIN. Thick sources
of vertical: dense packing (67 A^2), long/ordered tails (low q), extended tails.
Thin sources that preserve tilt: short q (1.0-1.2), loose packing (88 A^2) - both
partly unphysical. This is hydrophobic matching.

OPM annotation (pdb/1rkl.pdb): "1/2 of bilayer thickness: 12.4" -> 1RKL was
fitted to a ~24.8 A hydrophobic membrane, THINNER than DOPC (~27-30 A). Its ~35
deg tilt = matching a long helix to that thin membrane. In thicker DOPC, partial
straightening toward vertical is physically EXPECTED; full 35 deg retention may
not be the correct target for DOPC.

All physical levers now ruled out as clean fixes: compaction-q (unphysical),
target angular flip (worse), membrane density correction (worse). The remaining
question is whether the model's DOPC effective hydrophobic thickness (as felt by
the protein via SC-CGL/BB-CGL) is too large vs real DOPC, or whether full tilt
retention is simply not expected for 1RKL in DOPC.

## Update 5 (2026-07-16): full-res confirms bug is CG projection; 40ps route analyzed

Full-res (14-bead) 1RKL rerun, same box (110.8, area ~89), same 40ps clock, 50k
free production:
- TILT: holds ~44 -> 53 deg (mean last20% 56 deg), NOT vertical. CG collapsed to
  ~17 deg. => the tilt bug is purely the CG single-particle projection
  (SC-CGL/BB-CGL), NOT thickness/density (confirmed; user was right).

Alternate route analysis (run full-res DOPC at Upside 40ps + MARTINI x4 factor):
- Numerics: DOPC bond period 1.07 ps; Verlet stability limit 0.34 ps; 40ps is
  118x over the limit. Yet it does NOT blow up.
- Reason: at 40ps the strong stochastic thermostat dominates -> overdamped
  Boltzmann sampling. STRUCTURE is correct (all bonds mean=r0, std=thermal
  0.45-0.51 A; tilt held), but DYNAMICS are dead.
- Diffusion: D=0.0109 A^2/tu; x4 -> 0.044 um^2/s vs experimental 11.5
  (factor needed ~1000, not 4; even x56 only 0.6). So the MARTINI x4 mapping
  cannot recover physical diffusion at 40ps.
- VERDICT: full-res-at-40ps is a valid STRUCTURAL reference (holds tilt) but NOT
  a viable transport model. This is exactly why the CG model exists: removing the
  bonds + the CG mass-scale/GLE calibration is what makes 40ps give correct DOPC
  transport. User's intuition confirmed (with the nuance: it's overdamped, not
  a blow-up).

Path forward: the fix stays in the CG SC-CGL/BB-CGL projection. Full-res
trajectory (outputs/martini_1rkl_hybrid_full) is now the ground-truth reference
for comparing near-protein lipid tail accommodation vs the CG.

## Update 6 (2026-07-16): ROOT CAUSE — CG lipid orientations are frozen/over-ordered

Comparing full-res (ground truth, holds tilt) vs CG (drifts) near-protein lipid
geometry revealed the defect is orientational order:
- Full-res DOPC: <|n_z|>=0.815, S2=0.54, director tilt ~35 deg (fluid, disordered).
- CG CGL:        <|n_z|>=0.995, S2=0.985, tilt ~6 deg (near-crystalline).

The CG orientations are FROZEN: over the whole 200 tu run each lipid orientation
changes only ~6 deg from its initial placement (0.26 deg per 2 ns; autocorr 0.992;
decorrelation time ~125x the run). So the CGL directors are stuck at the ORDERED
initial placement and never sample the fluid disorder real DOPC has.

Cause (leading hypothesis): the design calibrates ONLY translation (mass-scale
0.05 + two-mode GLE) so diffusion is right at 40ps; ROTATION is left at full
inertia with a simple rotational thermostat -> orientational relaxation is far too
slow -> frozen. (Analogous to the full-res bonds being frozen at 40ps.) A frozen
picket-fence membrane cannot tilt/splay to accommodate a tilted protein, so the
protein aligns to the rigid normal-aligned lipids -> vertical drift.

This is a CGL-CGL / lipid-dynamics defect (matches user: CGL-CGL mis-calibrated,
"too fluid" translationally but frozen rotationally). The fix must make the CGL
orientations sample the correct DOPC orientational disorder (S2~0.54) at the 40ps
clock — physically (calibrate rotation like translation) and universally.

Tests in progress: (a) initialize CG membrane disordered (|n_z|~0.82) and check
if protein holds tilt; (b) map orientation dynamics to find how to unfreeze it.

## Update 7 (2026-07-16): CULMINATING DIAGNOSIS — rotation frozen + CGL-CGL can't hold fluid bilayer

Root-cause chain (all evidence-backed):
1. CGL ROTATION is uncalibrated. Translation is accelerated (mass x0.05 + 2-mode
   GLE, giving correct diffusion after the x56=14x4 factor). Rotation is left at
   FULL inertia (1883, unscaled) + an over-damped thermostat (tau=0.008,
   friction 125 vs natural rate 0.03). Rotational diffusion D_rot=T*tau/I~4e-6 ->
   reorientation ~1e5 T_up -> orientations move ~4 deg over a whole run = FROZEN
   at the ordered prep placement (S2=0.985) vs fluid DOPC S2~0.54 (full-res).
2. Frozen picket-fence membrane forces the embedded protein to align vertical.
3. Unfreezing rotation DOES let the coupling relax to the correct fluid disorder
   (rotcal inertia x0.05,tau5 -> S2 0.99->0.45, potential dropped) — so the CG
   coupling's EQUILIBRIUM order is correct (~0.5). BUT:
4. At the correct disorder the CG bilayer LOSES Z-COHESION: lipid z-spread
   doubles (8.3A full-res -> 16.2A CG-disordered) and the protein still goes
   vertical. The single-particle CGL coupling lacks the tail-interdigitation
   cohesion that holds a real disordered bilayer tight. So the CG bilayer is only
   tight when (incorrectly) over-ordered/frozen.

Consequence: the fix is a COUPLED recalibration, not a single knob:
 (a) Calibrate CGL rotation like translation (inertia scaled x0.05 for
     consistency with the mass scale + a proper rotational thermostat/GLE out of
     the over-damped regime) so orientations sample the correct fluid disorder at
     40ps WITHOUT breaking the x56 diffusion.
 (b) Recalibrate the CGL-CGL coupling so a FLUID bilayer (S2~0.5) is stable and
     TIGHT (z-spread ~8A, leaflet ~15A) matching full-res — i.e. add the
     z-cohesion at disorder that the single-particle currently lacks.
 (c) Verify protein holds tilt universally (1RKL + 1AFO), diffusion correct
     (x56~11.5), structure stable.

Ruled-out non-fixes (all tested): compaction-q, target angular flip, membrane
density/thickness, static disordered init (fell apart), gentle unfreeze at full
inertia (too slow to disorder). Full-res-at-40ps route: structurally valid but
diffusion ~1000x too slow (not a transport model).

## Update 8 (2026-07-16): part 1 implemented; part 2 quantified + target characterized

Part 1 (rotation calibration) IMPLEMENTED (physical, universal):
- py/martini_prepare_system.py apply_cgl_transport_metadata: rotational inertia
  now scaled by the same mass_scale (0.05) as the translational mass (it is
  I=sum m_i r_i^2 over the same beads; idempotent from stored base).
- rotational thermostat 0.008 -> 5.0 (system thermostat timescale):
  py/martini_prepare_system_lib.py DEFAULT_CGL_ROTATIONAL_THERMOSTAT_TIMESCALE
  and example/16.MARTINI/run_sim_hybrid.sh.
- Validated (matching rotcal run): membrane relaxes 0.99 -> S2 0.49 (correct
  fluid disorder). Rotation no longer frozen.

Part 2 (CGL-CGL disorder-robust cohesion) QUANTIFIED:
- At correct disorder, CG intra-leaflet z-spread=11.4 A vs full-res 2.9 A
  (~4x looser; z-confinement stiffness ~15x too weak). Diffusion side-effect:
  x56 -> 20 um^2/s (symptom of loose membrane; transport calibration to be
  re-validated after part 2, since the x56=11.5 was set on the frozen membrane).
- Physical target (full-res mechanism): leaflets held by TAIL INTERDIGITATION —
  tail-ends ~5 A from midplane, 42% of tails cross into the opposing leaflet;
  heads at |z|~18 A (head-head ~37 A = correct DOPC thickness). This
  interdigitation is flexible-tail-driven and DIRECTOR-INDEPENDENT.
- The single-particle CGL cohesion is director-TIED (rigid conformers), so it
  vanishes at disorder. Fix: add a source-derived, director-robust inter-lipid
  cohesion (flexible-tail interdigitation from two-molecule source), validated
  against full-res. Then re-validate diffusion (x56) + re-tune mass-scale if
  needed.

## Update 9 (2026-07-16): Option A validated physically; base pair grid has singular Coulomb

Option A viability CONFIRMED: the director-robust inter-leaflet interdigitation
cohesion is real and capturable.
- Best interdigitated (transverse-offset) config gives -15..-23 E_up per pair,
  nearly constant across director tilt 0/35/55 deg (interdig_proto2/3).
- The current pair build MISSES it: (i) it never samples the transverse offset
  between the two lipids (only azimuth/bead-frame/conformer), and (ii) tau=25 is a
  soft-min that averages the cohesive configs away with repulsive head-on ones.
- Recipe: at the cross-leaflet mask, sample transverse offsets and project at the
  physical T0 -> recovers ~-18 E_up director-robust cohesion.

BUT a fast patch is blocked by a base-grid numerical artifact:
- energy_grid_raw at inter-leaflet bins is dominated by SINGULAR Coulomb
  (NC3+/PO4- collapse): values -1e5..-1e13 E_up at cos=-1 bins (79/6860 bins
  < -1e6 E_up). The 1e-6 nm guard doesn't tame Coulomb -q^2/r.
- So min(base, cohesion) = base everywhere on the mask (base is far more
  "attractive" numerically), and the raw grid can't be cleanly patched.
- The runtime interaction_param (log1p, clipped, smoothed) tames this, so the
  sim is stable, but the raw grid is not a clean quantity to modify.

Conclusion: Part 2's clean fix is a careful pair-table BUILDER re-derivation that
must (a) fix/soft-core the singular Coulomb consistently, (b) add transverse-offset
sampling for director-robust cohesion, (c) rebuild dopc.h5 + re-validate pure
bilayer tightness at S2~0.5, hybrid tilt, and x56 diffusion. This is a multi-cycle
force-field re-derivation, not a quick term.

## Update 10 (2026-07-16): Option A hits the clustering/tau tension (user's worry confirmed)

Corrected analysis (no soft-core; model handles Coulomb via reference=-1.89 +
log1p clamping). Real gap: inter-leaflet pair well is shallow (-1.9 E_up) vs the
genuine interdigitation cohesion (~-18 E_up), because the build never samples the
transverse offset and tau=25 soft-min washes the sparse cohesive configs.

Tested Option A end-to-end: physical-T0 offset-sampled inter-leaflet cohesion at
the cross-leaflet mask ONLY (kept tau=25 laterally), patched onto the
rotation-calibrated 1RKL:
- GOOD: bilayer TIGHTENED, intra-leaflet z-spread 16 -> ~6 A (toward full-res).
- BAD: severe CLUSTERING (22/36 occupancy cells empty), protein over-tilted
  35 -> 67 deg (lay flat in broken membrane), diffusion blew up.

Root of the failure = fundamental single-particle tension:
- Capturing cohesion needs a SHARP projection (low tau), which UNDER-counts the
  excluded-volume entropy (picks the one interdigitated offset, ignores the many
  clashing ones) -> loses the short-range repulsive wall -> collapse/cluster.
- The SOFT tau=25 correctly counts that entropy (no clustering) but washes the
  cohesion -> shallow well -> tilt bug.
- No single projection temperature recovers both cohesion AND excluded volume,
  because the offset/tail-conformation DOFs carrying both were integrated out.
  (User's clustering worry confirmed empirically.)

Implication: within the single oriented-particle representation, holding a fluid
bilayer with correct cohesion appears to require the missing DOFs. Likely
resolution = minimal multi-bead lipid (Option C: e.g. 3 beads) that carries
explicit excluded volume + interdigitation natively, OR a fundamentally better
projection that separates the attractive and repulsive offset ensembles.

## Update 11 (2026-07-16): Option A exhausted — single-particle @40ps can't hold a fluid tilted bilayer

Tried all single-particle routes to disorder-robust cohesion; all fail on a
coupled trilemma (tight bilayer + fluid disorder + no clustering + tilt-holding):
- Offset sampling (sharp T0): captures cohesion but removes the r-wall
  (offset lets lipids reach the well at any r) -> collapse.
- Nested projection (offset T0 + orientation tau25): still monotonic-attractive,
  no wall -> collapse.
- Targeted well from full-res U(r) (keep wall, deepen r>=1.55 only): TIGHTENED
  ok (z-spread ~5), but RE-ORDERED the membrane (S2 0.45->0.7) because the
  cohesion sits at tails-facing (ordered) orientations, so the protein re-drifted
  vertical (35->13), and it still clustered moderately (8-9/36 empty).

Root (definitive): the single oriented particle couples cohesion, orientational
order, and excluded volume because the decoupling DOFs (tail conformation,
transverse offset) are integrated out. Full-res decouples them via flexible tails
(cohere regardless of head-tail director). No projection choice recovers all.

Deeper design tension: decoupling REQUIRES tail DOFs = bonds, but bonds are
incompatible with the 40ps step (118x the ~0.34ps velocity-Verlet stability
limit; even constrained, nonbonded/angle limits keep MARTINI near ~40fs). So a
multi-bead lipid (Option C) that would hold the bilayer cannot run at 40ps
without a dual-timestep integrator.

CONCLUSION: within the single-particle @40ps design, the protein-tilt bug is
intrinsic. Part 1 (rotation calibration) is a real, keepable fix; part 2
(disorder-robust cohesion) is not achievable in this representation. Resolutions:
(i) dual-timestep multi-bead lipids (lipids ~40fs, protein 40ps); (ii) a minimal
extra internal CGL DOF (a tail director separate from the head-tail orientation)
that provides director-robust cohesion; (iii) accept the limitation.

## Update 12 (2026-07-16): Option 1 (tail-director DOF) refuted by full-res data

Premise of Option 1: a tail director t, more ordered toward the midplane than the
head director n, could carry the cohesion without re-ordering n. Full-res data
refutes it:
- head->tail director n: S2=0.538; tail director t (GL->MEAN tail-end): S2=0.446;
  combined-10-bead PCA: S2=0.054 (that was the SPLAY axis, a wrong measure).
- Checked the two-tail structure explicitly (DOPC has sn-1 + sn-2 tails): PER-TAIL
  directors are EVEN MORE disordered -- sn-1 A S2=0.303, sn-2 B S2=0.342, per-tail
  PCA 0.17/0.27. The two tails splay 48+-24 deg. Averaging two tails splayed +-24
  about the normal FLATTERS the order (0.446) vs either real tail (0.30). So the
  two-tail structure is NOT the reason Option 1 fails; correcting for it
  STRENGTHENS the refutation -- there is no ordered tail director, per-tail or
  averaged. The 48 deg splay is itself a conformational mode (on top of per-tail
  bend + reach): tails interdigitate by flexing/splaying, not by pointing in any
  ordered direction. No number of rigid directors can represent that.
- n-t angle 16.6+-8.7 deg -> coupling k_nt ~ 47 kJ/mol/rad^2 (moderate). So any
  order the cohesion imposes on t leaks back into n (16 deg stiffness) and
  re-drives the protein vertical -> no decoupling.

Interpretation: in full-res the director-robust cohesion is carried by TAIL
CONFORMATIONAL FLEXIBILITY (tails bend/reach to interdigitate regardless of
orientation), NOT by any director. A second orientation DOF cannot represent
bending; a per-lipid scalar cannot represent the pairwise/directional reach
without a bilayer reference. => no single-particle DOF (orientation or scalar)
captures it.

DEFINITIVE: the tilt-preserving inter-leaflet cohesion is intrinsically a
flexible-tail conformational, pairwise effect. Single oriented-particle @40ps
cannot reproduce it (Options: compaction-q, angular, density, offset/nested/
targeted-well cohesion, and now tail-director -- all exhausted). Resolution:
Option 2 (explicit multi-bead tails, dual-timestep: lipids ~40fs, protein 40ps)
or accept the limitation. Part 1 (rotation calibration) remains a real fix.

## Update 13 (2026-07-16): tilt bug re-diagnosed via controlled runs -- NOT cohesion

User asked for a last attempt to save the single-particle CGL. Re-examined with a
CORRECT tilt metric (CA-atoms-only helix axis; earlier PCA was polluted by NA/CL
ions -> gave false readings). Clean facts:
- Bug reproduced: single-particle 1RKL drifts 35 -> 6 deg. Full-res reference
  HOLDS ~53 deg (not 35).
- The drift happens while the membrane stays TIGHT (z-spread ~8) and ORDERED
  (S2 0.985) throughout -> the tilt bug is NOT the cohesion trilemma (that was
  about fluid+tight). It is the protein-lipid coupling / membrane relaxation.

Controlled short runs from the 34.5 deg start (freeze/vary one thing):
  E1 free positions, q free, ordered dirs        : 35 -> 6   (drift)
  E2 free positions, q FROZEN@0.43, ordered dirs : 35 -> 7   (drift) => compaction
     collapse is a SYMPTOM, not the cause (refutes the compaction hypothesis).
  E3c FROZEN positions, ordered dirs             : 35 -> 32  (HOLDS)
  E3  FROZEN positions, DISORDERED dirs (S2 0.5) : 35 -> 55  (rises to full-res!)
Conclusions:
- The drift REQUIRES membrane positional relaxation (frozen positions hold tilt).
- DISORDERED directors give the protein the CORRECT tilt preference (~55 deg).
  Ordered directors + free positions -> the membrane relaxes into a vertical-
  driving picket-fence.
  E4 free positions, UNFROZEN rotation (disorder): tilt 35 -> ~27 (HELD, not 6)
     BUT membrane DISINTEGRATED: z-spread 6->30, S2 crashed 0.5->0.15.

FINAL (single-particle, proven by controlled runs): the ONE director DOF controls
BOTH cohesion and protein tilt, requiring OPPOSITE values:
  ordered directors  -> tight membrane BUT protein driven vertical (via free-
                        position relaxation into a picket fence)
  disordered dirs    -> correct protein tilt BUT membrane disintegrates
No single director value gives both. The pair potential's free equilibrium is
loose+disordered (z 30, S2 0.15); only the frozen prep placement keeps it
ordered+tight -> the pair cannot prefer a tight fluid (S2~0.5) bilayer (= the
Option A cohesion trilemma, reconfirmed sharply). The only escape (director-
AVERAGED protein coupling that ignores the artificial order) needs a bilayer-
normal depth reference = forbidden ("no bilayer info in CGL FF"). => single
particle cannot be saved; multi-bead (separate tail-cohesion beads vs the head
the protein sees) is the physical decoupling. This CONFIRMS the Option 2 plan.

## Update 14 (2026-07-16): ghost-tail concept probe -- hybrid feasible, but needs explicit tail PACKING

Cheap proxy (g1_local): free positions (real dynamics), disorder ONLY first-shell
directors (16/275 within 15A of protein), keep bulk ordered+tight.
Result: tilt 34.5 -> wanders 10-32 (ends ~11); bulk z-spread STAYS TIGHT 8.8.
- POSITIVE: local modification does NOT break the bulk (vs E4 global disorder ->
  z 30). Hybrid/local resolution is FEASIBLE -- can keep the calibrated bulk.
- POSITIVE: removing the local director grain stops the MONOTONIC vertical
  collapse (E1 clean ->6; proxy oscillates instead).
- BUT: director disorder alone does NOT pin the tilt (protein wanders, no
  restoring force to 35-55). The tilt-pinning force is explicit tail PACKING
  (excluded volume): E3 pinned 55 via frozen-position scaffold; full-res pins via
  flexible tails packing around the protein. The proxy had no tails -> no pin.

IMPLICATION: ghost tails need the actual tail BEADS (excluded volume + hydrophobic
LJ), not just local director disorder. Testing that requires the Brownian
integrator + multi-bead particle/LJ infra = the SAME core as full multi-bead
(Option 2 P1-P2). So the core build is required either way; ghost-tail (local
resolution) becomes a deployment optimization on top, with adaptive-boundary
complexity that full (uniform) multi-bead avoids. Recommend: build the multi-bead
core (solves tilt via explicit tail packing everywhere), keep ghost-tail as a
later speed optimization.

## Update 15 (2026-07-16): why hybrid (ghost) schemes are hard -- time-scale mismatch

Full-res runs at --time-step 0.004 --integrator v with FULL stiff bonds (1250
kJ/mol/nm^2) and is stable ONLY because 0.004 is an UNCALIBRATED mechanical step
worth ~fs of real time (dt*sqrt(k/m)~0.003 << 2). It is SLOW (many steps/real
time), not fast. The CGL's 0.004 is a CALIBRATED 40 ps/step (mass-scale 0.05 +
GLE). So the SAME dt means ~fs (full-res) vs 40 ps (CGL) of effective dynamics.
=> a CGL-bulk + full-res-shell hybrid has a ~10^4x effective-timescale mismatch
at the interface (membrane tears). Any ghost/hybrid must re-calibrate the local
multi-bead region to 40 ps/step anyway (= the same Brownian work as full 3-bead),
AND full-res is a WORSE 40 ps Brownian unit (light beads -> ~2.5x stability margin
vs 3-bead ~55x). => uniform 3-bead is the robust choice; ghost = later speed
optimization only. Confirms Option 2 uniform-3-bead plan.

## Update 16 (2026-07-16): the hard-LJ-core wall blocks single-step 40 ps for ALL multi-bead

Implemented the Brownian (overdamped) integrator (src/martini_brownian.{h,cpp};
/input/brownian atom_index+friction; hooked into deriv_engine integration_cycle,
skips those atoms in Verlet; registered in main.cpp; builds clean). It is CORRECT:
- Verlet on the full-res bilayer .up = stable (confirms the .up + plumbing).
- Brownian is stable at small dt.
BUT at the 40 ps step it NaNs, for BOTH models:
- full-res (factor 4, gamma=0.0215): stable only at dt<=0.0001 T_up (~1 ps) => needs
  ~40 sub-steps per 40 ps.
- 5-bead (factor 11.2, gamma=0.168, 8x more friction): ALSO NaNs at 40 ps.
Root cause (diagnosed, not a bug): the physical 40-ps diffusion step (~0.1-0.6 A)
occasionally pushes a bead into the raw MARTINI LJ core (r^-12), where the explicit
overdamped step overshoots -> NaN. Initial config is clean (min r/sigma=1.09); it's
the dynamics. The CGL survives at 40 ps ONLY because its pair table is SOFTENED
(log1p/reference-clamped = no hard core).

DEEP TENSION: the hard excluded volume that PINS the protein tilt is exactly what
DESTABILIZES the 40 ps overdamped step. Tilt-pinning and single-step-40ps-stability
are incompatible with explicit dynamics + a hard core.

Options for a working 40 ps multi-bead membrane:
 (a) DUAL-TIMESTEP: sub-step lipids (~1 ps, ~40 inner steps) + protein at 40 ps.
     Physical (hard core kept, tilt pinned), robust; ~40x lipid force evals slower.
     = the original Option 2(i) dual-timestep.
 (b) IMPLICIT overdamped integrator: single-step 40 ps, hard core, unconditionally
     stable. Major build (Newton solve per step).
 (c) DISPLACEMENT-CAPPED Brownian: cap |dx| per step (integrator safeguard, NOT a
     potential change). Single-step 40 ps, hard core, cheap; small dynamics bias.
 (d) SOFTENED core (CGL-style): 40 ps stable but risks weak tilt-pinning; "soft-core".
The Brownian integrator is the shared basis for (a)/(b)/(c). Correction factor 4
for full-res / 14*4/N for N-bead confirmed (user directive).

## Update 17 (2026-07-16): IMPLICIT Brownian integrator WORKS -- full-res stable+fluid at 40 ps

User chose the implicit integrator. Implemented in src/martini_brownian.cpp. Journey:
- Full linearly-implicit (I+mob*H)dx=b via matrix-free CG: STABLE (no NaN) but froze
  (dx~0) because the LJ concave regions make A indefinite (pAp<0) -> CG breaks.
- Diagonal Hutchinson implicit: NaN'd (noisy diag underestimates core curvature ->
  huge core drift not suppressed -> overshoot).
- WORKING scheme: curvature-ALONG-THE-FORCE implicit. Per bead, one matrix-free
  probe measures curvature in the direction of its own force (the direction the stiff
  drift acts); step = b/(1+mob*hdir), hdir clamped >=0; + a per-step displacement cap
  (1.5 A) as a guaranteed safety net. Cost ~1 extra force eval/step.
RESULT (full-res DOPC bilayer, dt=0.004=40 ps, factor-4 friction gamma=0.0215):
  STABLE (5000 steps, finite), thickness relaxes 62->47 A, tail |S2| 0.97->0.52
  (FLUID), 2 leaflets intact. => full-res runs at a single 40 ps overdamped step.
This is the goal: reuse the validated full-res FF at 40 ps (no CG parameterization
risk). Cap/curvature is an integrator scheme (not a force-field change).
NEXT: (1) calibrate lateral D * 4 = physical (tune gamma; check cap activation isn't
biasing D); (2) correct area/thickness (NPT or right area); (3) VALIDATE 1RKL tilt at
40 ps with full-res+implicit; (4) verify sampling vs reference (implicit is approximate).

## Update 18 (2026-07-16): SOLVED -- full-res + implicit Brownian @40 ps HOLDS the tilt

1RKL hybrid (protein Verlet + full-res 14-bead lipids on the implicit Brownian
integrator, factor-4 gamma=0.0215, dt=0.004=40 ps):
  tilt 34.5 -> 47.6 deg (RISING toward the full-res equilibrium ~53 deg);
  lipid z-spread stable 10.2 -> 9.5 (tight membrane); no NaN; 122 s / 4000 steps.
Compare: single-particle CGL 34.5 -> 6 (collapses vertical = the bug); full-res
reference (small dt) -> ~53. So the implicit-Brownian full-res REPRODUCES the
full-res tilt physics at a single 40 ps step, reusing the validated FF with ZERO
force-field compromises. THE TILT BUG IS SOLVED via this route.

Path that worked: full-res model (correct FF) + overdamped dynamics (physically
correct at 40 ps) + implicit step (curvature-along-force + displacement cap) to
survive the hard LJ core. This sidesteps ALL the CG parameterization struggles
(5-bead head amphipathy/thickness) by reusing full-res.

Remaining (refinement, not blockers): (1) calibrate lateral D*4=physical (tune
gamma); (2) confirm the 1.5A displacement cap rarely activates (no dynamics bias);
(3) longer run to confirm tilt settles ~53 and stability; (4) speed (implicit adds
1 probe compute/step; optimize if needed); (5) 1AFO + regression.

## Update 19 (2026-07-16): diffusion calibration reveals implicit over-damps lateral D

Implicit @40 ps holds the TILT (structure correct) but OVER-DAMPS lateral diffusion:
- directional full-suppression: D*4 = 0.08 um^2/s (140x too slow).
- FDT semi-implicit (suppress drift by (1+mob*h), scale noise by 1/sqrt(1+mob*h)):
  D*4 = 1.2 um^2/s (9x too slow), stable+fluid (S2 0.41).
Analysis: D = kT/(gamma + dt*h) -> as gamma->0, D -> kT/(dt*h) = a CAP independent of
gamma. With true diagonal h~15, cap ~4 um^2/s (x4), still ~3x below target 11.5.
=> the drift suppression that gives 40 ps stability also caps the lateral diffusion.
Cheap per-bead curvature (directional/diagonal) also OVERESTIMATES stiffness (caged
lipid stiff in all dirs) -> extra over-damping. Only a FULL mode-selective implicit
(Newton-Krylov: soft eigenmodes -> A~1 -> free diffusion) avoids it, but that froze on
the LJ's indefinite (concave) Hessian.

FUNDAMENTAL: correct factor-4 lateral diffusion AND single-step 40 ps stability are in
tension for the stiff MARTINI potential. Options:
 (a) SUB-STEP explicit overdamped (dt_inner~1 ps, ~40 inner steps): correct diffusion
     (no suppression) + stable; cost ~40x lipid force evals (dual-timestep).
 (b) FULL Newton-Krylov implicit with indefiniteness handling (MINRES+trust region):
     correct diffusion + single 40 ps; major build, uncertain.
 (c) ACCEPT over-damped diffusion (tilt/structure CORRECT, lateral dynamics ~3-9x slow;
     report an effective correction factor rather than the physical 4).
Tilt is SOLVED regardless; this is only the lateral-diffusion timescale.

## Update 20 (2026-07-16): factor-4 diffusion at single 40 ps step is fundamentally blocked

Exhausted the "clever" single-step routes to raise lateral D to the physical target
(x4=11.5), all hitting the SAME wall (the hard-LJ stability mechanism caps D):
- mass-scaling: N/A (overdamped is mass-free; the cap D=kT/(N dt h) is set by potential
  stiffness + dt, not mass; CGL's mass-scaling works only via SOFT potentials).
- implicit (drift-suppress + FDT noise): D*4 ~ 1.2 (suppression damps the diffusive mode).
- capped-explicit (full noise + displacement cap, no suppression): D*4 ~ 2.0-2.4,
  FLAT vs gamma (the cap itself limits D: lower gamma -> bigger noise -> hits cap),
  and bonds over-broaden (std 0.9-1.0 vs 0.45).
Ceiling ~2-5 um^2/s (x4), 2-5x below 11.5. ROOT: to diffuse the physical ~0.4 A/step
past a hard LJ core without overshoot needs a stability measure (suppression OR cap),
and that measure caps the lateral diffusion. Provably not fixable by mass/friction.

ONLY escapes (each relaxes a user constraint): (a) sub-step small dt [user: no];
(b) soften the LJ core [user: wary]; (c) accept ~2-9x slow [user: no]; (d) FULL
mode-selective implicit (regularized Newton-Krylov: free the soft diffusive eigenmode,
suppress stiff radial modes) -- the only route meeting all constraints, but a major
uncertain build and likely still ~1.5x off from the regularization.
TILT IS SOLVED regardless (55.5 deg, matches full-res); diffusion is the open item.

## Update 21 (2026-07-16): diffusion target confirmed + soft-core retracted (user was right)

Soft-core RETRACTED: it changes the physical diffusion (caging barriers set lateral
D via Kramers; softening lowers the barrier -> faster hopping). The 40 ps overshoot
is at the CONTACT/caging distance (r~sigma-step), so softening there twists transport,
not just numerics. Also breaks the MARTINI x4 calibration (tied to the hard-core FF).
User's skepticism correct.

TARGET CONFIRMED (cg_lipid_potentials.tex): 11.5 um^2/s = experimental DOPC at 303 K
(Dickson2014/Lipid14; other refs 6.5-9.5). So DOPC microscopic lateral D IS ~10, not
~1 (I was wrong). Factor for full-res = 4 (56=14x4, no lumping). So D_sim target =
11.5/4 = 2.875 A^2/T_up. My overdamped full-res gives 0.589 -> ~5x gap (real).
Question: gap in the INTEGRATOR (my overdamped throttle for 40ps stability) or the FF
(hard-core caging)? -> measuring native full-res (small-dt Verlet) diffusion now.
Time unit: 1 reduced T_up = tau0 = sqrt(m_up A^2/E_up) = 0.203 ps; dt=0.004 = 0.81 fs
mechanical. (CGL's "0.004 T_up = 40 ps" is the CALIBRATED time via mass-scale+GLE.)

## Update 22 (2026-07-16): CHEAP fix validated -- M~5 substep of overdamped hits target D

Gap quantified: lipid D is ~5x below target (D_sim 0.59 vs 2.9). It's CAP-limited, not
caging-fundamental: shrinking the step 5x (dt=0.0008) UN-CAPS it. At gamma=0.006:
 40 ps single step: D*4 = 2.0 (capped)
 dt/5 step:         D*4 = 8.07 (un-capped) -> near target 11.5; lower gamma to ~0.0042
So M~5 sub-steps of the CURRENT overdamped integrator (no SHAKE -- overdamped is stable
at any inner dt) + a gamma calibration reaches the target. This is the CHEAP path
(~2-5x slowdown), NOT the M~1000 native substep (that was for the inertial mechanism).
Gives the correct diffusion COEFFICIENT + ~8 ps dynamical resolution; M tunable for finer.
Implementation: RESPA/dual-clock -- lipid forces (martini_potential+dist/angle_spring)
inner M x via existing integrator-level mechanism (level 0 fast); protein core outer 1x
(level 1 slow). Bounded build (the level mechanism already exists in deriv_engine).
Also: neighbor/cell lists on the lipid nonbonded (currently all-pairs O(N^2)) speeds each
eval further. TILT already solved.

## Update 23 (2026-07-17): literature survey confirms sub-stepping (no single-step alternative)

Surveyed single-step (non-MTS) large-timestep methods for a HARD-CORE bilayer at 40 ps.
Controlling fact: hard dry-MARTINI bilayer explicit stability ~20-40 fs (Arnarez2015);
40 ps is ~1000x beyond. Classical stiff-integrator toolbox (SHAKE/LINCS constraints,
HMR, implicit/exponential/IMEX, LN Schlick, mollified impulse, BAOAB, preconditioning)
buys only 2-10x -> three orders of magnitude short. Where A-stable at 40 ps
(implicit/exponential), they damp/distort the cage-rattle modes that PRODUCE diffusion
-> fail on D. Verdict families:
- GLE / Mori-Zwanzig memory (Ceriotti2010, Li2015, Klippenstein2021): most principled
  single-force-eval route, reproduces correct D + T-dependence by construction -- BUT the
  memory kernel must be integrated on the sub-40-ps timescale -> re-introduces sub-stepping
  in disguise. Converges back to our choice.
- Overdamped/Brownian single 40-ps step: under-diffuses ~5x (our finding); matching D
  needs friction tuning = twisting.
- ML learned propagators (Score Dynamics 10-ps steps JCTC2024; Hamiltonian Flow Maps
  ICML2026): genuinely large steps, can preserve kinetics, but need per-system training
  data (generated by sub-stepped MD), non-universal, unproven for condensed bilayers.
  Timewarp = sampler, not kinetics.
- DPD/soft-CG: large steps only by softening the core (forbidden) + tuned transport.
- aMD/GaMD/hyperdynamics: accelerate rare-barrier crossings, not diffusion; lateral
  diffusion is not a rare-barrier process. MSM/kMC: surrogate from sub-stepped MD, not an
  online integrator for coupled protein+bilayer.
CONCLUSION: RESPA sub-stepping is physically correct and effectively unavoidable for a
hard-core bilayer at 40 ps. Our implementation is the right call.

## Update 24 (2026-07-17): sub-stepping calibrated + multi-T validated

FINAL working config (full-res bilayer, dense S=8.5 grid): RESPA sub-stepping M=10,
total-step cap 2.0 A, gamma_ref=0.0035 (E_up*T_up/A^2), t_ref=0.8647, Ea=10.3 E_up.
Calibration at 303 K: D*4 = 11.24 um^2/s (target 11.5, ratio 0.98) -- dead on.
Multi-T (SAME gamma_ref + Arrhenius gamma(T), NO per-T retuning):
  280 K: D*4=3.08 (law 4.39, ratio 0.70)
  303 K: D*4=11.24 (law 11.50, ratio 0.98)
  333 K: D*4=35.55 (law 33.51, ratio 1.06)
-> Arrhenius gamma(T) transfers across a 12x diffusion range; caging factor ~T-independent
(low-T 30% deviation = stronger low-T caging + noisy slow-D). gamma(T) trend is physical:
gamma DECREASES with T (D increases, Arrhenius). This IS "tuning gamma" -- user-accepted,
with a clear sensible T-trend, no lookup table.
Cost: M=10 => 10x lipid force-evals/outer step (honest cost of the dense hard-core cage;
protein core stays 1x/outer, level-0 optimization is future P6).
User accepts gamma-tuning; note single-step gamma-tuning is FUTILE (cap-limited, D flat vs
gamma: 2.35/1.98/2.00 at gamma=0.0215/0.010/0.006) -- sub-stepping is what gives gamma leverage.

## Update 25 (2026-07-17): CRITICAL -- tilt NOT held under correct (fast) diffusion

1RKL on-system, sub-stepping M=10 + calibrated gamma(T) (D*4=14.7, on target), stable:
tilt trajectory 34->24->18->15->20->20->29->28->26->13->6 deg (mean last 30% = 17.3),
vs full-res reference ~53 and initial ~35. The tilt DROPS toward vertical -- the ORIGINAL
BUG behavior.
KEY INSIGHT: the earlier "tilt solved" (34.5->55.5) was measured with the SLOW/over-damped
diffusion (implicit, D*4~0.1-1.2). With correct FAST diffusion the membrane can relax and
the tilt is lost. => tilt-holding was likely a SLOW-DIFFUSION ARTIFACT; the physics is NOT
actually settled once dynamics are correct.
LEADS: (a) SC-env (protein side chain <-> lipid) coupling may be too weak / a raw LJ rather
than a proper PMF (user flag #1) -> can't anchor the tilt when the membrane is mobile;
(b) reconsider what physically holds a TM tilt in full-res dry-MARTINI.
ACTION: PAUSE the CGL cleanup (removing CGL is premature if the tilt isn't held). Verify the
regression (longer run + compare to a small-dt reference trajectory over matched real time),
then diagnose (SC-env PMF first).

## Update 26 (2026-07-17): SC-env is already a PMF; tilt loss is likely the diffusion-enabled ratchet

SC-env investigation: the ACTIVE full-res path ALREADY uses a proper PMF for SC-env
(martini_sc_table_1body + _run_sc_task/_build_sc_table_group: rotamer-weighted, multi-bead-
averaged, tempered Boltzmann free energy, radial x angular(cos theta)). So user #1 is
SATISFIED for SC-env. BB-env = raw single-bead MARTINI LJ, which is the correct NATIVE
backbone treatment (4 atoms->1 bead standard MARTINI), not a multi-atom-collapse case.
=> tilt loss is NOT from SC-env.
HYPOTHESIS: tilt loss = the original RATCHET, now enabled by correct fast diffusion. With
fast lipid diffusion the membrane repacks to the protein's current tilt each step, walking
the coupling-energy minimum toward vertical. Slow/over-damped diffusion froze the membrane
-> tilt held (the 55.5 artifact). Full-res reference (~53) ran too little REAL time to show
the ratchet. So the tilt bug may live in the FF (protein-lipid coupling / membrane repack),
exposed only under correct dynamics -- OR ~20 deg is the true equilibrium and 35 was the
(non-equilibrium) start. DISCRIMINATOR: slow-diffusion control (should reproduce ~55) vs
fast (drops) confirms the diffusion-speed dependence; then it's a FF physics question.

## Update 27 (2026-07-17): PIVOTAL -- ~17 deg tilt is PHYSICALLY CORRECT; the "tilt bug" may be a misdiagnosis

Literature (WALP/KALP, Im-group PMFs, ss-NMR/fluorescence): a single WELL-MATCHED TM
alpha-helix in a thick fluid PC bilayer (DOPC ~27-29 A hydrophobic) sits NEARLY VERTICAL,
equilibrium tilt ~5-25 deg. ~17 deg is squarely correct. 30-40 deg is a POSITIVE-MISMATCH
signature (long helix / thin bilayer like DMPC), NOT a matched helix in DOPC.
- 1RKL = Ost4p: 36-res KINKED single TM helix, ~18 hydrophobic res ~ matched to DOPC ->
  should be nearly vertical. No experimental bilayer tilt (solved in micelles).
- 1AFO = glycophorin A dimer: each monomer ~15-20 deg from normal (nearly vertical);
  consistent with "inserted near-vertical".
- Tilt is a SLOW DOF (~us to converge), broad +/-10-15 deg distribution.
REFRAME: the sub-stepping (correct fast diffusion) relaxing 1RKL to ~17 deg is likely the
CORRECT physical equilibrium. The 55 deg "held" tilt (slow/over-damped diffusion, frozen
membrane) was the ARTIFACT. The original "protein tilts to vertical" bug may have been the
protein CORRECTLY relaxing a non-equilibrium ~35 deg start toward its matched-helix vertical
equilibrium -- i.e. NOT a bug. => the whole tilt-fix effort may have chased an artifact.
DECISIVE TEST (user-requested): standard wet MARTINI 2.2 + explicit water + DOPC reference
of 1RKL in GROMACS -> if it ALSO relaxes to ~vertical, ~17 deg is confirmed correct and the
full-res dry-MARTINI sub-stepping is right. Caveat: need us-scale + broad distribution.

## Update 28 (2026-07-17): RESOLVED -- wet MARTINI 2.2 confirms 1RKL is nearly vertical; tilt "bug" was a misdiagnosis

Wet MARTINI 2.2 ground truth (user's CHARMM-GUI system: 1RKL + 66 DOPC + 812 W + ions,
GROMACS 2026.3, min+eq+30 ns prod): 1RKL tilt relaxes to mean(last 50%) = 10.3 +/- 5.4 deg
(nearly vertical). Independent of Upside/dry-MARTINI.
THREE-WAY AGREEMENT that 1RKL is nearly vertical in DOPC:
  wet MARTINI 2.2 (explicit water):        ~10 deg
  dry-MARTINI full-res + sub-stepping:     ~17 deg
  literature (matched TM helix in DOPC):   ~5-25 deg
The 55 deg "held" tilt ONLY appeared with artificially frozen (slow/over-damped) lipids.
CONCLUSION: the original "protein tilts to vertical" was CORRECT physics (a matched single
TM helix relaxes toward vertical), NOT a bug. The tilt-fix effort chased a slow-diffusion
artifact. The full-res + sub-stepping model is VALIDATED: correct tilt AND correct diffusion.
VTF trajectory for the user's boss: ~/Downloads/charmm-gui-8427092652/gromacs/martini22_run/
1rkl_martini22_wet.vtf (+ step7_production.xtc/.tpr, step6.6_equilibration.gro).
=> CGL cleanup UNBLOCKED (physics validated).

## Update 29 (2026-07-17): CORRECTION -- dry-MARTINI+sub-stepping DISPERSES membrane proteins; earlier "validated" was premature

Measuring 1RKL tilt with a CONSISTENT backbone-only (BB) axis (user caught that wet was
measured with SC beads included, skewing the PCA) revealed:
- Initials MATCH (dry 34.2, wet 34.3) -- same insertion. Good.
- BUT both dry runs DISPERSE the protein-membrane system under sub-stepping:
  * tilt_impl.up (=stage_7.0 Upside prep, protein Verlet + lipid Brownian): protein z-gap
    0->32 A, lipid z-span ->130 A over the run.
  * CHARMM-GUI-geometry dry build (all-beads Brownian): protein z-gap 0->87 A, tilt 34->69.
  => NOT a bad config: our dry-MARTINI + sub-stepping does NOT stably hold a membrane protein.
  Pure bilayer WAS stable, so it's the protein+vacuum hybrid (likely dry-MARTINI lipid
  evaporation into the implicit-solvent vacuum + protein diffusing out under fast dynamics;
  wet MARTINI has water to hold both).
- My earlier "dry reproduces vertical / physics validated" (Update 28) was measured on these
  DISPERSING runs -> RETRACTED. The dry Upside model is NOT yet validated for membrane proteins.

SOLID GROUND TRUTH (wet MARTINI 2.2, CHARMM-GUI + GROMACS, stable, matches literature):
  1RKL single helix: BB tilt -> ~10 deg (vertical).
  1AFO GpA dimer: monomers 19-22 deg, dimer axis ~6 deg (vertical). Matches GpA NMR/MD.
  => the tilt "bug" (protein going vertical) IS the correct physics; the 55 deg was the
  slow/frozen-membrane artifact. That conclusion stands (from the wet ground truth + lit).

OPEN ENGINEERING PROBLEM: stabilize the dry-MARTINI hybrid under sub-stepping (anti-
evaporation box / membrane restraint / proper dry-MARTINI setup / check the SC-env & BB-env
coupling strength holds the protein without water) before it can simulate membrane proteins.
Figures: ~/Downloads/1rkl_tilt_comparison.png (wet valid + dry-dispersal shown),
1afo_tilt_martini22.png (wet). Wet VTFs saved under ~/Downloads/charmm-gui-*/gromacs/martini22_run/.

## Update 29 (2026-07-17): dry-MARTINI 1RKL apples-to-apples (Option 2) + M=10 correction

Built a dry-MARTINI 1RKL-in-DOPC run directly from the CHARMM-GUI structure (protein
elastic network + PROA_P.itp bonded, 12-bead martini_v2.0 DOPC, DRY v2.1 nonbonded matrix,
overdamped RESPA brownian) as an apples-to-apples dry-vs-wet comparison on the SAME structure.
- STRUCTURE: step5_charmm2gmx.pdb is PRE-minimization (tail clashes r<0.6A, E~6e11) and blows
  up. Must start from step6.6_equilibration.gro (the equilibrated frame the wet prod started
  from; box 47.955 = 69.7 A^2/lipid). Then minimize under the dry FF (E -3146 -> -6260).
- INTEGRATOR CORRECTION: Update 24's "M=10 stable, D*4=11.24" was on the SOFTENED CGL pair
  table. On the RAW hard-core LJ (full-res/dry), M=10 DISPERSES the bilayer (z-span 48->103,
  energy spikes): the overdamped step is stable near the LJ min (mob*K~0.45) but unstable on
  the repulsive wall r<4.5 (mob*K~15) -> 2A cap fires -> energy injection -> heating. Need
  n_substep~40 (dt_inner 1e-4) so the cap essentially never fires (matches Update 16's "~40
  substeps"). KEY: D = 2 kT dt / gamma is INDEPENDENT of M, so raising M preserves the
  diffusion calibration exactly and only fixes hard-core stability. At M=40: membrane intact
  (PO4 head-head ~39A, z-span plateaus ~65A), protein embedded (z-gap ~2A), stable 6000 steps.
- TILT: dry 1RKL BB tilt starts 33 deg, samples ~10-32 (mean ~28 over 120ns), consistent with
  a matched TM helix (wet MARTINI 2.2 ref ~10, dry sub-stepping ref ~17, lit 5-25). Confirms
  1RKL is near-vertical/moderately tilted in DOPC, NOT the artifactual 55 deg (frozen-membrane).

## Update 30 (2026-07-17): dispersal ROOT-CAUSED + FIXED

The dry-MARTINI hybrid dispersal has TWO sub-causes, both now understood:
1. DEFORMATION: the protein was VERLET at 40 ps, but its stiff protein-lipid LJ collisions
   (from the Brownian, fast-moving lipids) kicked the backbone apart (BB bonds 3.4->4.6 A,
   Rg 12.6->24). The lipids had the sub-step+cap; the Verlet protein did not. FIX: sub-step
   the protein too (put it in /input/brownian) -> resolved (protein stays folded, Rg ~14).
2. DRIFT/EXPULSION: even folded, the protein diffuses OUT of the dry (implicit-solvent)
   membrane because there is no water to confine its soluble parts -> no z-confinement well.
   High friction only slows it (kinetic, not thermodynamic: z-gap 87->27 with x15 friction).
   FIX: a light z-anchor (restraint_position, z-only k=10 on ~3 central BB beads) supplies the
   confinement water would give, while leaving the tilt free.

C++ change (martini_brownian.cpp): optional per-bead friction_scale (aligned with atom_index);
gamma_i = gamma(T)*scale_i. Protein scale ~15 (protein-appropriate ~1 um^2/s x4 diffusion),
lipid scale 1 (calibrated). Both sub-stepped (M=10).

RESULT (dry, all sub-stepped, protein x15 + z-anchor): protein FOLDED (Rg 14) + EMBEDDED
(z-gap 0-5) + tilt mean(2nd half) = 21 +/- 7 deg -> PHYSICAL near-vertical, consistent with
wet MARTINI 2.2 (~10) and literature (5-25). Dispersal FIXED.

Note: the z-anchor is a pragmatic stand-in for the missing implicit-solvent confinement of the
protein; the deeper physical fix is an implicit-membrane solvation term (z-dependent per-bead-
type transfer free energy). Recorded as future work.

## Update 31 (2026-07-17): implicit single-step timescale RE-TESTED on raw full-res -- capped ~30x low

User asked to measure the timescale the IMPLICIT integrator can reach on the full-res DOPC
bilayer. Re-added a GATED implicit branch to martini_brownian.cpp (curvature-ALONG-force
suppression via one matrix-free probe: drift/(1+mob*h), FDT noise/sqrt(1+mob*h), 2A cap),
ran single-step M=1 on the raw hard-core 14-bead bilayer (128 lipids), then REVERTED the
branch (production RESPA integrator clean, rebuild passes).
RESULT: stable (finite) at one 40 ps step, but lateral D*4 = 0.38/0.42/0.41/0.26 um^2/s at
gamma = 0.0005/0.0035/0.010/0.0215 -- essentially FLAT across a 40x gamma range. That is the
intrinsic-cap signature D = kT/(dt*h) (set by curvature h + dt, NOT gamma): the drift
suppression that buys 40 ps stability also throttles the diffusive mode. Ceiling ~0.4 um^2/s,
~30x below target 11.5, un-tunable. Over-damping also freezes/thickens the membrane (PO4
thickness 55-67A at low gamma = the frozen-membrane artifact). RESPA M=40 ref (short 4000-step
run) = 6.3 (calibrated value 11.2). CONFIRMS Updates 19-20: the implicit is a STRUCTURE-only
integrator (holds tilt, freezes diffusion), NOT viable for correct single-step dynamics.
Reaffirms Option 3 (RESPA sub-step lipids + Verlet protein via BB-env->PMF).

## Update 32 (2026-07-17): directional-coupled implicit + GLE reconsidered -- single-step mode-selection fails

User: (a) try coupling implicit+overdamp; (b) reconsider GLE+damping. Tested (a): a GATED
directional-coupled branch (suppress ONLY the radial/into-wall component implicitly, leave the 2
TANGENTIAL/diffusive components as FREE explicit overdamped), single-step gamma-scan on the raw
full-res bilayer, then reverted.
RESULT: D*4 = 6.7/32.0/17.2/6.2 at gamma=0.0005/0.0035/0.010/0.0215 -- NOT flat (cap removed, the
tangential-free idea DOES un-throttle D), BUT the membrane DISINTEGRATES (z-span 134-182A vs box
74A; PO4 "thick" 22-39A meaningless). The high D is dispersal, not in-membrane diffusion. Root:
the cage is MULTI-directional; freeing the tangential lets beads overshoot SIDE neighbors -> bilayer
comes apart. So single-step gives a clean either/or: full-suppress = stable+frozen (D*4~0.4);
radial-only = correct-ish D + disintegration. No middle.

GLE+damping (b) reconsidered (prior Update 23 only covered GLE-alone). Framing: TWO independent
barriers at 40 ps -- (1) STABILITY = conservative position overshoot into r^-12 wall (set by force*dt);
(2) DIFFUSION = transport coeff (set by friction/noise). GLE (colored friction, +/- extra damping) is
a barrier-2 method: it reproduces correct D + VACF + T-dependence by construction, but friction acts
on dissipation, NOT the conservative overshoot -> it does NOT touch barrier-1. Frequency-selective GLE
(high friction on stiff modes, low on the diffusive mode) is the TEMPORAL analog of the coupled
spatial mode-split; at a single 40 ps step it cannot resolve the frequencies (memory time ~ps << 40ps)
-> reduces to Markovian overdamped (=5x/cap) OR unstable aux vars. All three mode-selection routes fail
their own way: spatial(coupled)->disintegrate; frequency(GLE)->Markovian; Hessian(Newton-Krylov)->froze
on indefinite LJ Hessian. COMMON ROOT: separating stiff from diffusive modes REQUIRES resolving the
sub-40-ps structure = sub-stepping. CGL precedent: GLE ran single-step at 40 ps ONLY because CGL was
SOFT (barrier-1 gone). GLE+hard-core still needs sub-stepping. GLE's real value = a colored friction
ON TOP of the sub-stepping (better VACF); the current Markovian sub-stepping already gives correct D+T.
=> Option 3 stands (RESPA sub-step lipids + Verlet protein via BB-env PMF).

## Update 33 (2026-07-17): BB-env -> PMF implemented as smoothed prep-only grids

Option 3's protein-side half is done: BB-env (backbone<->lipid) is a thermally-averaged PMF, so a
single-step Verlet protein is not impulsively kicked by the RESPA-substep lipids. NO C++ change --
reuses MartiniPotential + its existing BB-proxy gradient projection; only prep grids change.
- BB-env carrier confirmed: martini_potential.cpp SKIPS N/CA/C/O carrier<->env pairs (609-623) and
  routes the interaction through the virtual BB proxy (atom_name "BB"), projecting its gradient onto
  the 4 carriers. So the smoothing target is exactly the (BB proxy, lipid) pairs.
- Classify by IDENTITY not type: the BB proxy and its N/CA/C/O carriers share the same MARTINI bead
  type+charge (P5/Qd/Qa...), and lipid beads reuse Na/Qa too. is_bb_proxy = (particle_class==PROTEIN
  & atom_names=="BB"); is_lipid = particle_class==LIPID. A boolean folded into the coeff key gives
  BB-env pairs their own coeff rows (marked by bb_env_coefficient_mask); everything else is untouched.
- Smoothing math: U_bb(r) = -kT ln <exp(-U_raw(|r zhat+delta|)/kT)>, delta~N(0,sigma_bb^2 I),
  sigma_bb=UPSIDE_BB_ENV_RMSF (def 1.0 A), kT=0.8647 E_up. Isotropy reduction -> Gauss-Hermite(delta_z)
  x Gauss-Laguerre(perp radius), order 64, deterministic + SMOOTH (max |2nd-diff| ~1.5e-4, no MC
  noise -> clean force spline). U_raw = exact build_particle_h5 LJ+Coulomb with max(r,0.1 sig) clamp.
- Properties verified: sigma_bb->0 reproduces raw to 1e-8 rel; core softened (charged BB-env row raw
  14748 -> PMF -0.93 E_up @3A); tail -> raw (residual is the physical Mori-Zwanzig term ~ -15 A
  sigma^2 / r^8 for LJ and the 1/r Jensen shift for Coulomb, ~1-2% at cutoff, vanishing with r). The
  raw grids carry NO tail-shift, so the PMF adds none (matches the codebase convention).
- Master parity: pure DOPC bilayer produces byte-identical coefficients/combined_energy_grids/
  coeff_indices vs pre-change (reconstructed pre-change lib, diffed); 0 BB-env rows.
- Smoke (1RKL hybrid, protein Verlet + lipids brownian M=10, 300 free steps): smoothed Rg stable
  12.8, hbonds 13-30; unsmoothed (RMSF=0.01) degrades -- hbonds 5-10, Rg drift 12.8->13.3,
  protein_potential +77..+456. So the smoothing measurably protects the backbone.

Lesson: the `--minimize` upside calls legitimately print `avg_kinetic_energy/1.5kT nan` and a
`inf us/systems/step` timing (no velocities, --duration 0). Do NOT read those as a blow-up; check the
actual DYNAMICS `step N / M ...` lines (finite temperature/Rg/potential) for stability.

## Update 34 (2026-07-17): drift CAUSE = the BB-env PMF (my twist); raw FF confines correctly. + dt was non-standard.

TWO findings resolving the "protein moves out of bilayer" problem:
1. TIMESTEP: the dryMARTINI workflow ran the protein at dt=0.004 (inherited from the abandoned CGL) =
   HALF the Upside standard 0.009 (04.HDX/`--time-step` default). A slowdown that partly masked the drift.
   Switched to dt=0.009. Lipid D=kT/gamma is dt-INDEPENDENT, so keep dt_inner=0.0001 by scaling M 40->90;
   gamma unchanged, D unchanged. Protein + lipid OUTER dt both = 0.009.
2. DRIFT CAUSE (isolated at dt=0.009, ~1320 steps ~130 ns): flex+RAW FF STAYS (z-offset drift +0.6 A,
   folded Rg 12.5, embedded all frames); flex+PMF DRIFTS (-5.4 A). => the BB-env PMF SOFTENING is the
   drift cause; the RAW dry-MARTINI FF confines the protein correctly (matches the rigid test + the user's
   prediction: correct FF -> correct position, no confining potential needed).

REINTERPRETATION: the "kick" that motivated the BB-env PMF (O2) was almost certainly an UNDER-RESOLVED-lipid
artifact -- the old M=10 runs let the displacement cap fire -> energy injection -> kicked the protein --
NOT the raw BB-env. At well-resolved M (dt_inner~0.0001), the raw BB-env folds AND confines with no kick.

FIX (satisfies all user constraints): REMOVE the BB-env PMF (O2) entirely. Production = raw dry-MARTINI FF
+ well-resolved lipids (M so dt_inner~0.0001) + protein single-step Verlet at standard dt=0.009. No PMF,
no potential, no slowing, no protein sub-stepping; position emerges from the correct FF. Verifying long-term
stability (flex_raw extension) + re-confirming lipid D at 0.009/M=90 before removing the O2 code.
OPEN: tilt ~40 deg (stable, above wet ~10); insertion tilt, us-slow to relax -- assess with long runs.

## Update 35 (2026-07-18): single-step design space exhaustively closed (incl. reflecting BD)

Systematically evaluated + (where testable) empirically ran every "on top of damping" route to
correct single-step lipid diffusion on the HARD core. All fail against one invariant: correct D
fixes a per-bead kick ~2 Ang (~sigma/2) that a hard r^-12 core cannot absorb in one step.
- soft-core alone: alpha~1 -> contact energy 0 -> collapse.
- soft-core + damping (joint alpha-gamma scan): mild alpha (0.3-0.5) stays ~15x too slow at any
  gamma; only heavy alpha(0.7, excluded volume gone) approaches D via softened barriers (Kramers
  twist). No "just enough" point.
- implicit (drift suppression): D*4~0.4, flat vs gamma (throttled diffusive mode).
- directional/coupled implicit (radial suppress, tangential free): un-throttles D but bilayer
  disintegrates (multi-directional cage).
- mass-scaling: rescales dt AND dynamics together -> no timescale compression; light mass ->
  unstable; heavy -> ballistic/slow. (CGL mass-scaling worked only because CGL was soft.)
- GLE/colored noise: reproduces D+T by construction but memory kernel needs sub-40ps integration
  -> sub-step in disguise / Markovian limit at single step.
- REFLECTING BD (last untested; ran it): NaN. The force probe needed to detect the overlap
  overflows on the deep r^-12 overlaps a 2 Ang single-step kick creates in a dense bilayer; and one
  reflection off the NET force cannot handle the multi-neighbor cage. Proper (geometric,
  per-neighbor, iterative) reflection = event-driven collision resolution = sub-stepping in disguise.
CONCLUSION: single-step correct-timescale on a hard-core bilayer is impossible without giving up
the hard core (collapse/twist), universality (ML surrogate + training data), or the full-res model
(CGL). Sub-stepping (M~90, hard core intact) is the physical answer; the productive move is to make
it cheap (neighbor lists [done]; but the level-0 inner refresh was later found INEFFECTIVE, see
Update 36), not to eliminate it.

## Update 36 (2026-07-18): level-0 lipid-only inner refresh gives NO speedup (reverted)

Implemented the P6 optimization: inner Brownian substeps 1..M-2 call compute(DerivMode,0) (skip
the protein FF, integrator_level 1, which is fixed during substeps); the last substep does a full
compute so the protein Verlet update stays exact. Proved BIT-IDENTICAL (1RKL hybrid, old vs new,
identical incl NaN onset -- a lipid bead's sens never reads a protein-FF node, so the level-0 force
equals the full force).

BUT differential timing (20 vs 200 steps, startup cancelled) = 0.99x: 352 ms/outer-step for BOTH.
The inner loop is dominated by martini_potential (the lipid nonbonded spline compute, level 0,
~3.9 ms x M=90, recomputed regardless because lipids move every substep); the protein FF is
negligible for an embedded TM helix. Force-splitting saves nothing here.

-> REVERTED (rebuilt binary bit-identical to pre-change). The multi-step cost is M x (lipid
nonbonded). The only real speed levers: (a) REDUCE M (min-M; confounded on toy bilayers, needs a
proper equilibrated run), or (b) make martini_potential cheaper (cell list already done). Do NOT
re-attempt force-splitting for the hybrid -- the protein FF is not the bottleneck.

## Update 37 (2026-07-18): min-M determined CLEANLY = 45 (2x speedup) via cap-firing scan

M=90 (dt_inner=1e-4) was conservative. Determined the minimum stable M using the DISPLACEMENT-CAP
FIRING RATE (per-bead, per-substep) as the signal -- it directly measures overdamped overshoot into
cores and is immune to the undulation/dispersal confounds that ruined the earlier thickness/frac
scans. Added a temporary cap counter; built + equilibrated a dense 128-lipid DOPC bilayer (SPACING
8.1 ~0.66 nm^2/lipid) at M=90 (cap_frac 6e-6, intact frac 1.00, potential -11136 stable), then
scanned M from the equilibrated state:

  M   dt_inner  cap_frac  thick  pot_last
  90  1.0e-4    0         30.2   -11136
  60  1.5e-4    0         30.3   -11033
  45  2.0e-4    3e-6      29.8   -11031   <- clean; structure/energy == M=90
  30  3.0e-4    1.8e-3    28.8   -10616   <- knee: cap jumps 600x, thinning begins
  20  4.5e-4    1.2e-2    25.4    -8768
  <=15          collapses (thickness 22->19, potential -> positive by M<=7)

min-M = 45 (dt_inner 2e-4) -> 2x speedup with D/structure/energy preserved (cap 3e-6 negligible).
M=30 (3x) is past the knee (measurable thinning + cohesion loss) -> not recommended. Temporary
counter REMOVED; binary bit-identical to pre-diagnostic. Production default still M=90; lowering to
45 is a 2x win (M is a STABILITY knob, D is M-independent) -- pending confirmation + a spot-check on
the real 1RKL hybrid (protein-lipid contacts may pack slightly tighter, but 45 has margin to the
M=35 knee).

## Update 38 (2026-07-18): min-M=45 does NOT transfer to the hybrid; the protein interface needs M~90

Per the "verify on the hybrid before changing the default" decision, tested M=90/75/60/50/45 on a
current-setup, equilibrated 1RKL hybrid (raw BB-env, dt=0.009, from md_1rkl.up's last frame). The
protein-lipid interface (BB-env/SC-env forces are far steeper than lipid-lipid) makes the lipids
cap-fire ~1.7% even at M=90 (vs ~0 on the pure bilayer), and reducing M degrades monotonically:

  M   cap_frac  fracL  pot(0->400 steps)
  90  1.7e-2    0.94   -13898 -> -14350   (stable, cohesive, protein embedded Z~1)
  75  2.5e-2    0.84   -13898 ->  -6894   (cohesion halved)
  60  3.5e-2    0.76   -13898 ->  -4490
  50  4.8e-2    0.72   -13898 ->  -4915
  45  5.7e-2    0.67   -13898 ->  +6867   (dispersed; protein flips midplane; potential positive)

CONCLUSION: for the protein-embedded system min-M ~ 90 -- the pure-bilayer 45 is optimistic because
it lacks the stiff protein interface. PRODUCTION KEEPS M=90 (no safe reduction with a protein
present). The bulk lipids far from the protein are over-resolved at 90 (fine at 45); a
spatially-varying M (fine near the protein, coarse in bulk) or an inertial-Langevin/BAOAB
single-step scheme are the only further-speedup routes -- future work. This is exactly why the
"verify on the hybrid first" step mattered: the toy-bilayer optimum would have been wrong for production.

## Update 39 (2026-07-18): LITERATURE -- inertial Langevin (g-JF/BAOAB) single-step is the way off sub-stepping

Literature review on "correct D, single step, hard cores, on top of damping, no multi-step":

- ROUTE (strongly supported): replace OVERDAMPED BD with INERTIAL LANGEVIN (real mass + friction),
  integrated with g-JF (Gronbech-Jensen & Farago, Mol.Phys. 2013) or BAOAB (Leimkuhler & Matthews
  2013), ONE step per outer dt; friction gamma sets D = kT/(m gamma) in the diffusive limit,
  step-count-independent. g-JF gives provably correct flat-potential D and exact harmonic Boltzmann
  for ANY gamma/dt within the Verlet stability limit; BAOAB adds high-friction/diffusive-limit accuracy.
- EXISTENCE PROOF: Dry Martini (Arnarez ... Marrink, JCTC 2015, 10.1021/ct500477k) runs exactly this
  -- GROMACS `sd` (leap-frog inertial Langevin) at ~20 fs, friction (tau_t ~ few ps) tuned to recover
  lipid diffusion (D-vs-friction table in its SI). Hard-core (LJ+Coulomb r^-12) implicit-solvent
  lipids, single ~20 fs step, no sub-stepping.
- MECHANISM (why inertial beats overdamped on hard cores): inertia velocity-limits the per-step
  displacement (~v_thermal*dt) and the core decelerates the bead before overlap; the overdamped
  displacement (dt/gamma)|grad U| instead diverges at contact. The single step is bounded only by the
  Verlet stability limit of the stiffest core (tens of fs for CG beads), NOT by needing inner substeps.
- RULED OUT: MALA / Metropolis-adjusted BD (stable but accept/reject destroys the dynamics -> wrong
  timescale); mollified / Langevin-impulse MTS (still sub-steps the fast force); DPD (large dt only
  because cores are SOFT -- violates our hard-core rule).

RECONCILIATION with Update 35 ("single-step space exhaustively closed"): that closure was for the
OVERDAMPED family (implicit, coupled, GLE, reflecting BD). The INERTIAL-Langevin route is different
and was NOT among them. It is corroborated in-house: our PROTEIN already runs inertial Verlet at
dt=0.009 stably against the same martini_potential hard cores (incl BB-env) that make the OVERDAMPED
lipids cap-fire -> the instability is a property of overdamped integration, not the cores. So making
the LIPIDS inertial Langevin (+ friction for D) should be single-step-stable at dt=0.009 too.

IMPLICATION: prototype inertial-Langevin lipids (g-JF or BAOAB, real CG mass, friction calibrated to
D_DOPC with the factor-4 + Arrhenius law) as a single-step replacement for overdamped-BD + M=90
sub-stepping. If D + Arrhenius + bilayer stability hold, this removes the M-loop (~up to 90x fewer
inner computes) AND removes the protein-interface cap-firing that forces M=90 (Update 38). Major
integrator change -> needs go-ahead + careful D/stability validation. Papers: dry Martini
10.1021/ct500477k; g-JF arXiv:1212.1244; BAOAB arXiv:1203.5428.

## Update 40 (2026-07-18): g-JF inertial-Langevin single-step -- STABLE, but D ~10x too slow (underdamped)

Prototyped g-JF (Gronbech-Jensen & Farago) inertial Langevin for the lipids, gated by /input/brownian
attr "langevin" (default 0 = overdamped, master parity). Bilayer, single step dt=0.009 (no RESPA):
- WORKS mechanically: STABLE (rc=0, potential -11272 cohesive, finite), correct thermostat
  (Tkin = 1.00 kT), bilayer structure preserved (frac 1.00, thickness 30.4 ~ overdamped 29.9).
  Single-step against the hard core IS achievable (as the literature predicted; dt=0.009 ~ 1.8 fs).
- BUT diffusion is ~8-10x too SLOW at the reused (overdamped-calibrated) friction. Matched-time PO4
  MSD (A^2): t=2.7/10.8/27 Tup -> g-JF 0.73/3.56/7.78 vs overdamped 24.4/36.2/60.1. g-JF over 360 Tup
  (26.5) is still below overdamped at 27 Tup (60) -> it is NOT catching up; effective D_gjf is genuinely lower.
- CAUSE: real CG mass (~6 m_up) + calibrated zeta -> tau=m/zeta ~ 57 ps -> the inertial lipid is
  UNDERDAMPED (rattles ballistically in its cage; momentum ramps from rest). g-JF gives the correct
  FREE-particle D=kT/zeta + correct Boltzmann (literature), but the effective diffusion in a CAGED
  bilayer is set by Kramers cage-escape, which in the underdamped regime is slower than the overdamped
  Smoluchowski dynamics at the same zeta. The overdamped ZETA calibration does NOT transfer.

VERDICT: single-step stability is real (mechanism validated), but correct TIMESCALE is not free --
g-JF needs its OWN friction re-calibration in the inertial/underdamped regime (as dry MARTINI did) +
the Arrhenius D(T) re-fit, complicated by the Kramers turnover (lowering zeta raises free D but can
slow cage escape further). Not a drop-in. Trade-off: overdamped BD = correct D, needs M=90 sub-steps;
g-JF single-step = fast + stable, D-wrong until re-calibrated. PRODUCTION STAYS on overdamped M=90.
g-JF kept as a GATED prototype (langevin attr, default off) for a future re-calibration effort.

## Update 41 (2026-07-18): g-JF calibration EXHAUSTED -- cannot reach correct D single-step

Followed Update 40 by trying to calibrate g-JF to the overdamped D (bilayer, single-step dt=0.009;
overdamped ref MSD@27=60):
- FRICTION scan (zeta x0.3..x30, real mass): MSD@27 stuck ~6.5 A^2 (cage-trapped), INSENSITIVE to zeta
  over 100x, ~10x below overdamped. Kramers turnover: low zeta -> underdamped/slow cage escape; high
  zeta -> low free D. Friction is NOT the knob at real mass.
- Terminal D (g-JF x1, 1800 Tup run): D=0.0055 -- ~50x below the overdamped bilayer-apparent (0.28);
  severely cage-trapped even at long times.
- MASS scan (real -> x0.0001, zeta x1): mass IS a knob (lighter -> more overdamped -> faster). MSD@27:
  x1=7, x0.1=14, x0.01=25, x0.001=41 (rising toward 60) -- BUT the single-step STABILITY WALL hits at
  x0.0003 (blow-up/nan: the overdamped-limit core overshoot returns as inertia vanishes). Lightest
  stable mass x0.001 -> ~41 (still ~1.5x low).
- MASS x0.001 + FRICTION tweak: MSD@27 peaks ~48 at zeta x0.5, then declines; never reaches 60.

VERDICT: g-JF single-step gives correct STABILITY + thermostat (Tkin=1.00) + structure (mechanism
works, as the literature claims), but CANNOT reproduce the calibrated lipid D. Correct D needs the
overdamped limit (light mass), which reintroduces the exact core-overshoot instability the
sub-stepping was built to avoid -> the correct-D-vs-single-step-stability tension is FUNDAMENTAL and
persists in the inertial regime. Best achievable (mass x0.001 at the stability edge + zeta tune) is
~0.8x the target D, with an unphysical 1000x mass + fragile stability + near-certain hybrid failure
(the protein interface needed M=90 even for overdamped). dry MARTINI avoids this only by accepting an
approximate D at a larger timestep. CONCLUSION: overdamped M=90 remains the correct-D production path.
g-JF retained as a gated prototype (langevin attr, off) documenting the exhausted route.

## Update 42 (2026-07-18): LITERATURE (round 2) -- correct-D + ~40ps single step on hard cores is fundamentally impossible

Second review targeting the SPECIFIC tension (correct D + single ~40 ps step + hard r^-12 core, no
sub-step, no softening, no surrogate). VERDICT: no published, transferable, spline-compatible method
achieves it. Correct caged long-time D REQUIRES resolving excluded-volume collisions at their natural
scale; you can change the FORM of that resolution but not remove it. Three equivalent forms (all same
collision-resolution cost, none a genuine big step):
- RESPA SUB-STEPPING (our current M=90): exact D; uses a displacement-cap safeguard.
- DYNAMIC MONTE CARLO (DMC), time-rescaled (Sanz-Marenduzzo JCP 2010; Cuetos-Patti PRE 2015;
  Corbella-Cuetos-Patti PCCP 2018): small trial moves + REJECT hard-core overlaps + rescale time by
  the acceptance ratio -> provably = overdamped BD in the small-move limit -> CORRECT D, UNCONDITIONALLY
  STABLE (rejection can never NaN), transferable, needs only U (spline-compatible, no analytic events).
  NOT a big single step -- it re-expresses M sub-steps as M stable rejection micro-moves. Wins over our
  RESPA+cap: guaranteed correct D + unconditional stability + NO arbitrary cap (rejection = exact
  excluded volume). Same cost.
- EVENT-DRIVEN / first-passage BD (Scala-Voigtmann JCP 2007; Donev FPKMC 2010; GFRD): exact D via
  analytic collision times -- needs ANALYTIC event times (incompatible with a spline table; our naive
  reflection NaN'd for exactly this) and DEGENERATES in a dense caged bilayer. Ruled out for us.

ONLY ROUTE TO A GENUINE ~40 ps SINGLE STEP: accept an APPROXIMATE D + TIME-RESCALE (the Martini
standard). Martini does not reproduce absolute D; it accepts it and multiplies reported time by a
calibrated factor (~4, but "3-8x, interpret with caution", process-dependent). So g-JF single-step
(stable, fast) + a calibrated time-rescale = the 40 ps/step goal, at the cost of approximate D. CAVEAT
(critical): the rescale factor is process- AND temperature-dependent -> a single global factor can
distort our Arrhenius D(T); must calibrate across the T-range and verify ~constant before relying on it.

BOTTOM LINE: exact-D + 40ps-single-step + hard-core is IMPOSSIBLE (fundamental). Choose (a) EXACT D via
collision resolution -- current M=90 sub-stepping, or DMC (cleaner: correct D, unconditionally stable,
no cap, same cost); or (b) 40 ps SINGLE STEP via g-JF + Martini-style time-rescaling (approximate D;
verify T-constant rescale). No option is simultaneously big-single-step + exact-D + hard-core.

ADDENDUM (2nd targeted review, 3 categories): the ARITHMETIC makes it concrete -- at DOPC D~1e-5 nm^2/ps
the net lateral move in 40 ps is only RMS ~0.4 A, but the bead undergoes ~1e3-1e4 sub-ps cage collisions
in that window; correct D is the Kramers balance of those collisions, so a 40 ps step must INTEGRATE the
collisions out, not resolve them -> only a coarse-time propagator (learned or externally-supplied D) can
do that. Also explicitly ruled out: GJ-2GJ / correct-kinetics Verlet thermostats (linear/equilibrium
guarantees only; Verlet stability limit unchanged), colored-noise GLE (sampling tool, same small step),
super-time-stepping RKC/RKL (multi-stage sub-stepping in disguise; linear-parabolic only), Mori-Zwanzig
memory kernels (coarse-grain DOFs not the time axis; sub-ps kernel -> fine step). The ONLY literal
"one step = long physical time" methods are ML transfer-operator propagators (ITO/Timewarp/TITO 2023-26)
-- forbidden here (non-transferable surrogate). FINAL: must relax exactly ONE of {single-40ps-step,
exact-D, hard-core, universal, spline}: drop the big step -> DMC/sub-step (exact D); or drop exact D ->
g-JF + calibrated effective-time rescale (big step). Design space closed (2 reviews + full empirical test).

## Update 43 (2026-07-18): Route C (g-JF + time-rescale) -- T-viable for the bilayer, but the rescale is inconsistent for the coupled hybrid

Tested Route C (g-JF single-step + calibrated effective-time rescale). Multi-T (280/303/333 K)
g-JF/overdamped MSD@27 ratio: R(T) = 7.36, 6.29, 8.23 (mean 7.3, spread 1.31x, non-monotonic ->
consistent with the ~20-30% run-to-run D noise + within Martini's own "3-8x, caution" tolerance). So a
single rescale R~7.3 preserves the Arrhenius D(T) at bilayer precision -> Route C is T-VIABLE for a
PURE BILAYER (correct effective D via a ~7.3x time-rescale).
Speedup: only ~6x -- single-step g-JF is ~45x fewer force evals/step than M=90, but the 7.3x-slower
per-step D needs ~7.3x more steps to sample the same diffusion -> net ~90/(7.3*2) ~ 6x.
HYBRID OBSTACLE (decisive): the rescale changes the LIPID effective clock by 7.3x relative to the
Verlet protein. In the coupled hybrid the g-JF lipids diffuse ~7.3x slower (in sim-time) than the
protein moves; a GLOBAL time-rescale then makes the protein 7.3x too slow (wrong for the trained
Upside FF), while a lipid-ONLY rescale corrupts the in-sim protein-lipid coupling (the rescale is a
post-hoc time relabel; it cannot fix the coupling that happens at the sim step). So a single rescale
cannot put protein + lipids on one correct clock -> Route C does NOT cleanly serve the membrane-protein
hybrid (the project goal). Route B (DMC) avoids this entirely: it gives the EXACT lipid D with no
rescale, preserving the current correct protein-lipid relative timescale (same as overdamped M=90),
while dropping the displacement cap.

## Update 44 (2026-07-18): all-atom 1RKL tilt reference (100 ns GPU) -- dry-MARTINI OVER-TILTS

All-atom 1RKL (CHARMM-GUI 8430021493) ran 100 ns on midway3 GPU (V100, 318 ns/day). Backbone tilt
(31 CA, PCA principal axis vs membrane normal z) relaxes from 34.7 deg (insertion) to ~16-22 deg
(2nd-half mean 22.5+/-6.2, last 20 ns 16.3+/-3.2) -- near-vertical/moderate, consistent with
wet-MARTINI ~10 deg + literature 5-25 deg. So the CORRECT 1RKL tilt is ~16-22 deg; the dry-MARTINI
hybrid's ~40 deg is an OVER-TILT (~15-25 deg too high). The all-atom relaxes 34->~20 in ~50 ns; the
dry-MARTINI hybrid did NOT relax similarly. Whether the dry over-tilt is (a) incomplete relaxation
(dry runs too short) or (b) a genuine dry-MARTINI FF over-tilt needs a longer dry 1RKL run to
distinguish. The all-atom ~16-22 deg is the ground truth to match. protein.xtc/gro saved in scratchpad.

METHOD-CONSISTENT 1RKL tilt comparison (all BB/CA-PCA principal-axis vs z):
  wet MARTINI 2.2 (explicit water, 30 ns): 10.3 +/- 5.6 deg  (relaxes 29.7 -> ~10)
  all-atom (CHARMM-GUI, 100 ns GPU):       ~16-22 deg (2nd-half 22.5, last 20ns 16.3, still dropping)
  literature (matched TM helix, DOPC):     5-25 deg
  dry-MARTINI hybrid (Upside):             ~40 deg  <-- OUTLIER / over-tilt (~20-30 deg too high)
Three independent refs agree 1RKL is near-vertical/moderate; the dry-MARTINI hybrid over-tilts. This
dovetails with the Route-C obstacle: too-slow lipid rearrangement around the protein -> the tilt DOF
cannot relax -> it stays stuck high. So the over-tilt may be a lipid-dynamics/coupling artifact
(protein-adjacent lipids not rearranging fast enough) rather than an FF error -- a longer dry run
and/or the interface-lipid dynamics are the next thing to check.

## Update 45 (2026-07-18): LITERATURE (round 3) -- global time-rescale is INVALID for the coupled hybrid; the real question is over-tilt KINETIC vs THERMODYNAMIC

Third targeted review (Route C for a COUPLED protein-lipid hybrid). VERDICT: a single global
effective-time rescale is PROVABLY INVALID for a multi-component system -- CG acceleration is
component/process-dependent (Martini factor 2-10x, "not uniform for all DOF"; "may be different in
different components"; "a single reduced timescale cannot describe both components" -- Marrink CSR 2013,
WIREs 2023; friction-parameterization JCP 2023; time-scale-consistent blends J.Phys.Mater. 2020;
protein over-aggregation Javanainen PLoS One 2017 / Alessandri JCTC 2021). A scalar rescale relabels
the shared clock but cannot repair a RELATIVE mismatch. The sanctioned fix is PER-COMPONENT
friction/memory matching (equalize the clocks BEFORE any rescale), not a global scalar. For us:
- Can't friction-match the lipids UP (7.3x slowdown is a hard-core Kramers cap -- dead knob).
- Only self-consistent Route-C salvage: add friction to the PROTEIN integrator to slow it to the
  lipids' 7.3x clock, then rescale x7.3. Legit in principle but (i) discards protein-integrator
  efficiency, (ii) needs protein friction tuned to exactly the lipid factor, (iii) works ONLY IF the
  mismatch is a clean uniform KINETIC factor.
- No root-fix integrator (correct overdamped-D + single-step + hard-core) exists (rounds 1-2).
=> Route C is inherently pure-bilayer-only; for the hybrid it distorts relative protein-lipid dynamics.
Salvages: (1) friction-match protein->lipid clock then rescale (sub-project, conditional); (2)
source-level lipid-D correction (MZ/GLE kernel or DMC/event-driven) = reintroduces sub-stepping;
(3) accept the artifact (bilayer-only). None a clean win.

CRITICAL PIVOT: the ~40 vs ~20 deg over-tilt is either KINETIC (lipids too slow to let the tilt relax
-> correct-D + longer run fixes it) OR THERMODYNAMIC (dry-MARTINI FF prefers ~40 -> NO integrator/
timescale fix helps; FF limitation). Must be determined BEFORE more integrator work. TEST: run
dry-MARTINI 1RKL (M=90, correct D) from BOTH a tilted (~40) and a vertical (~20) start; converge to the
same tilt = FF equilibrium (thermodynamic); each stays near its start = kinetic trapping. This decides
whether the over-tilt is a sampling problem (fixable) or an FF problem (fundamental) -- more decisive
than the C/B integrator choice.

MEASURED dry-MARTINI tilt (CA-PCA, same method) -- FROZEN DOF, not an equilibrium: flex_raw_009_ext
40.7+/-0.4, flex_raw_009_long 40.8+/-0.2, o3_long 45.6+/-0.8. The tilt sits at its INITIAL ~40 deg and
never moves; std 0.2-0.4 deg vs a real fluctuating tilt's +/-5-6 deg (all-atom, wet MARTINI). So the dry
tilt DOF is FROZEN (kinetic trap), NOT a thermodynamic equilibrium at 40 deg. The protein was inserted
at ~40 deg and cannot relax -- consistent with the interface-stiffness finding (annular lipids cap-fire
~1.7% even at M=90 -> they do not rearrange -> the helix cannot straighten). (The earlier "dry -> vertical"
was on dispersing/frozen confounded runs, retracted Update 29.)

WET-MARTINI EVIDENCE (points to KINETIC + rules out "MARTINI itself"): wet MARTINI 2.2 1RKL RELAXES the
tilt 29.7 -> ~10 deg over 30 ns (measured, BB-PCA), i.e. the standard MARTINI FF gets the tilt right and
the tilt DOF CAN relax. So the over-tilt is NOT a generic MARTINI problem -- it is specific to the
dry-MARTINI-in-Upside stack (dry lipid FF, Upside protein FF, integrator, implicit solvent). The
wet-relaxation strongly suggests the dry ~40 deg is a KINETIC failure-to-relax (slow dry-lipid dynamics
around the protein) rather than an FF equilibrium at 40 deg -- to be confirmed by the tilted-vs-vertical
convergence test above.

## Update 46 (2026-07-18): effective-time-factor DIRECTION -- code comment is inverted

Working out the manuscript physics exposed a sign error in the effective-time-factor wiring.
g-JF lipids diffuse ~R=7.3x TOO SLOW (D_exp = R * D_sim, measured R = od_MSD/gjf_MSD ~ 7.3 > 1).
The correct mapping is therefore physical_time = sim_time / R and D_physical = R * D_sim (time is
COMPRESSED, opposite sign to MARTINI's factor-4 which dilates time because MARTINI over-diffuses).
But run_sim_hybrid.sh + martini_prepare_system_lib.py document it as "physical_time = sim_time * this"
with UPSIDE_LIPID_TIME_FACTOR=7.3 -- that reads as time DILATION (as if the lipids were too fast).
CLEANUP FIX: keep the stored value = R (7.3, the diffusion multiplier D_phys = R*D_sim) but correct the
comment/semantics to "physical_time = sim_time / effective_time_factor" everywhere it appears (run
script, prep lib, any analysis that consumes the attribute). Verify no analysis multiplies time by R.

## Update 47 (2026-07-18): refined R calibration + temperature scope + R(t) lag-dependence

Final calibration (calib_final.py, clean windowed MSD, PO4 single-bead, COM-removed) on the clean binary:
- R@27 (matched diffusive-onset window): 280K=7.27, 303K=6.70, 333K=8.01; spread 1.20x.
- USER SCOPE: only protein-relevant temperatures matter -- room temp up to ~50-60 C, i.e. ~0.85-0.95 T_up
  ("around 1.0 T_up"). 280K (7 C) is below that -> DROPPED. Relevant range = 303K (room) + 333K (60 C),
  R ~ 6.7-8.0, i.e. R ~ 7. This confirms the documented ~7.3 for the range we care about.
- R(t) LAG-DEPENDENCE (at 303K): R = od_MSD/gjf_MSD declines monotonically 53.8(t=1) -> 6.7(t=27) ->
  5.1(t=81 T_up), NO clean plateau in the accessible window -- the overdamped reference cages (MSD
  saturates ~62) while g-JF creeps, so the ratio keeps shrinking. => a single scalar R is a
  timescale-anchored APPROXIMATION (matches the diffusion at the lag it is calibrated to), reinforcing
  the "approximate for coupled kinetics" caveat. Reported R is at the diffusive-onset window (~7).
- NOTE: t/T_up in the R(t) table is TIME (Upside time unit ~0.2 ps), NOT temperature -- do not conflate.
- Tex updated: dropped 280K, table = 303K(6.7)/333K(8.0), R~7 over room->60C, R(t) note added.

## Lessons

- Before building a corrective FF term to fix a symptom, rule out SETUP artifacts first --
  timestep and sub-step resolution above all. The BB-env PMF (O2) was built to fix a protein
  "kick" that was actually caused by a non-standard timestep (0.004 vs Upside standard 0.009,
  inherited from the abandoned CGL) + under-resolved lipids (M=10 -> cap-firing -> energy
  injection). At the standard, well-resolved setup the raw FF needed no PMF -- and the PMF
  itself then CAUSED the drift. Check that dt matches the standard examples (04.HDX) and that
  sub-stepping is converged BEFORE adding any corrective potential/softening.
- Validate tilt claims with the actual conservative energy landscape, not just
  the trajectory. The force field's static landscape was fine; the pathology was
  the coupled membrane relaxation (ratchet) enabled by too-cheap tail collapse.
- A "relief" derived from an isolated molecule can be non-physical in a packed
  phase: it over-rewards deformation because the neighbor/packing cost is
  omitted. Derive deformation free energies in the phase where they act.
- Prior `.bak` names (`target_angle_massscale`) show earlier attempts twisted
  parameters / mass to hide the angle drift. Do not repeat: the fix must restore
  the missing packing physics, not damp the dynamics.
- Removing whole Python functions programmatically: a naive "def line -> next
  column-0 line" scan MIS-cuts multi-line signatures whose closing `)` sits at
  column 0 (e.g. `def f(\n    a,\n) -> T:`), deleting only the header and orphaning
  the body. Use `ast` (node.lineno..node.end_lineno) for reliable spans, or verify
  with `py_compile` immediately after each bulk removal. (Hit this on
  prepare_system_lib.py / prepare_system.py during the CGL cleanup.)

## Cell-list neighbor search — performance lessons (2026-07-18)

Replaced the martini_potential O(N^2) rebuild scan with an O(N) cell list. Key facts
learned (rigorous measurement on the M=90 full-res DOPC bilayer):

- The cache-friendly O(N^2) scan is HARD to beat at small N. A naive linked-list cell
  list REGRESSED (277k vs scan 166-218k us/step at 1792 atoms) because linked-list
  traversal + random pos reads + per-candidate binary search are cache-hostile, while the
  scan is a tight sequential/vectorizable loop over contiguous pairs[]. Fix = counting-sort
  atoms into cells and read positions from a cell-sorted scratch array (sequential), and
  look up the pair index only AFTER the distance filter (few in-cutoff candidates, not all).
  This cut the 1792-atom rebuild 12.6ms -> 4.2ms, making cb2 match the scan; the O(N) win
  then shows at scale (7168 atoms: 1.7x/step vs scan; old cb10 unrunnable there).
- The RESPA step at the M=90 test size is dominated by the Brownian integrator + the fixed
  spline compute over in-cutoff pairs, NOT the rebuild. So the rebuild optimization alone
  gives little at 1792 atoms; the 1.28x there comes from lowering cache_buffer 10->2 (leaner
  active list) which the cheap O(N) rebuild makes affordable. "Several-fold" is a large-N
  property (O(N^2) rebuild), not a 1792-atom one.
- cache_buffer is a U-curve: too small (cb1) over-rebuilds; too large (cb10) fattens the
  90x/step compute list. With an O(N) rebuild the optimum is small (cb2-4). Kept 2.0.
- Bit-identity discipline: sort the active list ascending (== old scan np order) so
  energy/force accumulation order is preserved; compute dr with the SAME operand order
  (pos_a - pos_b, a<b) and minimum_image so the inclusion test is bit-identical. A skin
  pair (12-14 A) contributes zero, so a trajectory diff will NOT catch a skin set mismatch —
  gate the set equality directly with an in-code assert (old scan vs cell), then remove it.

## Update 48 (2026-07-18): unit conversion baked into h5 (P11) -- coulomb_k is dead in compute

Moved native->Upside unit conversion out of C++ into the Python h5-build. Key realizations:
- The main martini_potential node's `coulomb_k` is READ but NEVER used in compute -- the LJ+Coulomb
  energy is fully baked into `combined_energy_grids` by convert_stage (particles group already in eup:
  unique_eps_eup/unique_sig_ang/combined_energy_grids). So the 3 old conversion attrs on the config
  node were consumed only by the PYTHON softening builder (post-process), not the engine. That is why
  Gate C (bilayer) is BIT-IDENTICAL (-11039.306641, |diff|=0.000000): changing where coulomb_k is
  computed cannot move the potential. New scheme writes a single baked attr `coulomb_constant`.
- sc_table conversion: keep the whole build/read logic in native units and convert+rename ONLY at the
  final create_dataset step (grid_ang=grid_nm*10, *_energy_eup=*_kj_mol/2.914952774272). The C++ tail-
  subtraction shift is retained (it is a physics zeroing at the cutoff, applied identically before/after);
  (native[ig]-tail)/E == eup[ig]-eup_tail, so the float32 result is unchanged (bake was exact, maxabs 0.0).
- Dimensionless datasets (angular_profile, rotamer_angular_profile, cos_theta_grid) are NOT converted.
- The engine is now unit-agnostic: it reads grid_ang/*_energy_eup/coulomb_constant and does zero math.

## Update 49 (2026-07-19): hybrid production NaN root-caused = MIXED-INTEGRATOR coupling; fixed

Both 1rkl + 1afo NaN at stage-7 production (g-JF single-step lipids + Verlet protein). Faithful reproducer
(1afo.stage_7.0.up frame0 = -17409) NaNs at step 3500, matching the log. Diagnosis (each ruled out with a
test): NOT friction (alpha*dt/2m ~ 1e-6 at CG mass -> g-JF friction negligible; x1..x30 barely moved kT
1.96->1.79; higher -> overdamped-unstable); NOT the compute/force order (reorder so protein uses f_n:
kT 1.96->1.94, no effect); NOT a missing push-off (BB-env has 876960 lipid-protein pairs in the lipid
force; SC-env sc_table_1body writes grad to the env/lipid atoms); NOT P11 (heating predates it, go/no-go
was kT 1.06 at 1000 steps). ROOT CAUSE: the g-JF lipid uses the symmetric (a*f_n+f_n1)*dt/2 velocity
update; the Upside protein integrator does not -> they exchange the shared interface force inconsistently
-> energy injected at the hard interface -> hot lipids overcome the push-off, creep 4A->2A into the
protein (verified in the trajectory) -> r^-12 spike -> NaN. The pure bilayer is fine (single consistent
integrator). DISCRIMINATOR: all-atom g-JF (consistent) -> kT 1.96->1.10, 20000 steps no NaN, protein
folded (Rg 15.8->15.4).
FIX (prep-side, martini_prepare_system_lib.py, gated langevin=1): put ALL atoms on the g-JF path;
friction_scale=0 for PROTEIN (g-JF at zero friction == velocity-Verlet, thermostatted by the Upside OU),
=1 for lipids (g-JF unchanged -> lipid D/timescale preserved). Overdamped fallback (langevin=0) unchanged
(M-substepping resolves the interface). No C++ change. Verified rebuilt-binary: hybrid 5000 steps no NaN
kT 1.10. PROTEIN-ALONE check (chig): velocity-Verlet(fric=0)+OU vs Upside-default v = energy 17.0 vs 17.4,
kT 0.984 vs 0.986, Rg within 1 sigma -> protein core ensemble unchanged.
OPEN: DOPC lateral D -> ps/step. g-JF is sub-diffusive (cage-trapped) so the effective-time factor is
lag-dependent (R does not plateau, cf Update 47); measured ~9-15 ps/step (lag 27-90) vs the 40 ps/step
target. Awaiting user's convention (which D anchor; is total = 4 x effective_time_factor) before tuning
the factor. The fix does NOT change the lipid D.

## Update 50 (2026-07-19): C++ analytic LJ/Coulomb encoding removed (spline-only engine)

martini_potential eval is purely spline (combined_spline). Removed the dead analytic legacy: node-level
epsilon/sigma + coulomb_k reads (C++); the eps/sig/qi/qj/qq analytic coefficients in PairParam (the spline
is now indexed one-per-unique-coefficient-row, active-pair via coeff_is_zero[ip]); unified lj_cutoff/
coul_cutoff -> single cutoff attr. KEPT the coulomb_constant attr: the PYTHON soft-core builder (coulomb_
soften=1 in the minimization/npt_equil stages, which ARE live) reads it; the C++ no longer does. Bit-
identical (-11039.306641, trajectory maxabs 0.0). Engine is now unit-agnostic AND analytic-LJ/Coulomb-free.

## Update 51 (2026-07-19): friction_scale clamp defeated protein velocity-Verlet; fixed

Bug in martini_brownian.cpp load: `friction_scale[i] = (v > 0.f) ? v : 1.f` silently rewrote the protein's
intended 0 -> 1, so the "velocity-Verlet + OU" protein (Update 49) actually ran with full g-JF friction.
Confirmed the g-JF math collapses to EXACT velocity-Verlet at friction_scale=0 (b=1,a=1,bnoise=0; lines
137/145/153-162). It stayed stable only because gamma is tiny (~0.0035) so the spurious friction was
negligible -- but it was not the clean chosen integrator, and the load comment was stale (old "protein
uses larger friction" design). Fix: clamp `>= 0.f` (0 legit for g-JF; negatives still fall back to 1) +
rewrote the comment. langevin=0 (RESPA) never writes friction_scale so its `mob_base/sc` div-by-zero path
is unreachable. Needs rebuild + re-verify.

## Update 52 (2026-07-19): DOPC diffusion -- the 40fs/step + factor-4 + go-to-lag picture

Resolved the recurring "g-JF too fast vs too slow" flip: it was a time-UNIT confusion. The raw g-JF motion
has a fixed MSD/step; the D you report is whatever the declared time/step makes it. Analytic anchor (bilayer
eq_m90, g-JF, formal T_up=0.2029ps, dt=0.009 -> 1.826 fs/step):
  - At the FORMAL 1.826 fs/step: D_apparent ~3400 um2/s at lag 27 (÷4 -> 850), ~74x too fast vs DOPC 11.5.
  - Declaring each step = 40 fs (dry-MARTINI standard, R_eff=40/1.826=21.9) THEN factor-4 (total 160 fs/step):
    short lag (27) reads ~39 um2/s (too fast) but at lag ~270 D settles to ~11.3 um2/s ~= exp 11.5.
This is EXACTLY the user's "use 4, test D first (short lag = not working), then go to lag": the short-lag D
is ballistic/cage-rattle (too fast); the long-lag D (diffusive regime) matches DOPC at 40 fs/step + factor-4.
NOTE the target is almost certainly 40 FS/step (dry-MARTINI standard), not 40 ps (which would make g-JF
~250x too slow -- physically impossible with this displacement). CAVEAT: g-JF is sub-diffusive so verify a
genuine diffusive plateau exists near lag 270 (long-run in progress) before fixing effective_time_factor.

## Update 53 (2026-07-19): g-JF bilayer has NO diffusive plateau -- the lag-270 "match" is a transient crossing

Long g-JF bilayer run (eq_m90, 400000 single-steps, 3600 T_up, 2001 frames) MSD(lag):
  lag  27/54/108/216/432/864/1728 -> MSD 7.5/10.6/14.4/19.6/26.1/33.5/40.3 A^2
  local alpha (d lnMSD/d lnlag): 0.51 -> 0.44 -> 0.44 -> 0.42 -> 0.36 -> 0.27.
alpha FALLS toward 0 (caging), NOT toward 1 (Fickian). So the g-JF bilayer never reaches a diffusive
regime on this window; the earlier "D ~ 11.3 um2/s at lag ~270" (Update 52) is a COINCIDENTAL crossing on
a falling curve, not a real diffusion coefficient. Consequence: no single effective-time factor / factor-4
maps the g-JF to one DOPC D -- confirms + sharpens the lag-dependent-R problem (Updates 47, 52). Still
MOVING (MSD 7.5->40, not flat) so not a hard cage, but strongly sub-diffusive.
RESOLVED (overdamped M=45 control, same box, lag 27->54): overdamped D=0.425->0.271, alpha=0.35 -- ALSO
sub-diffusive, and NOT better than g-JF (alpha 0.51) at these lags (overdamped D is ~6x g-JF, consistent
with R~6-7). So the sub-diffusivity is NOT unique to the g-JF near-zero friction -- it is present with proper
friction too => a SYSTEM/regime property (small 128-lipid patch + the physical short-intermediate-time
sub-diffusive regime, real DOPC alpha~0.5-0.6). CONSEQUENCE: neither integrator yields a clean Fickian D on
this patch; a diffusion-based timescale would need a much larger membrane + far longer runs (to reach the
~10-100 ns Fickian crossover) for EITHER integrator. The "use overdamped for D" fallback is thus also
system-limited here, not a quick win.

## Update 54 (2026-07-19): REMD works correctly for lipids (subagent-verified); one bug fixed

Upside REMD is C++ single-process multi-System, COORDINATE-swap convention (temperatures pinned to slots;
only positions + replica_indices swap; ReplicaExchange in main.cpp). Verified correct for dry-MARTINI:
(1) per-slot T reaches the lipid integrator -- set_brownian_temperature(sys.temperature) is called at the
thermostat cadence (main.cpp ~1237), same as the protein OU; slot T is fixed so lipid noise/gamma always
use the right kT. (2) Lipid potential IS in the swap criterion (MartiniPotential is a potential_term; the
swap uses engine.potential = sum of all potential_term nodes). (3) NO momentum rescaling needed/appropriate
-- config-swap keeps momenta+T pinned to slots (Maxwell-Boltzmann at each slot's T); kinetic terms cancel
in the acceptance ratio. (4) Arrhenius gamma(T) only sets timescale, not the equilibrium exp(-U/kT), so it
is REMD-transferable. FIXED: main.cpp set_brownian_temperature was called at setup BEFORE sys->temperature
was initialized (uninitialized read; register already defaults rt->temperature=t_ref, and the run loop
re-sets it at nr=0, so it was harmless) -> moved the call to right after sys->set_temperature(). Rebuilt,
bilayer bit-identical (-11039.306641), hybrid stable. LATENT (not the active path, flagged not fixed): the
`mv` RESPA integrator overload does not call apply_brownian_step nor skip brownian_mask -> mv + brownian
lipids would silently mis-integrate the lipids. MARTINI uses `v`, so not triggered.

## Update 55 (2026-07-19): metadynamics feasibility (subagent-evaluated)

Deriv engine is a reverse-mode autodiff graph (DerivComputation nodes: CoordNode has output+sens;
PotentialNode writes grad into its input's sens). A metad bias is just a new PotentialNode reading a CV
CoordNode, adding a sum-of-Gaussians bias energy + writing dV/ds into cv.sens -> engine propagates to atoms.
Template: src/spring.cpp (Spring/WallSpring; commented AFM/MovingSpring shows the time-dependent pattern).
CORE (protein): VERY FEASIBLE ~1-2 days. Most CVs already exist as nodes (Distance/Angle/Dihedral/
GroupCenter/Sigmoid/Sum -> composite CVs incl native-contact Q via Distance->Sigmoid->Sum). Rg-as-CV, RMSD
need small new CoordNodes (Rg exists only as a print-only diagnostic). HYBRID (membrane): node is identical
(integrator-agnostic), but (a) membrane geometry CVs (thickness/APL/leaflet/tilt/insertion) mostly DON'T
exist -> new geometry nodes (tilt/insertion trivial via GroupCenter+VectorAngle/Distance1D; thickness/APL
+2-4 days each); (b) sub-diffusive lipids -> metad on collective LIPID CVs is fragile (deposition must be
slow vs the slow CV relaxation); protein-CV metad-in-membrane is well-posed (protein is velocity-Verlet+OU,
normal timescales). PLUMED: write_plumed (advanced_config.py) is a DEAD STUB -- no C++ node, not in build,
no libplumed anywhere. Native (well-tempered, bias on a grid evaluated by the existing deBoor spline evals)
fits the autodiff-graph + spline-table philosophy far better -> recommend native. Deposition timing must
key off the true round_num via a set_round hook (mirroring set_brownian_temperature) because compute(Deriv
Mode) runs ~4x/round (3 stages + the g-JF recompute). CONFIRMED: effective_time_factor is NEVER read by the
engine (register_brownian_for_engine ignores it) -- it is analysis-only metadata, by design (C++ = pure sim
time); so deposit Gaussians on the sim-time (round) schedule, convert to physical time only at analysis.

## Update 56 (2026-07-19): what CAN be claimed at the correct timescale -- static mechanics, not transport

Question: if lateral D can't anchor the bilayer timescale (sub-diffusive, findings 53), what can? Principle:
diffusion (and full reorientation, flip-flop) are CAGE-ESCAPE processes needing PHYSICAL friction; the g-JF
has ~none -> ALL long-time transport is wrong. Only LOCAL/inertial motions (set by physical CG mass + correct
FF curvature, friction-independent) are ~right. So no long-time dynamical observable anchors the g-JF time.
STATIC STRUCTURE measured (g-JF, 128-DOPC, T=0.8647, gjf_long.up 1500 frames) -- validates the FORCE FIELD:
  APL = 65.6 A^2 (exp ~67.4, -3%)  |  P-P thickness = 43.3 A (exp 37-39, +10-15% THICK)
  chain order P2 tail-avg = 0.30 (MARTINI DOPC 0.2-0.3, high end)  |  tilt 22.5+/-12.3 deg, director S=0.74
  -> bilayer is FLUID (not gel; S=0.74/tilt 22 is fluid-range) but somewhat OVER-ORDERED/THICK. FF-fidelity
     note (report, do NOT twist per CLAUDE.md). Excess thickness partly the area-thickness anticorrelation at
     the fixed (slightly-small) NVT box; a zero-tension NPgT equilibration would let the FF pick its APL/thick.
ROTATIONAL ACF (director PO4->mean(C5A,C5B), C2(t)=<P2(u0.ut)>): 1.00->0.997(1.8)->0.94(108)->0.80(1728),
  smooth monotonic, NO plateau -- decays SLOWLY = same cage-escape wall as lateral D (confirms the prediction:
  reorientation is also not a clean g-JF timescale anchor). BUT it is a clean curve -> CURVE-MATCHABLE to
  atomistic (option 4). CLAIMS AVAILABLE (ranked): (1) correct thermodynamics + nominal sampling time (best
  for REMD/metad, both are equilibrium-sampling); (2) correct STATIC MECHANICS -- kappa (needs bigger patch),
  K_A (needs NPgT), APL/thickness/S/tilt (have them, ~DOPC) = the strongest positive physical claim;
  (3) overdamped path for D specifically (friction calibrated to D); (4) curve-match MSD/F(q,t)/rot-ACF to
  atomistic for an effective-time factor without a clean D. NOT viable: any single-scalar transport timescale.

## Update 57 (2026-07-19): hydrophobic thickness right; K_A NOT measurable with the current barostat

Hydrophobic core d_c (C1A/C1B leaflet-leaflet) = 27.8 A = 2.78 nm ~= DOPC (2.7-2.9). GL-GL 35.2, PO4-PO4
43.3 -> the P-P excess (findings 56) is HEADGROUP projection (PO4 ~4 A above glycerol/side), NOT the
hydrophobic core -> the functionally-important (TM hydrophobic-matching) thickness is correct.
K_A ATTEMPT FAILED (tooling, not physics): zero-tension NPgT (barostat type=1 "PR"), g-JF, 200000 steps.
The box area barely fluctuated -- sd(A)=0.95 A^2 vs the ~81 A^2 expected at K_A~265 (85x suppressed) -> K_A
came out 1.9e6 mN/m (nonsense); pressure never reached target (P_xy 9.4e-5 vs 2.07e-5). ROOT CAUSE (box.cpp):
the "Parrinello-Rahman" path is NOT a true extended-Lagrangian barostat -- it applies an ad-hoc box-velocity
damping (0.95/step) + tight scale clamps [0.99,1.01], and Berendsen is clamped harder. Both are damped/clamped
RELAXATION schemes for equilibration; neither samples the NPT fluctuation spectrum -> K_A-from-fluctuations is
unavailable, and the weak coupling makes tension-area unreliable too (box didn't respond). Mean APL under the
barostat (65.5) == the NVT APL (65.6) -> the FF's zero-tension APL is ~65.5 (good, ~DOPC). To get K_A/kappa
would need a fluctuation-correct barostat (proper PR or an MC-barostat) or the NVT stress-fluctuation/Born
method (needs the per-frame virial exposed). Polymer-brush kappa = K_A*d_c^2/24 is then one step from a valid
K_A; since d_c is already correct, kappa tracks K_A. Flagged; not fixed (barostat rewrite is a real task).

## Update 58 (2026-07-19): "do all 4" -- MC barostat, temp-sweep, metadynamics node, curve-match, analysis

E1 MC BAROSTAT (box.cpp BarostatType::MonteCarlo, type=2): exact-NPT Metropolis, scales MOLECULE COMs
(reads /input/molecule_ids) so dU is intermolecular (area elasticity, not bond stiffness); one xy + one z
trial/call; accept exp(-b[dU+P dV])*(V'/V)^N_mol. RESTORES area fluctuations (sd(A) 0.95->29 vs the broken
damped-relaxation barostat). BUT K_A undersampled: block-Var grows with window (K_A 5658->4337->2693->1867
mN/m, still falling) -> slow area relaxation coupled to the sub-diffusive g-JF lipids; converged number needs
long sampling (a 500k-step run is running). Zero-tension APL relaxes 65.6->60 (FF MORE CONDENSED than the NVT
box implied -> consistent with too-thick/over-ordered). kappa tracks K_A. main.cpp passes sys.temperature.
E2 TEMP SWEEP (temp_sweep.py, g-JF, T_up 0.5-1.5) -- COMPLETE, all 8 rungs: STABLE at every rung (0.5-1.5);
D monotonic in T (Arrhenius), D@27 = 18.8/30.3/39.9/41.8/46.9/52.7/60.6/70.5 at T=0.5/0.7/0.86/0.9/1.0/1.1/
1.3/1.5; ~303K rung D@270~13 um2/s ~ DOPC. g-JF is stable + physically-trending across the whole REMD ladder
(supports task 2). NOTE thickness grows mildly with T here (40.9->47.5) -- a short-run/expansion detail.
E3 METADYNAMICS (src/metadynamics.{h,cpp} PotentialNode<1>; well-tempered Gaussians summed analytically,
dV/ds->cv.sens; deposit hook metadynamics::maybe_deposit in main.cpp; CMakeLists x2; advanced_config.
write_metadynamics; logs output/metadynamics_cv {2} + metadynamics_bias {grid_bins}). SMOKE PASS: biasing a
bead z drove the CV over 9.8 A (vs ~1 A confined), well-tempered bias 0->6.9 E_up, FES logged. Quantitative
FES-vs-unbiased-PMF: 30k-step run NOT converged (atom escaped the well; RMSD 1.2>kT) -- node correct, needs a
longer/bounded run. Bilayer bit-identical after ALL C++ changes (clamp fix + main.cpp REMD temp + MC barostat
+ metad) -> zero regression to master behavior.
E4 CURVE-MATCH (curve_match.py): CG g-JF DOPC director rot-ACF vs all-atom CHARMM36 (134 DOPC, 11.5 ns @
100 ps, ~/Downloads/charmm-gui-8430021493). POOR overlay (rms 0.14; shapes differ -- all-atom director carries
fast tail flexibility the stiff CG tail lacks). Confirms quantitatively: NO single effective-time factor maps
CG lipid time to physical (rotation vs diffusion vs curve-shape all disagree = mode-dependent clock).
E5 ANALYSIS (example/00.AnalysisScripts/helpers/): analyze_metadynamics.py (FES=-(g/(g-1))V_bias from the
logged bias grid; TESTED) + analyze_remd.py (per-slot demux, per-T observables, MBAR F(CV) at target T with
pymbar>=4 weights() + histogram fallback; TESTED on a real 4-replica REMD run) + README. REMD end-to-end
run PASSED (bl.run.*.up, 4 slots fixed T 0.8-1.1, <E> rises with T) -> REMD works with the g-JF lipid.

## Update 59 (2026-07-19): tex fully rewritten (coherent), thesis reframed

drymartini_upside_interface.tex rewritten end-to-end (not patched): thesis is now ML protein FF (Upside,
data-rich) + PHYSICAL membrane (dry-MARTINI), coupled via RAW cross-interactions -> NO membrane-protein
training data needed -> reaches the data-sparse regime (intro contrasts the prior group's ~40-protein
implicit-membrane model). Scope made honest + central: STATICS/thermodynamics rigorous (g-JF + new MC
barostat both Boltzmann/NPT-correct -> structure, mechanics, conformational ensembles, local stability);
DYNAMICS effective (new Sec: lipid sub-diffusive in BOTH integrators, mode-dependent clock, NO single
effective-time factor -- corrects the old "R~7 recovers D [VALIDATED]" claim to [EFFECTIVE]). HDX is the
flagship (new Sec): EX2 protection factor = equilibrium K_op = thermodynamic -> in the trustworthy domain;
break/reform RATES cancel; caveats = converge rare opens via REMD/metad + "open AND water-accessible"
criterion in a dry membrane. New Outlook: compose-by-physical-coupling generalizes to rare-data systems
(nanoparticles). Compiles clean (pdflatex, 260 KB).

## Update 60 (2026-07-19): secondary-structure failure and coupled clock, corrected diagnosis

The former `friction_scale=100` conclusion was a system-dependent thermostat fit, not a physical
protein--membrane drag model, and is superseded. The actual failure combined three independent defects:

1. O and the MARTINI BB centre are derived from N/CA/C, but the old force path regenerated them by mutating
   the integrator position array. A `martini_hybrid_position` CoordNode now constructs O/BB without mutating
   the integrated coordinates. Its reverse pass follows Upside's regenerated-O contract: it copies the
   14/54, 12/54, and 12/54 BB shares to persistent N/CA/C and discards derived-site sensitivities. O/BB are
   not integrated or thermostatted particles.
2. `/output/pos` used the raw integration array, so its O/BB entries were stale. Uncorrected DSSP therefore
   reported helix loss even when the N/CA/C helix core remained native. The logger now copies the evaluated
   hybrid coordinate node. A direct three-frame check found maximum O and BB mapping errors of
   `1.91e-6 A`; ordinary trajectory analysis and stage restarts now see the force-graph geometry.
3. The soluble-protein step `dt=0.009` under-resolves the hard coupled interface. Equal-duration 1RKL
   regressions (45 Upside units, all SC-env and BB-env interactions active) gave:

   - `dt=.009`: carrier kT mean/peak 1.160/2.539; H-bonds 32.79 -> min 7.12 -> 32.25.
   - `dt=.0045`: carrier kT mean/peak .872/1.101; H-bonds 32.79 -> min 20.23 -> 30.50.
   - `dt=.00225`: carrier kT mean/peak .945/1.480; H-bonds 32.79 -> min 22.43 -> 32.54.

The current logger rerun at `.00225` was bitwise identical for every real coordinate and momentum; only
the formerly stale virtual positions changed. Direct raw-output DSSP found 23 helical residues initially,
24 finally, and all 23 native helical sites retained at the endpoint. O/BB mapping maxima were
`1.15e-6/1.06e-6 A`, versus `10.76/11.38 A` with the old logger.

The hybrid default is therefore `dt=.00225`. Default step counts and frame strides are multiplied by four
to preserve the previous nominal duration and output cadence. The fix does not soften, disable, or cap any
SC-env or BB-env interaction; all runtime pair forces remain spline-table evaluations.

Fresh current-code 1AFO at `dt=.00225` for the same 45 units retained the native helix core at 45/45
(minimum 44/45), total DSSP helix 54 -> 54, H-bonds 93.44 -> 83.36, and Rg 15.79 -> 15.50 A. Its final
nine-unit block had carrier/lipid kT .917/.836. Saved O and BB mapping maxima were `1.90e-6` and
`1.21e-6 A`, and O/BB momenta were exactly zero. At `dt=.009`, the comparable final-block carrier kT,
H-bonds, and native core were 1.139, 68.51, and 43.8/45.

Thermostatting now follows one-bath-per-particle fluctuation--dissipation bookkeeping. Environment beads
use absolute g-JF friction `alpha=m/(4 ps / 0.202896643 ps)` (mass-6 alpha `0.304344964`); positive-friction
particles are excluded from global OU. Protein N/CA/C carriers use the zero-friction Verlet limit of the
same g-JF update and retain the standard Upside OU bath. The mechanical step is 0.457 fs. The often quoted
40 ps per Upside time unit is an empirical protein conformational acceleration, not the mixed-system clock.

This gives the bilayer the standard dry-MARTINI stochastic relaxation and lets the protein receive
stochastic drag through conservative collisions with that damped bath. It does not establish a quantitative
protein resistance tensor: exact parallel, normal, and rotational protein--membrane drag would require a
momentum-conserving pair thermostat calibrated to a conventional dry-MARTINI or atomistic reference. No
uncalibrated DPD/contact coefficient is introduced merely to make secondary structure survive.

## Update 61 (2026-07-20): literal 40 ps protein-step clock requires a synchronously coupled lipid micro-solver

The requested protein clock is one integration step = 40 ps physical, not 40 ps per Upside time unit. Code
inspection adds a critical scheduling fact: `DerivEngine::integration_cycle` advances the protein through
three consecutive stages, but the current lipid propagator runs only once on stage zero. Thus the mixed
system already has an explicit 3x schedule mismatch. With the standard MARTINI factor four
(`t_physical = 4 t_CG`), the lipid must represent 10 ps raw MARTINI dynamics per protein stage, or 30 ps per
three-stage cycle. For `D_DOPC = 11.5 um^2/s` at 303 K, the acceptance target is `0.184 A^2` lateral MSD per
stage and `0.552 A^2` per cycle. Equivalently, raw-CG D is 46 um^2/s. Multiplying again by four would
double-count the correction.

The current single-step g-JF path cannot represent this clock, and the existing overdamped fallback is not
adequate for protein friction: it advances the lipids through all inner Brownian force refreshes while the
protein carriers receive only the final `pos->sens`. The time-integrated fluctuating protein--lipid reaction
is discarded. The clean architecture is a mixed multirate SDE: on every inner step, update lipids with
constant-mobility overdamped Brownian dynamics and FDT noise, and advance protein momenta/coordinates
symmetrically from the simultaneous full conservative force. Every SC-env and BB-env reaction follows the
hybrid force-routing contract. The proposed joint canonical distribution is invariant; eliminating the lipid
coordinates gives the protein a colored memory/friction kernel without an empirical protein drag multiplier.

Calibrate only lipid mobility on a large pure-DOPC patch to the 0.184 A^2/update slope; choose the inner-step
count by convergence, not by diffusion. The existing hard-core tests suggest about 90 overdamped inner steps
for the hybrid, while an inertial dry-MARTINI reference spanning 10 ps would require 250--500 conventional
40--20 fs steps. Do not use the current displacement cap in the calibrated path. Validate a Fickian MSD
window plus APL, thickness, order, and compressibility, then measure protein lateral/normal/rotational
resistance and membrane force/torque correlations. Add a calibrated pair/contact GLE only if this emergent
friction misses those independent reference observables.

## Update 62 (2026-07-20): standard `.009` single-step constraint and molecular-versus-bead target

The user clarified that both components should use the same numerical timestep as the ordinary Upside
examples and that bilayer multistepping should be avoided if possible. `example/02.ReplicaExchangeSimulation`
does not override the engine default, so its timestep is `.009`; the ordinary Verlet path groups three `.009`
steps per outer bookkeeping cycle. The next implementation must therefore test a one-step all-particle path
at `.009` before considering any lipid subcycling.

The primary transport target remains whole-DOPC COM diffusion: at 303 K, 11.5 um^2/s and 40 physical ps per
`.009` step require `0.184 A^2` lateral MSD per numerical step. The MARTINI-equivalent statement is 46 um^2/s
over 10 raw-CG ps per step. For a bonded molecule, every bead's long-time trajectory is the COM trajectory
plus bounded internal motion, so every bead and the COM have the same asymptotic diffusion coefficient. A
separate sub-molecular fallback can only calibrate the *bare/free bead mobility*. Its initial g-JF friction is
`alpha_bead = kT_ref*dt/(D_phys*40 ps)`; at `kT=.8647`, `dt=.009`, and 11.5 um^2/s this is about 0.1692 in
Upside units. A free-draining molecular initial estimate is `alpha_bead/14`, but hard-core caging makes the
actual DOPC COM diffusion a measured, nonlinear quantity.

Lesson: do not assign physical time to outer logging cycles or infer a bead-specific long-time diffusion in
a bonded lipid. Verify the numerical timestep against the reference example, state whether the observable is
molecular COM diffusion or bare bead mobility, and apply the MARTINI factor four exactly once.

## Update 63 (2026-07-20): molecular-D scan fails; use an explicit protein--bilayer friction clock

The user refined the claim hierarchy: prefer the factor-four-corrected molecular DOPC lateral diffusion; if
40 ps/step is unreachable, define the kinetic mapping through the friction experienced by protein in the
bilayer. A one-step `.009` 1RKL scan used lipid frictions 0.3043 (old standard), 0.1692 (free-bead D target),
0.0121 (14-bead free-draining molecular estimate), 0.0035, and 0.001. Inferred whole-DOPC COM diffusion was
0.055, 0.084, 0.155, 0.168, and 0.175 um^2/s respectively, versus 11.5 target. The low-friction response has
already reached a Kramers/cage plateau about 66x too slow, while protein/lipid kinetic temperatures rise and
H-bonds degrade. Friction reduction cannot produce the molecular clock in a stable single step.

The selected fallback maps the known dry-MARTINI 4 ps raw stochastic relaxation rather than pretending that
the measured COM diffusion is correct. Factor four makes that 16 ps physical. With one `.009` integration
step representing 40 ps physical (10 ps raw), `tau_up=.009*(4/10)=.0036` and `alpha_i=m_i/tau_up`. The same
mass-scaled g-JF bath is applied to environment particles and to real N/CA/C protein carriers whose initial
positions lie within the existing 12 A spline interaction cutoff of a lipid bead; noncontact carriers retain
the standard Upside OU bath. This is FDT-consistent and gives the protein a literal bilayer-friction clock,
but it is explicitly not a claim that DOPC COM diffusion equals 11.5 um^2/s.

A 5000-step upper-bound test applying the mapped friction to every 1RKL carrier reduced mean/final protein
kinetic energy from 1.69/1.72 (lipid-only friction clock) to 1.29/1.21 and raised the minimum H-bond score
from 16.0 to 19.8. The production mask uses the physical 12 A interaction cutoff instead of damping
noncontact protein indiscriminately.

The final implementation removes the overdamped/substep branch and its displacement cap. Every physical
particle takes one g-JF update per `.009` step. H5 records the full clock contract and the exact contact
carrier indices. Preparation uses ordinary native dry-MARTINI damping; switching to the high production
friction before resolving packing overlaps froze the bath and failed, so the friction clock activates only
for `npt_prod`. Fixed and z-fixed constraints are enforced inside g-JF; regenerated O/BB sites remain
nonintegrated with exactly zero momentum.

Fresh final-code workflows passed for both proteins. In 5000-step stage-7 runs, mapped all-atom DSSP helix
counts stayed exactly 7 for 1RKL and 27 for 1AFO. Mean protein/lipid kinetic energies were 1.286/1.297 and
1.300/1.298, respectively, versus `3 kT/2 = 1.297`. Current short-window molecular DOPC estimates were
`7.13e-5` and `7.00e-5 um^2/s`, versus 11.5 target. Separate 50,000-step continuations on the same production
integrator also kept helix counts fixed and averaged 99.9%/100.2% of target kinetic energy; their molecular
diffusion estimates were `6.72e-5` and `6.04e-5 um^2/s`. SC-env and BB-env spline nodes remained present,
and no conservative potential was modified.

Lesson: when diffusion and friction clocks disagree, name the calibrated observable in H5 and in the paper.
Never report a friction-calibrated trajectory as having the target molecular lateral diffusion.

## Update 64 (2026-07-20): user-reported stage-7 freeze reopens the friction design

The user found that the delivered 1AFO and 1RKL stage-7 VTF trajectories are effectively frozen. The prior
acceptance gates over-weighted finite coordinates, secondary-structure retention, and canonical kinetic
energy. Those checks can all pass for a high-friction g-JF process whose momenta thermalize but whose
coordinates barely move. The production value `tau_up=.0036` gives `alpha*dt/m=2.5`, so mobility and actual
coordinate displacement must be measured directly before any friction-clock claim is accepted.

Lesson: never accept a kinetic calibration from temperature and structural stability alone. Every friction
change must gate on protein displacement/orientation and whole-molecule lipid MSD over the saved trajectory,
and the reported VTF must be inspected rather than only the source H5 statistics.

## Update 65 (2026-07-20): particle-mobility clock removes the freeze without melting the TM core

The saved broken trajectories quantify the regression. Drift-removed DOPC molecular COM motion was only
`0.0081 A` per saved frame and `0.17--0.24 A` net RMS over 1,001 frames. The cause is the rejected production
mapping `tau_up=.0036`: for a mass-6 bead, `alpha=1666.7` and `alpha*dt/(2m)=1.25`. g-JF still produces a
thermal-looking momentum distribution in this regime, so temperature was a false acceptance signal.

Three nonfreezing controls were compared for 5,000 steps. Native dry-MARTINI damping (`alpha=0.3043`, no
protein g-JF friction) and the previous single-step regime (`alpha_lipid=.0035`, `alpha_protein=.35`) both
moved DOPC but heated the interface and reduced 1RKL helicity. Direct bare-particle friction on lipid-contacting
protein carriers likewise destabilized 1RKL. The accepted model uses the user-approved sub-molecular fallback:
`D_raw=4*11.5=46 um^2/s`, `dt_raw=40/4=10 ps`,
`D_bead,up=D_raw*1e-4*dt_raw/.009=5.1111 A^2/U`, and
`alpha_bead=kT/D_bead,up=0.1691804`. Each environment bead receives this friction. A real protein N/CA/C
carrier receives the additive friction `n_contact*alpha_bead`, where `n_contact` is the number of DOPC beads
inside the existing 12 A spline cutoff. Counts are refreshed after stage handoff, minimization promotion, and
production continuation, then held fixed during each segment so the SDE does not silently acquire
position-dependent multiplicative noise.

Matched 50,000-step tests at the shared `.009` timestep passed without substeps. DOPC COM net RMS motion was
`3.43 A` (1AFO) and `3.72 A` (1RKL), compared with less than `0.24 A` in the broken outputs. Measured molecular
diffusion remained only `0.0132` and `0.0152 um^2/s`; this is reported as a failed molecular target, not hidden.
Total kinetic energy was within 1% of `3kT/2`; protein subset means were within 1.3%. 1AFO finished with all
54 helical residues. In 1RKL, residues 10--29 retained mostly 0.9--0.99 helix occupancy; the N-terminal
extension at residues 6--9 relaxed to coil, leaving a 19-residue transmembrane helix. Valid fresh workflows
for both proteins completed normal minimization, equilibration, burn-in, production, and VTF extraction with
finite moving coordinates and correct H5 friction/contact metadata. Preparation uses native damping in
softened stages and strong FDT overlap-settling damping in full-core stages 6.1/6.4--6.6; neither bath is a
kinetic claim. Production alone uses the particle clock. A deliberately abbreviated workflow that disabled
minimization exploded from unresolved packing overlaps and is not evidence about the friction model.

Lesson: derive fallback friction from the named observable, propagate it to protein through a documented
local physical approximation, and validate both coordinate transport and structural occupancy on the final
trajectory. Never skip hard-core preparation merely to shorten an end-to-end test.

## Update 66 (2026-07-20): HDX is an equilibrium structural proxy; the hybrid path is not production-ready

The `example/00.AnalysisScripts` uptake path does not read simulated elapsed time as exchange time. For each
amide donor, `get_protection_state.py` assigns a binary protected flag from backbone H-bond score, an Asp/Glu
side-chain-contact proxy, and backbone/side-chain burial (plus an optional membrane-surface term that the driver
does not enable). `4.calc_D_uptake.py` MBAR-reweights those flags to `p_protected`, calculates sequence-, pD-,
and temperature-dependent intrinsic `k_chem`, then applies the EX2-like approximation
`k_obs = k_chem * (1-p_protected)` and `D(t)=1-exp(-k_obs*t)` in experimental seconds. Thus a wrong friction or
40 ps mapping does not directly rescale the HDX time axis; it matters because it controls decorrelation and the
ability to sample opening/closing equilibria. The old frozen trajectories fail that requirement.

The stock extraction path cannot consume the hybrid coordinate layout. `mdtraj_upside.traj_from_upside` asserts
`n_coord = 3*n_res`, whereas the stage-7 hybrid files have 2,894 coordinates for 72-residue 1AFO and 4,098 for
31-residue 1RKL. The standalone `py/martini_protection_state.py` understands full hybrid coordinates but is not
called by `analysis.sh` or `3.get_protaction_states.sh`; no `PS.npy`, percentD, or HDX output currently exists
under `example/16.MARTINI`.

The current hybrid criterion also saturates membrane protection. It marks an amide protected whenever its N is
between the two global PO4 planes, independently of local hydration or defects. On the 11 stage-6.6 frames,
1RKL has 29/29 donors protected in every frame although 65.5% of donor-frame observations are H-bond-open. The
uptake code handles exact `p_protected=1` with an empirical floor `k_obs=k_chem/1000`, so such a prediction would
mostly report the geometric membrane rule and the 1000-fold floor, not a sampled protein free-energy signal.
For 1AFO the mean exchange-competent fraction is 0.132, but 65/67 donors never change state in those frames.
These preparation-frame numbers diagnose criterion behavior only; they are not production HDX estimates.

The corrected 50,000-step production length is `50000*40 ps = 2 us` under the requested literal
40-ps-per-numerical-step protein clock. The earlier 18-ns statement incorrectly multiplied the empirical clock
by the reduced numerical step a second time. Two microseconds is materially better sampling, but a single
trajectory is still not evidence of converged rare opening populations.
Quantitative hybrid HDX requires a wired hybrid extraction path, a defensible membrane/water-accessibility model,
replica/enhanced sampling, protection-probability convergence and effective-sample-size checks, and external
peptide-level validation. Until then, use any calculated profile only as a qualitative structural-protection map.

Lesson: before judging HDX from a dynamics timescale, inspect whether the estimator is kinetic or equilibrium.
Then verify coordinate-layout compatibility and whether binary protection rules saturate; a mathematically valid
uptake curve can still contain almost no trajectory-derived information.

## Update 67 (2026-07-20): 1RKL runs at the intended 303 K; endpoint DSSP hid reversible helix loss

The completed 1RKL production invocation uses `--temperature 0.8647`, and all 1,001 H5 frames store that value.
With `1 T_up=350.588235 K`, this is 303.15 K. The last 200 frames give protein and lipid mean kinetic energies
1.011 and 1.003 times `3 kT/2`, respectively. g-JF receives the same `kT` as the OU thermostat, and positive-
friction degrees of freedom are excluded from OU. There is no evidence for a missing factor of four or another
integrator-dependent temperature conversion; the MARTINI factor four changes time, not thermodynamic temperature.

Direct DSSP from the mapped N/CA/C/O coordinates confirms the user's visual correction. Residues 10--29 average
84.5% helical occupancy over the 2-us trajectory. They briefly fall to four helical residues near frame 543
(about 1.09 us on the declared clock), then recover to 19/20 in the final frame. The core CA RMSD from frame zero
is only 0.99 A on average and 1.58 A at maximum, so this is a reversible local hydrogen-bond/geometric disruption,
not global unfolding. The central alpha2 region (including residues 18--24 identified in the experimental study)
is the most persistent. Endpoint-only statements that the core “remained helical” were insufficient and are
superseded by occupancy and minimum-event reporting.

Correct kinetic temperature does not prove the nonlinear hard-core configurational distribution is converged at
finite timestep. The existing matched-duration study already found `.009` under-resolves the coupled interface,
whereas `.00225` retained the native helix. g-JF has favorable harmonic configurational statistics but no exact
finite-step guarantee for this nonlinear hybrid force field. Therefore lowering the thermostat temperature to
hide helix loss is not defensible. Match the experimental temperature (`T_up=T_exp/350.588235`) and establish the
largest acceptable timestep with `.009/.0045/.00225` structural-population comparisons. The older Upside membrane
examples use `T_up=0.80` (280.5 K); that is a different thermodynamic condition, not an integrator correction.

For HDX, the clean hybrid adapter is to extract N/CA/C coordinates through `/input/hybrid_bb_map/atom_indices`,
evaluate those coordinates with a separate protein-only `-HDX.up` analysis engine to recover the stock backbone
H-bond, Asp/Glu-contact, and protein-burial terms, and combine them with membrane accessibility evaluated from the
full hybrid frame. The current global-phosphate-plane Boolean is not quantitative: it permanently protects a
stable transmembrane donor and saturates 1RKL. Since dry MARTINI contains no water, local water accessibility must
come from a calibrated depth/water-activity function or an explicit-water reference; otherwise membrane-corrected
HDX should be reported only as a sensitivity analysis.

The adapter must save the full hybrid potential for MBAR, not the protein-only analysis-engine energy. It must
also save `T.npy` in Upside `kT` units. The `example/00.AnalysisScripts/README.md` statement that `T.npy` is Kelvin
is incorrect: `4.calc_D_uptake.py` constructs `beta=1/T` directly and expects values such as 0.85; only the
intrinsic-chemistry branch converts that target to Kelvin. Feeding 303 instead of 0.8647 would invalidate MBAR.

Lesson: validate secondary structure with per-residue occupancy, time blocks, and worst excursions, not endpoint
counts. A correct kinetic-temperature average cannot substitute for configurational timestep convergence.

## Update 68 (2026-07-20): maximize reuse by projecting a standard HDX trajectory view

The user corrected the architecture priority: the hybrid analysis must maximize reuse of `example/04.HDX` and
`example/00.AnalysisScripts`, not fork their formulas. The narrow adapter should create a protein-only HDX-view
H5 from each hybrid replica. Its `/input` comes from the ordinary `-HDX.up` config; each output group's positions
are N/CA/C mapped through `hybrid_bb_map/atom_indices[:, :3]`. Full coupled-system potential, Upside temperature,
and compatible H-bond logs are copied from the hybrid group. This makes `mdtraj_upside`,
`get_info_from_upside_traj.py`, and `get_protection_state.py` see their native `3*n_res` contract while MBAR still
uses the correct protein+bilayer Hamiltonian.

`example/04.HDX` remains the reference for REMD, MBAR-weighted protection probabilities, residue dG, and
denaturant dependence. `example/00.AnalysisScripts` is the maintained superset for experimental preprocessing,
intrinsic chemistry, peptide uptake, stability summaries, comparison, and plotting. New hybrid-specific code is
limited to the H5 projection and an optional, separately validated membrane-accessibility combination. The stock
protein-only `PS.npy` must be retained beside any combined PS so the membrane correction is observable and
reversible.

Lesson: when adapting a new trajectory representation to a mature analysis workflow, translate the data into the
existing contract. Do not duplicate the established estimator merely because the source coordinate layout differs.

## Update 69 (2026-07-20): hybrid HDX adapter implemented; current 1RKL fails quantitative gates

The implemented projector creates the standard `3*n_res` H5 view and preserves the full coupled observables.
On all 1,001 1RKL frames, projected N/CA/C coordinates are exactly equal to the hybrid-map selection, and
potential, H-bond, temperature, and time arrays are exactly equal to the source. The ordinary
`get_protection_state.py` path then produces `(1001,29)` PS. Without an external accessibility array,
`PS_protein.npy` and downstream `PS.npy` are exactly identical. A synthetic check confirmed
`P_protected,total = 1-(1-P_protected,protein)*P_water_accessible` exactly.

The full single-replica workflow now completes through uptake, stability, and summary. Two latent workflow
defects reached by this case were fixed without changing the estimator: scalar `T.npy` is accepted for MBAR,
and stability weights use a max-shift before exponentiation. The hybrid wrapper reads the sampled `T_up` from
replica 0 and sets all single-temperature targets to that value; the current trajectory is analyzed at
`0.8647`, not the legacy unsampled `1.14` default.

The generated numbers remain non-quantitative. Current molecular DOPC diffusion is `0.0152 um^2/s`, about 750
times below the `11.5 um^2/s` target. More decisively for this equilibrium estimator, first/second-half donor
protection differs by `0.078` on average and `0.621` at residue 10; the smallest estimated effective sample
count is only `3.74`. Dry MARTINI also lacks a calibrated water-accessibility field. The friction-clock mismatch
does not rescale experimental seconds in the EX2 calculation, but it severely limits configurational mixing, and
the protein's additive contact bath has not been validated as the physical bilayer hydrodynamic friction.

Lesson: a mechanically correct HDX pipeline does not validate its source ensemble. Gate trust on per-residue
block convergence and effective sample size, sampled-temperature matching, timestep convergence, membrane
accessibility calibration, and independent-replica or REMD agreement.

## Update 70 (2026-07-20): low temperature does not exclude a BB-env force regression

The user reports that 1RKL still loses part of its alpha helix at `T_up=0.80`, whereas older revisions such as
`b1041bb6640b34edf55fe371b5cf943d054ba187` did not. This supersedes any inference that the observed helix loss
can be explained primarily by the 303 K condition. The audit must compare the complete particle--backbone force
path across revisions: source table values and units, H5 type mapping, pair classification/exclusions, spline
derivative sign and scale, virtual-backbone geometry, and Jacobian propagation to N/CA/C/O carriers.

Lesson: when a structural regression persists after lowering temperature and an older revision is stable, do
not keep interpreting it as thermodynamic behavior. Establish trajectory provenance and bisect the full force
path against the known-good revision before changing temperature, friction, or structural restraints.

## Update 71 (2026-07-20): cited trajectory predates the `.80` temperature edit

The cited `martini_1rkl_hybrid_full/1rkl.stage_7.0` output records `--temperature 0.8647` and stores
`0.864700019` in every `/output/temperature` frame. Its `0.80` H5 value is the DOPC Brownian-friction reference,
not the runtime thermostat temperature. The current workflow script has since been edited to default
`TEMPERATURE=0.8`, so a fresh output is required before evaluating the user's low-temperature observation.

The force audit found no raw particle--BB spline regression relative to `b1041bb`: the canonical pair grids are
identical, and cutoff, type assignment, force sign/scale, pair exclusions, and Newton reaction are unchanged.
The material derivative-path regression was the full BB virtual-site Jacobian replacing the historical static
mapping, whose disposable O share leaves 70.37% of the translational BB reaction on persistent N/CA/C coordinates.
The older default also used `dt=.002`, whereas the cited current output used `.009`. Per user direction, neither
path will be changed until a fresh `.80` run separates temperature provenance from these other differences.

Lesson: verify the completed H5 invocation and `/output/temperature` before interpreting a script's current
defaults or a friction-reference temperature as the temperature of an already-generated trajectory.

## Update 72 (2026-07-20): one temperature controls the full hybrid system

`run_sim_hybrid.sh` previously allowed `TEMPERATURE` and `UPSIDE_DOPC_REFERENCE_TEMPERATURE_UP` to differ. The
former controls the actual Upside thermostat for the coupled system; the latter enters the DOPC mobility/friction
calibration. A mismatch does not create a second thermostat, but it calibrates friction at one `kT` and drives its
noise at another, changing the nominal diffusion. The workflow now assigns the DOPC reference directly from the
single authoritative `TEMPERATURE` and overwrites any independently supplied reference-temperature environment
value. `run.py` already enforces the same relationship.

Lesson: temperature-dependent friction calibration and runtime FDT noise must derive from one system-temperature
setting unless a deliberately temperature-dependent transport model is explicitly implemented and validated.

## Update 73 (2026-07-20): friction derivation and trust boundary consolidated

The interface manuscript now follows the implemented calibration directly. With $f_\mathrm{CG}=4$,
$D_\mathrm{target}=11.5~\mu\mathrm{m}^2/\mathrm{s}$, 40 ps declared protein time, and $dt=.009$, the raw
mapping is 10 ps and $46~\mu\mathrm{m}^2/\mathrm{s}$. This gives
$D_\mathrm{bead,up}=5.11111~\text{\AA}^2/t_\mathrm{up}$ and, at the unified $T_\mathrm{up}=.8$,
$\alpha_\mathrm{bead}=.15652$. Protein carrier friction is the static per-segment contact count times this bead
friction. This is an additive local FDT damping approximation; it omits correlated, anisotropic, hydrodynamic,
and memory effects and is not a measured protein resistance.

The evidence boundary is now explicit. One-temperature wiring, one-bath bookkeeping, metadata, virtual-site force
projection, and HDX representation mapping are implementation results. Bare-particle mobility is calibrated by
construction. Whole-DOPC diffusion fails by roughly three orders of magnitude; the 40 ps clock cannot be assigned
to collective trajectory kinetics. Bilayer/protein equilibrium populations remain conditional on timestep,
charged-cutoff, ensemble, and replica convergence. Current HDX is not quantitative because its source ensemble
fails block/ESS gates and dry MARTINI lacks calibrated water accessibility. The unequal kinetic clock does not
rescale EX2 seconds directly; it undermines sampling of the protected probability.

Lesson: organize validation by implementation, configurational distribution, and kinetics. Never promote an
equilibrium-form estimator to a trustworthy prediction before its finite-step ensemble and sampling convergence
have passed.

## Update 74 (2026-07-20): unified-temperature pre-fix trajectories (superseded structurally)

Fresh stage-7 files replaced the previously audited outputs while the manuscript was being revised. Both 1RKL and
1AFO now contain 1,001 frames at runtime $T_\mathrm{up}=.8$ with friction-reference $T_\mathrm{up}=.8$,
$D_\mathrm{bead,up}=5.11111$, and $\alpha_\mathrm{bead}=.15652$. Last-200-frame protein/lipid kinetic-energy
ratios to $3kT/2$ are 1.005/0.991 for 1RKL and 1.002/1.018 for 1AFO, so the unified thermostat works as specified.

These files still used the full-Jacobian BB reverse route and the unprotected handoff, so their structural
measurements are not results of the final interface. They were nevertheless sufficient to isolate temperature
and transport: drift-removed DOPC COM fits were window-dependent and far below target, with 1RKL giving
$0.0053$--$0.0112~\mu\mathrm{m}^2/\mathrm{s}$ and 1AFO
$0.0027$--$0.0103~\mu\mathrm{m}^2/\mathrm{s}$ over tested 20--1000 ns windows, versus 11.5 target. Their failed
1RKL DSSP result motivated the force-route and handoff audit and is replaced by Update 75.

Lesson: when outputs can be regenerated concurrently, re-read artifact provenance immediately before finalizing
result-specific documentation; a correct historical diagnosis can become stale during the same task.

## Update 75 (2026-07-20): BB routing follows regenerated-O semantics; pre-production is rigid

The user corrected the integration contract behind the historical reverse map. In `b1041bb`, a BB gradient was
copied with N/CA/C/O mass fractions 14/54, 12/54, 12/54, and 16/54. Upside overwrites O during the integration
cycle instead of integrating it, so the O share is intentionally disposable. The persistent backbone must
therefore receive only 38/54, or 70.37%, of the translational BB reaction. Although a full coordinate Jacobian
is the derivative of the reconstructed virtual site, it is not the requested Upside regenerated-O force-routing
contract. The current coordinate node now reconstructs O and BB in the forward pass, copies only the N/CA/C
mass-fraction shares in reverse, and discards sensitivity on derived O and BB sites.

Stage-resolved DSSP localized the original loss: all saved frames in stages 6.0--6.6 had 23 helical residues,
but the obsolete stage-7 trajectory began at 13 after an unprotected production-interface minimization. The
non-negotiable preparation rule is now stronger than a fixed minimization or finite positional spring. BB-env
and SC-env activate in `production_handoff`, while the complete protein remains one rigid body through both
minimization and burn-in. The stage is relabelled `production` only afterward, which removes the rigid group and
starts flexible sampling. Full corrected 50,000-step runs now complete this sequence. Handoff pair-distance
changes are at most 0.017 A for 1RKL and 0.0037 A for 1AFO. Corrected DSSP first/final/minimum/mean helical counts
are 23/23/16/21.88 for 1RKL and 54/51/48/50.77 for 1AFO; 1RKL residues 10--28 are individually helical in
96.5--100% of frames. The structural regression is resolved.

The engine and ordinary example infrastructure confirm the standard Upside numerical step is `.009`; `.09`
would be tenfold too large. The official Martini 2 parameter contract is a separate open issue: LJ is shifted or
force-switched from 0.9 to 1.2 nm and electrostatics are shifted/reaction-field treated to 1.2 nm. The current raw
hard cutoff is therefore not force-field-equivalent and should be corrected in a dedicated table-validation
phase rather than mixed into the reverse-map fix.

Lesson: do not infer the intended reverse force map from virtual-site calculus alone when the host integrator
deletes and regenerates one of the mapped atoms. Audit which force shares survive the integration cycle. Also,
"protected during minimization" is not equivalent to "rigid before production": enforce one rigid protein body
through every pre-production minimization and MD segment, and release it only at the explicit production handoff.
Once that construct passes a full run, rewrite the primary documentation around the final design and remove the
superseded debugging narrative rather than accumulating successive correction notes.
