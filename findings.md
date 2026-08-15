# Findings

## Health gates: CN_MAX removed from BOTH jobs; count broken bonds instead (2026-08-09)

Measured on a forced NP tear (270-0-0 at dt=0.01), recording when each candidate criterion first fires:

```
FIRST FIRING TIME PER CRITERION
   maxCN (>3.5 A)                   t=168.0
   count (>=5 bonds >2 A)           t=168.0
   potential (non-finite or >1e6)   never
   |coord| (>1e4 A)                 never
```

- **`maxCN` is redundant** — the count fires at the identical frame, so the worst-bond threshold buys
  nothing, and it is the fragile one (healthy max 2.659 A vs torn 3.93 A is a knife-edge; at 2.5 A it
  false-fired on a healthy chunk and cost a 6 h glpG block).
- **An energy test can never guard NP.** In the same tear the protein reached 431 broken bonds with the
  potential still finite at +3e5, so anything keyed on non-finite or large energy never fires.
  Conversely glpG's blow-up goes fully NaN within one 46-step interval, so there `isfinite` suffices.
  The two jobs need different detectors — another reason not to transfer settings between them.
- **The count is a robust discriminant, not a tuned knob:** healthy frames have 0-2 stretched bonds
  (verified: all six NP systems and healthy glpG give 0), a torn one has 279-431. Any cut between 3 and
  200 behaves identically.
- Both drivers now carry **no invented magnitudes** — `CN_MAX`, `POT_MAX`, `COORD_MAX` deleted.
  NP: non-finite positions OR >=5 stretched bonds. glpG: non-finite potential scanned over the whole
  chunk OR >=5 stretched bonds in the final frame.

## Core audit: spline table is CORRECT; one real defect — zero-force hole below 0.1*sigma (2026-08-09)

Audited the MARTINI core after the NP/glpG blow-ups. Most of it checks out; one genuine defect.

**Verified correct (do not re-litigate):**
- The spline table **equals the analytic dry-MARTINI form**: potential-shifted LJ + reaction-field
  Coulomb, both zero at 12 A. Max relative error **3.6e-11 across all 23 (eps, sigma, qq) triples**
  over r=2.5--11.9 A. This satisfies the "Spline Tables Must Be The Original Potential" assertion.
- Watch the convention: GROMACS `epsilon_rf = 0` means **infinity** (conducting boundary), so
  `(eps_rf - eps_r)/(2 eps_rf + eps_r) = 1/2`, NOT -1. Taking the 0 literally makes the charged triples
  look ~100% wrong at large r. That was my error, not the table's.
- No force caps, no NaN-skipping, no energy clamping in the pair-force path (`martini_potential.cpp`
  `eval_pair_force`). The `kMinDistance` values (1e-6 / 1e-8 A) are divide-by-zero epsilons for exactly
  coincident particles, not physics guards.

**REAL DEFECT — `py/martini_build_tables.py:532`: `r = max(r, 0.1 * sig)`**
- The builder floors the radius at `0.1*sigma` before evaluating LJ, so the tabulated potential is
  **constant and the force is exactly ZERO** below 0.47 A (sigma=4.7), 0.60 A (6.0), 0.62 A (6.2).
  Measured: first nonzero slope at r=0.589 A for the sigma=6.0 triple; true LJ there is 7.6e20 vs a
  tabulated 3.16e12.
- This is a hole in the repulsive core and a deviation from the published functional form.
- **UNRELATED to the NP/glpG blow-ups — tested, not assumed.** In healthy dynamics the max per-step
  displacement is ~0.045 A and the closest evaluated pair is 2.403 A, so reaching the hole means
  crossing **1.80 A (~40 steps)** of correctly-tabulated 1e5--1e8 E_up repulsion. Unreachable from a
  healthy state. It only becomes accessible after a cascade has already injected far-above-thermal
  velocities, by which point the configuration is destroyed regardless. If anything the true un-floored
  LJ (U~1e16 at 0.3 A) would eject harder. I earlier called it a "plausible amplifier" — that was
  speculation and is retracted.
- Still worth fixing on spline-fidelity grounds: tabulate the true form down to r->0 (large but finite
  in float64), or state and justify the floor rather than hard-coding it silently. Do not expect it to
  change blow-up behaviour.

## glpG REMD is NOT the same job as NP — stop transferring analysis between them (2026-08-09)

**NP and glpG are different simulations and must be analysed separately.** Conflating them cost a
healthy 6 h glpG block (a threshold borrowed from NP false-positived) and produced a wrong root-cause
writeup.

| | NP (`np_1AO6_prod`) | glpG (`remd_glpG-*`) |
|---|---|---|
| method | regular MD, 6 independent trajectories, single T=0.8647 | **REMD**, 48 replicas T=0.70-0.90 with configuration exchange |
| purpose | NP adsorption footprinting (K190) | **HDX** protection factors |
| integrator | **pure velocity-Verlet**, no `/input/brownian` | **MIXED**: 2736 atoms (ION 432 + LIPID 1674 + PROTEIN 630) overdamped **Brownian**; only 420 protein atoms velocity-Verlet |
| timestep | free; fixed by changing 0.005 -> 0.001 | **hard-locked 0.009** (`martini_brownian.cpp:100` throws on mismatch); friction tuned against it for lipid D=11.5 um2/s. A calibrated quantity — do not change |

- **The glpG mechanism is OPEN.** The earlier claim of "same root cause, dt 8.3x over the limit" applied
  NP's velocity-Verlet stability analysis to a system whose environment is overdamped Brownian. It does
  not transfer. The 8.3x number was independently wrong too (max over all 23 grids and mu=0.5 for pairs
  the engine never evaluates; per-pair it is ~0.56x of the harmonic limit, i.e. marginal at worst).
- What IS solid for glpG: the 79ALA blow-up was real; REMD exchange **recirculates** the destroyed
  configuration indefinitely so it contaminates whichever slot it occupies (1, later 2 configs walking
  the ladder, not "4 failed replicas"); and `isfinite` misses it because env coords reach +/-4.65e12 A
  while staying finite.
- Old, superseded reasoning retained below for context only:
- The blow-ups were hidden because **REMD exchange launders them**: `run_remd.py:59` reseeds from
  `output/pos[-1]` (the destroyed frame), then a swap brings a healthy configuration in from a
  neighbouring slot, so the trajectory looks finite from frame 1 and the bad slot walks (3->13->23->33).
- Two traps this cost me, both worth remembering: (a) the environment coordinates at the failed frame
  are `+/-4.65e12 A` -- numerically **finite**, physically destroyed, so `isfinite` passes; (b) it looks
  like a chunk-teardown logging artifact but is not -- potential/kinetic/hbond/rama are all NaN at that
  frame. Use an exploded-coordinate test (`abs(pos).max() > 1e4`) alongside peptide C-N.
- Only 79ALA showed NaN at scan time; the other three variants are equally over the limit and must not
  be assumed safe.
- Fix is a timestep reduction (0.0006 for NP-equivalent margin = 15x cost) and is **not yet applied** --
  it needs a cost decision, and it is worth first checking whether 2.64 A contacts mean the micelle is
  over-compressed. Full record, diagnostics and resume steps: `glpG_REMD_blowup_debug.md`.

## NP-1AO6 blow-up ROOT CAUSE: dt exceeds the MARTINI LJ-core stability limit (2026-08-08)

- Velocity-Verlet is stable only for `dt < 2/omega`. Taking `omega = sqrt(U''/mu)` from the **tabulated**
  grid curvature (not an analytic fit), the limit collapses as pairs approach:
  `dt_max = 0.0209 @ 4.0 A, 0.00808 @ 3.5 A, 0.00430 @ 3.2 A, 0.00278 @ 3.0 A, 0.00186 @ 2.84 A`.
- The closest protein-environment approach the system actually samples is **2.84 A** (1st pct 2.88,
  median 3.54; 40% of frames below 3.5 A). So **dt=0.005 was 2.7x over the limit**, and over it for any
  approach closer than ~3.25 A. Failure is therefore stochastic, which explains the 16x spread in onset
  across faces, why 5/6 survived while one tore, and why 0.009 -> 0.005 only slowed it.
- Proven by an A/B restart from the identical pre-tear frame (same seed, only dt differs):
  dt=0.005 -> 42 broken bonds and Epot -7988 -> +9130 by t=19.8, avg KE/1.5kT = 5.34;
  dt=0.001 -> 0 broken bonds, Epot -7988 -> -8146, avg KE/1.5kT = 1.006.
- **FIX: dt = 0.001** (`NP_DT` in `np_prod.sbatch`; `build_all.py` MD dict with step counts x5 to hold
  physical time). ~1,765 time units per 36 h block.
- Ruled out, each with evidence: PBC wrap (min-image distances identical to raw); bond-spring stiffness
  (`dt/dt_max = 0.035`, 29x margin -- the springs are NOT the constraint); anomalous close contacts
  (protein-ION/GOLD/MPA all steady ~4 A across the failure). The tearing bond was **68 A from the gold**,
  so this is not the "protein enters the gold LJ core" story assumed earlier.
- Do not diagnose the nonbonded term by summing the `pairs` dataset offline: that gives 1.58e14 against a
  logged -8020, because the engine does not evaluate intra-protein pairs (Upside's own FF handles those).
  The backbone O(i)-N(i+1) pairs sitting at 0.1-0.2 A are a **red herring** for the same reason.
- Underlying tension worth remembering: Upside protein atoms are mass 1.0 m_up vs MARTINI beads at 6.0,
  and `dt_max ~ sqrt(mu)`, so the light all-atom backbone is what forces the small timestep. Mass
  repartitioning would buy ~1.9x but changes kinetics -- not done.

## NP-1AO6 blow-up at dt=0.005 after findings-88 fix: backbone spring instability at large amplitude during unfolding (2026-08-11)

Job 53234804 ran with `NP_DT=0.005` after findings 88 removed O sites from the MARTINI pair table.
Runs 1 (90-0-0) and 4 (0-90-0) were destroyed at block 1, t≈250; runs 0,2,3,5 healthy.

**Root cause: backbone spring bonds at large oscillation amplitude during unfolding, not MARTINI LJ.**

Evidence and mechanism:
- Blow-up localised at residues 119–122 (MET-CYS-THR-ALA-PHE), which are 70+ Å from the NP surface.
  No MARTINI BB-env pair is close: nearest ion is 8.6+ Å throughout the failure window.
- At frames 920–935 (t=248–252), backbone bonds (CA-C, C-N, N-CA) oscillate at amplitudes 0.3–0.5 Å
  — 2–3× the thermal amplitude sqrt(kT/k)=sqrt(1.0/48)=0.14 Å — as protein unfolding transfers
  energy into backbone spring modes.
- The spring bond linear stability criterion (ω×dt = sqrt(48/0.5) × 0.005 = 0.049 << 2) is satisfied
  throughout. But accuracy at large amplitude requires much smaller dt: at amplitude 0.5 Å, the bond
  force reverses sign over a half-period of ≈0.32 t_u (1.2 frame intervals at 0.27 t_u per frame).
  With only 54 force evaluations per frame interval, the reversal is tracked too coarsely.
- At frame 935, a coincident alignment of spring forces (CA-C compressed + C-N stretched) produces
  ≈14.4 E_up/Å sustained force on C(121). Over 54 steps of dt=0.005, this injects 12 E_up of kinetic
  energy into C(121) (velocity jumps from 0.96 to 4.96 Å/t_u), stretching CA(121)-C(121) to 2.58 Å.
  A cascading wild oscillation follows; final blow-up at frame 944 (t=254.88).

**Key distinction from the Aug-2 blow-up (findings §96):**
- Aug-2: MARTINI LJ core instability at O-site approaches of 2.84 Å (ω×dt = 2.7× over the limit).
  Fix: remove O from pair table (findings 88).
- Aug-11: Backbone spring at large amplitude during unfolding. MARTINI is irrelevant. The linear
  stability criterion is met. The failure is from integration accuracy at large displacement.
- The old claim "the springs are NOT the constraint; dt/dt_max = 0.035" (findings §111) applies only
  to small-amplitude (thermal) oscillations. It does not hold when unfolding drives bonds to 3-7× thermal.

**Fix:** `NP_DT=0.001`. This gives 270 force evaluations per 0.27 t_u interval — sufficient to track
the force reversal as atoms move. Runs 1 and 4: `/output` (blown-up block 1) deleted; `/input/mom`
has `restart_valid=1` (set at block 1 start from block 0's last frame), so the driver restarts cleanly
from the block-0 endpoint. Resubmitted as job 53251911.

**Do not raise dt above 0.001 for NP.** Even though the MARTINI contact limit allows dt≈0.3, protein
unfolding creates a secondary constraint: backbone spring accuracy at large amplitude. The 50 t_u
validation test is too short to sample the full unfolding regime (failure occurs at t>250).

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

## Update 76 (2026-08-04): glpG-79HIS TM4 unfolding is a DDM-slab hydrophobic-mismatch artifact

The hybrid glpG-RKRK-79HIS DDM runs lose TM4 while the implicit-membrane runs keep it. The implicit result is
the correct one, and the hybrid failure is a preparation defect, not protein-force-field behaviour.

Construct mapping (verified against UniProt P09391 and the experimental HDX file, whose protein field reads
`GLPG_S201T_TM`): construct index + 66 = E. coli GlpG numbering. The base construct is already the catalytically
dead S201T (construct 135). `79HIS`/`79ALA` is WT H145 vs H145A; `S115T` is S181T. Helices from the input PDB:
TM1 25-47, L1 48-76, TM2 81-100, TM3 104-126, **TM4 134-153 (WT 200-219)**, TM5 160-176, TM6 184-208.

Measured on `glpG-RKRK-79HIS_local_T0.8647_prod.vtf` (401 frames, 210 residues, 425 DDM, 100.8 A box):
TM4 is the only helix that fails -- O(i)..N(i+4) helicity 0.90 in the first decile to 0.17 in the last, internal
RMSD 4.45 A, against 0.84-0.96 helicity and 1.6-2.5 A for TM1/TM2/TM3/TM5/TM6. TM4 is also the most buried
helix (22.6 other-helix CA within 12 A; 6-16 for the rest), consistent with the literature statement that five
helices surround a short TM4 that begins deep in the membrane at S201.

The failure runs through the TM4:TM6 interface -- the GxxxGxxxA/GxxxAxxG heterotypic packing whose disruption
(G261V) costs 28.1 C of Tm. CA-CA<10 A contacts fall 127 -> 50-60 by frame 60, with 0 -> 19-25 DDM beads having
three or more TM helices within 8 A. Detergent intercalation precedes the helix loss (interface eroded by frame
60, helicity collapse at ~frame 270).

Root cause: a DDM *lamellar* slab cannot solvate GlpG. The CHARMM-GUI template (`charmm-gui-8543403667`,
500 DDM in 100x100 A, 40 A^2/DDM) has a C1-C3 tail core of only 11 A; in the run it is 12.7-13.9 A against a
28-30 A protein TM belt, and 50% of TM backbone CA have a polar maltose bead as their nearest detergent bead.
Every residue that unwinds sits at or beyond the edge of that 14 A core (TM4 WT 208-219 at z = -5 to -21 versus
a core of -4.6..+9.3; also construct 50-57, 96-98, 176, 185); no residue inside the core unwinds. This is not
tunable: lamellar thickness = 2*V_tail/APL, and V(C12) ~ 324 A^3 means a 28 A core needs APL ~ 23 A^2, which a
maltose head (>= 40 A^2) cannot reach. DDM forms micelles for exactly this reason. Note also that
`--membrane-thickness-angstrom` does not set slab geometry -- it is only read for ion counting
(`py/martini_prepare_system.py:477-480`), so the 48.8/52 A values passed had no geometric effect.

Second, smaller defect in the diagnostic run: production began before the insertion void healed. Whole-DDM
deletion at the 5.276 A CG-envelope clearance emptied the first solvation shell -- 14 DDM beads in the 3-6 A
shell at frame 0 versus 213 at the end, mean nearest-DDM over TM CA 11.7 A -> 5.6 A -- and with min 500 /
eq 500x6 / burn-in 0 that shell refilled during the analysed segment, partly through the protein interior.

Contrast with the earlier `glpG-RKRK-79HIS_DDM_REMD_T0.70` system, which is a *micelle*: 235 DDM in a 186.8 A
box, xy r95 = 45 A against box/2 = 93 A, tail core 31.5 A matched to the 28 A belt. There TM4 only degrades
0.85 -> 0.69. Not a controlled comparison (T = 0.70 versus 0.8647), but the direction is unambiguous.

Experiment settles which result is right. In `glpG-RKRK-79HIS_rawuptake_HXMS.csv` (pD 9) the TM4 peptide 140-144
(ALMGY = WT 206-210) is 4.6% deuterated at 12 min, 17% at 1 h, 50% at 24 h, 97% at 20 d, i.e. dG_op ~ 10-11
kcal/mol -- matching the ~10-13 kcal/mol the implicit-membrane run reports there. TM3 109-118 shows essentially
no exchange even after 20 days, which is why it goes off-scale implicitly. A trajectory that unwinds TM4 within
a sub-microsecond segment is inconsistent with both. The hybrid REMD dG profile is globally 0-7 kcal/mol against
10-15+ experimentally, so every hybrid glpG HDX number inherits this defect.

Lesson: before trusting any hybrid membrane-protein run, gate the preparation on two geometric measurements
rather than on the requested parameters -- (1) the environment's tail-bead hydrophobic thickness measured *at the
protein* against the protein's own TM belt span, and (2) first-shell occupancy (nearest environment bead per TM
CA ~ 4.5-5.5 A) at the start of the analysed segment. A detergent whose packing parameter favours micelles must
be built as a micelle; asking CHARMM-GUI's Bilayer Builder for it yields a slab that is mostly headgroup.

## Update 77 (2026-08-04): DDM rebuilt as a micelle -- what the construction actually requires

Implementation of the findings-76 fix. `derive_environment_morphology` counts acyl chains as connected
components of the apolar (`C1`-`C5`) bond subgraph in the lipid ITP: one tail -> micelle, two or more ->
bilayer. DDM -> micelle, DOPC/POPC/POPG -> bilayer, detergent+lipid mixtures rejected. Deriving it from
topology rather than a name table or a CLI flag makes the unphysical combination unreachable, and matches the
module's existing rule that no lipid identity is hardcoded. Note `parse_lipid_from_itp` already returns
**0-based** bead indices; subtracting 1 again silently split DDM's single tail into two and reported it as a
bilayer.

`build_detergent_micelle` wraps a monolayer on the belt (tails in, heads out), belt span from the OPM
reference's own thickness REMARK. Three construction facts, each found by measurement:

* **Seed from the shell VOLUME, not a surface.** Seeding from convex-hull support points gave 32 molecules
  (653 A^2 per molecule against a 40 A^2 template) with 25 bald belt beads -- the hull has few vertices, so
  one monomer per vertex is a sparse shell. Filling the shell volume gives 186.
* **Fill innermost outward.** Random-order tail-tip seeding let molecules seeded far out block the contact
  layer: tails smeared over 63 A and coverage got *worse* despite 304 molecules. Sorting seeds by distance
  to the belt packs tails against the protein and pushes heads outward -- the actual PDC organization.
* **The template is not a density reference.** The CHARMM-GUI step5 assembly is PRE-minimization and contains
  bead pairs 0.24 A apart. Measuring an intermolecular contact distance from it produced a 0.237 A packing
  threshold, which accepted 861 mutually overlapping molecules. The force field's own contact distance
  (2^(1/6) sigma_max = 5.276 A for DDM) is the correct spacing, and the 6.0-6.6 ladder condenses from there.

Exposure must be defined the same way the coverage is enforced, or the gate is unsatisfiable. "A free cell
within contact range" wrongly marked *buried* beads exposed (the 79 flagged beads had burial 68-105 protein
neighbours within 12 A against a median of 67, i.e. core beads like TM4's). The correct test is that some
free outside cell at contact range has THAT bead as its own nearest protein bead -- anything shielded never
wins the nearest lookup. That cut 317 belt beads to 106 genuinely exposed, of which 105 are tail-covered.

**A packed-state thickness span cannot be gated on.** The first gate compared the environment's 5-95
percentile tail-bead z span with OPM's hydrophobic thickness. It fails DOPC (20.4 A), POPE/POPG (20.6 A) and
DDM (11.0 A) alike against a 22.9 A limit, because CHARMM-GUI templates are laterally compressed and a
clipped percentile is not the same quantity as a relaxed hydrophobic thickness. Gating on it would have
broken the validated lamellar path. What ships instead is `assert_environment_solvation` at the production
handoff, on equilibrated coordinates and per belt residue: hard-fail on vacuum (any belt site with no
environment bead within 2x the contact distance), and REPORT acyl-tail reach plus the local tail-core
thickness. The report is deliberately not a gate, because on a post-damage snapshot both recover -- detergent
floods the cavity a stripped helix leaves -- so they diagnose rather than decide. The contrast is still
unambiguous: micelle seed 45.1 A local tail core, 0 sites beyond tail reach; old lamellar replica 13.9 A,
2 sites; belt 28.2 A.

Validated 79HIS micelle: 186 DDM, aggregate r95 43 A in a 137 A cube (compact, not spread), tail-z span 45 A
over the 28.2 A belt, first-shell distance tightening 7.04 -> 5.76 A, all six TM helices 0.97-1.00 helicity
including TM4, box fixed throughout (barostat off, as a finite aggregate requires).

Lesson: when a build-time gate and the thing it guards are measured differently, the gate is worthless --
either measure both on the same state (equilibrated, local, per residue) or demote it to a report. And check
a new geometric criterion against the paths it must NOT break before shipping it, not only against the case
that motivated it.

The DOPC bilayer regression confirms the lamellar path survived the refactor: morphology resolves to
`bilayer`, the xy-only tensionless barostat is still live (box 200.2 -> 213.3 A laterally with z pinned at
123.697), tile/carve insertion runs, ions come out at 287 pairs -> 287 Na / 293 Cl neutralising the +6
protein, and the solvation gate passes with 0 bare sites, 0 beyond acyl-tail reach and a 23.0 A local tail
core against the 28.2 A belt -- correctly NOT warned, since that is ordinary lipid/protein adaptation rather
than mismatch. Ion parity was also asserted numerically: for the recorded slab inputs the refactored
excluded-volume form returns the same 65 pairs as the old `A_xy*(L_z - t) - V_protein` expression.

## Update 78 (2026-08-05): NP-1AO6 job 53032366 -- the tear is NOT a gold-core collision

Job 53032366 ended its own chain after its health check found `np.run.1` at peptide C-N 5.29 A (limit 2.5).
Diagnosed from the trajectory; the mechanism is different from the documented one.

* **Single-step onset.** For 2559 of 2974 frames every face was thermal (max C-N 1.67-1.88 A, zero bonds
  over 2.5 A). Between frame 2559 and 2560 -- one recorded frame, ~38 steps -- face 1 jumped to 14.29 A at
  res550. Not gradual stretching: res551's N moved 10.21 A in that interval while the median atom moved
  0.35 A (99th percentile 1.07 A), a ~30x outlier confined to res550-551 backbone plus the res550 `BB`
  virtual site.
* **Not an environment collision.** With the protein mask applied, res550/551 sit **~40 A from the nearest
  GOLD bead** and 36-38 A from MPA at onset; nearest ION 5.2-5.7 A; nearest non-local protein atom 2.4-2.8 A
  before and 3.5 A after. Nothing was in contact.
* **Gold contact elsewhere was stable all run.** Minimum protein-GOLD distance stayed 3.0-4.3 A from frame
  400 to the end (residues 76-92 and 467-494 in persistent contact) with zero bond breakage. So the earlier
  hypothesis -- backbone tears when a segment creeps into the gold LJ core -- does NOT explain this run.
* **Propagation.** By the last frame 22 bonds exceed 2.5 A, scattered in 1-3 residue fragments over
  res 526-575 (HSA domain III). CA(i)-CA(i+1) there averages 5.96 A against 4.00 A elsewhere and 3.8 A
  native, so the local backbone geometry is genuinely distorted, not just the C-N term. Rg stays
  27.6 -> 31.1 A and the other five faces ended at 1.73-1.83 A.

This is a different failure from the Aug-2 dt=0.009 blow-up, which was global (~350 bonds at 12-17 A, C-N
minima collapsed to 0.2-0.6 A, reaching NaN). Here one localized impulse appears in an otherwise healthy
system with no proximate partner, which points at the bonded/interface force path rather than a steric
event -- the `BB` virtual site moving with the break is the reason to look there first, since env force is
evaluated at that site and redistributed to N/CA/C. Unconfirmed from trajectory data alone; the next step is
to instrument per-term forces on res549-552 around the onset, not to lower dt again.

Two analysis traps hit while doing this, both worth remembering. Selecting `residue_ids == r` without the
`particle_class == "PROTEIN"` mask picks up GOLD/MPA/ION beads that share residue numbering, which reported a
protein-gold distance of exactly 0.00 A (a gold bead measured against itself). And taking "distance to gold"
at the *argmax-C-N* residue per frame is meaningless while the argmax wanders -- it gave 36-61 A for the
wrong reason. Pin the residue first, then track it.

Separately: all six faces carry **218 K+ and zero Cl-**, which matches the CURRENT builder -- `np_hybrid.py`
asserts counterions-only and 203 MPA(-1) + protein(-15) = -218 exactly. The `plan.md` note of "2053 pairs,
K+ 2053 / Cl- 1835, 3888 ions" describes a SUPERSEDED convention. The running systems are therefore not a
stale build, but the two conventions differ physically: counterions-only means no bulk ionic strength at all,
whereas the earlier intent was 0.15 M. Which is wanted is a modelling decision for the user, not a bug.

## Update 79 (2026-08-05): glpG micelle REMD -- one runaway replica takes the whole ladder down

The four micelle REMD jobs started 2026-08-05T00:43. After ~11.5 h, **two of four are dead** and two are
healthy:

| variant | state | first NaN |
|---|---|---|
| glpG-RKRK-79HIS | NaN, all 48 replicas | step 12190 |
| glpG-RKRK-79HIS_S115T | NaN, all 48 replicas | step 6624 |
| glpG-RKRK-79ALA | healthy at step 9338 | -- |
| glpG-RKRK-79ALA_S115T | healthy at step 10755 | -- |

**The mechanism is single-replica runaway plus REMD propagation, not a systematic instability.** At step
12144 in 79HIS, 47 of 48 replicas were perfectly healthy (hbonds 178-212, Rg 19.2-19.9, protein_potential
-2168 to -2549, total -8447 to -9440). At the same step `system 16` had **Rg 935,837 A and
protein_potential 9.6e15**, and `system 17` was already positive at +7581. One replica diverged; the next
exchange sweep carried a non-finite energy into the Metropolis test and poisoned all 48 within one step.
Temperature is not the trigger: systems 16/17 sit at T=0.77, mid-ladder, not at the 0.90 top.

Transient positive `protein_potential` excursions are common and normally recover -- the two surviving runs
peaked at +36,450 (79ALA, step 13064) and +41,012 (79ALA_S115T, step 720) with hbonds 195-202 and Rg
18.5-20.0 throughout, i.e. structurally fine. So the rare event is not the spike itself but a spike that
fails to recover. With 48 replicas x 4 variants x many chunks the chance of at least one runaway is high,
which is why half the batch died while the other half looks pristine.

Two separable defects, and the second is the one that turned a 1-replica problem into a 48-replica loss:

1. **Numerics:** a rare non-recovering force spike at dt=0.009. Cause not yet localized. The 400-step
   single-replica smoke test on the seed could never have caught this -- it needs ~10^4 steps across 48
   replicas, so seed smoke tests bound loading and setup, NOT stability.
2. **`run_remd.py` has no finiteness guard.** It exchanges on whatever energy it reads, so one NaN replica
   destroys the ladder, and because the driver reseeds `/input` from the last output frame the corrupted
   state is then written back and the chain self-resubmits onto it. This is the same "isfinite is not a
   health check" lesson already learned on the NP track (`health()` there checks peptide C-N), except the
   REMD driver never got the equivalent. It should test every replica for finiteness (and a sane Rg /
   sign-of-potential) before each exchange sweep, and drop or halt a bad replica rather than swap with it.

Both NaN jobs were stopped (STOP file then scancel 53036661, 53036664) so they could not keep resubmitting
onto corrupted state; 53036667 and 53036669 were left running, since their data is clean -- but they carry
the same exposure until the guard exists.

## Update 80 (2026-08-05): the spline table was NOT dry-MARTINI -- bare cutoff instead of reaction field

Hunting the glpG runaway (findings 79) by elimination. Everything below was checked against the prepared
79HIS micelle seed and the engine source, and is CLEAN:

* **Table indexing** -- `coefficient_indices` span 0-21 over a 22-row table, all rows used, none out of range.
* **Intra-protein exclusion** -- `skip_pair_if_intra_protein` does fire. It has to: the pair list contains
  550,725 PROTEIN-PROTEIN entries including 390 below 1 A, because the `BB` virtual site is coincident with
  the backbone atoms it is built from (closest 0.704 A, a backbone C against its own BB proxy).
* **Contact distances at production start** -- every CROSS-class minimum is healthy: ION-ION 5.64,
  ION-LIPID 4.12, ION-PROTEIN 5.01, LIPID-LIPID 4.11, LIPID-PROTEIN 3.26 A, and zero pairs below 3 A. The
  physically active protein-env set is BB-vs-env only (protein-env pairs are skipped unless the protein
  side is a derived BB proxy), and its minimum is 4.35 A with 3 pairs under 4.5 A.
* **Masses** -- env 6.0 m_up (72 Da, a MARTINI bead) and protein 1.0 m_up (12 Da, one heavy atom). Correct.
* **Pair list** -- a textbook Verlet scheme: built to `cutoff + cache_buffer` = 14 A, rebuilt when any atom
  moves more than `0.5*cache_buffer` = 1.0 A. Two atoms closing 2 A between rebuilds cannot cross the
  12 A cutoff undetected, so no tunnelling.
* **Spline evaluation** -- `LayeredClampedSpline1D::evaluate_value_and_deriv` clamps both ends (`x>=nx-1`
  and `x<=0` return the clamped value with zero derivative) and `x_bin` indexes `nx-1` segments safely.
  No out-of-bounds read.
* **Spline ringing** -- despite the grid spanning 12 orders of magnitude at 0.012 A spacing, only 3 of 81
  rows show a wrong-sign derivative inside the repulsive core and only at ~1e-2 dE/dx; overshoot at the
  minimum is 2.6e-5 E_up. Not a force-sign error.
* **Thermostat** -- and a documentation correction: `--thermostat-interval -1` does NOT mean NVE.
  `max(1., round(-1/(inner_step*dt)))` evaluates to **1**, so the thermostat fires EVERY step. The note in
  `progress.md` describing it as "effectively NVE" was wrong.

**The one real defect: the spline table was not the dry-MARTINI potential.** It stored bare LJ plus bare
`1/r` Coulomb, hard-truncated at 1.2 nm. Published dry Martini is run with `coulombtype = reaction-field`,
`epsilon_r = 15`, `epsilon_rf = 0` (infinite -> conducting) and `vdw-modifier = Potential-shift`, so both
terms reach the cutoff smoothly. The stored table therefore had a step at r_c of `k*qq/r_c` = **2.65 E_up
for a charged pair, about 3.4 kT** -- verified analytically (LJ -0.057 + Coulomb -2.648 = -2.705, matching
the stored -2.7049 exactly). This is the correction flagged as outstanding back in findings 75.

Fixed in `py/martini_build_tables.py`, in BOTH nonbonded builders (the particle-particle grids and the
SC-env/BB-env `_pair_energy_and_grads`, which had the same bare form including its analytic gradient).
`scratchpad/verify_table_matches_drymartini.py` now asserts equivalence by rebuilding the reference in
native kJ/mol + nm and converting:

| table | max deviation from the reference form | rows non-zero at the cutoff |
|---|---|---|
| old (`martini.h5.bak.pre-reactionfield`) | **3.95 E_up** | 81 / 81 |
| regenerated | 2.2e-11 E_up (round-off) | **0 / 81** |

Note the neutral pairs change by a constant only (forces bit-identical); charged pairs change by 3.1-3.9
E_up across the sampled range, because the reaction field is a genuine r-dependent term. So this alters
results for anything with charges -- ions, PC/PG headgroups, charged residues. DDM itself is neutral.

**Causation is NOT yet established, and one attempt to establish it failed.** A near-NVE drift comparison
was inconclusive: both tables drifted DOWN by a similar amount (-332 vs -418 E_up over 20,000 steps) with
identical kinetic ratios (1.033 vs 1.036). That measurement was dominated by ongoing structural relaxation,
not by conservation error, and 20 frames is far too coarse to separate them. Reasoning it through: a step in
ENERGY at the cutoff breaks conservation but applies no large FORCE, so it should cause gradual heating
rather than the single ~1e15 kick actually observed. An A/B reproduction (12 seeds per arm, 16,000 steps at
T=0.77, the temperature and step range where both losses occurred) is the test that decides it.

**A/B result: NULL, and underpowered by design.** 12 seeds x 16,000 steps at T=0.77 per arm: 0 diverged,
0 NaN, worst Rg 20.6 A (old) and 20.8 A (new) against a native ~19-20 A. So it distinguishes nothing.
The arithmetic I should have done first: the observed failure rate is ~2 events across 4 variants x 48
replicas x ~1e4 steps, i.e. ~2 per 2e6 replica-steps. Each arm here was 1.9e5 replica-steps -- about a 10%
chance of a single event even in the defective arm. A powered local A/B would need ~48 runs x 16,000 steps
per arm just to expect ~1 event, which is a poor way to spend the compute.

**Where that leaves the fix.** The table correction stands on its own: the user's requirement is that the
spline table BE the original dry-MARTINI potential after unit conversion, and it demonstrably was not
(3.95 E_up deviation, 81/81 rows discontinuous at the cutoff). That is fixed and verified independently of
the NaN question. But it must NOT be reported as the proven cause of the runaway -- it is an unproven
contributor, and the mechanism argues against it being the direct trigger (an energy step applies no large
force).

The efficient way to settle causation is the natural experiment we owe anyway: re-prep and re-run all four
glpG variants on the corrected table and compare the failure rate against the observed 2/4. If a replica
still diverges, capture it with dense output around the onset rather than trying to reproduce a 1-in-1e6
event locally.

Lesson: compute the expected event count BEFORE running an A/B on a rare stochastic failure. Two tests were
spent here without discriminating power -- one confounded by structural relaxation, one 10x under-exposed.

## Update 81 (2026-08-05): guards removed from the engine; NO GUARDS is now a standing rule

User directive: no guard, anywhere, for any numerical problem. Recorded as a **NO GUARDS** rule in
`AGENTS.md` alongside a companion rule that a spline table must BE the published potential, not a variant.

Factual note for the record: the `src/main.cpp` blow-up guard was **not** added in this session -- no C++
file was touched while diagnosing findings 79-80 (`git status src/` was clean, so the guard was already
committed on `martini-dev`). Master has neither it nor `martini_potential.cpp`.

Removed, and the engine rebuilt clean:
* `src/main.cpp` -- the blow-up guard that aborted the run on a non-finite potential or kinetic energy, plus
  its `blew_up` / `blew_up_ns` / `blew_up_round` state, the loop break, and the FATAL block.
  `compute_logged_kinetic_energy` stays: it is still used for ordinary logging.
* `src/martini_potential.cpp` -- three silent-skip masks, which were the harmful ones. The pair loop's
  `if(isfinite(pot) && isfinite(force_mag))` dropped non-finite pair contributions outright; the main.cpp
  comment even admitted this by noting its kinetic check existed to catch "diverging momenta even when
  non-finite pair forces are masked out of the potential". The SC-env path had the same pattern twice
  (`!isfinite(value) || !isfinite(dVdr) || !isfinite(dVdcoord)) continue`).
* Kept in `update_martini_node_boxes`: `!(scale_xy > 0.f) || !(scale_z > 0.f)`. A box length being positive
  is a domain precondition, not a masked numerical error. The `isfinite` clause beside it was removed.

**Operational consequence, stated plainly:** a divergence now propagates instead of being stopped. NaN will
enter the log, be exchanged between replicas, and be written into restarts. That is the intended trade --
the evidence survives instead of the run -- but it means a bad run must be caught by monitoring, and the
`glpG` 11-hour NaN propagation becomes expected behaviour rather than a stale-binary symptom (which makes
the question of whether the cluster binary carried the guard moot).

**Three existing checks are arguably in scope and were left alone pending a ruling**, because none masks a
numerical error: `assert_environment_solvation` (prep-time, fails a build whose belt faces vacuum), the
`np_hybrid.py` ion assertions (neutrality and salt composition of a built system), and `run_np_prod.py`'s
`health()` peptide C-N check, which ends a chain rather than propagate a torn system. The last of these is
the closest to a guard under the new rule and is the one most worth a decision.

## Update 82 (2026-08-05): CLC-ec1 (Robertson eLife 2021) -- blocked on a POPE/POPG bead-model mismatch

Scoped as equilibrium monomer + dimer with a lipid-solvation comparison (their dimerization dG needs
free-energy machinery this pipeline does not have, and metadynamics was already rejected on glpG).

What works, and is worth keeping:
* Their `step8_production.mdp` independently confirms Update 80: `coulombtype = reaction-field`,
  `epsilon_r = 15`, `epsilon_rf = 0`, `vdw-modifier = Potential-shift-verlet`, `rvdw = rcoulomb = 1.2`.
  Their `ref_t = 303.15 K` is exactly our 0.8647 T_up.
* Protein input: their deposited PDBs are coarse-grained, so 1OTS is needed instead. 1OTS is CLC-ec1 **plus
  a Fab complex** -- chains A/B are the transporter, C-F must be stripped. The OPM-oriented copy carries
  `1/2 of bilayer thickness: 14.9` (a 29.8 A belt) and contains only A/B. Monomer = chain A (444 residues),
  dimer = A+B (885). Multi-chain is supported (`derive_chain_break_metadata`), and OPM orientation came out
  perfect (444/444 and 885/885 Ca within 4 A, since the input already came from the OPM frame).
* Box sizing is solvable. Their 292 A box would be ~37,700 particles / 7e8 all-pairs entries (~12 GB), which
  the all-pairs list cannot carry. The prep floors the box at one template tile, so the tile has to shrink:
  cropping their bilayer to a 150 A tile from its densest region gave 719 lipids at 62.6 A^2/leaflet against
  their 61.7 bulk and a composition of exactly 2.00:1 -- i.e. representative, not the protein-shaped hole.
  That lands the dimer at 13,053 particles and the monomer at 10,848, both tractable.
* Hydrophobic match is good: packed tail span 28.2 A against the 29.8 A OPM belt.

**The blocker: POPE/POPG bead models disagree, so POPE/POPG has never actually run through this pipeline.**

| source | POPE beads | tails |
|---|---|---|
| Robertson `last_frame.pdb` | 12: NH3 PO4 GL1 GL2 C1A **D2A** C3A C4A / C1B C2B C3B C4B | 4 + 4, unsaturation in tail A |
| CHARMM-GUI Martini Maker | 12: identical to the above | 4 + 4 |
| our `dry_martini_v2.1_lipids.itp` | 13: NH3 PO4 GL1 GL2 C1A C2A C3A C4A / C1B C2B **D3B** C4B C5B | 4 + 5, unsaturation in tail B |

Different bead count, different tail lengths, different unsaturation position -- not a renaming. DOPC matches
exactly between CHARMM-GUI and our itp (14 beads, same names), which is why DOPC and DDM work today and
POPE/POPG silently never did. A 2:1 mixture built by substituting our own POPG headgroup (`GL0`, P4, q=0)
for POPE's (`NH3`, Qd, q=+1) is exact -- our POPE and POPG differ in atom 1 only -- but it does not help,
because the CHARMM-GUI template it starts from is already the 12-bead model.

Options, in preference order: (1) write a 13-bead POPE/POPG bilayer builder matching our itp, in the style of
`scratchpad/build_ddm_slab.py` for DDM -- the itp is the authority and the coordinates must be generated for
it; reusable for any dry-MARTINI lipid. (2) Add the 12-bead POPE/POPG to the itp -- rejected: it means
running a lipid the force field was not parameterised for, against the rule that a table must be the
published form. (3) Substitute DOPC -- rejected: the paper's result is specifically about POPE/POPG
composition, so the membrane cannot be swapped.

Lesson: a lipid being present in the ITP does not mean the pipeline can build it. Check that the bead names
in the membrane template and the ITP agree BEFORE sizing boxes and fetching structures -- one `diff` of the
two name lists would have found this in the first minute rather than after the orientation, cropping and
composition work.

## Update 83 (2026-08-05): POPE/POPG builder -- rigid-rod conformers cannot start at the target area

Wrote `scratchpad/build_popepopg_bilayer.py` to generate the 13-bead POPE/POPG bilayer our ITP needs
(Update 82). What is correct: geometry is read from the ITP and bond lengths are exact (worst deviation
0.0011 A), composition is exactly 2.00:1, the arrangement is a proper bilayer (tails inside heads), and the
tile loads back through `read_charmm_gui_membrane` with its `[ molecules ]` cross-check passing.

What is wrong: the starting packing is ~2.2x too tight. Measured against the paper's own equilibrated
bilayer, per-molecule nearest intermolecular bead distance (minimum-image in xy):

| | min | median | p90 | APL |
|---|---|---|---|---|
| Robertson equilibrated, 12-bead | 4.21 | 4.60 | 4.71 A | 61.7 A^2 |
| this build, pre-equilibration | 0.20 | 2.05 | 2.46 A | 69.4 A^2 |

Their median of 4.60 A is essentially the MARTINI sigma of 4.7 A, achieved at a TIGHTER area than ours. The
cause is structural, not a tunable: at 69.4 A^2 the lipid axes are 8.33 A apart while a rigid conformer with
both tails at fixed lateral offsets is ~6 A wide at every height, leaving ~2.3 A -- exactly the median
observed. A real lipid spreads its beads along z with tilt and disorder, so neighbours interleave rather than
confront at the same height. Four attempts to fix it inside the rigid-rod picture all failed: reducing the
tail splay 12 -> 5 deg and GL2 20 -> 60 deg (no change), first-bond-only splay plus jitter-only spin
(min 0.39 -> 1.72, median unchanged), z-stagger (made min WORSE, 1.61 -> 0.20), and wider spacing.

Correct route, and it doubles as the deliverable the user actually asked for: start deliberately LOOSE and
condense with a pure-bilayer run under the xy barostat -- the same two-step the DDM slab template used --
then use the equilibrated tile as the template for the protein systems. That equilibration IS the lipid
validation: compare the condensed APL, bilayer thickness and this packing statistic against Robertson's
61.7 A^2 / 4.60 A. Note the loose-start trick only works for a PURE bilayer, where the box is free to shrink;
in a protein system the box is pinned by the protein footprint, so the lipid count must already suit the
area, which is why the template has to be condensed first and measured second.

Lesson: I had the reference bilayer in hand the whole time and spent four iterations tuning against
intuition instead of measuring the target statistic on it first. When reproducing someone's system, compute
the comparison metric on THEIR data before touching your own builder -- it converts guesswork into a
one-line pass/fail.

## Update 85 (2026-08-05): dry-MARTINI membranes are NVT -- a barostat cannot find the area without solvent

User's objection, and it is right: with holes expected in a bilayer, a tensionless barostat is the wrong
ensemble. Three independent reasons converge, and the second is the one that settles it generally.

1. **A defect is a stress sink.** Lateral tension relaxes through a pore or an under-filled patch instead of
   through the lipid area, so zero measured tension stops marking the intact bilayer's equilibrium; the
   barostat compresses the intact regions past it. CLC-ec1 dimerization is *about* lipid-depleted states.
2. **Implicit solvent has no solvent virial.** Dry MARTINI carries no water, so the lateral pressure the
   barostat reads is missing the solvent contribution entirely and is not the physical one. Dry Martini
   (Arnarez et al., JCTC 2015) is specified NVT for exactly this reason.
3. Robertson et al. used semiisotropic Parrinello-Rahman legitimately because **their system is wet** --
   `tc-grps = membrane water_and_ions Protein`, `energygrps ... water_and_ions` in `step8_production.mdp`.
   The water is merely stripped from the frames they deposited (`compressed-x-grps = membrane Protein`),
   which is why the frames look solvent-free and invited the wrong inference.

Measured consequence, our own tile run both ways: under the tensionless xy barostat it condensed to
**APL 56.0 A^2 against their 61.7 (-9.4%)**, with tail core +1.6% and head-head +0.9% -- over-compressed and
correspondingly thicker/tighter, which is the signature. Not a pore artifact in that run (0% empty bins,
continuous midplane density), just the barostat in a small periodic tile.

So the area is now an **input matched to theirs**, taken from an equilibrated simulation of the same lipids at
the same temperature, and the validation moves to the area-sensitive structural observables measured AT that
area. Their `step8_production.mdp` also confirmed two things we already had right: **T = 303.15 K** (= 0.8647
T_up, our value) and `reaction-field / eps_r 15 / eps_rf 0 / Potential-shift / rvdw = rcoulomb = 1.2`, i.e.
independent confirmation of the corrected spline tables (Update 80).

Their deposited composition is 1814 POPE : 959 POPG = **1.892:1**, not the nominal 2:1 they describe -- worth
knowing before treating an exact 2.00 ratio as a fidelity requirement.

Lesson: before adopting an ensemble from a paper you are reproducing, check whether their system has the
degrees of freedom that make it valid. "They ran NPT" does not transfer to an implicit-solvent model, and the
frames they deposit may have had the solvent stripped, hiding the difference.

## Update 86 (2026-08-05): the bilayer prep squared up every tile, silently

`prepare_bilayer_structure` parsed the template's CRYST1 box into `bilayer_box`, used it only for tiling, and
then sized the actual box from the **lipid coordinate extent** with `force_square_xy=True`. For a rectangular
periodic tile that squares the box up to the longer edge: my 84.41 x 73.10 A tile became **84.26 x 84.26**,
a 15% stretch along y that both opens a vacuum stripe in the membrane and changes the area per lipid from the
61.70 A^2 it was built at to 71.0. It also loses the sub-bead gap between the outermost bead and the tile edge
even for square tiles (84.41 -> 84.26).

Fixed by using the CRYST1 record as the lateral box (`force_xy_box=target_xy`, `force_square_xy=False`) and
by letting `force_xy_box` be an (x, y) pair rather than a single edge. A periodic tile's box is the tile; it
is not recoverable from the coordinates, so a template without CRYST1 is now a hard error.

This is the same bug class as the earlier "box inflated to 83.3 A for a 75 A tile", which I had patched by
wrapping molecules into the tile so that extent ~= box. Wrapping made the symptom vanish for square tiles
while leaving the cause -- sizing a periodic box from a coordinate extent -- in place.

Lesson: when a fix works by making a wrong input look right (wrapping so the extent matches the box), the
cause is still there and will resurface on the first input that breaks the coincidence. Fix the derivation,
not the input.

## Update 88 (2026-08-10): the hybrid interface discards most of its own force -- Update 87's verdict is withdrawn

Update 87 concluded the glpG blow-ups were a timestep failure at protein-lipid contacts. **That verdict is
wrong.** It rested on omega*dt computed for the closest protein-lipid contacts, and every one of those
contacts turned out to be an O site -- a site whose force the engine throws away. The stiffness I measured
was the stiffness of pairs that exert no force at all.

`HybridPositionNode` (src/martini_hybrid.cpp:146) rebuilds both the O and the BB site from the three
carrier atoms (N, CA, C) every step, so neither is an independent degree of freedom. `propagate_deriv`
(line 189) then marks both as `derived`, which excludes them from the pass-through loop, and redistributes
only `bb_idx`'s gradient to the carriers. **The O site's sensitivity is neither passed through nor
redistributed: every O-environment force is silently discarded.** Separately, `map_backbone_sites` builds
BB from *four* weighted components (N, CA, C, O; `d < 4`, weights summing to 1) while `propagate_deriv`
redistributes over only *three* (`carrier < 3`), so BB's O weight, 16/54 = 0.2963, is dropped as well.

Measured at frame 67 of the 79HIS_S115T chunk, summing the analytic martini pair force per site:

| site | n | sum abs F (E_up/A) | max abs F |
|---|---|---|---|
| O (discarded) | 210 | 12537 | 2600 |
| N/CA/C (transmitted) | 630 | 5909 | 363 |
| BB (29.6% under-transmitted) | 210 | 482 | 15 |

More force is thrown away than delivered: ratio 2.12. The full pair force sums to 1e-12 as it must, but
with the O contribution removed the residual is **3590 E_up/A of net force**, i.e. spurious momentum
injected every step. The largest single discarded force is 2600 E_up/A (7577 kJ/mol/A) on THR115's O,
seven times the largest force the protein actually feels.

**Correction (same day):** I cited the engine's `--potential-deriv-agreement` (6.2% on the glpG chunk) as
corroboration. It is not evidence of anything. The metric is
`sqrt(sum (fd-analytic)^2 / sum fd^2)` with an FD step of 1e-3 A against a float-precision total potential
of order 1e4 E_up, so the per-component FD signal (~1e-3*F) is comparable to float round-off on the total
(~2e-3), and the number is round-off dominated at this system size. Measured on the *same* 1rkl file, the
pre-fix and post-fix binaries give **0.41960 and 0.41959** -- identical. Do not use this flag as a gate on
a system this large; it is a developer probe and its own help says so. The dropped-O-force defect stands on
the code path and on the analytically summed pair forces above, not on this number.

This unifies two things Update 87 left unexplained. Non-conservative forces do net work regardless of
step size, which is exactly why the `avg_kinetic_energy/1.5kT` excess measured 2.43 / 2.16 / 2.18 % across
a 4x dt range instead of falling ~4x per halving. And energy pumped into specific degrees of freedom is a
far better candidate for a blow-up localised on one residue's backbone than a marginal stability limit:
the O of a peptide unit absorbs a force its own neighbours never feel, so the local force balance is broken
and the springs tear.

A third finding, which needs a decision on intent rather than a bug fix. All five backbone sites
(N, CA, C, O, BB) appear in the martini pair list with **identical env pair counts (442260 each)**, the
same MARTINI type per residue, and full-strength epsilon (1.201 E_up = 3.50 kJ/mol, the same values the
lipid-lipid pairs use -- not divided down). Dry-MARTINI represents a residue backbone as ONE BB bead, and
N, CA, C, O are precisely the atoms that bead stands for, so representing them both as individual
sigma = 4.7 A beads and as a BB bead counts the backbone-environment interaction five times over. Either
only BB should be in the pair list, or the per-site epsilon should be partitioned; which one is a model
question. `_write_nonbonded_pairs` (py/martini_prepare_system_lib.py:2077) builds all unique i<j pairs with
only nrexcl=1 and ITP exclusions, so nothing ever restricted the protein side to BB.

In ENERGY the over-count is much milder than five-fold, because the four atom sites cancel: at frame 67
N/CA/C sum to -2295.7 E_up while O contributes **+2117.5** (its geometric placement lands it inside env
repulsive cores), giving a backbone-env total of -1368.9 against -1190.6 for BB alone -- a **1.15x**
over-count, a -182.7 E_up (-533 kJ/mol) shift if the four sites are removed. An earlier draft of this entry
called it fivefold; that is wrong for the energy.

The FORCE damage does not cancel, which is where the harm is: 12537 E_up/A discarded on O, 5909 E_up/A of
spurious-but-delivered force on N/CA/C, and 3590 E_up/A of net one-sided force per step. Note its
direction -- the env partner's sensitivity IS passed through while the protein's is dropped, so the
one-sided force acts ON THE ENVIRONMENT. That matches the blow-up's first measurable symptom, which was
`lipid_kinetic` reaching 2.487e16 at frame 68 while the protein kinetic energy was merely NaN.

Consequence to face squarely: every glpG-DDM and NP-1AO6 run made with this code has a partially
disconnected protein-environment interface and a non-conservative force. CLAUDE.md requires that BB-env
never be turned off, and roughly two-thirds of it has been off. Results from those runs should not be
interpreted until the gradient is fixed, and the mass and timestep questions cannot be answered
meaningfully before then, because fixing the gradient changes every force in the system.

Lesson: when the stiffest interaction in a diagnosis is also the most suspicious one, check that it is
connected to the dynamics at all before building a theory on its magnitude. I measured omega*dt for 210
sites for two rounds without asking whether their forces reached the protein -- the answer was that they
did not, and the real defect was one level below where I was looking.

## Update 89 (2026-08-10): example 16 was unrunnable, and two of its defects were silent

Running example 16 to check stability after the interface fix (findings 88) required repairing four
independent pre-existing faults. None were caused by the interface change; the example had simply been left
behind by the CHARMM-GUI prep redesign.

1. **Stale CLI in `run_sim_hybrid.sh`.** It passed `--bilayer-pdb`, `--protein-placement-mode`,
   `--protein-surface-gap`, `--protein-lipid-min-gap`, `--protein-lipid-cutoff-step` and
   `--protein-lipid-cutoff-max`, none of which exist any more. The membrane now comes from
   `--membrane-pdb` + `--lipid-name`; `parameters/dryMARTINI/DOPC.pdb` is itself a CHARMM-GUI product
   (CRYST1 50.09 x 50.09 x 85.0) and works directly, with its water and ions stripped on read.
2. **`str(None)` forwarded as an argument.** `prepare_workflow_hybrid_artifacts` forwarded
   `--membrane-top str(args.membrane_top)` and `--opm-reference str(args.opm_reference)` unconditionally, so
   an unset optional became the literal string `"None"` and the child opened a file named `None`.
   Optional arguments are now appended only when set.
3. **NumPy 2.** `ndarray.ptp()` was removed in NumPy 2.0 and the venv is on 2.4.6.
4. **The stage-7.0 solvation gate was unsatisfiable without `--opm-reference`.** `belt_half_thickness_a`
   fell back to `0.0`, so `|z - midplane| <= 0` selected nothing and the gate raised "found no backbone
   sites inside the belt" on a *correctly inserted* protein -- 26 of 31 BB sites were within 15 A of the
   midplane at the time it fired. The belt is the protein's published hydrophobic band, so it does not
   exist without an OPM reference: the fallback is now `None`, the gate reports that it was skipped and
   that solvation was NOT verified, and the genuine empty-belt error now names the half-thickness, the
   midplane, the closest BB offset and the BB span so it is diagnosable. Micelle mode still requires
   `--opm-reference` up front, as before.

Separately, **every VTF the workflow has ever written mislabelled the environment.** The metadata lookup was
hardcoded to `example/16.MARTINI/pdb/<id>.MARTINI.pdb`, but the workflow writes that file to
`<run_dir>/hybrid_prep/`, so it was never found and all 3612 lipid beads came out as `UNK` -- silently, with
positions intact, so the trajectory looked fine while being unselectable by lipid in a viewer.
`find_martini_metadata_pdb` now searches the run directories implied by the `.up` paths passed in, every
candidate derived from an explicit argument so no system can inherit another's labels.

Stability after the findings-88 interface fix, 20000 steps at dt = 0.009 and T = 0.8 T_up (280.5 K),
from the stage-7.0 coordinates of each build:

| | 1rkl | 1AFO |
|---|---|---|
| non-finite / abs(V) > 1e6 frames | 0 / 0 | 0 / 0 |
| potential | -18234 -> -19789 | -18208 -> -19412 |
| Rg(CA) | 12.76 -> 13.03 A | 15.86 -> 14.40 A |
| worst peptide C-N | 1.62 A, 0 above 2 A | 1.66 A, 0 above 2 A |
| KE/1.5kT mean | 1.021 | 1.032 |

For 1rkl the bilayer was checked directly: two leaflets with exactly 258/258 headgroups, peaks at +16.8 and
-15.4 A, thickness relaxing 34.6 -> 30.5 A, and zero lipids leaving. The large `|pos|max` in the
environment is unwrapped **ion** diffusion, not lipid loss -- worth knowing before reading it as evaporation.

The kinetic excess is still +2 to +3%, so the findings-88 prediction that it would collapse to 1.000 is
**refuted**: the one-sided O force was not its cause and the real one is still unidentified.

Lesson: a silent fallback is worse than a missing input. Both the `0.0` belt half-thickness and the
unfound metadata PDB produced confident, wrong behaviour -- one an error that accused the geometry, the
other a trajectory that looked complete. Neither said "the thing I needed was not there."

## Update 90 (2026-08-13): the NaN trigger is still unidentified; the LJ core table has a zero-force plateau

The four glpG-DDM micelle REMD variants all died at block 2--3 of 12. What the trajectories show, and
what each measurement rules in or out.

**The event.** Located the origin by scanning all 48 slots of `glpG-RKRK-79HIS_S115T` for the earliest
non-finite frame: slot 45 (T = 0.890978, near the top of the ladder), group `output`, frame 129. It then
appears one frame later in slot 38, then 32, 26, 20, 14, 8, 2, then back up through 3 and 9 -- the wreck
random-walking the ladder by exchange, destroying every slot it lands in. That is why 48/48 replicas end
up destroyed from a single event, and why a slot looks healthy again the frame after it went NaN.

It is a **single-frame catastrophe from a fully healthy state**. Slot 45 frame 128: potential -7070,
protein/lipid KE 1.51/1.37, max|pos| 489 A. Frame 129: positions non-finite, `lipid_kinetic` = 6.246e18.
Frames are ~60 steps apart, so the whole event fits inside 60 steps. `lipid_kinetic` blowing up while the
protein KE is merely NaN is the same signature as Update 88, but that cause is excluded below.

**Ruled out, each by measurement rather than argument.**

* *The Update-88 five-fold backbone over-count.* Both the local `fixcheck_1rkl` build and the actual
  cluster seeds carry **only BB** in `martini_potential.pairs` (31 and 210 sites, one per residue). The
  fix is in these runs. Corollary: the backbone `O` site being the closest protein atom to the
  environment (2.97--3.32 A, vs 4.0--4.4 A for BB, and all 8 closest protein-env contacts are O-to-tail)
  is *expected*, not a defect -- O carries no MARTINI interaction, so nothing repels it. An earlier draft
  of this entry read that proximity as the cause; that was wrong.
* *A stale pair list.* `cache_buffer` = 2.0 A and `pairlist_needs_rebuild` is called before every force
  evaluation, so the list is at most one step old and the skin covers 2 x 1.0 A of approach.
* *Minimum image.* `simulation_box::minimum_image` uses `roundf(dr/box)`, which is correct for
  arbitrarily large separations -- relevant because positions run to 489 A in a 137.4 A box (unwrapped
  ion diffusion, `--disable-recentering`), i.e. 3.5 box lengths, which a single-shift implementation
  would have got wrong.
* *The table build formula.* The tabulated grid reproduces its own analytic expression (potential-shifted
  LJ + reaction-field Coulomb, `roundf`-free) to a max relative deviation of **2.3e-12**.
* *Timestep instability in the sampled range.* At the tabulated forces, one-step displacement at
  dt = 0.009 for a 72 m_up lipid bead is **0.0058 A at r = 3 A** and **1.2e-4 A at r = 4 A**. Nothing
  marginal there.
* *Thermal access to the catastrophic region.* Reaching r = 3 A against 4*eps*[(sig/r)^12 - (sig/r)^6]
  costs ~450 E_up, about 500 kT at T = 0.89; r = 1 A costs ~6.5e8 E_up. Unreachable.
* *Reproduction from the last healthy state.* 4000 steps locally, single replica, from slot 45 frame 128
  at T = 0.890978 and dt = 0.009: stable, KE/1.5kT = 1.065. This test is **weak** and does not exonerate
  anything -- the trajectory stores no momenta (`run_remd.py` re-thermalises each chunk), so the exact
  microstate cannot be rebuilt, and the frame-129 configuration arrived from another slot anyway.

**What the table does contain.** `py/martini_build_tables.py:533` floors the radius,
`r = max(r, 0.1 * sig)`, over a grid that starts at `PARTICLES_R_MIN_A = 0.0`. Consequences, measured on
a real build (sigma up to 6.0 A, 1000 points, dr = 0.0120 A):

* 50 grid points are inside the floor, tabulating a **constant** 3.16e12 E_up -- so the force is
  **exactly zero for r < 0.1*sig** (0.47--0.60 A). A pair that gets inside feels no restoring force at
  all, which turns a recoverable close approach into free interpenetration.
* Grid maximum **6.18e12 E_up**, and the steepest tabulated force is **1.10e14 E_up/A at r = 0.48 A**,
  giving a one-step displacement of **1.2e8 A** for a lipid bead at dt = 0.009.

This violates two explicit rules: the spline table must equal the published functional form exactly, and
no arbitrary capping. It is a genuine defect and it is what converts a close approach into a 1e18
kinetic energy. **It is not, however, demonstrated to be the trigger** -- nothing measured explains how
a pair crosses from the thermally accessible region (>= 3 A, ~500 kT of margin) into it. Fixing the
floor is a change to the force-field table's domain and spacing, so it needs a decision rather than a
unilateral edit: removing the floor alone puts an infinity at r = 0, and simply raising `r_min` moves the
zero-force plateau to the new inner knot, because the clamped spline returns a constant with zero
derivative below its domain.

Open, and the honest state of it: the trigger is unknown. The next measurement that would settle it is
per-step instrumentation inside a running REMD ladder (max |force| and min interacting-pair distance per
step, per slot), because the 60-step frame interval and the missing momenta make the stored trajectories
unable to resolve the onset.

## Update 91 (2026-08-13): MBAR silently returns uniform weights for a hybrid coupled potential

`helpers/calc_hdx_ht.py` and `4.calc_D_uptake.py` built `beta[l] * cE0[k]` from the raw energies with no
reference subtraction. For the protein-only potentials they were written against, O(1e2--1e3), that is
fine. `00.AnalysisScripts/README.md` requires a hybrid trajectory's `Energy.npy` to be the full
coupled-system potential, which for a protein plus a DDM micelle is ~-7.6e3 E_up: `beta*U` reaches
-1.2e4, `exp(-u)` overflows, and the solver never leaves `f_k = 0`.

Measured on `glpG-RKRK-79HIS_S115T` (48 states, 423 frames): raw gave f_k spread **0.000**, neighbour
overlap **0.0000**, and **0/423** columns carrying weight at every target temperature except the bottom
rung; mean-subtracted gave f_k spread **71.07** rising monotonically, neighbour overlap 0.115--0.128, all
423 columns weighted, ESS 671--3378 of 20304, dG median 2.2--3.8 kcal/mol.

The failure mode is the dangerous part: `f_k = 0` makes every weight equal, so the estimator returns an
unweighted average over the whole pooled ladder while reporting it as a reweighted ensemble at one
temperature. It raises nothing. Tell-tales are a gradient norm of exactly `sqrt(n_rep-1) * n_frames`,
`max_delta` of exactly 0, and ESS of exactly `n_rep * n_frames`. Two of the four variants had already
produced plausible-looking dG plots this way.

Fixed in both files by referencing `cE0` to its pooled mean, which is exact (f_k shifts by -beta_k*C and
the exp(beta_target*C) on each weight cancels in the normalisation) and leaves the protein-only path
numerically identical. `03.TrajectoryAnalysis/2.mbar_meltingCurve_freeEnergy.py` and
`04.HDX/4.calc_HDX.py` have the same construction but are byte-identical to master and are only ever fed
protein-only energies, so they were left alone under master parity.

Also note for anyone reading those solver messages: passing the pymbar-3 style 3D `u_kln` under the
installed **pymbar 4.0.3** is *not* a bug -- 3D and 2D give identical `f_k` (max abs difference 0).

Lesson: when a solver reports a gradient that is an exact function of the array shape
(`sqrt(K-1) * N` here) rather than of the data, it has not solved anything. Check that before reading any
number downstream of it.

## Update 92 (2026-08-13): the NaN trigger — the LJ core table was force-free, and particles reached it

Findings 90 left the trigger unidentified and concluded the catastrophic core region was ~500 kT out of
reach. **That conclusion was wrong.** Per-step instrumentation (`UPSIDE_MARTINI_PAIR_DIAG`, jobs
53324867/53324868: 48 replicas, exchange disabled so a blow-up stays in the slot that made it, 340 000
steps, 9.2 h) measured the opposite.

| | 53324867 | 53324868 |
|---|---|---|
| reported approaches < 1 A | 4742 | 5032 |
| approaches < 0.6 A (inside the floored plateau) | 1446 | 1432 |
| closest approach | **0.0355 A** | **0.0466 A** |
| largest force delivered anywhere | 3.4e11 E_up/A | 1.0e11 E_up/A |

The floor was the defect. `martini_build_tables.py` evaluated the grid at `r = max(r, 0.1*sig)` on a
domain starting at r = 0, so below 0.1*sig (0.47-0.60 A) the tabulated potential was a **constant** and
therefore exerted **no force at all**. At 0.0355 A the true dry-MARTINI LJ force is **5.8e29 E_up/A**;
the table's largest delivered force anywhere in the run was 3.4e11 -- about **18 orders of magnitude too
weak**. The clearest single event: a pair at **0.0804 A while the whole box's maximum force was
7.28 E_up/A**, i.e. that pair felt nothing. 93% of sub-0.6 A events coincided with a >1e6 E_up/A force
elsewhere (a neighbour being launched); 3% had no force above 1e3 anywhere.

Offending pairs are **environment-environment** (LIPID-LIPID and LIPID-ION; PROTEIN 0-1049,
LIPID 1050-2723, ION 2724-3155), which is why `lipid_kinetic` is what explodes while the protein KE is
merely NaN -- the Update-88 signature, reached by a different route.

Why findings 90's barrier argument failed: it computed one-step displacement for an **inertial** particle
of mass 72 (0.0058 A at r = 3 A). But ION and LIPID are integrated as **overdamped Brownian**; only the
420 protein atoms are velocity-Verlet. An overdamped step is proportional to the force, so a large force
gives a large displacement that can overshoot *through* a partner, and once inside the force-free plateau
there is nothing to eject it. **Entry by overshoot is inferred, not measured**; the floor removing the
exit is measured. Lesson: check which integrator governs the particles before computing a stability
margin for them.

**Fix.** Two coupled changes, because the domain was declared in one place and assumed in another:
* `py/martini_build_tables.py`: floor removed, `PARTICLES_R_MIN_A = 0.0 -> 0.3`, grid built vectorised
  over the true potential, and a new assertion that every grid point equals the analytic
  potential-shifted-LJ + reaction-field form to 1e-12 relative (CLAUDE.md requires asserting this rather
  than assuming it). 0.3 A is far inside anything reachable: the core there is ~5e17 E_up/A.
* `src/martini_potential.cpp`: `r_min`/`r_max` were **hardcoded to [0, 12]** and ignored the
  `r_min_ang`/`r_max_ang` attributes the builder already wrote. Changing the builder's domain alone would
  have silently mapped every distance onto the wrong knot. It now reads the domain from the table and
  validates it.

Verified locally against a real glpG-DDM system with the corrected table patched in: initial potential
**-7811.76**, min pair distance 4.0349 A and max force 33.3441 E_up/A -- **identical to the old table**,
so the change is a no-op everywhere the system actually samples, while the core now delivers
8.5e15 E_up/A at 0.35 A and 7.3e13 at 0.5 A where the old table delivered ~0. Old .up files still read
correctly (their own attrs say r_min = 0), so the C++ change is backward compatible.

Residual risk, stated plainly: the clamped spline still returns a constant with zero derivative below its
domain, so a single overdamped step that lands under 0.3 A would still find no restoring force. The dead
zone is halved (0.6 -> 0.3 A) and now sits behind a 5e17 E_up/A wall, but it is not eliminated. Whether
the overdamped step size for ION/LIPID can produce such a jump is the open question.

## Update 93 (2026-08-14): three defects the first corrected-table run exposed

**1. The MD loop destroyed its own error messages.** `src/main.cpp:1252` integrates systems under
`#pragma omp parallel for` with no exception trapping. An exception cannot leave an OpenMP structured
block, so any `throw string(...)` from the engine mid-run called `std::terminate`: the local POPE/POPG
run died at step 155 640 reporting nothing but `libc++abi: terminating due to uncaught exception` and
`Abort trap: 6`. The setup loop at line 929 already documents this hazard and traps for exactly this
reason; the integration loop did not. Now traps per system, keeps the first message, and rethrows once
serial so the existing top-level handler prints it with the round number.

**2. `martini_hdx_membrane_accessibility.py` defaulted to one system's bead names.**
`--tail-bead-names` defaulted to `C1,C2,C3` -- the DDM tails. On POPE/POPG (`C1A..C4A`, `C1B..C5B`,
`D3B`) it found zero tail beads and raised. It raised loudly here, but that is luck: a lipid that happened
to contain a bead called C1 would have been scored on the wrong subset silently. The default is now empty
and the tails are detected from the trajectory by the MARTINI apolar naming pattern `^[CD]\d[A-Z]?$`,
which covers a single-tailed detergent and a two-tailed lipid alike; an explicit list still overrides.
Verified: auto-detection on POPE/POPG selects C1A,C1B,C2A,C2B,C3A,C4A,C4B,C5B,D3B (2511 beads) and gives
an accessibility array bit-identical to passing those names by hand.

**3. Removing the table floor did not stop core penetration -- it changed the consequence.** On the
corrected table the same system still reached **0.2078 A**, with 32 approaches under 1.0 A and 6 under
0.3 A, and forces up to **1.27e15 E_up/A**. Under the old floored table those pairs coasted through
force-free; now they receive an enormous impulse and the run dies. So findings 92's fix removed a real
defect (the table is now the published potential) but the **entry** mechanism is untouched, and the
overdamped Brownian overshoot named there as "inferred, not measured" is now the outstanding cause.
Note the residual dead zone below the new r_min = 0.3 A was entered, which is exactly the risk that entry
recorded.

Also fixed while chasing these: `inject_particles_table` restated the radial domain instead of copying it
from the source table (see findings 92), and the C++ hardcoded the 1000-point grid size in four places;
the grid geometry is now declared once by `py/martini_build_tables.py` and carried through.

Lesson: a diagnostic that cannot report is worse than no diagnostic. Two of these three were only findable
because the pair instrumentation was still on -- the crash itself said nothing at all.

## Update 94 (2026-08-15): the BB proxy is now a derived site, and the old placement is gone

Update 88 fixed the hybrid interface by propagating the O/BB sensitivities through a hand-written
placement Jacobian. That worked but was fragile, and the user's preference was the other route: build the
BB proxy out of nodes Upside already differentiates. That is now what happens, and the hand-rolled path
has been deleted rather than left as a fallback.

**What the proxy is now.** `HybridPositionNode` takes two arguments, `pos` and `infer_H_O`, and BB is the
mass-weighted centre of N/CA/C/O with weights [14,12,12,16]/54, where the O is *Upside's own derived
carbonyl O* rather than one rebuilt from a stored local frame. Because every term is then a linear
combination of node outputs, `propagate_deriv` is a constant-weight split: the N/CA/C shares go to
`pos.sens`, the O share goes into `infer_H_O`'s sensitivity, and infer_H_O's own chain rule carries it
back to CA, C and the next residue's N. That is where the frame's rotational dependence lives, so no
placement Jacobian is needed here at all -- 184 lines of frame/Jacobian helpers were removed.

**The counter-force is not special-cased.** `martini_potential` runs on the hybrid node's output, which is
a pass-through copy of `pos` with only the BB and O slots overwritten. So an environment particle's
sensitivity is returned to `pos.sens` unchanged (the `derived[]` test in `propagate_deriv`), while the
protein side is redistributed by weights that sum to 1. Newton's third law therefore holds by
construction. Measured on 1rkl: force remaining on the BB slots and the O slots is exactly **0**, and the
total force sums to **1.69e-09** of the largest single force.

**The C-terminal residue has no acceptor.** `infer_H_O` builds a carbonyl O from the *next* residue's N,
so it emits `n_res - 1` acceptors. That one BB is the N/CA/C mass centre with the weights renormalised,
and its O slot is left as it came in. This is not a guard -- there is genuinely no such site -- and no
MARTINI term reads that slot.

**Also fixed: the Python engine computed a different model than the binary.** `engine_c_library.cpp` never
called `load_masses_for_engine` / `register_fix_rigid_for_engine` / `register_stage_params_for_engine` /
`register_hybrid_for_engine`, which `main.cpp` does. Any analysis through `upside_engine` was therefore
evaluating a system with the hybrid interface inactive: **-12827 vs -18172 E_up** on the same
coordinates. Every `get_output`/`energy` result taken through the Python engine before this is suspect.

**Deployment has a compatibility cost, and it is not optional.** A two-argument node cannot load a config
that declares one, so every existing `.up` on the cluster became unloadable the moment the binary was
replaced. `py/martini_upgrade_hybrid_args.py` rewrites that attribute in place; it refuses a file with no
`infer_H_O` rather than producing a config that cannot load. Lesson: when a node's arity changes, the
migration is part of the change, not a follow-up -- and a config held open by a running job can only be
migrated in the gap before the next block, so the migration has to live in that job's submit script.

## Update 95 (2026-08-15): the NP footprint answers the albumin lysine claim, but only for one orientation

Testing Carlson et al. 2025 against the 1AO6 / 5 nm MPA-AuNP runs (12 262 frames sampled from ~1.33 M
across six orientations, CB-anchored contact at 8 A, burial as CB coordination within 10 A on the native
first frame). The paper's albumin claims are decreased labelling at **K12, K73, K190, K525** in water and
**K190, K525, K541** in TES, a common site "around K190 (Subunit II)", and a preference to "open and expose
its center". Construct index = HSA number - 5; all five map to lysines.

**The window decides the answer, so it has to be chosen before looking.** Albumin does not merely deform
here, it unravels: Rg goes 26.3 A -> 39 / 50 / 86 / 103 / 152 / 172 A across the six orientations. Once
spread, contact stops discriminating anything -- core and surface mean contact are 0.272 vs 0.261, i.e.
equal -- so a whole-trajectory average would "support" almost any site claim. Restricted to adsorbed and
still-compact frames (Rg < 1.25x native) the same measure separates properly: core 0.128 vs surface 0.245.
Only **393 of 12 262 frames (3.2%)** are in that window.

* **K190 is not supported.** In the compact window its contact is **0.000**, in every orientation
  separately. In the spread state it is 19th percentile of 58 lysines and *further* from the particle than
  most (72nd percentile by mean minimum distance). Its measured burial is 18.3, the **50th percentile** --
  so the "in the core" premise is itself weak for 1AO6 by CB coordination.
* **K525 and K541 are supported**: contact 0.545 (81st percentile) and 0.690 (91st) in the compact window.
  The paper already says of K525 that it "better agrees with the computational projects", and this agrees.
* **"Open and expose its center" is supported**, and not merely as a by-product of unravelling: even
  within the compact window core residues lose **9.19** CB neighbours against **0.37** at the surface.
* **The site evidence is one orientation.** The compact window is 71% orientation 2, and the contacting
  patch -- a contiguous C-terminal run K536/K538/K541/K545/K557/K560/K573/K574 plus THR566, CYS567,
  ALA569, GLU570/571, ASP563 -- scores 0.96-1.00 in orientation 2 and **0.00 in all five others**. This is
  one binding pose, not a demonstrated preference.

**Independent cross-check.** Orientation 2 is 180-0-0, and the group's earlier BSA(3V03) MARTINI/DPD run
recorded for that same orientation "LYS535 and LYS537 facing NP". A different protein and a different
model land on the same face, which is real corroboration of the pose rather than of the K190 claim.

**Caveat that limits all of it:** the experiment labels albumin whose CD still shows its secondary
structure, and these runs go far past that. Orientations 1 and 3 reach Rg 152 and 172 A in a **200 A box**,
so they self-interact through the periodic image and are unusable for any structural conclusion. Testing
site preference properly needs sampling of the adsorbed-but-folded state, which these runs barely visit.

Method note: contact is scored at CB because that is where `martini_sc_table_1body` anchors the
sidechain-environment term, and the reconstruction (Kabsch fit of the stored `affine_alignment` reference
onto N/CA/C, then the fixed CB offset) matches the engine to **1.5e-4 A**. A CB cutoff under-counts lysine
relative to short residues by ~6 A of sidechain, so rankings are taken lysine-against-lysine.

## Update 96 (2026-08-15): the flat HDX profile was an analysis ceiling, not undersampling

Retracts the diagnosis in the same session's earlier note that the glpG dG profile was limited by sampled
amide opening events. It was limited by a threshold in `calc_hdx_ht.py`.

`residue_dg_profile_with_error` has two conventions. `residue_dg_from_pf_jscripts` sends only
`mean_pf >= 1.0` to the +1000 sentinel; `residue_dg_from_pf_step6_plot` sends everything above
**0.99999**. Since dG = 0.001987 * temp_scale * ln(pf/(1-pf)), that second clip is a hard ceiling of
0.001987 * 298 * ln(99999) = **6.82 kcal/mol** at T = 0.85. The T-slice was calling it -- even though the
jscripts helper's own comment says it is "the legacy jscripts residue-dG mean calculation used for
DG_res_T_slice" -- so the profile could not exceed 6.8 no matter what the data said.

Measured on the POPE/POPG wildtype ladder at 10% of the run, same frames both ways:

| | T = 0.75 | T = 0.85 |
|---|---|---|
| clipped path, dG max | 5.94 | 6.31 |
| jscripts path, dG max | **17.35** | **19.58** |
| residues with `0.99999 <= pf < 1` (discarded by the clip) | 77 | 44 |
| residues with `pf == 1.0` exactly | **0** | **0** |

So every "never-exchanging" residue was manufactured by the clip, and the profile's compression to -5..+7
against a reference peaking at 10-20 was entirely this. Two things made the wrong diagnosis look right:
the sentinel count did fall with more sampling (65 -> 44 at T = 0.85, because more frames push more
residues below 0.99999), and the ceiling is close enough to a plausible statistical limit that
`0.5921 * ln(ESS)` lands in the same range. Neither is evidence for undersampling.

Why pf may legitimately sit far closer to 1 than 1e-5: reweighted MBAR weights span many orders of
magnitude, so a residue that is unprotected only in high-temperature frames has those frames
exponentially down-weighted at a low target temperature. The jackknife error bar is what reports how well
resolved that is; a clip just throws the residue away. Fixed by using the jscripts convention for the
T-slice, which also makes it consistent with the full-temperature profile that already used it.

**Residual, and genuinely sampling-limited:** 32 residues at T = 0.85 and 75 at T = 0.70 still reach
`pf == 1.0` in float64 because the unprotected frames underflow to zero weight at the far target
temperatures. Those are real resolution limits and do improve with more frames -- unlike the clip.

Lesson: when a computed quantity saturates at a suspiciously round value, check the code path for a
threshold before concluding anything about the physics. Recomputing the same frames both ways cost one
minute and overturned a diagnosis that had already redirected a weekend of compute.

## Update 97 (2026-08-15): the NaN cascade, caught with diagnostics on -- the sub-0.3 A plateau is the sink, not the source

The local POPE/POPG wildtype ladder lost replica 15 (T = 0.90, the top rung) at step ~145 700 of 1.3 M,
with the pair diagnostics running. This is the first capture of the whole sequence, and it **revises
findings 92/93**, which named the table's zero-force core as the trigger.

Causal order, from the diagnostic's own step counter:

| diag_step | min_dist | pair | max_force |
|---|---|---|---|
| 5205397 | **1.83 A** | 1191e-1438e (lipid-lipid) | 3.4e6 |
| 5205403 | 1.59 A | 1022p-3920e (BB proxy-lipid) | 1.2e7, pot 1.86e6 |
| 5206114 | 1.05 A | 1146e-2701e | 5.5e10 |
| 5206469 | **0.245 A** | 1155e-3033e | 1.3e12 |

The first sub-0.3 A approach arrives **1072 diag-steps after** the first enormous force. So the flat region
below the table's inner knot is where the cascade **ends up**, not what starts it. The trigger is two
**environment** beads reaching 1.83 A -- comfortably inside the valid table domain, where the tabulated
force is correct and simply huge. Removing the r = 0.1*sig floor (findings 92) was still right, but it was
never going to prevent this.

The escalation is visible in the logged scalars too: total potential -2.17e4 (frame 2332) -> **+5.4e3**
(2333) -> recovers -> **+1.98e6** (2428) -> NaN (2429). It oscillates in and out for ~96 frames before
diverging, which is why a single spot check can miss it.

**Two health checks that both fail on this, and the one that works.** `isfinite` passes for 96 frames after
the system is destroyed. The documented peptide C-N bond count also passes -- frame 2428 has 2 broken bonds,
because the protein is intact and the damage is in the lipids. What separates them is the **sign of the
total potential**: a condensed bilayer plus protein sits near -2.2e4, and a positive total means core
overlap. That is physics, not a tuned cut. `scratchpad/local_popg_79HIS/reseed.py` applies both, since they
catch different failures, and picked frame 2425 (-2.18e4) where finiteness alone would have picked the
destroyed frame 2428 and bond-count alone would have agreed with it.

**Diagnostic reading note.** After the blow-up the diagnostic prints `max_force 0` with pair indices `-1`,
because every force is NaN and `NaN > diag_max_force` is false, so nothing updates the accumulator. That is
a post-mortem signature, not a measurement of zero force. Likewise a repeated `min_dist 0.0010 (840p-1050e)`
is not a physical 0.001 A contact -- the same atoms are 83-96 A apart in the stored coordinates. Do not
read the diagnostic's post-NaN lines as evidence about the mechanism; only the pre-NaN escalation is.

Restart cost nothing: `reseed.py` rotates `/output` to `output_previous_<n>` and reseeds `/input/pos`, and
`py/martini_remd_concat.py` joins the segments at analysis time.
