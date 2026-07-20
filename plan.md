# Project Goal

dryMARTINI hybrid membrane in Upside. Full-resolution 14-bead DOPC + protein with a timestep
validated for the coupled hard-core interface. Single-particle CGL subsystem ABANDONED (removed).

CURRENT STABILITY BASELINE: every physical degree of freedom uses one g-JF force update at the standard
Upside numerical timestep `dt=.009`; there is no lipid substep or displacement cap. Preparation stages use
ordinary dry-MARTINI damping. Production maps the native 4 ps dry-MARTINI relaxation to a 16 ps physical
relaxation after the factor-four correction and therefore to `.4` numerical steps (`tau_up=.0036`). This
FDT-consistent friction is applied to all environment beads and to real protein N/CA/C carriers initially
inside the existing 12 A protein--environment spline range. Noncontact protein carriers retain the standard
Upside OU bath. Regenerated O and BB sites remain virtual and BB-env forces are mapped through the complete
coordinate Jacobian. BB-env and SC-env remain active and spline-table-only. This is a contact-local friction
clock; it is not a claim that molecular DOPC diffusion reaches the 40 ps/step target.

# Architecture & Key Decisions

- Master reference for diffs: /Users/yinhan/Documents/upside2-md-master (NO martini files ->
  entire martini subsystem is Clean-Slate scope; keep it impeccably clean, remove dead code).
- C++ (kept): martini_potential, dist_spring, angle_spring, martini_masses, martini_brownian,
  martini_fix_rigid, martini_hybrid, martini_stage_params, sc_table/particle tables. CGL removed.
- Integrator: one all-particle g-JF update per `.009` numerical timestep. Production uses the calibrated
  dry-MARTINI friction clock for environment beads and initially lipid-contacting protein carriers; carriers
  with zero g-JF friction use the ordinary Upside OU bath. Positive-friction particles are excluded from OU,
  so no degree of freedom is double-thermostatted.
- Unit contract: native dry-MARTINI -> Upside conversion happens ONCE at Python h5-build; runtime
  h5/config store Upside-unit values; C++ engine does ZERO unit conversion.
- H5 FF: parameters/ff_2.1/martini.h5 (/particles + /sc_table). No version numbers; back up + overwrite.
- Cell-list neighbor search (O(N)); bilayer baseline potential -11039.306641 (bit-identical gate).

# COMPLETED (P1-P11, condensed — see findings.md for detail)

- [x] C1-C6 Clean-slate excision of CGL (C++, Python, H5, ITP, MD scripts, docs); build passes.
- [x] BB-env reverted to raw LJ grid (PMF drove protein out of bilayer); standard timestep 0.009.
- [x] Cell-list pairlist (O(N), bit-identical); min-M study (hybrid needs M~90 overdamped).
- [x] g-JF single-step lipids implemented (findings 39-45); the later transport study rejected a
      universal effective-time factor because diffusion and rotation require different mappings.
- [x] P10 finalization: tex rewritten for the final method; master-parity cleanup (MASTER files
      minimal/gated; NEW martini files clean); convert_stage split; dead C++/Python removed.
- [x] P11 unit-baking: conversions moved to Python h5-build; C++ unit-agnostic; bit-identical gate.
- [x] Historical Update 49 stabilized the mixed integrator by putting all real particles on the g-JF
      schedule. Phase G supersedes its root-cause interpretation with the virtual-coordinate and
      timestep evidence below.
- [x] Update 50: C++ analytic LJ/Coulomb removed -> spline-only engine (bit-identical).

# COMPLETED PHASE (2026-07-20)

## H. Match DOPC transport at the standard Upside timestep (2026-07-20)

- [x] Verify the standard used by `example/02.ReplicaExchangeSimulation`: default single-force timestep
      `dt=.009`; its ordinary Verlet integration groups three such timesteps between outer bookkeeping calls.
      The hybrid must use the same `.009` numerical timestep for both protein and bilayer.
- [x] Define the physical clock per numerical timestep: 40 ps protein time and 10 ps raw MARTINI time after
      the factor four. For 11.5 um^2/s, the molecular target is `MSD_xy = 0.184 A^2` per `.009` step.
- [x] Define the fallback observable correctly: every bead in a bonded DOPC has the same asymptotic diffusion
      as its molecular COM. A distinct sub-molecular clock therefore means calibrated bare/free bead mobility,
      not a different long-time bead diffusion coefficient.
- [x] Test the existing single-step all-particle g-JF path at `.009`. Reducing friction through the molecular
      free-draining estimate and nearly to zero saturates DOPC COM diffusion near 0.175 um^2/s (66x below
      target) while heating the interface, so molecular diffusion cannot define 40 ps/step on this path.
- [x] Select the accepted fallback wording and observable: the dry-MARTINI friction clock felt by the protein
      in the bilayer. Map the 4 ps raw relaxation through the factor four (16 ps physical), and apply its exact
      mass-scaled g-JF friction to environment beads and protein carriers within the 12 A spline cutoff.
- [x] Complete the selected single-step scheme and explicit H5 clock/transport metadata, plus a reproducible
      molecular-COM and bare-bead calibration workflow rather than a hard-coded multiplier.
- [x] Verify integration wiring, kinetic temperature, molecular DOPC diffusion, and hybrid secondary structure
      over 50,000-step continuations of freshly prepared 1RKL and 1AFO systems.

### H Decisions

- The final observable is physical DOPC lateral diffusion after the MARTINI time correction, not a
  nominal metadata rescaling. The factor four must enter the clock/calibration exactly once.
- Equilibrium potentials, hard cores, SC-env, and BB-env are immutable. Only masses, dissipative/noise
  dynamics, and multirate integration may be calibrated.
- Use the same `.009` numerical timestep as the ordinary Upside examples and prefer the existing one-step
  all-particle g-JF architecture, which applies full SC-env/BB-env reactions through the virtual-site Jacobian.
- Molecular COM diffusion is primary and must always be reported. The accepted fallback is a friction clock,
  not a diffusion claim: raw dry-MARTINI relaxation 4 ps becomes 16 ps physical after the factor four, hence
  `tau_up = dt * 4/10 = .0036` at `.009` and `alpha_i = m_i/tau_up`. Protein carriers inside the same 12 A
  spline interaction range receive that FDT-consistent bath and are excluded from the ordinary OU bath.
- Multistep Brownian dynamics remains a documented last resort only. The existing overdamped fallback is not
  acceptable for protein friction because it discards intermediate protein--lipid reactions, and its
  displacement cap is prohibited.
- The protein contact set is classified once from the production starting structure. Making friction depend
  continuously on position would require a thermodynamically consistent multiplicative-noise formulation;
  a silent distance-switched Langevin coefficient is not used.

## G. Re-open secondary-structure and coupled-timescale fix (current user request)
- [x] Preserve all SC-env and BB-env interactions and the spline-only force-field path.
- [x] Reproduce and quantify the 1AFO/1RKL failure from the saved production outputs.
- [x] Audit the integrator/thermostat equations and configuration wiring; distinguish a numerical
      heat leak from a physical protein-bilayer friction model.
- [x] Replace the empirical `friction_scale=100` default because it is not derived or transferable;
      use a fluctuation-dissipation-consistent coupling that leaves the equilibrium potential intact.
- [x] Serialize regenerated O/BB coordinates from the hybrid coordinate node so trajectory analysis and
      stage restarts see the same geometry used by the force graph.
- [x] Select the largest converged outer timestep from matched-time 1RKL regressions; do not mask a
      discretization heat leak with extra protein friction.
- [x] Re-prepare both examples, verify fresh workflow wiring, and run current-code matched-time
      structure/temperature regressions from equilibrated production inputs.
- [x] Document the physically defensible time interpretation, its calibration procedure, and limits.

### Revised Decisions
- The empirical `friction_scale=100` protein thermostat is removed: its target-temperature crossover
  is system-dependent, so it cannot establish that the protein feels the membrane's physical drag.
- A single scalar "ps per step" must not be claimed for both components unless one observable-independent
  mapping is demonstrated. The final design must separate equilibrium correctness from kinetic calibration.
- Root cause found in the virtual backbone path: regenerated O and BB sites are not physical degrees of
  freedom. They must not be integrated or thermostatted, and the BB-env gradient must include the full
  derivative of the regenerated-O/BB position with respect to the real N/CA/C carriers. The initial
  finite-difference Jacobian still produced timestep-independent NVE drift and interface heating, so it
  is replaced by a forward-mode automatic derivative of the coordinate construction.
- The raw all-pairs filter is now explicit and stage-independent: production admits only the derived BB
  proxy as a protein--environment raw pair. SC-env remains fully active exactly once through its dedicated
  spline PMF; N/CA/C/O and fixed SC bookkeeping sites do not acquire discarded raw gradients.
- Saved `/output/pos` previously copied the raw integration array, whose virtual O/BB entries were stale.
  Logging now copies the already-evaluated `martini_hybrid_position` node instead; mutating the integration
  coordinates inside force evaluation is prohibited because it breaks the differentiable graph.
- The final fallback supersedes the earlier `.00225` mechanical-clock baseline: protein and bilayer now use
  the same standard `.009` numerical step, with contact-local production friction keeping both kinetic
  subsets at the target temperature without changing conservative forces.
- Particles with positive stochastic g-JF friction must be excluded from the global OU thermostat;
  zero-friction protein carriers remain in the standard Upside OU bath.
- The implemented fallback is a static contact classification with an FDT-consistent one-body Langevin bath,
  calibrated to the corrected dry-MARTINI relaxation clock. It is not a momentum-conserving pair thermostat
  and does not establish a full protein resistance tensor.

## A. Crash fix + clamp bug  -- DONE + verified
- [x] Clamp bug fixed (findings 51): friction_scale load `(v>0)?v:1` -> `>= 0` (0 = velocity-Verlet).
- [x] Rebuilt obj/; bilayer bit-identical (-11039.306641); hybrid stable (kT 1.118, Rg 15.65->15.64,
      protein folded, no NaN over 5000 steps with the protein now truly velocity-Verlet).

## B. Lipid diffusion check across T_up 0.5-1.5 (task 1)  [user active instruction]
- [x] Resolved the "too fast/too slow" unit confusion (findings 52): raw g-JF displacement is fixed.
- [!] LONG-RUN RESULT (findings 53): g-JF bilayer has NO diffusive plateau -- alpha falls 0.51->0.27
      (toward caging). Overdamped control (same box) is ALSO sub-diffusive (alpha 0.35) -> the caging is a
      SYSTEM/regime property, NOT the g-JF integrator. Neither integrator gives a clean Fickian D here; a
      diffusion timescale needs a much bigger patch + longer runs. No single-scalar transport timescale.
- [x] ALTERNATIVE anchors (findings 56-57): STATIC MECHANICS is the claim. Measured (g-JF DOPC, ~303K):
      APL 65.5 A^2 (exp 67.4), hydrophobic d_c 2.78 nm (exp 2.7-2.9, ON TARGET), chain order P2 0.30,
      tilt 22 deg / S 0.74 (fluid). Rotational ACF also cage-limited (same wall) but curve-matchable.
      K_A/kappa BLOCKED: barostat is a damped/clamped relaxation scheme, not fluctuation-correct (findings 57).
- [ ] OPTIONAL (user's call): further MC-barostat sampling for K_A and larger/longer transport runs.
      These are not blockers for REMD/metad equilibrium sampling.

## E. DO ALL 4 + analysis scripts (user directive 2026-07-19) -- DONE (findings 58)
- [x] E1 MC barostat (box.cpp type=2, COM-scaling, exact NPT Metropolis): WORKS -- area fluctuations
      restored (sd 0.95->29). K_A undersampled (block-Var grows, still falling from 1867; slow area relax ~
      sub-diffusive lipids) -> 500k-step run cooking for a converged number. APL relaxes 65.6->60 (FF condensed).
- [x] E2 temp sweep: STABLE every rung T_up 0.5-1.5; D monotonic in T (Arrhenius); ~303K rung D@270~13~DOPC.
- [x] E3 metadynamics node (src/metadynamics.{h,cpp} + main.cpp hook + CMakeLists + advanced_config): SMOKE
      PASS (drives CV, well-tempered deposition, FES logged). Quantitative FES-recovery needs a longer run.
- [x] E4 curve-match to all-atom CHARMM36 DOPC: POOR overlay -> confirms NO single effective-time factor
      (mode-dependent CG clock; rotation != diffusion).
- [x] E5 (user add-on) 00.AnalysisScripts: analyze_metadynamics.py + analyze_remd.py (both TESTED) + README.
      REMD end-to-end run PASSED (4 slots, g-JF lipid) -> task 2 confirmed at run level.
- ALL C++ (clamp fix, main.cpp REMD temp, MC barostat, metad) -> bilayer bit-identical, zero master regression.
- OPEN (need long sampling, not blockers): converged K_A number (running); converged metad FES-recovery.

## C. REMD for lipid (task 2)  -- DONE: works correctly (findings 54)
- [x] Verified: coord-swap REMD; per-slot T reaches lipid integrator; lipid energy in the criterion;
      no momentum rescaling needed. FIXED main.cpp uninitialized-temperature
      call (moved after set_temperature); rebuilt, bilayer bit-identical, hybrid stable. Latent mv-path
      gap flagged (not the active path).

## D. Metadynamics feasibility (task 3)  -- DONE: evaluated (findings 55)
- [x] Core: very feasible (~1-2d, new PotentialNode bias reading a CV node; most CVs exist). Hybrid:
      feasible but membrane geometry CVs must be built + collective-lipid CVs fragile (sub-diffusive).
      PLUMED stub is dead; native well-tempered (spline-grid bias) recommended. Design sketch in findings.

## F. Superseded diagnosis (replaced in findings 60)
- The earlier `friction_scale=100` workaround reduced apparent overheating but hid the actual problems:
  O/BB virtual sites were integrated and thermostatted, BB-env forces omitted the virtual-site chain
  rule, and lipids received both g-JF and OU noise. Phase G replaces that workaround at the source.
- ASIDE: long K_A (500k steps, MC barostat) still ~1886 mN/m (undersampled; slow area relax ~ sub-diffusive
  lipids) -> quantitative K_A needs far longer sampling / bigger patch (report as method+bound, not a number).

# Review

- Root causes resolved: virtual-coordinate graph/chain rule, stale virtual-site trajectories, double-bath
  bookkeeping, empirical protein friction, and unequal component clocks.
- Final fresh 5000-step `.009` workflows: mapped DSSP helix count remained 7/7 for 1RKL and 27/27 for 1AFO;
  mean protein/lipid kinetic energies were 1.286/1.297 and 1.300/1.298 against the 1.297 target.
- 50,000-step production continuations kept the same helix counts; mean kinetic energies were 99.9% and
  100.2% of target. O/BB momenta were exactly zero.
- Whole-DOPC COM diffusion is measured honestly: roughly `6e-5`--`7e-5 um^2/s`, not the 11.5 target. The
  final claim is the contact-local corrected dry-MARTINI friction clock.
- Build, Python compilation, shell syntax, diff whitespace, H5 metadata, constraint smokes, analyzer synthetic
  regression, fresh workflows, and manuscript compilation all pass.

# Known Errors / Blockers
- Molecular DOPC diffusion is not matched at 40 ps/step; the single-step friction scan saturates about 66x
  below target before the final high-friction clock is applied.
- The contact-local relaxation clock is not a quantitative hydrodynamic resistance tensor. Such a tensor
  would require independent parallel, normal, and rotational reference data and a momentum-conserving model.
- DO NOT break the full-res path or master behavior. Back up any h5 before overwrite. Git READ-ONLY.
