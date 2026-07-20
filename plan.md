# Project Goal

dryMARTINI hybrid membrane in Upside. Full-resolution 14-bead DOPC + protein at the standard
Upside step. Single-particle CGL subsystem ABANDONED (removed).

METHOD (final): protein = standard Upside velocity-Verlet (dt=0.009); lipids = g-JF inertial-
Langevin SINGLE step (langevin=1) with an effective-time factor mapping sim-time -> physical time.
BB-env pairs use the RAW dry-MARTINI LJ grid; SC-env uses the PMF. All pair interactions are
spline-table-only. Overdamped-Brownian RESPA (langevin=0, M=90) retained as the exact-D calibration
reference and default-off fallback (master parity).

# Architecture & Key Decisions

- Master reference for diffs: /Users/yinhan/Documents/upside2-md-master (NO martini files ->
  entire martini subsystem is Clean-Slate scope; keep it impeccably clean, remove dead code).
- C++ (kept): martini_potential, dist_spring, angle_spring, martini_masses, martini_brownian,
  martini_fix_rigid, martini_hybrid, martini_stage_params, sc_table/particle tables. CGL removed.
- Integrator: g-JF at friction_scale=0 reduces EXACTLY to velocity-Verlet -> protein rides the
  same g-JF path (friction_scale=0, thermostatted by OU); lipids friction_scale=1. This unifies the
  integrator across the hard protein-lipid interface (the fix for the mixed-integrator crash).
- Unit contract: native dry-MARTINI -> Upside conversion happens ONCE at Python h5-build; runtime
  h5/config store Upside-unit values; C++ engine does ZERO unit conversion.
- H5 FF: parameters/ff_2.1/martini.h5 (/particles + /sc_table). No version numbers; back up + overwrite.
- Cell-list neighbor search (O(N)); bilayer baseline potential -11039.306641 (bit-identical gate).

# COMPLETED (P1-P11, condensed — see findings.md for detail)

- [x] C1-C6 Clean-slate excision of CGL (C++, Python, H5, ITP, MD scripts, docs); build passes.
- [x] BB-env reverted to raw LJ grid (PMF drove protein out of bilayer); standard timestep 0.009.
- [x] Cell-list pairlist (O(N), bit-identical); min-M study (hybrid needs M~90 overdamped).
- [x] g-JF single-step lipids implemented + Route C chosen (findings 39-45): exact D impossible
      single-step on hard cores; approximate D via effective-time factor is the accepted path.
- [x] P10 finalization: tex rewritten for the final method; master-parity cleanup (MASTER files
      minimal/gated; NEW martini files clean); convert_stage split; dead C++/Python removed.
- [x] P11 unit-baking: conversions moved to Python h5-build; C++ unit-agnostic; bit-identical gate.
- [x] Update 49 crash fix: hybrid production NaN root-caused = mixed g-JF(lipid)/Verlet(protein)
      coupling injecting energy at the interface. FIX (prep, gated langevin=1): all atoms on g-JF,
      protein friction_scale=0 (= velocity-Verlet + OU), lipids =1. Verified stable 5k-20k steps.
- [x] Update 50: C++ analytic LJ/Coulomb removed -> spline-only engine (bit-identical).

# ACTIVE PHASE (2026-07-19)

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
- [ ] OPTIONAL (user's call): proper MC/PR barostat for K_A -> kappa; OR temp-sweep stability check; OR
      curve-match to atomistic for an effective time. Not blockers for REMD/metad (equilibrium sampling).

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
      no momentum rescaling needed; Arrhenius gamma is transferable. FIXED main.cpp uninitialized-temp
      call (moved after set_temperature); rebuilt, bilayer bit-identical, hybrid stable. Latent mv-path
      gap flagged (not the active path).

## D. Metadynamics feasibility (task 3)  -- DONE: evaluated (findings 55)
- [x] Core: very feasible (~1-2d, new PotentialNode bias reading a CV node; most CVs exist). Hybrid:
      feasible but membrane geometry CVs must be built + collective-lipid CVs fragile (sub-diffusive).
      PLUMED stub is dead; native well-tempered (spline-grid bias) recommended. Design sketch in findings.

## F. BUG: protein secondary structure melts in hybrid (findings 60)  -- FIXED (prep), re-prep needed
- [x] Root cause: PROTEIN ran at friction_scale=0 (velocity-Verlet, from the Update-51 clamp fix) ->
      no dissipation vs the noisy g-JF lipid bath across the stiff interface -> overheats to kT 1.4-2.0
      (vs target 0.86; lipids fine at 0.86) -> H-bonds melt (1afo 93->58, 1rkl 34->9).
- [x] Fix: prep gives PROTEIN its own Langevin friction, friction_scale=UPSIDE_PROTEIN_FRICTION_SCALE
      (default 100, was 0). Faithful scan: fs=0 kT1.58(melts), fs=100 kT0.87~target (H-bonds plateau ~60,
      folded, still fluctuating = HDX-relevant reversible opening), fs=300 over-cools. Supersedes the
      "velocity-Verlet + OU" protein design (that lacks dissipation).
- [x] HDX correctness preserved: friction is Boltzmann-correct (sets T + effective kinetics, NOT the
      equilibrium); once kT=target the EX2 protection factor K_op is correct. Use min friction that holds
      T (don't over-damp -> would under-sample opening); converge with REMD/metad; water-accessible criterion.
- [ ] 1rkl long verify (running); existing configs must be RE-PREPPED to get the fix (config value, no rebuild).
- ASIDE: long K_A (500k steps, MC barostat) still ~1886 mN/m (undersampled; slow area relax ~ sub-diffusive
  lipids) -> quantitative K_A needs far longer sampling / bigger patch (report as method+bound, not a number).

# Known Errors / Blockers
- g-JF lipid is sub-diffusive (cage-trapped) -> effective-time factor is lag-dependent (findings 47,52);
  the calibration anchors D at the diffusive-regime lag.
- DO NOT break the full-res path or master behavior. Back up any h5 before overwrite. Git READ-ONLY.
