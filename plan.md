# CURRENT PHASE (2026-08-05): NP full-ion rebuild + 6-face stability; glpG micelle REMD half-lost to NaN

## Project Goal

Two tracks running at once. (1) Rebuild the six 1AO6 + MPA-AuNP faces with **0.15 M bulk salt** instead of
counterions-only, because the counterion convention carried no ionic strength at all and leaves the
protein-NP electrostatics unscreened; establish whether all six are stable. (2) Deal with the glpG micelle
REMD batch, half of which died of a single-replica runaway that replica exchange spread to all 48 slots.

## Architecture & Key Decisions

- **NP ions**: `build_ions` now takes `salt_molar` + the solvent-accessible volume and spends the
  neutralizing excess FROM the salt budget rather than adding it on top (adding on top would inflate ionic
  strength by the entire 218-charge MPA+protein excess). `salt_molar = 0` reproduces counterions-only, so
  there is one config value and no dead toggle. Realized per face: **K+ 706 / Cl- 488, net +218**
  (neutralizing protein -15 + MPA -203), solvent 7.82e6 A^3 in the 200 A box.
  *Caveat to settle:* spending from the budget puts the achieved **ionic strength at 0.127 M, not 0.15 M**,
  because the fixed 218 excess is a large fraction of a 706-pair budget. Hitting I = 0.15 M exactly would
  need K+ 815 / Cl- 597. The old 284 A box reached I = 0.142 M for the same reason (bigger budget, same
  excess). User decision which target matters.
- The counterion builds and their trajectories are **preserved, not deleted**: local
  `scratchpad/NP-footprinting/counterion_backup/` (1.0 GB) and cluster
  `/home/yinhanw/project/NP-1AO6/prod_counterion/` (3.2 GB), so the salt-vs-no-salt comparison stays
  possible.
- **Stability cannot be judged by a short check.** `build_all.py`'s built-in probe is 2000 steps; the
  counterion failure appeared at ~step 97,000 of a 112,966-step chunk, and the glpG REMD seed passed a
  400-step smoke test then died at step 12,190. Seed smoke tests bound loading and setup, NOT stability.
  The 6-face probe therefore runs a full chunk-equivalent (112,966 steps at dt=0.005) per face.
- **glpG REMD needs a finiteness guard before any further submission.** `run_remd.py` exchanges on whatever
  energy it reads, so one diverged replica poisons all 48 in a single sweep, and because the driver reseeds
  `/input` from the last output frame the corruption is written back and self-resubmitted. The NP driver
  already does this correctly (`health()` on peptide C-N, ends the chain rather than propagating); the REMD
  driver never got the equivalent.

## Execution Phases

- [x] N1 Full-salt ion path in `np_hybrid.py` + `salt_molar` in `build_all.py`; arithmetic verified.
- [x] N2 All six faces rebuilt: finite, no blow-up, gold_shift 0.00, Rg 26.0-27.2 with drift <= 0.5 A.
- [x] N3 Counterion builds/trajectories preserved; six full-salt seeds uploaded and verified on the cluster
      (5174 atoms, K+ 706 / Cl- 488, net +218, finite); production submitted as **job 53080076**.
- [~] N4 6-face long stability probe: reached ~40-44k of 112,966 steps on ALL six with no tear, then stopped
      because the force-field table was found defective mid-probe. Valid as an old-table baseline only.
- [ ] N5 Report stability per face; decide the ionic-strength target (0.127 M achieved vs 0.15 M nominal).

## Force-field correction (2026-08-05) -- the reason everything is being rebuilt

The spline table was NOT the dry-MARTINI potential: bare LJ + bare 1/r hard-truncated at 1.2 nm, giving a
2.65 E_up (~3.4 kT) energy step at the cutoff for every charged pair. Published dry Martini uses
reaction-field electrostatics (eps_r 15, eps_rf 0) + potential-shifted LJ, both reaching zero at the cutoff.
Fixed in both nonbonded builders in `py/martini_build_tables.py`; verified by
`scratchpad/verify_table_matches_drymartini.py` (old 3.95 E_up deviation and 81/81 rows discontinuous ->
new 2.2e-11 and 0/81). Findings 80. NOT proven to be the NaN cause -- two attempts to show that were
uninformative (one confounded by relaxation, one 10x under-exposed); the powered test is the production
re-run's failure rate against the observed 2/4.

- [x] F1 Corrected both table builders; regenerated `martini.h5` (old backed up as
      `martini.h5.bak.pre-reactionfield`); equivalence asserted against the analytic form.
- [x] F2 All guards removed from the engine per the standing NO GUARDS rule now in AGENTS.md: the
      `main.cpp` blow-up abort and three silent-skip finiteness masks in `martini_potential.cpp`. Rebuilt
      clean locally; user rebuilt the cluster binary (13:08, 0 guard refs in cluster source).
      Consequence: divergences now propagate into logs, exchanges and restarts by design.
- [ ] F3 Rebuild every seed on the corrected table -- 4 glpG variants (`prepared_rf/`) + 6 NP faces -- then
      resubmit both chains. The defective-table builds are kept alongside for provenance.
- [x] G1 glpG micelle REMD: 79ALA (53036667) and 79ALA_S115T (53036669) healthy and left running; 79HIS
      (53036661) and 79HIS_S115T (53036664) NaN, stopped with STOP + scancel (findings 79).
- [ ] G2 Add the finiteness guard to `run_remd.py`, then re-run the two lost variants from the clean seeds
      (still held locally, md5-verified).
- [ ] G3 Localize the runaway: rerun 79HIS from its clean seed with dense output to catch the divergence.

## Known Errors / Blockers

- Two of four glpG micelle variants are lost to NaN. The two survivors carry the same exposure until G2
  lands: they have already shown transient protein_potential spikes of +36,450 and +41,012 that recovered.
- Do not compare the counterion and full-salt NP runs as if the ionic strength were the only difference
  until the 0.127-vs-0.15 M question above is settled.

---

# PREVIOUS PHASE (2026-08-04): rebuild the DDM environment as a micelle, then re-run the glpG REMD

## Project Goal

Replace the DDM *lamellar slab* environment with a physical detergent **micelle**, because DDM is a micelle
under the experimental HDX condition and the slab's 11--14 A hydrophobic core unfolds glpG TM4 (findings 76).
Then re-prepare all four glpG-RKRK variants and re-run the midway3 REMD jobs so the HDX dG profiles rest on a
physically correct environment.

## Architecture & Key Decisions

- **Morphology is derived from the ITP, not a CLI flag and not a name table.** `derive_environment_morphology`
  counts acyl chains as connected components of the apolar (MARTINI `C1`-`C5`) bond subgraph: one tail ->
  micelle, two or more -> bilayer. DDM -> micelle; DOPC/POPC/POPG -> bilayer; a detergent+lipid mixture is
  rejected as ambiguous. No user-facing knob can select the unphysical combination.
- **`build_detergent_micelle()`** replaces tile/crop/carve for micelle morphology. Inputs are all measured,
  never tuned: the oriented protein CG envelope (OPM frame, midplane z=0); the monomer conformer read from the
  existing template PDB; and the hydrophobic half-thickness from the OPM reference REMARK
  (`1/2 of bilayer thickness: 14.1` -> 28.2 A, matching glpG's measured 28.2 A belt span).
- The shell hugs the protein's hydrophobic belt only (tails inward, heads outward), leaving the polar caps
  bare -- a protein-detergent complex, not a closed vesicle. Monomers are seeded by random sequential
  adsorption over the shell VOLUME, **innermost first**, spaced by the force field's contact distance
  (2^(1/6) sigma_max). N therefore follows from belt geometry + the force field: 186 DDM for glpG, in the
  experimental PDC range. Two earlier attempts are recorded so they are not retried: seeding from convex-hull
  support points gave a sparse 32-molecule shell with bald patches, and seeding tail tips in random order
  smeared tails across the shell (63 A span, worse coverage).
- The CHARMM-GUI step5 assembly is **pre-minimization and contains hard clashes** (bead pairs 0.24 A apart), so
  it cannot serve as a packing-density reference. The force field's contact distance is used instead and the
  6.0-6.6 ladder condenses the aggregate.
- **No barostat for micelle morphology** (`npt_enable = 0` at every stage). A micelle in implicit solvent has no
  lateral tension, so the xy-only tensionless barostat is meaningless and would squeeze the box onto the
  micelle. Box = protein+micelle extent plus vacuum padding, fixed.
- Ions: solvent volume = box minus protein minus *measured micelle volume* (not `A_xy x thickness`), placed
  outside the aggregate. `estimate_salt_pairs` takes an excluded volume so both morphologies share one
  formula (bit-identical arithmetic for the bilayer). `place_ions` takes one `reject(trial)` predicate
  instead of `exclude_z`: the z-band for a bilayer, nearer-to-tail-than-head for a micelle.
- `--membrane-thickness-angstrom` is KEPT: it has a real job (count ions against the equilibrated thickness,
  not the compressed CHARMM-GUI start). Only its help text was misleading. It is rejected in micelle mode,
  where the aggregate volume is measured instead. Renaming it would touch the DOPC path for no gain.
- Micelle-frame analysis convention: a finite aggregate has no fixed normal (measured: asphericity fluctuates
  0.19-0.33, long axis drifts 6.4 -> 14.6 deg off box z, while the periodic slab is pinned at 0.21 / <3 deg).
  Depth and tilt must be measured against the aggregate's instantaneous short principal axis, never box z, or
  environment wobble reads as protein bending. HDX dG itself is normal-independent, so the wobble is harmless
  for the deliverable.
- **Gating, revised after measurement.** The build-time span check was DROPPED as a hard gate: a packed-state
  5-95 percentile tail span reads ~20 A for *every* lipid (DOPC 20.4, POPE/POPG 20.6, DDM 11.0) because
  CHARMM-GUI templates are compressed, so comparing it with OPM's relaxed thickness would have falsely failed
  the DOPC path too. What ships instead:
  * `assert_environment_solvation` runs at the production handoff on **equilibrated** coordinates and per belt
    residue -- the only like-for-like comparison. It HARD-FAILS on vacuum (any belt site with no environment
    bead within 2x the contact distance). That catches the diag run's unhealed insertion void (mean 11.7 A).
  * Acyl-tail reach and the **local** equilibrated tail-core thickness are REPORTED with a loud warning, not
    gated, because on a post-damage snapshot they recover: once a helix is stripped, detergent floods the
    cavity. Measured contrast is still stark -- micelle seed 45.1 A local tail core / 0 sites beyond tail
    reach, old lamellar replica 13.9 A / 2 sites, against a 28.2 A belt.
  * Build-time residual bald spots are reported, not fatal: a rigid conformer cannot always reach a deep
    recess, and refining the grid only surfaces more such sites (1.5 A -> 116 exposed, 1.25 A -> 122).
  Promoting the thickness check to a hard gate requires P7 first.
- The bilayer path must keep master parity. The refactor is arithmetically identity-preserving
  (`excluded = A_xy*thickness + V_protein`, so `box_volume - excluded` equals the old
  `A_xy*(box_z - thickness) - V_protein`; the `reject` predicate applies the same z-band test at the same
  point in the RNG stream), but P7 verifies it end-to-end.
- The DDM itp has no `[ position_restraints ]`, so the `lipidhead_fc` ladder is already inert for DDM; the
  6.0-6.6 LJ-softening ladder is retained as the annealing schedule.

## Execution Phases

- [x] P1 Tail-count-derived morphology + `build_detergent_micelle()`; barostat off for a micelle.
- [x] P2 Gates G1 (build-time span) and G2 (post-equilibration solvation) + summary fields.
- [x] P3 Local validation: short 79HIS micelle run -- 186 DDM, r95 43 A micelle, tail span 45 A over the
      28 A belt, first shell 7.04 -> 5.76 A, TM4 helicity 0.97-1.00, box fixed.
- [x] P4 All four variants prepared through stage 7.0; every one passed the solvation gate with 0 bare and
      0 sites beyond acyl-tail reach, local tail core 45.1-47.5 A against the 28.2 A belt.
- [x] P5a Cancel the slab REMD jobs (52970530-533) with STOP files first; NP job untouched.
- [x] P5b Free 74 GB of invalid slab trajectories; keep the slab seeds + hdx_results.
- [x] P5c New cluster base `glpG_DDM_micelle_REMD`; prep code synced to the cluster UPSIDE_HOME.
- [x] P5d Seeds uploaded (md5-verified local vs cluster) and four 48-replica REMD jobs submitted:
      **53036661** 79HIS, **53036664** 79HIS_S115T, **53036667** 79ALA, **53036669** 79ALA_S115T.
      All PENDING; cancelling the running slab jobs forfeited their queue position, so expect a wait.
- [ ] P6 Re-run the HDX analysis; compare per-residue dG against the pD 9 HXMS data.
- [x] P7 DOPC bilayer regression PASSED (rc=0, glpG in DOPC): morphology `bilayer`; xy-only barostat still
      live (box 200.2 -> 213.3 A in xy, z pinned at 123.697); tile/carve ran (64 lipids removed, 15176 lipid
      atoms kept); ions 287 pairs -> 287 Na / 293 Cl, neutralising the +6 protein; solvation gate 0 bare and
      0 beyond tail reach, local tail core 23.0 A vs 28.2 A belt with no warning (correct -- that is normal
      lipid/protein adaptation, not mismatch). Ion-count parity asserted numerically: the refactored
      excluded-volume formula returns 65 pairs for the recorded slab inputs, identical to the pre-refactor
      formula and to the logged `na_added`.

## Known Errors / Blockers

- Every glpG-DDM hybrid HDX dG result from the lamellar prep is invalid (findings 76). Do not compare the old
  and new profiles as if both were meaningful.
- Reaching midway3 needs the user to approve a Duo push; P5 cannot run unattended.

---

# PAUSED PARALLEL TRACK (2026-08-04): NP-1AO6 rebuild with corrected ions, then diagnose the blow-up

## Project Goal

Rebuild the six 1AO6 + 5 nm MPA-AuNP faces with the corrected ion/box prep, verify them through
pre-production locally, emit a PDB of the final configuration for review, then run a production-length
local segment to find out whether the backbone still tears at dt=0.005.

## Architecture & Key Decisions

- Ion prep corrected 2026-08-03 (see `findings.md`): fixed 284 A box for all six faces (was 232--284,
  orientation-dependent) and counterions spent FROM the salt budget as the reference does (was added on
  top). Result per face: 2053 pairs, K+ 2053 / Cl- 1835, 3888 ions, excess +218 cancels MPA+protein.
- `build_system` now asserts achieved ionic strength and exact neutrality from the PLACED ion counts and
  rejects a box too small for the orientation. All three guards negative-tested. This is the check that
  was missing while the ion count was silently re-derived on every rebuild.
- The ion fix is independent of the blow-up. Free salt lands at 0.134 M (I=0.142 M) under the reference
  convention because 218 counterions are a bigger slice of 284^3 than of the reference's 400^3; the
  one-line alternative (counterions on top) gives 4324 ions / I=0.158 M.
- Pre-production = 4-stage soft-core LJ ramp (0.20 -> 0, minimize + 1000 steps MD each) then a 2000-step
  stability MD. Production is dt=0.005, `--integrator v`, thermostat interval -1 (effectively NVE).
- Do not modify or disable SC-env / BB-env or any other physical interaction. Git stays read-only;
  `py/` stays byte-identical (all changes live in `scratchpad/NP-footprinting/`, gitignored).

## Execution Phases

- [x] Correct the ion/box prep and add build-time assertions (`np_hybrid.py`, `build_all.py`).
- [ ] Rebuild all six faces + pre-production locally (`build_all.py`).
- [ ] Emit a PDB of the final post-equilibration configuration for user review.
- [ ] Run a production-length local segment past t~105 to test for backbone tearing.
- [ ] Decide the integrator/dt change based on what that segment shows.

## Known Errors / Blockers

- **Ion convention is unsettled and must be decided.** `np_hybrid.py` currently enforces *counterions only*
  (hard assertion that no salt pair is placed), and the six running faces match it exactly: 218 K+, 0 Cl-,
  since 203 MPA(-1) + protein(-15) = -218. The earlier note in this plan of "2053 pairs, K+ 2053 / Cl- 1835,
  3888 ions at 0.15 M" is a SUPERSEDED convention. These are physically different systems -- counterions-only
  carries no bulk ionic strength -- so for a corona-footprinting study the choice matters. User decision.
- Judge health by Rg / peptide C-N / sign of the potential, never by isfinite: a torn backbone stays finite
  for a long time (the Aug-2 chain resubmitted for 25 h onto already-destroyed systems).
- **Tearing root cause still unresolved, but the gold-LJ-core hypothesis is now excluded** for job 53032366
  (findings 78): the onset was a single ~38-step impulse at res550-551 sitting ~40 A from the nearest gold
  bead, while residues in genuine 3.0-4.3 A gold contact all run never broke. Next step is to instrument
  per-term forces on res549-552 across the onset frame, NOT to lower dt again -- dt 0.009 -> 0.005 plus a
  soft-core ramp only delayed it, and this failure has no steric trigger to blame.
- Job 53032366 ended its own chain correctly and does NOT self-resubmit; 5 of 6 faces finished healthy
  (C-N 1.73-1.83 A). Nothing is queued for this track.

---

# ARCHIVED (2026-07-28): glpG 79HIS REMD HDX dG-by-residue audit — COMPLETE

The July 23 profile was invalid: the source REMD continuation began with NaN energies in all 48 slots,
the HDX extractor had no finite-frame filter, and the MBAR fit did not converge. The plot also used
zero-based donor IDs as residue labels and rendered saturated protection as off-scale +inf spikes joined
to finite points. Exact re-scoring of all 478 T=0.70 VTF frames gives no observed openings for 16 of
residues 30--50, so the helix is strongly protected and the old figure was the artifact. No physics or
force-field settings were changed.

---

# Project Goal

dryMARTINI hybrid membrane in Upside. Full-resolution 14-bead DOPC + protein with a timestep
validated for the coupled hard-core interface. Single-particle CGL subsystem ABANDONED (removed).

REJECTED BASELINE: production mapped a 4 ps relaxation time directly onto `.4` numerical steps
(`tau_up=.0036`). For a DOPC bead (`m=6`) this gave `alpha=1666.7` and
`alpha*dt/(2m)=1.25`, suppressing molecular translation while leaving momentum temperature apparently
normal. The saved stage-7 trajectories are invalid for kinetics.

CURRENT BASELINE: production uses the factor-four-corrected bare-particle mobility clock
(`alpha_bead=0.1691804`) plus additive lipid-contact friction on protein carriers. Both components use one
`.009` g-JF step with no lipid substeps. The molecular DOPC target remains unmatched and is reported.

# Architecture & Key Decisions

- Master reference for diffs: /Users/yinhan/Documents/upside2-md-master (NO martini files ->
  entire martini subsystem is Clean-Slate scope; keep it impeccably clean, remove dead code).
- C++ (kept): martini_potential, dist_spring, angle_spring, martini_masses, martini_brownian,
  martini_fix_rigid, martini_hybrid, martini_stage_params, sc_table/particle tables. CGL removed.
- Integrator: one all-particle g-JF update per `.009` numerical timestep. Production uses the requested
  particle-level fallback: a free MARTINI bead has the lateral mobility corresponding to 40 ps/protein step
  after applying the factor four exactly once (`10 ps` raw MARTINI time and `46 um^2/s` raw target). A protein
  N/CA/C carrier receives the additive FDT friction of the lipid beads inside the existing 12 A spline
  interaction range. Zero-contact carriers use the ordinary Upside OU bath. Positive-friction particles are
  excluded from OU, so no degree of freedom is double-thermostatted.
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

## L. Implement the hybrid-to-existing-HDX compatibility layer

- [x] Add a reusable projector that converts each full hybrid replica into the standard `3*n_res` HDX-view H5
      while retaining full-system potential, temperature, time, and H-bond observables.
- [x] Add an auditable protection-state combiner that preserves stock protein PS and optionally applies an
      externally calibrated per-frame water-accessibility array without changing the existing uptake formulas.
- [x] Update `example/00.AnalysisScripts` driver and steps 2--3 to support ordinary and `martini_hybrid` inputs
      through configuration, with the ordinary path behavior unchanged.
- [x] Add a 16.MARTINI wrapper/configuration for projecting the current stage-7 outputs and running the maintained
      `00.AnalysisScripts` workflow rather than duplicating HDX calculations.
- [x] Update both READMEs, including the correct Upside-unit `T.npy` contract and the scientific limitations of
      dry, uncalibrated membrane accessibility.
- [x] Test projection shape/content, stock protection extraction, downstream single-replica uptake/stability,
      ordinary-path shell syntax, and full-system energy/temperature preservation.
- [x] Reassess HDX trust specifically for unequal lipid/protein kinetic clocks and document the acceptance gates.

### L Decisions

- The compatibility layer is representation-only: mapped N/CA/C positions feed the stock protection extractor,
  while full hybrid potential, H-bond, temperature, and time arrays are copied exactly for downstream weighting.
- A single-temperature wrapper reads `T_up` from replica 0 and analyzes only at that sampled temperature. The
  ordinary workflow retains its existing default when no override is supplied.
- The current 1RKL result is not quantitatively trustworthy. Molecular DOPC diffusion is `0.0152 um^2/s` versus
  `11.5 um^2/s`, the worst first/second-half protection-probability difference is `0.621`, and dry MARTINI has no
  calibrated water accessibility. Unequal clocks do not rescale the EX2 uptake axis, but here they accompany
  demonstrably inadequate ensemble mixing.

# ACTIVE PHASE (2026-07-20)

## O. Consolidate the confirmed BB-routing and rigid-handoff design in the manuscript

- [x] Re-audit the current 1RKL artifact and extract only result claims supported by the corrected construct.
- [x] Rewrite the manuscript as one coherent design narrative: regenerated-O semantics, partial BB force routing,
      rigid pre-production handoff, explicit release to flexible production, and unchanged SC-env/BB-env physics.
- [x] Remove obsolete candidate-fix, unprotected-handoff, and pending-correction language throughout the TeX.
- [x] Preserve the separate transport and HDX trust boundaries; structural repair must not be presented as a
      calibration of molecular DOPC kinetics or membrane water accessibility.
- [x] Compile twice, inspect the rendered text for warnings and contradictions, and run repository diff checks.

### O Decisions

- The confirmed interface construct is the manuscript baseline, not an appended workaround: forward O/BB
  reconstruction, 14/54--12/54--12/54 reverse routing to persistent N/CA/C, and disposal of the regenerated-O
  share under the Upside integration cycle.
- The protein is one rigid body throughout every pre-production segment. BB-env and SC-env are active during the
  dedicated `production_handoff`; flexibility begins only at the explicit `production` transition.
- Resolving secondary-structure loss validates this force-routing and handoff construct. It does not by itself
  validate the 40 ps kinetic mapping, DOPC molecular diffusion, HDX populations, or membrane accessibility.

### O Review

- The manuscript now derives the forward O/BB construction and the persistent reverse sensitivities as explicit
  equations, then describes one rigid pre-production protocol through the flexible-production transition.
- Current corrected artifacts contain 1,001 frames at `T_up=.8`. DSSP gives 1RKL first/final/minimum/mean helical
  counts of `23/23/16/21.88`; residues 10--28 have 96.5--100% helical occupancy. The corresponding 1AFO counts
  are `54/51/48/50.77`.
- Two-pass `pdflatex`, warning/contradiction scans, `chktex`, and `git diff --check` pass.

## N. Rewrite friction calibration and evidence boundaries in the interface manuscript (completed)

- [x] Trace the implemented temperature, diffusion, friction, contact-count, and thermostat equations from the
      workflow/H5 builder so every manuscript formula matches code.
- [x] Replace obsolete or contradictory timescale language with one derivation from the factor-four-corrected
      bare-particle target to the g-JF friction and contact-local protein drag.
- [x] Rewrite the results, HDX, discussion, and conclusion claims into explicit validated, conditional, and
      unsupported domains, including current molecular diffusion and ensemble-convergence failures.
- [x] Remove duplicated patch-like prose and ensure the document reads as one methods narrative.
- [x] Compile the TeX and inspect warnings/errors; run repository diff checks.
- [x] Recompute trust-critical diagnostics after the fresh unified-temperature stage-7 outputs replaced the
      previously audited trajectories, then update all result-specific manuscript statements.

### N Decisions

- One system `TEMPERATURE` controls both runtime FDT noise and the temperature used to calibrate friction.
- The calibrated observable is free/bare MARTINI-particle mobility. Whole-DOPC molecular diffusion and a full
  membrane hydrodynamic resistance tensor are not calibrated.
- Correct thermostat statistics and a mechanically working HDX adapter are necessary but insufficient for a
  trustworthy HDX prediction. Current results remain qualitative until timestep, ensemble, replica/enhanced-
  sampling, and membrane-water-accessibility gates pass.

### N Review

- The manuscript derives $D_\mathrm{bead,up}=5.11111$ and
  $\alpha_\mathrm{bead}=0.15652$ at the unified default $T_\mathrm{up}=0.8$, documents the static contact-count
  approximation, and distinguishes preparation damping from production calibration.
- An evidence table now separates implementation checks, calibrated inputs, failed molecular transport, open
  equilibrium populations, and the mechanically correct but scientifically unvalidated HDX adapter.
- Two-pass pdflatex completes without warnings; arithmetic assertions, manuscript consistency checks, chktex,
  and git diff --check pass.
- The pre-fix 1RKL/1AFO artifacts stored runtime and calibration $T_\mathrm{up}=.8$, the derived
  $\alpha_\mathrm{bead}=.15652$, and 1,001 frames. Kinetic temperature passes, molecular diffusion remains more
  than $10^3$-fold too slow, and their structural result is superseded by the corrected-interface runs in O.

## M. Diagnose the 1RKL BB-env force regression against `b1041bb`

- [x] Check the cited stage-7 VTF/H5 provenance and identify affected residues rather than relying on the final
      frame. The cited output is an older `T_up=0.8647`, `dt=.009` run, not the requested `T_up=.80` run.
- [x] Diff the particle--backbone table generation, H5 schema/wiring, C++ pair-force evaluation, and virtual-BB
      derivative propagation against commit `b1041bb6640b34edf55fe371b5cf943d054ba187`.
- [x] Compare current and historical BB-env potential/force curves in common Upside units and test force sign,
      magnitude, residue/type assignment, pair exclusions, and Newton-pair reactions on representative inputs.
- [x] Isolate the first behavior-changing commit or code path and distinguish a force regression from timestep,
      thermostat, startup-hold, or trajectory-provenance effects.
- [x] Re-audit the fresh `T_up=.80` H5 invocation, temperature/friction metadata, and residue-wise DSSP
      timeline. The unified temperature is correct, but 1RKL still has substantial helix loss.
- [x] Compare the BB--particle interaction definition across the coarse and full-resolution historical paths:
      interaction centre, type assignment, radial table/cutoff treatment, pair admission, and virtual-site force
      projection.
- [x] Determine the Upside-correct BB--particle reverse force mapping under regenerated-O integration semantics;
      distinguish the persistent N/CA/C force from the force written to disposable O.
- [x] Restore the historical partial BB routing (14/54 N, 12/54 CA, 12/54 C; disposable 16/54 O share omitted)
      directly in the current coordinate-node architecture, while retaining the protected stage-7 handoff and
      never disabling SC-env or BB-env.
- [x] Make the entire stage-7 interface handoff conformationally rigid: activate BB-env/SC-env in a dedicated
      rigid `production_handoff` stage, equilibrate the membrane, and switch to flexible `production` only after
      burn-in. Remove the superseded stage-7 positional-spring/release controls.
- [x] Verify partial-routing source/H5 weights, stage-6.6-to-handoff rigid geometry, shell/Python syntax, full
      build, manuscript compilation, and repository diff.
- [x] Regenerate the complete 1RKL stage-7 trajectory with the corrected handoff and repeat the DSSP/timestep
      audit; the corrected production trajectory preserves the 1RKL helix.

### M Decisions

- There is one authoritative `TEMPERATURE` for the Upside thermostat and DOPC friction calibration.
- Raw particle--BB tables, type assignment, derivative sign/scale, and pair classes match `b1041bb`; the material
  regression was the reverse route. Persistent N/CA/C receive 14/54, 12/54, and 12/54 of BB sensitivity, while
  the regenerated-O share is disposable under the Upside integration cycle.
- The complete protein is one rigid body in every pre-production stage, including interface-active
  `production_handoff`; only the explicit `production` transition starts flexible protein dynamics.
- The standard numerical step remains `.009`. This structural correction does not alter the separate transport
  limitation or validate a 40 ps molecular clock.

### M Review

- The reverse pass now applies H5 weights 14/54, 12/54, and 12/54 to persistent N/CA/C. The stored 16/54 O
  weight remains part of the forward BB definition but its derived-site sensitivity is discarded.
- Full corrected 1RKL and 1AFO workflows exercised interface-active rigid handoff followed by flexible production.
  Maximum persistent-carrier pair-distance changes through the handoff were `0.017 A` and `0.0037 A`, respectively.
- Corrected DSSP results resolve the regression: 1RKL begins and ends production with 23 helical residues, while
  residues 10--28 remain helical in 96.5--100% of frames.
- C++ build, Python compilation, shell syntax, removed-control checks, two-pass TeX compilation, and
  `git diff --check` pass.

# PREVIOUSLY COMPLETED PHASE (2026-07-20)

## K. Reassess 1RKL temperature/secondary structure and define a hybrid HDX workflow

- [x] Read the production invocation and H5 temperature/friction metadata and convert the requested Upside
      temperature to Kelvin.
- [x] Measure protein/lipid kinetic temperatures independently and determine whether the new integrator needs
      a distinct temperature conversion.
- [x] Quantify residue-level DSSP occupancy, transient helix loss, recovery, and backbone RMSD from the final
      stage-7 trajectory rather than relying on endpoint counts.
- [x] Reconcile the result with prior timestep-resolution evidence and the temperature used by ordinary Upside
      membrane examples.
- [x] Specify the coordinate, protection-state, ensemble, and validation adapters needed to use full hybrid
      trajectories in the existing HDX analysis contract.

### K Decisions

- Production is at `T_up=0.8647 = 303.15 K`. Protein and lipid last-window kinetic energies are respectively
  1.1% and 0.3% above `3 kT/2`; there is no missing integrator-specific temperature scaling. Never lower the
  reported temperature merely to compensate for structural loss. Use `T_up=T_exp/350.588235` for HDX.
- DSSP confirms a real but reversible hydrogen-bond-network disruption: the residue 10--29 core averages 84.5%
  helix occupancy, falls to four helical residues near 1.09 us, and ends at 19/20; its CA RMSD stays below
  1.58 A. The central experimental alpha2 region remains the most stable. Endpoint-only validation previously
  hid this event.
- Correct kinetic temperature does not establish configurational timestep convergence for nonlinear hard-core
  forces. Existing matched-duration evidence found `.009` under-resolves the coupled interface while `.00225`
  preserves the helix. The next controlled test is a current-model `.009/.0045/.00225` comparison at the same
  physical temperature and duration. Do not relabel this as a temperature-conversion problem.
- A mechanically compatible HDX adapter should evaluate mapped hybrid N/CA/C coordinates with a separate
  protein-only `-HDX.up` engine, reuse the standard H-bond/side-chain/burial terms, and combine them with full-
  system membrane accessibility. The current global-PO4 binary slab is not quantitative because it saturates
  stable transmembrane donors.
- Hybrid MBAR must use the full coupled-system `output/potential`. Its `T.npy` must remain in Upside `kT` units
  (for example 0.8647), despite the current README incorrectly describing Kelvin; `4.calc_D_uptake.py` uses
  `beta=1/T` in Upside units and converts to Kelvin only for intrinsic chemistry.
- Dry MARTINI contains no water, so quantitative membrane accessibility must be supplied by a calibrated local
  depth/water-activity model or an explicit-water reference. Until that is validated, report H-bond/protein-
  burial HDX separately from the membrane-accessibility correction.
- Revised per user direction: do not create a parallel hybrid HDX implementation. Add one projection adapter
  that writes a protein-only HDX-view `.up` with mapped N/CA/C positions but copied full-hybrid potential,
  temperature, and H-bond logs. Then reuse `example/04.HDX/1.config.py`, the stock steps 2--3, and
  `example/00.AnalysisScripts` steps 4--6 unchanged. Apply any calibrated membrane-accessibility correction as
  an explicit postprocessing layer on the stock `PS.npy`, preserving the unmodified protein-only PS for audit.

# PREVIOUSLY COMPLETED PHASE (2026-07-20)

## J. Audit the HDX analysis workflow and reassess trust

- [x] Trace protection-state extraction from trajectory coordinates and identify whether it is an equilibrium
      population estimator, a kinetic opening estimator, or a structural proxy.
- [x] Trace protection factors, intrinsic exchange rates, peptide uptake, temperature mapping, and experimental
      comparison, including all fitted or empirical parameters.
- [x] Determine which outputs are invalidated by frozen trajectories and which could remain mathematically
      defined but statistically untrustworthy.
- [x] Run targeted synthetic/current-output checks where needed and issue a method-specific trust assessment.
- [x] Record conclusions and limitations in findings/progress without changing the HDX method unless a code
      defect is independently demonstrated.

### J Decisions

- The uptake calculation is an equilibrium EX2-like structural proxy, not direct MD-time kinetics. It averages a
  binary protected state, uses `k_obs = k_chem * (1 - p_protected)`, and evaluates uptake analytically in seconds.
  The 40 ps/Upside-time mapping therefore affects sampling efficiency/convergence, not the uptake time axis.
- The stock extractor and trajectory-observable scripts are incompatible with hybrid files because they assume
  exactly three coordinates per residue. Hybrid stage-7 inputs contain 2,894 (1AFO) and 4,098 (1RKL) particles,
  versus protein-only expectations of 216 and 93. The hybrid-aware extractor is not wired into `analysis.sh`.
- The available hybrid-aware membrane criterion is a geometric prior: any donor between global phosphate planes
  is protected. On stage-6.6 it makes all 29 1RKL donors protected in all 11 frames despite 65.5% H-bond-open
  observations; the downstream `mean_pf == 1` branch then imposes `k_obs = k_chem/1000`. This is saturation,
  not a resolved protection factor.
- No HDX/PS/percentD artifacts exist under `example/16.MARTINI`, so there is no produced hybrid HDX result to
  validate. Existing frozen outputs are unusable for ensemble estimates; corrected single 50,000-step runs span
  2 us on the stated 40 ps-per-numerical-step protein clock but still need convergence/replica evidence.
- Do not modify the HDX model during this audit. A quantitative hybrid workflow first needs an explicit mapped-
  coordinate path, robust membrane/water-accessibility treatment, and convergence/replica validation.

# PREVIOUSLY COMPLETED PHASE (2026-07-20)

## I. Resolve production freeze introduced by the friction-clock mapping

- [x] Quantify protein and DOPC coordinate mobility in the saved 1AFO and 1RKL stage-7 VTF/H5 trajectories.
- [x] Derive the actual coordinate mobility implied by the current mass-scaled g-JF friction and identify the
      mismatch between relaxation-clock metadata, numerical integration time, and physical 40 ps/step time.
- [x] Replace the freezing friction treatment with the particle-mobility/additive-contact model without
      changing SC-env, BB-env, hard cores, or spline tables.
- [x] Re-run fresh 1AFO and 1RKL workflows and gate on coordinate motion, DOPC COM MSD, secondary structure,
      kinetic temperature, virtual-site constraints, and build/static checks.
- [x] Update the manuscript and calibration claims to match the verified observable.

### I Decisions

- Do not restore the native-relaxation or earlier `friction_scale=100` regimes: both unfreeze DOPC but heat
  the interface and reduce 1RKL helicity in matched 5,000-step tests.
- Production bare-bead friction is derived from the stated observable, not fitted:
  `D_raw = 4 D_target`, `dt_raw = 40 ps / 4`,
  `D_up = D_raw * 1e-4 A^2/ps * dt_raw / .009`, and `alpha_bead = kT/D_up`.
  At the defaults this is `alpha_bead=0.1691804`.
- Protein drag is local and additive: `alpha_protein,i = n_contact,i * alpha_bead`, where `n_contact,i`
  counts DOPC beads within the already-defined 12 A spline cutoff in the production starting structure.
  This is a static, one-body approximation to interfacial drag, not a full hydrodynamic resistance tensor.
- The acceptance gate includes coordinate-space mobility and DSSP, not kinetic temperature alone. The
  5,000-step prototypes moved DOPC COM by 1.8--2.0 A RMS and retained 54/54 (1AFO) and 19/24 (1RKL) final
  helical residues, with last-window protein kinetic energy within 2% of `3kT/2`.
- A wiring-only fresh-workflow shortcut that disabled all minimization and used 50-step equilibration stages
  was invalid: it entered stage 6.0 at potential `+131697`, overheated immediately, and exploded when the
  full hard core activated. It is not a model result. Fresh-workflow verification must retain the normal
  minimization and 500-step soft-core schedule before testing the new production friction.
- Revised after fresh 1AFO exposed a real pre-production failure: stage 6.3 injects large kinetic energy as
  the core is sharpened, and the production-mobility bath is too weak to settle those overlaps before the
  full hard core. Use three explicit phases: native dry-MARTINI damping in softened early equilibration;
  strong FDT overlap-settling damping in full-core stages 6.1/6.4--6.6; and the derived particle/additive-
  contact clock only in production. Preparation damping is not a kinetic claim or analyzed trajectory.
  Conservative softening, restraints, and all interface interactions remain unchanged.

# PREVIOUSLY COMPLETED PHASE (2026-07-20)

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
  all-particle g-JF architecture, which applies SC-env/BB-env reactions through the hybrid coordinate node.
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
  freedom and must not be integrated or thermostatted. The coordinate node constructs both sites without
  mutating integration state, while its reverse pass follows the regenerated-O contract by routing the
  14/54, 12/54, and 12/54 BB shares to persistent N/CA/C and discarding the O share.
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

- The user-visible freeze is removed: 50,000-step DOPC COM net RMS motion is `3.43 A` (1AFO) and `3.72 A`
  (1RKL), versus `0.17--0.24 A` in the rejected outputs.
- Total kinetic energy is within 1% of target; 1AFO retains all 54 helical residues, and the physical 1RKL
  transmembrane core at residues 10--29 retains mostly 0.9--0.99 helix occupancy.
- Measured molecular DOPC diffusion is `0.0132--0.0152 um^2/s`, not 11.5. The verified claim is the
  factor-four-corrected bare-particle mobility plus an additive local protein-friction approximation.
- A valid fresh workflow completed minimization, equilibration, burn-in, production, and VTF extraction with
  finite moving coordinates and refreshed contact metadata for both 1AFO and 1RKL. Fresh 5,000-step DOPC
  COM net RMS motion was `1.45 A` and `1.22 A`, respectively.
- Build, Python compilation, shell syntax, diff whitespace, H5 metadata, fresh workflow, 50,000-step
  regressions, DSSP/mobility analysis, and manuscript compilation pass.

# Known Errors / Blockers
- The CHARMM-GUI DDM *lamellar* template cannot solvate glpG: its tail core is 11-14 A against a 28-30 A
  protein TM belt, which unfolds TM4 through the TM4:TM6 interface. DDM must be built as a micelle (or the
  lamellar environment must use a bilayer-forming lipid). All glpG-DDM hybrid HDX dG results produced from the
  redesigned CHARMM-GUI prep are invalid until this is fixed (findings 76).
- Molecular DOPC diffusion is not matched at 40 ps/step; the fallback is explicitly particle-level.
- The additive contact friction is not a quantitative hydrodynamic resistance tensor. Such a tensor
  would require independent parallel, normal, and rotational reference data and a momentum-conserving model.
- DO NOT break the full-res path or master behavior. Back up any h5 before overwrite. Git READ-ONLY.
