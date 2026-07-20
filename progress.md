# Progress

## P11 DONE -- unit conversion baked into h5, ZERO C++ conversion (2026-07-18)

Moved native dry-MARTINI -> Upside unit conversion out of the C++ runtime into the Python h5-build.
Runtime h5/config now stores Upside-unit values; C++ does no conversion.

Edits (all landed):
- py/martini_build_tables.py `_build_sc_table_group`: convert+rename ONLY at final create_dataset
  (grid_ang = grid_nm*10; *_energy_eup = *_energy_kj_mol/2.914952774272). Internal arrays stay native.
  angular_profile/rotamer_angular_profile/cos_theta_grid dimensionless -> unchanged.
- py/martini_prepare_system_lib.py: main node writes single `coulomb_constant = native*(L/E)` attr
  (dropped 3 conversion attrs); softening builder reads baked `coulomb_constant`; inject_stage7 reads/
  writes new *_eup/grid_ang names, dropped conversion-attr writes on the sc_table node.
- py/martini_prepare_system.py: SC_LIBRARY_REQUIRED_DATASETS renamed to new names.
- src/martini_potential.cpp (3 nodes): read grid_ang + *_energy_eup + coulomb_constant directly;
  removed all conversion members/reads and the /E, *L math (tail-subtraction shift kept). NOTE coulomb_k
  is read-but-unused in compute (Coulomb baked into combined_energy_grids); attr validates presence.
- Regenerated parameters/ff_2.1/martini.h5 (backup .bak, no version number). Rebuilt obj/ (make exit 0).

Gates: A py_compile+make 0; B bake exact (maxabs 0.0); C bilayer BIT-IDENTICAL -11039.306641
(|diff|=0.000000, old binary confirmed same baseline first); D hybrid rc==0, no NaN, potential[0]=1.9e20
(expected r^-12, not gated). See plan P11, findings Update 47.

## Route C BUILT -- g-JF single-step lipids + effective-time factor (2026-07-18)

User chose Route C (fast/elegant/no-multi-step, approximate D via a MARTINI-style effective-time factor)
after the literature (findings 39-45) closed exact-D single-step. Build:
- g-JF inertial-Langevin single-step lipid integrator: martini_brownian.cpp, gated /input/brownian attr
  langevin (default 0 = overdamped RESPA). Already built + verified (stable, Tkin=1.00, structure intact).
- Prep wiring (py/martini_prepare_system_lib.py): writes langevin (env UPSIDE_LIPID_LANGEVIN, default 0)
  + effective_time_factor (env UPSIDE_LIPID_TIME_FACTOR, default 1.0). Defaults keep overdamped ->
  master parity. Convention: physical_time = sim_time * effective_time_factor (applied in analysis).
- HYBRID go/no-go PASSED: md_1rkl.up with g-JF lipids + Verlet protein, 1000 steps -> stable, bilayer
  intact (frac 0.95, thick 39), protein embedded (protZ 0.4), Tkin 1.06. Inertia handles the stiff
  protein interface single-step (which forced M=90 for the overdamped lipids).
- Calibration: effective-time factor R ~ 7.3 (D_overdamped/D_gjf, MSD@27; first pass R(T)=7.36/6.29/8.23,
  T-const within noise); 2nd-seed refinement running.
- Enable: UPSIDE_LIPID_LANGEVIN=1 UPSIDE_LIPID_TIME_FACTOR=7.3. Caveat: single global factor is a
  MARTINI-class approximation of the protein<->lipid relative clock (protein core stays Verlet).
See plan P9, findings 39-45.

## g-JF inertial-Langevin single-step lipids -- prototyped; STABLE but D too slow (2026-07-18)

Implemented g-JF (Gronbech-Jensen & Farago) inertial Langevin for the lipids in martini_brownian.cpp,
gated by /input/brownian attr "langevin" (default 0 = overdamped; existing runs bit-identical, master
parity). Motivated by the literature review (findings Update 39) as the route off the M=90 RESPA
sub-stepping. Bilayer test, single step dt=0.009 (no sub-stepping): STABLE (potential cohesive,
finite), correct thermostat (Tkin=1.00 kT), bilayer structure preserved (frac 1.00, thickness ~30) --
so single-step against the hard core WORKS. BUT the lipid diffusion is ~8-10x too slow at the reused
overdamped-calibrated friction (matched-time MSD): the real CG mass makes the inertial regime
UNDERDAMPED (tau=m/zeta~57 ps), and Kramers cage-escape there is slower than the overdamped
Smoluchowski dynamics -- the overdamped zeta does not transfer. Correct D would need a dedicated
inertial-regime friction re-calibration (+ Arrhenius D(T) re-fit); not a drop-in. Production stays on
overdamped M=90. g-JF kept as a gated prototype. See findings Update 40, plan P8.

CALIBRATION EXHAUSTED (2026-07-18): friction scan (x0.3..x30) -> MSD friction-insensitive (cage-trapped,
~10x low); terminal D (1800 Tup) = 0.0055 (~50x low); mass scan -> mass is a knob but the single-step
stability wall hits at mass x0.0003 (blow-up), lightest stable x0.001 -> ~1.5x low; mass x0.001 +
friction tweak peaks ~0.8x target, never reaches it. So g-JF cannot reach the calibrated D single-step
(correct D needs the overdamped/light-mass limit, which reintroduces the core-overshoot instability --
fundamental tension). Overdamped M=90 stays production. See findings Update 41.

## min-M determined = 45 (2x speedup) via cap-firing scan (2026-07-18)

Determined the minimum stable lipid sub-step count M. Signal = displacement-cap firing rate
(temporary counter in martini_brownian.cpp, since removed) -- per-bead, immune to the undulation
confounds that ruined earlier scans. Built + equilibrated a dense 128-lipid DOPC bilayer (SPACING
8.1 ~0.66 nm^2/lipid) at M=90 (cap 6e-6, intact, potential -11136), then scanned M from the
equilibrated state: cap negligible + energy/thickness identical to M=90 down to M=45 (cap 3e-6);
sharp knee at M=30 (cap 1.8e-3, thinning + cohesion loss begins); collapses below. min-M = 45
(dt_inner 2e-4) -> 2x. M is a stability knob (D is M-independent) so 45 preserves physics. Temporary
counter removed; production binary bit-identical to pre-diagnostic. See findings Update 37.

VERIFIED ON THE HYBRID (2026-07-18): min-M=45 does NOT transfer. On an equilibrated 1RKL hybrid
(raw BB-env, dt=0.009) the protein interface makes lipids cap-fire ~1.7% even at M=90, and reducing
M degrades monotonically (M=75 halves cohesion; M=45 disperses + potential positive). Hybrid
min-M ~ 90 -> PRODUCTION DEFAULT STAYS M=90 (no change). Bulk lipids are over-resolved at 90 but the
protein interface is the binding constraint. See findings Update 38. (Cap counter was re-added for
this check, then removed; binary bit-identical to production.)

## P6 level-0 lipid-only inner refresh — INVESTIGATED, INEFFECTIVE, REVERTED (2026-07-18)

Tried skipping the protein FF in the lipid Brownian inner substeps: compute(DerivMode,0) for
substeps 1..M-2 (protein fixed during substeps), full compute on the last substep to keep the
protein Verlet update exact + build_integrator_levels for the brownian path (src/martini_brownian.cpp,
src/main.cpp). Proven BIT-IDENTICAL (1RKL hybrid, old vs new, identical incl NaN onset) but
differential timing = 0.99x (352 ms/outer-step both). Inner loop is dominated by martini_potential
(lipid nonbonded x M=90); protein FF is negligible for an embedded TM helix. REVERTED both edits;
rebuilt binary bit-identical to pre-change. Real speed lever = reduce M (min-M). See findings Update 36.

## All-atom 1RKL tilt reference moved to midway3 GPU (2026-07-18)

Offloaded the CHARMM-GUI 8430021493 GROMACS job to midway3 (job 52367423, gpu partition, Tesla
V100, 316 ns/day). Re-grompp'd on GROMACS 2025.3 from the text inputs (2026 .tpr/.cpt are
version-locked); 100 ns production at /project/trsosnic/yinhanw/allatom_1rkl_8430021493. Local Mac
run stopped (superseded). Connection via ssh ControlMaster (one Duo approval); ops notes in
~/.codex/AGENTS.md.

## Cell-list neighbor search for martini_potential — DONE (2026-07-18)

Replaced the O(N^2) pairlist rebuild scan with an O(N) cell list producing the IDENTICAL
active set; lowered default cache_buffer 10.0 -> 2.0.

Files: src/martini_potential.cpp (MartiniPotential: CSR adjacency + cache-friendly cell-list
rebuild); py/martini_prepare_system_lib.py (cache_buffer 10.0 -> 2.0).

Design:
- Construction: build_pair_adjacency() = CSR over non-zero-coeff pairs keyed by lower atom
  (row sorted by neighbor since pairs are stored first-asc/second-asc) -> binary-search lookup.
- rebuild: counting-sort atoms into a cell grid (edge >= active_cutoff) with cell-sorted
  position scratch (sequential/cache-friendly reads); PBC via wrapped coords + mod-wrapped
  27-cell stencil (neighbor_cells dedups collapsed cells for small grids); non-periodic uses
  coord-extent grid + clamped stencil. Per candidate: dedup by b>a, distance filter
  (minimum_image_scalar, bit-identical to scan), then lookup np. Sort active list ascending
  -> identical accumulation order -> bit-for-bit identical energy/forces.

Verify (ALL PASS):
- Identical active set: temporary in-code assert (cell set == old scan set, incl. skin pairs)
  = 0 mismatches over 900+ substeps, cb2 & cb10, periodic AND non-periodic paths.
- Bit-identical physics: full potential trajectory old-scan vs new-cell IDENTICAL at matched
  cache_buffer (seed 7), and cb2-vs-cb10 IDENTICAL (skin-independent), 120-step run identical.
- Stability: M=90 bilayer dt=0.009 integrator v, 120 steps finite/stable, membrane intact.
- Speed (contention-robust ratios; GROMACS running): M=90 (1792 atoms) new(cell,cb2) vs
  old(scan,cb10) = 1.28x. Rebuild is O(N): 4170 us@1792 -> 19026 us@7168. At 7168 atoms cb2
  (same compute), cell vs scan = 1.7x on the step (rebuild-dominated); old cb10 unrunnable at
  that size. Small-box gain is capped by the Brownian+spline compute floor (rebuild was never
  the 1792-atom bottleneck); the O(N) win grows with N. cache_buffer=2 is near-optimal
  (cb1 over-rebuilds; cb3-4 ~5% better; kept 2.0 for a lean list per request).

## Finalize: raw BB-env + standard timestep — DONE (2026-07-17)

Reverted the earlier BB-env PMF (it drove the embedded protein OUT of the bilayer; the RAW
dry-MARTINI FF confines it correctly) and moved to the Upside standard timestep.

(A) py/martini_prepare_system_lib.py — BB-env back to RAW:
- convert_stage: removed is_bb_proxy/is_lipid_bead flags, unique_coeff_is_bb_env, the is_bb_env
  coeff-key arg (back to 4-tuple (eps,sig,q1,q2)), the per-pair bb_env grouping (back to grouping
  by class_j only), and the bb_env_coefficient_mask dataset write.
- Removed bb_env_pmf_grid() + BB_ENV_PMF_QUAD_ORDER.
- inject_particles_table: removed the PMF branch + UPSIDE_BB_ENV_RMSF; every coefficient row now
  gets the raw matched LJ+Coulomb grid. C++ bb_env_interface_potential (energy reporting) untouched.

(B) Standard timestep + matched sub-stepping (dt_inner ~ 1e-4 -> M=90 at dt=0.009):
- run_sim_hybrid.sh: EQ_TIME_STEP 0.010->0.009, PROD_TIME_STEP 0.004->0.009, UPSIDE_LIPID_NSUBSTEP
  10->90 (+ comment fix).
- run.py: eq_time_step / prod_time_step -> 0.009.
- prepare_system_lib.py: n_substep env default 10->90 (run.py doesn't set the env var, so this keeps
  it consistent). run_sim_1rkl_full.sh / run_sim_1afo_full.sh are thin wrappers -> no timestep to change.

Verify (all PASS):
- Repo grep clean of is_bb_env/BB_ENV_RMSF/bb_env_pmf_grid/bb_env_coefficient_mask and of any bb_env
  in the dryMARTINI interface (py/martini_* + example/16.MARTINI). Remaining bb_env = core Upside
  bb_env.dat feature (unrelated) + C++ energy reporting.
- 1RKL hybrid prep OK: 22 coeff rows; /input/brownian = all 3850 LIPID beads, 0 protein, n_substep=90.
- BB-proxy<->lipid coeff grid == raw analytic LJ+Coulomb (charged row @3A = 14748 raw, NOT the PMF's
  -0.93); matches to float32 storage precision (max|d| only at r->0 core where LJ~1e12).
- 100-step MD dt=0.009 integrator v: finite/stable (initial clash 2.4e13 settles in ~10 steps to
  bounded ~-1e4; Rg 12.8; no NaN; exit 0).
- Pure DOPC bilayer: combined_energy_grids + coefficients + coefficient_indices BYTEWISE identical
  new-vs-old(PMF) code (reconstructed old lib in scratchpad). Master parity preserved.

## CGL Clean-Slate removal + FF/workflow consolidation — DONE (2026-07-17)

Removed the abandoned single-particle CGL (coarse-grained lipid) subsystem end to end.
KEEP path = full-resolution 14-bead DOPC + overdamped martini_brownian integrator.

C++ (src/):
- Deleted martini_cg_lipid.cpp; removed `namespace martini_cg_lipid` from martini.h
  (kept `#include "martini_brownian.h"`); removed CGL random streams from random.h
  (kept BROWNIAN_RANDOM_STREAM=7); removed the integrate_dynamic_orientation/compaction
  calls from all 3 integration_cycle overloads in deriv_engine.cpp (kept martini_brownian);
  removed is_cg_lipid_term_name + all register/set/apply/add/clear CGL wiring from main.cpp
  (kept martini_brownian); dropped martini_cg_lipid.cpp from both CMakeLists.
- `make upside` PASSES.

Python (py/):
- martini_cg_lipid_params.py DELETED; _pair_param + dopc_max_sigma_nm relocated into
  martini_itp_reader.py; importer fixed in prepare_system.py.
- prepare_system_lib.py: removed all coarse branches + CGL helpers/validators/build_cg_lipid_array/
  inject_cg_lipid_nodes; convert_stage is now unconditional full-res; kept /input/brownian
  injection, inject_particles_table, inject_stage7_sc_table_nodes, dist/angle/martini_potential.
- prepare_system.py: removed apply_cgl_transport_metadata, promote_cgl_*, cgl_wrapped_z_outliers,
  the --cg-lipid-*/--cgl-gle-*/--dopc-h5/--lipid-resolution args, all coarse branches; rewired to
  read parameters/ff_2.1/martini.h5 + itp from example/16.MARTINI/dryMARTINI_itp.
- martini_build_tables.py: removed build_dopc_h5 + 162 CGL/DOPC-CG derivation functions
  (13960->1113 lines); added build_martini_h5 (writes /particles + /sc_table in one file);
  build_particle_h5/build_sidechain_h5 unchanged.
- martini_gen_params.py: rewritten to build the single ff_2.1/martini.h5 via build_martini_h5.
- martini_extract_vtf.py: removed CGL rod-display fns + args; core VTF exporter intact.

H5 / ITP:
- parameters/ff_2.1/martini.h5 created (merged existing /particles + /sc_table). dopc.h5 + all
  its .bak deleted (backed up to scratchpad). Superseded dryMARTINI particle/sidechain h5
  removed (.bak retained).
- *.itp moved to example/16.MARTINI/dryMARTINI_itp/. DOPC.pdb stays in parameters/dryMARTINI/.

MD scripts / docs:
- run_sim_hybrid.sh full-only (single martini.h5 check, itp path). Deleted run_sim_1rkl.sh,
  run_sim_1afo.sh, run_sim_1rkl_outlipid.sh, run_sim_1afo_outlipid.sh. Kept + de-CGL'd
  run_sim_1rkl_full.sh, run_sim_1afo_full.sh. run.py + build_martini_h5_{m1,slurm}.sh updated.
- Deleted cg_lipid_potentials.tex(+aux/pdf/...), analyze_cgl_bilayer.py, analyze_hybrid_validation.py.
  readme.md trimmed.

Verification:
- Rebuild PASSED. All KEEP files grep-clean of cg_lipid/cgl/compose_vector6d/coarse/lipid_resolution.
- Full-res bilayer prep+run: 72 DOPC (1008 lipid beads on /input/brownian), potential
  32853 -> -193 E_up over 30 steps, stable, no NaN. build_particle_h5 smoke test OK with new itp.

Not changed (flagged): AGENTS.md line 16 keeps a CGL *principle* statement (governance doc,
analogous to CLAUDE.md) — left intentionally, not CGL implementation.

## 2026-07-18: finalize method (g-JF single-step) + master-based cleanup

TEX: drymartini_upside_interface.tex fully rewritten for the FINAL method -- g-JF single-step lipids
(Boltzmann-correct statics + hard core at the coarse step) + effective-time factor R~7.3
(t_phys = t_sim / R; g-JF under-diffuses, inverse sign of MARTINI factor-4). Overdamped RESPA kept as
the exact-D calibration reference (design-space corner i). Honest scope: rigorous for the bilayer,
approximate for coupled hybrid kinetics; tilt [OPEN]. Compiles clean (pdflatex).

CLEANUP (vs upside2-md-master): master-diff inventory via subagent. Findings + fixes:
- C1 MASTER PARITY (done, scratch-build verified; obj/ left intact for the running calibration):
  CMakeLists_M1 OpenMP restored (libomp present); main.cpp inner_step default back to 3 (1 only on the
  MARTINI --duration-steps path); end-of-run frame collect gated to duration-steps; cosmetics reverted
  (%.2f, help string, indent, one-line logger); deriv_engine.h integration_stage -> 10-arg w/ nullptr
  defaults. KEPT (required/used, not churn): h5_support HDF5-2.x rewrite (system HDF5 2.1.1), box virtual
  (used via dynamic_cast), CMakeLists_Other (clean additive). Finding 46: effective-time-factor comment
  was inverted -> fixed to physical_time = sim_time / R in run_sim_hybrid.sh + prep lib.
- C3/C4/C5 (in progress, delegated to a bounded subagent w/ verify-0-callers + scratch build): delete
  verified-dead C++ (martini_hybrid/fix_rigid/masses/stage_params/box + martini.h decls), MartiniNodeRegistrar
  dead guards, dead py functions, py shebang/future. Exclusions: RESPA dual-path, box virtual, live env channel.

CALIBRATION (running, b18x4dw4y): long 303K overdamped reference for R(t) plateau + R@27(T) T-check;
refine R vs the documented 7.3. Re-run on the clean binary after cleanup.

## Remove legacy analytic LJ/Coulomb encoding from MartiniPotential (2026-07-19)

All pair interactions already evaluate ONLY through the pre-fitted combined spline; removed the dead
analytic scaffolding. Files: src/martini_potential.cpp, py/martini_prepare_system_lib.py.
- C++: deleted node-level epsilon/sigma/coulomb_k fields + attr reads (coulomb_constant read at load
  actually THREW on eq_m90, which lacks that attr -> removal also fixes loadability); dropped
  eps/sig/qi/qj/qq from PairParam; replaced the map<(eps,sig,qq)> of splines with a vector<spline>
  indexed per unique coefficient row (param_table[ip].combined_spline = &combined_splines[ip]); added
  coeff_is_zero[] (all-4-raw-coeffs==0) for the active-pair check; unified lj_cutoff/coul_cutoff ->
  single `cutoff` attr (=12.0). Removed dead unique_lj_params + <set> include. Construction-time
  (eps,sig,qq)->grid-row matching KEPT (self-validating table lookup, not analytic eval).
- Python: stopped writing dead epsilon/sigma attrs (+ weighted_median/coeff_counts helpers) and dead
  lj_cutoff/coul_cutoff attrs. KEPT coulomb_constant: the soft-core builder (inject_particles_table,
  has_softening) reads it and IS live in npt_equil/npt_equil_reduced equilibration stages -- task's
  "softening never active" premise was incorrect (STAGE_PARAMS enables it), so removing it would break
  equil prep.
- GATE: fresh scratch build make=0; eq_m90 initial FF potential = -11039.306640625 (== baseline
  -11039.306641). Refactor proven bit-identical to the dead-reads-only intermediate across the full
  2-step trajectory (potential + positions, max abs diff 0.0). obj/ NOT rebuilt.
