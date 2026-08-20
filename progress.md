# Progress log

## Current phase (2026-08-10/11) — findings-88 fix deployed; NP dt corrected

### findings-88 fix — COMPLETE, production running
`HybridPositionNode::propagate_deriv` now routes the O-site sensitivity through the placement Jacobian
and BB's fourth weight is no longer dropped. `_write_nonbonded_pairs` restricts the protein side to BB
only (N/CA/C/O env pairs = 0). Seeds rebuilt for all four glpG variants and all six NP faces on the
corrected pair table and binary. Jobs:
- **glpG REMD**: jobs 53233848/49/51/52 (79HIS, 79HIS_S115T, 79ALA, 79ALA_S115T), PENDING as of 08-11.
- **NP**: job 53251911, RUNNING as of 08-11, dt=0.001, block 2.

**Still open:** avg_kinetic_energy/1.5kT is +2.1% above 1.000 after the fix. Second cause unidentified;
not blocking production. `--potential-deriv-agreement` is unusable as a gate on these systems (round-off
dominated at this potential scale).

### NP dt corrected — COMPLETE
Job 53234804 (dt=0.005) had runs 1 and 4 destroyed at t≈250 during protein unfolding. Root cause:
backbone spring bonds (k=48, ω=9.8 rad/t_u) are driven to 3–7× thermal amplitude by conformational
changes; at dt=0.005 the integration accuracy is insufficient to track force reversals at large amplitude.
Reverted to dt=0.001; /output from destroyed runs deleted; resubmitted as job 53251911. See findings.md
"NP-1AO6 blow-up at dt=0.005 after findings-88 fix".

## Environment notes (still valid)
- midway3 default shell is zsh (1-indexed arrays) — use explicit index mapping in upload loops.
- A micelle has no fixed normal (asphericity fluctuates 0.19–0.33); depth and tilt must be measured
  against the aggregate's instantaneous short principal axis, never box z.
- Cancelling a job is not always a no-op — a 9-second run was enough for reseed() to consume /output
  on all six NP configs and strand the next job.

## 2026-08-13
- **NaN investigation (glpG-DDM micelle REMD).** Located the blow-up origin (slot 45, T=0.890978,
  frame 129) and showed it propagates through the ladder by exchange, which is why 48/48 replicas die
  from one event. Ruled out by measurement: the Update-88 backbone over-count (seeds carry BB only),
  stale pair list, minimum image, the table-build formula (2.3e-12 vs analytic), timestep margin in the
  sampled range (one-step dx 0.0058 A at r=3 A), and thermal access to the core (~500 kT). Found a real
  defect in the table -- `r = max(r, 0.1*sig)` gives a constant 3.16e12 E_up plateau, i.e. exactly zero
  force below 0.47--0.60 A, with a peak tabulated force of 1.10e14 E_up/A. **Trigger still unidentified**;
  fixing the floor needs a domain/spacing decision. See findings 90.
- **Fixed: MBAR returned uniform weights for hybrid coupled potentials.** Referenced `cE0` to its pooled
  mean in `00.AnalysisScripts/helpers/calc_hdx_ht.py` and `4.calc_D_uptake.py`. Exact transformation;
  protein-only path unchanged. Left the two master-identical MBAR sites alone. See findings 91.
- Files modified: `example/00.AnalysisScripts/helpers/calc_hdx_ht.py`,
  `example/00.AnalysisScripts/4.calc_D_uptake.py`, `findings.md`, `plan.md`, `progress.md`.
- Delivered the four-variant HDX dG-vs-residue plots to ~/Downloads (all four MBAR solves converged
  after the fix); cancelled the stuck 79HIS job 53294511.

## 2026-08-15
- **BB proxy reworked onto Upside's own derived sites, and the old placement deleted.**
  `HybridPositionNode` now takes `pos` + `infer_H_O` and builds BB as the mass-weighted N/CA/C/O centre,
  so its derivative is a constant-weight split and no hand-written placement Jacobian exists any more
  (184 lines removed, plus the `sc_env_backbone_hold_steps` flag and the reference-frame datasets).
  Verified on 1rkl (20 000 steps): 0 non-finite, potential -18173 -> -19810, Rg 12.76 -> 12.43, worst
  C-N 1.53 A, KE/1.5kT 1.013, force on BB and O slots exactly 0, total force ratio 1.69e-09.
- **Fixed: `engine_c_library.cpp` never registered the MARTINI nodes**, so every analysis through
  `upside_engine` evaluated a different model than the binary (-12827 vs -18172 E_up).
- **Deployed to midway3** (7 files, rebuilt `upside` + `libupside.so`) and migrated the configs:
  new `py/martini_upgrade_hybrid_args.py` rewrites the one-argument `martini_hybrid_position` attribute.
  Applied to the four queued POPE/POPG seeds; wired into `np_prod.sbatch` for the six NP replicas, which
  a running job holds open. Confirmed the four queued glpG jobs need no resubmission (see remote_jobs.md).
- **Relaunched the local wildtype glpG REMD**, 16 replicas, 400 000 steps (~12.3 h at 32 400 steps/h),
  warm-started from the previous run's final frame. That run had died at ~153 k/200 k steps with its
  potential still falling monotonically, so it was all relaxation -- reusable as an initial condition,
  not as MBAR samples.
- **NP footprint analysis** (`scratchpad/np_footprint/`): per-residue CB-to-NP contact over ~2000 frames
  x 6 orientations drawn from ~1.3 M frames, to test the Carlson et al. 2025 core-histidine claim. CB is
  reconstructed the way `affine_alignment` does and matches the engine to 1.5e-4 A.
- Files modified: `src/martini_hybrid.cpp`, `src/martini_internal.h`, `src/engine_c_library.cpp`,
  `py/martini_upgrade_hybrid_args.py` (new), `py/martini_prepare_system{,_lib}.py`,
  `scratchpad/local_popg_79HIS/{run_local.sh,hdx.sh,warm_start.py}`, `scratchpad/np_footprint/*`,
  `findings.md` (94), `remote_jobs.md`, `progress.md`.
- **HDX pipeline validated before committing the weekend to it.** Dry-ran the whole chain against the
  first 5% of the new local run, into a scratch work dir at reduced parallelism. It completes, and the
  MBAR diagnostics measured directly: f_k spread 106.9 (degenerate case is exactly 0.000), neighbour
  overlap 0.093-0.268 (median 0.233, better than the 48-replica DDM ladder's 0.115-0.128), ESS
  632-1063 of 9120 pooled frames. So the thinned 16-replica ladder reweights fine; the limit is sampled
  opening events, which is what pins 65-112 of 203 amides at the p_f = 1 sentinel and compresses dG to
  -5..+7 against the reference's 10-20 kcal/mol peaks.
  Consequently relaunched the local run at **1 300 000 steps** (~40 h) warm-started from step 22 200,
  frame interval widened to 60 steps (~30 G total). Frames stream continuously, so it can be analysed at
  any point. `hdx.sh` is now parameterised (`HDX_WORK`, `HDX_START_FRAME`, `HDX_JOBS`, `HDX_OUT`) and its
  relative `--trimmed-dir` bug fixed.
- **NP footprint moved to Slurm** (53372772) after the login-node run lost a race with the production
  job's `reseed()` renaming `output` mid-read; `np_footprint.py` now reads only rotated blocks. It also
  keeps per-frame min-distance/burial/Rg instead of running sums, so the compact (experimentally
  comparable) window can be separated from the spread state after the fact.
- **Carlson et al. 2025 claim identified precisely** from the local draft: decreased labelling at K12,
  K73, K190, K525 (water) and K190/K525/K541 (TES), a common site "around K190 (Subunit II)", and the
  protein preferring to "open and expose its center". First (whole-trajectory) pass put K190 at the 19th
  percentile among 58 lysines, below its burial-matched background -- but that average mixes in a nearly
  unravelled chain, which is why the per-frame re-analysis was needed before drawing a conclusion.
- **Staged the cluster HDX path so it is ready when block 1 lands.** New `py/martini_remd_concat.py`
  joins a chained run's rotated chunk groups into one `/output` (required: `/output` alone is the last
  ~300 frames) and drops whole any chunk with a non-finite potential. Verified on a real 3-chunk file:
  ordering, boundary frames, all 16 datasets, and strided output identical to the naive slice.
  `hdx_cluster.sbatch` then runs the same pipeline as the local job. Consolidated `write_hybrid_energy.py`
  and `plot_ref_style.py` out of scratchpad into `example/00.AnalysisScripts/helpers/` so local and
  cluster share one copy, and re-verified the local dry run afterwards.
- **Near-miss caught: the cluster's analysis scripts had drifted.** Its `calc_hdx_ht.py` and
  `4.calc_D_uptake.py` still lacked the findings-91 reference subtraction, so a cluster HDX run would have
  silently returned uniform MBAR weights. Re-uploaded. Only the C++ build is kept current by install.sh.
- **Retracted the undersampling diagnosis; the flat dG profile was an analysis ceiling (findings 96).**
  The T-slice was using `residue_dg_from_pf_step6_plot`, which clips mean_pf at 0.99999 and so caps dG at
  6.82 kcal/mol at T=0.85 -- exactly where the profile sat. On the same frames the jscripts convention
  gives dG max 19.58, and the number of residues with pf == 1.0 exactly is zero, so every
  "never-exchanging" amide was manufactured by the clip. Fixed in `calc_hdx_ht.py` and re-deployed to the
  cluster. The 10%-of-run profile now spans -5..+22 with the six-helix band structure visible and is in
  ~/Downloads. Remaining off-scale residues (32 at T=0.85) are a real resolution limit and keep improving
  with sampling, so the 1.3 M-step run continues.
- **Made the local run shutdown-safe (machine off Sunday 06:00).** Resized it to 740 000 steps so it
  finishes on its own at ~05:20 Sunday rather than being killed mid-write. Restarts no longer cost data:
  `reseed.py` rotates `/output` to the next `output_previous_<n>` and restarts from that segment's last
  healthy frame, and `hdx.sh` now joins the segments with `py/martini_remd_concat.py` before projecting --
  reading `/output` alone would have silently analysed only the newest segment. Verified on the real
  3-segment files: 2448 + 697 - 300 discard = 2845 frames joined, matching exactly.
  `resume.sh [STEPS]` is the one-command recovery after the reboot (or after any blow-up).
- **Local ladder made self-healing for the unattended overnight window.** System 15 (T=0.90) blew up a
  second time at step ~226 800; `supervise.sh` now runs rounds against a wall-clock deadline and reseeds
  destroyed replicas between them. Two of my own data-handling faults fixed alongside (findings 99):
  `reseed.py | head -3` SIGPIPE'd the writer so only 3 of 16 replicas were reseeded at the 06:15 restart
  (no corruption -- every seam is continuous at 0.0000 A -- but 13 replicas re-ran ~41 000 steps), and the
  concat dropped destroyed chunks whole, which would have discarded 2332 good frames from run 15 to remove
  ~100 bad ones. Now filtered per frame on finite-and-negative potential.
- **Plot now reproduces the reference structure.** With the healed multi-segment data: 163-179 of 203
  residues resolve (was 128-171), off-scale down to 24-40, dG -4..+20. Peak/valley positions match the
  reference (valleys 50-70, 118-133, 150-162, 178-190). Remaining difference is spikiness and the
  off-scale count, which is genuine sampling refinement.
- **Local wildtype run stopped at the user's call** with 11 533 frames/replica (~692 k steps, 184 528
  pooled). Final HDX figure in ~/Downloads: 187/203 residues resolved at T=0.85, 16 off-scale, dG -2.4..+18.6.
- **Root-caused the residual flatness (findings 100/101).** Ruled out by measurement: integrator, H-bond
  assignment (agrees with DSSP to 8%), lipid voids (none), hydrophobic mismatch (-2.8 A), burial threshold.
  It is tertiary packing: helical-core CA-RMSD plateaus at 4.15 A (DDM 2.61 A), so both protection terms
  fail together 3.34% of frames -> dG 1.99, reproducing the observed median. Localised to three
  protein-protein terms MARTINI cannot supply (protein self-burial, backbone burial coupling, and the
  sidechain-backbone H-bond competition solved inside the rotamer node). The sidechain 1-body field was
  deliberately and correctly replaced by martini_sc_table_1body; membrane.h5 is correctly absent.
  Proposed remedy documented as plan.md RD1 -- NOT implemented, needs sign-off (changes rotamer arity).
- **RD1 built and running as a 3-arm experiment** (`scratchpad/rd1_env_test/`): `base` (as shipped), `env`
  (protein self-burial restored), `envfull` (all non-membrane Upside terms restored; differs from a
  standard config by one node, `placement_scalar`). Built with Upside's own writers; rotamer arg list
  extended, no C++ change needed since `RotamerSidechain` sums a variable-length prob_nodes vector.
  envfull pre-flight 3000 steps: KE/1.5kT 1.009, Rg 19.44->19.62, 0 broken, 0 non-finite.
  2 M steps each, finishing before the Sunday shutdown; `rd1_resume.sh` recovers from an interruption.
  First reading is uninformative (all three still at the shared start, envfull only 6 k steps).
- **Fixed a systematic 0.568 A CB placement error (findings 102).** `affine_alignment` centres each residue
  frame on the centroid of N/CA/C, so a placement in that frame must be centroid-relative;
  `martini_prepare_system_lib.CB_PLACEMENT` stored the raw CB. Every sidechain interaction site -- the
  anchor of the entire SC-env term -- was displaced by exactly the centroid. Found by comparing the hybrid
  config array-by-array against a standard Upside config of the same protein (the bonded and H-bond core
  came back bit-identical, which is what made the two differing arrays worth reading). Fixed to match
  upside_config to 5.8e-8; spline tables unaffected. Identical-configuration potential moves ~+680 E_up.
  RD1 arms rebuilt on the corrected placement, stability-tested, relaunched at 1.2 M steps.
- **RD1 rejected (findings 103).** On the CB-corrected placement, restoring the protein-protein
  environment terms changed helical-core RMSD by +0.10 A (self-burial) and -0.08 A (all terms) against a
  4.61 A control -- inside scatter, against a ~1.5 A target. The CB fix itself also did not improve fold
  fidelity, though it was a real bug worth fixing. Cause of the ~4.5 A deviation remains unidentified; the
  one untested node-level difference is the rotamer 1-body representation (fixed vs rama-dependent).
  Stopped the arms and relaunched the wildtype ladder on the CB-corrected seed alone (1.4 M steps,
  `scratchpad/local_popg_cbfix/`), which is the Tuesday deliverable path.
- **Verified the CB correction against ground truth, then applied it to the cluster.** Placing CB by each
  convention and comparing to the actual crystal CB atoms (187 residues): raw/as-shipped deviates 0.576 A
  mean, corrected deviates 0.047 A -- the residual being idealized-geometry error. Independent of any code
  reading. Then patched the four cluster seeds, cancelled 53349015/17/18/20 at 16 chunks of block 1, moved
  that data to `<variant>_pre_cbfix/` (102 G, reversible), and resubmitted as 53410263-66 from the
  corrected seeds. Nothing analysable was lost -- protection-state variance is zero that early in a chain.
- Local CB-corrected ladder resumed after the shutdown (146 k steps preserved, all 16 replicas readable);
  running 200 k more, sized to finish ~14:30 before the next expected shutdown.
- **Retracted findings 96 (see findings 104).** Master deliberately uses the clipped step-6 convention for
  the dG T-slice; I had removed it. The clip's 6.8 kcal/mol ceiling nearly coincides with the estimator's
  measured statistical limit (5.5 kcal/mol at ESS 10 697 / 85 760 frames), and removing it let dG reach
  19.6 where 100% of values above 8 were carried by <2 effective frames -- the source of the
  spike-within-a-helix pattern the user objected to. Reverted; the profile is now smooth up to the ceiling
  with unresolvable residues marked off-scale. Reference peak magnitudes of 10-20 are unreachable with a
  binary protection state (would need ESS ~1e14), so that figure used a different estimator.

### 2026-08-16 (late): HDX figure fixes 2 and 3 -- DSSP helices and helix N-terminal caps
Established by measuring the **unweighted** protection state that the excursions and the roughness are three
different things (findings 106), then that the "broken helices" were an annotation error (findings 107).
Applied the two annotation fixes:
* `plot_ref_style.py` rewritten with `--structure`: DSSP helices shaded, amides in the first four positions of
  a helix drawn as open markers (no $i-4$ acceptor, so fast exchange is expected there), and a per-region
  summary printed (helix interior / cap / loop medians, plus any helix-interior amide below 0).
* `scratchpad/local_popg_79HIS/hdx.sh` now reports one temperature (`HDX_REPORT_T`, default 0.85) and passes
  `--structure`, since the cold rungs contribute mostly estimator collapse.

Deliverable regenerated: `~/Downloads/glpG-RKRK-79HIS_POPEPOPG_CBfix_dG_vs_residue.png`.
T = 0.85, 180 of 203 amides resolved. Helix interior median 2.43, N-terminal cap 1.23, loop 0.07 kcal/mol.
**Only 3 helix-interior amides fall below 0: idx 98, 99, 150** -- the same three that were H-bonded in the
prepared structure and lost it in the run. That is the whole of the remaining fold defect, and it is Fix 4.

### 2026-08-17: the within-helix spikes are an analysis omission, not the force field (findings 113)
User pointed at `glpG_dG_T_slice_masterstyle.png` (= `local_popg_cbfix/hdx/results/..._DG_res_T_slice.png`):
inside a TM helix every amide should read `+inf` as a block, as the implicit-bilayer reference does, and
instead it reads 0.5--4 with isolated needles. Diagnosed on **unweighted** single-replica protection states,
so no MBAR or clipping is involved: the hybrid pipeline scores a lipid-facing TM amide as exchanging whenever
its backbone H-bond flickers, because `get_protection_state.py` only knows about *protein* burial and the
hybrid HDX topology has no membrane `surface` node for `--use-TM-region`. The replacement,
`py/martini_hdx_membrane_accessibility.py` -> `combine_hdx_protection.py --water-accessibility`, exists in the
repo and `hdx.sh` never calls it, so `PS.npy` == `PS_protein.npy` bit-for-bit.
Applying it: bilayer-embedded helical amides at `+inf` go 35/92 -> 79/92, contiguous runs 8,6,5,4 -> 17,16,15,14.
Figure: `~/Downloads/glpG_PS_lipidshield_diagnosis.png` (two panels, same trajectory).
Retracts findings 106 item 4 (measured the omission at 6 amides, ~7x too small, via a >50%-exchanging cutoff).
Residual is 13 of 92 deep amides at 2.3--4.2 kcal/mol -- the real, much smaller fold defect.
**Not yet wired in:** `--cutoff`/`--min-contacts` are uncalibrated and set where the bilayer edge falls;
they need fixing against the PO4 leaflet planes before any dG from this path is quoted.

**Fix applied and re-analysed (D8).** `martini_hdx_membrane_accessibility.py` rewritten with a
parameter-free criterion — radius = flat first minimum of the intermolecular tail-tail g(r) (7.000 A,
identical on all 16 replicas; a bare argmin flips between 6.88 and 7.12), threshold = 1 contact, which is
where a phosphate bead sits (median 0 tail neighbours there, vs ester 2, first tail 6, terminal tail 11).
`--cutoff`/`--min-contacts` removed. A slab between the leaflet planes was measured and rejected: it also
protects the amides lining glpG's polar interior, giving +inf runs of 51/46/42 residues and swallowing
34 of 55 loop donors; master ANDs its own slab with a lipid-facing `surface` test for the same reason.
`hdx.sh` now runs the accessibility step per replica and folds it in with `--water-accessibility`.
Full 16-replica re-analysis in `scratchpad/local_popg_cbfix/hdx_acc/` (11 068 frames/replica, 177 k pooled).
T-slice sentinels at T = 0.85 go 23 -> 94, contiguous runs 3,3,3,3,2,2 -> 16,15,14,12,10,8, and per helix
28-47 17/19, 160-174 13/15, 104-126 16/23, 184-206 16/23, 134-150 10/17, 81-102 11/22. The two interfacial
helices (49-56, 62-66) stay fully resolved, which is correct — they sit outside the boundary. Negative dG
excursions are gone. Figures: `~/Downloads/glpG-RKRK-79HIS_POPEPOPG_CBfix_dG_vs_residue_acc.png` and
`hdx_acc/results/..._DG_res_T_slice.png`. Caveat for any before/after comparison against `hdx/`: that run
had 5 399 frames/replica, so the two differ by sampling as well as by the membrane term — the clean
comparison is `~/Downloads/glpG_PS_lipidshield_diagnosis.png`, both panels from identical frames.

**Renderer switched to the unclipped convention.** `plot_ref_style.py` collapsed to a single mode: it draws
`dG_unclipped`, so an amide with p_f == 1 runs off the axis and a non-exchanging TM helix reads as one
continuous excursion instead of a capped plateau (the reference layout, `~/Downloads/79ALA_0.85.png`); the
`--unclipped` toggle and the censor-at-the-bound path are gone, and the axis is the reference's -20..30. The
DSSP shading and helix N-terminal open markers (findings 107) are kept, and the per-temperature resolution
limit stays as a dotted line, since values above it are lower bounds. This is only correct now that the
membrane term is in the protection state: the clip was masking the missing term, and findings 104 reverted an
earlier attempt to remove it for exactly that reason. At T = 0.85, 85 of 203 amides off scale (81 in helices,
runs 15/15/14/12/9/8), 118 resolved spanning -0.53 to 21.3, helix-interior median 4.44 against loop 1.63, and
residue 150 the only helix-interior amide below zero at any temperature but 0.90.
Deliverables: `~/Downloads/glpG-RKRK-79HIS_POPEPOPG_CBfix_dG_vs_residue_acc.png` (5 rungs) and
`~/Downloads/glpG_dG_vs_residue_3T.png` (0.75/0.80/0.85, the reference's own temperature set).
`calc_hdx_ht.py`'s `_DG_res_T_slice.png` deliberately keeps master's clipped step-6 convention (master parity).

**Figure styling matched to the reference (`~/Downloads/79HIS_0.85.png`).** The annotation added in the
findings-107 pass — DSSP helix shading, open markers on helix N-terminal amides, a dotted resolution-limit
line per temperature — produced a legend with eight entries in two columns that covered the top-left of the
plot. All of it is off the figure now: bare white axis, a light grey grid, one line and one marker per
temperature at the reference's weights (blue circle / red square / green triangle, lw 2, ms 8), the
reference's -20..30 axis, and a legend of temperatures only. It sits lower-left rather than the reference's
upper-left, because nothing in this profile falls below -1 kcal/mol while the reference's own corner is under
the excursions off helix 1 here. The helix bounds and the helix-interior / N-cap / loop medians are printed by
the script instead of drawn, so nothing was lost. `--temperatures` now defaults to the reference's
0.75,0.80,0.85. Deliverables re-rendered: `glpG_dG_vs_residue_3T.png` (reference's three rungs) and
`glpG-RKRK-79HIS_POPEPOPG_CBfix_dG_vs_residue_acc.png` (all five).

### 2026-08-17 (later): cluster HDX is empty because the cluster protein is frozen (findings 116)
Ran the findings-113 pipeline on the 48-replica cluster ladder (synced 3 analysis files, wired the membrane
term into `hdx_cluster.sbatch`, added `HDX_WORK`; 48 replicas process in 7 min). All four variants came out
degenerate — 188 of 203 amides off scale, resolved values to −53.9 kcal/mol. Root-caused through three wrong
answers of my own: MBAR is healthy (ESS 13.4%, overlap 0.25, better than local); the −2.5σ energy drift is real
but not the cause (the equilibrated tail alone gives the same result); the membrane term is not the cause
(protein-only p_f is already 1.0 for 169 of 203). The protein has **zero internal dynamics** — raw internal CA
RMSD 0.000 Å across chunks — because `/input/stage_parameters.current_stage` is `production_handoff`, and
`martini_hybrid.cpp:637-641` holds the protein rigid for any stage that is not exactly `production`. The local
seed says `production`. One attribute; fix is `set_stage_label(..., "production")` + restart, not applied.
Also confirmed the local run has **no** PBC problem: 0 CA–CA breaks >5 Å in the projected trajectory over
11 068 frames, all protein coordinates inside the box (max |x| 44.9 vs half-box 49.7). The PBC-crossing tail in
the VTF is a rendering artifact of `martini_extract_vtf.py`, which centres on the protein's circular-mean COM
and then wraps everything into [−L/2, L/2); with a 69.5 Å protein extent in a 99.3 Å box that centre is
ill-conditioned and can image part of the chain across. The analysis never sees it — Rg/RMSD/protection run on
unwrapped projected coordinates with no periodic box.

### 2026-08-17 (evening): cluster restart verified in flight, local run extended, poster moved to OneDrive
* **Cluster pull/rebuild is a no-op for the physics.** All 102 files in `src/` are byte-identical before and
  after (the apparent mismatch was macOS vs Linux `sort` collation); the pulled commit `a31d335` touches only
  analysis Python and docs. Nothing needs re-running on that account. The four glpG jobs were still PENDING at
  rebuild time so they start on the new binary; `np_1AO6_prod` (53411347) is stepping cleanly with 0 errors.
* **Local ladder extended for the poster figure.** Reseeded all 16 replicas (rotated to
  `output_previous_5`, restart from frame 5268) and relaunched at 700 k steps, ~21 h, into
  `scratchpad/local_popg_cbfix/`. First step healthy: total_potential -21574, Rg 20.9. Roughly doubles the
  11 068 frames/replica the current figure rests on. `hdx.sh` can be run against partial output at any time.
* **Poster updated and moved out of the repo** to `/Users/yinhan/OneDrive/hybrid_interface_poster_2026-08-18/`
  (OneDrive running, 7.1 M, rebuild verified from the new location; every hardcoded path repointed, including
  `render_system.pml` and both renderers' output paths). `make_figures.sh` now defaults to `hdx_acc` and takes
  `POSTER_HDX`; the `--unclipped`/censored figure calls are gone with the flag. Content changes: added "The
  protection state needs the lipid" (findings 113), updated frames to 177 088 and the T = 0.85 medians to
  4.44/4.26/1.63, retracted the ESS-censoring-as-novelty claim (findings 114), and **dropped the −3.4 kcal/mol
  per T_up temperature slope** — on the 83 residues resolved at every temperature the median slope is −0.3 with
  52% negative, so it is not systematic; the never-exchange count falling 119 → 85 of 203 is what survives.
  Added `README.md` there with the rebuild recipe, the overflow constraint (all three columns full), the pending
  refresh after the longer run, and the two claims that are easy to get wrong.

### 2026-08-17 (late): re-checked the poster's failure claims, rewrote that section
Re-measured every limitation on the current dataset rather than reasoning from the analysis change. **All hold**,
two are worse than the poster stated: helix-interior H-bond occupancy 0.636 (was 0.681), 34 of 108 donors below
0.5 occupancy (was 24), protein/lipid temperature excess +5.74%/+1.08%, occupancy-vs-friction r = +0.537.
The spike fix was an analysis fix; these are simulation properties, so none moved.

**It did change how they must be read.** With the membrane term, helix-interior protection is 0.989 against an
H-bond occupancy of 0.636, and 89.3% of broken-H-bond frames are still scored protected (was 79%). So the clean
profile is *weaker* evidence about the H-bond network than before, and that caveat now belongs on the poster.

Also fixed: `fig3` was unreproducible — `make_fidelity_figures.py` read a raw array from a **previous session's
temp directory** that no longer exists, so `make_figures.sh` would have failed. The array now lives beside the
figure as `raw_hbond_occupancy.npy` and `HDX` points at `hdx_acc`; fig3 and fig4 both regenerate.

Rewrote the limitations section at the user's request — it read as machine-written. Dropped the "two limitations
each with a defined next step" framing, the 5-row dt-scan table (now 3 rows) and the self-auditing prose; kept one
plain point ("The helices are too loose ... rank those numbers against each other, not against experiment") plus
the timestep/coupling split. Filled the freed space with fig4, which is real data supporting it. Fixed two real
errors found while editing: prose quoted 0.64 occupancy directly above a table reading 0.711 (long ladder vs
matched 30 k controls — different runs), and the box said "the profile above" for a figure in another column.

**Added the implicit-membrane comparison to the poster.** `79HIS_0.85.png` now sits beside this work's profile
as one panel (`make_comparison_figure.py` -> `fig1_dG_comparison.png`), replacing the hybrid-only figure. The
two are genuinely comparable in observable, construct, temperatures (0.75/0.80/0.85) and axis (-20..30), which is
why the panel works; they are **not** comparable in pipeline (the implicit figure came from legacy jscripts with
`--use-TM-region`) and its per-residue numbers are not available, only the PNG. RESULTS.md records that, and that
the three comparison statements on the poster are read off the figures rather than computed from shared arrays.
Each source panel is cropped to its ink bounding box before compositing -- two panels in one column leaves ~5 in
each, and the source margins were the difference between legible tick labels and not.

### 2026-08-17 (night): poster restructured -- Validation and stock-Upside sections out, NP methods panel in
* **Removed `Validation` and `Against stock Upside`** at the user's call: the implicit-vs-hybrid comparison now
  carries that role. The stock table's load-bearing numbers moved into the limitation box (occupancy 0.71->0.76
  and cooperativity 1.80x->1.53x under dt/4, against stock's 0.95 and 4.18x) so the "survives dt/4" claim keeps
  its evidence. Note the insertion / bilayer-integrity / no-numerical-failures statement is now unstated on the
  poster; it is validation of the system rather than of the result, and the comparison does not cover it.
* **NP panel added as a METHODS claim, not a result** (user's choice after I laid out the options). Two rendered
  states side by side -- albumin folded beside the 5 nm MPA-AuNP, and albumin spread over it late in the run,
  DSSP helix 74% -> 52% -- plus three sentences: the coupling needed no change for a different environment, and
  the same excess opening shows up here as albumin failing to adsorb folded.
  New, all reproducible: `export_np_snapshot.py` (cluster-side frame -> PDB, chain per particle class, chain-walk
  PBC unwrap since a spread albumin exceeds half the box), `render_np.pml`, `make_np_figure.py`.
  Two things worth remembering: PyMOL's `dss` finds **no** secondary structure on this CG backbone while
  mdtraj's DSSP finds 74% helix on the same coordinates, so SS is assigned explicitly from mdtraj; and `zoom`
  counts hidden atoms, so the 1194 ions had to be removed rather than hidden or the system rendered as a speck.
* `fig4_interface_fidelity` dropped from the poster -- it was interim filler for the gap the NP panel now fills,
  and its numbers are stated in the limitation box.

**NP re-run assessment (findings 95 re-tested on 74% more data, and a defect found).** The current footprint
(12 312 frames sampled from 2.32 M, against 1.33 M before) reproduces every finding-95 conclusion: K190 contact
**0.000, the lowest of all 58 lysines**; K525 0.542 (81st pct) and K541 0.678 (88th) supported; dominant
orientation still 71%. The adsorbed-and-compact window **shrank, 3.2% -> 1.9%**, i.e. more sampling makes it
worse. Separately, the NP configs still carry the **uncorrected CB placement** of findings 102 -- their
`placement_fixed_point_vector_only_CB` reads [0.0000, 0.9438, 1.2068] against the corrected
[-0.0198, 1.5117, 1.2068], the 0.568 A centroid offset -- which matters doubly because the footprint is
CB-anchored. `current_stage` is `production`, so no rigid-protein bug. A re-run therefore needs the CB fix
**and** a larger box (the box is near-vacuum, 5174 beads in 200^3 A, so enlarging it is nearly free), but neither
addresses the over-unfolding, which is most likely the same coupling defect as glpG's loose helices. Not
resubmitted -- doing so now would spend days to reproduce the same dead end.

**Poster, final structure.** Dropped the H-bond occupancy panel (`fig3`) at the user's call; its headline number
(occupancy 0.64 against the starting structure's 0.95) moved into the limitation box, which is now the only place
that claim is made, and the dt-scan numbers there are labelled as matched short controls so they cannot be read
against the ladder number. Restored one sentence of system validation (insertion, 38.6 A bilayer, no numerical
failures) since the implicit-vs-hybrid comparison validates the result and not the system. Column 3 is now:
interface checks -> the nanoparticle system (two rendered states + the footprint against Carlson et al.'s
lysines, `make_np_footprint_figure.py`, current data) -> the limitation box. The NP panel states plainly that it
is thesis work in which the simulation never reproduced the labelling experiment, and that the reason is the same
excess opening. RESULTS.md gains a Result 4 section with the re-measured numbers and the two defects to fix
before any NP re-run.

**Poster prose rewritten for register (user: "make it formal, professional, more human, not AI like").** Went
through all three columns, not just the limitation box. Removed the constructions that read as machine-written:
the "reaches the first ... gives up the second" schematic pairing, the triadic roadmap ("the profile, its
validation, and what remains to improve"), self-referential hedging ("which is why it is worth saying here"),
the grandiose "the ordering the physics requires", and the heavy em-dash chains throughout. Replaced with plain
declarative sentences of varied length in normal scientific register; "The main thing still wrong" became
"Principal limitation". Every number is unchanged.

Two real errors surfaced while rewriting: "177 088 pooled frames" was stated in both column 2 and column 3, and
the box's first paragraph asserted the occupancy claim twice over. Both fixed. Column 1 keeps a little slack
(0.029); columns 2 and 3 are at 0.014 and 0.013. `.pptx`, `.pdf` and preview re-rendered; `poster_content.py.bak3`
holds the pre-rewrite text.

**Real temperature units on the poster.** Upside's reduced T meant nothing to an outside audience. 1 T_up =
350.588 K, so the ladder is 245-316 K and T = 0.85 is exactly 298 K (which is why 0.85 is the reference rung).
The poster now gives the conversion once where the ladder is introduced, then uses Kelvin in the body text:
298 K for the reported profile, 245/316 K for the resolution bounds and the cooling trend, 263/280/298 K in the
figure caption. **Figure legends stay in reduced units on purpose** -- the implicit-membrane panel is a supplied
PNG that cannot be re-rendered, so relabelling only our half would make the comparison inconsistent; the caption
carries the mapping. RESULTS.md gains a temperature-units section with the full table.

**Poster cut to poster length; full text kept as a reading copy.** 1222 -> 795 words, a 35% cut, and the freed
space went into type rather than white space: body 24 -> 29 pt, headings 34 -> 38 pt, so it reads at three feet.
Paragraphs became fragments where they could, the three-conventions passage collapsed to one sentence, and the
"Interface checks" prose became bullets. No number or caveat was dropped -- only words.
The full version is preserved verbatim as `poster_content_reading.py` and rendered to `*_reading.pptx/.pdf`;
both renderers now take `POSTER_CONTENT` and suffix their output, so the two cannot overwrite each other.
Two defects the larger type exposed and I fixed: the heading "The same interface, a different environment" was
**clipped, not wrapped** (42 characters does not fit a 11.9 in column at 38 pt -- now "Same interface, new
environment"), and a long bullet ran into the gutter. Both are silent failures -- the OVERFLOW check only watches
column height -- so the preview has to be looked at. Noted in the README.

**NP panel reframed; NP campaign audited and found two generations behind (findings 117).** The poster no longer
claims the simulation failed to reproduce the labelling experiment -- that would be reporting a defective run as
a result. It now shows both states (prepared and adsorbed) with labels carrying no measurement, captioned
"Sampling in progress", and the footprint figure is off the poster. The freed space went to restoring the
stock-Upside benchmark table, which is what the limitation box argues from.
The audit: the NP configs carry `r_min_ang = 0.00` (pre-findings-92 LJ core table, force-free at short range) and
the uncorrected CB placement (findings 102). The 2026-08-15 `martini_upgrade_hybrid_args.py` migration fixed the
node arity only, and its success had been read as evidence the configs were current. Both defects bear directly
on adsorption, so findings 95 and the 2026-08-17 re-measurement are **superseded as tests of Carlson et al.**
Re-run therefore recommended, needing a rebuild (the LJ tables must be regenerated), a larger box, and shorter
replicate runs rather than six long ones -- not yet submitted, awaiting the user's call on cluster time.

**Jumper Figure 1 added to the poster (column 1, "The Upside core").** Figure 1 of PLOS Comput Biol 14, e1006342
-- the six-step Upside cycle -- downloaded from PLOS after checking the DOI resolves to that title and that the
licence is **CC BY 4.0**; the caption credits the authors and the licence, which is what CC BY requires. It earns
its space by replacing text: it shows O/H placement, side-chain decoration, the joint rotamer solve and the
pull-back of forces onto C/N/CA, so three of the five "Upside core" bullets went. It also shows exactly where the
environment couples (steps 4-5), which is the point.
Paid for the rest by moving References to the foot of column 3 and dropping reference 7 (Gronbech-Jensen & Farago,
the g-JF propagator) -- **it had become uncited** when the thermostat text was cut; found by extracting every
[n] from the rendered blocks and diffing against the list. DSSP renumbered 8 -> 7.

**NP campaign rebuilt and relaunched (53456240).** Killed 53411347, rebuilt all six orientations locally with
`build_all.py`, verified each carries `r_min_ang 0.30` + 2-arg `martini_hybrid_position` + corrected CB + stage
`production`, and proved one loads *and integrates* on the cluster binary (5000 steps, Rg steady 25.8-26.0)
**before** deleting anything. Then removed the 206 G of pre-fix trajectories; `prod/` is 970 M, old `.out` logs
kept, snapshots archived off-cluster. Relaunched with `NP_MAX_BLOCKS=3` rather than the default 8 to cap disk,
since the informative window is early -- extendable by resubmitting.
Two things worth remembering. (1) `np_hybrid.py` was a **third** stale generation: it wrote
`martini_hybrid_position` with one argument, so the rebuild failed to load until fixed at source; the running
configs had only ever been patched post-hoc. (2) `--frame-interval` is in **time units, not steps**: my first
cluster load test passed `500` with dt = 0.001, which is 500 000 steps, so only step 0 was recorded and the
kinetic-energy average over one frame printed `-nan`. That looked like a defect and was not one -- with
`--frame-interval 0.5` the run reports `avg_kinetic_energy/1.5kT 0.866`, below 1.0 only because the average
includes the cold start from a minimised configuration over one thermostat timescale.

**Rigid-protein fix verified by measurement (not by reading the gate).** Internal CA RMSD on the relaunched
53441123 grows 0.730 -> 1.124 -> 1.440 A over the first 100 frames, against 0.000-0.001 A across whole chunks in
the frozen run. The glpG campaign is producing real dynamics again.

### 2026-08-18: extended local dataset, poster finalised, NP site conclusions retracted
* **Local ladder finished** (22 836 frames/replica, 363 152 pooled, KE/1.5kT 1.015-1.023). Re-analysed into
  `hdx_acc2`. `~/Downloads/glpG_dG_vs_residue_3T.png` regenerated, and the poster's ΔG panels rebuilt from the
  new npz so figure and text agree. At 298 K: 123 of 203 resolved (was 118), −0.5 to 18.3 kcal/mol, 80 never
  exchanging (was 85), helix interior 4.49 / cap 3.51 / loop 1.71, resolution bound 5.02–6.46 kcal/mol.
  Helix-interior H-bond occupancy re-measured at **0.588**, not 0.64 — it keeps drifting down with sampling
  (0.681 → 0.636 → 0.588), so the helices loosen progressively rather than settling.
* **Poster:** removed "Against stock Upside" and "Interface checks" (user: both read as AI-written), grew the
  glpG system image to full column width, and removed a redundant five-temperature ΔG panel I had added — the
  user correctly spotted it was the same data as the comparison figure's right half.
  Added `fig7_melting.png` from the workflow's melting curve, rebuilt rather than reused because the shipped
  `_Tm_curve.png` has an **inverted** hydrogen-bond axis. It earns its place: MBAR reweighted matches the direct
  average to 1.7 bonds of 162 and 0.05 Å in Rg, which validates the reweighting behind every ΔG on the poster.
  Surveyed the rest of the workflow's outputs in RESULTS.md with reasons for not using them — notably the
  denaturant/m-value panels, which are properties of a synthetic perturbation, and the convergence KS test,
  whose p = 0.0 reflects autocorrelation rather than drift.
* **Slurm:** 53441123 running (19 h, 9 chunks, dynamics verified); 53441124/25/26 still queue-starved after 19 h,
  so there is no four-variant result today.
* **CoF panel (user request):** ran the shipped D-uptake/cooperativity workflow on the hybrid local result
  (`scratchpad/local_popg_cbfix/hdx_acc2`). It first returned an identical COF for all 63 peptides; cause was
  `4.calc_D_uptake.py`'s `legacy_T_range` default of 1.14 (~400 K), which puts k_chem at 4.9e5 s^-1 so every amide
  is exchanged before the first experimental time point (findings 120). With `legacy_T_range=0.85` the curves are
  63 distinct sigmoids and the uptake comparison reaches **mean R² = 0.693 over 57 peptides** — the first real
  experiment-facing number from the hybrid run.
  Fixed two latent shared-script incompatibilities to get there: `5.analyze_D_uptake.py` sliced the 0-d `T.npy`
  the shipped extractor writes, and `write_hybrid_energy.py` emitted a flat `Energy.npy`. Added a peptide-level
  `<pdb>_COF_peptides.csv` export beside the existing residue-level one.
  Poster: `fig9_cof_map.png` (HX-MS above hybrid, 63 peptides, ranked within each panel) replaces the H-bond
  landscape panel. The CoF agreement is weak and the panel says so — Spearman 0.23 over the 55 shared peptides,
  and the workflow's linear R² is not usable because COF normalisation inflates near-flat experimental curves.
  Also centred all poster images and enlarged the Jumper cycle figure (0.186 → 0.215 of poster width).
* **Poster restructured (2026-08-19)** around the model comparison the work is actually for. Added an Abstract box;
  a GlpG panel (`orient_glpg.py` rotates the built PDB onto the mean of its six TM helix axes, 9.1 deg rms tilt,
  and drops the tag; two PyMOL views); an HDX section carrying the stretched-exponential fit and both free-energy
  equations; and **fig10, per-amide dG for Upside's implicit membrane against the hybrid — r = 0.60, Spearman
  0.69 over the 134 amides both resolve.** Dropped the CoF panel (the group no longer uses cooperativity factors)
  and the "Why this is hard" section.
  Three data facts behind fig10, all recorded in RESULTS.md: both axes use the plain hydrogen-bonded-or-buried
  criterion so the comparison is of models not protection states; `--use-TM-region` cannot run on the implicit side
  because its legacy HDX topology has no `surface` node; and 31 of 48 implicit PS files are truncated at 4 MiB
  boundaries, so each is read to its last complete frame.
* **2026-08-19, cluster analysis:** ran `hdx_cluster.sbatch` with `HDX_LIVE=1` on the two wildtype POPE/POPG runs to
  get the plain (hydrogen-bonded-or-buried) protection state for the implicit-vs-hybrid comparison.
  **79HIS (48 replicas, 6081 frames each, ESS 37 492): r = 0.642, Spearman 0.711 over 144 amides** — up from
  0.600/0.690/134 with the local 16-replica set, and only 1 amide saturated against 27. **79ALA, the construct that
  matches the implicit dataset: r = 0.598 over 136** — statistically the same, so the single-residue difference was
  never the limitation; the implicit side is, with 58 of 203 amides censored. Poster panel repointed to the 79HIS
  cluster data, with the protected fractions cached in the poster folder so the figure needs no trajectory.
* **NP block 2 re-measured** (18 246 frames, 3.2% adsorbed-and-compact): reproduces findings 118 unchanged. None of
  Carlson et al.'s five lysines is contacted (K12 0.000, K73 0.005, K190 0.000, K525 0.000, K541 0.033); the patch is
  Lys313 0.33, Glu311 0.31, Asp314 0.27, Asp562 0.26, Lys560 0.26. Max contact 0.328, so still no preferred pose.
* **2026-08-19, poster headline panel rebuilt on matched data.** Found implicit-membrane REMD for the *same*
  construct in `~/cds3/glpG-run3-REMD/79HIS/1` (48 replicas x 11 500 frames) and reduced it with a new
  `popepopg_REMD/hdx_implicit.sbatch` to the plain protection state, ESS 167 977. Paired against the 48-replica
  hybrid: **Spearman 0.662, Pearson 0.546 over 155 amides**, medians 1.36 vs 0.92, everything else held fixed
  (construct, ladder, MBAR target, criterion). Deliberately *not* the highest r on offer — the mismatched-construct
  pairing gives 0.642 over 144 — because that one rests on implicit files 31 of whose 48 were truncated. The two
  implicit datasets agree with each other at r = 0.944, so provenance was the only thing at stake. Pearson is
  cutoff-dependent here (0.55 -> 0.66 as the ceiling region is excluded) while Spearman is not, so the panel leads
  with Spearman.
  Operational lesson: `~/cds3` is `/cds3/...`, which compute nodes cannot read — stage to `/project` first.
* **Integrator reviewed, not changed** (findings 123, per instruction). The hybrid bypasses Upside's three-stage
  Predescu integrator entirely: `integration_cycle` returns after one g-JF Langevin stage whenever `/input/brownian`
  exists, and that list includes the 630 protein backbone sites. `--integrator mv` cannot help — its level 1 is the
  *slow* set at `dt x inner_step`. Timestep left at 0.009.
* **2026-08-19, pooled the four implicit repeats.** `cds3/glpG-run3-REMD/79HIS` has four independent runs of the same
  48-replica ladder, not one; all four reduced (jobs 53539424, 53540593/94/95), each reweighted to 0.85 and the
  protected fractions averaged, pooled ESS 679 839. **Spearman 0.676, Pearson 0.554 over 170 amides**, censored on the
  implicit side down from 47 to 32. Both statistics improve at every fixed cutoff (below 5 kcal/mol: Spearman
  0.662 -> 0.685, Pearson 0.584 -> 0.615), so this is sampling rather than dataset selection. Poster updated.
* **2026-08-19, PI feedback merged (ver3 -> ver4).** Pulled the PI's hand-edited pptx content back into
  `poster_content.py` and regenerated. Removed the p_f equation, credited Peng and Faruk with Crossref-verified
  citations (noting [9] is a Biophys Soc meeting abstract), added a "What was hard" section, and rebuilt the GlpG
  figure from `system.pdb` so it carries no rotation at all — the earlier tilt was my own helix-axis averaging — with
  the leaflets' measured phosphate planes drawn as the two boundary lines. Also repaired three ver3 accidents: the
  abstract body text had been lost, the GlpG heading overlapped its figure, and the References heading was parked
  outside the column.
* **2026-08-19, two figure changes.** Reweighted both models to 0.75/0.80/0.85 on the cluster (`pf_multi_T.py`, one
  MBAR solve per dataset). The profile panel is now an **overlay** of hybrid on implicit in one axis at three
  temperatures (`make_profile_overlay_figure.py`), replacing the side-by-side pair whose left half was a data-less
  supplied PNG; the implicit model saturates across the TM segments while the hybrid follows it in the loops and
  reads lower in the helices. The scatter now shows all three temperatures: **rho = 0.709 (263 K, n=156), 0.739
  (280 K, n=162), 0.676 (298 K, n=170)** — the best agreement is not at the best-sampled temperature, which is
  another argument against the sampling explanation. "What was hard" gained the lipid-to-one-vector-particle dead
  end (PMF route works for side chains; tail bending and compression have no single-particle representation).
* **2026-08-19, ver4 -> ver5.** Extracted the hand-edited deck's structure back into `poster_content.py` (References
  to column two, What-was-hard and Future-directions to column three, bullets, 14 pt refs, author block beside the
  title, bottom-up-PMF item). Restored the abstract body and the "Measuring stability by exchange" body, both of
  which ver4 had lost under live headings. Both comparison figures are now generated at their printed size so the
  script's point sizes are the printed ones; `fig12` uses $\Delta G$; `fig10` split into three panels, one per
  temperature, each labelled with rho, R^2 and n.
