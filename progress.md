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
