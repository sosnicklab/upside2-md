# CURRENT PHASE (2026-09-08): finish ff3.0 on two hosts, then an FF-convergence arm test

ConDiv gly-sym retraining is running on **two** hosts at once and both reach `MAX_STEPS=500` on
Wednesday 2026-09-09 afternoon (job ids and live state in `remote_jobs.md`). They diverged at step
338, so they are independent stochastic trajectories of the same training and will yield two
different force fields. **Everything from here to glpG production must run unattended** — the user
is away from the Mac on Thursday 2026-09-10 and a Claude session only exists while that Mac is on.
The autonomy lives in Slurm scripts on midway2, not in a monitoring session.

## Revised decision 2026-09-08: the arm test comes back, as an FF-convergence test

The 12 h A/B arm test was skipped on 2026-09-07 because the local four-arm experiment settled the
coverage question (only trained-pair + coverage clears TM4). That left `planned_job.md` with an open
item: skipping it removed the **only pre-production check on the final force field**, so production
would have launched on an FF whose TM4 was never measured.

Two trainers finishing at once turns that liability into the check. The arm test is retargeted:

| | old arms | new arms |
|---|---|---|
| arm A / M | plain hybrid + trained pair only | coverage nodes + **midway2's** trained FF |
| arm B / R | coverage nodes + full trained FF | coverage nodes + **rockfish's** trained FF |

Both arms now carry the coverage recipe, which is settled on physics grounds (see "The intra-protein
split" below: coverage is the arm that implements the intended model, not merely the self-consistent
one). The only variable left is which training trajectory produced the tables. That makes the test
an n=2 convergence check on the deliverable:

* the two arms agree within scatter -> training converged, TM4 number is real, launch on midway2's
* the two arms disagree -> 500 steps is not convergence, and that is a finding, not a coin flip
* both fail the health gate -> launch nothing, exactly as before

**And it is not a formality.** Measured at step ~355 (`findings.md`): the two runs drifted the same
distance from ff_2.1 to within 0.1-1.7% but 23 degrees apart in direction, so they sit 33-40% of that
drift from each other and 15-17% rms apart on every trained table. A single ConDiv run does not pin
the force field to better than a third of what training changed. The corollary is a trap the `STEP`
stamp now blocks: the trajectory difference is the same size as ~85 steps of training, so pairing
force fields from different steps measures neither.

**midway2's FF is the declared primary**, so it wins any tie. Reason: the seeds, the install chain,
the production data and the baseline all live on midway2, and it is 2 steps ahead. This replaces the
old tie-break ("Arm B, the self-consistent choice"), which no longer discriminates now that both
arms carry the same node recipe.

**Rockfish's FF is an optional input.** midway2 cannot reach rockfish from inside a Slurm job, so
`parameters/ff_3.0_trained_rf/sidechain.h5` has to be delivered by hand. If it is absent when the
arm test builds, arm R is skipped and the chain runs a single-arm TM4 check, which is still the
pre-production check that was missing. That is logged loudly, never silently.

## The chain (each step submitted as a Slurm dependency of the previous)

1. `srun_mdw2.sh` — `STEPS_PER_JOB=130` steps per job, inside the 36 h wall.
2. `check_continue.sbatch` — resubmits until `MAX_STEPS=500`, then extracts the FF and hands off.
   Its step counter now derives the step from the newest **written** checkpoint's directory name
   (`XX*38 + YY + 1`), not from a file count. See Known Errors.
3. `extract_ff.py` → `parameters/ff_3.0_trained/{sidechain,environment}.h5`.
4. `compare_ff.py` — a free numerical diff of the two trained `sidechain.h5` files, logged before
   anything is launched. Costs nothing and is the first evidence about convergence.
5. `run_arm_test.sbatch` — build arm M (always) and arm R (if the rockfish FF was delivered), both
   from copies of the same seed with the coverage nodes injected, verify both, launch both 12 h.
6. `decide_arm.py` — health gate, then the winner, with midway2 primary on a tie.
7. `decide_and_launch.sbatch` — install the winner into the four real seeds and launch glpG
   production. **glpG production is idle** (all four variants COMPLETED 2026-09-05), so this is now
   a launch rather than a cancel-and-relaunch. See "The baseline data" below.
8. `push_progress.sh` — commit + push so the result is readable off-cluster.

## The baseline data, and the disk bound that decides it

`decide_and_launch.sbatch` was written when production was live, and its relaunch step does
`rm -f $V/$V.run.*.up`. Those files are no longer a half-finished run: they are the **completed
pre-ff3 baseline**, 28 replicas x ~825 MB x 4 variants = ~92 GB, and the derived HDX analysis in
`$V/hdx/` (53 GB) was computed from them.

**Measured growth rate** (`run.0` of two variants, 2026-09-08): a replica is
`280 MB + 43.7 KB x frames`, fitted from 79HIS (51 output groups, 13117 frames, 853 MB) against
79ALA (32 groups, 9730 frames, 705 MB). At `REMD_FRAMES_PER_CHUNK=300` that is ~13 MB per chunk per
replica, so **~364 MB per chunk per variant, ~1.6 GB/h across all four** at the ~1.1 chunk/h the
first campaign sustained.

**That rate is the binding constraint, not `REMD_MAX_BLOCKS=12`.** Quota is 1.45 T of a 1.64 T hard
group limit, so headroom is 195 GB:

| | GB |
|---|---|
| headroom today | 195 |
| 2-arm test, 12 h each | ~22 |
| glpG campaign, per day across 4 variants | ~38 |
| a full 12-block campaign (12 x 36 h) | **~690, which does not fit under any deletion plan** |

So two things follow, and they are separate decisions:

1. **The old raw ladder has to go for the campaign to start with room.** Keeping it leaves 103 GB,
   which the arm test plus ~2 days of campaign consumes, and then the campaign hits the quota
   mid-write, which is ENOSPC and a corrupt HDF5 file. The proposal is to **move `run.0`** (the
   T=0.70 rung, the one `BASELINE_TM_pre_ff3.txt` measured) to `$V/pre_ff3/`, **delete rungs 1-27**,
   and leave `$V/hdx/` untouched. `run.0` at 3.3 GB total keeps a re-measurable old-FF trajectory
   for any future TM4 / Rg / core-RMSD diff at 3.6% of the disk, and rungs 1-27 exist for
   multi-temperature MBAR whose output is already computed and kept.
2. **`REMD_MAX_BLOCKS=12` is not a reachable target and should be set deliberately.** Even at 284 GB
   reclaimed headroom the campaign runs ~7 days before the quota binds. The TM4 and HDX checks read
   the first blocks, so a lower bound is the honest setting; leaving it at 12 encodes an intention
   the disk cannot honour.

The `scancel mdw2_glpG*` step stays. It is a precondition (nothing may be writing those files when
they are moved), not a workaround, and it is a no-op today.

## NP: rebuilt from scratch, NOT in the Thursday chain

Decided 2026-09-07. NP is **out** of `decide_and_launch.sbatch` and will be rebuilt from scratch on
the trained FF afterwards, for two independent reasons:

* **Disk.** Each of the 6 replicas is ~41 GB carrying 98 accumulated `output_previous_*` groups
  (246 GB total). Resetting `block_count` and appending a fresh 8-block campaign adds ~330 GB and
  crosses the 1.64 T hard group quota, i.e. ENOSPC mid-write and a corrupt HDF5 file.
* **Physics.** `run_np_prod.py` reseeds `/input/pos` from `output/pos[-1]`, so the current
  coordinates are a partly unfolded, partly adsorbed structure produced by the OLD force field
  (Rg reached 230 Å). Patching in place would continue an old-FF trajectory under a new FF.

Build assets survive on the cluster in `NP-1AO6/`: `np_hybrid.py` (FCC gold core, MPA tethering,
ions, hybrid assembly — imported as a library), `build_k190.py`, `build_k190_proximal.py`,
`1AO6_A.pdb`, `orientation_map.txt`. Only the `build_all.py` wrapper was lost with the gitignored
scratchpad. `np_hybrid.py` hardcodes `parameters/ff_2.1/` and must be retargeted at
`parameters/ff_3.0_trained/` for `sidechain.h5`/`environment.h5` only — `hbond.h5`, `sheet` and
`martini.h5` were **not** trained and stay at ff_2.1.

**The rebuild must use whichever coverage recipe wins the arm test**, so NP and glpG finally share
one Hamiltonian. That is the whole point of doing it this way rather than patching.

Unresolved: `prod/orientation_map.txt` records cardinal Euler faces while `build_k190_proximal.py`
builds K190-proximal perturbations, and `remote_jobs.md` / `progress.md` disagree about which the
current replicas are. Must be settled from the files themselves before rebuilding.

## Decisions fixed in advance (cannot be asked on Thursday)

* **Health gate** (both arms): all potentials finite, Rg 15–25 Å, < 5 stretched peptide C–N bonds
  in the final frame, helix fraction computable. An arm that fails cannot win.
* **Winner**: higher `mean(TM1, TM4)` helix fraction. If the gap is < 0.05 (inside scatter),
  **midway2's force field wins**, because the seeds, the install chain, the production data and the
  baseline all live on midway2 and it finished 2 steps ahead. A tie is the *expected* result here
  and it is the good one: it means the two independent training trajectories converged.
* **Only one arm available** (the rockfish FF was not delivered in time): run it alone. A single-arm
  result is still the pre-production TM4 check, and it is reported as n=1 rather than as agreement.
* **Both arms fail** → launch nothing, keep backups, write a report. Never run production on a FF
  that failed its health check.
* Baseline for comparison: `popepopg_REMD_mdw2/BASELINE_TM_pre_ff3.txt`, captured 2026-09-05 from
  the pre-install trajectories with the same sign-validated script.

## The intra-protein split (the principle the arm test is really about)

dry-MARTINI must contribute **nothing** intra-protein — it does not reproduce correct protein
dynamics — so the entire intra-protein Hamiltonian is the Upside core's. Verified active, not
assumed: `/input/hybrid_control` attr `exclude_intra_protein_martini = 1` in the glpG seed, the glpG
live replica and the NP live replica; honoured at `src/martini_hybrid.cpp:384-385,790` and
`src/martini_potential.cpp:716-718`. `run_arm_test.sbatch` and `decide_and_launch.sbatch` assert it
before launching anything.

The consequence: glpG's intra-protein terms are `rama_map_pot`, `hbond_energy`/`protein_hbond`,
`rotamer/pair_interaction`, the bonded springs and `backbone_pairs` — and **nothing else**. Versus a
standard Upside core it is missing `hbond_coverage` and `hbond_coverage_hydrophobe`, both
intra-protein and neither substitutable by MARTINI. Backbone H-bond strength in glpG is therefore a
pure function of backbone dihedrals with no burial dependence at all (`src/hbond.cpp:585`).
So Arm B is the arm that actually implements the intended model, not merely the self-consistent one.

Also verified (documented past failure, checked by measurement not by reading the flag):
`current_stage = production` and `rama_map_potential std = 69.95` in the live replica, so the protein
is genuinely mobile despite `preprod_protein_mode = rigid_body`. A frozen protein would give 0.000
and would make the arm test meaningless.

Caveat on Arm B: 42 of glpG's 210 residues (all GLY and ALA) carry no sidechain bead in the hybrid,
so they contribute no coverage. GLY is physically correct to omit; ALA is a hybrid limitation. Arm B
is "coverage from ~80% of residues", not a full standard-Upside coverage term. Say so in any writeup.

## Why an arm test at all

glpG and NP are both hybrids but run **different Hamiltonians**, by build accident not design. The
in-tree hybrid prep never calls `write_rotamer_backbone`/`write_environment` and has no flag to do
so; NP's configs came from the RD1 `envfull` recipe via a gitignored script (`remote_jobs.md`,
"envfull+300Å rebuild"). Unifying them is the goal; the test decides in which direction.

RD1 (below) rejected these terms, but tested them with the **old** FF where coverage was not
co-adapted to the pair term. With a jointly-trained FF that argument does not carry over: 56% of the
trained parameters (coverage 22.4% + hydrophobe 33.6%) are tables Arm A discards, and they moved
*more* during training than the pair term did.

---

# Standing HDX objective (paused while the FF is retrained)

The wildtype (79HIS) ΔG-vs-residue plot is delivered. The other three variants need several more REMD
blocks before their ΔG is meaningful; that sampling restarts on whichever arm wins the test above, so
nothing accumulated before the FF install carries forward.

- [ ] **G2** Identify the second cause of the +2–3% `avg_kinetic_energy/1.5kT` excess (dt-independent,
      present in 1rkl/1AFO).
- Before ANY cluster-side HDX run, re-upload `calc_hdx_ht.py` and `4.calc_D_uptake.py` — the cluster
  copies drift from the repo and will silently return uniform MBAR weights (findings 91).
- `HDX_LIVE=1` is mandatory while REMD is running: the driver rotates `/output` to
  `output_previous_<n>` every chunk, so the live group can vanish mid-read.

Do NOT: change dt for glpG (hard-locked to 0.009; brownian friction is tuned against it), change masses,
widen `destroyed()` thresholds, or add any guard.

---

# Architecture & Key Decisions

- **Master reference for diffs**: `/Users/yinhan/Documents/upside2-md-master` (entire martini subsystem
  is Clean-Slate scope; keep it impeccably clean, remove dead code).
- **dryMARTINI interface** (C++): `martini_potential`, `martini_hybrid`, `martini_brownian`, `martini_masses`,
  `martini_fix_rigid`, `martini_stage_params`. CGL removed.
- **Integrator**: glpG uses one all-particle g-JF step per `.009` numerical timestep; NP uses pure
  velocity-Verlet at dt=0.001.
- **Unit contract**: native dry-MARTINI → Upside conversion happens ONCE at Python h5-build; runtime h5
  and config store Upside-unit values; C++ engine does ZERO unit conversion.
- **H5 FF**: `parameters/ff_2.1/martini.h5` (`/particles` + `/sc_table`). No version numbers; back up
  and overwrite.
- **Spline table**: must equal the published dry-MARTINI functional form (reaction-field electrostatics
  ε_r=15, ε_rf=0; potential-shifted LJ; both reach zero at 1.2 nm). Verified by
  `scratchpad/verify_table_matches_drymartini.py`.
- **No barostat for micelle morphology** (`npt_enable=0`). Bilayer path runs NVT at a target APL derived
  from the reference (Robertson et al. 61.7 Å² for POPE:POPG 2:1).
- **NP campaign**: six independent velocity-Verlet trajectories at T=0.8647, dt=0.001. Self-resubmits
  up to NP_MAX_BLOCKS=8. Gate: ≥5 peptide C–N bonds > 2.0 Å → DESTROYED, no resubmit.
- **glpG campaign**: 48-replica T-REMD (T 0.70–0.90, dt=0.009, dt HARD-LOCKED). Gate: non-finite
  potential OR ≥5 stretched bonds → DESTROYED, no resubmit.
- **REMD momentum rescaling**: coord_swap exchanges pos AND mom, rescaling each by sqrt(T_dest/T_src).
  A finiteness guard is NOT added (NO GUARDS rule; it would hide the defect it is meant to catch).
- **Protein presents BB only** to the MARTINI pair table. SC-env is active and must never be disabled.
  All intra-protein interactions are handled by the Upside core FF.
- **glpG environment morphology**: derived from ITP acyl-chain count (one tail → micelle; two or more →
  bilayer). DDM → micelle with 186 molecules wrapping the 28.2 Å hydrophobic belt; barostat off.
- **Bilayer path**: NVT at target APL; tile/carve geometry; xy-barostat kept for CHARMM-GUI-derived
  systems only, until a trusted target APL exists for those lipids.

## Resolved 2026-09-08: the outage detour is over, and it left two trainers

The midway2 GPFS outage (2026-09-07 ~17:42 CDT, training wedged at step 338) forced training onto
the Mac at step 269 and then onto rockfish at step 275. Storage came back 2026-09-08 morning and
midway2 resumed from its own step-338 checkpoint, so both hosts are now training toward 500 and the
Mac run is stopped. Details in `remote_jobs.md`.

What the detour cost and left behind:

* **69 steps were forfeited once and then were not** — midway2's 338 survived on disk, so its
  lineage never lost them; only the rockfish lineage carries the gap, which is why the two differ.
* **A real defect was fixed to enable local running**: Upside trapped on exit under clang whenever
  Monte Carlo was enabled, so every local worker reported `WORKER_FAIL` on valid data
  (`findings.md` §3.9). Verified bit-identical, so the fix is a keeper regardless of the detour.
* **The reference artifacts are stranded.** `scratchpad/ff3_retraining/` and `training/gly-sym/`
  live on the machine that ran the Mac trainer, not on this one; this checkout received the Sep 7-8
  work as `git pull` fast-forwards, which do not carry a gitignored directory. From here the
  cluster copies are the only readable copies of the chain scripts.

## Known Errors / Blockers

- **`check_continue.sbatch`'s step counter is 2 low, and it decides when to install.** `count_steps`
  counts `checkpoint.pkl` files under `run_output`, but `epoch_02_minibatch_09` and
  `epoch_02_minibatch_14` never wrote one, so the count trails the true step permanently. At true
  step 500 it reads 498, concludes training is unfinished, and submits a 4-node job for 2 more
  minibatches — overshooting to 502 and making the two hosts' force fields come from different
  steps, which would wreck the n=2 comparison. Fix: derive the step from the newest written
  checkpoint's directory name (`XX*38 + YY + 1`), the arithmetic `continue_mdw2.sbatch` and
  `continue_rf.sbatch` already use.
- avg_kinetic_energy/1.5kT is +2.1% above 1.000 after the findings-88 fix. dt-independent; present
  in 1rkl/1AFO. Second cause unidentified (G2 open).
- Molecular DOPC diffusion is not matched at the 40 ps/step clock (measured: 0.015 µm²/s vs 11.5 µm²/s
  target). Fallback is explicitly particle-level friction. Not a blocker for REMD equilibrium sampling.
- **NP albumin over-unfolds and does not reproduce the paper.** Six K190-proximal orientations (runs 0–5, block 3 running). Rg reaches 230.9 Å on run.3 (200 Å box — self-interaction through PBC). Only 3.2% of block-2 frames are adsorbed-and-compact, and none of the paper's five target lysines (K12, K73, K190, K525, K541) are contacted. The paper's central claim (K190 most protected) is contradicted. No footprint npz exists for block-3 data; run `np_footprint.py` after block 3 completes. Larger box required for meaningful structural conclusions on the spread state.
- R4 (CLC-ec1 monomer+dimer on the validated bilayer) is deferred; not scheduled.
- **NaN trigger unidentified (blocker for all glpG-DDM production).** Blow-up origin located and the
  propagation mechanism explained, but nothing measured accounts for a pair crossing from >= 3 A (~500 kT
  margin) into the catastrophic core region. Needs per-step instrumentation inside a running ladder;
  stored trajectories cannot resolve it (60-step frames, no momenta). See findings 90.
- **LJ core table floor: fixed, with a residual.** The `r = max(r, 0.1*sig)` floor is gone and `r_min` is
  0.3 Å, asserted against the analytic form (findings 92). Residual: the clamped spline still flattens
  below the 0.3 Å inner knot, and that domain **was entered** — 6 approaches under 0.3 Å on the corrected
  table (findings 93). Removing the floor changed the consequence, not the entry mechanism.
- **`martini_hdx_project.py` energy contract: settled in favour of the README.** `write_hybrid_energy.py`
  overwrites each replica's `Energy.npy` with the full coupled potential, referenced to the pooled mean.
  The projector's protein-only re-scoring is left alone for the non-hybrid path.


---

# RD1 (2026-08-15): restore the protein-protein terms — rejected then, being retested now

Kept because the current arm test is a direct successor and cites these numbers.

**Then:** three arms (`base` / `env` / `envfull`) from one identical T=0.70 starting configuration,
2 M steps each, on the CB-corrected placement. Helical-core CA-RMSD came out **base 4.61 ± 0.09,
env 4.71 ± 0.12, envfull 4.53 ± 0.08 Å** against a ~2.6 Å target — every arm inside run-to-run
scatter, so the terms were rejected and the production prep left unchanged (findings 103).
`envfull` restored `hbond_coverage`, `hbond_coverage_hydrophobe`, `environment_coverage_hb`,
`hb_environment_coverage_hn/oc`, `bb_sigmoid_coupling_environment`, `hbbb_coverage`,
`cat_pos_bb_coverage`, `placement_fixed_point_vector_scalar`. It was stable (`avg_kinetic_energy/1.5kT`
1.009, 0 broken bonds, 0 non-finite) and needed no C++ change, since `RotamerSidechain` sums a
variable-length `prob_nodes` vector.

**Why it is being retested:** RD1 used the OLD force field, where the coverage tables were not
co-adapted to the pair term. The retrained FF optimises pair/coverage/hydrophobe jointly, so
"restoring these terms does not help" does not carry over — see the CURRENT PHASE section. The new
Arm B is narrower than `envfull`: it restores only `hbond_coverage` + `hbond_coverage_hydrophobe`
(plus their prerequisite placement node), not the environment terms, which dry-MARTINI legitimately
replaces.

**Underlying measurement that motivated RD1** (findings 100/101): glpG's helical-core CA-RMSD plateaus
at 4.15 Å in POPE/POPG vs 2.61 Å in DDM, with backbone H-bond occupancy 0.844 where crystal geometry
scores ~1.0. Ruled out by measurement: integrator, H-bond assignment (agrees with DSSP to 8%), lipid
voids, hydrophobic mismatch, burial threshold.

**Risk carried into the new test:** `environment.h5` was trained against implicit solvent, so an
uncovered residue reads as water-exposed, and the hybrid has no `membrane.h5` to correct that. Failure
mode is over-burying the TM surface — watch for Rg collapsing below the crystal's 20.4 Å.
