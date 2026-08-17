# CURRENT PHASE (2026-08-16): deliver the glpG HDX ΔG-vs-residue plot before Tuesday 2026-08-18

The BB proxy rework (findings 94) is done, deployed and verified; the phase is now the deliverable. Target
shape is `~/Downloads/79HIS_0.90.png`: six broad peaks reaching ΔG 10–20 kcal/mol at residues ~35–48,
~85–100, ~110–125, ~135–148, ~165–175, ~190–200.

**Ladder is fine; the profile shape was limited by an analysis omission, and the fold defect is a small
residual (findings 113).** MBAR reweights correctly (f_k spread 106.9, neighbour overlap 0.093–0.268, better
than the 48-replica DDM ladder). The dG ceiling near 6.8 kcal/mol is master's deliberate resolution guard,
**not** a bug — findings 96 claimed otherwise and is retracted (findings 104).

The dominant cause of "helix amides at ΔG ≈ 2 with isolated needles instead of a contiguous +∞ block" is
that the protection state has **no membrane term**: `get_protection_state.py` scores only protein burial,
the hybrid HDX topology has no `surface` node for master's `--use-TM-region`, and the replacement
(`py/martini_hdx_membrane_accessibility.py` → `combine_hdx_protection.py --water-accessibility`) is never
called by `hdx.sh`. Measured unweighted on one replica: applying it takes bilayer-embedded helical amides at
+∞ from 35/92 to 79/92 and turns needles into 14–17-residue blocks, i.e. the reference figure's shape.
Findings 106 item 4 sized this at 6 amides and is retracted.

What remains is a genuine but much smaller fold defect: 13 of 92 deep helical amides stay finite at
2.3–4.2 kcal/mol, consistent with helical H-bond occupancy ~0.76 against stock's 0.950 (findings 108–112).
Its cause is the environment coupling, not the integrator — reducing dt fourfold removed the temperature
defect entirely and left occupancy and cooperativity unchanged (findings 112). Eight candidates eliminated
by measurement (integrator, timestep, H-bond assignment, lipid packing, hydrophobic mismatch, burial
threshold, CB placement, absent protein–protein terms — RD1, rejected, findings 103). The one untested
node-level difference from a standard Upside config is the rotamer one-body representation: fixed
(`placement_fixed_scalar`) versus rama-dependent (`placement_scalar`).

- [x] **D1** Rework the BB proxy onto `infer_H_O` + constant-weight split; delete the old placement path.
- [x] **D2** Deploy to midway3, rebuild, migrate every config (`py/martini_upgrade_hybrid_args.py`).
      All four glpG jobs and the six NP replicas confirmed loading on the new binary.
- [x] **D3** Validate the HDX pipeline end-to-end and measure MBAR overlap before committing the weekend.
- [x] **D4** Local wildtype REMD, 16 replicas. **Stopped 2026-08-15 at the user's call** with enough
      sampling: 11 533 frames/replica (~692 k steps), 184 528 pooled, across 4 segments, 0 NaN at stop.
      Final figure delivered: 187/203 residues resolved at T = 0.85, 16 off-scale, ΔG −2.4…+18.6, six-peak
      structure matching the reference's positions.
- [x] **D5** Fixed the 0.568 Å CB placement error (findings 102), verified against real crystal CB atoms
      (0.047 Å vs 0.576 Å), and deployed it everywhere.
- [ ] **D6** Local wildtype ladder on the CB-corrected seed — `scratchpad/local_popg_cbfix/`, the Tuesday
      deliverable. Resumed in segments around the user's shutdowns; ~346 k steps/replica after the segment
      ending 2026-08-16 14:30.
- [ ] **D7** Cluster REMD relaunched from CB-corrected seeds 2026-08-16 (53410263–66), back at block 1 of
      12. Gives the other three variants and a 48-rung cross-check; **not** a Tuesday contributor.

**Measured throughput, which settles which dataset the wildtype figure comes from.** Cluster calibration
gives 613.7 ms/step for 48 replicas = 12.8 ms per replica-step; the local machine does 9.00 steps/s for 16
replicas = 6.9 ms per replica-step. So the local run is ~1.9x faster per replica-step and, running 40 h,
reaches **1.3 M steps/replica against the cluster block's ~206 k — 6.3x longer trajectories**, at
comparable pooled frames (347 k vs 245 k). Since the sentinel problem is a shortage of sampled opening
events rather than of ladder rungs (overlap already verified adequate at 16 rungs), trajectory length is
what matters and the local run wins on it. An earlier note in this file calling the cluster "better data
by a wide margin" was wrong — it assumed 48 replicas meant proportionally more sampling.
- [ ] **D8** Membrane term wired into the HDX pipeline (findings 113). `hdx.sh` now runs
      `martini_hdx_membrane_accessibility.py` per replica off the joined hybrid trajectory and folds it in
      with `combine_hdx_protection.py --water-accessibility`. The criterion carries **no free parameter**:
      radius = flat first minimum of the intermolecular tail–tail g(r) (7.00 Å, identical on all 16
      replicas), threshold = 1 contact, which is where a phosphate bead sits (median 0 tail neighbours at
      that radius). `--cutoff`/`--min-contacts` removed. Full 16-replica re-analysis running into
      `scratchpad/local_popg_cbfix/hdx_acc/`; the old protein-only results are kept alongside in `hdx/`.
- [ ] **G2** Identify the second cause of the +2–3% avg_kinetic_energy/1.5kT excess (dt-independent).

Do NOT: change dt for glpG (hard-locked to 0.009; brownian friction tuned against it), change masses,
widen destroyed() thresholds, or add any guard.

## Secondary: NP footprint — answered, see findings 95

K190 is **not** favoured (0.000 contact whenever the protein is still compact, 50th-percentile burial);
K525/K541 are, and "opens and exposes its center" holds. But 71% of the compact window is one orientation,
so it is one binding pose, not a preference. Blocked from going further by over-unfolding, below.

---

# How to run, resume, and analyse

The wildtype HDX deliverable. `scratchpad/local_popg_cbfix/` is the CB-corrected ladder (findings 102);
`scratchpad/local_popg_79HIS/` holds the earlier dataset behind the delivered figure and all the scripts.

```bash
# resume the ladder after any stop (shutdown, blow-up, or a finished segment)
python3 scratchpad/local_popg_79HIS/reseed.py scratchpad/local_popg_cbfix/glpG-RKRK-79HIS.run.*.up
bash scratchpad/local_popg_79HIS/launch_cbfix_ladder.sh <steps>      # ~34 700 steps/h, 16 replicas

# figure from whatever has accumulated; joins all segments itself
HDX_DATA=scratchpad/local_popg_cbfix HDX_WORK=scratchpad/local_popg_cbfix/hdx \
  bash scratchpad/local_popg_79HIS/hdx.sh
```

`reseed.py` rotates each replica's `/output` to the next `output_previous_<n>` and restarts from its last
healthy frame, so segments accumulate and a resume loses nothing. Add `HDX_LIVE=1` to the analysis only if
a run might rotate its `/output` mid-analysis. Sizing a segment to a wall-clock deadline is the norm here:
pick `<steps>` = hours x 34 700 so it ends on its own rather than being killed.

Cluster HDX for all four variants (needs several blocks first — protection-state variance is zero early in
a chain), and the SSH socket needs a Duo push the user approves:

```bash
expect scratchpad/mdw3_master.exp                                   # USER approves Duo
B=/home/yinhanw/project/popepopg_REMD
for V in glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T; do
  sbatch --job-name=hdx_$V --output=$B/logs/hdx.%j.out --partition=caslake --account=pi-trsosnic \
    --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --time=04:00:00 \
    --export=ALL,PDB_ID=$V,HDX_LIVE=1 $B/hdx_cluster.sbatch
done
```

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

## Known Errors / Blockers

- avg_kinetic_energy/1.5kT is +2.1% above 1.000 after the findings-88 fix. dt-independent; present
  in 1rkl/1AFO. Second cause unidentified (G2 open).
- Molecular DOPC diffusion is not matched at the 40 ps/step clock (measured: 0.015 µm²/s vs 11.5 µm²/s
  target). Fallback is explicitly particle-level friction. Not a blocker for REMD equilibrium sampling.
- **NP albumin over-unfolds past the experimental regime.** Rg 26.3 → 39/50/86/103/152/172 Å across the
  six orientations, against a CD measurement showing secondary structure largely retained. Orientations 1
  and 3 (Rg 152/172 Å) exceed the **200 Å box** and self-interact through the periodic image, so they are
  unusable for structural conclusions. Only 3.2% of frames are adsorbed-and-compact, and 71% of those come
  from orientation 2 — which is why the footprint can identify a pose but not a site preference. Testing
  the paper's claim properly needs the adsorbed-but-folded state sampled: a larger box at minimum.
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

# Completed work (condensed)

**2026-08-10**: findings-88 fix (P1/P2); seeds rebuilt; all campaigns resubmitted.  
**2026-08-05**: NP ion convention settled (counterions-only, 218 K+); corrected spline table (reaction-field
  + potential-shifted LJ); glpG micelle REMD launched on fixed binary.  
**2026-08-04**: DDM environment rebuilt as a micelle (186 DDM, tail core 45 Å over 28 Å belt); DOPC
  bilayer regression passed; all four glpG variants submitted.  
**2026-07-20**: g-JF single-step lipid integrator; BB-env force regression found and fixed (partial BB
  routing 14/54, 12/54, 12/54 to N/CA/C; O share disposable); HDX compatibility layer; manuscript
  rewritten; cell-list pairlist; unit-baking to Python h5-build.

---

# Revised Decisions — proposed, NOT yet implemented

## RD1 (2026-08-15): restore the protein-protein environment terms -- TESTED AND REJECTED (findings 103)

**Status: REJECTED 2026-08-16. Both phases built and tested on the CB-corrected placement; neither
moved the helical-core RMSD beyond run-to-run scatter (base 4.61, env 4.71, envfull 4.53 A; target ~2.6).
The production prep is UNCHANGED and should stay that way. See findings 103.** (Architecture Change Rule:
this alters the rotamer solver's inputs, a core pattern).

**Why.** Measured (findings 100/101): glpG's helical-core CA-RMSD plateaus at 4.15 A in POPE/POPG against
2.61 A in DDM, and backbone H-bond occupancy is 0.844 where the crystal geometry scores ~1.0. Both
protection terms then fail together in 3.34% of frames, giving dG = 1.99 kcal/mol -- which reproduces the
observed HDX median exactly. Ruled out by measurement: integrator, H-bond assignment (agrees with DSSP to
8%), lipid voids (none), hydrophobic mismatch (-2.8 A), and the burial threshold.

**What is missing.** MARTINI supplies only protein-environment interactions, so three protein-protein
terms have no substitute: `sigmoid_coupling_environment` (many-body protein self-burial),
`bb_sigmoid_coupling_environment` + `hb_environment_coverage_hn/oc` (backbone burial coupling), and
`hbond_coverage` / `hbond_coverage_hydrophobe` (sidechain-backbone H-bond competition, solved inside the
rotamer node). The sidechain 1-body field, by contrast, was deliberately and correctly replaced by
`martini_sc_table_1body`. `membrane.h5` is correctly absent -- that is the implicit bilayer.

**The open design question.** Standard Upside feeds `rotamer` two coverage nodes; the hybrid feeds it one
MARTINI node. Restoring the terms means deciding how the two 1-body fields compose -- summed into one
input, or the C++ `rotamer` node extended to accept both. Composition must not double-count: Upside
coverage is over protein CB neighbours, MARTINI over environment beads, so they are disjoint by
construction, but this needs verifying in the C++ rather than assuming.

**Risk to watch.** `environment.h5` was trained against implicit solvent, so an uncovered residue reads as
water-exposed; master corrects that in membranes with `membrane.h5`, which the hybrid does not have. The
failure mode is over-burying the TM surface. Accept the change only if helical-core RMSD moves toward
2.6 A **without** Rg collapsing below the crystal's 20.4 A (DDM already sits at 18.6 A, i.e. compacted).

**Cheap prior.** The self-burial term disfavours the drifted state by only 8.0 E_up = 5.6 kcal/mol, so it
is probably not sufficient alone; the H-bond coverage coupling is the likelier dominant piece.

**Test as built** — `scratchpad/rd1_env_test/`, three arms from one identical starting configuration
(the finished production run's T=0.70 endpoint), 2 M steps each at T=0.70:

| arm | config | nodes restored |
|---|---|---|
| `base` | hybrid as shipped | none (control) |
| `env` | + phase 1 | `weighted_pos`, `environment_coverage_sc`, `sigmoid_coupling_environment` |
| `envfull` | + phase 2 | also `hbond_coverage`, `hbond_coverage_hydrophobe`, `environment_coverage_hb`, `hb_environment_coverage_hn/oc`, `bb_sigmoid_coupling_environment`, `hbbb_coverage`, `cat_pos_bb_coverage`, `placement_fixed_point_vector_scalar` |

`envfull` differs from a standard Upside config by exactly one node, `placement_scalar`, whose hybrid
counterpart `placement_fixed_scalar` is present (fixed vs dynamic rotamer probabilities — a rotamer
representation choice, not an environment term). `membrane.h5` stays out. Built by calling Upside's own
writers, not reimplementing them; `hbond_coverage`/`hbond_coverage_hydrophobe` are appended to the
rotamer's argument list, which needs no C++ change because `RotamerSidechain` sums a variable-length
`prob_nodes` vector.

Pre-flight (envfull, 3000 steps): `avg_kinetic_energy/1.5kT` 1.009, potential -24110 -> -24368, Rg
19.44 -> 19.62, 0 broken bonds, 0 non-finite. Stable.

Verdict via `compare_arms.py`, which encodes the acceptance rule. **Do not read it before the arms are at
comparable step counts** — all three start at core-RMSD ~4.3 and the original run needed ~150 k steps to
reorganise.

**If `envfull` wins**, the follow-on is to move the injection into `martini_prepare_system_lib.py`, rebuild
the production seeds, and relaunch the REMD ladder. That is roughly a day of sampling, so it only reaches
a figure by Tuesday if started Sunday night.
