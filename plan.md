# CURRENT PHASE (2026-08-15): deliver the glpG HDX ΔG-vs-residue plot before Monday 2026-08-17

The BB proxy rework (findings 94) is done, deployed and verified; the phase is now the deliverable. Target
shape is `~/Downloads/79HIS_0.90.png`: six broad peaks reaching ΔG 10–20 kcal/mol at residues ~35–48,
~85–100, ~110–125, ~135–148, ~165–175, ~190–200.

**The limit is sampled amide opening events, and this is measured, not assumed.** A dry run of the full
pipeline on the first 5% of the local run gives f_k spread 106.9 and neighbour overlap 0.093–0.268
(median 0.233, better than the 48-replica DDM ladder), so the thinned ladder reweights fine. What is
missing is events: 65–112 of 203 amides sit at the p_f = 1 sentinel and ΔG is compressed to −5…+7.

- [x] **D1** Rework the BB proxy onto `infer_H_O` + constant-weight split; delete the old placement path.
- [x] **D2** Deploy to midway3, rebuild, migrate every config (`py/martini_upgrade_hybrid_args.py`).
      All four glpG jobs and the six NP replicas confirmed loading on the new binary.
- [x] **D3** Validate the HDX pipeline end-to-end and measure MBAR overlap before committing the weekend.
- [ ] **D4** Local wildtype REMD, 16 replicas × 1.3 M steps (~40 h, Sunday PM). **Primary source for the
      wildtype figure.** Analyse at any point.
- [ ] **D5** Cluster REMD, 4 variants × 48 replicas, block 1 of 12 due Sunday ~13:00. Supplies the other
      three variants and a 48-rung ladder, not more sampling of the wildtype.

**Measured throughput, which settles which dataset the wildtype figure comes from.** Cluster calibration
gives 613.7 ms/step for 48 replicas = 12.8 ms per replica-step; the local machine does 9.00 steps/s for 16
replicas = 6.9 ms per replica-step. So the local run is ~1.9x faster per replica-step and, running 40 h,
reaches **1.3 M steps/replica against the cluster block's ~206 k — 6.3x longer trajectories**, at
comparable pooled frames (347 k vs 245 k). Since the sentinel problem is a shortage of sampled opening
events rather than of ladder rungs (overlap already verified adequate at 16 rungs), trajectory length is
what matters and the local run wins on it. An earlier note in this file calling the cluster "better data
by a wide margin" was wrong — it assumed 48 replicas meant proportionally more sampling.
- [ ] **G2** Identify the second cause of the +2–3% avg_kinetic_energy/1.5kT excess (dt-independent).

Do NOT: change dt for glpG (hard-locked to 0.009; brownian friction tuned against it), change masses,
widen destroyed() thresholds, or add any guard.

## Secondary: NP footprint — answered, see findings 95

K190 is **not** favoured (0.000 contact whenever the protein is still compact, 50th-percentile burial);
K525/K541 are, and "opens and exposes its center" holds. But 71% of the compact window is one orientation,
so it is one binding pose, not a preference. Blocked from going further by over-unfolding, below.

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
