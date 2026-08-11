# CURRENT PHASE (2026-08-10): hybrid interface force propagation fix (findings 88)

The glpG/NP instability root cause: `HybridPositionNode::propagate_deriv` discarded the O site's
sensitivity entirely and dropped BB's fourth (O) weight, while all five backbone sites (N, CA, C, O, BB)
sat in the MARTINI pair table. Fixed: (a) restrict the protein side of the pair table to BB; (b) route
BB's O-weight sensitivity through the placement Jacobian. Fix scope is the dryMARTINI interface only —
no Upside core change.

**Status 2026-08-10/11:** P1, P2, P4 done. Seeds rebuilt, both campaigns running.
G1 is unusable as a gate (round-off dominated at this potential scale — see findings 88).
G2 FAILED: kinetic excess still +2.1%; second cause unidentified, not blocking production.

- [x] **P1 (py, prep)** Restrict protein side of MARTINI pair list to BB only.
- [x] **P2 (C++, martini_hybrid.cpp)** Complete BB's gradient via O-placement Jacobian; propagate O sensitivity.
- [ ] **G2** Identify the second cause of the +2–3% avg_kinetic_energy/1.5kT excess (dt-independent; present
      after the fix in 1rkl/1AFO).
- [x] **P4** Rebuild all seeds on the corrected table; resubmit both campaigns.
- [ ] **G3** Localize the 79HIS runaway: rerun from clean seed with dense output to catch the divergence.
      (Low priority — both variants are now running on the fixed binary.)

Do NOT: change dt for glpG (hard-locked to 0.009; brownian friction tuned against it), change masses,
widen destroyed() thresholds, or add any guard.

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
- NP run.3 protein now spans 246 Å in a 200 Å box — can interact with its own periodic image in the
  most extended conformation.
- R4 (CLC-ec1 monomer+dimer on the validated bilayer) is deferred; not scheduled.

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
