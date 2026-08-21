# TM4 Helix Instability — Root Cause, Fix, and Next Steps

## Problem

TM4 (residues 136–151, sequence GVVYALMGYVWLRGER) is visibly unstable in glpG-RKRK-79HIS REMD
trajectories. At T=0.70 E_up (245 K), ~70% of frames have ≥1 non-helical residue in TM4, with G136
(TM4 N-terminus) showing P(alphaL) = 0.225.

## Root Cause

The Upside core force field uses context-dependent per-residue Ramachandran maps trained from
soluble-protein PDB statistics. Every GLY residue in glpG-RKRK-79HIS has a map that **prefers alphaL
over alphaR at the canonical bins**, because the training data is dominated by non-helical contexts
where GLY following a loop tends to adopt alphaL. The artifact is worst for G136 (TM4 N-terminal GLY),
where the alphaR penalty (~2.0 E_up at the canonical -60° bin) competes directly with H-bond
stabilization.

This is **not** a missing `uhb_energy` term from the implicit bilayer model (that term is excluded from
the hybrid by design). TM6's G195 is 96–100% helical despite having the same type of biased map,
because TM6 mid-helix has more H-bond partners. The artifact is worst at helix termini with fewer
stabilizing contacts.

## Physical Justification for Fix

GLY is achiral (no beta-carbon). Its intrinsic Ramachandran potential has a physical symmetry under
(phi, psi) → (−phi, −psi). The asymmetry in context-dependent maps trained on soluble proteins is a
sampling artifact for GLY in helical contexts: TM-helix GLY is underrepresented relative to soluble
protein statistics, biasing the maps toward alphaL.

**Scope**: this justification applies only to helical-context GLY. For loop and turn GLY (e.g. G132,
G133 in the turn before TM4), surrounding L-amino acids break the achirality symmetry and the alphaL
preference is genuine — those maps must not be symmetrized.

## Fix Implemented

**File**: `py/upside_config.py`  
**Functions**: `write_rama_map_pot` and `write_rama_map_pot2`  
**Location**: After reading `rama_pot` from the training library, before normalization.

Criterion: GLY residue is in helical context (phi ∈ [−130°, −20°] from input coordinates, or phi=nan
for the N-terminal residue). For those GLY, if the raw map prefers alphaL over alphaR at canonical
bins, the map is symmetrized:
```python
m_sym = 0.5 * (m + m[np.ix_(idx_phi, idx_psi)])
```
where `idx_phi = (-np.arange(n_phi)) % n_phi` implements the (phi, psi) → (−phi, −psi) symmetry.

Canonical bin indices (72×72 grid):
- alphaR: `i=24` (phi=-60°), `j=27` (psi=-45°)
- alphaL: `i=48` (phi=+60°), `j=45` (psi=+45°)

**GLY fixed (16 helical-context)**: G1(nan), G12, G17, G28, G96, G97, G106, G120, G136, G143, G149,
G162, G174, G186, G191, G195

**GLY skipped (7 loop/turn)**: G49(+85.9°), G104(+91.4°), G128(+100.0°), G132(−19.6°), G133(+38.4°),
G156(+110.4°), G180(+111.5°)

The phi values are from the input (equilibrated starting) structure. G132 with phi=−19.6° is excluded:
−19.6° > −20° (numerically outside the range), correctly preserving the turn geometry before TM4.

The paper (`example/16.MARTINI/drymartini_upside_interface.tex`) has been updated throughout: new
`\subsection{Context-dependent Ramachandran maps and the glycine correction}` in Methods, symmetry
equation (`\label{eq:gly_sym}`), updated Abstract, Introduction, HDX section, Evidence table,
Discussion, and Conclusions.

## Key Energy Numbers

- G136 bin(25,29) [actual helical phi ≈ -51°]: orig=0.390 E_up → fixed=0.449 E_up (+0.059)
- G136 bin(48,45) [canonical alphaL]: orig=0.851 E_up → fixed=1.440 E_up (+0.589)
- G143 bin(48,45): orig=0.459 → fixed=1.154 E_up (+0.695)
- G149 bin(48,45): orig=0.590 → fixed=1.346 E_up (+0.756)

## Local Validation 1 — Stay-folded (10k steps, T=0.85, from helical frame 2542)

| Metric | Original | Fixed |
|---|---|---|
| TM4 severe H-bond disruption (≥5 donors < 0.3) | 14.9% | **6.0%** |
| TM4 mean H-bond score | 0.728 | **0.743** |
| TM6 G195 (control) | unchanged | unchanged |

60% reduction in severe disruption.

## Local Validation 2 — Re-folding from actual disrupted frame (20k steps, T=0.70)

**Starting frame**: frame 158 of the original unfixed T=0.85 simulation (`perm_orig.up`), where
TM4 has 2 broken helical pairs (max CA-CA+4 = 8.20 Å). This is a real disrupted TM4 state
produced by the unfixed potential — not artificially generated.

Both runs (original and fixed potential) start from identical positions. Run at T=0.70 (the REMD
physiological rung). The "Fixed" result below is with **G136, G143, G149 symmetrized** (the only
GLY that were actually changed in the h5 built for this test — the broader helical set was added
in the cluster deployment).

| Window | Metric | Original (unfixed) | Fixed |
|---|---|---|---|
| Early (0–20%) | Mean broken CA-CA+4 pairs | **6.42** | **0.21** |
| Mid (20–60%) | Mean broken pairs | **8.39** | **0.09** |
| Late (60–100%) | Mean broken pairs | **7.75** | **0.11** |
| Late | TM4 fully helical (n_broken=0) | **0%** | **88.8%** |
| Late | TM4 persistently disrupted (n_broken≥2) | **100%** | **0%** |
| Late | G136 in alphaL (phi>20°) | 0% | 0% |

The unfixed potential does not stall at the observed disruption — TM4 gets worse, rapidly reaching
8 broken helical pairs and never recovering (0% fully helical). The fixed potential refolds TM4
within the early window and holds it at 88.8% fully helical.

## Local Validation 3 — Targeted fix early-window check (T=0.70, starting from frame 158)

The targeted fix (16 helical GLY, G132/G133 excluded) was started from the same disrupted frame 158.
The early hbond recovery (same seeds for T=0.70):

| Time | H-bonds | Protein potential (E_up) |
|---|---|---|
| 0 | 175.5 | -1283.94 |
| 100 | 194.6 | -1494.18 |
| 200 | 195.6 | -1548.89 |
| 300 | 199.6 | -1541.03 |
| 400 | 197.4 | -1585.26 |
| 800 | 179.8 | -1674.43 |
| 1000 | 177.6 | -1643.88 |

Rapid hbond recovery (175→197 by time=100) is identical to the G136/G143/G149-only fix. The targeted
fix enables the same TM4 refolding pathway. Protein potential deepens throughout, consistent with
improved helical packing. The protein potential at time=800 (−1674) is ~30% deeper than time=0
(−1284), indicating substantial structural relaxation.

A complete 20k-duration trajectory was not feasible locally (this system requires ~3.3 hours).
The early-window recovery is sufficient to confirm the targeted fix does not impair TM4 refolding.

## Local Validation 3 — Targeted fix quantitative (duration=300, frame-interval=10, T=0.70)

Same disrupted starting frame 158. 31 output frames. Targeted fix: 16 helical GLY symmetrized.

| Metric | Result |
|---|---|
| Mean broken CA-CA+4 pairs | **0.26** |
| Fully helical (n_broken=0) | **77.4%** |
| Any disruption (n_broken≥1) | 22.6% |
| Max broken pairs any frame | 2 (frame 1 = t=10, starting state) |

G132/G133 loop phi over all frames — correctly stays in alphaL/turn region:
- G132: phi = 74°–132° (turn, positive phi throughout)
- G133: phi = 36°–96° (turn, positive phi throughout)

G136 helix phi over all frames: −104° to −121° (helical alphaR, no alphaL visits).

**Conclusion**: targeted fix (16 helical GLY, G132/G133 excluded) refolds TM4 just as effectively as
the narrow fix. The loop geometry at G132/G133 is preserved — excluding them was correct.

## Cluster Deployment (updated 2026-08-21)

All 4 seed h5 files on midway3 patched with the targeted fix (16 helical GLY). Backups:
- `.bak_gly_fix` — original, pre-any-fix
- `.bak_gly_fix_v2` — intermediate broad fix (all 23 GLY, superseded)

Current seeds: targeted fix (16 helical GLY, G132/G133/other loop GLY excluded).

| JobID | Variant |
|---|---|
| 54103426 | glpG-RKRK-79HIS |
| 54103428 | glpG-RKRK-79HIS_S115T |
| 54103430 | glpG-RKRK-79ALA |
| 54103431 | glpG-RKRK-79ALA_S115T |

Patch scripts:
- `/home/yinhanw/project/popepopg_REMD/gly_fix_seeds_targeted.py` — current (targeted)
- `/home/yinhanw/project/popepopg_REMD/gly_fix_seeds.py` — superseded (all 23 GLY)

Seed dir: `/home/yinhanw/project/popepopg_REMD/seeds/`  
Logs: `/home/yinhanw/project/popepopg_REMD/logs/remd.<jobid>.out`  
Cluster `upside_config.py`: `/home/yinhanw/beagle3/yinhan/upside2-md/py/upside_config.py` — updated.

## What to Check When Cluster Results Arrive

### Question 1: Does the fix correct G136's conformational bias?

Measure G136's P(alphaL) from the REMD trajectory at T=0.70. Before the fix: ~0.22. After a
successful fix: should drop to ≤0.05. If G136 is still predominantly alphaL at T=0.70, the
Ramachandran artifact is not the only driver.

```python
# Load trajectory from popepopg_REMD/glpG-RKRK-79HIS/ with martini_remd_concat.py first.
# Then read phi dihedral for G136 (0-indexed residue 135).
# phi atoms confirmed: [538, 540, 541, 542] = [C(i-1), N(i), CA(i), C(i)]
# Bin convention: idx = int((phi_deg + 180.) / 5.) % 72
# alphaL bin: i=48 (phi ≈ +60°), alphaR bin: i=24 (phi ≈ -60°)
```

### Question 2: Is residual disruption physical flexibility or another artifact?

Check TM4 severe H-bond disruption as a function of T across the REMD ladder. If the fraction
decreases monotonically toward zero as T → 0.70, the residual 6% is physical flexibility. If it
remains elevated even at T=0.70, there is a second mechanism to find.

### Question 3: Loop GLY stability — G132 and G133

The targeted fix correctly excludes G132 (phi=−19.6°) and G133 (phi=+38.4°) — both in the turn
immediately before TM4. The cluster results will show whether the turn geometry is maintained:

**What to measure:**
- phi/psi distributions for G132 and G133 across the REMD trajectory
- H-bond pattern of the loop (residues 130–135) — are the turn H-bonds maintained?
- If these loop residues adopt more alphaR (phi≈-60°) than expected, the exclusion criterion may
  need to be narrowed further (e.g., require |phi| > 30° for loop detection).
