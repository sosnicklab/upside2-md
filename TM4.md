# TM4 Helix Instability — Root Cause, Fix, and Next Steps

## Problem

TM4 (residues 136–151, sequence GVVYALMGYVWLRGER) is visibly unstable in glpG-RKRK-79HIS REMD
trajectories. At T=0.70 E_up (245 K), ~70% of frames have ≥1 non-helical residue in TM4, with G136
(TM4 N-terminus) showing P(alphaL) = 0.225.

## Root Cause

The Upside core force field uses context-dependent per-residue Ramachandran maps trained from
soluble-protein PDB statistics. Every GLY residue in glpG-RKRK-79HIS has a map that **prefers alphaL
over alphaR at the canonical bins**, because in soluble-protein data, GLY at a helix N-terminus
following a loop is predominantly seen in non-helical conformations. The artifact is worst for G136
(TM4 N-terminal GLY following a loop), where the alphaR penalty (~2.0 E_up at the canonical -60° bin)
competes directly with H-bond stabilization.

This is **not** a missing `uhb_energy` term from the implicit bilayer model (that term is excluded from
the hybrid by design). TM6's G195 is 96–100% helical despite having the same type of biased map,
because TM6 mid-helix has more H-bond partners. The artifact is worst at helix termini with fewer
stabilizing contacts.

## Physical Justification for Fix

GLY is achiral (no beta-carbon). Its intrinsic Ramachandran potential has a physical symmetry under
(phi, psi) → (−phi, −psi). The asymmetry in the trained maps is entirely a sampling artifact.
Symmetrizing removes the artifact without introducing any new physics.

## Fix Implemented

**File**: `py/upside_config.py`  
**Functions**: `write_rama_map_pot` and `write_rama_map_pot2`  
**Location**: After reading `rama_pot` from the training library, before normalization.

For every GLY residue where the raw map prefers alphaL over alphaR (canonical bins), the map is
symmetrized:
```python
m_sym = 0.5 * (m + m[np.ix_(idx_phi, idx_psi)])
```
where `idx_phi = (-np.arange(n_phi)) % n_phi` implements the (phi, psi) → (−phi, −psi) symmetry in
the stored grid convention. The condition `m[i_alphaR, j_alphaR] > m[i_alphaL, j_alphaL]` fires for
**all 23 GLY residues** in glpG-RKRK-79HIS — equivalent to unconditional GLY symmetrization for this
protein.

Canonical bin indices (72×72 grid):
- alphaR: `i=24` (phi=-60°), `j=27` (psi=-45°)
- alphaL: `i=48` (phi=+60°), `j=45` (psi=+45°)

The paper (`example/16.MARTINI/drymartini_upside_interface.tex`) has been updated throughout: new
`\subsection{Context-dependent Ramachandran maps and the glycine correction}` in Methods, symmetry
equation (`\label{eq:gly_sym}`), updated Abstract, Introduction, HDX section, Evidence table,
Discussion, and Conclusions.

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
physiological rung).

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
within the early window and holds it at 88.8% fully helical. This is a clean pass/fail test:
the fix enables refolding from a state that the unfixed potential cannot escape.

## Key Energy Numbers

- G136 bin(25,29) [actual helical phi ≈ -51°]: orig=0.390 E_up → fixed=0.449 E_up (+0.059)
- G136 bin(48,45) [canonical alphaL]: orig=0.851 E_up → fixed=1.440 E_up (+0.589)
- G143 bin(48,45): orig=0.459 → fixed=1.154 E_up (+0.695)
- G149 bin(48,45): orig=0.590 → fixed=1.346 E_up (+0.756)

## Cluster Deployment (done 2026-08-21)

All 4 seed h5 files on midway3 were patched with the same symmetrization logic. Backups at
`seeds/<V>.up.bak_gly_fix`. New REMD jobs submitted:

| JobID | Variant |
|---|---|
| 54095760 | glpG-RKRK-79HIS |
| 54095761 | glpG-RKRK-79HIS_S115T |
| 54095762 | glpG-RKRK-79ALA |
| 54095763 | glpG-RKRK-79ALA_S115T |

Patch script: `/home/yinhanw/project/popepopg_REMD/gly_fix_seeds.py`  
Seed dir: `/home/yinhanw/project/popepopg_REMD/seeds/`  
Logs: `/home/yinhanw/project/popepopg_REMD/logs/remd.<jobid>.out`

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

### Question 3: Loop GLY stability — the main open risk

The symmetrization was applied to all 23 GLY residues, including loop residues G132 and G133 which
are in the turn immediately before TM4. Those residues may genuinely need their alphaL preference to
form the correct turn geometry that caps TM4 at its N-terminus. If the fix over-corrects them, a new
artifact is introduced at the TM4 cap that partially offsets the gain from fixing G136.

**What to measure:**
- phi/psi distributions for G132 and G133 across the REMD trajectory
- H-bond pattern of the loop (residues 130–135) — are the turn H-bonds maintained?
- Protein Rg vs T — if the fix destabilizes loops, Rg at T=0.70 should increase relative to the
  pre-fix REMD data

**If loop GLY residues are destabilized:** the fix must be scoped to TM residues only. At h5 build
time, identify which GLY residues are in TM helices (e.g., by phi occupancy in the input structure
or by secondary structure assignment) and apply the symmetrization only to those. The current
condition (energy criterion alone) is a proxy that works for glpG because all 23 GLY have the
artifact, but it is not guaranteed to be safe in general.

## If a Second Fix Is Needed

If G132/G133 are destabilized, the targeted fix is:

1. At h5 build time in `write_rama_map_pot`, detect TM-helix GLY by checking whether the input
   structure has the residue in a helical phi range (phi ∈ [-100°, -30°]) rather than relying on
   the energy criterion alone.
2. Apply symmetrization only to those GLY, leaving loop/turn GLY untouched.
3. Patch the 4 cluster seeds again (backup with `.bak_gly_fix_v2`), cancel, resubmit.

This is a refinement of the same physical argument — GLY achirality justifies symmetrization only
in the helical context where the soluble-protein training data is unrepresentative. In turns, GLY
genuinely populates alphaL (a real conformation for achiral glycine), so symmetrizing those maps
raises a real preference to match a phantom one.
