# Progress — glpG HDX: REMD-system fix + lipid de-hardcode (2026-07-25)

## Goal
Metad denatures glpG (both CVs) → deprioritized; REMD is the working HDX method.
Fix the REMD prep system, then run 36 h REMD for 4 glpG variants, then investigate metad.

## Done
- **Ion count = molarity × box volume** (removed the template `effective_volume_fraction`
  fudge; deleted `infer_effective_ion_volume_fraction_from_template`). No-op for DDM (795
  pairs / 1596 ions), verified on-cluster. Box kept at 3× (user decision); clustering accepted.
- **Lipid de-hardcode (full):** removed `_LIPID_SPEC` + `UPSIDE_MARTINI_LIPID` env gate;
  lipid identity now via `--lipid-name` CLI; topology discovered by moleculetype name across
  the input ITPs (no hardcoded filename); `--bilayer-pdb` defaults to `dryMARTINI/<lipid>.pdb`;
  `dopc_*`→`lipid_*`; `parse_dopc_from_itp`→`parse_lipid_from_itp(path,molname)`;
  `dopc_max_sigma_nm`→`lipid_max_sigma_nm`. Drivers pass `--lipid-name DDM`.
  (Left `UPSIDE_DOPC_*` diffusion env vars — kinetics calibration, not lipid-type selection.)
- **Clearance bug fixed:** `derive_lipid_contact_clearance_angstrom` reads the ACTIVE lipid's
  ITP → DDM corrected 6.96 Å (DOPC-derived, wrong) → 5.28 Å.
- Validated locally: compile clean, no leftover DOPC symbols, topology parity (DDM/DOPC bead
  types+bonds+angles identical), per-lipid clearance, CLI threading, structure build
  (795 pairs, charge-neutral). Synced 3 files to cluster (md5 match).

## Preps run LOCALLY (caslake congested), uploaded, REMD-only on Slurm
- 4 preps built locally with corrected code (clearance target 5.28 Å confirmed in logs; folded
  production hbonds~200, 795 ions). `.up` at `example/16.MARTINI/scratchpad/local_preps/<V>/checkpoints/`.
- Uploaded to cluster; load-tested + faithful reseed+50-step run on cluster obj/upside → folded/finite
  (locally-built `.up` is cluster-compatible).
- **REMD submitted (caslake, 48cpu exclusive, 35:45:00, no prep queue, self-resubmit ≤12 blocks):**
  79HIS=52647003, 79HIS_S115T=52647004, 79ALA=52647005, 79ALA_S115T=52647006 (pending, Resources).
- Monitor running to confirm first block folded once they start.

## DDM box-fill fix (2026-07-27) — CURRENT
User saw DDM "aggregated with gaps" in the REMD VTF. Root cause: the tiler filled only the
crop window (~3× span COM extent) at native density, but `set_box_from_lipid_xy` set the box to
the DDM **molecule extent** (~187 Å, since molecule tails reach ~30 Å past their COMs) → a sparse
tail-only outer ring; being under-dense, the detergent then contracts → gaps. NOT a viz artifact.
- User decision: keep the box large (≥3× span), add more DDM to fill it.
- Fix (dev-only `martini_prepare_system{,_lib}.py`, not in master): round the tiling window up to
  whole template tiles; add `force_xy_box` to `set_box_from_lipid_xy` so the box = the fill window
  (not molecule extent), COM-centered (tails wrap under PBC). Belt radius `r_lipid` computed
  min-image so PBC-split template molecules don't inflate it (was 67→315 Å box bug; now 2.7→189 Å).
- Verified (structure build, 79HIS): box 189×189 (4.6× span, ≥3×, ≈ old 187), DDM 621 mols (was
  235; "more DDM"), midplane density FLAT 0.074–0.078 to the edge (uniform, no gap/belt).
- In progress: full local prep (min→eq→40k burn-in→prod) to confirm no MD contraction, then re-prep
  all 4 variants locally, upload `.up`, re-run REMD.

## Metadynamics REMOVED from branch (2026-07-27)
Metad confirmed ineffective (denatures glpG; fundamental CV/observable mismatch). Removed from the
branch: deleted `src/metadynamics.{cpp,h}` + both CMake entries; removed the `#include` +
`maybe_deposit` call from main.cpp; removed `write_metadynamics` from advanced_config.py; deleted
`analyze_metadynamics.py`. KEPT (not metad): the NaN blow-up guard (REMD robustness), all
hybrid/MARTINI code, and the PLUMED interface (`write_plumed` is in master). Rebuilt clean; hybrid
runs folded (hbonds 194.5, Rg 20.4, finite) — metad code was dormant for the hybrid, so parity holds.
- **coord_operator.cpp Sum-node change REVERTED to master.** It is a genuine gradient bug
  (`weight[id_pos[i]]` — wrong index + OOB read when id_pos is a subset; inconsistent with
  compute_value's `weight[i]`), but LATENT: no real config builds a Sum node (verified 0 Sum nodes in
  stage_7.0.up; only the metad CV did). Reverting restores C++ core parity, zero behavioral change.
  Latent upstream bug — fix separately in master if the Sum node is ever used.
- **DDM "lipids bend a lot" = physical MARTINI, NOT a bug.** AngleSpring uses the exact MARTINI type-2
  form V=½k(cosθ−cosθ0)² (martini_potential.cpp:1975); k converted by /2.91495; stiff angles sit on
  θ0 (2-4-7: 152.5°/149.7°, 5-4-2: 70.9°/70°). The cosine angle has ZERO restoring force/curvature at
  θ0=180°, so the C12 tail angles (k=11, 77) are floppy by design; DDM itp has no dihedrals. Stiffening
  would violate the no-twisting rule. FF is faithful; CG detergent just looks floppy. No re-run needed.

## Prior state (context)
- REMD is the reliable HDX method.
- REMD crashed at 2/12 blocks (replica instability → NaN frames). ΔG re-run (stageBC, NaN-sanitized)
  jobs 52689487/88/89 running. Sum-node `coord_operator.cpp` fix kept LOCAL only, parity-proven.
  ← these REMD runs used the OLD (gappy) DDM box; superseded by the box-fill re-run.
