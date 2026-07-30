# Progress — dryMARTINI hybrid prep redesign (2026-07-30)

Approved plan: CHARMM-GUI-only bilayer + martinize + OPM orientation; xy-only tensionless MC-barostat
NPT; solvent-volume ions. Done so far:
- **box.cpp/box.h** (subagent): removed broken Berendsen/PR barostat + virial + shrink guards; MC
  barostat is now the sole path (COM-based, energy-based, freezes z via mc_dmax_z=0). Clean compile
  (+23/−313). master untouched.
- **py prep rewrite**: `read_charmm_gui_membrane` (strip water, drop ions, CRYST1 required, count
  cross-check); `martinize_protein_cg`; `orient_protein_to_opm` (Biopython BLOSUM62 + Kabsch + ≤4Å
  refine); rewrote `prepare_mixed_structure` (unified OPM orient → same rigid transform to backbone +
  CG envelope → tile → **overlap vs CG envelope** → box → solvent-volume ions in slabs); CLI
  `--charmm-gui-dir/--membrane-pdb/--membrane-top/--protein-cg-pdb/--opm-reference/--protein-pbc-margin`;
  removed retry loop + placement/surface-gap/min-gap/cutoff-step-max dead args.
- **barostat config** (py): equilibration stages write MC attrs (target_p=0 tensionless,
  compressibility_z=0, mc_dmax_z=0 xy-only, mc_dmax_xy, mc_seed); removed dead `type`/`barostat_type`
  threading + `stage_npt_targets`; production NVT (run.py prod_70_npt_enable 1→0).
- **ions**: `estimate_salt_pairs(box, molar, membrane_thickness_z, protein_volume)` = molarity ×
  solvent volume; `place_ions(..., exclude_z=)` restricts to solvent slabs.
- **Verified**: regression_lipid_box_fill PASSES (mode=bilayer byte-identical). Smoke test on 79HIS:
  OPM 163/186 Cα core, box 100×100×153.6 (was 200²), 446 DDM, **106 ion pairs (was 795), 0 ions in
  membrane**, TM in slab (63–91) + soluble domain in solvent (15–62), no water. martinize 487 beads OK.

- **P6 DONE (awaiting user review)**: full local run-hybrid-workflow on 79HIS completed rc=0
  (6.0→6.6→7.0, "Workflow Complete", no NaN). MC barostat: xy 100.0→100.8 across stages, **z frozen
  153.63, no collapse**. DDM relaxed 33→48.8 Å (physical DDM bilayer thickness — NPT did real work).
  Protein stayed folded (hbonds ~200, Rg 20.2). Staged to ~/Downloads: preproduction_6.6.pdb +
  .vtf + packed_initial.pdb. GATE: awaiting user OK before P8.
- Watch items for review: DDM thickening 33→49 Å (physical? eyeball belt); ions started 0-in-membrane
  but ~20/218 diffused into the DDM z-band during MD (near heads vs core?).

- **Review-feedback fixes (2026-07-30)**:
  - "ions xy outside bilayer" = my review PDB dumped Upside's *unwrapped* MD coords. Fixed with
    per-molecule min-image wrapping to the membrane center; verified 0/218 ions outside the DDM
    footprint, 0 in the hydrophobic core (near-membrane ions sit at the polar heads). Not a sim bug.
  - Ion concentration was exactly 0.15 M at prep, but (a) box z was oversized (stacked
    box_padding_z=20 + pbc_margin=15 = 70 A margin) and (b) the membrane relaxes 33->49 A, so the
    realized concentration drifted to ~0.17 M. Fix: box_padding_z default 0 (z-margin = pbc_margin
    only) -> box_z 153.6->123.7; and new `--membrane-thickness-angstrom` counts ions against the
    EQUILIBRATED thickness (48.8 A for DDM) so production realizes 0.1496 M. Ions 218->136.
  - readme.md rewritten as one coherent prep guide; run.py updated to the new CHARMM-GUI+OPM
    interface (was passing removed args).
  - Confirming local re-run (out2) in progress: prep confirms box 123.7 + 65 pairs; equilibration
    folded, no NaN.

- **Rigid-protein check (user request)**: CONFIRMED protein holds rigid through pre-production.
  Mechanism: hybrid_control preprod_protein_mode="rigid_body" -> martini_hybrid.cpp preprod_rigid;
  enforce_preprod_rigid_stage = preprod_rigid && (stage!="production") -> set_dynamic_rigid_groups
  for all 6.0-6.6 + 7.0 handoff/burn-in; released only at "production" relabel. Empirical: hbonds
  194.5 / Rg 20.4 exactly constant across 6.0-6.6, moves only in 7.0.
- **Confirming re-run (out2) COMPLETE**: box 100.7x100.7x123.7, 65 ion pairs, realized 0.154 M
  (equilibrated DDM t=52.1 A; used 48.8 for count -> ~3% high; use ~52 for exact 0.15), 0/136 ions
  outside DDM footprint, 0/136 in hydrophobic core, folded. Wrapped review PDB re-staged to ~/Downloads.

TODO: P8 (after user OK on the wrapped review PDB) re-run 4 glpG variants on cluster.

---
## Prior — box-fix DDM REMD: extract VTF+ΔG, repair disk-full corruption (2026-07-30)

Box-fix REMD (4 DDM variants) ran ~33 h (26/27 chunks) then died rc=1 on a full /project disk.
Raws intact (~2.28 TB). Extraction job 52780633 (raws kept) COMPLETED with mixed results:
- **79HIS**: VTF ✅ (227 MB) + ΔG ✅ — clean.
- **79ALA**: VTF ✅ (329 MB) + ΔG ✗ (`zero-size array` — NaN final chunk).
- **79HIS_S115T**: VTF ✗ + ΔG ✗ (`NaN in frame 0 of /output`).
- **79ALA_S115T**: VTF ✗ + ΔG ✅.
Root cause: the disk-full death left NaN in the last (interrupted) `/output` chunk of some
replicas, breaking the VTF extractor and the protection/MBAR. The ~26 completed chunks are fine.

Downloaded the 4 GOOD outputs to `~/Downloads` (79HIS VTF+ΔG, 79ALA VTF, 79ALA_S115T ΔG).
Submitted repair+re-extract **52828371** (caslake, 16 cpu, 14 h): per replica drop any output
group with non-finite pos (guaranteeing a valid `/output`), then re-do VTF+ΔG for the 3 affected
variants (79ALA, 79HIS_S115T, 79ALA_S115T). Raws NOT deleted — awaiting user verification.
No simulation physics or C++ changed.

## Prior phase summary — glpG HDX: REMD-system fix + lipid de-hardcode (2026-07-25)

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
