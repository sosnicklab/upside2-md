# Task Plan

## Project Goal
- Define and implement a hybrid dry-MARTINI + Upside workflow where protein dynamics are handled by Upside, environment is handled by dry MARTINI, and coupling starts only at production stage.

## Architecture & Key Decisions
- Data preparation includes packing OPM protein `pdb/1rkl.pdb` into MARTINI DOPC bilayer `pdb/bilayer.MARTINI.pdb`, with MARTINI 2.2 coarse-graining and 0.15 M NaCl ion placement.
- Protein backbone state is represented by Upside atoms (`N, CA, C, O`) as primary DOFs.
- MARTINI protein backbone beads (`BB`) are computed from COM of mapped `N, CA, C, O` atoms each cycle; no per-step rebuilding of atoms from BB.
- Dry-MARTINI forces on protein `BB`/mapped sidechain proxies are projected onto protein atoms/DOFs in Upside; Upside then applies its own force field update.
- Hybrid force exchange is disabled before production stage; pre-production keeps protein rigid to relax environment.
- All MARTINI interactions within the same protein are disabled to avoid double-counting protein internal physics already modeled by Upside.
- Rotamer mapping is used only for protein-environment coupling, not protein-protein or intra-protein sidechain interactions.
- Hybrid `.up` schema will add dedicated groups under `/input`: `hybrid_control`, `hybrid_bb_map`, `hybrid_sc_map`, and `hybrid_env_topology`.
- Add optional rigid-body drift control for Upside protein motion in production: RMSD alignment or equivalent COM+rotation removal mode, enabled by configuration and evaluated for physical side effects.

## Execution Phases
- [x] Phase 1: Finalize hybrid data schema for `.up` and packed-system metadata.
- [x] Phase 2: Implement data-preparation script (CG conversion, bilayer packing, ions, mapping export).
- [x] Phase 3: Implement Upside input parsing for hybrid mapping blocks and stage-gated runtime flags.
- [x] Phase 4: Implement hybrid coupling loop in Upside (MARTINI->Upside force projection, BB COM refresh, production-only activation).
- [x] Phase 5: Implement sidechain rotamer-environment mapping path with intra-protein MARTINI exclusions.
- [x] Phase 6: Implement/benchmark rigid-body drift mitigation mode (frame-to-frame alignment option) and diagnostics.
- [x] Phase 7: Validate correctness, stability, and energy/force accounting; document usage and limitations.

## Known Errors / Blockers
- Requires concrete integration points in Upside C++ step loop and force modules.
- Requires confirmation of available martinize2 command and MARTINI 2.2 assets in runtime environment.
- Rigid-body alignment scope confirmed: apply alignment to coupling coordinates only; preserve raw integrated coordinates in saved trajectories.
- Phase 4 assumes valid BB mapping indices into current position array (`bb_atom_index` preferred; residue/name inference fallback may fail on nonstandard topologies).
- Phase 5 currently uses one deterministic sidechain proxy entry per MARTINI sidechain bead (`rotamer_probability=1.0`); full Upside rotamer-probability coupling remains future work.
- Phase 6 alignment uses BB-frame rigid transform between consecutive steps (not full Kabsch over all protein atoms); validate stability impact in Phase 7.
- Phase 7 completed smoke/integrity validation; full long-horizon physical validation remains to be benchmarked on production trajectories.

## Revised Decisions
- Hybrid coupling starts only at production stage; pre-production protein remains rigid.
- Pre-production rigid hold includes entire protein plus only lipid `PO4` headgroup beads (not `NC3/GL1/GL2`).
- MARTINI intra-protein interactions are globally excluded for this hybrid mode.
- BB MARTINI positions are derived each cycle from COM of Upside `N,CA,C,O`; no inverse backbone reconstruction from BB.
- If enabled, frame-to-frame RMSD alignment is used only in coupling frame construction and does not overwrite primary trajectory coordinates.
