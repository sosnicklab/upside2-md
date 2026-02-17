# Task Plan

## Project Goal
- Define and implement a hybrid dry-MARTINI + Upside workflow where protein dynamics are handled by Upside, environment is handled by dry MARTINI, and coupling starts only at production stage.

## Architecture & Key Decisions
- Data preparation includes packing OPM protein `pdb/1rkl.pdb` into MARTINI DOPC bilayer `pdb/bilayer.MARTINI.pdb`, with MARTINI 2.2 coarse-graining and 0.15 M NaCl ion placement.
- Protein backbone state is represented by Upside atoms (`N, CA, C, O`) as primary DOFs.
- MARTINI protein backbone beads (`BB`) are computed from COM of mapped `N, CA, C, O` atoms each cycle; no per-step rebuilding of atoms from BB.
- Dry-MARTINI forces on protein `BB`/mapped sidechain proxies are projected onto protein atoms/DOFs in Upside; Upside then applies its own force field update.
- Hybrid force exchange is disabled before production stage; pre-production keeps protein rigid to relax environment.
- In production-stage hybrid mode, MARTINI intra-protein interactions are selectively filtered to avoid double-counting protein internal physics already modeled by Upside: no BB-BB, no SC-SC, and BB-SC only for same residue.
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
- [x] Phase 8: Audit pre-production rigidity and potential-output semantics in `run_sim_1rkl.sh` workflow artifacts.
- [x] Phase 9: Implement split potential logging (`Upside-only` vs `MARTINI`) and enforce pre-production NPT coupling to match CHARMM-GUI stage settings.

## Known Errors / Blockers
- Requires concrete integration points in Upside C++ step loop and force modules.
- Rigid-body alignment scope confirmed: apply alignment to coupling coordinates only; preserve raw integrated coordinates in saved trajectories.
- Phase 4 assumes valid BB mapping indices into current position array (`bb_atom_index` preferred; residue/name inference fallback may fail on nonstandard topologies).
- Phase 5 probabilistic sidechain coupling is implemented in C++; remaining work is physical calibration/stability validation for production trajectories.
- Phase 6 alignment uses BB-frame rigid transform between consecutive steps (not full Kabsch over all protein atoms); validate stability impact in Phase 7.
- Phase 7 completed smoke/integrity validation; full long-horizon physical validation remains to be benchmarked on production trajectories.
- Packed-system mapping export currently uses MARTINI protein index-space proxy targets (BB self-targets; SC proxies targeting residue BB). Full direct Upside `N,CA,C,O` target index export into production `.up` remains a separate integration task.
- Legacy static protein CG assets (`pdb/1rkl.MARTINI.pdb`, `pdb/1rkl_proa.itp`) are not reliable for dry-MARTINI compatibility; workflow mitigation is runtime MARTINI2 generation via `martinize.py -ff martini22`.
- Stage 6.0 HDF5 close-ID runtime failure is resolved (hybrid-control parsing bug in C++ fixed); remaining workflow issue is downstream VTF extraction mismatch when using `extract_martini_vtf.py` with a non-runtime structure source.
- Current `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up` has no `rotamer` or `placement*_point_vector_only` nodes, so live probabilistic rotamer weighting is inactive there and hybrid sidechain rows remain static (`rotamer_id=0` fallback).
- After removing SC placement/probability fallback, production files lacking `rotamer` + `placement*_point_vector_only` nodes will fail fast in hybrid mode and must be regenerated with those nodes present.
- Current stage-generation workflow (`prepare_martini.py` in `run_sim_1rkl.sh` / `test_prod_run_sim_1rkl.sh`) does not emit Upside `affine_alignment`/`placement_*`/`rotamer` nodes, so production-stage strict SC coupling fails unless scripts augment stage-7 files with these nodes.
- Stage-7 script augmentation now injects required `affine_alignment`/`placement_fixed_*`/`rotamer` nodes with `ff_2.1/sidechain.h5` compatibility; previous production SIGBUS (degenerate affine triplets) is fixed. Remaining blocker is production-stage physical instability (rapid potential blow-up/NaN) rather than a runtime memory crash.
- `output/potential` is a total-engine scalar (sum over active potential nodes), not a protein-only all-atom Upside-backbone channel; this can vary during pre-production even with rigid protein because environment terms evolve.

## Revised Decisions
- Hybrid coupling starts only at production stage; pre-production protein remains rigid.
- Pre-production rigid hold keeps the entire protein fully fixed in space (`x,y,z`) while lipid `PO4` headgroup beads are constrained only along `z` (`x,y` remain free); other lipid beads remain unconstrained.
- MARTINI intra-protein interactions in production hybrid mode follow selective rules:
  - no BB-BB (bonded or nonbonded),
  - no SC-SC (bonded or nonbonded),
  - BB-SC allowed only for same residue.
- BB MARTINI positions are derived each cycle from COM of Upside `N,CA,C,O`; no inverse backbone reconstruction from BB.
- Deterministic single-rotamer sidechain projection in C++ is removed; sidechain coupling uses per-rotamer probabilistic mapping rows from `hybrid_sc_map`.
- Sidechain probabilistic projection/force evaluation is implemented in Upside C++ runtime and is active only during production-stage hybrid coupling; preparation scripts remain unchanged for pre-production rigid setup.
- If enabled, frame-to-frame RMSD alignment is used only in coupling frame construction and does not overwrite primary trajectory coordinates.
- Packed-system hybrid mapping artifacts must keep all BB/SC projection indices in the same runtime coordinate space as `protein_membership`; validator/runtime now enforce protein-local bounds.
- Hybrid sidechain placement in production will be sourced from Upside rotamer placement coordinates (`placement_fixed_point_vector_only` / `placement_point_vector_only`) per rotamer state when available, with runtime fallback to `hybrid_sc_map` projected rows if placement nodes are absent.
- Probabilistic sidechain weights in production will use live Upside rotamer marginals from C++ runtime rotamer state when available, with fallback to static `hybrid_sc_map/rotamer_probability` only when rotamer marginals are unavailable.
- Hybrid runtime may add ordering dependencies so rotamer/placement sidechain states are computed before MARTINI coupling each step; this is a scheduling change only (no Python preparation changes).
- Sidechain-environment force transfer in hybrid production prefers explicit `hybrid_sc_map` projection targets/weights generated from martinize bonded topology, with residue-based MARTINI `BB` target lookup retained only as fallback.
- Production-stage MARTINI protein interaction filtering is runtime-enforced in C++ (based on `protein_membership`, `atom_roles`, and `residue_ids`) so non-production stages remain unaffected while protein is rigid.
- Production-stage MARTINI nonbonded interactions for non-protein particles are hard-sphere-like (repulsive-only) while protein-involving MARTINI interactions remain dry-MARTINI for hybrid coupling.
- Preparation exports per-`BB` backbone reference metadata (`N,CA,C,O` atom indices/coordinates plus comment strings) in `hybrid_bb_map` for auditability; this metadata does not override runtime coordinate arrays unless explicitly mapped into runtime index space.
- Runtime hybrid scheduling enforces `martini_potential` execution after all coordinate-node updates each integration sub-step, so MARTINI coupling always consumes the latest backbone-derived coordinates for that sub-step.
- `run_sim_1rkl.sh` defaults to MARTINI2.2 protein coarse-graining via local `martinize.py` (`-ff martini22`) and performs a preflight protein-ITP-vs-dry-FF mass compatibility check to prevent silent MARTINI3/dry-MARTINI mixing.
- Hybrid rotamer/placement node binding is discovered from the built `DerivEngine` graph at registration time (not by pre-scanning HDF5 `/input/potential`), to avoid HDF5 handle/path edge cases in hybrid-control parsing.
- Bilayer preparation for hybrid packing must tile/crop lipid coordinates in XY so lipid edges extend to at least one protein-length margin beyond protein edges on each side, and the XY box edges are set exactly to lipid edges.
- Ion placement for target salinity uses a template-calibrated effective ion-accessible volume fraction (derived from bilayer template ions/box) to avoid full-box overestimation in dry-MARTINI packed systems.
- Hybrid bilayer packing enforces square XY dimensions in the generated system by expanding the shorter side to the longer side after minimum protein-margin constraints are applied.
- Hybrid packed PDB output canonicalizes lipid residue naming to `DOPC` (accepting `DOP`/`DOPC` as inputs) and preserves 4-character residue-name formatting.
- Hybrid box sizing enforces a minimum `Z` length of `3x` the protein `Z` span in generated packed systems.
- Hybrid packing places the lipid bilayer midplane at the box `Z` center in generated PDBs.
- Pre-production rigid hold (protein + configured lipid headgroup roles, default `PO4`) is fixed in space, not only force-constrained.
- Pre-production fixed-in-space enforcement uses an integration-level fixed-atom mask (skip momentum/position updates for fixed atoms) plus force/momentum zeroing; direct per-step coordinate rewrites are avoided to prevent early-step NaN instability.
- Lipid-overlap removal for protein insertion uses a tighter default cutoff (`3.0 Å`) to avoid creating an artificially large annular cavity around the protein.
- NPT barostat updates in hybrid pre-production rigid mode remain enabled; scaling is applied only to unconstrained DOFs (protein fixed atoms are not scaled, z-fixed headgroups remain z-fixed).
- `run_sim_1rkl.sh` converts requested MD steps to Upside CLI time units (`duration`, `frame-interval`) to avoid unintended 100x step inflation.
- VTF export is stage-wise: modes are `mode 1` (all MARTINI particles) for non-production stages and `mode 2` (non-protein MARTINI + protein all-atom backbone backmapped from `hybrid_bb_map`) for production stage.
- Upside CLI supports `--duration-steps` for step-count-driven stage lengths, and `run_sim_1rkl.sh` uses it so stage step counts are independent of `dt`.
- Sidechain-to-backbone force transfer mapping will prefer explicit `hybrid_sc_map` projection targets/weights (prepared from martinize.py protein bonded topology/force constants) with BB-residue lookup retained only as runtime fallback.
- Minimization mode no longer forces a post-minimization stage switch to `production`; runtime now restores the pre-minimization stage label to prevent production-only hybrid/UpSide activation in pre-production stage files.
- Production-stage hard-sphere behavior is implemented in C++ MARTINI runtime as repulsive-only WCA-like nonbonded interactions for non-protein/non-protein pairs; this is controlled by `hybrid_control` and active only when hybrid is active.
- Sidechain placement/probability fallback is removed in production hybrid mode: SC mapping requires live Upside `rotamer` and `placement*_point_vector_only` nodes, and static `hybrid_sc_map` fallback is no longer used for SC position/probability evaluation.
- SC rotamer projection in production uses a residue-level rigid transform between reference and current rotamer placement frames so mapped SC geometry remains rigid across rotamer states.
- Production-stage workflow scripts will explicitly augment stage-7 `.up` files with `affine_alignment`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, and `rotamer` nodes (built from sidechain library + current MARTINI protein mapping) before running production.
- Pre-production rigid hold verification should use runtime hybrid metadata (`hybrid_control` + `protein_membership`) and trajectory displacement checks; presence of `/input/fix_rigid` datasets is not required because pre-production fixed masks are generated dynamically in C++ runtime.
- Potential logging will be split:
  - `output/potential`: non-MARTINI (Upside-side) potential contribution only.
  - `output/martini_potential`: MARTINI potential contribution, preserving visibility of the previously observed scalar behavior for MARTINI stages.
- Pre-production NPT coupling should follow CHARMM-GUI stage intent:
  - stages `6.0-6.1`: `ref_p = 0.0 bar`, `tau_p = 4.0`, `compressibility = 3e-4 bar^-1 (xy), 0.0 (z)`.
  - stages `6.2-6.6`: `ref_p = 1.0 bar` (semi-isotropic; `z` compressibility remains `0.0`), `tau_p = 4.0`, same compressibility.
  - convert pressure to Upside units using `1 bar = 0.000020659477 E_up/Angstrom^3`.
- Hybrid pre-production rigid hold and NPT are both enabled by scaling only unconstrained DOFs during barostat box/coordinate updates (fully fixed atoms remain fixed; z-fixed atoms keep fixed z).
